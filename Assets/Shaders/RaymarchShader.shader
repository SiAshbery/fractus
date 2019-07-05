Shader "Fractus/RaymarchShader"
{
    Properties
    {
        _MainTex ("Texture", 2D) = "white" {}
    }
    SubShader
    {
        // No culling or depth
        Cull Off ZWrite Off ZTest Always

        Pass
        {
            CGPROGRAM
            #pragma vertex vert
            #pragma fragment frag
            // prafgma target sets the target shader model that the HLSL will be compiled to.
            // The higher themore modern.
            // 3.0 supportsL:
            // DX9 shader model 3.0: derivative instructions, texture LOD sampling, 10 interpolators, more math/texture instructions allowed.
            // Not supported on DX11 feature level 9.x GPUs (e.g. most Windows Phone devices).
            // Might not be fully supported by some OpenGL ES 2.0 devices, depending on driver extensions present and features used.
            #pragma target 3.0

            #include "UnityCG.cginc"
            // Include distance functions library.
            #include "DistanceFunctions.cginc"
            
            sampler2D _MainTex;

            // Use the depth texture to render meshes in front of SDFs correctly.
            // When the ray hits the distance value of the depth texture before it has hit the surface of 
            // an SDF we know that the mesh is in front of the SDF, so we should stop marching rays and draw the mesh.
            uniform sampler2D _CameraDepthTexture;

            // To create a raymarching shader we need to get the origin position of the ray which is the camera's world space position.
            // Unity has a built in variable for this called _WorldSpaceCameraPos
            // In the on render image we'll send the camera position to a float 4 in this shader

            // We also need the ray direction based on the field-of-view and aspect ratio of the camera
            // In C# we will calculate a direction for each corner of the screen
            // and send those directions in a matrix to the shader

            // The CamFrustum will be calculated in eye space and we will need to convert that to woprld space
            // So we need a transformation matrix to send to the shader. And that will be the CamToWorldTransformation
            // Which is supplied directly from the camera
            uniform float4x4 _CamFrustum, _CamToWorld;
            // Set a max distance for rays to travel before they are abandoned
            uniform float _maxDistance;
            // Set a max distance for rays to travel before they are abandoned
            uniform int _maxIterations;
            //lighting
            uniform float3 _lightDir;
            uniform fixed4 _lightCol;
            // Set shape attributes in editor
            uniform float4 _sphere1, _box1, _recursiveTet1, _mandelBulb1, _mandelBox1;
            uniform float _recursiveTet1Offset;
            uniform int _recursiveTet1Iterations, _mandelBulb1Iterations, _mandelBox1Iterations;
            uniform float _mandelBulb1Power, _mandelBulb1Bailout, _mandelBox1Scale, _mandelBox1SphereRadius;
			// Set repetation interval for modulus
			uniform float3 _modInterval, _mandelBox1FoldLimit;

            struct appdata
            {
                float4 vertex : POSITION;
                float2 uv : TEXCOORD0;
            };

            struct v2f
            {
                float2 uv : TEXCOORD0;
                float4 vertex : SV_POSITION;
                // stores the ray direction
                float3 ray : TEXCOORD1;
            };

            v2f vert (appdata v)
            {
                v2f o;
                // ??
                half index = v.vertex.z;
                v.vertex.z = 0;

                o.vertex = UnityObjectToClipPos(v.vertex);
                o.uv = v.uv;

                // Set ray to corresponding cam frustum
                o.ray = _CamFrustum[(int)index].xyz;

                // normalize ray in z direction
                o.ray /= abs(o.ray.z);

                // convert the ray from eye space to world space
                o.ray = mul(_CamToWorld, o.ray);

                return o;
            }

            float distanceField(float3 p)
            {
				// p = position
				// Repeat along axis:
				// First argument is the position and is an inout which means it is changed by the function
				// Second argument is the size of the repeated chunk and should be double the size of the shape
				// I.E a Cube that is 2x2x2 would need a Mod size of 4 to fit perfectly into the repeation.
				float modX = pMod1(p.x, _modInterval.x);
				float modY = pMod1(p.y, _modInterval.y);
				float modZ = pMod1(p.z, _modInterval.z);
				float Sphere1 = sdSphere(p - _sphere1.xyz, _sphere1.w);
				float Box1 = sdBox(p - _box1.xyz, _box1.www);
                return opS(Sphere1,sdMandelBox(p - _mandelBox1.xyz, _mandelBox1Iterations, _mandelBox1Scale, _mandelBox1SphereRadius, _mandelBox1FoldLimit.xyz));
                //return sdMandelBulb(p - _mandelBulb1.xyz, _mandelBulb1Power, _mandelBulb1Bailout, _mandelBulb1Iterations);
				//return opS(Sphere1,Box1);
                //return sdRTet(p - _recursiveTet1.xyz, _recursiveTet1.w,_recursiveTet1Offset, _recursiveTet1Iterations);
            }

            float3 getNormal(float3 p)
            {
                // p = position
                // this function gets the normal value of the current position,
                // we can use this in our lighting model.

                // Getting the normal in a distance field is slightly different than doing so for a 3D mesh
                // The gradient of the distance field is the same as the normal at that point
                // The gradient is the dirivative of the field in the x,y and z directions.

                // To calculate the normal we offset the current position by a small number
                // We calculate the position + the offset and next the position - the offset
                // We subtract those results from each other for each axis
                
                // This means we have the calculate the distance field 6 extra times to get the normal

                const float2 offset = float2(0.001, 0.0);

                // calculate the normal for each of the axis of p
                // the offset calls with .xyy etc are to quickly create
                // a float3, the x of offset has a value of 0.001 and the y has 0.0
                // so this is a shorthand way of writing float3(offset.x, offset.y, offset.y)
                float3 n = float3(
                    distanceField(p + offset.xyy) - distanceField(p - offset.xyy),
                    distanceField(p + offset.yxy) - distanceField(p - offset.yxy),
                    distanceField(p + offset.yyx) - distanceField(p - offset.yyx));
                return normalize(n);
            }

            fixed4 raymarching(float3 ro, float3 rd, float depth)
            {
                // ro = ray origin, rd = ray direction
                fixed4 result = fixed4(0,0,0,1);

                // to implement sphere trasing we need to make a loop
                // so that it can march along the ray until we hit a surface

                // track distance traveled along ray direction
                float t = 0;

                for (int i = 0; i < _maxIterations; i++)
                {
                    // Drawn the environment if we have gone too far without a collision.
                    // or if the distance we have gone beyond the depth value for a mesh in the unity scene.
                    if (t > _maxDistance || t >= depth )
                    {
                        // for rays that miss the distance filed we need to set the alpha value to 0
                        // So that we can render the unity scene
                        result = fixed4(rd, 0);
                        break;
                    }

                    // get current position
                    // th origin vector added to the direction gives us our direction of travel which we multiply by how far we have gone to ge the current position.
                    float3 p = ro + rd * t;
                    // Check for a collision with a distance field
                    // d = closest distance to a surface
                    float d = distanceField(p);

                    // if the result of d is negative e.g. -1 We are instande an object
                    // if it is positive e.g. 1, we are outside (this is our distance from a surface)
                    // if it is 0 (or within a tolerance e.g. 0.01) we are at the surface.
                    if (d < 0.01)
                    {
                        // shading! calculating normals and a lambertian light model, fun!
                        // normals!
                        float3 n = getNormal(p);
                        // light!
                        // Lighting requires the dot product of the inversed lighting direction and the normal direction
                        float light = dot(-_lightDir, n);
                        result = fixed4(_lightCol.rgb * light,1);
                        break;
                    }

                    // If we have not met any break criteria, track the distance traveled
                    t += d;
                }

                return result;
            }

            fixed4 frag (v2f i) : SV_Target
            {
                //find the current depth of the univty scene for this ray.
                float depth = LinearEyeDepth(tex2D(_CameraDepthTexture, i.uv).r);
                depth *= length(i.ray);
                // We need to assign the colors of the MainTex so we can blend between the raymarching effects and the unity scene.
                fixed3 col = tex2D(_MainTex, i.uv);
                // set ray direction
                float3 rayDirection = normalize(i.ray.xyz);
                // Set the ray origin
                float3 rayOrigin = _WorldSpaceCameraPos;

                // get color values from raymarching function
                fixed4 result = raymarching(rayOrigin,rayDirection,depth);
                // Check to see if we have hit the distance field or not and show the raymarch shader or unity scene accordingly
                // To do so we will use the w value of the returning result, when there is a hit we set the alpha component to 1
                // for rays that miss we need to set the alpha value to 0
                // If the ray has hit something in the distance field the values of result.w will be 1
                // So 1.0 - 1.0 would be 0, so we don't render the col value.
                // We then add the result color value multiplied by result w which if it hasn't hit anything is 0 so that would in turn resulte as 0
                return fixed4(col * (1.0 - result.w) + result.xyz * result.w,1.0);
            }
            ENDCG
        }
    }
}
