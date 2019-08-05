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

            uniform sampler2D _patternTex;

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
			uniform float _accuracy;
            //lighting
            uniform float3 _lightDir;
            uniform fixed4 _lightCol;
			uniform fixed4 _mainCol, _fogColor, _glowColor;
			uniform float _lightIntensity, _fogDensity;
			uniform float _shadowIntensity, _specularHighlight, _specularIntensity, _glowSharpness;
			uniform float2 _shadowDistance;
			uniform float _shadowPenumbra;
            // Set shape attributes in editor
            uniform float4 _sphere1, _box1, _recursiveTet1, _apollonian1;
            uniform float3 _mandelBox1, _mandelBulb1;
            uniform float _recursiveTet1Offset;
            uniform int _recursiveTet1Iterations, _mandelBulb1Iterations, _mandelBox1Iterations, _apollonian1Iterations;
            uniform float _mandelBulb1Power, _mandelBulb1Bailout, _mandelBox1Scale, _mandelBox1SphereRadius, _apollonian1Scale;
			// Set repetation interval for modulus
			uniform float3 _modInterval, _mandelBox1FoldLimit, _apollonian1Size;

			uniform float _repeatX, _repeatY, _repeatZ;

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
				float time = _Time;
				// p = position
				// Repeat along axis:
				// First argument is the position and is an inout which means it is changed by the function
				// Second argument is the size of the repeated chunk and should be double the size of the shape
				// I.E a Cube that is 2x2x2 would need a Mod size of 4 to fit perfectly into the repeation.
				// float modX = pMod1(p.x, _modInterval.x);
				// float modY = pMod1(p.y, _modInterval.y);
				// float modZ = pMod1(p.z, _modInterval.z);
			
				//float Sphere1 = sdSphere(p - _sphere1.xyz, _sphere1.w);
				//float Box1 = sdBox(p - _box1.xyz, _box1.www);
				//float mandelBulb = sdMandelBulb(p - _mandelBulb1.xyz, _mandelBulb1Power, _mandelBulb1Bailout, _mandelBulb1Iterations);
                //float mandelBox = sdMandelBox(p - _mandelBox1.xyz, _mandelBox1Iterations, _mandelBox1Scale + abs(sin(time) * cos(time)), _mandelBox1SphereRadius, _mandelBox1FoldLimit.xyz);
				//_apollonian1Size.x += abs(sin(time) * cos(time)) * 0.002;
				//_apollonian1Size.z -= abs(sin(time) * cos(time)) * 0.002;
                float apollonian = sdApollonian(p - _apollonian1.xyz, _apollonian1Scale, _apollonian1Iterations, _apollonian1Size);
                //float sierpinskiTri = sdRTet(p - _recursiveTet1.xyz, _recursiveTet1.w,_recursiveTet1Offset, _recursiveTet1Iterations);
                //float julia = sdJulia(p);
				//return Sphere1;
				return apollonian;
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

			uniform float _aoStepSize, _aoIntensity;
			uniform int _aoIterations;

			float ambientOcclusion(float3 p, float3 n) {
				float step = _aoStepSize;
				float ao = 0.0;
				float dist;
				for (int i = 1; i <= _aoIterations; i++) {
					dist = step * i;
					ao += max(0.0, (dist - distanceField(p + n * dist)) / dist);
				}
				return (1.0 - ao * _aoIntensity);
			}

			float hardShadow(float3 ro, float3 rd, float mint, float maxt) {
				// t = distance travelled. 
				// We start at our minimum distance and we proceed along the inverse of the lighting direction
				// and continue checking our distance field. If we have a collision we know we have hit something and we
				// return 0 which will be multiplied by the light value to give a shadow.
				// if we get to our max distance(maxt) we have no shadows and return 1.0
				for (float t = mint; t < maxt;) {
					float h = distanceField(ro + rd * t);
					if (h < _accuracy) {
						return 0.0;
					}
					t += h;
				}
				return 1.0;
			}

			float softShadow(float3 ro, float3 rd, float mint, float maxt, float k) {
				// k = penumbra
				float result = 1.0;
				for (float t = mint; t < maxt;) {
					float h = distanceField(ro + rd * t);
					if (h < _accuracy) {
						return 0.0;
					}
					// instead of returning 0 for the rays that dont hit an object
                    // we return the value of the closest distance it came to an object
					result = min(result, k * h / t);
					t += h;
				}
				return result;
			}
            

            float3 applySimpleFog( in float3  rgb,       // original color of the pixel
               in float d // camera to point distance
               )
            {
                //float direction = normalize(d);
                //fog density
                //float b = 0.25;
                float fogAmount = 1.0 - exp( -(abs(d))*(_fogDensity) );
                // glsl uses mix, hlsl uses lerp
                return lerp( rgb, _fogColor, fogAmount );
            }

            float3 applyLitFog( in float3  rgb,      // original color of the pixel
               in float d, // camera to point distance
               in float3  rd)  // camera to point vector
            {
                float direction = normalize(rd);
                float fogAmount = direction - exp( -d*_fogDensity );
                float sunAmount = max( dot( rd, _lightDir ), 0.0 );
                float3 color  = lerp( _fogColor, // bluish
                           _lightCol, // yellowish
                           pow(_lightIntensity / 4,8.0) );
                return lerp( rgb, color, fogAmount );
            }

			float3 applySpecular(float3 col, float3 reflected, float shadow) {
				float specular = max(dot(reflected, _lightDir), 0.0);
				specular = pow(specular, _specularHighlight);
				col += specular * _lightCol * shadow * _specularIntensity;
				return col;
			}
			
			float3 applyGlow(float3 col, float minD) {
                if (_glowSharpness > 0) {
			        col += (1.0 - minD) * (1.0 - minD) * _glowColor;
                }
			    return col;
			}

            float4 texcube( sampler2D sam, in float3 p, in float3 n ) {
                float4 x = tex2D( sam, p.yz );
                float4 y = tex2D( sam, p.zx );
                float4 z = tex2D( sam, p.xy );
                return x*abs(n.x) + y*abs(n.y) + z*abs(n.z);
            }

            float3 applySurfaceCol(float3 col, float3 p, float3 n) {
                float3 result = texcube( _patternTex, 0.5*p, n ).xyz;
                //col = col*(0.2+0.8*texcube( iChannel2, 4.0*vec3(2.0,8.0,2.0)*pos, nor ).x);
                //vec3 verde = vec3(1.0,0.9,0.2);
                //verde *= texture( iChannel2, pos.xz ).xyz;
                //col = mix( col, 0.8*verde, hh );
        
                float vv = smoothstep( 0.0, 0.8, n.y )*smoothstep(0.0, 0.1, p.y-0.8 );
                //verde = vec3(0.2,0.45,0.1);
                //verde *= texture( iChannel2, 30.0*pos.xz ).xyz;
                //verde += 0.2*texture( iChannel2, 1.0*pos.xz ).xyz;
                vv *= smoothstep( 0.0, 0.5, tex2D( _patternTex, 0.1*p.xz + 0.01*n.x ).x );
                //col = mix( col, verde*1.1, vv );
                result = lerp(col,result,vv);
                return result;
            }

			float3 Shading(float3 p, float3 n, float t, float3 ray, float minD) {
				float3 result;

				float3 reflected = ray.xyz - 2.0 * dot(ray.xyz, n) * n;

				// Diffuse color
				float3 color = _mainCol.rgb;

                color = applySurfaceCol(color, p, n);
				// Directional light
				float3 light = (_lightCol * dot(-_lightDir, n) * 0.5 + 0.5) * _lightIntensity;
				// Shadows
				float shadow = softShadow(p, -_lightDir, _shadowDistance.x, _shadowDistance.y, _shadowPenumbra) * 0.5 + 0.5;
				shadow = max(0.0,pow(shadow, _shadowIntensity));
				float ao = ambientOcclusion(p, n);
               
				result = color * light * shadow * ao;
            
				if (_specularHighlight > 0) {
					result = applySpecular(result, reflected, shadow);
				}
				
				
				//result = applyGlow(result, minD);
				
                //result = applyTexture(result, p);
                float r2 = dot(p,p);
                float4 orb = float4(1000.0,1000.0,1000.0,1000.0);
                orb = min( orb, float4(abs(p),r2));
        
                //result += applySurfaceCol(result, orb, p);
				
				// Fog should always be calculated last
				if (_fogDensity > 0) {
                    result = applySimpleFog(result, t);
                }


				return result;


			}

  

            fixed4 raymarching(float3 ro, float3 rd, float depth, float2 uv)
            {
                // min distance to surface
	            float minD = 1.0;
                // ro = ray origin, rd = ray direction
                float3 p;

                int didHit = 0;

                //Screen space texture
                //fixed4 result = fixed4(tex2D(_patternTex, uv).xyz,1);

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
                        float alpha;
                     

                        // Inversed normalization of the minimum distance.
                        // The closer it is to the surface 
                        if (_glowSharpness > 0) {
                            alpha = max(0, 1 - (minD*_glowSharpness));
                        } else {
                            alpha = 0;
                        }

                       
                        result = fixed4(applyGlow(float3(0,0,0), minD),alpha);
                      
                        break;
                    }

                    // get current position
                    // th origin vector added to the direction gives us our direction of travel which we multiply by how far we have gone to ge the current position.
                    p = ro + rd * t;
                    // Check for a collision with a distance field
                    // d = closest distance to a surface
                    float d = distanceField(p);

                    // if the result of d is negative e.g. -1 We are instande an object
                    // if it is positive e.g. 1, we are outside (this is our distance from a surface)
                    // if it is 0 (or within a tolerance e.g. 0.01) we are at the surface.
                    if (d < _accuracy)
                    {
                        didHit = 1;
                        break;
                    }

                    // If we have not met any break criteria, track the distance traveled
                    t += d;
                    minD = min(minD, _glowSharpness * d / t);
                }

                if (didHit > 0) {      
                    // shading! calculating normals and a lambertian light model, fun!
                    // normals!
                    float3 n = getNormal(p);
                    // light!
                    // Lighting requires the dot product of the inversed lighting direction and the normal direction
                    float3 s = Shading(p, n, t, -rd, minD);
                    result = fixed4(s, 1);
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
                fixed4 result = raymarching(rayOrigin,rayDirection,depth, i.uv);
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
