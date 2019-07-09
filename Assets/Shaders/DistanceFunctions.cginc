// Sphere
// s: radius
float sdSphere(float3 p, float s)
{
	return length(p) - s;
}

// Box
// b: size of box in x/y/z
float sdBox(float3 p, float3 b)
{
	float3 d = abs(p) - b;
	return min(max(d.x, max(d.y, d.z)), 0.0) +
		length(max(d, 0.0));
}

// Octahedron
float sdOctahedron(float3 p, float s)
{
    p = abs(p);
    float m = p.x+p.y+p.z-s;
    float3 q;
         if( 3.0*p.x < m ) q = p.xyz;
    else if( 3.0*p.y < m ) q = p.yzx;
    else if( 3.0*p.z < m ) q = p.zxy;
    else return m*0.57735027;
    
    float k = clamp(0.5*(q.z-q.y+s),0.0,s); 
    return length(float3(q.x,q.y-s+k,q.z-k)); 
}

// OPERATORS

void sphereFold(inout float3 z, inout float dz, float3 radius) {
    float minRadius2 = radius.x;
    float fixedRadius2 = radius.y;
    float r = dot(z,z);
    float r2 = dot(z,z);
    if (r<minRadius2) { 
        // linear inner scaling
        float temp = (fixedRadius2/minRadius2);
         z = mul(z,temp);
        dz = mul(dz,temp);
    } else if (r2<fixedRadius2) { 
        // this is the actual sphere inversion
        float temp =(fixedRadius2/r2);
        z = mul(z,temp);
        dz = mul(dz,temp);
    }
}

void boxFold(inout float3 z, inout float dz, int foldingLimit) {
    z = clamp(z, -foldingLimit, foldingLimit) * 2.0 - z;
}

//Fractals

// Recursive tetrahedron
float sdRTet(float3 p, float scale, float offset, int iterations)
{
    int n = 0;
    while (n < iterations) {
       if(p.x+p.y<0) p.xy = -p.yx; // fold 1
       if(p.x+p.z<0) p.xz = -p.zx; // fold 2
       if(p.y+p.z<0) p.zy = -p.yz; // fold 3    
       p = p*scale - offset*(scale-1.0);
       n++;
    }
    return (length(p) ) * pow(scale, -float(n));
}

float sdMandelBulb(float3 p, float power, float bailout, float iterations) {
    float3 z = p;
    float dr = 0.5;
    float r = 0.0;
    for (int i = 0; i < iterations ; i++) {
        r = length(z);
        if (r>bailout) break;
        
        // convert to polar coordinates
        float theta = acos(z.z/r);
        float phi = atan2(z.y,z.x);
        dr = pow(r, power-1.0)*power*dr + 1.0;
        
        // scale and rotate the point
        float zr = pow(r,power);
        theta = theta*power;
        phi = phi*power;
        
        // convert back to cartesian coordinates
        z = zr*float3(sin(theta)*cos(phi), sin(phi)*sin(theta), cos(theta));
        z+=p;
    }
    return 0.5*log(r)*r/dr;
}

float sdMandelBox(float3 p, int iterations, float scale, float sphereRadius, float3 foldLimit)
{
    //int iterations = 12;
    //float scale = 10.0;
    // float3 offset = p;
    // float dr = 1.0;
    // for (int n = 0; n < iterations; n++) {
    //    boxFold(p,dr,foldLimit);       // Reflect
    //    sphereFold(p,dr,sphereRadius);    // Sphere Inversion
    //    
    //            p=scale*p + offset;  // Scale & Translate
    //            dr = dr*abs(scale)+1.0;
    //}
    //float r = length(p);
    //return r/abs(dr);



    float MR2 = sphereRadius * sphereRadius;
    // precomputed somewhere
    float4 scalevec = float4(scale, scale, scale, abs(scale)) / MR2;
    float C1 = abs(scale-1.0), C2 = pow(abs(scale), float(1-iterations));

    // distance estimate
    float4 pos = float4(p.xyz, 1.0), p0 = float4(p.xyz, 1.0);  // p.w is knighty's DEfactor
    for (int i=0; i<iterations; i++) {
    // foldLimit.xyz: good values to start from (1.0, 1.0, 2.0)
    pos.xyz = clamp(pos.xyz, -foldLimit.x, foldLimit.y) * foldLimit.z - pos.xyz;  // box fold: min3, max3, mad3
    float r2 = dot(pos.xyz, pos.xyz);  // dp3
    pos.xyzw = mul(pos.xyzw,clamp(max(MR2/r2, MR2), 0.0, 1.0));  // sphere fold: div1, max1.sat, mul4
    pos.xyzw = pos*scalevec + p0;  // mad4
  }
  return (length(pos.xyz) - C1) / pos.w - C2;
}

float sdApollonian( float3 p, float scale, int iterations, float3 size )
{
    //float scale = 1.0;

    //float4 orb = float4(1000.0,1000.0,1000.0,1000.0); 
    
    //for( int i=0; i<8;i++ )
    //{
    //    p = -1.0 + 2.0*frac(0.5*p+0.5);

    //    float r2 = dot(p,p);
        
    //    orb = min( orb, float4(abs(p),r2) );
        
    //    float k = s/r2;
    //    p     *= k;
    //    scale *= k;
    //}
    
    //return 0.25*abs(p.y)/scale;

    //float3 CSize = float3(1.0, 1.0, 1.3);
    p = p.xzy;
    // float scale = 0.5;
    for( int i=0; i < iterations ;i++ )
    {
        p = 2.0*clamp(p, -size, size) - p;
        float r2 = dot(p,p);
        //float r2 = dot(p,p+sin(p.z*.3)); //Alternate fractal
        float k = max((2.0)/(r2), 0.027);
        p     *= k;
        scale *= k;
    }
    float l = length(p.xy);
    float rxy = l - 4.0;
    float n = l * p.z;
    rxy = max(rxy, -(n) / 4.0);
    return (rxy) / abs(scale);
}

// BOOLEAN OPERATORS //

// Union
float opU(float d1, float d2)
{
	return min(d1, d2);
}

// Subtraction
float opS(float d1, float d2)
{
	return max(-d1, d2);
}

// Intersection
float opI(float d1, float d2)
{
	return max(d1, d2);
}

// Repeats a distance function along an axis.
// Mod Position Axis
float pMod1 (inout float p, float size)
{
	float halfsize = size * 0.5;
	float c = floor((p+halfsize)/size);
	p = fmod(p+halfsize,size)-halfsize;
	p = fmod(-p+halfsize,size)-halfsize;
	return c;
}
