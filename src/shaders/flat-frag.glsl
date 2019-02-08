#version 300 es
precision highp float;

uniform vec3 u_Eye, u_Ref, u_Up;
uniform vec2 u_Dimensions;
uniform float u_Time;

uniform highp float u_SlowFactor;
uniform vec4 u_Color;

in vec2 fs_Pos;
out vec4 out_Col;

const int MAX_MARCHING_STEPS = 255;
const float MIN_DIST = 0.0;
const float MAX_DIST = 100.0;
const float EPSILON = 0.0001;

float triangleWave(float x, float freq, float amplitude) {
  return abs(mod((x * freq), amplitude) - (0.5 * amplitude));
}

/**
 * SDF cominbation operations 
 */ 
float intersectOp(float distA, float distB) {
    return max(distA, distB);
}

float unionOp(float distA, float distB) {
    return min(distA, distB);
}

float differenceOp(float distA, float distB) {
    return max(distA, -distB);
}

/** 
 * Signed distance functions (SDF)
 */

float boxSDF(vec3 p, vec3 boxDim) {
    vec3 d = abs(p) - boxDim;
    float insideDistance = min(max(d.x, max(d.y, d.z)), 0.0);
    float outsideDistance = length(max(d, 0.0));
    return insideDistance + outsideDistance;
}

float sphereSDF(vec3 p, float r) {
    return length(p) - r;
}

float roundBoxSDF( in vec3 p, in vec3 b, in float r) {
    vec3 q = abs(p) - b;
    return min(max(q.x,max(q.y,q.z)),0.0) + length(max(q,0.0)) - r;
}

float roundConeSDF( in vec3 p, in float r1, float r2, float h ) {
    vec2 q = vec2( length(p.xz), p.y );
    
    float b = (r1-r2)/h;
    float a = sqrt(1.0-b*b);
    float k = dot(q,vec2(-b,a));
    
    if( k < 0.0 ) return length(q) - r1;
    if( k > a*h ) return length(q-vec2(0.0,h)) - r2;
        
    return dot(q, vec2(a,b) ) - r1;
}

float ellipsoidSDF( in vec3 p, in vec3 r ) {
    float k0 = length(p / r);
    float k1 = length(p / (r * r));
    return k0 * (k0 - 1.0) / k1;
}

float cylinderSDF(vec3 p, vec2 h) {
  vec2 d = abs(vec2(length(p.xz),p.y)) - h;
  return min(max(d.x,d.y),0.0) + length(max(d,0.0));
}

// arbitrary orientation
float cylinderSDF(vec3 p, vec3 a, vec3 b, float r) {
    vec3 pa = p - a;
    vec3 ba = b - a;
    float baba = dot(ba,ba);
    float paba = dot(pa,ba);
#if 0    
    float ibal = inversesqrt(baba);
    float x = length(pa-ba*paba*ibal*ibal) - r;
    float y = (abs(paba-baba*0.5)-baba*0.5)*ibal;
    return min(max(x,y),0.0) + length(max(vec2(x,y),0.0));
#else
    float x = length(pa*baba-ba*paba) - r*baba;
    float y = abs(paba-baba*0.5)-baba*0.5;
    float x2 = x*x;
    float y2 = y*y*baba;
    float d = (max(x,y)<0.0)?-min(x2,y2):(((x>0.0)?x2:0.0)+((y>0.0)?y2:0.0));
    return sign(d)*sqrt(abs(d))/baba;
#endif    
}

float capsuleSDF(vec3 p, vec3 a, vec3 b, float r) {
	vec3 pa = p-a, ba = b-a;
	float h = clamp( dot(pa,ba)/dot(ba,ba), 0.0, 1.0 );
	return length( pa - ba*h ) - r;
}

float equilateralTriangleSDF(  in vec2 p ) {
    const float k = 1.73205;//sqrt(3.0);
    p.x = abs(p.x) - 1.0;
    p.y = p.y + 1.0/k;
    if( p.x + k*p.y > 0.0 ) p = vec2( p.x - k*p.y, -k*p.x - p.y )/2.0;
    p.x += 2.0 - 2.0*clamp( (p.x+2.0)/2.0, 0.0, 1.0 );
    return -length(p)*sign(p.y);
}

float triPrismSDF( vec3 p, vec2 h ) {
    vec3 q = abs(p);
    float d1 = q.z-h.y;
    h.x *= 0.866025;
    float d2 = equilateralTriangleSDF(p.xy/h.x)*h.x;
    return length(max(vec2(d1,d2),0.0)) + min(max(d1,d2), 0.);
}

float sceneSDF(vec3 pos) {
    // Union all the shapes together so all will be rendered together

    // Boxular goop
    float res = unionOp(boxSDF(pos + vec3(5.0, 1.0 + triangleWave(u_Time / u_SlowFactor, 1.0, 1.0), 0.0), vec3(0.5, 0.5, 0.2)),
                        boxSDF(pos + vec3(5.0, triangleWave(u_Time / u_SlowFactor, 1.0, 1.0), 0.0), vec3(1.0, 1.0, 1.0)));
    res = unionOp(res,sphereSDF(pos + vec3(4.5, 0.75 + triangleWave(u_Time / u_SlowFactor, 1.0, 1.0), 1.0), 0.1)); // left eye
    res = unionOp(res,sphereSDF(pos + vec3(5.5, 0.75 + triangleWave(u_Time / u_SlowFactor, 1.0, 1.0), 1.0), 0.1)); // right eye
    res = unionOp(res, boxSDF(pos + vec3(5.2, 2.0 + triangleWave(u_Time / u_SlowFactor, 1.0, 1.0), 0.0), vec3(0.06, 1.5, 0.06))); // right leg
    res = unionOp(res, boxSDF(pos + vec3(4.8, 2.0 + triangleWave(u_Time / u_SlowFactor, 1.0, 1.0), 0.0), vec3(0.06, 1.5, 0.06))); // left leg

    // Antler goop
    res = unionOp(res, roundConeSDF(pos + vec3(-5.0, sin(u_Time / u_SlowFactor), 0.0), 1.5, 0.2, 2.2)); // head
    res = unionOp(res, roundBoxSDF(pos + vec3(-5.0, 1.8 + sin(u_Time / u_SlowFactor), 0.0), vec3(0.30, 0.5, 0.30), 0.2)); // body
    res = unionOp(res, ellipsoidSDF(pos + vec3(-5.6, 0.5 + sin(u_Time / u_SlowFactor), 1.5), vec3(0.1, 0.1, 0.05))); // left eye
    res = unionOp(res, ellipsoidSDF(pos + vec3(-4.4, 0.5 + sin(u_Time / u_SlowFactor), 1.5), vec3(0.1, 0.1, 0.05))); // right eye
    res = unionOp(res, triPrismSDF(pos + vec3(-5.0, 0.6 + sin(u_Time / u_SlowFactor), 1.3), vec2(0.15, 0.15))); // nose
    res = unionOp(res, capsuleSDF(pos + vec3(-4.3, -2.0 + sin(u_Time / u_SlowFactor), 0.0), vec3(0.9,0.,0.), vec3(0.0,0.8,0.0), .1)); // right antler
    res = unionOp(res, capsuleSDF(pos + vec3(-4.5, -2.2 + sin(u_Time / u_SlowFactor), 0.0), vec3(-0.9,0.,0.), vec3(0.0,0.8,0.0), .1)); // right antler
    res = unionOp(res, capsuleSDF(pos + vec3(-2.8, -2.5 + sin(u_Time / u_SlowFactor), 0.0), vec3(0.9,0.,0.), vec3(0.0,0.8,0.0), .1)); // right antler

    res = unionOp(res, capsuleSDF(pos + vec3(-6.55, -1.8 + sin(u_Time / u_SlowFactor), 0.0), vec3(-1.9,0.,0.), vec3(0.0,1.4,0.0), .1)); // left antler
    res = unionOp(res, capsuleSDF(pos + vec3(-6.3, -2.9 + sin(u_Time / u_SlowFactor), 0.0), vec3(0.7,0.,0.), vec3(0.0,0.5,0.0), .1)); // left antler

    res = unionOp(res, capsuleSDF(pos + vec3(-5.0, 2.4 + sin(u_Time / u_SlowFactor), 0.0), vec3(0.9,0.,0.), vec3(0.0,1.8,0.0), .1)); // arm 
    res = unionOp(res, capsuleSDF(pos + vec3(-5.0, 2.4 + sin(u_Time / u_SlowFactor), 0.0), vec3(-0.9,0.,0.), vec3(0.0,1.8,0.0), .1)); // arm

    res = unionOp(res, capsuleSDF(pos + vec3(-5.2, 5.4 + sin(u_Time / u_SlowFactor), 0.0), vec3(0.0,0.,0.), vec3(0.0,4.8,0.0), .07)); // arm
    res = unionOp(res, capsuleSDF(pos + vec3(-4.8, 5.4 + sin(u_Time / u_SlowFactor), 0.0), vec3(0.0,0.,0.), vec3(0.0,4.8,0.0), .07)); // arm 
    return res;
}

/**
 * BVH
 */ 
bool didIntersectBoundingBox(vec3 dir, vec3 origin, out float dist, vec3 min, vec3 max) { 
    dir = normalize(dir);
    float tmin = (min.x - origin.x) / dir.x; 
    float tmax = (max.x - origin.x) / dir.x; 

    if (tmin > tmax) {
      // Swap the values
      float temp = tmin;
      tmin = tmax;
      tmax = temp;
    }
    float tymin = (min.y - origin.y) / dir.y; 
    float tymax = (max.y - origin.y) / dir.y; 
    if (tymin > tymax) {
      // Swap the values
      float temp = tymin;
      tymin = tymax;
      tymax = temp;
    }
 
    if ((tmin > tymax) || (tymin > tmax)) {
      return false; 
    }
 
    if (tymin > tmin) {
      tmin = tymin; 
    }
 
    if (tymax < tmax) {
      tmax = tymax; 
    }
 
    float tzmin = (min.z - origin.z) / dir.z; 
    float tzmax = (max.z - origin.z) / dir.z; 
 
    if (tzmin > tzmax) {
      // Swap the values
      float temp = tzmin;
      tzmin = tzmax;
      tzmax = temp;
    }
 
    if ((tmin > tzmax) || (tzmin > tmax)) {
      return false; 
    }
 
    if (tzmin > tmin) {
      tmin = tzmin; 
    } 
 
    if (tzmax < tmax) {
      tmax = tzmax; 
    }
    return true; 
}

// NOTE: hard-coded values (calculated by hand on paper...)
bool antlerGoopBoundingBox(vec3 dir, vec3 origin, out float dist) {
  vec3 min = vec3(-6.5 , -7.9 + triangleWave(u_Time / u_SlowFactor, 1.0, 1.0), -10.80);
  vec3 max = vec3(10.5, 10.8 + triangleWave(u_Time / u_SlowFactor, 1.0, 1.0), 10.3);
  return didIntersectBoundingBox(dir, origin, dist, min, max);
}

bool boxGoopBoundingBox(vec3 dir, vec3 origin, out float dist) {
  vec3 center = vec3(5.0, 1.0 + triangleWave(u_Time / u_SlowFactor, 1.0, 1.0), 0.0);
  vec3 min = vec3(-10.5, -8.5 + triangleWave(u_Time / u_SlowFactor, 1.0, 1.0), -0.5);
  vec3 max = vec3(10.5, 5.0 + triangleWave(u_Time / u_SlowFactor, 1.0, 1.0), 0.5);
  return didIntersectBoundingBox(dir, origin, dist, min, max);
}

// Ray marching
float shortestDistanceToSurface(vec3 eye, vec3 marchingDirection, float start, float end) {
    // BVH
    float dist = 0.0;
    if (!antlerGoopBoundingBox(marchingDirection, eye, dist)) {
      return 100000.0;
    }
    if (!boxGoopBoundingBox(marchingDirection, eye, dist)) {
      return 100000.0;
    }
    float depth = start;
    for (int i = 0; i < MAX_MARCHING_STEPS; i++) {
        float dist = sceneSDF(eye + depth * marchingDirection);
        if (dist < EPSILON) {
			    return depth;
        }
        depth += dist;
        if (depth >= end) {
          return end;
        }
    }
    return end;
}

// Perform a gradient calculation to approximate the normal
vec3 estimateNormal(vec3 p) {
  return normalize(vec3(
      sceneSDF(vec3(p.x + EPSILON, p.y, p.z)) - sceneSDF(vec3(p.x - EPSILON, p.y, p.z)),
      sceneSDF(vec3(p.x, p.y + EPSILON, p.z)) - sceneSDF(vec3(p.x, p.y - EPSILON, p.z)),
      sceneSDF(vec3(p.x, p.y, p.z  + EPSILON)) - sceneSDF(vec3(p.x, p.y, p.z - EPSILON))
  ));
}

// Calculate the phong contribution to light intensity
vec3 phongContribForLight(vec3 k_d, vec3 k_s, float alpha, vec3 p, vec3 eye,
                          vec3 lightPos, vec3 lightIntensity) {
  vec3 N = estimateNormal(p);
  vec3 L = normalize(lightPos - p);
  vec3 V = normalize(eye - p);
  vec3 R = normalize(reflect(-L, N));
  
  float dotLN = dot(L, N);
  float dotRV = dot(R, V);
  
  if (dotLN < 0.0) {
      // Light not visible from this point on the surface
      return vec3(0.0, 0.0, 0.0);
  } 
  if (dotRV < 0.0) {
      // Light reflection in opposite direction as viewer, apply only diffuse component
      return lightIntensity * (k_d * dotLN);
  }
  return lightIntensity * (k_d * dotLN + k_s * pow(dotRV, alpha));
}

vec3 phongIllumination(vec3 k_a, vec3 k_d, vec3 k_s, float alpha, vec3 p, vec3 eye) {
  // Add ambient light so scene is not totally dark (scale for more/less)
  const vec3 ambientLight = 0.6 * vec3(1.0, 1.0, 1.0);
  vec3 color = ambientLight * k_a;
  
  // Add light 1 in front of the reindeer creature
  vec3 light1Pos = vec3(-6.0, 2.0, -6.0 );
  vec3 light1Intensity = vec3(0.5, 0.5, 0.5);
  
  color += phongContribForLight(k_d, k_s, alpha, p, eye, light1Pos, light1Intensity);

  // Add light 2 in front of the square creature
  vec3 light2Pos = vec3(2.0, 
                        2.0,
                        2.0);
  vec3 light2Intensity = vec3(0.4, 0.4, 0.4);
  
  color += phongContribForLight(k_d, k_s, alpha, p, eye, light2Pos, light2Intensity);    
  return color;
}

vec3 castRay() {
  float sx = fs_Pos.x;
  float sy = fs_Pos.y;

  float len = length(u_Ref - u_Eye);
  vec3 look = normalize(u_Ref - u_Eye);
  vec3 right = normalize(cross(look, u_Up));
  vec3 up = cross(right, look);
  float tan_fovy = tan(45.0 / 2.0); // FOV = 45 degrees
  float aspect = u_Dimensions.x / u_Dimensions.y;
  vec3 V = up * len * tan_fovy;
  vec3 H = right * len * aspect * tan_fovy;

  vec3 p = u_Ref + sx * H + sy * V;
  vec3 dir = normalize(p - u_Eye);

  return dir;
}

/**
 * Noise functions
 */
vec2 random2( vec2 p , vec2 seed) {
  return fract(sin(vec2(dot(p + seed, vec2(311.7, 127.1)), dot(p + seed, vec2(269.5, 183.3)))) * 85734.3545);
}

float worley(float x, float y, float rows, float cols) {
    float xPos = x * float(rows) / 20.0;
    float yPos = y * float(cols) / 20.0;

    float minDist = 60.0;
    vec2 minVec = vec2(0.0, 0.0);

    // Find closest point
    for (int i = -1; i < 2; i++) {
        for (int j = -1; j < 2; j++) {
            vec2 currGrid = vec2(floor(float(xPos)) + float(i), floor(float(yPos)) + float(j));
            vec2 currNoise = currGrid + random2(currGrid, vec2(2.0, 1.0));
            float currDist = distance(vec2(xPos, yPos), currNoise);
            if (currDist <= minDist) {
                minDist = currDist;
                minVec = currNoise;
            }
        }
    }
    return minDist;
}

void main() {
  // Get ray direction
  vec3 dir = castRay();

  // Ray march along ray
  vec3 colorFromScene = vec3(0, 0, 0);
  float dist = shortestDistanceToSurface(u_Eye, dir, MIN_DIST, MAX_DIST);

  // Lambert's Law for shading
  //vec3 normal = estimateNormal(vec3(fs_Pos, 1.0));
  // float diffuseTerm = dot(normalize(normal), normalize(fs_LightVec));
  // diffuseTerm = clamp(diffuseTerm, 0.0, 1.0);

  vec3 color = vec3(0.0, 0.0, 0.0);
  if (dist > MAX_DIST - EPSILON) {
    // Procedural polka dot sky
    float noise = worley(fs_Pos.x, fs_Pos.y, 100.0, 100.0);
    if (noise < 0.4) {
      // Polka dot
      out_Col = vec4(214.0 / 255.0, 201.0 / 255.0, 201.0 / 255.0, 1.0);
    } else {
      // Cloud
      out_Col = vec4(163.0 / 255.0, 186.0 / 255.0, 195.0 / 255.0, 1.0);
    }
		return;
  } else {
    // Get the closest point on the surface to the eyepoint along the view ray
    vec3 pClosest = u_Eye + dist * dir;
      
    vec3 K_a = vec3(0.2, 0.2, 0.2);
    vec3 K_d = vec3(u_Color[0] / 255.0, u_Color[1] / 255.0, u_Color[2] / 255.0);
    vec3 K_s = vec3(1.0, 1.0, 1.0);
    float shininess = 10.0;
      
    color = phongIllumination(K_a, K_d, K_s, shininess, pClosest, u_Eye);
  }

  out_Col = vec4(color, 1.0);
}
