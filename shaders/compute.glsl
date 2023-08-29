#version 460 core

// BEGIN OF LIBRARY INCLUDE

#define PI 3.14159265
#define TAU (2*PI)
#define PHI (sqrt(5.0)*0.5 + 0.5)

// Clamp to [0,1] - this operation is free under certain circumstances.
// For further information see
// http://www.humus.name/Articles/Persson_LowLevelThinking.pdf and
// http://www.humus.name/Articles/Persson_LowlevelShaderOptimization.pdf
#define saturate(x) clamp(x, 0.0, 1.0)

// Sign function that doesn't return 0
float sgn(float x) {
	return (x<0.0)?-1.0:1.0;
}

vec2 sgn(vec2 v) {
	return vec2((v.x<0.0)?-1.0:1.0, (v.y<0.0)?-1.0:1.0);
}

float square (float x) {
	return x*x;
}

vec2 square (vec2 x) {
	return x*x;
}

vec3 square (vec3 x) {
	return x*x;
}

float lengthSqr(vec3 x) {
	return dot(x, x);
}


// Maximum/minumum elements of a vector
float vmax(vec2 v) {
	return max(v.x, v.y);
}

float vmax(vec3 v) {
	return max(max(v.x, v.y), v.z);
}

float vmax(vec4 v) {
	return max(max(v.x, v.y), max(v.z, v.w));
}

float vmin(vec2 v) {
	return min(v.x, v.y);
}

float vmin(vec3 v) {
	return min(min(v.x, v.y), v.z);
}

float vmin(vec4 v) {
	return min(min(v.x, v.y), min(v.z, v.w));
}




////////////////////////////////////////////////////////////////
//
//             PRIMITIVE DISTANCE FUNCTIONS
//
////////////////////////////////////////////////////////////////
//
// Conventions:
//
// Everything that is a distance function is called fSomething.
// The first argument is always a point in 2 or 3-space called <p>.
// Unless otherwise noted, (if the object has an intrinsic "up"
// side or direction) the y axis is "up" and the object is
// centered at the origin.
//
////////////////////////////////////////////////////////////////

// Modified version of fSphere witch takes into account besides the radius of the sphere also it's position
float fSphereModified(vec3 p, vec4 r) {
	return length(p - r.xyz) - r.w;
}

float fSphere(vec3 p, float r) {
	return length(p) - r;
}


// Plane with normal n (n is normalized) at some distance from the origin
float fPlane(vec3 p, vec3 n, float distanceFromOrigin) {
	return dot(p, n) + distanceFromOrigin;
}

// Cheap Box: distance to corners is overestimated
float fBoxCheap(vec3 p, vec3 b) { //cheap box
	return vmax(abs(p) - b);
}

// Box: correct distance to corners
float fBox(vec3 p, vec3 b) {
	vec3 d = abs(p) - b;
	return length(max(d, vec3(0))) + vmax(min(d, vec3(0)));
}

// Same as above, but in two dimensions (an endless box)
float fBox2Cheap(vec2 p, vec2 b) {
	return vmax(abs(p)-b);
}

float fBox2(vec2 p, vec2 b) {
	vec2 d = abs(p) - b;
	return length(max(d, vec2(0))) + vmax(min(d, vec2(0)));
}


// Endless "corner"
float fCorner (vec2 p) {
	return length(max(p, vec2(0))) + vmax(min(p, vec2(0)));
}

// Blobby ball object. You've probably seen it somewhere. This is not a correct distance bound, beware.
float fBlob(vec3 p) {
	p = abs(p);
	if (p.x < max(p.y, p.z)) p = p.yzx;
	if (p.x < max(p.y, p.z)) p = p.yzx;
	float b = max(max(max(
		dot(p, normalize(vec3(1.0, 1.0, 1.0))),
		dot(p.xz, normalize(vec2(PHI+1.0, 1.0)))),
		dot(p.yx, normalize(vec2(1.0, PHI)))),
		dot(p.xz, normalize(vec2(1.0, PHI))));
	float l = length(p);
	return l - 1.5 - 0.2 * (1.5 / 2.0)* cos(min(sqrt(1.01 - b / l)*(PI / 0.25), PI));
}

// Cylinder standing upright on the xz plane
float fCylinder(vec3 p, float r, float height) {
	float d = length(p.xz) - r;
	d = max(d, abs(p.y) - height);
	return d;
}

// Capsule: A Cylinder with round caps on both sides
float fCapsule(vec3 p, float r, float c) {
	return mix(length(p.xz) - r, length(vec3(p.x, abs(p.y) - c, p.z)) - r, step(c, abs(p.y)));
}

// Distance to line segment between <a> and <b>, used for fCapsule() version 2below
float fLineSegment(vec3 p, vec3 a, vec3 b) {
	vec3 ab = b - a;
	float t = clamp((dot(p - a, ab) / dot(ab, ab)), 0.0, 1.0);
	return length((ab*t + a) - p);
}

// Capsule version 2: between two end points <a> and <b> with radius r
float fCapsule(vec3 p, vec3 a, vec3 b, float r) {
	return fLineSegment(p, a, b) - r;
}

// Torus in the XZ-plane
float fTorus(vec3 p, float smallRadius, float largeRadius) {
	return length(vec2(length(p.xz) - largeRadius, p.y)) - smallRadius;
}

// A circle line. Can also be used to make a torus by subtracting the smaller radius of the torus.
float fCircle(vec3 p, float r) {
	float l = length(p.xz) - r;
	return length(vec2(p.y, l));
}

// A circular disc with no thickness (i.e. a cylinder with no height).
// Subtract some value to make a flat disc with rounded edge.
float fDisc(vec3 p, float r) {
	float l = length(p.xz) - r;
	return l < 0.0 ? abs(p.y) : length(vec2(p.y, l));
}

// Hexagonal prism, circumcircle variant
float fHexagonCircumcircle(vec3 p, vec2 h) {
	vec3 q = abs(p);
	return max(q.y - h.y, max(q.x*sqrt(3.0)*0.5 + q.z*0.5, q.z) - h.x);
	//this is mathematically equivalent to this line, but less efficient:
	//return max(q.y - h.y, max(dot(vec2(cos(PI/3), sin(PI/3)), q.zx), q.z) - h.x);
}

// Hexagonal prism, incircle variant
float fHexagonIncircle(vec3 p, vec2 h) {
	return fHexagonCircumcircle(p, vec2(h.x*sqrt(3.0)*0.5, h.y));
}

// Cone with correct distances to tip and base circle. Y is up, 0 is in the middle of the base.
float fCone(vec3 p, float radius, float height) {
	vec2 q = vec2(length(p.xz), p.y);
	vec2 tip = q - vec2(0.0, height);
	vec2 mantleDir = normalize(vec2(height, radius));
	float mantle = dot(tip, mantleDir);
	float d = max(mantle, -q.y);
	float projected = dot(tip, vec2(mantleDir.y, -mantleDir.x));

	// distance to tip
	if ((q.y > height) && (projected < 0.0)) {
		d = max(d, length(tip));
	}

	// distance to base ring
	if ((q.x > radius) && (projected > length(vec2(height, radius)))) {
		d = max(d, length(q - vec2(radius, 0.0)));
	}
	return d;
}

//
// "Generalized Distance Functions" by Akleman and Chen.
// see the Paper at https://www.viz.tamu.edu/faculty/ergun/research/implicitmodeling/papers/sm99.pdf
//
// This set of constants is used to construct a large variety of geometric primitives.
// Indices are shifted by 1 compared to the paper because we start counting at Zero.
// Some of those are slow whenever a driver decides to not unroll the loop,
// which seems to happen for fIcosahedron und fTruncatedIcosahedron on nvidia 350.12 at least.
// Specialized implementations can well be faster in all cases.
//

const vec3 GDFVectors[19] = vec3[](
	normalize(vec3(1.0, 0.0, 0.0)),
	normalize(vec3(0.0, 1.0, 0.0)),
	normalize(vec3(0.0, 0.0, 1.0)),

	normalize(vec3(1.0, 1.0, 1.0)),
	normalize(vec3(-1.0, 1.0, 1.0)),
	normalize(vec3(1.0, -1.0, 1.0)),
	normalize(vec3(1.0, 1.0, -1.0)),

	normalize(vec3(0.0, 1.0, PHI+1.0)),
	normalize(vec3(0.0, -1.0, PHI+1.0)),
	normalize(vec3(PHI+1.0, 0.0, 1.0)),
	normalize(vec3(-PHI-1.0, 0.0, 1.0)),
	normalize(vec3(1.0, PHI+1.0, 0.0)),
	normalize(vec3(-1.0, PHI+1.0, 0.0)),

	normalize(vec3(0.0, PHI, 1.0)),
	normalize(vec3(0.0, -PHI, 1.0)),
	normalize(vec3(1.0, 0.0, PHI)),
	normalize(vec3(-1.0, 0.0, PHI)),
	normalize(vec3(PHI, 1.0, 0.0)),
	normalize(vec3(-PHI, 1.0, 0.0))
);

// Version with variable exponent.
// This is slow and does not produce correct distances, but allows for bulging of objects.
float fGDF(vec3 p, float r, float e, int begin, int end) {
	float d = 0.0;
	for (int i = begin; i <= end; ++i)
		d += pow(abs(dot(p, GDFVectors[i])), e);
	return pow(d, 1.0/e) - r;
}

// Version with without exponent, creates objects with sharp edges and flat faces
float fGDF(vec3 p, float r, int begin, int end) {
	float d = 0.0;
	for (int i = begin; i <= end; ++i)
		d = max(d, abs(dot(p, GDFVectors[i])));
	return d - r;
}

// Primitives follow:

float fOctahedron(vec3 p, float r, float e) {
	return fGDF(p, r, e, 3, 6);
}

float fDodecahedron(vec3 p, float r, float e) {
	return fGDF(p, r, e, 13, 18);
}

float fIcosahedron(vec3 p, float r, float e) {
	return fGDF(p, r, e, 3, 12);
}

float fTruncatedOctahedron(vec3 p, float r, float e) {
	return fGDF(p, r, e, 0, 6);
}

float fTruncatedIcosahedron(vec3 p, float r, float e) {
	return fGDF(p, r, e, 3, 18);
}

float fOctahedron(vec3 p, float r) {
	return fGDF(p, r, 3, 6);
}

float fDodecahedron(vec3 p, float r) {
	return fGDF(p, r, 13, 18);
}

float fIcosahedron(vec3 p, float r) {
	return fGDF(p, r, 3, 12);
}

float fTruncatedOctahedron(vec3 p, float r) {
	return fGDF(p, r, 0, 6);
}

float fTruncatedIcosahedron(vec3 p, float r) {
	return fGDF(p, r, 3, 18);
}


////////////////////////////////////////////////////////////////
//
//                DOMAIN MANIPULATION OPERATORS
//
////////////////////////////////////////////////////////////////
//
// Conventions:
//
// Everything that modifies the domain is named pSomething.
//
// Many operate only on a subset of the three dimensions. For those,
// you must choose the dimensions that you want manipulated
// by supplying e.g. <p.x> or <p.zx>
//
// <inout p> is always the first argument and modified in place.
//
// Many of the operators partition space into cells. An identifier
// or cell index is returned, if possible. This return value is
// intended to be optionally used e.g. as a random seed to change
// parameters of the distance functions inside the cells.
//
// Unless stated otherwise, for cell index 0, <p> is unchanged and cells
// are centered on the origin so objects don't have to be moved to fit.
//
//
////////////////////////////////////////////////////////////////



// Rotate around a coordinate axis (i.e. in a plane perpendicular to that axis) by angle <a>.
// Read like this: R(p.xz, a) rotates "x towards z".
// This is fast if <a> is a compile-time constant and slower (but still practical) if not.
void pR(inout vec2 p, float a) {
	p = cos(a)*p + sin(a)*vec2(p.y, -p.x);
}

// Shortcut for 45-degrees rotation
void pR45(inout vec2 p) {
	p = (p + vec2(p.y, -p.x))*sqrt(0.5);
}

// Repeat space along one axis. Use like this to repeat along the x axis:
// <float cell = pMod1(p.x,5);> - using the return value is optional.
float pMod1(inout float p, float size) {
	float halfsize = size*0.5;
	float c = floor((p + halfsize)/size);
	p = mod(p + halfsize, size) - halfsize;
	return c;
}

// Same, but mirror every second cell so they match at the boundaries
float pModMirror1(inout float p, float size) {
	float halfsize = size*0.5;
	float c = floor((p + halfsize)/size);
	p = mod(p + halfsize,size) - halfsize;
	p *= mod(c, 2.0)*2.0 - 1.0;
	return c;
}

// Repeat the domain only in positive direction. Everything in the negative half-space is unchanged.
float pModSingle1(inout float p, float size) {
	float halfsize = size*0.5;
	float c = floor((p + halfsize)/size);
	if (p >= 0)
		p = mod(p + halfsize, size) - halfsize;
	return c;
}

// Repeat only a few times: from indices <start> to <stop> (similar to above, but more flexible)
float pModInterval1(inout float p, float size, float start, float stop) {
	float halfsize = size*0.5;
	float c = floor((p + halfsize)/size);
	p = mod(p+halfsize, size) - halfsize;
	if (c > stop) { //yes, this might not be the best thing numerically.
		p += size*(c - stop);
		c = stop;
	}
	if (c <start) {
		p += size*(c - start);
		c = start;
	}
	return c;
}


// Repeat around the origin by a fixed angle.
// For easier use, num of repetitions is use to specify the angle.
float pModPolar(inout vec2 p, float repetitions) {
	float angle = 2.0*PI/repetitions;
	float a = atan(p.y, p.x) + angle/2.0;
	float r = length(p);
	float c = floor(a/angle);
	a = mod(a,angle) - angle/2.0;
	p = vec2(cos(a), sin(a))*r;
	// For an odd number of repetitions, fix cell index of the cell in -x direction
	// (cell index would be e.g. -5 and 5 in the two halves of the cell):
	if (abs(c) >= (repetitions/2)) c = abs(c);
	return c;
}

// Repeat in two dimensions
vec2 pMod2(inout vec2 p, vec2 size) {
	vec2 c = floor((p + size*0.5)/size);
	p = mod(p + size*0.5,size) - size*0.5;
	return c;
}

// Same, but mirror every second cell so all boundaries match
vec2 pModMirror2(inout vec2 p, vec2 size) {
	vec2 halfsize = size*0.5;
	vec2 c = floor((p + halfsize)/size);
	p = mod(p + halfsize, size) - halfsize;
	p *= mod(c,vec2(2.0))*2.0 - vec2(1);
	return c;
}

// Same, but mirror every second cell at the diagonal as well
vec2 pModGrid2(inout vec2 p, vec2 size) {
	vec2 c = floor((p + size*0.5)/size);
	p = mod(p + size*0.5, size) - size*0.5;
	p *= mod(c,vec2(2.0))*2.0 - vec2(1.0);
	p -= size/2;
	if (p.x > p.y) p.xy = p.yx;
	return floor(c/2.0);
}

// Repeat in three dimensions
vec3 pMod3(inout vec3 p, vec3 size) {
	vec3 c = floor((p + size*0.5)/size);
	p = mod(p + size*0.5, size) - size*0.5;
	return c;
}

// Mirror at an axis-aligned plane which is at a specified distance <dist> from the origin.
float pMirror (inout float p, float dist) {
	float s = sgn(p);
	p = abs(p)-dist;
	return s;
}

// Mirror in both dimensions and at the diagonal, yielding one eighth of the space.
// translate by dist before mirroring.
vec2 pMirrorOctant (inout vec2 p, vec2 dist) {
	vec2 s = sgn(p);
	pMirror(p.x, dist.x);
	pMirror(p.y, dist.y);
	if (p.y > p.x)
		p.xy = p.yx;
	return s;
}

// Reflect space at a plane
float pReflect(inout vec3 p, vec3 planeNormal, float offset) {
	float t = dot(p, planeNormal)+offset;
	if (t < 0.0) {
		p = p - (2.0*t)*planeNormal;
	}
	return sgn(t);
}


////////////////////////////////////////////////////////////////
//
//             OBJECT COMBINATION OPERATORS
//
////////////////////////////////////////////////////////////////
//
// We usually need the following boolean operators to combine two objects:
// Union: OR(a,b)
// Intersection: AND(a,b)
// Difference: AND(a,!b)
// (a and b being the distances to the objects).
//
// The trivial implementations are min(a,b) for union, max(a,b) for intersection
// and max(a,-b) for difference. To combine objects in more interesting ways to
// produce rounded edges, chamfers, stairs, etc. instead of plain sharp edges we
// can use combination operators. It is common to use some kind of "smooth minimum"
// instead of min(), but we don't like that because it does not preserve Lipschitz
// continuity in many cases.
//
// Naming convention: since they return a distance, they are called fOpSomething.
// The different flavours usually implement all the boolean operators above
// and are called fOpUnionRound, fOpIntersectionRound, etc.
//
// The basic idea: Assume the object surfaces intersect at a right angle. The two
// distances <a> and <b> constitute a new local two-dimensional coordinate system
// with the actual intersection as the origin. In this coordinate system, we can
// evaluate any 2D distance function we want in order to shape the edge.
//
// The operators below are just those that we found useful or interesting and should
// be seen as examples. There are infinitely more possible operators.
//
// They are designed to actually produce correct distances or distance bounds, unlike
// popular "smooth minimum" operators, on the condition that the gradients of the two
// SDFs are at right angles. When they are off by more than 30 degrees or so, the
// Lipschitz condition will no longer hold (i.e. you might get artifacts). The worst
// case is parallel surfaces that are close to each other.
//
// Most have a float argument <r> to specify the radius of the feature they represent.
// This should be much smaller than the object size.
//
// Some of them have checks like "if ((-a < r) && (-b < r))" that restrict
// their influence (and computation cost) to a certain area. You might
// want to lift that restriction or enforce it. We have left it as comments
// in some cases.
//
// usage example:
//
// float fTwoBoxes(vec3 p) {
//   float box0 = fBox(p, vec3(1));
//   float box1 = fBox(p-vec3(1), vec3(1));
//   return fOpUnionChamfer(box0, box1, 0.2);
// }
//
////////////////////////////////////////////////////////////////


// The "Chamfer" flavour makes a 45-degree chamfered edge (the diagonal of a square of size <r>):
float fOpUnionChamfer(float a, float b, float r) {
	return min(min(a, b), (a - r + b)*sqrt(0.5));
}

// Intersection has to deal with what is normally the inside of the resulting object
// when using union, which we normally don't care about too much. Thus, intersection
// implementations sometimes differ from union implementations.
float fOpIntersectionChamfer(float a, float b, float r) {
	return max(max(a, b), (a + r + b)*sqrt(0.5));
}

// Difference can be built from Intersection or Union:
float fOpDifferenceChamfer (float a, float b, float r) {
	return fOpIntersectionChamfer(a, -b, r);
}

// The "Round" variant uses a quarter-circle to join the two objects smoothly:
float fOpUnionRound(float a, float b, float r) {
	vec2 u = max(vec2(r - a,r - b), vec2(0.0));
	return max(r, min (a, b)) - length(u);
}

float fOpIntersectionRound(float a, float b, float r) {
	vec2 u = max(vec2(r + a,r + b), vec2(0.0));
	return min(-r, max (a, b)) + length(u);
}

float fOpDifferenceRound (float a, float b, float r) {
	return fOpIntersectionRound(a, -b, r);
}


// The "Columns" flavour makes n-1 circular columns at a 45 degree angle:
float fOpUnionColumns(float a, float b, float r, float n) {
	if ((a < r) && (b < r)) {
		vec2 p = vec2(a, b);
		float columnradius = r*sqrt(2.0)/((n-1.0)*2.0+sqrt(2.0));
		pR45(p);
		p.x -= sqrt(2.0)/2.0*r;
		p.x += columnradius*sqrt(2.0);
		if (mod(n,2.0) == 1.0) {
			p.y += columnradius;
		}
		// At this point, we have turned 45 degrees and moved at a point on the
		// diagonal that we want to place the columns on.
		// Now, repeat the domain along this direction and place a circle.
		pMod1(p.y, columnradius*2.0);
		float result = length(p) - columnradius;
		result = min(result, p.x);
		result = min(result, a);
		return min(result, b);
	} else {
		return min(a, b);
	}
}

float fOpDifferenceColumns(float a, float b, float r, float n) {
	a = -a;
	float m = min(a, b);
	//avoid the expensive computation where not needed (produces discontinuity though)
	if ((a < r) && (b < r)) {
		vec2 p = vec2(a, b);
		float columnradius = r*sqrt(2.0)/n/2.0;
		columnradius = r*sqrt(2.0)/((n-1.0)*2.0+sqrt(2.0));

		pR45(p);
		p.y += columnradius;
		p.x -= sqrt(2.0)/2.0*r;
		p.x += -columnradius*sqrt(2.0)/2.0;

		if (mod(n,2.0) == 1.0) {
			p.y += columnradius;
		}
		pMod1(p.y,columnradius*2.0);

		float result = -length(p) + columnradius;
		result = max(result, p.x);
		result = min(result, a);
		return -min(result, b);
	} else {
		return -m;
	}
}

float fOpIntersectionColumns(float a, float b, float r, float n) {
	return fOpDifferenceColumns(a,-b,r, n);
}

// The "Stairs" flavour produces n-1 steps of a staircase:
// much less stupid version by paniq
float fOpUnionStairs(float a, float b, float r, float n) {
	float s = r/n;
	float u = b-r;
	return min(min(a,b), 0.5 * (u + a + abs ((mod (u - a + s, 2.0 * s)) - s)));
}

// We can just call Union since stairs are symmetric.
float fOpIntersectionStairs(float a, float b, float r, float n) {
	return -fOpUnionStairs(-a, -b, r, n);
}

float fOpDifferenceStairs(float a, float b, float r, float n) {
	return -fOpUnionStairs(-a, b, r, n);
}


// Similar to fOpUnionRound, but more lipschitz-y at acute angles
// (and less so at 90 degrees). Useful when fudging around too much
// by MediaMolecule, from Alex Evans' siggraph slides
float fOpUnionSoft(float a, float b, float r) {
	float e = max(r - abs(a - b), 0.0);
	return min(a, b) - e*e*0.25/r;
}


// produces a cylindical pipe that runs along the intersection.
// No objects remain, only the pipe. This is not a boolean operator.
float fOpPipe(float a, float b, float r) {
	return length(vec2(a, b)) - r;
}

// first object gets a v-shaped engraving where it intersect the second
float fOpEngrave(float a, float b, float r) {
	return max(a, (a + r - abs(b))*sqrt(0.5));
}

// first object gets a capenter-style groove cut out
float fOpGroove(float a, float b, float ra, float rb) {
	return max(a, min(a + ra, rb - abs(b)));
}

// first object gets a capenter-style tongue attached
float fOpTongue(float a, float b, float ra, float rb) {
	return min(a, max(a - ra, abs(b) - rb));
}


// END OF LIBRARY INCLUDE

layout(local_size_x = 8, local_size_y = 8, local_size_z = 1) in;
layout(rgba32f, binding = 0) uniform image2D screen;

uniform vec2 u_resolution;
uniform float u_time;
uniform vec3 u_pos;
uniform vec3 u_dir;

// Maximum number of ray marching steps
#define MAX_STEPS 100 
// Maximum distance for ray marching
#define MAX_DIST 1000.
// Minimum distance to consider a surface hit
#define SURF_DIST .01  

struct Camera {
	vec3 position;
	vec3 direction;
	vec3 look;
};

mat3 getRotationMatrix(vec3 dir)
{
    float cosYaw = cos(dir.y);
    float sinYaw = sin(dir.y);
    float cosPitch = cos(dir.x);
    float sinPitch = sin(dir.x);
    
    mat3 yawMatrix = mat3(
        cosYaw, 0, -sinYaw,
        0, 1, 0,
        sinYaw, 0, cosYaw
    );
    
    mat3 pitchMatrix = mat3(
        1, 0, 0,
        0, cosPitch, sinPitch,
        0, -sinPitch, cosPitch
    );
    
    return yawMatrix * pitchMatrix;
}

float sdCylinder( vec3 p, vec3 c )
{
  return length(p.xz-c.xy)-c.z;
}

vec2 fOpUnion(vec2 res1, vec2 res2)
{
    return (res1.x < res2.x) ? res1 : res2;
}

vec2 fOpUnionID(vec2 res1, vec2 res2) {
    return (res1.x < res2.x) ? res1 : res2;
}

// Define the size of each voxel in the grid
#define VOXEL_SIZE 2.0

// Distance function for a single voxel
float sdVoxel(vec3 p, vec3 voxelCenter) {
    vec3 d = abs(p - voxelCenter) - vec3(VOXEL_SIZE * 0.5);
    return length(max(d, 0.0));
}

vec3 getMaterial(vec3 p, float id) {
    vec3 m;
	vec3 red = vec3(1.0, 0.0, 0.0);
	vec3 green = vec3(0.0, 1.0, 0.0);
	vec3 blue = vec3(0.0, 0.0, 1.0);
	vec3 unimat = vec3(0.3, 0.1, 0.4);
	
    switch (int(id)) {
        case 1:
			m = red; break;
        case 2:
			m = vec3(0.2 + 0.4 * mod(floor(p.x) + floor(p.z), 2.0)); break;
		case 3:
			m = green; break;
		case 4:
			m = blue; break;
		case 5:
			m = unimat; break;
		default:
            m = vec3(1.0);  // Default material
            break;
    }
    return m;
}

vec2 GetDist(vec3 p) {
    
	//pMod3(p, vec3(5));
	vec3 e = p;
	pMod3(e, vec3(5));
	
    float planeDist = fPlane(p, vec3(0.0, 1.0, 0.0), 1.0);
    float planeID = 2.0;
    vec2 plane = vec2(planeDist, planeID);
    // Sphere
	// To rotate all the p vector in ny axis by some specific angle (here, u_time)
	pR(p.zy, u_time/10);
    float sphereDist = fSphere(vec3(e.x, e.y, e.z), 1.0);
    float sphereID = 5.0;
    vec2 sphere = vec2(sphereDist, sphereID);
	
	float boxDist = fBoxCheap(vec3(p.x + 1.5, p.y, p.z), vec3(1.0, 200.0, 1.0));
	float boxID = 5.0;
	vec2 box = vec2(boxDist, boxID);
	
	float wallDist = fBox2(p.xy, vec2(1, 15));
	float wallID = 4.0;
	vec2 wall = vec2(wallDist, wallID);
	
	float uniDist = fOpUnionRound(sphereDist, boxDist, 1.0);
	float uniID = 5.0;
	vec2 unionObj = vec2(uniDist, uniID);
	
	float voxDist = sdVoxel(p, vec3(2.0, 0.0, 1.0));
	float voxID = 1.0;
	vec2 vox = vec2(voxDist, voxID);

    // Select the minimum distance between the sphere and the plane
    //vec2 d = vec2(min(sphereDist, planeDist), 0.0);
	//vec2 d = vec2(sphereDist, planeDist);
    //float d = sphereDist;
	vec2 d = fOpUnionID(sphere, plane);
	d = fOpUnionID(d, box);
	//d = fOpUnionID(d, wall);
	d = fOpUnionID(d, vox);
	d = fOpUnionID(d, vec2(uniDist, uniID));
	
	/*
	vec2 d = sphere;
    if (d.x > plane.x)
        d = plane;
    if (d.x > sphere.x)
        d = sphere;
	*/

    return d;
}

vec2 RayMarch(vec3 ro, vec3 rd) {
    vec2 dO = vec2(0.0);  // Current distance along the ray (how far from the origin you are) "Object"
    vec2 dS = vec2(0.0);

	 float epsilon = 0.001;
	
    // March the ray while the marching value is less than a predefined maximum  
    for (int i = 0; i < MAX_STEPS; i++) {
		// Current position along the ray will be the ray origin plus the ray direction multiplied with the current  closest distance estimated by the distance function
        vec3 p = ro + rd * (dO.x + epsilon);  // Compute the current position along the ray
        dS = GetDist(p);  // Get the distance to the nearest surface from the current position
        dO.x += dS.x;  // Increment the distance along the ray by the distance to the closest surface
        dO.y = dS.y;

        // Break the loop if the distance along the ray exceeds the maximum or if the distance to the surface is below the threshold
        if (dO.x > MAX_DIST || dS.x < SURF_DIST)
            break;
    }

    return dO;
}

vec3 GetNormal(vec3 p) {
    vec2 d = GetDist(p);  // Distance to the nearest surface from the point 'p'
    vec2 e = vec2(0.01, 0.0);  // Small offset for gradient computation

    // Estimate the surface normal using the gradient of the distance field
    vec3 n = d.x - vec3(
        GetDist(p - e.xyy).x,
        GetDist(p - e.yxy).x,
        GetDist(p - e.yyx).x
    );

    return normalize(n);
}

float GetLight(vec3 p) {
    vec3 lightPos = vec3(0.0, 5.0, 6.0);  // Position of the light source
    lightPos.xz += vec2(sin(u_time), cos(u_time)) * 2.;  // Move the light source based on time
    vec3 l = normalize(lightPos - p);  // Direction to the light source
    vec3 n = GetNormal(p);  // Surface normal at the point 'p'

    float dif = clamp(dot(n, l), 0., 1.);  // Diffuse lighting contribution based on the dot product of the normal and light direction
    float d = RayMarch(p + n * SURF_DIST * 2., l).x;  // Check for shadows by marching a ray slightly offset from the surface
    if (d < length(lightPos - p))  // If the distance to the surface is less than the distance to the light source, it's in shadow
        dif *= .1;

    return dif;
}
float getAmbientOcclusion(vec3 p, vec3 normal) {
    float occ = 0.0;
    float weight = 1.0;
    for (int i = 0; i < 8; i++) {
        float len = 0.01 + 0.02 * float(i * i);
        float dist = GetDist(p + normal * len).x;
        occ += (len - dist) * weight;
        weight *= 0.85;
    }
    return 1.0 - clamp(0.6 * occ, 0.0, 1.0);
}

vec3 GetSunLight(vec3 p, vec3 rd, vec3 color) {
    vec3 sunDirection = vec3(10.0, 55.0, -20.0);  // Direction of the sun
	vec3 l = normalize(sunDirection - p);
    vec3 n = GetNormal(p);  // Surface normal at the point 'p'
    vec3 v = -rd;
	vec3 r = reflect(-l, n);

	vec3 specColor = vec3(0.5);
	vec3 specular = specColor * pow(clamp(dot(r, v), 0.0, 1.0), 10.0);
    vec3 dif = color * clamp(dot(l, n), 0.0, 1.0);  // Diffuse lighting contribution based on the dot product of the normal and sun direction
	vec3 ambient = color *0.05;
	
	float d = RayMarch(p + n * SURF_DIST * 2., l).x;  // Check for shadows by marching a ray slightly offset from the surface
    if (d < length(sunDirection - p))  // If the distance to the surface is less than the distance to the light source, it's in shadow
        dif *= .1;
	
    // No shadow calculation for sun light, as it illuminates the whole scene
    float occ = getAmbientOcclusion(p, n);

    return dif + ambient + specular;
	//return dif; // to use if wanted a simpler lighting model
}

void main()
{
    ivec2 pixel_coords = ivec2(gl_GlobalInvocationID.xy);

    // These three lines of code corespond to my UV coordinates
    ivec2 dims = imageSize(screen);
    float x = (float(pixel_coords.x - 0.5 * dims.x) / dims.x); // transforms to [-1.0, 1.0]
    float y = (float(pixel_coords.y - 0.5 * dims.y) / dims.y); // transforms to [-1.0, 1.0]
	
	//vec3 ro = vec3(u_pos);
	//vec3 pos = ro;
	
	//vec3 dir = vec3(u_dir.y, u_dir.x, u_dir.z);
	
	vec2 uv = vec2(x, y);
	
	vec3 col = vec3(0.0);  // "Blank" color value
	
	Camera camera;
	
	camera.position = vec3(u_pos);
	camera.direction = vec3(u_dir.y, u_dir.x, u_dir.z);
	camera.look = normalize(vec3(uv.x, uv.y, 1.0));
	
	//vec3 rd = normalize(vec3(uv.x, uv.y, 1.0));  // Ray direction
	
	mat3 rotationMatrix = getRotationMatrix(camera.direction);
	camera.look = rotationMatrix * camera.look;
	
	//float aspectRatio = dims.x / dims.y;
	
	vec2 d = RayMarch(camera.position, camera.look);  // Perform ray marching to get the distance to the nearest surface

    vec3 p = camera.position + camera.look * d.x;  // Compute the intersection point on the surface
    
    vec3 material = getMaterial(p, d.y); // Get the material for each object

    col = GetSunLight(p, camera.look, material);  // Compute the diffuse lighting contribution
    //col = dif;  // Assign the diffuse lighting value to the color
	
	vec3 bgColor = vec3(0.2, 0.4, 0.6);  // Background color (adjust as desired)
    col = mix(col, bgColor, smoothstep(0.0, 1.0, d.x / MAX_DIST));  // Blend the color with the background based on the distance
	//col += bgColor - max(1.0 * rd.y, 0.0);
	// fog
	col += mix(col, bgColor, 1.0 - exp(-0.00008 * d.x * d.x));

	
	col = pow(col, vec3(0.4545));  // Apply gamma correction
	imageStore(screen, pixel_coords, vec4(col, 1.0));
}
