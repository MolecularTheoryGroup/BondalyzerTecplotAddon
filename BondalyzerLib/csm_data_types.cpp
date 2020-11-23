
#include <cmath>
#include <vector>
#include <string>

#include "TECADDON.h"

#include "CSM_DATA_TYPES.h"

#include <armadillo>
using namespace arma;

using std::string;
using std::to_string;

int ElementCoreElectronCount(int AtomicNumber) {
	if (AtomicNumber < 3)
		return 0;
	else if (AtomicNumber < 11)
		return 2;
	else if (AtomicNumber < 19)
		return 10;
	else if (AtomicNumber < 37)
		return 18;
	else if (AtomicNumber < 55)
		return 36;
	else if (AtomicNumber < 87)
		return 54;
	else
		return 86;
}

double dot2(const vec & v) 
{ return dot(v, v); }

vec const LogSpace(double const & low, double const & high, int n){
	return exp(linspace(log(low), log(high), n));
}

double DistSqr(vec const & A, vec const & B) { 
	return sum(square(B - A)); 
}
double Distance(vec const & A, vec const & B) { 
	return norm(A - B); 
}

double VectorAngle(vec3 const & A, vec3 const & B){
	return acos(dot(A, B) / (norm(A) * norm(B)));
}

vec3 const SphericalToCartesian(double const & r, double const & theta, double const & phi){
	vec3 out;
	double sinTheta = sin(theta);
	out << sinTheta * cos(phi)
		<< sinTheta * sin(phi)
		<< cos(theta);
	return out * r;
}

const mat44		RotationMatrix(double const & Angle, vec3 const & Axis){
	double L = norm(Axis);
	double LSqr = L * L;
	double L_sinAngle = L * sin(Angle), cosAngle = cos(Angle), OneMinusCosAngle = 1. - cosAngle;
	vec3 AxisSqr = square(Axis);

	mat44 Out;

	Out << (AxisSqr[0] + (AxisSqr[1] + AxisSqr[2]) * cosAngle) / LSqr <<
		(Axis[0] * Axis[1] * OneMinusCosAngle - Axis[2] * L_sinAngle) / LSqr <<
		(Axis[0] * Axis[2] * OneMinusCosAngle + Axis[1] * L_sinAngle) / LSqr <<
		0.0 << endr <<

		(Axis[0] * Axis[1] * OneMinusCosAngle + Axis[2] * L_sinAngle) / LSqr <<
		(AxisSqr[1] + (AxisSqr[0] + AxisSqr[2]) * cosAngle) / LSqr <<
		(Axis[1] * Axis[2] * OneMinusCosAngle - Axis[0] * L_sinAngle) / LSqr <<
		0.0 << endr <<

		(Axis[0] * Axis[2] * OneMinusCosAngle - Axis[1] * L_sinAngle) / LSqr <<
		(Axis[1] * Axis[2] * OneMinusCosAngle + Axis[0] * L_sinAngle) / LSqr <<
		(AxisSqr[2] + (AxisSqr[0] + AxisSqr[1]) * cosAngle) / LSqr <<
		0.0 << endr <<

		0.0 << 0.0 << 0.0 << 1.0;

// 	Out << (AxisSqr[0] + (AxisSqr[1] + AxisSqr[2]) * cos(Angle)) / LSqr <<
// 		(Axis[0] * Axis[1] * (1 - cos(Angle)) - Axis[2] * L * sin(Angle)) / LSqr <<
// 		(Axis[0] * Axis[2] * (1 - cos(Angle)) + Axis[1] * L * sin(Angle)) / LSqr <<
// 		0.0 << endr <<
// 
// 		(Axis[0] * Axis[1] * (1 - cos(Angle)) + Axis[2] * L * sin(Angle)) / LSqr <<
// 		(AxisSqr[1] + (AxisSqr[0] + AxisSqr[2]) * cos(Angle)) / LSqr <<
// 		(Axis[1] * Axis[2] * (1 - cos(Angle)) - Axis[0] * L * sin(Angle)) / LSqr <<
// 		0.0 << endr <<
// 
// 		(Axis[0] * Axis[2] * (1 - cos(Angle)) - Axis[1] * L * sin(Angle)) / LSqr <<
// 		(Axis[1] * Axis[2] * (1 - cos(Angle)) + Axis[0] * L * sin(Angle)) / LSqr <<
// 		(AxisSqr[2] + (AxisSqr[0] + AxisSqr[1]) * cos(Angle)) / LSqr <<
// 		0.0 << endr <<
// 
// 		0.0 << 0.0 << 0.0 << 1.0;

	return Out;
}

/*
 *	Rotates Point Angle radians clockwise around Axis
 *	Note that Point will be rotated around the origin, so
 *	remember to translate it to whereever it needs to be afterwards.
 */
vec3 const Rotate(vec3 const & Point, double const & Angle, vec3 const & Axis){
	mat44 RotMat = RotationMatrix(Angle, Axis);
	vec4 TmpVec4 = RotMat * join_cols(Point, ones<vec>(1));
	return vec3(TmpVec4.subvec(0, 2));
}

double const TriangleArea(vec3 const & A, vec3 const & B, vec3 const & C){
	vec3 AB = B - A, AC = C - A;
// 
	double M = norm(AB)*norm(AC);
	double Theta = acos(dot(AB, AC) / M);

	double b = 0.5 * M * sin(Theta);

	double a = 0.5 * sqrt(pow(AB[1]*AC[2] - AB[2]*AC[1], 2)
		+ pow(AB[2]*AC[0] - AB[0]*AC[2], 2)
		+ pow(AB[0]*AC[1] - AB[1]*AC[0], 2));

	return a;
}

/*
		 v2
		/|\
	   / | \
	e1/  |  \e3
	 /   |e2 \
  v1/____|____\v3
	\ e6 |    /
	 \   |   /
	e4\  |  /e5
	   \ | /
		\|/
		 v4

 * V{1-4} are the vertices and
 * e{1-6} are the edge lengths
*/

/*
 *	Calculate tetrahedron volume using squared edge lengths.
 *	Taken from    http://keisan.casio.com/exec/system/1329962711
 */
double const TetVolume(vector<double> const & e2){
	REQUIRE(e2.size() == 6);

	double v = 0.0069444444444444 * (
		  e2[0] * e2[4] * (e2[1] + e2[2] + e2[3] + e2[5] - e2[0] - e2[4])
		+ e2[1] * e2[5] * (e2[0] + e2[2] + e2[3] + e2[4] - e2[1] - e2[5])
		+ e2[2] * e2[3] * (e2[0] + e2[1] + e2[4] + e2[5] - e2[2] - e2[3])
		- e2[0] * e2[1] * e2[3]
		- e2[1] * e2[2] * e2[4]
		- e2[0] * e2[2] * e2[5]
		- e2[3] * e2[4] * e2[5]
		);

	return sqrt(v);
}

/*
 *	Calculate tetrahedron volume given vertices
 *	Taken from https://en.wikipedia.org/wiki/Tetrahedron#Volume
 */

/*
 *	Calculate total volume of all tets formed between
 *	the vertices of a hexahedron and an internal point.
 *	V is list of vertices, where 0 (i) is the internal point.
 */
double const TetVolume(vec3 const & a, vec3 const & b, vec3 const & c, vec3 const & d){
	return abs(dot(a - d, cross(b - d, c - d))) / 6.0;
}


/*
	5-------6	    e-------f
   /|      /|	   /|      /|
  / |  *0 / |	  / |  *i / |
 1-------2  |	 a-------b  |
 |  7----|--8	 |  g----|--h
 | /     | /	 | /     | /
 3-------4	     c-------d


	   5--------------6
	  /|             /|
	 / |            / |
	/  |           /  |
   /   |   *0     /   |
  1--------------2    |
  |    |         |    |
  |    7---------|----8
  |   /          |   /
  |  /           |  / 
  | /            | /
  3--------------4

 */
// 26 edges:
//	12 for the edges
//	6 for face crossing
//	8 for interior point to vertices
vector<vector<int> > EdgeInds = {
	{ 1, 2 },	// 0
	{ 1, 3 },	// 1
	{ 1, 5 },	// 2
	{ 2, 4 }, 	// 3
	{ 2, 6 },	// 4
	{ 3, 4 }, 	// 5
	{ 3, 7 },	// 6
	{ 5, 6 },	// 7
	{ 5, 7 },	// 8
	{ 4, 8 }, 	// 9
	{ 6, 8 }, 	// 10
	{ 7, 8 },	// 11

	{ 1, 4 },	// 12
	{ 1, 6 },	// 13
	{ 1, 7 },	// 14
	{ 2, 8 },	// 15
	{ 3, 8 },	// 16
	{ 5, 8 },	// 17

	{ 0, 1 },	// 18
	{ 0, 2 },	// 19
	{ 0, 3 },	// 20
	{ 0, 4 },	// 21
	{ 0, 5 },	// 22
	{ 0, 6 },	// 23
	{ 0, 7 },	// 24
	{ 0, 8 }	// 25
};
// indices of all the tets with 0 as vertex
vector<vector<int> > TetInds = {
	{ 0, 1, 2, 3 }, { 0, 2, 3, 4 },
	{ 0, 1, 3, 5 }, { 0, 3, 5, 7 },
	{ 0, 1, 2, 5 }, { 0, 2, 5, 6 },
	{ 0, 2, 4, 6 }, { 0, 4, 6, 8 },
	{ 0, 3, 4, 7 }, { 0, 4, 7, 8 },
	{ 0, 5, 6, 7 }, { 0, 6, 7, 8 }
};
double const HexahedronInternalPointTetVolume(vector<vec3> const & V){
	REQUIRE(V.size() == 9);

	double vol = 0.0;
	for (auto & vi : TetInds) vol += TetVolume(V[vi[0]], V[vi[1]], V[vi[2]], V[vi[3]]);

	return vol;
}

/*
 *	Calculate volume of parallelpiped using the lattice vector that defines it
 */
double const ParallepipedVolume(vector<vec3> const & LV){
	REQUIRE(LV.size() == 3);

	return  abs(dot(LV[0], cross(LV[1], LV[2])));
}

double const ParallepipedVolume(mat33 const & LV){
	return det(LV);
}

bool const ParallelpidedPointIsInternal(mat33 const & LV, vec3 const & Origin, vec3 const & Pt){
	
	vector<vec3> V(9);
	V[0] = Pt;
	V[1] = Origin;
	V[2] = Origin + LV.col(0);
	V[3] = Origin + LV.col(1);
	V[4] = V[2] + LV.col(1);
	V[5] = Origin + LV.col(2);
	V[6] = V[5] + LV.col(0);
	V[7] = V[5] + LV.col(1);
	V[8] = V[6] + LV.col(1);

	double PtTetVol = HexahedronInternalPointTetVolume(V);
	double FullVol = ParallepipedVolume(LV);

	return (PtTetVol <= FullVol);
}

/*
 *	Calculates the scalar triple product
 */
double const ScTP(vec3 const & a, vec3 const & b, vec3 const & c){
	return dot(a, cross(b, c));
}

/*
 *	Calculates the barycentric coordinates of a point p in a tetrahedron
 *	defined by a-b-c-d
 */
vec4 const BaryTet(vec3 const & a,
	vec3 const & b,
	vec3 const & c,
	vec3 const & d,
	vec3 const & p)
{
	vec3 vap = p - a;
	vec3 vbp = p - b;

	vec3 vab = b - a;
	vec3 vac = c - a;
	vec3 vad = d - a;

	vec3 vbc = c - b;
	vec3 vbd = d - b;
	// ScTP computes the scalar triple product
	double va6 = ScTP(vbp, vbd, vbc);
	double vb6 = ScTP(vap, vac, vad);
	double vc6 = ScTP(vap, vad, vab);
	double vd6 = ScTP(vap, vab, vac);
	double v6 = 1 / ScTP(vab, vac, vad);

	vec4 bary;
	bary << va6*v6 << vb6*v6 << vc6*v6 << vd6*v6;

	return bary;
}

string GetEdgeString(int ei, int ej) {
	return (ei < ej ? to_string(ei) + "," + to_string(ej) : to_string(ej) + "," + to_string(ei));
}