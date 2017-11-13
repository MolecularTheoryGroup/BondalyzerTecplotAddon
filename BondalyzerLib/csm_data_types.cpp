
#include <cmath>
#include <vector>

#include "TECADDON.h"

#include "CSM_DATA_TYPES.h"

#include <armadillo>
using namespace arma;

double DistSqr(const vec & A, const vec & B) { 
	return sum(square(B - A)); 
}
double Distance(const vec & A, const vec & B) { 
	return norm(A - B); 
}

double VectorAngle(const vec3 & A, const vec3 & B){
	return acos(dot(A, B) / (norm(A) * norm(B)));
}

const vec3 SphericalToCartesian(const double & r, const double & theta, const double & phi){
	vec3 out;
	double sinTheta = sin(theta);
	out << sinTheta * cos(phi)
		<< sinTheta * sin(phi)
		<< cos(theta);
	return out * r;
}

const mat44		RotationMatrix(const double & Angle, vec3 Axis){
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
const vec3 Rotate(const vec3 & Point, const double & Angle, vec3 Axis){
	mat44 RotMat = RotationMatrix(Angle, Axis);
	vec4 TmpVec4 = RotMat * join_cols(Point, ones<vec>(1));
	return vec3(TmpVec4.subvec(0, 2));
}

const double TriangleArea(const vec3 & A, const vec3 & B, const vec3 & C){
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
const double TetVolume(const vector<double> & e2){
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
const double TetVolume(const vec3 & a, const vec3 & b, const vec3 & c, const vec3 & d){
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
const double HexahedronInternalPointTetVolume(const vector<vec3> & V){
	REQUIRE(V.size() == 9);

	double vol = 0.0;
	for (auto & vi : TetInds) vol += TetVolume(V[vi[0]], V[vi[1]], V[vi[2]], V[vi[3]]);

	return vol;
}

/*
 *	Calculate volume of parallelpiped using the lattice vector that defines it
 */
const double ParallepipedVolume(const vector<vec3> & LV){
	REQUIRE(LV.size() == 3);

	return  abs(dot(LV[0], cross(LV[1], LV[2])));
}

const double ParallepipedVolume(const mat33 & LV){
	return det(LV);
}

const bool ParallelpidedPointIsInternal(const mat33 & LV, const vec3 & Origin, const vec3 & Pt){
	
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