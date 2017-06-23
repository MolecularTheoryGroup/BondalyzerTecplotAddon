
#include <cmath>
#include <vector>

#include "TECADDON.h"

#include "CSM_DATA_TYPES.h"

#include <armadillo>
using namespace arma;

double DistSqr(const vec & A, const vec & B) { return sum(square(B - A)); }
double Distance(const vec & A, const vec & B) { return norm(A - B); }

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