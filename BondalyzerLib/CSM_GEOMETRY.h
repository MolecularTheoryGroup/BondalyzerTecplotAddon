#pragma once
#ifndef CSMGEOMETRY_H_
#define CSMGEOMETRY_H_

#include <vector>
#include <armadillo>


using namespace arma;

using std::vector;

const bool CSMArrow(const vec3 & Origin,
	const vec3 & Dir,
	const double & Length,
	const double & Radius,
	const double & ArrowheadLengthRatio,
	const double & ArrowheadRadiusRatio,
	vector<vec3> & Nodes,
	vector<vector<int> > & ElemList);

#endif