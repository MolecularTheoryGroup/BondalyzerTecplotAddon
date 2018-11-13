#pragma once
#ifndef CSMGEOMETRY_H_
#define CSMGEOMETRY_H_

#include <vector>
#include <armadillo>

#define SAVE_INT_OBJECTS

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

const mat getPoints(const mat & Edge, const int & NumPts);
mat cubTrans(const mat & e1, const mat & e2, const mat & e3);
const mat cubTrans_func(const mat & abc,
	const vec & s,
	const vec & t,
	const vec & u);
const double cubJacobian(mat & abc, const double & s, const double & t, const double & u);
void rquad(const int & N, const double & k, vec & x, vec & w);
const vector<vec> tetraquad(const int & N, const mat & Verts);
const vector<vec> GetWeightsPoints(const int & N);

void IntegrateUsingIsosurfaces(vector<GradPath_c*> & GPs,
	const int & nPts,
	VolExtentIndexWeights_s & VolInfo,
	const vector<FieldDataPointer_c> & VarPtrs,
	vector<double> & IntVals,
	const bool ContainsBondPath = false);

const bool GetSortedParameterEdgeMidpoints(const vector<vector<int> > & TriElems,
	const vector<vec3> & NodeList,
	vector<vec3> & SortedEdgeMidpoints,
	vector<vector<int> > & PerimeterEdges);

void GetTriElementConnectivityList(const vector<vector<int> > * ElemListPtr,
	vector<vector<int> > & ElemConnectivity,
	const int & NumSharedCorners);

#endif