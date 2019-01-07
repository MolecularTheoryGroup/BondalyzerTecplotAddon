#pragma once
#ifndef CSMGEOMETRY_H_
#define CSMGEOMETRY_H_

#include <vector>
#include <map>
#include <set>
#include <armadillo>

#define SAVE_INT_OBJECTS

using namespace arma;

using std::vector;

bool CSMArrow(vec3 const & Origin,
	vec3 const & Dir,
	double const & Length,
	double const & Radius,
	double const & ArrowheadLengthRatio,
	double const & ArrowheadRadiusRatio,
	vector<vec3> & Nodes,
	vector<vector<int> > & ElemList);

mat getPoints(mat const & Edge, int NumPts);
mat cubTrans(mat const & e1, mat const & e2, mat const & e3);
mat cubTrans_func(mat const & abc,
	vec const & s,
	vec const & t,
	vec const & u);
double cubJacobian(mat & abc, double const & s, double const & t, double const & u);
void rquad(int N, double const & k, vec & x, vec & w);
vector<vec> tetraquad(int N, mat const & Verts);
vector<vec> GetWeightsPoints(int N);

void IntegrateUsingIsosurfaces(vector<GradPath_c*> & GPs,
	int nPts,
	VolExtentIndexWeights_s & VolInfo,
	vector<FieldDataPointer_c> const & VarPtrs,
	vector<double> & IntVals,
	bool const ContainsBondPath = false);

void NewIntegrateUsingIsosurfaces(vector<GradPath_c*> const & GPs,
	int nPts,
	VolExtentIndexWeights_s & VolInfo,
	vector<FieldDataPointer_c> const & VarPtrs,
	vector<double> & IntVals,
	vec3 const * BondCPPos = nullptr);

bool GetSortedParameterEdgeMidpoints(vector<vector<int> > const & TriElems,
	vector<vec3> const & NodeList,
	vector<vec3> & SortedEdgeMidpoints,
	vector<vector<int> > & PerimeterEdges);

void GetTriElementConnectivityList(vector<vector<int> > const * ElemListPtr,
	vector<vector<int> > & ElemConnectivity,
	int NumSharedCorners);

string GetEdgeString(int ei, int ej);

Boolean_t Vec3PathResample(vector<vec3> const & OldXYZList, int NumPoints, vector<vec3> & NewXYZList);

bool ProjectedPointToTriangleIsInterior(vec3 const & P0, vec3 & TP, vec3 const & T1, vec3 const & T2, vec3 const & T3);
double PointDistanceToTriangleSquared(vec3 const & P, vec3 & ClosestPoint, vec3 const & T1, vec3 const & T2, vec3 const & T3, bool RecordPoint = true);

void TriangleEdgeMidPointSubdivide(
	vector<vec3> & nodes, 
	vector<vector<int> > & elems, 
	int TriNum,
	vector<int> & newTriNums);
void TriangleEdgeMidPointSubdivideAroundNodes(
	vector<vec3> & nodes,
	vector<vector<int> > & tris,
	vector<std::set<int> > & trisOfNode,
	vector<int> const & nodeNums,
	vector<int> & newNodeNums,
	vector<int> & newTriNums);
void UpdateSubdivideSphericalTriangulationWithConstraintNodesAndSegments(
	vector<vec3> const & inNodes, // input sphere nodes
	vector<vector<int> > const & inTris, // input sphere triangular elements in terms of node indices
	vec3 const & sphereCenter,
	double const & sphereRadius,
	vector<vector<vec3> > & constraintNodes, // input 2d vector<vector<int>>; constraint nodes organized by segments
	vector<vec3> & outNodes, // new sphere nodes
	vector<vector<int> > & outTris, // new sphere elements in terms of outNodes indices
	std::map<int, int> & outNodeConstraintNodeIndices, // indicates which constraintNode each outNode corresponds to (-1 if none or corresponds to segment)
	std::map<int, int> & outNodeConstraintSegmentIndices); // indicates which constraint segment each outNode corresponds to (-1 if no correspondance to any constraint));

void TriangulatedSphereJiggleMesh(vector<vec3> & Nodes,
	vector<std::set<int> > const & NodeConnectivity,
	vector<bool> const & NodeIsConstrained,
	vec3 const & SphereCenter,
	double const & SphereRadius);

void GetMeshNodeConnectivity(vector<vector<int> > const & Elems,
	int NumNodes,
	vector<std::set<int> > & NodeConnectivity);

#endif