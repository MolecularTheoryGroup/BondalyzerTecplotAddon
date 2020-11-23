#include "TECADDON.h"
#if !defined (MSWIN)
#include <unistd.h>
#endif
#include "GLOBAL.h"

#include <vector>
#include <armadillo>
#include <string>
#include <map>
#include <set>
#include <unordered_set>
#include <queue>
#include <stdio.h>

#include <iostream>
#include <fstream>

#include "CSM_GEOMETRY.h"
#include "CSM_DATA_TYPES.h"
#include "CSM_GRAD_PATH.h"

#include "updateSphericalTriangulation.h"

//#define _GBADEBUG

using namespace arma;

using std::vector;

static int const ArrowRotationSteps = 10;

bool CSMArrow(vec3 const & Origin, 
	vec3 const & Dir, 
	double const & Length, 
	double const & Radius, 
	double const & ArrowheadLengthRatio,
	double const & ArrowheadRadiusRatio,
	vector<vec3> & Nodes,
	vector<vector<int> > & ElemList)
{
	int NumNodes = ArrowRotationSteps * 7 + 1; 
	// RotSteps for bottom of cylinder, then minor and major arrowhead, then 1 for center of bottom and RotSteps for tip
	// The bottom, minor, major, and tip nodes are doubled in order to produce proper shading on arrow
	int NumElems = ArrowRotationSteps * 6; // RotSteps triangles for bottom and top and RotSteps quads (two triangles) for cylinder sides and bottom of cap
	Nodes.resize(NumNodes); 
	ElemList.resize(NumElems, vector<int>(3));

	double InnerRadius = Radius * (1.0 - ArrowheadRadiusRatio);
	double RotStep = 2.0 * PI / static_cast<double>(ArrowRotationSteps);

	// First get normal vectors
	vec3 v1 = normalise(Dir);
	vec3 v2 = v1;
	v2[2] += 1;
	vec3 v3 = normalise(cross(v1, v2));

	v2 = normalise(cross(v1, v3));

	// Cylinder bottom points
	int NodeInd = 0;

	for (int l = 0; l < 2; ++l){
		for (int i = 0; i < ArrowRotationSteps; ++i){
			Nodes[NodeInd++] = Origin + Rotate(v2 * InnerRadius, static_cast<double>(i)* RotStep, v1);
		}
	}

	// Arrowhead inner points
	vec3 ArrowheadBottom = Origin + v1 * Length * (1.0 - ArrowheadLengthRatio);

	for (int l = 0; l < 2; ++l){
		for (int i = 0; i < ArrowRotationSteps; ++i){
			Nodes[NodeInd++] = ArrowheadBottom + Rotate(v2 * InnerRadius, static_cast<double>(i)* RotStep, v1);
		}
	}

	// Arrowhead outer points
	for (int l = 0; l < 2; ++l){
		for (int i = 0; i < ArrowRotationSteps; ++i){
			Nodes[NodeInd++] = ArrowheadBottom + Rotate(v2 * Radius, static_cast<double>(i)* RotStep, v1);
		}
	}

	// Arrow tip points
	for (int i = 0; i < ArrowRotationSteps; ++i){
		Nodes[NodeInd++] = Origin + v1 * Length;
	}

	//Cylinder Bottom center
	int BottomCenter = NodeInd++;
	Nodes[BottomCenter] = Origin;

	// Now make the connectivity list.
	// First the cylinder bottom
	int ElemInd = 0, ei = 0;

	for (int i = 0; i < ArrowRotationSteps; ++i){
		ElemList[ElemInd++] = { BottomCenter, ei, i < ArrowRotationSteps - 1 ? ei + 1 : 0 };
		ei++;
	}

// 	ei = 0;

	// Now the cylinder sides, with each side as a pair of triangles
	// This happens twice, first for the cylinder, then for the bottom of the cap
	for (int l = 0; l < 2; ++l){
		for (int i = 0; i < ArrowRotationSteps - 1; ++i){
			ElemList[ElemInd++] = { ei, ei + 1, ei + ArrowRotationSteps };
			ei++;
			ElemList[ElemInd++] = { ei, ei + ArrowRotationSteps, ei + ArrowRotationSteps - 1 };
		}
		ElemList[ElemInd++] = { ei, ei - ArrowRotationSteps + 1, ei + ArrowRotationSteps };
		ei++;
		ElemList[ElemInd++] = { ei - ArrowRotationSteps, ei, ei + ArrowRotationSteps - 1 };
		ei += ArrowRotationSteps;
	}

// 	ei -= ArrowRotationSteps;

// 	ei += ArrowRotationSteps;

	// Now the top of the arrowhead
	for (int i = 0; i < ArrowRotationSteps - 1; ++i){
		// 		ElemList[ElemInd++] = { Tip, ei, i < ArrowRotationSteps - 1 ? ei + 1 : ei - ArrowRotationSteps + 1 };
		ElemList[ElemInd++] = { ei, ei + 1, ei + ArrowRotationSteps};
		ei++;
	}
	ElemList[ElemInd++] = { ei, ei + 1 - ArrowRotationSteps, ei + ArrowRotationSteps };

	return true;
}

/*
 %Calculates evenly spaced points along an edge
 % Input:
 % edge - edge data given in an ordered numedgesX3 array
 % npts - number of points
 % Output:
 % pts - evenly space points
 */
mat getPoints(mat const & Edge, int NumPts){
	// initialize pts
	mat Pts(NumPts, 3);

	// calculate the diff between consecutive points on the edge
	mat EdgeDiff = Edge.tail_rows(Edge.n_rows - 1) - Edge.head_rows(Edge.n_rows - 1);

	// calculate the distance between consecutive points
	vec Dist = join_cols(mat(vector<double>({ 0. })), sqrt(square(EdgeDiff.col(0)) + square(EdgeDiff.col(1)) + square(EdgeDiff.col(2))));

	// get total length of edge
	for (int i = 1; i < Dist.n_elem; ++i)
		Dist(i) += Dist(i - 1);

	// calculate total distance from first point to each point on the edge

	// get npts evenly space distances from 0 to total distance
	vec PtDist = linspace(0, Dist(Dist.n_elem - 1), NumPts);

	// first point is first point on edge
	Pts.row(0) = Edge.row(0);
	// last point is last point on edge
	Pts.row(NumPts - 1) = Edge.row(Edge.n_rows - 1);

	// loop through remaining points
	int Ind = 0;
	for (int i = 1; i < NumPts - 1; ++i){
		// d is the distance to the point we want
		double d = PtDist(i);

		// find the closest point closer to the first point than the point we want
		for (int j = Ind; j < Dist.n_elem - 1; ++j){
			if (Dist(j) <= d && Dist(j + 1) > d){
				Ind = j;
				break;
			}
		}
		// calculate distance from closest point to point we want
		double Rem = d - Dist(Ind);

		// get dist from closest point before and closest point after point we want
		double Tot = Dist(Ind + 1) - Dist(Ind);

		// get ratio
		double Rat = Rem / Tot;

		// define the point we want
		Pts.row(i) = Edge.row(Ind) + (Edge.row(Ind + 1) - Edge.row(Ind)) * Rat;
	}

	return Pts;
}

mat cubTrans(mat const & e1, mat const & e2, mat const & e3){

	// Define the 20 nodes on the reference tet
	mat Rpt;
	Rpt << 0 << 0 << 0 << endr
		<< 1. / 3. << 0 << 0 << endr
		<< 2. / 3. << 0 << 0 << endr

		<< 1 << 0 << 0 << endr
		<< 0 << 1. / 3. << 0 << endr
		<< 0 << 2. / 3. << 0 << endr

		<< 0 << 1 << 0 << endr
		<< 0 << 0 << 1. / 3. << endr
		<< 0 << 0 << 2. / 3. << endr

		<< 0 << 0 << 1 << endr
		<< 2. / 3. << 1. / 3. << 0 << endr
		<< 1. / 3. << 2. / 3. << 0 << endr

		<< 2. / 3. << 0 << 1. / 3. << endr
		<< 1. / 3. << 0 << 2. / 3. << endr
		<< 0 << 2. / 3. << 1. / 3. << endr

		<< 0 << 1. / 3. << 2. / 3. << endr
		<< 1. / 3. << 1. / 3. << 0 << endr
		<< 1. / 3. << 0 << 1. / 3. << endr

		<< 0 << 1. / 3. << 1. / 3. << endr
		<< 1. / 3. << 1. / 3. << 1. / 3.;

	// 		mat Rpt({
	// 			{ 0, 0, 0 }, { 1. / 3., 0, 0 }, { 2. / 3., 0, 0 },
	// 			{ 1, 0, 0 }, { 0, 1. / 3., 0 }, { 0, 2. / 3., 0 },
	// 			{ 0, 1, 0 }, { 0, 0, 1. / 3. }, { 0, 0, 2. / 3. },
	// 			{ 0, 0, 1 }, { 2. / 3., 1. / 3., 0 }, { 1. / 3., 2. / 3., 0 },
	// 			{ 2. / 3., 0, 1. / 3. }, { 1. / 3., 0, 2. / 3. }, { 0, 2. / 3., 1. / 3. },
	// 			{ 0, 1. / 3., 2. / 3. }, { 1. / 3., 1. / 3., 0 }, { 1. / 3., 0, 1. / 3. },
	// 			{ 0, 1. / 3., 1. / 3. }, { 1. / 3., 1. / 3., 1. / 3. }
	// 		});

	// Get 4 equally spaced points on each edge
	mat pt1 = getPoints(e1, 4),
		pt2 = getPoints(e2, 4),
		pt3 = getPoints(e3, 4);

	// Get 3 equally space points on each edge
	// (second point is the midpoint we need to get face points)
	rowvec xx = getPoints(e1, 3).row(1),
		yy = getPoints(e2, 3).row(1),
		zz = getPoints(e3, 3).row(1);

	// Define the corresponding 20 points on the volume
	mat pt(20, 3);
	int pti = 0;
	for (int i = 0; i < 4; ++i)
		pt.row(pti++) = pt1.row(i);
	for (int i = 1; i < 4; ++i)
		pt.row(pti++) = pt2.row(i);
	for (int i = 1; i < 4; ++i)
		pt.row(pti++) = pt3.row(i);
	pt.row(pti++) = (pt1.row(3) * (2. / 3.)) + (pt2.row(3) * (1. / 3.));
	pt.row(pti++) = (pt1.row(3) * (1. / 3.)) + (pt2.row(3) * (2. / 3.));

	pt.row(pti++) = (pt1.row(3) * (2. / 3.)) + (pt3.row(3) * (1. / 3.));
	pt.row(pti++) = (pt1.row(3) * (1. / 3.)) + (pt3.row(3) * (2. / 3.));

	pt.row(pti++) = (pt2.row(3) * (2. / 3.)) + (pt3.row(3) * (1. / 3.));
	pt.row(pti++) = (pt2.row(3) * (1. / 3.)) + (pt3.row(3) * (2. / 3.));

	pt.row(pti++) = (xx * (1. / 2.)) + (yy * (1. / 2.));
	pt.row(pti++) = (xx * (1. / 2.)) + (zz * (1. / 2.));
	pt.row(pti++) = (yy * (1. / 2.)) + (zz * (1. / 2.));

	pt.row(pti++) = (pt1.row(3) * (1. / 3.)) + (pt2.row(3) * (1. / 3.)) + (pt3.row(3) * (1. / 3.));

	// Build the system matrix A
	vec C1 = Rpt.col(0), C2 = Rpt.col(1), C3 = Rpt.col(2);
	vec C12 = square(C1), C22 = square(C2), C32 = square(C3);

	mat A(Rpt.n_rows, Rpt.n_rows);
	pti = 0;

	A.col(pti++) = vec(Rpt.n_rows, fill::ones);
	A.col(pti++) = C1;
	A.col(pti++) = C2;
	A.col(pti++) = C3;
	A.col(pti++) = C12;
	A.col(pti++) = C22;
	A.col(pti++) = C32;
	A.col(pti++) = C1%C2;
	A.col(pti++) = C1%C3;
	A.col(pti++) = C2%C3;
	A.col(pti++) = pow(C1, 3);
	A.col(pti++) = pow(C2, 3);
	A.col(pti++) = pow(C3, 3);
	A.col(pti++) = C12%C2;
	A.col(pti++) = C12%C3;
	A.col(pti++) = C22%C1;
	A.col(pti++) = C22%C3;
	A.col(pti++) = C32%C1;
	A.col(pti++) = C32%C2;
	A.col(pti++) = C1%C2%C3;

	// Solve the systems to get the coefficients
	mat abc(pt.n_rows, 3);

	for (int i = 0; i < 3; ++i){
		abc.col(i) = solve(A, pt.col(i));
	}

	return abc.t();
}

/*
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %Transforms points from the reference tet to the volume.
 % Inputs:
 % s - s coordinates of the points on the reference tet
 % t - t coordinates of the points on the reference tet
 % u - u coordinates of the points on the reference tet
 % a - vector of coefficients for the (s,t,u) -> x mapping
 % b - vector of coefficients for the (s,t,u) -> y mapping
 % c - vector of coefficients for the (s,t,u) -> z mapping
 % Outputs:
 % x - x coordinates of the points on the volume
 % y - y coordinates of the points on the volume
 % z - z coordinates of the points on the volume

 % The mapping is given by
 % x=a1+a2s+a3t+a4u+a5s^2+a6t^2+a7u^2+a8st+a9su+a10tu+
 % a11s^3+a12t^3+a13u^3+a14s^2t+a15s^2u+a16t^2s+a17t^2u+
 % a18u^2s+a19u^2t+a20stu
 %
 % y and z transforms are defined the same way with b and c
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 */
mat cubTrans_func(mat const & abc,
	vec const & s,
	vec const & t,
	vec const & u)
{
	mat xyz(s.n_elem, 3);
	for (int i = 0; i < 3; ++i){
		xyz.col(i) =
			(s*abc(i, 1) + abc(i, 0))
			+ t*abc(i, 2)
			+ u*abc(i, 3)
			+ square(s)*abc(i, 4)
			+ square(t)*abc(i, 5)
			+ square(u)*abc(i, 6)
			+ s%t*abc(i, 7)
			+ s%u*abc(i, 8)
			+ t%u*abc(i, 9)
			+ pow(s, 3)*abc(i, 10)
			+ pow(t, 3)*abc(i, 11)
			+ pow(u, 3)*abc(i, 12)
			+ square(s) % t*abc(i, 13)
			+ square(s) % u*abc(i, 14)
			+ square(t) % s*abc(i, 15)
			+ square(t) % u*abc(i, 16)
			+ square(u) % s*abc(i, 17)
			+ square(u) % t*abc(i, 18)
			+ s%t%u*abc(i, 19);
	}
	return xyz;
}

/*
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % Computes the determinant of the Jacobian for the mapping.
 % Inputs:
 % s - s coordinates of the points on the reference tet
 % t - t coordinates of the points on the reference tet
 % u - u coordinates of the points on the reference tet
 % a - vector of coefficients for the (s,t,u) -> x mapping
 % b - vector of coefficients for the (s,t,u) -> y mapping
 % c - vector of coefficients for the (s,t,u) -> z mapping
 % Output:
 % J - determinant of the Jacobian at the point (s,t,u)
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 */
double cubJacobian(mat & abc, double const & s, double const & t, double const & u){
	// Get vectors of terms for the s,t,u partial derivatives of the mapping
	vector<vec> dstu({ vector<double>({ 0.,
		1.,
		0.,
		0.,
		s * 2,
		0.,
		0.,
		t,
		u,
		0.,
		(s*s) * 3,
		0.,
		0.,
		s*t * 2,
		s*u * 2,
		t*t,
		0.,
		u*u,
		0.,
		t*u }),
		vector<double>({ 0.,
		0.,
		1.,
		0.,
		0.,
		t * 2,
		0.,
		s,
		0.,
		u,
		0.,
		(t*t) * 3,
		0.,
		s*s,
		0.,
		t*s * 2,
		t*u * 2,
		0.,
		u*u,
		s*u }),
		vector<double>({ 0.,
		0.,
		0.,
		1.,
		0.,
		0.,
		u * 2,
		0.,
		s,
		t,
		0.,
		0.,
		(u*u) * 2,
		0.,
		s*s,
		0.,
		t*t,
		u*s * 2,
		u*t * 2,
		s*t }) });

	// Define Jacobian by multiplying by the coefficients
	mat33 J;
	for (unsigned int i = 0; i < 3; ++i){
		for (unsigned int j = 0; j < 3; ++j){
			J.at(i, j) = dot(abc.row(i), dstu[j]);
		}
	}

	// Get the determinant
	return det(J);
}

void rquad(int N, double const & k, vec & x, vec & w)
{
	double k1 = k + 1, k2 = k + 2;

	vec n = linspace<vec>(1, N, N);

	vec nnk = 2 * n + k;

	rowvec A = join_rows(mat(vector<double>({ k / k2 })), ((ones<vec>(N) * (k*k)) / (nnk % (nnk + 2))).t());

	n = n.tail(N - 1);
	nnk = nnk.tail(N - 1);

	double B1 = 4. * k1 / (k2*k2*(k + 3));

	vec nk = n + k, nnk2 = square(nnk);
	vec B = 4. * square(n%nk) / (square(nnk2) - nnk2);

	mat ab = join_rows(A.t(), join_cols(vec(vector<double>({ pow(2, k1) / k1, B1 })), B));

	vec s = sqrt(ab(span(1, N - 1), 1));

	mat V;
	vec X;

	eig_sym(X, V, diagmat(ab(span(0, N - 1), 0)) + diagmat(s, -1) + diagmat(s, 1));

	// Grid points
	x = (X + 1) / 2;

	// Quadrature weights
	w = pow(0.5, k1) * ab(0, 1) * square(V.row(0).t());

	return;
}

/*
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %
 % tetraquad.m - Gaussian Quadrature for a tetrahedron
 %
 % Construct Gauss points and weights for a tetrahedral domain with vertices
 % specified by the 4x3 matrix vert. Where each row contains the (x,y,z) for
 % a vertex.
 %
 % Sample usage:
 %
 % Suppose f(x,y,z)=x^2+y^2+z^2. Then let's integrate this over a regular
 % tetrahedron.
 %
 % >>vert=[1/sqrt(3) 0 0; -sqrt(3)/6,1/2,0;-sqrt(3)/6,-1/2,0;0 0 sqrt(6)/3];
 % >>[X,Y,Z,W]=tetraquad(4,vert);
 % >>F=X.^2+Y.^2+Z.^2;
 % >>Q=W'*F;
 %
 % Written by: Greg von Winckel
 % Contact: gregvw(at)math(dot)unm(dot)edu
 % http://math.unm.edu/~gregvw
 %
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 */
vector<vec> tetraquad(int N, mat const & Verts)
{

	vector<vec> q(3), w123(3);

	for (int i = 0; i < 3; ++i)
		rquad(N, 2 - i, q[i], w123[i]);

	int N2 = N*N;
	int N3 = N2*N;
	vec q1(N3), q2(N3), q3(N3);

	for (int i = 0; i < N; ++i){
		for (int j = 0; j < N; ++j){
			for (int k = 0; k < N; ++k){
				int Ind = i*N2 + j*N + k;
				q1[Ind] = q[0][j];
				q2[Ind] = q[1][k];
				q3[Ind] = q[2][i];
			}
		}
	}

	vec x = ones<vec>(N3) -q1,
		y = (ones<vec>(N3) -q2) % q1,
		z = q1 % q2 % q3;

	vec w = reshape(reshape(w123[1] * w123[0].t(), N2, 1) * w123[2].t(), N3, 1);
	// 		mat c = mat({ { 1, 0, 0, 0 }, { -1, 1, 0, 0 }, { -1, 0, 1, 0 }, { -1, 0, 0, 1 } }) * Verts;
	mat c = reshape(mat(vector<double>({ 1, 0, 0, 0, -1, 1, 0, 0, -1, 0, 1, 0, -1, 0, 0, 1 })), 4, 4).t() * Verts;
	vec W = w * fabs(det(c.rows(1, 3)));

	// % Change of coordinates 
	mat XYZ = join_rows(ones<vec>(N3), join_rows(x, join_rows(y, z))) * c;

	// 	vector<vector<double> > vXYZ(XYZ.n_rows, vector<double>(XYZ.n_cols));
	// 	vector<double> vW(W.n_rows);
	// 
	// 	for (int i = 0; i < XYZ.n_rows; ++i){
	// 		for (int j = 0; j < XYZ.n_cols; ++j){
	// 			vXYZ[i][j] = XYZ(i, j);
	// 		}
	// 	}
	// 	for (int i = 0; i < W.n_rows; ++i){
	// 		vW[i] = W(i);
	// 	}

	return vector<vec>({ XYZ.col(0), XYZ.col(1), XYZ.col(2), W });
}



vector<vec> GetWeightsPoints(int N){

	mat Verts;
	Verts << 0 << 0 << 0 << endr
		<< 1 << 0 << 0 << endr
		<< 0 << 1 << 0 << endr
		<< 0 << 0 << 1;
	// 		mat Verts = { { 0, 0, 0 }, { 1, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 } };

	return tetraquad(N, Verts);
}

/*
 *	Takes a vector of vec3 that define a polygon in 3space
 *	and triangulates them according to the two "halves" of 
 *	the polygon when split at the vertices of greatest
 *	euclidean distance.
 */
vector<vector<int> > TriangulatePolygon(vector<vec3> const & V){
	vector<vector<int> > T;
	double ChkDistSqr = -1.0, TmpDistSqr;
	int MaxDistInds[2];
	int NumVerts = V.size();

	/*
	*	First find the vertices of maximum distance
	*/
	for (int i = 0; i < NumVerts - 1; ++i){
		for (int j = i + 1; j < NumVerts; ++j){
			TmpDistSqr = DistSqr(V[i], V[j]);
			if (TmpDistSqr > ChkDistSqr){
				ChkDistSqr = TmpDistSqr;
				MaxDistInds[0] = i;
				MaxDistInds[1] = j;
			}
		}
	}

	/*
	*	Create lists of the "left" and "right" path indices
	*/
	int lSize = MaxDistInds[1] - MaxDistInds[0],
		rSize = NumVerts - MaxDistInds[1] + MaxDistInds[0];
	vector<int> L(lSize), R(rSize);
	for (int i = 0; i < lSize; ++i) L[i] = MaxDistInds[0] + i;
	for (int i = 0; i < rSize; ++i) R[rSize - i - 1] = (MaxDistInds[1] + i) % NumVerts;

#ifdef _GBADEBUG
	// save "left" and "right" paths as scatter"
	vector<vec3> vv;
	for (int i : L) vv.push_back(V[i]);
	// 		SaveVec3VecAsScatterZone(vv, to_string(Iter + 1) + "Left path", ColorIndex_t(0), { 1, 2, 3 });

	vv.clear();
	for (int i : R) vv.push_back(V[i]);
	// 		SaveVec3VecAsScatterZone(vv, to_string(Iter + 1) + "Right path", ColorIndex_t(1), { 1, 2, 3 });
	// for coloring triangles
#endif

	/*
	*	Stitch the paths to get the triangulation
	*/
	StitchPaths(L, R, V, T);

	return T;
}


/*
 *	Given a list of GP pointers, a list of indices specifying points along each GP,
 *	and a number specifying which GP to use to define a value of electron density,
 *	find points along the rest of the GPs with the same rho value.
 */

bool GetIsoRhoGPPoints(vector<GradPath_c const *> const & GPs, vector<unsigned int> & IndexList, int MinGPNum, vector<vec3> & NewVerts, double const & NewStepCutoffRatio, bool MultiStep = false){
	unsigned int NumGPs = GPs.size();
	NewVerts.resize(NumGPs);
	bool NotTerminated = true;
	/*
	*	Get rho value at the GP point
	*/
	double RhoVal = GPs[MinGPNum]->RhoAt(IndexList[MinGPNum]);

	/*
	*	Now linearly interpolate the same value on the rest of the GPs
	*/
	for (int i = 0; i < NumGPs && NotTerminated; ++i){
		if (i != MinGPNum){
			/*
			*	Step down the GP until a value lower than RhoVal
			*/
			if (MultiStep){
				while (GPs[i]->RhoAt(IndexList[i]) > RhoVal && NotTerminated){
					IndexList[i]++;
					NotTerminated = IndexList[i] < GPs[i]->GetCount() - 1;
				}
			}
			else{
				if (GPs[i]->RhoAt(IndexList[i]) > RhoVal){
					IndexList[i]++;
				}
			}
			NotTerminated = IndexList[i] < GPs[i]->GetCount() - 1;

			/*
			*	If end of GP has been reached then the end of all the GPs will be used.
			*/
			if (NotTerminated){
				/*
				*	If the interpolated point is very close to the next point on GPs[i], then the next iteration
				*	of the while loop will result in a very tiny polyhedron, so only use the interpolated point
				*	if it's considerably far from the next GP point.
				*	If it's too close to the next GP point, then just use the next GP point.
				*/
				double StepRatio = (RhoVal - GPs[i]->RhoAt(IndexList[i] - 1))
					/ (GPs[i]->RhoAt(IndexList[i]) - GPs[i]->RhoAt(IndexList[i] - 1));

				if (StepRatio < NewStepCutoffRatio){
					/*
					*	Get the interpolated position as offset from old position
					*/
					NewVerts[i] = GPs[i]->XYZAt(IndexList[i] - 1)
						+ (GPs[i]->XYZAt(IndexList[i]) - GPs[i]->XYZAt(IndexList[i] - 1))
						* StepRatio;
				}
				else{
					NewVerts[i] = GPs[i]->XYZAt(IndexList[i]);
				}
			}
		}
		else{
			NewVerts[i] = GPs[i]->XYZAt(IndexList[i]);
			IndexList[MinGPNum]++;
		}
	}

	return NotTerminated;
}

/*
 *	Given a list of GP pointers, a list of indices specifying points along each GP,
 *	and a number specifying which GP to use to define a value of electron density,
 *	find points along the rest of the GPs with the same rho value.
 *	This implementation doesn't interpolate between GP points.
 */

bool GetIsoRhoGPPoints2(vector<GradPath_c const *> const & GPs, vector<unsigned int> & IndexList, int MinGPNum, vector<vec3> & NewVerts) {
	unsigned int NumGPs = GPs.size();
	NewVerts.resize(NumGPs);
	bool NotTerminated = true;
	/*
	*	Get rho value at the GP point
	*/
	double RhoVal = GPs[MinGPNum]->RhoAt(IndexList[MinGPNum]);

	/*
	*	Now linearly interpolate the same value on the rest of the GPs
	*/
	for (int i = 0; i < NumGPs && NotTerminated; ++i) {
		if (i != MinGPNum) {

			int Ind;
			double Weight;
			vec3 Pt;
			NotTerminated = !GPs[i]->GetPointAtRhoValue(RhoVal, Pt, Ind, Weight);
			NewVerts[i] = GPs[i]->XYZAt(Ind);
			IndexList[i] = Ind;

			NotTerminated = IndexList[i] < GPs[i]->GetCount() - 1;
		}
		else {
			NewVerts[i] = GPs[i]->XYZAt(IndexList[i]);
			IndexList[MinGPNum]++;
		}
	}

	return NotTerminated;
}

/*
 *	Integrate volume and variables (for provided var pointers) over a tetrahedron specified
 *	according to vertex indices and a list of vertices.
 */
void IntTet(const vector<vec3 const*> & Vptr, 
	vector<int> const & Ind, 
	int nPts, 
	vector<FieldDataPointer_c> const & VarPtrs, 
	VolExtentIndexWeights_s & VolInfo, 
	vector<double> & IntVals){
	/*
	*	I'll do this in two ways;
	*	one that generates quadrature weights for a reference tetrahedron
	*	and then maps them onto each of these tets,
	*  and then another that directly generates quadrature weights directly
	*  for the points of these tets.
	*/

	/*
	*	Direct quadrature weights on each tet
	*/
	mat V(4, 3);
	for (int i = 0; i < Ind.size(); ++i){
		for (int j = 0; j < 3; ++j) V(i, j) = Vptr[Ind[i]]->at(j);
	}
	vector<vec> xyzw = tetraquad(nPts, V);

#ifdef _GBADEBUG
	vector<vec3> GaussQuadPts(xyzw[0].n_elem);
	for (int i = 0; i < GaussQuadPts.size(); ++i)
		GaussQuadPts[i] << xyzw[0][i] << xyzw[1][i] << xyzw[2][i];
	SaveVec3VecAsScatterZone(GaussQuadPts, "Gauss Quad Pts", Black_C, { 1, 2, 3 });
#endif


	/*
	*	Calculate integral
	*/
	for (int p = 0; p < xyzw[3].n_elem; ++p){
		vec3 Pt;
		Pt << xyzw[0][p] << xyzw[1][p] << xyzw[2][p];
		if (SetIndexAndWeightsForPoint(Pt, VolInfo)){
			for (int i = 0; i < VarPtrs.size(); ++i){
				IntVals[i] += ValByCurrentIndexAndWeightsFromRawPtr(VolInfo, VarPtrs[i]) * xyzw[3][p];
			}
			IntVals.back() += xyzw[3][p];
		}
	}
} 


/*
*	These are the indices of the two possible sets of tetrahedra that can result from the decomposition
*	of the intermediate trigonal prisms.
*	The first set results from (VI2,VI6) < (VI3,VI5) and the other from when that's not satisfied.
*/
vector<vector<vector<unsigned int> > > IntTetInds = {
	// 		{
	// 			{ 1, 2, 3, 6 },
	// 			{ 1, 2, 6, 5 },
	// 			{ 1, 5, 6, 4 }
	// 		},
	// 		{
	// 			{ 1, 2, 3, 5 },
	// 			{ 1, 5, 3, 6 },
	// 			{ 1, 5, 6, 4 }
	// 		}
	{
		{ 0, 1, 2, 5 },
		{ 0, 1, 5, 4 },
		{ 0, 4, 5, 3 }
	},
	{
		{ 0, 1, 2, 4 },
		{ 0, 4, 2, 5 },
		{ 0, 4, 5, 3 }
	}
};

/*
 *	An integration method for gradient bundles that takes advantage of the
 *	face that no gradient bundle has divergent surfaces.
 *	We step "down" the gradient paths that make up the gradient bundle according
 *	to the path structure (i.e. line segments) of the GPs. The line segment of shortest
 *	length is used to define a value of electron density, and then points with the same amount
 *	of density are interpolated on the other GPs, defining a polygon whose vertices all have the
 *	same value of charge density. The next polygon is found by taking another step down the GPs,
 *	and together they form a polyhedron. The "top" face of the polyhedron is triangulated and the
 *	triangulation is used to decompose the polyhedron into trigonal prisms which are decomposed
 *	further into tetrahedra, which are then integrated using Gauss quadrature rules.
 */
void IntegrateUsingIsosurfaces(vector<GradPath_c const *> & GPs,
	int nPts,
	VolExtentIndexWeights_s & VolInfo,
	vector<FieldDataPointer_c> const & VarPtrs,
	vector<double> & IntVals,
	bool const ContainsBondPath)
{
	REQUIRE(GPs.size() > 3);
	for (auto const * p : GPs){
		REQUIRE(p->IsMade() && p->GetCount() > 3);
	}

	if (IntVals.size() != VarPtrs.size() + 1)
		IntVals.resize(VarPtrs.size() + 1, 0);

	/*
	*	We're stepping down the GPs according to values of rho.
	*	We'll assume that they all either terminate at the same rho value or critical point.
	*	Even if they don't terminate at the same rho value, we'll assume that the termination of any GP will be at a sufficiently low rho value
	*		such that most properties of interest have converged.
	*	Also assuming that GPs go "downhill" from a GBA sphere around a nuclear CP.
	*	At end of each loop, confirm that there is at least one more point in each grad path, otherwise exit.
	*/
	bool NotTerminated = true;

	// Used in avoiding very small volume polyhedra.
	double NewStepCutoffRatio = 0.7;

	unsigned int NumGPs = GPs.size();

	/*
	*	If ContainsBondPath, then need to know which GPs correspond to the bond path, because they're
	*	degenerate until the bond point is reached.
	*	We'll use another vector of GP pointers to control which GPs are used. For GBs that do not
	*	contain a bond path, this doesn't change anything.
	*	
	*	This adds two cases to the function: a special case that happens when the GPs split after
	*	the bond point, and another typical case that happens for the remainder of the function
	*	once the full set of GPs is in use. 
	*	
	*	The special case is easy to deal with. The triangulation of the "top" polygon results in
	*	both trigonal prisms and tetrahedra, and the tetrahedra are easy to identify because they
	*	correspond to triangles in the triangulation with an edge corresponding to degenerate GPs.
	*	Triangles without a degenerate GP edge will form trigonal prisms.
	*	
	*	The second case requires no special treatment.
	*	
	*	First, get the indices of non-degenerate GPs.
	*/
	vector<GradPath_c const *> GPi = GPs;
	int BondPointNodeNum = -1;
	int DegenGPNumFull = -1, DegenGPNum;
	vector<vector<int> > T2;
	vector<bool> IsDegenerate;
	vector<unsigned int> OldIndexList;

 	if (ContainsBondPath){
 		/*
 		 *	First, get the list of unique grad paths and overwrite GPs 
 		 *	with it.
 		 */
 		IsDegenerate.resize(NumGPs, false);
 		GradPath_c const * DegenGP = nullptr;
 		for (int i = 0; i < NumGPs - 1; ++i){
 			for (int j = i + 1; j < NumGPs; ++j){
 				if (sum(GPi[i]->XYZAt(0) == GPi[j]->XYZAt(0)) == 3){
 					IsDegenerate[j] = true;
 					if (DegenGPNumFull < 0){
 						DegenGP = GPi[i];
 						DegenGPNumFull = i;
 					}
 				}
 			}
 		}
 		GPs.clear();
 		for (int i = 0; i < NumGPs; ++i){
 			if (!IsDegenerate[i]){
 				if (i == DegenGPNumFull) DegenGPNum = GPs.size();
 				GPs.push_back(GPi[i]);
 			}
 		}
 
 		/*
 		 *	Now we need to find where the degenerate GPs stop being degenerate,
 		 *	which means we need to find where the bond point occurs, but only for the
 		 *	single GP in the degenerate set that we use for the degenerate portion of
 		 *	the function (DegenCheckGPs[0]).
 		 *	We'll find the bond point by checking for where the GP makes a 90 degree turn,
 		 *	then confirm that's the last shared point among the degenerate GPs.
 		 *	This is done by finding the point where the neighboring line segments are most
 		 *	perpendicular, then checking that this point is perpendicular enough.
 		 *	
 		 *	First, step down the GP, checking for were the dot product of the vectors
 		 *	describing neighboring line segments becomes sufficiently close to 0.
 		 */
 		double MinValue = DBL_MAX, TmpValue;
 		for (int i = 2; i < DegenGP->GetCount(); ++i){
 			TmpValue = abs(dot(
 				normalise(DegenGP->XYZAt(i) - DegenGP->XYZAt(i - 1)),
 				normalise(DegenGP->XYZAt(i - 1) - DegenGP->XYZAt(i - 2))
 				));
 			if (TmpValue < MinValue){
 				MinValue = TmpValue;
 				BondPointNodeNum = i - 1;
 			}
 		}
 		if (BondPointNodeNum < 0 || MinValue > 0.1){
 			TecUtilDialogErrMsg("Failed to find perpendicular (enough) point on bond path GP");
 		}
 		/*
 		 *	With the bond point found, find the same point on the degenerate GPs
 		 *	and confirm that the next point is not degenerate.
 		 */
//  		bool BondPointFound = true;
//  		for (int i = 0; i < NumGPs && BondPointFound; ++i){
//  			if (IsDegenerate[i]){
//  				for (int j = 0; j < GPi[i]->GetCount(); ++j){
//  					if (sum(GPi[i]->XYZAt(j) == DegenGP->XYZAt(BondPointNodeNum)) == 3){
//  						BondPointFound = (sum(GPi[i]->XYZAt(j + 1) == DegenGP->XYZAt(BondPointNodeNum + 1)) == 0);
//  						break;
//  					}
//  				}
//  			}
//  		}
//  		if (!BondPointFound){
//  			TecUtilDialogErrMsg("Failed to find bond point");
//  		}
 
 		/*
 		 *	All done here.
 		 *	In the main function loop, once DegenGP has reached the bond point, we'll
 		 *	deal with the special case at the bond point and then return to the full
 		 *	list of GPs for the remainder of the function.
 		 */
 
 		NumGPs = GPs.size();
 	}



	// Maintain index numbers for each path
	vector<unsigned int> IndexList(NumGPs, 1);

	// For storing the "top" and "bottom" planes of intermediate polyhedrons
	vector<vector<vec3> > Planes(2, vector<vec3>(NumGPs));

	// Get initial plane
	for (int i = 0; i < NumGPs; ++i){
		Planes[0][i] = GPs[i]->XYZAt(0);
	}


	double ChkDistSqr, TmpDistSqr;

	/*
	 *	Because the triangulation of each plane can change as we 
	 *	move down the GB, we need to find a triangulation to use
	 *	throughout the entire process.
	 *	If we're in an open system then many of the GBs will terminate
	 *	at a threshold value of rho. In this case we can use the 
	 *	terminal plane to determine the triangulation.
	 *	If the GB terminates at a cage CP then the process of finding
	 *	a good plane to use for the triangulation becomes more complicated,
	 *	so we'll just punt and use the triangulation from the halfway plane
	 *	according to the number of points in the first GP.
	 *	
	 *	There's a special case for GBs that terminate at a point, but it's easy to
	 *	deal with. The triangulation
	 *	will result in tetrahedra that can be immediately integrated rather than
	 *	trigonal prisms.
	 */

	/*
	 * Check to see if the GB terminates at a point (i.e. cage CP) by
	 * testing for degenerate terminal points.
	 */
	double TermCheckDist = 0.01,
		AvgTermDist = 0.0;// (sum(GPs[0]->XYZAt(-1) == GPs[1]->XYZAt(-1)) == 3);
	double Denom = 0.0;
	for (int i = 0; i < GPs.size() - 1; ++i) {
		Denom += GPs.size() - (i + 1);
		for (int j = i + 1; j < GPs.size(); ++j) {
			AvgTermDist += Distance(GPs[i]->XYZAt(-1), GPs[j]->XYZAt(-1));
		}
	}
	AvgTermDist /= Denom;
	bool TermAtPoint = (AvgTermDist < TermCheckDist);

	vector<vector<int> > T;

 	if (ContainsBondPath){
 		/*
 		 *	For a bond path coincident GB, need to use the triangulation from the polygon
 		 *	immediately after the bond point passed.
 		 *	Then the triangulation used before the bond point is the mapping of the triangulation
 		 *	from after.
 		 */
 		vector<unsigned int> TmpIndices(GPi.size(), 1);
 
 		
 		vector<vec3> TmpPlane(GPi.size());
//  		if (TermAtPoint){
 			/*
 			*	Get the polygon triangulation immediately after the bond point is passed
 			*/
 			TmpIndices[DegenGPNum] = GPi[DegenGPNum]->RhoAt(BondPointNodeNum + 1);
 			if (!GetIsoRhoGPPoints(GPi, TmpIndices, DegenGPNum, TmpPlane, NewStepCutoffRatio, true)){
 				TecUtilDialogErrMsg("GP(s) terminated too close to bond point rho value");
 			}
//  		}
//  		else{
//  			/*
//  			 *	Get the polygon triangulation from the end of the GB
//  			 */
//  			for (int i = 0; i < GPi.size(); ++i){
//  				TmpPlane[i] = GPi[i]->XYZAt(-1);
//  			}
//  		}
 		T2 = TriangulatePolygon(TmpPlane);
 
 		/*
 		 * Now map the triangulation back to the degenerate GP case.
 		 * Do this by including only triangles that have an edge of degenerate nodes,
 		 * and then update those that don't, changing the degenerate indices to DegenGPNum.
 		 */
 		bool TriIsDegenerate;
 		for (auto const & t : T2){
 			bool EdgeIsDegenerate = false;
 			for (int i = 0; i < 3 && !EdgeIsDegenerate; ++i){
 				EdgeIsDegenerate = true;
 				for (int e = 0; e < 2; ++e){
 					EdgeIsDegenerate = (EdgeIsDegenerate && (IsDegenerate[t[(i + e) % 3]] || t[(i + e) % 3] == DegenGPNumFull));
 				}
 			}
 			if (!EdgeIsDegenerate){
 				T.push_back(t);
 				for (int & i : T.back()) if (IsDegenerate[i]) i = DegenGPNum;
 			}
 		}
 
 		/*
 		 *	The nondegenerate triangulation points to vertex indices of the full GP set,
 		 *	so need to update it for the reduced set.
 		 *	
 		 *	First, get the offset of each vertex index.
 		 */
 
 		vector<int> IndOffsets(GPi.size(), 0);
 		int vOffset = 0;
 		for (int i = 0; i < GPi.size(); ++i){
 			if (IsDegenerate[i]) vOffset++;
 			IndOffsets[i] = vOffset;
 		}
 		for (auto & t : T) for (int & i : t) i -= IndOffsets[i];
 	}
 	else{
		if (TermAtPoint){
			vector<unsigned int> TmpIndices(NumGPs, 1);
			TmpIndices[0] = GPs[0]->GetCount() / 2;
			if (!GetIsoRhoGPPoints(GPs, TmpIndices, 0, Planes[1], NewStepCutoffRatio)){
				TecUtilDialogErrMsg("GP(s) terminated too close to halfway value of first GP (GB terminating at cage point)");
			}
		}
		else{
			for (int i = 0; i < NumGPs; ++i){
				Planes[1][i] = GPs[i]->XYZAt(-1);
			}
		}

		/*
		*	Now triangulate Planes[1].
		*/
		T = TriangulatePolygon(Planes[1]);
 	}

#ifdef _GBADEBUG
	int Iter = 0;
	vector<vector<double> > IntValList;
#endif

	/*
	*	Main loop
	*/
	while (NotTerminated){
		vector<vec3 const*> Vptr;
		int MinGPNum = -1;

		
 		if (ContainsBondPath && T.size() != T2.size()){
			/*
			 *	Haven't reached the bond point yet, so keep track of the old indices in case
			 *	the bond point is reached.
			 */
			OldIndexList = IndexList;
 		}

		ChkDistSqr = DBL_MAX;

		/*
		*	Loop over GPs to find the GP with the shortest current line segment
		*/
		for (int i = 0; i < NumGPs; ++i){
			TmpDistSqr = DistSqr(Planes[0][i], GPs[i]->XYZAt(IndexList[i]));
			if (TmpDistSqr < ChkDistSqr){
				ChkDistSqr = TmpDistSqr;
				MinGPNum = i;
			}
		}

		/*
			*	Get corresponding isoRho points from the rest of the GPs and update the IndexList values accordingly.
			*/
		NotTerminated = GetIsoRhoGPPoints(GPs, IndexList, MinGPNum, Planes[1], NewStepCutoffRatio);

		/*
		*	If end of GP has been reached then the end of all the GPs will be used.
		*/
		if (!NotTerminated){
			for (int i = 0; i < NumGPs; ++i) Planes[1][i] = GPs[i]->XYZAt(-1);
		}

 		if (ContainsBondPath
 			&& NumGPs != GPi.size()
 			&& IndexList[DegenGPNum] >= BondPointNodeNum+1){
 			/*
 			*	The function just passed the bond point, so need to deal with
 			*	the special transition case and set things up so the rest of
 			*	the function works.
 			*
 			*	For the special case we need to switch to the post bond point triangulation
 			*	(T2) and the full set of gradient paths (GPi), and get a new set of indices
 			*	to bring the degenerate grad paths to where they need to be, since we were
 			*	working only with the degenerate grad path of lowest index up until this
 			*	iteration of the loop.
 			*	
 			*	This can result in a few iterations where zero volume tetrahedra are produced,
 			*	but that's not a big deal and not worth dealing with; just let them evaluate
 			*	to a zero volume integral.
 			*/
 			GPs = GPi;
 			NumGPs = GPs.size();
 			IndexList.resize(NumGPs, 1);
			int vi = 0;
			for (int i = 0; i < NumGPs; ++i){
				if (IsDegenerate[i]){
					IndexList[i] = OldIndexList[DegenGPNum];
				}
				else{
					IndexList[i] = OldIndexList[vi++];
				}
			}
			NotTerminated = GetIsoRhoGPPoints(GPs, IndexList, MinGPNum, Planes[1], NewStepCutoffRatio, true);
			T = T2;

			/*
			 *	Need to update Planes[0] so that the indices references in the
			 *	original triangulation, that included all the GPs, are correct.
			 */
 			vector<vec3> TmpPlane(NumGPs);
 			vi = 0;
 			for (int i = 0; i < NumGPs; ++i){
 				if (IsDegenerate[i]){
 					TmpPlane[i] = Planes[0][DegenGPNum];
 				}
 				else{
 					TmpPlane[i] = Planes[0][vi++];
 				}
 			}
 			Planes[0] = TmpPlane;
 		}


		Vptr.reserve(Planes[0].size() + Planes[1].size());
		for (auto & p : Planes) for (auto & v : p) Vptr.push_back(&v);


#ifdef _GBADEBUG
		// Save each plane as scatter
		for (int i = 0; i < 2; ++i){
			SaveVec3VecAsScatterZone(Planes[i], to_string(Iter + 1) + "Plane" + to_string(i + 1), ColorIndex_t(i % 7), { 1, 2, 3 });
		}
		int tColor = 0;
#endif

		for (auto const & t : T){
#ifdef _GBADEBUG
			// Save triangle as scatter
			vector<vec3> tv;
 			for (int i = 0; i < 3; ++i) tv.push_back(Planes[1][t[i]]);
				SaveVec3VecAsScatterZone(tv, to_string(Iter + 1) + "Tri " + to_string(tColor), ColorIndex_t(tColor % 7), { 1, 2, 3 });
				tColor++;
#endif
// 			if (ContainsBondPath && )
			/*
			*	Use the triangulation to form the smaller polyhedra (either prisms or tets).
			*	Only the vertex indices will be stored, where zero through (NumGPs-1) are the Planes[1]
			*	vertices and NumGPs through 2*NumGPs-1 are Planes[0] (so index 5 in Planes[0] would have
			*	index 5 + NumGPs in its polyhedron.
			*
			*	Each edge of the triangle will become a quadrilateral face of the resulting trigonal prism.
			*	The prism will be divided into three tetrahedra according to the vertex index-based scheme
			*	described in https://www.researchgate.net/profile/Julien_Dompierre/publication/221561839_How_to_Subdivide_Pyramids_Prisms_and_Hexahedra_into_Tetrahedra/links/0912f509c0b7294059000000/How-to-Subdivide-Pyramids-Prisms-and-Hexahedra-into-Tetrahedra.pdf?_sg%5B0%5D=d8hb2mSqOn8bMAMT9j7qeJSmt1LQMpm6e6wYw56TZlf5zUvn6hPON0p1o70xAeqVlOlgDlxK8dFWTrc4ENeghQ.2_L7UztfcTV4BIZrkyuj8t7fW8qYHgFr845e7-9IahcCYeAa-OArmNheBT5motwFWsgIMtEI_6kx5sftMTkU9Q&_sg%5B1%5D=uepkYtIkj8JbJ7neHScBUzt8DB3aVwHkOgcXVPlgizlmEQyIDIAjT_rZtzmx_5DUSy2pJLhQsM9-kHPAIu69Z4o8z3fe2VjWQBxSn6wOhois.2_L7UztfcTV4BIZrkyuj8t7fW8qYHgFr845e7-9IahcCYeAa-OArmNheBT5motwFWsgIMtEI_6kx5sftMTkU9Q&_iepl=
			*	Basically, the lowest index vertex has a new edge along the two quadrilateral faces in which
			*	it appears, then the remaining quadrilateral face is split again according to the indices of
			*	its vertices.
			*
			*	We'll start by adding a layer of indirection such that the prism's "first" vertex is always
			*	the vertex of lowest index.
			*/

			vector<unsigned int> VI(6);
			int MinCorner = INT_MAX;
			for (int i : t) MinCorner = MIN(MinCorner, i);
			for (int i = 0; i < 3; ++i){
				VI[i] = t[(MinCorner + i) % 3];
				VI[i+3] = VI[i] + NumGPs;
			}

#ifdef _GBADEBUG
			vector<vec3> hexv;
			for (auto const & i : VI){
				int PlaneNum = (i < NumGPs ? 0 : 1);
				hexv.push_back(Planes[PlaneNum][i % NumGPs]);
			}
			SaveVec3VecAsScatterZone(hexv, to_string(Iter + 1) + "Hex " + to_string(tColor + 1), ColorIndex_t(tColor % 7), { 1, 2, 3 });
#endif

			/*
			*	Which of the two sets of tetrahedra will result depends entirely on the indicices
			*	of the quadrilateral face opposite VI[0], so it's an easy code.
			*/
			int TetSetNum = int(MIN(VI[1], VI[5]) < MIN(VI[2], VI[4]));

			/*
			*	Now loop over the three tets to compute their integral
			*/
#ifdef _GBADEBUG
			// for coloring tets
			int tetColor = 0;
#endif // _DEBUG

			for (auto const & tet : IntTetInds[TetSetNum])
			{

				/*
				*	If checking for degenerate tet points, do it here and quit if point(s) are degenerate.
				*/
				bool TetDegen = false;
// 				if (ContainsBondPath || (TermAtPoint && !NotTerminated)) {
					for (int i = 0; i < 3 && !TetDegen; ++i) {
						for (int j = i + 1; j < 4 && !TetDegen; ++j) {
							if (sum(*Vptr[VI[tet[i]]] == *Vptr[VI[tet[j]]]) == 3)
								TetDegen = true;
						}
					}
// 				}

				if (TetDegen) 
					continue;
				
				vector<int> Ind(4);
				for (int i = 0; i < 4; ++i) Ind[i] = VI[tet[i]];

#ifdef _GBADEBUG
				// Save tet as scatter
				vector<vec3> tetv;
				for (auto const & i : tet){
					// 					int PlaneNum = (VI[i] < NumGPs ? 0 : 1);
					// 					tetv.push_back(Planes[PlaneNum][VI[i] % NumGPs]);
					tetv.push_back(*Vptr[VI[i]]);
				}
				SaveTetVec3VecAsFEZone(tetv, to_string(Iter + 1) + " Tet " + to_string(tColor) + " " + to_string(tetColor), ColorIndex_t(tetColor % 7), { 1, 2, 3 });
				tetColor++;
				IntValList.push_back(IntVals);
				int valSize = IntValList.size();
#endif // _DEBUG

				IntTet(Vptr, Ind, nPts, VarPtrs, VolInfo, IntVals);
			}
		}

		Planes[0] = Planes[1];

#ifdef _GBADEBUG
		Iter++;
// 		if (Iter > 20) break;
#endif
	}

	return;
}


/*
 *	Same as above, but assuming that gradient bundles have either 3 or 4 gradient
 *	paths.
 *	Used for GBA after adding the ring surfaces triangulation update code. No gradient 
 *	paths are seeded down the edges of triangular sphere elements, so there's always
 *	3 or 4 GPs.
 */
void NewIntegrateUsingIsosurfaces(vector<GradPath_c const *> const & GPsIn,
	int nPts,
	VolExtentIndexWeights_s & VolInfo,
	vector<FieldDataPointer_c> const & VarPtrs,
	vector<double> & IntVals,
	vec3 const * BondCPPos)
{
#ifdef _GBADEBUG
	std::ofstream DebugLog("C:\\Users\\Haiiro\\Documents\\Tecplot Addons\\debuglog.csv", std::ofstream::app);
	DebugLog << "\n\nNew run:\n\n";
#endif

	auto GPs = GPsIn;
	REQUIRE(GPs.size() == 3 || GPs.size() == 4);
	for (auto const * p : GPs) {
		REQUIRE(p->IsMade() && p->GetCount() > 3);
	}

	if (IntVals.size() != VarPtrs.size() + 1)
		IntVals.resize(VarPtrs.size() + 1, 0);

	/*
	*	We're stepping down the GPs according to values of rho.
	*	We'll assume that they all either terminate at the same rho value or critical point.
	*	Even if they don't terminate at the same rho value, we'll assume that the termination of any GP will be at a sufficiently low rho value
	*		such that most properties of interest have converged.
	*	Also assuming that GPs go "downhill" from a GBA sphere around a nuclear CP.
	*	At end of each loop, confirm that there is at least one more point in each grad path, otherwise exit.
	*/
	bool NotTerminated = true;

	// Used in avoiding very small volume polyhedra.
	double NewStepCutoffRatio = 0.7;

	unsigned int NumGPs = GPs.size();

	/*
	*	If ContainsBondPath, then need to know which GPs correspond to the bond path, because they're
	*	degenerate until the bond point is reached.
	*	We'll use another vector of GP pointers to control which GPs are used. For GBs that do not
	*	contain a bond path, this doesn't change anything.
	*
	*	This adds two cases to the function: a special case that happens when the GPs split after
	*	the bond point, and another typical case that happens for the remainder of the function
	*	once the full set of GPs is in use.
	*
	*	The special case is easy to deal with. The triangulation of the "top" polygon results in
	*	both trigonal prisms and tetrahedra, and the tetrahedra are easy to identify because they
	*	correspond to triangles in the triangulation with an edge corresponding to degenerate GPs.
	*	Triangles without a degenerate GP edge will form trigonal prisms.
	*
	*	The second case requires no special treatment.
	*
	*	First, get the indices of non-degenerate GPs.
	*/
	vector<GradPath_c const*> GPi = GPs;
	int BondPointNodeNum = -1;
	int DegenGPNumFull = -1, DegenGPNum = 0;
	vector<vector<int> > T2;
	vector<bool> IsDegenerate;
	vector<unsigned int> OldIndexList;

	bool ContainsBondPath = (BondCPPos != nullptr);

	if (ContainsBondPath) {
		// Find the pair of GPs that include the bond path by checking for
		// which GPs contain the bond point.
		vector<vec3> ClosestPtsToBondPt(NumGPs);
		vector<int> ClosestPtNumsToBondPt(NumGPs, -1);
		IsDegenerate.resize(NumGPs, false);
		GradPath_c const * DegenGP = nullptr;
		for (int i = 0; i < NumGPs - 1; ++i) {
			if (ClosestPtNumsToBondPt[i] < 0)
				ClosestPtsToBondPt[i] = GPs[i]->ClosestPoint(*BondCPPos, ClosestPtNumsToBondPt[i]);
			if (approx_equal(ClosestPtsToBondPt[i], *BondCPPos, "absdiff", 0.01)) {
				for (int j = i + 1; j < NumGPs; ++j) {
					if (ClosestPtNumsToBondPt[j] < 0)
						ClosestPtsToBondPt[j] = GPs[j]->ClosestPoint(*BondCPPos, ClosestPtNumsToBondPt[j]);
					if (approx_equal(ClosestPtsToBondPt[i], *BondCPPos, "absdiff", 0.01)){
						IsDegenerate[j] = true;
						if (DegenGPNumFull < 0) {
							DegenGP = GPi[i];
							DegenGPNumFull = i;
						}
						break;
					}
				}
				break;
			}
		}
		GPs.clear();
		for (int i = 0; i < NumGPs; ++i) {
			if (!IsDegenerate[i]) {
				if (i == DegenGPNumFull) DegenGPNum = GPs.size();
				GPs.push_back(GPi[i]);
			}
		}
		BondPointNodeNum = ClosestPtNumsToBondPt[DegenGPNumFull];

		NumGPs = GPs.size();
	}



	// Maintain index numbers for each path
	vector<unsigned int> IndexList(NumGPs, 1);

	// For storing the "top" and "bottom" planes of intermediate polyhedrons
	vector<vector<vec3> > Planes(2, vector<vec3>(NumGPs));

	// Get initial plane
	for (int i = 0; i < NumGPs; ++i) {
		Planes[0][i] = GPs[i]->XYZAt(0);
	}


	double ChkDistSqr, TmpDistSqr;

	/*
	 *	Because the triangulation of each plane can change as we
	 *	move down the GB, we need to find a triangulation to use
	 *	throughout the entire process.
	 *	If we're in an open system then many of the GBs will terminate
	 *	at a threshold value of rho. In this case we can use the
	 *	terminal plane to determine the triangulation.
	 *	If the GB terminates at a cage CP then the process of finding
	 *	a good plane to use for the triangulation becomes more complicated,
	 *	so we'll just punt and use the triangulation from the halfway plane
	 *	according to the number of points in the first GP.
	 *
	 *	There's a special case for GBs that terminate at a point, but it's easy to
	 *	deal with. The triangulation
	 *	will result in tetrahedra that can be immediately integrated rather than
	 *	trigonal prisms.
	 */

	 /*
	  * Check to see if the GB terminates at a point (i.e. cage CP) by
	  * testing for degenerate terminal points.
	  */
	double TermCheckDist = 0.01,
		AvgTermDist = 0.0;// (sum(GPs[0]->XYZAt(-1) == GPs[1]->XYZAt(-1)) == 3);
	double Denom = 0.0;
	for (int i = 0; i < GPs.size() - 1; ++i) {
		Denom += GPs.size() - (i + 1);
		for (int j = i + 1; j < GPs.size(); ++j) {
			AvgTermDist += Distance(GPs[i]->XYZAt(-1), GPs[j]->XYZAt(-1));
		}
	}
	AvgTermDist /= Denom;
	bool TermAtPoint = (AvgTermDist < TermCheckDist);

	vector<vector<int> > T;

	int MinGPNum = DegenGPNum;

	if (ContainsBondPath) {
		/*
		 *	For a bond path coincident GB, need to use the triangulation from the polygon
		 *	immediately after the bond point passed.
		 *	Then the triangulation used before the bond point is the mapping of the triangulation
		 *	from after.
		 */
		vector<unsigned int> TmpIndices(GPi.size(), 1);


		vector<vec3> TmpPlane(GPi.size());
		//  		if (TermAtPoint){
					/*
					*	Get the polygon triangulation immediately after the bond point is passed
					*/
		TmpIndices[DegenGPNum] = GPi[DegenGPNum]->RhoAt(BondPointNodeNum + 1);
// 		if (!GetIsoRhoGPPoints(GPi, TmpIndices, DegenGPNum, TmpPlane, NewStepCutoffRatio, true)) {
		if (!GetIsoRhoGPPoints2(GPi, TmpIndices, DegenGPNum, TmpPlane)) {
			TecUtilDialogErrMsg("GP(s) terminated too close to bond point rho value");
		}

		// Connect the minimum distance nodes to form the triangulation of the iso-plane
		// after the bond point.

		if (DistSqr(TmpPlane[0], TmpPlane[2]) < DistSqr(TmpPlane[1], TmpPlane[3])){
			T2 = {
				{ 0, 1, 2},
				{ 0, 2, 3}
			};
		}
		else{
			T2 = {
				{ 0, 1, 3},
				{ 1, 2, 3}
			};
		}

// 		T2 = TriangulatePolygon(TmpPlane);

		/*
		 * Now map the triangulation back to the degenerate GP case.
		 * Do this by including only triangles that have an edge of degenerate nodes,
		 * and then update those that don't, changing the degenerate indices to DegenGPNum.
		 */
		bool TriIsDegenerate;
		for (auto const & t : T2) {
			bool EdgeIsDegenerate = false;
			for (int i = 0; i < 3 && !EdgeIsDegenerate; ++i) {
				EdgeIsDegenerate = true;
				for (int e = 0; e < 2; ++e) {
					EdgeIsDegenerate = (EdgeIsDegenerate && (IsDegenerate[t[(i + e) % 3]] || t[(i + e) % 3] == DegenGPNumFull));
				}
			}
			if (!EdgeIsDegenerate) {
				T.push_back(t);
				for (int & i : T.back()) if (IsDegenerate[i]) i = DegenGPNum;
			}
		}

		/*
		 *	The nondegenerate triangulation points to vertex indices of the full GP set,
		 *	so need to update it for the reduced set.
		 *
		 *	First, get the offset of each vertex index.
		 */

		vector<int> IndOffsets(GPi.size(), 0);
		int vOffset = 0;
		for (int i = 0; i < GPi.size(); ++i) {
			if (IsDegenerate[i]) vOffset++;
			IndOffsets[i] = vOffset;
		}
		for (auto & t : T) for (int & i : t) i -= IndOffsets[i];
	}
	else {
		if (TermAtPoint) {
			vector<unsigned int> TmpIndices(NumGPs, 1);
			TmpIndices[0] = GPs[0]->GetCount() / 2;
// 			if (!GetIsoRhoGPPoints1(GPs, TmpIndices, 0, Planes[1], NewStepCutoffRatio)) {
			if (!GetIsoRhoGPPoints2(GPs, TmpIndices, 0, Planes[1])) {
				TecUtilDialogErrMsg("GP(s) terminated too close to halfway value of first GP (GB terminating at cage point)");
			}
		}
		else {
			for (int i = 0; i < NumGPs; ++i) {
				Planes[1][i] = GPs[i]->XYZAt(-1);
			}
		}

		/*
		*	Now triangulate Planes[1].
		*/
// 		T = TriangulatePolygon(Planes[1]);
		T = { {0,1,2} };
	}

#ifdef _GBADEBUG
	int Iter = 0;
	vector<vector<double> > IntValList;
#endif

	/*
	*	Main loop
	*/
	while (NotTerminated) {
		vector<vec3 const*> Vptr;


		if (ContainsBondPath && T.size() != T2.size()) {
			/*
			 *	Haven't reached the bond point yet, so keep track of the old indices in case
			 *	the bond point is reached.
			 */
			OldIndexList = IndexList;
		}

		ChkDistSqr = DBL_MAX;

// 		/*
// 		*	Loop over GPs to find the GP with the shortest current line segment
// 		*/
// 		for (int i = 0; i < NumGPs; ++i) {
// 
// 			TmpDistSqr = DistSqr(Planes[0][i], GPs[i]->XYZAt(IndexList[i]));
// 			if (TmpDistSqr < ChkDistSqr) {
// 				ChkDistSqr = TmpDistSqr;
// 				MinGPNum = i;
// 			}
// 		}

		/*
			*	Get corresponding isoRho points from the rest of the GPs and update the IndexList values accordingly.
			*/
// 		NotTerminated = GetIsoRhoGPPoints(GPs, IndexList, MinGPNum, Planes[1], NewStepCutoffRatio);
		int MinGPInd = IndexList[MinGPNum];
		NotTerminated = GetIsoRhoGPPoints2(GPs, IndexList, MinGPNum, Planes[1]);
		if (IndexList[MinGPNum] == MinGPInd && MinGPInd < GPs[MinGPNum]->GetCount() - 1) {
			IndexList[MinGPNum]++;
			Planes[1][MinGPNum] = GPs[MinGPNum]->XYZAt(IndexList[MinGPNum]);
		}

		/*
		*	If end of GP has been reached then the end of all the GPs will be used.
		*/
		if (!NotTerminated) {
			for (int i = 0; i < NumGPs; ++i) Planes[1][i] = GPs[i]->XYZAt(-1);
		}

		if (ContainsBondPath
			&& NumGPs != GPi.size()
			&& IndexList[DegenGPNum] >= BondPointNodeNum + 1) {
			/*
			*	The function just passed the bond point, so need to deal with
			*	the special transition case and set things up so the rest of
			*	the function works.
			*
			*	For the special case we need to switch to the post bond point triangulation
			*	(T2) and the full set of gradient paths (GPi), and get a new set of indices
			*	to bring the degenerate grad paths to where they need to be, since we were
			*	working only with the degenerate grad path of lowest index up until this
			*	iteration of the loop.
			*
			*	This can result in a few iterations where zero volume tetrahedra are produced,
			*	but that's not a big deal and not worth dealing with; just let them evaluate
			*	to a zero volume integral.
			*/
			GPs = GPi;
			NumGPs = GPs.size();
			IndexList.resize(NumGPs, 1);
			int vi = 0;
			for (int i = 0; i < NumGPs; ++i) {
				if (IsDegenerate[i]) {
					IndexList[i] = OldIndexList[DegenGPNum];
				}
				else {
					IndexList[i] = OldIndexList[vi++];
				}
			}
// 			NotTerminated = GetIsoRhoGPPoints(GPs, IndexList, MinGPNum, Planes[1], NewStepCutoffRatio, true);
			NotTerminated = GetIsoRhoGPPoints2(GPs, IndexList, MinGPNum, Planes[1]);
			T = T2;

			/*
			 *	Need to update Planes[0] so that the indices references in the
			 *	original triangulation, that included all the GPs, are correct.
			 */
			vector<vec3> TmpPlane(NumGPs);
			vi = 0;
			for (int i = 0; i < NumGPs; ++i) {
				if (IsDegenerate[i]) {
					TmpPlane[i] = Planes[0][DegenGPNum];
				}
				else {
					TmpPlane[i] = Planes[0][vi++];
				}
			}
			Planes[0] = TmpPlane;
		}


		Vptr.reserve(Planes[0].size() + Planes[1].size());
		for (auto & p : Planes) for (auto & v : p) Vptr.push_back(&v);


#ifdef _GBADEBUG
		// Save each plane as scatter
		for (int i = 0; i < 2; ++i) {
			 			SaveVec3VecAsScatterZone(Planes[i], to_string(Iter + 1) + "Plane" + to_string(i + 1), ColorIndex_t(i % 7), { 1, 2, 3 });
		}
		int tColor = 0;
#endif

		for (auto const & t : T) {
#ifdef _GBADEBUG
			// Save triangle as scatter
			vector<vec3> tv;
			 			for (int i = 0; i < 3; ++i) tv.push_back(Planes[1][t[i]]);
						 			SaveVec3VecAsScatterZone(tv, to_string(Iter + 1) + "Tri " + to_string(tColor), ColorIndex_t(tColor % 7), { 1, 2, 3 });
			tColor++;
#endif
			// 			if (ContainsBondPath && )
						/*
						*	Use the triangulation to form the smaller polyhedra (either prisms or tets).
						*	Only the vertex indices will be stored, where zero through (NumGPs-1) are the Planes[1]
						*	vertices and NumGPs through 2*NumGPs-1 are Planes[0] (so index 5 in Planes[0] would have
						*	index 5 + NumGPs in its polyhedron.
						*
						*	Each edge of the triangle will become a quadrilateral face of the resulting trigonal prism.
						*	The prism will be divided into three tetrahedra according to the vertex index-based scheme
						*	described in https://www.researchgate.net/profile/Julien_Dompierre/publication/221561839_How_to_Subdivide_Pyramids_Prisms_and_Hexahedra_into_Tetrahedra/links/0912f509c0b7294059000000/How-to-Subdivide-Pyramids-Prisms-and-Hexahedra-into-Tetrahedra.pdf?_sg%5B0%5D=d8hb2mSqOn8bMAMT9j7qeJSmt1LQMpm6e6wYw56TZlf5zUvn6hPON0p1o70xAeqVlOlgDlxK8dFWTrc4ENeghQ.2_L7UztfcTV4BIZrkyuj8t7fW8qYHgFr845e7-9IahcCYeAa-OArmNheBT5motwFWsgIMtEI_6kx5sftMTkU9Q&_sg%5B1%5D=uepkYtIkj8JbJ7neHScBUzt8DB3aVwHkOgcXVPlgizlmEQyIDIAjT_rZtzmx_5DUSy2pJLhQsM9-kHPAIu69Z4o8z3fe2VjWQBxSn6wOhois.2_L7UztfcTV4BIZrkyuj8t7fW8qYHgFr845e7-9IahcCYeAa-OArmNheBT5motwFWsgIMtEI_6kx5sftMTkU9Q&_iepl=
						*	Basically, the lowest index vertex has a new edge along the two quadrilateral faces in which
						*	it appears, then the remaining quadrilateral face is split again according to the indices of
						*	its vertices.
						*
						*	We'll start by adding a layer of indirection such that the prism's "first" vertex is always
						*	the vertex of lowest index.
						*/

			vector<unsigned int> VI(6);
			int MinCorner = INT_MAX;
			for (int i : t) MinCorner = MIN(MinCorner, i);
			for (int i = 0; i < 3; ++i) {
				VI[i] = t[(MinCorner + i) % 3];
				VI[i + 3] = VI[i] + NumGPs;
			}

#ifdef _GBADEBUG
			vector<vec3> hexv;
			for (auto const & i : VI) {
				int PlaneNum = (i < NumGPs ? 0 : 1);
				hexv.push_back(Planes[PlaneNum][i % NumGPs]);
			}
			SaveVec3VecAsScatterZone(hexv, to_string(Iter + 1) + "Hex " + to_string(tColor + 1), ColorIndex_t(tColor % 7), { 1, 2, 3 });
#endif

			/*
			*	Which of the two sets of tetrahedra will result depends entirely on the indicices
			*	of the quadrilateral face opposite VI[0], so it's an easy code.
			*/
			int TetSetNum = int(MIN(VI[1], VI[5]) < MIN(VI[2], VI[4]));

			/*
			*	Now loop over the three tets to compute their integral
			*/
#ifdef _GBADEBUG
			// for coloring tets
			int tetColor = 0;
#endif // _DEBUG

			for (auto const & tet : IntTetInds[TetSetNum])
			{

				/*
				*	If checking for degenerate tet points, do it here and quit if point(s) are degenerate.
				*/
				bool TetDegen = false;
				// 				if (ContainsBondPath || (TermAtPoint && !NotTerminated)) {
				for (int i = 0; i < 3 && !TetDegen; ++i) {
					for (int j = i + 1; j < 4 && !TetDegen; ++j) {
						if (sum(*Vptr[VI[tet[i]]] == *Vptr[VI[tet[j]]]) == 3)
							TetDegen = true;
					}
				}
				// 				}

				if (TetDegen)
					continue;

				vector<int> Ind(4);
				for (int i = 0; i < 4; ++i) Ind[i] = VI[tet[i]];

#ifdef _GBADEBUG
				// Save tet as scatter
				vector<vec3> tetv;
				for (auto const & i : tet) {
					// 					int PlaneNum = (VI[i] < NumGPs ? 0 : 1);
					// 					tetv.push_back(Planes[PlaneNum][VI[i] % NumGPs]);
					tetv.push_back(*Vptr[VI[i]]);
				}
#ifdef _GBADEBUG
				DebugLog << "Iter " << Iter + 1 << ",Tet " << tColor << " " << tetColor;
#endif // _DEBUG

				SaveTetVec3VecAsFEZone(tetv, to_string(Iter + 1) + " Tet " + to_string(tColor) + " " + to_string(tetColor), ColorIndex_t(tetColor % 7), { 1, 2, 3 });
				tetColor++;
				IntValList.push_back(IntVals);
				int valSize = IntValList.size();
				auto OldIndVals = IntVals;
#endif // _DEBUG

				IntTet(Vptr, Ind, nPts, VarPtrs, VolInfo, IntVals);

#ifdef _GBADEBUG
				vector<double> NewIntVals(IntVals.size(), 0.0);
				for (int i = 0; i < NewIntVals.size(); ++i)
					NewIntVals[i] = IntVals[i] - OldIndVals[i];

				for (double i : NewIntVals)
					DebugLog << "," << i;
				DebugLog << std::endl;
#endif
			}
		}

		Planes[0] = Planes[1];

#ifdef _GBADEBUG
		Iter++;
		// 		if (Iter > 20) break;
#endif
	}

#ifdef _GBADEBUG
	DebugLog.close();
#endif

	return;
}

/*
 *	Same as above, but simply stepping down GPs one point at a time without regard to value of rho.
 *	Used for GBA after adding the ring surfaces triangulation update code. No gradient
 *	paths are seeded down the edges of triangular sphere elements, so there's always
 *	3 or 4 GPs.
 */
vector<vector<vector<int> > > NewIntegration1Tets = { // These are tet mappings once a midpoint vertex has been added
		{// for 3 GP GBs. midpoint vertex is 6
			{0,1,2,6}, // first 6 for the 3 quad sides of the trigonal prism
			{1,2,3,6},
			{2,3,4,6},
			{3,4,5,6},
			{4,5,0,6},
			{5,0,1,6},
			{0,2,4,6}, // then 2 for the triangular top and bottom
			{1,3,5,6},
			{0,1,5,7}  // then an extra tet for the 3-4 GP transition GB, extra transition vertex is 7
		},
		{// for 4 GP GBs. midpoint vertex is 8
			{0,1,2,8}, // first 8 are wrapping around the "sides" of the rectangular prism
			{1,2,3,8},
			{2,3,4,8},
			{3,4,5,8},
			{4,5,6,8},
			{5,6,7,8},
			{6,7,0,8},
			{7,0,1,8},
			{0,2,4,8}, // then 4 more for the top and bottom
			{0,4,6,8},
			{1,3,5,8},
			{1,5,7,8}
		}
};
vector<vector<int> > NewIntegration1TetDegenCheckGPs = {
	{
		{0,1,1,2,2,1,-1,-1},
		{0,1,1,2,2,3,3,0,-1,-1,-1,-1}
	}
};
void NewIntegrateUsingIsosurfaces1(vector<GradPath_c const *> const & GPsIn,
	int nPts,
	VolExtentIndexWeights_s & VolInfo,
	vector<FieldDataPointer_c> const & VarPtrs,
	vector<double> & IntVals,
	vec3 const * BondCPPos,
	AddOn_pa * AddOnID)
{
#ifdef _GBADEBUG
	std::ofstream DebugLog("C:\\Users\\Haiiro\\Documents\\Tecplot Addons\\debuglog.csv", std::ofstream::app);
	DebugLog << "\n\nNew run:\n\n";
#endif

	auto GPs = GPsIn;
	REQUIRE(GPs.size() == 3 || GPs.size() == 4);
	for (auto const * p : GPs) {
		REQUIRE(p->IsMade() && p->GetCount() > 3);
	}

	if (IntVals.size() != VarPtrs.size() + 1)
		IntVals.resize(VarPtrs.size() + 1, 0);

	/*
	*	We're stepping down the GPs according to values of rho.
	*	We'll assume that they all either terminate at the same rho value or critical point.
	*	Even if they don't terminate at the same rho value, we'll assume that the termination of any GP will be at a sufficiently low rho value
	*		such that most properties of interest have converged.
	*	Also assuming that GPs go "downhill" from a GBA sphere around a nuclear CP.
	*	At end of each loop, confirm that there is at least one more point in each grad path, otherwise exit.
	*/
	bool NotTerminated = true;

	// Used in avoiding very small volume polyhedra.
	double NewStepCutoffRatio = 0.7;

	unsigned int NumGPs = GPs.size();

	int BondPointNodeNum = -1;

	bool ContainsBondPath = (BondCPPos != nullptr);

	if (ContainsBondPath){
		// Get GP point number that corresponds to the bond point.
		vec3 tmpvec = GPs[0]->ClosestPoint(*BondCPPos, BondPointNodeNum);
	}


	double ChkDistSqr, TmpDistSqr;

	/*
		*	Because the triangulation of each plane can change as we
		*	move down the GB, we need to find a triangulation to use
		*	throughout the entire process.
		*	If we're in an open system then many of the GBs will terminate
		*	at a threshold value of rho. In this case we can use the
		*	terminal plane to determine the triangulation.
		*	If the GB terminates at a cage CP then the process of finding
		*	a good plane to use for the triangulation becomes more complicated,
		*	so we'll just punt and use the triangulation from the halfway plane
		*	according to the number of points in the first GP.
		*
		*	There's a special case for GBs that terminate at a point, but it's easy to
		*	deal with. The triangulation
		*	will result in tetrahedra that can be immediately integrated rather than
		*	trigonal prisms.
		*/

		/*
		* Check to see if the GB originates or terminates at a point (i.e. cage CP) by
		* testing for degenerate terminal points.
		*/
// 	double TermCheckDist = 0.01,
// 		AvgTermDist = 0.0;
// 	double Denom = 0.0;
// 	for (int i = 0; i < GPs.size() - 1; ++i) {
// 		Denom += GPs.size() - (i + 1);
// 		for (int j = i + 1; j < GPs.size(); ++j) {
// 			AvgTermDist += Distance(GPs[i]->XYZAt(-1), GPs[j]->XYZAt(-1));
// 		}
// 	}
// 	AvgTermDist /= Denom;
// 	bool TermAtPoint = (AvgTermDist < TermCheckDist);

#ifdef _GBADEBUG
	vector<vector<double> > IntValList;
	high_resolution_clock::time_point Time1;
	bool UserQuit;
	if (AddOnID != nullptr) {
		Time1 = high_resolution_clock::now();
		UserQuit = !StatusUpdate(0, GPs[0]->GetCount(), "Integrating GB", *AddOnID, Time1);
	}
#endif

	/*
	*	Main loop
	*/
	int PtInd = 0;
	vector<vector<int> > IndLists(2, vector<int>(NumGPs, PtInd));
	vector<vec3 const*> Vptr;
	Vptr.reserve(6);
	vector<vector<int> > Tets;
	vector<bool> IsDegenerate;
	vector<int> DegenerateParent;
	while (NotTerminated) {
		for (int gpi = 0; gpi < NumGPs; ++gpi) {
			if (IndLists[1][gpi] < GPs[gpi]->GetCount() - 1)
				IndLists[1][gpi]++;
		}

		vector<vector<vec3> > Planes(2);
		for (int gpi = 0; gpi < NumGPs; ++gpi) {
			for (int i = 0; i < 2; ++i) {
				Planes[i].push_back(GPs[gpi]->XYZAt(IndLists[i][gpi]));
			}
		}

#ifdef _GBADEBUG
		if (AddOnID != nullptr)
			UserQuit = !StatusUpdate(PtInd, GPs[0]->GetCount(), "Integrating GB", *AddOnID, Time1);
		// Save each plane as scatter
		for (int i = 0; i < 2; ++i) {
// 			SaveVec3VecAsScatterZone(Planes[i], to_string(PtInd + 1) + " Plane " + to_string(i + 1), ColorIndex_t(i % 7), { 1, 2, 3 });
		}
		int tColor = 0;
#endif

// #ifdef _DEBUG
// 		// Save triangle as scatter
// 		vector<vec3> tv;
// 		for (int i = 0; i < 3; ++i) tv.push_back(Planes[1][t[i]]);
// 		SaveVec3VecAsScatterZone(tv, to_string(PtInd + 1) + " Tri " + to_string(tColor), ColorIndex_t(tColor % 7), { 1, 2, 3 });
// 		tColor++;
// #endif

#ifdef _GBADEBUG
		// for coloring tets
		int tetColor = 0;
#endif // _DEBUG

		Vptr.clear();
		vec3 MidPt = zeros(3);
		Tets.clear();
		IsDegenerate.clear();
		DegenerateParent.clear();

		if (!ContainsBondPath || (BondPointNodeNum > 0 && IndLists[0][0] <= BondPointNodeNum)){
			// 3 GP or transition GB case
			for (int gpi = 0; gpi < 3; ++gpi) {
				for (int i = 0; i < 2; ++i) {
					Vptr.push_back(&Planes[i][gpi]);
				}
			}

			// now add the midpoint vertex, excluding degenerate copies of vertices
			IsDegenerate.resize(Vptr.size(), false);
			DegenerateParent.resize(Vptr.size(), -1);
			for (int vi = 0; vi < Vptr.size() - 1; ++vi) {
				if (!IsDegenerate[vi]) {
					for (int vj = vi + 1; vj < Vptr.size(); ++vj) {
						if (approx_equal(*Vptr[vi], *Vptr[vj], "absdiff", 1e-10)) {
							IsDegenerate[vj] = true;
							DegenerateParent[vj] = vi;
						}
					}
				}
			}
			int NumNonDegenVerts = 0;
			for (int vi = 0; vi < Vptr.size(); ++vi) {
				if (!IsDegenerate[vi]) {
					MidPt += *Vptr[vi];
					NumNonDegenVerts++;
				}
			}
			MidPt /= (double)NumNonDegenVerts;
			Vptr.push_back(&MidPt);
			DegenerateParent.push_back(-1);

			Tets.insert(Tets.begin(), NewIntegration1Tets[0].begin(), NewIntegration1Tets[0].end() - 1);

			// for 3-4 GP transition GB, check if the extra transition
			// tet is degenerate and add if not
			if (ContainsBondPath && IndLists[0][0] == BondPointNodeNum
				&& IndLists[0][0] != IndLists[1][0]
				&& IndLists[0].back() != IndLists[1].back())
			{
				Vptr.push_back(&Planes[1][3]);
				DegenerateParent.push_back(-1);
				Tets.push_back(NewIntegration1Tets[0].back());
			}
		}
		else{
			// There is a bond path and we're now in the 4 GP section of the GB
			for (int gpi = 0; gpi < 4; ++gpi) {
				for (int i = 0; i < 2; ++i) {
					Vptr.push_back(&Planes[i][gpi]);
				}
			}

			// now add the midpoint vertex, excluding degenerate copies of vertices
			IsDegenerate.resize(Vptr.size(), false);
			DegenerateParent.resize(Vptr.size(), -1);
			for (int vi = 0; vi < Vptr.size() - 1; ++vi){
				if (!IsDegenerate[vi]) {
					for (int vj = vi + 1; vj < Vptr.size(); ++vj) {
						if (approx_equal(*Vptr[vi], *Vptr[vj], "absdiff", 1e-12)){
							IsDegenerate[vj] = true;
							DegenerateParent[vj] = vi;
						}
					}
				}
			}
			int NumNonDegenVerts = 0;
			for (int vi = 0; vi < Vptr.size(); ++vi) {
				if (!IsDegenerate[vi]) {
					MidPt += *Vptr[vi];
					NumNonDegenVerts++;
				}
			}
			MidPt /= (double)NumNonDegenVerts;
			Vptr.push_back(&MidPt);
			DegenerateParent.push_back(-1);
			Tets = NewIntegration1Tets[1];
		}

#ifdef _GBADEBUG
		vector<vec3> V;
		for (auto v : Vptr) V.push_back(*v);
		SaveVec3VecAsScatterZone(V, to_string(PtInd + 1) + " Vertices " + to_string(tColor + 1), ColorIndex_t(tColor % 7), { 1, 2, 3 });
#endif	

		for (auto const & tet : Tets)
		{
			bool TetIsDegen = false;
			for (int ti = 0; ti < 3 && !TetIsDegen; ++ti){
				for (int tj = ti + 1; tj < 4 && !TetIsDegen; ++tj){
					TetIsDegen = (DegenerateParent[tet[ti]] == tet[tj] || DegenerateParent[tet[tj]] == tet[ti]);
					if (!TetIsDegen && IndLists[0][0] == 0)
						TetIsDegen = approx_equal(*Vptr[tet[ti]], *Vptr[tet[tj]], "absdiff", 1e-12);
				}
			}
#ifdef _GBADEBUG
			// Save tet as scatter
			vector<vec3> tetv;
			for (auto const & i : tet) {
				tetv.push_back(*Vptr[i]);
			}
			DebugLog << "Iter " << PtInd + 1 << ",Tet " << tColor << " " << tetColor;

			if (TetIsDegen)
				SaveTetVec3VecAsFEZone(tetv, to_string(PtInd + 1) + " Skipped Tet " + to_string(tColor) + " " + to_string(tetColor), ColorIndex_t((tetColor % 7) + 1), { 1, 2, 3 });
			else
				SaveTetVec3VecAsFEZone(tetv, to_string(PtInd + 1) + " Tet " + to_string(tColor) + " " + to_string(tetColor), ColorIndex_t((tetColor % 7)+1), { 1, 2, 3 });
			tetColor++;
			IntValList.push_back(IntVals);
			int valSize = IntValList.size();
			auto OldIndVals = IntVals;
#endif // _DEBUG

			IntTet(Vptr, tet, nPts, VarPtrs, VolInfo, IntVals);

#ifdef _GBADEBUG
			vector<double> NewIntVals(IntVals.size(), 0.0);
			for (int i = 0; i < NewIntVals.size(); ++i)
				NewIntVals[i] = IntVals[i] - OldIndVals[i];

			for (double i : NewIntVals)
				DebugLog << "," << i;
			DebugLog << std::endl;
#endif
		}

		NotTerminated = false;
		for (int gpi = 0; gpi < NumGPs && !NotTerminated; ++gpi) {
			NotTerminated = IndLists[1][gpi] < GPs[gpi]->GetCount() - 1;
		}

		IndLists[0] = IndLists[1];

#ifdef _GBADEBUG
		PtInd++;
		// 		if (Iter > 20) break;
#endif
	}

#ifdef _GBADEBUG
	DebugLog.close();
#endif

	return;
}

/*
 * Integrates a hexahedron (trigonal prism) as three tetrahedrons.
 * "Bottom" triangle is indices (0,2,4), and the "top" tri is (1,3,5).
 */
void IntTrigonalPrism(vector<vec3 const*> const & InVptr, 
	int nPts, 
	vector<FieldDataPointer_c> const & VarPtrs, 
	VolExtentIndexWeights_s & VolInfo, 
	vector<double> & IntVals,
	int PtInd = 0)
{
	REQUIRE(InVptr.size() == 6);
	REQUIRE(!IntVals.empty() && IntVals.size() == VarPtrs.size() + 1);
	for (auto const & i : VarPtrs) REQUIRE(i.IsReady());

	auto Vptr = InVptr;
	vec3 MidPt = zeros(3);


	vector<bool> IsDegenerate(Vptr.size(), false);
	vector<int> DegenerateParent(Vptr.size(), -1);

	// Check for degenerate points
	for (int vi = 0; vi < Vptr.size() - 1; ++vi) {
		if (!IsDegenerate[vi]) {
			for (int vj = vi + 1; vj < Vptr.size(); ++vj) {
				if (approx_equal(*Vptr[vi], *Vptr[vj], "absdiff", 1e-10)) {
					IsDegenerate[vj] = true;
					DegenerateParent[vj] = vi;
				}
			}
		}
	}

	int NumNonDegenVerts = 0;
	for (int vi = 0; vi < Vptr.size(); ++vi) {
		if (!IsDegenerate[vi]) {
			MidPt += *Vptr[vi];
			NumNonDegenVerts++;
		}
	}

	MidPt /= (double)NumNonDegenVerts;
	Vptr.push_back(&MidPt);
	DegenerateParent.push_back(-1);
	IsDegenerate.push_back(false);

	vector<vector<int> > Tets(NewIntegration1Tets[0].cbegin(), NewIntegration1Tets[0].cend() - 1);


#ifdef _GBADEBUG
	// for coloring tets
	int tetColor = 0;
	vector<vec3> V;
	for (auto v : Vptr) V.push_back(*v);
	SaveVec3VecAsScatterZone(V, to_string(PtInd + 1) + " Vertices", ColorIndex_t(0), { 1, 2, 3 });
#endif	

	for (auto const & tet : Tets)
	{
		bool TetIsDegen = false;
		for (int ti = 0; ti < 3 && !TetIsDegen; ++ti) {
			for (int tj = ti + 1; tj < 4 && !TetIsDegen; ++tj) {
				TetIsDegen = (DegenerateParent[tet[ti]] == tet[tj] || DegenerateParent[tet[tj]] == tet[ti]);
// 				if (!TetIsDegen)
// 					TetIsDegen = approx_equal(*Vptr[tet[ti]], *Vptr[tet[tj]], "absdiff", 1e-12);
			}
		}
#ifdef _GBADEBUG
		// Save tet as scatter
		vector<vec3> tetv;
		for (auto const & i : tet) {
			tetv.push_back(*Vptr[i]);
		}

		if (TetIsDegen)
			SaveTetVec3VecAsFEZone(tetv, to_string(PtInd + 1) + " Skipped Tet " + to_string(tetColor), ColorIndex_t((tetColor % 7) + 1), { 1, 2, 3 });
		else
			SaveTetVec3VecAsFEZone(tetv, to_string(PtInd + 1) + " Tet " + to_string(tetColor), ColorIndex_t((tetColor % 7) + 1), { 1, 2, 3 });
		tetColor++;
#endif // _DEBUG

		IntTet(Vptr, tet, nPts, VarPtrs, VolInfo, IntVals);
	}
}

void NewIntegrateUsingIsosurfaces2(vector<GradPath_c const *> const & GPsIn,
	int nPts,
	VolExtentIndexWeights_s & VolInfo,
	vector<FieldDataPointer_c> const & VarPtrs,
	vector<double> & IntVals,
	vec3 const * BondCPPos,
	AddOn_pa * AddOnID)
{
#ifdef _GBADEBUG
	std::ofstream DebugLog("C:\\Users\\Haiiro\\Documents\\Tecplot Addons\\debuglog.csv", std::ofstream::app);
	DebugLog << "\n\nNew run:\n\n";
#endif

	auto GPs = GPsIn;
	ENSURE(GPs.size() >= 3);
	for (auto const * p : GPs) {
		ENSURE(p->IsMade() && p->GetCount() > 3);
	}

	if (IntVals.size() != VarPtrs.size() + 1)
		IntVals.resize(VarPtrs.size() + 1, 0);

	/*
	*	We're stepping down the GPs according to values of rho.
	*	We'll assume that they all either terminate at the same rho value or critical point.
	*	Even if they don't terminate at the same rho value, we'll assume that the termination of any GP will be at a sufficiently low rho value
	*		such that most properties of interest have converged.
	*	Also assuming that GPs go "downhill" from a GBA sphere around a nuclear CP.
	*	At end of each loop, confirm that there is at least one more point in each grad path, otherwise exit.
	*/
	bool NotTerminated = true;

	// Used in avoiding very small volume polyhedra.
	double NewStepCutoffRatio = 0.7;

	unsigned int NumGPs = GPs.size();

// 	/*
// 	  * Check to see if the GB terminates at a point (i.e. cage CP) by
// 	  * testing for degenerate terminal points.
// 	  */
// 	double TermCheckDist = 0.01,
// 		AvgTermDist = 0.0;// (sum(GPs[0]->XYZAt(-1) == GPs[1]->XYZAt(-1)) == 3);
// 	double Denom = 0.0;
// 	for (int i = 0; i < GPs.size() - 1; ++i) {
// 		Denom += GPs.size() - (i + 1);
// 		for (int j = i + 1; j < GPs.size(); ++j) {
// 			AvgTermDist += Distance(GPs[i]->XYZAt(-1), GPs[j]->XYZAt(-1));
// 		}
// 	}
// 	AvgTermDist /= Denom;
// 	bool TermAtPoint = (AvgTermDist < TermCheckDist);

	vector<vector<int> > T;

	bool ContainsBondPath = (BondCPPos != nullptr);

// 	if (!TermAtPoint){
// 		vector<vec3> V;
// 		V.reserve(NumGPs);
// 		for (auto const & p : GPs) V.push_back(p->XYZAt(-1));
// 		T = TriangulatePolygon(V);
// 	}
// 	else {

	// Get indices of maximum GB cross-sectional area.
	int TmpCount = GPs[0]->GetCount();
	vector<double> AreaAtInd(TmpCount, 0.0);
	vector<vector<vector<int> > > Triangulations(TmpCount);
	vector<vector<unsigned int> > TmpIndLists(TmpCount, vector<unsigned int>(GPs.size(), 0));
	vector<vector<vec3> > ThVertices(omp_get_num_procs(), vector<vec3>(GPs.size()));
	int StartInd = 0, EndInd = TmpCount;
	if (ContainsBondPath){
		vec3 ClosestPoint = GPs[0]->ClosestPoint(*BondCPPos, StartInd);
		if (TmpCount - StartInd > 30) {
			StartInd += 5;
		}
	}if (TmpCount - StartInd > 30) {
		EndInd -= 5;
	}
//#pragma omp parallel for
	for (int i = StartInd; i < EndInd; ++i) {
		int ThNum = omp_get_thread_num();
		TmpIndLists[i][0] = i;
		if (GetIsoRhoGPPoints(GPsIn, TmpIndLists[i], 0, ThVertices[ThNum], 0.5, true)){
			Triangulations[i] = TriangulatePolygon(ThVertices[ThNum]);
			for (auto const & t : Triangulations[i]){
				AreaAtInd[i] += TriArea(ThVertices[ThNum][t[0]], ThVertices[ThNum][t[1]], ThVertices[ThNum][t[2]]);
			}
		}
	}

	int MaxAreaInd = -1;
	double MaxArea = 0.0;
	for (auto i = 0; i < TmpCount; ++i){
		if (AreaAtInd[i] > MaxArea){
			MaxArea = AreaAtInd[i];
			MaxAreaInd = i;
		}
	}

	if (MaxAreaInd > 0){
		T = Triangulations[MaxAreaInd];
	}
	else{
		T = Triangulations[(StartInd + EndInd) / 2];
// #pragma omp critical(BadGBCrosSectionArea)
// 		{
// 			TecUtilDialogErrMsg("Failed to find max area cross-section of GB");
// 		}
	}
// 
// 	if (ContainsBondPath) {
// 		vector<vec3> V;
// 		V.reserve(NumGPs);
// 		int ClosestPointToBond = -1;
// 		for (auto const & p : GPs) {
// 			if (approx_equal(p->ClosestPoint(*BondCPPos, ClosestPointToBond), *BondCPPos, "absdiff", 1e-12)) {
// 				int CheckPoint = (ClosestPointToBond + p->GetCount()) / 2;
// 				for (auto const & q : GPs)
// 					V.push_back(q->XYZAt(CheckPoint > q->GetCount() ? q->GetCount() / 2 : CheckPoint));
// 			}
// 			V.push_back(p->XYZAt(p->GetCount() / 2));
// 		}
// 		T = TriangulatePolygon(V);
// 	}
// 	else {
// 		vector<vec3> V;
// 		V.reserve(NumGPs);
// 		int HalfPoint = GPs[0]->GetCount() / 2;
// 		for (auto const & p : GPs)
// 			V.push_back(p->XYZAt(HalfPoint > p->GetCount() ? p->GetCount() / 2 : HalfPoint));
// 		T = TriangulatePolygon(V);
// 	}
// 	}


	double ChkDistSqr, TmpDistSqr;

	/*
		*	Because the triangulation of each plane can change as we
		*	move down the GB, we need to find a triangulation to use
		*	throughout the entire process.
		*	If we're in an open system then many of the GBs will terminate
		*	at a threshold value of rho. In this case we can use the
		*	terminal plane to determine the triangulation.
		*	If the GB terminates at a cage CP then the process of finding
		*	a good plane to use for the triangulation becomes more complicated,
		*	so we'll just punt and use the triangulation from the halfway plane
		*	according to the number of points in the first GP.
		*
		*	There's a special case for GBs that terminate at a point, but it's easy to
		*	deal with. The triangulation
		*	will result in tetrahedra that can be immediately integrated rather than
		*	trigonal prisms.
		*/

		/*
		* Check to see if the GB originates or terminates at a point (i.e. cage CP) by
		* testing for degenerate terminal points.
		*/
		// 	double TermCheckDist = 0.01,
		// 		AvgTermDist = 0.0;
		// 	double Denom = 0.0;
		// 	for (int i = 0; i < GPs.size() - 1; ++i) {
		// 		Denom += GPs.size() - (i + 1);
		// 		for (int j = i + 1; j < GPs.size(); ++j) {
		// 			AvgTermDist += Distance(GPs[i]->XYZAt(-1), GPs[j]->XYZAt(-1));
		// 		}
		// 	}
		// 	AvgTermDist /= Denom;
		// 	bool TermAtPoint = (AvgTermDist < TermCheckDist);

#ifdef _GBADEBUG
	vector<vector<double> > IntValList;
	high_resolution_clock::time_point Time1;
	bool UserQuit;
	if (AddOnID != nullptr) {
		Time1 = high_resolution_clock::now();
		UserQuit = !StatusUpdate(0, GPs[0]->GetCount(), "Integrating GB", *AddOnID, Time1);
	}
#endif

	/*
	*	Main loop
	*/
	int PtInd = 0;
	std::deque<vector<vec3> > Planes;
	Planes.push_back(vector<vec3>());
	Planes.back().reserve(NumGPs);
	for (auto const & p : GPs) Planes.back().push_back(p->XYZAt(PtInd));
	std::deque<vector<int> > IndLists;
	IndLists.push_back(vector<int>(NumGPs, PtInd));
	vector<vec3 const*> Vptr;
	Vptr.reserve(6);
	vector<vector<int> > Tets;
	vector<bool> IsDegenerate;
	vector<int> DegenerateParent;
	while (NotTerminated) {
		IndLists.push_back(IndLists.back());
		for (int gpi = 0; gpi < NumGPs; ++gpi) {
			if (IndLists.back()[gpi] < GPs[gpi]->GetCount() - 1)
				IndLists.back()[gpi]++;
		}

		Planes.push_back(vector<vec3>());
		Planes.back().reserve(NumGPs);
		for (int gpi = 0; gpi < NumGPs; ++gpi) {
			Planes.back().push_back(GPs[gpi]->XYZAt(IndLists.back()[gpi]));
		}

#ifdef _GBADEBUG
		if (AddOnID != nullptr)
			UserQuit = !StatusUpdate(PtInd, GPs[0]->GetCount(), "Integrating GB", *AddOnID, Time1);
		// Save each plane as scatter
		for (int i = 0; i < 2; ++i) {
			 SaveVec3VecAsScatterZone(Planes[i], to_string(PtInd + 1) + " Plane " + to_string(i + 1), ColorIndex_t(i % 7), { 1, 2, 3 });
		}
		int tColor = 0;
#endif

		for (auto const & tri : T) {
#ifdef _GBADEBUG
			// Save triangle as scatter
			vector<vec3> tv;
			for (int i = 0; i < 3; ++i) tv.push_back(Planes[1][tri[i]]);
			SaveVec3VecAsScatterZone(tv, to_string(PtInd + 1) + " Tri " + to_string(tColor), ColorIndex_t(tColor % 7), { 1, 2, 3 });
			tColor++;
#endif
			Vptr.clear();
			for (auto i : tri){
				Vptr.push_back(&Planes.front()[i]);
				Vptr.push_back(&Planes.back()[i]);
			}

			IntTrigonalPrism(Vptr, nPts, VarPtrs, VolInfo, IntVals);
		}

		NotTerminated = false;
		for (int gpi = 0; gpi < NumGPs && !NotTerminated; ++gpi) {
			NotTerminated = IndLists[1][gpi] < GPs[gpi]->GetCount() - 1;
		}

		IndLists.pop_front();
		Planes.pop_front();

#ifdef _GBADEBUG
		PtInd++;
		// 		if (Iter > 20) break;
#endif
	}

#ifdef _GBADEBUG
	DebugLog.close();
#endif

	return;
}


void GetTriElementConnectivityList(vector<vector<int> > const * ElemListPtr,
	vector<vector<int> > & ElemConnectivity,
	int NumSharedCorners)
{
	REQUIRE(NumSharedCorners == 1 || NumSharedCorners == 2);
	int NumElems = ElemListPtr->size();
	ElemConnectivity.resize(NumElems);
	int NumNeighbors = (NumSharedCorners == 1 ? 12 : 3);
	for (auto & i : ElemConnectivity) i.reserve(NumNeighbors); // A triangle can have as many as 12 neighbors if NumSharedCorners == 1, else 3

#pragma omp parallel for schedule(dynamic)
	for (int i = 0; i < NumElems - 1; ++i) {
		for (int j = i + 1; j < NumElems; ++j) {
			int NumMatches = 0;
			for (int ci = 0; ci < 3 && NumMatches < NumSharedCorners; ++ci) for (int cj = 0; cj < 3 && NumMatches < NumSharedCorners; ++cj) {
				if (ElemListPtr->at(i).at(ci) == ElemListPtr->at(j).at(cj)) {
					NumMatches++;
				}
			}
			if (NumMatches >= NumSharedCorners) {
#pragma omp critical
				{
					ElemConnectivity[i].push_back(j);
					ElemConnectivity[j].push_back(i);
				}
			}
		}
	}
}


bool GetPerimeterEdges(vector<vector<int> > const & TriElems, vector<vec3> const & TriNodes, vector<vector<int> > & PerimeterEdges) {

	/*
	 *	We'll use a hash table (unordered map) to count occurrence of edges.
	 *	Any edge that only occurs once is a perimeter edge.
	 *	Also store the triangle the edge came from.
	 *	Since we only care about edges that occur once, those
	 *	that occur twice don't need the other triangle number recorded.
	 */
	bool IsOk = true;

	std::map<string,vector<int> > EdgeMap;
	for (int ti = 0; ti < TriElems.size(); ++ti){
		for (int e = 0; e < 3; ++e) {
			int e2 = (e + 1) % 3;
			string eStr = GetEdgeString(TriElems[ti][e], TriElems[ti][e2]);
			if (EdgeMap.count(eStr) == 0)
				EdgeMap[eStr] = { 1, ti };
			else
				EdgeMap[eStr][0]++;
		}
	}

	/*
	 *	Also need the edge-shared element connectivity
	 */
	vector<vector<int> > ElemConnectivity;
	GetTriElementConnectivityList(&TriElems, ElemConnectivity, 2);

	/*
	 *	We can now have a list of all the perimeter edges.
	 *	Now we can find the perimeter loop by starting at an edge and
	 *	continually finding an edge that shares a node with the current edge.
	 *	It's possible that there are "dangling" triangles with only one node
	 *	shared with neighboring triangles, which results in the edges of such
	 *	dangling triangles potentially being skipped when walking around the perimeter
	 *	because now there are multiple edges that share a node with the current edge.
	 *	We'll check for any dangling triangles by identifying nodes that appear more than two
	 *	times in the UnsortedEdgeList.
	 */
	std::map<int, int> NodeOccurenceCount;
	std::map<int, std::set<string> > NodeEdges;
	vector<vector<int> > UnsortedEdges;
	std::map<string, int> EdgeNums;
	vector<int> UnsortedEdgeTriNums;
	for (auto const e : EdgeMap) {
		if (e.second[0] == 1) {
			vector<int> edgeNodes = SplitStringInt(e.first);
			EdgeNums[e.first] = UnsortedEdges.size();
			UnsortedEdges.push_back(edgeNodes);
			UnsortedEdgeTriNums.push_back(e.second[1]);
			for (auto const & i : edgeNodes) {
				NodeOccurenceCount[i]++;
				NodeEdges[i].insert(e.first);
			}
		}
	}

	std::set<int> RepeatNodes;
	for (auto const & n : NodeOccurenceCount) {
		if (n.second > 2) {
			RepeatNodes.insert(n.first);
		}
	}


// 	vector<vector<int> > RepeatNodeEdgeOrders;
// 	vector<int> RepeatNodes;
// 	for (auto const & n : NodeOccurenceCount) {
// 		if (n.second > 2) {
// 		/*
// 		 *	There are dangling elements.
// 		 *	This changes two things: first we need to guarantee that we
// 		 *	don't start the perimeter walk from an edge that shares the repeated
// 		 *	node, and second we need to determine the correct way to walk
// 		 *	through the edges sharing the repeated node.
// 		 *	To figure out the correct traversal of the dangling element,
// 		 *	we'll just move according to minimum distance edge midpoint.
// 		 *	This assumes that the correct edge to move to next is the edge
// 		 *	whose midpoint is the closest to the midpoint of the current edge.
// 		 */
// 
// 			/*
// 			 *	First collect the edges that share the repeated node
// 			 */
// 			RepeatNodes.push_back(n.first);
// 			vector<int> RepeatNodeNeighbors, RepeatNodeEdgeNums;
// 			vector<const vector<int>* > RepeatNodeEdges;
// 			for (int e = 0; e < UnsortedEdges.size(); ++e) {
// 				for (int ei = 0; ei < 2; ++ei) {
// 					if (UnsortedEdges[e][ei] == n.first) {
// 						RepeatNodeEdges.push_back(&UnsortedEdges[e]);
// 						RepeatNodeEdgeNums.push_back(e);
// 						RepeatNodeNeighbors.push_back(UnsortedEdges[e][(ei + 1) % 2]);
// 						break;
// 					}
// 				}
// 				if (RepeatNodeEdges.size() == n.second)
// 					break;
// 			}
// 
// 			/*
// 			 *	Now find the two edges in this neighborhood that share another edge;
// 			 *	that shared edge is the third edge of the dangling element.
// 			 */
// 			for (int i = 0; i < RepeatNodeNeighbors.size() - 1 && RepeatNodeEdgeOrders.size() < RepeatNodes.size(); ++i) {
// 				for (int j = i + 1; j < RepeatNodeNeighbors.size() && RepeatNodeEdgeOrders.size() < RepeatNodes.size(); ++j) {
// 					for (int e = 0; e < UnsortedEdges.size() && RepeatNodeEdgeOrders.size() < RepeatNodes.size(); ++e) {
// 						for (int ei = 0; ei < 2 && RepeatNodeEdgeOrders.size() < RepeatNodes.size(); ++ei) {
// 							if (UnsortedEdges[e][ei] == RepeatNodeNeighbors[i] 
// 								&& UnsortedEdges[e][(ei + 1) % 2] == RepeatNodeNeighbors[j]) {
// 								/*
// 								 *	Now e is the third edge of the dangling element and we can know
// 								 *	the right way to traverse all these edges.
// 								 *	The traversal starts from one of edges that is not on the dangling element,
// 								 *	the second edge is the closer of the two dangling element edges that also
// 								 *	shares the repeat node,
// 								 *	the third element is the third edge of the dangling element,
// 								 *	the fourth is the remaining edge of the dangling element,
// 								 *	and the fifth is the only remaining edge.
// 								 */
// 								for (int k = 0; k < RepeatNodeNeighbors.size(); ++k) {
// 									if (k != i && k != j) {
// 										/*
// 										 *	This is one of the two edges that is not on the dangling element.
// 										 *	It will be the first edge in the local ordering around the dangling element.
// 										 *	Now check the distance from the midpoint of this edge to the midpoints of
// 										 *	the two edges on the dangling element.
// 										 *	The shorter of the distances will determine the next edge in the local ordering,
// 										 *	followed by the third edge 'e', then the other dangling element edge, then 
// 										 *	the remaining edge that shares the repeat node.
// 										 */
// 										RepeatNodeEdgeOrders.push_back({ RepeatNodeEdgeNums[k] });
// 										vec3 kVec = 0.5 * (TriNodes[RepeatNodeEdges[k]->at(0)] + TriNodes[RepeatNodeEdges[k]->at(1)]);
// 										vector<int> ij = { i,j };
// 										vector<double> ijDist;
// 										for (int ii : ij) {
// 											ijDist.push_back(DistSqr(kVec, 0.5 * (TriNodes[RepeatNodeEdges[ii]->at(0)] + TriNodes[RepeatNodeEdges[ii]->at(1)])));
// 										}
// 										if (ijDist[0] < ijDist[1]) {
// 											RepeatNodeEdgeOrders.back().push_back(RepeatNodeEdgeNums[i]);
// 											RepeatNodeEdgeOrders.back().push_back(e);
// 											RepeatNodeEdgeOrders.back().push_back(RepeatNodeEdgeNums[j]);
// 										}
// 										else {
// 											RepeatNodeEdgeOrders.back().push_back(RepeatNodeEdgeNums[j]);
// 											RepeatNodeEdgeOrders.back().push_back(e);
// 											RepeatNodeEdgeOrders.back().push_back(RepeatNodeEdgeNums[i]);
// 										}
// 										/*
// 										 *	Now all the edges except one from RepeatNodeEdgeNums has been added,
// 										 *	so search for the remaining edge and add it.
// 										 */
// 										for (int ii : RepeatNodeEdgeNums)
// 											if (std::find(RepeatNodeEdgeOrders.back().begin(), RepeatNodeEdgeOrders.back().end(), ii) == RepeatNodeEdgeOrders.back().end())
// 												RepeatNodeEdgeOrders.back().push_back(ii);
// 										break;
// 									}
// 								}
// 
// 							}
// 						}
// 					}
// 				}
// 			}
// 		}
// 		
// 		
// 	}

	/*
	 *	Now find a loop around the perimeter.
	 *	If the size of UnsortedEdges is zero then there is only one basin on the entire 
	 *	sphere.
	 */
	if (UnsortedEdges.size() > 0) {
		/*
		 *	If there are any dangling edges, make sure that we're not starting at one of
		 *	the edges that shares a repeated node of a dangling edge.
		 *	If we are, find one that's not in the neighborhood of any repeat node and swap.
		 */
		if (RepeatNodes.count(UnsortedEdges[0][1]) > 0) {
			for (auto & e : UnsortedEdges) {
				if (RepeatNodes.count(e[1]) == 0 && e[0] != UnsortedEdges[0][0] && e[1] != UnsortedEdges[0][0]) {
					string edgeStr = GetEdgeString(e[0], e[1]);
					if (NodeEdges[UnsortedEdges[0][1]].count(edgeStr) == 0) {
						vector<int> tmpVec = UnsortedEdges[0];
						UnsortedEdges[0] = e;
						e = tmpVec;
						break;
					}
				}
			}
		}

		PerimeterEdges.resize(UnsortedEdges.size(), vector<int>(2));
		vector<bool> EdgeAdded(UnsortedEdges.size(), false);
		int EdgeNum = 0;
		PerimeterEdges[EdgeNum][0] = UnsortedEdges[0][0];
		PerimeterEdges[EdgeNum][1] = UnsortedEdges[0][1];
		EdgeAdded[0] = true;
		int Iter = 0;
		while (EdgeNum < UnsortedEdges.size() - 1 && Iter <= UnsortedEdges.size()) {
			Iter++;
			for (int i = 0; i < UnsortedEdges.size(); ++i) {
				if (RepeatNodes.count(PerimeterEdges[EdgeNum][1]) > 0) {
					/*
					 *	The terminal node of the current edge is a repeat node, then
					 *	need to find the correct edge, of the 2 or more that share the
					 *	node.
					*/
					/*
					 *	Find the edges that include this node and that don't share an edge with
					 *	the triangle of the current edge.
					 */
					vector<vector<int> > tmpEdges;
					vector<int> tmpEdgeNums;
					string curEdgeStr = GetEdgeString(PerimeterEdges[EdgeNum][0], PerimeterEdges[EdgeNum][1]);
					int curEdgeTriNum = EdgeMap[curEdgeStr][1];
					for (auto & s : NodeEdges[PerimeterEdges[EdgeNum][1]]) {
						int tmpEdgeNum = EdgeNums[s];
						if (!EdgeAdded[tmpEdgeNum]) {
							int triNum = EdgeMap[s][1];
							if (std::find(ElemConnectivity[curEdgeTriNum].begin(), ElemConnectivity[curEdgeTriNum].end(), triNum) == ElemConnectivity[curEdgeTriNum].end()) {
								tmpEdges.push_back(SplitStringInt(s, ","));
								tmpEdgeNums.push_back(tmpEdgeNum);
							}
						}
					}

					/*
					 *	Now we've all all the candidate edges.
					 *	The correct next edge is assumed to be the one whose midpoint is closest
					 *	(Euclidean distance) to the midpoint of the current edge.
					 */
					vec3 curMidPt = (TriNodes[PerimeterEdges[EdgeNum][0]] + TriNodes[PerimeterEdges[EdgeNum][1]]) * 0.5;
					double minDist = DBL_MAX;
					int minInd = -1;
					for (int ei = 0; ei < tmpEdges.size(); ++ei) {
						vec3 midPt = (TriNodes[tmpEdges[ei][0]] + TriNodes[tmpEdges[ei][1]]) * 0.5;
						double tmpDist = DistSqr(curMidPt, midPt);
						if (tmpDist < minDist) {
							minDist = tmpDist;
							minInd = ei;
						}
					}
					if (minInd >= 0) {
						int ei = (UnsortedEdges[tmpEdgeNums[minInd]][0] == PerimeterEdges[EdgeNum][1] ? 0 : 1);
						EdgeNum++;
						EdgeAdded[tmpEdgeNums[minInd]] = true;
						PerimeterEdges[EdgeNum][0] = UnsortedEdges[tmpEdgeNums[minInd]][ei];
						PerimeterEdges[EdgeNum][1] = UnsortedEdges[tmpEdgeNums[minInd]][(ei + 1) % 2];
					}
				}
				if (!EdgeAdded[i]) {
					for (int ei = 0; ei < 2; ++ei) {
						if (UnsortedEdges[i][ei] == PerimeterEdges[EdgeNum][1]) {
							EdgeNum++;
							EdgeAdded[i] = true;
							PerimeterEdges[EdgeNum][0] = UnsortedEdges[i][ei];
							PerimeterEdges[EdgeNum][1] = UnsortedEdges[i][(ei + 1) % 2];
							break;
						}
					}
				}
			}
		}
		IsOk = (Iter < UnsortedEdges.size());
		if (IsOk) {
			for (int i = 0; i < UnsortedEdges.size(); ++i) {
				if (!EdgeAdded[i]) {
					for (int ei = 0; ei < 2; ++ei) {
						if (UnsortedEdges[i][ei] == PerimeterEdges[0][0]) {
							EdgeNum++;
							EdgeAdded[i] = true;
							PerimeterEdges[EdgeNum][0] = UnsortedEdges[i][ei];
							PerimeterEdges[EdgeNum][1] = UnsortedEdges[i][(ei + 1) % 2];
							break;
						}
					}
				}
			}
		}
	}
	return IsOk;
}

bool GetSortedParameterEdgeMidpoints(vector<vector<int> > const & TriElems, 
	vector<vec3> const & NodeList, 
	vector<vec3> & SortedEdgeMidpoints,
	vector<vector<int> > & PerimeterEdges)
{
	bool IsOk = GetPerimeterEdges(TriElems, NodeList, PerimeterEdges);

	SortedEdgeMidpoints.resize(PerimeterEdges.size());
	for (int i = 0; i < PerimeterEdges.size(); ++i) {
		SortedEdgeMidpoints[i] = 0.5 * (NodeList[PerimeterEdges[i][0]] + NodeList[PerimeterEdges[i][1]]);
	}
	return IsOk;
}

Boolean_t Vec3PathResample(vector<vec3> const & OldXYZList, int NumPoints, vector<vec3> & NewXYZList) {
	Boolean_t IsOk = NumPoints > 1;

	int OldCount = OldXYZList.size();

	// 	if (IsOk && NumPoints < OldCount){
	if (IsOk) {
		if (NumPoints == 2) {
			NewXYZList = { OldXYZList.front(),OldXYZList.back() };
		}
		else {
			NewXYZList.resize(NumPoints);

			double Length = 0;
			for (int i = 1; i < OldXYZList.size(); ++i)
				Length += Distance(OldXYZList[i], OldXYZList[i - 1]);

			double DelLength = Length / static_cast<double>(NumPoints - 1);

			double ArcLength = 0.0,
				ArcLengthI = 0.0,
				ArcLengthIm1 = 0.0;

			vec3 PtI, PtIm1;

			PtI = OldXYZList[0];

			NewXYZList[0] = PtI;

			int OldI = 0;

			for (int NewI = 1; NewI < NumPoints - 1; ++NewI) {
				ArcLength += DelLength;

				while (OldI < OldCount - 1 && ArcLengthI < ArcLength) {
					++OldI;

					ArcLengthIm1 = ArcLengthI;
					PtIm1 = PtI;

					PtI = OldXYZList[OldI];

					ArcLengthI += Distance(PtI, PtIm1);
				}

				double Ratio = (ArcLength - ArcLengthIm1) / (ArcLengthI - ArcLengthIm1);
				NewXYZList[NewI] = PtIm1 + (PtI - PtIm1) * Ratio;

				if (OldI >= OldCount) {
					while (NewI < NumPoints) {
						NewI++;
						if (NewI < NumPoints) {
							NewXYZList[NewI] = PtIm1 + (PtI - PtIm1) * Ratio;
						}
					}
				}
			}

			/*
			*	Add last point
			*/

			NewXYZList[NumPoints - 1] = OldXYZList[OldCount - 1];


			IsOk = NewXYZList.size() == NumPoints;
		}
	}
	else IsOk = FALSE;

	return IsOk;
}


bool ProjectedPointToTriangleIsInterior(vec3 const & P0, vec3 & TP, vec3 const & T1, vec3 const & T2, vec3 const & T3){
	// 3d method from http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.104.4264&rep=rep1&type=pdf

	// normal of triangle and its norm
	vec3 TN = cross(T2 - T1, T3 - T1);
	double TNLen = norm(TN);

	// first corner to P0 and its norm
	vec3 VP = P0 - T1;
	double VPLen = norm(VP);

	// alpha is the angle between TN and VP
	double CosAlpha = dot(VP, TN) / (VPLen * TNLen);

	// distance from P0 to triangle plane
	double PointTriDist = VPLen * CosAlpha;

	// vector from P0 to P1 on triangle
	vec3 PP = TN * (-PointTriDist) / TNLen;

	// point on triangle
	TP = P0 + PP;

	// for determining if interior or not
// 	vec3 V12 = T2 - T1,
// 		V21 = -V12,
// 		V13 = T3 - T1,
// 		V31 = -V13,
// 		V23 = T3 - T2,
// 		V32 = -V23;
// 
// 	vec3 V1 = V21 / norm(V21) + V31 / norm(V31),
// 		V2 = V32 / norm(V32) + V12 / norm(V12),
// 		V3 = V13 / norm(V13) + V23 / norm(V23);
// 
// 	vec3 T1P = T1 - TP,
// 		T2P = T2 - TP,
// 		T3P = T3 - TP;
// 
// 	double F1 = dot(cross(V1, -T1P), TN),
// 		F2 = dot(cross(V2, -T2P), TN),
// 		F3 = dot(cross(V3, -T3P), TN);


	bool IsInterior = dot(cross(T1 - TP, T2 - TP), TN) >= 0.0;

	return IsInterior;
}

double PointDistanceToTriangleSquared(vec3 const & P, vec3 & ClosestPoint, vec3 const & T1, vec3 const & T2, vec3 const & T3, bool RecordPoint) {
	/*
	 * from http://www.iquilezles.org/www/articles/distfunctions/distfunctions.htm
float dot2( in vec3 v ) { return dot(v,v); }

float udTriangle( in vec3 v1, in vec3 v2, in vec3 v3, in vec3 p )
{
	// prepare data
	vec3 v21 = v2 - v1; vec3 p1 = p - v1;
	vec3 v32 = v3 - v2; vec3 p2 = p - v2;
	vec3 v13 = v1 - v3; vec3 p3 = p - v3;
	vec3 nor = cross( v21, v13 );

	return sqrt( // inside/outside test
				 (sign(dot(cross(v21,nor),p1)) +
				  sign(dot(cross(v32,nor),p2)) +
				  sign(dot(cross(v13,nor),p3))<2.0)
				  ?
				  // 3 edges
				  min( min(
				  dot2(v21*clamp(dot(v21,p1)/dot2(v21),0.0,1.0)-p1),
				  dot2(v32*clamp(dot(v32,p2)/dot2(v32),0.0,1.0)-p2) ),
				  dot2(v13*clamp(dot(v13,p3)/dot2(v13),0.0,1.0)-p3) )
				  :
				  // 1 face
				  dot(nor,p1)*dot(nor,p1)/dot2(nor) );
}
*/
// 	vec3 v21 = T2 - T1,
// 		v32 = T3 - T2,
// 		v13 = T1 - T3,
// 		nor = cross(v21, v13);
// 
// 	vec3 p1 = P - T1,
// 		p2 = P - T2,
// 		p3 = P - T3;
// 
// 	double MinDistSqr = DBL_MAX;
// 
// 	// inside/outside test
// 	if (SIGN(dot(cross(v21, nor), p1))
// 		+ SIGN(dot(cross(v32, nor), p2))
// 		+ SIGN(dot(cross(v13, nor), p3))
// 		< 2.0)
// 	{
// 		// 3 edges
// 
// 		vec3 u = v21, v = -v13, w = p1;
// 		double gamma = (dot(cross(u, w), nor)) / dot2(nor),
// 			beta = (dot(cross(w, v), nor)) / dot2(nor),
// 			alpha = 1.0 - gamma - beta;
// 		ClosestPoint = (T1 * CLAMP(alpha, 0.0, 1.0)) + (T2 * CLAMP(beta, 0.0, 1.0)) + (T3 * CLAMP(gamma, 0.0, 1.0));
// 		MinDistSqr = DistSqr(ClosestPoint, P);
// 
// 
// // 		vec3 e1point = (v21 * CLAMP(dot(v21, p1) / dot2(v21), 0.0, 1.0)),
// // 			e2point = (v32 * CLAMP(dot(v32, p2) / dot2(v32), 0.0, 1.0)),
// // 			e3point = (v13 * CLAMP(dot(v13, p3) / dot2(v13), 0.0, 1.0));
// // 		double e1DistSqr = dot2(e1point - p1),
// // 			e2DistSqr = dot2(e2point - p2),
// // 			e3DistSqr = dot2(e3point - p3);
// // 
// // 		if (e1DistSqr < e2DistSqr)
// // 		{
// // 			
// // 			if (e1DistSqr < e3DistSqr)
// // 			{
// // 				ClosestPoint = e1point;
// // 				MinDistSqr = e1DistSqr;
// // 			}
// // 			else if (e1DistSqr == e3DistSqr)
// // 			{
// // 				ClosestPoint = T1;
// // 				MinDistSqr = e1DistSqr;
// // 			}
// // 			else if (e1DistSqr == e2DistSqr)
// // 			{
// // 				ClosestPoint = T2;
// // 				MinDistSqr = e1DistSqr;
// // 			}
// // 			else
// // 			{
// // 				ClosestPoint = e3point;
// // 				MinDistSqr = e3DistSqr;
// // 			}
// // 		}
// // 		else if (e2DistSqr < e1DistSqr)
// // 		{
// // 			
// // 			if (e2DistSqr < e3DistSqr)
// // 			{
// // 				ClosestPoint = e2point;
// // 				MinDistSqr = e2DistSqr;
// // 			}
// // 			else if (e2DistSqr == e3DistSqr)
// // 			{
// // 				ClosestPoint = T3;
// // 				MinDistSqr = e2DistSqr;
// // 			}
// // 			else if (e2DistSqr == e1DistSqr)
// // 			{
// // 				ClosestPoint = T2;
// // 				MinDistSqr = e2DistSqr;
// // 			}
// // 			else
// // 			{
// // 				ClosestPoint = e3point;
// // 				MinDistSqr = e3DistSqr;
// // 			}
// // 		}
// 	}
// 	else
// 	{
// 		// triangle face
// 
// 		nor = normalise(nor);
// 		ClosestPoint = P - (dot(p1, nor)) * nor;
// 		MinDistSqr = DistSqr(P, ClosestPoint);
// // 		MinDistSqr = std::pow(dot(nor, p1), 2) / dot2(nor);
// // 		if (RecordPoint) ClosestPoint = P + nor * (-sqrt(MinDistSqr)) / dot2(nor);
// 	}
// 
// 	return sqrt(MinDistSqr);
// 	


/*
 * from https://www.gamedev.net/forums/topic/552906-closest-point-on-triangle/
 */
double MinDistSqr;

vec3 edge0 = T2 - T1,
edge1 = T3 - T1,
v0 = T1 - P;

double a = dot2(edge0),
b = dot(edge0, edge1),
c = dot2(edge1),
d = dot(edge0, v0),
e = dot(edge1, v0),

Det = a * c - b * b,
s = b * e - c * d,
t = b * d - a * e;

if (s + t < Det)
{
	if (s < 0.0)
	{
		if (t < 0.0)
		{
			if (d < 0.0)
			{
				s = CLAMP(-d / a, 0.0, 1.0);
				t = 0.0;
			}
			else
			{
				s = 0.0;
				t = CLAMP(-e / c, 0.0, 1.0);
			}
		}
		else
		{
			s = 0.0;
			t = CLAMP(-e / c, 0.0, 1.0);
		}
	}
	else if (t < 0.0)
	{
		s = CLAMP(-d / a, 0.0, 1.0);
		t = 0.0;
	}
	else
	{
		double invDet = 1.0 / Det;
		s *= invDet;
		t *= invDet;
	}
}
else
{
	if (s < 0.0)
	{
		double tmp0 = b + d;
		double tmp1 = c + e;
		if (tmp1 > tmp0)
		{
			double numer = tmp1 - tmp0;
			double denom = a - 2 * b + c;
			s = CLAMP(numer / denom, 0.0, 1.0);
			t = 1 - s;
		}
		else
		{
			t = CLAMP(-e / c, 0.0, 1.0);
			s = 0.0;
		}
	}
	else if (t < 0.0)
	{
		if (a + d > b + e)
		{
			double numer = c + e - b - d;
			double denom = a - 2 * b + c;
			s = CLAMP(numer / denom, 0.0, 1.0);
			t = 1 - s;
		}
		else
		{
			s = CLAMP(-e / c, 0.0, 1.0);
			t = 0.0;
		}
	}
	else
	{
		double numer = c + e - b - d;
		double denom = a - 2 * b + c;
		s = CLAMP(numer / denom, 0.0, 1.0);
		t = 1.0 - s;
	}
}

ClosestPoint = T1 + s * edge0 + t * edge1;
return Distance(ClosestPoint, P);

}

void TriangleEdgeMidPointSubdivide(
	vector<vec3> & nodes,
	vector<vector<int> > & elems,
	int TriNum,
	vector<int> & newTriNums)
{
	vector<int> NewNodes(3);
	newTriNums.resize(4);
	newTriNums[0] = TriNum;
	int ni = nodes.size();
	int ti = elems.size();
	for (int i = 0; i < 3; ++i) {
		nodes.push_back(
			(
				nodes[elems[TriNum][i]]
				+ nodes[elems[TriNum][(i + 1) % 3]]
				) / 2.
		);
		NewNodes[i] = ni++;
		newTriNums[i + 1] = ti++;
	}
	elems.push_back(vector<int>({ elems[TriNum][1], NewNodes[1], NewNodes[0] }));
	elems.push_back(vector<int>({ elems[TriNum][2], NewNodes[2], NewNodes[1] }));
	elems.push_back(NewNodes);
	elems[TriNum] = vector<int>({ elems[TriNum][0], NewNodes[0], NewNodes[2] });
	return;
}


/*
 *	Perform midpoint subdivision for all triangular
 *	elements that contain nodeNum.
 */
void TriangleMidPointSubdivide(
	vector<vec3> & nodes,
	vector<vector<int> > & tris,
	vector<std::set<int> > & trisOfNode,
	int triNum,
	int & newNodeNum,
	vector<int> & newTriNums,
	vec3 * newNode = nullptr)
{
	newTriNums.clear();

	auto t = tris[triNum];

	for (auto ni : t){
		trisOfNode[ni].erase(triNum);
	}

	newNodeNum = nodes.size();
	if (newNode != nullptr){
		nodes.push_back(*newNode);
	}
	else {
		nodes.push_back((nodes[t[0]] + nodes[t[1]] + nodes[t[2]]) / 3.0);
	}

	newTriNums.push_back(tris.size());
	tris.push_back({ t[0], t[1], newNodeNum });
	newTriNums.push_back(tris.size());
	tris.push_back({ t[1], t[2], newNodeNum });
	newTriNums.push_back(triNum);
	tris[triNum] = { t[2], t[0], newNodeNum };

	// Update triOfNode
	trisOfNode.resize(nodes.size());
	for (auto const & ti : newTriNums) {
		for (auto ni : tris[ti]) {
			trisOfNode[ni].insert(ti);
		}
	}

}


/*
 *	Perform edge midpoint subdivision for all triangular
 *	elements that contain nodeNum.
 */
void TriangleEdgeMidPointSubdivideAroundNodes(
	vector<vec3> & nodes,
	vector<vector<int> > & tris,
	vector<std::set<int> > & trisOfNode,
	vector<int> const & nodeNums,
	vector<int> & newNodeNums,
	vector<int> & newTriNums)
{
	newNodeNums.clear();
	newTriNums.clear();
	// Get list of unique triangles, each with their edges,
	// and unique edges, with the node index of their midpoint.
	std::map<int, vector<std::pair<int, int> > > TriEdgeMap;
	std::map<std::pair<int, int>, int> EdgeMidpointNodeNumMap;
	for (auto ni : nodeNums) {
		for (auto ti : trisOfNode[ni]) {
			for (int ci = 0; ci < 3; ++ci) {
				int c1 = tris[ti][ci],
					c2 = tris[ti][(ci + 1) % 3];
				std::pair<int, int> e(c1 < c2 ? std::make_pair(c1, c2) : std::make_pair(c2, c1));
				TriEdgeMap[ti].push_back(e);
				if (EdgeMidpointNodeNumMap.count(e) == 0) {
					nodes.push_back((nodes[c1] + nodes[c2]) * 0.5);
					newNodeNums.push_back(nodes.size() - 1);
					EdgeMidpointNodeNumMap[e] = newNodeNums.back();
				}
			}
		}
	}

	// Update tris
	for (auto const & t : TriEdgeMap){
		// Remove old triangle index from trisOfNode
		for (auto ni : tris[t.first]) trisOfNode[ni].erase(t.first);

		// Three new corner elements made from the existing corner
		// and two connected edge midpoints.
		// Existing triangle is replaced with the new midpoint nodes.
		for (int ci = 0; ci < 3; ++ci){
			tris.push_back({
				tris[t.first][ci],
				EdgeMidpointNodeNumMap[t.second[ci]],
				EdgeMidpointNodeNumMap[t.second[(ci+2) % 3]]
				}
			);
			newTriNums.push_back(tris.size() - 1);
		}
		for (int ci = 0; ci < 3; ++ci){
			tris[t.first][ci] = EdgeMidpointNodeNumMap[t.second[ci]];
		}
		newTriNums.push_back(t.first);
	}

	// Update triOfNode
	trisOfNode.resize(nodes.size());
	for (auto const & ti : newTriNums){
		for (auto ni : tris[ti]){
			trisOfNode[ni].insert(ti);
		}
	}

}

/*
 * Updates a spherical triangulation using a set of constraint nodes (and segments)
 * such that the resulting spherical triangulation has nodes on each of the constraint nodes
 * and, if constraint segments are provided, no edges that cross constraint edges.
 */
void UpdateSubdivideSphericalTriangulationWithConstraintNodesAndSegments(
											vector<vec3> const & inNodes, // input sphere nodes
											vector<vector<int> > const & inTris, // input sphere triangular elements in terms of node indices
											vec3 const & sphereCenter,
											double const & sphereRadius,
											vector<vector<vec3> > & constraintNodes, // input 2d vector<vector<int>>; constraint nodes organized by segments
											vector<vec3> & outNodes, // new sphere nodes
											vector<vector<int> > & outTris, // new sphere elements in terms of outNodes indices
											std::map<int, int> & outNodeConstraintNodeIndices, // indicates which constraintNode each outNode corresponds to (-1 if none or corresponds to segment)
											std::map<int, int> & outNodeConstraintSegmentIndices) // indicates which constraint segment each outNode corresponds to (-1 if no correspondance to any constraint)
{
	outNodes = inNodes;
	outTris = inTris;
 
 	/*
 	 *	We'll keep a reverse element list; a list of triangles
 	 *	for each node.
 	 */
 	vector<std::set<int> > trisOfNode(inNodes.size());
 	for (int ti = 0; ti < inTris.size(); ++ti){
 		for (auto ni : inTris[ti]){
 			trisOfNode[ni].insert(ti);
 		}
 	}

	/*
	 *	The segments of constraint nodes can share endpoints.
	 *	We'll check for duplicates and only use unique constraint
	 *	node positions.
	 */
	std::set<std::pair<int,int> > dupConstraintNodes;
	bool CheckForCoincidentConstraintNodes = false;
	for (auto const & i : constraintNodes){
		if (i.size() > 1){
			CheckForCoincidentConstraintNodes = true;
			break;
		}
	}
	/*
	 	*	Check for nearly equal intersection path segment endpoints,
	 	*	i.e. where two paths intersect.
	 	*	These intersections should only occur at a bond-path-sphere intersection point,
	 	*	at which there is already a constrained sphere node, so when path intersections are
	 	*	found, move all the coincident points to the closest node on the sphere.
	 	*/
	if (CheckForCoincidentConstraintNodes) {
		/*
	 *	Get average sphere triangulation edge length so that we can
	 *	resample intersection paths at roughly the same spacing.
	 *	I'll do this by looping over triangles, which will double count every
	 *	edge so the average should still be valid.
	 */
		double SphereEdgeLenMean = 0.0;
		int SphereNumEdges = 0;
		for (auto const & t : outTris) {
			for (int i = 0; i < 3; ++i) {
				SphereEdgeLenMean += Distance(outNodes[t[i]], outNodes[t[(i + 1) % 3]]);
				SphereNumEdges++;
			}
		}
		SphereEdgeLenMean /= (double)SphereNumEdges;
		double CoincidentCheckEpsilon = SphereEdgeLenMean * 0.1;
		CoincidentCheckEpsilon *= CoincidentCheckEpsilon;
		for (int i = 0; i < constraintNodes.size() - 1; ++i) {
			vector<std::pair<int, int> > CoincidentPointIndices;
			for (int j = i + 1; j < constraintNodes.size(); ++j) {
				if (DistSqr(constraintNodes[i][0], constraintNodes[j][0]) <= CoincidentCheckEpsilon) {
					CoincidentPointIndices.push_back(std::make_pair(i, 0));
					CoincidentPointIndices.push_back(std::make_pair(j, 0));
				}
				else if (DistSqr(constraintNodes[i][0], constraintNodes[j].back()) <= CoincidentCheckEpsilon) {
					CoincidentPointIndices.push_back(std::make_pair(i, 0));
					CoincidentPointIndices.push_back(std::make_pair(j, constraintNodes[j].size() - 1));
				}
				else if (DistSqr(constraintNodes[i].back(), constraintNodes[j][0]) <= CoincidentCheckEpsilon) {
					CoincidentPointIndices.push_back(std::make_pair(i, constraintNodes[i].size() - 1));
					CoincidentPointIndices.push_back(std::make_pair(j, 0));
				}
				else if (DistSqr(constraintNodes[i].back(), constraintNodes[j].back()) <= CoincidentCheckEpsilon) {
					CoincidentPointIndices.push_back(std::make_pair(i, constraintNodes[i].size() - 1));
					CoincidentPointIndices.push_back(std::make_pair(j, constraintNodes[j].size() - 1));
				}
			}

			if (CoincidentPointIndices.size() > 0) {
				int SphereNodeNum = -1;
				for (auto const & p : CoincidentPointIndices) {
					dupConstraintNodes.insert(p);
					double MinNodeDistSqr = DBL_MAX;
					double TmpNodeDistSqr;
					int MinNodeIndex = -1;
					for (int ni = 0; ni < outNodes.size(); ++ni) {
						TmpNodeDistSqr = DistSqr(constraintNodes[p.first][p.second], outNodes[ni]);
						if (TmpNodeDistSqr < MinNodeDistSqr) {
							MinNodeDistSqr = TmpNodeDistSqr;
							MinNodeIndex = ni;
						}
					}
					if (MinNodeIndex >= 0) {
						// 					if (SphereNodeNum >= 0){
						// 						if (SphereNodeNum != MinNodeIndex) {
						// 							TecUtilDialogErrMsg("Coincident intersection path endpoints do not share same closest sphere node!");
						// 						}
						// 					}
						// 					else{
						// 						SphereNodeNum = MinNodeIndex;
						// 					}
						constraintNodes[p.first][p.second] = outNodes[MinNodeIndex];
					}
				}
			}
		}
	}
 
 	/*
 	 *	First match constraintNodes to closest inNodes.
 	 *	We'll loop over inNodes for each constraintNode to
 	 *	find the minimum Euclidean distance.
 	 *	If a single inNode is the closest node for two or more
 	 *	constraintNodes, then subdivide all the triangles
 	 *	it's a part of.
 	 *	We'll store which constraintNode each sphere node
 	 *	was closest to so that we can recheck all of the 
 	 *	constraintNodes the sphere node was closest to
 	 *	after subdividing.
 	 */
 	std::queue<std::pair<int,int> > checkConstraintSegNodeNums;
 	for (int segNum = 0; segNum < constraintNodes.size(); ++segNum) {
 		for (int nodeNum = 0; nodeNum < constraintNodes[segNum].size(); ++nodeNum) {
			std::pair<int, int> p(segNum, nodeNum);
 			checkConstraintSegNodeNums.push(p);
 		}
 	}

 	while (!checkConstraintSegNodeNums.empty()){
 		auto c = constraintNodes[checkConstraintSegNodeNums.front().first][checkConstraintSegNodeNums.front().second];
 
 		//  Get closest sphere node
 		double MinDistSqr = DBL_MAX, TmpDistSqr;
 		int minI = -1;
 		for (int ni = 0; ni < outNodes.size(); ++ni) {
 			TmpDistSqr = DistSqr(c, outNodes[ni]);
 			if (TmpDistSqr < MinDistSqr) {
 				minI = ni;
 				MinDistSqr = TmpDistSqr;
 			}
 		}
 		if (minI >= 0) {
 			if (outNodeConstraintNodeIndices.count(minI) > 0 && dupConstraintNodes.count(checkConstraintSegNodeNums.front()) == 0) {
 				// Sphere node is closest to two or more, so subdivide elements of node
// 				vector<int> newNodeNums;
// 				TriangleEdgeMidPointSubdivideAroundNodes(outNodes, outTris, trisOfNode, { minI }, newNodeNums, vector<int>());
// 
// 				// Project new points back to sphere radius
//  				for (auto ni : newNodeNums){
// 					outNodes[ni] = sphereCenter + normalise(outNodes[ni] - sphereCenter) * sphereRadius;
//  				}

				int triNum, newNodeNum;
				MinDistSqr = DBL_MAX;
				for (auto ti : trisOfNode[minI]) {
					TmpDistSqr = DistSqr(c, (outNodes[outTris[ti][0]] + outNodes[outTris[ti][1]] + outNodes[outTris[ti][2]]) / 3.0);
					if (TmpDistSqr < MinDistSqr){
						MinDistSqr = TmpDistSqr;
						triNum = ti;
					}
				}

				TriangleMidPointSubdivide(outNodes, outTris, trisOfNode, triNum, newNodeNum, vector<int>(), &c);
				outNodes[newNodeNum] = sphereCenter + normalise(outNodes[newNodeNum] - sphereCenter) * sphereRadius;
 
 				// Add both constraint nodes close to sphere node back onto queue
 				checkConstraintSegNodeNums.push(std::make_pair(outNodeConstraintSegmentIndices[minI],outNodeConstraintNodeIndices[minI]));
//  				checkConstraintSegNodeNums.push(checkConstraintSegNodeNums.front());
 
				// Reset constraint node index for sphere node
				outNodeConstraintNodeIndices.erase(minI);
				outNodeConstraintSegmentIndices.erase(minI);
 			}
 			else {
 				outNodeConstraintSegmentIndices[minI] = checkConstraintSegNodeNums.front().first;
 				outNodeConstraintNodeIndices[minI] = checkConstraintSegNodeNums.front().second;
 			}
 		}
 
 		checkConstraintSegNodeNums.pop();
 	}

	// Now for constraint node move the closest sphere node to it
//  	for (int ni = 0; ni < outNodes.size(); ++ni){
// 		if (outNodeConstraintNodeIndices.count(ni) > 0) {
// 			outNodes[ni] = constraintNodes[outNodeConstraintSegmentIndices[ni]][outNodeConstraintNodeIndices[ni]];
// 		}
//  	}

	std::map<int, int>::const_iterator ni, si = outNodeConstraintSegmentIndices.cbegin();
	for (ni = outNodeConstraintNodeIndices.cbegin(); ni != outNodeConstraintNodeIndices.cend(); ++ni){
		outNodes[ni->first] = constraintNodes[si->second][ni->second];
		si++;
	}
}

double TriangulatedSphereJiggleMesh(vector<vec3> & Nodes,
	vector<std::set<int> > const & NodeConnectivity,
	vector<bool> const & NodeIsConstrained,
	vec3 const & SphereCenter,
	double const & SphereRadius)
{
	int NumNodes = Nodes.size();
	vector<vec3> NewNodes(NumNodes);
	vector<double> NodeMovedDistSqr(NumNodes);

#pragma omp parallel for
	for (int ni = 0; ni < NumNodes; ++ni) {
		if (!NodeIsConstrained[ni]) {
			vec3 MidPt = zeros(3);
			for (auto nj : NodeConnectivity[ni])
				MidPt += Nodes[nj];

			MidPt /= (double)NodeConnectivity[ni].size();
			MidPt = SphereCenter + normalise(MidPt - SphereCenter) * SphereRadius;
			NewNodes[ni] = MidPt;
		}
		else
			NewNodes[ni] = Nodes[ni];

		NodeMovedDistSqr[ni] = DistSqr(Nodes[ni], NewNodes[ni]);
	}

	double MaxMovedDist = 0.;
	for (auto const & d : NodeMovedDistSqr)
		MaxMovedDist = MAX(MaxMovedDist, d);

	Nodes = NewNodes;

	return MaxMovedDist;
}

void GetMeshNodeConnectivity(vector<vector<int> > const & Elems,
	int NumNodes,
	vector<std::set<int> > & NodeConnectivity)
{
	NodeConnectivity.clear();
	NodeConnectivity.resize(NumNodes);

	for (auto const & e : Elems){
		for (int ci = 0; ci < e.size(); ++ci){
			for (int cj = 1; cj <= 2; ++cj){
				NodeConnectivity[e[ci]].insert(e[(ci + cj) % e.size()]);
			}
		}
	}
}

double TriArea(vec3 const & p1, vec3 const & p2, vec3 const & p3){
	vec3 ABC;

	ABC[0] = p1[2] * (p2[1] - p3[1])
		+ p2[2] * (p3[1] - p1[1])
		+ p3[2] * (p1[1] - p2[1]);

	ABC[1] = p1[0] * (p2[2] - p3[2])
		+ p2[0] * (p3[2] - p1[2])
		+ p3[0] * (p1[2] - p2[2]);

	ABC[2] = p1[1] * (p2[0] - p3[0])
		+ p2[1] * (p3[0] - p1[0])
		+ p3[1] * (p1[0] - p2[0]);

	return 0.5 * norm(ABC);
}

double TriPerimeter(vec3 const & p1, vec3 const & p2, vec3 const & p3){
	return Distance(p1, p2) + Distance(p2, p3) + Distance(p3, p1);
}

// Returns sum of differences between corner angles and 60 degrees,
// so gives the "distance" of the triangle from being equilateral.
// 0 means perfectly equilateral, and the upper limit is 240 degrees
double const MaxTriBadness = 240. / DEGPERRADIANS;
double const TriBadnessCheckAngle = 60. / DEGPERRADIANS;
double TriBadness(vec3 const & p1, vec3 const & p2, vec3 const & p3){
// 	Old implementation used the ratio of tri area to perimeter, but results
// 	in low area triangles looking more bad than they are.
// 	double a = TriArea(p1, p2, p3);
// 	double p = TriPerimeter(p1, p2, p3);
// 
// 	return (p > 0.0 ? a / p : 0.0);
// 	
// 	
	vec3 v1 = p2 - p1,
		v2 = p3 - p2,
		v3 = p1 - p3;
	double AngleDiffSum = abs(TriBadnessCheckAngle - VectorAngle(v1,-v2));
	AngleDiffSum += abs(TriBadnessCheckAngle - VectorAngle(v2,-v3));
	AngleDiffSum += abs(TriBadnessCheckAngle - VectorAngle(v3,-v1));

	return AngleDiffSum / MaxTriBadness;
}

/*
 *	Remove duplicate nodes (nodes with same position)
 *	and revise element list to reflect the removal of
 *	the nodes.
 */
void RemoveDupicateNodesFromMesh(vector<vec3> & NodeList, vector<vector<int> > & ElemList, vector<int> * OldToNewNodes) {
	vector<bool> NodeIsDuplicate(NodeList.size(), false);
	vector<int> NodeNums(NodeList.size());
	bool DupFound = false;

	for (int i = 0; i < NodeList.size(); ++i)
		NodeNums[i] = i;

	for (int i = 0; i < NodeList.size() - 1; ++i) {
		if (!NodeIsDuplicate[i]) {
			for (int j = i + 1; j < NodeList.size(); ++j) {
				if (!NodeIsDuplicate[j] && approx_equal(NodeList[i], NodeList[j], "absdiff", 1e-12)) {
					NodeIsDuplicate[j] = true;
					NodeNums[j] = NodeNums[i];
					DupFound = true;
				}
			}
		}
	}

	if (DupFound) {
		vector<vec3> NewXYZ;
		NewXYZ.reserve(NodeList.size());

		for (int n = 0; n < NodeIsDuplicate.size(); ++n) {
			if (!NodeIsDuplicate[n]) {
				NewXYZ.push_back(NodeList[n]);
				NodeNums[n] = NewXYZ.size() - 1;
				for (int ni = n + 1; ni < NodeIsDuplicate.size(); ++ni) {
					if (NodeIsDuplicate[ni] && NodeNums[ni] == n) {
						NodeNums[ni] = NodeNums[n];
					}
				}
			}
		}

		NodeList = NewXYZ;

		vector<vector<int> > NewElems;
		NewElems.reserve(ElemList.size());
		for (auto & e : ElemList) {
			for (int & ei : e) {
				ei = NodeNums[ei];
			}
			if (e[0] != e[1] && e[0] != e[2] && e[1] != e[2]) {
				NewElems.push_back(e);
			}
			else {
				int a = 1;
			}
		}

		ElemList = NewElems;

		if (OldToNewNodes != nullptr) {
			*OldToNewNodes = NodeNums;
		}
	}

	NodeList.shrink_to_fit();
	ElemList.shrink_to_fit();
}

void RemoveDuplicatePointsFromVec3Vec(vector<vec3> & Points, double tol, vector<int> * PtNumsOldToNew)
{
	if (Points.size() <= 1)
		return;

	vector<vec3> NewPoints;
	NewPoints.reserve(Points.size());

	vector<int> PtNums;
	vector<int> * PtNumPtr;
	if (PtNumsOldToNew != nullptr){
		PtNumPtr = PtNumsOldToNew;
	}
	else{
		PtNumPtr = &PtNums;
	}
	PtNumPtr->resize(Points.size());
	for (int i = 0; i < Points.size(); ++i)
		PtNumPtr->at(i) = i;

	for (int i = 0; i < Points.size() - 1; ++i){
		if (PtNumPtr->at(i) == i) {
			for (int j = i + 1; j < Points.size(); ++j) {
				if (approx_equal(Points[i], Points[PtNumPtr->at(j)], "absdiff", tol)) {
					PtNumPtr->at(j) = PtNumPtr->at(i);
				}
			}
		}
	}

	for (int i = 0; i < Points.size() - 1; ++i){
		if (PtNumPtr->at(i) == i) {
			NewPoints.push_back(Points[i]);
		}
	}

	NewPoints.push_back(Points.back());

	Points = NewPoints;
	Points.shrink_to_fit();

	return;
}

/*
 *	Distance between a point x0 and a line defined by x1 and x2.
 *	From http://mathworld.wolfram.com/Point-LineDistance3-Dimensional.html
 */
double PointLineDist(vec3 const & p, vec3 const & a, vec3 const & b) {
	double LineDist = Distance(b, a);
	if (LineDist <= 1e-12) {
		return Distance(p, a);
	}
	else {
		double numer = norm(cross(p - a, p - b));
		return numer / LineDist;
	}
}
double PointLineDist(tpcsm::Vec3 const & p, tpcsm::Vec3 const & a, tpcsm::Vec3 const & b) {
	double LineDist = (b - a).getNorm();
	if ( LineDist <= 1e-12) {
		return (a - p).getNorm();
	}
	else {
		double numer = ((p - a).cross(p - b)).getNorm();
		double denom = (b - a).getNorm();
		return numer / denom;
	}
}

vec3 ProjectPointToLine(vec3 const & p, vec3 const & a, vec3 const & b){
	vec3 ap = p - a,
		ab = b - a;
	return a + ab * dot(ap, ab) / dot(ab, ab);
}
tpcsm::Vec3 ProjectPointToLine(tpcsm::Vec3 const & p, tpcsm::Vec3 const & a, tpcsm::Vec3 const & b){
	tpcsm::Vec3 ap = p - a,
		ab = b - a;
	return a + ab * ap.dot(ab) / ab.dot(ab);
}

vec3 ClosestPointToPath(vector<vec3> const & Path, vec3 const & CheckPt, int & PtNum)
{
	vec3 ClosestPt;
	if (!Path.empty()) {
		// Simple sequential search down whole GP
		ClosestPt = Path[0];
		double MinSqrDist = DistSqr(Path[0], CheckPt);
		PtNum = 0;
		int Count = Path.size();
		for (int i = 1; i < Count - 1; ++i) {
			double TempSqrDist = DistSqr(Path[i], CheckPt);
			if (TempSqrDist < MinSqrDist) {
				MinSqrDist = TempSqrDist;
				ClosestPt = Path[i];
				PtNum = i;
			}
		}

		vec3 AB;
		double t;
		vector<bool> PointChecked(Count, false);
		do
		{
			PointChecked[PtNum] = true;
			AB = Path[PtNum + 1] - Path[PtNum];
			t = dot(CheckPt - Path[PtNum], AB) / dot(AB, AB);
			if (t > 1.0)
				if (PtNum < Count - 2 && !PointChecked[PtNum + 1])
					PtNum++;
				else
					t = 1.0;
			else if (t < 0.0)
				if (PtNum > 1 && !PointChecked[PtNum - 1])
					PtNum--;
				else
					t = 0.0;
		} while (t > 1.0 || t < 0.0);

		ClosestPt = Path[PtNum] + AB * t;
	}

	return ClosestPt;
}