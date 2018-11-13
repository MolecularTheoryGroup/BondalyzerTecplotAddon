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
#include <stdio.h>

#include "GSM_GEOMETRY.h"
#include "CSM_DATA_TYPES.h"
#include "CSM_GRAD_PATH.h"

using namespace arma;

using std::vector;

static const int ArrowRotationSteps = 10;

const bool CSMArrow(const vec3 & Origin, 
	const vec3 & Dir, 
	const double & Length, 
	const double & Radius, 
	const double & ArrowheadLengthRatio,
	const double & ArrowheadRadiusRatio,
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
const mat getPoints(const mat & Edge, const int & NumPts){
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

mat cubTrans(const mat & e1, const mat & e2, const mat & e3){

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
const mat cubTrans_func(const mat & abc,
	const vec & s,
	const vec & t,
	const vec & u)
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
const double cubJacobian(mat & abc, const double & s, const double & t, const double & u){
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

void rquad(const int & N, const double & k, vec & x, vec & w)
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
const vector<vec> tetraquad(const int & N, const mat & Verts)
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



const vector<vec> GetWeightsPoints(const int & N){

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
vector<vector<int> > TriangulatePolygon(const vector<vec3> & V){
	vector<vector<int> > T;
	double ChkDistSqr = DBL_MIN, TmpDistSqr;
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

#ifdef _DEBUG
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

const bool GetIsoRhoGPPoints(const vector<GradPath_c*> GPs, vector<unsigned int> & IndexList, const int & MinGPNum, vector<vec3> & NewVerts, const double & NewStepCutoffRatio, const bool & MultiStep = false){
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
 *	Integrate volume and variables (for provided var pointers) over a tetrahedron specified
 *	according to vertex indices and a list of vertices.
 */
void IntTet(const vector<vec3*> & Vptr, const vector<unsigned int> & Ind, const int & nPts, const vector<FieldDataPointer_c> & VarPtrs, VolExtentIndexWeights_s & VolInfo, vector<double> & IntVals){
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

#ifdef _DEBUG
	vector<vec3> GaussQuadPts(xyzw[0].n_elem);
	for (int i = 0; i < GaussQuadPts.size(); ++i)
		GaussQuadPts[i] << xyzw[0][i] << xyzw[1][i] << xyzw[2][i];
// 	SaveVec3VecAsScatterZone(GaussQuadPts, "Gauss Quad Pts", Black_C, { 1, 2, 3 });
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
void IntegrateUsingIsosurfaces(vector<GradPath_c*> & GPs,
	const int & nPts,
	VolExtentIndexWeights_s & VolInfo,
	const vector<FieldDataPointer_c> & VarPtrs,
	vector<double> & IntVals,
	const bool ContainsBondPath)
{
	REQUIRE(GPs.size() > 3);
	for (const auto * p : GPs){
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
	vector<GradPath_c*> GPi = GPs;
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
 		GradPath_c * DegenGP = NULL;
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
 		for (const auto & t : T2){
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

#ifdef _DEBUG
	int Iter = 0;
	vector<vector<double> > IntValList;
#endif

	/*
	*	Main loop
	*/
	while (NotTerminated){
		vector<vec3*> Vptr;
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


#ifdef _DEBUG
		// Save each plane as scatter
		for (int i = 0; i < 2; ++i){
// 			SaveVec3VecAsScatterZone(Planes[i], to_string(Iter + 1) + "Plane" + to_string(i + 1), ColorIndex_t(i % 7), { 1, 2, 3 });
		}
		int tColor = 0;
#endif

		for (const auto & t : T){
#ifdef _DEBUG
			// Save triangle as scatter
			vector<vec3> tv;
// 			for (int i = 0; i < 3; ++i) tv.push_back(Planes[1][t[i]]);
			// 			SaveVec3VecAsScatterZone(tv, to_string(Iter + 1) + "Tri " + to_string(tColor), ColorIndex_t(tColor % 7), { 1, 2, 3 });
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

#ifdef _DEBUG
			vector<vec3> hexv;
// 			for (const auto & i : VI){
// 				int PlaneNum = (i < NumGPs ? 0 : 1);
// 				hexv.push_back(Planes[PlaneNum][i % NumGPs]);
// 			}
// 			SaveVec3VecAsScatterZone(hexv, to_string(Iter + 1) + "Hex " + to_string(tColor + 1), ColorIndex_t(tColor % 7), { 1, 2, 3 });
#endif

			/*
			*	Which of the two sets of tetrahedra will result depends entirely on the indicices
			*	of the quadrilateral face opposite VI[0], so it's an easy code.
			*/
			int TetSetNum = int(MIN(VI[1], VI[5]) < MIN(VI[2], VI[4]));

			/*
			*	Now loop over the three tets to compute their integral
			*/
#ifdef _DEBUG
			// for coloring tets
			int tetColor = 0;
#endif // _DEBUG

			for (const auto & tet : IntTetInds[TetSetNum])
			{

				/*
				*	If checking for degenerate tet points, do it here and quit if point(s) are degenerate.
				*/
				bool TetDegen = false;
				if (ContainsBondPath || (TermAtPoint && !NotTerminated)) for (int i = 0; i < 3 && !TetDegen; ++i){
					for (int j = i + 1; j < 4 && !TetDegen; ++j){
						if (sum(*Vptr[VI[tet[i]]] == *Vptr[VI[tet[j]]]) == 3)
							TetDegen = true;
					}
				}

				if (TetDegen) 
					continue;
				
				vector<unsigned int> Ind(4);
				for (int i = 0; i < 4; ++i) Ind[i] = VI[tet[i]];

#ifdef _DEBUG
				// Save tet as scatter
				vector<vec3> tetv;
				for (const auto & i : tet){
					// 					int PlaneNum = (VI[i] < NumGPs ? 0 : 1);
					// 					tetv.push_back(Planes[PlaneNum][VI[i] % NumGPs]);
					tetv.push_back(*Vptr[VI[i]]);
				}
// 				SaveTetVec3VecAsFEZone(tetv, to_string(Iter + 1) + " Tet " + to_string(tColor) + " " + to_string(tetColor), ColorIndex_t(tetColor % 7), { 1, 2, 3 });
				tetColor++;
				IntValList.push_back(IntVals);
				int valSize = IntValList.size();
#endif // _DEBUG

				IntTet(Vptr, Ind, nPts, VarPtrs, VolInfo, IntVals);
			}
		}

		Planes[0] = Planes[1];

#ifdef _DEBUG
		Iter++;
// 		if (Iter > 20) break;
#endif
	}

	return;
}


void GetTriElementConnectivityList(const vector<vector<int> > * ElemListPtr,
	vector<vector<int> > & ElemConnectivity,
	const int & NumSharedCorners)
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

const string GetEdgeString(const int & ei, const int & ej) {
	return (ei < ej ? to_string(ei) + "," + to_string(ej) : to_string(ej) + "," + to_string(ei));
}

const bool GetPerimeterEdges(const vector<vector<int> > & TriElems, const vector<vec3> & TriNodes, vector<vector<int> > & PerimeterEdges) {

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
	for (const auto e : EdgeMap) {
		if (e.second[0] == 1) {
			vector<int> edgeNodes = SplitStringInt(e.first);
			EdgeNums[e.first] = UnsortedEdges.size();
			UnsortedEdges.push_back(edgeNodes);
			UnsortedEdgeTriNums.push_back(e.second[1]);
			for (const auto & i : edgeNodes) {
				if (NodeOccurenceCount.count(i) == 0) {
					NodeOccurenceCount[i] = 1;
					NodeEdges[i] = { e.first };
				}
				else {
					NodeOccurenceCount[i]++;
					NodeEdges[i].insert(e.first);
				}
			}
		}
	}

	std::set<int> RepeatNodes;
	for (const auto & n : NodeOccurenceCount) {
		if (n.second > 2) {
			RepeatNodes.insert(n.first);
		}
	}


// 	vector<vector<int> > RepeatNodeEdgeOrders;
// 	vector<int> RepeatNodes;
// 	for (const auto & n : NodeOccurenceCount) {
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
// 										for (const int & ii : ij) {
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
// 										for (const int & ii : RepeatNodeEdgeNums)
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

const bool GetSortedParameterEdgeMidpoints(const vector<vector<int> > & TriElems, 
	const vector<vec3> & NodeList, 
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