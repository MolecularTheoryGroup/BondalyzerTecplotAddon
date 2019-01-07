#include <vector>
#include <cmath>
#include <iostream>
#include <sstream>

#include <gsl/gsl_math.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_vector.h>

#include <omp.h>

#ifdef WIN32
#include <windows.h>
#endif

#include "meshgen2d_sphere.h"
#include "constants.h"

#include "parsearguments.h"
#include "screenoutput.h"

#define DEBUG 0

using std::vector;
using std::cout;
using std::cerr;
using std::endl;


struct OptParams{
	vector<vector<int> > AdjList;
	vector<vector<int> > TriAdjList;
	vector<char> AdjListCounts;
	vector<char> TriAdjListCounts;
	int NumPts, NumTri, OptPoint;
	int numCPU;

	point* p;
	triangle* t;
};

point const GSLVecToPoint(gsl_vector const *v){
	point tmp;
	for (int i = 0 ; i < 3 ; ++i)
		tmp[i] = gsl_vector_get(v, i);

	return tmp;
}

double CostFunc(gsl_vector const *v, void* params){

	OptParams* Params = (OptParams*)params;

	/*
	 * Get the lengths between the point and it's neighboring points
	 */
	vector<double> TAreas;
	TAreas.reserve(Params->AdjListCounts[Params->OptPoint]);
	//TAreas.reserve(TriAdjList[OptPoint].size());

	point CurPt = GSLVecToPoint(v);

	for (int PtNum : Params->AdjList[Params->OptPoint]){
		if (PtNum >= 0)
			TAreas.push_back(Params->p[PtNum].magsqr(CurPt));
	}

		/*
		 * get the areas of the local triangles (actually, the
		 * area squared, times two)
		 */
// 	for (int TNum : Params->TriAdjList[Params->OptPoint]){
// 
// 		if (TNum >= 0){
// 			double A, B, C;
// 
// 			if (Params->t[TNum].n1 == Params->OptPoint){
// 				A = CurPt.z * (Params->p[Params->t[TNum].n2].y - Params->p[Params->t[TNum].n3].y)
// 					+ Params->p[Params->t[TNum].n2].z * (Params->p[Params->t[TNum].n3].y - CurPt.y)
// 					+ Params->p[Params->t[TNum].n3].z * (CurPt.y - Params->p[Params->t[TNum].n2].y);
// 
// 				B = CurPt.x * (Params->p[Params->t[TNum].n2].z - Params->p[Params->t[TNum].n3].z)
// 					+ Params->p[Params->t[TNum].n2].x * (Params->p[Params->t[TNum].n3].z - CurPt.z)
// 					+ Params->p[Params->t[TNum].n3].x * (CurPt.z - Params->p[Params->t[TNum].n2].z);
// 
// 				C = CurPt.y * (Params->p[Params->t[TNum].n2].x - Params->p[Params->t[TNum].n3].x)
// 					+ Params->p[Params->t[TNum].n2].y * (Params->p[Params->t[TNum].n3].x - CurPt.x)
// 					+ Params->p[Params->t[TNum].n3].y * (CurPt.x - Params->p[Params->t[TNum].n2].x);
// 			}
// 			else if (Params->t[TNum].n2 == Params->OptPoint){
// 				A = Params->p[Params->t[TNum].n1].z * (CurPt.y - Params->p[Params->t[TNum].n3].y)
// 					+ CurPt.z * (Params->p[Params->t[TNum].n3].y - Params->p[Params->t[TNum].n1].y)
// 					+ Params->p[Params->t[TNum].n3].z * (Params->p[Params->t[TNum].n1].y - CurPt.y);
// 
// 				B = Params->p[Params->t[TNum].n1].x * (CurPt.z - Params->p[Params->t[TNum].n3].z)
// 					+ CurPt.x * (Params->p[Params->t[TNum].n3].z - Params->p[Params->t[TNum].n1].z)
// 					+ Params->p[Params->t[TNum].n3].x * (Params->p[Params->t[TNum].n1].z - CurPt.z);
// 
// 				C = Params->p[Params->t[TNum].n1].y * (CurPt.x - Params->p[Params->t[TNum].n3].x)
// 					+ CurPt.y * (Params->p[Params->t[TNum].n3].x - Params->p[Params->t[TNum].n1].x)
// 					+ Params->p[Params->t[TNum].n3].y * (Params->p[Params->t[TNum].n1].x - CurPt.x);
// 			}
// 			else {
// 				A = Params->p[Params->t[TNum].n1].z * (Params->p[Params->t[TNum].n2].y - CurPt.y)
// 					+ Params->p[Params->t[TNum].n2].z * (CurPt.y - Params->p[Params->t[TNum].n1].y)
// 					+ CurPt.z * (Params->p[Params->t[TNum].n1].y - Params->p[Params->t[TNum].n2].y);
// 
// 				B = Params->p[Params->t[TNum].n1].x * (Params->p[Params->t[TNum].n2].z - CurPt.z)
// 					+ Params->p[Params->t[TNum].n2].x * (CurPt.z - Params->p[Params->t[TNum].n1].z)
// 					+ CurPt.x * (Params->p[Params->t[TNum].n1].z - Params->p[Params->t[TNum].n2].z);
// 
// 				C = Params->p[Params->t[TNum].n1].y * (Params->p[Params->t[TNum].n2].x - CurPt.x)
// 					+ Params->p[Params->t[TNum].n2].y * (CurPt.x - Params->p[Params->t[TNum].n1].x)
// 					+ CurPt.y * (Params->p[Params->t[TNum].n1].x - Params->p[Params->t[TNum].n2].x);
// 			}
// 
// 			//TAreas.push_back(0.5 * sqrt(A*A + B*B + C*C));
// 			TAreas.push_back(A*A + B*B + C*C);
// 		}
// 	}

	/*
	 * get the average areas
	 */
	double Mean = 0;
	for (double const & i : TAreas)
		Mean += i;
	Mean /= (double)TAreas.size();

	/*
	 * get the variance of the areas of the triangles
	 */
	double Var = 0;
	for (double const & i : TAreas)
		Var += (i - Mean) * (i - Mean);
	Var /= (double)TAreas.size();

	return Var;
}

double TriAreaVar(OptParams * Params){
	double MaxVar = 0.0;
	vector<double> MaxList(Params->NumPts);

	omp_set_num_threads(Params->numCPU);

#pragma omp parallel for
	for (int PtNum = 0; PtNum < Params->NumPts; ++PtNum){
		vector<double> TAreas;
		TAreas.reserve(Params->TriAdjListCounts[PtNum]);
		for (int TriNum : Params->TriAdjList[PtNum]){
			if (TriNum >= 0){
				double A, B, C;

				A = Params->p[Params->t[TriNum].n1].z*(Params->p[Params->t[TriNum].n2].y - Params->p[Params->t[TriNum].n3].y)
					+ Params->p[Params->t[TriNum].n2].z*(Params->p[Params->t[TriNum].n3].y - Params->p[Params->t[TriNum].n1].y)
					+ Params->p[Params->t[TriNum].n3].z*(Params->p[Params->t[TriNum].n1].y - Params->p[Params->t[TriNum].n2].y);

				B = Params->p[Params->t[TriNum].n1].x*(Params->p[Params->t[TriNum].n2].z - Params->p[Params->t[TriNum].n3].z)
					+ Params->p[Params->t[TriNum].n2].x*(Params->p[Params->t[TriNum].n3].z - Params->p[Params->t[TriNum].n1].z)
					+ Params->p[Params->t[TriNum].n3].x*(Params->p[Params->t[TriNum].n1].z - Params->p[Params->t[TriNum].n2].z);

				C = Params->p[Params->t[TriNum].n1].y*(Params->p[Params->t[TriNum].n2].x - Params->p[Params->t[TriNum].n3].x)
					+ Params->p[Params->t[TriNum].n2].y*(Params->p[Params->t[TriNum].n3].x - Params->p[Params->t[TriNum].n1].x)
					+ Params->p[Params->t[TriNum].n3].y*(Params->p[Params->t[TriNum].n1].x - Params->p[Params->t[TriNum].n2].x);

				TAreas.push_back(A*A + B*B + C*C);
			}
		}

		/*
		 * get the average areas
		 */
		double Mean = 0;
		for (double const & i : TAreas)
			Mean += i;
		Mean /= (double)TAreas.size();

		/*
		 * get the variace of the areas of the triangles
		 */
		double Var = 0;
		for (double const & i : TAreas)
			Var += pow(i - Mean, 2);
		Var /= (double)TAreas.size();
		MaxList[PtNum] = Var;
	}

	for (int i = 0; i < Params->NumPts; ++i)
		if (MaxList[i] > MaxVar)
			MaxVar = MaxList[i];

	return MaxVar;
}

MeshStatus_e meshgen2D_sphere(double Radius,
	int Level, 
	vector<point> & ConstrainedVertices,
	vector<int> & MovedPointNums,
	point *& pIn, 
	triangle *& tIn,
	int **& eIn,
	int &NumPtsIn, 
	int &NumTriIn,
	int &NumEdgesIn)
{  

	//int** e;
	//int NumEdges;

	MeshStatus_e MeshStatus = SUCCESS;

	OptParams Params;

#ifdef WIN32
	SYSTEM_INFO sysinfo;
	GetSystemInfo(&sysinfo);

	Params.numCPU = sysinfo.dwNumberOfProcessors;
#else
	numCPU = 8;
#endif


	unsigned int MaxIter = 1000;
	double Tol = 1e-10;

	MeshCreateIcosa(Radius, Params.NumPts, Params.NumTri, Params.p, Params.t);

  // Subdivide mesh "level" number of times
  if (Level>0)
	{
		MeshSubdivide(Level, Radius, Params.NumPts, Params.NumTri, Params.p, Params.t);
	  Tol /= pow(10.0, Level);
	}

  /*
   * Tim:
   * Here is where I add the procedure to contrain vertices
   * and optimize the sphere mesh based on the constraints.
   */

  /*
   * Adding constrained vertices to system by moving point closest to
   * constraint to the constrained point location
   */
  MovedPointNums.clear();
  for (point cp : ConstrainedVertices){
	  double MinLenSqr = Radius * 10.0;
	  int MovedPointNum = -1;
	  for (int i = 0; i < Params.NumPts; ++i){
		  double TmpMin = Params.p[i].magsqr(cp);
		  if (TmpMin < MinLenSqr){
			  MinLenSqr = TmpMin;
			  MovedPointNum = i;
		  }
	  }
	  if (MovedPointNum >= 0){
		  bool NewPoint = true;
		  for (int pn : MovedPointNums)
			  if (MovedPointNum == pn){
				  NewPoint = false;
				  break;
			  }
		  if (NewPoint){
			  Params.p[MovedPointNum] = cp;
			  MovedPointNums.push_back(MovedPointNum);
		  }
		  else{
			  cerr << "Point " << MovedPointNum << " was moved more than once!";
			  return FAIL_NEEDREFINEMENT;
		  }
	  }
	  else{
		  cerr << "Min distance not found for constraint {" << cp.x << ", " << cp.y << ", " << cp.z << "}";
		  return FAIL_INVALIDCONSTRAINT;
	  }
  }

  /*
   * now perform the optimization:
   * 	iterate over each non-constrained point of the sphere mesh
   * 	minimizing a cost function based on local triangle area
   * 	variance. (really, we're using the variance of 2*area^2,
   * 	since we're only after having the same area and don't care
   * 	what the area is, this saves us some sqrt calculations.)
   */

	/*
	 * make the gsl error function pointer for the minimizer
	 */
	gsl_multimin_function ErrorFunc;
	ErrorFunc.n = 3;
	ErrorFunc.f = &CostFunc;
	ErrorFunc.params = (void*)&Params;

	/*
	 *	populate the point-triangle adjacency list
	 */

	Params.AdjList.resize(Params.NumPts, vector<int>(6, -1));
	Params.TriAdjList.resize(Params.NumPts, vector<int>(6, -1));
	Params.AdjListCounts.resize(Params.NumPts, 0);
	Params.TriAdjListCounts.resize(Params.NumPts, 0);
	for (int PtNum = 0; PtNum < Params.NumPts; ++PtNum){
		for (int ti = 0; ti < Params.NumTri; ++ti){
			for (int i = 0; i < 3; ++i){
				if (Params.t[ti][i] == PtNum){
					Params.TriAdjList[PtNum][Params.TriAdjListCounts[PtNum]++] = ti;
					for (int j = 0; j < 3; ++j){
						if (Params.t[ti][j] != PtNum){
							bool IsFound = false;
							for (int k = 0; !IsFound && k < Params.AdjListCounts[PtNum]; ++k){
								IsFound = (Params.AdjList[PtNum][k] == Params.t[ti][j]);
							}
							if (!IsFound)
								Params.AdjList[PtNum][Params.AdjListCounts[PtNum]++] = Params.t[ti][j];
						}
					}
				}
			}
		}
	}

	/*
	 *	Check that no two constrained points share a neighbor
	 */
	bool MeshFail = false;
	for (int ChkPt : MovedPointNums){
		for (int ChkNeighbor : Params.AdjList[ChkPt]){
			for (int ChkPt2 : MovedPointNums){
				if (ChkPt != ChkPt2){
					for (int ChkNeighbor2 : Params.AdjList[ChkPt2]){
						if (!MeshFail) MeshFail = (ChkNeighbor == ChkNeighbor2 && ChkNeighbor >= 0);
						if (MeshFail)
							int b = 1;
					}
				}
			}
		}
		if (MeshFail) break;
	}
	if (MeshFail){
		Params.AdjList.clear();
		Params.TriAdjList.clear();
		Params.AdjListCounts.clear();
		Params.TriAdjListCounts.clear();
		return FAIL_NEEDREFINEMENT;
	}

	

	if (ConstrainedVertices.size() > 0){

		vector<bool> PointIsConstrained(Params.NumPts, false);
		for (int i = 0; i < ConstrainedVertices.size(); ++i)
			PointIsConstrained[MovedPointNums[i]] = true;

		/*
		 *	First, do several iterations of jiggle-mesh,
		 *	moving each node to the midpoint of its
		 *	neighborhood
		 */

		for (int i = 0; i < 50; ++i){
			for (int PtNum = 0; PtNum < Params.NumPts; ++PtNum){
				if (!PointIsConstrained[PtNum]){

					double NewPoint[3] = { 0 };

					for (int i : Params.AdjList[PtNum]){
						if (i >= 0){
							for (int k = 0; k < 3; ++k){
								NewPoint[k] += Params.p[i][k];
							}
						}
					}
					for (int i = 0; i < 3; ++i){
						NewPoint[i] /= double(Params.AdjListCounts[PtNum]);
						Params.p[PtNum][i] = NewPoint[i];
					}

					/*
					* Project new point to sphere surface
					*/
					// Compute spherical coordinates
					double   rad = sqrt(pow(Params.p[PtNum].x, 2) + pow(Params.p[PtNum].y, 2) + pow(Params.p[PtNum].z, 2));
					double theta = acos(Params.p[PtNum].z / rad);
					double   phi = atan2(Params.p[PtNum].y, Params.p[PtNum].x);

					// Project point onto a sphere of radius "radius"
					Params.p[PtNum].x = Radius * sin(theta) * cos(phi);
					Params.p[PtNum].y = Radius * sin(theta) * sin(phi);
					Params.p[PtNum].z = Radius * cos(theta);
				}
			}
		}


		double OldVar = 0.0;
		double NewVar = TriAreaVar(&Params);
		double DVar = fabs(NewVar - OldVar);

		unsigned int GIter = 0;
		while (DVar >= Tol && GIter < MaxIter){
			GIter++;
			for (Params.OptPoint = 0; Params.OptPoint < Params.NumPts; ++Params.OptPoint){
				/*
				 * only optimize points that aren't contrained
				 */
				bool IsConstrained = false;

				if (!PointIsConstrained[Params.OptPoint]){

					/*
					 * need a step size and initial guess for the minimizer.
					 * arbitrarily pick a fraction of the average
					 * distance between the vertices of the triangles
					 * in the vertex neighborhood for step size,
					 * and pick the midpoint of the neighborhood for
					 * the guess.
					 */

					double DStep = 0;
					double G[3] = { 0 };
					for (int i : Params.AdjList[Params.OptPoint]){
						if (i >= 0){
							for (int k = 0; k < 3; ++k){
								G[k] += Params.p[i][k];
								DStep += pow(Params.p[i][k] - Params.p[Params.OptPoint][k], 2);
							}
						}
					}
					for (int i = 0; i < 3; ++i)
						G[i] /= double(Params.AdjListCounts[Params.OptPoint]);

					DStep = sqrt(DStep) / double(Params.AdjListCounts[Params.OptPoint]) * 0.4;

					gsl_vector* Step = gsl_vector_alloc(3);
					gsl_vector* Guess = gsl_vector_alloc(3);
					for (int i = 0; i < 3; ++i){
						gsl_vector_set(Step, i, DStep);
						gsl_vector_set(Guess, i, G[i]);
					}

					/*
					 * initialize the minimizer
					 */
					gsl_multimin_fminimizer_type const *T;
					T = gsl_multimin_fminimizer_nmsimplex2;

					gsl_multimin_fminimizer *s;
					s = gsl_multimin_fminimizer_alloc(T, 3);

					gsl_multimin_fminimizer_set(s, &ErrorFunc, Guess, Step);

					/*
					 * prepare and start iterations
					 */
					unsigned int Iter = 0;
					int Status = GSL_CONTINUE;

					while (Status == GSL_CONTINUE && Iter < MaxIter){
						++Iter;

						Status = gsl_multimin_fminimizer_iterate(s);

					}

					if (Status != GSL_SUCCESS || Iter >= MaxIter){
						cerr << "GSL minimizer failed on point " << Params.OptPoint << "\n";
					}

					gsl_vector* NewPoint = gsl_multimin_fminimizer_x(s);
					Params.p[Params.OptPoint] = GSLVecToPoint(NewPoint);

					/*
					 * Project new point to sphere surface
					 */
					// Compute spherical coordinates
					double   rad = sqrt(pow(Params.p[Params.OptPoint].x, 2) + pow(Params.p[Params.OptPoint].y, 2) + pow(Params.p[Params.OptPoint].z, 2));
					double theta = acos(Params.p[Params.OptPoint].z / rad);
					double   phi = atan2(Params.p[Params.OptPoint].y, Params.p[Params.OptPoint].x);

					// Project point onto a sphere of radius "radius"
					Params.p[Params.OptPoint].x = Radius * sin(theta) * cos(phi);
					Params.p[Params.OptPoint].y = Radius * sin(theta) * sin(phi);
					Params.p[Params.OptPoint].z = Radius * cos(theta);


					gsl_multimin_fminimizer_free(s);
					gsl_vector_free(Step);
					gsl_vector_free(Guess);
				}
			}
			OldVar = NewVar;
			NewVar = TriAreaVar(&Params);
			DVar = fabs(NewVar - OldVar);
			//		if (GIter > 29)
			//			break;
		}

		if (GIter == MaxIter && DVar >= Tol){
			MeshStatus = SUCCESS_NOTCONVERGED;
		}

		cout << GIter << " iterations total\n" << Tol << " > " << DVar << "\n";

	}

  // Compute element areas and dual-element areas
	double* area = new double[Params.NumTri];
	double* cdual = new double[Params.NumPts];
	MeshComputeAreas(Params.NumTri, Params.NumPts, Params.p, Params.t, area, cdual);

  // Store spherical coordinates of p in psphere
	point* psphere = new point[Params.NumPts];
	MeshSphereCoords(Params.NumPts, Params.p, psphere);

  // Create information about edges
  int numedges = 1;
  dTensor2* edge;
  iTensor2* tedge;
  iTensor2* eelem;
  int* bnd_node = new int[1];
  int* ghost_link = new int[1];
  edge = new dTensor2(numedges,4);
  tedge = new iTensor2(Params.NumTri, 3);
  eelem = new iTensor2(numedges,2);
  MeshEdgeData(Params.NumPts, Params.NumTri, Params.p, psphere, Params.t, area, numedges,
			   edge,tedge,eelem);

  int NumEdges = numedges;
  int** e;
  e = new int*[NumEdges];
  for (int i = 0; i < NumEdges; ++i){
	  e[i] = new int[2];
	  e[i][0] = e[i][1] = -1;
  }

  int EdgeNum = 0;
  for (int i = 0; i < Params.NumPts; ++i){
	  for (int j = 0; j < Params.AdjListCounts[i]; ++j){
		  bool IsFound = false;
		  for (int k = 0; k < NumEdges; ++k){
			  if (e[k][0] > 0 || e[k][1] > 0)
			  {
				  if ((e[k][0] == i && e[k][1] == Params.AdjList[i][j])
					  || (e[k][0] == Params.AdjList[i][j] && e[k][1] == i)){
					  IsFound = true;
					  break;
				  }
			  }
		  }
		  if (!IsFound){
			  e[EdgeNum][0] = i;
			  e[EdgeNum][1] = Params.AdjList[i][j];
			  EdgeNum++;
		  }
	  }
  }

  // Store all mesh info 
  mesh Mesh(Params.NumTri, Params.NumTri, Params.NumPts, Params.NumPts, 0, numedges, 0);
  MeshStore(psphere, Params.t, bnd_node, ghost_link, area, cdual, edge, tedge, eelem, Mesh);

  // Adjust edge information so that the element to the "left"
  // of the edge is always "upwind" of the unit normal to the
  // edge and the element to the "right" of the edge is always
  // "downwind" of the unit normal to the edge
  MeshOrientEdge(Mesh);

  // Output grid
//   ScreenOutput(Mesh);
  //Mesh.OutputMesh(outputdir);
  
//   delete t,p;

  pIn = Params.p;
  tIn = Params.t;
  eIn = e;

  Params.AdjList.clear();
  Params.TriAdjList.clear();
  Params.AdjListCounts.clear();
  Params.TriAdjListCounts.clear();

  NumPtsIn = Params.NumPts;
  NumTriIn = Params.NumTri;
  NumEdgesIn = numedges;

  return MeshStatus;
}
