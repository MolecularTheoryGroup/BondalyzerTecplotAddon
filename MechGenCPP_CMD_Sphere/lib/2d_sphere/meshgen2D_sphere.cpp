#include <vector>
#include <cmath>
#include <iostream>
#include <sstream>

#include <gsl/gsl_math.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_vector.h>

#include "mesh.h"
#include "constants.h"

#include "parsearguments.h"
#include "screenoutput.h"

#define DEBUG 0

using std::vector;
using std::cout;
using std::cerr;
using std::endl;

point* p;
triangle* t;
int numpts, numtri, OptPoint;
std::vector<int> TList;
std::vector<double> TAreas;

unsigned int MaxIter = 1000;
double Tol = 1e-9;

point GSLVecToPoint(const gsl_vector *v){
	point tmp;
	for (int i = 0 ; i < 3 ; ++i)
		tmp[i] = gsl_vector_get(v, i);

	return tmp;
}

double CostFunc(const gsl_vector *v, void* params){

	/*
	 * Get the lengths between the point and it's neighboring points
	 */

	vector<int> UsedPoints;
	UsedPoints.reserve(5);
	TAreas.clear();
	TAreas.reserve(6);

	point CurPt = GSLVecToPoint(v);

	for (vector<int>::iterator TNum = TList.begin() ; TNum != TList.end() ; TNum++){
		for (int TNode = 0 ; TNode < 3 ; ++TNode){
			if (t[*TNum][TNode] != OptPoint){
				bool NewPoint = true;
				for (int i = 0 ; i < UsedPoints.size() && NewPoint ; ++i)
					if (UsedPoints[i] == t[*TNum][TNode])
						NewPoint = false;

				if (NewPoint){
					UsedPoints.push_back(t[*TNum][TNode]);
					TAreas.push_back(p[t[*TNum][TNode]].magsqr(CurPt));
//					TAreas.push_back(p[t[*TNum][TNode]].magsqr(p[OptPoint]));
				}
			}
		}

		/*
		 * get the areas of the local triangles (actually, the
		 * area squared, times two)
		 */

//		double A,B,C;
//
//		if (t[*TNum].n1 == OptPoint){
//			A = CurPt.z * (p[t[*TNum].n2].y - p[t[*TNum].n3].y)
//					+ p[t[*TNum].n2].z * (p[t[*TNum].n3].y - CurPt.y)
//					+ p[t[*TNum].n3].z * (CurPt.y - p[t[*TNum].n2].y);
//
//			B = CurPt.x * (p[t[*TNum].n2].z - p[t[*TNum].n3].z)
//					+ p[t[*TNum].n2].x * (p[t[*TNum].n3].z - CurPt.z)
//					+ p[t[*TNum].n3].x * (CurPt.z - p[t[*TNum].n2].z);
//
//			C = CurPt.y * (p[t[*TNum].n2].x - p[t[*TNum].n3].x)
//					+ p[t[*TNum].n2].y * (p[t[*TNum].n3].x - CurPt.x)
//					+ p[t[*TNum].n3].y * (CurPt.x - p[t[*TNum].n2].x);
//		}
//		else if (t[*TNum].n2 == OptPoint){
//			A = p[t[*TNum].n1].z * (CurPt.y - p[t[*TNum].n3].y)
//					+ CurPt.z * (p[t[*TNum].n3].y - p[t[*TNum].n1].y)
//					+ p[t[*TNum].n3].z * (p[t[*TNum].n1].y - CurPt.y);
//
//			B = p[t[*TNum].n1].x * (CurPt.z - p[t[*TNum].n3].z)
//					+ CurPt.x * (p[t[*TNum].n3].z - p[t[*TNum].n1].z)
//					+ p[t[*TNum].n3].x * (p[t[*TNum].n1].z - CurPt.z);
//
//			C = p[t[*TNum].n1].y * (CurPt.x - p[t[*TNum].n3].x)
//					+ CurPt.y * (p[t[*TNum].n3].x - p[t[*TNum].n1].x)
//					+ p[t[*TNum].n3].y * (p[t[*TNum].n1].x - CurPt.x);
//		}
//		else {
//			A = p[t[*TNum].n1].z * (p[t[*TNum].n2].y - CurPt.y)
//					+ p[t[*TNum].n2].z * (CurPt.y - p[t[*TNum].n1].y)
//					+ CurPt.z * (p[t[*TNum].n1].y - p[t[*TNum].n2].y);
//
//			B = p[t[*TNum].n1].x * (p[t[*TNum].n2].z - CurPt.z)
//					+ p[t[*TNum].n2].x * (CurPt.z - p[t[*TNum].n1].z)
//					+ CurPt.x * (p[t[*TNum].n1].z - p[t[*TNum].n2].z);
//
//			C = p[t[*TNum].n1].y * (p[t[*TNum].n2].x - CurPt.x)
//					+ p[t[*TNum].n2].y * (CurPt.x - p[t[*TNum].n1].x)
//					+ CurPt.y * (p[t[*TNum].n1].x - p[t[*TNum].n2].x);
//		}
//
//		TAreas.push_back( 0.5 * sqrt(A*A + B*B + C*C) );
//		TAreas.push_back( A*A + B*B + C*C );
	}

	/*
	 * get the average areas
	 */
	double Mean = 0;
	for (vector<double>::iterator i = TAreas.begin() ; i != TAreas.end() ; i++)
		Mean += *i;
	Mean /= (double)TAreas.size();

	/*
	 * get the variace of the areas of the triangles
	 */
	double Var = 0;
	for (vector<double>::iterator i = TAreas.begin() ; i != TAreas.end() ; i++)
		Var += pow(*i - Mean, 2);
	Var /= (double)TAreas.size();

	return Var;
}

double TriAreaVar(){
	TAreas.clear();
	TAreas.reserve(numtri);

	for (int i = 0 ; i < numtri ; ++i){
		double A,B,C;

		A = p[t[i].n1].z*(p[t[i].n2].y-p[t[i].n3].y)+p[t[i].n2].z*(p[t[i].n3].y-p[t[i].n1].y)
				+p[t[i].n3].z*(p[t[i].n1].y-p[t[i].n2].y);

		B = p[t[i].n1].x*(p[t[i].n2].z-p[t[i].n3].z)+p[t[i].n2].x*(p[t[i].n3].z-p[t[i].n1].z)
				+p[t[i].n3].x*(p[t[i].n1].z-p[t[i].n2].z);

		C = p[t[i].n1].y*(p[t[i].n2].x-p[t[i].n3].x)+p[t[i].n2].y*(p[t[i].n3].x-p[t[i].n1].x)
				+p[t[i].n3].y*(p[t[i].n1].x-p[t[i].n2].x);

		TAreas.push_back( A*A + B*B + C*C );
	}

	/*
	 * get the average areas
	 */
	double Mean = 0;
	for (vector<double>::iterator i = TAreas.begin() ; i != TAreas.end() ; i++)
		Mean += *i;
	Mean /= (double)TAreas.size();

	/*
	 * get the variace of the areas of the triangles
	 */
	double Var = 0;
	for (vector<double>::iterator i = TAreas.begin() ; i != TAreas.end() ; i++)
		Var += pow(*i - Mean, 2);
	Var /= (double)TAreas.size();

	return Var;
}

int meshgen2D_sphere(int argc, char**argv)
{  
  // Output title information
  printf("\n");
  printf("   -------------------------------------------------   \n");
  printf("   |  MeshGenC++: An Unstructured Grid Generator   |   \n");
  printf("   |    C++ Version of Per-Olof Persson and        |   \n");
  printf("   |    Gilbert Strang's DistMesh algorithm for    |   \n");
  printf("   |    unstructured grid generation               |   \n");
  printf("   |                                               |   \n");
  printf("   |                   Written by                  |   \n");
  printf("   |    James A. Rossmanith and the DoGPack Team   |   \n");
  printf("   -------------------------------------------------   \n");
  printf("\n");

//  // Needed functions
//  void ParseArguments(int argc,
//		      char**argv,
//		      char*& outputdir);
//  void MeshInputData(double&, int&);
//  void MeshCreateIcosa(const double radius, int& numpts, int& numtri,
//		       point*& p, triangle*& t);
//  void MeshSubdivide(const int level, const double radius,
//		     int& numpts, int& numtri, point*& p, triangle*& t);
//  void MeshSphereCoords(const int, const point* p, point*& psphere);
//  void MeshEdgeData(const int& numpts, const int& numtri, const point p[],
//		    const point psphere[], const triangle t[], const double area[],
//		    int& numedges, dTensor2*& edge, iTensor2*& tedge,
//		    iTensor2*& eelem);
//  void MeshOrientEdge(mesh& Mesh);
//  void MeshComputeAreas(const int numtri, const int numpts,
//			const point p[], triangle t[],
//			double*& area, double*& cdual);
//  void ScreenOutput(const mesh& Mesh);
//  void MeshStore(point p[], triangle t[], int bnd_node[],
//		 int ghost_link[], double area[], double cdual[],
//		 dTensor2*& edge, iTensor2*& tedge,
//		 iTensor2*& eelem, mesh& Mesh);
  
  // Parse arguments -- sets directory to which output will be sent,
  //                    the default is "output"  
  char* outputdir = new char[12];
  outputdir[0] ='m';
  outputdir[1] ='e';
  outputdir[2] ='s';
  outputdir[3] ='h';
  outputdir[4] ='_';
  outputdir[5] ='o';
  outputdir[6] ='u';
  outputdir[7] ='t';
  outputdir[8] ='p';
  outputdir[9] ='u';
  outputdir[10]='t';
  outputdir[11]='\0';
  ParseArguments(argc,argv,outputdir);

  // Get data from input file
  double radius;
  int level;
  MeshInputData(radius,level);

  // Create icosahedral mesh
//  point* p;
//  triangle* t;
//  int numpts,numtri;
  MeshCreateIcosa(radius,numpts,numtri,p,t);

  // Subdivide mesh "level" number of times
  if (level>0)
	{
	  MeshSubdivide(level,radius,numpts,numtri,p,t);
	  Tol /= pow(10.0, level);
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
  vector<point> ConstrainedVertices;
  vector<int> MovedPointNums;
  ConstrainedVertices.push_back(point(0.413554, -0.603113, -0.682076));
  ConstrainedVertices.push_back(point(0.825853, 0.474478, -0.304693));
  ConstrainedVertices.push_back(point(-0.0148872, 0.98687, 0.160826));
  ConstrainedVertices.push_back(point(-0.637447, -0.0119685, -0.770401));
  ConstrainedVertices.push_back(point(0.0315743, -0.0234456, 0.999226));
  ConstrainedVertices.push_back(point(-0.567848, -0.810362, 0.144434));
  ConstrainedVertices.push_back(point(0.614949, -0.772711, -0.157342));
  ConstrainedVertices.push_back(point(-0.138839, -0.28429, 0.948632));
  ConstrainedVertices.push_back(point(-0.875728, 0.247566, 0.414503));
  ConstrainedVertices.push_back(point(-0.296182, 0.280499, -0.913015));
  for (vector<point>::iterator cp = ConstrainedVertices.begin() ; cp != ConstrainedVertices.end() ; cp++){
	  double MinLenSqr = radius * 10.0;
	  int MovedPointNum = -1;
	  for (int i = 0 ; i < numpts ; ++i){
		  double TmpMin = p[i].magsqr(*cp);
		  if (TmpMin < MinLenSqr){
			  MinLenSqr = TmpMin;
			  MovedPointNum = i;
		  }
	  }
	  if (MovedPointNum >= 0){
		  bool NewPoint = true;
		  for (vector<int>::iterator pn = MovedPointNums.begin() ; pn != MovedPointNums.end() ; pn++)
			  if (MovedPointNum == *pn){
				  NewPoint = false;
				  break;
			  }
		  if (NewPoint){
			  p[MovedPointNum] = *cp;
			  MovedPointNums.push_back(MovedPointNum);
		  }
		  else{
			  cerr << "Point " << MovedPointNum << " was moved more than once!";
			  return 1;
		  }
	  }
	  else{
		  cerr << "Min distance not found for constraint {" << cp->x << ", " << cp->y << ", " << cp->z << "}";
		  return 1;
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
	ErrorFunc.params = (void*)NULL;

	double OldVar = 0.0;
	double NewVar = TriAreaVar();
	double DVar = fabs(NewVar - OldVar);

	unsigned int GIter = 0;
	while (DVar >= Tol && GIter < MaxIter){
		GIter++;
		for (OptPoint = 0 ; OptPoint < numpts ; ++OptPoint){
			/*
			 * only optimize points that aren't contrained
			 */
			bool IsConstrained = false;
#if (DEBUG)
			cout << endl;
#endif

			for (vector<int>::iterator i = MovedPointNums.begin() ; i != MovedPointNums.end() && !IsConstrained ; i++)
				if (*i == OptPoint)
					IsConstrained = true;

			if (!IsConstrained){
				/*
				 * first, get the list of triangles that share the point
				 */
				TList.clear();
				TList.reserve(5); //there are at least 5 triangles per vertex
				for (int ti = 0 ; ti < numtri ; ++ti){
					if (t[ti].n1 == OptPoint || t[ti].n2 == OptPoint || t[ti].n3 == OptPoint)
						TList.push_back(ti);
				}
				TAreas.clear();
				TAreas.reserve(TList.size());

				/*
				 * need a step size and initial guess for the minimizer.
				 * arbitrarily pick a fraction of the average
				 * distance between the vertices of the triangles
				 * in the vertex neighborhood or step size,
				 * and pick the midpoint of the neighborhood for
				 * the guess.
				 */

				double DStep = 0;
				double G[3] = {0};
				vector<int> PList;
#if (DEBUG)
			vector<std::stringstream> ss(3);
			ss[0] << "x: ";
			ss[1] << "y: ";
			ss[2] << "z: ";
#endif
				for (vector<int>::iterator i = TList.begin() ; i != TList.end() ; ++i)
					for (int j = 0 ; j < 3 ; ++j){
						if (t[*i][j] != OptPoint){
							bool InList = false;
							for (int k = 0 ; k < PList.size() && !InList ; ++k)
								if (PList[k] == t[*i][j])
									InList = true;
							if (!InList){
								PList.push_back(t[*i][j]);
								for (int k = 0 ; k < 3 ; ++k){
	#if (DEBUG)
			ss[k] << p[t[*i][j]][k] << ", ";

	#endif
									G[k] += p[t[*i][j]][k];
									DStep += pow(p[t[*i][j]][k] - p[OptPoint][k], 2);
								}
							}
						}
					}
				for (int i = 0 ; i < 3 ; ++i)
					G[i] /= double(PList.size());

#if (DEBUG)
		cout << "point " << OptPoint << " before: {" << p[OptPoint].x << ", " << p[OptPoint].y << ", " << p[OptPoint].z << "}\n";
		for (int i = 0 ; i < 3 ; ++i)
			cout << ss[i].str() << endl;
		cout << "midpoint: ";
		for (int i = 0 ; i < 3 ; ++i)
			cout << G[i] << ", ";
		cout << endl;
#endif

//				for (int i = 0 ; i < 3 ; ++i)
//					p[OptPoint][i] = G[i];

				DStep = sqrt(DStep) / double(PList.size()) * 0.1;

				gsl_vector* Step = gsl_vector_alloc(3);
				gsl_vector* Guess = gsl_vector_alloc(3);
				for (int i = 0 ; i < 3 ; ++i){
					gsl_vector_set(Step, i, DStep);
//					gsl_vector_set(Guess, i, p[OptPoint][i]);
					gsl_vector_set(Guess, i, G[i]);
				}

				/*
				 * initialize the minimizer
				 */
				const gsl_multimin_fminimizer_type *T;
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

	#if (DEBUG)
			cout << "iteration " << Iter << ": Var = " << CostFunc(gsl_multimin_fminimizer_x(s), (void*) NULL) << "\n";
	#endif
				}

				if (Status != GSL_SUCCESS || Iter >= MaxIter){
					cerr << "GSL minimizer failed on point " << OptPoint << "\n";
				}

				gsl_vector* NewPoint = gsl_multimin_fminimizer_x(s);
				p[OptPoint] = GSLVecToPoint(NewPoint);
//				for (int i = 0 ; i < 3 ; ++i)
//					p[OptPoint][i] = gsl_vector_get(NewPoint, i);

				/*
				 * Project new point to sphere surface
				 */
				// Compute spherical coordinates
				double   rad = sqrt(pow(p[OptPoint].x,2)+pow(p[OptPoint].y,2)+pow(p[OptPoint].z,2));
				double theta = acos(p[OptPoint].z/rad);
				double   phi = atan2(p[OptPoint].y,p[OptPoint].x);

				// Project point onto a sphere of radius "radius"
				p[OptPoint].x = radius * sin(theta) * cos(phi);
				p[OptPoint].y = radius * sin(theta) * sin(phi);
				p[OptPoint].z = radius * cos(theta);

	#if (DEBUG)
			cout << "after: {" << p[OptPoint].x << ", " << p[OptPoint].y << ", " << p[OptPoint].z << "}\n";
	#endif

				gsl_multimin_fminimizer_free(s);
				gsl_vector_free(Step);
				gsl_vector_free(Guess);
			}
		}
		OldVar = NewVar;
		NewVar = TriAreaVar();
		DVar = fabs(NewVar - OldVar);
//		if (GIter > 29)
//			break;
	}

	cout << GIter << " iterations total\n" << Tol << " > " << DVar << "\n";



  // Compute element areas and dual-element areas
  double* area = new double[numtri];
  double* cdual = new double[numpts];
  MeshComputeAreas(numtri,numpts,p,t,area,cdual);

  // Store spherical coordinates of p in psphere
  point* psphere = new point[numpts];
  MeshSphereCoords(numpts,p,psphere);

  // Create information about edges
  int numedges = 1;
  dTensor2* edge;
  iTensor2* tedge;
  iTensor2* eelem;
  int* bnd_node = new int[1];
  int* ghost_link = new int[1];
  edge = new dTensor2(numedges,4);
  tedge = new iTensor2(numtri,3);
  eelem = new iTensor2(numedges,2);
  MeshEdgeData(numpts,numtri,p,psphere,t,area,numedges,
			   edge,tedge,eelem);

  // Store all mesh info 
  mesh Mesh(numtri,numtri,numpts,numpts,0,numedges,0);
  MeshStore(psphere,t,bnd_node,ghost_link,area,cdual,edge,tedge,eelem,Mesh);

  // Adjust edge information so that the element to the "left"
  // of the edge is always "upwind" of the unit normal to the
  // edge and the element to the "right" of the edge is always
  // "downwind" of the unit normal to the edge
  MeshOrientEdge(Mesh);

  // Output grid
  ScreenOutput(Mesh);
  Mesh.OutputMesh(outputdir);
  
  delete t,p;

  return 0;
}
