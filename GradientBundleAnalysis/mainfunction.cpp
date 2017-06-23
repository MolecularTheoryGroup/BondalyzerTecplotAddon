
#include "TECADDON.h"
#include "ADDGLBL.h"
#include "GUIDEFS.h"
#if !defined (MSWIN)
#include <unistd.h>
#endif

//#include <cstdlib>
#include <sstream>
#include <string>
#include <vector>

// for profiling
#include <chrono>
#include <ctime>
#include <ratio>

#include <omp.h>

#include "ZONEVARINFO.h"
#include "CSM_DATA_TYPES.h"
#include "meshgen2d_sphere.h"
#include "STREAMTRACEOPERATIONS.h"
#include "VIEWRESULTS.h"

#include "MAINFUNCTION.h"

#include <armadillo>
using namespace arma;

using std::string;
using std::to_string;
using std::stoi;
using std::stringstream;
using std::vector;

//for profiling
using std::chrono::high_resolution_clock;
using std::chrono::duration;
using std::chrono::duration_cast;





enum RadMode{
	ABSOLUTERADIUS = 1,
	MINCPDISTRATIO
};

void MainFunction(){

	SYSTEM_INFO sysinfo;
	GetSystemInfo(&sysinfo);

	int numCPU = sysinfo.dwNumberOfProcessors;
	omp_set_num_threads(numCPU);

	int MemoryRequired;

	TecUtilLockStart(AddOnID);
	TecUtilDrawGraphics(FALSE);
	TecUtilInterfaceSuspend(TRUE);

	high_resolution_clock::time_point Time1 = high_resolution_clock::now();
	high_resolution_clock::time_point Time2;

	Boolean_t IsOk = TRUE;


// 	high_resolution_clock::time_point Time1 = high_resolution_clock::now();
// 	high_resolution_clock::time_point Time2;

	EntIndex_t NumVars = TecUtilDataSetGetNumVars();

	int GroupNum = 0;

	/*
	 *	Get job parameters from dialog
	 */
	Boolean_t UseCutoff = TRUE;
	UseCutoff = TecGUIToggleGet(TGLOpenSys_TOG_D1);
	double CutoffVal = 0.001;
	TecGUITextFieldGetDouble(TFCutoff_TF_D1, &CutoffVal);
	EntIndex_t CutoffVarNum = VarNumByName(string("Electron Density"));

	/*
	 *	When checking if two streamtraces are straddling two IBs,
	 *	use the directions of their ends and the distance between
	 *	their ends to check.
	 *	IBCheckDistRatio is the fraction of the sphere radius that,
	 *	if the streamtrace endpoints are closer than,
	 *	they are considered to be at the same point and not in
	 *	different IBs.
	 *	IBCheckAngle is the angle, in degrees, that if greater than the
	 *	angle between the end segments of the streamtraces, they are
	 *	considered to be parallel and considered in the same IB.
	 */
	double IBCheckDistRatio = 0.2;
	TecGUITextFieldGetDouble(TFIBDist_TF_D1, &IBCheckDistRatio);
	double IBCheckAngle = 20.0;
	TecGUITextFieldGetDouble(TFIBAng_TF_D1, &IBCheckAngle);

	LgIndex_t * CPNums = NULL;
	LgIndex_t NumSelectedCPs = 1;
	TecGUIListGetSelectedItems(MLSelCPs_MLST_D1, &CPNums, &NumSelectedCPs);
	/*
	*	Set sphere radius and mesh refinement level
	*	(will be set by user in future)
	*/


	RadMode RadiusMode = MINCPDISTRATIO;
	RadiusMode = (RadMode)TecGUIRadioBoxGetToggle(RBRadMode_RADIO_D1);
	double UserRadius = 0.25;
	TecGUITextFieldGetDouble(TFRad_TF_D1, &UserRadius);
	int Level = 3;
	TecGUITextFieldGetLgIndex(TFLevel_TFS_D1, &Level);
	LgIndex_t NumSTPoints = 100;
	TecGUITextFieldGetLgIndex(TFSTPts_TF_D1, &NumSTPoints);


	EntIndex_t XYZVarNums[3];
	FieldData_pa CPXYZPtrs[3];
	EntIndex_t CPZoneNum = ZoneNumByName(string("Critical Points"));

	EntIndex_t OldNumZones = TecUtilDataSetGetNumZones();

	EntIndex_t CPTypeVarNum;
	FieldData_pa CPTypeVarPtr;

// 	vector<vector<LgIndex_t> > CPIndexList(4, vector<LgIndex_t>());

	LgIndex_t NumCPs[4];
	if (IsOk){
		AuxData_pa CPAuxData = TecUtilAuxDataZoneGetRef(CPZoneNum);
		if (VALID_REF(CPAuxData)){
			for (int i = 0; i < 4; ++i){
				char* TempCStr;
				AuxDataType_e ADTJunk;
				Boolean_t BJunk;
				IsOk = TecUtilAuxDataGetItemByName(CPAuxData, CCDataNumCPs[i].c_str(), (ArbParam_t*)&TempCStr, &ADTJunk, &BJunk);
				if (IsOk){
					NumCPs[i] = stoi(TempCStr);
				}
				TecUtilStringDealloc(&TempCStr);
			}
		}
	}

	LgIndex_t IJKMax[3];

	TecUtilZoneGetIJK(CPZoneNum, &IJKMax[0], &IJKMax[1], &IJKMax[2]);
	IsOk = (IJKMax[0] > 1 && IJKMax[1] == 1 && IJKMax[2] == 1);
	if (!IsOk){
		TecUtilDrawGraphics(TRUE);
		TecUtilInterfaceSuspend(FALSE);
		TecUtilDialogErrMsg("Critical point zone dimensions incorrect.");
		TecUtilLockFinish(AddOnID);
		return;
	}

	/*
	*	Get XYZ readable pointers to the crit points zone to get the
	*	position of the fifth atom (working with cubane, so this is the first carbon).
	*/

	if (CPZoneNum <= 0){
		TecUtilDrawGraphics(TRUE);
		TecUtilInterfaceSuspend(FALSE);
		TecUtilDialogErrMsg("Couldn't get critical points zone.");
		TecUtilLockFinish(AddOnID);
		return;
	}

	TecUtilAxisGetVarAssignments(&XYZVarNums[0], &XYZVarNums[1], &XYZVarNums[2]);
	for (int i = 0; i < 3 && IsOk; ++i)
		IsOk = (XYZVarNums[i] > 0);

	if (IsOk){
		for (int i = 0; i < 3 && IsOk; ++i){
			CPXYZPtrs[i] = TecUtilDataValueGetReadableRef(CPZoneNum, XYZVarNums[i]);
			IsOk = VALID_REF(CPXYZPtrs[i]);
		}
	}

	if (!IsOk){
		TecUtilDrawGraphics(TRUE);
		TecUtilInterfaceSuspend(FALSE);
		TecUtilDialogErrMsg("Couldn't get XYZ vars.");
		TecUtilLockFinish(AddOnID);
		return;
	}

	/*
	*	Get lists of nearfield ring and bond CP indices, and find
	*	the closest other CP
	*/
	if (IsOk){
		CPTypeVarNum = VarNumByName(string("CritPointType"));
		if (CPTypeVarNum > 0){
			CPTypeVarPtr = TecUtilDataValueGetReadableRef(CPZoneNum, CPTypeVarNum);
			IsOk = VALID_REF(CPTypeVarPtr);
		}
		else
			IsOk = FALSE;
	}

// 	int CPRanks[4] = { -3, -1, 1, 3 };
// 	for (int i = 1; i <= IJKMax[0] && IsOk; ++i){
// 		int CPType = TecUtilDataValueGetByRef(CPTypeVarPtr, i);
// 		for (int j = 0; j < 4; ++j){
// 			if (CPType == CPRanks[j]){
// 				CPIndexList[j].push_back(i);
// 			}
// 		}
// 	}

	for (int SelectCPNum = 0; SelectCPNum < NumSelectedCPs && IsOk; ++SelectCPNum){

		LgIndex_t CPNum = 43;

		double Radius = 0.25;

		stringstream ProgressStr;
		string CPName;

		if (IsOk){
			int RealCPNum = 0;
			string CPTypeStr;
			stringstream ss;
			char* TempCStr = TecGUIListGetString(MLSelCPs_MLST_D1, CPNums[SelectCPNum]);
			CPName = TempCStr;

			ProgressStr << "Processing " << CPName;
			if (NumSelectedCPs > 1){
				ProgressStr << " ... CP " << SelectCPNum + 1 << " of " << NumSelectedCPs;
			}

			ss << TempCStr;
			ss >> CPTypeStr >> CPNum;

			for (int i = 0; i < 4; ++i){
				if (CPTypeStr == RankStrs[i]){
					int RealCPNum = 0;
					for (int j = 0; j < i; ++j){
						RealCPNum += NumCPs[j];
					}
					CPNum = RealCPNum + CPNum;
					break;
				}
			}

			TecUtilStringDealloc(&TempCStr);
		}
		else{
			TecUtilDialogErrMsg("Didn't get CP number correctly");
		}

		if (IsOk){
			stringstream ss;
			ss << ProgressStr.str() << " ... Generating Mesh";
			TecUtilDrawGraphics(TRUE);
			TecUtilPleaseWait(ss.str().c_str(), TRUE);
			TecUtilDrawGraphics(FALSE);
		}

		int CPType;

		CPType = TecUtilDataValueGetByRef(CPTypeVarPtr, CPNum);
		vec3 CPPos;

		for (int i = 0; i < 3; ++i){
			CPPos[i] = TecUtilDataValueGetByRef(CPXYZPtrs[i], CPNum);
		}

		vector<point> IntersectionPoints;
		vector<vector<LgIndex_t> > IntZoneSaddleCPNodeNums;
		vector<vector<LgIndex_t> > AllIntCPNodeNums;
		vector<string> AllIntVolumeCPNames;

		vector<vec3> IntCPPos;

		vec3 MinVolPt, MaxVolPt;


		EntIndex_t NumZones = TecUtilDataSetGetNumZones();

		point* p = NULL;
		triangle* t = NULL;
		int** e = NULL;
		int NumPoints, NumTriangles, NumEdges;

		vector<int> MovedPointNums;

		MeshStatus_e MeshStatus = FAIL_INVALIDCONSTRAINT;
	
		

		/*
		 *	Get min and max XYZ for the system, just in case we need it.
		 */

		if (IsOk){
			EntIndex_t VolZoneNum = ZoneNumByName(string("Full Volume"));
			if (VolZoneNum <= 0){
				VolZoneNum = 1;
				// 			TecUtilDialogErrMsg("Couldn't get volume zone");
				//TecUtilLockFinish(AddOnID);
	// 			return;
			}
			for (int i = 0; i < 3; ++i){
				TecUtilDataValueGetMinMaxByZoneVar(VolZoneNum, XYZVarNums[i], &MinVolPt[i], &MaxVolPt[i]);
			}
		}


		double ClosestCPDist = 1e50;
		if (IsOk){
			for (int i = 1; i <= IJKMax[0] && IsOk; ++i){
				if (i != CPNum){
					vec3 OtherCP;
					for (int j = 0; j < 3; ++j){
						OtherCP[j] = TecUtilDataValueGetByRef(CPXYZPtrs[j], i);
					}
					ClosestCPDist = MIN(ClosestCPDist, Distance(OtherCP, CPPos));
				}
			}
		}

		if (RadiusMode == ABSOLUTERADIUS){
			Radius = UserRadius;
		}
		else if (RadiusMode == MINCPDISTRATIO){
			Radius = UserRadius * ClosestCPDist;
		}


		/*
		 *	Find all volumetric special gradient paths that will intersect
		 *	the sphere.
		 */
		TecUtilDataLoadBegin();
		for (EntIndex_t CurZoneNum = CPZoneNum + 1; CurZoneNum <= OldNumZones && IsOk; ++CurZoneNum){
			LgIndex_t TmpIJK[3];
			TecUtilZoneGetIJK(CurZoneNum, &TmpIJK[0], &TmpIJK[1], &TmpIJK[2]);
			if (TmpIJK[0] > 1 && TmpIJK[1] == 1 && TmpIJK[2] == 1){
				// If 1-d (line) zone.
				// Check that either initial or terminal point is within sphere radius
				FieldData_pa TmpXYZPtrs[3];
				for (int i = 0; i < 3 && IsOk; ++i){
					TmpXYZPtrs[i] = TecUtilDataValueGetReadableRef(CurZoneNum, XYZVarNums[i]);
					IsOk = VALID_REF(TmpXYZPtrs[i]);
				}
				Boolean_t ConnectsToCP = FALSE;
				Boolean_t StartsFrom1;
				LgIndex_t StartEndCPs[2];
				string StartEndTypes[2];
				Boolean_t ZoneOK = TRUE;
				AuxData_pa TempAuxData = TecUtilAuxDataZoneGetRef(CurZoneNum);
				ZoneOK = VALID_REF(TempAuxData);
				if (ZoneOK){
					for (int i = 0; i < 2 && ZoneOK; ++i){
						AuxDataType_e AuxDataTypeJunk;
						Boolean_t BoolJunk;
						char* TempCStr;
						if (ZoneOK){
							if (TecUtilAuxDataGetItemByName(TempAuxData, CCDataGPEndNums[i].c_str(), (ArbParam_t*)&TempCStr, &AuxDataTypeJunk, &BoolJunk)){
								StartEndCPs[i] = stoi(TempCStr);
							}
							TecUtilStringDealloc(&TempCStr);
						}
						if (ZoneOK){
							if (TecUtilAuxDataGetItemByName(TempAuxData, CCDataGPEndTypes[i].c_str(), (ArbParam_t*)&TempCStr, &AuxDataTypeJunk, &BoolJunk)){
								StartEndTypes[i] = TempCStr;
							}
							TecUtilStringDealloc(&TempCStr);
						}
					}
				}
				if (ZoneOK){
					for (int i = 0; i < 2; ++i){
						if (StartEndCPs[i] == CPNum){
							ConnectsToCP = TRUE;
							StartsFrom1 = (i == 0);
							for (int j = 1; j < 3; ++j){
								if (StartEndTypes[(i + 1) % 2] == RankStrs[j]){
									vec3 TmpCP;
									for (int k = 0; k < 3; ++k){
										TmpCP[k] = TecUtilDataValueGetByRef(CPXYZPtrs[k], StartEndCPs[(i + 1) % 2]);
									}
									IntZoneSaddleCPNodeNums.push_back(vector<LgIndex_t>());
									IntZoneSaddleCPNodeNums[IntZoneSaddleCPNodeNums.size() - 1].push_back(IntersectionPoints.size());
									IntZoneSaddleCPNodeNums[IntZoneSaddleCPNodeNums.size() - 1].push_back(CurZoneNum);
									IntZoneSaddleCPNodeNums[IntZoneSaddleCPNodeNums.size() - 1].push_back(StartEndCPs[(i + 1) % 2]);
									IntCPPos.push_back(TmpCP);

									if (TmpIJK[0] < NumSTPoints){
										NumSTPoints = TmpIJK[0];
										if (NumSTPoints % 2 != 0){
										 	NumSTPoints--;
										}
									}
								}
								int CPCount = 0;
								for (int k = 0; k < 4; ++k){
									if (StartEndTypes[(i + 1) % 2] == RankStrs[k]){
										AllIntVolumeCPNames.push_back(StartEndTypes[(i + 1) % 2] + string(" ") + to_string(StartEndCPs[(i + 1) % 2] - CPCount));
										AllIntCPNodeNums.push_back(vector<LgIndex_t>());
										AllIntCPNodeNums[AllIntCPNodeNums.size() - 1].push_back(IntersectionPoints.size());
										AllIntCPNodeNums[AllIntCPNodeNums.size() - 1].push_back(StartEndCPs[(i + 1) % 2]);
									}
									CPCount += NumCPs[k];
								}
							}
						}
					}
				}
// 				if (IsOk){
// 					for (LgIndex_t PtNum = 1; PtNum <= TmpIJK[0]; PtNum += TmpIJK[0] - 1){
// 						Vec3_c TmpPt;
// 						for (int i = 0; i < 3; ++i){
// 							TmpPt[i] = TecUtilDataValueGetByRef(TmpXYZPtrs[i], PtNum);
// 						}
// 						if (Distance(TmpPt, CPPos) < MAX(ClosestCPDist * 0.2, Radius)){
// 							ConnectsToCP = TRUE;
// 							StartsFrom1 = (PtNum == 1);
// 						}
// 					}
// 				}
				if (ConnectsToCP && IsOk){

					/*
					*	Check if the gradient path connects to a near field.
					*	bond or ring. If so, record the zone and CP number.
					*/
// 					for (LgIndex_t PtNum = 1; PtNum <= TmpIJK[0]; PtNum += TmpIJK[0] - 1){
// 						Vec3_c TmpPt;
// 						for (int i = 0; i < 3; ++i){
// 							TmpPt[i] = TecUtilDataValueGetByRef(TmpXYZPtrs[i], PtNum);
// 						}
// 						Boolean_t IsFound = FALSE;
// 						for (int i = 0; i < 4 && !IsFound; ++i){
// 							Vec3_c TmpCP;
// 							for (int j = 0; j < CPIndexList[i].size() && !IsFound; ++j){
// 								if (CPIndexList[i][j] != CPNum){
// 									for (int k = 0; k < 3; ++k){
// 										TmpCP[k] = TecUtilDataValueGetByRef(CPXYZPtrs[k], CPIndexList[i][j]);
// 									}
// 									if (Distance(TmpCP, TmpPt) < MAX(ClosestCPDist * 0.2, Radius)){
// 										if (i == 1 || i == 2){
// 											IntZoneSaddleCPNodeNums.push_back(vector<LgIndex_t>());
// 											IntZoneSaddleCPNodeNums[IntZoneSaddleCPNodeNums.size() - 1].push_back(IntersectionPoints.size());
// 											IntZoneSaddleCPNodeNums[IntZoneSaddleCPNodeNums.size() - 1].push_back(CurZoneNum);
// 											IntZoneSaddleCPNodeNums[IntZoneSaddleCPNodeNums.size() - 1].push_back(CPIndexList[i][j]);
// 											IntCPPos.push_back(TmpCP);
// 											IsFound = TRUE;
// 
// 											if (TmpIJK[0] < NumSTPoints){
// 												NumSTPoints = TmpIJK[0];
// 												if (NumSTPoints % 2 != 0){
// 													NumSTPoints--;
// 												}
// 											}
// 										}
// 										AllIntSaddleCPNodeNums.push_back(vector<LgIndex_t>());
// 										AllIntSaddleCPNodeNums[AllIntSaddleCPNodeNums.size() - 1].push_back(IntersectionPoints.size());
// 										AllIntSaddleCPNodeNums[AllIntSaddleCPNodeNums.size() - 1].push_back(CPIndexList[i][j]);
// 									}
// 								}
// 							}
// 						}
// 					}

					// An intersection exists, so need to find its location.
					vec3 Pt1, Pt2;
					LgIndex_t StartPt = 1, Step = 1;
					if (!StartsFrom1){
						StartPt = TmpIJK[0];
						Step = -1;
					}
					Boolean_t IntersectionFound = FALSE;
					double Dist1,  Dist2;
					for (int i = 0; i < 3; ++i){
						Pt1[i] = TecUtilDataValueGetByRef(TmpXYZPtrs[i], StartPt + Step);
					}
					Pt2 = Pt1;
					Dist2 = Dist1 = Distance(Pt1, CPPos);
					StartPt += Step;
					for (LgIndex_t PtNum = StartPt + Step; PtNum >= 1 && PtNum <= TmpIJK[0] && !IntersectionFound && IsOk; PtNum += Step){
						for (int i = 0; i < 3; ++i){
							Pt1[i] = TecUtilDataValueGetByRef(TmpXYZPtrs[i], PtNum);
						}
						Dist1 = Distance(Pt1, CPPos);
						if ((Dist1 < Radius && Dist2 >= Radius) ||
							(Dist1 >= Radius && Dist2 < Radius)){
							// Now Pt1 and Pt2 are straddling the sphere.
							// Interpolate to find the intersection point.
							double DistanceRatio = (Radius - Dist1) / (Dist2 - Dist1);
							point IntPt;
							for (int i = 0; i < 3; ++i){
								IntPt[i] = (Pt1[i] + DistanceRatio * (Pt2[i] - Pt1[i])) - CPPos[i];
							}

							/*
							* Project intersection point to sphere surface
							*/
							// Compute spherical coordinates
							double   rad = sqrt(pow(IntPt.x, 2) + pow(IntPt.y, 2) + pow(IntPt.z, 2));
							double theta = acos(IntPt.z / rad);
							double   phi = atan2(IntPt.y, IntPt.x);

							// Project point onto a sphere of radius "radius"
							IntPt.x = Radius * sin(theta) * cos(phi);
							IntPt.y = Radius * sin(theta) * sin(phi);
							IntPt.z = Radius * cos(theta);


							IntersectionFound = TRUE;
							/*
							 *	Make sure this intersection isn't too close to another
							 */

							Boolean_t IntersectionGood = TRUE;
							for (int i = 0; i < IntersectionPoints.size() && IntersectionGood; ++i){
								vec3 TmpPt1, TmpPt2;
								for (int j = 0; j < 3; ++j){
									TmpPt1[j] = IntPt[j];
									TmpPt2[j] = IntersectionPoints[i][j];
								}
								IntersectionGood = (Distance(TmpPt1, TmpPt2) >= 0.01);
							}
							
							if (IntersectionGood)
								IntersectionPoints.push_back(IntPt);


						}
						Pt2 = Pt1;
						Dist2 = Dist1;
					}
				}
			}
		}

		TecUtilDataLoadEnd();


		/*
		 *	Get optimized spherical mesh
		 */

		int MaxRefine = 2;
		int RefineNum = 0;
		while(MeshStatus != 0){
			TecUtilDataLoadBegin();/*
			 *	Get memory usage estimate for mesh
			 */
			NumTriangles = 20 * pow(4, Level);
			NumPoints = NumTriangles / 2 + 2;
			NumEdges = NumTriangles * double(3.0 / 2.0);
			/*
			*	point and triangle adjacency lists in meshgen.
			*	Typically have six neighbor point and triangles
			*	per point.
			*/
			int AdjListSize = sizeof(int) * 2 * 6 * NumPoints;
			MemoryRequired = (sizeof(point) * NumPoints 
				+ sizeof(triangle) * NumTriangles + AdjListSize
				+ sizeof(int) * 2 * NumEdges) / 1024;
			TecUtilMemoryChangeNotify(MemoryRequired);



			MeshStatus = meshgen2D_sphere(Radius, Level, IntersectionPoints, MovedPointNums, p, t, e, NumPoints, NumTriangles, NumEdges);



			/*
			 *	The point and triangle lists are kept. only
			 *	the adjacency lists are destroyed.
			 */
			MemoryRequired = (AdjListSize) / 1024;
			TecUtilMemoryChangeNotify(-MemoryRequired);

			TecUtilDataLoadEnd();
			if (MeshStatus == FAIL_NEEDREFINEMENT){
				//break;
				if (RefineNum >= MaxRefine){
					break;
				}
				else{
					//TecUtilDialogMessageBox("Mesh not fine enough for constrained points. Refining mesh.", MessageBoxType_Warning);
					TecUtilMemoryChangeNotify(-(sizeof(point) * NumPoints
						+ sizeof(triangle) * NumTriangles) / 1024);
					delete p, t;
					p = NULL;
					t = NULL;
					e = NULL;
					Level++;
					RefineNum++;
				}
			}
			else if (MeshStatus == FAIL_INVALIDCONSTRAINT){
				TecUtilDrawGraphics(TRUE);
				TecUtilInterfaceSuspend(FALSE);
				TecUtilDialogErrMsg("Constrained point too far from mesh");
				TecUtilLockFinish(AddOnID);
				return;
			}
			else if (MeshStatus == SUCCESS_NOTCONVERGED){
// 				TecUtilInterfaceSuspend(FALSE);
// 				TecUtilDialogMessageBox("Mesh did not converge.", MessageBoxType_Warning);
// 				TecUtilInterfaceSuspend(TRUE);
				break;
			}

		}

		IsOk = (MeshStatus == SUCCESS || MeshStatus == SUCCESS_NOTCONVERGED && NumPoints > 0 && NumTriangles > 0);
		if (!IsOk){
			TecUtilDrawGraphics(TRUE);
			TecUtilInterfaceSuspend(FALSE);
			TecUtilDialogErrMsg("Failed to make mesh.");
			TecUtilLockFinish(AddOnID);
			return;
		}

		if (IsOk){
			/*
			 *	Translate sphere to the CP position
			 */

			TecUtilDataLoadBegin();
			for (int i = 0; i < NumPoints; ++i){
				for (int j = 0; j < 3; ++j){
					p[i][j] += CPPos[j];
				}
			}
			TecUtilDataLoadEnd();

			/*
			*	Add triangle finite element zone for sphere
			*/

			IsOk = TecUtilDataSetAddZone("Sphere", NumPoints, NumTriangles, 0, ZoneType_FETriangle, NULL);
			EntIndex_t SurfZoneNum = TecUtilDataSetGetNumZones();
			FieldData_pa SurfXYZPtrs[3];

			for (int i = 0; i < 3 && IsOk; ++i){
				SurfXYZPtrs[i] = TecUtilDataValueGetWritableRef(SurfZoneNum, XYZVarNums[i]);
				IsOk = VALID_REF(SurfXYZPtrs[i]);
			}
			NodeMap_pa NodeMap;

			if (IsOk){
				NodeMap = TecUtilDataNodeGetWritableRef(SurfZoneNum);
				IsOk = VALID_REF(NodeMap);

				TecUtilDataLoadBegin();

				for (LgIndex_t NodeNum = 0; NodeNum < NumPoints; ++NodeNum){
					for (int i = 0; i < 3; ++i){
						TecUtilDataValueSetByRef(SurfXYZPtrs[i], NodeNum + 1, p[NodeNum][i]);
					}
				}
				for (LgIndex_t TriNum = 0; TriNum < NumTriangles; ++TriNum){
					for (int i = 0; i < 3; ++i){
						TecUtilDataNodeSetByRef(NodeMap, TriNum + 1, i + 1, t[TriNum][i] + 1);
					}
				}

				TecUtilDataLoadEnd();
			}
			if (IsOk){
				Set_pa TmpSet = TecUtilSetAlloc(FALSE);
				TecUtilSetAddMember(TmpSet, SurfZoneNum, FALSE);

				ArgList_pa CurrentArgList = TecUtilArgListAlloc();
				TecUtilArgListAppendString(CurrentArgList, SV_P1, SV_FIELDMAP);
				TecUtilArgListAppendString(CurrentArgList, SV_P2, SV_EFFECTS);
				TecUtilArgListAppendString(CurrentArgList, SV_P3, SV_USETRANSLUCENCY);
				TecUtilArgListAppendSet(CurrentArgList, SV_OBJECTSET, TmpSet);
				TecUtilArgListAppendArbParam(CurrentArgList, SV_IVALUE, FALSE);
				TecUtilStyleSetLowLevelX(CurrentArgList);
				TecUtilArgListDealloc(&CurrentArgList);

				TecUtilStateChanged(StateChange_ZonesAdded, (ArbParam_t)TmpSet);
				TecUtilZoneSetScatter(SV_SHOW, TmpSet, 0.0, FALSE);
				TecUtilZoneSetMesh(SV_SHOW, TmpSet, 0.0, TRUE);
				TecUtilZoneSetActive(TmpSet, AssignOp_PlusEquals);

				TecUtilSetDealloc(&TmpSet);
			}

			/*
			*	Save node and triangle numbers to the fe volume zone's aux data
			*/

			if (IsOk){
				AuxData_pa TempAuxData = TecUtilAuxDataZoneGetRef(SurfZoneNum);
				IsOk = VALID_REF(TempAuxData);
				if (IsOk){
					stringstream ss;
					ss << GBADataPrefix << GBASphereCPName;
					IsOk = TecUtilAuxDataSetItem(TempAuxData, ss.str().c_str(), (ArbParam_t)CPName.c_str(), AuxDataType_String, TRUE);
				}
				if (IsOk){
					stringstream ss;
					ss << GBADataPrefix << GBAZoneType;
					IsOk = TecUtilAuxDataSetItem(TempAuxData, ss.str().c_str(), (ArbParam_t)GBAZoneTypeSphereZone.c_str(), AuxDataType_String, TRUE);
				}
			}
		}

		/*
		*	Get a list of all unique nodes that neighbor nodes that
		*	were moved to GP intersections
		*/
		TecUtilDataLoadBegin();
		vector<vector<int> > ConstrainedNeighborNodesNum(MovedPointNums.size(), vector<int>());
		LgIndex_t NumNeighborNodes = 0;
		for (int i = 0; i < MovedPointNums.size(); ++i){
			ConstrainedNeighborNodesNum[i].reserve(6);
			for (int j = 0; j < NumTriangles; ++j){
				for (int k = 0; k < 3; ++k){
					if (t[j][k] == MovedPointNums[i]){
						for (int ii = 0; ii < 3; ++ii){
							if (t[j][ii] != MovedPointNums[i]){
								Boolean_t IsFound = FALSE;
								for (int jj = 0; jj < ConstrainedNeighborNodesNum[i].size() && !IsFound; ++jj){
									if (t[j][ii] == ConstrainedNeighborNodesNum[i][jj]){
										IsFound = TRUE;
									}
								}
								if (!IsFound)
									ConstrainedNeighborNodesNum[i].push_back(t[j][ii]);
							}
						}
					}
				}
			}
			NumNeighborNodes += ConstrainedNeighborNodesNum[i].size();
			for (int j = 0; j < IntZoneSaddleCPNodeNums.size(); ++j){
				if (i == IntZoneSaddleCPNodeNums[j][0]){
					IntZoneSaddleCPNodeNums[j].push_back(MovedPointNums[i]);
				}
			}
			for (int j = 0; j < AllIntCPNodeNums.size(); ++j){
				if (i == AllIntCPNodeNums[j][0]){
					AllIntCPNodeNums[j].push_back(MovedPointNums[i]);
				}
			}
		}

		TecUtilDataLoadEnd();

		// For storing edges that straddle IBs
		vector<vector<int> > InterIBEdges;

		/*
		 *	Prepare for streamtrace creation
		 */

		EntIndex_t GradXYZVarNums[3];

		if (IsOk){
			vector<string> TmpStrs = {
				"X Density Gradient",
				"Y Density Gradient",
				"Z Density Gradient"
			};
			for (int i = 0; i < 3 && IsOk; ++i){
				GradXYZVarNums[i] = VarNumByName(TmpStrs[i]);
				IsOk = (GradXYZVarNums[i] > 0);
			}
			if (!IsOk){
				TecUtilDrawGraphics(TRUE);
				TecUtilInterfaceSuspend(FALSE);
				TecUtilDialogErrMsg("Couldn't find gradient vector variables.");
				TecUtilStatusSuspend(FALSE);
				TecUtilStatusFinishPercentDone();
				TecUtilDataLoadEnd();
				TecUtilLockFinish(AddOnID);
				return;
			}
		}

		if (IsOk){

			// Enable system volume zone

			Set_pa TmpSet = TecUtilSetAlloc(FALSE);
			TecUtilSetAddMember(TmpSet, 1, FALSE);
			TecUtilZoneSetActive(TmpSet, AssignOp_PlusEquals);
			TecUtilSetDealloc(&TmpSet);

			//	Set streamtrace options 

			vector<string> TmpStrs1 = {
				SV_STREAMTRACELAYERS,
				SV_GLOBALTHREEDVECTOR,
				SV_GLOBALTHREEDVECTOR,
				SV_GLOBALTHREEDVECTOR,
				SV_STREAMATTRIBUTES
			};
			vector<string> TmpStrs2 = {
				SV_SHOW,
				SV_UVAR,
				SV_VVAR,
				SV_WVAR,
				SV_MAXSTEPS
			};
			vector<int> TmpValues = {
				1,
				GradXYZVarNums[0],
				GradXYZVarNums[1],
				GradXYZVarNums[2],
				10000
			};

			for (int i = 0; i < 5 && IsOk; ++i){
				if (IsOk){
					ArgList_pa argList = TecUtilArgListAlloc();
					TecUtilArgListAppendString(argList, SV_P1, TmpStrs1[i].data());
					TecUtilArgListAppendString(argList, SV_P2, TmpStrs2[i].data());
					TecUtilArgListAppendArbParam(argList, SV_IVALUE, TmpValues[i]);
					SetValueReturnCode_e ReturnCode = TecUtilStyleSetLowLevelX(argList);
					TecUtilArgListDealloc(&argList);
					IsOk = (ReturnCode == SetValueReturnCode_Ok || ReturnCode == SetValueReturnCode_DuplicateValue);
				}
			}

			if (IsOk){
				ArgList_pa argList = TecUtilArgListAlloc();
				TecUtilArgListAppendString(argList, SV_P1, SV_STREAMATTRIBUTES);
				TecUtilArgListAppendString(argList, SV_P2, SV_MINCELLFRACTION);
				TecUtilArgListAppendDouble(argList, SV_DVALUE, (double)0.00001);
				SetValueReturnCode_e ReturnCode = TecUtilStyleSetLowLevelX(argList);
				TecUtilArgListDealloc(&argList);
				IsOk = (ReturnCode == SetValueReturnCode_Ok || ReturnCode == SetValueReturnCode_DuplicateValue);
			}
		}


		/*
		 *	Seed streamtraces for all nodes of the sphere mesh.
		 *	If a node is one of the constraints that corresponds to 
		 *	a near-field ring or bond CP, then make sure the 
		 *	streamtrace terminates at the CP.
		 */



		EntIndex_t NumZonesBeforeVolumes = TecUtilDataSetGetNumZones();

		MemoryRequired = (sizeof(bool) * (2 * NumPoints + 2 * NumTriangles) + sizeof(EntIndex_t) * NumTriangles) / 1024;
		TecUtilMemoryChangeNotify(MemoryRequired);

		vector<bool> NodeHasSaddleCP(NumPoints, false),
			NodeIsConstrained(NumPoints, false),
			TriangleHasSaddleCP(NumTriangles, false),
			TriangleHasConstrainedNode(NumTriangles, false);
		for (int i = 0; i < IntZoneSaddleCPNodeNums.size(); ++i){
			NodeHasSaddleCP[IntZoneSaddleCPNodeNums[i][3]] = true;
		}
		for (int i = 0; i < MovedPointNums.size(); ++i){
			NodeIsConstrained[MovedPointNums[i]] = true;
		}

		vector<EntIndex_t> TriangleFEZoneNums(NumTriangles);

		/*
		 *	Start triangle-based seeding
		 */

		// Just do XYZ values for these volume zones.
		EntIndex_t ActualNumVars = NumVars;
		NumVars = 3;


		if (IsOk){

			LgIndex_t NumNodes;
			LgIndex_t NumElems;

			TecUtilDrawGraphics(TRUE);
			TecUtilPleaseWait(NULL, FALSE);
			TecUtilStatusStartPercentDone(ProgressStr.str().c_str(), FALSE, TRUE);
			TecUtilStatusSuspend(TRUE);
			TecUtilDrawGraphics(FALSE);


			for (int TriNum = 0; TriNum < NumTriangles && IsOk; ++TriNum){


				TecUtilDataLoadBegin();

				int STZoneNums[5] = { -1, -1, -1, -1, -1 };
				bool ConcatenateST[4] = { false, false, false, false };
				int DoubleCornerNum;
				int STCount = 0;
				vec3 DoubleCornerNodePos;

				int NewSTZoneNum = TecUtilDataSetGetNumZones();

				for (int CornerNum = 0; CornerNum < 3 && IsOk; ++CornerNum){
					if (!NodeHasSaddleCP[t[TriNum][CornerNum]]){
						if (NodeIsConstrained[t[TriNum][CornerNum]]){
							TriangleHasConstrainedNode[TriNum] = true;
						}

						IsOk = TecUtilStreamtraceAdd(1, Streamtrace_VolumeLine,
							StreamDir_Both, p[t[TriNum][CornerNum]].x,
							p[t[TriNum][CornerNum]].y,
							p[t[TriNum][CornerNum]].z,
							0, 0, 0);

						STZoneNums[CornerNum] = ++NewSTZoneNum;

						STCount++;
					}
				}

				if (STCount < 3){
					if (IsOk){
						IsOk = TecUtilCreateStreamZones(FALSE);
					}
					if (IsOk){
						IsOk = TecUtilStreamtraceDeleteAll();
					}

					if (IsOk){
						vector<vector<FieldData_pa>> STGetPtrs(2, vector<FieldData_pa>(3, NULL));
						vector<vector<FieldValueGetFunction_pf> > STGetFuncts(2, vector<FieldValueGetFunction_pf>(3, NULL));
						LgIndex_t TmpIJK[2][3];
						int CornerNum;

						for (int k = 0; k < 3; ++k){
							if (STZoneNums[k] <= 0){
								for (int i = 0; i < 2 && IsOk; ++i){
									int OtherCornerNum = (k + i + 1) % 3;
									if (STZoneNums[OtherCornerNum] > 0){
										TecUtilZoneGetIJK(STZoneNums[OtherCornerNum], &TmpIJK[i][0], &TmpIJK[i][1], &TmpIJK[i][2]);
										for (int j = 0; j < 3 && IsOk; ++j){
											STGetPtrs[i][j] = TecUtilDataValueGetReadableRef(STZoneNums[OtherCornerNum], XYZVarNums[j]);
											IsOk = (VALID_REF(STGetPtrs[i][j]));
											if (IsOk){
												STGetFuncts[i][j] = TecUtilDataValueRefGetGetFunc(STGetPtrs[i][j]);
												IsOk = (VALID_REF(STGetFuncts[i][j]));
											}
										}
									}
									else{
										TecUtilDrawGraphics(TRUE);
										TecUtilInterfaceSuspend(FALSE);
										TecUtilDialogErrMsg("Incorrect zone number");
										TecUtilStatusSuspend(FALSE);
										TecUtilStatusFinishPercentDone();
										TecUtilDataLoadEnd();
										TecUtilLockFinish(AddOnID);
										return;
									}
								}
								CornerNum = k;
								break;
							}
						}




						omp_set_num_threads(2);

						vec3 SeedPts[2];
						vec3 NodePos;
						for (int i = 0; i < 3; ++i)
							NodePos[i] = p[t[TriNum][CornerNum]][i];

						/*
						*	Get node's associated volume CP, volume GP.
						*/
						int VolCPNum = -1,
							VolGPNum = -1;
						vec3 VolCPPos;
						for (int i = 0; i < IntZoneSaddleCPNodeNums.size(); ++i){
							if (IntZoneSaddleCPNodeNums[i][3] == t[TriNum][CornerNum]){
								VolGPNum = IntZoneSaddleCPNodeNums[i][1];
								VolCPNum = IntZoneSaddleCPNodeNums[i][2];
								VolCPPos = IntCPPos[i];
							}
						}

						STZoneNums[4] = VolGPNum;

						TriangleHasSaddleCP[TriNum] = true;
						TriangleHasConstrainedNode[TriNum] = true;

						DoubleCornerNum = CornerNum;
						DoubleCornerNodePos = NodePos;
	#pragma omp parallel for
						for (int CornerI = 0; CornerI < 2; ++CornerI){
						

							/*
							*	For each neighboring node, need to get an orthogonal set
							*	of unit vectors to describe the positions of the neighbor
							*	node's streamtraces's closest point to the volume CP
							*	relative to the volume CP
							*/
							double TestDist = 100.0 * Distance(NodePos, VolCPPos);
							/*
							*	Find the point on the neighbor node's streamtrace that
							*	is closest to the volume CP
							*/

							vec3 ClosestPoint;
							double MinSqrDist = TestDist;
							if (IsOk){
								vec3 TmpPt;
								for (int STPtNum = 0; STPtNum < TmpIJK[CornerI][0]; ++STPtNum){
									for (int i = 0; i < 3; ++i){
										TmpPt[i] = STGetFuncts[CornerI][i](STGetPtrs[CornerI][i], STPtNum);
									}
									if (sum(TmpPt == zeros<vec>(3)) == 3){
										int a = 1;
									}
									double TmpSqrDist = MIN(MinSqrDist, DistSqr(TmpPt, VolCPPos));
									if (TmpSqrDist < MinSqrDist){
										ClosestPoint = TmpPt;
										MinSqrDist = TmpSqrDist;
									}
								}
							}
							/*
							*	v1 is the direction from the main node to the volume CP.
							*	v2 is the direction from main node to neighbor node, and is
							*	orthogonal to v1.
							*	v3 is orthogonal to v1 and v2.
							*/

							vec3 v1 = VolCPPos - NodePos,
								v2 = ClosestPoint - VolCPPos,
								v3 = cross(v1, v2);
							v1 = normalise( v1);
							v3 = normalise( v3);
							v2 = normalise( cross(v3, v1));

							vec3 SeedPt = v1 * (Radius * 0.1);
							double SeedTheta = PI / 2.0;

							vec4 TmpVec4 = join_cols(SeedPt, ones<vec>(1));
							TmpVec4 = RotationMatrix(SeedTheta, v3) * TmpVec4;
							SeedPts[CornerI] = vec3(TmpVec4) + VolCPPos;
						}

						StreamDir_e StreamDir;
						if (CPType < 0){
							StreamDir = StreamDir_Reverse;
						}
						else{
							StreamDir = StreamDir_Forward;
						}

						for (int CornerI = 0; CornerI < 2 && IsOk; ++CornerI){
							IsOk = TecUtilStreamtraceAdd(1, Streamtrace_VolumeLine, StreamDir, SeedPts[CornerI][0], SeedPts[CornerI][1], SeedPts[CornerI][2], 0, 0, 0);

							if (CornerI == 0){
								STZoneNums[CornerNum] = ++NewSTZoneNum;
								ConcatenateST[CornerNum] = true;
							}
							else{
								STZoneNums[3] = ++NewSTZoneNum;
								ConcatenateST[3] = true;
							}
						}

						if (!IsOk){
							TecUtilDrawGraphics(TRUE);
							TecUtilInterfaceSuspend(FALSE);
							TecUtilDialogErrMsg("Failed during saddle CP streamtraces.");
							TecUtilStatusSuspend(FALSE);
							TecUtilStatusFinishPercentDone();
							TecUtilDataLoadEnd();
							TecUtilLockFinish(AddOnID);
							return;
						}
					}
				}

				if (IsOk){
					IsOk = TecUtilCreateStreamZones(FALSE);
				}
				if (IsOk){
					IsOk = TecUtilStreamtraceDeleteAll();
				}

				/*
				 *	Streamtraces are made for this triangle.
				 *	Now need to create the FE_Quad zone that will
				 *	be this triangle's volume.
				 *	For now, assume both ends of triagonal tube are
				 *	capped with triangles.
				 */

				vector<vector<vector<double> > > ResampledStreamtraces;
				if (TriangleHasSaddleCP[TriNum]){
					NumNodes = 4 * NumSTPoints;
					NumElems = 2 + 4 * (NumSTPoints - 1);

					MemoryRequired = sizeof(double) * 4 * NumSTPoints * NumVars / 1024;
					TecUtilMemoryChangeNotify(MemoryRequired);
					ResampledStreamtraces.resize(4, vector<vector<double> >(NumSTPoints, vector<double>(NumVars, 0)));
				}
				else
				{
					NumNodes = 3 * NumSTPoints;
					NumElems = 2 + 3 * (NumSTPoints - 1);

					MemoryRequired = sizeof(double) * 3 * NumSTPoints * NumVars / 1024;
					TecUtilMemoryChangeNotify(MemoryRequired);
					ResampledStreamtraces.resize(3, vector<vector<double> >(NumSTPoints, vector<double>(NumVars, 0)));
				}
		

				if (IsOk){
					FieldData_pa CutoffVarFDPtrs[5];
					FieldValueGetFunction_pf CutoffVarGetFuncts[5];
				
					vector<vector<FieldData_pa>> STGetPtrs(5, vector<FieldData_pa>(NumVars, NULL));
					vector<vector<FieldValueGetFunction_pf> > STGetFuncts(5, vector<FieldValueGetFunction_pf>(NumVars, NULL));

					LgIndex_t STIJKMax[5][3];

					for (int i = 0; i < 5 && IsOk; ++i){
						if (STZoneNums[i] > 0){
							TecUtilZoneGetIJK(STZoneNums[i], &STIJKMax[i][0], &STIJKMax[i][1], &STIJKMax[i][2]);
							IsOk = (STIJKMax[i][0] > 1 && STIJKMax[i][1] == 1 && STIJKMax[i][2] == 1);
							for (int j = 0; j < 3 && IsOk; ++j){
								STGetPtrs[i][j] = TecUtilDataValueGetReadableRef(STZoneNums[i], XYZVarNums[j]);
								IsOk = (VALID_REF(STGetPtrs[i][j]));
								if (IsOk){
									STGetFuncts[i][j] = TecUtilDataValueRefGetGetFunc(STGetPtrs[i][j]);
									IsOk = (VALID_REF(STGetFuncts[i][j]));
								}
							}

							CutoffVarFDPtrs[i] = TecUtilDataValueGetReadableRef(STZoneNums[i], CutoffVarNum);
							IsOk = VALID_REF(CutoffVarFDPtrs[i]);
							if (IsOk){
								CutoffVarGetFuncts[i] = TecUtilDataValueRefGetGetFunc(CutoffVarFDPtrs[i]);
								IsOk = VALID_REF(CutoffVarGetFuncts[i]);
							}
						}
					}
					EntIndex_t VarNum = 3;
					for (int i = 1; i <= NumVars && IsOk; ++i){
						Boolean_t IsXYZVarNum = FALSE;
						for (int j = 0; j < 3 && !IsXYZVarNum; ++j){
							IsXYZVarNum = (XYZVarNums[j] == i);
						}
						if (!IsXYZVarNum){
							for (int j = 0; j < 5 && IsOk; ++j){
								if (STZoneNums[j] > 0){
									STGetPtrs[j][VarNum] = TecUtilDataValueGetReadableRef(STZoneNums[j], i);
									IsOk = (VALID_REF(STGetPtrs[j][VarNum]));
									if (IsOk){
										STGetFuncts[j][VarNum] = TecUtilDataValueRefGetGetFunc(STGetPtrs[j][VarNum]);
										IsOk = (VALID_REF(STGetFuncts[j][VarNum]));
									}
								}
							}
							VarNum++;
						}
					}
					if (!IsOk){
						TecUtilDrawGraphics(TRUE);
						TecUtilInterfaceSuspend(FALSE);
						TecUtilDialogErrMsg("Failed to get streamtrace zone pointers for resampling");
						TecUtilStatusSuspend(FALSE);
						TecUtilStatusFinishPercentDone();
						TecUtilDataLoadEnd();
						TecUtilLockFinish(AddOnID);
						return;
					}

					if (TriangleHasSaddleCP[TriNum]){
						omp_set_num_threads(4);
					}
					else{
						omp_set_num_threads(3);
					}

	#pragma omp parallel for
					for (int CornerNum = 0; CornerNum < 4; ++CornerNum){
						if (ConcatenateST[CornerNum]){
							/*
							*	Get length of streamtrace zones to be resampled and
							*	distance between each point of resampled zones.
							*	Here we also find the REAL max I value for the streamtraces,
							*	since it usually stops moving well before it reaches its
							*	max number of points.
							*/
							double DelLength;
							double STLength[2] = { 0.0 };

							int STNum = CornerNum;

							EntIndex_t GetSTZoneNums[2] = { STZoneNums[STNum], STZoneNums[4] };
							int PtrNums[2] = { STNum, 4 };
							Boolean_t StartsAtNode;

							/*
							*	Need to make sure that we start at the far end of one of the
							*	streamtraces so that when we concatenate the two, its one
							*	continuous line.
							*/

							Boolean_t StartFrom1[2] = { TRUE };
							double MinGap = 1e50;
							int HalfNumSTPoints[2];
							if (IsOk){
								vec3 EndPoints[4];
								for (int i = 0; i < 2; ++i){
									for (int j = 0; j < 3; ++j){
										EndPoints[2 * i][j] = STGetFuncts[PtrNums[i]][j](STGetPtrs[PtrNums[i]][j], 0);
										EndPoints[1 + 2 * i][j] = STGetFuncts[PtrNums[i]][j](STGetPtrs[PtrNums[i]][j], STIJKMax[PtrNums[i]][0] - 1);
									}
								}
								int GapNum = 1, MinGapNum = 1;
								for (int i = 0; i < 2; ++i){
									for (int j = 2; j < 4; ++j){
										double TmpGap = DistSqr(EndPoints[i], EndPoints[j]);
										if (TmpGap < MinGap){
											MinGap = TmpGap;
											MinGapNum = GapNum;
										}
										GapNum++;
									}
								}
								switch (MinGapNum){
									case 1:
										/*
										*	point 1 of the two streamtraces are next to eachother, so go in
										*	order last1 first1 first2 last2
										*/
										StartFrom1[0] = FALSE;
										StartFrom1[1] = TRUE;
										break;
									case 2:
										/*
										*	point 1 of the first streamtrace is next to last of second one, so go in
										*	order last1 first1 last2 first 2
										*/
										StartFrom1[0] = FALSE;
										StartFrom1[1] = FALSE;
										break;
									case 3:
										/*
										*	last point of first streamtrace next to first point of second, so go in
										*	order first1 last1 first2 last2
										*/
										StartFrom1[0] = TRUE;
										StartFrom1[1] = TRUE;
										break;
									default:
										/*
										*	last point of first streamtrace next to last point of second, so go in
										*	order first1 last1 last2 first2
										*/
										StartFrom1[0] = TRUE;
										StartFrom1[1] = FALSE;
										break;
								}
							}



							/*
							*	Both paths get their own LowI and HighI but we'll only use two
							*	of the four values.
							*/

							Boolean_t CutoffPtFound = FALSE;
							LgIndex_t LowI[2] = { 0, 0 };
							LgIndex_t HighI[2] = { STIJKMax[PtrNums[0]][0] - 1, STIJKMax[PtrNums[1]][0] - 1 };
							vec3 StartPos, EndPos;

							for (int PassNum = 0; PassNum < 2; ++PassNum){
								/*
								*	Figure out if the beginning or end of the streamtrace is closer to the seed point
								*/
								if (IsOk && PassNum == 1){
									vec3 BegPt, EndPt, EndPt2;
									for (int i = 0; i < 3; ++i){
										BegPt[i] = STGetFuncts[PtrNums[PassNum]][i](STGetPtrs[PtrNums[PassNum]][i], 0);
										EndPt[i] = STGetFuncts[PtrNums[PassNum]][i](STGetPtrs[PtrNums[PassNum]][i], STIJKMax[PtrNums[PassNum]][0] - 1);
									}
									EndPt2 = EndPt;
									int PtNum = STIJKMax[PtrNums[PassNum]][0] - 1;
									while (sum(EndPt == EndPt2) == 3 && PtNum > 1){
										PtNum--;
										EndPt2 = EndPt;
										for (int i = 0; i < 3; ++i){
											EndPt[i] = STGetFuncts[PtrNums[PassNum]][i](STGetPtrs[PtrNums[PassNum]][i], PtNum);
										}
									}
									STIJKMax[PtrNums[PassNum]][0] = PtNum + 1;
									EndPt = EndPt2;
									StartsAtNode = (Distance(CPPos, BegPt) < Distance(CPPos, EndPt));
								}

								/*
								*	Now need to find the intersection of the existing gradient path
								*	with the sphere, and use that point as the start point.
								*	Also need to find the intersection of the streamtrace
								*	with the user-defined isosurface and use that at the end point
								*/

								/*
									*	First, the intersection of the existing gradient path
									*	with the sphere.
									*/

								if (IsOk){
									vec3 PtI, PtIm1;
									int StartNum = 0;
									int Step = 1;
									if (PassNum == 1){
										if (!StartsAtNode){
											StartNum = STIJKMax[PtrNums[PassNum]][0] - 1;
											Step = -1;
										}
										for (int i = 0; i < 3; ++i)
											PtIm1[i] = STGetFuncts[PtrNums[PassNum]][i](STGetPtrs[PtrNums[PassNum]][i], StartNum);
										double Dist1, Dist2 = Distance(PtIm1, CPPos);

										StartNum += Step;

										for (LgIndex_t i = StartNum; i < STIJKMax[PtrNums[PassNum]][0] && i >= 0; i += Step){
											for (int j = 0; j < 3; ++j)
												PtI[j] = STGetFuncts[PtrNums[PassNum]][j](STGetPtrs[PtrNums[PassNum]][j], i);

											Dist1 = Distance(PtI, CPPos);
											if ((Dist1 < Radius && Dist2 >= Radius) ||
												(Dist1 >= Radius && Dist2 < Radius)){
												// Now Pt1 and Pt2 are straddling the sphere.
												// Interpolate to find the intersection point.
												double DistanceRatio = (Radius - Dist1) / (Dist2 - Dist1);

												StartPos = (PtI + ((PtIm1 - PtI) * DistanceRatio)) - CPPos;
												//TODO
												/*
												*	Here, if other variable values become desired, need to interpolate
												*	to find their values and store them so they can be used below
												*/

												/*
												* Project intersection point to sphere surface
												*/
												// Compute spherical coordinates
												double   rad = sqrt(pow(StartPos[0], 2) + pow(StartPos[1], 2) + pow(StartPos[2], 2));
												double theta = acos(StartPos[2] / rad);
												double   phi = atan2(StartPos[1], StartPos[0]);

												// Project point onto a sphere of radius "radius"
												StartPos[0] = Radius * sin(theta) * cos(phi);
												StartPos[1] = Radius * sin(theta) * sin(phi);
												StartPos[2] = Radius * cos(theta);

												StartPos = StartPos + CPPos;

												if (StartsAtNode){
													LowI[PassNum] = i;
												}
												else{
													/*
													*	We treat STIJKMax[0] as the upper bound in the below code,
													*	so adding one here lets us get to the value we want (i)
													*	by maintaining < STIJKMax[0]
													*/
													HighI[PassNum] = i + 1;
												}

												break;
											}
											Dist2 = Dist1;
											PtIm1 = PtI;
										}
									}

									/*
									*	Now find the intersection of the streamtrace
									*	with the isosurface (if one is used)
									*/

									if (IsOk && UseCutoff && PassNum == 0){
										StartNum = STIJKMax[PtrNums[PassNum]][0] - 1;
										Step = -1;
										if (!StartFrom1[PassNum]){
											StartNum = 0;
											Step = 1;
										}

										double RhoI, RhoIm1;
										for (int i = 0; i < 3; ++i)
											PtIm1[i] = STGetFuncts[PtrNums[PassNum]][i](STGetPtrs[PtrNums[PassNum]][i], StartNum);
										RhoIm1 = CutoffVarGetFuncts[PtrNums[PassNum]](CutoffVarFDPtrs[PtrNums[PassNum]], StartNum);

										StartNum += Step;

										for (LgIndex_t i = StartNum; i < STIJKMax[PtrNums[PassNum]][0] && i >= 0; i += Step){
											for (int j = 0; j < 3; ++j)
												PtI[j] = STGetFuncts[PtrNums[PassNum]][j](STGetPtrs[PtrNums[PassNum]][j], i);

											RhoI = CutoffVarGetFuncts[PtrNums[PassNum]](CutoffVarFDPtrs[PtrNums[PassNum]], i);

											if ((RhoI < CutoffVal && RhoIm1 >= CutoffVal) ||
												(RhoI >= CutoffVal && RhoIm1 < CutoffVal)){
												// Now Pt1 and Pt2 are straddling the sphere.
												// Interpolate to find the intersection point.
												double RhoRatio = (CutoffVal - RhoI) / (RhoIm1 - RhoI);

												EndPos = (PtI + ((PtIm1 - PtI) * RhoRatio));

												//TODO
												/*
												*	Here, if other variable values become desired, need to interpolate
												*	to find their values and store them so they can be used below
												*/

												if (!StartFrom1[PassNum]){
													HighI[PassNum] = i;
												}
												else{
													LowI[PassNum] = i;
												}

												CutoffPtFound = TRUE;

												break;
											}
											RhoIm1 = RhoI;
											PtIm1 = PtI;
										}
									}
								}
							}

							if (IsOk){
								for (int PassNum = 0; PassNum < 2; ++PassNum){
									vec3 PtI = MinVolPt - 10, PtIm1 = MinVolPt - 10;

									int StartNum = LowI[PassNum];
									int Step = 1;
									if (!StartFrom1[PassNum]){
										StartNum = HighI[PassNum];
										Step = -1;
									}

									if (PassNum == 1 && !StartFrom1[PassNum]){
										PtIm1 = StartPos;
									}
									else if (PassNum == 0 && CutoffPtFound){
										PtIm1 = EndPos;
									}
									else{
										for (int i = 0; i < 3; ++i)
											PtIm1[i] = STGetFuncts[PtrNums[PassNum]][i](STGetPtrs[PtrNums[PassNum]][i], StartNum);
									}

									for (LgIndex_t i = StartNum; i < HighI[PassNum] + 1 && i >= LowI[PassNum]; i += Step){
										for (int j = 0; j < 3; ++j)
											PtI[j] = STGetFuncts[PtrNums[PassNum]][j](STGetPtrs[PtrNums[PassNum]][j], i);
										STLength[PassNum] += Distance(PtI, PtIm1);
										PtIm1 = PtI;
									}

									if (PassNum == 1 && StartFrom1[PassNum]){
										STLength[PassNum] += Distance(PtI, StartPos);
									}
								}
							}

							double TotalLength = STLength[0] + STLength[1] + sqrt(MinGap);
							HalfNumSTPoints[0] = int(STLength[0] / TotalLength * (double)NumSTPoints);
							HalfNumSTPoints[1] = NumSTPoints - HalfNumSTPoints[0];
							DelLength = TotalLength / double(NumSTPoints + 1);


							/*
							*	Now step down streamtrace zone again, adding a point to the
							*	resampled zone each time a distance of DelLength is reached,
							*	interpolating values of all variables each time. First and
							*	last points will always be included.
							*/
							if (IsOk){
								vec3 PtI = MinVolPt - 10, PtIm1 = MinVolPt - 10;
								vector<double> ValsI(NumVars, MinVolPt[0] - 10);
								vector<double> ValsIm1(NumVars, MinVolPt[0] - 10);
								double ArcLength = 0.0,
									ArcLengthI = 0.0,
									ArcLengthIm1 = 0.0;
								int NewI = 0;

								for (int PassNum = 1; PassNum >= 0; --PassNum){

									int StartNum = LowI[PassNum];
									int Step = 1;
									if (PassNum == 0 && StartFrom1[PassNum] && CPType < 0){
										StartNum = HighI[PassNum];
										Step = -1;
									}
									else if (PassNum == 1 && !StartsAtNode){
										StartNum = HighI[PassNum];
										Step = -1;
									}

									int OldI = StartNum;

									//	Set first point and get values needed for first step in loop

									if (PassNum == 1){
										PtI = StartPos;
										for (int i = 0; i < 3; ++i){
											ResampledStreamtraces[STNum][NewI][i] = StartPos[i];
											ValsI[i] = PtI[i];
										}
									}
									else{
										for (int i = 0; i < NumVars; ++i){
											ValsI[i] = STGetFuncts[PtrNums[PassNum]][i](STGetPtrs[PtrNums[PassNum]][i], OldI);
											ResampledStreamtraces[STNum][NewI][i] = ValsI[i];
											if (i < 3){
												PtI[i] = ValsI[i];
											}
										}
									}

									NewI++;
									OldI += Step;
									int NewIEnd = HalfNumSTPoints[PassNum];
									if (PassNum == 0){
										NewIEnd += HalfNumSTPoints[1];
										ArcLengthI += Distance(PtI, PtIm1);
									}

#if (0)
									bool BreakBool = true;
									vector<int> ZoneList = { 389, 390, 391, 396, 398, 399 };
									for (vector<int>::iterator it = ZoneList.begin(); it != ZoneList.end() && BreakBool; it++)
										if (TriNum == *it)
											BreakBool = false;
#endif

									PtIm1 = PtI;
									ValsIm1 = ValsI;

									for (NewI; NewI < NewIEnd - int(PassNum == 0 && CutoffPtFound); ++NewI){
										ArcLength += DelLength;

										/*
										*	Move along old path until a segment of DelLength has
										*	been traversed, then interpolate to get new values to
										*	add to new streamtrace zone
										*/
										
										while (OldI < HighI[PassNum] + 1 && OldI >= LowI[PassNum] && ArcLengthI < ArcLength){
											ArcLengthIm1 = ArcLengthI;
											PtIm1 = PtI;
											ValsIm1 = ValsI;

											for (int i = 0; i < NumVars; ++i){
												ValsI[i] = STGetFuncts[PtrNums[PassNum]][i](STGetPtrs[PtrNums[PassNum]][i], OldI);
												if (i < 3){
													PtI[i] = ValsI[i];
												}
											}

											ArcLengthI += Distance(PtI, PtIm1);

											OldI += Step;
										}

										/*
										*	DelLength has been reached, so time to add a new point
										*	to the new streamtrace zone
										*/
										double Ratio = (ArcLength - ArcLengthIm1) / (ArcLengthI - ArcLengthIm1);
										for (int i = 0; i < NumVars; ++i){
											ResampledStreamtraces[STNum][NewI][i] = ValsIm1[i] + Ratio * (ValsI[i] - ValsIm1[i]);
										}
										if (OldI >= HighI[PassNum] + 1 || OldI < LowI[PassNum]){
											while (NewI < NewIEnd - int(PassNum == 0 && CutoffPtFound)){
												NewI++;
												if (NewI < NumSTPoints){
													for (int i = 0; i < NumVars; ++i){
														ResampledStreamtraces[STNum][NewI][i] = ValsIm1[i] + Ratio * (ValsI[i] - ValsIm1[i]);
													}
												}
											}
										}
									}
								}

								if (CutoffPtFound){
									for (int i = 0; i < 3; ++i){
										ResampledStreamtraces[STNum][NumSTPoints - 1][i] = EndPos[i];
									}
								}

							}
						}
						else if (STZoneNums[CornerNum] > 0){
							/*
							*	Get length of streamtrace zone to be resampled and
							*	distance between each point of resampled zone.
							*	Here we also find the REAL max I value for the streamtrace,
							*	since it usually stops moving well before it reaches its
							*	max number of points.
							*/

							double STLength = 0.0, DelLength;
							/*
							 *	Figure out if the beginning or end of the streamtrace is closer to the seed point
							 */
							Boolean_t StartsAtNode;
							if (IsOk){
								vec3 BegPt, EndPt, EndPt2;
								for (int i = 0; i < 3; ++i){
									BegPt[i] = STGetFuncts[CornerNum][i](STGetPtrs[CornerNum][i], 0);
									EndPt[i] = STGetFuncts[CornerNum][i](STGetPtrs[CornerNum][i], STIJKMax[CornerNum][0]-1);
								}
								EndPt2 = EndPt;
								int PtNum = STIJKMax[CornerNum][0] - 1;
								while (sum(EndPt == EndPt2) == 3 && PtNum > 1){
									PtNum--;
									EndPt2 = EndPt;
									for (int i = 0; i < 3; ++i){
										EndPt[i] = STGetFuncts[CornerNum][i](STGetPtrs[CornerNum][i], PtNum);
									}
								}
								STIJKMax[CornerNum][0] = PtNum + 1;
								EndPt = EndPt2;
								StartsAtNode = (Distance(CPPos, BegPt) < Distance(CPPos, EndPt));
							}

							/*
							 *	Now need to find the intersection of the streamtrace
							 *	with the sphere, and use that point as the start point.
							 *	Also need to find the intersection of the streamtrace
							 *	with the user-defined isosurface
							 */

							Boolean_t CutoffPtFound = FALSE;
							LgIndex_t LowI = 0, HighI = STIJKMax[CornerNum][0] - 1;
							vec3 StartPos, EndPos;

							if (IsOk){
								vec3 PtI, PtIm1;
								int StartNum = 0;
								int Step = 1;
								if (!StartsAtNode){
									StartNum = STIJKMax[CornerNum][0] - 1;
									Step = -1;
								}
								for (int i = 0; i < 3; ++i)
									PtIm1[i] = STGetFuncts[CornerNum][i](STGetPtrs[CornerNum][i], StartNum);
								double Dist1, Dist2 = Distance(PtIm1, CPPos);

								StartNum += Step;

								for (LgIndex_t i = StartNum; i < STIJKMax[CornerNum][0] && i >= 0; i += Step){
									for (int j = 0; j < 3; ++j)
										PtI[j] = STGetFuncts[CornerNum][j](STGetPtrs[CornerNum][j], i);

									Dist1 = Distance(PtI, CPPos);
									if ((Dist1 < Radius && Dist2 >= Radius) ||
										(Dist1 >= Radius && Dist2 < Radius)){
										// Now Pt1 and Pt2 are straddling the sphere.
										// Interpolate to find the intersection point.
										double DistanceRatio = (Radius - Dist1) / (Dist2 - Dist1);

										StartPos = (PtI + ((PtIm1 - PtI) * DistanceRatio)) - CPPos;
										//TODO
										/*
										 *	Here, if other variable values become desired, need to interpolate
										 *	to find their values and store them so they can be used below
										 */

										/*
										* Project intersection point to sphere surface
										*/
										// Compute spherical coordinates
										double   rad = sqrt(pow(StartPos[0], 2) + pow(StartPos[1], 2) + pow(StartPos[2], 2));
										double theta = acos(StartPos[2] / rad);
										double   phi = atan2(StartPos[1], StartPos[0]);

										// Project point onto a sphere of radius "radius"
										StartPos[0] = Radius * sin(theta) * cos(phi);
										StartPos[1] = Radius * sin(theta) * sin(phi);
										StartPos[2] = Radius * cos(theta);

										StartPos = StartPos + CPPos;

										if (StartsAtNode){
											LowI = i;
										}
										else{
											/*
											 *	We treat STIJKMax[0] as the upper bound in the below code,
											 *	so adding one here lets us get to the value we want (i)
											 *	by maintaining < STIJKMax[0]
											 */
											HighI = i + 1;
											STIJKMax[CornerNum][0] = i + 1;
										}

										break;
									}
									Dist2 = Dist1;
									PtIm1 = PtI;
								}

								/*
								 *	Now find the intersection with the isosurface (if one is used)
								 */
								if (UseCutoff){
									StartNum = STIJKMax[CornerNum][0] - 1;
									Step = -1;
									if (!StartsAtNode){
										StartNum = 0;
										Step = 1;
									}

									double RhoI, RhoIm1;
									for (int i = 0; i < 3; ++i)
										PtIm1[i] = STGetFuncts[CornerNum][i](STGetPtrs[CornerNum][i], StartNum);
									RhoIm1 = CutoffVarGetFuncts[CornerNum](CutoffVarFDPtrs[CornerNum], StartNum);

									StartNum += Step;

									for (LgIndex_t i = StartNum; i < STIJKMax[CornerNum][0] && i >= 0; i += Step){
										for (int j = 0; j < 3; ++j)
											PtI[j] = STGetFuncts[CornerNum][j](STGetPtrs[CornerNum][j], i);

										RhoI = CutoffVarGetFuncts[CornerNum](CutoffVarFDPtrs[CornerNum], i);

										if ((RhoI < CutoffVal && RhoIm1 >= CutoffVal) ||
											(RhoI >= CutoffVal && RhoIm1 < CutoffVal)){
											// Now Pt1 and Pt2 are straddling the sphere.
											// Interpolate to find the intersection point.
											double RhoRatio = (CutoffVal - RhoI) / (RhoIm1 - RhoI);

											EndPos = (PtI + ((PtIm1 - PtI) * RhoRatio));

											//TODO
											/*
											*	Here, if other variable values become desired, need to interpolate
											*	to find their values and store them so they can be used below
											*/

											if (StartsAtNode){
												HighI = i;
											}
											else{
												LowI = i;
											}

											CutoffPtFound = TRUE;

											break;
										}
										RhoIm1 = RhoI;
										PtIm1 = PtI;
									}
								}
							}

							if (IsOk){
								vec3 PtI = MinVolPt - 10, PtIm1 = MinVolPt - 10;

								int StartNum = LowI;
								int Step = 1;
								if (!StartsAtNode){
									StartNum = HighI;
									Step = -1;
								}

								PtIm1 = StartPos;
								for (LgIndex_t i = StartNum; i < HighI + 1 && i >= LowI; i += Step){
									for (int j = 0; j < 3; ++j)
										PtI[j] = STGetFuncts[CornerNum][j](STGetPtrs[CornerNum][j], i);
									STLength += Distance(PtI, PtIm1);
									PtIm1 = PtI;
								}
								DelLength = STLength / double(NumSTPoints + 1);
							}

							/*
							*	Now step down streamtrace zone again, adding a point to the
							*	resampled zone each time a distance of DelLength is reached,
							*	interpolating values of all variables each time. First and
							*	last points will always be included.
							*/
							if (IsOk){
								vec3 PtI = MinVolPt - 10, PtIm1 = MinVolPt - 10;

								int StartNum = LowI;
								int Step = 1;
								if (!StartsAtNode){
									StartNum = HighI;
									Step = -1;
								}

								vector<double> ValsI(NumVars, MinVolPt[0] - 10);
								vector<double> ValsIm1(NumVars, MinVolPt[0] - 10);
								double ArcLength = 0.0,
									ArcLengthI = 0.0,
									ArcLengthIm1 = 0.0;

								//	Set first point and get values needed for first step in loop

								if (StartsAtNode || CutoffPtFound){
									PtI = StartPos;
									for (int i = 0; i < 3; ++i){
										ResampledStreamtraces[CornerNum][0][i] = StartPos[i];
										ValsI[i] = PtI[i];
									}
								}
								else{
									for (int i = 0; i < 3; ++i){
										ValsI[i] = STGetFuncts[CornerNum][i](STGetPtrs[CornerNum][i], StartNum);
										PtI[i] = ValsI[i];
										ResampledStreamtraces[CornerNum][0][i] = PtI[i];
									}
								}
								PtIm1 = PtI;
								ValsIm1 = ValsI;
								LgIndex_t OldI = StartNum;

								//TecUtilDataLoadBegin();

								for (int NewI = 1; NewI < NumSTPoints - (int)CutoffPtFound; ++NewI){
									ArcLength += DelLength;

									/*
									*	Move along old path until a segment of DelLength has
									*	been traversed, then interpolate to get new values to
									*	add to new streamtrace zone
									*/
									while (OldI < HighI + 1 && OldI >= LowI && ArcLengthI < ArcLength){
										OldI += Step;

										ArcLengthIm1 = ArcLengthI;
										PtIm1 = PtI;
										ValsIm1 = ValsI;

										for (int i = 0; i < NumVars; ++i){
											ValsI[i] = STGetFuncts[CornerNum][i](STGetPtrs[CornerNum][i], OldI);
											if (i < 3){
												PtI[i] = ValsI[i];
											}
										}

										ArcLengthI += Distance(PtI, PtIm1);
									}

									/*
									*	DelLength has been reached, so time to add a new point
									*	to the new streamtrace zone
									*/
									double Ratio = (ArcLength - ArcLengthIm1) / (ArcLengthI - ArcLengthIm1);
									for (int i = 0; i < NumVars; ++i){
										ResampledStreamtraces[CornerNum][NewI][i] = ValsIm1[i] + Ratio * (ValsI[i] - ValsIm1[i]);
									}
									if (OldI >= HighI + 1 || OldI < LowI){
										while (NewI < NumSTPoints - (int)CutoffPtFound){
											NewI++;
											if (NewI < NumSTPoints){
												for (int i = 0; i < NumVars; ++i){
													ResampledStreamtraces[CornerNum][NewI][i] = ValsIm1[i] + Ratio * (ValsI[i] - ValsIm1[i]);
												}
											}
										}
									}
								}

								if (CutoffPtFound){
									for (int i = 0; i < 3; ++i){
										ResampledStreamtraces[CornerNum][NumSTPoints-1][i] = EndPos[i];
									}
								}

							}
						}
						/*
						 *	Make sure the streamtrace starts at the sphere
						 */

						vec3 BegPt, EndPt, NodePos;
						if (CornerNum < 3){
							for (int i = 0; i < 3; ++i){
								NodePos[i] = p[t[TriNum][CornerNum]][i];
								BegPt[i] = ResampledStreamtraces[CornerNum][0][i];
								EndPt[i] = ResampledStreamtraces[CornerNum][NumSTPoints - 1][i];
							}
							if (Distance(EndPt, NodePos) < Distance(BegPt, NodePos)){
								vector<vector<double> > TmpVector = ResampledStreamtraces[CornerNum];
								for (int i = 0; i < NumSTPoints; ++i){
									ResampledStreamtraces[CornerNum][i] = TmpVector[NumSTPoints - i - 1];
								}
							}
						}
						else if (TriangleHasSaddleCP[TriNum]){
							for (int i = 0; i < 3; ++i){
								BegPt[i] = ResampledStreamtraces[3][0][i];
								EndPt[i] = ResampledStreamtraces[3][NumSTPoints - 1][i];
							}
							NodePos = DoubleCornerNodePos;
							if (Distance(EndPt, NodePos) < Distance(BegPt, NodePos)){
								vector<vector<double> > TmpVector = ResampledStreamtraces[CornerNum];
								for (int i = 0; i < NumSTPoints; ++i){
									ResampledStreamtraces[CornerNum][i] = TmpVector[NumSTPoints - i - 1];
								}
							}
						}
					}
				}

				//	Delete the streamtrace zones
					Set_pa DeleteZoneSet = TecUtilSetAlloc(FALSE);
					for (int i = 0; i < 4 && IsOk; ++i){
						if (STZoneNums[i] > 0){
							IsOk = TecUtilSetAddMember(DeleteZoneSet, STZoneNums[i], FALSE);
						}
					}
					if (IsOk){
						bool DoDel = true;
#if (0)
						vector<int> NoDelList = { 389, 390, 391, 396, 398, 399 };
						for (vector<int>::iterator it = NoDelList.begin(); it != NoDelList.end() && DoDel; it++)
							if (TriNum == *it)
								DoDel = false;
#endif
						if (DoDel)
							IsOk = TecUtilDataSetDeleteZone(DeleteZoneSet);
						else{
							TecUtilZoneSetScatter(SV_SHOW, DeleteZoneSet, 0.0, FALSE);
							TecUtilZoneSetContour(SV_SHOW, DeleteZoneSet, 0.0, FALSE);
							TecUtilZoneSetMesh(SV_SHOW, DeleteZoneSet, 0.0, TRUE);
							TecUtilZoneSetActive(DeleteZoneSet, AssignOp_PlusEquals);
						}
						// 			if (IsOk){
						// 				TecUtilStateChanged(StateChange_ZonesDeleted, (ArbParam_t)DeleteZoneSet);
						// 			}
						TecUtilSetDealloc(&DeleteZoneSet);
					}
					else{
						TecUtilDrawGraphics(TRUE);
						TecUtilInterfaceSuspend(FALSE);
						TecUtilDialogErrMsg("Failed to delete streamtrace zone");
						TecUtilStatusSuspend(FALSE);
						TecUtilStatusFinishPercentDone();
						TecUtilDataLoadEnd();
						TecUtilLockFinish(AddOnID);
						return;
					}

				/*
				 *	Check the 3 or 4 resampled streamtraces to see if they lead
				 *	to different irreducible bundles
				 */

				if (IsOk){
					for (int i = 0; i < 3; ++i){
						Boolean_t IsFound = FALSE;
						for (int j = 0; j < InterIBEdges.size(); ++j){
							if ((InterIBEdges[j][0] == t[TriNum][i] && InterIBEdges[j][1] == t[TriNum][(i + 1) % 3])
								|| (InterIBEdges[j][1] == t[TriNum][i] && InterIBEdges[j][0] == t[TriNum][(i + 1) % 3])){
								IsFound = TRUE;
								break;
							}
						}
						if (!IsFound){
							vec3 PtsI[2], PtsIp1[2];
							int jj = 0;
							for (int j = NumSTPoints - 2; j < NumSTPoints; ++j){
								for (int k = 0; k < 3; ++k){
									PtsI[jj][k] = ResampledStreamtraces[i][j][k];
									PtsIp1[jj][k] = ResampledStreamtraces[(i + 1) % 3][j][k];
								}
								jj++;
							}
							double Dist = Distance(PtsI[1], PtsIp1[1]);
							vec3 VI, VIp1;
							VI = PtsI[1] - PtsI[0];
							VIp1 = PtsIp1[1] - PtsIp1[0];
							double Ang = acosf(dot(VI, VIp1) / (norm(VI) * norm(VIp1))) * 180.0 / PI;
							if (Ang > IBCheckAngle && Dist > IBCheckDistRatio){
								InterIBEdges.push_back(vector<int>(2));
								InterIBEdges[InterIBEdges.size() - 1][0] = t[TriNum][i];
								InterIBEdges[InterIBEdges.size() - 1][1] = t[TriNum][(i + 1) % 3];
							}
						}
					}
				}

				/*
				 *	Now the three corners have been resampled (and concatenated),
				 *	so can create the FEQuad volume
				 */

				if (IsOk){
					IsOk = TecUtilDataSetAddZone("FE Volume", NumNodes, NumElems, 0, ZoneType_FEQuad, NULL);
				}

				TriangleFEZoneNums[TriNum] = TecUtilDataSetGetNumZones();

				vector<FieldData_pa> SetFDPtrs(NumVars, NULL);
				vector<FieldValueSetFunction_pf> SetFuncPtrs(NumVars, NULL);
				if (IsOk){
					for (int i = 0; i < 3 && IsOk; ++i){
						SetFDPtrs[i] = TecUtilDataValueGetWritableRef(TriangleFEZoneNums[TriNum], XYZVarNums[i]);
						IsOk = (VALID_REF(SetFDPtrs[i]));
						if (IsOk){
							SetFuncPtrs[i] = TecUtilDataValueRefGetSetFunc(SetFDPtrs[i]);
							IsOk = (VALID_REF(SetFuncPtrs[i]));
						}
					}
					EntIndex_t VarNum = 3;
					for (int i = 1; i <= NumVars && IsOk; ++i){
						Boolean_t IsXYZVarNum = FALSE;
						for (int j = 0; j < 3 && !IsXYZVarNum; ++j){
							IsXYZVarNum = (XYZVarNums[j] == i);
						}
						if (!IsXYZVarNum){
							SetFDPtrs[VarNum] = TecUtilDataValueGetWritableRef(TriangleFEZoneNums[TriNum], i);
							IsOk = (VALID_REF(SetFDPtrs[VarNum]));
							if (IsOk){
								SetFuncPtrs[VarNum] = TecUtilDataValueRefGetSetFunc(SetFDPtrs[VarNum]);
								IsOk = (VALID_REF(SetFuncPtrs[VarNum]));
							}
							VarNum++;
						}
					}
					if (!IsOk){
						TecUtilDrawGraphics(TRUE);
						TecUtilInterfaceSuspend(FALSE);
						TecUtilDialogErrMsg("Failed to get pointers for creating FE volume");
						TecUtilStatusSuspend(FALSE);
						TecUtilStatusFinishPercentDone();
						TecUtilDataLoadEnd();
						TecUtilLockFinish(AddOnID);
						return;
					}

				}
			
			
				if (IsOk && TriangleHasSaddleCP[TriNum]){
					NodeMap_pa NodeMap;
					if (IsOk){
	// 					int STNum = 0;
	// 					for (int i = 0; i + STNum < 4; ++i){
	// 						for (int j = 0; j < NumSTPoints; ++j){
	// 							for (int k = 0; k < 3; ++k){
	// 								SetFuncPtrs[k](SetFDPtrs[k], i + STNum + 4 * j, ResampledStreamtraces[i][j][k]);
	// 							}
	// 						}
	// 						if (i == DoubleCornerNum){
	// 							STNum++;
	// 							for (int j = 0; j < NumSTPoints; ++j){
	// 								for (int k = 0; k < 3; ++k){
	// 									SetFuncPtrs[k](SetFDPtrs[k], i + STNum + 4 * j, ResampledStreamtraces[3][j][k]);
	// 								}
	// 							}
	// 						}
	// 					}
						int STNum[4];
						switch (DoubleCornerNum){
							case 0:
								STNum[0] = 1;
								STNum[1] = 2;
								STNum[2] = 3;
								STNum[3] = 0;
								break;
							case 1:
								STNum[0] = 2;
								STNum[1] = 0;
								STNum[2] = 3;
								STNum[3] = 1;
								break;
							case 2:
								STNum[0] = 0;
								STNum[1] = 1;
								STNum[2] = 3;
								STNum[3] = 2;
								break;
						}
						for (int i = 0; i < 4; ++i){
							for (int j = 0; j < NumSTPoints; ++j){
								for (int k = 0; k < 3; ++k){
									SetFuncPtrs[k](SetFDPtrs[k], i + 4 * j, ResampledStreamtraces[STNum[i]][j][k]);
								}
							}
						}

						NodeMap = TecUtilDataNodeGetWritableRef(TriangleFEZoneNums[TriNum]);
						IsOk = VALID_REF(NodeMap);
					}
					/*
					*	Remember, node numbers are 1-based
					*/
					int ei = 1;
					for (int i = 1; i <= 4; ++i){
						// 4th corner gets repeated
						TecUtilDataNodeSetByRef(NodeMap, ei, i, i);
					}
					ei++;
					for (int i = 0; i < NumSTPoints - 1 && ei < NumElems - 1; ++i){
						int ii = 1 + 4 * i;
						for (int j = 0; j < 3 && ei < NumElems - 1; ++j){
							int Vals[4] = { ii, ii + 1, ii + 5, ii + 4 };
							for (int k = 0; k < 4; ++k){
								TecUtilDataNodeSetByRef(NodeMap, ei, k + 1, Vals[k]);
							}
							ii++;
							ei++;
						}
						int Vals[4] = { ii, ii - 3, ii + 1, ii + 4 };
						for (int k = 0; k < 4; ++k){
							TecUtilDataNodeSetByRef(NodeMap, ei, k + 1, Vals[k]);
						}
						ei++;
					}
					for (int i = 1; i <= 4; ++i){
						// 4th corner gets repeated
						TecUtilDataNodeSetByRef(NodeMap, ei, i, i + 4 * (NumSTPoints - 1));
					}
				}
				else if (IsOk){
					NodeMap_pa NodeMap;
					if (IsOk){
						for (int i = 0; i < 3; ++i){
							for (int j = 0; j < NumSTPoints; ++j){
								for (int k = 0; k < 3; ++k){
									SetFuncPtrs[k](SetFDPtrs[k], i + 3 * j, ResampledStreamtraces[i][j][k]);
								}
							}
						}

						NodeMap = TecUtilDataNodeGetWritableRef(TriangleFEZoneNums[TriNum]);
						IsOk = VALID_REF(NodeMap);
					}

					/*
					*	Remember, node numbers are 1-based
					*/
					int ei = 1;
					for (int i = 1; i <= 4; ++i){
						// 4th corner gets repeated
						TecUtilDataNodeSetByRef(NodeMap, ei, i, MIN(i, 3));
					}
					ei++;
					for (int i = 0; i < NumSTPoints - 1 && ei < NumElems - 1; ++i){
						int ii = 1 + 3 * i;
						for (int j = 0; j < 2 && ei < NumElems - 1; ++j){
							int Vals[4] = { ii, ii + 1, ii + 4, ii + 3 };
							for (int k = 0; k < 4; ++k){
								TecUtilDataNodeSetByRef(NodeMap, ei, k + 1, Vals[k]);
							}
							ii++;
							ei++;
						}
						int Vals[4] = { ii, ii - 2, ii + 1, ii + 3 };
						for (int k = 0; k < 4; ++k){
							TecUtilDataNodeSetByRef(NodeMap, ei, k + 1, Vals[k]);
						}
						ei++;
					}
					for (int i = 1; i <= 4; ++i){
						// 4th corner gets repeated
						TecUtilDataNodeSetByRef(NodeMap, ei, i, MIN(i + 3 * (NumSTPoints - 1), NumNodes));
					}
				}

				ResampledStreamtraces.clear();
				TecUtilMemoryChangeNotify(-MemoryRequired);

				/*
				 *	Save node and triangle numbers to the fe volume zone's aux data
				 */

				if (IsOk){
					AuxData_pa TempAuxData = TecUtilAuxDataZoneGetRef(TriangleFEZoneNums[TriNum]);
					IsOk = VALID_REF(TempAuxData);
					if (IsOk){
						stringstream ss;
						ss << GBADataPrefix << GBASphereCPName;
						IsOk = TecUtilAuxDataSetItem(TempAuxData, ss.str().c_str(), (ArbParam_t)CPName.c_str(), AuxDataType_String, TRUE);
					}
					if (IsOk){
						stringstream ss;
						ss << GBADataPrefix << GBAZoneType;
						IsOk = TecUtilAuxDataSetItem(TempAuxData, ss.str().c_str(), (ArbParam_t)GBAZoneTypeFEVolumeZone.c_str(), AuxDataType_String, TRUE);
					}
					if (IsOk){
						stringstream ss;
						ss << GBADataPrefix << GBAElemNum;
						IsOk = TecUtilAuxDataSetItem(TempAuxData, ss.str().c_str(), (ArbParam_t)to_string(TriNum + 1).c_str(), AuxDataType_String, TRUE);
					}
					for (int i = 0; i < 3 && IsOk; ++i){
						stringstream ss;
						ss << GBADataPrefix << GBANodeNums[i];
						IsOk = TecUtilAuxDataSetItem(TempAuxData, ss.str().c_str(), (ArbParam_t)to_string(t[TriNum][i] + 1).c_str(), AuxDataType_String, TRUE);
					}
					if (IsOk && TriangleHasConstrainedNode[TriNum]){
						Boolean_t IsFound = FALSE;
						for (int i = 0; i < 3 && !IsFound; ++i){
							for (int j = 0; j < AllIntCPNodeNums.size() && !IsFound; ++j){
								if (t[TriNum][i] == AllIntCPNodeNums[j][2]){
									IsFound = TRUE;
									stringstream ss;
									ss << GBADataPrefix << GBAVolumeCPName;
									IsOk = TecUtilAuxDataSetItem(TempAuxData, ss.str().c_str(), (ArbParam_t)AllIntVolumeCPNames[j].c_str(), AuxDataType_String, TRUE);
								}
							}
						}
					}
				}


				TecUtilDrawGraphics(TRUE);
				if (!SetPercent(TriNum, NumTriangles, ProgressStr.str().c_str())){
					TecUtilDrawGraphics(TRUE);
					TecUtilInterfaceSuspend(FALSE);
					TecUtilStatusFinishPercentDone();
					TecUtilStatusSuspend(FALSE);
					TecUtilDataLoadEnd();
					TecUtilLockFinish(AddOnID);
					return;
				}
				TecUtilDrawGraphics(FALSE);

				TecUtilDataLoadEnd();

			}

			TecUtilDrawGraphics(TRUE);
			TecUtilStatusSuspend(FALSE);
			TecUtilStatusFinishPercentDone();
			TecUtilDrawGraphics(FALSE);

			omp_set_num_threads(numCPU);

			TecUtilDataLoadBegin();

			if (IsOk && InterIBEdges.size() > 0){
				IsOk = TecUtilDataSetAddZone("Inter-IB Edges", InterIBEdges.size(), 1, 1, ZoneType_Ordered, NULL);
				EntIndex_t IBZoneNum;
				FieldData_pa IBXYZFDPtrs[3];
				if (IsOk){
					IBZoneNum = TecUtilDataSetGetNumZones();
					for (int i = 0; i < 3 && IsOk; ++i){
						IBXYZFDPtrs[i] = TecUtilDataValueGetWritableRef(IBZoneNum, XYZVarNums[i]);
						IsOk = VALID_REF(IBXYZFDPtrs[i]);
					}
				}
				if (IsOk){
					for (int i = 0; i < InterIBEdges.size(); ++i){
						vec3 Pts[2];
						for (int j = 0; j < 2; ++j){
							for (int k = 0; k < 3; ++k){
								Pts[j][k] = p[InterIBEdges[i][j]][k];
							}
						}
						vec3 MidPt = (Pts[0] + Pts[1]) * 0.5 - CPPos;

						/*
						* Project midpoint to just above sphere surface
						*/
						// Compute spherical coordinates
						double   rad = norm(MidPt);
						double theta = acos(MidPt[2] / rad);
						double   phi = atan2(MidPt[1], MidPt[0]);

						// Project point onto a sphere of radius "radius"
						MidPt[0] = Radius * 1.01 * sin(theta) * cos(phi);
						MidPt[1] = Radius * 1.01 * sin(theta) * sin(phi);
						MidPt[2] = Radius * 1.01 * cos(theta);

						MidPt = MidPt + CPPos;

						for (int j = 0; j < 3; ++j){
							TecUtilDataValueSetByRef(IBXYZFDPtrs[j], i + 1, MidPt[j]);
						}
					}
					Set_pa IBSet = TecUtilSetAlloc(FALSE);
					TecUtilSetAddMember(IBSet, IBZoneNum, FALSE);

					TecUtilZoneSetScatter(SV_SHOW, IBSet, 0.0, TRUE);
					TecUtilZoneSetContour(SV_SHOW, IBSet, 0.0, FALSE);
					TecUtilZoneSetMesh(SV_SHOW, IBSet, 0.0, FALSE);

					TecUtilZoneSetScatter(SV_COLOR, IBSet, NULL, Black_C);
					TecUtilZoneSetScatter(SV_FRAMESIZE, IBSet, 0.5, NULL);
					TecUtilZoneSetScatterSymbolShape(SV_GEOMSHAPE, IBSet, GeomShape_Sphere);

					TecUtilZoneSetActive(IBSet, AssignOp_PlusEquals);

					TecUtilSetDealloc(&IBSet);
				}
				if (IsOk){
					AuxData_pa TempAuxData = TecUtilAuxDataZoneGetRef(IBZoneNum);
					IsOk = VALID_REF(TempAuxData);
					if (IsOk){
						stringstream ss;
						ss << GBADataPrefix << GBASphereCPName;
						IsOk = TecUtilAuxDataSetItem(TempAuxData, ss.str().c_str(), (ArbParam_t)CPName.c_str(), AuxDataType_String, TRUE);
					}
					if (IsOk){
						stringstream ss;
						ss << GBADataPrefix << GBAZoneType;
						IsOk = TecUtilAuxDataSetItem(TempAuxData, ss.str().c_str(), (ArbParam_t)GBAZoneTypeIBEdgeZone.c_str(), AuxDataType_String, TRUE);
					}
				}
			}
		
			EntIndex_t NumZonesAfterVolumes = TecUtilDataSetGetNumZones();

			Set_pa TmpSet = TecUtilSetAlloc(FALSE);
			vector<Set_pa> ShowSet(MovedPointNums.size(), NULL);
			for (int i = 0; i < MovedPointNums.size(); ++i){
				ShowSet[i] = TecUtilSetAlloc(FALSE);
			}

			for (int i = 0; i < NumTriangles && IsOk; ++i){
				if (TriangleHasConstrainedNode[i] || TriangleHasSaddleCP[i]){
					for (int j = 0; j < MovedPointNums.size() && IsOk; ++j){
						for (int k = 0; k < 3; ++k){
							if (t[i][k] == MovedPointNums[j]){
								IsOk = TecUtilSetAddMember(ShowSet[j], TriangleFEZoneNums[i], FALSE);
							}
						}
					}
				}
				IsOk = TecUtilSetAddMember(TmpSet, TriangleFEZoneNums[i], FALSE);
			}

			NumVars = ActualNumVars;

			TriangleFEZoneNums.clear();
			NodeHasSaddleCP.clear();
			NodeIsConstrained.clear();
			TriangleHasConstrainedNode.clear();
			TriangleHasSaddleCP.clear();
			MemoryRequired = (sizeof(bool) * (2*NumPoints + 2*NumTriangles) + sizeof(EntIndex_t) * NumTriangles) / 1024;
			TecUtilMemoryChangeNotify(-MemoryRequired);

			if (IsOk){
				// Grouping the nodes neighboring a GP node
				ColorIndex_t TmpColor = 0;
			
				TecUtilZoneSetScatter(SV_SHOW, TmpSet, 0.0, FALSE);
				TecUtilZoneSetContour(SV_SHOW, TmpSet, 0.0, TRUE);
				TecUtilZoneSetMesh(SV_SHOW, TmpSet, 0.0, TRUE);
				TecUtilZoneSetMesh(SV_COLOR, TmpSet, 0.0, Red_C);
				TecUtilZoneSetActive(TmpSet, AssignOp_MinusEquals);

				ArgList_pa CurrentArgList = TecUtilArgListAlloc();
				TecUtilArgListAppendString(CurrentArgList, SV_P1, SV_FIELDMAP);
				TecUtilArgListAppendString(CurrentArgList, SV_P2, SV_EFFECTS);
				TecUtilArgListAppendString(CurrentArgList, SV_P3, SV_USETRANSLUCENCY);
				TecUtilArgListAppendSet(CurrentArgList, SV_OBJECTSET, TmpSet);
				TecUtilArgListAppendArbParam(CurrentArgList, SV_IVALUE, TRUE);
				TecUtilStyleSetLowLevelX(CurrentArgList);
				TecUtilArgListClear(CurrentArgList);

				TecUtilArgListAppendString(CurrentArgList, SV_P1, SV_FIELDMAP);
				TecUtilArgListAppendString(CurrentArgList, SV_P2, SV_EFFECTS);
				TecUtilArgListAppendString(CurrentArgList, SV_P3, SV_TRANSLUCENCY);
				TecUtilArgListAppendSet(CurrentArgList, SV_OBJECTSET, TmpSet);
				TecUtilArgListAppendArbParam(CurrentArgList, SV_IVALUE, 80);
				TecUtilStyleSetLowLevelX(CurrentArgList);
				TecUtilArgListClear(CurrentArgList);

				// Setting unique colors for GP nodes and their neighbors
				for (int i = 0; i < MovedPointNums.size(); ++i){
					TecUtilZoneSetMesh(SV_COLOR, ShowSet[i], 0.0, i % 63);
					TecUtilZoneSetActive(ShowSet[i], AssignOp_MinusEquals);

					//	Modify group number
					CurrentArgList = TecUtilArgListAlloc();
					TecUtilArgListAppendString(CurrentArgList, SV_P1, SV_FIELDMAP);
					TecUtilArgListAppendString(CurrentArgList, SV_P2, SV_GROUP);
					TecUtilArgListAppendSet(CurrentArgList, SV_OBJECTSET, ShowSet[i]);
					TecUtilArgListAppendArbParam(CurrentArgList, SV_IVALUE, (LgIndex_t)(50 + GroupNum++));
					TecUtilStyleSetLowLevelX(CurrentArgList);
					TecUtilArgListDealloc(&CurrentArgList);

					TecUtilSetDealloc(&ShowSet[i]);
				}
				TecUtilSetDealloc(&TmpSet);
				TmpSet = TecUtilSetAlloc(FALSE);
				TecUtilSetAddMember(TmpSet, 1, FALSE);
				TecUtilZoneSetActive(TmpSet, AssignOp_MinusEquals);
			}
			TecUtilSetDealloc(&TmpSet);

			TecUtilDataLoadEnd();
		}

		delete p, t;
		for (int i = 0; i < NumEdges; ++i){
			delete[] e[i];
		}
		delete e;
		MemoryRequired = sizeof(point) * NumPoints
			+ sizeof(triangle) * NumTriangles
			+ sizeof(int) * 2 * NumEdges;
		TecUtilMemoryChangeNotify(-MemoryRequired / 1024);

		/*
		 *	Have tecplot unload data for everything created, but only
		 *	if there are more CPs to process.
		 */

		if (NumSelectedCPs > 1 && SelectCPNum < NumSelectedCPs - 1){
			TecUtilDrawGraphics(TRUE);
			TecUtilPleaseWait(string(ProgressStr.str() + string(" ... Unloading Data")).c_str(), TRUE);
			TecUtilDrawGraphics(FALSE);
			EntIndex_t NewNumZones = TecUtilDataSetGetNumZones() - NumZones;
			for (EntIndex_t i = NumZones + 1; i <= NewNumZones && IsOk; ++i){
				for (EntIndex_t j = 1; j <= NumVars && IsOk; ++j){
					IsOk = TecUtilDataValueUnload(i, j);
				}
			}
			TecUtilDrawGraphics(TRUE);
			TecUtilPleaseWait(NULL, FALSE);
			TecUtilDrawGraphics(FALSE);
		}

	}



	Time2 = high_resolution_clock::now();

	duration<double> TimeSpan = duration_cast<duration<double>>(Time2 - Time1);
	TecUtilDialogMessageBox(string(string("Runtime: ") + to_string(TimeSpan.count()) + string(" seconds.")).c_str(), MessageBoxType_Information);


	TecUtilArrayDealloc((void **)&CPNums);

	TecUtilDrawGraphics(TRUE);
	TecUtilInterfaceSuspend(FALSE);
	TecUtilLockFinish(AddOnID);
	return;
}

Boolean_t SetPercent(unsigned int CurrentNum, unsigned int TotalNum, const char* ProgresssText){
	unsigned int Percent = MIN((int)((double)CurrentNum / (double)TotalNum * 100.), 100);

	Boolean_t IsOk = TRUE;
	stringstream ss;
	ss << ProgresssText << "  (" << Percent << "% Complete)";

	TecUtilStatusSuspend(FALSE);
	TecUtilStatusSetPercentDoneText(ss.str().c_str());
	IsOk = TecUtilStatusCheckPercentDone(Percent);
	TecUtilStatusSuspend(TRUE);

	return IsOk;
}