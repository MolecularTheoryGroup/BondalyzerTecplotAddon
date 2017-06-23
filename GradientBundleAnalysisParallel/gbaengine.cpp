
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
#include <fstream>
#include <algorithm>

// for profiling
#include <chrono>
#include <ctime>
#include <ratio>

#include <omp.h>

#include "CSM_DATA_SET_INFO.h"
#include "CSM_DATA_TYPES.h"
#include "meshgen2d_sphere.h"
#include "VIEWRESULTS.h"
#include "CSM_FE_VOLUME.h"
#include "CSM_GRAD_PATH.h"
#include "CSM_CRIT_POINTS.h"
#include "CSM_GUI.h"

#include "GBAENGINE.h"

#include <armadillo>
using namespace arma;

using std::string;
using std::to_string;
using std::stoi;
using std::stringstream;
using std::vector;
using std::ofstream;
using std::ios;
using std::endl;

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


	TecUtilLockStart(AddOnID);
	TecUtilDrawGraphics(FALSE);

	Boolean_t DoTiming = TRUE;

	high_resolution_clock::time_point Time1, Time2;
	if (DoTiming){
		Time1 = high_resolution_clock::now();
	}


	Boolean_t IsOk = TRUE;
	int MemoryRequired;
	string TmpString;

	vector<string> GoodSGPEndCPTypes = { "Atom", "Bond" };

	EntIndex_t NumVars = TecUtilDataSetGetNumVars();

	/*
	*	Get min and max XYZ for the system, just in case we need it.
	*/


	vector<EntIndex_t> XYZVarNums(3);
	TecUtilAxisGetVarAssignments(&XYZVarNums[0], &XYZVarNums[1], &XYZVarNums[2]);
	for (int i = 0; i < 3 && IsOk; ++i)
		IsOk = (XYZVarNums[i] > 0);

	vec3 VolMinXYZ, VolMaxXYZ;
	vector<int> VolMaxIJK(3);
	EntIndex_t VolZoneNum = ZoneNumByName("Full Volume");
	if (IsOk){
		if (VolZoneNum <= 0){
			VolZoneNum = 1;
			// 			TecUtilDialogErrMsg("Couldn't get volume zone");
			//TecUtilLockFinish(AddOnID);
			// 			return;
		}
// 		for (int i = 0; i < 3; ++i){
// 			TecUtilDataValueGetMinMaxByZoneVar(VolZoneNum, XYZVarNums[i], &VolMinXYZ[i], &VolMaxXYZ[i]);
// 		}
		ZoneXYZVarGetMinMax_Ordered3DZone(XYZVarNums, VolZoneNum, VolMinXYZ, VolMaxXYZ);
		TecUtilZoneGetIJK(VolZoneNum, &VolMaxIJK[0], &VolMaxIJK[1], &VolMaxIJK[2]);
	}

	vector<FieldDataType_e> VarDataTypes(NumVars);
	vector<ValueLocation_e> VarLocations(NumVars, ValueLocation_Nodal);
	
	for (int i = 0; i < NumVars; ++i){
		VarDataTypes[i] = TecUtilDataValueGetType(VolZoneNum, i + 1);
	}


	/*
	*	Check to see if there's any existing GBA zones present.
	*	If there are, then the data types and data locations need
	*	to be consistent with those zones
	*	i.e. need to use cell-centered data for integrated variables.
	*
	*	Just check for a GBA sphere zone, since that's sufficient.
	*/

	for (int i = 1; i <= TecUtilDataSetGetNumZones(); ++i){
		if (AuxDataZoneHasItem(i, CSMAuxData.GBA.ZoneType)){
			for (int j = 3; j < NumVars; ++j){
				VarDataTypes[j] = TecUtilDataValueGetType(i, j + 1);
				VarLocations[j] = TecUtilDataValueGetLocation(i, j + 1);
			}
			break;
		}
	}

	int GroupNum = 0;

	/*
	 *	Get job parameters from dialog
	 */
	Boolean_t UseCutoff = TRUE;
	UseCutoff = TecGUIToggleGet(TGLOpenSys_TOG_T1_1);
	double CutoffVal = 0.001;
	TecGUITextFieldGetDouble(TFCutoff_TF_T1_1, &CutoffVal);
	EntIndex_t CutoffVarNum = VarNumByName(string("Electron Density"));
	LgIndex_t NumEdgeGPs = TecGUIScaleGetValue(SCNumEdgeGPs_SC_T1_1);

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
// 	double IBCheckDistRatio = 0.2;
// 	TecGUITextFieldGetDouble(TFIBDist_TF_T1_1, &IBCheckDistRatio);
// 	double IBCheckAngle = 20.0;
// 	TecGUITextFieldGetDouble(TFIBAng_TF_T1_1, &IBCheckAngle);

	LgIndex_t * CPNums = NULL;
	LgIndex_t NumSelectedCPs = 1;
	TecGUIListGetSelectedItems(MLSelCPs_MLST_T1_1, &CPNums, &NumSelectedCPs);
	/*
	*	Set sphere radius and mesh refinement level
	*	(will be set by user in future)
	*/


	RadMode RadiusMode = MINCPDISTRATIO;
	RadiusMode = (RadMode)TecGUIRadioBoxGetToggle(RBRadMode_RADIO_T1_1);
	double UserRadius = 0.25;
	TecGUITextFieldGetDouble(TFRad_TF_T1_1, &UserRadius);
	int Level = 3;
	TecGUITextFieldGetLgIndex(TFLevel_TFS_T1_1, &Level);
	LgIndex_t NumSTPoints = 100;
	TecGUITextFieldGetLgIndex(TFSTPts_TF_T1_1, &NumSTPoints);



	EntIndex_t OldNumZones = TecUtilDataSetGetNumZones();	

	if (!IsOk){
		TecUtilDrawGraphics(TRUE);
		TecUtilDialogErrMsg("Couldn't get XYZ vars.");
		TecUtilLockFinish(AddOnID);
		return;
	}

	vector<int> NumCPs;
	for (int SelectCPNum = 0; SelectCPNum < NumSelectedCPs && IsOk; ++SelectCPNum){

		LgIndex_t CPNum = 1;
		int CPType = -3;

		double Radius = 0.25;

		stringstream ProgressStr;
		string CPString, CPName;
		bool IsCP;

		int RealCPNum = 0;
		EntIndex_t CPTypeVarNum;
		EntIndex_t CPZoneNum;
		LgIndex_t CPIJKMax[3];

		vec3 CPPos;
		if (IsOk){
			CPPos = GetCoordsFromListItem(CPNums[SelectCPNum], MLSelCPs_MLST_T1_1, &CPString, &CPName, &CPNum, &CPZoneNum, &IsCP, &NumCPs);

			if (IsCP){

// 				CPZoneNum = ZoneNumByName("Critical Points");
				CPTypeVarNum = VarNumByName("CritPointType");

// 				LgIndex_t NumCPs[4];
// 				if (IsOk){
// 					for (int i = 0; i < 4; ++i){
// 						NumCPs[i] = stoi(AuxDataZoneGetItem(CPZoneNum, CCDataNumCPs[i]));
// 					}
// 				}
// 

				TecUtilZoneGetIJK(CPZoneNum, &CPIJKMax[0], &CPIJKMax[1], &CPIJKMax[2]);
				IsOk = (CPIJKMax[0] > 1 && CPIJKMax[1] == 1 && CPIJKMax[2] == 1);
				if (!IsOk){
					TecUtilDrawGraphics(TRUE);
					TecUtilDialogErrMsg("Critical point zone dimensions incorrect.");
					TecUtilLockFinish(AddOnID);
					return;
				}

				/*
				*	Get XYZ readable pointers to the crit points zone to get the
				*	position of the fifth atom (working with cubane, so this is the first carbon).
				*/

// 				if (CPZoneNum <= 0){
// 					TecUtilDrawGraphics(TRUE);
// 					TecUtilDialogErrMsg("Couldn't get critical points zone.");
// 					TecUtilLockFinish(AddOnID);
// 					return;
// 				}

				for (int i = 0; i < 4; ++i){
					if (CPName == RankStrs[i]){
						int RealCPNum = 0;
						for (int j = 0; j < i; ++j){
							RealCPNum += NumCPs[j];
						}
						CPNum = RealCPNum + CPNum;
						break;
					}
				}


				CPType = (int)TecUtilDataValueGetByZoneVar(CPZoneNum, CPTypeVarNum, CPNum);
			}

			ProgressStr << "Processing " << CPString;
			if (NumSelectedCPs > 1){
				ProgressStr << " ... CP " << SelectCPNum + 1 << " of " << NumSelectedCPs;
			}

// 			TecUtilStringDealloc(&TempCStr);
		}
		else{
			TecUtilDialogErrMsg("Didn't get CP number correctly");
		}

		if (IsOk){
			TmpString = ProgressStr.str() + string(" ... Step 1 of 4 ... Generating mesh");
			TecUtilDrawGraphics(TRUE);
			TecUtilPleaseWait(TmpString.c_str(), TRUE);
			TecUtilDrawGraphics(FALSE);
		}

// 		for (int i = 0; i < 3; ++i){
// 			CPPos[i] = TecUtilDataValueGetByZoneVar(CPZoneNum, XYZVarNums[i], CPNum);
// 		}

		vector<point> IntersectionPoints;
		vector<vector<LgIndex_t> > IntZoneSaddleCPNodeNums;
		vector<vector<LgIndex_t> > AllIntCPNodeNums;
		vector<string> AllIntVolumeCPNames;

		vector<vec3> IntCPPos;

		EntIndex_t NumZones = TecUtilDataSetGetNumZones();

		point* p = NULL;
		triangle* t = NULL;
		int** e = NULL;
		int NumPoints, NumTriangles, NumEdges;

		vector<int> MovedPointNums;
		string MovedPointCPTypes, MovedPointCPNames;

		MeshStatus_e MeshStatus = FAIL_INVALIDCONSTRAINT;


		double ClosestCPDist = 1e50;
		if (IsOk){
			if (IsCP){
				for (int i = 1; i <= CPIJKMax[0] && IsOk; ++i){
					if (i != CPNum && (int)TecUtilDataValueGetByZoneVar(CPZoneNum, CPTypeVarNum, i) != CPType){
						vec3 OtherCP;
						for (int j = 0; j < 3; ++j){
							OtherCP[j] = TecUtilDataValueGetByZoneVar(CPZoneNum, XYZVarNums[j], i);
						}
						ClosestCPDist = MIN(ClosestCPDist, Distance(OtherCP, CPPos));
					}
				}
			}
			else{
				for (int ZoneNum = 1; ZoneNum <= TecUtilDataSetGetNumZones(); ++ZoneNum){
					char* ZoneNameCStr;
					if (TecUtilZoneIsOrdered(ZoneNum) && TecUtilZoneGetName(ZoneNum, &ZoneNameCStr)){
						string ZoneName = ZoneNameCStr;
						TecUtilStringDealloc(&ZoneNameCStr);

						int ElemNum = SearchVectorForString(ElementSymbolList, ZoneName, false);
						if (ElemNum < 0)
							ElemNum = SearchVectorForString(ElementNameList, ZoneName);
						if (ElemNum >= 0){
							int NumPts;
							TecUtilZoneGetIJK(ZoneNum, &NumPts, NULL, NULL);
							for (int i = 1; i <= NumPts; ++i){
								if (ZoneNum != CPZoneNum || i != CPNum){
									vec3 OtherCP;
									for (int j = 0; j < 3; ++j){
										OtherCP[j] = TecUtilDataValueGetByZoneVar(ZoneNum, XYZVarNums[j], i);
									}
									ClosestCPDist = MIN(ClosestCPDist, Distance(OtherCP, CPPos));
								}
							}
						}
					}
				}
			}
		}

		if (RadiusMode == ABSOLUTERADIUS || ClosestCPDist == 1e50){
			Radius = UserRadius;
		}
		else if (RadiusMode == MINCPDISTRATIO){
			Radius = UserRadius * ClosestCPDist;
		}

		vector<int> XYZRhoVarNums = XYZVarNums;
		XYZRhoVarNums.push_back(CutoffVarNum);

		vector<GradPath_c> IntGPs;


		/*   
		 *	Find all volumetric special gradient paths that will intersect
		 *	the sphere.
		 */
		TecUtilDataLoadBegin();
		for (EntIndex_t CurZoneNum = 1; CurZoneNum <= OldNumZones && IsOk && IsCP; ++CurZoneNum){
			LgIndex_t TmpIJK[3];

			Boolean_t ZoneOK = TRUE;
			Boolean_t ZoneIsActive = TecUtilZoneIsActive(CurZoneNum);
			if (!ZoneIsActive){
				Set_pa TempSet = TecUtilSetAlloc(FALSE);
				TecUtilSetAddMember(TempSet, CurZoneNum, FALSE);
				TecUtilZoneSetActive(TempSet, AssignOp_PlusEquals);
				TecUtilSetDealloc(&TempSet);
			}
			ZoneOK = !AuxDataZoneHasItem(CurZoneNum, CSMAuxData.GBA.ZoneType);
			if (ZoneOK)
				TecUtilZoneGetIJK(CurZoneNum, &TmpIJK[0], &TmpIJK[1], &TmpIJK[2]);

			if (ZoneOK && TmpIJK[0] > 1 && TmpIJK[1] == 1 && TmpIJK[2] == 1){
				// If 1-d (line) zone that wasn't made by GBA.
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
				if (ZoneOK){
					for (int i = 0; i < 2 && ZoneOK; ++i){
						ZoneOK = AuxDataZoneGetItem(CurZoneNum, CCDataGPEndNums[i], TmpString);
						if (ZoneOK){
							StartEndCPs[i] = stoi(TmpString);
							ZoneOK = AuxDataZoneGetItem(CurZoneNum, CCDataGPEndTypes[i], TmpString);
							if (ZoneOK)
								ZoneOK = (SearchVectorForString(GoodSGPEndCPTypes, TmpString) >= 0);
							if (ZoneOK)
								StartEndTypes[i] = TmpString;
						}
					}
				}

				Boolean_t IntValid = TRUE;
				
				if (ZoneOK){
					Boolean_t IntRecorded = FALSE;
					for (int i = 0; i < 2 && IsOk; ++i){
						if (StartEndCPs[i] == CPNum){
							ConnectsToCP = TRUE;
							StartsFrom1 = (i == 0);
							for (int j = 1; j < 4 && IsOk && IntValid; ++j){
								if (StartEndTypes[(i + 1) % 2] == RankStrs[j]){
									GradPath_c TmpGP(CurZoneNum, XYZRhoVarNums, AddOnID);
									if (Distance(CPPos, TmpGP[0]) < Distance(CPPos, TmpGP[-1])){
										TmpGP.Reverse();
									}
									IsOk = TmpGP.Trim(CPPos, Radius);
									if (!IsOk){
										TecUtilDialogErrMsg("Failed to trim existing GP");
									}

									/*
									 *	Now check the intersection point against the previously found intersections
									 *	to make sure there are not two intersections at the same point
									 */

									for (const GradPath_c & OldGP : IntGPs){
										double IntIntDist = Distance(TmpGP[-1], OldGP[-1]);
										if (IntIntDist < 0.1){
											IntValid = FALSE;
										}
// 										double a = 1;
									}

									if (IntValid){

										IntGPs.push_back(TmpGP);

										vec3 TmpCP;
										for (int k = 0; k < 3 && IsOk; ++k){
											TmpCP[k] = TecUtilDataValueGetByZoneVar(CPZoneNum, XYZVarNums[k], StartEndCPs[(i + 1) % 2]);
										}
										IntZoneSaddleCPNodeNums.push_back(vector<LgIndex_t>());
										IntZoneSaddleCPNodeNums[IntZoneSaddleCPNodeNums.size() - 1].push_back(static_cast<int>(IntersectionPoints.size()));
										IntZoneSaddleCPNodeNums[IntZoneSaddleCPNodeNums.size() - 1].push_back(CurZoneNum);
										IntZoneSaddleCPNodeNums[IntZoneSaddleCPNodeNums.size() - 1].push_back(StartEndCPs[(i + 1) % 2]);
										IntCPPos.push_back(TmpCP);
									}

// 									if (TmpIJK[0] < NumSTPoints){
// 										NumSTPoints = TmpIJK[0];
// 										if (NumSTPoints % 2 != 0){
// 											NumSTPoints--;
// 										}
// 									}
								}
								if (IntValid){
									int CPCount = 0;
									for (int k = 0; k < 4; ++k){
										if (StartEndTypes[(i + 1) % 2] == RankStrs[k] && !IntRecorded){
											AllIntVolumeCPNames.push_back(StartEndTypes[(i + 1) % 2] + string(" ") + to_string(StartEndCPs[(i + 1) % 2] - CPCount));
											AllIntCPNodeNums.push_back(vector<LgIndex_t>());
											AllIntCPNodeNums[AllIntCPNodeNums.size() - 1].push_back(static_cast<int>(IntersectionPoints.size()));
											AllIntCPNodeNums[AllIntCPNodeNums.size() - 1].push_back(StartEndCPs[(i + 1) % 2]);
											Boolean_t ZoneFound = FALSE;
											if (StartEndTypes[(i + 1) % 2] == RankStrs[1]){
												// if the intersecting SGP connects to a bond point,
												// then find the CP on the other side of the bond path.

												// The other end of the bond path should be the zone before of after CurZoneNum, so check those first
												for (int CheckNum = -1; CheckNum < 2; CheckNum += 2){
													if (CurZoneNum + CheckNum <= NumZones && AuxDataZoneGetItem(CurZoneNum + CheckNum, CCDataGPEndNums[0], TmpString)){
														// the beginning point for bond paths is always the bond point, so only check CCDataGPEndNums[0]
														if (stoi(TmpString) == StartEndCPs[(i + 1) % 2]){
															// Matching bond path found, so get nuclear CP on other side
															int CheckNucCPNum = stoi(AuxDataZoneGetItem(CurZoneNum + CheckNum, CCDataGPEndNums[1]));
															// now search through the nuclear zones, if present, to see what type of element the CP is
															if (CheckNucCPNum > 0){
																double MinNucCPCheckDist = 1e150, TempDist;
																string ClosestElementName;
																int ClosestAtomNum = -1;
																for (int CheckZoneNum = 1; CheckZoneNum <= NumZones; ++CheckZoneNum){
																	if (CheckZoneNum != CurZoneNum
																		&& CheckZoneNum != CurZoneNum + CheckNum
																		&& TecUtilZoneIsOrdered(CheckZoneNum)
																		&& AuxDataZoneItemMatches(CheckZoneNum, DLZoneType, DLZoneTypeNuclearPositions)){
																		// need to now find the 
																		int NucIJK[3];
																		TecUtilZoneGetIJK(CheckZoneNum, &NucIJK[0], &NucIJK[1], &NucIJK[2]);
																		vec3 NucPt, NucCPPt;
																		for (int Dir = 0; Dir < 3; ++Dir){
																			NucCPPt[Dir] = TecUtilDataValueGetByZoneVar(CurZoneNum + CheckNum, XYZVarNums[Dir], CheckNucCPNum);
																		}
																		for (int NucCheckPt = 1; NucCheckPt <= NucIJK[0]; ++NucCheckPt){
																			for (int Dir = 0; Dir < 3; ++Dir){
																				NucPt[Dir] = TecUtilDataValueGetByZoneVar(CheckZoneNum, XYZVarNums[Dir], NucCheckPt);
																			}
																			TempDist = DistSqr(NucPt, NucCPPt);
																			if (TempDist < MinNucCPCheckDist){
																				MinNucCPCheckDist = TempDist;
																				ClosestElementName = AuxDataZoneGetItem(CheckZoneNum, DLZoneAtomicSpecies);
																				ClosestAtomNum = NucCheckPt;
																			}
																		}
																	}
																}
																if (ClosestAtomNum > 0){
																	MovedPointCPNames += ClosestElementName + " " + to_string(ClosestAtomNum) + ",";
																	ZoneFound = TRUE;
																}
															}
														}
													}
												}
											}
											MovedPointCPTypes += StartEndTypes[(i + 1) % 2] + ",";
											if (!ZoneFound)
												MovedPointCPNames += StartEndTypes[(i + 1) % 2] + string(" ") + to_string(StartEndCPs[(i + 1) % 2] - CPCount) + ",";
											IntRecorded = TRUE;
											break;
										}
										CPCount += NumCPs[k];
									}
								}
							}
							break;
						}
					}
				}
				if (ConnectsToCP && IsOk && IntValid){

					// An intersection exists, so need to find its location. 
					// INTERSECTION LOCATION IS THE SAME AS THE TRIMMED (LAST) POINT OF IntGPs[-1]
					// 
					point IntPt;
					for (int i = 0; i < 3; ++i){
						IntPt[i] = IntGPs[IntGPs.size() - 1][-1][i] - CPPos[i];
					}
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
// 						TODO: The distance cutoff of 0.1 is completely arbitrary. Make it less so.
						IntersectionGood = (Distance(TmpPt1, TmpPt2) >= 0.1);
					}

					if (IntersectionGood)
						IntersectionPoints.push_back(IntPt);

// 					Below is the (spurious) code to find the intersection of the GP with the sphere.
// 					This is unnecessary because the point was already found when the GP was trimmed above.

// 					vec3 Pt1, Pt2;
// 					LgIndex_t StartPt = 1, Step = 1;
// 					if (!StartsFrom1){
// 						StartPt = TmpIJK[0];
// 						Step = -1;
// 					}
// 					Boolean_t IntersectionFound = FALSE;
// 					double Dist1, Dist2;
// 					for (int i = 0; i < 3; ++i){
// 						Pt1[i] = TecUtilDataValueGetByRef(TmpXYZPtrs[i], StartPt + Step);
// 					}
// 					Pt2 = Pt1;
// 					Dist2 = Dist1 = Distance(Pt1, CPPos);
// 					StartPt += Step;
// 					for (LgIndex_t PtNum = StartPt + Step; PtNum >= 1 && PtNum <= TmpIJK[0] && !IntersectionFound && IsOk; PtNum += Step){
// 						for (int i = 0; i < 3; ++i){
// 							Pt1[i] = TecUtilDataValueGetByRef(TmpXYZPtrs[i], PtNum);
// 						}
// 						Dist1 = Distance(Pt1, CPPos);
// 						if ((Dist1 < Radius && Dist2 >= Radius) ||
// 							(Dist1 >= Radius && Dist2 < Radius)){
// 							// Now Pt1 and Pt2 are straddling the sphere.
// 							// Interpolate to find the intersection point.
// 							double DistanceRatio = (Radius - Dist1) / (Dist2 - Dist1);
// 							point IntPt;
// 							for (int i = 0; i < 3; ++i){
// 								IntPt[i] = (Pt1[i] + DistanceRatio * (Pt2[i] - Pt1[i])) - CPPos[i];
// 							}
// 
// 							/*
// 							* Project intersection point to sphere surface
// 							*/
// 							// Compute spherical coordinates
// 							double   rad = sqrt(pow(IntPt.x, 2) + pow(IntPt.y, 2) + pow(IntPt.z, 2));
// 							double theta = acos(IntPt.z / rad);
// 							double   phi = atan2(IntPt.y, IntPt.x);
// 
// 							// Project point onto a sphere of radius "radius"
// 							IntPt.x = Radius * sin(theta) * cos(phi);
// 							IntPt.y = Radius * sin(theta) * sin(phi);
// 							IntPt.z = Radius * cos(theta);
// 
// 
// 							IntersectionFound = TRUE;
// 							/*
// 							 *	Make sure this intersection isn't too close to another
// 							 */
// 
// 							Boolean_t IntersectionGood = TRUE;
// 							for (int i = 0; i < IntersectionPoints.size() && IntersectionGood; ++i){
// 								vec3 TmpPt1, TmpPt2;
// 								for (int j = 0; j < 3; ++j){
// 									TmpPt1[j] = IntPt[j];
// 									TmpPt2[j] = IntersectionPoints[i][j];
// 								}
// 								IntersectionGood = (Distance(TmpPt1, TmpPt2) >= 0.1);
// 							}
// 
// 							if (IntersectionGood)
// 								IntersectionPoints.push_back(IntPt);
// 
// 
// 						}
// 						Pt2 = Pt1;
// 						Dist2 = Dist1;
// 					}
				}
			}
			if (!ZoneIsActive && IsOk){
				Set_pa TmpSet = TecUtilSetAlloc(FALSE);
				IsOk = TecUtilSetAddMember(TmpSet, CurZoneNum, FALSE);
				if (IsOk)
					IsOk = (TecUtilZoneSetActive(TmpSet, AssignOp_MinusEquals) <= 1);
				TecUtilSetDealloc(&TmpSet);
				if (!IsOk)
					TecUtilDialogErrMsg("Failed to deactivate zone during intersection check");
			}
		}
		TecUtilDataLoadEnd();


		/*
		 *	Get optimized spherical mesh
		 */

		int MaxRefine = 2;
		int RefineNum = 0;
		TecUtilDataLoadBegin();
		while (MeshStatus != 0){
			/*
			 *	Get memory usage estimate for mesh
			 */
			NumTriangles = 20 * (int)pow(4, Level);
			NumPoints = NumTriangles / 2 + 2;
			NumEdges = int(NumTriangles * double(3.0 / 2.0));
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


			if (MeshStatus == FAIL_NEEDREFINEMENT){
				//break;
				if (RefineNum >= MaxRefine){
					break;
				}
				else{
					//TecUtilDialogMessageBox("Mesh not fine enough for constrained points. Refining mesh.", MessageBoxType_Warning);
					TecUtilMemoryChangeNotify(-((int)sizeof(point) * NumPoints
						+ (int)sizeof(triangle) * NumTriangles) / 1024);
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
				TecUtilDataLoadEnd();
				TecUtilDialogErrMsg("Constrained point too far from mesh");
				TecUtilLockFinish(AddOnID);
				return;
			}
			else if (MeshStatus == SUCCESS_NOTCONVERGED){
				// 				TecUtilDialogMessageBox("Mesh did not converge.", MessageBoxType_Warning);
				break;
			}

		}
		TecUtilDataLoadEnd();

		IsOk = (MeshStatus == SUCCESS || MeshStatus == SUCCESS_NOTCONVERGED && NumPoints > 0 && NumTriangles > 0);
		if (!IsOk){
			TecUtilDrawGraphics(TRUE);
			TecUtilDialogErrMsg("Failed to make mesh.");
			TecUtilLockFinish(AddOnID);
			return;
		}

		/*
		 *	Get a list of edge numbers for each triangle
		 */
		vector<vector<int> > TriangleEdgeNums(NumTriangles), NodeEdgeNums(NumPoints);
		for (auto & i : TriangleEdgeNums) i.reserve(3);
#pragma omp parallel for
		for (int TriNum = 0; TriNum < NumTriangles; ++TriNum){
			for (int tNode = 0; tNode < 3; ++tNode){
				bool EdgeFound = false;
				for (int EdgeNum = 0; EdgeNum < NumEdges && !EdgeFound; ++EdgeNum){
					for (int eNode = 0; eNode < 2; ++eNode){
						if (e[EdgeNum][eNode] == t[TriNum][tNode] && e[EdgeNum][(eNode + 1) % 2] == t[TriNum][(tNode + 1) % 3]){
							TriangleEdgeNums[TriNum].push_back(EdgeNum);
							EdgeFound = true;
						}
					}
				}
			}
		}
		for (auto & i : NodeEdgeNums) i.reserve(6);
#pragma omp parallel for
		for (int NodeNum = 0; NodeNum < NumPoints; ++NodeNum){
			for (int EdgeNum = 0; EdgeNum < NumEdges; ++EdgeNum){
				for (int eNode = 0; eNode < 2; ++eNode){
					if (e[EdgeNum][eNode] == NodeNum && std::find(NodeEdgeNums[NodeNum].begin(), NodeEdgeNums[NodeNum].end(), EdgeNum) == NodeEdgeNums[NodeNum].end()){
						NodeEdgeNums[NodeNum].push_back(EdgeNum);
					}
				}
			}
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

			if (IsOk){
				ArgList_pa ArgList = TecUtilArgListAlloc();
				TecUtilArgListAppendString(ArgList, SV_NAME, string("Sphere (" + CPString + ")").c_str());
				TecUtilArgListAppendInt(ArgList, SV_ZONETYPE, ZoneType_FETriangle);
				TecUtilArgListAppendInt(ArgList, SV_IMAX, NumPoints);
				TecUtilArgListAppendInt(ArgList, SV_JMAX, NumTriangles);
				TecUtilArgListAppendArray(ArgList, SV_VALUELOCATION, VarLocations.data());
				TecUtilArgListAppendArray(ArgList, SV_VARDATATYPE, VarDataTypes.data());
				IsOk = TecUtilDataSetAddZoneX(ArgList);
				TecUtilArgListDealloc(&ArgList);
			}

// 			IsOk = TecUtilDataSetAddZone("Sphere", NumPoints, NumTriangles, 0, ZoneType_FETriangle, NULL);
			EntIndex_t SurfZoneNum = TecUtilDataSetGetNumZones();
			FieldData_pa SurfXYZPtrs[3];

			for (int i = 0; i < 3 && IsOk; ++i){
				SurfXYZPtrs[i] = TecUtilDataValueGetWritableRef(SurfZoneNum, XYZVarNums[i]);
				IsOk = VALID_REF(SurfXYZPtrs[i]);
			}

			if (IsOk){
				NodeMap_pa NodeMap;
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
				IsOk = AuxDataZoneSetItem(SurfZoneNum, CSMAuxData.GBA.SphereCPName, CPString);
				if (IsOk){
					IsOk = AuxDataZoneSetItem(SurfZoneNum, CSMAuxData.GBA.ZoneType, CSMAuxData.GBA.ZoneTypeSphereZone);
				}
				if (IsOk){
					IsOk = AuxDataZoneSetItem(SurfZoneNum, CSMAuxData.GBA.SphereCPNum, to_string(CPNum));
				}
				if (IsOk){
					string MovedPointNumsStr = "";
					for (const auto & i : MovedPointNums){
						MovedPointNumsStr += to_string(i+1) + ","; // +1 to switch from base-0 to base-1 so the node nums match with resulting zone
					}
					IsOk = AuxDataZoneSetItem(SurfZoneNum, CSMAuxData.GBA.SphereConstrainedNodeNums, MovedPointNumsStr);
				}
				if (IsOk){
					IsOk = AuxDataZoneSetItem(SurfZoneNum, CSMAuxData.GBA.SphereConstrainedNodeIntersectCPTypes, MovedPointCPTypes);
				}
				if (IsOk){
					IsOk = AuxDataZoneSetItem(SurfZoneNum, CSMAuxData.GBA.SphereConstrainedNodeIntersectCPNames, MovedPointCPNames);
				}
				AuxDataZoneSetItem(SurfZoneNum, CSMAuxData.GBA.NumGBs, to_string(NumTriangles));
				AuxDataZoneSetItem(SurfZoneNum, CSMAuxData.GBA.GPsPerGB, to_string(3 + 3 * NumEdgeGPs));
				AuxDataZoneSetItem(SurfZoneNum, CSMAuxData.GBA.PointsPerGP, to_string(NumSTPoints));
				AuxDataZoneSetItem(SurfZoneNum, CSMAuxData.GBA.SourceZoneNum, to_string(CPZoneNum));
			}
		}

		/*
		*	Get a list of all unique nodes that neighbor nodes that
		*	were moved to GP intersections
		*/
// 		TecUtilDataLoadBegin();
		vector<vector<int> > ConstrainedNeighborNodesNum(MovedPointNums.size());
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
			NumNeighborNodes += (int)ConstrainedNeighborNodesNum[i].size();
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
// 		TecUtilDataLoadEnd();

		/*
		 *	Need to find the nodes closest to any SGP-sphere intersections for
		 *	non-saddle CPs
		 *	TODO: this is a hack job. stupid to check all the distances for such little information.
		 *	find a better way to do it.
		 */
		for (int i = 0; i < AllIntCPNodeNums.size(); ++i){
			if (AllIntCPNodeNums[i].size() < 3){
				vec3 Pt1, Pt2;
				for (int j = 0; j < 3; ++j){
					Pt1[j] = TecUtilDataValueGetByZoneVar(CPZoneNum, XYZVarNums[j], AllIntCPNodeNums[i][1]);
				}
				double MinDist = 1e200;
				int MinInd = -1;
				for (int j = 0; j < NumPoints; ++j){
					for (int k = 0; k < 3; ++k)
						Pt2[k] = p[j][k];
					double TmpDist = DistSqr(Pt1, Pt2);
					if (TmpDist < MinDist){
						MinDist = TmpDist;
						MinInd = j;
					}
				}
				if (MinInd >= 0)
					AllIntCPNodeNums[i].push_back(MinInd);
			}
		}



		/*
		 *	Get all necessary raw pointers for GradPath creation.
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

		}

		vector<FieldDataPointer_c> GradRawPtrs(3);
		FieldDataPointer_c RhoRawPtr;

// 		Data load begin tells Tecplot that I don't want it to reorganize memory until I call
// 		Data load end, so that the raw pointers fetched below remain valid.
		TecUtilDataLoadBegin();

		for (int i = 0; i < 4 && IsOk; ++i){
			if (i == 0){
				IsOk = RhoRawPtr.GetReadPtr(VolZoneNum, CutoffVarNum);
			}
			else{
				IsOk = GradRawPtrs[i - 1].GetReadPtr(VolZoneNum, GradXYZVarNums[i - 1]);
			}
		}

		EntIndex_t NumZonesBeforeVolumes = TecUtilDataSetGetNumZones();

		MemoryRequired = (sizeof(Boolean_t) * (2 * NumPoints + 2 * NumTriangles) + sizeof(EntIndex_t) * NumTriangles) / 1024;
		TecUtilMemoryChangeNotify(MemoryRequired);

		vector<Boolean_t> NodeHasSaddleCP(NumPoints, FALSE),
			NodeIsConstrained(NumPoints, FALSE),
			TriangleHasSaddleCP(NumTriangles, FALSE),
			TriangleHasConstrainedNode(NumTriangles, FALSE);
		for (int i = 0; i < IntZoneSaddleCPNodeNums.size(); ++i){
			NodeHasSaddleCP[IntZoneSaddleCPNodeNums[i][3]] = TRUE;
		}
		for (int i = 0; i < MovedPointNums.size(); ++i){
			NodeIsConstrained[MovedPointNums[i]] = TRUE;
		}

		
#pragma omp parallel for
		for (int i = 0; i < NumTriangles; ++i){
			for (int j = 0; j < 3; ++j){
				if (NodeHasSaddleCP[t[i][j]]){
					TriangleHasSaddleCP[i] = TRUE;
					TriangleHasConstrainedNode[i] = TRUE;
					break;
				}
				else if (NodeIsConstrained[t[i][j]]){
					TriangleHasConstrainedNode[i] = TRUE;
					break;
				}
			}
		}

		TecUtilMemoryChangeNotify((NumSTPoints * 4 * NumPoints * sizeof(double)) / 1024);

		/*
		*	Seed grad paths for all non-saddle-cp nodes
		*	and then for edges.
		*/


		vector<GradPath_c> GPsNonSaddle(NumPoints), GPsEdges(NumEdges * NumEdgeGPs);
		GPTerminate_e HowTerminate;
		if (UseCutoff)
			HowTerminate = GPTerminate_AtRhoValue;
		else
			HowTerminate = GPTerminate_AtBoundary;

		StreamDir_e StreamDir;
		if (CPType == -3)
			StreamDir = StreamDir_Reverse;
		else
			StreamDir = StreamDir_Forward;


		TecUtilDrawGraphics(TRUE);
		TecUtilPleaseWait(NULL, FALSE);
		TecUtilDrawGraphics(FALSE);


		if (IsOk){
			TecUtilDrawGraphics(TRUE);
			TecUtilPleaseWait(NULL, FALSE);
			TmpString = ProgressStr.str() + string(" ... Step 2 of 4 ... Seeding gradient paths");
			TecUtilStatusStartPercentDone(TmpString.c_str(), TRUE, TRUE);
			TecUtilStatusSuspend(TRUE);
			TecUtilDrawGraphics(FALSE);
		}

		Boolean_t UserQuit = FALSE;

		int TmpNumIterations = NumPoints / numCPU;
		int NumCompleted = 0;

#pragma omp parallel for schedule(dynamic)
		for (int PtNum = 0; PtNum < NumPoints; ++PtNum){
			if (omp_get_thread_num() == 0){
				TecUtilDrawGraphics(TRUE);
				UserQuit = !SetPercent(NumCompleted, TmpNumIterations, TmpString, AddOnID);
				NumCompleted++;
#pragma omp flush (UserQuit)
				TecUtilDrawGraphics(FALSE);
			}
#pragma omp flush (UserQuit)

			if (!NodeHasSaddleCP[PtNum] && !UserQuit){
				vec3 NodePos;
				for (int ii = 0; ii < 3; ++ii)
					NodePos[ii] = p[PtNum][ii];

				IsOk = GPsNonSaddle[PtNum].SetupGradPath(NodePos,
					StreamDir, 
					NumSTPoints,
					HowTerminate, 
					NULL, vector<FieldDataPointer_c>(), NULL, NULL,
					&CutoffVal, 
					VolMaxIJK,
					VolMaxXYZ,
					VolMinXYZ,
					GradRawPtrs, 
					RhoRawPtr);

				if (IsOk)
					IsOk = GPsNonSaddle[PtNum].Seed();

				if (IsOk && DistSqr(NodePos, GPsNonSaddle[PtNum].XYZAt(0)) > DistSqr(NodePos, GPsNonSaddle[PtNum].XYZAt(GPsNonSaddle[PtNum].GetCount() - 1))){
					IsOk = GPsNonSaddle[PtNum].Reverse();
				}
			}
		}

		vector<vector<int> > ConstrainedNeighborEdgeNodesNum = ConstrainedNeighborNodesNum;
		vector<vector<int> > ConstrainedNeighborhoodEdgeNums(ConstrainedNeighborNodesNum.size());

#pragma omp parallel for schedule(dynamic)
		for (int EdgeNum = 0; EdgeNum < NumEdges; ++EdgeNum){
			if (omp_get_thread_num() == 0){
				TecUtilDrawGraphics(TRUE);
				UserQuit = !SetPercent(NumCompleted, TmpNumIterations, TmpString, AddOnID);
				NumCompleted++;
#pragma omp flush (UserQuit)
				TecUtilDrawGraphics(FALSE);
			}
#pragma omp flush (UserQuit)
			// Get the first edge node and a vector to step down the edge
			vec3 eNodes[2], DelVec, SeedPt;
			for (int i = 0; i < 2; ++i){
				for (int j = 0; j < 3; ++j){
					eNodes[i][j] = p[e[EdgeNum][i]][j];
				}
			}
			DelVec = (eNodes[1] - eNodes[0]) / static_cast<double>(NumEdgeGPs + 1);
			for (int EdgeGPNum = 0; EdgeGPNum < NumEdgeGPs; ++EdgeGPNum){
				if (!UserQuit){
					SeedPt = eNodes[0] + DelVec * static_cast<double>(EdgeGPNum + 1);
					int GPInd = EdgeNum * NumEdgeGPs + EdgeGPNum;

					IsOk = GPsEdges[GPInd].SetupGradPath(SeedPt,
						StreamDir,
						NumSTPoints,
						HowTerminate,
						NULL, vector<FieldDataPointer_c>(), NULL, NULL,
						&CutoffVal,
						VolMaxIJK,
						VolMaxXYZ,
						VolMinXYZ,
						GradRawPtrs,
						RhoRawPtr);

					if (IsOk)
						IsOk = GPsEdges[GPInd].Seed();

					if (IsOk && DistSqr(SeedPt, GPsEdges[GPInd].XYZAt(0)) > DistSqr(SeedPt, GPsEdges[GPInd].XYZAt(GPsEdges[GPInd].GetCount() - 1))){
						IsOk = GPsEdges[GPInd].Reverse();
					}
				}
			}

			// Also update the ConstrainedNeighborEdgeNodesNum list so that the saddle GPs know their
			// new neighbors.
			for (int m = 0; m < MovedPointNums.size(); ++m){
				int ei;
				for (ei = 0; ei < 2; ++ei){
					if (MovedPointNums[m] == e[EdgeNum][ei]){
						break;
					}
				}
				if (ei < 2){
					for (int & n : ConstrainedNeighborEdgeNodesNum[m]){
						if (e[EdgeNum][(ei + 1) % 2] == n){
							if (ei == 0){
								n = EdgeNum * NumEdgeGPs;
							}
							else{
								n = EdgeNum * NumEdgeGPs + NumEdgeGPs - 1;
							}
							break;
						}
					}
					break;
				}
			}
		}

		if (UserQuit){
			TecUtilDrawGraphics(TRUE);
			TecUtilStatusSuspend(FALSE);
			TecUtilStatusFinishPercentDone();

			MemoryRequired = sizeof(point) * NumPoints
				+ sizeof(triangle) * NumTriangles
				+ sizeof(int) * 2 * NumEdges;
			TecUtilMemoryChangeNotify(-MemoryRequired / 1024);
			delete p, t;
			for (int i = 0; i < NumEdges; ++i){
				delete[] e[i];
			}
			delete e;

			TecUtilMemoryChangeNotify((NumSTPoints * 4 * NumPoints * sizeof(double)) / 1024);
			GPsNonSaddle.clear();

			TecUtilDataLoadEnd();
			TecUtilLockFinish(AddOnID);
			return;
		}


		/*
		 *	Seed grad paths for saddle-cp nodes
		 */
		TecUtilMemoryChangeNotify((NumSTPoints * 4 * NumPoints * sizeof(double)) / 1024);
		TecUtilMemoryChangeNotify((NumSTPoints * 4 * IntZoneSaddleCPNodeNums.size() * sizeof(double)) / 1024);
		vector<vector<GradPath_c> > GPsSaddle(IntZoneSaddleCPNodeNums.size()), GPsSaddleEdges(IntZoneSaddleCPNodeNums.size());
		vector<vector<int> > GPsSaddleClosestPtNum(IntZoneSaddleCPNodeNums.size()), GPsNonSaddleClosestPtNums(IntZoneSaddleCPNodeNums.size());
		
#ifdef _DEBUG
		vector<vector<int> > EdgeList(NumEdges, vector<int>(2));
#pragma omp parallel for
		for (int i = 0; i < NumEdges; ++i){
			for (int j = 0; j < 2; ++j){
				EdgeList[i][j] = e[i][j];
			}
		}
#endif // _DEBUG


#ifndef _DEBUG
#pragma omp parallel for schedule(dynamic)
#endif // !_DEBUG
		for (int i = 0; i < IntZoneSaddleCPNodeNums.size(); ++i){
			/*
			 *	Get the point's index in the list of moved points so
			 *	neighbor information can be looked up.
			 *	Its index in MovedPointNums is the same as its index
			 *	in ConstrainedNeighborNodesNum.
			 */
			Boolean_t IsFound = FALSE;
			for (int j = 0; j < MovedPointNums.size() && !IsFound; ++j){
				if (IntZoneSaddleCPNodeNums[i][3] == MovedPointNums[j]){
					IsFound = TRUE;

					GPsSaddle[i].resize(ConstrainedNeighborNodesNum[j].size());
					GPsSaddleEdges[i].resize(ConstrainedNeighborNodesNum[j].size() * NumEdgeGPs);
					ConstrainedNeighborhoodEdgeNums.reserve(ConstrainedNeighborNodesNum[j].size());
					GPsSaddleClosestPtNum[i].resize(ConstrainedNeighborNodesNum[j].size());
					GPsNonSaddleClosestPtNums[i].resize(ConstrainedNeighborNodesNum[j].size());

					vec3 NodePos;
					for (int ii = 0; ii < 3; ++ii)
						NodePos[ii] = p[MovedPointNums[j]][ii];

					vec3 VolCPPos = IntCPPos[i];

					vector<vec3> SeedPts(ConstrainedNeighborNodesNum[j].size());

					for (int k = 0; k < ConstrainedNeighborNodesNum[j].size() && IsOk; ++k){
						/*
						*	For each neighboring node, need to get an orthogonal set
						*	of unit vectors to describe the positions of the neighbor
						*	node's streamtraces's closest point to the volume CP
						*	relative to the volume CP.
						*/
						/*
						 *	Get the closest point between the saddle CP and
						 *	the neighbor node's grad path.
						 */
// 						vec3 ClosestPoint = GPsNonSaddle[ConstrainedNeighborNodesNum[j][k]].ClosestPoint(VolCPPos, GPsNonSaddleClosestPtNums[i][k]);
						vec3 ClosestPoint = GPsEdges[ConstrainedNeighborEdgeNodesNum[j][k]].ClosestPoint(VolCPPos, GPsNonSaddleClosestPtNums[i][k]);
						/*
						*	v1 is the direction from the main node to the volume CP.
						*	v2 is the direction from main node to neighbor node, and is
						*	orthogonal to v1.
						*	v3 is orthogonal to v1 and v2.
						*/

						vec3 v1 = VolCPPos - NodePos,
							v2 = ClosestPoint - VolCPPos,
							v3 = cross(v1, v2);
						v1 = normalise(v1);
						v3 = normalise(v3);
						v2 = normalise(cross(v3, v1));

						SeedPts[k] = v1 * (Radius * 0.1);
						double SeedTheta = PI / 2.0;
						vec4 TmpVec4 = join_cols(SeedPts[k], ones<vec>(1));
						TmpVec4 = RotationMatrix(SeedTheta, v3) * TmpVec4;

						SeedPts[k] = vec3(TmpVec4.subvec(0, 2)) + VolCPPos;

						IsOk = GPsSaddle[i][k].SetupGradPath(SeedPts[k], StreamDir, NumSTPoints, HowTerminate,
							NULL, vector<FieldDataPointer_c>(), NULL, NULL, &CutoffVal,
							VolMaxIJK, VolMaxXYZ, VolMinXYZ,
							GradRawPtrs,
							RhoRawPtr);

						if (IsOk){
							IsOk = GPsSaddle[i][k].Seed();
// 							GPsSaddle[i][k].SaveAsOrderedZone("Saddle GP " + CPName + " Node " + to_string(ConstrainedNeighborNodesNum[j][k]), Green_C);
						}

						if (IsOk){
							GPsSaddle[i][k].ConcatenateResample(IntGPs[i], NumSTPoints, GPsSaddleClosestPtNum[i][k]);
// 							GPsSaddle[i][k].SaveAsOrderedZone("Full Saddle GP " + CPName + " Node " + to_string(ConstrainedNeighborNodesNum[j][k]), Purple_C);
						}

						if (IsOk && DistSqr(NodePos, GPsSaddle[i][k].XYZAt(0)) > DistSqr(NodePos, GPsSaddle[i][k].XYZAt(GPsSaddle[i][k].GetCount() - 1))){
							IsOk = GPsSaddle[i][k].Reverse();
						}
// 						else{
// 							GPsSaddleClosestPtNum[i][k] = GPsSaddle[i][k].GetCount() - GPsSaddleClosestPtNum[i][k] - 1;
// 						}


					}

					/*
					*	Now make edge GPs between this saddle GP and its clockwise neighboring saddle GP
					*/
					for (int k = 0; k < ConstrainedNeighborNodesNum[j].size() && IsOk; ++k){
						// always make saddle edge GPs in direction of first to second edge node
						// (where the edge is the edge opposite the sphere saddle node)
						// 

						int StartNodeNum, EndNodeNum;
						vector<int> n;
						for (int kk = k + 1; kk < ConstrainedNeighborNodesNum[j].size(); ++kk){
							bool NodeFound = false;
							n = vector<int>({ ConstrainedNeighborNodesNum[j][k], ConstrainedNeighborNodesNum[j][kk] });
							for (int EdgeNum = 0; EdgeNum < NodeEdgeNums[n[0]].size() && !NodeFound; ++EdgeNum){
								if (std::find(ConstrainedNeighborhoodEdgeNums[j].begin(), ConstrainedNeighborhoodEdgeNums[j].end(), NodeEdgeNums[n[0]][EdgeNum]) == ConstrainedNeighborhoodEdgeNums[j].end()){
									for (int ei = 0; ei < 2; ++ei){
										if (e[NodeEdgeNums[n[0]][EdgeNum]][ei] == n[0] && e[NodeEdgeNums[n[0]][EdgeNum]][(ei + 1) % 2] == n[1]) {
											NodeFound = true;
											ConstrainedNeighborhoodEdgeNums[j].push_back(NodeEdgeNums[n[0]][EdgeNum]);
											if (ei == 0){
												n = vector<int>({ k, kk });
											}
											else{
												n = vector<int>({ kk, k });
											}
											break;
										}
									}
								}
							}
							if (NodeFound){
								vec3 RotVec = SeedPts[n[0]] - VolCPPos, v1 = normalise(RotVec), v2 = normalise(SeedPts[n[1]] - VolCPPos);
								vec3 RotAxis = normalise(cross(v1, v2));
								vec4 TmpVec4a = join_cols(RotVec, ones<vec>(1));
								double Alpha = acos(dot(v1, v2));
								double DelAlpha = Alpha / static_cast<double>(NumEdgeGPs + 1);

								for (int ei = 0; ei < NumEdgeGPs; ++ei){
									vec4 TmpVec4b = RotationMatrix(DelAlpha * static_cast<double>(ei + 1), RotAxis) * TmpVec4a;
									vec3 SeedPt = vec3(TmpVec4b.subvec(0, 2)) + VolCPPos;

									int GPInd = (ConstrainedNeighborhoodEdgeNums[j].size() - 1) * NumEdgeGPs + ei;

									IsOk = GPsSaddleEdges[i][GPInd].SetupGradPath(SeedPt, StreamDir, NumSTPoints, HowTerminate,
										NULL, vector<FieldDataPointer_c>(), NULL, NULL, &CutoffVal,
										VolMaxIJK, VolMaxXYZ, VolMinXYZ,
										GradRawPtrs,
										RhoRawPtr);

									if (IsOk){
										IsOk = GPsSaddleEdges[i][GPInd].Seed();
										// 							GPsSaddle[i][k].SaveAsOrderedZone("Saddle GP " + CPName + " Node " + to_string(ConstrainedNeighborNodesNum[j][k]), Green_C);
									}

									if (IsOk){
										GPsSaddleEdges[i][GPInd].ConcatenateResample(IntGPs[i], NumSTPoints, GPsSaddleClosestPtNum[i][k]);
										// 							GPsSaddle[i][k].SaveAsOrderedZone("Full Saddle GP " + CPName + " Node " + to_string(ConstrainedNeighborNodesNum[j][k]), Purple_C);
									}

									if (IsOk && DistSqr(NodePos, GPsSaddleEdges[i][GPInd].XYZAt(0)) > DistSqr(NodePos, GPsSaddleEdges[i][GPInd].XYZAt(GPsSaddleEdges[i][GPInd].GetCount() - 1))){
										IsOk = GPsSaddleEdges[i][GPInd].Reverse();
									}
								}
							}
						}
					}


					break;
				}
			}

			if (!IsFound){
				TecUtilDialogErrMsg("Couldn't find point in list");
			}
		}
		TecUtilDataLoadEnd();

// #pragma omp parallel for
// 		for (int i = 0; i < GPsNonSaddle.size(); ++i){
// 			if (GPsNonSaddle[i].IsMade()){
// 				GPsNonSaddle[i].PointPrepend(CPPos, 0.0);
// 			}
// 		}
// 
// #pragma omp parallel for
// 		for (int i = 0; i < GPsSaddle.size(); ++i){
// 			for (int j = 0; j < GPsSaddle[i].size(); ++j){
// 				if (GPsSaddle[i][j].IsMade()){
// 					GPsSaddle[i][j].PointPrepend(CPPos, 0.0);
// 				}
// 			}
// 		}
// 
// 		for (int i = 0; i < GPsNonSaddle.size(); ++i){
// 			if (GPsNonSaddle[i].IsMade()){
// 				GPsNonSaddle[i].SaveAsOrderedZone("GP " + CPName + " " + "Node " + to_string(i + 1));
// 				AuxDataZoneSetItem(GPsNonSaddle[i].GetZoneNum(), GBASphereCPName, CPName);
// 				AuxDataZoneSetItem(GPsNonSaddle[i].GetZoneNum(), GBASphereCPNum, to_string(CPNum));
// 				AuxDataZoneSetItem(GPsNonSaddle[i].GetZoneNum(), GBAGPNodeNum, to_string(i + 1));
// 				AuxDataZoneSetItem(GPsNonSaddle[i].GetZoneNum(), GBAZoneType, GBAZoneTypeGradPath);
// 			}
// 		}
// 
// 		for (int i = 0; i < GPsSaddle.size(); ++i){
// 			for (int j = 0; j < GPsSaddle[i].size(); ++j){
// 				if (GPsSaddle[i][j].IsMade()){
// 					GPsSaddle[i][j].SaveAsOrderedZone("GP " + CPName + " " + "Node " + to_string(MovedPointNums[i] + 1) + "." + to_string(j + 1));
// 					AuxDataZoneSetItem(GPsSaddle[i][j].GetZoneNum(), GBASphereCPName, CPName);
// 					AuxDataZoneSetItem(GPsSaddle[i][j].GetZoneNum(), GBASphereCPNum, to_string(CPNum));
// 					AuxDataZoneSetItem(GPsSaddle[i][j].GetZoneNum(), GBAGPNodeNum, to_string(MovedPointNums[i] + 1));
// 					AuxDataZoneSetItem(GPsSaddle[i][j].GetZoneNum(), GBAZoneType, GBAZoneTypeGradPath);
// 				}
// 			}
// 		}


		/*
		 *	Now have all gradient paths.
		 *	The concatenated gradient paths are indexed identically
		 *	to the SaddleNeighborNodesNum list.
		 *	
		 *	Now find the sphere edges whose grad paths straddle
		 *	an IB.
		 */


// 		TecUtilMemoryChangeNotify((NumEdges * (sizeof(Boolean_t) + sizeof(vec3))) / 1024);
// 
// 		vector<Boolean_t> EdgeStraddlesIB(NumEdges, FALSE);
// 		vector<vec3> EdgeMidpoints(NumEdges);
// 
// 
// 		TecUtilDataLoadBegin();
// #pragma omp parallel for schedule(dynamic)
// 		for (int EdgeNum = 0; EdgeNum < NumEdges; ++EdgeNum){
// 			if (NodeIsConstrained[e[EdgeNum][0]] || NodeIsConstrained[e[EdgeNum][1]]){
// 				/*
// 				 *	WTF right?
// 				 *	Well, if one of the edge's nodes is a saddle-point node then
// 				 *	this finds the GP that was made for that node to connect to 
// 				 *	the edge's other node.
// 				 */
// 				for (int i = 0; i < 2; ++i){
// 					if (NodeIsConstrained[e[EdgeNum][i]]){
// 						for (int j = 0; j < MovedPointNums.size(); ++j){
// 							if (e[EdgeNum][i] == MovedPointNums[j]){
// 								for (int k = 0; k < ConstrainedNeighborNodesNum[j].size(); ++k){
// 									if (e[EdgeNum][(i + 1) % 2] == ConstrainedNeighborNodesNum[j][k]){
// 										for (int ii = 0; ii < IntZoneSaddleCPNodeNums.size(); ++ii){
// 											if (IntZoneSaddleCPNodeNums[ii][3] == MovedPointNums[j]){
// 												EdgeStraddlesIB[EdgeNum] = GPsStraddleIB(GPsNonSaddle[e[EdgeNum][(i + 1) % 2]], GPsSaddle[ii][k], IBCheckAngle, IBCheckDistRatio);
// 											}
// 										}
// 									}
// 								}
// 							}
// 						}
// 						break;
// 					}
// 				}
// 			}
// 			else
// 			{
// 				EdgeStraddlesIB[EdgeNum] = GPsStraddleIB(GPsNonSaddle[e[EdgeNum][0]], GPsNonSaddle[e[EdgeNum][1]], IBCheckAngle, IBCheckDistRatio);
// 			}
// 
// 			if (EdgeStraddlesIB[EdgeNum]){
// 				vec3 Pts[2];
// 
// 				for (int i = 0; i < 2; ++i)
// 					for (int j = 0; j < 3; ++j)
// 						Pts[i][j] = p[e[EdgeNum][i]][j];
// 
// 				EdgeMidpoints[EdgeNum] = (Pts[0] + Pts[1]) * 0.5;
// 			}
// 		}
// 
// 		int IBEdgeCount = 0;
// 		for (int i = 0; i < NumEdges; ++i)
// 			IBEdgeCount += (int)EdgeStraddlesIB[i];
// 
// 		if (IsOk && IBEdgeCount > 0){
// 			IsOk = TecUtilDataSetAddZone("Inter-IB Edges", IBEdgeCount, 1, 1, ZoneType_Ordered, NULL);
// 			EntIndex_t IBZoneNum;
// 			FieldData_pa IBXYZFDPtrs[3];
// 			if (IsOk){
// 				IBZoneNum = TecUtilDataSetGetNumZones();
// 				for (int i = 0; i < 3 && IsOk; ++i){
// 					IBXYZFDPtrs[i] = TecUtilDataValueGetWritableRef(IBZoneNum, XYZVarNums[i]);
// 					IsOk = VALID_REF(IBXYZFDPtrs[i]);
// 				}
// 			}
// 			if (IsOk){
// 				IBEdgeCount = 1;
// 				for (int i = 0; i < NumEdges; ++i){
// 					if (EdgeStraddlesIB[i]){
// 
// 						for (int j = 0; j < 3; ++j){
// 							TecUtilDataValueSetByRef(IBXYZFDPtrs[j], IBEdgeCount, EdgeMidpoints[i][j]);
// 						}
// 						IBEdgeCount++;
// 					}
// 				}
// 				Set_pa IBSet = TecUtilSetAlloc(FALSE);
// 				TecUtilSetAddMember(IBSet, IBZoneNum, FALSE);
// 
// 				TecUtilZoneSetScatter(SV_SHOW, IBSet, 0.0, TRUE);
// 				TecUtilZoneSetContour(SV_SHOW, IBSet, 0.0, FALSE);
// 				TecUtilZoneSetMesh(SV_SHOW, IBSet, 0.0, FALSE);
// 
// 				TecUtilZoneSetScatter(SV_COLOR, IBSet, NULL, Black_C);
// 				TecUtilZoneSetScatter(SV_FRAMESIZE, IBSet, 0.5, NULL);
// 				TecUtilZoneSetScatterSymbolShape(SV_GEOMSHAPE, IBSet, GeomShape_Sphere);
// 
// 				TecUtilZoneSetActive(IBSet, AssignOp_PlusEquals);
// 
// 				TecUtilSetDealloc(&IBSet);
// 			}
// 			if (IsOk){
// 				IsOk = AuxDataZoneSetItem(IBZoneNum, GBASphereCPName, CPName);
// 				if (IsOk)
// 					IsOk = AuxDataZoneSetItem(IBZoneNum, GBAZoneType, GBAZoneTypeIBEdgeZone);
// 			}
// 		}
// 		TecUtilDataLoadEnd();
// 
// 
// 		EdgeMidpoints.clear();
// 		EdgeStraddlesIB.clear();
// 		TecUtilMemoryChangeNotify((NumEdges * (sizeof(Boolean_t) + sizeof(vec3))) / 1024);

		/*
		 *	Now put together the grad paths into fe volume
		 *	objects.
		 */

// 		int NumSaddleFEZones = 0;
// 		for (int i = 0; i < GPsSaddle.size(); ++i)
// 			NumSaddleFEZones += static_cast<int>(GPsSaddle[i].size());

		TecUtilMemoryChangeNotify((NumSTPoints * 4 * 4 * NumTriangles * sizeof(double)) / 1024);
		vector<FEVolume_c> FEVolumes(NumTriangles);


		TmpString = ProgressStr.str() + string(" ... Step 3 of 4 ... Making volumes");
		TmpNumIterations = NumTriangles / numCPU;
		NumCompleted = 0;

		vector<vector<int> > GPSaddleNodeNums1(NumTriangles), GPSaddleNodeNums2(NumTriangles), GPNonSaddleNodeNums(NumTriangles);
		for (int i = 0; i < NumTriangles; ++i){
			GPSaddleNodeNums1[i].reserve(3);
			GPSaddleNodeNums2[i].reserve(3);
			GPNonSaddleNodeNums[i].reserve(3);
		}


		TecUtilDataLoadBegin();
#ifndef _DEBUG
#pragma omp parallel for schedule(dynamic)
#endif
		for (int TriNum = 0; TriNum < NumTriangles; ++TriNum){
			if (omp_get_thread_num() == 0){
				TecUtilDrawGraphics(TRUE);
				UserQuit = !SetPercent(NumCompleted, TmpNumIterations, TmpString, AddOnID);
				NumCompleted++;
#pragma omp flush (UserQuit)
				TecUtilDrawGraphics(FALSE);
			}
#pragma omp flush (UserQuit)
			if (!UserQuit){
				vector<GradPath_c*> GPPtrs;
				if (TriangleHasSaddleCP[TriNum]){
					GPPtrs.reserve(3 + 3 * NumEdgeGPs + 1);
					/*
					 *	Need to find which GPs to use
					 *	for the FE volume.
					 */
					int Count = 0;
					int NonConstrainedNodeNums[2];
					int ConstrainedGPNums[2] = { -1, -1 };
					int ConstrainedNodeNum;
					int SaddleNum;
					int ConstrainedCornerNum;
					for (int i = 0; i < 3; ++i){
						if (NodeIsConstrained[t[TriNum][i]]){
							ConstrainedNodeNum = t[TriNum][i];
							ConstrainedCornerNum = i;
							for (int j = 0; j < 2; ++j)
								NonConstrainedNodeNums[j] = t[TriNum][(i + j + 1) % 3];
						}
					}

					Count = 0;
					Boolean_t IsFound = FALSE;

					for (int i = 0; i < IntZoneSaddleCPNodeNums.size(); ++i){
						if (IntZoneSaddleCPNodeNums[i][3] == ConstrainedNodeNum){
							SaddleNum = i;
							for (int j = 0; j < MovedPointNums.size(); ++j){
								if (IntZoneSaddleCPNodeNums[i][3] == MovedPointNums[j]){
									for (int k = 0; k < ConstrainedNeighborNodesNum[j].size(); ++k){
										for (int ii = 0; ii < 2; ++ii){
											if (ConstrainedNeighborNodesNum[j][k] == NonConstrainedNodeNums[ii]){
												ConstrainedGPNums[ii] = k;
												Count++;
												break;
											}
										}
										if (Count >= 2) IsFound = TRUE;
									}
									break;
								}
							}
							break;
						}
					}
					if (IsFound){
						vector<int> Nodes({
							NonConstrainedNodeNums[0],
							NonConstrainedNodeNums[1],
							ConstrainedNodeNum
						});
						for (int n = 0; n < Nodes.size(); ++n){
							if (n < Nodes.size() - 1){
								GPPtrs.push_back(&GPsNonSaddle[Nodes[n]]);
							}
							else{
								GPPtrs.push_back(&GPsSaddle[SaddleNum][ConstrainedGPNums[1]]);
								bool EdgeFound = false;
								for (int ei = 0; ei < ConstrainedNeighborhoodEdgeNums[SaddleNum].size() && !EdgeFound; ++ei){
									if (std::find(TriangleEdgeNums[TriNum].begin(), TriangleEdgeNums[TriNum].end(), ConstrainedNeighborhoodEdgeNums[SaddleNum][ei]) != TriangleEdgeNums[TriNum].end()){
										if (e[ConstrainedNeighborhoodEdgeNums[SaddleNum][ei]][0] == NonConstrainedNodeNums[1]){
											for (int i = 0; i < NumEdgeGPs; ++i){
												GPPtrs.push_back(&GPsSaddleEdges[SaddleNum][ei * NumEdgeGPs + i]);
											}
										}
										else{
											for (int i = NumEdgeGPs - 1; i >= 0; --i){
												GPPtrs.push_back(&GPsSaddleEdges[SaddleNum][ei * NumEdgeGPs + i]);
											}
										}
									}
								}
								GPPtrs.push_back(&GPsSaddle[SaddleNum][ConstrainedGPNums[0]]);
							}
							bool EdgeFound = false;
							for (int ei = 0; ei < 3 && !EdgeFound; ++ei){
								for (int i = 0; i < 2; ++i){
									if (e[TriangleEdgeNums[TriNum][ei]][i] == Nodes[n] && e[TriangleEdgeNums[TriNum][ei]][(i + 1) % 2] == Nodes[(n + 1) % Nodes.size()]) {
										if (i == 0){
											for (int j = 0; j < NumEdgeGPs; ++j){
												GPPtrs.push_back(&GPsEdges[TriangleEdgeNums[TriNum][ei] * NumEdgeGPs + j]);
											}
										}
										else{
											for (int j = NumEdgeGPs - 1; j >= 0; --j){
												GPPtrs.push_back(&GPsEdges[TriangleEdgeNums[TriNum][ei] * NumEdgeGPs + j]);
											}
										}
										EdgeFound = true;
										break;
									}
								}
							}
						}
// 						GPPtrs = {&GPsNonSaddle[NonConstrainedNodeNums[0]],
// 							&GPsNonSaddle[NonConstrainedNodeNums[1]],
// 							&GPsSaddle[ConstrainedNodeNum][ConstrainedGPNums[1]],
// 							&GPsSaddle[ConstrainedNodeNum][ConstrainedGPNums[0]]};
							
						/*IsOk = FEVolumes[TriNum].Make(GPsNonSaddle[NonConstrainedNodeNums[0]],
							GPsNonSaddle[NonConstrainedNodeNums[1]],
							GPsSaddle[ConstrainedNodeNum][ConstrainedGPNums[1]],
							GPsSaddle[ConstrainedNodeNum][ConstrainedGPNums[0]]);*/
						if (IsOk){
							for (int i = 0; i < 2; ++i){
								GPNonSaddleNodeNums[TriNum].push_back(NonConstrainedNodeNums[i]);
								GPSaddleNodeNums1[TriNum].push_back(ConstrainedNodeNum);
								GPSaddleNodeNums2[TriNum].push_back(ConstrainedGPNums[i]);
							}
						}
					}
					else{
						TecUtilDialogErrMsg("Couldn't find GPs to make FE volume");
					}
				}
				else{
					GPPtrs.reserve(3 + 3 * NumEdgeGPs);
					for (int n = 0; n < 3; ++n){
						GPPtrs.push_back(&GPsNonSaddle[t[TriNum][n]]);
						bool EdgeFound = false;
						for (int ei = 0; ei < 3 && !EdgeFound; ++ei){
							for (int i = 0; i < 2; ++i){
								if (e[TriangleEdgeNums[TriNum][ei]][i] == t[TriNum][n] && e[TriangleEdgeNums[TriNum][ei]][(i + 1) % 2] == t[TriNum][(n + 1) % 3]) {
									if (i == 0){
										for (int j = 0; j < NumEdgeGPs; ++j){
											GPPtrs.push_back(&GPsEdges[TriangleEdgeNums[TriNum][ei] * NumEdgeGPs + j]);
										}
									}
									else{
										for (int j = NumEdgeGPs - 1; j >= 0; --j){
											GPPtrs.push_back(&GPsEdges[TriangleEdgeNums[TriNum][ei] * NumEdgeGPs + j]);
										}
									}
									EdgeFound = true;
									break;
								}
							}
						}
					}
// 					GPPtrs = {&GPsNonSaddle[t[TriNum][0]],
// 						&GPsNonSaddle[t[TriNum][1]],
// 						&GPsNonSaddle[t[TriNum][2]]};
					/*IsOk = FEVolumes[TriNum].Make(GPsNonSaddle[t[TriNum][0]],
						GPsNonSaddle[t[TriNum][1]],
						GPsNonSaddle[t[TriNum][2]]);*/
					if (IsOk){
						for (int i = 0; i < 3; ++i){
							GPNonSaddleNodeNums[TriNum].push_back(t[TriNum][i]);
						}
					}
				}
				if (IsOk){
					IsOk = FEVolumes[TriNum].Make(GPPtrs);
				}
			}
		}

		/*
		 *	Save neighboring element numbers to GP zones
		 */
// 		for (int TriNum = 0; TriNum < NumTriangles; ++TriNum){
// 			if (FEVolumes[TriNum].GetNumSides() == 4){
// 				vector<GradPathBase_c> TmpGPs;
// 				for (int i = 0; i < 2; ++i){
// 					/*
// 					*	Temp code for writing out GPs as files
// 					*	for Charles so that he can test his integration code.
// 					*	TODO: get rid of this code!
// 					*/
// 					int SBegPt[] = { 0, GPsSaddleClosestPtNum[GPSaddleNodeNums1[TriNum][i]][GPSaddleNodeNums2[TriNum][i]] };
// 					int SEndPt[] = { GPsSaddleClosestPtNum[GPSaddleNodeNums1[TriNum][i]][GPSaddleNodeNums2[TriNum][i]], -1 };
// 					int NBegPt[] = { 0, GPsNonSaddleClosestPtNums[GPSaddleNodeNums1[TriNum][i]][GPSaddleNodeNums2[TriNum][i]] };
// 					int NEndPt[] = { GPsNonSaddleClosestPtNums[GPSaddleNodeNums1[TriNum][i]][GPSaddleNodeNums2[TriNum][i]], -1 };
// 
// 
// 					for (int j = 0; j < 2; ++j){
// 						GradPathBase_c Tmp;
// 						string b = "_";
// 						if (j == 1)
// 							b += "B_";
// 						Tmp = GPsNonSaddle[GPNonSaddleNodeNums[TriNum][i]].SubGP(NBegPt[j], NEndPt[j]);
// 						Tmp.SaveAsCSV(string("Elem_" + to_string(TriNum) + b + to_string(i) + ".csv"));
// 						Tmp.SaveAsOrderedZone(string("Elem_" + to_string(TriNum) + b + to_string(i)), Red_C);
// 						TmpGPs.push_back(Tmp);
// 
// 						if (i == 0 || j == 1){
// 							Tmp = GPsSaddle[GPSaddleNodeNums1[TriNum][i]][GPSaddleNodeNums2[TriNum][i]].SubGP(SBegPt[j], SEndPt[j]);
// 							Tmp.SaveAsCSV(string("Elem_" + to_string(TriNum) + b + to_string(2 + i) + ".csv"));
// 							Tmp.SaveAsOrderedZone(string("Elem_" + to_string(TriNum) + b + to_string(2 + i)), Blue_C);
// 							TmpGPs.push_back(Tmp);
// 						}
// 					}
// 
// 					AuxDataZoneSetItem(GPsNonSaddle[GPNonSaddleNodeNums[TriNum][i]].GetZoneNum(), GBAGPClosestPtNumToCP, to_string(GPsNonSaddleClosestPtNums[GPSaddleNodeNums1[TriNum][i]][GPSaddleNodeNums2[TriNum][i]]));
// 					AuxDataZoneSetItem(GPsSaddle[GPSaddleNodeNums1[TriNum][i]][GPSaddleNodeNums2[TriNum][i]].GetZoneNum(), GBAGPClosestPtNumToCP, to_string(GPsSaddleClosestPtNum[GPSaddleNodeNums1[TriNum][i]][GPSaddleNodeNums2[TriNum][i]]));
// 				}
// 			}
// 			else{
// 				/*
// 					 *	Temp code for writing out triplets of GPs as triplets of files
// 					 *	for Charles so that he can test his integration code.
// 					 *	TODO: get rid of this code!
// 					 */
// 				for (int i = 0; i < 3; ++i){
// 					GradPath_c Tmp = GPsNonSaddle[t[TriNum][i]];
// 					Tmp.SaveAsCSV(string("Elem_" + to_string(TriNum) + '_' + to_string(i) + ".csv"));
// 					Tmp.SaveAsOrderedZone(string("Elem_" + to_string(TriNum) + '_' + to_string(i)));
// 				}
// 			}
// 			for (int i = 0; i < GPSaddleNodeNums1[TriNum].size(); ++i){
// 				string TmpStr;
// 				if (AuxDataZoneGetItem(GPsSaddle[GPSaddleNodeNums1[TriNum][i]][GPSaddleNodeNums2[TriNum][i]].GetZoneNum(), GBAGPElemNums, TmpStr))
// 					AuxDataZoneSetItem(GPsSaddle[GPSaddleNodeNums1[TriNum][i]][GPSaddleNodeNums2[TriNum][i]].GetZoneNum(), GBAGPElemNums, TmpStr + "," + to_string(TriNum + 1));
// 				else
// 					AuxDataZoneSetItem(GPsSaddle[GPSaddleNodeNums1[TriNum][i]][GPSaddleNodeNums2[TriNum][i]].GetZoneNum(), GBAGPElemNums, to_string(TriNum + 1));
// 			}
// 			for (int i = 0; i < GPNonSaddleNodeNums[TriNum].size(); ++i){
// 				string TmpStr;
// 				if (AuxDataZoneGetItem(GPsNonSaddle[GPNonSaddleNodeNums[TriNum][i]].GetZoneNum(), GBAGPElemNums, TmpStr))
// 					AuxDataZoneSetItem(GPsNonSaddle[GPNonSaddleNodeNums[TriNum][i]].GetZoneNum(), GBAGPElemNums, TmpStr + "," + to_string(TriNum + 1));
// 				else
// 					AuxDataZoneSetItem(GPsNonSaddle[GPNonSaddleNodeNums[TriNum][i]].GetZoneNum(), GBAGPElemNums, to_string(TriNum + 1));
// 			}
// 		}

		if (UserQuit){
			TecUtilDrawGraphics(TRUE);
			TecUtilStatusSuspend(FALSE);
			TecUtilStatusFinishPercentDone();

			MemoryRequired = sizeof(point) * NumPoints
				+ sizeof(triangle) * NumTriangles
				+ sizeof(int) * 2 * NumEdges;
			TecUtilMemoryChangeNotify(-MemoryRequired / 1024);
			delete p, t;
			for (int i = 0; i < NumEdges; ++i){
				delete[] e[i];
			}
			delete e;

			TecUtilMemoryChangeNotify(-int(NumSTPoints * 4 * (NumPoints + IntZoneSaddleCPNodeNums.size()) * sizeof(double)) / 1024);
			GPsNonSaddle.clear();
			GPsSaddle.clear();

			TecUtilMemoryChangeNotify((NumSTPoints * 4 * 4 * NumTriangles * sizeof(double)) / 1024);
			FEVolumes.clear();

			TecUtilDataLoadEnd();
			TecUtilLockFinish(AddOnID);
			return;
		}
		TecUtilDataLoadEnd();


		TecUtilMemoryChangeNotify(-int(NumSTPoints * 4 * (NumPoints + IntZoneSaddleCPNodeNums.size()) * sizeof(double)) / 1024);
		GPsNonSaddle.clear();
		GPsSaddle.clear();

		/*
		 *	Test volume integration of all FE zones
		 */

// 		double IntSum = 0;
// 
// 		for (const auto & i : FEVolumes){
// 			IntSum += i.IntVolume(20, CPPos);
// 		}
// 
// 		TecUtilDialogMessageBox(string("Volume is " + to_string(IntSum)).c_str(), MessageBoxType_Information);

		/*
		 *	Now have structure of all zones in memory,
		 *	so make all FE volume zones.
		 */

		
		TmpString = ProgressStr.str() + string(" ... Step 4 of 4 ... Saving volumes");

		TecUtilDataLoadBegin();
		vector<EntIndex_t> TriangleFEZoneNums(NumTriangles, -1);
		for (int TriNum = 0; TriNum < NumTriangles && IsOk; ++TriNum){
			TecUtilDrawGraphics(TRUE);
			if (!SetPercent(TriNum, NumTriangles, TmpString, AddOnID)){
				TecUtilDrawGraphics(TRUE);
				TecUtilStatusSuspend(FALSE);
				TecUtilStatusFinishPercentDone();

				MemoryRequired = sizeof(point) * NumPoints
					+ sizeof(triangle) * NumTriangles
					+ sizeof(int) * 2 * NumEdges;
				TecUtilMemoryChangeNotify(-MemoryRequired / 1024);
				delete p, t;
				for (int i = 0; i < NumEdges; ++i){
					delete[] e[i];
				}
				delete e;

				TecUtilMemoryChangeNotify((NumSTPoints * 4 * 4 * NumTriangles * sizeof(double)) / 1024);
				FEVolumes.clear();

				TecUtilDataLoadEnd();
				TecUtilLockFinish(AddOnID);
				return;
			}
			TecUtilDrawGraphics(FALSE);

			//TriangleFEZoneNums[TriNum] = FEVolumes[TriNum].SaveAsFEZone(VarDataTypes, XYZVarNums, CutoffVarNum);
			TriangleFEZoneNums[TriNum] = FEVolumes[TriNum].SaveAsTriFEZone("Gradient Bundle " + to_string(TriNum + 1) + " (Zone " + to_string(TecUtilDataSetGetNumZones() + 1) + ")", VarDataTypes, VarLocations, XYZVarNums);

			if (TriangleFEZoneNums[TriNum] > 0){
				/*
				*	Save node and triangle numbers to the fe volume zone's aux data
				*/

				if (IsOk)
					IsOk = AuxDataZoneSetItem(TriangleFEZoneNums[TriNum], CSMAuxData.GBA.SphereCPName, CPString);
				if (IsOk)
					IsOk = AuxDataZoneSetItem(TriangleFEZoneNums[TriNum], CSMAuxData.GBA.ZoneType, CSMAuxData.GBA.ZoneTypeFEVolumeZone);
				if (IsOk)
					IsOk = AuxDataZoneSetItem(TriangleFEZoneNums[TriNum], CSMAuxData.GBA.ElemNum, to_string(TriNum + 1));

				for (int i = 0; i < 3 && IsOk; ++i)
					IsOk = AuxDataZoneSetItem(TriangleFEZoneNums[TriNum], CSMAuxData.GBA.NodeNums[i], to_string(t[TriNum][i] + 1));

				if (IsOk && TriangleHasConstrainedNode[TriNum]){
					Boolean_t IsFound = FALSE;
					for (int i = 0; i < 3 && !IsFound; ++i){
						for (int j = 0; j < AllIntCPNodeNums.size() && !IsFound; ++j){
// 							if (NodeIsConstrained[t[TriNum][i]]){
							if (t[TriNum][i] == AllIntCPNodeNums[j][2]){
								IsFound = TRUE;
								IsOk = AuxDataZoneSetItem(TriangleFEZoneNums[TriNum], CSMAuxData.GBA.VolumeCPName, AllIntVolumeCPNames[j]);
							}
						}
					}
				}
			}
		}
		TecUtilDataLoadEnd();


		TecUtilMemoryChangeNotify((NumSTPoints * 4 * 4 * NumTriangles * sizeof(double)) / 1024);
		FEVolumes.clear();

		MemoryRequired = sizeof(point) * NumPoints
			+ sizeof(triangle) * NumTriangles
			+ sizeof(int) * 2 * NumEdges;
		TecUtilMemoryChangeNotify(-MemoryRequired / 1024);
		delete p, t;
		for (int i = 0; i < NumEdges; ++i){
			delete[] e[i];
		}
		delete e;


		TmpString = ProgressStr.str() + string(" ... Unloading Data");
		EntIndex_t NewNumZones = TecUtilDataSetGetNumZones() - NumZones;
		for (EntIndex_t i = NumZones + 1; i <= NewNumZones && IsOk; ++i){
			TecUtilDrawGraphics(TRUE);
			SetPercent(i, NewNumZones, TmpString, AddOnID);
			TecUtilDrawGraphics(FALSE);
			for (EntIndex_t j = 0; j < 3 && IsOk; ++j){
				IsOk = TecUtilDataValueUnload(i, XYZVarNums[j]);
			}
		}
		TecUtilDrawGraphics(TRUE);
		TecUtilStatusSuspend(FALSE);
		TecUtilStatusFinishPercentDone();
		TecUtilDrawGraphics(FALSE);

	}

	if (DoTiming){
		Time2 = high_resolution_clock::now();
		duration<double> TimeSpan = duration_cast<duration<double>>(Time2 - Time1);
// 		TecUtilDialogMessageBox(string(string("Runtime: ") + to_string(TimeSpan.count()) + string(" seconds.")).c_str(), MessageBoxType_Information);
	}

	Set_pa VolSet = TecUtilSetAlloc(FALSE);
	TecUtilSetAddMember(VolSet, VolZoneNum, FALSE);
	TecUtilZoneSetActive(VolSet, AssignOp_MinusEquals);
	TecUtilSetDealloc(&VolSet);

	TecUtilArrayDealloc((void **)&CPNums);

	TecGUITabSetCurrentPage(TAB1_TB_D1, 3);

	TecUtilDrawGraphics(TRUE);
	TecUtilStatusSuspend(FALSE);
	TecUtilLockFinish(AddOnID);
	return;
}



void GradPathTest(){

	TecUtilLockStart(AddOnID);

	SYSTEM_INFO sysinfo;
	GetSystemInfo(&sysinfo);

	int numCPU = sysinfo.dwNumberOfProcessors;
	omp_set_num_threads(numCPU);

	/*
	 *	Getting all the pointers and system info needed to construct the grad path object
	 */
	EntIndex_t VolZoneNum = ZoneNumByName(string("Full Volume"));

	vector<int> XYZVarNums(3);
	TecUtilAxisGetVarAssignments(&XYZVarNums[0], &XYZVarNums[1], &XYZVarNums[2]);


	EntIndex_t RhoVarNum = VarNumByName(string("Electron Density"));
	vector<int> IntVarList = { RhoVarNum };

	vector<int> ZoneNumList;
	int NumZones = TecUtilDataSetGetNumZones();
	ZoneNumList.reserve(NumZones);

	TecUtilDataLoadBegin();

	/*
	 *	Compute the total rho in the whole system.
	 *	Simply add rho from every point.
	 */
// 	double SysTotal = 0;
// 	int IJK[3];
// 	TecUtilZoneGetIJK(VolZoneNum, &IJK[0], &IJK[1], &IJK[2]);
// 	for (int i = 1; i <= IJK[0]; ++i)
// 		for (int j = 1; j <= IJK[1]; ++j)
// 			for (int k = 1; k <= IJK[2]; ++k)
// 				SysTotal += TecUtilDataValueGetByZoneVar(VolZoneNum, RhoVarNum, IndexFromIJK(i, j, k, IJK[0], IJK[1]));
// 	TecUtilDialogMessageBox(to_string(SysTotal).c_str(), MessageBoxType_Information);

	/*
	 *	Integrate rho over all volumes for Atom 1,
	 *	including sphere zone.
	 *	Try at multiple resolutions for subcell sampling
	 *	and save results and times in Out.txt.
	 */
	vector<FEVolume_c> VolumeList;
	VolumeList.reserve(NumZones);
	string TmpString = "Atom 1";
	for (int i = 1; i <= NumZones; ++i)
		if (TecUtilZoneIsFiniteElement(i)
			&& AuxDataZoneItemMatches(i, CSMAuxData.GBA.SphereCPName, TmpString))
		{
			VolumeList.push_back(FEVolume_c(i, VolZoneNum, XYZVarNums, IntVarList));
		}


	TecUtilStatusStartPercentDone("Integrating", TRUE, TRUE);
	TecUtilStatusSuspend(TRUE);

	vector<int> ResList = { 1, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24 };
	vector<double> TotalList;
	vector<double> TimeList;

	std::ofstream OutFile("Out.csv", std::ios::out | std::ios::app);
	if (!OutFile.is_open())
		TecUtilDialogErrMsg("Failed to open output file");


	OutFile << "Resolution,Time,Integral\n";

	for (int ResNum = 0; ResNum < 3; ResNum++){
	//for (int ResNum = 0; ResNum < ResList.size(); ResNum++){
		high_resolution_clock::time_point Time1, Time2;
		Time1 = high_resolution_clock::now();

		stringstream ProgressBase;
		ProgressBase << "Integrating at res " << ResList[ResNum] << " (" << ResNum + 1 << " of " << ResList.size() << "), zone ";

		double Total = 0.0;
#pragma omp parallel for
		for (int i = 0; i < VolumeList.size(); ++i){
			if (omp_get_thread_num() == 0)
				SetPercent(i * numCPU, static_cast<int>(VolumeList.size()), string(ProgressBase.str() + to_string(i * numCPU) + string(" of ") + to_string(VolumeList.size())), AddOnID);
			//if (VolumeList[i].GetZoneNum() == 132)
// 			VolumeList[i].DoIntegration(ResList[ResNum], FALSE);
		}

		for (int i = 0; i < VolumeList.size(); ++i)
			Total += VolumeList[i].GetIntResults()[0];

		Time2 = high_resolution_clock::now();
		duration<double> TimeSpan = duration_cast<duration<double>>(Time2 - Time1);
		TotalList.push_back(Total);
		TimeList.push_back(TimeSpan.count());

		OutFile << ResList[ResNum] << "," << TimeList[ResNum] << "," << TotalList[ResNum] << std::endl;
		//TecUtilDialogMessageBox(to_string(Total * Ang3PerBohr3).c_str(), MessageBoxType_Information);
	}

	TecUtilStatusSuspend(FALSE);
	TecUtilStatusFinishPercentDone();

	OutFile.close();

	TecUtilDataLoadEnd();

	return;

	Set_pa VolSet = TecUtilSetAlloc(FALSE);
	TecUtilSetAddMember(VolSet, VolZoneNum, FALSE);
	TecUtilZoneSetActive(VolSet, AssignOp_PlusEquals);


	//EntIndex_t XYZVarNums[3];
	TecUtilAxisGetVarAssignments(&XYZVarNums[0], &XYZVarNums[1], &XYZVarNums[2]);

	vector<int> MaxIJK(3);
	TecUtilZoneGetIJK(VolZoneNum, &MaxIJK[0], &MaxIJK[1], &MaxIJK[2]);

	vec3 MaxXYZ, MinXYZ;
	for (int i = 0; i < 3; ++i)
		TecUtilVarGetMinMax(XYZVarNums[i], &MinXYZ[i], &MaxXYZ[i]);

	string GradVarNames[3] = { "X Den", "Y Den", "Z Den" };
	EntIndex_t GradVarNums[3];
	for (int i = 0; i < 3; ++i)
		GradVarNums[i] = VarNumByName(GradVarNames[i]);

	vector<FieldDataPointer_c> GradPtrs(3);
	FieldDataPointer_c RhoPtr;

	for (int i = 0; i < 3; ++i){
		GradPtrs[i].GetReadPtr(VolZoneNum, GradVarNums[i]);
	}
	RhoPtr.GetReadPtr(VolZoneNum, RhoVarNum);

	//Vec3_c StartPoint(-6.4627935942777901, 1.5801340272202937, -0.40415061637454408);
	vec3 StartPoint;
	StartPoint << -6.77934 << -2.2642 << 0.46312;

	double TermRhoValue = 1e-3;

	GradPath_c GP(StartPoint,
		StreamDir_Reverse, 100,
		GPTerminate_AtBoundary,
		NULL, vector<FieldDataPointer_c>(), NULL, NULL, &TermRhoValue,
		MaxIJK, MaxXYZ, MinXYZ,
		GradPtrs,
		RhoPtr);

	TecUtilDataLoadBegin();

	if (GP.IsReady()){
		GP.Seed();
	}
// 	if (GP.IsMade()){
// 		GP.Resample(100);
// 		GP.Reverse();
// 	}

// 	EntIndex_t ConcatZoneNum = 26;
// 
// 	vector<EntIndex_t> VarNumsVec(4);
// 	for (int i = 0; i < 3; ++i)
// 		VarNumsVec[i] = XYZVarNums[i];
// 	VarNumsVec[3] = RhoVarNum;
// 
// 	GradPath_c GP2 = GradPath_c(ConcatZoneNum, VarNumsVec);
// 	GP2.Trim(StartPoint, 1);
// 	GP2.SaveAsOrderedZone();

	//GP += GradPath_c(ConcatZoneNum, VarNumsVec);

	GP.SaveAsOrderedZone("Test GradPath");


	TecUtilDataLoadEnd();

	TecUtilZoneSetActive(VolSet, AssignOp_MinusEquals);
	TecUtilSetDealloc(&VolSet);

	TecUtilDialogMessageBox("Finished", MessageBoxType_Information);

	TecUtilLockFinish(AddOnID);
}