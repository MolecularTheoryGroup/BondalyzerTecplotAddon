
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
#include <set>
#include <map>

// for profiling
#include <chrono>
#include <ctime>
#include <ratio>

#include <omp.h>

#include "Set.h"
#include "StyleValue.h"

#include "CSM_DATA_SET_INFO.h"
#include "CSM_DATA_TYPES.h"
#include "meshgen2d_sphere.h"
#include "VIEWRESULTS.h"
#include "CSM_FE_VOLUME.h"
#include "CSM_GRAD_PATH.h"
#include "CSM_CRIT_POINTS.h"
#include "CSM_GUI.h"
#include "ZONEVARINFO.h"
#include "CSM_GEOMETRY.h"

#include "GBAENGINE.h"

#include "spdlog/spdlog.h"
#include "spdlog/sinks/rotating_file_sink.h" // support for rotating file logging

#include <armadillo>
using namespace arma;
using namespace tecplot::toolbox;

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

const double DivergentGPMaxTerminalDotProduct = cos(170. * PI / 180.);

const string LogName = "GBA_log";
const int LogSizeMB = 5;
const int LogNumLogs = 10;

enum RadMode{
	ABSOLUTERADIUS = 1,
	MINCPDISTRATIO
};

const string GetNucleusNameForCP(const vec3 & CPPt, const vector<int> XYZVarNums = { 1,2,3 }, const double & DistCutoff = 0.1)
{
	REQUIRE(XYZVarNums.size() == 3);
	for (const auto & i : XYZVarNums) REQUIRE(i > 0);
	REQUIRE(DistCutoff > 0);

	double MinDist = DBL_MAX, TmpDist;
	int MinZoneNum = -1, MinPtNum = -1;

	TecUtilDataLoadBegin();
	int NumZones = TecUtilDataSetGetNumZones();
	for (int z = 1; z <= NumZones; ++z) {
		if (AuxDataZoneItemMatches(z, CSMAuxData.DL.ZoneType, CSMAuxData.DL.ZoneTypeNuclearPositions)) {
			FieldVecPointer_c XYZPtr;
			if (XYZPtr.GetReadPtr(z, XYZVarNums)) {
				for (int i = 0; i < XYZPtr.Size(); ++i) {
					TmpDist = DistSqr(XYZPtr[i], CPPt);
					if (TmpDist < MinDist) {
						MinZoneNum = z;
						MinPtNum = i;
						MinDist = TmpDist;
					}
				}
			}
			else {
				TecUtilDialogErrMsg("Failed to get XYZ vec pointer for zone when finding nucleus name");
			}
		}
	}
	TecUtilDataLoadEnd();
	string outStr;
	if (MinZoneNum > 0 && sqrt(MinDist) <= DistCutoff){
		char* cstr;
		if (TecUtilZoneGetName(MinZoneNum, &cstr)){
			outStr = cstr + to_string(MinPtNum + 1);
		}
		else{
			TecUtilDialogErrMsg("Failed to get zone name when finding nucleus name");
		}
		TecUtilStringDealloc(&cstr);
	}
	return outStr;
}

void MainFunction(){

	/*
	 *	Setup log
	 */
	//auto Logger = spdlog::rotating_logger_mt(LogName, "C:\\Users\\Haiiro\\Desktop\\logs\\" + LogName + ".txt", 1024 * 1024 * LogSizeMB, LogNumLogs);


	SYSTEM_INFO sysinfo;
	GetSystemInfo(&sysinfo);

	int numCPU = sysinfo.dwNumberOfProcessors;
	omp_set_num_threads(numCPU);


	TecUtilLockStart(AddOnID);

	Boolean_t DoTiming = TRUE;

	high_resolution_clock::time_point Time1, Time2;
	if (DoTiming){
		Time1 = high_resolution_clock::now();
	}


	Boolean_t IsOk = TRUE;
	int MemoryRequired;
	string TmpString;

	vector<string> GoodSGPEndCPTypes = { "Nuclear", "Atom", "Bond", "Ring", "Cage" };

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


	bool SaveGPs = TecGUIToggleGet(TGLsGP_TOG_T1_1);
	bool SaveGBs = TecGUIToggleGet(TGLsGB_TOG_T1_1);


	VolExtentIndexWeights_s VolInfo;
	GetVolInfo(VolZoneNum, XYZVarNums, FALSE, VolInfo);


	EntIndex_t OldNumZones = TecUtilDataSetGetNumZones();	

	/*
	*	information for integration
	*/
	vector<string> AtomNameList;
	vector<string> IntVarNameList;
	vector<int> IntVarNumList;
	int IntResolution = TecGUIScaleGetValue(SCPrecise_SC_T1_1);
	vector<FieldDataPointer_c> IntVarPtrs;

	Boolean_t DoIntegration = TecGUIToggleGet(TGLInt_TOG_T1_1);

	if (DoIntegration){
		AtomNameList = ListGetSelectedStrings(MLSelCPs_MLST_T1_1);
		IntVarNameList = ListGetSelectedStrings(MLSelVars_MLST_T1_1);
		IntVarNumList = ListGetSelectedItemNums(MLSelVars_MLST_T1_1);

		IntVarPtrs.resize(IntVarNumList.size());
		for (int i = 0; i < IntVarPtrs.size(); ++i){
			IsOk = IntVarPtrs[i].GetReadPtr(VolZoneNum, IntVarNumList[i]);
		}
	}
	

	if (!IsOk){
		TecUtilDrawGraphics(TRUE);
		TecUtilDialogErrMsg("Couldn't get XYZ vars.");
		TecUtilLockFinish(AddOnID);
		return;
	}

	EntIndex_t CPTypeVarNum;
	CPTypeVarNum = VarNumByName("CritPointType");

	Set oldSphereZonesToDelete;
	bool deleteOldSphereZones = false, deleteOldSphereZonesAsked = false;

	vector<int> NumCPs;
	for (int SelectCPNum = 0; SelectCPNum < NumSelectedCPs && IsOk; ++SelectCPNum){

		LgIndex_t CPNum = -1;
		int CPType = -3;

		double Radius = 0.25;

		stringstream ProgressStr;
		string CPString, CPName, CPTypeStr;
		bool IsCP;

		int RealCPNum = 0;
		EntIndex_t CPZoneNum;
		LgIndex_t CPIJKMax[3];

		vec3 CPPos;
		string NucleusName;
		if (IsOk){
			CPPos = GetCoordsFromListItem(CPNums[SelectCPNum], MLSelCPs_MLST_T1_1, &CPString, &CPName, &CPNum, &CPZoneNum, &IsCP, &NumCPs);
			NucleusName = GetNucleusNameForCP(CPPos);
			if (IsCP){

// 				CPZoneNum = ZoneNumByName("Critical Points");

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
					if (CPName == RankStrs[i] || i == 0 && CPName == "Atom"){
						int RealCPNum = 0;
						for (int j = 0; j < i; ++j){
							RealCPNum += NumCPs[j];
						}
						CPNum = RealCPNum + CPNum;
						break;
					}
				}


				CPType = (int)TecUtilDataValueGetByZoneVar(CPZoneNum, CPTypeVarNum, CPNum);
				for (int i = 0; i < CPTypeList.size(); ++i){
					if (CPType == CPTypeList[i]){
						CPTypeStr = CPNameList[i];
					}
				}
			}

			/*
			 *	Delete any existing sphere for the current CP if so chosen by the user
			 */
			for (int z = 1; z <= TecUtilDataSetGetNumZones(); ++z){
				if (AuxDataZoneItemMatches(z, CSMAuxData.GBA.SphereCPName, CPString) && (!deleteOldSphereZonesAsked || deleteOldSphereZones)){
					string tmpStr = "An existing sphere zone was found for " + CPString + ". Would you like to erase it (and all other conflicting sphere zones) before continuing?";
					if (!deleteOldSphereZonesAsked){
						deleteOldSphereZones = TecUtilDialogMessageBox(tmpStr.c_str(), MessageBox_YesNo);
						deleteOldSphereZonesAsked = true;
					}
					if (deleteOldSphereZones) {
						oldSphereZonesToDelete += z;
						for (int zi = z + 1; zi <= TecUtilDataSetGetNumZones(); zi++){
							if (AuxDataZoneItemMatches(zi, CSMAuxData.GBA.SourceZoneNum, to_string(z)))
								oldSphereZonesToDelete += zi;
						}
					}
					break;
				}
			}

			ProgressStr << "Processing " << (NucleusName == "" ? CPString : NucleusName);
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
			StatusLaunch(TmpString, AddOnID, FALSE);
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
		string MovedPointCPTypes, MovedPointCPNames, MovedPointCPNamesTotalCount;
		vector<bool> MovedCPIsBond;

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

		if (RadiusMode == ABSOLUTERADIUS || ClosestCPDist == 1e50 || ClosestCPDist <= 0){
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
				TecUtilZoneSetActive(Set(CurZoneNum).getRef(), AssignOp_PlusEquals);
// 				Set_pa TempSet = TecUtilSetAlloc(FALSE);
// 				TecUtilSetAddMember(TempSet, CurZoneNum, FALSE);
// 				TecUtilZoneSetActive(TempSet, AssignOp_PlusEquals);
// 				TecUtilSetDealloc(&TempSet);
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
				if (!ZoneOK){
					ZoneOK = TRUE;
					for (int i = 0; i < 2 && ZoneOK; ++i){
						ZoneOK = AuxDataZoneGetItem(CurZoneNum, CSMAuxData.CC.GPEndNumStrs[i], TmpString);
						if (ZoneOK){
							StartEndCPs[i] = stoi(TmpString);
							ZoneOK = AuxDataZoneGetItem(CurZoneNum, CSMAuxData.CC.GPEndTypes[i], TmpString);
							if (ZoneOK)
								ZoneOK = (SearchVectorForString(GoodSGPEndCPTypes, TmpString) >= 0);
							if (ZoneOK)
								StartEndTypes[i] = TmpString;
						}
					}
				}

				Boolean_t IntValid = ZoneOK && (CPTypeStr == "" || StartEndTypes[0] == CPTypeStr || StartEndTypes[1] == CPTypeStr) && StartEndTypes[0] != StartEndTypes[1];
				
				if (ZoneOK){
					Boolean_t IntRecorded = FALSE;
					for (int i = 0; i < 2 && IsOk; ++i){
						if (StartEndCPs[i] == CPNum){
							ConnectsToCP = TRUE;
							StartsFrom1 = (i == 0);
// 							for (int j = 2; j < 4 && IsOk && IntValid; ++j){
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
										IntZoneSaddleCPNodeNums.back().push_back(static_cast<int>(IntersectionPoints.size()));
										IntZoneSaddleCPNodeNums.back().push_back(CurZoneNum);
										IntZoneSaddleCPNodeNums.back().push_back(StartEndCPs[(i + 1) % 2]);
										IntCPPos.push_back(TmpCP);
									}

									// 									if (TmpIJK[0] < NumSTPoints){
									// 										NumSTPoints = TmpIJK[0];
									// 										if (NumSTPoints % 2 != 0){
									// 											NumSTPoints--;
									// 										}
									// 									}
									if (IntValid){
										int CPCount = 0;
										for (int k = 0; k < RankStrs.size(); ++k){
											if (StartEndTypes[(i + 1) % 2] == RankStrs[k] && !IntRecorded){
												AllIntVolumeCPNames.push_back(StartEndTypes[(i + 1) % 2] + string(" ") + to_string(StartEndCPs[(i + 1) % 2] - CPCount));
												AllIntCPNodeNums.push_back(vector<LgIndex_t>());
												AllIntCPNodeNums.back().push_back(static_cast<int>(IntersectionPoints.size()));
												AllIntCPNodeNums.back().push_back(StartEndCPs[(i + 1) % 2]);
												Boolean_t ZoneFound = FALSE;
// 												if (StartEndTypes[(i + 1) % 2] == RankStrs[1]){
// 													// if the intersecting SGP connects to a bond point,
// 													// then find the CP on the other side of the bond path.
// 
// 													// The other end of the bond path should be the zone before of after CurZoneNum, so check those first
// 													for (int CheckNum = -1; CheckNum < 2; CheckNum += 2){
// 														if (CurZoneNum + CheckNum <= NumZones && AuxDataZoneGetItem(CurZoneNum + CheckNum, CCDataGPEndNums[0], TmpString)){
// 															// the beginning point for bond paths is always the bond point, so only check CCDataGPEndNums[0]
// 															if (stoi(TmpString) == StartEndCPs[(i + 1) % 2]){
// 																// Matching bond path found, so get nuclear CP on other side
// 																int CheckNucCPNum = stoi(AuxDataZoneGetItem(CurZoneNum + CheckNum, CCDataGPEndNums[1]));
// 																// now search through the nuclear zones, if present, to see what type of element the CP is
// 																if (CheckNucCPNum > 0){
// 																	double MinNucCPCheckDist = DBL_MAX, TempDist;
// 																	string ClosestElementName;
// 																	int ClosestAtomNum = -1;
// 																	for (int CheckZoneNum = 1; CheckZoneNum <= NumZones; ++CheckZoneNum){
// 																		if (CheckZoneNum != CurZoneNum
// 																			&& CheckZoneNum != CurZoneNum + CheckNum
// 																			&& TecUtilZoneIsOrdered(CheckZoneNum)
// 																			&& AuxDataZoneItemMatches(CheckZoneNum, DLZoneType, DLZoneTypeNuclearPositions)){
// 																			// need to now find the 
// 																			int NucIJK[3];
// 																			TecUtilZoneGetIJK(CheckZoneNum, &NucIJK[0], &NucIJK[1], &NucIJK[2]);
// 																			vec3 NucPt, NucCPPt;
// 																			for (int Dir = 0; Dir < 3; ++Dir){
// 																				NucCPPt[Dir] = TecUtilDataValueGetByZoneVar(CurZoneNum + CheckNum, XYZVarNums[Dir], CheckNucCPNum);
// 																			}
// 																			for (int NucCheckPt = 1; NucCheckPt <= NucIJK[0]; ++NucCheckPt){
// 																				for (int Dir = 0; Dir < 3; ++Dir){
// 																					NucPt[Dir] = TecUtilDataValueGetByZoneVar(CheckZoneNum, XYZVarNums[Dir], NucCheckPt);
// 																				}
// 																				TempDist = DistSqr(NucPt, NucCPPt);
// 																				if (TempDist < MinNucCPCheckDist){
// 																					MinNucCPCheckDist = TempDist;
// 																					ClosestElementName = AuxDataZoneGetItem(CheckZoneNum, DLZoneAtomicSpecies);
// 																					ClosestAtomNum = NucCheckPt;
// 																				}
// 																			}
// 																		}
// 																	}
// 																	if (ClosestAtomNum > 0){
// 																		MovedPointCPNames += ClosestElementName + " " + to_string(ClosestAtomNum) + ",";
// 																		ZoneFound = TRUE;
// 																	}
// 																}
// 															}
// 														}
// 													}
// 												}
												MovedPointCPTypes += StartEndTypes[(i + 1) % 2] + ",";
												MovedCPIsBond.push_back(StartEndTypes[(i + 1) % 2] == RankStrs[1]);
												if (!ZoneFound) {
													MovedPointCPNames += StartEndTypes[(i + 1) % 2] + string(" ") + to_string(StartEndCPs[(i + 1) % 2] - CPCount) + ",";
													MovedPointCPNamesTotalCount += StartEndTypes[(i + 1) % 2] + string(" ") + to_string(StartEndCPs[(i + 1) % 2]) + ",";
												}
												IntRecorded = TRUE;
												break;
											}
// 											CPCount += NumCPs[k == 0 ? k : k - 1];
											CPCount += NumCPs[k];
										}
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
						IntPt[i] = IntGPs.back()[-1][i] - CPPos[i];
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

		FESurface_c Sphere;
		EntIndex_t SurfZoneNum;

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

			vector<vec3> P(NumPoints);
			for (int i = 0; i < NumPoints; ++i) for (int j = 0; j < 3; ++j) P[i][j] = p[i][j];
			vector<vector<int> > T(NumTriangles, vector<int>(3));
			for (int i = 0; i < NumTriangles; ++i) for (int j = 0; j < 3; ++j) T[i][j] = t[i][j];

			Sphere.MakeFromNodeElemList(P, T);
			vector<ValueLocation_e> DataLocs(TecUtilDataSetGetNumVars(), ValueLocation_Nodal);
			vector<FieldDataType_e> DataTypes(TecUtilDataSetGetNumVars(), FieldDataType_Double);
			vector<string> IntVarPrefixes = { "I: ","IN: ","INS: " };
			for (int i = 0; i < DataLocs.size(); ++i) {
				char * cstr;
				TecUtilVarGetName(i+1, &cstr);
				string str = cstr;
				TecUtilStringDealloc(&cstr);
				for (const auto & s : IntVarPrefixes) {
					if (str.length() >= s.length()){
						int j = str.compare(0, s.length(), s);
						if (j == 0){
							DataLocs[i] = ValueLocation_CellCentered;
							break;
						}
					}
				}
			}
			SurfZoneNum = Sphere.SaveAsTriFEZone(NucleusName != "" ? NucleusName + " Sphere" : "Sphere (" + CPString + ")", DataTypes, DataLocs, XYZVarNums);
			Set TmpSet(SurfZoneNum);
			TecUtilZoneSetMesh(SV_SHOW, TmpSet.getRef(), 0.0, FALSE);
			TecUtilZoneSetContour(SV_SHOW, TmpSet.getRef(), 0.0, TRUE);
			StyleValue styleValue;
			styleValue.set((Boolean_t)1, SV_FIELDLAYERS, SV_SHOWCONTOUR);

// 			if (IsOk){
// 				ArgList_pa ArgList = TecUtilArgListAlloc();
// 				TecUtilArgListAppendString(ArgList, SV_NAME, string("Sphere (" + CPString + ")").c_str());
// 				TecUtilArgListAppendInt(ArgList, SV_ZONETYPE, ZoneType_FETriangle);
// 				TecUtilArgListAppendInt(ArgList, SV_IMAX, NumPoints);
// 				TecUtilArgListAppendInt(ArgList, SV_JMAX, NumTriangles);
// 				TecUtilArgListAppendArray(ArgList, SV_VALUELOCATION, VarLocations.data());
// 				TecUtilArgListAppendArray(ArgList, SV_VARDATATYPE, VarDataTypes.data());
// 				IsOk = TecUtilDataSetAddZoneX(ArgList);
// 				TecUtilArgListDealloc(&ArgList);
// 			}
// 
// // 			IsOk = TecUtilDataSetAddZone("Sphere", NumPoints, NumTriangles, 0, ZoneType_FETriangle, NULL);
// 			EntIndex_t SurfZoneNum = TecUtilDataSetGetNumZones();
// 			FieldData_pa SurfXYZPtrs[3];
// 
// 			for (int i = 0; i < 3 && IsOk; ++i){
// 				SurfXYZPtrs[i] = TecUtilDataValueGetWritableRef(SurfZoneNum, XYZVarNums[i]);
// 				IsOk = VALID_REF(SurfXYZPtrs[i]);
// 			}
// 
// 			if (IsOk){
// 				NodeMap_pa NodeMap;
// 				NodeMap = TecUtilDataNodeGetWritableRef(SurfZoneNum);
// 				IsOk = VALID_REF(NodeMap);
// 
// 
// 				TecUtilDataLoadBegin();
// 				for (LgIndex_t NodeNum = 0; NodeNum < NumPoints; ++NodeNum){
// 					for (int i = 0; i < 3; ++i){
// 						TecUtilDataValueSetByRef(SurfXYZPtrs[i], NodeNum + 1, p[NodeNum][i]);
// 					}
// 				}
// 				for (LgIndex_t TriNum = 0; TriNum < NumTriangles; ++TriNum){
// 					for (int i = 0; i < 3; ++i){
// 						TecUtilDataNodeSetByRef(NodeMap, TriNum + 1, i + 1, t[TriNum][i] + 1);
// 					}
// 				}
// 				TecUtilDataLoadEnd();
// 
// 			}
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
				AuxDataZoneSetItem(SurfZoneNum, CSMAuxData.GBA.SphereCPName, CPString);
				AuxDataZoneSetItem(SurfZoneNum, CSMAuxData.GBA.ZoneType, CSMAuxData.GBA.ZoneTypeSphereZone);
				AuxDataZoneSetItem(SurfZoneNum, CSMAuxData.GBA.SphereCPNum, to_string(CPNum));
				string MovedPointNumsStr = "";
				for (const auto & i : MovedPointNums){
					MovedPointNumsStr += to_string(i+1) + ","; // +1 to switch from base-0 to base-1 so the node nums match with resulting zone
				}
				AuxDataZoneSetItem(SurfZoneNum, CSMAuxData.GBA.SphereConstrainedNodeNums, MovedPointNumsStr);
				AuxDataZoneSetItem(SurfZoneNum, CSMAuxData.GBA.SphereConstrainedNodeIntersectCPTypes, MovedPointCPTypes);
				AuxDataZoneSetItem(SurfZoneNum, CSMAuxData.GBA.SphereConstrainedNodeIntersectCPNames, MovedPointCPNames);
				AuxDataZoneSetItem(SurfZoneNum, CSMAuxData.GBA.SphereConstrainedNodeIntersectCPTotalOffsetNames, MovedPointCPNamesTotalCount);
				AuxDataZoneSetItem(SurfZoneNum, CSMAuxData.GBA.NumGBs, to_string(NumTriangles));
				AuxDataZoneSetItem(SurfZoneNum, CSMAuxData.GBA.GPsPerGB, to_string(3 + 3 * NumEdgeGPs));
				AuxDataZoneSetItem(SurfZoneNum, CSMAuxData.GBA.PointsPerGP, to_string(NumSTPoints));
				AuxDataZoneSetItem(SurfZoneNum, CSMAuxData.GBA.SourceZoneNum, to_string(CPZoneNum));
				AuxDataZoneSetItem(SurfZoneNum, CSMAuxData.GBA.SourceNucleusName, NucleusName);
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
				StatusDrop(AddOnID);
				TecUtilDialogErrMsg("Couldn't find gradient vector variables.");
				TecUtilDataLoadEnd();
				TecUtilLockFinish(AddOnID);
				return;
			}
		}

// 		if (IsOk){
// 
// 			// Enable system volume zone
// 
// 			Set_pa TmpSet = TecUtilSetAlloc(FALSE);
// 			TecUtilSetAddMember(TmpSet, 1, FALSE);
// 			TecUtilZoneSetActive(TmpSet, AssignOp_PlusEquals);
// 			TecUtilSetDealloc(&TmpSet);
// 
// 		}

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
// 			NodeHasSaddleCP[IntZoneSaddleCPNodeNums[i][3]] = TRUE;
			NodeHasSaddleCP[IntZoneSaddleCPNodeNums[i][3]] = MovedCPIsBond[i];
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


		StatusDrop(AddOnID);


		if (IsOk){
			TmpString = ProgressStr.str() + string(" ... Step 2 of 4 ... Seeding gradient paths");
			StatusLaunch(TmpString, AddOnID, TRUE);
		}

		Boolean_t UserQuit = FALSE;

		int TmpNumIterations = NumPoints + NumEdges * NumEdgeGPs;
		int NumCompleted = 0;

		high_resolution_clock::time_point StatusStartTime = high_resolution_clock::now();

// #ifndef _DEBUG
#pragma omp parallel for schedule(dynamic)
// #endif
		for (int PtNum = 0; PtNum < NumPoints; ++PtNum){
			if (omp_get_thread_num() == 0){
				UserQuit = !StatusUpdate(NumCompleted, TmpNumIterations, TmpString, AddOnID, StatusStartTime);
#pragma omp flush (UserQuit)
			}
#pragma omp flush (UserQuit)
#pragma omp atomic
			NumCompleted++;

			if (!NodeHasSaddleCP[PtNum] && !UserQuit){
				vec3 NodePos;
				for (int ii = 0; ii < 3; ++ii)
					NodePos[ii] = p[PtNum][ii];

// 				IsOk = GPsNonSaddle[PtNum].SetupGradPath(NodePos,
// 					StreamDir, 
// 					NumSTPoints,
// 					HowTerminate, 
// 					NULL, vector<FieldDataPointer_c>(), NULL, NULL,
// 					&CutoffVal, 
// 					VolMaxIJK,
// 					VolMaxXYZ,
// 					VolMinXYZ,
// 					GradRawPtrs, 
// 					RhoRawPtr);

				IsOk = GPsNonSaddle[PtNum].SetupGradPath(NodePos, 
					StreamDir, 
					NumSTPoints, 
					GPType_Classic, 
					HowTerminate, 
					NULL, NULL, NULL, 
					&CutoffVal, 
					VolInfo, 
					vector<FieldDataPointer_c>(), 
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

// #ifndef _DEBUG
#pragma omp parallel for schedule(dynamic)
// #endif
		for (int EdgeNum = 0; EdgeNum < NumEdges; ++EdgeNum){
			if (omp_get_thread_num() == 0){
				UserQuit = !StatusUpdate(NumCompleted, TmpNumIterations, TmpString, AddOnID, StatusStartTime);
#pragma omp flush (UserQuit)
			}
#pragma omp flush (UserQuit)
#pragma omp atomic
			NumCompleted++;
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

// 					IsOk = GPsEdges[GPInd].SetupGradPath(SeedPt,
// 						StreamDir,
// 						NumSTPoints,
// 						HowTerminate,
// 						NULL, vector<FieldDataPointer_c>(), NULL, NULL,
// 						&CutoffVal,
// 						VolMaxIJK,
// 						VolMaxXYZ,
// 						VolMinXYZ,
// 						GradRawPtrs,
// 						RhoRawPtr);

					IsOk = GPsEdges[GPInd].SetupGradPath(SeedPt,
						StreamDir,
						NumSTPoints,
						GPType_Classic,
						HowTerminate,
						NULL, NULL, NULL,
						&CutoffVal,
						VolInfo,
						vector<FieldDataPointer_c>(),
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
			StatusDrop(AddOnID);

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
		vector<vector<GradPath_c> > GPsSaddle(IntZoneSaddleCPNodeNums.size(),vector<GradPath_c>(6)), GPsSaddleEdges(IntZoneSaddleCPNodeNums.size(), vector<GradPath_c>(6*NumEdgeGPs));
		vector<vector<int> > GPsSaddleClosestPtNum(IntZoneSaddleCPNodeNums.size(), vector<int>(6)), GPsNonSaddleClosestPtNums(IntZoneSaddleCPNodeNums.size(), vector<int>(6));
		ConstrainedNeighborhoodEdgeNums.reserve(6);
		int NumIntZoneSaddleCPNodeNums = IntZoneSaddleCPNodeNums.size();
		
#ifdef _DEBUG
		vector<vector<int> > EdgeList(NumEdges, vector<int>(2));
#pragma omp parallel for
		for (int i = 0; i < NumEdges; ++i){
			for (int j = 0; j < 2; ++j){
				EdgeList[i][j] = e[i][j];
			}
		}
#endif // _DEBUG


// #ifndef _DEBUG
 #pragma omp parallel for schedule(dynamic)
// #endif // !_DEBUG
		for (int i = 0; i < NumIntZoneSaddleCPNodeNums; ++i){
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

					//GPsSaddle[i].resize(ConstrainedNeighborNodesNum[j].size());
					//GPsSaddleEdges[i].resize(ConstrainedNeighborNodesNum[j].size() * NumEdgeGPs);
					//ConstrainedNeighborhoodEdgeNums.reserve(ConstrainedNeighborNodesNum[j].size());
					//GPsSaddleClosestPtNum[i].resize(ConstrainedNeighborNodesNum[j].size());
					//GPsNonSaddleClosestPtNums[i].resize(ConstrainedNeighborNodesNum[j].size());

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

// 						IsOk = GPsSaddle[i][k].SetupGradPath(SeedPts[k], StreamDir, NumSTPoints, HowTerminate,
// 							NULL, vector<FieldDataPointer_c>(), NULL, NULL, &CutoffVal,
// 							VolMaxIJK, VolMaxXYZ, VolMinXYZ,
// 							GradRawPtrs,
// 							RhoRawPtr);

						IsOk = GPsSaddle[i][k].SetupGradPath(SeedPts[k],
							StreamDir,
							NumSTPoints,
							GPType_Classic,
							HowTerminate,
							NULL, NULL, NULL,
							&CutoffVal,
							VolInfo,
							vector<FieldDataPointer_c>(),
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

// 									IsOk = GPsSaddleEdges[i][GPInd].SetupGradPath(SeedPt, StreamDir, NumSTPoints, HowTerminate,
// 										NULL, vector<FieldDataPointer_c>(), NULL, NULL, &CutoffVal,
// 										VolMaxIJK, VolMaxXYZ, VolMinXYZ,
// 										GradRawPtrs,
// 										RhoRawPtr);

									IsOk = GPsSaddleEdges[i][GPInd].SetupGradPath(SeedPt,
										StreamDir,
										NumSTPoints,
										GPType_Classic,
										HowTerminate,
										NULL, NULL, NULL,
										&CutoffVal,
										VolInfo,
										vector<FieldDataPointer_c>(),
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
// 

		if (SaveGPs){
			Set ZoneSet;
			for (int i = 0; i < GPsNonSaddle.size(); ++i){
				if (GPsNonSaddle[i].IsMade()){
					GPsNonSaddle[i].SaveAsOrderedZone("GP " + CPName + " " + "Node " + to_string(i + 1));
					ZoneSet += GPsNonSaddle[i].GetZoneNum();
					AuxDataZoneSetItem(GPsNonSaddle[i].GetZoneNum(), GBASphereCPName, CPString);
					AuxDataZoneSetItem(GPsNonSaddle[i].GetZoneNum(), GBASphereCPNum, to_string(CPNum));
					AuxDataZoneSetItem(GPsNonSaddle[i].GetZoneNum(), GBAGPNodeNum, to_string(i + 1));
					AuxDataZoneSetItem(GPsNonSaddle[i].GetZoneNum(), GBAZoneType, GBAZoneTypeGradPath);
				}
			}

			for (int i = 0; i < GPsSaddle.size(); ++i){
				for (int j = 0; j < GPsSaddle[i].size(); ++j){
					if (GPsSaddle[i][j].IsMade()){
						GPsSaddle[i][j].SaveAsOrderedZone("GP " + CPName + " " + "Node " + to_string(MovedPointNums[i] + 1) + "." + to_string(j + 1));
						ZoneSet += GPsSaddle[i][j].GetZoneNum();
						AuxDataZoneSetItem(GPsSaddle[i][j].GetZoneNum(), GBASphereCPName, CPString);
						AuxDataZoneSetItem(GPsSaddle[i][j].GetZoneNum(), GBASphereCPNum, to_string(CPNum));
						AuxDataZoneSetItem(GPsSaddle[i][j].GetZoneNum(), GBAGPNodeNum, to_string(MovedPointNums[i] + 1));
						AuxDataZoneSetItem(GPsSaddle[i][j].GetZoneNum(), GBAZoneType, GBAZoneTypeGradPath);
					}
				}
			}

			for (int i = 0; i < GPsEdges.size(); ++i){
				if (GPsEdges[i].IsMade()){
					GPsEdges[i].SaveAsOrderedZone("GP " + CPName + " " + "EdgeGP " + to_string(i + 1), Blue_C);
					ZoneSet += GPsEdges[i].GetZoneNum();
					AuxDataZoneSetItem(GPsEdges[i].GetZoneNum(), GBASphereCPName, CPString);
					AuxDataZoneSetItem(GPsEdges[i].GetZoneNum(), GBASphereCPNum, to_string(CPNum));
					AuxDataZoneSetItem(GPsEdges[i].GetZoneNum(), GBAGPNodeNum, to_string(i + 1));
					AuxDataZoneSetItem(GPsEdges[i].GetZoneNum(), GBAZoneType, GBAZoneTypeGradPath);
				}
			}

			for (int i = 0; i < GPsSaddleEdges.size(); ++i){
				for (int j = 0; j < GPsSaddleEdges[i].size(); ++j){
					if (GPsSaddleEdges[i][j].IsMade()){
						GPsSaddleEdges[i][j].SaveAsOrderedZone("GP " + CPName + " " + "Saddle Edge " + to_string(MovedPointNums[i] + 1) + "." + to_string(j + 1));
						ZoneSet += GPsSaddleEdges[i][j].GetZoneNum();
						AuxDataZoneSetItem(GPsSaddleEdges[i][j].GetZoneNum(), GBASphereCPName, CPString);
						AuxDataZoneSetItem(GPsSaddleEdges[i][j].GetZoneNum(), GBASphereCPNum, to_string(CPNum));
						AuxDataZoneSetItem(GPsSaddleEdges[i][j].GetZoneNum(), GBAGPNodeNum, to_string(MovedPointNums[i] + 1));
						AuxDataZoneSetItem(GPsSaddleEdges[i][j].GetZoneNum(), GBAZoneType, GBAZoneTypeGradPath);
					}
				}
			}
			TecUtilZoneSetActive(ZoneSet.getRef(), AssignOp_MinusEquals);
		}


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
		vector<FESurface_c> FEVolumes(NumTriangles);
		vector<bool> IsDivergentGB(NumTriangles, false);

		TmpString = ProgressStr.str() + string(" ... Step 3 of 4 ... Integrating volumes");
		TmpNumIterations = NumTriangles;
// #ifndef _DEBUG
// 			TmpNumIterations /= numCPU;
// #endif // !_DEBUG

		vector<vector<double> > IntVals(NumTriangles, vector<double>(IntVarNumList.size() + 1, 0));
			
		NumCompleted = 0;

		vector<vector<int> > GPSaddleNodeNums1(NumTriangles), GPSaddleNodeNums2(NumTriangles), GPNonSaddleNodeNums(NumTriangles);
		for (int i = 0; i < NumTriangles; ++i){
			GPSaddleNodeNums1[i].reserve(3);
			GPSaddleNodeNums2[i].reserve(3);
			GPNonSaddleNodeNums[i].reserve(3);
		}

// 		TecUtilDialogMessageBox("Before making FEVolume", MessageBoxType_Information);

		vector<VolExtentIndexWeights_s> ThreadVolInfo(omp_get_num_procs(), VolInfo);

		StatusStartTime = high_resolution_clock::now();

		TecUtilDataLoadBegin();
// #ifndef _DEBUG
#pragma omp parallel for schedule(dynamic)
// #endif
		for (int TriNum = 0; TriNum < NumTriangles; ++TriNum){
// 		for (int TriNum = 24; TriNum < 25; ++TriNum){
			if (omp_get_thread_num() == 0){
				UserQuit = !StatusUpdate(NumCompleted, TmpNumIterations, TmpString, AddOnID, StatusStartTime);
#pragma omp flush (UserQuit)
			}
#pragma omp flush (UserQuit)
			if (!UserQuit){
				vector<GradPath_c*> GPPtrs;
				if (TriangleHasSaddleCP[TriNum]){
					GPPtrs.reserve(4 + 4 * NumEdgeGPs);
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
								if (GPsSaddle[SaddleNum][ConstrainedGPNums[1]].IsMade()) {
									GPPtrs.push_back(&GPsSaddle[SaddleNum][ConstrainedGPNums[1]]);
									bool EdgeFound = false;
									for (int ei = 0; ei < ConstrainedNeighborhoodEdgeNums[SaddleNum].size() && !EdgeFound; ++ei) {
										if (std::find(TriangleEdgeNums[TriNum].begin(), TriangleEdgeNums[TriNum].end(), ConstrainedNeighborhoodEdgeNums[SaddleNum][ei]) != TriangleEdgeNums[TriNum].end()) {
											if (e[ConstrainedNeighborhoodEdgeNums[SaddleNum][ei]][0] == NonConstrainedNodeNums[1]) {
												for (int i = 0; i < NumEdgeGPs; ++i) {
													if (GPsSaddleEdges[SaddleNum][ei * NumEdgeGPs + i].IsMade())
														GPPtrs.push_back(&GPsSaddleEdges[SaddleNum][ei * NumEdgeGPs + i]);
												}
											}
											else {
												for (int i = NumEdgeGPs - 1; i >= 0; --i) {
													if (GPsSaddleEdges[SaddleNum][ei * NumEdgeGPs + i].IsMade())
														GPPtrs.push_back(&GPsSaddleEdges[SaddleNum][ei * NumEdgeGPs + i]);
												}
											}
										}
									}
									GPPtrs.push_back(&GPsSaddle[SaddleNum][ConstrainedGPNums[0]]);
								}
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
							
						/*IsOk = FEVolumes[TriNum].MakeGradientBundle(GPsNonSaddle[NonConstrainedNodeNums[0]],
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
					/*IsOk = FEVolumes[TriNum].MakeGradientBundle(GPsNonSaddle[t[TriNum][0]],
						GPsNonSaddle[t[TriNum][1]],
						GPsNonSaddle[t[TriNum][2]]);*/
					if (IsOk){
						for (int i = 0; i < 3; ++i){
							GPNonSaddleNodeNums[TriNum].push_back(t[TriNum][i]);
						}
					}
				}
				if (IsOk){
					/*
					 *	As a workaround until ring surface divergent GBs are handled correctly,
					 *	Use the old integration method for divergent GBs, which are determined by
					 *	checking for a dot product of ~-1 between the last line segments of two
					 *	GPs.
					 */
					vector<vec3*> TerminalVecs(GPPtrs.size(), NULL);
					for (int i = 0; i < GPPtrs.size() - 1 && !IsDivergentGB[TriNum]; ++i){
						if (TerminalVecs[i] == NULL){
							TerminalVecs[i] = new vec3; 
							*TerminalVecs[i] = normalise(GPPtrs[i]->XYZAt(-1) - GPPtrs[i]->XYZAt(-2));
						}
						for (int j = i + 1; j < GPPtrs.size() && !IsDivergentGB[TriNum]; ++j){
							if (TerminalVecs[j] == NULL){
								TerminalVecs[j] = new vec3;
								*TerminalVecs[j] = normalise(GPPtrs[j]->XYZAt(-1) - GPPtrs[j]->XYZAt(-2));
							}
							IsDivergentGB[TriNum] = (dot(*TerminalVecs[i], *TerminalVecs[j]) < DivergentGPMaxTerminalDotProduct);
						}
					}
					for (auto * i : TerminalVecs) if (i != NULL) delete i;
					if (SaveGBs || IsDivergentGB[TriNum]){
						// 					if (TriNum < 3) TecUtilDialogMessageBox("Before making FEVolume", MessageBoxType_Information);
						IsOk = FEVolumes[TriNum].MakeFromGPs(GPPtrs, true, true);
						// 					if (TriNum < 3) TecUtilDialogMessageBox("After making FEVolume", MessageBoxType_Information);
					}
					if (!IsDivergentGB[TriNum] && DoIntegration)
						IntegrateUsingIsosurfaces(GPPtrs, IntResolution, ThreadVolInfo[omp_get_thread_num()], IntVarPtrs, IntVals[TriNum], TriangleHasSaddleCP[TriNum]);
				}
#pragma omp atomic
				NumCompleted++;
			}
		}


		/*
		 *	Now collect and save integration values
		 */
		NumZones = TecUtilDataSetGetNumZones();
		vector<string> NewVarNames;
		vector<int> NewVarNums;
		vector<FieldDataType_e> ZoneDataTypes(NumZones, FieldDataType_Bit);

		if (DoIntegration) {

			for (int i = 1; i <= NumZones; ++i) {
				if (AuxDataZoneHasItem(i, CSMAuxData.GBA.ZoneType)) {
					ZoneDataTypes[i - 1] = FieldDataType_Double;
				}
			}
			vector<ValueLocation_e> ZoneDataLocs(NumZones, ValueLocation_CellCentered);
			// 	vector<ValueLocation_e> ZoneDataLocs(NumZones, ValueLocation_Nodal);
			vector<string> VarNameAppends = { "I", "N", "S" };
			for (int i = 0; i < IntVarNameList.size() && IsOk; ++i) {
				for (int j = 0; j < 3; ++j) {
					string TmpStr;
					for (int k = 0; k <= j; ++k)
						TmpStr += VarNameAppends[k];
					NewVarNames.push_back(TmpStr + ": " + IntVarNameList[i]);
					NewVarNums.push_back(VarNumByName(NewVarNames.back()));
					if (NewVarNums[i * 3 + j] < 0) {
						ArgList_pa ArgList = TecUtilArgListAlloc();

						IsOk = TecUtilArgListAppendString(ArgList, SV_NAME, NewVarNames.back().c_str());
						if (IsOk)
							IsOk = TecUtilArgListAppendArray(ArgList, SV_VARDATATYPE, ZoneDataTypes.data());
						if (IsOk)
							IsOk = TecUtilArgListAppendArray(ArgList, SV_VALUELOCATION, ZoneDataLocs.data());

						if (IsOk)
							IsOk = TecUtilDataSetAddVarX(ArgList);

						if (IsOk)
							NewVarNums.back() = TecUtilDataSetGetNumVars();

						if (IsOk) {
							Set_pa TmpSet = TecUtilSetAlloc(FALSE);
							IsOk = TecUtilSetAddMember(TmpSet, NewVarNums.back(), FALSE);
							if (IsOk)
								TecUtilStateChanged(StateChange_VarsAdded, (ArbParam_t)TmpSet);
							TecUtilSetDealloc(&TmpSet);
						}

						TecUtilArgListDealloc(&ArgList);
					}
				}
			}
			// 			TecUtilDialogMessageBox((VectorToString(NewVarNames) + ": " + VectorToString(NewVarNums)).c_str(), MessageBoxType_Information);
			if (IsOk) {
				for (int j = 0; j < 3; ++j) {
					string TmpStr;
					for (int k = 0; k <= j; ++k)
						TmpStr += VarNameAppends[k];
					TmpStr += ": Volume";
					if (std::find(NewVarNames.begin(), NewVarNames.end(), TmpStr) == NewVarNames.end()) {
						NewVarNames.push_back(TmpStr);
						NewVarNums.push_back(VarNumByName(NewVarNames.back()));
						if (NewVarNums.back() < 0) {
							ArgList_pa ArgList = TecUtilArgListAlloc();

							IsOk = TecUtilArgListAppendString(ArgList, SV_NAME, NewVarNames.back().c_str());
							if (IsOk)
								IsOk = TecUtilArgListAppendArray(ArgList, SV_VARDATATYPE, ZoneDataTypes.data());
							if (IsOk)
								IsOk = TecUtilArgListAppendArray(ArgList, SV_VALUELOCATION, ZoneDataLocs.data());

							if (IsOk)
								IsOk = TecUtilDataSetAddVarX(ArgList);

							if (IsOk)
								NewVarNums.back() = TecUtilDataSetGetNumVars();

							if (IsOk) {
								Set_pa TmpSet = TecUtilSetAlloc(FALSE);
								IsOk = TecUtilSetAddMember(TmpSet, NewVarNums.back(), FALSE);
								if (IsOk)
									TecUtilStateChanged(StateChange_VarsAdded, (ArbParam_t)TmpSet);
								TecUtilSetDealloc(&TmpSet);
							}

							TecUtilArgListDealloc(&ArgList);
						}
					}
				}
			}

			// 			TecUtilDialogMessageBox((VectorToString(NewVarNames) + ": " + VectorToString(NewVarNums)).c_str(), MessageBoxType_Information);
			vector<string> AuxDataNames(IntVarNameList.size());
			for (int i = 0; i < IntVarNameList.size(); ++i)
				AuxDataNames[i] = RemoveStringChar(IntVarNameList[i], " ");
			AuxDataNames.push_back(string("Volume"));

			Sphere = FESurface_c(Sphere.GetZoneNum(), VolZoneNum, XYZVarNums, IntVarNumList, true);

			/*
			 *	Save any divergent GBs as zones in order to integrate them, then integrate them.
			 *	Also integrate the sphere here, but only if divergent GBs were found.
			 */
			bool DivergentGBsExist = false;
			for (int i = 0; i < NumTriangles && !DivergentGBsExist; ++i) {
				DivergentGBsExist = IsDivergentGB[i];
			}
			if (DivergentGBsExist) {
				TecUtilDataLoadBegin();
				vector<FESurface_c> TempGBs(NumTriangles);
				int NumVolumes = 0;
				Set TmpSet;
				for (int i = 0; i < NumTriangles; ++i) {
					if (IsDivergentGB[i] && FEVolumes[i].IsMade()) {
						DivergentGBsExist = true;
						int TmpZoneNum = FEVolumes[i].SaveAsTriFEZone(XYZVarNums, "Temp GB" + to_string(i + 1));
						TempGBs[i] = FESurface_c(TmpZoneNum, VolZoneNum, XYZVarNums, IntVarNumList);
						TmpSet += TmpZoneNum;
						NumVolumes++;
					}
				}
				NumVolumes++; // +1 for the Sphere integration

				TecUtilZoneSetActive(TmpSet.getRef(), AssignOp_PlusEquals);
				TecUtilRedrawAll(TRUE);

				/*
				 *	Now integrate all the saved zones
				 */
				int NumComplete = 0;
				StatusStartTime = high_resolution_clock::now();
				// #ifndef _DEBUG
#pragma omp parallel for schedule(dynamic)
// #endif
				for (int i = 0; i < NumTriangles; ++i) {
					if (omp_get_thread_num() == 0)
					{
						stringstream ss;
						ss << "Integrating other volumes " << NumComplete + 1;
						if (AtomNameList.size() > 1)
							ss << " of " << NumVolumes;
						UserQuit = !StatusUpdate(NumComplete, NumVolumes, ss.str(), AddOnID, StatusStartTime);
#pragma omp flush (UserQuit)
					}
#pragma omp flush (UserQuit)
					if (!UserQuit){
						if (TempGBs[i].IsMade()) {
							TempGBs[i].DoIntegrationNew(IntResolution - 1, TRUE);
							IntVals[i] = TempGBs[i].GetIntResults();
#pragma omp atomic
							NumComplete++;
						}
						if (!Sphere.IntResultsReady() && omp_get_thread_num() == 0) {
							Sphere.DoIntegrationNew(IntResolution - 1, TRUE);
#pragma omp atomic
							NumComplete++;
						}
					}

				}

				/*
				 *	Now delete the zones
				 */
				Set DelZones;
				for (const auto & i : TempGBs) {
					if (i.IsMade())
						DelZones += i.GetZoneNum();
				}
				TecUtilDataSetDeleteZone(DelZones.getRef());



				TecUtilDataLoadEnd();
			}
			else {
				stringstream ss;
				ss << "Integrating sphere";
				TecUtilPleaseWait(ss.str().c_str(), TRUE);
				Sphere.DoIntegrationNew(IntResolution - 1, TRUE);
				TecUtilPleaseWait(ss.str().c_str(), FALSE);
			}



			vector<double> SphereTriangleAreas;
			vector<vector<double> > SphereElemIntVals = Sphere.GetTriSphereIntValsByElem(&SphereTriangleAreas);
			for (int i = 0; i < NumTriangles; ++i) {
				for (int j = 0; j < SphereElemIntVals[i].size(); ++j) {
					SphereElemIntVals[i][j] += IntVals[i][j];
				}
			}
			vector<vector<double> > NormalizedValues = SphereElemIntVals;
			vector<double> TotalList(SphereElemIntVals[0].size(), 0.0),
				TotalNormlizedList(SphereElemIntVals[0].size(), 0.0),
				IntScaleFactors(SphereElemIntVals[0].size());
			for (int i = 0; i < SphereElemIntVals.size(); ++i) {
				for (int j = 0; j < SphereElemIntVals[i].size(); ++j) {
					NormalizedValues[i][j] /= SphereTriangleAreas[i];
					TotalNormlizedList[j] += NormalizedValues[i][j];
					TotalList[j] += SphereElemIntVals[i][j];
				}
			}
			for (int i = 0; i < SphereElemIntVals[0].size(); ++i) {
				IntScaleFactors[i] = TotalList[i] / TotalNormlizedList[i];
			}
			for (int i = 0; i < SphereElemIntVals.size(); ++i) {
				vector<double> TmpVec = SphereElemIntVals[i];
				SphereElemIntVals[i] = vector<double>();
				SphereElemIntVals[i].reserve(3 * TmpVec.size());
				for (int j = 0; j < TmpVec.size(); ++j) {
					SphereElemIntVals[i].push_back(TmpVec[j]);
					SphereElemIntVals[i].push_back(NormalizedValues[i][j]);
					SphereElemIntVals[i].push_back(NormalizedValues[i][j] * IntScaleFactors[j]);
				}
			}

			if (std::find(IntVarNameList.begin(), IntVarNameList.end(), "Volume") == IntVarNameList.end())
				IntVarNameList.push_back("Volume");
			AuxDataZoneSetItem(Sphere.GetZoneNum(), CSMAuxData.GBA.AtomicBasinIntegrationVariables, VectorToString(IntVarNameList, ","));
			AuxDataZoneSetItem(Sphere.GetZoneNum(), CSMAuxData.GBA.AtomicBasinIntegrationValues, VectorToString(TotalList, ","));
			AuxDataZoneSetItem(Sphere.GetZoneNum(), CSMAuxData.GBA.SphereRadius, to_string(Radius));
			AuxDataZoneSetItem(Sphere.GetZoneNum(), CSMAuxData.GBA.SphereRadiusCPDistRatio, to_string(UserRadius));


			TecUtilDataLoadBegin();

			// 		if (SphereElemIntVals[0].size() != NewVarNums.size()) {
			// 			TecUtilDialogErrMsg("Discrepency between number of integrated variables to be saved and the number of variable numbers to use when saving.");
			// 		}

			vector<FieldDataPointer_c> Ptrs(SphereElemIntVals[0].size());
			for (int j = 0; j < Ptrs.size(); ++j) {
				Ptrs[j].GetWritePtr(Sphere.GetZoneNum(), NewVarNums[j]);
			}

			// Write cell-centered integration values to sphere and FE zones
#pragma omp parallel for
			for (int e = 0; e < NumTriangles; ++e) {
				for (int j = 0; j < SphereElemIntVals[e].size(); ++j) {
					Ptrs[j].Write(e, SphereElemIntVals[e][j]);
				}
			}
			for (auto & i : Ptrs)
				i.Close();



			TecUtilDataLoadEnd();
		}

		if (UserQuit){
			StatusDrop(AddOnID);

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


		if (SaveGBs){
			TmpString = ProgressStr.str() + string(" ... Step 4 of 4 ... Saving volumes");
			TecUtilDataLoadBegin();
			vector<EntIndex_t> TriangleFEZoneNums(NumTriangles, -1);
			for (int TriNum = 0; TriNum < NumTriangles && IsOk; ++TriNum){
				if (!StatusUpdate(TriNum, NumTriangles, TmpString, AddOnID)){
					StatusDrop(AddOnID);

					MemoryRequired = sizeof(point)* NumPoints
						+ sizeof(triangle)* NumTriangles
						+ sizeof(int)* 2 * NumEdges;
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

				//TriangleFEZoneNums[TriNum] = FEVolumes[TriNum].SaveAsFEZone(VarDataTypes, XYZVarNums, CutoffVarNum);

				// 			if (TriNum < 3) TecUtilDialogMessageBox("before saving FEVolume", MessageBoxType_Information);
				//TriangleFEZoneNums[TriNum] = FEVolumes[TriNum].SaveAsTriFEZone("Gradient Bundle " + to_string(TriNum + 1) + " (Zone " + to_string(TecUtilDataSetGetNumZones() + 1) + ")", VarDataTypes, VarLocations, XYZVarNums);
				TriangleFEZoneNums[TriNum] = FEVolumes[TriNum].SaveAsTriFEZone("Gradient Bundle " + to_string(TriNum + 1) + " (Zone " + to_string(TecUtilDataSetGetNumZones() + 1) + ")", vector<FieldDataType_e>(), vector<ValueLocation_e>(), XYZVarNums);
				// 			if (TriNum < 3) TecUtilDialogMessageBox("after saving FEVolume", MessageBoxType_Information);

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
		}

		


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
			StatusUpdate(i, NewNumZones, TmpString, AddOnID);
			TecUtilDrawGraphics(FALSE);
			for (EntIndex_t j = 0; j < 3 && IsOk; ++j){
				IsOk = TecUtilDataValueUnload(i, XYZVarNums[j]);
			}
		}
		StatusDrop(AddOnID);

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

	if (deleteOldSphereZones)
		TecUtilDataSetDeleteZone(oldSphereZonesToDelete.getRef());

	TecUtilArrayDealloc((void **)&CPNums);

	TecGUITabSetCurrentPage(TAB1_TB_D1, 2);

// 	TecUtilDialogMessageBox("Finished making GBs", MessageBoxType_Information);

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
	vector<FESurface_c> VolumeList;
	VolumeList.reserve(NumZones);
	string TmpString = "Atom 1";
	for (int i = 1; i <= NumZones; ++i)
		if (TecUtilZoneIsFiniteElement(i)
			&& AuxDataZoneItemMatches(i, CSMAuxData.GBA.SphereCPName, TmpString))
		{
			VolumeList.push_back(FESurface_c(i, VolZoneNum, XYZVarNums, IntVarList));
		}


	StatusLaunch("Integrating", AddOnID, TRUE);

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
				StatusUpdate(i * numCPU, static_cast<int>(VolumeList.size()), string(ProgressBase.str() + to_string(i * numCPU) + string(" of ") + to_string(VolumeList.size())), AddOnID);
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

	StatusDrop(AddOnID);

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

void ElemWiseGradPath(const int & StartingElem, 
						const int & DirNum, // 0 for downhill, 1 for uphill
						const vector<vector<int> > & ElemConnectivity,
						const vector<vector<double> > & ElemIntVals, 
						const int & BasinVarNum, 
						vector<vector<int> > & TerminalMinMaxIndices)
{
	int CurrentElem = StartingElem;
	int NextElem = CurrentElem;
	vector<int> IntermediateElems;
	IntermediateElems.reserve(100);
	/*
	 *	Starting at StartingElem, find the neighboring element NextElem with the least or greatest value
	 *	(depending on the value of DirNum), then step to that element.
	 *	If NextElem is CurrentElem then we've arrived at the critical element,
	 *	or if TerminalMinMaxIndices[DirNum][CurrentElem] >= 0 then we've reach an element that has already
	 *	been checked, hence the rest of the path has already been found and recorded.
	 *	We save all the intermediate elements' indices, so the last index recorded is the index of the
	 *	min or max.
	 */
	int iter = 0;
	while (iter < ElemConnectivity.size()) {
		IntermediateElems.push_back(CurrentElem);
		/*
		 *	Check the neighboring elements' values and update NextElem when the lower (or higher) value is found,
		 *	which will change the comparison for the next iteration because NextElem is used to lookup the "old" value.
		 */
		for (int i = 0; i < ElemConnectivity[CurrentElem].size(); ++i) {
			if ((DirNum == 0 && ElemIntVals[BasinVarNum][ElemConnectivity[CurrentElem][i]] < ElemIntVals[BasinVarNum][NextElem])
				|| (DirNum == 1 && ElemIntVals[BasinVarNum][ElemConnectivity[CurrentElem][i]] > ElemIntVals[BasinVarNum][NextElem]))
			{
				NextElem = ElemConnectivity[CurrentElem][i];
			}
		}
		if (NextElem == CurrentElem) 
			break;
		else if (TerminalMinMaxIndices[NextElem][DirNum] >= 0) {
			IntermediateElems.push_back(TerminalMinMaxIndices[NextElem][DirNum]);
			break;
		}
		else 
			CurrentElem = NextElem;

		iter++;
	}
	for (const int & i : IntermediateElems) TerminalMinMaxIndices[i][DirNum] = IntermediateElems.back();

	REQUIRE(iter < ElemConnectivity.size());

	return;
}

void GetSphereBasinIntegrations(const FESurface_c & Sphere,
									const vector<vector<double> > & ElemIntVals,
									const int & BasinVarNum,
									vector<vector<int> > & MinMaxIndices,
									vector<vector<vector<double> > > & BasinIntVals,
									vector<vector<vector<vec3> > > & BasinNodes,
									vector<vector<vector<vector<int> > > > & BasinElems,
									vector<vector<vector<int> > > & SphereNodeNums,
									vector<vector<vector<int> > > & SphereElemNums)
{
	REQUIRE(Sphere.IsMade());

	/*
		First, use the node connectivity list already inside the Sphere object
		to generate an element connectivity list.
		Two elements are neighbors if they share an edge.
	*/
	auto * NodeConnectivityPtr = Sphere.GetConnectivityListPtr();
	auto * ElemListPtr = Sphere.GetElemListPtr();
	auto * XYZListPtr = Sphere.GetXYZListPtr();

	REQUIRE(NodeConnectivityPtr != NULL && ElemListPtr != NULL && XYZListPtr != NULL);

	int NumElems = ElemListPtr->size();

	/*
	 *	ElemIntVals is here number of integrated variables X number of sphere elements,
	 *	which is the transpose of its shape when returned from Sphere.GetTriSphereIntValsByElem()
	 */
	int NumIntVars = ElemIntVals.size();
	REQUIRE(BasinVarNum >= 0 && BasinVarNum < NumIntVars);
	REQUIRE(ElemIntVals[0].size() == NumElems);

	vector<vector<int> > ElemConnectivity;
	GetTriElementConnectivityList(ElemListPtr, ElemConnectivity, 1);

	/*
	 *	Now I can start doing element-wise gradient paths to find all the local minima and maxima.
	 *	For each element I'll store the indices of the elements you terminate at by
	 *	going "uphill" or "downhill" through whatever values are selected (e.g. condensed charge density).
	 *	Initialize the whole to -1 and do this in parallel.
	 *	While tracing a gradient path, keep track of all the elements you traversed to get to the
	 *	terminal element, then set the terminal index for all the intermediate elements too.
	 */

	/*
	 *	There can be noise that results in spurious local min/max elements.
	 *	We'll loop, applying smoothing to the integration values used to find the basins,
	 *	and when the number of min/max elems found converges, assume that's correct
	 *	(the actual min/max elem numbers might change due to the smoothing, so only use
	 *	number of min/max elems found to determine convergence.
	 *	When convergence in achieved, use the older results (with one less round of
	 *	smoothing) because they're less altered.
	 *	Limit to MaxNumSmoothing rounds of smoothing.
	 */
	int MaxNumSmoothing = 20;
	int NumConvergedIter = 5;
	vector<int> NumMinMax(2,-1), OldNumMinMax(2,INT_MAX);
	MinMaxIndices.resize(2);
	vector<vector<int> > OldMinMaxIndices(2);
	vector<vector<int> > TerminalMinMaxIndices, OldTerminalMinMaxIndices;
	vector<vector<double> > SmoothElemIntVals, TmpSmoothElemIntVals;
	SmoothElemIntVals.push_back(ElemIntVals[BasinVarNum]);

	int Iter = 0, convegedIter = 0;
	while (Iter < MaxNumSmoothing && 
		convegedIter < NumConvergedIter) {
		Iter++;
		OldNumMinMax = NumMinMax;
		OldMinMaxIndices = MinMaxIndices;
		MinMaxIndices.assign(2, vector<int>());
		OldTerminalMinMaxIndices = TerminalMinMaxIndices;
		TerminalMinMaxIndices.assign(ElemConnectivity.size(), vector<int>(2, -1));
		/*
		 *	Do a simple smoothing on the ElemIntVals based on average neighborhood value
		 */
		TmpSmoothElemIntVals = SmoothElemIntVals;
#pragma omp parallel for
		for (int i = 0; i < NumElems; ++i) {
			for (const int & j : ElemConnectivity[i]) {
				SmoothElemIntVals[0][i] += TmpSmoothElemIntVals[0][j];
			}
			SmoothElemIntVals[0][i] /= double(ElemConnectivity[i].size() + 1);
		}

#pragma omp parallel for schedule(dynamic)
		for (int i = 0; i < NumElems; ++i) {
			for (int DirNum = 0; DirNum < 2; ++DirNum)
				ElemWiseGradPath(i, DirNum, ElemConnectivity, SmoothElemIntVals, 0, TerminalMinMaxIndices);
		}

		/*
		 *	Now all the terminal elements have been identified.
		 *	Get a list of the unique min and max element indices.
		 */
#pragma omp parallel for
		for (int Dir = 0; Dir < 2; ++Dir) {
			for (const auto & i : TerminalMinMaxIndices) {
				if (std::find(MinMaxIndices[Dir].begin(), MinMaxIndices[Dir].end(), i[Dir]) == MinMaxIndices[Dir].end())
					MinMaxIndices[Dir].push_back(i[Dir]);
			}
			/*
			 *	Sort the min and max indices by index number
			 *	(though it may make more sense to sort by value of
			 *	each terminal element's BasinVarNum)
			 */
			std::sort(MinMaxIndices[Dir].begin(), MinMaxIndices[Dir].end());
			NumMinMax[Dir] = MinMaxIndices[Dir].size();
		}
		if (NumMinMax == OldNumMinMax) convegedIter++;
		else convegedIter = 0;
	}
	

	/*
	 *	Now get the integrals of all the variables for the basins found.
	 *	These are stored in a 3d irregular vector whose first dimension is 
	 *	of length 2 (min and max), second dimension is of length the number 
	 *	of mins and maxes respectively, and the third dimension is of length
	 *	number of integration variables.
	 *	Also save lists of nodes and elements for each min/max.
	 */
	BasinIntVals.resize(2);
	BasinNodes.resize(2);
	BasinElems.resize(2);
	SphereElemNums.resize(2);
	SphereNodeNums.resize(2);
	for (int i = 0; i < 2; ++i) {
		BasinIntVals[i] = vector<vector<double> >(MinMaxIndices[i].size(), vector<double>(NumIntVars, 0));
		BasinNodes[i] = vector<vector<vec3> >(MinMaxIndices[i].size());
		BasinElems[i] = vector<vector<vector<int> > >(MinMaxIndices[i].size());
		SphereElemNums[i] = vector<vector<int> >(MinMaxIndices[i].size());
		for (auto & j : SphereElemNums[i]) j.reserve(int(double(NumElems) / double(MinMaxIndices[i].size())));
		SphereNodeNums[i] = vector<vector<int> >(MinMaxIndices[i].size());
		for (auto & j : SphereNodeNums[i]) j.reserve(int(double(NumElems) / double(MinMaxIndices[i].size())));
	}

	for (int Dir = 0; Dir < 2; ++Dir) {
 #pragma omp parallel for schedule(dynamic)
		for (int te = 0; te < MinMaxIndices[Dir].size(); ++te) {
			vector<int> NewNodeNums(XYZListPtr->size(), -1);
			for (int e = 0; e < NumElems; ++e) {
				if (TerminalMinMaxIndices[e][Dir] == MinMaxIndices[Dir][te]) {
					for (int v = 0; v < NumIntVars; ++v) {
						BasinIntVals[Dir][te][v] += ElemIntVals[v][e];
					}
					BasinElems[Dir][te].push_back(vector<int>());
					BasinElems[Dir][te].back().reserve(3);
					for (const int & ei : ElemListPtr->at(e)) {
						if (NewNodeNums[ei] < 0) {
							NewNodeNums[ei] = BasinNodes[Dir][te].size();
							BasinNodes[Dir][te].push_back(XYZListPtr->at(ei));

							if (std::find(SphereNodeNums[Dir][te].begin(), SphereNodeNums[Dir][te].end(), ei) == SphereNodeNums[Dir][te].end())
								SphereNodeNums[Dir][te].push_back(ei);
						}
						BasinElems[Dir][te].back().push_back(NewNodeNums[ei]);
					}
					SphereElemNums[Dir][te].push_back(e);
				}
			}
		}
	}

	return;
}

void FindSphereBasins() {
	/*
	 *	Get the selected Sphere and integration variable
	 */

	int SphereZoneNum, NumIntVars = 0, IntNum = 0;
	string SphereName;
	NumIntVars = TecGUIListGetItemCount(SLSelVar_SLST_T3_1);
	vector<string> IntVarNames(NumIntVars);

	if (NumIntVars > 0 && TecGUIListGetItemCount(SLSelSphere_SLST_T3_1) > 0) {
		IntNum = TecGUIListGetSelectedItem(SLSelVar_SLST_T3_1);
	}

	string BasinDefineVarName;
	if (IntNum > 0) {
		char* cstr = TecGUIListGetString(SLSelVar_SLST_T3_1, IntNum);
		BasinDefineVarName = cstr;
		TecUtilStringDealloc(&cstr);
	}

	vector<int> SphereZoneNums;
	vector<string> SphereZoneNames, NuclearNames;
	for (int z = 1; z <= TecUtilDataSetGetNumZones(); ++z) {
		if (AuxDataZoneHasItem(z, CSMAuxData.GBA.AtomicBasinIntegrationValues)) {
			SphereZoneNums.push_back(z);
			char *cstr;
			TecUtilZoneGetName(z, &cstr);
			SphereZoneNames.push_back(cstr);
			NuclearNames.push_back(AuxDataZoneGetItem(z, CSMAuxData.GBA.SourceNucleusName));
			TecUtilStringDealloc(&cstr);
		}
	}

// 	/*
// 	 *	Right now it's not handled correctly if a user runs this twice on a single 
// 	 *	data set, so delete any preexisting zones that would be created here
// 	 *	to prevent a conflict later.
// 	 */
// 	Set oldZoneSet;
// 	for (int z = 1; z <= TecUtilDataSetGetNumZones(); ++z) {
// 		if (AuxDataZoneHasItem(z, CSMAuxData.GBA.CondensedBasinInfo)) {
// 			string tmpStr;
// 			if (AuxDataZoneGetItem(z, CSMAuxData.GBA.SphereCPName, tmpStr)
// 				&& std::find(SphereZoneNames.begin(), SphereZoneNames.end(), tmpStr) != SphereZoneNames.end()) {
// 				oldZoneSet += z;
// 			}
// 		}
// 	}
// 	if (!oldZoneSet.isEmpty()) {
// 		TecUtilDataSetDeleteZone(oldZoneSet.getRef());
// 	}

	high_resolution_clock::time_point startTime = high_resolution_clock::now();
	if (IntNum > 0 && SphereZoneNums.size() > 0) {
		StatusLaunch("Finding gradient bundles...(please wait)", AddOnID);
		for (int z = 0; z < SphereZoneNums.size(); ++z) {
// 			z = 1;
			SphereName = (NuclearNames[z] != "" ? NuclearNames[z] : SphereZoneNames[z]);
			StatusUpdate(z, SphereZoneNums.size(), "Finding gradient bundles for " + SphereZoneNames[z] + " (" + to_string(z + 1) + " of " + to_string(SphereZoneNums.size()) + " )", AddOnID, startTime);
			SphereZoneNum = SphereZoneNums[z];
			REQUIRE(SphereZoneNum > 0 && TecUtilZoneGetType(SphereZoneNum) == ZoneType_FETriangle);
			FESurface_c Sphere(SphereZoneNum, ZoneNumByName("Full Volume"), { 1,2,3 }, vector<int>(), true);
			/*
				*	Get all the integration values from the Sphere Elements
				*/
			int NumElems = Sphere.GetElemListPtr()->size();
			vector<vector<double> > ElemIntVals(NumIntVars, vector<double>(NumElems));
			for (int i = 0; i < NumIntVars; ++i) {
				int IntVarNum = VarNumByName(TecGUIListGetString(SLSelVar_SLST_T3_1, i + 1));
				REQUIRE(IntVarNum > 0);
				FieldData_pa CCRef = TecUtilDataValueGetReadableCCRef(SphereZoneNum, IntVarNum);
				REQUIRE(VALID_REF(CCRef));
				TecUtilDataValueArrayGetByRef(CCRef, 1, NumElems, ElemIntVals[i].data());
				IntVarNames[i] = TecGUIListGetString(SLSelVar_SLST_T3_1, i + 1);
			}

			/*
				*	Get all the basins.
				*	First prepare the arrays to receive the data.
				*	We receive one array with the min/max node indices on the sphere (2 x number of min/max),
				*	one with the integration values inside each (2 x number of min/max x number of int vars),
				*	one with the nodes of each (2 x number of min/max x number of nodes),
				*	and one with the elements of each (2 x number of min/max x number of elements).
				*/
			vector<vector<int> > MinMaxIndices;
			vector<vector<vector<double> > > BasinIntVals;
			vector<vector<vector<vec3> > > BasinNodes;
			vector<vector<vector<vector<int> > > > BasinElems;
			vector<vector<vector<int> > > SphereElemNums, SphereNodeNums;

			GetSphereBasinIntegrations(Sphere, ElemIntVals, IntNum - 1, MinMaxIndices, BasinIntVals, BasinNodes, BasinElems, SphereNodeNums, SphereElemNums);

			/*
				*	Make a surface from gradient paths seeded at the basin parameter edge midpoints.
				*	Start by collecting the seed points for all the gradient paths.
				*/

				/*
				*	Surfaces of basins that include sphere nodes corresponding to bond paths or
				*	ring surfaces need special treatment.
				*	Collect the sphere node numbers with correspondence to CPs in the volume zone
				*	and the types of CPs.
				*/

			vector<int> SphereConstrainedNodeNums, SphereConstrainedGPNodeNums, SphereConstrainedNodeCPNums;
			vector<string> SphereConstrainedNodeNames, SphereConstrainedNodeNamesTotalCount, SphereConstrainedNodeTypeStrs;
			vector<GradPath_c> SphereConstrainedNodeGPs;
			int SphereNuclearCPNum, NumConstrainedGPs = 0;
			string TmpStr;
			if (AuxDataZoneGetItem(Sphere.GetZoneNum(), CSMAuxData.GBA.SphereConstrainedNodeNums, TmpStr)) {
				SphereConstrainedNodeNums = SplitStringInt(TmpStr, ",");
				for (int & i : SphereConstrainedNodeNums) i--;// moving from base 1 to base 0
			}
			else {
				TecUtilDialogErrMsg("Failed to get constrained node numbers");
			}


			if (AuxDataZoneGetItem(Sphere.GetZoneNum(), CSMAuxData.GBA.SphereCPNum, TmpStr)) {
				SphereNuclearCPNum = stoi(TmpStr);
			}
			else {
				TecUtilDialogErrMsg("Failed to get constrained node numbers");
			}

			if (AuxDataZoneGetItem(Sphere.GetZoneNum(), CSMAuxData.GBA.SphereConstrainedNodeIntersectCPNames, TmpStr)) {
				SphereConstrainedNodeNames = SplitString(TmpStr, ",");
				vector<string> strVec;
				for (const auto & i : SphereConstrainedNodeNames) {
					if (i.length() > 0) {
						strVec.push_back(i);
					}
				}
				SphereConstrainedNodeNames = strVec;
			}
			if (AuxDataZoneGetItem(Sphere.GetZoneNum(), CSMAuxData.GBA.SphereConstrainedNodeIntersectCPTotalOffsetNames, TmpStr)) {
				SphereConstrainedNodeNamesTotalCount = SplitString(TmpStr, ",");
				vector<string> strVec;
				for (int i = 0; i < SphereConstrainedNodeNamesTotalCount.size(); ++i) {
					if (SphereConstrainedNodeNamesTotalCount[i].length() > 0) {
						strVec.push_back(SphereConstrainedNodeNamesTotalCount[i]);
						vector<string> TmpStrVec = SplitString(SphereConstrainedNodeNamesTotalCount[i]);
						if (TmpStrVec.size() > 1) {
							SphereConstrainedNodeTypeStrs.push_back(TmpStrVec[0]);
							SphereConstrainedNodeCPNums.push_back(stoi(TmpStrVec[1]));
							if (SphereConstrainedNodeTypeStrs.back() == CPNameList[1]) {
								NumConstrainedGPs++;
								SphereConstrainedGPNodeNums.push_back(SphereConstrainedNodeNums[i]);
							}
						}
					}
				}
				SphereConstrainedNodeNamesTotalCount = strVec;
			}
			else {
				TecUtilDialogErrMsg("Failed to get constrained node numbers");
			}

			/*
				*	Get all necessary raw pointers for GradPath creation.
				*/
			EntIndex_t CutoffVarNum = VarNumByName(string("Electron Density"));
			int VolZoneNum = ZoneNumByName("Full Volume");
			EntIndex_t GradXYZVarNums[3];
			Boolean_t IsOk = TRUE;
			if (IsOk) {
				vector<string> TmpStrs = {
					"X Density Gradient",
					"Y Density Gradient",
					"Z Density Gradient"
				};
				for (int i = 0; i < 3 && IsOk; ++i) {
					GradXYZVarNums[i] = VarNumByName(TmpStrs[i]);
					IsOk = (GradXYZVarNums[i] > 0);
				}
				if (!IsOk) {
					StatusDrop(AddOnID);
					TecUtilDialogErrMsg("Couldn't find gradient vector variables.");
					TecUtilDataLoadEnd();
					TecUtilLockFinish(AddOnID);
					return;
				}
			}

// 			if (IsOk) {
// 
// 				// Enable system volume zone
// 
// 				Set_pa TmpSet = TecUtilSetAlloc(FALSE);
// 				TecUtilSetAddMember(TmpSet, 1, FALSE);
// 				TecUtilZoneSetActive(TmpSet, AssignOp_PlusEquals);
// 				TecUtilSetDealloc(&TmpSet);
// 
// 			}

			vector<FieldDataPointer_c> GradRawPtrs(3);
			FieldDataPointer_c RhoRawPtr;

			// 		Data load begin tells Tecplot that I don't want it to reorganize memory until I call
			// 		Data load end, so that the raw pointers fetched below remain valid.
			TecUtilDataLoadBegin();

			for (int i = 0; i < 4 && IsOk; ++i) {
				if (i == 0) {
					IsOk = RhoRawPtr.GetReadPtr(VolZoneNum, CutoffVarNum);
				}
				else {
					IsOk = GradRawPtrs[i - 1].GetReadPtr(VolZoneNum, GradXYZVarNums[i - 1]);
				}
			}

			double TermRhoValue;
			vector<vector<FESurface_c> > WedgeSurfaces(2);
			vector<VolExtentIndexWeights_s> VolInfo(omp_get_num_procs());
			GetVolInfo(VolZoneNum, { 1,2,3 }, FALSE, VolInfo[0]);
			for (int i = 1; i < VolInfo.size(); ++i) VolInfo[i] = VolInfo[0];
			TecGUITextFieldGetDouble(TFCutoff_TF_T1_1, &TermRhoValue);
			int NumBasins = 0;
			for (int i = 0; i < 2; ++i) {
				NumBasins += MinMaxIndices[i].size();
				WedgeSurfaces[i].resize(MinMaxIndices[i].size());
			}
			int NumGPPts;
			TecGUITextFieldGetLgIndex(TFSTPts_TF_T1_1, &NumGPPts);

			for (const auto & i : SphereConstrainedNodeCPNums) {
				/*
					*	Now get the bond path half for this node
					*/
				for (int zi = 1; zi <= TecUtilDataSetGetNumZones(); ++zi) {
					if (AuxDataZoneItemMatches(zi, CSMAuxData.CC.GPEndNumStrs[0], to_string(i))
						&& AuxDataZoneItemMatches(zi, CSMAuxData.CC.GPEndNumStrs[1], to_string(SphereNuclearCPNum))
						&& AuxDataZoneItemMatches(zi, CSMAuxData.CC.ZoneSubType, CSMAuxData.CC.ZoneSubTypeBondPathSegment)) {
						SphereConstrainedNodeGPs.push_back(GradPath_c(zi, { 1,2,3,CutoffVarNum }, AddOnID));
						SphereConstrainedNodeGPs.back().Reverse();
						break;
					}
				}
			}

			REQUIRE(NumConstrainedGPs == SphereConstrainedNodeGPs.size());

//  		#ifndef _DEBUG
#pragma omp parallel for schedule(dynamic)
// #endif
			for (int b = 0; b < NumBasins; ++b) {
// 				break;
				int i, j;
				if (b < MinMaxIndices[0].size()) {
					i = 0;
					j = b;
				}
				else {
					i = 1;
					j = b - MinMaxIndices[0].size();
				}

// 				i = 0; j = 2;

				int ThreadNum = omp_get_thread_num();
				vector<vec3> SeedPts;
				vector<vector<int> > BasinPerimeterEdges;
				if (GetSortedParameterEdgeMidpoints(BasinElems[i][j], BasinNodes[i][j], SeedPts, BasinPerimeterEdges)) {
					int NumSeedPts = SeedPts.size();
					/*
						*	Check to see if any bond path-sphere intersections appear in the perimeter.
						*	If they do, need to add them to the surface.
						*/
					std::set<int> ConstrainedNodeNums(SphereConstrainedGPNodeNums.begin(), SphereConstrainedGPNodeNums.end());
					bool ContainsConstrainedNode = false;
					for (int e = 0; e < BasinPerimeterEdges.size() && !ContainsConstrainedNode; ++e) {
						for (const int & ei : BasinPerimeterEdges[e]) {
							if (ConstrainedNodeNums.count(SphereNodeNums[i][j][ei]) > 0) {
								ContainsConstrainedNode = true;
								break;
							}
						}
					}
					vector<GradPath_c> BasinParameterGPs(NumSeedPts);
					vector<GradPath_c*> GPPtrs;
					for (int p = 0; p < NumSeedPts; ++p) {
						BasinParameterGPs[p].SetupGradPath(
							SeedPts[p],
							StreamDir_Both,
							NumGPPts,
							GPType_Classic,
							GPTerminate_AtRhoValue,
							NULL, NULL, NULL,
							&TermRhoValue,
							VolInfo[ThreadNum],
							vector<FieldDataPointer_c>(),
							GradRawPtrs,
							RhoRawPtr);
						GPPtrs.push_back(&BasinParameterGPs[p]);
					}
#pragma omp parallel for
					for (int p = 0; p < NumSeedPts; ++p) {
						BasinParameterGPs[p].Seed();
						BasinParameterGPs[p].Reverse();
					}
					if (ContainsConstrainedNode) {
						/*
							*	There are one or more constrained nodes that need to be added
							*	to the surface.
							*/
						auto * XYZListPtr = Sphere.GetXYZListPtr();
						vector<GradPath_c*> TmpGPPtrs;
						for (int p = 0; p < NumSeedPts; ++p) {
							TmpGPPtrs.push_back(GPPtrs[p]);
							for (int n = 0; n < SphereConstrainedGPNodeNums.size(); ++n) {
								if (SphereConstrainedGPNodeNums[n] == SphereNodeNums[i][j][BasinPerimeterEdges[p][1]]) {
									/*
										*	Now need to take the bond path segment for the constrained node
										*	and add new gradient paths in the interatomic surface to be stitched with the
										*	neighboring edge midpoint gradient paths.
										*	First, get the closest points on the neighboring edge midpoint gradient
										*	paths to the bond path and use that information to find the
										*	angle between them in the interatomic plane and then the
										*	seed points for the new gradient paths in the interatomic plane.
										*/
									vector<vec3> TmpSeedPts;
									for (int gpi = 0; gpi < 2; ++gpi) {
										vec3 ClosestPt = GPPtrs[(p + gpi) % GPPtrs.size()]->ClosestPoint(SphereConstrainedNodeGPs[n][-1]);

										/*
											*	Now place the closest point in the interatomic plane
											*/
										vec3 v1 = normalise(SphereConstrainedNodeGPs[n][-2] - SphereConstrainedNodeGPs[n][-1]),
											v3 = ClosestPt - SphereConstrainedNodeGPs[n][-1],
											v2 = normalise(cross(v1, v3));
										double tmpLen = norm(v3);
										vec3 tmpv3 = normalise(cross(v1, v2));
										if (dot(tmpv3, v3) < 0)
											v3 = -tmpv3;
										else
											v3 = tmpv3;

										TmpSeedPts.push_back(v3);
									}

									/*
										*	Now we've got the two vectors to use to seed the neighboring gradient paths
										*/
									for (const auto & seedPt : TmpSeedPts) {
										TmpGPPtrs.push_back(new GradPath_c);
										TmpGPPtrs.back()->SetupGradPath(SphereConstrainedNodeGPs[n][-1] + (seedPt * 0.5),
											StreamDir_Reverse,
											NumGPPts,
											GPType_Classic,
											GPTerminate_AtRhoValue,
											NULL, NULL, NULL,
											&TermRhoValue,
											VolInfo[ThreadNum],
											vector<FieldDataPointer_c>(),
											GradRawPtrs,
											RhoRawPtr);
										TmpGPPtrs.back()->Seed();
										*TmpGPPtrs.back() = SphereConstrainedNodeGPs[n] + *TmpGPPtrs.back();
									}

								}
							}
						}
						GPPtrs = TmpGPPtrs;
					}
					if (GPPtrs.size() > 2) {
						WedgeSurfaces[i][j].MakeFromGPs(GPPtrs, true);
					}
				}
			}

			/*
				*	Now save the basins and wedge surfaces as zones.
				*/

			/*
			 *	Place nodes of basin a bit further away from the nuclear
			 *	CP than the sphere so that the resulting surfaces aren't
			 *	coincident with the sphere.
			 */
			string CPZoneNumStr, CPIndStr, RadStr;
			vec3 CPPos;
			double SphereRadius = -1;
			if (AuxDataZoneGetItem(SphereZoneNum, CSMAuxData.GBA.SourceZoneNum, CPZoneNumStr)
				&& AuxDataZoneGetItem(SphereZoneNum, CSMAuxData.GBA.SphereCPNum, CPIndStr)
				&& AuxDataZoneGetItem(SphereZoneNum, CSMAuxData.GBA.SphereRadius, RadStr)){
				vector<int> XYZVarNums = { 1,2,3 };
				int CPZoneNum = stoi(CPZoneNumStr);
				int CPInd = stoi(CPIndStr);
				SphereRadius = stod(RadStr);
				for (int i = 0; i < 3; ++i) {
					CPPos[i] = TecUtilDataValueGetByZoneVar(CPZoneNum, XYZVarNums[i], CPInd);
				}
			}



			Set ZoneSet;
			vector<string> MinMax = { "Min","Max" };
			for (int i = 0; i < 2; ++i) {
				for (int b = 0; b < MinMaxIndices[i].size(); ++b) {
					if (SphereRadius > 0) {
						// #pragma omp parallel for
						for (int vi = 0; vi < BasinNodes[i][b].size(); ++vi) {
							BasinNodes[i][b][vi] += normalise(BasinNodes[i][b][vi] - CPPos) * 0.001;
						}
					}
					FESurface_c Basin(BasinNodes[i][b], BasinElems[i][b]);
					REQUIRE(Basin.IsMade());
					int SphereConstrainedNodeNum = -1;
					for (int n = 0; n < SphereConstrainedNodeNums.size(); ++n) {
						if (((i == 0 && SphereConstrainedNodeTypeStrs[n] == RankStrs[3])
							|| (i == 1 && SphereConstrainedNodeTypeStrs[n] == RankStrs[1]))
							&& std::find(SphereNodeNums[i][b].begin(), SphereNodeNums[i][b].end(), SphereConstrainedNodeNums[n]) != SphereNodeNums[i][b].end())
						{
							SphereConstrainedNodeNum = n;
							break;
						}
					}

					int ZoneNum = Basin.SaveAsTriFEZone({ 1,2,3 }, SphereName + ": " + MinMax[i] + " basin (node " + to_string(MinMaxIndices[i][b]) + ") " + BasinDefineVarName);
					if (ZoneNum > 0) {
						if (i == 1) ZoneSet += ZoneNum;
						/*
							*	Set Aux data
							*/
						AuxDataZoneSetItem(ZoneNum, CSMAuxData.GBA.SourceZoneNum, to_string(SphereZoneNum));
						AuxDataZoneSetItem(ZoneNum, CSMAuxData.GBA.SphereCPName, SphereZoneNames[z]);
						AuxDataZoneSetItem(ZoneNum, CSMAuxData.GBA.SourceNucleusName, NuclearNames[z]);
						AuxDataZoneSetItem(ZoneNum, CSMAuxData.GBA.CondensedBasinDefiningVariable, BasinDefineVarName);
						AuxDataZoneSetItem(ZoneNum, CSMAuxData.GBA.IntVarNames, VectorToString(IntVarNames, ","));
						AuxDataZoneSetItem(ZoneNum, CSMAuxData.GBA.IntVarVals, VectorToString(BasinIntVals[i][b], ","));
						AuxDataZoneSetItem(ZoneNum, CSMAuxData.GBA.ZoneType, (i == 0 ? CSMAuxData.GBA.ZoneTypeCondensedRepulsiveBasin : CSMAuxData.GBA.ZoneTypeCondensedAttractiveBasin));
						AuxDataZoneSetItem(ZoneNum, CSMAuxData.GBA.CondensedBasinSphereElements, VectorToString(SphereElemNums[i][b], ","));
						AuxDataZoneSetItem(ZoneNum, CSMAuxData.GBA.CondensedBasinSphereNodes, VectorToString(SphereNodeNums[i][b], ","));
						AuxDataZoneSetItem(ZoneNum, CSMAuxData.GBA.CondensedBasinCentralElementIndex, to_string(MinMaxIndices[i][b]));

						if (SphereConstrainedNodeNum >= 0) {
							AuxDataZoneSetItem(ZoneNum, CSMAuxData.GBA.CondensedBasinName, SphereConstrainedNodeNames[SphereConstrainedNodeNum]);
							AuxDataZoneSetItem(ZoneNum, CSMAuxData.GBA.CondensedBasinIsSpecialGradientBundle, "true");
						}
						else {
							AuxDataZoneSetItem(ZoneNum, CSMAuxData.GBA.CondensedBasinName, SphereName + ": " + MinMax[i] + " basin (node " + to_string(MinMaxIndices[i][b]) + ")");
						}
						AuxDataZoneSetItem(ZoneNum, CSMAuxData.GBA.CondensedBasinInfo, SphereName + ": " + MinMax[i] + " (node " + to_string(MinMaxIndices[i][b]) + ")");

						Set TmpSet(ZoneNum);
						TecUtilZoneSetMesh(SV_SHOW, TmpSet.getRef(), 0.0, TRUE);
						TecUtilZoneSetMesh(SV_COLOR, TmpSet.getRef(), 0.0, ColorIndex_t((b % 7) + 1));
						TecUtilZoneSetMesh(SV_LINEPATTERN, TmpSet.getRef(), 0.0, (i == 0 ? LinePattern_Dotted : LinePattern_Solid));
						TecUtilZoneSetShade(SV_SHOW, TmpSet.getRef(), 0.0, FALSE);
						TecUtilZoneSetShade(SV_COLOR, TmpSet.getRef(), 0.0, ColorIndex_t((b % 7) + 1));
						TecUtilZoneSetContour(SV_SHOW, TmpSet.getRef(), 0.0, FALSE);
						TecUtilZoneSetEdgeLayer(SV_SHOW, TmpSet.getRef(), 0.0, FALSE);


					}

  					if (WedgeSurfaces[i][b].IsMade()) {
  						ZoneNum = WedgeSurfaces[i][b].SaveAsTriFEZone({ 1,2,3 }, SphereName + ": " + MinMax[i] + " wedge (node " + to_string(MinMaxIndices[i][b]) + ") " + BasinDefineVarName);
  						if (ZoneNum > 0) {
  							if (i == 1) ZoneSet += ZoneNum;
  							/*
  								*	Set Aux data
  								*/
  							AuxDataZoneSetItem(ZoneNum, CSMAuxData.GBA.SourceZoneNum, to_string(SphereZoneNum));
  							AuxDataZoneSetItem(ZoneNum, CSMAuxData.GBA.SphereCPName, SphereZoneNames[z]);
  							AuxDataZoneSetItem(ZoneNum, CSMAuxData.GBA.SourceNucleusName, NuclearNames[z]);
  							AuxDataZoneSetItem(ZoneNum, CSMAuxData.GBA.CondensedBasinDefiningVariable, BasinDefineVarName);
  							AuxDataZoneSetItem(ZoneNum, CSMAuxData.GBA.IntVarNames, VectorToString(IntVarNames, ","));
  							AuxDataZoneSetItem(ZoneNum, CSMAuxData.GBA.IntVarVals, VectorToString(BasinIntVals[i][b], ","));
  							AuxDataZoneSetItem(ZoneNum, CSMAuxData.GBA.ZoneType, (i == 0 ? CSMAuxData.GBA.ZoneTypeCondensedRepulsiveBasinWedge : CSMAuxData.GBA.ZoneTypeCondensedAttractiveBasinWedge));
  							if (SphereConstrainedNodeNum >= 0) {
  								AuxDataZoneSetItem(ZoneNum, CSMAuxData.GBA.CondensedBasinName, SphereConstrainedNodeNames[SphereConstrainedNodeNum]);
								AuxDataZoneSetItem(ZoneNum, CSMAuxData.GBA.CondensedBasinIsSpecialGradientBundle, "true");
  							}
  							else {
  								AuxDataZoneSetItem(ZoneNum, CSMAuxData.GBA.CondensedBasinName, SphereName + ": " + MinMax[i] + " wedge (node " + to_string(MinMaxIndices[i][b]) + ")");
  							}
  							AuxDataZoneSetItem(ZoneNum, CSMAuxData.GBA.CondensedBasinInfo, SphereName + ": " + MinMax[i] + " (node " + to_string(MinMaxIndices[i][b]) + ")");
  
  							Set TmpSet(ZoneNum);
  							TecUtilZoneSetMesh(SV_SHOW, TmpSet.getRef(), 0.0, FALSE);
  							TecUtilZoneSetMesh(SV_COLOR, TmpSet.getRef(), 0.0, ColorIndex_t((b % 7) + 1));
  							TecUtilZoneSetMesh(SV_LINEPATTERN, TmpSet.getRef(), 0.0, (i == 0 ? LinePattern_Dotted : LinePattern_Solid));
  							TecUtilZoneSetShade(SV_SHOW, TmpSet.getRef(), 0.0, TRUE);
  							TecUtilZoneSetShade(SV_COLOR, TmpSet.getRef(), 0.0, ColorIndex_t((b % 7) + 1));
  							TecUtilZoneSetContour(SV_SHOW, TmpSet.getRef(), 0.0, FALSE);
  							TecUtilZoneSetEdgeLayer(SV_SHOW, TmpSet.getRef(), 0.0, TRUE);
  						}
  					}
				}
// 				break;
			}

			TecUtilDataLoadEnd();

// 			TecUtilZoneSetActive(ZoneSet.getRef(), AssignOp_PlusEquals);

// 			break;
		}
		StatusDrop(AddOnID);
	}

	GBAResultViewerPopulateGBs();
}