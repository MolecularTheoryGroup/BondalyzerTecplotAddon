
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
#include <queue>

// for profiling
#include <chrono>
#include <ctime>
#include <ratio>

#include <omp.h>

#include "Set.h"
#include "ArgList.h"
#include "StyleValue.h"

#include "CSM_DATA_SET_INFO.h"
#include "CSM_DATA_TYPES.h"
#include "meshgen2d_sphere.h"
#include "VIEWRESULTS.h"
#include "CSM_FE_VOLUME.h"
#include "CSM_GRAD_PATH.h"
#include "CSM_CRIT_POINTS.h"
#include "CSM_GBAGUI.h"
#include "CSM_GUI.h"
#include "ZONEVARINFO.h"
#include "CSM_GEOMETRY.h"
#include "CSM_CALC_VARS.h"

#include "updateSphericalTriangulation.h"

#include "GBAENGINE.h"

#include "spdlog/spdlog.h"
#include "spdlog/sinks/rotating_file_sink.h" // support for rotating file logging

#include <armadillo>
using namespace arma;
using namespace tecplot::toolbox;


//  #define SUPERDEBUG
// #define DEBUG_SAVEZONES

#ifdef _DEBUG
#endif

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

double const DivergentGPMinTerminalDotProduct = cos(150. / DEGPERRADIANS);

vector<vector<int> > const ElemSubdivisionNewElemIndices = { {0,3,5},{1,4,3},{2,5,4},{3,4,5} }; //counterclockwise
// vector<vector<int> > const ElemSubdivisionNewElemIndices = { {0,5,3},{1,3,4},{2,4,5},{3,5,4} }; //clockwise

string const LogName = "GBA_log";
int const LogSizeMB = 5;
int const LogNumLogs = 10;

vector<string> const NETypeStrs = { "C", "R", "B", "RB" };

enum GBATriangulatedSphereNodeElemType_e
{
	NETypeInvalid = -1,
	NETypeC, // goes to cage cp, system boundary, or rho cutoff
	NETypeR, // goes along a ring surface to a ring cp, then to a cage cp or rho cutoff
	NETypeB, // goes along bond path
	NETypeRB // B type on same elem as R type
};

enum RadMode{
	ABSOLUTERADIUS = 1,
	MINCPDISTRATIO
};

string GetNucleusNameForCP(vec3 const & CPPt, vector<int> const XYZVarNums = { 1,2,3 }, double const & DistCutoff = 0.1)
{
	REQUIRE(XYZVarNums.size() == 3);
	for (auto const & i : XYZVarNums) REQUIRE(i > 0);
	REQUIRE(DistCutoff > 0);

	double MinDist = DBL_MAX, TmpDist;
	int MinZoneNum = -1, MinPtNum = -1;

	TecUtilDataLoadBegin();
	int NumZones = TecUtilDataSetGetNumZones();
	for (int z = 1; z <= NumZones; ++z) {
		if (AuxDataZoneItemMatches(z, CSMAuxData.DL.ZoneType, CSMAuxData.DL.ZoneTypeNuclearPositions)) {
			FieldVecPointer_c XYZPtr;
			if (XYZPtr.InitializeReadPtr(z, XYZVarNums)) {
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
			outStr = string(cstr) + " " + to_string(MinPtNum + 1);
		}
		else{
			TecUtilDialogErrMsg("Failed to get zone name when finding nucleus name");
		}
		TecUtilStringDealloc(&cstr);
	}
	return outStr;
}

void NewMainFunction() {
	SYSTEM_INFO sysinfo;
	GetSystemInfo(&sysinfo);

	int numCPU = sysinfo.dwNumberOfProcessors;
	omp_set_num_threads(numCPU);

	Boolean_t IsOk = TRUE;

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
	if (IsOk) {
		if (VolZoneNum <= 0) {
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

	for (int i = 0; i < NumVars; ++i) {
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

	for (int i = 1; i <= TecUtilDataSetGetNumZones(); ++i) {
		if (AuxDataZoneHasItem(i, CSMAuxData.GBA.ZoneType)) {
			for (int j = 3; j < NumVars; ++j) {
				VarDataTypes[j] = TecUtilDataValueGetType(i, j + 1);
				VarLocations[j] = TecUtilDataValueGetLocation(i, j + 1);
			}
			break;
		}
	}

	int GroupNum = 0;

	LgIndex_t * CPNums = nullptr;
	LgIndex_t NumSelectedCPs = 1;
	TecGUIListGetSelectedItems(MLSelCPs_MLST_T1_1, &CPNums, &NumSelectedCPs);

	Set oldSphereZonesToDelete;
	bool deleteOldSphereZones = false, deleteOldSphereZonesAsked = false;


	/*
 *	Delete any existing sphere for the current CP if so chosen by the user
 */
	for (int SelectCPNum = 0; SelectCPNum < NumSelectedCPs && IsOk; ++SelectCPNum)
	{
		int CPType = -3;
		string CPString, CPName;
		int junkZoneNum;
		LgIndex_t CPNum = -1;
		bool IsCP;
		vector<int> NumCPs;
		GetCoordsFromListItem(CPNums[SelectCPNum], MLSelCPs_MLST_T1_1, &CPString, &CPName, &CPNum, &junkZoneNum, &IsCP, &NumCPs);
		for (int z = 1; z <= TecUtilDataSetGetNumZones(); ++z) {
			if (AuxDataZoneItemMatches(z, CSMAuxData.GBA.SphereCPName, CPString) && (!deleteOldSphereZonesAsked || deleteOldSphereZones)) {
				string tmpStr = "An existing sphere zone was found for " + CPString + ". Would you like to erase it (and all other conflicting sphere zones) before continuing?";
				if (!deleteOldSphereZonesAsked) {
					deleteOldSphereZones = TecUtilDialogMessageBox(tmpStr.c_str(), MessageBox_YesNo);
					deleteOldSphereZonesAsked = true;
				}
				if (deleteOldSphereZones) {
					oldSphereZonesToDelete += z;
					for (int zi = z + 1; zi <= TecUtilDataSetGetNumZones(); zi++) {
						if (AuxDataZoneItemMatches(zi, CSMAuxData.GBA.SphereCPName, CPString))
							oldSphereZonesToDelete += zi;
					}
				}
				break;
			}
		}
	}


	/*
	 *	Get job parameters from dialog
	 */
	double CutoffVal = 0.001;
	TecGUITextFieldGetDouble(TFCutoff_TF_T1_1, &CutoffVal);
	EntIndex_t RhoVarNum = VarNumByName("Electron Density");

	vector<int> GradVarNums, HessVarNums;
	GradVarNums.push_back(VarNumByName("X Density", true));
	if (GradVarNums.back() > 0) {
		for (int i = 1; i < 3; ++i)
			GradVarNums.push_back(GradVarNums[0] + i);
	}
	else {
		GradVarNums.clear();
	}

	HessVarNums.push_back(VarNumByName("XX Density", true));
	if (HessVarNums.back() > 0) {
		for (int i = 1; i < 6; ++i)
			HessVarNums.push_back(HessVarNums[0] + i);
	}
	else
		HessVarNums.clear();

	// Grad and Hess needed, so calculate them if they werent already present
	bool DeleteGradVars = false, DeleteHessVars = false;
	if (GradVarNums.empty()){// || HessVarNums.empty()){
		CalcVarsOptions_s opt;
		opt.AddOnID = AddOnID;
		opt.CalcForAllZones = FALSE;
		opt.CalcZoneNum = VolZoneNum;
		opt.RhoVarNum = RhoVarNum;
		opt.HasGrad = (!GradVarNums.empty());
		if (opt.HasGrad) {
			opt.GradVarNums = GradVarNums;
		}
		else{
			DeleteGradVars = true;
			opt.CalcVarList = { CalcGradientVectors };
		}
// 		DeleteHessVars = true;
// 		opt.CalcVarList.push_back(CalcHessian);
		opt.IsPeriodic = FALSE;

		CalcVars(opt);

		if (GradVarNums.empty()){
			GradVarNums.push_back(VarNumByName("X Density", true));
			for (int i = 1; i < 3; ++i)
				GradVarNums.push_back(GradVarNums[0] + i);

			if (GradVarNums.empty())
				return;
		}

// 		if (HessVarNums.empty()){
// 			HessVarNums.push_back(VarNumByName("XX Density", true));
// 			for (int i = 1; i < 6; ++i)
// 				HessVarNums.push_back(HessVarNums[0] + i);
// 
// 			if (HessVarNums.empty())
// 				return;
// 		}
	}

	Boolean_t RadialSphereApprx = TecGUIToggleGet(TGLNoSphInt_TOG_T1_1);

	double EdgeGPCheckDistance = GBADefaultEdgeGPSpacing;
	EdgeGPCheckDistance = (double)TecGUIScaleGetValue(SCEGPDist_SC_T1_1) / 10.;
	double RingBondEdgeGPCheckDistanceFactor = 0.6;
	int MaxEdgeGPs = 300;
	bool MinEdgeGPs = false;

	bool TestRun = TecGUIToggleGet(TGLSphTest_TOG_T1_1);

#ifdef _DEBUG
// 	MaxEdgeGPs = -1;
#endif


	/*
	*	Set sphere radius and mesh refinement level
	*/

	int NumGBsPerElectron = GBADefaultGBPerE;
	// 	TecGUITextFieldGetLgIndex(TFGBPerE_TF_T1_1, &NumGBsPerElectron);

		// Set number of gradient bundles around bond path intersection nodes
	int NumAngularGBsAroundBPIntersectionNodes = GBADefaultBPAngularGBs;
	// 	TecGUITextFieldGetLgIndex(SCBPGBs_SC_T1_1, &NumAngularGBsAroundBPIntersectionNodes);
	NumAngularGBsAroundBPIntersectionNodes = TecGUIScaleGetValue(SCBPGBs_SC_T1_1);
	NumAngularGBsAroundBPIntersectionNodes = MAX(1, NumAngularGBsAroundBPIntersectionNodes);
	double BPNodeCheckAngle = 2.0 * PI / (double)NumAngularGBsAroundBPIntersectionNodes;

	double TargetDensityPerGB = 1.0 / (double)NumGBsPerElectron;

	int MaxGBSubdivisionLevel = GBADefaultMaxGBSubdivisionLevel;
	// 	TecGUITextFieldGetLgIndex(TFGBMaxSD_TF_T1_1, &MaxGBSubdivisionLevel);

	int NumberOfPreBondPathElemSubdivision = GBAMaxSubdivisionLevel;
// 	// 	TecGUITextFieldGetLgIndex(SCBPGBInit_SC_T1_1, &NumberOfPreBondPathElemSubdivision);
// 	NumberOfPreBondPathElemSubdivision = TecGUIScaleGetValue(SCBPGBInit_SC_T1_1);


	GPResampleMethod_e ResampleMethod = GPResampleMethod_Linear;

	TecUtilDataLoadBegin();

	FieldDataPointer_c RhoPtr;
	vector<FieldDataPointer_c> GradPtrs, HessPtrs;
// 	CSMGuiUnlock();
	GetReadPtrsForZone(VolZoneNum, RhoVarNum, GradVarNums, HessVarNums, RhoPtr, GradPtrs, HessPtrs);
// 	CSMGuiLock();

	int CPZoneNum = ZoneNumByName(CSMZoneName.CriticalPoints);

	EntIndex_t CPTypeVarNum;
	CPTypeVarNum = VarNumByName(CSMVarName.CritPointType);

// 	CritPoints_c CPs(CPZoneNum, XYZVarNums, CPTypeVarNum, RhoVarNum);
	VolExtentIndexWeights_s VolInfo;
	GetVolInfo(VolZoneNum, XYZVarNums, FALSE, VolInfo);

	vector<VolExtentIndexWeights_s> ThVolInfo(numCPU, VolInfo);

	MultiRootParams_s MR;
	MR.VolInfo = &VolInfo;
	MR.IsPeriodic = FALSE;
	MR.HasGrad = (GradPtrs.size() == 3);
	MR.HasHess = (HessPtrs.size() == 6);
	MR.RhoPtr = &RhoPtr;
	MR.GradPtrs = &GradPtrs;
	MR.HessPtrs = &HessPtrs;
	MR.BasisVectors = &VolInfo.BasisVectors;

	CritPoints_c CPs(CPZoneNum, XYZVarNums, CPTypeVarNum, RhoVarNum, &MR);
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





	RadMode RadiusMode = MINCPDISTRATIO;
	RadiusMode = (RadMode)TecGUIRadioBoxGetToggle(RBRadMode_RADIO_T1_1);
	double UserRadius = 0.25;
	TecGUITextFieldGetDouble(TFRad_TF_T1_1, &UserRadius);
	int Level = TecGUIScaleGetValue(SCMinGBs_SC_T1_1);
	int MaxSubdivisions = TecGUIScaleGetValue(SCBPGBInit_SC_T1_1);
	int SubdivisionTightness = TecGUIScaleGetValue(SCSDtight_SC_T1_1);

	LgIndex_t NumGPPoints = 100;
	TecGUITextFieldGetLgIndex(TFSTPts_TF_T1_1, &NumGPPoints);

	double SphereJiggleMeshMaxMovedNodeDistTol = 1e-8;


	bool SaveGPs = TecGUIToggleGet(TGLsGP_TOG_T1_1);
	bool SaveGBs = TecGUIToggleGet(TGLsGB_TOG_T1_1);


	EntIndex_t OldNumZones = TecUtilDataSetGetNumZones();

	/*
		*	information for integration
		*/
	vector<string> AtomNameList;
	vector<string> IntVarNameList, BaseIntVarNameList;
	vector<int> IntVarNumList;
	int IntResolution = TecGUIScaleGetValue(SCPrecise_SC_T1_1);
	vector<FieldDataPointer_c> IntVarPtrs;

	int DensityVarNumInIntList;

	
	AtomNameList = ListGetSelectedStrings(MLSelCPs_MLST_T1_1);
	IntVarNameList = ListGetSelectedStrings(MLSelVars_MLST_T1_1);
	BaseIntVarNameList = IntVarNameList;
	IntVarNumList = ListGetSelectedItemNums(MLSelVars_MLST_T1_1);
	for (auto & i : IntVarNumList)
		i += 3;
	for (int i = 0; i < IntVarNumList.size(); ++i){
		if (IntVarNumList[i] == RhoVarNum){
			DensityVarNumInIntList = i;
			break;
		}
	}

	IntVarPtrs.resize(IntVarNumList.size());

	if (!IsOk) {
		TecUtilDialogErrMsg("Couldn't get XYZ vars."); 
		TecUtilDataLoadEnd();
		return;
	}


	StatusLaunch("GBA: Working", AddOnID, FALSE);


	EntIndex_t CPNuclearZoneNum = ZoneNumByName(CSMZoneName.CPType[0]);

	

	vector<int> NumCPs;
	for (int SelectCPNum = 0; SelectCPNum < NumSelectedCPs && IsOk; ++SelectCPNum)
	{
		LgIndex_t CPNum = -1;
		int CPType = -3;

		int SphereZoneNum = -1;

		double Radius = 0.25;

		vector<vector<double> > IntVals, SphereIntVals;
		vector<double> SphereTotalIntVals;

		stringstream ProgressStr1;
		string CPString, CPName, CPTypeStr;
		bool IsCP;

		int RealCPNum = 0;
		LgIndex_t CPIJKMax[3];

		vec3 CPPos;
		string NucleusName;
		if (IsOk)
		{
			int junkZoneNum;
			CPPos = GetCoordsFromListItem(CPNums[SelectCPNum], MLSelCPs_MLST_T1_1, &CPString, &CPName, &CPNum, &junkZoneNum, &IsCP, &NumCPs);
			NucleusName = GetNucleusNameForCP(CPPos);
			CPType = CPType_Nuclear;
			CPNum = NuclearNameToCPNum[CPString];
			CPPos = CPs.GetXYZ(CPNum);
			CPNum++;
			CPTypeStr = CPNameList[0];
// 			if (IsCP)
// 			{
// 			TecUtilZoneGetIJK(CPZoneNum, &CPIJKMax[0], &CPIJKMax[1], &CPIJKMax[2]);
// 			IsOk = (CPIJKMax[0] > 1 && CPIJKMax[1] == 1 && CPIJKMax[2] == 1);
// 			if (!IsOk) {
// 				TecUtilDialogErrMsg("Critical point zone dimensions incorrect.");
// 				TecUtilDataLoadEnd();
// 				return;
// 			}

			/*
			*	Get XYZ readable pointers to the crit points zone to get the
			*	position of the fifth atom (working with cubane, so this is the first carbon).
			*/

// 			for (int i = 0; i < 4; ++i) {
// 				if (CPName == RankStrs[i] || i == 0 && CPName == "Atom") {
// 					int RealCPNum = 0;
// 					for (int j = 0; j < i; ++j) {
// 						RealCPNum += NumCPs[j];
// 					}
// 					CPNum = RealCPNum + CPNum;
// 					break;
// 				}
// 			}


// 			CPType = (int)TecUtilDataValueGetByZoneVar(CPZoneNum, CPTypeVarNum, CPNum);
// 			for (int i = 0; i < CPTypeList.size(); ++i) {
// 				if (CPType == CPTypeList[i]) {
// 					CPTypeStr = CPNameList[i];
// 				}
// 			}	
// 			}


			ProgressStr1 << "Processing " << (NucleusName.empty() ? CPString : NucleusName) << ": ";
			if (NumSelectedCPs > 1) {
				ProgressStr1 << "(CP " << SelectCPNum + 1 << " of " << NumSelectedCPs << "): ";
			}

			// 			TecUtilStringDealloc(&TempCStr);
		}
		else {
			TecUtilDialogErrMsg("Didn't get CP number correctly");
		}

		string TmpString;

		if (IsOk) {
			TmpString = ProgressStr1.str() + "Generating mesh";
			StatusUpdate(1, 1000, TmpString, AddOnID);
		}

		// Get closest CP if necessary
		if (RadiusMode == MINCPDISTRATIO) {
			double MinDistSqr = DBL_MAX;
			for (int i = CPTypeNum_Bond; i <= CPTypeNum_Cage; ++i) {
				for (int j = 0; j < CPs.NumCPs(i); ++j) {
					MinDistSqr = MIN(MinDistSqr, DistSqr(CPPos, CPs.GetXYZ(i, j)));
				}
			}
			Radius = sqrt(MinDistSqr) * UserRadius;
		}
		else
			Radius = UserRadius;


		double GPNCPTermRadius = MAX(0.001, Radius * 0.01);
		double GPTermRadius = MIN(0.01, GPNCPTermRadius);
		double RTypeGPTermRadiusFactor = 1.0;
		double RTypeTermRadius = GPTermRadius * RTypeGPTermRadiusFactor;


		/*
		 *	Get sphere intersecting bond paths
		 */

		vector<vec3> IntersectingBondPathPoints;
		vector<GradPath_c> IntersectingBondPaths;

		vector<int> XYZRhoVarNums = XYZVarNums;
		XYZRhoVarNums.push_back(RhoVarNum);

		for (int z = 1; z <= TecUtilDataSetGetNumZones(); ++z) {
			if (AuxDataZoneItemMatches(z, CSMAuxData.CC.ZoneSubType, CSMAuxData.CC.ZoneSubTypeBondPathSegment)) {
				vec3 IntPoint;
				GradPath_c TmpGP(z, XYZRhoVarNums, AddOnID);

				// if approximating the sphere interior
				if (RadialSphereApprx && TmpGP.Trim(CPPos, Radius)) {
					IntersectingBondPaths.push_back(TmpGP);
					IntersectingBondPathPoints.push_back(TmpGP.XYZAt(-1));
				}
				// if including the sphere interior
				else if (!RadialSphereApprx && TmpGP.GetSphereIntersectionPoint(CPPos, Radius, IntPoint) >= 0) {
					IntersectingBondPathPoints.push_back(IntPoint);
					IntersectingBondPaths.push_back(TmpGP);
				}

			}
		}

		/*
		 *	Get optimized spherical mesh
		 */

		vector<int> MovedPointNums;
		vector<point> IntersectionPoints;
		for (auto const & p1 : IntersectingBondPathPoints){
			auto p = p1 - CPPos;
			IntersectionPoints.emplace_back(p[0], p[1], p[2]);
		}

		int NumNodes, NumElems, NumEdges;
		point * meshP;
		triangle * meshT;

		int MaxRefine = 2;
		int RefineNum = 0;
		MeshStatus_e MeshStatus = FAIL_INVALIDCONSTRAINT;
		while (MeshStatus != 0) {
			int ** meshE;

			MeshStatus = meshgen2D_sphere(Radius, Level, IntersectionPoints, MovedPointNums, meshP, meshT, meshE, NumNodes, NumElems, NumEdges);


			if (MeshStatus == FAIL_NEEDREFINEMENT) {
				//break;
				if (RefineNum >= MaxRefine) {
					break;
				}
				else {
					Level++;
					RefineNum++;
				}
			}
			else if (MeshStatus == FAIL_INVALIDCONSTRAINT) {
				TecUtilDialogErrMsg("Constrained point too far from mesh");
				TecUtilDataLoadEnd();
				StatusDrop(AddOnID);
				return;
			}
			else if (MeshStatus == SUCCESS_NOTCONVERGED) {
				// 				TecUtilDialogMessageBox("Mesh did not converge.", MessageBoxType_Warning);
				break;
			}

		}

		IsOk = (MeshStatus == SUCCESS || MeshStatus == SUCCESS_NOTCONVERGED && NumNodes > 0 && NumElems > 0);
		if (!IsOk) {
			TecUtilDialogErrMsg("Failed to make mesh."); 
			TecUtilDataLoadEnd();
			StatusDrop(AddOnID);
			return;
		}

		tpcsm::Vec3 SphereCenterVec(CPPos[0], CPPos[1], CPPos[2]);

		vector<tpcsm::Vec3> PreNodes(NumNodes);
		vector<TriNodes> PreElems(NumElems);

		for (int i = 0; i < NumNodes; ++i)
			PreNodes[i] = tpcsm::Vec3(meshP[i][0], meshP[i][1], meshP[i][2]) + SphereCenterVec;
		for (int i = 0; i < NumElems; ++i)
			PreElems[i] = TriNodes(meshT[i][0], meshT[i][1], meshT[i][2]);

#ifdef DEBUG_SAVEZONES
		{
			vector<vec3> SphereNodes;
			vector<vector<int> > SphereElems;
			SphereNodes.reserve(PreNodes.size());
			for (auto const & i : PreNodes)
				SphereNodes.push_back({ i.x(), i.y(), i.z() });
			NumNodes = SphereNodes.size();
			SphereElems.reserve(PreElems.size());
			for (auto const & i : PreElems)
				SphereElems.push_back({ i[0], i[1], i[2] });
			FESurface_c TmpSphere(SphereNodes, SphereElems);
			TmpSphere.SaveAsTriFEZone({ 1,2,3 }, "Initial sphere mesh");
			TecUtilZoneSetActive(Set(TecUtilDataSetGetNumZones()).getRef(), AssignOp_MinusEquals);
		}
#endif

		/*
		 *	Get average sphere triangulation edge length so that we can
		 *	resample intersection paths at roughly the same spacing.
		 *	I'll do this by looping over triangles, which will double count every
		 *	edge so the average should still be valid.
		 */
		std::set<Edge> EdgeList;
		for (auto const & t : PreElems){
			for (int ci = 0; ci < 3; ++ci){
				EdgeList.insert(MakeEdge(t[ci], t[(ci + 1) % 3]));
			}
		}
		double SphereEdgeLenMean = 0.0;
// 		int SphereNumEdges = 0;

		for (auto const & e : EdgeList){
			SphereEdgeLenMean += (PreNodes[e.first] - PreNodes[e.second]).getNorm();
		}
		SphereEdgeLenMean /= (double)EdgeList.size();

// 		for (auto const & t : PreElems) {
// 			for (int i = 0; i < 3; ++i) {
// 				SphereEdgeLenMean += (PreNodes[t[i]] - PreNodes[t[(i + 1) % 3]]).getNorm();
// 				SphereNumEdges++;
// 			}
// 		}
// 		SphereEdgeLenMean /= (double)SphereNumEdges;

		double CoincidentCheckEpsilon = SphereEdgeLenMean * 0.2;
		double CoincidentCheckEpsilonSqr = CoincidentCheckEpsilon * CoincidentCheckEpsilon;

		vector<vector<vec3> > AllIntSurfPoints;
		vector<int> IntSurfZoneNums;
		vector<FESurface_c> IntSurfs;
		vector<vector<GradPath_c> > IntSurfRingCagePaths;


		int NumCompleted = 0;
		bool UserQuit = false;
		high_resolution_clock::time_point Time1, Time2;

		/*
		 *	Get sphere intersecting ring surfaces
		 */
		vector<int> RingSurfZoneNums;
		int NumZones = TecUtilDataSetGetNumZones();
#ifdef SUPERDEBUG
// 		vector<int> ZoneNumVec = { 100,118,260,293,319,429 };
// 		vector<int> ZoneNumVec = { 302,346,377,405 };
		vector<int> ZoneNumVec = { 238 };
		ZoneNumVec.resize(NumZones);
		for (int z = 1; z <= NumZones; ++z)
			ZoneNumVec[z - 1] = z;

		for (int z : ZoneNumVec){
 #else
 		for (int z = 1; z <= NumZones; ++z) {
 #endif
  			if (AuxDataZoneItemMatches(z, CSMAuxData.CC.ZoneSubType, CSMAuxData.CC.ZoneSubTypeRSSegment)) {
  				for (int ci = 0; ci < 3; ++ci) {
  					if (AuxDataZoneItemMatches(z, CSMAuxData.CC.ZFSCornerCPTypes[ci], CPNameList[CPTypeNum_Nuclear])
  						&& AuxDataZoneItemMatches(z, CSMAuxData.CC.ZFSCornerCPNumStrs[ci], to_string(CPNum)))
  					{
  						RingSurfZoneNums.push_back(z);
  						break;
  					}
  				}
  			}
//  			if (AuxDataZoneItemMatches(z, CSMAuxData.CC.ZoneSubType, CSMAuxData.CC.ZoneSubTypeRS)) {
// 				RingSurfZoneNums.push_back(z);
// 			}
		}

		std::map<int, int> RingCPNumToRingLineNum;

		SphereEdgeLenMean *= 0.85;

		NumCompleted = 0;
		TmpString = ProgressStr1.str() + "Finding ring surface intersections";
		StatusLaunch(TmpString, AddOnID, TRUE, TRUE, TRUE);
		Time1 = high_resolution_clock::now();
		for (int z : RingSurfZoneNums) {
			UserQuit = !StatusUpdate(NumCompleted, RingSurfZoneNums.size(), TmpString, AddOnID, Time1);
			if (UserQuit){
				if (TecUtilDialogMessageBox("Cancel run? All progress will be lost if you press \"Yes.\"", MessageBoxType_YesNo)) {
					break;
				}
				else {
					StatusLaunch(TmpString, AddOnID, TRUE, TRUE, TRUE);
					StatusUpdate(NumCompleted, RingSurfZoneNums.size(), TmpString, AddOnID, Time1);
					UserQuit = false;
				}
			}
			NumCompleted++;
			if (!UserQuit) {
				FESurface_c Surf(z, XYZVarNums);
				vector<vec3> IntPoints = Surf.GetSphereIntersectionPath(CPPos, Radius);
				if (IntPoints.size() > 1) {
					/*
						 * Now project each point back to the sphere radius
						 */
					for (auto & p : IntPoints) {
						p = CPPos + normalise(p - CPPos) * Radius;
					}
					vector<vec3> ResampledIntPoints;
					double PathLen = 0;
					for (int i = 0; i < IntPoints.size() - 1; ++i)
						PathLen += Distance(IntPoints[i], IntPoints[i + 1]);
					if (Vec3PathResample(IntPoints, MAX(int(PathLen / (SphereEdgeLenMean) ) + 1, 2), ResampledIntPoints)) {
#ifdef DEBUG_SAVEZONES
						SaveVec3VecAsScatterZone(IntPoints, "intersection with zone" + to_string(z));
						SetZoneStyle({ TecUtilDataSetGetNumZones() }, ZoneStyle_Path, (z % 60) + 1, 0.5);
						TecUtilZoneSetActive(Set(TecUtilDataSetGetNumZones()).getRef(), AssignOp_MinusEquals);
#endif
						IntPoints = ResampledIntPoints;

						/*
							 * Now project each point back to the sphere radius
							 */
						for (auto & p : IntPoints) {
							p = CPPos + normalise(p - CPPos) * Radius;
						}
						IntSurfs.push_back(Surf);
						IntSurfs.back().GeneratePointElementDistanceCheckData();

						/*
						 *	Get the ring-cage paths originating at the ring point.
						 *	Repeat the generation of the pair of ring path segments
						 *	with increasing distance of the seed point in each direction
						 *	until they terminate in different directions.
						 *
						 */
						int RingCPNum = -1;
						for (int i = 0; i < 3; ++i) {
							if (AuxDataZoneItemMatches(z, CSMAuxData.CC.ZFSCornerCPTypes[i], CPNameList[CPTypeNum_Ring])) {
								RingCPNum = stoi(AuxDataZoneGetItem(z, CSMAuxData.CC.ZFSCornerCPNumStrs[i])) - 1;
								break;
							}
						}
						if (RingCPNum >= 0) {
							if (AuxDataZoneItemMatches(z, CSMAuxData.CC.ZoneSubType, CSMAuxData.CC.ZoneSubTypeRS))
								RingCPNum = CPs.GetTypeNumOffsetFromTotOffset(RingCPNum)[1];
							vec3 RingCPPos, PrincDir;
							RingCPPos = CPs.GetXYZ(CPTypeNum_Ring, RingCPNum);
							vec3 EigVals;
							mat33 EigVecs;
							CalcEigenSystemForPoint(RingCPPos, EigVals, EigVecs, MR);
							PrincDir = EigVecs.col(0);

							int CPGlobalNum = CPs.GetTotOffsetFromTypeNumOffset(CPTypeNum_Ring, RingCPNum);
							double CPRho = ValAtPointByPtr(RingCPPos, VolInfo, RhoPtr);

							RingCPNumToRingLineNum[RingCPNum] = IntSurfRingCagePaths.size();
							IntSurfRingCagePaths.emplace_back();

							double iter = 1;
							double offset = 0.02;
							double dPdt = 1;
							while (dPdt > 0.5) {
								double dir = -1.0;
								IntSurfRingCagePaths.back() = vector<GradPath_c>(2);
								for (auto & GP : IntSurfRingCagePaths.back()) {
									vec3 SeedPt = RingCPPos + (PrincDir * offset * iter * dir);
									GP.SetupGradPath(
										SeedPt,
										StreamDir_Reverse, NumGPPoints * 2, GPType_Classic, GPTerminate_AtCP, nullptr,
										&CPs, &GPTermRadius, &CutoffVal, VolInfo, HessPtrs, GradPtrs, RhoPtr);
									GP.SetTerminalCPTypeNum(CPTypeNum_Cage);
									GP.SetStartEndCPNum(CPGlobalNum, 0);
									GP.Seed(false);
									GP.PointPrepend(RingCPPos, CPRho);


									dir += 2.0;
								}
								dPdt = dot(IntSurfRingCagePaths.back()[0].XYZAt(-1) - IntSurfRingCagePaths.back()[0].XYZAt(-2), IntSurfRingCagePaths.back()[1].XYZAt(-1) - IntSurfRingCagePaths.back()[1].XYZAt(-2));
								iter++;
							}
						}
						else {
							TecUtilDialogErrMsg("Failed to find ring CP for sphere-intersecting SZFS zone");
						}

#ifdef DEBUG_SAVEZONES
						SaveVec3VecAsScatterZone(IntPoints, "intersection with zone resampled" + to_string(z));
						SetZoneStyle({ TecUtilDataSetGetNumZones() }, ZoneStyle_Path, (z % 60) + 1, 0.5);
						TecUtilZoneSetActive(Set(TecUtilDataSetGetNumZones()).getRef(), AssignOp_MinusEquals);
#endif

						AllIntSurfPoints.push_back(IntPoints);
						IntSurfZoneNums.push_back(z);
					}
				}
			}
		}
		if (UserQuit) {
			TecUtilDataLoadEnd();
			StatusDrop(AddOnID);
			return;
		}

#ifdef DEBUG_SAVEZONES

		{
			vector<vec3> tmp;
			for (auto const & p : AllIntSurfPoints) {
				tmp.insert(tmp.end(), p.cbegin(), p.cend());
			}

			SaveVec3VecAsScatterZone(tmp, "constraint nodes before removing dups", White_C);
			SetZoneStyle({ TecUtilDataSetGetNumZones() }, ZoneStyle_Path, Blue_C, 0.8);
			TecUtilZoneSetScatter(SV_SHOW, Set(TecUtilDataSetGetNumZones()).getRef(), 0, TRUE);
			TecUtilZoneSetActive(Set(TecUtilDataSetGetNumZones()).getRef(), AssignOp_MinusEquals);
		}
#endif

		// Store map of edges to the constraint segment indices,
		// so that when a gradient bundle is split we can easily 
		// know if it lies on a constraint (ring surface intersection)
		// and which constraint.
		std::map<Edge, int> EdgeIntSurfNums;

		vector<vec3> SphereNodes;
		vector<vector<int> > SphereElems;
		vector<Index_t> OutConstrainedSegmentIndices;
		vector<vec3> constraintNodes;
		vector<Edge> constraintSegments;
		vector<tpcsm::Vec3> tmpConstraintNodes;
		vector<int> constraintSegmentSurfNums;
		
		// Store type of each node and element
		vector<GBATriangulatedSphereNodeElemType_e> NodeTypes, ElemTypes;

		// Update mesh for ring surface intersections if present
		if (!AllIntSurfPoints.empty()) 
		{

			// Move constraint nodes sufficiently close to bond path nodes to the bond path node
			for (auto & pList : AllIntSurfPoints) {
				for (auto & p : pList) {
					for (auto const & b : IntersectingBondPathPoints) {
						if (DistSqr(b, p) < CoincidentCheckEpsilonSqr) {
							p = b;
						}
					}
				}
			}

 			/*
 			 *	Check for nearly equal intersection path segment endpoints,
 			 *	i.e. where two paths intersect.
 			 *	These intersections should only occur at a bond-path-sphere intersection point,
 			 *	at which there is already a constrained sphere node, so when path intersections are
 			 *	found, move all the coincident points to the closest node on the sphere.
 			 *	
 			 *	Correction (2/6/2019): Some ring surfaces can converge well before reaching a common bond path
 			 *	intersection, so there can be points that are too close together basically anywhere
 			 *	between two ring surface intersection paths.
 			 *	So we'll check all points.
 			 *	This can be seen in simple aluminum.
 			 *	
 			 *	Correction 2: If two intersection paths converge before the common bond point then 
 			 *	their points may collide, breaking the triangulation update later.
 			 *	Need to look for colliding edges instead, then the would-be collisions
 			 *	
 			 *	Correction 3 (2/7/2019): This is a much larger problem than previously thought.
 			 *	In the event that two constraint segments run parallel and overlap,
 			 *	it means that the ring surface intersections share an edge on the sphere,
 			 *	and there is no one-to-one mapping of ring surface intersection constraint segments 
			 *	to sphere mesh edges.
			 *	There may even be cases where three or more ring surface intersections share duplicate edges
			 *	sufficiently close to a shared bond path intersection, though I haven't seen it yet.
 			 *	
 			 *	The problem is three-fold.
 			 *	1. simplify the duplicate constraint points and overlapping edges
 			 *	2. store ring surface identify and edge gps to for
 			 *		edges with more than one intersecting ring surface
 			 *	3. when gradient bundles are formed, detect which set of edge gps
 			 *		corresponds to which edge-coincident elements
 			 *	
 			 *	1) It's easy enough to identify and simplify the duplicate constraint nodes as well as to
 			 *	treat overlaping constraint edges, even when neither vertex of neither edge are coincident.
 			 *	
 			 *	2) Currently, I store which edges correspond to which ring surface intersections.
			 *	I could instead store a vector<int>, and in the event there are two or more constraint 
			 *	segments corresponding to the same sphere mesh edge, simply put them in the list.
			 *	Where are is only one constraining segment, nothing else needs to change.
			 *	
 			 *	Similarly, when edge gradient paths are produced, a vector<grad_path> is stored for each edge
 			 *	in a hash table according to edge string.
 			 *	It's also easy to instead store a 2d grad_path vector whose indices match the indices of stored
 			 *	ring surface indentifiers.
 			 *	
 			 *	3) Both involve identifying which set of edge GPs corresponds to one edge-coincident element
 			 *		or another by comparing edge gps in each set with the node gp opposite the edge in each 
 			 *		element. 
 			 *		I should be able to only compare to one edge gp in each set, but it may be safer
 			 *		to take the average path distance to all edge gps in each set.
 			 *		Then the question is when this should take place.
 			 *		I've thought to two approaches for this.
 			 *		i) only modify the specific code in the gradient bundle generation step, so all of this would
 			 *			happen inside the loop over non-bond path elements, then again in the bond path element
 			 *			loop
 			 *		ii) perform the identification for all elements before gradient bundle generation, storing the 
 			 *			results in a hash table to be used for non-bond path and bond path elements later on.
 			 *	
 			 *	Might be a second before I implement this. I'm currently very busy....
 			 */


			// Move ALL coincident constraint points to a common point.
// 			{
// 				double dupCheckDistSqr = CoincidentCheckEpsilonSqr * 0.001;
// 				for (int i = 0; i < AllIntSurfPoints.size() - 1; ++i) {
// 					for (int j = i + 1; j < AllIntSurfPoints.size(); ++j) {
// 						int ClosePtNum;
// 						for (int ji = 0; ji < AllIntSurfPoints[j].size(); ++ji) {
// 							auto & pj = AllIntSurfPoints[j][ji];
// 							if (DistSqr(ClosestPointToPath(AllIntSurfPoints[i], pj, ClosePtNum), pj) <= dupCheckDistSqr) {
// 								AllIntSurfPoints[j][ji] = AllIntSurfPoints[i][ClosePtNum];
// 							}
// 						}
// 					}
// 				}
// 			}

			// Move only coincident constraint points that are the terminal points of their
			// respective ring surface intersection paths.
			for (int i = 0; i < AllIntSurfPoints.size() - 1; ++i) {
				vector<std::pair<int, int> > CoincidentPointIndices;
				for (int j = i + 1; j < AllIntSurfPoints.size(); ++j) {
					if (DistSqr(AllIntSurfPoints[i][0], AllIntSurfPoints[j][0]) <= CoincidentCheckEpsilonSqr) {
						CoincidentPointIndices.emplace_back(i, 0);
						CoincidentPointIndices.emplace_back(j, 0);
					}
					else if (DistSqr(AllIntSurfPoints[i][0], AllIntSurfPoints[j].back()) <= CoincidentCheckEpsilonSqr) {
						CoincidentPointIndices.emplace_back(i, 0);
						CoincidentPointIndices.emplace_back(j, AllIntSurfPoints[j].size() - 1);
					}
					else if (DistSqr(AllIntSurfPoints[i].back(), AllIntSurfPoints[j][0]) <= CoincidentCheckEpsilonSqr) {
						CoincidentPointIndices.emplace_back(i, AllIntSurfPoints[i].size() - 1);
						CoincidentPointIndices.emplace_back(j, 0);
					}
					else if (DistSqr(AllIntSurfPoints[i].back(), AllIntSurfPoints[j].back()) <= CoincidentCheckEpsilonSqr) {
						CoincidentPointIndices.emplace_back(i, AllIntSurfPoints[i].size() - 1);
						CoincidentPointIndices.emplace_back(j, AllIntSurfPoints[j].size() - 1);
					}
				}

				if (!CoincidentPointIndices.empty()) {
					for (int ci = 0; ci < CoincidentPointIndices.size(); ci += 2) {
						auto p1 = CoincidentPointIndices[ci],
							p2 = CoincidentPointIndices[ci + 1];

						AllIntSurfPoints[p2.first][p2.second] = AllIntSurfPoints[p1.first][p1.second];
					}
				}
 			}

			for (int j = 0; j < AllIntSurfPoints.size(); ++j) {
				auto IntPoints = AllIntSurfPoints[j];
				tmpConstraintNodes.emplace_back(IntPoints[0][0], IntPoints[0][1], IntPoints[0][2]);
				for (int i = 1; i < IntPoints.size(); ++i) {
					tmpConstraintNodes.emplace_back(IntPoints[i][0], IntPoints[i][1], IntPoints[i][2]);
					constraintSegments.emplace_back(tmpConstraintNodes.size() - 2, tmpConstraintNodes.size() - 1);
					constraintSegmentSurfNums.push_back(j);
				}
			}


 			/*
 			 *	Need to remove duplicate constraint nodes and update edge connectivity to reflect the change.
 			 */
 			vector<int> nodeIndices(tmpConstraintNodes.size());
 			for (int i = 0; i < nodeIndices.size(); ++i)
 				nodeIndices[i] = i;
 			vector<tpcsm::Vec3> tmpNodes;
 			tmpNodes.reserve(tmpConstraintNodes.size());
			double dupCheckDistSqr = CoincidentCheckEpsilonSqr * 0.001;
 			for (int ni = 0; ni < tmpConstraintNodes.size() - 1; ++ni) {
 				bool isFound = false;
 				for (int nj = ni + 1; nj < tmpConstraintNodes.size(); ++nj) {
 					if ((tmpConstraintNodes[nodeIndices[ni]] - tmpConstraintNodes[nodeIndices[nj]]).getNormSquared() <= dupCheckDistSqr) {
 						/*
 						 *	duplicate node found, so search through edges, changing nj to ni when found.
 						 */
 						for (auto & e : constraintSegments) {
 							if (e.first == nj)
 								e.first = ni;
 
 							if (e.second == nj)
 								e.second = ni;
 						}
 						nodeIndices[nj] = nodeIndices[ni];
 						isFound = true;
 					}
 				}
 			}
 			/*
 			 * Now need to make the new list of nodes and update the edges again to reflect this.
 			 */
 			for (int ni = 0; ni < tmpConstraintNodes.size(); ++ni) {
 				if (nodeIndices[ni] == ni) {
 					/*
 					 *	Node is not duplicate
 					 */
 					tmpNodes.push_back(tmpConstraintNodes[ni]);
 					for (auto & e : constraintSegments) {
 						if (e.first == ni)
 							e.first = tmpNodes.size() - 1;
 
 						if (e.second == ni)
 							e.second = tmpNodes.size() - 1;
 					}
 				}
 			}
 			tmpConstraintNodes = tmpNodes;
 
 			/*
 			 *	Now move nodes on the sphere that are sufficiently close to a constraint node
 			 *	to eliminate the ugly triangles that result when two points are very close.
 			 */
 

			vector<std::pair<int, double> > MinDistNodes(PreNodes.size(), std::make_pair(-1, DBL_MAX));
			vector<bool> PreNodeMoved(PreNodes.size(), false), 
				ConstraintNodeMoved(tmpConstraintNodes.size(), false);
			for (int pi = 0; pi < tmpConstraintNodes.size(); ++pi) {
				auto p = tmpConstraintNodes[pi];
				double MinDistSqr = DBL_MAX, TmpDistSqr;
				int MinI;
				for (int ni = 0; ni < PreNodes.size(); ++ni) {
					auto n = PreNodes[ni];
					TmpDistSqr = (p - n).getNormSquared();
					if (TmpDistSqr < MinDistSqr) {
						MinDistSqr = TmpDistSqr;
						MinI = ni;
					}
				}
				if (MinDistNodes[MinI].first < 0 || (MinDistNodes[MinI].first >= 0 && MinDistSqr < MinDistNodes[MinI].second)){
					if (MinDistSqr < MinDistNodes[MinI].second){
						MinDistNodes[MinI] = std::make_pair(pi, MinDistSqr);
					}
				}
			}
			for (int ni = 0; ni < PreNodes.size(); ++ni){
				if (!PreNodeMoved[ni] && MinDistNodes[ni].second <= CoincidentCheckEpsilonSqr) {
					PreNodes[ni] = tmpConstraintNodes[MinDistNodes[ni].first];
					PreNodeMoved[ni] = true;
					ConstraintNodeMoved[MinDistNodes[ni].first] = true;
				}
			}

#ifdef DEBUG_SAVEZONES
			SphereNodes.reserve(PreNodes.size());
			for (auto const & i : PreNodes)
				SphereNodes.push_back({ i.x(), i.y(), i.z() });
			NumNodes = SphereNodes.size();
			SphereElems.reserve(PreElems.size());
			for (auto const & i : PreElems)
				SphereElems.push_back({ i[0], i[1], i[2] });
			FESurface_c TmpSphere(SphereNodes, SphereElems);
			TmpSphere.SaveAsTriFEZone({ 1,2,3 }, "Before triangulation update after node moving to constraint nodes");
			TecUtilZoneSetActive(Set(TecUtilDataSetGetNumZones()).getRef(), AssignOp_MinusEquals);
			SphereNodes.clear();
			SphereElems.clear();

			{
				vector<vec3> tmp;
				for (auto const & i : tmpConstraintNodes)
					tmp.push_back({ i.x(),i.y(),i.z() });

				SaveVec3VecAsScatterZone(tmp, "constraint nodes", White_C);
				SetZoneStyle({ TecUtilDataSetGetNumZones() }, ZoneStyle_Path, Blue_C, 0.8);
				TecUtilZoneSetScatter(SV_SHOW, Set(TecUtilDataSetGetNumZones()).getRef(), 0, TRUE);
				TecUtilZoneSetActive(Set(TecUtilDataSetGetNumZones()).getRef(), AssignOp_MinusEquals);
			}
#endif
// 			TecUtilDataLoadEnd();
// 			StatusDrop(AddOnID);
// 			return;

			/*
			 *	Now do the same check but for edges, that is, check for edges that are sufficiently
			 *	close to a constraint point. If found, split the edge and its two elements.	
			 */


			 // Keep track of sorted distances from edges to constraint points in two directions.
			// That is, for each edge, keep a sorted list (std::map) that maps distances to constraint nodes,
			// and for each constraint node, keep a sorted list that maps distances to edges.
			// Using both of these, we can make sure that the closest node/edge combinations are found and
			// used.
			// SphereEdgeTriNumsMap[edge].first = maps distance to constraint point to the point's index 
			// SphereEdgeTriNumsMap[edge].second = set of element indices (always two elements)
			CoincidentCheckEpsilon *= 0.5;
			std::map<Edge, std::pair<std::map<double,int>,std::set<int> > > SphereEdgeConstraintNodeDistanceTriNumsMap;
			std::vector<std::map<double, Edge> > MinDistEdges(tmpConstraintNodes.size());
			for (int ti = 0; ti < PreElems.size(); ++ti) {
				auto * t = &PreElems[ti];
				SphereEdgeConstraintNodeDistanceTriNumsMap[MakeEdge(t->v1(), t->v2())].second.insert(ti);
				SphereEdgeConstraintNodeDistanceTriNumsMap[MakeEdge(t->v1(), t->v3())].second.insert(ti);
				SphereEdgeConstraintNodeDistanceTriNumsMap[MakeEdge(t->v2(), t->v3())].second.insert(ti);
			}

			// Keep track of min distance from each constraint point to any edge
			for (auto & e : SphereEdgeConstraintNodeDistanceTriNumsMap){
				auto x1 = &PreNodes[e.first.first], x2 = &PreNodes[e.first.second];
				for (int pi = 0; pi < tmpConstraintNodes.size(); ++pi){
					auto x0 = &tmpConstraintNodes[pi];
					if (!ConstraintNodeMoved[pi]){
						double TmpDist = PointLineDist(*x0, *x1, *x2);
						if (TmpDist <= CoincidentCheckEpsilon) {
							double EdgeLen = (*x2 - *x1).getNormSquared();
							if ((*x0 - *x1).getNormSquared() < EdgeLen && (*x0 - *x2).getNormSquared() < EdgeLen) {
								MinDistEdges[pi][TmpDist] = e.first;
								e.second.first[TmpDist] = pi;
							}
						}
					}
				}
			}
			
			// Now loop over constraint nodes and if the closest edge is within tolerance
			// (they're all within tolerance),
			// split the edge, but only if the node is also the closest constraint node to the
			// edge. 
			// In the event that an edge is the closest to a node but there are nodes closer to the
			// edge, then remove that edge from the distance-edge map for the current node.
			// With this criteria, the list of closest nodes held by each edge take precedent over
			// the list of closest edges held by each node, and it is possible that all of the close-
			// enough edges to a node are "taken" by other nodes and that no edge will be split for
			// a particular node.
			// Put another way, a constraint node will only be used to split an edge if it is the closest
			// edge to it and it is the closest constraint node to the edge.
			for (int pi = 0; pi < tmpConstraintNodes.size(); ++pi){
				auto p = &MinDistEdges[pi];
				while (!p->empty()) {
					auto pe = p->cbegin()->second;
					auto e = SphereEdgeConstraintNodeDistanceTriNumsMap[pe];
					double EdgeLen = (PreNodes[pe.first] - PreNodes[pe.second]).getNormSquared();
					if (PointLineDist(tmpConstraintNodes[pi], PreNodes[pe.first], PreNodes[pe.second]) <= CoincidentCheckEpsilon 
						&& (tmpConstraintNodes[pi] - PreNodes[pe.first]).getNormSquared() < EdgeLen && (tmpConstraintNodes[pi] - PreNodes[pe.second]).getNormSquared() < EdgeLen
						&& e.first.cbegin()->second == pi)
					{
						// Split the edge
						int ni = PreNodes.size();
						PreNodes.push_back(tmpConstraintNodes[pi]);
						PreNodeMoved.push_back(true);

						// Update elements.
						for (auto const & ti : e.second){
							// Find opposite corner
							for (int ci = 0; ci < 3; ++ci){
								int c1 = PreElems[ti][ci],
									c2 = PreElems[ti][(ci + 1) % 3];
								if (p->cbegin()->second == MakeEdge(c1, c2)){
									// Split element
									int c3 = PreElems[ti][(ci + 2) % 3];
									PreElems[ti] = TriNodes(c1, c3, ni);
									PreElems.emplace_back(c2, ni, c3);
								}
							}
						}

						// Now loop over all other edges in SphereEdgeTriNumsMap and update any edges of
						// the new or modified elements to reflect their new triangles
						for (auto & e2 : SphereEdgeConstraintNodeDistanceTriNumsMap){
							for (int const & ti : e.second){
								if (e2.second.second.count(ti)){
									e2.second.second.clear();
									for (int tj = 0; tj < PreElems.size() && e2.second.second.size() < 2; ++tj){
										for (int ci = 0; ci < 3; ++ci){
											int n1 = PreElems[tj][ci],
												n2 = PreElems[tj][(ci + 1) % 3];
											if (MakeEdge(n1,n2) == e2.first){
												e2.second.second.insert(tj);
												break;
											}
										}
									}
								}
							}
						}

						break;
					}
					p->erase(p->begin());
				}
			}

#ifdef DEBUG_SAVEZONES
			SphereNodes.reserve(PreNodes.size());
			for (auto const & i : PreNodes)
				SphereNodes.push_back({ i.x(), i.y(), i.z() });
			NumNodes = SphereNodes.size();
			SphereElems.reserve(PreElems.size());
			for (auto const & i : PreElems)
				SphereElems.push_back({ i[0], i[1], i[2] });
			TmpSphere = FESurface_c(SphereNodes, SphereElems);
			TmpSphere.SaveAsTriFEZone({ 1,2,3 }, "Before triangulation update after edge splitting");
			TecUtilZoneSetActive(Set(TecUtilDataSetGetNumZones()).getRef(), AssignOp_MinusEquals);
			SphereNodes.clear();
			SphereElems.clear();
#endif

			/*
			 *	Do the same check again, but this time checking sphere nodes to see if they are close enough to
			 *	a constraint segment.
			 *	If so, project the node to the segment and add it to tmpConstraintNodes and split the segment.
			 */
// 			CoincidentCheckEpsilon *= 2.0;
			PreNodeMoved.resize(PreNodes.size(), false); // non-destructive resizing

			int oldNumConstraintSegments = 0;
			while (oldNumConstraintSegments < constraintSegments.size()) {
				vector<std::map<double, int> > MinDistConstraintSegmentToPreNodes(constraintSegments.size());
				vector<std::map<double, int> > MinDistPreNodeToConstraintSegment(PreNodes.size());
				for (int ei = 0; ei < constraintSegments.size(); ++ei) {
					auto e = constraintSegments[ei];
					for (int ni = 0; ni < PreNodes.size(); ++ni) {
						if (!PreNodeMoved[ni]) {
							double TmpDist = PointLineDist(PreNodes[ni], tmpConstraintNodes[e.first], tmpConstraintNodes[e.second]);
							if (TmpDist <= CoincidentCheckEpsilon) {
								double EdgeLen = (tmpConstraintNodes[e.first] - tmpConstraintNodes[e.second]).getNorm();
								if ((PreNodes[ni] - tmpConstraintNodes[e.second]).getNorm() < EdgeLen
									&& (tmpConstraintNodes[e.first] - PreNodes[ni]).getNorm() < EdgeLen)
								{
									MinDistConstraintSegmentToPreNodes[ei][TmpDist] = ni;
									MinDistPreNodeToConstraintSegment[ni][TmpDist] = ei;
								}
							}
						}
					}
				}

				oldNumConstraintSegments = constraintSegments.size();
				for (int ei = 0; ei < oldNumConstraintSegments; ++ei) {
					double EdgeLen = (tmpConstraintNodes[constraintSegments[ei].first] - tmpConstraintNodes[constraintSegments[ei].second]).getNorm();
					while (!MinDistConstraintSegmentToPreNodes[ei].empty()) {
						auto MinDist = MinDistConstraintSegmentToPreNodes[ei].cbegin();
						if (!PreNodeMoved[MinDist->second] && MinDistPreNodeToConstraintSegment[MinDist->second].cbegin()->second == ei
							&& PointLineDist(PreNodes[MinDist->second], tmpConstraintNodes[constraintSegments[ei].first], tmpConstraintNodes[constraintSegments[ei].second]))
						{
							if ((PreNodes[MinDist->second] - tmpConstraintNodes[constraintSegments[ei].second]).getNorm() < EdgeLen
								&& (tmpConstraintNodes[constraintSegments[ei].first] - PreNodes[MinDist->second]).getNorm() < EdgeLen)
							{
								// Project node to constraint segment
								PreNodeMoved[MinDist->second] = true;
								PreNodes[MinDist->second] = ProjectPointToLine(PreNodes[MinDist->second], tmpConstraintNodes[constraintSegments[ei].first], tmpConstraintNodes[constraintSegments[ei].second]);
								PreNodes[MinDist->second] = SphereCenterVec + (PreNodes[MinDist->second] - SphereCenterVec).normalize() * Radius;

								// Add node to tmpConstaintNodes and split the segment
								int newNodeNum = tmpConstraintNodes.size();
								constraintSegmentSurfNums.push_back(constraintSegmentSurfNums[ei]);
								tmpConstraintNodes.push_back(PreNodes[MinDist->second]);
								constraintSegments.push_back(MakeEdge(constraintSegments[ei].first, newNodeNum));
								constraintSegments[ei] = MakeEdge(constraintSegments[ei].second, newNodeNum);

								break;
							}
						}
						MinDistConstraintSegmentToPreNodes[ei].erase(MinDist);
					}
				}
			}

#ifdef DEBUG_SAVEZONES
			SphereNodes.reserve(PreNodes.size());
			for (auto const & i : PreNodes)
				SphereNodes.push_back({ i.x(), i.y(), i.z() });
			NumNodes = SphereNodes.size();
			SphereElems.reserve(PreElems.size());
			for (auto const & i : PreElems)
				SphereElems.push_back({ i[0], i[1], i[2] });
			TmpSphere = FESurface_c(SphereNodes, SphereElems);
			TmpSphere.SaveAsTriFEZone({ 1,2,3 }, "Before triangulation update after moving nodes to constraint segments");
			TecUtilZoneSetActive(Set(TecUtilDataSetGetNumZones()).getRef(), AssignOp_MinusEquals);
			SphereNodes.clear();
			SphereElems.clear();
#endif

// 			TecUtilDataLoadEnd();
// 			StatusDrop(AddOnID);
// 			return;


			/*
			 *	Update triangulation to include ring surface intersections as edges
			 */

			vector<tpcsm::Vec3> OutXYZs;
			vector<TriNodes> OutTris;
			char const * StatusMessage;


			if (tpcsm::updateSphericalTriangulation(PreNodes,
				PreElems,
				SphereCenterVec,
				Radius,
				tmpConstraintNodes,
				constraintSegments,
				false,
				OutXYZs,
				OutConstrainedSegmentIndices,
				OutTris,
				StatusMessage))
			{
				SphereNodes.reserve(OutXYZs.size());
				for (auto const & i : OutXYZs)
					SphereNodes.push_back({i.x(), i.y(), i.z()});

				NumNodes = SphereNodes.size();

				SphereElems.reserve(OutTris.size());
				for (auto const & i : OutTris)
					SphereElems.push_back({ i[0], i[1], i[2] });

				NumElems = SphereElems.size();

				constraintNodes.reserve(tmpConstraintNodes.size());
				for (auto const & i : tmpConstraintNodes)
					constraintNodes.push_back({ i.x(), i.y(), i.z() });
			}
			else{
				TecUtilDialogErrMsg(string("Failed to update mesh with ring surface intersections: " + string(StatusMessage)).c_str());
				continue;
			}
		}
		else {
			SphereNodes.reserve(PreNodes.size());
			for (auto const & i : PreNodes)
				SphereNodes.push_back({i.x(), i.y(), i.z()});

			SphereElems.reserve(PreElems.size());
			for (auto const & i : PreElems)
				SphereElems.push_back({i[0], i[1], i[2]});

			OutConstrainedSegmentIndices.resize(SphereNodes.size(), -1);
		}

		NumNodes = SphereNodes.size();
		NumElems = SphereElems.size();

#ifdef DEBUG_SAVEZONES
  		FESurface_c TmpSphere(SphereNodes, SphereElems);
		TmpSphere.SaveAsTriFEZone({ 1,2,3 }, "after triangulation update 1");
		TecUtilZoneSetActive(Set(TecUtilDataSetGetNumZones()).getRef(), AssignOp_MinusEquals);
#endif

		// Check for and remove any duplicate nodes
// 		vector<int> SphereNodesOldToNew;
// 		RemoveDupicateNodesFromMesh(SphereNodes, SphereElems, &SphereNodesOldToNew);

		NumNodes = SphereNodes.size();
		NumElems = SphereElems.size();

// 		TmpSphere = FESurface_c(SphereNodes, SphereElems);
// 		TmpSphere.SaveAsTriFEZone(XYZVarNums, "after triangulation update 1 minus dup nodes");

		// The triangulation update can sometimes just ignore constraint nodes.
		// The forgotten nodes still appear in SphereNodes, but they might not be
		// in any elements, so they're spurious.
		// We only care about such nodes if they happen to also be constraint nodes,
		// in which case they'll be discovered in the following loop.
		// To identify them as having been ignored, we'll prepare a list
		// of all the edges in the mesh, then the constraint edges found
		// below can be checked for occurrence in the mesh.
		// An ignored node will necessary fall along an existing edge,
		// so that edge can be found and split at the ignored node and
		// the neighboring elements updated to reflect the existence of
		// the ignored node.
		std::map<Edge, std::set<int> > SphereEdgeTriNumsMap;
		std::set<int> SphereNodeIndSet;
		if (!AllIntSurfPoints.empty())
		{
			for (int ti = 0; ti < SphereElems.size(); ++ti) {
				auto * t = &SphereElems[ti];
				for (int ci = 0; ci < 3; ++ci) {
					Edge e = MakeEdge((*t)[ci], (*t)[(ci + 1) % 3]);
					SphereEdgeTriNumsMap[e].insert(ti);
					SphereNodeIndSet.insert((*t)[ci]);
				}
			}
		}

		// Make sphere
		FESurface_c Sphere(SphereNodes, SphereElems);
		std::set<int> ConstrainedEdgeNodes;
		if (!AllIntSurfPoints.empty())
		{
			// Store constraingSegments into the EdgeConstraint map
			for (int i = 0; i < constraintSegments.size(); ++i) {
				vec3 p1 = { tmpConstraintNodes[constraintSegments[i].first].x(), tmpConstraintNodes[constraintSegments[i].first].y(), tmpConstraintNodes[constraintSegments[i].first].z() },
					p2 = { tmpConstraintNodes[constraintSegments[i].second].x(), tmpConstraintNodes[constraintSegments[i].second].y(), tmpConstraintNodes[constraintSegments[i].second].z() };
				int n1 = Sphere.GetClosestNodeToPoint(p1),
					n2 = Sphere.GetClosestNodeToPoint(p2);
				Edge e = MakeEdge(n1, n2);
				for (auto ni : { n1,n2 }) {
					OutConstrainedSegmentIndices[ni] = i;
					ConstrainedEdgeNodes.insert(ni);
				}
				if (SphereEdgeTriNumsMap.count(e))
					EdgeIntSurfNums[e] = constraintSegmentSurfNums[i];
				else {
					// This constrained edge is not a current edge in the mesh.
					// Loop over the edges of the mesh to find the one for which
					// this edge is half.
					// Do this by identifying the node that is missing from the sphere,
					// then finding the edge on the sphere that the missing node lies
					// within.
					int newNode = -1;
					if (!SphereNodeIndSet.count(n1))
						newNode = n1;
					else if (!SphereNodeIndSet.count(n2))
						newNode = n2;
					if (newNode >= 0) {
						SphereNodeIndSet.insert(newNode);
						double MinDotPdt = DBL_MAX;
						Edge MinEdge;
						for (const auto & ei : SphereEdgeTriNumsMap) {
							double TmpDotPdt = dot(normalise(SphereNodes[ei.first.first] - SphereNodes[newNode]), normalise(SphereNodes[ei.first.second] - SphereNodes[newNode]));
							if (TmpDotPdt < MinDotPdt) {
								MinDotPdt = TmpDotPdt;
								MinEdge = ei.first;
							}
						}

						// Add ring surface info to new edges
						EdgeIntSurfNums[MakeEdge(MinEdge.first, newNode)] = constraintSegmentSurfNums[i];
						EdgeIntSurfNums[MakeEdge(MinEdge.second, newNode)] = constraintSegmentSurfNums[i];

						// Correct edge found. Now split the two triangles sharing that edge.
						for (auto ti : SphereEdgeTriNumsMap[MinEdge]) {
							// Find corner opposite the edge
							for (int ci = 0; ci < 3; ++ci) {
								int c1 = (ci + 1) % 3;
								if (MakeEdge(SphereElems[ti][ci], SphereElems[ti][c1]) == MinEdge) {
									int c2 = (ci + 2) % 3;
									SphereElems.push_back({ SphereElems[ti][c2], SphereElems[ti][ci], newNode });
									SphereElems[ti] = { SphereElems[ti][c2], newNode, SphereElems[ti][c1] };
								}
							}
						}
					}
					else{
						// In this case, the updated triangulation split an existing edge that intersected a constraint segment without keeping track of the split.
						// The node is already in the sphere and is along the current constraint segment.
						// A later check will correct this.
						// Here, just add the edge as is.
						EdgeIntSurfNums[e] = constraintSegmentSurfNums[i];
					}
				}
			}

#ifdef DEBUG_SAVEZONES
			TmpSphere = FESurface_c(SphereNodes, SphereElems);
			TmpSphere.SaveAsTriFEZone(XYZVarNums, "after triangulation update lost constraint nodes found");
			TecUtilZoneSetActive(Set(TecUtilDataSetGetNumZones()).getRef(), AssignOp_MinusEquals);
#endif
		}

		if (!AllIntSurfPoints.empty()) {
			// Another pass to remove tris that have an edge that cuts interior to the sphere.
			// This happens when the triangulation update code makes an element from three nodes
			// along a single ring surface intersection path.
			// To check, loop over every element and compare against every ring surface intersection
			// path.
			// 
			double CheckTotalMinDistance = 1e-10;
			for (int ti = 0; ti < NumElems; ++ti) {
				bool BadElem = false;
				for (auto const & p : AllIntSurfPoints) {
					// For each element, loop over each intersection path
					for (int pi = 2; pi < p.size() && !BadElem; ++pi) {
						// Iterate down triplets of points of the intersection path by triplets starting at 2, so the 
						// first triplet is (0,1,2).
						// Maintain the closest distance of each element node to any of the three
						// points being checked on the intersection path.
						vec3 MinDistances = { DBL_MAX, DBL_MAX, DBL_MAX };
						for (int ci = 0; ci < 3; ++ci) {
							// For each corner of the element, find the minimum distance to any of the points
							// being checked on the intersection path.
							for (int cj = pi - 2; cj <= pi; ++cj) { // pi starts as 2, so the lowest cj will be 0
								MinDistances[ci] = MIN(MinDistances[ci], Distance(SphereNodes[SphereElems[ti][ci]], p[cj]));
							}
						}
						// If the sum of minimum distances is very small then the element is assumed to be bad and is removed from the sphere.
						if (sum(MinDistances) <= CheckTotalMinDistance) {
							SphereElems[ti].clear();
							BadElem = true;
						}
					}
					if (BadElem)
						break;
				}
			}

			vector<vector<int> > TmpElems;
			for (auto const & e : SphereElems) {
				if (!e.empty()) {
					TmpElems.push_back(e);
				}
			}

			SphereElems = TmpElems;
			NumElems = SphereElems.size();

#ifdef DEBUG_SAVEZONES
			TmpSphere = FESurface_c(SphereNodes, SphereElems);
			TmpSphere.SaveAsTriFEZone({ 1,2,3 }, "after triangulation update remove constraint segment tris");
			TecUtilZoneSetActive(Set(TecUtilDataSetGetNumZones()).getRef(), AssignOp_MinusEquals);
#endif
		}

		if (!AllIntSurfPoints.empty())
		{
			// The updated mesh may contain new nodes along constraint edges,
			// that is, a constraint edge provided to the triangulation update
			// had a new node inserted along it.
			// These new nodes will necessarily appear along a straight line connecting
			// the nodes of the edge on which they were created, so we can find the 
			// constraint edge that sandwiches the new node by checking the dot product
			// between the new node and edge nodes.
			// The vectors pointing from the edge nodes to the new node should be
			// pointing in the exact opposite direction, so the edge that gives a normalised
			// dot product of -1 is the correct edge.
			// Need to find those and add their respective edges to EdgeIntSurfNums.
			// Find the right constrained edge by finding the minimum dot product of
			// the normalized vectors from edge nodes to each constraint node.
			for (int ni = 0; ni < NumNodes; ++ni) {
				if (OutConstrainedSegmentIndices[ni] >= 0 && !ConstrainedEdgeNodes.count(ni)) {
					double MinDotPdt = DBL_MAX, TmpDotPdt;
					Edge MinDotPdtEdge;
					for (const auto & EdgeIntSurfNum : EdgeIntSurfNums) {
						vec3 v1 = normalise(SphereNodes[EdgeIntSurfNum.first.first] - SphereNodes[ni]);
						vec3 v2 = normalise(SphereNodes[EdgeIntSurfNum.first.second] - SphereNodes[ni]);
						double TmpDotPdt = dot(v1, v2);
						if (TmpDotPdt < MinDotPdt) {
							MinDotPdt = TmpDotPdt;
							MinDotPdtEdge = EdgeIntSurfNum.first;
						}
					}
					int SurfNum = EdgeIntSurfNums[MinDotPdtEdge];
					EdgeIntSurfNums[MakeEdge(MinDotPdtEdge.first, ni)] = SurfNum;
					EdgeIntSurfNums[MakeEdge(MinDotPdtEdge.second, ni)] = SurfNum;
					ConstrainedEdgeNodes.insert(ni);
					EdgeIntSurfNums.erase(MinDotPdtEdge);
				}
			}
		}


		Sphere = FESurface_c(SphereNodes, SphereElems);
		/*
		 *	Find nodes on new sphere that correspond to the bond path intersections
		 *	with the original sphere.
		 */
		vector<int> IntBondPathNodes(IntersectingBondPathPoints.size());
		std::unordered_map<int, GradPath_c> IntBondPathSegments; // stored according to intersecting node index
		for (int bi = 0; bi < IntersectingBondPathPoints.size(); ++bi) {
			double ClosestNodeDist;
			IntBondPathNodes[bi] = Sphere.GetClosestNodeToPoint(IntersectingBondPathPoints[bi], &ClosestNodeDist);
			IntBondPathSegments[IntBondPathNodes[bi]] = IntersectingBondPaths[bi];
		}

// 		if (!AllIntSurfPoints.empty())
// 		{
// 			// The updated may also contain new nodes along constraint edges that,
// 			// unlike those being check for above, are present in one or more elements
// 			// but are also absent from other elements that should have been bisected
// 			// by the new node.
// 			// As above, these will necessarily lie along existing edges.
// 			// To find them, loop over all R or RB type elements, and for each constrained
// 			// edge of an element, check all constrained nodes to see if any lie interior
// 			// to the edge.
// 			for (int ti = 0; ti < SphereElems.size(); ++ti) {
// 				bool SplitElement = false;
// 				for (int ci = 0; ci < 3 && !SplitElement; ++ci) {
// 					int n1 = SphereElems[ti][ci];
// 					int n2 = SphereElems[ti][(ci + 1) % 3];
// 					if ((OutConstrinedSegmentIndices[n1] >= 0 && OutConstrinedSegmentIndices[n1] == OutConstrinedSegmentIndices[n2])
// 						|| (OutConstrinedSegmentIndices[n1] >= 0 || IntBondPathSegments.count(n2))
// 						|| (OutConstrinedSegmentIndices[n2] >= 0 || IntBondPathSegments.count(n1)))
// 					{
// 						// Constrained edge found, now loop over all constrained nodes and check
// 						// if each node is interior to the edge by checking the dot product of 
// 						// the vectors pointing from each edge node to the check node.
// 						for (int ni = 0; ni < NumNodes; ++ni) {
// 							if (ni != n1 && ni != n2 && OutConstrinedSegmentIndices[ni] >= 0) {
// 								vec3 v1 = normalise(SphereNodes[ni] - SphereNodes[n1]);
// 								vec3 v2 = normalise(SphereNodes[ni] - SphereNodes[n2]);
// 								double DotPdt = dot(v1, v2);
// 								if (DotPdt < -0.95) {
// 									// Partially missing node found, so split the edge/element
// 									// Get opposite corner of element
// 									int n3 = SphereElems[ti][(ci + 2) % 3];
// 									if (ni != n3) {
// 										// Add first half as new element
// 										SphereElems.push_back({ n1,ni,n3 });
// 										// Replace old element with new second half
// 										SphereElems[ti] = { ni,n2,n3 };
// 										// Add edges to EdgeIntSurfNums
// 										Edge OldEdge = MakeEdge(n1, n2),
// 											NewEdge1 = MakeEdge(n1, ni),
// 											NewEdge2 = MakeEdge(ni, n2);
// 										int IntSurfNum;
// 										if (EdgeIntSurfNums.count(OldEdge)) {
// 											IntSurfNum = EdgeIntSurfNums[OldEdge];
// 										}
// 										else if (EdgeIntSurfNums.count(NewEdge1)) {
// 											IntSurfNum = EdgeIntSurfNums[NewEdge1];
// 										}
// 										else if (EdgeIntSurfNums.count(NewEdge2)) {
// 											IntSurfNum = EdgeIntSurfNums[NewEdge2];
// 										}
// 										else {
// 											IntSurfNum = constraintSegmentSurfNums[OutConstrinedSegmentIndices[ni]];
// 										}
// 										EdgeIntSurfNums[MakeEdge(n1, ni)] = IntSurfNum;
// 										EdgeIntSurfNums[MakeEdge(ni, n2)] = IntSurfNum;
// 										EdgeIntSurfNums.erase(OldEdge);
// 										SplitElement = true;
// 										ti--;
// 										break;
// 									}
// 								}
// 							}
// 						}
// 					}
// 				}
// 			}
// #ifdef DEBUG_SAVEZONES
// 			TmpSphere = FESurface_c(SphereNodes, SphereElems);
// 			TmpSphere.SaveAsTriFEZone(XYZVarNums, "after triangulation update nodes along constraint edges found");
// #endif
// 		}

		NumElems = SphereElems.size();

		

		if (!AllIntSurfPoints.empty())
		{
			// Pass over elems/nodes to make sure all constrained edges are in EdgeIntSurfNums.
			// This is to check for edges of R-B type nodes, but we'll include all constrained
			// edges  in the check.
			for (int ti = 0; ti < NumElems; ++ti) {
				for (int ci = 0; ci < 3; ++ci) {
					int n1 = SphereElems[ti][ci];
					int n2 = SphereElems[ti][(ci + 1) % 3];
					Edge TmpEdge = MakeEdge(n1, n2);
					if (!EdgeIntSurfNums.count(TmpEdge))
					{
						if (OutConstrainedSegmentIndices[n1] >= 0 && OutConstrainedSegmentIndices[n1] == OutConstrainedSegmentIndices[n2]) {
							EdgeIntSurfNums[TmpEdge] = constraintSegmentSurfNums[OutConstrainedSegmentIndices[n1]];
						}
						else if (OutConstrainedSegmentIndices[n2] >= 0 && IntBondPathSegments.count(n1)) {
							EdgeIntSurfNums[TmpEdge] = constraintSegmentSurfNums[OutConstrainedSegmentIndices[n2]];
						}
						else if (OutConstrainedSegmentIndices[n1] >= 0 && IntBondPathSegments.count(n2)) {
							EdgeIntSurfNums[TmpEdge] = constraintSegmentSurfNums[OutConstrainedSegmentIndices[n1]];
						}
					}
				}
			}
		}

		/*
		*	Now we need to seed the gradient paths from the nodes and along the
		*	edges of the triangular mesh.
		*	***actually, don't put gradient paths along the edges. It introduces
		*		unnecessary complexity, makes the algorithm much more expensive,
		*		and only captures a little bit more information about the curvature
		*		of the zero flux surfaces formed from stitching the gradient
		*		paths together.
		*		This also simplifies the tetrahedral decomposition used in the
		*		integration code.
		*		The information lost by excluding edges can be recovered while only using
		*		the nodes of triangles by simply refining the triangulation.
		*
		*	There are three types of nodes, according to the treatment of the paths
		*	they form:
		*		C type:		Paths that terminate at a cage CP or at a rho cutoff
		*					require no special treatment.
		*		R type:		Paths that coincide with a ring surface become will be seeded
		*					as surface paths and terminate at a ring point, then
		*					then split into the two halves of the ring path away from the
		*					ring point to the cage.
		*		B type:		Paths that coincide with a bond path will be replaced with the
		*					bond path, then a new gradient path will be seeded away from
		*					the bond point in the bond plane (interatomic surface) for each
		*					of the nodes that share an edge with the node of the original
		*					path.
		*
		*	We'll store all the gradient paths in a single vector and keep a new list of triangles
		*	where each triangle is defined by the indices of the gradient paths that correspond to
		*	its nodes.
		*	R, and B type nodes will result in new gradient paths being added to the list.
		*/


		// reinterpolate rho values and force to be monotomic for bond path segments
		for (auto & b : IntBondPathSegments){
			b.second.MakeRhoValuesMonotonic(&VolInfo, &RhoPtr);
			if (b.second.RhoAt(0) < b.second.RhoAt(-1)){
				b.second.Reverse();
			}
		}

		// Get all node types;
		NodeTypes.resize(NumNodes, NETypeInvalid);
// #pragma omp parallel for

		for (int ni = 0; ni < NumNodes; ++ni) {
			if (IntBondPathSegments.count(ni))
				NodeTypes[ni] = NETypeB;
			else if (OutConstrainedSegmentIndices[ni] >= 0)
				NodeTypes[ni] = NETypeR;
			else {
				NodeTypes[ni] = NETypeC;
				for (int i = 0; i < constraintNodes.size() && NodeTypes[ni] != NETypeR; ++i) {

					if (DistSqr(SphereNodes[ni], constraintNodes[i]) < CoincidentCheckEpsilonSqr) {
						for (int j = 0; j < constraintSegments.size(); ++j) {
							if (constraintSegments[j].first == i || constraintSegments[j].second == i) {
								OutConstrainedSegmentIndices[ni] = j;
								NodeTypes[ni] = NETypeR;
								break;
							}
						}
					}
				}
			}
		}

		// For B type nodes:
		// Get node connectivity and prepare a list of
		// GP lists so that all GPs coming from a B type
		// node have the same number of points in the
		// bond path segment.
		// 
// 		vector<std::set<int> > NodeConnectivity;
// 		GetMeshNodeConnectivity(*Sphere.GetElemListPtr(), NumNodes, NodeConnectivity);

		
// 		TecUtilDataLoadEnd();
// 		StatusDrop(AddOnID);
// 		return;
 		

		ElemTypes.resize(NumElems, NETypeInvalid);
// #pragma omp parallel for
		for (int ti = 0; ti < NumElems; ++ti)
			for (int ni : SphereElems[ti])
				ElemTypes[ti] = MAX(ElemTypes[ti], NodeTypes[ni]);

//  		if (!AllIntSurfPoints.empty()){
//  			// Another pass to remove tris (elements) that are too irregular (low area vs perimeter).
//  			// First get the edge-to-triangle map.
//  			SphereEdgeTriNumsMap.clear();
//  			for (int ti = 0; ti < SphereElems.size(); ++ti) {
//  				auto * t = &SphereElems[ti];
//  				for (int ci = 0; ci < 3; ++ci) {
//  					Edge e = MakeEdge((*t)[ci], (*t)[(ci + 1) % 3]);
//  					SphereEdgeTriNumsMap[e].insert(ti);
//  				}
//  			}
//  
//  			// Now see which elements will be removed
//  			double TriBadnessHighCheck = 0.45;
//  			vector<bool> RemoveElem(NumElems, false);
// 			int NumRemovedElems = 0;
//  			for (int ti = 0; ti < NumElems; ++ti) {
// 				if (TriBadnessHighCheck < TriBadness(SphereNodes[SphereElems[ti][0]], SphereNodes[SphereElems[ti][1]], SphereNodes[SphereElems[ti][2]])) {
// 					RemoveElem[ti] = true;
// 					NumRemovedElems++;
// 				}
//  			}
//  
// 			if (NumRemovedElems > 0) {
// 				// Now remove the bad elements. 
// 				// There are two cases of bad tris:
// 				// 1. A tri with one really short edge.
// 				//		In this case, we'll collapse the short edge, removing one node and the two elements that shared the edge,
// 				//		since the other sharing element must also be bad, and updating neighbor elements to reflect the removed
// 				//		node.
// 				//		If the edge to be removed is a constrained edge, then collapse it 
// 				// 2. A tri with one long and two shorter edges, such that the node between the shorter edges is very close to the 
// 				//		long edge. In this case we'll collapse the bad element itself, removing the node between the shorter edges
// 				//		and updating neighbor elements to reflect the node removal.
// 				// 
// 				vector<int> NodesOldToNew(NumNodes),
// 					ElemsOldToNew(NumElems);
// 				for (int ni = 0; ni < NumNodes; ++ni) NodesOldToNew[ni] = ni;
// 				double EdgeLengthToPerimeterLowLimit = 0.1;
// 				for (int ti = 0; ti < NumElems; ++ti) {
// 					if (ElemTypes[ti] < NETypeB && RemoveElem[ti] && SphereElems[ti].size() == 3
// 						&& NodesOldToNew[SphereElems[ti][0]] >= 0 && NodesOldToNew[SphereElems[ti][1]] >= 0 && NodesOldToNew[SphereElems[ti][2]] >= 0
// 						&& NodesOldToNew[SphereElems[ti][0]] != NodesOldToNew[SphereElems[ti][1]] && NodesOldToNew[SphereElems[ti][0]] != NodesOldToNew[SphereElems[ti][2]] && NodesOldToNew[SphereElems[ti][1]] != NodesOldToNew[SphereElems[ti][2]]
// 						&& TriBadnessHighCheck < TriBadness(SphereNodes[NodesOldToNew[SphereElems[ti][0]]], SphereNodes[NodesOldToNew[SphereElems[ti][1]]], SphereNodes[NodesOldToNew[SphereElems[ti][2]]])) {
// 						std::map<double, Edge> EdgeLengths;
// 						for (int ci = 0; ci < 3; ++ci) {
// 							int cj = (ci + 1) % 3;
// 							EdgeLengths[Distance(SphereNodes[SphereElems[ti][ci]], SphereNodes[SphereElems[ti][cj]])] = MakeEdge(ci, cj);
// 						}
// 						double Perimeter = TriPerimeter(SphereNodes[SphereElems[ti][0]], SphereNodes[SphereElems[ti][1]], SphereNodes[SphereElems[ti][2]]);
// 						double SmallEdgeLenOverPerimeter = EdgeLengths.cbegin()->first / Perimeter;
// 						// We have edge lengths sorted low to high in EdgeLengths.
// 						// If the ratio of the shortest edge length to Perimeter is low enough, 
// 						// Then it's a type 1 bad tri, otherwise it's a type 2.
// 						if (Perimeter > 0.0 && SmallEdgeLenOverPerimeter < EdgeLengthToPerimeterLowLimit) {
// 							// Type 1 bad triangle, and EdgeLengths.cbegin().second is the edge to be removed
// 							// Remove the two elements on this edge, then remove EdgeLengths.cbegin().second.second 
// 							// (the higher index node of the edge).
// 
// 							for (int tj : SphereEdgeTriNumsMap[EdgeLengths.cbegin()->second]) {
// 								SphereElems[tj].clear();
// 							}
// 							if (NodeTypes[EdgeLengths.cbegin()->second.second] == NETypeC)
// 								NodesOldToNew[EdgeLengths.cbegin()->second.second] = EdgeLengths.cbegin()->second.first;
// 							else
// 								NodesOldToNew[EdgeLengths.cbegin()->second.first] = EdgeLengths.cbegin()->second.second;
// 							EdgeIntSurfNums.erase(EdgeLengths.cbegin()->second);
// 						}
// 						else {
// 							// Type 2 bad tri or zero length perimeter.
// 							// Either way, remove the element and the node shared between the two
// 							// shorter edges in EdgeLenghts.
// 							auto e1 = EdgeLengths.cbegin(),
// 								e2 = e1;
// 							e2++;
// 							bool NodeFound = false;
// 							// 						if (!EdgeIntSurfNums.count(e1->second) && !EdgeIntSurfNums.count(e2->second)) {
// 							for (int ei = 0; ei < 2 && !NodeFound; ++ei) for (int ej = 0; ej < 2 && !NodeFound; ++ej) {
// 								int n1 = (ei == 0 ? e1->second.first : e1->second.second),
// 									n2 = (ej == 0 ? e2->second.first : e2->second.second);
// 								if (n1 == n2) {
// 									NodeFound = true;
// 									if (NodeTypes[n1] == NETypeC) {
// 										NodesOldToNew[n1] = -1;
// 										SphereElems[ti].clear();
// 									}
// 								}
// 							}
// 							// 						}
// 						}
// 					}
// 				}
// 
// 				// Now update everything to reflect the removals. This includes:
// 				// SphereElems, SphereNodes, NodeTypes, ElemTypes, NumElems, NumNodes,
// 				// OutConstrinedSegmentIndices
// 				// 
// 				vector<vec3> TmpNodes;
// 				TmpNodes.reserve(NumNodes);
// 				vector<vector<int> > TmpElems;
// 				vector<int> TmpNodesOldToNew(NumNodes, -1),
// 					TmpElemsOldToNew(NumElems, -1),
// 					TmpConstrainedSegmentIndices;
// 
// 				vector< GBATriangulatedSphereNodeElemType_e> TmpNodeTypes, TmpElemTypes;
// 				TmpNodeTypes.reserve(NumNodes);
// 				TmpElemTypes.reserve(NumElems);
// 				TmpConstrainedSegmentIndices.reserve(NumNodes);
// 
// 				for (int ni = 0; ni < NumNodes; ++ni) {
// 					if (NodesOldToNew[ni] == ni){
// 						TmpNodesOldToNew[ni] = TmpNodes.size();
// 						TmpNodes.push_back(SphereNodes[ni]);
// 						TmpNodeTypes.push_back(NodeTypes[ni]);
// 						TmpConstrainedSegmentIndices.push_back(OutConstrinedSegmentIndices[ni]);
// 					}
// 				}
// 
// 				for (int ti = 0; ti < NumElems; ++ti){
// 					if (!SphereElems[ti].empty()){
// 						auto e = &SphereElems[ti];
// 						for (int & ni : *e){
// 							ni = TmpNodesOldToNew[NodesOldToNew[ni]];
// 						}
// 						if (e->at(0) >= 0 && e->at(1) >= 0 && e->at(2) >= 0
// 							&& e->at(0) != e->at(1) && e->at(0) != e->at(2) && e->at(1) != e->at(2))
// 						{
// 							TmpElemsOldToNew[ti] = TmpElems.size();
// 							TmpElems.push_back(*e);
// 							TmpElemTypes.push_back(ElemTypes[ti]);
// 						}
// 					}
// 				}
// 
// 				TmpNodes.shrink_to_fit();
// 				TmpNodeTypes.shrink_to_fit();
// 				TmpElemTypes.shrink_to_fit();
// 				TmpConstrainedSegmentIndices.shrink_to_fit();
// 
// 				SphereNodes = TmpNodes;
// 				SphereElems = TmpElems;
// 
// 				NumNodes = SphereNodes.size();
// 				NumElems = SphereElems.size();
// 
// 				NodeTypes = TmpNodeTypes;
// 				ElemTypes = TmpElemTypes;
// 				OutConstrinedSegmentIndices = TmpConstrainedSegmentIndices;
// 
// #ifdef DEBUG_SAVEZONES
// 				TmpSphere = FESurface_c(SphereNodes, SphereElems);
// 				TmpSphere.SaveAsTriFEZone({ 1,2,3 }, "after triangulation update remove bad tris");
// #endif
// 			}
// 
//  		}

		// Store intersecting ring surface indices for each node (-1 if none)
		vector<int> NodeIntSurfNums(NumNodes, -1);
		// #pragma omp parallel for
		for (int ni = 0; ni < NumNodes; ++ni) {
			if (OutConstrainedSegmentIndices[ni] >= 0) {
				NodeIntSurfNums[ni] = constraintSegmentSurfNums[OutConstrainedSegmentIndices[ni]];
			}
		}



		// Perform element subdivision based on the minimum normalized arc length difference 
		// between the element midpoint and the nearest ring surface and bond path 
		// intersection points.
		// The normalization is according to the maximum bond path or ring surface intersection achieved by any element on the sphere.
		// The number of subdivisions is determined by a minimum element solid angle threshold specified by the user.
		// Elements near bond path intersections will have this solid angle when the subdivision is complete, 
		// while elements along ring surface intersections will have larger solid angles than this (less subdivision).
		// If an element is sufficiently far from a ring surface intersecting sphere node, then only its distance
		// to the nearest bond path intersecting node will be used.
		// In this way, the same amount of subdivision will occur around bond path intersections despite
		// the presence of ring surface intersections.
		// This fraction is determined by a linear combination of the
		// minimum normalized bond path and ring surface intersection arc length differences, d_b and d_r respectively.
		// Some parameters are defined to control the behavior of the linear combination:
		//	* d_r_c, double in (0,1), a cutoff value for d_r
		//	* w_b, double in (0,1), the weight given to d_b if d_r < d_r_c
		//	
		
 		 
		
		// Each time you subdivide the area of the new elements is 1/4 that of the original, so ElemSubdivisionLevel will
		// be incremented by 4 for edge-midpoint subdivision (1 elem -> 4 elements), and by 2 for the elements that are split as the result of a neighbor element's subdivision.
		vector<int> ElemSubdivisionLevel(NumElems, 0); 
		std::queue<int> ElemsToSubdivide;
		MaxSubdivisions = TecGUIScaleGetValue(SCBPGBInit_SC_T1_1);
		if (MaxSubdivisions > 0) {
			// First determine the maximum subdivision level (assuming uniform elements)
// 		int MaxSubdivisions = 6;// int(double(MaxNumElems) / double(NumElems));
			MaxSubdivisions++;
			vec DistCutoffVals = (LogSpace(1, pow(2, SubdivisionTightness+5), MaxSubdivisions + 1) - 1.0) / (pow(2, SubdivisionTightness+5) - 1);
			int MaxAreaSubdivisions = MaxSubdivisions * 4;

			double SphereArea = pow(Radius, 2.0) * 4.0 * PI;
			double MinElementArea = 1.0 / double(NumElems) * SphereArea / pow(4.0,MaxSubdivisions-1);
			double d_r_c = 0.3;
			double w_b = 0.2;

			// Outer loop MaxSubdivisions times
			for (int si = 0; si < MaxSubdivisions; ++si) {

				// Get maximum bond path and ring surface arc distances to any element
				double MaxBPDistance = -1.0;
				double MaxRSDistance = -1.0;
				double MinBPDistance = DBL_MAX;
				double MinRSDistance = DBL_MAX;
				for (int ti = 0; ti < NumElems; ++ti) {
					for (int ni = 0; ni < NumNodes; ++ni) {
						if (NodeTypes[ni] >= NETypeR) {
							double tmpDist = VectorAngle(SphereNodes[ni] - CPPos, TriangleElemMidPoint(SphereNodes, SphereElems, ti) - CPPos);
							if (NodeTypes[ni] >= NETypeB) {
								MaxBPDistance = MAX(MaxBPDistance, tmpDist);
								MinBPDistance = MIN(MinBPDistance, tmpDist);
							}
							else {
								MaxRSDistance = MAX(MaxRSDistance, tmpDist);
								MinRSDistance = MIN(MinRSDistance, tmpDist);
							}
						}
					}
				}
				double MaxArcDistance = MAX(MaxBPDistance, MaxRSDistance);

				if (MaxArcDistance > 0) {
					if (MaxRSDistance > 0.0){
						MaxRSDistance -= MinRSDistance;
					}
					MaxBPDistance -= MinBPDistance;
					// Now loop over elements, adding them to the subdivision queue if necessary
					// Get maximum bond path and ring surface arc distances to any element
					for (int ti = 0; ti < NumElems; ++ti) {
						// only look at element if its area isn't already too small
						double ElemArea = TriangleElemArea(SphereNodes, SphereElems, ti);
						if (ElemArea / MinElementArea > 2.0) {
							double d_b = DBL_MAX;
							double d_r = DBL_MAX;
							for (int ni = 0; ni < NumNodes; ++ni) {
								if (NodeTypes[ni] >= NETypeR) {
									double tmpDist = VectorAngle(SphereNodes[ni] - CPPos, TriangleElemMidPoint(SphereNodes, SphereElems, ti) - CPPos);
									if (NodeTypes[ni] >= NETypeB) {
										d_b = MIN(d_b, tmpDist);
									}
									else {
										d_r = MIN(d_r, tmpDist);
									}
								}
							}
							double d;
							d_b -= MinBPDistance;
							if (MaxRSDistance > 0.0){
								d_r -= MinRSDistance;
// 								d_b = (1.0 - (d_b / MaxBPDistance));
// 								d_r = (1.0 - (d_r / MaxRSDistance));
								d_b /= MaxBPDistance;
								d_r /= MaxRSDistance;

  								if (d_r > d_b){//d_r_c) {
  									d = d_b;
  								}
  								else {
  									d = w_b * d_b + (1.0 - w_b) * d_r;
  								}
//  								d = MIN(d_r, d_b);

// 								d = (1 - d_b) * d_b +  d_b * d_r;
							}
							else{
								d = d_b / MaxBPDistance;
							}


// 							int SubdivisionLevel = int(round(double(MaxSubdivisions) * d));
							int SubdivisionLevel = MaxSubdivisions;
							for (int i = 1; i < MaxSubdivisions; ++i){
								if (d <= DistCutoffVals[i]){
									SubdivisionLevel = i;
									break;
								}
							}
							SubdivisionLevel = (MaxSubdivisions - SubdivisionLevel) * 4;
							if (SubdivisionLevel > ElemSubdivisionLevel[ti]) {
								ElemsToSubdivide.push(ti);
							}
						}
					}

					if (ElemsToSubdivide.empty()) {
						break;
					}

					// subdivide any elements in the queue
					while (!ElemsToSubdivide.empty()) {
						if (!ElemsToSubdivide.empty()) {
							std::map<Edge, int> NewEdgeNodes;
							while (!ElemsToSubdivide.empty()) {
								int ti = ElemsToSubdivide.front();
								ElemSubdivisionLevel[ti] += 4;
								auto oldElem = SphereElems[ti];

								// Premake edge objects for element
								vector<Edge> ElemEdges(3);
								for (int ci = 0; ci < 3; ++ci)
									ElemEdges[ci] = MakeEdge(oldElem[ci], oldElem[(ci + 1) % 3]);

								// Make new nodes at midpoints of edges, using existing nodes if present.
								// Also add new node types.
								for (int ci = 0; ci < 3; ++ci) {
									if (!NewEdgeNodes.count(ElemEdges[ci])) {
										int NewNodeNum = SphereNodes.size();
										oldElem.push_back(NewNodeNum);
										NewEdgeNodes[ElemEdges[ci]] = NewNodeNum;
										SphereNodes.emplace_back((SphereNodes[ElemEdges[ci].first] + SphereNodes[ElemEdges[ci].second]) * 0.5);
										SphereNodes.back() = CPPos + normalise(SphereNodes.back() - CPPos) * Radius;
										if (EdgeIntSurfNums.count(ElemEdges[ci])) {
											NodeTypes.push_back(NETypeR);
											int EdgeIntSurfNum = EdgeIntSurfNums[ElemEdges[ci]];
											NodeIntSurfNums.push_back(EdgeIntSurfNum);
											EdgeIntSurfNums[MakeEdge(ElemEdges[ci].first, NewNodeNum)] = EdgeIntSurfNum;
											EdgeIntSurfNums[MakeEdge(ElemEdges[ci].second, NewNodeNum)] = EdgeIntSurfNum;
											EdgeIntSurfNums.erase(ElemEdges[ci]);
										}
										else {
											NodeTypes.push_back(NETypeC);
											NodeIntSurfNums.push_back(-1);
										}
									}
									else {
										oldElem.push_back(NewEdgeNodes[ElemEdges[ci]]);
									}
								}

								// Make new elements
								vector<int> newElems(4);
								for (int ei = 0; ei < 3; ++ei) {
									newElems[ei] = SphereElems.size();
									SphereElems.push_back({
										oldElem[ElemSubdivisionNewElemIndices[ei][0]],
										oldElem[ElemSubdivisionNewElemIndices[ei][1]],
										oldElem[ElemSubdivisionNewElemIndices[ei][2]]
										});
									ElemSubdivisionLevel.push_back(ElemSubdivisionLevel[ti]);
									ElemTypes.push_back(NETypeInvalid);
								}
								SphereElems[ti] = {
										oldElem[ElemSubdivisionNewElemIndices[3][0]],
										oldElem[ElemSubdivisionNewElemIndices[3][1]],
										oldElem[ElemSubdivisionNewElemIndices[3][2]]
								};
								newElems[3] = ti;
								ElemTypes[ti] = NETypeInvalid;
								for (int ei = 0; ei < 4; ++ei) {
									for (int ni : SphereElems[newElems[ei]])
										ElemTypes[newElems[ei]] = MAX(ElemTypes[newElems[ei]], NodeTypes[ni]);
								}

								ElemsToSubdivide.pop();
							}

							// Now update triangles that were not subdivided but share an edge with a triangle
							// that was subdivided. Loop over all elements and check if any of its edges are in
							// the NewEdgeNodes map. If an edge was split, then split the triangle in two
							// by creating a new edge to the existing node for the split edge.

							for (int ti = 0; ti < SphereElems.size(); ++ti) {
								bool TriSplit = false;
								for (int ci = 0; ci < 3 && !TriSplit; ++ci) {
									Edge e = MakeEdge(SphereElems[ti][ci], SphereElems[ti][(ci + 1) % 3]);
									if (NewEdgeNodes.count(e)) {
										// This edge was already split for the triangle on its other side.
										// Split the current triangle into two triangles by introducing a new
										// edge between the node at the center of the current edge and the
										// opposite corner of the triangle.
										int EdgeNodeNum = NewEdgeNodes[e];
										int OppositeCornerNode = SphereElems[ti][(ci + 2) % 3];
										SphereElems.push_back({ e.first, OppositeCornerNode, EdgeNodeNum });
										SphereElems[ti] = { e.second, EdgeNodeNum, OppositeCornerNode };

										ElemTypes.push_back(NETypeInvalid);
										ElemSubdivisionLevel[ti] += 2;
										ElemSubdivisionLevel.push_back(ElemSubdivisionLevel[ti]);
										ElemTypes[ti] = NETypeInvalid;

										for (int ei = 0; ei < 3; ++ei)
											ElemTypes[ti] = MAX(ElemTypes[ti], NodeTypes[SphereElems[ti][ei]]);
										for (int ei = 0; ei < 3; ++ei)
											ElemTypes.back() = MAX(ElemTypes.back(), NodeTypes[SphereElems.back()[ei]]);

										TriSplit = true;
									}
								}

								if (TriSplit) // remaining edges on triangle still need to be checked, so repeat this triangle
									ti--;
							}
						}



						NumNodes = SphereNodes.size();
						NumElems = SphereElems.size();

						// Jiggle Mesh
			// Do a few iterations of jiggle mesh to make the mesh
			// less ugly.
						vector<std::set<int> > NodeConnectivity;
						GetMeshNodeConnectivity(SphereElems, NumNodes, NodeConnectivity);

						vector<bool> NodeIsConstrained(NumNodes, false);
						for (int ni = 0; ni < NumNodes; ++ni) {
							if (NodeTypes[ni] > NETypeC) {
								NodeIsConstrained[ni] = true;
								if (MaxSubdivisions > 0 && NodeTypes[ni] >= NETypeB) {
									for (int nj : NodeConnectivity[ni])
										NodeIsConstrained[nj] = true;
								}
							}
						}

// 						while (TriangulatedSphereJiggleMesh(SphereNodes, NodeConnectivity, NodeIsConstrained, CPPos, Radius) > SphereJiggleMeshMaxMovedNodeDistTol) {}

						for (int i = 0; i < 10 && TriangulatedSphereJiggleMesh(SphereNodes, NodeConnectivity, NodeIsConstrained, CPPos, Radius) > SphereJiggleMeshMaxMovedNodeDistTol; ++i) {}
					}
				}
				
			}
		}

		//// Old explicit subdivision around bond points:
		//// 
		//// Perform midpoint subdivision for elements with bond path nodes before 
		//// Splitting them according to the angular constraint.
		//// This decreases the size of any resulting artifacts in the contours
		//// of condensed charge density.

		//std::queue<int> ElemsToSubdivide;
		//vector<int> ElemSubdivisionLevel(NumElems, 0);
		//for (int si = 0; si < NumberOfPreBondPathElemSubdivision; ++si) {
		//	for (int ti = 0; ti < NumElems; ++ti)
		//		if (ElemTypes[ti] >= NETypeB)
		//			ElemsToSubdivide.push(ti);

		//	if (!ElemsToSubdivide.empty()) {
		//		std::map<Edge, int> NewEdgeNodes;
		//		while (!ElemsToSubdivide.empty()) {
		//			int ti = ElemsToSubdivide.front();
		//			ElemSubdivisionLevel[ti] += 4;
		//			auto oldElem = SphereElems[ti];

		//			// Premake edge objects for element
		//			vector<Edge> ElemEdges(3);
		//			for (int ci = 0; ci < 3; ++ci)
		//				ElemEdges[ci] = MakeEdge(oldElem[ci], oldElem[(ci + 1) % 3]);

		//			// Make new nodes at midpoints of edges, using existing nodes if present.
		//			// Also add new node types.
		//			for (int ci = 0; ci < 3; ++ci) {
		//				if (!NewEdgeNodes.count(ElemEdges[ci])) {
		//					int NewNodeNum = SphereNodes.size();
		//					oldElem.push_back(NewNodeNum);
		//					NewEdgeNodes[ElemEdges[ci]] = NewNodeNum;
		//					SphereNodes.emplace_back((SphereNodes[ElemEdges[ci].first] + SphereNodes[ElemEdges[ci].second]) * 0.5);
		//					SphereNodes.back() = CPPos + normalise(SphereNodes.back() - CPPos) * Radius;
		//					if (EdgeIntSurfNums.count(ElemEdges[ci])) {
		//						NodeTypes.push_back(NETypeR);
		//						int EdgeIntSurfNum = EdgeIntSurfNums[ElemEdges[ci]];
		//						NodeIntSurfNums.push_back(EdgeIntSurfNum);
		//						EdgeIntSurfNums[MakeEdge(ElemEdges[ci].first, NewNodeNum)] = EdgeIntSurfNum;
		//						EdgeIntSurfNums[MakeEdge(ElemEdges[ci].second, NewNodeNum)] = EdgeIntSurfNum;
		//						EdgeIntSurfNums.erase(ElemEdges[ci]);
		//					}
		//					else {
		//						NodeTypes.push_back(NETypeC);
		//						NodeIntSurfNums.push_back(-1);
		//					}
		//				}
		//				else {
		//					oldElem.push_back(NewEdgeNodes[ElemEdges[ci]]);
		//				}
		//			}

		//			// Make new elements
		//			vector<int> newElems(4);
		//			for (int ei = 0; ei < 3; ++ei) {
		//				newElems[ei] = SphereElems.size();
		//				SphereElems.push_back({
		//					oldElem[ElemSubdivisionNewElemIndices[ei][0]],
		//					oldElem[ElemSubdivisionNewElemIndices[ei][1]],
		//					oldElem[ElemSubdivisionNewElemIndices[ei][2]]
		//					});
		//				ElemSubdivisionLevel.push_back(ElemSubdivisionLevel[ti]);
		//				ElemTypes.push_back(NETypeInvalid);
		//			}
		//			SphereElems[ti] = {
		//					oldElem[ElemSubdivisionNewElemIndices[3][0]],
		//					oldElem[ElemSubdivisionNewElemIndices[3][1]],
		//					oldElem[ElemSubdivisionNewElemIndices[3][2]]
		//			};
		//			newElems[3] = ti;
		//			ElemTypes[ti] = NETypeInvalid;
		//			for (int ei = 0; ei < 4; ++ei) {
		//				for (int ni : SphereElems[newElems[ei]])
		//					ElemTypes[newElems[ei]] = MAX(ElemTypes[newElems[ei]], NodeTypes[ni]);
		//			}

		//			ElemsToSubdivide.pop();
		//		}

		//		// Now update triangles that were not subdivided but share an edge with a triangle
		//		// that was subdivided. Loop over all elements and check if any of its edges are in
		//		// the NewEdgeNodes map. If an edge was split, then split the triangle in two
		//		// by creating a new edge to the existing node for the split edge.

		//		for (int ti = 0; ti < SphereElems.size(); ++ti) {
		//			bool TriSplit = false;
		//			for (int ci = 0; ci < 3 && !TriSplit; ++ci) {
		//				Edge e = MakeEdge(SphereElems[ti][ci], SphereElems[ti][(ci + 1) % 3]);
		//				if (NewEdgeNodes.count(e)) {
		//					// This edge was already split for the triangle on its other side.
		//					// Split the current triangle into two triangles by introducing a new
		//					// edge between the node at the center of the current edge and the
		//					// opposite corner of the triangle.
		//					int EdgeNodeNum = NewEdgeNodes[e];
		//					int OppositeCornerNode = SphereElems[ti][(ci + 2) % 3];
		//					SphereElems.push_back({ e.first, OppositeCornerNode, EdgeNodeNum });
		//					SphereElems[ti] = { e.second, EdgeNodeNum, OppositeCornerNode };

		//					ElemTypes.push_back(NETypeInvalid);
		//					ElemSubdivisionLevel[ti] += 2;
		//					ElemSubdivisionLevel.push_back(ElemSubdivisionLevel[ti]);
		//					ElemTypes[ti] = NETypeInvalid;

		//					for (int ei = 0; ei < 3; ++ei)
		//						ElemTypes[ti] = MAX(ElemTypes[ti], NodeTypes[SphereElems[ti][ei]]);
		//					for (int ei = 0; ei < 3; ++ei)
		//						ElemTypes.back() = MAX(ElemTypes.back(), NodeTypes[SphereElems.back()[ei]]);

		//					TriSplit = true;
		//				}
		//			}

		//			if (TriSplit) // remaining edges on triangle still need to be checked, so repeat this triangle
		//				ti--;
		//		}
		//	}



		//	NumNodes = SphereNodes.size();
		//	NumElems = SphereElems.size();
		//}

		for (int bi : IntBondPathNodes) {
			// Need elems around each bond path node, the edges of
			// those elements, and the first edge
			// neighbors of those elements.

			vector<int> NewElems;
			for (int ti = 0; ti < NumElems; ++ti) {
				for (int ci = 0; ci < 3; ++ci) {
					if (SphereElems[ti][ci] == bi) {
						// Check to see if and how many new nodes are necessary
						// by comparing the angle between element edges to the 
						// user-specified angle.

						Edge FarEdge = MakeEdge(SphereElems[ti][(ci + 1) % 3], SphereElems[ti][(ci + 2) % 3]);
						double TriAngle = VectorAngle(SphereNodes[FarEdge.first] - SphereNodes[SphereElems[ti][ci]], SphereNodes[FarEdge.second] - SphereNodes[SphereElems[ti][ci]]);
						int NumNewNodes = floorf(TriAngle / BPNodeCheckAngle);
						if (NumNewNodes > 0) {
							// Make new nodes along far edge opposite to bond path node
							double StepDist = Distance(SphereNodes[FarEdge.first], SphereNodes[FarEdge.second]) / double(NumNewNodes + 1);
							vec3 StepVec = normalise(SphereNodes[FarEdge.second] - SphereNodes[FarEdge.first]) * StepDist;
							vector<int> EdgeNodeNums = { FarEdge.first };
							for (int ni = 1; ni <= NumNewNodes; ++ni) {
								EdgeNodeNums.push_back(SphereNodes.size());
								SphereNodes.emplace_back(SphereNodes[FarEdge.first] + StepVec * (double)ni);
								SphereNodes.back() = CPPos + normalise(SphereNodes.back() - CPPos) * Radius;
							}
							EdgeNodeNums.push_back(FarEdge.second);

							ElemSubdivisionLevel[ti] += NumNewNodes;

							// Split bond path element
							for (int ni = 2; ni < EdgeNodeNums.size(); ++ni) {
								NewElems.push_back(SphereElems.size());
								SphereElems.push_back({ SphereElems[ti][ci], EdgeNodeNums[ni], EdgeNodeNums[ni - 1] });
								ElemSubdivisionLevel.push_back(ElemSubdivisionLevel[ti]);
							}
							SphereElems[ti] = { SphereElems[ti][ci], EdgeNodeNums[1], EdgeNodeNums[0] };
							NewElems.push_back(ti);

							// Now find element that shares the edge that was split
							// and split that element along the same new nodes.
							bool NeighborFound = false;
							for (int tj = 0; tj < NumElems && !NeighborFound; ++tj) {
								if (tj != ti) {
									for (int cj = 0; cj < 3; ++cj) {
										if (MakeEdge(SphereElems[tj][cj], SphereElems[tj][(cj + 1) % 3]) == FarEdge) {
											// Split neighbor element
											ElemSubdivisionLevel[tj] += NumNewNodes;
											int FarCorner = (cj + 2) % 3;
											for (int ni = 2; ni < EdgeNodeNums.size(); ++ni) {
												NewElems.push_back(SphereElems.size());
												SphereElems.push_back({ SphereElems[tj][FarCorner], EdgeNodeNums[ni], EdgeNodeNums[ni - 1] });
												ElemSubdivisionLevel.push_back(ElemSubdivisionLevel[tj]);
											}
											SphereElems[tj] = { SphereElems[tj][FarCorner], EdgeNodeNums[1], EdgeNodeNums[0] };
											NewElems.push_back(tj);
											NeighborFound = true;
											break;
										}
									}
								}
							}
						}
						break;
					}
				}
			}

			NumNodes = SphereNodes.size();
			NumElems = SphereElems.size();
			NodeTypes.resize(NumNodes, NETypeC);
			ElemTypes = vector<GBATriangulatedSphereNodeElemType_e>(NumElems, NETypeInvalid);

			for (int ti = 0; ti < NumElems; ++ti)
				for (int ni : SphereElems[ti])
					ElemTypes[ti] = MAX(ElemTypes[ti], NodeTypes[ni]);
		}

		

		NumNodes = SphereNodes.size();
		NumElems = SphereElems.size();
		OutConstrainedSegmentIndices.resize(NumNodes, -1);
		NodeIntSurfNums.resize(NumNodes, -1);

		vector<std::set<int> > NodeConnectivity;
		GetMeshNodeConnectivity(SphereElems, NumNodes, NodeConnectivity);

		vector<bool> NodeIsConstrained(NumNodes, false);
		for (int ni = 0; ni < NumNodes; ++ni) {
			if (NodeTypes[ni] > NETypeC) {
				NodeIsConstrained[ni] = true;
				if (MaxSubdivisions > 0 && NodeTypes[ni] >= NETypeB) {
					for (int nj : NodeConnectivity[ni])
						NodeIsConstrained[nj] = true;
				}
			}
		}

// 		while (TriangulatedSphereJiggleMesh(SphereNodes, NodeConnectivity, NodeIsConstrained, CPPos, Radius) > SphereJiggleMeshMaxMovedNodeDistTol) {}
		for (int i = 0; i < 10 && TriangulatedSphereJiggleMesh(SphereNodes, NodeConnectivity, NodeIsConstrained, CPPos, Radius) > SphereJiggleMeshMaxMovedNodeDistTol; ++i) {}

		if (TestRun) {
			auto TmpSphere = FESurface_c(SphereNodes, SphereElems);
			TmpSphere.SaveAsTriFEZone({ 1,2,3 }, "after bond path intersection element subdivision");
			TecUtilZoneSetActive(Set(TecUtilDataSetGetNumZones()).getRef(), AssignOp_PlusEquals);
			continue;
// 			StatusDrop(AddOnID);
// 			return;
		}

#ifdef DEBUG_SAVEZONES
		TmpSphere = FESurface_c(SphereNodes, SphereElems);
		TmpSphere.SaveAsTriFEZone({ 1,2,3 }, "after bond path intersection element subdivision");
		TecUtilZoneSetActive(Set(TecUtilDataSetGetNumZones()).getRef(), AssignOp_MinusEquals);
#endif

		bool RingPresent = false;
		for (int i = 0; i < NodeTypes.size() && !RingPresent; ++i){
			RingPresent = (NodeTypes[i] == NETypeR);
		}

		/*
		 *	Setup and seed all type C and R paths.
		 *	R paths will be combined with the ring path segments found above.
		 *	We need R and C paths in order to get B and RB paths.
		 */
		vector<GradPath_c> GradPaths;
		vector<GradPath_c> InnerGradPaths;

		/*
		 *	GPs are also seeded along the edges of triangular elements whose
		 *	node GPs experience too great a separation at some point down the GPs.
		 *	Keep the edge GPs inside GradPaths (above) and maintain indices of the edge
		 *	GPs in a <edge, vector<int>> map.
		 */
		std::map<Edge, vector<int> > EdgeGPMap, RingEdgeGPMap;

		// Start the adaptive GBA loop.
		// This will continue until there are no more gradient bundles that
		// need to be subdivided.
		// 
		vector<int> DivergentGBs,
			NotMonotonicallyDecreasingGPs,
			NotMadeGPs,
			InvalidGBs;

		vector<vector<GradPath_c const*> > GPPtrList;
		vector<vector<int> > SphereElemGPInds(SphereElems);

		int GBAIter = 0;
		double TotalDensity = 0.0;
		int NumElementsToSubdivide = -1;

		vector<int> ElemTodo;

		do
		{
			NumCompleted = 0;
			GBAIter++;


#ifdef SUPERDEBUG
// 			ElemTodo = { 196 }; // base 0 (tecplot is base 1)
  			ElemTodo = { 0 };
//         			ElemTodo.clear();
//  // 				for (int i = 0; i < NumElems; ++i) {
// 			for (int i = 0; i < NumElems; ++i) {
// 				if (ElemTypes[i] >= NETypeR)
// 					ElemTodo.push_back(i);
// 			}

//				// Get ring and bond type elements and their n-th nearest neighbors
//   				int MaxBFSDepth = 4; // if 3, gets ring/bond elems and their 4 nearest neighbors
//   				vector<vector<int> > SphereElemConnectivity;
//   				GetTriElementConnectivityList(&SphereElems, SphereElemConnectivity, 1);
//   				std::set<int> ElemTodoSet;
//   				auto CompNEType = NETypeR;
//   
//   				for (int i = 0; i < NumElems; ++i) {
//   					if (ElemTypes[i] >= CompNEType) {
//   						ElemTodoSet.insert(i);
//   						// BFS with depth limit
//   						std::queue<int> ToCheck;
//   						vector<bool> Visited(SphereElems.size(), false);
//   						int CurDepth = 0;
//   						ToCheck.push(i);
//   						ToCheck.push(-1); //Depth level marker
//   						Visited[i] = true;
//   						while (!ToCheck.empty()) {
//   							if (ToCheck.front() == -1) {
//   								// Moving down another depth
//   								ToCheck.push(-1);
//   								if (++CurDepth > MaxBFSDepth)
//   									break;
//   							}
//   							else {
//   								int e = ToCheck.front();
//   // 								if (ElemTypes[e] >= CompNEType) {
//   									ElemTodoSet.insert(e);
//   // 								}
//   								for (int n : SphereElemConnectivity[e]) {
//   									if (!Visited[n]) {
//   										ToCheck.push(n);
//   										Visited[n] = true;
//   									}
//   								}
//   							}
//   							ToCheck.pop();
//   						}
//   					}
//   				}
//   
//   				ElemTodo = vector<int>(ElemTodoSet.begin(), ElemTodoSet.end());
   					
//   			ElemTodo.clear();
//   			for (int ti = 0; ti < NumElems; ++ti) {
//   				if (ElemTypes[ti] >= NETypeB) {
// 					for (int i = 0; i < 3; ++i) {
// 						if (NodeTypes[SphereElems[ti][i]] == NETypeR) {
// 							ElemTodo.push_back(ti);
// 							break;
// 						}
// 					}
//   				}
//   			}

#endif

			// This sizes GradPaths for the first iteration and 
			// destroys the composite paths (for R and B type nodes)
			// from the previous iteration so that they can be found
			// again for the new nodes/elements.
			// This probably results in remaking some composite GPs
			// for elements that weren't subdivided, but since most
			// subdivision will occur closer to R and B type nodes,
			// this shouldn't be a lot of wasted work.
			GradPaths.resize(NumNodes);

			// Subdivide any elements in the queue
			if (!ElemsToSubdivide.empty()) {
				std::map<Edge,int> NewEdgeNodes;
				while (!ElemsToSubdivide.empty()) {
					int ti = ElemsToSubdivide.front();
					ElemsToSubdivide.pop();
					if (ElemSubdivisionLevel[ti] >= MaxGBSubdivisionLevel)
						continue;
					ElemSubdivisionLevel[ti] += 4;

					auto oldElem = SphereElems[ti];

					// Premake edge objects for element
					vector<Edge> ElemEdges(3);
					for (int ci = 0; ci < 3; ++ci)
						ElemEdges[ci] = MakeEdge(oldElem[ci], oldElem[(ci + 1) % 3]);

					// Make new nodes at midpoints of edges, using existing nodes if present.
					// Also add new node types.
					for (int ci = 0; ci < 3; ++ci){
						if (!NewEdgeNodes.count(ElemEdges[ci])) {
							int NewNodeNum = SphereNodes.size();
							oldElem.push_back(NewNodeNum);
							NewEdgeNodes[ElemEdges[ci]] = NewNodeNum;
							SphereNodes.emplace_back((SphereNodes[ElemEdges[ci].first] + SphereNodes[ElemEdges[ci].second]) * 0.5);
							SphereNodes.back() = CPPos + normalise(SphereNodes.back() - CPPos) * Radius;
							if (EdgeIntSurfNums.count(ElemEdges[ci])) {
								NodeTypes.push_back(NETypeR);
								int EdgeIntSurfNum = EdgeIntSurfNums[ElemEdges[ci]];
								NodeIntSurfNums.push_back(EdgeIntSurfNum);
								EdgeIntSurfNums[MakeEdge(ElemEdges[ci].first, NewNodeNum)] = EdgeIntSurfNum;
								EdgeIntSurfNums[MakeEdge(ElemEdges[ci].second, NewNodeNum)] = EdgeIntSurfNum;
								EdgeIntSurfNums.erase(ElemEdges[ci]);
							}
							else {
								NodeTypes.push_back(NETypeC);
								NodeIntSurfNums.push_back(-1);
							}
						}
						else{
							oldElem.push_back(NewEdgeNodes[ElemEdges[ci]]);
						}
					}

					// Make new elements
					vector<int> newElems(4);
					for (int ei = 0; ei < 3; ++ei) {
						newElems[ei] = SphereElems.size();
						SphereElems.push_back({ 
							oldElem[ElemSubdivisionNewElemIndices[ei][0]], 
							oldElem[ElemSubdivisionNewElemIndices[ei][1]],
							oldElem[ElemSubdivisionNewElemIndices[ei][2]] 
							});
						ElemSubdivisionLevel.push_back(ElemSubdivisionLevel[ti]);
						ElemTypes.push_back(NETypeInvalid);
					}
					SphereElems[ti] = {
							oldElem[ElemSubdivisionNewElemIndices[3][0]],
							oldElem[ElemSubdivisionNewElemIndices[3][1]],
							oldElem[ElemSubdivisionNewElemIndices[3][2]]
						};
					newElems[3] = ti;
					ElemTypes[ti] = NETypeInvalid;
					IntVals[ti] = vector<double>(IntVarNumList.size() + 1, 0.0);
					for (int ei = 0; ei < 4; ++ei) {
						for (int ni : SphereElems[newElems[ei]])
							ElemTypes[newElems[ei]] = MAX(ElemTypes[newElems[ei]], NodeTypes[ni]);
					}
				}

				// Now update triangles that were not subdivided but share an edge with a triangle
				// that was subdivided. Loop over all elements and check if any of its edges are in
				// the NewEdgeNodes map. If an edge was split, then split the triangle in two
				// by creating a new edge to the existing node for the split edge.

				IntVals.resize(SphereElems.size(), vector<double>(IntVarNumList.size() + 1, 0.0));
				
				for (int ti = 0; ti < SphereElems.size(); ++ti){
					bool TriSplit = false;
					for (int ci = 0; ci < 3 && !TriSplit; ++ci){
						Edge e = MakeEdge(SphereElems[ti][ci], SphereElems[ti][(ci + 1) % 3]);
						if (NewEdgeNodes.count(e)){
							// This edge was already split for the triangle on its other side.
							// Split the current triangle into two triangles by introducing a new
							// edge between the node at the center of the current edge and the
							// opposite corner of the triangle.
							ElemSubdivisionLevel[ti] += 2;
							int EdgeNodeNum = NewEdgeNodes[e];
							int OppositeCornerNode = SphereElems[ti][(ci + 2) % 3];
							SphereElems.push_back({ e.first, OppositeCornerNode, EdgeNodeNum });
							ElemSubdivisionLevel.push_back(ElemSubdivisionLevel[ti]);
							SphereElems[ti] = { e.second, EdgeNodeNum, OppositeCornerNode };
							
							ElemTypes.push_back(NETypeInvalid);
							IntVals.emplace_back(IntVarNumList.size() + 1, 0.0);
							ElemTypes[ti] = NETypeInvalid;
							IntVals[ti] = vector<double>(IntVarNumList.size() + 1, 0.0);
							
							for (int ei = 0; ei < 3; ++ei)
								ElemTypes[ti] = MAX(ElemTypes[ti], NodeTypes[SphereElems[ti][ei]]);
							for (int ei = 0; ei < 3; ++ei)
								ElemTypes.back() = MAX(ElemTypes.back(), NodeTypes[SphereElems.back()[ei]]);

							TriSplit = true;
						}
					}

					if (TriSplit) // remaining edges on triangle still need to be checked, so repeat this triangle
						ti--;
				}
			}



			NumNodes = SphereNodes.size();
			NumElems = SphereElems.size();

			// Jiggle Mesh
			// Do a few iterations of jiggle mesh to make the mesh
			// less ugly.
			vector<std::set<int> > NodeConnectivity;
			GetMeshNodeConnectivity(SphereElems, NumNodes, NodeConnectivity);

 			vector<bool> NodeIsConstrained(NumNodes, false);
 			for (int ni = 0; ni < NumNodes; ++ni) {
 				if (NodeTypes[ni] > NETypeC){
 					NodeIsConstrained[ni] = true;
 // 					if (NodeTypes[ni] >= NETypeB){
 // 						for (int nj : NodeConnectivity[ni])
 // 							NodeIsConstrained[nj] = true;
 // 					}
 				}
 			}
 
//  			while(TriangulatedSphereJiggleMesh(SphereNodes, NodeConnectivity, NodeIsConstrained, CPPos, Radius) > SphereJiggleMeshMaxMovedNodeDistTol){}
			for (int i = 0; i < 10 && TriangulatedSphereJiggleMesh(SphereNodes, NodeConnectivity, NodeIsConstrained, CPPos, Radius) > SphereJiggleMeshMaxMovedNodeDistTol; ++i) {}


#ifdef DEBUG_SAVEZONES
 			TmpSphere = FESurface_c(SphereNodes, SphereElems);
			TmpSphere.SaveAsTriFEZone({ 1,2,3 }, "after triangulation update jiggle mesh");
			TecUtilZoneSetActive(Set(TecUtilDataSetGetNumZones()).getRef(), AssignOp_MinusEquals);

			vector<string> s = { "B", "R", "C" };
			vector<ColorIndex_t> c = { Red_C, Green_C, Blue_C };
			vector<GBATriangulatedSphereNodeElemType_e> TypeList = { NETypeB, NETypeR, NETypeC };
			for (int i = 0; i < 3; ++i) {
				vector<vec3> tmpVec;
				for (int ni = 0; ni < NumNodes; ++ni) {
					if (NodeTypes[ni] == TypeList[i])
						tmpVec.push_back(SphereNodes[ni]);
				}
				SaveVec3VecAsScatterZone(tmpVec, s[i] + "-type nodes", c[i]);
				TecUtilZoneSetActive(Set(TecUtilDataSetGetNumZones()).getRef(), AssignOp_MinusEquals);

				tmpVec.clear();
				for (int ti = 0; ti < NumElems; ++ti) {
					if (ElemTypes[ti] == TypeList[i]) {
						vec3 midPt = zeros(3);
						for (int ci = 0; ci < 3; ++ci) {
							midPt += SphereNodes[SphereElems[ti][ci]];
						}
						tmpVec.push_back(midPt / 3.);
					}
				}
				SaveVec3VecAsScatterZone(tmpVec, s[i] + "-type elems", c[i]);
				TecUtilZoneSetActive(Set(TecUtilDataSetGetNumZones()).getRef(), AssignOp_MinusEquals);
			}

			// 		TecUtilDataLoadEnd();
			// 		StatusDrop(AddOnID);
			// 		return;
#endif

  			GradPaths.resize(NumNodes);
  			IntVals.resize(NumElems, vector<double>(IntVarNumList.size() + 1, 0.0));
			if (RadialSphereApprx) {
				SphereIntVals.resize(NumElems, vector<double>(IntVarNumList.size() + 1, 0.0));
				SphereTotalIntVals.resize(IntVarNumList.size() + 1, 0.0);
			}
			
// 			GradPaths = vector<GradPath_c>(NumNodes);
// 			IntVals = vector<vector<double> >(NumElems, vector<double>(IntVarNumList.size() + 1, 0.0));

			OutConstrainedSegmentIndices.resize(NumNodes, -1);
			NodeIntSurfNums.resize(NumNodes, -1);


			SphereElemGPInds = SphereElems;

			stringstream ProgressStr2;

			ProgressStr2 << ProgressStr1.str() << " Pass " << GBAIter;

			if (GBAIter > 1)
				ProgressStr2 << " (" << TotalDensity << " e): ";
			else
				ProgressStr2 << ": ";

			Time1 = high_resolution_clock::now();
			TmpString = ProgressStr2.str() + "Make GPs";

			int NumToDo = 0;
			for (int ni = 0; ni < NumNodes; ++ni)
				NumToDo += (int)!GradPaths[ni].IsMade();
			NumCompleted = 0;

			/*
			 *	Changing it so grad paths don't go into the sphere
			 */

// 			if (RadialSphereApprx)
				InnerGradPaths.resize(NumNodes);

#ifdef SUPERDEBUG
			std::set<int> NodesTodo;
			for (int ti : ElemTodo) {
				for (int ni : SphereElems[ti])
					NodesTodo.insert(ni);
			}
			vector<int> NodesTodoVec(NodesTodo.begin(), NodesTodo.end());
			NumToDo = NodesTodoVec.size();

// #pragma omp parallel for schedule(dynamic)
// 			for (int nii = 0; nii < NumToDo; ++nii) {
// 				int ni = NodesTodoVec[nii];
// #else
// #pragma omp parallel for schedule(dynamic)
// 			for (int ni = 0; ni < NumNodes; ++ni) {
#endif
// 				if (!GradPaths[ni].IsMade()) {
// 					if (!RadialSphereApprx) {
// 						InnerGradPaths[ni].SetupGradPath(SphereNodes[ni], StreamDir_Forward, NumGPPoints * 2, GPType_Classic, GPTerminate_AtCP, nullptr, &CPs, &GPNCPTermRadius, &CutoffVal, ThVolInfo[0], HessPtrs, GradPtrs, RhoPtr);
// 						InnerGradPaths[ni].SetTerminalCPTypeNum(CPTypeNum_Nuclear);
// // 					}
// 					if (NodeTypes[ni] == NETypeC) {
// 						GradPaths[ni].SetupGradPath(SphereNodes[ni], StreamDir_Reverse, NumGPPoints * 2, GPType_Classic, GPTerminate_AtCP, nullptr, &CPs, &GPTermRadius, &CutoffVal, ThVolInfo[0], HessPtrs, GradPtrs, RhoPtr);
// 						GradPaths[ni].SetTerminalCPTypeNum(CPTypeNum_Cage);
// 					}
// 					else if (NodeTypes[ni] == NETypeR) {
// 						GradPaths[ni].SetupGradPath(SphereNodes[ni], StreamDir_Reverse, NumGPPoints * 2, GPType_Classic, GPTerminate_AtCP, nullptr, &CPs, &RTypeTermRadius, &CutoffVal, ThVolInfo[0], HessPtrs, GradPtrs, RhoPtr, &IntSurfs[NodeIntSurfNums[ni]]);
// 						GradPaths[ni].SetTerminalCPTypeNum(CPTypeNum_Ring);
// 					}
// 					else continue;
// 				}
// 				if (!GradPaths[ni].IsMade()) {
// 					if (NodeTypes[ni] == NETypeC) {
// 						if (RadialSphereApprx)
// 							GradPaths[ni].SetupGradPath(SphereNodes[ni], StreamDir_Reverse, NumGPPoints, GPType_Classic, GPTerminate_AtCP, nullptr, &CPs, &GPTermRadius, &CutoffVal, ThVolInfo[0], HessPtrs, GradPtrs, RhoPtr);
// 						else
// 							GradPaths[ni].SetupGradPath(SphereNodes[ni], StreamDir_Both, NumGPPoints, GPType_Classic, GPTerminate_AtCP, nullptr, &CPs, &GPTermRadius, &CutoffVal, ThVolInfo[0], HessPtrs, GradPtrs, RhoPtr);
// 
// 						GradPaths[ni].SetTerminalCPTypeNum(CPTypeNum_Cage);
// 					}
// 					else if (NodeTypes[ni] == NETypeR) {
// 						GradPaths[ni].SetupGradPath(SphereNodes[ni], StreamDir_Reverse, NumGPPoints, GPType_Classic, GPTerminate_AtCP, nullptr, &CPs, &RTypeTermRadius, &CutoffVal, ThVolInfo[0], HessPtrs, GradPtrs, RhoPtr, &IntSurfs[NodeIntSurfNums[ni]]);
// 						GradPaths[ni].SetTerminalCPTypeNum(CPTypeNum_Ring);
// 
// 						if (!RadialSphereApprx) {
// 							InnerGradPaths[ni].SetupGradPath(SphereNodes[ni], StreamDir_Forward, NumGPPoints, GPType_Classic, GPTerminate_AtCP, nullptr, &CPs, &GPTermRadius, &CutoffVal, ThVolInfo[0], HessPtrs, GradPtrs, RhoPtr);
// 							InnerGradPaths[ni].SetTerminalCPTypeNum(CPTypeNum_Nuclear);
// 						}
// 					}
// 					else continue;
// 				}
// 			}

	#ifdef SUPERDEBUG
#pragma omp parallel for schedule(dynamic)
			for (int nii = 0; nii < NodesTodoVec.size(); ++nii) {
				int ni = NodesTodoVec[nii];
	#else
	 #pragma omp parallel for schedule(dynamic)
 			for (int ni = 0; ni < NumNodes; ++ni) {
	#endif
				int ThNum = omp_get_thread_num();
				if (ThNum == 0) {
					UserQuit = !StatusUpdate(NumCompleted, NumToDo, TmpString, AddOnID, Time1);
	#pragma omp flush (UserQuit)
				}
	#pragma omp flush (UserQuit)
				if (UserQuit) {
#pragma omp critical
					{
						if (ThNum == 0) {
							if (TecUtilDialogMessageBox("Cancel run? All progress will be lost if you press \"Yes.\"", MessageBoxType_YesNo)) {
								UserQuit = true;
							}
							else {
								StatusLaunch(TmpString, AddOnID, TRUE, TRUE, TRUE);
								StatusUpdate(NumCompleted, NumToDo, TmpString, AddOnID, Time1);
								UserQuit = false;
							}
// 						}
// 						if (omp_get_thread_num() == 0){
// 							if (!TecUtilDialogMessageBox("Cancel run? All progress will be lost if you press \"Yes.\"", MessageBoxType_YesNo)) {
// // 								StatusLaunch(TmpString, AddOnID, TRUE, TRUE, TRUE);
// // 								StatusUpdate(NumCompleted, NumToDo, TmpString, AddOnID, Time1);
// 								UserQuit = false;
// 							}
#pragma omp flush (UserQuit)
						}
#pragma omp flush (UserQuit)
					}
				}

				if (!UserQuit && !GradPaths[ni].IsMade()) {
#pragma omp atomic
					NumCompleted++;
					// 					if (!RadialSphereApprx) {
					InnerGradPaths[ni].SetupGradPath(SphereNodes[ni], StreamDir_Forward, NumGPPoints * 2, GPType_Classic, GPTerminate_AtCP, nullptr, &CPs, &GPNCPTermRadius, &CutoffVal, ThVolInfo[0], HessPtrs, GradPtrs, RhoPtr);
					InnerGradPaths[ni].SetTerminalCPTypeNum(CPTypeNum_Nuclear);
					// 					}
					if (NodeTypes[ni] == NETypeC) {
						GradPaths[ni].SetupGradPath(SphereNodes[ni], StreamDir_Reverse, NumGPPoints * 2, GPType_Classic, GPTerminate_AtCP, nullptr, &CPs, &GPTermRadius, &CutoffVal, ThVolInfo[0], HessPtrs, GradPtrs, RhoPtr);
						GradPaths[ni].SetTerminalCPTypeNum(CPTypeNum_Cage);
					}
					else if (NodeTypes[ni] == NETypeR) {
						GradPaths[ni].SetupGradPath(SphereNodes[ni], StreamDir_Reverse, NumGPPoints * 2, GPType_Classic, GPTerminate_AtCP, nullptr, &CPs, &RTypeTermRadius, &CutoffVal, ThVolInfo[0], HessPtrs, GradPtrs, RhoPtr, &IntSurfs[NodeIntSurfNums[ni]]);
						GradPaths[ni].SetTerminalCPTypeNum(CPTypeNum_Ring);
					}


// 					if (RadialSphereApprx){
// 						if (NodeTypes[ni] <= NETypeR) {
// 							GradPaths[ni].Seed();
// 						}
// 						else continue;
  					InnerGradPaths[ni].Seed(false);
  					InnerGradPaths[ni].RemoveKinks();
  					if (RadialSphereApprx){// && InnerGradPaths[ni].IsMade() && InnerGradPaths[ni].RhoAt(0) < InnerGradPaths[ni].RhoAt(-1)) {
  						InnerGradPaths[ni].Reverse();
  					}

						// just do straight lines for inner grad paths
// 					vector<vec3> Pts(NumGPPoints);
// 					vec3 StepVec = (CPPos - SphereNodes[ni]) / (double(NumGPPoints-1));
// 					for (int i = 0; i < NumGPPoints; ++i){
// 						Pts[i] = SphereNodes[ni] + StepVec * double(i);
// 					}
// 					InnerGradPaths[ni] = GradPath_c(Pts);
// 					InnerGradPaths[ni].ReinterpolateRhoValuesFromVolume(&ThVolInfo[omp_get_thread_num()], &RhoPtr);

// 					if (RadialSphereApprx) {// && InnerGradPaths[ni].IsMade() && InnerGradPaths[ni].RhoAt(0) < InnerGradPaths[ni].RhoAt(-1)) {
// 						InnerGradPaths[ni].Reverse();
// 					}

// 					}
// 					else {
						if (NodeTypes[ni] <= NETypeR) {
							GradPaths[ni].Seed(RadialSphereApprx);
							GradPaths[ni].RemoveKinks();
							
							if (!RadialSphereApprx) {
								GradPaths[ni] = ConcatenateResample({ GradPaths[ni],InnerGradPaths[ni] }, NumGPPoints * 2, {}, ResampleMethod);
							}
						}
						else continue;
// 						if (NodeTypes[ni] == NETypeC) {
// 							GradPaths[ni].Seed();
// 						}
// 						else if (NodeTypes[ni] == NETypeR) {
// 							GradPaths[ni].Seed();
// 							GradPaths[ni].RemoveKinks();
// 							InnerGradPaths[ni].Seed();
// 							GradPaths[ni] = ConcatenateResample({ GradPaths[ni],InnerGradPaths[ni].SubGP(1, -1) }, NumGPPoints, {}, ResampleMethod);
// 						}
// 						else continue;
// 					}

					if (GradPaths[ni].IsMade() && GradPaths[ni].RhoAt(0) < GradPaths[ni].RhoAt(-1)) GradPaths[ni].Reverse();
				}
 			}

// 			int NumGPPointsOuter;
			int NumGPPointsSphereInterior = NumGPPoints / 10;
// 			NumGPPoints = (NumGPPoints * 9) / 10;
// 			int NumGPPointsPathSegment = NumGPPoints / (RingPresent ? 3 : 2);

			if (RadialSphereApprx) {
				/*
				 *	Resample inner and outer grad paths as if they were connected, using the lengths of the outer grad
				 *	paths to determine the number of points for their corresponding inner grad paths.
				 *	Outer path lengths should be
				 */

// 				double AvgOuterGradPathLength = 0.0, AvgInnerGradPathLength = 0.0, MaxOuterGradPathLength = DBL_MIN;
// 				int NumMadeGradPaths = 0, NumMadeInnerGradPaths = 0;
// 
// 				for (int i = 0; i < GradPaths.size(); ++i) {
// 					if (GradPaths[i].IsMade()) {
// 						NumMadeGradPaths++;
// 						double OuterLength = GradPaths[i].GetLength();
// 						if (OuterLength > MaxOuterGradPathLength) {
// 							MaxOuterGradPathLength = OuterLength;
// 						}
// 						AvgOuterGradPathLength += OuterLength;
// 					}
// 				}
// 				for (int i = 0; i < InnerGradPaths.size(); ++i) {
// 					if (InnerGradPaths[i].IsMade()) {
// 						AvgInnerGradPathLength += InnerGradPaths[i].GetLength();
// 						NumMadeInnerGradPaths++;
// 					}
// 				}
// 
// 				REQUIRE(NumMadeGradPaths > 0);
// 				REQUIRE(NumMadeInnerGradPaths > 0);
// 
// 				AvgOuterGradPathLength /= double(NumMadeGradPaths);
// 				AvgInnerGradPathLength /= double(NumMadeInnerGradPaths);
// 				// 				NumGPPointsOuter = int(AvgOuterGradPathLength / (AvgInnerGradPathLength + AvgOuterGradPathLength) * double(NumGPPoints));
// 				NumGPPointsOuter = int(MaxOuterGradPathLength / (AvgInnerGradPathLength + MaxOuterGradPathLength) * double(NumGPPoints));
// 				int InnerPoints = NumGPPoints - NumGPPointsOuter;
				

				for (int i = 0; i < GradPaths.size(); ++i) {
					if (GradPaths[i].IsMade()) {
						GradPaths[i].Resample(NumGPPoints, ResampleMethod);
					}
				}
				for (int i = 0; i < InnerGradPaths.size(); ++i) {
					if (InnerGradPaths[i].IsMade()) {
						InnerGradPaths[i].Resample(NumGPPointsSphereInterior, ResampleMethod);
					}
				}

				// 				int MinInnerGPNumPoints = INT_MAX;
				// 
				// 				for (const auto &  g : InnerGradPaths){
				// 					if (g.IsMade()) {
				// 						MinInnerGPNumPoints = MIN(MinInnerGPNumPoints, g.GetCount());
				// 					}
				// 				}
				// 				for (auto & g : InnerGradPaths){
				// 					if (g.IsMade() && g.GetCount() > MinInnerGPNumPoints){
				// 						g.Resample(MinInnerGPNumPoints, ResampleMethod);
				// 					}
				// 				}
			}
			else {
				for (int i = 0; i < GradPaths.size(); ++i) {
					if (GradPaths[i].IsMade()) {
						GradPaths[i].Resample(NumGPPoints, ResampleMethod);
					}
				}
			}
			

			if (UserQuit) {
				TecUtilDataLoadEnd();
				StatusDrop(AddOnID);
				return;
			}
			Time1 = high_resolution_clock::now();
			TmpString = ProgressStr2.str() + "Make ring surface GPs";

#ifdef SUPERDEBUG
			NumToDo = ElemTodo.size();
			NumCompleted = 0;
			for (int ti : ElemTodo) {
				UserQuit = !StatusUpdate(NumCompleted++, NumToDo, TmpString, AddOnID, Time1);
				if (UserQuit){
					if (TecUtilDialogMessageBox("Cancel run? All progress will be lost if you press \"Yes.\"", MessageBoxType_YesNo)) {
						break;
					}
					else {
						StatusLaunch(TmpString, AddOnID, TRUE, TRUE, TRUE);
						StatusUpdate(NumCompleted, NumToDo, TmpString, AddOnID, Time1);
						UserQuit = false;
					}
				}
#else
			NumToDo = 0;
			for (int ti = 0; ti < NumElems; ++ti)
				NumToDo += int(ElemTypes[ti] >= NETypeR && IntVals[ti][DensityVarNumInIntList] == 0.0);
			NumCompleted = 0;

			for (int ti = 0; ti < NumElems && !UserQuit; ++ti) {
				UserQuit = !StatusUpdate(NumCompleted, NumToDo, TmpString, AddOnID, Time1);
				if (UserQuit){
					if (TecUtilDialogMessageBox("Cancel run? All progress will be lost if you press \"Yes.\"", MessageBoxType_YesNo)) {
						break;
					}
					else {
						StatusLaunch(TmpString, AddOnID, TRUE, TRUE, TRUE);
						StatusUpdate(NumCompleted, NumToDo, TmpString, AddOnID, Time1);
						UserQuit = false;
					}
				}
				if (ElemTypes[ti] >= NETypeR && IntVals[ti][DensityVarNumInIntList] == 0.0)
					NumCompleted++;
#endif
				if (IntVals[ti][DensityVarNumInIntList] == 0.0) {
					/*
						*	For each corner, check if it's type R or B.
						*	If type R then just need to update the triangle
						*	to point to the correct R path.
						*	If B, need to construct GPs to finish the paths,
						*	using information from the neighboring GPs.
						*/
					auto * tri = &SphereElemGPInds[ti];
					for (int c1 = 0; c1 < 3; ++c1) {
						/*
							*	First pass over triangle to check for R type nodes
							*	to update them to the correct side.
							*	There could be 0, 1, or 2 R type nodes on the triangle.
							*	If 1, can use the guaranteed C type node to find the correct
							*	ring path, by comparing distances of the two possible ring path
							*	segment terminuses to that of the C type node path.
							*	If 2, can check both ring paths for each node to find the pair
							*  whose end points are closest together.
							*  Also if 2 then we'll do them both together, so can break after
							*  they're found.
							*/


						if (NodeTypes[SphereElems[ti][c1]] == NETypeR && SphereElems[ti][c1] == tri->at(c1)) {

							if (NodeTypes[SphereElems[ti][(c1 + 1) % 3]] == NETypeR && NodeTypes[SphereElems[ti][(c1 + 2) % 3]] == NETypeR) {
								/*
									*	All three nodes of the triangle are type R, meaning the
									*	gradient bundle has edges going to two different cages
									*	(It's impossible for a triangle to go to three different
									*	rings, as that would mean there's a bond point interior
									*	to the element, so it has to be two).
									*
									*	In this case we'll just look for the cage CP (last point
									*	of ring path segment) that's closest to the ring cps (midpoint of the
									*	end points of the ring surface paths) for all three nodes.
									*	This seems like an OK test since the cage point needs to be
									*	present (otherwise couldn't have two rings this close to each
									*	other) and is the common terminus of the correct paths, where
									*	the other two cage points are not common among the correct paths.
									*	(Even if the presumed cage point isn't in the system, the end
									*	points of the ring path segments both falling to where the cage 
									*	*would* be still provide an approximate location.)
									*/


									// Get midpoint of ring points;
								vec3 MidPt = { 0,0,0 };
								for (int ni : SphereElems[ti]) {
									MidPt += GradPaths[ni].XYZAt(-1);
								}
								MidPt /= SphereElems[ti].size();

								// Get which ring path segments terminate closer to the ring cp midpoint
								// and add that path.
								for (int ni = 0; ni < SphereElems[ti].size(); ++ni) {
									double MinDistSqr = DBL_MAX, TmpDistSqr;
									int MinPathNum;
									int SurfNum = NodeIntSurfNums[SphereElems[ti][ni]];
									for (int ii = 0; ii < 2; ++ii) {
										TmpDistSqr = DistSqr(IntSurfRingCagePaths[SurfNum][ii].XYZAt(-1), MidPt);
										if (TmpDistSqr < MinDistSqr) {
											MinDistSqr = TmpDistSqr;
											MinPathNum = ii;
										}
									}
									// 									GradPaths.push_back(ConcatenateResample({ GradPaths[SphereElems[ti][ni]], IntSurfRingCagePaths[SurfNum][MinPathNum] }, NumGPPoints, { -1, RNodeRingPathSegmentLength[SphereElems[ti][ni]] }, ResampleMethod));
									GradPaths.push_back(ConcatenateResample({ GradPaths[SphereElems[ti][ni]], IntSurfRingCagePaths[SurfNum][MinPathNum] }, NumGPPoints, {}, ResampleMethod));
									tri->at(ni) = GradPaths.size() - 1;
								}


								break;
							}

							for (int ei = 1; ei < 3; ++ei) {
								int c2 = (c1 + ei) % 3;

								string logStr = "elem " + to_string(ti) + ", node " + to_string(SphereElems[ti][c1]) + "-" + to_string(SphereElems[ti][c2]) + ": ";
								/*
									*	Check to see if this corner is R type
									*/
								if (ElemTypes[ti] == NETypeB && NodeTypes[SphereElems[ti][c2]] == NETypeR) {
									/*
										* Both corners c1 and c2 are R type.
										* They each correspond to a pair of ring path segments.
										* Because these are two corners of a triangle, two of the four
										* ring path segments must terminate at (or near) the same cage CP,
										* while the other two go off to different cage CPs.
										* We'll find the correct ring path segments by finding the pair with the nearest
										* end points.
										*/
									int minI, minJ, SurfNumI, SurfNumJ;
									SurfNumI = NodeIntSurfNums[SphereElems[ti][c1]];
									SurfNumJ = NodeIntSurfNums[SphereElems[ti][c2]];
									double MinDistSqr = DBL_MAX;
									for (int ii = 0; ii < 2; ++ii) {
										for (int jj = 0; jj < 2; ++jj) {
											double TmpDistSqr = DistSqr(IntSurfRingCagePaths[SurfNumI][ii].XYZAt(-1),
												IntSurfRingCagePaths[SurfNumJ][jj].XYZAt(-1));
											if (TmpDistSqr < MinDistSqr) {
												MinDistSqr = TmpDistSqr;
												minI = ii;
												minJ = jj;
											}
										}
									}
									/*
										*	Now the correct R paths for this triangle are those specified
										*	by minI and minJ
										*/


									tri->at(c1) = GradPaths.size();
									// 									GradPaths.push_back(ConcatenateResample({ GradPaths[SphereElems[ti][c1]], IntSurfRingCagePaths[SurfNumI][minI] }, NumGPPoints, { -1, RNodeRingPathSegmentLength[SphereElems[ti][c1]] }, ResampleMethod));
									GradPaths.push_back(ConcatenateResample({ GradPaths[SphereElems[ti][c1]], IntSurfRingCagePaths[SurfNumI][minI] }, NumGPPoints, {}, ResampleMethod));

									tri->at(c2) = GradPaths.size();
									// 									GradPaths.push_back(ConcatenateResample({ GradPaths[SphereElems[ti][c2]], IntSurfRingCagePaths[SurfNumJ][minJ] }, NumGPPoints, { -1, RNodeRingPathSegmentLength[SphereElems[ti][c2]] }, ResampleMethod));
									GradPaths.push_back(ConcatenateResample({ GradPaths[SphereElems[ti][c2]], IntSurfRingCagePaths[SurfNumJ][minJ] }, NumGPPoints, {}, ResampleMethod));

									break;
								}
								else if (NodeTypes[SphereElems[ti][c2]] == NETypeC) {
									/*
										*	Just need to check the two ring paths to see which terminates
										*	closer to the neighbor GP
										*/
									double MinDistSqr = DBL_MAX;
									int SurfNum = NodeIntSurfNums[SphereElems[ti][c1]];
									int minI;
									for (int ii = 0; ii < 2; ++ii) {
										double TmpDistSqr = DistSqr(IntSurfRingCagePaths[SurfNum][ii].XYZAt(-1), GradPaths[SphereElems[ti][c2]].XYZAt(-1));
										if (TmpDistSqr < MinDistSqr) {
											MinDistSqr = TmpDistSqr;
											minI = ii;
										}
									}

									tri->at(c1) = GradPaths.size();
									// 									GradPaths.push_back(ConcatenateResample({ GradPaths[SphereElems[ti][c1]], IntSurfRingCagePaths[SurfNum][minI] }, NumGPPoints, { -1, RNodeRingPathSegmentLength[SphereElems[ti][c1]] }, ResampleMethod));
									GradPaths.push_back(ConcatenateResample({ GradPaths[SphereElems[ti][c1]], IntSurfRingCagePaths[SurfNum][minI] }, NumGPPoints, {}, ResampleMethod));

									break;
								}
							}
						}
					}
				}
			}


			if (UserQuit) {
				TecUtilDataLoadEnd();
				StatusDrop(AddOnID);
				return;
			}

			std::map<Edge, std::pair<Edge, Edge> > GPIndEdgeToNodeIndEdgeMap;

			if (MaxEdgeGPs > 0) {
				/*
				 *	With all C and R type GPs found, can now create/update edge GPs by comparing
				 *	the maximum achieved iso-rho separation distance between GPs.
				 *	Triangular element edges that coincide with ring surface intersections
				 *	can be skipped because we're already assuming the behavior of rho along
				 *	the ring surface.
				 *	Only need to check edges that are of node type C-C, C-R, and C-B.
				 *	Doing this now because SphereElems == SphereElemGPInds at this point,
				 *	so edges made here can be used in subsequent post-subdivision iterations.
				 *
				 *	Do one pass in serial to get the unique list of edges of the right node types,
				 *	then a parallel pass to get GP distances for those edges,
				 *	then another serial pass to setup edge GPs,
				 *	then another parallel pass to make the new GPs.
				 *
				 *	An odd number of edge GPs will always be used so that no edge GPs need be
				 *	remade upon element subdivision, i.e. a GP will always lie at the edge
				 *	midpoint that will become a new node upon subdivision.
				 *
				 *	Here we make the assumption that the separation between a type C node GP between
				 *	either a type R or B node GP is well-approximated using only the partial R or B node
				 *	GP, that is, using only the partial R and B node GPs available at this point.
				 */


				Time1 = high_resolution_clock::now();
				TmpString = ProgressStr2.str() + "Preparing edge GPs";

				 // Get unique list of edges to check
				vector<std::pair<Edge, double> > EdgeDistList;
				std::set<Edge> TmpEdgeSet;
				std::map<Edge, std::pair<Edge, Edge> > GPIndEdgeToNodeIndEdgeMap;
				EdgeDistList.reserve(NumNodes * 2);

#ifdef SUPERDEBUG
				for (auto ti : ElemTodo) {
#else
				for (int ti = 0; ti < NumElems; ++ti) {
#endif
					auto tri = &SphereElems[ti];
					for (int ci = 0; ci < 3; ++ci) {
						int c2 = (ci + 1) % 3;
						Edge TmpEdge = MakeEdge(tri->at(ci), tri->at(c2), false);
						Edge TmpEdge2 = MakeEdge(SphereElemGPInds[ti][ci], SphereElemGPInds[ti][c2]);
						if (!EdgeGPMap.count(TmpEdge2))
						{
							TmpEdgeSet.insert(TmpEdge2);
						}
						GPIndEdgeToNodeIndEdgeMap[TmpEdge2] = std::make_pair(MakeEdge(SphereElemGPInds[ti][ci], SphereElemGPInds[ti][c2], false), MakeEdge(tri->at(ci), tri->at(c2), false));
					}
				}

				for (auto const & e : TmpEdgeSet) EdgeDistList.emplace_back(e, 0.0);

				if (!EdgeDistList.empty()) {
					// Get separation distances for edge GPs
					NumToDo = EdgeDistList.size();
					NumCompleted = 0;
#ifdef SUPERDEBUG
#pragma omp parallel for schedule(dynamic)
#endif
					for (int ei = 0; ei < NumToDo; ++ei) {
						int ThNum = omp_get_thread_num();
						if (ThNum == 0) {
							UserQuit = !StatusUpdate(NumCompleted, NumToDo, TmpString, AddOnID, Time1);
#pragma omp flush (UserQuit)
						}
#pragma omp flush (UserQuit)
						if (UserQuit) {
#pragma omp critical
							{
								if (omp_get_thread_num() == 0 && !TecUtilDialogMessageBox("Cancel run? All progress will be lost if you press \"Yes.\"", MessageBoxType_YesNo)) {
									StatusLaunch(TmpString, AddOnID, TRUE, TRUE, TRUE);
									StatusUpdate(NumCompleted, NumToDo, TmpString, AddOnID, Time1);
									UserQuit = false;
#pragma omp flush (UserQuit)
								}
#pragma omp flush (UserQuit)
							}
						}
						if (!UserQuit) {
#pragma omp atomic
							NumCompleted++;
							auto e = &EdgeDistList[ei];
							auto e2 = GPIndEdgeToNodeIndEdgeMap[e->first];
							auto en = e2.second;
							auto egp = e2.first;
							// If a R or B type node is present, use that as the basis for the distance check
							// because the distance check function here only uses the points of the path
							// that calls the method.
							// Otherwise, just use the path with the shortest number of points.
							// 
							vec3 TmpVec;
							int TmpInt;
							if (NodeTypes[en.first] >= NETypeB) {
								GradPath_c const * BPPtr = &IntBondPathSegments[en.first];
								BPPtr->GetMaxSeparationMidpointFromOtherGP(GradPaths[egp.second], TmpVec, TmpInt, TmpInt, TmpVec, e->second, MAX(1,BPPtr->GetCount() / 20));
							}
							else if (NodeTypes[en.second] >= NETypeB) {
								GradPath_c const * BPPtr = &IntBondPathSegments[en.second];
								BPPtr->GetMaxSeparationMidpointFromOtherGP(GradPaths[egp.first], TmpVec, TmpInt, TmpInt, TmpVec, e->second, MAX(1, BPPtr->GetCount() / 20));
							}
							else {
								GradPaths[e->first.first].GetMaxSeparationMidpointFromOtherGPRhoBased(GradPaths[e->first.second], TmpVec, TmpInt, TmpInt, TmpVec, e->second, MAX(1, GradPaths[e->first.first].GetCount() / 20));
							}
						}
					}


					if (UserQuit) {
						TecUtilDataLoadEnd();
						StatusDrop(AddOnID);
						return;
					}

					// Prepare edge GPs

					GradPath_c TmpGP;
// 					if (RadialSphereApprx)
						TmpGP.SetupGradPath(vec3(), StreamDir_Reverse, NumGPPoints * 2, GPType_Classic, GPTerminate_AtCP, nullptr, &CPs, &GPTermRadius, &CutoffVal, ThVolInfo[0], HessPtrs, GradPtrs, RhoPtr);
// 					else
// 						TmpGP.SetupGradPath(vec3(), StreamDir_Both, NumGPPoints * 2, GPType_Classic, GPTerminate_AtCP, nullptr, &CPs, &GPTermRadius, &CutoffVal, ThVolInfo[0], HessPtrs, GradPtrs, RhoPtr);
					TmpGP.SetTerminalCPTypeNum(CPTypeNum_Cage);

					int NumGPsOld = GradPaths.size();

					Time1 = high_resolution_clock::now();
					TmpString = ProgressStr2.str() + "Make Edge GPs";

					// 				NumEdges = EdgeDistList.size();
					NumEdges = 0;
					for (auto const & e : EdgeDistList) {
						auto e2 = GPIndEdgeToNodeIndEdgeMap[e.first];
						auto en = e2.second;
						double TmpCheckDistance = EdgeGPCheckDistance;
						if (NodeTypes[en.first] >= NETypeR || NodeTypes[en.second] >= NETypeR)
						{
							TmpCheckDistance *= RingBondEdgeGPCheckDistanceFactor;
						}
						if (e.second >= TmpCheckDistance) {
							NumEdges++;
						}
					}
					NumCompleted = 0;
					bool BadGP = false;
					int BadEdge1, BadEdge2, BadNode;
					int NumEdgeGPs = 0, TotNumEdgeGPs;
					int TotNumEdges = EdgeDistList.size();
					if (MinEdgeGPs) {
						NumToDo = EdgeDistList.size();
					}
					else {
						NumToDo = 0;
						for (auto const & e : EdgeDistList) {
							auto e2 = GPIndEdgeToNodeIndEdgeMap[e.first];
							auto en = e2.second;
							double TmpCheckDistance = EdgeGPCheckDistance;
							if (NodeTypes[en.first] >= NETypeR || NodeTypes[en.second] >= NETypeR)
							{
								TmpCheckDistance *= RingBondEdgeGPCheckDistanceFactor;
							}
							if (e.second >= TmpCheckDistance) {
								NumToDo++;
							}
						}
					}

					vector<vector<int> > EdgeGPNumListVec(TotNumEdges);

					// #pragma omp parallel for schedule(dynamic)
					for (int ei = 0; ei < TotNumEdges; ++ei) {
						// 					int ThNum = omp_get_thread_num();
						// 					NumCompleted = NumEdgeGPs;

						// 					TotNumEdgeGPs = NumEdges * int((double)NumEdgeGPs / double(ei + 1));
						// 					if (ThNum == 0) {
						UserQuit = !StatusUpdate(NumCompleted, NumToDo, TmpString, AddOnID, Time1);
						if (UserQuit){
							if (!TecUtilDialogMessageBox("Cancel run? All progress will be lost if you press \"Yes.\"", MessageBoxType_YesNo)) {
								StatusLaunch(TmpString, AddOnID, TRUE, TRUE, TRUE);
								StatusUpdate(NumCompleted, NumToDo, TmpString, AddOnID, Time1);
								UserQuit = false;
							}
						}
						if (UserQuit)
							break;

						auto e = EdgeDistList[ei];
						auto e2 = GPIndEdgeToNodeIndEdgeMap[e.first];
						auto en = e2.second;
						auto NodeEdge = MakeEdge(en);
						auto egp = e2.first;
						
						double TmpCheckDistance = EdgeGPCheckDistance;
						if (NodeTypes[en.first] >= NETypeR || NodeTypes[en.second] >= NETypeR)
						{
							TmpCheckDistance *= RingBondEdgeGPCheckDistanceFactor;
						}

// 						if ((NodeTypes[en.first] == NETypeB && NodeTypes[en.second] != NETypeB)
// 							|| (NodeTypes[en.first] != NETypeB && NodeTypes[en.second] == NETypeB)
// 							|| (NodeTypes[en.first] == NETypeR && NodeTypes[en.second] == NETypeC)
// 							|| (NodeTypes[en.first] == NETypeC && NodeTypes[en.second] == NETypeR)
// 							)
// 						{
// 							TmpCheckDistance *= RingBondEdgeGPCheckDistanceFactor;
// 						}
// 						else if (
// 							NodeTypes[en.first] == NETypeR
// 							&& NodeTypes[en.second] == NETypeR
// 							&& NodeIntSurfNums[en.first] == NodeIntSurfNums[en.second])
// 						{
// 							TmpCheckDistance *= RingBondEdgeGPCheckDistanceFactor * 0.5;
// 						}

						if (!UserQuit && (MinEdgeGPs || e.second >= TmpCheckDistance)) {
							NumCompleted++;
							std::list<std::pair<vec3, GradPath_c> > EdgeGPs;
							vector<int> EdgeGPNums;
							if (NodeTypes[en.first] == NETypeB) {
								EdgeGPs.emplace_back(SphereNodes[en.first], IntBondPathSegments[en.first]);
							}
							else {
								EdgeGPs.emplace_back(SphereNodes[en.first], GradPaths[egp.first]);
							}
							if (!EdgeGPs.back().second.IsMade()) {
#pragma omp critical(BadGP1)
								{
									BadGP = true;
									BadNode = e.first.first;
									BadEdge1 = e.first.first;
									BadEdge2 = e.first.second;
									TecUtilDialogErrMsg(string("Bad GP in edge GP generation. edge (" + to_string(BadEdge1) + "," + to_string(BadEdge2) + "), node " + to_string(BadNode)).c_str());
								}
							}
							if (NodeTypes[en.second] == NETypeB) {
								EdgeGPs.emplace_back(SphereNodes[en.second], IntBondPathSegments[en.second]);
							}
							else {
								EdgeGPs.emplace_back(SphereNodes[en.second], GradPaths[egp.second]);
							}

							if (!EdgeGPs.back().second.IsMade()) {
#pragma omp critical(BadGP2)
								{
									BadGP = true;
									BadNode = e.first.second;
									BadEdge1 = e.first.first;
									BadEdge2 = e.first.second;
									TecUtilDialogErrMsg(string("Bad GP in edge GP generation. edge (" + to_string(BadEdge1) + "," + to_string(BadEdge2) + "), node " + to_string(BadNode)).c_str());
								}
							}

							for (auto & e : EdgeGPs){
								e.second.SetVolInfo(VolInfo);
								e.second.SetRhoPtr(RhoPtr);
							}

							FESurface_c const * SurfPtr = nullptr;
							GradPath_c const * GPPtr = nullptr;
							if ((// Edge corresponds to a ring surface intersection, so get pointers to the surface
								// and the correct ring path to concatenate the surface path edge paths with.
								NodeTypes[en.first] == NETypeR
								&& NodeTypes[en.second] == NETypeR
								&& EdgeIntSurfNums.count(NodeEdge))
								// 								&& NodeIntSurfNums[en.first] == NodeIntSurfNums[en.second])
								|| (// Edge corresponds to a ring surface intersection, and the second edge node is a bond path type.
									NodeTypes[en.first] == NETypeR
									&& NodeTypes[en.second] == NETypeB)
								)
							{

								// 								int SurfNum = NodeIntSurfNums[en.first];
								int SurfNum = EdgeIntSurfNums[NodeEdge];
								SurfPtr = &IntSurfs[SurfNum];

								// Get the ring path to finish off the new edge GPs.
								// Ring surface GPs have already been made, so just check
								// distance of terminal ring paths points to terminal
								// points of the ring surface GPs.
								if (DistSqr(GradPaths[egp.first][-1], IntSurfRingCagePaths[SurfNum][0][-1]) < DistSqr(GradPaths[egp.first][-1], IntSurfRingCagePaths[SurfNum][1][-1])) {
									GPPtr = &IntSurfRingCagePaths[SurfNum][0];
								}
								else {
									GPPtr = &IntSurfRingCagePaths[SurfNum][1];
								}
							}
							else if (
								// Edge corresponds to a ring surface intersection, and the first edge node is a bond path type.
								NodeTypes[en.second] == NETypeR
								&& NodeTypes[en.first] == NETypeB)
							{
								int SurfNum = NodeIntSurfNums[en.second];
								SurfPtr = &IntSurfs[SurfNum];

								// Get the ring path to finish off the new edge GPs.
								// Ring surface GPs have already been made, so just check
								// distance of terminal ring paths points to terminal
								// points of the ring surface GPs.
								if (DistSqr(GradPaths[egp.second][-1], IntSurfRingCagePaths[SurfNum][0][-1]) < DistSqr(GradPaths[egp.second][-1], IntSurfRingCagePaths[SurfNum][1][-1])) {
									GPPtr = &IntSurfRingCagePaths[SurfNum][0];
								}
								else {
									GPPtr = &IntSurfRingCagePaths[SurfNum][1];
								}
							}

							auto Pt1 = EdgeGPs.cbegin();
							auto Pt2 = Pt1;
							Pt2++;
							double TmpDist;
							vec3 TmpVec1, TmpVec2;
							int TmpInd1, TmpInd2;
							auto GP = TmpGP;
							auto GPInner = GP;
							GP.SetDir(StreamDir_Reverse);
							GPInner.SetDir(StreamDir_Forward);
							GPInner.SetTermPoint(CPPos);
							GPInner.SetGPTermType(GPTerminate_AtPoint);
							GPInner.SetTerminalCPTypeNum(CPTypeNum_Nuclear);
							if (SurfPtr != nullptr) {
								GP.SetSurfPtr(SurfPtr);
								GP.SetTerminalCPTypeNum(CPTypeNum_Ring);
							}
							else {
								GP.SetTerminalCPTypeNum(CPTypeNum_Cage);
							}
							
							bool DidConverge = true;
							int Iter = 0;
							while (!UserQuit && Pt2 != EdgeGPs.cend() && Iter <= MaxEdgeGPs && DidConverge)
							{
								UserQuit = !StatusUpdate(NumCompleted, NumToDo, TmpString + ": " + to_string(Iter + 1), AddOnID, Time1);
								if (UserQuit) {
									if (!TecUtilDialogMessageBox("Cancel run? All progress will be lost if you press \"Yes.\"", MessageBoxType_YesNo)) {
										StatusLaunch(TmpString, AddOnID, TRUE, TRUE, TRUE);
										StatusUpdate(NumCompleted, NumToDo, TmpString + ": " + to_string(Iter + 1), AddOnID, Time1);
										UserQuit = false;
									}
								}
								if (UserQuit)
									break;
								// Always check of GP2 to GP1
								// If one GP is a bond path segment, need to use that as the basis for the distance check
								// 
								if (approx_equal(Pt1->first, Pt2->first, "absdiff", 1e-10)) {
									DidConverge = !EdgeGPs.size() > 2;
								}
								Iter++;
								if (DidConverge) {
									GradPath_c const *GP1 = &Pt1->second,
										*GP2 = &Pt2->second;
									bool RBPathPresent = false;
									if (NodeTypes[en.second] >= NETypeR && sum(Pt2->first == SphereNodes[en.second]) == 3) {
										GP1 = &Pt2->second;
										GP2 = &Pt1->second;
										RBPathPresent = true;
									}
									if (!RBPathPresent && NodeTypes[en.first] >= NETypeR && sum(Pt1->first == SphereNodes[en.first]) == 3)
										RBPathPresent = true;
									bool DistCheckSuccessful;
									if (RBPathPresent)
										DistCheckSuccessful = GP1->GetMaxSeparationMidpointFromOtherGP(*GP2, TmpVec1, TmpInd1, TmpInd2, TmpVec2, TmpDist);
									else
										DistCheckSuccessful = GP1->GetMaxSeparationMidpointFromOtherGPRhoBased(*GP2, TmpVec1, TmpInd1, TmpInd2, TmpVec2, TmpDist);
									if (DistCheckSuccessful) {
										if (TmpDist > TmpCheckDistance || (Iter == 1 && MinEdgeGPs)) {
											vec3 Midpt;
											vec3 EdgeMidpt = (Pt1->first + Pt2->first) * 0.5;
											vector<GradPath_c> GPList;
											vector<int> GPNumPointsList;
											if (NodeTypes[en.first] == NodeTypes[en.second]
												&& (NodeTypes[en.first] == NETypeC || EdgeIntSurfNums.count(en))) {
												// no strong deviation along cage-cage edges or ring-ring (if they're both on
												// the same ring, so simply seed an edge 
												// GP from the midpoint of Pt1 and Pt2 on the edge
												GP.Clear();
												GP.SetStartPoint(EdgeMidpt);
												GP.Seed(false);
												GPList = { GP };
												GPNumPointsList.push_back(-1);

												if (!RadialSphereApprx) {
													GPInner.Clear();
													GPInner.SetStartPoint(EdgeMidpt);
													GPInner.Seed(false);
													GPList.push_back(GPInner);
													GPNumPointsList.push_back(-1);
												}
											}
											else {
												double EdgeDist = Distance(Pt1->first, Pt2->first);

												if (GP1->GetSeparationMidpointAtDistFromOtherGP(*GP2, Midpt, TmpInd1, TmpInd2, TmpVec2, MAX(EdgeDist * 2.0, TmpDist * 0.1), 0.0)) {
													// all other types of edge have one path that makes one or more 90-degree turns
													// while the other path makes fewer, so we'll use that region of path deviation
													// to define the new GPs
													GP.Clear();
													GP.SetStartPoint(Midpt);
													GP.Seed(false);
													GPList = { GP };
													GPNumPointsList.push_back(-1);
													if (RadialSphereApprx) {
// 														GPInner.Clear();
// 														GPInner.SetStartPoint(Midpt);
// 														GPInner.Seed(false);
// 														GPInner.Trim(CPPos, Radius);
// 														GPList.push_back(GPInner);
// 														GPNumPointsList.push_back(-1);
// 														EdgeMidpt = GPInner[-1];

														double TermDist = Distance(EdgeMidpt, Pt1->first) * 3;
														GradPath_c BackGP({ GP1,GP2 }, Midpt, StreamDir_Forward, { std::make_pair(0, MIN(TmpInd1 + 5, GP1->GetCount() - 1)), std::make_pair(0, MIN(TmpInd2 + 5, GP2->GetCount() - 1)) }, &EdgeMidpt, &TermDist, GPTerminate_AtPointRadius);
// 														BackGP = BackGP.SubGP(0, -2);
														vec3 endPt = CPPos + normalise(BackGP[-1] - CPPos) * Radius;
														endPt = ProjectPointToLine(endPt, SphereNodes[en.first], SphereNodes[en.second]);
														EdgeMidpt = endPt;
														BackGP.PointAppend(endPt, GP1->RhoAt(0));
														GPList.push_back(BackGP);
														GPNumPointsList.push_back(-1);
													}
													else {
// 														GPInner.Clear();
// 														GPInner.SetStartPoint(Midpt);
// 														GPInner.Seed(false);
// 														GPList.push_back(GPInner);
// 														GPNumPointsList.push_back(-1);
// 														GPInner.Trim(CPPos, Radius);
// 														EdgeMidpt = GPInner[-1];

														GradPath_c BackGP({ GP1,GP2 }, Midpt, StreamDir_Forward, { std::make_pair(0, MIN(TmpInd1 + 5, GP1->GetCount() - 1)), std::make_pair(0, MIN(TmpInd2 + 5, GP2->GetCount() - 1)) }, &CPPos, &GPTermRadius);
														GPList.push_back(BackGP);
														GPNumPointsList.push_back(-1);
														BackGP.Trim(CPPos, Radius);
														EdgeMidpt = BackGP[-1];
													}
												}
											}

											if (SurfPtr != nullptr) {
												GPList.push_back(*GPPtr);
												GPNumPointsList.push_back(-1);
											}

											GP = ConcatenateResample(GPList, NumGPPoints, {}, ResampleMethod);
											GP.MakeRhoValuesMonotonic(&VolInfo, &RhoPtr);
											EdgeGPs.emplace(Pt2, EdgeMidpt, GP);
											Pt2--;
											continue;
										}
									}
									Pt1++;
									Pt2++;
								}
								else if (NodeTypes[en.second] != NETypeB && NodeTypes[en.first] != NETypeB) {
									// Deviating edge found, which means that not all ring surfaces were found prior to
									// running GBA. Ideally, there will be a ring CP in proximity to the midpoint
									// between the deviating Pt1 and Pt2 gradient paths.
									// Seed a GP from the current midpoint.
									// If there is a ring CP present, then connect the GP to the ring CP right before
									// its closest pass and then cap it with the ring paths from the ring CP.
									// If no ring CP is present then the new GP will deviate to one side or the other,
									// and we'll find the deviation midpoint between the new GP and the one going
									// the other direction, and do a binary GP search to get one going in the
									// opposite direction, and use those to cap the midpoint GP at its closest point
									// to the deviation point.
								}
							}

							if (UserQuit)
								break;

							EdgeGPs.pop_back();
							EdgeGPs.pop_front();
							for (auto gp : EdgeGPs) {
								if (gp.second.RhoAt(0) < gp.second.RhoAt(-1))
									gp.second.Reverse();
#pragma omp critical(UpdateGradPaths)
								{
									EdgeGPNums.push_back(GradPaths.size());
									GradPaths.push_back(gp.second);
								}
							}
#pragma omp critical(UpdateEdgeGPNumList)
							EdgeGPNumListVec[ei] = EdgeGPNums;
							// 						NumEdgeGPs += EdgeGPNums.size();
						}
					}

					if (UserQuit) {
						TecUtilDataLoadEnd();
						StatusDrop(AddOnID);
						return;
					}

					for (int ei = 0; ei < EdgeGPNumListVec.size(); ++ei) {
						if (!EdgeGPNumListVec[ei].empty()) {
							EdgeGPMap[EdgeDistList[ei].first] = EdgeGPNumListVec[ei];
						}
					}


					// 				for (auto const & i : EdgeGPMap) {
					// 					auto e2 = NodeToGPIndMap[i.first];
					// 					GradPaths[e2.first].SaveAsOrderedZone("GPs for edge GP creation");
					// 					GradPaths[e2.second].SaveAsOrderedZone("GPs for edge GP creation");
					// 					for (auto j : i.second) GradPaths[j].SaveAsOrderedZone("edge (" + to_string(i.first.first) + "," + to_string(i.first.second) + ") gp " + to_string(j));
					// 				}

				}
			}

			// For all B type nodes, get the average of the index of the closest
			// points to the bond point for all neighboring node GPs.
			// (This now doesn't do anything, as the GPs are "aligned" according to rho value.
			// 
			std::map<int, int> BNodeBondPathSegmentLength;
#ifdef SUPERDEBUG
			for (int nii = 0; nii < NodesTodoVec.size(); ++nii) {
				int ni = NodesTodoVec[nii];
#else
			for (int ni = 0; ni < NumNodes; ++ni) {
#endif
				if (NodeTypes[ni] == NETypeB) {
					// Get bond point
// 					vec3 BondCP = IntBondPathSegments[ni].XYZAt(0);
// 
// 					// Get index of closest point to bond point for neighboring node GPs
// 
//  					vector<double> ClosestPointNums;
//  					for (int nj : NodeConnectivity[ni]) {
//  						// If edge GPs exist between ni and nj, then use the closest GP to ni
//  						// on the edge for this check.
//  						GradPath_c const * GPPtr = nullptr;
//  						Edge TmpEdge = MakeEdge(ni, nj);
//  						if (EdgeGPMap.count(TmpEdge)){
//  							if (ni < nj)
//  								GPPtr = &GradPaths[EdgeGPMap[TmpEdge].front()];
//  							else
//  								GPPtr = &GradPaths[EdgeGPMap[TmpEdge].back()];
//  						}
//  						if (GPPtr == nullptr){
//  							GPPtr = &GradPaths[nj];
//  						}
//  						int PointNum;
//  						if (GPPtr != nullptr && GPPtr->IsMade()) {
//  							GPPtr->ClosestPoint(BondCP, PointNum);
//  							ClosestPointNums.push_back(PointNum);
//  						}
//  					}
//  
//  					int AvgClosestPoint = mean(vec(ClosestPointNums));
//  					BNodeBondPathSegmentLength[ni] = AvgClosestPoint;
// 					
					// Different method: Want the same point spacing on the bond path as 
					// its neighboring paths, so get the average point spacing of neighboring
					// paths then see which point on the bond path should be used in order
					// to achieve the same point spacing
					// 
					vector<double> PointSpacings;
					for (int nj : NodeConnectivity[ni]){
						// If edge GPs exist between ni and nj, then use the closest GP to ni
						// on the edge for this check.
						GradPath_c const * GPPtr = nullptr;
						Edge TmpEdge = MakeEdge(ni, nj);
						if (EdgeGPMap.count(TmpEdge)) {
							if (ni < nj)
								GPPtr = &GradPaths[EdgeGPMap[TmpEdge].front()];
							else
								GPPtr = &GradPaths[EdgeGPMap[TmpEdge].back()];
						}
						if (GPPtr == nullptr) {
							GPPtr = &GradPaths[nj];
						}
						if (GPPtr != nullptr && GPPtr->IsMade()) {
							PointSpacings.push_back(GPPtr->GetLength() / (double)GPPtr->GetCount());
						}
					}
					double AvgPointSpacing = mean(vec(PointSpacings));
					int BPNumPoints = IntBondPathSegments[ni].GetLength() / AvgPointSpacing;
					BNodeBondPathSegmentLength[ni] = BPNumPoints;
				}
			}

			Time1 = high_resolution_clock::now();
			TmpString = ProgressStr2.str() + "Make bond path GPs";

#ifdef SUPERDEBUG
			vector<int> BTypeElemTodo;
			for (int ti : ElemTodo){
				if (ElemTypes[ti] >= NETypeB){
					BTypeElemTodo.push_back(ti);
				}
			}
			NumToDo = BTypeElemTodo.size();
			NumCompleted = 0;
			for (int ti : BTypeElemTodo) {
				UserQuit = !StatusUpdate(NumCompleted++, NumToDo, TmpString, AddOnID, Time1);
				if (UserQuit){
					if (TecUtilDialogMessageBox("Cancel run? All progress will be lost if you press \"Yes.\"", MessageBoxType_YesNo)) {
						break;
					}
					else {
						StatusLaunch(TmpString, AddOnID, TRUE, TRUE, TRUE);
						StatusUpdate(NumCompleted, NumToDo, TmpString, AddOnID, Time1);
						UserQuit = false;
					}
				}
#else
			NumToDo = 0;
			for (int ti = 0; ti < NumElems; ++ti)
				NumToDo += int(ElemTypes[ti] >= NETypeB);
			NumCompleted = 0;

			for (int ti = 0; ti < NumElems && !UserQuit; ++ti) {
				UserQuit = !StatusUpdate(NumCompleted, NumToDo, TmpString, AddOnID, Time1);
				if (UserQuit){
					if (TecUtilDialogMessageBox("Cancel run? All progress will be lost if you press \"Yes.\"", MessageBoxType_YesNo)) {
						break;
					}
					else {
						StatusLaunch(TmpString, AddOnID, TRUE, TRUE, TRUE);
						StatusUpdate(NumCompleted, NumToDo, TmpString, AddOnID, Time1);
						UserQuit = false;
					}
				}
				if (ElemTypes[ti] >= NETypeB)
					NumCompleted++;
	#endif
				/*
				 * Second pass over triangle corners to setup the B type
				 * nodes now that we have the C and R type paths.
				 */
				auto * tri = &SphereElemGPInds[ti];
				for (int c1 = 0; c1 < 3; ++c1) {
					if (NodeTypes[SphereElems[ti][c1]] == NETypeB) {
						/*
						 *	Type B or RB found, so need to  determine where to
						 *	seed gradient paths in the interatomic surface for
						 *	the bond path. We'll get the GPs for the
						 *	other two corners of the triangle and find the closest point
						 *	to the bond point and use that to determine the seed direction.
						 */
						GradPath_c const * BPPtr = &IntBondPathSegments[SphereElems[ti][c1]];
						vec3 BCP = BPPtr->XYZAt(-1);
						double BCPRho = BPPtr->RhoAt(-1);
						int BondCPNum;
						for (int bi = 0; bi < CPs.NumBonds(); ++bi) {
							if (sum(CPs.GetXYZ(CPTypeNum_Bond, bi) == BCP) == 3) {
								BondCPNum = CPs.GetTotOffsetFromTypeNumOffset(CPTypeNum_Bond, bi);
								break;
							}
						}


						int NewGPNums[2];

						GradPath_c NewGPs[2];

						for (int ei = 1; ei <= 2; ++ei) {
							int c2 = (c1 + ei) % 3;


							string logStr = "elem " + to_string(ti) + ", node " + to_string(SphereElems[ti][c1]) + "-" + to_string(SphereElems[ti][c2]) + ": ";

							GradPath_c const * GPPtr;
// 							Edge TmpEdge = MakeEdge(SphereElems[ti][c1], SphereElems[ti][c2]);
							Edge TmpEdge = MakeEdge(SphereElemGPInds[ti][c1], SphereElemGPInds[ti][c2]);
							if (EdgeGPMap.count(TmpEdge)){
								if (DistSqr(BCP, GradPaths[EdgeGPMap[TmpEdge].front()].ClosestPoint(BCP)) < DistSqr(BCP, GradPaths[EdgeGPMap[TmpEdge].back()].ClosestPoint(BCP)))
									GPPtr = &GradPaths[EdgeGPMap[TmpEdge].front()];
								else
									GPPtr = &GradPaths[EdgeGPMap[TmpEdge].back()];
							}
							else{
								GPPtr = &GradPaths[SphereElems[ti][c2]];
							}

//  							vec3 EigenValues;
//  							mat33 EigenVectors;
//  							CalcEigenSystemForPoint(BPPtr->XYZAt(0), EigenValues, EigenVectors, MR);

							int ClosestPointNum;
// 							vec3 ClosestPt = GradPaths[tri->at(c2)].ClosestPoint(BPPtr->XYZAt(0), ClosestPointNum);
							vec3 ClosestPt = GPPtr->ClosestPoint(BCP, ClosestPointNum);
							vec3 v1 = BPPtr->XYZAt(-5) - BPPtr->XYZAt(-1);
							//  								vec3 v1 = EigenVectors.col(2);
// 							vec3 v1 = CPs.GetPrincDir(BondCPNum);
							vec3 v2 = ClosestPt - BCP;
							vec3 v3 = cross(v1, v2);

							v2 = normalise(cross(v1, v3));

							vec3 SeedPt = BCP + v2 * 0.01;

							if (DistSqr(SeedPt, ClosestPt) > DistSqr(BCP, ClosestPt))
								SeedPt = BCP + v2 * (-0.01);

							FESurface_c const * SurfPtr = nullptr;
							GradPath_c const * RPPtr = nullptr;
							if (NodeTypes[SphereElems[ti][c2]] == NETypeR) {
								ElemTypes[ti] = NETypeRB;
								int SurfNum = NodeIntSurfNums[SphereElems[ti][c2]];
								SurfPtr = &IntSurfs[SurfNum];

								/*
								 * Need to get the correct ring path segment to finish off the GP
								 */
								double MinDistSqr = DBL_MAX, TmpDistSqr;
								int MinI = -1;
								for (int ii = 0; ii < 2; ++ii) {
									TmpDistSqr = DistSqr(IntSurfRingCagePaths[SurfNum][ii].XYZAt(-1), GradPaths[tri->at(c2)].XYZAt(-1));
									if (TmpDistSqr < MinDistSqr) {
										MinDistSqr = TmpDistSqr;
										MinI = ii;
									}
								}
								if (MinI >= 0) {
									RPPtr = &IntSurfRingCagePaths[SurfNum][MinI];
								}
								else{
									TecUtilDialogErrMsg("Failed to find ring cage path to complete atom-bond-ring-cage path");
								}
							}

							GradPath_c GP(SeedPt, StreamDir_Reverse, NumGPPoints * 2, GPType_Classic, GPTerminate_AtCP, nullptr, &CPs, &GPTermRadius, &CutoffVal, VolInfo, HessPtrs, GradPtrs, RhoPtr);
							GP.SetSurfPtr(SurfPtr, 200);
							if (NodeTypes[SphereElems[ti][c2]] == NETypeR)
								GP.SetTerminalCPTypeNum(CPTypeNum_Ring);
							else
								GP.SetTerminalCPTypeNum(CPTypeNum_Cage);

							GP.SetStartEndCPNum(BondCPNum, 0);
							GP.Seed();
							if (SurfPtr != nullptr)
								GP.ProjectPathToSurface();
							GP.MakeRhoValuesMonotonic(&VolInfo, &RhoPtr);

							NewGPs[ei - 1] = GP;
							NewGPs[ei - 1].PointPrepend(BCP, BCPRho);

							vector<GradPath_c> GPList = { GP };
							vector<int> GPNumPointList = { -1 };

							string tmpStr;
							if (NodeTypes[SphereElems[ti][c2]] == NETypeR) {
								int TmpInt;
								NewGPs[ei - 1] = ConcatenateResample({ NewGPs[ei - 1], *RPPtr }, NumGPPoints, {}, ResampleMethod);
								NewGPs[ei - 1].MakeRhoValuesMonotonic(&VolInfo, &RhoPtr);
								GPList.push_back(*RPPtr);
								GPNumPointList.push_back(-1);
							}
							GPList.push_back(*BPPtr);
							GPNumPointList.push_back(BNodeBondPathSegmentLength[SphereElems[ti][c1]]);

// 							GP = ConcatenateResample(GPList, NumGPPoints, GPNumPointList, ResampleMethod);
							GP = ConcatenateResample(GPList, NumGPPoints, {}, ResampleMethod);

							if (GP.RhoAt(0) < GP.RhoAt(-1))
								GP.Reverse();

							int GPInd = GradPaths.size();
							GradPaths.push_back(GP);

							NewGPNums[ei - 1] = GPInd;
						}

						tri->insert(tri->begin() + c1, NewGPNums[1]);
						tri->at(c1 + 1) = NewGPNums[0];

						// Now reorder tri so that the new GPs are the first and last in the tri
						auto tmptri = *tri;
						for (int i = 0; i < tmptri.size(); ++i){
							tri->at(i) = tmptri[(c1 + i + 1) % tmptri.size()];
						}

						


						// Now add edge GPs if present between opposite corners of the element
						int c2 = (c1 + 1) % 3,
							c3 = (c1 + 2) % 3;
// 						Edge TmpEdge = MakeEdge(SphereElems[ti][c2], SphereElems[ti][c3]);
// 						Edge TmpEdge = MakeEdge(SphereElemGPInds[ti][c2], SphereElemGPInds[ti][c3]);
						int BondPathPtNum = BNodeBondPathSegmentLength[SphereElems[ti][c1]];
						if (MaxEdgeGPs > 0){
							bool BadGP = false;
							int BadEdge1, BadEdge2, BadNode;

							GradPath_c GP(vec3(), StreamDir_Reverse, NumGPPoints * 2, GPType_Classic, GPTerminate_AtCP, nullptr, &CPs, &GPTermRadius, &CutoffVal, VolInfo, HessPtrs, GradPtrs, RhoPtr);
							GP.SetTerminalCPTypeNum(CPTypeNum_Cage);
							vector<int> GPNumPointList = { -1, BondPathPtNum };
							std::list<std::pair<vec3, GradPath_c> > EdgeGPs;
							// Move bounding seed points 1% closer to eachother
							vec3 SeedPts[2];
							for (int i = 0; i < 2; ++i) {
								SeedPts[i] = NewGPs[i][1];
// 								NewGPs[i].SaveAsOrderedZone("initial GP " + to_string(i + 1));
							}

// 							SaveVec3VecAsScatterZone({ SeedPts[0],SeedPts[1] }, "initial seed pts");

							vec3 TmpVec = SeedPts[1] - SeedPts[0];
							SeedPts[0] += TmpVec * 0.005;
							SeedPts[1] -= TmpVec * 0.005;

// 							SaveVec3VecAsScatterZone({ SeedPts[0],SeedPts[1] }, "shifted initial seed pts");
 							EdgeGPs.emplace_back(SeedPts[1], NewGPs[1]);
// 							if (!EdgeGPs.back().second.IsMade()) {
// #pragma omp critical(BadGP1)
// 								{
// 									BadGP = true;
// 									BadNode = TmpEdge.first;
// 									BadEdge1 = TmpEdge.first;
// 									BadEdge2 = TmpEdge.second;
// 									TecUtilDialogErrMsg(string("Bad GP in bond path edge GP generation. edge (" + to_string(BadEdge1) + "," + to_string(BadEdge2) + "), node " + to_string(BadNode)).c_str());
// 								}
// 							}
 							EdgeGPs.emplace_back(SeedPts[0], NewGPs[0]);
// 							if (!EdgeGPs.back().second.IsMade()) {
// #pragma omp critical(BadGP2)
// 								{
// 									BadGP = true;
// 									BadNode = TmpEdge.first;
// 									BadEdge1 = TmpEdge.first;
// 									BadEdge2 = TmpEdge.second;
// 									TecUtilDialogErrMsg(string("Bad GP in bond path edge GP generation. edge (" + to_string(BadEdge1) + "," + to_string(BadEdge2) + "), node " + to_string(BadNode)).c_str());
// 								}
// 							}
							auto Pt1 = EdgeGPs.begin();
							auto Pt2 = Pt1;
							Pt2++;
							double TmpDist;
							vec3 TmpVec1, TmpVec2;
							vec3 Midpt;
							int TmpInd1, TmpInd2;
							int Iter = 0;
							double TmpCheckDistance = EdgeGPCheckDistance * RingBondEdgeGPCheckDistanceFactor;
							double TmpTermDist = 0.05;
							while (!UserQuit && Pt2 != EdgeGPs.end() && Iter <= MaxEdgeGPs && !approx_equal(Pt1->first, Pt2->first, "absdiff", 1e-10)) {
								// 								if (Pt1->second.GetMaxSeparationMidpointFromOtherGPRhoBased(Pt2->second, TmpVec1, TmpInd1, TmpInd2, TmpVec2, TmpDist)) {
// 								if (Pt1->second.GetMaxSeparationMidpointFromOtherGP(Pt2->second, TmpVec1, TmpInd1, TmpInd2, TmpVec2, TmpDist)) {
								UserQuit = !StatusUpdate(NumCompleted, NumToDo, TmpString + ": " + to_string(Iter + 1), AddOnID, Time1);
								if (UserQuit){
									if (TecUtilDialogMessageBox("Cancel run? All progress will be lost if you press \"Yes.\"", MessageBoxType_YesNo)) {
										break;
									}
									else {
										StatusLaunch(TmpString, AddOnID, TRUE, TRUE, TRUE);
										StatusUpdate(NumCompleted, NumToDo, TmpString + ": " + to_string(Iter + 1), AddOnID, Time1);
										UserQuit = false;
									}
								}
								Iter++;
								if (Pt1->second.GetSeparationMidpointAtDistFromOtherGP(Pt2->second, Midpt, TmpInd1, TmpInd2, TmpVec2, TmpCheckDistance, 0.01) || (Iter == 1 && MinEdgeGPs && Pt1->second.GetMaxSeparationMidpointFromOtherGP(Pt2->second, Midpt, TmpInd1, TmpInd2, TmpVec2, TmpDist) && Pt1->second.GetSeparationMidpointAtDistFromOtherGP(Pt2->second, Midpt, TmpInd1, TmpInd2, TmpVec2, TmpDist * 0.3, 0.01))){
// 									if (TmpDist > TmpCheckDistance){
// 										Midpt = (Pt1->first + Pt2->first) * 0.5;
// 										SaveVec3VecAsScatterZone({ Midpt }, "midpoint " + to_string(Iter));
										Pt1->second.SetTermPointRadius(0.05);
// 										TecUtilDialogMessageBox(string("bond edge gp before surf path iter " + to_string(Iter)).c_str(), MessageBoxType_Information);
										GradPath_c TmpGP({ &Pt1->second,&Pt2->second }, Midpt, StreamDir_Forward, { std::make_pair(0, MIN(TmpInd1 + 10, Pt1->second.GetCount() - 1)), std::make_pair(0, MIN(TmpInd2 + 10, Pt2->second.GetCount() - 1)) }, &BCP, &TmpTermDist);
										if (!TmpGP.IsMade()) {
											TecUtilDialogErrMsg("Failed to fill gradient bundle surface GPs along inter-atomic surface.");
											break;
										}
// 										TmpGP.SaveAsOrderedZone("midpoint gp backward " + to_string(Iter));
										GP.Clear();
										GP.SetStartPoint(Midpt);
// 										TecUtilDialogMessageBox(string("bond edge gp before free path iter " + to_string(Iter)).c_str(), MessageBoxType_Information);
										GP.Seed(false);
										GP.RemoveKinks();
// 										GP.SaveAsOrderedZone("midpoint gp forward " + to_string(Iter));
										GP = ConcatenateResample({ GP, TmpGP }, NumGPPoints, {}, ResampleMethod);
										GP.MakeRhoValuesMonotonic(&VolInfo, &RhoPtr);
										if (DistSqr(BCP, GP[0]) > DistSqr(BCP, GP[-1]))
											GP.Reverse();
// 										GP.SaveAsOrderedZone("midpoint gp both " + to_string(Iter));
										EdgeGPs.emplace(Pt2, Midpt, GP);
										// 										EdgeGPs.emplace(Pt2, Midpt, ConcatenateResample({ GP, *BPPtr }, NumGPPoints, GPNumPointList, ResampleMethod));
// 										TecUtilDialogMessageBox(string("bond edge gp iter " + to_string(Iter)).c_str(), MessageBoxType_Information);
										Pt2--;
										continue;
// 									}
								}
// 								else{
// 									TecUtilDialogErrMsg("Failed to get bong path GB edge GPs");
// 								}
								Pt1++;
								Pt2++;
							}
							if (UserQuit)
								break;
							EdgeGPs.pop_back();
							EdgeGPs.pop_front();
							double Dist1 = 1, Dist2 = 2;
							if (EdgeGPs.size() > 1) {
								NewGPs[1].GetMaxSeparationMidpointFromOtherGP(EdgeGPs.front().second, TmpVec1, TmpInd1, TmpInd2, TmpVec2, Dist1);
								NewGPs[0].GetMaxSeparationMidpointFromOtherGP(EdgeGPs.front().second, TmpVec1, TmpInd1, TmpInd2, TmpVec2, Dist2);
// 								GradPaths[NewGPNums[1]].GetMaxSeparationMidpointFromOtherGP(EdgeGPs.front().second, TmpVec1, TmpInd1, TmpInd2, TmpVec2, Dist1);
// 								GradPaths[NewGPNums[0]].GetMaxSeparationMidpointFromOtherGP(EdgeGPs.front().second, TmpVec1, TmpInd1, TmpInd2, TmpVec2, Dist2);
							}
							if (Dist1 < Dist2){
								for (auto gp = EdgeGPs.begin(); gp != EdgeGPs.end(); ++gp){
									tri->push_back(GradPaths.size());
									GradPaths.push_back(ConcatenateResample({ gp->second, *BPPtr }, NumGPPoints, {}, ResampleMethod));
// 									GradPaths.push_back(ConcatenateResample({ gp->second, *BPPtr }, NumGPPoints, GPNumPointList, ResampleMethod));
								}
							}
							else{
								for (auto gp = EdgeGPs.rbegin(); gp != EdgeGPs.rend(); ++gp) {
									tri->push_back(GradPaths.size());
									GradPaths.push_back(ConcatenateResample({ gp->second, *BPPtr }, NumGPPoints, {}, ResampleMethod));
// 									GradPaths.push_back(ConcatenateResample({ gp->second, *BPPtr }, NumGPPoints, GPNumPointList, ResampleMethod));
								}
							}
						}

						// Now reorder the element in the original SphereElems list.
						// tri is ordered so that the first and last elements both correspond
						// to the B type node, so reorder the original element so that the B type node
						// comes first.
						for (int ci = 0; ci < 3; ++ci) {
							if (NodeTypes[SphereElems[ti][ci]] == NETypeB) {
								vector<int> NewElem;
								for (int cj = 0; cj < 3; ++cj) {
									NewElem.push_back(SphereElems[ti][(ci + cj) % 3]);
								}
								SphereElems[ti] = NewElem;
							}
						}

						break;
					}
				}
			}

			if (UserQuit) {
				TecUtilDataLoadEnd();
				StatusDrop(AddOnID);
				return;
			}

			// Check that GPs are monotomically decreasing in rho value.
			// If one isn't, recalculate rho values along the path, as the 
			// resampled paths might be fixed this way.
			// If still not monotomically decreasing, fix that by sorting the rho
			// values: that sounds stupid, but...
			// Also save the GP for inspection.
			NumElems = SphereElemGPInds.size();
			int NumGPs = GradPaths.size();
			vector<bool> GBIsDivergent(NumElems, false),
				GPNotMonotomicallyDecreasing(NumGPs, false),
				GPNotMade(NumGPs, false),
				GPFixedAfterRecalculatingRho(NumGPs, true),
				GBIsValid(NumElems, true);
			vector<GradPath_c> TmpGPs(NumGPs);
			vector<int> BadPtNums(NumGPs);

#pragma omp parallel for schedule(dynamic)
			for (int ni = 0; ni < NumGPs; ++ni){
				if (GradPaths[ni].IsMade()) {
					GradPaths[ni].MakeRhoValuesMonotonic(&ThVolInfo[omp_get_thread_num()]);
					if (GradPaths[ni].RhoAt(0) < GradPaths[ni].RhoAt(-1)) {
						GradPaths[ni].Reverse();
					}
				}
			}

#if 1
			/*
			 *	Before moving on to integration, align all grad paths according to iso values of 
			 *	log(rho). 
			 *	This is to make sure that gradient bundle tetrahedral decomposition results in the most
			 *	regular tetrahedra possible, which can be an issue as some gradient paths are much longer
			 *	than others. 
			 *	
			 *	There will be one or more sets of rho values to use, for each cage CP coincident
			 *	with the bader atom, and one for GPs truncated at the system boundary in open systems.
			 *	Each bond path or ring surface-coincident gradient bundle has a unique set of gradient paths,
			 *	and so all gradient paths that terminate at a given cage can have its points resampled according
			 *	to the iso-rho values for that cage (or system boundary).
			 *	Rho values at bond and ring CPs will be included in the list of iso-rho values so that bond and
			 *	ring-coincident gradient paths keep their cp-coincident points.
			 *	
			 *	Start by collecting the set of unique terminating rho values within some tolerance.
			 *	When the resampling occurs, the path endpoints won't change, so it doesn't matter that we 
			 *	actually capture the rho values at the end, just that we can effectively categorize them.
			 *	Use a map of <double> -> set<double> to store the min rho values and point towards
			 *	a set of rho isovalues.
			 *	First we'll load the path terminal rho values into the map if it's not already in it,
			 *	then compute the isorho values for the unique set.
			 *	
			 *	The only "valid" terminal rho values are either the rho value used for GP truncation
			 *	or the rho values of cage CPs.
			 *	We'll only include GPs that have one of these terminal rho values
			 */

			std::set<double> ValidRhoVals({ CutoffVal });
			
			for (int i = 0; i < CPs.NumCages(); ++i){
				ValidRhoVals.insert(CPs.GetRho(CPTypeNum_Cage, i));
			}

			std::map<double, std::set<double> > RhoMinMapToIsoRhoVals;
			vector<double> GradPathMinRho(GradPaths.size(), -1.0);
			vector<bool> GradPathChecked(GradPaths.size(), false);
			double RhoHighValue;
			double RhoRelTol = 1e-1; // relative tolorance of 10%
			if (RadialSphereApprx){
				vec SphereNodeRhoVals(SphereNodes.size());
				for (int i = 0; i < SphereNodes.size(); ++i){
					SphereNodeRhoVals[i] = ValAtPointByPtr(SphereNodes[i], VolInfo, RhoPtr);
				}
				RhoHighValue = mean(SphereNodeRhoVals);

			}
			else{
				RhoHighValue = ValAtPointByPtr(CPPos, VolInfo, RhoPtr);
			}

//  			std::set<double> BondRingCPRhoVals;
//  			vector<vector<int> > GradPathCPCoincidentInds(GradPaths.size());
//  			for (int ti = 0; ti < SphereElemGPInds.size(); ++ti) {
//  				// Include edge GPs for the element if present
//  				auto GPInds = SphereElemGPInds[ti];
//  				for (int c1 = 0; c1 < 3; ++c1) {
//  					int ni = SphereElems[ti][c1];
//  					if (NodeTypes[ni] == NETypeR)
//  						ni = SphereElemGPInds[ti][c1];
//  
//  					int nj = SphereElems[ti][(c1 + 1) % 3];
//  					if (NodeTypes[nj] == NETypeR)
//  						nj = SphereElemGPInds[ti][(c1 + 1) % 3];
//  
//  					Edge TmpEdge = MakeEdge(ni, nj);
//  					if (EdgeGPMap.count(TmpEdge)) {
//  						auto EdgeGPInds = EdgeGPMap[TmpEdge];
//  						if (!EdgeGPInds.empty()) {
//  							GPInds.insert(GPInds.end(), EdgeGPInds.cbegin(), EdgeGPInds.cend());
//  						}
//  					}
//  				}
//  				for (int gi : GPInds) {
//  					if (!GradPathChecked[gi]){
//  						GradPathChecked[gi] = true;
//  						if (GradPaths[gi].IsMade()) {
//  							double TermRhoVal = GradPaths[gi].RhoAt(-1);
//  							double ClosestValidTermVal = DBL_MIN;
//  							for (auto v : ValidRhoVals){
//  								if (abs(v - TermRhoVal) < abs(ClosestValidTermVal - TermRhoVal)){
//  									ClosestValidTermVal = v;
//  								}
//  							}
//  							if (abs(TermRhoVal - ClosestValidTermVal) / TermRhoVal < RhoRelTol){
//  								GradPathMinRho[gi] = ClosestValidTermVal;
//  								RhoMinMapToIsoRhoVals.emplace(ClosestValidTermVal, std::set<double>());
//  
//  								if (ElemTypes[ti] >= NETypeR) {
//  									GradPathCPCoincidentInds[gi] = GradPaths[gi].GetCPCoincidentPoints(&CPs);
//  									for (int cpind : GradPathCPCoincidentInds[gi]) {
//  										auto pathrho = GradPaths[gi].RhoAt(cpind);
//  										bool InSet = false;
//  										for (auto cprho : BondRingCPRhoVals) {
//  											if (abs(pathrho - cprho) / pathrho < RhoRelTol) {
//  												InSet = true;
//  												break;
//  											}
//  										}
//  
//  										if (!InSet) {
//  											BondRingCPRhoVals.insert(pathrho);
//  										}
//  									}
//  								}
//  							}
//  						}
//  					}
//  				}
//  			}
//  			vector<int> GPsToDo;
//  			GPsToDo.reserve(GradPaths.size());
//  			for (int i = 0; i < GradPathMinRho.size(); ++i) {
//  				if (GradPathMinRho[i] > 0) {
//  					GPsToDo.push_back(i);
//  				}
//  			}
//  
//  			
//  
//  			/*
//  			  *	Now generate log values for each <min,max> and add bond/ring rho values.
//  			  */
//  			for (auto & m : RhoMinMapToIsoRhoVals) {
//  				auto RhoVals = conv_to<vector<double> >::from(LogSpace(m.first, RhoHighValue * 0.999, NumGPPoints - BondRingCPRhoVals.size()));
//  				m.second = std::set<double>(RhoVals.begin(), RhoVals.end());
//  				for (auto cprho : BondRingCPRhoVals) {
//  					m.second.insert(cprho);
//  				}
//  			}
//  
//  			/*
//  			 * Now resample all Gradpaths.
//  			 * First store the isorho vals in descending order as vectors that can be passed to the resample function.
//  			 */
//  			std::map<double,vector<double> > RhoMinMaxMapToIsoRhoValsVec;
//  			for (auto const & m : RhoMinMapToIsoRhoVals) {
//  				RhoMinMaxMapToIsoRhoValsVec[m.first] = vector<double>(m.second.crbegin(), m.second.crend());
//  			}
//  
//  			NumToDo = GPsToDo.size();
//  #pragma omp parallel for
//  			for (int gi = 0; gi < NumToDo; ++gi) {
//  				GradPaths[GPsToDo[gi]].ResampleByRhoVals(RhoMinMaxMapToIsoRhoValsVec[GradPathMinRho[GPsToDo[gi]]], GradPathCPCoincidentInds[GPsToDo[gi]]);
//  			}

			std::set<CPType_e> SaddleCPTypes = { CPType_Bond, CPType_Ring };
// 			std::set<double> BondRingCPRhoVals;
// 			for (int ti = 0; ti < SphereElemGPInds.size(); ++ti) {
// 				// Include edge GPs for the element if present
// 				auto GPInds = SphereElemGPInds[ti];
// 				for (int c1 = 0; c1 < 3; ++c1) {
// 					int ni = SphereElems[ti][c1];
// 					if (NodeTypes[ni] == NETypeR)
// 						ni = SphereElemGPInds[ti][c1];
// 
// 					int nj = SphereElems[ti][(c1 + 1) % 3];
// 					if (NodeTypes[nj] == NETypeR)
// 						nj = SphereElemGPInds[ti][(c1 + 1) % 3];
// 
// 					Edge TmpEdge = MakeEdge(ni, nj);
// 					if (EdgeGPMap.count(TmpEdge)) {
// 						auto EdgeGPInds = EdgeGPMap[TmpEdge];
// 						if (!EdgeGPInds.empty()) {
// 							GPInds.insert(GPInds.end(), EdgeGPInds.cbegin(), EdgeGPInds.cend());
// 						}
// 					}
// 				}
// 				for (int gi : GPInds) {
// 					if (!GradPathChecked[gi]) {
// 						GradPathChecked[gi] = true;
// 						if (GradPaths[gi].IsMade()) {
// 							if (ElemTypes[ti] > NETypeR) {
// 								for (int cpind : GradPaths[gi].GetCPCoincidentPoints(&CPs, SaddleCPTypes)) {
// 									auto pathrho = GradPaths[gi].RhoAt(cpind);
// 									bool InSet = false;
// 									for (auto cprho : BondRingCPRhoVals) {
// 										if (abs(pathrho - cprho) / pathrho < RhoRelTol) {
// 											InSet = true;
// 											break;
// 										}
// 									}
// 
// 									if (!InSet) {
// 										BondRingCPRhoVals.insert(pathrho);
// 									}
// 								}
// 							}
// 						}
// 					}
// 				}
// 			}

			std::map<int, int> TerminalCPIndToMaxLenGPIndMap;
			for (int ti = 0; ti < SphereElemGPInds.size(); ++ti) {
				// Include edge GPs for the element if present
				if (ElemTypes[ti] > NETypeR) {
					auto GPInds = SphereElemGPInds[ti];
					for (int c1 = 0; c1 < 3; ++c1) {
						int ni = SphereElems[ti][c1];
						if (NodeTypes[ni] == NETypeR)
							ni = SphereElemGPInds[ti][c1];

						int nj = SphereElems[ti][(c1 + 1) % 3];
						if (NodeTypes[nj] == NETypeR)
							nj = SphereElemGPInds[ti][(c1 + 1) % 3];

						Edge TmpEdge = MakeEdge(ni, nj);
						if (EdgeGPMap.count(TmpEdge)) {
							auto EdgeGPInds = EdgeGPMap[TmpEdge];
							if (!EdgeGPInds.empty()) {
								GPInds.insert(GPInds.end(), EdgeGPInds.cbegin(), EdgeGPInds.cend());
							}
						}
					}
					for (int gi : GPInds) {
						if (!GradPathChecked[gi]) {
							GradPathChecked[gi] = true;
							if (GradPaths[gi].IsMade()) {
								double TermCPInd = GradPaths[gi].GetStartEndCPNum(1);
								auto TermCPType = CPs.GetTypeFromTotOffset(TermCPInd);
								if (TermCPType == CPType_Cage || TermCPType == CPType_Invalid) {
									double PathLen = GradPaths[gi].GetLength();
									if (TerminalCPIndToMaxLenGPIndMap.count(TermCPInd)) {
										if (PathLen > GradPaths[TerminalCPIndToMaxLenGPIndMap[TermCPInd]].GetLength()) {
											TerminalCPIndToMaxLenGPIndMap[TermCPInd] = gi;
										}
									}
									else {
										TerminalCPIndToMaxLenGPIndMap[TermCPInd] = gi;
									}
								}
							}
						}
					}
				}
			}

			std::map<int, vector<double> > TerminalCPIndToRhoValsMap;
			for (auto r : TerminalCPIndToMaxLenGPIndMap){
				auto GPRhoPtr = GradPaths[r.second].GetRhoListPtr();
				auto RhoVals = vector<double>(GPRhoPtr->cbegin(), GPRhoPtr->cend());
				// Check for duplicate rho vals in list, and if present shift them over, protecting CP-coincident points
				auto CPPts = GradPaths[r.second].GetCPCoincidentPoints(&CPs, SaddleCPTypes);
				auto CPPtsSet = std::set<int>(CPPts.begin(), CPPts.end());
				for (int i = 2; i < RhoVals.size() - 1; ++i) {
					if (abs(RhoVals[i] - RhoVals[i-1]) < 1e-8){
						if (!CPPtsSet.count(i)){
							RhoVals[i] = (RhoVals[i] + RhoVals[i + 1]) * 0.5;
						}
						else if (!CPPtsSet.count(i-1)){
							RhoVals[i-1] = (RhoVals[i-1] + RhoVals[i-2]) * 0.5;
						}
					}
				}
				GradPaths[r.second].ReplaceRhoList(RhoVals);
				TerminalCPIndToRhoValsMap[r.first] = RhoVals;
			}

			NumGPs = GradPaths.size();
			
#ifndef _DEBUG
#pragma omp parallel for schedule(dynamic)
#endif
			for (int gi = 0; gi < NumGPs; ++gi) {
				if (GradPaths[gi].IsMade()) {
					double TermCPInd = GradPaths[gi].GetStartEndCPNum(1);
					auto TermCPType = CPs.GetTypeFromTotOffset(TermCPInd);
					if ((TermCPType == CPType_Cage || TermCPType == CPType_Invalid) && gi != TerminalCPIndToMaxLenGPIndMap[TermCPInd]) {
						auto CPPts = GradPaths[gi].GetCPCoincidentPoints(&CPs, SaddleCPTypes);
						auto RhoVals = TerminalCPIndToRhoValsMap[TermCPInd];
						if (!RhoVals.empty()) {
							for (auto cpInd : CPPts) {
								double cpVal = GradPaths[gi].RhoAt(cpInd);
								int MinInd = -1;
								double MinDiff = DBL_MAX;
								for (int i = 0; i < RhoVals.size(); ++i) {
									double tmpDiff = abs(log(RhoVals[i]) - log(cpVal));
									if (tmpDiff < MinDiff) {
										MinDiff = tmpDiff;
										MinInd = i;
									}
								}
								RhoVals[MinInd] = cpVal;
							}

							// Resample twice to better converge at the right points.
	// 						for (int i = 0; i < 2; ++i) {
							GradPaths[gi].ResampleByRhoVals(RhoVals, GradPaths[gi].GetCPCoincidentPoints(&CPs, SaddleCPTypes));
							GradPaths[gi].ReinterpolateRhoValuesFromVolume(&ThVolInfo[omp_get_thread_num()]);
							// 						}
						}
					}
				}
			}

// 			NumGPs = GradPaths.size();
// #pragma omp parallel for schedule(dynamic)
// 			for (int gi = 0; gi < NumGPs; ++gi) {
// 				if (GradPaths[gi].IsMade()) {
// 					GradPaths[gi].ReinterpolateRhoValuesFromVolume(&ThVolInfo[omp_get_thread_num()]);
// 					if (GradPaths[gi].RhoAt(0) < GradPaths[gi].RhoAt(-1)) {
// 						GradPaths[gi].Reverse();
// 					}
// 					if (abs(GradPaths[gi].RhoAt(0) - RhoHighValue) / RhoHighValue < RhoRelTol) {
// 						auto CPPts = GradPaths[gi].GetCPCoincidentPoints(&CPs, SaddleCPTypes);
// 						auto RhoVals = conv_to<vector<double> >::from(LogSpace(GradPaths[gi].RhoAt(-1), RhoHighValue, NumGPPoints));
// 						for (auto cpInd : CPPts){
// 							double cpVal = GradPaths[gi].RhoAt(cpInd);
// 							int MinInd = -1;
// 							double MinDiff = DBL_MAX;
// 							for (int i = 0; i < RhoVals.size(); ++i){
// 								double tmpDiff = abs(log(RhoVals[i]) - log(cpVal));
// 								if (tmpDiff < MinDiff){
// 									MinDiff = tmpDiff;
// 									MinInd = i;
// 								}
// 							}
// 							RhoVals[MinInd] = cpVal;
// 						}
// 						std::reverse(RhoVals.begin(), RhoVals.end());
// 
// 						// Resample twice to better converge at the right points.
// 						for (int i = 0; i < 2; ++i) {
// 							GradPaths[gi].ResampleByRhoVals(RhoVals, GradPaths[gi].GetCPCoincidentPoints(&CPs, SaddleCPTypes));
// 							GradPaths[gi].ReinterpolateRhoValuesFromVolume(&ThVolInfo[omp_get_thread_num()]);
// 						}
// 					}
// 				}
// 			}
			 

#endif

#if 0
			/*
			 *	Before moving on to integration (and GB creation), align the path representation
			 *	of all paths with that of bond path GPs.
			 *	Do this by making a FIFO queue<std::pair<GradPath*, GradPath*> >, where pair.first points to 
			 *	the GP to be aligned according to the nodes of pair.second, and the queue is filled using a 
			 *	breadth-first search starting from the set of bond GPs.
			 *	So the bond GPs' nearest neighbors will be the first to be aligned, and will use the bond
			 *	GPs as the alignment source.
			 *	The BFS will have to use the element connectivity because the GPs no longer map one-to-one 
			 *	to the sphere nodes, and the edge GP information so that the queue
			 *	will include all nodal and edge GPs, including those that have been concatenated to bond path and ring line 
			 *	segments.
			 *	When the queue is formed, we simply go through the queue, aligning each pair.first GP according
			 *	to the pair.second GP.
			 *	
			 *	When complete, this will guarantee that the tetrahedra that are eventually integrated are as regular as possible.
			 */
			std::queue<std::pair<int, std::pair<int, int> > > GPIndQueue; // GB and source/target GP indices for each GP
			vector<GradPath_c*> FinishedGPs;

				/* 
				* First, populate the queue with the GPs corresponding to the edges added 
				* after the bond point for bond path nodes.
				* 
				* Becuase all the bond-path GPs need to also be aligned, for each bond path
				* we'll check 
				* 
				* For bond path elements, the order these edges are added is inconsequential
				* (I'm pretty sure), but for bond-ring elements, need to start with the edge
				* closest to bond-ring-cage GP, since it's the path with two right angles to
				* which the others need to be aligned.
				*/

			for (int ti = 0; ti < NumElems; ++ti){
				if (ElemTypes[ti] > NETypeR) {
					// We know that the first node of the element is the bond node.
						
					int NumEdgePaths = SphereElemGPInds[ti].size() - 4;
					if (NumEdgePaths > 0) {
						if (ElemTypes[ti] == NETypeB || (NodeTypes[SphereElems[ti][1]] == NETypeR && NodeTypes[SphereElems[ti][2]] == NETypeR)) {
							// If the other two nodes of the element are BOTH ring type,
							// then we'll split them so each half is aligned to its closer
							// bond-ring-cage path.
							// This is also the behavior to follow for bond-only elements.
							// 

							// First one half, using the element's fourth node as its first
							// alignment source GP.
							// For odd number of edge paths, this pass gets one more pass
							// compared to the next pass.
							for (int i = 0; i < NumEdgePaths / 2 + 1; ++i) {
								if (GradPaths[SphereElemGPInds[ti][4 + i]].IsMade() && GradPaths[SphereElemGPInds[ti][3 + i]].IsMade()) {
									GPIndQueue.push(std::make_pair(ti, std::make_pair(SphereElemGPInds[ti][4 + i], SphereElemGPInds[ti][3 + i])));
								}
							}

							// Now the other half, in reverse order
							for (int i = 0; i < NumEdgePaths / 2; ++i) {
								int e1 = SphereElemGPInds[ti].size() - i - 1;
								int e2 = (e1 + 1) % SphereElemGPInds[ti].size();
								if (GradPaths[SphereElemGPInds[ti][e1]].IsMade() && GradPaths[SphereElemGPInds[ti][e2]].IsMade()) {
									GPIndQueue.push(std::make_pair(ti, std::make_pair(SphereElemGPInds[ti][e1], SphereElemGPInds[ti][e2])));
								}
							}
						}
						else if (NodeTypes[SphereElems[ti][1]] == NETypeR) {
							// If ONLY the second node of the element is ring type, then add
							// the post-bond point edges in reverse order (relative to their
							// order in SphereElemGPInds).
							for (int i = 0; i < NumEdgePaths; ++i) {
								int e1 = SphereElemGPInds[ti].size() - i - 1;
								int e2 = (e1 + 1) % SphereElemGPInds[ti].size();
								if (GradPaths[SphereElemGPInds[ti][e1]].IsMade() && GradPaths[SphereElemGPInds[ti][e2]].IsMade()) {
									GPIndQueue.push(std::make_pair(ti, std::make_pair(SphereElemGPInds[ti][e1], SphereElemGPInds[ti][e2])));
								}
							}
						}
						else if (NodeTypes[SphereElems[ti][2]] == NETypeR) {
							// If ONLY the third node is ring type, then add the post-bond
							// point edges in original order.
							for (int i = 0; i < NumEdgePaths; ++i) {
								if (GradPaths[SphereElemGPInds[ti][4 + i]].IsMade() && GradPaths[SphereElemGPInds[ti][3 + i]].IsMade()) {
									GPIndQueue.push(std::make_pair(ti, std::make_pair(SphereElemGPInds[ti][4 + i], SphereElemGPInds[ti][3 + i])));
								}
							}
						}
					}
				}
			}

			/*
				* Now start the breadth-first search of elements (since we are
				* no longer able to generate node connectivity, as we're not
				* interested in nodes but rather GPs), using a FIFO element queue first
				* populated with the bond path elements.
				* The search will grow out from each bond path on the sphere simultaneously,
				* guaranteeing that, to the extent that the sphere meshing is distributed in
				* a physically reasonable way, we'll be aligning each GP according to the
				* "closest" bond path GPs.
				* 
				* Because we're searching by elements, but are thinking in terms of nodal
				* connectivity, it's less straightforward to keep track of which
				* nodes and edges have already been visited.
				* We'll keep track of nodes and edges separately, so that we can discern
				* between edges that have been added, and those for which only one, the other,
				* or both nodes have been added but not the corresponding edge GPs.
				* By only operating on edges where one or more nodes have already been processed,
				* we'll effectively be moving the BFS out from bond path nodes according to 
				* nodal connectivity.
				* 
				* How an edge is identified depends on the type of element:
				*	*	For cage- and ring-type elements, the edge can be identified simply by 
				*		the SphereElemGPInds[ti] indices. 
				*		This still recovers shared edges with nodes of type cage-cage and cage-ring, 
				*		and preserves the correct mapping of node and edge GPs along ring surfaces.
				*	
				*	*	For bond type elements, the bond node itself is no longer pointed to by
				*		SphereElemGPInds[ti], which instead points to the two new gradient paths
				*		made by concatenating the existing bond path (and ring line if bond-ring element).
				*		To recover shared edges with a bond type node, they will be identified using
				*		the SphereElems[ti][0] node index, which then reveals the correspondance to
				*		the other elements that share the same bond node.
				*		Because we're only using edges as a means of checking whether one has been 
				*		visited already in the search, no information is lost here, and we still use
				*		SphereElemGPInds[ti] for identifying whether nodes have been visited and
				*		recovering GPs.
				*		
				* 
				* The search will move out according to sphere elements, but we'll be concerned
				* more with edges.
				* There are four possibilities when encountering an edge:
				* 1. One node in the edge will have already been added, and the corresponding edge GPs will
				*		be added to the GP queue in order of increasing distance from the already added node.
				* 2. Both nodes in the edge have been added, but the GPs corresponding to the edge have not.
				*		This is what happens when the searches starting from two different bond paths collide.
				*		Here we'll split the corresponding edge GPs, adding have in order of increasing distance
				*		from their respective node.
				* 3. Both nodes and the corresponding edge GPs have been added: ignore and move on.
				* 4. Neither nodes have been added: ignore and they'll be added in a subsequent iteration.
				*/

			// Prepare queue and sets for workspace
			std::queue<int> ElemIndQueue;
			std::set<Edge> EdgeWasVisited;
			std::set<int> NodeWasVisited, ElemWasVisited;

			// First add the bond(-ring) type elements to the queue
			// and bond-type nodes to the node set, setting up the
			// initial conditions for the BFS to progress.
			for (int ti = 0; ti < NumElems; ++ti){
				if (ElemTypes[ti] > NETypeR){
					// First and last nodes of element points to the bond node
					NodeWasVisited.insert(SphereElemGPInds[ti][0]);
					NodeWasVisited.insert(SphereElemGPInds[ti].back());
					ElemWasVisited.insert(ti);
					if (SphereElemGPInds[ti].size() > 3){
						bool ElemMade = true;
						for (int i = 0; i < SphereElemGPInds[ti].size() && ElemMade; ++i) {
							ElemMade = GradPaths[SphereElemGPInds[ti][i]].IsMade();
						}
						if (ElemMade) {
							ElemIndQueue.push(ti);
						}
					}
				}
			}

			// Get element connectivity
			vector<vector<int> > SphereElemConnectivity;
			GetTriElementConnectivityList(&SphereElems, SphereElemConnectivity, 1);

			// Start the element BFS
			while (!ElemIndQueue.empty()){
				// Get front element
				int ti = ElemIndQueue.front();
				ElemIndQueue.pop();

				// Add neighboring elems
				for (int tj : SphereElemConnectivity[ti]){
					if (!ElemWasVisited.count(tj)){
						ElemWasVisited.insert(tj);
						bool ElemMade = true;
						for (int i = 0; i < SphereElemGPInds[tj].size() && ElemMade; ++i){
							ElemMade = GradPaths[SphereElemGPInds[tj][i]].IsMade();
						}
						if (ElemMade) {
							ElemIndQueue.push(tj);
						}
					}
				}

				// Check elem edges for cases 1, 2, or 3, as indicated in the above comments
				for (int ci = 0; ci < 3; ++ci){
					int cj = (ci + 1) % 3;
					Edge e = MakeEdge(SphereElemGPInds[ti][ci], SphereElemGPInds[ti][cj]);
					Edge egp = MakeEdge(SphereElemGPInds[ti][ci], SphereElemGPInds[ti][cj]);
					if (ElemTypes[ti] > NETypeR) {
						cj = (ci + 1);
						e = MakeEdge(
							NodeTypes[SphereElems[ti][ci]] == NETypeB ? SphereElems[ti][ci] : SphereElemGPInds[ti][ci],
							NodeTypes[SphereElems[ti][cj % 3]] == NETypeB ? SphereElems[ti][cj % 3] : SphereElemGPInds[ti][cj]);
					}
					if (!EdgeWasVisited.count(e) && GradPaths[egp.first].IsMade() && GradPaths[egp.second].IsMade()){
						// Edge wasn't visited. See which case from comments above; 1, 2, or 4
						if (!NodeWasVisited.count(SphereElemGPInds[ti][ci]) && !NodeWasVisited.count(SphereElemGPInds[ti][cj])){
							// case 4; continue
							continue;
						}
						else {
							EdgeWasVisited.insert(e);
							auto EdgeGPInds = EdgeGPMap[egp];
							int NumEdgeGPs = EdgeGPInds.size();
							GradPath_c *GP1 = &GradPaths[SphereElemGPInds[ti][ci]], *GP2 = &GradPaths[SphereElemGPInds[ti][cj]];
							if (GP1->IsMade() && GP2->IsMade()) {
								if (NumEdgeGPs > 0) {
									double DistFront = 0.0, DistBack = 0.0;
									int ind1, ind2;
									vec3 vec1, vec2;
									if (GP1->GetMaxSeparationMidpointFromOtherGPRhoBased(GradPaths[EdgeGPInds.front()], vec1, ind1, ind2, vec2, DistFront)
										&& GP1->GetMaxSeparationMidpointFromOtherGPRhoBased(GradPaths[EdgeGPInds.back()], vec1, ind1, ind2, vec2, DistBack))
									{
										bool ciCloser = (DistFront < DistBack);
										if (!NodeWasVisited.count(SphereElemGPInds[ti][ci])) {
											// SphereElemGPInds[ti][ci] hasn't been added but SphereElemGPInds[ti][cj] has, so add edge
											// GPs and SphereElemGPInds[ti][ci] to the gp queue starting from SphereElemGPInds[ti][cj].
											if (ciCloser) {
												// Add in reverse order
												for (int i = 0; i < NumEdgeGPs; ++i) {
													if (i == 0) {
														GPIndQueue.push(std::make_pair(ti, std::make_pair(EdgeGPInds.back(), SphereElemGPInds[ti][cj])));
													}
													else {
														GPIndQueue.push(std::make_pair(ti, std::make_pair(EdgeGPInds[NumEdgeGPs - i - 1], EdgeGPInds[NumEdgeGPs - i])));
													}
												}
												// Then add ci
												GPIndQueue.push(std::make_pair(ti, std::make_pair(SphereElemGPInds[ti][ci], EdgeGPInds.front())));
											}
											else {
												// Add in forward order
												for (int i = 0; i < NumEdgeGPs; ++i) {
													if (i == 0) {
														GPIndQueue.push(std::make_pair(ti, std::make_pair(EdgeGPInds.front(), SphereElemGPInds[ti][cj])));
													}
													else {
														GPIndQueue.push(std::make_pair(ti, std::make_pair(EdgeGPInds[i], EdgeGPInds[i - 1])));
													}
												}
												// Then add ci
												GPIndQueue.push(std::make_pair(ti, std::make_pair(SphereElemGPInds[ti][ci], EdgeGPInds.back())));
											}
											NodeWasVisited.insert(SphereElemGPInds[ti][ci]);
										}
										if (!NodeWasVisited.count(SphereElemGPInds[ti][cj])) {
											// SphereElemGPInds[ti][cj] hasn't been added but SphereElemGPInds[ti][ci] has, so add edge
											// GPs and SphereElemGPInds[ti][cj] to the gp queue starting from SphereElemGPInds[ti][ci].
											if (ciCloser) {
												// Add in forward order
												for (int i = 0; i < NumEdgeGPs; ++i) {
													if (i == 0) {
														GPIndQueue.push(std::make_pair(ti, std::make_pair(EdgeGPInds.front(), SphereElemGPInds[ti][ci])));
													}
													else {
														GPIndQueue.push(std::make_pair(ti, std::make_pair(EdgeGPInds[i], EdgeGPInds[i - 1])));
													}
												}
												// Then add ci
												GPIndQueue.push(std::make_pair(ti, std::make_pair(SphereElemGPInds[ti][cj], EdgeGPInds.back())));
											}
											else {
												// Add in reverse order
												for (int i = 0; i < NumEdgeGPs; ++i) {
													if (i == 0) {
														GPIndQueue.push(std::make_pair(ti, std::make_pair(EdgeGPInds.back(), SphereElemGPInds[ti][ci])));
													}
													else {
														GPIndQueue.push(std::make_pair(ti, std::make_pair(EdgeGPInds[NumEdgeGPs - i - 1], EdgeGPInds[NumEdgeGPs - i])));
													}
												}
												// Then add ci
												GPIndQueue.push(std::make_pair(ti, std::make_pair(SphereElemGPInds[ti][cj], EdgeGPInds.front())));
											}
											NodeWasVisited.insert(SphereElemGPInds[ti][cj]);
										}
										else {
											// both ci and cj have been added already, so split the edge gps as we did before.
											if (ciCloser) {
												// ci get's low half of edge gps
												for (int i = 0; i < NumEdgeGPs / 2 + 1; ++i) {
													if (i == 0) {
														GPIndQueue.push(std::make_pair(ti, std::make_pair(EdgeGPInds[i], SphereElemGPInds[ti][ci])));
													}
													else {
														GPIndQueue.push(std::make_pair(ti, std::make_pair(EdgeGPInds[i], EdgeGPInds[i - 1])));
													}
												}
												// cj get's upper half of edge gps
												for (int i = 0; i < NumEdgeGPs / 2; ++i) {
													if (i == 0) {
														GPIndQueue.push(std::make_pair(ti, std::make_pair(EdgeGPInds[NumEdgeGPs - i - 1], SphereElemGPInds[ti][cj])));
													}
													else {
														GPIndQueue.push(std::make_pair(ti, std::make_pair(EdgeGPInds[NumEdgeGPs - i - 1], EdgeGPInds[NumEdgeGPs - i - 2])));
													}
												}
											}
											else {
												// cj get's low half of edge gps
												for (int i = 0; i < NumEdgeGPs / 2 + 1; ++i) {
													if (i == 0) {
														GPIndQueue.push(std::make_pair(ti, std::make_pair(EdgeGPInds[i], SphereElemGPInds[ti][cj])));
													}
													else {
														GPIndQueue.push(std::make_pair(ti, std::make_pair(EdgeGPInds[i], EdgeGPInds[i - 1])));
													}
												}
												// ci get's upper half of edge gps
												for (int i = 0; i < NumEdgeGPs / 2; ++i) {
													if (i == 0) {
														GPIndQueue.push(std::make_pair(ti, std::make_pair(EdgeGPInds[NumEdgeGPs - i - 1], SphereElemGPInds[ti][ci])));
													}
													else {
														GPIndQueue.push(std::make_pair(ti, std::make_pair(EdgeGPInds[NumEdgeGPs - i - 1], EdgeGPInds[NumEdgeGPs - i - 2])));
													}
												}
											}
										}
									}
								}
								else {
									// No edge GPs between node GPs, so only need to look for the case where
									// one node GP was added and the other wasn't, and the unadded GP will
									// point to the added as its parent.
									// Don't care about the cases where neither or both GPs are already added.
									if (!NodeWasVisited.count(SphereElemGPInds[ti][ci]) && NodeWasVisited.count(SphereElemGPInds[ti][cj])) {
										// Add ci to queue pointing to cj
										GPIndQueue.push(std::make_pair(ti, std::make_pair(SphereElemGPInds[ti][ci], SphereElemGPInds[ti][cj])));
										NodeWasVisited.insert(SphereElemGPInds[ti][ci]);
									}
									else if (NodeWasVisited.count(SphereElemGPInds[ti][ci]) && !NodeWasVisited.count(SphereElemGPInds[ti][cj])) {
										// Add cj to queue pointing to ci
										GPIndQueue.push(std::make_pair(ti, std::make_pair(SphereElemGPInds[ti][cj], SphereElemGPInds[ti][ci])));
										NodeWasVisited.insert(SphereElemGPInds[ti][cj]);
									}
								}
							}
						}
					}// else case 3; continue
				}

					
			}

			/*
				*	The GP queue is now full of <GradPath*, GradPath*> pairs, where
				*	pair.first will have it's points aligned with those of pair.second,
				*	in the optimal order such that the placement of nodes along
				*	bond-path GPs are applied to all other GPs.
				*	Rather than having misalignment occur at bond and ring type nodes,
				*	it will now occur at the maximum "distance" (in terms of mesh edge steps)
				*	from bond path nodes.
				*	
				*	At this point we simply move through the GP queue, aligning paths as we go,
				*	and the ordering is such that the source (pair.second) path is guaranteed
				*	to have had its points already aligned with its "parent" GP.
				*/
			vector<vector<int> > GPConstraintInds(GradPaths.size());
			while (!GPIndQueue.empty()){
				auto GPInds = GPIndQueue.front();
				GPIndQueue.pop();

// 					auto sourceGP = GradPaths[GPInds.second.second];
// 					std::stringstream ss;
// 					ss << "GB " << GPInds.first << ": alignment source gp " << GPInds.second.second << " --> " << GPInds.second.first;
// 					sourceGP.SaveAsOrderedZone(ss.str().c_str(), Green_C);

				if (GPConstraintInds[GPInds.second.second].empty()) {
					GPConstraintInds[GPInds.second.second] = GradPaths[GPInds.second.second].GetCPCoincidentPoints(&CPs);
				}
				GPConstraintInds[GPInds.second.first] = GPConstraintInds[GPInds.second.second];

// 					vector<vec3> constraintPoints;
// 					for (int i : GPConstraintInds[GPInds.second.second]){
// 						constraintPoints.push_back(sourceGP[i]);
// 					}
// 					ss.str("");
// 					ss << "GB " << GPInds.first << ": constraint points for gp " << GPInds.second.second;
// 					SaveVec3VecAsScatterZone(constraintPoints, ss.str().c_str());
// 
// 					ss.str("");
// 					ss << "GB " << GPInds.first << ": target gp before " << GPInds.second.first << " <-- " << GPInds.second.second;
// 					GradPaths[GPInds.second.first].SaveAsOrderedZone(ss.str().c_str(), Cyan_C);

				vector<int> ProtectedPointInds = GradPaths[GPInds.second.first].GetCPCoincidentPoints(&CPs);
				GradPaths[GPInds.second.first].AlignToOtherPath(GradPaths[GPInds.second.second], GPConstraintInds[GPInds.second.second], ProtectedPointInds);
				GradPaths[GPInds.second.first].ReinterpolateRhoValuesFromVolume(&VolInfo, &RhoPtr);

// 					ss.str("");
// 					ss << "GB " << GPInds.first << ": target gp after " << GPInds.second.first << " <-- " << GPInds.second.second;
// 					GradPaths[GPInds.second.first].SaveAsOrderedZone(ss.str().c_str(), Blue_C);
			}
#endif

			/*
			 * We now have all the gradient paths completed in a single vector,
			 * and a representation of the gradient bundles in terms of the
			 * gradient path indices in the correct order.
			 * Now simply loop over each triangle, creating a gradient bundle for it.
			 */

// 			if (SaveGPs)
				GPPtrList.resize(NumElems);


			GPNotMonotomicallyDecreasing = vector<bool>(GradPaths.size(), false);

			TmpString = ProgressStr2.str() + "Integrating";

			for (int i = 0; i < IntVarPtrs.size(); ++i) {
				IsOk = IntVarPtrs[i].InitializeReadPtr(VolZoneNum, IntVarNumList[i]);
			}

			Time1 = high_resolution_clock::now();

	#ifdef SUPERDEBUG
			NumToDo = ElemTodo.size();
			NumCompleted = 0;
// 			for (int ti : ElemTodo) {
// #pragma omp parallel for schedule(dynamic)
			for (int tii = 0; tii < NumToDo; ++tii){
				int ti = ElemTodo[tii];
#else
			NumToDo = 0;
			for (int ti = 0; ti < NumElems; ++ti)
				NumToDo += int(IntVals[ti][DensityVarNumInIntList] == 0.0);
			NumCompleted = 0;
			Time1 = high_resolution_clock::now();
	#pragma omp parallel for schedule(dynamic)
			for (int ti = 0; ti < NumElems; ++ti) {
	#endif
				if (omp_get_thread_num() == 0) {
					UserQuit = !StatusUpdate(NumCompleted, NumToDo, TmpString, AddOnID, Time1);
	#pragma omp flush (UserQuit)
				}
	#pragma omp flush (UserQuit)
				if (UserQuit) {
#pragma omp critical
					{
						if (omp_get_thread_num() == 0 && !TecUtilDialogMessageBox("Cancel run? All progress will be lost if you press \"Yes.\"", MessageBoxType_YesNo)) {
							StatusLaunch(TmpString, AddOnID, TRUE, TRUE, TRUE);
							StatusUpdate(NumCompleted, NumToDo, TmpString, AddOnID, Time1);
							UserQuit = false;
#pragma omp flush (UserQuit)
						}
#pragma omp flush (UserQuit)
					}
				}

				if (!UserQuit && IntVals[ti][DensityVarNumInIntList] == 0.0) {
#pragma omp atomic
					NumCompleted++;

					vector<GradPath_c const *> GPPtrs;
// 					for (auto ni : SphereElemGPInds[ti]) {
// 						GradPath_c const *GP1 = &GradPaths[ni];
// 						GPPtrs.push_back(GP1);
// 
// 						if (GP1->IsMade()) {
// 							for (int i = 0; i < GP1->GetCount() - 1; ++i) {
// 								if (GP1->RhoAt(i) < GP1->RhoAt(i + 1)) {
// 									GPNotMonotomicallyDecreasing[ni] = true;
// 									GBIsValid[ti] = false;
// 									break;
// 								}
// 							}
// 						}
// 						else {
// 							GPNotMade[ni] = true;
// 							GBIsValid[ti] = false;
// 						}
// 					}

// 					if (ElemTypes[ti] != NETypeB) {
// 						for (int tj = 0; tj < 3; ++tj) {
// 							int ni = SphereElems[ti][tj];
// 							auto const GP1 = &GradPaths[SphereElemGPInds[ti][tj]];
// 
// 							GPPtrs.push_back(GP1);
// 
// 							// now add and edge GPs between GP1 and the next GP
// 							int nj = SphereElems[ti][(tj + 1) % 3];
// 
// 							Edge TmpEdge = MakeEdge(ni, nj);
// 							if (EdgeGPMap.count(TmpEdge)) {
// 								auto EdgeGPInds = EdgeGPMap[TmpEdge];
// 								if (ni < nj) {
// 									for (auto ei = EdgeGPInds.cbegin(); ei != EdgeGPInds.cend(); ei++) {
// 										GPPtrs.push_back(&GradPaths[*ei]);
// 									}
// 								}
// 								else {
// 									for (auto ei = EdgeGPInds.crbegin(); ei != EdgeGPInds.crend(); ei++) {
// 										GPPtrs.push_back(&GradPaths[*ei]);
// 									}
// 								}
// 							}
// 						}
// 					}
// 					else{
// 						for (int tj = 0; tj < 3; ++tj) {
// 							int ni = SphereElems[ti][tj];
// 							auto const GP1 = &GradPaths[SphereElemGPInds[ti][tj]];
// 
// 							GPPtrs.push_back(GP1);
// 
// 							// now add and edge GPs between GP1 and the next GP
// 							int nj = SphereElems[ti][(tj + 1) % 3];
// 
// 							Edge TmpEdge = MakeEdge(ni, nj);
// 							if (EdgeGPMap.count(TmpEdge)) {
// 								auto EdgeGPInds = EdgeGPMap[TmpEdge];
// 								if (ni < nj) {
// 									for (auto ei = EdgeGPInds.cbegin(); ei != EdgeGPInds.cend(); ei++) {
// 										GPPtrs.push_back(&GradPaths[*ei]);
// 									}
// 								}
// 								else {
// 									for (auto ei = EdgeGPInds.crbegin(); ei != EdgeGPInds.crend(); ei++) {
// 										GPPtrs.push_back(&GradPaths[*ei]);
// 									}
// 								}
// 							}
// 						}
// 
// 						for (int tj = 3; tj < SphereElemGPInds[ti].size(); ++tj){
// 							GPPtrs.push_back(&GradPaths[SphereElemGPInds[ti][tj]]);
// 						}
// 					}

					for (int tj = 0; tj < 3; ++tj) {
						int ni = SphereElems[ti][tj];
						if (NodeTypes[ni] == NETypeR)
							ni = SphereElemGPInds[ti][tj];
						auto const GP1 = &GradPaths[SphereElemGPInds[ti][tj]];

						GPPtrs.push_back(GP1);

						// now add and edge GPs between GP1 and the next GP
						int nj = SphereElems[ti][(tj + 1) % 3];
						if (NodeTypes[nj] == NETypeR)
							nj = SphereElemGPInds[ti][(tj + 1) % 3];

						Edge TmpEdge = MakeEdge(ni, nj);

						if (EdgeGPMap.count(TmpEdge)) {

							auto egp = GPIndEdgeToNodeIndEdgeMap[TmpEdge].first;

							auto EdgeGPInds = EdgeGPMap[TmpEdge];
							double DistFront = 0.0, DistBack = 0.0;
							int ind1, ind2;
							vec3 vec1, vec2;
							if (!EdgeGPInds.empty() && GP1->GetMaxSeparationMidpointFromOtherGPRhoBased(GradPaths[EdgeGPInds.front()], vec1, ind1, ind2, vec2, DistFront)
								&& GP1->GetMaxSeparationMidpointFromOtherGPRhoBased(GradPaths[EdgeGPInds.back()], vec1, ind1, ind2, vec2, DistBack))
							{
								if (DistFront < DistBack) {
									for (auto ei = EdgeGPInds.cbegin(); ei != EdgeGPInds.cend(); ei++) {
										if (*ei >= GradPaths.size()) {
#pragma omp critical(ForCheck)
											{
												TecUtilDialogErrMsg("edge GP index out of bounds (GB creation forward)");
											}
										}
										GPPtrs.push_back(&GradPaths[*ei]);
									}
								}
								else{
									for (auto ei = EdgeGPInds.crbegin(); ei != EdgeGPInds.crend(); ei++) {
										if (*ei >= GradPaths.size()) {
#pragma omp critical(BackCheck)
											{
												TecUtilDialogErrMsg("edge GP index out of bounds (GB creation backward)");
											}
										}
										GPPtrs.push_back(&GradPaths[*ei]);
									}
								}
							}
// 							if (ni == egp.first) {
// // 								if (ni < nj) {
// 								for (auto ei = EdgeGPInds.cbegin(); ei != EdgeGPInds.cend(); ei++) {
// 									if (*ei >= GradPaths.size()){
// #pragma omp critical(ForCheck)
// 										{
// 											TecUtilDialogErrMsg("edge GP index out of bounds (GB creation forward)");
// 										}
// 									}
// 									GPPtrs.push_back(&GradPaths[*ei]);
// 								}
// 							}
// 							else {
// 								for (auto ei = EdgeGPInds.crbegin(); ei != EdgeGPInds.crend(); ei++) {
// 									if (*ei >= GradPaths.size()) {
// #pragma omp critical(BackCheck)
// 										{
// 											TecUtilDialogErrMsg("edge GP index out of bounds (GB creation backward)");
// 										}
// 									}
// 									GPPtrs.push_back(&GradPaths[*ei]);
// 								}
// 							}
						}
					}

					// This loop only applies to B type elements, where the edge GPs
					// going from the bond point are referenced in SphereElemGPInds.
					for (int tj = 3; tj < SphereElemGPInds[ti].size(); ++tj) {
						if (SphereElemGPInds[ti][tj] >= GradPaths.size()) {
#pragma omp critical(BondCheck)
							{
								TecUtilDialogErrMsg("edge GP index out of bounds (GB creation forward)");
							}
						}
						GPPtrs.push_back(&GradPaths[SphereElemGPInds[ti][tj]]);
					}

					// Check that no GPs are deviating and that rho for GPs is monotonically decreasing
					for (int gpi = 0; gpi < GPPtrs.size(); ++gpi) {
						GradPath_c const *GP1 = GPPtrs[gpi], *GP2 = GPPtrs[(gpi + 1) % GPPtrs.size()];
						if (GP1->IsMade() && GP2->IsMade()) {
							if (Distance(GP1->XYZAt(-1), GP2->XYZAt(-1)) > 0.5 && dot(GP1->XYZAt(-1) - GP1->XYZAt(-2), GP2->XYZAt(-1) - GP2->XYZAt(-2)) < DivergentGPMinTerminalDotProduct) {
								GBIsDivergent[ti] = true;
								GBIsValid[ti] = false;
							}
						}
					}

// 					if (SaveGPs && GPPtrList[ti].empty())
					if (GPPtrList[ti].empty())
						GPPtrList[ti] = GPPtrs;

					if (GBIsValid[ti]) {
						if (ElemTypes[ti] >= NETypeB) {
							vec3 BondCPPos;
							for (int ni : SphereElems[ti]) {
								if (NodeTypes[ni] >= NETypeB) {
									BondCPPos = IntBondPathSegments[ni].RhoAt(0) > IntBondPathSegments[ni].RhoAt(-1) ? IntBondPathSegments[ni].XYZAt(0) : IntBondPathSegments[ni].XYZAt(-1);
								}

							}
  							NewIntegrateUsingIsosurfaces2(GPPtrs, IntResolution, ThVolInfo[omp_get_thread_num()], IntVarPtrs, IntVals[ti], &BondCPPos, &AddOnID);
						}
						else {
  							NewIntegrateUsingIsosurfaces2(GPPtrs, IntResolution, ThVolInfo[omp_get_thread_num()], IntVarPtrs, IntVals[ti], nullptr, &AddOnID);
						}

					}


					if (RadialSphereApprx) {
						/* 
						 * Now do sphere integration if RadialSphereApprx.
						 * Much simpler as there are no edge GPs
						 */

						vector<const GradPath_c *> InnerGPPtrs(3);
						for (int tj = 0; tj < 3; ++tj){
							InnerGPPtrs[tj] = &InnerGradPaths[SphereElems[ti][tj]];
						}

						NewIntegrateUsingIsosurfaces2(InnerGPPtrs, IntResolution, ThVolInfo[omp_get_thread_num()], IntVarPtrs, SphereIntVals[ti], nullptr, &AddOnID);
					}
				}
			}

			


			for (auto & IntVarPtr : IntVarPtrs) {
				IntVarPtr.Close();
			}

			DivergentGBs.clear();
			NotMonotonicallyDecreasingGPs.clear();
			NotMadeGPs.clear();
			InvalidGBs.clear();

			for (int ni = 0; ni < GradPaths.size(); ++ni){
				if (GPNotMade[ni]) NotMadeGPs.push_back(ni);
				if (GPNotMonotomicallyDecreasing[ni]) NotMonotonicallyDecreasingGPs.push_back(ni);
			}
			for (int ti = 0; ti < NumElems; ++ti){
				if (GBIsDivergent[ti]) DivergentGBs.push_back(ti);
				if (!GBIsValid[ti]) InvalidGBs.push_back(ti);
			}

			if (UserQuit) {
				TecUtilDataLoadEnd();
				StatusDrop(AddOnID);
				return;
			}

			if (RadialSphereApprx) {
				Sphere = FESurface_c(SphereNodes, SphereElems);
				vec SphereTriangleAreas(Sphere.TriSphereElemSolidAngles());
				double TotalArea = sum(SphereTriangleAreas);
				vec SphereTriangleAreaFactors = SphereTriangleAreas / TotalArea;
				SphereTotalIntVals = vector<double>(SphereIntVals[0].size(), 0.0);
				for (int ti = 0; ti < SphereElems.size(); ++ti) {
					for (int iVar = 0; iVar < SphereTotalIntVals.size(); ++iVar) {
						SphereTotalIntVals[iVar] += SphereIntVals[ti][iVar];
					}
				}
				// Don't actually put the values back into the GBs; the sphere total and GB solid angles are saved, which can be used to get these values if desired
				for (int ti = 0; ti < SphereElems.size(); ++ti){
					for (int iVar = 0; iVar < SphereTotalIntVals.size(); ++iVar) {
						IntVals[ti][iVar] += SphereTotalIntVals[iVar] * SphereTriangleAreaFactors[ti];
					}
				}
			}


			TotalDensity = 0.0;
			for (auto const & elemIntVals : IntVals)
				TotalDensity += elemIntVals[DensityVarNumInIntList];


			for (int ti = 0; ti < SphereElems.size(); ++ti){
				if (IntVals[ti][DensityVarNumInIntList] > TargetDensityPerGB
					&& ElemSubdivisionLevel[ti] < MaxGBSubdivisionLevel)
					ElemsToSubdivide.push(ti);
			}
			
// 			if (!ElemsToSubdivide.empty() && SaveGBs){
// 				Sphere = FESurface_c(SphereNodes, SphereElems);
// 				Sphere.SaveAsTriFEZone(XYZVarNums, "Pass " + to_string(GBAIter) + ": sphere: " + to_string(ElemsToSubdivide.size()) + " elems subdivided");
// 
// 				auto tmpQueue = ElemsToSubdivide;
// 
// 				while (!tmpQueue.empty()){
// 					int ti = tmpQueue.front();
// 					if (GradientBundles[ti].IsMade()){
// 						GradientBundles[ti].SaveAsTriFEZone(XYZVarNums, "Pass " + to_string(GBAIter) + ": elem " + to_string(ti) + ": rho = " + to_string(IntVals[ti][DensityVarNumInIntList]));
// 					}
// 					tmpQueue.pop();
// 				}
// 			}

			if (NumElementsToSubdivide > 0 && ElemsToSubdivide.size() > NumElementsToSubdivide)
				NumElementsToSubdivide = 0;
			else
				NumElementsToSubdivide = ElemsToSubdivide.size();

			break;
		} while (NumElementsToSubdivide > 0);

		if (!DivergentGBs.empty()){
// 			for (int ti : DivergentGBs){
// 				if (GradientBundles[ti].IsMade()){
// 					GradientBundles[ti].SaveAsTriFEZone(XYZVarNums, "Divergent GB for elem " + to_string(ti));
// 				}
// 			}
			string s = IntVectorToRangeString(DivergentGBs);
			TecUtilDialogErrMsg(string("Divergent gradient bundles found for element(s): " + s).c_str());
		}

		if (!InvalidGBs.empty()) {
			// 			for (int ti : DivergentGBs){
			// 				if (GradientBundles[ti].IsMade()){
			// 					GradientBundles[ti].SaveAsTriFEZone(XYZVarNums, "Divergent GB for elem " + to_string(ti));
			// 				}
			// 			}
			string s = IntVectorToRangeString(InvalidGBs);
			TecUtilDialogErrMsg(string("Invalid gradient bundles found for element(s): " + s).c_str());
		}

		if (!NotMadeGPs.empty()) {
			string s = IntVectorToRangeString(NotMadeGPs);
			TecUtilDialogErrMsg(string("Gradient paths not made for gradient path(s): " + s).c_str());
		}

		if (!NotMonotonicallyDecreasingGPs.empty()) {
  			for (int ni : NotMonotonicallyDecreasingGPs){
  				if (GradPaths[ni].IsMade()){
  					GradPaths[ni].SaveAsOrderedZone("Not monotomically decreasing GP " + to_string(ni));
  				}
  			}
			string s = IntVectorToRangeString(NotMonotonicallyDecreasingGPs);
			TecUtilDialogErrMsg(string("Gradient paths not monotonically decreasing for gradient path(s): " + s).c_str());
		}

		if (std::find(IntVarNameList.begin(), IntVarNameList.end(), "Volume") == IntVarNameList.end())
			IntVarNameList.emplace_back("Volume");

		int VolumeNumInVarList = IntVals[0].size() - 1;

		int TeNumInVarList = -1;
		for (int vi = 0; vi < IntVarNameList.size(); ++vi) {
			auto v = IntVarNameList[vi];
			if (v.find("inetic") != string::npos) {
				TeNumInVarList = vi;
				IntVarNameList.emplace_back("Te Per Electron");
				break;
			}
		}

		IntVarNameList.emplace_back("Solid Angle (alpha)");
		IntVarNameList.emplace_back("Volume Fraction");
		string ElemName = SplitString(NucleusName, " ", true, true).front();
		int AtomicNumber = SearchVectorForString(ElementSymbolList, ElemName, false) + 1;
		if (AtomicNumber > 0) {
			IntVarNameList.emplace_back("Area Bader Charge");
			IntVarNameList.emplace_back("Valence Density");
			IntVarNameList.emplace_back("Valence Bader Charge");
			if (TeNumInVarList >= 0)
				IntVarNameList.emplace_back("Te Per Valence Electron");
			IntVarNameList.emplace_back("Area Condensed Charge");

			IntVarNameList.emplace_back("Volume Bader Charge");
// 			IntVarNameList.emplace_back("Volume Valence Density");
// 			IntVarNameList.emplace_back("Volume Valence Bader Charge");
// 			if (TeNumInVarList >= 0)
// 				IntVarNameList.emplace_back("Volume Te Per Valence Electron");
			IntVarNameList.emplace_back("Volume Condensed Charge");
		}
// 		IntVarNameList.emplace_back("P over N");
// 		IntVarNameList.emplace_back("P/(N alpha)"); // shannon condensed entropy
// 		IntVarNameList.emplace_back("S[P]"); // shannon condensed entropy
// // 		IntVarNameList.emplace_back("S[N alpha P]"); // shannon condensed entropy
// // 		IntVarNameList.emplace_back("S[N P]"); // shannon condensed entropy
// // 		IntVarNameList.emplace_back("S[alpha P]"); // shannon condensed entropy
// // 		IntVarNameList.emplace_back("S[P / (N alpha)]"); // shannon condensed entropy
// // 		IntVarNameList.emplace_back("S[P / N]"); // shannon condensed entropy
// 		IntVarNameList.emplace_back("S[P / alpha]"); // shannon condensed entropy
// 		IntVarNameList.emplace_back("S[alpha]"); // shannon condensed entropy
// 		IntVarNameList.emplace_back("S[P] / S[alpha]"); // shannon condensed entropy over max entropy
// 		IntVarNameList.emplace_back("S[alpha] - S[P]"); // shannon condensed entropy over max entropy

		/*
		 *	Now collect and save integration values
		 */
		NumZones = TecUtilDataSetGetNumZones();
		vector<string> NewVarNames;
		vector<int> NewVarNums;
		vector<FieldDataType_e> ZoneDataTypes(NumZones, FieldDataType_Bit);


		for (int i = 1; i <= NumZones; ++i) {
			if (AuxDataZoneHasItem(i, CSMAuxData.GBA.ZoneType)) {
				ZoneDataTypes[i - 1] = FieldDataType_Float;
			}
		}
		vector<ValueLocation_e> ZoneDataLocs(NumZones, ValueLocation_CellCentered);
		vector<string> VarNameAppends = { "I", "N", "S" };
		for (int i = 0; i < IntVarNameList.size() && IsOk; ++i) {
			for (int j = 0; j < 3; ++j) {
				string TmpStr;
				for (int k = 0; k <= j; ++k)
					TmpStr += VarNameAppends[k];
				NewVarNames.push_back(TmpStr + ": " + IntVarNameList[i]);
				NewVarNums.push_back(VarNumByName(NewVarNames.back()));
				if (NewVarNums[i * 3 + j] < 0) {
					ArgList args;
					args.appendString(SV_NAME, NewVarNames.back().c_str());
					args.appendArray(SV_VARDATATYPE, ZoneDataTypes.data());
					args.appendArray(SV_VALUELOCATION, ZoneDataLocs.data());
					IsOk = TecUtilDataSetAddVarX(args.getRef());

					if (IsOk) {
						NewVarNums.back() = TecUtilDataSetGetNumVars();
						TecUtilStateChanged(StateChange_VarsAdded, (ArbParam_t)Set(NewVarNums.back()).getRef());
					}
				}
			}
		}
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
						ArgList args;
						args.appendString(SV_NAME, NewVarNames.back().c_str());
						args.appendArray(SV_VARDATATYPE, ZoneDataTypes.data());
						args.appendArray(SV_VALUELOCATION, ZoneDataLocs.data());
						IsOk = TecUtilDataSetAddVarX(args.getRef());

						if (IsOk) {
							NewVarNums.back() = TecUtilDataSetGetNumVars();
							TecUtilStateChanged(StateChange_VarsAdded, (ArbParam_t)Set(NewVarNums.back()).getRef());
						}
					}
				}
			}
		}

		// add vars for gradient bundle probabilities based on gradient bundle partition function and their individual contributions to the partition function
		if (IsOk && TeNumInVarList >= 0) {
			string TmpStr = "I: Gradient bundle partition function contribution (from IN Te)";
			if (std::find(NewVarNames.begin(), NewVarNames.end(), TmpStr) == NewVarNames.end()) {
				NewVarNames.push_back(TmpStr);
				NewVarNums.push_back(VarNumByName(NewVarNames.back()));
				if (NewVarNums.back() < 0) {
					ArgList args;
					args.appendString(SV_NAME, NewVarNames.back().c_str());
					args.appendArray(SV_VARDATATYPE, ZoneDataTypes.data());
					args.appendArray(SV_VALUELOCATION, ZoneDataLocs.data());
					IsOk = TecUtilDataSetAddVarX(args.getRef());

					if (IsOk) {
						NewVarNums.back() = TecUtilDataSetGetNumVars();
						TecUtilStateChanged(StateChange_VarsAdded, (ArbParam_t)Set(NewVarNums.back()).getRef());
					}
				}
			}

			TmpStr = "I: Gradient bundle probability (from IN Te)";
			if (std::find(NewVarNames.begin(), NewVarNames.end(), TmpStr) == NewVarNames.end()) {
				NewVarNames.push_back(TmpStr);
				NewVarNums.push_back(VarNumByName(NewVarNames.back()));
				if (NewVarNums.back() < 0) {
					ArgList args;
					args.appendString(SV_NAME, NewVarNames.back().c_str());
					args.appendArray(SV_VARDATATYPE, ZoneDataTypes.data());
					args.appendArray(SV_VALUELOCATION, ZoneDataLocs.data());
					IsOk = TecUtilDataSetAddVarX(args.getRef());

					if (IsOk) {
						NewVarNums.back() = TecUtilDataSetGetNumVars();
						TecUtilStateChanged(StateChange_VarsAdded, (ArbParam_t)Set(NewVarNums.back()).getRef());
					}
				}
			}
		}


		vector<string> AuxDataNames(IntVarNameList.size());
		for (int i = 0; i < IntVarNameList.size(); ++i)
			AuxDataNames[i] = RemoveStringChar(IntVarNameList[i], " ");
		AuxDataNames.emplace_back("Volume");

		vector<ValueLocation_e> DataLocs(TecUtilDataSetGetNumVars(), ValueLocation_Nodal);
		vector<FieldDataType_e> DataTypes(TecUtilDataSetGetNumVars(), FieldDataType_Float);
		vector<string> IntVarPrefixes = { "I: ","IN: ","INS: " };
		for (int i = 0; i < DataLocs.size(); ++i) {
			char * cstr;
			TecUtilVarGetName(i + 1, &cstr);
			string str = cstr;
			TecUtilStringDealloc(&cstr);
			for (auto const & s : IntVarPrefixes) {
				if (str.length() >= s.length()) {
					int j = str.compare(0, s.length(), s);
					if (j == 0) {
						DataLocs[i] = ValueLocation_CellCentered;
						break;
					}
				}
			}
		}

		Sphere = FESurface_c(SphereNodes, SphereElems);
		SphereZoneNum = Sphere.SaveAsTriFEZone(!NucleusName.empty() ? NucleusName + " Sphere" : "Sphere (" + CPString + ")", DataTypes, DataLocs, XYZVarNums);

		Set TmpSet(SphereZoneNum);
		TecUtilZoneSetMesh(SV_SHOW, TmpSet.getRef(), 0.0, FALSE);
		TecUtilZoneSetContour(SV_SHOW, TmpSet.getRef(), 0.0, TRUE);
		StyleValue styleValue;
		styleValue.set((Boolean_t)1, SV_FIELDLAYERS, SV_SHOWCONTOUR);

		if (IsOk)
		{
			ArgList_pa CurrentArgList = TecUtilArgListAlloc();
			TecUtilArgListAppendString(CurrentArgList, SV_P1, SV_FIELDMAP);
			TecUtilArgListAppendString(CurrentArgList, SV_P2, SV_EFFECTS);
			TecUtilArgListAppendString(CurrentArgList, SV_P3, SV_USETRANSLUCENCY);
			TecUtilArgListAppendSet(CurrentArgList, SV_OBJECTSET, TmpSet.getRef());
			TecUtilArgListAppendArbParam(CurrentArgList, SV_IVALUE, FALSE);
			TecUtilStyleSetLowLevelX(CurrentArgList);
			TecUtilArgListDealloc(&CurrentArgList);

			TecUtilStateChanged(StateChange_ZonesAdded, (ArbParam_t)TmpSet.getRef());
			TecUtilZoneSetScatter(SV_SHOW, TmpSet.getRef(), 0.0, FALSE);
			TecUtilZoneSetMesh(SV_SHOW, TmpSet.getRef(), 0.0, TRUE);
			TecUtilZoneSetActive(TmpSet.getRef(), AssignOp_PlusEquals);
		}

		/*
		*	Save node and triangle numbers to the fe volume zone's aux data
		*/

		if (IsOk)
		{
			// Get intersecting special gradient paths for identifying condensed basins
			vector<string> CheckPathTypes = {
				CSMAuxData.CC.ZoneSubTypeBondPathSegment,
				CSMAuxData.CC.ZoneSubTypeCageNuclearPath,
				CSMAuxData.CC.ZoneSubTypeRingNuclearPath
			};
			stringstream IntersectionNodeNums, IntersectionCPTypes, IntersectionCPNames, IntersectionCPNamesTotalCount;
			vec3 IntPoint;
			vector<int> XYZRhoVarNums = XYZVarNums;
			XYZRhoVarNums.push_back(RhoVarNum);
			for (int z = 1; z < TecUtilDataSetGetNumZones(); ++z) {
				if (AuxDataZoneItemMatches(z, CSMAuxData.CC.ZoneType, CSMAuxData.CC.ZoneTypeGP))
				{
					bool match = false;
					for (auto const & s : CheckPathTypes) {
						if (AuxDataZoneItemMatches(z, CSMAuxData.CC.ZoneSubType, s)) {
							match = true;
							break;
						}
					}
					if (match) {
						GradPath_c GP(z, XYZRhoVarNums, AddOnID);
						if (GP.GetSphereIntersectionPoint(CPPos, Radius, IntPoint) >= 0) {
							for (int i = 0; i < 2; ++i) {
								if (AuxDataZoneItemMatches(z, CSMAuxData.CC.GPEndTypes[i], CPNameList[CPTypeNum_Nuclear])) {
									int j = (i + 1) % 2;
									int NodeNum = Sphere.GetClosestNodeToPoint(IntPoint);
									int IntCPNum = stoi(AuxDataZoneGetItem(z, CSMAuxData.CC.GPEndNumStrs[j]));
									vector<int> CPTypeNumAndOffset = CPs.GetTypeNumOffsetFromTotOffset(IntCPNum-1);
									if (CPTypeNumAndOffset[0] < 0 || CPTypeNumAndOffset[1] < 0){
										TecUtilDialogErrMsg("Failed to get CP information for intersecting special gradient path");
									}
									else {
										IntersectionNodeNums << NodeNum << ",";
										IntersectionCPNamesTotalCount << CPNameList[CPTypeNumAndOffset[0]] << " " << IntCPNum << ",";
										IntersectionCPTypes << CPNameList[CPTypeNumAndOffset[0]] << ",";
										IntersectionCPNames << CPNameList[CPTypeNumAndOffset[0]] << " " << CPTypeNumAndOffset[1] + 1 << ",";
									}
								}
							}
						}
					}
				}
			}


			AuxDataZoneSetItem(SphereZoneNum, CSMAuxData.GBA.SphereCPName, CPString);
			AuxDataZoneSetItem(SphereZoneNum, CSMAuxData.GBA.ZoneType, CSMAuxData.GBA.ZoneTypeSphereZone);
			AuxDataZoneSetItem(SphereZoneNum, CSMAuxData.GBA.SphereCPNum, to_string(CPNum));
			AuxDataZoneSetItem(SphereZoneNum, CSMAuxData.GBA.SphereConstrainedNodeNums, IntersectionNodeNums.str());
			AuxDataZoneSetItem(SphereZoneNum, CSMAuxData.GBA.SphereConstrainedNodeIntersectCPTypes, IntersectionCPTypes.str());
			AuxDataZoneSetItem(SphereZoneNum, CSMAuxData.GBA.SphereConstrainedNodeIntersectCPNames, IntersectionCPNames.str());
			AuxDataZoneSetItem(SphereZoneNum, CSMAuxData.GBA.SphereConstrainedNodeIntersectCPTotalOffsetNames, IntersectionCPNamesTotalCount.str());
			AuxDataZoneSetItem(SphereZoneNum, CSMAuxData.GBA.NumGBs, to_string(NumElems));
			AuxDataZoneSetItem(SphereZoneNum, CSMAuxData.GBA.PointsPerGP, to_string(NumGPPoints));
			AuxDataZoneSetItem(SphereZoneNum, CSMAuxData.GBA.SourceZoneNum, to_string(CPNuclearZoneNum));
			AuxDataZoneSetItem(SphereZoneNum, CSMAuxData.GBA.SourceNucleusName, NucleusName);
			AuxDataZoneSetItem(SphereZoneNum, CSMAuxData.GBA.SphereRadius, to_string(Radius));
			AuxDataZoneSetItem(SphereZoneNum, CSMAuxData.GBA.SphereRadiusCPDistRatio, to_string(UserRadius));
			AuxDataZoneSetItem(SphereZoneNum, CSMAuxData.GBA.RadialSphereApproximation, RadialSphereApprx ? "Yes" : "No");
			AuxDataZoneSetItem(SphereZoneNum, CSMAuxData.GBA.SphereElemSymbol, ElemName);
			AuxDataZoneSetItem(SphereZoneNum, CSMAuxData.GBA.SphereOrigin, string(to_string(CPPos[0]) + "," + to_string(CPPos[1]) + "," + to_string(CPPos[2])));
		}

		Sphere = FESurface_c(SphereZoneNum, VolZoneNum, XYZVarNums, IntVarNumList, true);
// 		Sphere.DoIntegrationNew(IntResolution - 1, TRUE);
		

		vec SphereTriangleAreas(Sphere.TriSphereElemSolidAngles());
		double TotalArea = sum(SphereTriangleAreas);
		vec SphereTriangleAreaFactors = SphereTriangleAreas / TotalArea;

		/*
		 * Add derived properties if applicable
		 */

		TotalDensity = 0.0;
		double TotalVolume = 0.0;
		for (auto const & elemIntVals : IntVals) {
			TotalDensity += elemIntVals[DensityVarNumInIntList];
			TotalVolume += elemIntVals[VolumeNumInVarList];
		}

		if (TeNumInVarList >= 0) {
			for (int i = 0; i < NumElems; ++i) {
				IntVals[i].push_back(IntVals[i][TeNumInVarList] / MAX(IntVals[i][DensityVarNumInIntList], 1e-100)); // Te per electron
			}
		}

		for (int i = 0; i < NumElems; ++i) {
			IntVals[i].push_back(SphereTriangleAreaFactors[i]); // Solid angle (area on unit sphere)
			IntVals[i].push_back(IntVals[i][VolumeNumInVarList] / MAX(TotalVolume, 1e-100));  // Volume fraction
		}

		int VolumeFractionNumInVarList = IntVals[0].size() - 1;

		/*
		 * Now add the condensed valence density, bader charge, and valence bader charge
		 */
		if (AtomicNumber > 0) {
			AuxDataZoneSetItem(SphereZoneNum, CSMAuxData.GBA.SphereAtomicNumber, to_string(AtomicNumber));

			double CoreElectronCount = ElementCoreElectronCount(AtomicNumber),
				ValenceElectronCount = TotalDensity - CoreElectronCount;

			AuxDataZoneSetItem(SphereZoneNum, CSMAuxData.GBA.SphereAtomicCoreECount, to_string(CoreElectronCount));

			for (int i = 0; i < NumElems; ++i) {
				IntVals[i].push_back(((double)AtomicNumber * SphereTriangleAreaFactors[i]) - IntVals[i][DensityVarNumInIntList]); // Area Bader charge
				IntVals[i].push_back(IntVals[i][DensityVarNumInIntList] - (CoreElectronCount * SphereTriangleAreaFactors[i])); // Valence charge
				IntVals[i].push_back((ValenceElectronCount * SphereTriangleAreaFactors[i]) - IntVals[i].back()); // Valence Bader charge
				if (TeNumInVarList >= 0) {
					IntVals[i].push_back(IntVals[i][TeNumInVarList] / MAX(IntVals[i][IntVals[i].size() - 2], 1e-100)); // Te per valence electron
				}
				IntVals[i].push_back((TotalDensity * SphereTriangleAreaFactors[i]) - IntVals[i][DensityVarNumInIntList]); // Area condensed charge

				IntVals[i].push_back(((double)AtomicNumber * IntVals[i][VolumeFractionNumInVarList]) - IntVals[i][DensityVarNumInIntList]); // Volume Bader charge
// 				IntVals[i].push_back(IntVals[i][DensityVarNumInIntList] - (CoreElectronCount * IntVals[i][VolumeFractionNumInVarList])); // Volume Valence
// 				IntVals[i].push_back((ValenceElectronCount * IntVals[i][VolumeFractionNumInVarList]) - IntVals[i].back()); // Volume Valence Bader charge
// 				if (TeNumInVarList >= 0) {
// 					IntVals[i].push_back(IntVals[i][TeNumInVarList] / MAX(IntVals[i][IntVals.size() - 2], 1e-100)); // Te per Volume valence electron
// 				}
				IntVals[i].push_back((TotalDensity * IntVals[i][VolumeFractionNumInVarList]) - IntVals[i][DensityVarNumInIntList]); // Volume condensed charge
			}
		}
// 		double OneOverTotalDensity = 1.0 / TotalDensity;
// 		for (int i = 0; i < NumElems; ++i) {
// 			double val, prob = IntVals[i][DensityVarNumInIntList] * OneOverTotalDensity;
// 			IntVals[i].push_back(prob); // P over N
// 			val = prob / SphereTriangleAreaFactors[i]; // P/(N alpha)
// 			IntVals[i].push_back(val);
// 			IntVals[i].push_back(-prob * log(prob)); // S[P] shannon condensed entropy
// // 			double val = prob * SphereTriangleAreaFactors[i] * (double)NumElems; // S[N alpha P]
// // 			IntVals[i].push_back(-val * log(val)); 
// // 			val = prob * (double)NumElems; // S[N P]
// // 			IntVals[i].push_back(-val * log(val)); 
// // 			val = prob * SphereTriangleAreaFactors[i]; // S[alpha P]
// // 			IntVals[i].push_back(-val * log(val)); 
// // 			val = prob / (SphereTriangleAreaFactors[i] * (double)NumElems); // S[P / (N alpha)]
// // 			IntVals[i].push_back(-val * log(val)); 
// // 			val = prob / (double)NumElems; // S[P / N]
// // 			IntVals[i].push_back(-val * log(val)); 
// 			val = prob / SphereTriangleAreaFactors[i]; // S[P / alpha]
// 			IntVals[i].push_back(-val * log(val));
// 			val = SphereTriangleAreaFactors[i]; // S[alpha]
// 			IntVals[i].push_back(-val * log(val));
// 			IntVals[i].push_back((-prob * log(prob)) / (-val * log(val)));
// 			IntVals[i].push_back((-val * log(val)) - (-prob * log(prob)));
// 		}

		/*
		 *	Get the normalized and normalized+scaled values
		 */
		vector<vector<double> > NormalizedValues = IntVals;
		vector<double> TotalList(IntVals[0].size(), 0.0),
			TotalNormlizedList(IntVals[0].size(), 0.0),
			IntScaleFactors(IntVals[0].size());

		for (int i = 0; i < NumElems; ++i) {
			for (int j = 0; j < IntVals[i].size(); ++j) {
// 				NormalizedValues[i][j] /= SphereTriangleAreas[i];
				NormalizedValues[i][j] /= SphereTriangleAreaFactors[i];
				TotalNormlizedList[j] += NormalizedValues[i][j];
				TotalList[j] += IntVals[i][j];
			}
		}

		for (int i = 0; i < IntVals[0].size(); ++i) {
			IntScaleFactors[i] = TotalList[i] / TotalNormlizedList[i];
		}

		for (int i = 0; i < NumElems; ++i) {
			auto TmpVec = IntVals[i];
			IntVals[i] = vector<double>();
			IntVals[i].reserve(3 * TmpVec.size());
			for (int j = 0; j < TmpVec.size(); ++j) {
				IntVals[i].push_back(TmpVec[j]);
				IntVals[i].push_back(NormalizedValues[i][j]);
				IntVals[i].push_back(NormalizedValues[i][j] * IntScaleFactors[j]);
			}
		}

		// Calculate Te gradient bundle partition function, each gradient bundle's contribution to it, and each gradient bundle's resulting probability.
		double TePartitionFunction = 0.0;
		double TePartitionFunctionBeta = 1.0;
		double TotalProb = 0.0;
		if (TeNumInVarList >= 0) {
			int TeNumInINSVarList = (TeNumInVarList) * 3 + 1; // points to the area-normalized "IN: Te" values.
			// First partition function terms and total partition function
			for (int i = 0; i < NumElems; ++i) {
				IntVals[i].push_back(exp(-TePartitionFunctionBeta * IntVals[i][TeNumInINSVarList]));
				TePartitionFunction += IntVals[i].back();
			}
			// now probabilities
			for (int i = 0; i < NumElems; ++i){
				IntVals[i].push_back(IntVals[i].back() / TePartitionFunction);
				TotalProb += IntVals[i].back();
			}
		}

		IntVarNameList.emplace_back("Volume");
		TotalList.push_back(TotalVolume);
		IntVarNameList.emplace_back("Area");
		TotalList.push_back(TotalArea);
		if (TeNumInVarList >= 0) {
			IntVarNameList.emplace_back("TePartitionFunction");
			TotalList.push_back(TePartitionFunction);

			IntVarNameList.emplace_back("TeGBProbability");
			TotalList.push_back(TotalProb);
		}
		
		AuxDataZoneSetItem(SphereZoneNum, CSMAuxData.GBA.AtomicBasinIntegrationVariables, VectorToString(IntVarNameList, ","));
		AuxDataZoneSetItem(SphereZoneNum, CSMAuxData.GBA.AtomicBasinIntegrationValues, VectorToString(TotalList, ","));

		AuxDataZoneSetItem(SphereZoneNum, CSMAuxData.GBA.GPRhoCutoff, to_string(CutoffVal));
		AuxDataZoneSetItem(SphereZoneNum, CSMAuxData.GBA.SphereMinNumGPs, to_string(Level));
		AuxDataZoneSetItem(SphereZoneNum, CSMAuxData.GBA.SphereSubdivisionLevel, to_string(MaxSubdivisions));
		AuxDataZoneSetItem(SphereZoneNum, CSMAuxData.GBA.SphereSubdivisionTightness, to_string(SubdivisionTightness));
		AuxDataZoneSetItem(SphereZoneNum, CSMAuxData.GBA.SphereNumBondPathCoincidentGBs, to_string(NumAngularGBsAroundBPIntersectionNodes));
		AuxDataZoneSetItem(SphereZoneNum, CSMAuxData.GBA.GBSurfaceGPMaxSpacing, to_string(EdgeGPCheckDistance));

		if (RadialSphereApprx) {
			AuxDataZoneSetItem(SphereZoneNum, CSMAuxData.GBA.IntVarSphereInteriorVals, VectorToString(SphereTotalIntVals, ","));
// 			vector<double> TotalIntTotals = SphereTotalIntVals; // Saving this since sphere integration isn't represented in the GB intVals
// 			for (int i = 0; i < SphereTotalIntVals.size(); ++i){
// 				TotalIntTotals[i] += TotalList[i];
// 			}
// 			AuxDataZoneSetItem(SphereZoneNum, CSMAuxData.GBA.AtomicBasinIntegrationTotalValues, VectorToString(TotalIntTotals, ","));

			// some stats for the sphere interior region to gauge 1) whether it was in the radial region and 2) to gauge grid error in the nuclear region
			vector<double> SphereMins(SphereTotalIntVals.size()),
				SphereMaxes(SphereTotalIntVals.size()),
				SphereMeans(SphereTotalIntVals.size()),
				SphereStdDevs(SphereTotalIntVals.size()),
				SphereRanges(SphereTotalIntVals.size());
			
			for (int iVar = 0; iVar < SphereTotalIntVals.size(); ++iVar){
				vector<double> tmpVector(SphereIntVals.size());
				for (int ti = 0; ti < SphereIntVals.size(); ++ti){
					tmpVector[ti] = SphereIntVals[ti][iVar];
				}
				vec tmpVec(tmpVector);
				SphereMins[iVar] = min(tmpVec);
				SphereMaxes[iVar] = max(tmpVec);
				SphereMeans[iVar] = mean(tmpVec);
				SphereStdDevs[iVar] = stddev(tmpVec);
				SphereRanges[iVar] = SphereMaxes[iVar] - SphereMins[iVar];
			}
			AuxDataZoneSetItem(SphereZoneNum, CSMAuxData.GBA.IntVarSphereInteriorMaxes, VectorToString(SphereMaxes, ","));
			AuxDataZoneSetItem(SphereZoneNum, CSMAuxData.GBA.IntVarSphereInteriorMins, VectorToString(SphereMins, ","));
			AuxDataZoneSetItem(SphereZoneNum, CSMAuxData.GBA.IntVarSphereInteriorMeans, VectorToString(SphereMeans, ","));
			AuxDataZoneSetItem(SphereZoneNum, CSMAuxData.GBA.IntVarSphereInteriorStdDevs, VectorToString(SphereStdDevs, ","));
			AuxDataZoneSetItem(SphereZoneNum, CSMAuxData.GBA.IntVarSphereInteriorRanges, VectorToString(SphereRanges, ","));
		}



		vector<FieldDataPointer_c> Ptrs(IntVals[0].size());
		for (int j = 0; j < Ptrs.size(); ++j) {
			Ptrs[j].InitializeWritePtr(SphereZoneNum, NewVarNums[j]);
		}

		// Write cell-centered integration values to sphere and FE zones
		for (int e = 0; e < NumElems; ++e) {
			for (int j = 0; j < IntVals[e].size(); ++j) {
				Ptrs[j].Write(e, IntVals[e][j]);
			}
		}
		for (auto & i : Ptrs)
			i.Close();

		vector<bool> SaveGP;

		int NumColors = 63,
			ColorShift = 1;
		ColorIndex_t ColorIter = 0;

		if (SaveGBs) {
			Time1 = high_resolution_clock::now();
			TmpString = ProgressStr1.str() + "Saving GBs";
			if (SaveGPs)
				TmpString += " & GPs";
#ifdef SUPERDEBUG
			for (int ti : ElemTodo) {
				{
#else
			for (int ti = 0; ti < GPPtrList.size(); ++ti) {
				if (StatusUpdate(ti, GPPtrList.size() * (2 ? RadialSphereApprx : 1), TmpString, AddOnID, Time1)) {
#endif
					FESurface_c GradientBundle;
					GradientBundle.MakeFromGPs(GPPtrList[ti], true, true);
					GradientBundle.SaveAsTriFEZone(CPString + ": Gradient Bundle " + to_string(ti + 1) + " " + NETypeStrs[ElemTypes[ti]] + "-Type (Zone " + to_string(TecUtilDataSetGetNumZones() + 1) + ")", vector<FieldDataType_e>(TecUtilDataSetGetNumVars(), FieldDataType_Float), vector<ValueLocation_e>(GradientBundle.GetXYZListPtr()->size(), ValueLocation_Nodal), XYZVarNums, RhoVarNum);
					
					ColorIndex_t GBColor = -1;
					for (auto g : GPPtrList[ti]) {
						if (g->GetStartEndCPNum(1) >= 0) {
							GBColor = ColorIndex_t((g->GetStartEndCPNum(1) % NumColors) + 1);
							break;
						}
					}
					if (GBColor == -1) {
						GBColor = ColorIndex_t((++ColorIter % NumColors) + 1);
					}

					if (GradientBundle.IsMade() && GradientBundle.GetZoneNum() > 0) {
						TecUtilZoneSetActive(Set(GradientBundle.GetZoneNum()).getRef(), AssignOp_MinusEquals);

						TecUtilZoneSetShade(SV_SHOW, Set(GradientBundle.GetZoneNum()).getRef(), 0.0, TRUE);
						TecUtilZoneSetShade(SV_COLOR, Set(GradientBundle.GetZoneNum()).getRef(), 0.0, GBColor);

						AuxDataZoneSetItem(GradientBundle.GetZoneNum(), CSMAuxData.GBA.SourceZoneNum, to_string(SphereZoneNum));
						AuxDataZoneSetItem(GradientBundle.GetZoneNum(), CSMAuxData.GBA.SphereCPName, CPString);
						AuxDataZoneSetItem(GradientBundle.GetZoneNum(), CSMAuxData.GBA.ZoneType, CSMAuxData.GBA.ZoneTypeDGB);
						AuxDataZoneSetItem(GradientBundle.GetZoneNum(), CSMAuxData.GBA.ElemNum, to_string(ti + 1));
						AuxDataZoneSetItem(GradientBundle.GetZoneNum(), CSMAuxData.GBA.SourceNucleusName, NucleusName);

						for (int i = 0; i < 3 && IsOk; ++i)
							AuxDataZoneSetItem(GradientBundle.GetZoneNum(), CSMAuxData.GBA.NodeNums[i], to_string(SphereElems[ti][i] + 1));


						// Save integration results to gradient bundle
						for (int vi = 0; vi < NewVarNums.size(); ++vi) {
							FieldData_pa FDPtr = TecUtilDataValueGetWritableNativeRef(GradientBundle.GetZoneNum(), NewVarNums[vi]);
							if (VALID_REF(FDPtr)) {

								
								if (TecUtilDataValueGetType(GradientBundle.GetZoneNum(), NewVarNums[vi]) == FieldDataType_Double) {
									vector<double> tmpDblVec(GradientBundle.GetXYZListPtr()->size(), IntVals[ti][vi]);
									TecUtilDataValueArraySetByRef(FDPtr, 1, tmpDblVec.size(), tmpDblVec.data());
								}
								else {
									vector<float> tmpDblVec(GradientBundle.GetXYZListPtr()->size(), IntVals[ti][vi]);
									TecUtilDataValueArraySetByRef(FDPtr, 1, tmpDblVec.size(), tmpDblVec.data());
								}
							}
							else {
								TecUtilDialogErrMsg("Failed to get write pointer for gradient bundle");
							}

// 
// 							if (ti >= 5) {
// 								TecUtilDialogErrMsg(string(to_string(IntVals[ti-5][vi]) + " " + to_string(TecUtilDataValueGetByZoneVar(GradientBundles[ti-5].GetZoneNum(), NewVarNums[vi], 1))).c_str());
// 							}

// 							for (int i = 0; i < GradientBundle.GetXYZListPtr()->size(); ++i){
// 								TecUtilDataValueSetByZoneVar(GradientBundle.GetZoneNum(), NewVarNums[vi], i + 1, IntVals[ti][vi]);
// 							}

// 							FieldDataPointer_c tmpPtr;
// 							tmpPtr.InitializeWritePtr(GradientBundle.GetZoneNum(), NewVarNums[vi]);
// 							for (int i = 0; i < tmpPtr.Size(); ++i){
// 								tmpPtr.Write(i, IntVals[ti][vi]);
// 							}
						}
					}
					else {
						TecUtilDialogErrMsg(string("Failed to save GB " + to_string(ti + 1)).c_str());
					}

					if (SaveGPs) {
// 						for (int i = 0; i < SphereElemGPInds[ti].size(); ++i) {
// 							int ni = SphereElemGPInds[ti][i];
// 							if (GradPaths[ni].IsMade()) {
// 								GradPaths[ni].SaveAsOrderedZone("GP " + CPName + ": Element " + to_string(ti + 1) + ", Node " + to_string(ni + 1));
// 								TecUtilZoneSetActive(Set(GradPaths[ni].GetZoneNum()).getRef(), AssignOp_MinusEquals);
// 								AuxDataZoneSetItem(GradPaths[ni].GetZoneNum(), CSMAuxData.GBA.SourceZoneNum, to_string(SphereZoneNum));
// 								AuxDataZoneSetItem(GradPaths[ni].GetZoneNum(), CSMAuxData.GBA.SphereCPName, CPString);
// 								AuxDataZoneSetItem(GradPaths[ni].GetZoneNum(), CSMAuxData.GBA.SphereCPNum, to_string(CPNum));
// 								AuxDataZoneSetItem(GradPaths[ni].GetZoneNum(), CSMAuxData.GBA.GPNodeNum, to_string(ni + 1));
// 								AuxDataZoneSetItem(GradPaths[ni].GetZoneNum(), CSMAuxData.GBA.ZoneType, CSMAuxData.GBA.ZoneTypeGradPath);
// 								AuxDataZoneSetItem(GradPaths[ni].GetZoneNum(), CSMAuxData.GBA.SourceNucleusName, NucleusName);
// // 								AuxDataZoneSetItem(GradPaths[ni].GetZoneNum(), "Debug.NodeTypeSurfaceNumber", to_string(NodeTypes[SphereElems[ti][i]]) + ", " + to_string(NodeIntSurfNums[SphereElems[ti][ni]]));
// 							}
// 						}
						for (int i = 0; i < GPPtrList[ti].size(); ++i) {
							auto gp = *GPPtrList[ti][i];
							if (gp.IsMade()) {
								gp.SaveAsOrderedZone(CPString + ": GP " + CPName + ": Element " + to_string(ti + 1) + "." + to_string(i + 1) + " with GB");
								if (gp.GetZoneNum() > 0) {
									auto GPColor = ColorIndex_t((i % NumColors) + 1);
									TecUtilZoneSetActive(Set(gp.GetZoneNum()).getRef(), AssignOp_MinusEquals);
#ifdef SUPERDEBUG
									TecUtilZoneSetMesh(SV_COLOR, Set(gp.GetZoneNum()).getRef(), 0.0, GPColor);
									TecUtilZoneSetMesh(SV_SHOW, Set(gp.GetZoneNum()).getRef(), 0.0, TRUE);
									TecUtilZoneSetScatter(SV_COLOR, Set(gp.GetZoneNum()).getRef(), 0.0, GPColor);
									TecUtilZoneSetScatter(SV_FRAMESIZE, Set(gp.GetZoneNum()).getRef(), 0.4, 0);
									TecUtilZoneSetScatter(SV_FILLMODE, Set(gp.GetZoneNum()).getRef(), 0.0, FillMode_UseLineColor);
									TecUtilZoneSetScatterSymbolShape(SV_GEOMSHAPE, Set(gp.GetZoneNum()).getRef(), GeomShape_Circle);
#else
									TecUtilZoneSetMesh(SV_COLOR, Set(gp.GetZoneNum()).getRef(), 0.0, GBColor);
#endif
									AuxDataZoneSetItem(gp.GetZoneNum(), CSMAuxData.GBA.SourceZoneNum, to_string(SphereZoneNum));
									AuxDataZoneSetItem(gp.GetZoneNum(), CSMAuxData.GBA.SphereCPName, CPString);
									AuxDataZoneSetItem(gp.GetZoneNum(), CSMAuxData.GBA.SphereCPNum, to_string(CPNum));
									AuxDataZoneSetItem(gp.GetZoneNum(), CSMAuxData.GBA.ZoneType, CSMAuxData.GBA.ZoneTypeGradPath);
									AuxDataZoneSetItem(gp.GetZoneNum(), CSMAuxData.GBA.SourceNucleusName, NucleusName);
									// 								AuxDataZoneSetItem(GradPaths[ni].GetZoneNum(), "Debug.NodeTypeSurfaceNumber", to_string(NodeTypes[SphereElems[ti][i]]) + ", " + to_string(NodeIntSurfNums[SphereElems[ti][ni]]));
								}
								else{
									TecUtilDialogErrMsg(string("Failed to save GP " + to_string(i + 1) + " of GB " + to_string(ti + 1)).c_str());
								}
							}
						}
					}
				}
			}
		}
		else if (SaveGPs) {
			SaveGP.resize(GradPaths.size(), false);
			int NumToSave = 0;

#ifdef SUPERDEBUG
			for (int ti : ElemTodo)
#else
			for (int ti = 0; ti < GPPtrList.size(); ++ti)
#endif
				for (int ni : SphereElemGPInds[ti]) {
					NumToSave++;
					SaveGP[ni] = true;
				}

			TmpString = ProgressStr1.str() + "Saving GPs";
			Time1 = high_resolution_clock::now();

			int NumComplete = 0;
			for (int ni = 0; ni < GradPaths.size(); ++ni) {
				if (StatusUpdate(NumComplete, NumToSave + InnerGradPaths.size(), TmpString, AddOnID, Time1) && SaveGP[ni] && GradPaths[ni].IsMade()) {
					ColorIndex_t GPColor = -1;
					if (GradPaths[ni].GetStartEndCPNum(1) >= 0){
						GPColor = ColorIndex_t((GradPaths[ni].GetStartEndCPNum(1) % NumColors) + 1);
					}
					else{
						GPColor = ColorIndex_t((++ColorIter % NumColors) + 1);
					}

					NumComplete++;
					string TypeStr;
					if (ni < NodeTypes.size())
						TypeStr = " " + NETypeStrs[NodeTypes[ni]] + "-Type";
					GradPaths[ni].SaveAsOrderedZone(CPString + ": GP " + CPName + ": Node " + to_string(ni + 1) + TypeStr);
					if (GradPaths[ni].GetZoneNum() > 0) {
						TecUtilZoneSetActive(Set(GradPaths[ni].GetZoneNum()).getRef(), AssignOp_MinusEquals);
						TecUtilZoneSetMesh(SV_COLOR, Set(GradPaths[ni].GetZoneNum()).getRef(), 0.0, GPColor);
						AuxDataZoneSetItem(GradPaths[ni].GetZoneNum(), CSMAuxData.GBA.SourceZoneNum, to_string(Sphere.GetZoneNum()));
						AuxDataZoneSetItem(GradPaths[ni].GetZoneNum(), CSMAuxData.GBA.SphereCPName, CPString);
						AuxDataZoneSetItem(GradPaths[ni].GetZoneNum(), CSMAuxData.GBA.SphereCPNum, to_string(CPNum));
						AuxDataZoneSetItem(GradPaths[ni].GetZoneNum(), CSMAuxData.GBA.GPNodeNum, to_string(ni + 1));
						AuxDataZoneSetItem(GradPaths[ni].GetZoneNum(), CSMAuxData.GBA.ZoneType, CSMAuxData.GBA.ZoneTypeGradPath);
						AuxDataZoneSetItem(GradPaths[ni].GetZoneNum(), CSMAuxData.GBA.SourceNucleusName, NucleusName);
					}
					else {
						TecUtilDialogErrMsg(string("Failed to save GP " + to_string(ni + 1)).c_str());
					}
				}
			}
		}

		if (RadialSphereApprx && false) {
			if (SaveGBs) {
#ifdef SUPERDEBUG
				for (int ti : ElemTodo) {
					{
#else
				for (int ti = 0; ti < GPPtrList.size(); ++ti) {
					if (StatusUpdate(ti + GPPtrList.size(), GPPtrList.size() * (2 ? RadialSphereApprx : 1), TmpString, AddOnID, Time1)) {
#endif
						vector<const GradPath_c *> InnerGPPtrs(3);
						for (int tj = 0; tj < 3; ++tj) {
							InnerGPPtrs[tj] = &InnerGradPaths[SphereElems[ti][tj]];
						}
						FESurface_c GradientBundle;
						GradientBundle.MakeFromGPs(InnerGPPtrs, true, true);
						GradientBundle.SaveAsTriFEZone(CPString + ": Sphere Gradient Bundle " + to_string(ti + 1) + " " + NETypeStrs[ElemTypes[ti]] + "-Type (Zone " + to_string(TecUtilDataSetGetNumZones() + 1) + ")", vector<FieldDataType_e>(TecUtilDataSetGetNumVars(), FieldDataType_Float), vector<ValueLocation_e>(GradientBundle.GetXYZListPtr()->size(), ValueLocation_Nodal), XYZVarNums, RhoVarNum);

						ColorIndex_t GBColor = -1;
						for (auto g : GPPtrList[ti]) {
							if (g->GetStartEndCPNum(1) >= 0) {
								GBColor = ColorIndex_t((g->GetStartEndCPNum(1) % NumColors) + 1);
								break;
							}
						}
						if (GBColor == -1) {
							GBColor = ColorIndex_t((++ColorIter % NumColors) + 1);
						}

						if (GradientBundle.IsMade() && GradientBundle.GetZoneNum() > 0) {
							TecUtilZoneSetActive(Set(GradientBundle.GetZoneNum()).getRef(), AssignOp_MinusEquals);

							TecUtilZoneSetShade(SV_SHOW, Set(GradientBundle.GetZoneNum()).getRef(), 0.0, TRUE);
							TecUtilZoneSetShade(SV_COLOR, Set(GradientBundle.GetZoneNum()).getRef(), 0.0, GBColor);

							AuxDataZoneSetItem(GradientBundle.GetZoneNum(), CSMAuxData.GBA.SourceZoneNum, to_string(SphereZoneNum));
							AuxDataZoneSetItem(GradientBundle.GetZoneNum(), CSMAuxData.GBA.SphereCPName, CPString);
							AuxDataZoneSetItem(GradientBundle.GetZoneNum(), CSMAuxData.GBA.ZoneType, CSMAuxData.GBA.ZoneTypeDGB);
							AuxDataZoneSetItem(GradientBundle.GetZoneNum(), CSMAuxData.GBA.ElemNum, to_string(ti + 1));
							AuxDataZoneSetItem(GradientBundle.GetZoneNum(), CSMAuxData.GBA.SourceNucleusName, NucleusName);

							for (int i = 0; i < 3 && IsOk; ++i)
								AuxDataZoneSetItem(GradientBundle.GetZoneNum(), CSMAuxData.GBA.NodeNums[i], to_string(SphereElems[ti][i] + 1));


							// Save integration results to gradient bundle
							for (int vi = 0; vi < NewVarNums.size(); ++vi) {
								FieldData_pa FDPtr = TecUtilDataValueGetWritableNativeRef(GradientBundle.GetZoneNum(), NewVarNums[vi]);
								if (VALID_REF(FDPtr)) {
									if (TecUtilDataValueGetType(GradientBundle.GetZoneNum(), NewVarNums[vi]) == FieldDataType_Double) {
										vector<double> tmpDblVec(GradientBundle.GetXYZListPtr()->size(), IntVals[ti][vi]);
										TecUtilDataValueArraySetByRef(FDPtr, 1, tmpDblVec.size(), tmpDblVec.data());
									}
									else {
										vector<float> tmpDblVec(GradientBundle.GetXYZListPtr()->size(), IntVals[ti][vi]);
										TecUtilDataValueArraySetByRef(FDPtr, 1, tmpDblVec.size(), tmpDblVec.data());
									}
								}
								else {
									TecUtilDialogErrMsg("Failed to get write pointer for gradient bundle");
								}

								// 
								// 							if (ti >= 5) {
								// 								TecUtilDialogErrMsg(string(to_string(IntVals[ti-5][vi]) + " " + to_string(TecUtilDataValueGetByZoneVar(GradientBundles[ti-5].GetZoneNum(), NewVarNums[vi], 1))).c_str());
								// 							}

								// 							for (int i = 0; i < GradientBundle.GetXYZListPtr()->size(); ++i){
								// 								TecUtilDataValueSetByZoneVar(GradientBundle.GetZoneNum(), NewVarNums[vi], i + 1, IntVals[ti][vi]);
								// 							}

								// 							FieldDataPointer_c tmpPtr;
								// 							tmpPtr.InitializeWritePtr(GradientBundle.GetZoneNum(), NewVarNums[vi]);
								// 							for (int i = 0; i < tmpPtr.Size(); ++i){
								// 								tmpPtr.Write(i, IntVals[ti][vi]);
								// 							}
							}
						}
						else {
							TecUtilDialogErrMsg(string("Failed to save GB " + to_string(ti + 1)).c_str());
						}

						if (SaveGPs) {
							// 						for (int i = 0; i < SphereElemGPInds[ti].size(); ++i) {
							// 							int ni = SphereElemGPInds[ti][i];
							// 							if (GradPaths[ni].IsMade()) {
							// 								GradPaths[ni].SaveAsOrderedZone("GP " + CPName + ": Element " + to_string(ti + 1) + ", Node " + to_string(ni + 1));
							// 								TecUtilZoneSetActive(Set(GradPaths[ni].GetZoneNum()).getRef(), AssignOp_MinusEquals);
							// 								AuxDataZoneSetItem(GradPaths[ni].GetZoneNum(), CSMAuxData.GBA.SourceZoneNum, to_string(SphereZoneNum));
							// 								AuxDataZoneSetItem(GradPaths[ni].GetZoneNum(), CSMAuxData.GBA.SphereCPName, CPString);
							// 								AuxDataZoneSetItem(GradPaths[ni].GetZoneNum(), CSMAuxData.GBA.SphereCPNum, to_string(CPNum));
							// 								AuxDataZoneSetItem(GradPaths[ni].GetZoneNum(), CSMAuxData.GBA.GPNodeNum, to_string(ni + 1));
							// 								AuxDataZoneSetItem(GradPaths[ni].GetZoneNum(), CSMAuxData.GBA.ZoneType, CSMAuxData.GBA.ZoneTypeGradPath);
							// 								AuxDataZoneSetItem(GradPaths[ni].GetZoneNum(), CSMAuxData.GBA.SourceNucleusName, NucleusName);
							// // 								AuxDataZoneSetItem(GradPaths[ni].GetZoneNum(), "Debug.NodeTypeSurfaceNumber", to_string(NodeTypes[SphereElems[ti][i]]) + ", " + to_string(NodeIntSurfNums[SphereElems[ti][ni]]));
							// 							}
							// 						}
							for (int i = 0; i < 3; ++i) {
								auto gp = InnerGradPaths[SphereElems[ti][i]];
								if (gp.IsMade()) {
									gp.SaveAsOrderedZone(CPString + ": Sphere GP " + CPName + ": Element " + to_string(ti + 1) + "." + to_string(i + 1) + " with GB");
									if (gp.GetZoneNum() > 0) {
										TecUtilZoneSetActive(Set(gp.GetZoneNum()).getRef(), AssignOp_MinusEquals);
										TecUtilZoneSetMesh(SV_COLOR, Set(gp.GetZoneNum()).getRef(), 0.0, GBColor);
										AuxDataZoneSetItem(gp.GetZoneNum(), CSMAuxData.GBA.SourceZoneNum, to_string(SphereZoneNum));
										AuxDataZoneSetItem(gp.GetZoneNum(), CSMAuxData.GBA.SphereCPName, CPString);
										AuxDataZoneSetItem(gp.GetZoneNum(), CSMAuxData.GBA.SphereCPNum, to_string(CPNum));
										AuxDataZoneSetItem(gp.GetZoneNum(), CSMAuxData.GBA.ZoneType, CSMAuxData.GBA.ZoneTypeGradPath);
										AuxDataZoneSetItem(gp.GetZoneNum(), CSMAuxData.GBA.SourceNucleusName, NucleusName);
										// 								AuxDataZoneSetItem(GradPaths[ni].GetZoneNum(), "Debug.NodeTypeSurfaceNumber", to_string(NodeTypes[SphereElems[ti][i]]) + ", " + to_string(NodeIntSurfNums[SphereElems[ti][ni]]));
									}
									else {
										TecUtilDialogErrMsg(string("Failed to save GP " + to_string(i + 1) + " of GB " + to_string(ti + 1)).c_str());
									}
								}
							}
						}
					}
				}
			}
			else if (SaveGPs) {
				SaveGP.resize(GradPaths.size(), false);
				int NumToSave = 0;

	#ifdef SUPERDEBUG
				for (int ti : ElemTodo)
	#else
				for (int ti = 0; ti < GPPtrList.size(); ++ti)
	#endif
					for (int ni : SphereElemGPInds[ti]) {
						NumToSave++;
						SaveGP[ni] = true;
					}

				TmpString = ProgressStr1.str() + "Saving GPs";
				Time1 = high_resolution_clock::now();

				int NumComplete = 0;
				for (int ni = 0; ni < GradPaths.size(); ++ni) {
					if (StatusUpdate(NumComplete + NumToSave, NumToSave + InnerGradPaths.size(), TmpString, AddOnID, Time1) && SaveGP[ni] && InnerGradPaths[ni].IsMade()) {
						ColorIndex_t GPColor = -1;
						if (InnerGradPaths[ni].GetStartEndCPNum(1) >= 0) {
							GPColor = ColorIndex_t((InnerGradPaths[ni].GetStartEndCPNum(1) % NumColors) + 1);
						}
						else {
							GPColor = ColorIndex_t((++ColorIter % NumColors) + 1);
						}

						NumComplete++;
						string TypeStr;
						if (ni < NodeTypes.size())
							TypeStr = " " + NETypeStrs[NodeTypes[ni]] + "-Type";
						InnerGradPaths[ni].SaveAsOrderedZone(CPString + ": Sphere GP " + CPName + ": Node " + to_string(ni + 1) + TypeStr);
						if (InnerGradPaths[ni].GetZoneNum() > 0) {
							TecUtilZoneSetActive(Set(InnerGradPaths[ni].GetZoneNum()).getRef(), AssignOp_MinusEquals);
							TecUtilZoneSetMesh(SV_COLOR, Set(InnerGradPaths[ni].GetZoneNum()).getRef(), 0.0, GPColor);
							AuxDataZoneSetItem(InnerGradPaths[ni].GetZoneNum(), CSMAuxData.GBA.SourceZoneNum, to_string(Sphere.GetZoneNum()));
							AuxDataZoneSetItem(InnerGradPaths[ni].GetZoneNum(), CSMAuxData.GBA.SphereCPName, CPString);
							AuxDataZoneSetItem(InnerGradPaths[ni].GetZoneNum(), CSMAuxData.GBA.SphereCPNum, to_string(CPNum));
							AuxDataZoneSetItem(InnerGradPaths[ni].GetZoneNum(), CSMAuxData.GBA.GPNodeNum, to_string(ni + 1));
							AuxDataZoneSetItem(InnerGradPaths[ni].GetZoneNum(), CSMAuxData.GBA.ZoneType, CSMAuxData.GBA.ZoneTypeGradPath);
							AuxDataZoneSetItem(InnerGradPaths[ni].GetZoneNum(), CSMAuxData.GBA.SourceNucleusName, NucleusName);
						}
						else {
							TecUtilDialogErrMsg(string("Failed to save GP " + to_string(ni + 1)).c_str());
						}
					}
				}
			}
		}

		/*
		 * Find and save topological cage wedges to zones
		 */
		{
			std::map<int, std::set<int> > CageNumToElemNums;
			auto Elems = Sphere.GetElemListPtr();
			
			vector<bool> ElemVisited(Elems->size(), false);
			Sphere.GenerateElemConnectivity(2);
			auto ElemConnectivity = Sphere.GetElemConnectivityListPtr();

			int FarFieldCageIter = 0;

			for (int ei = 0; ei < Elems->size(); ++ei) {
				if (!ElemVisited[ei]){
					ElemVisited[ei] = true;

					std::set<int> CageWedgeElemNums;
					CageWedgeElemNums.insert(ei);
					std::queue<int> q;
					q.push(ei);

					while (!q.empty()){
						// check neighbor edges to see if ring surface intersection
						auto eq = q.front();
						q.pop();
						for (auto const & ej : ElemConnectivity->at(eq)){
							if (!ElemVisited[ej]){
								ElemVisited[ej] = true;
								bool SharedEdgeFound = false;
								for (int ci = 0; ci < 3 && !SharedEdgeFound; ++ci){
									auto edge_i = MakeEdge(Elems->at(eq)[ci], Elems->at(eq)[(ci + 1) % 3]);
									for (int cj = 0; cj < 3 && !SharedEdgeFound; ++cj){
										auto edge_j = MakeEdge(Elems->at(ej)[cj], Elems->at(ej)[(cj + 1) % 3]);
										if (edge_i == edge_j) {
											SharedEdgeFound = true;
											if (!EdgeIntSurfNums.count(edge_i)) {
												// shared edge is not a ring surface intersection, so add neighbor
												q.push(ej);
												CageWedgeElemNums.insert(ej);
											}
											else{
												ElemVisited[ej] = false;
											}
										}
									}
								}
							}
						}
					}

					// check whether to keep set or not
					if (CageWedgeElemNums.size() < Elems->size()){
						// Fewer elements in cage wedge basin than in sphere, meaning there
						// must be another region completely separated by ring surface intersections on the sphere
						// Get cage CP number for wedge. by finding common cage CP to the ring surface intersections
					
						// get ring surface edges
						std::set<Edge> RingSurfEdges;
						for (auto const & ei : CageWedgeElemNums){
							for (int ci = 0; ci < 3; ++ci){
								auto edge = MakeEdge(Elems->at(ei)[ci], Elems->at(ei)[(ci + 1) % 3]);
								if (EdgeIntSurfNums.count(edge)){
									RingSurfEdges.insert(edge);
								}
							}
						}

						// get list of unique indices specifying ring-cage paths
						std::set<std::pair<int, int> > RingCagePathInds;
						for (auto const & e : RingSurfEdges){
							int SurfNum = EdgeIntSurfNums[e];
							RingCagePathInds.emplace(SurfNum, 0);
							RingCagePathInds.emplace(SurfNum, 1);
						}
						auto RingCagePathIndsCopy = RingCagePathInds;

						// Now loop over the RingSurfEdges and remove from RingCagePathInds any pair that isn't present in all RingSurfEdges
						for (auto const & e : RingSurfEdges){
							int SurfNum = EdgeIntSurfNums[e];
							for (auto const & ri : RingCagePathInds){
								if (RingCagePathIndsCopy.count(ri) && Distance(IntSurfRingCagePaths[SurfNum][0][-1], IntSurfRingCagePaths[ri.first][ri.second][-1]) > 0.05
									&& Distance(IntSurfRingCagePaths[SurfNum][1][-1], IntSurfRingCagePaths[ri.first][ri.second][-1]) > 0.05){
									RingCagePathIndsCopy.erase(ri);
								}
							}
							if (RingCagePathIndsCopy.size() == 1){
								break;
							}
						}

						// Get the cage CP
						if (!RingCagePathIndsCopy.empty()) {
							auto RCPathInd = *RingCagePathIndsCopy.cbegin();
							auto RingCagePath = &IntSurfRingCagePaths[RCPathInd.first][RCPathInd.second];
							int TermCPNum = RingCagePath->GetStartEndCPNum(1);
							if (TermCPNum < 0 || CageNumToElemNums.count(TermCPNum) || CPs.GetTypeFromTotOffset(TermCPNum) == CPType_CageFF) {
								CageNumToElemNums[--FarFieldCageIter] = CageWedgeElemNums;
							}
							else {
								CageNumToElemNums[TermCPNum] = CageWedgeElemNums;
							}
						}
						else{
							CageNumToElemNums[--FarFieldCageIter] = CageWedgeElemNums;
						}
					}
				}
			}

// 			Now we have all topological cage wedge element sets; must be two or more
			if (CageNumToElemNums.size() > 1) {
// 				For each cage wedge set save as a zone like a condensed basin zone,
// 					but don't need to actually compute integration totals as those will be included in data exports
// 					Zone name will be(using python f - string notation) f"{sphere name}: Cage wedge {cage cp number}",
// 						where cage cp number is the absolute cp number of the cage cp of present,
// 						or if far field it could be a negative iterator starting at - 1, or maybe just add "FF" before a positive iterator.
// 						Whatever feels best
// 				
// 				If there's exactly two "cages" but only one is local, as in cubane, then need to confirm that we got the association between wedges and cages correct with a distance comparison of the average wedge element position and the local cage CP
				if (CageNumToElemNums.size() == 2) {
					auto it = CageNumToElemNums.cbegin();
					auto it2 = it;
					it2++;
					if ((it->first >= 0 && it2->first < 0) || (it->first < 0 && it2->first >= 0)) {
						bool firstHasCP = (it->first >= 0);

						Sphere.GenerateElemMidpoints();
						auto ElemMidpoints = Sphere.GetElemMidpointsPtr();

						vec3 avgXYZ1 = { 0,0,0 }, avgXYZ2 = { 0,0,0 };
						it = CageNumToElemNums.cbegin();
						for (auto const & ei : it->second) {
							avgXYZ1 += ElemMidpoints->at(ei);
						}
						avgXYZ1 /= it->second.size();

						vec3 cageCPXYZ;
						if (it->first >= 0) {
							cageCPXYZ = CPs.GetXYZ(it->first);
						}

						it++;
						for (auto const & ei : it->second) {
							avgXYZ2 += ElemMidpoints->at(ei);
						}
						avgXYZ2 /= it->second.size();

						if (it->first >= 0) {
							cageCPXYZ = CPs.GetXYZ(it->first);
						}

						if ((firstHasCP && DistSqr(avgXYZ1, cageCPXYZ) > DistSqr(avgXYZ2, cageCPXYZ))
							|| (!firstHasCP && DistSqr(avgXYZ2, cageCPXYZ) > DistSqr(avgXYZ1, cageCPXYZ))){
							auto tmpmap = CageNumToElemNums;
							it = tmpmap.cbegin();
							it2 = it;
							it2++;

							CageNumToElemNums[it->first] = it2->second;
							CageNumToElemNums[it2->first] = it->second;
						}
					}
				}
						
				int ColorInd = 0;
				for (auto const & c : CageNumToElemNums) {
					std::map<int, int> OldToNewNodeNum;
					std::set<int> SphereNodeNums;
					vector<vector<int> > NewElems;
					for (auto const & ei : c.second){
						NewElems.push_back(Elems->at(ei));
					}

					vector<vec3> NewNodes;
					NewNodes.reserve(Sphere.GetNumNodes());
					for (auto & ei : NewElems) {
						for (auto & ni : ei){
							if (!OldToNewNodeNum.count(ni)) {
								OldToNewNodeNum[ni] = NewNodes.size();
								SphereNodeNums.insert(ni);
								NewNodes.push_back(CPPos + (SphereNodes[ni] - CPPos) * 1.01);
							}
							ni = OldToNewNodeNum[ni];
						}
					}
 					
					FESurface_c NewSurf(NewNodes, NewElems);
					string ZoneName = NucleusName +": ", CageName = "Cage ";
					int CPOffset;
					if (c.first >= 0) {
						CPOffset = CPs.GetTypeNumOffsetFromTotOffset(c.first)[1];
						CageName += to_string(CPOffset + 1);
					}
					else{
						CageName += "FF " + to_string(-c.first);
					}
					ZoneName += CageName;
					int ZoneNum = NewSurf.SaveAsTriFEZone(XYZVarNums, ZoneName);
					// 					Save same aux data as condensed basins to be used when exporting
					if (ZoneNum > 0) {
						/*
							*	Set Aux data
							*/
						AuxDataZoneSetItem(ZoneNum, CSMAuxData.GBA.SourceZoneNum, to_string(SphereZoneNum));
						AuxDataZoneSetItem(ZoneNum, CSMAuxData.GBA.SphereCPName, CPString);
						AuxDataZoneSetItem(ZoneNum, CSMAuxData.GBA.SourceNucleusName, NucleusName);
						AuxDataZoneSetItem(ZoneNum, CSMAuxData.GBA.ZoneType, CSMAuxData.GBA.ZoneTypeTopoCageWedge);
						AuxDataZoneSetItem(ZoneNum, CSMAuxData.GBA.CondensedBasinSphereElements, VectorToString(vector<int>(c.second.begin(), c.second.end()), ","));
						AuxDataZoneSetItem(ZoneNum, CSMAuxData.GBA.CondensedBasinSphereNodes, VectorToString(vector<int>(SphereNodeNums.begin(), SphereNodeNums.end()), ","));

						AuxDataZoneSetItem(ZoneNum, CSMAuxData.GBA.CondensedBasinName, CageName);
						AuxDataZoneSetItem(ZoneNum, CSMAuxData.GBA.CondensedBasinInfo, ZoneName);

						Set TmpSet(ZoneNum);
						TecUtilZoneSetMesh(SV_SHOW, TmpSet.getRef(), 0.0, FALSE);
						TecUtilZoneSetMesh(SV_COLOR, TmpSet.getRef(), 0.0, ColorIndex_t((ColorInd++ % 7) + 1));
						TecUtilZoneSetShade(SV_SHOW, TmpSet.getRef(), 0.0, TRUE);
						TecUtilZoneSetShade(SV_COLOR, TmpSet.getRef(), 0.0, ColorIndex_t((ColorInd % 7) + 1));
						TecUtilZoneSetContour(SV_SHOW, TmpSet.getRef(), 0.0, FALSE);
						TecUtilZoneSetEdgeLayer(SV_SHOW, TmpSet.getRef(), 0.0, FALSE);


					}
				}
			}

		}

		IntVarNameList = BaseIntVarNameList;
	}

	if (deleteOldSphereZones)
		TecUtilDataSetDeleteZone(oldSphereZonesToDelete.getRef());

	Set DelVars;
	if (DeleteGradVars){
		for (int i : GradVarNums)
			DelVars += i;
	}

	if (DeleteHessVars){
		for (int i : HessVarNums)
			DelVars += i;
	}

	if (!DelVars.isEmpty()){
		TecUtilDataSetDeleteVar(DelVars.getRef());
	}

	TecUtilDataLoadEnd();
	StatusDrop(AddOnID);
	return;
}


void GenerateCondensedVariables(int SphereZoneNum, int CondensedRhoVarNum) {
	int NumZones = TecUtilDataSetGetNumZones();
	int NumVars = TecUtilDataSetGetNumVars();
	REQUIRE(CondensedRhoVarNum > 0 && CondensedRhoVarNum <= NumVars);
	REQUIRE(SphereZoneNum > 0 && SphereZoneNum <= NumZones);
	REQUIRE(AuxDataZoneItemMatches(SphereZoneNum, CSMAuxData.GBA.ZoneType, CSMAuxData.GBA.ZoneTypeSphereZone));
	Boolean_t IsOk = TRUE;

	// Get sphere
	FESurface_c Sphere(SphereZoneNum, { 1,2,3 });

	// Get list of condensed variables
	string SearchStr = "I: ";
	vector<int> CondensedVarNums;
	vector<string> CondensedVarBaseNames;
	for (int i = 1; i <= NumVars; ++i){
		char* TmpCStr;
		if  (TecUtilVarGetName(i, &TmpCStr)){
			string TmpStr = TmpCStr;
			TecUtilStringDealloc(&TmpCStr);
			if (TmpStr.length() > 3 && TmpStr.substr(0, 3) == SearchStr){
				CondensedVarNums.push_back(i);
				CondensedVarBaseNames.push_back(StringRemoveSubString(TmpStr, SearchStr));
			}
		}
	}

	// Create new variables if necessary

	vector<string> NewVarNames;
	vector<int> NewVarNums;
	vector<string> VarNamePrepends = { "I: ", "INS: " };
	vector<ValueLocation_e> ZoneDataLocs(NumZones, ValueLocation_CellCentered);
	vector<FieldDataType_e> ZoneDataTypes(NumZones, FieldDataType_Bit);
	for (int i = 1; i <= NumZones; ++i) {
		if (AuxDataZoneHasItem(i, CSMAuxData.GBA.ZoneType)) {
			ZoneDataTypes[i - 1] = FieldDataType_Float;
		}
	}

	vector<string> ExtraVarBaseNames = { "Valence Density","Bader Charge","Valence Bader Charge", "Volume" };
	vector<string> FullVarBaseNames = CondensedVarBaseNames;
	FullVarBaseNames.insert(FullVarBaseNames.end(), ExtraVarBaseNames.begin(), ExtraVarBaseNames.end());

	vector<int> ExtraVarNums;

	for (string const & BaseVarName : FullVarBaseNames) {
		if (!IsOk) break;
		for (string const & Prefix : VarNamePrepends) {
			if (!IsOk) break;
			NewVarNames.push_back(Prefix + BaseVarName);
			NewVarNums.push_back(VarNumByName(NewVarNames.back()));
			if (NewVarNums.back() < 0) {
				ArgList args;
				args.appendString(SV_NAME, NewVarNames.back().c_str());
				args.appendArray(SV_VARDATATYPE, ZoneDataTypes.data());
				args.appendArray(SV_VALUELOCATION, ZoneDataLocs.data());
				IsOk = TecUtilDataSetAddVarX(args.getRef());

				if (IsOk) {
					NewVarNums.back() = TecUtilDataSetGetNumVars();
					TecUtilStateChanged(StateChange_VarsAdded, (ArbParam_t)Set(NewVarNums.back()).getRef());
				}
			}

		}
	}


	vec SphereTriangleAreas(Sphere.TriSphereElemSolidAngles());
	double TotalArea = sum(SphereTriangleAreas);
	vec SphereTriangleAreaFactors = SphereTriangleAreas / TotalArea;

	FieldDataPointer_c CondensedRhoReadPtr;
	
	TecUtilDataLoadBegin();
	CondensedRhoReadPtr.InitializeReadPtr(SphereZoneNum, CondensedRhoVarNum);


	TecUtilDataLoadEnd();
}

// void MainFunction(){
// 
// 	/*
// 	 *	Setup log
// 	 */
// 	//auto Logger = spdlog::rotating_logger_mt(LogName, "C:\\Users\\Haiiro\\Desktop\\logs\\" + LogName + ".txt", 1024 * 1024 * LogSizeMB, LogNumLogs);
// 
// 
// 	SYSTEM_INFO sysinfo;
// 	GetSystemInfo(&sysinfo);
// 
// 	int numCPU = sysinfo.dwNumberOfProcessors;
// 	omp_set_num_threads(numCPU);
// 
// 
// 	TecUtilLockStart(AddOnID);
// 	CSMGuiLock();
// 
// 	Boolean_t DoTiming = TRUE;
// 
// 	high_resolution_clock::time_point Time1, Time2;
// 	if (DoTiming){
// 		Time1 = high_resolution_clock::now();
// 	}
// 
// 
// 	Boolean_t IsOk = TRUE;
// 	int MemoryRequired;
// 	string TmpString;
// 
// 	vector<string> GoodSGPEndCPTypes = { "Nuclear", "Atom", "Bond", "Ring", "Cage" };
// 
// 	EntIndex_t NumVars = TecUtilDataSetGetNumVars();
// 
// 	/*
// 	*	Get min and max XYZ for the system, just in case we need it.
// 	*/
// 
// 
// 	vector<EntIndex_t> XYZVarNums(3);
// 	TecUtilAxisGetVarAssignments(&XYZVarNums[0], &XYZVarNums[1], &XYZVarNums[2]);
// 	for (int i = 0; i < 3 && IsOk; ++i)
// 		IsOk = (XYZVarNums[i] > 0);
// 
// 
// 	vec3 VolMinXYZ, VolMaxXYZ;
// 	vector<int> VolMaxIJK(3);
// 	EntIndex_t VolZoneNum = ZoneNumByName("Full Volume");
// 	if (IsOk){
// 		if (VolZoneNum <= 0){
// 			VolZoneNum = 1;
// 			// 			TecUtilDialogErrMsg("Couldn't get volume zone");
// 			//TecUtilLockFinish(AddOnID);
// 			// 			return;
// 		}
// // 		for (int i = 0; i < 3; ++i){
// // 			TecUtilDataValueGetMinMaxByZoneVar(VolZoneNum, XYZVarNums[i], &VolMinXYZ[i], &VolMaxXYZ[i]);
// // 		}
// 		ZoneXYZVarGetMinMax_Ordered3DZone(XYZVarNums, VolZoneNum, VolMinXYZ, VolMaxXYZ);
// 		TecUtilZoneGetIJK(VolZoneNum, &VolMaxIJK[0], &VolMaxIJK[1], &VolMaxIJK[2]);
// 	}
// 
// 	vector<FieldDataType_e> VarDataTypes(NumVars);
// 	vector<ValueLocation_e> VarLocations(NumVars, ValueLocation_Nodal);
// 	
// 	for (int i = 0; i < NumVars; ++i){
// 		VarDataTypes[i] = TecUtilDataValueGetType(VolZoneNum, i + 1);
// 	}
// 
// 
// 	/*
// 	*	Check to see if there's any existing GBA zones present.
// 	*	If there are, then the data types and data locations need
// 	*	to be consistent with those zones
// 	*	i.e. need to use cell-centered data for integrated variables.
// 	*
// 	*	Just check for a GBA sphere zone, since that's sufficient.
// 	*/
// 
// 	for (int i = 1; i <= TecUtilDataSetGetNumZones(); ++i){
// 		if (AuxDataZoneHasItem(i, CSMAuxData.GBA.ZoneType)){
// 			for (int j = 3; j < NumVars; ++j){
// 				VarDataTypes[j] = TecUtilDataValueGetType(i, j + 1);
// 				VarLocations[j] = TecUtilDataValueGetLocation(i, j + 1);
// 			}
// 			break;
// 		}
// 	}
// 
// 	int GroupNum = 0;
// 
// 	/*
// 	 *	Get job parameters from dialog
// 	 */
// 	Boolean_t UseCutoff = TRUE;
// 	UseCutoff = TecGUIToggleGet(TGLOpenSys_TOG_T1_1);
// 	double CutoffVal = 0.001;
// 	TecGUITextFieldGetDouble(TFCutoff_TF_T1_1, &CutoffVal);
// 	EntIndex_t CutoffVarNum = VarNumByName(string("Electron Density"));
// 	LgIndex_t NumEdgeGPs = TecGUIScaleGetValue(SCNumEdgeGPs_SC_T1_1);
// 
// 	/*
// 	 *	When checking if two streamtraces are straddling two IBs,
// 	 *	use the directions of their ends and the distance between
// 	 *	their ends to check.
// 	 *	IBCheckDistRatio is the fraction of the sphere radius that,
// 	 *	if the streamtrace endpoints are closer than,
// 	 *	they are considered to be at the same point and not in
// 	 *	different IBs.
// 	 *	IBCheckAngle is the angle, in degrees, that if greater than the
// 	 *	angle between the end segments of the streamtraces, they are
// 	 *	considered to be parallel and considered in the same IB.
// 	 */
// // 	double IBCheckDistRatio = 0.2;
// // 	TecGUITextFieldGetDouble(TFIBDist_TF_T1_1, &IBCheckDistRatio);
// // 	double IBCheckAngle = 20.0;
// // 	TecGUITextFieldGetDouble(TFIBAng_TF_T1_1, &IBCheckAngle);
// 
// 	LgIndex_t * CPNums = nullptr;
// 	LgIndex_t NumSelectedCPs = 1;
// 	TecGUIListGetSelectedItems(MLSelCPs_MLST_T1_1, &CPNums, &NumSelectedCPs);
// 	/*
// 	*	Set sphere radius and mesh refinement level
// 	*	(will be set by user in future)
// 	*/
// 
// 	
// 
// 
// 	RadMode RadiusMode = MINCPDISTRATIO;
// 	RadiusMode = (RadMode)TecGUIRadioBoxGetToggle(RBRadMode_RADIO_T1_1);
// 	double UserRadius = 0.25;
// 	TecGUITextFieldGetDouble(TFRad_TF_T1_1, &UserRadius);
// 	int Level = 3;
// 	TecGUITextFieldGetLgIndex(TFLevel_TFS_T1_1, &Level);
// 	LgIndex_t NumSTPoints = 100;
// 	TecGUITextFieldGetLgIndex(TFSTPts_TF_T1_1, &NumSTPoints);
// 
// 
// 	bool SaveGPs = TecGUIToggleGet(TGLsGP_TOG_T1_1);
// 	bool SaveGBs = TecGUIToggleGet(TGLsGB_TOG_T1_1);
// 
// 
// 	VolExtentIndexWeights_s VolInfo;
// 	GetVolInfo(VolZoneNum, XYZVarNums, FALSE, VolInfo);
// 
// 
// 	EntIndex_t OldNumZones = TecUtilDataSetGetNumZones();	
// 
// 	/*
// 	*	information for integration
// 	*/
// 	vector<string> AtomNameList;
// 	vector<string> IntVarNameList;
// 	vector<int> IntVarNumList;
// 	int IntResolution = TecGUIScaleGetValue(SCPrecise_SC_T1_1);
// 	vector<FieldDataPointer_c> IntVarPtrs;
// 
// 	Boolean_t DoIntegration = TecGUIToggleGet(TGLInt_TOG_T1_1);
// 
// 	if (DoIntegration){
// 		AtomNameList = ListGetSelectedStrings(MLSelCPs_MLST_T1_1);
// 		IntVarNameList = ListGetSelectedStrings(MLSelVars_MLST_T1_1);
// 		IntVarNumList = ListGetSelectedItemNums(MLSelVars_MLST_T1_1);
// 
// 		IntVarPtrs.resize(IntVarNumList.size());
// 		for (int i = 0; i < IntVarPtrs.size(); ++i){
// 			IsOk = IntVarPtrs[i].InitializeReadPtr(VolZoneNum, IntVarNumList[i]);
// 		}
// 	}
// 	
// 
// 	if (!IsOk){
// 		TecUtilDrawGraphics(TRUE);
// 		TecUtilDialogErrMsg("Couldn't get XYZ vars.");
// 		TecUtilLockFinish(AddOnID);
// 		return;
// 	}
// 
// 	EntIndex_t CPTypeVarNum;
// 	CPTypeVarNum = VarNumByName("CritPointType");
// 
// 	Set oldSphereZonesToDelete;
// 	bool deleteOldSphereZones = false, deleteOldSphereZonesAsked = false;
// 
// 	vector<int> NumCPs;
// 	for (int SelectCPNum = 0; SelectCPNum < NumSelectedCPs && IsOk; ++SelectCPNum){
// 
// 		LgIndex_t CPNum = -1;
// 		int CPType = -3;
// 
// 		double Radius = 0.25;
// 
// 		stringstream ProgressStr;
// 		string CPString, CPName, CPTypeStr;
// 		bool IsCP;
// 
// 		int RealCPNum = 0;
// 		EntIndex_t CPZoneNum;
// 		LgIndex_t CPIJKMax[3];
// 
// 		vec3 CPPos;
// 		string NucleusName;
// 		if (IsOk){
// 			CPPos = GetCoordsFromListItem(CPNums[SelectCPNum], MLSelCPs_MLST_T1_1, &CPString, &CPName, &CPNum, &CPZoneNum, &IsCP, &NumCPs);
// 			NucleusName = GetNucleusNameForCP(CPPos);
// 			if (IsCP){
// 
// // 				CPZoneNum = ZoneNumByName("Critical Points");
// 
// // 				LgIndex_t NumCPs[4];
// // 				if (IsOk){
// // 					for (int i = 0; i < 4; ++i){
// // 						NumCPs[i] = stoi(AuxDataZoneGetItem(CPZoneNum, CCDataNumCPs[i]));
// // 					}
// // 				}
// // 
// 
// 				TecUtilZoneGetIJK(CPZoneNum, &CPIJKMax[0], &CPIJKMax[1], &CPIJKMax[2]);
// 				IsOk = (CPIJKMax[0] > 1 && CPIJKMax[1] == 1 && CPIJKMax[2] == 1);
// 				if (!IsOk){
// 					TecUtilDrawGraphics(TRUE);
// 					TecUtilDialogErrMsg("Critical point zone dimensions incorrect.");
// 					TecUtilLockFinish(AddOnID);
// 					return;
// 				}
// 
// 				/*
// 				*	Get XYZ readable pointers to the crit points zone to get the
// 				*	position of the fifth atom (working with cubane, so this is the first carbon).
// 				*/
// 
// // 				if (CPZoneNum <= 0){
// // 					TecUtilDrawGraphics(TRUE);
// // 					TecUtilDialogErrMsg("Couldn't get critical points zone.");
// // 					TecUtilLockFinish(AddOnID);
// // 					return;
// // 				}
// 
// 				for (int i = 0; i < 4; ++i){
// 					if (CPName == RankStrs[i] || i == 0 && CPName == "Atom"){
// 						int RealCPNum = 0;
// 						for (int j = 0; j < i; ++j){
// 							RealCPNum += NumCPs[j];
// 						}
// 						CPNum = RealCPNum + CPNum;
// 						break;
// 					}
// 				}
// 
// 
// 				CPType = (int)TecUtilDataValueGetByZoneVar(CPZoneNum, CPTypeVarNum, CPNum);
// 				for (int i = 0; i < CPTypeList.size(); ++i){
// 					if (CPType == CPTypeList[i]){
// 						CPTypeStr = CPNameList[i];
// 					}
// 				}
// 			}
// 
// 			/*
// 			 *	Delete any existing sphere for the current CP if so chosen by the user
// 			 */
// 			for (int z = 1; z <= TecUtilDataSetGetNumZones(); ++z){
// 				if (AuxDataZoneItemMatches(z, CSMAuxData.GBA.SphereCPName, CPString) && (!deleteOldSphereZonesAsked || deleteOldSphereZones)){
// 					string tmpStr = "An existing sphere zone was found for " + CPString + ". Would you like to erase it (and all other conflicting sphere zones) before continuing?";
// 					if (!deleteOldSphereZonesAsked){
// 						deleteOldSphereZones = TecUtilDialogMessageBox(tmpStr.c_str(), MessageBox_YesNo);
// 						deleteOldSphereZonesAsked = true;
// 					}
// 					if (deleteOldSphereZones) {
// 						oldSphereZonesToDelete += z;
// 						for (int zi = z + 1; zi <= TecUtilDataSetGetNumZones(); zi++){
// 							if (AuxDataZoneItemMatches(zi, CSMAuxData.GBA.SourceZoneNum, to_string(z)))
// 								oldSphereZonesToDelete += zi;
// 						}
// 					}
// 					break;
// 				}
// 			}
// 
// 			ProgressStr << "Processing " << (NucleusName == "" ? CPString : NucleusName);
// 			if (NumSelectedCPs > 1){
// 				ProgressStr << " ... CP " << SelectCPNum + 1 << " of " << NumSelectedCPs;
// 			}
// 
// // 			TecUtilStringDealloc(&TempCStr);
// 		}
// 		else{
// 			TecUtilDialogErrMsg("Didn't get CP number correctly");
// 		}
// 
// 		if (IsOk){
// 			TmpString = ProgressStr.str() + string(" ... Step 1 of 4 ... Generating mesh");
// 			StatusLaunch(TmpString, AddOnID, FALSE);
// 		}
// 
// // 		for (int i = 0; i < 3; ++i){
// // 			CPPos[i] = TecUtilDataValueGetByZoneVar(CPZoneNum, XYZVarNums[i], CPNum);
// // 		}
// 
// 		vector<point> IntersectionPoints;
// 		vector<vector<LgIndex_t> > IntZoneSaddleCPNodeNums;
// 		vector<vector<LgIndex_t> > AllIntCPNodeNums;
// 		vector<string> AllIntVolumeCPNames;
// 
// 		vector<vec3> IntCPPos;
// 
// 		EntIndex_t NumZones = TecUtilDataSetGetNumZones();
// 
// 		point* p = nullptr;
// 		triangle* t = nullptr;
// 		int** e = nullptr;
// 		int NumPoints, NumTriangles, NumEdges;
// 
// 		vector<int> MovedPointNums;
// 		string MovedPointCPTypes, MovedPointCPNames, MovedPointCPNamesTotalCount;
// 		vector<bool> MovedCPIsBond;
// 
// 		MeshStatus_e MeshStatus = FAIL_INVALIDCONSTRAINT;
// 
// 
// 		double ClosestCPDist = 1e50;
// 		if (IsOk){
// 			if (IsCP){
// 				for (int i = 1; i <= CPIJKMax[0] && IsOk; ++i){
// 					if (i != CPNum && (int)TecUtilDataValueGetByZoneVar(CPZoneNum, CPTypeVarNum, i) != CPType){
// 						vec3 OtherCP;
// 						for (int j = 0; j < 3; ++j){
// 							OtherCP[j] = TecUtilDataValueGetByZoneVar(CPZoneNum, XYZVarNums[j], i);
// 						}
// 						ClosestCPDist = MIN(ClosestCPDist, Distance(OtherCP, CPPos));
// 					}
// 				}
// 			}
// 			else{
// 				for (int ZoneNum = 1; ZoneNum <= TecUtilDataSetGetNumZones(); ++ZoneNum){
// 					char* ZoneNameCStr;
// 					if (TecUtilZoneIsOrdered(ZoneNum) && TecUtilZoneGetName(ZoneNum, &ZoneNameCStr)){
// 						string ZoneName = ZoneNameCStr;
// 						TecUtilStringDealloc(&ZoneNameCStr);
// 
// 						int ElemNum = SearchVectorForString(ElementSymbolList, ZoneName, false);
// 						if (ElemNum < 0)
// 							ElemNum = SearchVectorForString(ElementNameList, ZoneName);
// 						if (ElemNum >= 0){
// 							int NumPts;
// 							TecUtilZoneGetIJK(ZoneNum, &NumPts, nullptr, nullptr);
// 							for (int i = 1; i <= NumPts; ++i){
// 								if (ZoneNum != CPZoneNum || i != CPNum){
// 									vec3 OtherCP;
// 									for (int j = 0; j < 3; ++j){
// 										OtherCP[j] = TecUtilDataValueGetByZoneVar(ZoneNum, XYZVarNums[j], i);
// 									}
// 									ClosestCPDist = MIN(ClosestCPDist, Distance(OtherCP, CPPos));
// 								}
// 							}
// 						}
// 					}
// 				}
// 			}
// 		}
// 
// 		if (RadiusMode == ABSOLUTERADIUS || ClosestCPDist == 1e50 || ClosestCPDist <= 0){
// 			Radius = UserRadius;
// 		}
// 		else if (RadiusMode == MINCPDISTRATIO){
// 			Radius = UserRadius * ClosestCPDist;
// 		}
// 
// 		vector<int> XYZRhoVarNums = XYZVarNums;
// 		XYZRhoVarNums.push_back(CutoffVarNum);
// 
// 		vector<GradPath_c> IntGPs;
// 
// 
// 		/*   
// 		 *	Find all volumetric special gradient paths that will intersect
// 		 *	the sphere.
// 		 */
// 		TecUtilDataLoadBegin();
// 		for (EntIndex_t CurZoneNum = 1; CurZoneNum <= OldNumZones && IsOk && IsCP; ++CurZoneNum){
// 			LgIndex_t TmpIJK[3];
// 
// 			Boolean_t ZoneOK = TRUE;
// 			Boolean_t ZoneIsActive = TecUtilZoneIsActive(CurZoneNum);
// 			if (!ZoneIsActive){
// 				TecUtilZoneSetActive(Set(CurZoneNum).getRef(), AssignOp_PlusEquals);
// // 				Set_pa TempSet = TecUtilSetAlloc(FALSE);
// // 				TecUtilSetAddMember(TempSet, CurZoneNum, FALSE);
// // 				TecUtilZoneSetActive(TempSet, AssignOp_PlusEquals);
// // 				TecUtilSetDealloc(&TempSet);
// 			}
// 			ZoneOK = !AuxDataZoneHasItem(CurZoneNum, CSMAuxData.GBA.ZoneType);
// 			if (ZoneOK)
// 				TecUtilZoneGetIJK(CurZoneNum, &TmpIJK[0], &TmpIJK[1], &TmpIJK[2]);
// 
// 			if (ZoneOK && TmpIJK[0] > 1 && TmpIJK[1] == 1 && TmpIJK[2] == 1){
// 				// If 1-d (line) zone that wasn't made by GBA.
// 				// Check that either initial or terminal point is within sphere radius
// 				FieldData_pa TmpXYZPtrs[3];
// 				for (int i = 0; i < 3 && IsOk; ++i){
// 					TmpXYZPtrs[i] = TecUtilDataValueGetReadableRef(CurZoneNum, XYZVarNums[i]);
// 					IsOk = VALID_REF(TmpXYZPtrs[i]);
// 				}
// 				Boolean_t ConnectsToCP = FALSE;
// 				Boolean_t StartsFrom1;
// 				LgIndex_t StartEndCPs[2];
// 				string StartEndTypes[2];
// 				for (int i = 0; i < 2 && ZoneOK; ++i){
// 					ZoneOK = AuxDataZoneGetItem(CurZoneNum, CCDataGPEndNums[i], TmpString);
// 					if (ZoneOK){
// 						StartEndCPs[i] = stoi(TmpString);
// 						ZoneOK = AuxDataZoneGetItem(CurZoneNum, CCDataGPEndTypes[i], TmpString);
// 						if (ZoneOK)
// 							ZoneOK = (SearchVectorForString(GoodSGPEndCPTypes, TmpString) >= 0);
// 						if (ZoneOK)
// 							StartEndTypes[i] = TmpString;
// 					}
// 				}
// 				if (!ZoneOK){
// 					ZoneOK = TRUE;
// 					for (int i = 0; i < 2 && ZoneOK; ++i){
// 						ZoneOK = AuxDataZoneGetItem(CurZoneNum, CSMAuxData.CC.GPEndNumStrs[i], TmpString);
// 						if (ZoneOK){
// 							StartEndCPs[i] = stoi(TmpString);
// 							ZoneOK = AuxDataZoneGetItem(CurZoneNum, CSMAuxData.CC.GPEndTypes[i], TmpString);
// 							if (ZoneOK)
// 								ZoneOK = (SearchVectorForString(GoodSGPEndCPTypes, TmpString) >= 0);
// 							if (ZoneOK)
// 								StartEndTypes[i] = TmpString;
// 						}
// 					}
// 				}
// 
// 				Boolean_t IntValid = ZoneOK && (CPTypeStr == "" || StartEndTypes[0] == CPTypeStr || StartEndTypes[1] == CPTypeStr) && StartEndTypes[0] != StartEndTypes[1];
// 				
// 				if (ZoneOK){
// 					Boolean_t IntRecorded = FALSE;
// 					for (int i = 0; i < 2 && IsOk; ++i){
// 						if (StartEndCPs[i] == CPNum){
// 							ConnectsToCP = TRUE;
// 							StartsFrom1 = (i == 0);
// // 							for (int j = 2; j < 4 && IsOk && IntValid; ++j){
// 							for (int j = 1; j < 4 && IsOk && IntValid; ++j){
// 								if (StartEndTypes[(i + 1) % 2] == RankStrs[j]){
// 									GradPath_c TmpGP(CurZoneNum, XYZRhoVarNums, AddOnID);
// 									if (Distance(CPPos, TmpGP[0]) < Distance(CPPos, TmpGP[-1])){
// 										TmpGP.Reverse();
// 									}
// 									IsOk = TmpGP.Trim(CPPos, Radius);
// 									if (!IsOk){
// 										TecUtilDialogErrMsg("Failed to trim existing GP");
// 									}
// 
// 									/*
// 									 *	Now check the intersection point against the previously found intersections
// 									 *	to make sure there are not two intersections at the same point
// 									 */
// 
// 									for (GradPath_c const & OldGP : IntGPs){
// 										double IntIntDist = Distance(TmpGP[-1], OldGP[-1]);
// 										if (IntIntDist < 0.1){
// 											IntValid = FALSE;
// 										}
// 										// 										double a = 1;
// 									}
// 
// 									if (IntValid){
// 
// 										IntGPs.push_back(TmpGP);
// 
// 										vec3 TmpCP;
// 										for (int k = 0; k < 3 && IsOk; ++k){
// 											TmpCP[k] = TecUtilDataValueGetByZoneVar(CPZoneNum, XYZVarNums[k], StartEndCPs[(i + 1) % 2]);
// 										}
// 										IntZoneSaddleCPNodeNums.push_back(vector<LgIndex_t>());
// 										IntZoneSaddleCPNodeNums.back().push_back(static_cast<int>(IntersectionPoints.size()));
// 										IntZoneSaddleCPNodeNums.back().push_back(CurZoneNum);
// 										IntZoneSaddleCPNodeNums.back().push_back(StartEndCPs[(i + 1) % 2]);
// 										IntCPPos.push_back(TmpCP);
// 									}
// 
// 									// 									if (TmpIJK[0] < NumSTPoints){
// 									// 										NumSTPoints = TmpIJK[0];
// 									// 										if (NumSTPoints % 2 != 0){
// 									// 											NumSTPoints--;
// 									// 										}
// 									// 									}
// 									if (IntValid){
// 										int CPCount = 0;
// 										for (int k = 0; k < RankStrs.size(); ++k){
// 											if (StartEndTypes[(i + 1) % 2] == RankStrs[k] && !IntRecorded){
// 												AllIntVolumeCPNames.push_back(StartEndTypes[(i + 1) % 2] + string(" ") + to_string(StartEndCPs[(i + 1) % 2] - CPCount));
// 												AllIntCPNodeNums.push_back(vector<LgIndex_t>());
// 												AllIntCPNodeNums.back().push_back(static_cast<int>(IntersectionPoints.size()));
// 												AllIntCPNodeNums.back().push_back(StartEndCPs[(i + 1) % 2]);
// 												Boolean_t ZoneFound = FALSE;
// // 												if (StartEndTypes[(i + 1) % 2] == RankStrs[1]){
// // 													// if the intersecting SGP connects to a bond point,
// // 													// then find the CP on the other side of the bond path.
// // 
// // 													// The other end of the bond path should be the zone before of after CurZoneNum, so check those first
// // 													for (int CheckNum = -1; CheckNum < 2; CheckNum += 2){
// // 														if (CurZoneNum + CheckNum <= NumZones && AuxDataZoneGetItem(CurZoneNum + CheckNum, CCDataGPEndNums[0], TmpString)){
// // 															// the beginning point for bond paths is always the bond point, so only check CCDataGPEndNums[0]
// // 															if (stoi(TmpString) == StartEndCPs[(i + 1) % 2]){
// // 																// Matching bond path found, so get nuclear CP on other side
// // 																int CheckNucCPNum = stoi(AuxDataZoneGetItem(CurZoneNum + CheckNum, CCDataGPEndNums[1]));
// // 																// now search through the nuclear zones, if present, to see what type of element the CP is
// // 																if (CheckNucCPNum > 0){
// // 																	double MinNucCPCheckDist = DBL_MAX, TempDist;
// // 																	string ClosestElementName;
// // 																	int ClosestAtomNum = -1;
// // 																	for (int CheckZoneNum = 1; CheckZoneNum <= NumZones; ++CheckZoneNum){
// // 																		if (CheckZoneNum != CurZoneNum
// // 																			&& CheckZoneNum != CurZoneNum + CheckNum
// // 																			&& TecUtilZoneIsOrdered(CheckZoneNum)
// // 																			&& AuxDataZoneItemMatches(CheckZoneNum, DLZoneType, DLZoneTypeNuclearPositions)){
// // 																			// need to now find the 
// // 																			int NucIJK[3];
// // 																			TecUtilZoneGetIJK(CheckZoneNum, &NucIJK[0], &NucIJK[1], &NucIJK[2]);
// // 																			vec3 NucPt, NucCPPt;
// // 																			for (int Dir = 0; Dir < 3; ++Dir){
// // 																				NucCPPt[Dir] = TecUtilDataValueGetByZoneVar(CurZoneNum + CheckNum, XYZVarNums[Dir], CheckNucCPNum);
// // 																			}
// // 																			for (int NucCheckPt = 1; NucCheckPt <= NucIJK[0]; ++NucCheckPt){
// // 																				for (int Dir = 0; Dir < 3; ++Dir){
// // 																					NucPt[Dir] = TecUtilDataValueGetByZoneVar(CheckZoneNum, XYZVarNums[Dir], NucCheckPt);
// // 																				}
// // 																				TempDist = DistSqr(NucPt, NucCPPt);
// // 																				if (TempDist < MinNucCPCheckDist){
// // 																					MinNucCPCheckDist = TempDist;
// // 																					ClosestElementName = AuxDataZoneGetItem(CheckZoneNum, DLZoneAtomicSpecies);
// // 																					ClosestAtomNum = NucCheckPt;
// // 																				}
// // 																			}
// // 																		}
// // 																	}
// // 																	if (ClosestAtomNum > 0){
// // 																		MovedPointCPNames += ClosestElementName + " " + to_string(ClosestAtomNum) + ",";
// // 																		ZoneFound = TRUE;
// // 																	}
// // 																}
// // 															}
// // 														}
// // 													}
// // 												}
// 												MovedPointCPTypes += StartEndTypes[(i + 1) % 2] + ",";
// 												MovedCPIsBond.push_back(StartEndTypes[(i + 1) % 2] == RankStrs[1]);
// 												if (!ZoneFound) {
// 													MovedPointCPNames += StartEndTypes[(i + 1) % 2] + string(" ") + to_string(StartEndCPs[(i + 1) % 2] - CPCount) + ",";
// 													MovedPointCPNamesTotalCount += StartEndTypes[(i + 1) % 2] + string(" ") + to_string(StartEndCPs[(i + 1) % 2]) + ",";
// 												}
// 												IntRecorded = TRUE;
// 												break;
// 											}
// // 											CPCount += NumCPs[k == 0 ? k : k - 1];
// 											CPCount += NumCPs[k];
// 										}
// 									}
// 								}
// 							}
// 							break;
// 						}
// 					}
// 				}
// 				if (ConnectsToCP && IsOk && IntValid){
// 
// 					// An intersection exists, so need to find its location. 
// 					// INTERSECTION LOCATION IS THE SAME AS THE TRIMMED (LAST) POINT OF IntGPs[-1]
// 					// 
// 					point IntPt;
// 					for (int i = 0; i < 3; ++i){
// 						IntPt[i] = IntGPs.back()[-1][i] - CPPos[i];
// 					}
// 					/*
// 					*	Make sure this intersection isn't too close to another
// 					*/
// 
// 					Boolean_t IntersectionGood = TRUE;
// 					for (int i = 0; i < IntersectionPoints.size() && IntersectionGood; ++i){
// 						vec3 TmpPt1, TmpPt2;
// 						for (int j = 0; j < 3; ++j){
// 							TmpPt1[j] = IntPt[j];
// 							TmpPt2[j] = IntersectionPoints[i][j];
// 						}
// // 						TODO: The distance cutoff of 0.1 is completely arbitrary. Make it less so.
// 						IntersectionGood = (Distance(TmpPt1, TmpPt2) >= 0.1);
// 					}
// 
// 					if (IntersectionGood)
// 						IntersectionPoints.push_back(IntPt);
// 
// // 					Below is the (spurious) code to find the intersection of the GP with the sphere.
// // 					This is unnecessary because the point was already found when the GP was trimmed above.
// 
// // 					vec3 Pt1, Pt2;
// // 					LgIndex_t StartPt = 1, Step = 1;
// // 					if (!StartsFrom1){
// // 						StartPt = TmpIJK[0];
// // 						Step = -1;
// // 					}
// // 					Boolean_t IntersectionFound = FALSE;
// // 					double Dist1, Dist2;
// // 					for (int i = 0; i < 3; ++i){
// // 						Pt1[i] = TecUtilDataValueGetByRef(TmpXYZPtrs[i], StartPt + Step);
// // 					}
// // 					Pt2 = Pt1;
// // 					Dist2 = Dist1 = Distance(Pt1, CPPos);
// // 					StartPt += Step;
// // 					for (LgIndex_t PtNum = StartPt + Step; PtNum >= 1 && PtNum <= TmpIJK[0] && !IntersectionFound && IsOk; PtNum += Step){
// // 						for (int i = 0; i < 3; ++i){
// // 							Pt1[i] = TecUtilDataValueGetByRef(TmpXYZPtrs[i], PtNum);
// // 						}
// // 						Dist1 = Distance(Pt1, CPPos);
// // 						if ((Dist1 < Radius && Dist2 >= Radius) ||
// // 							(Dist1 >= Radius && Dist2 < Radius)){
// // 							// Now Pt1 and Pt2 are straddling the sphere.
// // 							// Interpolate to find the intersection point.
// // 							double DistanceRatio = (Radius - Dist1) / (Dist2 - Dist1);
// // 							point IntPt;
// // 							for (int i = 0; i < 3; ++i){
// // 								IntPt[i] = (Pt1[i] + DistanceRatio * (Pt2[i] - Pt1[i])) - CPPos[i];
// // 							}
// // 
// // 							/*
// // 							* Project intersection point to sphere surface
// // 							*/
// // 							// Compute spherical coordinates
// // 							double   rad = sqrt(pow(IntPt.x, 2) + pow(IntPt.y, 2) + pow(IntPt.z, 2));
// // 							double theta = acos(IntPt.z / rad);
// // 							double   phi = atan2(IntPt.y, IntPt.x);
// // 
// // 							// Project point onto a sphere of radius "radius"
// // 							IntPt.x = Radius * sin(theta) * cos(phi);
// // 							IntPt.y = Radius * sin(theta) * sin(phi);
// // 							IntPt.z = Radius * cos(theta);
// // 
// // 
// // 							IntersectionFound = TRUE;
// // 							/*
// // 							 *	Make sure this intersection isn't too close to another
// // 							 */
// // 
// // 							Boolean_t IntersectionGood = TRUE;
// // 							for (int i = 0; i < IntersectionPoints.size() && IntersectionGood; ++i){
// // 								vec3 TmpPt1, TmpPt2;
// // 								for (int j = 0; j < 3; ++j){
// // 									TmpPt1[j] = IntPt[j];
// // 									TmpPt2[j] = IntersectionPoints[i][j];
// // 								}
// // 								IntersectionGood = (Distance(TmpPt1, TmpPt2) >= 0.1);
// // 							}
// // 
// // 							if (IntersectionGood)
// // 								IntersectionPoints.push_back(IntPt);
// // 
// // 
// // 						}
// // 						Pt2 = Pt1;
// // 						Dist2 = Dist1;
// // 					}
// 				}
// 			}
// 			if (!ZoneIsActive && IsOk){
// 				Set_pa TmpSet = TecUtilSetAlloc(FALSE);
// 				IsOk = TecUtilSetAddMember(TmpSet, CurZoneNum, FALSE);
// 				if (IsOk)
// 					IsOk = (TecUtilZoneSetActive(TmpSet, AssignOp_MinusEquals) <= 1);
// 				TecUtilSetDealloc(&TmpSet);
// 				if (!IsOk)
// 					TecUtilDialogErrMsg("Failed to deactivate zone during intersection check");
// 			}
// 		}
// 		TecUtilDataLoadEnd();
// 
// 
// 		/*
// 		 *	Get optimized spherical mesh
// 		 */
// 
// 		int MaxRefine = 2;
// 		int RefineNum = 0;
// 		TecUtilDataLoadBegin();
// 		while (MeshStatus != 0){
// 			/*
// 			 *	Get memory usage estimate for mesh
// 			 */
// 			NumTriangles = 20 * (int)pow(4, Level);
// 			NumPoints = NumTriangles / 2 + 2;
// 			NumEdges = int(NumTriangles * double(3.0 / 2.0));
// 			/*
// 			*	point and triangle adjacency lists in meshgen.
// 			*	Typically have six neighbor point and triangles
// 			*	per point.
// 			*/
// 			int AdjListSize = sizeof(int) * 2 * 6 * NumPoints;
// 			MemoryRequired = (sizeof(point) * NumPoints
// 				+ sizeof(triangle) * NumTriangles + AdjListSize
// 				+ sizeof(int) * 2 * NumEdges) / 1024;
// 			TecUtilMemoryChangeNotify(MemoryRequired);
// 
// 			MeshStatus = meshgen2D_sphere(Radius, Level, IntersectionPoints, MovedPointNums, p, t, e, NumPoints, NumTriangles, NumEdges);
// 
// 			/*
// 			 *	The point and triangle lists are kept. only
// 			 *	the adjacency lists are destroyed.
// 			 */
// 			MemoryRequired = (AdjListSize) / 1024;
// 			TecUtilMemoryChangeNotify(-MemoryRequired);
// 
// 
// 			if (MeshStatus == FAIL_NEEDREFINEMENT){
// 				//break;
// 				if (RefineNum >= MaxRefine){
// 					break;
// 				}
// 				else{
// 					//TecUtilDialogMessageBox("Mesh not fine enough for constrained points. Refining mesh.", MessageBoxType_Warning);
// 					TecUtilMemoryChangeNotify(-((int)sizeof(point) * NumPoints
// 						+ (int)sizeof(triangle) * NumTriangles) / 1024);
// 					delete p, t;
// 					p = nullptr;
// 					t = nullptr;
// 					e = nullptr;
// 					Level++;
// 					RefineNum++;
// 				}
// 			}
// 			else if (MeshStatus == FAIL_INVALIDCONSTRAINT){
// 				TecUtilDrawGraphics(TRUE);
// 				TecUtilDataLoadEnd();
// 				TecUtilDialogErrMsg("Constrained point too far from mesh");
// 				TecUtilLockFinish(AddOnID);
// 				return;
// 			}
// 			else if (MeshStatus == SUCCESS_NOTCONVERGED){
// 				// 				TecUtilDialogMessageBox("Mesh did not converge.", MessageBoxType_Warning);
// 				break;
// 			}
// 
// 		}
// 		TecUtilDataLoadEnd();
// 
// 		IsOk = (MeshStatus == SUCCESS || MeshStatus == SUCCESS_NOTCONVERGED && NumPoints > 0 && NumTriangles > 0);
// 		if (!IsOk){
// 			TecUtilDrawGraphics(TRUE);
// 			TecUtilDialogErrMsg("Failed to make mesh.");
// 			TecUtilLockFinish(AddOnID);
// 			return;
// 		}
// 
// 		/*
// 		 *	Get a list of edge numbers for each triangle
// 		 */
// 		vector<vector<int> > TriangleEdgeNums(NumTriangles), NodeEdgeNums(NumPoints);
// 		for (auto & i : TriangleEdgeNums) i.reserve(3);
// #pragma omp parallel for
// 		for (int TriNum = 0; TriNum < NumTriangles; ++TriNum){
// 			for (int tNode = 0; tNode < 3; ++tNode){
// 				bool EdgeFound = false;
// 				for (int EdgeNum = 0; EdgeNum < NumEdges && !EdgeFound; ++EdgeNum){
// 					for (int eNode = 0; eNode < 2; ++eNode){
// 						if (e[EdgeNum][eNode] == t[TriNum][tNode] && e[EdgeNum][(eNode + 1) % 2] == t[TriNum][(tNode + 1) % 3]){
// 							TriangleEdgeNums[TriNum].push_back(EdgeNum);
// 							EdgeFound = true;
// 						}
// 					}
// 				}
// 			}
// 		}
// 		for (auto & i : NodeEdgeNums) i.reserve(6);
// #pragma omp parallel for
// 		for (int NodeNum = 0; NodeNum < NumPoints; ++NodeNum){
// 			for (int EdgeNum = 0; EdgeNum < NumEdges; ++EdgeNum){
// 				for (int eNode = 0; eNode < 2; ++eNode){
// 					if (e[EdgeNum][eNode] == NodeNum && std::find(NodeEdgeNums[NodeNum].begin(), NodeEdgeNums[NodeNum].end(), EdgeNum) == NodeEdgeNums[NodeNum].end()){
// 						NodeEdgeNums[NodeNum].push_back(EdgeNum);
// 					}
// 				}
// 			}
// 		}
// 
// 		FESurface_c Sphere;
// 		EntIndex_t SurfZoneNum;
// 
// 		if (IsOk) {
// 			/*
// 			 *	Translate sphere to the CP position
// 			 */
// 
// 			TecUtilDataLoadBegin();
// 			for (int i = 0; i < NumPoints; ++i) {
// 				for (int j = 0; j < 3; ++j) {
// 					p[i][j] += CPPos[j];
// 				}
// 			}
// 			TecUtilDataLoadEnd();
// 
// 			/*
// 			*	Add triangle finite element zone for sphere
// 			*/
// 
// 			vector<vec3> P(NumPoints);
// 			for (int i = 0; i < NumPoints; ++i) for (int j = 0; j < 3; ++j) P[i][j] = p[i][j];
// 			vector<vector<int> > T(NumTriangles, vector<int>(3));
// 			for (int i = 0; i < NumTriangles; ++i) for (int j = 0; j < 3; ++j) T[i][j] = t[i][j];
// 
// 			IsOk = Sphere.MakeFromNodeElemList(P, T);
// 		}
// 
// 
// 		if (IsOk){
// 			vector<ValueLocation_e> DataLocs(TecUtilDataSetGetNumVars(), ValueLocation_Nodal);
// 			vector<FieldDataType_e> DataTypes(TecUtilDataSetGetNumVars(), FieldDataType_Double);
// 			vector<string> IntVarPrefixes = { "I: ","IN: ","INS: " };
// 			for (int i = 0; i < DataLocs.size(); ++i) {
// 				char * cstr;
// 				TecUtilVarGetName(i+1, &cstr);
// 				string str = cstr;
// 				TecUtilStringDealloc(&cstr);
// 				for (auto const & s : IntVarPrefixes) {
// 					if (str.length() >= s.length()){
// 						int j = str.compare(0, s.length(), s);
// 						if (j == 0){
// 							DataLocs[i] = ValueLocation_CellCentered;
// 							break;
// 						}
// 					}
// 				}
// 			}
// 			SurfZoneNum = Sphere.SaveAsTriFEZone(NucleusName != "" ? NucleusName + " Sphere" : "Sphere (" + CPString + ")", DataTypes, DataLocs, XYZVarNums);
// 
// 			Set TmpSet(SurfZoneNum);
// 			TecUtilZoneSetMesh(SV_SHOW, TmpSet.getRef(), 0.0, FALSE);
// 			TecUtilZoneSetContour(SV_SHOW, TmpSet.getRef(), 0.0, TRUE);
// 			StyleValue styleValue;
// 			styleValue.set((Boolean_t)1, SV_FIELDLAYERS, SV_SHOWCONTOUR);
// 
// // 			if (IsOk){
// // 				ArgList_pa ArgList = TecUtilArgListAlloc();
// // 				TecUtilArgListAppendString(ArgList, SV_NAME, string("Sphere (" + CPString + ")").c_str());
// // 				TecUtilArgListAppendInt(ArgList, SV_ZONETYPE, ZoneType_FETriangle);
// // 				TecUtilArgListAppendInt(ArgList, SV_IMAX, NumPoints);
// // 				TecUtilArgListAppendInt(ArgList, SV_JMAX, NumTriangles);
// // 				TecUtilArgListAppendArray(ArgList, SV_VALUELOCATION, VarLocations.data());
// // 				TecUtilArgListAppendArray(ArgList, SV_VARDATATYPE, VarDataTypes.data());
// // 				IsOk = TecUtilDataSetAddZoneX(ArgList);
// // 				TecUtilArgListDealloc(&ArgList);
// // 			}
// // 
// // // 			IsOk = TecUtilDataSetAddZone("Sphere", NumPoints, NumTriangles, 0, ZoneType_FETriangle, nullptr);
// // 			EntIndex_t SurfZoneNum = TecUtilDataSetGetNumZones();
// // 			FieldData_pa SurfXYZPtrs[3];
// // 
// // 			for (int i = 0; i < 3 && IsOk; ++i){
// // 				SurfXYZPtrs[i] = TecUtilDataValueGetWritableRef(SurfZoneNum, XYZVarNums[i]);
// // 				IsOk = VALID_REF(SurfXYZPtrs[i]);
// // 			}
// // 
// // 			if (IsOk){
// // 				NodeMap_pa NodeMap;
// // 				NodeMap = TecUtilDataNodeGetWritableRef(SurfZoneNum);
// // 				IsOk = VALID_REF(NodeMap);
// // 
// // 
// // 				TecUtilDataLoadBegin();
// // 				for (LgIndex_t NodeNum = 0; NodeNum < NumPoints; ++NodeNum){
// // 					for (int i = 0; i < 3; ++i){
// // 						TecUtilDataValueSetByRef(SurfXYZPtrs[i], NodeNum + 1, p[NodeNum][i]);
// // 					}
// // 				}
// // 				for (LgIndex_t TriNum = 0; TriNum < NumTriangles; ++TriNum){
// // 					for (int i = 0; i < 3; ++i){
// // 						TecUtilDataNodeSetByRef(NodeMap, TriNum + 1, i + 1, t[TriNum][i] + 1);
// // 					}
// // 				}
// // 				TecUtilDataLoadEnd();
// // 
// // 			}
// 			if (IsOk){
// 				Set_pa TmpSet = TecUtilSetAlloc(FALSE);
// 				TecUtilSetAddMember(TmpSet, SurfZoneNum, FALSE);
// 
// 				ArgList_pa CurrentArgList = TecUtilArgListAlloc();
// 				TecUtilArgListAppendString(CurrentArgList, SV_P1, SV_FIELDMAP);
// 				TecUtilArgListAppendString(CurrentArgList, SV_P2, SV_EFFECTS);
// 				TecUtilArgListAppendString(CurrentArgList, SV_P3, SV_USETRANSLUCENCY);
// 				TecUtilArgListAppendSet(CurrentArgList, SV_OBJECTSET, TmpSet);
// 				TecUtilArgListAppendArbParam(CurrentArgList, SV_IVALUE, FALSE);
// 				TecUtilStyleSetLowLevelX(CurrentArgList);
// 				TecUtilArgListDealloc(&CurrentArgList);
// 
// 				TecUtilStateChanged(StateChange_ZonesAdded, (ArbParam_t)TmpSet);
// 				TecUtilZoneSetScatter(SV_SHOW, TmpSet, 0.0, FALSE);
// 				TecUtilZoneSetMesh(SV_SHOW, TmpSet, 0.0, TRUE);
// 				TecUtilZoneSetActive(TmpSet, AssignOp_PlusEquals);
// 
// 				TecUtilSetDealloc(&TmpSet);
// 			}
// 
// 			/*
// 			*	Save node and triangle numbers to the fe volume zone's aux data
// 			*/
// 
// 			if (IsOk){
// 				AuxDataZoneSetItem(SurfZoneNum, CSMAuxData.GBA.SphereCPName, CPString);
// 				AuxDataZoneSetItem(SurfZoneNum, CSMAuxData.GBA.ZoneType, CSMAuxData.GBA.ZoneTypeSphereZone);
// 				AuxDataZoneSetItem(SurfZoneNum, CSMAuxData.GBA.SphereCPNum, to_string(CPNum));
// 				string MovedPointNumsStr = "";
// 				for (auto const & i : MovedPointNums){
// 					MovedPointNumsStr += to_string(i+1) + ","; // +1 to switch from base-0 to base-1 so the node nums match with resulting zone
// 				}
// 				AuxDataZoneSetItem(SurfZoneNum, CSMAuxData.GBA.SphereConstrainedNodeNums, MovedPointNumsStr);
// 				AuxDataZoneSetItem(SurfZoneNum, CSMAuxData.GBA.SphereConstrainedNodeIntersectCPTypes, MovedPointCPTypes);
// 				AuxDataZoneSetItem(SurfZoneNum, CSMAuxData.GBA.SphereConstrainedNodeIntersectCPNames, MovedPointCPNames);
// 				AuxDataZoneSetItem(SurfZoneNum, CSMAuxData.GBA.SphereConstrainedNodeIntersectCPTotalOffsetNames, MovedPointCPNamesTotalCount);
// 				AuxDataZoneSetItem(SurfZoneNum, CSMAuxData.GBA.NumGBs, to_string(NumTriangles));
// 				AuxDataZoneSetItem(SurfZoneNum, CSMAuxData.GBA.GPsPerGB, to_string(3 + 3 * NumEdgeGPs));
// 				AuxDataZoneSetItem(SurfZoneNum, CSMAuxData.GBA.PointsPerGP, to_string(NumSTPoints));
// 				AuxDataZoneSetItem(SurfZoneNum, CSMAuxData.GBA.SourceZoneNum, to_string(CPZoneNum));
// 				AuxDataZoneSetItem(SurfZoneNum, CSMAuxData.GBA.SourceNucleusName, NucleusName);
// 				AuxDataZoneSetItem(Sphere.GetZoneNum(), CSMAuxData.GBA.SphereRadius, to_string(Radius));
// 				AuxDataZoneSetItem(Sphere.GetZoneNum(), CSMAuxData.GBA.SphereRadiusCPDistRatio, to_string(UserRadius));
// 			}
// 		}
// 
// 		/*
// 		*	Get a list of all unique nodes that neighbor nodes that
// 		*	were moved to GP intersections
// 		*/
// // 		TecUtilDataLoadBegin();
// 		vector<vector<int> > ConstrainedNeighborNodesNum(MovedPointNums.size());
// 		LgIndex_t NumNeighborNodes = 0;
// 		for (int i = 0; i < MovedPointNums.size(); ++i){
// 			ConstrainedNeighborNodesNum[i].reserve(6);
// 			for (int j = 0; j < NumTriangles; ++j){
// 				for (int k = 0; k < 3; ++k){
// 					if (t[j][k] == MovedPointNums[i]){
// 						for (int ii = 0; ii < 3; ++ii){
// 							if (t[j][ii] != MovedPointNums[i]){
// 								Boolean_t IsFound = FALSE;
// 								for (int jj = 0; jj < ConstrainedNeighborNodesNum[i].size() && !IsFound; ++jj){
// 									if (t[j][ii] == ConstrainedNeighborNodesNum[i][jj]){
// 										IsFound = TRUE;
// 									}
// 								}
// 								if (!IsFound)
// 									ConstrainedNeighborNodesNum[i].push_back(t[j][ii]);
// 							}
// 						}
// 					}
// 				}
// 			}
// 			NumNeighborNodes += (int)ConstrainedNeighborNodesNum[i].size();
// 			for (int j = 0; j < IntZoneSaddleCPNodeNums.size(); ++j){
// 				if (i == IntZoneSaddleCPNodeNums[j][0]){
// 					IntZoneSaddleCPNodeNums[j].push_back(MovedPointNums[i]);
// 				}
// 			}
// 			for (int j = 0; j < AllIntCPNodeNums.size(); ++j){
// 				if (i == AllIntCPNodeNums[j][0]){
// 					AllIntCPNodeNums[j].push_back(MovedPointNums[i]);
// 				}
// 			}
// 		}
// // 		TecUtilDataLoadEnd();
// 
// 		/*
// 		 *	Need to find the nodes closest to any SGP-sphere intersections for
// 		 *	non-saddle CPs
// 		 */
// 		for (int i = 0; i < AllIntCPNodeNums.size(); ++i){
// 			if (AllIntCPNodeNums[i].size() < 3){
// 				vec3 Pt1, Pt2;
// 				for (int j = 0; j < 3; ++j){
// 					Pt1[j] = TecUtilDataValueGetByZoneVar(CPZoneNum, XYZVarNums[j], AllIntCPNodeNums[i][1]);
// 				}
// 				double MinDist = 1e200;
// 				int MinInd = -1;
// 				for (int j = 0; j < NumPoints; ++j){
// 					for (int k = 0; k < 3; ++k)
// 						Pt2[k] = p[j][k];
// 					double TmpDist = DistSqr(Pt1, Pt2);
// 					if (TmpDist < MinDist){
// 						MinDist = TmpDist;
// 						MinInd = j;
// 					}
// 				}
// 				if (MinInd >= 0)
// 					AllIntCPNodeNums[i].push_back(MinInd);
// 			}
// 		}
// 
// 
// 
// 		/*
// 		 *	Get all necessary raw pointers for GradPath creation.
// 		 */
// 
// 		EntIndex_t GradXYZVarNums[3];
// 
// 		if (IsOk){
// 			vector<string> TmpStrs = {
// 				"X Density Gradient",
// 				"Y Density Gradient",
// 				"Z Density Gradient"
// 			};
// 			for (int i = 0; i < 3 && IsOk; ++i){
// 				GradXYZVarNums[i] = VarNumByName(TmpStrs[i]);
// 				IsOk = (GradXYZVarNums[i] > 0);
// 			}
// 			if (!IsOk){
// 				StatusDrop(AddOnID);
// 				TecUtilDialogErrMsg("Couldn't find gradient vector variables.");
// 				TecUtilDataLoadEnd();
// 				TecUtilLockFinish(AddOnID);
// 				return;
// 			}
// 		}
// 
// // 		if (IsOk){
// // 
// // 			// Enable system volume zone
// // 
// // 			Set_pa TmpSet = TecUtilSetAlloc(FALSE);
// // 			TecUtilSetAddMember(TmpSet, 1, FALSE);
// // 			TecUtilZoneSetActive(TmpSet, AssignOp_PlusEquals);
// // 			TecUtilSetDealloc(&TmpSet);
// // 
// // 		}
// 
// 		vector<FieldDataPointer_c> GradRawPtrs(3);
// 		FieldDataPointer_c RhoRawPtr;
// 
// // 		Data load begin tells Tecplot that I don't want it to reorganize memory until I call
// // 		Data load end, so that the raw pointers fetched below remain valid.
// 		TecUtilDataLoadBegin();
// 
// 		for (int i = 0; i < 4 && IsOk; ++i){
// 			if (i == 0){
// 				IsOk = RhoRawPtr.InitializeReadPtr(VolZoneNum, CutoffVarNum);
// 			}
// 			else{
// 				IsOk = GradRawPtrs[i - 1].InitializeReadPtr(VolZoneNum, GradXYZVarNums[i - 1]);
// 			}
// 		}
// 
// 		EntIndex_t NumZonesBeforeVolumes = TecUtilDataSetGetNumZones();
// 
// 		MemoryRequired = (sizeof(Boolean_t) * (2 * NumPoints + 2 * NumTriangles) + sizeof(EntIndex_t) * NumTriangles) / 1024;
// 		TecUtilMemoryChangeNotify(MemoryRequired);
// 
// 		vector<Boolean_t> NodeHasSaddleCP(NumPoints, FALSE),
// 			NodeIsConstrained(NumPoints, FALSE),
// 			TriangleHasSaddleCP(NumTriangles, FALSE),
// 			TriangleHasConstrainedNode(NumTriangles, FALSE);
// 		for (int i = 0; i < IntZoneSaddleCPNodeNums.size(); ++i){
// // 			NodeHasSaddleCP[IntZoneSaddleCPNodeNums[i][3]] = TRUE;
// 			NodeHasSaddleCP[IntZoneSaddleCPNodeNums[i][3]] = MovedCPIsBond[i];
// 		}
// 		for (int i = 0; i < MovedPointNums.size(); ++i){
// 			NodeIsConstrained[MovedPointNums[i]] = TRUE;
// 		}
// 
// 
// 		TecUtilMemoryChangeNotify((NumSTPoints * 4 * NumPoints * sizeof(double)) / 1024);
// 
// 		/*
// 		*	Seed grad paths for all non-saddle-cp nodes
// 		*	and then for edges.
// 		*/
// 
// 
// 		vector<GradPath_c> GPsNonSaddle(NumPoints), GPsEdges(NumEdges * NumEdgeGPs);
// 		GPTerminate_e HowTerminate;
// 		if (UseCutoff)
// 			HowTerminate = GPTerminate_AtRhoValue;
// 		else
// 			HowTerminate = GPTerminate_AtBoundary;
// 
// 		StreamDir_e StreamDir;
// 		if (CPType == -3)
// 			StreamDir = StreamDir_Reverse;
// 		else
// 			StreamDir = StreamDir_Forward;
// 
// 
// 		StatusDrop(AddOnID);
// 
// 
// 		if (IsOk){
// 			TmpString = ProgressStr.str() + string(" ... Step 2 of 4 ... Seeding gradient paths");
// 			StatusLaunch(TmpString, AddOnID, TRUE);
// 		}
// 
// 		Boolean_t UserQuit = FALSE;
// 
// 		int TmpNumIterations = NumPoints + NumEdges * NumEdgeGPs;
// 		int NumCompleted = 0;
// 
// 		high_resolution_clock::time_point StatusStartTime = high_resolution_clock::now();
// 
// // #ifndef _DEBUG
// #pragma omp parallel for schedule(dynamic)
// // #endif
// 		for (int PtNum = 0; PtNum < NumPoints; ++PtNum){
// 			if (omp_get_thread_num() == 0){
// 				UserQuit = !StatusUpdate(NumCompleted, TmpNumIterations, TmpString, AddOnID, StatusStartTime);
// #pragma omp flush (UserQuit)
// 			}
// #pragma omp flush (UserQuit)
// #pragma omp atomic
// 			NumCompleted++;
// 
// 			if (!NodeHasSaddleCP[PtNum] && !UserQuit){
// 				vec3 NodePos;
// 				for (int ii = 0; ii < 3; ++ii)
// 					NodePos[ii] = p[PtNum][ii];
// 
// // 				IsOk = GPsNonSaddle[PtNum].SetupGradPath(NodePos,
// // 					StreamDir, 
// // 					NumSTPoints,
// // 					HowTerminate, 
// // 					nullptr, vector<FieldDataPointer_c>(), nullptr, nullptr,
// // 					&CutoffVal, 
// // 					VolMaxIJK,
// // 					VolMaxXYZ,
// // 					VolMinXYZ,
// // 					GradRawPtrs, 
// // 					RhoRawPtr);
// 
// 				IsOk = GPsNonSaddle[PtNum].SetupGradPath(NodePos, 
// 					StreamDir, 
// 					NumSTPoints, 
// 					GPType_Classic, 
// 					HowTerminate, 
// 					nullptr, nullptr, nullptr, 
// 					&CutoffVal, 
// 					VolInfo, 
// 					vector<FieldDataPointer_c>(), 
// 					GradRawPtrs, 
// 					RhoRawPtr);
// 
// 				if (IsOk)
// 					IsOk = GPsNonSaddle[PtNum].Seed();
// 
// 				if (IsOk && DistSqr(NodePos, GPsNonSaddle[PtNum].XYZAt(0)) > DistSqr(NodePos, GPsNonSaddle[PtNum].XYZAt(GPsNonSaddle[PtNum].GetCount() - 1))){
// 					IsOk = GPsNonSaddle[PtNum].Reverse();
// 				}
// 			}
// 		}
// 
// 		vector<vector<int> > ConstrainedNeighborEdgeNodesNum = ConstrainedNeighborNodesNum;
// 		vector<vector<int> > ConstrainedNeighborhoodEdgeNums(ConstrainedNeighborNodesNum.size());
// 
// // #ifndef _DEBUG
// #pragma omp parallel for schedule(dynamic)
// // #endif
// 		for (int EdgeNum = 0; EdgeNum < NumEdges; ++EdgeNum){
// 			if (omp_get_thread_num() == 0){
// 				UserQuit = !StatusUpdate(NumCompleted, TmpNumIterations, TmpString, AddOnID, StatusStartTime);
// #pragma omp flush (UserQuit)
// 			}
// #pragma omp flush (UserQuit)
// #pragma omp atomic
// 			NumCompleted++;
// 			// Get the first edge node and a vector to step down the edge
// 			vec3 eNodes[2], DelVec, SeedPt;
// 			for (int i = 0; i < 2; ++i){
// 				for (int j = 0; j < 3; ++j){
// 					eNodes[i][j] = p[e[EdgeNum][i]][j];
// 				}
// 			}
// 			DelVec = (eNodes[1] - eNodes[0]) / static_cast<double>(NumEdgeGPs + 1);
// 			for (int EdgeGPNum = 0; EdgeGPNum < NumEdgeGPs; ++EdgeGPNum){
// 				if (!UserQuit){
// 					SeedPt = eNodes[0] + DelVec * static_cast<double>(EdgeGPNum + 1);
// 					int GPInd = EdgeNum * NumEdgeGPs + EdgeGPNum;
// 
// // 					IsOk = GPsEdges[GPInd].SetupGradPath(SeedPt,
// // 						StreamDir,
// // 						NumSTPoints,
// // 						HowTerminate,
// // 						nullptr, vector<FieldDataPointer_c>(), nullptr, nullptr,
// // 						&CutoffVal,
// // 						VolMaxIJK,
// // 						VolMaxXYZ,
// // 						VolMinXYZ,
// // 						GradRawPtrs,
// // 						RhoRawPtr);
// 
// 					IsOk = GPsEdges[GPInd].SetupGradPath(SeedPt,
// 						StreamDir,
// 						NumSTPoints,
// 						GPType_Classic,
// 						HowTerminate,
// 						nullptr, nullptr, nullptr,
// 						&CutoffVal,
// 						VolInfo,
// 						vector<FieldDataPointer_c>(),
// 						GradRawPtrs,
// 						RhoRawPtr);
// 
// 					if (IsOk)
// 						IsOk = GPsEdges[GPInd].Seed();
// 
// 					if (IsOk && DistSqr(SeedPt, GPsEdges[GPInd].XYZAt(0)) > DistSqr(SeedPt, GPsEdges[GPInd].XYZAt(GPsEdges[GPInd].GetCount() - 1))){
// 						IsOk = GPsEdges[GPInd].Reverse();
// 					}
// 				}
// 			}
// 
// 			// Also update the ConstrainedNeighborEdgeNodesNum list so that the saddle GPs know their
// 			// new neighbors.
// 			for (int m = 0; m < MovedPointNums.size(); ++m){
// 				int ei;
// 				for (ei = 0; ei < 2; ++ei){
// 					if (MovedPointNums[m] == e[EdgeNum][ei]){
// 						break;
// 					}
// 				}
// 				if (ei < 2){
// 					for (int & n : ConstrainedNeighborEdgeNodesNum[m]){
// 						if (e[EdgeNum][(ei + 1) % 2] == n){
// 							if (ei == 0){
// 								n = EdgeNum * NumEdgeGPs;
// 							}
// 							else{
// 								n = EdgeNum * NumEdgeGPs + NumEdgeGPs - 1;
// 							}
// 							break;
// 						}
// 					}
// 					break;
// 				}
// 			}
// 		}
// 
// 		if (UserQuit){
// 			StatusDrop(AddOnID);
// 
// 			MemoryRequired = sizeof(point) * NumPoints
// 				+ sizeof(triangle) * NumTriangles
// 				+ sizeof(int) * 2 * NumEdges;
// 			TecUtilMemoryChangeNotify(-MemoryRequired / 1024);
// 			delete p, t;
// 			for (int i = 0; i < NumEdges; ++i){
// 				delete[] e[i];
// 			}
// 			delete e;
// 
// 			TecUtilMemoryChangeNotify((NumSTPoints * 4 * NumPoints * sizeof(double)) / 1024);
// 			GPsNonSaddle.clear();
// 
// 			TecUtilDataLoadEnd();
// 			TecUtilLockFinish(AddOnID);
// 			return;
// 		}
// 
// 
// 		/*
// 		 *	Seed grad paths for saddle-cp nodes
// 		 */
// 		TecUtilMemoryChangeNotify((NumSTPoints * 4 * NumPoints * sizeof(double)) / 1024);
// 		TecUtilMemoryChangeNotify((NumSTPoints * 4 * IntZoneSaddleCPNodeNums.size() * sizeof(double)) / 1024);
// 		vector<vector<GradPath_c> > GPsSaddle(IntZoneSaddleCPNodeNums.size(),vector<GradPath_c>(6)), GPsSaddleEdges(IntZoneSaddleCPNodeNums.size(), vector<GradPath_c>(6*NumEdgeGPs));
// 		vector<vector<int> > GPsSaddleClosestPtNum(IntZoneSaddleCPNodeNums.size(), vector<int>(6)), GPsNonSaddleClosestPtNums(IntZoneSaddleCPNodeNums.size(), vector<int>(6));
// 		ConstrainedNeighborhoodEdgeNums.reserve(6);
// 		int NumIntZoneSaddleCPNodeNums = IntZoneSaddleCPNodeNums.size();
// 		
// #ifdef _DEBUG
// 		vector<vector<int> > EdgeList(NumEdges, vector<int>(2));
// #pragma omp parallel for
// 		for (int i = 0; i < NumEdges; ++i){
// 			for (int j = 0; j < 2; ++j){
// 				EdgeList[i][j] = e[i][j];
// 			}
// 		}
// #endif // _DEBUG
// 
// 
// // #ifndef _DEBUG
//  #pragma omp parallel for schedule(dynamic)
// // #endif // !_DEBUG
// 		for (int i = 0; i < NumIntZoneSaddleCPNodeNums; ++i){
// 			/*
// 			 *	Get the point's index in the list of moved points so
// 			 *	neighbor information can be looked up.
// 			 *	Its index in MovedPointNums is the same as its index
// 			 *	in ConstrainedNeighborNodesNum.
// 			 */
// 			Boolean_t IsFound = FALSE;
// 			for (int j = 0; j < MovedPointNums.size() && !IsFound; ++j){
// 				if (IntZoneSaddleCPNodeNums[i][3] == MovedPointNums[j]){
// 					IsFound = TRUE;
// 
// 					//GPsSaddle[i].resize(ConstrainedNeighborNodesNum[j].size());
// 					//GPsSaddleEdges[i].resize(ConstrainedNeighborNodesNum[j].size() * NumEdgeGPs);
// 					//ConstrainedNeighborhoodEdgeNums.reserve(ConstrainedNeighborNodesNum[j].size());
// 					//GPsSaddleClosestPtNum[i].resize(ConstrainedNeighborNodesNum[j].size());
// 					//GPsNonSaddleClosestPtNums[i].resize(ConstrainedNeighborNodesNum[j].size());
// 
// 					vec3 NodePos;
// 					for (int ii = 0; ii < 3; ++ii)
// 						NodePos[ii] = p[MovedPointNums[j]][ii];
// 
// 					vec3 VolCPPos = IntCPPos[i];
// 
// 					vector<vec3> SeedPts(ConstrainedNeighborNodesNum[j].size());
// 
// 					for (int k = 0; k < ConstrainedNeighborNodesNum[j].size() && IsOk; ++k){
// 						/*
// 						*	For each neighboring node, need to get an orthogonal set
// 						*	of unit vectors to describe the positions of the neighbor
// 						*	node's streamtraces's closest point to the volume CP
// 						*	relative to the volume CP.
// 						*/
// 						/*
// 						 *	Get the closest point between the saddle CP and
// 						 *	the neighbor node's grad path.
// 						 */
// // 						vec3 ClosestPoint = GPsNonSaddle[ConstrainedNeighborNodesNum[j][k]].ClosestPoint(VolCPPos, GPsNonSaddleClosestPtNums[i][k]);
// 						vec3 ClosestPoint = GPsEdges[ConstrainedNeighborEdgeNodesNum[j][k]].ClosestPoint(VolCPPos, GPsNonSaddleClosestPtNums[i][k]);
// 						/*
// 						*	v1 is the direction from the main node to the volume CP.
// 						*	v2 is the direction from main node to neighbor node, and is
// 						*	orthogonal to v1.
// 						*	v3 is orthogonal to v1 and v2.
// 						*/
// 
// 						vec3 v1 = VolCPPos - NodePos,
// 							v2 = ClosestPoint - VolCPPos,
// 							v3 = cross(v1, v2);
// 						v1 = normalise(v1);
// 						v3 = normalise(v3);
// 						v2 = normalise(cross(v3, v1));
// 
// 						SeedPts[k] = v1 * (Radius * 0.1);
// 						double SeedTheta = PI / 2.0;
// 						vec4 TmpVec4 = join_cols(SeedPts[k], ones<vec>(1));
// 						TmpVec4 = RotationMatrix(SeedTheta, v3) * TmpVec4;
// 
// 						SeedPts[k] = vec3(TmpVec4.subvec(0, 2)) + VolCPPos;
// 
// // 						IsOk = GPsSaddle[i][k].SetupGradPath(SeedPts[k], StreamDir, NumSTPoints, HowTerminate,
// // 							nullptr, vector<FieldDataPointer_c>(), nullptr, nullptr, &CutoffVal,
// // 							VolMaxIJK, VolMaxXYZ, VolMinXYZ,
// // 							GradRawPtrs,
// // 							RhoRawPtr);
// 
// 						IsOk = GPsSaddle[i][k].SetupGradPath(SeedPts[k],
// 							StreamDir,
// 							NumSTPoints,
// 							GPType_Classic,
// 							HowTerminate,
// 							nullptr, nullptr, nullptr,
// 							&CutoffVal,
// 							VolInfo,
// 							vector<FieldDataPointer_c>(),
// 							GradRawPtrs,
// 							RhoRawPtr);
// 
// 						if (IsOk){
// 							IsOk = GPsSaddle[i][k].Seed();
// // 							GPsSaddle[i][k].SaveAsOrderedZone("Saddle GP " + CPName + " Node " + to_string(ConstrainedNeighborNodesNum[j][k]), Green_C);
// 						}
// 
// 						if (IsOk){
// 							GPsSaddle[i][k].ConcatenateResample(IntGPs[i], NumSTPoints, GPsSaddleClosestPtNum[i][k]);
// // 							GPsSaddle[i][k].SaveAsOrderedZone("Full Saddle GP " + CPName + " Node " + to_string(ConstrainedNeighborNodesNum[j][k]), Purple_C);
// 						}
// 
// 						if (IsOk && DistSqr(NodePos, GPsSaddle[i][k].XYZAt(0)) > DistSqr(NodePos, GPsSaddle[i][k].XYZAt(GPsSaddle[i][k].GetCount() - 1))){
// 							IsOk = GPsSaddle[i][k].Reverse();
// 						}
// // 						else{
// // 							GPsSaddleClosestPtNum[i][k] = GPsSaddle[i][k].GetCount() - GPsSaddleClosestPtNum[i][k] - 1;
// // 						}
// 
// 
// 					}
// 
// 					/*
// 					*	Now make edge GPs between this saddle GP and its clockwise neighboring saddle GP
// 					*/
// 					for (int k = 0; k < ConstrainedNeighborNodesNum[j].size() && IsOk; ++k){
// 						// always make saddle edge GPs in direction of first to second edge node
// 						// (where the edge is the edge opposite the sphere saddle node)
// 						// 
// 
// 						int StartNodeNum, EndNodeNum;
// 						vector<int> n;
// 						for (int kk = k + 1; kk < ConstrainedNeighborNodesNum[j].size(); ++kk){
// 							bool NodeFound = false;
// 							n = vector<int>({ ConstrainedNeighborNodesNum[j][k], ConstrainedNeighborNodesNum[j][kk] });
// 							for (int EdgeNum = 0; EdgeNum < NodeEdgeNums[n[0]].size() && !NodeFound; ++EdgeNum){
// 								if (std::find(ConstrainedNeighborhoodEdgeNums[j].begin(), ConstrainedNeighborhoodEdgeNums[j].end(), NodeEdgeNums[n[0]][EdgeNum]) == ConstrainedNeighborhoodEdgeNums[j].end()){
// 									for (int ei = 0; ei < 2; ++ei){
// 										if (e[NodeEdgeNums[n[0]][EdgeNum]][ei] == n[0] && e[NodeEdgeNums[n[0]][EdgeNum]][(ei + 1) % 2] == n[1]) {
// 											NodeFound = true;
// 											ConstrainedNeighborhoodEdgeNums[j].push_back(NodeEdgeNums[n[0]][EdgeNum]);
// 											if (ei == 0){
// 												n = vector<int>({ k, kk });
// 											}
// 											else{
// 												n = vector<int>({ kk, k });
// 											}
// 											break;
// 										}
// 									}
// 								}
// 							}
// 							if (NodeFound){
// 								vec3 RotVec = SeedPts[n[0]] - VolCPPos, v1 = normalise(RotVec), v2 = normalise(SeedPts[n[1]] - VolCPPos);
// 								vec3 RotAxis = normalise(cross(v1, v2));
// 								vec4 TmpVec4a = join_cols(RotVec, ones<vec>(1));
// 								double Alpha = acos(dot(v1, v2));
// 								double DelAlpha = Alpha / static_cast<double>(NumEdgeGPs + 1);
// 
// 								for (int ei = 0; ei < NumEdgeGPs; ++ei){
// 									vec4 TmpVec4b = RotationMatrix(DelAlpha * static_cast<double>(ei + 1), RotAxis) * TmpVec4a;
// 									vec3 SeedPt = vec3(TmpVec4b.subvec(0, 2)) + VolCPPos;
// 
// 									int GPInd = (ConstrainedNeighborhoodEdgeNums[j].size() - 1) * NumEdgeGPs + ei;
// 
// // 									IsOk = GPsSaddleEdges[i][GPInd].SetupGradPath(SeedPt, StreamDir, NumSTPoints, HowTerminate,
// // 										nullptr, vector<FieldDataPointer_c>(), nullptr, nullptr, &CutoffVal,
// // 										VolMaxIJK, VolMaxXYZ, VolMinXYZ,
// // 										GradRawPtrs,
// // 										RhoRawPtr);
// 
// 									IsOk = GPsSaddleEdges[i][GPInd].SetupGradPath(SeedPt,
// 										StreamDir,
// 										NumSTPoints,
// 										GPType_Classic,
// 										HowTerminate,
// 										nullptr, nullptr, nullptr,
// 										&CutoffVal,
// 										VolInfo,
// 										vector<FieldDataPointer_c>(),
// 										GradRawPtrs,
// 										RhoRawPtr);
// 
// 									if (IsOk){
// 										IsOk = GPsSaddleEdges[i][GPInd].Seed();
// 										// 							GPsSaddle[i][k].SaveAsOrderedZone("Saddle GP " + CPName + " Node " + to_string(ConstrainedNeighborNodesNum[j][k]), Green_C);
// 									}
// 
// 									if (IsOk){
// 										GPsSaddleEdges[i][GPInd].ConcatenateResample(IntGPs[i], NumSTPoints, GPsSaddleClosestPtNum[i][k]);
// 										// 							GPsSaddle[i][k].SaveAsOrderedZone("Full Saddle GP " + CPName + " Node " + to_string(ConstrainedNeighborNodesNum[j][k]), Purple_C);
// 									}
// 
// 									if (IsOk && DistSqr(NodePos, GPsSaddleEdges[i][GPInd].XYZAt(0)) > DistSqr(NodePos, GPsSaddleEdges[i][GPInd].XYZAt(GPsSaddleEdges[i][GPInd].GetCount() - 1))){
// 										IsOk = GPsSaddleEdges[i][GPInd].Reverse();
// 									}
// 								}
// 							}
// 						}
// 					}
// 
// 
// 					break;
// 				}
// 			}
// 
// 			if (!IsFound){
// 				TecUtilDialogErrMsg("Couldn't find point in list");
// 			}
// 		}
// 		TecUtilDataLoadEnd();
// 
// // #pragma omp parallel for
// // 		for (int i = 0; i < GPsNonSaddle.size(); ++i){
// // 			if (GPsNonSaddle[i].IsMade()){
// // 				GPsNonSaddle[i].PointPrepend(CPPos, 0.0);
// // 			}
// // 		}
// // 
// // #pragma omp parallel for
// // 		for (int i = 0; i < GPsSaddle.size(); ++i){
// // 			for (int j = 0; j < GPsSaddle[i].size(); ++j){
// // 				if (GPsSaddle[i][j].IsMade()){
// // 					GPsSaddle[i][j].PointPrepend(CPPos, 0.0);
// // 				}
// // 			}
// // 		}
// // 
// // 
// 
// 		if (SaveGPs){
// 			Set ZoneSet;
// 			for (int i = 0; i < GPsNonSaddle.size(); ++i){
// 				if (GPsNonSaddle[i].IsMade()){
// 					GPsNonSaddle[i].SaveAsOrderedZone("GP " + CPName + " " + "Node " + to_string(i + 1));
// 					ZoneSet += GPsNonSaddle[i].GetZoneNum();
// 					AuxDataZoneSetItem(GPsNonSaddle[i].GetZoneNum(), GBASphereCPName, CPString);
// 					AuxDataZoneSetItem(GPsNonSaddle[i].GetZoneNum(), GBASphereCPNum, to_string(CPNum));
// 					AuxDataZoneSetItem(GPsNonSaddle[i].GetZoneNum(), GBAGPNodeNum, to_string(i + 1));
// 					AuxDataZoneSetItem(GPsNonSaddle[i].GetZoneNum(), GBAZoneType, GBAZoneTypeGradPath);
// 				}
// 			}
// 
// 			for (int i = 0; i < GPsSaddle.size(); ++i){
// 				for (int j = 0; j < GPsSaddle[i].size(); ++j){
// 					if (GPsSaddle[i][j].IsMade()){
// 						GPsSaddle[i][j].SaveAsOrderedZone("GP " + CPName + " " + "Node " + to_string(MovedPointNums[i] + 1) + "." + to_string(j + 1));
// 						ZoneSet += GPsSaddle[i][j].GetZoneNum();
// 						AuxDataZoneSetItem(GPsSaddle[i][j].GetZoneNum(), GBASphereCPName, CPString);
// 						AuxDataZoneSetItem(GPsSaddle[i][j].GetZoneNum(), GBASphereCPNum, to_string(CPNum));
// 						AuxDataZoneSetItem(GPsSaddle[i][j].GetZoneNum(), GBAGPNodeNum, to_string(MovedPointNums[i] + 1));
// 						AuxDataZoneSetItem(GPsSaddle[i][j].GetZoneNum(), GBAZoneType, GBAZoneTypeGradPath);
// 					}
// 				}
// 			}
// 
// 			for (int i = 0; i < GPsEdges.size(); ++i){
// 				if (GPsEdges[i].IsMade()){
// 					GPsEdges[i].SaveAsOrderedZone("GP " + CPName + " " + "EdgeGP " + to_string(i + 1), Blue_C);
// 					ZoneSet += GPsEdges[i].GetZoneNum();
// 					AuxDataZoneSetItem(GPsEdges[i].GetZoneNum(), GBASphereCPName, CPString);
// 					AuxDataZoneSetItem(GPsEdges[i].GetZoneNum(), GBASphereCPNum, to_string(CPNum));
// 					AuxDataZoneSetItem(GPsEdges[i].GetZoneNum(), GBAGPNodeNum, to_string(i + 1));
// 					AuxDataZoneSetItem(GPsEdges[i].GetZoneNum(), GBAZoneType, GBAZoneTypeGradPath);
// 				}
// 			}
// 
// 			for (int i = 0; i < GPsSaddleEdges.size(); ++i){
// 				for (int j = 0; j < GPsSaddleEdges[i].size(); ++j){
// 					if (GPsSaddleEdges[i][j].IsMade()){
// 						GPsSaddleEdges[i][j].SaveAsOrderedZone("GP " + CPName + " " + "Saddle Edge " + to_string(MovedPointNums[i] + 1) + "." + to_string(j + 1));
// 						ZoneSet += GPsSaddleEdges[i][j].GetZoneNum();
// 						AuxDataZoneSetItem(GPsSaddleEdges[i][j].GetZoneNum(), GBASphereCPName, CPString);
// 						AuxDataZoneSetItem(GPsSaddleEdges[i][j].GetZoneNum(), GBASphereCPNum, to_string(CPNum));
// 						AuxDataZoneSetItem(GPsSaddleEdges[i][j].GetZoneNum(), GBAGPNodeNum, to_string(MovedPointNums[i] + 1));
// 						AuxDataZoneSetItem(GPsSaddleEdges[i][j].GetZoneNum(), GBAZoneType, GBAZoneTypeGradPath);
// 					}
// 				}
// 			}
// 			TecUtilZoneSetActive(ZoneSet.getRef(), AssignOp_MinusEquals);
// 		}
// 
// 
// 		/*
// 		 *	Now have all gradient paths.
// 		 *	The concatenated gradient paths are indexed identically
// 		 *	to the SaddleNeighborNodesNum list.
// 		 *	
// 		 *	Now find the sphere edges whose grad paths straddle
// 		 *	an IB.
// 		 */
// 
// 
// // 		TecUtilMemoryChangeNotify((NumEdges * (sizeof(Boolean_t) + sizeof(vec3))) / 1024);
// // 
// // 		vector<Boolean_t> EdgeStraddlesIB(NumEdges, FALSE);
// // 		vector<vec3> EdgeMidpoints(NumEdges);
// // 
// // 
// // 		TecUtilDataLoadBegin();
// // #pragma omp parallel for schedule(dynamic)
// // 		for (int EdgeNum = 0; EdgeNum < NumEdges; ++EdgeNum){
// // 			if (NodeIsConstrained[e[EdgeNum][0]] || NodeIsConstrained[e[EdgeNum][1]]){
// // 				/*
// // 				 *	WTF right?
// // 				 *	Well, if one of the edge's nodes is a saddle-point node then
// // 				 *	this finds the GP that was made for that node to connect to 
// // 				 *	the edge's other node.
// // 				 */
// // 				for (int i = 0; i < 2; ++i){
// // 					if (NodeIsConstrained[e[EdgeNum][i]]){
// // 						for (int j = 0; j < MovedPointNums.size(); ++j){
// // 							if (e[EdgeNum][i] == MovedPointNums[j]){
// // 								for (int k = 0; k < ConstrainedNeighborNodesNum[j].size(); ++k){
// // 									if (e[EdgeNum][(i + 1) % 2] == ConstrainedNeighborNodesNum[j][k]){
// // 										for (int ii = 0; ii < IntZoneSaddleCPNodeNums.size(); ++ii){
// // 											if (IntZoneSaddleCPNodeNums[ii][3] == MovedPointNums[j]){
// // 												EdgeStraddlesIB[EdgeNum] = GPsStraddleIB(GPsNonSaddle[e[EdgeNum][(i + 1) % 2]], GPsSaddle[ii][k], IBCheckAngle, IBCheckDistRatio);
// // 											}
// // 										}
// // 									}
// // 								}
// // 							}
// // 						}
// // 						break;
// // 					}
// // 				}
// // 			}
// // 			else
// // 			{
// // 				EdgeStraddlesIB[EdgeNum] = GPsStraddleIB(GPsNonSaddle[e[EdgeNum][0]], GPsNonSaddle[e[EdgeNum][1]], IBCheckAngle, IBCheckDistRatio);
// // 			}
// // 
// // 			if (EdgeStraddlesIB[EdgeNum]){
// // 				vec3 Pts[2];
// // 
// // 				for (int i = 0; i < 2; ++i)
// // 					for (int j = 0; j < 3; ++j)
// // 						Pts[i][j] = p[e[EdgeNum][i]][j];
// // 
// // 				EdgeMidpoints[EdgeNum] = (Pts[0] + Pts[1]) * 0.5;
// // 			}
// // 		}
// // 
// // 		int IBEdgeCount = 0;
// // 		for (int i = 0; i < NumEdges; ++i)
// // 			IBEdgeCount += (int)EdgeStraddlesIB[i];
// // 
// // 		if (IsOk && IBEdgeCount > 0){
// // 			IsOk = TecUtilDataSetAddZone("Inter-IB Edges", IBEdgeCount, 1, 1, ZoneType_Ordered, nullptr);
// // 			EntIndex_t IBZoneNum;
// // 			FieldData_pa IBXYZFDPtrs[3];
// // 			if (IsOk){
// // 				IBZoneNum = TecUtilDataSetGetNumZones();
// // 				for (int i = 0; i < 3 && IsOk; ++i){
// // 					IBXYZFDPtrs[i] = TecUtilDataValueGetWritableRef(IBZoneNum, XYZVarNums[i]);
// // 					IsOk = VALID_REF(IBXYZFDPtrs[i]);
// // 				}
// // 			}
// // 			if (IsOk){
// // 				IBEdgeCount = 1;
// // 				for (int i = 0; i < NumEdges; ++i){
// // 					if (EdgeStraddlesIB[i]){
// // 
// // 						for (int j = 0; j < 3; ++j){
// // 							TecUtilDataValueSetByRef(IBXYZFDPtrs[j], IBEdgeCount, EdgeMidpoints[i][j]);
// // 						}
// // 						IBEdgeCount++;
// // 					}
// // 				}
// // 				Set_pa IBSet = TecUtilSetAlloc(FALSE);
// // 				TecUtilSetAddMember(IBSet, IBZoneNum, FALSE);
// // 
// // 				TecUtilZoneSetScatter(SV_SHOW, IBSet, 0.0, TRUE);
// // 				TecUtilZoneSetContour(SV_SHOW, IBSet, 0.0, FALSE);
// // 				TecUtilZoneSetMesh(SV_SHOW, IBSet, 0.0, FALSE);
// // 
// // 				TecUtilZoneSetScatter(SV_COLOR, IBSet, nullptr, Black_C);
// // 				TecUtilZoneSetScatter(SV_FRAMESIZE, IBSet, 0.5, nullptr);
// // 				TecUtilZoneSetScatterSymbolShape(SV_GEOMSHAPE, IBSet, GeomShape_Sphere);
// // 
// // 				TecUtilZoneSetActive(IBSet, AssignOp_PlusEquals);
// // 
// // 				TecUtilSetDealloc(&IBSet);
// // 			}
// // 			if (IsOk){
// // 				IsOk = AuxDataZoneSetItem(IBZoneNum, GBASphereCPName, CPName);
// // 				if (IsOk)
// // 					IsOk = AuxDataZoneSetItem(IBZoneNum, GBAZoneType, GBAZoneTypeIBEdgeZone);
// // 			}
// // 		}
// // 		TecUtilDataLoadEnd();
// // 
// // 
// // 		EdgeMidpoints.clear();
// // 		EdgeStraddlesIB.clear();
// // 		TecUtilMemoryChangeNotify((NumEdges * (sizeof(Boolean_t) + sizeof(vec3))) / 1024);
// 
// 		/*
// 		 *	Now put together the grad paths into fe volume
// 		 *	objects.
// 		 */
// 
// // 		int NumSaddleFEZones = 0;
// // 		for (int i = 0; i < GPsSaddle.size(); ++i)
// // 			NumSaddleFEZones += static_cast<int>(GPsSaddle[i].size());
// 
// 		TecUtilMemoryChangeNotify((NumSTPoints * 4 * 4 * NumTriangles * sizeof(double)) / 1024);
// 		vector<FESurface_c> FEVolumes(NumTriangles);
// 		vector<bool> IsDivergentGB(NumTriangles, false);
// 
// 		TmpString = ProgressStr.str() + string(" ... Step 3 of 4 ... Integrating volumes");
// 		TmpNumIterations = NumTriangles;
// // #ifndef _DEBUG
// // 			TmpNumIterations /= numCPU;
// // #endif // !_DEBUG
// 
// 		vector<vector<double> > IntVals(NumTriangles, vector<double>(IntVarNumList.size() + 1, 0));
// 			
// 		NumCompleted = 0;
// 
// 		vector<vector<int> > GPSaddleNodeNums1(NumTriangles), GPSaddleNodeNums2(NumTriangles), GPNonSaddleNodeNums(NumTriangles);
// 		for (int i = 0; i < NumTriangles; ++i){
// 			GPSaddleNodeNums1[i].reserve(3);
// 			GPSaddleNodeNums2[i].reserve(3);
// 			GPNonSaddleNodeNums[i].reserve(3);
// 		}
// 
// // 		TecUtilDialogMessageBox("Before making FEVolume", MessageBoxType_Information);
// 
// 		vector<VolExtentIndexWeights_s> ThreadVolInfo(omp_get_num_procs(), VolInfo);
// 
// 		StatusStartTime = high_resolution_clock::now();
// 
// 		TecUtilDataLoadBegin();
// // #ifndef _DEBUG
// #pragma omp parallel for schedule(dynamic)
// // #endif
// 		for (int TriNum = 0; TriNum < NumTriangles; ++TriNum){
// // 		for (int TriNum = 24; TriNum < 25; ++TriNum){
// 			if (omp_get_thread_num() == 0){
// 				UserQuit = !StatusUpdate(NumCompleted, TmpNumIterations, TmpString, AddOnID, StatusStartTime);
// #pragma omp flush (UserQuit)
// 			}
// #pragma omp flush (UserQuit)
// 			if (!UserQuit){
// 				vector<GradPath_c*> GPPtrs;
// 				if (TriangleHasSaddleCP[TriNum]){
// 					GPPtrs.reserve(4 + 4 * NumEdgeGPs);
// 					/*
// 					 *	Need to find which GPs to use
// 					 *	for the FE volume.
// 					 */
// 					int Count = 0;
// 					int NonConstrainedNodeNums[2];
// 					int ConstrainedGPNums[2] = { -1, -1 };
// 					int ConstrainedNodeNum;
// 					int SaddleNum;
// 					int ConstrainedCornerNum;
// 					for (int i = 0; i < 3; ++i){
// 						if (NodeIsConstrained[t[TriNum][i]]){
// 							ConstrainedNodeNum = t[TriNum][i];
// 							ConstrainedCornerNum = i;
// 							for (int j = 0; j < 2; ++j)
// 								NonConstrainedNodeNums[j] = t[TriNum][(i + j + 1) % 3];
// 						}
// 					}
// 
// 					Count = 0;
// 					Boolean_t IsFound = FALSE;
// 
// 					for (int i = 0; i < IntZoneSaddleCPNodeNums.size(); ++i){
// 						if (IntZoneSaddleCPNodeNums[i][3] == ConstrainedNodeNum){
// 							SaddleNum = i;
// 							for (int j = 0; j < MovedPointNums.size(); ++j){
// 								if (IntZoneSaddleCPNodeNums[i][3] == MovedPointNums[j]){
// 									for (int k = 0; k < ConstrainedNeighborNodesNum[j].size(); ++k){
// 										for (int ii = 0; ii < 2; ++ii){
// 											if (ConstrainedNeighborNodesNum[j][k] == NonConstrainedNodeNums[ii]){
// 												ConstrainedGPNums[ii] = k;
// 												Count++;
// 												break;
// 											}
// 										}
// 										if (Count >= 2) IsFound = TRUE;
// 									}
// 									break;
// 								}
// 							}
// 							break;
// 						}
// 					}
// 					if (IsFound){
// 						vector<int> Nodes({
// 							NonConstrainedNodeNums[0],
// 							NonConstrainedNodeNums[1],
// 							ConstrainedNodeNum
// 						});
// 						for (int n = 0; n < Nodes.size(); ++n){
// 							if (n < Nodes.size() - 1){
// 								GPPtrs.push_back(&GPsNonSaddle[Nodes[n]]);
// 							}
// 							else{
// 								if (GPsSaddle[SaddleNum][ConstrainedGPNums[1]].IsMade()) {
// 									GPPtrs.push_back(&GPsSaddle[SaddleNum][ConstrainedGPNums[1]]);
// 									bool EdgeFound = false;
// 									for (int ei = 0; ei < ConstrainedNeighborhoodEdgeNums[SaddleNum].size() && !EdgeFound; ++ei) {
// 										if (std::find(TriangleEdgeNums[TriNum].begin(), TriangleEdgeNums[TriNum].end(), ConstrainedNeighborhoodEdgeNums[SaddleNum][ei]) != TriangleEdgeNums[TriNum].end()) {
// 											if (e[ConstrainedNeighborhoodEdgeNums[SaddleNum][ei]][0] == NonConstrainedNodeNums[1]) {
// 												for (int i = 0; i < NumEdgeGPs; ++i) {
// 													if (GPsSaddleEdges[SaddleNum][ei * NumEdgeGPs + i].IsMade())
// 														GPPtrs.push_back(&GPsSaddleEdges[SaddleNum][ei * NumEdgeGPs + i]);
// 												}
// 											}
// 											else {
// 												for (int i = NumEdgeGPs - 1; i >= 0; --i) {
// 													if (GPsSaddleEdges[SaddleNum][ei * NumEdgeGPs + i].IsMade())
// 														GPPtrs.push_back(&GPsSaddleEdges[SaddleNum][ei * NumEdgeGPs + i]);
// 												}
// 											}
// 										}
// 									}
// 									GPPtrs.push_back(&GPsSaddle[SaddleNum][ConstrainedGPNums[0]]);
// 								}
// 							}
// 							bool EdgeFound = false;
// 							for (int ei = 0; ei < 3 && !EdgeFound; ++ei){
// 								for (int i = 0; i < 2; ++i){
// 									if (e[TriangleEdgeNums[TriNum][ei]][i] == Nodes[n] && e[TriangleEdgeNums[TriNum][ei]][(i + 1) % 2] == Nodes[(n + 1) % Nodes.size()]) {
// 										if (i == 0){
// 											for (int j = 0; j < NumEdgeGPs; ++j){
// 												GPPtrs.push_back(&GPsEdges[TriangleEdgeNums[TriNum][ei] * NumEdgeGPs + j]);
// 											}
// 										}
// 										else{
// 											for (int j = NumEdgeGPs - 1; j >= 0; --j){
// 												GPPtrs.push_back(&GPsEdges[TriangleEdgeNums[TriNum][ei] * NumEdgeGPs + j]);
// 											}
// 										}
// 										EdgeFound = true;
// 										break;
// 									}
// 								}
// 							}
// 						}
// // 						GPPtrs = {&GPsNonSaddle[NonConstrainedNodeNums[0]],
// // 							&GPsNonSaddle[NonConstrainedNodeNums[1]],
// // 							&GPsSaddle[ConstrainedNodeNum][ConstrainedGPNums[1]],
// // 							&GPsSaddle[ConstrainedNodeNum][ConstrainedGPNums[0]]};
// 							
// 						/*IsOk = FEVolumes[TriNum].MakeGradientBundle(GPsNonSaddle[NonConstrainedNodeNums[0]],
// 							GPsNonSaddle[NonConstrainedNodeNums[1]],
// 							GPsSaddle[ConstrainedNodeNum][ConstrainedGPNums[1]],
// 							GPsSaddle[ConstrainedNodeNum][ConstrainedGPNums[0]]);*/
// 						if (IsOk){
// 							for (int i = 0; i < 2; ++i){
// 								GPNonSaddleNodeNums[TriNum].push_back(NonConstrainedNodeNums[i]);
// 								GPSaddleNodeNums1[TriNum].push_back(ConstrainedNodeNum);
// 								GPSaddleNodeNums2[TriNum].push_back(ConstrainedGPNums[i]);
// 							}
// 						}
// 					}
// 					else{
// 						TecUtilDialogErrMsg("Couldn't find GPs to make FE volume");
// 					}
// 				}
// 				else{
// 					GPPtrs.reserve(3 + 3 * NumEdgeGPs);
// 					for (int n = 0; n < 3; ++n){
// 						GPPtrs.push_back(&GPsNonSaddle[t[TriNum][n]]);
// 						bool EdgeFound = false;
// 						for (int ei = 0; ei < 3 && !EdgeFound; ++ei){
// 							for (int i = 0; i < 2; ++i){
// 								if (e[TriangleEdgeNums[TriNum][ei]][i] == t[TriNum][n] && e[TriangleEdgeNums[TriNum][ei]][(i + 1) % 2] == t[TriNum][(n + 1) % 3]) {
// 									if (i == 0){
// 										for (int j = 0; j < NumEdgeGPs; ++j){
// 											GPPtrs.push_back(&GPsEdges[TriangleEdgeNums[TriNum][ei] * NumEdgeGPs + j]);
// 										}
// 									}
// 									else{
// 										for (int j = NumEdgeGPs - 1; j >= 0; --j){
// 											GPPtrs.push_back(&GPsEdges[TriangleEdgeNums[TriNum][ei] * NumEdgeGPs + j]);
// 										}
// 									}
// 									EdgeFound = true;
// 									break;
// 								}
// 							}
// 						}
// 					}
// // 					GPPtrs = {&GPsNonSaddle[t[TriNum][0]],
// // 						&GPsNonSaddle[t[TriNum][1]],
// // 						&GPsNonSaddle[t[TriNum][2]]};
// 					/*IsOk = FEVolumes[TriNum].MakeGradientBundle(GPsNonSaddle[t[TriNum][0]],
// 						GPsNonSaddle[t[TriNum][1]],
// 						GPsNonSaddle[t[TriNum][2]]);*/
// 					if (IsOk){
// 						for (int i = 0; i < 3; ++i){
// 							GPNonSaddleNodeNums[TriNum].push_back(t[TriNum][i]);
// 						}
// 					}
// 				}
// 				if (IsOk){
// 					/*
// 					 *	As a workaround until ring surface divergent GBs are handled correctly,
// 					 *	Use the old integration method for divergent GBs, which are determined by
// 					 *	checking for a dot product of ~-1 between the last line segments of two
// 					 *	GPs.
// 					 */
// 					vector<vec3*> TerminalVecs(GPPtrs.size(), nullptr);
// 					for (int i = 0; i < GPPtrs.size() - 1 && !IsDivergentGB[TriNum]; ++i){
// 						if (TerminalVecs[i] == nullptr){
// 							TerminalVecs[i] = new vec3; 
// 							*TerminalVecs[i] = normalise(GPPtrs[i]->XYZAt(-1) - GPPtrs[i]->XYZAt(-2));
// 						}
// 						for (int j = i + 1; j < GPPtrs.size() && !IsDivergentGB[TriNum]; ++j){
// 							if (TerminalVecs[j] == nullptr){
// 								TerminalVecs[j] = new vec3;
// 								*TerminalVecs[j] = normalise(GPPtrs[j]->XYZAt(-1) - GPPtrs[j]->XYZAt(-2));
// 							}
// 							IsDivergentGB[TriNum] = abs(dot(*TerminalVecs[i], *TerminalVecs[j]) < DivergentGPMaxTerminalDotProduct);
// 						}
// 					}
// 					for (auto * i : TerminalVecs) if (i != nullptr) delete i;
// 					if (SaveGBs || IsDivergentGB[TriNum]){
// 						// 					if (TriNum < 3) TecUtilDialogMessageBox("Before making FEVolume", MessageBoxType_Information);
// 						IsOk = FEVolumes[TriNum].MakeFromGPs(GPPtrs, true, true);
// 						// 					if (TriNum < 3) TecUtilDialogMessageBox("After making FEVolume", MessageBoxType_Information);
// 					}
// 					if (!IsDivergentGB[TriNum] && DoIntegration)
// 						IntegrateUsingIsosurfaces(GPPtrs, IntResolution, ThreadVolInfo[omp_get_thread_num()], IntVarPtrs, IntVals[TriNum], TriangleHasSaddleCP[TriNum]);
// 				}
// #pragma omp atomic
// 				NumCompleted++;
// 			}
// 		}
// 
// 
// 		/*
// 		 *	Now collect and save integration values
// 		 */
// 		NumZones = TecUtilDataSetGetNumZones();
// 		vector<string> NewVarNames;
// 		vector<int> NewVarNums;
// 		vector<FieldDataType_e> ZoneDataTypes(NumZones, FieldDataType_Bit);
// 
// 		if (DoIntegration) {
// 
// 			for (int i = 1; i <= NumZones; ++i) {
// 				if (AuxDataZoneHasItem(i, CSMAuxData.GBA.ZoneType)) {
// 					ZoneDataTypes[i - 1] = FieldDataType_Double;
// 				}
// 			}
// 			vector<ValueLocation_e> ZoneDataLocs(NumZones, ValueLocation_CellCentered);
// 			// 	vector<ValueLocation_e> ZoneDataLocs(NumZones, ValueLocation_Nodal);
// 			vector<string> VarNameAppends = { "I", "N", "S" };
// 			for (int i = 0; i < IntVarNameList.size() && IsOk; ++i) {
// 				for (int j = 0; j < 3; ++j) {
// 					string TmpStr;
// 					for (int k = 0; k <= j; ++k)
// 						TmpStr += VarNameAppends[k];
// 					NewVarNames.push_back(TmpStr + ": " + IntVarNameList[i]);
// 					NewVarNums.push_back(VarNumByName(NewVarNames.back()));
// 					if (NewVarNums[i * 3 + j] < 0) {
// 						ArgList_pa ArgList = TecUtilArgListAlloc();
// 
// 						IsOk = TecUtilArgListAppendString(ArgList, SV_NAME, NewVarNames.back().c_str());
// 						if (IsOk)
// 							IsOk = TecUtilArgListAppendArray(ArgList, SV_VARDATATYPE, ZoneDataTypes.data());
// 						if (IsOk)
// 							IsOk = TecUtilArgListAppendArray(ArgList, SV_VALUELOCATION, ZoneDataLocs.data());
// 
// 						if (IsOk)
// 							IsOk = TecUtilDataSetAddVarX(ArgList);
// 
// 						if (IsOk)
// 							NewVarNums.back() = TecUtilDataSetGetNumVars();
// 
// 						if (IsOk) {
// 							Set_pa TmpSet = TecUtilSetAlloc(FALSE);
// 							IsOk = TecUtilSetAddMember(TmpSet, NewVarNums.back(), FALSE);
// 							if (IsOk)
// 								TecUtilStateChanged(StateChange_VarsAdded, (ArbParam_t)TmpSet);
// 							TecUtilSetDealloc(&TmpSet);
// 						}
// 
// 						TecUtilArgListDealloc(&ArgList);
// 					}
// 				}
// 			}
// 			// 			TecUtilDialogMessageBox((VectorToString(NewVarNames) + ": " + VectorToString(NewVarNums)).c_str(), MessageBoxType_Information);
// 			if (IsOk) {
// 				for (int j = 0; j < 3; ++j) {
// 					string TmpStr;
// 					for (int k = 0; k <= j; ++k)
// 						TmpStr += VarNameAppends[k];
// 					TmpStr += ": Volume";
// 					if (std::find(NewVarNames.begin(), NewVarNames.end(), TmpStr) == NewVarNames.end()) {
// 						NewVarNames.push_back(TmpStr);
// 						NewVarNums.push_back(VarNumByName(NewVarNames.back()));
// 						if (NewVarNums.back() < 0) {
// 							ArgList_pa ArgList = TecUtilArgListAlloc();
// 
// 							IsOk = TecUtilArgListAppendString(ArgList, SV_NAME, NewVarNames.back().c_str());
// 							if (IsOk)
// 								IsOk = TecUtilArgListAppendArray(ArgList, SV_VARDATATYPE, ZoneDataTypes.data());
// 							if (IsOk)
// 								IsOk = TecUtilArgListAppendArray(ArgList, SV_VALUELOCATION, ZoneDataLocs.data());
// 
// 							if (IsOk)
// 								IsOk = TecUtilDataSetAddVarX(ArgList);
// 
// 							if (IsOk)
// 								NewVarNums.back() = TecUtilDataSetGetNumVars();
// 
// 							if (IsOk) {
// 								Set_pa TmpSet = TecUtilSetAlloc(FALSE);
// 								IsOk = TecUtilSetAddMember(TmpSet, NewVarNums.back(), FALSE);
// 								if (IsOk)
// 									TecUtilStateChanged(StateChange_VarsAdded, (ArbParam_t)TmpSet);
// 								TecUtilSetDealloc(&TmpSet);
// 							}
// 
// 							TecUtilArgListDealloc(&ArgList);
// 						}
// 					}
// 				}
// 			}
// 
// 			// 			TecUtilDialogMessageBox((VectorToString(NewVarNames) + ": " + VectorToString(NewVarNums)).c_str(), MessageBoxType_Information);
// 			vector<string> AuxDataNames(IntVarNameList.size());
// 			for (int i = 0; i < IntVarNameList.size(); ++i)
// 				AuxDataNames[i] = RemoveStringChar(IntVarNameList[i], " ");
// 			AuxDataNames.push_back(string("Volume"));
// 
// 			Sphere = FESurface_c(Sphere.GetZoneNum(), VolZoneNum, XYZVarNums, IntVarNumList, true);
// 
// 
// 			/*
// 			 *	Save any divergent GBs as zones in order to integrate them, then integrate them.
// 			 *	Also integrate the sphere here, but only if divergent GBs were found.
// 			 */
// 			bool DivergentGBsExist = false;
// 			for (int i = 0; i < NumTriangles && !DivergentGBsExist; ++i) {
// 				DivergentGBsExist = IsDivergentGB[i];
// 			}
// 			if (DivergentGBsExist) {
// 				TecUtilDataLoadBegin();
// 				vector<FESurface_c> TempGBs(NumTriangles);
// 				int NumVolumes = 0;
// 				Set TmpSet;
// 				for (int i = 0; i < NumTriangles; ++i) {
// 					if (IsDivergentGB[i] && FEVolumes[i].IsMade()) {
// 						DivergentGBsExist = true;
// 						int TmpZoneNum = FEVolumes[i].SaveAsTriFEZone(XYZVarNums, "Temp GB" + to_string(i + 1));
// 						TempGBs[i] = FESurface_c(TmpZoneNum, VolZoneNum, XYZVarNums, IntVarNumList);
// 						TmpSet += TmpZoneNum;
// 						NumVolumes++;
// 					}
// 				}
// 				NumVolumes++; // +1 for the Sphere integration
// 
// 				TecUtilZoneSetActive(TmpSet.getRef(), AssignOp_PlusEquals);
// 				TecUtilRedrawAll(TRUE);
// 
// 				/*
// 				 *	Now integrate all the saved zones
// 				 */
// 				int NumComplete = 0;
// 				StatusStartTime = high_resolution_clock::now();
// 				// #ifndef _DEBUG
// #pragma omp parallel for schedule(dynamic)
// // #endif
// 				for (int i = 0; i < NumTriangles; ++i) {
// 					if (omp_get_thread_num() == 0)
// 					{
// 						stringstream ss;
// 						ss << "Integrating other volumes " << NumComplete + 1;
// 						if (AtomNameList.size() > 1)
// 							ss << " of " << NumVolumes;
// 						UserQuit = !StatusUpdate(NumComplete, NumVolumes, ss.str(), AddOnID, StatusStartTime);
// #pragma omp flush (UserQuit)
// 					}
// #pragma omp flush (UserQuit)
// 					if (!UserQuit){
// 						if (TempGBs[i].IsMade()) {
// 							TempGBs[i].DoIntegrationNew(IntResolution - 1, TRUE);
// 							IntVals[i] = TempGBs[i].GetIntResults();
// #pragma omp atomic
// 							NumComplete++;
// 						}
// 						if (!Sphere.IntResultsReady() && omp_get_thread_num() == 0) {
// 							Sphere.DoIntegrationNew(IntResolution - 1, TRUE);
// #pragma omp atomic
// 							NumComplete++;
// 						}
// 					}
// 
// 				}
// 
// 
// 				/*
// 				 *	Now delete the zones
// 				 */
// 				Set DelZones;
// 				for (auto const & i : TempGBs) {
// 					if (i.IsMade())
// 						DelZones += i.GetZoneNum();
// 				}
// 				TecUtilDataSetDeleteZone(DelZones.getRef());
// 
// 
// 
// 				TecUtilDataLoadEnd();
// 			}
// 			else {
// 				stringstream ss;
// 				ss << "Integrating sphere";
// 				TecUtilPleaseWait(ss.str().c_str(), TRUE);
// 				Sphere.DoIntegrationNew(IntResolution - 1, TRUE);
// 				TecUtilPleaseWait(ss.str().c_str(), FALSE);
// 			}
// 
// 
// 
// 			vector<double> SphereTriangleAreas;
// 			vector<vector<double> > SphereElemIntVals = Sphere.GetTriSphereIntValsByElem(&SphereTriangleAreas);
// 			for (int i = 0; i < NumTriangles; ++i) {
// 				for (int j = 0; j < SphereElemIntVals[i].size(); ++j) {
// 					SphereElemIntVals[i][j] += IntVals[i][j];
// 				}
// 			}
// 			vector<vector<double> > NormalizedValues = SphereElemIntVals;
// 			vector<double> TotalList(SphereElemIntVals[0].size(), 0.0),
// 				TotalNormlizedList(SphereElemIntVals[0].size(), 0.0),
// 				IntScaleFactors(SphereElemIntVals[0].size());
// 			for (int i = 0; i < SphereElemIntVals.size(); ++i) {
// 				for (int j = 0; j < SphereElemIntVals[i].size(); ++j) {
// 					NormalizedValues[i][j] /= SphereTriangleAreas[i];
// 					TotalNormlizedList[j] += NormalizedValues[i][j];
// 					TotalList[j] += SphereElemIntVals[i][j];
// 				}
// 			}
// 			for (int i = 0; i < SphereElemIntVals[0].size(); ++i) {
// 				IntScaleFactors[i] = TotalList[i] / TotalNormlizedList[i];
// 			}
// 			for (int i = 0; i < SphereElemIntVals.size(); ++i) {
// 				vector<double> TmpVec = SphereElemIntVals[i];
// 				SphereElemIntVals[i] = vector<double>();
// 				SphereElemIntVals[i].reserve(3 * TmpVec.size());
// 				for (int j = 0; j < TmpVec.size(); ++j) {
// 					SphereElemIntVals[i].push_back(TmpVec[j]);
// 					SphereElemIntVals[i].push_back(NormalizedValues[i][j]);
// 					SphereElemIntVals[i].push_back(NormalizedValues[i][j] * IntScaleFactors[j]);
// 				}
// 			}
// 
// 			if (std::find(IntVarNameList.begin(), IntVarNameList.end(), "Volume") == IntVarNameList.end())
// 				IntVarNameList.push_back("Volume");
// 			AuxDataZoneSetItem(Sphere.GetZoneNum(), CSMAuxData.GBA.AtomicBasinIntegrationVariables, VectorToString(IntVarNameList, ","));
// 			AuxDataZoneSetItem(Sphere.GetZoneNum(), CSMAuxData.GBA.AtomicBasinIntegrationValues, VectorToString(TotalList, ","));
// 
// 
// 			TecUtilDataLoadBegin();
// 
// 			// 		if (SphereElemIntVals[0].size() != NewVarNums.size()) {
// 			// 			TecUtilDialogErrMsg("Discrepency between number of integrated variables to be saved and the number of variable numbers to use when saving.");
// 			// 		}
// 
// 			vector<FieldDataPointer_c> Ptrs(SphereElemIntVals[0].size());
// 			for (int j = 0; j < Ptrs.size(); ++j) {
// 				Ptrs[j].InitializeWritePtr(Sphere.GetZoneNum(), NewVarNums[j]);
// 			}
// 
// 			// Write cell-centered integration values to sphere and FE zones
// #pragma omp parallel for
// 			for (int e = 0; e < NumTriangles; ++e) {
// 				for (int j = 0; j < SphereElemIntVals[e].size(); ++j) {
// 					Ptrs[j].Write(e, SphereElemIntVals[e][j]);
// 				}
// 			}
// 			for (auto & i : Ptrs)
// 				i.Close();
// 
// 
// 
// 			TecUtilDataLoadEnd();
// 		}
// 
// 		if (UserQuit){
// 			StatusDrop(AddOnID);
// 
// 			MemoryRequired = sizeof(point) * NumPoints
// 				+ sizeof(triangle) * NumTriangles
// 				+ sizeof(int) * 2 * NumEdges;
// 			TecUtilMemoryChangeNotify(-MemoryRequired / 1024);
// 			delete p, t;
// 			for (int i = 0; i < NumEdges; ++i){
// 				delete[] e[i];
// 			}
// 			delete e;
// 
// 			TecUtilMemoryChangeNotify(-int(NumSTPoints * 4 * (NumPoints + IntZoneSaddleCPNodeNums.size()) * sizeof(double)) / 1024);
// 			GPsNonSaddle.clear();
// 			GPsSaddle.clear();
// 
// 			TecUtilMemoryChangeNotify((NumSTPoints * 4 * 4 * NumTriangles * sizeof(double)) / 1024);
// 			FEVolumes.clear();
// 
// 			TecUtilDataLoadEnd();
// 			TecUtilLockFinish(AddOnID);
// 			return;
// 		}
// 		TecUtilDataLoadEnd();
// 
// 
// 		TecUtilMemoryChangeNotify(-int(NumSTPoints * 4 * (NumPoints + IntZoneSaddleCPNodeNums.size()) * sizeof(double)) / 1024);
// 		GPsNonSaddle.clear();
// 		GPsSaddle.clear();
// 
// 		/*
// 		 *	Test volume integration of all FE zones
// 		 */
// 
// // 		double IntSum = 0;
// // 
// // 		for (auto const & i : FEVolumes){
// // 			IntSum += i.IntVolume(20, CPPos);
// // 		}
// // 
// // 		TecUtilDialogMessageBox(string("Volume is " + to_string(IntSum)).c_str(), MessageBoxType_Information);
// 
// 		/*
// 		 *	Now have structure of all zones in memory,
// 		 *	so make all FE volume zones.
// 		 */
// 
// 
// 		if (SaveGBs){
// 			TmpString = ProgressStr.str() + string(" ... Step 4 of 4 ... Saving volumes");
// 			TecUtilDataLoadBegin();
// 			vector<EntIndex_t> TriangleFEZoneNums(NumTriangles, -1);
// 			for (int TriNum = 0; TriNum < NumTriangles && IsOk; ++TriNum){
// 				if (!StatusUpdate(TriNum, NumTriangles, TmpString, AddOnID)){
// 					StatusDrop(AddOnID);
// 
// 					MemoryRequired = sizeof(point)* NumPoints
// 						+ sizeof(triangle)* NumTriangles
// 						+ sizeof(int)* 2 * NumEdges;
// 					TecUtilMemoryChangeNotify(-MemoryRequired / 1024);
// 					delete p, t;
// 					for (int i = 0; i < NumEdges; ++i){
// 						delete[] e[i];
// 					}
// 					delete e;
// 
// 					TecUtilMemoryChangeNotify((NumSTPoints * 4 * 4 * NumTriangles * sizeof(double)) / 1024);
// 					FEVolumes.clear();
// 
// 					TecUtilDataLoadEnd();
// 					TecUtilLockFinish(AddOnID);
// 					return;
// 				}
// 
// 				//TriangleFEZoneNums[TriNum] = FEVolumes[TriNum].SaveAsFEZone(VarDataTypes, XYZVarNums, CutoffVarNum);
// 
// 				// 			if (TriNum < 3) TecUtilDialogMessageBox("before saving FEVolume", MessageBoxType_Information);
// 				//TriangleFEZoneNums[TriNum] = FEVolumes[TriNum].SaveAsTriFEZone("Gradient Bundle " + to_string(TriNum + 1) + " (Zone " + to_string(TecUtilDataSetGetNumZones() + 1) + ")", VarDataTypes, VarLocations, XYZVarNums);
// 				TriangleFEZoneNums[TriNum] = FEVolumes[TriNum].SaveAsTriFEZone("Gradient Bundle " + to_string(TriNum + 1) + " (Zone " + to_string(TecUtilDataSetGetNumZones() + 1) + ")", vector<FieldDataType_e>(), vector<ValueLocation_e>(), XYZVarNums);
// 				// 			if (TriNum < 3) TecUtilDialogMessageBox("after saving FEVolume", MessageBoxType_Information);
// 
// 				if (TriangleFEZoneNums[TriNum] > 0){
// 					/*
// 					*	Save node and triangle numbers to the fe volume zone's aux data
// 					*/
// 
// 					if (IsOk)
// 						IsOk = AuxDataZoneSetItem(TriangleFEZoneNums[TriNum], CSMAuxData.GBA.SphereCPName, CPString);
// 					if (IsOk)
// 						IsOk = AuxDataZoneSetItem(TriangleFEZoneNums[TriNum], CSMAuxData.GBA.ZoneType, CSMAuxData.GBA.ZoneTypeFEVolumeZone);
// 					if (IsOk)
// 						IsOk = AuxDataZoneSetItem(TriangleFEZoneNums[TriNum], CSMAuxData.GBA.ElemNum, to_string(TriNum + 1));
// 
// 					for (int i = 0; i < 3 && IsOk; ++i)
// 						IsOk = AuxDataZoneSetItem(TriangleFEZoneNums[TriNum], CSMAuxData.GBA.NodeNums[i], to_string(t[TriNum][i] + 1));
// 
// 					if (IsOk && TriangleHasConstrainedNode[TriNum]){
// 						Boolean_t IsFound = FALSE;
// 						for (int i = 0; i < 3 && !IsFound; ++i){
// 							for (int j = 0; j < AllIntCPNodeNums.size() && !IsFound; ++j){
// 								// 							if (NodeIsConstrained[t[TriNum][i]]){
// 								if (t[TriNum][i] == AllIntCPNodeNums[j][2]){
// 									IsFound = TRUE;
// 									IsOk = AuxDataZoneSetItem(TriangleFEZoneNums[TriNum], CSMAuxData.GBA.VolumeCPName, AllIntVolumeCPNames[j]);
// 								}
// 							}
// 						}
// 					}
// 				}
// 			}
// 			TecUtilDataLoadEnd();
// 			TecUtilMemoryChangeNotify((NumSTPoints * 4 * 4 * NumTriangles * sizeof(double)) / 1024);
// 			FEVolumes.clear();
// 		}
// 
// 		
// 
// 
// 		MemoryRequired = sizeof(point) * NumPoints
// 			+ sizeof(triangle) * NumTriangles
// 			+ sizeof(int) * 2 * NumEdges;
// 		TecUtilMemoryChangeNotify(-MemoryRequired / 1024);
// 		delete p, t;
// 		for (int i = 0; i < NumEdges; ++i){
// 			delete[] e[i];
// 		}
// 		delete e;
// 
// 
// 		TmpString = ProgressStr.str() + string(" ... Unloading Data");
// 		EntIndex_t NewNumZones = TecUtilDataSetGetNumZones() - NumZones;
// 		for (EntIndex_t i = NumZones + 1; i <= NewNumZones && IsOk; ++i){
// 			TecUtilDrawGraphics(TRUE);
// 			StatusUpdate(i, NewNumZones, TmpString, AddOnID);
// 			TecUtilDrawGraphics(FALSE);
// 			for (EntIndex_t j = 0; j < 3 && IsOk; ++j){
// 				IsOk = TecUtilDataValueUnload(i, XYZVarNums[j]);
// 			}
// 		}
// 		StatusDrop(AddOnID);
// 
// 	}
// 
// 	if (DoTiming){
// 		Time2 = high_resolution_clock::now();
// 		duration<double> TimeSpan = duration_cast<duration<double>>(Time2 - Time1);
// // 		TecUtilDialogMessageBox(string(string("Runtime: ") + to_string(TimeSpan.count()) + string(" seconds.")).c_str(), MessageBoxType_Information);
// 	}
// 
// 	Set_pa VolSet = TecUtilSetAlloc(FALSE);
// 	TecUtilSetAddMember(VolSet, VolZoneNum, FALSE);
// 	TecUtilZoneSetActive(VolSet, AssignOp_MinusEquals);
// 	TecUtilSetDealloc(&VolSet);
// 
// 	if (deleteOldSphereZones)
// 		TecUtilDataSetDeleteZone(oldSphereZonesToDelete.getRef());
// 
// 	TecUtilArrayDealloc((void **)&CPNums);
// 
// 	TecGUITabSetCurrentPage(TAB1_TB_D1, 2);
// 
// // 	TecUtilDialogMessageBox("Finished making GBs", MessageBoxType_Information);
// 
// 	TecUtilLockFinish(AddOnID);
// 	return;
// }



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
		GradPtrs[i].InitializeReadPtr(VolZoneNum, GradVarNums[i]);
	}
	RhoPtr.InitializeReadPtr(VolZoneNum, RhoVarNum);

	//Vec3_c StartPoint(-6.4627935942777901, 1.5801340272202937, -0.40415061637454408);
	vec3 StartPoint;
	StartPoint << -6.77934 << -2.2642 << 0.46312;

	double TermRhoValue = 1e-3;

	GradPath_c GP(StartPoint,
		StreamDir_Reverse, 100,
		GPTerminate_AtBoundary,
		nullptr, vector<FieldDataPointer_c>(), nullptr, nullptr, &TermRhoValue,
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

void ElemWiseGradPath(int StartingElem, 
						int DirNum, // 0 for downhill, 1 for uphill
						vector<vector<int> > const & ElemConnectivity,
						vector<vector<double> > const & ElemIntVals, 
						int BasinVarNum, 
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
	for (int i : IntermediateElems) TerminalMinMaxIndices[i][DirNum] = IntermediateElems.back();

	REQUIRE(iter < ElemConnectivity.size());

	return;
}



void GetSphereBasinIntegrations(FESurface_c const & Sphere,
									vector<vector<double> > const & ElemIntVals,
									int BasinVarNum,
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
	auto * NodeConnectivityPtr = Sphere.GetNodeConnectivityListPtr();
	auto * ElemListPtr = Sphere.GetElemListPtr();
	auto * XYZListPtr = Sphere.GetXYZListPtr();

	REQUIRE(NodeConnectivityPtr != nullptr && ElemListPtr != nullptr && XYZListPtr != nullptr);

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
	int NumConvergedIter = 3;
	vector<int> NumMinMax(2,-1), OldNumMinMax(2,INT_MAX);
	MinMaxIndices.resize(2);
	vector<vector<int> > OldMinMaxIndices(2);
	vector<vector<int> > TerminalMinMaxIndices, OldTerminalMinMaxIndices;
	vector<vector<double> > SmoothElemIntVals, TmpSmoothElemIntVals;
	SmoothElemIntVals.push_back(ElemIntVals[BasinVarNum]);

	TerminalMinMaxIndices.assign(ElemConnectivity.size(), vector<int>(2, -1));
	OldTerminalMinMaxIndices = TerminalMinMaxIndices;
	MinMaxIndices.assign(2, vector<int>());

	int Iter = 0, ConvegedIter[] = { 0,0 };
	while (Iter < MaxNumSmoothing && (ConvegedIter[0] < NumConvergedIter || ConvegedIter[1] < NumConvergedIter)) {
		Iter++;
	// 		OldNumMinMax = NumMinMax;
	// 		OldMinMaxIndices = MinMaxIndices;
	// 		MinMaxIndices.assign(2, vector<int>());
	// 		OldTerminalMinMaxIndices = TerminalMinMaxIndices;
	// 		TerminalMinMaxIndices.assign(ElemConnectivity.size(), vector<int>(2, -1));
		for (int Dir = 0; Dir < 2; ++Dir) {
			if (ConvegedIter[Dir] <= 1) {
				OldNumMinMax[Dir] = NumMinMax[Dir];
				OldMinMaxIndices[Dir] = MinMaxIndices[Dir];
				for (int i = 0; i < ElemConnectivity.size(); ++i)
					OldTerminalMinMaxIndices[i][Dir] = TerminalMinMaxIndices[i][Dir];
			}
		}
		MinMaxIndices.assign(2, vector<int>());
		TerminalMinMaxIndices.assign(ElemConnectivity.size(), vector<int>(2, -1));

// #pragma omp parallel for schedule(dynamic)
		for (int i = 0; i < NumElems; ++i) {
			for (int Dir = 0; Dir < 2; ++Dir) {
				if (ConvegedIter[Dir] < NumConvergedIter) {
					ElemWiseGradPath(i, Dir, ElemConnectivity, SmoothElemIntVals, 0, TerminalMinMaxIndices);
				}
			}
		}

		/*
		 *	Now all the terminal elements have been identified.
		 *	Get a list of the unique min and max element indices.
		 */
// #pragma omp parallel for
		for (int Dir = 0; Dir < 2; ++Dir) {
			if (ConvegedIter[Dir] < NumConvergedIter) {
				for (auto const & i : TerminalMinMaxIndices) {
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
		}
		for (int Dir = 0; Dir < 2; ++Dir){
			if (NumMinMax[Dir] == OldNumMinMax[Dir]) ConvegedIter[Dir]++;
		}
		if (Iter < MaxNumSmoothing && (NumMinMax[0] != OldNumMinMax[0] || NumMinMax[1] != OldNumMinMax[1])) {
			/*
			 *	Do a simple smoothing on the ElemIntVals based on average neighborhood value
			 */
			TmpSmoothElemIntVals = SmoothElemIntVals;
#pragma omp parallel for
			for (int i = 0; i < NumElems; ++i) {
				for (int j : ElemConnectivity[i]) {
					SmoothElemIntVals[0][i] += TmpSmoothElemIntVals[0][j];
				}
				SmoothElemIntVals[0][i] /= double(ElemConnectivity[i].size() + 1);
			}

			if (NumMinMax[0] != OldNumMinMax[0])
				ConvegedIter[0] = 0;
			else if (NumMinMax[1] != OldNumMinMax[1])
				ConvegedIter[1] = 0;
		}
	}
	if (Iter > 1){
		TerminalMinMaxIndices = OldTerminalMinMaxIndices;
		MinMaxIndices = OldMinMaxIndices;
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
//  #pragma omp parallel for schedule(dynamic)
		for (int te = 0; te < MinMaxIndices[Dir].size(); ++te) {
			vector<int> NewNodeNums(XYZListPtr->size(), -1);
			for (int e = 0; e < NumElems; ++e) {
 				if (TerminalMinMaxIndices[e][Dir] == MinMaxIndices[Dir][te]) {
					for (int v = 0; v < NumIntVars; ++v) {
						BasinIntVals[Dir][te][v] += ElemIntVals[v][e];
					}
					BasinElems[Dir][te].push_back(vector<int>());
					BasinElems[Dir][te].back().reserve(3);
					for (int ei : ElemListPtr->at(e)) {
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

 	/*
 	 *	Delete any preexisting zones that would be created here
 	 *	to prevent a conflict later.
 	 */
 	Set oldZoneSet;
 	for (int z = 1; z <= TecUtilDataSetGetNumZones(); ++z) {
 		if (AuxDataZoneHasItem(z, CSMAuxData.GBA.CondensedBasinInfo)
			&& AuxDataZoneItemMatches(z, CSMAuxData.GBA.CondensedBasinDefiningVariable, BasinDefineVarName)) {
 			string tmpStr;
 			if (AuxDataZoneGetItem(z, CSMAuxData.GBA.SphereCPName, tmpStr)
 				&& std::find(SphereZoneNames.begin(), SphereZoneNames.end(), tmpStr) != SphereZoneNames.end()) {
 				oldZoneSet += z;
 			}
 		}
 	}

	/*
				*	Get all necessary raw pointers for GradPath creation.
				*/
	EntIndex_t CutoffVarNum = VarNumByName(string("Electron Density"));
	int VolZoneNum = ZoneNumByName("Full Volume");
	EntIndex_t GradXYZVarNums[3];
	vector<int> GradVarNums;
	Boolean_t IsOk = TRUE;
	if (IsOk) {
		vector<string> TmpStrs = {
			"X Density Gradient",
			"Y Density Gradient",
			"Z Density Gradient"
		};
		for (int i = 0; i < 3 && IsOk; ++i) {
			GradXYZVarNums[i] = VarNumByName(TmpStrs[i]);
			if (GradXYZVarNums[i] <= 0) {
				GradVarNums.clear();
				break;
			}
			GradVarNums.push_back(GradXYZVarNums[i]);
		}
	}

	vector<VolExtentIndexWeights_s> VolInfo(omp_get_num_procs());
	GetVolInfo(VolZoneNum, { 1,2,3 }, FALSE, VolInfo[0]);
	for (int i = 1; i < VolInfo.size(); ++i) VolInfo[i] = VolInfo[0];

	bool DeleteGradVars = false;
	if (GradVarNums.empty() && VolInfo[0].NumPts() < 600000000) {
		CalcVarsOptions_s opt;
		opt.AddOnID = AddOnID;
		opt.CalcForAllZones = FALSE;
		opt.CalcZoneNum = VolZoneNum;
		opt.RhoVarNum = CutoffVarNum;
		opt.HasGrad = (!GradVarNums.empty());
		if (opt.HasGrad) {
			opt.GradVarNums = GradVarNums;
		}
		else {
			DeleteGradVars = true;
			opt.CalcVarList = { CalcGradientVectors };
		}
		opt.IsPeriodic = FALSE;

		CalcVars(opt);

		if (GradVarNums.empty()) {
			GradVarNums.push_back(VarNumByName("X Density", true));
			for (int i = 1; i < 3; ++i) {
				GradVarNums.push_back(GradVarNums[0] + i);
				GradXYZVarNums[i - 1] = GradVarNums.back();
			}

			if (GradVarNums.empty())
				return;
		}
	}

	StatusLaunch("Finding gradient bundles...(please wait)", AddOnID);

	// Resize spheres to original size
	ResizeSpheres(1.0, TRUE, FALSE);

	high_resolution_clock::time_point startTime = high_resolution_clock::now();
	if (IntNum > 0 && SphereZoneNums.size() > 0) {
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
			vector<float> TmpVals(NumElems);
			for (int i = 0; i < NumIntVars; ++i) {
				int IntVarNum = VarNumByName(TecGUIListGetString(SLSelVar_SLST_T3_1, i + 1));
				REQUIRE(IntVarNum > 0);
				FieldData_pa CCRef = TecUtilDataValueGetReadableCCRef(SphereZoneNum, IntVarNum);
				REQUIRE(VALID_REF(CCRef));
				if (TecUtilDataValueGetType(SphereZoneNum, IntVarNum) == FieldDataType_Float) {
					TecUtilDataValueArrayGetByRef(CCRef, 1, NumElems, TmpVals.data());
					ElemIntVals[i] = vector<double>(TmpVals.begin(), TmpVals.end());
				}
				else {
					TecUtilDataValueArrayGetByRef(CCRef, 1, NumElems, ElemIntVals[i].data());
				}
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
// 				for (int & i : SphereConstrainedNodeNums) i--;// moving from base 1 to base 0
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
				for (auto const & i : SphereConstrainedNodeNames) {
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
					IsOk = RhoRawPtr.InitializeReadPtr(VolZoneNum, CutoffVarNum);
				}
				else {
					if (GradXYZVarNums[i] <= 0) {
						GradRawPtrs.clear();
						break;
					}
					IsOk = GradRawPtrs[i - 1].InitializeReadPtr(VolZoneNum, GradXYZVarNums[i - 1]);
				}
			}

			double TermRhoValue;
			vector<vector<FESurface_c> > WedgeSurfaces(2);
			
			TecGUITextFieldGetDouble(TFCutoff_TF_T1_1, &TermRhoValue);
			int NumBasins = 0;
			for (int i = 0; i < 2; ++i) {
				NumBasins += MinMaxIndices[i].size();
				WedgeSurfaces[i].resize(MinMaxIndices[i].size());
			}
			int NumGPPts;
			TecGUITextFieldGetLgIndex(TFSTPts_TF_T1_1, &NumGPPts);

			for (auto const & i : SphereConstrainedNodeCPNums) {
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
				&& AuxDataZoneGetItem(SphereZoneNum, CSMAuxData.GBA.SphereRadius, RadStr)) {
				vector<int> XYZVarNums = { 1,2,3 };
				int CPZoneNum = stoi(CPZoneNumStr);
				int CPInd = stoi(CPIndStr);
				SphereRadius = stod(RadStr);
				for (int i = 0; i < 3; ++i) {
					CPPos[i] = TecUtilDataValueGetByZoneVar(CPZoneNum, XYZVarNums[i], CPInd);
				}
			}

			// Assign each special gradient path sphere intersection to the condensed basin that most
			// "points" in the direction of the intersection point.
			// First, precompute the average vector for each basin pointing from the sphere center
			// to each of the nodes of the basin.// Store the mapping of basins to the constrained nodes to which they point.

			vector<vector<std::pair<int, double> > > CondensedBasinToConstrainedNodeNums(BasinNodes.size());
			vector<std::map<double, std::pair<int, int> > > ConstrainedNodeToCondensedBasin(SphereConstrainedNodeNums.size());

			vector<vector<vec3> > MinMaxBasinAverageDirVecs(BasinNodes.size());
			for (int i = 0; i < BasinNodes.size(); ++i) {
				MinMaxBasinAverageDirVecs[i] = vector<vec3>(BasinNodes[i].size());
				CondensedBasinToConstrainedNodeNums[i].resize(BasinNodes[i].size(), std::make_pair(-1,-1));
			}

// #pragma omp parallel for schedule(dynamic)
			for (int b = 0; b < NumBasins; ++b) {
				int i, j;
				if (b < MinMaxIndices[0].size()) {
					i = 0;
					j = b;
				}
				else {
					i = 1;
					j = b - MinMaxIndices[0].size();
				}
				MinMaxBasinAverageDirVecs[i][j] = { 0,0,0 };
// #pragma omp parallel for
				for (int ti = 0; ti < BasinElems[i][j].size(); ++ti){
					vec3 WeightedDir = zeros(3);
					for (auto const & ni : BasinElems[i][j][ti])
						WeightedDir += (BasinNodes[i][j][ni] - CPPos);
					WeightedDir *= abs(ElemIntVals[IntNum][SphereElemNums[i][j][ti]]);
// #pragma omp critical(UpdateMinMaxBasinAverageDirVecs)
					{
						MinMaxBasinAverageDirVecs[i][j] += WeightedDir;
					}
				}
// #pragma omp parallel for
// 				for (int ni = 0; ni < BasinNodes[i][j].size(); ++ni)
// 					MinMaxBasinAverageDirVecs[i][j] += (BasinNodes[i][j][ni] - CPPos);
				MinMaxBasinAverageDirVecs[i][j] = normalise(MinMaxBasinAverageDirVecs[i][j]);
			}

			// Now loop over the sphere's constrained nodes and finding the basin that 
			// points most in the direction of the constrained node.
			// Want to match each constrained node with the most closely aligned condensed basin,
			// but the best basin for a particular constrained node might be even more aligned with some other
			// constrained node, in which case we'd use a ranked system so that the basin goes with its
			// best fit constrained node, and the other constrained node will go to its next best fit.
			// To accomplish this we'll keep a sorted list for each of the constrained nodes so each one knows
			// its own ranking, and also record the best fit constrained node for each basin that will be used
			// to store the actual correspondence between basin/node.
			// 
			auto * XYZListPtr = Sphere.GetXYZListPtr();
			for (int n = 0; n < SphereConstrainedNodeNums.size(); ++n) {
				vec3 const * constrinedVec = &XYZListPtr->at(SphereConstrainedNodeNums[n]);
				for (int b = 0; b < NumBasins; ++b) {
					int i, j;
					if (b < MinMaxIndices[0].size()) {
						i = 0;
						j = b;
					}
					else {
						i = 1;
						j = b - MinMaxIndices[0].size();
					}
					if (std::find(SphereNodeNums[i][j].begin(), SphereNodeNums[i][j].end(), SphereConstrainedNodeNums[n]) != SphereNodeNums[i][j].end())
					{
						double TmpDotPdt = dot(MinMaxBasinAverageDirVecs[i][j], normalise(*constrinedVec - CPPos));
						ConstrainedNodeToCondensedBasin[n][TmpDotPdt] = std::make_pair(i, j);
						if (TmpDotPdt > CondensedBasinToConstrainedNodeNums[i][j].second) {
							CondensedBasinToConstrainedNodeNums[i][j].second = TmpDotPdt;
							CondensedBasinToConstrainedNodeNums[i][j].first = n;
						}
					}
				}
			}

			// Now CondensedBasinToConstrainedNodeNums points each condensed basin to its
			// most aligned constrained node, but it could be that two condensed basins are
			// pointing to the same node.
			// Meanwhile, ConstrainedNodeToCondensedBasin gives each constrained node a ranked
			// list of condensed basins according to best alignment.
			// We'll loop over ConstrainedNodeToCondensedBasin and for each constrained node
			// see if its most aligned condensed basin points back at it. If so, move on, if not,
			// start moving down the constrained node's ranked list until a condensed basin is found
			// that points back at the constrained node.
			// Then loop over CondensedBasinToConstrainedNodeNums to do the same check in reverse.
			// This time, if the constrained node doesn't point back to the condensed basin, then
			// it means there's another condensed basin that was better aligned with the constrained
			// node and the condensed basin shouldn't point to any constrained node.

			for (int n = 0; n < ConstrainedNodeToCondensedBasin.size(); ++n) {
				auto cn = ConstrainedNodeToCondensedBasin[n].rbegin();
				while (CondensedBasinToConstrainedNodeNums[cn->second.first][cn->second.second].first != n){
					ConstrainedNodeToCondensedBasin[n].erase(cn->first);
					if (ConstrainedNodeToCondensedBasin[n].empty()) {
// 						TecUtilDialogErrMsg("Failed to find condensed basin corresponding to SGP sphere intersection");
						break;
					}
					else
						cn = ConstrainedNodeToCondensedBasin[n].rbegin();
				}
			}

			// Now update CondensedBasinToConstrainedNodeNums
			for (int i = 0; i < CondensedBasinToConstrainedNodeNums.size(); ++i){
				for (int j = 0; j < CondensedBasinToConstrainedNodeNums[i].size(); ++j){
					auto * cb = &CondensedBasinToConstrainedNodeNums[i][j];
					if (cb->first >= 0) {
						auto cn = ConstrainedNodeToCondensedBasin[cb->first].rbegin();
						if (cn->second.first != i || cn->second.second != j) {
							cb->first = -1;
							cb->second = -1;
						}
					}
				}
			}

//  		#ifndef _DEBUG
// #pragma omp parallel for schedule(dynamic)
// #endif
			for (int b = 0; b < NumBasins; ++b) {
				if (!StatusUpdate(z, SphereZoneNums.size(), "Finding gradient bundles for " + SphereZoneNames[z] + " (" + to_string(z + 1) + " of " + to_string(SphereZoneNums.size()) + "; basin " + to_string(b + 1) + " of " + to_string(NumBasins) + ")", AddOnID, startTime)){
					break;
				}
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

				vector<vec3> SeedPts;
				vector<vector<int> > BasinPerimeterEdges;
				if (GetSortedPerimeterEdgeMidpoints(BasinElems[i][j], BasinNodes[i][j], SeedPts, BasinPerimeterEdges)) {
					/*
					 *	Prepare workspace for wedge cap
					 */
// 					vector<vector<int> > CapElems = BasinElems[i][j];
// 					vector<std::pair<int,int> > PerimeterEdgeElems(BasinPerimeterEdges.size());
// 					vector<vec3> CapNodes(BasinNodes[i][j].size());
// 					vector<bool> NodeIsPerimeter(CapNodes.size(), false);
// 					// Below map only used so we can get perimeter nodes
// 					std::map<Edge, std::pair<int,int> > EdgeToElemMap;
// 					for (int ei = 0; ei < CapElems.size(); ++ei) {
// 						for (int ej = 0; ej < 3; ++ej)
// 							EdgeToElemMap[MakeEdge(CapElems[ei][ej], CapElems[ei][(ej + 1) % 3])] = std::make_pair(ei,ej);
// 					}
// 					for (int ei = 0; ei < BasinPerimeterEdges.size(); ++ei){
// 						PerimeterEdgeElems[ei] = EdgeToElemMap[MakeEdge(BasinPerimeterEdges[ei][0], BasinPerimeterEdges[ei][1])];
// 						for (auto ej : BasinPerimeterEdges[ei])
// 							NodeIsPerimeter[ej] = true;
// 					}
// 
// 
// 					vector<GradPath_c> ThreadGPs(1);
// 					ThreadGPs[0].SetupGradPath(
// 						vec3(),
// 						StreamDir_Both,
// 						NumGPPts,
// 						GPType_Classic,
// 						GPTerminate_AtRhoValue,
// 						nullptr, nullptr, nullptr,
// 						&TermRhoValue,
// 						VolInfo[ThreadNum],
// 						vector<FieldDataPointer_c>(),
// 						GradRawPtrs,
// 						RhoRawPtr);
// 					ThreadGPs = vector<GradPath_c>(omp_get_num_procs(), ThreadGPs[0]);
// 
// 					// get cap points: terminal points of all non-perimeter basin nodes
// 					int NumBasinNodes = BasinNodes[i][j].size();
// #pragma omp parallel for schedule(dynamic)
// 					for (int ni = 0; ni < NumBasinNodes; ++ni){
// 						GradPath_c * GPPtr = &ThreadGPs[omp_get_thread_num()];
// 						GPPtr->Clear();
// 						GPPtr->SetStartPoint(BasinNodes[i][j][ni]);
// 						GPPtr->Seed(false);
// 						CapNodes[ni] = GPPtr->XYZAt(-1);
// 					}

					int NumSeedPts = SeedPts.size();
					/*
						*	Check to see if any bond path-sphere intersections appear in the perimeter.
						*	If they do, need to add them to the surface.
						*/
					std::set<int> ConstrainedNodeNums(SphereConstrainedGPNodeNums.begin(), SphereConstrainedGPNodeNums.end());
					bool ContainsConstrainedNode = false;
					for (int e = 0; e < BasinPerimeterEdges.size() && !ContainsConstrainedNode; ++e) {
						for (int ei : BasinPerimeterEdges[e]) {
							if (ConstrainedNodeNums.count(SphereNodeNums[i][j][ei]) > 0) {
								ContainsConstrainedNode = true;
								break;
							}
						}
					}
					vector<GradPath_c> BasinParimeterGPs(NumSeedPts);
					vector<GradPath_c*> GPPtrs;
					for (int p = 0; p < NumSeedPts; ++p) {
						BasinParimeterGPs[p].SetupGradPath(
							SeedPts[p],
							StreamDir_Both,
							NumGPPts,
							GPType_Classic,
							GPTerminate_AtRhoValue,
							nullptr, nullptr, nullptr,
							&TermRhoValue,
							VolInfo[0],
							vector<FieldDataPointer_c>(),
							GradRawPtrs,
							RhoRawPtr);
						GPPtrs.push_back(&BasinParimeterGPs[p]);
					}
					bool UserQuit = false;
					int NumCompleted = 0;
#ifndef _DEBUG
#pragma omp parallel for schedule(dynamic)
#endif
					for (int p = 0; p < NumSeedPts; ++p) {
						BasinParimeterGPs[p].Seed();
						BasinParimeterGPs[p].Reverse();
					}
					if (UserQuit){
						break;
					}
					if (ContainsConstrainedNode) {
						/*
							*	There are one or more constrained nodes that need to be added
							*	to the surface.
							*/
// 						auto * XYZListPtr = Sphere.GetXYZListPtr();
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
									for (auto const & seedPt : TmpSeedPts) {
										TmpGPPtrs.push_back(new GradPath_c);
										TmpGPPtrs.back()->SetupGradPath(SphereConstrainedNodeGPs[n][-1] + (seedPt * 0.5),
											StreamDir_Reverse,
											NumGPPts,
											GPType_Classic,
											GPTerminate_AtRhoValue,
											nullptr, nullptr, nullptr,
											&TermRhoValue,
											VolInfo[0],
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
						vector<GradPath_c const *> ConstGPPtrs;
						for (auto * g : GPPtrs)
							ConstGPPtrs.push_back(const_cast<GradPath_c const *>(g));
						WedgeSurfaces[i][j].MakeFromGPs(ConstGPPtrs, true);
					}
				}
			}

			/*
				*	Now save the basins and wedge surfaces as zones.
				*/

			



			Set ZoneSet;
			vector<string> MinMax = { "Min","Max" };
			for (int i = 0; i < 2; ++i) {
				for (int b = 0; b < MinMaxIndices[i].size(); ++b) {
					if (SphereRadius > 0) {
						for (int vi = 0; vi < BasinNodes[i][b].size(); ++vi) {
							BasinNodes[i][b][vi] += normalise(BasinNodes[i][b][vi] - CPPos) * 0.001;
						}
					}
					FESurface_c Basin(BasinNodes[i][b], BasinElems[i][b]);
					REQUIRE(Basin.IsMade());

// 					int SphereConstrainedNodeNum = -1;
// 					for (int n = 0; n < SphereConstrainedNodeNums.size(); ++n) {
//  						if (((i == 0 && SphereConstrainedNodeTypeStrs[n] == RankStrs[3])
//  							|| (i == 1 && SphereConstrainedNodeTypeStrs[n] == RankStrs[1]))
//  							&& std::find(SphereNodeNums[i][b].begin(), SphereNodeNums[i][b].end(), SphereConstrainedNodeNums[n]) != SphereNodeNums[i][b].end())
//  						{
//  							SphereConstrainedNodeNum = n;
//  							break;
//  						}
// 					}
// 					This is now taken care of in a much more robust way above.

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

						if (CondensedBasinToConstrainedNodeNums[i][b].first >= 0) {
							AuxDataZoneSetItem(ZoneNum, CSMAuxData.GBA.CondensedBasinName, SphereConstrainedNodeNames[CondensedBasinToConstrainedNodeNums[i][b].first]);
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
  							if (CondensedBasinToConstrainedNodeNums[i][b].first >= 0) {
  								AuxDataZoneSetItem(ZoneNum, CSMAuxData.GBA.CondensedBasinName, SphereConstrainedNodeNames[CondensedBasinToConstrainedNodeNums[i][b].first]);
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


		if (!oldZoneSet.isEmpty()) {
			TecUtilDataSetDeleteZone(oldZoneSet.getRef());
		}
	}

	Set DelVars;
	if (DeleteGradVars) {
		for (int i : GradVarNums)
			DelVars += i;
	}

	if (!DelVars.isEmpty()) {
		TecUtilDataSetDeleteVar(DelVars.getRef());
	}

	GBAResultViewerPopulateGBs();
}

struct CreateGBProbeClientData_s
{
	VolExtentIndexWeights_s VolInfo;
	int VolZoneNum;
	vector<int> XYZVarNums;
	int RhoVarNum;
	vector<int> GradVarNums;
	double SeedRadius;
	int NumGPs;
	int PointsPerGP;
	bool SaveGPs;
	int GroupNum;
	bool CapSurf;
};

CreateGBProbeClientData_s ClientDataStruct;

void CreateCircularGBsProbeCB(Boolean_t WasSuccessful, Boolean_t IsNearestPoint, ArbParam_t ClientData)
{
	TecUtilLockStart(AddOnID);
	if (WasSuccessful) {

// 		CreateGBProbeClientData_s *ClientDataStruct = (CreateGBProbeClientData_s*)ClientData;

		int ZoneNum = TecUtilProbeFieldGetZone();
		if (AuxDataZoneItemMatches(ZoneNum, CSMAuxData.GBA.ZoneType, CSMAuxData.GBA.ZoneTypeSphereZone)) {
			CSMGUILock();

			string SphereName = AuxDataZoneGetItem(ZoneNum, CSMAuxData.GBA.SourceNucleusName);
			if (SphereName == "") SphereName = "GBA Sphere";

			vec3 SphereOrigin;
			if (!GetSphereOrigin(ZoneNum, SphereOrigin)) return;

			vec3 SeedOrigin;
			for (int i = 0; i < 3; ++i)
				SeedOrigin[i] = TecUtilProbeFieldGetValue(ClientDataStruct.XYZVarNums[i]);

// 			SaveVec3VecAsScatterZone({ SeedOrigin }, "seed origin");

			// Now get the vector to rotate around the seed origin
			vec3 RotAxis = normalise(SphereOrigin - SeedOrigin);
			vec3 v2 = RotAxis;
			v2[abs(v2).index_min()] += 1.0;
			v2 = normalise(v2);
			vec3 RotVec = normalise(cross(v2,RotAxis)) * ClientDataStruct.SeedRadius;
// 			for (int i = 0; i < 3; ++i){
// 				v2[i] += 100.0;
// 				RotVec = normalise(cross(RotAxis, v2));
// 				if (abs(dot(RotVec, RotAxis)) > 0.9)
// 					break;
// 			}
// 			RotVec *= ClientDataStruct.SeedRadius;

			// Get read pointers
			FieldDataPointer_c RhoPtr;
			vector<FieldDataPointer_c> GradPtrs;

			TecUtilDataLoadBegin();
			GetReadPtrsForZone(ClientDataStruct.VolZoneNum, ClientDataStruct.RhoVarNum, ClientDataStruct.GradVarNums, vector<int>(), RhoPtr, GradPtrs, vector<FieldDataPointer_c>());

			int NumThreads = omp_get_num_procs();

			vector<GradPath_c> GPs(ClientDataStruct.NumGPs);
			vector<GradPath_c const*> GPPtrs(ClientDataStruct.NumGPs);
			vector<vec3> SeedPts;
			SeedPts.reserve(ClientDataStruct.NumGPs + 1);

			for (int i = 0; i < ClientDataStruct.NumGPs; ++i){
				// Get seed point (assume we don't need to project it back to sphere)
				vec3 SeedPt = SeedOrigin + Rotate( RotVec, TWOPI * (double)i / (double)ClientDataStruct.NumGPs, RotAxis);
// 				SaveVec3VecAsScatterZone({ SeedPt }, "seed origin");
				SeedPts.push_back(SeedPt);
				GPs[i].SetupGradPath(SeedPt, StreamDir_Both, ClientDataStruct.PointsPerGP, GPType_Classic, GPTerminate_AtBoundary, nullptr, nullptr, nullptr, nullptr, ClientDataStruct.VolInfo, vector<FieldDataPointer_c>(), GradPtrs, RhoPtr);
				GPPtrs[i] = &GPs[i];
			}
			SeedPts.push_back(SeedOrigin);

#pragma omp parallel for schedule(dynamic)
			for (int i = 0; i < ClientDataStruct.NumGPs; ++i){
				GPs[i].Seed();
			}

			int NumVars = TecUtilDataSetGetNumVars();

			if (ClientDataStruct.SaveGPs){
				for (int i = 0; i < GPs.size(); ++i){
					GPs[i].SaveAsOrderedZone("Circular GB Path " + to_string(i + 1));
				}
			}

			FESurface_c Surf;
			Surf.MakeFromGPs(GPPtrs, true, ClientDataStruct.CapSurf, true);
			int SurfZoneNum = Surf.SaveAsTriFEZone("GB (" + SphereName + ", r=" + DoubleToString(ClientDataStruct.SeedRadius) + ")", vector<FieldDataType_e>(NumVars, FieldDataType_Float), vector<ValueLocation_e>(NumVars, ValueLocation_Nodal), ClientDataStruct.XYZVarNums);

			if (SurfZoneNum > 0) {

				// Write rho values to each node on the GB

				FieldVecPointer_c SurfXYZPtr;
				SurfXYZPtr.InitializeReadPtr(SurfZoneNum, ClientDataStruct.XYZVarNums);

				FieldDataPointer_c SurfRhoPtr;
				SurfRhoPtr.InitializeWritePtr(SurfZoneNum, ClientDataStruct.RhoVarNum);

				int NumNodes = SurfXYZPtr.Size();
				vector<VolExtentIndexWeights_s> ThVolInfo(omp_get_num_procs(), ClientDataStruct.VolInfo);
#pragma omp parallel for
				for (int i = 0; i < NumNodes; ++i) {
					SurfRhoPtr.Write(i, ValAtPointByPtr(SurfXYZPtr[i], ThVolInfo[omp_get_thread_num()], RhoPtr));
				}

				// Write any GBA integrated values to ALL nodes on the GB
				vector<string> SearchStrs = { "I: ","IN: ","INS: " };
				for (int v = 1; v <= NumVars; ++v) {
					char *TmpCStr;
					if (TecUtilVarGetName(v, &TmpCStr)) {
						string TmpStr = TmpCStr;
						TecUtilStringDealloc(&TmpCStr);
						for (auto const & s : SearchStrs) {
							if (TmpStr.find(s) == 0) {
								double VarVal = TecUtilProbeFieldGetValue(v);
								FieldData_pa TmpPtr = TecUtilDataValueGetWritableNativeRef(SurfZoneNum, v);
								

								if (TecUtilDataValueGetType(SurfZoneNum, v) == FieldDataType_Double) {
									vector<double> ValVec(NumNodes, VarVal);
									TecUtilDataValueArraySetByRef(TmpPtr, 1, NumNodes, ValVec.data());
								}
								else {
									vector<float> ValVec(NumNodes, VarVal);
									TecUtilDataValueArraySetByRef(TmpPtr, 1, NumNodes, ValVec.data());
								}

								break;
							}
						}
					}
				}

				TecUtilDataLoadEnd();
				// 			int SurfZoneNum = Surf.SaveAsTriFEZone(ClientDataStruct.XYZVarNums, "GB (" + SphereName + ", r=" + DoubleToString(ClientDataStruct.SeedRadius) + ")");
				Set SurfZoneSet(SurfZoneNum);
				SetZoneStyle({ SurfZoneNum }, ZoneStyle_Surface);
				TecUtilZoneSetActive(SurfZoneSet.getRef(), AssignOp_PlusEquals);

				TecUtilZoneSetMesh(SV_SHOW, SurfZoneSet.getRef(), 0.0, FALSE);
				TecUtilZoneSetShade(SV_SHOW, SurfZoneSet.getRef(), 0.0, FALSE);
				TecUtilZoneSetScatter(SV_SHOW, SurfZoneSet.getRef(), 0.0, FALSE);
				TecUtilZoneSetEdgeLayer(SV_SHOW, SurfZoneSet.getRef(), 0.0, FALSE);

				TecUtilZoneSetContour(SV_SHOW, SurfZoneSet.getRef(), 0.0, TRUE);
				TecUtilZoneSetContour(SV_CONTOURTYPE, SurfZoneSet.getRef(), 0.0, ContourType_Flood);
				TecUtilZoneSetContour(SV_FLOODCOLORING, SurfZoneSet.getRef(), 0.0, ContourColoring_Group8);

				/*
				 *	Change group number of zones
				 */
				if (ClientDataStruct.GroupNum > 1) {
					TecUtilMacroExecuteCommand(string("$!FIELDMAP [" + to_string(SurfZoneNum) + "] GROUP = " + to_string(ClientDataStruct.GroupNum)).c_str());
				}

				// Save disc of seed points as zone
				vector<vector<int> > DiscElems(ClientDataStruct.NumGPs,vector<int>(3));
				for (int i = 0; i < ClientDataStruct.NumGPs; ++i) {
					DiscElems[i] = { i, (i + 1) % ClientDataStruct.NumGPs, ClientDataStruct.NumGPs };
					SeedPts[i] = SphereOrigin + (SeedPts[i] - SphereOrigin) * 1.000001;
				}
				SeedPts.back() = SphereOrigin + (SeedPts.back() - SphereOrigin) * 1.000001;


				FESurface_c SeedDisc(SeedPts, DiscElems);
				int DiscZoneNum = SeedDisc.SaveAsTriFEZone("GB Seed Disc (" + SphereName + ", r=" + DoubleToString(ClientDataStruct.SeedRadius) + ")", vector<FieldDataType_e>(NumVars, FieldDataType_Float), vector<ValueLocation_e>(NumVars, ValueLocation_Nodal), ClientDataStruct.XYZVarNums);
				if (DiscZoneNum > 0){
					// Write any GBA integrated values to ALL nodes on the GB seed disc
					for (int v = 1; v <= NumVars; ++v) {
						char *TmpCStr;
						if (TecUtilVarGetName(v, &TmpCStr)) {
							string TmpStr = TmpCStr;
							TecUtilStringDealloc(&TmpCStr);
							for (auto const & s : SearchStrs) {
								if (TmpStr.find(s) == 0) {
									double VarVal = TecUtilProbeFieldGetValue(v);
									FieldData_pa TmpPtr = TecUtilDataValueGetWritableNativeRef(DiscZoneNum, v);

									if (TecUtilDataValueGetType(DiscZoneNum, v) == FieldDataType_Double) {
										vector<double> ValVec(ClientDataStruct.NumGPs + 1, VarVal);
										TecUtilDataValueArraySetByRef(TmpPtr, 1, ClientDataStruct.NumGPs + 1, ValVec.data());
									}
									else {
										vector<float> ValVec(ClientDataStruct.NumGPs + 1, VarVal);
										TecUtilDataValueArraySetByRef(TmpPtr, 1, ClientDataStruct.NumGPs + 1, ValVec.data());
									}

									break;
								}
							}
						}
					}

					Set DiscZoneSet(DiscZoneNum);
					TecUtilZoneSetActive(DiscZoneSet.getRef(), AssignOp_MinusEquals);

					TecUtilZoneSetMesh(SV_SHOW, DiscZoneSet.getRef(), 0.0, FALSE);
					TecUtilZoneSetScatter(SV_SHOW, DiscZoneSet.getRef(), 0.0, FALSE);
					TecUtilZoneSetEdgeLayer(SV_SHOW, DiscZoneSet.getRef(), 0.0, TRUE);

					TecUtilZoneSetContour(SV_SHOW, DiscZoneSet.getRef(), 0.0, TRUE);
					TecUtilZoneSetContour(SV_CONTOURTYPE, DiscZoneSet.getRef(), 0.0, ContourType_Flood);
					TecUtilZoneSetContour(SV_FLOODCOLORING, DiscZoneSet.getRef(), 0.0, ContourColoring_Group8);

					TecUtilZoneSetShade(SV_SHOW, DiscZoneSet.getRef(), 0.0, FALSE);
					TecUtilZoneSetShade(SV_COLOR, DiscZoneSet.getRef(), 0.0, Red_C);
				}
			}
			else
				TecUtilDataLoadEnd();

			CSMGUIUnlock();
		}
	}
	TecUtilLockFinish(AddOnID);
}

void STDCALL CreateCircularGBsReturnUserInfo(bool const GuiSuccess,
	vector<GuiField_c> const & Fields,
	vector<GuiField_c> const PassthroughFields) {
	if (!GuiSuccess) return;

	TecUtilLockStart(AddOnID);

	Boolean_t IsPeriodic;

// 	CreateGBProbeClientData_s ClientDataStruct;
	ClientDataStruct.XYZVarNums.resize(3);

	int fNum = 0;

	ClientDataStruct.VolZoneNum = Fields[fNum++].GetReturnInt();
	fNum++;
	for (int i = 0; i < 3; ++i) ClientDataStruct.XYZVarNums[i] = i + Fields[fNum].GetReturnInt();
	fNum++;
	fNum++;
	ClientDataStruct.RhoVarNum = Fields[fNum++].GetReturnInt();
	fNum++;


	if (Fields[fNum++].GetReturnBool()) {
		ClientDataStruct.GradVarNums.resize(3);
		for (int i = 0; i < 3; ++i) ClientDataStruct.GradVarNums[i] = i + Fields[fNum].GetReturnInt();
	}
	fNum++;
	fNum++;

	ClientDataStruct.SeedRadius = Fields[fNum++].GetReturnDouble();
	ClientDataStruct.NumGPs = Fields[fNum++].GetReturnInt();
	ClientDataStruct.PointsPerGP = Fields[fNum++].GetReturnInt();
	ClientDataStruct.SaveGPs = Fields[fNum++].GetReturnBool();
	ClientDataStruct.GroupNum = Fields[fNum++].GetReturnInt();
	ClientDataStruct.CapSurf = Fields[fNum++].GetReturnBool();

	GetVolInfo(ClientDataStruct.VolZoneNum, ClientDataStruct.XYZVarNums, FALSE, ClientDataStruct.VolInfo);

	ArgList Args;
	Args.appendString(SV_STATUSLINETEXT, "Click at center of gradient bundle on a GBA sphere");
// 	Args.appendArbParam(SV_CLIENTDATA, ArbParam_t(&ClientDataStruct));
	Args.appendFunction(SV_CALLBACKFUNCTION, CreateCircularGBsProbeCB);

 	TecUtilProbeInstallCallbackX(Args.getRef());

	TecUtilLockFinish(AddOnID);
}

void CreateCircularGBsGetUserInfo() {

	vector<GuiField_c> Fields = {
		GuiField_c(Gui_ZoneSelect, "Volume zone", CSMZoneName.FullVolume.substr(0,10)),
		GuiField_c(Gui_VertSep)
	};

	// 	vector<string> XYZStr = { "X", "Y", "Z" }, HessXYZStr = { "XX", "XY", "XZ", "YY", "YZ", "ZZ" };

	Fields.emplace_back(Gui_VarSelect, "X", "X");
	// 	for (int i = 0; i < 3; ++i) Fields.push_back(GuiField_c(Gui_VarSelect, XYZStr[i], XYZStr[i]));
	Fields.emplace_back(Gui_VertSep);

	Fields.emplace_back(Gui_VarSelect, "Electron Density", CSMVarName.Dens);
	Fields.emplace_back(Gui_VertSep);

	int iTmp = Fields.size();
	Fields.emplace_back(Gui_ToggleEnable, "Density gradient vector variables present");
	// 	for (int i = 0; i < 3; ++i){
	// 		Fields[iTmp].AppendSearchString(to_string(Fields.size()));
	// 		Fields.push_back(GuiField_c(Gui_VarSelect, XYZStr[i], CSMVarName.DensGradVec[i]));
	// 		if (i < 2) Fields[iTmp].AppendSearchString(",");
	// 	}
	Fields[iTmp].AppendSearchString(to_string(Fields.size()));
	Fields.emplace_back(Gui_VarSelect, CSMVarName.DensGradVec[0], CSMVarName.DensGradVec[0]);
	Fields.emplace_back(Gui_VertSep);

	Fields.emplace_back(Gui_Double, "GB seed radius", "0.1");
	Fields.emplace_back(Gui_Int, "Number of GPs", "30");
	Fields.emplace_back(Gui_Int, "Points per GP", "300");
	Fields.emplace_back(Gui_Toggle, "Save GPs", "0");
	Fields.emplace_back(Gui_Int, "Dest. group number", "2");
	Fields.emplace_back(Gui_Toggle, "Cap gradient bundle", "1");


	CSMGui("Circular gradient bundle", Fields, CreateCircularGBsReturnUserInfo, AddOnID);
}


/*
 *	This is a function that lets a user select two points on a GBA sphere, then, using the sphere center as the third point, make a slice through the sphere, interpolate condensed property values along the intersection of the slice and the sphere, and make a 1d plot from the results.
 *	Here's the basic process:
 *	1. Given the two selected points and the sphere center, construct orthonormal axes for the plane.
 *	2. Find all sphere edges that intersect the plane and get interpolated condensed properties for the intersections.
 *		(http://www.geomalgorithms.com/a05-_intersect-1.html)
 *	3. Map all the intersections onto x' and y' (the two in-plane directions for the plane). Just take dot product with x' and y' vectors.
 *	4. Get theta values based on the x' and y' values.
 *	5. Sort based on theta value.
 *	6. Shift theta values so that the two selected points will appear centered in the resulting plot.
 *	7. [If gradient bundles have been found, then we want to also plot integration values for applicable condensed maximum and minumum basins.]
 *		Every element on the sphere belongs to one maximum and one minimum basin, but we only care about having these values for the maxima and minima of our resulting plot
 *		a. For each condensed variable in the 1d plot data, for each max (min), need to find the maximum (minimum) basin it belongs to, then store the total condensed integration value for that variable.
 *		b. The resulting values will be stored in two ways; as zones for making a 1d plot in Tecplot, and as AuxData to the main zone of the 1d plot so that it's available in an exporded .dat file for the first zone of the 1d plot.
 *		c. In the resulting plot, we'll plot the total integration values in tecplot as scatter points at the same theta value as the max/min they correspond to.
 *			In Mathematica, where we want to plot this type of 1d plot for multiple systems simultaneously, I'll have to think of a better way to show the information, since multiple colored lines with multiple similarly colored points could get messy quick.
 */
vector<vec3> Make1dPlotProbePoints;
void Make1dPlotFromSphereSliceProbeCB(Boolean_t WasSuccessful, Boolean_t IsNearestPoint, ArbParam_t ClientData)
{

}