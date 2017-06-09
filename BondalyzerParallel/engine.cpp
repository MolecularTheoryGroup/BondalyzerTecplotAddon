/*
*****************************************************************
*****************************************************************
*******                                                  ********
*******    (C) Copyright 1988-2013 by TECPLOT INC.       *******
*******                                                  ********
*****************************************************************
*****************************************************************
*/

#include <vector>
#include <string>
#include <string.h>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <ctime>
#include <cstdlib>
#include<algorithm>

#include <omp.h>

// #include <windows.h>

#include "TECADDON.h"
#include "ADDGLBL.h"
#include "ADDONVER.h"
#include "GUIDEFS.h"
#include "ADKUTIL.h"

#include "CSM_DATA_TYPES.h"
#include "CSM_DATA_SET_INFO.h"
#include "CSM_CRIT_POINTS.h"
#include "CSM_CALC_VARS.h"
#include "CSM_GRAD_PATH.h"
#include "CSM_FE_VOLUME.h"

#include "KFc.h"

#include "ENGINE.h"

#include <armadillo>
using namespace arma;


using std::string;
using std::to_string;
using std::vector;
using std::stringstream;
using std::ofstream;
using std::ios;
using std::endl;
using std::setprecision;




const string DS_ZoneName_Delim = "_-_";
const string T41Ext = ".t41";


void RefineActiveZones(){
	/*
	 *	Get number of times to refine
	 */

	char* UserText;
	if (!TecUtilDialogGetSimpleText("Enter number times to refine zones", "2", &UserText))
		return;
	int NumRefine = std::stoi(UserText);
	TecUtilStringDealloc(&UserText);

	vector<FEVolume_c> Vols;


	TecUtilDialogLaunchPercentDone("Refining zones", TRUE);

	for (int i = 1; i <= TecUtilDataSetGetNumZones(); ++i){
		FEVolume_c Vol(i, vector<int>({ 1, 2, 3 }));
		if (Vol.IsMade()){
			Vols.push_back(Vol);
		}
	}

	bool UserQuit = false;
	int NumVols = Vols.size();
	int NumVolsPerThread = int(double(NumVols) / double(omp_get_num_procs()));
	int VolNum = 0;
#ifndef _DEBUG
#pragma omp parallel for
#endif
	for (int i = 0; i < NumVols; ++i){
		if (omp_get_thread_num() == 0){
			UserQuit = !TecUtilDialogCheckPercentDone(MIN(100, MAX(int(double(++VolNum) / double(NumVols) * 100.0), 1)));
#pragma omp flush (UserQuit)
		}
		else{
			++VolNum;
		}
#pragma omp flush (UserQuit)
		if (!UserQuit){
			for (int r = 0; r < NumRefine; ++r)
				Vols[i].Refine();
		}
	}

	for (int i = 0; i < NumVols; ++i){
		char* ZoneName;
		TecUtilZoneGetName(Vols[i].GetZoneNum(), &ZoneName);
		Vols[i].SaveAsTriFEZone(vector<int>({ 1, 2, 3 }), string(ZoneName) + " (refined)");
		TecUtilStringDealloc(&ZoneName);
	}

	TecUtilDialogDropPercentDone();
}

const bool FEZoneDFS(const int & NodeNum,
				const NodeMap_pa & NodeMap,
				const NodeToElemMap_pa & NodeToElemMap,
				const int & ZoneNum,
				vector<bool> & IsVisited){
	bool IsOk = true;
	int NumElems = TecUtilDataNodeToElemMapGetNumElems(NodeToElemMap, NodeNum + 1);
	for (int e = 1; e <= NumElems && IsOk; ++e){
		int ei = TecUtilDataNodeToElemMapGetElem(NodeToElemMap, NodeNum + 1, e);
		for (int n = 0; n < 3 && IsOk; ++n){
			int ni = TecUtilDataNodeGetByRef(NodeMap, ei, n + 1) - 1;
			if (!IsVisited[ni]){
				IsVisited[ni] = true;
				IsOk = FEZoneDFS(ni, NodeMap, NodeToElemMap, ZoneNum, IsVisited);
			}
		}
	}

	return IsOk;
}

void GetClosedIsoSurfaceFromPoints(){
	int IsoZoneNum = MAX(1, ZoneNumByName("Iso"));
	int	PointZoneNum = MAX(1, ZoneNumByName("Critical Points"));
	vector<int> PointINums;

	int NumZones = TecUtilDataSetGetNumZones();

	char *UserInput;

	while (true){
		if (!TecUtilDialogGetSimpleText("Generate closed isosurface. Enter Isosurface zone number and Point zone number, separated by commas", string(to_string(IsoZoneNum) + "," + to_string(PointZoneNum)).c_str(), &UserInput)){
			return;
		}

		vector<string> Answers = SplitString(string(UserInput), ',');
		TecUtilStringDealloc(&UserInput);
		if (Answers.size() != 2){
			TecUtilDialogErrMsg("Expected 2 numbers. Try again.");
			IsoZoneNum = MAX(1, ZoneNumByName("Iso"));
			PointZoneNum = MAX(1, ZoneNumByName("Critical Points"));
		}
		else{
			IsoZoneNum = stoi(Answers[0]);
			PointZoneNum = stoi(Answers[1]);
			break;
		}
	}

	for (auto i : vector<int>{ IsoZoneNum, PointZoneNum }){
		if (i < 1 || i > NumZones){
			TecUtilDialogErrMsg("Invalid zone number provided. Quitting.");
			return;
		}
	}

	int PointIJK[3], IsoIJK[3];
	TecUtilZoneGetIJK(IsoZoneNum, &IsoIJK[0], &IsoIJK[1], &IsoIJK[2]); // {NumNodes, NumElems, NumNodesPerElem}
	TecUtilZoneGetIJK(PointZoneNum, &PointIJK[0], &PointIJK[1], &PointIJK[2]);

	if (!TecUtilZoneIsFiniteElement(IsoZoneNum)){
		TecUtilDialogErrMsg("Isosurface zone is not finite element. Quitting.");
		return;
	}
	else if (!TecUtilZoneIsOrdered(PointZoneNum) || PointIJK[1] > 1 || PointIJK[2] > 1){
		TecUtilDialogErrMsg("Point zone is not i-ordered. Quitting.");
		return;
	}

	while (true){
		if (!TecUtilDialogGetSimpleText("Enter I-number(s) of point(s) around which to extract closed isosurfaces, separated by commas", "1", &UserInput)){
			return;
		}

		vector<string> Answers = SplitString(string(UserInput), ',');
		TecUtilStringDealloc(&UserInput);

		for (const auto & i : Answers){
			PointINums.push_back(stoi(i));
			if (PointINums.back() < 1 || PointINums.back() > PointIJK[0]){
				TecUtilDialogErrMsg("Invalid point i number. Try again.");
				PointINums.clear();
				break;
			}
		}

		if (PointINums.size() > 0) break;
	}

	TecUtilDataLoadBegin();

	vector<FieldDataPointer_c> IsoReadPtrs(TecUtilDataSetGetNumVars());
	for (int i = 0; i < IsoReadPtrs.size(); ++i){
		if (!IsoReadPtrs[i].GetReadPtr(IsoZoneNum, i + 1)){
			TecUtilDialogErrMsg("Failed to get isosurface read pointer. Quitting.");
			return;
		}
	}

	vector<int> NodeNums;

	for (const int & PointINum : PointINums){
		vec3 Pt;
		for (int i = 0; i < 3; ++i) Pt[i] = TecUtilDataValueGetByZoneVar(PointZoneNum, i + 1, PointINum);

		// Get closest node in isosurface to point
		// Assume xyz vars are {1,2,3}


		double MinDistSqr = DBL_MAX, TmpDistSqr;
		int NodeNum = -1;
		vec3 TmpVec;
		for (int NodeNum = 0; NodeNum < IsoIJK[0]; ++NodeNum){
			for (int i = 0; i < 3; ++i) TmpVec[i] = IsoReadPtrs[i][NodeNum];
			TmpDistSqr = DistSqr(Pt, TmpVec);
			if (TmpDistSqr < MinDistSqr){
				MinDistSqr = TmpDistSqr;
				NodeNum = NodeNum;
			}
		}

		if (NodeNum < 0){
			TecUtilDialogErrMsg("Failed to find node of minimum distance. Quitting.");
			return;
		}

		NodeNums.push_back(NodeNum);
	}

	GetClosedIsoSurface(IsoZoneNum, IsoReadPtrs, NodeNums);

	TecUtilDataLoadEnd();
}

void GetClosedIsoSurfaceFromNodes(){
	int IsoZoneNum = MAX(1, ZoneNumByName("Iso"));
	vector<int> NodeNums;

	int NumZones = TecUtilDataSetGetNumZones();

	char *UserInput;

	while (true){
		if (!TecUtilDialogGetSimpleText("Generate closed isosurface. Enter Isosurface zone number:", string(to_string(IsoZoneNum)).c_str(), &UserInput)){
			return;
		}

		vector<string> Answers = SplitString(string(UserInput), ',');
		TecUtilStringDealloc(&UserInput);
		if (Answers.size() != 1){
			TecUtilDialogErrMsg("Expected 1 number. Try again.");
			IsoZoneNum = MAX(1, ZoneNumByName("Iso"));
		}
		else{
			IsoZoneNum = stoi(Answers[0]);
			break;
		}
	}

	if (IsoZoneNum < 1 || IsoZoneNum > NumZones){
		TecUtilDialogErrMsg("Invalid zone number provided. Quitting.");
		return;
	}

	int IsoIJK[3];
	TecUtilZoneGetIJK(IsoZoneNum, &IsoIJK[0], &IsoIJK[1], &IsoIJK[2]); // {NumNodes, NumElems, NumNodesPerElem}

	if (!TecUtilZoneIsFiniteElement(IsoZoneNum)){
		TecUtilDialogErrMsg("Isosurface zone is not finite element. Quitting.");
		return;
	}

	while (true){
		if (!TecUtilDialogGetSimpleText("Enter I-number(s) of point(s) around which to extract closed isosurfaces, separated by commas", "1", &UserInput)){
			return;
		}

		vector<string> Answers = SplitString(string(UserInput), ',');
		TecUtilStringDealloc(&UserInput);

		for (const auto & i : Answers){
			NodeNums.push_back(stoi(i) - 1);
			if (NodeNums.back() < 1 || NodeNums.back() > IsoIJK[0]){
				TecUtilDialogErrMsg("Invalid point i number. Try again.");
				NodeNums.clear();
				break;
			}
		}

		if (NodeNums.size() > 0) break;
	}

	TecUtilDataLoadBegin();

	vector<FieldDataPointer_c> IsoReadPtrs(TecUtilDataSetGetNumVars());
	for (int i = 0; i < IsoReadPtrs.size(); ++i){
		if (!IsoReadPtrs[i].GetReadPtr(IsoZoneNum, i + 1)){
			TecUtilDialogErrMsg("Failed to get isosurface read pointer. Quitting.");
			return;
		}
	}

	GetClosedIsoSurface(IsoZoneNum, IsoReadPtrs, NodeNums);

	TecUtilDataLoadEnd();
}

void GetClosedIsoSurface(const int & IsoZoneNum, const vector<FieldDataPointer_c> & IsoReadPtrs, const vector<int> & NodeNums){

	int IsoIJK[3];
	TecUtilZoneGetIJK(IsoZoneNum, &IsoIJK[0], &IsoIJK[1], &IsoIJK[2]); // {NumNodes, NumElems, NumNodesPerElem}

	char *IsoZoneName;
	TecUtilZoneGetName(IsoZoneNum, &IsoZoneName);

	Set_pa ZoneSet = TecUtilSetAlloc(TRUE);

	TecUtilDialogLaunchPercentDone("Generating isosurface subzone", FALSE);
	int PtNum = 1;

	NodeMap_pa NodeMap = TecUtilDataNodeGetReadableRef(IsoZoneNum);
	NodeToElemMap_pa NodeToElemMap = TecUtilDataNodeToElemMapGetReadableRef(IsoZoneNum);

	for (const int & NodeNum : NodeNums){
		if (!TecUtilDialogCheckPercentDone(0)){
			TecUtilDialogDropPercentDone();
			return;
		}
		TecUtilDialogSetPercentDoneText(string("Generating isosurface subzone: " + to_string(PtNum++) + " of " + to_string(NodeNums.size())).c_str());


		/*
		 * Now run a depth first search through the isosurface to get the connected component
		 * of the minimum distance node.
		 */
		vector<bool> NodeVisited(IsoIJK[0], false);

		if (!VALID_REF(NodeMap)){
			TecUtilDialogErrMsg("Isosurface node pointer(s) invalid. Quitting.");
			return;
		}

		FEZoneDFS(NodeNum, NodeMap, NodeToElemMap, IsoZoneNum, NodeVisited);

		/*
		 * Collect nodes/elements of single found connected component
		 */

		int NumNewNodes = 0;
		for (const auto & i : NodeVisited) NumNewNodes += (int)i;

		if (NumNewNodes <= 0){
			continue;
		}

		vector<int> NodeNumsNewToOld(NumNewNodes);
		vector<vector<int> > Elems(IsoIJK[1]);
		for (auto & i : Elems) i.reserve(3);

		int NodeNumElems;
		int NewNodeNum = -1;

		for (int n = 0; n < IsoIJK[0]; ++n){
			if (NodeVisited[n]){
				NodeNumsNewToOld[++NewNodeNum] = n;
				NodeNumElems = TecUtilDataNodeToElemMapGetNumElems(NodeToElemMap, n + 1);
				for (int e = 1; e <= NodeNumElems; ++e) Elems[TecUtilDataNodeToElemMapGetElem(NodeToElemMap, n + 1, e) - 1].push_back(NewNodeNum);
			}
		}

		vector<FieldDataPointer_c> IsoWritePtrs(IsoReadPtrs.size());

		vector<vector<int> > NewElems;
		NewElems.reserve(IsoIJK[1]);
		for (const auto & i : Elems){
			if (i.size() == 3) NewElems.push_back(i);
			else if (i.size() > 0 && i.size() != 3){
				TecUtilDialogErrMsg("Element with less/more than 3 nodes found, and I didn't plan for that.... Quitting.");
				return;
			}
		}

		/*
		 * Make new zone with single connected component
		 */

		if (!TecUtilDataSetAddZone(string("CP " + to_string(NodeNum) + ": " + IsoZoneName).c_str(), NumNewNodes, NewElems.size(), 3, ZoneType_FETriangle, NULL)){
			TecUtilDialogErrMsg("Failed to make new iso zone. Quitting.");
			return;
		}

		int NewZoneNum = TecUtilDataSetGetNumZones();

		for (int i = 0; i < IsoReadPtrs.size(); ++i){
			if (!IsoWritePtrs[i].GetWritePtr(NewZoneNum, i + 1)){
				TecUtilDialogErrMsg("Failed to get write pointer(s) for new iso zone. Quitting.");
				return;
			}
		}

		for (int n = 0; n < NumNewNodes; ++n) for (int i = 0; i < IsoWritePtrs.size(); ++i) IsoWritePtrs[i].Write(n, IsoReadPtrs[i][NodeNumsNewToOld[n]]);

		NodeMap_pa NewNodeMap = TecUtilDataNodeGetWritableRef(NewZoneNum);
		if (!VALID_REF(NewNodeMap)){
			TecUtilDialogErrMsg("Failed to get node map pointer for new iso zone. Quitting.");
			return;
		}

		for (int e = 0; e < NewElems.size(); ++e) for (int ei = 0; ei < 3; ++ei) TecUtilDataNodeSetByRef(NewNodeMap, e + 1, ei + 1, NewElems[e][ei] + 1);

		TecUtilSetAddMember(ZoneSet, NewZoneNum, TRUE);

	}
	
	if (TecUtilSetGetMemberCount(ZoneSet) > 0){
		TecUtilZoneSetMesh(SV_SHOW, ZoneSet, 0.0, FALSE);
		TecUtilZoneSetScatter(SV_SHOW, ZoneSet, 0.0, FALSE);
		TecUtilZoneSetShade(SV_SHOW, ZoneSet, 0.0, FALSE);
		TecUtilZoneSetContour(SV_SHOW, ZoneSet, 0.0, TRUE);
		TecUtilZoneSetActive(ZoneSet, AssignOp_PlusEquals);
	}

	TecUtilStringDealloc(&IsoZoneName);
	TecUtilSetDealloc(&ZoneSet);
	TecUtilDialogDropPercentDone();
}

void GetInfoFromUserForBondalyzer(BondalyzerSteps_e CalcType){
	int NumZones = TecUtilDataSetGetNumZones();
	int NumVars = TecUtilDataSetGetNumVars();

	int RhoVarNum = MAX(1, VarNumByName("Electron Density"));
	vector<int> XYZVarNums(3), GradVarNums(3), HessVarNums(6);
	vector<string> XYZStr = { "X", "Y", "Z" }, HessXYZStr = { "XX", "XY", "XZ", "YY", "YZ", "ZZ" };
	for (int i = 0; i < 3; ++i){
		XYZVarNums[i] = MAX(1, VarNumByName(XYZStr[i]));
		GradVarNums[i] = MAX(1, VarNumByName(CSMVarName.DensGradVec[i]));
	}
	for (int i = 0; i < 6; ++i)
		HessVarNums[i] = MAX(1, VarNumByName(CSMVarName.DensHessTensor[i]));
	int VolZoneNum = MAX(1, ZoneNumByName("Full Volume"));
	Boolean_t IsPeriodic;

	while (true){
		char *c;
		if (TecUtilDialogGetSimpleText("Enter volume zone number to search for critical points", to_string(VolZoneNum).c_str(), &c)){
			if (StringIsInt(c)){
				VolZoneNum = atoi(c);
				TecUtilStringDealloc(&c);
				if (0 < VolZoneNum && VolZoneNum <= NumZones) break;
			}
			else TecUtilDialogErrMsg("Enter integer.");
		}
		else return;
	}
	while (true){
		char *c;
		if (TecUtilDialogGetSimpleText("Periodic boundary condition?\ny for yes, n for no.", "n", &c)){
			if (*c == 'y') IsPeriodic = TRUE;
			else if (*c == 'n') IsPeriodic = FALSE;
			else continue;
			break;
		}
		else return;
	}

	if (!TecUtilDialogGetVariables("Select electron density variable",
		"Electron density", NULL, NULL, &RhoVarNum, NULL, NULL)
		|| !TecUtilDialogGetVariables("Select the spatial X, Y, and Z variables",
		"X", "Y", "Z", &XYZVarNums[0], &XYZVarNums[1], &XYZVarNums[2])) return;

	if (TecUtilDialogMessageBox("Are electron density gradient variables present?", MessageBoxType_YesNo)){
		if (!TecUtilDialogGetVariables("Select the density gradient X, Y, and Z variables",
			"X gradient", "Y gradient", "Z gradient",
			&GradVarNums[0], &GradVarNums[1], &GradVarNums[2])) return;
	}
	else GradVarNums = vector<int>();

	if (GradVarNums.size() != 0 && TecUtilDialogMessageBox("Are electron density hessian variables present?", MessageBoxType_YesNo)){
		if (!TecUtilDialogGetVariables("Select the density hessian XX, XY, and XZ variables",
			"XX", "XY", "XZ",
			&HessVarNums[0], &HessVarNums[1], &HessVarNums[2])
			|| !TecUtilDialogGetVariables("Select the density hessian YY, YZ, and ZZ variables",
			"YY", "YZ", "ZZ",
			&HessVarNums[3], &HessVarNums[4], &HessVarNums[5])) return;
	}
	else HessVarNums = vector<int>();

	if (CalcType == CRITICALPOINTS){
		FindCritPoints(VolZoneNum, XYZVarNums, RhoVarNum, GradVarNums, HessVarNums, IsPeriodic);
	}
	else if (CalcType == BONDPATHS || CalcType == RINGLINES){
		int CPType;
		if (CalcType == BONDPATHS) CPType = BONDCP;
		else CPType = RINGCP;

		int CPZoneNum = MAX(1, ZoneNumByName("Critical Points", false, true));
		int CPTypeVarNum = MAX(1, VarNumByName(CPTypeVarName));

		while (true){
			char *c;
			if (TecUtilDialogGetSimpleText("Enter critical points zone number", to_string(CPZoneNum).c_str(), &c)){
				if (StringIsInt(c)){
					CPZoneNum = atoi(c);
					TecUtilStringDealloc(&c);
					if (0 < CPZoneNum && CPZoneNum <= NumZones) break;
				}
				else TecUtilDialogErrMsg("Enter integer.");
			}
			else return;
		}

		if (!TecUtilDialogGetVariables("Select the critical point type variable",
			"CP type", NULL, NULL,
			&CPTypeVarNum, NULL, NULL)) return;

		FindBondRingLines(VolZoneNum, CPZoneNum, CPType, XYZVarNums, RhoVarNum, GradVarNums, HessVarNums, IsPeriodic);
	}
}

const Boolean_t GetReadPtrsForZone(const int & ZoneNum,
	const int & RhoVarNum,
	const vector<int> & GradVarNums,
	const vector<int> & HessVarNums,
	FieldDataPointer_c & RhoPtr,
	vector<FieldDataPointer_c> & GradPtrs,
	vector<FieldDataPointer_c> & HessPtrs)
{
	int NumZones = TecUtilDataSetGetNumZones();
	int NumVars = TecUtilDataSetGetNumVars();

	REQUIRE(ZoneNum > 0 && ZoneNum <= NumZones);
	REQUIRE(RhoVarNum > 0 && RhoVarNum <= NumVars);
	REQUIRE(GradVarNums.size() == 3 || GradVarNums.size() == 0);
	for (const auto & i : GradVarNums) REQUIRE(i > 0 && i <= NumVars);
	REQUIRE(HessVarNums.size() == 6 || HessVarNums.size() == 0);
	for (const auto & i : HessVarNums) REQUIRE(i > 0 && i <= NumVars);

	RhoPtr.Close();
	GradPtrs.resize(GradVarNums.size());
	HessPtrs.resize(HessVarNums.size());

	if (!RhoPtr.GetReadPtr(ZoneNum, RhoVarNum)) {
		TecUtilDialogErrMsg(string("Failed to get rho read pointer (zone " + to_string(ZoneNum) + ", var " + to_string(RhoVarNum) + ")").c_str());
		return FALSE;
	}
	for (int i = 0; i < GradVarNums.size(); ++i) if (!GradPtrs[i].GetReadPtr(ZoneNum, GradVarNums[i])){
		TecUtilDialogErrMsg(string("Failed to get gradient read pointer (zone " + to_string(ZoneNum) + ", var " + to_string(GradVarNums[i]) + ")").c_str());
		return FALSE;
	}
	for (int i = 0; i < HessVarNums.size(); ++i) if (!HessPtrs[i].GetReadPtr(ZoneNum, HessVarNums[i])){
		TecUtilDialogErrMsg(string("Failed to get hessian read pointer (zone " + to_string(ZoneNum) + ", var " + to_string(HessVarNums[i]) + ")").c_str());
		return FALSE;
	}

	return TRUE;
}

const Boolean_t GetVolInfo(const int & VolZoneNum,
	const vector<int> & XYZVarNums,
	const Boolean_t & IsPeriodic,
	VolExtentIndexWeights_s & VolInfo)
{
	TecUtilZoneGetIJK(VolZoneNum, &VolInfo.MaxIJK[0], &VolInfo.MaxIJK[1], &VolInfo.MaxIJK[2]);
	for (int i = 0; i < 3; ++i) REQUIRE(VolInfo.MaxIJK[i] >= 3); // if less than 3 points, can't do numerical gradients (if necessary)

	VolInfo.AddOnID = AddOnID;
	ZoneXYZVarGetMinMax_Ordered3DZone(XYZVarNums, VolZoneNum, VolInfo.MinXYZ, VolInfo.MaxXYZ);
	VolInfo.DelXYZ = GetDelXYZ_Ordered3DZone(XYZVarNums, VolZoneNum);
	VolInfo.IsPeriodic = IsPeriodic;

	return TRUE;
}


const Boolean_t FindCritPoints(const int & VolZoneNum,
								const vector<int> & XYZVarNums,
								const int & RhoVarNum,
								const vector<int> & GradVarNums,
								const vector<int> & HessVarNums,
								const Boolean_t & IsPeriodic)
{
	TecUtilLockStart(AddOnID);

	FieldDataPointer_c RhoPtr;
	vector<FieldDataPointer_c> GradPtrs, HessPtrs;

	TecUtilDataLoadBegin();

	TecUtilPleaseWait("Reading data into memory...", TRUE);

	if (!GetReadPtrsForZone(VolZoneNum, 
		RhoVarNum, GradVarNums, HessVarNums,
		RhoPtr, GradPtrs, HessPtrs)){
		TecUtilDialogErrMsg("Failed to get read pointer(s)");
		return FALSE;
	}

	VolExtentIndexWeights_s VolInfo;

	if (!GetVolInfo(VolZoneNum, XYZVarNums, IsPeriodic, VolInfo)){
		TecUtilDialogErrMsg("Failed to get volume zone info");
		return FALSE;
	}

	CritPoints_c VolCPs;

	TecUtilPleaseWait(NULL, FALSE);

	if (FindCPs(VolCPs, VolInfo, DefaultCellSpacing, DefaultRhoCutoff, IsPeriodic, RhoPtr, GradPtrs, HessPtrs)){
		VolCPs.SaveAsOrderedZone(XYZVarNums, RhoVarNum, TRUE);
	}

	TecUtilDataLoadEnd();

	TecUtilLockFinish(AddOnID);

	TecUtilDialogMessageBox("Finished", MessageBoxType_Information);

	return TRUE;
}



/*
 *	Create bond paths/ring lines from a source CP zone(s)
 */
const Boolean_t FindBondRingLines(const int & VolZoneNum,
	const int & CPZoneNum,
	const char & CPType,
	const vector<int> & XYZVarNums,
	const int & RhoVarNum,
	const vector<int> & GradVarNums,
	const vector<int> & HessVarNums,
	const Boolean_t & IsPeriodic)
{
	TecUtilLockStart(AddOnID);

	int NumZones = TecUtilDataSetGetNumZones();
	int NumVars = TecUtilDataSetGetNumVars();

	REQUIRE(CPZoneNum > 0 && CPZoneNum <= NumZones);
	REQUIRE(XYZVarNums.size() == 3);
	for (const auto & i : XYZVarNums) REQUIRE(i > 0 && i <= NumVars);

	int TypeInd = std::find(CPTypeList, CPTypeList + 6, CPType) - CPTypeList;
	if (TypeInd >= 6){
		TecUtilDialogErrMsg("Invalid CP type specified");
		return FALSE;
	}

	int CPTypeVarNum = VarNumByName(CPTypeVarName);
	if (CPTypeVarNum < 0){
		TecUtilDialogErrMsg("Failed to find CP Type variable");
		return FALSE;
	}

	FieldDataPointer_c RhoPtr;
	vector<FieldDataPointer_c> GradPtrs, HessPtrs;

	TecUtilDataLoadBegin();

	if (!GetReadPtrsForZone(VolZoneNum,
		RhoVarNum, GradVarNums, HessVarNums,
		RhoPtr, GradPtrs, HessPtrs)){
		TecUtilDialogErrMsg("Failed to get read pointer(s)");
		return FALSE;
	}

	VolExtentIndexWeights_s VolInfo;

	if (!GetVolInfo(VolZoneNum, XYZVarNums, IsPeriodic, VolInfo)){
		TecUtilDialogErrMsg("Failed to get volume zone info");
		return FALSE;
	}

	vector<int> iJunk(2);
	int NumCPs;

	TecUtilZoneGetIJK(CPZoneNum, &NumCPs, &iJunk[0], &iJunk[1]);
	for (const auto & i : iJunk) if (i > 1){
		TecUtilDialogErrMsg("CP zone is not i-ordered");
		return FALSE;
	}

	vec3 EigVals;
	mat33 EigVecs;

	MultiRootParams_s MR;
	MR.CalcType = GPType_Classic;
	MR.VolInfo = &VolInfo;
	MR.IsPeriodic = IsPeriodic;
	MR.HasGrad = (GradPtrs.size() == 3);
	MR.HasHess = (HessPtrs.size() == 6);
	MR.RhoPtr = &RhoPtr;
	MR.GradPtrs = &GradPtrs;
	MR.HessPtrs = &HessPtrs;
	vector<vec3> BV(3);
	BV[0] << 1. << 0. << 0.;
	BV[1] << 0. << 1. << 0.;
	BV[2] << 0. << 0. << 1.;
	MR.BasisVectors = &BV;

	int PrincipalDirInd = 2;
	int EndCPNumforName = 0;
	StreamDir_e GPDir = StreamDir_Forward;
	ColorIndex_t PathColor = Black_C;
	if (CPType == RINGCP){
		PrincipalDirInd = 0;
		GPDir = StreamDir_Reverse;
		PathColor = Green_C;
		EndCPNumforName++;
	}

	CritPoints_c CPs(CPZoneNum, XYZVarNums, RhoVarNum, CPTypeVarNum);

	vec3 StartPoint;

	const double StartPointOffset = 0.1 * CPs.GetMinCPDist();
	const int NumGPPts = 100;
	double TermRadius = 0.2;
	double RhoCutoff = DefaultRhoCutoff;

// 	vector<GradPath_c> GPs(CPs.NumCPs(TypeInd) * 2);
	vector<GradPath_c> GPs;

// 	int gpNum = 0;
	for (int cpNum = 0; cpNum < CPs.NumCPs(TypeInd); ++cpNum){
		CalcEigenSystemForPoint(CPs.GetXYZ(TypeInd, cpNum), EigVals, EigVecs, MR);
		for (int d = -1; d < 2; d += 2){
			StartPoint = CPs.GetXYZ(TypeInd, cpNum) + normalise(EigVecs.col(PrincipalDirInd)) * StartPointOffset * static_cast<double>(d);

// 			GPs[gpNum].SetStartEndCPNum(CPs.GetTotOffsetFromTypeNumOffset(TypeInd, cpNum), 0);
// 			GPs[gpNum++].SetupGradPath(StartPoint, GPDir, NumGPPts, GPType_Classic, GPTerminate_AtCP, NULL, &CPs, &TermRadius, &RhoCutoff, VolInfo, HessPtrs, GradPtrs, RhoPtr);

			GPs.push_back(GradPath_c(StartPoint, GPDir, NumGPPts, GPType_Classic, GPTerminate_AtCP, NULL, &CPs, &TermRadius, &RhoCutoff, VolInfo, HessPtrs, GradPtrs, RhoPtr));
			GPs.back().SetStartEndCPNum(CPs.GetTotOffsetFromTypeNumOffset(TypeInd, cpNum), 0);

		}
	}

	int NumGPs = GPs.size();

#ifndef _DEBUG
#pragma omp parallel for schedule(dynamic)
#endif
	for (int iGP = 0; iGP < NumGPs; ++iGP){
		GPs[iGP].Seed(false);
		GPs[iGP].PointPrepend(CPs.GetXYZ(GPs[iGP].GetStartEndCPNum()[0]), CPs.GetRho(GPs[iGP].GetStartEndCPNum()[0]));
		GPs[iGP].Resample(NumGPPts);
		if (CPType == RINGCP) GPs[iGP].Reverse();
	}

	for (int iGP = 0; iGP < NumGPs; iGP += 2){
		GradPath_c GP = GPs[iGP];
		int OldStartCPNum = GP.GetStartEndCPNum()[EndCPNumforName];
		GP += GPs[iGP + 1];
		vector<int> StartEndCPNums = GP.GetStartEndCPNum();
		vector<vector<int> > StartEndCPTypeAndOffset(2);
		StartEndCPTypeAndOffset[0] = CPs.GetTypeNumOffsetFromTotOffset(OldStartCPNum);
		string Name = "GP: " + CPNameList[StartEndCPTypeAndOffset[0][0]] + " " + to_string(StartEndCPTypeAndOffset[0][1] + 1) + " (";
		for (int i = 0; i < 2; ++i){
			if (StartEndCPNums[i] >= 0){
				StartEndCPTypeAndOffset[i] = CPs.GetTypeNumOffsetFromTotOffset(StartEndCPNums[i]);
				Name += CPNameList[StartEndCPTypeAndOffset[i][0]][0] + to_string(StartEndCPTypeAndOffset[i][1] + 1);
			}
			else Name += "FF";

			if (i == 0) Name += "-";
		}

		Name += ")";

		GP.SaveAsOrderedZone(Name, vector<FieldDataType_e>(), XYZVarNums, RhoVarNum, TRUE, PathColor);
		for (int i = 0; i < 2; ++i){
			AuxDataZoneSetItem(GP.GetZoneNum(), CSMAuxData.CC.GPEndNums[i], to_string(StartEndCPNums[i] + 1));
			if (StartEndCPNums[i] >= 0) AuxDataZoneSetItem(GP.GetZoneNum(), CSMAuxData.CC.GPEndTypes[i], CPNameList[StartEndCPTypeAndOffset[i][0]]);
			else AuxDataZoneSetItem(GP.GetZoneNum(), CSMAuxData.CC.GPEndTypes[i], "FF");
		}
		if (CPType == RINGCP) AuxDataZoneSetItem(GP.GetZoneNum(), CSMAuxData.CC.ZoneSubType, CSMAuxData.CC.ZoneSubTypeRingLine);
		else AuxDataZoneSetItem(GP.GetZoneNum(), CSMAuxData.CC.ZoneSubType, CSMAuxData.CC.ZoneSubTypeBondPath);
	}

	for (auto & GP : GPs){
		vector<int> StartEndCPNums = GP.GetStartEndCPNum();
		vector<vector<int> > StartEndCPTypeAndOffset(2);
		string Name = "GP: ";
		for (int i = 0; i < 2; ++i){
			if (StartEndCPNums[i] >= 0){
				StartEndCPTypeAndOffset[i] = CPs.GetTypeNumOffsetFromTotOffset(StartEndCPNums[i]);
				Name += CPNameList[StartEndCPTypeAndOffset[i][0]] + " " + to_string(StartEndCPTypeAndOffset[i][1] + 1);
			}
			else Name += "FF";

			if (i == 0) Name += " to ";
		}

		GP.SaveAsOrderedZone(Name, vector<FieldDataType_e>(), XYZVarNums, RhoVarNum, FALSE, PathColor);
		for (int i = 0; i < 2; ++i){
			AuxDataZoneSetItem(GP.GetZoneNum(), CSMAuxData.CC.GPEndNums[i], to_string(StartEndCPNums[i] + 1));
			if (StartEndCPNums[i] >= 0) AuxDataZoneSetItem(GP.GetZoneNum(), CSMAuxData.CC.GPEndTypes[i], CPNameList[StartEndCPTypeAndOffset[i][0]]);
			else AuxDataZoneSetItem(GP.GetZoneNum(), CSMAuxData.CC.GPEndTypes[i], "FF");
		}
		if (CPType == RINGCP) AuxDataZoneSetItem(GP.GetZoneNum(), CSMAuxData.CC.ZoneSubType, CSMAuxData.CC.ZoneSubTypeRingLine);
		else AuxDataZoneSetItem(GP.GetZoneNum(), CSMAuxData.CC.ZoneSubType, CSMAuxData.CC.ZoneSubTypeBondPath);
	}

	TecUtilDataLoadEnd();

	TecUtilLockFinish(AddOnID);

	TecUtilDialogMessageBox("Finished", MessageBoxType_Information);

	return TRUE;
}