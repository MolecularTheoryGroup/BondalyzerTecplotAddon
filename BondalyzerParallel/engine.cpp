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
#include <algorithm>

#include <omp.h>

// #include <windows.h>

#include "TECADDON.h"
#include "ADDGLBL.h"
#include "ADDONVER.h"
#include "GUIDEFS.h"
#include "ADKUTIL.h"

#include <gsl/gsl_errno.h>
#include <gsl/gsl_min.h>

#include "CSM_DATA_TYPES.h"
#include "CSM_DATA_SET_INFO.h"
#include "CSM_CRIT_POINTS.h"
#include "CSM_CALC_VARS.h"
#include "CSM_GRAD_PATH.h"
#include "CSM_FE_VOLUME.h"
#include "CSM_GUI.h"

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


struct MinFuncParams_GPLengthInPlane{
	StreamDir_e Direction = StreamDir_Invalid;
	int NumGPPoints = -1;
	GPType_e GPType = GPType_Invalid;
	GPTerminate_e GPTerminate = GPTerminate_Invalid;
	vec3 * TermPoint = NULL;
	CritPoints_c * CPs = NULL;
	double * TermPointRadius = NULL;
	double * TermValue = NULL;
	VolExtentIndexWeights_s * VolInfo = NULL;
	const vector<FieldDataPointer_c> * HessPtrs = NULL;
	const vector<FieldDataPointer_c> * GradPtrs = NULL;
	const FieldDataPointer_c * RhoPtr = NULL;

	GradPath_c GP;

	vec3 StartPointOrigin;
	vec3 RotVec;
	vec3 RotAxis;
	int StartCPNum = -1;
	int EndCPNum = -1;
};


const string DS_ZoneName_Delim = "_-_";
const string T41Ext = ".t41";

const static int NumCircleGPs = 50;
const static double SmallAngleFactor = 0.5;
const static int MaxIter_GPLengthInPlane = 100;

//DEBUG
#ifdef _DEBUG
vector<vector<double> > MinFunc_GPAngleLengths;
#endif

/*
 *	This is the convergence criteria for a one-dimensional GSL minimization
 *	based on gradient path length.
 */
const static double Tolerance_GPLengthInPlane = 1e-4;

/*
 *	If the terminal ends of two gradient paths diverge away from eachother
 *	such that the inner (dot) product of their terminal ends is less than
 *	below, they are considered divergent and to be terminating at different 
 *	far field CPs (if they were near field, then we would know which CP and
 *	don't need to resort to checking for divergence).
 */
const static double Tolerance_DivergentGPInnerProd = -0.9;

BondalyzerSteps_e CurrentCalcType;


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
				const int & NumNodesPerElem,
				vector<bool> & IsVisited){
	bool IsOk = true;
	IsVisited[NodeNum] = true;
	int NumElems = TecUtilDataNodeToElemMapGetNumElems(NodeToElemMap, NodeNum + 1);
	for (int e = 1; e <= NumElems && IsOk; ++e){
		int ei = TecUtilDataNodeToElemMapGetElem(NodeToElemMap, NodeNum + 1, e);
		for (int n = 0; n < NumNodesPerElem && IsOk; ++n){
			int ni = TecUtilDataNodeGetByRef(NodeMap, ei, n + 1) - 1;
			if (!IsVisited[ni]){
				IsVisited[ni] = true;
				IsOk = FEZoneDFS(ni, NodeMap, NodeToElemMap, NumNodesPerElem, IsVisited);
			}
		}
	}

	return IsOk;
}

// void FEZoneDFS(const int & NodeNum, const vector<vector<int> > & Elems, vector<bool> & IsVisited){
// 	IsVisited[NodeNum] = true;
// 	for ()
// }

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
			if (NodeNums.back() < 0 || NodeNums.back() >= IsoIJK[0]){
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

void GetAllClosedIsoSurfaces(){
	int IsoZoneNum = MAX(1, ZoneNumByName("Iso"));
	vector<int> NodeNums;

	int NumZones = TecUtilDataSetGetNumZones();

	char *UserInput;

	while (true){
		if (!TecUtilDialogGetSimpleText("Generate closed isosurfaces. Enter Isosurface zone number:", string(to_string(IsoZoneNum)).c_str(), &UserInput)){
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

void GetClosedIsoSurface(const int & IsoZoneNum, const vector<FieldDataPointer_c> & IsoReadPtrs, vector<int> & NodeNums){

	int IsoIJK[3];
	TecUtilZoneGetIJK(IsoZoneNum, &IsoIJK[0], &IsoIJK[1], &IsoIJK[2]); // {NumNodes, NumElems, NumNodesPerElem}

	ZoneType_e IsoZoneType = TecUtilZoneGetType(IsoZoneNum);

	char *IsoZoneName;
	TecUtilZoneGetName(IsoZoneNum, &IsoZoneName);

	Set_pa ZoneSet = TecUtilSetAlloc(TRUE);

	TecUtilDialogLaunchPercentDone("Generating isosurface subzone", FALSE);
	int PtNum = 1;

	NodeMap_pa NodeMap = TecUtilDataNodeGetReadableRef(IsoZoneNum);
	NodeToElemMap_pa NodeToElemMap = TecUtilDataNodeToElemMapGetReadableRef(IsoZoneNum);

	if (NodeNums.size() == 0){
		NodeNums.resize(IsoIJK[0]);
		for (int n = 0; n < IsoIJK[0]; ++n) NodeNums[n] = n;
	}

// 	 Construct internal representation of the IsoSurface.
// 	 Only need the connectivity information.
// 	 NodeConnectivity will be a list of connected nodes for each node
	vector<vector<int> > Elems(IsoIJK[1]);
	for (auto & e : Elems) e.reserve(IsoIJK[2]);
 
	int NodesPerElem = TecUtilDataNodeGetNodesPerElem(NodeMap);
	if (NodesPerElem != IsoIJK[2]) TecUtilDialogErrMsg("Number of nodes per element doesn't match!"); //sanity check
 
	// Get the full list of elements, each element a list of node numbers.
	// Everything is being stored in base-0.
	for (int e = 0; e < IsoIJK[1]; ++e) for (int n = 0; n < IsoIJK[2]; ++n) Elems[e].push_back(TecUtilDataNodeGetByRef(NodeMap, e + 1, n + 1) - 1);

	vector<bool> TotalNodeVisited(IsoIJK[0], false);

	int SubZoneNum = 1;

	for (const int & NodeNum : NodeNums){
		if (TotalNodeVisited[NodeNum]) continue;

		if (!TecUtilDialogCheckPercentDone(0)){
			TecUtilDialogDropPercentDone();
			return;
		}
		if (NodeNums.size() != IsoIJK[0])
			TecUtilDialogSetPercentDoneText(string("Generating isosurface subzone: " + to_string(PtNum++) + " of " + to_string(NodeNums.size())).c_str());
		else
			TecUtilDialogSetPercentDoneText(string("Generating isosurface subzone: " + to_string(PtNum++)).c_str());


		/*
		 * Now run a depth first search through the isosurface to get the connected component
		 * of the minimum distance node.
		 */
		vector<bool> NodeVisited(IsoIJK[0], false);

		if (!VALID_REF(NodeMap)){
			TecUtilDialogErrMsg("Isosurface node pointer(s) invalid. Quitting.");
			return;
		}

		FEZoneDFS(NodeNum, NodeMap, NodeToElemMap, IsoIJK[2], NodeVisited);

		for (int n = 0; n < NodeVisited.size(); ++n){
			TotalNodeVisited[n] = (TotalNodeVisited[n] || NodeVisited[n]);
		}

		/*
		 * Collect nodes/elements of single found connected component
		 */

		int NumNewNodes = 0;
		for (const auto & i : NodeVisited) NumNewNodes += (int)i;

		if (NumNewNodes <= 0){
			continue;
		}

		vector<int> NodeNumsNewToOld(NumNewNodes);
		vector<int> NodeNumsOldToNew(IsoIJK[0]);

		int NodeNumElems;
		int NewNodeNum = -1;

		for (int n = 0; n < IsoIJK[0]; ++n){
			if (NodeVisited[n]){
				NodeNumsNewToOld[++NewNodeNum] = n;
				NodeNumsOldToNew[n] = NewNodeNum;
			}
		}

		vector<FieldDataPointer_c> IsoWritePtrs(IsoReadPtrs.size());

		vector<vector<int> > NewElems;
		NewElems.reserve(IsoIJK[1]);
		vector<bool> ElemAdded(IsoIJK[1], false);

		for (int e = 0; e < IsoIJK[1]; ++e) if (!ElemAdded[e]) for (const auto & n : Elems[e]) if (NodeVisited[n])
		{
			ElemAdded[e] = true;
			NewElems.push_back(Elems[e]);
			for (auto & ne : NewElems.back()){
				ne = NodeNumsOldToNew[ne];
			}
			break;
		}

		/*
		 * Make new zone with single connected component
		 */

		if (!TecUtilDataSetAddZone(string(IsoZoneName + string(": Subzone ") + to_string(SubZoneNum++)).c_str(), NumNewNodes, NewElems.size(), IsoIJK[2], IsoZoneType, NULL)){
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

		for (int e = 0; e < NewElems.size(); ++e){
			for (int ei = 0; ei < IsoIJK[2]; ++ei){
				TecUtilDataNodeSetByRef(NewNodeMap, e + 1, ei + 1, NewElems[e][ei] + 1);
			}
		}

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

void BondalyzerReturnUserInfo(const bool GuiSuccess, const vector<GuiField_c> & Fields){
	if (!GuiSuccess) return;

	TecUtilLockStart(AddOnID);

	int VolZoneNum, RhoVarNum;
	vector<int> XYZVarNums(3), GradVarNums, HessVarNums;
	Boolean_t IsPeriodic;

	int fNum = 0;

	VolZoneNum = Fields[fNum++].GetReturnInt();
	IsPeriodic = Fields[fNum++].GetReturnBool();
	fNum++;
	for (int i = 0; i < 3; ++i) XYZVarNums[i] = Fields[fNum++].GetReturnInt();
	fNum++;
	RhoVarNum = Fields[fNum++].GetReturnInt();
	fNum++;


	if (Fields[fNum++].GetReturnBool()){
		GradVarNums.resize(3);
		for (int i = 0; i < 3; ++i) GradVarNums[i] = Fields[fNum++].GetReturnInt();
	}
	else fNum += 3;
	fNum++;

	if (Fields[fNum++].GetReturnBool()){
		HessVarNums.resize(6);
		for (int i = 0; i < 6; ++i) HessVarNums[i] = Fields[fNum++].GetReturnInt();
	}
	else fNum += 6;

	fNum++;

	if (CurrentCalcType == BondalyzerSteps_CriticalPoints){
		double CellSpacing = Fields[fNum].GetReturnDouble();
		FindCritPoints(VolZoneNum, XYZVarNums, RhoVarNum, GradVarNums, HessVarNums, IsPeriodic, CellSpacing);
	}
	else if (CurrentCalcType == BondalyzerSteps_BondPaths || CurrentCalcType == BondalyzerSteps_RingLines
		|| CurrentCalcType == BondalyzerSteps_InteratomicSurfaces || CurrentCalcType == BondalyzerSteps_RingSurfaces){

		int CPTypeVarNum = Fields[fNum++].GetReturnInt(), AllCPsZoneNum = Fields[fNum++].GetReturnInt(), SelectCPsZoneNum = Fields[fNum].GetReturnInt();
		vector<int> SelectedCPs = Fields[fNum].GetReturnIntVec();

		CPType_e CPType;
		if (CurrentCalcType == BondalyzerSteps_BondPaths || CurrentCalcType == BondalyzerSteps_InteratomicSurfaces) CPType = CPType_BondCP;
		else if (CurrentCalcType == BondalyzerSteps_RingLines || CurrentCalcType == BondalyzerSteps_RingSurfaces) CPType = CPType_RingCP;

		if (CurrentCalcType == BondalyzerSteps_BondPaths || CurrentCalcType == BondalyzerSteps_RingLines)
			FindBondRingLines(VolZoneNum, AllCPsZoneNum, SelectedCPs, CPTypeVarNum, CPType, XYZVarNums, RhoVarNum, GradVarNums, HessVarNums, IsPeriodic);
		else if (CurrentCalcType == BondalyzerSteps_InteratomicSurfaces || CurrentCalcType == BondalyzerSteps_RingSurfaces)
			FindBondRingSurfaces(VolZoneNum, AllCPsZoneNum, SelectedCPs, CPTypeVarNum, CPType, XYZVarNums, RhoVarNum, GradVarNums, HessVarNums, IsPeriodic);
	}

	TecUtilDialogMessageBox("Finished", MessageBoxType_Information);

	TecUtilLockFinish(AddOnID);
}

void BondalyzerGetUserInfo(BondalyzerSteps_e CalcType){
	// new implementation with CSM_Gui

	CurrentCalcType = CalcType;

	vector<GuiField_c> Fields = {
		GuiField_c(Gui_ZoneSelect, "Volume zone", "Full Volume"),
		GuiField_c(Gui_Toggle, "Periodic system"),
		GuiField_c(Gui_VertSep)
	};

	vector<string> XYZStr = { "X", "Y", "Z" }, HessXYZStr = { "XX", "XY", "XZ", "YY", "YZ", "ZZ" };

	for (int i = 0; i < 3; ++i) Fields.push_back(GuiField_c(Gui_VarSelect, XYZStr[i], XYZStr[i]));
	Fields.push_back(GuiField_c(Gui_VertSep));

	Fields.push_back(GuiField_c(Gui_VarSelect, "Electron density", CSMVarName.Dens));
	Fields.push_back(GuiField_c(Gui_VertSep));

	Fields.push_back(GuiField_c(Gui_ToggleEnable, "Density gradient vector variables present"));
	int iTmp = Fields.size() - 1;
	for (int i = 0; i < 3; ++i){
		Fields.push_back(GuiField_c(Gui_VarSelect, XYZStr[i], CSMVarName.DensGradVec[i]));
		Fields[iTmp].AppendSearchString(to_string(Fields.size() - 1));
		if (i < 2) Fields[iTmp].AppendSearchString(",");
	}
	Fields.push_back(GuiField_c(Gui_VertSep));

	Fields.push_back(GuiField_c(Gui_ToggleEnable, "Density Hessian variables present"));
	iTmp = Fields.size() - 1;
	for (int i = 0; i < HessXYZStr.size(); ++i){
		Fields.push_back(GuiField_c(Gui_VarSelect, HessXYZStr[i], CSMVarName.DensHessTensor[i].c_str()));
		Fields[iTmp].AppendSearchString(to_string(Fields.size() - 1));
		if (i < 5) Fields[iTmp].AppendSearchString(",");
	}

	int VolZoneNum, RhoVarNum;
	vector<int> XYZVarNums(3), GradVarNums, HessVarNums;
	Boolean_t IsPeriodic;

	string Title = "Find " + BondalyzerStepGUITitles[CalcType];

	if (CalcType == BondalyzerSteps_CriticalPoints){
		Fields.push_back(GuiField_c(Gui_VertSep));
		Fields.push_back(GuiField_c(Gui_Double, "CP search grid spacing", to_string(DefaultCellSpacing)));
	}
	else if (CalcType > BondalyzerSteps_CriticalPoints){

		Fields.push_back(GuiField_c(Gui_VertSep));
// 		Fields.push_back(GuiField_c(Gui_ZoneSelect, "Critical points zone", "Critical Points"));
		Fields.push_back(GuiField_c(Gui_VarSelect, "CP type variable", CSMVarName.CritPointType));
		string SearchString = "Critical Points";
		Fields.push_back(GuiField_c(Gui_ZoneSelect, "All critical points zone", SearchString));

		int CPTypeNum;
		if (CalcType == BondalyzerSteps_BondPaths || CalcType == BondalyzerSteps_InteratomicSurfaces) CPTypeNum = 1;
		else if (CalcType == BondalyzerSteps_RingLines || CalcType == BondalyzerSteps_RingSurfaces) CPTypeNum = 2;
		int CPZoneNum = ZoneNumByName(SearchString + ": " + CPNameList[CPTypeNum]);
		if (CPZoneNum > 0) SearchString += ": " + CPNameList[CPTypeNum];

		Fields.push_back(GuiField_c(Gui_ZonePointSelectMulti, "Source critical point(s)", SearchString));
	}

	CSMGui(Title, Fields, BondalyzerReturnUserInfo, AddOnID);
}

const int DeletePointsFromIOrderedZone(const int & ZoneNum, const vector<int> & DelPointNums){
	int NumZones = TecUtilDataSetGetNumZones();
	int NumVars = TecUtilDataSetGetNumVars();

	REQUIRE(ZoneNum > 0 && ZoneNum <= NumZones);

// 	if (DelPointNums.size() == 0){
// 		TecUtilZoneCopy(ZoneNum, 1, 0, 1, 1, 0, 1, 1, 0, 1);
// 		return TecUtilDataSetGetNumZones();
// 	}

	TecUtilDataLoadBegin();

	vector<FieldDataPointer_c> VarReadPtrs(NumVars), VarWritePtrs(NumVars);

	int IJK[3];
	TecUtilZoneGetIJK(ZoneNum, &IJK[0], &IJK[1], &IJK[2]);
	vector<FieldDataType_e> FDTypes(NumVars);

	for (int v = 0; v < NumVars; ++v){
		VarReadPtrs[v].GetReadPtr(ZoneNum, v + 1);
		FDTypes[v] = VarReadPtrs[v].GetDFType();
	}

	char* ZoneName;
	TecUtilZoneGetName(ZoneNum, &ZoneName);
	int NewZoneNum;

	if (TecUtilDataSetAddZone(ZoneName, IJK[0] - DelPointNums.size(), IJK[1], IJK[2], ZoneType_Ordered, FDTypes.data())){
		NewZoneNum = TecUtilDataSetGetNumZones();
		for (int v = 0; v < NumVars; ++v){
			VarWritePtrs[v].GetWritePtr(NewZoneNum, v + 1);
		}
		int NewPtNum = 0;
		for (int OldPtNum = 0; OldPtNum < IJK[0]; ++OldPtNum){
			if (!std::binary_search(DelPointNums.begin(), DelPointNums.end(), OldPtNum + 1)){
				for (int v = 0; v < NumVars; ++v) VarWritePtrs[v].Write(NewPtNum, VarReadPtrs[v][OldPtNum]);
				NewPtNum++;
			}
		}


	}
	else TecUtilDialogErrMsg("Failed to make new zone when deleting CPs");
	TecUtilStringDealloc(&ZoneName);
	TecUtilDataLoadEnd();

	return NewZoneNum;
}

void DeleteCPsReturnUserInfo(const bool GuiSuccess, const vector<GuiField_c> & Fields){
	if (!GuiSuccess) return;

	TecUtilLockStart(AddOnID);

	bool DeleteFromOtherCPZones = Fields[1].GetReturnBool();
	vector<int> CPNums = Fields[0].GetReturnIntVec();
	int ZoneNum = Fields[0].GetReturnInt();
	vector<int> OtherCPZoneNums;
	vector<vector<int> > OtherCPNums;

	if (DeleteFromOtherCPZones){
		TecUtilDataLoadBegin();

		vector<FieldData_pa> CPZoneXYZRefs, OtherCPZoneXYZRefs[3], OtherCPZoneTypeRefs;
		FieldData_pa CPZoneCPTypeRef;
		int XYZVarNums[3] = { -1, -1, -1 };
		int CPTypeVarNum = VarNumByName(CSMVarName.CritPointType);
		TecUtilAxisGetVarAssignments(&XYZVarNums[0], &XYZVarNums[1], &XYZVarNums[2]);
		if (XYZVarNums[0] < 0 || XYZVarNums[1] < 0 || XYZVarNums[2] < 0) for (int i = 1; i <= 3; ++i) XYZVarNums[i - 1] = i;

		for (int i = 0; i < 3; ++i) CPZoneXYZRefs.push_back(TecUtilDataValueGetReadableNativeRef(ZoneNum, XYZVarNums[i]));
		CPZoneCPTypeRef = TecUtilDataValueGetReadableNativeRef(ZoneNum, CPTypeVarNum);

		for (int z = 1; z <= TecUtilDataSetGetNumZones(); ++z) if (z != ZoneNum 
			&& TecUtilZoneIsOrdered(z) 
			&& AuxDataZoneItemMatches(z, CSMAuxData.CC.ZoneType, CSMAuxData.CC.ZoneTypeCPs))
		{
			OtherCPZoneNums.push_back(z);
			OtherCPNums.push_back(vector<int>());
			OtherCPZoneTypeRefs.push_back(TecUtilDataValueGetReadableNativeRef(z, CPTypeVarNum));
			for (int i = 0; i < 3; ++i) OtherCPZoneXYZRefs[i].push_back(TecUtilDataValueGetReadableNativeRef(z, XYZVarNums[i]));
		}

		for (const int & CP : CPNums){
			vec3 CPPos;
			for (int i = 0; i < 3; ++i) CPPos[i] = TecUtilDataValueGetByRef(CPZoneXYZRefs[i], CP);
			int CPType = TecUtilDataValueGetByRef(CPZoneCPTypeRef, CP);

			for (int z = 0; z < OtherCPZoneNums.size(); ++z){
				int IJK[3];
				TecUtilZoneGetIJK(OtherCPZoneNums[z], &IJK[0], &IJK[1], &IJK[2]);
				for (int oCP = 1; oCP <= IJK[0]; ++oCP){
					if (TecUtilDataValueGetByRef(OtherCPZoneTypeRefs[z], oCP) == CPType){
						bool IsMatch = true;
						for (int i = 0; i < 3 && IsMatch; ++i) IsMatch = (TecUtilDataValueGetByRef(OtherCPZoneXYZRefs[i][z], oCP) == CPPos[i]);
						if (IsMatch) OtherCPNums[z].push_back(oCP);
					}

				}
			}
		}

		TecUtilDataLoadEnd();
	}

	Set_pa ZoneNumsSet = TecUtilSetAlloc(TRUE);
	vector<int> NewZoneNums, OldZoneNums;


	if (DeleteFromOtherCPZones){
		bool ZoneRun = false;
		for (int z = 0; z < OtherCPZoneNums.size(); ++z){
			if (!ZoneRun && OtherCPZoneNums[z] > ZoneNum){
				ZoneRun = true;
				int NewZoneNum = DeletePointsFromIOrderedZone(ZoneNum, CPNums);
				if (NewZoneNum > 0){
					TecUtilSetAddMember(ZoneNumsSet, NewZoneNum, TRUE);
					NewZoneNums.push_back(NewZoneNum);
					OldZoneNums.push_back(ZoneNum);
				}
			}
			int NewZoneNum = DeletePointsFromIOrderedZone(OtherCPZoneNums[z], OtherCPNums[z]);
			if (NewZoneNum > 0){
				TecUtilSetAddMember(ZoneNumsSet, NewZoneNum, TRUE);
				NewZoneNums.push_back(NewZoneNum);
				OldZoneNums.push_back(OtherCPZoneNums[z]);
			}
		}
		if (!ZoneRun){
			int NewZoneNum = DeletePointsFromIOrderedZone(ZoneNum, CPNums);
			if (NewZoneNum > 0){
				TecUtilSetAddMember(ZoneNumsSet, NewZoneNum, TRUE);
				NewZoneNums.push_back(NewZoneNum);
				OldZoneNums.push_back(ZoneNum);
			}
		}
	}
	else{
		int NewZoneNum = DeletePointsFromIOrderedZone(ZoneNum, CPNums);
		if (NewZoneNum > 0){
			TecUtilSetAddMember(ZoneNumsSet, NewZoneNum, TRUE);
			NewZoneNums.push_back(NewZoneNum);
			OldZoneNums.push_back(ZoneNum);
		}
	}

	TecUtilStateChanged(StateChange_ZonesAdded, (ArbParam_t)ZoneNumsSet);

	TecUtilZoneSetScatter(SV_SHOW, ZoneNumsSet, 0.0, TRUE);
	TecUtilZoneSetScatter(SV_FRAMESIZE, ZoneNumsSet, 1, FALSE);
	TecUtilZoneSetScatterSymbolShape(SV_GEOMSHAPE, ZoneNumsSet, GeomShape_Sphere);
	TecUtilZoneSetMesh(SV_SHOW, ZoneNumsSet, 0.0, FALSE);
	TecUtilZoneSetContour(SV_SHOW, ZoneNumsSet, 0.0, FALSE);
	TecUtilZoneSetShade(SV_SHOW, ZoneNumsSet, 0.0, FALSE);
	TecUtilZoneSetVector(SV_SHOW, ZoneNumsSet, 0.0, FALSE);

	TecUtilZoneSetActive(ZoneNumsSet, AssignOp_PlusEquals);

	int TmpZoneNum;
	TecUtilSetForEachMember(TmpZoneNum, ZoneNumsSet){ SetCPZone(TmpZoneNum); }

	TecUtilSetClear(ZoneNumsSet);
	TecUtilSetAddMember(ZoneNumsSet, ZoneNum, TRUE);
	for (const int & z : OtherCPZoneNums) TecUtilSetAddMember(ZoneNumsSet, z, TRUE);

	TecUtilDataSetDeleteZone(ZoneNumsSet);

	TecUtilSetDealloc(&ZoneNumsSet);
	TecUtilLockFinish(AddOnID);
}

void DeleteCPsGetUserInfo(){
	vector<GuiField_c> Fields = {
		GuiField_c(Gui_ZonePointSelectMulti, "Critical point(s)", "Critical Points"),
		GuiField_c(Gui_Toggle, "Also delete from other CP zone(s)", "1")
	};

	CSMGui("Delete critical point(s)", Fields, DeleteCPsReturnUserInfo, AddOnID);
}

const Boolean_t GetReadPtrsForZone(const int & ZoneNum,
	const int & RhoVarNum,
	const vector<int> & GradVarNums,
	const vector<int> & HessVarNums,
	FieldDataPointer_c & RhoPtr,
	vector<FieldDataPointer_c> & GradPtrs,
	vector<FieldDataPointer_c> & HessPtrs)
{
	TecUtilPleaseWait("Loading data", TRUE);

	Boolean_t IsOk = TRUE;

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
		IsOk = FALSE;
	}
	for (int i = 0; i < GradVarNums.size() && IsOk; ++i) if (!GradPtrs[i].GetReadPtr(ZoneNum, GradVarNums[i])){
		TecUtilDialogErrMsg(string("Failed to get gradient read pointer (zone " + to_string(ZoneNum) + ", var " + to_string(GradVarNums[i]) + ")").c_str());
		IsOk = FALSE;
	}
	for (int i = 0; i < HessVarNums.size() && IsOk; ++i) if (!HessPtrs[i].GetReadPtr(ZoneNum, HessVarNums[i])){
		TecUtilDialogErrMsg(string("Failed to get hessian read pointer (zone " + to_string(ZoneNum) + ", var " + to_string(HessVarNums[i]) + ")").c_str());
		IsOk = FALSE;
	}

	TecUtilPleaseWait("Loading data", FALSE);

	return IsOk;
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
								const Boolean_t & IsPeriodic,
								const double & CellSpacing)
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

	if (FindCPs(VolCPs, VolInfo, CellSpacing, DefaultRhoCutoff, IsPeriodic, RhoPtr, GradPtrs, HessPtrs)){
		VolCPs.SaveAsOrderedZone(XYZVarNums, RhoVarNum, TRUE);
	}

	TecUtilDataLoadEnd();

	TecUtilLockFinish(AddOnID);

	return TRUE;
}



/*
 *	Create bond paths/ring lines from a source CP zone(s)
 */
const Boolean_t FindBondRingLines(const int & VolZoneNum,
	const int & AllCPsZoneNum,
	const vector<int> & SelectedCPNums,
	const int & CPTypeVarNum,
	const CPType_e & CPType,
	const vector<int> & XYZVarNums,
	const int & RhoVarNum,
	const vector<int> & GradVarNums,
	const vector<int> & HessVarNums,
	const Boolean_t & IsPeriodic)
{
	TecUtilLockStart(AddOnID);

	int NumZones = TecUtilDataSetGetNumZones();
	int NumVars = TecUtilDataSetGetNumVars();

	REQUIRE(AllCPsZoneNum > 0 && AllCPsZoneNum <= NumZones);
	REQUIRE(XYZVarNums.size() == 3);
	for (const auto & i : XYZVarNums) REQUIRE(i > 0 && i <= NumVars);

	int TypeInd = std::find(CPTypeList, CPTypeList + 6, CPType) - CPTypeList;
	if (TypeInd >= 6){
		TecUtilDialogErrMsg("Invalid CP type specified");
		return FALSE;
	}

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

	TecUtilZoneGetIJK(AllCPsZoneNum, &NumCPs, &iJunk[0], &iJunk[1]);
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

	int EndCPNumforName = 0;
	StreamDir_e GPDir = StreamDir_Forward;
	ColorIndex_t PathColor = Black_C;
	if (CPType == CPType_RingCP){
		GPDir = StreamDir_Reverse;
		PathColor = Green_C;
		EndCPNumforName++;
	}


	CritPoints_c CPs(AllCPsZoneNum, XYZVarNums, RhoVarNum, CPTypeVarNum, &MR);

	if (SelectedCPNums.size() > CPs.NumCPs(TypeInd) || SelectedCPNums.back() - 1 > CPs.NumCPs(TypeInd)){
		TecUtilDialogErrMsg("Discrepancy between selected CPs and selected CP zone");
		return FALSE;
	}

	vec3 StartPoint;

	const double StartPointOffset = 0.1 * CPs.GetMinCPDist();
	const int NumGPPts = 100;
	double TermRadius = 0.2;
	double RhoCutoff = DefaultRhoCutoff;

	vector<GradPath_c> GPs;

	for (int cpNum = 0; cpNum < SelectedCPNums.size(); ++cpNum){
		for (int d = -1; d < 2; d += 2){
			StartPoint = CPs.GetXYZ(TypeInd, SelectedCPNums[cpNum] - 1) + CPs.GetPrincDir(TypeInd, SelectedCPNums[cpNum] - 1) * StartPointOffset * static_cast<double>(d);

			GPs.push_back(GradPath_c(StartPoint, GPDir, NumGPPts, GPType_Classic, GPTerminate_AtCP, NULL, &CPs, &TermRadius, &RhoCutoff, VolInfo, HessPtrs, GradPtrs, RhoPtr));
			GPs.back().SetStartEndCPNum(CPs.GetTotOffsetFromTypeNumOffset(TypeInd, SelectedCPNums[cpNum] - 1), 0);

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
		if (CPType == CPType_RingCP) GPs[iGP].Reverse();
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
		AuxDataZoneSetItem(GP.GetZoneNum(), CSMAuxData.CC.ZoneType, CSMAuxData.CC.ZoneTypeGP);
		if (CPType == CPType_RingCP) AuxDataZoneSetItem(GP.GetZoneNum(), CSMAuxData.CC.ZoneSubType, CSMAuxData.CC.ZoneSubTypeRingLine);
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
		AuxDataZoneSetItem(GP.GetZoneNum(), CSMAuxData.CC.ZoneType, CSMAuxData.CC.ZoneTypeGP);
		if (CPType == CPType_RingCP) AuxDataZoneSetItem(GP.GetZoneNum(), CSMAuxData.CC.ZoneSubType, CSMAuxData.CC.ZoneSubTypeRingLine);
		else AuxDataZoneSetItem(GP.GetZoneNum(), CSMAuxData.CC.ZoneSubType, CSMAuxData.CC.ZoneSubTypeBondPath);
	}

	TecUtilDataLoadEnd();

	TecUtilLockFinish(AddOnID);

	return TRUE;
}

double MinFunc_GPLength_InPlane(double alpha, void * params){
	MinFuncParams_GPLengthInPlane * GPParams = reinterpret_cast<MinFuncParams_GPLengthInPlane*>(params);

	GPParams->GP = GradPath_c(
		GPParams->StartPointOrigin + Rotate(GPParams->RotVec, alpha, GPParams->RotAxis),
		GPParams->Direction,
		GPParams->NumGPPoints,
		GPParams->GPType,
		GPParams->GPTerminate,
		GPParams->TermPoint,
		GPParams->CPs,
		GPParams->TermPointRadius,
		GPParams->TermValue,
		*GPParams->VolInfo,
		*GPParams->HessPtrs,
		*GPParams->GradPtrs,
		*GPParams->RhoPtr);

	GPParams->GP.SetStartEndCPNum(GPParams->StartCPNum, 0);

	GPParams->GP.Seed(false);

	if (!GPParams->GP.IsMade())
		TecUtilDialogErrMsg("Gradient path in MinFunc_GPLength_InPlane failed");
	if (GPParams->EndCPNum >= 0 && GPParams->GP.GetStartEndCPNum(1) != GPParams->EndCPNum)
		TecUtilDialogErrMsg("Gradient path in MinFunc_GPLength_InPlane terminated at wrong CP");

	GPParams->GP.PointPrepend(GPParams->StartPointOrigin, 0);
#ifdef _DEBUG
	MinFunc_GPAngleLengths.push_back(vector<double>({ alpha, GPParams->GP.GetLength() }));
#endif

	return GPParams->GP.GetLength();
}

/*
*	Create interatomic (bond) and ring surfaces from a source CP zone(s)
*	Also creates bond-ring, and bond-cage (ring-nuclear) paths.
*/
const Boolean_t FindBondRingSurfaces(const int & VolZoneNum,
	const int & AllCPsZoneNum,
	const vector<int> & SelectedCPNums,
	const int & CPTypeVarNum,
	const CPType_e & CPType,
	const vector<int> & XYZVarNums,
	const int & RhoVarNum,
	const vector<int> & GradVarNums,
	const vector<int> & HessVarNums,
	const Boolean_t & IsPeriodic)
{
	TecUtilLockStart(AddOnID);

	int NumZones = TecUtilDataSetGetNumZones();
	int NumVars = TecUtilDataSetGetNumVars();

	REQUIRE(AllCPsZoneNum > 0 && AllCPsZoneNum <= NumZones);
	REQUIRE(XYZVarNums.size() == 3);
	for (const auto & i : XYZVarNums) REQUIRE(i > 0 && i <= NumVars);

	int TypeInd = std::find(CPTypeList, CPTypeList + 6, CPType) - CPTypeList;
	if (TypeInd >= 6){
		TecUtilDialogErrMsg("Invalid CP type specified");
		return FALSE;
	}

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

	TecUtilZoneGetIJK(AllCPsZoneNum, &NumCPs, &iJunk[0], &iJunk[1]);
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


	CritPoints_c CPs(AllCPsZoneNum, XYZVarNums, RhoVarNum, CPTypeVarNum, &MR);

	if (SelectedCPNums.size() > CPs.NumCPs(TypeInd) || SelectedCPNums.back() - 1 > CPs.NumCPs(TypeInd)){
		TecUtilDialogErrMsg("Discrepancy between selected CPs and selected CP zone");
		return FALSE;
	}

	int EndCPNumforName = 1;
	StreamDir_e GPDir = StreamDir_Reverse;
	ColorIndex_t PathColor = Black_C;
	CPType_e PrimaryTermType = CPType_CageCP,
		SecondaryTermType = CPType_RingCP;
	if (CPType == CPType_RingCP){
		GPDir = StreamDir_Forward;
		PathColor = Green_C;
		EndCPNumforName--;
		PrimaryTermType = CPType_NuclearCP;
		SecondaryTermType = CPType_BondCP;
	}

	vec3 StartPoint;

	const double StartPointOffset = 0.1 * CPs.GetMinCPDist();
	const int NumGPPts = 100;
	double TermRadius = 0.2;
	double RhoCutoff = DefaultRhoCutoff;

	NumCPs = SelectedCPNums.size();
// 	NumCPs = CPs.NumCPs(TypeInd);

	double AngleStep = 2. * PI / static_cast<double>(NumCircleGPs);

	/*
	 * Find the CP to CP connections by seeding GPs around each CP along a circle of constant radius
	 * in the plane normal to the CP principal direction.
	 * For ring (interatomic) surfaces:
	 * Most of these GPs will terminate at cage (nuclear) CPs, but will switch as we seed with an increasing
	 * angle around the CP either from one cage (nuclear) CP to another, or terminate at the intermediate bond ring (bond) CP.
	 * If the GPs transition directly from one cage (nuclear) CP to another, then take the half-angle between the two GPs to
	 * search for the intermediate ring (bond) CP.
	 * Use each GP terminating at a cage (nuclear) or ring (bond) CP, as the starting guess for a minimum length GP connecting
	 * the source CP to the terminating CP.
	 * This will be a one-dimensional minimization, for which GSL has (presumedly) fast algorithms.
	 */

	vector<vector<GradPath_c> > GPsPerCP(NumCPs), SGPsPerCP(NumCPs);

#ifndef _DEBUG
// #pragma omp parallel for schedule(dynamic)
#endif
	for (int iCP = 0; iCP < NumCPs; ++iCP){
		int cpNum = SelectedCPNums[iCP] - 1;
		MinFuncParams_GPLengthInPlane GPParams;
		GPParams.Direction = GPDir;
		GPParams.NumGPPoints = NumGPPts;
		GPParams.GPType = GPType_Classic;
		GPParams.GPTerminate = GPTerminate_AtCP;
		GPParams.TermPoint = NULL;
		GPParams.CPs = &CPs;
		GPParams.TermPointRadius = &TermRadius;
		GPParams.TermValue = &RhoCutoff;
		GPParams.VolInfo = &VolInfo;
		GPParams.HessPtrs = &HessPtrs;
		GPParams.GradPtrs = &GradPtrs;
		GPParams.RhoPtr = &RhoPtr;

		GPParams.StartPointOrigin = CPs.GetXYZ(TypeInd, cpNum);
		GPParams.RotVec = CPs.GetEigVecs(TypeInd, cpNum).col(1);
		GPParams.RotAxis = CPs.GetPrincDir(TypeInd, cpNum);

		GPParams.StartCPNum = CPs.GetTotOffsetFromTypeNumOffset(TypeInd, cpNum);

		vec3 StartVec = normalise(CPs.GetEigVecs(TypeInd, cpNum).col(1)) * StartPointOffset;
		vector<GradPath_c> GPs;
		vector<double> GPLengths(NumCircleGPs);
		vector<CPType_e> GPEndCPTypes(NumCircleGPs);
		vector<int> GPEndCPNums(NumCircleGPs);

		/*
		 *	Seed GPs around the CP
		 */
		for (int gpNum = 0; gpNum < NumCircleGPs; ++gpNum){
			double RotAngle = static_cast<double>(gpNum) * AngleStep;

			StartPoint = CPs.GetXYZ(TypeInd, cpNum) + Rotate(StartVec, RotAngle, CPs.GetPrincDir(TypeInd, cpNum));

			GPs.push_back(GradPath_c(StartPoint, GPDir, NumGPPts, GPType_Classic, GPTerminate_AtCP, NULL, &CPs, &TermRadius, &RhoCutoff, VolInfo, HessPtrs, GradPtrs, RhoPtr));
			GPs.back().SetStartEndCPNum(CPs.GetTotOffsetFromTypeNumOffset(TypeInd, cpNum), 0);
		}

#ifndef _DEBUG
// #pragma omp parallel for schedule(dynamic)
#endif
		for (int gpNum = 0; gpNum < NumCircleGPs; ++gpNum){
			GPs[gpNum].Seed(false);
			GPs[gpNum].PointPrepend(CPs.GetXYZ(TypeInd, cpNum), CPs.GetRho(TypeInd, cpNum));
			GPLengths[gpNum] = GPs[gpNum].GetLength();
			GPEndCPNums[gpNum] = GPs[gpNum].GetStartEndCPNum(1);
			GPEndCPTypes[gpNum] = CPs.GetTypeFromTotOffset(GPEndCPNums[gpNum]);
		}

		// DEBUG
// 		for (int gpNum = 0; gpNum < NumCircleGPs; ++gpNum){
// 			GradPath_c GP = GPs[gpNum];
// 			vector<int> StartEndCPNums = GP.GetStartEndCPNum();
// 			vector<vector<int> > StartEndCPTypeAndOffset(2);
// 			string Name = "GP: ";
// 			for (int i = 0; i < 2; ++i){
// 				if (StartEndCPNums[i] >= 0){
// 					StartEndCPTypeAndOffset[i] = CPs.GetTypeNumOffsetFromTotOffset(StartEndCPNums[i]);
// 					Name += CPNameList[StartEndCPTypeAndOffset[i][0]] + " " + to_string(StartEndCPTypeAndOffset[i][1] + 1);
// 				}
// 				else Name += "FF";
// 
// 				if (i == 0) Name += " to ";
// 			}
// 
// 			GP.SaveAsOrderedZone(Name, vector<FieldDataType_e>(), XYZVarNums, RhoVarNum, FALSE, PathColor);
// 			for (int i = 0; i < 2; ++i){
// 				AuxDataZoneSetItem(GP.GetZoneNum(), CSMAuxData.CC.GPEndNums[i], to_string(StartEndCPNums[i] + 1));
// 				if (StartEndCPNums[i] >= 0) AuxDataZoneSetItem(GP.GetZoneNum(), CSMAuxData.CC.GPEndTypes[i], CPNameList[StartEndCPTypeAndOffset[i][0]]);
// 				else AuxDataZoneSetItem(GP.GetZoneNum(), CSMAuxData.CC.GPEndTypes[i], "FF");
// 			}
// 			AuxDataZoneSetItem(GP.GetZoneNum(), CSMAuxData.CC.ZoneType, CSMAuxData.CC.ZoneTypeGP);
// 			AuxDataZoneSetItem(GP.GetZoneNum(), "Seed angle", to_string(static_cast<double>(gpNum)* AngleStep).c_str());
// 		}
// 		return TRUE;

		/*
		 *	Use the seeded gradient paths to get starting guesses for special gradient paths.
		 *	Start with gradient paths terminating at the PrimaryTermType (nuclear or cage).
		 *	1.	Find the minimum length GP terminating at a given CP
		 *	2.	Test a small angle different to the left/right of that GP to see in which direction
		 *		the minimum length GP will lie.
		 *	3.	Then we know which two GPs will bracket the minimum length GP, and can start the
		 *		one-dimensional minimization.
		 */


		vector<vec3> MinFunc_AlphaLowMidHigh;
		vector<int> TermCPNums;
		
		vector<bool> GPUsed(NumCircleGPs, false);
		for (int gpNum = 0; gpNum < NumCircleGPs; ++gpNum){
			if (!GPUsed[gpNum] || (GPEndCPTypes[gpNum] == PrimaryTermType
				&& GPEndCPTypes[(gpNum + 1) % NumCircleGPs] == PrimaryTermType
				&& GPEndCPNums[gpNum] != GPEndCPNums[(gpNum + 1) % NumCircleGPs])
				|| (GPEndCPTypes[gpNum] == CPType_Invalid
				&& GPEndCPTypes[(gpNum + 1) % NumCircleGPs] == CPType_Invalid
				&& dot(GPs[gpNum][-1] - GPs[gpNum][-2], GPs[(gpNum + 1) % NumCircleGPs][-1] - GPs[(gpNum + 1) % NumCircleGPs][-2]) < Tolerance_DivergentGPInnerProd)){
				int TermCPNum = -1;
				int MinGPNum;
				double GPLen = -1.;

				double AlphaLower, AlphaUpper, AlphaGuess;
				double AngleFactor;
				if ((GPEndCPTypes[gpNum] == PrimaryTermType
					&& GPEndCPTypes[(gpNum + 1) % NumCircleGPs] == PrimaryTermType
					&& GPEndCPNums[gpNum] != GPEndCPNums[(gpNum + 1) % NumCircleGPs])
					|| (GPEndCPTypes[gpNum] == CPType_Invalid
					&& GPEndCPTypes[(gpNum + 1) % NumCircleGPs] == CPType_Invalid
					&& dot(GPs[gpNum][-1] - GPs[gpNum][-2], GPs[(gpNum + 1) % NumCircleGPs][-1] - GPs[(gpNum + 1) % NumCircleGPs][-2]) < Tolerance_DivergentGPInnerProd)){
					/*
					*	The gpNum'th gradient path and the (gpNum + 1)'th gradient path diverge to two CPs.
					*	This means there is a (ring) saddle CP to be found between these gradient paths.
					*	We'll use a binary search to get at the CP sandwiched between these two.
					*/
					int CPNumLower, CPNumUpper;
					vec3 TermVecLower, TermVecUpper;

					if (GPEndCPTypes[gpNum] == PrimaryTermType){
						CPNumLower = GPEndCPNums[gpNum];
						CPNumUpper = GPEndCPNums[(gpNum + 1) % NumCircleGPs];
					}
					else{
						TermVecLower = GPs[gpNum][-1] - GPs[gpNum][-2];
						TermVecUpper = GPs[(gpNum + 1) % NumCircleGPs][-1] - GPs[(gpNum + 1) % NumCircleGPs][-2];
					}

					AlphaLower = static_cast<double>(gpNum)* AngleStep;
					AlphaUpper = static_cast<double>(gpNum + 1) * AngleStep;

					while (TermCPNum < 0){
						AlphaGuess = (AlphaLower + AlphaUpper) * 0.5;
						StartPoint = CPs.GetXYZ(TypeInd, cpNum) + Rotate(StartVec, AlphaGuess, CPs.GetPrincDir(TypeInd, cpNum));
						GradPath_c GP(StartPoint, GPDir, NumGPPts, GPType_Classic, GPTerminate_AtCP, NULL, &CPs, &TermRadius, &RhoCutoff, VolInfo, HessPtrs, GradPtrs, RhoPtr);
						GP.SetStartEndCPNum(CPs.GetTotOffsetFromTypeNumOffset(TypeInd, cpNum), 0);
						GP.Seed(false);
						if (CPs.GetTypeFromTotOffset(GP.GetStartEndCPNum(1)) == SecondaryTermType){
							TermCPNum = GP.GetStartEndCPNum(1);
							GP.PointPrepend(CPs.GetXYZ(TypeInd, cpNum), CPs.GetRho(TypeInd, cpNum));
							GPLen = GP.GetLength();
						}
						else{
							bool GoLeft;
							if (GPEndCPTypes[gpNum] == PrimaryTermType){
								if (GP.GetStartEndCPNum(1) == CPNumLower) GoLeft = false;
								else if (GP.GetStartEndCPNum(1) == CPNumUpper) GoLeft = true;
								else{
									TecUtilDialogErrMsg("GP terminated at unexpected CP when performing binary search.");
									break;
								}
							}
							else{
								vec3 TmpTermVec = GP[-1] - GP[-2];
								if (dot(TmpTermVec, TermVecLower) > 0.9) GoLeft = false;
								else if (dot(TmpTermVec, TermVecUpper) > 0.9) GoLeft = true;
								else{
									TecUtilDialogErrMsg("GP terminated not parallel to either bounding GPs when performing binary search.");
									break;
								}
							}
							if (GoLeft) AlphaUpper = AlphaGuess;
							else AlphaLower = AlphaGuess;
						}
					}
				}
				else if (GPEndCPTypes[gpNum] == PrimaryTermType || GPEndCPTypes[gpNum] == SecondaryTermType){
					/*
					*	Found a GP terminating at a CP, so now get minimum length GP
					*	terminating at the same CP
					*/
					TermCPNum = GPs[gpNum].GetStartEndCPNum(1);
					MinGPNum = gpNum;
					for (int i = gpNum; i < gpNum + NumCircleGPs; ++i){
						int ii = i % NumCircleGPs;
						if (!GPUsed[ii] && GPs[ii].GetStartEndCPNum(1) == TermCPNum){
							GPUsed[ii] = true;
							if (GPLengths[ii] < GPLengths[MinGPNum]){
								MinGPNum = ii;
							}
						}
					}

					AlphaGuess = static_cast<double>(MinGPNum) * AngleStep;
					GPLen = GPLengths[MinGPNum];
				}

				if (TermCPNum >= 0){
					/*
					*	We found a GP that terminates at the saddle CP, so now have to find the actual lower and upper
					*	bounds on the minimization.
					*	The following is modified from above: find direction of minimum, then step in that direction to find
					*	a better guess at the minimum.
					*/
					int MinGPDir; // will take the values of {-1, 1}
					vector<double> GPDirLengths(2, -1);
					/*
					*	Now find the angular direction in which the minimum length GP will lie.
					*/
					for (int i = 0; i < 2; ++i){
						MinGPDir = (i == 0) ? -1 : 1;
						double RotAngle;
						AngleFactor = 1.0;
						int Iter = 0;
						while (Iter < 50){
							Iter++;
							RotAngle = AlphaGuess + (AngleStep * AngleFactor * SmallAngleFactor * static_cast<double>(MinGPDir));
							StartPoint = CPs.GetXYZ(TypeInd, cpNum) + Rotate(StartVec, RotAngle, CPs.GetPrincDir(TypeInd, cpNum));
							GradPath_c GP(StartPoint, GPDir, NumGPPts, GPType_Classic, GPTerminate_AtCP, NULL, &CPs, &TermRadius, &RhoCutoff, VolInfo, HessPtrs, GradPtrs, RhoPtr);
							GP.SetStartEndCPNum(CPs.GetTotOffsetFromTypeNumOffset(TypeInd, cpNum), 0);
							GP.Seed(false);
							if (GP.GetStartEndCPNum(1) == TermCPNum){
								GP.PointPrepend(CPs.GetXYZ(TypeInd, cpNum), CPs.GetRho(TypeInd, cpNum));
								GPDirLengths[i] = GP.GetLength();
								break;
							}
							else{
								AngleFactor *= SmallAngleFactor;
							}
						}
						if (GPDirLengths[i] > 0){
							if (i == 0) AlphaLower = RotAngle;
							else AlphaUpper = RotAngle;
							if (GPDirLengths[i] < GPLen){
								break;
							}
						}
					}

					if (GPDirLengths[0] < 0 && GPDirLengths[1] < 0){
						TecUtilDialogErrMsg("Failed to find direction of minimum length GP.");
						continue;
					}

					if (GPDirLengths[0] > 0 && GPDirLengths[1] > 0 && GPDirLengths[0] > GPLen && GPDirLengths[1] > GPLen){
						/*
						*	This is the special case where the GPs seeded to either side of MinGPNum both
						*	terminated at the correct CP and had a greater length than GPLenghts[MinGPNum],
						*	meaning that they serve as the upper and lower bounds of the minimization.
						*	In this case the proper values are already saved to AlphaLower and AlphaUpper.
						*/
					}
					else{
						/*
						*	Now need to find a minimum length GP terminating at the saddle CP
						*	so we can get an upper and lower angular bound on the CP.
						*/
						if (MinGPDir < 0) AlphaUpper = AlphaGuess;
						else AlphaLower = AlphaGuess;
						double RotAngle = AlphaGuess;
						while (true){
							RotAngle += (AngleStep * AngleFactor * SmallAngleFactor * static_cast<double>(MinGPDir));
							StartPoint = CPs.GetXYZ(TypeInd, cpNum) + Rotate(StartVec, RotAngle, CPs.GetPrincDir(TypeInd, cpNum));
							GradPath_c GP(StartPoint, GPDir, NumGPPts, GPType_Classic, GPTerminate_AtCP, NULL, &CPs, &TermRadius, &RhoCutoff, VolInfo, HessPtrs, GradPtrs, RhoPtr);
							GP.SetStartEndCPNum(CPs.GetTotOffsetFromTypeNumOffset(TypeInd, cpNum), 0);
							GP.Seed(false);
							if (GP.GetStartEndCPNum(1) == TermCPNum){
								GP.PointPrepend(CPs.GetXYZ(TypeInd, cpNum), CPs.GetRho(TypeInd, cpNum));
								double TmpGPLen = GP.GetLength();
								if (TmpGPLen > GPLen){
									if (MinGPDir < 0){
										AlphaLower = RotAngle;
									}
									else{
										AlphaUpper = RotAngle;
									}
									break;
								}
								else AlphaGuess = RotAngle;
								GPLen = TmpGPLen;
							}
							else{
								AngleFactor *= SmallAngleFactor;
								RotAngle = AlphaGuess;
							}
						}


						AlphaGuess = (AlphaLower + AlphaUpper) * 0.5; // guess is the midpoint between the two angles
					}

					MinFunc_AlphaLowMidHigh.push_back(vec3());
					MinFunc_AlphaLowMidHigh.back() << AlphaLower << AlphaGuess << AlphaUpper;

					TermCPNums.push_back(TermCPNum);
				}
			}
		}

		if (MinFunc_AlphaLowMidHigh.size() != TermCPNums.size()){
			TecUtilDialogErrMsg("Number of starting angle values doesn't match number of terminating CPs");
			continue;
		}

		int NumSGPs = MinFunc_AlphaLowMidHigh.size();
		SGPsPerCP[iCP].resize(NumSGPs);
#ifndef _DEBUG
// #pragma omp parallel for schedule(dynamic)
#endif
		for (int i = 0; i < NumSGPs; ++i){
			/*
			*	Initialize and run the GSL one-dimensional minimization.
			*	https://www.gnu.org/software/gsl/manual/html_node/One-dimensional-Minimization.html#One-dimensional-Minimization
			*/

			MinFuncParams_GPLengthInPlane TmpGPParams = GPParams;
			TmpGPParams.EndCPNum = TermCPNums[i];

			//DEBUG
#ifdef _DEBUG
			MinFunc_GPAngleLengths.clear();
#endif

			gsl_function F;
			F.function = &MinFunc_GPLength_InPlane;
			F.params = reinterpret_cast<void*>(&TmpGPParams);

			const gsl_min_fminimizer_type *T;
			T = gsl_min_fminimizer_brent;

			gsl_min_fminimizer *s;
			s = gsl_min_fminimizer_alloc(T);

			// DEBUG check the GP lengths that will be found using the low, guess, and high angle values
#ifdef _DEBUG
			vector<double> TmpLens(3);
			for (int j = 0; j < 3; ++j){
				StartPoint = CPs.GetXYZ(TypeInd, cpNum) + Rotate(StartVec, MinFunc_AlphaLowMidHigh[i][j], CPs.GetPrincDir(TypeInd, cpNum));
				GradPath_c GP(StartPoint, GPDir, NumGPPts, GPType_Classic, GPTerminate_AtCP, NULL, &CPs, &TermRadius, &RhoCutoff, VolInfo, HessPtrs, GradPtrs, RhoPtr);

				GP.SetStartEndCPNum(CPs.GetTotOffsetFromTypeNumOffset(TypeInd, cpNum), 0);
				GP.Seed(false);
				GP.PointPrepend(CPs.GetXYZ(TypeInd, cpNum), CPs.GetRho(TypeInd, cpNum));

				TmpLens[j] = GP.GetLength();
			}
#endif

			gsl_min_fminimizer_set(s, &F, MinFunc_AlphaLowMidHigh[i][1], MinFunc_AlphaLowMidHigh[i][0], MinFunc_AlphaLowMidHigh[i][2]);

			int Status;
			int Iter = 0;

			double OldLen = reinterpret_cast<MinFuncParams_GPLengthInPlane*>(F.params)->GP.GetLength();
			double NewLen;

			/*
			*	Perform the minimization
			*/
			do
			{
				Iter++;
				Status = gsl_min_fminimizer_iterate(s);

				if (Status == GSL_EBADFUNC || Status == GSL_FAILURE){
					TecUtilDialogErrMsg("GSL one-dimensional minimization failed");
					break;
				}

// 				if (reinterpret_cast<MinFuncParams_GPLengthInPlane*>(F.params)->GP.GetStartEndCPNum(1) != TermCPNums[i]){
// 					TecUtilDialogErrMsg("Gradient path in GSL one-dimensional minimization terminated at wrong CP");
// 				}

// 				MinFunc_AlphaLowMidHigh[i][1] = gsl_min_fminimizer_x_minimum(s);
// 				MinFunc_AlphaLowMidHigh[i][0] = gsl_min_fminimizer_x_lower(s);
// 				MinFunc_AlphaLowMidHigh[i][2] = gsl_min_fminimizer_x_upper(s);
// 
// 				Status = gsl_min_test_interval(MinFunc_AlphaLowMidHigh[i][0], MinFunc_AlphaLowMidHigh[i][2], Tolerance_GPLengthInPlane, 0.0);

  				NewLen = reinterpret_cast<MinFuncParams_GPLengthInPlane*>(F.params)->GP.GetLength();
  
  				Status = (abs(OldLen - NewLen) < Tolerance_GPLengthInPlane) ? GSL_SUCCESS : GSL_CONTINUE;

				OldLen = NewLen;

			} while (Status == GSL_CONTINUE && Iter < MaxIter_GPLengthInPlane);

			gsl_min_fminimizer_free(s);

			if (Status != GSL_SUCCESS || Iter >= MaxIter_GPLengthInPlane){
				TecUtilDialogErrMsg("GSL one-dimensional minimization did not converge");
// 				continue;
			}

			/*
			*	Minimizating was successful, so save the resulting SGP
			*/
			SGPsPerCP[iCP][i] = reinterpret_cast<MinFuncParams_GPLengthInPlane*>(F.params)->GP;
// 			SGPsPerCP[iCP][i].SetupGradPath(GPParams.StartPointOrigin + Rotate(GPParams.RotVec, MinFunc_AlphaLowMidHigh[i][1], GPParams.RotAxis),
// 				GPParams.Direction,
// 				GPParams.NumGPPoints,
// 				GPParams.GPType,
// 				GPParams.GPTerminate,
// 				GPParams.TermPoint,
// 				GPParams.CPs,
// 				GPParams.TermPointRadius,
// 				GPParams.TermValue,
// 				*GPParams.VolInfo,
// 				*GPParams.HessPtrs,
// 				*GPParams.GradPtrs,
// 				*GPParams.RhoPtr);
// 			SGPsPerCP[iCP][i].SetStartEndCPNum(CPs.GetTotOffsetFromTypeNumOffset(TypeInd, cpNum), 0);
// 			SGPsPerCP[iCP][i].Seed();
// 			SGPsPerCP[iCP][i].PointPrepend(CPs.GetXYZ(TypeInd, cpNum), CPs.GetRho(TypeInd, cpNum));

			/*
			 */

			SGPsPerCP[iCP][i].Resample(NumGPPts);
			if (TypeInd == CPType_BondCP) SGPsPerCP[iCP][i].Reverse();
		}

	}

	/*
	 * Save special gradient paths
	 */

	for (int iCP = 0; iCP < NumCPs; ++iCP){
		for (auto & GP : SGPsPerCP[iCP]){
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
			AuxDataZoneSetItem(GP.GetZoneNum(), CSMAuxData.CC.ZoneType, CSMAuxData.CC.ZoneTypeGP);
		}
	}

	//DEBUG
	TecUtilDataLoadEnd();
	TecUtilLockFinish(AddOnID);
	return TRUE;


	vector<GradPath_c> GPs;

	/*
	 *	Shotgun approach.
	 *	First seed GPs around each CP along a constant radius circle in the 
	 *	plane normal to the CP principal direction.
	 *	These GPs will identify all bond-cage, ring-nuclear, and bond-ring connections.
	 */

	for (int cpNum = 0; cpNum < NumCPs; ++cpNum){
		vec3 StartVec = normalise(CPs.GetEigVecs(TypeInd, cpNum).col(1)) * StartPointOffset;

		for (int gpNum = 0; gpNum < NumCircleGPs; ++gpNum){
			double RotAngle = static_cast<double>(gpNum)* AngleStep;

			StartPoint = CPs.GetXYZ(TypeInd, cpNum) + Rotate(StartVec, RotAngle, CPs.GetPrincDir(TypeInd, cpNum));

			GPs.push_back(GradPath_c(StartPoint, GPDir, NumGPPts, GPType_Classic, GPTerminate_AtCP, NULL, &CPs, &TermRadius, &RhoCutoff, VolInfo, HessPtrs, GradPtrs, RhoPtr));
			GPs.back().SetStartEndCPNum(CPs.GetTotOffsetFromTypeNumOffset(TypeInd, cpNum), 0);

		}
// 		break;
	}

	int NumGPs = GPs.size();

#ifndef _DEBUG
#pragma omp parallel for schedule(dynamic)
#endif
	for (int iGP = 0; iGP < NumGPs; ++iGP){
		GPs[iGP].Seed(false);
		GPs[iGP].PointPrepend(CPs.GetXYZ(GPs[iGP].GetStartEndCPNum()[0]), CPs.GetRho(GPs[iGP].GetStartEndCPNum()[0]));
		GPs[iGP].Resample(NumGPPts);
		if (CPType == CPType_BondCP) GPs[iGP].Reverse();
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
		AuxDataZoneSetItem(GP.GetZoneNum(), CSMAuxData.CC.ZoneType, CSMAuxData.CC.ZoneTypeGP);
	}

	TecUtilDataLoadEnd();

	TecUtilLockFinish(AddOnID);

	return TRUE;
}