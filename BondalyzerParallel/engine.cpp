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
#include <cmath>

#include <omp.h>

// #include <windows.h>

#include "TECADDON.h"
#include "ADDGLBL.h"
#include "ADDONVER.h"
#include "GUIDEFS.h"
#include "ADKUTIL.h"

#include "Set.h"
#include "ArgList.h"
#include "StyleValue.h"

#include <gsl/gsl_errno.h>
#include <gsl/gsl_min.h>

#include "CSM_DATA_TYPES.h"
#include "CSM_DATA_SET_INFO.h"
#include "CSM_CRIT_POINTS.h"
#include "CSM_CALC_VARS.h"
#include "CSM_GRAD_PATH.h"
#include "CSM_FE_VOLUME.h"
#include "CSM_GUI.h"
#include "CSM_GEOMETRY.h"

#include "KFc.h"

#include "ENGINE.h"

#include <armadillo>
using namespace arma;
using namespace tecplot::toolbox;


using std::string;
using std::to_string;
using std::vector;
using std::stringstream;
using std::ofstream;
using std::ios;
using std::endl;
using std::setprecision;

#define PI2 PI * 2.

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

struct MinFuncParams_RCSFuncInPlane{
	const FieldDataPointer_c* RCSFuncPtr = NULL;
	VolExtentIndexWeights_s* VolInfo = NULL;
	vec3 StartPointOrigin;
	vec3 RotVec;
	vec3 RotAxis;
};


const string DS_ZoneName_Delim = "_-_";
const string T41Ext = ".t41";

const static int MinNumCircleCheckGPs = 360;
const static int MinNumCircleGPs = 120;
const static int MinNumRCSFuncCircleCheckPts = 3600;
const static int DefaultRCSNumCheckRadii = 128;
const static int DefaultRCSNumCheckRadiiConverged = 4;
const static double DefaultRCSCheckRadiiLow = 0.2;
const static double DefaultRCSCheckRadiiHigh = 0.8;
const static double DefaultRCSCheckRadiiHighMultiplier = 1.4; // if <= 1 near field SGPs are found, increase the radius for RCSFunc-based search
const static int DefaultRCSFuncCheckNumPts = 3; 
const static int DefaultRCSFuncConvergence = 32;
const static int DefaultRCSSGPAngleCheck = 15; // If RCS-based SGP is less than this angle [degrees] from an existing SGP it is discarded
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
const static double Tolerance_GPLengthInPlane = 1e-5;

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

	vector<FESurface_c> Vols;

	TecUtilDrawGraphics(FALSE);
	TecUtilInterfaceSuspend(TRUE);
	TecUtilWorkAreaSuspend(TRUE);

	StatusLaunch("Refining zones", AddOnID, TRUE);

	for (int i = 1; i <= TecUtilDataSetGetNumZones(); ++i){
		FESurface_c Vol(i, vector<int>({ 1, 2, 3 }));
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
			UserQuit = !StatusUpdate(++VolNum, NumVols, "Refining zones", AddOnID);
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

	StatusDrop(AddOnID);

	TecUtilDrawGraphics(TRUE);
	TecUtilInterfaceSuspend(FALSE);
	TecUtilWorkAreaSuspend(FALSE);
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

		vector<string> Answers = SplitString(string(UserInput), ",");
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

		vector<string> Answers = SplitString(string(UserInput), ",");
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

		vector<string> Answers = SplitString(string(UserInput), ",");
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

		vector<string> Answers = SplitString(string(UserInput), ",");
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

		vector<string> Answers = SplitString(string(UserInput), ",");
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

	string StatusStr = "Generating isosurface subzone";
	StatusLaunch(StatusStr, AddOnID, FALSE);
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
		if (!StatusUpdate(0, 1, StatusStr + ": " + to_string(PtNum++) + (NodeNums.size() != IsoIJK[0] ? " of " + to_string(NodeNums.size()) : ""), AddOnID)){
			StatusDrop(AddOnID);
			return;
		}


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
	StatusDrop(AddOnID);
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
	else if (CurrentCalcType >= BondalyzerSteps_BondPaths){

		int RidgeFuncVarNum;

		if (CurrentCalcType >= BondalyzerSteps_InteratomicSurfaces){
			if (Fields[fNum++].GetReturnBool())
				RidgeFuncVarNum = Fields[fNum++].GetReturnInt();
			else fNum++;
		}

		int CPTypeVarNum = Fields[fNum++].GetReturnInt(), 
			AllCPsZoneNum = Fields[fNum++].GetReturnInt(), 
			SelectCPsZoneNum = Fields[fNum].GetReturnInt();
		vector<int> SelectedCPs = Fields[fNum].GetReturnIntVec();

		CPType_e CPType;
		if (CurrentCalcType == BondalyzerSteps_BondPaths || CurrentCalcType == BondalyzerSteps_InteratomicSurfaces) CPType = CPType_BondCP;
		else if (CurrentCalcType == BondalyzerSteps_RingLines || CurrentCalcType == BondalyzerSteps_RingSurfaces) CPType = CPType_RingCP;

		if (CurrentCalcType == BondalyzerSteps_BondPaths || CurrentCalcType == BondalyzerSteps_RingLines)
			FindBondRingLines(VolZoneNum, AllCPsZoneNum, SelectedCPs, CPTypeVarNum, CPType, XYZVarNums, RhoVarNum, GradVarNums, HessVarNums, IsPeriodic);
		else if (CurrentCalcType == BondalyzerSteps_InteratomicSurfaces || CurrentCalcType == BondalyzerSteps_RingSurfaces)
			FindBondRingSurfaces(VolZoneNum, AllCPsZoneNum, SelectedCPs, CPTypeVarNum, CPType, XYZVarNums, RhoVarNum, GradVarNums, HessVarNums, RidgeFuncVarNum, IsPeriodic);
	}

	TecUtilDialogMessageBox("Finished", MessageBoxType_Information);

	TecUtilLockFinish(AddOnID);
}

void BondalyzerGetUserInfo(BondalyzerSteps_e CalcType){
	// new implementation with CSM_Gui

	CurrentCalcType = CalcType;

	vector<GuiField_c> Fields = {
		GuiField_c(Gui_ZoneSelect, "Volume zone", CSMZoneName.FullVolume.substr(0,10)),
		GuiField_c(Gui_Toggle, "Periodic system"),
		GuiField_c(Gui_VertSep)
	};

	vector<string> XYZStr = { "X", "Y", "Z" }, HessXYZStr = { "XX", "XY", "XZ", "YY", "YZ", "ZZ" };

	for (int i = 0; i < 3; ++i) Fields.push_back(GuiField_c(Gui_VarSelect, XYZStr[i], XYZStr[i]));
	Fields.push_back(GuiField_c(Gui_VertSep));

	Fields.push_back(GuiField_c(Gui_VarSelect, "Electron density", CSMVarName.Dens));
	Fields.push_back(GuiField_c(Gui_VertSep));

	int iTmp = Fields.size();
	Fields.push_back(GuiField_c(Gui_ToggleEnable, "Density gradient vector variables present"));
	for (int i = 0; i < 3; ++i){
		Fields[iTmp].AppendSearchString(to_string(Fields.size()));
		Fields.push_back(GuiField_c(Gui_VarSelect, XYZStr[i], CSMVarName.DensGradVec[i]));
		if (i < 2) Fields[iTmp].AppendSearchString(",");
	}
	Fields.push_back(GuiField_c(Gui_VertSep));

	iTmp = Fields.size();
	Fields.push_back(GuiField_c(Gui_ToggleEnable, "Density Hessian variables present"));
	for (int i = 0; i < HessXYZStr.size(); ++i){
		Fields[iTmp].AppendSearchString(to_string(Fields.size()));
		Fields.push_back(GuiField_c(Gui_VarSelect, HessXYZStr[i], CSMVarName.DensHessTensor[i]));
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
	else if (CalcType >= BondalyzerSteps_BondPaths){
		Fields.push_back(GuiField_c(Gui_VertSep));
		if (CalcType >= BondalyzerSteps_InteratomicSurfaces){
			Fields.push_back(GuiField_c(Gui_ToggleEnable, "1-RCS function present", to_string(Fields.size() + 1)));
			Fields.push_back(GuiField_c(Gui_VarSelect, "1-ridge function variable", CSMVarName.EberlyFunc[0]));
		}

// 		Fields.push_back(GuiField_c(Gui_ZoneSelect, "Critical points zone", "Critical Points"));
		Fields.push_back(GuiField_c(Gui_VarSelect, "CP type variable", CSMVarName.CritPointType));
		string SearchString = CSMZoneName.CriticalPoints;
		Fields.push_back(GuiField_c(Gui_ZoneSelect, "All critical points zone", SearchString));

		int CPTypeNum;
		if (CalcType == BondalyzerSteps_BondPaths || CalcType == BondalyzerSteps_InteratomicSurfaces) CPTypeNum = 1;
		else if (CalcType == BondalyzerSteps_RingLines || CalcType == BondalyzerSteps_RingSurfaces) CPTypeNum = 2;
		int CPZoneNum = ZoneNumByName(CSMZoneName.CPType[CPTypeNum]);
		if (CPZoneNum > 0) SearchString = CSMZoneName.CPType[CPTypeNum];

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
		FDTypes[v] = VarReadPtrs[v].FDType();
	}

	char* ZoneName;
	TecUtilZoneGetName(ZoneNum, &ZoneName);
	int NewZoneNum = -1;

	if (IJK[0] - DelPointNums.size() > 0 && TecUtilDataSetAddZone(ZoneName, IJK[0] - DelPointNums.size(), IJK[1], IJK[2], ZoneType_Ordered, FDTypes.data())){
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

	int fNum = 0;

	vector<int> CPNums = Fields[fNum].GetReturnIntVec();
	int ZoneNum = Fields[fNum++].GetReturnInt(),
		CPTypeVarNum = Fields[fNum++].GetReturnInt();
	bool DeleteFromOtherCPZones = Fields[fNum++].GetReturnBool();
	vector<int> OtherCPZoneNums;
	vector<vector<int> > OtherCPNums;

	if (DeleteFromOtherCPZones){
		TecUtilDataLoadBegin();

		vector<FieldData_pa> CPZoneXYZRefs, OtherCPZoneXYZRefs[3], OtherCPZoneTypeRefs;
		FieldData_pa CPZoneCPTypeRef;
		int XYZVarNums[3] = { -1, -1, -1 };
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
					NewZoneNums.push_back(NewZoneNum);
					OldZoneNums.push_back(ZoneNum);
				}
			}
			int NewZoneNum = DeletePointsFromIOrderedZone(OtherCPZoneNums[z], OtherCPNums[z]);
			if (NewZoneNum > 0){
				NewZoneNums.push_back(NewZoneNum);
				OldZoneNums.push_back(OtherCPZoneNums[z]);
			}
		}
		if (!ZoneRun){
			int NewZoneNum = DeletePointsFromIOrderedZone(ZoneNum, CPNums);
			if (NewZoneNum > 0){
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

// 	TecUtilStateChanged(StateChange_ZonesAdded, (ArbParam_t)ZoneNumsSet);
// 
// 	TecUtilZoneSetScatter(SV_SHOW, ZoneNumsSet, 0.0, TRUE);
// 	TecUtilZoneSetScatter(SV_FRAMESIZE, ZoneNumsSet, 1, FALSE);
// 	TecUtilZoneSetScatterSymbolShape(SV_GEOMSHAPE, ZoneNumsSet, GeomShape_Sphere);
// 	TecUtilZoneSetMesh(SV_SHOW, ZoneNumsSet, 0.0, FALSE);
// 	TecUtilZoneSetContour(SV_SHOW, ZoneNumsSet, 0.0, FALSE);
// 	TecUtilZoneSetShade(SV_SHOW, ZoneNumsSet, 0.0, FALSE);
// 	TecUtilZoneSetVector(SV_SHOW, ZoneNumsSet, 0.0, FALSE);
// 
// 	TecUtilZoneSetActive(ZoneNumsSet, AssignOp_PlusEquals);

	for (const int & i : NewZoneNums){
		SetCPZone(i);

		if (DeleteFromOtherCPZones && AuxDataZoneItemMatches(i, CSMAuxData.CC.ZoneSubType, CSMAuxData.CC.ZoneTypeCPsAll)){
			Set_pa TmpSet = TecUtilSetAlloc(FALSE);
			TecUtilSetAddMember(TmpSet, i, FALSE);
			TecUtilZoneSetActive(TmpSet, AssignOp_MinusEquals);
			TecUtilSetDealloc(&TmpSet);
		}

	}


	TecUtilSetClear(ZoneNumsSet);
	TecUtilSetAddMember(ZoneNumsSet, ZoneNum, TRUE);
	for (const int & z : OtherCPZoneNums) TecUtilSetAddMember(ZoneNumsSet, z, TRUE);

	TecUtilDataSetDeleteZone(ZoneNumsSet);

	TecUtilSetDealloc(&ZoneNumsSet);
	TecUtilLockFinish(AddOnID);
}

void DeleteCPsGetUserInfo(){
	vector<GuiField_c> Fields = {
		GuiField_c(Gui_ZonePointSelectMulti, "Critical point(s)", CSMZoneName.CriticalPoints),
		GuiField_c(Gui_VarSelect, "CP type", CSMVarName.CritPointType),
		GuiField_c(Gui_Toggle, "Also delete from other CP zone(s)", "1")
	};

	CSMGui("Delete critical point(s)", Fields, DeleteCPsReturnUserInfo, AddOnID);
}

void ExtractCPsReturnUserInfo(const bool GuiSuccess, const vector<GuiField_c> & Fields){
	if (!GuiSuccess) return;

	TecUtilLockStart(AddOnID);

	int fNum = 0;
	vector<int> CPListNums = Fields[fNum].GetReturnIntVec();
	vector<string> CPListNames = Fields[fNum].GetReturnStringVec();
	int ZoneNum = Fields[fNum++].GetReturnInt();

	int NumZones = TecUtilDataSetGetNumZones();

	Set_pa ZoneNumsSet = TecUtilSetAlloc(TRUE);
	vector<int> NewZoneNums;
	vector<vector<int> > CPNums(CPNameList.size(), vector<int>()), CPTypeNums(CPNameList.size(), vector<int>());

	for (int iCP = 0; iCP < CPListNums.size(); ++iCP){
		vector<string> CPName = SplitString(CPListNames[iCP], " ");
		int ind = VectorGetElementNum(CPNameList, CPName[0]);
		CPNums[ind].push_back(CPListNums[iCP]);
		CPTypeNums[ind].push_back(stoi(CPName[1]));
	}

	Set_pa OldZoneSet = TecUtilSetAlloc(FALSE);
	int CPNumTypes = 0;
	for (const auto & t : CPTypeNums) if (t.size() > 0) CPNumTypes++;

	if (CPNumTypes > 1 && AuxDataZoneItemMatches(ZoneNum, CSMAuxData.CC.ZoneSubType, CSMAuxData.CC.ZoneTypeCPsAll)){
		vector<int> AllCPNums;
		for (const auto & i : CPNums) AllCPNums.insert(AllCPNums.end(), i.cbegin(), i.cend());
		int NumCPs, iJunk;
		TecUtilZoneGetIJK(ZoneNum, &NumCPs, &iJunk, &iJunk);
		TecUtilSetAddMember(OldZoneSet, ZoneNum, FALSE);
		/*
		 *	I already made a function for making new zones with certain points deleted, so just use that
		 */
		vector<int> DeleteNums;
		for (int i = 1; i <= NumCPs; ++i){
			if (VectorGetElementNum(AllCPNums, i) < 0) DeleteNums.push_back(i);
		}
		NewZoneNums.push_back(DeletePointsFromIOrderedZone(ZoneNum, DeleteNums));

		
	}


	for (int t = 0; t < CPNameList.size(); ++t){
		if (CPNums[t].size() > 0){
			int SubZoneNum = -1;
			for (int z = ZoneNum; z <= NumZones && SubZoneNum < 0; ++z){
				if (AuxDataZoneItemMatches(z, CSMAuxData.CC.ZoneSubType, CSMAuxData.CC.CPSubTypes[t])) SubZoneNum = z;
			}
			if (SubZoneNum > 0){
				int NumCPs, iJunk;
				TecUtilZoneGetIJK(SubZoneNum, &NumCPs, &iJunk, &iJunk);
				TecUtilSetAddMember(OldZoneSet, SubZoneNum, FALSE);
				vector<int> DeleteNums;
				for (int i = 1; i <= NumCPs; ++i){
					if (VectorGetElementNum(CPTypeNums[t], i) < 0) DeleteNums.push_back(i);
				}
				NewZoneNums.push_back(DeletePointsFromIOrderedZone(SubZoneNum, DeleteNums));
			}
		}
	}

	for (const int & i : NewZoneNums) SetCPZone(i);
	TecUtilZoneSetActive(OldZoneSet, AssignOp_MinusEquals);
	TecUtilSetDealloc(&OldZoneSet);

	TecUtilLockFinish(AddOnID);
}

void ExtractCPsGetUserInfo(){
	vector<GuiField_c> Fields = {
		GuiField_c(Gui_ZonePointSelectMulti, "Critical point(s)", CSMZoneName.CriticalPoints)
	};

	CSMGui("Extract critical point(s)", Fields, ExtractCPsReturnUserInfo, AddOnID);
}

void VarNameFindReplaceReturnUserInfo(const bool GuiSuccess, const vector<GuiField_c> & Fields){
	if (!GuiSuccess) return;

	TecUtilLockStart(AddOnID);

	int FNum = 0;

	string FindStr = Fields[FNum++].GetReturnString(),
		ReplaceStr = Fields[FNum++].GetReturnString();

	vector<int> SearchVars = Fields[FNum++].GetReturnIntVec();

	Set_pa ChangedVars = TecUtilSetAlloc(TRUE);

	for (const int & i : SearchVars){
		char* TmpCStr;

		if (TecUtilVarGetName(i, &TmpCStr)){
			TecUtilVarRename(i, StringReplaceSubString(TmpCStr, FindStr, ReplaceStr).c_str());
			TecUtilSetAddMember(ChangedVars, i, TRUE);
		}

		TecUtilStringDealloc(&TmpCStr);
	}

	TecUtilStateChanged(StateChange_VarsAltered, (ArbParam_t)ChangedVars);
	TecUtilSetDealloc(&ChangedVars);

	TecUtilLockFinish(AddOnID);
}

void VarNameFindReplaceGetUserInfo(){
	vector<GuiField_c> Fields = {
		GuiField_c(Gui_String, "Find string", ""),
		GuiField_c(Gui_String, "Replace string", ""),
		GuiField_c(Gui_VarSelectMulti, "Variables", "")
	};

	CSMGui("Variable name: find replace", Fields, VarNameFindReplaceReturnUserInfo, AddOnID);
}

void ZoneNameFindReplaceReturnUserInfo(const bool GuiSuccess, const vector<GuiField_c> & Fields){
	if (!GuiSuccess) return;

	TecUtilLockStart(AddOnID);

	int FNum = 0;

	string FindStr = Fields[FNum++].GetReturnString(),
		ReplaceStr = Fields[FNum++].GetReturnString();

	vector<int> SearchZones = Fields[FNum++].GetReturnIntVec();

	Set_pa ChangedZones = TecUtilSetAlloc(TRUE);

	for (const int & i : SearchZones){
		char* TmpCStr;

		if (TecUtilZoneGetName(i, &TmpCStr)){
			TecUtilZoneRename(i, StringReplaceSubString(TmpCStr, FindStr, ReplaceStr).c_str());
			TecUtilSetAddMember(ChangedZones, i, TRUE);
		}

		TecUtilStringDealloc(&TmpCStr);
	}

	TecUtilStateChanged(StateChange_ZonesAdded, (ArbParam_t)ChangedZones);
	TecUtilSetDealloc(&ChangedZones);

	TecUtilLockFinish(AddOnID);
}

void ZoneNameFindReplaceGetUserInfo(){
	vector<GuiField_c> Fields = {
		GuiField_c(Gui_String, "Find string", ""),
		GuiField_c(Gui_String, "Replace string", ""),
		GuiField_c(Gui_ZoneSelectMulti, "Zones", "")
	};

	CSMGui("Zone name: find replace", Fields, ZoneNameFindReplaceReturnUserInfo, AddOnID);
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

	StatusLaunch("Reading data into memory...", AddOnID, FALSE);

	VolExtentIndexWeights_s VolInfo;
	VolInfo.AddOnID = AddOnID;

	if (!GetVolInfo(VolZoneNum, XYZVarNums, IsPeriodic, VolInfo)){
		TecUtilDialogErrMsg("Failed to get volume zone info");
		return FALSE;
	}

	if (!GetReadPtrsForZone(VolZoneNum, 
		RhoVarNum, GradVarNums, HessVarNums,
		RhoPtr, GradPtrs, HessPtrs)){
		TecUtilDialogErrMsg("Failed to get read pointer(s)");
		return FALSE;
	}

	CSMGuiLock();

	CritPoints_c VolCPs;

	StatusDrop(AddOnID);

	double RhoCutoff = DefaultRhoCutoff;

	if (FindCPs(VolCPs, VolInfo, CellSpacing, RhoCutoff, IsPeriodic, RhoPtr, GradPtrs, HessPtrs)){
		VolCPs.SaveAsOrderedZone(XYZVarNums, RhoVarNum, TRUE);
	}

	TecUtilDataLoadEnd();

	CSMGuiUnlock();
	TecUtilLockFinish(AddOnID);

	return TRUE;
}


const string MakeStringFromCPNums(const vector<int> & CPNums, const CritPoints_c & CPs, vector<vector<int> > & CPTypeOffsets){
	string Out = " (";

	CPTypeOffsets.clear();

	int NumFarField = 0;
	vector<vector<int> > CPsByType(CPNameList.size());

	for (const auto & CP : CPNums){
		if (CP < 0) NumFarField++;
		else{
			vector<int> CPTypeNumOffset = CPs.GetTypeNumOffsetFromTotOffset(CP);
			CPsByType[CPTypeNumOffset[0]].push_back(CPTypeNumOffset[1]);
		}
	}

	for (auto & v : CPsByType) sort(v.begin(), v.end());
	int CornerNum = 0;
	for (int i = 0; i < CPsByType.size(); ++i){
		for (int j = 0; j < CPsByType[i].size(); ++j){
			CornerNum++;
			Out += CPNameList[i][0] + to_string(CPsByType[i][j] + 1);
			if (CornerNum < CPNums.size()) Out += "-";

			CPTypeOffsets.push_back(vector<int>({ i, CPsByType[i][j] }));
		}
	}

	for (int i = 0; i < NumFarField; ++i){
		CornerNum++;
		Out += "FF";
		if (CornerNum < CPNums.size()) Out += "-";

		CPTypeOffsets.push_back(vector<int>({ -1, -1 }));
	}

	Out += ")";

	return Out;
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

	int TypeInd = std::find(CPTypeList, std::end(CPTypeList), CPType) - CPTypeList;
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
	VolInfo.AddOnID = AddOnID;

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
	mat33 I = eye<mat>(3, 3);
	MR.BasisVectors = &I;


	CritPoints_c CPs(AllCPsZoneNum, XYZVarNums, CPTypeVarNum, RhoVarNum, &MR);

	int EndCPNumforName = 0;
	StreamDir_e GPDir = StreamDir_Forward;
	ColorIndex_t PathColor = Black_C;
	vector<CPType_e> MinDistTypes = { CPType_NuclearCP, CPType_BondCP };
	if (CPType == CPType_RingCP){
		MinDistTypes = { CPType_NuclearCP, CPType_RingCP };
		GPDir = StreamDir_Reverse;
		PathColor = Green_C;
		EndCPNumforName++;
	}



	if (SelectedCPNums.size() > CPs.NumCPs(TypeInd) || SelectedCPNums.back() - 1 > CPs.NumCPs(TypeInd)){
		TecUtilDialogErrMsg("Discrepancy between selected CPs and selected CP zone");
		return FALSE;
	}


	CSMGuiLock();

	vec3 StartPoint;


	const double StartPointOffset = 0.05 * CPs.GetMinCPDist(MinDistTypes);
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
		if (GPDir == StreamDir_Reverse) GPs[iGP].Reverse();
	}

	for (int iGP = 0; iGP < NumGPs; iGP += 2){
		GradPath_c GP = GPs[iGP];
		vector<int> CPNums;
		CPNums.push_back(GP.GetStartEndCPNum()[EndCPNumforName]);
		int OldStartCPNum = GP.GetStartEndCPNum()[EndCPNumforName];
		GP += GPs[iGP + 1];
		vector<int> StartEndCPNums = GP.GetStartEndCPNum();
		vector<vector<int> > StartEndCPTypeAndOffset(2);
		StartEndCPTypeAndOffset[0] = CPs.GetTypeNumOffsetFromTotOffset(OldStartCPNum);
		string Name = (CPType == CPType_BondCP ? CSMZoneName.BondPath : CSMZoneName.RingLine) + CPNameList[StartEndCPTypeAndOffset[0][0]] + " " + to_string(StartEndCPTypeAndOffset[0][1] + 1);
		Name += MakeStringFromCPNums(StartEndCPNums, CPs, StartEndCPTypeAndOffset);
// 		for (int i = 0; i < 2; ++i){
// 			if (StartEndCPNums[i] >= 0){
// 				StartEndCPTypeAndOffset[i] = CPs.GetTypeNumOffsetFromTotOffset(StartEndCPNums[i]);
// 				Name += CPNameList[StartEndCPTypeAndOffset[i][0]][0] + to_string(StartEndCPTypeAndOffset[i][1] + 1);
// 			}
// 			else Name += "FF";
// 
// 			if (i == 0) Name += "-";
// 		}
// 
// 		Name += ")";

		if (StartEndCPNums[0] < 0){
			int TmpInt = StartEndCPNums[0];
			StartEndCPNums[0] = StartEndCPNums[1];
			StartEndCPNums[1] = TmpInt;
		}
		GP.SaveAsOrderedZone(Name, vector<FieldDataType_e>(), XYZVarNums, RhoVarNum, TRUE, PathColor);
		for (int i = 0; i < 2; ++i){
			AuxDataZoneSetItem(GP.GetZoneNum(), CSMAuxData.CC.GPEndNumStrs[i], to_string(StartEndCPNums[i] + 1));
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
		string Name = CSMZoneName.SpecialGradientPath;
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
			AuxDataZoneSetItem(GP.GetZoneNum(), CSMAuxData.CC.GPEndNumStrs[i], to_string(StartEndCPNums[i] + 1));
			if (StartEndCPNums[i] >= 0) AuxDataZoneSetItem(GP.GetZoneNum(), CSMAuxData.CC.GPEndTypes[i], CPNameList[StartEndCPTypeAndOffset[i][0]]);
			else AuxDataZoneSetItem(GP.GetZoneNum(), CSMAuxData.CC.GPEndTypes[i], "FF");
		}
		AuxDataZoneSetItem(GP.GetZoneNum(), CSMAuxData.CC.ZoneType, CSMAuxData.CC.ZoneTypeGP);
		if (CPType == CPType_RingCP) AuxDataZoneSetItem(GP.GetZoneNum(), CSMAuxData.CC.ZoneSubType, CSMAuxData.CC.ZoneSubTypeRingLine);
		else AuxDataZoneSetItem(GP.GetZoneNum(), CSMAuxData.CC.ZoneSubType, CSMAuxData.CC.ZoneSubTypeBondPath);
	}

	TecUtilDataLoadEnd();


	CSMGuiUnlock();
	TecUtilLockFinish(AddOnID);

	return TRUE;
}

double MinFunc_GPLength_InPlane(double alpha, void * params){
	MinFuncParams_GPLengthInPlane * GPParams = reinterpret_cast<MinFuncParams_GPLengthInPlane*>(params);

	vec3 StartPoint = GPParams->StartPointOrigin + Rotate(GPParams->RotVec, alpha, GPParams->RotAxis);

	GPParams->GP = GradPath_c(
		StartPoint,
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

template <typename T>
const T PeriodicDistance1D(T a, T b, const T & low, const T & high){
	T span = high - low;
	T halfSpan = span / static_cast<T>(2);

	T rawDist = a - b;
	rawDist = rawDist > 0 ? rawDist : -rawDist;

	if (rawDist > halfSpan){
		if (b > a){
			b -= span;
		}
		else a -= span;
	}

	T dist = a - b;
	dist = dist > 0 ? dist : -dist;

	return dist;
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
	const int & RCSFuncVarNum,
	const Boolean_t & IsPeriodic)
{
	TecUtilLockStart(AddOnID);

	int NumZones = TecUtilDataSetGetNumZones();
	int NumVars = TecUtilDataSetGetNumVars();

	REQUIRE(AllCPsZoneNum > 0 && AllCPsZoneNum <= NumZones);
	REQUIRE(XYZVarNums.size() == 3);
	for (const auto & i : XYZVarNums) REQUIRE(i > 0 && i <= NumVars);

	int TypeInd = std::find(CPTypeList, std::end(CPTypeList), CPType) - CPTypeList;
	if (TypeInd >= 6){
		TecUtilDialogErrMsg("Invalid CP type specified");
		return FALSE;
	}

	if (CPTypeVarNum < 0){
		TecUtilDialogErrMsg("Failed to find CP Type variable");
		return FALSE;
	}

	FieldDataPointer_c RhoPtr, RCSFuncPtr;
	vector<FieldDataPointer_c> GradPtrs, HessPtrs;

	TecUtilDataLoadBegin();

	if (RCSFuncVarNum > 0 && !RCSFuncPtr.GetReadPtr(VolZoneNum, RCSFuncVarNum)){
		TecUtilDialogErrMsg("Failed to get read pointers");
		return FALSE;
	}

	if (!GetReadPtrsForZone(VolZoneNum,
		RhoVarNum, GradVarNums, HessVarNums,
		RhoPtr, GradPtrs, HessPtrs)){
		TecUtilDialogErrMsg("Failed to get read pointer(s)");
		return FALSE;
	}

	VolExtentIndexWeights_s VolInfo;
	VolInfo.AddOnID = AddOnID;

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
	mat33 I = eye<mat>(3, 3);
	MR.BasisVectors = &I;


	CritPoints_c CPs(AllCPsZoneNum, XYZVarNums, CPTypeVarNum, RhoVarNum, &MR);

	if (SelectedCPNums.size() > CPs.NumCPs(TypeInd) || SelectedCPNums.back() - 1 > CPs.NumCPs(TypeInd)){
		TecUtilDialogErrMsg("Discrepancy between selected CPs and selected CP zone");
		return FALSE;
	}

	CSMGuiLock();

	vector<CPType_e> CPTypesForMinDistCheck;
	int EndCPNumforName = 0;
	StreamDir_e GPDir = StreamDir_Reverse;
	ColorIndex_t PathColor = Black_C;
	CPType_e PrimaryTermType = CPType_CageCP,
		SecondaryTermType = CPType_RingCP;
	if (CPType == CPType_RingCP){
		GPDir = StreamDir_Forward;
		PathColor = Green_C;
		EndCPNumforName++;
		PrimaryTermType = CPType_NuclearCP;
		SecondaryTermType = CPType_BondCP;
		CPTypesForMinDistCheck = { CPType_BondCP, CPType_RingCP };
	}
	else{
		if (CPs.NumBonds() <= 1) CPTypesForMinDistCheck.push_back(CPType_NuclearCP);
		
		CPTypesForMinDistCheck.push_back(CPType_BondCP);

		if (CPs.NumRings() > 0) CPTypesForMinDistCheck.push_back(CPType_RingCP);
	}



	vec3 StartPoint;

	const int NumGPPts = 100;
	double TermRadius = 0.2;
	double RhoCutoff = DefaultRhoCutoff;
	int NumCircleCheckGPs = MinNumCircleCheckGPs;
	int NumCircleGPs = MinNumCircleGPs;

	NumCPs = SelectedCPNums.size();
// 	NumCPs = CPs.NumCPs(TypeInd);


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


	double StartPointOffset;
	double AngleStep;

	string StatusStr = "Finding " + string(CPType == CPType_BondCP ? "interatomic" : "ring") + " surfaces";
	StatusLaunch(StatusStr, AddOnID, TRUE);

#ifndef _DEBUG
// #pragma omp parallel for schedule(dynamic)
#endif
	for (int iCP = 0; iCP < NumCPs; ++iCP){
		if (!StatusUpdate(iCP + 1, NumCPs, StatusStr + ": (" + to_string(iCP+1) + " of " + to_string(NumCPs) + ")", AddOnID)) break;
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

		StartPointOffset = 0.05 * CPs.GetMinCPDist(CPTypesForMinDistCheck);
		vec3 StartVec = normalise(CPs.GetEigVecs(TypeInd, cpNum).col(1)) * StartPointOffset;

		GPParams.StartPointOrigin = CPs.GetXYZ(TypeInd, cpNum);
		GPParams.RotVec = StartVec;
		GPParams.RotAxis = CPs.GetPrincDir(TypeInd, cpNum);

		GPParams.StartCPNum = CPs.GetTotOffsetFromTypeNumOffset(TypeInd, cpNum);

		AngleStep = PI2 / static_cast<double>(NumCircleCheckGPs);

		vector<GradPath_c> GPs;
		vector<double> GPLengths(NumCircleCheckGPs);
		vector<CPType_e> GPEndCPTypes(NumCircleCheckGPs);
		vector<int> GPEndCPNums(NumCircleCheckGPs);

		/*
		 *	Seed GPs around the CP
		 */
		for (int gpNum = 0; gpNum < NumCircleCheckGPs; ++gpNum){
			double RotAngle = static_cast<double>(gpNum) * AngleStep;

			StartPoint = GPParams.StartPointOrigin + Rotate(StartVec, RotAngle, GPParams.RotAxis);

			GPs.push_back(GradPath_c(StartPoint, GPDir, NumGPPts, GPType_Classic, GPTerminate_AtCP, NULL, &CPs, &TermRadius, &RhoCutoff, VolInfo, HessPtrs, GradPtrs, RhoPtr));
			GPs.back().SetStartEndCPNum(GPParams.StartCPNum, 0);
		}

#ifndef _DEBUG
#pragma omp parallel for schedule(dynamic)
#endif
		for (int gpNum = 0; gpNum < NumCircleCheckGPs; ++gpNum){
			GPs[gpNum].Seed(false);
			GPs[gpNum].PointPrepend(GPParams.StartPointOrigin, CPs.GetRho(TypeInd, cpNum));
			GPLengths[gpNum] = GPs[gpNum].GetLength();
			GPEndCPNums[gpNum] = GPs[gpNum].GetStartEndCPNum(1);
			GPEndCPTypes[gpNum] = CPs.GetTypeFromTotOffset(GPEndCPNums[gpNum]);
		}

		// DEBUG
// 		for (int gpNum = 0; gpNum < NumCircleCheckGPs; ++gpNum){
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
// 
// 		StatusDrop(AddOnID);
// 		CSMGuiUnlock();
//  		return TRUE;

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
		vector<double> SGPSeedAngles;
		vector<int> TermCPNums;
		
		vector<bool> GPUsed(NumCircleCheckGPs, false);
		for (int gpNum = 0; gpNum < NumCircleCheckGPs; ++gpNum){
			if (!GPUsed[gpNum] || (GPEndCPTypes[gpNum] == PrimaryTermType
				&& GPEndCPTypes[(gpNum + 1) % NumCircleCheckGPs] == PrimaryTermType
				&& GPEndCPNums[gpNum] != GPEndCPNums[(gpNum + 1) % NumCircleCheckGPs])
				|| (GPEndCPTypes[gpNum] == CPType_Invalid
				&& GPEndCPTypes[(gpNum + 1) % NumCircleCheckGPs] == CPType_Invalid
				&& dot(GPs[gpNum][-1] - GPs[gpNum][-2], GPs[(gpNum + 1) % NumCircleCheckGPs][-1] - GPs[(gpNum + 1) % NumCircleCheckGPs][-2]) < Tolerance_DivergentGPInnerProd)){

				for (int CheckNum = 0; CheckNum < 2; ++CheckNum){
					int TermCPNum = -1;
					int MinGPNum;
					double GPLen = -1.;
					double AlphaLower, AlphaUpper, AlphaGuess;
					double AngleFactor;
					if (CheckNum == 0 && (GPEndCPTypes[gpNum] == PrimaryTermType || GPEndCPTypes[gpNum] == SecondaryTermType)
						&& find(TermCPNums.begin(), TermCPNums.end(), GPEndCPNums[gpNum]) == TermCPNums.end()){
						/*
						 *	Only process if this TermCPNum isn't in the list already
						 */
						/*
						*	Found a GP terminating at a CP, so now get minimum length GP
						*	terminating at the same CP
						*/
						TermCPNum = GPs[gpNum].GetStartEndCPNum(1);
						MinGPNum = gpNum;
						for (int i = gpNum; i < gpNum + NumCircleCheckGPs; ++i){
							int ii = i % NumCircleCheckGPs;
							if (!GPUsed[ii] && GPs[ii].GetStartEndCPNum(1) == TermCPNum){
								GPUsed[ii] = true;
								if (GPLengths[ii] < GPLengths[MinGPNum]){
									MinGPNum = ii;
								}
							}
						}

						AlphaGuess = static_cast<double>(MinGPNum)* AngleStep;
						GPLen = GPLengths[MinGPNum];
					}
					else if (CheckNum == 1 && (GPEndCPTypes[gpNum] == PrimaryTermType
						&& GPEndCPTypes[(gpNum + 1) % NumCircleCheckGPs] == PrimaryTermType
						&& GPEndCPNums[gpNum] != GPEndCPNums[(gpNum + 1) % NumCircleCheckGPs])
						|| (GPEndCPTypes[gpNum] == CPType_Invalid
						&& GPEndCPTypes[(gpNum + 1) % NumCircleCheckGPs] == CPType_Invalid
						&& dot(GPs[gpNum][-1] - GPs[gpNum][-2], GPs[(gpNum + 1) % NumCircleCheckGPs][-1] - GPs[(gpNum + 1) % NumCircleCheckGPs][-2]) < Tolerance_DivergentGPInnerProd)){
						/*
						*	The gpNum'th gradient path and the (gpNum + 1)'th gradient path diverge to two CPs.
						*	This means there is a (ring) saddle CP to be found between these gradient paths.
						*	We'll use a binary search to get at the CP sandwiched between these two.
						*/
						int CPNumLower, CPNumUpper;
						vec3 TermVecLower, TermVecUpper;

						if (GPEndCPTypes[gpNum] == PrimaryTermType){
							CPNumLower = GPEndCPNums[gpNum];
							CPNumUpper = GPEndCPNums[(gpNum + 1) % NumCircleCheckGPs];
						}
						else{
							TermVecLower = GPs[gpNum][-1] - GPs[gpNum][-2];
							TermVecUpper = GPs[(gpNum + 1) % NumCircleCheckGPs][-1] - GPs[(gpNum + 1) % NumCircleCheckGPs][-2];
						}

						AlphaLower = static_cast<double>(gpNum)* AngleStep;
						AlphaUpper = static_cast<double>(gpNum + 1) * AngleStep;

						while (TermCPNum < 0){
							AlphaGuess = (AlphaLower + AlphaUpper) * 0.5;
							StartPoint = GPParams.StartPointOrigin + Rotate(StartVec, AlphaGuess, GPParams.RotAxis);
							GradPath_c GP(StartPoint, GPDir, NumGPPts, GPType_Classic, GPTerminate_AtCP, NULL, &CPs, &TermRadius, &RhoCutoff, VolInfo, HessPtrs, GradPtrs, RhoPtr);
							GP.SetStartEndCPNum(GPParams.StartCPNum, 0);
							GP.Seed(false);
							if (CPs.GetTypeFromTotOffset(GP.GetStartEndCPNum(1)) == SecondaryTermType){
								TermCPNum = GP.GetStartEndCPNum(1);
								GP.PointPrepend(GPParams.StartPointOrigin, CPs.GetRho(TypeInd, cpNum));
								GPLen = GP.GetLength();
							}
							else{
								bool GoLeft;
								if (GPEndCPTypes[gpNum] == PrimaryTermType){
									if (GP.GetStartEndCPNum(1) == CPNumLower) GoLeft = false;
									else if (GP.GetStartEndCPNum(1) == CPNumUpper) GoLeft = true;
									else{
										TecUtilDialogErrMsg(string("GP terminated at unexpected CP when performing binary search.\n\ngpNum = " + to_string(gpNum)).c_str());
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
								StartPoint = GPParams.StartPointOrigin + Rotate(StartVec, RotAngle, GPParams.RotAxis);
								GradPath_c GP(StartPoint, GPDir, NumGPPts, GPType_Classic, GPTerminate_AtCP, NULL, &CPs, &TermRadius, &RhoCutoff, VolInfo, HessPtrs, GradPtrs, RhoPtr);
								GP.SetStartEndCPNum(GPParams.StartCPNum, 0);
								GP.Seed(false);
								if (GP.GetStartEndCPNum(1) == TermCPNum){
									GP.PointPrepend(GPParams.StartPointOrigin, CPs.GetRho(TypeInd, cpNum));
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
								StartPoint = GPParams.StartPointOrigin + Rotate(StartVec, RotAngle, GPParams.RotAxis);
								GradPath_c GP(StartPoint, GPDir, NumGPPts, GPType_Classic, GPTerminate_AtCP, NULL, &CPs, &TermRadius, &RhoCutoff, VolInfo, HessPtrs, GradPtrs, RhoPtr);
								GP.SetStartEndCPNum(GPParams.StartCPNum, 0);
								GP.Seed(false);
								if (GP.GetStartEndCPNum(1) == TermCPNum){
									GP.PointPrepend(GPParams.StartPointOrigin, CPs.GetRho(TypeInd, cpNum));
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
		}

		if (MinFunc_AlphaLowMidHigh.size() != TermCPNums.size()){
			TecUtilDialogErrMsg("Number of starting angle values doesn't match number of terminating CPs");
			continue;
		}

		int NumSGPs = MinFunc_AlphaLowMidHigh.size();
		SGPsPerCP[iCP].resize(NumSGPs);
		vector<vector<GradPath_c> > SaddleCPSupplementGPs(NumSGPs, vector<GradPath_c>(2));
		SGPSeedAngles.resize(NumSGPs);
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
			vector<int> TmpEndCPNums(3);
			vector<vec3> TmpStartPoints(3);
			for (int j = 0; j < 3; ++j){
				TmpStartPoints[j] = GPParams.StartPointOrigin + Rotate(StartVec, MinFunc_AlphaLowMidHigh[i][j], GPParams.RotAxis);
				GradPath_c GP(TmpStartPoints[j], GPDir, NumGPPts, GPType_Classic, GPTerminate_AtCP, NULL, &CPs, &TermRadius, &RhoCutoff, VolInfo, HessPtrs, GradPtrs, RhoPtr);

				GP.SetStartEndCPNum(GPParams.StartCPNum, 0);
				GP.Seed(false);
				GP.PointPrepend(GPParams.StartPointOrigin, CPs.GetRho(TypeInd, cpNum));

				TmpLens[j] = GP.GetLength();
				TmpEndCPNums[j] = GP.GetStartEndCPNum(1);
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

// 			SGPsPerCP[iCP][i] = reinterpret_cast<MinFuncParams_GPLengthInPlane*>(F.params)->GP;
			SGPsPerCP[iCP][i].SetupGradPath(GPParams.StartPointOrigin + Rotate(StartVec, MinFunc_AlphaLowMidHigh[i][1], GPParams.RotAxis), GPDir, NumGPPts, GPType_Classic, GPTerminate_AtCP, NULL, &CPs, &TermRadius, &RhoCutoff, VolInfo, HessPtrs, GradPtrs, RhoPtr);
			SGPsPerCP[iCP][i].SetStartEndCPNum(GPParams.StartCPNum, 0);
			SGPsPerCP[iCP][i].Seed(false);
			SGPsPerCP[iCP][i].PointPrepend(GPParams.StartPointOrigin, CPs.GetRho(TypeInd, cpNum));
			SGPsPerCP[iCP][i].Resample(NumGPPts);

			SGPSeedAngles[i] = MinFunc_AlphaLowMidHigh[i][1];

			if (GPDir == StreamDir_Reverse) SGPsPerCP[iCP][i].Reverse();

			if (CPs.GetTypeFromTotOffset(TermCPNums[i]) == SecondaryTermType){
				/*
				*	Need to find the direction normal to the SGP at the bond/ring CP and in the plane of the ring/bond.
				*	Use unit vectors for 1) the last two points of the SGP and 2) a vector going from the originating
				*	bond/ring CP to another arbitrary point in the plane, to get three orthonormal unit vectors to use.
				*/
				vec3 pt1, pt2, pt3, v1, v2, v3;
				vec3 CPPos = CPs.GetXYZ(TermCPNums[i]);
				double RotAngle = SGPSeedAngles[i] + PI * 0.25;
				pt3 = GPParams.StartPointOrigin + Rotate(StartVec, RotAngle, GPParams.RotAxis);

				if (GPDir == StreamDir_Forward){
					pt1 = SGPsPerCP[iCP][i][-1];
					pt2 = SGPsPerCP[iCP][i][-2];
				}
				else{
					pt1 = SGPsPerCP[iCP][i][0];
					pt2 = SGPsPerCP[iCP][i][1];
				}

// 				v1 = normalise(pt2 - pt1);
				v1 = normalise(GPParams.StartPointOrigin - CPPos);
				v2 = pt3 - pt2;
				v3 = normalise(cross(v2, v1));
				v2 = normalise(cross(v3, v1));

				/*
				*	One of the GPs we seed will be the first of a surface segment and pointing in the positive rotation direction,
				*	and the other the last of a segment pointing in the negative rotation direction,
				*	so need to store them for later.
				*/
				vec3 PtUp = CPPos + v2 * 0.05 * CPs.GetMinCPDist(CPTypesForMinDistCheck);
				vec3 PtDn = CPPos - v2 * 0.05 * CPs.GetMinCPDist(CPTypesForMinDistCheck);

				if (dot(v2, pt3 - GPParams.StartPointOrigin) > 0){
					/*
					 *	v2 points in the positive rotation direction
					 */

					SaddleCPSupplementGPs[i][0].SetupGradPath(PtUp, GPDir, NumGPPts,
						GPType_Classic, GPTerminate_AtRhoValue, NULL,
						&CPs, &TermRadius, &RhoCutoff, VolInfo,
						HessPtrs, GradPtrs, RhoPtr);
					SaddleCPSupplementGPs[i][1].SetupGradPath(PtDn, GPDir, NumGPPts,
						GPType_Classic, GPTerminate_AtRhoValue, NULL,
						&CPs, &TermRadius, &RhoCutoff, VolInfo,
						HessPtrs, GradPtrs, RhoPtr);
					for (GradPath_c & GP : SaddleCPSupplementGPs[i]){
						GP.SetStartEndCPNum(TermCPNums[i], 0);
						GP.Seed(true);
						//DEBUG
						GP.PointPrepend(CPPos, 0.0);
// 						GP.SaveAsOrderedZone("Supp GP", Red_C);
					}
				}
				else{
					/*
					*	v2 points in the positive rotation direction
					*/

					SaddleCPSupplementGPs[i][0].SetupGradPath(PtDn, GPDir, NumGPPts,
					GPType_Classic, GPTerminate_AtRhoValue, NULL,
					&CPs, &TermRadius, &RhoCutoff, VolInfo,
					HessPtrs, GradPtrs, RhoPtr);
					SaddleCPSupplementGPs[i][1].SetupGradPath(PtUp, GPDir, NumGPPts,
						GPType_Classic, GPTerminate_AtRhoValue, NULL,
						&CPs, &TermRadius, &RhoCutoff, VolInfo,
						HessPtrs, GradPtrs, RhoPtr);
					for (GradPath_c & GP : SaddleCPSupplementGPs[i]){
						GP.SetStartEndCPNum(TermCPNums[i], 0);
						GP.Seed(true);
						//DEBUG
						GP.PointPrepend(CPPos, 0.0);
// 						GP.SaveAsOrderedZone("Supp GP", Red_C);
					}
				}
			}
		}

		/*
		 * All near-field SGPs found. Now look for far-field SGPs using the RCS function
		 */

		if (RCSFuncPtr.IsReady()){
			/*
			 *	Get the RCS function values at points around the CP 
			 *	in the normal plane at constant radius
			 */
			int NumRCSFuncCircleCheckPts = MinNumRCSFuncCircleCheckPts;
			AngleStep = PI2 / static_cast<double>(NumRCSFuncCircleCheckPts);

			int RCSNumCheckRadii = DefaultRCSNumCheckRadii;
			double RCSCheckRadiiLow = DefaultRCSCheckRadiiLow;
			double RCSCheckRadiiHigh = DefaultRCSCheckRadiiHigh;
			double RCSCheckRadiiHighMultiplier = DefaultRCSCheckRadiiHighMultiplier;
			bool LongRangeRCSCheck = (SGPsPerCP[iCP].size() <= 1);
			if (LongRangeRCSCheck){
				RCSCheckRadiiHigh *= RCSCheckRadiiHighMultiplier;
				StartPointOffset = RCSCheckRadiiHigh * CPs.GetMinCPDist(CPTypesForMinDistCheck);
			}
			else StartPointOffset = RCSCheckRadiiLow * CPs.GetMinCPDist(CPTypesForMinDistCheck);

			double RCSCheckRadiiStep = (RCSCheckRadiiHigh - RCSCheckRadiiLow) / static_cast<double>(DefaultRCSNumCheckRadii);

			int iter = 0;
			int NumNeighborPtsToCkeck = DefaultRCSFuncCheckNumPts;
			int RCSNumCheckRadiiConverged = DefaultRCSNumCheckRadiiConverged;
			int RCSFuncConvergence = DefaultRCSFuncConvergence;
			int RCSSGPAngleCheck = DefaultRCSSGPAngleCheck;
			double RCSSGPRadCheck = PI2 / static_cast<double>(RCSSGPAngleCheck);
			vector<vector<int> > MaxPtNums, MinPtNums;

			bool IsConverged = false;

			while (iter < RCSNumCheckRadii && !IsConverged){
				iter++;

				MaxPtNums.push_back(vector<int>());
				MinPtNums.push_back(vector<int>());

				vec3 StartVec = normalise(GPParams.RotVec) * StartPointOffset;

				vector<double> RCSFuncCheckVals(NumRCSFuncCircleCheckPts);

				for (int i = 0; i < NumRCSFuncCircleCheckPts; ++i){
					double RotAngle = AngleStep * static_cast<double>(i);
					StartPoint = GPParams.StartPointOrigin + Rotate(StartVec, RotAngle, GPParams.RotAxis);
					RCSFuncCheckVals[i] = ValAtPointByPtr(StartPoint, VolInfo, RCSFuncPtr);
				}

				/*
				 *	Find the points that are local max/min
				 */
				for (int i = 0; i < NumRCSFuncCircleCheckPts; ++i){
					bool IsMax = true;
					bool IsMin = true;

					for (int n = 0; n < NumNeighborPtsToCkeck && (IsMax || IsMin); ++n){
						vector<int> im1p1 = {
							i - n - 1,
							i - n,
							i + n,
							i + n + 1
						};
						for (int & j : im1p1) j = (j + NumRCSFuncCircleCheckPts) % NumRCSFuncCircleCheckPts;

						IsMax = (IsMax && RCSFuncCheckVals[im1p1[0]] <= RCSFuncCheckVals[im1p1[1]] && RCSFuncCheckVals[im1p1[2]] >= RCSFuncCheckVals[im1p1[3]]);
						IsMin = (IsMin && RCSFuncCheckVals[im1p1[0]] >= RCSFuncCheckVals[im1p1[1]] && RCSFuncCheckVals[im1p1[2]] <= RCSFuncCheckVals[im1p1[3]]);
					}

					if (IsMax) MaxPtNums.back().push_back(i);
					else if (IsMin) MinPtNums.back().push_back(i);
				}

				if (MaxPtNums.back().size() > 0 && MaxPtNums.back().size() == MinPtNums.back().size()){
					/*
					 * Check that the sets of min/max points have converged.
					 * This check that the most recent NumConvergedRCSFuncCheckRadii of
					 * min/max point numbers have converged within the tolerance specified by
					 * RCSFuncConvergence.
					 * Convergence means that
					 * 1.	the same number of max/min points were found for
					 *		NumConvergedRCSFuncCheckRadii different radii, and
					 * 2.	the difference in point number between the same max/min
					 *		for two subsequent radii is less or equal to
					 *		RCSFuncConvergence
					 */
					if (MaxPtNums.size() >= RCSNumCheckRadiiConverged && MinPtNums.size() >= RCSNumCheckRadiiConverged){
						IsConverged = true;
						for (int l = 1; l < RCSNumCheckRadiiConverged && IsConverged; ++l){
							IsConverged = (
								MaxPtNums[MaxPtNums.size() - l].size() == MaxPtNums[MaxPtNums.size() - l - 1].size()
								&& MinPtNums[MinPtNums.size() - l].size() == MinPtNums[MinPtNums.size() - l - 1].size()
								&& MaxPtNums[MaxPtNums.size() - l].size() == MinPtNums[MinPtNums.size() - l].size()
								);
						}
						for (int l = 1; l < RCSNumCheckRadiiConverged && IsConverged; ++l){
							int DiffSum = 0;
							for (int p = 0; p < MaxPtNums.back().size() && IsConverged; ++p){
								IsConverged = abs((MaxPtNums[MaxPtNums.size() - l][p] - MaxPtNums[MaxPtNums.size() - l - 1][p])) <= RCSFuncConvergence
									&& abs((MinPtNums[MinPtNums.size() - l][p] - MinPtNums[MinPtNums.size() - l - 1][p])) <= RCSFuncConvergence;
							}
						}
					}
					if (IsConverged){
						/*
						 *	Now that the min/max points have converged, take the average point num for each max/min
						 *	and use that as the starting guess for the minimization.
						 */
						vector<double> AvgMaxNum(MaxPtNums.back().size()), AvgMinNum(MinPtNums.back().size());
						for (int p = 0; p < MaxPtNums.back().size() && IsConverged; ++p){
							for (int l = 1; l <= RCSNumCheckRadiiConverged && IsConverged; ++l){
								AvgMaxNum[p] += MaxPtNums[MaxPtNums.size() - l][p];
								AvgMinNum[p] += MinPtNums[MinPtNums.size() - l][p];
							}
							AvgMaxNum[p] /= static_cast<double>(RCSNumCheckRadiiConverged);
							AvgMinNum[p] /= static_cast<double>(RCSNumCheckRadiiConverged);
						}

						/*
						 *	Now check the angles found from the RCS function against the angles
						 *	for the near field SGPs found above.
						 *	If the angles are too close, compare the lengths of the two
						 *	GPs and keep the shorter one.
						 */

						for (const double & iRCS : AvgMinNum){
							double aRCS = iRCS * AngleStep;
							StartPoint = GPParams.StartPointOrigin + Rotate(GPParams.RotVec, aRCS, GPParams.RotAxis);
							GradPath_c GP(StartPoint, GPDir, NumGPPts, GPType_Classic, GPTerminate_AtCP, NULL, &CPs, &TermRadius, &RhoCutoff, VolInfo, HessPtrs, GradPtrs, RhoPtr);
							GP.SetStartEndCPNum(CPs.GetTotOffsetFromTypeNumOffset(TypeInd, cpNum), 0);
							GP.Seed(false);
							GP.PointPrepend(GPParams.StartPointOrigin, CPs.GetRho(TypeInd, cpNum));
							GP.Resample(NumGPPts);
							if (GPDir == StreamDir_Reverse) GP.Reverse();
							bool GPFound = false;
							for (int iSGP = 0; iSGP < SGPsPerCP[iCP].size() && !GPFound; ++iSGP){
								double aSGP = SGPSeedAngles[iSGP];
								if (PeriodicDistance1D(aSGP, aRCS, 0., PI2) < RCSSGPRadCheck){
									GPFound = true;
									if (GP.GetLength() < SGPsPerCP[iCP][iSGP].GetLength()){
										SGPsPerCP[iCP][iSGP] = GP;
										SGPSeedAngles[iSGP] = aRCS;
									}
								}
							}
							if (!GPFound){
								SGPsPerCP[iCP].push_back(GP);
								SGPSeedAngles.push_back(aRCS);
							}
						}
					}

				}

				if (LongRangeRCSCheck) StartPointOffset -= RCSCheckRadiiStep;
				else StartPointOffset += RCSCheckRadiiStep;
			}
		}

		vector<vector<GradPath_c> > SurfaceGPs;
		vector<vector<int> > SurfaceGPCPNums;
		NumSGPs = SGPsPerCP[iCP].size();

		if (NumSGPs > 1){
			/*
			 *	All SGPs found. Now fill in remaining space with regular GPs in order to
			 *	construct surfaces.
			 *	First get all seed angles shifted to be between 0 and 2pi, then sort them.
			 */
			for (double & a : SGPSeedAngles){
				if (a < 0) a += PI2;
				else if (a >= PI2) a -= PI2;
			}
			vector<double> OldSeedAngles = SGPSeedAngles;
			std::sort(SGPSeedAngles.begin(), SGPSeedAngles.end());

			/*
			 *	Also want to know the mapping from the old order (of SGPs) to the sorted order.
			 *	So SGPSeedAngles[0] == OldSeedAngles[SGPNumsNewToOld[0]]
			 */
			vector<int> SGPNumsNewToOld(NumSGPs);
			for (int n = 0; n < NumSGPs; ++n){
				for (int o = 0; o < NumSGPs; ++o){
					if (SGPSeedAngles[n] == OldSeedAngles[o]){
						SGPNumsNewToOld[n] = o;
						break;
					}
				}
			}

			/*
			 *	Get all the GPs (and SGPs) for the surface.
			 *	Now find out how many GPs need to go in between each pair of SGPs in order to satisfy
			 *	NumCircleGPs.
			 */
			double GPSepAngle = PI2 / static_cast<double>(NumCircleGPs);
			double TwoPi = PI2;

			for (int i = 0; i < NumSGPs; ++i){
				int ip1 = (i + 1) % NumSGPs;
				double AngleDistance = SGPSeedAngles[ip1] - SGPSeedAngles[i];
				if (i == NumSGPs - 1) AngleDistance = PI2 + AngleDistance;
				double AngleRatio = AngleDistance / TwoPi;
				int NumSepGPs = static_cast<int>(static_cast<double>(NumCircleGPs) * AngleRatio);
				double SepAngleStep = AngleDistance / static_cast<double>(NumSepGPs);

				SurfaceGPs.push_back(vector<GradPath_c>());
				SurfaceGPCPNums.push_back(vector<int>({ GPParams.StartCPNum }));
				/*
				 *	For SGPs terminating at bond/ring CPs, need to seed two
				 *	GPs "above" and "below" in order to have a complete surface.
				 */
				int TermCPNum = -1;
				if (SGPNumsNewToOld[i] < TermCPNums.size()) TermCPNum = TermCPNums[SGPNumsNewToOld[i]];
				if (TermCPNum < 0 || CPs.GetTypeFromTotOffset(TermCPNum) == PrimaryTermType){
					SurfaceGPs.back().push_back(SGPsPerCP[iCP][SGPNumsNewToOld[i]]);
					SurfaceGPCPNums.back().push_back(SGPsPerCP[iCP][SGPNumsNewToOld[i]].GetStartEndCPNum(EndCPNumforName));
				}
				else{
					GradPath_c GP = SGPsPerCP[iCP][SGPNumsNewToOld[i]];
// #ifdef _DEBUG
// 					GP.SaveAsOrderedZone("Test SGP", Yellow_C);
// 					SaddleCPSupplementGPs[SGPNumsNewToOld[i]][0].SaveAsOrderedZone("Test Supp GP", Black_C);
// #endif
					GP.ConcatenateResample(SaddleCPSupplementGPs[SGPNumsNewToOld[i]][0], NumGPPts);
					SurfaceGPCPNums.back().push_back(SaddleCPSupplementGPs[SGPNumsNewToOld[i]][0].GetStartEndCPNum(0));

// #ifdef _DEBUG
// 					GP.SaveAsOrderedZone("Concatenated GP front", Blue_C);
// #endif

					SurfaceGPs.back().push_back(GP);
				}

				/*
				 *	Now add the filler GPs
				 */
				for (int j = 1; j < NumSepGPs; ++j){
					double RotAngle = SGPSeedAngles[i] + SepAngleStep * static_cast<double>(j);
					StartPoint = GPParams.StartPointOrigin + Rotate(StartVec, RotAngle, GPParams.RotAxis);
					GradPath_c GP(StartPoint, GPDir, NumGPPts, GPType_Classic, GPTerminate_AtRhoValue, NULL, &CPs, &TermRadius, &RhoCutoff, VolInfo, HessPtrs, GradPtrs, RhoPtr);
					GP.SetStartEndCPNum(CPs.GetTotOffsetFromTypeNumOffset(TypeInd, cpNum), 0);
					GP.Seed(false);
					GP.PointPrepend(GPParams.StartPointOrigin, CPs.GetRho(TypeInd, cpNum));
					GP.Resample(NumGPPts);
					if (GPDir == StreamDir_Reverse) GP.Reverse();
// #ifdef _DEBUG
// 					GP.SaveAsOrderedZone("Filler GP");
// #endif
					SurfaceGPs.back().push_back(GP);
				}
				TermCPNum = -1;
				if (SGPNumsNewToOld[(i + 1) % SGPNumsNewToOld.size()] < TermCPNums.size()) TermCPNum = TermCPNums[SGPNumsNewToOld[(i + 1) % SGPNumsNewToOld.size()]];
				if (TermCPNum < 0 || CPs.GetTypeFromTotOffset(TermCPNum) == PrimaryTermType){
					SurfaceGPs.back().push_back(SGPsPerCP[iCP][SGPNumsNewToOld[(i + 1) % SGPNumsNewToOld.size()]]);
					SurfaceGPCPNums.back().push_back(SGPsPerCP[iCP][SGPNumsNewToOld[(i + 1) % SGPNumsNewToOld.size()]].GetStartEndCPNum(EndCPNumforName));
				}
				else{
					GradPath_c GP = SGPsPerCP[iCP][SGPNumsNewToOld[(i + 1) % SGPNumsNewToOld.size()]];
// #ifdef _DEBUG
// 					GP.SaveAsOrderedZone("Test SGP", Yellow_C);
// 					SaddleCPSupplementGPs[SGPNumsNewToOld[(i + 1) % SGPNumsNewToOld.size()]][1].SaveAsOrderedZone("Test Supp GP", Black_C);
// #endif

					GP.ConcatenateResample(SaddleCPSupplementGPs[SGPNumsNewToOld[(i + 1) % SGPNumsNewToOld.size()]][1], NumGPPts);
					SurfaceGPCPNums.back().push_back(SaddleCPSupplementGPs[SGPNumsNewToOld[(i + 1) % SGPNumsNewToOld.size()]][0].GetStartEndCPNum(0));

// #ifdef _DEBUG
// 					GP.SaveAsOrderedZone("Concatenated GP back", Blue_C);
// #endif

					SurfaceGPs.back().push_back(GP);
				}
			}
		}
		else if (NumSGPs > 0){
			/*
			 *	This is a special case where there is only 1 SGP.
			 *	This means the SGP necessarily connects a bond to a ring (or vice versa)
			 *	so we'll simply make two SZFSs, each with half of the circle.
			 */
			AngleStep = PI2 / static_cast<double>(NumCircleGPs);
			for (int iHalf = 0; iHalf < 2; ++iHalf){
				SurfaceGPCPNums.push_back(vector<int>({ GPParams.StartCPNum }));
				SurfaceGPCPNums.back().push_back(SGPsPerCP[iCP][0].GetStartEndCPNum(EndCPNumforName));
				SurfaceGPs.push_back(vector<GradPath_c>());

				if (iHalf == 0){
					GradPath_c GP = SGPsPerCP[iCP][0];
					//DEBUG
					// 					GP.SaveAsOrderedZone("Test SGP", Yellow_C);
					// 					SaddleCPSupplementGPs[SGPNumsNewToOld[i]][0].SaveAsOrderedZone("Test Supp GP", Black_C);
					//END DEBUG
					GP.ConcatenateResample(SaddleCPSupplementGPs[0][0], NumGPPts);
// 					SurfaceGPCPNums.back().push_back(SaddleCPSupplementGPs[0][0].GetStartEndCPNum(0));

					//DEBUG
					// 					GP.SaveAsOrderedZone("Concatenated GP front", Blue_C);
					//END DEBUG

					SurfaceGPs.back().push_back(GP);
				}

				int jLow = 1 + iHalf * ((NumCircleGPs - 1) / 2);
				int jHigh = jLow + ((NumCircleGPs - 1) / 2);
				for (int j = jLow; j <= jHigh; ++j){
					double RotAngle = SGPSeedAngles[0] + AngleStep * static_cast<double>(j);
					StartPoint = GPParams.StartPointOrigin + Rotate(StartVec, RotAngle, GPParams.RotAxis);
					GradPath_c GP(StartPoint, GPDir, NumGPPts, GPType_Classic, GPTerminate_AtRhoValue, NULL, &CPs, &TermRadius, &RhoCutoff, VolInfo, HessPtrs, GradPtrs, RhoPtr);
					GP.SetStartEndCPNum(CPs.GetTotOffsetFromTypeNumOffset(TypeInd, cpNum), 0);
					GP.Seed(false);
					GP.PointPrepend(GPParams.StartPointOrigin, CPs.GetRho(TypeInd, cpNum));
					GP.Resample(NumGPPts);
					if (GPDir == StreamDir_Reverse) GP.Reverse();
					SurfaceGPs.back().push_back(GP);
				}

				if (iHalf == 1){
					GradPath_c GP = SGPsPerCP[iCP][0];
					//DEBUG
					// 					GP.SaveAsOrderedZone("Test SGP", Yellow_C);
					// 					SaddleCPSupplementGPs[SGPNumsNewToOld[i]][0].SaveAsOrderedZone("Test Supp GP", Black_C);
					//END DEBUG
					GP.ConcatenateResample(SaddleCPSupplementGPs[0][1], NumGPPts);
					// 					SurfaceGPCPNums.back().push_back(SaddleCPSupplementGPs[0][0].GetStartEndCPNum(0));

					//DEBUG
					// 					GP.SaveAsOrderedZone("Concatenated GP front", Blue_C);
					//END DEBUG

					SurfaceGPs.back().push_back(GP);
				}
			}
		}
		else{
			SurfaceGPs.push_back(vector<GradPath_c>());
			AngleStep = PI2 / static_cast<double>(NumCircleGPs);
			for (int j = 0; j < NumCircleGPs; ++j){
				double RotAngle = AngleStep * static_cast<double>(j);
				StartPoint = GPParams.StartPointOrigin + Rotate(StartVec, RotAngle, GPParams.RotAxis);
				GradPath_c GP(StartPoint, GPDir, NumGPPts, GPType_Classic, GPTerminate_AtRhoValue, NULL, &CPs, &TermRadius, &RhoCutoff, VolInfo, HessPtrs, GradPtrs, RhoPtr);
				GP.SetStartEndCPNum(CPs.GetTotOffsetFromTypeNumOffset(TypeInd, cpNum), 0);
				GP.Seed(false);
				GP.PointPrepend(GPParams.StartPointOrigin, CPs.GetRho(TypeInd, cpNum));
				GP.Resample(NumGPPts);
				if (GPDir == StreamDir_Reverse) GP.Reverse();
				SurfaceGPs.back().push_back(GP);
			}
		}

		/*
		* Save special gradient paths
		*/

		AngleStep = 2. * PI / static_cast<int>(NumCircleGPs);

		for (auto & GP : SGPsPerCP[iCP]){
			vector<int> StartEndCPNums = GP.GetStartEndCPNum();
			vector<vector<int> > StartEndCPTypeAndOffset(2);
			string Name = "SGP_";
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
				AuxDataZoneSetItem(GP.GetZoneNum(), CSMAuxData.CC.GPEndNumStrs[i], to_string(StartEndCPNums[i] + 1));
				if (StartEndCPNums[i] >= 0) AuxDataZoneSetItem(GP.GetZoneNum(), CSMAuxData.CC.GPEndTypes[i], CPNameList[StartEndCPTypeAndOffset[i][0]]);
				else AuxDataZoneSetItem(GP.GetZoneNum(), CSMAuxData.CC.GPEndTypes[i], "FF");
			}
			AuxDataZoneSetItem(GP.GetZoneNum(), CSMAuxData.CC.ZoneType, CSMAuxData.CC.ZoneTypeGP);
		}

		/*
		 *	Now make the FE surfaces.
		 *	One for each SGP and one big total surface
		 */
		vector<FESurface_c> Surfaces;
		FESurface_c WholeSurface;
		vector<vector<GradPath_c*> > SurfaceGPPtrs;
		bool SingleSurface = (SurfaceGPs.size() == 1);
		string ZoneNameBase = "SZFS_",
			WholeSurfaceZoneNameBase = (CPType == CPType_BondCP ? "IAS_" : "RS_"),
			CPNameBase = CPNameList[TypeInd] + " " + to_string(cpNum + 1);
		if (SingleSurface) ZoneNameBase = WholeSurfaceZoneNameBase;

		for (int i = 0; i < SurfaceGPs.size(); ++i){
			auto * GPVec = &SurfaceGPs[i];
			SurfaceGPPtrs.push_back(vector<GradPath_c*>());
			for (auto & GP : *GPVec){
				SurfaceGPPtrs.back().push_back(&GP);
			}

// 			SurfaceGPPtrs.back().push_back(&SurfaceGPs[(i + 1) % SurfaceGPs.size()][0]);

// 			for (GradPath_c * GP : SurfaceGPPtrs.back()){
// 				GP->SaveAsOrderedZone("TestGP", Green_C);
// 			}

			Surfaces.push_back(FESurface_c());
			Surfaces.back().MakeFromGPs(SurfaceGPPtrs.back(), SingleSurface);

			//DEBUG
// 			if (i == 2){
// 				for (int ii = 0; ii < SurfaceGPs[i].size() - 1; ++ii){
// 					vector<GradPath_c*> TmpGPPtrs;
// 					TmpGPPtrs.push_back(&SurfaceGPs[i][ii]);
// 					TmpGPPtrs.push_back(&SurfaceGPs[i][ii + 1]);
// 					FESurface_c Surf;
// 					Surf.MakeFromGPs(TmpGPPtrs, false);
// 					FESurface_c::SetZoneStyle(Surf.SaveAsTriFEZone(XYZVarNums, string("TmpSurface " + to_string(ii+1))),
// 						AssignOp_MinusEquals);
// 				}
// 			}

			string SZFSStr = CPNameBase;
			vector<vector<int> > CornerCPTypeOffsets;

			WholeSurface += Surfaces.back();

			if (!SingleSurface){
				SZFSStr += MakeStringFromCPNums(SurfaceGPCPNums[i], CPs, CornerCPTypeOffsets);
				int SurfZoneNum = FESurface_c::SetZoneStyle(Surfaces.back().SaveAsTriFEZone(XYZVarNums, ZoneNameBase + SZFSStr),
					AssignOp_MinusEquals);

				for (int c = 0; c < CornerCPTypeOffsets.size(); ++c){
					AuxDataZoneSetItem(SurfZoneNum, CSMAuxData.CC.ZFSCornerCPNumStrs[c], to_string(CornerCPTypeOffsets[c][1] + 1));
					if (CornerCPTypeOffsets[c][0] >= 0) AuxDataZoneSetItem(SurfZoneNum, CSMAuxData.CC.ZFSCornerCPTypes[c], CPNameList[CornerCPTypeOffsets[c][0]]);
					else AuxDataZoneSetItem(SurfZoneNum, CSMAuxData.CC.ZFSCornerCPTypes[c], "FF");
				}
				AuxDataZoneSetItem(SurfZoneNum, CSMAuxData.CC.ZoneType, CSMAuxData.CC.ZoneTypeSZFS);
				AuxDataZoneSetItem(SurfZoneNum, CSMAuxData.CC.ZoneSubType, (CPType == CPType_RingCP ? CSMAuxData.CC.ZoneSubTypeRS : CSMAuxData.CC.ZoneSubTypeIAS));
			}
		}

		if (WholeSurface.IsMade()){
			int SurfZoneNum = FESurface_c::SetZoneStyle(WholeSurface.SaveAsTriFEZone(XYZVarNums, WholeSurfaceZoneNameBase + CPNameBase),
				AssignOp_PlusEquals, TRUE, FALSE);

			AuxDataZoneSetItem(SurfZoneNum, CSMAuxData.CC.ZFSCornerCPNumStrs[0], to_string(CPs.GetTotOffsetFromTypeNumOffset(TypeInd, cpNum) + 1));
			AuxDataZoneSetItem(SurfZoneNum, CSMAuxData.CC.ZFSCornerCPTypes[0], CPNameList[TypeInd]);
			AuxDataZoneSetItem(SurfZoneNum, CSMAuxData.CC.ZoneType, CSMAuxData.CC.ZoneTypeSZFS);
			AuxDataZoneSetItem(SurfZoneNum, CSMAuxData.CC.ZoneSubType, (CPType == CPType_RingCP ? CSMAuxData.CC.ZoneSubTypeRS : CSMAuxData.CC.ZoneSubTypeIAS));
		}
	}

// 	/*
// 	 * Save special gradient paths
// 	 */
// 
// 	AngleStep = 2. * PI / static_cast<int>(NumCircleGPs);
// 
// 	for (int iCP = 0; iCP < NumCPs; ++iCP){
// 		for (auto & GP : SGPsPerCP[iCP]){
// 			vector<int> StartEndCPNums = GP.GetStartEndCPNum();
// 			vector<vector<int> > StartEndCPTypeAndOffset(2);
// 			string Name = "SGP_";
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
// 		}
// 	}

	StatusDrop(AddOnID);

	//DEBUG
	TecUtilDataLoadEnd();
	
	CSMGuiUnlock();
	TecUtilLockFinish(AddOnID);
	return TRUE;


	vector<GradPath_c> GPs;

	/*
	 *	Shotgun approach.
	 *	First seed GPs around each CP along a constant radius circle in the 
	 *	plane normal to the CP principal direction.
	 *	These GPs will identify all bond-cage, ring-nuclear, and bond-ring connections.
	 */

	StartPointOffset = 0.05 * CPs.GetMinCPDist();

	for (int cpNum = 0; cpNum < NumCPs; ++cpNum){
		vec3 StartVec = normalise(CPs.GetEigVecs(TypeInd, cpNum).col(1)) * StartPointOffset;

		for (int gpNum = 0; gpNum < NumCircleCheckGPs; ++gpNum){
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
		string Name = "GP_";
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
			AuxDataZoneSetItem(GP.GetZoneNum(), CSMAuxData.CC.GPEndNumStrs[i], to_string(StartEndCPNums[i] + 1));
			if (StartEndCPNums[i] >= 0) AuxDataZoneSetItem(GP.GetZoneNum(), CSMAuxData.CC.GPEndTypes[i], CPNameList[StartEndCPTypeAndOffset[i][0]]);
			else AuxDataZoneSetItem(GP.GetZoneNum(), CSMAuxData.CC.GPEndTypes[i], "FF");
		}
		AuxDataZoneSetItem(GP.GetZoneNum(), CSMAuxData.CC.ZoneType, CSMAuxData.CC.ZoneTypeGP);
	}

	TecUtilDataLoadEnd();

	
	CSMGuiUnlock();
	TecUtilLockFinish(AddOnID);

	return TRUE;
}

void ConnectCPsReturnUserInfo(const bool GuiSuccess, const vector<GuiField_c> & Fields){
	TecUtilLockStart(AddOnID);

	vector<int> XYZVarNums(3);

	int fNum = 0;

	for (int i = 0; i < 3; ++i) XYZVarNums[i] = Fields[fNum++].GetReturnInt();
	fNum++;

	int NumPts = Fields[fNum++].GetReturnInt(),
		CPTypeVarNum = Fields[fNum++].GetReturnInt(),
		CPZoneNum = Fields[fNum].GetReturnInt();
	vector<int> CPNums = Fields[fNum].GetReturnIntVec();

	TecUtilDataLoadBegin();

	CritPoints_c CPs(CPZoneNum, XYZVarNums, CPTypeVarNum);

	Set_pa ZoneSet = TecUtilSetAlloc(FALSE);

	for (int i = 0; i < CPNums.size() - 1; ++i){
		for (int j = i + 1; j < CPNums.size(); ++j){
			vec3 a = CPs.GetXYZ(CPNums[i]-1),
				s = (CPs.GetXYZ(CPNums[j]-1) - a) / double(NumPts - 1);

			vector<vector<int> > StartEndCPTypeAndOffset(2);
			string ZoneName = "Line " + MakeStringFromCPNums({ CPNums[i]-1, CPNums[j]-1 }, CPs, StartEndCPTypeAndOffset);
			Boolean_t IsOk = TecUtilDataSetAddZone(ZoneName.c_str(), NumPts, 1, 1, ZoneType_Ordered, NULL);
			FieldVecPointer_c XYZPtr;
			IsOk = IsOk && XYZPtr.GetWritePtr(TecUtilDataSetGetNumZones(), XYZVarNums);
			if (IsOk){
				for (int n = 0; n < NumPts; ++n) XYZPtr.Write(n, a + s * double(n));
			}
			TecUtilSetAddMember(ZoneSet, TecUtilDataSetGetNumZones(), FALSE);
		}
	}

	TecUtilDataLoadEnd();

	TecUtilZoneSetScatter(SV_SHOW, ZoneSet, 0.0, FALSE);
	TecUtilZoneSetMesh(SV_SHOW, ZoneSet, 0.0, TRUE);
	TecUtilZoneSetMesh(SV_COLOR, ZoneSet, 0.0, Black_C);
	TecUtilZoneSetActive(ZoneSet, AssignOp_PlusEquals);

	TecUtilLockFinish(AddOnID);
}

void ConnectCPsGetUserInfo(){
	vector<GuiField_c> Fields = {
		GuiField_c(Gui_VarSelect, "X", "X"),
		GuiField_c(Gui_VarSelect, "Y", "Y"),
		GuiField_c(Gui_VarSelect, "Z", "Z"),
		GuiField_c(Gui_VertSep),
		GuiField_c(Gui_Int, "Number of points in lines", "100"),
		GuiField_c(Gui_VarSelect, "CP Type", CSMVarName.CritPointType),
		GuiField_c(Gui_ZonePointSelectMulti, "CPs", CSMZoneName.CriticalPoints)
	};

	CSMGui("Connect CPs with lines", Fields, ConnectCPsReturnUserInfo, AddOnID);
}

void DrawEigenvectotArrorsReturnUserInfo(const bool GuiSuccess, const vector<GuiField_c> & Fields){
	TecUtilLockStart(AddOnID);

	int fNum = 0;

	bool IsOk = true;

	bool InverseLength = Fields[fNum++].GetReturnBool();
	int SkipNum = abs(Fields[fNum++].GetReturnInt());
	bool AllZones = Fields[fNum++].GetReturnBool(),
		AllSurfaces = Fields[fNum++].GetReturnBool(),
		AllLines = Fields[fNum++].GetReturnBool(),
		OnlySelectedZones = Fields[fNum++].GetReturnBool();
	vector<int> SelectedZones = Fields[fNum++].GetReturnIntVec();
	int RhoVarNum = Fields[fNum++].GetReturnInt();

	vector<int> XYZVarNums(3), EigValVarNums(3);
	vector<vector<int> > EigVecVarNums(3, vector<int>(3));

	for (int i = 0; i < 3; ++i) XYZVarNums[i] = Fields[fNum++].GetReturnInt();
	for (int i = 0; i < 3; ++i) EigValVarNums[i] = Fields[fNum++].GetReturnInt();
	for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j) EigVecVarNums[i][j] = Fields[fNum++].GetReturnInt();


	int NumZones = TecUtilDataSetGetNumZones();
	int NumVars = TecUtilDataSetGetNumVars();
	vector<int> IJK(3);

	vector<FieldDataPointer_c> RhoPtrs;
	vector<FieldVecPointer_c> XYZPtrs;
	vector<vector<FieldVecPointer_c> > EigVecPtrs;
	vector<vector<FieldDataPointer_c> > EigValPtrs;

	TecUtilDataLoadBegin();

	if (!OnlySelectedZones && AllZones || AllSurfaces || AllLines ){
		SelectedZones.clear();
		for (int z = 1; z <= NumZones; ++z){
			bool RunZone = true;
			ZoneType_e t = TecUtilZoneGetType(z);
			TecUtilZoneGetIJK(z, &IJK[0], &IJK[1], &IJK[2]);
			if (ZoneType_Ordered != t || IJK[2] > 1){
				// If finite element or (1,2)d ordered zone
				RunZone = (
					AllZones
					|| (AllLines && (ZoneType_FELineSeg == t || (ZoneType_Ordered == t && IJK[1] == 1)))
					|| (AllSurfaces && (ZoneType_FETriangle == t || ZoneType_FEQuad == t || ZoneType_FEPolygon == t || ZoneType_Ordered))
					);
			}
			if (RunZone){
				SelectedZones.push_back(z);
			}
		}
	}

	int NumSelectedZones = SelectedZones.size();

	for (const int & z : SelectedZones){
		RhoPtrs.push_back(FieldDataPointer_c());
		if (!RhoPtrs.back().GetReadPtr(z, RhoVarNum)){
			TecUtilDialogErrMsg("Failed to get rho data pointer");
			return;
		}
		XYZPtrs.push_back(FieldVecPointer_c());
		if (!XYZPtrs.back().GetReadPtr(z, XYZVarNums)){
			TecUtilDialogErrMsg("Failed to get XYZ data pointer");
			return;
		}
		EigValPtrs.push_back(vector<FieldDataPointer_c>());
		for (const int & v : EigValVarNums){
			EigValPtrs.back().push_back(FieldDataPointer_c());
			if (!EigValPtrs.back().back().GetReadPtr(z, v)){
				TecUtilDialogErrMsg("Failed to get eigenvalue data pointer");
				return;
			}
		}
		EigVecPtrs.push_back(vector<FieldVecPointer_c>());
		for (const vector<int> & e : EigVecVarNums) {
			EigVecPtrs.back().push_back(FieldVecPointer_c());
			if (!EigVecPtrs.back().back().GetReadPtr(z, e)){
				TecUtilDialogErrMsg("Failed to get eigenvector data pointer");
				return;
			}
		}
	}

	vector<vector<FESurface_c> > ArrowZones(NumSelectedZones, vector<FESurface_c>(3));

#ifndef _DEBUG
#pragma omp parallel for
#endif
	for (int z = 0; z < NumSelectedZones; ++z){
		int NumPoints = EigValPtrs[z][0].Size();
		int NumArrows = NumPoints / SkipNum;
		vector<vector<FESurface_c> > AllArrows(3);
		vector<vec3> Nodes;
		vector<vector<int> > Elements;
		for (int n = 0; n < NumArrows; ++n){
			int ni = MIN(n * SkipNum, NumPoints - 1);
			for (int ei = 0; ei < 3; ++ei){
				if (CSMArrow(XYZPtrs[z][ni], EigVecPtrs[z][ei][ni], log(fabs(EigValPtrs[z][ei][ni]) / RhoPtrs[z][ni] + 1), 0.1, 0.2, 0.5, Nodes, Elements)){
					AllArrows[ei].push_back(FESurface_c(Nodes, Elements));
					if (ei == 0) AllArrows[ei].back().SaveAsTriFEZone(XYZVarNums, "TestArrow");
				}
			}
// 			TecUtilDataLoadEnd();
// 			TecUtilLockFinish(AddOnID);
// 			return;
		}
		for (int ei = 0; ei < 3; ++ei) for (int n = 0; n < AllArrows[ei].size(); ++n){
			ArrowZones[z][ei] += AllArrows[ei][n];
		}
	}


	TecUtilDataLoadEnd();

	vector<Set_pa> ArrowZoneSets(3);
	for (Set_pa & s : ArrowZoneSets) s = TecUtilSetAlloc(FALSE);
	Set_pa AllZonesSet = TecUtilSetAlloc(FALSE);
	for (int z = 0; z < ArrowZones.size(); ++z){
		char *BaseZoneNameCStr;
		TecUtilZoneGetName(XYZPtrs[z].ZoneNum(), &BaseZoneNameCStr);

		string NewZoneName = BaseZoneNameCStr + string(" (Eigvec ");
		TecUtilStringDealloc(&BaseZoneNameCStr);

		for (int ei = 0; ei < 3; ++ei){
			ArrowZones[z][ei].SaveAsTriFEZone(XYZVarNums, (NewZoneName + to_string(ei + 1) + " arrow)").c_str());
			TecUtilSetAddMember(AllZonesSet, ArrowZones[z][ei].GetZoneNum(), FALSE);
			TecUtilSetAddMember(ArrowZoneSets[ei], ArrowZones[z][ei].GetZoneNum(), FALSE);
		}
	}

	TecUtilSetDealloc(&AllZonesSet);
	for (Set_pa & s : ArrowZoneSets) TecUtilSetDealloc(&s);

	TecUtilLockFinish(AddOnID);
}

void DrawEigenvectorArrowsGetUserInfo(){
	vector<GuiField_c> Fields = {
		GuiField_c(Gui_Toggle, "Use inverse eigenvalue magnitudes for length", "0"),
		GuiField_c(Gui_Int, "Draw arrow every n points", "1"),
		GuiField_c(Gui_Toggle, "All:", "0"),
		GuiField_c(Gui_Toggle, "Surfaces"),
		GuiField_c(Gui_Toggle, "Lines", "0"),
		GuiField_c(Gui_Toggle, "Selected", "0"),
		GuiField_c(Gui_ZoneSelectMulti, "Zones", ""),
		GuiField_c(Gui_VarSelect, CSMVarName.Dens, CSMVarName.Dens)
	};

	for (const string & s : { "X", "Y", "Z" }) Fields.push_back(GuiField_c(Gui_VarSelect, s, s));
	for (const string & s : CSMVarName.EigVals) Fields.push_back(GuiField_c(Gui_VarSelect, s, s));
	for (const string & s : CSMVarName.EigVecs) Fields.push_back(GuiField_c(Gui_VarSelect, s, s));

	CSMGui("Connect CPs with lines", Fields, DrawEigenvectotArrorsReturnUserInfo, AddOnID);
}

void TestFunction(){

	TecUtilZoneSetActive(Set(1).getRef(), AssignOp_Equals);
		
}


const bool ExtractRadiusContourLinesToIOrderedPoints(const vector<int> & ZoneNums, 
	const vec3 & Origin, const double & Radius, const vector<int> & XYZVarNums, vector<GradPath_c> & ContourLines)
{
	TecUtilLockStart(AddOnID);

	int NumZones = TecUtilDataSetGetNumZones();

	// First define a RadiusSquared variable for the specified zones

	Set ZoneSet;
	ArgList Args;

	for (const int & z : ZoneNums){
		REQUIRE(z > 0 && z <= NumZones);
		ZoneSet += z;
	}

	string TempVarName = "TempRadiusSqrVar";
	int TempVarNum = VarNumByName(TempVarName);

	if (TempVarNum <= 0){
		// Need to make new temp variable
		vector<FieldDataType_e> FDTypes(NumZones, FieldDataType_Bit);
		for (const int & z : ZoneNums) FDTypes[z - 1] = FieldDataType_Double;
		Args.clear();
		Args.appendString(SV_NAME, TempVarName);
		Args.appendArray(SV_VARDATATYPE, FDTypes.data());
		Args.appendInt(SV_DEFERVARCREATION, TRUE);

		TecUtilDataSetAddVarX(Args.getRef());
	}

	TempVarNum = VarNumByName(TempVarName);
	if (TempVarNum <= 0){
		TecUtilDialogErrMsg("Failed to make temp radius variable");
		return false;
	}

	string DataAlterEqn = "V" + to_string(TempVarNum) + "=";
	for (int i = 0; i < 3; ++i)
		DataAlterEqn += "(V" + to_string(XYZVarNums[i]) + " - " + to_string(Origin[i]) + ")**2" + (i<2 ? "+" : "");

	TecUtilDataAlter(DataAlterEqn.c_str(), ZoneSet.getRef(),
		1, 0, 1,
		1, 0, 1,
		1, 0, 1,
		FieldDataType_Double);


	// Then loop through the zones, for each extracting a contour at the specified radius to a FELineSeg zone.
	// This will be done in a new frame to prevent messing with the contours in the current frame
	
	// Get frame id of current frame
	UniqueID_t OldFrame = TecUtilFrameGetUniqueID(), NewFrame;

	// Make new frame and get its id
	TecUtilFrameCreateNew(TRUE, 1, 1, 1, 1);
	NewFrame = TecUtilFrameGetUniqueID();

	// Set new frame plot type to 3d cartesian
	Args.clear();
	Args.appendString(SV_P1, SV_PLOTTYPE);
	Args.appendArbParam(SV_IVALUE, PlotType_Cartesian3D);
	TecUtilStyleSetLowLevelX(Args.getRef());

	// Set the first contour group to the new variable

	Args.clear();
	Args.appendInt(SV_CONTOURGROUP, 1);
	Args.appendInt(SV_VAR, TempVarNum);
	TecUtilContourSetVariableX(Args.getRef());

	// Set the first contour group to have a single contour at the specified radius

	double RadiusSqr = Radius * Radius;
	Args.clear();
	Args.appendInt(SV_CONTOURLEVELACTION, ContourLevelAction_New);
	Args.appendInt(SV_NUMVALUES, 1);
	Args.appendArray(SV_RAWDATA, (void*)&RadiusSqr);
	TecUtilContourLevelX(Args.getRef());

	// Turn on contours and set group for zones
	TecUtilZoneSetContour(SV_SHOW, ZoneSet.getRef(), 0.0, TRUE);
	StyleValue Style;
	Style.set((SmInteger_t)1, ZoneSet, SV_FIELDMAP, SV_CONTOUR, SV_LINECONTOURGROUP);
	Style.set(ContourType_Lines, ZoneSet, SV_FIELDMAP, SV_CONTOUR, SV_CONTOURTYPE); 
	Style.set(TRUE, SV_FIELDLAYERS, SV_SHOWCONTOUR);

// 	TecUtilArgListDealloc(&argList);

	// Now actually loop through zones to extract contour lines
	ZoneSet.clear();

	for (const int & z : ZoneNums){
		TecUtilZoneSetActive(Set(z).getRef(), AssignOp_Equals);

		Args.clear();
		Args.appendInt(SV_CONTLINECREATEMODE, ContLineCreateMode_OneZonePerIndependentPolyline);
		TecUtilCreateContourLineZonesX(Args.getRef());

		int NewZoneNum = TecUtilDataSetGetNumZones();

		if (NewZoneNum - NumZones == 1){
			ContourLines.push_back(GradPath_c(NewZoneNum, XYZVarNums, AddOnID));

			ZoneSet += NewZoneNum;
			NumZones = NewZoneNum;
		}
	}

	TecUtilDataSetDeleteVar(Set(TempVarNum).getRef());
	TecUtilDataSetDeleteZone(ZoneSet.getRef());
	TecUtilFrameDeleteActive();

	TecUtilLockFinish(AddOnID);

	return true;
}

void ExtractRSIntersectionsReturnUserInfo(const bool GuiSuccess, const vector<GuiField_c> & Fields){
	int fNum = 0;

	int VolZoneNum = Fields[fNum++].GetReturnInt(),
		CPTypeVarNum = Fields[fNum++].GetReturnInt();
	int AllCPZoneNum = Fields[fNum++].GetReturnInt();
	vector<int> CPList = Fields[fNum].GetReturnIntVec();
	int CPZoneNum = Fields[fNum++].GetReturnInt();
	double Radius = Fields[fNum++].GetReturnDouble();
	int ResampleNumPoints = Fields[fNum++].GetReturnInt();
	vector<int> XYZVarNums;
	for (int i = 0; i < 3; ++i) XYZVarNums.push_back(Fields[fNum++].GetReturnInt());

	char* ZoneName;
	TecUtilZoneGetName(CPZoneNum, &ZoneName);
	int CPZoneTypeInd = VectorGetElementNum(CSMZoneName.CPType, string(ZoneName));
	// Will be -1 if CPs selected from total CP zone

	if (CSMZoneName.CriticalPoints == ZoneName && CPZoneTypeInd != 0 && CPZoneTypeInd != 3){
		// Check that the user selected points from the total, nuclear, or cage critical points zones
		TecUtilDialogErrMsg("Must select nuclear or cage critical points");
		TecUtilStringDealloc(&ZoneName);
		return;
	}
	TecUtilStringDealloc(&ZoneName);

	TecUtilLockStart(AddOnID);

	TecUtilDrawGraphics(FALSE);

	int NumZones = TecUtilDataSetGetNumZones();

	CPType_e CPZoneType = CPTypeList[CPZoneTypeInd];

	CritPoints_c CPs(AllCPZoneNum, XYZVarNums, CPTypeVarNum);

	vector<string> SurfaceSearchStrings = { CSMAuxData.CC.ZoneSubTypeRS, CSMAuxData.CC.ZoneSubTypeIAS };
	vector<string> PathSearchStrings = { CSMAuxData.CC.ZoneSubTypeBondPath, CSMAuxData.CC.ZoneSubTypeRingLine };

	for (const int & CPNum1 : CPList){
		vector<int> TypeNumOffset = (CPZoneTypeInd < 0 ? CPs.GetTypeNumOffsetFromTotOffset(CPNum1 - 1) : vector<int>({ CPZoneTypeInd, CPNum1 - 1 }));
		CPType_e CPType = CPTypeList[TypeNumOffset[0]];
		if (CPType != CPType_NuclearCP && CPType != CPType_CageCP) continue;

		int CPTotOffset = CPs.GetTotOffsetFromTypeNumOffset(TypeNumOffset[0], TypeNumOffset[1]) + 1;

		// Get zone numbers of intersecting ring (interatomic) surfaces
		string SurfaceSearchString = (CPType == CPType_NuclearCP ? CSMAuxData.CC.ZoneSubTypeRS : CSMAuxData.CC.ZoneSubTypeIAS);
		vector<int> IntersectingZoneNums;
		for (int z = 1; z <= NumZones; ++z){
			if (AuxDataZoneItemMatches(z, CSMAuxData.CC.ZoneSubType, SurfaceSearchStrings[CPType == CPType_NuclearCP ? 0 : 1])){
				for (const string & CornerStr : CSMAuxData.CC.ZFSCornerCPNumStrs){
					int CornerNum = stoi(AuxDataZoneGetItem(z, CornerStr));
					if (CornerNum == CPTotOffset){
						IntersectingZoneNums.push_back(z);
						break;
					}
				}
			}
			else if (AuxDataZoneItemMatches(z, CSMAuxData.CC.ZoneSubType, PathSearchStrings[CPType == CPType_NuclearCP ? 0 : 1])){
				for (const string & CornerStr : CSMAuxData.CC.GPEndNumStrs){
					int CornerNum = stoi(AuxDataZoneGetItem(z, CornerStr));
					if (CornerNum == CPTotOffset){
						IntersectingZoneNums.push_back(z);
						break;
					}
				}
			}
		}

		// Extract intersection lines as gradient paths
		vector<GradPath_c> Intersections;
		if (IntersectingZoneNums.size() > 0){
			ExtractRadiusContourLinesToIOrderedPoints(IntersectingZoneNums, CPs.GetXYZ(TypeNumOffset[0], TypeNumOffset[1]), Radius, XYZVarNums, Intersections);
		}

		// Save intersections as zones
		for (int iGP = 0; iGP < Intersections.size(); ++iGP){
			Intersections[iGP].Resample(ResampleNumPoints);
			Intersections[iGP].SaveAsOrderedZone(string("Intersection of CP " + to_string(CPNum1) + " w/ Zone " + to_string(IntersectingZoneNums[iGP])));
		}
	}

	TecUtilDrawGraphics(TRUE);

	TecUtilLockFinish(AddOnID);
}

void ExtractRSIntersectionsGetUserInfo(){
	vector<GuiField_c> Fields({
		GuiField_c(Gui_ZoneSelect,				CSMZoneName.FullVolume,			CSMZoneName.FullVolume),
		GuiField_c(Gui_VarSelect,				"CP type variable",				CSMVarName.CritPointType),
		GuiField_c(Gui_ZoneSelect,				"Critical points zone",			CSMZoneName.CriticalPoints),
		GuiField_c(Gui_ZonePointSelectMulti,	"Nuclear CP zone",				CSMZoneName.CPType[0]),
		GuiField_c(Gui_Double,					"Intersection Radius",			"0.2"),
		GuiField_c(Gui_Int,						"n points per intersection",	"20")
	});

	for (const string & s : { "X", "Y", "Z" }) Fields.push_back(GuiField_c(Gui_VarSelect, s, s));

	CSMGui("Extract Ring Surface Sphere Intersections", Fields, ExtractRSIntersectionsReturnUserInfo, AddOnID);
}

void GBAGetUserInfo(){
	
}