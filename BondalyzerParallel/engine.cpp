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
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multimin.h>

#include "CSM_DATA_TYPES.h"
#include "CSM_DATA_SET_INFO.h"
#include "CSM_CRIT_POINTS.h"
#include "CSM_CALC_VARS.h"
#include "CSM_GRAD_PATH.h"
#include "CSM_FE_VOLUME.h"
#include "CSM_GUI.h"
#include "CSM_GEOMETRY.h"

#include "updateSphericalTriangulation.h"

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
	vector<FieldDataPointer_c> const * HessPtrs = NULL;
	vector<FieldDataPointer_c> const * GradPtrs = NULL;
	FieldDataPointer_c const * RhoPtr = NULL;

	GradPath_c GP;

	vec3 StartPointOrigin;
	vec3 RotVec;
	vec3 RotAxis;
	int StartCPNum = -1;
	int EndCPPosition = -1;
	int EndCPNum = -1;
};

struct MinFuncParams_GPLengthSpherical{
	StreamDir_e Direction = StreamDir_Invalid;
	int NumGPPoints = -1;
	GPType_e GPType = GPType_Invalid;
	GPTerminate_e GPTerminate = GPTerminate_Invalid;
	vec3 * TermPoint = NULL;
	CritPoints_c * CPs = NULL;
	double * TermPointRadius = NULL;
	double * TermValue = NULL;
	VolExtentIndexWeights_s * VolInfo = NULL;
	vector<FieldDataPointer_c> const * HessPtrs = NULL;
	vector<FieldDataPointer_c> const * GradPtrs = NULL;
	FieldDataPointer_c const * RhoPtr = NULL;

	GradPath_c GP;

	vec3 StartPointOrigin;
	double r;
	int StartCPNum = -1;
	int EndCPNum = -1;
};


string const DS_ZoneName_Delim = "_-_";
string const T41Ext = ".t41";

static int const MinNumCircleCheckGPs = 1080;
static int const MinNumCircleGPs = 120;
static int const MinNumRCSFuncCircleCheckPts = 3600;
static int const DefaultRCSNumCheckRadii = 128;
static int const DefaultRCSNumCheckRadiiConverged = 4;
static double const DefaultRCSCheckRadiiLow = 0.2;
static double const DefaultRCSCheckRadiiHigh = 0.8;
static double const DefaultRCSCheckRadiiHighMultiplier = 1.4; // if <= 1 near field SGPs are found, increase the radius for RCSFunc-based search
static int const DefaultRCSFuncCheckNumPts = 3; 
static int const DefaultRCSFuncConvergence = 32;
static int const DefaultRCSSGPAngleCheck = 15; // If RCS-based SGP is less than this angle [degrees] from an existing SGP it is discarded
static double const SmallAngleFactor = 0.5;
static int const MaxIter_GPLengthInPlane = 100;
static double const SurfRCSMinAngleCheck = 1e-6;

//DEBUG
#ifdef _DEBUG
vector<vector<double> > MinFunc_GPAngleLengths;
#endif

/*
 *	This is the convergence criteria for a one-dimensional GSL minimization
 *	based on gradient path length.
 */
static double const Tolerance_GPLengthInPlane = 1e-5;

/*
 *	If the terminal ends of two gradient paths diverge away from eachother
 *	such that the inner (dot) product of their terminal ends is less than
 *	below, they are considered divergent and to be terminating at different 
 *	far field CPs (if they were near field, then we would know which CP and
 *	don't need to resort to checking for divergence).
 */
static double const Tolerance_DivergentGPInnerProd = -0.9;

BondalyzerCalcType_e CurrentCalcType;


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

	CSMGuiLock();

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
			UserQuit = !StatusUpdate(VolNum++, NumVols, "Refining zones", AddOnID);
#pragma omp flush (UserQuit)
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

	CSMGuiUnlock();
}

bool const FEZoneDFS(int NodeNum,
				NodeMap_pa const & NodeMap,
				NodeToElemMap_pa const & NodeToElemMap,
				int NumNodesPerElem,
				vector<bool> & IsVisited){
	bool IsOk = true;
	IsVisited[NodeNum] = true;
	int NumElems = TecUtilDataNodeToElemMapGetNumElems(NodeToElemMap, NodeNum + 1);
	for (int e = 1; e <= NumElems && IsOk; ++e){
		int ei = TecUtilDataNodeToElemMapGetElem(NodeToElemMap, NodeNum + 1, e);
		for (int n = 0; n < NumNodesPerElem && IsOk; ++n){
			int ni = TecUtilDataNodeGetByRef(NodeMap, ei, n + 1) - 1;
			if (!IsVisited[ni]){
				IsOk = FEZoneDFS(ni, NodeMap, NodeToElemMap, NumNodesPerElem, IsVisited);
			}
		}
	}

	return IsOk;
}

// void FEZoneDFS(int NodeNum, vector<vector<int> > const & Elems, vector<bool> & IsVisited){
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

		for (auto const & i : Answers){
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

	for (int PointINum : PointINums){
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

		for (auto const & i : Answers){
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

	TecUtilLockStart(AddOnID);

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

	TecUtilLockFinish(AddOnID);
}

void GetClosedIsoSurface(int IsoZoneNum, vector<FieldDataPointer_c> const & IsoReadPtrs, vector<int> & NodeNums){

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

	for (int NodeNum : NodeNums){
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
		for (auto const & i : NodeVisited) NumNewNodes += (int)i;

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

		for (int e = 0; e < IsoIJK[1]; ++e) if (!ElemAdded[e]) for (auto const & n : Elems[e]) if (NodeVisited[n])
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

		vector<int> NewZoneIJK(3);
		TecUtilZoneGetIJK(NewZoneNum, &NewZoneIJK[0], &NewZoneIJK[1], &NewZoneIJK[2]);

		for (int n = 0; n < NumNewNodes; ++n) for (int i = 0; i < IsoWritePtrs.size(); ++i) IsoWritePtrs[i].Write(n, IsoReadPtrs[i][NodeNumsNewToOld[n]]);

		NodeMap_pa NewNodeMap = TecUtilDataNodeGetWritableRef(NewZoneNum);
		if (!VALID_REF(NewNodeMap)){
			TecUtilDialogErrMsg("Failed to get node map pointer for new iso zone. Quitting.");
			return;
		}

		for (int e = 0; e < NewElems.size(); ++e){
			int e1 = e + 1;
			for (int ei = 0; ei < IsoIJK[2]; ++ei){
				int ei1 = ei + 1;
				int n = NewElems[e][ei % NewElems[e].size()] + 1;
				TecUtilDataNodeSetByZone((EntIndex_t)NewZoneNum, (LgIndex_t)e1, (LgIndex_t)ei1, (NodeMap_t)n);
// 				TecUtilDataNodeSetByRef(NewNodeMap, e1, ei1, n);
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

void BondalyzerReturnUserInfo(bool const GuiSuccess, 
	vector<GuiField_c> const & Fields,
	vector<GuiField_c> const PassthroughFields){
	if (!GuiSuccess) return;

	TecUtilLockStart(AddOnID);

	int VolZoneNum, RhoVarNum;
	vector<int> XYZVarNums(3), GradVarNums, HessVarNums;
	Boolean_t IsPeriodic;

	int fNum = 0;

	VolZoneNum = Fields[fNum++].GetReturnInt();
	IsPeriodic = Fields[fNum++].GetReturnBool();
	fNum++;
	for (int i = 0; i < 3; ++i) XYZVarNums[i] = i + Fields[fNum].GetReturnInt();
	fNum++;
	fNum++;
	RhoVarNum = Fields[fNum++].GetReturnInt();
	fNum++;


	if (Fields[fNum++].GetReturnBool()){
		GradVarNums.resize(3);
		for (int i = 0; i < 3; ++i) GradVarNums[i] = i + Fields[fNum].GetReturnInt();
// 		for (int i = 0; i < 3; ++i) GradVarNums[i] = Fields[fNum++].GetReturnInt();
	}
	fNum++;
// 	else fNum += 3;
	fNum++;

	if (Fields[fNum++].GetReturnBool()){
		HessVarNums.resize(6);
		for (int i = 0; i < 6; ++i) HessVarNums[i] = i + Fields[fNum].GetReturnInt();
// 		for (int i = 0; i < 6; ++i) HessVarNums[i] = Fields[fNum++].GetReturnInt();
	}
	fNum++;
// 	else fNum += 6;

	fNum++;

	if (CurrentCalcType == BondalyzerCalcType_CriticalPoints){
		double CellSpacing = Fields[fNum].GetReturnDouble();
		FindCritPoints(VolZoneNum, XYZVarNums, RhoVarNum, GradVarNums, HessVarNums, IsPeriodic, CellSpacing);
	}
	else if (CurrentCalcType >= BondalyzerCalcType_BondPaths && CurrentCalcType < BondalyzerCalcType_GBA){

		int RidgeFuncVarNum = 0;

		if (CurrentCalcType >= BondalyzerCalcType_InteratomicSurfaces){
			if (Fields[fNum++].GetReturnBool())
				RidgeFuncVarNum = Fields[fNum++].GetReturnInt();
			else fNum++;
		}

		int CPTypeVarNum = Fields[fNum++].GetReturnInt();
		vector<int> OtherCPZoneNums = Fields[fNum++].GetReturnIntVec();
		int SelectCPsZoneNum = Fields[fNum].GetReturnInt();
		vector<int> SelectedCPs = Fields[fNum++].GetReturnIntVec();

		CPType_e CPType;
		if (CurrentCalcType == BondalyzerCalcType_BondPaths || CurrentCalcType == BondalyzerCalcType_InteratomicSurfaces) CPType = CPType_Bond;
		else if (CurrentCalcType == BondalyzerCalcType_RingLines || CurrentCalcType == BondalyzerCalcType_RingSurfaces) CPType = CPType_Ring;

		if (CurrentCalcType == BondalyzerCalcType_BondPaths || CurrentCalcType == BondalyzerCalcType_RingLines)
			FindBondRingLines(VolZoneNum, OtherCPZoneNums, SelectCPsZoneNum, SelectedCPs, CPTypeVarNum, CPType, XYZVarNums, RhoVarNum, GradVarNums, HessVarNums, IsPeriodic);
		else if (CurrentCalcType == BondalyzerCalcType_CageNuclearPaths)
			FindCageNuclearPaths(VolZoneNum, OtherCPZoneNums, SelectCPsZoneNum, SelectedCPs, CPTypeVarNum, XYZVarNums, RhoVarNum, GradVarNums, HessVarNums, IsPeriodic);
		else if (CurrentCalcType == BondalyzerCalcType_InteratomicSurfaces || CurrentCalcType == BondalyzerCalcType_RingSurfaces)
			FindBondRingSurfaces(VolZoneNum, OtherCPZoneNums, SelectCPsZoneNum, SelectedCPs, CPTypeVarNum, CPType, XYZVarNums, RhoVarNum, GradVarNums, HessVarNums, RidgeFuncVarNum, IsPeriodic);
	}
	else if (CurrentCalcType == BondalyzerCalcType_GBA){
		fNum = 0;
		vector<int> CPList = PassthroughFields[fNum].GetReturnIntVec();
		int CPZoneNum = PassthroughFields[fNum++].GetReturnInt(),
			CPTypeVarNum = PassthroughFields[fNum++].GetReturnInt();
		fNum++;
		double SphereRadius = PassthroughFields[fNum++].GetReturnDouble();
		int RadiusMode = PassthroughFields[fNum++].GetReturnInt(),
			RefinementLevel = PassthroughFields[fNum++].GetReturnInt();
		fNum++;
		double RhoCutoff = PassthroughFields[fNum++].GetReturnDouble();
		int NumGBEdgePoints = PassthroughFields[fNum++].GetReturnInt(),
			GPNumPoints = PassthroughFields[fNum++].GetReturnInt();

		GBA_Generation(VolZoneNum, 
			CPZoneNum, 
			CPTypeVarNum, 
			XYZVarNums, 
			RhoVarNum, 
			GradVarNums, 
			HessVarNums, 
			IsPeriodic,
			CPList, 
			SphereRadius,
			RadiusMode, 
			RefinementLevel, 
			RhoCutoff, 
			NumGBEdgePoints, 
			GPNumPoints);
	}
	else if (CurrentCalcType == BondalyzerCalcType_Batch){
		bool UseCPZones = Fields[fNum++].GetReturnBool();
		vector<int> CPZoneNums = Fields[fNum++].GetReturnIntVec();
		int CPTypeVarNum = Fields[fNum++].GetReturnInt();
		double CellSpacing = Fields[fNum++].GetReturnDouble();
		int RidgeFuncVarNum = 0;
		if (Fields[fNum++].GetReturnBool())
			RidgeFuncVarNum = Fields[fNum++].GetReturnInt();
		else fNum++;
		fNum++;
		vector<bool> DoCalcs;
		for (int i = (int)BondalyzerCalcType_BondPaths; i <= (int)BondalyzerCalcType_RingSurfaces; ++i) DoCalcs.push_back(Fields[fNum++].GetReturnBool());

		int NumZones = TecUtilDataSetGetNumZones();

		if (!UseCPZones && FindCritPoints(VolZoneNum, XYZVarNums, RhoVarNum, GradVarNums, HessVarNums, IsPeriodic, CellSpacing)){
			CPZoneNums.clear();
			for (int i = NumZones + 1; i <= TecUtilDataSetGetNumZones(); ++i) CPZoneNums.push_back(i);
			CPTypeVarNum = VarNumByName(CSMVarName.CritPointType);
		}

		BondalyzerBatch(VolZoneNum, CPZoneNums, CPTypeVarNum, XYZVarNums, RhoVarNum, GradVarNums, HessVarNums, IsPeriodic, RidgeFuncVarNum, DoCalcs);
	}

	TecUtilDialogMessageBox("Finished", MessageBoxType_Information);

	TecUtilLockFinish(AddOnID);
}

void BondalyzerGetUserInfo(BondalyzerCalcType_e CalcType, vector<GuiField_c> const PassthroughFields){
	// new implementation with CSM_Gui

	CurrentCalcType = CalcType;

	vector<GuiField_c> Fields = {
		GuiField_c(Gui_ZoneSelect, "Volume zone", CSMZoneName.FullVolume.substr(0,10)),
		GuiField_c(Gui_Toggle, "Periodic system"),
		GuiField_c(Gui_VertSep)
	};

// 	vector<string> XYZStr = { "X", "Y", "Z" }, HessXYZStr = { "XX", "XY", "XZ", "YY", "YZ", "ZZ" };

	Fields.push_back(GuiField_c(Gui_VarSelect, "X", "X"));
// 	for (int i = 0; i < 3; ++i) Fields.push_back(GuiField_c(Gui_VarSelect, XYZStr[i], XYZStr[i]));
	Fields.push_back(GuiField_c(Gui_VertSep));

	Fields.push_back(GuiField_c(Gui_VarSelect, "Electron Density", CSMVarName.Dens));
	Fields.push_back(GuiField_c(Gui_VertSep));

	int iTmp = Fields.size();
	Fields.push_back(GuiField_c(Gui_ToggleEnable, "Density gradient vector variables present"));
// 	for (int i = 0; i < 3; ++i){
// 		Fields[iTmp].AppendSearchString(to_string(Fields.size()));
// 		Fields.push_back(GuiField_c(Gui_VarSelect, XYZStr[i], CSMVarName.DensGradVec[i]));
// 		if (i < 2) Fields[iTmp].AppendSearchString(",");
// 	}
	Fields[iTmp].AppendSearchString(to_string(Fields.size()));
	Fields.push_back(GuiField_c(Gui_VarSelect, CSMVarName.DensGradVec[0], CSMVarName.DensGradVec[0]));
	Fields.push_back(GuiField_c(Gui_VertSep));

	iTmp = Fields.size();
	Fields.push_back(GuiField_c(Gui_ToggleEnable, "Density Hessian variables present"));
// 	for (int i = 0; i < HessXYZStr.size(); ++i){
// 		Fields[iTmp].AppendSearchString(to_string(Fields.size()));
// 		Fields.push_back(GuiField_c(Gui_VarSelect, HessXYZStr[i], CSMVarName.DensHessTensor[i]));
// 		if (i < 5) Fields[iTmp].AppendSearchString(",");
// 	}
	Fields[iTmp].AppendSearchString(to_string(Fields.size()));
	Fields.push_back(GuiField_c(Gui_VarSelect, CSMVarName.DensHessTensor[0], CSMVarName.DensHessTensor[0]));

	Fields.push_back(GuiField_c(Gui_VertSep));

	string Title = "Find " + BondalyzerStepGUITitles[CalcType];
	if (CalcType == BondalyzerCalcType_GBA) Title = BondalyzerStepGUITitles[CalcType];

	if (CalcType == BondalyzerCalcType_CriticalPoints){
		Fields.push_back(GuiField_c(Gui_Double, "CP search grid spacing", to_string(DefaultCellSpacing)));
	}
	else if (CalcType >= BondalyzerCalcType_BondPaths && CalcType < BondalyzerCalcType_GBA){
		if (CalcType >= BondalyzerCalcType_InteratomicSurfaces){
			Fields.push_back(GuiField_c(Gui_ToggleEnable, "1-RCS function present", to_string(Fields.size() + 1)));
			Fields.push_back(GuiField_c(Gui_VarSelect, "1-ridge function variable", CSMVarName.EberlyFunc[0]));
		}

// 		Fields.push_back(GuiField_c(Gui_ZoneSelect, "Critical points zone", "Critical Points"));
		Fields.push_back(GuiField_c(Gui_VarSelect, "CP type variable", CSMVarName.CritPointType));
		string SearchString = CSMZoneName.CriticalPoints;
// 		if (CalcType == BondalyzerCalcType_BondPaths || CalcType == BondalyzerCalcType_InteratomicSurfaces) SearchString += CSMZoneName.CPType[2];
// 		else if (CalcType == BondalyzerCalcType_RingLines || CalcType == BondalyzerCalcType_RingSurfaces) SearchString += CSMZoneName.CPType[1];


		Fields.push_back(GuiField_c(Gui_ZoneSelectMulti, "Other critical points zones", SearchString));

		int CPTypeNum;
		if (CalcType == BondalyzerCalcType_BondPaths || CalcType == BondalyzerCalcType_InteratomicSurfaces) CPTypeNum = 1;
		else if (CalcType == BondalyzerCalcType_RingLines || CalcType == BondalyzerCalcType_RingSurfaces) CPTypeNum = 2;
		else if (CalcType == BondalyzerCalcType_CageNuclearPaths) CPTypeNum = 3;
		int CPZoneNum = ZoneNumByName(CSMZoneName.CPType[CPTypeNum]);
		if (CPZoneNum > 0) SearchString = CSMZoneName.CPType[CPTypeNum];
		else SearchString = CSMZoneName.CriticalPoints;

		Fields.push_back(GuiField_c(Gui_ZonePointSelectMulti, "Source critical point(s)", SearchString));
	}
	else if (CalcType == BondalyzerCalcType_Batch){
		int TogNum = Fields.size();
		Fields.push_back(GuiField_c(Gui_ToggleEnable, "Use existing critical points zone(s)?"));
		Fields[TogNum].AppendSearchString(to_string(Fields.size()));
		Fields.push_back(GuiField_c(Gui_ZoneSelectMulti, "Critical points", CSMZoneName.CriticalPoints));
		Fields[TogNum].AppendSearchString("," + to_string(Fields.size()));
		Fields.push_back(GuiField_c(Gui_VarSelect, "CP type variable", CSMVarName.CritPointType));
		Fields.push_back(GuiField_c(Gui_Double, "CP search grid spacing", to_string(DefaultCellSpacing)));
		Fields.push_back(GuiField_c(Gui_ToggleEnable, "1-RCS function present", to_string(Fields.size() + 1)));
		Fields.push_back(GuiField_c(Gui_VarSelect, "1-ridge function variable", CSMVarName.EberlyFunc[0]));
		Fields.push_back(GuiField_c(Gui_Label, "Calculation steps"));
		for (int i = (int)BondalyzerCalcType_BondPaths; i <= (int)BondalyzerCalcType_RingSurfaces; ++i) Fields.push_back(GuiField_c(Gui_Toggle, BondalyzerStepGUITitles[i], "1"));
	}

	CSMGui(Title, Fields, BondalyzerReturnUserInfo, AddOnID);
}

int const DeletePointsFromIOrderedZone(int ZoneNum, vector<int> const & DelPointNums){
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

void DeleteCPsReturnUserInfo(bool const GuiSuccess, 
	vector<GuiField_c> const & Fields,
	vector<GuiField_c> const PassthroughFields){
	if (!GuiSuccess) return;

	TecUtilLockStart(AddOnID);

	int fNum = 0;

	vector<int> CPNums = Fields[fNum].GetReturnIntVec();
	int ZoneNum = Fields[fNum++].GetReturnInt(),
		CPTypeVarNum = Fields[fNum++].GetReturnInt();
	bool DeleteFromOtherCPZones = Fields[fNum++].GetReturnBool();
	vector<int> OtherCPZoneNums, UserOtherCPZoneNums = Fields[fNum++].GetReturnIntVec();
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

// 		for (int z = 1; z <= TecUtilDataSetGetNumZones(); ++z)
		for (int const z : UserOtherCPZoneNums)
		if (z != ZoneNum 
			&& TecUtilZoneIsOrdered(z) 
			&& AuxDataZoneItemMatches(z, CSMAuxData.CC.ZoneType, CSMAuxData.CC.ZoneTypeCPs))
		{
			OtherCPZoneNums.push_back(z);
			OtherCPNums.push_back(vector<int>());
			OtherCPZoneTypeRefs.push_back(TecUtilDataValueGetReadableNativeRef(z, CPTypeVarNum));
			for (int i = 0; i < 3; ++i) OtherCPZoneXYZRefs[i].push_back(TecUtilDataValueGetReadableNativeRef(z, XYZVarNums[i]));
		}

		for (int CP : CPNums){
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

	for (int i : NewZoneNums){
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
	for (int z : OtherCPZoneNums) TecUtilSetAddMember(ZoneNumsSet, z, TRUE);

	TecUtilDataSetDeleteZone(ZoneNumsSet);

	TecUtilSetDealloc(&ZoneNumsSet);
	TecUtilLockFinish(AddOnID);
}

void DeleteCPsGetUserInfo(){
	vector<GuiField_c> Fields = {
		GuiField_c(Gui_ZonePointSelectMulti, "Critical point(s)", CSMZoneName.CriticalPoints),
		GuiField_c(Gui_VarSelect, "CP type", CSMVarName.CritPointType),
		GuiField_c(Gui_Toggle, "Also delete from other CP zone(s)", "0"),
		GuiField_c(Gui_ZoneSelectMulti, "Other CP zone(s)", CSMZoneName.CriticalPoints)
	};

	CSMGui("Delete critical point(s)", Fields, DeleteCPsReturnUserInfo, AddOnID);
}

void ExtractCPsReturnUserInfo(bool const GuiSuccess, 
	vector<GuiField_c> const & Fields,
	vector<GuiField_c> const PassthroughFields){
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
	for (auto const & t : CPTypeNums) if (t.size() > 0) CPNumTypes++;

	if (CPNumTypes > 1 && AuxDataZoneItemMatches(ZoneNum, CSMAuxData.CC.ZoneSubType, CSMAuxData.CC.ZoneTypeCPsAll)){
		vector<int> AllCPNums;
		for (auto const & i : CPNums) AllCPNums.insert(AllCPNums.end(), i.cbegin(), i.cend());
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

	for (int i : NewZoneNums) SetCPZone(i);
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

void VarNameFindReplaceReturnUserInfo(bool const GuiSuccess, 
	vector<GuiField_c> const & Fields,
	vector<GuiField_c> const PassthroughFields){
	if (!GuiSuccess) return;

	TecUtilLockStart(AddOnID);

	int FNum = 0;

	string FindStr = Fields[FNum++].GetReturnString(),
		ReplaceStr = Fields[FNum++].GetReturnString();

	vector<int> SearchVars = Fields[FNum++].GetReturnIntVec();

	Set_pa ChangedVars = TecUtilSetAlloc(TRUE);

	for (int i : SearchVars){
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

void ZoneNameFindReplaceReturnUserInfo(bool const GuiSuccess, 
	vector<GuiField_c> const & Fields,
	vector<GuiField_c> const PassthroughFields){
	if (!GuiSuccess) return;

	TecUtilLockStart(AddOnID);

	int FNum = 0;

	string FindStr = Fields[FNum++].GetReturnString(),
		ReplaceStr = Fields[FNum++].GetReturnString();

	vector<int> SearchZones = Fields[FNum++].GetReturnIntVec();

	Set_pa ChangedZones = TecUtilSetAlloc(TRUE);

	for (int i : SearchZones){
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

Boolean_t const GetReadPtrsForZone(int ZoneNum,
	int RhoVarNum,
	vector<int> const & GradVarNums,
	vector<int> const & HessVarNums,
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
	for (auto const & i : GradVarNums) REQUIRE(i > 0 && i <= NumVars);
	REQUIRE(HessVarNums.size() == 6 || HessVarNums.size() == 0);
	for (auto const & i : HessVarNums) REQUIRE(i > 0 && i <= NumVars);

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


Boolean_t FindCritPoints(int VolZoneNum,
								vector<int> const & XYZVarNums,
								int RhoVarNum,
								vector<int> const & GradVarNums,
								vector<int> const & HessVarNums,
								Boolean_t IsPeriodic,
								double const & CellSpacing)
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


string const MakeStringFromCPNums(vector<int> const & CPNums, CritPoints_c const & CPs, vector<vector<int> > & CPTypeOffsets){
	string Out = " (";

	CPTypeOffsets.clear();

	int NumFarField = 0;
	vector<vector<int> > CPsByType(CPNameList.size());

	for (auto const & CP : CPNums){
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
Boolean_t FindBondRingLines(int VolZoneNum,
	vector<int> const & OtherCPZoneNums,
	int SelectedCPZoneNum,
	vector<int> SelectedCPNums,
	int CPTypeVarNum,
	CPType_e const & CPType,
	vector<int> const & XYZVarNums,
	int RhoVarNum,
	vector<int> const & GradVarNums,
	vector<int> const & HessVarNums,
	Boolean_t IsPeriodic)
{
	TecUtilLockStart(AddOnID);

	int NumZones = TecUtilDataSetGetNumZones();
	int NumVars = TecUtilDataSetGetNumVars();

	REQUIRE(XYZVarNums.size() == 3);
	for (auto const & i : XYZVarNums) REQUIRE(i > 0 && i <= NumVars);
	for (auto const & i : OtherCPZoneNums) REQUIRE(i > 0 && i <= NumZones);

	int TypeInd = VectorGetElementNum(CPTypeList, CPType);
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
	string CPZoneCheckString;

	char* tmpName;
	TecUtilZoneGetName(SelectedCPZoneNum, &tmpName);
	CPZoneCheckString = tmpName;
	TecUtilStringDealloc(&tmpName);

	for (int z : OtherCPZoneNums){
		TecUtilZoneGetIJK(z, &NumCPs, &iJunk[0], &iJunk[1]);
		for (auto const & i : iJunk) if (i > 1){
			TecUtilDialogErrMsg("CP zone is not i-ordered");
			return FALSE;
		}
// 		TecUtilZoneGetName(z, &tmpName);
// 		if (CPZoneCheckString == tmpName){
// 			TecUtilDialogErrMsg(string("Duplicate CP zone specified: " + string(tmpName)).c_str());
// 			return FALSE;
// 		}
// 		TecUtilStringDealloc(&tmpName);
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
	MR.BasisVectors = &VolInfo.BasisNormalized;


	CritPoints_c AllCPs(SelectedCPZoneNum, XYZVarNums, CPTypeVarNum, RhoVarNum, &MR);

	for (int z : OtherCPZoneNums){
		AllCPs += CritPoints_c(z, XYZVarNums, CPTypeVarNum, RhoVarNum, &MR);
	}

	int EndCPNumforName = 0;
	StreamDir_e GPDir = StreamDir_Forward;
	ColorIndex_t PathColor = Black_C;
	vector<CPType_e> MinDistTypes = { CPType_Nuclear, CPType_Bond };
	if (CPType == CPType_Ring){
		MinDistTypes = { CPType_Nuclear, CPType_Ring };
		GPDir = StreamDir_Reverse;
		PathColor = Green_C;
		EndCPNumforName++;
	}

	if (SelectedCPNums.size() == 0){
		for (int i = 0; i < AllCPs.NumCPs(TypeInd); ++i) SelectedCPNums.push_back(i+1);
	}

	CSMGuiLock();

	vec3 StartPoint;

	double const StartPointOffset = 0.1 * AllCPs.GetMinCPDist(MinDistTypes);
	int const NumGPPts = 100;
	double TermRadius = 0.2;
	double RhoCutoff = DefaultRhoCutoff;

	vector<GradPath_c> GPs;

	for (int iCP = 0; iCP < SelectedCPNums.size(); ++iCP){
		int CPInd = SelectedCPNums[iCP] - 1;
		for (int d = -1; d < 2; d += 2){
			StartPoint = AllCPs.GetXYZ(TypeInd, CPInd) + AllCPs.GetPrincDir(TypeInd, CPInd) * StartPointOffset * static_cast<double>(d);

			GPs.push_back(GradPath_c(StartPoint, GPDir, NumGPPts, GPType_Classic, GPTerminate_AtCP, NULL, &AllCPs, &TermRadius, &RhoCutoff, VolInfo, HessPtrs, GradPtrs, RhoPtr));
			GPs.back().SetStartEndCPNum(AllCPs.GetTotOffsetFromTypeNumOffset(TypeInd, CPInd), 0);

		}
	}

	int NumGPs = GPs.size();

#ifndef _DEBUG
#pragma omp parallel for schedule(dynamic)
#endif
	for (int iGP = 0; iGP < NumGPs; ++iGP){
		GPs[iGP].Seed(false);
		GPs[iGP].PointPrepend(AllCPs.GetXYZ(GPs[iGP].GetStartEndCPNum()[0]), AllCPs.GetRho(GPs[iGP].GetStartEndCPNum()[0]));
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
		StartEndCPTypeAndOffset[0] = AllCPs.GetTypeNumOffsetFromTotOffset(OldStartCPNum);
		int CPNumForType = StartEndCPTypeAndOffset[0][1];
		string Name = (CPType == CPType_Bond ? CSMZoneName.BondPath : CSMZoneName.RingLine) + CPNameList[StartEndCPTypeAndOffset[0][0]] + " " + to_string(StartEndCPTypeAndOffset[0][1] + 1);
		Name += MakeStringFromCPNums(StartEndCPNums, AllCPs, StartEndCPTypeAndOffset);
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
			if (StartEndCPNums[i] >= 0) 
				AuxDataZoneSetItem(GP.GetZoneNum(), CSMAuxData.CC.GPEndTypes[i], CPNameList[StartEndCPTypeAndOffset[i][0]]);
			else 
				AuxDataZoneSetItem(GP.GetZoneNum(), CSMAuxData.CC.GPEndTypes[i], "FF");
		}
		AuxDataZoneSetItem(GP.GetZoneNum(), CSMAuxData.CC.GPMiddleCPNum, to_string(CPNums[0] + 1));
		AuxDataZoneSetItem(GP.GetZoneNum(), CSMAuxData.CC.GPMiddleCPNumForType, to_string(CPNumForType + 1));
		AuxDataZoneSetItem(GP.GetZoneNum(), CSMAuxData.CC.GPMiddleCPType, (CPType == CPType_Bond ? CPNameList[1] : CPNameList[2]));
		AuxDataZoneSetItem(GP.GetZoneNum(), CSMAuxData.CC.ZoneType, CSMAuxData.CC.ZoneTypeGP);
		if (CPType == CPType_Ring) 
			AuxDataZoneSetItem(GP.GetZoneNum(), CSMAuxData.CC.ZoneSubType, CSMAuxData.CC.ZoneSubTypeRingLine);
		else 
			AuxDataZoneSetItem(GP.GetZoneNum(), CSMAuxData.CC.ZoneSubType, CSMAuxData.CC.ZoneSubTypeBondPath);
	}

	for (auto & GP : GPs){
		vector<int> StartEndCPNums = GP.GetStartEndCPNum();
		vector<vector<int> > StartEndCPTypeAndOffset(2);
		string Name = CSMZoneName.SpecialGradientPath;
		for (int i = 0; i < 2; ++i){
			if (StartEndCPNums[i] >= 0){
				StartEndCPTypeAndOffset[i] = AllCPs.GetTypeNumOffsetFromTotOffset(StartEndCPNums[i]);
				Name += CPNameList[StartEndCPTypeAndOffset[i][0]] + " " + to_string(StartEndCPTypeAndOffset[i][1] + 1);
			}
			else Name += "FF";

			if (i == 0) Name += " to ";
		}

		GP.SaveAsOrderedZone(Name, vector<FieldDataType_e>(), XYZVarNums, RhoVarNum, FALSE, PathColor);
		for (int i = 0; i < 2; ++i){
			AuxDataZoneSetItem(GP.GetZoneNum(), CSMAuxData.CC.GPEndNumStrs[i], to_string(StartEndCPNums[i] + 1));
			if (StartEndCPNums[i] >= 0) 
				AuxDataZoneSetItem(GP.GetZoneNum(), CSMAuxData.CC.GPEndTypes[i], CPNameList[StartEndCPTypeAndOffset[i][0]]);
			else 
				AuxDataZoneSetItem(GP.GetZoneNum(), CSMAuxData.CC.GPEndTypes[i], "FF");
		}
		AuxDataZoneSetItem(GP.GetZoneNum(), CSMAuxData.CC.ZoneType, CSMAuxData.CC.ZoneTypeGP);
		if (CPType == CPType_Ring) 
			AuxDataZoneSetItem(GP.GetZoneNum(), CSMAuxData.CC.ZoneSubType, CSMAuxData.CC.ZoneSubTypeRingLineSegment);
		else 
			AuxDataZoneSetItem(GP.GetZoneNum(), CSMAuxData.CC.ZoneSubType, CSMAuxData.CC.ZoneSubTypeBondPathSegment);
	}

	TecUtilDataLoadEnd();


	CSMGuiUnlock();
	TecUtilLockFinish(AddOnID);

	return TRUE;
}

/*
 *	Create N equidistributed points on a sphere using the algorithm 
 *	described by Deserno at https://www.cmu.edu/biolphys/deserno/pdf/sphere_equi.pdf
 *	Turns out that the number of points on the sphere is only (close to) the
 *	specified number if r == 1.
 */
void EquidistributedSpherePoints(
	int N, 
	double const & r,
	vector<vec3> & PointVec,
	vector<double> & ThetaVec,
	vector<double> & PhiVec)
{
	REQUIRE(N > 10 && r > 0);

	PointVec.clear();
	ThetaVec.clear();
	PhiVec.clear();
	PointVec.reserve(N);
	ThetaVec.reserve(N);
	PhiVec.reserve(N);

// 	double a = 4. * PI * r * r / (double)N;
	double a = 4. * PI / (double)N;
	double d = sqrt(a);
	double M_theta = round(PI / d);
	double d_theta = PI / M_theta;
	double d_phi = a / d_theta;

	for (double m = 0; m < M_theta; ++m){
		double theta = PI * (m + 0.5) / M_theta;
		double M_phi = round(PI2 * sin(theta) / d_phi);
		for (double n = 0; n < M_phi; ++n){
			double phi = PI2 * n / M_phi;
			PointVec.push_back(SphericalToCartesian(r, theta, phi));
			ThetaVec.push_back(theta);
			PhiVec.push_back(phi);
		}
	}
}

double MinFunc_GPLength_Spherical(gsl_vector const * ThetaPhi, void * params){
	MinFuncParams_GPLengthSpherical * GPParams = reinterpret_cast<MinFuncParams_GPLengthSpherical*>(params);

	GPParams->GP = GradPath_c(
		GPParams->StartPointOrigin + SphericalToCartesian(GPParams->r, gsl_vector_get(ThetaPhi, 0), gsl_vector_get(ThetaPhi, 1)),
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

#ifdef _DEBUG
	if (!GPParams->GP.IsMade())
		TecUtilDialogErrMsg("Gradient path in MinFunc_GPLength_Spherical failed");
// 	if (GPParams->EndCPNum >= 0 && GPParams->GP.GetStartEndCPNum(1) != GPParams->EndCPNum)
// 		TecUtilDialogErrMsg("Gradient path in MinFunc_GPLength_Spherical terminated at wrong CP");
#endif

	if (GPParams->Direction != StreamDir_Both)
		GPParams->GP.PointPrepend(GPParams->StartPointOrigin, 0);

	return GPParams->GP.GetLength();
}

/*
*	Create cage-nuclear paths from a source CP zone(s)
*/
Boolean_t FindCageNuclearPaths(int VolZoneNum,
	vector<int> const & OtherCPZoneNums,
	int SelectedCPZoneNum,
	vector<int> SelectedCPNums,
	int CPTypeVarNum,
	vector<int> const & XYZVarNums,
	int RhoVarNum,
	vector<int> const & GradVarNums,
	vector<int> const & HessVarNums,
	Boolean_t IsPeriodic)
{
	TecUtilLockStart(AddOnID);

	int NumZones = TecUtilDataSetGetNumZones();
	int NumVars = TecUtilDataSetGetNumVars();

	REQUIRE(XYZVarNums.size() == 3);
	for (auto const & i : XYZVarNums) REQUIRE(i > 0 && i <= NumVars);
	for (auto const & i : OtherCPZoneNums) REQUIRE(i > 0 && i <= NumZones);

	int TypeInd = 3;

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
	string CPZoneCheckString;

	char* tmpName;
	TecUtilZoneGetName(SelectedCPZoneNum, &tmpName);
	CPZoneCheckString = tmpName;
	TecUtilStringDealloc(&tmpName);

	for (int z : OtherCPZoneNums){
		TecUtilZoneGetIJK(z, &NumCPs, &iJunk[0], &iJunk[1]);
		for (auto const & i : iJunk) if (i > 1){
			TecUtilDialogErrMsg("CP zone is not i-ordered");
			return FALSE;
		}
// 		TecUtilZoneGetName(z, &tmpName);
// 		if (CPZoneCheckString == tmpName){
// 			TecUtilDialogErrMsg(string("Duplicate CP zone specified: " + string(tmpName)).c_str());
// 			return FALSE;
// 		}
// 		TecUtilStringDealloc(&tmpName);
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
	MR.BasisVectors = &VolInfo.BasisNormalized;


	CritPoints_c AllCPs(SelectedCPZoneNum, XYZVarNums, CPTypeVarNum, RhoVarNum, &MR);

	for (int z : OtherCPZoneNums){
		AllCPs += CritPoints_c(z, XYZVarNums, CPTypeVarNum, RhoVarNum, &MR);
	}

	int EndCPNumforName = 0;
	StreamDir_e GPDir = StreamDir_Both;
	ColorIndex_t PathColor = Red_C;
	vector<CPType_e> MinDistTypes = { CPType_Ring, CPType_Nuclear }; //ring CPs are the closest to cage CPs, and don't want to seed past a neighboring ring CP, so use that distance.

	if (SelectedCPNums.size() == 0){
		for (int i = 0; i < AllCPs.NumCPs(TypeInd); ++i) SelectedCPNums.push_back(i + 1);
	}

	CSMGuiLock();

	vec3 StartPoint;

	int const NumSpherePoints = 5000;
	int const NumGPPts = 100;
	double TermRadius = 0.2;
	double RhoCutoff = DefaultRhoCutoff;

	for (int iCP = 0; iCP < SelectedCPNums.size(); ++iCP){
		double StartPointOffset = 0.3 * AllCPs.GetMinCPDist(TypeInd, SelectedCPNums[iCP] - 1, MinDistTypes);

		vector<GradPath_c> GPs;
		vector<vec3> SeedPoints;
		vector<double> SeedTheta, SeedPhi;

		EquidistributedSpherePoints(NumSpherePoints, StartPointOffset, SeedPoints, SeedTheta, SeedPhi);

		for (vec3 & Pt : SeedPoints){
			Pt += AllCPs.GetXYZ(TypeInd, SelectedCPNums[iCP] - 1);
			GPs.push_back(GradPath_c(Pt, GPDir, NumGPPts, GPType_Classic, GPTerminate_AtCP, NULL, &AllCPs, &TermRadius, &RhoCutoff, VolInfo, HessPtrs, GradPtrs, RhoPtr));
			GPs.back().SetStartEndCPNum(AllCPs.GetTotOffsetFromTypeNumOffset(TypeInd, SelectedCPNums[iCP] - 1), 0);
		}

		int NumGPs = GPs.size();
		vector<double> GPLengths(NumGPs);
		vector<CPType_e> GPEndCPTypes(NumGPs);
		vector<int> GPEndCPNums(NumGPs);

#ifndef _DEBUG
#pragma omp parallel for schedule(dynamic)
#endif
		for (int iGP = 0; iGP < NumGPs; ++iGP){
			GPs[iGP].Seed(false);
// 			GPs[iGP].PointPrepend(AllCPs.GetXYZ(GPs[iGP].GetStartEndCPNum()[0]), AllCPs.GetRho(GPs[iGP].GetStartEndCPNum()[0]));
			GPLengths[iGP] = GPs[iGP].GetLength();
			GPEndCPNums[iGP] = GPs[iGP].GetStartEndCPNum(1);
			GPEndCPTypes[iGP] = AllCPs.GetTypeFromTotOffset(GPEndCPNums[iGP]);
		}

		vector<bool> GPUsed(GPs.size(), false);
		vector<int> EndCPNums;
		vector<double> StartTheta, StartPhi;

		for (int iGP = 0; iGP < NumGPs; ++iGP){
			if (!GPUsed[iGP] && GPEndCPTypes[iGP] == CPType_Nuclear){
				int MinInd = iGP;
				for (int jGP = iGP + 1; jGP < iGP + NumGPs - 1; ++jGP){
					int ind = jGP % NumGPs;
					if (!GPUsed[ind] && GPEndCPNums[ind] == GPEndCPNums[iGP]){
						if (GPLengths[ind] < GPLengths[MinInd]) MinInd = ind;
						GPUsed[ind] = true;
					}
				}
				EndCPNums.push_back(GPEndCPNums[MinInd]);
				StartTheta.push_back(SeedTheta[MinInd]);
				StartPhi.push_back(SeedPhi[MinInd]);
			}
			GPUsed[iGP] = true;
		}

		int NumSGPs = EndCPNums.size();
		vector<GradPath_c> SGPs(NumSGPs);

#ifndef _DEBUG
#pragma omp parallel for schedule(dynamic)
#endif
		for (int iGP = 0; iGP < NumSGPs; ++iGP){
			MinFuncParams_GPLengthSpherical GPParams;
			GPParams.Direction = GPDir;
			GPParams.NumGPPoints = NumGPPts;
			GPParams.GPType = GPType_Classic;
			GPParams.GPTerminate = GPTerminate_AtCP;
			GPParams.TermPoint = NULL;
			GPParams.CPs = &AllCPs;
			GPParams.TermPointRadius = &TermRadius;
			GPParams.TermValue = &RhoCutoff;
			GPParams.VolInfo = &VolInfo;
			GPParams.HessPtrs = &HessPtrs;
			GPParams.GradPtrs = &GradPtrs;
			GPParams.RhoPtr = &RhoPtr;

			GPParams.StartPointOrigin = AllCPs.GetXYZ(TypeInd, SelectedCPNums[iCP] - 1);
			GPParams.r = StartPointOffset;
			GPParams.StartCPNum = AllCPs.GetTotOffsetFromTypeNumOffset(TypeInd, SelectedCPNums[iCP] - 1);
			GPParams.EndCPNum = EndCPNums[iGP];

			gsl_multimin_fminimizer_type const *T = gsl_multimin_fminimizer_nmsimplex2;
			gsl_multimin_fminimizer *s = NULL;
			gsl_vector *ss, *x;
			gsl_multimin_function minex_func;

			size_t iter = 0;
			int status;
			double size;

			/* Starting point */
			x = gsl_vector_alloc(2);
			gsl_vector_set(x, 0, StartTheta[iGP]);
			gsl_vector_set(x, 1, StartPhi[iGP]);

			/* Set initial step sizes to 0.5 degree (0.0087 radians) */
			ss = gsl_vector_alloc(2);
			gsl_vector_set_all(ss, 0.0087);

			/* Initialize method and iterate */
			minex_func.n = 2;
			minex_func.f = MinFunc_GPLength_Spherical;
			minex_func.params = reinterpret_cast<void*>(&GPParams);

			s = gsl_multimin_fminimizer_alloc(T, 2);
			gsl_multimin_fminimizer_set(s, &minex_func, x, ss);

			do
			{
				iter++;
				status = gsl_multimin_fminimizer_iterate(s);

				if (status)
					break;

				size = gsl_multimin_fminimizer_size(s);
				status = gsl_multimin_test_size(size, 1e-6);
			} while (status == GSL_CONTINUE && iter < 100);

			if (status == GSL_SUCCESS){
				SGPs[iGP] = reinterpret_cast<MinFuncParams_GPLengthSpherical*>(minex_func.params)->GP;
				SGPs[iGP].Resample(NumGPPts);
			}

			gsl_vector_free(x);
			gsl_vector_free(ss);
			gsl_multimin_fminimizer_free(s);
		}

		for (auto & GP : SGPs){
			vector<int> StartEndCPNums = GP.GetStartEndCPNum();
			vector<vector<int> > StartEndCPTypeAndOffset(2);
			string Name = CSMZoneName.SpecialGradientPath;
			for (int i = 0; i < 2; ++i){
				if (StartEndCPNums[i] >= 0){
					StartEndCPTypeAndOffset[i] = AllCPs.GetTypeNumOffsetFromTotOffset(StartEndCPNums[i]);
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
			AuxDataZoneSetItem(GP.GetZoneNum(), CSMAuxData.CC.ZoneSubType, CSMAuxData.CC.ZoneSubTypeCageNuclearPath);
		}
	}

// 	int NumGPs = GPs.size();
// 
// #ifndef _DEBUG
// #pragma omp parallel for schedule(dynamic)
// #endif
// 	for (int iGP = 0; iGP < NumGPs; ++iGP){
// 		GPs[iGP].Seed(false);
// 		GPs[iGP].PointPrepend(AllCPs.GetXYZ(GPs[iGP].GetStartEndCPNum()[0]), AllCPs.GetRho(GPs[iGP].GetStartEndCPNum()[0]));
// 		GPs[iGP].Resample(NumGPPts);
// 		if (GPDir == StreamDir_Reverse) GPs[iGP].Reverse();
// 	}
// 
// 	for (int iGP = 0; iGP < NumGPs; iGP += 2){
// 		GradPath_c GP = GPs[iGP];
// 		vector<int> CPNums;
// 		CPNums.push_back(GP.GetStartEndCPNum()[EndCPNumforName]);
// 		int OldStartCPNum = GP.GetStartEndCPNum()[EndCPNumforName];
// 		GP += GPs[iGP + 1];
// 		vector<int> StartEndCPNums = GP.GetStartEndCPNum();
// 		vector<vector<int> > StartEndCPTypeAndOffset(2);
// 		StartEndCPTypeAndOffset[0] = AllCPs.GetTypeNumOffsetFromTotOffset(OldStartCPNum);
// 		string Name = (CPType == CPType_BondCP ? CSMZoneName.BondPath : CSMZoneName.RingLine) + CPNameList[StartEndCPTypeAndOffset[0][0]] + " " + to_string(StartEndCPTypeAndOffset[0][1] + 1);
// 		Name += MakeStringFromCPNums(StartEndCPNums, AllCPs, StartEndCPTypeAndOffset);
// 		// 		for (int i = 0; i < 2; ++i){
// 		// 			if (StartEndCPNums[i] >= 0){
// 		// 				StartEndCPTypeAndOffset[i] = CPs.GetTypeNumOffsetFromTotOffset(StartEndCPNums[i]);
// 		// 				Name += CPNameList[StartEndCPTypeAndOffset[i][0]][0] + to_string(StartEndCPTypeAndOffset[i][1] + 1);
// 		// 			}
// 		// 			else Name += "FF";
// 		// 
// 		// 			if (i == 0) Name += "-";
// 		// 		}
// 		// 
// 		// 		Name += ")";
// 
// 		if (StartEndCPNums[0] < 0){
// 			int TmpInt = StartEndCPNums[0];
// 			StartEndCPNums[0] = StartEndCPNums[1];
// 			StartEndCPNums[1] = TmpInt;
// 		}
// 		GP.SaveAsOrderedZone(Name, vector<FieldDataType_e>(), XYZVarNums, RhoVarNum, TRUE, PathColor);
// 		for (int i = 0; i < 2; ++i){
// 			AuxDataZoneSetItem(GP.GetZoneNum(), CSMAuxData.CC.GPEndNumStrs[i], to_string(StartEndCPNums[i] + 1));
// 			if (StartEndCPNums[i] >= 0) AuxDataZoneSetItem(GP.GetZoneNum(), CSMAuxData.CC.GPEndTypes[i], CPNameList[StartEndCPTypeAndOffset[i][0]]);
// 			else AuxDataZoneSetItem(GP.GetZoneNum(), CSMAuxData.CC.GPEndTypes[i], "FF");
// 		}
// 		AuxDataZoneSetItem(GP.GetZoneNum(), CSMAuxData.CC.ZoneType, CSMAuxData.CC.ZoneTypeGP);
// 		if (CPType == CPType_RingCP) AuxDataZoneSetItem(GP.GetZoneNum(), CSMAuxData.CC.ZoneSubType, CSMAuxData.CC.ZoneSubTypeRingLine);
// 		else AuxDataZoneSetItem(GP.GetZoneNum(), CSMAuxData.CC.ZoneSubType, CSMAuxData.CC.ZoneSubTypeBondPath);
// 	}
// 
	

	TecUtilDataLoadEnd();


	CSMGuiUnlock();
	TecUtilLockFinish(AddOnID);

	return TRUE;
}

/*
*	Shotgun GPs around selected nuclear/cage CPs and save them as ordered zones.
*/
Boolean_t const GradientPathsOnSphere(int NumSphereGPs,
	int VolZoneNum,
	vector<int> const & OtherCPZoneNums,
	int SelectedCPZoneNum,
	vector<int> const & SelectedCPNums,
	int CPTypeVarNum,
	vector<int> const & XYZVarNums,
	int RhoVarNum,
	vector<int> const & GradVarNums,
	vector<int> const & HessVarNums,
	Boolean_t IsPeriodic)
{
	TecUtilLockStart(AddOnID);

	int NumZones = TecUtilDataSetGetNumZones();
	int NumVars = TecUtilDataSetGetNumVars();

	REQUIRE(XYZVarNums.size() == 3);
	for (auto const & i : XYZVarNums) REQUIRE(i > 0 && i <= NumVars);
	for (auto const & i : OtherCPZoneNums) REQUIRE(i > 0 && i <= NumZones);

	char *ZoneName;
	TecUtilZoneGetName(SelectedCPZoneNum, &ZoneName);
	int TypeInd = VectorGetElementNum(CSMZoneName.CPType, string(ZoneName));
	if ((TypeInd >= 0 && TypeInd != 0 && TypeInd != 3) || (TypeInd < 0 && string(ZoneName) != CSMZoneName.CriticalPoints)){
		TecUtilDialogErrMsg("Please select nuclear, cage, or \"Critical Points\" zone");
		return FALSE;
	}
	TecUtilStringDealloc(&ZoneName);

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
	string CPZoneCheckString;

	char* tmpName;
	TecUtilZoneGetName(SelectedCPZoneNum, &tmpName);
	CPZoneCheckString = tmpName;
	TecUtilStringDealloc(&tmpName);

	for (int z : OtherCPZoneNums){
		TecUtilZoneGetIJK(z, &NumCPs, &iJunk[0], &iJunk[1]);
		for (auto const & i : iJunk) if (i > 1){
			TecUtilDialogErrMsg("CP zone is not i-ordered");
			return FALSE;
		}
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
	MR.BasisVectors = &VolInfo.BasisNormalized;


	CritPoints_c AllCPs(SelectedCPZoneNum, XYZVarNums, CPTypeVarNum, RhoVarNum, &MR);

	for (int z : OtherCPZoneNums){
		AllCPs += CritPoints_c(z, XYZVarNums, CPTypeVarNum, RhoVarNum, &MR);
	}

	int EndCPNumforName = 0;
	StreamDir_e GPDir = StreamDir_Both;
	ColorIndex_t PathColor = Red_C;
	
	CSMGuiLock();

	vec3 StartPoint;

	int const NumSpherePoints = NumSphereGPs;
	int const NumGPPts = 100;
	double TermRadius = 0.2;
	double RhoCutoff = DefaultRhoCutoff;

	Set NewZones;

	for (int iCP = 0; iCP < SelectedCPNums.size(); ++iCP){
		vector<int> CPTypeIndOffset;
		if (TypeInd < 0) CPTypeIndOffset = AllCPs.GetTypeNumOffsetFromTotOffset(SelectedCPNums[iCP] - 1);
		else CPTypeIndOffset = { TypeInd, SelectedCPNums[iCP] - 1 };
		if (CPTypeIndOffset[0] != 0 && CPTypeIndOffset[0] != 3) continue;


		double StartPointOffset = 0.3 * AllCPs.GetMinCPDist(CPTypeIndOffset[0], CPTypeIndOffset[1]);

		vector<GradPath_c> GPs;
		vector<vec3> SeedPoints;
		vector<double> SeedTheta, SeedPhi;

		EquidistributedSpherePoints(NumSpherePoints, StartPointOffset, SeedPoints, SeedTheta, SeedPhi);

		for (vec3 & Pt : SeedPoints){
			Pt += AllCPs.GetXYZ(CPTypeIndOffset[0], CPTypeIndOffset[1]);
			GPs.push_back(GradPath_c(Pt, GPDir, NumGPPts, GPType_Classic, GPTerminate_AtCP, NULL, &AllCPs, &TermRadius, &RhoCutoff, VolInfo, HessPtrs, GradPtrs, RhoPtr));
			GPs.back().SetStartEndCPNum(AllCPs.GetTotOffsetFromTypeNumOffset(CPTypeIndOffset[0], CPTypeIndOffset[1]), (CPTypeIndOffset[0] == 0 ? 1 : 0));
		}

		int NumGPs = GPs.size();
		vector<int> GPEndCPNums(NumGPs);

#ifndef _DEBUG
#pragma omp parallel for schedule(dynamic)
#endif
		for (int iGP = 0; iGP < NumGPs; ++iGP){
			GPs[iGP].Seed(false);
			GPEndCPNums[iGP] = (CPTypeIndOffset[0] == 0 ? GPs[iGP].GetStartEndCPNum(0) : GPs[iGP].GetStartEndCPNum(1));
		}

		vector<bool> GPUsed(GPs.size(), false);
		vector<int> EndCPNums;

		for (int iGP = 0; iGP < NumGPs; ++iGP){
			if (!GPUsed[iGP]){
				EndCPNums.push_back(GPEndCPNums[iGP]);
				for (int jGP = iGP + 1; jGP < iGP + NumGPs - 1; ++jGP){
					int ind = jGP % NumGPs;
					if (!GPUsed[ind] && GPEndCPNums[ind] == GPEndCPNums[iGP]){
						GPUsed[ind] = true;
					}
				}
			}
			GPUsed[iGP] = true;
		}

		std::sort(EndCPNums.begin(), EndCPNums.end());
		if (EndCPNums.size() > 1 && EndCPNums[0] < 0){
			vector<int> NewEndCPNums(EndCPNums.begin() + 1, EndCPNums.end());
			NewEndCPNums.push_back(EndCPNums[0]);
			EndCPNums = NewEndCPNums;
		}

		int iGP = 0;
		for (int jCP = 0; jCP < EndCPNums.size(); ++jCP){
			int ColorInd = jCP;
			if ((ColorIndex_t)jCP >= White_C) ColorInd++;
			ColorInd = ColorInd % 64;
			for (int iGP = 0; iGP < GPs.size(); ++iGP){
				if (GPEndCPNums[iGP] == EndCPNums[jCP]){
					vector<int> StartEndCPNums = GPs[iGP].GetStartEndCPNum();
					vector<vector<int> > StartEndCPTypeAndOffset(2);
					string Name = CSMZoneName.SpecialGradientPath;
					for (int i = 0; i < 2; ++i){
						if (StartEndCPNums[i] >= 0){
							StartEndCPTypeAndOffset[i] = AllCPs.GetTypeNumOffsetFromTotOffset(StartEndCPNums[i]);
							Name += CPNameList[StartEndCPTypeAndOffset[i][0]] + " " + to_string(StartEndCPTypeAndOffset[i][1] + 1);
						}
						else Name += "FF";

						if (i == 0) Name += " to ";
					}

					NewZones += GPs[iGP].SaveAsOrderedZone(Name, vector<FieldDataType_e>(), XYZVarNums, RhoVarNum, FALSE, (ColorIndex_t)ColorInd);
					for (int i = 0; i < 2; ++i){
						AuxDataZoneSetItem(GPs[iGP].GetZoneNum(), CSMAuxData.CC.GPEndNumStrs[i], to_string(StartEndCPNums[i] + 1));
						if (StartEndCPNums[i] >= 0) AuxDataZoneSetItem(GPs[iGP].GetZoneNum(), CSMAuxData.CC.GPEndTypes[i], CPNameList[StartEndCPTypeAndOffset[i][0]]);
						else AuxDataZoneSetItem(GPs[iGP].GetZoneNum(), CSMAuxData.CC.GPEndTypes[i], "FF");
					}
					AuxDataZoneSetItem(GPs[iGP].GetZoneNum(), CSMAuxData.CC.ZoneType, CSMAuxData.CC.ZoneTypeGP);
					AuxDataZoneSetItem(GPs[iGP].GetZoneNum(), CSMAuxData.CC.ZoneSubType, CSMAuxData.CC.ZoneSubTypeCageNuclearPath);
				}
			}
		}
	}

	TecUtilZoneSetActive(NewZones.getRef(), AssignOp_PlusEquals);

	TecUtilDataLoadEnd();


	CSMGuiUnlock();
	TecUtilLockFinish(AddOnID);

	return TRUE;
}

void GradientPathsOnSphereReturnUserInfo(bool const GuiSuccess,
	vector<GuiField_c> const & Fields,
	vector<GuiField_c> const PassthroughFields){
	if (!GuiSuccess) return;

	TecUtilLockStart(AddOnID);

	int VolZoneNum, RhoVarNum;
	vector<int> XYZVarNums(3), GradVarNums, HessVarNums;
	Boolean_t IsPeriodic;

	int fNum = 0;

	VolZoneNum = Fields[fNum++].GetReturnInt();
	IsPeriodic = Fields[fNum++].GetReturnBool();
	fNum++;
	for (int i = 0; i < 3; ++i) XYZVarNums[i] = i + Fields[fNum].GetReturnInt();
	fNum++;
	fNum++;
	RhoVarNum = Fields[fNum++].GetReturnInt();
	fNum++;


	if (Fields[fNum++].GetReturnBool()){
		GradVarNums.resize(3);
		for (int i = 0; i < 3; ++i) GradVarNums[i] = i + Fields[fNum].GetReturnInt();
	}
	fNum++;
	fNum++;

	if (Fields[fNum++].GetReturnBool()){
		HessVarNums.resize(6);
		for (int i = 0; i < 6; ++i) HessVarNums[i] = i + Fields[fNum].GetReturnInt();
	}
	fNum += 2;

	int CPTypeVarNum = Fields[fNum++].GetReturnInt();
	vector<int> OtherCPZoneNums = Fields[fNum++].GetReturnIntVec();
	int SelectCPsZoneNum = Fields[fNum].GetReturnInt();
	vector<int> SelectedCPs = Fields[fNum++].GetReturnIntVec();
	int NumGPs = Fields[fNum++].GetReturnInt();

	GradientPathsOnSphere(NumGPs, VolZoneNum, OtherCPZoneNums, SelectCPsZoneNum, SelectedCPs, CPTypeVarNum, XYZVarNums, RhoVarNum, GradVarNums, HessVarNums, IsPeriodic);

	TecUtilDialogMessageBox("Finished", MessageBoxType_Information);

	TecUtilLockFinish(AddOnID);
}

void GradientPathsOnSphereGetUserInfo(){

	vector<GuiField_c> Fields = {
		GuiField_c(Gui_ZoneSelect, "Volume zone", CSMZoneName.FullVolume.substr(0, 10)),
		GuiField_c(Gui_Toggle, "Periodic system"),
		GuiField_c(Gui_VertSep)
	};

	Fields.push_back(GuiField_c(Gui_VarSelect, "X", "X"));

	Fields.push_back(GuiField_c(Gui_VertSep));

	Fields.push_back(GuiField_c(Gui_VarSelect, "Electron Density", CSMVarName.Dens));
	Fields.push_back(GuiField_c(Gui_VertSep));

	int iTmp = Fields.size();
	Fields.push_back(GuiField_c(Gui_ToggleEnable, "Density gradient vector variables present"));

	Fields[iTmp].AppendSearchString(to_string(Fields.size()));
	Fields.push_back(GuiField_c(Gui_VarSelect, CSMVarName.DensGradVec[0], CSMVarName.DensGradVec[0]));
	Fields.push_back(GuiField_c(Gui_VertSep));

	iTmp = Fields.size();
	Fields.push_back(GuiField_c(Gui_ToggleEnable, "Density Hessian variables present"));

	Fields[iTmp].AppendSearchString(to_string(Fields.size()));
	Fields.push_back(GuiField_c(Gui_VarSelect, CSMVarName.DensHessTensor[0], CSMVarName.DensHessTensor[0]));

	Fields.push_back(GuiField_c(Gui_VertSep));

	Fields.push_back(GuiField_c(Gui_VarSelect, "CP type variable", CSMVarName.CritPointType));
	string SearchString = CSMZoneName.CriticalPoints;

	Fields.push_back(GuiField_c(Gui_ZoneSelectMulti, "Other critical points zones", SearchString));

	Fields.push_back(GuiField_c(Gui_ZonePointSelectMulti, "Source critical point(s)", SearchString));
 
	Fields.push_back(GuiField_c(Gui_Int, "Number of gradient paths", "2048"));

	CSMGui("GPs around cage/nuclear CPs", Fields, GradientPathsOnSphereReturnUserInfo, AddOnID);
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

	if (GPParams->Direction != StreamDir_Both)
		GPParams->GP.SetStartEndCPNum(GPParams->StartCPNum, 0);

	GPParams->GP.Seed(false);

	if (!GPParams->GP.IsMade())
		TecUtilDialogErrMsg("Gradient path in MinFunc_GPLength_InPlane failed");
// 	if (GPParams->EndCPNum >= 0 && GPParams->EndCPPosition >= 0 && GPParams->GP.GetStartEndCPNum(GPParams->EndCPPosition) != GPParams->EndCPNum)
// 		TecUtilDialogErrMsg("Gradient path in MinFunc_GPLength_InPlane terminated at wrong CP");
// 	if (GPParams->StartCPNum >= 0 && GPParams->EndCPPosition >= 0 && GPParams->GP.GetStartEndCPNum((GPParams->EndCPPosition + 1) % 2) != GPParams->StartCPNum)
// 		TecUtilDialogErrMsg("Gradient path in MinFunc_GPLength_InPlane originated at wrong CP");

	if (GPParams->Direction != StreamDir_Both)
		GPParams->GP.PointPrepend(GPParams->StartPointOrigin, 0);

#ifdef _DEBUG
	MinFunc_GPAngleLengths.push_back(vector<double>({ alpha, GPParams->GP.GetLength() }));
#endif

	return GPParams->GP.GetLength();
}

template <typename T>
T const PeriodicDistance1D(T a, T b, T const & low, T const & high){
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
	if (dist < 0) dist *= -1.;

	return dist;
}

/*
*	Create interatomic (bond) and ring surfaces from a source CP zone(s)
*	Also creates bond-ring, and bond-cage (ring-nuclear) paths.
*/
Boolean_t FindBondRingSurfaces(int VolZoneNum,
	vector<int> const & OtherCPZoneNums,
	int SelectedCPZoneNum,
	vector<int> SelectedCPNums,
	int CPTypeVarNum,
	CPType_e const & CPType,
	vector<int> const & XYZVarNums,
	int RhoVarNum,
	vector<int> const & GradVarNums,
	vector<int> const & HessVarNums, 
	int RCSFuncVarNum,
	Boolean_t IsPeriodic)
{
	TecUtilLockStart(AddOnID);

	int NumZones = TecUtilDataSetGetNumZones();
	int NumVars = TecUtilDataSetGetNumVars();

	REQUIRE(XYZVarNums.size() == 3);
	for (auto const & i : XYZVarNums) REQUIRE(i > 0 && i <= NumVars);
	for (auto const & i : OtherCPZoneNums) REQUIRE(i > 0 && i <= NumZones);

	int TypeInd = VectorGetElementNum(CPTypeList, CPType);
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
	string CPZoneCheckString;

	char* tmpName;
	TecUtilZoneGetName(SelectedCPZoneNum, &tmpName);
	CPZoneCheckString = tmpName;
	TecUtilStringDealloc(&tmpName);

	for (int z : OtherCPZoneNums){
		TecUtilZoneGetIJK(z, &NumCPs, &iJunk[0], &iJunk[1]);
		for (auto const & i : iJunk) if (i > 1){
			TecUtilDialogErrMsg("CP zone is not i-ordered");
			return FALSE;
		}
// 		TecUtilZoneGetName(z, &tmpName);
// 		if (CPZoneCheckString == tmpName){
// 			TecUtilDialogErrMsg(string("Duplicate CP zone specified: " + string(tmpName)).c_str());
// 			return FALSE;
// 		}
// 		TecUtilStringDealloc(&tmpName);
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
	MR.BasisVectors = &VolInfo.BasisNormalized;


	CritPoints_c AllCPs(SelectedCPZoneNum, XYZVarNums, CPTypeVarNum, RhoVarNum, &MR);

	for (int z : OtherCPZoneNums){
		AllCPs += CritPoints_c(z, XYZVarNums, CPTypeVarNum, RhoVarNum, &MR);
	}

	AllCPs.RemoveSpuriousCPs();

	if (SelectedCPNums.size() == 0){
		for (int i = 0; i < AllCPs.NumCPs(TypeInd); ++i) SelectedCPNums.push_back(i + 1);
	}

	CSMGuiLock();

	vector<CPType_e> CPTypesForMinDistCheck;
	int EndCPPosition = 1;
	StreamDir_e GPDir = StreamDir_Reverse;
	ColorIndex_t PathColor = Black_C;
	CPType_e PrimaryTermType = CPType_Cage,
		SecondaryTermType = CPType_Ring;
	if (CPType == CPType_Ring){
		GPDir = StreamDir_Forward;
		PathColor = Green_C;
// 		EndCPPosition++;
		PrimaryTermType = CPType_Nuclear;
		SecondaryTermType = CPType_Bond;
		CPTypesForMinDistCheck = { CPType_Bond, CPType_Ring };
	}
	else{
		if (AllCPs.NumBonds() <= 1) CPTypesForMinDistCheck.push_back(CPType_Nuclear);
		
		CPTypesForMinDistCheck.push_back(CPType_Bond);

		if (AllCPs.NumRings() > 0) CPTypesForMinDistCheck.push_back(CPType_Ring);
	}



	vec3 StartPoint;

	int const NumGPPts = 100;
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

	string StatusStr = "Finding " + string(CPType == CPType_Bond ? "interatomic" : "ring") + " surfaces";
	StatusLaunch(StatusStr, AddOnID, TRUE);

#ifndef _DEBUG
// #pragma omp parallel for schedule(dynamic)
#endif
	for (int iCP = 0; iCP < NumCPs; ++iCP){
		if (!StatusUpdate(iCP, NumCPs, StatusStr + ": (" + to_string(iCP+1) + " of " + to_string(NumCPs) + ")", AddOnID)) break;
		int cpNum = SelectedCPNums[iCP] - 1;
		MinFuncParams_GPLengthInPlane GPParams;
		GPParams.Direction = GPDir;
		GPParams.NumGPPoints = NumGPPts;
		GPParams.GPType = GPType_Classic;
		GPParams.GPTerminate = GPTerminate_AtCP;
		GPParams.TermPoint = NULL;
		GPParams.CPs = &AllCPs;
		GPParams.TermPointRadius = &TermRadius;
		GPParams.TermValue = &RhoCutoff;
		GPParams.VolInfo = &VolInfo;
		GPParams.HessPtrs = &HessPtrs;
		GPParams.GradPtrs = &GradPtrs;
		GPParams.RhoPtr = &RhoPtr;

		StartPointOffset = 0.1 * AllCPs.GetMinCPDist(TypeInd, cpNum);
		vec3 StartVec = normalise(AllCPs.GetEigVecs(TypeInd, cpNum).col(1)) * StartPointOffset;

		GPParams.StartPointOrigin = AllCPs.GetXYZ(TypeInd, cpNum);
		GPParams.RotVec = StartVec;
		GPParams.RotAxis = AllCPs.GetPrincDir(TypeInd, cpNum);

		GPParams.StartCPNum = AllCPs.GetTotOffsetFromTypeNumOffset(TypeInd, cpNum);

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

			GPs.push_back(GradPath_c(StartPoint, GPDir, NumGPPts, GPType_Classic, GPTerminate_AtCP, NULL, &AllCPs, &TermRadius, &RhoCutoff, VolInfo, HessPtrs, GradPtrs, RhoPtr));
			GPs.back().SetStartEndCPNum(GPParams.StartCPNum, 0);
		}

#ifndef _DEBUG
#pragma omp parallel for schedule(dynamic)
#endif
		for (int gpNum = 0; gpNum < NumCircleCheckGPs; ++gpNum){
			GPs[gpNum].Seed(false);
			if (GPDir != StreamDir_Both) GPs[gpNum].PointPrepend(GPParams.StartPointOrigin, AllCPs.GetRho(TypeInd, cpNum));
			GPLengths[gpNum] = GPs[gpNum].GetLength();
			GPEndCPNums[gpNum] = GPs[gpNum].GetStartEndCPNum(EndCPPosition);
			GPEndCPTypes[gpNum] = AllCPs.GetTypeFromTotOffset(GPEndCPNums[gpNum]);
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
						TermCPNum = GPEndCPNums[gpNum];
// 						TermCPNum = GPs[gpNum].GetStartEndCPNum(1);
						MinGPNum = gpNum;
						for (int i = gpNum; i < gpNum + NumCircleCheckGPs; ++i){
							int ii = i % NumCircleCheckGPs;
// 							if (!GPUsed[ii] && GPs[ii].GetStartEndCPNum(1) == TermCPNum){
							if (!GPUsed[ii] && GPEndCPNums[ii] == TermCPNum){
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

						AlphaLower = static_cast<double>(gpNum) * AngleStep;
						AlphaUpper = static_cast<double>(gpNum + 1) * AngleStep;

						int Iter = 0;
						while (TermCPNum < 0 && Iter < 15){
							Iter++;
							AlphaGuess = (AlphaLower + AlphaUpper) * 0.5;
							StartPoint = GPParams.StartPointOrigin + Rotate(StartVec, AlphaGuess, GPParams.RotAxis);
							GradPath_c GP(StartPoint, GPDir, NumGPPts, GPType_Classic, GPTerminate_AtCP, NULL, &AllCPs, &TermRadius, &RhoCutoff, VolInfo, HessPtrs, GradPtrs, RhoPtr);
							GP.SetStartEndCPNum(GPParams.StartCPNum, 0);
							GP.Seed(false);
							if (AllCPs.GetTypeFromTotOffset(GP.GetStartEndCPNum(EndCPPosition)) == SecondaryTermType){
								TermCPNum = GP.GetStartEndCPNum(EndCPPosition);
								if (GPDir != StreamDir_Both) GP.PointPrepend(GPParams.StartPointOrigin, AllCPs.GetRho(TypeInd, cpNum));
								GPLen = GP.GetLength();

// 								AlphaLower = (static_cast<double>(gpNum) * AngleStep + AlphaGuess) * 0.5;
// 								AlphaUpper = (static_cast<double>(gpNum + 1) * AngleStep + AlphaGuess) * 0.5;
							}
							else{
								bool GoLeft;
								if (GPEndCPTypes[gpNum] == PrimaryTermType){
									if (GP.GetStartEndCPNum(EndCPPosition) == CPNumLower) GoLeft = false;
									else if (GP.GetStartEndCPNum(EndCPPosition) == CPNumUpper) GoLeft = true;
									else{
// 										TecUtilDialogErrMsg(string("GP terminated at unexpected CP when performing binary search.\n\ngpNum = " + to_string(gpNum)).c_str());
										break;
									}
								}
								else{
									vec3 TmpTermVec = GP[-1] - GP[-2];
									if (dot(TmpTermVec, TermVecLower) > 0.9) GoLeft = false;
									else if (dot(TmpTermVec, TermVecUpper) > 0.9) GoLeft = true;
									else{
// 										TecUtilDialogErrMsg("GP terminated not parallel to either bounding GPs when performing binary search.");
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
								GradPath_c GP(StartPoint, GPDir, NumGPPts, GPType_Classic, GPTerminate_AtCP, NULL, &AllCPs, &TermRadius, &RhoCutoff, VolInfo, HessPtrs, GradPtrs, RhoPtr);
								GP.SetStartEndCPNum(GPParams.StartCPNum, 0);
								GP.Seed(false);
								if (GP.GetStartEndCPNum(EndCPPosition) == TermCPNum){
									if (GPDir != StreamDir_Both) GP.PointPrepend(GPParams.StartPointOrigin, AllCPs.GetRho(TypeInd, cpNum));
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
							int Iter = 0;
							while (Iter < 50){
								Iter++;
								RotAngle += (AngleStep * AngleFactor * SmallAngleFactor * static_cast<double>(MinGPDir));
								StartPoint = GPParams.StartPointOrigin + Rotate(StartVec, RotAngle, GPParams.RotAxis);
								GradPath_c GP(StartPoint, GPDir, NumGPPts, GPType_Classic, GPTerminate_AtCP, NULL, &AllCPs, &TermRadius, &RhoCutoff, VolInfo, HessPtrs, GradPtrs, RhoPtr);
								GP.SetStartEndCPNum(GPParams.StartCPNum, 0);
								GP.Seed(false);
								if (GP.GetStartEndCPNum(EndCPPosition) == TermCPNum){
									if (GPDir != StreamDir_Both) GP.PointPrepend(GPParams.StartPointOrigin, AllCPs.GetRho(TypeInd, cpNum));
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

							if (Iter >= 50) AlphaLower = AlphaUpper;

							AlphaGuess = (AlphaLower + AlphaUpper) * 0.5; // guess is the midpoint between the two angles
						}

						if (AlphaUpper - AlphaLower >= SurfRCSMinAngleCheck){

							MinFunc_AlphaLowMidHigh.push_back(vec3());
							MinFunc_AlphaLowMidHigh.back() << AlphaLower << AlphaGuess << AlphaUpper;

							TermCPNums.push_back(TermCPNum);
						}
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

			gsl_min_fminimizer_type const *T;
			T = gsl_min_fminimizer_brent;

			gsl_min_fminimizer *s;
			s = gsl_min_fminimizer_alloc(T);

			// DEBUG check the GP lengths that will be found using the low, guess, and high angle values
// #ifdef _DEBUG
			vector<double> TmpLens(3);
			vector<int> TmpEndCPNums(3);
			vector<vec3> TmpStartPoints(3);
			for (int j = 0; j < 3; ++j){
				TmpStartPoints[j] = GPParams.StartPointOrigin + Rotate(StartVec, MinFunc_AlphaLowMidHigh[i][j], GPParams.RotAxis);
				GradPath_c GP(TmpStartPoints[j], GPDir, NumGPPts, GPType_Classic, GPTerminate_AtCP, NULL, &AllCPs, &TermRadius, &RhoCutoff, VolInfo, HessPtrs, GradPtrs, RhoPtr);

				GP.SetStartEndCPNum(GPParams.StartCPNum, 0);
				GP.Seed(false);
				if (GPDir != StreamDir_Both) GP.PointPrepend(GPParams.StartPointOrigin, AllCPs.GetRho(TypeInd, cpNum));

				TmpLens[j] = GP.GetLength();
				TmpEndCPNums[j] = GP.GetStartEndCPNum(EndCPPosition);
			}
			if (TmpLens[1] > TmpLens[0] || TmpLens[1] > TmpLens[2]){
				// initial guess and bounds for 1-d minimization is wrong, so skip minimization and shrink SGP vectors as necessary
				continue;
			}
// #endif

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
			SGPsPerCP[iCP][i].SetupGradPath(GPParams.StartPointOrigin + Rotate(StartVec, MinFunc_AlphaLowMidHigh[i][1], GPParams.RotAxis), GPDir, NumGPPts, GPType_Classic, GPTerminate_AtCP, NULL, &AllCPs, &TermRadius, &RhoCutoff, VolInfo, HessPtrs, GradPtrs, RhoPtr);
			SGPsPerCP[iCP][i].SetStartEndCPNum(GPParams.StartCPNum, 0);
			SGPsPerCP[iCP][i].Seed(false);
			if (GPDir != StreamDir_Both) SGPsPerCP[iCP][i].PointPrepend(GPParams.StartPointOrigin, AllCPs.GetRho(TypeInd, cpNum));
			SGPsPerCP[iCP][i].Resample(NumGPPts);

			SGPSeedAngles[i] = MinFunc_AlphaLowMidHigh[i][1];

			if (GPDir == StreamDir_Reverse) SGPsPerCP[iCP][i].Reverse();

			if (AllCPs.GetTypeFromTotOffset(TermCPNums[i]) == SecondaryTermType){
				/*
				*	Need to find the direction normal to the SGP at the bond/ring CP and in the plane of the ring/bond.
				*	Use unit vectors for 1) the last two points of the SGP and 2) a vector going from the originating
				*	bond/ring CP to another arbitrary point in the plane, to get three orthonormal unit vectors to use.
				*/
				vec3 pt1, pt2, pt3, v1, v2, v3;
				vec3 CPPos = AllCPs.GetXYZ(TermCPNums[i]);
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
				vec3 PtUp = CPPos + v2 * 0.05 * AllCPs.GetMinCPDist(CPTypesForMinDistCheck);
				vec3 PtDn = CPPos - v2 * 0.05 * AllCPs.GetMinCPDist(CPTypesForMinDistCheck);

				if (dot(v2, pt3 - GPParams.StartPointOrigin) > 0){
					/*
					 *	v2 points in the positive rotation direction
					 */

					SaddleCPSupplementGPs[i][0].SetupGradPath(PtUp, GPDir, NumGPPts,
						GPType_Classic, GPTerminate_AtRhoValue, NULL,
						&AllCPs, &TermRadius, &RhoCutoff, VolInfo,
						HessPtrs, GradPtrs, RhoPtr);
					SaddleCPSupplementGPs[i][1].SetupGradPath(PtDn, GPDir, NumGPPts,
						GPType_Classic, GPTerminate_AtRhoValue, NULL,
						&AllCPs, &TermRadius, &RhoCutoff, VolInfo,
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
					&AllCPs, &TermRadius, &RhoCutoff, VolInfo,
					HessPtrs, GradPtrs, RhoPtr);
					SaddleCPSupplementGPs[i][1].SetupGradPath(PtUp, GPDir, NumGPPts,
						GPType_Classic, GPTerminate_AtRhoValue, NULL,
						&AllCPs, &TermRadius, &RhoCutoff, VolInfo,
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
				StartPointOffset = RCSCheckRadiiHigh * AllCPs.GetMinCPDist(CPTypesForMinDistCheck);
			}
			else StartPointOffset = RCSCheckRadiiLow * AllCPs.GetMinCPDist(CPTypesForMinDistCheck);

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

						for (double const & iRCS : AvgMinNum){
							double aRCS = iRCS * AngleStep;
							StartPoint = GPParams.StartPointOrigin + Rotate(GPParams.RotVec, aRCS, GPParams.RotAxis);
							GradPath_c GP(StartPoint, GPDir, NumGPPts, GPType_Classic, GPTerminate_AtCP, NULL, &AllCPs, &TermRadius, &RhoCutoff, VolInfo, HessPtrs, GradPtrs, RhoPtr);
							GP.SetStartEndCPNum(AllCPs.GetTotOffsetFromTypeNumOffset(TypeInd, cpNum), 0);
							GP.Seed(false);
							GP.PointPrepend(GPParams.StartPointOrigin, AllCPs.GetRho(TypeInd, cpNum));
							GP.Resample(NumGPPts);
							if (GPDir == StreamDir_Reverse) GP.Reverse();
							bool GPFound = false;
							for (int iSGP = 0; iSGP < SGPsPerCP[iCP].size() && !GPFound; ++iSGP){
								if (SGPsPerCP[iCP][iSGP].IsMade()){
									double aSGP = SGPSeedAngles[iSGP];
									if (PeriodicDistance1D(aSGP, aRCS, 0., PI2) < RCSSGPRadCheck){
										GPFound = true;
										if (GP.GetLength() < SGPsPerCP[iCP][iSGP].GetLength()){
											SGPsPerCP[iCP][iSGP] = GP;
											SGPSeedAngles[iSGP] = aRCS;
										}
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

		{
			vector<GradPath_c> GoodSGPs;
			for (auto & sgp : SGPsPerCP[iCP]) if (sgp.IsMade()) GoodSGPs.push_back(sgp);
			SGPsPerCP[iCP] = GoodSGPs;
		}
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
				if (TermCPNum < 0 || AllCPs.GetTypeFromTotOffset(TermCPNum) == PrimaryTermType){
					SurfaceGPs.back().push_back(SGPsPerCP[iCP][SGPNumsNewToOld[i]]);
					SurfaceGPCPNums.back().push_back(SGPsPerCP[iCP][SGPNumsNewToOld[i]].GetStartEndCPNum(EndCPPosition));
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
					GradPath_c GP(StartPoint, GPDir, NumGPPts, GPType_Classic, GPTerminate_AtRhoValue, NULL, &AllCPs, &TermRadius, &RhoCutoff, VolInfo, HessPtrs, GradPtrs, RhoPtr);
					GP.SetStartEndCPNum(AllCPs.GetTotOffsetFromTypeNumOffset(TypeInd, cpNum), 0);
					GP.Seed(false);
					GP.PointPrepend(GPParams.StartPointOrigin, AllCPs.GetRho(TypeInd, cpNum));
					GP.Resample(NumGPPts);
					if (GPDir == StreamDir_Reverse) GP.Reverse();
// #ifdef _DEBUG
// 					GP.SaveAsOrderedZone("Filler GP");
// #endif
					SurfaceGPs.back().push_back(GP);
				}
				TermCPNum = -1;
				if (SGPNumsNewToOld[(i + 1) % SGPNumsNewToOld.size()] < TermCPNums.size()) TermCPNum = TermCPNums[SGPNumsNewToOld[(i + 1) % SGPNumsNewToOld.size()]];
				if (TermCPNum < 0 || AllCPs.GetTypeFromTotOffset(TermCPNum) == PrimaryTermType){
					SurfaceGPs.back().push_back(SGPsPerCP[iCP][SGPNumsNewToOld[(i + 1) % SGPNumsNewToOld.size()]]);
					SurfaceGPCPNums.back().push_back(SGPsPerCP[iCP][SGPNumsNewToOld[(i + 1) % SGPNumsNewToOld.size()]].GetStartEndCPNum(EndCPPosition));
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
				SurfaceGPCPNums.back().push_back(SGPsPerCP[iCP][0].GetStartEndCPNum(EndCPPosition));
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
					GradPath_c GP(StartPoint, GPDir, NumGPPts, GPType_Classic, GPTerminate_AtRhoValue, NULL, &AllCPs, &TermRadius, &RhoCutoff, VolInfo, HessPtrs, GradPtrs, RhoPtr);
					GP.SetStartEndCPNum(AllCPs.GetTotOffsetFromTypeNumOffset(TypeInd, cpNum), 0);
					GP.Seed(false);
					GP.PointPrepend(GPParams.StartPointOrigin, AllCPs.GetRho(TypeInd, cpNum));
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
				GradPath_c GP(StartPoint, GPDir, NumGPPts, GPType_Classic, GPTerminate_AtRhoValue, NULL, &AllCPs, &TermRadius, &RhoCutoff, VolInfo, HessPtrs, GradPtrs, RhoPtr);
				GP.SetStartEndCPNum(AllCPs.GetTotOffsetFromTypeNumOffset(TypeInd, cpNum), 0);
				GP.Seed(false);
				GP.PointPrepend(GPParams.StartPointOrigin, AllCPs.GetRho(TypeInd, cpNum));
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
					StartEndCPTypeAndOffset[i] = AllCPs.GetTypeNumOffsetFromTotOffset(StartEndCPNums[i]);
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
			WholeSurfaceZoneNameBase = (CPType == CPType_Bond ? "IAS_" : "RS_"),
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
			Surfaces.back().MakeFromGPs(SurfaceGPPtrs.back(), SingleSurface, false);

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
				SZFSStr += MakeStringFromCPNums(SurfaceGPCPNums[i], AllCPs, CornerCPTypeOffsets);
				int SurfZoneNum = FESurface_c::SetZoneStyle(Surfaces.back().SaveAsTriFEZone(XYZVarNums, ZoneNameBase + SZFSStr),
					AssignOp_MinusEquals);

				for (int c = 0; c < CornerCPTypeOffsets.size(); ++c){
					AuxDataZoneSetItem(SurfZoneNum, CSMAuxData.CC.ZFSCornerCPNumStrs[c], to_string(CornerCPTypeOffsets[c][1] + 1));
					if (CornerCPTypeOffsets[c][0] >= 0) AuxDataZoneSetItem(SurfZoneNum, CSMAuxData.CC.ZFSCornerCPTypes[c], CPNameList[CornerCPTypeOffsets[c][0]]);
					else AuxDataZoneSetItem(SurfZoneNum, CSMAuxData.CC.ZFSCornerCPTypes[c], "FF");
				}
				AuxDataZoneSetItem(SurfZoneNum, CSMAuxData.CC.ZoneType, CSMAuxData.CC.ZoneTypeSZFS);
				AuxDataZoneSetItem(SurfZoneNum, CSMAuxData.CC.ZoneSubType, (CPType == CPType_Ring ? CSMAuxData.CC.ZoneSubTypeRSSegment : CSMAuxData.CC.ZoneSubTypeIASSegment));
			}
		}

		if (WholeSurface.IsMade()){
			int SurfZoneNum = FESurface_c::SetZoneStyle(WholeSurface.SaveAsTriFEZone(XYZVarNums, WholeSurfaceZoneNameBase + CPNameBase),
				AssignOp_PlusEquals, TRUE, FALSE);

			AuxDataZoneSetItem(SurfZoneNum, CSMAuxData.CC.ZFSCornerCPNumStrs[0], to_string(AllCPs.GetTotOffsetFromTypeNumOffset(TypeInd, cpNum) + 1));
			AuxDataZoneSetItem(SurfZoneNum, CSMAuxData.CC.ZFSCornerCPTypes[0], CPNameList[TypeInd]);
			AuxDataZoneSetItem(SurfZoneNum, CSMAuxData.CC.ZoneType, CSMAuxData.CC.ZoneTypeSZFS);
			AuxDataZoneSetItem(SurfZoneNum, CSMAuxData.CC.ZoneSubType, (CPType == CPType_Ring ? CSMAuxData.CC.ZoneSubTypeRS : CSMAuxData.CC.ZoneSubTypeIAS));
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

	StartPointOffset = 0.05 * AllCPs.GetMinCPDist();

	for (int cpNum = 0; cpNum < NumCPs; ++cpNum){
		vec3 StartVec = normalise(AllCPs.GetEigVecs(TypeInd, cpNum).col(1)) * StartPointOffset;

		for (int gpNum = 0; gpNum < NumCircleCheckGPs; ++gpNum){
			double RotAngle = static_cast<double>(gpNum)* AngleStep;

			StartPoint = AllCPs.GetXYZ(TypeInd, cpNum) + Rotate(StartVec, RotAngle, AllCPs.GetPrincDir(TypeInd, cpNum));

			GPs.push_back(GradPath_c(StartPoint, GPDir, NumGPPts, GPType_Classic, GPTerminate_AtCP, NULL, &AllCPs, &TermRadius, &RhoCutoff, VolInfo, HessPtrs, GradPtrs, RhoPtr));
			GPs.back().SetStartEndCPNum(AllCPs.GetTotOffsetFromTypeNumOffset(TypeInd, cpNum), 0);

		}
// 		break;
	}

	int NumGPs = GPs.size();

#ifndef _DEBUG
#pragma omp parallel for schedule(dynamic)
#endif
	for (int iGP = 0; iGP < NumGPs; ++iGP){
		GPs[iGP].Seed(false);
		GPs[iGP].PointPrepend(AllCPs.GetXYZ(GPs[iGP].GetStartEndCPNum()[0]), AllCPs.GetRho(GPs[iGP].GetStartEndCPNum()[0]));
		GPs[iGP].Resample(NumGPPts);
		if (CPType == CPType_Bond) GPs[iGP].Reverse();
	}

	for (auto & GP : GPs){
		vector<int> StartEndCPNums = GP.GetStartEndCPNum();
		vector<vector<int> > StartEndCPTypeAndOffset(2);
		string Name = "GP_";
		for (int i = 0; i < 2; ++i){
			if (StartEndCPNums[i] >= 0){
				StartEndCPTypeAndOffset[i] = AllCPs.GetTypeNumOffsetFromTotOffset(StartEndCPNums[i]);
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
		AuxDataZoneSetItem(GP.GetZoneNum(), CSMAuxData.CC.ZoneSubType, (CPType == CPType_Bond ? CSMAuxData.CC.ZoneSubTypeBondCagePath : CSMAuxData.CC.ZoneSubTypeRingNuclearPath));
	}

	TecUtilDataLoadEnd();

	
	CSMGuiUnlock();
	TecUtilLockFinish(AddOnID);

	return TRUE;
}

void ConnectCPsReturnUserInfo(bool const GuiSuccess, 
	vector<GuiField_c> const & Fields,
	vector<GuiField_c> const PassthroughFields){
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

void DrawEigenvectotArrorsReturnUserInfo(bool const GuiSuccess, 
	vector<GuiField_c> const & Fields,
	vector<GuiField_c> const PassthroughFields){
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

	for (int z : SelectedZones){
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
		for (int v : EigValVarNums){
			EigValPtrs.back().push_back(FieldDataPointer_c());
			if (!EigValPtrs.back().back().GetReadPtr(z, v)){
				TecUtilDialogErrMsg("Failed to get eigenvalue data pointer");
				return;
			}
		}
		EigVecPtrs.push_back(vector<FieldVecPointer_c>());
		for (vector<int> const & e : EigVecVarNums) {
			EigVecPtrs.back().push_back(FieldVecPointer_c());
			if (!EigVecPtrs.back().back().GetReadPtr(z, e)){
				TecUtilDialogErrMsg("Failed to get eigenvector data pointer");
				return;
			}
		}
	}

	vector<vector<FESurface_c> > ArrowZones(NumSelectedZones, vector<FESurface_c>(3));

#ifndef _DEBUG
#pragma omp parallel for schedule(dynamic)
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

	for (string const & s : { "X", "Y", "Z" }) Fields.push_back(GuiField_c(Gui_VarSelect, s, s));
	for (string const & s : CSMVarName.EigVals) Fields.push_back(GuiField_c(Gui_VarSelect, s, s));
	for (string const & s : CSMVarName.EigVecs) Fields.push_back(GuiField_c(Gui_VarSelect, s, s));

	CSMGui("Connect CPs with lines", Fields, DrawEigenvectotArrorsReturnUserInfo, AddOnID);
}


bool const ExtractRadiusContourLinesToIOrderedPoints(vector<int> & ZoneNums, 
	vec3 const & Origin, double const & Radius, vector<int> const & XYZVarNums, vector<GradPath_c> & ContourLines)
{
	TecUtilLockStart(AddOnID);

	int NumZones = TecUtilDataSetGetNumZones();

	// First define a RadiusSquared variable for the specified zones

	Set ZoneSet;
	ArgList Args;

	for (int z : ZoneNums){
		REQUIRE(z > 0 && z <= NumZones);
		ZoneSet += z;
	}

	string TempVarName = "TempRadiusSqrVar";
	int TempVarNum = VarNumByName(TempVarName);

	if (TempVarNum <= 0){
		// Need to make new temp variable
		vector<FieldDataType_e> FDTypes(NumZones, FieldDataType_Bit);
		for (int z : ZoneNums) FDTypes[z - 1] = FieldDataType_Double;
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

	vector<int> IntZoneNums;

	for (int z : ZoneNums){
		TecUtilZoneSetActive(Set(z).getRef(), AssignOp_Equals);

		Args.clear();
		Args.appendInt(SV_CONTLINECREATEMODE, ContLineCreateMode_OneZonePerIndependentPolyline);
		TecUtilCreateContourLineZonesX(Args.getRef());

		int NewZoneNum = TecUtilDataSetGetNumZones();

		if (NewZoneNum - NumZones == 1){
			ContourLines.push_back(GradPath_c(NewZoneNum, XYZVarNums, AddOnID));

			ZoneSet += NewZoneNum;
			NumZones = NewZoneNum;
			IntZoneNums.push_back(z);
		}
	}

	ZoneNums = IntZoneNums;

// 	TecUtilDataSetDeleteVar(Set(TempVarNum).getRef());
// 	TecUtilDataSetDeleteZone(ZoneSet.getRef());
	TecUtilFrameDeleteActive();

	TecUtilLockFinish(AddOnID);

	return true;
}

void ExtractRSIntersectionsReturnUserInfo(bool const GuiSuccess,
	vector<GuiField_c> const & Fields,
	vector<GuiField_c> const PassthroughFields){
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

	ExtractSurfaceSphereIntersections(VolZoneNum, CPTypeVarNum, AllCPZoneNum, CPList, CPZoneNum, Radius, ResampleNumPoints, XYZVarNums, vector<GradPath_c>());
}

void ExtractSurfaceSphereIntersections(
	int VolZoneNum,
	int CPTypeVarNum,
	int AllCPZoneNum,
	vector<int> const & CPList,
	int CPZoneNum,
	double const & Radius,
	int ResampleNumPoints,
	vector<int> const & XYZVarNums,
	vector<GradPath_c> & Intersections)
{
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

	for (int CPNum1 : CPList){
		vector<int> TypeNumOffset = (CPZoneTypeInd < 0 ? CPs.GetTypeNumOffsetFromTotOffset(CPNum1 - 1) : vector<int>({ CPZoneTypeInd, CPNum1 - 1 }));
		CPType_e CPType = CPTypeList[TypeNumOffset[0]];
		if (CPType != CPType_Nuclear && CPType != CPType_Cage) continue;

		int CPTotOffset = CPs.GetTotOffsetFromTypeNumOffset(TypeNumOffset[0], TypeNumOffset[1]) + 1;

		// Get zone numbers of intersecting ring (interatomic) surfaces
		string SurfaceSearchString = (CPType == CPType_Nuclear ? CSMAuxData.CC.ZoneSubTypeRS : CSMAuxData.CC.ZoneSubTypeIAS);
		vector<int> IntersectingZoneNums;
		for (int z = 1; z <= NumZones; ++z){
			if (AuxDataZoneItemMatches(z, CSMAuxData.CC.ZoneSubType, SurfaceSearchStrings[CPType == CPType_Nuclear ? 0 : 1])){
				for (string const & CornerStr : CSMAuxData.CC.ZFSCornerCPNumStrs){
					int CornerNum = stoi(AuxDataZoneGetItem(z, CornerStr));
					if (CornerNum == CPTotOffset){
						IntersectingZoneNums.push_back(z);
						break;
					}
				}
			}
			else if (AuxDataZoneItemMatches(z, CSMAuxData.CC.ZoneSubType, PathSearchStrings[CPType == CPType_Nuclear ? 0 : 1])){
				for (string const & CornerStr : CSMAuxData.CC.GPEndNumStrs){
					int CornerNum = stoi(AuxDataZoneGetItem(z, CornerStr));
					if (CornerNum == CPTotOffset){
						IntersectingZoneNums.push_back(z);
						break;
					}
				}
			}
		}

		// Extract intersection lines as gradient paths
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

	for (string const & s : { "X", "Y", "Z" }) Fields.push_back(GuiField_c(Gui_VarSelect, s, s));

	CSMGui("Extract Ring Surface Sphere Intersections", Fields, ExtractRSIntersectionsReturnUserInfo, AddOnID);
}



void GBAReturnUserInfo(bool const GuiSuccess,
	vector<GuiField_c> const & Fields,
	vector<GuiField_c> const PassthroughFields)
{
	if (!GuiSuccess) return;

	BondalyzerGetUserInfo(BondalyzerCalcType_GBA, Fields);
}

void GBAGetUserInfo()
{
	vector<GuiField_c> Fields({
		GuiField_c(Gui_ZonePointSelectMulti, "Select CPs", CSMZoneName.CPType[0]),
		GuiField_c(Gui_VarSelect, "CP type variable", CSMVarName.CritPointType),
		GuiField_c(Gui_Label, "Sphere mesh:"),
		GuiField_c(Gui_Double, "Radius", "0.2"),
		GuiField_c(Gui_Radio, "Radius is", "Absolute,Fraction of min CP dist."),
		GuiField_c(Gui_Int, "Icosahedron refinement", "3"),
		GuiField_c(Gui_Label, "Gradient bundles (GBs):"),
		GuiField_c(Gui_Double, "Rho cutoff value", "0.001"),
		GuiField_c(Gui_Int, "# GB edge points", "3"),
		GuiField_c(Gui_Int, "# points per gradient path", "50")
	});

	CSMGui("Gradient bundle analysis", Fields, GBAReturnUserInfo, AddOnID);
}

Boolean_t GBA_Generation(
	int VolZoneNum,
	int CPZoneNum,
	int CPTypeVarNum,
	vector<int> const & XYZVarNums,
	int RhoVarNum,
	vector<int> const & GradVarNums,
	vector<int> const & HessVarNums,
	Boolean_t IsPeriodic,
	vector<int> const & SelectedCPNums,
	double const & SphereRadius,
	int RadiusMode,
	int RefinementLevel,
	double const & RhoCutoff,
	int NumGBEdgePoints,
	int GPNumPoints
	)
{
	TecUtilLockStart(AddOnID);

	int NumZones = TecUtilDataSetGetNumZones();
	int NumVars = TecUtilDataSetGetNumVars();

	REQUIRE(CPZoneNum > 0 && CPZoneNum <= NumZones);
	REQUIRE(XYZVarNums.size() == 3);
	for (auto const & i : XYZVarNums) REQUIRE(i > 0 && i <= NumVars);

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

	TecUtilZoneGetIJK(CPZoneNum, &NumCPs, &iJunk[0], &iJunk[1]);
	for (auto const & i : iJunk) if (i > 1){
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
	MR.BasisVectors = &VolInfo.BasisNormalized;


	CritPoints_c CPs(CPZoneNum, XYZVarNums, CPTypeVarNum, RhoVarNum, &MR);

	CSMGuiLock();

	// This is where the magic happens!

	CSMGuiUnlock();

	TecUtilLockFinish(AddOnID);
}

/*
 *	Given a list of SelectedCPNums in zone SelectedCPZoneNum, find,
 *	by comparing XYZ coordinates, the same CPs in zone AllCPsZoneNum.
 */
Boolean_t CPNumbersMapBetweenZones(int AllCPsZoneNum,
	int SelectedCPZoneNum,
	vector<int> const & XYZVarNums,
	vector<int> const & SelectedCPNums,
	vector<int> & MappedCPNums)
{
	TecUtilLockStart(AddOnID);

	int NumZones = TecUtilDataSetGetNumZones();
	int NumVars = TecUtilDataSetGetNumVars();

	REQUIRE(AllCPsZoneNum > 0 && AllCPsZoneNum <= NumZones);
	REQUIRE(SelectedCPZoneNum > 0 && SelectedCPZoneNum <= NumZones);
	REQUIRE(XYZVarNums.size() == 3);
	for (auto const & i : XYZVarNums) REQUIRE(i > 0 && i <= NumVars);

	FieldVecPointer_c SelXYZ, AllXYZ;

	VERIFY(SelXYZ.GetReadPtr(SelectedCPZoneNum, XYZVarNums) && AllXYZ.GetReadPtr(AllCPsZoneNum, XYZVarNums));

	for (int i : SelectedCPNums) REQUIRE(i <= SelXYZ.Size()); // i are base-1 indices

	MappedCPNums.resize(SelectedCPNums.size());

	Boolean_t AllFound = TRUE;

	for (int iCP = 0; iCP < SelectedCPNums.size(); ++iCP){
		// First assume that the CPs are already mapped, so first check the same cp index in the full cp zone
		int SelInd = SelectedCPNums[iCP] - 1; // moving from base-1 to base-0
		if (SelInd >= AllXYZ.Size() || sum(SelXYZ[SelInd] == AllXYZ[SelInd]) != 3){
			// CP didn't match between zones, so start searching.
			// For first CP, start sequential search from first CP.
			// For rest of CPs, start search from the CP following the previous mapped CP number.
			
			Boolean_t CPFound = FALSE;
			int StartSearchInd = (iCP > 0 ? MappedCPNums[iCP - 1] + 1 : 0);
			for (int jCP = StartSearchInd; jCP < StartSearchInd + AllXYZ.Size() && !CPFound; ++jCP){
				int jInd = jCP % AllXYZ.Size();
				if (sum(SelXYZ[SelInd] == AllXYZ[jInd]) == 3){
					MappedCPNums[iCP] = jInd + 1; // moving from base-0 to base-1
					CPFound = TRUE;
				}
			}

			if (AllFound && !CPFound) AllFound = FALSE;
		}
	}

	TecUtilLockFinish(AddOnID);

	return AllFound;
}

Boolean_t BondalyzerBatch(int VolZoneNum,
	vector<int> const & CPZoneNums,
	int CPTypeVarNum,
	vector<int> const & XYZVarNums,
	int RhoVarNum,
	vector<int> const & GradVarNums,
	vector<int> const & HessVarNums,
	Boolean_t IsPeriodic,
	int RidgeFuncVarNum,
	vector<bool> const & CalcSteps)
{
	// Check and run each calculation step. (minus 1 because critical points isn't included)
	if (CalcSteps[(int)BondalyzerCalcType_BondPaths - 1]){
		FindBondRingLines(VolZoneNum, CPZoneNums, CPZoneNums[0], vector<int>(), CPTypeVarNum, CPType_Bond, XYZVarNums, RhoVarNum, GradVarNums, HessVarNums, IsPeriodic);
	}
	if (CalcSteps[(int)BondalyzerCalcType_RingLines - 1]){
		FindBondRingLines(VolZoneNum, CPZoneNums, CPZoneNums[0], vector<int>(), CPTypeVarNum, CPType_Ring, XYZVarNums, RhoVarNum, GradVarNums, HessVarNums, IsPeriodic);
	}
	if (CalcSteps[(int)BondalyzerCalcType_CageNuclearPaths - 1]){
		FindCageNuclearPaths(VolZoneNum, CPZoneNums, CPZoneNums[0], vector<int>(), CPTypeVarNum, XYZVarNums, RhoVarNum, GradVarNums, HessVarNums, IsPeriodic);
	}
	if (CalcSteps[(int)BondalyzerCalcType_InteratomicSurfaces - 1]){
		FindBondRingSurfaces(VolZoneNum, CPZoneNums, CPZoneNums[0], vector<int>(), CPTypeVarNum, CPType_Bond, XYZVarNums, RhoVarNum, GradVarNums, HessVarNums, RidgeFuncVarNum, IsPeriodic);
	}
	if (CalcSteps[(int)BondalyzerCalcType_RingSurfaces - 1]){
		FindBondRingSurfaces(VolZoneNum, CPZoneNums, CPZoneNums[0], vector<int>(), CPTypeVarNum, CPType_Ring, XYZVarNums, RhoVarNum, GradVarNums, HessVarNums, RidgeFuncVarNum, IsPeriodic);
	}

	return TRUE;
}

void CombineCPZonesReturnUserInfo(bool const GuiSuccess,
	vector<GuiField_c> const & Fields,
	vector<GuiField_c> const PassthroughFields)
{
	int fNum = 0;
	vector<int> ZoneNums = Fields[fNum++].GetReturnIntVec();
	bool DeleteSourceZones = Fields[fNum++].GetReturnBool();
	int RhoVarNum = Fields[fNum++].GetReturnInt();
	int CPTypeVarNum = Fields[fNum++].GetReturnInt();
	vector<int> XYZVarNums(3);
	for (int i = 0; i < 3; ++i) XYZVarNums[i] = Fields[fNum].GetReturnInt() + i;


	TecUtilLockStart(AddOnID);

	CritPoints_c CPs;

	Set CPZones;

	for (int z : ZoneNums) {
		CPZones += z;
		CPs += CritPoints_c(z, XYZVarNums, CPTypeVarNum, RhoVarNum);
	}

	CPs.RemoveSpuriousCPs();

	CPs.SaveAsOrderedZone(XYZVarNums, RhoVarNum, TRUE);

	if (DeleteSourceZones) TecUtilZoneDelete(CPZones.getRef());

	TecUtilLockFinish(AddOnID);
}

void CombineCPZonesGetUserInfo(){
	vector<GuiField_c> Fields = {
		GuiField_c(Gui_ZoneSelectMulti, "CP Zones"),
		GuiField_c(Gui_Toggle, "Delete source CP zones?", "1"),
		GuiField_c(Gui_VarSelect, "Electron density variable", CSMVarName.Dens),
		GuiField_c(Gui_VarSelect, "CP type variable", CSMVarName.CritPointType),
		GuiField_c(Gui_VarSelect, "X", "X")
	};

	CSMGui("Combine CP zones", Fields, CombineCPZonesReturnUserInfo, AddOnID);
}


void MakeSurfaceFromPathZonesReturnUserInfo(bool const GuiSuccess,
	vector<GuiField_c> const & Fields,
	vector<GuiField_c> const PassthroughFields)
{
	int fNum = 0;
	bool UseStreamtraces = Fields[fNum++].GetReturnBool();
	vector<int> ZoneNums = Fields[fNum++].GetReturnIntVec();
	int ResampleNumPoints = Fields[fNum++].GetReturnInt();
	int GroupNum = Fields[fNum++].GetReturnInt();
	bool ConnectFirstAndLastPaths = Fields[fNum++].GetReturnBool();
	bool DelSourceZones = Fields[fNum++].GetReturnBool();
	vector<int> XYZVarNums(3);
	for (int i = 0; i < 3; ++i) XYZVarNums[i] = Fields[fNum++].GetReturnInt();

	TecUtilLockStart(AddOnID);

	if (UseStreamtraces && TecUtilStreamtraceGetCount() > 1 && TecUtilStreamtracesAreActive()) {
		ZoneNums.clear();
		int startZoneNum = TecUtilDataSetGetNumZones();
		if (TecUtilCreateStreamZones(FALSE)){
			for (int z = startZoneNum+1; z <= TecUtilDataSetGetNumZones(); ++z){
				ZoneNums.push_back(z);
			}
		}
		if (ZoneNums.back() - ZoneNums[0] < 2){
			TecUtilDialogErrMsg("Not enough streamtrace zones were created");
			TecUtilLockFinish(AddOnID);
			return;
		}
		DelSourceZones = true;
	}



	if (ResampleNumPoints > 0){
		/*
		 *	First, read in all zones as GPs to resample.
		 */
		vector<GradPath_c> GPs;
		vector<GradPath_c*> GPPtrs;
		vector<int> XYZRhoVarNums = XYZVarNums;
		XYZRhoVarNums.push_back(4);
		for (int z : ZoneNums){
			GPs.push_back(GradPath_c(z, XYZRhoVarNums, AddOnID));
			GPs.back().Resample(ResampleNumPoints);
		}
		for (auto & gp : GPs)
			GPPtrs.push_back(&gp);

		/*
		 *	Now make and save surface
		 */
		FESurface_c Surf;
		Surf.MakeFromGPs(GPPtrs, ConnectFirstAndLastPaths);

		if (Surf.IsMade()){
			Surf.SaveAsTriFEZone(XYZVarNums, "Stitched path surface");
		}
	}
	else {
		/*
		 *	First, need to read in all the points of all the zones to a single vector of vec3
		 *	while recording the indices in that vector of each individual path.
		 */
		vector<vec3> P;
		vector<vector<int> > IndList;
		int Ind = 0;
		for (int z : ZoneNums) {
			REQUIRE(TecUtilZoneIsOrdered(z));
			int IJK[3];
			TecUtilZoneGetIJK(z, &IJK[0], &IJK[1], &IJK[2]);
			REQUIRE(IJK[0] > 1 && IJK[1] == 1 && IJK[2] == 1);

			/*
			 *	Preallocate space in P assuming all paths are same number of points.
			 */
			if (P.size() == 0) {
				P.reserve(IJK[0] * ZoneNums.size());
			}
			IndList.push_back(vector<int>(IJK[0]));

			FieldVecPointer_c XYZVarPtr;
			if (XYZVarPtr.GetReadPtr(z, XYZVarNums)) {
				for (int i = 0; i < IJK[0]; ++i) {
					P.push_back(XYZVarPtr[i]);
					IndList.back()[i] = Ind++;
				}
			}
		}
		if (ConnectFirstAndLastPaths)
			IndList.push_back(IndList[0]);

		/*
		 *	Second, construct the triangular elements by stitching neighboring paths together.
		 */
		vector<vector<int> > Elem;
		for (int i = 0; i < IndList.size() - 1; ++i) {
			StitchPaths(IndList[i], IndList[i + 1], P, Elem);
		}

		/*
		 * Now create and save surface
		 */
		FESurface_c Surf(P, Elem);
		if (Surf.IsMade()) {
			Surf.SaveAsTriFEZone(XYZVarNums, "Stitched path surface");
		}
	}

	/*
	 *	Deactivate source zones and activate the new zone
	 */
	Set sourceZones;
	for (int z : ZoneNums) sourceZones += z;
	TecUtilZoneSetActive(sourceZones.getRef(), AssignOp_MinusEquals);
	TecUtilZoneSetActive(Set(TecUtilDataSetGetNumZones()).getRef(), AssignOp_PlusEquals);

	/*
	 *	Delete source zones
	 */
	if (DelSourceZones) {
		Set delZones;
		for (int z : ZoneNums) delZones += z;
		TecUtilDataSetDeleteZone(delZones.getRef());
		ZoneNums.clear();
	}

	/*
	 *	Change group number of zones
	 */
	if (GroupNum > 1) {
		ZoneNums.push_back(TecUtilDataSetGetNumZones());
		string ZoneStr = "[" + to_string(ZoneNums[0]);
		for (int z = 1; z < ZoneNums.size(); ++z) ZoneStr += "," + to_string(ZoneNums[z]);
		ZoneStr += "]";
		TecUtilMacroExecuteCommand(string("$!FIELDMAP " + ZoneStr + " GROUP = " + to_string(GroupNum)).c_str());
	}

	if (UseStreamtraces){
		TecUtilStreamtraceDeleteAll();
	}

	

	TecUtilLockFinish(AddOnID);
}

void MakeSurfaceFromPathZonesGetUserInfo() {
	vector<GuiField_c> Fields = {
		GuiField_c(Gui_Toggle, "Use current set of streamtraces"),
		GuiField_c(Gui_ZoneSelectMulti, "Path Zones"),
		GuiField_c(Gui_Int, "Resample with at most # points (0 for no resample)", "500"),
		GuiField_c(Gui_Int, "Group Number", "1"),
		GuiField_c(Gui_Toggle, "Connect first and last paths?"),
		GuiField_c(Gui_Toggle, "Delete source zones after surface creation?", "1"),
		GuiField_c(Gui_VarSelect, "X", "X"),
		GuiField_c(Gui_VarSelect, "Y", "Y"),
		GuiField_c(Gui_VarSelect, "Z", "Z")
	};

	CSMGui("Make surface from path zones", Fields, MakeSurfaceFromPathZonesReturnUserInfo, AddOnID);
}


void MakeSliceFromPointSelectionReturnUserInfo(bool const GuiSuccess,
	vector<GuiField_c> const & Fields,
	vector<GuiField_c> const PassthroughFields){
	if (!GuiSuccess) return;

	TecUtilLockStart(AddOnID);

	int fNum = 1;

	vector<int> SelectedCPNums = Fields[fNum].GetReturnIntVec();
	vector<string> SelectedCPNames = Fields[fNum].GetReturnStringVec();
	int CPZoneNum = Fields[fNum++].GetReturnInt();
	vector<int> VolumeZoneNums = Fields[fNum++].GetReturnIntVec();
	int SelectedPlane = Fields[fNum++].GetReturnInt();

	if (SelectedCPNums.size() < 1 || SelectedCPNums.size() > 3){
		TecUtilDialogErrMsg("Select 1-3 points");
		return;
	}

// 	CSMGuiLock();
	TecUtilLockStart(AddOnID);

	Set VolumeZones;
	Set_pa ActiveZones = TecUtilSetAlloc(FALSE);
	TecUtilZoneGetActive(&ActiveZones);
	for (int z : VolumeZoneNums) VolumeZones += z;
	TecUtilZoneSetActive(VolumeZones.getRef(), AssignOp_Equals);
	

	vector<int> XYZVarNums = { -1, -1, -1 };
	TecUtilAxisGetVarAssignments(&XYZVarNums[0], &XYZVarNums[1], &XYZVarNums[2]);
	if (XYZVarNums[0] < 0 || XYZVarNums[1] < 0 || XYZVarNums[2] < 0) for (int i = 1; i <= 3; ++i) XYZVarNums[i - 1] = i;

	FieldVecPointer_c XYZPtr;

	TecUtilDataLoadBegin();

	XYZPtr.GetReadPtr(CPZoneNum, XYZVarNums);

	vec3 Pt = zeros<vec>(3);
	string ZoneName = "Slc: ";
	for (int i = 0; i < SelectedCPNames.size(); ++i){
		ZoneName += SelectedCPNames[i][0];
		ZoneName += SelectedCPNames[i].back();
		if (i < SelectedCPNames.size() - 1) ZoneName += "-";
	}
	ZoneName += " (Zone " + to_string(CPZoneNum) + ")";

	for (int & i : SelectedCPNums) Pt += XYZPtr[--i];
	Pt /= (double)SelectedCPNums.size();

	vec3 Normal = zeros<vec>(3);

	if (SelectedCPNums.size() == 1){
		Normal[SelectedPlane - 1] = 1;
	}
	else if (SelectedCPNums.size() == 2){
		Normal = XYZPtr[SelectedCPNums[1]] - XYZPtr[SelectedCPNums[0]];
	}
	else if (SelectedCPNums.size() == 3){
		Normal = cross(XYZPtr[SelectedCPNums[1]] - XYZPtr[SelectedCPNums[0]], XYZPtr[SelectedCPNums[2]] - XYZPtr[SelectedCPNums[0]]);
	}

	TecUtilDataLoadEnd();

	TecUtilCreateSliceZoneFromPlane(SliceSource_VolumeZones,
		Pt[0], Pt[1], Pt[2],
		Normal[0], Normal[1], Normal[2]);

	Set NewZone(TecUtilDataSetGetNumZones());
	TecUtilZoneSetMesh(SV_SHOW, NewZone.getRef(), 0.0, FALSE);
	TecUtilZoneSetScatter(SV_SHOW, NewZone.getRef(), 0.0, FALSE);
	TecUtilZoneSetShade(SV_SHOW, NewZone.getRef(), 0.0, FALSE);
	TecUtilZoneSetContour(SV_SHOW, NewZone.getRef(), 0.0, TRUE);

	TecUtilZoneRename(TecUtilDataSetGetNumZones(), ZoneName.c_str());

	TecUtilSetAddMember(ActiveZones, TecUtilDataSetGetNumZones(), FALSE);
	TecUtilZoneSetActive(ActiveZones, AssignOp_Equals);

	TecUtilSetDealloc(&ActiveZones);

	TecUtilLockFinish(AddOnID);
// 	CSMGuiUnlock();
}

void MakeSliceFromPointSelectionGetUserInfo(){
	vector<GuiField_c> Fields = {
		GuiField_c(Gui_Label, "Define plane by 2 (midpoint) or more (avg. normal) points"),
		GuiField_c(Gui_ZonePointSelectMulti, "Critical point(s)", CSMZoneName.CriticalPoints),
		GuiField_c(Gui_ZoneSelectMulti, "Volume zone(s)", CSMZoneName.FullVolume.substr(0,10)),
		GuiField_c(Gui_Radio, "Plane (if single point)", "X,Y,Z")
	};

	CSMGui("Make slice from selected points", Fields, MakeSliceFromPointSelectionReturnUserInfo, AddOnID);
}



void TestFunction(){
	TecUtilLockStart(AddOnID);
	/*
	*	Quick test of tecplot toolbox set class
	*/
	// 	TecUtilZoneSetActive(Set(1).getRef(), AssignOp_Equals);

	/*
	*	Grad path test for adaptive step size
	*/
// 	vec3 StartPoint;
// 	StartPoint << 0.1 << 0.1 << 0.1;
// 
// 	VolExtentIndexWeights_s VolInfo;
// 	GetVolInfo(1, { 1, 2, 3 }, FALSE, VolInfo);
// 
// 	vector<FieldDataPointer_c> GradPtrs(3), HessPtrs(6);
// 	FieldDataPointer_c RhoPtr;
// 
// 	int VarNum = 4;
// 
// 	RhoPtr.GetReadPtr(1, VarNum++);
// 
// 	for (int i = 0; i < 3; ++i) GradPtrs[i].GetReadPtr(1, VarNum++);
// 	for (int i = 0; i < 6; ++i) HessPtrs[i].GetReadPtr(1, VarNum++);
// 
// 	double RhoCutoff = 1e-3;
// 	GradPath_c GP(StartPoint, StreamDir_Forward, 200, GPType_Classic, GPTerminate_AtRhoValue, NULL, NULL, NULL, &RhoCutoff, VolInfo, HessPtrs, GradPtrs, RhoPtr);
// 
// 	GP.Seed(true);
// 
// 	GP.SaveAsOrderedZone("test GP");
// 	

	/*
	 *	Testing the creation and initialzation of a fe_quad zone, because this part is crashing in the isosurface extraction tool.
	 *	This works, so wonder why it's crashing in the other function...
	 */

// 	FieldVecPointer_c XYZPtr;
// 	vector<vec3> Nodes(4);
// 	for (int i = 0; i < 2; ++i) for (int j = 0; j < 2; ++j) Nodes[i * 2 + j] << i << j << 0;
// 	vector<int> Elems = { 1, 2, 4, 3 };
// 
// 	TecUtilDataSetAddZone("test zone", 4, 1, 4, ZoneType_FEQuad, NULL);
// 	int ZoneNum = TecUtilDataSetGetNumZones();
// 
// 	TecUtilDataLoadBegin();
// 
// 	XYZPtr.GetWritePtr(ZoneNum, { 1, 2, 3 });
// 
// 	for (int i = 0; i < 4; ++i){
// 		TecUtilDataNodeSetByZone(ZoneNum, 1, i + 1, Elems[i]);
// 		XYZPtr.Write(i, Nodes[i]);
// 	}
// 	
// 	TecUtilDataLoadEnd();
// 	

	/*
	 *	Test rolling average-based time remaining on status
	 */
// 
// 	high_resolution_clock::time_point StartTime = high_resolution_clock::now();
// 
// 	StatusLaunch("test", AddOnID, TRUE, TRUE);
// 
// 	int N = 200;
// 
// 	for (int i = 0; i < N; ++i){
// 		StatusUpdate(i, N, "test", AddOnID, StartTime);
// 		Sleep(100 * (i % 5));
// 	}
// 
// 	StatusDrop(AddOnID);
// 	

	/*
	 *	Test isoRho value based integration
	 */

	/*
	 *	First on a typical GB, starting at the sphere end terminating not at a point.
	 */
// 	vector<int> GPZoneNums = { 55, 191, 60, 126, 30, 124 };
// 	vector<GradPath_c*> GPPtrs;
// 	for (int i : GPZoneNums){
// 		GPPtrs.push_back(new GradPath_c);
// 		*GPPtrs.back() = GradPath_c(i, { 1, 2, 3, 4 }, AddOnID);
// 	}
// 	VolExtentIndexWeights_s VolInfo;
// 	GetVolInfo(1, { 1, 2, 3 }, FALSE, VolInfo);
// 	vector<FieldDataPointer_c> VarPtrs(1);
// 	VarPtrs[0].GetReadPtr(1, 4);
// 	vector<double> IntVals;
// 	IntegrateUsingIsosurfaces(GPPtrs, 3, VolInfo, VarPtrs, IntVals);
// 
// 	/*
// 	 *	Get the sphere so we can add the sphere-interior portion of the integral
// 	 */
// 	FESurface_c Sphere(22, 1, { 1, 2, 3 }, { 4 });
// 	Sphere.DoIntegrationNew(2, TRUE);
// 	vector<double> SphereTriangleAreas;
// 	vector<vector<double> > SphereElemIntVals = Sphere.GetTriSphereIntValsByElem(&SphereTriangleAreas);
// 
// 	for (int i = 0; i < IntVals.size(); ++i){
// 		IntVals[i] += SphereElemIntVals[62][i];
// 
// 		// 		TecUtilDialogMessageBox(to_string(IntVals[i]).c_str(), MessageBoxType_Information);
// 	}
// 
// 	/*
// 	 *	Now a bond-path-coincident GB, again not terminating at a point.
// 	 */
// 
// 	GPZoneNums = { 73, 278, 74, 94, 48, 177, 47, 93 };
// 	GPPtrs.clear();
// 	for (int i : GPZoneNums){
// 		GPPtrs.push_back(new GradPath_c);
// 		*GPPtrs.back() = GradPath_c(i, { 1, 2, 3, 4 }, AddOnID);
// 	}
// 	IntegrateUsingIsosurfaces(GPPtrs, 2, VolInfo, VarPtrs, IntVals, true);
// 
// 	for (int i = 0; i < IntVals.size(); ++i){
// 		IntVals[i] += SphereElemIntVals[20][i];
// 
// 		// 		TecUtilDialogMessageBox(to_string(IntVals[i]).c_str(), MessageBoxType_Information);
// 	}
// 	

	/*
	 *	Test project to triangle function
	 */
// vec3 OP, TP;
// vector<vec3> tri;
// int numIter = 3;
// 
// vector<vec3> axes = { {1,0,0},{0,1,0},{0,0,1} };
// for (auto const & i : axes) {
// 	SaveVec3VecAsScatterZone({ i * 2, {0,0,0} }, string("axis"), Green_C);
// 	TecUtilZoneSetScatter(SV_SHOW, tecplot::toolbox::Set(TecUtilDataSetGetNumZones()).getRef(), 0, FALSE);
// 	TecUtilZoneSetMesh(SV_SHOW, tecplot::toolbox::Set(TecUtilDataSetGetNumZones()).getRef(), 0, TRUE);
// 	TecUtilZoneSetMesh(SV_LINETHICKNESS, tecplot::toolbox::Set(TecUtilDataSetGetNumZones()).getRef(), 0.2, 0);
// 	TecUtilZoneSetMesh(SV_COLOR, tecplot::toolbox::Set(TecUtilDataSetGetNumZones()).getRef(), 0, Black_C);
// 	TecUtilZoneSetActive(tecplot::toolbox::Set(TecUtilDataSetGetNumZones()).getRef(), AssignOp_MinusEquals);
// }
// 
// tri = { { 1,0,0 }, { 0,1,0 }, { 0,0,1 } };
// for (int i = 0; i < numIter; ++i) {
// 	OP = { i * 1. / (double)numIter, 0, 0 };
// 	for (int j = 0; j < numIter; ++j) {
// 		OP[1] += 1. / (double)numIter;
// 		for (int k = 0; k < numIter; ++k) {
// 			OP[2] += 1. / (double)numIter;
// 
// 			SaveVec3VecAsScatterZone({ tri[0], tri[1], tri[2], tri[0] }, string("triangle"));
// 			TecUtilZoneSetScatter(SV_SHOW, tecplot::toolbox::Set(TecUtilDataSetGetNumZones()).getRef(), 0, FALSE);
// 			TecUtilZoneSetMesh(SV_SHOW, tecplot::toolbox::Set(TecUtilDataSetGetNumZones()).getRef(), 0, TRUE);
// 			TecUtilZoneSetMesh(SV_LINETHICKNESS, tecplot::toolbox::Set(TecUtilDataSetGetNumZones()).getRef(), 0.2, 0);
// 			TecUtilZoneSetMesh(SV_COLOR, tecplot::toolbox::Set(TecUtilDataSetGetNumZones()).getRef(), 0, Red_C);
// 			TecUtilZoneSetActive(tecplot::toolbox::Set(TecUtilDataSetGetNumZones()).getRef(), AssignOp_MinusEquals);
// 
// 			double testProj = PointDistanceToTriangleSquared(OP, TP, tri[0], tri[1], tri[2]);
// 			SaveVec3VecAsScatterZone({ OP }, string("old point"), Red_C);
// 			TecUtilZoneSetScatter(SV_FRAMESIZE, tecplot::toolbox::Set(TecUtilDataSetGetNumZones()).getRef(), 1, 0);
// 			TecUtilZoneSetActive(tecplot::toolbox::Set(TecUtilDataSetGetNumZones()).getRef(), AssignOp_MinusEquals);
// 			SaveVec3VecAsScatterZone({ TP }, string("new point"), Green_C);
// 			TecUtilZoneSetScatter(SV_FRAMESIZE, tecplot::toolbox::Set(TecUtilDataSetGetNumZones()).getRef(), 1, 0);
// 			TecUtilZoneSetActive(tecplot::toolbox::Set(TecUtilDataSetGetNumZones()).getRef(), AssignOp_MinusEquals);
// 
// 			SaveVec3VecAsScatterZone({ OP, TP }, string("old to new point"), Green_C);
// 			TecUtilZoneSetScatter(SV_SHOW, tecplot::toolbox::Set(TecUtilDataSetGetNumZones()).getRef(), 0, FALSE);
// 			TecUtilZoneSetMesh(SV_SHOW, tecplot::toolbox::Set(TecUtilDataSetGetNumZones()).getRef(), 0, TRUE);
// 			TecUtilZoneSetMesh(SV_LINETHICKNESS, tecplot::toolbox::Set(TecUtilDataSetGetNumZones()).getRef(), 0.2, 0);
// 			TecUtilZoneSetMesh(SV_COLOR, tecplot::toolbox::Set(TecUtilDataSetGetNumZones()).getRef(), 0, Blue_C);
// 			TecUtilZoneSetActive(tecplot::toolbox::Set(TecUtilDataSetGetNumZones()).getRef(), AssignOp_MinusEquals);
// 		}
// 	}
// }
// 
// tri = { { 0,0,0 }, { 0,1,0 }, { 1,0,0 } };
// for (int i = 0; i < numIter; ++i) {
// 	OP = { i * 1. / (double)numIter, 0, 0 };
// 	for (int j = 0; j < numIter; ++j) {
// 		OP[1] += 1. / (double)numIter;
// 		for (int k = 0; k < numIter; ++k) {
// 			OP[2] += 1. / (double)numIter;
// 
// 			SaveVec3VecAsScatterZone({ tri[0], tri[1], tri[2], tri[0] }, string("triangle"));
// 			TecUtilZoneSetScatter(SV_SHOW, tecplot::toolbox::Set(TecUtilDataSetGetNumZones()).getRef(), 0, FALSE);
// 			TecUtilZoneSetMesh(SV_SHOW, tecplot::toolbox::Set(TecUtilDataSetGetNumZones()).getRef(), 0, TRUE);
// 			TecUtilZoneSetMesh(SV_LINETHICKNESS, tecplot::toolbox::Set(TecUtilDataSetGetNumZones()).getRef(), 0.2, 0);
// 			TecUtilZoneSetMesh(SV_COLOR, tecplot::toolbox::Set(TecUtilDataSetGetNumZones()).getRef(), 0, Red_C);
// 			TecUtilZoneSetActive(tecplot::toolbox::Set(TecUtilDataSetGetNumZones()).getRef(), AssignOp_MinusEquals);
// 
// 			double testProj = PointDistanceToTriangleSquared(OP, TP, tri[0], tri[1], tri[2]);
// 			SaveVec3VecAsScatterZone({ OP }, string("old point"), Red_C);
// 			TecUtilZoneSetScatter(SV_FRAMESIZE, tecplot::toolbox::Set(TecUtilDataSetGetNumZones()).getRef(), 1, 0);
// 			TecUtilZoneSetActive(tecplot::toolbox::Set(TecUtilDataSetGetNumZones()).getRef(), AssignOp_MinusEquals);
// 			SaveVec3VecAsScatterZone({ TP }, string("new point"), Green_C);
// 			TecUtilZoneSetScatter(SV_FRAMESIZE, tecplot::toolbox::Set(TecUtilDataSetGetNumZones()).getRef(), 1, 0);
// 			TecUtilZoneSetActive(tecplot::toolbox::Set(TecUtilDataSetGetNumZones()).getRef(), AssignOp_MinusEquals);
// 
// 			SaveVec3VecAsScatterZone({ OP, TP }, string("old to new point"), Green_C);
// 			TecUtilZoneSetScatter(SV_SHOW, tecplot::toolbox::Set(TecUtilDataSetGetNumZones()).getRef(), 0, FALSE);
// 			TecUtilZoneSetMesh(SV_SHOW, tecplot::toolbox::Set(TecUtilDataSetGetNumZones()).getRef(), 0, TRUE);
// 			TecUtilZoneSetMesh(SV_LINETHICKNESS, tecplot::toolbox::Set(TecUtilDataSetGetNumZones()).getRef(), 0.2, 0);
// 			TecUtilZoneSetMesh(SV_COLOR, tecplot::toolbox::Set(TecUtilDataSetGetNumZones()).getRef(), 0, Blue_C);
// 			TecUtilZoneSetActive(tecplot::toolbox::Set(TecUtilDataSetGetNumZones()).getRef(), AssignOp_MinusEquals);
// 		}
// 	}
// }
// 
// return;

	/*
	 *	Test surface grad path
	 */
// 	int SurfZoneNum = 50;
// 	int RhoVarNum = 4;
// 	VolExtentIndexWeights_s VolInfo;
// 	GetVolInfo(1, { 1,2,3 }, FALSE, VolInfo);
// 	FieldDataPointer_c RhoPtr;
// 	GetReadPtrsForZone(1, RhoVarNum, vector<int>(), vector<int>(), RhoPtr, vector<FieldDataPointer_c>(), vector<FieldDataPointer_c>());
// 
// 	FESurface_c Surf(SurfZoneNum, { 1,2,3 });
// // 	Surf.GenerateElemConnectivity();
// // 	Surf.GenerateElemMidpoints();
// 	Surf.GeneratePointElementDistanceCheckData();
// 	double TermVal = 0.001;
// 
// 	GradPath_c SurfGP({ -0.649,-0.672,-0.166 }, StreamDir_Reverse, 1000, GPType_Classic, GPTerminate_AtRhoValue, NULL, NULL, NULL, &TermVal, VolInfo, vector<FieldDataPointer_c>(), vector<FieldDataPointer_c>(), RhoPtr, &Surf);
// 
// 	SurfGP.Seed(false);
// 	if (SurfGP.IsMade()) SurfGP.SaveAsOrderedZone("SurfGradPath", Green_C);
// 
// 	TecUtilLockFinish(AddOnID);
// 	return;

	/*
	 *	test surface-sphere intersection function for file Z:\Tecplot\StorageWorkspace\2018_11_GBS-sphere-triangulationupdate\Al43.lpk
	 */
	int SphereZoneNum = 364;
	int CPZoneNum = 3;

	vec3 SphereCenter;
	for (int i = 0; i < 3; ++i){
		SphereCenter[i] = TecUtilDataValueGetByZoneVar(CPZoneNum, i + 1, 2);
	}

	double SphereRadius = stod(AuxDataZoneGetItem(SphereZoneNum, CSMAuxData.GBA.SphereRadius));

	FESurface_c Sphere(SphereZoneNum, { 1,2,3 });
	vector<vec3> const * SphereXYZPtr = Sphere.GetXYZListPtr();
	vector<vector<int> > const * SphereElemListPtr = Sphere.GetElemListPtr();

	/*
	 *	Get average sphere triangulation edge length so that we can 
	 *	resample intersection paths at roughly the same spacing.
	 *	I'll do this by looping over triangles, which will double count every 
	 *	edge so the average should still be valid.
	 */
	double SphereEdgeLenMean = 0.0;
	int SphereNumEdges = 0;
	for (auto const & t : *SphereElemListPtr){
		for (int i = 0; i < 3; ++i){
			SphereEdgeLenMean += Distance(SphereXYZPtr->at(t[i]), SphereXYZPtr->at(t[(i + 1) % 3]));
			SphereNumEdges++;
		}
	}
	SphereEdgeLenMean /= (double)SphereNumEdges;

	vector<vector<vec3> > AllIntPoints;
	vector<int> IntZoneNums;
	vector<FESurface_c> IntSurfs;

	TecUtilPleaseWait("finding ring surface intersections", TRUE);

	for (int z = 1; z <= TecUtilDataSetGetNumZones(); ++z){
		if (AuxDataZoneItemMatches(z, CSMAuxData.CC.ZoneSubType, CSMAuxData.CC.ZoneSubTypeRSSegment)){
			FESurface_c Surf(z, { 1,2,3 });
			vector<vec3> IntPoints = Surf.GetSphereIntersectionPath(SphereCenter, SphereRadius);
			if (IntPoints.size() > 1) {
				vector<vec3> ResampledIntPoints;
				double PathLen = 0;
				for (int i = 0; i < IntPoints.size() - 1; ++i)
					PathLen += Distance(IntPoints[i], IntPoints[i + 1]);
				if (Vec3PathResample(IntPoints, MAX(int(PathLen / SphereEdgeLenMean) + 1, 4), ResampledIntPoints)){
					IntPoints = ResampledIntPoints;
// 					Surf.GenerateElemConnectivity();
// 					Surf.GenerateElemMidpoints();
					Surf.GeneratePointElementDistanceCheckData();
					IntSurfs.push_back(Surf);

					/* 
					 * Now project each point back to the sphere radius
					 */
					for (auto & p : IntPoints){
						p = SphereCenter + normalise(p - SphereCenter) * SphereRadius;
					}

// 					SaveVec3VecAsScatterZone(IntPoints, string("Intersection of zone " + to_string(z)));
					AllIntPoints.push_back(IntPoints);
					IntZoneNums.push_back(z);
				}
			}
// 			break;
		}
	}

	/*
	 *	Check for nearly equal intersection path segment endpoints,
	 *	i.e. where two paths intersect.
	 *	These intersections should only occur at a bond-path-sphere intersection point,
	 *	at which there is already a constrained sphere node, so when path intersections are
	 *	found, move all the coincident points to the closest node on the sphere.
	 */
	double CoincidentCheckEpsilon = SphereEdgeLenMean * 0.1;
	CoincidentCheckEpsilon *= CoincidentCheckEpsilon;
	for (int i = 0; i < AllIntPoints.size() - 1; ++i) {
		vector<std::pair<int, int> > CoincidentPointIndices;
		for (int j = i+1; j < AllIntPoints.size(); ++j){
			if (DistSqr(AllIntPoints[i][0], AllIntPoints[j][0]) <= CoincidentCheckEpsilon){
				CoincidentPointIndices.push_back(std::pair<int, int>(i, 0));
				CoincidentPointIndices.push_back(std::pair<int, int>(j, 0));
			}
			else if (DistSqr(AllIntPoints[i][0], AllIntPoints[j].back()) <= CoincidentCheckEpsilon) {
				CoincidentPointIndices.push_back(std::pair<int, int>(i, 0));
				CoincidentPointIndices.push_back(std::pair<int, int>(j, AllIntPoints[j].size() - 1));
			}
			else if (DistSqr(AllIntPoints[i].back(), AllIntPoints[j][0]) <= CoincidentCheckEpsilon) {
				CoincidentPointIndices.push_back(std::pair<int, int>(i, AllIntPoints[i].size() - 1));
				CoincidentPointIndices.push_back(std::pair<int, int>(j, 0));
			}
			else if (DistSqr(AllIntPoints[i].back(), AllIntPoints[j].back()) <= CoincidentCheckEpsilon) {
				CoincidentPointIndices.push_back(std::pair<int, int>(i, AllIntPoints[i].size() - 1));
				CoincidentPointIndices.push_back(std::pair<int, int>(j, AllIntPoints[j].size() - 1));
			}
		}

		if (CoincidentPointIndices.size() > 0){
			int SphereNodeNum = -1;
			for (auto const & p : CoincidentPointIndices){
				double MinNodeDistSqr = DBL_MAX;
				double TmpNodeDistSqr;
				int MinNodeIndex = -1;
				for (int ni = 0; ni < SphereXYZPtr->size(); ++ni){
					TmpNodeDistSqr = DistSqr(AllIntPoints[p.first][p.second], SphereXYZPtr->at(ni));
					if (TmpNodeDistSqr < MinNodeDistSqr){
						MinNodeDistSqr = TmpNodeDistSqr;
						MinNodeIndex = ni;
					}
				}
				if (MinNodeIndex >= 0){
// 					if (SphereNodeNum >= 0){
// 						if (SphereNodeNum != MinNodeIndex) {
// 							TecUtilDialogErrMsg("Coincident intersection path endpoints do not share same closest sphere node!");
// 						}
// 					}
// 					else{
// 						SphereNodeNum = MinNodeIndex;
// 					}
					AllIntPoints[p.first][p.second] = SphereXYZPtr->at(MinNodeIndex);
				}
			}
		}

		SaveVec3VecAsScatterZone(AllIntPoints[i], string("Intersection of zone " + to_string(IntZoneNums[i])));
	}
	SaveVec3VecAsScatterZone(AllIntPoints.back(), string("Intersection of zone " + to_string(IntZoneNums.back())));

// 	return;

	vector<tpcsm::Vec3> initialNodeXYZs(SphereXYZPtr->size());
	for (int i = 0; i < initialNodeXYZs.size(); ++i)
		initialNodeXYZs[i] = tpcsm::Vec3(SphereXYZPtr->at(i)[0], SphereXYZPtr->at(i)[1], SphereXYZPtr->at(i)[2]);

	vector<TriNodes> initialTriangleNodes(SphereElemListPtr->size());
	for (int i = 0; i < initialTriangleNodes.size(); ++i)
		initialTriangleNodes[i] = TriNodes(SphereElemListPtr->at(i)[0], SphereElemListPtr->at(i)[1], SphereElemListPtr->at(i)[2]);

	tpcsm::Vec3 SphereCenterVec(SphereCenter[0], SphereCenter[1], SphereCenter[2]);

	vector<tpcsm::Vec3> constraintNodes;
	vector<Edge> constraintSegments;
	vector<int> constraintSegmentSurfNums;
	for (int j = 0; j < AllIntPoints.size(); ++j) {
		vector<vec3> IntPoints = AllIntPoints[j];
		constraintNodes.push_back(tpcsm::Vec3(IntPoints[0][0], IntPoints[0][1], IntPoints[0][2]));
		for (int i = 1; i < IntPoints.size(); ++i) {
			constraintNodes.push_back(tpcsm::Vec3(IntPoints[i][0], IntPoints[i][1], IntPoints[i][2]));
			constraintSegments.push_back(Edge(constraintNodes.size() - 2, constraintNodes.size() - 1));
			constraintSegmentSurfNums.push_back(j);
		}
	}

	/*
	 *	Need to remove duplicate constraint nodes and update edge connectivity to reflect the change.
	 */
	vector<int> nodeIndices(constraintNodes.size());
	for (int i = 0; i < nodeIndices.size(); ++i)
		nodeIndices[i] = i;
	vector<tpcsm::Vec3> tmpNodes;
	tmpNodes.reserve(constraintNodes.size());
	for (int ni = 0; ni < constraintNodes.size() - 1; ++ni){
		bool isFound = false;
		for (int nj = ni + 1; nj < constraintNodes.size(); ++nj){
			if (constraintNodes[nodeIndices[ni]] == constraintNodes[nodeIndices[nj]]){
				/*
				 *	duplicate node found, so search through edges, changing nj to ni when found.
				 */
				for (auto & e : constraintSegments){
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
	for (int ni = 0; ni < constraintNodes.size(); ++ni){
		if (nodeIndices[ni] == ni){
			/*
			 *	Node is not duplicate
			 */
			tmpNodes.push_back(constraintNodes[ni]);
			for (auto & e : constraintSegments){
				if (e.first == ni)
					e.first = tmpNodes.size() - 1;

				if (e.second == ni)
					e.second = tmpNodes.size() - 1;
			}
		}
	}
	constraintNodes = tmpNodes;

	/*
	 *	Now move nodes on the sphere that are sufficiently close to a constraint node
	 *	to eliminate the ugly triangles that result when two points are very close.
	 */
	for (auto & n : initialNodeXYZs){
		for (auto & p : constraintNodes){
			if ((p - n).getNormSquared() <= CoincidentCheckEpsilon){
				n = p;
				break;
			}
		}
	}


	/*
	 *	Update triangulation to include ring surface intersections as edges
	 */
	TecUtilPleaseWait("updating triangulation", TRUE);

	vector<tpcsm::Vec3> OutXYZs;
	vector<Index_t> OutConstrinedSegmentIndices;
	vector<TriNodes> OutTris;
	const char* StatusMessage;

	if (tpcsm::updateSphericalTriangulation(initialNodeXYZs,
		initialTriangleNodes,
		SphereCenterVec,
		SphereRadius,
		constraintNodes,
		constraintSegments,
		false,
		OutXYZs,
		OutConstrinedSegmentIndices,
		OutTris,
		StatusMessage)){

		vector<vec3> NewNodes(OutXYZs.size());
		for (int i = 0; i < NewNodes.size(); ++i)
			NewNodes[i] << OutXYZs[i].x() << OutXYZs[i].y() << OutXYZs[i].z();

		vector<vector<int> > NewElems(OutTris.size());
		for (int i = 0; i < OutTris.size(); ++i)
			NewElems[i] = { OutTris[i].v1(), OutTris[i].v2(), OutTris[i].v3() };

		FESurface_c NewSphere(NewNodes, NewElems);

		if (NewSphere.IsMade()){
			NewSphere.SaveAsTriFEZone({ 1,2,3 }, "UpdatedSphere");
		}

		/*
		 *	seed surface gradient paths for all the nodes that 
		 *	fall on constraint edges.
		 */

		TecUtilPleaseWait("seeding surface GPs on ring surfaces", TRUE);

		VolExtentIndexWeights_s VolInfo;
  		GetVolInfo(1, { 1,2,3 }, FALSE, VolInfo);
  		FieldDataPointer_c RhoPtr;
  		GetReadPtrsForZone(1, 4, vector<int>(), vector<int>(), RhoPtr, vector<FieldDataPointer_c>(), vector<FieldDataPointer_c>());

		for (int i = 0; i < OutConstrinedSegmentIndices.size(); ++i){
			if (OutConstrinedSegmentIndices[i] >= 0) {
				int ConstrainedSegNum = i,
					ConstrainedSegSurfNum = constraintSegmentSurfNums[OutConstrinedSegmentIndices[i]];
				if (ConstrainedSegSurfNum < IntSurfs.size()){
					FESurface_c * SurfPtr = &IntSurfs[ConstrainedSegSurfNum];
					GradPath_c GP(NewNodes[i], StreamDir_Both, 500, GPType_Classic, GPTerminate_AtBoundary, NULL, NULL, NULL, NULL, VolInfo, vector<FieldDataPointer_c>(), vector<FieldDataPointer_c>(), RhoPtr, SurfPtr);
					GP.Seed();
					if (GP.IsMade()){
						GP.SaveAsOrderedZone(string("SurfGP node " + to_string(i) + " cSegNum " + to_string(i) + " surfNum " + to_string(ConstrainedSegSurfNum) + " surfZoneNum " + to_string(SurfPtr->GetZoneNum())), Green_C);
					}
				}
			}
		}
	}

	TecUtilPleaseWait("finding ring surface intersections", FALSE);
	

	TecUtilLockFinish(AddOnID);
}

struct GradientPathToolData_s{
	VolExtentIndexWeights_s VolInfo;
	FieldDataPointer_c RhoPtr;
	vector<FieldDataPointer_c> GradPtrs, HessPtrs;
	StreamDir_e Dir;

// 	GradientPathToolData_s(VolExtentIndexWeights_s &VolInfo_in,
// 		FieldVecPointer_c & RhoPtr_in,
// 		vector<FieldVecPointer_c> & GradPtrs_in, 
// 		vector<FieldVecPointer_c> & HessPtrs_in)
// 	{
// 		VolInfo = VolInfo_in;
// 		RhoPtr = RhoPtr_in;
// 		GradPtrs = GradPtrs_in;
// 		HessPtrs = HessPtrs_in;
// 	}
};

GradientPathToolData_s GPToolData;

void STDCALL GradientPathToolProbeCB(Boolean_t WasSuccessful,
	Boolean_t isNearestPoint,
	ArbParam_t ClientData){

// 	EntIndex_t ZoneNum = TecUtilProbeFieldGetZone();
// 	if (ZoneNum == TECUTILBADZONENUMBER){
// 		TecUtilDialogErrMsg("Failed to get zone number!");
// 		return;
// 	}

// 	EntIndex_t ind = -1;
// 	if (isNearestPoint) ind = TecUtilProbeGetPointIndex();

// 	GradientPathToolData_s *GPToolData = reinterpret_cast<GradientPathToolData_s*>(ClientData);

	vec3 Pt;
	for (int i = 0; i < 3; ++i) Pt[i] = TecUtilProbeFieldGetValue(i+1);

	TecUtilLockStart(AddOnID);

	CSMGuiLock();

	
	double TermValue = 0;
	GradPath_c GP(Pt, GPToolData.Dir, 1000, GPType_Classic, GPTerminate_AtRhoValue, NULL, NULL, NULL, &TermValue, GPToolData.VolInfo, GPToolData.HessPtrs, GPToolData.GradPtrs, GPToolData.RhoPtr);
	GP.Seed(false);
	GP.SaveAsOrderedZone();

	TecUtilRedraw(TRUE);


	CSMGuiUnlock();


	TecUtilLockFinish(AddOnID);
}

void GradientPathToolReturnUserInfo(bool const GuiSuccess,
	vector<GuiField_c> const & Fields,
	vector<GuiField_c> const PassthroughFields){
	if (!GuiSuccess) return;

	TecUtilLockStart(AddOnID);

	int fNum = 0;
	int VolZoneNum = Fields[fNum++].GetReturnInt();
	StreamDir_e StreamDir = StreamDir_e(Fields[fNum++].GetReturnInt()-1);
	int RhoVarNum = Fields[fNum++].GetReturnInt();
	GPToolData.Dir = StreamDir;
	GetVolInfo(VolZoneNum, { 1, 2, 3 }, FALSE, GPToolData.VolInfo);

	Boolean_t IsOk = GPToolData.RhoPtr.GetReadPtr(VolZoneNum, RhoVarNum);

	if (IsOk && Fields[fNum++].GetReturnBool()){
		GPToolData.GradPtrs.resize(3);
		for (int i = 0; i < 3; ++i){
			if (!GPToolData.GradPtrs[i].GetReadPtr(VolZoneNum, Fields[fNum].GetReturnInt() + i)){
				GPToolData.GradPtrs.clear();
				break;
			}
		}
	}
	else GPToolData.GradPtrs.clear();
	fNum++;
	if (IsOk && Fields[fNum++].GetReturnBool()){
		GPToolData.HessPtrs.resize(6);
		for (int i = 0; i < 6; ++i){
			if (!GPToolData.HessPtrs[i].GetReadPtr(VolZoneNum, Fields[fNum].GetReturnInt() + i)){
				GPToolData.HessPtrs.clear();
				break;
			}
		}
	}
	else GPToolData.HessPtrs.clear();

	ArgList ProbeArgs;
	ProbeArgs.appendFunction(SV_CALLBACKFUNCTION, GradientPathToolProbeCB);
	ProbeArgs.appendString(SV_STATUSLINETEXT, "Select point at which to seed gradient path");
// 	ProbeArgs.appendArbParam(SV_CLIENTDATA, reinterpret_cast<ArbParam_t>(&GPToolData));

	TecUtilProbeInstallCallbackX(ProbeArgs.getRef());
	TecUtilProbeAllowCOBs();

	TecUtilLockFinish(AddOnID);
}

void GradientPathToolGetUserInfo(){
	vector<GuiField_c> Fields = {
		GuiField_c(Gui_ZoneSelect, "Full volume zone", CSMZoneName.FullVolume),
		GuiField_c(Gui_Radio, "Gradient path direction", "Forward,Reverse,Both"),
		GuiField_c(Gui_VarSelect, "Electron density variable number", CSMVarName.Dens)
	};

	int iTmp = Fields.size();
	Fields.push_back(GuiField_c(Gui_ToggleEnable, "Density gradient vector variables present"));
	Fields[iTmp].AppendSearchString(to_string(Fields.size()));
	Fields.push_back(GuiField_c(Gui_VarSelect, "X Gradient variable number", CSMVarName.DensGradX));
	iTmp = Fields.size();
	Fields.push_back(GuiField_c(Gui_ToggleEnable, "Density Hessian tensor variables present"));
	Fields[iTmp].AppendSearchString(to_string(Fields.size()));
	Fields.push_back(GuiField_c(Gui_VarSelect, "XX Hessian variable number", CSMVarName.DensHessXX));

	CSMGui("Gradient path tool", Fields, GradientPathToolReturnUserInfo, AddOnID);
}