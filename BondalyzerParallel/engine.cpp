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
#include <map>
#include <set>
#include <list>
#include <deque>
#include <queue>
#include <io.h>

#include <omp.h>


#include <chrono>
#include <ctime>
#include <ratio>

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

#include <gsl/gsl_bspline.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_rng.h>

#include <gsl/gsl_deriv.h>

#include "CSM_DATA_TYPES.h"
#include "CSM_DATA_LOADER_FUNCTIONS.h"
#include "CSM_DATA_SET_INFO.h"
#include "CSM_CRIT_POINTS.h"
#include "CSM_CALC_VARS.h"
#include "CSM_GRAD_PATH.h"
#include "CSM_FE_VOLUME.h"
#include "CSM_GUI.h"
#include "CSM_GEOMETRY.h"

#include "updateSphericalTriangulation.h"
#include "meshgen2d_sphere.h"

#include "KFc.h"

#include "ENGINE.h"

#include <armadillo>
#include "gsl/gsl_vector_double.h"
using namespace arma;
using namespace tecplot::toolbox;


#define DEBUG_PRINTPATHS
// #define DEBUG_SAVESCREENSHOTS

using std::string;
using std::to_string;
using std::stoi;
using std::vector;
using std::stringstream;
using std::ofstream;
using std::ios;
using std::endl;
using std::setprecision;

#define PI2 PI * 2.
#define ONEOVERPI 1. / PI


//for profiling
using std::chrono::high_resolution_clock;
using std::chrono::duration;
using std::chrono::duration_cast;

struct MinFuncParams_GPLengthBetweenPoints {
	GradPath_c * GP = nullptr;

	std::pair<vec3 const *,vec3 const *> * StartPoints = nullptr;
};

struct MinFuncParams_GPLengthInPlane{
	StreamDir_e Direction = StreamDir_Invalid;
	int NumGPPoints = -1;
	GPType_e GPType = GPType_Invalid;
	GPTerminate_e GPTerminate = GPTerminate_Invalid;
	vec3 * TermPoint = nullptr;
	CritPoints_c * CPs = nullptr;
	double * TermPointRadius = nullptr;
	double * TermValue = nullptr;
	VolExtentIndexWeights_s * VolInfo = nullptr;
	vector<FieldDataPointer_c> const * HessPtrs = nullptr;
	vector<FieldDataPointer_c> const * GradPtrs = nullptr;
	FieldDataPointer_c const * RhoPtr = nullptr;

	GradPath_c GP;

	vec3 StartPointOrigin;
	vec3 RotVec;
	vec3 RotAxis;
	int StartCPNum = -1;
	int EndCPPosition = -1;
	int EndCPNum = -1;
};

struct MinFuncParams_GPLengthInPlane2 {
	GradPath_c * GP;

	vec3 const * Origin;
	vec3 const * RotVec;
	vec3 const * RotAxis;
};

struct MinFuncParams_GPLengthSpherical{
	StreamDir_e Direction = StreamDir_Invalid;
	int NumGPPoints = -1;
	GPType_e GPType = GPType_Invalid;
	GPTerminate_e GPTerminate = GPTerminate_Invalid;
	vec3 * TermPoint = nullptr;
	CritPoints_c * CPs = nullptr;
	double * TermPointRadius = nullptr;
	double * TermValue = nullptr;
	VolExtentIndexWeights_s * VolInfo = nullptr;
	vector<FieldDataPointer_c> const * HessPtrs = nullptr;
	vector<FieldDataPointer_c> const * GradPtrs = nullptr;
	FieldDataPointer_c const * RhoPtr = nullptr;

	GradPath_c GP;

	vec3 StartPointOrigin;
	double r;
	int StartCPNum = -1;
	int EndCPNum = -1;
};


string const DS_ZoneName_Delim = "_-_";
string const T41Ext = ".t41";

static int const MinNumCircleCheckGPs = 180;
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
static int const DefaultSaddleGPLengthCheckNum = 15; // Number of GPs used to find a minimum length GP to a saddle CP between two deviating GPs
static double const DefaultSaddleGPLengthCheckWeightRange = 0.9;

//DEBUG
#ifdef _DEBUG
vector<vector<double> > MinFunc_GPAngleLengths;
#endif

/*
 *	This is the convergence criteria for a one-dimensional GSL minimization
 *	based on gradient path length.
 */
static double const Tolerance_GPLengthInPlane = 1e-8;

/*
 *	If the terminal ends of two gradient paths diverge away from eachother
 *	such that the inner (dot) product of their terminal ends is less than
 *	below, they are considered divergent and to be terminating at different 
 *	far field CPs (if they were near field, then we would know which CP and
 *	don't need to resort to checking for divergence).
 */
static double const Tolerance_DivergentGPInnerProd = -0.9;

BondalyzerCalcType_e CurrentCalcType;

enum GBATriangulatedSphereNodeElemType_e
{
	NETypeInvalid = -1,
	NETypeC, // goes to cage cp, system boundary, or rho cutoff
	NETypeR, // goes along a ring surface to a ring cp, then to a cage cp or rho cutoff
	NETypeB, // goes along bond path
	NETypeRB // B type on same elem as R type
};


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

	CSMGUILock();

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

	CSMGUIUnlock();
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

void FEZoneBFS(int NodeNum,
	NodeMap_pa const & NodeMap,
	NodeToElemMap_pa const & NodeToElemMap,
	int NumNodesPerElem,
	vector<bool> & IsVisited)
{
	std::queue<int> q;
	q.push(NodeNum);
	IsVisited[NodeNum] = true;

	while (!q.empty()){
		NodeNum = q.front();
		q.pop();
		int NumElems = TecUtilDataNodeToElemMapGetNumElems(NodeToElemMap, NodeNum + 1);
		for (int e = 1; e <= NumElems; ++e) {
			int ei = TecUtilDataNodeToElemMapGetElem(NodeToElemMap, NodeNum + 1, e);
			for (int n = 0; n < NumNodesPerElem; ++n) {
				int ni = TecUtilDataNodeGetByRef(NodeMap, ei, n + 1) - 1;
				if (!IsVisited[ni]) {
					q.push(ni);
					IsVisited[ni] = true;
				}
			}
		}
	}
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
		if (!IsoReadPtrs[i].InitializeReadPtr(IsoZoneNum, i + 1)){
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
		if (!IsoReadPtrs[i].InitializeReadPtr(IsoZoneNum, i + 1)){
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
		if (!IsoReadPtrs[i].InitializeReadPtr(IsoZoneNum, i + 1)){
			TecUtilDialogErrMsg("Failed to get isosurface read pointer. Quitting.");
			return;
		}
	}

	GetClosedIsoSurface(IsoZoneNum, IsoReadPtrs, NodeNums);

	TecUtilDataLoadEnd();

	TecUtilLockFinish(AddOnID);
}

void GetClosedIsoSurfaceProbeCallback(Boolean_t WasSuccessful, Boolean_t IsNearestPoint, ArbParam_t ClientData){
	TecUtilLockStart(AddOnID);
	if (WasSuccessful){
		// first check that zone is isosurface
		int ZoneNum = TecUtilProbeFieldGetZone();
		if (ZoneNum > 0 && TecUtilZoneIsFiniteElement(ZoneNum)) {
			TecUtilDataLoadBegin();
			vector<FieldDataPointer_c> IsoReadPtrs(TecUtilDataSetGetNumVars());
			for (int i = 0; i < IsoReadPtrs.size(); ++i) {
				if (!IsoReadPtrs[i].InitializeReadPtr(ZoneNum, i + 1)) {
					TecUtilDialogErrMsg("Failed to get isosurface read pointer. Quitting.");
					return;
				}
			}

			int NodeNum = -1;
			if (IsNearestPoint){
				NodeNum = TecUtilProbeGetPointIndex();
			}
			else{
				// get nearest node to point.
				vector<int> XYZVarNums = { 1,2,3 };
				vec3 ProbePoint;
				for (int i = 0; i < 3; ++i)
					ProbePoint[i] = TecUtilProbeFieldGetValue(XYZVarNums[i]);
				FieldVecPointer_c XYZPtr;
				XYZPtr.InitializeReadPtr(ZoneNum, XYZVarNums);
				double MinDistSqr = DBL_MAX;
				int MinDistNode = -1;
				for (int i = 0; i < IsoReadPtrs[0].Size(); ++i){
					double TmpDistSqr = DistSqr(ProbePoint, XYZPtr[i]);
					if (TmpDistSqr < MinDistSqr){
						MinDistSqr = TmpDistSqr;
						MinDistNode = i;
					}
				}
				if (MinDistNode < 0){
					TecUtilDialogErrMsg("Failed to find minimum distance node to selected point");
				}
				else{
					NodeNum = MinDistNode;
				}
			}
			
			vector<int> NodeNums = { NodeNum };

			GetClosedIsoSurface(ZoneNum, IsoReadPtrs, NodeNums);

			TecUtilDataLoadEnd();
		}
	}
	TecUtilLockFinish(AddOnID);
}

void GetClosedIsoSurfacesInstallProbeCallback() {
	ArgList args;
	args.appendFunction(SV_CALLBACKFUNCTION, GetClosedIsoSurfaceProbeCallback);
	args.appendString(SV_STATUSLINETEXT, "Select extracted isosurface component to extract");
	TecUtilProbeInstallCallbackX(args.getRef());
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

// 		FEZoneDFS(NodeNum, NodeMap, NodeToElemMap, IsoIJK[2], NodeVisited);
		FEZoneBFS(NodeNum, NodeMap, NodeToElemMap, IsoIJK[2], NodeVisited);

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

		for (int e = 0; e < IsoIJK[1]; ++e) {
			if (!ElemAdded[e]) {
				for (auto const & n : Elems[e]) {
					if (NodeVisited[n])
					{
						ElemAdded[e] = true;
						NewElems.push_back(Elems[e]);
						for (auto & ne : NewElems.back()) {
							ne = NodeNumsOldToNew[ne];
						}
						break;
					}
				}
			}
		}

		/*
		 * Make new zone with single connected component
		 */

		if (!TecUtilDataSetAddZone(string(IsoZoneName + string(": Subzone ") + to_string(SubZoneNum++)).c_str(), NumNewNodes, NewElems.size(), IsoIJK[2], IsoZoneType, nullptr)){
			TecUtilDialogErrMsg("Failed to make new iso zone. Quitting.");
			return;
		}

		int NewZoneNum = TecUtilDataSetGetNumZones();

		for (int i = 0; i < IsoReadPtrs.size(); ++i){
			if (!IsoWritePtrs[i].InitializeWritePtr(NewZoneNum, i + 1)){
				TecUtilDialogErrMsg("Failed to get write pointer(s) for new iso zone. Quitting.");
				return;
			}
		}

		vector<int> NewZoneIJK(3);
		TecUtilZoneGetIJK(NewZoneNum, &NewZoneIJK[0], &NewZoneIJK[1], &NewZoneIJK[2]);

// 		for (int n = 0; n < NumNewNodes; ++n) for (int i = 0; i < IsoWritePtrs.size(); ++i) IsoWritePtrs[i].Write(n, IsoReadPtrs[i][NodeNumsNewToOld[n]]);
		for (int i = 0; i < IsoWritePtrs.size(); ++i) {
			if (IsoWritePtrs[i].FDType() != FieldDataType_Bit) {
				for (int n = 0; n < NumNewNodes; ++n) {
					IsoWritePtrs[i].Write(n, IsoReadPtrs[i][NodeNumsNewToOld[n]]);
				}
			}
		}

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

	bool PrecalcVars = Fields[fNum++].GetReturnBool();

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
		double CellSpacing = Fields[fNum++].GetReturnDouble();
		int ConvergedIterations = Fields[fNum++].GetReturnInt();
		bool DoDataAgitation = Fields[fNum++].GetReturnBool();
		double AgitationFactor = Fields[fNum++].GetReturnDouble();
		AgitationFactor = CLAMP(AgitationFactor, 0.0, 0.5);
		int AgitationNumIter = Fields[fNum++].GetReturnInt();
		AgitationNumIter = MAX(AgitationNumIter, 0);
		if (!DoDataAgitation){
			AgitationNumIter = 0;
		}
		FindCritPoints(VolZoneNum, XYZVarNums, RhoVarNum, GradVarNums, HessVarNums, IsPeriodic, CellSpacing, ConvergedIterations, PrecalcVars, AgitationFactor, AgitationNumIter);
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

		double RhoCutoff = Fields[fNum++].GetReturnDouble();
		int NumGPPts = Fields[fNum++].GetReturnInt();

		CPType_e CPType;
		if (CurrentCalcType == BondalyzerCalcType_BondPaths || CurrentCalcType == BondalyzerCalcType_InteratomicSurfaces) CPType = CPType_Bond;
		else if (CurrentCalcType == BondalyzerCalcType_RingLines || CurrentCalcType == BondalyzerCalcType_RingSurfaces) CPType = CPType_Ring;

		if (CurrentCalcType == BondalyzerCalcType_BondPaths || CurrentCalcType == BondalyzerCalcType_RingLines)
			FindBondRingLines(VolZoneNum, OtherCPZoneNums, SelectCPsZoneNum, SelectedCPs, CPTypeVarNum, CPType, XYZVarNums, RhoVarNum, GradVarNums, HessVarNums, IsPeriodic, PrecalcVars, NumGPPts);
		else if (CurrentCalcType == BondalyzerCalcType_CageNuclearPaths)
			FindCageNuclearPaths(VolZoneNum, OtherCPZoneNums, SelectCPsZoneNum, SelectedCPs, CPTypeVarNum, XYZVarNums, RhoVarNum, GradVarNums, HessVarNums, IsPeriodic, PrecalcVars);
		else if (CurrentCalcType == BondalyzerCalcType_InteratomicSurfaces || CurrentCalcType == BondalyzerCalcType_RingSurfaces) {
			// 			 FindBondRingSurfaces(VolZoneNum, OtherCPZoneNums, SelectCPsZoneNum, SelectedCPs, CPTypeVarNum, CPType, XYZVarNums, RhoVarNum, GradVarNums, HessVarNums, RidgeFuncVarNum, IsPeriodic);
			bool DebugMode = Fields[fNum].GetReturnBool();
			FindBondRingSurfaces2(VolZoneNum, OtherCPZoneNums, SelectCPsZoneNum, SelectedCPs, CPTypeVarNum, CPType, XYZVarNums, RhoVarNum, GradVarNums, HessVarNums, RidgeFuncVarNum, IsPeriodic, DebugMode, RhoCutoff, NumGPPts, PrecalcVars);
		}
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
		int ConvergedIterations = Fields[fNum++].GetReturnInt();
		bool DoDataAgitation = Fields[fNum++].GetReturnBool();
		double AgitationFactor = Fields[fNum++].GetReturnDouble();
		AgitationFactor = CLAMP(AgitationFactor, 0.0, 0.5);
		int AgitationNumIter = Fields[fNum++].GetReturnInt();
		AgitationNumIter = MAX(AgitationNumIter, 0);
		if (!DoDataAgitation) {
			AgitationNumIter = 0;
		}
		int RidgeFuncVarNum = 0;
		if (Fields[fNum++].GetReturnBool())
			RidgeFuncVarNum = Fields[fNum++].GetReturnInt();
		else fNum++;

		double RhoCutoff = Fields[fNum++].GetReturnDouble();
		int NumGPPts = Fields[fNum++].GetReturnInt();

		fNum++;
		vector<bool> DoCalcs;
		for (int i = (int)BondalyzerCalcType_BondPaths; i <= (int)BondalyzerCalcType_RingSurfaces; ++i) DoCalcs.push_back(Fields[fNum++].GetReturnBool());

		int NumZones = TecUtilDataSetGetNumZones();

		if (!UseCPZones && FindCritPoints(VolZoneNum, XYZVarNums, RhoVarNum, GradVarNums, HessVarNums, IsPeriodic, CellSpacing, ConvergedIterations, true, AgitationFactor, AgitationNumIter)){
			CPZoneNums.clear();
			for (int i = NumZones + 1; i <= TecUtilDataSetGetNumZones(); ++i) CPZoneNums.push_back(i);
			CPTypeVarNum = VarNumByName(CSMVarName.CritPointType);
		}

		BondalyzerBatch(VolZoneNum, CPZoneNums, CPTypeVarNum, XYZVarNums, RhoVarNum, GradVarNums, HessVarNums, IsPeriodic, RidgeFuncVarNum, DoCalcs, NumGPPts);
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

	Fields.push_back(GuiField_c(Gui_Toggle, "Precalculate missing variables", "1"));

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
		Fields.push_back(GuiField_c(Gui_Int, "Search converged iterations", to_string(1)));
		Fields.push_back(GuiField_c(Gui_Toggle, "Try to resolve spurious CPs using data agitation", "0"));
		Fields.push_back(GuiField_c(Gui_Double, "Spurious CPs: Data agitation factor (fraction of rho value):", to_string(0.0001)));
		Fields.push_back(GuiField_c(Gui_Int, "Spurious CPs: Number agitation cycles:", to_string(1000)));

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
		Fields.push_back(GuiField_c(Gui_Double, "GP rho cutoff", to_string(DefaultRhoCutoff)));
		Fields.push_back(GuiField_c(Gui_Int, "GP Number of points", "500"));
	}
	else if (CalcType == BondalyzerCalcType_Batch){
		int TogNum = Fields.size();
		Fields.push_back(GuiField_c(Gui_ToggleEnable, "Use existing critical points zone(s)?"));
		Fields[TogNum].AppendSearchString(to_string(Fields.size()));
		Fields.push_back(GuiField_c(Gui_ZoneSelectMulti, "Critical points", CSMZoneName.CriticalPoints));
		Fields[TogNum].AppendSearchString("," + to_string(Fields.size()));
		Fields.push_back(GuiField_c(Gui_VarSelect, "CP type variable", CSMVarName.CritPointType));
		Fields.push_back(GuiField_c(Gui_Double, "CP search grid spacing", to_string(DefaultCellSpacing)));
		Fields.push_back(GuiField_c(Gui_Int, "Search converged iterations", to_string(0)));
		Fields.push_back(GuiField_c(Gui_Toggle, "Try to resolve spurious CPs using data agitation", "0"));
		Fields.push_back(GuiField_c(Gui_Double, "Spurious CPs: Data agitation factor (fraction of rho value):", to_string(0.0001)));
		Fields.push_back(GuiField_c(Gui_Int, "Spurious CPs: Number agitation cycles:", to_string(1000)));
		Fields.push_back(GuiField_c(Gui_ToggleEnable, "1-RCS function present", to_string(Fields.size() + 1)));
		Fields.push_back(GuiField_c(Gui_VarSelect, "1-ridge function variable", CSMVarName.EberlyFunc[0]));
		Fields.push_back(GuiField_c(Gui_Double, "GP rho cutoff", to_string(DefaultRhoCutoff)));
		Fields.push_back(GuiField_c(Gui_Int, "GP Number of points", "500"));
		Fields.push_back(GuiField_c(Gui_Label, "Calculation steps"));
		for (int i = (int)BondalyzerCalcType_BondPaths; i <= (int)BondalyzerCalcType_RingSurfaces; ++i) Fields.push_back(GuiField_c(Gui_Toggle, BondalyzerStepGUITitles[i], "1"));
	}

	if (CalcType == BondalyzerCalcType_InteratomicSurfaces || CalcType == BondalyzerCalcType_RingSurfaces)
		Fields.push_back(GuiField_c(Gui_Toggle, "Interactive mode", "0"));

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
		VarReadPtrs[v].InitializeReadPtr(ZoneNum, v + 1);
		FDTypes[v] = VarReadPtrs[v].FDType();
	}

	char* ZoneName;
	TecUtilZoneGetName(ZoneNum, &ZoneName);
	int NewZoneNum = -1;

	if (IJK[0] - DelPointNums.size() > 0 && TecUtilDataSetAddZone(ZoneName, IJK[0] - DelPointNums.size(), IJK[1], IJK[2], ZoneType_Ordered, FDTypes.data())){
		NewZoneNum = TecUtilDataSetGetNumZones();
		for (int v = 0; v < NumVars; ++v){
			VarWritePtrs[v].InitializeWritePtr(NewZoneNum, v + 1);
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
	CSMGUILock();

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
	}


	TecUtilSetClear(ZoneNumsSet);
	TecUtilSetAddMember(ZoneNumsSet, ZoneNum, TRUE);
	for (int z : OtherCPZoneNums) {
		TecUtilSetAddMember(ZoneNumsSet, z, TRUE);
	}

	TecUtilDataSetDeleteZone(ZoneNumsSet);


	int MinZoneNumOld = TecUtilDataSetGetNumZones();
	for (int z : OldZoneNums) {
		MinZoneNumOld = MIN(MinZoneNumOld, z);
	}
	int MinZoneNumNew = TecUtilDataSetGetNumZones();
	for (int z : NewZoneNums) {
		MinZoneNumNew = MIN(MinZoneNumNew, z - TecUtilSetGetMemberCount(ZoneNumsSet));
	}
	SetZoneNum(MinZoneNumNew, MinZoneNumOld);

	TecUtilSetDealloc(&ZoneNumsSet);
	CSMGUIUnlock();
	TecUtilLockFinish(AddOnID);
}

void DeleteCPsGetUserInfo(){
	string OtherCPZoneSearchString = CSMZoneName.CriticalPoints;
	for (auto s : CSMZoneName.CPType){
		OtherCPZoneSearchString += "," + s;
	}
	vector<GuiField_c> Fields = {
		GuiField_c(Gui_ZonePointSelectMulti, "Critical point(s)", CSMZoneName.CriticalPoints),
		GuiField_c(Gui_VarSelect, "CP type", CSMVarName.CritPointType),
		GuiField_c(Gui_Toggle, "Also delete from other CP zone(s)", "1"),
		GuiField_c(Gui_ZoneSelectMulti, "Other CP zone(s)", OtherCPZoneSearchString)
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




Boolean_t FindCritPoints(int VolZoneNum,
								vector<int> const & XYZVarNums,
								int RhoVarNum,
								vector<int> & GradVarNums,
								vector<int> & HessVarNums,
								Boolean_t IsPeriodic,
								double const & CellSpacing,
								int ConvergedIterations,
								bool PrecalcVars,
	double AgitationFactor,
	int AgitationMaxNumIter)
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

	// Grad and Hess needed, so calculate them if they werent already present
	bool DeleteGradVars = false, DeleteHessVars = false;
	if (PrecalcVars && (GradVarNums.empty() || HessVarNums.empty())) {
		CalcVarsOptions_s opt;
		opt.AddOnID = AddOnID;
		opt.CalcForAllZones = FALSE;
		opt.CalcZoneNum = VolZoneNum;
		opt.RhoVarNum = RhoVarNum;
		opt.HasGrad = (!GradVarNums.empty());
		if (opt.HasGrad) {
			opt.GradVarNums = GradVarNums;
		}
		else {
			DeleteGradVars = true;
			opt.CalcVarList = { CalcGradientVectors };
		}
		DeleteHessVars = true;
		opt.CalcVarList.push_back(CalcHessian);
		opt.IsPeriodic = FALSE;

		CalcVars(opt);

		if (GradVarNums.empty()) {
			GradVarNums.push_back(VarNumByName("X Density", true));
			for (int i = 1; i < 3; ++i)
				GradVarNums.push_back(GradVarNums[0] + i);

			if (GradVarNums.empty())
				return FALSE;
		}

		if (HessVarNums.empty()) {
			HessVarNums.push_back(VarNumByName("XX Density", true));
			for (int i = 1; i < 6; ++i)
				HessVarNums.push_back(HessVarNums[0] + i);

			if (HessVarNums.empty())
				return FALSE;
		}
	}

	if (!GetReadPtrsForZone(VolZoneNum, 
		RhoVarNum, GradVarNums, HessVarNums,
		RhoPtr, GradPtrs, HessPtrs)){
		TecUtilDialogErrMsg("Failed to get read pointer(s)");
		return FALSE;
	}

	CSMGUILock();

	CritPoints_c VolCPs;

	StatusDrop(AddOnID);

	double RhoCutoff = DefaultRhoCutoff;

// 	if (FindCPs(VolCPs, VolInfo, CellSpacing, RhoCutoff, IsPeriodic, RhoPtr, GradPtrs, HessPtrs)){
// 		auto CPZoneNums = VolCPs.SaveAsOrderedZone(XYZVarNums, RhoVarNum, TRUE, VolZoneNum);
// 	}

	int NumCPs = 0;
	int Iter = 0;
	int NumMatch = 0;
	int CheckMatch = ConvergedIterations;
	int MaxIter = 5;
	
	Boolean_t UserQuit;

	double Spacing = CellSpacing;
	do 
	{
		CritPoints_c TmpCPs;
		UserQuit = !FindCPs(TmpCPs, VolInfo, Spacing, RhoCutoff, IsPeriodic, RhoPtr, GradPtrs, HessPtrs, AgitationFactor, AgitationMaxNumIter);
		if (UserQuit)
			break;

		VolCPs += TmpCPs;

		VolCPs.RemoveSpuriousCPs(Spacing);

		if (NumCPs == VolCPs.NumCPs())
			NumMatch++;
		else
			NumMatch = 0;

		NumCPs = VolCPs.NumCPs();

		Spacing *= 0.97;

	} while (NumMatch < CheckMatch && ++Iter < MaxIter);

	if (!UserQuit)
		auto CPZoneNums = VolCPs.SaveAsOrderedZone(XYZVarNums, RhoVarNum, TRUE, VolZoneNum);

	Set DelVars;
	if (DeleteGradVars) {
		for (int i : GradVarNums)
			DelVars += i;
	}

	if (DeleteHessVars) {
		for (int i : HessVarNums)
			DelVars += i;
	}

	if (!DelVars.isEmpty()) {
		TecUtilDataSetDeleteVar(DelVars.getRef());
	}

	TecUtilDataLoadEnd();

	CSMGUIUnlock();
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
	Boolean_t IsPeriodic,
	bool PrecalcVars,
	int NumGPPts)
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
		MinDistTypes = { CPType_Cage, CPType_Ring };
		GPDir = StreamDir_Reverse;
		PathColor = Green_C;
		EndCPNumforName++;
	}

	if (SelectedCPNums.size() == 0){
		for (int i = 0; i < AllCPs.NumCPs(TypeInd); ++i) SelectedCPNums.push_back(i+1);
	}

	CSMGUILock();

	vec3 StartPoint;

	double const StartPointOffset = 0.01 * AllCPs.GetMinCPDist(MinDistTypes);
// 	int const NumGPPts = 300;
	double TermRadius = 0.1 * AllCPs.GetMinCPDist(MinDistTypes);
	double RhoCutoff = DefaultRhoCutoff;

	CPTypeNum_e TermTypeNum = (CPType == CPType_Bond ? CPTypeNum_Nuclear : CPTypeNum_Cage);

	vector<GradPath_c> GPs;

	for (int iCP = 0; iCP < SelectedCPNums.size(); ++iCP){
		int CPInd = SelectedCPNums[iCP] - 1;
		for (int d = -1; d < 2; d += 2){
			StartPoint = AllCPs.GetXYZ(TypeInd, CPInd) + AllCPs.GetPrincDir(TypeInd, CPInd) * StartPointOffset * static_cast<double>(d);

			GPs.push_back(GradPath_c(StartPoint, GPDir, NumGPPts, GPType_Classic, GPTerminate_AtCP, nullptr, &AllCPs, &TermRadius, &RhoCutoff, VolInfo, HessPtrs, GradPtrs, RhoPtr));
			GPs.back().SetTerminalCPTypeNum(TermTypeNum);
			GPs.back().SetStartEndCPNum(AllCPs.GetTotOffsetFromTypeNumOffset(TypeInd, CPInd), 0);

		}
	}

	int NumGPs = GPs.size();

	vector<VolExtentIndexWeights_s> Th_VolInfo(omp_get_num_procs(), VolInfo);

#ifndef _DEBUG
#pragma omp parallel for schedule(dynamic)
#endif
	for (int iGP = 0; iGP < NumGPs; ++iGP){
		GPs[iGP].Seed(false);
		GPs[iGP].PointPrepend(AllCPs.GetXYZ(GPs[iGP].GetStartEndCPNum()[0]), AllCPs.GetRho(GPs[iGP].GetStartEndCPNum()[0]));
		GPs[iGP].Resample(NumGPPts);
		GPs[iGP].MakeRhoValuesMonotonic(&Th_VolInfo[omp_get_thread_num()], &RhoPtr);
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

		AuxDataZoneSetItem(GP.GetZoneNum(), CSMAuxData.GBA.SourceZoneNum, to_string(ZoneFinalSourceZoneNum(SelectedCPZoneNum, true)));
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

		AuxDataZoneSetItem(GP.GetZoneNum(), CSMAuxData.GBA.SourceZoneNum, to_string(ZoneFinalSourceZoneNum(SelectedCPZoneNum, true)));
	}

	TecUtilDataLoadEnd();


	CSMGUIUnlock();
	TecUtilLockFinish(AddOnID);

	return TRUE;
}

/*
 *	Create bond paths/ring lines from a source CP zone(s)
 */
Boolean_t DrawRepresentationQuadrics(int VolZoneNum,
	int SelectedCPZoneNum,
	vector<int> SelectedCPNums,
	int CPTypeVarNum,
	vector<int> const & XYZVarNums,
	int RhoVarNum,
	vector<int> const & GradVarNums,
	vector<int> const & HessVarNums,
	double SizeFactor)
{
	TecUtilLockStart(AddOnID);

	int NumZones = TecUtilDataSetGetNumZones();
	int NumVars = TecUtilDataSetGetNumVars();

	REQUIRE(XYZVarNums.size() == 3);
	for (auto const & i : XYZVarNums) REQUIRE(i > 0 && i <= NumVars);

	if (CPTypeVarNum < 0) {
		TecUtilDialogErrMsg("Failed to find CP Type variable");
		return FALSE;
	}

	FieldDataPointer_c RhoPtr;
	vector<FieldDataPointer_c> GradPtrs, HessPtrs;

	TecUtilDataLoadBegin();

	if (!GetReadPtrsForZone(VolZoneNum,
		RhoVarNum, GradVarNums, HessVarNums,
		RhoPtr, GradPtrs, HessPtrs)) {
		TecUtilDialogErrMsg("Failed to get read pointer(s)");
		return FALSE;
	}

	VolExtentIndexWeights_s VolInfo;
	VolInfo.AddOnID = AddOnID;

	if (!GetVolInfo(VolZoneNum, XYZVarNums, FALSE, VolInfo)) {
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

	vec3 EigVals;
	mat33 EigVecs;

	MultiRootParams_s MR;
	MR.CalcType = GPType_Classic;
	MR.VolInfo = &VolInfo;
	MR.IsPeriodic = FALSE;
	MR.HasGrad = (GradPtrs.size() == 3);
	MR.HasHess = (HessPtrs.size() == 6);
	MR.RhoPtr = &RhoPtr;
	MR.GradPtrs = &GradPtrs;
	MR.HessPtrs = &HessPtrs;
	MR.BasisVectors = &VolInfo.BasisNormalized;

	CritPoints_c AllCPs(SelectedCPZoneNum, XYZVarNums, CPTypeVarNum, RhoVarNum, &MR);

	if (SelectedCPNums.size() == 0) {
		for (int i = 0; i < AllCPs.NumCPs(); ++i) SelectedCPNums.push_back(i + 1);
	}

	CSMGUILock();

	int NumHyperboloidPts = 120;
	int SphereSubdivision = 3;

	double AngleStep = TWOPI / double(NumHyperboloidPts);

	Set NewZoneNums, DeleteZoneNums;

	for (int iCP = 0; iCP < SelectedCPNums.size(); ++iCP) {
		int CPInd = SelectedCPNums[iCP] - 1;
		
		auto CPPos = AllCPs.GetXYZ(CPInd);
		auto EigenVecs = AllCPs.GetEigVecs(CPInd);
		auto EigenVals = AllCPs.GetEigVals(CPInd);
		auto CPType = AllCPs.GetTypeFromTotOffset(CPInd);
		int TypeInd = VectorGetElementNum(CPTypeList, CPType);
		int TypeOffset = AllCPs.GetTypeNumOffsetFromTotOffset(CPInd)[1];
		auto DistToOtherCP = AllCPs.GetMinCPDist(CPInd);
		auto r = DistToOtherCP * SizeFactor;
		if (DistToOtherCP < 0.0){
			DistToOtherCP = 0.5;
		}

		stringstream ZoneName;
		FESurface_c Surf;
		vector<vec3> Nodes;
		vector<vector<int> > Elems;

		ZoneName << CPNameList[TypeInd] << " " << to_string(TypeOffset + 1) << " quadric: ";

		if (CPType == CPType_Bond || CPType == CPType_Ring) {
			vec3 RotAxis, v1, v2 = EigenVecs.col(1);;
			double th, phi;

			if (CPType == CPType_Bond) {
				RotAxis = EigenVecs.col(2);
				v1 = EigenVecs.col(0);
				th = atan(sqrt(abs(EigenVals[0] / EigenVals[2])));
				phi = atan(sqrt(abs(EigenVals[1] / EigenVals[2])));
			}
			else {
				RotAxis = EigenVecs.col(0);
				v1 = EigenVecs.col(2);
				th = atan(sqrt(abs(EigenVals[0] / EigenVals[2])));
				phi = atan(sqrt(abs(EigenVals[1] / EigenVals[2])));
			}

			ZoneName << setprecision(4) << "theta = " << th * 180. * ONEOVERPI << "; phi = " << phi * 180. * ONEOVERPI;

			Nodes.resize(2*NumHyperboloidPts+9);

			for (int i = 0; i < NumHyperboloidPts; ++i){
				double alpha = double(i) * AngleStep;
				vec3 NewNode = Rotate(v1, alpha, RotAxis);
				vec3 NormVec = Rotate(v2, alpha, RotAxis);
				double beta = th * pow(dot(NewNode, v1), 2) + phi * pow(dot(NewNode, v2), 2);
				NewNode = Rotate(NewNode, beta, NormVec) * r;
				Nodes[i] = CPPos + NewNode;
				Nodes[i + NumHyperboloidPts] = CPPos - NewNode;
			}

			// Now make tri surface
			
			for (int i = 0; i < NumHyperboloidPts; ++i){
				Elems.push_back(vector<int>({ i, (i + 1) % NumHyperboloidPts, 2 * NumHyperboloidPts }));
				Elems.push_back(vector<int>({ i + NumHyperboloidPts, (i + 1) % NumHyperboloidPts + NumHyperboloidPts, 2 * NumHyperboloidPts }));
			}

			int ii = 2 * NumHyperboloidPts, CPii = ii++;
			Nodes[CPii] = CPPos;

			vec3 NewNode = Rotate(v1, th, v2) * r;
			Nodes[ii++] = CPPos + NewNode;
			Nodes[ii++] = CPPos - NewNode;
			Elems.push_back(vector<int>({ ii - 1,ii - 2,CPii }));

			NewNode = Rotate(v1, -th, v2) * r;
			Nodes[ii++] = CPPos + NewNode;
			Nodes[ii++] = CPPos - NewNode;
			Elems.push_back(vector<int>({ ii - 1,ii - 2,CPii }));

			NewNode = Rotate(v2, phi, v1) * r;
			Nodes[ii++] = CPPos + NewNode;
			Nodes[ii++] = CPPos - NewNode;
			Elems.push_back(vector<int>({ ii - 1,ii - 2,CPii }));

			NewNode = Rotate(v2, -phi, v1) * r;
			Nodes[ii++] = CPPos + NewNode;
			Nodes[ii++] = CPPos - NewNode;
			Elems.push_back(vector<int>({ ii - 1,ii - 2,CPii }));
		}
		else {
			// Make triangulated sphere, then loop over each point, adjusting its distance from the center based on its alignment with eigenvectors
			int NumNodes, NumElems, NumEdges;
			point * meshP;
			triangle * meshT;
			int ** meshE;
			vector<int> MovedPointNums;
			vector<point> IntersectionPoints;
			auto MeshStatus = meshgen2D_sphere(1.0, SphereSubdivision, IntersectionPoints, MovedPointNums, meshP, meshT, meshE, NumNodes, NumElems, NumEdges);

			vec3 v1, v2 = EigenVecs.col(1) , v3;
			double r1 = r, r2, r3;
			if (CPType == CPType_Nuclear){
				v1 = EigenVecs.col(0);
				v3 = EigenVecs.col(2);
				r2 = EigenVals[1] / EigenVals[0] * r;
				r3 = EigenVals[2] / EigenVals[0] * r;
			}
			else {
				v1 = EigenVecs.col(2);
				v3 = EigenVecs.col(0);
				r2 = EigenVals[1] / EigenVals[2] * r;
				r3 = EigenVals[0] / EigenVals[2] * r;
			}
			ZoneName << setprecision(4) << "EigenValues = (" << EigenVals[0] << "; " << EigenVals[1] << "; " << EigenVals[2] << ")";

			Nodes.resize(NumNodes);
			Nodes.reserve(NumNodes + 3 * NumHyperboloidPts);
			Elems.resize(NumElems, vector<int>(3));

			for (int i = 0; i < NumNodes; ++i) {
				auto & n = Nodes[i];
				n = { meshP[i].x, meshP[i].y, meshP[i].z };
				double rr = r1 * pow(dot(v1, n), 2)
					+ r2 * pow(dot(v2, n), 2)
					+ r3 * pow(dot(v3, n), 2);
				n *= rr;
				n = CPPos + n;
			}
			for (int i = 0; i < NumElems; ++i) {
				Elems[i] = { meshT[i].n1, meshT[i].n2, meshT[i].n3 };
			}

			for (int i = 0; i < NumHyperboloidPts; ++i){
				double alpha = double(i) * AngleStep;
				auto n = Rotate(v1, alpha, v2);
				double rr = r1 * pow(dot(v1, n), 2)
					+ r2 * pow(dot(v2, n), 2)
					+ r3 * pow(dot(v3, n), 2);
				n *= rr;
				n = CPPos + n;
				Nodes.push_back(n);
				Elems.push_back(vector<int>({ NumNodes + i, NumNodes + i, NumNodes + ((i + 1) % NumHyperboloidPts) }));
			}
			NumNodes = Nodes.size();

			for (int i = 0; i < NumHyperboloidPts; ++i) {
				double alpha = double(i) * AngleStep;
				auto n = Rotate(v1, alpha, v3);
				double rr = r1 * pow(dot(v1, n), 2)
					+ r2 * pow(dot(v2, n), 2)
					+ r3 * pow(dot(v3, n), 2);
				n *= rr;
				n = CPPos + n;
				Nodes.push_back(n);
				Elems.push_back(vector<int>({ NumNodes + i, NumNodes + i, NumNodes + ((i + 1) % NumHyperboloidPts) }));
			}
			NumNodes = Nodes.size();

			for (int i = 0; i < NumHyperboloidPts; ++i) {
				double alpha = double(i) * AngleStep;
				auto n = Rotate(v2, alpha, v1);
				double rr = r1 * pow(dot(v1, n), 2)
					+ r2 * pow(dot(v2, n), 2)
					+ r3 * pow(dot(v3, n), 2);
				n *= rr;
				n = CPPos + n;
				Nodes.push_back(n);
				Elems.push_back(vector<int>({ NumNodes + i, NumNodes + i, NumNodes + ((i + 1) % NumHyperboloidPts) }));
			}
			NumNodes = Nodes.size();

		}

		int OldZoneNum = ZoneNumByName(ZoneName.str());
		if (OldZoneNum > 0) {
			DeleteZoneNums += OldZoneNum;
		}

		Surf.MakeFromNodeElemList(Nodes, Elems);
		int ZoneNum = Surf.SaveAsTriFEZone(XYZVarNums, ZoneName.str());
		SetZoneStyle({ ZoneNum }, ZoneStyle_Surface);
		TecUtilZoneSetShade(SV_COLOR, Set(ZoneNum).getRef(), 0.0,
			CPType == CPType_Nuclear ? White_C
			: CPType == CPType_Bond ? Red_C
			: CPType == CPType_Ring ? Green_C
			: Cyan_C);
		TecUtilZoneSetEdgeLayer(SV_COLOR, Set(ZoneNum).getRef(), 0.0,
			CPType == CPType_Nuclear ? (ColorIndex_t)9 // gray
			: CPType == CPType_Bond ? Red_C
			: CPType == CPType_Ring ? Green_C
			: Cyan_C);
		NewZoneNums += ZoneNum;

	}

	TecUtilZoneSetActive(NewZoneNums.getRef(), AssignOp_PlusEquals);
	TecUtilZoneSetContour(SV_SHOW, NewZoneNums.getRef(), 0.0, FALSE);
	TecUtilZoneSetMesh(SV_SHOW, NewZoneNums.getRef(), 0.0, FALSE);
	TecUtilZoneSetShade(SV_SHOW, NewZoneNums.getRef(), 0.0, TRUE);
	TecUtilZoneSetEdgeLayer(SV_SHOW, NewZoneNums.getRef(), 0.0, TRUE);
	TecUtilZoneSetEdgeLayer(SV_LINETHICKNESS, NewZoneNums.getRef(), 0.4, TRUE);
	StyleValue styleValue;
	styleValue.set((Boolean_t)1, NewZoneNums, SV_FIELDMAP, SV_EFFECTS, SV_USETRANSLUCENCY);
	styleValue.set((SmInteger_t)2, NewZoneNums, SV_FIELDMAP, SV_EFFECTS, SV_SURFACETRANSLUCENCY);
	styleValue.set((Boolean_t)1, SV_FIELDLAYERS, SV_SHOWSHADE);
	if (!DeleteZoneNums.isEmpty()){
		TecUtilZoneDelete(DeleteZoneNums.getRef());
	}

	TecUtilDataLoadEnd();


	CSMGUIUnlock();
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

// 	GPParams->GP = GradPath_c(
// 		GPParams->StartPointOrigin + SphericalToCartesian(GPParams->r, gsl_vector_get(ThetaPhi, 0), gsl_vector_get(ThetaPhi, 1)),
// 		GPParams->Direction,
// 		GPParams->NumGPPoints,
// 		GPParams->GPType,
// 		GPParams->GPTerminate,
// 		GPParams->TermPoint,
// 		GPParams->CPs,
// 		GPParams->TermPointRadius,
// 		GPParams->TermValue,
// 		*GPParams->VolInfo,
// 		*GPParams->HessPtrs,
// 		*GPParams->GradPtrs,
// 		*GPParams->RhoPtr);

	GPParams->GP.Clear();
	GPParams->GP.SetStartPoint(GPParams->StartPointOrigin + SphericalToCartesian(GPParams->r, gsl_vector_get(ThetaPhi, 0), gsl_vector_get(ThetaPhi, 1)));

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
	Boolean_t IsPeriodic,
	bool PrecalcVars)
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

	CSMGUILock();

	vec3 StartPoint;

	int const NumSpherePoints = 5000;
	int const NumGPPts = 300;
	double TermRadius = 1e-4;
	double RhoCutoff = DefaultRhoCutoff;

	high_resolution_clock::time_point Time1;
	Time1 = high_resolution_clock::now();

	string StatusStr = "Finding cage-nuclear paths";
	StatusLaunch(StatusStr, AddOnID, TRUE);

	for (int iCP = 0; iCP < SelectedCPNums.size(); ++iCP) {
		if (!StatusUpdate(iCP, SelectedCPNums.size(), StatusStr + ": (" + to_string(iCP + 1) + " of " + to_string(SelectedCPNums.size()) + ")", AddOnID, Time1)) break;
		double StartPointOffset = 0.3 * AllCPs.GetMinCPDist(TypeInd, SelectedCPNums[iCP] - 1, MinDistTypes);

		vector<GradPath_c> GPs;
		vector<vec3> SeedPoints;
		vector<double> SeedTheta, SeedPhi;

		EquidistributedSpherePoints(NumSpherePoints, StartPointOffset, SeedPoints, SeedTheta, SeedPhi);

		for (vec3 & Pt : SeedPoints){
			Pt += AllCPs.GetXYZ(TypeInd, SelectedCPNums[iCP] - 1);
			GPs.push_back(GradPath_c(Pt, GPDir, NumGPPts, GPType_Classic, GPTerminate_AtCP, nullptr, &AllCPs, &TermRadius, &RhoCutoff, VolInfo, HessPtrs, GradPtrs, RhoPtr));
			GPs.back().SetTerminalCPTypeNum(CPTypeNum_Nuclear);
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
			GPParams.TermPoint = nullptr;
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

			GPParams.GP = GPs[0];

			gsl_multimin_fminimizer_type const *T = gsl_multimin_fminimizer_nmsimplex2;
			gsl_multimin_fminimizer *s = nullptr;
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
			AuxDataZoneSetItem(GP.GetZoneNum(), CSMAuxData.GBA.SourceZoneNum, to_string(SelectedCPZoneNum));
		}
	}


	StatusDrop(AddOnID);

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


	CSMGUIUnlock();
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
	
	CSMGUILock();

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
			GPs.push_back(GradPath_c(Pt, GPDir, NumGPPts, GPType_Classic, GPTerminate_AtCP, nullptr, &AllCPs, &TermRadius, &RhoCutoff, VolInfo, HessPtrs, GradPtrs, RhoPtr));
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


	CSMGUIUnlock();
	TecUtilLockFinish(AddOnID);

	return TRUE;
}

/*
*	Shotgun GPs around selected nuclear/cage CPs and save them as ordered zones.
*/
Boolean_t const SimpleSurfacesAroundSaddles(int NumGPs,
	int VolZoneNum,
	vector<int> const & OtherCPZoneNums,
	int SelectedCPZoneNum,
	vector<int> const & SelectedCPNums,
	int CPTypeVarNum,
	vector<int> const & XYZVarNums,
	int RhoVarNum,
	vector<int> & GradVarNums,
	vector<int> const & HessVarNums,
	Boolean_t IsPeriodic,
	double RhoCutoff,
	int NumGPPts,
	bool SavePaths,
	bool SaveSurfaces,
	bool EvenPointSpacing, 
	double SeedOffset,
	bool PrependSaddlePoint,
	bool PrecalcVars)
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
	if ((TypeInd >= 0 && TypeInd != 1 && TypeInd != 2) || (TypeInd < 0 && string(ZoneName) != CSMZoneName.CriticalPoints)) {
		TecUtilDialogErrMsg("Please select bond, ring, or \"Critical Points\" zone");
		return FALSE;
	}
	TecUtilStringDealloc(&ZoneName);

	if (CPTypeVarNum < 0) {
		TecUtilDialogErrMsg("Failed to find CP Type variable");
		return FALSE;
	}

	FieldDataPointer_c RhoPtr;
	vector<FieldDataPointer_c> GradPtrs, HessPtrs;

	TecUtilDataLoadBegin();

	if (!GetReadPtrsForZone(VolZoneNum,
		RhoVarNum, GradVarNums, HessVarNums,
		RhoPtr, GradPtrs, HessPtrs)) {
		TecUtilDialogErrMsg("Failed to get read pointer(s)");
		return FALSE;
	}

	// Grad and Hess needed, so calculate them if they weren't already present
	bool DeleteGradVars = false;
	if (PrecalcVars && GradVarNums.empty()) {
		CalcVarsOptions_s opt;
		opt.AddOnID = AddOnID;
		opt.CalcForAllZones = FALSE;
		opt.CalcZoneNum = VolZoneNum;
		opt.RhoVarNum = RhoVarNum;
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
			for (int i = 1; i < 3; ++i)
				GradVarNums.push_back(GradVarNums[0] + i);

			if (GradVarNums.empty())
				return FALSE;
		}
	}

	VolExtentIndexWeights_s VolInfo;
	VolInfo.AddOnID = AddOnID;

	if (!GetVolInfo(VolZoneNum, XYZVarNums, IsPeriodic, VolInfo)) {
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

	for (int z : OtherCPZoneNums) {
		TecUtilZoneGetIJK(z, &NumCPs, &iJunk[0], &iJunk[1]);
		for (auto const & i : iJunk) if (i > 1) {
			TecUtilDialogErrMsg("CP zone is not i-ordered");
			return FALSE;
		}
	}


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

	for (int z : OtherCPZoneNums) {
		AllCPs += CritPoints_c(z, XYZVarNums, CPTypeVarNum, RhoVarNum, &MR);
	}

	int EndCPNumforName = 0;
	ColorIndex_t PathColor = Red_C;

	CSMGUILock();

	vec3 StartPoint;

	int const NumSpherePoints = NumGPs;
// 	int const NumGPPts = 100;
// 	double RhoCutoff = DefaultRhoCutoff;

	Set NewZones;

	StatusLaunch("Working...", AddOnID, TRUE);

	int numIter = SelectedCPNums.size() * NumGPs;
	int iter = 0;
	auto Time1 = high_resolution_clock::now();

	for (int iCP = 0; iCP < SelectedCPNums.size(); ++iCP) {
		vector<int> CPTypeIndOffset;
		if (TypeInd < 0) CPTypeIndOffset = AllCPs.GetTypeNumOffsetFromTotOffset(SelectedCPNums[iCP] - 1);
		else CPTypeIndOffset = { TypeInd, SelectedCPNums[iCP] - 1 };
		if (CPTypeIndOffset[0] != 1 && CPTypeIndOffset[0] != 2) continue;

		vec3 EigVals = AllCPs.GetEigVals(SelectedCPNums[iCP] - 1);
		mat33 EigVecs = AllCPs.GetEigVecs(SelectedCPNums[iCP] - 1);

		double StartPointOffset = SeedOffset;

		vec3 RotAxis,
			RotVec = normalise(EigVecs.col(1)) * StartPointOffset;
		StreamDir_e GPDir;

		if (AllCPs.GetTypeFromTotOffset(SelectedCPNums[iCP] - 1) == CPType_Bond){
			RotAxis = normalise(EigVecs.col(2));
			GPDir = StreamDir_Reverse;
		}
		else {
			RotAxis = normalise(EigVecs.col(0));
			GPDir = StreamDir_Forward;
		}


		vector<GradPath_c> GPs;
		double RotStep = TWOPI / double(NumGPs);
		double TermRadius = 0.05;

		vector<GradPath_c const *> GPPtrs;
		for (int iGP = 0; iGP < NumGPs; ++iGP){
			vec3 StartPoint = AllCPs.GetXYZ(CPTypeIndOffset[0], CPTypeIndOffset[1]) + Rotate(RotVec, RotStep * double(iGP), RotAxis);
			GPs.push_back(GradPath_c(StartPoint, GPDir, NumGPPts * 5, GPType_Classic, GPTerminate_AtCP, nullptr, &AllCPs, &TermRadius, &RhoCutoff, VolInfo, HessPtrs, GradPtrs, RhoPtr));
			GPs.back().SetStartEndCPNum(AllCPs.GetTotOffsetFromTypeNumOffset(CPTypeIndOffset[0], CPTypeIndOffset[1]), 0);
		}
		int gpIter = 0;
		bool UserQuit = false;

#ifndef _DEBUG
#pragma omp parallel for schedule(dynamic)
#endif
		for (int iGP = 0; iGP < NumGPs; ++iGP) {
			int ThNum = omp_get_thread_num();
			if (ThNum == 0) {
				UserQuit = !StatusUpdate(iter, numIter, "CP " + to_string(iCP + 1) + " of " + to_string(SelectedCPNums.size()) + ": GP " + to_string(gpIter + 1) + " of " + to_string(NumGPs), AddOnID, Time1, false);
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
							StatusLaunch("CP " + to_string(iCP + 1) + " of " + to_string(SelectedCPNums.size()) + ": GP " + to_string(gpIter + 1) + " of " + to_string(NumGPs), AddOnID, TRUE, TRUE, TRUE);
							StatusUpdate(iter, numIter, "CP " + to_string(iCP + 1) + " of " + to_string(SelectedCPNums.size()) + ": GP " + to_string(gpIter + 1) + " of " + to_string(NumGPs), AddOnID, Time1, false);
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
			if (!UserQuit){
#pragma omp atomic
				iter++;
				gpIter++;

				GPs[iGP].Seed(false);
				if (PrependSaddlePoint) {
					GPs[iGP].PointPrepend(AllCPs.GetXYZ(CPTypeIndOffset[0], CPTypeIndOffset[1]), AllCPs.GetRho(CPTypeIndOffset[0], CPTypeIndOffset[1]));
				}
				GPs[iGP].Resample(NumGPPts, EvenPointSpacing ? GPResampleMethod_Linear : GPResampleMethod_Adaptive);
			}
		}
		if (UserQuit){
			TecUtilDataLoadEnd();
			CSMGUIUnlock();
			StatusDrop(AddOnID);
			return FALSE;
		}

		for (int i = 0; i < GPs.size(); ++i) {
			auto & GP = GPs[i];
			if (SavePaths) {
				int zoneNum = GP.SaveAsOrderedZone("Simple Surface " + CPNameList[CPTypeIndOffset[0]] + " " + to_string(CPTypeIndOffset[1]) + " GP " + to_string(i+1));
				NewZones += zoneNum;
				SetZoneStyle({ zoneNum }, ZoneStyle_Path);
			}
			GPPtrs.push_back(&GP);
		}

		if (SaveSurfaces) {
			FESurface_c OutSurf;
			OutSurf.MakeFromGPs(GPPtrs, true, false, false);

			int zoneNum = OutSurf.SaveAsTriFEZone({ 1,2,3 }, "Simple Surface " + CPNameList[CPTypeIndOffset[0]] + " " + to_string(CPTypeIndOffset[1]));
			NewZones += zoneNum;

			SetZoneStyle({ zoneNum }, ZoneStyle_Surface);
		}
	}

	if (!NewZones.isEmpty()) {
		TecUtilZoneSetActive(NewZones.getRef(), AssignOp_PlusEquals);
	}

	StatusDrop(AddOnID);

	TecUtilDataLoadEnd();


	CSMGUIUnlock();
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

void DrawRepresentationQuadricsReturnUserInfo(bool const GuiSuccess,
	vector<GuiField_c> const & Fields,
	vector<GuiField_c> const PassthroughFields) {
	if (!GuiSuccess) return;

	TecUtilLockStart(AddOnID);

	int VolZoneNum, RhoVarNum;
	vector<int> XYZVarNums(3), GradVarNums, HessVarNums;
	Boolean_t IsPeriodic;

	int fNum = 0;

	VolZoneNum = Fields[fNum++].GetReturnInt();
	for (int i = 0; i < 3; ++i) XYZVarNums[i] = i + Fields[fNum].GetReturnInt();
	fNum++;
	RhoVarNum = Fields[fNum++].GetReturnInt();


	if (Fields[fNum++].GetReturnBool()) {
		GradVarNums.resize(3);
		for (int i = 0; i < 3; ++i) GradVarNums[i] = i + Fields[fNum].GetReturnInt();
	}
	fNum++;

	if (Fields[fNum++].GetReturnBool()) {
		HessVarNums.resize(6);
		for (int i = 0; i < 6; ++i) HessVarNums[i] = i + Fields[fNum].GetReturnInt();
	}
	fNum++;

	int CPTypeVarNum = Fields[fNum++].GetReturnInt();
	int SelectCPsZoneNum = Fields[fNum].GetReturnInt();
	vector<int> SelectedCPs = Fields[fNum++].GetReturnIntVec();
	double SizeFactor = Fields[fNum++].GetReturnDouble();

	DrawRepresentationQuadrics(VolZoneNum, SelectCPsZoneNum, SelectedCPs, CPTypeVarNum, XYZVarNums, RhoVarNum, GradVarNums, HessVarNums, SizeFactor);

	TecUtilDialogMessageBox("Finished", MessageBoxType_Information);

	TecUtilLockFinish(AddOnID);
}

void DrawRepresentationQuadricsGetUserInfo() {

	vector<GuiField_c> Fields = {
		GuiField_c(Gui_ZoneSelect, "Volume zone", CSMZoneName.FullVolume.substr(0, 10)),
		GuiField_c(Gui_VarSelect, "X", "X"),
		GuiField_c(Gui_VarSelect, "Electron Density", CSMVarName.Dens)
	};

	int iTmp = Fields.size();
	Fields.push_back(GuiField_c(Gui_ToggleEnable, "Density gradient vector variables present"));

	Fields[iTmp].AppendSearchString(to_string(Fields.size()));
	Fields.push_back(GuiField_c(Gui_VarSelect, CSMVarName.DensGradVec[0], CSMVarName.DensGradVec[0]));

	iTmp = Fields.size();
	Fields.push_back(GuiField_c(Gui_ToggleEnable, "Density Hessian variables present"));

	Fields[iTmp].AppendSearchString(to_string(Fields.size()));
	Fields.push_back(GuiField_c(Gui_VarSelect, CSMVarName.DensHessTensor[0], CSMVarName.DensHessTensor[0]));

	Fields.push_back(GuiField_c(Gui_VarSelect, "CP type variable", CSMVarName.CritPointType));
	string SearchString = CSMZoneName.CriticalPoints;

	Fields.push_back(GuiField_c(Gui_ZonePointSelectMulti, "Source critical point(s)", SearchString));

	Fields.push_back(GuiField_c(Gui_Double, "Size of quadric (as factor of closest CP distance)", "0.3"));

	CSMGui("Draw Representation Quadrics", Fields, DrawRepresentationQuadricsReturnUserInfo, AddOnID);
}

void SimpleSurfacesAroundSaddlesReturnUserInfo(bool const GuiSuccess,
	vector<GuiField_c> const & Fields,
	vector<GuiField_c> const PassthroughFields) {
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

	bool PrecalcVars = Fields[fNum++].GetReturnBool();


	if (Fields[fNum++].GetReturnBool()) {
		GradVarNums.resize(3);
		for (int i = 0; i < 3; ++i) GradVarNums[i] = i + Fields[fNum].GetReturnInt();
	}
	fNum++;
	fNum++;

	if (Fields[fNum++].GetReturnBool()) {
		HessVarNums.resize(6);
		for (int i = 0; i < 6; ++i) HessVarNums[i] = i + Fields[fNum].GetReturnInt();
	}
	fNum += 2;

	int CPTypeVarNum = Fields[fNum++].GetReturnInt();
	vector<int> OtherCPZoneNums = Fields[fNum++].GetReturnIntVec();
	int SelectCPsZoneNum = Fields[fNum].GetReturnInt();
	vector<int> SelectedCPs = Fields[fNum++].GetReturnIntVec();
	int NumGPs = Fields[fNum++].GetReturnInt();
	int PointsPerPath = Fields[fNum++].GetReturnInt();
	double RhoCutoff = Fields[fNum++].GetReturnDouble();

	double SeedOffset = Fields[fNum++].GetReturnDouble();

	bool EvenPointSpacing = Fields[fNum++].GetReturnBool();
	bool PrependSaddlePoint = Fields[fNum++].GetReturnBool();
	bool SavePaths = Fields[fNum++].GetReturnBool();
	bool SaveSurfaces = Fields[fNum++].GetReturnBool();

	SimpleSurfacesAroundSaddles(NumGPs, VolZoneNum, OtherCPZoneNums, SelectCPsZoneNum, SelectedCPs, CPTypeVarNum, XYZVarNums, RhoVarNum, GradVarNums, HessVarNums, IsPeriodic, RhoCutoff, NumGPs, SavePaths, SaveSurfaces, EvenPointSpacing, SeedOffset, PrependSaddlePoint, PrecalcVars);

	TecUtilDialogMessageBox("Finished", MessageBoxType_Information);

	TecUtilLockFinish(AddOnID);
}

void SimpleSurfacesAroundSaddlesGetUserInfo() {
	vector<GuiField_c> Fields = {
		GuiField_c(Gui_ZoneSelect, "Volume zone", CSMZoneName.FullVolume.substr(0, 10)),
		GuiField_c(Gui_Toggle, "Periodic system"),
		GuiField_c(Gui_VertSep)
	};

	Fields.push_back(GuiField_c(Gui_VarSelect, "X", "X"));

	Fields.push_back(GuiField_c(Gui_VertSep));

	Fields.push_back(GuiField_c(Gui_VarSelect, "Electron Density", CSMVarName.Dens));
	Fields.push_back(GuiField_c(Gui_VertSep));

	Fields.push_back(GuiField_c(Gui_Toggle, "Precalculate missing variables", "1"));

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

	Fields.push_back(GuiField_c(Gui_Int, "Number of gradient paths", "360"));

	Fields.push_back(GuiField_c(Gui_Int, "Number of points per path", "500"));

	Fields.push_back(GuiField_c(Gui_Double, "Path rho cutoff", "0.001"));

	Fields.push_back(GuiField_c(Gui_Double, "Seed point offset from saddle point", "0.01"));

	Fields.push_back(GuiField_c(Gui_Toggle, "Even point spacing along gradient paths", "0"));

	Fields.push_back(GuiField_c(Gui_Toggle, "Prepend saddle CP as first point on paths", "1"));

	Fields.push_back(GuiField_c(Gui_Toggle, "Save gradient paths as zones", "0"));

	Fields.push_back(GuiField_c(Gui_Toggle, "Save gradient paths as surface", "1"));

	CSMGui("Simple surfaces (and GPs) around bond/ring CPs", Fields, SimpleSurfacesAroundSaddlesReturnUserInfo, AddOnID);
}

void ExportPathDataReturnUserInfo(bool const GuiSuccess,
	vector<GuiField_c> const & Fields,
	vector<GuiField_c> const PassthroughFields) {
	if (!GuiSuccess) return;

	TecUtilLockStart(AddOnID);

	vector<int> ExportZoneNums;
	vector<int> XYZVarNums(3);
	int RhoVarNum;
	int fNum = 0;

	ExportZoneNums = Fields[fNum++].GetReturnIntVec();
	for (int i = 0; i < 3; ++i) XYZVarNums[i] = i + Fields[fNum].GetReturnInt();
	RhoVarNum = Fields[++fNum].GetReturnInt();
	XYZVarNums.push_back(RhoVarNum);

	TecUtilLockStart(AddOnID);
	TecUtilPleaseWait("Exporting... Please wait.", TRUE);
	CSMGUILock();

	char* FileNameCStr;

	vector<string> ColumnHeadings = {
		"Zone number",
		"Zone name",
		"Number of points",
		"Length",
		"Total curvature (sum of angle changes along path) [radians]",
		"Average curvature k (total curvature over length) [radians / length]",
		"Net curvature (angle difference between endpoints; first and last 1% of path) [radians]",
		"Average net curvature (net curvature over length) [radians / length]",
		"Total torsion (sum of angle changes of binormal along path) [radians]",
		"Average torsion (total torsion over length) [radians / length]",
		"Total curvature-scaled torsion (sum of angle changes of binormal along path) [radians]",
		"Average curvature-scaled torsion (total torsion over length) [radians / length]"
	};

	if (TecUtilDialogGetFileName(SelectFileOption_WriteFile, &FileNameCStr, "Comma separated values", "Bondalyzer exported path data.csv", "*.csv")) {
		string FileName = FileNameCStr;
		TecUtilStringDealloc(&FileNameCStr);

		ofstream OutFile(FileName.c_str(), std::ios::trunc);
		if (OutFile.is_open()){
			OutFile << StringJoin(ColumnHeadings, ",");
			OutFile << endl;
			for (auto zi : ExportZoneNums){
				int ijk[3];
				TecUtilZoneGetIJK(zi, &ijk[0], &ijk[1], &ijk[2]);
				if (TecUtilZoneIsOrdered(zi) && ijk[0] > 2 && ijk[1] == 1 && ijk[2] == 1) {
					GradPath_c GP(zi, XYZVarNums, AddOnID);
					double l = GP.GetLength(),
						k = GP.ComputeTotalCurvature(),
						t = GP.ComputeTotalTorsion(),
						t1 = GP.ComputeTotalTorsion(true);

					int startInd = MAX(1, GP.GetIndAtLength(0.01 * l)),
						endInd = MIN(GP.GetCount() - 2, GP.GetIndAtLength(0.99 * l));
					double k1 = VectorAngleMagnitude(GP[startInd] - GP[0], GP[-1] - GP[-endInd - 1]);

					char* ZoneNameCStr;
					TecUtilZoneGetName(zi, &ZoneNameCStr);
					OutFile << zi
						<< "," << ZoneNameCStr
						<< "," << GP.GetCount()
						<< "," << l
						<< "," << k
						<< "," << k / l
						<< "," << k1
						<< "," << k1 / l
						<< "," << t
						<< "," << t / l
						<< "," << t1
						<< "," << t1 / l
						<< endl;
					TecUtilStringDealloc(&ZoneNameCStr);
				}
			}

			OutFile.close();

			TecUtilDialogMessageBox("Finished", MessageBox_Information);
		}
		else{
			TecUtilDialogErrMsg("Failed to open file for writing! (perhaps it's open in another program?)");
		}

		
	}

	CSMGUIUnlock();
	TecUtilPleaseWait("Exporting... Please wait.", FALSE);

	TecUtilLockFinish(AddOnID);
}

void ExportPathDataGetUserInfo() {
	vector<GuiField_c> Fields = {
		GuiField_c(Gui_ZoneSelectMulti, "Zones to export", ""),
		GuiField_c(Gui_VarSelect, "X", "X"),
		GuiField_c(Gui_VarSelect, "Electron Charge Density", "Electron Density")
	};


	CSMGui("Export gradient path data", Fields, ExportPathDataReturnUserInfo, AddOnID);
}

double MinFunc_GPLength_InPlane(double alpha, void * params) {
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

double MinFunc_GPLength_InPlane2(double alpha, void * params) {
	MinFuncParams_GPLengthInPlane2 * GPParams = reinterpret_cast<MinFuncParams_GPLengthInPlane2*>(params);

	vec3 StartPoint = *GPParams->Origin + Rotate(*GPParams->RotVec, alpha, *GPParams->RotAxis);
	GPParams->GP->Clear();
	GPParams->GP->SetStartPoint(StartPoint);

	GPParams->GP->Seed(false);

	if (!GPParams->GP->IsMade())
		TecUtilDialogErrMsg("Gradient path in MinFunc_GPLength_InPlane failed");
	
	GPParams->GP->PointPrepend(*GPParams->Origin, 0);

	return GPParams->GP->GetLength();
}

double MinFunc_GPLength_BetweenPoints(double Path1Weight, void * params){
	MinFuncParams_GPLengthBetweenPoints * GPParams = reinterpret_cast<MinFuncParams_GPLengthBetweenPoints*>(params);

	GPParams->GP->Clear();

	vec3 StartPoint = *GPParams->StartPoints->first * Path1Weight + *GPParams->StartPoints->second * (1.0 - Path1Weight);
	GPParams->GP->SetStartPoint(StartPoint);

	if (GPParams->GP->GetSurfPtr() != nullptr) {
		GradPath_c BackGP = *GPParams->GP;
		BackGP.SetStartPoint(StartPoint);
		BackGP.SetDir(BackGP.GetDir() == StreamDir_Forward ? StreamDir_Reverse : StreamDir_Forward);
		BackGP.SetGPTermType(GPTerminate_AtPoint);
		BackGP.SetTermPointRadius(1e-4);

		BackGP.Seed(false);

// 		GPParams->GP->SetStartPoint(BackGP[0]);
		GPParams->GP->SetSurfPtr(nullptr);
		GPParams->GP->Seed(false);
		GPParams->GP->SetSurfPtr(BackGP.GetSurfPtr());

		*GPParams->GP += BackGP;
	}
	else
		GPParams->GP->Seed(false);


	if (!GPParams->GP->IsMade())
		TecUtilDialogErrMsg("Gradient path in MinFunc_GPLength_InPlane failed");

	double Length = GPParams->GP->GetLength();

	return Length;
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

bool MinLengthGPInPlane(GradPath_c & GP, 
	vec3 const & Origin, 
	vec3 const & RotVec, 
	vec3 const & RotAxis, 
	double const & aLow, 
	double const & aHigh,
	double & aMin,
	double & lenMin)
{
	MinFuncParams_GPLengthInPlane2 GPParams;
	GPParams.GP = &GP;
	GPParams.Origin = &Origin;
	GPParams.RotVec = &RotVec;
	GPParams.RotAxis = &RotAxis;

	gsl_function F;
	F.function = &MinFunc_GPLength_InPlane2;
	F.params = reinterpret_cast<void*>(&GPParams);

	gsl_min_fminimizer_type const *T;
	T = gsl_min_fminimizer_brent;

	gsl_min_fminimizer *s;
	s = gsl_min_fminimizer_alloc(T);

	gsl_min_fminimizer_set(s, &F, (aLow + aHigh) * 0.5, aLow, aHigh);

	int Iter = 0;
	int Status;

	double OldLen = DBL_MAX, NewLen;

	do
	{
		Status = gsl_min_fminimizer_iterate(s);
		if (Status == GSL_EBADFUNC || Status == GSL_FAILURE) {
			TecUtilDialogErrMsg("GSL one-dimensional minimization failed");
			break;
		}
		NewLen = GPParams.GP->GetLength();

		Status = (abs(OldLen - NewLen) < Tolerance_GPLengthInPlane) ? GSL_SUCCESS : GSL_CONTINUE;

		OldLen = NewLen;
	} while (Status == GSL_CONTINUE && ++Iter < MaxIter_GPLengthInPlane);

	if (Status == GSL_SUCCESS) {
		aMin = s->x_minimum;
		lenMin = s->f_minimum;
		vec3 StartPoint = Origin + Rotate(RotVec, aMin, RotAxis);
		GP.Clear();
		GP.SetStartPoint(StartPoint);
		GP.Seed(false);
		GP.PointPrepend(Origin, 0.0);
	}

	gsl_min_fminimizer_free(s);

	return (Status == GSL_SUCCESS);
}

bool MinLengthGPBetweenPointsCheckStartPoints(GradPath_c & GP, std::pair<vec3 const *, vec3 const *> & StartPoints, double & bMin)
{
	MinFuncParams_GPLengthBetweenPoints GPParams;
	GPParams.GP = &GP;
	GPParams.StartPoints = &StartPoints;

	double lowLen = MinFunc_GPLength_BetweenPoints(0.0, reinterpret_cast<void*>(&GPParams)),
		highLen = MinFunc_GPLength_BetweenPoints(1.0, reinterpret_cast<void*>(&GPParams)),
		minLen = MinFunc_GPLength_BetweenPoints(bMin, reinterpret_cast<void*>(&GPParams));

	return (minLen < lowLen && minLen < highLen);
}

bool MinLengthGPBetweenPoints(GradPath_c & GP, std::pair<vec3 const *,vec3 const *> & StartPoints, double & bMin)
{
	MinFuncParams_GPLengthBetweenPoints GPParams;
	GPParams.GP = &GP;
	GPParams.StartPoints = &StartPoints;

// 	double lowLen = MinFunc_GPLength_BetweenPoints(0.0, reinterpret_cast<void*>(&GPParams)),
// 		highLen = MinFunc_GPLength_BetweenPoints(1.0, reinterpret_cast<void*>(&GPParams)),
// 		minLen = MinFunc_GPLength_BetweenPoints(bMin, reinterpret_cast<void*>(&GPParams));
// 
// 	while(minLen > lowLen || minLen > highLen){
// 		if (minLen > lowLen){
// 			bMin = (bMin + 1.0) * 0.5;
// 		}
// 		else{
// 			bMin = bMin * 0.5;
// 		}
// 		minLen = MinFunc_GPLength_BetweenPoints(bMin, reinterpret_cast<void*>(&GPParams));
// 	}

	gsl_function F;
	F.function = &MinFunc_GPLength_BetweenPoints;
	F.params = reinterpret_cast<void*>(&GPParams);

	gsl_min_fminimizer_type const *T;
	T = gsl_min_fminimizer_brent;

	gsl_min_fminimizer *s;
	s = gsl_min_fminimizer_alloc(T);

	gsl_min_fminimizer_set(s, &F, bMin, 0.0, 1.0);

	int Iter = 0;
	int Status;

	double OldLen = DBL_MAX, NewLen;
	double OldWeight = 0.5, NewWeight;

	do 
	{
		Status = gsl_min_fminimizer_iterate(s);
		if (Status == GSL_EBADFUNC || Status == GSL_FAILURE) {
			TecUtilDialogErrMsg("GSL one-dimensional minimization failed");
			break;
		}
		NewLen = GPParams.GP->GetLength();
		NewWeight = s->x_minimum;

// 		Status = gsl_min_test_interval(OldLen, NewLen, Tolerance_GPLengthInPlane, Tolerance_GPLengthInPlane);
// 		if (Status != GSL_SUCCESS)
// 			gsl_min_test_interval(OldWeight, NewWeight, Tolerance_GPLengthInPlane, Tolerance_GPLengthInPlane);

		Status = (abs(OldLen - NewLen) < Tolerance_GPLengthInPlane) ? GSL_SUCCESS : GSL_CONTINUE;

		OldLen = NewLen;
		OldWeight = NewWeight;
	} while (Status == GSL_CONTINUE && ++Iter < MaxIter_GPLengthInPlane);

	if (Status == GSL_SUCCESS) {
		bMin = s->x_minimum;
		vec3 minStartPoint = *StartPoints.first * bMin + *StartPoints.second * (1.0 - bMin);
		GP.Clear();
		GP.SetStartPoint(minStartPoint);
		GP.Seed(false);
	}

	gsl_min_fminimizer_free(s);

	return (Status == GSL_SUCCESS);
}

double MidPtStdDevScalar(std::deque<vec3> const & MidPtDeque) {
	vec3 mean = zeros<vec>(3);
	double stddev = 0.0;
	for (auto const & i : MidPtDeque)
		mean += i;
	mean /= (double)MidPtDeque.size();

	for (auto const & i : MidPtDeque)
		stddev += Distance(mean, i);
	stddev /= (double)MidPtDeque.size();

	return stddev;
};


GradPath_c MinLengthGPFromGPTripletWithCommonTerminus(vector<GradPath_c> LeftMidRightGPs, int MaxNumIter, double MinGPLengthTolerence){
	bool IsOk = LeftMidRightGPs.size() == 3 && (MaxNumIter > 0 || MinGPLengthTolerence > 0);
	bool MinFound = false;
	if (IsOk) {
		for (auto gp : LeftMidRightGPs){
			IsOk = IsOk && gp.IsMade() && gp.GetStartEndCPNum(1) >= 0;
		}
		for (int i = 0; i < 2 && IsOk; ++i) {
			IsOk = (LeftMidRightGPs[i].GetStartEndCPNum(1) == LeftMidRightGPs[i + 1].GetStartEndCPNum(1));
		}
		MinFound = IsOk && LeftMidRightGPs[0].GetLength() > LeftMidRightGPs[1].GetLength() && LeftMidRightGPs[1].GetLength() < LeftMidRightGPs[2].GetLength();
	}

	if (!MinFound ){
		return GradPath_c();
	}

	FESurface_c Surf;
	vector<GradPath_c const *> GPPtrs;
	for (auto const & gp : LeftMidRightGPs) {
		GPPtrs.push_back(&gp);
	}
	Surf.MakeFromGPs(GPPtrs);
	Surf.GeneratePointElementDistanceCheckData();

	GradPath_c *OutGP;
	GradPath_c BaseGP = LeftMidRightGPs[1];
	BaseGP.Clear();
	BaseGP.SetSurfPtr(&Surf);
	BaseGP.SetDir(StreamDir_Both);
	BaseGP.SetTerminalCPTypeNums({ CPTypeNum_Nuclear, CPTypeNum_Bond, CPTypeNum_Ring, CPTypeNum_Cage });

	std::list<GradPath_c> GPList(LeftMidRightGPs.begin(), LeftMidRightGPs.end());

	auto Left = GPList.begin();
	auto Center = Left;
	Center++;
	auto Right = Center;
	Right++;

	int Iter = 0;
	double OldLength = DBL_MAX, NewLength = 0, LengthAbsDiff = abs(OldLength - NewLength);
	bool IsConverged = false;

	while (IsOk && MinFound && (MaxNumIter > 0 && Iter < MaxNumIter) && !IsConverged){
		Iter++;
		MinFound = false;
		while (!MinFound && Right != GPList.end()){
			MinFound = Left->GetLength() > Center->GetLength() && Center->GetLength() < Right->GetLength();
			if (!MinFound){
				Left++; 
				Center++; 
				Right++;
			}
		}
		

		if (MinFound){
			// Seed two paths between Left and Center and between Center and Right, then position the iterators so that the paths checked on the next loop are Left--Left-Center--Right.
			OutGP = &*Center;
			NewLength = Center->GetLength();
			LengthAbsDiff = abs(OldLength - NewLength);
			IsConverged = (MinGPLengthTolerence > 0 && LengthAbsDiff < MinGPLengthTolerence);

			if (!IsConverged) {
				vector<std::pair<GradPath_c*, GradPath_c*> > GPPairs;
				GPPairs.emplace_back(&*Center, &*Right);
				GPPairs.emplace_back(&*Left, &*Center);
				bool FirstPass = true;
				for (auto GPPair : GPPairs) {
					GradPath_c TmpGP = BaseGP;
					vec3 SeedPt, Pt2;
					int Ind1, Ind2;
					double SepDist;
					IsOk = GPPair.first->GetMaxSeparationMidpointFromOtherGP(*GPPair.second, SeedPt, Ind1, Ind2, Pt2, SepDist);
					if (IsOk) {
						TmpGP.SetStartPoint(SeedPt);
						TmpGP.Seed();
						IsOk = TmpGP.IsMade() && TmpGP.GetStartEndCPNum() == GPPair.first->GetStartEndCPNum();
					}
					if (IsOk) {
						if (FirstPass) {
							GPList.insert(Right, TmpGP);
							FirstPass = false;
						}
						else {
							GPList.insert(Center, TmpGP);
						}
					}

					if (!IsOk) {
						break;
					}
				}
				if (IsOk) {
					Center = Left;
					Center++;
					Right = Center;
					Right++;

					// Center now points to the new path, and Right to the old Center
					OldLength = NewLength;
				}
			}
		}
	}
	/*
	if (IsOk && LeftMidRightGPs.front()->GetMaxSeparationMidpointFromOtherGP(*LeftMidRightGPs.back(), GPMidPt, Path1Ind, Path2Ind, Path2Pt, MaxDist)
			&& ((!Path1->IsSGP && !Path2->IsSGP && MaxDist > CheckDist)
				|| MaxDist > CheckDist * 0.8))
	{
		// 				GradPath_c TmpGP({ &Path1->GP,&Path2->GP }, GPMidPt, (CPType == CPType_Ring ? StreamDir_Reverse : StreamDir_Forward), { std::make_pair(0, MIN(Path1Ind + NeighborGPSurfaceIndOffset, Path1->GP.GetCount() - 1)), std::make_pair(0, MIN(Path2Ind + NeighborGPSurfaceIndOffset, Path2->GP.GetCount() - 1)) }, &CPPosition);
		GradPath_c TmpGP({ &Path1->GP,&Path2->GP }, GPMidPt, (CPType == CPType_Ring ? StreamDir_Reverse : StreamDir_Forward), { std::make_pair(0, MAX(Path1->GP.GetIndAtLength(Path1->GP.GetLength(Path1Ind) + Path1->GP.GetLength() * NeighborGPSurfaceIndLengthOffsetFactor), Path1Ind + NeighborGPSurfaceIndOffset)), std::make_pair(0, MAX(Path2->GP.GetIndAtLength(Path2->GP.GetLength(Path2Ind) + Path2->GP.GetLength() * NeighborGPSurfaceIndLengthOffsetFactor), Path2Ind + NeighborGPSurfaceIndOffset)) }, &CPPosition);

		GPWithInfo GP;
		GP.CPType = Path1->CPType;
		GP.CPNum = Path1->CPNum;
		GP.Identifier = (Path1->Identifier + Path2->Identifier) * 0.5;
		// 				GP.GP = GradPath_c({ &Path1->GP, &Path2->GP }, 0.5, Path1->GP.RhoAt(0), 1);
		GP.GP = Path1->GP;
		GP.GP.Clear();
		GP.GP.SetStartPoint(TmpGP[0]);
		GP.GP.SetDir(GPDir);
		GP.GP.SetSurfPtr(nullptr);
		GP.GP.SetTerminalCPTypeNum(GPTerminalCPTypeNum);
		GP.GP.SetGPTermType(GPTerminate_AtRhoValue);
		GP.GP.Seed(false);

		int tmpInd = 0;
		while ((CPType == CPType_Ring && GP.GP.RhoAt(tmpInd) < TmpGP.RhoAt(0)) || (CPType == CPType_Bond && GP.GP.RhoAt(tmpInd) > TmpGP.RhoAt(0)))
			tmpInd++;
		if (tmpInd > 0)
			GP.GP = GP.GP.SubGP(tmpInd, -1);

		GP.GP += TmpGP;
		GP.GP.Resample(NumGPPts + 5);
		GP.GP.RemoveKinks();
		GP.GP.Resample(NumGPPts);
		GP.GP.ReinterpolateRhoValuesFromVolume(&VolInfo, &RhoPtr);

		if (DistSqr(GP.GP[0], CPPosition) > DistSqr(GP.GP[-1], CPPosition))
			// Need to reverse GP
			GP.GP.Reverse();
		// 				if (CPType == CPType_Ring && GP.GP.RhoAt(0) > GP.GP.RhoAt(-1)) {
		// 					// Need to reverse GP
		// 					GP.GP.Reverse();
		// 				}
		// 				
#ifdef DEBUG_PRINTPATHS
		SubPassNum++;
		DebugPathStringBase = "Path1 id:" + to_string(Path1->Identifier) + ", CP:" + to_string(Path1->CPNum) + "\nPath2 id:" + to_string(Path2->Identifier) + ", CP:" + to_string(Path2->CPNum) + "\nFiller path";
		DoPrintPaths = DoPrintPaths && DebugSaveGPWaitForUser(GP.GP, XYZVarNums, RhoVarNum, DebugPathStringBase, false, PassNum, SubPassNum, 1);
#endif

		GPList.insert(Path2, GP);
		Path2--;
	}
	}*/

	return *OutGP;
}


/*
*	Create interatomic (bond) and ring surfaces from a source CP zone(s)
*	Also creates bond-ring, and bond-cage (ring-nuclear) paths.
*/
Boolean_t FindBondRingSurfaces2(int VolZoneNum,
	vector<int> const & OtherCPZoneNums,
	int SelectedCPZoneNum,
	vector<int> SelectedCPNums,
	int CPTypeVarNum,
	CPType_e const & CPType,
	vector<int> const & XYZVarNums,
	int RhoVarNum,
	vector<int> GradVarNums,
	vector<int> HessVarNums,
	int RCSFuncVarNum,
	Boolean_t IsPeriodic,
	bool DebugMode,
	double RhoCutoff,
	int NumGPPts,
	bool PrecalcVars)
{
	TecUtilLockStart(AddOnID);

	int NumZones = TecUtilDataSetGetNumZones();
	int NumVars = TecUtilDataSetGetNumVars();

	REQUIRE(XYZVarNums.size() == 3);
	for (auto const & i : XYZVarNums) REQUIRE(i > 0 && i <= NumVars);
	for (auto const & i : OtherCPZoneNums) REQUIRE(i > 0 && i <= NumZones);

	int TypeInd = VectorGetElementNum(CPTypeList, CPType);
	if (TypeInd >= 6) {
		TecUtilDialogErrMsg("Invalid CP type specified");
		return FALSE;
	}

	if (CPTypeVarNum < 0) {
		TecUtilDialogErrMsg("Failed to find CP Type variable");
		return FALSE;
	}

	// Grad and Hess needed, so calculate them if they weren't already present
	bool DeleteGradVars = false;
	if (PrecalcVars && GradVarNums.empty()) {
		CalcVarsOptions_s opt;
		opt.AddOnID = AddOnID;
		opt.CalcForAllZones = FALSE;
		opt.CalcZoneNum = VolZoneNum;
		opt.RhoVarNum = RhoVarNum;
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
			for (int i = 1; i < 3; ++i)
				GradVarNums.push_back(GradVarNums[0] + i);

			if (GradVarNums.empty())
				return FALSE;
		}
	}

	FieldDataPointer_c RhoPtr, RCSFuncPtr;
	vector<FieldDataPointer_c> GradPtrs, HessPtrs;

	vector<int> XYZRhoVarNums = XYZVarNums;
	XYZRhoVarNums.push_back(RhoVarNum);

	TecUtilDataLoadBegin();

	if (RCSFuncVarNum > 0 && !RCSFuncPtr.InitializeReadPtr(VolZoneNum, RCSFuncVarNum)) {
		TecUtilDialogErrMsg("Failed to get read pointers");
		return FALSE;
	}

	if (!GetReadPtrsForZone(VolZoneNum,
		RhoVarNum, GradVarNums, HessVarNums,
		RhoPtr, GradPtrs, HessPtrs)) {
		TecUtilDialogErrMsg("Failed to get read pointer(s)");
		return FALSE;
	}

	VolExtentIndexWeights_s VolInfo;
	VolInfo.AddOnID = AddOnID;

	if (!GetVolInfo(VolZoneNum, XYZVarNums, IsPeriodic, VolInfo)) {
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

	for (int z : OtherCPZoneNums) {
		TecUtilZoneGetIJK(z, &NumCPs, &iJunk[0], &iJunk[1]);
		for (auto const & i : iJunk) if (i > 1) {
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

	bool UseTotalCPOffset = AllCPs.NumAtoms() > 0;

	for (int z : OtherCPZoneNums) {
		AllCPs += CritPoints_c(z, XYZVarNums, CPTypeVarNum, RhoVarNum, &MR);
	}

	AllCPs.RemoveSpuriousCPs();

	if (SelectedCPNums.size() == 0) {
		for (int i = 0; i < AllCPs.NumCPs(TypeInd); ++i) SelectedCPNums.push_back(i + 1);
	}

	CSMGUILock();

	vector<CPType_e> CPTypesForMinDistCheck;
	CPTypeNum_e GPTerminalCPSaddleTypeNum = CPTypeNum_Ring,
		GPTerminalCPTypeNum = CPTypeNum_Cage;
	vector<CPTypeNum_e> GPTerminalCPTypeNums;
	int EndCPPosition = 1;
	StreamDir_e GPDir = StreamDir_Reverse;
	ColorIndex_t PathColor = Black_C;
	CPType_e PrimaryCPTermType = CPType_Cage,
		SaddleCPTermType = CPType_Ring;
	if (CPType == CPType_Ring) {
		GPDir = StreamDir_Forward;
		PathColor = Green_C;
		// 		EndCPPosition++;
		PrimaryCPTermType = CPType_Nuclear;
		SaddleCPTermType = CPType_Bond;
		CPTypesForMinDistCheck = { CPType_Bond, CPType_Ring };
		GPTerminalCPTypeNum = CPTypeNum_Nuclear;
		GPTerminalCPSaddleTypeNum = CPTypeNum_Bond;
	}
	else {
		if (AllCPs.NumBonds() <= 1) CPTypesForMinDistCheck.push_back(CPType_Nuclear);

		CPTypesForMinDistCheck.push_back(CPType_Bond);

		if (AllCPs.NumRings() > 0) CPTypesForMinDistCheck.push_back(CPType_Ring);
	}
 	GPTerminalCPTypeNums = { GPTerminalCPTypeNum, GPTerminalCPSaddleTypeNum };
//	GPTerminalCPTypeNums = { GPTerminalCPSaddleTypeNum };

// 	int const NumGPPts = 800;
	// Increase number of GP points when working, then resample it to specified number
	double NumGPPtsWorkingFactor = 5.0;
	NumGPPts *= NumGPPtsWorkingFactor;
	double TermRadius = 0.02;
	double SaddleTermRadius = TermRadius * 5.0;
	double CoincidentPointCheckDistSqr = 1e-12;
	double CoincidentPointCheckAngle = 1e-40;
	double GPMidPointConvergenceTolerance = 1e-6;
	double GPMidPointUseCurvatureMidpointTolerance = 0.05;
	int GPMidPointConvergenceNumCheckPts = 10;
// 	double RhoCutoff = DefaultRhoCutoff;
	int NumCircleCheckGPs = MinNumCircleCheckGPs;
	int NumCircleGPs = MinNumCircleGPs;
	double NeighborGPSurfaceIndLengthOffsetFactor = 0.001;
	int NeighborGPSurfaceIndOffset = 1;



	double StartPointOffset;
	double AngleStep;
	AngleStep = PI2 / static_cast<double>(NumCircleCheckGPs);
	double GPDeviationCheckAngleFactor = DefaultGPDeviationCheckAngleFactor;
	double GPDeviationCheckDistFactor = DefaultGPDeviationCheckDistFactor;
	double GPDeviationStartPtLengthFactor = DefaultGPDeviationStartPtLengthFactor;


	double DelXYZNorm = norm(VolInfo.PointSpacingV123);

	NumCPs = SelectedCPNums.size();

	struct GPWithInfo{
		GradPath_c GP = GradPath_c();
		double Identifier = -1.0;
		int CPNum = -1;
		int SaddleCPNum = -1;
		CPType_e CPType = CPType_Invalid;
		bool IsSGP = false;
		vector<GradPath_c> SaddleEndGPs;
	};

	auto DebugSaveGPList = [](std::list<GPWithInfo> & GPList, int PassNum, vector<int> XYZVarNums, int RhoVarNum, CritPoints_c const & AllCPs, bool SaveAllGPs = true, bool ActivateZones = false)
	{
		string PassStr = "Pass " + to_string(PassNum) + ": ";
		Set ZoneSet;
		for (auto & p : GPList) {
			GradPath_c * tmpGP = &p.GP;
			if (SaveAllGPs || p.GP.GetZoneNum() <= 0) {
				vector<int> StartEndCPNums = tmpGP->GetStartEndCPNum();
				vector<vector<int> > StartEndCPTypeAndOffset(2);
				string Name = PassStr + (p.IsSGP ? "SGP: " : "GP: ");
				for (int i = 0; i < 2; ++i) {
					if (StartEndCPNums[i] >= 0) {
						StartEndCPTypeAndOffset[i] = AllCPs.GetTypeNumOffsetFromTotOffset(StartEndCPNums[i]);
						Name += CPNameList[StartEndCPTypeAndOffset[i][0]] + " " + to_string(StartEndCPTypeAndOffset[i][1] + 1);
					}
					else Name += "FF";

					if (i == 0) Name += " to ";
				}

				tmpGP->SaveAsOrderedZone(Name, vector<FieldDataType_e>(), XYZVarNums, RhoVarNum, FALSE, ColorIndex_t(1 + (PassNum % 5)));
				for (int i = 0; i < 2; ++i) {
					AuxDataZoneSetItem(tmpGP->GetZoneNum(), CSMAuxData.CC.GPEndNumStrs[i], to_string(StartEndCPNums[i] + 1));
					if (StartEndCPNums[i] >= 0) AuxDataZoneSetItem(tmpGP->GetZoneNum(), CSMAuxData.CC.GPEndTypes[i], CPNameList[StartEndCPTypeAndOffset[i][0]]);
					else AuxDataZoneSetItem(tmpGP->GetZoneNum(), CSMAuxData.CC.GPEndTypes[i], "FF");
				}
				AuxDataZoneSetItem(tmpGP->GetZoneNum(), CSMAuxData.CC.ZoneType, CSMAuxData.CC.ZoneTypeGP);
				AuxDataZoneSetItem(tmpGP->GetZoneNum(), "Seed angle", to_string(p.Identifier));
				AuxDataZoneSetItem(tmpGP->GetZoneNum(), "Length", to_string(p.GP.GetLength()));
				ZoneSet += tmpGP->GetZoneNum();
			}
		}
		if (ActivateZones)
			TecUtilZoneSetActive(ZoneSet.getRef(), AssignOp_PlusEquals);

#ifdef DEBUG_SAVESCREENSHOTS
		TecUtilRedrawAll(TRUE);
		TecUtilMacroExecuteCommand("$!EXPORTSETUP EXPORTFORMAT = JPEG");
		TecUtilMacroExecuteCommand("$!EXPORTSETUP IMAGEWIDTH = 1000");
		TecUtilMacroExecuteCommand("$!EXPORTSETUP USESUPERSAMPLEANTIALIASING = NO");

		TecUtilMacroExecuteCommand(string("$!EXPORTSETUP EXPORTFNAME = \'C:\\Users\\Haiiro\\Desktop\\pass " + to_string(PassNum) + " full.jpg\'").c_str());
		TecUtilMacroExecuteCommand("$!EXPORT EXPORTREGION = CURRENTFRAME");

// 		sprintf(buf, "$!EXPORTSETUP EXPORTFORMAT = PNG\n$!EXPORTSETUP IMAGEWIDTH = 2000\n$!EXPORTSETUP USESUPERSAMPLEANTIALIASING = YES\n$!EXPORTSETUP EXPORTFNAME = \'W:\\Tecplot\\StorageWorkspace\\2019_GBAAlgotithm\\RSScreens\\pass %d full.png\'\n$!EXPORT \n  EXPORTREGION = CURRENTFRAME", PassNum);
// 		TecUtilMacroExecuteCommand(buf);
// 		delete buf;
#endif
	};

	auto DebugSaveGPWaitForUser = [](GradPath_c GP, vector<int> XYZVarNums, int RhoVarNum, string statusStr = "", bool DeleteGP = false, int PassNum = -1, int SubPassNum = -1, int SubSubPassNum = -1)
	{
#ifdef DEBUG_PRINTPATHS
		CSMGUIUnlock();
		int ZoneNum = GP.SaveAsOrderedZone("Temp GP Zone");
		Boolean_t DoContinue;
		if (ZoneNum > 0) {
			Set tmpSet(ZoneNum);
			// 		tmpSet += ZoneNum;
			SetZoneStyle({ ZoneNum }, ZoneStyle_Path, Red_C, 0.5);
			TecUtilZoneSetActive(tmpSet.getRef(), AssignOp_PlusEquals);
			TecUtilRedrawAll(TRUE);
			TecUtilRedrawAll(TRUE);
			
#ifndef DEBUG_SAVESCREENSHOTS
			DoContinue = TecUtilDialogMessageBox(string(statusStr + "\nWaiting to continue or quit...").c_str(), MessageBox_YesNo);
#else
			DoContinue = TRUE;
			if (PassNum >= 0 && SubPassNum >= 0) {
				TecUtilMacroExecuteCommand("$!EXPORTSETUP EXPORTFORMAT = JPEG");
				TecUtilMacroExecuteCommand("$!EXPORTSETUP IMAGEWIDTH = 1000");
				TecUtilMacroExecuteCommand("$!EXPORTSETUP USESUPERSAMPLEANTIALIASING = NO");

				TecUtilMacroExecuteCommand(string("$!EXPORTSETUP EXPORTFNAME = \'C:\\Users\\Haiiro\\Desktop\\pass " + to_string(PassNum) + "-" + to_string(SubPassNum) + "-" + to_string(SubSubPassNum) + ".jpg\'").c_str());
				TecUtilMacroExecuteCommand("$!EXPORT EXPORTREGION = CURRENTFRAME");

// 				sprintf(buf, "$!EXPORTSETUP EXPORTFORMAT = PNG\n$!EXPORTSETUP IMAGEWIDTH = 2000\n$!EXPORTSETUP USESUPERSAMPLEANTIALIASING = YES\n$!EXPORTSETUP EXPORTFNAME = \'W:\\Tecplot\\StorageWorkspace\\2019_GBAAlgotithm\\RSScreens\\pass %d-%04d-%d.png\'\n$!EXPORT \n  EXPORTREGION = CURRENTFRAME", PassNum, SubPassNum, SubSubPassNum);
// 				TecUtilMacroExecuteCommand(buf);
// 				delete buf;
			}
#endif // !DEBUG_SAVESCREENSHOTS

			if (DeleteGP && !tmpSet.isEmpty())
				TecUtilZoneDelete(tmpSet.getRef());

			SetZoneStyle({ ZoneNum }, ZoneStyle_Path);
		}
		else{
			DoContinue = TecUtilDialogMessageBox(string(statusStr + "\nFAILED TO PRINT GP\nWaiting to continue or quit...").c_str(), MessageBox_YesNo);
		}
		CSMGUILock();

		return DoContinue;
#else
		return true;
#endif
	};

	auto GPIncrement = [](std::list<GPWithInfo>::iterator & GP, std::list<GPWithInfo> & GPList)
	{
		if (GP == GPList.end()){
			GP = GPList.begin();
// 			GP++;
		}
		else {
			GP++;
			if (GP == GPList.end()) {
				GP = GPList.begin();
				// 			GP++;
			}
		}
	};
	auto GPDecrement = [](std::list<GPWithInfo>::iterator & GP, std::list<GPWithInfo> & GPList)
	{
		if (GP == GPList.begin()) {
			GP = GPList.end();
			GP--;// --;
		}
		else {
			GP--;
		}
	};


	high_resolution_clock::time_point Time1;
	Time1 = high_resolution_clock::now();

	string StatusStr = "Finding " + string(CPType == CPType_Bond ? "interatomic" : "ring") + " surfaces";
	StatusLaunch(StatusStr, AddOnID, TRUE);

	// Get set of CPs connected by GPs (First neighbor only)
	AllCPs.GenerateCPGraph();

	int StepNum = 0;
	vector<string> StepStrings = {
		": Seeding initial GPs",
		": Finding saddle-saddle paths",
		": Filling sparse regions",
		": Preparing saddle-terminal CP paths",
		": Finding saddle-terminal CP paths"
	};
	int NumSteps = StepStrings.size();
	for (int i = 0; i < NumSteps; ++i)
		StepStrings[i] += " (Step " + to_string(i + 1) + " of " + to_string(NumSteps) + ")";
	int StatusTotalNum = NumCPs * NumSteps;

	for (int iCP = 0; iCP < NumCPs; ++iCP) {
		int DebugNumZonesBefore = TecUtilDataSetGetNumZones();

		string CPStatusStr = StatusStr + ": (" + to_string(iCP + 1) + " of " + to_string(NumCPs) + ")";
		StepNum = 0;
		if (!StatusUpdate(iCP * NumSteps + StepNum, StatusTotalNum, CPStatusStr + StepStrings[StepNum++], AddOnID, Time1)) break;
		int CPNum = SelectedCPNums[iCP] - 1;

		if (UseTotalCPOffset) {
			vector<int> TypeOffset = AllCPs.GetTypeNumOffsetFromTotOffset(CPNum);
			CPNum = TypeOffset[1];
		}

		double MinCPDist = AllCPs.GetMinCPDist(TypeInd, CPNum, CPTypesForMinDistCheck);

		TermRadius = 0.01 * MinCPDist;
		SaddleTermRadius = TermRadius * 50.0;

		StartPointOffset = MIN(0.05 * MinCPDist, DelXYZNorm);
		double SourceCPTermRadius = 0.05 * StartPointOffset;
		vec3 StartVec = normalise(AllCPs.GetEigVecs(TypeInd, CPNum).col(1)) * StartPointOffset;
		vec3 CPPosition = AllCPs.GetXYZ(TypeInd, CPNum);
		double CPRho = AllCPs.GetRho(TypeInd, CPNum);
		vec3 RotAxis = AllCPs.GetPrincDir(TypeInd, CPNum);
		int CurrentCPNum = AllCPs.GetTotOffsetFromTypeNumOffset(TypeInd, CPNum);

		std::list<GPWithInfo> GPList;
		vector<GradPath_c> GPs, BackGPs;
		GPs.reserve(NumCircleCheckGPs);
		BackGPs.reserve(NumCircleCheckGPs);
		vector<double> GPLengths(NumCircleCheckGPs);
		vector<CPType_e> GPEndCPTypes(NumCircleCheckGPs);
		vector<int> GPEndCPNums(NumCircleCheckGPs);

		vector<vec3> SeedDiscPts;
		SeedDiscPts.reserve(NumCircleCheckGPs + 1);
		SeedDiscPts.push_back(CPPosition);

		/*
		 *	Seed GPs around the CP.
		 *	These GPs will only terminate at nuclear/cage CPs.
		 *	Bond/ring CPs and missing nuclear/cage CPs will be found using these results.
		 */
		TermRadius *= 15.0;
		for (int gpNum = 0; gpNum < NumCircleCheckGPs; ++gpNum) {
			double RotAngle = static_cast<double>(gpNum) * AngleStep;

			vec3 SeedPt = Rotate(StartVec, RotAngle, RotAxis);
			SeedDiscPts.push_back(CPPosition + SeedPt * 1.5);
			SeedPt += CPPosition;

			GPs.push_back(GradPath_c(SeedPt, GPDir, NumGPPts, GPType_Classic, GPTerminate_AtCP, nullptr, &AllCPs, &StartPointOffset, &RhoCutoff, VolInfo, HessPtrs, GradPtrs, RhoPtr));
			GPs.back().SetTerminalCPTypeNums(GPTerminalCPTypeNums);
			GPs.back().SetStartEndCPNum(CurrentCPNum, 0);

			BackGPs.push_back(GradPath_c(SeedPt, (GPDir == StreamDir_Forward ? StreamDir_Reverse : StreamDir_Forward), NumGPPts, GPType_Classic, GPTerminate_AtCP, nullptr, &AllCPs, &SourceCPTermRadius, &RhoCutoff, VolInfo, HessPtrs, GradPtrs, RhoPtr));
			BackGPs.back().SetGPTermType(GPTerminate_AtPoint);
			BackGPs.back().SetTermPoint(CPPosition);
			BackGPs.back().SetTermPointRadius(1e-4);
		}
		TermRadius /= 15.0;

		// Prepare FESurface for the disc used to seed GPs. 
		// This surface will be used to keep the BackGPs constrained to the 
		// planer approximation being used local to the origin saddle CP.
		vector<vector<int> > SeedDiscTris;
		SeedDiscTris.reserve(NumCircleCheckGPs);
		for (int i = 1; i < NumCircleCheckGPs; ++i)
			SeedDiscTris.push_back({ 0, i, i + 1 });
		SeedDiscTris.push_back({ 0,NumCircleCheckGPs,1 });
		FESurface_c SeedDisc(SeedDiscPts, SeedDiscTris);
		SeedDisc.GeneratePointElementDistanceCheckData();
		for (auto & i : BackGPs)
			i.SetSurfPtr(&SeedDisc);

		if (DebugMode){
			TecUtilDialogMessageBox("Seeding initial GPs to define surface", MessageBoxType_Information);
		}

// 		SeedDisc.SaveAsTriFEZone(XYZVarNums, "seed disc");
		// 		StatusDrop(AddOnID);
		// 		CSMGuiUnlock();
		// 		return TRUE;

// #ifndef _DEBUG
#pragma omp parallel for schedule(dynamic)
// #endif
		for (int gpNum = 0; gpNum < NumCircleCheckGPs; ++gpNum) {
			GPs[gpNum].Seed(false);
			// 			GPs[gpNum].PointPrepend(CPPosition, AllCPs.GetRho(TypeInd, cpNum));
			BackGPs[gpNum].Seed(false);
			int BackGPEndCPNum = BackGPs[gpNum].GetStartEndCPNum(EndCPPosition);
			CPType_e BackGPEndType = AllCPs.GetTypeFromTotOffset(BackGPEndCPNum);
			if (BackGPEndType != CPType) {
				// GP back to CP at center of surface didn't terminate at the CP, so
				// find the closest point on the GP to the CP and terminate it there.
				int ClosestPtInd;
				vec3 ClosestPt = BackGPs[gpNum].ClosestPoint(CPPosition, ClosestPtInd);
				BackGPs[gpNum] = BackGPs[gpNum].SubGP(0, ClosestPtInd - 1);
				BackGPs[gpNum].PointAppend(CPPosition, CPRho);
				BackGPs[gpNum].SetStartEndCPNum(CurrentCPNum, 1);
			}

			GPs[gpNum] += BackGPs[gpNum];

			GPs[gpNum].Resample(NumGPPts);

			GPs[gpNum].ReinterpolateRhoValuesFromVolume();

			GPs[gpNum].MakeRhoValuesMonotonic();

// 			if ((CPType == CPType_Ring && GPs[gpNum].RhoAt(0) > GPs[gpNum].RhoAt(-1))
// 				|| (CPType == CPType_Bond && GPs[gpNum].RhoAt(0) < GPs[gpNum].RhoAt(-1))) {
// 				GPs[gpNum].Reverse();
// 			}
			if (DistSqr(GPs[gpNum][0], CPPosition) > DistSqr(GPs[gpNum][-1], CPPosition))
				// Need to reverse GP
				GPs[gpNum].Reverse();

			GPLengths[gpNum] = GPs[gpNum].GetLength();
			GPEndCPNums[gpNum] = GPs[gpNum].GetStartEndCPNum(EndCPPosition);
			GPEndCPTypes[gpNum] = AllCPs.GetTypeFromTotOffset(GPEndCPNums[gpNum]);
		}

		// Reseed GPs with smaller temrinal radius

// #ifndef _DEBUG
#pragma omp parallel for schedule(dynamic)
// #endif
		for (int gpNum = 0; gpNum < NumCircleCheckGPs; ++gpNum) {
			if (GPEndCPTypes[gpNum] != SaddleCPTermType) {
				GPs[gpNum].Clear();
				GPs[gpNum].SetTerminalCPTypeNum(GPTerminalCPTypeNum);
				GPs[gpNum].SetTermPointRadius(TermRadius);
				GPs[gpNum].Seed(false);

				GPs[gpNum] += BackGPs[gpNum];

				GPs[gpNum].Resample(NumGPPts);

				GPs[gpNum].ReinterpolateRhoValuesFromVolume();

				GPs[gpNum].MakeRhoValuesMonotonic();
			}
		}


		// Add GPs to map, but only if not saddle-terminating
		for (int gpNum = 0; gpNum < NumCircleCheckGPs; ++gpNum) {
			if (GPEndCPTypes[gpNum] != SaddleCPTermType) {
				double RotAngle = static_cast<double>(gpNum) * AngleStep;
				GPList.push_back(GPWithInfo());
// 				GPList.back().Identifier = RotAngle;
				GPList.back().Identifier = double(gpNum * 10000);
				GPList.back().GP = GPs[gpNum];
				GPList.back().CPNum = GPEndCPNums[gpNum];
				GPList.back().CPType = GPEndCPTypes[gpNum];
			}
		}
		GPs.clear();
		BackGPs.clear();
		GPLengths.clear();
		GPEndCPNums.clear();
		GPEndCPTypes.clear();

		// Copy first GP to end of list and set its angle to 2Pi.
		// a circle of GPs.
		GPList.insert(GPList.end(), GPList.front());
		GPList.back().Identifier = PI2 + GPList.front().Identifier;

		// DEBUG
		int PassNum = -1;
		bool DoContinue = true;
		bool DoPrintPaths = DebugMode;

#ifdef DEBUG_PRINTPATHS
		int SubPassNum = 0;
		if (DebugMode) {
			DebugSaveGPList(GPList, ++PassNum, XYZVarNums, RhoVarNum, AllCPs, true, true);
		}
#endif
//       
//  		DoContinue = false;
//  		if (!DoContinue) break;
//  		


		// Doing a removal to paths too close to their neighbors before the saddle-saddle step


		int Path1Ind, Path2Ind;
		vec3 Path2Pt, GPMidPt, GPMidPointOld;
		double Weight;
		double CheckDist = 0.2 * MinCPDist;;
		double MaxDist;

		auto Path1 = GPList.begin();
		auto Path2 = Path1;
		auto Path3 = Path2;


		string DebugPathStringBase;


		// Pass to remove GPs that are consistently too close to their neighbors.
		// If the paths to the left and right of a test path are too close, then remove
		// the test path.

// 		Path1 = GPList.begin();
// 		Path2 = Path1;
// 		Path2++;
// 		Path3 = Path2;
// 		Path3++;
// 		CheckDist = 0.01 * MinCPDist;
// 		while (Path2 != GPList.end()) {
// 			if (!Path1->IsSGP && !Path2->IsSGP
// 				&& (Path1->CPNum == Path2->CPNum
// 					|| (Path1->CPType == CPType_Invalid
// 						&& Path2->CPType == CPType_Invalid))
// 				&& Path1->GP.GetMaxSeparationMidpointFromOtherGP(Path2->GP, GPMidPt, Path1Ind, Path2Ind, Path2Pt, MaxDist)
// 				&& MaxDist < CheckDist)
// 			{
// 				// Paths are too close, replace them with a path at their midpoint
// 
// 				GPWithInfo GP;
// 				GP.CPType = Path1->CPType;
// 				GP.CPNum = Path1->CPNum;
// 				GP.Angle = (Path1->Angle + Path2->Angle) * 0.5;
// 				// 				GP.GP = GradPath_c({ &Path1->GP, &Path2->GP }, 0.5, Path1->GP.RhoAt(0), 1);
// 				GP.GP = Path1->GP;
// 				GP.GP.Clear();
// 				GP.GP.SetStartPoint(GPMidPt);
// 				GP.GP.SetDir(GPDir);
// 				GP.GP.SetSurfPtr(nullptr);
// 				GP.GP.SetGPTermType(GPTerminate_AtPoint);
// 				GP.GP.SetTermPoint(AllCPs.GetXYZ(GP.CPNum));
// 				GP.GP.SetTermPointRadius(1e-4);
// 				GP.GP.Seed(false);
// 
// 				GradPath_c TmpGP({ &Path1->GP,&Path2->GP }, GPMidPt, (CPType == CPType_Ring ? StreamDir_Reverse : StreamDir_Forward), { std::make_pair(0, MAX(Path1->GP.GetIndAtLength(Path1->GP.GetLength(Path1Ind) + Path1->GP.GetLength() * NeighborGPSurfaceIndLengthOffsetFactor), Path1Ind + NeighborGPSurfaceIndOffset)), std::make_pair(0, MAX(Path2->GP.GetIndAtLength(Path2->GP.GetLength(Path2Ind) + Path2->GP.GetLength() * NeighborGPSurfaceIndLengthOffsetFactor), Path2Ind + NeighborGPSurfaceIndOffset)) }, &CPPosition);
// 				GP.GP += TmpGP;
// 				GP.GP.RemoveKinks();
// 				GP.GP.Resample(NumGPPts);
// 				GP.GP.ReinterpolateRhoValuesFromVolume(&VolInfo, &RhoPtr);
// 
// 				if (CPType == CPType_Ring && GP.GP.RhoAt(0) > GP.GP.RhoAt(-1)) {
// 					// Need to reverse GP
// 					GP.GP.Reverse();
// 				}
// 
// 				GPList.insert(Path2, GP);
// 
// 				Path3 = Path2;
// 				Path3--;
// 				GPList.erase(Path2);
// 				GPList.erase(Path1);
// 				Path1 = Path3;
// 				Path2 = Path3;
// 				Path2++;
// 			}
// 			else {
// 				Path1++;
// 				Path2++;
// 			}
// 		}
// 		
		if (DebugMode) {
			TecUtilRedrawAll(TRUE);
			TecUtilDialogMessageBox("Removing GPs too close together (GP triplet check)", MessageBoxType_Information);
		}

		Path1 = GPList.begin();
		Path2 = Path1;
		Path2++;
		Path3 = Path2;
		Path3++;
		CheckDist = 0.04 * MinCPDist;
		vector<int> IndList;
		while (DoContinue &&  Path1 != GPList.end()) {
			DoContinue = StatusUpdate(iCP * NumSteps + StepNum, StatusTotalNum, CPStatusStr + StepStrings[StepNum], AddOnID, Time1);
			if ((!Path2->IsSGP
				&& (Path1->CPNum == Path2->CPNum
					&& Path1->CPNum == Path3->CPNum)
				|| (Path1->CPType == CPType_Invalid
					&& Path2->CPType == CPType_Invalid
					&& Path3->CPType == CPType_Invalid))
				//  				&& Path1->SaddleCPNum < 0 && Path3->SaddleCPNum < 0
				&& Path2->GP.GetMaxSeparationFromNeighboringGPs(vector<GradPathBase_c const *>({ &Path1->GP, &Path3->GP }), IndList, Path2Ind) < CheckDist)
				//  				&& Path1->GP.GetMaxSeparationMidpointFromOtherGP(Path3->GP, GPMidPt, Path1Ind, Path2Ind, Path2Pt, MaxDist)
				//  				&& MaxDist < CheckDist)
			{
				// Paths are too close, so remove the path between them.

#ifdef DEBUG_PRINTPATHS
				SubPassNum++;
				DebugPathStringBase = "Path1 id:" + to_string(Path1->Identifier) + ", CP:" + to_string(Path1->CPNum) + "\nPath2 id:" + to_string(Path2->Identifier) + ", CP:" + to_string(Path2->CPNum) + "\nPath3 id:" + to_string(Path3->Identifier) + ", CP:" + to_string(Path3->CPNum) + "\nRemoving path 2";
// 				DoPrintPaths = DoPrintPaths && DebugSaveGPWaitForUser(Path2->GP, XYZVarNums, RhoVarNum, DebugPathStringBase, true, PassNum, SubPassNum, 1);
#endif
				GPList.erase(Path2);
				Path2 = Path3;
				GPIncrement(Path3, GPList);
			}
			else {
				Path1++;
				GPIncrement(Path2, GPList);
				GPIncrement(Path3, GPList);
			}
		}


		// DEBUG
#ifdef DEBUG_PRINTPATHS
		SubPassNum = 0;
		DoPrintPaths = DebugMode;
		if (DoPrintPaths) {
			{
// 				int DebugNumZonesAfter = TecUtilDataSetGetNumZones();
// 				Set DebugDeleteZoneSet;
// 				for (int i = DebugNumZonesBefore + 1; i <= DebugNumZonesAfter; ++i)
// 					DebugDeleteZoneSet += i;
// 
// 				if (!DebugDeleteZoneSet.isEmpty())
// 					TecUtilZoneDelete(DebugDeleteZoneSet.getRef());
			}
			DebugSaveGPList(GPList, ++PassNum, XYZVarNums, RhoVarNum, AllCPs, true, true);
		}
#endif
		// 
		//  		DoQuit = true;
		if (!DoContinue) break;

		if (DebugMode) {
			TecUtilRedrawAll(TRUE);
			TecUtilDialogMessageBox("Removing GPs too close together (GP neighbor check)", MessageBoxType_Information);
		}

		Path1 = GPList.begin();
		Path2 = Path1;
		Path2++;
		CheckDist = 0.015 * MinCPDist;
		while (DoContinue && Path1 != GPList.end()) {
			DoContinue = StatusUpdate(iCP * NumSteps + StepNum, StatusTotalNum, CPStatusStr + StepStrings[StepNum], AddOnID, Time1);
			if ((!Path2->IsSGP
				&& Path1->CPNum == Path2->CPNum
				|| (Path1->CPType == CPType_Invalid
					&& Path2->CPType == CPType_Invalid))
				&& Path1->SaddleCPNum < 0 && Path3->SaddleCPNum < 0
				&& Path1->GP.GetMaxSeparationMidpointFromOtherGP(Path2->GP, GPMidPt, Path1Ind, Path2Ind, Path2Pt, MaxDist)
				&& MaxDist < CheckDist)
			{
#ifdef DEBUG_PRINTPATHS
				SubPassNum++;
				DebugPathStringBase = "Path1 id:" + to_string(Path1->Identifier) + ", CP:" + to_string(Path1->CPNum) + "\nPath2 id:" + to_string(Path2->Identifier) + ", CP:" + to_string(Path2->CPNum) + "\nRemoving path 2";
// 				DoPrintPaths = DoPrintPaths && DebugSaveGPWaitForUser(Path2->GP, XYZVarNums, RhoVarNum, DebugPathStringBase, true, PassNum, SubPassNum, 1);
#endif
				// Paths are too close, so remove the path between them.
				GPList.erase(Path2);
				Path2 = Path1;
				GPIncrement(Path2, GPList);
			}
			else {
				Path1++;
				GPIncrement(Path2, GPList);
			}
		}


		// DEBUG
#ifdef DEBUG_PRINTPATHS
		SubPassNum = 0;
		DoPrintPaths = DebugMode;
		if (DoPrintPaths) {
			{
				int DebugNumZonesAfter = TecUtilDataSetGetNumZones();
				Set DebugDeleteZoneSet;
				for (int i = DebugNumZonesBefore + 1; i <= DebugNumZonesAfter; ++i)
					DebugDeleteZoneSet += i;

				if (!DebugDeleteZoneSet.isEmpty())
					TecUtilZoneDelete(DebugDeleteZoneSet.getRef());
			}
			DebugSaveGPList(GPList, ++PassNum, XYZVarNums, RhoVarNum, AllCPs, true, true);
		}
#endif


// 		int Path1Ind, Path2Ind;
// 		vec3 Path2Pt, GPMidPt, GPMidPointOld;
// 		double Weight;
// 		double CheckDist = 0.2 * MinCPDist;;
// 		double MaxDist;
// 		auto Path1 = GPList.begin();
// 		auto Path2 = Path1;
// 		auto Path3 = Path2;


		CheckDist = 0.2 * MinCPDist;

		Path1 = GPList.begin();
		Path2 = Path1;
		Path3 = Path2;


		// Pass to get bond-ring paths and any missing bond-cage or ring-nuclear paths.
		// Step through the GPMap (around the bond/ring point) to check for GPs that 
		// terminate at different CPs. Only nuclear/cage-terminal GPs were made, so
		// any deviation should be between CPs that don't already share a bond/ring path.
		// When a deviation is found, start bisecting the area between the GPs once they
		// achieve a cutoff deviation angle.

		Path1 = GPList.begin();
		Path2 = Path1;
		Path2++;
		std::deque<vec3> GPMidPtDeque;
		std::set<int> FoundSaddleCPNums;


		DoContinue = StatusUpdate(iCP * NumSteps + StepNum, StatusTotalNum, CPStatusStr + StepStrings[StepNum], AddOnID, Time1);

		if (DebugMode) {
			TecUtilRedrawAll(TRUE);
			TecUtilDialogMessageBox("Finding Saddle-saddle paths", MessageBoxType_Information);
		}


		while (DoContinue && Path2 != GPList.end()) {
			TermRadius = 0.01 * MinCPDist;
			DoContinue = StatusUpdate(iCP * NumSteps + StepNum, StatusTotalNum, CPStatusStr + StepStrings[StepNum] + " found " + to_string(FoundSaddleCPNums.size()) + " saddle(s)", AddOnID, Time1);
			if (Path2->IsSGP);
			{
				// It's possible that a bond CP will, against the odds, be found when none of its
				// neighboring nuclear CPs have yet been discovered. This makes it harder to ensure
				// that the corresponding bond-nuclear paths are ordered correctly.
				// So if Path2 is a SGP, might have to switch the order of Path2's SaddleEndSGPs based
				// on the SaddleEndGPs end-point distance to the terminal CPs of Path1 and Path3
				if (Path2->SaddleEndGPs.size() == 2 && Path2->SaddleEndGPs[0].IsMade() && Path2->SaddleEndGPs[1].IsMade()) {
					Path3 = Path2;
					Path3++;
					if (DistSqr(Path2->SaddleEndGPs[0][-1], Path3->GP[-1]) < DistSqr(Path2->SaddleEndGPs[1][-1], Path1->GP[-1])){
						Path2->SaddleEndGPs = { Path2->SaddleEndGPs[1], Path2->SaddleEndGPs[0] };
					}
				}

			}
			if (((!Path1->IsSGP && !Path2->IsSGP)
				|| (Path1->IsSGP && Path1->SaddleEndGPs[1].GetStartEndCPNum(1) != Path2->CPNum)
				|| (Path2->IsSGP && Path2->SaddleEndGPs[0].GetStartEndCPNum(1) != Path1->CPNum))
				&& ((Path1->CPNum >= 0 && Path2->CPNum >= 0 && Path1->CPNum != Path2->CPNum)
					|| (Path1->GP.GetMaxSeparationMidpointFromOtherGP(Path2->GP, GPMidPt, Path1Ind, Path2Ind, Path2Pt, MaxDist) && MaxDist > CheckDist))
				&& Path1->GP.GetDeviationMidpointAsTerminalAngleAndDistanceFactorFromOtherGP(Path2->GP, GPDeviationCheckAngleFactor, GPDeviationCheckDistFactor, GPDeviationStartPtLengthFactor, GPMidPt, Path1Ind, Path2Ind, Path2Pt, MinCPDist * 0.05))
			{
				// Paths deviate. Need to seed a new GP between them once their
				// deviation angle is path the check angle.
				GPWithInfo GP;

				// DEBUG
				GP.Identifier = (Path1->Identifier + Path2->Identifier) * 0.5;

#ifdef DEBUG_PRINTPATHS
				SubPassNum++;
				if (DebugMode) {
					DebugPathStringBase = "Path1 id:" + to_string(Path1->Identifier) + ", CP:" + to_string(Path1->CPNum) + "\nPath2 id:" + to_string(Path2->Identifier) + ", CP:" + to_string(Path2->CPNum);
				}
#endif

				double MidPtRhoVal = ValAtPointByPtr(GPMidPt, VolInfo, RhoPtr);
				double AngleDiff = abs(Path1->Identifier - Path2->Identifier);

				// 				if (AngleDiff < CoincidentPointCheckAngle){

				// 				double TmpDistSqr = DistSqr(Path1->GP.XYZAt(Path1Ind), Path2Pt);
				// 				double TmpDistSqr = DistSqr(Path1->GP.ClosestPoint(GPMidPt), Path2->GP.ClosestPoint(GPMidPt));
				// 				double TmpCheckDistSqr = CoincidentPointCheckDistSqr;
				// 				if (Path1Ind < Path1->GP.GetCount() - 2)
				// 					TmpCheckDistSqr = MAX(TmpCheckDistSqr, DistSqr(Path1->GP.XYZAt(Path1Ind), Path1->GP.XYZAt(Path1Ind + 1)));
				// 				else if (Path1Ind > 0)
				// 					TmpCheckDistSqr = MAX(TmpCheckDistSqr, DistSqr(Path1->GP.XYZAt(Path1Ind), Path1->GP.XYZAt(Path1Ind - 1)));
				// 
				// 				if (Path2Ind < Path2->GP.GetCount() - 2)
				// 					TmpCheckDistSqr = MAX(TmpCheckDistSqr, DistSqr(Path2->GP.XYZAt(Path2Ind), Path2->GP.XYZAt(Path2Ind + 1)));
				// 				else if (Path2Ind > 0)
				// 					TmpCheckDistSqr = MAX(TmpCheckDistSqr, DistSqr(Path2->GP.XYZAt(Path2Ind), Path2->GP.XYZAt(Path2Ind - 1)));

				GPMidPtDeque.push_back(GPMidPt);
				while (GPMidPtDeque.size() > GPMidPointConvergenceNumCheckPts)
					GPMidPtDeque.pop_front();

				bool IsConverged = (GPMidPtDeque.size() >= 3 && MidPtStdDevScalar(GPMidPtDeque) < GPMidPointConvergenceTolerance || AngleDiff < CoincidentPointCheckAngle);
				// 
				// 				if (IsConverged)
				// 					break;

				if (IsConverged) {
					GPMidPt = (Path1->GP.ClosestMaxCurvaturePoint(GPMidPt, Path1Ind) + Path2->GP.ClosestMaxCurvaturePoint(GPMidPt, Path2Ind)) * 0.5;
					GPMidPtDeque.push_back(GPMidPt);
					while (GPMidPtDeque.size() > GPMidPointConvergenceNumCheckPts)
						GPMidPtDeque.pop_front();
				}

				// 				double MidPtStdDev = MidPtStdDevScalar(GPMidPtDeque);

				// 				if (TmpDistSqr < TmpCheckDistSqr && AngleDiff < CoincidentPointCheckAngle) {
				if (IsConverged && !AllCPs.HasEdge(Path1->CPNum, Path2->CPNum, 2)) {
// 					if (DebugMode) {
// 						TecUtilDialogMessageBox("Finding Saddle-saddle paths: deviating: converged", MessageBoxType_Information);
// 					}
					// Path1 and Path2 are essentially the same, deviating at a saddle
					// that wasn't found or isn't near-field, at GPMidPt.
					// Use the Path1-Path2 midpoint GP as the saddle GP, and use the SubGPs
					// Path1 and Path2 as the end GPs.
					// 
// 					GP.GP = GradPath_c({ &Path1->GP,&Path2->GP }, 0.5, MidPtRhoVal, -1);
//  					GP.GP = GradPath_c({ &Path1->GP,&Path2->GP }, GPMidPt, (CPType == CPType_Ring ? StreamDir_Reverse : StreamDir_Forward), { std::make_pair(0, MIN(Path1Ind + NeighborGPSurfaceIndOffset, Path1->GP.GetCount() - 1)), std::make_pair(0, MIN(Path2Ind + NeighborGPSurfaceIndOffset, Path2->GP.GetCount() - 1)) }, &CPPosition);
					GP.GP = GradPath_c({ &Path1->GP,&Path2->GP }, 
						GPMidPt, 
						(CPType == CPType_Ring ? StreamDir_Reverse : StreamDir_Forward), 
						{ std::make_pair(0, MAX(Path1->GP.GetIndAtLength(Path1->GP.GetLength(Path1Ind) + Path1->GP.GetLength() * NeighborGPSurfaceIndLengthOffsetFactor), Path1Ind + NeighborGPSurfaceIndOffset)), std::make_pair(0, MAX(Path2->GP.GetIndAtLength(Path2->GP.GetLength(Path2Ind) + Path2->GP.GetLength() * NeighborGPSurfaceIndLengthOffsetFactor), Path2Ind + NeighborGPSurfaceIndOffset)) }, 
						&CPPosition);
					GP.GP.ReinterpolateRhoValuesFromVolume(&VolInfo, &RhoPtr);
					GP.Identifier = (Path1->Identifier + Path2->Identifier) * 0.5;
					GP.IsSGP = true;
					GP.CPType = SaddleCPTermType;
					GP.SaddleEndGPs.push_back(Path1->GP.SubGP(MIN(Path1Ind + 2, Path1->GP.GetCount() - 2), -1));
					GP.SaddleEndGPs.push_back(Path2->GP.SubGP(MIN(Path2Ind + 2, Path2->GP.GetCount() - 2), -1));

					DoPrintPaths = DoPrintPaths && DebugSaveGPWaitForUser(GP.GP, XYZVarNums, RhoVarNum, DebugPathStringBase + "\nTerminal CP: " + to_string(GP.CPNum) + ", type: " + to_string(GP.CPType) + "\nConverged non-saddle path", PassNum, SubPassNum, 1);
					DoPrintPaths = DoPrintPaths && DebugSaveGPWaitForUser(GP.SaddleEndGPs[0], XYZVarNums, RhoVarNum, DebugPathStringBase + "\nTerminal CP: " + to_string(GP.CPNum) + ", type: " + to_string(GP.CPType) + "\nConverged non-saddle path: saddle GP 1", PassNum, SubPassNum, 2);
					DoPrintPaths = DoPrintPaths && DebugSaveGPWaitForUser(GP.SaddleEndGPs[1], XYZVarNums, RhoVarNum, DebugPathStringBase + "\nTerminal CP: " + to_string(GP.CPNum) + ", type: " + to_string(GP.CPType) + "\nConverged non-saddle path: saddle GP 2", PassNum, SubPassNum, 3);

					if (DistSqr(GP.GP[0], CPPosition) > DistSqr(GP.GP[-1], CPPosition))
						// Need to reverse GP
						GP.GP.Reverse();
// 					if (CPType == CPType_Ring && GP.GP.RhoAt(0) > GP.GP.RhoAt(-1)) {
// 						// Need to reverse GP
// 						GP.GP.Reverse();
// 					}


					GPList.insert(Path2, GP);

					// Remove the paths to the left and right of the new path.
					Path2--; // Points to new GP in GPList.
					Path3 = Path2;
					Path3--;
					GPList.erase(Path3);
					Path3 = Path2;
					Path3++;
					GPList.erase(Path3);

					// Resume the search to the right of Path2
					Path1 = Path2;
					Path2++;
				}
				else {
					// Path1 and Path2 deviate, so seed a GP at the deviation midpoint and see where it goes.

					// Get closest bond point to deviation point, then use the closest CP distance to that point to adjust the termination radius
					int CloseCPInd;
					double CloseCPDist;
					vec ClosePt = AllCPs.ClosestPoint(GPMidPt, CloseCPInd, CloseCPDist);
					double TmpMinDist = AllCPs.GetMinCPDist(CloseCPInd);
					double TmpSaddleTermRadius = TmpMinDist * 0.1;
					TermRadius = TmpSaddleTermRadius;

					GP.GP.SetupGradPath(GPMidPt, GPDir, NumGPPts, GPType_Classic, GPTerminate_AtCP, nullptr, &AllCPs, &TermRadius, &RhoCutoff, VolInfo, HessPtrs, GradPtrs, RhoPtr);
					GP.GP.SetTerminalCPTypeNums(GPTerminalCPTypeNums);
					GP.GP.SetStartEndCPNum(CurrentCPNum, 0);
					GP.GP.Seed(false);

					//DEBUG
// 					GP.GP.SaveAsOrderedZone();

					GP.CPNum = GP.GP.GetStartEndCPNum(EndCPPosition);
					GP.CPType = AllCPs.GetTypeFromTotOffset(GP.CPNum);

 					// Do distance check to see if path skimmed saddle point
// 					int ClosePtNum = GP.GP.GetCount() - 1;
// 					CloseCPInd = -1;
// 					double ClosePtDist;
//   					for (int jCP = 0; jCP < AllCPs.NumCPs(GPTerminalCPSaddleTypeNum); ++jCP){
// 						int TmpClosePtNum;
// 						auto TmpDist = Distance(GP.GP.ClosestPoint(AllCPs.GetXYZ(GPTerminalCPSaddleTypeNum, jCP), TmpClosePtNum), AllCPs.GetXYZ(GPTerminalCPSaddleTypeNum, jCP));
// 						if (TmpDist <= TmpSaddleTermRadius) {
// 							int TmpCPNum = AllCPs.GetTotOffsetFromTypeNumOffset(GPTerminalCPSaddleTypeNum, jCP);
// 							if (!FoundSaddleCPNums.count(TmpCPNum) && TmpClosePtNum < ClosePtNum) {
// 								CloseCPInd = TmpCPNum;
// 								ClosePtNum = TmpClosePtNum;
// 								ClosePtDist = TmpDist;
// 							}
//   						}
//   					}
// 					if (CloseCPInd >= 0){
// 						GP.CPNum = CloseCPInd;
// 						GP.CPType = SaddleCPTermType;
// 						GP.GP = GP.GP.SubGP(0, ClosePtNum);
// 						GP.GP.PointAppend(AllCPs.GetXYZ(CloseCPInd), AllCPs.GetRho(CloseCPInd));
// 					}
					

					bool Reiterate = false;

					if (IsConverged && GP.CPType != SaddleCPTermType) {
// 						if (DebugMode) {
// 							TecUtilDialogMessageBox("Finding Saddle-saddle paths: deviating: converged at non-saddle", MessageBoxType_Information);
// 						}
						// The deviation check has stalled, but we know there exists a bond/ring path
						// connection between the Path1 and Path2 terminal CPs. Stitch the existing
						// bond/ring path together with Path1 and Path2 to form a surface that
						// can be used to seed a GP from the saddle CP back into the central CP.
						// So, gather the two bond/ring path segments, then make the surface, seed
						// back to the central CP, save to GPList, and move on in the search.



						// First get the bond/ring paths
						vector<int> CPNums;
						AllCPs.HasEdge(Path1->CPNum, Path2->CPNum, 2, &CPNums);
						GP.SaddleCPNum = CPNums[1];
						GP.CPNum = GP.SaddleCPNum;
						GP.IsSGP = true;
						GP.CPType = SaddleCPTermType;
						GP.SaddleEndGPs.resize(2);

						vec3 SaddleCPPos = AllCPs.GetXYZ(GP.SaddleCPNum);

						DoPrintPaths = DoPrintPaths && DebugSaveGPWaitForUser(GP.GP, XYZVarNums, RhoVarNum, DebugPathStringBase + "\nTerminal CP: " + to_string(GP.CPNum) + ", type: " + to_string(GP.CPType) + "\nConverged at saddle", PassNum, SubPassNum, 1);

						for (int i = 0; i < 2; ++i) {
							int ZoneNum;
							auto TmpPath = (i == 0 ? Path1 : Path2);
							AllCPs.HasEdge(TmpPath->CPNum, GP.SaddleCPNum, &ZoneNum);
							GP.SaddleEndGPs[i] = GradPath_c(ZoneNum, XYZRhoVarNums, AddOnID);

							if (DistSqr(GP.SaddleEndGPs[i][0], SaddleCPPos) > DistSqr(GP.SaddleEndGPs[i][-1], SaddleCPPos))
								// Need to reverse GP
								GP.SaddleEndGPs[i].Reverse();
// 							if (CPType == CPType_Ring && GP.SaddleEndGPs[i].RhoAt(0) > GP.SaddleEndGPs[i].RhoAt(-1))
// 								// Need to reverse GP
// 								GP.SaddleEndGPs[i].Reverse();
// 								
						}

						// Now check that the end GPs terminate at the terminuses of Path1
						// and Path2 respectively.
						if (DistSqr(GP.SaddleEndGPs[0].XYZAt(-1), Path1->GP.XYZAt(-1)) > DistSqr(GP.SaddleEndGPs[0].XYZAt(-1), Path2->GP.XYZAt(-1)))
							// Need to swap end GPs because they're going in the wrong directions.
							GP.SaddleEndGPs = vector<GradPath_c>(GP.SaddleEndGPs.rbegin(), GP.SaddleEndGPs.rend());

						DoPrintPaths = DoPrintPaths && DebugSaveGPWaitForUser(GP.SaddleEndGPs[0], XYZVarNums, RhoVarNum, DebugPathStringBase + "\nTerminal CP: " + to_string(GP.CPNum) + ", type: " + to_string(GP.CPType) + "\nConverged at saddle: saddle GP 1", PassNum, SubPassNum, 2);
						DoPrintPaths = DoPrintPaths && DebugSaveGPWaitForUser(GP.SaddleEndGPs[1], XYZVarNums, RhoVarNum, DebugPathStringBase + "\nTerminal CP: " + to_string(GP.CPNum) + ", type: " + to_string(GP.CPType) + "\nConverged at saddle: saddle GP 2", PassNum, SubPassNum, 3);

						// Get seed point from saddle CP, which is normal to the CPs principal curvature direction
						vec3 v1 = AllCPs.GetPrincDir(GP.SaddleCPNum);
						vec3 v2 = ((Path1->GP.ClosestPoint(SaddleCPPos) + Path2->GP.ClosestPoint(SaddleCPPos)) * 0.5) - SaddleCPPos;
						vec3 v3 = cross(v1, v2);
						v2 = normalise(cross(v1, v3));
						vec3 SeedPt = SaddleCPPos + v2 * 0.01;
						if (DistSqr(SaddleCPPos, CPPosition) < DistSqr(SeedPt, CPPosition))
							SeedPt = SaddleCPPos - v2 * 0.01;

						// Make surface from the SaddleEndGPs Path1, and Path2
						auto SaddlePath = GP.SaddleEndGPs[0] + GP.SaddleEndGPs[1];
						// 						if (DistSqr(SaddlePath[0], Path2->GP[0]) > DistSqr(SaddlePath[0], Path1->GP[0]))
						// 							SaddlePath.Reverse();
						int SaddlePathCPInd;
						SaddlePath.ClosestPoint(SaddleCPPos, SaddlePathCPInd);
						FESurface_c TmpSurf;
						TmpSurf.MakeFromGPs(vector<GradPath_c const *>({ &Path1->GP, &Path2->GP }), false, false, false, SaddlePath.GetXYZListPtr(), Path1Ind, Path2Ind, SaddlePathCPInd);
// 						TmpSurf.SaveAsTriFEZone(XYZVarNums, "test surface");
						TmpSurf.GeneratePointElementDistanceCheckData();

						// Seed GP back to the surface origin CP
						GradPath_c GPSaddle(SeedPt, GPDir == StreamDir_Forward ? StreamDir_Reverse : StreamDir_Forward, NumGPPts, GPType_Classic, GPTerminate_AtPoint, &CPPosition, &AllCPs, &TermRadius, &RhoCutoff, VolInfo, HessPtrs, GradPtrs, RhoPtr, &TmpSurf);
						GPSaddle.Seed(false);
						GPSaddle.PointPrepend(SaddleCPPos, AllCPs.GetRho(GP.SaddleCPNum));
						GPSaddle.SetStartEndCPNum(GP.SaddleCPNum, 0);
						if (sum(GPSaddle[-1] == CPPosition) != 3) {
							GPSaddle.PointAppend(CPPosition, CPRho);
							GPSaddle.SetStartEndCPNum(CurrentCPNum, 1);
						}
						GPSaddle.ReinterpolateRhoValuesFromVolume(&VolInfo, &RhoPtr);
						GPSaddle.RemoveKinks();
						if (DistSqr(GPSaddle[0], CPPosition) > DistSqr(GPSaddle[-1], CPPosition))
							// Need to reverse GP
							GPSaddle.Reverse();
// 						if ((CPType == CPType_Ring && GPSaddle.RhoAt(0) > GPSaddle.RhoAt(-1)) || (CPType == CPType_Bond && GPSaddle.RhoAt(0) < GPSaddle.RhoAt(-1)))
// 							GPSaddle.Reverse();
						GP.GP = GPSaddle;

						GPList.insert(Path2, GP);
						Path1 = Path2;
						Path2++;

						// Midpoint has converged and there is a known saddle CP connecting
						// the terminal CPs of Path1 and Path2.
						// Start a binary search between the Points at Path1Ind and Path2Ind
						// until the saddle CP is found.

// 						int Iter = 0;
// 						int MaxIter = 50;
// 						double Low = 0.0, Mid = 0.5, High = 1.0;
// 						if (GP.CPNum == Path1->CPNum)
// 							Low = Mid;
// 						else if (GP.CPNum == Path2->CPNum)
// 							High = Mid;
// 
// 						while (GP.CPType != SaddleCPTermType && Iter++ < MaxIter){
// 							Mid = (Low + High) * 0.5;
// 							GPMidPt = Path1->GP.XYZAt(Path1Ind) * (1.0 - Mid) + Path2->GP.XYZAt(Path2Ind) * Mid;
// 							GP.GP.Clear();
// 							GP.GP.SetStartPoint(GPMidPt);
// 							GP.GP.Seed(false);
// 							GP.CPNum = GP.GP.GetStartEndCPNum(EndCPPosition);
// 							GP.CPType = AllCPs.GetTypeFromTotOffset(GP.CPNum);
// 							if (GP.CPType != SaddleCPTermType){
// 								if (GP.CPNum == Path1->CPNum)
// 									Low = Mid;
// 								else if (GP.CPNum == Path2->CPNum)
// 									High = Mid;
// 								else
// 									TecUtilDialogErrMsg("Close range binary saddle search GP terminated at invalid CP");
// 							}
// 							else{
// 								GP.Angle = Path1->Angle * (1.0 - Mid) + Path2->Angle * Mid;
// 							}
// 						}
// 
// 						if (Iter >= MaxIter) {
// 							// Close range binary search for saddle CP failed, but we know that one
// 							// both exists and that we're as close as were going to get to it using
// 							// the deviation checking method, so we'll instead stitch the existing
// 							// bond/ring path together with Path1 and Path2 to form a surface that
// 							// can be used to seed a GP from the saddle CP back into the central CP.
// 							// So, gather the two bond/ring path segments, then make the surface, seed
// 							// back to the central CP, save to GPList, and move on in the search.
// 							
// 							TecUtilDialogErrMsg("Close range binary saddle search GP failed");
// 							GPMidPt = (Path1->GP.XYZAt(Path1Ind) + Path2->GP.XYZAt(Path2Ind)) * 0.5;
// 							GP.GP.Clear();
// 							GP.GP.SetStartPoint(GPMidPt);
// 							GP.GP.Seed(false);
// 							GP.CPNum = GP.GP.GetStartEndCPNum(EndCPPosition);
// 							GP.CPType = AllCPs.GetTypeFromTotOffset(GP.CPNum);
// 						}

					}
					// The GP could have terminated at either of the terminal CPs or at a third
					// CP, which could be a nuclear/cage or bond/ring CP.
					// If it terminates at one of the current terminal CPs then adjust the Path1
					// and Path2 iterators to divide again.
					// If the GP terminates at a third CP that does not share an edge with both of the terminal
					// CPs in the CP graph, then we continue the process.
					// If it terminates at a bond/ring CP that shares an edge with both terminal CPs,
					// then start a 1-d minimization to get the minimum length GP as a function of the
					// weight used to place the seed point between the two deviating GPs.
					else if (GP.CPType == SaddleCPTermType) {
// 						if (DebugMode) {
// 							TecUtilDialogMessageBox("Finding Saddle-saddle paths: deviating: term at saddle", MessageBoxType_Information);
// 						}
						// New path terminates at a saddle.
						// 
						// Before we're satisfied with this saddle-saddle path, want to make sure it wasn't too 
						// far of a long-shot, that is, that we're closer to the terminal saddle point than to the originating
						// saddle point.
						// If we're closer to the originating than the terminal saddle point, by some factor, then redo this 
						// such that it doesn't terminate at the saddle point.
						// Take the distance between the path seed point and the originating and terminating saddle points to
						// make the decision.
						// Just using a straight-line distance here rather than partial path lengths.
						// 
						double DistStartToSeed = Distance(GPMidPt, CPPosition),
							DistSeedToEnd = Distance(GPMidPt, GP.GP[-1]);

						double DistSeedToStartEndRatioCutoff = 0.2; // this is the maximum value allowed of (DistSeedToEnd / DistStartToSeed), so if set to 0.25, then DistStartToSeed must be 4 times larger than DistSeedToEnd, that is, we must be much closer to the terminating saddle than to the originating one.

						double DistSeedToStartEndRatio;

						if (DistStartToSeed > 1e-4){
							DistSeedToStartEndRatio = DistSeedToEnd / DistStartToSeed;
						}
						else {
							DistSeedToStartEndRatio = DistSeedToStartEndRatioCutoff + 1.;
						}

						if (DistSeedToStartEndRatio > DistSeedToStartEndRatioCutoff
							|| FoundSaddleCPNums.count(GP.CPNum)) {
							// Cutoff is exeeded, so redo the path so that it won't terminate at the saddle
							GP.GP.Clear();
							GP.GP.SetupGradPath(GPMidPt, GPDir, NumGPPts, GPType_Classic, GPTerminate_AtCP, nullptr, &AllCPs, &TermRadius, &RhoCutoff, VolInfo, HessPtrs, GradPtrs, RhoPtr);
							GP.GP.SetTerminalCPTypeNum(GPTerminalCPTypeNum);
							GP.GP.SetStartEndCPNum(CurrentCPNum, 0);
							GP.GP.Seed(false);

							//DEBUG
		// 					GP.GP.SaveAsOrderedZone();

							GP.CPNum = GP.GP.GetStartEndCPNum(EndCPPosition);
							GP.CPType = AllCPs.GetTypeFromTotOffset(GP.CPNum);

							Reiterate = true;
						}
						else{


							// We'll find the shortest length GP by starting a 1-D minimization of path length
							// as a function of the weight 'b' used to position the seed point between the points
							// on Path1 and Path2.
							// First need to find a lower and upper bound on 'b' for the minimization.
							// Do a binary search on either side of the min to get larger values;
							// 

							GradPath_c GPSaddle(vec3(), GPDir, NumGPPts, GPType_Classic, GPTerminate_AtCP, nullptr, &AllCPs, &TermRadius, &RhoCutoff, VolInfo, HessPtrs, GradPtrs, RhoPtr);
							GPSaddle.SetTerminalCPTypeNum(GPTerminalCPSaddleTypeNum);
							GPSaddle.SetStartEndCPNum(CurrentCPNum, 0);
							int Iter = 0, MaxIter = 10;
							int OuterIter = 0, MaxOuterIter = 10;
							double CheckLen = GP.GP.GetLength();
							vec3 Pts[2];
							double Lengths[] = { -1.0,-1.0 };
							while (OuterIter++ < MaxOuterIter && (Lengths[0] < CheckLen || Lengths[1] < CheckLen)) {
								for (int i = 0; i < 2; ++i) {
									vec3 LeftPt = (i == 0 ? Path1->GP.XYZAt(Path1Ind) : Path2Pt);
									vec3 RightPt = GPMidPt;
									double tmpLength = CheckLen - 1.0;
									Iter = 0;
									do {
										Pts[i] = (LeftPt + RightPt) * 0.5;
										GPSaddle.Clear();
										GPSaddle.SetStartPoint(Pts[i]);
										GPSaddle.Seed(false);
										if (GPSaddle.GetStartEndCPNum(EndCPPosition) != GP.CPNum) {
											LeftPt = Pts[i];
											tmpLength = CheckLen - 1.0;
										}
										else if (GPSaddle.GetLength() < CheckLen) {
											RightPt = Pts[i];
										}
										else {
											tmpLength = GPSaddle.GetLength();
										}
										Lengths[i] = tmpLength;
									} while (Iter++ < MaxIter && Lengths[i] < CheckLen);
									if (Iter >= MaxIter)
										break;
								}
								GPSaddle.Clear();
								GPSaddle.SetStartPoint((Pts[0] + Pts[1]) * 0.5);
								GPSaddle.Seed(false);
								CheckLen = GPSaddle.GetLength();
							}
							std::pair<vec3 const *, vec3 const *> StartPoints(&Pts[0], &Pts[1]);
							double bMin = 0.5;
							if (Iter >= MaxIter || !MinLengthGPBetweenPointsCheckStartPoints(GPSaddle, StartPoints, bMin)) {
								// Failed to get a lower and upper bound for the minimization, so remake the GP so that 
								// it won't terminate at the saddle CP, then we can iterate and try again.
								GP.GP.Clear();
								GP.GP.SetStartPoint(GPMidPt);
								GP.GP.SetTerminalCPTypeNum(GPTerminalCPTypeNum);
								GP.GP.Seed(false);

								GP.CPNum = GP.GP.GetStartEndCPNum(EndCPPosition);
								GP.CPType = AllCPs.GetTypeFromTotOffset(GP.CPNum);
								// 
								DoPrintPaths = DoPrintPaths && DebugSaveGPWaitForUser(GP.GP, XYZVarNums, RhoVarNum, DebugPathStringBase + "\nTerminal CP: " + to_string(GP.CPNum) + ", type: " + to_string(GP.CPType) + "\n Failed to get min-length saddle-saddle path", PassNum, SubPassNum, 1);

								Reiterate = true;
							}
							else {
								// 							FESurface_c TmpSurf;
								// 							vector<GradPath_c> SubGPs = { Path1->GP.SubGP(0, Path1->GP.GetIndAtLength(Path1->GP.GetLength(Path1Ind) + Path1->GP.GetLength() * NeighborGPSurfaceIndLengthOffsetFactor)), Path2->GP.SubGP(0,  Path2->GP.GetIndAtLength(Path2->GP.GetLength(Path2Ind) + Path2->GP.GetLength() * NeighborGPSurfaceIndLengthOffsetFactor)) };
								// 							vector<GradPath_c const *> GPPtrs = { &SubGPs[0], &SubGPs[1] };
								// 							TmpSurf.MakeFromGPs(GPPtrs);
								// 							TmpSurf.GeneratePointElementDistanceCheckData();
								// 							GPSaddle.SetSurfPtr(&TmpSurf);
								// 							GPSaddle.SetTermPoint(CPPosition);

								if (MinLengthGPBetweenPoints(GPSaddle, StartPoints, bMin)) {
									FoundSaddleCPNums.insert(GP.CPNum);
									// Minimum length GP found.
									// Combine with midpoint GP from Path1 and Path2 and put into GPList.
									vec3 bMinPt = Pts[0] * bMin + Pts[1] * (1.0 - bMin);
									//  								MidPtRhoVal = ValAtPointByPtr(bMinPt, VolInfo, RhoPtr);
									 // 								GP.GP = GradPath_c({ &Path1->GP,&Path2->GP }, bMin, MidPtRhoVal, -1);
									 // 								GP.GP = GradPath_c({ &Path1->GP,&Path2->GP }, bMin, Path1->GP.RhoAt(Path1Ind), -1);
									 // 								GP.GP.ReinterpolateRhoValuesFromVolume(&VolInfo, &RhoPtr);
									// 								GP.GP = GradPath_c({ &Path1->GP,&Path2->GP }, bMinPt, (CPType == CPType_Ring ? StreamDir_Reverse : StreamDir_Forward), { std::make_pair(0, MIN(Path1Ind + NeighborGPSurfaceIndOffset, Path1->GP.GetCount() - 1)), std::make_pair(0, MIN(Path2Ind + NeighborGPSurfaceIndOffset, Path2->GP.GetCount() - 1)) }, &CPPosition);
									GP.GP = GradPath_c({ &Path1->GP,&Path2->GP }, bMinPt, (CPType == CPType_Ring ? StreamDir_Reverse : StreamDir_Forward), { std::make_pair(0, MAX(Path1->GP.GetIndAtLength(Path1->GP.GetLength(Path1Ind) + Path1->GP.GetLength() * NeighborGPSurfaceIndLengthOffsetFactor), Path1Ind + NeighborGPSurfaceIndOffset)), std::make_pair(0, MAX(Path2->GP.GetIndAtLength(Path2->GP.GetLength(Path2Ind) + Path2->GP.GetLength() * NeighborGPSurfaceIndLengthOffsetFactor), Path2Ind + NeighborGPSurfaceIndOffset)) }, &CPPosition);
									// Need to trim points off GP.GP to make sure that the resulting rho values
									// when combined with GPSaddle are monotonic down the path.
	// 								GPSaddle.Resample(NumGPPts, GPResampleMethod_Linear);
									int tmpInd = 0;
									while ((CPType == CPType_Ring && GP.GP.RhoAt(tmpInd) > GPSaddle.RhoAt(0)) || (CPType == CPType_Bond && GP.GP.RhoAt(tmpInd) < GPSaddle.RhoAt(0)))
										tmpInd++;

									if (tmpInd > 0)
										GP.GP = GP.GP.SubGP(tmpInd, -1);
									GP.GP += GPSaddle;
									// 								GP.GP = GPSaddle;

									GP.GP.ReinterpolateRhoValuesFromVolume(&VolInfo, &RhoPtr);

									GP.GP.Resample(NumGPPts + 5);
									GP.GP.RemoveKinks();

									if (DistSqr(GP.GP[0], CPPosition) > DistSqr(GP.GP[-1], CPPosition))
										// Need to reverse GP
										GP.GP.Reverse();

									// 								if ((CPType == CPType_Ring && GP.GP.RhoAt(0) > GP.GP.RhoAt(-1))
									// 									|| (CPType == CPType_Bond && GP.GP.RhoAt(0) < GP.GP.RhoAt(-1))) {
									// 									GP.GP.Reverse();
									// 								}
									GP.Identifier = Path1->Identifier * bMin + Path2->Identifier * (1.0 - bMin);
									GP.SaddleCPNum = GP.GP.GetStartEndCPNum(EndCPPosition);
									GP.CPType = AllCPs.GetTypeFromTotOffset(GP.SaddleCPNum);
									GP.IsSGP = true;


									DoPrintPaths = DoPrintPaths && DebugSaveGPWaitForUser(GP.GP, XYZVarNums, RhoVarNum, DebugPathStringBase + "\nTerminal CP: " + to_string(GP.CPNum) + ", type: " + to_string(GP.CPType) + "\nTerminated at saddle", PassNum, SubPassNum, 2);

									// Find the end GPs to combine with the saddle GP at a later step.
									// The saddle-min/max CP paths could go to the terminal CPs of Path1 and Path2,
									// or they may go to another CP, in which case there are more CPs to be found
									// in this loop.
									GP.SaddleEndGPs.resize(2);
									int ZoneNums[2];
									vec3 SeedOrigin = AllCPs.GetXYZ(GP.SaddleCPNum);
									int Iter = 1;
									bool IsEdgeGP[] = { false,false };

									// EDIT: First, check for existing paths then seed new ones if necessary

									vector<std::list<GPWithInfo>::iterator> TmpPaths({ Path1,Path2 });

									for (int i = 0; i < 2; ++i) {
										if (AllCPs.HasEdge(TmpPaths[i]->CPNum, GP.SaddleCPNum, &ZoneNums[i])) {
											IsEdgeGP[i] = true;
											GP.SaddleEndGPs[i] = GradPath_c(ZoneNums[i], XYZRhoVarNums, AddOnID);
											vec3 SaddleCPPos = AllCPs.GetXYZ(GP.SaddleCPNum);
											if (DistSqr(GP.SaddleEndGPs[i][0], SaddleCPPos) > DistSqr(GP.SaddleEndGPs[i][-1], SaddleCPPos))
												// Need to reverse GP
												GP.SaddleEndGPs[i].Reverse();
											// 											if (CPType == CPType_Ring && GP.SaddleEndGPs[i].RhoAt(0) > GP.SaddleEndGPs[i].RhoAt(-1)) {
											// 												// Need to reverse GP
											// 												GP.SaddleEndGPs[i].Reverse();
											// 											}
											DoPrintPaths = DoPrintPaths && DebugSaveGPWaitForUser(GP.SaddleEndGPs[i], XYZVarNums, RhoVarNum, DebugPathStringBase + "\nTerminal CP: " + to_string(GP.CPNum) + ", type: " + to_string(GP.CPType) + "\nTerminated at saddle: existing saddle GP", PassNum, SubPassNum, 3);
										}
									}

									while (!GP.SaddleEndGPs[0].IsMade() || !GP.SaddleEndGPs[1].IsMade()) {
										vec3 PrinceDir = AllCPs.GetPrincDir(GP.SaddleCPNum);
										bool DoLoop = false;
										for (int i = 0; i < 2; ++i) {
											// 										if (AllCPs.HasEdge(TmpPaths[i]->CPNum, GP.SaddleCPNum, &ZoneNums[i])) {
											// 											IsEdgeGP[i] = true;
											// 											GP.SaddleEndGPs[i] = GradPath_c(ZoneNums[i], XYZRhoVarNums, AddOnID);
											// 											vec3 SaddleCPPos = AllCPs.GetXYZ(GP.SaddleCPNum);
											// 											if (DistSqr(GP.SaddleEndGPs[i][0], SaddleCPPos) > DistSqr(GP.SaddleEndGPs[i][-1], SaddleCPPos))
											// 												// Need to reverse GP
											// 												GP.SaddleEndGPs[i].Reverse();
											// // 											if (CPType == CPType_Ring && GP.SaddleEndGPs[i].RhoAt(0) > GP.SaddleEndGPs[i].RhoAt(-1)) {
											// // 												// Need to reverse GP
											// // 												GP.SaddleEndGPs[i].Reverse();
											// // 											}
											// 											DoPrintPaths = DoPrintPaths && DebugSaveGPWaitForUser(GP.SaddleEndGPs[i], XYZVarNums, RhoVarNum, DebugPathStringBase + "\nTerminal CP: " + to_string(GP.CPNum) + ", type: " + to_string(GP.CPType) + "\nTerminated at saddle: existing saddle GP", PassNum, SubPassNum, 3);
											// 										}
											// 										else 
											if (!GP.SaddleEndGPs[i].IsMade())
											{
												DoLoop = true;
												auto * sGP = &GP.SaddleEndGPs[i];
												vec3 OffsetVec = PrinceDir * 0.01 * (double)Iter;
												vec3 SeedPt = SeedOrigin + OffsetVec;
												if (GP.SaddleEndGPs[(i + 1) % 2].IsMade() && DistSqr(GP.SaddleEndGPs[(i + 1) % 2][-1], SeedPt) < DistSqr(GP.SaddleEndGPs[(i + 1) % 2][-1], SeedOrigin)) {
													SeedPt = SeedOrigin - OffsetVec;
												}

												sGP->SetupGradPath(SeedPt, GPDir, NumGPPts, GPType_Classic, GPTerminate_AtCP, nullptr, &AllCPs, &TermRadius, &RhoCutoff, VolInfo, HessPtrs, GradPtrs, RhoPtr);
												sGP->SetTerminalCPTypeNum(GPTerminalCPTypeNum);
												sGP->SetStartEndCPNum(GP.SaddleCPNum, 0);
												sGP->Seed(false);
												sGP->PointPrepend(SeedOrigin, AllCPs.GetRho(GP.SaddleCPNum));

												if (GP.SaddleEndGPs[(i + 1) % 2].IsMade() && DistSqr(GP.SaddleEndGPs[(i + 1) % 2].XYZAt(-1), sGP->XYZAt(-1)) < TermRadius * TermRadius) {
													// if (i == 1 && DistSqr(GP.SaddleEndGPs[0].XYZAt(-1), sGP->XYZAt(-1)) < TermRadius * TermRadius) {
		// 											if ((i == 1 && DistSqr(GP.SaddleEndGPs[0].XYZAt(-1), sGP->XYZAt(-1)) < TermRadius * TermRadius)
		// 												|| (DistSqr(TmpPaths[i]->GP.XYZAt(-1), sGP->XYZAt(-1)) > DistSqr(TmpPaths[(i+1)%2]->GP.XYZAt(-1), sGP->XYZAt(-1)))) {
													SeedPt = SeedOrigin - OffsetVec;
													sGP->Clear();
													sGP->SetStartPoint(SeedPt);
													sGP->SetStartEndCPNum(GP.SaddleCPNum, 0);
													sGP->Seed(false);
													sGP->PointPrepend(SeedOrigin, AllCPs.GetRho(GP.SaddleCPNum));
												}
												DoPrintPaths = DoPrintPaths && DebugSaveGPWaitForUser(GP.SaddleEndGPs[i], XYZVarNums, RhoVarNum, DebugPathStringBase + "\nTerminal CP: " + to_string(GP.CPNum) + ", type: " + to_string(GP.CPType) + "\nTerminated at saddle: new saddle GP", PassNum, SubPassNum, 4);
											}
											// 										PrinceDir = -PrinceDir;
										}
										if (DoLoop && DistSqr(GP.SaddleEndGPs[0].XYZAt(-1), GP.SaddleEndGPs[1].XYZAt(-1)) < CoincidentPointCheckDistSqr) {
											// Gradpaths went to same terminus, so repeat with larger seed point offset
											for (int i = 0; i < 2; ++i) {
												if (!IsEdgeGP[i])
													GP.SaddleEndGPs[i].Clear();
											}
										}
										Iter++;
									}

									// Now check that the end GPs terminate at the terminuses of Path1
									// and Path2 respectively.
	// 								if (DistSqr(GP.SaddleEndGPs[0].XYZAt(-1), Path1->GP.XYZAt(-1)) > DistSqr(GP.SaddleEndGPs[0].XYZAt(-1), Path2->GP.XYZAt(-1))) {
									if (DistSqr(GP.SaddleEndGPs[0].XYZAt(-1), Path1->GP.XYZAt(-1)) > DistSqr(GP.SaddleEndGPs[1].XYZAt(-1), Path1->GP.XYZAt(-1))) {
										// Need to swap end GPs because they're going in the wrong directions.
										GP.SaddleEndGPs = vector<GradPath_c>(GP.SaddleEndGPs.rbegin(), GP.SaddleEndGPs.rend());
									}

									for (auto & i : GP.SaddleEndGPs)
										i.ReinterpolateRhoValuesFromVolume(&VolInfo, &RhoPtr);

									GP.GP.ReinterpolateRhoValuesFromVolume(&VolInfo, &RhoPtr);

									if (DistSqr(GP.GP[0], CPPosition) > DistSqr(GP.GP[-1], CPPosition))
										// Need to reverse GP
										GP.GP.Reverse();
									// 								if (CPType == CPType_Ring && GP.GP.RhoAt(0) > GP.GP.RhoAt(-1)) {
									// 									// Need to reverse GP
									// 									GP.GP.Reverse();
									// 								}

									GPList.insert(Path2, GP);

									// Update the GPList iterators to continue the search
									// If the first end GP doesn't share its terminus with Path1 then
									// resume between the saddle GP and Path1.
									// Else if the second end GP doesn't share its terminus with Path2
									// then resume between the saddle GP and Path2.
									// Else the saddle end GPs share terminuses with Path1 and Path2 so
									// resume the search to the right of Path2.
									if (DistSqr(GP.SaddleEndGPs[0].XYZAt(-1), Path1->GP.XYZAt(-1)) > CoincidentPointCheckDistSqr) {
										Path2--;
									}
									else if (DistSqr(GP.SaddleEndGPs[1].XYZAt(-1), Path2->GP.XYZAt(-1)) > CoincidentPointCheckDistSqr) {
										Path1++;
									}
									else {
										Path1 = Path2;
										Path2++;
									}
								}
							}
						}
					}
					else {
						Reiterate = true;
					}
					if (Reiterate) {
						// GP doesn't go into a saddle point, in which case we still want to keep the new GP
						// to alleviate the deviation, so insert into GPList and resume loop between Path1 and
						// the new GP.
// 						GP.GP += GradPath_c({ &Path1->GP,&Path2->GP }, 0.5, MidPtRhoVal, -1);

// 						GradPath_c TmpGP({ &Path1->GP,&Path2->GP }, GPMidPt, (CPType == CPType_Ring ? StreamDir_Reverse : StreamDir_Forward), { std::make_pair(0, MIN(Path1Ind + NeighborGPSurfaceIndOffset, Path1->GP.GetCount() - 1)), std::make_pair(0, MIN(Path2Ind + NeighborGPSurfaceIndOffset, Path2->GP.GetCount() - 1)) }, &CPPosition);
// 						if (Path1->CPType != CPType_Nuclear && Path2->CPType != CPType_Nuclear){
// 							TecUtilDialogMessageBox("iterating between two missing CPs: before back GP", MessageBoxType_Information);
// 						}
						GradPath_c TmpGP({ &Path1->GP,&Path2->GP }, GPMidPt, (CPType == CPType_Ring ? StreamDir_Reverse : StreamDir_Forward), { std::make_pair(0, MAX(Path1->GP.GetIndAtLength(Path1->GP.GetLength(Path1Ind) + Path1->GP.GetLength() * NeighborGPSurfaceIndLengthOffsetFactor), Path1Ind + NeighborGPSurfaceIndOffset)), std::make_pair(0, MAX(Path2->GP.GetIndAtLength(Path2->GP.GetLength(Path2Ind) + Path2->GP.GetLength() * NeighborGPSurfaceIndLengthOffsetFactor), Path2Ind + NeighborGPSurfaceIndOffset)) }, &CPPosition);

						bool debugPrint = false;
						if (debugPrint) {
							FESurface_c TmpSurf;
							vector<GradPath_c> SubGPs = { Path1->GP.SubGP(0,MAX(Path1->GP.GetIndAtLength(Path1->GP.GetLength(Path1Ind) + Path1->GP.GetLength() * NeighborGPSurfaceIndLengthOffsetFactor), Path1Ind + NeighborGPSurfaceIndOffset)), Path2->GP.SubGP(0, MAX(Path2->GP.GetIndAtLength(Path2->GP.GetLength(Path2Ind) + Path2->GP.GetLength() * NeighborGPSurfaceIndLengthOffsetFactor), Path2Ind + NeighborGPSurfaceIndOffset)) };
							// 							for (auto & i : SubGPs)
							// 								i.Reverse();
							TmpSurf.MakeFromGPs(vector<GradPath_c const *>({ &SubGPs[0], &SubGPs[1] }));
							TmpSurf.SaveAsTriFEZone(XYZVarNums, "back GP surface");
						}

// 						if (Path1->CPType != CPType_Nuclear && Path2->CPType != CPType_Nuclear) {
// 							TecUtilDialogMessageBox("iterating between two missing CPs: before forward GP", MessageBoxType_Information);
// 						}
						GP.GP.SetTerminalCPTypeNum(GPTerminalCPTypeNum);
						GP.GP.SetTermPointRadius(TermRadius);
						GP.GP.Clear();
						GP.GP.Seed(false);

						GP.CPNum = GP.GP.GetStartEndCPNum(EndCPPosition);
						GP.CPType = AllCPs.GetTypeFromTotOffset(GP.CPNum);

// 						if (Path1->CPType != CPType_Nuclear && Path2->CPType != CPType_Nuclear) {
// 							TecUtilDialogMessageBox("iterating between two missing CPs: before combining GPs", MessageBoxType_Information);
// 						}

						int tmpInd = 0;
						while (tmpInd < GP.GP.GetCount() && (CPType == CPType_Ring && GP.GP.RhoAt(tmpInd) < TmpGP.RhoAt(0)) || (CPType == CPType_Bond && GP.GP.RhoAt(tmpInd) > TmpGP.RhoAt(0)))
							tmpInd++;
						if (tmpInd > 0 && tmpInd < GP.GP.GetCount())
							GP.GP = GP.GP.SubGP(tmpInd, -1);
						GP.GP += TmpGP;

// 						if (Path1->CPType != CPType_Nuclear && Path2->CPType != CPType_Nuclear) {
// 							TecUtilDialogMessageBox("iterating between two missing CPs: before treating GP", MessageBoxType_Information);
// 						}

						GP.Identifier = (Path1->Identifier + Path2->Identifier) * 0.5;
						GP.GP.Resample(NumGPPts + 5);
						GP.GP.RemoveKinks();
						GP.GP.Resample(NumGPPts);
						GP.GP.ReinterpolateRhoValuesFromVolume(&VolInfo, &RhoPtr);
						if (DistSqr(GP.GP[0], CPPosition) > DistSqr(GP.GP[-1], CPPosition))
							// Need to reverse GP
							GP.GP.Reverse();
						GPList.insert(Path2, GP);

// 						if (Path1->CPType != CPType_Nuclear && Path2->CPType != CPType_Nuclear) {
// 							TecUtilDialogMessageBox(string("iterating between two missing CPs: before printing GP\nGP made: " + to_string(GP.GP.IsMade())).c_str(), MessageBoxType_Information);
// 						}

						DoPrintPaths = DoPrintPaths && DebugSaveGPWaitForUser(GP.GP, XYZVarNums, RhoVarNum, DebugPathStringBase + "\nTerminal CP: " + to_string(GP.CPNum) + ", type: ", PassNum, SubPassNum, 1);

						// recheck region between new GP and Path1 by setting Path2 to the new GP
						Path2--;
					}
				}
			}
			else {

// 			if (DebugMode) {
// 				TecUtilDialogMessageBox("Finding Saddle-saddle paths: not deviating", MessageBoxType_Information);
// 			}
				Path1++;
				Path2++;
			}
		}

// 		if (DebugMode) {
// 			TecUtilDialogMessageBox("Finished with saddle-saddle paths", MessageBoxType_Information);
// 		}
// 		
// 		return TRUE;

		// DEBUG;
#ifdef DEBUG_PRINTPATHS
		SubPassNum = 0;
		DoPrintPaths = DebugMode;
		if (DoPrintPaths) {
			{
// 				int DebugNumZonesAfter = TecUtilDataSetGetNumZones();
// 				Set DebugDeleteZoneSet;
// 				for (int i = DebugNumZonesBefore + 1; i <= DebugNumZonesAfter; ++i)
// 					DebugDeleteZoneSet += i;
// 				if (!DebugDeleteZoneSet.isEmpty())
// 					TecUtilZoneDelete(DebugDeleteZoneSet.getRef());
			}
			DebugSaveGPList(GPList, ++PassNum, XYZVarNums, RhoVarNum, AllCPs, true, true);
		}
#endif
// 
// 		DoContinue = true;
		if (!DoContinue) break;

		if (DebugMode) {
			TecUtilRedrawAll(TRUE);
			TecUtilDialogMessageBox("Finding terminal paths for Saddle-saddle paths", MessageBoxType_Information);
		}

		vector<GPWithInfo> SGPVector;

		// GPList now contains paths that terminate at saddle points and
		// also has the left and right paths to concatenate with the saddle
		// paths in order to construct distinct, non-overlapping special zero
		// flux surface segments.
		// Do a second pass over GPList to actually perform the concatenation
		// of the saddle paths with there end GPs.
		// (Assume correctness in the end GPs; i.e. that the end GPs for each
		// saddle path do indeed share terminuses with the paths that neighbor
		// the saddle path.)
		Path1 = GPList.begin();
		while (DoContinue && Path1 != GPList.end()) {
			DoContinue = StatusUpdate(iCP * NumSteps + StepNum, StatusTotalNum, CPStatusStr + StepStrings[StepNum], AddOnID, Time1);
			if (Path1->IsSGP) {
#ifdef DEBUG_PRINTPATHS
				SubPassNum = 0;
				DebugPathStringBase = "Path1 id:" + to_string(Path1->Identifier) + ", CP:" + to_string(Path1->CPNum);
#endif
				// Found a saddle path.
				// Add two new paths that are a combination of the saddle path with
				// each of its end GPs that share terminuses with the paths around the
				// saddle path.

				// Point Path2 to the path before the saddle path
				Path2 = Path1;
				Path2--;

				// flip order of Path1's saddle GPs if necessary
				if (DistSqr(Path1->SaddleEndGPs[0].XYZAt(-1), Path2->GP.XYZAt(-1)) > DistSqr(Path1
					->SaddleEndGPs[1].XYZAt(-1), Path2->GP.XYZAt(-1))) {
					// Need to swap end GPs because they're going in the wrong directions.
					Path1->SaddleEndGPs = vector<GradPath_c>(Path1->SaddleEndGPs.rbegin(), Path1->SaddleEndGPs.rend());
				}

// 				for (int i = 0; i < 2; ++i)
// 					Path1->SaddleEndGPs[i].SaveAsOrderedZone("saddle path " + to_string(i));

				// Make new GP from saddle path and the left end GP and copy relevant info
				GPWithInfo GP;
				GP.GP = Path1->GP;
				GP.GP.ConcatenateResample(Path1->SaddleEndGPs[0], NumGPPts, GPResampleMethod_Linear);
				// 				TmpGP.GP = Path1->GP + Path1->SaddleEndGPs[0];
				GP.CPNum = Path2->CPNum;
				GP.SaddleCPNum = Path1->SaddleCPNum;
				GP.Identifier = Path1->Identifier;
				GP.IsSGP = true;
				GP.CPType = Path2->CPType;
				GP.GP.ReinterpolateRhoValuesFromVolume(&VolInfo, &RhoPtr);

				if (DistSqr(GP.GP[0], CPPosition) > DistSqr(GP.GP[-1], CPPosition))
					// Need to reverse GP
					GP.GP.Reverse();
// 				if (CPType == CPType_Ring && GP.GP.RhoAt(0) > GP.GP.RhoAt(-1)) {
// 					// Need to reverse GP
// 					GP.GP.Reverse();
// 				}

				DoPrintPaths = DoPrintPaths && DebugSaveGPWaitForUser(GP.GP, XYZVarNums, RhoVarNum, DebugPathStringBase + "\nLeft concatenated saddle GP", PassNum, SubPassNum, 1);

				// Insert into GPList before Path1
				GPList.insert(Path1, GP);

				// Point Path2 to the path after the saddle path
				Path2 = Path1;
				Path2++;
				// Update TmpGP for the right end GP
				GP.GP = Path1->GP;
				GP.GP.ConcatenateResample(Path1->SaddleEndGPs[1], NumGPPts, GPResampleMethod_Linear);
				// 				TmpGP.GP = Path1->GP + Path1->SaddleEndGPs[1];
				GP.CPNum = Path2->CPNum;
				GP.CPType = Path2->CPType;
				GP.GP.ReinterpolateRhoValuesFromVolume(&VolInfo, &RhoPtr);

				if (DistSqr(GP.GP[0], CPPosition) > DistSqr(GP.GP[-1], CPPosition))
					// Need to reverse GP
					GP.GP.Reverse();
// 				if (CPType == CPType_Ring && GP.GP.RhoAt(0) > GP.GP.RhoAt(-1)) {
// 					// Need to reverse GP
// 					GP.GP.Reverse();
// 				}

				DoPrintPaths = DoPrintPaths && DebugSaveGPWaitForUser(GP.GP, XYZVarNums, RhoVarNum, DebugPathStringBase + "\nRight concatenated saddle GP", PassNum, SubPassNum, 2);
				// Insert into GPList after Path1 (before Path2)
				GPList.insert(Path2, GP);

				// Save then delete the saddle path
				SGPVector.push_back(*Path1);
				GPList.erase(Path1);

				// The erasure invalidates the Path1 iterator, so start it back up at Path2
				// (the path after the saddle path).
				Path1 = Path2;
			}
			else
				Path1++;
		}

		// DEBUG
#ifdef DEBUG_PRINTPATHS
		SubPassNum = 0;
		DoPrintPaths = DebugMode;
		if (DoPrintPaths) {
			{
// 				int DebugNumZonesAfter = TecUtilDataSetGetNumZones();
// 				Set DebugDeleteZoneSet;
// 				for (int i = DebugNumZonesBefore + 1; i <= DebugNumZonesAfter; ++i)
// 					DebugDeleteZoneSet += i;
// 
// 				if (!DebugDeleteZoneSet.isEmpty())
// 					TecUtilZoneDelete(DebugDeleteZoneSet.getRef());
			}
			DebugSaveGPList(GPList, ++PassNum, XYZVarNums, RhoVarNum, AllCPs, true, true);
		}
#endif
//   
//   		DoContinue = true;
		if (!DoContinue) break;

		// Pass to remove GPs that are consistently too close to their neighbors.
		// If the paths to the left and right of a test path are too close, then remove
		// the test path.

// 		Path1 = GPList.begin();
// 		Path2 = Path1;
// 		Path2++;
// 		Path3 = Path2;
// 		Path3++;
// 		CheckDist = 0.01 * MinCPDist;
// 		while (Path2 != GPList.end()) {
// 			if (!Path1->IsSGP && !Path2->IsSGP
// 				&& (Path1->CPNum == Path2->CPNum
// 					|| (Path1->CPType == CPType_Invalid
// 						&& Path2->CPType == CPType_Invalid))
// 				&& Path1->GP.GetMaxSeparationMidpointFromOtherGP(Path2->GP, GPMidPt, Path1Ind, Path2Ind, Path2Pt, MaxDist)
// 				&& MaxDist < CheckDist)
// 			{
// 				// Paths are too close, replace them with a path at their midpoint
// 
// 				GPWithInfo GP;
// 				GP.CPType = Path1->CPType;
// 				GP.CPNum = Path1->CPNum;
// 				GP.Angle = (Path1->Angle + Path2->Angle) * 0.5;
// 				// 				GP.GP = GradPath_c({ &Path1->GP, &Path2->GP }, 0.5, Path1->GP.RhoAt(0), 1);
// 				GP.GP = Path1->GP;
// 				GP.GP.Clear();
// 				GP.GP.SetStartPoint(GPMidPt);
// 				GP.GP.SetDir(GPDir);
// 				GP.GP.SetSurfPtr(nullptr);
// 				GP.GP.SetGPTermType(GPTerminate_AtPoint);
// 				GP.GP.SetTermPoint(AllCPs.GetXYZ(GP.CPNum));
// 				GP.GP.SetTermPointRadius(1e-4);
// 				GP.GP.Seed(false);
// 
// 				GradPath_c TmpGP({ &Path1->GP,&Path2->GP }, GPMidPt, (CPType == CPType_Ring ? StreamDir_Reverse : StreamDir_Forward), { std::make_pair(0, MAX(Path1->GP.GetIndAtLength(Path1->GP.GetLength(Path1Ind) + Path1->GP.GetLength() * NeighborGPSurfaceIndLengthOffsetFactor), Path1Ind + NeighborGPSurfaceIndOffset)), std::make_pair(0, MAX(Path2->GP.GetIndAtLength(Path2->GP.GetLength(Path2Ind) + Path2->GP.GetLength() * NeighborGPSurfaceIndLengthOffsetFactor), Path2Ind + NeighborGPSurfaceIndOffset)) }, &CPPosition);
// 				GP.GP += TmpGP;
// 				GP.GP.RemoveKinks();
// 				GP.GP.Resample(NumGPPts);
// 				GP.GP.ReinterpolateRhoValuesFromVolume(&VolInfo, &RhoPtr);
// 
// 				if (CPType == CPType_Ring && GP.GP.RhoAt(0) > GP.GP.RhoAt(-1)) {
// 					// Need to reverse GP
// 					GP.GP.Reverse();
// 				}
// 
// 				GPList.insert(Path2, GP);
// 
// 				Path3 = Path2;
// 				Path3--;
// 				GPList.erase(Path2);
// 				GPList.erase(Path1);
// 				Path1 = Path3;
// 				Path2 = Path3;
// 				Path2++;
// 			}
// 			else {
// 				Path1++;
// 				Path2++;
// 			}
// 		}
// 		
		if (DebugMode) {
			TecUtilRedrawAll(TRUE);
			TecUtilDialogMessageBox("Removing GPs too close together (GP triplet check)", MessageBoxType_Information);
		}

		Path1 = GPList.begin();
		Path2 = Path1;
		Path2++;
		Path3 = Path2;
		Path3++;
		CheckDist = 0.04 * MinCPDist;
// 		vector<int> IndList;
		IndList.clear();
		while (DoContinue &&  Path1 != GPList.end()) {
			DoContinue = StatusUpdate(iCP * NumSteps + StepNum, StatusTotalNum, CPStatusStr + StepStrings[StepNum], AddOnID, Time1);
			if ((!Path2->IsSGP
				&& (Path1->CPNum == Path2->CPNum
					&& Path1->CPNum == Path3->CPNum)
				|| (Path1->CPType == CPType_Invalid
					&& Path2->CPType == CPType_Invalid
					&& Path3->CPType == CPType_Invalid))
				//  				&& Path1->SaddleCPNum < 0 && Path3->SaddleCPNum < 0
				&& Path2->GP.GetMaxSeparationFromNeighboringGPs(vector<GradPathBase_c const *>({ &Path1->GP, &Path3->GP }), IndList, Path2Ind) < CheckDist)
				//  				&& Path1->GP.GetMaxSeparationMidpointFromOtherGP(Path3->GP, GPMidPt, Path1Ind, Path2Ind, Path2Pt, MaxDist)
				//  				&& MaxDist < CheckDist)
			{
				// Paths are too close, so remove the path between them.

#ifdef DEBUG_PRINTPATHS
				SubPassNum++;
				DebugPathStringBase = "Path1 id:" + to_string(Path1->Identifier) + ", CP:" + to_string(Path1->CPNum) + "\nPath2 id:" + to_string(Path2->Identifier) + ", CP:" + to_string(Path2->CPNum) + "\nPath3 id:" + to_string(Path3->Identifier) + ", CP:" + to_string(Path3->CPNum) + "\nRemoving path 2";
				DoPrintPaths = DoPrintPaths && DebugSaveGPWaitForUser(Path2->GP, XYZVarNums, RhoVarNum, DebugPathStringBase, true, PassNum, SubPassNum, 1);
#endif
				GPList.erase(Path2);
				Path2 = Path3;
				GPIncrement(Path3, GPList);
			}
			else {
				Path1++;
				GPIncrement(Path2, GPList);
				GPIncrement(Path3, GPList);
			}
		}


		// DEBUG
#ifdef DEBUG_PRINTPATHS
		SubPassNum = 0;
		DoPrintPaths = DebugMode;
		if (DoPrintPaths) {
			{
// 				int DebugNumZonesAfter = TecUtilDataSetGetNumZones();
// 				Set DebugDeleteZoneSet;
// 				for (int i = DebugNumZonesBefore + 1; i <= DebugNumZonesAfter; ++i)
// 					DebugDeleteZoneSet += i;
// 
// 				if (!DebugDeleteZoneSet.isEmpty())
// 					TecUtilZoneDelete(DebugDeleteZoneSet.getRef());
			}
			DebugSaveGPList(GPList, ++PassNum, XYZVarNums, RhoVarNum, AllCPs, true, true);
		}
#endif
// 
//  		DoQuit = true;
		if (!DoContinue) break;

		if (DebugMode) {
			TecUtilRedrawAll(TRUE);
			TecUtilDialogMessageBox("Removing GPs too close together (GP neighbor check)", MessageBoxType_Information);
		}

		Path1 = GPList.begin();
		Path2 = Path1;
		Path2++;
		CheckDist = 0.015 * MinCPDist;
		while (DoContinue && Path1 != GPList.end()) {
			DoContinue = StatusUpdate(iCP * NumSteps + StepNum, StatusTotalNum, CPStatusStr + StepStrings[StepNum], AddOnID, Time1);
			if ((!Path2->IsSGP
				&& Path1->CPNum == Path2->CPNum
				|| (Path1->CPType == CPType_Invalid
					&& Path2->CPType == CPType_Invalid))
				  	&& Path1->SaddleCPNum < 0 && Path3->SaddleCPNum < 0
				&& Path1->GP.GetMaxSeparationMidpointFromOtherGP(Path2->GP, GPMidPt, Path1Ind, Path2Ind, Path2Pt, MaxDist)
				&& MaxDist < CheckDist)
			{
#ifdef DEBUG_PRINTPATHS
				SubPassNum++;
				DebugPathStringBase = "Path1 id:" + to_string(Path1->Identifier) + ", CP:" + to_string(Path1->CPNum) + "\nPath2 id:" + to_string(Path2->Identifier) + ", CP:" + to_string(Path2->CPNum) + "\nRemoving path 2";
				DoPrintPaths = DoPrintPaths && DebugSaveGPWaitForUser(Path2->GP, XYZVarNums, RhoVarNum, DebugPathStringBase, true, PassNum, SubPassNum, 1);
#endif
				// Paths are too close, so remove the path between them.
				GPList.erase(Path2);
				Path2 = Path1;
				GPIncrement(Path2, GPList);
			}
			else {
				Path1++;
				GPIncrement(Path2, GPList);
			}
		}


		// DEBUG
#ifdef DEBUG_PRINTPATHS
		SubPassNum = 0;
		DoPrintPaths = DebugMode;
		if (DoPrintPaths) {
			{
// 				int DebugNumZonesAfter = TecUtilDataSetGetNumZones();
// 				Set DebugDeleteZoneSet;
// 				for (int i = DebugNumZonesBefore + 1; i <= DebugNumZonesAfter; ++i)
// 					DebugDeleteZoneSet += i;
// 
// 				if (!DebugDeleteZoneSet.isEmpty())
// 					TecUtilZoneDelete(DebugDeleteZoneSet.getRef());
			}
			DebugSaveGPList(GPList, ++PassNum, XYZVarNums, RhoVarNum, AllCPs, true, true);
		}
#endif
//     
//      		DoContinue = true;
		if (!DoContinue) break;

		// Another pass to add GPs between pairs that experience too high a value
		// of maximum separation.
		Path1 = GPList.begin();
		Path2 = Path1;
		Path2++;
		CheckDist = 0.3 * MinCPDist;

		DoContinue = StatusUpdate(iCP * NumSteps + StepNum, StatusTotalNum, CPStatusStr + StepStrings[StepNum++], AddOnID, Time1);

		// Remake path identifiers so that it can be useful as a way to cap the depth to which new paths are added.
		double tmpDbl = 0.0;
		for (auto & GP : GPList){
			GP.Identifier = tmpDbl;
			tmpDbl += 10.0;
		}

		if (DebugMode) {
			TecUtilRedrawAll(TRUE);
			TecUtilDialogMessageBox("Filling sparse regions", MessageBoxType_Information);
		}
		while (DoContinue && Path1 != GPList.end()) {
			DoContinue = StatusUpdate(iCP * NumSteps + StepNum, StatusTotalNum, CPStatusStr + StepStrings[StepNum], AddOnID, Time1);
			if (
// 				((Path1->CPNum == Path2->CPNum)
// 				|| (Path1->CPType == CPType_Invalid && Path2->CPType == CPType_Invalid))
// 				&& 
				(!Path1->IsSGP || !Path2->IsSGP)
// 				&& Path1->GP.GetSeparationMidpointAtDistFromOtherGP(Path2->GP, GPMidPt, Path1Ind, Path2Ind, Path2Pt, CheckDist, GPDeviationStartPtLengthFactor))
 				&& Path2->Identifier - Path1->Identifier > CoincidentPointCheckAngle
				&& Path1->GP.GetMaxSeparationMidpointFromOtherGP(Path2->GP, GPMidPt, Path1Ind, Path2Ind, Path2Pt, MaxDist)
				// && Path1->GP.GetMaxSeparationMidpointFromOtherGPRhoBased(Path2->GP, GPMidPt, Path1Ind, Path2Ind, Path2Pt, MaxDist, 1)
				&& ((!Path1->IsSGP && !Path2->IsSGP && MaxDist > CheckDist)
					|| MaxDist > CheckDist * 0.8))
			{
				// 				GradPath_c TmpGP({ &Path1->GP,&Path2->GP }, GPMidPt, (CPType == CPType_Ring ? StreamDir_Reverse : StreamDir_Forward), { std::make_pair(0, MIN(Path1Ind + NeighborGPSurfaceIndOffset, Path1->GP.GetCount() - 1)), std::make_pair(0, MIN(Path2Ind + NeighborGPSurfaceIndOffset, Path2->GP.GetCount() - 1)) }, &CPPosition);
				GradPath_c TmpGP({ &Path1->GP,&Path2->GP }, GPMidPt, (CPType == CPType_Ring ? StreamDir_Reverse : StreamDir_Forward), { std::make_pair(0, MAX(Path1->GP.GetIndAtLength(Path1->GP.GetLength(Path1Ind) + Path1->GP.GetLength() * NeighborGPSurfaceIndLengthOffsetFactor), Path1Ind + NeighborGPSurfaceIndOffset)), std::make_pair(0, MAX(Path2->GP.GetIndAtLength(Path2->GP.GetLength(Path2Ind) + Path2->GP.GetLength() * NeighborGPSurfaceIndLengthOffsetFactor), Path2Ind + NeighborGPSurfaceIndOffset)) }, &CPPosition);

				GPWithInfo GP;
				GP.CPType = Path1->CPType;
				GP.CPNum = Path1->CPNum;
				GP.Identifier = (Path1->Identifier + Path2->Identifier) * 0.5;
				// 				GP.GP = GradPath_c({ &Path1->GP, &Path2->GP }, 0.5, Path1->GP.RhoAt(0), 1);
				GP.GP = Path1->GP;
				GP.GP.Clear();
				GP.GP.SetStartPoint(TmpGP[0]);
				GP.GP.SetDir(GPDir);
				GP.GP.SetSurfPtr(nullptr);
				GP.GP.SetTerminalCPTypeNum(GPTerminalCPTypeNum);
				GP.GP.SetGPTermType(GPTerminate_AtRhoValue);
				GP.GP.Seed(false);

				int tmpInd = 0;
				while ((CPType == CPType_Ring && GP.GP.RhoAt(tmpInd) < TmpGP.RhoAt(0)) || (CPType == CPType_Bond && GP.GP.RhoAt(tmpInd) > TmpGP.RhoAt(0)))
					tmpInd++;
				if (tmpInd > 0)
					GP.GP = GP.GP.SubGP(tmpInd, -1);

				GP.GP += TmpGP;
				GP.GP.Resample(NumGPPts + 5);
				GP.GP.RemoveKinks();
				GP.GP.Resample(NumGPPts);
				GP.GP.ReinterpolateRhoValuesFromVolume(&VolInfo, &RhoPtr);

				if (DistSqr(GP.GP[0], CPPosition) > DistSqr(GP.GP[-1], CPPosition))
					// Need to reverse GP
					GP.GP.Reverse();
// 				if (CPType == CPType_Ring && GP.GP.RhoAt(0) > GP.GP.RhoAt(-1)) {
// 					// Need to reverse GP
// 					GP.GP.Reverse();
// 				}
// 				
#ifdef DEBUG_PRINTPATHS
				SubPassNum++;
				DebugPathStringBase = "Path1 id:" + to_string(Path1->Identifier) + ", CP:" + to_string(Path1->CPNum) + "\nPath2 id:" + to_string(Path2->Identifier) + ", CP:" + to_string(Path2->CPNum) + "\nFiller path";
				DoPrintPaths = DoPrintPaths && DebugSaveGPWaitForUser(GP.GP, XYZVarNums, RhoVarNum, DebugPathStringBase, false, PassNum, SubPassNum, 1);
#endif

				GPList.insert(Path2, GP);
				Path2--;
			}
			else {
				Path1++;
				// 				Path2++;
// 				GPIncrement(Path1, GPList);
				GPIncrement(Path2, GPList);
			}
		}


		// DEBUG
#ifdef DEBUG_PRINTPATHS
		SubPassNum = 0;
		DoPrintPaths = DebugMode;
		if (DoPrintPaths) {
			{
// 				int DebugNumZonesAfter = TecUtilDataSetGetNumZones();
// 				Set DebugDeleteZoneSet;
// 				for (int i = DebugNumZonesBefore + 1; i <= DebugNumZonesAfter; ++i)
// 					DebugDeleteZoneSet += i;
// 
// 				if (!DebugDeleteZoneSet.isEmpty())
// 					TecUtilZoneDelete(DebugDeleteZoneSet.getRef());
			}
			DebugSaveGPList(GPList, ++PassNum, XYZVarNums, RhoVarNum, AllCPs, true, true);
		}
#endif
// 
// 		DoQuit = true;
		if (!DoContinue) break;

// 
// 		// Another pass to make sure that there are at least MinNumGPsPerTerminalCP GPs terminating at each
// 		// of the nuclear/cage CPs in the surface.
		int MinNumGPsPerTerminalCP = 3;
		Path1 = GPList.begin();
		int CheckCPNum = Path1->CPNum;
		auto EndPath = GPList.begin();

		DoContinue = StatusUpdate(iCP * NumSteps + StepNum, StatusTotalNum, CPStatusStr + StepStrings[StepNum++], AddOnID, Time1);

		if (DebugMode) {
			TecUtilRedrawAll(TRUE);
			TecUtilDialogMessageBox("Preparing saddle-terminal path search", MessageBoxType_Information);
		}

		while (++Path1 != GPList.end() && Path1->CPNum == CheckCPNum) {}

		EndPath = Path1;
		GPDecrement(EndPath, GPList);

		if (Path1 != GPList.end()) {
			CheckCPNum = Path1->CPNum;
			Path2 = Path1;
			int CurrentCPCount = 0;
			int iter = 0;
			while (false) {
// 				while (DoContinue && Path1 != EndPath) {
				if (iter == 0){
					GPIncrement(EndPath, GPList);
				}
				iter++;
				DoContinue = StatusUpdate(iCP * NumSteps + StepNum, StatusTotalNum, CPStatusStr + StepStrings[StepNum], AddOnID, Time1);
				if (Path1->CPNum == CheckCPNum) {
					CurrentCPCount++;
				}
				else {
					if (CurrentCPCount < MinNumGPsPerTerminalCP) {
						// Add GPs between all existing GPs, then iterate to check to see
						// if CurrentCPCount is high enough.
						FESurface_c TmpSurf;
						TmpSurf.MakeFromGPs(vector<GradPath_c const *>({ &Path1->GP, &Path2->GP }));
// 						TmpSurf.SaveAsTriFEZone(XYZVarNums, "test surface");
						TmpSurf.GeneratePointElementDistanceCheckData();

						auto StartPath = Path2;
						Path3 = Path1;
						Path1 = Path2;
						GPIncrement(Path2, GPList);

						GradPath_c TmpGP = Path1->GP;
						TmpGP.SetDir(GPDir == StreamDir_Forward ? StreamDir_Reverse : StreamDir_Forward);
						TmpGP.SetTermPoint(CPPosition);
						TmpGP.SetGPTermType(GPTerminate_AtPoint);
						TmpGP.SetTermPointRadius(1e-4);
						TmpGP.SetSurfPtr(&TmpSurf);

						while (Path2 != Path3) {
							GPWithInfo GP;
							Path1->GP.GetMaxSeparationMidpointFromOtherGP(Path2->GP, GPMidPt, Path1Ind, Path2Ind, Path2Pt, MaxDist);
							GP.GP.SetupGradPath(GPMidPt, StreamDir_Both, NumGPPts, GPType_Classic, GPTerminate_AtCP, nullptr, &AllCPs, &TermRadius, &RhoCutoff, VolInfo, HessPtrs, GradPtrs, RhoPtr, &TmpSurf);
// 							GP.GP.SetupGradPath(GPMidPt, GPDir, NumGPPts, GPType_Classic, GPTerminate_AtCP, nullptr, &AllCPs, &TermRadius, &RhoCutoff, VolInfo, HessPtrs, GradPtrs, RhoPtr);
							GP.GP.SetTerminalCPTypeNum(GPTerminalCPTypeNum);
							GP.GP.SetStartEndCPNum(CurrentCPNum, 0);
							GP.GP.Seed(false);

							GP.CPNum = GP.GP.GetStartEndCPNum(EndCPPosition);
							GP.CPType = AllCPs.GetTypeFromTotOffset(GP.CPNum);

// 							TmpGP.Clear();
// 							TmpGP.SetStartPoint(GPMidPt);
// 							TmpGP.Seed(false);
// 
// 							int tmpInd = 0;
// 							while ((CPType == CPType_Ring && GP.GP.RhoAt(tmpInd) < TmpGP.RhoAt(0)) || (CPType == CPType_Bond && GP.GP.RhoAt(tmpInd) > TmpGP.RhoAt(0)))
// 								tmpInd++;
// 							if (tmpInd > 0)
// 								GP.GP = GP.GP.SubGP(tmpInd, -1);
// 							GP.GP += TmpGP;
							GP.Identifier = (Path1->Identifier + Path2->Identifier) * 0.5;
							GP.GP.RemoveKinks();
							GP.GP.Resample(NumGPPts);
							GP.GP.ReinterpolateRhoValuesFromVolume(&VolInfo, &RhoPtr);

							if (DistSqr(GP.GP[0], CPPosition) > DistSqr(GP.GP[-1], CPPosition))
								// Need to reverse GP
								GP.GP.Reverse();

// 							if (CPType == CPType_Ring && GP.GP.RhoAt(0) > GP.GP.RhoAt(-1))
// 								// Need to reverse GP
// 								GP.GP.Reverse();
// 								
#ifdef DEBUG_PRINTPATHS
							DebugPathStringBase = "Path1 id:" + to_string(Path1->Identifier) + ", CP:" + to_string(Path1->CPNum) + "\nPath2 id:" + to_string(Path2->Identifier) + ", CP:" + to_string(Path2->CPNum) + "\nFiller path for saddle-min minimization";
							DoPrintPaths = DoPrintPaths && DebugSaveGPWaitForUser(GP.GP, XYZVarNums, RhoVarNum, DebugPathStringBase);
#endif

							GPList.insert(Path2, GP);

							Path1 = Path2;
							GPIncrement(Path2, GPList);
						}
						Path1 = StartPath;
						Path2 = Path1;
					}
					else {
						Path2 = Path1;
					}
					CurrentCPCount = 0;
					CheckCPNum = Path1->CPNum;
				}
				GPIncrement(Path1, GPList);
			}
		}

		// DEBUG
#ifdef DEBUG_PRINTPATHS
		DoPrintPaths = DebugMode;
		if (DoPrintPaths) {
			{
// 				int DebugNumZonesAfter = TecUtilDataSetGetNumZones();
// 				Set DebugDeleteZoneSet;
// 				for (int i = DebugNumZonesBefore + 1; i <= DebugNumZonesAfter; ++i)
// 					DebugDeleteZoneSet += i;
// 
// 				if (!DebugDeleteZoneSet.isEmpty())
// 					TecUtilZoneDelete(DebugDeleteZoneSet.getRef());
			}
			DebugSaveGPList(GPList, ++PassNum, XYZVarNums, RhoVarNum, AllCPs, true, true);
		}
#endif
//  
//  		DoContinue = true;
		if (!DoContinue) break;


		// GPList now contains GPs that go to all CPs in the interatomic/ring surface,
		// Including the minimum length saddle-saddle GPs.
		// Do another pass to find minimum length saddle-min/max GPs.
		// Use three iterators Path1/2/3, and check for when Path2's length is 
		// less than 1 and 3.
		// (Recall that the last GP in GPlist is the same as the first)
		// 

		Path1 = GPList.begin();
		Path2 = Path1;
		Path2++;
		Path3 = Path2;
		Path3++;

		DoContinue = StatusUpdate(iCP * NumSteps + StepNum, StatusTotalNum, CPStatusStr + StepStrings[StepNum++], AddOnID, Time1);

		if (DebugMode) {
			TecUtilRedrawAll(TRUE);
			TecUtilDialogMessageBox("Finding saddle-terminal paths", MessageBoxType_Information);
		}

		std::set<int> UsedCPNums;

		while (DoContinue && Path1 != GPList.end()) {
			DoContinue = StatusUpdate(iCP * NumSteps + StepNum, StatusTotalNum, CPStatusStr + StepStrings[StepNum], AddOnID, Time1);
			// if path2 is a special gradient path, then it means we're at one
			// of the saddle path pairs made in the previous pass, at which we 
			// don't need to check for minimum, so move on.
			if (!Path2->IsSGP
				&& Path2->CPType != CPType_Invalid
				&& Path2->CPNum >= 0
				&& Path1->CPNum == Path2->CPNum
				&& Path2->CPNum == Path3->CPNum
				&& !UsedCPNums.count(Path2->CPNum)
				&& ((Path2->GP.GetLength() < Path1->GP.GetLength()
					&& Path2->GP.GetLength() < Path3->GP.GetLength())
					|| (Path1->SaddleCPNum >= 0 && Path3->SaddleCPNum >= 0))
				)
			{
#ifdef DEBUG_PRINTPATHS
				SubPassNum++;
					DebugPathStringBase = "Path1 id:" + to_string(Path1->Identifier) + ", CP:" + to_string(Path1->CPNum) + ", Length:" + to_string(Path1->GP.GetLength()) + "\nPath2 id:" + to_string(Path2->Identifier) + ", CP:" + to_string(Path2->CPNum) + ", Length:" + to_string(Path2->GP.GetLength()) + "\nPath3 id:" + to_string(Path3->Identifier) + ", CP:" + to_string(Path3->CPNum) + ", Length:" + to_string(Path3->GP.GetLength()) + "\nMinimum length paths to be used as min length starting point";
					DoPrintPaths = DoPrintPaths && DebugSaveGPWaitForUser(Path2->GP, XYZVarNums, RhoVarNum, DebugPathStringBase);
#endif 

				UsedCPNums.insert(Path2->CPNum);

				// Start bypass of finding actual minimum length GP
// 				Path2->IsSGP = true;
// 				SGPVector.push_back(*Path2);
// 				GPIncrement(Path3, GPList);
// 				GPIncrement(Path2, GPList);
// 				Path1++;
// 				continue;
				// end bypass

 				GPWithInfo GP = *Path2;
 				GP.GP.Clear();
 				GP.GP.SetTermPoint(Path2->GP.XYZAt(-1));
 				GP.GP.SetGPTermType(GPTerminate_AtPoint);
 				GP.IsSGP = true;
 
				// Check to see that we have the actual minimum for the terminal CP
				double MinLength = Path2->GP.GetLength();
				auto MinPath = Path2;
				int TmpCPNum = Path2->CPNum;
				while (Path2->CPNum == TmpCPNum && Path2->SaddleCPNum < 0) {
					GPDecrement(Path2, GPList);
				}

				GPIncrement(Path2, GPList);
				while (Path2->CPNum == TmpCPNum && Path2->SaddleCPNum < 0)
				{
					if (Path2->GP.GetLength() < MinLength){
						MinPath = Path2;
						MinLength = Path2->GP.GetLength();
					}
					GPIncrement(Path2, GPList);
				}

				Path2 = MinPath;
				Path2->GP.SetTermPointRadius(SourceCPTermRadius);
				Path1 = Path2;
				GPDecrement(Path1, GPList);
				Path3 = Path2;
				GPIncrement(Path3, GPList);

#ifdef DEBUG_PRINTPATHS
				DebugPathStringBase = "Path1 id:" + to_string(Path1->Identifier) + ", CP:" + to_string(Path1->CPNum) + ", Length:" + to_string(Path1->GP.GetLength()) + "\nPath2 id:" + to_string(Path2->Identifier) + ", CP:" + to_string(Path2->CPNum) + ", Length:" + to_string(Path2->GP.GetLength()) + "\nPath3 id:" + to_string(Path3->Identifier) + ", CP:" + to_string(Path3->CPNum) + ", Length:" + to_string(Path3->GP.GetLength()) + "\nVerified minimum length paths to be used as min length starting point";
				DoPrintPaths = DoPrintPaths && DebugSaveGPWaitForUser(Path2->GP, XYZVarNums, RhoVarNum, DebugPathStringBase, false, PassNum, SubPassNum, 1);
#endif 

				// Begin implementation of simple binary search instead of GSL 1d minimizer
				GradPath_c MinLengthGP = MinLengthGPFromGPTripletWithCommonTerminus({ Path1->GP, Path2->GP, Path3->GP }, 5, -1);
				if (MinLengthGP.IsMade()){
					// We were already close to the minimum length GP, so just replace Path2 with the min length GP
// 					MinLengthGP.Resample(NumGPPts);
					Path2->GP = MinLengthGP;
					Path2->IsSGP = true;
					SGPVector.push_back(*Path2);

#ifdef DEBUG_PRINTPATHS
					DebugPathStringBase = "Path1 id:" + to_string(Path1->Identifier) + ", CP:" + to_string(Path1->CPNum) + ", Length:" + to_string(Path1->GP.GetLength()) + "\nPath2 id:" + to_string(Path2->Identifier) + ", CP:" + to_string(Path2->CPNum) + ", Length:" + to_string(Path2->GP.GetLength()) + "\nPath3 id:" + to_string(Path3->Identifier) + ", CP:" + to_string(Path3->CPNum) + ", Length:" + to_string(Path3->GP.GetLength()) + "\n Minimum length path";
					DoPrintPaths = DoPrintPaths && DebugSaveGPWaitForUser(Path2->GP, XYZVarNums, RhoVarNum, DebugPathStringBase, false, PassNum, SubPassNum, 1);
#endif 
				}
				// End implementation

				// Minimum length Path2 found. Now move Path1 backwards and Path3 forwards
				// to get wider bounds for the minimizer.
// 				int NumBufferGPs = 5;
// 				int LeftNumBufferGPs = 0,
// 					RightNumBufferGPs = 0;
// 				double TmpLength = MinLength;
// 				for (int i = 0; i < NumBufferGPs && Path1->CPNum == Path2->CPNum && Path1->SaddleCPNum < 0; ++i) {
// 					double CheckDist = Path1->GP.GetLength();
// 					if (CheckDist > TmpLength) {
// 						TmpLength = CheckDist;
// 						GPDecrement(Path1, GPList);
// 						LeftNumBufferGPs++;
// 					}
// 					else
// 						break;
// 				}
// 				TmpLength = MinLength;
// 				for (int i = 0; i < NumBufferGPs && Path3->CPNum == Path2->CPNum && Path3->SaddleCPNum < 0; ++i) {
// 					double CheckLength = Path3->GP.GetLength();
// 					if (CheckLength > TmpLength) {
// 						TmpLength = CheckLength;
// 						GPIncrement(Path3, GPList);
// 						RightNumBufferGPs++;
// 					}
// 					else
// 						break;
// 				}
// 				if (LeftNumBufferGPs != RightNumBufferGPs){
// 					int Diff = abs(LeftNumBufferGPs - RightNumBufferGPs);
// 					if (LeftNumBufferGPs > RightNumBufferGPs){
// 						for (int i = 0; i < Diff; ++i) {
// 							GPIncrement(Path1, GPList);
// 						}
// 					}
// 					else{
// 						for (int i = 0; i < Diff; ++i) {
// 							GPDecrement(Path3, GPList);
// 						}
// 					}
// 				}
// 				vector<GradPathBase_c const *> GPPtrs;
// 				auto Path4 = Path1;
// 				while (Path4 != Path3) {
// 					GPIncrement(Path4, GPList);
// 					GPPtrs.push_back(&(Path4++)->GP);
// 				}
// 				GPPtrs.push_back(&Path4->GP);
// 
// #ifdef DEBUG_PRINTPATHS
// 				for (auto GPPtr : GPPtrs) {
// 					DebugPathStringBase = "CP:" + to_string(Path2->CPNum) + "\nPath2 id:" + to_string(Path2->Identifier) + "\nbuffer paths for minimum length saddle-max path search";
// 					DoPrintPaths = DoPrintPaths && DebugSaveGPWaitForUser(*GPPtr, XYZVarNums, RhoVarNum, DebugPathStringBase, false, PassNum, SubPassNum, 2);
// 				}
// #endif
// 					
// // 					vec3 midPt, Pt2;
// 				double bMin;// , MaxDist;
// 				vector<int> IndList;
// 
// 				double MaxDist = Path2->GP.GetMaxSeparationFromNeighboringGPs(vector<GradPathBase_c const *>({ &Path1->GP, &Path3->GP }), IndList, Path2Ind);
// 				int Path1Ind = IndList.front(), 
// 					Path3Ind = IndList.back();
// 
// //  					Path1->GP.GetMaxSeparationMidpointFromOtherGP(Path3->GP, midPt, Path1Ind, Path3Ind, Pt2, MaxDist);
//  
//   				FESurface_c TmpSurf;
// //  					vector<GradPath_c> SubGPs = { Path1->GP.SubGP(0,MAX(Path1->GP.GetIndAtLength(Path1->GP.GetLength(Path1Ind) + Path1->GP.GetLength() * NeighborGPSurfaceIndLengthOffsetFactor), Path1Ind + NeighborGPSurfaceIndOffset)), Path3->GP.SubGP(0, MAX(Path3->GP.GetIndAtLength(Path3->GP.GetLength(Path3Ind) + Path3->GP.GetLength() * NeighborGPSurfaceIndLengthOffsetFactor), Path3Ind + NeighborGPSurfaceIndOffset)) };
// // 				TmpSurf.MakeFromGPs(vector<GradPath_c const *>({ &Path1->GP, &Path3->GP }), false, false, false);
// 				TmpSurf.MakeFromGPs(GPPtrs);
//  				TmpSurf.GeneratePointElementDistanceCheckData();
//  				GP.GP.SetSurfPtr(&TmpSurf);
//  				GP.GP.SetTermPoint(CPPosition);
//  				GP.GP.SetTermPointRadius(1e-4);
//  
//  				vec3 Pt1 = Path1->GP[Path1Ind],
// 					Pt3 = Path3->GP[Path3Ind];
//  				std::pair<vec3 const *, vec3 const *> StartPoints(&Pt1, &Pt3);
// 
// 				// If Path1 (or Path3) contain a saddle point, then we can't use it for the 
// 				// length minimization because the GPs might deviate the wrong way off the 
// 				// saddle CP.
// 				// In that case, we need to adjust StartPoints in order to make sure that 
// 				// the saddle path isn't used and that the StartPoints still correspond
// 				// to GPs of greater length than that of Path2.
// 				MinFuncParams_GPLengthBetweenPoints Params;
// 				Params.GP = &GP.GP;
// 				Params.StartPoints = &StartPoints;
// 				vec3 AB = Pt3 - Pt1;
// 				bMin = 1.0 - (dot(Path2->GP[Path2Ind] - Pt1, AB) / dot(AB, AB));
// // 				double CheckLength = MinFunc_GPLength_BetweenPoints(0.5, reinterpret_cast<void*>(&Params));
// // 				double MinMaxWeights[] = { 0.0, 1.0 };
// // // 				double EndLengths[] = { 
// // // 					Path1->SaddleCPNum < 0 ? Path1->GP.GetLength() : -1, 
// // // 					Path3->SaddleCPNum < 0 ? Path3->GP.GetLength() : -1 };
// // 				double EndLengths[] = {
// // 					Path1->SaddleCPNum < 0 ? MinFunc_GPLength_BetweenPoints(MinMaxWeights[1], reinterpret_cast<void*>(&Params)) : -1,
// // 					Path3->SaddleCPNum < 0 ? MinFunc_GPLength_BetweenPoints(MinMaxWeights[0], reinterpret_cast<void*>(&Params)) : -1 };
// // 				while (EndLengths[0] < CheckLength || EndLengths[1] < CheckLength)
// // 				{
// // 					for (int i = 0; i < 2; ++i) {
// // 						if (EndLengths[i] < CheckLength) {
// // 							auto TmpPath = Path1;
// // 							double StaticWeight = 0.0, Mid, DynamicWeight = 0.5;
// // 							if (i == 1) {
// // 								TmpPath = Path3;
// // 								StaticWeight = 1.0;
// // 							}
// // 							// Start a binary search to find a point 
// // 							double TmpLength = 0.0;
// // 							do
// // 							{
// // 								Mid = (StaticWeight + DynamicWeight) * 0.5;
// // 								TmpLength = MinFunc_GPLength_BetweenPoints(Mid, reinterpret_cast<void*>(&Params));
// // 
// // 								DynamicWeight = Mid;
// // 							} while (TmpLength < CheckLength);
// // 							MinMaxWeights[i] = DynamicWeight;
// // 							EndLengths[i] = TmpLength;
// // 						}
// // 					}
// // 					vec3 *PtPtrs[] = { &Pt1, &Pt3 };
// // 					for (int i = 0; i < 2; ++i)
// // 						*PtPtrs[i] = Path1->GP[Path1Ind] * (1.0 - MinMaxWeights[i]) + Path3->GP[Path3Ind] * MinMaxWeights[i];
// // 
// // 					CheckLength = MinFunc_GPLength_BetweenPoints(0.5, reinterpret_cast<void*>(&Params));
// // 				}
// 
// // 					MinFuncParams_GPLengthBetweenPoints Params;
// // 					Params.GP = &GP.GP;
// // 					Params.StartPoints = &StartPoints;
// // 
// // 					MinFunc_GPLength_BetweenPoints(0.0, reinterpret_cast<void *>(&Params));
// // 					Params.GP->SaveAsOrderedZone("new path 1", Blue_C);
// // 
// // 					MinFunc_GPLength_BetweenPoints(1.0, reinterpret_cast<void *>(&Params));
// // 					Params.GP->SaveAsOrderedZone("new path 1", Blue_C);
// // 
// // 					MinFunc_GPLength_BetweenPoints(0.5, reinterpret_cast<void *>(&Params));
// // 					Params.GP->SaveAsOrderedZone("new path mid", Blue_C);
// 
//    				if (MinLengthGPBetweenPoints(GP.GP, StartPoints, bMin)) {
//    					// 						GP.GP += GradPath_c({ &Path1->GP, &Path3->GP }, bMin, Path1->GP.RhoAt(Path1Ind), -1);
// //       					GP.GP.Clear();
//       				vec3 SeedPt = Pt1 * bMin + Pt3 * (1.0 - bMin);
// //       					GP.GP.SetStartPoint(SeedPt);
// //       					GP.GP.Seed(false);
// //       				GradPath_c TmpGP({ &Path1->GP,&Path3->GP }, SeedPt, (CPType == CPType_Ring ? StreamDir_Reverse : StreamDir_Forward), { std::make_pair(0, MAX(Path1->GP.GetIndAtLength(Path1->GP.GetLength(Path1Ind) + Path1->GP.GetLength() * NeighborGPSurfaceIndLengthOffsetFactor), Path1Ind + NeighborGPSurfaceIndOffset)), std::make_pair(0, MAX(Path3->GP.GetIndAtLength(Path3->GP.GetLength(Path3Ind) + Path3->GP.GetLength() * NeighborGPSurfaceIndLengthOffsetFactor), Path3Ind + NeighborGPSurfaceIndOffset)) }, &CPPosition);
// 					GradPath_c TmpGP = GP.GP;
// 					TmpGP.Clear();
// 					TmpGP.SetDir(CPType == CPType_Ring ? StreamDir_Reverse : StreamDir_Forward);
// 					TmpGP.SetTermPoint(CPPosition);
// 					TmpGP.SetStartPoint(SeedPt);
// 					TmpGP.Seed(false);
// 					if (sum(TmpGP[-1] == CPPosition) != 3) {
// 						TmpGP.PointAppend(CPPosition, CPRho);
// 						TmpGP.SetStartEndCPNum(CurrentCPNum, 1);
// 					}
// 
//       				GP.GP += TmpGP;
// 					GP.GP.RemoveKinks();
//    					GP.Identifier = Path1->Identifier * bMin + Path3->Identifier * (1.0 - bMin);
// 					GP.GP.Resample(NumGPPts);
// 					GP.GP.ReinterpolateRhoValuesFromVolume(&VolInfo, &RhoPtr);
// 
// 					if (DistSqr(GP.GP[0], CPPosition) > DistSqr(GP.GP[-1], CPPosition))
// 						// Need to reverse GP
// 						GP.GP.Reverse();
// #ifdef DEBUG_PRINTPATHS
// 					DebugPathStringBase = "Path1 id:" + to_string(Path1->Identifier) + ", CP:" + to_string(Path1->CPNum) + "\nPath2 id:" + to_string(Path2->Identifier) + ", CP:" + to_string(Path2->CPNum) + "\nMinimum length saddle-max path";
// 					DoPrintPaths = DoPrintPaths && DebugSaveGPWaitForUser(GP.GP, XYZVarNums, RhoVarNum, DebugPathStringBase, false, PassNum, SubPassNum, 3);
// #endif
// 
// 					*Path2 = GP;
// 					SGPVector.push_back(GP);
// 
// 					Path2 = MinPath;
// 					Path1 = MinPath;
// 					GPDecrement(Path1, GPList);
// 					Path3 = MinPath;
// 					GPIncrement(Path3, GPList);
//    				}
			}
			// Increment the iterators, 
			GPIncrement(Path3, GPList);
			GPIncrement(Path2, GPList);
			Path1++;
		}

		// DEBUG
//        	DebugSaveGPList(GPList, ++PassNum, XYZVarNums, RhoVarNum, AllCPs);
//    
//  		DoQuit = true;
		if (!DoContinue) break;

		int DebugNumZonesAfter = TecUtilDataSetGetNumZones();
		if (DebugNumZonesAfter != DebugNumZonesBefore){
#ifdef DEBUG_PRINTPATHS
			if (!DebugMode || TecUtilDialogMessageBox("Delete debug printed zones?", MessageBoxType_YesNo)) 
			{
				Set DebugDeleteZoneSet;
				for (int i = DebugNumZonesBefore + 1; i <= DebugNumZonesAfter; ++i)
					DebugDeleteZoneSet += i;

				if (!DebugDeleteZoneSet.isEmpty())
					TecUtilZoneDelete(DebugDeleteZoneSet.getRef());
			}
#else
			{
				Set DebugDeleteZoneSet;
				for (int i = DebugNumZonesBefore + 1; i <= DebugNumZonesAfter; ++i)
					DebugDeleteZoneSet += i;

				if (!DebugDeleteZoneSet.isEmpty())
					TecUtilZoneDelete(DebugDeleteZoneSet.getRef());
			}
#endif
		}

		if (DebugMode) {
			TecUtilRedrawAll(TRUE);
			TecUtilDialogMessageBox("Saving SGPs and SZFSs", MessageBoxType_Information);
		}

		// Now resample GPs to user-specified number of points
		NumGPPts /= NumGPPtsWorkingFactor;
		for (auto & GP : GPList){
			GP.GP.Resample(NumGPPts, GP.IsSGP ? GPResampleMethod_Linear : GPResampleMethod_Adaptive);
		}

		// Now the full GPList is made, so make SGPs, special zero flux surface segments
		// and the full interatomic/ring surface.
		// First SGPs
		Set ZoneSet;
		for (auto & SGP : SGPVector){
			auto * GP = &SGP.GP;
			vector<int> StartEndCPNums = GP->GetStartEndCPNum();
			vector<vector<int> > StartEndCPTypeAndOffset(2);
			string Name = "SGP_";
			for (int i = 0; i < 2; ++i) {
				if (StartEndCPNums[i] >= 0) {
					StartEndCPTypeAndOffset[i] = AllCPs.GetTypeNumOffsetFromTotOffset(StartEndCPNums[i]);
					Name += CPNameList[StartEndCPTypeAndOffset[i][0]] + " " + to_string(StartEndCPTypeAndOffset[i][1] + 1);
				}
				else Name += "FF";

				if (i == 0) Name += " to ";
			}

			GP->Resample(NumGPPts, GPResampleMethod_Linear);

			ZoneSet += GP->SaveAsOrderedZone(Name, vector<FieldDataType_e>(), XYZVarNums, RhoVarNum, FALSE, PathColor);
			for (int i = 0; i < 2; ++i) {
				AuxDataZoneSetItem(GP->GetZoneNum(), CSMAuxData.CC.GPEndNumStrs[i], to_string(StartEndCPNums[i] + 1));
				if (StartEndCPNums[i] >= 0) AuxDataZoneSetItem(GP->GetZoneNum(), CSMAuxData.CC.GPEndTypes[i], CPNameList[StartEndCPTypeAndOffset[i][0]]);
				else AuxDataZoneSetItem(GP->GetZoneNum(), CSMAuxData.CC.GPEndTypes[i], "FF");
			}
			AuxDataZoneSetItem(GP->GetZoneNum(), CSMAuxData.CC.ZoneType, CSMAuxData.CC.ZoneTypeGP);
			AuxDataZoneSetItem(GP->GetZoneNum(), CSMAuxData.GBA.SourceZoneNum, to_string(ZoneFinalSourceZoneNum(SelectedCPZoneNum, true)));
			if (CPType == CPType_Ring) {
				if (StartEndCPNums[1] >= 0) {
					if (StartEndCPTypeAndOffset[1][0] == CPTypeNum_Bond)
						AuxDataZoneSetItem(GP->GetZoneNum(), CSMAuxData.CC.ZoneSubType, CSMAuxData.CC.ZoneSubTypeRingBondPath);
					else if (StartEndCPTypeAndOffset[1][0] == CPTypeNum_Nuclear)
						AuxDataZoneSetItem(GP->GetZoneNum(), CSMAuxData.CC.ZoneSubType, CSMAuxData.CC.ZoneSubTypeRingNuclearPath);
				}
			}
			else if (CPType == CPType_Bond) {
				if (StartEndCPNums[1] >= 0) {
					if (StartEndCPTypeAndOffset[1][0] == CPTypeNum_Ring)
						AuxDataZoneSetItem(GP->GetZoneNum(), CSMAuxData.CC.ZoneSubType, CSMAuxData.CC.ZoneSubTypeRingBondPath);
					else if (StartEndCPTypeAndOffset[1][0] == CPTypeNum_Cage)
						AuxDataZoneSetItem(GP->GetZoneNum(), CSMAuxData.CC.ZoneSubType, CSMAuxData.CC.ZoneSubTypeBondCagePath);
				}
			}
		}
		if (!ZoneSet.isEmpty()) {
			TecUtilZoneSetActive(ZoneSet.getRef(), AssignOp_PlusEquals);
		}

		// We'll use two iterators again, Path1 and Path2.
		// Step Path1 through GPList until a saddle-terminal GP pair is found,
		// then step Path2 starting at (Path1 + 1) until a min/max-terminal
		// GP is found.
		vector<FESurface_c> Surfaces;
		FESurface_c WholeSurface;
		Path1 = GPList.begin();
		EndPath = GPList.end();

		string ZoneNameBase = "SZFS_",
			WholeSurfaceZoneNameBase = (CPType == CPType_Bond ? "IAS_" : "RS_"),
			CPNameBase = CPNameList[TypeInd] + " " + to_string(CPNum + 1);

		if (SGPVector.empty()) {
			vector<GradPath_c const *> GPPtrs;
			GPPtrs.reserve(GPList.size());
			for (auto const & p : GPList)
				GPPtrs.push_back(&p.GP);
			GPPtrs.pop_back();
			WholeSurface.MakeFromGPs(GPPtrs, true);
		}
		else{
			while (Path1 != EndPath) {
				// Loop until a SGP is found.
				if (Path1->IsSGP) {
					if (EndPath == GPList.end() && Path1 != GPList.begin())
						EndPath = Path1;
					// Now collect all the GPs between this SGP and the next.
					// The GPs will be combined into a surface.
					Path2 = Path1;
					GPIncrement(Path2, GPList);
					if (!Path2->IsSGP) {
						// Path1 is the beginning of an SZFS.
						// (If this loop is skipped then it's the end of one
						// and we had to move one to the right.)
						Path2 = Path1;
					}
					else {
						Path1 = Path2;
					}
					vector<GradPath_c const *> GPPtrs;
					do
					{
						GPPtrs.push_back(&Path2->GP);
						GPIncrement(Path2, GPList);
					} while (!Path2->IsSGP && Path2 != EndPath);
					GPPtrs.push_back(&Path2->GP);

					vector<int> CPNums = { CurrentCPNum,
						(Path1->SaddleCPNum >= 0 ? Path1->SaddleCPNum : Path1->CPNum),
						(Path2->SaddleCPNum >= 0 ? Path2->SaddleCPNum : Path2->CPNum) };

					// We've got the GPs for the surface.
					// Make the surface.

					Surfaces.push_back(FESurface_c());
					Surfaces.back().MakeFromGPs(GPPtrs);
					WholeSurface += Surfaces.back();

					string SZFSStr = CPNameBase;
					vector<vector<int> > CornerCPTypeOffsets;

					SZFSStr += MakeStringFromCPNums(CPNums, AllCPs, CornerCPTypeOffsets);
					int SurfZoneNum = FESurface_c::SetZoneStyle(Surfaces.back().SaveAsTriFEZone(XYZVarNums, ZoneNameBase + SZFSStr),
						AssignOp_MinusEquals);

					for (int c = 0; c < CornerCPTypeOffsets.size(); ++c) {
						AuxDataZoneSetItem(SurfZoneNum, CSMAuxData.CC.ZFSCornerCPNumStrs[c], to_string(CornerCPTypeOffsets[c][1] + 1));
						if (CornerCPTypeOffsets[c][0] >= 0) AuxDataZoneSetItem(SurfZoneNum, CSMAuxData.CC.ZFSCornerCPTypes[c], CPNameList[CornerCPTypeOffsets[c][0]]);
						else AuxDataZoneSetItem(SurfZoneNum, CSMAuxData.CC.ZFSCornerCPTypes[c], "FF");
					}
					AuxDataZoneSetItem(SurfZoneNum, CSMAuxData.CC.ZoneType, CSMAuxData.CC.ZoneTypeSZFS);
					AuxDataZoneSetItem(SurfZoneNum, CSMAuxData.CC.ZoneSubType, (CPType == CPType_Ring ? CSMAuxData.CC.ZoneSubTypeRSSegment : CSMAuxData.CC.ZoneSubTypeIASSegment));

					// Set Path1 to continue the search.
					Path1 = Path2;
				}
				else {
					GPIncrement(Path1, GPList);
				}
			}
		}

		// Save whole surface
		int SurfZoneNum = FESurface_c::SetZoneStyle(WholeSurface.SaveAsTriFEZone(XYZVarNums, WholeSurfaceZoneNameBase + CPNameBase),
			AssignOp_MinusEquals, TRUE, FALSE);

		AuxDataZoneSetItem(SurfZoneNum, CSMAuxData.CC.ZFSCornerCPNumStrs[0], to_string(AllCPs.GetTotOffsetFromTypeNumOffset(TypeInd, CPNum) + 1));
		AuxDataZoneSetItem(SurfZoneNum, CSMAuxData.CC.ZFSCornerCPTypes[0], CPNameList[TypeInd]);
		AuxDataZoneSetItem(SurfZoneNum, CSMAuxData.CC.ZoneType, CSMAuxData.CC.ZoneTypeSZFS);
		AuxDataZoneSetItem(SurfZoneNum, CSMAuxData.CC.ZoneSubType, (CPType == CPType_Ring ? CSMAuxData.CC.ZoneSubTypeRS : CSMAuxData.CC.ZoneSubTypeIAS));

		// All done with this ring/interatomic surface.
	}

	Set DelVars;
	if (DeleteGradVars) {
		for (int i : GradVarNums)
			DelVars += i;
	}

	if (!DelVars.isEmpty()) {
		TecUtilDataSetDeleteVar(DelVars.getRef());
	}

	StatusDrop(AddOnID);

	//DEBUG
	TecUtilDataLoadEnd();

	CSMGUIUnlock();
	TecUtilLockFinish(AddOnID);
	return TRUE;
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
	Boolean_t IsPeriodic,
	bool PrecalcVars)
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

	if (RCSFuncVarNum > 0 && !RCSFuncPtr.InitializeReadPtr(VolZoneNum, RCSFuncVarNum)){
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

	CSMGUILock();

	vector<CPType_e> CPTypesForMinDistCheck;
	vector<CPTypeNum_e> GPTerminalCPTypeNums;
	CPTypeNum_e GPTerminalSaddleTypeNum = CPTypeNum_Ring;
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
		GPTerminalCPTypeNums = { CPTypeNum_Nuclear, CPTypeNum_Bond };
		GPTerminalSaddleTypeNum = CPTypeNum_Bond;
	}
	else{
		if (AllCPs.NumBonds() <= 1) CPTypesForMinDistCheck.push_back(CPType_Nuclear);
		
		CPTypesForMinDistCheck.push_back(CPType_Bond);

		if (AllCPs.NumRings() > 0) CPTypesForMinDistCheck.push_back(CPType_Ring);

		GPTerminalCPTypeNums = { CPTypeNum_Ring, CPTypeNum_Cage };
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
		GPParams.TermPoint = nullptr;
		GPParams.CPs = &AllCPs;
		GPParams.TermPointRadius = &TermRadius;
		GPParams.TermValue = &RhoCutoff;
		GPParams.VolInfo = &VolInfo;
		GPParams.HessPtrs = &HessPtrs;
		GPParams.GradPtrs = &GradPtrs;
		GPParams.RhoPtr = &RhoPtr;

		StartPointOffset = MIN(0.1 * AllCPs.GetMinCPDist(TypeInd, cpNum), norm(VolInfo.PointSpacingV123));
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

			GPs.push_back(GradPath_c(StartPoint, GPDir, NumGPPts, GPType_Classic, GPTerminate_AtCP, nullptr, &AllCPs, &TermRadius, &RhoCutoff, VolInfo, HessPtrs, GradPtrs, RhoPtr));
			GPs.back().SetTerminalCPTypeNums(GPTerminalCPTypeNums);
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
// 		for (int gpNum = 0; gpNum < NumCircleCheckGPs; ++gpNum) {
// 			GradPath_c GP = GPs[gpNum];
// 			vector<int> StartEndCPNums = GP.GetStartEndCPNum();
// 			vector<vector<int> > StartEndCPTypeAndOffset(2);
// 			string Name = "GP: ";
// 			for (int i = 0; i < 2; ++i) {
// 				if (StartEndCPNums[i] >= 0) {
// 					StartEndCPTypeAndOffset[i] = AllCPs.GetTypeNumOffsetFromTotOffset(StartEndCPNums[i]);
// 					Name += CPNameList[StartEndCPTypeAndOffset[i][0]] + " " + to_string(StartEndCPTypeAndOffset[i][1] + 1);
// 				}
// 				else Name += "FF";
// 
// 				if (i == 0) Name += " to ";
// 			}
// 
// 			GP.SaveAsOrderedZone(Name, vector<FieldDataType_e>(), XYZVarNums, RhoVarNum, FALSE, PathColor);
// 			for (int i = 0; i < 2; ++i) {
// 				AuxDataZoneSetItem(GP.GetZoneNum(), CSMAuxData.CC.GPEndNumStrs[i], to_string(StartEndCPNums[i] + 1));
// 				if (StartEndCPNums[i] >= 0) AuxDataZoneSetItem(GP.GetZoneNum(), CSMAuxData.CC.GPEndTypes[i], CPNameList[StartEndCPTypeAndOffset[i][0]]);
// 				else AuxDataZoneSetItem(GP.GetZoneNum(), CSMAuxData.CC.GPEndTypes[i], "FF");
// 			}
// 			AuxDataZoneSetItem(GP.GetZoneNum(), CSMAuxData.CC.ZoneType, CSMAuxData.CC.ZoneTypeGP);
// 			AuxDataZoneSetItem(GP.GetZoneNum(), "Seed angle", to_string(static_cast<double>(gpNum)* AngleStep).c_str());
// 		}

// 		StatusDrop(AddOnID);
// 		CSMGuiUnlock();
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
						while (TermCPNum < 0 && Iter < 100 && AlphaUpper - AlphaLower > 1e-12){
							Iter++;
							AlphaGuess = (AlphaLower + AlphaUpper) * 0.5;
							StartPoint = GPParams.StartPointOrigin + Rotate(StartVec, AlphaGuess, GPParams.RotAxis);
							GradPath_c GP(StartPoint, GPDir, NumGPPts, GPType_Classic, GPTerminate_AtCP, nullptr, &AllCPs, &TermRadius, &RhoCutoff, VolInfo, HessPtrs, GradPtrs, RhoPtr);
							GP.SetTerminalCPTypeNum(GPTerminalSaddleTypeNum);
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
							while (Iter < 100 && AngleFactor >= 1e-12){
								Iter++;
								RotAngle = AlphaGuess + (AngleStep * AngleFactor * SmallAngleFactor * static_cast<double>(MinGPDir));
								StartPoint = GPParams.StartPointOrigin + Rotate(StartVec, RotAngle, GPParams.RotAxis);
								GradPath_c GP(StartPoint, GPDir, NumGPPts, GPType_Classic, GPTerminate_AtCP, nullptr, &AllCPs, &TermRadius, &RhoCutoff, VolInfo, HessPtrs, GradPtrs, RhoPtr);
								GP.SetTerminalCPTypeNum(GPTerminalSaddleTypeNum);
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
							while (Iter < 100 && AngleFactor >= 1e-12){
								Iter++;
								RotAngle += (AngleStep * AngleFactor * SmallAngleFactor * static_cast<double>(MinGPDir));
								StartPoint = GPParams.StartPointOrigin + Rotate(StartVec, RotAngle, GPParams.RotAxis);
								GradPath_c GP(StartPoint, GPDir, NumGPPts, GPType_Classic, GPTerminate_AtCP, nullptr, &AllCPs, &TermRadius, &RhoCutoff, VolInfo, HessPtrs, GradPtrs, RhoPtr);
								GP.SetTerminalCPTypeNum(GPTerminalSaddleTypeNum);
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
				GradPath_c GP(TmpStartPoints[j], GPDir, NumGPPts, GPType_Classic, GPTerminate_AtCP, nullptr, &AllCPs, &TermRadius, &RhoCutoff, VolInfo, HessPtrs, GradPtrs, RhoPtr);
				GP.SetTerminalCPTypeNums(GPTerminalCPTypeNums);
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
			SGPsPerCP[iCP][i].SetupGradPath(GPParams.StartPointOrigin + Rotate(StartVec, MinFunc_AlphaLowMidHigh[i][1], GPParams.RotAxis), GPDir, NumGPPts, GPType_Classic, GPTerminate_AtCP, nullptr, &AllCPs, &TermRadius, &RhoCutoff, VolInfo, HessPtrs, GradPtrs, RhoPtr);
			SGPsPerCP[iCP][i].SetTerminalCPTypeNums(GPTerminalCPTypeNums);
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
				vec3 PtUp = CPPos + v2 * MIN(0.05 * AllCPs.GetMinCPDist(CPTypesForMinDistCheck), norm(VolInfo.PointSpacingV123));
				vec3 PtDn = CPPos - v2 * MIN(0.05 * AllCPs.GetMinCPDist(CPTypesForMinDistCheck), norm(VolInfo.PointSpacingV123));

				if (dot(v2, pt3 - GPParams.StartPointOrigin) > 0){
					/*
					 *	v2 points in the positive rotation direction
					 */

					SaddleCPSupplementGPs[i][0].SetupGradPath(PtUp, GPDir, NumGPPts,
						GPType_Classic, GPTerminate_AtRhoValue, nullptr,
						&AllCPs, &TermRadius, &RhoCutoff, VolInfo,
						HessPtrs, GradPtrs, RhoPtr);
					SaddleCPSupplementGPs[i][1].SetupGradPath(PtDn, GPDir, NumGPPts,
						GPType_Classic, GPTerminate_AtRhoValue, nullptr,
						&AllCPs, &TermRadius, &RhoCutoff, VolInfo,
						HessPtrs, GradPtrs, RhoPtr);
					for (GradPath_c & GP : SaddleCPSupplementGPs[i]){
						GP.SetTerminalCPTypeNums(GPTerminalCPTypeNums);
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
					GPType_Classic, GPTerminate_AtRhoValue, nullptr,
					&AllCPs, &TermRadius, &RhoCutoff, VolInfo,
					HessPtrs, GradPtrs, RhoPtr);
					SaddleCPSupplementGPs[i][1].SetupGradPath(PtUp, GPDir, NumGPPts,
						GPType_Classic, GPTerminate_AtRhoValue, nullptr,
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
							GradPath_c GP(StartPoint, GPDir, NumGPPts, GPType_Classic, GPTerminate_AtCP, nullptr, &AllCPs, &TermRadius, &RhoCutoff, VolInfo, HessPtrs, GradPtrs, RhoPtr);
							GP.SetTerminalCPTypeNums(GPTerminalCPTypeNums);
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
					GradPath_c GP(StartPoint, GPDir, NumGPPts, GPType_Classic, GPTerminate_AtRhoValue, nullptr, &AllCPs, &TermRadius, &RhoCutoff, VolInfo, HessPtrs, GradPtrs, RhoPtr);
					GP.SetTerminalCPTypeNums(GPTerminalCPTypeNums);
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
					GradPath_c GP(StartPoint, GPDir, NumGPPts, GPType_Classic, GPTerminate_AtRhoValue, nullptr, &AllCPs, &TermRadius, &RhoCutoff, VolInfo, HessPtrs, GradPtrs, RhoPtr);
					GP.SetTerminalCPTypeNums(GPTerminalCPTypeNums);
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
				GradPath_c GP(StartPoint, GPDir, NumGPPts, GPType_Classic, GPTerminate_AtRhoValue, nullptr, &AllCPs, &TermRadius, &RhoCutoff, VolInfo, HessPtrs, GradPtrs, RhoPtr);
				GP.SetTerminalCPTypeNums(GPTerminalCPTypeNums);
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
		Set ZoneSet;
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

			ZoneSet += GP.SaveAsOrderedZone(Name, vector<FieldDataType_e>(), XYZVarNums, RhoVarNum, FALSE, PathColor);
			for (int i = 0; i < 2; ++i){
				AuxDataZoneSetItem(GP.GetZoneNum(), CSMAuxData.CC.GPEndNumStrs[i], to_string(StartEndCPNums[i] + 1));
				if (StartEndCPNums[i] >= 0) AuxDataZoneSetItem(GP.GetZoneNum(), CSMAuxData.CC.GPEndTypes[i], CPNameList[StartEndCPTypeAndOffset[i][0]]);
				else AuxDataZoneSetItem(GP.GetZoneNum(), CSMAuxData.CC.GPEndTypes[i], "FF");
			}
			AuxDataZoneSetItem(GP.GetZoneNum(), CSMAuxData.CC.ZoneType, CSMAuxData.CC.ZoneTypeGP);
			if (CPType == CPType_Ring){
				if (StartEndCPNums[1] >= 0) {
					if (StartEndCPTypeAndOffset[1][0] == CPTypeNum_Bond)
						AuxDataZoneSetItem(GP.GetZoneNum(), CSMAuxData.CC.ZoneSubType, CSMAuxData.CC.ZoneSubTypeRingBondPath);
					else if (StartEndCPTypeAndOffset[1][0] == CPTypeNum_Nuclear)
						AuxDataZoneSetItem(GP.GetZoneNum(), CSMAuxData.CC.ZoneSubType, CSMAuxData.CC.ZoneSubTypeRingNuclearPath);
				}
			}
			else if (CPType == CPType_Bond) {
				if (StartEndCPNums[1] >= 0) {
					if (StartEndCPTypeAndOffset[1][0] == CPTypeNum_Ring)
						AuxDataZoneSetItem(GP.GetZoneNum(), CSMAuxData.CC.ZoneSubType, CSMAuxData.CC.ZoneSubTypeRingBondPath);
					else if (StartEndCPTypeAndOffset[1][0] == CPTypeNum_Cage)
						AuxDataZoneSetItem(GP.GetZoneNum(), CSMAuxData.CC.ZoneSubType, CSMAuxData.CC.ZoneSubTypeBondCagePath);
				}
			}
		}
		if (!ZoneSet.isEmpty()){
			TecUtilZoneSetActive(ZoneSet.getRef(), AssignOp_PlusEquals);
		}

		/*
		 *	Now make the FE surfaces.
		 *	One for each SGP and one big total surface
		 */
		vector<FESurface_c> Surfaces;
		FESurface_c WholeSurface;
		vector<vector<GradPath_c const*> > SurfaceGPPtrs;
		bool SingleSurface = (SurfaceGPs.size() == 1);
		string ZoneNameBase = "SZFS_",
			WholeSurfaceZoneNameBase = (CPType == CPType_Bond ? "IAS_" : "RS_"),
			CPNameBase = CPNameList[TypeInd] + " " + to_string(cpNum + 1);
		if (SingleSurface) ZoneNameBase = WholeSurfaceZoneNameBase;

		for (int i = 0; i < SurfaceGPs.size(); ++i){
			auto * GPVec = &SurfaceGPs[i];
			SurfaceGPPtrs.push_back(vector<GradPath_c const *>());
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
				AssignOp_MinusEquals, TRUE, FALSE);

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
	
	CSMGUIUnlock();
	TecUtilLockFinish(AddOnID);
	return TRUE;


// 	vector<GradPath_c> GPs;
// 
// 	/*
// 	 *	Shotgun approach.
// 	 *	First seed GPs around each CP along a constant radius circle in the 
// 	 *	plane normal to the CP principal direction.
// 	 *	These GPs will identify all bond-cage, ring-nuclear, and bond-ring connections.
// 	 */
// 
// 	StartPointOffset = 0.05 * AllCPs.GetMinCPDist();
// 
// 	for (int cpNum = 0; cpNum < NumCPs; ++cpNum){
// 		vec3 StartVec = normalise(AllCPs.GetEigVecs(TypeInd, cpNum).col(1)) * StartPointOffset;
// 
// 		for (int gpNum = 0; gpNum < NumCircleCheckGPs; ++gpNum){
// 			double RotAngle = static_cast<double>(gpNum)* AngleStep;
// 
// 			StartPoint = AllCPs.GetXYZ(TypeInd, cpNum) + Rotate(StartVec, RotAngle, AllCPs.GetPrincDir(TypeInd, cpNum));
// 
// 			GPs.push_back(GradPath_c(StartPoint, GPDir, NumGPPts, GPType_Classic, GPTerminate_AtCP, nullptr, &AllCPs, &TermRadius, &RhoCutoff, VolInfo, HessPtrs, GradPtrs, RhoPtr));
// 			GPs.back().SetStartEndCPNum(AllCPs.GetTotOffsetFromTypeNumOffset(TypeInd, cpNum), 0);
// 
// 		}
// // 		break;
// 	}
// 
// 	int NumGPs = GPs.size();
// 
// #ifndef _DEBUG
// #pragma omp parallel for schedule(dynamic)
// #endif
// 	for (int iGP = 0; iGP < NumGPs; ++iGP){
// 		GPs[iGP].Seed(false);
// 		GPs[iGP].PointPrepend(AllCPs.GetXYZ(GPs[iGP].GetStartEndCPNum()[0]), AllCPs.GetRho(GPs[iGP].GetStartEndCPNum()[0]));
// 		GPs[iGP].Resample(NumGPPts);
// 		if (CPType == CPType_Bond) GPs[iGP].Reverse();
// 	}
// 
// 	for (auto & GP : GPs){
// 		vector<int> StartEndCPNums = GP.GetStartEndCPNum();
// 		vector<vector<int> > StartEndCPTypeAndOffset(2);
// 		string Name = "GP_";
// 		for (int i = 0; i < 2; ++i){
// 			if (StartEndCPNums[i] >= 0){
// 				StartEndCPTypeAndOffset[i] = AllCPs.GetTypeNumOffsetFromTotOffset(StartEndCPNums[i]);
// 				Name += CPNameList[StartEndCPTypeAndOffset[i][0]] + " " + to_string(StartEndCPTypeAndOffset[i][1] + 1);
// 			}
// 			else Name += "FF";
// 
// 			if (i == 0) Name += " to ";
// 		}
// 
// 		GP.SaveAsOrderedZone(Name, vector<FieldDataType_e>(), XYZVarNums, RhoVarNum, FALSE, PathColor);
// 		for (int i = 0; i < 2; ++i){
// 			AuxDataZoneSetItem(GP.GetZoneNum(), CSMAuxData.CC.GPEndNumStrs[i], to_string(StartEndCPNums[i] + 1));
// 			if (StartEndCPNums[i] >= 0) AuxDataZoneSetItem(GP.GetZoneNum(), CSMAuxData.CC.GPEndTypes[i], CPNameList[StartEndCPTypeAndOffset[i][0]]);
// 			else AuxDataZoneSetItem(GP.GetZoneNum(), CSMAuxData.CC.GPEndTypes[i], "FF");
// 		}
// 		AuxDataZoneSetItem(GP.GetZoneNum(), CSMAuxData.CC.ZoneType, CSMAuxData.CC.ZoneTypeGP);
// 		AuxDataZoneSetItem(GP.GetZoneNum(), CSMAuxData.CC.ZoneSubType, (CPType == CPType_Bond ? CSMAuxData.CC.ZoneSubTypeBondCagePath : CSMAuxData.CC.ZoneSubTypeRingNuclearPath));
// 	}
// 
// 	TecUtilDataLoadEnd();
// 
// 	
// 	CSMGuiUnlock();
// 	TecUtilLockFinish(AddOnID);
// 
// 	return TRUE;
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
			Boolean_t IsOk = TecUtilDataSetAddZone(ZoneName.c_str(), NumPts, 1, 1, ZoneType_Ordered, nullptr);
			FieldVecPointer_c XYZPtr;
			IsOk = IsOk && XYZPtr.InitializeWritePtr(TecUtilDataSetGetNumZones(), XYZVarNums);
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
		if (!RhoPtrs.back().InitializeReadPtr(z, RhoVarNum)){
			TecUtilDialogErrMsg("Failed to get rho data pointer");
			return;
		}
		XYZPtrs.push_back(FieldVecPointer_c());
		if (!XYZPtrs.back().InitializeReadPtr(z, XYZVarNums)){
			TecUtilDialogErrMsg("Failed to get XYZ data pointer");
			return;
		}
		EigValPtrs.push_back(vector<FieldDataPointer_c>());
		for (int v : EigValVarNums){
			EigValPtrs.back().push_back(FieldDataPointer_c());
			if (!EigValPtrs.back().back().InitializeReadPtr(z, v)){
				TecUtilDialogErrMsg("Failed to get eigenvalue data pointer");
				return;
			}
		}
		EigVecPtrs.push_back(vector<FieldVecPointer_c>());
		for (vector<int> const & e : EigVecVarNums) {
			EigVecPtrs.back().push_back(FieldVecPointer_c());
			if (!EigVecPtrs.back().back().InitializeReadPtr(z, e)){
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

	CSMGui("Draw eigenvector arrows", Fields, DrawEigenvectotArrorsReturnUserInfo, AddOnID);
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
		for (int z : ZoneNums) FDTypes[z - 1] = FieldDataType_Float;
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
		FieldDataType_Float);


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

	CSMGUILock();

	// This is where the magic happens!

	CSMGUIUnlock();

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

	VERIFY(SelXYZ.InitializeReadPtr(SelectedCPZoneNum, XYZVarNums) && AllXYZ.InitializeReadPtr(AllCPsZoneNum, XYZVarNums));

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
	vector<bool> const & CalcSteps,
	int NumGPPts)
{
	// Check and run each calculation step. (minus 1 because critical points isn't included)
	if (CalcSteps[(int)BondalyzerCalcType_BondPaths - 1]){
		FindBondRingLines(VolZoneNum, CPZoneNums, CPZoneNums[0], vector<int>(), CPTypeVarNum, CPType_Bond, XYZVarNums, RhoVarNum, GradVarNums, HessVarNums, IsPeriodic, true, NumGPPts);
	}
	if (CalcSteps[(int)BondalyzerCalcType_RingLines - 1]){
		FindBondRingLines(VolZoneNum, CPZoneNums, CPZoneNums[0], vector<int>(), CPTypeVarNum, CPType_Ring, XYZVarNums, RhoVarNum, GradVarNums, HessVarNums, IsPeriodic, true, NumGPPts);
	}
	if (CalcSteps[(int)BondalyzerCalcType_CageNuclearPaths - 1]){
		FindCageNuclearPaths(VolZoneNum, CPZoneNums, CPZoneNums[0], vector<int>(), CPTypeVarNum, XYZVarNums, RhoVarNum, GradVarNums, HessVarNums, IsPeriodic);
	}
	if (CalcSteps[(int)BondalyzerCalcType_InteratomicSurfaces - 1]){
// 		FindBondRingSurfaces(VolZoneNum, CPZoneNums, CPZoneNums[0], vector<int>(), CPTypeVarNum, CPType_Bond, XYZVarNums, RhoVarNum, GradVarNums, HessVarNums, RidgeFuncVarNum, IsPeriodic);
		FindBondRingSurfaces2(VolZoneNum, CPZoneNums, CPZoneNums[0], vector<int>(), CPTypeVarNum, CPType_Bond, XYZVarNums, RhoVarNum, GradVarNums, HessVarNums, RidgeFuncVarNum, IsPeriodic);
	}
	if (CalcSteps[(int)BondalyzerCalcType_RingSurfaces - 1]){
// 		FindBondRingSurfaces(VolZoneNum, CPZoneNums, CPZoneNums[0], vector<int>(), CPTypeVarNum, CPType_Ring, XYZVarNums, RhoVarNum, GradVarNums, HessVarNums, RidgeFuncVarNum, IsPeriodic);
		FindBondRingSurfaces2(VolZoneNum, CPZoneNums, CPZoneNums[0], vector<int>(), CPTypeVarNum, CPType_Ring, XYZVarNums, RhoVarNum, GradVarNums, HessVarNums, RidgeFuncVarNum, IsPeriodic);
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

	if (DeleteSourceZones && !CPZones.isEmpty()) TecUtilZoneDelete(CPZones.getRef());

	CPs.RemoveSpuriousCPs();

	CPs.SaveAsOrderedZone(XYZVarNums, RhoVarNum, TRUE);

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
	int VolZoneNum = Fields[fNum++].GetReturnInt();
	vector<int> ZoneNums = Fields[fNum++].GetReturnIntVec();
	int ResampleNumPoints = Fields[fNum++].GetReturnInt();
	int GroupNum = Fields[fNum++].GetReturnInt();
	bool ConnectFirstAndLastPaths = Fields[fNum++].GetReturnBool();
	bool CapSurface = Fields[fNum++].GetReturnBool();
	double RhoCutoff = Fields[fNum++].GetReturnDouble();
	bool DelSourceZones = Fields[fNum++].GetReturnBool();
	bool PairwiseSurfaces = Fields[fNum++].GetReturnBool();
	int RhoVarNum = Fields[fNum++].GetReturnInt();
	vector<int> XYZVarNums(3);
	for (int i = 0; i < 3; ++i) XYZVarNums[i] = i + Fields[fNum].GetReturnInt();

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

	if (PairwiseSurfaces){
		for (int zi = 0; zi < ZoneNums.size() + (ConnectFirstAndLastPaths ? 0 : -1); ++zi) {
			FESurface_c Surf;

			if (ResampleNumPoints > 0) {
				/*
				 *	First, read in all zones as GPs to resample.
				 */
				
				vector<GradPath_c const *> GPPtrs;
				vector<int> XYZRhoVarNums = XYZVarNums;
				if (RhoVarNum <= 0)
					RhoVarNum = VarNumByName("Electron Density", true);
				if (RhoVarNum > 0)
					XYZRhoVarNums.push_back(RhoVarNum);
				else
					XYZRhoVarNums.push_back(4);
				
				vector<GradPath_c> GPs({ GradPath_c(ZoneNums[zi], XYZRhoVarNums, AddOnID), GradPath_c(ZoneNums[(zi + 1) % ZoneNums.size()], XYZRhoVarNums, AddOnID) });
// 				for (int z : ZoneNums) {
// 					GPs.push_back(GradPath_c(z, XYZRhoVarNums, AddOnID));
// 					if (CapSurface) {
// 						GPs.back().TruncateAtRhoValue(RhoCutoff);
// 					}
// 					GPs.back().Resample(ResampleNumPoints);
// 					if (CapSurface && GPs.back().RhoAt(0) < GPs.back().RhoAt(-1)) {
// 						GPs.back().Reverse();
// 					}
// 				}
				for (auto & gp : GPs)
					GPPtrs.push_back(&gp);

				/*
				 *	Now make and save surface
				 */
				Surf.MakeFromGPs(GPPtrs, false, CapSurface);
			}
			else {
				/*
				 *	First, need to read in all the points of all the zones to a single vector of vec3
				 *	while recording the indices in that vector of each individual path.
				 */
				vector<vec3> P;
				vector<vector<int> > IndList;
				int Ind = 0;
				for (int zj = zi; zj <= zi + 1; ++zj) {
					int z = ZoneNums[zj % ZoneNums.size()];
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
					if (XYZVarPtr.InitializeReadPtr(z, XYZVarNums)) {
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
				Surf = FESurface_c(P, Elem);
				// 		if (Surf.IsMade()) {
				// 			Surf.SaveAsTriFEZone(XYZVarNums, "Stitched path surface");
				// 		}
			}

			if (Surf.IsMade()) {
				int SurfZoneNum = Surf.SaveAsTriFEZone(XYZVarNums, "Stitched path surface");
				if (SurfZoneNum > 0) {
					FieldDataPointer_c RhoPtr;

					VolExtentIndexWeights_s VolInfo;

					if (GetVolInfo(VolZoneNum, XYZVarNums, FALSE, VolInfo)) {
						TecUtilDataLoadBegin();
						RhoPtr.InitializeReadPtr(VolZoneNum, RhoVarNum);
						// Write rho values to each node on the GB

						FieldVecPointer_c SurfXYZPtr;
						SurfXYZPtr.InitializeReadPtr(SurfZoneNum, XYZVarNums);

						FieldDataPointer_c SurfRhoPtr;
						SurfRhoPtr.InitializeWritePtr(SurfZoneNum, RhoVarNum);

						int NumNodes = SurfXYZPtr.Size();
						vector<VolExtentIndexWeights_s> ThVolInfo(omp_get_num_procs(), VolInfo);
#pragma omp parallel for
						for (int i = 0; i < NumNodes; ++i) {
							SurfRhoPtr.Write(i, ValAtPointByPtr(SurfXYZPtr[i], ThVolInfo[omp_get_thread_num()], RhoPtr));
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
						TecUtilZoneSetContour(SV_FLOODCOLORING, SurfZoneSet.getRef(), 0.0, ContourColoring_Group1);
					}
				}
			}
		}
	}
	else {
		FESurface_c Surf;

		if (ResampleNumPoints > 0) {
			/*
			 *	First, read in all zones as GPs to resample.
			 */
			vector<GradPath_c> GPs;
			vector<GradPath_c const *> GPPtrs;
			vector<int> XYZRhoVarNums = XYZVarNums;
			if (RhoVarNum <= 0)
				RhoVarNum = VarNumByName("Electron Density", true);
			if (RhoVarNum > 0)
				XYZRhoVarNums.push_back(RhoVarNum);
			else
				XYZRhoVarNums.push_back(4);
			for (int z : ZoneNums) {
				GPs.push_back(GradPath_c(z, XYZRhoVarNums, AddOnID));
				if (CapSurface) {
					GPs.back().TruncateAtRhoValue(RhoCutoff);
				}
				GPs.back().Resample(ResampleNumPoints);
				if (CapSurface && GPs.back().RhoAt(0) < GPs.back().RhoAt(-1)) {
					GPs.back().Reverse();
				}
			}
			for (auto & gp : GPs)
				GPPtrs.push_back(&gp);

			/*
			 *	Now make and save surface
			 */
			Surf.MakeFromGPs(GPPtrs, ConnectFirstAndLastPaths, CapSurface);
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
				if (XYZVarPtr.InitializeReadPtr(z, XYZVarNums)) {
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
			Surf = FESurface_c(P, Elem);
			// 		if (Surf.IsMade()) {
			// 			Surf.SaveAsTriFEZone(XYZVarNums, "Stitched path surface");
			// 		}
		}

		if (Surf.IsMade()) {
			int SurfZoneNum = Surf.SaveAsTriFEZone(XYZVarNums, "Stitched path surface");
			if (SurfZoneNum > 0) {
				FieldDataPointer_c RhoPtr;

				VolExtentIndexWeights_s VolInfo;

				if (GetVolInfo(VolZoneNum, XYZVarNums, FALSE, VolInfo)) {
					TecUtilDataLoadBegin();
					RhoPtr.InitializeReadPtr(VolZoneNum, RhoVarNum);
					// Write rho values to each node on the GB

					FieldVecPointer_c SurfXYZPtr;
					SurfXYZPtr.InitializeReadPtr(SurfZoneNum, XYZVarNums);

					FieldDataPointer_c SurfRhoPtr;
					SurfRhoPtr.InitializeWritePtr(SurfZoneNum, RhoVarNum);

					int NumNodes = SurfXYZPtr.Size();
					vector<VolExtentIndexWeights_s> ThVolInfo(omp_get_num_procs(), VolInfo);
#pragma omp parallel for
					for (int i = 0; i < NumNodes; ++i) {
						SurfRhoPtr.Write(i, ValAtPointByPtr(SurfXYZPtr[i], ThVolInfo[omp_get_thread_num()], RhoPtr));
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
					TecUtilZoneSetContour(SV_FLOODCOLORING, SurfZoneSet.getRef(), 0.0, ContourColoring_Group1);
				}
			}
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
		GuiField_c(Gui_ZoneSelect, "Volume Zone", "Full Volume"),
		GuiField_c(Gui_ZoneSelectMulti, "Path Zones"),
		GuiField_c(Gui_Int, "Resample with at most # points (0 for no resample)", "500"),
		GuiField_c(Gui_Int, "Group Number", "1"),
		GuiField_c(Gui_Toggle, "Connect first and last paths?"),
		GuiField_c(Gui_Toggle, "Cap surface (if closed surface)?"),
		GuiField_c(Gui_Double, "Rho value at which to cap", "0.001"),
		GuiField_c(Gui_Toggle, "Delete source zones after surface creation?", "1"),
		GuiField_c(Gui_Toggle, "One surface for each pair?", "0"),
		GuiField_c(Gui_VarSelect, "Electron Density", "Electron Density"),
		GuiField_c(Gui_VarSelect, "X", "X")
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
		TecUtilLockFinish(AddOnID);
		return;
	}

	CSMGUILock();

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

	XYZPtr.InitializeReadPtr(CPZoneNum, XYZVarNums);

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

	CSMGUIUnlock();

	TecUtilZoneRename(TecUtilDataSetGetNumZones(), ZoneName.c_str());

	TecUtilSetAddMember(ActiveZones, TecUtilDataSetGetNumZones(), FALSE);
	TecUtilZoneSetActive(ActiveZones, AssignOp_Equals);

	TecUtilSetDealloc(&ActiveZones);

	TecUtilLockFinish(AddOnID);
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



// void TestFunction(){
// 	TecUtilLockStart(AddOnID);
// 	/*
// 	*	Quick test of tecplot toolbox set class
// 	*/
// 	// 	TecUtilZoneSetActive(Set(1).getRef(), AssignOp_Equals);
// 
// 	/*
// 	*	Grad path test for adaptive step size
// 	*/
// // 	vec3 StartPoint;
// // 	StartPoint << 0.1 << 0.1 << 0.1;
// // 
// // 	VolExtentIndexWeights_s VolInfo;
// // 	GetVolInfo(1, { 1, 2, 3 }, FALSE, VolInfo);
// // 
// // 	vector<FieldDataPointer_c> GradPtrs(3), HessPtrs(6);
// // 	FieldDataPointer_c RhoPtr;
// // 
// // 	int VarNum = 4;
// // 
// // 	RhoPtr.InitializeReadPtr(1, VarNum++);
// // 
// // 	for (int i = 0; i < 3; ++i) GradPtrs[i].InitializeReadPtr(1, VarNum++);
// // 	for (int i = 0; i < 6; ++i) HessPtrs[i].InitializeReadPtr(1, VarNum++);
// // 
// // 	double RhoCutoff = 1e-3;
// // 	GradPath_c GP(StartPoint, StreamDir_Forward, 200, GPType_Classic, GPTerminate_AtRhoValue, nullptr, nullptr, nullptr, &RhoCutoff, VolInfo, HessPtrs, GradPtrs, RhoPtr);
// // 
// // 	GP.Seed(true);
// // 
// // 	GP.SaveAsOrderedZone("test GP");
// // 	
// 
// 	/*
// 	 *	Testing the creation and initialzation of a fe_quad zone, because this part is crashing in the isosurface extraction tool.
// 	 *	This works, so wonder why it's crashing in the other function...
// 	 */
// 
// // 	FieldVecPointer_c XYZPtr;
// // 	vector<vec3> Nodes(4);
// // 	for (int i = 0; i < 2; ++i) for (int j = 0; j < 2; ++j) Nodes[i * 2 + j] << i << j << 0;
// // 	vector<int> Elems = { 1, 2, 4, 3 };
// // 
// // 	TecUtilDataSetAddZone("test zone", 4, 1, 4, ZoneType_FEQuad, nullptr);
// // 	int ZoneNum = TecUtilDataSetGetNumZones();
// // 
// // 	TecUtilDataLoadBegin();
// // 
// // 	XYZPtr.InitializeWritePtr(ZoneNum, { 1, 2, 3 });
// // 
// // 	for (int i = 0; i < 4; ++i){
// // 		TecUtilDataNodeSetByZone(ZoneNum, 1, i + 1, Elems[i]);
// // 		XYZPtr.Write(i, Nodes[i]);
// // 	}
// // 	
// // 	TecUtilDataLoadEnd();
// // 	
// 
// 	/*
// 	 *	Test rolling average-based time remaining on status
// 	 */
// // 
// // 	high_resolution_clock::time_point StartTime = high_resolution_clock::now();
// // 
// // 	StatusLaunch("test", AddOnID, TRUE, TRUE);
// // 
// // 	int N = 200;
// // 
// // 	for (int i = 0; i < N; ++i){
// // 		StatusUpdate(i, N, "test", AddOnID, StartTime);
// // 		Sleep(100 * (i % 5));
// // 	}
// // 
// // 	StatusDrop(AddOnID);
// // 	
// 
// 	/*
// 	 *	Test isoRho value based integration
// 	 */
// 
// 	/*
// 	 *	First on a typical GB, starting at the sphere end terminating not at a point.
// 	 */
// // 	vector<int> GPZoneNums = { 55, 191, 60, 126, 30, 124 };
// // 	vector<GradPath_c*> GPPtrs;
// // 	for (int i : GPZoneNums){
// // 		GPPtrs.push_back(new GradPath_c);
// // 		*GPPtrs.back() = GradPath_c(i, { 1, 2, 3, 4 }, AddOnID);
// // 	}
// // 	VolExtentIndexWeights_s VolInfo;
// // 	GetVolInfo(1, { 1, 2, 3 }, FALSE, VolInfo);
// // 	vector<FieldDataPointer_c> VarPtrs(1);
// // 	VarPtrs[0].InitializeReadPtr(1, 4);
// // 	vector<double> IntVals;
// // 	IntegrateUsingIsosurfaces(GPPtrs, 3, VolInfo, VarPtrs, IntVals);
// // 
// // 	/*
// // 	 *	Get the sphere so we can add the sphere-interior portion of the integral
// // 	 */
// // 	FESurface_c Sphere(22, 1, { 1, 2, 3 }, { 4 });
// // 	Sphere.DoIntegrationNew(2, TRUE);
// // 	vector<double> SphereTriangleAreas;
// // 	vector<vector<double> > SphereElemIntVals = Sphere.GetTriSphereIntValsByElem(&SphereTriangleAreas);
// // 
// // 	for (int i = 0; i < IntVals.size(); ++i){
// // 		IntVals[i] += SphereElemIntVals[62][i];
// // 
// // 		// 		TecUtilDialogMessageBox(to_string(IntVals[i]).c_str(), MessageBoxType_Information);
// // 	}
// // 
// // 	/*
// // 	 *	Now a bond-path-coincident GB, again not terminating at a point.
// // 	 */
// // 
// // 	GPZoneNums = { 73, 278, 74, 94, 48, 177, 47, 93 };
// // 	GPPtrs.clear();
// // 	for (int i : GPZoneNums){
// // 		GPPtrs.push_back(new GradPath_c);
// // 		*GPPtrs.back() = GradPath_c(i, { 1, 2, 3, 4 }, AddOnID);
// // 	}
// // 	IntegrateUsingIsosurfaces(GPPtrs, 2, VolInfo, VarPtrs, IntVals, true);
// // 
// // 	for (int i = 0; i < IntVals.size(); ++i){
// // 		IntVals[i] += SphereElemIntVals[20][i];
// // 
// // 		// 		TecUtilDialogMessageBox(to_string(IntVals[i]).c_str(), MessageBoxType_Information);
// // 	}
// // 	
// 
// 	/*
// 	 *	Test project to triangle function
// 	 */
// // vec3 OP, TP;
// // vector<vec3> tri;
// // int numIter = 3;
// // 
// // vector<vec3> axes = { {1,0,0},{0,1,0},{0,0,1} };
// // for (auto const & i : axes) {
// // 	SaveVec3VecAsScatterZone({ i * 2, {0,0,0} }, string("axis"), Green_C);
// // 	TecUtilZoneSetScatter(SV_SHOW, tecplot::toolbox::Set(TecUtilDataSetGetNumZones()).getRef(), 0, FALSE);
// // 	TecUtilZoneSetMesh(SV_SHOW, tecplot::toolbox::Set(TecUtilDataSetGetNumZones()).getRef(), 0, TRUE);
// // 	TecUtilZoneSetMesh(SV_LINETHICKNESS, tecplot::toolbox::Set(TecUtilDataSetGetNumZones()).getRef(), 0.2, 0);
// // 	TecUtilZoneSetMesh(SV_COLOR, tecplot::toolbox::Set(TecUtilDataSetGetNumZones()).getRef(), 0, Black_C);
// // 	TecUtilZoneSetActive(tecplot::toolbox::Set(TecUtilDataSetGetNumZones()).getRef(), AssignOp_MinusEquals);
// // }
// // 
// // tri = { { 1,0,0 }, { 0,1,0 }, { 0,0,1 } };
// // for (int i = 0; i < numIter; ++i) {
// // 	OP = { i * 1. / (double)numIter, 0, 0 };
// // 	for (int j = 0; j < numIter; ++j) {
// // 		OP[1] += 1. / (double)numIter;
// // 		for (int k = 0; k < numIter; ++k) {
// // 			OP[2] += 1. / (double)numIter;
// // 
// // 			SaveVec3VecAsScatterZone({ tri[0], tri[1], tri[2], tri[0] }, string("triangle"));
// // 			TecUtilZoneSetScatter(SV_SHOW, tecplot::toolbox::Set(TecUtilDataSetGetNumZones()).getRef(), 0, FALSE);
// // 			TecUtilZoneSetMesh(SV_SHOW, tecplot::toolbox::Set(TecUtilDataSetGetNumZones()).getRef(), 0, TRUE);
// // 			TecUtilZoneSetMesh(SV_LINETHICKNESS, tecplot::toolbox::Set(TecUtilDataSetGetNumZones()).getRef(), 0.2, 0);
// // 			TecUtilZoneSetMesh(SV_COLOR, tecplot::toolbox::Set(TecUtilDataSetGetNumZones()).getRef(), 0, Red_C);
// // 			TecUtilZoneSetActive(tecplot::toolbox::Set(TecUtilDataSetGetNumZones()).getRef(), AssignOp_MinusEquals);
// // 
// // 			double testProj = PointDistanceToTriangleSquared(OP, TP, tri[0], tri[1], tri[2]);
// // 			SaveVec3VecAsScatterZone({ OP }, string("old point"), Red_C);
// // 			TecUtilZoneSetScatter(SV_FRAMESIZE, tecplot::toolbox::Set(TecUtilDataSetGetNumZones()).getRef(), 1, 0);
// // 			TecUtilZoneSetActive(tecplot::toolbox::Set(TecUtilDataSetGetNumZones()).getRef(), AssignOp_MinusEquals);
// // 			SaveVec3VecAsScatterZone({ TP }, string("new point"), Green_C);
// // 			TecUtilZoneSetScatter(SV_FRAMESIZE, tecplot::toolbox::Set(TecUtilDataSetGetNumZones()).getRef(), 1, 0);
// // 			TecUtilZoneSetActive(tecplot::toolbox::Set(TecUtilDataSetGetNumZones()).getRef(), AssignOp_MinusEquals);
// // 
// // 			SaveVec3VecAsScatterZone({ OP, TP }, string("old to new point"), Green_C);
// // 			TecUtilZoneSetScatter(SV_SHOW, tecplot::toolbox::Set(TecUtilDataSetGetNumZones()).getRef(), 0, FALSE);
// // 			TecUtilZoneSetMesh(SV_SHOW, tecplot::toolbox::Set(TecUtilDataSetGetNumZones()).getRef(), 0, TRUE);
// // 			TecUtilZoneSetMesh(SV_LINETHICKNESS, tecplot::toolbox::Set(TecUtilDataSetGetNumZones()).getRef(), 0.2, 0);
// // 			TecUtilZoneSetMesh(SV_COLOR, tecplot::toolbox::Set(TecUtilDataSetGetNumZones()).getRef(), 0, Blue_C);
// // 			TecUtilZoneSetActive(tecplot::toolbox::Set(TecUtilDataSetGetNumZones()).getRef(), AssignOp_MinusEquals);
// // 		}
// // 	}
// // }
// // 
// // tri = { { 0,0,0 }, { 0,1,0 }, { 1,0,0 } };
// // for (int i = 0; i < numIter; ++i) {
// // 	OP = { i * 1. / (double)numIter, 0, 0 };
// // 	for (int j = 0; j < numIter; ++j) {
// // 		OP[1] += 1. / (double)numIter;
// // 		for (int k = 0; k < numIter; ++k) {
// // 			OP[2] += 1. / (double)numIter;
// // 
// // 			SaveVec3VecAsScatterZone({ tri[0], tri[1], tri[2], tri[0] }, string("triangle"));
// // 			TecUtilZoneSetScatter(SV_SHOW, tecplot::toolbox::Set(TecUtilDataSetGetNumZones()).getRef(), 0, FALSE);
// // 			TecUtilZoneSetMesh(SV_SHOW, tecplot::toolbox::Set(TecUtilDataSetGetNumZones()).getRef(), 0, TRUE);
// // 			TecUtilZoneSetMesh(SV_LINETHICKNESS, tecplot::toolbox::Set(TecUtilDataSetGetNumZones()).getRef(), 0.2, 0);
// // 			TecUtilZoneSetMesh(SV_COLOR, tecplot::toolbox::Set(TecUtilDataSetGetNumZones()).getRef(), 0, Red_C);
// // 			TecUtilZoneSetActive(tecplot::toolbox::Set(TecUtilDataSetGetNumZones()).getRef(), AssignOp_MinusEquals);
// // 
// // 			double testProj = PointDistanceToTriangleSquared(OP, TP, tri[0], tri[1], tri[2]);
// // 			SaveVec3VecAsScatterZone({ OP }, string("old point"), Red_C);
// // 			TecUtilZoneSetScatter(SV_FRAMESIZE, tecplot::toolbox::Set(TecUtilDataSetGetNumZones()).getRef(), 1, 0);
// // 			TecUtilZoneSetActive(tecplot::toolbox::Set(TecUtilDataSetGetNumZones()).getRef(), AssignOp_MinusEquals);
// // 			SaveVec3VecAsScatterZone({ TP }, string("new point"), Green_C);
// // 			TecUtilZoneSetScatter(SV_FRAMESIZE, tecplot::toolbox::Set(TecUtilDataSetGetNumZones()).getRef(), 1, 0);
// // 			TecUtilZoneSetActive(tecplot::toolbox::Set(TecUtilDataSetGetNumZones()).getRef(), AssignOp_MinusEquals);
// // 
// // 			SaveVec3VecAsScatterZone({ OP, TP }, string("old to new point"), Green_C);
// // 			TecUtilZoneSetScatter(SV_SHOW, tecplot::toolbox::Set(TecUtilDataSetGetNumZones()).getRef(), 0, FALSE);
// // 			TecUtilZoneSetMesh(SV_SHOW, tecplot::toolbox::Set(TecUtilDataSetGetNumZones()).getRef(), 0, TRUE);
// // 			TecUtilZoneSetMesh(SV_LINETHICKNESS, tecplot::toolbox::Set(TecUtilDataSetGetNumZones()).getRef(), 0.2, 0);
// // 			TecUtilZoneSetMesh(SV_COLOR, tecplot::toolbox::Set(TecUtilDataSetGetNumZones()).getRef(), 0, Blue_C);
// // 			TecUtilZoneSetActive(tecplot::toolbox::Set(TecUtilDataSetGetNumZones()).getRef(), AssignOp_MinusEquals);
// // 		}
// // 	}
// // }
// // 
// // return;
// 
// 	/*
// 	 *	Test surface grad path
// 	 */
// // 	int SurfZoneNum = 50;
// // 	int RhoVarNum = 4;
// // 	VolExtentIndexWeights_s VolInfo;
// // 	GetVolInfo(1, { 1,2,3 }, FALSE, VolInfo);
// // 	FieldDataPointer_c RhoPtr;
// // 	GetReadPtrsForZone(1, RhoVarNum, vector<int>(), vector<int>(), RhoPtr, vector<FieldDataPointer_c>(), vector<FieldDataPointer_c>());
// // 
// // 	FESurface_c Surf(SurfZoneNum, { 1,2,3 });
// // // 	Surf.GenerateElemConnectivity();
// // // 	Surf.GenerateElemMidpoints();
// // 	Surf.GeneratePointElementDistanceCheckData();
// // 	double TermVal = 0.001;
// // 
// // 	GradPath_c SurfGP({ -0.649,-0.672,-0.166 }, StreamDir_Reverse, 1000, GPType_Classic, GPTerminate_AtRhoValue, nullptr, nullptr, nullptr, &TermVal, VolInfo, vector<FieldDataPointer_c>(), vector<FieldDataPointer_c>(), RhoPtr, &Surf);
// // 
// // 	SurfGP.Seed(false);
// // 	if (SurfGP.IsMade()) SurfGP.SaveAsOrderedZone("SurfGradPath", Green_C);
// // 
// // 	TecUtilLockFinish(AddOnID);
// // 	return;
// // 	
// 
// 	/*
// 	 *	test split string with blank removal
// 	 */
// // vector<string> out = SplitString("381,357,375,,345,944,1329,", ",", false, true);
// // out = SplitString("381,357,375,,345,944,1329,", ",", false, false);
// // out = SplitString("381,357,375,,345,944,1329,", ",", true);
// // vector<int> outi = SplitStringInt("381,357,375,,345,944,1329,");
// // outi = SplitStringInt("381,357,375,,345,944,1329,");
// // outi = SplitStringInt("381,357,375,,345,944,1329,");
// // vector<double> outd = SplitStringDbl("381,357,375,,345,944,1329,");
// // outd = SplitStringDbl("381,357,375,,345,944,1329,");
// // outd = SplitStringDbl("381,357,375,,345,944,1329,");
// 
// 
// 
// // 	
// // 
// // 	/*
// // 	 *	test surface-sphere intersection function for file Z:\Tecplot\StorageWorkspace\2018_11_GBS-sphere-triangulationupdate\Al43.lpk
// // 	 */
// // 	int SphereZoneNum = 364;
// // 	int CPZoneNum = 3;
// // 
// // 	vec3 SphereCenter;
// // 	for (int i = 0; i < 3; ++i){
// // 		SphereCenter[i] = TecUtilDataValueGetByZoneVar(CPZoneNum, i + 1, 2);
// // 	}
// // 
// // 	double SphereRadius = stod(AuxDataZoneGetItem(SphereZoneNum, CSMAuxData.GBA.SphereRadius));
// // 
// // 	/*
// // 	*	Test gradpath_c sphere intersection code
// // 	*/
// // // 	vec3 IntPoint;
// // // 	bool DoesIntersect = GradPath_c(21, { 1,2,3,4 }, AddOnID).GetSphereIntersectionPoint(SphereCenter, SphereRadius, IntPoint);
// // // 	SaveVec3VecAsScatterZone({ IntPoint });
// // // 	DoesIntersect = GradPath_c(22, { 1,2,3,4 }, AddOnID).GetSphereIntersectionPoint(SphereCenter, SphereRadius, IntPoint);
// // // 	return;
// // 	/*
// // 	 *	end test gradpath_c sphere intersection code
// // 	 */
// // 
// // 	FESurface_c OldSphere(SphereZoneNum, { 1,2,3 });
// // 	vector<vec3> const * SphereXYZPtr = OldSphere.GetXYZListPtr();
// // 	vector<vector<int> > const * SphereElemListPtr = OldSphere.GetElemListPtr();
// // 
// // 
// // // 	auto tmpNodes = *SphereXYZPtr;
// // // 	auto tmpElems = *SphereElemListPtr;
// // // 	vector<std::set<int> > trisOfNode(tmpNodes.size());
// // // 	for (int ti = 0; ti < tmpElems.size(); ++ti) {
// // // 		for (auto ni : tmpElems[ti]) {
// // // 			trisOfNode[ni].insert(ti);
// // // 		}
// // // 	}
// // // 	vector<int> newNodeNums, newElemNums;
// // // 
// // // 	TriangleEdgeMidPointSubdivideAroundNodes(tmpNodes, tmpElems, trisOfNode, { 0 }, newNodeNums, newElemNums);
// // // 	FESurface_c tmpSphere(tmpNodes, tmpElems);
// // // 	tmpSphere.SaveAsTriFEZone({ 1,2,3 }, "refined sphere");
// // // 	TecUtilZoneSetActive(Set(tmpSphere.GetZoneNum()).getRef(), AssignOp_Equals);
// // // 
// // // 	TecUtilLockFinish(AddOnID);
// // // 	return;
// // 
// // 
// // 	/*
// // 	 *	Get average sphere triangulation edge length so that we can 
// // 	 *	resample intersection paths at roughly the same spacing.
// // 	 *	I'll do this by looping over triangles, which will double count every 
// // 	 *	edge so the average should still be valid.
// // 	 */
// // 	double SphereEdgeLenMean = 0.0;
// // 	int SphereNumEdges = 0;
// // 	for (auto const & t : *SphereElemListPtr){
// // 		for (int i = 0; i < 3; ++i){
// // 			SphereEdgeLenMean += Distance(SphereXYZPtr->at(t[i]), SphereXYZPtr->at(t[(i + 1) % 3]));
// // 			SphereNumEdges++;
// // 		}
// // 	}
// // 	SphereEdgeLenMean /= (double)SphereNumEdges;
// // 
// // 	vector<vector<vec3> > AllIntSurfPoints;
// // 	vector<int> IntSurfZoneNums;
// // 	vector<FESurface_c> IntSurfs;
// // 	vector<vector<GradPath_c> > IntSurfRingCagePaths;
// // 
// // 	CritPoints_c CPs(CPZoneNum, { 1,2,3 }, 5, 4);
// // 	VolExtentIndexWeights_s VolInfo;
// // 	GetVolInfo(1, { 1,2,3 }, FALSE, VolInfo);
// // 	FieldDataPointer_c RhoPtr;
// // 	vector<FieldDataPointer_c> GradPtrs, HessPtrs;
// // 	GetReadPtrsForZone(1, 4, { 6,7,8 }, { 9,10,11,12,13,14 }, RhoPtr, GradPtrs, HessPtrs);
// // 
// // 	MultiRootParams_s MR;
// // 	MR.VolInfo = &VolInfo;
// // 	MR.IsPeriodic = FALSE;
// // 	MR.HasGrad = TRUE;
// // 	MR.HasHess = TRUE;
// // 	MR.RhoPtr = &RhoPtr;
// // 	MR.GradPtrs = &GradPtrs;
// // 	MR.HessPtrs = &HessPtrs;
// // 	MR.BasisVectors = &VolInfo.BasisVectors;
// // 
// // 
// // 	vector<VolExtentIndexWeights_s> ThVolInfo(omp_get_num_procs(), VolInfo);
// // 
// // 
// // 	/*
// // 	 *	Get sphere intersecting ring surfaces and bond paths
// // 	 */
// // 
// // 	vector<vec3> IntersectingBondPathPoints;
// // 	vector<GradPath_c> IntersectingBondPaths;
// // 
// // // 	for (int z : {182, 249}){
// // 	for (int z = 1; z <= TecUtilDataSetGetNumZones(); ++z){
// // 		if (AuxDataZoneItemMatches(z, CSMAuxData.CC.ZoneSubType, CSMAuxData.CC.ZoneSubTypeRSSegment)){
// // 			FESurface_c Surf(z, { 1,2,3 });
// // 			vector<vec3> IntPoints = Surf.GetSphereIntersectionPath(SphereCenter, SphereRadius);
// // 			if (IntPoints.size() > 1) {
// // 				vector<vec3> ResampledIntPoints;
// // 				double PathLen = 0;
// // 				for (int i = 0; i < IntPoints.size() - 1; ++i)
// // 					PathLen += Distance(IntPoints[i], IntPoints[i + 1]);
// // 				if (Vec3PathResample(IntPoints, MAX(int(PathLen / SphereEdgeLenMean) + 1, 4), ResampledIntPoints)){
// // 					IntPoints = ResampledIntPoints;
// // // 					Surf.GenerateElemConnectivity();
// // // 					Surf.GenerateElemMidpoints();
// // 					Surf.GeneratePointElementDistanceCheckData();
// // 					IntSurfs.push_back(Surf);
// // 
// // 					/*
// // 					 *	Get the ring-cage paths originating at the ring point.
// // 					 *	Repeat the generation of the pair of ring path segments
// // 					 *	with increasing distance of the seed point in each direction
// // 					 *	until they terminate in different directions.
// // 					 *	
// // 					 */
// // 					int CPNum = -1;
// // 					for (int i = 0; i < 3; ++i) {
// // 						if (AuxDataZoneItemMatches(z, CSMAuxData.CC.ZFSCornerCPTypes[i], "Ring")) {
// // 							CPNum = stoi(AuxDataZoneGetItem(z, CSMAuxData.CC.ZFSCornerCPNumStrs[i])) - 1;
// // 							break;
// // 						}
// // 					}
// // 					if (CPNum >= 0) {
// // 						vec3 CPPos, PrincDir;
// // 						CPPos = CPs.GetXYZ(CPTypeNum_Ring, CPNum);
// // 						vec3 EigVals;
// // 						mat33 EigVecs;
// // 						CalcEigenSystemForPoint(CPPos, EigVals, EigVecs, MR);
// // 						PrincDir = EigVecs.col(0);
// // 						double TermRadius = 0.1;
// // 
// // 						int CPGlobalNum = CPs.GetTotOffsetFromTypeNumOffset(CPTypeNum_Ring, CPNum);
// // 						double CPRho = ValAtPointByPtr(CPPos, VolInfo, RhoPtr);
// // 
// // 						IntSurfRingCagePaths.push_back(vector<GradPath_c>(2));
// // 
// // 						double iter = 1;
// // 						double offset = 0.02;
// // 						double dPdt = 1;
// // 						while (dPdt > 0) {
// // 							double dir = -1.0;
// // 							for (auto & GP : IntSurfRingCagePaths.back()) {
// // 								vec3 SeedPt = CPPos + (PrincDir * offset * iter * dir);
// // 								GP.SetupGradPath(
// // 									SeedPt,
// // 									StreamDir_Reverse, 100, GPType_Classic, GPTerminate_AtCP, nullptr,
// // 									&CPs, &TermRadius, nullptr, VolInfo, HessPtrs, GradPtrs, RhoPtr);
// // 								GP.SetStartEndCPNum(CPGlobalNum, 0);
// // 								GP.Seed();
// // 								GP.PointPrepend(CPPos, CPRho);
// // 
// // 								GP.SaveAsOrderedZone("ring path for zone " + to_string(z) + " with offset of " + to_string(offset * iter * dir), Green_C);
// // 								TecUtilZoneSetActive(Set(TecUtilDataSetGetNumZones()).getRef(), AssignOp_MinusEquals);
// // 
// // 								dir += 2.0;
// // 							}
// // 							dPdt = dot(IntSurfRingCagePaths.back()[0].XYZAt(-1) - IntSurfRingCagePaths.back()[0].XYZAt(-2), IntSurfRingCagePaths.back()[1].XYZAt(-1) - IntSurfRingCagePaths.back()[1].XYZAt(-2));
// // 							iter++;
// // 						}
// // 					}
// // 					else {
// // 						TecUtilDialogErrMsg("Failed to find ring CP for sphere-intersecting SZFS zone");
// // 					}
// // 
// // 					/* 
// // 					 * Now project each point back to the sphere radius
// // 					 */
// // 					for (auto & p : IntPoints){
// // 						p = SphereCenter + normalise(p - SphereCenter) * SphereRadius;
// // 					}
// // 
// // // 					SaveVec3VecAsScatterZone(IntPoints, string("Intersection of zone " + to_string(z)));
// // 					AllIntSurfPoints.push_back(IntPoints);
// // 					IntSurfZoneNums.push_back(z);
// // 				}
// // 			}
// // 		}
// // 		else if (AuxDataZoneItemMatches(z, CSMAuxData.CC.ZoneSubType, CSMAuxData.CC.ZoneSubTypeBondPathSegment)){
// // 			vec3 IntPoint;
// // 			GradPath_c TmpGP(z, { 1,2,3,4 }, AddOnID);
// // 			if (TmpGP.GetSphereIntersectionPoint(SphereCenter, SphereRadius, IntPoint)){
// // 				IntersectingBondPathPoints.push_back(IntPoint);
// // 				IntersectingBondPaths.push_back(TmpGP);
// // 			}
// // 		}
// // 	}
// // 
// // // 	/*
// // // 	 *	12/26/2018
// // // 	 *	Craig's triangulation update code has a problem where some of the 
// // // 	 *	ring surface intersection points aren't being recorded properly
// // // 	 *	so that their corresponding gradient bundles can be created
// // // 	 *	correctly.
// // // 	 *	Rather than try to fix that code or manually identify intersection
// // // 	 *	nodes, I'll just workaround the need for the triangulation update
// // // 	 *	in the first place.
// // // 	 *	In order to know which sphere nodes are on ring surface-sphere 
// // // 	 *	intersections I just need to make sure that there's at least
// // // 	 *	one preexisting node on the sphere that can be moved to each
// // // 	 *	of the intersection points, so I just need to move the closest 
// // // 	 *	sphere node to each intersection node. If a single sphere node
// // // 	 *	is closest to two or more intersection nodes, subdivide the
// // // 	 *	neighborhood of the sphere node until each intersection node
// // // 	 *	gets its own sphere node.
// // // 	 */
// // // 
// // // 	/*
// // // 	 *	First, perform this process for the bond path intersections,
// // // 	 *  because a sphere point that moves to a bond path is allowed 
// // // 	 *  to be associated with more than one ring surface intersection
// // // 	 *  point but not more than one bond path intersection point.
// // // 	 */
// // // 
// // //  	vector<vec3> NewNodes;
// // //  	vector<vector<int> > NewElems;
// // // 	std::map<int,int> outRingConstraintNodeIndices, outRingConstraintSegmentIndices, outBondConstraintNodeIndices;
// // //  	vector<vector<vec3> > InIntersectingBondPathPoints;
// // //  	for (auto p : IntersectingBondPathPoints){
// // //  		InIntersectingBondPathPoints.push_back({ p });
// // //  	}
// // //  	
// // //  	UpdateSubdivideSphericalTriangulationWithConstraintNodesAndSegments(*SphereXYZPtr, *SphereElemListPtr, SphereCenter, SphereRadius, InIntersectingBondPathPoints, NewNodes, NewElems, outBondConstraintNodeIndices, outRingConstraintSegmentIndices);
// // //  
// // // 	vector<vec3> outConstrainedNodes;
// // //  	for (auto const & i : outBondConstraintNodeIndices) {
// // // 		outConstrainedNodes.push_back(NewNodes[i.first]);
// // //  	}
// // //  	SaveVec3VecAsScatterZone(outConstrainedNodes, "bond nodes", Red_C);
// // // 
// // // 	TecUtilDialogMessageBox("before updating triangulation for ring surfaces", MessageBoxType_Information);
// // // 
// // // 	UpdateSubdivideSphericalTriangulationWithConstraintNodesAndSegments(NewNodes, NewElems, SphereCenter, SphereRadius, AllIntSurfPoints, NewNodes, NewElems, outRingConstraintNodeIndices, outRingConstraintSegmentIndices);
// // // 
// // // 	outConstrainedNodes.clear();
// // // 	for (auto const & i : outRingConstraintNodeIndices) {
// // // 		outConstrainedNodes.push_back(NewNodes[i.first]);
// // // 	}
// // // 	SaveVec3VecAsScatterZone(outConstrainedNodes, "ring nodes", Green_C);
// // //  
// // //  	FESurface_c NewSphere(NewNodes, NewElems);
// // //  
// // //  	if (NewSphere.IsMade()) {
// // //  		NewSphere.SaveAsTriFEZone({ 1,2,3 }, "UpdatedSphere");
// // //  		TecUtilZoneSetActive(Set(TecUtilDataSetGetNumZones()).getRef(), AssignOp_Equals);
// // //  	}
// // //  	else {
// // //  		TecUtilDialogErrMsg("Failed to create updated triangulation");
// // //  	}
// // // 
// // // 	TecUtilLockFinish(AddOnID);
// // // 	return;
// // // 
// //  	/*
// //  	 *	Check for nearly equal intersection path segment endpoints,
// //  	 *	i.e. where two paths intersect.
// //  	 *	These intersections should only occur at a bond-path-sphere intersection point,
// //  	 *	at which there is already a constrained sphere node, so when path intersections are
// //  	 *	found, move all the coincident points to the closest node on the sphere.
// //  	 */
// //  	double CoincidentCheckEpsilon = SphereEdgeLenMean * 0.2;
// //  	CoincidentCheckEpsilon *= CoincidentCheckEpsilon;
// //  	for (int i = 0; i < AllIntSurfPoints.size() - 1; ++i) {
// //  		vector<std::pair<int, int> > CoincidentPointIndices;
// //  		for (int j = i+1; j < AllIntSurfPoints.size(); ++j){
// //  			if (DistSqr(AllIntSurfPoints[i][0], AllIntSurfPoints[j][0]) <= CoincidentCheckEpsilon){
// //  				CoincidentPointIndices.push_back(std::pair<int, int>(i, 0));
// //  				CoincidentPointIndices.push_back(std::pair<int, int>(j, 0));
// //  			}
// //  			else if (DistSqr(AllIntSurfPoints[i][0], AllIntSurfPoints[j].back()) <= CoincidentCheckEpsilon) {
// //  				CoincidentPointIndices.push_back(std::pair<int, int>(i, 0));
// //  				CoincidentPointIndices.push_back(std::pair<int, int>(j, AllIntSurfPoints[j].size() - 1));
// //  			}
// //  			else if (DistSqr(AllIntSurfPoints[i].back(), AllIntSurfPoints[j][0]) <= CoincidentCheckEpsilon) {
// //  				CoincidentPointIndices.push_back(std::pair<int, int>(i, AllIntSurfPoints[i].size() - 1));
// //  				CoincidentPointIndices.push_back(std::pair<int, int>(j, 0));
// //  			}
// //  			else if (DistSqr(AllIntSurfPoints[i].back(), AllIntSurfPoints[j].back()) <= CoincidentCheckEpsilon) {
// //  				CoincidentPointIndices.push_back(std::pair<int, int>(i, AllIntSurfPoints[i].size() - 1));
// //  				CoincidentPointIndices.push_back(std::pair<int, int>(j, AllIntSurfPoints[j].size() - 1));
// //  			}
// //  		}
// //  
// //  		if (CoincidentPointIndices.size() > 0){
// //  			int SphereNodeNum = -1;
// //  			for (auto const & p : CoincidentPointIndices){
// //  				double MinNodeDistSqr = DBL_MAX;
// //  				double TmpNodeDistSqr;
// //  				int MinNodeIndex = -1;
// //  				for (int ni = 0; ni < SphereXYZPtr->size(); ++ni){
// //  					TmpNodeDistSqr = DistSqr(AllIntSurfPoints[p.first][p.second], SphereXYZPtr->at(ni));
// //  					if (TmpNodeDistSqr < MinNodeDistSqr){
// //  						MinNodeDistSqr = TmpNodeDistSqr;
// //  						MinNodeIndex = ni;
// //  					}
// //  				}
// //  				if (MinNodeIndex >= 0){
// //  // 					if (SphereNodeNum >= 0){
// //  // 						if (SphereNodeNum != MinNodeIndex) {
// //  // 							TecUtilDialogErrMsg("Coincident intersection path endpoints do not share same closest sphere node!");
// //  // 						}
// //  // 					}
// //  // 					else{
// //  // 						SphereNodeNum = MinNodeIndex;
// //  // 					}
// //  					AllIntSurfPoints[p.first][p.second] = SphereXYZPtr->at(MinNodeIndex);
// //  				}
// //  			}
// //  		}
// //  
// //  		SaveVec3VecAsScatterZone(AllIntSurfPoints[i], string("Intersection of zone " + to_string(IntSurfZoneNums[i])), White_C);
// //  // 		SetZoneStyle({}, ZoneStyle_Path, Blue_C, 0.2);
// //  	}
// //  	SaveVec3VecAsScatterZone(AllIntSurfPoints.back(), string("Intersection of zone " + to_string(IntSurfZoneNums.back())), White_C);
// //  // 	SetZoneStyle({}, ZoneStyle_Path, Blue_C, 0.2);
// //  
// //  // 	return;
// //  
// //  	vector<tpcsm::Vec3> initialNodeXYZs(SphereXYZPtr->size());
// //  	for (int i = 0; i < initialNodeXYZs.size(); ++i)
// //  		initialNodeXYZs[i] = tpcsm::Vec3(SphereXYZPtr->at(i)[0], SphereXYZPtr->at(i)[1], SphereXYZPtr->at(i)[2]);
// //  
// //  	vector<TriNodes> initialTriangleNodes(SphereElemListPtr->size());
// //  	for (int i = 0; i < initialTriangleNodes.size(); ++i)
// //  		initialTriangleNodes[i] = TriNodes(SphereElemListPtr->at(i)[0], SphereElemListPtr->at(i)[1], SphereElemListPtr->at(i)[2]);
// //  
// //  	tpcsm::Vec3 SphereCenterVec(SphereCenter[0], SphereCenter[1], SphereCenter[2]);
// //  
// //  	vector<tpcsm::Vec3> constraintNodes;
// //  	vector<Edge> constraintSegments;
// //  	vector<int> constraintSegmentSurfNums;
// //  	for (int j = 0; j < AllIntSurfPoints.size(); ++j) {
// //  		vector<vec3> IntPoints = AllIntSurfPoints[j];
// //  		constraintNodes.push_back(tpcsm::Vec3(IntPoints[0][0], IntPoints[0][1], IntPoints[0][2]));
// //  		for (int i = 1; i < IntPoints.size(); ++i) {
// //  			constraintNodes.push_back(tpcsm::Vec3(IntPoints[i][0], IntPoints[i][1], IntPoints[i][2]));
// //  			constraintSegments.push_back(Edge(constraintNodes.size() - 2, constraintNodes.size() - 1));
// //  			constraintSegmentSurfNums.push_back(j);
// //  		}
// //  	}
// //  
// //  	/*
// //  	 *	Need to remove duplicate constraint nodes and update edge connectivity to reflect the change.
// //  	 */
// //  	vector<int> nodeIndices(constraintNodes.size());
// //  	for (int i = 0; i < nodeIndices.size(); ++i)
// //  		nodeIndices[i] = i;
// //  	vector<tpcsm::Vec3> tmpNodes;
// //  	tmpNodes.reserve(constraintNodes.size());
// //  	for (int ni = 0; ni < constraintNodes.size() - 1; ++ni){
// //  		bool isFound = false;
// //  		for (int nj = ni + 1; nj < constraintNodes.size(); ++nj){
// //  			if (constraintNodes[nodeIndices[ni]] == constraintNodes[nodeIndices[nj]]){
// //  				/*
// //  				 *	duplicate node found, so search through edges, changing nj to ni when found.
// //  				 */
// //  				for (auto & e : constraintSegments){
// //  					if (e.first == nj)
// //  						e.first = ni;
// //  
// //  					if (e.second == nj)
// //  						e.second = ni;
// //  				}
// //  				nodeIndices[nj] = nodeIndices[ni];
// //  				isFound = true;
// //  			}
// //  		}
// //  	}
// //  	/* 
// //  	 * Now need to make the new list of nodes and update the edges again to reflect this.
// //  	 */
// //  	for (int ni = 0; ni < constraintNodes.size(); ++ni){
// //  		if (nodeIndices[ni] == ni){
// //  			/*
// //  			 *	Node is not duplicate
// //  			 */
// //  			tmpNodes.push_back(constraintNodes[ni]);
// //  			for (auto & e : constraintSegments){
// //  				if (e.first == ni)
// //  					e.first = tmpNodes.size() - 1;
// //  
// //  				if (e.second == ni)
// //  					e.second = tmpNodes.size() - 1;
// //  			}
// //  		}
// //  	}
// //  	constraintNodes = tmpNodes;
// //  
// //  	/*
// //  	 *	Now move nodes on the sphere that are sufficiently close to a constraint node
// //  	 *	to eliminate the ugly triangles that result when two points are very close.
// //  	 */
// //  // 	for (auto & n : initialNodeXYZs){
// //  // 		for (auto & p : constraintNodes){
// //  // 			if ((p - n).getNormSquared() <= CoincidentCheckEpsilon){
// //  // 				n = p;
// //  // 				break;
// //  // 			}
// //  // 		}
// //  // 	}
// //  
// //  	vector<bool> NodeMoved(initialNodeXYZs.size(), false);
// //  	for (auto & p : constraintNodes) {
// //  		double MinDistSqr = DBL_MAX, TmpDistSqr;
// //  		int MinI;
// //  		for (int ni = 0; ni < initialNodeXYZs.size(); ++ni) {
// //  			auto n = initialNodeXYZs[ni];
// //  			TmpDistSqr = (p - n).getNormSquared();
// //  			if (TmpDistSqr < MinDistSqr) {
// //  				MinDistSqr = TmpDistSqr;
// //  				MinI = ni;
// //  			}
// //  		}
// // 		if (!NodeMoved[MinI]) {
// // 			NodeMoved[MinI] = true;
// // 			initialNodeXYZs[MinI] = p;
// // 		}
// //  	}
// // 
// //  
// //  	/*
// //  	 *	Update triangulation to include ring surface intersections as edges
// //  	 */
// //  
// //  	vector<tpcsm::Vec3> OutXYZs;
// //  	vector<Index_t> OutConstrinedSegmentIndices;
// //  	vector<TriNodes> OutTris;
// //  	char const * StatusMessage;
// //  
// // 	if (tpcsm::updateSphericalTriangulation(initialNodeXYZs,
// // 		initialTriangleNodes,
// // 		SphereCenterVec,
// // 		SphereRadius,
// // 		constraintNodes,
// // 		constraintSegments,
// // 		false,
// // 		OutXYZs,
// // 		OutConstrinedSegmentIndices,
// // 		OutTris,
// // 		StatusMessage))
// // 	{
// // 
// // 		vector<vec3> SphereNodes(OutXYZs.size());
// // 		for (int i = 0; i < SphereNodes.size(); ++i)
// // 			SphereNodes[i] << OutXYZs[i].x() << OutXYZs[i].y() << OutXYZs[i].z();
// // 
// // 		vector<vector<int> > SphereElems(OutTris.size());
// // 		for (int i = 0; i < OutTris.size(); ++i)
// // 			SphereElems[i] = { OutTris[i].v1(), OutTris[i].v2(), OutTris[i].v3() };
// // 
// // 		FESurface_c Sphere(SphereNodes, SphereElems);
// // 
// // 		if (Sphere.IsMade()) {
// // 			Sphere.SaveAsTriFEZone({ 1,2,3 }, "UpdatedSphere");
// // 			TecUtilZoneSetActive(Set(TecUtilDataSetGetNumZones()).getRef(), AssignOp_MinusEquals);
// // 		}
// // 		else {
// // 			TecUtilDialogErrMsg("Failed to create updated triangulation");
// // 		}
// // 
// // 		/*
// // 		 *	seed surface gradient paths for all the nodes that
// // 		 *	fall on constraint edges.
// // 		 */
// // 
// // 		 // 
// // 		 // 		for (int i = 0; i < OutConstrinedSegmentIndices.size(); ++i){
// // 		 // 			if (OutConstrinedSegmentIndices[i] >= 0) {
// // 		 // 				int ConstrainedSegNum = i,
// // 		 // 					ConstrainedSegSurfNum = constraintSegmentSurfNums[OutConstrinedSegmentIndices[i]];
// // 		 // 				if (ConstrainedSegSurfNum < IntSurfs.size()){
// // 		 // 					FESurface_c * SurfPtr = &IntSurfs[ConstrainedSegSurfNum];
// // 		 // 					GradPath_c GP(NewNodes[i], StreamDir_Both, 500, GPType_Classic, GPTerminate_AtBoundary, nullptr, nullptr, nullptr, nullptr, VolInfo, vector<FieldDataPointer_c>(), vector<FieldDataPointer_c>(), RhoPtr, SurfPtr);
// // 		 // 					GP.Seed();
// // 		 // 					if (GP.IsMade()){
// // 		 // 						GP.SaveAsOrderedZone(string("SurfGP node " + to_string(i) + " cSegNum " + to_string(i) + " surfNum " + to_string(ConstrainedSegSurfNum) + " surfZoneNum " + to_string(SurfPtr->GetZoneNum())), Green_C);
// // 		 // 					}
// // 		 // 				}
// // 		 // 			}
// // 		 // 		}
// // 		 // 		
// // 
// // 				/*
// // 				 *	Now we need to seed the gradient paths from the nodes and along the
// // 				 *	edges of the triangular mesh.
// // 				 *	***actually, don't put gradient paths along the edges. It introduces
// // 				 *		unnecessary complexity, makes the algorithm much more expensive,
// // 				 *		and only captures a little bit more information about the curvature
// // 				 *		of the zero flux surfaces formed from stitching the gradient
// // 				 *		paths together.
// // 				 *		This also simplifies the tetrahedral decomposition used in the
// // 				 *		integration code.
// // 				 *		The information lost by excluding edges can be recovered while only using
// // 				 *		the nodes of triangles by simply refining the triangulation.
// // 				 *
// // 				 *	There are three types of nodes, according to the treatment of the paths
// // 				 *	they form:
// // 				 *		C type:		Paths that terminate at a cage CP or at a rho cutoff
// // 				 *					require no special treatment.
// // 				 *		R type:		Paths that coincide with a ring surface become will be seeded
// // 				 *					as surface paths and terminate at a ring point, then
// // 				 *					then split into the two halves of the ring path away from the
// // 				 *					ring point to the cage.
// // 				 *		B type:		Paths that coincide with a bond path will be replaced with the
// // 				 *					bond path, then a new gradient path will be seeded away from
// // 				 *					the bond point in the bond plane (interatomic surface) for each
// // 				 *					of the nodes that share an edge with the node of the original
// // 				 *					path.
// // 				 *
// // 				 *	We'll store all the gradient paths in a single vector and keep a new list of triangles
// // 				 *	where each triangle is defined by the indices of the gradient paths that correspond to
// // 				 *	its nodes.
// // 				 *	R, and B type nodes will result in new gradient paths being added to the list.
// // 				 */
// // 
// // 				 /*
// // 				  *	Find nodes on new sphere that correspond to the bond path intersections
// // 				  *	with the original sphere.
// // 				  */
// // 		std::unordered_map<int, GradPath_c> IntBondPathSegments; // stored according to intersecting node index
// // 		for (int bi = 0; bi < IntersectingBondPathPoints.size(); ++bi) {
// // 			IntBondPathSegments[Sphere.GetClosestNodeToPoint(IntersectingBondPathPoints[bi])] = IntersectingBondPaths[bi];
// // 		}
// // 
// // 		int NumNodes = SphereNodes.size(),
// // 			NumElems = SphereElems.size();
// // 
// // 		// Store type of each node
// // 		vector<GBATriangulatedSphereNodeElemType_e> NodeTypes(NumNodes, NETypeInvalid), ElemTypes(NumElems, NETypeInvalid);
// // 
// // 
// // 		// First get all node types;
// // //  #pragma omp parallel for
// // 		for (int ni = 0; ni < NumNodes; ++ni) {
// // 			if (IntBondPathSegments.count(ni) > 0)
// // 				NodeTypes[ni] = NETypeB;
// // 			else if (OutConstrinedSegmentIndices[ni] >= 0)
// // 				NodeTypes[ni] = NETypeR;
// // 			else {
// // 				NodeTypes[ni] = NETypeC;
// // 				for (int i = 0; i < constraintNodes.size(); ++i) {
// // 
// // 					if ((constraintNodes[i] - OutXYZs[ni]).getNormSquared() < CoincidentCheckEpsilon) {
// // 						for (int j = 0; j < constraintSegments.size(); ++j) {
// // 							if (constraintSegments[j].first == i || constraintSegments[j].second == i) {
// // 								OutConstrinedSegmentIndices[ni] = j;
// // 								NodeTypes[ni] = NETypeR;
// // 								break;
// // 							}
// // 						}
// // 					}
// // 					if (NodeTypes[ni] == NETypeR)
// // 						break;
// // 				}
// // 			}
// // 		}
// // 
// // 		// Get all element types
// // #pragma omp parallel for
// // 		for (int ti = 0; ti < NumElems; ++ti) {
// // 			for (int ni : SphereElems[ti]) {
// // 				ElemTypes[ti] = MAX(ElemTypes[ti], NodeTypes[ni]);
// // 			}
// // 		}
// // 
// // 		// #if defined _DEBUG
// // 		vector<int> BNodes, RNodes, CNodes;
// // 		for (int i = 0; i < NumNodes; ++i) {
// // 			if (NodeTypes[i] == NETypeB) BNodes.push_back(i);
// // 			else if (NodeTypes[i] == NETypeR) RNodes.push_back(i);
// // 			else if (NodeTypes[i] == NETypeC) CNodes.push_back(i);
// // 		}
// // 		int TotalNodes = BNodes.size() + RNodes.size() + CNodes.size();
// // 
// // 		TecUtilZoneSetActive(Set(Sphere.GetZoneNum()).getRef(), AssignOp_Equals);
// // 
// // 		vector<vector<int> > v = { BNodes, RNodes, CNodes };
// // 		vector<string> s = { "B nodes", "R nodes", "C nodes" };
// // 		vector<ColorIndex_t> c = { Red_C, Green_C, Red_C };
// // 		for (int i = 0; i < 3; ++i) {
// // 			auto vi = v[i];
// // 			vector<vec3> Vec3Vec;
// // 			for (auto i : vi) Vec3Vec.push_back(SphereNodes[i]);
// // 			SaveVec3VecAsScatterZone(Vec3Vec, s[i], c[i]);
// // 		}
// // 
// // 		// #endif
// // 	   //  		TecUtilLockFinish(AddOnID);
// // 	   //  		return;
// // 
// // 	   //  #define SUPERDEBUG
// // #ifdef SUPERDEBUG
// // 		vector<int> ElemTodo = { 6675 };
// // 		std::set<int> NodesTodo;
// // 		for (int ti : ElemTodo) {
// // 			for (int ni : SphereElems[ti])
// // 				NodesTodo.insert(ni);
// // 		}
// // 		vector<int> NodesTodoVec(NodesTodo.begin(), NodesTodo.end());
// // #endif
// // 
// // 		/*
// // 		 *	Setup and seed all type C and R paths.
// // 		 *	R paths will be combined with the ring path segments found above.
// // 		 *	We need R and C paths in order to get B and RB paths.
// // 		 */
// // 		vector<GradPath_c> GradPaths(NumNodes);
// // 		double TermRadius = 0.05;
// // 		double TermVal = 0.001;
// // 
// // #ifdef SUPERDEBUG
// // 		for (int nii = 0; nii < NodesTodoVec.size(); ++nii) {
// // 			int ni = NodesTodoVec[nii];
// // #else
// // #pragma omp parallel for schedule(dynamic)
// // 		for (int ni = 0; ni < NumNodes; ++ni) {
// // #endif
// // 			int ThNum = omp_get_thread_num();
// // 			if (NodeTypes[ni] == NETypeC) {
// // 				GradPaths[ni].SetupGradPath(SphereNodes[ni], StreamDir_Both, 100, GPType_Classic, GPTerminate_AtCP, nullptr, &CPs, &TermRadius, &TermVal, ThVolInfo[ThNum], HessPtrs, GradPtrs, RhoPtr);
// // 				GradPaths[ni].Seed();
// // 				// GradPaths[ni].SaveAsOrderedZone(string("C path for node " + to_string(ni)), Blue_C);
// // 			}
// // 			else if (NodeTypes[ni] == NETypeR) {
// // 				int SurfNum = constraintSegmentSurfNums[OutConstrinedSegmentIndices[ni]];
// // 				GradPaths[ni].SetupGradPath(SphereNodes[ni], StreamDir_Reverse, 100, GPType_Classic, GPTerminate_AtCP, nullptr, &CPs, &TermRadius, &TermVal, ThVolInfo[ThNum], HessPtrs, GradPtrs, RhoPtr, &IntSurfs[SurfNum]);
// // 				GradPaths[ni].Seed();
// // 				GradPath_c GP(SphereNodes[ni], StreamDir_Forward, 100, GPType_Classic, GPTerminate_AtCP, nullptr, &CPs, &TermRadius, &TermVal, ThVolInfo[ThNum], HessPtrs, GradPtrs, RhoPtr);
// // 				GP.Seed();
// // 				GradPaths[ni] += GP.SubGP(1, -1);
// // 
// // 				// 				GradPaths[ni].SaveAsOrderedZone(string("R path for node " + to_string(ni)), Cyan_C);
// // 			}
// // 			else continue;
// // 
// // 			if (GradPaths[ni].IsMade() && GradPaths[ni].RhoAt(0) < GradPaths[ni].RhoAt(-1)) GradPaths[ni].Reverse();
// // 
// // 			if (GradPaths[ni].RhoAt(0) <= GradPaths[ni].RhoAt(1))
// // 				GradPaths[ni].PointPrepend(GradPaths[ni].XYZAt(0), GradPaths[ni].RhoAt(0) + 1);
// // 		}
// // 		// 		for (int ni = 0; ni < GradPaths.size(); ++ni){
// // 		// 			if (GradPaths[ni].IsMade()){
// // 		// 				if (NodeTypes[ni] <= NodeTypeR){
// // 		// 					GradPaths[ni].SaveAsOrderedZone(string(NodeTypes[ni] == NodeTypeC ? "C" : "R") + string(" path for node ") + to_string(ni), NodeTypes[ni] == NodeTypeC ? Blue_C : Green_C);
// // 		// 				}
// // 		// 			}
// // 		// 		}
// // 
// // 			   /*
// // 				*	Now we can get the B paths.
// // 				*	We'll do this by looping over triangles because the selection of B  paths
// // 				*	is determined by the paths of the other two nodes on the triangle.
// // 				*	Start with serial pass to growNewElemGradPathInd, then another in parallel to
// // 				*	seed all the paths.
// // 				*/
// // 		std::unordered_map<int, vec3> BNodeSeedPts;
// // 		std::unordered_map<int, int> BNodeBondCPNums;
// // 		std::unordered_map<int, GradPath_c const*> BNodeBondPathPtrs;// , BNodeRingPathPtrs;
// //  // 		std::unordered_map<int, FESurface_c const *> BNodeRingSurfPtrs;
// // 
// // 		vector<vector<int> > NewElemGradPathInd(SphereElems);
// // #ifdef SUPERDEBUG
// // 		for (int ti : ElemTodo) {
// // #else
// // 		for (int ti = 0; ti < SphereElems.size(); ++ti) {
// // #endif
// // 			/*
// // 			 *	For each corner, check if it's type R or B.
// // 			 *	If type R then just need to update the triangle
// // 			 *	to point to the correct R path.
// // 			 *	If B, need to construct GPs to finish the paths,
// // 			 *	using information from the neighboring GPs.
// // 			 */
// // 			auto * tri = &NewElemGradPathInd[ti];
// // 			for (int c1 = 0; c1 < 3; ++c1) {
// // 				/*
// // 				 *	First pass over triangle to check for R type nodes
// // 				 *	to update them to the correct side.
// // 				 *	There could be 0, 1, or 2 R type nodes on the triangle.
// // 				 *	If 1, can use the guaranteed C type node to find the correct
// // 				 *	ring path.
// // 				 *	If 2, can check both ring paths for each node to find the pair
// // 				 *  whose end points are closest together.
// // 				 *  Also if 2 then we'll do them both together, so can break after
// // 				 *  they're found.
// // 				 */
// // 
// // 
// // 				if (NodeTypes[SphereElems[ti].at(c1)] == NETypeR && SphereElems[ti][c1] == tri->at(c1)) {
// // 
// // 					for (int ei = 1; ei < 3; ++ei) {
// // 						int c2 = (c1 + ei) % 3;
// // 
// // 						string logStr = "elem " + to_string(ti) + ", node " + to_string(SphereElems[ti].at(c1)) + "-" + to_string(SphereElems[ti].at(c2)) + ": ";
// // 						/*
// // 						 *	Check to see if this corner is R type
// // 						 */
// // 						if (ElemTypes[ti] == NETypeB && NodeTypes[SphereElems[ti].at(c2)] == NETypeR) {
// // 							/*
// // 							 * Both corners c1 and c2 are R type.
// // 							 * They each correspond to a pair of ring path segments.
// // 							 * Because these are two corners of a triangle, two of the four
// // 							 * ring path segments must terminate at (or near) the same cage CP,
// // 							 * while the other two go off to different cage CPs.
// // 							 * We'll find the correct ring path segments by finding the pair with the nearest
// // 							 * end points.
// // 							 */
// // 							int minI, minJ, SurfNumI, SurfNumJ;
// // 							SurfNumI = constraintSegmentSurfNums[OutConstrinedSegmentIndices[SphereElems[ti].at(c1)]];
// // 							SurfNumJ = constraintSegmentSurfNums[OutConstrinedSegmentIndices[SphereElems[ti].at(c2)]];
// // 							double MinDistSqr = DBL_MAX;
// // 							for (int ii = 0; ii < 2; ++ii) {
// // 								for (int jj = 0; jj < 2; ++jj) {
// // 									double TmpDistSqr = DistSqr(IntSurfRingCagePaths[SurfNumI][ii].XYZAt(-1),
// // 										IntSurfRingCagePaths[SurfNumJ][jj].XYZAt(-1));
// // 									if (TmpDistSqr < MinDistSqr) {
// // 										MinDistSqr = TmpDistSqr;
// // 										minI = ii;
// // 										minJ = jj;
// // 										// 										IntSurfRingCagePaths[SurfNumI][minI].SaveAsOrderedZone(logStr + "R neighbor surfI min path");
// // 										// 										TecUtilZoneSetActive(Set(TecUtilDataSetGetNumZones()).getRef(), AssignOp_MinusEquals);
// // 										// 										IntSurfRingCagePaths[SurfNumJ][minJ].SaveAsOrderedZone(logStr + "R neighbor surfJ min path");
// // 										// 										TecUtilZoneSetActive(Set(TecUtilDataSetGetNumZones()).getRef(), AssignOp_MinusEquals);
// // 									}
// // 									else {
// // 										// 										IntSurfRingCagePaths[SurfNumI][ii].SaveAsOrderedZone(logStr + "R neighbor surfI not min path");
// // 										// 										TecUtilZoneSetActive(Set(TecUtilDataSetGetNumZones()).getRef(), AssignOp_MinusEquals);
// // 										// 										IntSurfRingCagePaths[SurfNumJ][jj].SaveAsOrderedZone(logStr + "R neighbor surfJ not min path");
// // 										// 										TecUtilZoneSetActive(Set(TecUtilDataSetGetNumZones()).getRef(), AssignOp_MinusEquals);
// // 									}
// // 								}
// // 							}
// // 							/*
// // 							 *	Now the correct R paths for this triangle are those specified
// // 							 *	by minI and minJ
// // 							 */
// // 
// // 
// // 							 // 							IntSurfRingCagePaths[SurfNumI][minI].SaveAsOrderedZone(logStr + "R neighbor c1 ring path");
// // 							 // 							TecUtilZoneSetActive(Set(TecUtilDataSetGetNumZones()).getRef(), AssignOp_MinusEquals);
// // 							 // 							IntSurfRingCagePaths[SurfNumJ][minJ].SaveAsOrderedZone(logStr + "R neighbor c2 ring path");
// // 							 // 							TecUtilZoneSetActive(Set(TecUtilDataSetGetNumZones()).getRef(), AssignOp_MinusEquals);
// // 
// // 
// // 							GradPaths.push_back(ConcatenateResample({ GradPaths[SphereElems[ti].at(c1)], IntSurfRingCagePaths[SurfNumI][minI] }, 100));
// // 							tri->at(c1) = GradPaths.size() - 1;
// // 							//  							GradPaths.back().SaveAsOrderedZone(logStr + "R neighbor c1 new path");
// // 							//  							TecUtilZoneSetActive(Set(TecUtilDataSetGetNumZones()).getRef(), AssignOp_MinusEquals);
// // 
// // 							GradPaths.push_back(ConcatenateResample({ GradPaths[SphereElems[ti].at(c2)], IntSurfRingCagePaths[SurfNumJ][minJ] }, 100));
// // 							tri->at(c2) = GradPaths.size() - 1;
// // 							//  							GradPaths.back().SaveAsOrderedZone(logStr + "R neighbor c2 new path");
// // 							//  							TecUtilZoneSetActive(Set(TecUtilDataSetGetNumZones()).getRef(), AssignOp_MinusEquals);
// // 
// // 							break;
// // 						}
// // 						else if (NodeTypes[SphereElems[ti].at(c2)] == NETypeR && NodeTypes[SphereElems[ti].at((c2 + 1) % 3)] == NETypeR) {
// // 							/*
// // 							 *	All three nodes of the triangle are type R, meaning the
// // 							 *	gradient bundle has edges going to two different cages
// // 							 *	(It's impossible for a triangle to go to three different
// // 							 *	rings, so it has to be two).
// // 							 *
// // 							 *	In this case we'll just look for the cage CP (last point
// // 							 *	of ring path segment) that's closest to the ring cps (midpoint of the
// // 							 *	end points of the ring surface paths) for all three nodes.
// // 							 *	This seems like an OK test since the cage point needs to be
// // 							 *	present (otherwise couldn't have two rings this close to each
// // 							 *	other) and is the common terminus of the correct paths, where
// // 							 *	the other two cage points are not common among the correct paths.
// // 							 */
// // 
// // 							 // Get midpoint of ring points;
// // 							vec3 MidPt;
// // 							for (int ni : SphereElems[ti]) {
// // 								MidPt += GradPaths[ni].XYZAt(-1);
// // 							}
// // 							MidPt /= SphereElems[ti].size();
// // 
// // 							// Get which ring path segments terminate closer to the ring cp midpoint
// // 							// and add that path.
// // 							for (int ni = 0; ni < SphereElems[ti].size(); ++ni) {
// // 								double MinDistSqr = DBL_MAX, TmpDistSqr;
// // 								int MinPathNum;
// // 								int SurfNum = constraintSegmentSurfNums[OutConstrinedSegmentIndices[SphereElems[ti].at(ni)]];
// // 								for (int ii = 0; ii < 2; ++ii) {
// // 									TmpDistSqr = DistSqr(IntSurfRingCagePaths[SurfNum][ii].XYZAt(-1), MidPt);
// // 									if (TmpDistSqr < MinDistSqr) {
// // 										MinDistSqr = TmpDistSqr;
// // 										MinPathNum = ii;
// // 									}
// // 								}
// // 								GradPaths.push_back(ConcatenateResample({ GradPaths[SphereElems[ti][ni]], IntSurfRingCagePaths[SurfNum][MinPathNum] }, 100));
// // 								tri->at(ni) = GradPaths.size() - 1;
// // 							}
// // 
// // 							break;
// // 						}
// // 						else if (NodeTypes[SphereElems[ti].at(c2)] == NETypeC) {
// // 							/*
// // 							 *	Just need to check the two ring paths to see which terminates
// // 							 *	closer to the neighbor GP
// // 							 */
// // 							double MinDistSqr = DBL_MAX;
// // 							int SurfNum = constraintSegmentSurfNums[OutConstrinedSegmentIndices[SphereElems[ti].at(c1)]];
// // 							int minI;
// // 							for (int ii = 0; ii < 2; ++ii) {
// // 								double TmpDistSqr = DistSqr(IntSurfRingCagePaths[SurfNum][ii].XYZAt(-1), GradPaths[SphereElems[ti].at(c2)].XYZAt(-1));
// // 								if (TmpDistSqr < MinDistSqr) {
// // 									MinDistSqr = TmpDistSqr;
// // 									minI = ii;
// // 								}
// // 							}
// // 
// // 							//  							IntSurfRingCagePaths[SurfNum][minI].SaveAsOrderedZone(logStr + "C neighbor c2 ring path");
// // 							//  							TecUtilZoneSetActive(Set(TecUtilDataSetGetNumZones()).getRef(), AssignOp_MinusEquals);
// // 
// // 							GradPaths.push_back(ConcatenateResample({ GradPaths[SphereElems[ti].at(c1)], IntSurfRingCagePaths[SurfNum][minI] }, 100));
// // 							tri->at(c1) = GradPaths.size() - 1;
// // 							//  							GradPaths.back().SaveAsOrderedZone(logStr + "C neighbor c2 new path");
// // 							//  							TecUtilZoneSetActive(Set(TecUtilDataSetGetNumZones()).getRef(), AssignOp_MinusEquals);
// // 
// // 							break;
// // 						}
// // 						else continue;
// // 					}
// // 				}
// // 			}
// // 		}
// // 
// // #ifdef SUPERDEBUG
// // 		for (int ti : ElemTodo) {
// // #else
// // 		for (int ti = 0; ti < SphereElems.size(); ++ti) {
// // #endif
// // 			/*
// // 			 * Second pass over triangle corners to setup the B type
// // 			 * nodes now that we have the C and R type paths.
// // 			 */
// // 			auto * tri = &NewElemGradPathInd[ti];
// // 			for (int c1 = 0; c1 < 3; ++c1) {
// // 				if (NodeTypes[SphereElems[ti].at(c1)] == NETypeB) {
// // 					/*
// // 					 *	Type B or RB found, so need to  determine where to
// // 					 *	seed gradient paths in the interatomic surface for
// // 					 *	the bond path. We'll get the GPs for the
// // 					 *	other two corners of the triangle and find the closest point
// // 					 *	to the bond point and use that to determine the seed direction.
// // 					 */
// // 					GradPath_c const * BPPtr = &IntBondPathSegments[SphereElems[ti].at(c1)];
// // 					int BondCPNum;
// // 					for (int bi = 0; bi < CPs.NumBonds(); ++bi) {
// // 						if (sum(CPs.GetXYZ(CPTypeNum_Bond, bi) == BPPtr->XYZAt(0)) == 3) {
// // 							BondCPNum = CPs.GetTotOffsetFromTypeNumOffset(CPTypeNum_Bond, bi);
// // 							break;
// // 						}
// // 					}
// // 
// // 					int NewGPNums[2];
// // 
// // 					for (int ei = 1; ei <= 2; ++ei) {
// // 						int c2 = (c1 + ei) % 3;
// // 
// // 
// // 						string logStr = "elem " + to_string(ti) + ", node " + to_string(SphereElems[ti].at(c1)) + "-" + to_string(SphereElems[ti].at(c2)) + ": ";
// // 
// // 						//  						GradPath_c(*BPPtr).SaveAsOrderedZone(logStr + "bond path used", Red_C);
// // 						//  						TecUtilZoneSetActive(Set(TecUtilDataSetGetNumZones()).getRef(), AssignOp_MinusEquals);
// // 						//  						GradPath_c(GradPaths[tri->at(c2)]).SaveAsOrderedZone(logStr + "neighbor path used");
// // 						//  						TecUtilZoneSetActive(Set(TecUtilDataSetGetNumZones()).getRef(), AssignOp_MinusEquals);
// // 
// // 						vec3 EigenValues;
// // 						mat33 EigenVectors;
// // 						CalcEigenSystemForPoint(BPPtr->XYZAt(0), EigenValues, EigenVectors, MR);
// // 
// // 						vec3 ClosestPt = GradPaths[tri->at(c2)].ClosestPoint(BPPtr->XYZAt(0)),
// // 							// 							v1 = BPPtr->XYZAt(1) - BPPtr->XYZAt(0),
// // 							v1 = EigenVectors.col(2) + BPPtr->XYZAt(0),
// // 							v2 = ClosestPt - BPPtr->XYZAt(0),
// // 							v3 = cross(v1, v2);
// // 
// // 						// 						for (auto v : { v1,v2,v3 }){
// // 						// 							SaveVec3VecAsScatterZone({
// // 						// 								ClosestPt,
// // 						// 								ClosestPt + normalise(v)
// // 						// 								},
// // 						// 								logStr + "v123 first pass");
// // 						// 							SetZoneStyle();
// // 						// 						}
// // 						// 
// // 						// 						SaveVec3VecAsScatterZone({ ClosestPt }, logStr + "closest point on neighbor path");
// // 						// 						SetZoneStyle();
// // 
// // 						// 						v1 = normalise(v1);
// // 						// 						v3 = normalise(v3);
// // 						v2 = normalise(cross(v1, v3));
// // 
// // 						// 						for (auto v : { v1,v2,v3 }) {
// // 						// 							SaveVec3VecAsScatterZone({
// // 						// 								ClosestPt,
// // 						// 								ClosestPt + normalise(v)
// // 						// 								},
// // 						// 								logStr + "v123 second pass");
// // 						// 							SetZoneStyle();
// // 						// 						}
// // 
// // 						vec3 SeedPt = BPPtr->XYZAt(0) + v2 * 0.01;
// // 						if (DistSqr(SeedPt, ClosestPt) > DistSqr(BPPtr->XYZAt(0), ClosestPt))
// // 							SeedPt = BPPtr->XYZAt(0) + v2 * (-0.01);
// // 
// // 						// 						vec3 SeedPt = v1 * 0.01;
// // 						// 						double SeedTheta = PI * 0.5;
// // 						// 						vec4 TmpVec4 = join_cols(SeedPt, ones<vec>(1));
// // 						// 						TmpVec4 = RotationMatrix(SeedTheta, v3) * TmpVec4;
// // 						// 						SeedPt = vec3(TmpVec4.subvec(0, 2)) + BPPtr->XYZAt(0);
// // 
// // 						// 						SaveVec3VecAsScatterZone({ SeedPt }, logStr + "seed point");
// // 						// 						SetZoneStyle();
// // 
// // 						FESurface_c const * SurfPtr = nullptr;
// // 						GradPath_c const * RPPtr = nullptr;
// // 						if (NodeTypes[SphereElems[ti][c2]] == NETypeR) {
// // 							ElemTypes[ti] = NETypeRB;
// // 							int SurfNum = constraintSegmentSurfNums[OutConstrinedSegmentIndices[SphereElems[ti][c2]]];
// // 							SurfPtr = &IntSurfs[SurfNum];
// // 							/*
// // 							 * Need to get the correct ring path segment to finish off the GP
// // 							 */
// // 							for (int ii = 0; ii < 2; ++ii) {
// // 								if (approx_equal(IntSurfRingCagePaths[SurfNum][ii].XYZAt(-1), GradPaths[tri->at(c2)].XYZAt(-1), "absdiff", 0.01)) {
// // 									RPPtr = &IntSurfRingCagePaths[SurfNum][ii];
// // 									break;
// // 								}
// // 							}
// // 							if (RPPtr == nullptr) {
// // 								TecUtilDialogErrMsg("Failed to find ring cage path to complete atom-bond-ring-cage path");
// // 							}
// // 
// // 							//  							GradPath_c(*RPPtr).SaveAsOrderedZone(logStr + "ring path used", Cyan_C);
// // 							//  							TecUtilZoneSetActive(Set(TecUtilDataSetGetNumZones()).getRef(), AssignOp_MinusEquals);
// // 
// // 							//  							FESurface_c(*SurfPtr).SaveAsTriFEZone({ 1,2,3 }, logStr + "ring surface used");
// // 						}
// // 
// // 						GradPath_c GP(SeedPt, StreamDir_Reverse, 100, GPType_Classic, GPTerminate_AtCP, nullptr, &CPs, &TermRadius, &TermVal, VolInfo, HessPtrs, GradPtrs, RhoPtr, SurfPtr);
// // 						// 						GradPath_c GP(SeedPt, StreamDir_Reverse, 100, GPType_Classic, GPTerminate_AtCP, nullptr, &CPs, &TermRadius, &TermVal, VolInfo, HessPtrs, GradPtrs, RhoPtr, nullptr);
// // 						GP.SetStartEndCPNum(BondCPNum, 0);
// // 						GP.Seed();
// // 						//  						GradPath_c(GP).SaveAsOrderedZone(logStr + "new path", Green_C);
// // 						//  						TecUtilZoneSetActive(Set(TecUtilDataSetGetNumZones()).getRef(), AssignOp_MinusEquals);
// // 
// // 						vector<GradPath_c> GPList = { GP.SubGP(1, -1) };
// // 
// // 						string tmpStr;
// // 						if (NodeTypes[SphereElems[ti][c2]] == NETypeR) {
// // 							GPList.push_back(*RPPtr);
// // 							//  							tmpStr += "+ ring-cage";
// // 							//  							ConcatenateResample(GPList, 100).SaveAsOrderedZone(logStr + "new " + tmpStr + " path", Cyan_C);
// // 							//  							TecUtilZoneSetActive(Set(TecUtilDataSetGetNumZones()).getRef(), AssignOp_MinusEquals);
// // 						}
// // 						GPList.push_back(*BPPtr);
// // 						// 						tmpStr += "+ bond-nuclear";
// // 
// // 						GP = ConcatenateResample(GPList, 100);
// // 						//  						GP.SaveAsOrderedZone(logStr + "new " + tmpStr + " path", Red_C);
// // 						//  						TecUtilZoneSetActive(Set(TecUtilDataSetGetNumZones()).getRef(), AssignOp_MinusEquals);
// // 
// // 						if (GP.RhoAt(0) < GP.RhoAt(-1))
// // 							GP.Reverse();
// // 
// // 						GradPaths.push_back(GP);
// // 						int GPInd = GradPaths.size() - 1;
// // 						NewGPNums[ei - 1] = GPInd;
// // 						// 						if (ei == 1)
// // 						// 							tri->at(c1) = GPInd;
// // 						// 						else
// // 						// 							tri->insert(tri->begin() + c1, GPInd);
// // 
// // 						// 							GP.SaveAsOrderedZone(string("B path for elem " + to_string(ti) + ", node " + to_string(tri->back())), Red_C);
// // 
// // 						BNodeBondCPNums[GPInd] = BondCPNum;
// // 						BNodeSeedPts[GPInd] = SeedPt;
// // 						BNodeBondPathPtrs[GPInd] = BPPtr;
// // 						// 						BNodeRingPathPtrs[GPInd] = RPPtr;
// // 						// 						BNodeRingSurfPtrs[GPInd] = SurfPtr;
// // 					}
// // 
// // 					tri->insert(tri->begin() + c1, NewGPNums[1]);
// // 					tri->at(c1 + 1) = NewGPNums[0];
// // 
// // 					// 					for (auto ni : NewElemGradPathInd[ti]){
// // 					// 						if (GradPaths[ni].IsMade()) {
// // 					// 							GradPaths[ni].SaveAsOrderedZone("elem " + to_string(ti), Blue_C);
// // 					// 							TecUtilZoneSetActive(Set(TecUtilDataSetGetNumZones()).getRef(), AssignOp_MinusEquals);
// // 					// 						}
// // 					// 						else
// // 					// 							TecUtilDialogErrMsg(string("grad path(s) for elem " + to_string(ti) + " not made").c_str());
// // 					// 					}
// // 
// // 					break;
// // 				}
// // 			}
// // 		}
// // 
// // 		// 		NumNodes = GradPaths.size();
// // 		// #pragma omp parallel for schedule(dynamic)
// // 		// 		for (int ni = 0; ni < NumNodes; ++ni){
// // 		// 			if (BNodeSeedPts.count(ni) > 0 && !GradPaths[ni].IsMade()){
// // 		// 				int ThreadNum = omp_get_thread_num();
// // 		// 				GradPaths[ni].SetupGradPath(BNodeSeedPts[ni], StreamDir_Reverse, 100, GPType_Classic, GPTerminate_AtCP, nullptr, &CPs, &TermRadius, &TermVal, ThVolInfo[ThreadNum], HessPtrs, GradPtrs, RhoPtr, BNodeRingSurfPtrs[ni]);
// // 		// 				GradPaths[ni].SetStartEndCPNum(BNodeBondCPNums[ni], 0);
// // 		// 				GradPaths[ni].Seed();
// // 		// 				GradPath_c const *BPPtr = BNodeBondPathPtrs[ni],
// // 		// 					*RPPtr = BNodeRingPathPtrs[ni];
// // 		// 				GradPaths[ni].PointPrepend(BPPtr->XYZAt(0), BPPtr->RhoAt(0));
// // 		// 				GradPaths[ni] += *BPPtr;
// // 		// 				if (RPPtr != nullptr)
// // 		// 					GradPaths[ni] += *RPPtr;
// // 		// 				if (GradPaths[ni].RhoAt(0) < GradPaths[ni].RhoAt(-1))
// // 		// 					GradPaths[ni].Reverse();
// // 		// 			}
// // 		// 		}
// // 
// // 
// // 		// 		TecUtilLockFinish(AddOnID);
// // 		// 		return;
// // 
// // 			   /*
// // 				* We now have all the gradient paths completed in a single vector,
// // 				* and a representation of the gradient bundles in terms of the
// // 				* gradient path indices in the correct order.
// // 				* Now simply loop over each triangle, creating a gradient bundle for it.
// // 				*/
// // 
// // 		vector<FESurface_c> GradientBundles(NewElemGradPathInd.size());
// // 		NumElems = NewElemGradPathInd.size();
// // 
// // 		vector<vector<GradPath_c*> > GPPtrList(NewElemGradPathInd.size());
// // 
// // 		double const DivergentGPMaxTerminalDotProduct = cos(170. * PI / 180.);
// // 
// // // 		int NumGPs = GradPaths.size();
// // // #pragma omp parallel for schedule(dynamic)
// // // 		for (int ni = 0; ni < NumGPs; ++ni) {
// // // 			int ThNum = omp_get_thread_num();
// // // 			GradPaths[ni].ReinterpolateRhoValuesFromVolume(&ThVolInfo[ThNum]);
// // // 		}
// //  
// //  #ifdef SUPERDEBUG
// //  		for (int ti : ElemTodo) {
// //  #else
// //  #pragma omp parallel for schedule(dynamic)
// //  		for (int ti = 0; ti < NumElems; ++ti) {
// //  #endif
// //  			vector<GradPath_c*> GPPtrs;
// // 			for (auto ni : NewElemGradPathInd[ti]) {
// // 				GPPtrs.push_back(&GradPaths[ni]);
// // 			}
// //  			GradientBundles[ti].MakeFromGPs(GPPtrs, true, true);
// //  
// //  			GPPtrList[ti] = GPPtrs;
// //  		}
// //  
// //  		vector<string> ElemTypeStrs = { "C","R","Bo","RB" };
// //  #ifdef SUPERDEBUG
// //  		for (int ti : ElemTodo) {
// //  #else
// //  		for (int ti = 0; ti < NumElems; ++ti) {
// //  #endif
// //  			string tmpStr = ElemTypeStrs[ElemTypes[ti]] + "_elem " + to_string(ti) + ", ";
// //  //  			for (int ni = 0; ni < GPPtrList[ti].size(); ++ni) {
// //  //  				if (GPPtrList[ti][ni]->IsMade()) {
// //  //  					GPPtrList[ti][ni]->SaveAsOrderedZone(tmpStr + "corner " + to_string(ni));
// //  //  					SetZoneStyle({}, ZoneStyle_Path, Green_C, 0.4);
// //  //  					TecUtilZoneSetActive(Set(TecUtilDataSetGetNumZones()).getRef(), AssignOp_MinusEquals);
// //  //  				}
// //  //  			}
// //  			if (GradientBundles[ti].IsMade()) {
// //  				GradientBundles[ti].SaveAsTriFEZone({ 1,2,3 }, tmpStr + "GB");
// // 				SetZoneStyle({}, ZoneStyle_Surface, ElemTypes[ti] >= NETypeR ? ColorIndex_t(ti % 10) : White_C);
// //  				TecUtilZoneSetActive(Set(TecUtilDataSetGetNumZones()).getRef(), AssignOp_MinusEquals);
// // 				TecUtilZoneSetContour(SV_SHOW, Set(TecUtilDataSetGetNumZones()).getRef(), 0.0, FALSE);
// //  			}
// // 
// // 			// Check that no GPs are deviating and that rho for GPs is monotonically decreasing
// // 			for (int gpi = 0; gpi < GPPtrList[ti].size(); ++gpi) {
// // 				GradPath_c *GP1 = GPPtrList[ti][gpi], *GP2 = GPPtrList[ti][(gpi + 1) % GPPtrList[ti].size()];
// // 				if (GP1->IsMade() && GP2->IsMade()){
// // 					if (Distance(GP1->XYZAt(-1), GP2->XYZAt(-1)) > 0.5 && abs(dot(GP1->XYZAt(-1) - GP1->XYZAt(-1), GP2->XYZAt(-1) - GP2->XYZAt(-1))) < DivergentGPMaxTerminalDotProduct){
// // 						TecUtilDialogErrMsg(string("deviating GPs for element " + to_string(ti)).c_str());
// // 					}
// // 				}
// // 				else{
// // 					TecUtilDialogErrMsg(string("unmade GPs for element " + to_string(ti)).c_str());
// // 				}
// // 				for (int i = 0; i < GP1->GetCount() - 1; ++i){
// // 					if (GP1->RhoAt(i) < GP1->RhoAt(i+1)){
// // 						string tmpStr = "rho not monotonically decreasing for GP" + to_string(gpi) + " of element " + to_string(ti) + ": ";
// // 						tmpStr += to_string(GP1->RhoAt(i)) + " at point " + to_string(i) + ", " + to_string(GP1->RhoAt(i + 1)) + " at point " + to_string(i + 1);
// // 						GP1->SaveAsOrderedZone(tmpStr);
// // // 						TecUtilDialogErrMsg(tmpStr.c_str());
// // 						break;
// // 					}
// // 				}
// // 			}
// //  		}
// //  	}
// //  	
// // 	
// 
// 
// 	/*
// 	 *	Test GSL differentiation for grad and hessian
// 	 */
// 
// 	
// 
// 	TecUtilLockFinish(AddOnID);
// 	return;
// }
// 

struct GSLPartialDerivativeParams
{
	GSLPartialDerivativeParams(){}
	GSLPartialDerivativeParams(vec3 const * _xyz, FieldDataPointer_c const * _FuncFDPointer, VolExtentIndexWeights_s * _VolInfo)
	{
		xyz = _xyz;
		FuncFDPointer = _FuncFDPointer;
		VolInfo = _VolInfo;
	}

	int firstDir = -1;
	int secondDir = -1;
	vec3 const * xyz = nullptr;
	FieldDataPointer_c const * FuncFDPointer = nullptr;
	VolExtentIndexWeights_s * VolInfo = nullptr;
};

vec3 GradAtPoint(vec3 const & Pt, FieldDataPointer_c const & FDPtr, VolExtentIndexWeights_s & VolInfo);
double GradDirAtPoint(vec3 const & Pt, int dir, FieldDataPointer_c const & FDPtr, VolExtentIndexWeights_s & VolInfo);

double FunctionAtPointForGSLDeriv(double x, void * params)
{
	auto p = reinterpret_cast<GSLPartialDerivativeParams*>(params);

	vec3 pt = *(p->xyz);
	pt[p->firstDir] = x;

	double out = ValAtPointByPtr(pt, *(p->VolInfo), *(p->FuncFDPointer));

	return out;
}

double FunctionGradientAtPointForGSLDeriv(double x, void * params)
{
	// returns the firstDir component of the derivative w.r.t.
	// secondDir.
	auto p = reinterpret_cast<GSLPartialDerivativeParams*>(params);

	vec3 pt = *(p->xyz);
	pt[p->secondDir] = x;

	double out = GradDirAtPoint(pt, p->firstDir, *(p->FuncFDPointer), *(p->VolInfo));

	return out;
}

double GradDirAtPoint(vec3 const & Pt, int dir, FieldDataPointer_c const & FDPtr, VolExtentIndexWeights_s & VolInfo)
{
	REQUIRE(dir >= 0 && dir < 3);

	double Grad;
	GSLPartialDerivativeParams p(&Pt, &FDPtr, &VolInfo);

	gsl_function F;
	F.function = &FunctionAtPointForGSLDeriv;

	double h = 1e-8, abserror;
	double boundaryCheckDist = h * 5;

	p.firstDir = dir;
	F.params = reinterpret_cast<void*>(&p);
	if (VolInfo.IsPeriodic
		|| (Pt[dir] - VolInfo.MinXYZ[dir] > boundaryCheckDist
			&& VolInfo.MaxXYZ[dir] - Pt[dir] > boundaryCheckDist
			))
	{
		// Interior of system, so use central divided difference
		gsl_deriv_central(&F, Pt[dir], h, &Grad, &abserror);
	}
	else if (Pt[dir] - VolInfo.MinXYZ[dir] < boundaryCheckDist)
	{
		// Against lower bound of system, so use forward divided difference
		gsl_deriv_forward(&F, Pt[dir], h, &Grad, &abserror);
	}
	else
	{
		// Against upper bound of system, so use backward divided difference
		gsl_deriv_backward(&F, Pt[dir], h, &Grad, &abserror);
	}

	return Grad;
}

vec3 GradAtPoint(vec3 const & Pt, FieldDataPointer_c const & FDPtr, VolExtentIndexWeights_s & VolInfo)
{
	vec3 Grad;

	for (int i = 0; i < 3; ++i) {
		Grad[i] = GradDirAtPoint(Pt, i, FDPtr, VolInfo);
	}

	return Grad;
}

double HessDirAtPoint(vec3 const & Pt, int idir, int jdir, vector<FieldDataPointer_c> const & FDPtrs, VolExtentIndexWeights_s & VolInfo)
{
	// If 3 FDPtrs provided, assume that's a precalculated gradient,
	// or assume for the function if 1 provided.
	REQUIRE(FDPtrs.size() == 1 || FDPtrs.size() == 3);

	double Out;

	gsl_function F;

	GSLPartialDerivativeParams p;
	p.xyz = &Pt;
	p.VolInfo = &VolInfo;

	if (FDPtrs.size() == 3){
		// Use the idir component of the provided gradient as the function to 
		// differentiate in the jdir direction.
		F.function = &FunctionAtPointForGSLDeriv;
		p.FuncFDPointer = &FDPtrs[idir];
		p.firstDir = jdir;
		F.params = reinterpret_cast<void*>(&p);

		double h = 1e-8, abserror;
		double boundaryCheckDist = h * 3;

		if (VolInfo.IsPeriodic
			|| (Pt[jdir] - VolInfo.MinXYZ[jdir] > boundaryCheckDist
				&& VolInfo.MaxXYZ[jdir] - Pt[jdir] > boundaryCheckDist
				))
		{
			// Interior of system, so use central divided difference
			gsl_deriv_central(&F, Pt[jdir], h, &Out, &abserror);
		}
		else if (Pt[jdir] - VolInfo.MinXYZ[jdir] < boundaryCheckDist)
		{
			// Against lower bound of system, so use forward divided difference
			gsl_deriv_forward(&F, Pt[jdir], h, &Out, &abserror);
		}
		else
		{
			// Against upper bound of system, so use backward divided difference
			gsl_deriv_backward(&F, Pt[jdir], h, &Out, &abserror);
		}
	}
	else{
		double h = 1e-3, abserror;
		double boundaryCheckDist = h * 3;

		if (idir == jdir) {
			// For second derivative, avoid calculating the first derivative and 
			// just use a second-order approximation.
			
			vec3 Pt1 = Pt;
			vector<double> fVals;

			if (VolInfo.IsPeriodic
				|| (Pt[jdir] - VolInfo.MinXYZ[jdir] > boundaryCheckDist
					&& VolInfo.MaxXYZ[jdir] - Pt[jdir] > boundaryCheckDist
					))
			{
				// Interior of system, so use central divided difference
				for (int i = -1; i <= 1; ++i){
					Pt1[idir] = Pt[idir] + (double)i * h;
					fVals.push_back(ValAtPointByPtr(Pt1, VolInfo, FDPtrs[0]));
				}
				Out = (fVals[2] - 2. * fVals[1] + fVals[0]) / (h * h);
			}
			else if (Pt[jdir] - VolInfo.MinXYZ[jdir] < boundaryCheckDist)
			{
				// Against lower bound of system, so use forward divided difference
				for (int i = 0; i <= 2; ++i) {
					Pt1[idir] = Pt[idir] + (double)i * h;
					fVals.push_back(ValAtPointByPtr(Pt1, VolInfo, FDPtrs[0]));
				}
				Out = (fVals[2] - 2. * fVals[1] + fVals[0]) / (h * h);
			}
			else
			{
				// Against upper bound of system, so use backward divided difference
				for (int i = -2; i <= 0; ++i) {
					Pt1[idir] = Pt[idir] + (double)i * h;
					fVals.push_back(ValAtPointByPtr(Pt1, VolInfo, FDPtrs[0]));
				}
				Out = (fVals[2] - 2. * fVals[1] + fVals[0]) / (h * h);
			}
		}
		else {
			// Differentiate the idir component of the GSL derivative of the function 
			// in the jdir direction.
			F.function = &FunctionGradientAtPointForGSLDeriv;
			p.FuncFDPointer = &FDPtrs[0];
			p.firstDir = idir;
			p.secondDir = jdir;
			F.params = reinterpret_cast<void*>(&p);

			if (VolInfo.IsPeriodic
				|| (Pt[jdir] - VolInfo.MinXYZ[jdir] > boundaryCheckDist
					&& VolInfo.MaxXYZ[jdir] - Pt[jdir] > boundaryCheckDist
					))
			{
				// Interior of system, so use central divided difference
				gsl_deriv_central(&F, Pt[jdir], h, &Out, &abserror);
			}
			else if (Pt[jdir] - VolInfo.MinXYZ[jdir] < boundaryCheckDist)
			{
				// Against lower bound of system, so use forward divided difference
				gsl_deriv_forward(&F, Pt[jdir], h, &Out, &abserror);
			}
			else
			{
				// Against upper bound of system, so use backward divided difference
				gsl_deriv_backward(&F, Pt[jdir], h, &Out, &abserror);
			}
		}
	}

	return Out;
}

mat33 HessAtPoint(vec3 const & Pt, vector<FieldDataPointer_c> const & FDPtrs, VolExtentIndexWeights_s & VolInfo)
{
	mat33 Hess;

	for (int i = 0; i < 3; ++i){
		for (int j = 0; j < 3; ++j){
			if (i <= j){
				Hess.at(i, j) = HessDirAtPoint(Pt, i, j, FDPtrs, VolInfo);
			}
			else{
				Hess.at(i, j) = Hess.at(j, i);
			}
		}
	}

	return Hess;
}

void TestGSLPartialDifferentiation(){
	/*
	 *	Test GSL partial differentiation.
	 *
	 *	GSL only does 1d functions, so in order to do gradient we'll need
	 *	to construct a function for GSL where the 'params' contain not only
	 *	all the information necessary to calculate rho, but which dimension
	 *	to differentiate in. To get the partial w.r.t., e.g., x, we'll provide
	 *	the y and z coordinates of the point in 'params' so that the GSL routine
	 *	can only find the x derivative.
	 *
	 *	The same process can then be used to get the hessian by repeating the process
	 *	for each of the three gradient terms.
	 */

	 // Calc gradient of rho
	vec3 Pt = { 2,2,0 };
	// 	Pt = { 0,0,0 };

	int RhoVarVum = 4;
	int VolZoneNum = 1;
	vector<int> XYZVarNums = { 1,2,3 };

	VolExtentIndexWeights_s VolInfo;
	GetVolInfo(VolZoneNum, XYZVarNums, FALSE, VolInfo);

	FieldDataPointer_c RhoPtr;
	RhoPtr.InitializeReadPtr(VolZoneNum, RhoVarVum);

	vec3 grad = GradAtPoint(Pt, RhoPtr, VolInfo);

	// calc Hessian of rho
	mat33 hess = HessAtPoint(Pt, { RhoPtr }, VolInfo);

	// vals Hessian from precalculated gradient
	int GradXVarNum = 5;
	vector<FieldDataPointer_c> GradPtrs(3);
	for (int i = 0; i < 3; ++i)
		GradPtrs[i].InitializeReadPtr(1, GradXVarNum + i);
	mat33 hess2 = HessAtPoint(Pt, GradPtrs, VolInfo);

	mat33 hessdiff = hess2 - hess;

	return;
}

void TestGradPathPointFromRhoValue(){
	int GPZoneNum = 243, OtherGPZoneNum = 184;
	double rhoVal = 0.14;
	vec3 CPPos = { 0.58140290519360649,3.0061544952456658,-2.4207129874246442 };
	double TermRadius = 0.05;
	double TermVal = 0.001;

	GradPath_c GP(GPZoneNum, { 1,2,3,4 }, AddOnID), GP2(OtherGPZoneNum, { 1,2,3,4 }, AddOnID);

	vec3 Pt, Pt2;
	int Ind1, Ind2;
	double Weight;

	GP.GetDeviationMidpointAsTerminalAngleAndDistanceFactorFromOtherGP(GP2, DefaultGPDeviationCheckAngleFactor, DefaultGPDeviationCheckDistFactor, DefaultGPDeviationStartPtLengthFactor, Pt, Ind1, Ind2, Pt2);

	VolExtentIndexWeights_s VolInfo;
	GetVolInfo(1, { 1,2,3 }, FALSE, VolInfo);

	FieldDataPointer_c RhoPtr;
	RhoPtr.InitializeReadPtr(1, 4);

	GP.SetGPTermType(GPTerminate_AtPoint);
	GP.SetGPType(GPType_Classic);
	GP.SetNumGPPoints(1000);
	GP.SetTermPoint(CPPos);
	GP.SetVolInfo(VolInfo);
	GP.SetRhoPtr(RhoPtr);
	GP.SetTermPointRadius(TermRadius);
	GP.SetTermValue(TermVal);

	GradPath_c MidGP({ &GP,&GP2 }, Pt, StreamDir_Reverse, { std::make_pair(0, MIN(Ind1 + 10, GP.GetCount() - 1)), std::make_pair(0, MIN(Ind2 + 10, GP2.GetCount() - 1)) }, &CPPos);

	MidGP.SaveAsOrderedZone("Mid path");

	return;
}

void TestAdaptivePathResampling(){
	GradPath_c GP(20, { 1,2,3,4 }, AddOnID);

	TecUtilZoneSetActive(Set(20).getRef(), AssignOp_MinusEquals);

	GP.Resample(100);

	GP.SaveAsOrderedZone("resampled to 400 points", Blue_C);
}

void TestMaxGPDistFromNeighborGPs(){
	GradPath_c GP1(840, { 1,2,3,4 }, AddOnID);
	GradPath_c GP2(842, { 1,2,3,4 }, AddOnID);
	GradPath_c GP3(844, { 1,2,3,4 }, AddOnID);

	vector<int> IndList;
	int Ind;
	double MaxDist = GP2.GetMaxSeparationFromNeighboringGPs(vector<GradPathBase_c const *>({ &GP1, &GP3 }), IndList, Ind);

	return;
}

void GBAExtractSphereContourPoints(){

}

void TestFunction() {
	// 	TestGradPathPointFromRhoValue();
	// 	TestAdaptivePathResampling();
	// 	bool DoQuit = false;
	// 	while (!DoQuit) {
	// 		TestMaxGPDistFromNeighborGPs();
	// 	}
	// 	ArgList args;
	// 	args.appendInt(SV_CONTLINECREATEMODE, ContLineCreateMode_OneZonePerIndependentPolyline);
	// 	args.appendInt(SV_CONTLINECREATEMODE, ContLineCreateMode_OneZonePerContourLevel);
	// 	TecUtilCreateContourLineZonesX(args.getRef());
	// 	
		// Kludge INS values on GBA spheres (i.e. if an element on the sphere has a ridiculous value, set
		// it to the average of it's neighborhood.
// 	int NumZones = TecUtilDataSetGetNumZones();
// 	int NumVars = TecUtilDataSetGetNumVars();
// 	for (int z = 1; z <= NumZones; ++z) {
// 		if (AuxDataZoneItemMatches(z, CSMAuxData.GBA.ZoneType, CSMAuxData.GBA.ZoneTypeSphereZone)) {
// 			FESurface_c Sphere(z, { 1,2,3 });
// 			Sphere.GenerateElemConnectivity();
// 			int NumElems = Sphere.GetElemListPtr()->size();
// 			FieldDataPointer_c Ptr;
// 
// 			char * tmpCStr;
// 			if (TecUtilZoneGetName(z, &tmpCStr)) {
// 				string ZoneName = tmpCStr;
// 				TecUtilStringDealloc(&tmpCStr);
// 
// 				vector<vector<int> > const * ElemConnList = Sphere.GetElemConnectivityListPtr();
// 
// 				for (int v = 1; v <= NumVars; ++v) {
// 					if (TecUtilVarGetName(v, &tmpCStr)) {
// 						string VarName = tmpCStr;
// 						TecUtilStringDealloc(&tmpCStr);
// 						if (VarName.find("I: ") != string::npos
// 							|| VarName.find("INS: ") != string::npos) {
// 							TecUtilDataLoadBegin();
// 							Ptr.InitializeWritePtr(z, v);
// 							vec Vals(NumElems);
// 							for (int i = 0; i < NumElems; ++i)
// 								Vals[i] = Ptr[i];
// 
// 							// Get outer fences of the data to determine outliers
// 							vec ValsSorted = sort(Vals);
// 							double Median = median(ValsSorted),
// 								Q1 = ValsSorted[NumElems / 4],
// 								Q3 = ValsSorted[3 * NumElems / 4],
// 								InterquartileRange = Q3 - Q1,
// 								OuterFenceDist = InterquartileRange * 1000.0;
// 
// 
// 							int NumChanged = 0;
// 							vector<int> ChangedElemNums;
// 
// 							vector<bool> ElemIsOutlier(NumElems, false);
// 							for (int i = 0; i < NumElems; ++i) {
// 								ElemIsOutlier[i] = (abs(Vals[i] - Median) > OuterFenceDist);
// 
// 								for (int i = 0; i < NumElems; ++i) {
// 									if (ElemIsOutlier[i]) {
// 										NumChanged++;
// 										ChangedElemNums.push_back(i + 1);
// 
// 										int NumNeighbors = 0;
// 										vector<double> NeighborVals;
// 										for (int j = 0; j < ElemConnList->at(i).size(); ++j)
// 											if (!ElemIsOutlier[ElemConnList->at(i)[j]])
// 												NeighborVals.push_back(Vals[ElemConnList->at(i)[j]]);
// 										// 									TecUtilDataValueSetByZoneVar(z, v, i + 1, mean(NeighborVals));
// 										Ptr.Write(i, mean(vec(NeighborVals)));
// 									}
// 								}
// 								Ptr.Close();
// 								TecUtilDataLoadEnd();
// 								TecUtilStateChanged(StateChange_VarsAltered, (ArbParam_t)Set(v).getRef());
// 
// 								if (NumChanged > 0) {
// 									string ChangedString = "";
// 									for (int i = 0; i < ChangedElemNums.size(); ++i) {
// 										ChangedString += to_string(ChangedElemNums[i]);
// 										if (i < ChangedElemNums.size() - 1) ChangedString += ", ";
// 									}
// 									TecUtilDialogMessageBox(string("Zone: " + ZoneName
// 										+ "\nVar: " + VarName
// 										+ "\nMin/Median/Max = " + to_string(ValsSorted[0]) + "/" + to_string(Median) + "/" + to_string(ValsSorted[NumElems - 1])
// 										+ "\nQ1/Q3/IQR = " + to_string(Q1) + "/" + to_string(Q3) + "/" + to_string(InterquartileRange)
// 										+ "\nNumElemsChanged: " + to_string(NumChanged)
// 										+ "\nChanged elems: " + ChangedString).c_str(), MessageBox_Information);
// 								}
// 							}
// 						}
// 					}
// 				}
// 			}
// 		}
// 	}
// 	
// 	// testing gradient bundle integration on preexisting GPs.
// 	vector<vector<int> > GPStartEndZoneNums = {
// 		{
// 			3091,3113,38 //Bond-ring-ring
// 		},
// 		{
// 			3237,3265,38 //Bond-ring-cage
// 		},
// 		{
// 			15724,15735 //ring-ring-ring
// 		},
// 		{
// 			3214,3228 //ring-ring-cage
// 		},
// 		{
// 			3267,3288 //ring-cage-cage
// 		},
// 		{
// 			3975,3980 // tet cage-cage-cage 
// 		},
// 		{
// 			8267,8272 // oct cage-cage-cage early max area ind
// 		},
// 		{
// 			9002,9007 // oct cage-cage-cage late max area ind
// 		}
// 	};
// 
// 	VolExtentIndexWeights_s VolInfo;
// 	GetVolInfo(TecUtilDataSetGetNumZones(), { 1,2,3 }, FALSE, VolInfo);
// 
// 	FieldDataPointer_c RhoPtr;
// 	RhoPtr.InitializeReadPtr(TecUtilDataSetGetNumZones(), 4);
// 
// 	CritPoints_c CPs(2, { 1,2,3 }, 5, 4);
// 
// 	for (auto zStartEnd : GPStartEndZoneNums){
// 
// 		vector<GradPath_c> GPs;
// 		vector<GradPath_c const *> GPPtrs;
// 		for (int z = zStartEnd[0]; z <= zStartEnd[1]; ++z){
// 			GPs.push_back(GradPath_c(z, { 1,2,3,4 }, AddOnID));
// 		}
// 
// 		GPPtrs.resize(GPs.size());
// 		for (int i = 0; i < GPs.size(); ++i)
// 			GPPtrs[i] = &GPs[i];
// 
// 		FESurface_c GB;
// 		GB.MakeFromGPs(GPPtrs, true, true);
// 
// 		vector<double> IntVals(2, 0.0);
// 
// 		vec3 BondCP;
// 		vec3 * BondCPPtr = nullptr;
// 		if (zStartEnd.size() > 2) {
// 			BondCP = CPs.GetXYZ(zStartEnd[2]);
// 			BondCPPtr = &BondCP;
// 		}
// 
// 		NewIntegrateUsingIsosurfaces2(GPPtrs, 2, VolInfo, { RhoPtr }, IntVals, BondCPPtr, &AddOnID);
// 	}
// 	

	// Testing triangle badness
// vec3 p1 = { 1,1,1 },
// p2 = { 1,1,2 },
// p3 = { 2,1,1 };
// 
// double badness = TriBadness(p1, p2, p3);
// 
// p1 = { 1,1,1 };
// p2 = { 1,1,2 };
// p3 = { 1.05,1,1 };
// 
// badness = TriBadness(p1, p2, p3);
// 
// p1 = { 0,0,0 };
// p2 = { 1,0,0 };
// p3 = { 0.5,0.1,0 };
// 
// double dist = PointLineDist(p3, p1, p2);
// 
// p3[1] += 0.2;
// 
// dist = PointLineDist(p3, p1, p2);
// 
// tpcsm::Vec3 x0(0, 0, 0), x1(1, 0, 0), x2(0.5, 0.1, 0);
// dist = PointLineDist(x2, x0, x1);
// 
// x2.y() += .2;
// 
// dist = PointLineDist(x2, x0, x1);
// 
// vec3 p = ProjectPointToLine(p3, p1, p2);

// testing symmetry mirroring of ijk-ordered rectilinear volumes

int ZoneNum = 1;
string OriginStr = "0,0,0";
vec3 Origin = vec(SplitStringDbl(OriginStr));
vector<int> XYZVarNums = { 1,2,3 };

Set DeleteZones(ZoneNum);

VolExtentIndexWeights_s VolInfo;
GetVolInfo(ZoneNum, XYZVarNums, FALSE, VolInfo);
ZoneNum = VolumeZoneMirrorPlane(ZoneNum, 0, Origin, VolInfo, XYZVarNums);
DeleteZones += ZoneNum;
GetVolInfo(ZoneNum, XYZVarNums, FALSE, VolInfo);
VolumeZoneMirrorPlane(ZoneNum++, 1, Origin, VolInfo, XYZVarNums);
DeleteZones += ZoneNum;
GetVolInfo(ZoneNum, { 1,2,3 }, FALSE, VolInfo);
VolumeZoneMirrorPlane(ZoneNum, 2, Origin, VolInfo, XYZVarNums);

return;
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

	CSMGUILock();

	
	GradPath_c GP(Pt, GPToolData.Dir, 1000, GPType_Classic, GPTerminate_AtBoundary, nullptr, nullptr, nullptr, nullptr, GPToolData.VolInfo, GPToolData.HessPtrs, GPToolData.GradPtrs, GPToolData.RhoPtr);
	GP.Seed(false);
	GP.SaveAsOrderedZone();

	TecUtilRedraw(TRUE);


	CSMGUIUnlock();


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

	Boolean_t IsOk = GPToolData.RhoPtr.InitializeReadPtr(VolZoneNum, RhoVarNum);

	if (IsOk && Fields[fNum++].GetReturnBool()){
		GPToolData.GradPtrs.resize(3);
		for (int i = 0; i < 3; ++i){
			if (!GPToolData.GradPtrs[i].InitializeReadPtr(VolZoneNum, Fields[fNum].GetReturnInt() + i)){
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
			if (!GPToolData.HessPtrs[i].InitializeReadPtr(VolZoneNum, Fields[fNum].GetReturnInt() + i)){
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

/*
 *	Duplicates an IJK-ordered volume zone representing a
 *	rectilinear volume by mirroring it across its `PlaneNum`th basis vector
 *	with respect to `Origin`.
 *	Returns zone number of newly created zone.
 *	
 *	The input zone is assumed to be such that the origin lies at, or offset from,
 *	the point at I=J=K=1. This means if `PlaneNum` is 0, then new points are added
 *	for I less than 1 going in the negative direction, (0, -1, -2, -3) ect., and 
 *	the new values at I=J=K=1 will be the value at I=Mx, J=K=1.
 */
int VolumeZoneMirrorPlane(int ZoneNum, int PlaneNum, vec3 Origin, VolExtentIndexWeights_s VolInfo, vector<int> XYZVarNums)
{
	REQUIRE(ZoneNum > 0 && ZoneNum <= TecUtilDataSetGetNumZones());
	REQUIRE(PlaneNum >= 0 && PlaneNum <= 2);
	REQUIRE(XYZVarNums.size() == 3);
	for (auto XYZVarNum : XYZVarNums){
		REQUIRE(XYZVarNum > 0 && XYZVarNum <= TecUtilDataSetGetNumVars());
	}
	
	auto OriginFraction = VolInfo.XYZ_to_Fractional(Origin);

	auto NewIJK = VolInfo.MaxIJK, IJK = VolInfo.MaxIJK;
	vector<int> Offset = { 0,0,0 };
	NewIJK[PlaneNum] *= 2;
	if (abs(OriginFraction[PlaneNum]) < 1e-5){
		NewIJK[PlaneNum] -= 1;
		Offset[PlaneNum]++;
	}

	char *ZoneNameCStr;
	TecUtilZoneGetName(ZoneNum, &ZoneNameCStr);
	string ZoneName = ZoneNameCStr;
	TecUtilStringDealloc(&ZoneNameCStr);
	ZoneName += " " + to_string(PlaneNum) + "-mirrored";

	vector<FieldDataType_e> VarTypes;
	VarTypes.reserve(TecUtilDataSetGetNumVars());
	vector<FieldDataPointer_c> ReadPtrs(TecUtilDataSetGetNumVars()), WritePtrs(TecUtilDataSetGetNumVars());

	TecUtilDataLoadBegin();

	for (int v = 1; v <= TecUtilDataSetGetNumVars(); ++v){
		ReadPtrs[v - 1].InitializeReadPtr(ZoneNum, v);
		VarTypes.push_back(ReadPtrs[v - 1].FDType());
	}

	int NewZoneNum = -1;
	if (TecUtilDataSetAddZone(ZoneName.c_str(), NewIJK[0], NewIJK[1], NewIJK[2], ZoneType_Ordered, VarTypes.data())){
		NewZoneNum = TecUtilDataSetGetNumZones();
	}
	else{
		TecUtilDialogErrMsg("Failed to create new zone for symmetry mirroring!");
		return -1;
	}

	for (int v = 1; v <= TecUtilDataSetGetNumVars(); ++v) {
		WritePtrs[v - 1].InitializeWritePtr(NewZoneNum, v);
	}

	FieldVecPointer_c XYZPtr(vector<FieldDataPointer_c>({ ReadPtrs[XYZVarNums[0] - 1], ReadPtrs[XYZVarNums[1] - 1], ReadPtrs[XYZVarNums[2] - 1] })),
		NewXYZPtr(vector<FieldDataPointer_c>({ WritePtrs[XYZVarNums[0] - 1], WritePtrs[XYZVarNums[1] - 1], WritePtrs[XYZVarNums[2] - 1] }));

#ifndef _DEBUG
#pragma omp parallel for
#endif
	for (int i = 1; i <= IJK[0]; ++i){
		for (int j = 1; j <= IJK[1]; ++j){
			for (int k = 1; k <= IJK[2]; ++k){
				int ReadIndex = IndexFromIJK(i, j, k, IJK[0], IJK[1], IJK[2], FALSE) - 1;
				vec3 PtFrac = VolInfo.XYZ_to_Fractional(XYZPtr[ReadIndex]);
				vec3 PtNewFrac = PtFrac;
				PtNewFrac[PlaneNum] = OriginFraction[PlaneNum] - (PtFrac - OriginFraction)[PlaneNum];
				vec3 PtNew = VolInfo.Fractional_to_XYZ(PtNewFrac);
				
				vector<int> PtOldIJK = { i,j,k }, 
					PtNewIJK = PtOldIJK;

				if (Offset[PlaneNum] && PtOldIJK[PlaneNum] == 1){
					// a point on the symmetry plane, so just write it to its new position
					PtOldIJK[PlaneNum] = IJK[PlaneNum];
					int PtOldIndex = IndexFromIJK(PtOldIJK[0], PtOldIJK[1], PtOldIJK[2], NewIJK[0], NewIJK[1], NewIJK[2], FALSE) - 1;
					for (int vi = 0; vi < ReadPtrs.size(); ++vi) {
						WritePtrs[vi].Write(PtOldIndex, ReadPtrs[vi][ReadIndex]);
					}
				}
				else {
					PtNewIJK[PlaneNum] = IJK[PlaneNum] + 1 - PtOldIJK[PlaneNum];// -Offset[PlaneNum];
					PtOldIJK[PlaneNum] += IJK[PlaneNum] - Offset[PlaneNum];
					int PtOldIndex = IndexFromIJK(PtOldIJK[0], PtOldIJK[1], PtOldIJK[2], NewIJK[0], NewIJK[1], NewIJK[2], FALSE) - 1,
						PtNewIndex = IndexFromIJK(PtNewIJK[0], PtNewIJK[1], PtNewIJK[2], NewIJK[0], NewIJK[1], NewIJK[2], FALSE) - 1;
					// First copy the original point to its old and new positions
					for (int vi = 0; vi < ReadPtrs.size(); ++vi) {
						WritePtrs[vi].Write(PtOldIndex, ReadPtrs[vi][ReadIndex]);
						if (ReadPtrs[vi].VarNum() != XYZVarNums[0] && ReadPtrs[vi].VarNum() != XYZVarNums[1] && ReadPtrs[vi].VarNum() != XYZVarNums[2]) {
							WritePtrs[vi].Write(PtNewIndex, ReadPtrs[vi][ReadIndex]);
						}
					}
					// Then write the new XYZ values
					NewXYZPtr.Write(PtNewIndex, PtNew);
				}
			}
		}
	}

// 	for (int i = 1 + Offset[0]; i <= IJK[0]; ++i) {
// 		for (int j = 1 + Offset[1]; j <= IJK[1]; ++j) {
// 			for (int k = 1 + Offset[2]; k <= IJK[2]; ++k) {
// 				int ReadIndex = IndexFromIJK(i, j, k, IJK[0], IJK[1], IJK[2], FALSE) - 1;
// 				vec3 PtFrac = VolInfo.XYZ_to_Fractional(XYZPtr[ReadIndex]);
// 				vec3 PtNewFrac = PtFrac;
// 				PtNewFrac[PlaneNum] = OriginFraction[PlaneNum] - (PtFrac - OriginFraction)[PlaneNum];
// 				vec3 PtNew = VolInfo.Fractional_to_XYZ(PtNewFrac);
// 
// 				vector<int> PtOldIJK = { i,j,k },
// 					PtNewIJK = PtOldIJK;
// 				PtNewIJK[PlaneNum] = IJK[PlaneNum] + 1 - PtOldIJK[PlaneNum] - Offset[PlaneNum];
// 				PtOldIJK[PlaneNum] += IJK[PlaneNum] - Offset[PlaneNum];
// 				int PtOldIndex = IndexFromIJK(PtOldIJK[0], PtOldIJK[1], PtOldIJK[2], NewIJK[0], NewIJK[1], NewIJK[2], FALSE) - 1,
// 					PtNewIndex = IndexFromIJK(PtNewIJK[0], PtNewIJK[1], PtNewIJK[2], NewIJK[0], NewIJK[1], NewIJK[2], FALSE) - 1;
// 				// First copy the original point to its old and new positions
// 				for (int vi = 0; vi < ReadPtrs.size(); ++vi) {
// 					WritePtrs[vi].Write(PtOldIndex, ReadPtrs[vi][ReadIndex]);
// 					WritePtrs[vi].Write(PtNewIndex, ReadPtrs[vi][ReadIndex]);
// 				}
// 				// Then write the new XYZ values
// 				NewXYZPtr.Write(PtNewIndex, PtNew);
// 			}
// 		}
// 	}

	TecUtilDataLoadEnd();

	return NewZoneNum;
}

void SymmetryMirrorReturnUserInfo(bool const GuiSuccess,
	vector<GuiField_c> const & Fields,
	vector<GuiField_c> const PassthroughFields) {
	if (!GuiSuccess) return;


	int fNum = 0;

	int VolZoneNum = Fields[fNum++].GetReturnInt();
	vector<bool> MirrorXYZ(3);
	vector<int> XYZVarNums(3);
	for (int i = 0; i < 3; ++i){
		MirrorXYZ[i] = Fields[fNum++].GetReturnBool();
		XYZVarNums[i] = Fields[fNum++].GetReturnInt();
	}
// 	int XVarNum = Fields[fNum++].GetReturnInt();
	string OriginStr = Fields[fNum++].GetReturnString();
	auto OriginVec = SplitStringDbl(OriginStr);
	if (OriginVec.size() != 3){
		TecUtilDialogErrMsg("Provide the symmetry origin as three comma-delimited numbers e.g. \'0.5,-0.5,1.5\'");
		SymmetryMirrorGetUserInfo();
		return;
	}
	vec3 Origin = vec(OriginVec);
// 	vector<int> XYZVarNums = { 1,2,3 };
// 	for (int i = 0; i < 3; ++i){
// 		XYZVarNums[i] = XVarNum + i;
// 	}

	int ZoneNum = VolZoneNum;

	Set DeleteZones;

	TecUtilLockStart(AddOnID);
	CSMGUILock();

	VolExtentIndexWeights_s VolInfo;
	for (int i = 0; i < 3; ++i) {
		if (MirrorXYZ[i]) {
			GetVolInfo(ZoneNum, XYZVarNums, FALSE, VolInfo);
			DeleteZones += ZoneNum;
			ZoneNum = VolumeZoneMirrorPlane(ZoneNum, i, Origin, VolInfo, XYZVarNums);
		}
	}
// 	GetVolInfo(ZoneNum, XYZVarNums, FALSE, VolInfo);
// 	ZoneNum = VolumeZoneMirrorPlane(ZoneNum, 1, Origin, VolInfo, XYZVarNums);
// 	DeleteZones += ZoneNum;
// 	GetVolInfo(ZoneNum, { 1,2,3 }, FALSE, VolInfo);
// 	ZoneNum = VolumeZoneMirrorPlane(ZoneNum, 2, Origin, VolInfo, XYZVarNums);

	TecUtilZoneSetActive(Set(ZoneNum).getRef(), AssignOp_PlusEquals);

	char *ZoneNameCStr;
	TecUtilZoneGetName(VolZoneNum, &ZoneNameCStr);
	TecUtilZoneRename(ZoneNum, ZoneNameCStr);
	TecUtilStringDealloc(&ZoneNameCStr);

	if (!DeleteZones.isEmpty()) {
		TecUtilDataSetDeleteZone(DeleteZones.getRef());
	}

	SetZoneNum();

	CSMGUIUnlock();
	TecUtilLockFinish(AddOnID);

	return;
}

void SymmetryMirrorGetUserInfo() {
	vector<GuiField_c> Fields = {
		GuiField_c(Gui_ZoneSelect, "Zone to mirror", CSMZoneName.FullVolume),
		GuiField_c(Gui_Toggle, "Mirror in \"X\" direction", "1"),
		GuiField_c(Gui_VarSelect, "X variable", "X"),
		GuiField_c(Gui_Toggle, "Mirror in \"Y\" direction", "1"),
		GuiField_c(Gui_VarSelect, "Y variable", "Y"),
		GuiField_c(Gui_Toggle, "Mirror in \"Z\" direction", "1"),
		GuiField_c(Gui_VarSelect, "Z variable", "Z"),
		GuiField_c(Gui_String, "Symmetry origin", "0,0,0")
	};

	CSMGui("Symmetry mirror volume zone", Fields, SymmetryMirrorReturnUserInfo, AddOnID);
}


void TranslationalCopyReturnUserInfo(bool const GuiSuccess,
	vector<GuiField_c> const & Fields,
	vector<GuiField_c> const PassthroughFields) {
	if (!GuiSuccess) return;

	// Get system info and user input

	int fNum = 0;

	int VolZoneNum = Fields[fNum++].GetReturnInt();
	vector<int> MirrorXYZ(3);
	vector<int> XYZVarNums(3);
	XYZVarNums[0] = Fields[fNum++].GetReturnInt();
	XYZVarNums[1] = XYZVarNums[0] + 1;
	XYZVarNums[2] = XYZVarNums[0] + 2;
	fNum += 2;
	for (int i = 0; i < 3; ++i) {
		MirrorXYZ[i] = Fields[fNum++].GetReturnInt();
		MirrorXYZ[i] = MAX(abs(MirrorXYZ[i]), 1);
	}

	if (MirrorXYZ == vector<int>({ 1,1,1 }) || TecUtilZoneGetType(VolZoneNum) != ZoneType_Ordered){
		return;
	}

	int ZoneNum = VolZoneNum;

	Set DeleteZones(ZoneNum);

	TecUtilLockStart(AddOnID);
	TecUtilDataLoadBegin();

	VolExtentIndexWeights_s VolInfo;
	GetVolInfo(ZoneNum, XYZVarNums, FALSE, VolInfo);

	char *ZoneNameCStr;
	TecUtilZoneGetName(ZoneNum, &ZoneNameCStr);
	string ZoneName = ZoneNameCStr;
	TecUtilStringDealloc(&ZoneNameCStr);
	ZoneName += " (" + to_string(MirrorXYZ[0]) + " " + to_string(MirrorXYZ[1]) + " " + to_string(MirrorXYZ[2]) + ") supercell";

	// Get old nuclear positions
	vector<AtomGroup_s> OldAtomGroups;
	vector<int> OldAtomZoneNums;
	for (int z = 1; z <= TecUtilDataSetGetNumZones(); ++z) {
		if (AuxDataZoneItemMatches(z, CSMAuxData.DL.ZoneType, CSMAuxData.DL.ZoneTypeNuclearPositions)) {
			char *zName;
			int I, J, K;
			TecUtilZoneGetName(z, &zName);
			TecUtilZoneGetIJK(z, &I, &J, &K);

			OldAtomGroups.emplace_back(zName);
			OldAtomGroups.back().Count = I;
			OldAtomZoneNums.push_back(z);
			for (int i = 0; i < 3; ++i) {
				OldAtomGroups.back().Positions[i].resize(I);
				FieldData_pa TmpRef = TecUtilDataValueGetReadableNativeRef(z, XYZVarNums[i]);
				TecUtilDataValueArrayGetByRef(TmpRef, 1, I, OldAtomGroups.back().Positions[i].data());
			}
		}
	}

	// Setup new zone and get pointers to data

	vector<int> TotIJK(3);
	int TotNumPoints = 1;
	for (int i = 0; i < 3; ++i){
		TotIJK[i] = VolInfo.MaxIJK[i] * MirrorXYZ[i] - (MirrorXYZ[i] - 1);
		TotNumPoints *= TotIJK[i];
	}

	vector<FieldDataType_e> VarTypes;
	VarTypes.reserve(TecUtilDataSetGetNumVars());
	vector<FieldDataPointer_c> ReadPtrs(TecUtilDataSetGetNumVars()), WritePtrs(TecUtilDataSetGetNumVars());


	TecUtilPleaseWait("Preparing to copy volume zone...", TRUE);

	for (int v = 1; v <= TecUtilDataSetGetNumVars(); ++v) {
		ReadPtrs[v - 1].InitializeReadPtr(ZoneNum, v);
		VarTypes.push_back(ReadPtrs[v - 1].FDType());
	}

	int NewZoneNum = -1;
	if (TecUtilDataSetAddZone(ZoneName.c_str(), TotIJK[0], TotIJK[1], TotIJK[2], ZoneType_Ordered, VarTypes.data())) {
		NewZoneNum = TecUtilDataSetGetNumZones();
	}
	else {
		TecUtilDialogErrMsg("Failed to create new zone for translational mirroring!");
		TecUtilDataLoadEnd();
		TecUtilLockFinish(AddOnID);
		return;
	}

	for (int v = 1; v <= TecUtilDataSetGetNumVars(); ++v) {
		WritePtrs[v - 1].InitializeWritePtr(NewZoneNum, v);
	}

	FieldVecPointer_c XYZPtr(vector<FieldDataPointer_c>({ ReadPtrs[XYZVarNums[0] - 1], ReadPtrs[XYZVarNums[1] - 1], ReadPtrs[XYZVarNums[2] - 1] })),
		NewXYZPtr(vector<FieldDataPointer_c>({ WritePtrs[XYZVarNums[0] - 1], WritePtrs[XYZVarNums[1] - 1], WritePtrs[XYZVarNums[2] - 1] }));

	// Copy data to new zone

	TecUtilPleaseWait("", FALSE);
	CSMGUILock();
	auto StartTime = StatusLaunch("Copying volume zone...", AddOnID, TRUE);


	Boolean_t TaskQuit = FALSE;
	auto IJK = VolInfo.MaxIJK;
	for (int ii = 0; ii < MirrorXYZ[0] && !TaskQuit; ++ii) {
		for (int i = 1; i <= IJK[0] && !TaskQuit; ++i) {
			if (MirrorXYZ[0] > 1 && ii < MirrorXYZ[0] - 1 && i == IJK[0]){
				continue;
			}
			if (!StatusUpdate(i + IJK[0] * ii, TotIJK[0], "Copying volume zone...", AddOnID, StartTime)) {
				TaskQuit = TRUE;
			}
			for (int jj = 0; jj < MirrorXYZ[1]; ++jj) {
				for (int j = 1; j <= IJK[1] && !TaskQuit; ++j) {
					if (MirrorXYZ[1] > 1 && jj < MirrorXYZ[1] - 1 && j == IJK[1]) {
						continue;
					}
					for (int kk = 0; kk < MirrorXYZ[2]; ++kk) {
						for (int k = 1; k <= IJK[2]; ++k) {
							if (MirrorXYZ[2] > 1 && kk < MirrorXYZ[2] - 1 && k == IJK[2]) {
								continue;
							}
							int ReadIndex = IndexFromIJK(i, j, k, IJK[0], IJK[1], IJK[2], TRUE) - 1;
							int WriteIndex = IndexFromIJK(i + ii * IJK[0] - ii, j + jj * IJK[1] - jj, k + kk * IJK[2] - kk,
								TotIJK[0], TotIJK[1], TotIJK[2], TRUE) - 1;
							
							vec3 PtFrac = VolInfo.XYZ_to_Fractional(XYZPtr[ReadIndex]);
							vec3 PtNewFrac = PtFrac + vec3({ double(ii), double(jj), double(kk) });
							vec3 PtNew = VolInfo.Fractional_to_XYZ(PtNewFrac);

							// First copy the original point to its old and new positions
							for (int vi = 0; vi < ReadPtrs.size(); ++vi) {
								if (ReadPtrs[vi].VarNum() != XYZVarNums[0] && ReadPtrs[vi].VarNum() != XYZVarNums[1] && ReadPtrs[vi].VarNum() != XYZVarNums[2]) {
									WritePtrs[vi].Write(WriteIndex, ReadPtrs[vi][ReadIndex]);
								}
							}
							// Then write the new XYZ values
							NewXYZPtr.Write(WriteIndex, PtNew);
						}
					}
				}
			}
		}
	}

	StatusDrop(AddOnID);

	// Close data pointers
	for (auto & p : ReadPtrs){
		p.Close();
	}
	for (auto & p : WritePtrs) {
		p.Close();
	}
	XYZPtr.Close();
	NewXYZPtr.Close();

	// copy and save nuclear coordinates
	auto NewAtomGroups = OldAtomGroups;
	for (auto & a : NewAtomGroups){
		for (int ii = 0; ii < MirrorXYZ[0] && !TaskQuit; ++ii) {
			for (int jj = 0; jj < MirrorXYZ[1]; ++jj) {
				for (int kk = 0; kk < MirrorXYZ[2]; ++kk) {
					for (int pti = 0; pti < a.Count; ++pti) {
						vec3 OldPoint;
						for (int i = 0; i < 3; ++i) {
							OldPoint[i] = a.Positions[i][pti];
						}
						vec3 PtFrac = VolInfo.XYZ_to_Fractional(OldPoint);
						vec3 PtNewFrac = PtFrac + vec3({ double(ii), double(jj), double(kk) });
						vec3 PtNew = VolInfo.Fractional_to_XYZ(PtNewFrac);
						for (int i = 0; i < 3; ++i) {
							a.Positions[i].push_back(PtNew[i]);
						}
					}
				}
			}
		}
		a.Count = a.Positions[0].size();
	}
	vector<AtomColor_s> AtomColorList;
	PopulateAtomColorList(AtomColorList);

	CreateAtomZonesFromAtomGroupList(NewAtomGroups, { "X","Y","Z" }, vector<FieldDataType_e>(), 50);

	TecUtilZoneSetActive(Set(NewZoneNum).getRef(), AssignOp_PlusEquals);

	for (auto z : OldAtomZoneNums){
		DeleteZones += z;
	}

	if (!DeleteZones.isEmpty()) {
		TecUtilDataSetDeleteZone(DeleteZones.getRef());
	}

	SetZoneNum();

	TecUtilDataLoadEnd();
	CSMGUIUnlock();
	TecUtilLockFinish(AddOnID);

	return;
}

void TranslationalCopyGetUserInfo() {
	vector<GuiField_c> Fields = {
		GuiField_c(Gui_ZoneSelect, "Zone to mirror", CSMZoneName.FullVolume),
		GuiField_c(Gui_VarSelect, "X variable", "X"),
		GuiField_c(Gui_Label, "Define the supercell lattice factors"),
		GuiField_c(Gui_Label, "((1,1,1) describes the current cell)"),
		GuiField_c(Gui_Int, "X supercell factor", "1"),
		GuiField_c(Gui_Int, "Y supercell factor", "1"),
		GuiField_c(Gui_Int, "Z supercell factor", "1")
	};

	CSMGui("Translational mirror volume zone", Fields, TranslationalCopyReturnUserInfo, AddOnID);
}

void ImportNuclearCoordinatesFromXYZ(){

	char *XYZCStr, *FileNameCStr;

// 	TecUtilDialogGetSimpleText("Enter \"XYZ\" string describing nuclear coordinates", "ElementSymbol X Y Z", &XYZCStr);
	TecUtilDialogGetFileName(SelectFileOption_ReadSingleFile, &FileNameCStr, "XYZ file", "Nuclear coordinates.xyz", "*.xyz");
	std::ifstream SS(FileNameCStr);
	TecUtilStringDealloc(&FileNameCStr);

// 	stringstream SS;
// 	SS.str(XYZCStr);
// 	TecUtilStringDealloc(&XYZCStr);

	vector<AtomGroup_s> AtomGroupList;
	std::map<string, vector<vec3> > AtomNameToPositions;
	vector<EntIndex_t> XYZVarNums;
	EntIndex_t vNum = 1;
	for (auto v : { "X","Y","Z" }){
		XYZVarNums.push_back(MAX(VarNumByName(v), 1));
	}

	TecUtilDialogGetVariables("Select variables to use for \"XYZ\" values", "Column 1", "Column 2", "Column 3", &XYZVarNums[0], &XYZVarNums[1], &XYZVarNums[2]);



	while (SS.good()){
		string Name;
		vec3 Pos;
		int i;
		SS >> Name;
		for (i = 0; i < 3 && SS.good(); ++i){
			SS >> Pos[i];
		}
		if (!SS.good() && i < 3){
			TecUtilDialogErrMsg("Failed to parse atom XYZ list");
			return;
		}
		if (AtomNameToPositions.count(Name)){
			AtomNameToPositions[Name].push_back(Pos);
		}
		else{
			AtomNameToPositions[Name] = { Pos };
		}
	}

	SS.close();

	for (auto & Atom : AtomNameToPositions){
		AtomGroupList.emplace_back(Atom.first);
		AtomGroupList.back().Count = Atom.second.size();
		for (auto & p : Atom.second){
			for (int dir = 0; dir < 3; ++dir){
				AtomGroupList.back().Positions[dir].push_back(p[dir]);
			}
		}
	}

	vector<AtomColor_s> AtomColorList;
	PopulateAtomColorList(AtomColorList);

	vector<string> XYZNameList(3);
	for (int i = 0; i < 3; ++i){
		char *TmpCStr;
		TecUtilVarGetName(XYZVarNums[i], &TmpCStr);
		XYZNameList[i] = TmpCStr;
		TecUtilStringDealloc(&TmpCStr);
	}

	vector<FieldDataType_e> VarDataTypes(TecUtilDataSetGetNumVars(), FieldDataType_Bit);
	for (auto i : XYZVarNums){
		VarDataTypes[i - 1] = FieldDataType_Float;
	}

	CreateAtomZonesFromAtomGroupList(AtomGroupList, XYZNameList, VarDataTypes, 50);
}

/*
 * 		vector<string> VarNameList = {
			"P over N",
			"P/(N alpha)",
			"S[P]",
			"S[P / alpha]",
			"S[alpha]",
			"S[P] / S[alpha]",
			"S[alpha] - S[P]"
	};
 */

vector<vec> CalculateShannonEntropiesFromDensitySolidAngleValues(vec & DensityVec, vec & SolidAngleVec){
	REQUIRE(DensityVec.size() == SolidAngleVec.size() && DensityVec.size() > 0);
	vector<vec> Out;

	DensityVec = normalise(DensityVec);
	SolidAngleVec = normalise(SolidAngleVec);

	vec val, 
		prob = DensityVec;

	// Add values to output vector according to order of VarNameList above
	Out.push_back(prob); // Pr = P over N
	val = prob / SolidAngleVec; // P/(N alpha)
	Out.push_back(val);
	Out.emplace_back(-prob * log(prob)); // S[Pr] shannon condensed entropy
	val = prob / SolidAngleVec; // S[Pr / alpha]
	Out.emplace_back(-val * log(val));
	val = SolidAngleVec; // S[alpha]
	Out.emplace_back(-val * log(val));
	Out.emplace_back((-prob * log(prob)) / (-val * log(val)));
	Out.emplace_back((-prob * log(prob)) - (-val * log(val)));

	return Out;
}

vector<vector<double> > CalculateNormalizedAndScaledVectors(vec const & InVals, vec const & InWeights) {
	vector<vector<double> > OutVals = { arma::conv_to<vector<double> >::from(InVals) };

	vec NormalisedValues = InVals / InWeights;
	double NormSum = sum(NormalisedValues);
	vec ScaledValues = NormalisedValues;
	if (NormSum != 0.0){
		ScaledValues *= (sum(InVals) / NormSum);
	}

	return { arma::conv_to<vector<double> >::from(InVals), arma::conv_to<vector<double> >::from(NormalisedValues), arma::conv_to<vector<double> >::from(ScaledValues) };
}

vector<string> CalculateShannonEntropyVarNameList = {
		"Normalized condensed density | Pr = P/N",
		"Further normalized by solid angle | Pr/alpha",
		"Shannon entropy of Pr | S[Pr]",
		"Shannon entropy of normalized Pr | S[Pr/alpha]",
		"Maximum (radial) Shannon entropy | S[alpha]",
		"Ratio | S[P] / S[alpha]",
		"Shannon entropy of interaction | S[P] - S[alpha]"
};

enum CalculateShannonEntropiesRegionTypes_e{
	Region_AtomicBasin = 0,
	Region_AllAtomicBasins,
	Region_MoleculesMG,
	Region_MoleculesOS,
	Region_TopologicalCages,
	Region_SpecialGradientBundles,
	Region_CondensedBasins,
	Region_ActiveDGBs,
	Region_NumRegions
};

void CalculateAndSaveShannonEntropiesForSphereNameElements(
	string const & VarNamePrefix,
	vector<std::pair<string, vector<int> > > & SphereNameAndElements, // empty vector<int> means use all elements
	vector<bool> const & DoVar, // list of ints of indices of VarNameList (if new vars are added, need to update here, VarNameList, and in the other Function)
	std::map<string, FESurface_c> & SpheresByName,
	std::map<string, FieldDataPointer_c> & SphereNameToVarPtr,
	std::map<string, vec> & SphereNameToElemSolidAngles)
{
	// prepare new variables if necessary (double type for sphere zones and bit for others)
// 	vector<string> CalculateShannonEntropyVarNameList = {
// 		": P over N",
// 		": P/(N alpha)",
// 		": S[P]",
// 		": S[P / alpha]",
// 		": S[alpha]",
// 		": S[P] / S[alpha]",
// 		": S[P] - S[alpha]"
// 	};

	vector<FieldDataType_e> FDTypes(TecUtilDataSetGetNumZones(), FieldDataType_Bit);
	vector<ValueLocation_e> ValueLocations(TecUtilDataSetGetNumZones(), ValueLocation_Nodal);
	for (int zi = 1; zi < TecUtilDataSetGetNumZones(); ++zi){
		if (AuxDataZoneItemMatches(zi, CSMAuxData.GBA.ZoneType, CSMAuxData.GBA.ZoneTypeSphereZone)){
			FDTypes[zi - 1] = FieldDataType_Float;
			ValueLocations[zi - 1] = ValueLocation_CellCentered;
		}
	}

	vector<int> VarNums(CalculateShannonEntropyVarNameList.size() * 3, -1);
	vector<string> VarPrefixes = { "I","N","S" };
	for (int vi = 0; vi < CalculateShannonEntropyVarNameList.size(); ++vi){
		if (DoVar[vi]) {
			for (int vj = 0; vj < 3; ++vj) {
				string VarName = "";
				for (int vk = 0; vk <= vj; ++vk){
					VarName += VarPrefixes[vk];
				}
				VarName += ": " + VarNamePrefix + CalculateShannonEntropyVarNameList[vi];
				VarNums[(vi * 3) + vj] = VarNumByName(VarName);
				if (VarNums[(vi * 3) + vj] <= 0) {
					ArgList Args;
					Args.appendString(SV_NAME, VarName);
					Args.appendArray(SV_VARDATATYPE, FDTypes.data());
					Args.appendArray(SV_VALUELOCATION, ValueLocations.data());

					if (TecUtilDataSetAddVarX(Args.getRef())) {
						VarNums[(vi * 3) + vj] = TecUtilDataSetGetNumVars();
					}
					else {
						TecUtilDialogErrMsg(("Failed to make new variable: " + VarNamePrefix + CalculateShannonEntropyVarNameList[vi]).c_str());
						return;
					}
				}
			}
		}
	}

	// Prepare values and calculate output values
	vector<double> DensityVals, SolidAngleVals;
	for (auto & s : SphereNameAndElements) {
		auto VarPtr = SphereNameToVarPtr[s.first];
		auto SolidAngleVec = SphereNameToElemSolidAngles[s.first];
		if (s.second.empty()){
			s.second.resize(VarPtr.Size());
			for (int ei = 0; ei < VarPtr.Size(); ++ei) {
				s.second[ei] = ei;
			}
		}

		vector<double> tmpVals(s.second.size());
		for (int ei = 0; ei < s.second.size(); ++ei) {
			tmpVals[ei] = VarPtr[s.second[ei]];
		}
		DensityVals.insert(DensityVals.end(), tmpVals.begin(), tmpVals.end());

		for (int ei = 0; ei < s.second.size(); ++ei) {
			tmpVals[ei] = SolidAngleVec[s.second[ei]];
		}

		SolidAngleVals.insert(SolidAngleVals.end(), tmpVals.begin(), tmpVals.end());
	}

	// Now calculate the actual shannon entropy values for the specified spheres and elements


	vector<vec> EntropyValues = CalculateShannonEntropiesFromDensitySolidAngleValues(vec(DensityVals), vec(SolidAngleVals));
	vector<vector<double> > EntropyValuesFull;
	for (auto v : EntropyValues){
		auto NewVals = CalculateNormalizedAndScaledVectors(v, vec(SolidAngleVals));
		for (auto w : NewVals){
			EntropyValuesFull.push_back(w);
		}
	}

	// Then another loop over the spheres to write the new values
	int ei = 0;
	for (auto & s : SphereNameAndElements){
		vector<FieldDataPointer_c> VarPtrs(EntropyValuesFull.size());
		for (int vi = 0; vi < VarPtrs.size(); ++vi){
			if (VarNums[vi] > 0) {
				VarPtrs[vi].InitializeWritePtr(SpheresByName[s.first].GetZoneNum(), VarNums[vi]);
			}
		}

		for (int ej = 0; ej < s.second.size(); ++ej){
			for (int vi = 0; vi < VarPtrs.size(); ++vi) {
				if (VarNums[vi] > 0) {
					VarPtrs[vi].Write(s.second[ej], EntropyValuesFull[vi][ei]);
				}
			}
			ei++;
		}
	}
}

void CollectGradientBundlesByRegionType(int PVarNum, 
	int CPZoneNum, 
	vector<int> XYZVarNums, 
	int CPTypeVarNum,
	std::map<string, FESurface_c> & SpheresByName,
	std::map<string, FieldDataPointer_c> & SphereNameToVarPtr,
	std::map<string, vec> & SphereNameToElemSolidAngles,
	std::map<string, std::set<string> > & MoleculeMGNamesToSphereNames,
	std::map<string, std::set<string> > & MoleculeOSNamesToSphereNames,
	std::map<string, vector<int> > & CondensedBasinNameToSphereElems,
	std::map<string, string> & CondensedBasinNameToSphereName,
	std::map<string, std::set<string> > & SpecialGradientBundleNameToCondensedBasinNames,
	std::map<string, vector<int> > & SphereNameToActiveDGBElemNums,
	std::map<string, std::map<string, vector<int> > > & TopoCageNameToSphereNameToElems,
	std::map<string, std::map<string, double> > & CondensedBasinNameToBoundaryIntVarNamesToIntVals)
{
	REQUIRE(PVarNum > 0 && PVarNum <= TecUtilDataSetGetNumVars());
	REQUIRE(CPZoneNum > 0 && CPZoneNum <= TecUtilDataSetGetNumZones());

	CritPoints_c CPs(CPZoneNum, XYZVarNums, CPTypeVarNum);
	CritPoints_c CPsOneSkeleton = CPs;

	CPs.GenerateCPGraph(CSMAuxData.CC.MolecularGraphGPSubTypes);
	CPsOneSkeleton.GenerateCPGraph(CSMAuxData.CC.OneSkeletonGPSubTypes);

	//1. Get pointers to condensed density for all spheres
	// 
	// 
// 	std::map<string, FESurface_c> SpheresByName;
// 	std::map<string, FieldDataPointer_c> SphereNameToVarPtr;
// 	std::map<string, vec> SphereNameToElemSolidAngles;
// 
// 
// 	TecUtilDataLoadBegin();

	for (int zi = 1; zi <= TecUtilDataSetGetNumZones(); ++zi) {
		if (AuxDataZoneItemMatches(zi, CSMAuxData.GBA.ZoneType, CSMAuxData.GBA.ZoneTypeSphereZone)) {
			// sphere zone
			string SphereName = AuxDataZoneGetItem(zi, CSMAuxData.GBA.SourceNucleusName);
			FieldDataPointer_c VarPtr;
			VarPtr.InitializeReadPtr(zi, PVarNum);
			SphereNameToVarPtr[SphereName] = VarPtr;

			FESurface_c Sphere(zi, XYZVarNums);
			auto ElemAreas = Sphere.TriSphereElemSolidAngles();
			SphereNameToElemSolidAngles[SphereName] = ElemAreas;
			SpheresByName[SphereName] = Sphere;
		}
	}

	// Get "molecules" by searching for edges in the CP graph between the spheres present.
	// Because this is using map and set, duplicates are automatically ignored.
// 	std::map<string, std::set<string> > MoleculeMGNamesToSphereNames, MoleculeOSNamesToSphereNames;
	for (auto & si : SpheresByName) {
		std::set<string> MoleculeNameMG, MoleculeNameOS;
		int CPi = stoi(AuxDataZoneGetItem(si.second.GetZoneNum(), CSMAuxData.GBA.SphereCPNum))-1;
		for (auto & sj : SpheresByName) {
			int CPj = stoi(AuxDataZoneGetItem(sj.second.GetZoneNum(), CSMAuxData.GBA.SphereCPNum))-1;
			if (CPi != CPj) {
				if (CPs.HasEdge(CPi, CPj, -1)) {
					MoleculeNameMG.insert(si.first);
					MoleculeNameOS.insert(si.first);
					MoleculeNameMG.insert(sj.first);
					MoleculeNameOS.insert(sj.first);
				}
				else if (CPsOneSkeleton.HasEdge(CPi, CPj, -1)) {
					MoleculeNameOS.insert(si.first);
					MoleculeNameOS.insert(sj.first);
				}
			}
		}
		if (!MoleculeNameMG.empty()) {
			string MoleculeName = StringJoin(vector<string>(MoleculeNameMG.begin(), MoleculeNameMG.end()), " -- ");
			MoleculeMGNamesToSphereNames[MoleculeName] = MoleculeNameMG;
		}
		if (!MoleculeNameOS.empty()) {
			string MoleculeName = StringJoin(vector<string>(MoleculeNameOS.begin(), MoleculeNameOS.end()), " -- ");
			MoleculeOSNamesToSphereNames[MoleculeName] = MoleculeNameOS;
		}
	}

	//2. Get max/min basin information
	//	1. keep map for max/min basins <CondensedBasinInfo, BasinElementsVals>
	//	2. another map for special gradient bundles <CondensedBasinName, BasinElements> that will include the elements of all max/min basins that contribute to them
	//	3. another map for molecular ensembles defined in two ways
// 			1. by the molecular graph(i.e.bond paths)
// 			2. by the 1 - skeleton(i.e.by any special gradient paths)

// 	std::map<string, vector<int> > CondensedBasinNameToSphereElems;
// 	std::map<string, string> CondensedBasinNameToSphereName;
// 	std::map<string, std::set<string> > SpecialGradientBundleNameToCondensedBasinNames;
	for (int zi = 1; zi <= TecUtilDataSetGetNumZones(); ++zi) {
		if (AuxDataZoneItemMatches(zi, CSMAuxData.GBA.ZoneType, CSMAuxData.GBA.ZoneTypeCondensedAttractiveBasin) || AuxDataZoneItemMatches(zi, CSMAuxData.GBA.ZoneType, CSMAuxData.GBA.ZoneTypeCondensedRepulsiveBasin)) {
			// condensed basin
			string CondensedBasinName = AuxDataZoneGetItem(zi, CSMAuxData.GBA.CondensedBasinInfo) + " from " + AuxDataZoneGetItem(zi, CSMAuxData.GBA.CondensedBasinDefiningVariable);
			string SphereName = AuxDataZoneGetItem(zi, CSMAuxData.GBA.SourceNucleusName);
			auto ElemList = SplitStringInt(AuxDataZoneGetItem(zi, CSMAuxData.GBA.CondensedBasinSphereElements)); // 0-based
			CondensedBasinNameToSphereName[CondensedBasinName] = SphereName;
			CondensedBasinNameToSphereElems[CondensedBasinName] = ElemList;

			// Get boundary element integration values
			string BoundaryIntValStr, IntVarNames;
			if (AuxDataZoneGetItem(zi, CSMAuxData.GBA.CondensedBasinBoundaryIntVals, BoundaryIntValStr) && AuxDataZoneGetItem(zi, CSMAuxData.GBA.IntVarNames, IntVarNames)) {
				auto BoundaryIntValList = SplitStringDbl(BoundaryIntValStr);
				auto BoundaryIntNameList = SplitString(IntVarNames, ",");
				std::map<string, double> IntVarNamesToVals;
				for (int vi = 0; vi < BoundaryIntValList.size(); ++vi){
					IntVarNamesToVals[BoundaryIntNameList[vi]] = BoundaryIntValList[vi];
				}
				CondensedBasinNameToBoundaryIntVarNamesToIntVals[CondensedBasinName] = IntVarNamesToVals;
			}

			// Add to special gradient bundle if present
			if (AuxDataZoneItemMatches(zi, CSMAuxData.GBA.CondensedBasinIsSpecialGradientBundle, "true")) {
				string SpecialGradientBundleName = AuxDataZoneGetItem(zi, CSMAuxData.GBA.CondensedBasinName) + " from " + AuxDataZoneGetItem(zi, CSMAuxData.GBA.CondensedBasinDefiningVariable);
				SpecialGradientBundleNameToCondensedBasinNames[SpecialGradientBundleName].insert(CondensedBasinName);
			}
		}
	}

	// get topological cages
	for (int zi = 1; zi <= TecUtilDataSetGetNumZones(); ++zi) {
		if (AuxDataZoneItemMatches(zi, CSMAuxData.GBA.ZoneType, CSMAuxData.GBA.ZoneTypeTopoCageWedge)) {
			// condensed basin
			string CageName = AuxDataZoneGetItem(zi, CSMAuxData.GBA.CondensedBasinName);
			string SphereName = AuxDataZoneGetItem(zi, CSMAuxData.GBA.SourceNucleusName);
			auto ElemList = SplitStringInt(AuxDataZoneGetItem(zi, CSMAuxData.GBA.CondensedBasinSphereElements)); // 0-based
			TopoCageNameToSphereNameToElems[CageName][SphereName] = ElemList;
		}
	}

	// get  all active DGB zones
	std::map < string, std::set<int> > SphereNameToActiveDGBElemNumSet;

	for (int zi = 1; zi <= TecUtilDataSetGetNumZones(); ++zi) {
		if (AuxDataZoneItemMatches(zi, CSMAuxData.GBA.ZoneType, CSMAuxData.GBA.ZoneTypeDGB)) {
			SphereNameToActiveDGBElemNumSet[AuxDataZoneGetItem(zi, CSMAuxData.GBA.SourceNucleusName)].insert(stoi(AuxDataZoneGetItem(zi, CSMAuxData.GBA.ElemNum)) - 1);
		}
	}
	for (auto const & s : SphereNameToActiveDGBElemNumSet) {
		SphereNameToActiveDGBElemNums[s.first] = vector<int>(s.second.cbegin(), s.second.cend());
	}
}

/*
 *	Calculate Shannon entropies for every max/min condensed basin, every special gradient bundle
 *	(unions of max/min basins), and for entire molecules (either the whole system or molecules as
 *	defined by the molecular graph).
 */
void CalculateShannonEntropies(int PVarNum, int CPZoneNum, vector<int> XYZVarNums, int CPTypeVarNum,
	vector<bool> const & DoRegion, vector<bool> const & DoVar) {

	std::map<string, FESurface_c> SpheresByName;
	std::map<string, FieldDataPointer_c> SphereNameToVarPtr;
	std::map<string, vec> SphereNameToElemSolidAngles;
	std::map<string, std::set<string> > MoleculeMGNamesToSphereNames, MoleculeOSNamesToSphereNames;
	std::map<string, vector<int> > CondensedBasinNameToSphereElems, SphereNameToActiveDGBElemNums;
	std::map<string, std::map<string, vector<int> > > TopoCageNameToSphereNameToElems;
	std::map<string, string> CondensedBasinNameToSphereName;
	std::map<string, std::set<string> > SpecialGradientBundleNameToCondensedBasinNames;
	std::map<string, std::map<string, double> > CondensedBasinNameToBoundaryIntVals;


	TecUtilDataLoadBegin();

	CollectGradientBundlesByRegionType(PVarNum, CPZoneNum, XYZVarNums, CPTypeVarNum, SpheresByName, SphereNameToVarPtr, SphereNameToElemSolidAngles, MoleculeMGNamesToSphereNames, MoleculeOSNamesToSphereNames, CondensedBasinNameToSphereElems, CondensedBasinNameToSphereName, SpecialGradientBundleNameToCondensedBasinNames, SphereNameToActiveDGBElemNums, TopoCageNameToSphereNameToElems, CondensedBasinNameToBoundaryIntVals);

	//3. For each max/min basin
	//	1. compute the Shannon entropy for the basin
	//	2. save it in a “Condensed Basin Shannon Entropy” variable
	//4. For each special gradient path
	//	1. compute the Shannon entropy for the gradient bundle
	//	2. save it to a new variable named for the gradient bundle (e.g. “bond 1”)
	//		1. only spheres/basins contributing to the special gradient bundle need have double precision, the rest can have boolean to save space.
	//5. Compute the molecular ensemble Shannon entropy
	//	1. using all spheres (molecular boundaries be damned)
	//	2. save it to a common variable for all spheres

	vector<std::pair<string, vector<int> > > SphereNameAndElements;

	// First for atomic basins and the entire system
	if (DoRegion[Region_AtomicBasin] || DoRegion[Region_AllAtomicBasins]) {
		string VarNamePrefix = "Atomic basin: ";
		for (auto & s : SpheresByName) {
			SphereNameAndElements.push_back(std::make_pair(s.first, vector<int>()));

			// The individual sphere
			if (DoRegion[Region_AtomicBasin]) {
				CalculateAndSaveShannonEntropiesForSphereNameElements(VarNamePrefix, vector<std::pair<string, vector<int> > >({ std::make_pair(s.first, vector<int>()) }), DoVar, SpheresByName, SphereNameToVarPtr, SphereNameToElemSolidAngles);
			}
		}

		// Then the whole system, if more than one atom
		if (DoRegion[Region_AllAtomicBasins]) {
			if (SphereNameAndElements.size() > 1) {
				VarNamePrefix = "All atomic basins: ";
				CalculateAndSaveShannonEntropiesForSphereNameElements(VarNamePrefix, SphereNameAndElements, DoVar, SpheresByName, SphereNameToVarPtr, SphereNameToElemSolidAngles);
			}
		}
	}
	if (DoRegion[Region_MoleculesMG]) {
		for (auto & m : MoleculeMGNamesToSphereNames) {
			SphereNameAndElements.clear();
			for (auto & s : m.second) {
				SphereNameAndElements.push_back(std::make_pair(s, vector<int>()));
			}
			CalculateAndSaveShannonEntropiesForSphereNameElements("Molecule (mol. graph): ", SphereNameAndElements, DoVar, SpheresByName, SphereNameToVarPtr, SphereNameToElemSolidAngles);
		}
	}

	if (DoRegion[Region_MoleculesOS]) {
		for (auto & m : MoleculeOSNamesToSphereNames) {
			SphereNameAndElements.clear();
			for (auto & s : m.second) {
				SphereNameAndElements.push_back(std::make_pair(s, vector<int>()));
			}
			CalculateAndSaveShannonEntropiesForSphereNameElements("Molecule (1-skeleton): ", SphereNameAndElements, DoVar, SpheresByName, SphereNameToVarPtr, SphereNameToElemSolidAngles);
		}
	}

	if (DoRegion[Region_ActiveDGBs] && !SphereNameToActiveDGBElemNums.empty()) {
		string OutPath = "";
		SphereNameAndElements.clear();
		for (auto & s : SphereNameToActiveDGBElemNums) {
			SphereNameAndElements.push_back(std::make_pair(s.first, s.second));
		}
		CalculateAndSaveShannonEntropiesForSphereNameElements("Active dGBs: ", SphereNameAndElements, DoVar, SpheresByName, SphereNameToVarPtr, SphereNameToElemSolidAngles);
	}

	if (DoRegion[Region_TopologicalCages]) {
		// Then special gradient bundles
		string OutPath = "";
		for (auto & cage : TopoCageNameToSphereNameToElems) {
			SphereNameAndElements.clear();
			for (auto & cage_sphere : cage.second) {
				SphereNameAndElements.push_back(cage_sphere);
			}
			CalculateAndSaveShannonEntropiesForSphereNameElements(cage.first + ": ", SphereNameAndElements, DoVar, SpheresByName, SphereNameToVarPtr, SphereNameToElemSolidAngles);
		}
	}

	if (DoRegion[Region_SpecialGradientBundles]) {
		// Then special gradient bundles
		for (auto & sgb : SpecialGradientBundleNameToCondensedBasinNames) {
			SphereNameAndElements.clear();
			for (auto & cb : sgb.second) {
				SphereNameAndElements.push_back(std::make_pair(CondensedBasinNameToSphereName[cb], CondensedBasinNameToSphereElems[cb]));
			}
			CalculateAndSaveShannonEntropiesForSphereNameElements(sgb.first + ": ", SphereNameAndElements, DoVar, SpheresByName, SphereNameToVarPtr, SphereNameToElemSolidAngles);
		}
	}

	if (DoRegion[Region_CondensedBasins]) {
		// Then all condensed basins
		for (auto & cb : CondensedBasinNameToSphereName) {
			CalculateAndSaveShannonEntropiesForSphereNameElements(cb.first + ": ", vector<std::pair<string, vector<int> > >({ std::make_pair(cb.second, CondensedBasinNameToSphereElems[cb.first]) }), DoVar, SpheresByName, SphereNameToVarPtr, SphereNameToElemSolidAngles);
		}
	}


	TecUtilDataLoadEnd();
}


void CalculateShannonEntropyReturnUserInfo(bool const GuiSuccess,
	vector<GuiField_c> const & Fields,
	vector<GuiField_c> const PassthroughFields) {
	if (!GuiSuccess) return;


	int fNum = 0;

	int PVarNum = Fields[fNum++].GetReturnInt();
	int CPZoneNum = Fields[fNum++].GetReturnInt();
	int CPTypeVarNum = Fields[fNum++].GetReturnInt();
	int XVarNum = Fields[fNum++].GetReturnInt();

	vector<int> XYZVarNums(3);
	for (int i = 0; i < 3; ++i){
		XYZVarNums[i] = XVarNum + i;
	}

	vector<bool> DoVar(CalculateShannonEntropyVarNameList.size()), DoRegion(Region_NumRegions);
	fNum += 2;
	for (int i = 0; i < DoVar.size(); ++i){
		DoVar[i] = Fields[fNum++].GetReturnBool();
	}
	fNum += 2;
	for (int i = 0; i < DoRegion.size(); ++i){
		DoRegion[i] = Fields[fNum++].GetReturnBool();
	}
	

	TecUtilLockStart(AddOnID);
	TecUtilPleaseWait("Calculating Shannon entropy of condensed density... Please wait.", TRUE);
	CSMGUILock();

	CalculateShannonEntropies(PVarNum, CPZoneNum, XYZVarNums, CPTypeVarNum, DoRegion, DoVar);

	CSMGUIUnlock();
	TecUtilPleaseWait("Calculating Shannon entropy of condensed density... Please wait.", FALSE);
	TecUtilDialogMessageBox("Finished", MessageBox_Information);

	TecUtilLockFinish(AddOnID);

	return;
}

void CalculateShannonEntropyGetUserInfo() {
	vector<GuiField_c> Fields = {
		GuiField_c(Gui_VarSelect, "Condensed charge density (input function)", "I: " + CSMVarName.Dens),
		GuiField_c(Gui_ZoneSelect, "Critical points zone", CSMZoneName.CriticalPoints),
		GuiField_c(Gui_VarSelect, "Critical point type variable", CSMVarName.CritPointType),
		GuiField_c(Gui_VarSelect, "X variable", "X"),
		GuiField_c(Gui_Label, " "),
		GuiField_c(Gui_Label, "Calculate the following variables:")
	};

	for (auto v : CalculateShannonEntropyVarNameList){
		Fields.emplace_back(Gui_Toggle, v, "0");
	}
	Fields.back().SetSearchString("1"); // check only the last box

	Fields.emplace_back(Gui_Label, " ");
	Fields.emplace_back(Gui_Label, "For the following regions:");
	vector<std::pair<string, string> > RegionTypeList = {
		std::make_pair("every atomic basin (atomic chart)", "1"),
		std::make_pair("the system molecular ensemble (all atomic basins together)", "1"),
		std::make_pair("every molecular ensemble (by molecular graph)", "1"),
		std::make_pair("every molecular ensemble (by 1-skeleton)", "1"),
		std::make_pair("every topological cage", "1"),
		std::make_pair("every special gradient bundle (bonds etc.)", "0"),
		std::make_pair("every max/min condensed basin", "0"),
		std::make_pair("active dGB or sphere element zones", "1")
	};
	for (auto r : RegionTypeList) {
		Fields.emplace_back(Gui_Toggle, r.first, r.second);
	}

	CSMGui("Calculate Shannon entropy of condensed density", Fields, CalculateShannonEntropyReturnUserInfo, AddOnID);
}

void ExportGBADataForRegions(
	string const & FileName,
	string const & DatasetName,
	string const & RegionName,
	string const & OutDir,
	string & OutPath,
	vector<std::pair<string, vector<int> > > & SphereNameAndElements, // empty vector<int> means use all elements
	vector<string> const & IntVarNamesIn,
	vector<int> const & IntVarNumsIn,
	bool IncludeAllGBs,
	std::map<string, FESurface_c> & SpheresByName,
	std::map<string, FieldDataPointer_c> & SphereNameToVarPtr,
	std::map<string, vec> & SphereNameToElemSolidAngles,
	string const & AllGBOutPathSuffix = "",
	std::map<string, std::map<string,double> > SphereNameToBoundaryIntTotals = std::map<string, std::map<string, double> >())
{
	// Get alphabetically sorted list of variables so files with mismatched variable sets can be merged more easily
	std::map<string, int> SortedVarNames;
	for (int vi = 0; vi < IntVarNamesIn.size(); ++vi){
		SortedVarNames[IntVarNamesIn[vi]] = IntVarNumsIn[vi];
	}

	vector<string> IntVarNames;
	vector<int> IntVarNums;
	for (auto const & v : SortedVarNames){
		IntVarNames.push_back(v.first);
		IntVarNums.push_back(v.second);
	}
	

	// update elements for whole spheres
	for (auto & s : SphereNameAndElements){
		if (s.second.empty()) {
			s.second.resize(SphereNameToVarPtr[s.first].Size());
			for (int ei = 0; ei < s.second.size(); ++ei) {
				s.second[ei] = ei;
			}
		}
	}

	// Get data
	std::map<string, vector<vec> > SphereNameToAllVals;
	std::map<string, vector<double> > SphereNameToValTotals;
	vector<double> VarTotals(IntVarNums.size(), 0.0), BoundaryVarTotals(IntVarNums.size(), 0.0);
	for (auto s : SphereNameAndElements){
		int SphereZoneNum = SphereNameToVarPtr[s.first].ZoneNum();
		for (int vi = 0; vi < IntVarNums.size(); ++vi) {
			int v = IntVarNums[vi];
			FieldDataPointer_c Ptr;
			Ptr.InitializeReadPtr(SphereZoneNum, v);

			vec Vals(s.second.size());
			for (int ei = 0; ei < s.second.size(); ++ei){
				Vals[ei] = Ptr[s.second[ei]];
			}
			SphereNameToAllVals[s.first].push_back(Vals);
			double total = sum(Vals);
			VarTotals[vi] += total;
			SphereNameToValTotals[s.first].push_back(total);
			if (!SphereNameToBoundaryIntTotals.empty()){
				BoundaryVarTotals[vi] += SphereNameToBoundaryIntTotals[s.first][IntVarNames[vi]];
			}
		}
	}

	// open output file
	bool AppendToFile = (OutPath != "");
	if (!AppendToFile){
		OutPath = OutDir + "/" + StringRemoveSubString(FileName,":") + ".csv";
	}
	ofstream OutFile(OutPath.c_str(), AppendToFile ? std::ios::app : std::ios::trunc);
	ofstream OutFileGBs;
	if (IncludeAllGBs){
		OutFileGBs = ofstream((OutDir + "/" + StringRemoveSubString(FileName + AllGBOutPathSuffix,":") + "_dGBs.csv").c_str(), AppendToFile ? std::ios::app : std::ios::trunc);
	}
	if (OutFile.is_open()){
		// write headers (region name, sphere name (or full), var names)
		if (!AppendToFile) {
			OutFile << "Dataset name,Region name,Atom name";
			for (auto const & vn : IntVarNames) {
				OutFile << "," << vn;
			}
			if (!SphereNameToBoundaryIntTotals.empty()){
				for (auto const & vn : IntVarNames) {
					OutFile << "," << vn << " (boundary)";
				}
				for (auto const & vn : IntVarNames) {
					OutFile << "," << vn << " (with boundary error)";
				}
			}
			OutFile << endl;
		}
		if (SphereNameToValTotals.size() > 1) {
			// write the full values
			OutFile << DatasetName << "," << RegionName << ",full";
			for (auto const & total : VarTotals) {
				OutFile << "," << std::setprecision(8) << total;
			}
			if (!SphereNameToBoundaryIntTotals.empty()) {
				for (auto const & total : BoundaryVarTotals) {
					OutFile << "," << std::setprecision(8) << total;
				}
				for (int vi = 0; vi < BoundaryVarTotals.size(); ++vi){
					auto & val = VarTotals[vi];
					auto & bval = BoundaryVarTotals[vi];
					OutFile << "," << std::setprecision(8) << val - 0.5 * bval << " +/- " << 0.5 * abs(bval);
				}
			}
			OutFile << endl;
		}
		// then values for each sphere
		for (auto const & s : SphereNameToValTotals){
			OutFile << DatasetName << "," << RegionName << "," << s.first;
			for (auto const & val : s.second){
				OutFile << "," << std::setprecision(8) << val;
			}
			if (!SphereNameToBoundaryIntTotals.empty()) {
				auto & bvals = SphereNameToBoundaryIntTotals[s.first];
				for (auto const & varName : IntVarNames) {
					auto & val = bvals[varName];
					OutFile << "," << std::setprecision(8) << val;
				}
				for (int vi = 0; vi < IntVarNames.size(); ++vi) {
					auto & val = s.second[vi];
					auto & bval = bvals[IntVarNames[vi]];
					OutFile << "," << std::setprecision(8) << val - 0.5 * bval << " +/- " << 0.5 * abs(bval);
				}
			}
			OutFile << endl;
		}

		OutFile.close();

		if (IncludeAllGBs && OutFileGBs.is_open()) {
			// Now the individual dGBs
// 			if (!AppendToFile) {
				// write headers (sphere name, element number, var names)
				OutFileGBs << "Dataset name,Region name,Atom name,dGB index (starts at 0)";
				for (auto const & vn : IntVarNames) {
					OutFileGBs << "," << vn;
				}
				OutFileGBs << endl;
// 			}

			// values for each dGB of each sphere
			for (auto const & s : SphereNameAndElements) {
				auto *SphereElemVals = &SphereNameToAllVals[s.first];
				for (int ei = 0; ei < s.second.size(); ++ei) {
					OutFileGBs << DatasetName << "," << RegionName << "," << s.first << "," << s.second[ei];
					for (int vi = 0; vi < IntVarNums.size(); ++vi){
						OutFileGBs << "," << std::setprecision(8) << SphereElemVals->at(vi)[ei];
					}
					OutFileGBs << endl;
				}
			}

			OutFileGBs.close();
		}
	}
}

/*
 *	Export GBA integration data for a number of regions.
 */
void ExportGBAData(int PVarNum, int CPZoneNum, vector<int> XYZVarNums, int CPTypeVarNum,
	vector<bool> const & DoRegion, bool IncludeAllDGBs, vector<bool> IncludeINS) {

	std::map<string, FESurface_c> SpheresByName;
	std::map<string, FieldDataPointer_c> SphereNameToVarPtr;
	std::map<string, vec> SphereNameToElemSolidAngles;
	std::map<string, std::set<string> > MoleculeMGNamesToSphereNames, MoleculeOSNamesToSphereNames;
	std::map<string, vector<int> > CondensedBasinNameToSphereElems, SphereNameToActiveDGBElemNums;
	std::map<string, std::map<string, vector<int> > > TopoCageNameToSphereNameToElems;
	std::map<string, string> CondensedBasinNameToSphereName;
	std::map<string, std::set<string> > SpecialGradientBundleNameToCondensedBasinNames;
	std::map<string, std::map<string,double> > CondensedBasinNameToBoundaryIntVarNamesToIntVals;


	TecUtilDataLoadBegin();

	// first get the export directory
	char* FolderNameCStr;
	if (TecUtilDialogGetFolderName("Select folder to save files", &FolderNameCStr)) {

		CollectGradientBundlesByRegionType(PVarNum, CPZoneNum, XYZVarNums, CPTypeVarNum, SpheresByName, SphereNameToVarPtr, SphereNameToElemSolidAngles, MoleculeMGNamesToSphereNames, MoleculeOSNamesToSphereNames, CondensedBasinNameToSphereElems, CondensedBasinNameToSphereName, SpecialGradientBundleNameToCondensedBasinNames, SphereNameToActiveDGBElemNums, TopoCageNameToSphereNameToElems, CondensedBasinNameToBoundaryIntVarNamesToIntVals);

		vector<int> IntVarNums;
		vector<string> IntVarNames;
		vector<string> INSStrs = { "I: ", "IN: ", "INS: " };
		vector<string> IntCheckStrs;
		for (int i = 0; i < 3; ++i) {
			if (IncludeINS[i]) {
				IntCheckStrs.push_back(INSStrs[i]);
			}
		}
		IntCheckStrs.push_back(" Integration");
		for (int VarNum = 1; VarNum <= TecUtilDataSetGetNumVars(); ++VarNum) {
			char *VarName, *CheckStr;
			if (TecUtilVarGetName(VarNum, &VarName)) {
				for (string const & Str : IntCheckStrs) {
					CheckStr = std::strstr(VarName, Str.c_str());
					if (CheckStr != nullptr) {
						string TmpStr = VarName;
						std::replace(TmpStr.begin(), TmpStr.end(), ',', '.');
						IntVarNames.push_back(TmpStr);
						IntVarNums.push_back(VarNum);
						break;
					}
					// 								break; // Only want the "I:" integration values
				}
				TecUtilStringDealloc(&VarName);
			}
		}

		// export files for each region


		string OutDir = FolderNameCStr;
		TecUtilStringDealloc(&FolderNameCStr);

		char *DataSetNameCStr;
		int junkint;
		TecUtilDataSetGetInfo(&DataSetNameCStr, &junkint, &junkint);
		string DataSetName = DataSetNameCStr;
		TecUtilStringDealloc(&DataSetNameCStr);


		vector<std::pair<string, vector<int> > > SphereNameAndElements;

		// First for atomic basins and the entire system
		if (DoRegion[Region_AtomicBasin] || DoRegion[Region_AllAtomicBasins]) {
			string OutPath = "";
			for (auto & s : SpheresByName) {
				SphereNameAndElements.push_back(std::make_pair(s.first, vector<int>()));

				// The individual sphere
				if (DoRegion[Region_AtomicBasin]) {
					ExportGBADataForRegions(DataSetName + " Atomic basins", DataSetName, s.first, OutDir, OutPath, vector<std::pair<string, vector<int> > >({ std::make_pair(s.first, vector<int>()) }), IntVarNames, IntVarNums, IncludeAllDGBs, SpheresByName, SphereNameToVarPtr, SphereNameToElemSolidAngles, "_" + s.first);
// 					CalculateAndSaveShannonEntropiesForSphereNameElements(VarNamePrefix, vector<std::pair<string, vector<int> > >({ std::make_pair(s.first, vector<int>()) }), DoVar, SpheresByName, SphereNameToVarPtr, SphereNameToElemSolidAngles);
				}
			}

			// Then the whole system, if more than one atom
			if (DoRegion[Region_AllAtomicBasins]) {
				if (SphereNameAndElements.size() > 1) {
					OutPath = "";
					ExportGBADataForRegions(DataSetName + " All atomic basins", DataSetName, "All atomic basins", OutDir, OutPath, SphereNameAndElements, IntVarNames, IntVarNums, IncludeAllDGBs, SpheresByName, SphereNameToVarPtr, SphereNameToElemSolidAngles);
// 					CalculateAndSaveShannonEntropiesForSphereNameElements(VarNamePrefix, SphereNameAndElements, DoVar, SpheresByName, SphereNameToVarPtr, SphereNameToElemSolidAngles);
				}
			}
		}
		if (DoRegion[Region_MoleculesMG]) {
			string OutPath = "";
			for (auto & m : MoleculeMGNamesToSphereNames) {
				SphereNameAndElements.clear();
				for (auto & s : m.second) {
					SphereNameAndElements.push_back(std::make_pair(s, vector<int>()));
				}
				ExportGBADataForRegions(DataSetName + " Molecules (mol. graph) ", DataSetName, m.first, OutDir, OutPath, SphereNameAndElements, IntVarNames, IntVarNums, IncludeAllDGBs, SpheresByName, SphereNameToVarPtr, SphereNameToElemSolidAngles, "_" + m.first);
// 				CalculateAndSaveShannonEntropiesForSphereNameElements(m.first + ": ", SphereNameAndElements, DoVar, SpheresByName, SphereNameToVarPtr, SphereNameToElemSolidAngles);
			}
		}

		if (DoRegion[Region_MoleculesOS]) {
			string OutPath = "";
			for (auto & m : MoleculeOSNamesToSphereNames) {
				SphereNameAndElements.clear();
				for (auto & s : m.second) {
					SphereNameAndElements.push_back(std::make_pair(s, vector<int>()));
				}
				ExportGBADataForRegions(DataSetName + " Molecules (1-skeleton) ", DataSetName, m.first, OutDir, OutPath, SphereNameAndElements, IntVarNames, IntVarNums, IncludeAllDGBs, SpheresByName, SphereNameToVarPtr, SphereNameToElemSolidAngles, "_" + m.first);
// 				CalculateAndSaveShannonEntropiesForSphereNameElements(m.first + ": ", SphereNameAndElements, DoVar, SpheresByName, SphereNameToVarPtr, SphereNameToElemSolidAngles);
			}
		}

		if (DoRegion[Region_TopologicalCages]) {
			// Then special gradient bundles
			string OutPath = "";
			for (auto & cage : TopoCageNameToSphereNameToElems) {
				SphereNameAndElements.clear();
				for (auto & cage_sphere : cage.second) {
					SphereNameAndElements.push_back(cage_sphere);
				}
				ExportGBADataForRegions(DataSetName + " Topological cages", DataSetName, cage.first, OutDir, OutPath, SphereNameAndElements, IntVarNames, IntVarNums, IncludeAllDGBs, SpheresByName, SphereNameToVarPtr, SphereNameToElemSolidAngles, "_" + cage.first);
			}
		}

		if (DoRegion[Region_SpecialGradientBundles]) {
			// Then special gradient bundles
			string OutPath = "";
			for (auto & sgb : SpecialGradientBundleNameToCondensedBasinNames) {
				SphereNameAndElements.clear();
				std::map<string, std::map<string,double> > SphereNameToBoundaryIntTotals;
				for (auto & cb : sgb.second) {
					SphereNameAndElements.push_back(std::make_pair(CondensedBasinNameToSphereName[cb], CondensedBasinNameToSphereElems[cb]));
					if (CondensedBasinNameToBoundaryIntVarNamesToIntVals.count(cb)){
						SphereNameToBoundaryIntTotals[CondensedBasinNameToSphereName[cb]] = CondensedBasinNameToBoundaryIntVarNamesToIntVals[cb];
					}
				}
				ExportGBADataForRegions(DataSetName + " Special gradient bundles", DataSetName, sgb.first, OutDir, OutPath, SphereNameAndElements, IntVarNames, IntVarNums, IncludeAllDGBs, SpheresByName, SphereNameToVarPtr, SphereNameToElemSolidAngles, "_" + sgb.first, SphereNameToBoundaryIntTotals);
// 				CalculateAndSaveShannonEntropiesForSphereNameElements(sgb.first + ": ", SphereNameAndElements, DoVar, SpheresByName, SphereNameToVarPtr, SphereNameToElemSolidAngles);
			}
		}

		if (DoRegion[Region_CondensedBasins]) {
			// Then all condensed basins
			string OutPath = "";
			for (auto & cb : CondensedBasinNameToSphereName) {
				std::map<string, std::map<string, double> > SphereNameToBoundaryIntTotals;
				if (CondensedBasinNameToBoundaryIntVarNamesToIntVals.count(cb.first)) {
					SphereNameToBoundaryIntTotals[cb.second] = CondensedBasinNameToBoundaryIntVarNamesToIntVals[cb.first];
				}
				ExportGBADataForRegions(DataSetName + " Condensed basins", DataSetName, cb.first, OutDir, OutPath, vector<std::pair<string, vector<int> > >({ std::make_pair(cb.second, CondensedBasinNameToSphereElems[cb.first]) }), IntVarNames, IntVarNums, IncludeAllDGBs, SpheresByName, SphereNameToVarPtr, SphereNameToElemSolidAngles, "_" + cb.first, SphereNameToBoundaryIntTotals);
// 				CalculateAndSaveShannonEntropiesForSphereNameElements(cb.first + ": ", vector<std::pair<string, vector<int> > >({ std::make_pair(cb.second, CondensedBasinNameToSphereElems[cb.first]) }), DoVar, SpheresByName, SphereNameToVarPtr, SphereNameToElemSolidAngles);
			}
		}

		if (DoRegion[Region_ActiveDGBs] && !SphereNameToActiveDGBElemNums.empty()) {
			string OutPath = "";
			SphereNameAndElements.clear();
			for (auto & s : SphereNameToActiveDGBElemNums) {
				SphereNameAndElements.push_back(std::make_pair(s.first, s.second));
			}
			ExportGBADataForRegions(DataSetName + " Active dGBs", DataSetName, "active dGB", OutDir, OutPath, SphereNameAndElements, IntVarNames, IntVarNums, IncludeAllDGBs, SpheresByName, SphereNameToVarPtr, SphereNameToElemSolidAngles);
		}
	}


	TecUtilDataLoadEnd();
}

void ExportGBADataReturnUserInfo(bool const GuiSuccess,
	vector<GuiField_c> const & Fields,
	vector<GuiField_c> const PassthroughFields) {
	if (!GuiSuccess) return;


	int fNum = 0;

	int PVarNum = Fields[fNum++].GetReturnInt();
	int CPZoneNum = Fields[fNum++].GetReturnInt();
	int CPTypeVarNum = Fields[fNum++].GetReturnInt();
	int XVarNum = Fields[fNum++].GetReturnInt();

	vector<int> XYZVarNums(3);
	for (int i = 0; i < 3; ++i) {
		XYZVarNums[i] = XVarNum + i;
	}
	fNum++;

	bool IncludeAllDGBs = Fields[fNum++].GetReturnBool();

	fNum++;

	vector<bool> IncludeINS(3);
	for (int i = 0; i < 3; ++i){
		IncludeINS[i] = Fields[fNum++].GetReturnBool();
	}

	fNum += 2;
	vector<bool> DoRegion(Region_NumRegions);
	for (int i = 0; i < DoRegion.size(); ++i) {
		DoRegion[i] = Fields[fNum++].GetReturnBool();
	}


	TecUtilLockStart(AddOnID);
	TecUtilPleaseWait("Exporting... Please wait.", TRUE);
	CSMGUILock();

	ExportGBAData(PVarNum, CPZoneNum, XYZVarNums, CPTypeVarNum, DoRegion, IncludeAllDGBs, IncludeINS);

	CSMGUIUnlock();
	TecUtilPleaseWait("Exporting... Please wait.", FALSE);
	TecUtilDialogMessageBox("Finished", MessageBox_Information);

	TecUtilLockFinish(AddOnID);

	return;
}

void ExportGBADataGetUserInfo() {
	vector<GuiField_c> Fields = {
		GuiField_c(Gui_VarSelect, "Condensed charge density (input function)", "I: " + CSMVarName.Dens),
		GuiField_c(Gui_ZoneSelect, "Critical points zone", CSMZoneName.CriticalPoints),
		GuiField_c(Gui_VarSelect, "Critical point type variable", CSMVarName.CritPointType),
		GuiField_c(Gui_VarSelect, "X variable", "X"),
		GuiField_c(Gui_Label, " "),
		GuiField_c(Gui_Toggle, "Export individual dGB integration values", "0"),
		GuiField_c(Gui_Label, " "),
		GuiField_c(Gui_Toggle, "Export \"I: \" integration values", "1"),
		GuiField_c(Gui_Toggle, "Export \"IN: \" (normalized by solid angle) integration values", "0"),
		GuiField_c(Gui_Toggle, "Export \"INS: \" (normalized and scaled) integration values", "0"),
		GuiField_c(Gui_Label, " "),
		GuiField_c(Gui_Label, "Export data for the following types of regions:")
	};

	vector<std::pair<string, string> > RegionTypeList = {
		std::make_pair("Atomic basin (atomic chart)", "1"),
		std::make_pair("The system molecular ensemble (all atomic basins together)", "0"),
		std::make_pair("Molecular ensemble (by molecular graph)", "0"),
		std::make_pair("Molecular ensemble (by 1-skeleton)", "0"),
		std::make_pair("Topological cages", "0"),
		std::make_pair("Special gradient bundle (bonds etc.)", "1"),
		std::make_pair("Max/min condensed basin", "1"),
		std::make_pair("Active dGB or sphere element zones", "1")
	};
	for (auto r : RegionTypeList) {
		Fields.emplace_back(Gui_Toggle, r.first, r.second);
	}
	Fields.emplace_back(Gui_Label, " ");
	Fields.emplace_back(Gui_Label, "(be sure to create a new folder)");

	CSMGui("Export gradient bundle integration data", Fields, ExportGBADataReturnUserInfo, AddOnID);
}


void MapVolumeZoneVarsToOtherZonesReturnUserInfo(bool const GuiSuccess,
	vector<GuiField_c> const & Fields,
	vector<GuiField_c> const PassthroughFields) {
	if (!GuiSuccess) return;


	int fNum = 0;

	int SourceZoneNum = Fields[fNum++].GetReturnInt();
	int XVarNum = Fields[fNum++].GetReturnInt();
	auto CopyVarNums = Fields[fNum++].GetReturnIntVec(),
		DestZoneNums = Fields[fNum++].GetReturnIntVec();

	vector<int> XYZVarNums = { XVarNum, XVarNum + 1, XVarNum + 2 };

	TecUtilLockStart(AddOnID);

	MapAllVarsToAllZones(SourceZoneNum, XYZVarNums, CopyVarNums, DestZoneNums, AddOnID);

	TecUtilLockFinish(AddOnID);

	return;
}

void MapVolumeZoneVarsToOtherZonesGetUserInfo() {
	vector<GuiField_c> Fields = {
		GuiField_c(Gui_ZoneSelect, "Source zone", "Full Volume"),
		GuiField_c(Gui_VarSelect, "X", "X"),
		GuiField_c(Gui_VarSelectMulti, "Variables to map", ""),
		GuiField_c(Gui_ZoneSelectMulti, "Destination zones", "")
	};

	CSMGui("Map volume zone variables to other zones", Fields, MapVolumeZoneVarsToOtherZonesReturnUserInfo, AddOnID);
}

void GenerateGUIBondsReturnUserInfo(bool const GuiSuccess,
	vector<GuiField_c> const & Fields,
	vector<GuiField_c> const PassthroughFields) {
	if (!GuiSuccess) return;


	int fNum = 0;

	auto SourceZoneNums = Fields[fNum++].GetReturnIntVec();
	auto XVarNum = Fields[fNum++].GetReturnInt();
	auto BondCutoff = Fields[fNum++].GetReturnDouble();
	auto DoVDW = Fields[fNum++].GetReturnBool();
	auto VDWCutoff = Fields[fNum++].GetReturnDouble();
	auto MinLineThickness = Fields[fNum++].GetReturnDouble();
	auto MaxLineThickness = Fields[fNum++].GetReturnDouble();
	auto ReplaceOldBonds = Fields[fNum++].GetReturnBool();
	auto HeteroBondsOnly = Fields[fNum++].GetReturnBool();

	vector<int> XYZVarNums = { XVarNum, XVarNum + 1, XVarNum + 2 };

	TecUtilLockStart(AddOnID);

	DrawBondLinesBetweenNuclearPositions(SourceZoneNums, XYZVarNums, BondCutoff, DoVDW ? VDWCutoff : -1.0, MinLineThickness, MaxLineThickness, ReplaceOldBonds, HeteroBondsOnly);

	TecUtilLockFinish(AddOnID);

	return;
}

void GenerateGUIBondsGetUserInfo() {
	vector<GuiField_c> Fields = {
		GuiField_c(Gui_ZoneSelectMulti, "Nuclear position zones to use", NuclearPositionsZoneNameBase),
		GuiField_c(Gui_VarSelect, "X variable", "X"),
		GuiField_c(Gui_Double, "\"Bond\" cutoff distance", "3.0"),
		GuiField_c(Gui_Toggle, "Include intermolecular bonds", "0"),
		GuiField_c(Gui_Double, "Intermolecular bond cutoff distance", "8.0"),
		GuiField_c(Gui_Double, "Minimum bond line thickness", "0.4"),
		GuiField_c(Gui_Double, "Maximum bond line thickness", "0.8"),
		GuiField_c(Gui_Toggle, "Replace existing GUI bonds where collisions occur?", "1"),
		GuiField_c(Gui_Toggle, "Heteronuclear bonds only?", "0")
	};

	CSMGui("Generate GUI bonds", Fields, GenerateGUIBondsReturnUserInfo, AddOnID);
}


void ResizeSpheresReturnUserInfo(bool const GuiSuccess,
	vector<GuiField_c> const & Fields,
	vector<GuiField_c> const PassthroughFields) {
	if (!GuiSuccess) return;

	TecUtilLockStart(AddOnID);

	int fNum = 0;

	auto ZoneNums = Fields[fNum++].GetReturnIntVec();
	auto ScaleByVar = Fields[fNum++].GetReturnBool();
	auto ScaleVarNum = Fields[fNum++].GetReturnInt();
	auto ScaleFactor = Fields[fNum++].GetReturnDouble();
	auto LogScale = Fields[fNum++].GetReturnBool();
	auto SizeFactor = Fields[fNum++].GetReturnDouble();
	auto RelativeSize = Fields[fNum++].GetReturnBool();

	for (auto z : ZoneNums) {
		ResizeSphere(z, SizeFactor, !RelativeSize, ScaleByVar, ScaleVarNum, ScaleFactor, LogScale);
	}

	TecUtilRedrawAll(TRUE);

	TecUtilLockFinish(AddOnID);
}

void ResizeSpheresGetUserInfo() {

	vector<GuiField_c> Fields = {
		GuiField_c(Gui_ZoneSelectMulti, "Sphere zones", CSMZoneName.FullVolume.substr(0, 10)),
		GuiField_c(Gui_Toggle, "Scale radius with variable", "0"),
		GuiField_c(Gui_VarSelect, "Radius scale variable", "INS: Electron Density"),
		GuiField_c(Gui_Double, "Radius scale factor", "1.0"),
		GuiField_c(Gui_Toggle, "Log scaling", "1.0"),
		GuiField_c(Gui_Double, "New sphere size", "1.0"),
		GuiField_c(Gui_Toggle, "Relative to current size", "0")
	};

	CSMGui("Resize sphere zones", Fields, ResizeSpheresReturnUserInfo, AddOnID);
}


void FindContourCurvaturePoints(vector<int> ZoneNums, vector<int> const & XYZVarNums, bool UseCurrentContours){
	if (UseCurrentContours){
		int OldNumZones = TecUtilDataSetGetNumZones();
		TecUtilMacroExecuteCommand("$!CREATECONTOURLINEZONES CONTLINECREATEMODE = ONEZONEPERINDEPENDENTPOLYLINE");
		int NewNumZones = TecUtilDataSetGetNumZones();
		ZoneNums.clear();
		ZoneNums.reserve(NewNumZones - OldNumZones - 1);
		for (int z = OldNumZones + 1; z <= NewNumZones; ++z){
			ZoneNums.push_back(z);
		}
	}
	vector<vec3> ZGCPoints, MaxCurvaturePoints, MinCurvaturePoints;
	for (auto const & z : ZoneNums){
		bool IsOk = TecUtilZoneIsOrdered(z);
		int IJK[3];
		TecUtilZoneGetIJK(z, &IJK[0], &IJK[1], &IJK[2]);
		IsOk = IsOk && IJK[1] == 1 && IJK[2] == 1;
		FieldVecPointer_c XYZPtr;
		XYZPtr.InitializeReadPtr(z, XYZVarNums);
				// these zones are i-ordered and describe a closed path, so the last point and the first point should be the same
		IsOk = IsOk && approx_equal(XYZPtr[0], XYZPtr[-1], "absdiff", 1e-14);

		if (IsOk){

			vector<vec3> ZGCPointsLocal, MaxCurvaturePointsLocal, MinCurvaturePointsLocal;

			int N = XYZPtr.Size();
			auto XYZRhoVarNums = XYZVarNums;
			XYZRhoVarNums.push_back(XYZVarNums[0]);
			GradPath_c Path(z, XYZRhoVarNums, AddOnID);
			Path = Path.SubGP(0, -2);
			Path.Resample((3 * N) / 4, GPResampleMethod_Linear);
			N = Path.GetCount();

			// first get contour plane defined by three "distant" points
			vec3 vPlane = cross(Path.XYZAt(N / 3) - Path.XYZAt(0), Path.XYZAt(N * 2 / 3) - Path.XYZAt(0));

			// setup basis spline fit
// 			const size_t n_buffer = MIN(10, N - 2);
// 			const size_t k = 4;
// 			const size_t n = N + 2 * n_buffer;
// 			const size_t ncoeffs = n / 8;
// 			const size_t nbreak = ncoeffs - 2;
// 			size_t i, j;
// 			gsl_bspline_workspace *bw, *xw, *yw, *zw;
// 			gsl_vector *B;
// 			double dy;
// 			gsl_rng *r;
// 			gsl_vector *c, *xc, *yc, *zc, *w, *wx, *wy, *wz;
// 			gsl_vector *arclen, *curvature, *x, *y, *z;
// 			gsl_matrix *X, *cov, *xcov, *ycov, *zcov;
// 			gsl_multifit_linear_workspace *mw;
// 			double chisq, Rsq, dof, tss;
// 			vector<double> arclencheck(n);
// 
// 			gsl_rng_env_setup();
// 			r = gsl_rng_alloc(gsl_rng_default);
// 
// 			/* allocate a cubic bspline workspace */
// 
// 			bw = gsl_bspline_alloc(k, nbreak);
// 			xw = gsl_bspline_alloc(k, nbreak);
// 			yw = gsl_bspline_alloc(k, nbreak);
// 			zw = gsl_bspline_alloc(k, nbreak);
// 			B = gsl_vector_alloc(ncoeffs);
// 
// 			arclen = gsl_vector_alloc(n);
// 			curvature = gsl_vector_alloc(n);
// 			x = gsl_vector_alloc(n);
// 			y = gsl_vector_alloc(n);
// 			z = gsl_vector_alloc(n);
// 			X = gsl_matrix_alloc(n, ncoeffs);
// 			c = gsl_vector_alloc(ncoeffs);
// 			xc = gsl_vector_alloc(ncoeffs);
// 			yc = gsl_vector_alloc(ncoeffs);
// 			zc = gsl_vector_alloc(ncoeffs);
// 			w = gsl_vector_alloc(n);
// 			wx = gsl_vector_alloc(n);
// 			wy = gsl_vector_alloc(n);
// 			wz = gsl_vector_alloc(n);
// 			cov = gsl_matrix_alloc(ncoeffs, ncoeffs);
// 			xcov = gsl_matrix_alloc(ncoeffs, ncoeffs);
// 			ycov = gsl_matrix_alloc(ncoeffs, ncoeffs);
// 			zcov = gsl_matrix_alloc(ncoeffs, ncoeffs);
// 			mw = gsl_multifit_linear_alloc(n, ncoeffs);



			// Now move along the contour checking the sign of curvature. 
			// Each time the sign flips, add a new point
			// 
			int LastSign;
			bool FirstLoop = true;
			vector<double> C(N), CSign(N), ArcLen(N);
			ArcLen[0] = 0.0;
			double StartArcLength;
			for (int ni = 0; ni < N; ++ni) {
				int nj = (ni + 1) % N,
					nk = (ni + 2) % N;

				vec3 v1 = Path.XYZAt(nj) - Path.XYZAt(ni),
					v2 = Path.XYZAt(nk) - Path.XYZAt(nj),
					vNormal = cross(v1, v2);
				int TmpSign = SIGN(dot(vPlane, vNormal));
				C[nj] = VectorAngle(v1, v2) / ((norm(v1) + norm(v2)) * 0.5);
				CSign[nj] = C[nj] * double(TmpSign);
				if (nj != 0) {
					ArcLen[nj] = ArcLen[nj - 1] + norm(v1);
				}
				else{
					StartArcLength = ArcLen[N-1] + norm(v1);
				}
				// 				if (FirstLoop){
				// 					LastSign = TmpSign;
				// 					FirstLoop = false;
				// 				}
				// 				else if (TmpSign != LastSign){
				// 					ZGCPointsLocal.push_back((Path.XYZAt(ni) + Path.XYZAt(nj)) * 0.5);
				// 					LastSign = TmpSign;
				// 				}
			}

			/* this is the data to be fitted */

// 			for (i = 0; i < n; ++i)
// 			{
// 				int ii = i - n_buffer;
// 				if (ii < 0) {
// 					auto pt = Path.XYZAt(N + ii);
// 					gsl_vector_set(x, i, pt[0]);
// 					gsl_vector_set(y, i, pt[1]);
// 					gsl_vector_set(z, i, pt[2]);
// 					gsl_vector_set(arclen, i, ArcLen[N + ii] - StartArcLength);
// 					gsl_vector_set(curvature, i, CSign[N + ii]);
// 				}
// 				else if (ii >= N){
// 					auto pt = Path.XYZAt(ii % N);
// 					gsl_vector_set(x, i, pt[0]);
// 					gsl_vector_set(y, i, pt[1]);
// 					gsl_vector_set(z, i, pt[2]);
// 					gsl_vector_set(arclen, i, ArcLen[ii % N] + StartArcLength);
// 					gsl_vector_set(curvature, i, CSign[ii % N]);
// 				}
// 				else {
// 					auto pt = Path.XYZAt(ii);
// 					gsl_vector_set(x, i, pt[0]);
// 					gsl_vector_set(y, i, pt[1]);
// 					gsl_vector_set(z, i, pt[2]);
// 					gsl_vector_set(arclen, i, ArcLen[ii]);
// 					gsl_vector_set(curvature, i, CSign[ii]);
// 
// 					double sigma = CSign[ii] * 0.1;
// 					gsl_vector_set(w, i, 1.0 / (sigma * sigma));
// 					sigma = pt[0] * 0.1;
// 					gsl_vector_set(wx, i, 1.0 / (sigma * sigma));
// 					sigma = pt[1] * 0.1;
// 					gsl_vector_set(wy, i, 1.0 / (sigma * sigma));
// 					sigma = pt[2] * 0.1;
// 					gsl_vector_set(wz, i, 1.0 / (sigma * sigma));
// 				}
// 				
// 				arclencheck[i] = gsl_vector_get(arclen, i);
// 			}
// 
// 			/* use uniform breakpoints on [0, 15] */
// 			gsl_bspline_knots_uniform(ArcLen[N - n_buffer] - StartArcLength, StartArcLength + ArcLen[n_buffer], bw);
// 			gsl_bspline_knots_uniform(ArcLen[N - n_buffer] - StartArcLength, StartArcLength + ArcLen[n_buffer], xw);
// 			gsl_bspline_knots_uniform(ArcLen[N - n_buffer] - StartArcLength, StartArcLength + ArcLen[n_buffer], yw);
// 			gsl_bspline_knots_uniform(ArcLen[N - n_buffer] - StartArcLength, StartArcLength + ArcLen[n_buffer], zw);
// 
// 			/* construct the fit matrix X */
// 			for (i = 0; i < n; ++i)
// 			{
// 				double xi = gsl_vector_get(arclen, i);
// 
// 				/* compute B_j(xi) for all j */
// 				gsl_bspline_eval(xi, B, bw);
// 
// 				/* fill in row i of X */
// 				for (j = 0; j < ncoeffs; ++j)
// 				{
// 					double Bj = gsl_vector_get(B, j);
// 					gsl_matrix_set(X, i, j, Bj);
// 				}
// 			}
// 
// 			/* do the fit */
// 			gsl_multifit_wlinear(X, w, curvature, c, cov, &chisq, mw);
// 			gsl_multifit_wlinear(X, wx, x, xc, xcov, &chisq, mw);
// 			gsl_multifit_wlinear(X, wy, y, yc, ycov, &chisq, mw);
// 			gsl_multifit_wlinear(X, wz, z, zc, zcov, &chisq, mw);
// 
// 			dof = n - ncoeffs;
// 			tss = gsl_stats_wtss(w->data, 1, curvature->data, 1, curvature->size);
// 			Rsq = 1.0 - chisq / tss;
// 
// 
// 			// Now generate values of first and second derivative of curvature
// 			// Where first deriv goes to 0 we have max/min, and where second deriv goes to zero we have inflection points
// 
// 			int CheckNumPoints = N * 10;
// 			vector<double> Csmooth(CheckNumPoints), Csmoothcheck(CheckNumPoints), Csmooth1(CheckNumPoints), Cprime(CheckNumPoints), Cdoubleprime(CheckNumPoints);
// 			double xStep = (arclencheck.back() - arclencheck.front() - 0.5) / double(CheckNumPoints - 1);
// 			gsl_matrix *dB;
// 			dB = gsl_matrix_alloc(nbreak + k - 2, 3);
// 			gsl_bspline_deriv_workspace *dw;
// 			dw = gsl_bspline_deriv_alloc(k);
// 			vector<vec3> outputs(CheckNumPoints);
// 			for (i = 0; i < CheckNumPoints; ++i){
// 				double xi = arclencheck.front() + double(i) * xStep, yi, yerr;
// 				// get smoothed xyz values
// 				gsl_bspline_eval(xi, B, xw);
// 				gsl_multifit_linear_est(B, xc, xcov, &outputs[i][0], &yerr);
// 				gsl_bspline_eval(xi, B, yw);
// 				gsl_multifit_linear_est(B, yc, ycov, &outputs[i][1], &yerr);
// 				gsl_bspline_eval(xi, B, zw);
// 				gsl_multifit_linear_est(B, zc, zcov, &outputs[i][2], &yerr);
// 
// 				for (int dir = 0; dir < 3; ++dir){
// 					if (isnan(outputs[i][dir]) || abs(outputs[i][dir]) > 1e10){
// 						outputs[i][dir] = 0.0;
// 					}
// 				}
// 				
// 				// smoothed curvature
// 				gsl_bspline_eval(xi, B, bw);
// 				gsl_multifit_linear_est(B, c, cov, &Csmooth[i], &yerr);
// 
// 				gsl_bspline_deriv_eval(xi, 2, dB, bw, dw);
// 				// smoothed curvature again, but from the deriv function as the 0-th derivative (sanity check)
// 				for (j = 0; j < nbreak + k - 2; ++j) {
// 					gsl_vector_set(B, j, gsl_matrix_get(dB, j, 0));
// 				}
// 				gsl_multifit_linear_est(B, c, cov, &Csmoothcheck[i], &yerr);
// 
// 				// smoothed first deriv of curvature
// 				for (j = 0; j < nbreak + k - 2; ++j) {
// 					gsl_vector_set(B, j, gsl_matrix_get(dB, j, 1));
// 				}
// 				gsl_multifit_linear_est(B, c, cov, &Cprime[i], &yerr);
// 
// 				// smoothed second deriv
// 				for (j = 0; j < nbreak + k - 2; ++j) {
// 					gsl_vector_set(B, j, gsl_matrix_get(dB, j, 2));
// 				}
// 				gsl_multifit_linear_est(B, c, cov, &Cdoubleprime[i], &yerr);
// 			}
// 			gsl_bspline_deriv_free(dw);
// 			gsl_matrix_free(dB);
// 
// 			int newzonenum = SaveVec3VecAsScatterZone(outputs, "smoothed curve");
// 			SetZoneStyle({ newzonenum }, ZoneStyle_Path, Green_C, 0.5);

			// Now find in

// 			fprintf(stderr, "chisq/dof = %e, Rsq = %f\n",
// 				chisq / dof, Rsq);

// 			printf("\n\n");

			/* output the smoothed curve */
// 			{
// 				double xi, yi, yerr;
// 
// 				for (xi = 0.0; xi < 15.0; xi += 0.1)
// 				{
// 					gsl_bspline_eval(xi, B, bw);
// 					gsl_multifit_linear_est(B, c, cov, &yi, &yerr);
// 					printf("%f %f\n", xi, yi);
// 				}
// 			}
// 
// 			gsl_rng_free(r);
// 			gsl_bspline_free(bw);
// 			gsl_bspline_free(xw);
// 			gsl_bspline_free(yw);
// 			gsl_bspline_free(zw);
// 			gsl_vector_free(B);
// 			gsl_vector_free(arclen);
// 			gsl_vector_free(curvature);
// 			gsl_vector_free(x);
// 			gsl_vector_free(y);
// 			gsl_vector_free(z);
// 			gsl_matrix_free(X);
// 			gsl_vector_free(c);
// 			gsl_vector_free(w);
// 			gsl_matrix_free(cov);
// 			gsl_vector_free(xc);
// 			gsl_vector_free(wx);
// 			gsl_matrix_free(xcov);
// 			gsl_vector_free(yc);
// 			gsl_vector_free(wy);
// 			gsl_matrix_free(ycov);
// 			gsl_vector_free(zc);
// 			gsl_vector_free(wz);
// 			gsl_matrix_free(zcov);
// 			gsl_multifit_linear_free(mw);

			vector<double> AvgLocalCurvature(N, 0.0);
			int avgCheckPts = 5;
			for (int ni = 0; ni < N; ++ni) {
				for (int i = 0; i < avgCheckPts * 2 + 1; ++i){
					int ii = ni - avgCheckPts + i;
					if (ii < 0) {
						ii = N + ii;
					}
					else if (ii >= N){
						ii = ii % N;
					}
					AvgLocalCurvature[ni] += CSign[ii];
				}
				AvgLocalCurvature[ni] /= 2 * avgCheckPts + 1;
			}

			// Now get min/max curvature points
			vector<int> ZGCPointInds;
			for (int ni = 0; ni < N; ++ni) {
				bool IsInflect = ABS(AvgLocalCurvature[ni]) < 0.1;

				for (int i = 1; i < 2 && IsInflect; ++i) {
					int njp = (ni + i) % N;
					if (njp == 0) njp++;

					int njm = ni - i;
					if (njm < 0) {
						njm = N + njm;
					}

// 					IsInflect = (CSign[njm] > 0.0 && CSign[njp] < 0.0) || (CSign[njm] < 0.0 && CSign[njp] > 0.0);
					IsInflect = (ABS(AvgLocalCurvature[njm]) > ABS(AvgLocalCurvature[ni]) && ABS(AvgLocalCurvature[njp]) > ABS(AvgLocalCurvature[ni]));
				}
				if (IsInflect) {
					ZGCPointsLocal.push_back(Path.XYZAt(ni));
					ZGCPointInds.push_back(ni);
				}
			}

// 			// combine nearby inflection points separated by flat regions
// 			for (auto const & i : ZGCPointInds) {
// 				ZGCPointsLocal.push_back(Path[i]);
// 			}
			for (int i = 0; i < ZGCPointInds.size(); ++i){
				int j = (i + 1) % ZGCPointInds.size(),
					ni = ZGCPointInds[i],
					nj = ZGCPointInds[j];
				if (j == i){
					break;
				}
				double TmpDist = Distance(Path[ni], Path[nj]);
				if (TmpDist < 0.5){
					double curvesum = 0.0, tmplength = Distance(Path[MIN(ni, nj)], Path[(MIN(ni, nj) + 1) % N]);
					if (abs(ni - nj) == 1){
						ZGCPointsLocal[i] == (Path[ni] + Path[nj]) * 0.5;
						ZGCPointsLocal.erase(ZGCPointsLocal.begin() + j);
						ZGCPointInds.erase(ZGCPointInds.begin() + j);
						i--;
					}
					else {
						for (int ii = MIN(ni, nj) + 1; ii < MAX(ni, nj); ++ii) {
							curvesum += CSign[ii];
							tmplength += Distance(Path[ii], Path[(ii + 1) % N]);
						}
						if (curvesum < 0.1) {
							tmplength *= 0.5;
							int ii = ni;
							double tmplength2 = 0.0;
							while (tmplength2 < tmplength) {
								tmplength2 += Distance(Path[ii], Path[(ii + 1) % N]);
								ii++;
							}
							ZGCPointsLocal[i] = (Path[ii] + Path[ii-1]) * 0.5;
							ZGCPointInds[i] = ii;
							ZGCPointsLocal.erase(ZGCPointsLocal.begin() + j);
							ZGCPointInds.erase(ZGCPointInds.begin() + j);
							i--;
						}
					}
				}
			}

			if (ZGCPointsLocal.size() % 2 == 0) {
				ZGCPoints.insert(ZGCPoints.end(), ZGCPointsLocal.begin(), ZGCPointsLocal.end());
			}
			int NumCheckPts = N / 4;
			vector<vec3> maxLocal, minLocal;
			double minmaxcheckscale = 1.;
			for (int ni = 0; ni < N; ++ni){

				bool IsMax = true;
				for (int i = 1; i < NumCheckPts && IsMax; ++i){
					int njp = (ni + i) % N;

					int njm = ni - i;
					if (njm < 0){
						njm = N + njm;
					}

					IsMax = (AvgLocalCurvature[njm] < AvgLocalCurvature[ni] * minmaxcheckscale && AvgLocalCurvature[ni] * minmaxcheckscale > AvgLocalCurvature[njp]);
// 					IsMax = (C[njm] < C[ni] * minmaxcheckscale && C[ni] * minmaxcheckscale > C[njp]);
				}
				if (IsMax){
					maxLocal.push_back(Path.XYZAt(ni));
				}

				bool IsMin = true;
				for (int i = 1; i < NumCheckPts && IsMin; ++i) {
					int njp = (ni + i) % N;

					int njm = ni - i;
					if (njm < 0) {
						njm = N + njm;
					}

					IsMin = (AvgLocalCurvature[njm] * minmaxcheckscale > AvgLocalCurvature[ni] && AvgLocalCurvature[ni] < AvgLocalCurvature[njp] * minmaxcheckscale);
// 					IsMin = (C[njm] * minmaxcheckscale > C[ni] && C[ni] < C[njp] * minmaxcheckscale);
				}
				if (IsMin) {
					minLocal.push_back(Path.XYZAt(ni));
				}
			}

			if (!maxLocal.empty() && !minLocal.empty()){
				MaxCurvaturePoints.insert(MaxCurvaturePoints.end(), maxLocal.begin(), maxLocal.end());
				MinCurvaturePoints.insert(MinCurvaturePoints.end(), minLocal.begin(), minLocal.end());
			}

		}
	}

	if (!ZGCPoints.empty()){
		SaveVec3VecAsScatterZone(ZGCPoints, "Curvature inflection points", Purple_C);
	}
	if (!MaxCurvaturePoints.empty()) {
		SaveVec3VecAsScatterZone(MaxCurvaturePoints, "Maximum curvature points", Red_C);
	}
	if (!MinCurvaturePoints.empty()) {
		SaveVec3VecAsScatterZone(MinCurvaturePoints, "Minimum curvature points", Blue_C);
	}

	if (UseCurrentContours && !ZoneNums.empty()){
		Set DeleteZones;
		for(auto const & z : ZoneNums){
			DeleteZones += z;
		}
		TecUtilZoneDelete(DeleteZones.getRef());
	}
	TecUtilRedrawAll(TRUE);
}

void FindContourCurvaturePointsReturnUserInfo(bool const GuiSuccess,
	vector<GuiField_c> const & Fields,
	vector<GuiField_c> const PassthroughFields)
{
	int fNum = 0;
	bool UseExistingContours = Fields[fNum++].GetReturnBool();
	vector<int> ZoneNums = Fields[fNum++].GetReturnIntVec();
	vector<int> XYZVarNums(3);
	for (int i = 0; i < 3; ++i) XYZVarNums[i] = i + Fields[fNum].GetReturnInt();

	TecUtilLockStart(AddOnID);

	FindContourCurvaturePoints(ZoneNums, XYZVarNums, UseExistingContours);

	TecUtilLockFinish(AddOnID);
}

void FindContourCurvaturePointsGetUserInfo() {
	vector<GuiField_c> Fields = {
		GuiField_c(Gui_Toggle, "Use current set of contours"),
		GuiField_c(Gui_ZoneSelectMulti,  "Or using extracted contours zones:"),
		GuiField_c(Gui_VarSelect, "X", "X")
	};

	CSMGui("Find contour curvature max/min/inflection points", Fields, FindContourCurvaturePointsReturnUserInfo, AddOnID);
}