
#include <vector>
#include <string>
#include <sstream>
#include <algorithm>
#include <cmath>
#include <omp.h>

//#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_eigen.h>

#include "TECADDON.h"
#include "CSM_DATA_TYPES.h"
#include "CSM_CRIT_POINTS.h"
#include "CSM_DATA_SET_INFO.h"
#include "CSM_VOL_EXTENT_INDEX_WEIGHTS.h"

#include "CSM_CALC_VARS.h"

#include <armadillo>
using namespace arma;

using std::vector;
using std::string;
using std::to_string;
using std::stringstream;


/*
*	These indices are used for second-order derivative calculations.
*/
static vector<vector<int> > const ValInds = { { 0, 1, 3, 4 }, { 1, 3 }, { 2, 3, 4 }, { 0, 1, 2 } };


void CalcGradGradMagForDataset(Boolean_t IsPeriodic, AddOn_pa const & AddOnID){

	vector<string> GradXYZMagStr = CSMVarName.DensGradVec;
	GradXYZMagStr.push_back(CSMVarName.DensGradMag);

	vector<EntIndex_t> GradXYZMagVarNum(GradXYZMagStr.size(), -1);

	TecUtilLockStart(AddOnID);

	EntIndex_t RhoVarNum = VarNumByNameList(vector<string>({
		CSMVarName.Dens,
		"Rho",
		"rho"
	}));

	for (int i = 0; i < GradXYZMagStr.size(); ++i)
		GradXYZMagVarNum[i] = VarNumByName(GradXYZMagStr[i]);

	Boolean_t HasGrad = (GradXYZMagVarNum[0] > 0 && GradXYZMagVarNum[1] > 0 && GradXYZMagVarNum[2] > 0);
	Boolean_t HasGradMag = GradXYZMagVarNum[3] > 0;

	if (HasGrad && HasGradMag){
		TecUtilLockFinish(AddOnID);
		return;
	}

	if (RhoVarNum <= 0 && (!HasGrad || !HasGradMag)){
		TecUtilDialogErrMsg("Couldn't find electron density variable, so did not compute gradient or gradient magnitude!");
		TecUtilLockFinish(AddOnID);
		return;
	}

	EntIndex_t ZoneNum = ZoneNumByName(string("Full Volume"));
	vector<LgIndex_t> MaxIJK(3, -1);
	if (ZoneNum > 0){
		TecUtilZoneGetIJK(ZoneNum, &MaxIJK[0], &MaxIJK[1], &MaxIJK[2]);
	}
	if (ZoneNum <= 0 || MaxIJK[0] <= 0 || MaxIJK[1] <= 0 || MaxIJK[2] <= 0){
		TecUtilDialogErrMsg("Failed to get zone information. Didn't calculation gradient or gradient magnitude.");
		TecUtilLockFinish(AddOnID);
		return;
	}

	vec3 DelXYZ;
	vector<EntIndex_t> XYZVarNums(3, -1);
	TecUtilAxisGetVarAssignments(&XYZVarNums[0], &XYZVarNums[1], &XYZVarNums[2]);
	if (XYZVarNums[0] <= 0 || XYZVarNums[1] <= 0 || XYZVarNums[2] <= 0){
		vector<vector<string> > TmpStrs = {
			{ "X", "Y", "Z" }, { "I", "J", "K" }
		};
		Boolean_t VarsFound = FALSE;
		for (int i = 0; i < 2 && !VarsFound; ++i){
			for (int j = 0; j < 3; ++j)
				XYZVarNums[j] = TecUtilVarGetNumByName(TmpStrs[i][j].c_str());
			VarsFound = (XYZVarNums[0] > 0 && XYZVarNums[1] > 0 && XYZVarNums[2] > 0);
		}
		if (!VarsFound){
			TecUtilDialogErrMsg("Failed to get x, y, z variable numbers. Didn't calculate gradient or gradient magnitude");
			TecUtilLockFinish(AddOnID);
			return;
		}
	}

	/*
	*	Using this PointSpacingV123 for the distance in the derivative approximations
	*	assumes regular spacing in the data!!!
	*/
	// 	for (int i = 0; i < 3; ++i){
	// 		double MinMax[2];
	// 		TecUtilVarGetMinMax(XYZVarNums[i], &MinMax[0], &MinMax[1]);
	// 		PointSpacingV123[i] = (MinMax[1] - MinMax[0]) / (MaxIJK[i] - 1);
	// 	}
	DelXYZ = GetDelXYZ_Ordered3DZone(XYZVarNums, ZoneNum);

	vector<FieldDataPointer_c> GradXYZMagRawPtr(GradXYZMagVarNum.size());

	TecUtilDataLoadBegin();

	FieldDataPointer_c RhoRawPtr;
	RhoRawPtr.InitializeReadPtr(ZoneNum, RhoVarNum);

	EntIndex_t NumZones = TecUtilDataSetGetNumZones();
	vector<FieldDataType_e> DataType;
	DataType.resize(NumZones, FieldDataType_Float);

	vector<ValueLocation_e> DataLoc(NumZones, ValueLocation_Nodal);

	if (!HasGrad){
		for (int i = 0; i < 3; ++i){
			if (GradXYZMagVarNum[i] <= 0){
				ArgList_pa Args = TecUtilArgListAlloc();
				TecUtilArgListAppendString(Args, SV_NAME, GradXYZMagStr[i].c_str());
				TecUtilArgListAppendArray(Args, SV_VARDATATYPE, DataType.data());
				TecUtilArgListAppendArray(Args, SV_VALUELOCATION, DataLoc.data());
				TecUtilArgListAppendInt(Args, SV_SHAREVARWITHALLZONES, FALSE);

				if (TecUtilDataSetAddVarX(Args)){
					Set_pa NewVar = TecUtilSetAlloc(FALSE);
					GradXYZMagVarNum[i] = TecUtilDataSetGetNumVars();
					TecUtilSetAddMember(NewVar, GradXYZMagVarNum[i], FALSE);
					TecUtilStateChanged(StateChange_VarsAdded, reinterpret_cast<ArbParam_t>(NewVar));
					TecUtilSetDealloc(&NewVar);
				}
				TecUtilArgListDealloc(&Args);
			}

			GradXYZMagRawPtr[i].InitializeWritePtr(ZoneNum, GradXYZMagVarNum[i]);
		}

		if (!CalcGradForRegularVar(MaxIJK, DelXYZ, IsPeriodic, RhoRawPtr, GradXYZMagRawPtr, "gradient vector", AddOnID)){
			Set_pa Vars = TecUtilSetAlloc(FALSE);
			for (int i = 0; i < 3; ++i)
				TecUtilSetAddMember(Vars, GradXYZMagVarNum[i], FALSE);
			TecUtilDataSetDeleteVar(Vars);
			TecUtilSetDealloc(&Vars);
			TecUtilLockFinish(AddOnID);
			return;
		}
	}

	TecUtilDataLoadEnd();

	TecUtilDataLoadBegin();

	if (!HasGradMag){
		for (int i = 0; i < 4; ++i){
			if (GradXYZMagVarNum[i] <= 0){
				ArgList_pa Args = TecUtilArgListAlloc();
				TecUtilArgListAppendString(Args, SV_NAME, GradXYZMagStr[i].c_str());
				TecUtilArgListAppendArray(Args, SV_VARDATATYPE, DataType.data());
				TecUtilArgListAppendArray(Args, SV_VALUELOCATION, DataLoc.data());
				TecUtilArgListAppendInt(Args, SV_SHAREVARWITHALLZONES, FALSE);

				if (TecUtilDataSetAddVarX(Args)){
					Set_pa NewVar = TecUtilSetAlloc(FALSE);
					GradXYZMagVarNum[i] = TecUtilDataSetGetNumVars();
					TecUtilSetAddMember(NewVar, GradXYZMagVarNum[i], FALSE);
					TecUtilStateChanged(StateChange_VarsAdded, reinterpret_cast<ArbParam_t>(NewVar));
					TecUtilSetDealloc(&NewVar);
				}
				TecUtilArgListDealloc(&Args);
			}

			if (i < 3){
				GradXYZMagRawPtr[i].InitializeReadPtr(ZoneNum, GradXYZMagVarNum[i]);
			}
			else{
				GradXYZMagRawPtr[i].InitializeWritePtr(ZoneNum, GradXYZMagVarNum[i]);
			}
		}

		if (!CalcMagForRegularVectorVar(MaxIJK, GradXYZMagRawPtr, GradXYZMagRawPtr[3], "gradient magnitude", AddOnID)){
			Set_pa Vars = TecUtilSetAlloc(FALSE);
			TecUtilSetAddMember(Vars, GradXYZMagVarNum[3], FALSE);
			TecUtilDataSetDeleteVar(Vars);
			TecUtilSetDealloc(&Vars);
			TecUtilLockFinish(AddOnID);
			return;
		}
	}

	TecUtilDataLoadEnd();
	TecUtilLockFinish(AddOnID);
}

Boolean_t CalcGradForRegularVar(vector<int> const & IJKMax,
	vec3 const & DelXYZ,
	Boolean_t IsPeriodic,
	FieldDataPointer_c const & VarReadPtr,
	vector<FieldDataPointer_c> const & VarWritePtrs,
	string const & VarName,
	AddOn_pa const & AddOnID)
{
	Boolean_t TaskQuit = FALSE;

	TaskQuit = !(IJKMax.size() == 3
		&& VarWritePtrs.size() >= 3);
	for (int i = 0; i < 3 && !TaskQuit; ++i){
		TaskQuit = !VarWritePtrs[i].IsReady();
	}

	SYSTEM_INFO sysinfo;
	GetSystemInfo(&sysinfo);

	int numCPU = sysinfo.dwNumberOfProcessors;
	omp_set_num_threads(numCPU);

	vec3 OneOverDelXYZx2 = 1./(DelXYZ * 2.), OneOverDelXYZx12 = 1./(DelXYZ * 12.);
	vector<vector<double> > Vals(numCPU, vector<double>(5));
	vector<vector<int> > DirInd(numCPU, vector<int>(5));

	int StatusKMax = IJKMax[2] / numCPU;

	string StatusStr = string("Calculating ") + VarName;

	StatusLaunch(StatusStr.c_str(), AddOnID, TRUE);
	vector<vec3> OutValues(numCPU);

	if (!TaskQuit){
#ifndef DEBUG
#pragma omp parallel
#endif
		for (LgIndex_t kk = 1; kk <= IJKMax[2]; ++kk){
			int ThreadNum = omp_get_thread_num();
			if (ThreadNum == 0 && !StatusUpdate(kk-1, StatusKMax, StatusStr, AddOnID)){
				TaskQuit = TRUE;
#pragma omp flush (TaskQuit)
			}
#pragma omp flush (TaskQuit)
			for (LgIndex_t jj = 1; jj <= IJKMax[1] && !TaskQuit; ++jj){
				for (LgIndex_t ii = 1; ii <= IJKMax[0]; ++ii){
					int Index = IndexFromIJK(ii, jj, kk, IJKMax[0], IJKMax[1], IJKMax[2], IsPeriodic) - 1;
					CalcGradForNode(ii, jj, kk,
						OneOverDelXYZx2, OneOverDelXYZx12, Vals[ThreadNum], DirInd[ThreadNum], 0,
						IJKMax, IsPeriodic, VarReadPtr, OutValues[ThreadNum]);
					for (int i = 0; i < 3; ++i)
						VarWritePtrs[i].Write(Index, OutValues[ThreadNum][i]);
				}
			}
		}
	}

	StatusDrop(AddOnID);

	return !TaskQuit;
}

Boolean_t CalcMagForRegularVectorVar(vector<int> const & IJKMax,
	vector<FieldDataPointer_c> const & VarReadPtrs,
	FieldDataPointer_c const & VarWritePtr,
	string const & VarName,
	AddOn_pa const & AddOnID)
{
	Boolean_t TaskQuit = FALSE;

	TaskQuit = !(IJKMax.size() == 3
		&& VarReadPtrs.size() >= 3);

	SYSTEM_INFO sysinfo;
	GetSystemInfo(&sysinfo);

	int numCPU = sysinfo.dwNumberOfProcessors;
	omp_set_num_threads(numCPU);

	int StatusKMax = IJKMax[2] / numCPU;

	string StatusStr = string("Calculating ") + VarName;

	StatusLaunch(StatusStr.c_str(), AddOnID, TRUE);

	if (!TaskQuit){
#pragma omp parallel for
		for (LgIndex_t kk = 1; kk <= IJKMax[2]; ++kk){
			if (omp_get_thread_num() == 0 && !StatusUpdate(kk-1, StatusKMax, StatusStr, AddOnID)){
				TaskQuit = TRUE;
#pragma omp flush (TaskQuit)
			};
#pragma omp flush (TaskQuit)
			vec3 Vals;
			for (LgIndex_t jj = 1; jj <= IJKMax[1] && !TaskQuit; ++jj){
				for (LgIndex_t ii = 1; ii <= IJKMax[0]; ++ii){

					LgIndex_t Index = IndexFromIJK(ii, jj, kk, IJKMax[0], IJKMax[1], IJKMax[2], FALSE) - 1;

					for (int i = 0; i < 3; ++i)
						Vals[i] = VarReadPtrs[i][Index];

					VarWritePtr.Write(Index, norm(Vals));
				}
			}
		}
	}

	StatusDrop(AddOnID);

	return !TaskQuit;
}

void CalcGradForNode(int ii,
	int jj,
	int kk,
	vec3 const & OneOverDelXYZx2,
	vec3 const & OneOverDelXYZx12,
	vector<double> & Vals,
	vector<int> & DirInd,
	int StartDir,
	vector<int> const & IJKMax,
	Boolean_t IsPeriodic,
	FieldDataPointer_c const & VarReadPtr,
	vec3 & OutValues)
{
	for (int Dir = 0; Dir < 3; ++Dir){
		/*
		*	Use a combination of four potential methods for approximating
		*	derivative:
		*		0. High accuracy centered divided difference
		*		1. Centered divided difference
		*		2. High accuracy forward divided difference
		*		3. High accuracy backward divided difference
		*	Which is used depends on how close the point is to the system
		*	boundary (only for non-periodic systems, periodic always get
		*	high accuracy centered).
		*	If not periodic:
		*		Points on the boundary get the high accuracy forward or
		*			backward method
		*		Points 1 away from the boundary get the regular centered
		*		Points 2 or more away from the boundary get the centered
		*			high accuracy method
		*/

		int Method = 0;
		/*
		*	The 5 elements of DirInd of I,J,K are for the
		*	minus 2, minus 1, plus 0, plus 1, plus 2
		*	respectively.
		*/

		switch (Dir){
			case 0:
				for (int i = 0; i < 5; ++i)
					DirInd[i] = ii + i - 2;
				break;
			case 1:
				for (int i = 0; i < 5; ++i)
					DirInd[i] = jj + i - 2;
				break;
			case 2:
				for (int i = 0; i < 5; ++i)
					DirInd[i] = kk + i - 2;
				break;
		}

		if (!IsPeriodic && (DirInd[2] <= 2 || DirInd[2] >= IJKMax[Dir] - 1)){
			if (DirInd[2] == 1){
				Method = 2;
			}
			else if (DirInd[2] == IJKMax[Dir]){
				Method = 3;
			}
			else{ // point is 1 away from boundary
				Method = 1;
			}
		}

		/*
		*	Get the values for the current of the variable at the points found.
		*
		*	Elements of Vals[] correspond to the
		*	minus 2, minus 1, plus 0, plus 1, plus 2
		*	points in the current direction.
		*/

		switch (Dir){
			case 0:
				for (int i : ValInds[Method])
					Vals[i] = VarReadPtr[IndexFromIJK(DirInd[i], jj, kk, IJKMax[0], IJKMax[1], IJKMax[2], IsPeriodic) - 1];
				break;
			case 1:
				for (int i : ValInds[Method])
					Vals[i] = VarReadPtr[IndexFromIJK(ii, DirInd[i], kk, IJKMax[0], IJKMax[1], IJKMax[2], IsPeriodic) - 1];
				break;
			case 2:
				for (int i : ValInds[Method])
					Vals[i] = VarReadPtr[IndexFromIJK(ii, jj, DirInd[i], IJKMax[0], IJKMax[1], IJKMax[2], IsPeriodic) - 1];
				break;
		}

		if (OneOverDelXYZx12[Dir] != 0){
			switch (Method){
				case 0: // High accuracy centered divided difference
					OutValues[Dir] = (-Vals[4] + 8.0 * (Vals[3] - Vals[1]) + Vals[0]) * OneOverDelXYZx12[Dir];
					break;
				case 1: // Centered divided difference
					OutValues[Dir] = (Vals[3] - Vals[1]) * OneOverDelXYZx2[Dir];
					break;
				case 2: // High accuracy forward divided difference
					OutValues[Dir] = (-Vals[4] + 4.0 * Vals[3] - 3.0 * Vals[2]) * OneOverDelXYZx2[Dir];
					break;
				case 3: // High accuracy backward divided difference
					OutValues[Dir] = (3.0 * Vals[2] - 4.0 * Vals[1] + Vals[0]) * OneOverDelXYZx2[Dir];
					break;
			}
		}
		else OutValues[Dir] = 0.;

		/*
		*	Here's the previous implementation, which uses a
		*	combination of basic central/forward/backward
		*	divided difference methods instead of the higher
		*	accuracy versions above.
		*/
		// 						double DelMult = 2.0;
		// 						LgIndex_t PM1[3], PP1[3], Pt[3];
		// 						Pt[0] = PM1[0] = PP1[0] = ii;
		// 						Pt[1] = PM1[1] = PP1[1] = jj;
		// 						Pt[2] = PM1[2] = PP1[2] = kk;
		// 
		// 						PM1[Dir]--;
		// 						PP1[Dir]++;
		// 
		// 						if (!IsPeriodic){
		// 							Boolean_t BackForFiniteDiff = FALSE;
		// 							LgIndex_t IJK[3] = { ii, jj, kk };
		// 							for (int i = 0; i < 3; ++i){
		// 								if (IJK[i] == 1){
		// 									PM1[i] = Pt[i];
		// 									BackForFiniteDiff = TRUE;
		// 								}
		// 								else if (IJK[i] == IJKMax[i]){
		// 									PP1[i] = Pt[i];
		// 									BackForFiniteDiff = TRUE;
		// 								}
		// 							}
		// 							if (BackForFiniteDiff)
		// 								DelMult--;
		// 						}
		// 
		// 						LgIndex_t Point = IndexFromIJK(Pt[0], Pt[1], Pt[2], IJKMax[0], IJKMax[1], IJKMax[2], IsPeriodic) - 1;
		// 						LgIndex_t PtMinus = IndexFromIJK(PM1[0], PM1[1], PM1[2], IJKMax[0], IJKMax[1], IJKMax[2], IsPeriodic) - 1;
		// 						LgIndex_t  PtPlus = IndexFromIJK(PP1[0], PP1[1], PP1[2], IJKMax[0], IJKMax[1], IJKMax[2], IsPeriodic) - 1;
		// 
		// 						float PlusRho = VarReadPtr[PtPlus];
		// 						float MinusRho = VarReadPtr[PtMinus];
		// 						float DRho = (PlusRho - MinusRho) / (PointSpacingV123[Dir] * DelMult);
		// 						OutValues.at(Dir, Point) = DRho;
	}
}



void CalcGradForPoint(vec3 const & Point,
	vec3 const & DelXYZ,
	VolExtentIndexWeights_s & VolInfo,
	mat33 const & DirVects,
	int StartDir,
	Boolean_t IsPeriodic,
	vec & OutValues,
	FieldDataPointer_c const & VarReadPtr,
	GPType_e const & CalcType,
	void * Params)
{
	double Vals[5];
	vec3 Points[5];
	vec3 OneOverDelXYZx2 = 1./(DelXYZ * 2.), OneOverDelXYZx12 = 1./(DelXYZ * 12.);

	int NumDirs = sqrt(DirVects.size());

	for (int Dir = StartDir; Dir < NumDirs; ++Dir){
		/*
		*	Use a combination of four potential methods for approximating
		*	derivative:
		*		0. High accuracy centered divided difference
		*		1. Centered divided difference
		*		2. High accuracy forward divided difference
		*		3. High accuracy backward divided difference
		*	Which is used depends on how close the point is to the system
		*	boundary (only for non-periodic systems, periodic always get
		*	high accuracy centered).
		*	If not periodic:
		*		Points on the boundary get the high accuracy forward or
		*			backward method
		*		Points 1 away from the boundary get the regular centered
		*		Points 2 or more away from the boundary get the centered
		*			high accuracy method
		*/

		int Method = 0;
		/*
		*	The 5 elements of DirInd of I,J,K are for the
		*	minus 2, minus 1, plus 0, plus 1, plus 2
		*	respectively.
		*/

		for (int i = 0; i < 5; ++i)
			Points[i] = Point + DelXYZ % (DirVects.col(Dir) * static_cast<double>(i - 2));

		if (!IsPeriodic && (sum(Points[0] < VolInfo.MinXYZ) > 0 || sum(Points[4] > VolInfo.MaxXYZ) > 0)){
			if (sum(Points[0] < VolInfo.MinXYZ) > 0){
				if (sum(Points[1] < VolInfo.MinXYZ) > 0){ // High accuracy forward divided difference
					Method = 2;
				}
			}
			else if (sum(Points[4] > VolInfo.MaxXYZ) > 0){ // High accuracy backward divided difference
				if (sum(Points[3] > VolInfo.MaxXYZ) > 0){
					Method = 3;
				}
			}
			else{ // Centered divided difference
				Method = 1;
			}
		}
		else{ // High accuracy centered divided difference
			if (IsPeriodic){
				for (int i = 0; i < 5; ++i){
					for (int j = 0; j < 3; ++j){
						if (Points[i][j] < VolInfo.MinXYZ[j])
							Points[i][j] = VolInfo.MaxXYZ[j] - (VolInfo.MinXYZ[j] - Points[i][j]);
						else if (Points[i][j] > VolInfo.MaxXYZ[j])
							Points[i][j] = VolInfo.MinXYZ[j] + (Points[i][j] - VolInfo.MaxXYZ[j]);
					}
				}
			}
		}

		/*
		*	Get the values for the current direction of the variable at the points found.
		*
		*	Elements of Vals[] correspond to the
		*	minus 2, minus 1, plus 0, plus 1, plus 2
		*	points in the current direction.
		*/

		for (int i : ValInds[Method]){
			switch (CalcType){
				case GPType_NormalPlaneEberlyCP:
					Vals[i] = Eberly1RidgeFunction(Points[i], 0.0, TRUE, *reinterpret_cast<MultiRootParams_s*>(Params));
					break;
				case GPType_NEB:
					Vals[i] = NEBForceFunction(Points[i], *reinterpret_cast<MultiRootParams_s*>(Params));
				default:
					SetIndexAndWeightsForPoint(Points[i], VolInfo);
					Vals[i] = ValByCurrentIndexAndWeightsFromRawPtr(VolInfo, VarReadPtr);
					break;
			}
		}

		if (OneOverDelXYZx12[Dir]){
			switch (Method){
				case 0: // High accuracy centered divided difference
					OutValues[Dir] = (-Vals[4] + 8.0 * (Vals[3] - Vals[1]) + Vals[0]) * OneOverDelXYZx12[Dir];
					break;
				case 1: // Centered divided difference
					OutValues[Dir] = (Vals[3] - Vals[1]) * OneOverDelXYZx2[Dir];
					break;
				case 2: // High accuracy forward divided difference
					OutValues[Dir] = (-Vals[4] + 4.0 * Vals[3] - 3.0 * Vals[2]) * OneOverDelXYZx2[Dir];
					break;
				case 3: // High accuracy backward divided difference
					OutValues[Dir] = (3.0 * Vals[2] - 4.0 * Vals[1] + Vals[0]) * OneOverDelXYZx2[Dir];
					break;
			}
		}
		else OutValues[Dir] = 0.;
	}
}

void CalcHessForNode(int ii,
	int jj,
	int kk,
	vec3 const & OneOverDelXYZx2,
	vec3 const & OneOverDelXYZx12,
	vector<double> & Vals,
	vector<int> & DirInd,
	vector<int> const & IJKMax,
	Boolean_t IsPeriodic,
	FieldDataPointer_c const & ScalarReadPtr,
	vector<FieldDataPointer_c> const & GradReadPtrs,
	mat33 & Hess)
{
	vector<vec3> GradValues(3);
	vector<double> TmpVals(5);
	vector<int> TmpInd(5);

	Boolean_t HasGrad = TRUE;
	for (int i = 0; i < 3 && HasGrad; ++i)
		HasGrad = GradReadPtrs[i].IsReady();

	if (HasGrad){
		for (int i = 0; i < 3; ++i){
			CalcGradForNode(ii,
				jj,
				kk,
				OneOverDelXYZx2,
				OneOverDelXYZx12,
				Vals,
				DirInd,
				i,
				IJKMax,
				IsPeriodic,
				GradReadPtrs[i],
				GradValues[i]);
			for (int j = 0; j < 3; ++j){
				/*
				*	This if statement makes is so that the symmetrical elements of the
				*	hessian are simply copied from the top-right of the matrix rather
				*	than being recomputed.
				*/
				if (j >= i)
					Hess.at(i, j) = GradValues[i][j];
				else
					Hess.at(i, j) = GradValues[j][i];
			}
		}
	}
	else{
		vec3 Grad[5];
		for (int iDir = 0; iDir < 3; ++iDir){
			/*
			*	Use a combination of four potential methods for approximating
			*	derivative:
			*		0. High accuracy centered divided difference
			*		1. Centered divided difference
			*		2. High accuracy forward divided difference
			*		3. High accuracy backward divided difference
			*	Which is used depends on how close the point is to the system
			*	boundary (only for non-periodic systems, periodic always get
			*	high accuracy centered).
			*	If not periodic:
			*		Points on the boundary get the high accuracy forward or
			*			backward method
			*		Points 1 away from the boundary get the regular centered
			*		Points 2 or more away from the boundary get the centered
			*			high accuracy method
			*/

			int Method = 0;
			/*
			*	The 5 elements of DirInd of I,J,K are for the
			*	minus 2, minus 1, plus 0, plus 1, plus 2
			*	respectively.
			*/

			switch (iDir){
				case 0:
					for (int i = 0; i < 5; ++i)
						DirInd[i] = ii + i - 2;
					break;
				case 1:
					for (int i = 0; i < 5; ++i)
						DirInd[i] = jj + i - 2;
					break;
				case 2:
					for (int i = 0; i < 5; ++i)
						DirInd[i] = kk + i - 2;
					break;
			}

			if (!IsPeriodic && (DirInd[2] <= 2 || DirInd[2] >= IJKMax[iDir] - 1)){
				if (DirInd[2] == 1){
					Method = 2;
				}
				else if (DirInd[2] == IJKMax[iDir]){
					Method = 3;
				}
				else{ // point is 1 away from boundary
					Method = 1;
				}
			}

			/*
			*	Get the values for the current of the variable at the points found.
			*
			*	Elements of Vals[] correspond to the
			*	minus 2, minus 1, plus 0, plus 1, plus 2
			*	points in the current direction.
			*/

			switch (iDir){
				case 0:
					for (int i : ValInds[Method])
						CalcGradForNode(DirInd[i], jj, kk, OneOverDelXYZx2, OneOverDelXYZx12, TmpVals, TmpInd, iDir, IJKMax, IsPeriodic, ScalarReadPtr, Grad[i]);
					break;
				case 1:
					for (int i : ValInds[Method])
						CalcGradForNode(ii, DirInd[i], kk, OneOverDelXYZx2, OneOverDelXYZx12, TmpVals, TmpInd, iDir, IJKMax, IsPeriodic, ScalarReadPtr, Grad[i]);
					break;
				case 2:
					for (int i : ValInds[Method])
						CalcGradForNode(ii, jj, DirInd[i], OneOverDelXYZx2, OneOverDelXYZx12, TmpVals, TmpInd, iDir, IJKMax, IsPeriodic, ScalarReadPtr, Grad[i]);
					break;
			}

			for (int jDir = 0; jDir < 3; ++jDir){
				if (jDir >= iDir){
					if (OneOverDelXYZx12[iDir] != 0){
						switch (Method){
							case 0: // High accuracy centered divided difference
								Hess.at(iDir, jDir) = (-Grad[4][jDir] + 8.0 * (Grad[3][jDir] - Grad[1][jDir]) + Grad[0][jDir]) * OneOverDelXYZx12[iDir];
								break;
							case 1: // Centered divided difference
								Hess.at(iDir, jDir) = (Grad[3][jDir] - Grad[1][jDir]) * OneOverDelXYZx2[iDir];
								break;
							case 2: // High accuracy forward divided difference
								Hess.at(iDir, jDir) = (-Grad[4][jDir] + 4.0 * Grad[3][jDir] - 3.0 * Grad[2][jDir]) * OneOverDelXYZx2[iDir];
								break;
							case 3: // High accuracy backward divided difference
								Hess.at(iDir, jDir) = (3.0 * Grad[2][jDir] - 4.0 * Grad[1][jDir] + Grad[0][jDir]) * OneOverDelXYZx2[iDir];
								break;
						}
					}
					else Hess.at(iDir, jDir) = 0.;
				}
				else{
					Hess.at(iDir, jDir) = Hess.at(jDir, iDir);
				}
			}
		}
	}
}

void CalcHessForPoint(vec3 const & Point,
	vec3 const & DelXYZ,
	VolExtentIndexWeights_s & VolInfo,
	mat33 const & DirVects,
	Boolean_t IsPeriodic,
	mat & OutValues,
	FieldDataPointer_c const & VarReadPtr,
	GPType_e const & CalcType,
	void * Params)
{
	vec3 Points[5], Grad[5];
	vec3 OneOverDelXYZx2 = 1./(DelXYZ * 2.), OneOverDelXYZx12 = 1./(DelXYZ * 12.);

	int Rank = sqrt(DirVects.size());

	for (int iDir = 0; iDir < Rank; ++iDir){
		/*
		*	Use a combination of four potential methods for approximating
		*	derivative:
		*		1. High accuracy centered divided difference
		*		2. Centered divided difference
		*		3. High accuracy forward divided difference
		*		4. High accuracy backward divided difference
		*	Which is used depends on how close the point is to the system
		*	boundary (only for non-periodic systems, periodic always get
		*	high accuracy centered).
		*	If not periodic:
		*		Points on the boundary get the high accuracy forward or
		*			backward method
		*		Points 1 away from the boundary get the regular centered
		*		Points 2 or more away from the boundary get the centered
		*			high accuracy method
		*/

		int Method = 0;
		/*
		*	The 5 elements of DirInd of I,J,K are for the
		*	minus 2, minus 1, plus 0, plus 1, plus 2
		*	respectively.
		*/

		for (int i = 0; i < 5; ++i)
			Points[i] = Point + DelXYZ % DirVects.col(iDir) * static_cast<double>(i - 2);

		if (!IsPeriodic && (sum(Points[0] < VolInfo.MinXYZ) > 0 || sum(Points[4] > VolInfo.MaxXYZ) > 0)){
			if (sum(Points[0] < VolInfo.MinXYZ) > 0){
				if (sum(Points[1] < VolInfo.MinXYZ) > 0){ // High accuracy forward divided difference
					Method = 2;
				}
			}
			else if (sum(Points[4] > VolInfo.MaxXYZ) > 0){ // High accuracy backward divided difference
				if (sum(Points[3] > VolInfo.MaxXYZ) > 0){
					Method = 3;
				}
			}
			else{ // Centered divided difference
				Method = 1;
			}
		}
		else{ // High accuracy centered divided difference
			if (IsPeriodic){
				for (int i = 0; i < 5; ++i){
					for (int j = 0; j < 3; ++j){
						if (Points[i][j] < VolInfo.MinXYZ[j])
							Points[i][j] = VolInfo.MaxXYZ[j] - (VolInfo.MinXYZ[j] - Points[i][j]);
						else if (Points[i][j] > VolInfo.MaxXYZ[j])
							Points[i][j] = VolInfo.MinXYZ[j] + (Points[i][j] - VolInfo.MaxXYZ[j]);
					}
				}
			}
		}

		/*
		*	Get the values for the current direction of the variable at the points found.
		*
		*	Elements of Grad[] correspond to the
		*	minus 2, minus 1, plus 0, plus 1, plus 2
		*	points in the current direction.
		*/

		for (int i : ValInds[Method]){
			CalcGradForPoint(Points[i], DelXYZ, VolInfo, DirVects, iDir, IsPeriodic, Grad[i], VarReadPtr, CalcType, Params);
		}

		for (int jDir = 0; jDir < Rank; ++jDir){
			/*
			*	This if statement makes is so that the symmetrical elements of the
			*	hessian are simply copied from the top-right of the matrix rather
			*	than being recomputed.
			*/
			if (jDir >= iDir){
				if (DelXYZ[iDir] != 0){
					switch (Method){
						case 0: // High accuracy centered divided difference
							OutValues.at(iDir, jDir) = (-Grad[4][jDir] + 8.0 * (Grad[3][jDir] - Grad[1][jDir]) + Grad[0][jDir]) * OneOverDelXYZx12[iDir];
							break;
						case 1: // Centered divided difference
							OutValues.at(iDir, jDir) = (Grad[3][jDir] - Grad[1][jDir]) * OneOverDelXYZx2[iDir];
							break;
						case 2: // High accuracy forward divided difference
							OutValues.at(iDir, jDir) = (-Grad[4][jDir] + 4.0 * Grad[3][jDir] - 3.0 * Grad[2][jDir]) * OneOverDelXYZx2[iDir];
							break;
						case 3: // High accuracy backward divided difference
							OutValues.at(iDir, jDir) = (3.0 * Grad[2][jDir] - 4.0 * Grad[1][jDir] + Grad[0][jDir]) * OneOverDelXYZx2[iDir];
							break;
					}
				}
				else OutValues.at(iDir, jDir) = 0.;
			}
			else{
				OutValues.at(iDir, jDir) = OutValues.at(jDir, iDir);
			}
		}
	}
}

void CalcHessFor3DPoint(vec3 const & Point,
	vec3 const & DelXYZ,
	VolExtentIndexWeights_s & VolInfo,
	Boolean_t IsPeriodic,
	mat33 & Hess,
	vector<FieldDataPointer_c> const & VarReadPtrs,
	GPType_e const & CalcType,
	void * Params)
{
	vector<vec3> GradValues(3);

	for (int i = 0; i < 3; ++i){
		CalcGradForPoint(Point,
			DelXYZ,
			VolInfo,
			*reinterpret_cast<MultiRootParams_s*>(Params)->BasisVectors,
			i,
			IsPeriodic,
			GradValues[i],
			VarReadPtrs[i],
			GPType_Invalid,
			Params);
		for (int j = 0; j < 3; ++j){
			/*
			*	This if statement makes is so that the symmetrical elements of the
			*	hessian are simply copied from the top-right of the matrix rather
			*	than being recomputed.
			*/
			if (j >= i)
				Hess.at(i, j) = GradValues[i][j];
			else
				Hess.at(i, j) = GradValues[j][i];
		}
	}
}

void CalcHessForDataSet(Boolean_t IsPeriodic, AddOn_pa const & AddOnID){
	TecUtilLockStart(AddOnID);

	EntIndex_t VolZoneNum = ZoneNumByName(CSMZoneName.FullVolume);
	if (VolZoneNum < 0){
		TecUtilDialogErrMsg("Failed to get volume zone");
		return;
	}

	EntIndex_t RhoVarNum = VarNumByNameList({ "Rho", CSMVarName.Dens });
	if (RhoVarNum < 0){
		TecUtilDialogErrMsg("Failed to get density (rho) variable");
		return;
	}

	Boolean_t HasGrad = TRUE, HasHess = TRUE;

	vector<EntIndex_t> GradVarNums(3);
	for (int i = 0; i < 3 && HasGrad; ++i){
		GradVarNums[i] = VarNumByName(CSMVarName.DensGradVec[i]);
		HasGrad = GradVarNums[i] > 0;
	}

	if (!HasGrad){
		CalcGradGradMagForDataset(IsPeriodic, AddOnID);

		for (int i = 0; i < 3; ++i){
			GradVarNums[i] = VarNumByName(CSMVarName.DensGradVec[i]);
			HasGrad = GradVarNums[i] > 0;
		}

		if (!HasGrad){
			TecUtilDialogErrMsg("Failed to get gradient variables");
			return;
		}
	}

	EntIndex_t GradMagVarNum = VarNumByName(CSMVarName.DensGradMag);
	if (GradMagVarNum < 0){
		TecUtilDialogErrMsg("Failed to get gradient magnitude variable");
		return;
	}

	vector<EntIndex_t> HessVarNums(6);
	for (int i = 0; i < 6 && HasHess; ++i){
		HessVarNums[i] = VarNumByName(CSMVarName.DensHessTensor[i]);
		HasHess = HessVarNums[i] > 0;
	}

	Boolean_t PeriodicBC = FALSE;

	vector<EntIndex_t> XYZVarNums(3, -1);
	TecUtilAxisGetVarAssignments(&XYZVarNums[0], &XYZVarNums[1], &XYZVarNums[2]);
	if (XYZVarNums[0] <= 0 || XYZVarNums[1] <= 0 || XYZVarNums[2] <= 0){
		vector<vector<string> > TmpStrs = {
			{ "X", "Y", "Z" }, { "I", "J", "K" }
		};
		Boolean_t VarsFound = FALSE;
		for (int i = 0; i < 2 && !VarsFound; ++i){
			for (int j = 0; j < 3; ++j)
				XYZVarNums[j] = TecUtilVarGetNumByName(TmpStrs[i][j].c_str());
			VarsFound = (XYZVarNums[0] > 0 && XYZVarNums[1] > 0 && XYZVarNums[2] > 0);
		}
		if (!VarsFound){
			TecUtilDialogErrMsg("Failed to get x, y, z variable numbers.");
			TecUtilLockFinish(AddOnID);
			return;
		}
	}

	VolExtentIndexWeights_s VolInfo;
	VolInfo.AddOnID = AddOnID;
	VolInfo.IsPeriodic = PeriodicBC;
	GetVolInfo(VolZoneNum, XYZVarNums, IsPeriodic, VolInfo);

	/*
	*	Using this PointSpacingV123 for the distance in the derivative approximations
	*	assumes regular spacing in the data!!!
	*/

	SYSTEM_INFO sysinfo;
	GetSystemInfo(&sysinfo);

	int numCPU = sysinfo.dwNumberOfProcessors;
	omp_set_num_threads(numCPU);

	FieldDataPointer_c RhoPtr;
	vector<FieldDataPointer_c> GradPtrs(3);
	vector<FieldDataPointer_c> HessPtrs(6);

	Boolean_t IsOk = TRUE;

	TecUtilDataLoadBegin();

	RhoPtr.InitializeReadPtr(VolZoneNum, RhoVarNum);

	for (int i = 0; i < 3 && IsOk; ++i){
		IsOk = GradPtrs[i].InitializeReadPtr(VolZoneNum, GradVarNums[i]);
	}

	if (HasHess){
		TecUtilDialogErrMsg("1 or more Hessian variables already exist. Delete them and try again.");
		return;
	}

	mat33 I = eye<mat>(3, 3);
	vector<MultiRootParams_s> TmpParams(numCPU);
	for (int i = 0; i < numCPU; ++i){
		TmpParams[i].CalcType = GPType_NormalPlaneEberlyCP;
		// 		TmpParams[i].BasisVectors = &VolInfo.BasisVectors;
		TmpParams[i].BasisVectors = &I;
		TmpParams[i].HasHess = HasHess;
		TmpParams[i].IsPeriodic = PeriodicBC;
		TmpParams[i].RhoPtr = &RhoPtr;
		TmpParams[i].GradPtrs = &GradPtrs;
		TmpParams[i].HessPtrs = &HessPtrs;
		TmpParams[i].VolInfo = new VolExtentIndexWeights_s;
		*TmpParams[i].VolInfo = VolInfo;
	}

	/*
	*	Calculate the hessian for each point (node) in the system
	*/

	EntIndex_t NumZones = TecUtilDataSetGetNumZones();
	vector<FieldDataType_e> DataType;
	DataType.resize(NumZones, FieldDataType_Float);

	vector<ValueLocation_e> DataLoc(NumZones, ValueLocation_Nodal);

	for (int i = 0; i < 6; ++i){
		EntIndex_t NewVarNum;

		ArgList_pa Args = TecUtilArgListAlloc();
		TecUtilArgListAppendString(Args, SV_NAME, CSMVarName.DensHessTensor[i].c_str());
		TecUtilArgListAppendArray(Args, SV_VARDATATYPE, DataType.data());
		TecUtilArgListAppendArray(Args, SV_VALUELOCATION, DataLoc.data());
		TecUtilArgListAppendInt(Args, SV_SHAREVARWITHALLZONES, FALSE);

		if (TecUtilDataSetAddVarX(Args)){
			Set_pa NewVar = TecUtilSetAlloc(FALSE);
			NewVarNum = TecUtilDataSetGetNumVars();
			TecUtilSetAddMember(NewVar, NewVarNum, FALSE);
			TecUtilStateChanged(StateChange_VarsAdded, reinterpret_cast<ArbParam_t>(NewVar));
			TecUtilSetDealloc(&NewVar);
		}
		TecUtilArgListDealloc(&Args);

		HessPtrs[i].InitializeWritePtr(VolZoneNum, NewVarNum);
	}

	string const TmpStr = "Calculating Hessian of rho";
	StatusLaunch(TmpStr.c_str(), AddOnID, TRUE);

	int NumCompleted = 0, NumTotal = VolInfo.MaxIJK[2] / numCPU;

	int HessIndices[3][3] = {
		{ 0, 1, 2 },
		{ 1, 3, 4 },
		{ 2, 4, 5 }
	};

#pragma omp parallel for
	for (int k = 1; k <= VolInfo.MaxIJK[2]; ++k){
		int ThreadNum = omp_get_thread_num();
		LgIndex_t Ind;
		vec3 Pos;
		mat33 TmpHess;
		if (ThreadNum == 0){
			if (!StatusUpdate(NumCompleted, NumTotal, TmpStr, AddOnID)){
				IsOk = FALSE;
#pragma omp flush (IsOk)
			}
			else
				NumCompleted++;
		}
#pragma omp flush (IsOk)
		for (int j = 1; j <= VolInfo.MaxIJK[1] && IsOk; ++j){
			for (int i = 1; i <= VolInfo.MaxIJK[0]; ++i){
				Ind = IndexFromIJK(i, j, k, VolInfo.MaxIJK[0], VolInfo.MaxIJK[1]) - 1;
				Pos = VolInfo.MinXYZ + VolInfo.PointSpacingV123 % vec3(vector<double>({ i - 1., j - 1., k - 1. }).data());
// 				Pos = VolInfo.MinXYZ + VolInfo.BasisVectors * vec3(vector<double>({ double(i - 1), double(j - 1), double(k - 1) }).data());
				CalcHessFor3DPoint(Pos,
					VolInfo.PointSpacingV123,
					*TmpParams[ThreadNum].VolInfo,
					IsPeriodic,
					TmpHess,
					GradPtrs,
					GPType_Invalid,
					reinterpret_cast<MultiRootParams_s*>(&TmpParams[ThreadNum]));
				for (int ii = 0; ii < 3; ++ii){
					for (int jj = ii; jj < 3; ++jj){
						HessPtrs[HessIndices[ii][jj]].Write(Ind, TmpHess.at(ii, jj));
					}
				}
			}
		}
	}

	for (auto & i : TmpParams) delete i.VolInfo;

	StatusDrop(AddOnID);

	TecUtilDataLoadEnd();

	TecUtilLockFinish(AddOnID);
}

Boolean_t CalcEigenSystemForNode(int ii,
	int jj,
	int kk,
	vec3 & EigenValues,
	mat33 & EigenVectors,
	MultiRootParams_s & RootParams)
{
	Boolean_t IsOk = TRUE;

	mat33 Hessian;

	// 	/*
	// 	*	Prepare gsl data structures for the Hessian and eigen vector matrices
	// 	*	and eigen value vector.
	// 	*/
	// 	gsl_matrix * Hess = gsl_matrix_alloc(3, 3);
	// 	gsl_matrix * EigVecs = gsl_matrix_alloc(3, 3);
	// 	gsl_vector * EigVals = gsl_vector_alloc(3);

	/*
	*	Populate the Hessian matrix.
	*/

	if (RootParams.HasHess){
		/*
		*	Analytical Hessian available, so use that.
		*/

		int HessIndices[3][3] = {
			{ 0, 1, 2 },
			{ 1, 3, 4 },
			{ 2, 4, 5 }
		};

		int Index = IndexFromIJK(ii, jj, kk, RootParams.VolInfo->MaxIJK[0], RootParams.VolInfo->MaxIJK[1], RootParams.VolInfo->MaxIJK[2], RootParams.IsPeriodic) - 1;

		for (int i = 0; i < 3; ++i)
		for (int j = 0; j < 3; ++j)
			Hessian.at(i, j) = RootParams.HessPtrs->at(HessIndices[i][j])[Index];
		// 				gsl_matrix_set(Hess, i, j, RootParams.HessPtrs->at(HessIndices[i][j])[Index]);
	}
	else{
		/*
		*	No analytical Hessian, so need to find derivative numerically.
		*	Need to do it manually, since the GSL solver doesn't know not to
		*	go beyond the bounds of the system.
		*/
		// 		mat33 OutHess;
		vector<double> Vals(5);
		vector<int> DirInd(5);
		CalcHessForNode(ii,
			jj,
			kk,
			RootParams.VolInfo->PointSpacingV123 * 2,
			RootParams.VolInfo->PointSpacingV123 * 12,
			Vals, DirInd,
			RootParams.VolInfo->MaxIJK,
			RootParams.IsPeriodic,
			*RootParams.RhoPtr,
			*RootParams.GradPtrs,
			Hessian);

		// 		for (int i = 0; i < 3; ++i)
		// 			for (int j = 0; j < 3; ++j)
		// 				gsl_matrix_set(Hess, i, j, OutHess[i][j]);
	}

	// 	// Setup the GSL eigensystem workspace
	// 	gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc(3);
	// 
	// 	// Solve the eigensystem
	// 	gsl_eigen_symmv(Hess, EigVals, EigVecs, w);
	// 	/*
	// 	*	Simultaneously sort the eigen values and vectors in ascending order according
	// 	*	to the eigenvalues.
	// 	*/
	// 	gsl_eigen_symmv_sort(EigVals, EigVecs, GSL_EIGEN_SORT_VAL_ASC);
	// 
	// 	/*
	// 	*	Store the results in the CSM vector and matrix provided
	// 	*/
	// 	for (int i = 0; i < 3; ++i){
	// 		EigenValues[i] = gsl_vector_get(EigVals, i);
	// 		for (int j = 0; j < 3; ++j)
	// 			EigenVectors[j][i] = gsl_matrix_get(EigVecs, i, j);
	// 	}
	// 
	// 	/*
	// 	*	Clear the workspace
	// 	*/
	// 	gsl_eigen_symmv_free(w);
	// 	gsl_matrix_free(Hess);
	// 	gsl_matrix_free(EigVecs);
	// 	gsl_vector_free(EigVals);
	// 	
	eig_sym(EigenValues, EigenVectors, Hessian);

	EigenVectors = mat33(normalise(EigenVectors, 2, 1));
// 	EigenVectors = mat33(normalise(EigenVectors));

	// 	EigenVectors = mat33(EigenVectors.t());

	return IsOk;
}

Boolean_t CalcEigenSystemForPoint(vec3 & Point,
	vec & EigenValues,
	mat & EigenVectors,
	MultiRootParams_s & RootParams)
{
	Boolean_t IsOk = TRUE;

	mat33 Hessian;

	if (RootParams.AgitationFactor < 0 || RootParams.VolInfo2 == nullptr || RootParams.RhoPtr2 == nullptr) {
		if (RootParams.Index < 0 && !SetIndexAndWeightsForPoint(Point, *RootParams.VolInfo))
			// 	if (!SetIndexAndWeightsForPoint(Point, *RootParams.VolInfo))
			return FALSE;

		int NumDirs = sqrt(RootParams.BasisVectors->size());


		// 	/*
		// 	*	Prepare gsl data structures for the Hessian and eigen vector matrices
		// 	*	and eigen value vector.
		// 	*/
		// 	gsl_matrix * Hess = gsl_matrix_alloc(NumDirs, NumDirs);
		// 	gsl_matrix * EigVecs = gsl_matrix_alloc(NumDirs, NumDirs);
		// 	gsl_vector * EigVals = gsl_vector_alloc(NumDirs);

		/*
		*	Populate the Hessian matrix.
		*/

		if (RootParams.HasHess && NumDirs >= 3) {
			/*
			*	Analytical Hessian available, so use that.
			*/

			int HessIndices[3][3] = {
				{ 0, 1, 2 },
				{ 1, 3, 4 },
				{ 2, 4, 5 }
			};

			for (int i = 0; i < 3; ++i)
				for (int j = 0; j < 3; ++j)
					if (RootParams.Index >= 0)
						Hessian.at(i, j) = RootParams.HessPtrs->at(HessIndices[i][j])[RootParams.Index];
			// 					gsl_matrix_set(Hess, i, j, RootParams.HessPtrs->at(HessIndices[i][j])[RootParams.Index]);
					else
						Hessian.at(i, j) = ValByCurrentIndexAndWeightsFromRawPtr(*RootParams.VolInfo, RootParams.HessPtrs->at(HessIndices[i][j]));
			// 					gsl_matrix_set(Hess, i, j, ValByCurrentIndexAndWeightsFromRawPtr(*RootParams.VolInfo, RootParams.HessPtrs->at(HessIndices[i][j])));
		}
		else {
			/*
			*	No analytical Hessian, so need to find derivative numerically.
			*	Need to do it manually, since the GSL solver doesn't know not to
			*	go beyond the bounds of the system.
			*/
			if (!RootParams.HasGrad) {
				CalcHessForPoint(Point,
					RootParams.VolInfo->PointSpacingV123,
					*RootParams.VolInfo,
					*RootParams.BasisVectors,
					RootParams.IsPeriodic,
					Hessian,
					*RootParams.RhoPtr,
					GPType_Invalid,
					reinterpret_cast<void*>(&RootParams));
			}
			else {
				if (NumDirs >= 3) {
					if (RootParams.HasGrad) {
						CalcHessFor3DPoint(Point,
							RootParams.VolInfo->PointSpacingV123,
							*RootParams.VolInfo,
							RootParams.IsPeriodic,
							Hessian,
							*RootParams.GradPtrs,
							GPType_Invalid,
							reinterpret_cast<void*>(&RootParams));
					}


					// 			for (int i = 0; i < NumDirs; ++i)
					// 				for (int j = 0; j < NumDirs; ++j)
					// 					gsl_matrix_set(Hess, i, j, OutHess[i][j]);
				}
				else {
					CalcHessForPoint(Point,
						RootParams.VolInfo->PointSpacingV123,
						*RootParams.VolInfo,
						*RootParams.BasisVectors,
						RootParams.IsPeriodic,
						Hessian,
						RootParams.GradPtrs->at(0),
						GPType_Invalid,
						reinterpret_cast<void*>(&RootParams));

					// 			for (int i = 0; i < NumDirs; ++i)
					// 				for (int j = 0; j < NumDirs; ++j)
					// 					gsl_matrix_set(Hess, i, j, OutHess[i][j]);
				}
			}

		}
	}
	else{
		if (RootParams.Index < 0 && !SetIndexAndWeightsForPoint(Point, *RootParams.VolInfo2))
			// 	if (!SetIndexAndWeightsForPoint(Point, *RootParams.VolInfo))
			return FALSE;

		int NumDirs = sqrt(RootParams.VolInfo2->BasisVectors.size());

		/*
		*	Populate the Hessian matrix.
		*/

		CalcHessForPoint(Point,
			RootParams.VolInfo2->PointSpacingV123,
			*RootParams.VolInfo2,
			RootParams.VolInfo2->BasisVectors,
			FALSE,
			Hessian,
			*RootParams.RhoPtr2,
			GPType_Invalid,
			reinterpret_cast<void*>(&RootParams));
	}

	// 	// Setup the GSL eigensystem workspace
	// 	gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc(NumDirs);
	// 
	// 	// Solve the eigensystem
	// 	gsl_eigen_symmv(Hess, EigVals, EigVecs, w);
	// 	/*
	// 	*	Simultaneously sort the eigen values and vectors in ascending order according
	// 	*	to the eigenvalues.
	// 	*/
	// 	gsl_eigen_symmv_sort(EigVals, EigVecs, GSL_EIGEN_SORT_VAL_ASC);
	// 
	// 	/*
	// 	*	Store the results in the CSM vector and matrix provided
	// 	*/
	// 	for (int i = 0; i < NumDirs; ++i){
	// 		EigenValues[i] = gsl_vector_get(EigVals, i);
	// 		for (int j = 0; j < NumDirs; ++j)
	// 			EigenVectors[j][i] = gsl_matrix_get(EigVecs, i, j);
	// 	}
	// 
	// 	/*
	// 	*	Clear the workspace
	// 	*/
	// 	gsl_eigen_symmv_free(w);
	// 	gsl_matrix_free(Hess);
	// 	gsl_matrix_free(EigVecs);
	// 	gsl_vector_free(EigVals);
	// 	
	eig_sym(EigenValues, EigenVectors, Hessian, "std");

	EigenVectors = mat33(normalise(EigenVectors, 2, 1));
// 	EigenVectors = mat33(normalise(EigenVectors));

	return IsOk;
}

void CalcEigenSystemForDataSet(Boolean_t IsPeriodic, AddOn_pa const & AddOnID){
	TecUtilLockStart(AddOnID);

	EntIndex_t VolZoneNum = ZoneNumByName(CSMZoneName.FullVolume);
	if (VolZoneNum < 0){
		TecUtilDialogErrMsg("Failed to get volume zone");
		return;
	}

	EntIndex_t RhoVarNum = VarNumByNameList({ "Rho", CSMVarName.Dens });

	vector<EntIndex_t> GradVarNums(3);
	for (int i = 0; i < 3; ++i){
		GradVarNums[i] = VarNumByName(CSMVarName.DensGradVec[i]);
		if (GradVarNums[i] < 0){
			TecUtilDialogErrMsg("Failed to get gradient variable");
			return;
		}
	}

	EntIndex_t GradMagVarNum = VarNumByName(CSMVarName.DensGradMag);
	if (GradMagVarNum < 0){
		TecUtilDialogErrMsg("Failed to get gradient magnitude variable");
		return;
	}

	Boolean_t HasHess = TRUE;
	vector<EntIndex_t> HessVarNums(6);
	for (int i = 0; i < 6 && HasHess; ++i){
		HessVarNums[i] = VarNumByName(CSMVarName.DensHessTensor[i]);
		HasHess = HessVarNums[i] > 0;
	}

	Boolean_t PeriodicBC = FALSE;

	vector<EntIndex_t> XYZVarNums(3, -1);
	TecUtilAxisGetVarAssignments(&XYZVarNums[0], &XYZVarNums[1], &XYZVarNums[2]);
	if (XYZVarNums[0] <= 0 || XYZVarNums[1] <= 0 || XYZVarNums[2] <= 0){
		vector<vector<string> > TmpStrs = {
			{ "X", "Y", "Z" }, { "I", "J", "K" }
		};
		Boolean_t VarsFound = FALSE;
		for (int i = 0; i < 2 && !VarsFound; ++i){
			for (int j = 0; j < 3; ++j)
				XYZVarNums[j] = TecUtilVarGetNumByName(TmpStrs[i][j].c_str());
			VarsFound = (XYZVarNums[0] > 0 && XYZVarNums[1] > 0 && XYZVarNums[2] > 0);
		}
		if (!VarsFound){
			TecUtilDialogErrMsg("Failed to get x, y, z variable numbers.");
			TecUtilLockFinish(AddOnID);
			return;
		}
	}

	VolExtentIndexWeights_s VolInfo;
	VolInfo.AddOnID = AddOnID;
	VolInfo.IsPeriodic = PeriodicBC;
	GetVolInfo(VolZoneNum, XYZVarNums, PeriodicBC, VolInfo);

	SYSTEM_INFO sysinfo;
	GetSystemInfo(&sysinfo);

	int numCPU = sysinfo.dwNumberOfProcessors;
	omp_set_num_threads(numCPU);

	FieldDataPointer_c RhoPtr;
	vector<FieldDataPointer_c> GradPtrs(3);
	vector<FieldDataPointer_c> HessPtrs;

	Boolean_t IsOk = TRUE;

	TecUtilDataLoadBegin();

	RhoPtr.InitializeReadPtr(VolZoneNum, RhoVarNum);

	for (int i = 0; i < 3 && IsOk; ++i){
		IsOk = GradPtrs[i].InitializeReadPtr(VolZoneNum, GradVarNums[i]);
	}

	if (HasHess){
		HessPtrs.resize(6);
		for (int i = 0; i < 6 && IsOk; ++i)
			IsOk = HessPtrs[i].InitializeReadPtr(VolZoneNum, HessVarNums[i]);
	}

	vector<MultiRootParams_s> TmpParams(numCPU);
	mat33 I = eye<mat>(3, 3);
	for (int i = 0; i < numCPU; ++i){
		TmpParams[i].CalcType = GPType_NormalPlaneEberlyCP;
		TmpParams[i].BasisVectors = &I;
		TmpParams[i].HasHess = HasHess;
		TmpParams[i].IsPeriodic = PeriodicBC;
		TmpParams[i].RhoPtr = &RhoPtr;
		TmpParams[i].GradPtrs = &GradPtrs;
		TmpParams[i].HessPtrs = &HessPtrs;
		TmpParams[i].VolInfo = new VolExtentIndexWeights_s;
		*TmpParams[i].VolInfo = VolInfo;
	}

	vector<string> VarNames = {
		"Eigenvector 1-1",
		"Eigenvector 1-2",
		"Eigenvector 1-3",
		"Eigenvector 2-1",
		"Eigenvector 2-2",
		"Eigenvector 2-3",
		"Eigenvector 3-1",
		"Eigenvector 3-2",
		"Eigenvector 3-3",
		"Eigenvalue 1",
		"Eigenvalue 2",
		"Eigenvalue 3"
	};
	vector<FieldDataPointer_c> VarPtrs(VarNames.size());

	EntIndex_t NumZones = TecUtilDataSetGetNumZones();
	vector<FieldDataType_e> DataType;
	DataType.resize(NumZones, FieldDataType_Float);

	vector<ValueLocation_e> DataLoc(NumZones, ValueLocation_Nodal);

	EntIndex_t NewVarNum;

	for (int i = 0; i < VarNames.size(); ++i){
		if (VarNumByName(VarNames[i]) > 0){
			TecUtilDialogErrMsg("1 or more eigenvector or eigenvalue variables exist. Delete them and try again.");
			return;
		}

		ArgList_pa Args = TecUtilArgListAlloc();
		TecUtilArgListAppendString(Args, SV_NAME, VarNames[i].c_str());
		TecUtilArgListAppendArray(Args, SV_VARDATATYPE, DataType.data());
		TecUtilArgListAppendArray(Args, SV_VALUELOCATION, DataLoc.data());
		TecUtilArgListAppendInt(Args, SV_SHAREVARWITHALLZONES, FALSE);

		if (TecUtilDataSetAddVarX(Args)){
			Set_pa NewVar = TecUtilSetAlloc(FALSE);
			NewVarNum = TecUtilDataSetGetNumVars();
			TecUtilSetAddMember(NewVar, NewVarNum, FALSE);
			TecUtilStateChanged(StateChange_VarsAdded, reinterpret_cast<ArbParam_t>(NewVar));
			TecUtilSetDealloc(&NewVar);
		}
		TecUtilArgListDealloc(&Args);

		VarPtrs[i].InitializeWritePtr(VolZoneNum, NewVarNum);
	}

	string const TmpStr = "Calculating Eigen system";
	StatusLaunch(TmpStr.c_str(), AddOnID, TRUE);

	int NumCompleted = 0, NumTotal = VolInfo.MaxIJK[2] / numCPU;

#pragma omp parallel for
	for (int k = 1; k <= VolInfo.MaxIJK[2]; ++k){
		int ThreadNum = omp_get_thread_num();
		LgIndex_t Ind;
		vec3 Pos;
		vec3 EigenValues;
		mat33 EigenVectors;
		if (ThreadNum == 0){
			if (!StatusUpdate(NumCompleted, NumTotal, TmpStr, AddOnID)){
				IsOk = FALSE;
#pragma omp flush (IsOk)
			}
			else
				NumCompleted++;
		}
#pragma omp flush (IsOk)
		for (int j = 1; j <= VolInfo.MaxIJK[1] && IsOk; ++j){
			for (int i = 1; i <= VolInfo.MaxIJK[0]; ++i){
				Ind = IndexFromIJK(i, j, k, VolInfo.MaxIJK[0], VolInfo.MaxIJK[1]) - 1;
				Pos = VolInfo.MinXYZ + VolInfo.PointSpacingV123 % vec3(vector<double>({ i - 1., j - 1., k - 1. }).data());
				CalcEigenSystemForPoint(Pos,
					EigenValues,
					EigenVectors,
					TmpParams[ThreadNum]);
				for (int ii = 0; ii < 3; ++ii){
					for (int jj = 0; jj < 3; ++jj){
						VarPtrs[3 * ii + jj].Write(Ind, EigenVectors.at(ii, jj));
					}
				}
				for (int ii = 0; ii < 3; ++ii){
					VarPtrs[9 + ii].Write(Ind, EigenValues[ii]);
				}
			}
		}
	}

	for (auto & i : TmpParams) delete i.VolInfo;

	StatusDrop(AddOnID);

	TecUtilDataLoadEnd();

	TecUtilLockFinish(AddOnID);
}

void CalcEigenvecDotGradForDataSet(Boolean_t IsPeriodic,
	Boolean_t NormalizeGrad,
	AddOn_pa const & AddOnID)
{
	TecUtilLockStart(AddOnID);

	EntIndex_t VolZoneNum = ZoneNumByName(CSMZoneName.FullVolume);
	if (VolZoneNum < 0){
		TecUtilDialogErrMsg("Failed to get volume zone");
		return;
	}

	EntIndex_t RhoVarNum = VarNumByNameList({ "Rho", CSMVarName.Dens });

	vector<EntIndex_t> GradVarNums(3);
	for (int i = 0; i < 3; ++i){
		GradVarNums[i] = VarNumByName(CSMVarName.DensGradVec[i]);
		if (GradVarNums[i] < 0){
			TecUtilDialogErrMsg("Failed to get gradient variable");
			return;
		}
	}

	EntIndex_t GradMagVarNum = VarNumByName(CSMVarName.DensGradMag);
	if (GradMagVarNum < 0){
		TecUtilDialogErrMsg("Failed to get gradient magnitude variable");
		return;
	}

	Boolean_t HasHess = TRUE;
	vector<EntIndex_t> HessVarNums(6);
	for (int i = 0; i < 6 && HasHess; ++i){
		HessVarNums[i] = VarNumByName(CSMVarName.DensHessTensor[i]);
		HasHess = HessVarNums[i] > 0;
	}

	Boolean_t PeriodicBC = FALSE;

	vec3 DelXYZ;
	vector<EntIndex_t> XYZVarNums(3, -1);
	TecUtilAxisGetVarAssignments(&XYZVarNums[0], &XYZVarNums[1], &XYZVarNums[2]);
	if (XYZVarNums[0] <= 0 || XYZVarNums[1] <= 0 || XYZVarNums[2] <= 0){
		vector<vector<string> > TmpStrs = {
			{ "X", "Y", "Z" }, { "I", "J", "K" }
		};
		Boolean_t VarsFound = FALSE;
		for (int i = 0; i < 2 && !VarsFound; ++i){
			for (int j = 0; j < 3; ++j)
				XYZVarNums[j] = TecUtilVarGetNumByName(TmpStrs[i][j].c_str());
			VarsFound = (XYZVarNums[0] > 0 && XYZVarNums[1] > 0 && XYZVarNums[2] > 0);
		}
		if (!VarsFound){
			TecUtilDialogErrMsg("Failed to get x, y, z variable numbers.");
			TecUtilLockFinish(AddOnID);
			return;
		}
	}

	VolExtentIndexWeights_s VolInfo;
	VolInfo.AddOnID = AddOnID;
	VolInfo.IsPeriodic = PeriodicBC;
	GetVolInfo(VolZoneNum, XYZVarNums, PeriodicBC, VolInfo);

	SYSTEM_INFO sysinfo;
	GetSystemInfo(&sysinfo);

	int numCPU = sysinfo.dwNumberOfProcessors;
	omp_set_num_threads(numCPU);

	FieldDataPointer_c RhoPtr;
	vector<FieldDataPointer_c> GradPtrs(3);
	vector<FieldDataPointer_c> HessPtrs;

	Boolean_t IsOk = TRUE;

	TecUtilDataLoadBegin();

	RhoPtr.InitializeReadPtr(VolZoneNum, RhoVarNum);

	for (int i = 0; i < 3 && IsOk; ++i){
		IsOk = GradPtrs[i].InitializeReadPtr(VolZoneNum, GradVarNums[i]);
	}

	if (HasHess){
		HessPtrs.resize(6);
		for (int i = 0; i < 6 && IsOk; ++i)
			IsOk = HessPtrs[i].InitializeReadPtr(VolZoneNum, HessVarNums[i]);
	}

	vector<MultiRootParams_s> TmpParams(numCPU);
	mat33 I = eye<mat>(3, 3);
	for (int i = 0; i < numCPU; ++i){
		TmpParams[i].CalcType = GPType_NormalPlaneEberlyCP;
		TmpParams[i].BasisVectors = &I;
		TmpParams[i].HasHess = HasHess;
		TmpParams[i].IsPeriodic = PeriodicBC;
		TmpParams[i].RhoPtr = &RhoPtr;
		TmpParams[i].GradPtrs = &GradPtrs;
		TmpParams[i].HessPtrs = &HessPtrs;
		TmpParams[i].VolInfo = new VolExtentIndexWeights_s;
		*TmpParams[i].VolInfo = VolInfo;
	}

	vector<string> VarNames = {
		"Gradient dot Eigenvector 1",
		"Gradient dot Eigenvector 2",
		"Gradient dot Eigenvector 3"
	};
	vector<FieldDataPointer_c> VarPtrs(VarNames.size());

	EntIndex_t NumZones = TecUtilDataSetGetNumZones();
	vector<FieldDataType_e> DataType;
	DataType.resize(NumZones, FieldDataType_Float);

	vector<ValueLocation_e> DataLoc(NumZones, ValueLocation_Nodal);

	EntIndex_t NewVarNum;

	for (int i = 0; i < 3; ++i){
		if (VarNumByName(VarNames[i]) > 0){
			TecUtilDialogErrMsg("1 or more grad-dot-eigenvector variables exist. Delete them and try again.");
			return;
		}

		ArgList_pa Args = TecUtilArgListAlloc();
		TecUtilArgListAppendString(Args, SV_NAME, VarNames[i].c_str());
		TecUtilArgListAppendArray(Args, SV_VARDATATYPE, DataType.data());
		TecUtilArgListAppendArray(Args, SV_VALUELOCATION, DataLoc.data());
		TecUtilArgListAppendInt(Args, SV_SHAREVARWITHALLZONES, FALSE);

		if (TecUtilDataSetAddVarX(Args)){
			Set_pa NewVar = TecUtilSetAlloc(FALSE);
			NewVarNum = TecUtilDataSetGetNumVars();
			TecUtilSetAddMember(NewVar, NewVarNum, FALSE);
			TecUtilStateChanged(StateChange_VarsAdded, reinterpret_cast<ArbParam_t>(NewVar));
			TecUtilSetDealloc(&NewVar);
		}
		TecUtilArgListDealloc(&Args);

		VarPtrs[i].InitializeWritePtr(VolZoneNum, NewVarNum);
	}

	string const TmpStr = "Calculating dot products of Hessian eigenvectors with gradient of rho";
	StatusLaunch(TmpStr.c_str(), AddOnID, TRUE);

	int NumCompleted = 0, NumTotal = VolInfo.MaxIJK[2] / numCPU;

#pragma omp parallel for
	for (int k = 1; k <= VolInfo.MaxIJK[2]; ++k){
		int ThreadNum = omp_get_thread_num();
		LgIndex_t Ind;
		vec3 Pos;
		vec3 DotProducts;
		if (ThreadNum == 0){
			if (!StatusUpdate(NumCompleted, NumTotal, TmpStr, AddOnID)){
				IsOk = FALSE;
#pragma omp flush (IsOk)
			}
			else
				NumCompleted++;
		}
#pragma omp flush (IsOk)
		for (int j = 1; j <= VolInfo.MaxIJK[1] && IsOk; ++j){
			for (int i = 1; i <= VolInfo.MaxIJK[0]; ++i){
				Ind = IndexFromIJK(i, j, k, VolInfo.MaxIJK[0], VolInfo.MaxIJK[1]) - 1;
				Pos = VolInfo.MinXYZ + VolInfo.PointSpacingV123 % vec3(vector<double>({ i - 1., j - 1., k - 1. }).data());
				CalcEigenvecDotGradForPoint(Pos, DotProducts, NormalizeGrad, TmpParams[ThreadNum]);
				for (int ii = 0; ii < 3; ++ii){
					VarPtrs[ii].Write(Ind, DotProducts[ii]);
				}
			}
		}
	}

	for (auto & i : TmpParams) delete i.VolInfo;

	StatusDrop(AddOnID);

	TecUtilDataLoadEnd();

	TecUtilLockFinish(AddOnID);
}

void CalcEigenvecDotGradForPoint(vec3 Point,
	vec3 & DotProducts,
	vec3 & EigenValues,
	mat33 & EigenVectors,
	Boolean_t NormalizeGrad,
	MultiRootParams_s & RootParams)
{
	vec3 Gradient;

	SetIndexAndWeightsForPoint(Point, *RootParams.VolInfo);

	CalcEigenSystemForPoint(Point,
		EigenValues,
		EigenVectors,
		RootParams);

	for (int i = 0; i < 3; ++i)
		Gradient[i] = ValByCurrentIndexAndWeightsFromRawPtr(*RootParams.VolInfo, RootParams.GradPtrs->at(i));

	if (NormalizeGrad)
		Gradient = normalise(Gradient);
	else{
		for (int i = 0; i < 3; ++i)
			EigenVectors[i] *= EigenValues[i];
	}


	for (int i = 0; i < 3; ++i)
		DotProducts[i] = dot(Gradient, EigenVectors.col(i));
}

void CalcEigenvecDotGradForPoint(vec3 Point,
	vec3 & DotProducts,
	Boolean_t NormalizeGrad,
	MultiRootParams_s & RootParams)
{
	vec3 EigenValues;
	mat33 EigenVectors;

	CalcEigenvecDotGradForPoint(Point, DotProducts, EigenValues, EigenVectors, NormalizeGrad, RootParams);
}

/*
*	One Eberly definition for a point on a 1-ridge is a point such that
*	two of the dot (inner) products of the Hessian's eigen vectors with
*	the gradient are zero, and dot product of the third eigen vector with
*	the gradient is one.
*	Because the eigen vector whose corresponding eigen value has the most
*	positive value (compared to the other two eigen values) is guaranteed to
*	be the eigenvector pointing in the direction of the gradient, I only
*	need to check the dot product of this eigen vector with the gradient.
*	If the dot product is 1, then the other two dot products, of the other
*	two eigen vectors with the gradient, are guaranteed to be zero, because
*	the eigen vectors are orthonormal.
*	This function is being used in a multidimensional root finder, so it will
*	return 1 - x, where x is the dot product of the most positive eigen value's
*	eigen vector with the gradient.
*/
double Eberly1RidgeFunction(vec3 & Point,
	double const & RhoCutoff,
	Boolean_t NormalizeGrad,
	MultiRootParams_s & RootParams)
{
	double Val = 0.0;

	SetIndexAndWeightsForPoint(Point, *RootParams.VolInfo);

	double RhoValue = ValByCurrentIndexAndWeightsFromRawPtr(*RootParams.VolInfo, *RootParams.RhoPtr);

	if (RhoValue > RhoCutoff){
		vec3 DotPdts;

		CalcEigenvecDotGradForPoint(Point, DotPdts, NormalizeGrad, RootParams);

		vector<double> DotProducts(3);

		for (int i = 0; i < 3; ++i)
			DotProducts[i] = abs(DotPdts[i]);

		std::sort(DotProducts.begin(), DotProducts.end());

		for (int i = 0; i < 2; ++i)
			Val += DotProducts[i];
	}
	else Val = 1e100;

	return Val;
}

double Eberly2RidgeFunction(vec3 & Point,
	double const & RhoCutoff,
	Boolean_t NormalizeGrad,
	MultiRootParams_s & RootParams)
{
	double Val = 0.0;

	SetIndexAndWeightsForPoint(Point, *RootParams.VolInfo);

	double RhoValue = ValByCurrentIndexAndWeightsFromRawPtr(*RootParams.VolInfo, *RootParams.RhoPtr);

	if (RhoValue > RhoCutoff){
		vec3 DotPdts;

		CalcEigenvecDotGradForPoint(Point, DotPdts, NormalizeGrad, RootParams);

		vector<double> DotProducts(3);

		for (int i = 0; i < 3; ++i)
			DotProducts[i] = abs(DotPdts[i]);

		std::sort(DotProducts.begin(), DotProducts.end());

		for (int i = 0; i < 1; ++i)
			Val += DotProducts[i];
	}
	else Val = 1e100;

	return Val;
}

double NEBForceFunction(vec3 & Point,
	MultiRootParams_s & RootParams)
{
	double Val = 0.0;

	double DispTerm, GradTerm;

	vec2 Grad;

	CalcGradForPoint(Point, RootParams.VolInfo->PointSpacingV123,
		*RootParams.VolInfo,
		*RootParams.BasisVectors,
		0,
		RootParams.IsPeriodic,
		Grad,
		*RootParams.RhoPtr,
		GPType_NormalPlaneEberlyCP,
		reinterpret_cast<void*>(&RootParams));

	DispTerm = RootParams.KDisp *Distance(Point, (*RootParams.EquilPos));
	GradTerm = RootParams.KGrad * norm(Grad);

	Val = DispTerm + GradTerm;

	return Val;
}

void CalcEberlyFunctions(Boolean_t IsPeriodic, AddOn_pa const & AddOnID, double const & RhoCutoff){
	TecUtilLockStart(AddOnID);

	EntIndex_t VolZoneNum = ZoneNumByName(CSMZoneName.FullVolume);
	if (VolZoneNum < 0){
		TecUtilDialogErrMsg("Failed to get volume zone");
		return;
	}

	EntIndex_t RhoVarNum = VarNumByNameList({ "Rho", CSMVarName.Dens });

	vector<EntIndex_t> GradVarNums(3);
	for (int i = 0; i < 3; ++i){
		GradVarNums[i] = VarNumByName(CSMVarName.DensGradVec[i]);
		if (GradVarNums[i] < 0){
			TecUtilDialogErrMsg("Failed to get gradient variable");
			return;
		}
	}

	EntIndex_t GradMagVarNum = VarNumByName(CSMVarName.DensGradMag);
	if (GradMagVarNum < 0){
		TecUtilDialogErrMsg("Failed to get gradient magnitude variable");
		return;
	}

	Boolean_t HasHess = TRUE;
	vector<EntIndex_t> HessVarNums(6);
	for (int i = 0; i < 6 && HasHess; ++i){
		HessVarNums[i] = VarNumByName(CSMVarName.DensHessTensor[i]);
		HasHess = HessVarNums[i] > 0;
	}

	Boolean_t PeriodicBC = FALSE;

	LgIndex_t  IMax, JMax, KMax;

	vec3 DelXYZ;
	vector<EntIndex_t> XYZVarNums(3, -1);
	TecUtilAxisGetVarAssignments(&XYZVarNums[0], &XYZVarNums[1], &XYZVarNums[2]);
	if (XYZVarNums[0] <= 0 || XYZVarNums[1] <= 0 || XYZVarNums[2] <= 0){
		vector<vector<string> > TmpStrs = {
			{ "X", "Y", "Z" }, { "I", "J", "K" }
		};
		Boolean_t VarsFound = FALSE;
		for (int i = 0; i < 2 && !VarsFound; ++i){
			for (int j = 0; j < 3; ++j)
				XYZVarNums[j] = TecUtilVarGetNumByName(TmpStrs[i][j].c_str());
			VarsFound = (XYZVarNums[0] > 0 && XYZVarNums[1] > 0 && XYZVarNums[2] > 0);
		}
		if (!VarsFound){
			TecUtilDialogErrMsg("Failed to get x, y, z variable numbers.");
			TecUtilLockFinish(AddOnID);
			return;
		}
	}

	vector<LgIndex_t> MaxIJK(3);
	VolExtentIndexWeights_s VolInfo;
	VolInfo.AddOnID = AddOnID;
	VolInfo.IsPeriodic = PeriodicBC;
	GetVolInfo(VolZoneNum, XYZVarNums, PeriodicBC, VolInfo);

	SYSTEM_INFO sysinfo;
	GetSystemInfo(&sysinfo);

	int numCPU = sysinfo.dwNumberOfProcessors;
	omp_set_num_threads(numCPU);

	FieldDataPointer_c RhoPtr;
	vector<FieldDataPointer_c> GradPtrs(3);
	vector<FieldDataPointer_c> HessPtrs;

	Boolean_t IsOk = TRUE;

	TecUtilDataLoadBegin();

	RhoPtr.InitializeReadPtr(VolZoneNum, RhoVarNum);

	for (int i = 0; i < 3 && IsOk; ++i){
		IsOk = GradPtrs[i].InitializeReadPtr(VolZoneNum, GradVarNums[i]);
	}

	if (HasHess){
		HessPtrs.resize(6);
		for (int i = 0; i < 6 && IsOk; ++i)
			IsOk = HessPtrs[i].InitializeReadPtr(VolZoneNum, HessVarNums[i]);
	}

	vector<MultiRootParams_s> TmpParams(numCPU);
	mat33 I = eye<mat>(3, 3);
	for (auto & i : TmpParams){
		i.CalcType = GPType_NormalPlaneEberlyCP;
		i.BasisVectors = &I;
		i.HasHess = HasHess;
		i.IsPeriodic = PeriodicBC;
		i.RhoPtr = &RhoPtr;
		i.GradPtrs = &GradPtrs;
		i.HessPtrs = &HessPtrs;
		i.VolInfo = new VolExtentIndexWeights_s;
		*i.VolInfo = VolInfo;
	}

	/*
	*	Find sum of smallest two dot products of eigenvectors with gradient
	*	and save value at each point in the volume.
	*/
	{
		EntIndex_t NumZones = TecUtilDataSetGetNumZones();
		vector<FieldDataType_e> DataType;
		DataType.resize(NumZones, FieldDataType_Float);

		vector<ValueLocation_e> DataLoc(NumZones, ValueLocation_Nodal);

		EntIndex_t NewVarNum;

		ArgList_pa Args = TecUtilArgListAlloc();
		if (RhoCutoff > 0)
			TecUtilArgListAppendString(Args, SV_NAME, "Eberly 1-ridge with cutoff");
		else
			TecUtilArgListAppendString(Args, SV_NAME, "Eberly 1-ridge");
		TecUtilArgListAppendArray(Args, SV_VARDATATYPE, DataType.data());
		TecUtilArgListAppendArray(Args, SV_VALUELOCATION, DataLoc.data());
		TecUtilArgListAppendInt(Args, SV_SHAREVARWITHALLZONES, FALSE);

		if (TecUtilDataSetAddVarX(Args)){
			Set_pa NewVar = TecUtilSetAlloc(FALSE);
			NewVarNum = TecUtilDataSetGetNumVars();
			TecUtilSetAddMember(NewVar, NewVarNum, FALSE);
			TecUtilStateChanged(StateChange_VarsAdded, reinterpret_cast<ArbParam_t>(NewVar));
			TecUtilSetDealloc(&NewVar);
		}
		TecUtilArgListDealloc(&Args);

		FieldDataPointer_c VarPtr;
		VarPtr.InitializeWritePtr(VolZoneNum, NewVarNum);

		string const TmpStr = "Calculate Eberly 1-ridge function";
		StatusLaunch(TmpStr.c_str(), AddOnID, TRUE);

		int NumCompleted = 0, NumTotal = VolInfo.MaxIJK[2] / numCPU;

#pragma omp parallel for
		for (int k = 1; k <= VolInfo.MaxIJK[2]; ++k){
			if (omp_get_thread_num() == 0){
				if (!StatusUpdate(NumCompleted, NumTotal, TmpStr, AddOnID)){
					IsOk = FALSE;
#pragma omp flush (IsOk)
				}
				else
					NumCompleted++;
			}
#pragma omp flush (IsOk)
			LgIndex_t Ind;
			vec3 Pos;
			for (int j = 1; j <= VolInfo.MaxIJK[1] && IsOk; ++j){
				for (int i = 1; i <= VolInfo.MaxIJK[0]; ++i){
					Ind = IndexFromIJK(i, j, k, VolInfo.MaxIJK[0], VolInfo.MaxIJK[1]) - 1;
					Pos = VolInfo.MinXYZ + VolInfo.PointSpacingV123 % vec3(vector<double>({ i - 1., j - 1., k - 1. }).data());
					VarPtr.Write(Ind, Eberly1RidgeFunction(Pos, RhoCutoff, TRUE, TmpParams[omp_get_thread_num()]));
				}
			}
		}

		StatusDrop(AddOnID);
	}

	/*
	*	Again, but for 2-ridges
	*/
	{
		EntIndex_t NumZones = TecUtilDataSetGetNumZones();
		vector<FieldDataType_e> DataType;
		DataType.resize(NumZones, FieldDataType_Float);

		vector<ValueLocation_e> DataLoc(NumZones, ValueLocation_Nodal);

		EntIndex_t NewVarNum;

		ArgList_pa Args = TecUtilArgListAlloc();
		if (RhoCutoff > 0)
			TecUtilArgListAppendString(Args, SV_NAME, "Eberly 2-ridge with cutoff");
		else
			TecUtilArgListAppendString(Args, SV_NAME, "Eberly 2-ridge");
		TecUtilArgListAppendArray(Args, SV_VARDATATYPE, DataType.data());
		TecUtilArgListAppendArray(Args, SV_VALUELOCATION, DataLoc.data());
		TecUtilArgListAppendInt(Args, SV_SHAREVARWITHALLZONES, FALSE);

		if (TecUtilDataSetAddVarX(Args)){
			Set_pa NewVar = TecUtilSetAlloc(FALSE);
			NewVarNum = TecUtilDataSetGetNumVars();
			TecUtilSetAddMember(NewVar, NewVarNum, FALSE);
			TecUtilStateChanged(StateChange_VarsAdded, reinterpret_cast<ArbParam_t>(NewVar));
			TecUtilSetDealloc(&NewVar);
		}
		TecUtilArgListDealloc(&Args);

		FieldDataPointer_c VarPtr;
		VarPtr.InitializeWritePtr(VolZoneNum, NewVarNum);

		string const TmpStr = "Calculate Eberly 2-ridge function";
		StatusLaunch(TmpStr.c_str(), AddOnID, TRUE);

		int NumCompleted = 0, NumTotal = VolInfo.MaxIJK[2] / numCPU;

#pragma omp parallel for
		for (int k = 1; k <= VolInfo.MaxIJK[2]; ++k){
			if (omp_get_thread_num() == 0){
				if (!StatusUpdate(NumCompleted, NumTotal, TmpStr, AddOnID)){
					IsOk = FALSE;
#pragma omp flush (IsOk)
				}
				else
					NumCompleted++;
			}
#pragma omp flush (IsOk)
			LgIndex_t Ind;
			vec3 Pos;
			for (int j = 1; j <= VolInfo.MaxIJK[1] && IsOk; ++j){
				for (int i = 1; i <= VolInfo.MaxIJK[0]; ++i){
					Ind = IndexFromIJK(i, j, k, VolInfo.MaxIJK[0], VolInfo.MaxIJK[1]) - 1;
					Pos = VolInfo.MinXYZ + VolInfo.PointSpacingV123 % vec3(vector<double>({ i - 1., j - 1., k - 1. }).data());
					VarPtr.Write(Ind, Eberly2RidgeFunction(Pos, RhoCutoff, TRUE, TmpParams[omp_get_thread_num()]));
				}
			}
		}

		for (auto & i : TmpParams) delete i.VolInfo;

		StatusDrop(AddOnID);
	}

	TecUtilDataLoadEnd();

	TecUtilLockFinish(AddOnID);
}

void MapAllVarsToAllZones(AddOn_pa const & AddOnID)
{
	SYSTEM_INFO sysinfo;
	GetSystemInfo(&sysinfo);

	int numCPU = sysinfo.dwNumberOfProcessors;
	omp_set_num_threads(numCPU);

	TecUtilLockStart(AddOnID);

	EntIndex_t NumVars = TecUtilDataSetGetNumVars();
	EntIndex_t NumZones = TecUtilDataSetGetNumZones();

	EntIndex_t VolZoneNum = ZoneNumByName("Full");

	vector<EntIndex_t> XYZVarNums(3, -1);
	TecUtilAxisGetVarAssignments(&XYZVarNums[0], &XYZVarNums[1], &XYZVarNums[2]);
	if (XYZVarNums[0] <= 0 || XYZVarNums[1] <= 0 || XYZVarNums[2] <= 0){
		vector<vector<string> > TmpStrs = {
			{ "X", "Y", "Z" }, { "I", "J", "K" }
		};
		Boolean_t VarsFound = FALSE;
		for (int i = 0; i < 2 && !VarsFound; ++i){
			for (int j = 0; j < 3; ++j)
				XYZVarNums[j] = TecUtilVarGetNumByName(TmpStrs[i][j].c_str());
			VarsFound = (XYZVarNums[0] > 0 && XYZVarNums[1] > 0 && XYZVarNums[2] > 0);
		}
		if (!VarsFound){
			TecUtilDialogErrMsg("Failed to get x, y, z variable numbers.");
			TecUtilLockFinish(AddOnID);
			return;
		}
	}

// 	char tmpStr[100];
// 	sprintf(tmpStr, "Vol zone = %d, x/y/z var nums = %d/%d/%d", VolZoneNum, XYZVarNums[0], XYZVarNums[1], XYZVarNums[2]);
// 	TecUtilDialogMessageBox(tmpStr, MessageBoxType_Information);

	VolExtentIndexWeights_s VolInfo;
	VolInfo.AddOnID = AddOnID;
	GetVolInfo(VolZoneNum, XYZVarNums, FALSE, VolInfo);
	VolInfo.IsPeriodic = FALSE;


	vector<VolExtentIndexWeights_s> ThreadVolInfo(numCPU, VolInfo);

	Boolean_t IsOk = TRUE;


	vector<LgIndex_t> ZoneMaxIJK(3);

	vector<FieldDataPointer_c> VolReadPtrs(NumVars), ZoneWritePtrs(NumVars);// , XYZReadPtrs(3);
	TecUtilDataLoadBegin();

	for (int i = 1; i <= NumVars && IsOk; ++i){
		Boolean_t VarIsValid = TRUE;
		for (auto Beg = XYZVarNums.cbegin(), End = XYZVarNums.cend(); Beg != End && VarIsValid; Beg++)
			VarIsValid = (i != *Beg);
		if (VarIsValid)
			VolReadPtrs[i - 1].InitializeReadPtr(VolZoneNum, i);
	}

	for (auto & ChkPtr : VolReadPtrs){
		if (ChkPtr.IsReady()) {
			bool HasUniqueVals = false;

			for (int i = 0; i < MIN(4, ChkPtr.Size() - 1) && !HasUniqueVals; ++i) {
				for (int j = i + 1; j < MIN(5, ChkPtr.Size()) && !HasUniqueVals; ++j) {
					HasUniqueVals = (ChkPtr[j] != ChkPtr[i]);
				}
			}

			if (!HasUniqueVals) ChkPtr.Close();
		}
	}

	string StatusStr = "Mapping all variables to each zone";
	StatusLaunch(StatusStr.c_str(), AddOnID, TRUE);

	for (EntIndex_t ZoneNum = 1; ZoneNum <= NumZones && IsOk; ++ZoneNum){
		if (ZoneNum != VolZoneNum){
			IsOk = StatusUpdate(ZoneNum-1, NumZones, StatusStr, AddOnID);
			if (IsOk){
				// 				for (int i = 0; i < 3 && IsOk; ++i)
				// 					IsOk = XYZReadPtrs[i].InitializeReadPtr(ZoneNum, XYZVarNums[i]);
				for (int i = 1; i <= NumVars && IsOk; ++i)
					IsOk = ZoneWritePtrs[i - 1].InitializeWritePtr(ZoneNum, i);

				TecUtilZoneGetIJK(ZoneNum, &ZoneMaxIJK[0], &ZoneMaxIJK[1], &ZoneMaxIJK[2]);

				if (TecUtilZoneIsOrdered(ZoneNum)){
#ifndef DEBUG
#pragma omp parallel for
#endif
					for (int ii = 1; ii <= ZoneMaxIJK[0]; ++ii){
						int ThreadNum = omp_get_thread_num();
						int Index;
						vec3 Point;
						double ReadVal, WriteVal;

						for (int jj = 1; jj <= ZoneMaxIJK[1]; ++jj){
							for (int kk = 1; kk <= ZoneMaxIJK[2]; ++kk){
								Index = IndexFromIJK(ii, jj, kk, ZoneMaxIJK[0], ZoneMaxIJK[1]) - 1;

								for (int Dir = 0; Dir < 3; ++Dir)
									// 									Point[Dir] = XYZReadPtrs[Dir][Index];
									Point[Dir] = ZoneWritePtrs[Dir][Index];
								// 									Point[Dir] = TecUtilDataValueGetByZoneVar(ZoneNum, XYZVarNums[Dir], Index + 1);

								SetIndexAndWeightsForPoint(Point, ThreadVolInfo[ThreadNum]);

								for (int Var = 1; Var <= NumVars; ++Var){
									if (VolReadPtrs[Var - 1].IsReady()){
										ReadVal = ValByCurrentIndexAndWeightsFromRawPtr(ThreadVolInfo[ThreadNum], VolReadPtrs[Var - 1]);
										WriteVal = ZoneWritePtrs[Var - 1][Index];
										// 									if (WriteVal != 0)
										ZoneWritePtrs[Var - 1].Write(Index, ReadVal);
									}
								}
							}
						}
					}
				}
				else if (TecUtilZoneIsFiniteElement(ZoneNum)){
#ifndef DEBUG
#pragma omp parallel for
#endif
					for (int ii = 1; ii <= ZoneMaxIJK[0]; ++ii){
						int ThreadNum = omp_get_thread_num();
						vec3 Point;
						double ReadVal, WriteVal;

						for (int Dir = 0; Dir < 3; ++Dir)
							// 							Point[Dir] = XYZReadPtrs[Dir][ii - 1];
							Point[Dir] = ZoneWritePtrs[Dir][ii - 1];

						SetIndexAndWeightsForPoint(Point, ThreadVolInfo[ThreadNum]);

						for (int Var = 1; Var <= NumVars; ++Var){
							if (VolReadPtrs[Var - 1].IsReady()){
								ReadVal = ValByCurrentIndexAndWeightsFromRawPtr(ThreadVolInfo[ThreadNum], VolReadPtrs[Var - 1]);
								WriteVal = ZoneWritePtrs[Var - 1][ii - 1];
								// 							if (WriteVal != 0)
								ZoneWritePtrs[Var - 1].Write(ii - 1, ReadVal);
							}
						}
					}
				}
			}
		}
	}

	StatusDrop(AddOnID);

	TecUtilDataLoadEnd();

	TecUtilLockFinish(AddOnID);
}

void MapAllVarsToAllZones(int SourceZoneNum, vector<int> const & XYZVarNums, vector<int> const & CopyVarNums, vector<int> const & DestZoneNums, AddOn_pa const & AddOnID)
{
	REQUIRE(SourceZoneNum >= 1 && SourceZoneNum <= TecUtilDataSetGetNumZones());
	for (auto const & i : XYZVarNums){
		REQUIRE(i >= 1 && i <= TecUtilDataSetGetNumVars());
	}
	for (auto const & i : CopyVarNums) {
		REQUIRE(i >= 1 && i <= TecUtilDataSetGetNumVars());
	}
	for (auto const & i : DestZoneNums) {
		REQUIRE(i >= 1 && i <= TecUtilDataSetGetNumZones());
	}

	SYSTEM_INFO sysinfo;
	GetSystemInfo(&sysinfo);

	int numCPU = sysinfo.dwNumberOfProcessors;
	omp_set_num_threads(numCPU);

	TecUtilLockStart(AddOnID);

	// 	char tmpStr[100];
	// 	sprintf(tmpStr, "Vol zone = %d, x/y/z var nums = %d/%d/%d", VolZoneNum, XYZVarNums[0], XYZVarNums[1], XYZVarNums[2]);
	// 	TecUtilDialogMessageBox(tmpStr, MessageBoxType_Information);

	VolExtentIndexWeights_s VolInfo;
	VolInfo.AddOnID = AddOnID;
	int VolZoneNum = SourceZoneNum;
	GetVolInfo(VolZoneNum, XYZVarNums, FALSE, VolInfo);
	VolInfo.IsPeriodic = FALSE;


	vector<VolExtentIndexWeights_s> ThreadVolInfo(numCPU, VolInfo);

	Boolean_t IsOk = TRUE;


	vector<LgIndex_t> ZoneMaxIJK(3);

	vector<FieldDataPointer_c> VolReadPtrs;
	VolReadPtrs.reserve(CopyVarNums.size());
	TecUtilDataLoadBegin();

	for (int const & i : CopyVarNums) {
		VolReadPtrs.push_back(FieldDataPointer_c());
		VolReadPtrs.back().InitializeReadPtr(VolZoneNum, i);
	}

	for (auto & ChkPtr : VolReadPtrs) {
		if (ChkPtr.IsReady()) {
			bool HasUniqueVals = false;

			for (int i = 0; i < MIN(4, ChkPtr.Size() - 1) && !HasUniqueVals; ++i) {
				for (int j = i + 1; j < MIN(5, ChkPtr.Size()) && !HasUniqueVals; ++j) {
					HasUniqueVals = (ChkPtr[j] != ChkPtr[i]);
				}
			}

			if (!HasUniqueVals) ChkPtr.Close();
		}
	}

	string StatusStr = "Mapping all variables to each zone";
	StatusLaunch(StatusStr.c_str(), AddOnID, TRUE);

	int zoneIter = 0;
	for (int const ZoneNum : DestZoneNums) {
		zoneIter++;
		if (ZoneNum != VolZoneNum) {
			IsOk = StatusUpdate(zoneIter, DestZoneNums.size(), StatusStr, AddOnID);
			if (IsOk) {
				vector<FieldDataPointer_c> ZoneWritePtrs;
				ZoneWritePtrs.reserve(DestZoneNums.size());
				FieldVecPointer_c DestXYZ;
				DestXYZ.InitializeReadPtr(ZoneNum, XYZVarNums);

				for (int const & i : CopyVarNums) {
					ZoneWritePtrs.push_back(FieldDataPointer_c());
					IsOk = ZoneWritePtrs.back().InitializeWritePtr(ZoneNum, i);
				}

				TecUtilZoneGetIJK(ZoneNum, &ZoneMaxIJK[0], &ZoneMaxIJK[1], &ZoneMaxIJK[2]);

				if (TecUtilZoneIsOrdered(ZoneNum)) {
#ifndef DEBUG
#pragma omp parallel for
#endif
					for (int ii = 1; ii <= ZoneMaxIJK[0]; ++ii) {
						int ThreadNum = omp_get_thread_num();
						int Index;
						double ReadVal, WriteVal;

						for (int jj = 1; jj <= ZoneMaxIJK[1]; ++jj) {
							for (int kk = 1; kk <= ZoneMaxIJK[2]; ++kk) {
								Index = IndexFromIJK(ii, jj, kk, ZoneMaxIJK[0], ZoneMaxIJK[1]) - 1;

								SetIndexAndWeightsForPoint(DestXYZ[Index], ThreadVolInfo[ThreadNum]);

								for (int Var = 0; Var < ZoneWritePtrs.size(); ++Var) {
									ReadVal = ValByCurrentIndexAndWeightsFromRawPtr(ThreadVolInfo[ThreadNum], VolReadPtrs[Var]);
									WriteVal = ZoneWritePtrs[Var][Index];
									// 									if (WriteVal != 0)
									ZoneWritePtrs[Var].Write(Index, ReadVal);
								}
							}
						}
					}
				}
				else if (TecUtilZoneIsFiniteElement(ZoneNum)) {
#ifndef DEBUG
#pragma omp parallel for
#endif
					for (int ii = 1; ii <= ZoneMaxIJK[0]; ++ii) {
						int ThreadNum = omp_get_thread_num();
						double ReadVal, WriteVal;

						SetIndexAndWeightsForPoint(DestXYZ[ii-1], ThreadVolInfo[ThreadNum]);

						for (int Var = 0; Var < ZoneWritePtrs.size(); ++Var) {
							ReadVal = ValByCurrentIndexAndWeightsFromRawPtr(ThreadVolInfo[ThreadNum], VolReadPtrs[Var]);
							WriteVal = ZoneWritePtrs[Var][ii - 1];
							// 							if (WriteVal != 0)
							ZoneWritePtrs[Var].Write(ii - 1, ReadVal);
						}
					}
				}
			}
		}
	}

	StatusDrop(AddOnID);

	TecUtilDataLoadEnd();

	TecUtilLockFinish(AddOnID);
}


void CalcVars(CalcVarsOptions_s & Opt)
{
	SYSTEM_INFO sysinfo;
	GetSystemInfo(&sysinfo);
	int numCPU = sysinfo.dwNumberOfProcessors;
	omp_set_num_threads(numCPU);


	TecUtilLockStart(Opt.AddOnID);

	/*
	*	Get dataset info
	*/
	EntIndex_t NumZones, NumVars;
	Boolean_t IsOk = TecUtilDataSetGetInfo(nullptr, &NumZones, &NumVars);
	if (!IsOk){
		TecUtilDialogErrMsg("Failed to get dataset information");
		TecUtilLockFinish(Opt.AddOnID);
		return;
	}
	vector<ValueLocation_e> DataLoc(NumZones, ValueLocation_Nodal);
	vector<FieldDataType_e> DataType;
	DataType.resize(NumZones, FieldDataType_Float);


	/*
	*	Get XYZ variable numbers
	*/
	vector<EntIndex_t> XYZVarNums(3, -1);
	TecUtilAxisGetVarAssignments(&XYZVarNums[0], &XYZVarNums[1], &XYZVarNums[2]);
	if (XYZVarNums[0] <= 0 || XYZVarNums[1] <= 0 || XYZVarNums[2] <= 0){
		vector<vector<string> > TmpStrs = {
			{ "X", "Y", "Z" }, { "I", "J", "K" }
		};
		IsOk = FALSE;
		for (int i = 0; i < 2 && !IsOk; ++i){
			for (int j = 0; j < 3; ++j)
				XYZVarNums[j] = TecUtilVarGetNumByName(TmpStrs[i][j].c_str());
			IsOk = (XYZVarNums[0] > 0 && XYZVarNums[1] > 0 && XYZVarNums[2] > 0);
		}
		if (!IsOk){
			TecUtilDialogErrMsg("Failed to get x, y, z variable numbers.");
			TecUtilLockFinish(Opt.AddOnID);
			return;
		}
	}


	/*
	*	Get full volume zone number, meaning the largest of the volume zones, assumed to define the boundaries of the system.
	*/

	EntIndex_t VolZoneNum = -1;
	vec3 MaxXYZ = ones(3) * (-DBL_MAX), MinXYZ = ones(3) * DBL_MAX;
	for (int i = 1; i <= NumZones; ++i){
		if (TecUtilZoneIsOrdered(i)){
			int IJK[3];
			TecUtilZoneGetIJK(i, &IJK[0], &IJK[1], &IJK[2]);
			if (IJK[0] > 1 && IJK[1] > 1 && IJK[2] > 1){
				vector<vec3> TmpMinMaxXYZ = ZoneXYZVarGetMinMax_Ordered3DZone(XYZVarNums, i);
				if (sum(TmpMinMaxXYZ[1] > MaxXYZ) == 3 && sum(TmpMinMaxXYZ[0] < MinXYZ) == 3){
					VolZoneNum = i;
					MaxXYZ = TmpMinMaxXYZ[1];
					MinXYZ = TmpMinMaxXYZ[0];
				}
			}
		}
	}

	IsOk = (VolZoneNum > 0);
	if (!IsOk){
		TecUtilDialogErrMsg("Failed to get volume zone");
		TecUtilLockFinish(Opt.AddOnID);
		return;
	}

	/*
	*	Get system information
	*/
	VolExtentIndexWeights_s VolInfo;
	GetVolInfo(VolZoneNum, XYZVarNums, Opt.IsPeriodic, VolInfo);
	VolInfo.AddOnID = Opt.AddOnID;
	VolInfo.IsPeriodic = Opt.IsPeriodic;

	vec3 OneOverDelXYZx2 = VolInfo.PointSpacingV123 * 2,
		OneOverDelXYZx12 = VolInfo.PointSpacingV123 * 12;

	int HessIndices[3][3] = {
		{ 0, 1, 2 },
		{ 1, 3, 4 },
		{ 2, 4, 5 }
	};

	/*
	*	All necessary variable names
	*/
	vector<string> GradVarNames = CSMVarName.DensGradVec,
		HessVarNames = CSMVarName.DensHessTensor,
		EigSysVarNames = CSMVarName.EigSys,
		EigVecDotGradVarNames = CSMVarName.GradDotEigVecs,
		EberlyFuncVarNames = CSMVarName.EberlyFunc;
	string LapVarName = CSMVarName.DensLap,
		GaussCurvatureVarName = CSMVarName.GaussCurvature,
		GradMagVarName = CSMVarName.DensGradMag,
		EigenRankVarName = "Rank of Eigenvalues";
	for (int i = 0; i < 2; ++i){
		if (Opt.EberlyUseCutoff[i]){
			EberlyFuncVarNames[i] += string(" cutoff at " + to_string(Opt.EberlyCutoff[i]));
		}
	}

	/*
	*	For keeping track of which variables are already present in dataset
	*/
	Boolean_t HasGradMag = FALSE,
		HasEigSys = FALSE,
		HasLap = FALSE,
		HasGaussCurvature = FALSE,
		HasEigVecDotGrad = FALSE,
		HasEberlyFunc[2] = { FALSE },
		HasEigenRank = FALSE;

	EntIndex_t GradMagVarNum = -1,
		LapVarNum = -1,
		GaussCurvatureVarNum = -1,
		EigenRankVarNum = -1;

	vector<EntIndex_t> EigSysVarNums(12, -1),
		EigVecDotGradVarNums(3, -1),
		EberlyFuncVarNums(2, -1);

	/*
	*	Raw pointers for all variables
	*/
	FieldDataPointer_c RhoPtr,
		GradMagPtr,
		LapPtr,
		GaussCurvaturePtr,
		EigenRankPtr;;

	vector<FieldDataPointer_c> GradPtrs(3),
		HessPtrs(6),
		EigSysPtrs(12),
		EigVecDotGradPtrs(3),
		EberlyFuncPtrs(2);

	/*
	*	Then make multiple copies of some workspace variables so that
	*	each thread gets its own copy
	*/
	vector<vector<double> > ThreadVals(numCPU, vector<double>(5)),
		ThreadDotPtds(numCPU, vector<double>(3));
	vector<mat33> ThreadHessVec(numCPU);
	vector<vector<vec3> > ThreadPoints(numCPU, vector<vec3>(5));
	vector<vector<int> > ThreadDirInd(numCPU, vector<int>(5));
	vector<vec3> ThreadGrad(numCPU),
		ThreadEigVals(numCPU),
		ThreadCurPoint(numCPU);
	vector<mat33> ThreadHess(numCPU),
		ThreadEigVecs(numCPU);
	vector<double> TmpVal(numCPU);

	vector<MultiRootParams_s> Params(numCPU);
	vector<VolExtentIndexWeights_s> VI(numCPU, VolInfo);
	mat33 I = eye<mat>(3, 3);
	for (int i = 0; i < numCPU; ++i){
		Params[i].BasisVectors = &I;
		Params[i].CalcType = GPType_Invalid;
		Params[i].IsPeriodic = Opt.IsPeriodic;
		Params[i].VolInfo = &VI[i];
	}

	Boolean_t ShowRequiredVarsMessage = TRUE;

	/*
	*	Loop over each variable that needs to be calculated
	*/
	int IterNum = 0;
	for (auto CalcVar = Opt.CalcVarList.cbegin(), End = Opt.CalcVarList.cend(); CalcVar != End && IsOk; CalcVar++){
		/*
		*	First, see if variables are already present in dataset.
		*  If present, they'll be overwritten if the user requested
		*  to calculate them, or they'll be used in the calculation
		*  of other variables if necessary.
		*/
		if (*CalcVar >= CalcGradientVectors){
			if (!Opt.HasGrad && IterNum > 0){
				Opt.HasGrad = TRUE;
				for (int i = 0; i < 3 && Opt.HasGrad; ++i){
					Opt.GradVarNums[i] = VarNumByName(GradVarNames[i]);
					Opt.HasGrad = (Opt.GradVarNums[i] > 0);
				}
			}
		}
		if (*CalcVar == CalcGradientMagnitude){
			GradMagVarNum = VarNumByName(GradMagVarName);
			HasGradMag = (GradMagVarNum > 0);
		}
		if (*CalcVar >= CalcHessian){
			if (!Opt.HasHess && IterNum > 0){
				Opt.HasHess = TRUE;
				for (int i = 0; i < 3 && Opt.HasHess; ++i){
					Opt.HessVarNums[i] = VarNumByName(HessVarNames[i]);
					Opt.HasHess = (Opt.HessVarNums[i] > 0);
				}
			}
		}
		if (*CalcVar >= CalcEigenSystem){
			if (!HasEigSys){
				HasEigSys = TRUE;
				for (int i = 0; i < 12 && HasEigSys; ++i){
					EigSysVarNums[i] = VarNumByName(EigSysVarNames[i]);
					HasEigSys = (EigSysVarNums[i] > 0);
				}
			}
		}
		if (*CalcVar == CalcLaplacian){
			if (!HasLap){
				LapVarNum = VarNumByName(LapVarName);
				HasLap = (LapVarNum > 0);
			}
		}
		if (*CalcVar == CalcGaussianCurvature){
			if (!HasGaussCurvature){
				GaussCurvatureVarNum = VarNumByName(GaussCurvatureVarName);
				HasGaussCurvature = (GaussCurvatureVarNum > 0);
			}
		}
		if (*CalcVar == CalcEigenVectorsDotGradient){
			if (!HasEigVecDotGrad){
				HasEigVecDotGrad = TRUE;
				for (int i = 0; i < 3 && HasEigVecDotGrad; ++i){
					EigVecDotGradVarNums[i] = VarNumByName(EigVecDotGradVarNames[i]);
					HasEigVecDotGrad = (EigVecDotGradVarNums[i] > 0);
				}
			}
		}
		if (*CalcVar == CalcEberly1Ridge){
			if (!HasEberlyFunc[0]){
				EberlyFuncVarNums[0] = VarNumByName(EberlyFuncVarNames[0]);
				HasEberlyFunc[0] = (EberlyFuncVarNums[0] > 0);
			}
		}
		if (*CalcVar == CalcEberly2Ridge){
			if (!HasEberlyFunc[1]){
				EberlyFuncVarNums[1] = VarNumByName(EberlyFuncVarNames[1]);
				HasEberlyFunc[1] = (EberlyFuncVarNums[1] > 0);
			}
		}
		if (*CalcVar == CalcEigenRank){
			if (!HasEigenRank){
				EigenRankVarNum = VarNumByName(EigenRankVarName);
				HasEigenRank = (EigenRankVarNum > 0);
			}
		}

		Boolean_t CheckForRequiredVars = TRUE;
		if (!Opt.CalcForAllZones){
			int IJK[3];
			ZoneType_e ZoneType = TecUtilZoneGetType(Opt.CalcZoneNum);
			TecUtilZoneGetIJK(Opt.CalcZoneNum, &IJK[0], &IJK[1], &IJK[2]);
			if (TecUtilZoneIsFiniteElement(Opt.CalcZoneNum)
				|| (ZoneType == ZoneType_Ordered && IJK[2] <= 1)){
				CheckForRequiredVars = TRUE;
			}
			else
				CheckForRequiredVars = FALSE;
		}

		Boolean_t HasRequiredVars = TRUE;
		if (CheckForRequiredVars){
			switch (*CalcVar){
				case CalcGradientVectors:
					HasRequiredVars = Opt.HasGrad;
					break;
				case CalcGradientMagnitude:
					HasRequiredVars = Opt.HasGrad;
					break;
				case CalcHessian:
					HasRequiredVars = Opt.HasHess;
					break;
				case CalcEigenSystem:
					HasRequiredVars = Opt.HasHess;
					break;
				case CalcLaplacian:
					HasRequiredVars = (Opt.HasHess || HasEigSys);
					break;
				case CalcGaussianCurvature:
					HasRequiredVars = (Opt.HasHess || HasEigSys);
					break;
				case CalcEigenVectorsDotGradient:
					HasRequiredVars = (Opt.HasHess || HasEigSys) && Opt.HasGrad;
					break;
				case CalcEberly1Ridge:
					HasRequiredVars = (Opt.HasHess || HasEigSys) && Opt.HasGrad;
					break;
				case CalcEberly2Ridge:
					HasRequiredVars = (Opt.HasHess || HasEigSys) && Opt.HasGrad;
					break;
				case CalcEigenRank:
					HasRequiredVars = (Opt.HasHess || HasEigSys);
					break;
				default:
					break;
			}

			if (!HasRequiredVars && ShowRequiredVarsMessage){
				TecUtilDialogMessageBox("Required variables are not present for non-volume zones, so calculation will only be for volume zones.", MessageBoxType_Information);
				ShowRequiredVarsMessage = FALSE;
				// 				TecUtilDialogErrMsg("The selected variable(s) cannot be calculated for ridge zones because required source variable(s) are missing");
				// 				TecUtilLockFinish(Opt.AddOnID);
				// 				return;
			}
		}

		/*
		*	Here we create variables if they're to be calculated and don't
		*	already exist in the dataset.
		*/
		switch (*CalcVar){
			case CalcGradientVectors:
				for (int i = 0; i < 3; ++i){
					if (Opt.GradVarNums[i] <= 0){
						ArgList_pa Args = TecUtilArgListAlloc();
						TecUtilArgListAppendString(Args, SV_NAME, GradVarNames[i].c_str());
						TecUtilArgListAppendArray(Args, SV_VARDATATYPE, DataType.data());
						TecUtilArgListAppendArray(Args, SV_VALUELOCATION, DataLoc.data());
						TecUtilArgListAppendInt(Args, SV_SHAREVARWITHALLZONES, FALSE);

						if (TecUtilDataSetAddVarX(Args)){
							Set_pa NewVar = TecUtilSetAlloc(FALSE);
							Opt.GradVarNums[i] = TecUtilDataSetGetNumVars();
							TecUtilSetAddMember(NewVar, Opt.GradVarNums[i], FALSE);
							TecUtilStateChanged(StateChange_VarsAdded, reinterpret_cast<ArbParam_t>(NewVar));
							TecUtilSetDealloc(&NewVar);
						}
						TecUtilArgListDealloc(&Args);
					}
				}
				break;
			case CalcGradientMagnitude:
				if (GradMagVarNum <= 0){
					ArgList_pa Args = TecUtilArgListAlloc();
					TecUtilArgListAppendString(Args, SV_NAME, GradMagVarName.c_str());
					TecUtilArgListAppendArray(Args, SV_VARDATATYPE, DataType.data());
					TecUtilArgListAppendArray(Args, SV_VALUELOCATION, DataLoc.data());
					TecUtilArgListAppendInt(Args, SV_SHAREVARWITHALLZONES, FALSE);

					if (TecUtilDataSetAddVarX(Args)){
						Set_pa NewVar = TecUtilSetAlloc(FALSE);
						GradMagVarNum = TecUtilDataSetGetNumVars();
						TecUtilSetAddMember(NewVar, GradMagVarNum, FALSE);
						TecUtilStateChanged(StateChange_VarsAdded, reinterpret_cast<ArbParam_t>(NewVar));
						TecUtilSetDealloc(&NewVar);
					}
					TecUtilArgListDealloc(&Args);
				}
				break;
			case CalcHessian:
				for (int i = 0; i < 6; ++i){
					if (Opt.HessVarNums[i] <= 0){
						ArgList_pa Args = TecUtilArgListAlloc();
						TecUtilArgListAppendString(Args, SV_NAME, HessVarNames[i].c_str());
						TecUtilArgListAppendArray(Args, SV_VARDATATYPE, DataType.data());
						TecUtilArgListAppendArray(Args, SV_VALUELOCATION, DataLoc.data());
						TecUtilArgListAppendInt(Args, SV_SHAREVARWITHALLZONES, FALSE);

						if (TecUtilDataSetAddVarX(Args)){
							Set_pa NewVar = TecUtilSetAlloc(FALSE);
							Opt.HessVarNums[i] = TecUtilDataSetGetNumVars();
							TecUtilSetAddMember(NewVar, Opt.HessVarNums[i], FALSE);
							TecUtilStateChanged(StateChange_VarsAdded, reinterpret_cast<ArbParam_t>(NewVar));
							TecUtilSetDealloc(&NewVar);
						}
						TecUtilArgListDealloc(&Args);
					}
				}
				break;
			case CalcEigenSystem:
				for (int i = 0; i < 12; ++i){
					if (EigSysVarNums[i] <= 0){
						ArgList_pa Args = TecUtilArgListAlloc();
						TecUtilArgListAppendString(Args, SV_NAME, EigSysVarNames[i].c_str());
						TecUtilArgListAppendArray(Args, SV_VARDATATYPE, DataType.data());
						TecUtilArgListAppendArray(Args, SV_VALUELOCATION, DataLoc.data());
						TecUtilArgListAppendInt(Args, SV_SHAREVARWITHALLZONES, FALSE);

						if (TecUtilDataSetAddVarX(Args)){
							Set_pa NewVar = TecUtilSetAlloc(FALSE);
							EigSysVarNums[i] = TecUtilDataSetGetNumVars();
							TecUtilSetAddMember(NewVar, EigSysVarNums[i], FALSE);
							TecUtilStateChanged(StateChange_VarsAdded, reinterpret_cast<ArbParam_t>(NewVar));
							TecUtilSetDealloc(&NewVar);
						}
						TecUtilArgListDealloc(&Args);
					}
				}
				break;
			case CalcLaplacian:
				if (LapVarNum <= 0){
					ArgList_pa Args = TecUtilArgListAlloc();
					TecUtilArgListAppendString(Args, SV_NAME, LapVarName.c_str());
					TecUtilArgListAppendArray(Args, SV_VARDATATYPE, DataType.data());
					TecUtilArgListAppendArray(Args, SV_VALUELOCATION, DataLoc.data());
					TecUtilArgListAppendInt(Args, SV_SHAREVARWITHALLZONES, FALSE);

					if (TecUtilDataSetAddVarX(Args)){
						Set_pa NewVar = TecUtilSetAlloc(FALSE);
						LapVarNum = TecUtilDataSetGetNumVars();
						TecUtilSetAddMember(NewVar, LapVarNum, FALSE);
						TecUtilStateChanged(StateChange_VarsAdded, reinterpret_cast<ArbParam_t>(NewVar));
						TecUtilSetDealloc(&NewVar);
					}
					TecUtilArgListDealloc(&Args);
				}
				break;
			case CalcGaussianCurvature:
				if (GaussCurvatureVarNum <= 0){
					ArgList_pa Args = TecUtilArgListAlloc();
					TecUtilArgListAppendString(Args, SV_NAME, GaussCurvatureVarName.c_str());
					TecUtilArgListAppendArray(Args, SV_VARDATATYPE, DataType.data());
					TecUtilArgListAppendArray(Args, SV_VALUELOCATION, DataLoc.data());
					TecUtilArgListAppendInt(Args, SV_SHAREVARWITHALLZONES, FALSE);

					if (TecUtilDataSetAddVarX(Args)){
						Set_pa NewVar = TecUtilSetAlloc(FALSE);
						GaussCurvatureVarNum = TecUtilDataSetGetNumVars();
						TecUtilSetAddMember(NewVar, GaussCurvatureVarNum, FALSE);
						TecUtilStateChanged(StateChange_VarsAdded, reinterpret_cast<ArbParam_t>(NewVar));
						TecUtilSetDealloc(&NewVar);
					}
					TecUtilArgListDealloc(&Args);
				}
				break;
			case CalcEigenVectorsDotGradient:
				for (int i = 0; i < 3; ++i){
					if (EigVecDotGradVarNums[i] <= 0){
						ArgList_pa Args = TecUtilArgListAlloc();
						TecUtilArgListAppendString(Args, SV_NAME, EigVecDotGradVarNames[i].c_str());
						TecUtilArgListAppendArray(Args, SV_VARDATATYPE, DataType.data());
						TecUtilArgListAppendArray(Args, SV_VALUELOCATION, DataLoc.data());
						TecUtilArgListAppendInt(Args, SV_SHAREVARWITHALLZONES, FALSE);

						if (TecUtilDataSetAddVarX(Args)){
							Set_pa NewVar = TecUtilSetAlloc(FALSE);
							EigVecDotGradVarNums[i] = TecUtilDataSetGetNumVars();
							TecUtilSetAddMember(NewVar, EigVecDotGradVarNums[i], FALSE);
							TecUtilStateChanged(StateChange_VarsAdded, reinterpret_cast<ArbParam_t>(NewVar));
							TecUtilSetDealloc(&NewVar);
						}
						TecUtilArgListDealloc(&Args);
					}
				}
				break;
			case CalcEberly1Ridge:
				if (EberlyFuncVarNums[0] <= 0){
					ArgList_pa Args = TecUtilArgListAlloc();
					TecUtilArgListAppendString(Args, SV_NAME, EberlyFuncVarNames[0].c_str());
					TecUtilArgListAppendArray(Args, SV_VARDATATYPE, DataType.data());
					TecUtilArgListAppendArray(Args, SV_VALUELOCATION, DataLoc.data());
					TecUtilArgListAppendInt(Args, SV_SHAREVARWITHALLZONES, FALSE);

					if (TecUtilDataSetAddVarX(Args)){
						Set_pa NewVar = TecUtilSetAlloc(FALSE);
						EberlyFuncVarNums[0] = TecUtilDataSetGetNumVars();
						TecUtilSetAddMember(NewVar, EberlyFuncVarNums[0], FALSE);
						TecUtilStateChanged(StateChange_VarsAdded, reinterpret_cast<ArbParam_t>(NewVar));
						TecUtilSetDealloc(&NewVar);
					}
					TecUtilArgListDealloc(&Args);
				}
				break;
			case CalcEberly2Ridge:
				if (EberlyFuncVarNums[1] <= 0){
					ArgList_pa Args = TecUtilArgListAlloc();
					TecUtilArgListAppendString(Args, SV_NAME, EberlyFuncVarNames[1].c_str());
					TecUtilArgListAppendArray(Args, SV_VARDATATYPE, DataType.data());
					TecUtilArgListAppendArray(Args, SV_VALUELOCATION, DataLoc.data());
					TecUtilArgListAppendInt(Args, SV_SHAREVARWITHALLZONES, FALSE);

					if (TecUtilDataSetAddVarX(Args)){
						Set_pa NewVar = TecUtilSetAlloc(FALSE);
						EberlyFuncVarNums[1] = TecUtilDataSetGetNumVars();
						TecUtilSetAddMember(NewVar, EberlyFuncVarNums[1], FALSE);
						TecUtilStateChanged(StateChange_VarsAdded, reinterpret_cast<ArbParam_t>(NewVar));
						TecUtilSetDealloc(&NewVar);
					}
					TecUtilArgListDealloc(&Args);
				}
				break;
			case CalcEigenRank:
				if (EigenRankVarNum <= 0){
					ArgList_pa Args = TecUtilArgListAlloc();
					TecUtilArgListAppendString(Args, SV_NAME, EigenRankVarName.c_str());
					TecUtilArgListAppendArray(Args, SV_VARDATATYPE, DataType.data());
					TecUtilArgListAppendArray(Args, SV_VALUELOCATION, DataLoc.data());
					TecUtilArgListAppendInt(Args, SV_SHAREVARWITHALLZONES, FALSE);

					if (TecUtilDataSetAddVarX(Args)){
						Set_pa NewVar = TecUtilSetAlloc(FALSE);
						EigenRankVarNum = TecUtilDataSetGetNumVars();
						TecUtilSetAddMember(NewVar, EigenRankVarNum, FALSE);
						TecUtilStateChanged(StateChange_VarsAdded, reinterpret_cast<ArbParam_t>(NewVar));
						TecUtilSetDealloc(&NewVar);
					}
					TecUtilArgListDealloc(&Args);
				}
				break;
			default:
				TecUtilDialogErrMsg("Invalid variable used");
				IsOk = FALSE;
				break;
		}

		/*
		*	Now actually start calculating variables
		*/

		vector<int> MaxIJK(3);
		vec3 BlankPoint;

		string StatusStr = string("Calculating ") + CalcVarTypeNames[*CalcVar];
		StatusLaunch(StatusStr.c_str(), VolInfo.AddOnID, TRUE);

		/*
		*	Loop over every zone and calculate the necessary variables.
		*	If the user only wants the volume zone, then this loop will terminate after
		*	the first iteration.
		*
		*	There are different functions used depending on if the calculation is for
		*	the full volume zone or not. The node-based functions are faster and
		*	have less numerical error.
		*/
		for (EntIndex_t ZoneNum = 1; ZoneNum <= NumZones && IsOk; ++ZoneNum){

			Boolean_t IsVolZone = TecUtilZoneIsOrdered(ZoneNum);
			if (IsOk && (Opt.CalcForAllZones || ZoneNum == Opt.CalcZoneNum)){
				vector<int> IJK(3);
				if (IsVolZone){
					TecUtilZoneGetIJK(ZoneNum, &IJK[0], &IJK[1], &IJK[2]);
					IsVolZone = (IJK[2] > 1);
				}
				if (IsVolZone){
					VolInfo.MaxIJK = IJK;
					if (ZoneNum == VolZoneNum){
						GetVolInfo(ZoneNum, XYZVarNums, Opt.IsPeriodic, VolInfo);
					}
					else{
						GetVolInfo(ZoneNum, XYZVarNums, FALSE, VolInfo);
					}

					OneOverDelXYZx2 = 1./(VolInfo.PointSpacingV123 * 2);
					OneOverDelXYZx12 = 1./(VolInfo.PointSpacingV123 * 12);

					VI = vector<VolExtentIndexWeights_s>(numCPU, VolInfo);
					for (int i = 0; i < numCPU; ++i){
						Params[i].VolInfo = &VI[i];
						
					}
				}
			}

			if (IsOk && (Opt.CalcForAllZones || ZoneNum == Opt.CalcZoneNum))
			{
				/*
				*	Check if variables are present for the current zone.
				*	We already know if they are or aren't in the dataset,
				*	but we need to check that they have non-zero values
				*	for the current zone, otherwise we'll be doing expensive
				*	calculations using all zeros!
				*/
				if (*CalcVar >= CalcGradientVectors){
					if (Opt.HasGrad){
						Opt.HasGrad = FALSE;
						for (int i = 0; i < 3 && !Opt.HasGrad; ++i){
							FieldDataPointer_c TmpPtr;
							if (TmpPtr.InitializeReadPtr(ZoneNum, Opt.GradVarNums[i])){
								int jMax = TmpPtr.Size();
#pragma omp parallel for
								for (int j = 0; j < jMax; ++j){
									if (!Opt.HasGrad && TmpPtr[j] != 0.0)
										Opt.HasGrad = TRUE;
								}
							}
						}
					}
				}
				if (*CalcVar == CalcGradientMagnitude){
					if (HasGradMag){
						HasGradMag = FALSE;
						FieldDataPointer_c TmpPtr;
						if (TmpPtr.InitializeReadPtr(ZoneNum, GradMagVarNum)){
							int jMax = TmpPtr.Size();
#pragma omp parallel for
							for (int j = 0; j < jMax; ++j){
								if (!HasGradMag && TmpPtr[j] != 0.0)
									HasGradMag = TRUE;
							}
						}
					}
					GradMagVarNum = VarNumByName(GradMagVarName);
					HasGradMag = (GradMagVarNum > 0);
				}
				if (*CalcVar >= CalcHessian){
					if (Opt.HasHess){
						Opt.HasHess = FALSE;
						for (int i = 0; i < 6 && !Opt.HasHess; ++i){
							FieldDataPointer_c TmpPtr;
							if (TmpPtr.InitializeReadPtr(ZoneNum, Opt.HessVarNums[i])){
								int jMax = TmpPtr.Size();
#pragma omp parallel for
								for (int j = 0; j < jMax; ++j){
									if (!Opt.HasHess && TmpPtr[j] != 0.0)
										Opt.HasHess = TRUE;
								}
							}
						}
					}
				}
				if (*CalcVar >= CalcEigenSystem){
					if (HasEigSys){
						HasEigSys = FALSE;
						for (int i = 0; i < 12 && !HasEigSys; ++i){
							FieldDataPointer_c TmpPtr;
							if (TmpPtr.InitializeReadPtr(ZoneNum, EigSysVarNums[i])){
								int jMax = TmpPtr.Size();
#pragma omp parallel for
								for (int j = 0; j < jMax; ++j){
									if (!HasEigSys && TmpPtr[j] != 0.0)
										HasEigSys = TRUE;
								}
							}
						}
					}
				}
				if (*CalcVar == CalcLaplacian){
					if (HasLap){
						HasLap = FALSE;
						FieldDataPointer_c TmpPtr;
						if (TmpPtr.InitializeReadPtr(ZoneNum, LapVarNum)){
							int jMax = TmpPtr.Size();
#pragma omp parallel for
							for (int j = 0; j < jMax; ++j){
								if (!HasLap && TmpPtr[j] != 0.0)
									HasLap = TRUE;
							}
						}
					}
				}
				if (*CalcVar == CalcGaussianCurvature){
					if (HasGaussCurvature){
						HasGaussCurvature = FALSE;
						FieldDataPointer_c TmpPtr;
						if (TmpPtr.InitializeReadPtr(ZoneNum, GaussCurvatureVarNum)){
							int jMax = TmpPtr.Size();
#pragma omp parallel for
							for (int j = 0; j < jMax; ++j){
								if (!HasGaussCurvature && TmpPtr[j] != 0.0)
									HasGaussCurvature = TRUE;
							}
						}
					}
				}
				if (*CalcVar == CalcEigenVectorsDotGradient){
					if (HasEigVecDotGrad){
						HasEigVecDotGrad = FALSE;
						for (int i = 0; i < 3 && !HasEigVecDotGrad; ++i){
							FieldDataPointer_c TmpPtr;
							if (TmpPtr.InitializeReadPtr(ZoneNum, EigVecDotGradVarNums[i])){
								int jMax = TmpPtr.Size();
#pragma omp parallel for
								for (int j = 0; j < jMax; ++j){
									if (!HasEigVecDotGrad && TmpPtr[j] != 0.0)
										HasEigVecDotGrad = TRUE;
								}
							}
						}
					}
				}
				if (*CalcVar == CalcEberly1Ridge){
					if (HasEberlyFunc[0]){
						HasEberlyFunc[0] = FALSE;
						FieldDataPointer_c TmpPtr;
						if (TmpPtr.InitializeReadPtr(ZoneNum, EberlyFuncVarNums[0])){
							int jMax = TmpPtr.Size();
#pragma omp parallel for
							for (int j = 0; j < jMax; ++j){
								if (!HasEberlyFunc[0] && TmpPtr[j] != 0.0)
									HasEberlyFunc[0] = TRUE;
							}
						}
					}
				}
				if (*CalcVar == CalcEberly2Ridge){
					if (HasEberlyFunc[1]){
						HasEberlyFunc[1] = FALSE;
						FieldDataPointer_c TmpPtr;
						if (TmpPtr.InitializeReadPtr(ZoneNum, EberlyFuncVarNums[1])){
							int jMax = TmpPtr.Size();
#pragma omp parallel for
							for (int j = 0; j < jMax; ++j){
								if (!HasEberlyFunc[1] && TmpPtr[j] != 0.0)
									HasEberlyFunc[1] = TRUE;
							}
						}
					}
				}
				if (*CalcVar == CalcEigenRank){
					if (HasEigenRank){
						HasEigenRank = FALSE;
						FieldDataPointer_c TmpPtr;
						if (TmpPtr.InitializeReadPtr(ZoneNum, EigenRankVarNum)){
							int jMax = TmpPtr.Size();
#pragma omp parallel for
							for (int j = 0; j < jMax; ++j){
								if (!HasEigenRank && TmpPtr[j] != 0.0)
									HasEigenRank = TRUE;
							}
						}
					}
				}

				if (CheckForRequiredVars){
					switch (*CalcVar){
						case CalcGradientVectors:
							HasRequiredVars = Opt.HasGrad;
							break;
						case CalcGradientMagnitude:
							HasRequiredVars = Opt.HasGrad;
							break;
						case CalcHessian:
							HasRequiredVars = Opt.HasHess;
							break;
						case CalcEigenSystem:
							HasRequiredVars = Opt.HasHess;
							break;
						case CalcLaplacian:
							HasRequiredVars = (Opt.HasHess || HasEigSys);
							break;
						case CalcGaussianCurvature:
							HasRequiredVars = (Opt.HasHess || HasEigSys);
							break;
						case CalcEigenVectorsDotGradient:
							HasRequiredVars = (Opt.HasHess || HasEigSys) && Opt.HasGrad;
							break;
						case CalcEberly1Ridge:
							HasRequiredVars = (Opt.HasHess || HasEigSys) && Opt.HasGrad;
							break;
						case CalcEberly2Ridge:
							HasRequiredVars = (Opt.HasHess || HasEigSys) && Opt.HasGrad;
							break;
						case CalcEigenRank:
							HasRequiredVars = (Opt.HasHess || HasEigSys);
							break;
						default:
							break;
					}

					if (!HasRequiredVars && ShowRequiredVarsMessage){
						TecUtilDialogMessageBox("Required variables are not present for non-volume zones, so calculation will only be for volume zones.", MessageBoxType_Information);
						ShowRequiredVarsMessage = FALSE;
						// 				TecUtilDialogErrMsg("The selected variable(s) cannot be calculated for ridge zones because required source variable(s) are missing");
						// 				TecUtilLockFinish(Opt.AddOnID);
						// 				return;
					}
				}


				TecUtilDataLoadBegin();


				if (Opt.CalcForAllZones)
					IsOk = StatusUpdate(ZoneNum-1, NumZones, StatusStr, Opt.AddOnID, high_resolution_clock::time_point(), false);
				/*
				*	Get all the read pointers
				*/
				// 				if (*CalcVar > CalcEigVecDotGrad && HasEigVecDotGrad){
				// 					for (int i = 0; i < 3 && IsOk; ++i)
				// 						IsOk = EigVecDotGradPtrs[i].InitializeReadPtr(ZoneNum, EigVecDotGradVarNums[i]);
				// 				}
				// 				else 
				if (*CalcVar > CalcEigenSystem && HasEigSys){
					for (int i = 0; i < 12 && IsOk; ++i)
						IsOk = EigSysPtrs[i].InitializeReadPtr(ZoneNum, EigSysVarNums[i]);
				}
				else if (*CalcVar > CalcHessian && Opt.HasHess){
					for (int i = 0; i < 6 && IsOk; ++i)
						IsOk = HessPtrs[i].InitializeReadPtr(ZoneNum, Opt.HessVarNums[i]);
				}
				else if (*CalcVar > CalcGradientVectors && Opt.HasGrad){
					for (int i = 0; i < 3 && IsOk; ++i)
						IsOk = GradPtrs[i].InitializeReadPtr(ZoneNum, Opt.GradVarNums[i]);
				}
				else if (IsOk && *CalcVar > CalcInvalidVar)
					IsOk = RhoPtr.InitializeReadPtr(ZoneNum, Opt.RhoVarNum);

				if (*CalcVar > CalcEigenVectorsDotGradient){
					if (IsOk)
						IsOk = RhoPtr.InitializeReadPtr(ZoneNum, Opt.RhoVarNum);
				}

				if (*CalcVar >= CalcEigenVectorsDotGradient){
					if (Opt.HasGrad)
					for (int i = 0; i < 3 && IsOk; ++i)
						IsOk = GradPtrs[i].InitializeReadPtr(ZoneNum, Opt.GradVarNums[i]);
					else if (IsOk)
						IsOk = RhoPtr.InitializeReadPtr(ZoneNum, Opt.RhoVarNum);
				}

				/*
				*	Get the write pointers
				*/
				if (IsOk){
					switch (*CalcVar){
						case CalcGradientVectors:
							for (int i = 0; i < 3 && IsOk; ++i)
								IsOk = GradPtrs[i].InitializeWritePtr(ZoneNum, Opt.GradVarNums[i]);
							break;
						case CalcGradientMagnitude:
							if (IsOk)
								IsOk = GradMagPtr.InitializeWritePtr(ZoneNum, GradMagVarNum);
							break;
						case CalcHessian:
							for (int i = 0; i < 6 && IsOk; ++i)
								IsOk = HessPtrs[i].InitializeWritePtr(ZoneNum, Opt.HessVarNums[i]);
							break;
						case CalcEigenSystem:
							for (int i = 0; i < 12 && IsOk; ++i)
								IsOk = EigSysPtrs[i].InitializeWritePtr(ZoneNum, EigSysVarNums[i]);
							break;
						case CalcLaplacian:
							if (IsOk)
								IsOk = LapPtr.InitializeWritePtr(ZoneNum, LapVarNum);
							break;
						case CalcGaussianCurvature:
							if (IsOk)
								IsOk = GaussCurvaturePtr.InitializeWritePtr(ZoneNum, GaussCurvatureVarNum);
							break;
						case CalcEigenVectorsDotGradient:
							for (int i = 0; i < 3 && IsOk; ++i)
								IsOk = EigVecDotGradPtrs[i].InitializeWritePtr(ZoneNum, EigVecDotGradVarNums[i]);
							break;
						case CalcEberly1Ridge:
							if (IsOk)
								IsOk = EberlyFuncPtrs[0].InitializeWritePtr(ZoneNum, EberlyFuncVarNums[0]);
							break;
						case CalcEberly2Ridge:
							if (IsOk)
								IsOk = EberlyFuncPtrs[1].InitializeWritePtr(ZoneNum, EberlyFuncVarNums[1]);
							break;
						case CalcEigenRank:
							if (IsOk)
								IsOk = EigenRankPtr.InitializeWritePtr(ZoneNum, EigenRankVarNum);
							break;
						default:
							TecUtilDialogErrMsg("Invalid variable used");
							IsOk = FALSE;
							break;
					}

					for (auto & i : Params){
						i.RhoPtr = &RhoPtr;
						i.GradPtrs = &GradPtrs;
						i.HessPtrs = &HessPtrs;
						i.HasHess = Opt.HasHess;
					}

					TecUtilZoneGetIJK(ZoneNum, &MaxIJK[0], &MaxIJK[1], &MaxIJK[2]);
					Boolean_t IsOrdered = TecUtilZoneIsOrdered(ZoneNum);
					if (!IsOrdered){
						/*
						*	For a FE zone, MaxIJK[0] is the actual number of points in the system.
						*	Can still use the triple loop, but need to set MaxIJK[1] and MaxIJK[2] to 1
						*/
						MaxIJK[1] = MaxIJK[2] = 1;
					}

					int NumCompleted = 0;
					int NumTotal = MaxIJK[0] / numCPU;

					if (HasRequiredVars || (IsOrdered && MaxIJK[2] > 1)){

						/*
						*	Loop over every point in the zone and calculate the necessary variable
						*/

						high_resolution_clock::time_point Time1;
						Time1 = high_resolution_clock::now();
#ifndef _DEBUG
#pragma omp parallel for
#endif
						for (int ii = 1; ii <= MaxIJK[0]; ++ii){
							int ThreadNum = omp_get_thread_num();

							if (ThreadNum == 0 && !Opt.CalcForAllZones){
								if (StatusUpdate(NumCompleted, NumTotal, StatusStr, Opt.AddOnID, Time1, false))
									NumCompleted++;
								else{
									IsOk = FALSE;
#pragma omp flush (IsOk)
								}

							}
#pragma omp flush (IsOk)
							for (int jj = 1; jj <= MaxIJK[1] && IsOk; ++jj){
								for (int kk = 1; kk <= MaxIJK[2] && IsOk; ++kk){
									if (ZoneNum == VolZoneNum){
										Params[ThreadNum].Index = IndexFromIJK(ii, jj, kk,
											VolInfo.MaxIJK[0], VolInfo.MaxIJK[1], VolInfo.MaxIJK[2],
											VolInfo.IsPeriodic) - 1;
									}
									else{
										if (IsOrdered)
											Params[ThreadNum].Index = IndexFromIJK(ii, jj, kk, MaxIJK[0], MaxIJK[1]) - 1;
										else
											Params[ThreadNum].Index = ii - 1;
									}

									ThreadCurPoint[ThreadNum] = VolInfo.MinXYZ + VolInfo.BasisNormalized * (vec3(vector<double>({ double(ii) - 1, double(jj) - 1, double(kk) - 1 }).data()) % VolInfo.PointSpacingV123);

									/*
									*	Get necessary information for calculating the variable
									*/
									if (*CalcVar > CalcEigenVectorsDotGradient){
										if (*CalcVar == CalcEigenRank
											|| (*CalcVar == CalcEberly1Ridge && (!Opt.EberlyUseCutoff[0] || RhoPtr[Params[ThreadNum].Index] >= Opt.EberlyCutoff[0]))
											|| (*CalcVar == CalcEberly2Ridge && (!Opt.EberlyUseCutoff[1] || RhoPtr[Params[ThreadNum].Index] >= Opt.EberlyCutoff[1])))
										{
											// 										if (HasEigVecDotGrad)
											// 											for (int i = 0; i < 3; ++i)
											// 												ThreadDotPtds[ThreadNum][i] = abs(EigVecDotGradPtrs[i][Index]);
											if (HasEigSys)
											for (int i = 0; i < 3; ++i)
											for (int j = 0; j < 3; ++j)
												ThreadEigVecs[ThreadNum].at(i, j) = EigSysPtrs[3 * i + j][Params[ThreadNum].Index];
											else{
												if (IsVolZone)
												{
													CalcEigenSystemForNode(ii,
														jj,
														kk,
														ThreadEigVals[ThreadNum],
														ThreadEigVecs[ThreadNum],
														Params[ThreadNum]);
												}
												else
												{
													CalcEigenSystemForPoint(ThreadCurPoint[ThreadNum],
														ThreadEigVals[ThreadNum],
														ThreadEigVecs[ThreadNum],
														Params[ThreadNum]);
												}
											}

											if (Opt.HasGrad)
											for (int i = 0; i < 3; ++i)
												ThreadGrad[ThreadNum][i] = GradPtrs[i][Params[ThreadNum].Index];
											else
											if (IsVolZone)
											{
												CalcGradForNode(ii,
													jj,
													kk,
													OneOverDelXYZx2,
													OneOverDelXYZx12,
													ThreadVals[ThreadNum],
													ThreadDirInd[ThreadNum],
													0,
													VolInfo.MaxIJK,
													VolInfo.IsPeriodic,
													RhoPtr,
													ThreadGrad[ThreadNum]);
											}
											else
											{
												CalcGradForPoint(ThreadCurPoint[ThreadNum],
													VolInfo.PointSpacingV123,
													*Params[ThreadNum].VolInfo,
													*Params[ThreadNum].BasisVectors,
													0,
													Opt.IsPeriodic,
													ThreadGrad[ThreadNum],
													RhoPtr,
													GPType_Invalid,
													reinterpret_cast<void*>(&Params[ThreadNum]));
											}

											// 											double GradMag = norm(ThreadGrad[ThreadNum]);
											// 											if (GradMag < 0.01)
											// 												ThreadGrad[ThreadNum] /= GradMag;
											// 											if (0){
											// 												ThreadGrad[ThreadNum] = normalise(ThreadGrad[ThreadNum]);
											// 											}
											// 											else{
											// 												for (int i = 0; i < 3; ++i)
											// 													ThreadEigVecs[ThreadNum][i] *= ThreadEigVals[ThreadNum][i];
											// 											}

											for (int i = 0; i < 3; ++i)
												ThreadDotPtds[ThreadNum][i] = abs(dot(ThreadEigVecs[ThreadNum].col(i), ThreadGrad[ThreadNum]));

											std::sort(ThreadDotPtds[ThreadNum].begin(), ThreadDotPtds[ThreadNum].end());
										}
										else
											ThreadDotPtds[ThreadNum] = { 0.0, 0.0, 0.0 };
									}
									else if (*CalcVar > CalcEigenSystem && HasEigSys){
										for (int i = 0; i < 3; ++i)
										for (int j = 0; j < 3; ++j)
											ThreadEigVecs[ThreadNum].at(i, j) = EigSysPtrs[3 * i + j][Params[ThreadNum].Index];

										for (int i = 0; i < 3; ++i)
											ThreadEigVals[ThreadNum][i] = EigSysPtrs[9 + i][Params[ThreadNum].Index];
									}
									else if (*CalcVar > CalcHessian){
										if (IsVolZone)
										{
											CalcEigenSystemForNode(ii,
												jj,
												kk,
												ThreadEigVals[ThreadNum],
												ThreadEigVecs[ThreadNum],
												Params[ThreadNum]);
										}
										else
										{
											CalcEigenSystemForPoint(ThreadCurPoint[ThreadNum],
												ThreadEigVals[ThreadNum],
												ThreadEigVecs[ThreadNum],
												Params[ThreadNum]);
										}

									}
									else if (*CalcVar > CalcGradientMagnitude){
										if (IsVolZone)
										{
											CalcHessForNode(ii,
												jj,
												kk,
												OneOverDelXYZx2,
												OneOverDelXYZx12,
												ThreadVals[ThreadNum],
												ThreadDirInd[ThreadNum],
												VolInfo.MaxIJK,
												VolInfo.IsPeriodic,
												RhoPtr,
												GradPtrs,
												ThreadHess[ThreadNum]);
										}
										else{
											if (Opt.HasGrad)
											{
												CalcHessFor3DPoint(ThreadCurPoint[ThreadNum],
													VolInfo.PointSpacingV123,
													*Params[ThreadNum].VolInfo,
													Opt.IsPeriodic,
													ThreadHessVec[ThreadNum],
													GradPtrs,
													GPType_Invalid,
													reinterpret_cast<void*>(&Params[ThreadNum]));
											}
											else
											{
												CalcHessForPoint(ThreadCurPoint[ThreadNum],
													VolInfo.PointSpacingV123,
													*Params[ThreadNum].VolInfo,
													*Params[ThreadNum].BasisVectors,
													Opt.IsPeriodic,
													ThreadHessVec[ThreadNum],
													RhoPtr,
													GPType_Invalid,
													reinterpret_cast<void*>(&Params[ThreadNum]));
											}
										}
									}
									else if (*CalcVar > CalcGradientVectors && Opt.HasGrad){
										for (int i = 0; i < 3; ++i)
											ThreadGrad[ThreadNum][i] = GradPtrs[i][Params[ThreadNum].Index];
									}
									else if (IsOk && *CalcVar > CalcInvalidVar){
										if (IsVolZone)
										{
											CalcGradForNode(ii,
												jj,
												kk,
												OneOverDelXYZx2,
												OneOverDelXYZx12,
												ThreadVals[ThreadNum],
												ThreadDirInd[ThreadNum],
												0,
												VolInfo.MaxIJK,
												VolInfo.IsPeriodic,
												RhoPtr,
												ThreadGrad[ThreadNum]);
										}
										else
										{
											CalcGradForPoint(ThreadCurPoint[ThreadNum],
												VolInfo.PointSpacingV123,
												*Params[ThreadNum].VolInfo,
												*Params[ThreadNum].BasisVectors,
												0,
												Opt.IsPeriodic,
												ThreadGrad[ThreadNum],
												RhoPtr,
												GPType_Invalid,
												reinterpret_cast<void*>(&Params[ThreadNum]));
										}
									}

									if (*CalcVar == CalcEigenVectorsDotGradient){
										if (Opt.HasGrad)
										for (int i = 0; i < 3; ++i)
											ThreadGrad[ThreadNum][i] = GradPtrs[i][Params[ThreadNum].Index];
										else{
											if (IsVolZone)
											{
												CalcGradForNode(ii,
													jj,
													kk,
													OneOverDelXYZx2,
													OneOverDelXYZx12,
													ThreadVals[ThreadNum],
													ThreadDirInd[ThreadNum],
													0,
													VolInfo.MaxIJK,
													VolInfo.IsPeriodic,
													RhoPtr,
													ThreadGrad[ThreadNum]);
											}
											else
											{
												CalcGradForPoint(ThreadCurPoint[ThreadNum],
													VolInfo.PointSpacingV123,
													*Params[ThreadNum].VolInfo,
													*Params[ThreadNum].BasisVectors,
													0,
													Opt.IsPeriodic,
													ThreadGrad[ThreadNum],
													RhoPtr,
													GPType_Invalid,
													reinterpret_cast<void*>(&Params[ThreadNum]));
											}
										}

										for (int i = 0; i < 3; ++i)
// 											ThreadDotPtds[ThreadNum][i] = dot(ThreadEigVecs[ThreadNum].col(i), ThreadGrad[ThreadNum]);
											ThreadDotPtds[ThreadNum][i] = dot(ThreadEigVecs[ThreadNum].row(i), ThreadGrad[ThreadNum]);
									}

									switch (*CalcVar){
										case CalcGradientVectors:
											for (int i = 0; i < 3; ++i)
												GradPtrs[i].Write(Params[ThreadNum].Index, ThreadGrad[ThreadNum][i]);
											break;
										case CalcGradientMagnitude:
											GradMagPtr.Write(Params[ThreadNum].Index, norm(ThreadGrad[ThreadNum]));
											break;
										case CalcHessian:
											for (int i = 0; i < 3; ++i)
											for (int j = i; j < 3; ++j)
												HessPtrs[HessIndices[i][j]].Write(Params[ThreadNum].Index, ThreadHess[ThreadNum].at(i, j));
											break;
										case CalcEigenSystem:
											for (int i = 0; i < 3; ++i)
											for (int j = 0; j < 3; ++j)
// 												EigSysPtrs[3 * i + j].Write(Params[ThreadNum].Index, ThreadEigVecs[ThreadNum].at(i, j));
												EigSysPtrs[3 * j + i].Write(Params[ThreadNum].Index, ThreadEigVecs[ThreadNum].at(i, j));

											for (int i = 0; i < 3; ++i)
												EigSysPtrs[9 + i].Write(Params[ThreadNum].Index, ThreadEigVals[ThreadNum][i]);
											break;
										case CalcLaplacian:
											TmpVal[ThreadNum] = 0.0;
											for (int i = 0; i < 3; ++i)
												TmpVal[ThreadNum] += ThreadEigVals[ThreadNum][i];

											LapPtr.Write(Params[ThreadNum].Index, TmpVal[ThreadNum]);
											break;
										case CalcGaussianCurvature:
											TmpVal[ThreadNum] = ThreadEigVals[ThreadNum][0];
											for (int i = 1; i < 3; ++i)
												TmpVal[ThreadNum] *= ThreadEigVals[ThreadNum][i];

											GaussCurvaturePtr.Write(Params[ThreadNum].Index, TmpVal[ThreadNum]);
											break;
										case CalcEigenVectorsDotGradient:
											for (int i = 0; i < 3; ++i)
												EigVecDotGradPtrs[i].Write(Params[ThreadNum].Index, ThreadDotPtds[ThreadNum][i]);
											break;
										case CalcEberly1Ridge:
											TmpVal[ThreadNum] = 0.0;
											for (int i = 0; i < 2; ++i)
												TmpVal[ThreadNum] += ThreadDotPtds[ThreadNum][i];

											EberlyFuncPtrs[0].Write(Params[ThreadNum].Index, TmpVal[ThreadNum]);
											break;
										case CalcEberly2Ridge:
											EberlyFuncPtrs[1].Write(Params[ThreadNum].Index, ThreadDotPtds[ThreadNum][0]);
											break;
										case CalcEigenRank:
											TmpVal[ThreadNum] = 0;
											for (int i = 0; i < 3; ++i){
												if (ThreadEigVals[ThreadNum][i] >= 0)
													++TmpVal[ThreadNum];
												else
													--TmpVal[ThreadNum];
											}
											EigenRankPtr.Write(Params[ThreadNum].Index, TmpVal[ThreadNum]);
											break;
										default:
											break;
									}
								}
							}
						}
					}
				}
				TecUtilDataLoadEnd();
			}

			/*
			*	Close all pointers
			*/
			RhoPtr.Close();
			GradMagPtr.Close();
			LapPtr.Close();
			GaussCurvaturePtr.Close();
			EigenRankPtr.Close();

			for (auto & i : GradPtrs) i.Close();
			for (auto & i : HessPtrs) i.Close();
			for (auto & i : EigSysPtrs) i.Close();
			for (auto & i : EigVecDotGradPtrs) i.Close();
			for (auto & i : EberlyFuncPtrs) i.Close();

			/*
			*	If volume zone calculation, quit loop.
			*/
			if (ZoneNum == Opt.CalcZoneNum && !Opt.CalcForAllZones)
				break;
		}

		StatusDrop(VolInfo.AddOnID);
		// 		TecUtilDialogDropPercentDone();

		// 		switch (*CalcVar){
		// 			case CalcGradVec:
		// 				break;
		// 			case CalcGradMag:
		// 				break;
		// 			case CalcHess:
		// 				break;
		// 			case CalcEigSys:
		// 				break;
		// 			case CalcLap:
		// 				break;
		// 			case CalcGaussianCurvature:
		// 				break;
		// 			case CalcEigVecDotGrad:
		// 				break;
		// 			case CalcEberly1Ridge:
		// 				break;
		// 			case CalcEberly2Ridge:
		// 				break;
		// 			case CalcEigenRank:
		// 				break;
		// 			default:
		// 				break;
		// 		}

		IterNum++;
	}

	TecUtilLockFinish(Opt.AddOnID);
	return;
}

vector<double> LogLevels(double Min, double Max)
{
	vector<double> LevelVector;
	double BaseValue;
	if (Max <= Min)
	{
		LevelVector.push_back(0.0);
	}
	else
	{
		int LevelInt = 0, MultInt = 1, ExpInt = 0;
		double TempValue = 0;
		Boolean_t BreakLoop = FALSE;
		TempValue = 1.0;

		if (Min == 0.0)
			BaseValue = 0.01;
		else{
			while (!BreakLoop)
			{
				if (TempValue > Min) TempValue *= 0.1;
				else if (TempValue * 10 <= Min) TempValue *= 10;
				else
				{
					BaseValue = TempValue;
					BreakLoop = TRUE;
				}
			}
		}

		BreakLoop = FALSE;

		while (!BreakLoop)
		{
			TempValue = MultInt * BaseValue * pow(10.0, ExpInt);
			if (TempValue >= Max)
			{
				TempValue = Max;
				BreakLoop = TRUE;
			}
			LevelVector.push_back(TempValue);

			if (MultInt >= 9)
			{
				MultInt = 1;
				ExpInt++;
			}
			else MultInt++;
		}
	}

	return LevelVector;
}

void GaussianBlur(Boolean_t IsPeriodic,
	AddOn_pa const & AddOnID,
	EntIndex_t const & ZoneNum,
	EntIndex_t const & VarNum,
	string const & NewVarName,
	double const & Sigma)
{
	FieldDataPointer_c Var;
	Var.InitializeReadPtr(ZoneNum, VarNum);
	if (!Var.IsReady()) return;

	int rSig = std::ceil(Sigma * 2.57);

	vector<FieldDataType_e> FD(TecUtilDataSetGetNumZones());
	for (int i = 0; i < FD.size(); ++i)
		FD[i] = TecUtilDataValueGetType(i + 1, VarNum);

	if (VarNumByName(NewVarName) > 0){
		Set_pa VarList = TecUtilSetAlloc(FALSE);
		TecUtilSetAddMember(VarList, VarNumByName(NewVarName), FALSE);
		TecUtilDataSetDeleteVar(VarList);
		TecUtilSetDealloc(&VarList);
	}

	if (!TecUtilDataSetAddVar(NewVarName.c_str(), FD.data())) return;

	FieldDataPointer_c NewVar;
	NewVar.InitializeWritePtr(ZoneNum, TecUtilDataSetGetNumVars());
	if (!NewVar.IsReady()) return;

	vector<int> IJK = NewVar.MaxIJK();

#pragma omp parallel for
	for (int i = 1; i <= IJK[0]; ++i){
		for (int j = 1; j <= IJK[1]; ++j){
			for (int k = 1; k <= IJK[2]; ++k){
				double Val = 0;
				double wSum = 0;

				for (int iz = k - rSig; iz < k + rSig + 1; ++iz){
					for (int iy = j - rSig; iy < j + rSig + 1; ++iy){
						for (int ix = i - rSig; ix < i + rSig + 1; ++ix){
							int x = MIN(IJK[0], MAX(1, ix));
							int y = MIN(IJK[1], MAX(1, iy));
							int z = MIN(IJK[2], MAX(1, iz));

							double dSqr = pow(ix - i, 2) + pow(iy - j, 2) + pow(iz - k, 2);
							double w = std::exp(-dSqr / (2. * Sigma * Sigma)) / pow(PI * 2. * Sigma * Sigma, 2);

							Val += Var[IndexFromIJK(x, y, z, IJK[0], IJK[1]) - 1] * w;
							wSum += w;
						}
					}
				}

				if (wSum != 0)
					NewVar.Write(IndexFromIJK(i, j, k, IJK[0], IJK[1]) - 1, Val / wSum);
			}
		}
	}

	TecUtilDialogMessageBox("Done!", MessageBoxType_Information);
}

// source channel, target channel, width, height, radius
// function gaussBlur_1(scl, tcl, w, h, r) {
// 	var rs = Math.ceil(r * 2.57);     // significant radius
// 	for (var i = 0; i < h; i++)
// 		for (var j = 0; j < w; j++) {
// 			var val = 0, wsum = 0;
// 			for (var iy = i - rs; iy < i + rs + 1; iy++)
// 				for (var ix = j - rs; ix < j + rs + 1; ix++) {
// 					var x = Math.min(w - 1, Math.max(0, ix));
// 					var y = Math.min(h - 1, Math.max(0, iy));
// 					var dsq = (ix - j)*(ix - j) + (iy - i)*(iy - i);
// 					var wght = Math.exp(-dsq / (2 * r*r)) / (Math.PI * 2 * r*r);
// 					val += scl[y*w + x] * wght;  wsum += wght;
// 				}
// 			tcl[i*w + j] = Math.round(val / wsum);
// 		}
// }
// 

vec3 Transform2dTo3d(vec2 const & TwoPt,
	mat33 const & BasisVectors,
	vec3 const & Origin)
{
	vec3 ThreePt(Origin);
	for (int i = 0; i < 2; ++i){
		ThreePt += BasisVectors.col(i) * TwoPt[i];
	}

	return ThreePt;
}