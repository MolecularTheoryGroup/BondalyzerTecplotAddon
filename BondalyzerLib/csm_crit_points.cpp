#include <vector>
#include <string>

#include "omp.h"

#include "TECADDON.h"

#include "CSM_DATA_SET_INFO.h"
#include "CSM_VOL_EXTENT_INDEX_WEIGHTS.h"
#include "CSM_CALC_VARS.h"

#include "CSM_CRIT_POINTS.h"

using std::vector;
using std::string;
using std::to_string;


/*
*	Begin CritPoints_c methods
*/

/*
	*	Constructors/destructors
	*/

CritPoints_c::CritPoints_c()
{
	m_TotNumCPs = 0;
	for (int i = 0; i < 6; ++i)
		m_NumCPs[i] = 0;

	m_MinCPDist = -1;
	m_MinCPDistFound = FALSE;
	m_RhoCutoff = -1;
	m_Dimensions = -1;
}

CritPoints_c::CritPoints_c(const double & RhoCutoff, const int & NumDimensions){
	m_TotNumCPs = 0;
	for (int i = 0; i < 6; ++i)
		m_NumCPs[i] = 0;

	if (RhoCutoff >= 0)
		m_RhoCutoff = RhoCutoff;

	if (NumDimensions > 0)
		m_Dimensions = NumDimensions;

	m_MinCPDist = -1;
	m_MinCPDistFound = FALSE;
}

CritPoints_c::CritPoints_c(const vector<CritPoints_c> & CPLists){
	m_TotNumCPs = 0;
	for (int i = 0; i < 6; ++i)
		m_NumCPs[i] = 0;

	for (auto Beg = CPLists.cbegin(), End = CPLists.cend(); Beg != End; Beg++)
		this->Append(*Beg);

	m_MinCPDist = -1;
	m_MinCPDistFound = FALSE;
}

CritPoints_c::CritPoints_c(const int & CPZoneNum, 
	const vector<int> & XYZVarNums, 
	const int & RhoVarNum, 
	const int & CPTypeVarNum,
	MultiRootParams_s *MR) : CritPoints_c()
{
	int NumZones = TecUtilDataSetGetNumZones();
	int NumVars = TecUtilDataSetGetNumVars();

	REQUIRE(0 < CPZoneNum && CPZoneNum <= NumZones);
	REQUIRE(0 < RhoVarNum && RhoVarNum <= NumVars);
	REQUIRE(0 < CPTypeVarNum && CPTypeVarNum <= NumVars);
	REQUIRE(XYZVarNums.size() == 3);
	for (const auto & i : XYZVarNums) REQUIRE(0 < i && i < NumVars);

	FieldDataPointer_c RhoPtr, CPTypePtr;
	vector<FieldDataPointer_c> XYZPtrs(3);

	TecUtilDataLoadBegin();

	if (!RhoPtr.GetReadPtr(CPZoneNum, RhoVarNum)){
		TecUtilDialogErrMsg("Failed to get CP rho pointer");
		return;
	}
	if (!CPTypePtr.GetReadPtr(CPZoneNum, CPTypeVarNum)){
		TecUtilDialogErrMsg("Failed to get CP CPType pointer");
		return;
	}
	for (int i = 0; i < 3; ++i) if (!XYZPtrs[i].GetReadPtr(CPZoneNum, XYZVarNums[i])){
		TecUtilDialogErrMsg("Failed to get CP XYZ pointer");
		return;
	}

	for (int iCP = 0; iCP < RhoPtr.GetSize(); ++iCP){
		int ind = std::find(CPTypeList, CPTypeList + 6, CPTypePtr[iCP]) - CPTypeList;

		m_Rho[ind].push_back(RhoPtr[iCP]);
		m_XYZ[ind].push_back(vec3());
		for (int d = 0; d < 3; ++d) m_XYZ[ind].back()[d] = XYZPtrs[d][iCP];
	}

	if (MR != NULL){
		vec3 EigVals;
		mat33 EigVecs;

		for (int t = 0; t < 4; ++t){ // only calculate principal directions for near field CPs
			for (auto p : m_XYZ[t]){
				CalcEigenSystemForPoint(p, EigVals, EigVecs, *MR);
				m_PrincDir[t].push_back(normalise(EigVecs.col(CPPrincDirInds[t])));
				m_EigVals[t].push_back(EigVals);
				m_EigVecs[t].push_back(EigVecs);
			}
		}
	}

	TecUtilDataLoadEnd();

	for (int i = 0; i < 6; ++i){
		m_NumCPs[i] = m_Rho[i].size();
		m_TotNumCPs += m_Rho[i].size();
	}
}

CritPoints_c::~CritPoints_c()
{
}

/*
	*	Operator overloads
	*/

CritPoints_c & CritPoints_c::operator+=(const CritPoints_c & rhs){
	this->Append(rhs);

	return *this;
}

/*
	*	Getter methods
	*/

const double CritPoints_c::GetMinCPDist(){
	if (m_MinCPDistFound)
		return m_MinCPDist;
	else if (FindMinCPDist())
		return m_MinCPDist;
	else
		return -1;
}

vec3 CritPoints_c::GetXYZ(const int & TotOffset) const{
	vector<int> TypeNumOffset = GetTypeNumOffsetFromTotOffset(TotOffset);
	if (TypeNumOffset[0] >= 0)
		return GetXYZ(TypeNumOffset[0], TypeNumOffset[1]);

	return vec3();
}

const double CritPoints_c::GetRho(const int & TotOffset) const{
	vector<int> TypeNumOffset = GetTypeNumOffsetFromTotOffset(TotOffset);
	if (TypeNumOffset[0] >= 0)
		return GetRho(TypeNumOffset[0], TypeNumOffset[1]);

	return -1;
}

vec3 CritPoints_c::GetPrincDir(const int & TotOffset) const{
	vector<int> TypeNumOffset = GetTypeNumOffsetFromTotOffset(TotOffset);
	if (TypeNumOffset[0] >= 0)
		return GetPrincDir(TypeNumOffset[0], TypeNumOffset[1]);

	return vec3();
}
vec3 CritPoints_c::GetEigVals(const int & TotOffset) const{
	vector<int> TypeNumOffset = GetTypeNumOffsetFromTotOffset(TotOffset);
	if (TypeNumOffset[0] >= 0)
		return GetEigVals(TypeNumOffset[0], TypeNumOffset[1]);

	return vec3();
}
mat33 CritPoints_c::GetEigVecs(const int & TotOffset) const{
	vector<int> TypeNumOffset = GetTypeNumOffsetFromTotOffset(TotOffset);
	if (TypeNumOffset[0] >= 0)
		return GetEigVecs(TypeNumOffset[0], TypeNumOffset[1]);

	return mat33();
}

const Boolean_t CritPoints_c::IsValid() const{
	Boolean_t IsOk = TRUE;
	size_t TmpInt;

	for (int i = 0; i < 6 && IsOk; ++i){
		TmpInt = m_Rho[i].size();
		IsOk = (TmpInt == m_XYZ[i].size() && TmpInt == m_PrincDir[i].size());
	}

	return IsOk;
}

/*
	*	Setter methods
	*/

const Boolean_t CritPoints_c::AddPoint(const double & Rho,
	const vec3 & Pos,
	const vec3 & PrincDir,
	const char & Type)
{
	Boolean_t IsOk = FALSE;

	for (int i = 0; i < 6 && !IsOk; ++i){
		if (Type == CPTypeList[i]){
			m_Rho[i].push_back(Rho);
			m_XYZ[i].push_back(Pos);
			m_PrincDir[i].push_back(PrincDir);
			m_NumCPs[i]++;
			m_TotNumCPs++;
			IsOk = TRUE;
		}
	}

	return IsOk;
}

void CritPoints_c::Append(const CritPoints_c & rhs)
{
	for (int i = 0; i < 6; ++i){
		m_Rho[i].insert(m_Rho[i].end(), rhs.m_Rho[i].cbegin(), rhs.m_Rho[i].cend());
		m_XYZ[i].insert(m_XYZ[i].end(), rhs.m_XYZ[i].cbegin(), rhs.m_XYZ[i].cend());
		m_PrincDir[i].insert(m_PrincDir[i].end(), rhs.m_PrincDir[i].cbegin(), rhs.m_PrincDir[i].cend());
		m_NumCPs[i] += rhs.m_NumCPs[i];
	}
	m_TotNumCPs += rhs.m_TotNumCPs;
}


/*
	*	Mutators and other methods
	*/

const Boolean_t CritPoints_c::FindMinCPDist(){
	Boolean_t IsOk = (m_TotNumCPs > 0);
	int MinI = -1, MinJ = -1;
	double TmpDbl;

	m_MinCPDist = DBL_MAX;

	if (IsOk){
		for (int i = 0; i < m_TotNumCPs; ++i){
			for (int j = i + 1; j < m_TotNumCPs; ++j){
				TmpDbl = DistSqr(GetXYZ(i), GetXYZ(j));
				if (TmpDbl < m_MinCPDist){
					m_MinCPDist = TmpDbl;
					MinI = i;
					MinJ = j;
				}
			}
		}
		IsOk = (MinI >= 0 && MinJ >= 0);
	}

	if (IsOk){
		m_MinCPDist = Distance(GetXYZ(MinI), GetXYZ(MinJ));
		m_MinCPDistFound = TRUE;
	}

	return IsOk;
}

vector<int> CritPoints_c::GetTypeNumOffsetFromTotOffset(const int & TotOffset) const{
	vector<int> TypeNumAndOffset = { -1, -1 };

	int CPCount = 0;
	for (int i = 0; i < 6; ++i){
		if (TotOffset < CPCount + m_NumCPs[i]){
			TypeNumAndOffset[0] = i;
			TypeNumAndOffset[1] = TotOffset - CPCount;
			break;
		}
		CPCount += m_NumCPs[i];
	}

	return TypeNumAndOffset;
}

const int CritPoints_c::GetTotOffsetFromTypeNumOffset(const int & TypeNum, const int & TypeOffset) const
{
	int TotOffset = 0;
	for (int i = 0; i < TypeNum; ++i) TotOffset += m_NumCPs[i];
	TotOffset += TypeOffset;

	return TotOffset;
}


/*
	*	Private methods
	*/


/*
 *	Mutators and other methods
 */

const vector<int> CritPoints_c::SaveAsOrderedZone(const vector<int> & XYZVarNum, const int & RhoVarNum, const Boolean_t & SaveCPTypeZones){
	for (const auto & i : XYZVarNum) REQUIRE(i > 0 && i <= TecUtilDataSetGetNumVars());

	vector<int> NewZoneNums;
	/*
	 *	Make CP type variable if it doesn't exist yet.
	 */
	int CPTypeVarNum = VarNumByName(CPTypeVarName);
	FieldDataPointer_c CPTypePtr;
	if (CPTypeVarNum <= 0){
		/*
		 *	There is no CP type variable yet, which means there are also
		 *	no CP zones, so the data type for all current zones can be 
		 *	bit for the CP type variable.
		 */
		vector<FieldDataType_e> DataTypes(TecUtilDataSetGetNumZones(), FieldDataType_Bit);
		if (!TecUtilDataSetAddVar(CPTypeVarName.c_str(), DataTypes.data())){
			TecUtilDialogErrMsg("Failed to create CP type variable for CP zone");
			return{ -1 };
		}
		CPTypeVarNum = TecUtilDataSetGetNumVars();
	}
	vector<FieldDataType_e> DataTypes(TecUtilDataSetGetNumVars());
	for (int i = 0; i < DataTypes.size() - 1; ++i)
		DataTypes[i] = TecUtilDataValueGetType(1, i + 1);
	DataTypes.back() = FieldDataType_Int16;

	if (!TecUtilDataSetAddZone("Critical Points", NumCPs(), 1, 1, ZoneType_Ordered, DataTypes.data())){
		TecUtilDialogErrMsg("Failed to create CP zone");
		return{ -1 };
	}
	NewZoneNums.push_back(TecUtilDataSetGetNumZones());

	AuxDataZoneSetItem(NewZoneNums.back(), CSMAuxData.CC.ZoneType, CSMAuxData.CC.ZoneTypeCPs);

	Set_pa CPZoneSet = TecUtilSetAlloc(TRUE);
	TecUtilSetAddMember(CPZoneSet, NewZoneNums.back(), TRUE);

	TecUtilZoneSetActive(CPZoneSet, AssignOp_PlusEquals);

	vector<FieldDataPointer_c> XYZPtrs(3);
	FieldDataPointer_c RhoPtr;

	TecUtilDataLoadBegin();

	for (int i = 0; i < 3; ++i) if (!XYZPtrs[i].GetWritePtr(NewZoneNums.back(), XYZVarNum[i])){
		TecUtilDialogErrMsg("Failed to get XYZ pointers for CP zone");
		return{ -1 };
	}
	if (!RhoPtr.GetWritePtr(NewZoneNums.back(), RhoVarNum)){
		TecUtilDialogErrMsg("Failed to get rho pointer for CP zone");
		return{ -1 };
	}
	if (!CPTypePtr.GetWritePtr(NewZoneNums.back(), CPTypeVarNum)){
		TecUtilDialogErrMsg("Failed to get rho pointer for CP zone");
		return{ -1 };
	}

	int CPNum = 0;

	for (int t = 0; t < 6; ++t){
		AuxDataZoneSetItem(NewZoneNums.back(), CSMAuxData.CC.NumCPs[t], to_string(NumCPs(t)));
		for (int ti = 0; ti < NumCPs(t); ++ti){
			for (int d = 0; d < 3; ++d){
				XYZPtrs[d].Write(CPNum, GetXYZ(t, ti)[d]);
			}
			CPTypePtr.Write(CPNum++, CPTypeList[t]);
		}
	}

	TecUtilDataLoadEnd();

	if (SaveCPTypeZones){
		TecUtilZoneSetActive(CPZoneSet, AssignOp_MinusEquals);
		Set_pa CPTypeZoneSet = TecUtilSetAlloc(TRUE);
		for (int t = 0; t < 6; ++t){
			if (NumCPs(t) > 0){
// 				DataTypes.push_back(FieldDataType_Int16);
				if (!TecUtilDataSetAddZone(("Critical Points: " + CPNameList[t]).c_str(), NumCPs(t), 1, 1, ZoneType_Ordered, DataTypes.data())){
					TecUtilDialogErrMsg("Failed to create CP type zone");
					return{ -1 };
				}
				NewZoneNums.push_back(TecUtilDataSetGetNumZones());

				AuxDataZoneSetItem(NewZoneNums.back(), CSMAuxData.CC.ZoneType, CSMAuxData.CC.ZoneTypeCPs);
				AuxDataZoneSetItem(NewZoneNums.back(), CSMAuxData.CC.ZoneSubType, CSMAuxData.CC.CPSubTypes[t]);
				AuxDataZoneSetItem(NewZoneNums.back(), CSMAuxData.CC.NumCPs[t], to_string(NumCPs(t)));

				TecUtilSetAddMember(CPZoneSet, NewZoneNums.back(), TRUE);
				TecUtilSetClear(CPTypeZoneSet);
				TecUtilSetAddMember(CPTypeZoneSet, NewZoneNums.back(), TRUE);
				TecUtilZoneSetScatter(SV_COLOR, CPTypeZoneSet, 0.0, CPColorList[t]);

				TecUtilDataLoadBegin();

				for (auto & i : XYZPtrs) i.Close();
				for (int i = 0; i < 3; ++i) if (!XYZPtrs[i].GetWritePtr(NewZoneNums.back(), XYZVarNum[i])){
					TecUtilDialogErrMsg("Failed to get xyz pointers for CP  type zone");
					return{ -1 };
				}
				RhoPtr.Close();
				if (!RhoPtr.GetWritePtr(NewZoneNums.back(), RhoVarNum)){
					TecUtilDialogErrMsg("Failed to get rho pointer for CP type zone");
					return{ -1 };
				}
				FieldDataPointer_c CPTypeZonePtr;
				if (!CPTypeZonePtr.GetWritePtr(NewZoneNums.back(), CPTypeVarNum)){
					TecUtilDialogErrMsg("Failed to get cp type pointer for CP type zone");
					return{ -1 };
				}
				for (int ti = 0; ti < NumCPs(t); ++ti){
// 					string str = to_string(CPTypeList[t]) + " at ";
					for (int d = 0; d < 3; ++d){
						XYZPtrs[d].Write(ti, GetXYZ(t, ti)[d]);
// 						str += to_string(GetXYZ(t, ti)[d]) + ", ";
					}
					CPTypeZonePtr.Write(ti, CPTypeList[t]);
// 					TecUtilDialogMessageBox(str.c_str(), MessageBoxType_Information);
				}

				CPTypeZonePtr.Close();

				TecUtilZoneSetActive(CPTypeZoneSet, AssignOp_PlusEquals);

				TecUtilDataLoadEnd();
			}
		}
		TecUtilZoneSetContour(SV_SHOW, CPZoneSet, 0.0, FALSE);
	}
	else{
		TecUtilZoneSetActive(CPZoneSet, AssignOp_PlusEquals);
		TecUtilZoneSetContour(SV_SHOW, CPZoneSet, 0.0, TRUE);
	}


	TecUtilZoneSetMesh(SV_SHOW, CPZoneSet, 0.0, FALSE);
	TecUtilZoneSetShade(SV_SHOW, CPZoneSet, 0.0, FALSE);
	TecUtilZoneSetScatter(SV_SHOW, CPZoneSet, 0.0, TRUE);
	TecUtilZoneSetScatter(SV_FRAMESIZE, CPZoneSet, 1, FALSE);
	TecUtilZoneSetScatterSymbolShape(SV_GEOMSHAPE, CPZoneSet, GeomShape_Sphere);

	return NewZoneNums;
}


/*
*	End CritPoints_c methods
*/

/*
	*	Functions for the GSL MultiRoots root finder
	*/

/*
	*	Function to return the actual function (gradient) value
	*/

int F3D(const gsl_vector * pos, void * params, gsl_vector * GradValues){

	MultiRootParams_s * RootParams = reinterpret_cast<MultiRootParams_s*>(params);

	vec3 Point(pos->data);

	if (!SetIndexAndWeightsForPoint(Point, *RootParams->VolInfo))
		return GSL_ESANITY;

	// 	for (int i = 0; i < 3; ++i)
	// 		Point[i] = ValByCurrentIndexAndWeightsFromRawPtr(*RootParams->VolInfo, RootParams->GradPtrs->at(i));
	// 
	// 	if (RootParams->HasHess)
	// 		Point *= 0.1;

	// 	if (RootParams->HasHess){
	// 		for (int i = 0; i < 3; ++i)
	// 			gsl_vector_set(GradValues, i, ValByCurrentIndexAndWeightsFromRawPtr(*RootParams->VolInfo, RootParams->GradPtrs->at(i)));
	// 	}
	// 	else{
	
	if (RootParams->HasGrad) for (int i = 0; i < 3; ++i){
		gsl_vector_set(GradValues, i, ValByCurrentIndexAndWeightsFromRawPtr(*RootParams->VolInfo, RootParams->GradPtrs->at(i)));
	}
	else{
		vec3 Grad;
		CalcGradForPoint(Point, RootParams->VolInfo->DelXYZ, *RootParams->VolInfo, *RootParams->BasisVectors, 0, RootParams->IsPeriodic, Grad, *RootParams->RhoPtr, GPType_Invalid, params);
		for (int i = 0; i < 3; ++i){
			gsl_vector_set(GradValues, i, Grad[i]);
		}
	}
	// 	}

	return GSL_SUCCESS;
}

/*
	*	Function to return the derivatives (jacobian matrix of gradient, ie Hessian of rho)
	*/

int DF3D(const gsl_vector * pos, void * params, gsl_matrix * Jacobian){

	MultiRootParams_s * RootParams = reinterpret_cast<MultiRootParams_s*>(params);

	vec3 Point(pos->data);

	if (!SetIndexAndWeightsForPoint(Point, *RootParams->VolInfo))
		return GSL_ESANITY;

	if (RootParams->HasHess){
		/*
			*	Analytical Hessian available, so use that.
			*/

		int HessIndices[3][3] = {
			{ 0, 1, 2 },
			{ 1, 3, 4 },
			{ 2, 4, 5 }
		};

		for (int i = 0; i < 3; ++i){
			for (int j = 0; j < 3; ++j){
				if (j >= i)
					gsl_matrix_set(Jacobian, i, j, ValByCurrentIndexAndWeightsFromRawPtr(*RootParams->VolInfo, RootParams->HessPtrs->at(HessIndices[j][i])));
				else
					gsl_matrix_set(Jacobian, i, j, gsl_matrix_get(Jacobian, j, i));
			}
		}
	}
	else{
		/*
			*	No analytical Hessian, so need to find derivative numerically.
			*	Need to do it manually, since the GSL solver doesn't know not to
			*	go beyond the bounds of the system.
			*/
		mat33 Hess;
		if (RootParams->HasGrad){
			CalcHessFor3DPoint(Point,
				RootParams->VolInfo->DelXYZ,
				*RootParams->VolInfo,
				RootParams->IsPeriodic,
				Hess,
				*RootParams->GradPtrs,
				GPType_Invalid,
				params);
		}
		else{
			CalcHessForPoint(Point,
				RootParams->VolInfo->DelXYZ,
				*RootParams->VolInfo,
				*RootParams->BasisVectors,
				RootParams->IsPeriodic,
				Hess,
				*RootParams->RhoPtr,
				GPType_Invalid,
				params);
		}

		for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j) gsl_matrix_set(Jacobian, i, j, Hess.at(i, j));
	}

	return GSL_SUCCESS;
}

/*
	*	Function to return both grad and dgrad (function and derivatives)
	*/

int FDF3D(const gsl_vector * pos, void * params, gsl_vector * GradValues, gsl_matrix * Jacobian){

	int Status = F3D(pos, params, GradValues);

	if (Status == GSL_SUCCESS)
		Status = DF3D(pos, params, Jacobian);

	return Status;
}


const Boolean_t CritPointInCell(const vector<int> & IJK,
	vec3 & Point,
	vec3 & PrincDir,
	double & RhoValue,
	const double & RhoCutoff,
	char & Type,
	MultiRootParams_s & RootParams,
	MultiRootObjects_s & MR)
{
	// 	TecUtilDialogMessageBox(string("CritPointInCell, rho ptr good " + to_string(RootParams.RhoPtr->IsReady())).c_str(), MessageBoxType_Information);

	Boolean_t CPInCell = TRUE;

	Type = 0;
	int Status = GSL_SUCCESS;
	int Iter = 0;

	vec3 MinCellXYZ, MaxCellXYZ;

	for (int i = 0; i < 3; ++i){
		MinCellXYZ[i] = RootParams.VolInfo->DelXYZ[i] * static_cast<double>(IJK[i]) + RootParams.VolInfo->MinXYZ[i];
		MaxCellXYZ[i] = MinCellXYZ[i] + RootParams.VolInfo->DelXYZ[i];
		Point[i] = MinCellXYZ[i] + RootParams.VolInfo->DelXYZ[i] * 0.5;
		gsl_vector_set(MR.pos, i, Point[i]);
	}

	// 	TecUtilDialogMessageBox("start point set", MessageBoxType_Information);

	CPInCell = SetIndexAndWeightsForPoint(Point, *RootParams.VolInfo);
	// 	TecUtilDialogMessageBox("initial point indexed", MessageBoxType_Information);
	if (CPInCell){
		RhoValue = ValByCurrentIndexAndWeightsFromRawPtr(*RootParams.VolInfo, *RootParams.RhoPtr);
		CPInCell = (RhoValue >= RhoCutoff);
	}

	// 	TecUtilDialogMessageBox("initial point checked", MessageBoxType_Information);

	if (CPInCell){
		gsl_multiroot_fdfsolver_set(MR.s, &MR.Func, MR.pos);


		// 		TecUtilDialogMessageBox("solver set", MessageBoxType_Information);

		do
		{
			++Iter;

			// 			string str = "iteration " + to_string(Iter);
			// 			TecUtilDialogMessageBox(str.c_str(), MessageBoxType_Information);

			Status = gsl_multiroot_fdfsolver_iterate(MR.s);

			if (Status != GSL_CONTINUE && Status != GSL_SUCCESS)
				break;

			Point = MR.s->x->data;

			Status = gsl_multiroot_test_residual(MR.s->f, 1e-7);

			if (Iter > CheckPosIter || Iter >= MaxCPIter){
				CPInCell = sum(Point >= MinCellXYZ) == 3 && sum(Point <= MaxCellXYZ) == 3;
			}

		} while (CPInCell && Status == GSL_CONTINUE && Iter < MaxCPIter);

		Point = MR.s->x->data;

		CPInCell = sum(Point >= MinCellXYZ) == 3 && sum(Point <= MaxCellXYZ) == 3;

		RhoValue = 0.0;
		if (CPInCell){
			CPInCell = SetIndexAndWeightsForPoint(Point, *RootParams.VolInfo);
			RhoValue = ValByCurrentIndexAndWeightsFromRawPtr(*RootParams.VolInfo, *RootParams.RhoPtr);
		}

		if (CPInCell && RhoValue >= RhoCutoff){
			vec3 EigVals;
			mat33 EigVecs;

			CalcEigenSystemForPoint(Point,
				EigVals,
				EigVecs,
				RootParams);

			for (int i = 0; i < 3; ++i){
				if (EigVals[i] > 0)
					Type++;
				else
					Type--;
			}
			if (Type == ATOMCP || Type == RINGCP)
				PrincDir = EigVecs.row(0).t();
			else
				PrincDir = EigVecs.row(2).t();
		}
	}

	return CPInCell;
}

const Boolean_t CritPointInCell(
	const vec3 & CellMinXYZ,
	const vec3 & CellMaxXYZ,
	vec3 & Point,
	vec3 & PrincDir,
	double & RhoValue,
	const double & RhoCutoff,
	char & Type,
	MultiRootParams_s & RootParams,
	MultiRootObjects_s & MR)
{
	// 	TecUtilDialogMessageBox(string("CritPointInCell, rho ptr good " + to_string(RootParams.RhoPtr->IsReady())).c_str(), MessageBoxType_Information);

	Boolean_t CPInCell = TRUE;

	Type = 0;
	int Status = GSL_SUCCESS;
	int Iter = 0;

	Point = (CellMaxXYZ + CellMinXYZ) * 0.5;
	for (int i = 0; i < 3; ++i) gsl_vector_set(MR.pos, i, Point[i]);

	// 	TecUtilDialogMessageBox("start point set", MessageBoxType_Information);

	CPInCell = SetIndexAndWeightsForPoint(Point, *RootParams.VolInfo);
	// 	TecUtilDialogMessageBox("initial point indexed", MessageBoxType_Information);
	if (CPInCell){
		RhoValue = ValByCurrentIndexAndWeightsFromRawPtr(*RootParams.VolInfo, *RootParams.RhoPtr);
		CPInCell = (RhoValue >= RhoCutoff);
	}

	// 	TecUtilDialogMessageBox("initial point checked", MessageBoxType_Information);

	if (CPInCell){
		gsl_multiroot_fdfsolver_set(MR.s, &MR.Func, MR.pos);


		// 		TecUtilDialogMessageBox("solver set", MessageBoxType_Information);

		do
		{
			++Iter;

			// 			string str = "iteration " + to_string(Iter);
			// 			TecUtilDialogMessageBox(str.c_str(), MessageBoxType_Information);

			Status = gsl_multiroot_fdfsolver_iterate(MR.s);

			if (Status != GSL_CONTINUE && Status != GSL_SUCCESS)
				break;

			Point = MR.s->x->data;

			Status = gsl_multiroot_test_residual(MR.s->f, 1e-7);

			if (Iter > CheckPosIter || Iter >= MaxCPIter){
				CPInCell = sum(Point >= CellMinXYZ) == 3 && sum(Point <= CellMaxXYZ) == 3;
			}

		} while (CPInCell && Status == GSL_CONTINUE && Iter < MaxCPIter);

		Point = MR.s->x->data;

		CPInCell = sum(Point >= CellMinXYZ) == 3 && sum(Point <= CellMaxXYZ) == 3;

		RhoValue = 0.0;
		if (CPInCell){
			CPInCell = SetIndexAndWeightsForPoint(Point, *RootParams.VolInfo);
			RhoValue = ValByCurrentIndexAndWeightsFromRawPtr(*RootParams.VolInfo, *RootParams.RhoPtr);
		}

		if (CPInCell && RhoValue >= RhoCutoff){
			vec3 EigVals;
			mat33 EigVecs;

			CalcEigenSystemForPoint(Point,
				EigVals,
				EigVecs,
				RootParams);

			for (int i = 0; i < 3; ++i){
				if (EigVals[i] > 0)
					Type++;
				else
					Type--;
			}
			if (Type == ATOMCP || Type == RINGCP)
				PrincDir = EigVecs.row(0).t();
			else
				PrincDir = EigVecs.row(2).t();
		}
	}

	return CPInCell;
}



/*
*	Function for searching a subzone (ordered IJK)
*	for critical points.
*/
const Boolean_t FindCPs(CritPoints_c & CPs,
	VolExtentIndexWeights_s VolInfo,
	const Boolean_t & IsPeriodic,
	const vector<int> & StartIJK,
	const vector<int> & EndIJK,
	const FieldDataPointer_c & RhoPtr,
	const vector<FieldDataPointer_c> & GradXYZPtrs,
	const vector<FieldDataPointer_c> & HessPtrs)
{
	Boolean_t IsOk = (GradXYZPtrs.size() == 3
		&& StartIJK.size() == 3 && EndIJK.size() == 3
		&& (HessPtrs.size() == 0 || HessPtrs.size() == 6));

	if (!IsOk) return IsOk;

	vec3 TmpPoint, PrincDir;
	double TmpRho;
	char TmpType;

	MultiRootParams_s RootParams;
	RootParams.VolInfo = &VolInfo;
	RootParams.IsPeriodic = IsPeriodic;
	RootParams.RhoPtr = &RhoPtr;
	RootParams.GradPtrs = &GradXYZPtrs;
	RootParams.HessPtrs = &HessPtrs;

// 	TecUtilDialogMessageBox(string("FindCPs, RhoPtr.IsReady() = " + to_string(RhoPtr.IsReady())).c_str(), MessageBoxType_Information);
// 
// 	TecUtilDialogMessageBox(string("FindCPs, RootParams.RhoPtr.IsReady() = " + to_string(RootParams.RhoPtr->IsReady())).c_str(), MessageBoxType_Information);

	vector<vec3> BV(3);
	BV[0] << 1 << 0 << 0;
	BV[1] << 0 << 1 << 0;
	BV[2] << 0 << 0 << 1;
	RootParams.BasisVectors = &BV;

	RootParams.HasHess = HessPtrs.size() == 6;

	for (int i = 0; i < 6 && RootParams.HasHess; ++i)
		RootParams.HasHess = RootParams.HessPtrs->at(i).IsReady();

	MultiRootObjects_s MR;

	MR.Func = { &F3D, &DF3D, &FDF3D, 3, &RootParams };
	MR.pos = gsl_vector_alloc(3);
	MR.T = gsl_multiroot_fdfsolver_newton;
	MR.s = gsl_multiroot_fdfsolver_alloc(MR.T, 3);

	string StatusStr = "Finding critical points";
	int Range = EndIJK[2] - StartIJK[2];
	int kNum = 0;


	vector<int> IJK(3);

	if (StartIJK[2] == 1)
		TecUtilDialogLaunchPercentDone(StatusStr.c_str(), TRUE);

	for (IJK[2] = StartIJK[2]; IJK[2] <= EndIJK[2] && IsOk; ++IJK[2]){
		kNum++;
		if (StartIJK[2] == 1 && !SetPercent(kNum, Range, StatusStr, VolInfo.AddOnID)){
			IsOk = FALSE;
		}
		for (IJK[1] = StartIJK[1]; IJK[1] <= EndIJK[1] && IsOk; ++IJK[1]){
			for (IJK[0] = StartIJK[0]; IJK[0] <= EndIJK[0] && IsOk; ++IJK[0]){
// 				string str = "{i,j,k} = {";
// 				for (int i = 0; i < 3; ++i) str += to_string(IJK[i]) + ", ";
// 				TecUtilDialogMessageBox(str.c_str(), MessageBoxType_Information);
				if (CritPointInCell(IJK, TmpPoint, PrincDir, TmpRho, CPs.GetRhoCutoff(), TmpType, RootParams, MR)
					&& TmpType != 0)
				{
					IsOk = CPs.AddPoint(TmpRho, TmpPoint, PrincDir, TmpType);
				}
			}
		}
	}

	if (StartIJK[2] == 1)
		TecUtilDialogDropPercentDone();

	if (MR.s != NULL)
		gsl_multiroot_fdfsolver_free(MR.s);
	// 	if (MR.pos != NULL)
	// 		gsl_vector_free(MR.pos);

	return IsOk;
}



/*
*	Function for searching a subzone (ordered IJK)
*	for critical points using aribrary cells.
*/
const Boolean_t FindCPs(CritPoints_c & CPs,
	const VolExtentIndexWeights_s & VolInfo,
	const double & CellSpacing,
	const double & RhoCutoff,
	const Boolean_t & IsPeriodic,
	const FieldDataPointer_c & RhoPtr,
	const vector<FieldDataPointer_c> & GradXYZPtrs,
	const vector<FieldDataPointer_c> & HessPtrs)
{
	Boolean_t IsOk = ((GradXYZPtrs.size() == 3 || GradXYZPtrs.size() == 0)
		&& (HessPtrs.size() == 0 || HessPtrs.size() == 6));

	if (!IsOk) return IsOk;

	int NumThreads = omp_get_num_procs();

	vector<VolExtentIndexWeights_s> VolInfoList(NumThreads, VolInfo);

	vector<vec3> BV(3);
	BV[0] << 1 << 0 << 0;
	BV[1] << 0 << 1 << 0;
	BV[2] << 0 << 0 << 1;

	vector<MultiRootParams_s> RootParams(NumThreads);
	for (int r = 0; r < NumThreads; ++r){
		RootParams[r].VolInfo = &VolInfoList[r];
		RootParams[r].IsPeriodic = IsPeriodic;
		RootParams[r].RhoPtr = &RhoPtr;
		RootParams[r].GradPtrs = &GradXYZPtrs;
		RootParams[r].HessPtrs = &HessPtrs;

		RootParams[r].BasisVectors = &BV;

		RootParams[r].HasGrad = GradXYZPtrs.size() == 3;
		RootParams[r].HasHess = HessPtrs.size() == 6;


		for (int i = 0; i < 3 && RootParams[r].HasGrad; ++i)
			RootParams[r].HasGrad = RootParams[r].GradPtrs->at(i).IsReady();
		for (int i = 0; i < 6 && RootParams[r].HasHess; ++i)
			RootParams[r].HasHess = RootParams[r].HessPtrs->at(i).IsReady();
	}

	vector<MultiRootObjects_s> MR(NumThreads);

	for (int m = 0; m < NumThreads; ++m){
		MR[m].Func = { &F3D, &DF3D, &FDF3D, 3, &RootParams[m] };
		MR[m].pos = gsl_vector_alloc(3);
		MR[m].T = gsl_multiroot_fdfsolver_newton;
		MR[m].s = gsl_multiroot_fdfsolver_alloc(MR[m].T, 3);
	}

	vector<vec3> TmpPoint(NumThreads), PrincDir(NumThreads), CellMinXYZ(NumThreads), CellMaxXYZ(NumThreads);
	vector<double> TmpRho(NumThreads);
	vector<char> TmpType(NumThreads);

	string StatusStr = "Finding critical points";

	vector<CritPoints_c> ThreadCPs(NumThreads);

	vector<int> NumPtsXYZ(3);
	vec3 ExtentXYZ = VolInfo.MaxXYZ - VolInfo.MinXYZ;
	for (int d = 0; d < 3; ++d) NumPtsXYZ[d] = ExtentXYZ[d] / CellSpacing;

	int NumThreadPts = NumPtsXYZ[2] / NumThreads;
	int ThreadPtNum = 0;

	vector<int> IJK(3);

	TecUtilDialogLaunchPercentDone(StatusStr.c_str(), TRUE);

#pragma omp parallel for schedule(dynamic)
	for (int zi = 0; zi < NumPtsXYZ[2]; ++zi){
		int ThreadNum = omp_get_thread_num();
		if (ThreadNum == 0 && !SetPercent(ThreadPtNum++, NumThreadPts, StatusStr, VolInfo.AddOnID)){
			IsOk = FALSE;
#pragma omp flush (IsOk)
		}
#pragma omp flush (IsOk)
		CellMinXYZ[ThreadNum][2] = VolInfo.MinXYZ[2] + static_cast<double>(zi)* CellSpacing;
		if (zi < NumPtsXYZ[2] - 1) CellMaxXYZ[ThreadNum][2] = CellMinXYZ[ThreadNum][2] + CellSpacing;
		else CellMaxXYZ[ThreadNum][2] = VolInfo.MaxXYZ[2];

		CellMinXYZ[ThreadNum][1] = VolInfo.MinXYZ[1];

		for (int yi = 0; yi < NumPtsXYZ[1] && IsOk; ++yi){
			if (yi < NumPtsXYZ[1] - 1) CellMaxXYZ[ThreadNum][1] = CellMinXYZ[ThreadNum][1] + CellSpacing;
			else CellMaxXYZ[ThreadNum][1] = VolInfo.MaxXYZ[1];

			CellMinXYZ[ThreadNum][0] = VolInfo.MinXYZ[0];

			for (int xi = 0; xi < NumPtsXYZ[0] && IsOk; ++xi){
				if (xi < NumPtsXYZ[0] - 1) CellMaxXYZ[ThreadNum][0] = CellMinXYZ[ThreadNum][0] + CellSpacing;
				else CellMaxXYZ[ThreadNum][0] = VolInfo.MaxXYZ[0];

				if (CritPointInCell(CellMinXYZ[ThreadNum], 
					CellMaxXYZ[ThreadNum], 
					TmpPoint[ThreadNum], 
					PrincDir[ThreadNum],
					TmpRho[ThreadNum], 
					RhoCutoff, 
					TmpType[ThreadNum], 
					RootParams[ThreadNum],
					MR[ThreadNum])){
					ThreadCPs[ThreadNum].AddPoint(TmpRho[ThreadNum], TmpPoint[ThreadNum], PrincDir[ThreadNum], TmpType[ThreadNum]);
				}

				CellMinXYZ[ThreadNum][0] += CellSpacing;
			}

			CellMinXYZ[ThreadNum][1] += CellSpacing;
		}
	}

	TecUtilDialogDropPercentDone();

	for (auto & m : MR)	if (m.s != NULL) gsl_multiroot_fdfsolver_free(m.s);
	// 	if (MR.pos != NULL)
	// 		gsl_vector_free(MR.pos);

	if (IsOk){
		CPs = CritPoints_c(ThreadCPs);
	}

	return IsOk;
}

