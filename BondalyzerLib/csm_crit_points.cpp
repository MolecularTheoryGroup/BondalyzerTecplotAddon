#include <vector>
#include <string>

#include "omp.h"

#include "TECADDON.h"

#include "CSM_DATA_SET_INFO.h"
#include "CSM_VOL_EXTENT_INDEX_WEIGHTS.h"
#include "CSM_CALC_VARS.h"
#include "CSM_GRAD_PATH.h"

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
	const int & CPTypeVarNum,
	const int & RhoVarNum,
	MultiRootParams_s *MR) : CritPoints_c()
{
	int NumZones = TecUtilDataSetGetNumZones();
	int NumVars = TecUtilDataSetGetNumVars();

	REQUIRE(0 < CPZoneNum && CPZoneNum <= NumZones);
	REQUIRE(RhoVarNum <= NumVars);
	REQUIRE(0 < CPTypeVarNum && CPTypeVarNum <= NumVars);
	REQUIRE(XYZVarNums.size() == 3);
	for (const auto & i : XYZVarNums) REQUIRE(0 < i && i < NumVars);

	FieldDataPointer_c RhoPtr, CPTypePtr;
	vector<FieldDataPointer_c> XYZPtrs(3);

	TecUtilDataLoadBegin();

	if (RhoVarNum >= 0 && !RhoPtr.GetReadPtr(CPZoneNum, RhoVarNum)){
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

	for (int iCP = 0; iCP < CPTypePtr.Size(); ++iCP){
		int CPType = CPTypePtr[iCP];
		int ind = VectorGetElementNum(CPTypeList, (CPType_e)CPType);

		m_Rho[ind].push_back(RhoPtr.IsReady() ? RhoPtr[iCP] : 0.0);
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

const double CritPoints_c::GetMinCPDist(const vector<CPType_e> & CPTypes){
	if (m_MinCPDistFound && CPTypes == m_MinDistCPTypes)
		return m_MinCPDist;
	else if (FindMinCPDist(CPTypes))
		return m_MinCPDist;
	else
		return -1;
}


const double CritPoints_c::GetMinCPDist(const int & CPTypeInd, const int & CPOffset, const vector<CPType_e> & CPTypes){
	double MinDistSqr = DBL_MAX;

	vector<int> TypeIndList;
	for (const auto & i : CPTypes) TypeIndList.push_back(VectorGetElementNum(CPTypeList, i));

	for (const int & t : TypeIndList){
		for (const vec3 & pt : m_XYZ[t]){
			if (t != CPTypeInd || sum(GetXYZ(CPTypeInd, CPOffset) == pt) < 3){
				double TmpDistSqr = DistSqr(GetXYZ(CPTypeInd, CPOffset), pt);
				if (TmpDistSqr < MinDistSqr) MinDistSqr = TmpDistSqr;
			}
		}
	}

	if (MinDistSqr != DBL_MAX) return sqrt(MinDistSqr);
	else return -1;
}

const double CritPoints_c::GetMinCPDist(const int & CPTotOffset, const vector<CPType_e> & CPTypes){
	vector<int> TypeOffset = GetTypeNumOffsetFromTotOffset(CPTotOffset);

	return GetMinCPDist(TypeOffset[0], TypeOffset[1], CPTypes);
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

	RemoveSpuriousCPs();
}


/*
*	Mutators and other methods
*/

const Boolean_t CritPoints_c::FindMinCPDist(const vector<CPType_e> & CPTypes){
	Boolean_t IsOk = (m_TotNumCPs > 0);

	vector<int> TypeIndList;
	m_MinDistCPTypes = CPTypes;
	for (const auto & i : CPTypes) TypeIndList.push_back(VectorGetElementNum(CPTypeList, i));

	int MinI = -1, MinJ = -1, MinII = -1, MinJJ = -1;
	double TmpDbl;

	m_MinCPDist = DBL_MAX;

	if (IsOk){
		for (int i = 0; i < TypeIndList.size() - 1; ++i){
			for (int j = i; j < TypeIndList.size(); ++j){
				for (int ii = 0; ii < NumCPs(TypeIndList[i]); ++ii){
					for (int jj = 0; jj < NumCPs(TypeIndList[j]); ++jj){
						if (i != j || ii != jj){
							TmpDbl = DistSqr(GetXYZ(TypeIndList[i], ii), GetXYZ(TypeIndList[j], jj));
							if (TmpDbl < m_MinCPDist){
								m_MinCPDist = TmpDbl;
								MinI = i;
								MinII = ii;
								MinJ = j;
								MinJJ = jj;
							}
						}
					}
				}
			}
		}
		IsOk = (MinI >= 0 && MinJ >= 0 && MinII >= 0 && MinJJ >= 0);

		// 		for (int i = 0; i < m_TotNumCPs; ++i){
		// 			for (int j = i + 1; j < m_TotNumCPs; ++j){
		// 				TmpDbl = DistSqr(GetXYZ(i), GetXYZ(j));
		// 				if (TmpDbl < m_MinCPDist){
		// 					m_MinCPDist = TmpDbl;
		// 					MinI = i;
		// 					MinJ = j;
		// 				}
		// 			}
		// 		}
		// 		IsOk = (MinI >= 0 && MinJ >= 0);
	}

	if (IsOk){
		// 		m_MinCPDist = Distance(GetXYZ(MinI), GetXYZ(MinJ));
		m_MinCPDist = Distance(GetXYZ(TypeIndList[MinI], MinII), GetXYZ(TypeIndList[MinJ], MinJJ));
		m_MinCPDistFound = TRUE;
	}

	return IsOk;
}

void CritPoints_c::RemoveSpuriousCPs(const double & CheckDist){
	double CheckDistSqr = CheckDist * CheckDist;


	// Cycle through all CP types and use a distance (squared) criteria 
	// to check for duplicates.
	// CP types have priority according to type index, so if
	// a bond is too close to a nuclear CP, the bond CP is removed.
	// Duplicate CPs of same type will be replaced with a single CP
	// with the average (midpoint) values of the duplicates.
	// 
	// Do one pass to fix the spurious groups within each CP type
	for (int t = 0; t < 6; ++t){
		if (m_XYZ[t].size() > 0){
			// First do check within same type, replacing spurious 
			// groups with the "central" point for that group.

			// For maintaining a list of "groups" of spurious neighboring CPs that
			vector<vector<int> > DuplicateNeighbors(m_XYZ[t].size());

			// For marking each CP as duplicate or not
			vector<bool> IsSpurious(m_XYZ[t].size(), false);

			/*
			*	Each CP i gets checked against the [i+1, max].
			*	The first CP in a spurious neighbor group is the one that
			*	keeps the full list of CPs in the group.
			*	This is done by using the first position in DuplicateNeighbors[i]
			*	to point to the "parent" CP for that group, which for the first CP
			*	in the group is itself.
			*	All subsequent CPs in that group are then stored in the DuplicateNeighbors[]
			*	location indicated by its parent.
			*/

			for (int i = 0; i < m_XYZ[t].size(); ++i){
				if (!IsSpurious[i]) DuplicateNeighbors[i].push_back(i);
				for (int j = i + 1; j < m_XYZ[t].size(); ++j){
					if (DistSqr(m_XYZ[t][i], m_XYZ[t][j]) <= CheckDistSqr){
						IsSpurious[j] = true;
						DuplicateNeighbors[DuplicateNeighbors[i][0]].push_back(j);
						DuplicateNeighbors[j].push_back(DuplicateNeighbors[i][0]);
					}
				}
			}

			vector<vec3> NewXYZ, NewPD, NewEigVals;
			vector<mat33> NewEigVecs;
			vector<double> NewRho;
			for (int i = 0; i < m_XYZ[t].size(); ++i){
				if (!IsSpurious[i]){
					NewXYZ.push_back(zeros(3));
					if (i < m_PrincDir[t].size()) NewPD.push_back(zeros(3));
					if (i < m_EigVals[t].size()) NewEigVals.push_back(zeros(3));
					if (i < m_EigVecs[t].size()) NewEigVecs.push_back(zeros(3, 3));
					if (i < m_Rho[t].size()) NewRho.push_back(0);

					for (int & j : DuplicateNeighbors[i]){
						NewXYZ.back() += m_XYZ[t][j];
						if (j < m_PrincDir[t].size()) NewPD.back() += m_PrincDir[t][j];
						if (j < m_EigVals[t].size()) NewEigVals.back() += m_EigVals[t][j];
						if (j < m_EigVecs[t].size()) NewEigVecs.back() += m_EigVecs[t][j];
						if (j < m_Rho[t].size()) NewRho.back() += m_Rho[t][j];
					}

					if (DuplicateNeighbors[i].size() > 1){
						NewXYZ.back() /= (double)DuplicateNeighbors[i].size();
						if (NewPD.size() > 0) NewPD.back() /= (double)DuplicateNeighbors[i].size();
						if (NewEigVals.size() > 0) NewEigVals.back() /= (double)DuplicateNeighbors[i].size();
						if (NewEigVecs.size() > 0) NewEigVecs.back() /= (double)DuplicateNeighbors[i].size();
						if (NewRho.size() > 0) NewRho.back() /= (double)DuplicateNeighbors[i].size();
					}
				}
			}

			m_NumCPs[t] = NewXYZ.size();
			m_TotNumCPs += m_NumCPs[t];

			m_XYZ[t] = NewXYZ;
			m_PrincDir[t] = NewPD;
			m_EigVals[t] = NewEigVals;
			m_EigVecs[t] = NewEigVecs;
			m_Rho[t] = NewRho;
		}
	}

	// Now do a second pass where CPs are checked against the CPs of other types

	m_TotNumCPs = 0;

	// For marking each CP as duplicate or not
	vector<vector<bool> > IsSpurious(6);
	for (int t = 0; t < 6; ++t) IsSpurious[t].resize(m_XYZ[t].size(), false);

	for (int ti = 0; ti < 6; ++ti){
		if (m_XYZ[ti].size() > 0){
			for (int tj = ti + 1; tj < 6; ++tj){

				for (int i = 0; i < m_XYZ[ti].size(); ++i){
					if (!IsSpurious[ti][i]){
						for (int j = 0; j < m_XYZ[tj].size(); ++j){
							if (!IsSpurious[tj][j] && DistSqr(m_XYZ[ti][i], m_XYZ[tj][j]) <= CheckDistSqr){
								// Spurious CP found. The j cp is always removed because i cp has higher priority
								IsSpurious[tj][j] = true;
							}
						}
					}
				}
			}

			vector<vec3> NewXYZ, NewPD, NewEigVals;
			vector<mat33> NewEigVecs;
			vector<double> NewRho;
			for (int i = 0; i < m_XYZ[ti].size(); ++i){
				if (!IsSpurious[ti][i]){
					NewXYZ.push_back(m_XYZ[ti][i]);
					if (i < m_PrincDir[ti].size()) NewPD.push_back(m_PrincDir[ti][i]);
					if (i < m_EigVals[ti].size()) NewEigVals.push_back(m_EigVals[ti][i]);
					if (i < m_EigVecs[ti].size()) NewEigVecs.push_back(m_EigVecs[ti][i]);
					if (i < m_Rho[ti].size()) NewRho.push_back(m_Rho[ti][i]);
				}
			}

			m_NumCPs[ti] = NewXYZ.size();
			m_TotNumCPs += m_NumCPs[ti];

			m_XYZ[ti] = NewXYZ;
			m_PrincDir[ti] = NewPD;
			m_EigVals[ti] = NewEigVals;
			m_EigVecs[ti] = NewEigVecs;
			m_Rho[ti] = NewRho;
		}
	}



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

const CPType_e CritPoints_c::GetTypeFromTotOffset(const int & TotOffset) const{
	CPType_e Type = CPType_Invalid;
	if (TotOffset < 0) return Type;

	int CPCount = 0;
	for (int i = 0; i < 6; ++i){
		if (TotOffset < CPCount + m_NumCPs[i]){
			Type = CPTypeList[i];
			break;
		}
		CPCount += m_NumCPs[i];
	}
	return Type;
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
	int CPTypeVarNum = VarNumByName(CSMVarName.CritPointType);
	FieldDataPointer_c CPTypePtr;
	if (CPTypeVarNum <= 0){
		/*
		*	There is no CP type variable yet, which means there are also
		*	no CP zones, so the data type for all current zones can be
		*	bit for the CP type variable.
		*/
		vector<FieldDataType_e> DataTypes(TecUtilDataSetGetNumZones(), FieldDataType_Double);
		if (!TecUtilDataSetAddVar(CSMVarName.CritPointType.c_str(), DataTypes.data())){
			TecUtilDialogErrMsg("Failed to create CP type variable for CP zone");
			return{ -1 };
		}
		CPTypeVarNum = TecUtilDataSetGetNumVars();
	}
	vector<FieldDataType_e> DataTypes(TecUtilDataSetGetNumVars(), FieldDataType_Double);
// 	for (int i = 0; i < DataTypes.size(); ++i)
// 		DataTypes[i] = TecUtilDataValueGetType(1, i + 1);
// 	DataTypes[CPTypeVarNum - 1] = FieldDataType_Int16;

	if (!TecUtilDataSetAddZone(CSMZoneName.CriticalPoints.c_str(), NumCPs(), 1, 1, ZoneType_Ordered, DataTypes.data())){
		TecUtilDialogErrMsg("Failed to create CP zone");
		return{ -1 };
	}
	NewZoneNums.push_back(TecUtilDataSetGetNumZones());

	AuxDataZoneSetItem(NewZoneNums.back(), CSMAuxData.CC.ZoneType, CSMAuxData.CC.ZoneTypeCPs);
	AuxDataZoneSetItem(NewZoneNums.back(), CSMAuxData.CC.ZoneSubType, CSMAuxData.CC.ZoneTypeCPsAll);

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
	if (RhoVarNum > 0 && !RhoPtr.GetWritePtr(NewZoneNums.back(), RhoVarNum)){
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
			RhoPtr.Write(CPNum, m_Rho[t][ti]);
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
				if (!TecUtilDataSetAddZone(CSMZoneName.CPType[t].c_str(), NumCPs(t), 1, 1, ZoneType_Ordered, DataTypes.data())){
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
				if (RhoVarNum > 0 && !RhoPtr.GetWritePtr(NewZoneNums.back(), RhoVarNum)){
					TecUtilDialogErrMsg("Failed to get rho pointer for CP type zone");
					return{ -1 };
				}
				FieldDataPointer_c CPTypeZonePtr;
				if (!CPTypeZonePtr.GetWritePtr(NewZoneNums.back(), CPTypeVarNum)){
					TecUtilDialogErrMsg("Failed to get cp type pointer for CP type zone");
					return{ -1 };
				}
				for (int ti = 0; ti < NumCPs(t); ++ti){
					for (int d = 0; d < 3; ++d){
						XYZPtrs[d].Write(ti, GetXYZ(t, ti)[d]);
					}
					CPTypeZonePtr.Write(ti, CPTypeList[t]);
					RhoPtr.Write(ti, m_Rho[t][ti]);
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

void SetCPZone(const int & ZoneNum){
	REQUIRE(0 < ZoneNum && ZoneNum <= TecUtilDataSetGetNumZones());

	Set_pa TmpSet = TecUtilSetAlloc(TRUE);
	TecUtilSetAddMember(TmpSet, ZoneNum, TRUE);

	char* ZoneName;
	TecUtilZoneGetName(ZoneNum, &ZoneName);

	int IJK[3];
	TecUtilZoneGetIJK(ZoneNum, &IJK[0], &IJK[1], &IJK[2]);

	int CPTypeVarNum = VarNumByName(CSMVarName.CritPointType);
	if (CPTypeVarNum < 0){
		TecUtilDialogErrMsg("Failed to find CP type variable when setting CP zone information");
		return;
	}

	AuxDataZoneSetItem(ZoneNum, CSMAuxData.CC.ZoneType, CSMAuxData.CC.ZoneTypeCPs);

	bool CPTypeFound = false;
	for (int t = 0; t < 6; ++t){
		if (!CPTypeFound && CSMZoneName.CPType[t] == ZoneName){
			AuxDataZoneSetItem(ZoneNum, CSMAuxData.CC.ZoneSubType, CSMAuxData.CC.CPSubTypes[t]);
			TecUtilZoneSetScatter(SV_COLOR, TmpSet, 0.0, CPColorList[t]);
			CPTypeFound = true;
		}

		int NumCPs = 0;
		for (int i = 0; i < IJK[0]; ++i){
			if ((int)TecUtilDataValueGetByZoneVar(ZoneNum, CPTypeVarNum, i + 1) == (int)CPTypeList[t]) NumCPs++;
		}
		if (NumCPs > 0 || CPTypeFound){
			AuxDataZoneSetItem(ZoneNum, CSMAuxData.CC.NumCPs[t], to_string(NumCPs));
		}
		if (CPTypeFound) break;
	}
	if (!CPTypeFound){
		TecUtilZoneSetScatter(SV_COLOR, TmpSet, 0.0, Black_C);
		TecUtilZoneSetActive(TmpSet, AssignOp_MinusEquals);
	}

	TecUtilStateChanged(StateChange_ZonesAdded, (ArbParam_t)TmpSet);

	TecUtilZoneSetScatter(SV_SHOW, TmpSet, 0.0, TRUE);
	TecUtilZoneSetScatter(SV_FRAMESIZE, TmpSet, 1, FALSE);
	TecUtilZoneSetScatterSymbolShape(SV_GEOMSHAPE, TmpSet, GeomShape_Sphere);
	TecUtilZoneSetMesh(SV_SHOW, TmpSet, 0.0, FALSE);
	TecUtilZoneSetContour(SV_SHOW, TmpSet, 0.0, FALSE);
	TecUtilZoneSetShade(SV_SHOW, TmpSet, 0.0, FALSE);
	TecUtilZoneSetVector(SV_SHOW, TmpSet, 0.0, FALSE);

	TecUtilZoneSetActive(TmpSet, AssignOp_PlusEquals);

	TecUtilSetDealloc(&TmpSet);
	TecUtilStringDealloc(&ZoneName);
}

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
		CalcGradForPoint(Point, RootParams->VolInfo->DelXYZ, *RootParams->VolInfo, eye<mat>(3, 3), 0, RootParams->IsPeriodic, Grad, *RootParams->RhoPtr, GPType_Invalid, params);
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
				eye<mat>(3, 3),
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
			if (Type == CPType_Nuclear || Type == CPType_Ring)
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

	vec3 CheckPt;
	vec3 CellMinCheck = RootParams.VolInfo->BasisInverse * (CellMinXYZ - RootParams.VolInfo->MinXYZ),
		CellMaxCheck = RootParams.VolInfo->BasisInverse * (CellMaxXYZ - RootParams.VolInfo->MinXYZ);

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

			Status = gsl_multiroot_test_residual(MR.s->f, 1e-12);

			if (Iter > CheckPosIter || Iter >= MaxCPIter){
				CheckPt = RootParams.VolInfo->BasisInverse * (Point - RootParams.VolInfo->MinXYZ);
				CPInCell = sum(CheckPt >= CellMinCheck) == 3 && sum(CheckPt <= CellMaxCheck) == 3;
			}

		} while (CPInCell && Status == GSL_CONTINUE && Iter < MaxCPIter);

		if (CPInCell){
			Point = MR.s->x->data;

			CheckPt = RootParams.VolInfo->BasisInverse * (Point - RootParams.VolInfo->MinXYZ);
			CPInCell = sum(CheckPt >= CellMinCheck) == 3 && sum(CheckPt <= CellMaxCheck) == 3;
		}

		if (CPInCell){
			RhoValue = 0.0;
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
			if (Type == CPType_Nuclear || Type == CPType_Ring)
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

	mat33 I = eye<mat>(3, 3);
	RootParams.BasisVectors = &I;

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
		StatusLaunch(StatusStr.c_str(), VolInfo.AddOnID, TRUE);

	for (IJK[2] = StartIJK[2]; IJK[2] <= EndIJK[2] && IsOk; ++IJK[2]){
		kNum++;
		if (StartIJK[2] == 1 && !StatusUpdate(kNum-1, Range, StatusStr, VolInfo.AddOnID)){
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
		StatusDrop(VolInfo.AddOnID);

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
	double & RhoCutoff,
	const Boolean_t & IsPeriodic,
	FieldDataPointer_c & RhoPtr,
	vector<FieldDataPointer_c> & GradXYZPtrs,
	vector<FieldDataPointer_c> & HessPtrs)
{
	Boolean_t IsOk = ((GradXYZPtrs.size() == 3 || GradXYZPtrs.size() == 0)
		&& (HessPtrs.size() == 0 || HessPtrs.size() == 6));

	if (!IsOk) return IsOk;

	int NumThreads = omp_get_num_procs();

	vector<VolExtentIndexWeights_s> VolInfoList(NumThreads, VolInfo);

	// 	vector<vec3> BV(3);
	// 	BV[0] << 1 << 0 << 0;
	// 	BV[1] << 0 << 1 << 0;
	// 	BV[2] << 0 << 0 << 1;

	vector<MultiRootParams_s> RootParams(NumThreads);
	mat33 I = eye<mat>(3, 3);
	for (int r = 0; r < NumThreads; ++r){
		RootParams[r].CalcType = GPType_Classic;
		RootParams[r].VolInfo = &VolInfoList[r];
		RootParams[r].IsPeriodic = IsPeriodic;
		RootParams[r].RhoPtr = &RhoPtr;
		RootParams[r].GradPtrs = &GradXYZPtrs;
		RootParams[r].HessPtrs = &HessPtrs;

		RootParams[r].BasisVectors = &I;

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
	for (int d = 0; d < 3; ++d) NumPtsXYZ[d] = VolInfo.BasisExtent[d] / CellSpacing;

	vector<int> StartPt(3, 0), EndPt = NumPtsXYZ;

	if (!VolInfo.IsPeriodic){
		for (auto & i : StartPt) i += 1;
		for (auto & i : EndPt) i -= 1;
	}

	int NumThreadPts = NumPtsXYZ[2];
#ifndef _DEBUG
	NumThreadPts /= NumThreads;
#endif
	int ThreadPtNum = 0;

	StatusLaunch((StatusStr + " (Loading data...)").c_str(), VolInfo.AddOnID, TRUE);

	mat33 LatticeVector = VolInfo.BasisNormalized * CellSpacing;
	mat33 LatticeInverse = VolInfo.BasisNormalized * CellSpacing;

	cube RhoVals(NumPtsXYZ[0], NumPtsXYZ[1], NumPtsXYZ[2]);
#ifndef _DEBUG
#pragma omp parallel for schedule(dynamic)
#endif
	for (int zi = 0; zi < NumPtsXYZ[2]; ++zi){
		int ThreadNum = omp_get_thread_num();
		for (int yi = 0; yi < NumPtsXYZ[1]; ++yi){
			for (int xi = 0; xi < NumPtsXYZ[0]; ++xi){
				vec3 iXYZ;
				iXYZ << xi << yi << zi;
				vec3 Pt = VolInfo.MinXYZ + LatticeVector * iXYZ;
				double CellRho = RhoPtr.At(Pt, *RootParams[ThreadNum].VolInfo);
				RhoVals(xi, yi, zi) = CellRho;
			}
		}
	}

#ifndef _DEBUG
#pragma omp parallel for schedule(dynamic)
#endif
	for (int zi = StartPt[2]; zi < EndPt[2]; ++zi){
		int ThreadNum = omp_get_thread_num();
		if (ThreadNum == 0 && !StatusUpdate(ThreadPtNum++, NumThreadPts, StatusStr, VolInfo.AddOnID)){
			IsOk = FALSE;
#pragma omp flush (IsOk)
		}
#pragma omp flush (IsOk)
		// 		CellMinXYZ[ThreadNum][2] = VolInfo.MinXYZ[2] + static_cast<double>(zi)* CellSpacing;
		// 		if (zi < NumPtsXYZ[2] - 1) CellMaxXYZ[ThreadNum][2] = CellMinXYZ[ThreadNum][2] + CellSpacing;
		// 		else CellMaxXYZ[ThreadNum][2] = VolInfo.MaxXYZ[2];
		// 
		// 		CellMinXYZ[ThreadNum][1] = VolInfo.MinXYZ[1];

		for (int yi = StartPt[1]; yi < EndPt[1] && IsOk; ++yi){
			// 			if (yi < NumPtsXYZ[1] - 1) CellMaxXYZ[ThreadNum][1] = CellMinXYZ[ThreadNum][1] + CellSpacing;
			// 			else CellMaxXYZ[ThreadNum][1] = VolInfo.MaxXYZ[1];
			// 
			// 			CellMinXYZ[ThreadNum][0] = VolInfo.MinXYZ[0];

			for (int xi = StartPt[0]; xi < EndPt[0] && IsOk; ++xi){
				// 				if (xi < NumPtsXYZ[0] - 1) CellMaxXYZ[ThreadNum][0] = CellMinXYZ[ThreadNum][0] + CellSpacing;
				// 				else CellMaxXYZ[ThreadNum][0] = VolInfo.MaxXYZ[0];

				vec3 iXYZ;
				iXYZ << xi << yi << zi;
				CellMinXYZ[ThreadNum] = VolInfo.MinXYZ + LatticeVector * iXYZ;
				for (int d = 0; d < 3; ++d) {
					if (iXYZ[d] < NumPtsXYZ[d] - 1){
						CellMaxXYZ[ThreadNum][d] = CellMinXYZ[ThreadNum][d];
						for (int di = 0; di < 3; ++di) CellMaxXYZ[ThreadNum][d] += LatticeVector.at(d, di);
					}
					else{
						CellMaxXYZ[ThreadNum][d] = VolInfo.MinXYZ[d];
						for (int di = 0; di < 3; ++di) CellMaxXYZ[ThreadNum][d] += VolInfo.BasisVectors.at(d, di);
					}
				}

				/*
				*	Check for local Min/Max
				*/
				bool IsMaxMin;
				vector<double> Signs = { 1, -1 };
				// 				double CellRho = RhoPtr.At(CellMinXYZ[ThreadNum], *RootParams[ThreadNum].VolInfo);
				double CellRho = RhoVals(xi, yi, zi);
				for (int s = 0; s < 2; ++s){
					IsMaxMin = true;
					vec3 CompPt;
					for (int xj = xi - 1; xj <= xi + 1 && IsMaxMin; ++xj){
						int xk = xj;
						if (xk < 0){
							if (VolInfo.IsPeriodic) xk = NumPtsXYZ[0] - 1;
							else continue;
						}
						else if (xk >= NumPtsXYZ[0]){
							if (VolInfo.IsPeriodic) xk = 0;
							else continue;
						}
						for (int yj = yi - 1; yj <= yi + 1 && IsMaxMin; ++yj){
							int yk = yj;
							if (yk < 0){
								if (VolInfo.IsPeriodic) yk = NumPtsXYZ[1] - 1;
								else continue;
							}
							else if (yk >= NumPtsXYZ[1]){
								if (VolInfo.IsPeriodic) yk = 0;
								else continue;
							}
							for (int zj = zi - 1; zj <= zi + 1 && IsMaxMin; ++zj){
								int zk = zj;
								if (zk < 0){
									if (VolInfo.IsPeriodic) zk = NumPtsXYZ[2] - 1;
									else continue;
								}
								else if (zk >= NumPtsXYZ[2]){
									if (VolInfo.IsPeriodic) zk = 0;
									else continue;
								}
								if (xi != xk || yi != yk || zi != zk){
									// 									iXYZ << xj << yj << zj;
									// 									CompPt = VolInfo.MinXYZ + LatticeVector * iXYZ;
									// 									double CompVal = RhoPtr.At(CompPt, *RootParams[ThreadNum].VolInfo);
									double CompVal = RhoVals(xk, yk, zk);
									IsMaxMin = IsMaxMin && (Signs[s] * CellRho >= Signs[s] * CompVal);
								}
							}
						}
					}
					if (IsMaxMin){
						GradPath_c GP(
							CellMinXYZ[ThreadNum],
							(StreamDir_e)s,
							100,
							GPType_Classic,
							GPTerminate_AtRhoValue,
							NULL,
							&ThreadCPs[ThreadNum],
							NULL,
							&RhoCutoff,
							*RootParams[ThreadNum].VolInfo,
							*RootParams[ThreadNum].HessPtrs,
							*RootParams[ThreadNum].GradPtrs,
							*RootParams[ThreadNum].RhoPtr);

						// 						const vector<FieldDataPointer_c> * junkPtrs = NULL;
						// 						GradPath_c GP(
						// 							CellMinXYZ[ThreadNum],
						// 							(StreamDir_e)s,
						// 							100,
						// 							GPType_Classic,
						// 							GPTerminate_AtRhoValue,
						// 							NULL,
						// 							&ThreadCPs[ThreadNum],
						// 							NULL,
						// 							&RhoCutoff,
						// 							*RootParams[ThreadNum].VolInfo,
						// 							*junkPtrs,
						// 							*junkPtrs,
						// 							*RootParams[ThreadNum].RhoPtr);

						GP.Seed(false);

						if (GP.IsMade()){
#ifdef _DEBUG
							GP.SaveAsOrderedZone();
#endif

							CompPt = GP[-1];
							vec3 EigVals, PrincDir;
							mat33 EigVecs;
							char Type = 0;

							CalcEigenSystemForPoint(CompPt,
								EigVals,
								EigVecs,
								RootParams[ThreadNum]);

							for (int i = 0; i < 3; ++i){
								if (EigVals[i] > 0)
									Type++;
								else
									Type--;
							}
							if (Type == CPType_Nuclear || Type == CPType_Ring)
								PrincDir = EigVecs.row(0).t();
							else
								PrincDir = EigVecs.row(2).t();

							IsMaxMin = (Type == CPType_Nuclear || Type == CPType_Cage);

							// 						if (IsMaxMin)
							// 						ThreadCPs[ThreadNum].AddPoint(GP.RhoAt(-1), CompPt, PrincDir, Type);
							ThreadCPs[ThreadNum].AddPoint(GP.RhoAt(-1), CompPt, PrincDir, (s == 0 ? CPType_Nuclear : CPType_Cage));
						}

						break;
					}
				}
				/*
				*	Rigorous check using Newton-Raphson method
				*/
				if (!IsMaxMin && CritPointInCell(CellMinXYZ[ThreadNum],
					CellMaxXYZ[ThreadNum],
					TmpPoint[ThreadNum],
					PrincDir[ThreadNum],
					TmpRho[ThreadNum],
					RhoCutoff,
					TmpType[ThreadNum],
					RootParams[ThreadNum],
					MR[ThreadNum]))
				{
					ThreadCPs[ThreadNum].AddPoint(TmpRho[ThreadNum], TmpPoint[ThreadNum], PrincDir[ThreadNum], TmpType[ThreadNum]);
				}

				// 				CellMinXYZ[ThreadNum][0] += CellSpacing;
			}

			// 			CellMinXYZ[ThreadNum][1] += CellSpacing;
		}
	}

	StatusDrop(VolInfo.AddOnID);

	for (auto & m : MR)	if (m.s != NULL) gsl_multiroot_fdfsolver_free(m.s);
	// 	if (MR.pos != NULL)
	// 		gsl_vector_free(MR.pos);

	if (IsOk){
		CPs = CritPoints_c(ThreadCPs);
		CPs.RemoveSpuriousCPs(SpuriousCPDistanceRatioOfSearchGrid * CellSpacing);
	}

	return IsOk;
}

