#include "TECADDON.h"

#include "Set.h"
#include "ArgList.h"

#include <armadillo>
#include <string>

#include "CSM_FIELD_DATA_POINTER.h"
#include "CSM_VOL_EXTENT_INDEX_WEIGHTS.h"

using namespace arma;
using namespace tecplot::toolbox;

using std::to_string;

/*
*	Begin methods for FieldDataPointer_c
*/
Boolean_t FieldDataPointer_c::operator==(FieldDataPointer_c const & rhs) const{
	return (m_VoidPtr == rhs.m_VoidPtr

		&& m_ReadBitPtr == rhs.m_ReadBitPtr
		&& m_ReadBytePtr == rhs.m_ReadBytePtr
		&& m_ReadInt16Ptr == rhs.m_ReadInt16Ptr
		&& m_ReadInt32Ptr == rhs.m_ReadInt32Ptr
		&& m_ReadFltPtr == rhs.m_ReadFltPtr
		&& m_ReadDblPtr == rhs.m_ReadDblPtr

		&& m_WriteBitPtr == rhs.m_WriteBitPtr
		&& m_WriteBytePtr == rhs.m_WriteBytePtr
		&& m_WriteInt16Ptr == rhs.m_WriteInt16Ptr
		&& m_WriteInt32Ptr == rhs.m_WriteInt32Ptr
		&& m_WriteFltPtr == rhs.m_WriteFltPtr
		&& m_WriteDblPtr == rhs.m_WriteDblPtr

		&& m_Size == rhs.m_Size
		&& m_Zone == rhs.m_Zone
		&& m_Var == rhs.m_Var

		&& m_IsReady == rhs.m_IsReady
		&& m_IsReadPtr == rhs.m_IsReadPtr
		&& m_FDType == rhs.m_FDType);
}
FieldDataPointer_c & FieldDataPointer_c::operator=(FieldDataPointer_c const & rhs){
	if (this == &rhs) return *this;

	m_VoidPtr = rhs.m_VoidPtr;

	m_ReadBitPtr = rhs.m_ReadBitPtr;
	m_ReadBytePtr = rhs.m_ReadBytePtr;
	m_ReadInt16Ptr = rhs.m_ReadInt16Ptr;
	m_ReadInt32Ptr = rhs.m_ReadInt32Ptr;
	m_ReadFltPtr = rhs.m_ReadFltPtr;
	m_ReadDblPtr = rhs.m_ReadDblPtr;

	m_WriteBitPtr = rhs.m_WriteBitPtr;
	m_WriteBytePtr = rhs.m_WriteBytePtr;
	m_WriteInt16Ptr = rhs.m_WriteInt16Ptr;
	m_WriteInt32Ptr = rhs.m_WriteInt32Ptr;
	m_WriteFltPtr = rhs.m_WriteFltPtr;
	m_WriteDblPtr = rhs.m_WriteDblPtr;

	m_Size = rhs.m_Size;
	m_Zone = rhs.m_Zone;
	m_Var = rhs.m_Var;

	m_IsReadPtr = rhs.m_IsReadPtr;
	m_IsReady = rhs.m_IsReady;
	m_FDType = rhs.m_FDType;

	return *this;
}
double FieldDataPointer_c::operator[](unsigned int i) const{
	double Val = 0.0;
	REQUIRE(m_IsReady && i < m_Size && i >= 0);
	if (m_IsReadPtr){
		switch (m_FDType){
			case FieldDataType_Float:
				Val = static_cast<double>(m_ReadFltPtr[i]);
				break;
			case FieldDataType_Double:
				Val = static_cast<double>(m_ReadDblPtr[i]);
				break;
			case FieldDataType_Int16:
				Val = static_cast<double>(m_ReadInt16Ptr[i]);
				break;
			case FieldDataType_Int32:
				Val = static_cast<double>(m_ReadInt32Ptr[i]);
				break;
			case FieldDataType_Byte:
				Val = static_cast<double>(m_ReadBytePtr[i]);
				break;
			case FieldDataType_Bit:
				Val = static_cast<double>(m_ReadBitPtr[i]);
				break;
			default:
				break;
		}
	}
	else{
		switch (m_FDType){
			case FieldDataType_Float:
				Val = static_cast<double>(m_WriteFltPtr[i]);
				break;
			case FieldDataType_Double:
				Val = static_cast<double>(m_WriteDblPtr[i]);
				break;
			case FieldDataType_Int16:
				Val = static_cast<double>(m_WriteInt16Ptr[i]);
				break;
			case FieldDataType_Int32:
				Val = static_cast<double>(m_WriteInt32Ptr[i]);
				break;
			case FieldDataType_Byte:
				Val = static_cast<double>(m_WriteBytePtr[i]);
				break;
			case FieldDataType_Bit:
				Val = static_cast<double>(m_WriteBitPtr[i]);
				break;
			default:
				break;
		}
	}

	return Val;
}

double FieldDataPointer_c::At(vec3 & Pt, VolExtentIndexWeights_s & VolInfo) const{
	if (SetIndexAndWeightsForPoint(Pt, VolInfo)){
		return ValByCurrentIndexAndWeightsFromRawPtr(VolInfo, *this);
	}
	else{
		TecUtilDialogErrMsg("Failed to interpolate value with data pointer");
		return -1;
	}
}

Boolean_t FieldDataPointer_c::Write(unsigned int i, double const & Val) const{
	REQUIRE(m_IsReady && !m_IsReadPtr && i < m_Size && i >= 0);

	switch (m_FDType){
		case FieldDataType_Float:
			m_WriteFltPtr[i] = static_cast<float_t>(Val);
			break;
		case FieldDataType_Double:
			m_WriteDblPtr[i] = static_cast<double_t>(Val);
			break;
		case FieldDataType_Int16:
			m_WriteInt16Ptr[i] = static_cast<Int16_t>(Val);
			break;
		case FieldDataType_Int32:
			m_WriteInt32Ptr[i] = static_cast<Int32_t>(Val);
			break;
		case FieldDataType_Byte:
			m_WriteBytePtr[i] = static_cast<Byte_t>(Val);
			break;
		case FieldDataType_Bit:
			m_WriteBitPtr[i] = static_cast<bool>(Val);
			break;
		default:
			break;
	}

	return TRUE;
}
Boolean_t FieldDataPointer_c::GetReadPtr(int ZoneNum, int VarNum){
	m_IsReady = (
		ZoneNum >= 1 && ZoneNum <= TecUtilDataSetGetNumZones()
		&& VarNum >= 1 && VarNum <= TecUtilDataSetGetNumVars()
		);

	if (m_IsReady){
		TecUtilDataValueGetReadableRawPtr(ZoneNum, VarNum, &m_VoidPtr, &m_FDType);
		m_IsReady = (m_VoidPtr != nullptr && m_FDType != FieldDataType_Invalid);
	}
	if (m_IsReady){
		m_Zone = ZoneNum;
		m_Var = VarNum;
		m_IsReadPtr = TRUE;

		m_ValueLocation = TecUtilDataValueGetLocation(ZoneNum, VarNum);

		TecUtilZoneGetIJK(ZoneNum, &m_MaxIJK[0], &m_MaxIJK[1], &m_MaxIJK[2]);

		m_ZoneType = TecUtilZoneGetType(ZoneNum);

		if (m_ZoneType == ZoneType_Ordered){
			m_Size = m_MaxIJK[0] * m_MaxIJK[1] * m_MaxIJK[2];
		}
		else if (m_ZoneType != ZoneType_Invalid){
			if (m_ValueLocation == ValueLocation_Nodal)
				m_Size = m_MaxIJK[0];
			else
				m_Size = m_MaxIJK[1];
		}

		switch (m_FDType){
			case FieldDataType_Float:
				m_ReadFltPtr = reinterpret_cast<const float_t*>(m_VoidPtr);
				break;
			case FieldDataType_Double:
				m_ReadDblPtr = reinterpret_cast<const double_t*>(m_VoidPtr);
				break;
			case FieldDataType_Int16:
				m_ReadInt16Ptr = reinterpret_cast<const Int16_t*>(m_VoidPtr);
				break;
			case FieldDataType_Int32:
				m_ReadInt32Ptr = reinterpret_cast<const Int32_t*>(m_VoidPtr);
				break;
			case FieldDataType_Byte:
				m_ReadBytePtr = reinterpret_cast<const Byte_t*>(m_VoidPtr);
				break;
			case FieldDataType_Bit:
				m_ReadBitPtr = reinterpret_cast<const bool*>(m_VoidPtr);
				break;
			default:
				m_IsReady = FALSE;
				break;
		}
	}

// 	REQUIRE(m_IsReady);

	return m_IsReady;
}
Boolean_t FieldDataPointer_c::GetWritePtr(int ZoneNum, int VarNum){
	m_IsReady = (
		ZoneNum >= 1 && ZoneNum <= TecUtilDataSetGetNumZones()
		&& VarNum >= 1 && VarNum <= TecUtilDataSetGetNumVars()
		);

	if (m_IsReady){
		TecUtilDataValueGetWritableRawPtr(ZoneNum, VarNum, &m_VoidPtr, &m_FDType);
		m_IsReady = (m_VoidPtr != nullptr && m_FDType != FieldDataType_Invalid);
	}
	if (m_IsReady){
		m_Zone = ZoneNum;
		m_Var = VarNum;
		m_IsReadPtr = FALSE;

		m_ValueLocation = TecUtilDataValueGetLocation(ZoneNum, VarNum);

		TecUtilZoneGetIJK(ZoneNum, &m_MaxIJK[0], &m_MaxIJK[1], &m_MaxIJK[2]);

		m_ZoneType = TecUtilZoneGetType(ZoneNum);

		if (m_ZoneType == ZoneType_Ordered){
			m_Size = m_MaxIJK[0] * m_MaxIJK[1] * m_MaxIJK[2];
		}
		else if (m_ZoneType != ZoneType_Invalid){
			if (m_ValueLocation == ValueLocation_Nodal)
				m_Size = m_MaxIJK[0];
			else
				m_Size = m_MaxIJK[1];
		}

		switch (m_FDType){
			case FieldDataType_Float:
				m_WriteFltPtr = reinterpret_cast<float_t*>(m_VoidPtr);
				break;
			case FieldDataType_Double:
				m_WriteDblPtr = reinterpret_cast<double_t*>(m_VoidPtr);
				break;
			case FieldDataType_Int16:
				m_WriteInt16Ptr = reinterpret_cast<Int16_t*>(m_VoidPtr);
				break;
			case FieldDataType_Int32:
				m_WriteInt32Ptr = reinterpret_cast<Int32_t*>(m_VoidPtr);
				break;
			case FieldDataType_Byte:
				m_WriteBytePtr = reinterpret_cast<Byte_t*>(m_VoidPtr);
				break;
			case FieldDataType_Bit:
				m_WriteBitPtr = reinterpret_cast<bool*>(m_VoidPtr);
				break;
			default:
				m_IsReady = FALSE;
				break;
		}
	}

	REQUIRE(m_IsReady);

	return m_IsReady;
}

void FieldDataPointer_c::Close(){
	if (m_IsReady){
		if (!m_IsReadPtr && m_Zone > 0 && m_Var > 0){
			ArgList Args;
			Args.appendInt(SV_STATECHANGE, StateChange_VarsAltered);
			Args.appendSet(SV_ZONELIST, Set(m_Zone));
			Args.appendSet(SV_VARLIST, Set(m_Var));
			TecUtilStateChangedX(Args.getRef());
		}
		m_IsReady = FALSE;
	}
}


/*
 *	Begin methods for FieldVecPointer_c
 */


FieldVecPointer_c::FieldVecPointer_c(vector<FieldDataPointer_c> & InVecs){
	REQUIRE(InVecs.size() == 3);
	for (int i = 1; i < 3; ++i)
		REQUIRE(InVecs[i].Size() == InVecs[i - 1].Size());
	for (int i = 0; i < 3; ++i) Ptrs[i] = InVecs[i];
}

Boolean_t FieldVecPointer_c::operator == (FieldVecPointer_c const & rhs) const{
	for (int i = 0; i < 3; ++i) if (!(Ptrs[i] == rhs.Ptrs[i])) return FALSE;

	return TRUE;
}
FieldVecPointer_c & FieldVecPointer_c::operator=(FieldVecPointer_c const & rhs){
	if (rhs == *this)
		return *this;

	for (int i = 0; i < 3; ++i) Ptrs[i] = rhs.Ptrs[i];

	return *this;
}
vec3 FieldVecPointer_c::operator[](unsigned int i) const{
	vec3 a;
	for (int d = 0; d < 3; ++d) a[d] = Ptrs[d][i];

	return a;
}
Boolean_t FieldVecPointer_c::Write(unsigned int i, vec3 const & Vec) const{
	for (int d = 0; d < 3; ++d) if (!Ptrs[d].Write(i, Vec[d])) return FALSE;

	return TRUE;
}
Boolean_t FieldVecPointer_c::GetReadPtr(int ZoneNum, vector<int> const & VarNums){
	Boolean_t IsOk = TRUE;

	for (int i = 0; i < 3 && IsOk; ++i){
		IsOk = (Ptrs[i].GetReadPtr(ZoneNum, VarNums[i]) 
			&& (
			Ptrs[i].FDType() == FieldDataType_Double
						|| Ptrs[i].FDType() == FieldDataType_Float)
			);
		if (i > 0) IsOk = (IsOk && Ptrs[i].ValueLocation() == Ptrs[i - 1].ValueLocation()
			&& Ptrs[i].Size() == Ptrs[i - 1].Size());
	}

	if (!IsOk) Close();

	return IsOk;
}
Boolean_t FieldVecPointer_c::GetWritePtr(int ZoneNum, vector<int> const & VarNums){
	Boolean_t IsOk = TRUE;
	
	for (int i = 0; i < 3 && IsOk; ++i){
		IsOk = (Ptrs[i].GetWritePtr(ZoneNum, VarNums[i]) 
			&& (
			Ptrs[i].FDType() == FieldDataType_Double
						|| Ptrs[i].FDType() == FieldDataType_Float)
			);
		if (i > 0) IsOk = (IsOk && Ptrs[i].ValueLocation() == Ptrs[i - 1].ValueLocation());
	}

	if (!IsOk) Close();

	return IsOk;
}
void FieldVecPointer_c::Close(){
	for (int i = 0; i < 3; ++i) Ptrs[i].Close();
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
	for (int i = 0; i < GradVarNums.size() && IsOk; ++i) if (!GradPtrs[i].GetReadPtr(ZoneNum, GradVarNums[i])) {
		TecUtilDialogErrMsg(string("Failed to get gradient read pointer (zone " + to_string(ZoneNum) + ", var " + to_string(GradVarNums[i]) + ")").c_str());
		IsOk = FALSE;
	}
	for (int i = 0; i < HessVarNums.size() && IsOk; ++i) if (!HessPtrs[i].GetReadPtr(ZoneNum, HessVarNums[i])) {
		TecUtilDialogErrMsg(string("Failed to get hessian read pointer (zone " + to_string(ZoneNum) + ", var " + to_string(HessVarNums[i]) + ")").c_str());
		IsOk = FALSE;
	}

	TecUtilPleaseWait("Loading data", FALSE);

	return IsOk;
}