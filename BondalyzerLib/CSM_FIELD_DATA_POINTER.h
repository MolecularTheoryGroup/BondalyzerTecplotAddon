#pragma once
#ifndef CSMFIELDDATAPOINTER_H_
#define CSMFIELDDATAPOINTER_H_

#include <vector>
#include <armadillo>

using std::vector;
using namespace arma;

struct VolExtentIndexWeights_s;

class FieldDataPointer_c{
public:
	FieldDataPointer_c(){}
	~FieldDataPointer_c(){ Close(); }

	const Boolean_t operator==(const FieldDataPointer_c & rhs) const;
	FieldDataPointer_c & operator=(const FieldDataPointer_c & rhs);
	const double operator[](const unsigned int & i) const;
	const double At(vec3 & Pt, VolExtentIndexWeights_s & VolInfo) const;
	const Boolean_t Write(const unsigned int & i, const double & Val) const;
	const Boolean_t GetReadPtr(const int & ZoneNum, const int & VarNum);
	const Boolean_t GetWritePtr(const int & ZoneNum, const int & VarNum);
	void Close();

	const Boolean_t IsReady() const { return m_IsReady; }
	const unsigned int Size() const { return m_Size; }
	const vector<int> MaxIJK() const { return vector<int>(m_MaxIJK, m_MaxIJK + 3); }
	const int VarNum() const { return m_Var; }
	const int ZoneNum() const { return m_Zone; }
	const ZoneType_e ZoneType() const { return m_ZoneType; }
	const bool ZoneIsOrdered() const { return (m_ZoneType == ZoneType_Invalid); }
	const FieldDataType_e FDType() const { return m_FDType; }
	const ValueLocation_e ValueLocation() const { return m_ValueLocation; }
private:
	void* m_VoidPtr = NULL;

	const bool* m_ReadBitPtr = NULL;
	const Byte_t* m_ReadBytePtr = NULL;
	const Int16_t* m_ReadInt16Ptr = NULL;
	const Int32_t* m_ReadInt32Ptr = NULL;
	const float_t* m_ReadFltPtr = NULL;
	const double_t* m_ReadDblPtr = NULL;

	bool* m_WriteBitPtr = NULL;
	Byte_t* m_WriteBytePtr = NULL;
	Int16_t* m_WriteInt16Ptr = NULL;
	Int32_t* m_WriteInt32Ptr = NULL;
	float_t* m_WriteFltPtr = NULL;
	double_t* m_WriteDblPtr = NULL;

	int m_MaxIJK[3];
	unsigned int m_Size = 0;

	int m_Zone = -1;
	int m_Var = -1;

	Boolean_t m_IsReady = FALSE;
	Boolean_t m_IsReadPtr;
	FieldDataType_e m_FDType = FieldDataType_Invalid;
	ValueLocation_e m_ValueLocation = ValueLocation_Invalid;
	ZoneType_e m_ZoneType = ZoneType_Invalid;
};

class FieldVecPointer_c{
public:
	FieldVecPointer_c(){}
	FieldVecPointer_c(vector<FieldDataPointer_c> & InVecs);
	~FieldVecPointer_c(){ Close(); }

	const Boolean_t operator==(const FieldVecPointer_c & rhs) const;
	FieldVecPointer_c & operator=(const FieldVecPointer_c & rhs);
	const vec3 operator[](const unsigned int & i) const;
	const Boolean_t Write(const unsigned int & i, const vec3 & Vec) const;
	const Boolean_t GetReadPtr(const int & ZoneNum, const vector<int> & VarNums);
	const Boolean_t GetWritePtr(const int & ZoneNum, const vector<int> & VarNums);
	void Close();

	const Boolean_t IsReady() const { return Ptrs[0].IsReady(); }
	const unsigned int Size() const { return Ptrs[0].Size(); }
	const vector<int> MaxIJK() const { return Ptrs[0].MaxIJK(); }
	const int VarNum() const { return Ptrs[0].VarNum(); }
	const int ZoneNum() const { return Ptrs[0].ZoneNum(); }
	const ZoneType_e ZoneType() const { return Ptrs[0].ZoneType(); }
	const bool ZoneIsOrdered() const { return Ptrs[0].ZoneIsOrdered(); }
	const ValueLocation_e ValueLocation() const { return Ptrs[0].ValueLocation(); }
private:
	FieldDataPointer_c Ptrs[3];
};


#endif // !CSMFIELDDATAPOINTER_H_
