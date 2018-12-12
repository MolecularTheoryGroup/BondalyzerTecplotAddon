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

	Boolean_t operator==(FieldDataPointer_c const & rhs) const;
	FieldDataPointer_c & operator=(FieldDataPointer_c const & rhs);
	double operator[](unsigned int i) const;
	double At(vec3 & Pt, VolExtentIndexWeights_s & VolInfo) const;
	Boolean_t Write(unsigned int i, double const & Val) const;
	Boolean_t GetReadPtr(int ZoneNum, int VarNum);
	Boolean_t GetWritePtr(int ZoneNum, int VarNum);
	void Close();

	Boolean_t IsReady() const { return m_IsReady; }
	unsigned int Size() const { return m_Size; }
	vector<int> MaxIJK() const { return vector<int>(m_MaxIJK, m_MaxIJK + 3); }
	int VarNum() const { return m_Var; }
	int ZoneNum() const { return m_Zone; }
	ZoneType_e ZoneType() const { return m_ZoneType; }
	bool ZoneIsOrdered() const { return (m_ZoneType == ZoneType_Invalid); }
	FieldDataType_e FDType() const { return m_FDType; }
	ValueLocation_e ValueLocation() const { return m_ValueLocation; }
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

	Boolean_t operator==(FieldVecPointer_c const & rhs) const;
	FieldVecPointer_c & operator=(FieldVecPointer_c const & rhs);
	vec3 operator[](unsigned int i) const;
	Boolean_t Write(unsigned int i, vec3 const & Vec) const;
	Boolean_t GetReadPtr(int ZoneNum, vector<int> const & VarNums);
	Boolean_t GetWritePtr(int ZoneNum, vector<int> const & VarNums);
	void Close();

	Boolean_t IsReady() const { return Ptrs[0].IsReady(); }
	unsigned int Size() const { return Ptrs[0].Size(); }
	vector<int> MaxIJK() const { return Ptrs[0].MaxIJK(); }
	int VarNum() const { return Ptrs[0].VarNum(); }
	int ZoneNum() const { return Ptrs[0].ZoneNum(); }
	ZoneType_e ZoneType() const { return Ptrs[0].ZoneType(); }
	bool ZoneIsOrdered() const { return Ptrs[0].ZoneIsOrdered(); }
	ValueLocation_e ValueLocation() const { return Ptrs[0].ValueLocation(); }
private:
	FieldDataPointer_c Ptrs[3];
};


#endif // !CSMFIELDDATAPOINTER_H_
