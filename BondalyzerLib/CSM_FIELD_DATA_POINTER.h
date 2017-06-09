#pragma once
#ifndef CSMFIELDDATAPOINTER_H_
#define CSMFIELDDATAPOINTER_H_

#include <vector>

using std::vector;


class FieldDataPointer_c{
public:
	FieldDataPointer_c(){}
	~FieldDataPointer_c(){ Close(); }

	const Boolean_t operator==(const FieldDataPointer_c & rhs) const;
	FieldDataPointer_c & operator=(const FieldDataPointer_c & rhs);
	const double operator[](const unsigned int & i) const;
	const Boolean_t Write(const unsigned int & i, const double & Val) const;
	const Boolean_t GetReadPtr(const int & ZoneNum, const int & VarNum);
	const Boolean_t GetWritePtr(const int & ZoneNum, const int & VarNum);
	void Close();

	const Boolean_t IsReady() const { return m_IsReady; }
	const unsigned int GetSize() const { return m_Size; }
	const vector<int> GetMaxIJK() const { return vector<int>(m_MaxIJK, m_MaxIJK + 3); }
	const int GetVarNum() const { return m_Var; }
	const int GetZoneNum() const { return m_Zone; }
	const FieldDataType_e GetDFType() const { return m_FDType; }
	const ValueLocation_e GetValueLocation() const { return m_ValueLocation; }
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
};


#endif // !CSMFIELDDATAPOINTER_H_
