#pragma once


#include <vector>
#include <string>

#include "CSM_DATA_TYPES.h"

using std::vector;
using std::string;

#include <armadillo>
using namespace arma;


struct VolExtentIndexWeights_s{
	vector<int> MaxIJK;
	vec3 MaxXYZ;
	vec3 MinXYZ;
	vec3 PointSpacingV123;

	mat33 BasisVectors;
	mat33 BasisNormalized;
	mat33 BasisInverse;
	vec3 BasisExtent;
	Boolean_t IsPeriodic;
	Boolean_t IsRectilinear;

	AddOn_pa AddOnID;

	int Index[8];
	double Weights[8];

	VolExtentIndexWeights_s(){ 
		MaxIJK.resize(3, -1); 
	}
	VolExtentIndexWeights_s(VolExtentIndexWeights_s const & rhs) { *this = rhs; }

	VolExtentIndexWeights_s & operator=(VolExtentIndexWeights_s const & rhs);
	Boolean_t operator==(VolExtentIndexWeights_s const & rhs) const;

	vec3 XYZ_to_Fractional(vec3 const & XYZ) const {
		return this->BasisInverse * (XYZ - this->MinXYZ);
	}

	vec3 Fractional_to_XYZ(vec3 const & UVW) const {
		return this->BasisVectors * UVW + this->MinXYZ;
	}

	bool PointIsInterior(vec3 const & r){
		vec3 CheckPt = this->XYZ_to_Fractional(r);
		for (int i = 0; i < 3; ++i){
			if (CheckPt[i] < 0.0 || CheckPt[i] > 1.0){
				return false;
			}
		}
		return true;
	}


	bool IsReady() const { return (MaxIJK[0] > 0 && MaxIJK[1] > 0 && MaxIJK[2] > 0); }
};
Boolean_t GetVolInfo(int VolZoneNum,
	vector<int> const & XYZVarNums,
	Boolean_t IsPeriodic,
	VolExtentIndexWeights_s & VolInfo);

Boolean_t InitializeVolInfo(vector<unsigned long long> const & MaxIJK,
	vec3 const & MinXYZ,
	vec3 const & MaxXYZ,
	mat33 const & BasisVectors,
	vec3 const & PointSpacingV123,
	VolExtentIndexWeights_s & VolInfo);

Boolean_t SetIndexAndWeightsForPoint(vec3 Point, VolExtentIndexWeights_s & SysInfo);
vector<int> GetIJKForPoint(vec3 & Point, VolExtentIndexWeights_s & VolZoneInfo);
void GetCellCornerIndices(int CornerNum, int & i, int & j, int & k);

double ValByCurrentIndexAndWeightsFromRawPtr(VolExtentIndexWeights_s const & VolZoneInfo, FieldDataPointer_c const & FDPtr);

double ValAtPointByPtr(vec3 & Point, VolExtentIndexWeights_s & VolZoneInfo, FieldDataPointer_c const & FDPtr);
