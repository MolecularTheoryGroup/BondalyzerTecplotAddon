#pragma once
#ifndef CSMVOLEXTENTINDEXWEIGHTS_H_
#define CSMVOLEXTENTINDEXWEIGHTS_H_


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
	vec3 DelXYZ;
	mat33 BasisVectors;
	mat33 BasisNormalized;
	mat33 BasisInverse;
	vec3 BasisExtent;
	Boolean_t IsPeriodic;

	AddOn_pa AddOnID;

	int Index[8];
	double Weights[8];

	VolExtentIndexWeights_s(){ 
		MaxIJK.resize(3); 
	}
	VolExtentIndexWeights_s & operator=(VolExtentIndexWeights_s const & rhs);
	Boolean_t operator==(VolExtentIndexWeights_s const & rhs) const;

	VolExtentIndexWeights_s(VolExtentIndexWeights_s const & rhs){ *this = rhs; }
};
Boolean_t GetVolInfo(int VolZoneNum,
	vector<int> const & XYZVarNums,
	Boolean_t IsPeriodic,
	VolExtentIndexWeights_s & VolInfo);

Boolean_t SetIndexAndWeightsForPoint(vec3 Point, VolExtentIndexWeights_s & SysInfo);
vector<int> GetIJKForPoint(vec3 & Point, VolExtentIndexWeights_s & VolZoneInfo);
void GetCellCornerIndices(int CornerNum, int & i, int & j, int & k);

double ValByCurrentIndexAndWeightsFromRawPtr(VolExtentIndexWeights_s const & VolZoneInfo, FieldDataPointer_c const & FDPtr);

double ValAtPointByPtr(vec3 & Point, VolExtentIndexWeights_s & VolZoneInfo, FieldDataPointer_c const & FDPtr);

#endif