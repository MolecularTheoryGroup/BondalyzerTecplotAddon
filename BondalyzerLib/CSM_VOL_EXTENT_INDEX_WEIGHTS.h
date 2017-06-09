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
	Boolean_t IsPeriodic;

	AddOn_pa AddOnID;

	int Index[8];
	double Weights[8];

	VolExtentIndexWeights_s(){ MaxIJK.resize(3); }
	VolExtentIndexWeights_s & operator=(const VolExtentIndexWeights_s & rhs);
	const Boolean_t operator==(const VolExtentIndexWeights_s & rhs) const;

	VolExtentIndexWeights_s(const VolExtentIndexWeights_s & rhs){ *this = rhs; }
};

const Boolean_t SetIndexAndWeightsForPoint(vec3 & Point, VolExtentIndexWeights_s & SysInfo);
const vector<int> GetIJKForPoint(vec3 & Point, VolExtentIndexWeights_s & VolZoneInfo);
void GetCellCornerIndices(const int & CornerNum, int & i, int & j, int & k);

const double ValByCurrentIndexAndWeightsFromRawPtr(const VolExtentIndexWeights_s & VolZoneInfo, const FieldDataPointer_c & FDPtr);

#endif