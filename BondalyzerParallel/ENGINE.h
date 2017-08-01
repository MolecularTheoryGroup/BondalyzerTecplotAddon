/*
*****************************************************************
*****************************************************************
*******                                                  ********
*******    (C) Copyright 1988-2013 by TECPLOT INC.       *******
*******                                                  ********
*****************************************************************
*****************************************************************
*/

#ifndef  ENGINE_H_
#define ENGINE_H_ /* Only include once */

#include <fstream>
#include <vector>

#include "CSM_CRIT_POINTS.h"

using std::vector;
using std::string;

enum BondalyzerSteps_e{
	BondalyzerSteps_CriticalPoints,
	BondalyzerSteps_BondPaths,
	BondalyzerSteps_RingLines,
	BondalyzerSteps_InteratomicSurfaces,
	BondalyzerSteps_RingSurfaces,
	BondalyzerSteps_BondBundleSurfaces,
	BondalyzerSteps_RingBundleSurfaces,

	BondalyzerSteps_Invalid
};

const static vector<string> BondalyzerStepGUITitles = {
	"critical points",
	"bond paths",
	"ring lines",
	"interatomic surfaces",
	"ring surfaces",
	"bond bundle surfaces",
	"ring bundle surfaces"
};


void RefineActiveZones();
void GetClosedIsoSurfaceFromPoints();
void GetClosedIsoSurfaceFromNodes();
void GetAllClosedIsoSurfaces();
class FieldDataPointer_c;
void GetClosedIsoSurface(const int & IsoZoneNum, const std::vector<FieldDataPointer_c> & IsoReadPtrs, std::vector<int> & NodeNums);

void BondalyzerGetUserInfo(BondalyzerSteps_e CalcType);

const Boolean_t FindCritPoints(const int & VolZoneNum,
	const vector<int> & XYZVarNums,
	const int & RhoVarNum,
	const vector<int> & GradVarNums,
	const vector<int> & HessVarNums,
	const Boolean_t & IsPeriodic,
	const double & CellSpacing);

void DeleteCPsGetUserInfo();

const Boolean_t FindBondRingLines(const int & VolZoneNum,
	const int & AllCPsZoneNum,
	const vector<int> & SelectedCPNums,
	const int & CPTypeVarNum,
	const CPType_e & CPType,
	const vector<int> & XYZVarNums,
	const int & RhoVarNum,
	const vector<int> & GradVarNums,
	const vector<int> & HessVarNums,
	const Boolean_t & IsPeriodic);

const Boolean_t FindBondRingSurfaces(const int & VolZoneNum,
	const int & CPZoneNum,
	const vector<int> & SelectedCPNums,
	const int & CPTypeVarNum,
	const CPType_e & CPType,
	const vector<int> & XYZVarNums,
	const int & RhoVarNum,
	const vector<int> & GradVarNums,
	const vector<int> & HessVarNums,
	const Boolean_t & IsPeriodic);


#endif /* ENGINE_H_ */
