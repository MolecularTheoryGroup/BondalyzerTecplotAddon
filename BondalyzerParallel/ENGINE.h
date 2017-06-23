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
#include<vector>

using std::vector;
using std::string;

enum BondalyzerSteps_e{
	CRITICALPOINTS,
	BONDPATHS,
	RINGLINES,
	INTERATOMICSURFACES,
	RINGSURFACES,
	BONDBUNDLESURFACES,
	RINGBUNDLESURFACES
};

const static vector<string> BondalyzerStepGUITitles = {
	"critical points",
	"bond paths",
	"ring lines",
	"interatomic surfaces",
	"ring surfaces"
};

const static int NumCircleGPs = 100;

void RefineActiveZones();
void GetClosedIsoSurfaceFromPoints();
void GetClosedIsoSurfaceFromNodes();
void GetAllClosedIsoSurfaces();
class FieldDataPointer_c;
void GetClosedIsoSurface(const int & IsoZoneNum, const std::vector<FieldDataPointer_c> & IsoReadPtrs, std::vector<int> & NodeNums);

void GetInfoFromUserForBondalyzer(BondalyzerSteps_e CalcType);

const Boolean_t FindCritPoints(const int & VolZoneNum,
	const vector<int> & XYZVarNums,
	const int & RhoVarNum,
	const vector<int> & GradVarNums,
	const vector<int> & HessVarNums,
	const Boolean_t & IsPeriodic);

const Boolean_t FindBondRingLines(const int & VolZoneNum,
	const int & CPZoneNum,
	const int & CPTypeVarNum,
	const char & CPType,
	const vector<int> & XYZVarNums,
	const int & RhoVarNum,
	const vector<int> & GradVarNums,
	const vector<int> & HessVarNums,
	const Boolean_t & IsPeriodic);

const Boolean_t FindBondRingSurfaces(const int & VolZoneNum,
	const int & CPZoneNum,
	const int & CPTypeVarNum,
	const char & CPType,
	const vector<int> & XYZVarNums,
	const int & RhoVarNum,
	const vector<int> & GradVarNums,
	const vector<int> & HessVarNums,
	const Boolean_t & IsPeriodic);


#endif /* ENGINE_H_ */
