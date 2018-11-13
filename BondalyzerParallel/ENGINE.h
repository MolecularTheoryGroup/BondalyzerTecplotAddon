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
#include "CSM_GRAD_PATH.h"
#include "CSM_GUI.h"

using std::vector;
using std::string;

enum BondalyzerCalcType_e{
	BondalyzerCalcType_CriticalPoints,
	BondalyzerCalcType_BondPaths,
	BondalyzerCalcType_RingLines,
	BondalyzerCalcType_CageNuclearPaths,
	BondalyzerCalcType_InteratomicSurfaces,
	BondalyzerCalcType_RingSurfaces,
	BondalyzerCalcType_BondBundleSurfaces,
	BondalyzerCalcType_RingBundleSurfaces,

	BondalyzerCalcType_GBA,

	BondalyzerCalcType_Batch,

	BondalyzerCalcType_Invalid
};

const static vector<string> BondalyzerStepGUITitles = {
	"critical points",
	"bond paths",
	"ring lines",
	"cage-nuclear paths",
	"interatomic surfaces",
	"ring surfaces",
	"bond bundle surfaces",
	"ring bundle surfaces",
	"gradient bundle analysis",
	"batch analysis"
};


void RefineActiveZones();
void GetClosedIsoSurfaceFromPoints();
void GetClosedIsoSurfaceFromNodes();
void GetAllClosedIsoSurfaces();
void ConnectCPsGetUserInfo();
void DrawEigenvectorArrowsGetUserInfo();
class FieldDataPointer_c;
void GetClosedIsoSurface(const int & IsoZoneNum, const std::vector<FieldDataPointer_c> & IsoReadPtrs, std::vector<int> & NodeNums);
void MakeSurfaceFromPathZonesGetUserInfo();
void MakeSliceFromPointSelectionGetUserInfo();

void GradientPathsOnSphereGetUserInfo();

void BondalyzerGetUserInfo(BondalyzerCalcType_e CalcType, const vector<GuiField_c> PassthroughFields = vector<GuiField_c>());

void VarNameFindReplaceGetUserInfo();
void ZoneNameFindReplaceGetUserInfo();

void GradientPathToolGetUserInfo();

const Boolean_t FindCritPoints(const int & VolZoneNum,
	const vector<int> & XYZVarNums,
	const int & RhoVarNum,
	const vector<int> & GradVarNums,
	const vector<int> & HessVarNums,
	const Boolean_t & IsPeriodic,
	const double & CellSpacing);

void DeleteCPsGetUserInfo();
void ExtractCPsGetUserInfo();
void CombineCPZonesGetUserInfo();

const Boolean_t BondalyzerBatch(const int & VolZoneNum,
	const vector<int> & CPZoneNums,
	const int & CPTypeVarNum,
	const vector<int> & XYZVarNums,
	const int & RhoVarNum,
	const vector<int> & GradVarNums,
	const vector<int> & HessVarNums,
	const Boolean_t & IsPeriodic,
	const int & RidgeFuncVarNum,
	const vector<bool> & CalcSteps);

const Boolean_t FindBondRingLines(const int & VolZoneNum,
	const vector<int> & OtherCPZoneNums,
	const int & SelectedCPZoneNum,
	vector<int> SelectedCPNums,
	const int & CPTypeVarNum,
	const CPType_e & CPType,
	const vector<int> & XYZVarNums,
	const int & RhoVarNum,
	const vector<int> & GradVarNums,
	const vector<int> & HessVarNums,
	const Boolean_t & IsPeriodic);

const Boolean_t FindCageNuclearPaths(const int & VolZoneNum,
	const vector<int> & OtherCPZoneNums,
	const int & SelectedCPZoneNum,
	vector<int> SelectedCPNums,
	const int & CPTypeVarNum,
	const vector<int> & XYZVarNums,
	const int & RhoVarNum,
	const vector<int> & GradVarNums,
	const vector<int> & HessVarNums,
	const Boolean_t & IsPeriodic);

const Boolean_t FindBondRingSurfaces(const int & VolZoneNum,
	const vector<int> & OtherCPZoneNums,
	const int & SelectedCPZoneNum,
	vector<int> SelectedCPNums,
	const int & CPTypeVarNum,
	const CPType_e & CPType,
	const vector<int> & XYZVarNums,
	const int & RhoVarNum,
	const vector<int> & GradVarNums,
	const vector<int> & HessVarNums,
	const int & RCSFuncVarNum,
	const Boolean_t & IsPeriodic);

void ExtractRadiusContourLinesToIOrderedPoints(const vector<int> & ZoneNums,
	const vec3 & Origin, 
	const double & Radius,
	const vector<int> & XYZVarNums);

void ExtractRSIntersectionsGetUserInfo();
void ExtractSurfaceSphereIntersections(
	const int & VolZoneNum,
	const int & CPTypeVarNum,
	const int & AllCPZoneNum,
	const vector<int> & CPList,
	const int & CPZoneNum,
	const double & Radius,
	const int & ResampleNumPoints,
	const vector<int> & XYZVarNums,
	vector<GradPath_c> & Intersections);

const Boolean_t GBA_Generation(
	const int & VolZoneNum,
	const int & CPZoneNum,
	const int & CPTypeVarNum,
	const vector<int> & XYZVarNums,
	const int & RhoVarNum,
	const vector<int> & GradVarNums,
	const vector<int> & HessVarNums,
	const Boolean_t & IsPeriodic,
	const vector<int> & SelectedCPNums,
	const double & SphereRadius,
	const int & RadiusMode,
	const int & RefinementLevel,
	const double & RhoCutoff,
	const int & NumGBEdgePoints,
	const int & GPNumPoints
	);

const Boolean_t CPNumbersMapBetweenZones(const int & AllCPsZoneNum,
	const int & SelectedCPZoneNum,
	const vector<int> & XYZVarNums,
	const vector<int> & SelectedCPNums,
	vector<int> & MappedCPNums);

void TestFunction();


#endif /* ENGINE_H_ */
