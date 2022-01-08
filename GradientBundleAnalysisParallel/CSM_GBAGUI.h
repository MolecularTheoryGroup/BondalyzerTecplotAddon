#pragma once

#include <string>
#include <vector>

using std::vector;
using std::string;


/*
*	Default values for GUI elements
*/
int const DefaultCP = 1;
string const DefaultIntVarStr = "Electron Density";
// double const DefaultRhoCutoff = 0.001;
double const GBADefaultRadialSphereApprxRadius = 0.2;
double const GBADefaultSeedRadius = 0.6;
double const GBADefaultHSeedRadius = 0.2;
int const GBADefaultSphereMeshRefinementLevel = 4;
int const GBAMinSphereRefinementLevel = 2;
int const GBAMaxSphereRefinementLevel = 7;
int const GBADefaultIntegrationPrecision = 3;
int const GBAMinSubdivisionTightness = 1;
int const GBAMaxSubdivisionTightness = 6;
int const GBADefaultSubdivisionTightness = 2;
int const DefaultNumGPPts = 300;
double const DefaultIBDist = 0.05;
double const DefaultIBAng = 20;
int const GBADefaultGBPerE = 200;
int const GBADefaultBPAngularGBs = 1;
int const GBADefaultMaxBPAngularGBs = 1;
int const GBADefaultMaxGBSubdivisionLevel = 0;
int const GBAMaxMaxSubdivisionLevel = 7;
int const GBAMaxSubdivisionLevel = 1;
int const GBADefaultMaxNumElems = 100000;
double const GBADefaultMinEdgeGPSpacing = 0.1;
double const GBADefaultMaxEdgeGPSpacing = 2.0;
double const GBADefaultEdgeGPSpacing = 0.5;
int const GBADefaultNumINSSmoothing = 1;
int const GBADefaultNumINSSmoothingHAtoms = 4;

ColorIndex_t const GBABondWedgeColor = Red_C;
ColorIndex_t const GBALoneWedgeColor = Blue_C;


extern std::map<string, int> NuclearNameToCPNum;

void ListDeselect(LgIndex_t ListID);
vector<string> ListGetSelectedStrings(LgIndex_t  ListID);
vector<int> ListGetSelectedItemNums(LgIndex_t ListID);
void ListPopulateWithVarNames(LgIndex_t ListID, int StartVarNum = 1);
vec3 GetCoordsFromListItem(LgIndex_t ItemIndex,
	LgIndex_t ListIndex,
	string* ItemFullString = nullptr,
	string* ItemNameString = nullptr,
	int* ItemNumber = nullptr,
	int* ZoneNumber = nullptr,
	bool* IsCP = nullptr,
	vector<int>* NumberOfCPs = nullptr);

