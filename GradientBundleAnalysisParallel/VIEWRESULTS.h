#pragma once

#include <string>
#include <vector>

using std::string;
using std::vector;


vector<string> const RankStrs = { "Nuclear", "Bond", "Ring", "Cage" };


void GBAResultViewerSelectSphere();
void GBAResultViewerSelectIntVar();
void GBAResultViewerSelectGB();
void GBAResultViewerSelectCondensedGBs();
void GBAResultViewerToggleSphere();
void GBAResultViewerDeleteSphere();
void GBAResultViewerActivateAllGB();

void GBAResultViewerPopulateGBs();

void SortCPNameList(vector<string> & StrList);

void ToggleFEVolumesProbeInstallCB();
void STDCALL ToggleFEVolumesProbeCB(Boolean_t WasSuccessful,
	Boolean_t isNearestPoint,
	ArbParam_t ClientData);
void STDCALL SelectGBsInRegionProbeCB(Boolean_t WasSuccessful,
	Boolean_t isNearestPoint,
	ArbParam_t ClientData);

void ExportGBAData();
bool GetSphereOrigin(int SphereZoneNum, vec3 & Origin);
void ResizeSpheres(double const & SizeFactor, Boolean_t AllSpheres, Boolean_t AbsoluteRadius);

void ResultsVarListReload();