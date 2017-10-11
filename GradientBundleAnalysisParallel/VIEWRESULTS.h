
#ifndef VIEWRESULTS_H_
#define VIEWRESULTS_H_

#include <string>
#include <vector>

using std::string;
using std::vector;


const vector<string> RankStrs = { "Nuclear", "Atom", "Bond", "Ring", "Cage" };


void GBAResultViewerSelectSphere();
void GBAResultViewerSelectIntVar();
void GBAResultViewerSelectGB();
void GBAResultViewerToggleSphere();
void GBAResultViewerDeleteSphere();
void GBAResultViewerActivateAllGB();

void SortCPNameList(vector<string> & StrList);

void ToggleFEVolumesProbeInstallCB();
void STDCALL ToggleFEVolumesProbeCB(Boolean_t WasSuccessful,
	Boolean_t isNearestPoint,
	ArbParam_t ClientData);

#endif