
#ifndef VIEWRESULTS_H_
#define VIEWRESULTS_H_

#include <string>
#include <vector>

using std::string;
using std::vector;


const vector<string> RankStrs = { "Atom", "Bond", "Ring", "Cage" };
const string GBADataPrefix = "GBA.";
const string GBAElemNum = "Elem";
const vector<string> GBANodeNums = {
	"Node1",
	"Node2",
	"Node3"
};
const string GBASphereCPName = "SphereCP";
const string GBAZoneType = "ZoneType";
const string GBAZoneTypeSphereZone = "SphereMeshZone";
const string GBAZoneTypeFEVolumeZone = "FEVolumeZone";
const string GBAZoneTypeIBEdgeZone = "IBEdgeZone";
const string GBAVolumeCPName = "VolumeCP";

/*
*	Aux data names for bondalyzed files
*/
// For CP zone
const string CCDataNumCPs[4] = {
	"CompChem.NumCrtPtAtom",
	"CompChem.NumCrtPtBond",
	"CompChem.NumCrtPtRing",
	"CompChem.NumCrtPtCage"
};
// For gradiant path zone
const string CCDataGPEndNums[2] = {
	"CompChem.BegCrtPtNum",
	"CompChem.EndCrtPtNum"
};
const string CCDataGPEndTypes[2] = {
	"CompChem.BegCrtPtType",
	"CompChem.EndCrtPtType"
};

void GBAResultViewerPrepareGUI();
void GBAResultViewerSelectSphere();
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