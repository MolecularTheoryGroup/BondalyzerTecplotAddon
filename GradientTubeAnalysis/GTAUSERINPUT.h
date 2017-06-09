#ifndef GTAUSERINPUT_H_
#define GTAUSERINPUT_H_

#include "CSM_DATA_TYPES.h"
#include <vector>

#include <armadillo>
using namespace arma;




#define MINREQPOINTS 2
#define MAXPOINTS 3

struct CSM_PlaneSelect_pa
{
	CSM_PlaneSelect_pa(Set_pa* SPtr = NULL, MouseButtonMode_e* MBPtr = NULL, std::vector<vec3>* PPtr = NULL);

	Set_pa* SetPtr;
	MouseButtonMode_e* MouseButtonModePtr;
	std::vector<vec3>* PointsPtr;
};

void ProbeTest_MenuCB(Set_pa *ActiveZones,
	MouseButtonMode_e *OldMouseMode,
	std::vector<vec3> *PlanePoints);
void STDCALL ProbeTest_ProbeCB(Boolean_t WasSuccessful,
								Boolean_t isNearestPoint,
								ArbParam_t ClientData);
Boolean_t GetActiveZonesSetScatterZones(Set_pa *ActiveZones);
void ActivateOldZones(Set_pa *ActiveZones);
void ClearCPList(std::vector<vec3> *PlanePoints);
void QuitProbing(Set_pa *ActiveZones,
	MouseButtonMode_e *OldMouseMode,
	std::vector<vec3> *PlanePoints,
	bool RestoreActiveZones);
void QuitButton_CB(Set_pa *ActiveZones,
	MouseButtonMode_e *OldMouseMode,
	std::vector<vec3> *PlanePoints,
	bool RestoreActiveZones);

void GTAPrepareGui();
void PopulateVarList();

#endif