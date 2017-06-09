
#include <cstdlib>
#include <sstream>
#include <string>
#include <vector>

#include "TECADDON.h"
#include "ADDGLBL.h"
#include "GUIDEFS.h"
#if !defined (MSWIN)
#include <unistd.h>
#endif
#include "ArgList.h"
#include "CSM_DATA_TYPES.h"
#include "ZONEVARINFO.h"

#include "gtauserinput.h"

#include <armadillo>
using namespace arma;





CSM_PlaneSelect_pa::CSM_PlaneSelect_pa(Set_pa* SPtr, MouseButtonMode_e* MBPtr, std::vector<vec3>* PPtr){
	SetPtr = SPtr;
	MouseButtonModePtr = MBPtr;
	PointsPtr = PPtr;
}

void ProbeTest_MenuCB(Set_pa* ActiveZones, 
						MouseButtonMode_e* OldMouseMode, 
						std::vector<vec3>* PlanePoints)
{
	TecUtilLockOn();
	if (TecUtilDataSetIsAvailable() && TecUtilFrameGetPlotType() != PlotType_Sketch){
		if (TecGUIListGetItemCount(MLCPs_MLST_S1) < 3){
			GetActiveZonesSetScatterZones(ActiveZones);
			*OldMouseMode = TecUtilMouseGetCurrentMode();
			if (TecGUIListGetItemCount(MLCPs_MLST_S1) == 0 && PlanePoints->size() > 0){
				PlanePoints->clear();
				PlanePoints->reserve(MAXPOINTS);
			}
			CSM_PlaneSelect_pa* PlaneClientData = new CSM_PlaneSelect_pa(ActiveZones, OldMouseMode, PlanePoints);
			ArgList_pa ProbeArgs = TecUtilArgListAlloc();
			TecUtilArgListAppendFunction(ProbeArgs, SV_CALLBACKFUNCTION, ProbeTest_ProbeCB);
			TecUtilArgListAppendString(ProbeArgs, SV_STATUSLINETEXT, "Select CPs to define a plane");
			TecUtilArgListAppendArbParam(ProbeArgs, SV_CLIENTDATA, reinterpret_cast<ArbParam_t>(PlaneClientData));
			//tecplot::toolbox::ArgList ProbeArgs;
			//ProbeArgs.appendFunction(SV_CALLBACKFUNCTION, ProbeTest_ProbeCB);
			//ProbeArgs.appendString(SV_STATUSLINETEXT, "Select CPs to define a plane");
			//ProbeArgs.appendArbParam(SV_CLIENTDATA, reinterpret_cast<ArbParam_t>(PlaneClientData));
			//if (!TecUtilProbeInstallCallbackX(ProbeArgs.getRef())){
			if (!TecUtilProbeInstallCallbackX(ProbeArgs)){
				TecUtilDialogErrMsg("Failed to install probe callback function!");
			}
			//TecUtilProbeInstallCallback(ProbeTest_ProbeCB,
			//	"Click to run the probe test.");
			TecUtilArgListDealloc(&ProbeArgs);
			TecGUIDialogLaunch(Dialog1Manager);
			TecGUIButtonSetText(PBSelectCPs_BTN_S1, "Done");
		}
		else{
			QuitProbing(ActiveZones, OldMouseMode, PlanePoints, true);
			std::stringstream ss;
			ss << "There are already " << MAXPOINTS << " selected. Clear selection and try again.";
			TecUtilDialogErrMsg(ss.str().c_str());
		}
	}
	else
		TecUtilDialogErrMsg("To execute this add-on Tecplot must have\n"
		"a data set and the frame must be in XY,\n"
		"2D, or 3D.");
	TecUtilLockOff();
}

void STDCALL ProbeTest_ProbeCB(Boolean_t WasSuccessful, 
								Boolean_t isNearestPoint,
								ArbParam_t ClientData){
	if (!TecGUIDialogIsUp(Dialog1Manager))
		TecGUIDialogLaunch(Dialog1Manager);

	if (!isNearestPoint){
		TecUtilDialogErrMsg("You didn't hold the control key when you clicked...");
		if (!TecGUIDialogIsUp(Dialog1Manager))
			TecGUIDialogLaunch(Dialog1Manager);
		return;
	}

	EntIndex_t ZoneNum = TecUtilProbeFieldGetZone();
	if (ZoneNum == TECUTILBADZONENUMBER){
		TecUtilDialogErrMsg("Failed to get zone number!");
		return;
	}

	TecUtilLockStart(AddOnID);

	CSM_PlaneSelect_pa* PlanePtr = reinterpret_cast<CSM_PlaneSelect_pa*>(ClientData);
	std::vector<vec3>* PlanePoints = PlanePtr->PointsPtr;
	MouseButtonMode_e* OldMouseMode = PlanePtr->MouseButtonModePtr;
	Set_pa* ActiveZones = PlanePtr->SetPtr; 

	EntIndex_t XYZVarNums[3] = { 0 };

	TecUtilAxisGetVarAssignments(&XYZVarNums[0], &XYZVarNums[1], &XYZVarNums[2]);
	std::vector<std::string> XYZCStr = { "X", "Y", "Z" };

	for (int i = 0; i < 3; ++i){
		char* TmpVarNameCStr = NULL;
		if (TecUtilVarGetName(XYZVarNums[i], &TmpVarNameCStr))
			if (strncmp(TmpVarNameCStr, XYZCStr[i].c_str(),XYZCStr[i].length()) != 0){
			TecUtilDialogErrMsg("Cannot use systems without XYZ variables.");
			QuitProbing(ActiveZones, OldMouseMode, PlanePoints, true);
			return;
			}
	}

	int NumZoneDims = 0;
	LgIndex_t IJKMax[3] = { 0 };
	TecUtilZoneGetIJK(ZoneNum, &IJKMax[0], &IJKMax[1], &IJKMax[2]);

	for (int i = 0; i < 3; ++i){
		if (XYZVarNums[i] == TECUTILBADVARNUMBER){
			TecUtilDialogErrMsg("Failed to get variable assignments!");
			QuitProbing(ActiveZones, OldMouseMode, PlanePoints, true);
			return;
		}
		if (IJKMax[i] > 1) NumZoneDims++;
	}
	if (NumZoneDims > 1){
		std::stringstream ss;
		ss << "Probed zone is " << NumZoneDims << "-dimensional. Should be 1-D!";
		TecUtilDialogErrMsg(ss.str().c_str());
		QuitProbing(ActiveZones, OldMouseMode, PlanePoints, true);
		return;
	}

	AuxData_pa AuxDataRef = TecUtilAuxDataZoneGetRef(ZoneNum);
	if (!VALID_REF(AuxDataRef)){
		TecUtilDialogErrMsg("Failed to get aux data ref!");
		QuitProbing(ActiveZones, OldMouseMode, PlanePoints, true);
		return;
	}
	std::string AuxNameBase = "CompChem.NumCrtPt";
	std::string CPTypeNames[4] = { "Atom", "Bond", "Ring", "Cage" };
	/*
		CP offsets are the indices of the first of a CP type.
		*/
	EntIndex_t CPOffsets[4];
	CPOffsets[0] = 1;
	EntIndex_t TotalOffset = 1;
	for (int i = 0; i < 3; ++i){
		std::stringstream ss;
		ss << AuxNameBase << CPTypeNames[i];
		char* TmpCStr = NULL;
		Boolean_t bJunk;
		TecUtilAuxDataGetStrItemByName(AuxDataRef, ss.str().c_str(), &TmpCStr, &bJunk);
		if (TmpCStr != NULL){
			TotalOffset += atoi(TmpCStr);
			CPOffsets[i + 1] = TotalOffset;
		}
		else{
			TecUtilDialogErrMsg("Failed to get CP offsets!");
			QuitProbing(ActiveZones, OldMouseMode, PlanePoints, true);
			return;
		}
		TecUtilStringDealloc(&TmpCStr);
	}

	EntIndex_t CPNum = TecUtilProbeGetPointIndex();
	EntIndex_t CPNumInType = 0;
	std::string CPType;
	for (int i = 0; i < 4; ++i){
		if (i == 3 || (CPNum >= CPOffsets[i] && CPNum < CPOffsets[i + 1])){
			CPType = CPTypeNames[i];
			CPNumInType = CPNum - CPOffsets[i] + 1;
			break;
		}
	}

	vec3	PointLoc;

	for (int i = 0; i < 3; ++i)
		PointLoc[i] = TecUtilProbeFieldGetValue(XYZVarNums[i]);

	bool IsDuplicate = false;
	for (int i = 0; i < PlanePoints->size() && !IsDuplicate; ++i){
		if (sum(PointLoc == PlanePoints->at(i)) == 3){
			TecUtilDialogErrMsg("Point already selected! Please select a different point.");
			IsDuplicate = true;
		}
	}
	
	if (!IsDuplicate){
		PlanePoints->push_back(PointLoc);

		std::stringstream ss;
		ss << CPType << " " << CPNumInType;
		TecGUIListAppendItem(MLCPs_MLST_S1, ss.str().c_str());
	}

	if (TecGUIListGetItemCount(MLCPs_MLST_S1) >= MAXPOINTS){
		QuitProbing(ActiveZones, OldMouseMode, PlanePoints, true);
#if(GTADEBUG)
		for (int i = 0; i < PlanePoints->size(); ++i){
			std::stringstream ss;
			ss << i+1 << ". { X Y Z } = {";
			for (int j = 0; j < 3; ++j)
				ss << " " << PlanePoints->at(i)[j] << " ";
			
			ss << "}";
			TecUtilDialogMessageBox(ss.str().c_str(), MessageBoxType_Information);
		}
#endif
	}

	TecUtilLockFinish(AddOnID);
}

Boolean_t GetActiveZonesSetScatterZones(Set_pa* ActiveZones){
	/*
		Get set of active zones
	*/

	if (ActiveZones != NULL)
		TecUtilSetDealloc(ActiveZones);
	TecUtilZoneGetActive(ActiveZones);
	Boolean_t IsOk = !TecUtilSetIsEmpty(*ActiveZones);

	if (IsOk){
		/*
			Turn off all zones
		*/

		//	Prepare arglist
		ArgList_pa CurrentArgList = TecUtilArgListAlloc();
		Set_pa EmptySet = TecUtilSetAlloc(TRUE);

		//	Set argument list for active field maps
		TecUtilArgListAppendString(CurrentArgList, SV_P1, SV_ACTIVEFIELDMAPS);
		TecUtilArgListAppendArbParam(CurrentArgList, SV_IVALUE, (ArbParam_t)EmptySet);
		SetValueReturnCode_e ReturnCode = TecUtilStyleSetLowLevelX(CurrentArgList);

		TecUtilArgListDealloc(&CurrentArgList);

		IsOk = (ReturnCode == SetValueReturnCode_Ok || ReturnCode == SetValueReturnCode_DuplicateValue);
	}

	if (IsOk){
		/*
			Get set of CP zones
		*/
		Set_pa CPZones = TecUtilSetAlloc(TRUE);
// 		std::string CheckName("Critical Points Zone 1");
		EntIndex_t NumZones = TecUtilDataSetGetNumZones();

		for (EntIndex_t i = 1; i < NumZones && IsOk; ++i){
			/*
			 *	Instead of finding the CP zone based on its name, which was having some 
			 *	issues for an unknown reason, will instead assume that the CP zone is
			 *	the first zone that is a 1-D i-ordered zone.
			 */
			LgIndex_t IJK[3];
			TecUtilZoneGetIJK(i, &IJK[0], &IJK[1], &IJK[2]);
			if (IJK[0] > 1 && IJK[1] <= 1 && IJK[2] <= 1){
				IsOk = TecUtilSetAddMember(CPZones, i, TRUE);
				break;
			}
// 			ZoneName_t CurrentName;
// 			IsOk = TecUtilZoneGetName(i, &CurrentName);
// 			if (IsOk && !strncmp(CheckName.c_str(), CurrentName, CheckName.length())){
// 				IsOk = TecUtilSetAddMember(CPZones, i, TRUE);
// 				break;
// 			}
		}

		if (IsOk && !TecUtilSetIsEmpty(CPZones)){
			/*
				Turn on CP zones and enable scatter
			*/

			//	Enable scatter
			TecUtilZoneSetScatter(SV_SHOW, CPZones, 0.0, TRUE);

			SetValueReturnCode_e ReturnCode = TecUtilZoneSetActive(CPZones, AssignOp_Equals);

			if ((ReturnCode != SetValueReturnCode_Ok && ReturnCode != SetValueReturnCode_DuplicateValue)){
				TecUtilDialogErrMsg("Failed to activate only CP zones!");
				IsOk = FALSE;
			}
		}

		TecUtilSetDealloc(&CPZones);
	}

	return IsOk;
}

void ActivateOldZones(Set_pa* ActiveZones){
	if (!TecUtilSetIsEmpty(*ActiveZones)){
		SetValueReturnCode_e ReturnCode = TecUtilZoneSetActive(*ActiveZones, AssignOp_Equals);
	}
}

void ClearCPList(std::vector<vec3>* PlanePoints){
	TecGUIListDeleteAllItems(MLCPs_MLST_S1);
	PlanePoints->clear();
	PlanePoints->reserve(MAXPOINTS);
}

void QuitProbing(Set_pa* ActiveZones, 
					MouseButtonMode_e* OldMouseMode, 
					std::vector<vec3>* PlanePoints,
					bool RestoreActiveZones){
	if (TecGUIDialogIsUp(Dialog1Manager)){
		TecGUIDialogDrop(Dialog1Manager);
	}
	TecGUIButtonSetText(PBSelectCPs_BTN_S1, "Select CPs");
	if (RestoreActiveZones)
		ActivateOldZones(ActiveZones);
	if (TecUtilMouseGetCurrentMode() == *OldMouseMode){
		TecUtilMouseSetMode(MouseButtonMode_RotateRollerBall);
		TecUtilMouseSetMode(MouseButtonMode_Probe);
	}
	else{
		TecUtilMouseSetMode(MouseButtonMode_RotateRollerBall);
		TecUtilMouseSetMode(*OldMouseMode);
	}
}

void QuitButton_CB(Set_pa* ActiveZones, 
					MouseButtonMode_e* OldMouseMode, 
					std::vector<vec3>* PlanePoints,
					bool RestoreActiveZones){
	QuitProbing(ActiveZones, OldMouseMode, PlanePoints, false);
	ClearCPList(PlanePoints);
	if (TecGUIDialogIsUp(Dialog1Manager)){
		PlanePoints->clear();
	}
	TecGUISidebarActivate(TECGUITECPLOTSIDEBAR);
}


void GTAPrepareGui(){
	TecGUIToggleSet(TGLDeleteZon_TOG_S1, FALSE);
	TecGUIToggleSet(TGLHideZones_TOG_S1, TRUE);
	TecGUIToggleSet(TGLShowOnly_TOG_S1, FALSE);
	TecGUIRadioBoxSetToggle(RBCutoffDi_RADIO_S1, (LgIndex_t)LessThan);
	TecGUITextFieldSetString(TFCutoff_TF_S1, "1e-3");
	TecGUITextFieldSetString(TFRadOffVal_TF_S1, "0.025");
	TecGUITextFieldSetString(TFSlcBegAng_TF_S1, "0");
	TecGUITextFieldSetString(TFSlcEndAng_TF_S1, "360");
	TecGUITextFieldSetString(TFStBegAng_TF_S1, "0");
	TecGUITextFieldSetString(TFStEndAng_TF_S1, "360");
	TecGUIOptionMenuSet(OptSliceAngl_OPT_S1, 17);
	TecGUIOptionMenuSet(OptStreamAng_OPT_S1, 8);

	PopulateVarList();
}


void PopulateVarList(){
	TecGUIListDeleteAllItems(MLSelectVar_MLST_S1);
	TecGUIOptionMenuDeleteAllItems(OPTCutoffVar_OPT_S1);
	if (TecUtilDataSetIsAvailable()){
		EntIndex_t NumVars = TecUtilDataSetGetNumVars();

		if (NumVars > 0){
			std::vector<std::string> RhoStrs = { "Electron Density", "Rho", "rho" };
			EntIndex_t RhoVarNum = VarNumByNameList(RhoStrs);
			for (int i = 1; i <= NumVars; ++i){
				char* VarName = NULL;
				if (TecUtilVarGetName(i, &VarName)){
					if (/*strncmp(VarName, "X", 1) && strncmp(VarName, "Y", 1) && strncmp(VarName, "Z", 1) &&*/ VarName != NULL)
					{
						TecGUIListAppendItem(MLSelectVar_MLST_S1, VarName);
						TecGUIOptionMenuAppendItem(OPTCutoffVar_OPT_S1, VarName);
					}
					TecUtilStringDealloc(&VarName);
				}
			}
			if (RhoVarNum > 0){
				TecGUIListSetSelectedItem(MLSelectVar_MLST_S1, RhoVarNum);
				TecGUIOptionMenuSet(OPTCutoffVar_OPT_S1, RhoVarNum);
			}
		}
		else{
			TecGUIListAppendItem(MLSelectVar_MLST_S1, "No Variables...");
			TecGUIOptionMenuAppendItem(OPTCutoffVar_OPT_S1, "No Variables...");
		}
	}
	else{
		TecGUIListAppendItem(MLSelectVar_MLST_S1, "No Data Set...");
		TecGUIOptionMenuAppendItem(OPTCutoffVar_OPT_S1, "No Data Set...");
	}
}

// void PopulateVarList(){
// 	TecGUIOptionMenuDeleteAllItems(OptSelectVar_OPT_S1);
// 	if (TecUtilDataSetIsAvailable()){
// 		EntIndex_t NumVars = TecUtilDataSetGetNumVars();
// 		if (NumVars > 0){
// 			std::vector<std::string> RhoStrs = { "Electron Density", "Rho", "rho" };
// 			EntIndex_t RhoVarNum = VarNumByNameList(RhoStrs);
// 			for (int i = 1; i <= NumVars; ++i){
// 				char* VarName = NULL;
// 				if (TecUtilVarGetName(i, &VarName)){
// 					if (VarName != NULL)
// 						TecGUIOptionMenuAppendItem(OptSelectVar_OPT_S1, VarName);
// 					TecUtilStringDealloc(&VarName);
// 				}
// 			}
// 			TecGUIOptionMenuSet(OptSelectVar_OPT_S1, RhoVarNum);
// 		}
// 		else
// 			TecGUIOptionMenuAppendItem(OptSelectVar_OPT_S1, "No variables...");
// 	}
// 	else
// 		TecGUIOptionMenuAppendItem(OptSelectVar_OPT_S1, "No data set...");
// }



