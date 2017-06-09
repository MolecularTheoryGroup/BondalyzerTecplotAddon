#include "TECADDON.h"
#include "ADDGLBL.h"
#if !defined (MSWIN)
#include <unistd.h>
#endif
#include "GUIDEFS.h"

#include <string>
#include <sstream>
#include <vector>

#include "ZONEVARINFO.h"
#include "MAINFUNCTION.h"
#include "VIEWRESULTS.h"

#include "GUICB.h"

using std::string;
using std::stringstream;
using std::to_string;
using std::vector;

/*
*	Default values for GUI elements
*/
const int DefaultCP = 1;
const string DefaultIntVarStr = "Electron Density";
const Boolean_t DefaultIntegrate = FALSE;
const Boolean_t DefaultSystemIsOpen = TRUE;
const double DefaultRhoCutoff = 0.001;
const double DefaultRadius = 0.2;
const int DefaultLevel = 3;
const int MinLevel = 0;
const int MaxLevel = 10;
const int DefaultSTPts = 100;
const double DefaultIBDist = 0.05;
const double DefaultIBAng = 20;

/*
 *	Results viewer sidebar functions
 */

/**
*/
static void Sidebar1Activate_CB(void)
{
	/*  <<< This function is called when sidebar "Sidebar" is activated >>> */
	TecUtilLockStart(AddOnID);

	GBAResultViewerPrepareGUI();

	TecUtilLockFinish(AddOnID);
}


/**
*/
static void Sidebar1Deactivate_CB(void)
{
	/*   <<< This function is called when sidebar "Sidebar" is deactivated >>> */
}


/**
*/
static void SLSelSphere_SLST_S1_CB(const LgIndex_t *I)
{
	TecUtilLockStart(AddOnID);
	TRACE1("Single selection list (SLSelSphere_SLST_S1) item selected,  Item is: %d\n", *I);
	GBAResultViewerSelectSphere();
	TecUtilLockFinish(AddOnID);
}


/**
*/
static void TGLSphereVis_TOG_S1_CB(const LgIndex_t *I)
{
	TecUtilLockStart(AddOnID);
	TRACE1("Toggle (TGLSphere_TOG_S1) Value Changed,  New value is: %d\n", *I);
	GBAResultViewerToggleSphere();
	TecUtilLockFinish(AddOnID);
}

/**
*/
static void BTNSphereDel_BTN_S1_CB(void)
{
	TecUtilLockStart(AddOnID);
	TRACE("Delete Button Pushed\n");
	GBAResultViewerDeleteSphere();
	TecUtilLockFinish(AddOnID);
}


/**
*/
static void MLSelGB_MLST_S1_CB(const LgIndex_t *I)
{
	TecUtilLockStart(AddOnID);
	TRACE1("Multi selection list (MLSelVol_MLST_S1) item selected,  First Item is: %d\n", *I);
	GBAResultViewerSelectGB();
	TecUtilLockFinish(AddOnID);
}


/**
*/
static void BTNTogMode_BTN_S1_CB(void)
{
	TecUtilLockStart(AddOnID);
	TRACE("      Toggle Mode Button Pushed\n");
	ToggleFEVolumesProbeInstallCB();
	TecUtilLockFinish(AddOnID);
}

/**
*/
static void BTNAllGB_BTN_S1_CB(void)
{
	TecUtilLockStart(AddOnID);
	TRACE("         Activate All Button Pushed\n");
	GBAResultViewerActivateAllGB();
	TecUtilLockFinish(AddOnID);
}


/**
*/
static void BTNExitRes_BTN_S1_CB(void)
{
	TecUtilLockStart(AddOnID);
	TRACE("Exit Results Viewer Button Pushed\n");
	TecGUISidebarActivate(TECGUITECPLOTSIDEBAR);
	TecUtilLockFinish(AddOnID);
}


/*
 *	Process system dialog functions
 */

/**
*/
static void Dialog1CloseButton_CB(void)
{
	TecUtilLockStart(AddOnID);
	GBAProcessSystemDeleteCPLabels();
	TecGUIDialogDrop(Dialog1Manager);
	TecUtilLockFinish(AddOnID);
}

/**
*/
static void BTNRun_BTN_D1_CB(void)
{
	TecUtilLockStart(AddOnID);
	TRACE("       Run\n Button Pushed\n");
	int NumSelected, *IJunk;
	TecGUIListGetSelectedItems(MLSelCPs_MLST_D1, &IJunk, &NumSelected);
	TecUtilArrayDealloc((void**)&IJunk);
	if (NumSelected > 0){
		GBAProcessSystemDeleteCPLabels();
		TecGUIDialogDrop(Dialog1Manager);
		MainFunction();
		GBAResultViewerPrepareGUI();
	}
	else{
		TecUtilDialogMessageBox("Please select a critical point.", MessageBoxType_Warning);
		if (!TecGUIDialogIsUp(Dialog1Manager)){
			TecGUIDialogLaunch(Dialog1Manager);
		}
	}
	TecUtilLockFinish(AddOnID);
}


/**
*/
static void Dialog1Init_CB(void)
{
	TecUtilLockStart(AddOnID);
	/* <<< Add init code (if necessary) here>>> */
	TecUtilLockFinish(AddOnID);
}


/**
*/
static void MLSelCPs_MLST_D1_CB(const LgIndex_t *I)
{
	TecUtilLockStart(AddOnID);
	TRACE1("Multi selection list (MLSelCPs_MLST_D1) item selected,  First Item is: %d\n", *I);
	GBAProcessSystemLabelSelectedCPs();
	if (!TecGUIDialogIsUp(Dialog1Manager))
		TecGUIDialogLaunch(Dialog1Manager);
	TecUtilLockFinish(AddOnID);
}


/**
*/
static void TGLInt_TOG_D1_CB(const LgIndex_t *I)
{
	TecUtilLockStart(AddOnID);
	TRACE1("Toggle (TGLInt_TOG_D1) Value Changed,  New value is: %d\n", *I);
	TecGUIListDeleteAllItems(MLSelVars_MLST_D1);
	if (TecGUIToggleGet(TGLInt_TOG_D1))
		GBAProcessSystemPopulateIntVarList();
	TecUtilLockFinish(AddOnID);
}


/**
*/
static void MLSelVars_MLST_D1_CB(const LgIndex_t *I)
{
	TecUtilLockStart(AddOnID);
	TRACE1("Multi selection list (MLSelVars_MLST_D1) item selected,  First Item is: %d\n", *I);
	TecUtilLockFinish(AddOnID);
}


/**
*/
static LgIndex_t  TFCutoff_TF_D1_CB(const char *S)
{
	LgIndex_t IsOk = 1;
	TecUtilLockStart(AddOnID);
	TRACE1("Text field (TFCutoff_TF_D1) Value Changed,  New value is: %s\n", S);
	TecUtilLockFinish(AddOnID);
	return (IsOk);
}


/**
*/
static void TGLOpenSys_TOG_D1_CB(const LgIndex_t *I)
{
	TecUtilLockStart(AddOnID);
	TRACE1("Toggle (TGLOpenSys_TOG_D1) Value Changed,  New value is: %d\n", *I);
	TecUtilLockFinish(AddOnID);
}


/**
*/
static LgIndex_t  TFLevel_TFS_D1_ValueChanged_CB(const char *S)
{
	LgIndex_t IsOk = 1;
	TecUtilLockStart(AddOnID);
	TRACE1("Spin Control Text field (TFLevel_TFS_D1) Value Changed,  New value is: %s\n", S);
	LgIndex_t Value;
	if (TecGUITextFieldGetLgIndex(TFLevel_TFS_D1, &Value)){
		if (Value < MinLevel)
			TecGUITextFieldSetString(TFLevel_TFS_D1, to_string(MinLevel).c_str());
		else if (Value > MaxLevel)
			TecGUITextFieldSetString(TFLevel_TFS_D1, to_string(MaxLevel).c_str());
	}
	else{
		TecGUITextFieldSetString(TFLevel_TFS_D1, to_string(DefaultLevel).c_str());
	}
	GBAProcessSystemUpdateNumTriangles();
	TecUtilLockFinish(AddOnID);
	return (IsOk);
}


/**
*/
static void TFLevel_TFS_D1_ButtonUp_CB(void)
{
	TecUtilLockStart(AddOnID);
	TRACE0("Spin control (TFLevel_TFS_D1) up callback called.\n");
	TecGUISpinTextFieldIncLgIndex(TFLevel_TFS_D1, 1, MinLevel, MaxLevel);
	GBAProcessSystemUpdateNumTriangles();
	TecUtilLockFinish(AddOnID);
}


/**
*/
static void TFLevel_TFS_D1_ButtonDown_CB(void)
{
	TecUtilLockStart(AddOnID);
	TRACE0("Spin control (TFLevel_TFS_D1) down callback called.\n");
	TecGUISpinTextFieldIncLgIndex(TFLevel_TFS_D1, -1, MinLevel, MaxLevel);
	GBAProcessSystemUpdateNumTriangles();
	TecUtilLockFinish(AddOnID);
}


/**
*/
static LgIndex_t  TFRad_TF_D1_CB(const char *S)
{
	LgIndex_t IsOk = 1;
	TecUtilLockStart(AddOnID);
	TRACE1("Text field (TFRad_TF_D1) Value Changed,  New value is: %s\n", S); 
	double Value;
	if (!TecGUITextFieldGetDouble(TFRad_TF_D1, &Value)){
		TecGUITextFieldSetString(TFRad_TF_D1, to_string(DefaultRadius).c_str());
	}
	TecUtilLockFinish(AddOnID);
	return (IsOk);
}


/**
*/
static LgIndex_t  TFIBDist_TF_D1_CB(const char *S)
{
	LgIndex_t IsOk = 1;
	TecUtilLockStart(AddOnID);
	TRACE1("Text field (TFIBDist_TF_D1) Value Changed,  New value is: %s\n", S);
	double Value;
	if (!TecGUITextFieldGetDouble(TFIBDist_TF_D1, &Value)){
		TecGUITextFieldSetString(TFIBDist_TF_D1, to_string(DefaultRadius).c_str());
	}
	TecUtilLockFinish(AddOnID);
	return (IsOk);
}


/**
*/
static LgIndex_t  TFIBAng_TF_D1_CB(const char *S)
{
	LgIndex_t IsOk = 1;
	TecUtilLockStart(AddOnID);
	TRACE1("Text field (TFIBAng_TF_D1) Value Changed,  New value is: %s\n", S); 
	double Value;
	if (!TecGUITextFieldGetDouble(TFIBAng_TF_D1, &Value)){
		TecGUITextFieldSetString(TFIBAng_TF_D1, to_string(DefaultRadius).c_str());
	}
	TecUtilLockFinish(AddOnID);
	return (IsOk);
}

/**
*/
static LgIndex_t  TFSTPts_TF_D1_CB(const char *S)
{
	LgIndex_t IsOk = 1;
	TecUtilLockStart(AddOnID);
	TRACE1("Text field (TFSTPts_TF_D1) Value Changed,  New value is: %s\n", S);
	int Value;
	Boolean_t ValRead = TecGUITextFieldGetLgIndex(TFSTPts_TF_D1, &Value);
	if (!ValRead){
		TecGUITextFieldSetString(TFSTPts_TF_D1, to_string(DefaultSTPts).c_str());
	}
	else if (Value % 2 != 0){
		TecGUITextFieldSetString(TFSTPts_TF_D1, to_string(--Value).c_str());
	}
	TecUtilLockFinish(AddOnID);
	return (IsOk);
}

/**
*/
static void RBRadMode_RADIO_D1_CB(const LgIndex_t *I)
{
	TecUtilLockStart(AddOnID);
	TRACE1("RadioBox (RBRadMode_RADIO_D1) Value Changed,  New value is: %d\n", *I);
	TecUtilLockFinish(AddOnID);
}

Boolean_t GBAProcessSystemPrepareGUI(){
	/*
	*	Prepare CP list
	*/
	TecGUIListDeleteAllItems(MLSelCPs_MLST_D1);
	EntIndex_t CPZoneNum = ZoneNumByName(string("Critical Points"));
	EntIndex_t CPVarNum = VarNumByName(string("CritPointType"));
	if (CPZoneNum < 0 || CPVarNum < 0){
		TecUtilDialogErrMsg("Couldn't find critical points zone or CP type variable.");
		return FALSE;
	}
	int Ranks[4] = { -3, -1, 1, 3 };
	vector<string> RankStrs = { "Atom", "Bond", "Ring", "Cage" };
	LgIndex_t IJK[3]; LgIndex_t NumCPs[4] = { 0, 0, 0 , 0};
	TecUtilZoneGetIJK(CPZoneNum, &IJK[0], &IJK[1], &IJK[2]);
	for (int i = 1; i <= IJK[0]; ++i){
		int CPType = TecUtilDataValueGetByZoneVar(CPZoneNum, CPVarNum, i);
		for (int j = 0; j < 4; ++j){
			if (CPType == Ranks[j]){
				NumCPs[j]++;
				if (CPType == -3 || CPType == 3){
					stringstream ss;
					int CPOffset = 0;
					for (int k = 0; k < j; ++k)
						CPOffset += NumCPs[k];
					ss << RankStrs[j] << " " << i - CPOffset;
					TecGUIListAppendItem(MLSelCPs_MLST_D1, ss.str().c_str());
				}
				break;
			}
		}
	}
	/*
	*	Prepare Var list
	*/
	TecGUIListDeleteAllItems(MLSelVars_MLST_D1);
	TecGUIToggleSet(TGLInt_TOG_D1, DefaultIntegrate);
	if (DefaultIntegrate){
		GBAProcessSystemPopulateIntVarList();
	}

	/*
	*	System boundary options
	*/
	TecGUIToggleSet(TGLOpenSys_TOG_D1, DefaultSystemIsOpen);
	TecGUITextFieldSetString(TFCutoff_TF_D1, to_string(DefaultRhoCutoff).c_str());

	/*
	*	Mesh parameters
	*/
	TecGUITextFieldSetString(TFRad_TF_D1, to_string(DefaultRadius).c_str());
	TecGUIRadioBoxSetToggle(RBRadMode_RADIO_D1, 2);
	TecGUITextFieldSetString(TFLevel_TFS_D1, to_string(DefaultLevel).c_str());
	GBAProcessSystemUpdateNumTriangles();
	TecGUITextFieldSetString(TFSTPts_TF_D1, to_string(DefaultSTPts).c_str());

	/*
	*	IB detection paraneters
	*/

	TecGUITextFieldSetString(TFIBDist_TF_D1, to_string(DefaultIBDist).c_str());
	TecGUITextFieldSetString(TFIBAng_TF_D1, to_string(DefaultIBAng).c_str());

	return TRUE;
}

void GBAProcessSystemPopulateIntVarList(){
	EntIndex_t NumVars = TecUtilDataSetGetNumVars();
	for (int i = 1; i <= NumVars; ++i){
		char * VarName;
		TecUtilVarGetName(i, &VarName);
		TecGUIListAppendItem(MLSelVars_MLST_D1, VarName);
		TecUtilStringDealloc(&VarName);
	}
}

void GBAProcessSystemUpdateNumTriangles(){
	LgIndex_t Level;
	TecGUITextFieldGetLgIndex(TFLevel_TFS_D1, &Level);
	LgIndex_t NumTriangles = 20 * pow(4, Level);
	stringstream ss;
	ss << NumTriangles << " Triangles";
	TecGUILabelSetText(LBLNumTri_LBL_D1, ss.str().c_str());
}

void GBAProcessSystemLabelSelectedCPs(){
	TecUtilLockStart(AddOnID);
	LgIndex_t * SelectedCPs;
	LgIndex_t NumSelected;
	TecGUIListGetSelectedItems(MLSelCPs_MLST_D1, &SelectedCPs, &NumSelected);
	if (NumSelected > 0){
		TecUtilDrawGraphics(FALSE);
		TecUtilPickDeselectAll();
		for (int i = 0; i < CPLabelIDs.size(); ++i){
			TecUtilPickText(CPLabelIDs[i]);
		}
		TecUtilPickClear();
		CPLabelIDs.clear();

		EntIndex_t CPZoneNum = ZoneNumByName(string("Critical Points"));
		EntIndex_t CPVarNum = VarNumByName(string("CritPointType"));
		if (CPZoneNum < 0 || CPVarNum < 0){
			TecUtilDialogErrMsg("Couldn't find critical points zone or CP type variable.");
			return;
		}
		int Ranks[4] = { -3, -1, 1, 3 };
		LgIndex_t IJK[3]; LgIndex_t NumCPs[4] = { 0, 0, 0, 0 };
		TecUtilZoneGetIJK(CPZoneNum, &IJK[0], &IJK[1], &IJK[2]);
		for (int i = 1; i <= IJK[0]; ++i){
			int CPType = TecUtilDataValueGetByZoneVar(CPZoneNum, CPVarNum, i);
			for (int j = 0; j < 4; ++j){
				if (CPType == Ranks[j]){
					NumCPs[j]++;
					break;
				}
			}
		}

		EntIndex_t XYZVarNums[3] = { -1, -1, -1 };
		TecUtilAxisGetVarAssignments(&XYZVarNums[0], &XYZVarNums[1], &XYZVarNums[2]);
		if (XYZVarNums[0] < 0 || XYZVarNums[1] < 0 || XYZVarNums[2] < 0){
			for (int i = 1; i <= 3; ++i){
				XYZVarNums[i - 1] = i;
			}
		}
		double TextOffset = 0.1;
		vector<string> RankStrs = { "Atom", "Bond", "Ring", "Cage" };
		for (int i = 0; i < NumSelected; ++i){
			char * CPName = TecGUIListGetString(MLSelCPs_MLST_D1, SelectedCPs[i]);
			int CPNum;
			string CPTypeStr;
			stringstream ss;
			ss << CPName;
			ss >> CPTypeStr >> CPNum;
			for (int i = 0; i < 4; ++i){
				if (CPTypeStr == RankStrs[i]){
					int RealCPNum = 0;
					for (int j = 0; j < i; ++j){
						RealCPNum += NumCPs[j];
					}
					CPNum = RealCPNum + CPNum;
					break;
				}
			}

			double CPXYZ[3];
			for (int j = 0; j < 3; ++j){
				CPXYZ[j] = TecUtilDataValueGetByZoneVar(CPZoneNum, XYZVarNums[j], CPNum);
			}

			TecUtilSetupTransformations();
			double XYZGrid[3];
			TecUtilConvert3DPositionToGrid(CPXYZ[0], CPXYZ[1], CPXYZ[2], &XYZGrid[0], &XYZGrid[1], &XYZGrid[2]);


			Text_ID LabelID = TecUtilTextCreate(CoordSys_Grid, XYZGrid[0] + TextOffset, XYZGrid[1] + TextOffset, Units_Frame, 2, CPName);
			TecUtilTextBoxSetType(LabelID, TextBox_Filled);
			TecUtilTextBoxSetFillColor(LabelID, White_C);

			Boolean_t IsValid = TecUtilTextIsValid(LabelID);

			CPLabelIDs.push_back(LabelID);

			TecUtilStringDealloc(&CPName);
		}
		TecUtilDrawGraphics(TRUE);
	}
	TecUtilArrayDealloc((void **)&SelectedCPs);
	TecUtilLockFinish(AddOnID);
}

void GBAProcessSystemDeleteCPLabels(){
	TecUtilLockStart(AddOnID);
	MouseButtonMode_e MouseMode = TecUtilMouseGetCurrentMode();
	if (CPLabelIDs.size() > 0){
		TecUtilDrawGraphics(FALSE);
		TecUtilPickDeselectAll();
		for (int i = 0; i < CPLabelIDs.size(); ++i){
			TecUtilPickText(CPLabelIDs[i]);
		}
		TecUtilPickClear();
		CPLabelIDs.clear();
		TecUtilDrawGraphics(TRUE);
		//TecUtilRedraw(TRUE);
	}
	TecUtilMouseSetMode(MouseMode);
	TecUtilLockFinish(AddOnID);
}



#include "guibld.cpp"
