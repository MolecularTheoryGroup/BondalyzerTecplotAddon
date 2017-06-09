#include "TECADDON.h"
#include "ADDGLBL.h"
#if !defined (MSWIN)
#include <unistd.h>
#endif
#include "GUIDEFS.h"

#include <string>
#include <sstream>
#include <vector>
#include <cmath>
#include <cctype>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <algorithm>

#include <armadillo>

#include "CSM_DATA_SET_INFO.h"
#include "CSM_DATA_TYPES.h"
#include "GBAENGINE.h"
#include "INTEGRATE.h"
#include "VIEWRESULTS.h"
#include "CSM_GUI.h"

#include "GUICB.h"

using std::string;
using std::stringstream;
using std::to_string;
using std::vector;
using std::ofstream;

using namespace arma;


vector<Text_ID> CPLabelIDs;

LgIndex_t GBATabNumber = 1;
bool GBAIsOpen = false;
string GBAOldNumContours = "30";
const int GBAMaxNumContours = 1000;
const int GBAMinNumContours = 5;

int GBALogLin = 1;
int GBACntSrc = 1;
int GBAResultViewerSelectIntVarNum = 1;
int GBAResultViewerSelectSphereNum = 1;


/*
 *	Integration sidebar functions
 */


/**
*/
static void TGLIntVolInt_TOG_T2_1_CB(const LgIndex_t *I)
{
	TecUtilLockStart(AddOnID);
	TRACE1("Toggle (TGLIntVolInt_TOG_T2_1) Value Changed,  New value is: %d\n", *I);
	TecUtilLockFinish(AddOnID);
}


/**
*/
static void MLIntSelVar_MLST_T2_1_CB(const LgIndex_t *I)
{
	TecUtilLockStart(AddOnID);
	TRACE1("Multi selection list (MLIntSelVar_MLST_T2_1) item selected,  First Item is: %d\n", *I);
	TecUtilLockFinish(AddOnID);
}


/**
*/
static void BTNExitInt_BTN_T2_1_CB(void)
{
	TecUtilLockStart(AddOnID);
	TRACE("Exit Results Viewer Button Pushed\n");
	TecGUISidebarActivate(TECGUITECPLOTSIDEBAR);
	TecUtilLockFinish(AddOnID);
}


/**
*/
static void MLIntSelSph_MLST_T2_1_CB(const LgIndex_t *I)
{
	TecUtilLockStart(AddOnID);
	TRACE1("Multi selection list (MLIntSelSph_MLST_T2_1) item selected,  First Item is: %d\n", *I);
	TecUtilLockFinish(AddOnID);
}


/**
*/
static void BTNIntegrate_BTN_T2_1_CB(void)
{
	TecUtilLockStart(AddOnID);
	TRACE("Integrate                   Button Pushed\n");
	PrepareIntegration(TRUE);
	TecUtilLockFinish(AddOnID);
}


/**
*/
static void SCIntPrecise_SC_T2_1_CB(const LgIndex_t *I)
{
	TecUtilLockStart(AddOnID);
	TRACE1("Scale (SCIntPrecise_SC_T2_1) Value Changed,  New value is: %d\n", *I);
	TecGUILabelSetText(LBLIntPrecis_LBL_T2_1, IntPrecisionLabels[*I].c_str());
	IntPrecise = *I;
	TecUtilLockFinish(AddOnID);
}


/**
*/
static void SCIntPrecise_SCD_T2_1_CB(const LgIndex_t *I)
{
	TecUtilLockStart(AddOnID);
	TRACE1("Scale (SCIntPrecise_SCD_T2_1) Value Changed on drag,  New value is: %d\n", *I);
	TecGUILabelSetText(LBLIntPrecis_LBL_T2_1, IntPrecisionLabels[*I].c_str());
	IntPrecise = *I;
	TecUtilLockFinish(AddOnID);
}

/*
 *	Results viewer sidebar functions
 */



/**
*/
static void SLSelSphere_SLST_T3_1_CB(const LgIndex_t *I)
{
	TecUtilLockStart(AddOnID);
	TRACE1("Single selection list (SLSelSphere_SLST_T3_1) item selected,  Item is: %d\n", *I);
	int NewGBAResultViewerSelectSphereNum = TecGUIListGetSelectedItem(SLSelSphere_SLST_T3_1);
	if (NewGBAResultViewerSelectSphereNum != GBAResultViewerSelectSphereNum){
		GBAResultViewerSelectSphereNum = NewGBAResultViewerSelectSphereNum;
		GBAResultViewerSelectSphere();
		if (GBACntSrc == 1)
			GBAResultViewerSelectIntVar();
		GBAReloadDialog();
	}
	TecUtilLockFinish(AddOnID);
}


/**
*/
static void TGLSphereVis_TOG_T3_1_CB(const LgIndex_t *I)
{
	TecUtilLockStart(AddOnID);
	TRACE1("Toggle (TGLSphere_TOG_T3_1) Value Changed,  New value is: %d\n", *I);
	GBAResultViewerToggleSphere();
	GBAReloadDialog();
	TecUtilLockFinish(AddOnID);
}

/**
*/
static void BTNSphereDel_BTN_T3_1_CB(void)
{
	TecUtilLockStart(AddOnID);
	TRACE("Delete Button Pushed\n");
	GBAResultViewerDeleteSphere();
	GBAReloadDialog();
	TecUtilLockFinish(AddOnID);
}


/**
*/
static void MLSelGB_MLST_T3_1_CB(const LgIndex_t *I)
{
	TecUtilLockStart(AddOnID);
	TRACE1("Multi selection list (MLSelVol_MLST_T3_1) item selected,  First Item is: %d\n", *I);
	GBAResultViewerSelectGB();
	GBAReloadDialog();
	TecUtilLockFinish(AddOnID);
}

/**
*/
static void SLSelVar_SLST_T3_1_CB(const LgIndex_t *I)
{
	TecUtilLockStart(AddOnID);
	TRACE1("Single selection list (SLSelVar_SLST_T3_1) item selected,  Item is: %d\n", *I);
	int NewGBAResultViewerSelectIntVarNum = TecGUIListGetSelectedItem(SLSelVar_SLST_T3_1);
	if (NewGBAResultViewerSelectIntVarNum != GBAResultViewerSelectIntVarNum){
		GBAResultViewerSelectIntVarNum = NewGBAResultViewerSelectIntVarNum;
		GBAResultViewerSelectIntVar();
		GBAReloadDialog();
	}
	TecUtilLockFinish(AddOnID);
}


/**
*/
static void BTNTogMode_BTN_T3_1_CB(void)
{
	TecUtilLockStart(AddOnID);
	TRACE("      Toggle Mode Button Pushed\n");
	ToggleFEVolumesProbeInstallCB();
	GBAReloadDialog();
	TecUtilLockFinish(AddOnID);
}

/**
*/
static void BTNAllGB_BTN_T3_1_CB(void)
{
	TecUtilLockStart(AddOnID);
	TRACE("         Activate All Button Pushed\n");
	GBAResultViewerActivateAllGB();
	GBAReloadDialog();
	TecUtilLockFinish(AddOnID);
}


/**
*/
static void BTNExitRes_BTN_T3_1_CB(void)
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
	GBAIsOpen = false;
	GBAProcessSystemDeleteCPLabels();
	TecGUIDialogDrop(Dialog1Manager);
	TecUtilLockFinish(AddOnID);
}

/**
*/
static void BTNRun_BTN_T1_1_CB(void)
{
	TecUtilLockStart(AddOnID);
	TRACE("       Run\n Button Pushed\n");
	int NumSelected, *IJunk;
	TecGUIListGetSelectedItems(MLSelCPs_MLST_T1_1, &IJunk, &NumSelected);
	TecUtilArrayDealloc((void**)&IJunk);
	if (NumSelected > 0){
		GBAProcessSystemDeleteCPLabels();
		TecGUIDialogDrop(Dialog1Manager);
		MainFunction();
		GBAResultViewerPrepareGUI();
		GBATabNumber = 3;
		if (TecGUIToggleGet(TGLInt_TOG_T1_1))
			PrepareIntegration(FALSE);
	}
	else{
		TecUtilDialogMessageBox("Please select a critical point.", MessageBoxType_Warning);
	}
	GBAReloadDialog();
	TecUtilLockFinish(AddOnID);
}


/**
*/
static void Dialog1Init_CB(void)
{
	TecUtilLockStart(AddOnID);
	GBAIsOpen = true;
	TecGUITabSetCurrentPage(TAB1_TB_D1, GBATabNumber);
	GBAProcessSystemPrepareGUI();
	/* <<< Add init code (if necessary) here>>> */
	TecUtilLockFinish(AddOnID);
}


/**
*/
static void MLSelCPs_MLST_T1_1_CB(const LgIndex_t *I)
{
	TecUtilLockStart(AddOnID);
	TRACE1("Multi selection list (MLSelCPs_MLST_T1_1) item selected,  First Item is: %d\n", *I);
	GBAProcessSystemLabelSelectedCPs();
	GBAReloadDialog();
	TecUtilLockFinish(AddOnID);
}


/**
*/
static void TGLInt_TOG_T1_1_CB(const LgIndex_t *I)
{
	TecUtilLockStart(AddOnID);
	TRACE1("Toggle (TGLInt_TOG_T1_1) Value Changed,  New value is: %d\n", *I);
	TecGUIListDeleteAllItems(MLSelVars_MLST_T1_1);
	if (TecGUIToggleGet(TGLInt_TOG_T1_1))
		ListPopulateWithVarNames(MLSelVars_MLST_T1_1);
	TecUtilLockFinish(AddOnID);
}


/**
*/
static void MLSelVars_MLST_T1_1_CB(const LgIndex_t *I)
{
	TecUtilLockStart(AddOnID);
	TRACE1("Multi selection list (MLSelVars_MLST_T1_1) item selected,  First Item is: %d\n", *I);
	TecUtilLockFinish(AddOnID);
}


/**
*/
static LgIndex_t  TFCutoff_TF_T1_1_CB(const char *S)
{
	LgIndex_t IsOk = 1;
	TecUtilLockStart(AddOnID);
	TRACE1("Text field (TFCutoff_TF_T1_1) Value Changed,  New value is: %s\n", S);
	TecUtilLockFinish(AddOnID);
	return (IsOk);
}


/**
*/
static void TGLOpenSys_TOG_T1_1_CB(const LgIndex_t *I)
{
	TecUtilLockStart(AddOnID);
	TRACE1("Toggle (TGLOpenSys_TOG_T1_1) Value Changed,  New value is: %d\n", *I);
	TecUtilLockFinish(AddOnID);
}


/**
*/
static LgIndex_t  TFLevel_TFS_T1_1_ValueChanged_CB(const char *S)
{
	LgIndex_t IsOk = 1;
	TecUtilLockStart(AddOnID);
	TRACE1("Spin Control Text field (TFLevel_TFS_T1_1) Value Changed,  New value is: %s\n", S);
	LgIndex_t Value;
	if (TecGUITextFieldGetLgIndex(TFLevel_TFS_T1_1, &Value)){
		if (Value < MinLevel)
			TecGUITextFieldSetString(TFLevel_TFS_T1_1, to_string(MinLevel).c_str());
		else if (Value > MaxLevel)
			TecGUITextFieldSetString(TFLevel_TFS_T1_1, to_string(MaxLevel).c_str());
	}
	else{
		TecGUITextFieldSetString(TFLevel_TFS_T1_1, to_string(DefaultLevel).c_str());
	}
	GBAProcessSystemUpdateNumTriangles();
	TecUtilLockFinish(AddOnID);
	return (IsOk);
}


/**
*/
static void TFLevel_TFS_T1_1_ButtonUp_CB(void)
{
	TecUtilLockStart(AddOnID);
	TRACE0("Spin control (TFLevel_TFS_T1_1) up callback called.\n");
	TecGUISpinTextFieldIncLgIndex(TFLevel_TFS_T1_1, 1, MinLevel, MaxLevel);
	GBAProcessSystemUpdateNumTriangles();
	TecUtilLockFinish(AddOnID);
}


/**
*/
static void TFLevel_TFS_T1_1_ButtonDown_CB(void)
{
	TecUtilLockStart(AddOnID);
	TRACE0("Spin control (TFLevel_TFS_T1_1) down callback called.\n");
	TecGUISpinTextFieldIncLgIndex(TFLevel_TFS_T1_1, -1, MinLevel, MaxLevel);
	GBAProcessSystemUpdateNumTriangles();
	TecUtilLockFinish(AddOnID);
}

/**
*/
static void SCPrecise_SC_T1_1_CB(const LgIndex_t *I)
{
	TecUtilLockStart(AddOnID);
	TRACE1("Scale (SCPrecise_SC_T1_1) Value Changed,  New value is: %d\n", *I);
	TecGUILabelSetText(LBLPrecise_LBL_T1_1, IntPrecisionLabels[*I].c_str());
	IntPrecise = *I;
	TecUtilLockFinish(AddOnID);
}


/**
*/
static void SCPrecise_SCD_T1_1_CB(const LgIndex_t *I)
{
	TecUtilLockStart(AddOnID);
	TRACE1("Scale (SCPrecise_SCD_T1_1) Value Changed on drag,  New value is: %d\n", *I);
	TecGUILabelSetText(LBLPrecise_LBL_T1_1, IntPrecisionLabels[*I].c_str());
	IntPrecise = *I;
	TecUtilLockFinish(AddOnID);
}


/**
*/
static void TGLVolInt_TOG_T1_1_CB(const LgIndex_t *I)
{
	TecUtilLockStart(AddOnID);
	TRACE1("Toggle (TGLVolInt_TOG_T1_1) Value Changed,  New value is: %d\n", *I);
	TecUtilLockFinish(AddOnID);
}


/**
*/
static LgIndex_t  TFRad_TF_T1_1_CB(const char *S)
{
	LgIndex_t IsOk = 1;
	TecUtilLockStart(AddOnID);
	TRACE1("Text field (TFRad_TF_T1_1) Value Changed,  New value is: %s\n", S); 
	double Value;
	if (!TecGUITextFieldGetDouble(TFRad_TF_T1_1, &Value)){
		TecGUITextFieldSetString(TFRad_TF_T1_1, to_string(DefaultRadius).c_str());
	}
	TecUtilLockFinish(AddOnID);
	return (IsOk);
}

/**
*/
static void SCNumEdgeGPs_SC_T1_1_CB(const LgIndex_t *I)
{
	TecUtilLockStart(AddOnID);
	TRACE1("Scale (SCNumEdgeGPs_SC_T1_1) Value Changed,  New value is: %d\n", *I);
	GBAProcessSystemUpdateNumGPsPerGB();
	TecUtilLockFinish(AddOnID);
}


/**
*/
static void SCNumEdgeGPs_SCD_T1_1_CB(const LgIndex_t *I)
{
	TecUtilLockStart(AddOnID);
	TRACE1("Scale (SCNumEdgeGPs_SCD_T1_1) Value Changed on drag,  New value is: %d\n", *I);
	GBAProcessSystemUpdateNumGPsPerGB();
	TecUtilLockFinish(AddOnID);
}


/**
*/
// static LgIndex_t  TFIBDist_TF_T1_1_CB(const char *S)
// {
// 	LgIndex_t IsOk = 1;
// 	TecUtilLockStart(AddOnID);
// 	TRACE1("Text field (TFIBDist_TF_T1_1) Value Changed,  New value is: %s\n", S);
// 	double Value;
// 	if (!TecGUITextFieldGetDouble(TFIBDist_TF_T1_1, &Value)){
// 		TecGUITextFieldSetString(TFIBDist_TF_T1_1, to_string(DefaultRadius).c_str());
// 	}
// 	TecUtilLockFinish(AddOnID);
// 	return (IsOk);
// }


/**
*/
// static LgIndex_t  TFIBAng_TF_T1_1_CB(const char *S)
// {
// 	LgIndex_t IsOk = 1;
// 	TecUtilLockStart(AddOnID);
// 	TRACE1("Text field (TFIBAng_TF_T1_1) Value Changed,  New value is: %s\n", S); 
// 	double Value;
// 	if (!TecGUITextFieldGetDouble(TFIBAng_TF_T1_1, &Value)){
// 		TecGUITextFieldSetString(TFIBAng_TF_T1_1, to_string(DefaultRadius).c_str());
// 	}
// 	TecUtilLockFinish(AddOnID);
// 	return (IsOk);
// }

/**
*/
static LgIndex_t  TFSTPts_TF_T1_1_CB(const char *S)
{
	LgIndex_t IsOk = 1;
	TecUtilLockStart(AddOnID);
	TRACE1("Text field (TFSTPts_TF_T1_1) Value Changed,  New value is: %s\n", S);
	int Value;
	Boolean_t ValRead = TecGUITextFieldGetLgIndex(TFSTPts_TF_T1_1, &Value);
	if (!ValRead){
		TecGUITextFieldSetString(TFSTPts_TF_T1_1, to_string(DefaultSTPts).c_str());
	}
	else if (Value % 2 != 0){
		TecGUITextFieldSetString(TFSTPts_TF_T1_1, to_string(--Value).c_str());
	}
	TecUtilLockFinish(AddOnID);
	return (IsOk);
}

/**
*/
static void RBRadMode_RADIO_T1_1_CB(const LgIndex_t *I)
{
	TecUtilLockStart(AddOnID);
	TRACE1("RadioBox (RBRadMode_RADIO_T1_1) Value Changed,  New value is: %d\n", *I);
	TecUtilLockFinish(AddOnID);
}

Boolean_t GBAProcessSystemPrepareGUI(){
	/*
	*	Prepare CP list
	*/
	int NumSelected = -1;
	int *SelectedNums;
	if (TecGUIListGetItemCount(MLSelCPs_MLST_T1_1) > 0)
		TecGUIListGetSelectedItems(MLSelCPs_MLST_T1_1, &SelectedNums, &NumSelected);
	TecGUIListDeleteAllItems(MLSelCPs_MLST_T1_1);

	/*
	*	Add any atom zones from the data import (e.g. H, C, O, Ag)
	*/

	for (int ZoneNum = 1; ZoneNum <= TecUtilDataSetGetNumZones(); ++ZoneNum){
		char* ZoneNameCStr;
		if (TecUtilZoneIsOrdered(ZoneNum) && TecUtilZoneGetName(ZoneNum, &ZoneNameCStr)){
			string ZoneName = ZoneNameCStr;
			TecUtilStringDealloc(&ZoneNameCStr);

			int ElemNum = SearchVectorForString(ElementSymbolList, ZoneName, false);
			if (ElemNum < 0)
				ElemNum = SearchVectorForString(ElementNameList, ZoneName);
			if (ElemNum >= 0){
				int NumPts;
				TecUtilZoneGetIJK(ZoneNum, &NumPts, NULL, NULL);
				for (int Counter = 1; Counter <= NumPts; ++Counter){
					TecGUIListAppendItem(MLSelCPs_MLST_T1_1, string(ZoneName + " " + to_string(Counter)).c_str());
				}
			}
		}
	}

	/*
	 *	Add CPs from Bondalyzer
	 */

	EntIndex_t CPZoneNum = ZoneNumByName(string("Critical Points"), true);
	EntIndex_t CPVarNum = VarNumByName(string("CritPointType"));
// 	if (CPZoneNum < 0 || CPVarNum < 0){
// 		TecUtilDialogErrMsg("Couldn't find critical points zone or CP type variable.");
// 		return FALSE;
// 	}
	if (CPZoneNum > 0 && CPVarNum > 0){
		int Ranks[4] = { -3, -1, 1, 3 };
		vector<string> RankStrs = { "Atom", "Bond", "Ring", "Cage" };
		LgIndex_t IJK[3]; LgIndex_t NumCPs[4] = { 0, 0, 0, 0 };
		TecUtilZoneGetIJK(CPZoneNum, &IJK[0], &IJK[1], &IJK[2]);
		for (int i = 1; i <= IJK[0]; ++i){
			int CPType = (int)TecUtilDataValueGetByZoneVar(CPZoneNum, CPVarNum, i);
			for (int j = 0; j < 4; ++j){
				if (CPType == Ranks[j]){
					NumCPs[j]++;
					if (CPType == -3 || CPType == 3){
						stringstream ss;
						int CPOffset = 0;
						for (int k = 0; k < j; ++k)
							CPOffset += NumCPs[k];
						ss << RankStrs[j] << " " << i - CPOffset;
						TecGUIListAppendItem(MLSelCPs_MLST_T1_1, ss.str().c_str());
					}
					break;
				}
			}
		}
	}

	if (NumSelected > 0)
		TecGUIListSetSelectedItems(MLSelCPs_MLST_T1_1, SelectedNums, NumSelected);



	/*
	*	Prepare Var list
	*/

	TecGUIScaleShowNumericDisplay(SCPrecise_SC_T1_1, TRUE);
	TecGUIScaleSetLimits(SCPrecise_SC_T1_1, 0, 4, 0);
	TecGUIScaleSetValue(SCPrecise_SC_T1_1, IntPrecise);
	TecGUILabelSetText(LBLPrecise_LBL_T1_1, IntPrecisionLabels[IntPrecise].c_str());
// 	TecGUIScaleSetValue(SCPrecise_SC_T1_1, GetRecommendedIntPrecision());

	TecGUIScaleShowNumericDisplay(SCNumEdgeGPs_SC_T1_1, FALSE);
	TecGUIScaleSetLimits(SCNumEdgeGPs_SC_T1_1, 1, 30, 0);
	TecGUIScaleSetValue(SCNumEdgeGPs_SC_T1_1, 4);
	GBAProcessSystemUpdateNumGPsPerGB();

	NumSelected = -1;
	if (TecGUIListGetItemCount(MLSelVars_MLST_T1_1) > 0)
		TecGUIListGetSelectedItems(MLSelVars_MLST_T1_1, &SelectedNums, &NumSelected);

	TecGUIListDeleteAllItems(MLSelVars_MLST_T1_1);
	TecGUIToggleSet(TGLInt_TOG_T1_1, DefaultIntegrate);
	TecGUIToggleSet(TGLVolInt_TOG_T1_1, DefaultVolIntegrate);
	if (DefaultIntegrate){
		ListPopulateWithVarNames(MLSelVars_MLST_T1_1);
	}
	if (NumSelected > 0)
		TecGUIListSetSelectedItems(MLSelVars_MLST_T1_1, SelectedNums, NumSelected);

	/*
	*	System boundary options
	*/
	TecGUIToggleSet(TGLOpenSys_TOG_T1_1, DefaultSystemIsOpen);
	TecGUITextFieldSetString(TFCutoff_TF_T1_1, to_string(DefaultRhoCutoff).c_str());

	/*
	*	Mesh parameters
	*/
	TecGUITextFieldSetString(TFRad_TF_T1_1, to_string(DefaultRadius).c_str());
	TecGUIRadioBoxSetToggle(RBRadMode_RADIO_T1_1, 2);
	TecGUITextFieldSetString(TFLevel_TFS_T1_1, to_string(DefaultLevel).c_str());
	GBAProcessSystemUpdateNumTriangles();
	TecGUITextFieldSetString(TFSTPts_TF_T1_1, to_string(DefaultSTPts).c_str());

	/*
	*	IB detection paraneters
	*/

// 	TecGUITextFieldSetString(TFIBDist_TF_T1_1, to_string(DefaultIBDist).c_str());
// 	TecGUITextFieldSetString(TFIBAng_TF_T1_1, to_string(DefaultIBAng).c_str());
// 	
	/*
	 *	Options in results tab
	 */

	TecGUITextFieldSetString(TFNumContours_TF_T3_1, GBAOldNumContours.c_str());
	TecGUIRadioBoxSetToggle(RBLogLin_RADIO_T3_1, GBALogLin);
	TecGUIRadioBoxSetToggle(RBCntSrc_RADIO_T3_1, GBACntSrc);
	TecGUIToggleSet(TGLExInt_TOG_T3_1, TRUE);
	TecGUIToggleSet(TGLExGBs_TOG_T3_1, TRUE);
	TecGUIToggleSet(TGLShowMesh_TOG_T3_1, FALSE);

	return TRUE;
}



void GBAProcessSystemUpdateNumTriangles(){
	LgIndex_t Level;
	TecGUITextFieldGetLgIndex(TFLevel_TFS_T1_1, &Level);
	LgIndex_t NumTriangles = 20 * (int)pow(4, Level);
	stringstream ss;
	ss << NumTriangles << " Triangles";
	TecGUILabelSetText(LBLNumTri_LBL_T1_1, ss.str().c_str());
}

void GBAProcessSystemUpdateNumGPsPerGB(){
	LgIndex_t Level;
	Level = TecGUIScaleGetValue(SCNumEdgeGPs_SC_T1_1);
	LgIndex_t NumGPs = 3 * (1 + Level);
	stringstream ss;
	ss << NumGPs;
	TecGUILabelSetText(LBLGPperGB_LBL_T1_1, ss.str().c_str());
}

void GBAProcessSystemLabelSelectedCPs(){
	TecUtilLockStart(AddOnID);
	LgIndex_t * SelectedCPs;
	LgIndex_t NumSelected;
	TecGUIListGetSelectedItems(MLSelCPs_MLST_T1_1, &SelectedCPs, &NumSelected);
	if (NumSelected > 0){
		TecUtilDrawGraphics(FALSE);
		TecUtilPickDeselectAll();
		for (int i = 0; i < CPLabelIDs.size(); ++i){
			TecUtilPickText(CPLabelIDs[i]);
		}
		TecUtilPickClear();
		CPLabelIDs.clear();

		double TextOffset = 0.1;
		vector<int> NumCPs;
		for (int i = 0; i < NumSelected; ++i){
			string CPName;
// 			
			vec3 CPXYZ = GetCoordsFromListItem(SelectedCPs[i], MLSelCPs_MLST_T1_1, &CPName, NULL, NULL, NULL, NULL, &NumCPs);
			vec3 XYZGrid;
			TecUtilConvert3DPositionToGrid(CPXYZ[0], CPXYZ[1], CPXYZ[2], &XYZGrid[0], &XYZGrid[1], &XYZGrid[2]);


			Text_ID LabelID = TecUtilTextCreate(CoordSys_Grid, XYZGrid[0] + TextOffset, XYZGrid[1] + TextOffset, Units_Frame, 2, CPName.c_str());
			TecUtilTextBoxSetType(LabelID, TextBox_Filled);
			TecUtilTextBoxSetFillColor(LabelID, White_C);

			Boolean_t IsValid = TecUtilTextIsValid(LabelID);

			CPLabelIDs.push_back(LabelID);
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
		TecUtilRedraw(TRUE);
	}
	TecUtilMouseSetMode(MouseMode);
	TecUtilLockFinish(AddOnID);
}

void PrepareIntegration(Boolean_t IntegratingFromIntTab){
	LgIndex_t AtomListID, VarListID, TGLID, ScaleID;
	if (IntegratingFromIntTab){
		AtomListID = MLIntSelSph_MLST_T2_1;
		VarListID = MLIntSelVar_MLST_T2_1;
		TGLID = TGLIntVolInt_TOG_T2_1;
		ScaleID = SCIntPrecise_SC_T2_1;
	}
	else{
		AtomListID = MLSelCPs_MLST_T1_1;
		VarListID = MLSelVars_MLST_T1_1;
		TGLID = TGLVolInt_TOG_T1_1;
		ScaleID = SCPrecise_SC_T1_1;
	}

	vector<string> AtomNameList = ListGetSelectedStrings(AtomListID);
	vector<string> IntVarNameList = ListGetSelectedStrings(VarListID);
	vector<int> IntVarNumList = ListGetSelectedItemNums(VarListID);

	PerformIntegration(AtomNameList, IntVarNameList, IntVarNumList,
		TecGUIToggleGet(TGLID), TecGUIScaleGetValue(ScaleID));
	/*for (int i = 1; i < 5; ++i){
		PerformIntegration(AtomNameList, IntVarNameList, IntVarNumList,
			TecGUIToggleGet(TGLID), i);
	}*/


	TecGUITabSetCurrentPage(TAB1_TB_D1, 3);
}

/**
*/
static void TAB1_TBA_D1_CB(const LgIndex_t *I)
{
	TecUtilLockStart(AddOnID);
	TRACE0("Activate callback for tab (TAB1_TBA_D1) called\n");
	GBATabNumber = *I;
	switch (*I){
		case 1:
			GBAProcessSystemPrepareGUI();
			break;
		case 2:
			GBAIntegrationPrepareGUI();
			break;
		case 3:
			GBAResultViewerPrepareGUI();
			break;
		default:
			break;
	}
	GBAReloadDialog();
	TecUtilLockFinish(AddOnID);
}


/**
*/
static void TAB1_TBD_D1_CB(const LgIndex_t *I)
{
	TecUtilLockStart(AddOnID);
	TRACE0("Deactivate callback for tab (TAB1_TBD_D1) called\n");
	TecUtilLockFinish(AddOnID);
}


void GBAReloadDialog(){
	TecUtilLockStart(AddOnID);
	if (GBAIsOpen && !TecGUIDialogIsUp(Dialog1Manager)){
		TecGUIDialogLaunch(Dialog1Manager);
	}
	TecUtilLockFinish(AddOnID);
}

/**
*/
static void RBCntSrc_RADIO_T3_1_CB(const LgIndex_t *I)
{
	TecUtilLockStart(AddOnID);
	TRACE1("RadioBox (RBCntSrc_RADIO_T3_1) Value Changed,  New value is: %d\n", *I);
	GBACntSrc = TecGUIRadioBoxGetToggle(RBCntSrc_RADIO_T3_1);
	GBAResultViewerSelectIntVar();
	TecUtilLockFinish(AddOnID);
}

/**
*/
static void RBLogLin_RADIO_T3_1_CB(const LgIndex_t *I)
{
	TecUtilLockStart(AddOnID);
	TRACE1("RadioBox (RBLogLin_RADIO_T3_1) Value Changed,  New value is: %d\n", *I);
	GBALogLin = TecGUIRadioBoxGetToggle(RBLogLin_RADIO_T3_1);
	GBAResultViewerSelectIntVar();
	TecUtilLockFinish(AddOnID);
}

/**
*/
static LgIndex_t  TFNumContours_TF_T3_1_CB(const char *S)
{
	LgIndex_t IsOk = 1;
	TecUtilLockStart(AddOnID);
	string TFStr = S;
	bool IsNumber = true;
	for (const char & i : TFStr){
		IsNumber = (std::isdigit(i));
		if (!IsNumber) break;
	}
	if (IsNumber){
		GBAOldNumContours = to_string(MAX(GBAMinNumContours, MIN(GBAMaxNumContours, stoi(TFStr))));
	}
	TecGUITextFieldSetString(TFNumContours_TF_T3_1, GBAOldNumContours.c_str());
	TRACE1("Text field (TFNumContours_TF_T3_1) Value Changed,  New value is: %s\n", S);
	TecUtilLockFinish(AddOnID);
	return (IsOk);
}

/*
*	Result viewer sidebar functions
*/
void GBAResultViewerPrepareGUI(){

	TecUtilLockStart(AddOnID);

	/*
	*	Clear lists and stuff
	*/

	TecGUIListDeleteAllItems(SLSelSphere_SLST_T3_1);
	TecGUIListDeleteAllItems(MLSelGB_MLST_T3_1);
	TecGUIListDeleteAllItems(SLSelVar_SLST_T3_1);
	TecGUIToggleSet(TGLSphereVis_TOG_T3_1, FALSE);
	/*
	*	First, populate the list of spheres.
	*	Get a total list, then load them alphabetically
	*	with atoms first and cages second.
	*/

	Boolean_t IsOk = TRUE;

	vector<string> SphereCPNameList;

	EntIndex_t NumZones = TecUtilDataSetGetNumZones();

	string TmpStr1, TmpStr2;

	TmpStr1 = CSMAuxData.GBA.ZoneType;
	TmpStr2 = CSMAuxData.GBA.SphereCPName;

	for (EntIndex_t ZoneNum = 1; ZoneNum <= NumZones; ++ZoneNum){
		if (AuxDataZoneItemMatches(ZoneNum, TmpStr1, CSMAuxData.GBA.ZoneTypeSphereZone)){
			SphereCPNameList.push_back(AuxDataZoneGetItem(ZoneNum, TmpStr2));
		}
	}

	IsOk = (SphereCPNameList.size() > 0);
	if (IsOk){
		/*
		*	Sort list of spheres and select first one
		*/
		SortCPNameList(SphereCPNameList);
		for (const string & it : SphereCPNameList){
			TecGUIListAppendItem(SLSelSphere_SLST_T3_1, it.c_str());
		}
		TecGUIListSetSelectedItem(SLSelSphere_SLST_T3_1, GBAResultViewerSelectSphereNum);
		GBAResultViewerSelectSphere();
	}

	TecUtilLockFinish(AddOnID);
}

/**
*/
static void BTNExport_BTN_T3_1_CB(void)
{
	TecUtilLockStart(AddOnID);
	TRACE("Export to CSV Button Pushed\n");
	if ((TecGUIListGetItemCount(SLSelSphere_SLST_T3_1) > 0 && TecGUIToggleGet(TGLExGBs_TOG_T3_1)) || (TecGUIListGetItemCount(SLSelVar_SLST_T3_1) > 0 && TecGUIToggleGet(TGLExInt_TOG_T3_1))){
		char* FolderNameCStr;
		if (TecUtilDialogGetFolderName("Select folder to save files", &FolderNameCStr)){
			string FolderName = FolderNameCStr;
			TecUtilStringDealloc(&FolderNameCStr);
			if (TecGUIListGetItemCount(SLSelVar_SLST_T3_1) > 0 && TecGUIToggleGet(TGLExInt_TOG_T3_1)){
				if (TecGUIListGetItemCount(SLSelVar_SLST_T3_1) > 0){
					vector<int> IntVarNums;
					vector<string> IntVarNames;
					vector<string> IntCheckStrs = { "I: ", "IN: ", "INS: ", " Integration" };
					for (int VarNum = 1; VarNum <= TecUtilDataSetGetNumVars(); ++VarNum){
						char *VarName, *CheckStr;
						if (TecUtilVarGetName(VarNum, &VarName)){
							for (const string & Str : IntCheckStrs){
								CheckStr = std::strstr(VarName, Str.c_str());
								if (CheckStr != NULL){
									string TmpStr = VarName;
									std::replace(TmpStr.begin(), TmpStr.end(), ',', '.');
									IntVarNames.push_back(TmpStr);
									IntVarNums.push_back(VarNum);
									break;
								}
							}
							TecUtilStringDealloc(&VarName);
						}
					}
					for (int ZoneNum = 1; ZoneNum <= TecUtilDataSetGetNumZones(); ++ZoneNum){
						if (AuxDataZoneItemMatches(ZoneNum, CSMAuxData.GBA.ZoneType, CSMAuxData.GBA.ZoneTypeSphereZone)){
							char* ZoneName;
							if (TecUtilZoneGetName(ZoneNum, &ZoneName)){
								vector<int> IJK(3);
								TecUtilZoneGetIJK(ZoneNum, &IJK[0], &IJK[1], &IJK[2]);

								ofstream OutFile(string(FolderName + "/Zone_" + to_string(ZoneNum) + "_" + string(ZoneName) + "_IntegrationResults.csv").c_str(), std::ios::trunc);
								if (OutFile.is_open()){
									OutFile << "Zone," << ZoneName << "\nNumber of gradient bundles (GBs)," << IJK[1];
									string TmpStr;
									if (AuxDataZoneGetItem(ZoneNum, CSMAuxData.GBA.GPsPerGB, TmpStr))
										OutFile << "\nGradient paths (GPs) per GB," << TmpStr;
									if (AuxDataZoneGetItem(ZoneNum, CSMAuxData.GBA.PointsPerGP, TmpStr))
										OutFile << "\nPoints per GP," << TmpStr;
									if (AuxDataZoneGetItem(ZoneNum, CSMAuxData.GBA.IntPrecision, TmpStr))
										OutFile << "\nIntegration precision," << TmpStr;
									if (AuxDataZoneGetItem(ZoneNum, CSMAuxData.GBA.IntWallTime, TmpStr))
										OutFile << "\nIntegration runtime total [seconds]," << TmpStr;
									OutFile << "\n\nIntegration totals\nVariable name,Variable number,Value\n";
									
									vector<FieldDataPointer_c> Ptrs(IntVarNums.size());
									for (int i = 0; i < IntVarNums.size(); ++i)
										Ptrs[i].GetReadPtr(ZoneNum, IntVarNums[i]);

									for (int i = 0; i < Ptrs.size(); ++i){
										OutFile << IntVarNames[i] << "," << i+1 << ",";
										double Total = 0.0;
										for (int j = 0; j < Ptrs[i].GetSize(); ++j){
											Total += Ptrs[i][j];
										}
										OutFile << std::setprecision(16) << std::scientific << Total << '\n';
									}

									OutFile << "\nZone Name,Zone#,GB#";
									for (const string & i : IntVarNames)
										OutFile << "," << i;
									OutFile << '\n';

									for (int i = 0; i < Ptrs[0].GetSize(); ++i){
										char* GBZoneName;
										TecUtilZoneGetName(ZoneNum + i + 1, &GBZoneName);
										OutFile << GBZoneName << "," << ZoneNum + i + 1 << "," << i + 1;
										TecUtilStringDealloc(&GBZoneName);
										for (const auto & j : Ptrs)
											OutFile << std::setprecision(16) << std::scientific << "," << j[i];
										OutFile << "\n";
									}

									OutFile.close();
								}
								else
									TecUtilDialogMessageBox(string("Failed to open output file: " + string(FolderName + "/" + string(ZoneName) + "_IntegrationResults.csv")).c_str(), MessageBoxType_Error);
								TecUtilStringDealloc(&ZoneName);
							}
						}
					}
				}
			}
		}
	}
	TecUtilLockFinish(AddOnID);
}

/**
*/
static void TGLExGBs_TOG_T3_1_CB(const LgIndex_t *I)
{
	TecUtilLockStart(AddOnID);
	TRACE1("Toggle (TGLExGBs_TOG_T3_1) Value Changed,  New value is: %d\n", *I);
	TecUtilLockFinish(AddOnID);
}


/**
*/
static void TGLExInt_TOG_T3_1_CB(const LgIndex_t *I)
{
	TecUtilLockStart(AddOnID);
	TRACE1("Toggle (TGLExInt_TOG_T3_1) Value Changed,  New value is: %d\n", *I);
	TecUtilLockFinish(AddOnID);
}

/**
*/
static void TGLShowMesh_TOG_T3_1_CB(const LgIndex_t *I)
{
	TecUtilLockStart(AddOnID);
	TRACE1("Toggle (TGLShowMesh_TOG_T3_1) Value Changed,  New value is: %d\n", *I);
	TecUtilDrawGraphics(FALSE);

	Set_pa ZoneSet = TecUtilSetAlloc(FALSE);
	for (int ZoneNum = 1; ZoneNum <= TecUtilDataSetGetNumZones(); ++ZoneNum){
		if (AuxDataZoneHasItem(ZoneNum, CSMAuxData.GBA.ZoneType))
			TecUtilSetAddMember(ZoneSet, ZoneNum, FALSE);
	}

	TecUtilZoneSetMesh(SV_SHOW, ZoneSet, 0.0, TecGUIToggleGet(TGLShowMesh_TOG_T3_1));
	TecUtilZoneSetMesh(SV_COLOR, ZoneSet, 0.0, Custom2_C);

	TecUtilSetDealloc(&ZoneSet);

	TecUtilDrawGraphics(TRUE);
	TecUtilLockFinish(AddOnID);
}

#include "guibld.cpp"
