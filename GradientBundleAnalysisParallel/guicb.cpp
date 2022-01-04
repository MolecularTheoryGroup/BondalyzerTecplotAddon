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
#include <map>
#include <set>

#include <armadillo>

#include "ArgList.h"
#include "CSM_DATA_SET_INFO.h"
#include "CSM_DATA_TYPES.h"
#include "CSM_CRIT_POINTS.h"
#include "CSM_GUI.h"
#include "ENGINE.h"
#include "GBAENGINE.h"
#include "INTEGRATE.h"
#include "VIEWRESULTS.h"
#include "CSM_GBAGUI.h"

#include "GUICB.h"

using std::string;
using std::stringstream;
using std::to_string;
using std::vector;
using std::ofstream;
using std::make_pair;

using namespace arma;
using namespace tecplot::toolbox;

bool GBA_REINIT_DIALOG = true;

vector<Text_ID> CPLabelIDs;

LgIndex_t GBATabNumber = 1;
bool GBAIsOpen = false;
string GBAOldNumContours = "15";
int const GBAMaxNumContours = 1000;
int const GBAMinNumContours = 5;

int GBALogLin = 1;
int GBACntSrc = 2;
int GBAResultViewerSelectIntVarNum = 1;
int GBAResultViewerSelectSphereNum = 1;

int PrecisionOffset = 2;


/*
 *	Tooltips for GBA dialog
 */

// void GBAProcessSystemSetToolTips() {
// 	const vector<std::pair<int, string> >
// 		ToolTips = {
// 			make_pair(SCPrecise_SC_T1_1, "Determines the number of integration points that will be used in the Gaussian quadrature integration of gradient bundles. N^2 points will be used for the specified N."),
// 			make_pair(LBLPrecise_LBL_T1_1, "Determines the number of integration points that will be used in the Gaussian quadrature integration of gradient bundles. N^2 points will be used for the specified N."),
// 			make_pair(LBLPreci1_LBL_T1_1, "Determines the number of integration points that will be used in the Gaussian quadrature integration of gradient bundles. N^2 points will be used for the specified N."),
// 
// 			make_pair(TGLsGP_TOG_T1_1, "Save gradient paths used to create gradient bundles as zones"),
// 			make_pair(TGLsGB_TOG_T1_1, "Save gradient bundles as zones"),
// 
// 			make_pair(TFCutoff_TF_T1_1, "Value of electron density at which to truncate gradient bundles"),
// 			make_pair(LBLCutoff_LBL_T1_1, "Value of electron density at which to truncate gradient bundles"),
// 
// 			make_pair(TFSTPts_TF_T1_1, "Number of points used to represent the gradient paths that make up gradient bundle edges"),
// 			make_pair(LBLSTPts_LBL_T1_1, "Number of points used to represent the gradient paths that make up gradient bundle edges"),
// 
// 			make_pair(TFRad_TF_T1_1, "Specify the radius of spheres placed around the specified nuclear critical points according to the option selected below"),
// 			make_pair(LBLRad_LBL_T1_1, "Specify the radius of spheres placed around the specified nuclear critical points according to the option selected below"),
// 			make_pair(RBRadMode_RADIO_T1_1, "If \"Absolute\" is specified then this number is in the same units as the X, Y, and Z variables (i.e. bohr). If Fraction of CP dist then it specifies a factor of the distance to the closest critical point for each of the specified nuclear critical points, and can result in a different radius sphere for each nuclear critical point analyzed."),
// 
// 			make_pair(TFLevel_TFS_T1_1, "Specify the number of times to subdivide the icosahedron used to generate the initial sphere nodes"),
// 			make_pair(LBLNumTri_LBL_T1_1, "Specify the number of times to subdivide the icosahedron used to generate the initial sphere nodes"),
// 			make_pair(LBLLevel_LBL_T1_1, "Specify the number of times to subdivide the icosahedron used to generate the initial sphere nodes"),
// 
// 			make_pair(TFGBPerE_TF_T1_1, "Specify the minimum number of gradient bundles per electron to use in subdivision of gradient bundles (controls amount of charge density per gradient bundle)"),
// 			make_pair(LBLGBPerE_LBL_T1_1, "Specify the minimum number of gradient bundles per electron to use in subdivision of gradient bundles (controls amount of charge density per gradient bundle)"),
// 
// 			make_pair(TFBPGBInit_TF_T1_1, "Specify the number of passes of edge midpoint subdivision for gradient bundles that coincide with bond paths. A gradient bundle is turned into four gradient bundles with each subdivision pass."),
// 			make_pair(LBLBPGBi_LBL_T1_1, "Specify the number of passes of edge midpoint subdivision for gradient bundles that coincide with bond paths. A gradient bundle is turned into four gradient bundles with each subdivision pass."),
// 
// 			make_pair(TFBPGBs_TF_T1_1, "Specify a target number of gradient bundles N that should be coincident with each bond path. Target number of gradient bundles is achieved by splitting gradient bundles around a bond path intersection like putting more slices into a pie such that the resulting angle between gradient bundle edges is not greater than 2*pi / N."),
// 			make_pair(LBLBPGBa_LBL_T1_1, "Specify a target number of gradient bundles N that should be coincident with each bond path. Target number of gradient bundles is achieved by splitting gradient bundles around a bond path intersection like putting more slices into a pie such that the resulting angle between gradient bundle edges is not greater than 2*pi / N."),
// 
// 			make_pair(TFGBMaxSD_TF_T1_1, "Specify maximum number of gradient bundles possible by subdividing a single gradient bundle (e.g. One gradient bundle is subdivided into four each pass, so a value of 4 will limit each gradient bundle to a single pass of subdivision)"),
// 			make_pair(LBLGBMaxSD_LBL_T1_1, "Specify maximum number of gradient bundles possible by subdividing a single gradient bundle (e.g. One gradient bundle is subdivided into four each pass, so a value of 4 will limit each gradient bundle to a single pass of subdivision)")
// 	};
// 
// 	for (auto const & i : ToolTips)
// 		TecGUISetToolTip(i.first, i.second.c_str());
// }

// void GBAResultViewerSetToolTips() {
// 	const vector<std::pair<int, string> >
// 		ToolTips = {
// 					make_pair(SLSelSphere_SLST_T3_1, "Selecting a sphere recalculates contour values for the specified property according to the selected sphere if \"Make contour values from all sphere()\" is selected"),
// 
// 					make_pair(SLSelVar_SLST_T3_1, "Select integrated property to use to generate contour values (\"I: \" properties are the actual integration results, \"IN: \" have been normalized by dividing by the area of the triangular element from which each gradient bundle was made, and \"INS :\" have then been scaled so that the total integrated value is the same as for the \"I: \")"),
// 
// 					make_pair(BTNSmooth_BTN_T3_1, "Perform gaussian smoothing on the selected integrated property (for all spheres or just for the selected sphere). Five passes of smoothing are performed with a sigma coefficient of 0.5.")
// 	};
// 
// 	for (auto const & i : ToolTips)
// 		TecGUISetToolTip(i.first, i.second.c_str());
// }

/*
 *	Integration sidebar functions
 */


// /**
// */
// static void TGLIntVolInt_TOG_T2_1_CB(LgIndex_t const *I)
// {
// 	TecUtilLockStart(AddOnID);
// 	TRACE1("Toggle (TGLIntVolInt_TOG_T2_1) Value Changed,  New value is: %d\n", *I);
// 	TecUtilLockFinish(AddOnID);
// }
// 
// 
// /**
// */
// static void MLIntSelVar_MLST_T2_1_CB(LgIndex_t const *I)
// {
// 	TecUtilLockStart(AddOnID);
// 	TRACE1("Multi selection list (MLIntSelVar_MLST_T2_1) item selected,  First Item is: %d\n", *I);
// 	TecUtilLockFinish(AddOnID);
// }
// 
// 
// /**
// */
// static void BTNExitInt_BTN_T2_1_CB(void)
// {
// 	TecUtilLockStart(AddOnID);
// 	TRACE("Exit Results Viewer Button Pushed\n");
// 	TecGUISidebarActivate(TECGUITECPLOTSIDEBAR);
// 	TecUtilLockFinish(AddOnID);
// }
// 
// 
// /**
// */
// static void MLIntSelSph_MLST_T2_1_CB(LgIndex_t const *I)
// {
// 	TecUtilLockStart(AddOnID);
// 	TRACE1("Multi selection list (MLIntSelSph_MLST_T2_1) item selected,  First Item is: %d\n", *I);
// 	TecUtilLockFinish(AddOnID);
// }
// 
// 
// /**
// */
// static void BTNIntegrate_BTN_T2_1_CB(void)
// {
// 	TecUtilLockStart(AddOnID);
// 	TRACE("Integrate                   Button Pushed\n");
// 	PrepareIntegration(TRUE);
// 	TecUtilLockFinish(AddOnID);
// }
// 
// 
// /**
// */
// static void SCIntPrecise_SC_T2_1_CB(LgIndex_t const *I)
// {
// 	TecUtilLockStart(AddOnID);
// 	TRACE1("Scale (SCIntPrecise_SC_T2_1) Value Changed,  New value is: %d\n", *I);
// 	TecGUILabelSetText(LBLIntPrecis_LBL_T2_1, IntPrecisionLabels[*I].c_str());
// 	IntPrecise = *I;
// 	TecUtilLockFinish(AddOnID);
// }
// 
// 
// /**
// */
// static void SCIntPrecise_SCD_T2_1_CB(LgIndex_t const *I)
// {
// 	TecUtilLockStart(AddOnID);
// 	TRACE1("Scale (SCIntPrecise_SCD_T2_1) Value Changed on drag,  New value is: %d\n", *I);
// 	TecGUILabelSetText(LBLIntPrecis_LBL_T2_1, IntPrecisionLabels[*I].c_str());
// 	IntPrecise = *I;
// 	TecUtilLockFinish(AddOnID);
// }

/*
 *	Results viewer sidebar functions
 */



/**
*/
static void SLSelSphere_SLST_T3_1_CB(LgIndex_t const *I)
{
	TecUtilLockStart(AddOnID);
	TRACE1("Single selection list (SLSelSphere_SLST_T3_1) item selected,  Item is: %d\n", *I);
	int NewGBAResultViewerSelectSphereNum = TecGUIListGetSelectedItem(SLSelSphere_SLST_T3_1);
	if (NewGBAResultViewerSelectSphereNum != GBAResultViewerSelectSphereNum){
		GBAResultViewerSelectSphereNum = NewGBAResultViewerSelectSphereNum;
		CSMGUILock();
		GBAResultViewerSelectSphere();
		if (GBACntSrc == 1)
			GBAResultViewerSelectIntVar();
		CSMGUIUnlock();
		GBAReloadDialog();
	}
	TecUtilLockFinish(AddOnID);
}


/**
*/
static void TGLSphereVis_TOG_T3_1_CB(LgIndex_t const *I)
{
	TecUtilLockStart(AddOnID);
	TRACE1("Toggle (TGLSphere_TOG_T3_1) Value Changed,  New value is: %d\n", *I);
	CSMGUILock();
	GBAResultViewerToggleSphere();
	CSMGUIUnlock();
	GBAReloadDialog();
	TecUtilLockFinish(AddOnID);
}

/**
*/
static void BTNSphereDel_BTN_T3_1_CB(void)
{
	TecUtilLockStart(AddOnID);
	TRACE("Delete Button Pushed\n");
	CSMGUILock();
	GBAResultViewerDeleteSphere();
	CSMGUIUnlock();
	GBAReloadDialog();
	TecUtilLockFinish(AddOnID);
}

/**
 */
static void TGLSmoothAll_TOG_T3_1_CB(const LgIndex_t *I)
{
	TecUtilLockStart(AddOnID);
	TRACE1("Toggle (TGLSmoothAll_TOG_T3_1) Value Changed,  New value is: %d\n", *I);
	TecUtilLockFinish(AddOnID);
}

void SmoothGBAIntResults()
{
	int NumSmoothingPasses = 5;
	double SmoothingCoefficient = 0.5;

	LgIndex_t *SelectedNums;
	LgIndex_t NumSelected;

	TecGUIListGetSelectedItems(SLSelVar_SLST_T3_1, &SelectedNums, &NumSelected);

	if (NumSelected > 0) {
		char *IntVarNameCStr = TecGUIListGetString(SLSelVar_SLST_T3_1, TecGUIListGetSelectedItem(SLSelVar_SLST_T3_1));
		string IntVarNameStr = IntVarNameCStr;
		TecUtilStringDealloc(&IntVarNameCStr);

		if (TecGUIToggleGet(TGLINSV_TOG_T3_1)){
			IntVarNameStr = "INS: " + IntVarNameStr;
		}

		char* SphereNameCStr = TecGUIListGetString(SLSelSphere_SLST_T3_1, TecGUIListGetSelectedItem(SLSelSphere_SLST_T3_1));
		string SphereNameStr = SphereNameCStr;
		TecUtilStringDealloc(&SphereNameCStr);

		// Get the variable number
		int VarNum = TecUtilVarGetNumByName(IntVarNameStr.c_str());

		vector<int> SphereZoneNums;

		Boolean_t AllSpheres = TecGUIToggleGet(TGLSmoothAll_TOG_T3_1);

		for (int ZoneNum = 1; ZoneNum <= TecUtilDataSetGetNumZones(); ++ZoneNum) {
			if (AuxDataZoneItemMatches(ZoneNum, CSMAuxData.GBA.ZoneType, CSMAuxData.GBA.ZoneTypeSphereZone)
				&& (AllSpheres || AuxDataZoneItemMatches(ZoneNum, CSMAuxData.GBA.SphereCPName, SphereNameStr)
					|| AuxDataZoneItemMatches(ZoneNum, CSMAuxData.GBA.SourceNucleusName, SphereNameStr))) {
				SphereZoneNums.push_back(ZoneNum);
				if (!AllSpheres)
					break;
			}
		}

		for (int SphereZoneNum : SphereZoneNums) {
			if (SphereZoneNum > 0) {

				string MacroCmd = "$!SMOOTH ZONE = "
					+ to_string(SphereZoneNum)
					+ " VAR = " + to_string(VarNum)
					+ " NUMSMOOTHPASSES = " + to_string(NumSmoothingPasses)
					+ " SMOOTHWEIGHT = " + to_string(SmoothingCoefficient)
					+ " SMOOTHBNDRYCOND = FIXED";

				TecUtilMacroExecuteCommand(MacroCmd.c_str());
			}
		}
	}
}

/**
 */
static void BTNSmooth_BTN_T3_1_CB(void)
{
	TecUtilLockStart(AddOnID);
	TRACE("Smooth Button Pushed\n");

	CSMGUILock();
	SmoothGBAIntResults();
	CSMGUIUnlock();
	GBAReloadDialog();

	TecUtilLockFinish(AddOnID);
}


/**
*/
static void MLSelGB_MLST_T3_1_CB(LgIndex_t const *I)
{
	TecUtilLockStart(AddOnID);
	TRACE1("Multi selection list (MLSelVol_MLST_T3_1) item selected,  First Item is: %d\n", *I);
// 	GBAResultViewerSelectGB();
	CSMGUILock();
	GBAResultViewerSelectCondensedGBs();
	CSMGUIUnlock();
	GBAReloadDialog();
	TecUtilLockFinish(AddOnID);
}

/**
*/
static void SLSelVar_SLST_T3_1_CB(LgIndex_t const *I)
{
	TecUtilLockStart(AddOnID);
	TRACE1("Single selection list (SLSelVar_SLST_T3_1) item selected,  Item is: %d\n", *I);
	int NewGBAResultViewerSelectIntVarNum = TecGUIListGetSelectedItem(SLSelVar_SLST_T3_1);
	if (NewGBAResultViewerSelectIntVarNum != GBAResultViewerSelectIntVarNum){
		GBAResultViewerSelectIntVarNum = NewGBAResultViewerSelectIntVarNum;
		CSMGUILock();
		GBAResultViewerSelectIntVar();
		CSMGUIUnlock();
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
	CSMGUILock();
	GBAResultViewerActivateAllGB();
	CSMGUIUnlock();
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
	CSMGUILock();
	GBAProcessSystemDeleteCPLabels();
	CSMGUIUnlock();
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
	if (NumSelected > 0) {
		TecGUIDialogDrop(Dialog1Manager);
		CSMGUILock();
		GBAProcessSystemDeleteCPLabels();
		// 		MainFunction();
		NewMainFunction();
		CSMGUIUnlock();
		GBAResultViewerPrepareGUI();
// 		GBATabNumber = 3;
// 		Boolean_t DoIntegration = TecGUIToggleGet(TGLInt_TOG_T1_1);
// 		int NumSelectedVars, *SelectedVarNums;
// 		TecGUIListGetSelectedItems(MLSelVars_MLST_T1_1, &SelectedVarNums, &NumSelectedVars);
// 		DoIntegration = (NumSelectedVars > 0);
// 		if (DoIntegration)
// 			PrepareIntegration(FALSE);
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
static void MLSelCPs_MLST_T1_1_CB(LgIndex_t const *I)
{
	TecUtilLockStart(AddOnID);
	TRACE1("Multi selection list (MLSelCPs_MLST_T1_1) item selected,  First Item is: %d\n", *I);
	CSMGUILock();
	GBAProcessSystemLabelSelectedCPs();
	CSMGUIUnlock();
	GBAReloadDialog();
	TecUtilLockFinish(AddOnID);
}


/**
*/
static void MLSelVars_MLST_T1_1_CB(LgIndex_t const *I)
{
	TecUtilLockStart(AddOnID);
	TRACE1("Multi selection list (MLSelVars_MLST_T1_1) item selected,  First Item is: %d\n", *I);
	TecUtilLockFinish(AddOnID);
}


/**
*/
static LgIndex_t  TFCutoff_TF_T1_1_CB(char const *S)
{
	LgIndex_t IsOk = 1;
	TecUtilLockStart(AddOnID);
	TRACE1("Text field (TFCutoff_TF_T1_1) Value Changed,  New value is: %s\n", S);
	TecUtilLockFinish(AddOnID);
	return (IsOk);
}


/**
*/
static void TGLOpenSys_TOG_T1_1_CB(LgIndex_t const *I)
{
	TecUtilLockStart(AddOnID);
	TRACE1("Toggle (TGLOpenSys_TOG_T1_1) Value Changed,  New value is: %d\n", *I);
	TecUtilLockFinish(AddOnID);
}

// 
// /**
// */
// static LgIndex_t  TFLevel_TFS_T1_1_ValueChanged_CB(char const *S)
// {
// 	LgIndex_t IsOk = 1;
// 	TecUtilLockStart(AddOnID);
// 	TRACE1("Spin Control Text field (TFLevel_TFS_T1_1) Value Changed,  New value is: %s\n", S);
// 	LgIndex_t Value;
// 	if (TecGUITextFieldGetLgIndex(TFLevel_TFS_T1_1, &Value)){
// 		if (Value < GBAMinSphereRefinementLevel)
// 			TecGUITextFieldSetString(TFLevel_TFS_T1_1, to_string(GBAMinSphereRefinementLevel).c_str());
// 		else if (Value > GBAMaxSphereRefinementLevel)
// 			TecGUITextFieldSetString(TFLevel_TFS_T1_1, to_string(GBAMaxSphereRefinementLevel).c_str());
// 	}
// 	else{
// 		TecGUITextFieldSetString(TFLevel_TFS_T1_1, to_string(GBADefaultSphereMeshRefinementLevel).c_str());
// 	}
// 	GBAProcessSystemUpdateNumTriangles();
// 	TecUtilLockFinish(AddOnID);
// 	return (IsOk);
// }
// 
// 
// /**
// */
// static void TFLevel_TFS_T1_1_ButtonUp_CB(void)
// {
// 	TecUtilLockStart(AddOnID);
// 	TRACE0("Spin control (TFLevel_TFS_T1_1) up callback called.\n");
// 	TecGUISpinTextFieldIncLgIndex(TFLevel_TFS_T1_1, 1, GBAMinSphereRefinementLevel, GBAMaxSphereRefinementLevel);
// 	GBAProcessSystemUpdateNumTriangles();
// 	TecUtilLockFinish(AddOnID);
// }
// 
// 
// /**
// */
// static void TFLevel_TFS_T1_1_ButtonDown_CB(void)
// {
// 	TecUtilLockStart(AddOnID);
// 	TRACE0("Spin control (TFLevel_TFS_T1_1) down callback called.\n");
// 	TecGUISpinTextFieldIncLgIndex(TFLevel_TFS_T1_1, -1, GBAMinSphereRefinementLevel, GBAMaxSphereRefinementLevel);
// 	GBAProcessSystemUpdateNumTriangles();
// 	TecUtilLockFinish(AddOnID);
// }

/**
 */
static void SCMinGBs_SC_T1_1_CB(const LgIndex_t *I)
{
	TecUtilLockStart(AddOnID);
	TRACE1("Scale (SCMinGBs_SC_T1_1) Value Changed,  New value is: %d\n", *I);
	GBAProcessSystemUpdateNumTriangles();
	TecUtilLockFinish(AddOnID);
}


/**
 */
static void SCMinGBs_SCD_T1_1_CB(const LgIndex_t *I)
{
	TecUtilLockStart(AddOnID);
	TRACE1("Scale (SCMinGBs_SCD_T1_1) Value Changed on drag,  New value is: %d\n", *I);
	GBAProcessSystemUpdateNumTriangles();
	TecUtilLockFinish(AddOnID);
}

/**
*/
static void SCPrecise_SC_T1_1_CB(LgIndex_t const *I)
{
	TecUtilLockStart(AddOnID);
	TRACE1("Scale (SCPrecise_SC_T1_1) Value Changed,  New value is: %d\n", *I);
	TecGUILabelSetText(LBLPrecise_LBL_T1_1, IntPrecisionLabels[*I-PrecisionOffset].c_str());
	IntPrecise = *I;
	TecUtilLockFinish(AddOnID);
}


/**
*/
static void SCPrecise_SCD_T1_1_CB(LgIndex_t const *I)
{
	TecUtilLockStart(AddOnID);
	TRACE1("Scale (SCPrecise_SCD_T1_1) Value Changed on drag,  New value is: %d\n", *I);
	TecGUILabelSetText(LBLPrecise_LBL_T1_1, IntPrecisionLabels[*I-PrecisionOffset].c_str());
	IntPrecise = *I;
	TecUtilLockFinish(AddOnID);
}


/**
*/
static void TGLVolInt_TOG_T1_1_CB(LgIndex_t const *I)
{
	TecUtilLockStart(AddOnID);
	TRACE1("Toggle (TGLVolInt_TOG_T1_1) Value Changed,  New value is: %d\n", *I);
	TecUtilLockFinish(AddOnID);
}


/**
*/
static LgIndex_t  TFRad_TF_T1_1_CB(char const *S)
{
	LgIndex_t IsOk = 1;
	TecUtilLockStart(AddOnID);
	TRACE1("Text field (TFRad_TF_T1_1) Value Changed,  New value is: %s\n", S); 
	double Value;
	if (!TecGUITextFieldGetDouble(TFRad_TF_T1_1, &Value)){
		TecGUITextFieldSetString(TFRad_TF_T1_1, to_string(GBADefaultRadialSphereApprxRadius).c_str());
	}
	TecUtilLockFinish(AddOnID);
	return (IsOk);
}

/**
 */
static LgIndex_t  TFSdRad_TF_T1_1_CB(const char *S)
{
	LgIndex_t IsOk = 1;
	TecUtilLockStart(AddOnID);
	TRACE1("Text field (TFSdRad_TF_T1_1) Value Changed,  New value is: %s\n", S);
	double Value;
	if (!TecGUITextFieldGetDouble(TFSdRad_TF_T1_1, &Value)) {
		TecGUITextFieldSetString(TFSdRad_TF_T1_1, to_string(GBADefaultSeedRadius).c_str());
	}
	TecUtilLockFinish(AddOnID);
	return (IsOk);
}


/**
*/
// static LgIndex_t  TFIBDist_TF_T1_1_CB(char const *S)
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
// static LgIndex_t  TFIBAng_TF_T1_1_CB(char const *S)
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
static LgIndex_t  TFSTPts_TF_T1_1_CB(char const *S)
{
	LgIndex_t IsOk = 1;
	TecUtilLockStart(AddOnID);
	TRACE1("Text field (TFSTPts_TF_T1_1) Value Changed,  New value is: %s\n", S);
	int Value;
	Boolean_t ValRead = TecGUITextFieldGetLgIndex(TFSTPts_TF_T1_1, &Value);
	if (!ValRead){
		TecGUITextFieldSetString(TFSTPts_TF_T1_1, to_string(DefaultNumGPPts).c_str());
	}
	else if (Value % 2 != 0){
		TecGUITextFieldSetString(TFSTPts_TF_T1_1, to_string(--Value).c_str());
	}
	TecUtilLockFinish(AddOnID);
	return (IsOk);
}

/**
*/
static void RBRadMode_RADIO_T1_1_CB(LgIndex_t const *I)
{
	TecUtilLockStart(AddOnID);
	TRACE1("RadioBox (RBRadMode_RADIO_T1_1) Value Changed,  New value is: %d\n", *I);
	TecUtilLockFinish(AddOnID);
}

Boolean_t GBAProcessSystemPrepareGUI(){
// 	GBAProcessSystemSetToolTips();

	/*
	*	Prepare CP list
	*/
	int NumSelected = -1, 
		OldAtomCount = TecGUIListGetItemCount(MLSelCPs_MLST_T1_1), 
		OldVarCount = TecGUIListGetItemCount(MLSelVars_MLST_T1_1);
	int *SelectedNums;
	if (OldAtomCount > 0) {
		if (GBA_REINIT_DIALOG) {
			TecGUIListGetSelectedItems(MLSelCPs_MLST_T1_1, &SelectedNums, &NumSelected);
		}
		TecGUIListDeleteAllItems(MLSelCPs_MLST_T1_1);
	}

	/*
	 *	Get nuclear CPs according to nuclear name.
	 */
	vector<int> XYZVarNums;
	for (auto const & i : { "X","Y","Z" }){
		XYZVarNums.push_back(VarNumByName(i));
		if (XYZVarNums.back() <= 0){
			XYZVarNums = { 1,2,3 };
			break;
		}
	}
	int CPTypeVarNum = VarNumByName(CSMVarName.CritPointType);
	if (CPTypeVarNum > 0) {
		CritPoints_c NuclearCPs;
		for (int ZoneNum = 1; ZoneNum <= TecUtilDataSetGetNumZones(); ++ZoneNum) {
			if (AuxDataZoneItemMatches(ZoneNum, CSMAuxData.CC.ZoneSubType, CSMAuxData.CC.CPSubTypes[0])) {
				NuclearCPs += CritPoints_c(ZoneNum, { 1,2,3 }, CPTypeVarNum);
			}
		}

		NuclearNameToCPNum.clear();
		std::set<int> CPIndSet;
		double DistCutoff = 0.2;

		for (int ZoneNum = 1; ZoneNum <= TecUtilDataSetGetNumZones(); ++ZoneNum) {
			char* ZoneNameCStr;
			if (AuxDataZoneItemMatches(ZoneNum, CSMAuxData.DL.ZoneType, CSMAuxData.DL.ZoneTypeNuclearPositions) && TecUtilZoneGetName(ZoneNum, &ZoneNameCStr)) {
				string ZoneName = ZoneNameCStr;
				TecUtilStringDealloc(&ZoneNameCStr);

				int ElemNum = SearchVectorForString(ElementSymbolList, ZoneName, false);
				if (ElemNum < 0)
					ElemNum = SearchVectorForString(ElementNameList, ZoneName);
				if (ElemNum >= 0) {
					int NumPts;
					TecUtilZoneGetIJK(ZoneNum, &NumPts, nullptr, nullptr);
					for (int Counter = 1; Counter <= NumPts; ++Counter) {
						vec3 Pt;
						for (int i = 0; i < 3; ++i){
							Pt[i] = TecUtilDataValueGetByZoneVar(ZoneNum, XYZVarNums[i], Counter);
						}
						int CPInd;
						double CPDist;
						vec3 ClosestCPPt = NuclearCPs.ClosestPoint(Pt, CPInd, CPDist);
						if (CPDist < DistCutoff){
							string CPName = ZoneName + " " + to_string(Counter);
							NuclearNameToCPNum[CPName] = CPInd;
							CPIndSet.insert(CPInd);
							TecGUIListAppendItem(MLSelCPs_MLST_T1_1, CPName.c_str());
						}
					}
				}
			}
		}

		// Now add nuclear CPs that weren't sufficiently close to one of the nuclear positions
		for (int i = 0; i < NuclearCPs.NumAtoms(); ++i){
			if (!CPIndSet.count(i)){
				string CPName = CPNameList[0] + " " + to_string(i + 1);
				NuclearNameToCPNum[CPName] = i;
				TecGUIListAppendItem(MLSelCPs_MLST_T1_1, CPName.c_str());
			}
		}
	}
// 	else {
// 		/*
// 		 *	Add CPs from Bondalyzer
// 		 */
// 
// 		EntIndex_t CPZoneNum = ZoneNumByName(string("Critical Points"), false);
// 		EntIndex_t CPVarNum = VarNumByName(string("CritPointType"));
// 		// 	if (CPZoneNum < 0 || CPVarNum < 0){
// 		// 		TecUtilDialogErrMsg("Couldn't find critical points zone or CP type variable.");
// 		// 		return FALSE;
// 		// 	}
// 		if (CPZoneNum > 0 && CPVarNum > 0) {
// 			int Ranks[4] = { -3, -1, 1, 3 };
// 			LgIndex_t IJK[3]; LgIndex_t NumCPs[4] = { 0, 0, 0, 0 };
// 			TecUtilZoneGetIJK(CPZoneNum, &IJK[0], &IJK[1], &IJK[2]);
// 			for (int i = 1; i <= IJK[0]; ++i) {
// 				int CPType = (int)TecUtilDataValueGetByZoneVar(CPZoneNum, CPVarNum, i);
// 				for (int j = 0; j < 4; ++j) {
// 					if (CPType == Ranks[j]) {
// 						NumCPs[j]++;
// 						if (CPType == -3 || CPType == 3) {
// 							stringstream ss;
// 							int CPOffset = 0;
// 							for (int k = 0; k < j; ++k)
// 								CPOffset += NumCPs[k];
// 							ss << CPNameList[j] << " " << i - CPOffset;
// 							TecGUIListAppendItem(MLSelCPs_MLST_T1_1, ss.str().c_str());
// 						}
// 						break;
// 					}
// 				}
// 			}
// 		}
// 	}

	if (NumSelected > 0 && OldAtomCount == TecGUIListGetItemCount(MLSelCPs_MLST_T1_1)) {
		TecGUIListSetSelectedItems(MLSelCPs_MLST_T1_1, SelectedNums, NumSelected);
	}



	/*
	*	Prepare Var list
	*/

	TecGUIScaleShowNumericDisplay(SCPrecise_SC_T1_1, TRUE);
	TecGUIScaleSetLimits(SCPrecise_SC_T1_1, PrecisionOffset, PrecisionOffset + IntPrecisionLabels.size() - 1, 0);
	TecGUIScaleSetValue(SCPrecise_SC_T1_1, GBADefaultIntegrationPrecision);
	TecGUILabelSetText(LBLPrecise_LBL_T1_1, IntPrecisionLabels[GBADefaultIntegrationPrecision - PrecisionOffset].c_str());
// 	TecGUIScaleSetValue(SCPrecise_SC_T1_1, GetRecommendedIntPrecision());


	NumSelected = -1;
	if (GBA_REINIT_DIALOG && TecGUIListGetItemCount(MLSelVars_MLST_T1_1) > 0)
		TecGUIListGetSelectedItems(MLSelVars_MLST_T1_1, &SelectedNums, &NumSelected);

	if (TecGUIListGetItemCount(MLSelVars_MLST_T1_1) > 0) {
		TecGUIListDeleteAllItems(MLSelVars_MLST_T1_1);
	}
	TecGUIToggleSet(TGLsGP_TOG_T1_1, FALSE);
	TecGUIToggleSet(TGLsGB_TOG_T1_1, FALSE);
// 	TecGUIToggleSet(TGLVolInt_TOG_T1_1, DefaultVolIntegrate);
	ListPopulateWithVarNames(MLSelVars_MLST_T1_1, 4);
	if (NumSelected > 0 && OldVarCount == TecGUIListGetItemCount(MLSelVars_MLST_T1_1))
		TecGUIListSetSelectedItems(MLSelVars_MLST_T1_1, SelectedNums, NumSelected);
	else{
		int SelectNums[] = { 1 };
		TecGUIListSetSelectedItems(MLSelVars_MLST_T1_1, SelectNums, 1);
	}
// 	TecGUITextFieldSetString(TFGBPerE_TF_T1_1, to_string(GBADefaultGBPerE).c_str());
// 	TecGUITextFieldSetString(TFBPGBs_TF_T1_1, to_string(GBADefaultBPAngularGBs).c_str());
// 	TecGUITextFieldSetString(TFGBMaxSD_TF_T1_1, to_string(GBADefaultMaxGBSubdivisionLevel).c_str());
// 	TecGUITextFieldSetString(TFBPGBInit_TF_T1_1, to_string(GBADefaultNumberOfPreBondPathElemSubdivision).c_str());
// 	
	TecGUIScaleSetLimits(SCBPGBInit_SC_T1_1, 0, GBAMaxMaxSubdivisionLevel, 0);
	TecGUIScaleSetValue(SCBPGBInit_SC_T1_1, GBAMaxSubdivisionLevel);
	TecGUILabelSetText(LBLBPGBInit_LBL_T1_1, to_string(GBAMaxSubdivisionLevel).c_str());

	TecGUIScaleSetLimits(SCSDtight_SC_T1_1, GBAMinSubdivisionTightness, GBAMaxSubdivisionTightness, 0);
	TecGUIScaleSetValue(SCSDtight_SC_T1_1, GBADefaultSubdivisionTightness);
	TecGUILabelSetText(LBLSDtight_LBL_T1_1, to_string(GBADefaultSubdivisionTightness).c_str());

	TecGUIScaleSetLimits(SCBPGBs_SC_T1_1, 1, GBADefaultMaxBPAngularGBs, 0);
	TecGUIScaleSetValue(SCBPGBs_SC_T1_1, GBADefaultBPAngularGBs);
	TecGUILabelSetText(LBLBPGBs_LBL_T1_1, to_string(GBADefaultBPAngularGBs).c_str());

// 	double VolumeGridSpacing = norm(GetDelXYZ_Ordered3DZone(XYZVarNums, ZoneNumByName("Full Volume", false, true)));
	TecGUIScaleSetLimits(SCEGPDist_SC_T1_1, GBADefaultMinEdgeGPSpacing * 10, GBADefaultMaxEdgeGPSpacing * 10, 1);
	TecGUIScaleSetValue(SCEGPDist_SC_T1_1, GBADefaultEdgeGPSpacing * 10);
	TecGUILabelSetText(LBLEGPDistTy_LBL_T1_1, DoubleToString(GBADefaultEdgeGPSpacing, 1).c_str());

	TecGUIToggleSet(TGLNoSphInt_TOG_T1_1, TRUE);
	

	TecGUIToggleSet(TGLSmoothAll_TOG_T3_1, TRUE);

	/*
	*	System boundary options
	*/
	TecGUITextFieldSetString(TFCutoff_TF_T1_1, to_string(1e-3).c_str());

	/*
	*	Mesh parameters
	*/
	TecGUITextFieldSetString(TFRad_TF_T1_1, to_string(GBADefaultRadialSphereApprxRadius).c_str());
	TecGUITextFieldSetString(TFSdRad_TF_T1_1, to_string(GBADefaultSeedRadius).c_str());
	TecGUITextFieldSetString(TFHSdR_TF_T1_1, to_string(GBADefaultHSeedRadius).c_str());
	TecGUIRadioBoxSetToggle(RBRadMode_RADIO_T1_1, 2);
	TecGUIScaleSetLimits(SCMinGBs_SC_T1_1, GBAMinSphereRefinementLevel, GBAMaxSphereRefinementLevel, 0);
	TecGUIScaleSetValue(SCMinGBs_SC_T1_1, GBADefaultSphereMeshRefinementLevel);
	GBAProcessSystemUpdateNumTriangles();
	TecGUITextFieldSetString(TFSTPts_TF_T1_1, to_string(DefaultNumGPPts).c_str());

	/*
	*	IB detection parameters
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
	TecGUITextFieldSetString(TFGBSub_TF_T3_1, "3");
	TecGUIToggleSet(TGLShowMesh_TOG_T3_1, FALSE);
	if (TecGUIListGetItemCount(SLSelVar_SLST_T3_1) > 0) {
		string itemStr = TecGUIListGetString(SLSelVar_SLST_T3_1, TecGUIListGetSelectedItem(SLSelVar_SLST_T3_1));
		bool IsFound = false;
		for (auto & s : { "I: ", "IN: ", "INS: " }) {
			if (itemStr.find(s) != string::npos) {
				TecGUIToggleSet(TGLINSV_TOG_T3_1, FALSE);
				IsFound = true;
				break;
			}
		}
		if (!IsFound) {
			TecGUIToggleSet(TGLINSV_TOG_T3_1, TRUE);
		}
	}
	else {
		TecGUIToggleSet(TGLINSV_TOG_T3_1, TRUE);
	}

	return TRUE;
}

string humancount(long int N, long int d = 1){
	long int K = 1000,
		M = K * 1000,
		G = M * 1000;

	stringstream stream;

	if (N < K) {
		stream << N;
	}
	else if (N < M) {
		stream << std::fixed << std::setprecision(d) << double(N) / double(K) << "k";
	}
	else if (N < G) {
		stream << std::fixed << std::setprecision(d) << double(N) / double(M) << "M";
	}
	else {
		stream << std::fixed << std::setprecision(d) << double(N) / double(G) << "G";
	}

	return stream.str();
}

void GBAProcessSystemUpdateNumTriangles() {
	stringstream ss;
	LgIndex_t Level = TecGUIScaleGetValue(SCMinGBs_SC_T1_1);
	long int NumTriangles = 20 * (int)pow(4, Level);
	LgIndex_t MaxLevel = Level + TecGUIScaleGetValue(SCBPGBInit_SC_T1_1);
	if (MaxLevel == Level){
		ss << "appx. " << humancount(NumTriangles) << " GBs";
	}
	else {
		long int MaxNumTriangles = double(NumTriangles) + (20 * (int)pow(4, MaxLevel)) / pow(2,TecGUIScaleGetValue(SCSDtight_SC_T1_1));
		ss << humancount(NumTriangles) << "-" << humancount(MaxNumTriangles) << " GBs";
	}
	TecGUILabelSetText(LBLNumTri_LBL_T1_1, ss.str().c_str());
}


void GBAProcessSystemLabelSelectedCPs(){
	LgIndex_t * SelectedCPs;
	LgIndex_t NumSelected;
	TecGUIListGetSelectedItems(MLSelCPs_MLST_T1_1, &SelectedCPs, &NumSelected);
	if (NumSelected > 0){
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
			vec3 CPXYZ = GetCoordsFromListItem(SelectedCPs[i], MLSelCPs_MLST_T1_1, &CPName, nullptr, nullptr, nullptr, nullptr, &NumCPs);
			vec3 XYZGrid;
			TecUtilConvert3DPositionToGrid(CPXYZ[0], CPXYZ[1], CPXYZ[2], &XYZGrid[0], &XYZGrid[1], &XYZGrid[2]);


			Text_ID LabelID = TecUtilTextCreate(CoordSys_Grid, XYZGrid[0] + TextOffset, XYZGrid[1] + TextOffset, Units_Frame, 2, CPName.c_str());
			TecUtilTextBoxSetType(LabelID, TextBox_Filled);
			TecUtilTextBoxSetFillColor(LabelID, White_C);

			Boolean_t IsValid = TecUtilTextIsValid(LabelID);

			CPLabelIDs.push_back(LabelID);
		}
	}
	TecUtilArrayDealloc((void **)&SelectedCPs);
}

void GBAProcessSystemDeleteCPLabels(){
	MouseButtonMode_e MouseMode = TecUtilMouseGetCurrentMode();

	Text_ID Text;

	Text = TecUtilTextGetBase();
	while (Text != TECUTILBADID)
	{
		TecUtilTextDelete(Text);
		Text = TecUtilTextGetBase();
	}

// 	if (CPLabelIDs.size() > 0){
// 		TecUtilPickDeselectAll();
// 		for (int i = 0; i < CPLabelIDs.size(); ++i){
// 			TecUtilPickText(CPLabelIDs[i]);
// 		}
// 		TecUtilPickClear();
// 		CPLabelIDs.clear();
// 	}
	if (TecUtilMouseIsValidMode(MouseMode))
		TecUtilMouseSetMode(MouseMode);
}

void PrepareIntegration(Boolean_t IntegratingFromIntTab){
	LgIndex_t AtomListID, VarListID, TGLID, ScaleID;
	Boolean_t AllGBs;
// 	if (IntegratingFromIntTab){
// 		AtomListID = MLIntSelSph_MLST_T2_1;
// 		VarListID = MLIntSelVar_MLST_T2_1;
// 		TGLID = TGLIntVolInt_TOG_T2_1;
// 		ScaleID = SCIntPrecise_SC_T2_1;
// 		AllGBs = TecGUIToggleGet(TGLIntAct2_TOG_T2_1);
// 	}
// 	else{
		AtomListID = MLSelCPs_MLST_T1_1;
		VarListID = MLSelVars_MLST_T1_1;
// 		TGLID = TGLVolInt_TOG_T1_1;
		ScaleID = SCPrecise_SC_T1_1;
		AllGBs = TRUE;
// 	}

	vector<string> AtomNameList = ListGetSelectedStrings(AtomListID);
	vector<string> IntVarNameList = ListGetSelectedStrings(VarListID);
	vector<int> IntVarNumList = ListGetSelectedItemNums(VarListID);

	PerformIntegration(AtomNameList, IntVarNameList, IntVarNumList,
		TecGUIToggleGet(TGLID), TecGUIScaleGetValue(ScaleID), AllGBs);
	/*for (int i = 1; i < 5; ++i){
		PerformIntegration(AtomNameList, IntVarNameList, IntVarNumList,
			TecGUIToggleGet(TGLID), i);
	}*/


	TecGUITabSetCurrentPage(TAB1_TB_D1, 3);
}

/**
 */
static void TGLNoSphInt_TOG_T1_1_CB(const LgIndex_t *I)
{
	TecUtilLockStart(AddOnID);
	TRACE1("Toggle (TGLNoSphInt_TOG_T1_1) Value Changed,  New value is: %d\n", *I);
	TecUtilLockFinish(AddOnID);
}


/**
 */
static void SCBPGBInit_SC_T1_1_CB(const LgIndex_t *I)
{
	TecUtilLockStart(AddOnID);
	TRACE1("Scale (SCBPGBInit_SC_T1_1) Value Changed,  New value is: %d\n", *I);
	GBAProcessSystemUpdateNumTriangles();
	TecGUILabelSetText(LBLBPGBInit_LBL_T1_1, to_string(*I).c_str());
	TecUtilLockFinish(AddOnID);
}


/**
 */
static void SCBPGBInit_SCD_T1_1_CB(const LgIndex_t *I)
{
	TecUtilLockStart(AddOnID);
	TRACE1("Scale (SCBPGBInit_SCD_T1_1) Value Changed on drag,  New value is: %d\n", *I);
	GBAProcessSystemUpdateNumTriangles();
	TecGUILabelSetText(LBLBPGBInit_LBL_T1_1, to_string(*I).c_str());
	TecUtilLockFinish(AddOnID);
}


/**
 */
static void SCBPGBs_SC_T1_1_CB(const LgIndex_t *I)
{
	TecUtilLockStart(AddOnID);
	TRACE1("Scale (SCBPGBs_SC_T1_1) Value Changed,  New value is: %d\n", *I);
	TecGUILabelSetText(LBLBPGBs_LBL_T1_1, to_string(*I).c_str());
	TecUtilLockFinish(AddOnID);
}


/**
 */
static void SCBPGBs_SCD_T1_1_CB(const LgIndex_t *I)
{
	TecUtilLockStart(AddOnID);
	TRACE1("Scale (SCBPGBs_SCD_T1_1) Value Changed on drag,  New value is: %d\n", *I);
	TecGUILabelSetText(LBLBPGBs_LBL_T1_1, to_string(*I).c_str());
	TecUtilLockFinish(AddOnID);
}


/**
 */
static void SCEGPDist_SC_T1_1_CB(const LgIndex_t *I)
{
	TecUtilLockStart(AddOnID);
	TRACE1("Scale (SCEGPDist_SC_T1_1) Value Changed,  New value is: %d\n", *I);
	TecGUILabelSetText(LBLEGPDistTy_LBL_T1_1, DoubleToString(double(*I)/10., 1).c_str());
	TecUtilLockFinish(AddOnID);
}


/**
 */
static void SCEGPDist_SCD_T1_1_CB(const LgIndex_t *I)
{
	TecUtilLockStart(AddOnID);
	TRACE1("Scale (SCEGPDist_SCD_T1_1) Value Changed on drag,  New value is: %d\n", *I);
	TecGUILabelSetText(LBLEGPDistTy_LBL_T1_1, DoubleToString(double(*I) / 10., 1).c_str());
	TecUtilLockFinish(AddOnID);
}

/**
 */
// static LgIndex_t  TFGBPerE_TF_T1_1_CB(const char *S)
// {
// 	LgIndex_t IsOk = 1;
// 	TecUtilLockStart(AddOnID);
// 	TRACE1("Text field (TFGBPerE_TF_T1_1) Value Changed,  New value is: %s\n", S);
// 	if (StringIsInt(S) && stoi(string(S)) > 0){
// 		TecGUITextFieldSetString(TFGBPerE_TF_T1_1, S);
// 	}
// 	TecUtilLockFinish(AddOnID);
// 	return (IsOk);
// }

/**
 */
// static LgIndex_t  TFBPGBInit_TF_T1_1_CB(const char *S)
// {
// 	LgIndex_t IsOk = 1;
// 	TecUtilLockStart(AddOnID);
// 	TRACE1("Text field (TFBPGBInit_TF_T1_1) Value Changed,  New value is: %s\n", S);
// 	if (StringIsInt(S)) {
// 		int s = CLAMP(stoi(string(S)), 0, GBADefaultMaxNumberOfPreBondPathElemSubdivision);
// 		TecGUITextFieldSetString(TFBPGBInit_TF_T1_1, to_string(s).c_str());
// 	}
// 	TecUtilLockFinish(AddOnID);
// 	return (IsOk);
// }

/**
 */
// static LgIndex_t  TFBPGBs_TF_T1_1_CB(const char *S)
// {
// 	LgIndex_t IsOk = 1;
// 	TecUtilLockStart(AddOnID);
// 	TRACE1("Text field (TFBPGBs_TF_T1_1) Value Changed,  New value is: %s\n", S);
// 	if (StringIsInt(S) && stoi(string(S)) > 0) {
// 		TecGUITextFieldSetString(TFBPGBs_TF_T1_1, S);
// 	}
// 	TecUtilLockFinish(AddOnID);
// 	return (IsOk);
// }

/**
 */
// static LgIndex_t  TFGBMaxSD_TF_T1_1_CB(const char *S)
// {
// 	LgIndex_t IsOk = 1;
// 	TecUtilLockStart(AddOnID);
// 	TRACE1("Text field (TFGBMaxSD_TF_T1_1) Value Changed,  New value is: %s\n", S);
// 	if (StringIsInt(S) && stoi(string(S)) >= -1) {
// 		TecGUITextFieldSetString(TFGBMaxSD_TF_T1_1, S);
// 	}
// 	TecUtilLockFinish(AddOnID);
// 	return (IsOk);
// }

/**
*/
static void TAB1_TBA_D1_CB(LgIndex_t const *I)
{
	TecUtilLockStart(AddOnID);
	TRACE0("Activate callback for tab (TAB1_TBA_D1) called\n");
	GBATabNumber = *I;
	switch (*I){
		case 1:
			GBAProcessSystemPrepareGUI();
			break;
		case 2:
// 			GBAIntegrationPrepareGUI();
// 			break;
// 		case 3:
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
static void TAB1_TBD_D1_CB(LgIndex_t const *I)
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
static void RBCntSrc_RADIO_T3_1_CB(LgIndex_t const *I)
{
	TecUtilLockStart(AddOnID);
	TRACE1("RadioBox (RBCntSrc_RADIO_T3_1) Value Changed,  New value is: %d\n", *I);
	GBACntSrc = TecGUIRadioBoxGetToggle(RBCntSrc_RADIO_T3_1);
	CSMGUILock();
	GBAResultViewerSelectIntVar();
	CSMGUIUnlock();
	TecUtilLockFinish(AddOnID);
}

/**
*/
static void RBLogLin_RADIO_T3_1_CB(LgIndex_t const *I)
{
	TecUtilLockStart(AddOnID);
	TRACE1("RadioBox (RBLogLin_RADIO_T3_1) Value Changed,  New value is: %d\n", *I);
	GBALogLin = TecGUIRadioBoxGetToggle(RBLogLin_RADIO_T3_1);
	CSMGUILock();
	GBAResultViewerSelectIntVar();
	CSMGUIUnlock();
	TecUtilLockFinish(AddOnID);
}

/**
*/
static LgIndex_t  TFNumContours_TF_T3_1_CB(char const *S)
{
	LgIndex_t IsOk = 1;
	TecUtilLockStart(AddOnID);
	string TFStr = S;
	bool IsNumber = true;
	for (char i : TFStr){
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

// 	TecUtilLockStart(AddOnID);

// 	GBAResultViewerSetToolTips();

	/*
	*	Clear lists and stuff
	*/

	TecGUIListDeleteAllItems(SLSelSphere_SLST_T3_1);
	TecGUIListDeleteAllItems(MLSelGB_MLST_T3_1);
	TecGUIListDeleteAllItems(SLSelVar_SLST_T3_1);
	TecGUIToggleSet(TGLSphereVis_TOG_T3_1, FALSE);
	TecGUITextFieldSetString(TFGrpNum_TF_T3_1, "10");
	/*
	*	First, populate the list of spheres.
	*	Get a total list, then load them alphabetically
	*	with atoms first and cages second.
	*/

	Boolean_t IsOk = TRUE;

	std::set<string> SphereCPNameList;

	EntIndex_t NumZones = TecUtilDataSetGetNumZones();

	for (EntIndex_t ZoneNum = 1; ZoneNum <= NumZones; ++ZoneNum){
		if (AuxDataZoneItemMatches(ZoneNum, CSMAuxData.GBA.ZoneType, CSMAuxData.GBA.ZoneTypeSphereZone)){
			string SphereName = AuxDataZoneGetItem(ZoneNum, CSMAuxData.GBA.SourceNucleusName);
			if (SphereName == "")
				SphereName = AuxDataZoneGetItem(ZoneNum, CSMAuxData.GBA.SphereCPName);
			SphereCPNameList.insert(SphereName);
		}
	}

	if (TecGUIListGetItemCount(SLSelVar_SLST_T3_1) > 0){
		string itemStr = TecGUIListGetString(SLSelVar_SLST_T3_1, TecGUIListGetSelectedItem(SLSelVar_SLST_T3_1));
		bool IsFound = false;
		for (auto & s : { "I: ", "IN: ", "INS: " }){
			if (itemStr.find(s) != string::npos){
				TecGUIToggleSet(TGLINSV_TOG_T3_1, FALSE);
				IsFound = true;
				break;
			}
		}
		if (!IsFound){
			TecGUIToggleSet(TGLINSV_TOG_T3_1, TRUE);
		}
	}
	else {
		TecGUIToggleSet(TGLINSV_TOG_T3_1, TRUE);
	}
	ResultsVarListReload();

	TecGUIScaleSetLimits(SCRad_SC_T3_1, 1, 50, 1);
	TecGUIScaleSetValue(SCRad_SC_T3_1, 10);
	TecGUIToggleSet(TGLRadAbs_TOG_T3_1, FALSE);
	TecGUIToggleSet(TGLRadAll_TOG_T3_1, TRUE);
	TecGUILabelSetText(LBLRadLab_LBL_T3_1, "1.0");

	IsOk = (SphereCPNameList.size() > 0);
	if (IsOk){
		/*
		*	Sort list of spheres and select first one
		*/
// 		SortCPNameList(SphereCPNameList);
		for (string const & it : SphereCPNameList){
			TecGUIListAppendItem(SLSelSphere_SLST_T3_1, it.c_str());
		}
		TecGUIListSetSelectedItem(SLSelSphere_SLST_T3_1, GBAResultViewerSelectSphereNum);
		GBAResultViewerSelectSphere();

		GBAResultViewerPopulateGBs();
	}

// 	TecUtilLockFinish(AddOnID);
}

/**
*/
static void BTNExport_BTN_T3_1_CB(void)
{
	TecUtilLockStart(AddOnID);
	TRACE("Export to CSV Button Pushed\n");
// 	ExportGBADataGetUserInfo();
// 	ExportGBAData();
// 	
	TecUtilDialogMessageBox("Please use the export tool located in the MTG_Utilities menu.", MessageBoxType_Information);
	TecUtilLockFinish(AddOnID);
}

/**
*/
static void TGLExGBs_TOG_T3_1_CB(LgIndex_t const *I)
{
	TecUtilLockStart(AddOnID);
	TRACE1("Toggle (TGLExGBs_TOG_T3_1) Value Changed,  New value is: %d\n", *I);
	TecUtilLockFinish(AddOnID);
}


/**
*/
static void TGLExInt_TOG_T3_1_CB(LgIndex_t const *I)
{
	TecUtilLockStart(AddOnID);
	TRACE1("Toggle (TGLExInt_TOG_T3_1) Value Changed,  New value is: %d\n", *I);
	TecUtilLockFinish(AddOnID);
}

/**
*/
static void TGLShowMesh_TOG_T3_1_CB(LgIndex_t const *I)
{
	TecUtilLockStart(AddOnID);
	TRACE1("Toggle (TGLShowMesh_TOG_T3_1) Value Changed,  New value is: %d\n", *I);
	CSMGUILock();

	Set_pa ZoneSet = TecUtilSetAlloc(FALSE);
	for (int ZoneNum = 1; ZoneNum <= TecUtilDataSetGetNumZones(); ++ZoneNum){
		if (AuxDataZoneHasItem(ZoneNum, CSMAuxData.GBA.ZoneType))
			TecUtilSetAddMember(ZoneSet, ZoneNum, FALSE);
	}

	TecUtilZoneSetMesh(SV_SHOW, ZoneSet, 0.0, TecGUIToggleGet(TGLShowMesh_TOG_T3_1));
	TecUtilZoneSetMesh(SV_COLOR, ZoneSet, 0.0, Red_C);

	TecUtilSetDealloc(&ZoneSet);
	
	CSMGUIUnlock();
	TecUtilLockFinish(AddOnID);
}

/**
*/
static void BTNSelGB_BTN_T3_1_CB(void)
{
	TecUtilLockStart(AddOnID);
	TRACE("Select GB Button Pushed\n");
	/*
	 *	Check to see if probe is enabled, and switch if it is
	 */
	if (TecUtilMouseGetCurrentMode() == MouseButtonMode_Probe){
		TecUtilMouseSetMode(Mouse_RotateRollerBall);
	}
	ArgList args;
	args.appendFunction(SV_CALLBACKFUNCTION, SelectGBsInRegionProbeCB);
	args.appendString(SV_STATUSLINETEXT, "Select an element interior to the boundary of active GBs");
	TecUtilProbeInstallCallbackX(args.getRef());
	TecUtilLockFinish(AddOnID);
}

/**
*/
static LgIndex_t  TFGrpNum_TF_T3_1_CB(char const *S)
{
	LgIndex_t IsOk = 1;
	TecUtilLockStart(AddOnID);
	TRACE1("Text field (TFGrpNum_TF_T3_1) Value Changed,  New value is: %s\n", S);
	if (!StringIsInt(S)){
		TecGUITextFieldSetString(TFGrpNum_TF_T3_1, TecGUITextFieldGetString(TFGrpNum_TF_T3_1));
	}
	TecUtilLockFinish(AddOnID);
	return (IsOk);
}

// /**
// */
// static void TGLIntAct2_TOG_T2_1_CB(LgIndex_t const *I)
// {
// 	TecUtilLockStart(AddOnID);
// 	TRACE1("Toggle (TGLIntAct2_TOG_T2_1) Value Changed,  New value is: %d\n", *I);
// 	TecUtilLockFinish(AddOnID);
// }
// 
/**
*/
static void TGLsGP_TOG_T1_1_CB(LgIndex_t const *I)
{
	TecUtilLockStart(AddOnID);
	TRACE1("Toggle (TGLsGP_TOG_T1_1) Value Changed,  New value is: %d\n", *I);
	TecUtilLockFinish(AddOnID);
}

/**
*/
static void TGLsGB_TOG_T1_1_CB(LgIndex_t const *I)
{
	TecUtilLockStart(AddOnID);
	TRACE1("Toggle (TGLsGB_TOG_T1_1) Value Changed,  New value is: %d\n", *I);
	TecUtilLockFinish(AddOnID);
}

/**
 */
static void BTNFndBas_BTN_T3_1_CB(void)
{
	TecUtilLockStart(AddOnID);
	TRACE("Find basins Button Pushed\n");
	CSMGUILock();
	FindSphereBasins();
	CSMGUIUnlock();
	GBAResultViewerPrepareGUI();
	TecUtilLockFinish(AddOnID);
}

void ResizeSpheresCallback(bool DoResize = true) {
	double NewVal = (double)(TecGUIScaleGetValue(SCRad_SC_T3_1)) / 10.0;
	TecGUILabelSetText(LBLRadLab_LBL_T3_1, DoubleToString(NewVal, 2).c_str());
	if (DoResize) {
		ResizeSpheres(NewVal, TecGUIToggleGet(TGLRadAll_TOG_T3_1), TecGUIToggleGet(TGLRadAbs_TOG_T3_1));
		GBAReloadDialog();
	}
}

/**
 */
static void SCRad_SC_T3_1_CB(const LgIndex_t *I)
{
	TecUtilLockStart(AddOnID);
	TRACE1("Scale (SCRad_SC_T3_1) Value Changed,  New value is: %d\n", *I);
	ResizeSpheresCallback();
	TecUtilLockFinish(AddOnID);
}


/**
 */
static void SCRad_SCD_T3_1_CB(const LgIndex_t *I)
{
	TecUtilLockStart(AddOnID);
	TRACE1("Scale (SCRad_SCD_T3_1) Value Changed on drag,  New value is: %d\n", *I);
	ResizeSpheresCallback(false);
	TecUtilLockFinish(AddOnID);
}


/**
 */
static void TGLRadAll_TOG_T3_1_CB(const LgIndex_t *I)
{
	TecUtilLockStart(AddOnID);
	TRACE1("Toggle (TGLRadAll_TOG_T3_1) Value Changed,  New value is: %d\n", *I);
	TecUtilLockFinish(AddOnID);
}


/**
 */
static void TGLRadAbs_TOG_T3_1_CB(const LgIndex_t *I)
{
	TecUtilLockStart(AddOnID);
	TRACE1("Toggle (TGLRadAbs_TOG_T3_1) Value Changed,  New value is: %d\n", *I);
	ResizeSpheresCallback();
	TecUtilLockFinish(AddOnID);
}


/**
 */
static void SCSDtight_SC_T1_1_CB(const LgIndex_t *I)
{
	TecUtilLockStart(AddOnID);
	TRACE1("Scale (SCSDtight_SC_T1_1) Value Changed,  New value is: %d\n", *I);
	GBAProcessSystemUpdateNumTriangles();
	TecGUILabelSetText(LBLSDtight_LBL_T1_1, to_string(*I).c_str());
	TecUtilLockFinish(AddOnID);
}


/**
 */
static void SCSDtight_SCD_T1_1_CB(const LgIndex_t *I)
{
	TecUtilLockStart(AddOnID);
	TRACE1("Scale (SCSDtight_SCD_T1_1) Value Changed on drag,  New value is: %d\n", *I);
	GBAProcessSystemUpdateNumTriangles();
	TecGUILabelSetText(LBLSDtight_LBL_T1_1, to_string(*I).c_str());
	TecUtilLockFinish(AddOnID);
}


/**
 */
static void TGLSphTest_TOG_T1_1_CB(const LgIndex_t *I)
{
	TecUtilLockStart(AddOnID);
	TRACE1("Toggle (TGLSphTest_TOG_T1_1) Value Changed,  New value is: %d\n", *I);
	TecUtilLockFinish(AddOnID);
}

/**
 */
static void TGL_TOG_T3_1_CB(const LgIndex_t *I)
{
	TecUtilLockStart(AddOnID);
	TRACE1("Toggle (TGL_TOG_T3_1) Value Changed,  New value is: %d\n", *I);
	TecUtilLockFinish(AddOnID);
}

/**
 */
static void TGL_TOG_T1_1_CB(const LgIndex_t *I)
{
	TecUtilLockStart(AddOnID);
	TRACE1("Toggle (TGL_TOG_T1_1) Value Changed,  New value is: %d\n", *I);
	TecUtilLockFinish(AddOnID);
}

/**
 */
static LgIndex_t  TFHSdR_TF_T1_1_CB(const char *S)
{
	LgIndex_t IsOk = 1;
	TecUtilLockStart(AddOnID);
	TRACE1("Text field (TFHSdR_TF_T1_1) Value Changed,  New value is: %s\n", S);
	TecUtilLockFinish(AddOnID);
	return (IsOk);
}

/**
 */
static void TGLINSV_TOG_T3_1_CB(const LgIndex_t *I)
{
	TecUtilLockStart(AddOnID);
	TRACE1("Toggle (TGLINSV_TOG_T3_1) Value Changed,  New value is: %d\n", *I);
	ResultsVarListReload();
	TecUtilLockFinish(AddOnID);
}

/**
 */
static LgIndex_t  TFGBSub_TF_T3_1_CB(const char *S)
{
	LgIndex_t IsOk = 1;
	TecUtilLockStart(AddOnID);
	TRACE1("Text field (TFGBSub_TF_T3_1) Value Changed,  New value is: %s\n", S);
	TecUtilLockFinish(AddOnID);
	return (IsOk);
}

/**
 */
static void TGLRCSf_TOG_T3_1_CB(const LgIndex_t *I)
{
	TecUtilLockStart(AddOnID);
	TRACE1("Toggle (TGLRCSf_TOG_T3_1) Value Changed,  New value is: %d\n", *I);
	TecUtilLockFinish(AddOnID);
}

#include "guibld.cpp"
