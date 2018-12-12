#include "TECADDON.h"
#include "ADDGLBL.h"
#if !defined (MSWIN)
#include <unistd.h>
#endif
#include "GUIDEFS.h"

#include "CSM_DATA_TYPES.h"
#include "GTAUSERINPUT.h"
#include "GTAENGINE.h"

#include "GUICB.h"

#include <armadillo>
using namespace arma;




Set_pa ActiveZones;
MouseButtonMode_e OldMouseMode;
std::vector<vec3> PlanePoints;

/**
*/
static void Dialog1HelpButton_CB(void)
{
	TecUtilLockStart(AddOnID);
	TecUtilDialogMessageBox("On-line Help not available for this dialog.",
		MessageBox_Information);
	TecUtilLockFinish(AddOnID);
}


/**
*/
static void Dialog1CloseButton_CB(void)
{
	TecUtilLockStart(AddOnID); 
	QuitProbing(&ActiveZones, &OldMouseMode, &PlanePoints, true);
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
static void Sidebar1Activate_CB(void)
{
	/*  <<< This function is called when sidebar "Gradient Tube Analysis" is activated >>> */
}


/**
*/
static void Sidebar1Deactivate_CB(void)
{
	/*   <<< This function is called when sidebar "Gradient Tube Analysis" is deactivated >>> */
}


/**
*/
static void PBSelectCPs_BTN_S1_CB(void)
{
	TecUtilLockStart(AddOnID);
	TRACE("Select CPs Button Pushed\n");
	if (!TecGUIDialogIsUp(Dialog1Manager))
		ProbeTest_MenuCB(&ActiveZones, &OldMouseMode, &PlanePoints);
	else 
		QuitProbing(&ActiveZones, &OldMouseMode, &PlanePoints, true);
	TecUtilLockFinish(AddOnID);
}


/**
*/
static void PBClear_BTN_S1_CB(void)
{
	TecUtilLockStart(AddOnID);
	TRACE("Clear Selection Button Pushed\n");
	ClearCPList(&PlanePoints);
	TecUtilLockFinish(AddOnID);
}


/**
*/
static void PBRun_BTN_S1_CB(void)
{
	TecUtilLockStart(AddOnID);
	TRACE("Run Button Pushed\n");
	if (GTARunGTA(&PlanePoints))
		PopulateVarList();
	TecUtilLockFinish(AddOnID);
}


/**
*/
static void PBQuit_BTN_S1_CB(void)
{
	TecUtilLockStart(AddOnID);
	TRACE("Quit Button Pushed\n");
	if (TecGUIDialogIsUp(Dialog1Manager))
		QuitButton_CB(&ActiveZones, &OldMouseMode, &PlanePoints, true);
	else
		QuitButton_CB(&ActiveZones, &OldMouseMode, &PlanePoints, false);

	TecUtilLockFinish(AddOnID);
}


/**
*/
static void MLCPs_MLST_S1_CB(LgIndex_t const *I)
{
	TecUtilLockStart(AddOnID);
	TRACE1("Multi selection list (MLCPs_MLST_S1) item selected,  First Item is: %d\n", *I);
	TecUtilLockFinish(AddOnID);
}

// char *OptSelectVar_OPT_S1_List = "Option 1,Option 2,Option 3";
// /**
// */
// static void OptSelectVar_OPT_S1_CB(LgIndex_t const *I)
// {
// 	TecUtilLockStart(AddOnID);
// 	TRACE1("Option Menu (OptSelectVar_OPT_S1) value changed,  New value is: %d\n", *I);
// 	TecUtilLockFinish(AddOnID);
// }
 

char *OPTCutoffVar_OPT_S1_List = "Option 1,Option 2,Option 3";

/**
*/
static void OPTCutoffVar_OPT_S1_CB(LgIndex_t const *I)
{
	TecUtilLockStart(AddOnID);
	TRACE1("Option Menu (OPTCutoffVar_OPT_S1) value changed,  New value is: %d\n", *I);
	TecUtilLockFinish(AddOnID);
}

/**
*/
static void TGLDeleteZon_TOG_S1_CB(LgIndex_t const *I)
{
	TecUtilLockStart(AddOnID);
	TRACE1("Toggle (TGLDeleteZon_TOG_S1) Value Changed,  New value is: %d\n", *I);
	TecUtilLockFinish(AddOnID);
}

/**
*/
static void TGLHideZones_TOG_S1_CB(LgIndex_t const *I)
{
	TecUtilLockStart(AddOnID);
	TRACE1("Toggle (TGLHideZones_TOG_S1) Value Changed,  New value is: %d\n", *I);
	TecUtilLockFinish(AddOnID);
}

/**
*/
static void TGLShowOnly_TOG_S1_CB(LgIndex_t const *I)
{
	TecUtilLockStart(AddOnID);
	TRACE1("Toggle (TGLShowOnly_TOG_S1) Value Changed,  New value is: %d\n", *I);
	TecUtilLockFinish(AddOnID);
}

/**
*/
static LgIndex_t  TFCutoff_TF_S1_CB(char const *S)
{
	LgIndex_t IsOk = 1;
	TecUtilLockStart(AddOnID);
	TRACE1("Text field (TFCutoff_TF_S1) Value Changed,  New value is: %s\n", S);
	TecUtilLockFinish(AddOnID);
	return (IsOk);
}

/**
*/
static LgIndex_t  TFRadOffVal_TF_S1_CB(char const *S)
{
	LgIndex_t IsOk = 1;
	TecUtilLockStart(AddOnID);
	TRACE1("Text field (TFRadOffVal_TF_S1) Value Changed,  New value is: %s\n", S);
	TecUtilLockFinish(AddOnID);
	return (IsOk);
}


/**
*/
static void RBCutoffDi_RADIO_S1_CB(LgIndex_t const *I)
{
	TecUtilLockStart(AddOnID);
	TRACE1("RadioBox (RBCuroffDi_RADIO_S1) Value Changed,  New value is: %d\n", *I);
	TecUtilLockFinish(AddOnID);
}


char *OptSliceAngl_OPT_S1_List = "0.01,0.05,0.1,0.25,0.5,1,2,2.5,3,5,10,15,30,45,60,90,180";



/**
*/
static void OptSliceAngl_OPT_S1_CB(LgIndex_t const *I)
{
	TecUtilLockStart(AddOnID);
	TRACE1("Option Menu (OptSliceAngl_OPT_S1) value changed,  New value is: %d\n", *I);
	TecUtilLockFinish(AddOnID);
}


char *OptStreamAng_OPT_S1_List = "0.01,0.05,0.1,0.25,0.5,1,2,2.5,3,5,10,15,30,45,60,90";



/**
*/
static void OptStreamAng_OPT_S1_CB(LgIndex_t const *I)
{
	TecUtilLockStart(AddOnID);
	TRACE1("Option Menu (OptStreamAng_OPT_S1) value changed,  New value is: %d\n", *I);
	TecUtilLockFinish(AddOnID);
}

/**
*/
static LgIndex_t  TFStBegAng_TF_S1_CB(char const *S)
{
	LgIndex_t IsOk = 1;
	TecUtilLockStart(AddOnID);
	TRACE1("Text field (TFStBegAng_TF_S1) Value Changed,  New value is: %s\n", S);
	TecUtilLockFinish(AddOnID);
	return (IsOk);
}


/**
*/
static LgIndex_t  TFSlcBegAng_TF_S1_CB(char const *S)
{
	LgIndex_t IsOk = 1;
	TecUtilLockStart(AddOnID);
	TRACE1("Text field (TFSlcBegAng_TF_S1) Value Changed,  New value is: %s\n", S);
	TecUtilLockFinish(AddOnID);
	return (IsOk);
}


/**
*/
static LgIndex_t  TFSlcEndAng_TF_S1_CB(char const *S)
{
	LgIndex_t IsOk = 1;
	TecUtilLockStart(AddOnID);
	TRACE1("Text field (TFSlcEndAng_TF_S1) Value Changed,  New value is: %s\n", S);
	TecUtilLockFinish(AddOnID);
	return (IsOk);
}


/**
*/
static LgIndex_t  TFStEndAng_TF_S1_CB(char const *S)
{
	LgIndex_t IsOk = 1;
	TecUtilLockStart(AddOnID);
	TRACE1("Text field (TFStEndAng_TF_S1) Value Changed,  New value is: %s\n", S);
	TecUtilLockFinish(AddOnID);
	return (IsOk);
}


/**
*/
static void MLSelectVar_MLST_S1_CB(LgIndex_t const *I)
{
	TecUtilLockStart(AddOnID);
	TRACE1("Multi selection list (MLSelectVar_MLST_S1) item selected,  First Item is: %d\n", *I);
	TecUtilLockFinish(AddOnID);
}


/**
*/
static void SCProgress_SC_S1_CB(LgIndex_t const *I)
{
	TecUtilLockStart(AddOnID);
	TRACE1("Scale (SCProgress_SC_S1) Value Changed,  New value is: %d\n", *I);
	TecUtilLockFinish(AddOnID);
}

/**
*/
static void SCProgress_SCD_S1_CB(LgIndex_t const *I)
{
	TecUtilLockStart(AddOnID);
	TRACE1("Scale (SCProgress_SCD_S1) Value Changed on drag,  New value is: %d\n", *I);
	TecUtilLockFinish(AddOnID);
}




#include "guibld.cpp"
