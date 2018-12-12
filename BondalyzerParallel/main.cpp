/*
*****************************************************************
*****************************************************************
*******                                                  ********
*******    (C) Copyright 1988-2013 by TECPLOT INC.       *******
*******                                                  ********
*****************************************************************
*****************************************************************
*/

#include "TECADDON.h"
#include "ADDGLBL.h"
#include "ADDONVER.h"
#include "GUIDEFS.h"
// #include "SETCHEMBONDSTYLE.h"
#include "CSM_CALC_VARS.h"
#include "CSM_GUI.h"
#include "ENGINE.h"
#include <string.h>
#include <string>
#include <sstream>

using std::string;
using std::to_string;



AddOn_pa AddOnID;



/**
 * This function is called when the
 * $!ADDONCOMMAND macro command is
 * processed.
 */
static Boolean_t STDCALL MacroCommandCallback(char *MacroCommandString,  /* IN */
                                              char **ErrMsg)             /* OUT (only if returning FALSE) */
{

    Boolean_t IsOk = TRUE;

    /*
     * MacroCommandString is the add-on macro command string needing processing.
     *
     * *ErrMsg is an error message string which must be allocated and set by this
     * function if and only if the return value is FALSE.
     */

    TecUtilLockStart(AddOnID);

    /*
     * TODO: Process the macro command.
     *
     * Example:
     *
     * $!ADDONCOMMAND ADDONID='Scalar Gradient Topology' COMMAND='MYCOMMAND'
     */

    if (!strcmp(MacroCommandString, "MYCOMMAND")) /* For example */
    {
        /* IsOk = ProcessMacroCommand_MYCOMMAND(); */
    }

    if (!IsOk)
    {
        /*
         * Some kind of error, so inform the user about it.
         */

        *ErrMsg = TecUtilStringAlloc(1000, "String for Error Message");
        strcpy(*ErrMsg, "Error processing macro command");
    }
    else
    {
        /* Ignore the *ErrMsg parameter */
    }

    TecUtilLockFinish(AddOnID);
    return (IsOk);
}

/**
 */
static void STDCALL StateChangeCallback(StateChange_e StateChange)
{

    switch (StateChange)
    {
            /*
             * This function will be called by Tecplot
             * each time a state change occurs.
             *
             *
             * NOTE:
             *
             * Some State changes also have some supplemental "state"
             * information that can be retrieved if you desire.
             * Comments in the case statements below identify these
             * state changes.  To retrieve the supplemental information
             * use the functions TecUtilStateChangeGetXXXXX. You may
             * only call these functions during the scope of this
             * callback.  Once control returns from this call the
             * supplemental information will become unaccessible.
             *
             */

            /*   State Change                Supplemental information */
        case StateChange_VarsAltered:     /* set of altered variables */
        case StateChange_VarsAdded:       /* set of added variables */
        case StateChange_ZonesDeleted:    /* set of deleted zones */
        case StateChange_ZonesAdded:      /* set of added zones */
        case StateChange_NodeMapsAltered: /* set of node maps altered */
		case StateChange_MouseModeUpdate: /* the new mouse mode */
		{
											  MouseButtonMode_e m = TecUtilMouseGetCurrentMode();
											  if (m != MouseButtonMode_Invalid && m != MouseButtonMode_Select && m != MouseButtonMode_Probe){
												  //CSMGUIDeleteCPLabels(&AddOnID);
												  CSMGuiLabelSelectedPoints(&AddOnID);
											  }
		}
			break;
        case StateChange_Style:           /* Style Parameters P1,P2,P3,P4,P5,P6 */
        case StateChange_View:            /* View action (View_e) */
        case StateChange_Streamtrace:     /* Streamtrace action (Streamtrace_e) */
        case StateChange_AuxDataAltered:  /* Name, Auxiliary Location (AuxDataLocation_e), Var/Zone/or Map Num */
        case StateChange_AuxDataAdded:    /* Name, Auxiliary Location (AuxDataLocation_e), Var/Zone/or Map Num */
        case StateChange_AuxDataDeleted:  /* Name, Auxiliary Location (AuxDataLocation_e), Var/Zone/or Map Num */
        case StateChange_VarsDeleted:     /* set of deleted variables (zero based set) */
        case StateChange_VariableLockOn:  /* Locker name, Variable Num, VarLockMode */
        case StateChange_VariableLockOff: /* Unlocker name, Variable Num */
        case StateChange_DataSetLockOn:   /* Locker name */
        case StateChange_DataSetLockOff:  /* Unlocker name */

            /* State changes which do not have any supplemental "state" information. */
        case StateChange_TecplotIsInitialized:/* Tecplot is finished initializing */
        case StateChange_FrameDeleted:        /* A frame was delete */
        case StateChange_NewTopFrame:         /* A new frame has become the current frame */
        case StateChange_Text:                /* One or more text elements has changed */
        case StateChange_Geom:                /* One or more geometry elements has changed */
        case StateChange_DataSetReset:        /* A new dataset has been loaded */
        case StateChange_NewLayout:           /* The current layout has been cleared and reset */
        case StateChange_CompleteReset:       /* Anything could have happened */
        case StateChange_LineMapAssignment:   /* A line mapping definition has been altered (includes zone and axis information) */
        case StateChange_ContourLevels:       /* The contour levels have been altered */
        case StateChange_ModalDialogLaunch:   /* A modal dialog has been launched */
        case StateChange_ModalDialogDismiss:  /* A modal dialog has been dismissed */
        case StateChange_QuitTecplot:         /* Tecplot is about to exit */
        case StateChange_ZoneName:            /* The name of a zone has been altered */
        case StateChange_VarName:             /* The name of a variable has been altered */
        case StateChange_LineMapName:           /* The name of an X-Y mapping has been altered */
        case StateChange_LineMapAddDeleteOrReorder: /* The set of existing X-Y mappings has been altered */
        case StateChange_ColorMap:            /* The color mapping has been altered */
        case StateChange_ContourVar:          /* The contour variable has been reassigned */
        case StateChange_NewAxisVariables:    /* The axis variables have been reassigned */
        case StateChange_PickListCleared:     /* All picked objects are unpicked */
        case StateChange_PickListGroupSelect: /* A group of objects has been added to the pick list */
        case StateChange_PickListSingleSelect:/* A single object has been added to or removed from the pick list */
        case StateChange_PickListStyle:       /* An action has been performed on all of the objects in the pick list */
        case StateChange_DataSetFileName:     /* The current data set has been saved to a file */
        case StateChange_DataSetTitle:        /* The current data set title has been changed */
        case StateChange_DrawingInterrupted:  /* The user has interrupted the drawing */
        case StateChange_ImageExported:       /* An image frame was exported */


            /* Version 9 and later Note: If you are using modeless dialogs, you should
               trap the following state changes and take appropriate
               action when print preview is launched and dismissed.

               Usually you will either disable or close your dialog
               when print preview is launched. */

        case StateChange_PrintPreviewLaunch:  /* Modeless dialogs should close or disable themselves */
        case StateChange_PrintPreviewDismiss: /* Modeless dialogs can re-launch or enable themselves */


        case StateChange_SuspendInterface:    /* Replaces StateChange_DrawGraphicsOn */
        case StateChange_UnsuspendInterface:  /* Replaces StateChange_DrawGraphicsOff */
        {
            /* TODO: Add code to handle state changes.... */
        } break;
        default: break;
    } /* end switch */
}

static void STDCALL GradientPathToolCallback(void)
{
	TecUtilLockStart(AddOnID);
	if (TecUtilDataSetIsAvailable())
	{
		GradientPathToolGetUserInfo();
	}
	else
	{
		TecUtilDialogErrMsg("No data set in current frame.");
	}
	TecUtilLockFinish(AddOnID);
}


static void STDCALL CalcVarsMenuCallback(void)
{
	TecUtilLockStart(AddOnID);
	if (TecUtilDataSetIsAvailable())
	{
		BuildDialog2(MAINDIALOGID);
		TecGUIDialogLaunch(Dialog2Manager);
	}
	else
	{
		TecUtilDialogErrMsg("No data set in current frame.");
	}
	TecUtilLockFinish(AddOnID);
}

static void STDCALL GaussianBlurMenuCallback(void)
{
	const Boolean_t		UseDialog = TRUE;
	TecUtilLockStart(AddOnID);
	if (TecUtilDataSetIsAvailable())
	{
		char* ZoneVarSigma;
		if (TecUtilDialogGetSimpleText("Enter ZoneNum,VarNum,Sigma value for gaussian blur", "1,4,3", &ZoneVarSigma)){

			std::stringstream SS;
			SS << ZoneVarSigma;
			TecUtilStringDealloc(&ZoneVarSigma);
			vector<string> Strs;
			string Str;

			while (getline(SS, Str, ',')){
				Strs.push_back(Str);
			}

			if (Strs.size() == 3){
				int ZoneNum = stoi(Strs[0]);
				int VarNum = stoi(Strs[1]);
				double Sigma = stod(Strs[2]);

				char* CStr;
				TecUtilVarGetName(VarNum, &CStr);
				Str = CStr;
				Str += "_Blur_" + to_string(Sigma);

				GaussianBlur(FALSE, AddOnID, ZoneNum, VarNum, Str, Sigma);
			}
		}
	}
	else
	{
		TecUtilDialogErrMsg("No data set in current frame.");
	}
	TecUtilLockFinish(AddOnID);
}

static void STDCALL ZoneNameFindReplaceMenuCallback(void)
{
	TecUtilLockStart(AddOnID);
	if (TecUtilDataSetIsAvailable())
	{
		ZoneNameFindReplaceGetUserInfo();
	}
	else
	{
		TecUtilDialogErrMsg("No data set in current frame.");
	}

	TecUtilLockFinish(AddOnID);
}

static void STDCALL VarNameFindReplaceMenuCallback(void)
{
	TecUtilLockStart(AddOnID);
	if (TecUtilDataSetIsAvailable())
	{
		VarNameFindReplaceGetUserInfo();
	}
	else
	{
		TecUtilDialogErrMsg("No data set in current frame.");
	}

	TecUtilLockFinish(AddOnID);
}

static void STDCALL MapVarsToZonesMenuCallback(void)
{
	TecUtilLockStart(AddOnID);
	if (TecUtilDataSetIsAvailable())
	{
		CSMGuiLock();

		MapAllVarsToAllZones(AddOnID);

		CSMGuiUnlock();
	}
	else
	{
		TecUtilDialogErrMsg("No data set in current frame.");
	}

	TecUtilLockFinish(AddOnID);
}

static void STDCALL RefineActiveZonesCallback(void)
{
	TecUtilLockStart(AddOnID);
	if (TecUtilDataSetIsAvailable())
	{
		CSMGuiLock();

		RefineActiveZones();

		CSMGuiUnlock();
	}
	else
	{
		TecUtilDialogErrMsg("No data set in current frame.");
	}
	TecUtilLockFinish(AddOnID);
}

static void STDCALL MakeSliceFomeCPsCallback(void)
{
	TecUtilLockStart(AddOnID);
	if (TecUtilDataSetIsAvailable())
	{
		MakeSliceFromPointSelectionGetUserInfo();
	}
	else
	{
		TecUtilDialogErrMsg("No data set in current frame.");
	}
	TecUtilLockFinish(AddOnID);
}

static void STDCALL MakeSurfaceFromPathZonesMenuCallback(void)
{
	TecUtilLockStart(AddOnID);
	if (TecUtilDataSetIsAvailable())
	{
		MakeSurfaceFromPathZonesGetUserInfo();
	}
	else
	{
		TecUtilDialogErrMsg("No data set in current frame.");
	}

	TecUtilLockFinish(AddOnID);
}

static void STDCALL GradientPathsOnSphereMenuCallback(void)
{
	TecUtilLockStart(AddOnID);
	if (TecUtilDataSetIsAvailable())
	{
		GradientPathsOnSphereGetUserInfo();
	}
	else
	{
		TecUtilDialogErrMsg("No data set in current frame.");
	}

	TecUtilLockFinish(AddOnID);
}

static void STDCALL GetClosedIsoSurfaceFromPointsCallback(void)
{
	TecUtilLockStart(AddOnID);
	if (TecUtilDataSetIsAvailable())
	{
		CSMGuiLock();

		GetClosedIsoSurfaceFromPoints();

		CSMGuiUnlock();
	}
	else
	{
		TecUtilDialogErrMsg("No data set in current frame.");
	}
	TecUtilLockFinish(AddOnID);
}

static void STDCALL GetClosedIsoSurfaceFromNodesCallback(void)
{
	TecUtilLockStart(AddOnID);
	if (TecUtilDataSetIsAvailable())
	{
		CSMGuiLock();

		GetClosedIsoSurfaceFromNodes();

		CSMGuiUnlock();
	}
	else
	{
		TecUtilDialogErrMsg("No data set in current frame.");
	}
	TecUtilLockFinish(AddOnID);
}

static void STDCALL GetAllClosedIsoSurfacesCallback(void)
{
	TecUtilLockStart(AddOnID);
	if (TecUtilDataSetIsAvailable())
	{
		CSMGuiLock();

		GetAllClosedIsoSurfaces();

		CSMGuiUnlock();
	}
	else
	{
		TecUtilDialogErrMsg("No data set in current frame.");
	}
	TecUtilLockFinish(AddOnID);
}

static void STDCALL ConnectCPsCallback(void)
{
	TecUtilLockStart(AddOnID);
	if (TecUtilDataSetIsAvailable())
	{
		CSMGuiLock();

		ConnectCPsGetUserInfo();

		CSMGuiUnlock();
	}
	else
	{
		TecUtilDialogErrMsg("No data set in current frame.");
	}
	TecUtilLockFinish(AddOnID);
}

static void STDCALL BondalyzerBatchCallback(void)
{
	TecUtilLockStart(AddOnID);
	if (TecUtilDataSetIsAvailable())
	{
		BondalyzerGetUserInfo(BondalyzerCalcType_Batch);
	}
	else
	{
		TecUtilDialogErrMsg("No data set in current frame.");
	}
	TecUtilLockFinish(AddOnID);
}

static void STDCALL FindCritPointsCallback(void)
{
	TecUtilLockStart(AddOnID);
	if (TecUtilDataSetIsAvailable())
	{
		BondalyzerGetUserInfo(BondalyzerCalcType_CriticalPoints);
	}
	else
	{
		TecUtilDialogErrMsg("No data set in current frame.");
	}
	TecUtilLockFinish(AddOnID);
}

static void STDCALL DeleteCritPointsCallback(void)
{
	TecUtilLockStart(AddOnID);
	if (TecUtilDataSetIsAvailable())
	{
		DeleteCPsGetUserInfo();
	}
	else
	{
		TecUtilDialogErrMsg("No data set in current frame.");
	}
	TecUtilLockFinish(AddOnID);
}

static void STDCALL ExtractCritPointsCallback(void)
{
	TecUtilLockStart(AddOnID);
	if (TecUtilDataSetIsAvailable())
	{
		ExtractCPsGetUserInfo();
	}
	else
	{
		TecUtilDialogErrMsg("No data set in current frame.");
	}
	TecUtilLockFinish(AddOnID);
}

static void STDCALL CombineCritPointZonesCallback(void)
{
	TecUtilLockStart(AddOnID);
	if (TecUtilDataSetIsAvailable())
	{
		CombineCPZonesGetUserInfo();
	}
	else
	{
		TecUtilDialogErrMsg("No data set in current frame.");
	}
	TecUtilLockFinish(AddOnID);
}

static void STDCALL FindBondPathsCallback(void)
{
	TecUtilLockStart(AddOnID);
	if (TecUtilDataSetIsAvailable())
	{
		BondalyzerGetUserInfo(BondalyzerCalcType_BondPaths);
	}
	else
	{
		TecUtilDialogErrMsg("No data set in current frame.");
	}
	TecUtilLockFinish(AddOnID);
}

static void STDCALL FindRingLinesCallback(void)
{
	TecUtilLockStart(AddOnID);
	if (TecUtilDataSetIsAvailable())
	{
		BondalyzerGetUserInfo(BondalyzerCalcType_RingLines);
	}
	else
	{
		TecUtilDialogErrMsg("No data set in current frame.");
	}
	TecUtilLockFinish(AddOnID);
}

static void STDCALL FindCageNuclearPathsCallback(void)
{
	TecUtilLockStart(AddOnID);
	if (TecUtilDataSetIsAvailable())
	{
		BondalyzerGetUserInfo(BondalyzerCalcType_CageNuclearPaths);
	}
	else
	{
		TecUtilDialogErrMsg("No data set in current frame.");
	}
	TecUtilLockFinish(AddOnID);
}

static void STDCALL FindBondSurfacesCallback(void)
{
	TecUtilLockStart(AddOnID);
	if (TecUtilDataSetIsAvailable())
	{
		BondalyzerGetUserInfo(BondalyzerCalcType_InteratomicSurfaces);
	}
	else
	{
		TecUtilDialogErrMsg("No data set in current frame.");
	}
	TecUtilLockFinish(AddOnID);
}

static void STDCALL FindRingSurfacesCallback(void)
{
	TecUtilLockStart(AddOnID);
	if (TecUtilDataSetIsAvailable())
	{
		BondalyzerGetUserInfo(BondalyzerCalcType_RingSurfaces);
	}
	else
	{
		TecUtilDialogErrMsg("No data set in current frame.");
	}
	TecUtilLockFinish(AddOnID);
}

static void STDCALL DrawEigenvectorArrowsCallback(void)
{
	TecUtilLockStart(AddOnID);
	if (TecUtilDataSetIsAvailable())
	{
		CSMGuiLock();

		DrawEigenvectorArrowsGetUserInfo();

		CSMGuiUnlock();
	}
	else
	{
		TecUtilDialogErrMsg("No data set in current frame.");
	}
	TecUtilLockFinish(AddOnID);
}

static void STDCALL TestFunctionCallback(void)
{
	TecUtilLockStart(AddOnID);
	if (TecUtilDataSetIsAvailable())
	{
		CSMGuiLock();

		TestFunction();

		CSMGuiUnlock();
	}
	else
	{
		TecUtilDialogErrMsg("No data set in current frame.");
	}
	TecUtilLockFinish(AddOnID);
}


static void STDCALL ExtractRSIntersectionsCallback(void)
{
	TecUtilLockStart(AddOnID);
	if (TecUtilDataSetIsAvailable())
	{
		CSMGuiLock();

		ExtractRSIntersectionsGetUserInfo();

		CSMGuiUnlock();
	}
	else
	{
		TecUtilDialogErrMsg("No data set in current frame.");
	}
	TecUtilLockFinish(AddOnID);
}

/**
 * When Tecplot first loads an add-on, it makes a
 * call to initialize the add-on. This function
 * must be named InitTecAddOn, as shown below.
 */
EXPORTFROMADDON void STDCALL InitTecAddOn(void)
{


    /*
     * NOTE:  TecUtilLockOn MUST be used for InitTecAddOn instead
     *        of TecUtilLockStart because AddonID has yet to be
     *        established.  TecUtilLockOn is in effect an "anonymous"
     *        locking of Tecplot (old style).
     */

    TecUtilLockOn();

    /*
     * The function TecUtilAddOnRegister() is the
     * only function that is REQUIRED to be called from
     * the initialization function.
     *
     * The information you give Tecplot by calling
     * this function will show up in the Help/About Add-ons
     * dialog box.
     */

    /*
     * Note that if your add-on requires a specific version of Tecplot,
     * you would check for that here using TecUtilGetTecplotVersion()
     */

    AddOnID = TecUtilAddOnRegister(110,
                                   ADDON_NAME,
                                   "V" ADDON_VERSION"(" TecVersionId") " ADDON_DATE,
                                   "Tecplot");

    /*
     * Initialize the Tecplot GUI Builder libraries.
     */
    InitTGB();


    TecUtilMacroAddCommandCallback(ADDON_NAME,
                                   MacroCommandCallback);
    {
        ArgList_pa ArgList;
        ArgList = TecUtilArgListAlloc();
        TecUtilArgListAppendFunction(ArgList, SV_CALLBACKFUNCTION, (void const *)StateChangeCallback);
        TecUtilArgListAppendInt(ArgList,      SV_STATECHANGEMODE,        StateChangeMode_v100);
        TecUtilArgListAppendInt(ArgList,      SV_STATECHANGECALLBACKAPI, StateChangeCallbackAPI_ChangeOnly);
        TecUtilStateChangeAddCallbackX(ArgList);
        TecUtilArgListDealloc(&ArgList);
    }

	int MenuNum = 1;

	TecUtilMenuAddOption("MTG_Utilities",
		string("Calculate variables").c_str(),
		'\0',
		CalcVarsMenuCallback);

	TecUtilMenuAddOption("MTG_Utilities",
		string("Gradient path tool").c_str(),
		'\0',
		GradientPathToolCallback);

// 	TecUtilMenuAddOption("MTG_Bondalyzer",
// 		string(to_string(MenuNum++) + ". Gaussian blur").c_str(),
// 		'\0',
// 	GaussianBlurMenuCallback);

	TecUtilMenuAddOption("MTG_Utilities",
		string("Map volume zone variables to all zones").c_str(),
		'\0',
		MapVarsToZonesMenuCallback);

	TecUtilMenuAddOption("MTG_Utilities",
		string("Refine active zones").c_str(),
		'\0',
		RefineActiveZonesCallback);

	TecUtilMenuAddOption("MTG_Utilities",
		string("Gradient path shutgun around nuclear/cage CPs").c_str(),
		'\0',
		GradientPathsOnSphereMenuCallback);

	TecUtilMenuAddOption("MTG_Utilities",
		string("Make surface from path zones").c_str(),
		'\0',
		MakeSurfaceFromPathZonesMenuCallback);

	TecUtilMenuAddOption("MTG_Utilities",
		string("Create slice from CPs").c_str(),
		'\0',
		MakeSliceFomeCPsCallback);

	TecUtilMenuAddOption("MTG_Utilities",
		string("Get closed isosurface from points").c_str(),
		'\0',
		GetClosedIsoSurfaceFromPointsCallback);

	TecUtilMenuAddOption("MTG_Utilities",
		string("Get closed isosurface from nodes").c_str(),
		'\0',
		GetClosedIsoSurfaceFromNodesCallback);

	TecUtilMenuAddOption("MTG_Utilities",
		string("Get all closed isosurface").c_str(),
		'\0',
		GetAllClosedIsoSurfacesCallback);

	TecUtilMenuAddOption("MTG_Utilities",
		string("Connect CPs with lines").c_str(),
		'\0',
		ConnectCPsCallback);

	TecUtilMenuAddOption("MTG_Utilities",
		string("Draw arrows to indicate eigenvectors").c_str(),
		'\0',
		DrawEigenvectorArrowsCallback);

	TecUtilMenuAddOption("MTG_Utilities",
		string("Zone name: find and replace").c_str(),
		'\0',
		ZoneNameFindReplaceMenuCallback);

	TecUtilMenuAddOption("MTG_Utilities",
		string("Variable name: find and replace").c_str(),
		'\0',
		VarNameFindReplaceMenuCallback);

	TecUtilMenuAddOption("MTG_Bondalyzer",
		string("Batch analysis").c_str(),
		'\0',
		BondalyzerBatchCallback);

	TecUtilMenuAddOption("MTG_Bondalyzer",
		string("1. Find critical points").c_str(),
		'\0',
		FindCritPointsCallback);


	TecUtilMenuAddOption("MTG_Bondalyzer",
		string("1a. Delete critical point(s)").c_str(),
		'\0',
		DeleteCritPointsCallback);

	TecUtilMenuAddOption("MTG_Bondalyzer",
		string("1b. Extract critical point(s)").c_str(),
		'\0',
		ExtractCritPointsCallback);

	TecUtilMenuAddOption("MTG_Bondalyzer",
		string("1c. Combine critical point zones").c_str(),
		'\0',
		CombineCritPointZonesCallback);

	TecUtilMenuAddOption("MTG_Bondalyzer",
		string("2. Find bond paths").c_str(),
		'\0',
		FindBondPathsCallback);

	TecUtilMenuAddOption("MTG_Bondalyzer",
		string("3. Find ring lines").c_str(),
		'\0',
		FindRingLinesCallback);

	TecUtilMenuAddOption("MTG_Bondalyzer",
		string("4. Find cage-nuclear paths").c_str(),
		'\0',
		FindCageNuclearPathsCallback);

	TecUtilMenuAddOption("MTG_Bondalyzer",
		string("5. Find interatomic surfaces").c_str(),
		'\0',
		FindBondSurfacesCallback);

	TecUtilMenuAddOption("MTG_Bondalyzer",
		string("6. Find ring surfaces").c_str(),
		'\0',
		FindRingSurfacesCallback);

	TecUtilMenuAddOption("MTG_Utilities",
		string("Extract surface/CP intersections").c_str(),
		'\0',
		ExtractRSIntersectionsCallback);

	TecUtilMenuAddOption("MTG_Utilities",
		string("Test function").c_str(),
		'\0',
		TestFunctionCallback);


    /*
     * See note on TecUtilLockOn at start of this function.
     */
    TecUtilLockOff();
}

