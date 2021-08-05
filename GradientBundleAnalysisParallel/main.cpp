
#include <string.h>

#include "TECADDON.h"
#include "ADDGLBL.h"
#include "GUIDEFS.h"
#include "GUICB.h"

#include "CSM_GUI.h"

#include "GBAENGINE.h"
#include "VIEWRESULTS.h"
#include "INTEGRATE.h"

AddOn_pa AddOnID;



/**
 * This function is called when the
 * $!EXTENDEDCOMMAND macro command is
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
	 * $!EXTENDEDCOMMAND COMMANDPROCESSORID='General Purpose Sample' COMMAND='MYCOMMAND'
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
	TecUtilLockStart(AddOnID);
	if (TecGUIDialogIsUp(Dialog1Manager) &&
		((StateChange == StateChange_QuitTecplot) ||
		(TecUtilFrameGetPlotType() != PlotType_XYLine)))
	{
		TecGUIDialogDrop(Dialog1Manager);
	}




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
// 			TecUtilDialogMessageBox("vars altered", MessageBoxType_Information);
			break;
		case StateChange_VarsAdded:       /* set of added variables */
// 			TecUtilDialogMessageBox("vars added", MessageBoxType_Information);
			break;
		case StateChange_ZonesDeleted:    /* set of deleted zones */
// 			TecUtilDialogMessageBox("zones deleted", MessageBoxType_Information);
			break;
		case StateChange_ZonesAdded:      /* set of added zones */
// 			TecUtilDialogMessageBox("zones added", MessageBoxType_Information);
			break;
		case StateChange_NodeMapsAltered: /* set of node maps altered */
// 			TecUtilDialogMessageBox("node maps altered", MessageBoxType_Information);
			break;
		case StateChange_MouseModeUpdate: /* the new mouse mode */
			CSMGUILock();
			GBAProcessSystemDeleteCPLabels();
			CSMGUIUnlock();
// 			TecUtilDialogMessageBox("mouse mode update", MessageBoxType_Information);
			break;
		case StateChange_Style:           /* Style Parameters P1,P2,P3,P4,P5,P6 */
// 			TecUtilDialogMessageBox("style", MessageBoxType_Information);
// 			GBAReloadDialog();
			break;
		case StateChange_View:            /* View action (View_e) */
// 			TecUtilDialogMessageBox("view", MessageBoxType_Information);
// 			GBAReloadDialog();
			break;
		case StateChange_Streamtrace:     /* Streamtrace action (Streamtrace_e) */
// 			TecUtilDialogMessageBox("streamtrace", MessageBoxType_Information);
			break;
		case StateChange_AuxDataAltered:  /* Name, Auxiliary Location (AuxDataLocation_e), Var/Zone/or Map Num */
// 			TecUtilDialogMessageBox("aux data altered", MessageBoxType_Information);
			break;
		case StateChange_AuxDataAdded:    /* Name, Auxiliary Location (AuxDataLocation_e), Var/Zone/or Map Num */
// 			TecUtilDialogMessageBox("aux data added", MessageBoxType_Information);
			break;
		case StateChange_AuxDataDeleted:  /* Name, Auxiliary Location (AuxDataLocation_e), Var/Zone/or Map Num */
// 			TecUtilDialogMessageBox("aux data deleted", MessageBoxType_Information);
			break;
		case StateChange_VarsDeleted:     /* set of deleted variables (zero based set) */
// 			TecUtilDialogMessageBox("vars deleted", MessageBoxType_Information);
			break;
		case StateChange_VariableLockOn:  /* Locker name, Variable Num, VarLockMode */
// 			TecUtilDialogMessageBox("variable lock on", MessageBoxType_Information);
			break;
		case StateChange_VariableLockOff: /* Unlocker name, Variable Num */
// 			TecUtilDialogMessageBox("variable lock off", MessageBoxType_Information);
			break;
		case StateChange_DataSetLockOn:   /* Locker name */
// 			TecUtilDialogMessageBox("data set lock on", MessageBoxType_Information);
			break;
		case StateChange_DataSetLockOff:  /* Unlocker name */
// 			TecUtilDialogMessageBox("data set lock off", MessageBoxType_Information);
			break;

			/* State changes which do not have any supplemental "state" information. */
		case StateChange_TecplotIsInitialized:/* Tecplot is finished initializing */
// 			TecUtilDialogMessageBox("tecplot initialized", MessageBoxType_Information);
			break;
		case StateChange_FrameDeleted:        /* A frame was delete */
// 			TecUtilDialogMessageBox("frame deleted", MessageBoxType_Information);
			break;
		case StateChange_NewTopFrame:         /* A new frame has become the current frame */
// 			TecUtilDialogMessageBox("new top frame", MessageBoxType_Information);
			break;
		case StateChange_Text:                /* One or more text elements has changed */
// 			TecUtilDialogMessageBox("text", MessageBoxType_Information);
// 			GBAReloadDialog();
			break;
		case StateChange_Geom:                /* One or more geometry elements has changed */
// 			TecUtilDialogMessageBox("geom", MessageBoxType_Information);
			break;
		case StateChange_DataSetReset:        /* A new dataset has been loaded */
			// 			TecUtilDialogMessageBox("dataset reset", MessageBoxType_Information);
			break;
		case StateChange_NewLayout:           /* The current layout has been cleared and reset */
// 			TecUtilDialogMessageBox("new layout", MessageBoxType_Information);
			break;
		case StateChange_CompleteReset:       /* Anything could have happened */
			// 			TecUtilDialogMessageBox("complete reset", MessageBoxType_Information);
			break;
		case StateChange_LineMapAssignment:   /* A line mapping definition has been altered (includes zone and axis information) */
// 			TecUtilDialogMessageBox("line map assignment", MessageBoxType_Information);
			break;
		case StateChange_ContourLevels:       /* The contour levels have been altered */
// 			TecUtilDialogMessageBox("contour levels", MessageBoxType_Information);
// 			GBAReloadDialog();
			break;
		case StateChange_ModalDialogLaunch:   /* A modal dialog has been launched */
// 			TecUtilDialogMessageBox("molal dialog launch", MessageBoxType_Information);
			break;
		case StateChange_ModalDialogDismiss:  /* A modal dialog has been dismissed */
// 			TecUtilDialogMessageBox("modal dialog dismiss", MessageBoxType_Information);
			break;
		case StateChange_QuitTecplot:         /* Tecplot is about to exit */
// 			TecUtilDialogMessageBox("quit tecplot", MessageBoxType_Information);
			break;
		case StateChange_ZoneName:            /* The name of a zone has been altered */
// 			TecUtilDialogMessageBox("zone name", MessageBoxType_Information);
			break;
		case StateChange_VarName:             /* The name of a variable has been altered */
// 			TecUtilDialogMessageBox("var name", MessageBoxType_Information);
			break;
		case StateChange_LineMapName:           /* The name of an X-Y mapping has been altered */
// 			TecUtilDialogMessageBox("line map name", MessageBoxType_Information);
			break;
		case StateChange_LineMapAddDeleteOrReorder: /* The set of existing X-Y mappings has been altered */
// 			TecUtilDialogMessageBox("line map add delete or reorder", MessageBoxType_Information);
			break;
		case StateChange_ColorMap:            /* The color mapping has been altered */
// 			TecUtilDialogMessageBox("color map", MessageBoxType_Information);
			break;
		case StateChange_ContourVar:          /* The contour variable has been reassigned */
// 			TecUtilDialogMessageBox("contour var", MessageBoxType_Information);
			break;
		case StateChange_NewAxisVariables:    /* The axis variables have been reassigned */
// 			TecUtilDialogMessageBox("new axis variables", MessageBoxType_Information);
			break;
		case StateChange_PickListCleared:     /* All picked objects are unpicked */
// 			TecUtilDialogMessageBox("pick list cleared", MessageBoxType_Information);
			break;
		case StateChange_PickListGroupSelect: /* A group of objects has been added to the pick list */
// 			TecUtilDialogMessageBox("pick list group select", MessageBoxType_Information);
			break;
		case StateChange_PickListSingleSelect:/* A single object has been added to or removed from the pick list */
// 			TecUtilDialogMessageBox("pick list single select", MessageBoxType_Information);
			break;
		case StateChange_PickListStyle:       /* An action has been performed on all of the objects in the pick list */
// 			TecUtilDialogMessageBox("pick list style", MessageBoxType_Information);
			break;
		case StateChange_DataSetFileName:     /* The current data set has been saved to a file */
// 			TecUtilDialogMessageBox("data set file name", MessageBoxType_Information);
			break;
		case StateChange_DataSetTitle:        /* The current data set title has been changed */
// 			TecUtilDialogMessageBox("data set title", MessageBoxType_Information);
			break;
		case StateChange_DrawingInterrupted:  /* The user has interrupted the drawing */
// 			TecUtilDialogMessageBox("drawing interrupted", MessageBoxType_Information);
			break;
		case StateChange_ImageExported:       /* An image frame was exported */
// 			TecUtilDialogMessageBox("image exported", MessageBoxType_Information);
			break;
		case StateChange_PageDeleted:         /* A page was deleted */
// 			TecUtilDialogMessageBox("page deleted", MessageBoxType_Information);
			break;
		case StateChange_NewTopPage:          /* A different page was made the top page */
// 			TecUtilDialogMessageBox("new top page", MessageBoxType_Information);
			break;


			/* Version 9 and later Note: If you are using modeless dialogs, you should
			   trap the following state changes and take appropriate
			   action when print preview is launched and dismissed.

			   Usually you will either disable or close your dialog
			   when print preview is launched. */

		case StateChange_PrintPreviewLaunch:  /* Modeless dialogs should close or disable themselves */
// 			TecUtilDialogMessageBox("print preview launch", MessageBoxType_Information);
			break;
		case StateChange_PrintPreviewDismiss: /* Modeless dialogs can re-launch or enable themselves */
// 			TecUtilDialogMessageBox("print preview dismiss", MessageBoxType_Information);
			break;


		case StateChange_SuspendInterface:    /* Replaces StateChange_DrawGraphicsOn */
// 			TecUtilDialogMessageBox("suspend interface", MessageBoxType_Information);
			break;
		case StateChange_UnsuspendInterface:  /* Replaces StateChange_DrawGraphicsOff */
// 			TecUtilDialogMessageBox("unsuspend interface", MessageBoxType_Information);
			break;
		{
			/* TODO: Add code to handle state changes.... */
		} break;
		default: break;
	} /* end switch */
	TecUtilLockFinish(AddOnID);
}

/**
 */
static void STDCALL MenuCallback(void)
{
	TecUtilLockStart(AddOnID);
	BuildDialog1(MAINDIALOGID);
	if (TecUtilDataSetIsAvailable() && GBAProcessSystemPrepareGUI())
		TecGUIDialogLaunch(Dialog1Manager);
	else
		TecUtilDialogErrMsg("Load a bondalyzed file first.");
	TecUtilLockFinish(AddOnID);
}

/**
 */
static void STDCALL CreateCircularGBsMenuCallback(void)
{
	TecUtilLockStart(AddOnID);
	BuildDialog1(MAINDIALOGID);
	if (TecUtilDataSetIsAvailable())
		CreateCircularGBsGetUserInfo();
	else
		TecUtilDialogErrMsg("Load a bondalyzed file first.");
	TecUtilLockFinish(AddOnID);
}

static void STDCALL GPTestMenuCB(void)
{
	TecUtilLockStart(AddOnID);
	if (TecUtilDataSetIsAvailable()){
		GradPathTest();
	}
	else
		TecUtilDialogErrMsg("Load a bondalyzed file first.");
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
								   "Tecplot, Inc.");

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

	TecUtilMenuAddOption("MTG_Bondalyzer",
						 "Gradient Bundle Analysis",
						 '\0',
						 MenuCallback);
	TecUtilMenuAddOption("MTG_Utilities",
						"Create Circular GBs on GBA sphere",
						'\0',
						CreateCircularGBsMenuCallback);
// 	TecUtilMenuAddOption("MTG_Bondalyzer",
// 						"GP Test",
// 						'\0',
// 						GPTestMenuCB);


	/*
	 * See note on TecUtilLockOn at start of this function.
	 */
	TecUtilLockOff();
}

