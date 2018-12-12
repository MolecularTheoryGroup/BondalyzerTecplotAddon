#if defined MSWIN
#include "STDAFX.h"
#endif
#include "TECADDON.h"
#include "ADDGLBL.h"
#if !defined (MSWIN)
#include <unistd.h>
#endif
#include "GUIDEFS.h"

#include <sstream>

#include "CSM_DATA_SET_INFO.h"

#include "LoadData.h"

using std::stringstream;



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
	TecGUIDialogDrop(Dialog1Manager);
	TecUtilLockFinish(AddOnID);
}


/**
*/
static void Dialog1Init_CB(void)
{
	TecUtilLockStart(AddOnID);
	/* <<< Add init code (if necessary) here>>> */
// 	TecUtilLockFinish(AddOnID);
}


/**
*/
static void MLT41LoadVa_MLST_D1_CB(LgIndex_t const *I)
{
	TecUtilLockStart(AddOnID);
	TRACE1("Multi selection list (MLT41LoadVa_MLST_D1) item selected,  First Item is: %d\n", *I);
	TecUtilLockFinish(AddOnID);
}


/**
*/
static void Dialog1CancelButton_CB(void)
{
	QuitT41Load();
	TecGUIDialogDrop(Dialog1Manager);
	/* Modal Dialogs must call TecUtilLockStart prior to coming */
	/* up and then call TecUtilLockFinish when the Ok or Cancel */
	/* button is pressed.  Only TecUtilLockFinish is supplied here. */
	TecUtilLockFinish(AddOnID);
}


/**
*/
static void Dialog1OkButton_CB(void)
{
	TecGUIDialogDrop(Dialog1Manager);
	/* Modal Dialogs must call TecUtilLockStart prior to coming */
	/* up and then call TecUtilLockFinish when the Ok or Cancel */
	/* button is pressed.  Only TecUtilLockFinish is supplied here. */
	SelVarsLoadData();
	TecUtilLockFinish(AddOnID);
}

/**
*/
static void Dialog2HelpButton_CB(void)
{
	TecUtilLockStart(AddOnID);
	TecUtilDialogMessageBox("On-line Help not available for this dialog.",
		MessageBox_Information);
	TecUtilLockFinish(AddOnID);
}


/**
*/
static void Dialog2CancelButton_CB(void)
{
	TecGUIDialogDrop(Dialog2Manager);
	/* Modal Dialogs must call TecUtilLockStart prior to coming */
	/* up and then call TecUtilLockFinish when the Ok or Cancel */
	/* button is pressed.  Only TecUtilLockFinish is supplied here. */
	TecUtilLockFinish(AddOnID);
}


/**
*/
static void Dialog2OkButton_CB(void)
{
	TecGUIDialogDrop(Dialog2Manager);
	/* Modal Dialogs must call TecUtilLockStart prior to coming */
	/* up and then call TecUtilLockFinish when the Ok or Cancel */
	/* button is pressed.  Only TecUtilLockFinish is supplied here. */
	NumCellsLoadData();
	TecUtilLockFinish(AddOnID);
}


/**
*/
static void Dialog2Init_CB(void)
{
	/* Modal Dialogs must call TecUtilLockStart prior to coming */
	/* up and then call TecUtilLockFinish when the Ok or Cancel */
	/* button is pressed. */
	TecUtilLockStart(AddOnID);
	/* <<< Add init code (if necessary) here>>> */
	TecGUITextFieldSetLgIndex(XNC_TFS_D2, 1, FALSE);
	TecGUITextFieldSetLgIndex(YNC_TFS_D2, 1, FALSE);
	TecGUITextFieldSetLgIndex(ZNC_TFS_D2, 1, FALSE);
}


/**
*/
static LgIndex_t  XNC_TFS_D2_ValueChanged_CB(char const *S)
{
	LgIndex_t IsOk = 1;
	TecUtilLockStart(AddOnID);
	TRACE1("Spin Control Text field (XNC_TFS_D2) Value Changed,  New value is: %s\n", S);
	SpinValueChangedInt(XNC_TFS_D2);
	TecUtilLockFinish(AddOnID);
	return (IsOk);
}


/**
*/
static void XNC_TFS_D2_ButtonUp_CB(void)
{
	TecUtilLockStart(AddOnID);
	TRACE0("Spin control (XNC_TFS_D2) up callback called.\n");
	SpinButtonInt(XNC_TFS_D2, 1);
	TecUtilLockFinish(AddOnID);
}


/**
*/
static void XNC_TFS_D2_ButtonDown_CB(void)
{
	TecUtilLockStart(AddOnID);
	TRACE0("Spin control (XNC_TFS_D2) down callback called.\n");
	SpinButtonInt(XNC_TFS_D2, -1);
	TecUtilLockFinish(AddOnID);
}


/**
*/
static LgIndex_t  YNC_TFS_D2_ValueChanged_CB(char const *S)
{
	LgIndex_t IsOk = 1;
	TecUtilLockStart(AddOnID);
	TRACE1("Spin Control Text field (YNC_TFS_D2) Value Changed,  New value is: %s\n", S);
	SpinValueChangedInt(YNC_TFS_D2);
	TecUtilLockFinish(AddOnID);
	return (IsOk);
}


/**
*/
static void YNC_TFS_D2_ButtonUp_CB(void)
{
	TecUtilLockStart(AddOnID);
	TRACE0("Spin control (YNC_TFS_D2) up callback called.\n");
	SpinButtonInt(YNC_TFS_D2, 1);
	TecUtilLockFinish(AddOnID);
}


/**
*/
static void YNC_TFS_D2_ButtonDown_CB(void)
{
	TecUtilLockStart(AddOnID);
	TRACE0("Spin control (YNC_TFS_D2) down callback called.\n");
	SpinButtonInt(YNC_TFS_D2, -1);
	TecUtilLockFinish(AddOnID);
}


/**
*/
static LgIndex_t  ZNC_TFS_D2_ValueChanged_CB(char const *S)
{
	LgIndex_t IsOk = 1;
	TecUtilLockStart(AddOnID);
	TRACE1("Spin Control Text field (ZNC_TFS_D2) Value Changed,  New value is: %s\n", S);
	SpinValueChangedInt(ZNC_TFS_D2);
	TecUtilLockFinish(AddOnID);
	return (IsOk);
}


/**
*/
static void ZNC_TFS_D2_ButtonUp_CB(void)
{
	TecUtilLockStart(AddOnID);
	TRACE0("Spin control (ZNC_TFS_D2) up callback called.\n");
	SpinButtonInt(ZNC_TFS_D2, 1);
	TecUtilLockFinish(AddOnID);
}


/**
*/
static void ZNC_TFS_D2_ButtonDown_CB(void)
{
	TecUtilLockStart(AddOnID);
	TRACE0("Spin control (ZNC_TFS_D2) down callback called.\n");
	SpinButtonInt(ZNC_TFS_D2, -1);
	TecUtilLockFinish(AddOnID);
}





/**
*/
static void ExportBTN_BTN_T1_1_CB(void)
{
	TecUtilLockStart(AddOnID);
	TRACE("Export densf script Button Pushed\n");
	MakeDensfScriptForZones();
	TecUtilLockFinish(AddOnID);
}


/**
*/
static void ImportBTN_BTN_T1_1_CB(void)
{
	TecUtilLockStart(AddOnID);
	TRACE("Import Tape41 files Button Pushed\n");
	ImportAdditionalTape41Files(TRUE);
	TecUtilLockFinish(AddOnID);
}


/**
*/
static void ZoneList_MLST_T1_1_CB(LgIndex_t const *I)
{
	TecUtilLockStart(AddOnID);
	TRACE1("Multi selection list (ZoneList_MLST_T1_1) item selected,  First Item is: %d\n", *I);
	TecUtilLockFinish(AddOnID);
}


/**
*/
static void Dialog3CloseButton_CB(void)
{
	TecUtilLockStart(AddOnID);
	TecGUIDialogDrop(Dialog3Manager);
	TecUtilLockFinish(AddOnID);
}


/**
*/
static void Dialog3Init_CB(void)
{
	TecUtilLockStart(AddOnID);

	EntIndex_t NumVars, NumZones;

	Boolean_t IsOk = TecUtilDataSetGetInfo(NULL, &NumZones, &NumVars);

	if (IsOk){
		TecGUIListDeleteAllItems(ZoneList_MLST_T1_1);
	}

	EntIndex_t VolZoneNum = ZoneNumByName("Full");

	for (EntIndex_t ZoneNum = 1; ZoneNum <= NumZones && IsOk; ++ZoneNum){
		if (ZoneNum != VolZoneNum){
			ZoneName_t ZoneName;
			if (TecUtilZoneGetName(ZoneNum, &ZoneName)){
				TecGUIListAppendItem(ZoneList_MLST_T1_1, ZoneName);
			}
			TecUtilStringDealloc(&ZoneName);
		}
	}

	stringstream InfoLabel;

	InfoLabel << "This tool is for replacing data values in any zone" << '\n' <<
		"other than the full volume zone with analytical" << '\n' <<
		"values calculated using the 'densf' utility." << '\n' <<
		'\n' <<
		"To use : " << '\n' <<
		"1. Select the zones for which you'd like analytical" << '\n' <<
		"   values in tab 1, then press 'Export shell script'" << '\n' <<
		"2. Run the shell script using Terminal on a Mac or " << '\n' <<
		"   Linux machine (i.e. 'bash script.sh /path/to/tape21') " << '\n' <<
		"   from the same directory as the script file, using " << '\n' <<
		"   the full path to the Tape21 file used to make the " << '\n' <<
		"   Tecplot file as as the only argument to the script." << '\n' <<
		"   This will make a Tape41 file for each zone that was" << '\n' <<
		"   selected." << '\n' <<
		"3. Once the Tape41 files are made, move them back to the" << '\n' <<
		"   windows machine." << '\n' <<
		"4. Press 'Import Tape41 files' and select ALL of the Tape41" << '\n' <<
		"   files. They will be imported, overwriting the variable" << '\n' <<
		"   values for the zones that were exported.";



	TecGUILabelSetText(InfoText_LBL_T2_1, InfoLabel.str().c_str());

	TecUtilLockFinish(AddOnID);
}


/**
*/
static void TAB1_TBA_D3_CB(LgIndex_t const *I)
{
	TecUtilLockStart(AddOnID);
	TRACE0("Activate callback for tab (TAB1_TBA_D3) called\n");
	TecUtilLockFinish(AddOnID);
}


/**
*/
static void TAB1_TBD_D3_CB(LgIndex_t const *I)
{
	TecUtilLockStart(AddOnID);
	TRACE0("Deactivate callback for tab (TAB1_TBD_D3) called\n");
	TecUtilLockFinish(AddOnID);
}




#include "guibld.cpp"