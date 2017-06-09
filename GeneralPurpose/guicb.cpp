#include "TECADDON.h"
#include "ADDGLBL.h"
#if !defined (MSWIN)
#include <unistd.h>
#endif
#include "GUIDEFS.h"
#include "PROBETEST.h"




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
static void Dialog2CloseButton_CB(void)
{
	TecUtilLockStart(AddOnID);
	TecGUIDialogDrop(Dialog2Manager);
	TecUtilLockFinish(AddOnID);
}


/**
*/
static void Dialog2Init_CB(void)
{
	TecUtilLockStart(AddOnID);
	/* <<< Add init code (if necessary) here>>> */
	TecUtilLockFinish(AddOnID);
}


/**
*/
static void StartProbe_CB(void)
{
	TecUtilLockStart(AddOnID);
	TRACE("Click and then select point Button Pushed\n");
	ProbeTest_MenuCB();
	TecUtilLockFinish(AddOnID);
}




#include "guibld.cpp"
