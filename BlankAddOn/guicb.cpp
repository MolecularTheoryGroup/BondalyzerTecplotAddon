#include "TECADDON.h"
#include "ADDGLBL.h"
#if !defined (MSWIN)
#include <unistd.h>
#endif
#include "GUIDEFS.h"

#include "MAINFUNCTION.h"

#include "GUICB.h"


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
	TecUtilLockFinish(AddOnID);
}


/**
*/
static void BTN1_BTN_D1_CB(void)
{
	TecUtilLockStart(AddOnID);
	TRACE("Push Button Button Pushed\n");
	MainFunction();
	TecUtilLockFinish(AddOnID);
}



#include "guibld.cpp"
