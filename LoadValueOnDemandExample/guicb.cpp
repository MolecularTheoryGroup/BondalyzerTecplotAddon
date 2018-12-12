#include "ADDGLBL.h"
#if !defined (MSWIN)
#include <unistd.h>
#endif
#include "GUIDEFS.h"
#include "core/LoadData.h"

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
static void Dialog1CancelButton_CB(void)
{
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

    FILE     *MyFile;
    Boolean_t IsOk = TRUE;
    Boolean_t Result;
    char     *FileName;

    /* First check for a valid file name. */
    FileName = TecGUITextFieldGetString(FileName_TF_D1);
    if (VALID_NON_ZERO_LEN_STR(FileName))
    {
        MyFile = fopen(FileName, "rb");
        if (MyFile == NULL)
        {
            TecUtilDialogErrMsg("Bad file name from OK button.");
            IsOk = FALSE;
        }
        else
        {
            fclose(MyFile);
        }
    }
    else
    {
        TecUtilDialogErrMsg("Enter a file name.");
        IsOk = FALSE;
    }


    if (IsOk == TRUE)
    {
        StringList_pa  Instructions;
        Instructions = TecUtilStringListAlloc();
        TecUtilStringListAppendString(Instructions, "STANDARDSYNTAX");
        TecUtilStringListAppendString(Instructions, "1.0");
        TecUtilStringListAppendString(Instructions, "FILENAME");
        TecUtilStringListAppendString(Instructions, FileName);
        TecUtilStringDealloc(&FileName);
        Result = LoaderCallback(Instructions);
        TecUtilStringListDealloc(&Instructions);
        if (Result)
        {
            TecGUIDialogDrop(Dialog1Manager);
            TecUtilLockFinish(AddOnID);
        }
        else
        {
            TecUtilDialogErrMsg("Error loading the file.");
        }
    }
}


/**
*/
static void Dialog1Init_CB(void)
{
    /* Modal Dialogs must call TecUtilLockStart prior to coming */
    /* up and then call TecUtilLockFinish when the Ok or Cancel */
    /* button is pressed. */
    TecUtilLockStart(AddOnID);
    /* <<< Add init code (if necessary) here>>> */
}


/**
*/
static LgIndex_t  FileName_TF_D1_CB(char const *S)
{
    LgIndex_t IsOk = 1;
    TecUtilLockStart(AddOnID);
    TRACE1("Text field (FileName_TF_D1) Value Changed,  New value is: %s\n", S);
    TecUtilLockFinish(AddOnID);
    return (IsOk);
}


/**
*/
static void Browse_BTN_D1_CB(void)
{
    char *FName = NULL;
    TecUtilLockStart(AddOnID);
    if (TecUtilDialogGetFileName(SelectFileOption_ReadSingleFile,
                                 &FName,
                                 "Binary",
                                 "",
                                 "*.*"))
    {
        TecGUITextFieldSetString(FileName_TF_D1, FName);
        TecUtilStringDealloc(&FName);
    }

    TRACE("... Button Pushed\n");
    TecUtilLockFinish(AddOnID);
}





/**
*/
static void CreateFile_BTN_D1_CB(void)
{
    TecUtilLockStart(AddOnID);
    BuildDialog2(Dialog1Manager);
    TecGUIDialogLaunch(Dialog2Manager);

    TRACE("Make Test File Button Pushed\n");
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
static LgIndex_t  XStart_TF_D2_CB(char const *S)
{
    double XStart;
    LgIndex_t IsOk = 1;
    TecUtilLockStart(AddOnID);
    TRACE1("Text field (XStart_TF_D2) Value Changed,  New value is: %s\n", S);
    if (TecGUITextFieldGetDouble(XStart_TF_D2, &XStart) != TRUE)
    {
        TecUtilDialogErrMsg("Invalid number.");
        TecGUITextFieldSetDouble(XStart_TF_D2, 1.3, "%f");
    }


    TecUtilLockFinish(AddOnID);
    return (IsOk);
}


/**
*/
static LgIndex_t  YStart_TF_D2_CB(char const *S)
{
    double YStart;
    LgIndex_t IsOk = 1;
    TecUtilLockStart(AddOnID);
    TRACE1("Text field (YStart_TF_D2) Value Changed,  New value is: %s\n", S);
    if (TecGUITextFieldGetDouble(YStart_TF_D2, &YStart) != TRUE)
    {
        TecUtilDialogErrMsg("Invalid number.");
        TecGUITextFieldSetDouble(YStart_TF_D2, 0.6, "%2.1f");
    }


    TecUtilLockFinish(AddOnID);
    return (IsOk);
}


/**
*/
static LgIndex_t  XDelta_TF_D2_CB(char const *S)
{
    double XDelta;
    LgIndex_t IsOk = 1;
    TecUtilLockStart(AddOnID);
    TRACE1("Text field (XDelta_TF_D2) Value Changed,  New value is: %s\n", S);
    if (TecGUITextFieldGetDouble(XDelta_TF_D2, &XDelta) != TRUE)
    {
        TecUtilDialogErrMsg("Invalid number.");
        TecGUITextFieldSetDouble(XDelta_TF_D2, 1.2, "%2.1f");
    }

    if ((XDelta <= 0) || (XDelta >= 99))
    {
        TecUtilDialogErrMsg("Small positive number required.");
        TecGUITextFieldSetDouble(XDelta_TF_D2, 1.2, "%2.1f");
    }
    TecUtilLockFinish(AddOnID);
    return (IsOk);
}


/**
*/
static LgIndex_t  YDelta_TF_D2_CB(char const *S)
{
    double YDelta;
    LgIndex_t IsOk = 1;
    TecUtilLockStart(AddOnID);
    TRACE1("Text field (YDelta_TF_D2) Value Changed,  New value is: %s\n", S);
    if (TecGUITextFieldGetDouble(YDelta_TF_D2, &YDelta) != TRUE)
    {
        TecUtilDialogErrMsg("Invalid number.");
        TecGUITextFieldSetDouble(YDelta_TF_D2, 1.2, "%2.1f");
    }

    if ((YDelta <= 0) || (YDelta >= 99))
    {
        TecUtilDialogErrMsg("Small positive number required.");
        TecGUITextFieldSetDouble(YDelta_TF_D2, 2.1, "%2.1f");
    }
    TecUtilLockFinish(AddOnID);
    return (IsOk);
}



/**
*/
static void Dialog2OkButton_CB(void)
{
    int    IMax   = 5;
    int    JMax   = 5;
    double XStart = 0.1;
    double YStart = 1.3;
    double XDelta = 1.0;
    double YDelta = 2.1;

    TecGUITextFieldGetLgIndex(IMax_TF_D2, &IMax);
    TecGUITextFieldGetLgIndex(JMax_TF_D2, &JMax);
    TecGUITextFieldGetDouble(XStart_TF_D2, &XStart);
    TecGUITextFieldGetDouble(YStart_TF_D2, &YStart);
    TecGUITextFieldGetDouble(XDelta_TF_D2, &XDelta);
    TecGUITextFieldGetDouble(YDelta_TF_D2, &YDelta);

    std::string testFileName = CreateTestFile(IMax, JMax,
                                              XStart, YStart,
                                              XDelta, YDelta);

    if (testFileName.empty())
    {
        TecUtilDialogErrMsg("Unable to open file for writing.");
    }
    else
    {
        TecGUITextFieldSetString(FileName_TF_D1,
                                 testFileName.c_str());
    }

    TecGUIDialogDrop(Dialog2Manager);
    /* Modal Dialogs must call TecUtilLockStart prior to coming */
    /* up and then call TecUtilLockFinish when the Ok or Cancel */
    /* button is pressed.  Only TecUtilLockFinish is supplied here. */
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
}


/**
*/
static LgIndex_t  IMax_TF_D2_CB(char const *S)
{
    int IMax;
    LgIndex_t IsOk = 1;
    TecUtilLockStart(AddOnID);
    TRACE1("Text field (IMax_TF_D2) Value Changed,  New value is: %s\n", S);
    if (TecGUITextFieldGetLgIndex(IMax_TF_D2, &IMax) != TRUE)
    {
        TecUtilDialogErrMsg("Invalid number.");
        TecGUITextFieldSetLgIndex(IMax_TF_D2, 10, TRUE);
    }

    if ((IMax <= 0) || (IMax >= 99))
    {
        TecUtilDialogErrMsg("Small positive integer required.");
        TecGUITextFieldSetLgIndex(IMax_TF_D2, 10, TRUE);
    }
    TecUtilLockFinish(AddOnID);
    return (IsOk);
}



/**
*/
static LgIndex_t  JMax_TF_D2_CB(char const *S)
{
    int JMax;
    LgIndex_t IsOk = 1;
    TecUtilLockStart(AddOnID);
    TRACE1("Text field (JMax_TF_D2) Value Changed,  New value is: %s\n", S);
    if (TecGUITextFieldGetLgIndex(JMax_TF_D2, &JMax) != TRUE)
    {
        TecUtilDialogErrMsg("Small positive integer required.");
        TecGUITextFieldSetLgIndex(JMax_TF_D2, 10, TRUE);
    }

    if ((JMax <= 0) || (JMax >= 99))
    {
        TecUtilDialogErrMsg("Invalid number.");
        TecGUITextFieldSetLgIndex(JMax_TF_D2, 10, TRUE);
    }
    TecUtilLockFinish(AddOnID);
    return (IsOk);
}





#include "guibld.cpp"


