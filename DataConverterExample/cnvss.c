/*

  Converter Example Add-on:

  Loads a text file in the following form:

  Var 1  Var 2  ...  Var N
  data1  data2  ...  data N
  .
  .
  .
  data1  data2  ...  data N

  [Blank lines are ignored]

  Notes:

  1. There are no restrictions on the position
     of the data, except that the variable names
     must be in quotes and must be first.

  2. Data points must be separated with spaces.

*/


#include "TECADDON.h"
AddOn_pa AddOnID;

#include "ADDONVER.h"
#define MAX_VARS  100     /* max number of variables allowed */

#ifndef ASSERT
#include <assert.h>
#define ASSERT(exp) assert(exp)
#endif




/*
*
* Exported
*
*/
Boolean_t STDCALL ConvSS(char *OldName, char *TempName, char **MessageString);

/*
* private to
* this file
*/

static void STDCALL StateChangeCallback(StateChange_e WhichState,
                                        ArbParam_t NotUsed);

static Boolean_t DoConvSS(FILE *f, char *TempName, char **MessageString);


/*------------------------------------------------------------
*
* get_token
*
*------------------------------------------------------------*/
#define MAX_TOKEN_LEN 5000
static char _token[MAX_TOKEN_LEN]; /* global buffer for tokens */

static Boolean_t get_token(FILE *f)
{
    /* returns FALSE if no more tokens */

    int index = 0;
    char c;
    Boolean_t StopRightQuote;

    ASSERT(f);
    /*
     * Note that f is assumed to
     * have been opened in binary
     * mode.
     */

    /*
     * skip whitespace
     */
    while (fread(&c, sizeof(char), 1, f) == 1
           && (c == ' ' || c == ',' || c == '\t' || c == '\n' || c == '\r'))
    {
        /* keep going */
    }

    if (!feof(f))
    {
        /*
         * now we're sitting on a non-whitespace character
         */

        StopRightQuote = (c == '"');
        if (StopRightQuote)
        {
            _token[index++] = c;
            fread(&c, sizeof(char), 1, f);
        }



        do
        {
            if (index == MAX_TOKEN_LEN - 1)
                break; /* ouch, we really shouldn't have lines > 5000 */

            if (feof(f))
                break;

            if (StopRightQuote)
            {
                if (c == '"')
                {
                    _token[index++] = c;
                    break;
                }
            }
            else
            {
                /* note that a space or comma may terminate the token */
                if (c == ',' || c == ' ' || c == '\t' || c == '\n' || c == '\r')
                    break;
            }

            _token[index++] = c;
            fread(&c, sizeof(char), 1, f);
        }
        while (1);
    }

    _token[index] = '\0';

    return (strlen(_token) > 0);
}


/*
* Tecplot will explicitly call
* InitTecAddOn() when the add-on
* is first loaded. The only thing
* we *have* to do is call
* TecUtilAddOnRegisterInfo().
*/


EXPORTFROMADDON void STDCALL InitTecAddOn(void)
{
    /*
     * Remember: Tecplot must be locked
     * before calling any TecUtil functions...
     */
    TecUtilLockOn();
    AddOnID = TecUtilAddOnRegister(100, ADDON_NAME, ADDON_VERSION" ("TecVersionId") "ADDON_DATE,
                                   "Tecplot, Inc.");

    if (TecUtilGetTecplotVersion() < MinTecplotVersionAllowed)
    {
        char buffer[256];
        sprintf(buffer, "Add-on \"%s\" requires Tecplot version %s or greater", ADDON_NAME, TecVersionId);
        TecUtilDialogErrMsg(buffer);
    }

    else
    {

        TecUtilImportAddConverter(ConvSS, ADDON_NAME, "*.txt");
    }

    /*
     * Note that the number of
     * TecUtilLockFinish(AddOnID);  * TecUtilLockOff()'s must
     * equal the number of
     * TecUtilLockStart(AddOnID);  * TecUtilLockOn()'s
     */

    TecUtilLockOff();
}


Boolean_t STDCALL ConvSS(char *OldName, char *TempName, char **MessageString)
{
    Boolean_t IsOk = TRUE;
    FILE *f;

    TecUtilLockStart(AddOnID);
    /*
     * We'll free this ourselves
     * if there's no error
     */

    *MessageString = TecUtilStringAlloc(1000, "MessageString for CNVSS");
    strcpy(*MessageString, "Error reading data set");

    /*
     * Try to open the file...
     */

    f = fopen(OldName, "rb");

    if (!f)
    {
        strcpy(*MessageString, "Cannot open input file");
        IsOk = FALSE;
    }

    if (IsOk)
        IsOk = DoConvSS(f, TempName, MessageString);

    fclose(f);

    if (IsOk)
        TecUtilStringDealloc(MessageString);

    TecUtilLockFinish(AddOnID);
    return IsOk;
}

static fpos_t _DataStartPos;

static void GetVars(FILE* f, StringList_pa sl)
{
    char  c;
    char  buffer[5000];
    char  *Line = buffer;
    char  Var[100];
    int   Index = 0;
    char Delimiter = ' ';

    /* read up to the first newline */

    do
    {
        if (fread(&c, sizeof(char), 1, f) < 1)
            break;

        if (c != '\r' && c != '\n' && c != '\0')
            buffer[Index++] = c;
        else
            break;
    }
    while (1);

    buffer[Index] = '\0';
    /* now get the variable names */

    Index = 0;

    while (*Line)
    {
        Index = 0;
        if (*Line == '"')
        {
            /* skip to next double quote */
            Line++;
            while (*Line && *Line != '"')
                Var[Index++] = *Line++;
        }

        else
        {
            /* just read to the next delimiter */
            while (*Line && *Line != Delimiter)
                Var[Index++] = *Line++;
        }

        Var[Index] = '\0';
        TecUtilStringListAppendString(sl, Var);

        /* now skip to the next non-delimiter char */
        while (*Line && *Line != Delimiter)
            Line++;

        fgetpos(f, &_DataStartPos);

        /* skip to next non-delimiter char */
        while (*Line && (*Line == Delimiter || *Line == ' '))
            Line++;
    }
}


static Boolean_t DoConvSS(FILE *f, char *TempName, char **MessageString)
{
    Boolean_t IsOk = TRUE;
    StringList_pa sl_var = TecUtilStringListAlloc(); /* variable list */
    long FileSize;
    int i;
    int NumValues;
    int NumVars;
    int imax;
    Boolean_t PercentDoneLaunched = FALSE;


    /*
     * compute the file size
     */

    fseek(f, 0, SEEK_END);
    FileSize = ftell(f);
    rewind(f);

    /*
     * First,
     * we need to read
     * all of the variables,
     */


    GetVars(f, sl_var);

    /*
     * now we have all of the variables.
     * There must be at least one
     */

    if (IsOk && TecUtilStringListGetCount(sl_var) < 1)
    {
        strcpy(*MessageString, "No variables defined");
        IsOk = FALSE;
    }

    /*
     * count the number of data points
     */

    if (IsOk)
    {
        int Debug = 0;
        int VIsDouble = 1;
        char var_names[5000];
        NumValues = 0;
        NumVars = TecUtilStringListGetCount(sl_var);

        TecUtilDialogLaunchPercentDone("Scanning File...", TRUE);
        PercentDoneLaunched = TRUE;

        while (get_token(f))
        {
            TecUtilDialogCheckPercentDone((int)(100.0*ftell(f) / FileSize));
            NumValues++;
        }
        TecUtilDialogDropPercentDone();
        PercentDoneLaunched = FALSE;


        fsetpos(f, &_DataStartPos); /* rewind to the start of the data */

        /*
         * compute the number of data points
         */

        if (NumValues < 1)
        {
            strcpy(*MessageString, "Invalid or corrupt data file.");
            IsOk = FALSE;
        }

        if (IsOk)
        {

            int jmax = 1, kmax = 1;
            char *s;

            imax = NumValues / NumVars;
            strcpy(var_names, "");
            for (i = 1; i <= NumVars && IsOk; i++)
            {
                s = TecUtilStringListGetString(sl_var, i);
                strcat(var_names, s);
                if (i < NumVars)
                    strcat(var_names, ",");
                TecUtilStringDealloc(&s);
            }

            if (TecUtilTecIni("Converted Dataset", var_names, TempName, ".", &Debug, &VIsDouble) != 0)
            {
                strcpy(*MessageString, "Could not create datset");
                IsOk = FALSE;
            }

            if (IsOk && TecUtilTecZne("Zone 1", &imax, &jmax, &kmax, "POINT", NULL) != 0)
            {
                strcpy(*MessageString, "Could not add zone");
                IsOk = FALSE;
            }

        }

        /* now add the data */

        if (IsOk)
        {
            double *LineValues = (double*) calloc(NumValues, sizeof(double));

            TecUtilDialogLaunchPercentDone("Adding Data...", TRUE);
            PercentDoneLaunched = TRUE;

            for (i = 0; i < NumValues; i++)
            {
                get_token(f);
                LineValues[i] = atof(_token);
            }

            if (TecUtilTecDat(&NumValues, (void*)LineValues, &VIsDouble) != 0)
            {
                strcpy(*MessageString, "Error loading data");
                IsOk = FALSE;
            }

            free(LineValues);
        }
    }

    if (PercentDoneLaunched)
        TecUtilDialogDropPercentDone();

    if (IsOk && TecUtilTecEnd() != 0)
    {
        IsOk = FALSE;
        strcpy(*MessageString, "Invalid or corrupt data file");
    }

    TecUtilStringListDealloc(&sl_var);
    return IsOk;
}
