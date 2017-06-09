#include "LoadData.h"
#include "../Lock.h"

std::string CreateTestFile(LgIndex_t IMax,
                           LgIndex_t JMax,
                           double    XStart,
                           double    YStart,
                           double    XDelta,
                           double    YDelta)
{
    char   IMaxString[4];
    char   JMaxString[4];
    char   XStartString[8];
    char   YStartString[8];
    char   XDeltaString[8];
    char   YDeltaString[8];

    sprintf(IMaxString, "%d", IMax);
    sprintf(JMaxString, "%d", JMax);
    sprintf(XStartString, "%.2f", XStart);
    sprintf(YStartString, "%.2f", YStart);
    sprintf(XDeltaString, "%.2f", XDelta);
    sprintf(YDeltaString, "%.2f", YDelta);

    std::string FileNameString = "X";
    FileNameString += IMaxString;
    FileNameString += "_";
    FileNameString += XStartString;
    FileNameString += "_";
    FileNameString += XDeltaString;
    FileNameString += "_";
    FileNameString += "Y";
    FileNameString += JMaxString;
    FileNameString += "_";
    FileNameString += YStartString;
    FileNameString += "_";
    FileNameString += YDeltaString;
    FileNameString += ".bin";

    FILE* MyFile = fopen(FileNameString.c_str(), "wb");
    bool isOk = (MyFile != NULL);
    if (isOk)
    {
        isOk = (fwrite(&IMax, sizeof(int), 1, MyFile) == 1 &&
                fwrite(&JMax, sizeof(int), 1, MyFile) == 1);

        double XVarValue = XStart;
        for (LgIndex_t PointIndex = 0; PointIndex < IMax && isOk; ++PointIndex)
        {
            isOk = (fwrite(&XVarValue, sizeof(double), 1, MyFile) == 1);
            XVarValue += XDelta;
        }

        double YVarValue = YStart;
        for (LgIndex_t PointIndex = 0; PointIndex < JMax && isOk; ++PointIndex)
        {
            isOk = (fwrite(&YVarValue, sizeof(double), 1, MyFile) == 1);
            YVarValue += YDelta;
        }

        fclose(MyFile);
    }

    if (!isOk)
    {
        FileNameString.clear();
    }

    return FileNameString;
}

static Boolean_t STDCALL LoadOnDemandVarLoad(FieldData_pa FieldData)
{
    REQUIRE(FieldData != NULL);
    REQUIRE(VALID_REF((ClientDataValues_s *)TecUtilDataValueGetClientData(FieldData)));

    Lock lockStart;

    ClientDataValues_s      *MyClientData = (ClientDataValues_s *)TecUtilDataValueGetClientData(FieldData);
    int                      NumValues = MyClientData->NumValues;
    MyClientData->StagingData = (double *) malloc(sizeof(double) * NumValues);

    Boolean_t IsOk = (MyClientData->StagingData != NULL);
    if (IsOk)
    {
        FILE *MyFile = fopen(MyClientData->FileName.c_str(), "rb");
        IsOk = (MyFile != NULL);
        if (IsOk)
        {
            /* Go to the position in the file for the first value of the variable. */
            IsOk = (fseek(MyFile, MyClientData->SeekOffset, SEEK_SET) == 0);
            if (IsOk)
            {
                LgIndex_t PointIndex;
                for (PointIndex = 0; PointIndex < NumValues; PointIndex++)
                {
                    double Value;
                    int NumValuesRead = (int)fread(&Value, sizeof(double), 1, MyFile);
                    if (NumValuesRead == 1)
                    {
                        MyClientData->StagingData[PointIndex] = Value;
                    }
                    else
                    {
                        IsOk = FALSE;
                    }
                }
            }
            fclose(MyFile);
        }
    }

    return IsOk;
}

/**
 */
static Boolean_t STDCALL LoadOnDemandVarUnload(FieldData_pa FieldData)
{
    Lock lockStart;

    REQUIRE(FieldData != NULL);
    REQUIRE(VALID_REF((ClientDataValues_s *)TecUtilDataValueGetClientData(FieldData)));

    ClientDataValues_s *MyClientData = (ClientDataValues_s *)TecUtilDataValueGetClientData(FieldData);
    free(MyClientData->StagingData);
    MyClientData->StagingData = NULL; /* set to NULL to prevent potential of multiple free attempts from within LoadOnDemandVarCleanup */

    return TRUE;
}

/**
*/
static void STDCALL LoadOnDemandVarCleanup(FieldData_pa FieldData)
{
    Lock lockStart;

    REQUIRE(FieldData != NULL);
    REQUIRE(VALID_REF((ClientDataValues_s *)TecUtilDataValueGetClientData(FieldData)));

    ClientDataValues_s *MyClientData = (ClientDataValues_s *)TecUtilDataValueGetClientData(FieldData);
    free(MyClientData->StagingData);
    MyClientData->StagingData = NULL;
    delete MyClientData;
    MyClientData = NULL;
}


static double STDCALL GetFirstStagedData(const FieldData_pa FieldData, LgIndex_t  PointIndex)
{
    Lock lockStart;

    REQUIRE(FieldData != NULL);
    REQUIRE(VALID_REF((ClientDataValues_s *)TecUtilDataValueGetClientData(FieldData)));

    ClientDataValues_s *MyClientData  = (ClientDataValues_s *)TecUtilDataValueGetClientData(FieldData);
    double             *StagingData   = MyClientData->StagingData;

    /* Set x values. */
    int    StageIndex = PointIndex % MyClientData->IMax;
    double Value      = StagingData[StageIndex];

    return Value;
}



static double STDCALL GetSecondStagedData(const FieldData_pa FieldData, LgIndex_t  PointIndex)
{
    Lock lockStart;

    REQUIRE(FieldData != NULL);
    REQUIRE(VALID_REF((ClientDataValues_s *)TecUtilDataValueGetClientData(FieldData)));

    ClientDataValues_s      *MyClientData = (ClientDataValues_s *)TecUtilDataValueGetClientData(FieldData);
    double                  *StagingData = MyClientData->StagingData;
    double                   D_IMax = (double) MyClientData->IMax;

    /* Set y values. */
    double                   Quotient = (double)PointIndex / D_IMax;
    int                      StageIndex = (int) Quotient;
    double                   Value = StagingData[StageIndex];
    return Value;
}

/* This loader requires the data file to fit the following specifications:
 * Header: IMax and JMax
 * Data: The data values are double and are contiguous for each variable
 * (not "point" format.)
 */
Boolean_t STDCALL LoaderCallback(StringList_pa Instructions) /* IN */
{
    char    *Name = NULL;
    char    *ValueString = NULL;
    char     FileName[120];
    int      IMax;
    int      JMax;
    int      InstrIndex;
    FILE    *MyFile;

    Boolean_t IsOk = FALSE;
    LgIndex_t Count;
    REQUIRE(Instructions != NULL);

    Lock lockStart;
    /* Check that Instructions are in STANDARDSYNTAX.
    */
    Count = TecUtilStringListGetCount(Instructions);
    Name = TecUtilStringListGetString(Instructions, 1);
    if (strcmp(Name, "STANDARDSYNTAX") == 0)
    {
        IsOk = TRUE;
    }
    TecUtilStringDealloc(&Name);

    /* Get pairs of instructions and convert the value string to the form
     * needed to load the data.
     */
    InstrIndex = 3;
    while (IsOk && InstrIndex < Count)
    {
        Name  = TecUtilStringListGetString(Instructions, InstrIndex);
        ValueString = TecUtilStringListGetString(Instructions, InstrIndex + 1);
        InstrIndex += 2;

        if (strcmp(Name, "FILENAME") == 0)
        {
            strcpy(FileName, ValueString);
            REQUIRE(VALID_NON_ZERO_LEN_STR(FileName));
        }
        /*   else if (Str_ustrcmp(Name, "Function") == 0)
             {
             Function = atoi(ValueString);
             }
         * First time through, just making a grid, no function options.
         */
        else
        {
            break;
        }

        TecUtilStringDealloc(&Name);
        TecUtilStringDealloc(&ValueString);
    }


    /* Get header info : IMax, JMax */
    MyFile = fopen(FileName, "rb");
    if (MyFile != NULL)
    {
        IsOk = (fread(&IMax, sizeof(int), 1, MyFile) == 1 &&
                fread(&JMax, sizeof(int), 1, MyFile) == 1);
        fclose(MyFile);
    }
    else
        IsOk = FALSE;

    if (IsOk == TRUE)
    {
        EntIndex_t       VarIndex;
        FieldDataType_e *VarDataTypes;
        ArgList_pa       ArgList;
        LgIndex_t        SeekOffset = 2 * sizeof(int);
        StringList_pa    VarNames = TecUtilStringListAlloc();

        VarDataTypes = (FieldDataType_e *)TecUtilStringAlloc(sizeof(FieldDataType_Double) * 2,
                                                             "Var Data Types");
        for (VarIndex = 0; VarIndex < 2; VarIndex++)
            VarDataTypes[VarIndex] = FieldDataType_Double;

        /* For the first test, just one zone, two
         * variables in calculated data set.
         */
        TecUtilStringListAppendString(VarNames, "X"); // first variable
        TecUtilStringListAppendString(VarNames, "Y"); // second variable

        if (TecUtilDataSetCreate("Calculated Data Set", VarNames, TRUE))
        {

            ArgList = TecUtilArgListAlloc();
            TecUtilArgListAppendString(ArgList, SV_NAME,                "Calculated Zone");
            TecUtilArgListAppendInt(ArgList, SV_IMAX,                IMax);
            TecUtilArgListAppendInt(ArgList, SV_JMAX,                JMax);
            TecUtilArgListAppendInt(ArgList, SV_DEFERVARCREATION,    TRUE);
            TecUtilArgListAppendArray(ArgList, SV_VARDATATYPE, (void *)VarDataTypes);
            IsOk = TecUtilDataSetAddZoneX(ArgList);
            if (VarDataTypes)
                TecUtilStringDealloc((char**)(void*)&VarDataTypes);

            TecUtilArgListDealloc(&ArgList);
            if (IsOk)
            {
                for (VarIndex = 1; VarIndex <= 2; VarIndex++)
                {
                    ClientDataValues_s *ClientData = new ClientDataValues_s;
                    ClientData->FileName  = FileName;
                    ClientData->IMax      = IMax;
                    ClientData->JMax      = JMax;
                    if (VarIndex == 1)
                    {
                        ClientData->NumValues = IMax;
                        ClientData->SeekOffset = SeekOffset;
                        ClientData->StagingData = NULL;
                        IsOk = TecUtilDataValueCustomLOD((EntIndex_t)1,
                                                         VarIndex,
                                                         LoadOnDemandVarLoad,
                                                         LoadOnDemandVarUnload,
                                                         LoadOnDemandVarCleanup,
                                                         GetFirstStagedData,
                                                         NULL,
                                                         (ArbParam_t)ClientData);
                    }
                    else if (VarIndex == 2)
                    {
                        ClientData->NumValues = JMax;
                        ClientData->SeekOffset = SeekOffset + IMax * (sizeof(double));
                        ClientData->StagingData = NULL; /*(double *) malloc(sizeof(double) * JMax); */
                        IsOk = TecUtilDataValueCustomLOD((EntIndex_t)1,
                                                         VarIndex,
                                                         LoadOnDemandVarLoad,
                                                         LoadOnDemandVarUnload,
                                                         LoadOnDemandVarCleanup,
                                                         GetSecondStagedData,
                                                         NULL,
                                                         (ArbParam_t)ClientData);
                    }
                }
            }
        }
        TecUtilDataSetDefVarLoadFinish(IsOk);
        TecUtilFrameSetPlotType(PlotType_Cartesian2D);
        TecUtilRedraw(TRUE);
        TecUtilImportSetLoaderInstr(LoaderName.c_str(), Instructions);
        TecUtilStringListDealloc(&VarNames);
    }

    return (IsOk);
}
