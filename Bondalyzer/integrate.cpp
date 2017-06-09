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
#include "tecint.h"
#include "INTEGRATE.h"

static Boolean_t integrateScalarOverVolumeZone(EntIndex_t   zoneNum,
                                               EntIndex_t   scalarVarNum,
                                               const char * auxDataName,
                                               double *     integrationValue)
{
    Boolean_t isOk;

    EntIndex_t xVarNum = TecUtilVarGetNumByAssignment('X');
    EntIndex_t yVarNum = TecUtilVarGetNumByAssignment('Y');
    EntIndex_t zVarNum = TecUtilVarGetNumByAssignment('Z');

    Set_pa zoneSet = TecUtilSetAlloc(TRUE);
    TecUtilSetAddMember(zoneSet, zoneNum, TRUE);

    IndexRange_t iRange;
    iRange.Min  = 1;
    iRange.Max  = 0;
    iRange.Skip = 1;
    IndexRange_t jRange;
    jRange.Min  = 1;
    jRange.Max  = 0;
    jRange.Skip = 1;
    IndexRange_t kRange;
    kRange.Min  = 1;
    kRange.Max  = 0;
    kRange.Skip = 1;
    IntegrationResult_pa integrationResult;
    IntegrationResultAlloc(&integrationResult);
    ReturnStatus_e returnStatus = Integrate(integrationResult,
                                            VariableOption_ScalarIntegral,
                                            FALSE/*axisymmetric*/,
                                            SymmetryVar_X/*symmetryVar:NotUsed*/,
                                            0.0/*symmetryValue:NotUsed*/,
                                            scalarVarNum,
                                            xVarNum,
                                            yVarNum,
                                            zVarNum,
                                            IntegrateOver_Volume,
                                            zoneSet,
                                            iRange,
                                            jRange,
                                            kRange,
                                            FALSE /*absolute*/,
                                            FALSE /*excludeBlanked*/);
    if (returnStatus != ReturnStatus_OK)
    {
        TecUtilDialogErrMsg("Integrate() returned error.");
        isOk = FALSE;
    }
    if (isOk && IntegrationResultGetCount(integrationResult) != 2)   // one for zoneNum and one for the total
    {
        TecUtilDialogErrMsg("Integrate() returned wrong number of values.");
        isOk = FALSE;
    }

    if (isOk)
    {
        if (!IntegrationResultGetValue(integrationResult, zoneNum, 1, 1, 1, integrationValue))
        {
            TecUtilDialogErrMsg("Cannot get result of Integrate().");
            isOk = FALSE;
        }
    }

    if (isOk)
    {
        char strValue[200];
        sprintf_s(strValue, "%g", *integrationValue);
        AuxData_pa auxDataRef = TecUtilAuxDataZoneGetRef(zoneNum);
        if (auxDataRef != NULL)
        {
            if (!TecUtilAuxDataSetStrItem(auxDataRef,
                                          auxDataName,
                                          strValue,
                                          TRUE))
            {
                TecUtilDialogErrMsg("Cannot set aux data for zoneNum.");
                isOk = FALSE;
            }

        }
        else
        {
            TecUtilDialogErrMsg("Cannot get zoneNum aux-data ref.");
            isOk = FALSE;
        }
    }
    IntegrationResultDealloc(&integrationResult);
    TecUtilSetDealloc(&zoneSet);
    return isOk;
}


/*
 * The keyVolumeZone will not be used in totals
 */
Boolean_t integrateScalarOverVolumeZones(EntIndex_t   scalarVarNum,
                                         char const * auxDataName,
                                         Set_pa       zoneSet)
{
    Boolean_t isOk = TRUE;
    Boolean_t integratedAtLeastOneZone = FALSE;
    EntIndex_t numZones;
    TecUtilDataSetGetInfo(NULL, &numZones, NULL);
    double integrationTotal = 0.0;
    EntIndex_t minZoneNumIntegrated = MAXINT32;
    EntIndex_t maxZoneNumIntegrated = 0;
    SetIndex_t zoneNum;
    TecUtilSetForEachMember(zoneNum, zoneSet)
    {
        ZoneType_e zoneType = TecUtilZoneGetType((EntIndex_t)zoneNum);
        Boolean_t doIntegration;
        if (zoneType == ZoneType_Ordered)
        {
            EntIndex_t iMax, jMax, kMax;
            TecUtilZoneGetInfo((EntIndex_t)zoneNum, &iMax, &jMax, &kMax, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);
            doIntegration = (iMax > 1 && jMax > 1 && kMax > 1);
        }
        else
        {
            doIntegration = (zoneType == ZoneType_FETetra || zoneType == ZoneType_FEBrick);
            if (zoneType == ZoneType_FEPolyhedron)
            {
                TecUtilDialogErrMsg("Polyhedral zones not available for integration");
                isOk = FALSE;
            }
        }
        if (doIntegration)
        {
            double integrationValue;
            isOk = integrateScalarOverVolumeZone((EntIndex_t)zoneNum, scalarVarNum, auxDataName, &integrationValue);
            if (isOk)
                integrationTotal += integrationValue;
            integratedAtLeastOneZone = TRUE;
            maxZoneNumIntegrated = MAX(maxZoneNumIntegrated, (EntIndex_t)zoneNum);
            minZoneNumIntegrated = MIN(minZoneNumIntegrated, (EntIndex_t)zoneNum);
        }
    }

    if (isOk)
    {
        if (maxZoneNumIntegrated == 0)
            TecUtilDialogMessageBox("No volume zones to integrate.", MessageBox_Warning);
        else
        {
            // add total to data set aux data
            char strValue[200];
            sprintf_s(strValue, "%g", integrationTotal);
            char dataSetAuxDataName[200];

            if (minZoneNumIntegrated == maxZoneNumIntegrated)
                sprintf_s(dataSetAuxDataName, "%sZone%d", auxDataName, minZoneNumIntegrated);
            else
                sprintf_s(dataSetAuxDataName, "%sZones%d_%d", auxDataName, minZoneNumIntegrated, maxZoneNumIntegrated);

            AuxData_pa auxDataRef = TecUtilAuxDataDataSetGetRef();
            if (auxDataRef == NULL)
            {
                TecUtilDialogErrMsg("Cannot get data-set aux-data ref.");
                isOk = FALSE;
            }
            else if (!TecUtilAuxDataSetStrItem(auxDataRef,
                                               dataSetAuxDataName,
                                               strValue,
                                               FALSE))
            {
                TecUtilDialogErrMsg("Cannot set aux data for data set.");
                isOk = FALSE;
            }
        }
    }

    return isOk;
}
