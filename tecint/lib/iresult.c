
#if !defined ADDON
#define ADDON
#endif /* ADDON */
#include "TECADDON.h"

#include "tecint.h"
#include "iresult.h"

typedef struct SingleResult_s
{
    EntIndex_t Zone;
    LgIndex_t  I;
    LgIndex_t  J;
    LgIndex_t  K;
    double     Value;
} SingleResult_t;

typedef struct IntegrationResult_s
{
    Set_pa           ZoneSet;
    IntegrateOver_e  IntegrateOver;
    int              ResultCount;
    int              MaxResultCount;
    int              NextResultNum;
    SingleResult_t  *ResultList;
} IntegrationResult_t;

Boolean_t IntegrationResultAlloc(IntegrationResult_pa *IntegrationResult)
{
    Boolean_t IsOk = TRUE;

    if (!IntegrationResult)
    {
        TecUtilDialogErrMsg("Bad argument passed to IntegrationResultAlloc().");
        IsOk = FALSE;
    }

    if (IsOk)
    {
        *IntegrationResult = malloc(sizeof(IntegrationResult_t));
        if (!(*IntegrationResult))
        {
            TecUtilDialogErrMsg("Out of memory while allocating integration result.");
            IsOk = FALSE;
        }
    }

    if (IsOk)
    {
        (*IntegrationResult)->ZoneSet = TecUtilSetAlloc(TRUE);
        if (!(*IntegrationResult)->ZoneSet)
        {
            free(*IntegrationResult);
            *IntegrationResult = NULL;
            IsOk = FALSE;
        }
    }

    if (IsOk)
    {
        (*IntegrationResult)->IntegrateOver  = IntegrateOver_Volume;
        (*IntegrationResult)->ResultList     = NULL;
        (*IntegrationResult)->ResultCount    = 0;
        (*IntegrationResult)->MaxResultCount = 0;
        (*IntegrationResult)->NextResultNum  = 0;
    }

    return(IsOk);
}

void IntegrationResultDealloc(IntegrationResult_pa *IntegrationResult)
{
    TecUtilSetDealloc(&(*IntegrationResult)->ZoneSet);
    free((*IntegrationResult)->ResultList);
    free(*IntegrationResult);
    *IntegrationResult = NULL;
}

int IntegrationResultGetCount(IntegrationResult_pa IntegrationResult)
{
    return(IntegrationResult->ResultCount);
}

/*
 * Assign irrelevant indices to unity to enable bsearching.
 */
static void AssignSingleResultIndices(SingleResult_t *SingleResult,
                                      IntegrateOver_e IntegrateOver,
                                      EntIndex_t      Zone,
                                      LgIndex_t       I,
                                      LgIndex_t       J,
                                      LgIndex_t       K)
{
    SingleResult->Zone = Zone;
    SingleResult->I = 1;
    SingleResult->J = 1;
    SingleResult->K = 1;

    switch (IntegrateOver)
    {
        case IntegrateOver_Volume:
            break;
        case IntegrateOver_IPlanes:
            SingleResult->I = I;
            break;
        case IntegrateOver_JPlanes:
            SingleResult->J = J;
            break;
        case IntegrateOver_KPlanes:
            SingleResult->K = K;
            break;
        case IntegrateOver_ILines:
            SingleResult->J = J;
            SingleResult->K = K;
            break;
        case IntegrateOver_JLines:
            SingleResult->I = I;
            SingleResult->K = K;
            break;
        case IntegrateOver_KLines:
            SingleResult->I = I;
            SingleResult->J = J;
            break;
    }

    if (I == TECINT_TOTAL)
    {
        SingleResult->I = I;
        SingleResult->J = 1;
        SingleResult->K = 1;
    }

    if (Zone == TECINT_TOTAL)
    {
        SingleResult->I = 1;
        SingleResult->J = 1;
        SingleResult->K = 1;
    }
}

static int CompareSingleResults(const void *Ptr1, const void *Ptr2)
{
    SingleResult_t *Result1 = (SingleResult_t *)Ptr1;
    SingleResult_t *Result2 = (SingleResult_t *)Ptr2;

    if (Result1->Zone > Result2->Zone)
    {
        return(1);
    }
    else if (Result1->Zone < Result2->Zone)
    {
        return(-1);
    }
    else if (Result1->K > Result2->K)
    {
        return(1);
    }
    else if (Result1->K < Result2->K)
    {
        return(-1);
    }
    else if (Result1->J > Result2->J)
    {
        return(1);
    }
    else if (Result1->J < Result2->J)
    {
        return(-1);
    }
    else if (Result1->I > Result2->I)
    {
        return(1);
    }
    else if (Result1->I < Result2->I)
    {
        return(-1);
    }
    else
    {
        return(0);
    }
}

Boolean_t IntegrationResultGetValue(IntegrationResult_pa  IntegrationResult,
                                    EntIndex_t            Zone,
                                    LgIndex_t             I,
                                    LgIndex_t             J,
                                    LgIndex_t             K,
                                    double               *Value)
{
    Boolean_t IsOk = TRUE;

    if (IntegrationResult->ResultCount == 0)
    {
        IsOk = FALSE;
    }
    else if (!TecUtilSetIsMember(IntegrationResult->ZoneSet, Zone))
    {
        IsOk = FALSE;
    }
    else
    {
        SingleResult_t  Key;
        SingleResult_t *SingleResult;

        AssignSingleResultIndices(&Key, IntegrationResult->IntegrateOver, Zone, I, J, K);
        SingleResult = bsearch(&Key,
                               IntegrationResult->ResultList,
                               (size_t)IntegrationResult->ResultCount,
                               sizeof(SingleResult_t),
                               CompareSingleResults);
        if (!SingleResult)
        {
            IsOk = FALSE;
        }
        else
        {
            *Value = SingleResult->Value;
        }
    }

    return(IsOk);
}

Boolean_t IntegrationResultGetFirstValue(IntegrationResult_pa  IntegrationResult,
                                         EntIndex_t           *Zone,
                                         LgIndex_t            *I,
                                         LgIndex_t            *J,
                                         LgIndex_t            *K,
                                         double               *Value)
{
    Boolean_t IsOk = TRUE;

    IntegrationResult->NextResultNum = 0;
    IsOk = IntegrationResultGetNextValue(IntegrationResult,
                                         Zone,
                                         I,
                                         J,
                                         K,
                                         Value);

    return(IsOk);
}

Boolean_t IntegrationResultGetNextValue(IntegrationResult_pa  IntegrationResult,
                                        EntIndex_t           *Zone,
                                        LgIndex_t            *I,
                                        LgIndex_t            *J,
                                        LgIndex_t            *K,
                                        double               *Value)
{
    Boolean_t IsOk = TRUE;

    if (IntegrationResult->NextResultNum == IntegrationResult->ResultCount)
    {
        IsOk = FALSE;
    }
    else
    {
        SingleResult_t *SingleResult = &IntegrationResult->ResultList[IntegrationResult->NextResultNum];
        *Zone  = SingleResult->Zone;
        *I     = SingleResult->I;
        *J     = SingleResult->J;
        *K     = SingleResult->K;
        *Value = SingleResult->Value;
        IntegrationResult->NextResultNum++;
    }

    return(IsOk);
}

Boolean_t IntegrationResultGetNumberedValue(IntegrationResult_pa  IntegrationResult,
                                            int                   WhichValue,
                                            EntIndex_t           *Zone,
                                            LgIndex_t            *I,
                                            LgIndex_t            *J,
                                            LgIndex_t            *K,
                                            double               *Value)
{
    Boolean_t IsOk = TRUE;

    if (WhichValue <= 0 || IntegrationResult->ResultCount < WhichValue)
    {
        IsOk = FALSE;
    }
    else
    {
        SingleResult_t *SingleResult = &IntegrationResult->ResultList[WhichValue - 1];
        *Zone  = SingleResult->Zone;
        *I     = SingleResult->I;
        *J     = SingleResult->J;
        *K     = SingleResult->K;
        *Value = SingleResult->Value;
        IntegrationResult->NextResultNum = WhichValue;
    }

    return(IsOk);
}

void IntegrationResultInitialize(IntegrationResult_pa IntegrationResult,
                                 IntegrateOver_e      IntegrateOver)
{
    TecUtilSetClear(IntegrationResult->ZoneSet);
    IntegrationResult->IntegrateOver = IntegrateOver;
    if (IntegrationResult->MaxResultCount > 0)
    {
        free(IntegrationResult->ResultList);
    }
    IntegrationResult->ResultList     = NULL;
    IntegrationResult->ResultCount    = 0;
    IntegrationResult->MaxResultCount = 0;
    IntegrationResult->NextResultNum  = 0;
}

extern Boolean_t IntegrationResultAddValue(IntegrationResult_pa IntegrationResult,
                                           EntIndex_t           Zone,
                                           LgIndex_t            I,
                                           LgIndex_t            J,
                                           LgIndex_t            K,
                                           double               Value)
{
    Boolean_t IsOk = TRUE;

    if (Zone > 0)
    {
        IsOk = TecUtilSetAddMember(IntegrationResult->ZoneSet, Zone, TRUE);
    }

    if (IsOk)
    {
        if (IntegrationResult->MaxResultCount == 0)
        {
            IntegrationResult->ResultList = malloc(sizeof(SingleResult_t));
            if (!IntegrationResult->ResultList)
            {
                TecUtilDialogErrMsg("Out of memory while storing integration results.");
                IsOk = FALSE;
            }
            else
            {
                IntegrationResult->MaxResultCount = 1;
            }
        }
        else if (IntegrationResult->ResultCount == IntegrationResult->MaxResultCount)
        {
            SingleResult_t *NewResultList = realloc(IntegrationResult->ResultList,
                                                    2 * IntegrationResult->MaxResultCount * sizeof(SingleResult_t));
            if (!NewResultList)
            {
                TecUtilDialogErrMsg("Out of memory while storing integration results.");
                IsOk = FALSE;
            }
            else
            {
                IntegrationResult->ResultList = NewResultList;
                IntegrationResult->MaxResultCount *= 2;
            }
        }
    }

    if (IsOk)
    {
        SingleResult_t *SingleResult = &IntegrationResult->ResultList[IntegrationResult->ResultCount];
        AssignSingleResultIndices(SingleResult, IntegrationResult->IntegrateOver, Zone, I, J, K);
        SingleResult->Value = Value;
        IntegrationResult->ResultCount++;
    }

    return(IsOk);
}

