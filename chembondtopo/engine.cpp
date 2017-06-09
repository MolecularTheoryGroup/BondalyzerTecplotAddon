#include "TECADDON.h"
#include "ADDGLBL.h"
#include "ADDONVER.h"
#include "TASSERT.h"
#include "GUIDEFS.h"
#include "ADKUTIL.h"
#include "ALLOC.h"
#include "ARRLIST.h"
#include "CRITPOINTS.h"
#include "BUNDLES.h"
#include "ENGINE.h"
#include <string.h>



Boolean_t RecordMacroAndJournal(Strand_t StrandToAve,
                                EntIndex_t VarToAve)
{
    Boolean_t IsOk = TRUE;

    char *     DataSetTitle;
    EntIndex_t NumZones, NumVars;

    IsOk = TecUtilDataSetGetInfo(&DataSetTitle, &NumZones, &NumVars);
    TecUtilStringDealloc(&DataSetTitle);

    REQUIRE(VarToAve    > 0 && VarToAve    <= NumVars);
    REQUIRE(StrandToAve > 0 && StrandToAve <= TecUtilDataSetGetMaxStrandID());

    if (IsOk)
    {
        char  S[1000];
        char *MacroString;

        sprintf(S,             "VarToAve = %d, ", VarToAve);
        sprintf(&S[strlen(S)], "StrandToAve = %d ", StrandToAve);

        S[strlen(S)] = '\0';
        MacroString = TecUtilStringAlloc((int)strlen(S), "MacroString");
        strcpy(MacroString, S);

        if (TecUtilDataSetJournalIsValid())
            IsOk = TecUtilDataSetAddJournalCommand(ADDON_NAME, MacroString, NULL);
        if (TecUtilMacroIsRecordingActive())
            IsOk = TecUtilMacroRecordAddOnCommand(ADDON_NAME, MacroString);
        TecUtilStringDealloc(&MacroString);
    }

    ENSURE(VALID_BOOLEAN(IsOk));
    return IsOk;
}








/**
 * Find and return the GradPath zone number with the specified beginning
 * and ending critical point numbers. Search throught the Tepclot zones.
 *
 * param BegCPNum
 *     Number of the beginning critcial point number for the GradPath.
 *
 * param EndCPNum
 *     Number of the ending critcial point number for the GradPath.
 *
 * return
 *     Zone number if it works, zero otherwise.
 */

EntIndex_t   LineZoneGetByBegEndCP(LgIndex_t BegCrtPtNum,
                                   LgIndex_t EndCrtPtNum)
{
    EntIndex_t   Result = 0;
    Boolean_t    IsFound = FALSE;
    Boolean_t    IsOk = TRUE;
    EntIndex_t   NumZones, NumVars;
    EntIndex_t   ChrgDensVarNum = 0;
    LgIndex_t    BegCrtPtNumTry, EndCrtPtNumTry;
    FieldData_pa XVarFDPtr = NULL;
    FieldData_pa YVarFDPtr = NULL;
    FieldData_pa ZVarFDPtr = NULL;
    FieldData_pa CDVarFDPtr = NULL;

    REQUIRE(TecUtilDataSetIsAvailable());

    /* Get num zones in dataset */
    if (IsOk)
    {
        char *DatasetTitle = NULL;
        TecUtilDataSetGetInfo(&DatasetTitle, &NumZones, &NumVars);
        TecUtilStringDealloc(&DatasetTitle);
    }

    /* Query the Auxiliary data for the zone to fill some of the structure */
    if (IsOk)
    {
        EntIndex_t ii;

        for (ii = 1; !IsFound && ii <= NumZones; ii++)
        {
            AuxData_pa AuxDataRef = TecUtilAuxDataZoneGetRef(ii);
            if (AuxDataRef != NULL)
            {
                char      *Value;
                Boolean_t  Retain;

                /* Get BegCrtPtNum */
                if (TecUtilAuxDataGetStrItemByName(AuxDataRef, "CompChem.BegCrtPtNum",
                                                   &Value, &Retain))
                {
                    if (Value != NULL)
                    {
                        /* Parse the value */
                        if (sscanf(Value, "%d", &BegCrtPtNumTry) == 1)
                        {
                            BegCrtPtNumTry--;
                        }
                        else
                        {
                            IsOk = FALSE;
                        }
                        // release the allocated string copy
                        TecUtilStringDealloc(&Value);
                    }
                    else
                    {
                        IsOk = FALSE;
                    }
                }

                /* Get EndCrtPtNum */
                if (TecUtilAuxDataGetStrItemByName(AuxDataRef, "CompChem.EndCrtPtNum",
                                                   &Value, &Retain))
                {
                    if (Value != NULL)
                    {
                        /* Parse the value */
                        if (sscanf(Value, "%d", &EndCrtPtNumTry) == 1)
                        {
                            EndCrtPtNumTry--;
                        }
                        else
                        {
                            IsOk = FALSE;
                        }
                        // release the allocated string copy
                        TecUtilStringDealloc(&Value);
                    }
                    else
                    {
                        IsOk = FALSE;
                    }
                }
            }
            else  /* ZoneAuxDataRef is Null */
                IsOk = FALSE;

            if (IsOk && (BegCrtPtNum == BegCrtPtNumTry) &&
                (EndCrtPtNum == EndCrtPtNumTry))
            {
                Result = ii;
                IsFound = TRUE;
            }
        }
    }


    ENSURE(Result >= 0);
    return Result;
}








Boolean_t ExtractTopoInZonesAuxData(Bundles_pa *   BundlesPtr,
                                    CritPoints_pa *CritPointsPtr)
{
    Boolean_t  IsOk = TRUE;
    Boolean_t  IsStop = FALSE;
    Boolean_t  CrtPtZoneFound = FALSE;
    EntIndex_t CritPointZoneNum = 0;
    Boolean_t  ShowStatusBar = TRUE;
    EntIndex_t NumZones, NumVars;
    EntIndex_t ChrgDensVarNum = 4;  // Default for old files
    // EntIndex_t TypeVarNum = 8;   // Default for old files
    EntIndex_t UVarNum = 5;  // Default for old files
    EntIndex_t VVarNum = 6;  // Default for old files
    EntIndex_t WVarNum = 7;  // Default for old files
    EntIndex_t  ii;
    ArrList_pa Atom1ForBond = NULL;
    ArrList_pa Atom2ForBond = NULL;
    ArrList_pa Atom1ForBondZoneNum = NULL;
    ArrList_pa Atom2ForBondZoneNum = NULL;
    CritPoints_pa CritPoints = NULL;
    Bundles_pa    Bundles = NULL;

    EntIndex_t TypeVarNum = TECUTILSETNOTMEMBER;

    REQUIRE(TecUtilDataSetIsAvailable());
    REQUIRE(CritPoints == NULL || CritPointsIsValid(CritPoints));
    REQUIRE(Bundles == NULL || BundlesIsValid(Bundles));

    /* If CritPoints or Bundles exist, dealloc them */
    if (CritPoints != NULL) CritPointsDealloc(&CritPoints);
    if (Bundles != NULL) BundlesDealloc(&Bundles);

    /* Get num zones in dataset */
    if (IsOk)
    {
        char *DatasetTitle = NULL;
        TecUtilDataSetGetInfo(&DatasetTitle, &NumZones, &NumVars);
        TecUtilStringDealloc(&DatasetTitle);
    }

    REQUIRE(NumZones > 2);


    // Read ChrgDensVarNum, UVarNum, VVarNum, WVarNum, and TypeVarNum from aux data
    if (IsOk)
    {
        AuxData_pa DSAuxDataRef = TecUtilAuxDataDataSetGetRef();
        if (DSAuxDataRef != NULL)
        {
            char      *Value;
            Boolean_t  Retain;

            // ChrgDensVarNum
            if (TecUtilAuxDataGetStrItemByName(DSAuxDataRef, "CompChem.ChrgDensVarNum",
                                               &Value, &Retain))
            {
                if (Value != NULL)
                {
                    /* Parse the value */
                    ChrgDensVarNum = atoi(Value);  // Convert to zero-based

                    // release the allocated string copy
                    TecUtilStringDealloc(&Value);
                }
                else
                {
                    IsOk = FALSE;
                }
            }

            // UVarNum
            if (TecUtilAuxDataGetStrItemByName(DSAuxDataRef, "CompChem.GradXVarNum",
                                               &Value, &Retain))
            {
                if (Value != NULL)
                {
                    /* Parse the value */
                    UVarNum = atoi(Value);  // Convert to zero-based

                    // release the allocated string copy
                    TecUtilStringDealloc(&Value);
                }
                else
                {
                    IsOk = FALSE;
                }
            }

            // VVarNum
            if (TecUtilAuxDataGetStrItemByName(DSAuxDataRef, "CompChem.GradYVarNum",
                                               &Value, &Retain))
            {
                if (Value != NULL)
                {
                    /* Parse the value */
                    VVarNum = atoi(Value);  // Convert to zero-based

                    // release the allocated string copy
                    TecUtilStringDealloc(&Value);
                }
                else
                {
                    IsOk = FALSE;
                }
            }

            // WVarNum
            if (TecUtilAuxDataGetStrItemByName(DSAuxDataRef, "CompChem.GradZVarNum",
                                               &Value, &Retain))
            {
                if (Value != NULL)
                {
                    /* Parse the value */
                    WVarNum = atoi(Value);  // Convert to zero-based

                    // release the allocated string copy
                    TecUtilStringDealloc(&Value);
                }
                else
                {
                    IsOk = FALSE;
                }
            }

            // TypeVarNum
            if (TecUtilAuxDataGetStrItemByName(DSAuxDataRef, "CompChem.TypeVarNum",
                                               &Value, &Retain))
            {
                if (Value != NULL)
                {
                    /* Parse the value */
                    TypeVarNum = atoi(Value);  // Convert to zero-based

                    // release the allocated string copy
                    TecUtilStringDealloc(&Value);
                }
                else
                {
                    IsOk = FALSE;
                }
            }

        }
    }


    // Find the critical point zone, determine the number of atoms, bond
    // rings and cages.
    for (ii = 1; IsOk && !CrtPtZoneFound && ii <= NumZones; ii++)
    {
        /* Set Zone aux data to completely transfer necessary information to Tecplot */
        AuxData_pa ZnAuxDataRef = TecUtilAuxDataZoneGetRef(ii);
        if (ZnAuxDataRef != NULL)
        {
            char      *Value;
            Boolean_t  Retain;
            if (TecUtilAuxDataGetStrItemByName(ZnAuxDataRef, "CompChem.ZoneType",
                                               &Value, &Retain))
            {
                if (strcmp(Value, "CriticalPoints") == 0)
                {
                    CrtPtZoneFound = TRUE;
                    CritPointZoneNum = ii;
                }
                if (Value != NULL) TecUtilStringDealloc(&Value);
            }

            // If IsFound, get the number of critical points
            if (CrtPtZoneFound)
            {
                // CritPoints = CritPointsGetFromTPZone(CritPointZoneNum, 4, 8, 5, 6, 7);
                CritPoints = CritPointsGetFromTPZone(CritPointZoneNum, ChrgDensVarNum, TypeVarNum,
                                                     UVarNum, VVarNum, WVarNum);
                /* TODO: Automate the variable numbers
                CritPoints_pa CritPointsGetFromTPZone(EntIndex_t TPZoneNum,
                                                      EntIndex_t ChrgDensVarNum,
                                                      EntIndex_t TypeVarNum,
                                                      EntIndex_t UVarNum,
                                                      EntIndex_t VVarNum,
                                                      EntIndex_t WVarNum)
                                        */
            }
        }
    }


    // Build Bond-Atom list by searhing through the zone
    // aux data. We know how many bonds there are, and that there
    // should be two atoms per bond.
    for (ii = 1; ii <= NumZones; ii++)
    {
        // If zone is a Bond/Ring/Cage zone it has auxiliary
        // data items named
        //     "CompChem.BegCrtPtType" == "Bond"
        //     "CompChem.EndCrtPtType" == "Atom"
        AuxData_pa AuxDataRef = TecUtilAuxDataZoneGetRef(ii);
        if (AuxDataRef != NULL)
        {
            char *BegTypeString;
            ArbParam_t StrValue;
            AuxDataType_e Type;
            Boolean_t     Retain;
            if (TecUtilAuxDataGetItemByName(AuxDataRef, "CompChem.BegCrtPtType",
                                            &StrValue, &Type, &Retain))
            {
                BegTypeString = (char *)StrValue;
                if (strcmp(BegTypeString, "Bond") == 0)
                {
                    char *EndTypeString;

                    // A beginning type of "Bond" was found, look for end type of "Cage"
                    if (TecUtilAuxDataGetItemByName(AuxDataRef, "CompChem.EndCrtPtType",
                                                    &StrValue, &Type, &Retain))
                    {
                        EndTypeString = (char *)StrValue;
                        if (strcmp(EndTypeString, "Atom") == 0)
                        {
                            // Beg and End types of "Bond" and "Atom" were found
                            char         *Value;
                            LgIndex_t BondCrtPtNum, AtomCrtPtNum;

                            // Get the "Bond", "Ring", and "Cage" critical point numbers
                            // Bond number
                            if (TecUtilAuxDataGetStrItemByName(AuxDataRef, "CompChem.BegCrtPtNum",
                                                               &Value, &Retain))
                            {
                                if (Value != NULL)
                                {
                                    /* Parse the value */
                                    BondCrtPtNum = atoi(Value) - 1;  // Convert to zero-based

                                    // release the allocated string copy
                                    TecUtilStringDealloc(&Value);
                                }
                                else
                                {
                                    IsOk = FALSE;
                                }
                            }
                            // Atom number
                            if (IsOk && TecUtilAuxDataGetStrItemByName(AuxDataRef, "CompChem.EndCrtPtNum",
                                                                       &Value, &Retain))
                            {
                                if (Value != NULL)
                                {
                                    /* Parse the value */
                                    AtomCrtPtNum = atoi(Value) - 1;  // Convert to zero-based

                                    // release the allocated string copy
                                    TecUtilStringDealloc(&Value);
                                }
                                else
                                {
                                    IsOk = FALSE;
                                }
                            }
                            // Now build the arraylists of two Atoms that connect to the Bond
                            if (IsOk)
                            {
                                ArrListItem_u Item;
                                LgIndex_t BondNum = BondCrtPtNum - CritPointsGetBegOffset(CritPoints, (char)(-1));


                                // Is Atom1 set?
                                if (Atom1ForBond == NULL)
                                {
                                    int jj;
                                    Atom1ForBond = ArrListAlloc(CritPoints->NumCrtPtsM1, ArrListType_Long);
                                    Atom2ForBond = ArrListAlloc(CritPoints->NumCrtPtsM1, ArrListType_Long);
                                    // Initialize to -1 to identify unset values
                                    for (jj = 0; jj < CritPoints->NumCrtPtsM1; jj++)
                                    {
                                        Item.Long = -1;
                                        IsOk = ArrListSetItem(Atom1ForBond, jj, Item);
                                        IsOk = ArrListSetItem(Atom2ForBond, jj, Item);
                                    }
                                    CHECK(ArrListIsValid(Atom1ForBond));
                                    CHECK(ArrListIsValid(Atom2ForBond));
                                }

                                if (ArrListGetItem(Atom1ForBond, BondNum).Long == -1)
                                {
                                    Item.Long = AtomCrtPtNum;
                                    IsOk = ArrListSetItem(Atom1ForBond, BondNum, Item);
                                }
                                else if (ArrListGetItem(Atom2ForBond, BondNum).Long == -1)
                                {
                                    Item.Long = AtomCrtPtNum;
                                    IsOk = ArrListSetItem(Atom2ForBond, BondNum, Item);
                                }
                            }
                        }
                        // release the allocated string copy
                        if (EndTypeString != NULL) TecUtilStringDealloc(&EndTypeString);
                    }
                }
                // release the allocated string copy
                if (BegTypeString != NULL) TecUtilStringDealloc(&BegTypeString);
            }
        }
    }



    /* Search through zones looking for bond/ring/cage zones or
     * bond/far-field surface zones. Create two bundles for each
     * bond/ring/cage zone or bond/far-field surface zone using the
     * two atoms per bond list created earlier.
     */
    for (ii = 1; ii <= NumZones; ii++)
    {
        // If zone is a Bond/Ring/Cage zone it has auxiliary
        // data items named
        //     "CompChem.BegCrtPtType" == "Bond"
        //     "CompChem.EndCrtPtType" == "Cage"
        //     "CompChem.ThrdCrtPtType" == "Ring"
        AuxData_pa AuxDataRef = TecUtilAuxDataZoneGetRef(ii);
        if (AuxDataRef != NULL)
        {
            char *BegTypeString;
            ArbParam_t StrValue;
            AuxDataType_e Type;
            Boolean_t     Retain;
            if (TecUtilAuxDataGetItemByName(AuxDataRef, "CompChem.BegCrtPtType",
                                            &StrValue, &Type, &Retain))
            {
                BegTypeString = (char *)StrValue;
                if (strcmp(BegTypeString, "Bond") == 0)
                {
                    Boolean_t IsFarField = FALSE;
                    LgIndex_t BondCrtPtNum;
                    LgIndex_t CageCrtPtNum = -1;
                    LgIndex_t RingCrtPtNum = -1;
                    LgIndex_t Ring2CrtPtNum = -1;
                    char *EndTypeString = NULL;
                    char *IsFarFieldString = NULL;
                    char *ThrdTypeString = NULL;
                    char EndCPType = 0;

                    if (IsOk)
                    {
                        char *Value;

                        // Bond number
                        if (TecUtilAuxDataGetStrItemByName(AuxDataRef, "CompChem.BegCrtPtNum",
                                                           &Value, &Retain))
                        {
                            if (Value != NULL)
                            {
                                /* Parse the value */
                                if (sscanf(Value, "%d", &BondCrtPtNum) == 1)
                                {
                                    BondCrtPtNum--;
                                }
                                else
                                {
                                    IsOk = FALSE;
                                }
                                // release the allocated string copy
                                TecUtilStringDealloc(&Value);
                            }
                            else
                            {
                                IsOk = FALSE;
                            }
                        }
                    }


                    // A beginning type of "Bond" was found, read the IsFarField aux data
                    if (TecUtilAuxDataGetItemByName(AuxDataRef, "CompChem.IsFarField",
                                                    &StrValue, &Type, &Retain))
                    {
                        IsFarFieldString = (char *)StrValue;
                        if (strcmp(IsFarFieldString, "TRUE") == 0) IsFarField = TRUE;
                    }


                    // A beginning type of "Bond" was found, look for end type of "Cage"
                    if (TecUtilAuxDataGetItemByName(AuxDataRef, "CompChem.EndCrtPtType",
                                                    &StrValue, &Type, &Retain))
                    {
                        EndTypeString = (char *)StrValue;
                        if (strcmp(EndTypeString, "Cage") == 0)
                            EndCPType = 3;
                        else if (strcmp(EndTypeString, "Ring") == 0)
                            EndCPType = 1;

                        if (EndCPType == 3 || EndCPType == 1)
                        {
                            char         *Value;
                            // Cage or Ring number
                            if (IsOk && TecUtilAuxDataGetStrItemByName(AuxDataRef, "CompChem.EndCrtPtNum",
                                                                       &Value, &Retain))
                            {
                                if (Value != NULL)
                                {
                                    LgIndex_t CageOrRingCrtPtNum;
                                    /* Parse the value */
                                    if (sscanf(Value, "%d", &CageOrRingCrtPtNum) == 1)
                                    {
                                        if (EndCPType == 3)
                                            CageCrtPtNum = CageOrRingCrtPtNum - 1;
                                        else if (EndCPType == 1)
                                            RingCrtPtNum = CageOrRingCrtPtNum - 1;
                                    }
                                    else
                                    {
                                        IsOk = FALSE;
                                    }
                                    // release the allocated string copy
                                    TecUtilStringDealloc(&Value);
                                }
                                else
                                {
                                    IsOk = FALSE;
                                }
                            }
                        }
                    }



                    // Beg and End types of "Bond" and "Cage" (or "Ring") were found, look for thrd type of "Ring"
                    if (TecUtilAuxDataGetItemByName(AuxDataRef, "CompChem.ThrdCrtPtType",
                                                    &StrValue, &Type, &Retain))
                    {
                        ThrdTypeString = (char *)StrValue;
                        if (strcmp(ThrdTypeString, "Ring") == 0)
                        {
                            char         *Value;

                            // Ring number
                            if (IsOk && TecUtilAuxDataGetStrItemByName(AuxDataRef, "CompChem.ThrdCrtPtNum",
                                                                       &Value, &Retain))
                            {
                                if (Value != NULL)
                                {
                                    /* Parse the value */
                                    if (EndCPType == 3)
                                    {
                                        if (sscanf(Value, "%d", &RingCrtPtNum) == 1)
                                            RingCrtPtNum--;
                                        else
                                            IsOk = FALSE;
                                    }
                                    else if (EndCPType == 1)
                                    {
                                        if (sscanf(Value, "%d", &Ring2CrtPtNum) == 1)
                                            Ring2CrtPtNum--;
                                        else
                                            IsOk = FALSE;
                                    }

                                    // release the allocated string copy
                                    TecUtilStringDealloc(&Value);
                                }
                                else
                                {
                                    IsOk = FALSE;
                                }
                            }
                        }
                    }


                    // Now find the two Atoms that connect to the Bond, and create two Bundles
                    if (IsOk)

                    {
                        ArrListItem_u Item;
                        LgIndex_t Atom1, Atom2;
                        LgIndex_t BondNum = BondCrtPtNum - CritPointsGetBegOffset(CritPoints, (char)(-1));
                        LgIndex_t RingNum = -1;
                        LgIndex_t CageNum = -1;
                        LgIndex_t Ring2Num = -1;
                        LgIndex_t CageOrRingNum;
                        char      FarField = 0;
                        char      CageOrRingType = 3;  // Default to Cage

                        if (IsFarField) FarField = 1;

                        if (RingCrtPtNum > -1) RingNum = RingCrtPtNum - CritPointsGetBegOffset(CritPoints, (char)(1));
                        if (CageCrtPtNum > -1) CageNum = CageCrtPtNum - CritPointsGetBegOffset(CritPoints, (char)(3));
                        if (Ring2CrtPtNum > -1) Ring2Num = Ring2CrtPtNum - CritPointsGetBegOffset(CritPoints, (char)(1));

                        // Set CageOrRingNum appropriately
                        CageOrRingNum = CageNum;
                        if (CageNum == -1 && Ring2Num >= 0)
                        {
                            CageOrRingNum  = Ring2Num;
                            CageOrRingType = 1;
                        }

                        // Need either a far-field bundle, or Cage and Ring numbers defined
                        if (BondNum >= 0 && (IsFarField || (CageNum >= 0 && RingNum >= 0)))
                        {
                            if (Bundles == NULL) Bundles = BundlesAlloc();

                            Item = ArrListGetItem(Atom1ForBond, BondNum);
                            Atom1 = Item.Long;

                            // Create Bundle
                            IsOk = BundlesAppendAtEnd(Bundles, Atom1, BondNum, RingNum, CageOrRingNum,
                                                      CageOrRingType, FarField);

                            Item = ArrListGetItem(Atom2ForBond, BondNum);
                            Atom2 = Item.Long;
                            // Create Bundle
                            IsOk = BundlesAppendAtEnd(Bundles, Atom2, BondNum, RingNum, CageOrRingNum,
                                                      CageOrRingType, FarField);
                        }

                    }
                    // release the allocated string copy
                    if (ThrdTypeString != NULL) TecUtilStringDealloc(&ThrdTypeString);
                    // release the allocated string copy
                    if (EndTypeString != NULL) TecUtilStringDealloc(&EndTypeString);
                }
                // release the allocated string copy
                if (BegTypeString != NULL) TecUtilStringDealloc(&BegTypeString);
            }
        }
    }




    // if(ShowStatusBar) TecUtilStatusFinishPercentDone();

    *CritPointsPtr = CritPoints;
    *BundlesPtr = Bundles;

    ENSURE(VALID_BOOLEAN(IsOk));
    return(IsOk);
}







/**
 * Find if the item is in the long Array List.
 *
 * param ArrList
 *     An array list to search.
 * param Item
 *     Item you are search for.
 * return
 *     TRUE if it is in the list, FALSE otherwise.
 */
Boolean_t ItemIsInLongArrList(ArrList_pa    ArrList,
                              LgIndex_t     Item)
{
    Boolean_t InList = FALSE;
    int       ii, Count;

    REQUIRE(ArrListIsValid(ArrList));

    Count = ArrListGetCount(ArrList);

    for (ii = 0; !InList && ii < Count; ii++)
    {
        ArrListItem_u ItemInList = ArrListGetItem(ArrList, ii);
        if (Item == ItemInList.Long) InList = TRUE;
    }

    ENSURE(VALID_BOOLEAN(InList));
    return InList;
}




/**
 * Find and return the list of Atoms connected by any of the specified
 * list of Bonds.
 *
 * param IsAll
 *     Return a list of all atoms, ignore other inputs.
 *
 * param BondList
 *     Array list containing the Bond numbers which have a connection
 *     to the desired atoms.
 *
 * return
 *     Array list of atom numbers, NULL otherwise.
 */
ArrList_pa ListOfAtomsForBonds(Boolean_t  IsAll,
                               ArrList_pa BondList,
                               void      *BundlesVoid)
{
    ArrList_pa Result = NULL;
    Boolean_t IsOk = TRUE;
    Bundles_pa Bundles = (Bundles_pa)BundlesVoid;
    LgIndex_t NumBundles = BundlesGetCount(Bundles);
    LgIndex_t Atom, Bond, Ring, CageOrRing;
    char CageOrRingType, BundleType;
    int ii;

    for (ii = 0; IsOk && ii < NumBundles; ii++)
    {
        IsOk = BundlesGetFromOffset(Bundles, ii, &Atom, &Bond, &Ring, &CageOrRing, &CageOrRingType, &BundleType);
        if (IsOk && (IsAll || ItemIsInLongArrList(BondList, Bond)))
        {

            if (Result == NULL)
            {
                ArrListItem_u Item;
                Result = ArrListAlloc(20, ArrListType_Long);
                Item.Long = Atom;
                IsOk = ArrListAppendItem(Result, Item);
            }
            else
            {
                int jj;
                LgIndex_t AtomListCount = 0;
                Boolean_t IsFound = FALSE;

                AtomListCount = ArrListGetCount(Result);
                for (jj = 0; !IsFound && jj < AtomListCount; jj++)
                {
                    ArrListItem_u Item;
                    Item = ArrListGetItem(Result, jj);
                    if (Item.Long == Atom) IsFound = TRUE;  // Duplicate

                    // Found it, insert the atom
                    if (Item.Long > Atom)
                    {
                        IsFound = TRUE;
                        Item.Long = Atom;
                        IsOk = ArrListInsertItem(Result, jj, Item);
                    }
                }
                // Must be larger than the existing items in the list - append
                if (!IsFound)
                {
                    ArrListItem_u Item;
                    Item.Long = Atom;
                    IsOk = ArrListAppendItem(Result, Item);
                }
            }
        }
    }

    ENSURE(VALID_BOOLEAN(IsOk));
    ENSURE(Result == NULL || ArrListIsValid(Result));
    return Result;
}




/**
 * For ArrLists of type Long, insert an item into the array list in
 * ascending (sorted) order.
 *
 * param ArrListLong
 *     ArrList in which to insert item.
 * param LongItem
 *     Long integer to insert.
 *
 * return
 *     TRUE if it worked, NULL otherwise.
 */
Boolean_t ArrListLongSortedInsert(ArrList_pa ArrListLong,
                                  long       LongItem)
{
    int jj;
    Boolean_t IsOk = TRUE;
    LgIndex_t Count = 0;
    Boolean_t IsFound = FALSE;
    ArrListItem_u Item;

    REQUIRE(ArrListIsValid(ArrListLong));
    REQUIRE(ArrListGetType(ArrListLong) == ArrListType_Long);

    // Insert the long interger into the list in order
    Count = ArrListGetCount(ArrListLong);
    for (jj = 0; IsOk && !IsFound && jj < Count; jj++)
    {
        Item = ArrListGetItem(ArrListLong, jj);
        if (Item.Long == LongItem) IsFound = TRUE;  // Duplicate

        // Found it, insert the CPIn critical point number
        if (Item.Long > LongItem)
        {
            IsFound = TRUE;
            Item.Long = LongItem;
            IsOk = ArrListInsertItem(ArrListLong, jj, Item);
        }
    }
    // Must be larger than the existing items in the list - append
    if (!IsFound)
    {
        Item.Long = LongItem;
        IsOk = ArrListAppendItem(ArrListLong, Item);
    }

    ENSURE(VALID_BOOLEAN(IsOk));
    return IsOk;
}




/**
 * Find and return the list of CPTypeOut critical points connected to
 * any of the specified CPTypeIn critical points.
 *
 * param CPTypeIn
 *     Type of the input critical points (CPInList).
 * param CPTypeOut
 *     Type of the desired output critical point list.
 *
 * param CPInList
 *     Array list containing the CPTypeIn critical point numbers which
 *     have a connection to the desired CPTypeOut critical points.
 *
 * return
 *     Array list of CPTypeOut critical point numbers, NULL otherwise.
 */
ArrList_pa ListOfCPOutsConnetedToCPIns(char       CPTypeIn,
                                       char       CPTypeOut,
                                       ArrList_pa CPInList,
                                       void      *BundlesVoid)
{
    ArrList_pa Result = NULL;
    Boolean_t IsOk = TRUE;
    Bundles_pa Bundles = (Bundles_pa)BundlesVoid;
    LgIndex_t NumBundles = BundlesGetCount(Bundles);
    LgIndex_t Atom, Bond, Ring, CageOrRing;
    char CageOrRingCPType, BundleType;
    int ii;

    REQUIRE(CPTypeIn == -3 || CPTypeIn == -1 || CPTypeIn == 1 || CPTypeIn == 3);
    REQUIRE(CPTypeOut == -3 || CPTypeOut == -1 || CPTypeOut == 1 || CPTypeOut == 3);
    REQUIRE(ArrListIsValid(CPInList));
    REQUIRE(BundlesIsValid(Bundles));


    // Search through Bundles looking for CPTypeOut critical points connnected to
    // the CPTypeIn critical points in CPInList.
    for (ii = 0; IsOk && ii < NumBundles; ii++)
    {
        LgIndex_t CPInNum[2];
        int jj;

        CPInNum[0] = -1;
        CPInNum[1] = -1;

        IsOk = BundlesGetFromOffset(Bundles, ii, &Atom, &Bond, &Ring, &CageOrRing, &CageOrRingCPType, &BundleType);

        // Determine the CPTypeIn critical point number
        switch ((int)CPTypeIn)
        {
            case -3:
                CPInNum[0] = Atom;
                break;
            case -1:
                CPInNum[0] = Bond;
                break;
            case 1:
                CPInNum[0] = Ring;
                if (CageOrRingCPType == 1)
                    CPInNum[1] = CageOrRing;
                break;
            case 3:
                if (CageOrRingCPType == 3)
                    CPInNum[0] = CageOrRing;
                break;
        }

        for (jj = 0; IsOk && jj < 2; jj++)
        {
            if (CPInNum[jj] > -1 && ItemIsInLongArrList(CPInList, CPInNum[jj]))
            {
                LgIndex_t CPOutNum = -1;

                // Determine the CPTypeOut critical point number
                switch ((int)CPTypeOut)
                {
                    case -3:
                        CPOutNum = Atom;
                        break;
                    case -1:
                        CPOutNum = Bond;
                        break;
                    case 1:
                        CPOutNum = Ring;
                        break;
                    case 3:
                        if (CageOrRingCPType == 3)
                            CPOutNum = CageOrRing;
                        break;
                }

                // First entry into CPOutList
                if (Result == NULL)
                {
                    ArrListItem_u Item;
                    Result = ArrListAlloc(20, ArrListType_Long);
                    if (CPOutNum >= 0)
                    {
                        Item.Long = CPOutNum;
                        IsOk = ArrListAppendItem(Result, Item);
                    }
                }
                else
                {
                    LgIndex_t CPOutListCount = 0;
                    Boolean_t IsFound = FALSE;

                    // Insert the CPOutType critical points into the list in order
                    if (CPOutNum >= 0)
                    {
                        IsOk = ArrListLongSortedInsert(Result, CPOutNum);
                    }

                    // Special case for Ring output when CageOrRingType is Ring
                    if (CPTypeOut == 1 && CageOrRingCPType == 1 && CageOrRing > -1)
                    {
                        IsOk = ArrListLongSortedInsert(Result, CageOrRing);
                    }
                }
            }
        }
    }

    ENSURE(VALID_BOOLEAN(IsOk));
    ENSURE(Result == NULL || ArrListIsValid(Result));
    return Result;
}







/**
 * Find and return the list of zones (zone numbers) for bond lines
 * containing CPType critical points in specified CPInList.
 *
 * param CPType
 *     Type of the input critical points (CPInList) - either Bond or Atom.
 *
 * param CPInList
 *     Array list containing the CPTypeIn (Bond or Atom) critical point
 *     numbers which have a connection to the desired bond lines.
 *
 * return
 *     Array list of bond line zone numbers if it works, NULL otherwise.
 */
ArrList_pa ZoneListOfBondLinesForCritPoint(char          CritPointType,
                                           ArrList_pa    CritPointInList,
                                           CritPoints_pa CritPoints,
                                           void         *BundlesVoid)
{
    ArrList_pa Result = ArrListAlloc(20, ArrListType_Long);
    Bundles_pa Bundles = (Bundles_pa)BundlesVoid;
    Boolean_t  IsOk = TRUE;
    EntIndex_t NumZones;
    EntIndex_t zn;
    char       CPType = CritPointType;
    ArrList_pa CPInList = CritPointInList;

    REQUIRE(CPType == -3 || CPType == -1 || CPType == 1 || CPType == 3);  // Only Bond or Atom
    REQUIRE(ArrListIsValid(CPInList));
    REQUIRE(CritPointsIsValid(CritPoints));
    REQUIRE(BundlesIsValid(Bundles));

    /* Get num zones in dataset */
    if (IsOk)
    {
        char *DatasetTitle = NULL;
        EntIndex_t NumVars;
        TecUtilDataSetGetInfo(&DatasetTitle, &NumZones, &NumVars);
        TecUtilStringDealloc(&DatasetTitle);
    }

    // For Rings or Cages, must first build a list of Bonds that share Bundles with the Rings or Cages
    if (CPType == 1 || CPType == 3)
    {
        LgIndex_t bd;
        LgIndex_t NumBundles = BundlesGetCount(Bundles);

        CPInList = ArrListAlloc(60, ArrListType_Long);

        for (bd = 0; bd < NumBundles; bd++)
        {
            LgIndex_t Atom, Bond, Ring, CageOrRing;
            char  CageOrRingCPType, BundleType;
            IsOk = BundlesGetFromOffset(Bundles,  bd, &Atom, &Bond, &Ring, &CageOrRing, &CageOrRingCPType, &BundleType);
            if (CPType == 1)
            {
                LgIndex_t cpn;
                LgIndex_t NumCPInList = ArrListGetCount(CritPointInList);
                for (cpn = 0; cpn < NumCPInList; cpn++)
                {
                    ArrListItem_u Item = ArrListGetItem(CritPointInList, cpn);
                    LgIndex_t RingInList = Item.Long;
                    if (RingInList == Ring || (CageOrRingCPType == 1 && CageOrRing == RingInList))
                    {
                        IsOk = ArrListLongSortedInsert(CPInList, Bond);
                    }
                }
            }
            else if (CPType == 3)
            {
                LgIndex_t cpn;
                LgIndex_t NumCPInList = ArrListGetCount(CritPointInList);
                for (cpn = 0; cpn < NumCPInList; cpn++)
                {
                    ArrListItem_u Item = ArrListGetItem(CritPointInList, cpn);
                    LgIndex_t CageInList = Item.Long;
                    if (CageInList == CageOrRing && CageOrRingCPType == 3)
                    {
                        IsOk = ArrListLongSortedInsert(CPInList, Bond);
                    }
                }
            }
        }
        CPType = -1;
    }


    // Build Bond-Atom list by searhing through the zone
    // aux data. We know how many bonds there are, and that there
    // should be two atoms per bond.
    for (zn = 1; zn <= NumZones; zn++)
    {
        // If zone is a Bond-Atom zone it has auxiliary
        // data items named
        //     "CompChem.BegCrtPtType" == "Bond"
        //     "CompChem.EndCrtPtType" == "Atom"
        AuxData_pa AuxDataRef = TecUtilAuxDataZoneGetRef(zn);
        if (AuxDataRef != NULL)
        {
            char *BegTypeString;
            ArbParam_t StrValue;
            AuxDataType_e Type;
            Boolean_t     Retain;
            if (TecUtilAuxDataGetItemByName(AuxDataRef, "CompChem.BegCrtPtType",
                                            &StrValue, &Type, &Retain))
            {
                BegTypeString = (char *)StrValue;
                if (strcmp(BegTypeString, "Bond") == 0)
                {
                    char *EndTypeString;

                    // A beginning type of "Bond" was found, look for end type of "Cage"
                    if (TecUtilAuxDataGetItemByName(AuxDataRef, "CompChem.EndCrtPtType",
                                                    &StrValue, &Type, &Retain))
                    {
                        EndTypeString = (char *)StrValue;
                        if (strcmp(EndTypeString, "Atom") == 0)
                        {
                            // Beg and End types of "Bond" and "Atom" were found
                            char         *Value;

                            // Get the "Bond" or "Atom" critical point number and compare to
                            // That in CPInList
                            //
                            if ((int)CPType == -1)  // Bond
                            {

                                // Bond number
                                if (TecUtilAuxDataGetStrItemByName(AuxDataRef, "CompChem.BegCrtPtNum",
                                                                   &Value, &Retain))
                                {
                                    LgIndex_t BondCrtPtNum;
                                    LgIndex_t BondListCount = ArrListGetCount(CPInList);
                                    LgIndex_t BondNum;
                                    int       jj;

                                    if (Value != NULL)
                                    {
                                        /* Parse the value */
                                        BondCrtPtNum = atoi(Value) - 1;  // Convert to zero-based

                                        // release the allocated string copy
                                        TecUtilStringDealloc(&Value);
                                    }
                                    else
                                    {
                                        IsOk = FALSE;
                                    }

                                    BondNum = BondCrtPtNum - CritPointsGetBegOffset(CritPoints, (char)(-1));

                                    // If it is in the CPInList, add zone number to Result list
                                    for (jj = 0; jj < BondListCount; jj++)
                                    {
                                        ArrListItem_u Item = ArrListGetItem(CPInList, jj);
                                        LgIndex_t BondInList = Item.Long;
                                        if (BondInList == BondNum)
                                        {
                                            ArrListItem_u Item;
                                            Item.Long = zn;
                                            IsOk = ArrListAppendItem(Result, Item);
                                        }
                                    }
                                }
                            }
                            else if ((int)CPType == -3)  // Atom
                            {
                                // Atom number
                                if (IsOk && TecUtilAuxDataGetStrItemByName(AuxDataRef, "CompChem.EndCrtPtNum",
                                                                           &Value, &Retain))
                                {
                                    LgIndex_t AtomCrtPtNum;
                                    LgIndex_t AtomListCount = ArrListGetCount(CPInList);
                                    int       jj;

                                    if (Value != NULL)
                                    {
                                        /* Parse the value */
                                        AtomCrtPtNum = atoi(Value) - 1;  // Convert to zero-based

                                        // release the allocated string copy
                                        TecUtilStringDealloc(&Value);
                                    }
                                    else
                                    {
                                        IsOk = FALSE;
                                    }

                                    // If it is in the CPInList, add zone number to Result list
                                    for (jj = 0; jj < AtomListCount; jj++)
                                    {
                                        ArrListItem_u Item = ArrListGetItem(CPInList, jj);
                                        LgIndex_t AtomInList = Item.Long;
                                        if (AtomInList == AtomCrtPtNum)
                                        {
                                            ArrListItem_u Item;
                                            Item.Long = zn;
                                            IsOk = ArrListAppendItem(Result, Item);
                                        }
                                    }
                                }
                            }
                        }
                        // release the allocated string copy
                        if (EndTypeString != NULL) TecUtilStringDealloc(&EndTypeString);
                    }
                }
                // release the allocated string copy
                if (BegTypeString != NULL) TecUtilStringDealloc(&BegTypeString);
            }
        }
    }

    // For Rings or Cages, dealloc the temporary list of Bonds
    if (CPType == 1 || CPType == 3)
        ArrListDealloc(&CPInList);

    ENSURE(ArrListIsValid(Result));
    return Result;
}









/**
 * Find and return the list of zones (zone numbers) for edge lines of bundles
 * containing CPType critical points in specified CPInList.
 *
 * param GetEdges, GetSurfaces
 *     Include zones that are edges (GetEdges) or surfaces (GetSurfaces) of
 *     the selected bundles.
 * param CPType
 *     Type of the input critical points (CPInList) - either Bond or Atom.
 *
 * param CPInList
 *     Array list containing the CPTypeIn (Bond or Atom) critical point
 *     numbers which have a connection to the desired bond lines.
 * param CritPoints
 *     Structure containing all critical points and related information
 * param BundlesVoid
 *     Void pointer to structure defining all bundles
 *
 * return
 *     Array list of bond line zone numbers if it works, NULL otherwise.
 */
ArrList_pa ZoneListOfBundleZonesForCritPoint(Boolean_t     GetEdges,
                                             Boolean_t     GetSurfaces,
                                             char          CPType,
                                             ArrList_pa    CPList,
                                             CritPoints_pa CritPoints,
                                             void         *BundlesVoid)
{
    ArrList_pa Result = ArrListAlloc(20, ArrListType_Long);
    Bundles_pa Bundles = (Bundles_pa)BundlesVoid;
    ArrList_pa BundleList = ArrListAlloc(60, ArrListType_Long);
    Boolean_t  IsOk = TRUE;
    EntIndex_t NumZones;
    EntIndex_t zn;

    REQUIRE(CPType == -3 || CPType == -1 || CPType == 1 || CPType == 3);  // Only Bond or Atom
    REQUIRE(ArrListIsValid(CPList));
    REQUIRE(CritPointsIsValid(CritPoints));
    REQUIRE(BundlesIsValid(Bundles));

    /* Get num zones in dataset */
    if (IsOk)
    {
        char *DatasetTitle = NULL;
        EntIndex_t NumVars;
        TecUtilDataSetGetInfo(&DatasetTitle, &NumZones, &NumVars);
        TecUtilStringDealloc(&DatasetTitle);
    }

    // First build a list of Bundles numbers with the relevant Critical Point
    if (IsOk && (GetEdges || GetSurfaces))
    {
        LgIndex_t bd;
        LgIndex_t NumBundles = BundlesGetCount(Bundles);

        for (bd = 0; bd < NumBundles; bd++)
        {
            LgIndex_t CPNumBundle = -1;
            LgIndex_t Atom, Bond, Ring, CageOrRing;
            char CageOrRingCPType, BundleType;
            IsOk = BundlesGetFromOffset(Bundles,  bd, &Atom, &Bond, &Ring, &CageOrRing, &CageOrRingCPType, &BundleType);
            switch (CPType)
            {
                case -3:
                    CPNumBundle = Atom;
                    break;
                case -1:
                    CPNumBundle = Bond;
                    break;
                case  1:
                    CPNumBundle = Ring;
                    break;
                case  3:
                    if (CageOrRingCPType == 3)
                        CPNumBundle = CageOrRing;
                    break;
                default:
                    CPNumBundle = -1;
            }


            if (CPNumBundle >= 0)
            {
                LgIndex_t cpn;
                LgIndex_t NumCPInList = ArrListGetCount(CPList);
                for (cpn = 0; cpn < NumCPInList; cpn++)
                {
                    ArrListItem_u Item = ArrListGetItem(CPList, cpn);
                    LgIndex_t CPInList = Item.Long;
                    if (CPInList == CPNumBundle)
                    {
                        IsOk = ArrListLongSortedInsert(BundleList, bd);
                    }
                }
            }
            // Special case for far-field bundles with second ring
            if (CPType == Ring && CageOrRingCPType == 1 && CageOrRing > -1)
            {
                LgIndex_t cpn;
                LgIndex_t NumCPInList = ArrListGetCount(CPList);
                for (cpn = 0; cpn < NumCPInList; cpn++)
                {
                    ArrListItem_u Item = ArrListGetItem(CPList, cpn);
                    LgIndex_t CPInList = Item.Long;
                    if (CPInList == CageOrRing)
                    {
                        IsOk = ArrListLongSortedInsert(BundleList, bd);
                    }
                }
            }
        }
    }


    // Build ZoneList by searching through the zone aux data and comparing beginning
    // and ending critical points to the critical point pairs of bundle edge lines.
    if (IsOk && (GetEdges || GetSurfaces))
    {
        for (zn = 1; IsOk && zn <= NumZones; zn++)
        {
            AuxData_pa AuxDataRef = TecUtilAuxDataZoneGetRef(zn);

            // Skip if not 1-D zone (JMax = 1)
            LgIndex_t IMax, JMax, KMax;
            ZoneType_e ZoneType = TecUtilZoneGetType(zn);
            Boolean_t ZnIsOrdered = TecUtilZoneIsOrdered(zn);
            Boolean_t IsLine, IsSurface;

            TecUtilZoneGetInfo(zn, &IMax, &JMax, &KMax, NULL, NULL, NULL, NULL, NULL,
                               NULL, NULL, NULL, NULL, NULL);


            IsLine    = ZnIsOrdered && IMax > 1 && JMax <= 1 && KMax <= 1;

            IsSurface = (ZnIsOrdered && IMax > 1 && JMax > 1 && KMax <= 1) ||
                        ZoneType == ZoneType_FETriangle || ZoneType == ZoneType_FEQuad;

            if (AuxDataRef != NULL && (IsLine || IsSurface))
            {
                char *BegTypeString = NULL;
                char *EndTypeString = NULL;
                char *ThrdTypeString = NULL;
                ArbParam_t StrValue;
                AuxDataType_e Type;
                Boolean_t     Retain;
                char BegCrtPtType = 0;
                char EndCrtPtType = 0;
                char ThrdCrtPtType = 0;
                LgIndex_t BegCrtPtNum = -1;
                LgIndex_t EndCrtPtNum = -1;
                LgIndex_t ThrdCrtPtNum = -1;

                if (TecUtilAuxDataGetItemByName(AuxDataRef, "CompChem.BegCrtPtType",
                                                &StrValue, &Type, &Retain))
                {
                    char *Value;

                    BegTypeString = (char *)StrValue;
                    if (strcmp(BegTypeString, "Atom") == 0)
                        BegCrtPtType = -3;
                    else if (strcmp(BegTypeString, "Bond") == 0)
                        BegCrtPtType = -1;
                    else if (strcmp(BegTypeString, "Ring") == 0)
                        BegCrtPtType =  1;
                    else if (strcmp(BegTypeString, "Cage") == 0)
                        BegCrtPtType =  3;

                    // Get beginning critical point number
                    if (TecUtilAuxDataGetStrItemByName(AuxDataRef, "CompChem.BegCrtPtNum",
                                                       &Value, &Retain))
                    {
                        if (Value != NULL)
                        {
                            /* Parse the value */
                            BegCrtPtNum = atoi(Value) - 1;  // Convert to zero-based

                            // release the allocated string copy
                            TecUtilStringDealloc(&Value);
                        }
                        else
                        {
                            IsOk = FALSE;
                        }
                    }
                }
                // release the allocated string
                if (BegTypeString != NULL) TecUtilStringDealloc(&BegTypeString);


                if (TecUtilAuxDataGetItemByName(AuxDataRef, "CompChem.EndCrtPtType",
                                                &StrValue, &Type, &Retain))
                {
                    char *Value;

                    EndTypeString = (char *)StrValue;
                    if (strcmp(EndTypeString, "Atom") == 0)
                        EndCrtPtType = -3;
                    else if (strcmp(EndTypeString, "Bond") == 0)
                        EndCrtPtType = -1;
                    else if (strcmp(EndTypeString, "Ring") == 0)
                        EndCrtPtType =  1;
                    else if (strcmp(EndTypeString, "Cage") == 0)
                        EndCrtPtType =  3;

                    // Get beginning critical point number
                    if (TecUtilAuxDataGetStrItemByName(AuxDataRef, "CompChem.EndCrtPtNum",
                                                       &Value, &Retain))
                    {
                        if (Value != NULL)
                        {
                            /* Parse the value */
                            EndCrtPtNum = atoi(Value) - 1;  // Convert to zero-based

                            // release the allocated string copy
                            TecUtilStringDealloc(&Value);
                        }
                        else
                        {
                            IsOk = FALSE;
                        }
                    }

                }
                // release the allocated string
                if (EndTypeString != NULL) TecUtilStringDealloc(&EndTypeString);


                if (IsSurface && TecUtilAuxDataGetItemByName(AuxDataRef, "CompChem.ThrdCrtPtType",
                                                             &StrValue, &Type, &Retain))
                {
                    char *Value;

                    ThrdTypeString = (char *)StrValue;
                    if (strcmp(ThrdTypeString, "Atom") == 0)
                        ThrdCrtPtType = -3;
                    else if (strcmp(ThrdTypeString, "Bond") == 0)
                        ThrdCrtPtType = -1;
                    else if (strcmp(ThrdTypeString, "Ring") == 0)
                        ThrdCrtPtType =  1;
                    else if (strcmp(ThrdTypeString, "Cage") == 0)
                        ThrdCrtPtType =  3;

                    // Get beginning critical point number
                    if (TecUtilAuxDataGetStrItemByName(AuxDataRef, "CompChem.ThrdCrtPtNum",
                                                       &Value, &Retain))
                    {
                        if (Value != NULL)
                        {
                            /* Parse the value */
                            ThrdCrtPtNum = atoi(Value) - 1;  // Convert to zero-based

                            // release the allocated string copy
                            TecUtilStringDealloc(&Value);
                        }
                        else
                        {
                            IsOk = FALSE;
                        }
                    }

                }
                // release the allocated string
                if (ThrdTypeString != NULL) TecUtilStringDealloc(&ThrdTypeString);


                // Include zone in list is an edge or surface of a selected bundle
                if (IsOk && BegCrtPtType != 0)
                {
                    // IsLine: Include zone in list if it is an edge of a selected bundle
                    if (IsLine && GetEdges)
                    {
                        Boolean_t IsFound = FALSE;
                        LgIndex_t bd;
                        LgIndex_t NumBundles = ArrListGetCount(BundleList);
                        for (bd = 0; IsOk && !IsFound && bd < NumBundles; bd++)
                        {
                            ArrListItem_u Item = ArrListGetItem(BundleList, bd);
                            LgIndex_t BundleNum = Item.Long;

                            LgIndex_t BegCrtPtNumDel = BegCrtPtNum - CritPointsGetBegOffset(CritPoints, BegCrtPtType);
                            LgIndex_t EndCrtPtNumDel = -1;

                            if (EndCrtPtType != 0)
                                EndCrtPtNumDel = EndCrtPtNum - CritPointsGetBegOffset(CritPoints, EndCrtPtType);

                            if (BundlesIsEdge(Bundles, BundleNum, BegCrtPtType, BegCrtPtNumDel, EndCrtPtType, EndCrtPtNumDel))
                            {
                                IsFound = TRUE;
                                Item.Long = zn;
                                IsOk = ArrListAppendItem(Result, Item);
                            }
                        }
                    }
                    // IsSurface: Include zone in list if it is an surface of a selected bundle
                    else if (IsSurface && GetSurfaces)
                    {
                        Boolean_t IsFound = FALSE;
                        LgIndex_t bd;
                        LgIndex_t NumBundles = ArrListGetCount(BundleList);
                        for (bd = 0; IsOk && !IsFound && bd < NumBundles; bd++)
                        {
                            ArrListItem_u Item = ArrListGetItem(BundleList, bd);
                            LgIndex_t BundleNum = Item.Long;

                            LgIndex_t BegCrtPtNumDel = BegCrtPtNum - CritPointsGetBegOffset(CritPoints, BegCrtPtType);
                            LgIndex_t EndCrtPtNumDel = -1;
                            LgIndex_t ThrdCrtPtNumDel = -1;

                            if (EndCrtPtType != 0)
                                EndCrtPtNumDel = EndCrtPtNum - CritPointsGetBegOffset(CritPoints, EndCrtPtType);
                            if (ThrdCrtPtType != 0)
                                ThrdCrtPtNumDel = ThrdCrtPtNum - CritPointsGetBegOffset(CritPoints, ThrdCrtPtType);

                            if (BundlesIsSurface(Bundles, BundleNum, BegCrtPtType, BegCrtPtNumDel,
                                                 EndCrtPtType, EndCrtPtNumDel, ThrdCrtPtType, ThrdCrtPtNumDel))
                            {
                                IsFound = TRUE;
                                Item.Long = zn;
                                IsOk = ArrListAppendItem(Result, Item);
                            }
                        }
                    }
                }
            }
        }
    }

    // Dealloc the temporary list of Bundles
    if (BundleList != NULL)
        ArrListDealloc(&BundleList);

    ENSURE(ArrListIsValid(Result));
    return Result;
}




