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
#include "ADDGLBL.h"
#include "ADDONVER.h"
#include "TASSERT.h"
#include "GUIDEFS.h"
#include "ENGINE.h"
#include "ADKUTIL.h"
#include "ALLOC.h"
#include "ARRLIST.h"
#include "NORMALS.h"
#include "SURFELEMMAP.h"
#include "ZONEVARINFO.h"
#include "CRITPOINTS.h"
#include "GRADPATH.h"
#include <string.h>
#include "BUNDLES.h"



/**
 * Determine if the Bundles handle is sane.
 *
 * param Bundles
 *     Bundles structure in question.
 *
 * return
 *     TRUE if the Bundles structure is valid, otherwise FALSE.
 */
Boolean_t BundlesIsValid(Bundles_pa Bundles)
{
    Boolean_t IsValid = FALSE;

    IsValid = (VALID_REF(Bundles) &&
               VALID_REF(Bundles->Atom) && ArrListIsValid(Bundles->Atom) &&
               VALID_REF(Bundles->Bond) && ArrListIsValid(Bundles->Bond) &&
               VALID_REF(Bundles->Ring) && ArrListIsValid(Bundles->Ring) &&
               VALID_REF(Bundles->Cage) && ArrListIsValid(Bundles->Cage));

    /* Require the same count for each array list in Bundles structure. */
    if (IsValid)
    {
        LgIndex_t Count = ArrListGetCount(Bundles->Atom);
        IsValid = (ArrListGetCount(Bundles->Ring) == Count);
        if (IsValid) IsValid = (ArrListGetCount(Bundles->Bond) == Count);
        if (IsValid) IsValid = (ArrListGetCount(Bundles->Cage) == Count);
    }

    ENSURE(VALID_BOOLEAN(IsValid));
    return IsValid;
}




/**
 * Deallocates the Bundles handle and set the handle to NULL.
 *
 * note
 *     Item destruction is the responsibility of the caller.
 *
 * param
 *     Reference to a Bundles handle.
 */
void BundlesDealloc(Bundles_pa *Bundles)
{
    REQUIRE(VALID_REF(Bundles));
    REQUIRE(BundlesIsValid(*Bundles) || *Bundles == NULL);

    if (*Bundles != NULL)
    {
        /* release the ArrList's */
        if ((*Bundles)->Atom != NULL) ArrListDealloc(&((*Bundles)->Atom));
        if ((*Bundles)->Bond != NULL) ArrListDealloc(&((*Bundles)->Bond));
        if ((*Bundles)->Ring != NULL) ArrListDealloc(&((*Bundles)->Ring));
        if ((*Bundles)->Cage != NULL) ArrListDealloc(&((*Bundles)->Cage));

        /* Notify Tecplot of memory usage change */
        if ((*Bundles)->MemUsageReported != 0)
        {
            TecUtilMemoryChangeNotify(- (Int64_t)((*Bundles)->MemUsageReported));
            (*Bundles)->MemUsageReported = 0;
        }

        /* release the list structure itself */
        FREE_ITEM(*Bundles, "Bundles structure");
        *Bundles = NULL;
    }

    ENSURE(*Bundles == NULL);
}





/**
 * Empties the Bundles structure of all Bundles.
 *
 * param Bundles
 *     Bundles to clear.
 */
void BundlesClear(Bundles_pa Bundles)
{
    REQUIRE(BundlesIsValid(Bundles));

    ArrListClear(Bundles->Atom);
    ArrListClear(Bundles->Bond);
    ArrListClear(Bundles->Ring);
    ArrListClear(Bundles->Cage);

    ENSURE(BundlesIsValid(Bundles) && BundlesGetCount(Bundles) == 0);
}





/**
 * Removes a point from the Bundles array. The members following the item
 * removed are shifted down accordingly to fill the vacated space.
 *
 *
 * param Bundles
 *     Bundles structure containing the bundle to remove.
 * param ItemOffset
 *     Offset to the point.
 *
 * return
 *     TRUE if successful, FALSE otherwise.
 */
Boolean_t BundlesRemoveAtOffset(Bundles_pa Bundles,
                                LgIndex_t  PointOffset)
{
    Boolean_t IsOk = TRUE;

    REQUIRE(BundlesIsValid(Bundles));
    REQUIRE(0 <= PointOffset && PointOffset <= BundlesGetCount(Bundles) - 1);

    /* Remove the items for the array lists */
    ArrListRemoveItem(Bundles->Atom, PointOffset);
    ArrListRemoveItem(Bundles->Bond, PointOffset);
    ArrListRemoveItem(Bundles->Ring, PointOffset);
    ArrListRemoveItem(Bundles->Cage, PointOffset);

    IsOk = BundlesIsValid(Bundles);

    ENSURE(BundlesIsValid(Bundles));
    ENSURE(VALID_BOOLEAN(IsOk));
    return IsOk;
}



/**
 * Allocates a Bundles handle with a suitable default
 * capacity for the contained ArrList's.
 *
 * return
 *     Bundles handle if sufficient memory was available,
 *     otherewise a handle to NULL.
 */
Bundles_pa BundlesAlloc()
{
    Bundles_pa Result = NULL;

    Result = ALLOC_ITEM(Bundles_s, "Bundles structure");
    if (Result != NULL)
    {
        Result->MemUsageReported = 0;
        Result->Atom    = ArrListAlloc(60, ArrListType_Long);
        Result->Bond    = ArrListAlloc(60, ArrListType_Long);
        Result->Ring    = ArrListAlloc(60, ArrListType_Long);
        Result->Cage    = ArrListAlloc(60, ArrListType_Long);

        /* If it failed to allocate any of the array lists, clean-up and exit. */
        if (Result->Atom == NULL || Result->Bond == NULL || Result->Ring == NULL ||
            Result->Cage  == NULL)
        {
            if (Result->Atom != NULL) ArrListDealloc(&(Result->Atom));
            if (Result->Bond != NULL) ArrListDealloc(&(Result->Bond));
            if (Result->Ring != NULL) ArrListDealloc(&(Result->Ring));
            if (Result->Cage != NULL) ArrListDealloc(&(Result->Cage));
            FREE_ITEM(Result, "Bundles structure");
            Result = NULL;
        }
    }

    ENSURE(BundlesIsValid(Result) || Result == NULL);
    return Result;
}




/**
 * Gets the number of bundles currently in the Bundles structure.
 *
 * param
 *     Bundles structure in question.
 *
 * return
 *     Number of bundles in the Bundles structure.
 */
LgIndex_t BundlesGetCount(Bundles_pa Bundles)
{
    LgIndex_t Result = 0;

    REQUIRE(BundlesIsValid(Bundles));

    Result = ArrListGetCount(Bundles->Atom);

    ENSURE(Result >= 0);
    return Result;
}






/**
 * Places bundle components at the specified offset. If the
 * offset is beyond the end of the array lists, they are sized
 * accordingly and the intervening array values between the last
 * node of the original state and the last node of the new state are
 * assigned 0.
 *
 * note
 *     If components already exists at the specified location they are
 *     replaced.
 *
 * param Bundles
 *     Bundles target in which to set the coordinates.
 * param PointOffset
 *     Offset of the point/node.
 * param Atom, Bond, Ring, Cage
 *     Components to set at the specified offset.
 *
 * return
 *     TRUE if sufficient memory permitted the operation, otherwise FALSE.
 */

Boolean_t BundlesSetAtOffset(Bundles_pa Bundles,
                             LgIndex_t  PointOffset,
                             LgIndex_t  Atom,
                             LgIndex_t  Bond,
                             LgIndex_t  Ring,
                             LgIndex_t  Cage)
{
    Boolean_t IsOk = TRUE;
    ArrListItem_u Item;


    REQUIRE(BundlesIsValid(Bundles));
    REQUIRE(PointOffset >= 0);

    Item.Long = Atom;
    IsOk = ArrListSetItem(Bundles->Atom, PointOffset, Item);

    if (IsOk)
    {
        Item.Long = Bond;
        IsOk = ArrListSetItem(Bundles->Bond, PointOffset, Item);
    }

    if (IsOk)
    {
        Item.Long = Ring;
        IsOk = ArrListSetItem(Bundles->Ring, PointOffset, Item);
    }

    if (IsOk)
    {
        Item.Long = Cage;
        IsOk = ArrListSetItem(Bundles->Cage, PointOffset, Item);
    }

    ENSURE(BundlesIsValid(Bundles));

    ENSURE(VALID_BOOLEAN(IsOk));
    return IsOk;
}





/**
 * Inserts critical point components at the specified offset. The arrays
 * will be expanded to accomodate the additional value.
 *
 *
 * param Bundles
 *     Bundles target in which to set the coordinates.
 * param PointOffset
 *     Offset at which to insert the point/node.
 * param X,Y,Z
 *     Coordinates to set at the specified offset.
 * param ChrgDens, Type
 *     Charge density and critical point type to set at the specified
 *     offset.
 * param PrincDirX, ...Y, ...Z
 *     Principle direction components to set at the specified offset.
 *
 * return
 *     TRUE if sufficient memory permitted the operation, otherwise FALSE.
 */

Boolean_t BundlesInsertAtOffset(Bundles_pa Bundles,
                                LgIndex_t  PointOffset,
                                LgIndex_t  Atom,
                                LgIndex_t  Bond,
                                LgIndex_t  Ring,
                                LgIndex_t  Cage)
{
    Boolean_t IsOk = TRUE;
    ArrListItem_u Item;


    REQUIRE(BundlesIsValid(Bundles));
    REQUIRE(0 <= PointOffset && PointOffset <= BundlesGetCount(Bundles));

    Item.Long = Atom;
    IsOk = ArrListInsertItem(Bundles->Atom, PointOffset, Item);

    if (IsOk)
    {
        Item.Long = Bond;
        IsOk = ArrListInsertItem(Bundles->Bond, PointOffset, Item);
    }

    if (IsOk)
    {
        Item.Long = Ring;
        IsOk = ArrListInsertItem(Bundles->Ring, PointOffset, Item);
    }
    if (IsOk)
    {
        Item.Long = Cage;
        IsOk = ArrListInsertItem(Bundles->Cage, PointOffset, Item);
    }

    ENSURE(BundlesIsValid(Bundles));

    ENSURE(VALID_BOOLEAN(IsOk));
    return IsOk;
}





/**
 * Appends the bundle components to Bundles structure. The array lists
 * will be expanded to accommodate the additional items.
 *
 * param Bundles
 *     Bundles target to which the bundle components are to be appended.
 * param Atom, Bond, Ring, Cage
 *     Components of the bundle to append to the Bundles structure.
 *
 * return
 *     TRUE if sufficient memory permitted the operation, otherwise FALSE.
 */
Boolean_t BundlesAppendAtEnd(Bundles_pa  Bundles,
                             LgIndex_t   Atom,
                             LgIndex_t   Bond,
                             LgIndex_t   Ring,
                             LgIndex_t   Cage)
{
    Boolean_t IsOk = TRUE;
    LgIndex_t Count;

    REQUIRE(BundlesIsValid(Bundles));

    Count = BundlesGetCount(Bundles);

    IsOk = BundlesInsertAtOffset(Bundles, Count, Atom, Bond, Ring, Cage);

    ENSURE(BundlesIsValid(Bundles));

    ENSURE(VALID_BOOLEAN(IsOk));
    return IsOk;
}






/**
 * Gets the components of the Bundles structure for the
 * bundle at the specified offset.
 *
 * param Bundles
 *     Bundles structure containing the desired item.
 * param PointOffset
 *     Offset to the components of the Bundles.
 * param *Atom, *Bond, *Ring, *Cage
 *     Pointers to components of the bundle
 *
 * return
 *     TRUE if it works, FALSE otherwise.
 */

Boolean_t BundlesGetFromOffset(const Bundles_pa  Bundles,
                               const LgIndex_t   PointOffset,
                               LgIndex_t  *Atom,
                               LgIndex_t  *Bond,
                               LgIndex_t  *Ring,
                               LgIndex_t  *Cage)
{
    Boolean_t IsOk = TRUE;
    ArrListItem_u Item;

    REQUIRE(BundlesIsValid(Bundles));
    REQUIRE(PointOffset >= 0 && PointOffset < BundlesGetCount(Bundles));
    REQUIRE(VALID_REF(Atom));
    REQUIRE(VALID_REF(Bond));
    REQUIRE(VALID_REF(Ring));
    REQUIRE(VALID_REF(Cage));

    Item = ArrListGetItem(Bundles->Atom, PointOffset);
    *Atom = Item.Long;

    Item = ArrListGetItem(Bundles->Bond, PointOffset);
    *Bond = Item.Long;

    Item = ArrListGetItem(Bundles->Ring, PointOffset);
    *Ring = Item.Long;

    Item = ArrListGetItem(Bundles->Cage, PointOffset);
    *Cage = Item.Long;

    ENSURE(VALID_REF(Atom));
    ENSURE(VALID_REF(Bond));
    ENSURE(VALID_REF(Ring));
    ENSURE(VALID_REF(Cage));

    ENSURE(VALID_BOOLEAN(IsOk));
    return IsOk;
}







/**
 * Gets all irreducible bundles with the specified Atom and Cage.
 *
 * param Bundles
 *     Bundles structure containing the desired item.
 * params Atom, Cage
 *     Atom and Cage number used to search Bundles.
 * param BundlesAtomCage
 *     Bundles structure containing the subset of all Bundles that
 *     contain the specified Atom and Cage.
 *
 *     NOTE: BundlesAtomCage must be valid (allocated) before calling
 *       this function. Any information contained in BundlesAtomCage
 *       when this function is called is overwritten.
 *
 * return
 *     TRUE if it works, FALSE otherwise.
 */

Boolean_t BundlesGetBondRingFromAtomCage(const Bundles_pa    Bundles,
                                         const LgIndex_t     Atom,
                                         const LgIndex_t     Cage,
                                         Bundles_pa          BundlesAtomCage)
{
    Boolean_t IsOk = TRUE;
    LgIndex_t NumBundles = 0;
    LgIndex_t b;

    REQUIRE(BundlesIsValid(Bundles));
    REQUIRE(Atom >= 0);
    REQUIRE(Cage >= 0);
    REQUIRE(BundlesIsValid(BundlesAtomCage));

    BundlesClear(BundlesAtomCage);

    // Search through bundles, adding those that match criteria to BundlesAtomCage
    NumBundles = BundlesGetCount(Bundles);
    for (b = 0; IsOk && b < NumBundles; b++)
    {
        LgIndex_t AtomTry, Bond, Ring, CageTry;
        IsOk = BundlesGetFromOffset(Bundles, b, &AtomTry, &Bond, &Ring, &CageTry);
        if (IsOk && AtomTry == Atom && CageTry == Cage)
            IsOk = BundlesAppendAtEnd(BundlesAtomCage, Atom, Bond, Ring, Cage);
    }

    ENSURE(BundlesIsValid(BundlesAtomCage));

    ENSURE(VALID_BOOLEAN(IsOk));
    return IsOk;
}

/*
 * Create two bond bundles (one for each atom) for this Bond-Ring-Cage combination.
 * Used one-based numbers for Atom, Bond, Ring, and Cage so that negative numbers
 * can denote far-field Rings and Cages.
 *
 * Parameter
 *   BaseZoneNum
 *   BondCrtPtNum
 *   CritPoints
 *   CircleGradPath
 *   GPNumBeg
 *   GPNumEnd
 *
 * Returns TRUE if it works, FALSE otherwise.
 */
Boolean_t  CreateBundlesForBond(EntIndex_t        BaseZoneNum,
                                CritPoints_pa     CritPoints,
                                CircleGradPath_pa CircleGradPath,
                                LgIndex_t         CGPNumBeg,
                                LgIndex_t         CGPNumEnd,
                                Bundles_pa        Bundles)
{
    Boolean_t   IsOk = TRUE;
    EntIndex_t  NumZones = 0;
    LgIndex_t   NumGradPaths = CGPNumEnd - CGPNumBeg + 1;
    GradPath_pa ConnectPathBeg = NULL;
    GradPath_pa ConnectPathEnd = NULL;

    // New
    LgIndex_t  FirstCrtPtNum,      LastCrtPtNum;
    char       FirstCrtPtType,     LastCrtPtType;
    // End New
    LgIndex_t  BegCrtPtNum,        EndCrtPtNum,        ThrdCrtPtNum;
    char       BegCrtPtType,       EndCrtPtType,       ThrdCrtPtType;
    LgIndex_t  BegCrtPtTypeOffset, EndCrtPtTypeOffset, ThrdCrtPtTypeOffset;
    LgIndex_t  CGPNumBegP1, CGPNumEndM1;
    LgIndex_t  CircleGradPathCount = CircleGradPathGetCount(CircleGradPath);

    FieldData_pa XVarFDPtr = NULL;
    FieldData_pa YVarFDPtr = NULL;
    FieldData_pa ZVarFDPtr = NULL;
    FieldData_pa CDVarFDPtr = NULL;
    FieldData_pa TypeVarFDPtr = NULL;

    REQUIRE(TecUtilDataSetIsAvailable());

    /* Get num zones in dataset */
    if (IsOk)
        IsOk = TecUtilDataSetGetInfo(NULL, &NumZones, NULL);

    REQUIRE(BaseZoneNum > 0 && BaseZoneNum <= NumZones);
    REQUIRE(CritPointsIsValid(CritPoints));
    REQUIRE(CircleGradPathIsValid(CircleGradPath));
    REQUIRE(CGPNumBeg >= 0 && CGPNumBeg < CircleGradPathCount);
    REQUIRE(CGPNumEnd >= 0 && CGPNumEnd <= 2*CircleGradPathCount);

    /* Handle case where surface crosses branch cut */
    CGPNumBegP1 = CGPNumBeg + 1;
    //tmp if (CGPNumBegP1 >= CircleGradPathCount) CGPNumBegP1 -= CircleGradPathCount;

    CGPNumEndM1 = CGPNumEnd - 1;
    //tmp if (CGPNumEndM1 < 0) CGPNumEndM1 += CircleGradPathCount;

    //tmp if (NumGradPaths < 0) NumGradPaths += CircleGradPathCount;


    /* Find the three critical points in the Bond-Ring-Cage surface */
    if (IsOk)
    {
        GradPath_pa GradPath1 = CircleGradPathGetGP(CircleGradPath, CGPNumBeg);
        GradPath_pa GradPath2 = CircleGradPathGetGP(CircleGradPath, CGPNumEnd);

        BegCrtPtNum = GradPath1->BeginCrtPtNum;
        CHECK(BegCrtPtNum == GradPath2->BeginCrtPtNum);
        BegCrtPtType = CritPointsGetType(CritPoints, BegCrtPtNum);
        BegCrtPtTypeOffset = CritPointsGetBegOffset(CritPoints, BegCrtPtType);

        // TODO: In Progress: Refactoring determination of surface corner points
        FirstCrtPtNum  = GradPath1->EndCrtPtNum;
        if (FirstCrtPtNum >= 0)
            FirstCrtPtType = CritPointsGetType(CritPoints, FirstCrtPtNum);

        LastCrtPtNum   = GradPath2->EndCrtPtNum;
        if (LastCrtPtNum >= 0)
            LastCrtPtType  = CritPointsGetType(CritPoints, LastCrtPtNum);

        if (NumGradPaths > 2)
        {
            GradPath_pa GradPathMid;
            LgIndex_t   CGPNumMid = (CGPNumBeg + CGPNumEnd) / 2;

            /* Handle group that crosses branch line */
            /* tmp
            if (CGPNumEnd < CGPNumBeg)
              {
                CGPNumMid = (CGPNumBeg + CGPNumEnd + CircleGradPathCount)/2;
                if (CGPNumMid >= CircleGradPathCount) CGPNumMid -= CircleGradPathCount;
              } */

            GradPathMid = CircleGradPathGetGP(CircleGradPath, CGPNumMid);
            EndCrtPtNum = GradPathMid->EndCrtPtNum;
            if (EndCrtPtNum >= 0)
            {
                EndCrtPtType = CritPointsGetType(CritPoints, EndCrtPtNum);
                EndCrtPtTypeOffset = CritPointsGetBegOffset(CritPoints, EndCrtPtType);

                /* Determine the third critical point, the one in the middle of the dog-leg line */
                if (FirstCrtPtNum == EndCrtPtNum)
                {
                    ThrdCrtPtNum  = LastCrtPtNum;
                    ThrdCrtPtType = LastCrtPtType;
                    ThrdCrtPtTypeOffset = CritPointsGetBegOffset(CritPoints, ThrdCrtPtType);
                }
                else if (LastCrtPtNum == EndCrtPtNum)
                {
                    ThrdCrtPtNum  = FirstCrtPtNum;
                    ThrdCrtPtType = FirstCrtPtType;
                    ThrdCrtPtTypeOffset = CritPointsGetBegOffset(CritPoints, ThrdCrtPtType);
                }
                else
                    IsOk = FALSE;

            }
            else
            {
                EndCrtPtType = 0;
                EndCrtPtTypeOffset = -1;
            }
        }
        else
        {
            if (FirstCrtPtType == -3 || FirstCrtPtType == 3)
                EndCrtPtNum = FirstCrtPtNum;
            else
                EndCrtPtNum = LastCrtPtNum;
            EndCrtPtType = CritPointsGetType(CritPoints, EndCrtPtNum);
            EndCrtPtTypeOffset = CritPointsGetBegOffset(CritPoints, EndCrtPtType);
        }
    }

    /* find the Atom Crititcal Point numbers and append bundle to Bundles */
    if (IsOk)
    {
        LgIndex_t ii;
        LgIndex_t Count = 0;
        LgIndex_t Atom, Bond, Ring, Cage;
        LgIndex_t NumAtoms = CritPointsGetEndOffset(CritPoints, -3);

        CHECK(BegCrtPtType == -1);
        Bond = BegCrtPtNum - BegCrtPtTypeOffset;
        if (CritPointsGetType(CritPoints, FirstCrtPtNum) == 1 ||
            CritPointsGetType(CritPoints, FirstCrtPtNum) == 11)
        {
            if (CritPointsGetType(CritPoints, FirstCrtPtNum) == 1)
                Ring = FirstCrtPtNum - CritPointsGetBegOffset(CritPoints, 1) + 1;
            else  // Far-Fied Ring
                Ring = -1 - (FirstCrtPtNum - CritPointsGetBegOffset(CritPoints, 11));

            if (CritPointsGetType(CritPoints, LastCrtPtNum) == 3)
                Cage = LastCrtPtNum - CritPointsGetBegOffset(CritPoints, 3) + 1;
            else  // Far-Field Cage
                Cage = -1 - (LastCrtPtNum - CritPointsGetBegOffset(CritPoints, 13));
        }
        else  // FirstCrtPtNum == Cage
        {
            if (CritPointsGetType(CritPoints, FirstCrtPtNum) == 3)
                Cage = FirstCrtPtNum - CritPointsGetBegOffset(CritPoints, 3) + 1;
            else  // Far-Field Cage
                Cage = -1 - (FirstCrtPtNum - CritPointsGetBegOffset(CritPoints, 13));

            if (CritPointsGetType(CritPoints, LastCrtPtNum) == 1)
                Ring = LastCrtPtNum - CritPointsGetBegOffset(CritPoints, 1) + 1;
            else  // Far-Fied Ring
                Ring = -1 - (LastCrtPtNum - CritPointsGetBegOffset(CritPoints, 11));
        }

        for (ii = 0; IsOk && Count < 2 && ii < NumAtoms; ii++)
        {
            GradPath_pa BondAtomPath = GradPathGetByBegEndCP(BegCrtPtNum, ii);
            if (BondAtomPath != NULL)
            {
                Atom = BondAtomPath->EndCrtPtNum;

                // Save one-based numbers for Atom, Bond, Ring, and Cage
                Atom++;
                Bond++;

                IsOk = BundlesAppendAtEnd(Bundles, Atom, Bond, Ring, Cage);
                Count++;
            }
        }
    }

    REQUIRE(BundlesIsValid(Bundles));
    REQUIRE(VALID_BOOLEAN(IsOk));
    return IsOk;
}

								   /*
 * Create list of Cages (ArrList) that connect directly to the specified Atom.
 *
 * Parameter
 *   Atom
 *   Bundles
 *   CritPoints
 *
 * Returns AtomCageList if it works, NULL otherwise.
 */
ArrList_pa  CreateAtomCageList(LgIndex_t        Atom,
                               Bundles_pa       Bundles,
                               CritPoints_pa    CritPoints)
{
    ArrList_pa AtomCageList = NULL;
    Boolean_t  IsOk = TRUE;
    int ii;

    REQUIRE(BundlesIsValid(Bundles));
    REQUIRE(CritPointsIsValid(CritPoints));
    REQUIRE(Atom >= 0 && Atom < CritPointsGetEndOffset(CritPoints, -3));

    AtomCageList = ArrListAlloc(20, ArrListType_Long);

    if (AtomCageList != NULL)
    {
        /* Search through Bundles looking for this Atom */
        for (ii = 0; ii < BundlesGetCount(Bundles); ii++)
        {
            LgIndex_t AtomInBundle, Bond, Ring, Cage;
            IsOk = BundlesGetFromOffset(Bundles, ii, &AtomInBundle, &Bond, &Ring, &Cage);

            if (IsOk && AtomInBundle == Atom)
            {
                LgIndex_t NumCages = ArrListGetCount(AtomCageList);

                /* Insert cage (if unique) at appropriate location (increasing value) */
                if (NumCages > 0)
                {
                    int jj;
                    Boolean_t Found = FALSE;
                    Boolean_t Duplicate = FALSE;

                    for (jj = 0; !Found && !Duplicate && jj < NumCages; jj++)
                    {
                        ArrListItem_u Item;
                        LgIndex_t     CageInList;

                        Item = ArrListGetItem(AtomCageList, jj);
                        CageInList = Item.Long;

                        if (Cage == CageInList) Duplicate = TRUE;
                        if (Cage <= CageInList) Found = TRUE;
                    }

                    /* Insert at appropriate location. Ignore duplicates */
                    if (!Duplicate)
                    {
                        ArrListItem_u Item;
                        Item.Long = Cage;
                        if (Found)
                            IsOk = ArrListInsertItem(AtomCageList, jj - 1, Item);
                        else
                            IsOk = ArrListAppendItem(AtomCageList, Item);
                    }
                }
                else // First Item, so append
                {
                    ArrListItem_u Item;
                    Item.Long = Cage;
                    IsOk = ArrListAppendItem(AtomCageList, Item);
                }
            }
        }
    }

    if (!IsOk) ArrListDealloc(&AtomCageList);

    ENSURE(AtomCageList == NULL || ArrListIsValid(AtomCageList));

    return AtomCageList;
}

/*
 * Create list of Rings and Bonds (ArrList) that connect directly to the
 * specified Atom and Cage.
 *
 * Parameter
 *   Atom
 *   Cage
 *   Bundles
 *   CritPoints
 *
 * Returns AtomCageList if it works, NULL otherwise.
 */
ArrList_pa  CreateRingBondListForAtomCage(LgIndex_t        Atom,
                                          LgIndex_t        Cage,
                                          Bundles_pa       Bundles,
                                          CritPoints_pa    CritPoints)
{
    ArrList_pa BondRingList = NULL;
    Boolean_t  IsOk = TRUE;
    int ii;

    REQUIRE(BundlesIsValid(Bundles));
    REQUIRE(CritPointsIsValid(CritPoints));
    REQUIRE(Atom >= 0 && Atom < CritPointsGetEndOffset(CritPoints, -3));
    REQUIRE(Cage >= CritPointsGetBegOffset(CritPoints, 3) && Cage < CritPointsGetEndOffset(CritPoints, 3));

    BondRingList = ArrListAlloc(20, ArrListType_Long);

    if (BondRingList != NULL)
    {
        /* Search through Bundles looking for this Atom/Cage pair */
        for (ii = 0; ii < BundlesGetCount(Bundles); ii++)
        {
            LgIndex_t AtomInBundle, Bond, Ring, CageInBundle;
            IsOk = BundlesGetFromOffset(Bundles, ii, &AtomInBundle, &Bond, &Ring, &CageInBundle);

            if (IsOk && AtomInBundle == Atom && CageInBundle == Cage)
            {
                LgIndex_t NumBondRings = ArrListGetCount(BondRingList);

                /* Insert cage (if unique) at appropriate location (increasing value) */
                if (NumBondRings > 0)
                {
                    int jj;
                    Boolean_t BondDuplicate = FALSE;
                    Boolean_t RingDuplicate = FALSE;

                    for (jj = 0; jj < NumBondRings; jj++)
                    {
                        ArrListItem_u Item;
                        LgIndex_t     BondRingInList;

                        Item = ArrListGetItem(BondRingList, jj);
                        BondRingInList = Item.Long;

                        if (Bond == BondRingInList) BondDuplicate = TRUE;
                        if (Ring == BondRingInList) RingDuplicate = TRUE;
                    }

                    /* Append Bond. Ignore duplicates */
                    if (!BondDuplicate)
                    {
                        ArrListItem_u Item;
                        Item.Long = Bond;
                        IsOk = ArrListAppendItem(BondRingList, Item);
                    }
                    /* Append Ring. Ignore duplicates */
                    if (!RingDuplicate)
                    {
                        ArrListItem_u Item;
                        Item.Long = Ring;
                        IsOk = ArrListAppendItem(BondRingList, Item);
                    }
                }
                else // First Item, so append
                {
                    ArrListItem_u Item;
                    Item.Long = Bond;
                    IsOk = ArrListAppendItem(BondRingList, Item);
                    Item.Long = Ring;
                    IsOk = ArrListAppendItem(BondRingList, Item);
                }
            }
        }
    }

    if (!IsOk) ArrListDealloc(&BondRingList);

    ENSURE(BondRingList == NULL || ArrListIsValid(BondRingList));

    return BondRingList;
}

/*
 * Create list of Cages (ArrList) that connect directly to the specified Atom and
 * specified Bond.
 *
 * Parameter
 *   Atom
 *   Bundles
 *   CritPoints
 *
 * Returns AtomBondCageList if it works, NULL otherwise.
 */
ArrList_pa  CreateAtomBondCageList(LgIndex_t        Atom,
								   LgIndex_t        Bond,
                                   Bundles_pa       Bundles,
                                   CritPoints_pa    CritPoints)
{
    ArrList_pa AtomBondCageList = NULL;
    Boolean_t  IsOk = TRUE;
    int ii;

    REQUIRE(BundlesIsValid(Bundles));
    REQUIRE(CritPointsIsValid(CritPoints));
    REQUIRE(Atom > 0 && Atom <= CritPointsGetEndOffset(CritPoints, -3));
    REQUIRE(Bond > 0 && Bond <= CritPointsGetEndOffset(CritPoints, -1) - CritPointsGetEndOffset(CritPoints, -3) );

    AtomBondCageList = ArrListAlloc(20, ArrListType_Long);

    if (AtomBondCageList != NULL)
    {
        /* Search through Bundles looking for this Atom */
        for (ii = 0; ii < BundlesGetCount(Bundles); ii++)
        {
            LgIndex_t AtomInBundle, BondInBundle, Ring, Cage;
            IsOk = BundlesGetFromOffset(Bundles, ii, &AtomInBundle, &BondInBundle, &Ring, &Cage);

            if (IsOk && AtomInBundle == Atom && BondInBundle == Bond)
            {
                LgIndex_t NumCages = ArrListGetCount(AtomBondCageList);

                /* Insert cage (if unique) at appropriate location (increasing value) */
                if (NumCages > 0)
                {
                    int jj;
                    Boolean_t Found = FALSE;
                    Boolean_t Duplicate = FALSE;

                    for (jj = 0; !Found && !Duplicate && jj < NumCages; jj++)
                    {
                        ArrListItem_u Item;
                        LgIndex_t     CageInList;

                        Item = ArrListGetItem(AtomBondCageList, jj);
                        CageInList = Item.Long;

                        if (Cage == CageInList) Duplicate = TRUE;
                        if (Cage <= CageInList) Found = TRUE;
                    }

                    /* Insert at appropriate location. Ignore duplicates */
                    if (!Duplicate)
                    {
                        ArrListItem_u Item;
                        Item.Long = Cage;
                        if (Found)
                            IsOk = ArrListInsertItem(AtomBondCageList, jj - 1, Item);
                        else
                            IsOk = ArrListAppendItem(AtomBondCageList, Item);
                    }
                }
                else // First Item, so append
                {
                    ArrListItem_u Item;
                    Item.Long = Cage;
                    IsOk = ArrListAppendItem(AtomBondCageList, Item);
                }
            }
        }
    }

    if (!IsOk) ArrListDealloc(&AtomBondCageList);

    ENSURE(AtomBondCageList == NULL || ArrListIsValid(AtomBondCageList));

    return AtomBondCageList;
}


