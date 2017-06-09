#include "TECADDON.h"
#include "ADDGLBL.h"
#include "ADDONVER.h"
#include "TASSERT.h"
#include "GUIDEFS.h"
#include "ADKUTIL.h"
#include "ALLOC.h"
#include "ARRLIST.h"
#include "BUNDLES.h"
#include "CRITPOINTS.h"
#include "ENGINE.h"
#include <string.h>



/**
 * Determine if the Bundles handle is sane.
 *
 * param Bundles
 *     Bundles structure in question.
 *
 * return
 *     TRUE if the Bundles structure is valid, otherwise FALSE.
 */
Boolean_t BundlesIsValid(const Bundles_pa Bundles)
{
    Boolean_t IsValid = FALSE;

    IsValid = (VALID_REF(Bundles) &&
               VALID_REF(Bundles->Atom) && ArrListIsValid(Bundles->Atom) &&
               VALID_REF(Bundles->Bond) && ArrListIsValid(Bundles->Bond) &&
               VALID_REF(Bundles->Ring) && ArrListIsValid(Bundles->Ring) &&
               VALID_REF(Bundles->CageOrRing2) && ArrListIsValid(Bundles->CageOrRing2) &&
               VALID_REF(Bundles->CageOrRingCPType) && ArrListIsValid(Bundles->CageOrRingCPType) &&
               VALID_REF(Bundles->BundleType) && ArrListIsValid(Bundles->BundleType));

    /* Require the same count for each array list in Bundles structure. */
    if (IsValid)
    {
        LgIndex_t Count = ArrListGetCount(Bundles->Atom);
        IsValid = (ArrListGetCount(Bundles->Ring) == Count);
        if (IsValid) IsValid = (ArrListGetCount(Bundles->Bond) == Count);
        if (IsValid) IsValid = (ArrListGetCount(Bundles->CageOrRing2) == Count);
        if (IsValid) IsValid = (ArrListGetCount(Bundles->CageOrRingCPType) == Count);
        if (IsValid) IsValid = (ArrListGetCount(Bundles->BundleType) == Count);
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
        if ((*Bundles)->CageOrRing2 != NULL) ArrListDealloc(&((*Bundles)->CageOrRing2));
        if ((*Bundles)->CageOrRingCPType != NULL) ArrListDealloc(&((*Bundles)->CageOrRingCPType));
        if ((*Bundles)->BundleType != NULL) ArrListDealloc(&((*Bundles)->BundleType));

        /* Notify Tecplot of memory usage change */
        TecUtilMemoryChangeNotify(- (Int64_t)((*Bundles)->MemUsageReported));
        (*Bundles)->MemUsageReported = 0;

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
    ArrListClear(Bundles->CageOrRing2);
    ArrListClear(Bundles->CageOrRingCPType);
    ArrListClear(Bundles->BundleType);

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
    ArrListRemoveItem(Bundles->CageOrRing2, PointOffset);
    ArrListRemoveItem(Bundles->CageOrRingCPType, PointOffset);
    ArrListRemoveItem(Bundles->BundleType, PointOffset);

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
        Result->CageOrRing2 = ArrListAlloc(60, ArrListType_Long);
        Result->CageOrRingCPType  = ArrListAlloc(60, ArrListType_Char);
        Result->BundleType  = ArrListAlloc(60, ArrListType_Char);

        /* If it failed to allocate any of the array lists, clean-up and exit. */
        if (Result->Atom == NULL || Result->Bond == NULL || Result->Ring == NULL ||
            Result->CageOrRing2  == NULL || Result->CageOrRingCPType == NULL ||
            Result->BundleType == NULL)
        {
            if (Result->Atom != NULL) ArrListDealloc(&(Result->Atom));
            if (Result->Bond != NULL) ArrListDealloc(&(Result->Bond));
            if (Result->Ring != NULL) ArrListDealloc(&(Result->Ring));
            if (Result->CageOrRing2 != NULL) ArrListDealloc(&(Result->CageOrRing2));
            if (Result->CageOrRingCPType != NULL) ArrListDealloc(&(Result->CageOrRingCPType));
            if (Result->BundleType != NULL) ArrListDealloc(&(Result->BundleType));
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
                             LgIndex_t  CageOrRing2,
                             char       CageOrRingCPType,
                             char       BundleType)
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
        Item.Long = CageOrRing2;
        IsOk = ArrListSetItem(Bundles->CageOrRing2, PointOffset, Item);
    }

    if (IsOk)
    {
        Item.Char = CageOrRingCPType;
        IsOk = ArrListSetItem(Bundles->CageOrRingCPType, PointOffset, Item);
    }

    if (IsOk)
    {
        Item.Char = BundleType;
        IsOk = ArrListSetItem(Bundles->BundleType, PointOffset, Item);
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
                                LgIndex_t  CageOrRing2,
                                char       CageOrRingCPType,
                                char       BundleType)
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
        Item.Long = CageOrRing2;
        IsOk = ArrListInsertItem(Bundles->CageOrRing2, PointOffset, Item);
    }
    if (IsOk)
    {
        Item.Char = CageOrRingCPType;
        IsOk = ArrListInsertItem(Bundles->CageOrRingCPType, PointOffset, Item);
    }
    if (IsOk)
    {
        Item.Char = BundleType;
        IsOk = ArrListInsertItem(Bundles->BundleType, PointOffset, Item);
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
                             LgIndex_t   CageOrRing2,
                             char        CageOrRingCPType,
                             char        BundleType)
{
    Boolean_t IsOk = TRUE;
    LgIndex_t Count;

    REQUIRE(BundlesIsValid(Bundles));

    Count = BundlesGetCount(Bundles);

    IsOk = BundlesInsertAtOffset(Bundles, Count, Atom, Bond, Ring, CageOrRing2, CageOrRingCPType, BundleType);

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
                               LgIndex_t  *CageOrRing2,
                               char       *CageOrRingCPType,
                               char       *BundleType)
{
    Boolean_t IsOk = TRUE;
    ArrListItem_u Item;

    REQUIRE(BundlesIsValid(Bundles));
    REQUIRE(PointOffset >= 0 && PointOffset < BundlesGetCount(Bundles));
    REQUIRE(VALID_REF(Atom));
    REQUIRE(VALID_REF(Bond));
    REQUIRE(VALID_REF(Ring));
    REQUIRE(VALID_REF(CageOrRing2));
    REQUIRE(VALID_REF(CageOrRingCPType));
    REQUIRE(VALID_REF(BundleType));

    Item = ArrListGetItem(Bundles->Atom, PointOffset);
    *Atom = Item.Long;

    Item = ArrListGetItem(Bundles->Bond, PointOffset);
    *Bond = Item.Long;

    Item = ArrListGetItem(Bundles->Ring, PointOffset);
    *Ring = Item.Long;

    Item = ArrListGetItem(Bundles->CageOrRing2, PointOffset);
    *CageOrRing2 = Item.Long;

    Item = ArrListGetItem(Bundles->CageOrRingCPType, PointOffset);
    *CageOrRingCPType = Item.Char;

    Item = ArrListGetItem(Bundles->BundleType, PointOffset);
    *BundleType = Item.Char;

    ENSURE(VALID_REF(Atom));
    ENSURE(VALID_REF(Bond));
    ENSURE(VALID_REF(Ring));
    ENSURE(VALID_REF(CageOrRing2));
    ENSURE(VALID_REF(CageOrRingCPType));
    ENSURE(VALID_REF(BundleType));

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
        LgIndex_t AtomTry, Bond, Ring, CageTry, CageOrRingTry;
        char CageOrRingType, BundleType;

        IsOk = BundlesGetFromOffset(Bundles, b, &AtomTry, &Bond, &Ring, &CageOrRingTry,
                                    &CageOrRingType, &BundleType);
        if (CageOrRingType == 3)
        {
            CageTry = CageOrRingTry;

            if (IsOk && AtomTry == Atom && CageTry == Cage)
                IsOk = BundlesAppendAtEnd(BundlesAtomCage, Atom, Bond, Ring, CageTry,
                                          CageOrRingType, BundleType);
        }
    }

    ENSURE(BundlesIsValid(BundlesAtomCage));

    ENSURE(VALID_BOOLEAN(IsOk));
    return IsOk;
}








/**
 * Check to see if the critical point pair is an edge of the specified bundle.
 *
 * param Bundles
 *     Bundles structure containing the desired item.
 * params BundleNum
 *     Bundle number used to be searched.
 * param BegCrtPtType, BegCrtPtNum, EndCrtPtType, EndCrtPtNum
 *     Beginning and ending critical point types and numbers
 *
 * return
 *     TRUE if critical points match an edge of the bundle, FALSE otherwise.
 */
Boolean_t BundlesIsEdge(const Bundles_pa    Bundles,
                        const LgIndex_t     BundleNum,
                        const char          BegCrtPtType,
                        const LgIndex_t     BegCrtPtNum,
                        const char          EndCrtPtType,
                        const LgIndex_t     EndCrtPtNum)
{
    Boolean_t IsMatch = FALSE;
    Boolean_t IsOk = TRUE;

    LgIndex_t Atom, Bond, Ring, Cage, Ring2, CageOrRing;
    char CageOrRingType, BundleType;

    REQUIRE(BundlesIsValid(Bundles));
    REQUIRE(BundleNum >= 0 && BundleNum < BundlesGetCount(Bundles));
    REQUIRE(BegCrtPtType == -3 || BegCrtPtType == -1 || BegCrtPtType == 1 || BegCrtPtType == 3);
    REQUIRE(BegCrtPtNum >= 0);
    REQUIRE(EndCrtPtNum < 0 || EndCrtPtType == -3 || EndCrtPtType == -1 || EndCrtPtType == 1 || EndCrtPtType == 3);

    IsOk = BundlesGetFromOffset(Bundles, BundleNum, &Atom, &Bond, &Ring, &CageOrRing,
                                &CageOrRingType, &BundleType);
    if (CageOrRingType == 3)
    {
        Cage = CageOrRing;
        Ring2 = -1;
    }
    else if (CageOrRingType == 1)
    {
        Cage = -1;
        Ring2 = CageOrRing;
    }

    if (IsOk)
    {
        switch (BegCrtPtType)
        {
            case -3:
                IsMatch = BegCrtPtNum == Atom;
                break;
            case -1:
                IsMatch = BegCrtPtNum == Bond;
                break;
            case  1:
                IsMatch = BegCrtPtNum == Ring || BegCrtPtNum == Ring2;
                break;
            case  3:
                IsMatch = BegCrtPtNum == Cage;
                break;
            default:
                IsMatch = FALSE;
                break;
        }

        if (IsMatch)
        {
            switch (EndCrtPtType)
            {
                case -3:
                    IsMatch = EndCrtPtNum == Atom;
                    break;
                case -1:
                    IsMatch = EndCrtPtNum == Bond;
                    break;
                case  1:
                    IsMatch = EndCrtPtNum == Ring || EndCrtPtNum == Ring2;
                    break;
                case  3:
                    IsMatch = EndCrtPtNum == Cage;
                    break;
            }
        }
    }

    ENSURE(IsOk == TRUE);
    ENSURE(VALID_BOOLEAN(IsMatch));
    return IsMatch;
}










/**
 * Check to see if the critical point pair is an edge of the specified bundle.
 *
 * param Bundles
 *     Bundles structure containing the desired item.
 * params BundleNum
 *     Bundle number used to be searched.
 * param BegCrtPtType, BegCrtPtNum, EndCrtPtType, EndCrtPtNum, ThrdCrtPtType, ThrdCrtPtNum
 *     Beginning, ending, and third critical point types and numbers
 *
 * return
 *     TRUE if critical points match a surface of the bundle, FALSE otherwise.
 */
Boolean_t BundlesIsSurface(const Bundles_pa    Bundles,
                           const LgIndex_t     BundleNum,
                           const char          BegCrtPtType,
                           const LgIndex_t     BegCrtPtNum,
                           const char          EndCrtPtType,
                           const LgIndex_t     EndCrtPtNum,
                           const char          ThrdCrtPtType,
                           const LgIndex_t     ThrdCrtPtNum)
{
    Boolean_t IsMatch = FALSE;
    Boolean_t IsOk = TRUE;

    LgIndex_t Atom, Bond, Ring, Ring2, Cage, CageOrRing;
    char CageOrRingType, BundleType;

    REQUIRE(BundlesIsValid(Bundles));
    REQUIRE(BundleNum >= 0 && BundleNum < BundlesGetCount(Bundles));
    REQUIRE(BegCrtPtType == -3 || BegCrtPtType == -1 || BegCrtPtType == 1 || BegCrtPtType == 3);
    REQUIRE(BegCrtPtNum >= 0);
    REQUIRE(EndCrtPtNum < 0 || EndCrtPtType == -3 || EndCrtPtType == -1 || EndCrtPtType == 1 || EndCrtPtType == 3);
    REQUIRE(ThrdCrtPtNum < 0 || ThrdCrtPtType == -3 || ThrdCrtPtType == -1 || ThrdCrtPtType == 1 || ThrdCrtPtType == 3);

    IsOk = BundlesGetFromOffset(Bundles, BundleNum, &Atom, &Bond, &Ring, &CageOrRing,
                                &CageOrRingType, &BundleType);

    if (CageOrRingType == 3)
    {
        Ring2 = -1;
        Cage  = CageOrRing;
    }
    else if (CageOrRingType == 1)
    {
        Ring2 = CageOrRing;
        Cage = -1;
    }

    if (IsOk)
    {
        switch (BegCrtPtType)
        {
            case -3:
                IsMatch = BegCrtPtNum == Atom;
                break;
            case -1:
                IsMatch = BegCrtPtNum == Bond;
                break;
            case  1:
                IsMatch = BegCrtPtNum == Ring;
                break;
            case  3:
                IsMatch = BegCrtPtNum == Cage;
                break;
            default:
                IsMatch = FALSE;
                break;
        }

        if (IsMatch)
        {
            switch (EndCrtPtType)
            {
                case -3:
                    IsMatch = EndCrtPtNum == Atom;
                    break;
                case -1:
                    IsMatch = EndCrtPtNum == Bond;
                    break;
                case  1:
                    IsMatch = EndCrtPtNum == Ring || EndCrtPtNum == Ring2;
                    break;
                case  3:
                    IsMatch = EndCrtPtNum == Cage;
                    break;
            }
        }

        if (IsMatch)
        {
            switch (ThrdCrtPtType)
            {
                case -3:
                    IsMatch = ThrdCrtPtNum == Atom;
                    break;
                case -1:
                    IsMatch = ThrdCrtPtNum == Bond;
                    break;
                case  1:
                    IsMatch = ThrdCrtPtNum == Ring || EndCrtPtNum == Ring2;
                    break;
                case  3:
                    IsMatch = ThrdCrtPtNum == Cage;
                    break;
            }
        }
    }

    ENSURE(IsOk == TRUE);
    ENSURE(VALID_BOOLEAN(IsMatch));
    return IsMatch;
}



