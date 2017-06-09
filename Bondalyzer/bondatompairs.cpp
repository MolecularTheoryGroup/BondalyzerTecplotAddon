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
#include "BONDATOMPAIRS.h"
#include <string.h>



/**
 * Determine if the BondAtomPairs handle is sane.
 *
 * param BondAtomPairs
 *     BondAtomPairs structure in question.
 *
 * return
 *     TRUE if the BondAtomPairs structure is valid, otherwise FALSE.
 */
Boolean_t BondAtomPairsIsValid(BondAtomPairs_pa BondAtomPairs)
{
    Boolean_t IsValid = FALSE;

    IsValid = (VALID_REF(BondAtomPairs) &&
               VALID_REF(BondAtomPairs->Atom1) && ArrListIsValid(BondAtomPairs->Atom1) &&
               VALID_REF(BondAtomPairs->Atom2) && ArrListIsValid(BondAtomPairs->Atom2));

    /* Require the same count for each array list in BondAtomPairs structure. */
    if (IsValid)
    {
        LgIndex_t Count = ArrListGetCount(BondAtomPairs->Atom1);
        IsValid = (ArrListGetCount(BondAtomPairs->Atom2) == Count);
    }

    ENSURE(VALID_BOOLEAN(IsValid));
    return IsValid;
}




/**
 * Deallocates the BondAtomPairs handle and set the handle to NULL.
 *
 * note
 *     Item destruction is the responsibility of the caller.
 *
 * param
 *     Reference to a BondAtomPairs handle.
 */
void BondAtomPairsDealloc(BondAtomPairs_pa *BondAtomPairs)
{
    REQUIRE(VALID_REF(BondAtomPairs));
    REQUIRE(BondAtomPairsIsValid(*BondAtomPairs) || *BondAtomPairs == NULL);

    if (*BondAtomPairs != NULL)
    {
        /* release the ArrList's */
        if ((*BondAtomPairs)->Atom1 != NULL) ArrListDealloc(&((*BondAtomPairs)->Atom1));
        if ((*BondAtomPairs)->Atom2 != NULL) ArrListDealloc(&((*BondAtomPairs)->Atom2));

        /* Notify Tecplot of memory usage change */
        TecUtilMemoryChangeNotify(- (Int64_t)((*BondAtomPairs)->MemUsageReported));
        (*BondAtomPairs)->MemUsageReported = 0;

        /* release the list structure itself */
        FREE_ITEM(*BondAtomPairs, "BondAtomPairs structure");
        *BondAtomPairs = NULL;
    }

    ENSURE(*BondAtomPairs == NULL);
}





/**
 * Empties the BondAtomPairs structure of all Bond Atom Pairs.
 *
 * param BondAtomPairs
 *     BondAtomPairs to clear.
 */
void BondAtomPairsClear(BondAtomPairs_pa BondAtomPairs)
{
    REQUIRE(BondAtomPairsIsValid(BondAtomPairs));

    ArrListClear(BondAtomPairs->Atom1);
    ArrListClear(BondAtomPairs->Atom2);

    ENSURE(BondAtomPairsIsValid(BondAtomPairs) && BondAtomPairsGetCount(BondAtomPairs) == 0);
}





/**
 * Removes a point from the BondAtomPairs array. The members following the item
 * removed are shifted down accordingly to fill the vacated space.
 *
 *
 * param BondAtomPairs
 *     BondAtomPairs structure containing the bundle to remove.
 * param ItemOffset
 *     Offset to the point.
 *
 * return
 *     TRUE if successful, FALSE otherwise.
 */
Boolean_t BondAtomPairsRemoveAtOffset(BondAtomPairs_pa BondAtomPairs,
                                      LgIndex_t  PointOffset)
{
    Boolean_t IsOk = TRUE;

    REQUIRE(BondAtomPairsIsValid(BondAtomPairs));
    REQUIRE(0 <= PointOffset && PointOffset <= BondAtomPairsGetCount(BondAtomPairs) - 1);

    /* Remove the items for the array lists */
    ArrListRemoveItem(BondAtomPairs->Atom1, PointOffset);
    ArrListRemoveItem(BondAtomPairs->Atom2, PointOffset);

    IsOk = BondAtomPairsIsValid(BondAtomPairs);

    ENSURE(BondAtomPairsIsValid(BondAtomPairs));
    ENSURE(VALID_BOOLEAN(IsOk));
    return IsOk;
}



/**
 * Allocates a BondAtomPairs handle with a suitable default
 * capacity for the contained ArrList's.
 *
 * return
 *     BondAtomPairs handle if sufficient memory was available,
 *     otherewise a handle to NULL.
 */
BondAtomPairs_pa BondAtomPairsAlloc()
{
    BondAtomPairs_pa Result = NULL;

    Result = ALLOC_ITEM(BondAtomPairs_s, "BondAtomPairs structure");
    if (Result != NULL)
    {
        Result->MemUsageReported = 0;
        Result->Atom1   = ArrListAlloc(60, ArrListType_Long);
        Result->Atom2   = ArrListAlloc(60, ArrListType_Long);

        /* If it failed to allocate any of the array lists, clean-up and exit. */
        if (Result->Atom1 == NULL || Result->Atom2 == NULL)
        {
            if (Result->Atom1 != NULL) ArrListDealloc(&(Result->Atom1));
            if (Result->Atom2 != NULL) ArrListDealloc(&(Result->Atom1));
            FREE_ITEM(Result, "BondAtomPairs structure");
            Result = NULL;
        }
    }

    ENSURE(BondAtomPairsIsValid(Result) || Result == NULL);
    return Result;
}




/**
 * Gets the number of BondAtomPairs currently in the BondAtomPairs structure.
 *
 * param
 *     BondAtomPairs structure in question.
 *
 * return
 *     Number of bond atom pairs in the BondAtomPairs structure.
 */
LgIndex_t BondAtomPairsGetCount(BondAtomPairs_pa BondAtomPairs)
{
    LgIndex_t Result = 0;

    REQUIRE(BondAtomPairsIsValid(BondAtomPairs));

    Result = ArrListGetCount(BondAtomPairs->Atom1);

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
 * param BondAtomPairs
 *     BondAtomPairs target in which to set the coordinates.
 * param PointOffset
 *     Offset of the point/node.
 * param Atom, Bond, Ring, Cage
 *     Components to set at the specified offset.
 *
 * return
 *     TRUE if sufficient memory permitted the operation, otherwise FALSE.
 */

Boolean_t BondAtomPairsSetAtOffset(BondAtomPairs_pa BondAtomPairs,
                                   LgIndex_t  PointOffset,
                                   LgIndex_t   Atom1,
                                   LgIndex_t   Atom2)
{
    Boolean_t IsOk = TRUE;
    ArrListItem_u Item;


    REQUIRE(BondAtomPairsIsValid(BondAtomPairs));
    REQUIRE(PointOffset >= 0);

    Item.Long = Atom1;
    IsOk = ArrListSetItem(BondAtomPairs->Atom1, PointOffset, Item);

    if (IsOk)
    {
        Item.Long = Atom2;
        IsOk = ArrListSetItem(BondAtomPairs->Atom2, PointOffset, Item);
    }

    ENSURE(BondAtomPairsIsValid(BondAtomPairs));

    ENSURE(VALID_BOOLEAN(IsOk));
    return IsOk;
}





/**
 * Inserts critical point components at the specified offset. The arrays
 * will be expanded to accomodate the additional value.
 *
 *
 * param BondAtomPairs
 *     BondAtomPairs target in which to set the coordinates.
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

Boolean_t BondAtomPairsInsertAtOffset(BondAtomPairs_pa BondAtomPairs,
                                      LgIndex_t   PointOffset,
                                      LgIndex_t   Atom1,
                                      LgIndex_t   Atom2)
{
    Boolean_t IsOk = TRUE;
    ArrListItem_u Item;


    REQUIRE(BondAtomPairsIsValid(BondAtomPairs));
    REQUIRE(0 <= PointOffset && PointOffset <= BondAtomPairsGetCount(BondAtomPairs));

    Item.Long = Atom1;
    IsOk = ArrListInsertItem(BondAtomPairs->Atom1, PointOffset, Item);

    if (IsOk)
    {
        Item.Long = Atom2;
        IsOk = ArrListInsertItem(BondAtomPairs->Atom2, PointOffset, Item);
    }

    ENSURE(BondAtomPairsIsValid(BondAtomPairs));

    ENSURE(VALID_BOOLEAN(IsOk));
    return IsOk;
}





/**
 * Appends the bundle components to BondAtomPairs structure. The array lists
 * will be expanded to accommodate the additional items.
 *
 * param BondAtomPairs
 *     BondAtomPairs target to which the bundle components are to be appended.
 * param Atom, Bond, Ring, Cage
 *     Components of the bundle to append to the BondAtomPairs structure.
 *
 * return
 *     TRUE if sufficient memory permitted the operation, otherwise FALSE.
 */
Boolean_t BondAtomPairsAppendAtEnd(BondAtomPairs_pa  BondAtomPairs,
                                   LgIndex_t   Atom1,
                                   LgIndex_t   Atom2)
{
    Boolean_t IsOk = TRUE;
    LgIndex_t Count;

    REQUIRE(BondAtomPairsIsValid(BondAtomPairs));

    Count = BondAtomPairsGetCount(BondAtomPairs);

    IsOk = BondAtomPairsInsertAtOffset(BondAtomPairs, Count, Atom1, Atom2);

    ENSURE(BondAtomPairsIsValid(BondAtomPairs));

    ENSURE(VALID_BOOLEAN(IsOk));
    return IsOk;
}






/**
 * Gets the components of the BondAtomPairs structure for the
 * bundle at the specified offset.
 *
 * param BondAtomPairs
 *     BondAtomPairs structure containing the desired item.
 * param PointOffset
 *     Offset to the components of the BondAtomPairs.
 * param *Atom, *Bond, *Ring, *Cage
 *     Pointers to components of the bundle
 *
 * return
 *     TRUE if it works, FALSE otherwise.
 */

Boolean_t BondAtomPairsGetFromOffset(BondAtomPairs_pa  BondAtomPairs,
                                     LgIndex_t   PointOffset,
                                     LgIndex_t  *Atom1,
                                     LgIndex_t  *Atom2)
{
    Boolean_t IsOk = TRUE;
    ArrListItem_u Item;

    REQUIRE(BondAtomPairsIsValid(BondAtomPairs));
    REQUIRE(PointOffset >= 0 && PointOffset < BondAtomPairsGetCount(BondAtomPairs));
    REQUIRE(VALID_REF(Atom1));
    REQUIRE(VALID_REF(Atom2));

    Item = ArrListGetItem(BondAtomPairs->Atom1, PointOffset);
    *Atom1 = Item.Long;

    Item = ArrListGetItem(BondAtomPairs->Atom2, PointOffset);
    *Atom2 = Item.Long;

    ENSURE(VALID_REF(Atom1));
    ENSURE(VALID_REF(Atom2));

    ENSURE(VALID_BOOLEAN(IsOk));
    return IsOk;
}





