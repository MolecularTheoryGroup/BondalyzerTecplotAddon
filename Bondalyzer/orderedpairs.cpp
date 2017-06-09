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
#include "ORDEREDPAIRS.h"
#include <string.h>



/**
 * Determine if the OrderedPairs handle is sane.
 *
 * param OrderedPairs
 *     OrderedPairs structure in question.
 *
 * return
 *     TRUE if the OrderedPairs structure is valid, otherwise FALSE.
 */
Boolean_t OrderedPairsIsValid(const OrderedPairs_pa OrderedPairs)
{
    Boolean_t IsValid = FALSE;

    IsValid = (VALID_REF(OrderedPairs) &&
               VALID_REF(OrderedPairs->Item1) && ArrListIsValid(OrderedPairs->Item1) &&
               VALID_REF(OrderedPairs->Item2) && ArrListIsValid(OrderedPairs->Item2));

    /* Require the same count for each array list in OrderedPairs structure. */
    if (IsValid)
    {
        LgIndex_t Count = ArrListGetCount(OrderedPairs->Item1);
        IsValid = (ArrListGetCount(OrderedPairs->Item2) == Count);
    }

    ENSURE(VALID_BOOLEAN(IsValid));
    return IsValid;
}




/**
 * Deallocates the OrderedPairs handle and set the handle to NULL.
 *
 * note
 *     Item destruction is the responsibility of the caller.
 *
 * param
 *     Reference to a OrderedPairs handle.
 */
void OrderedPairsDealloc(OrderedPairs_pa *OrderedPairs)
{
    REQUIRE(VALID_REF(OrderedPairs));
    REQUIRE(OrderedPairsIsValid(*OrderedPairs) || *OrderedPairs == NULL);

    if (*OrderedPairs != NULL)
    {
        /* release the ArrList's */
        if ((*OrderedPairs)->Item1 != NULL) ArrListDealloc(&((*OrderedPairs)->Item1));
        if ((*OrderedPairs)->Item2 != NULL) ArrListDealloc(&((*OrderedPairs)->Item2));

        /* release the list structure itself */
        FREE_ITEM(*OrderedPairs, "OrderedPairs structure");
        *OrderedPairs = NULL;
    }

    ENSURE(*OrderedPairs == NULL);
}





/**
 * Empties the OrderedPairs structure of all Bond Item Pairs.
 *
 * param OrderedPairs
 *     OrderedPairs to clear.
 */
void OrderedPairsClear(OrderedPairs_pa OrderedPairs)
{
    REQUIRE(OrderedPairsIsValid(OrderedPairs));

    ArrListClear(OrderedPairs->Item1);
    ArrListClear(OrderedPairs->Item2);

    ENSURE(OrderedPairsIsValid(OrderedPairs) && OrderedPairsGetCount(OrderedPairs) == 0);
}





/**
 * Removes a point from the OrderedPairs array. The members following the item
 * removed are shifted down accordingly to fill the vacated space.
 *
 *
 * param OrderedPairs
 *     OrderedPairs structure containing the bundle to remove.
 * param ItemOffset
 *     Offset to the point.
 *
 * return
 *     TRUE if successful, FALSE otherwise.
 */
Boolean_t OrderedPairsRemoveAtOffset(OrderedPairs_pa OrderedPairs,
                                     const LgIndex_t PointOffset)
{
    Boolean_t IsOk = TRUE;

    REQUIRE(OrderedPairsIsValid(OrderedPairs));
    REQUIRE(0 <= PointOffset && PointOffset <= OrderedPairsGetCount(OrderedPairs) - 1);

    /* Remove the items for the array lists */
    ArrListRemoveItem(OrderedPairs->Item1, PointOffset);
    ArrListRemoveItem(OrderedPairs->Item2, PointOffset);

    IsOk = OrderedPairsIsValid(OrderedPairs);

    ENSURE(OrderedPairsIsValid(OrderedPairs));
    ENSURE(VALID_BOOLEAN(IsOk));
    return IsOk;
}



/**
 * Allocates a OrderedPairs handle with a suitable default
 * capacity for the contained ArrList's.
 *
 * return
 *     OrderedPairs handle if sufficient memory was available,
 *     otherewise a handle to NULL.
 */
OrderedPairs_pa OrderedPairsAlloc()
{
    OrderedPairs_pa Result = NULL;

    Result = ALLOC_ITEM(OrderedPairs_s, "OrderedPairs structure");
    if (Result != NULL)
    {
        Result->Item1   = ArrListAlloc(60, ArrListType_Long);
        Result->Item2   = ArrListAlloc(60, ArrListType_Long);

        /* If it failed to allocate any of the array lists, clean-up and exit. */
        if (Result->Item1 == NULL || Result->Item2 == NULL)
        {
            if (Result->Item1 != NULL) ArrListDealloc(&(Result->Item1));
            if (Result->Item2 != NULL) ArrListDealloc(&(Result->Item1));
            FREE_ITEM(Result, "OrderedPairs structure");
            Result = NULL;
        }
    }

    ENSURE(OrderedPairsIsValid(Result) || Result == NULL);
    return Result;
}




/**
 * Gets the number of OrderedPairs currently in the OrderedPairs structure.
 *
 * param
 *     OrderedPairs structure in question.
 *
 * return
 *     Number of bond Item pairs in the OrderedPairs structure.
 */
LgIndex_t OrderedPairsGetCount(const OrderedPairs_pa OrderedPairs)
{
    LgIndex_t Result = 0;

    REQUIRE(OrderedPairsIsValid(OrderedPairs));

    Result = ArrListGetCount(OrderedPairs->Item1);

    ENSURE(Result >= 0);
    return Result;
}






/**
 * Places pair components at the specified offset. If the
 * offset is beyond the end of the array lists, they are sized
 * accordingly and the intervening array values between the last
 * node of the original state and the last node of the new state are
 * assigned 0.
 *
 * note
 *     If components already exists at the specified location they are
 *     replaced.
 *
 * param OrderedPairs
 *     OrderedPairs target in which to set the coordinates.
 * param PointOffset
 *     Offset of the point/node.
 * param Item1, Item2
 *     Components to set at the specified offset.
 *
 * return
 *     TRUE if sufficient memory permitted the operation, otherwise FALSE.
 */

Boolean_t OrderedPairsSetAtOffset(OrderedPairs_pa OrderedPairs,
                                  const LgIndex_t PointOffset,
                                  const LgIndex_t Item1,
                                  const LgIndex_t Item2)
{
    Boolean_t IsOk = TRUE;
    ArrListItem_u Item;


    REQUIRE(OrderedPairsIsValid(OrderedPairs));
    REQUIRE(PointOffset >= 0);

    Item.Long = Item1;
    IsOk = ArrListSetItem(OrderedPairs->Item1, PointOffset, Item);

    if (IsOk)
    {
        Item.Long = Item2;
        IsOk = ArrListSetItem(OrderedPairs->Item2, PointOffset, Item);
    }

    ENSURE(OrderedPairsIsValid(OrderedPairs));

    ENSURE(VALID_BOOLEAN(IsOk));
    return IsOk;
}





/**
 * Inserts pair components at the specified offset. The arrays
 * will be expanded to accomodate the additional value.
 *
 *
 * param OrderedPairs
 *     OrderedPairs target in which to set the coordinates.
 * param PointOffset
 *     Offset at which to insert the pair.
 * param Item1, Item2
 *     Components of the pair to inser.
 *
 * return
 *     TRUE if sufficient memory permitted the operation, otherwise FALSE.
 */

Boolean_t OrderedPairsInsertAtOffset(OrderedPairs_pa OrderedPairs,
                                     const LgIndex_t PointOffset,
                                     const LgIndex_t Item1,
                                     const LgIndex_t Item2)
{
    Boolean_t IsOk = TRUE;
    ArrListItem_u Item;


    REQUIRE(OrderedPairsIsValid(OrderedPairs));
    REQUIRE(0 <= PointOffset && PointOffset <= OrderedPairsGetCount(OrderedPairs));

    Item.Long = Item1;
    IsOk = ArrListInsertItem(OrderedPairs->Item1, PointOffset, Item);

    if (IsOk)
    {
        Item.Long = Item2;
        IsOk = ArrListInsertItem(OrderedPairs->Item2, PointOffset, Item);
    }

    ENSURE(OrderedPairsIsValid(OrderedPairs));

    ENSURE(VALID_BOOLEAN(IsOk));
    return IsOk;
}





/**
 * Appends the pair components to OrderedPairs structure. The array lists
 * will be expanded to accommodate the additional items.
 *
 * param OrderedPairs
 *     OrderedPairs target to which the pair components are to be appended.
 * param Item1, Item2
 *     Components of the pair to append to the OrderedPairs structure.
 *
 * return
 *     TRUE if sufficient memory permitted the operation, otherwise FALSE.
 */
Boolean_t OrderedPairsAppendAtEnd(OrderedPairs_pa  OrderedPairs,
                                  const LgIndex_t  Item1,
                                  const LgIndex_t  Item2)
{
    Boolean_t IsOk = TRUE;
    LgIndex_t Count;

    REQUIRE(OrderedPairsIsValid(OrderedPairs));

    Count = OrderedPairsGetCount(OrderedPairs);

    IsOk = OrderedPairsInsertAtOffset(OrderedPairs, Count, Item1, Item2);

    ENSURE(OrderedPairsIsValid(OrderedPairs));

    ENSURE(VALID_BOOLEAN(IsOk));
    return IsOk;
}






/**
 * Gets the components of the OrderedPairs structure for the
 * pair at the specified offset.
 *
 * param OrderedPairs
 *     OrderedPairs structure containing the desired item.
 * param PointOffset
 *     Offset to the components of the OrderedPairs.
 * param *Item1, *Item2
 *     Pointers to components of the pair
 *
 * return
 *     TRUE if it works, FALSE otherwise.
 */

Boolean_t OrderedPairsGetFromOffset(const OrderedPairs_pa  OrderedPairs,
                                    const LgIndex_t        PointOffset,
                                    LgIndex_t       *Item1,
                                    LgIndex_t       *Item2)
{
    Boolean_t IsOk = TRUE;
    ArrListItem_u Item;

    REQUIRE(OrderedPairsIsValid(OrderedPairs));
    REQUIRE(PointOffset >= 0 && PointOffset < OrderedPairsGetCount(OrderedPairs));
    REQUIRE(VALID_REF(Item1));
    REQUIRE(VALID_REF(Item2));

    Item = ArrListGetItem(OrderedPairs->Item1, PointOffset);
    *Item1 = Item.Long;

    Item = ArrListGetItem(OrderedPairs->Item2, PointOffset);
    *Item2 = Item.Long;

    ENSURE(VALID_REF(Item1));
    ENSURE(VALID_REF(Item2));

    ENSURE(VALID_BOOLEAN(IsOk));
    return IsOk;
}








/**
 * Return the offset to a specified pair.
 *
 * param OrderedPairs
 *     OrderedPairs structure containing the desired item.
 * param Item1, Item2
 *     Pair to compare with ordered-pair list
 *
 * return
 *     Offset if the pair exists, -1 if the pair is not in the list.
 */

LgIndex_t OrderedPairsGetOffsetFromPair(const OrderedPairs_pa OrderedPairs,
                                        const LgIndex_t       Item1,
                                        const LgIndex_t       Item2)
{
    LgIndex_t PairOffset = -1;
    LgIndex_t Count;
    LgIndex_t TestItem1, TestItem2;
    LgIndex_t ii;

    REQUIRE(OrderedPairsIsValid(OrderedPairs));

    Count = OrderedPairsGetCount(OrderedPairs);

    for (ii = 0; PairOffset < 0 && ii < Count; ii++)
    {
        Boolean_t IsOk = OrderedPairsGetFromOffset(OrderedPairs, ii, &TestItem1, &TestItem2);
        if (TestItem1 == Item1 && TestItem2 == Item2) PairOffset = ii;
    }

    ENSURE(PairOffset >= -1 && PairOffset < Count);
    return PairOffset;
}








/**
 * Check to see if a specified pair is unique (i.e. not already in the
 * ordered pair set).
 *
 * param OrderedPairs
 *     OrderedPairs structure containing the desired item.
 * param Item1, Item2
 *     Pair to compare with ordered-pair list
 *
 * return
 *     TRUE if the pair is unique, FALSE if the pair is in the list.
 */

Boolean_t OrderedPairsIsUniquePair(const OrderedPairs_pa  OrderedPairs,
                                   const LgIndex_t        Item1,
                                   const LgIndex_t        Item2)
{
    Boolean_t IsUnique = TRUE;
    LgIndex_t Count;
    LgIndex_t TestItem1, TestItem2;
    LgIndex_t ii;

    REQUIRE(OrderedPairsIsValid(OrderedPairs));

    Count = OrderedPairsGetCount(OrderedPairs);

    for (ii = 0; IsUnique && ii < Count; ii++)
    {
        Boolean_t IsOk = OrderedPairsGetFromOffset(OrderedPairs, ii, &TestItem1, &TestItem2);
        if (TestItem1 == Item1 && TestItem2 == Item2) IsUnique = FALSE;
    }

    ENSURE(VALID_BOOLEAN(IsUnique));
    return IsUnique;
}

