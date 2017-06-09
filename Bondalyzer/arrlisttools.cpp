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
#include "ALLOC.h"
#include "ARRLIST.h"
#include "ARRLISTTOOLS.h"


/**
 * Check to see if an item is already in a Long ArrList.
 *
 * param ArrList
 *     Array list to which the item is to be compared.
 * param LongArrListItem
 *     Item to be compared.
 *
 * return
 *     TRUE if item is unique in ArrList, FALSE if it is already in the list.
 */
Boolean_t ArrListIsUniqueLongItem(const ArrList_pa    ArrList,
                                  const LgIndex_t     LongArrListItem)
{
    Boolean_t IsUnique = TRUE;
    ArrListItem_u Item;
    int ii;
    LgIndex_t NumItems = ArrListGetCount(ArrList);

    REQUIRE(ArrListIsValid(ArrList));
    REQUIRE(ArrListGetType(ArrList) == ArrListType_Long);  // Limitation for now

    // Check to see if it is unique (can't use binary search - not sorted)
    for (ii = 0; IsUnique && ii < NumItems; ii++)
    {
        Item = ArrListGetItem(ArrList, ii);
        if (LongArrListItem == Item.Long)
            IsUnique = FALSE;
    }

    ENSURE(VALID_BOOLEAN(IsUnique));
    return IsUnique;
}




/**
 * Append an Item to the node list only if the item doesn't already
 * exist in the array list.
 *
 * param ArrList
 *     Array list to which the item is to be appended.
 * param LongArrListItem
 *     Item to be appended (if unique).
 *
 * return
 *     Offet of item in list if it was unique or appended, -1 if failed.
 */
LgIndex_t ArrListAppendUniqueLongItem(ArrList_pa    ArrList,
                                      LgIndex_t     LongArrListItem)
{
    LgIndex_t Offset = -1;
    Boolean_t IsUnique = TRUE;
    Boolean_t IsOk = TRUE;
    ArrListItem_u Item;
    ArrListItem_u ItemOld;
    int ii;
    LgIndex_t NumItems = ArrListGetCount(ArrList);

    REQUIRE(ArrListIsValid(ArrList));
    REQUIRE(ArrListGetType(ArrList) == ArrListType_Long);  // Limitation for now

    // Check to see if it is unique (can't use binary search - not sorted)
    for (ii = 0; IsUnique && ii < NumItems; ii++)
    {
        ItemOld = ArrListGetItem(ArrList, ii);
        if (LongArrListItem == ItemOld.Long)
        {
            IsUnique = FALSE;
            Offset = ii;
        }
    }

    if (IsUnique)
    {
        Item.Long = LongArrListItem;
        IsOk = ArrListAppendItem(ArrList, Item);
        Offset = ArrListGetCount(ArrList) - 1;
    }

    // Unique but append failed (unlikely)
    CHECK(IsOk);
    if (!IsOk) Offset = -1;

    ENSURE(ArrListIsValid(ArrList));
    ENSURE(Offset >= -1 || Offset < ArrListGetCount(ArrList));
    return Offset;
}


