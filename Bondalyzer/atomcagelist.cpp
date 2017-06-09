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
#include "ATOMCAGELIST.h"
#include <string.h>


#if 0 // not currently used

/**
 * Determine if the AtomCageList handle is sane.
 *
 * param AtomCageList
 *     AtomCageList structure in question.
 *
 * return
 *     TRUE if the AtomCageList structure is valid, otherwise FALSE.
 */
Boolean_t AtomCageListIsValid(AtomCageList_pa AtomCageList)
{
    Boolean_t IsValid = FALSE;

    IsValid = (VALID_REF(AtomCageList) &&
               VALID_REF(AtomCageList->Cages) && ArrListIsValid(AtomCageList->Cages));

    ENSURE(VALID_BOOLEAN(IsValid));
    return IsValid;
}




/**
 * Deallocates the AtomCageList handle and set the handle to NULL.
 *
 * note
 *     Item destruction is the responsibility of the caller.
 *
 * param
 *     Reference to a AtomCageList handle.
 */
void AtomCageListDealloc(AtomCageList_pa *AtomCageList)
{
    REQUIRE(VALID_REF(AtomCageList));
    REQUIRE(AtomCageListIsValid(*AtomCageList) || *AtomCageList == NULL);

    if (*AtomCageList != NULL)
    {
        /* release the ArrList's */
        if ((*AtomCageList)->Cages != NULL) ArrListDealloc(&((*AtomCageList)->Cages));

        /* release the list structure itself */
        FREE_ITEM(*AtomCageList, "AtomCageList structure");
        *AtomCageList = NULL;
    }

    ENSURE(*AtomCageList == NULL);
}





/**
 * Empties the AtomCageList structure of all Bond Atom Pairs.
 *
 * param AtomCageList
 *     AtomCageList to clear.
 */
void AtomCageListClear(AtomCageList_pa AtomCageList)
{
    REQUIRE(AtomCageListIsValid(AtomCageList));

    ArrListClear(AtomCageList->Cages);

    ENSURE(AtomCageListIsValid(AtomCageList) && AtomCageListGetCount(AtomCageList) == 0);
}





/**
 * Removes a point from the AtomCageList array. The members following the item
 * removed are shifted down accordingly to fill the vacated space.
 *
 *
 * param AtomCageList
 *     AtomCageList structure containing the bundle to remove.
 * param ItemOffset
 *     Offset to the point.
 *
 * return
 *     TRUE if successful, FALSE otherwise.
 */
Boolean_t AtomCageListRemoveAtOffset(AtomCageList_pa AtomCageList,
                                     LgIndex_t      PointOffset)
{
    Boolean_t IsOk = TRUE;

    REQUIRE(AtomCageListIsValid(AtomCageList));
    REQUIRE(0 <= PointOffset && PointOffset <= AtomCageListGetCount(AtomCageList) - 1);

    /* Remove the items for the array lists */
    ArrListRemoveItem(AtomCageList->Cages, PointOffset);

    IsOk = AtomCageListIsValid(AtomCageList);

    ENSURE(AtomCageListIsValid(AtomCageList));
    ENSURE(VALID_BOOLEAN(IsOk));
    return IsOk;
}



/**
 * Allocates a AtomCageList handle with a suitable default
 * capacity for the contained ArrList's.
 *
 * return
 *     AtomCageList handle if sufficient memory was available,
 *     otherewise a handle to NULL.
 */
AtomCageList_pa AtomCageListAlloc()
{
    AtomCageList_pa Result = NULL;

    Result = ALLOC_ITEM(AtomCageList_s, "AtomCageList structure");
    if (Result != NULL)
    {
        Result->MemUsageReported = 0;
        Result->Cages   = ArrListAlloc(60, ArrListType_Long);

        /* If it failed to allocate any of the array lists, clean-up and exit. */
        if (Result->Cages == NULL)
        {
            if (Result->Cages != NULL) ArrListDealloc(&(Result->Cages));
            FREE_ITEM(Result, "AtomCageList structure");
            Result = NULL;
        }
    }

    ENSURE(AtomCageListIsValid(Result) || Result == NULL);
    return Result;
}




/**
 * Gets the number of AtomCageList currently in the AtomCageList structure.
 *
 * param
 *     AtomCageList structure in question.
 *
 * return
 *     Number of bond atom pairs in the AtomCageList structure.
 */
LgIndex_t AtomCageListGetCount(AtomCageList_pa AtomCageList)
{
    LgIndex_t Result = 0;

    REQUIRE(AtomCageListIsValid(AtomCageList));

    Result = ArrListGetCount(AtomCageList->Cages);

    ENSURE(Result >= 0);
    return Result;
}






/**
 * Places Cage numbers at the specified offset. If the
 * offset is beyond the end of the array lists, they are sized
 * accordingly and the intervening array values between the last
 * node of the original state and the last node of the new state are
 * assigned 0.
 *
 * note
 *     If components already exists at the specified location they are
 *     replaced.
 *
 * param AtomCageList
 *     AtomCageList target in which to set the coordinates.
 * param PointOffset
 *     Offset of the point/node.
 * param Cage
 *     Cage number to set at the specified offset.
 *
 * return
 *     TRUE if sufficient memory permitted the operation, otherwise FALSE.
 */

Boolean_t AtomCageListSetAtOffset(AtomCageList_pa AtomCageList,
                                  LgIndex_t  PointOffset,
                                  LgIndex_t  Cage)
{
    Boolean_t IsOk = TRUE;
    ArrListItem_u Item;


    REQUIRE(AtomCageListIsValid(AtomCageList));
    REQUIRE(PointOffset >= 0);

    Item.Long = Cage;
    IsOk = ArrListSetItem(AtomCageList->Cages, PointOffset, Item);

    ENSURE(AtomCageListIsValid(AtomCageList));

    ENSURE(VALID_BOOLEAN(IsOk));
    return IsOk;
}





/**
 * Inserts critical point components at the specified offset. The arrays
 * will be expanded to accomodate the additional value.
 *
 *
 * param AtomCageList
 *     AtomCageList target in which to set the coordinates.
 * param PointOffset
 *     Offset at which to insert the point/node.
 * param Cage
 *     Cage number to inser at PointOffset.
 *
 * return
 *     TRUE if sufficient memory permitted the operation, otherwise FALSE.
 */

Boolean_t AtomCageListInsertAtOffset(AtomCageList_pa AtomCageList,
                                     LgIndex_t   PointOffset,
                                     LgIndex_t   Cage)
{
    Boolean_t IsOk = TRUE;
    ArrListItem_u Item;


    REQUIRE(AtomCageListIsValid(AtomCageList));
    REQUIRE(0 <= PointOffset && PointOffset <= AtomCageListGetCount(AtomCageList));

    Item.Long = Cage;
    IsOk = ArrListInsertItem(AtomCageList->Cages, PointOffset, Item);

    ENSURE(AtomCageListIsValid(AtomCageList));

    ENSURE(VALID_BOOLEAN(IsOk));
    return IsOk;
}





/**
 * Appends the bundle components to AtomCageList structure. The array list
 * will be expanded to accommodate the additional items.
 *
 * param AtomCageList
 *     AtomCageList target to which the bundle components are to be appended.
 * param Cage
 *     Components of the bundle to append to the AtomCageList structure.
 *
 * return
 *     TRUE if sufficient memory permitted the operation, otherwise FALSE.
 */
Boolean_t AtomCageListAppendAtEnd(AtomCageList_pa  AtomCageList,
                                  LgIndex_t        Cage)
{
    Boolean_t IsOk = TRUE;
    LgIndex_t Count;

    REQUIRE(AtomCageListIsValid(AtomCageList));

    Count = AtomCageListGetCount(AtomCageList);

    IsOk = AtomCageListInsertAtOffset(AtomCageList, Count, Cage);

    ENSURE(AtomCageListIsValid(AtomCageList));

    ENSURE(VALID_BOOLEAN(IsOk));
    return IsOk;
}






/**
 * Gets the components of the AtomCageList structure for the
 * bundle at the specified offset.
 *
 * param AtomCageList
 *     AtomCageList structure containing the desired item.
 * param PointOffset
 *     Offset to the components of the AtomCageList.
 *
 * return
 *     Cage number if it works, -1 otherwise.
 */

LgIndex_t AtomCageListGetFromOffset(AtomCageList_pa  AtomCageList,
                                    LgIndex_t   PointOffset)
{
    LgIndex_t Cage = -1;
    Boolean_t IsOk = TRUE;
    ArrListItem_u Item;

    REQUIRE(AtomCageListIsValid(AtomCageList));
    REQUIRE(PointOffset >= 0 && PointOffset < AtomCageListGetCount(AtomCageList));

    Item = ArrListGetItem(AtomCageList->Cages, PointOffset);
    Cage = Item.Long;

    ENSURE(Cage >= -1);
    return Cage;
}





#endif
