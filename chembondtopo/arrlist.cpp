/*
*****************************************************************
*****************************************************************
*******                                                  ********
******      (C) Copyright 1989-2004  by TECPLOT INC.     ********
*******               All Rights Reserved.               ********
*******                                                  ********
*****************************************************************
*****************************************************************
*/
#ifdef IRISX
#include <sys/types.h>
extern pid_t __vfork(void);
#endif

#include "MASTER.h"


#define ARRLISTMODULE
#include "GLOBAL.h"
#include "TASSERT.h"
#include "ALLOC.h"
#include "ARRLIST.h"

/* default capacity if none is given and also used for block sizing */
#define DEFAULT_CAPACITY 32


/*
 * ABSTRACT:
 *
 * This general purpose list uses an array implementation. The array item
 * is guarenteed to be large enough to hold the intrinsic 'C' types such
 * as char, int, long, double, and any pointer. This list has no notion
 * of item allocation/deallocation so it is the client's responsibility to
 * release any referenced memory that is no longer needed.
 */


/* private ArrList structure */
typedef struct _ArrList_s
{
    ArrListType_e   Type;     /* type of array items */
    EntIndex_t      ItemSize; /* byte size of an individual item */
    char            *Array;   /* byte array for holding the items */
    LgIndex_t       Count;    /* number of items in the array */
    LgIndex_t       Capacity; /* maximum holding capacity of the array */
} ArrList_s;


/**
 * Copies the private array items from the specified source to the target
 * location. The buffers may overlap.
 *
 * note
 *     Originally this function was a macro that called memmove
 *     directly:
 *
 *         #define CopyArrayItems(TargetArray, TargetOffset, \
 *                                SourceArray, SourceOffset, \
 *                                Count, ItemSize) \
 *                     (memmove(&((TargetArray)[(TargetOffset)*ItemSize]), \
 *                              &((SourceArray)[(SourceOffset)*ItemSize]), \
 *                              Count*ItemSize))
 *
 * This however proved troublesome as some machines replaced the memmove
 * with a call to memcpy in the linker. The memcpy function does not support
 * overlapping moves so I could not use it. This function should be just
 * about as fast however so it is no big deal.
 *
 * param TargetArray
 *     Base address of the target array to receive the items.
 * param TargetOffset
 *     Target offset of the first item.
 * param SourceArray
 *     Base address of the source array supplying the items.
 * param SourceOffset
 *     Source offset of the first item.
 * param Count
 *     Number of items to copy.
 * param ItemSize
 *     Item size in bytes.
 */
static void CopyArrayItems(char       *TargetArray,
                           LgIndex_t  TargetOffset,
                           char       *SourceArray,
                           LgIndex_t  SourceOffset,
                           LgIndex_t  Count,
                           EntIndex_t ItemSize)
{
    LgIndex_t Index = 0;
    LgIndex_t CharCount = 0;
    char      *TargetPtr = NULL;
    char      *SourcePtr = NULL;

    REQUIRE(VALID_REF(TargetArray));
    REQUIRE(TargetOffset >= 0);
    REQUIRE(VALID_REF(SourceArray));
    REQUIRE(SourceOffset >= 0);
    REQUIRE(&TargetArray[TargetOffset] != &SourceArray[SourceOffset]);
    REQUIRE(Count >= 1);
    REQUIRE(1 <= ItemSize && ItemSize <= sizeof(ArrListItem_u));

    CharCount = Count * ItemSize;
    if (&TargetArray[TargetOffset] < &SourceArray[SourceOffset])
    {
        TargetPtr = &TargetArray[TargetOffset * ItemSize];
        SourcePtr = &SourceArray[SourceOffset * ItemSize];
        for (Index = 0; Index < CharCount; Index++, TargetPtr++, SourcePtr++)
            *TargetPtr = *SourcePtr;
    }
    else
    {
        TargetPtr = &TargetArray[TargetOffset * ItemSize + CharCount - 1];
        SourcePtr = &SourceArray[SourceOffset * ItemSize + CharCount - 1];
        for (Index = 0; Index < CharCount; Index++, TargetPtr--, SourcePtr--)
            *TargetPtr = *SourcePtr;
    }
}


/**
 * Adjusts the capacity request as necessary to minimize memory reallocations
 * for large lists. The adjusted capacity will be at least as big as requested
 * however it may be larger if it is determined that the space requirement is
 * growing faster.
 *
 * param ArrList
 *     Current capacity used as a helpful hint for the adjustment algorythm.
 * param RequestedCapacity
 *     Capacity request.
 *
 * return
 *     Adjusted capacity that is at least as large as the request.
 */
static LgIndex_t AdjustCapacityRequest(ArrList_pa ArrList,
                                       LgIndex_t  RequestedCapacity)
{
    LgIndex_t BLOCK_SIZE = DEFAULT_CAPACITY;
    LgIndex_t Result = 0;

    REQUIRE(ArrListIsValid(ArrList));
    REQUIRE(RequestedCapacity > ArrList->Capacity);

    /* compute a reasonable block size adjust the request to fit within */
    BLOCK_SIZE = MAX(BLOCK_SIZE, ArrList->Capacity / 2);
    Result = ((RequestedCapacity - 1) / BLOCK_SIZE + 1) * BLOCK_SIZE;

    ENSURE(Result >= RequestedCapacity);
    return Result;
}


/**
 * Enlarge the list capacity to accommodate, at a minimum, the requested
 * capacity.
 *
 * param ArrList
 *     Current capacity used as a helpful hint for the adjustment algorythm.
 * param RequestedCapacity
 *     Capacity request.
 *
 * return
 *     TRUE if the list could be enlarged, otherwise FALSE.
 */
static Boolean_t EnlargeListCapacity(ArrList_pa ArrList,
                                     LgIndex_t  RequestedCapacity)
{
    Boolean_t IsOk = FALSE;
    LgIndex_t AdjustedCapacity = 0;
    char      *EnlargedArray = NULL;

    REQUIRE(ArrListIsValid(ArrList));
    REQUIRE(RequestedCapacity > ArrList->Capacity);

    AdjustedCapacity = AdjustCapacityRequest(ArrList, RequestedCapacity);
    EnlargedArray = ALLOC_ARRAY(AdjustedCapacity * ArrList->ItemSize,
                                char, "array list");
    if (EnlargedArray == NULL && RequestedCapacity < AdjustedCapacity)
    {
        /* try again with minimum capacity request */
        AdjustedCapacity = RequestedCapacity;
        EnlargedArray = ALLOC_ARRAY(AdjustedCapacity * ArrList->ItemSize,
                                    char, "array list");
    }
    IsOk = (EnlargedArray != NULL);
    if (IsOk)
    {
        /* copy the items to the new list, release the old list, */
        /* and record the new list in the list structure         */
        if (ArrList->Array != NULL)
        {
            if (ArrList->Count != 0)
                CopyArrayItems(EnlargedArray, 0,
                               ArrList->Array, 0,
                               ArrList->Count,
                               ArrList->ItemSize);
            FREE_ARRAY(ArrList->Array, "array list");
        }
        ArrList->Array = EnlargedArray;
        ArrList->Capacity = AdjustedCapacity;
    }

    ENSURE(ArrListIsValid(ArrList));
    ENSURE(VALID_BOOLEAN(IsOk));
    return IsOk;
}


/**
 * Gets the size of an individual element.
 *
 * param Type
 *     Array list element type.
 *
 * return
 *     Element size corresponding to the type.
 */
static EntIndex_t GetElementSize(ArrListType_e Type)
{
    EntIndex_t Result = 0;

    REQUIRE(VALID_ENUM(Type, ArrListType_e));

    switch (Type)
    {
        case ArrListType_UnsignedChar:
            Result = sizeof(unsigned char);
            break;
        case ArrListType_UnsignedShort:
            Result = sizeof(unsigned short);
            break;
        case ArrListType_UnsignedInt:
            Result = sizeof(unsigned int);
            break;
        case ArrListType_UnsignedLong:
            Result = sizeof(unsigned long);
            break;
        case ArrListType_Char:
            Result = sizeof(char);
            break;
        case ArrListType_Short:
            Result = sizeof(short);
            break;
        case ArrListType_Int:
            Result = sizeof(int);
            break;
        case ArrListType_Long:
            Result = sizeof(long);
            break;
        case ArrListType_Double:
            Result = sizeof(double);
            break;
        case ArrListType_UnsignedCharPtr:
            Result = sizeof(unsigned char *);
            break;
        case ArrListType_UnsignedShortPtr:
            Result = sizeof(unsigned short *);
            break;
        case ArrListType_UnsignedIntPtr:
            Result = sizeof(unsigned int *);
            break;
        case ArrListType_UnsignedLongPtr:
            Result = sizeof(unsigned long *);
            break;
        case ArrListType_CharPtr:
            Result = sizeof(char *);
            break;
        case ArrListType_ShortPtr:
            Result = sizeof(short *);
            break;
        case ArrListType_IntPtr:
            Result = sizeof(int *);
            break;
        case ArrListType_LongPtr:
            Result = sizeof(long *);
            break;
        case ArrListType_DoublePtr:
            Result = sizeof(double *);
            break;
        case ArrListType_VoidPtr:
            Result = sizeof(void *);
            break;
        case ArrListType_FunctionPtr:
            Result = sizeof(void (*)());
            break;
        case ArrListType_Any: /* allows a mixed bag of items */
            Result = sizeof(ArrListItem_u);
            break;

        default:
            CHECK(FALSE);
            break;
    }

    ENSURE(1 <= Result && Result <= sizeof(ArrListItem_u));
    return Result;
}


/**
 * Determine if the list handle is sane.
 *
 * param ArrList
 *     Array list in question.
 *
 * return
 *     TRUE if the array list is valid, otherwise FALSE.
 */
Boolean_t ArrListIsValid(ArrList_pa ArrList)
{
    Boolean_t IsValid = FALSE;

    IsValid = (VALID_REF(ArrList) &&
               VALID_ENUM(ArrList->Type, ArrListType_e) &&
               (1 <= ArrList->ItemSize &&
                ArrList->ItemSize <= sizeof(ArrListItem_u)) &&
               (0 <= ArrList->Count &&
                ArrList->Count <= ArrList->Capacity));

    ENSURE(VALID_BOOLEAN(IsValid));
    return IsValid;
}


/**
 * Gets the specified array list's type.
 *
 * param ArrList
 *     Array list of which the type is desired.
 *
 * return
 *     Array list type.
 */
ArrListType_e ArrListGetType(ArrList_pa ArrList)
{
    ArrListType_e Result = ArrListType_Invalid;

    REQUIRE(ArrListIsValid(ArrList));

    Result = ArrList->Type;

    ENSURE(VALID_ENUM(Result, ArrListType_e));
    return Result;
}


/**
 * Allocates an array list handle with the estimated capacity
 * or a suitable default if an estimate is unavailable.
 *
 * param EstimatedCapacity
 *     Clients best guess at the estimated capacity need. If
 *     an estimate is not available zero the zero should be
 *     used to get the default capacity.
 * param Type
 *     Type of array list being allocated.
 *
 * return
 *     Array list handle if sufficient memory was available,
 *     otherewise a handle to NULL.
 */
ArrList_pa ArrListAlloc(LgIndex_t     EstimatedCapacity,
                        ArrListType_e Type)
{
    ArrList_pa Result = NULL;

    REQUIRE(EstimatedCapacity >= 0);
    REQUIRE(VALID_ENUM(Type, ArrListType_e));

    Result = ALLOC_ITEM(ArrList_s, "ArrList structure");
    if (Result != NULL)
    {
        if (EstimatedCapacity == 0)
            EstimatedCapacity = DEFAULT_CAPACITY;

        Result->Type = Type;
        Result->ItemSize = GetElementSize(Type);
        Result->Array = ALLOC_ARRAY(EstimatedCapacity * Result->ItemSize,
                                    char, "array list");
        Result->Count = 0;
        if (Result->Array != NULL)
        {
            Result->Capacity = EstimatedCapacity;
        }
        else
        {
            /* estimated capacity allocation failed so clean up */
            Result->Capacity = 0;
            ArrListDealloc(&Result);
        }
    }

    ENSURE(ArrListIsValid(Result) || Result == NULL);
    ENSURE(IMPLICATION(Result != NULL, Result->Capacity >= EstimatedCapacity));
    return Result;
}


/**
 * Deallocates the list handle and set the handle to NULL.
 *
 * note
 *     Item destruction is the responsibility of the caller.
 *
 * param
 *     Reference to an array list handle.
 */
void ArrListDealloc(ArrList_pa *ArrList)
{
    REQUIRE(VALID_REF(ArrList));
    REQUIRE(ArrListIsValid(*ArrList) || *ArrList == NULL);

    if (*ArrList != NULL)
    {
        /* release the list */
        if ((*ArrList)->Capacity != 0)
            FREE_ARRAY((*ArrList)->Array, "array list");

        /* release the list structure itself */
        FREE_ITEM(*ArrList, "ArrList structure");
        *ArrList = NULL;
    }

    ENSURE(*ArrList == NULL);
}


/**
 * Gets the number of items currently maintained by the list.
 *
 * param
 *     Array list in question.
 *
 * return
 *     Number of items maintained by the list.
 */
LgIndex_t ArrListGetCount(ArrList_pa ArrList)
{
    LgIndex_t Result = 0;

    REQUIRE(ArrListIsValid(ArrList));

    Result = ArrList->Count;

    ENSURE(Result >= 0);
    return Result;
}


/**
 * Empties the array list of all items.
 *
 * note
 *     It is the clients responsibility to deallocate items prior to this
 *     call if necessary.
 *
 * param ArrList
 *     Array list to clear.
 */
void ArrListClear(ArrList_pa ArrList)
{
    REQUIRE(ArrListIsValid(ArrList));

    ArrList->Count = 0;

    ENSURE(ArrListIsValid(ArrList) && ArrList->Count == 0);
}


/**
 * Clears 'Count' items from the array list. The members following the
 * items cleared are shifted down accordingly to fill the vacated space.
 *
 * param ArrList
 *     Array list containing the items to clear.
 * param ItemOffset
 *     Offset to the first item to clear in the list.
 * param Count
 *     Number of items to clear.
 */
void ArrListClearItems(ArrList_pa ArrList,
                       LgIndex_t  ItemOffset,
                       LgIndex_t  Count)
{
    REQUIRE(ArrListIsValid(ArrList));
    REQUIRE(0 <= ItemOffset && ItemOffset <= ArrList->Count - 1);
    REQUIRE(1 <= Count && ItemOffset + Count <= ArrList->Count);

    /* if we cleared the items from the middle of the array then     */
    /* shift the end items down by 'Count' to fill the vacated space */
    if (ItemOffset + Count < ArrList->Count)
        CopyArrayItems(ArrList->Array, ItemOffset,
                       ArrList->Array, ItemOffset + Count,
                       ArrList->Count - (ItemOffset + Count),
                       ArrList->ItemSize);

    /* update the count but leave the capacity alone */
    ArrList->Count -= Count;

    ENSURE(ArrListIsValid(ArrList));
}


/**
 * Clears an item from the array list. The members following the item
 * cleared are shifted down accordingly to fill the vacated space.
 *
 * param ArrList
 *     Array list containing the item to clear.
 * param ItemOffset
 *     Offset to the item in the list.
 */
void ArrListClearItem(ArrList_pa ArrList,
                      LgIndex_t  ItemOffset)
{
    REQUIRE(ArrListIsValid(ArrList));
    REQUIRE(0 <= ItemOffset && ItemOffset <= ArrList->Count - 1);

    ArrListClearItems(ArrList, ItemOffset, 1);

    ENSURE(ArrListIsValid(ArrList));
}


/**
 * Removes 'Count' items from the array list beginning at the specified
 * item offset. The members following the items removed are shifted down
 * accordingly to fill the vacated space.
 *
 * param ArrList
 *     Array list containing the items to remove.
 * param ItemOffset
 *     Offset to the first item to remove in the list.
 * param Count
 *     Number of items to remove.
 *
 * return
 *     Array list handle referring to the removed items if sufficient
 *     memory was available, otherewise a handle to NULL.
 */
ArrList_pa ArrListRemoveItems(ArrList_pa ArrList,
                              LgIndex_t  ItemOffset,
                              LgIndex_t  Count)
{
    ArrList_pa Result = NULL;

    REQUIRE(ArrListIsValid(ArrList));
    REQUIRE(0 <= ItemOffset && ItemOffset <= ArrList->Count - 1);
    REQUIRE(1 <= Count && ItemOffset + Count <= ArrList->Count);

    /* get a copy of the items and clear them from the source */
    Result = ArrListGetItems(ArrList, ItemOffset, Count);
    if (Result != NULL)
        ArrListClearItems(ArrList, ItemOffset, Count);

    ENSURE(ArrListIsValid(ArrList));
    ENSURE(ArrListIsValid(Result) || Result == NULL);
    return Result;
}


/**
 * Removes an item from the array list. The members following the item
 * removed are shifted down accordingly to fill the vacated space.
 *
 * note
 *     It is the clients responsibility to deallocate the removed item
 *     if necessary.
 *
 * param ArrList
 *     Array list containing the item to remove.
 * param ItemOffset
 *     Offset to the item in the list.
 *
 * return
 *     Item removed from the array list.
 */
ArrListItem_u ArrListRemoveItem(ArrList_pa ArrList,
                                LgIndex_t  ItemOffset)
{
    ArrListItem_u Result;

    REQUIRE(ArrListIsValid(ArrList));
    REQUIRE(0 <= ItemOffset && ItemOffset <= ArrList->Count - 1);

    /* record the orginal item */
    CopyArrayItems((char *)&Result, 0,
                   ArrList->Array, ItemOffset,
                   1, ArrList->ItemSize);

    /* clear the item from the array */
    ArrListClearItems(ArrList, ItemOffset, 1);

    ENSURE(ArrListIsValid(ArrList));
    return Result;
}


/**
 * Inserts copies of the items from the source list to the target list at
 * the specified offset. The target list will expand to accommodate the
 * additional items. The source list remains unchanged.
 *
 * param Target
 *     Array list receiving the source items.
 * param ItemOffset
 *     Offset at which to insert the source list items.
 * param Source
 *     Array list supplying the source items.
 *
 * return
 *     TRUE if sufficient memory permitted the operation, otherwise FALSE.
 */
Boolean_t ArrListInsert(ArrList_pa Target,
                        LgIndex_t  ItemOffset,
                        ArrList_pa Source)
{
    Boolean_t IsOk = TRUE;

    REQUIRE(ArrListIsValid(Target));
    REQUIRE(0 <= ItemOffset && ItemOffset <= Target->Count);
    REQUIRE(ArrListIsValid(Source));
    REQUIRE(Target != Source);
    REQUIRE(Target->Type == Source->Type);

    if (Source->Count != 0)
    {
        /* if necessary enlarge the target list to accommodate the request */
        if (Target->Count + Source->Count > Target->Capacity)
            IsOk = EnlargeListCapacity(Target, Target->Count + Source->Count);

        if (IsOk)
        {
            /* shift all items in the target list ahead of the  */
            /* insert position up by the number of items in the */
            /* source list to make room for the new items       */
            if (ItemOffset < Target->Count)
                CopyArrayItems(Target->Array, ItemOffset + Source->Count,
                               Target->Array, ItemOffset,
                               Target->Count - ItemOffset,
                               Target->ItemSize);

            /* insert the items and update the count */
            CopyArrayItems(Target->Array, ItemOffset,
                           Source->Array, 0,
                           Source->Count, Source->ItemSize);
            Target->Count += Source->Count;
        }
    }

    ENSURE(ArrListIsValid(Target));
    ENSURE(VALID_BOOLEAN(IsOk));
    return IsOk;
}


/**
 * Inserts the item into the array list at the specified offset.
 * The list will be expanded to accommodate the additonal item.
 *
 * param ArrList
 *     Array list target in which to insert the item.
 * param ItemOffset
 *     Offset at which to insert the item.
 * param Item
 *     Item to insert at the specified offset.
 *
 * return
 *     TRUE if sufficient memory permitted the operation, otherwise FALSE.
 */
Boolean_t ArrListInsertItem(ArrList_pa    ArrList,
                            LgIndex_t     ItemOffset,
                            ArrListItem_u Item)
{
    Boolean_t IsOk = TRUE;

    REQUIRE(ArrListIsValid(ArrList));
    REQUIRE(0 <= ItemOffset && ItemOffset <= ArrList->Count);

    /* if necessary enlarge the list to accommodate the request */
    if (ArrList->Count + 1 > ArrList->Capacity)
        IsOk = EnlargeListCapacity(ArrList, ArrList->Count + 1);

    if (IsOk)
    {
        /* shift all items in the target list ahead of the insert */
        /* position up by one to make room for the new item       */
        if (ItemOffset < ArrList->Count)
            CopyArrayItems(ArrList->Array, ItemOffset + 1,
                           ArrList->Array, ItemOffset,
                           ArrList->Count - ItemOffset,
                           ArrList->ItemSize);

        /* insert the item and update the count */
        CopyArrayItems(ArrList->Array, ItemOffset,
                       (char *)&Item, 0,
                       1, ArrList->ItemSize);
        ArrList->Count++;
    }

    ENSURE(ArrListIsValid(ArrList));
    ENSURE(VALID_BOOLEAN(IsOk));
    return IsOk;
}


/**
 * Gets copies of 'Count' items from the array list beginning at the
 * specified item offset.
 *
 * param ArrList
 *     Array list containing the items to copy.
 * param ItemOffset
 *     Offset to the first item to copy from the list.
 * param Count
 *     Number of items to copy.
 *
 * return
 *     Array list handle referring to the copied items if sufficient
 *     memory was available, otherewise a handle to NULL.
 */
ArrList_pa ArrListGetItems(ArrList_pa ArrList,
                           LgIndex_t  ItemOffset,
                           LgIndex_t  Count)
{
    ArrList_pa Result = NULL;

    REQUIRE(ArrListIsValid(ArrList));
    REQUIRE(0 <= ItemOffset && ItemOffset <= ArrList->Count - 1);
    REQUIRE(1 <= Count && ItemOffset + Count <= ArrList->Count);

    Result = ArrListAlloc(Count, ArrList->Type);
    if (Result != NULL)
    {
        /* copy the original items into the result */
        CopyArrayItems(Result->Array, 0,
                       ArrList->Array, ItemOffset,
                       Count, ArrList->ItemSize);
        Result->Count = Count;
    }

    ENSURE(ArrListIsValid(ArrList));
    ENSURE(ArrListIsValid(Result) || Result == NULL);
    return Result;
}


/**
 * Gets the item at the specified offset in the list.
 *
 * param ArrList
 *     Array list containing the desired item.
 * param ItemOffset
 *     Offset to the item in the list.
 *
 * return
 *     The requested item.
 */
ArrListItem_u ArrListGetItem(ArrList_pa ArrList,
                             LgIndex_t  ItemOffset)
{
    ArrListItem_u Result;

    REQUIRE(ArrListIsValid(ArrList));
    REQUIRE(0 <= ItemOffset && ItemOffset <= ArrList->Count - 1);

    CopyArrayItems((char *)&Result, 0,
                   ArrList->Array, ItemOffset,
                   1, ArrList->ItemSize);

    return Result;
}


/**
 * Places the item at the specified offset. If the offset is beyond the
 * end of the list it is sized accordingly and the intervening items
 * between the last item of the original state and the last item of the
 * new state are assigned 0.
 *
 * note
 *     If an item already exists at the specified location it is replaced
 *     therefore item destruction is the responsibility of the caller.
 *
 * param ArrList
 *     Array list target in which to set the item.
 * param ItemOffset
 *     Offset of the item.
 * param Item
 *     Item to set at the specified offset.
 *
 * return
 *     TRUE if sufficient memory permitted the operation, otherwise FALSE.
 */
Boolean_t ArrListSetItem(ArrList_pa    ArrList,
                         LgIndex_t     ItemOffset,
                         ArrListItem_u Item)
{
    Boolean_t IsOk = TRUE;

    REQUIRE(ArrListIsValid(ArrList));
    REQUIRE(ItemOffset >= 0);

    /* if necessary enlarge the list to accommodate the request */
    if (ItemOffset + 1 > ArrList->Capacity)
        IsOk = EnlargeListCapacity(ArrList, ItemOffset + 1);

    if (IsOk)
    {
        if (ItemOffset + 1 > ArrList->Count)
        {
            /* fill intervening items between the original last item */
            /* and the new last item with zeros; update the count    */
            if (ItemOffset > ArrList->Count)
                memset(&ArrList->Array[ArrList->Count*ArrList->ItemSize],
                       0, (ItemOffset - ArrList->Count)*ArrList->ItemSize);
            ArrList->Count = ItemOffset + 1;
        }
        CopyArrayItems(ArrList->Array, ItemOffset,
                       (char *)&Item, 0,
                       1, ArrList->ItemSize);
    }

    ENSURE(ArrListIsValid(ArrList));
    ENSURE(VALID_BOOLEAN(IsOk));
    return IsOk;
}


/**
 * Appends the item to the list. The list will be expanded
 * to accommodate the additional item.
 *
 * param ArrList
 *     Array list target to which the item is to be appended.
 * param Item
 *     Item to append to the array list.
 *
 * return
 *     TRUE if sufficient memory permitted the operation, otherwise FALSE.
 */
Boolean_t ArrListAppendItem(ArrList_pa    ArrList,
                            ArrListItem_u Item)
{
    Boolean_t IsOk = TRUE;

    REQUIRE(ArrListIsValid(ArrList));

    IsOk = ArrListInsertItem(ArrList, ArrList->Count, Item);

    ENSURE(ArrListIsValid(ArrList));
    ENSURE(VALID_BOOLEAN(IsOk));
    return IsOk;
}


/**
 * Appends copies of the items from the source list to the target list.
 * The source list remains unchanged.
 *
 * param Target
 *     Array list receiving the source items.
 * param Source
 *     Array list supplying the source items.
 *
 * return
 *     TRUE if sufficient memory permitted the operation, otherwise FALSE.
 */
Boolean_t ArrListAppend(ArrList_pa Target,
                        ArrList_pa Source)
{
    Boolean_t IsOk = TRUE;

    REQUIRE(ArrListIsValid(Target));
    REQUIRE(ArrListIsValid(Source));
    REQUIRE(Target != Source);
    REQUIRE(Target->Type == Source->Type);

    IsOk = ArrListInsert(Target, Target->Count, Source);

    ENSURE(ArrListIsValid(Target));
    ENSURE(VALID_BOOLEAN(IsOk));
    return IsOk;
}


/**
 * Copies the items of the array list.
 *
 * return
 *     Handle to a duplicate of the specified array list if sufficient
 *     memory permitted the operation, otherwise NULL.
 */
ArrList_pa ArrListCopy(ArrList_pa ArrList)
{
    ArrList_pa Result = NULL;

    REQUIRE(ArrListIsValid(ArrList));

    Result = ArrListGetItems(ArrList, 0, ArrList->Count);

    ENSURE(Result == NULL ||
           (ArrListIsValid(Result) && Result->Count == ArrList->Count));
    return Result;
}


/**
 * Creates a native 'C' array containing copies of the items held in the
 * source array list. The size of the native array is always one larger
 * than the Count and is assigned a value zero.
 *
 * param Source
 *     Array list containing the items of interest.
 * param Target
 *     Reference to a pointer capable of referencing a native 'C'
 *     array of the specified elements.
 * param Count
 *     Number of items returned in the native Target array not including
 *     the zero terminator.
 *
 * note
 *     The array returned in Target is allocated and it is the
 *     client's responsibility to release it when it is no longer
 *     needed. If sufficient
 *
 * return
 *     TRUE if sufficient memory permitted the operation, otherwise FALSE.
 */
Boolean_t ArrListToNative(ArrList_pa Source,
                          void       **Target,
                          LgIndex_t  *Count)
{
    Boolean_t IsOk = TRUE;

    REQUIRE(ArrListIsValid(Source));
    REQUIRE(VALID_REF(Target));

    *Target = (void *)ALLOC_ARRAY((Source->Count + 1) * Source->ItemSize,
                                  char, "native array");
    IsOk = (*Target != NULL);
    if (IsOk)
    {
        if (Source->Count > 0)
            CopyArrayItems((char *)*Target, 0,
                           Source->Array, 0,
                           Source->Count,
                           Source->ItemSize);

        /* zero terminate the end of the array */
        memset(&((char *)*Target)[Source->Count*Source->ItemSize],
               0, Source->ItemSize);

        *Count = Source->Count;
    }
    else
    {
        *Count = 0;
    }

    ENSURE(VALID_REF(*Target) || *Target == NULL);
    ENSURE(IMPLICATION(IsOk, *Count == Source->Count));
    ENSURE(VALID_BOOLEAN(IsOk));
    return IsOk;
}


/**
 * Creates an array list containing copies of the items held in the
 * native 'C' array.
 *
 * param Source
 *     Native 'C' array containing the items of interest.
 * param Count
 *     Number of items contained in the native 'C' array.
 * param Type
 *     Type of items contained in the native 'C' array.
 *
 * return
 *     Array list handle containing copies of the items held in the
 *     native 'C' array if sufficient memory was available, otherewise
 *     a handle to NULL.
 */
ArrList_pa ArrListFromNative(void         *Source,
                             LgIndex_t     Count,
                             ArrListType_e Type)
{
    ArrList_pa Result = NULL;

    REQUIRE(VALID_REF(Source));
    REQUIRE(Count >= 0);
    REQUIRE(VALID_ENUM(Type, ArrListType_e));

    Result = ArrListAlloc(Count, Type);
    if (Result != NULL && Count > 0)
    {
        CopyArrayItems(Result->Array, 0,
                       (char *)Source, 0,
                       Count,
                       Result->ItemSize);
        Result->Count = Count;
    }

    ENSURE(ArrListIsValid(Result) || Result == NULL);
    return Result;
}
