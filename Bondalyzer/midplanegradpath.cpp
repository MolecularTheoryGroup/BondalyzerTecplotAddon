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
#include "ALLOC.h"
#include "ARRLIST.h"
#include "ODERUNGEKUTTA.h"
#include "NORMALS.h"
#include "SURFELEMMAP.h"
#include "ZONEVARINFO.h"
#include "CRITPOINTS.h"
#include "GRADPATH.h"
#include "BUNDLES.h"
// #include "NORMALS.h"
#include "MIDPLANEGRADPATH.h"
#include "SURFACEFIT.h"
#include "ENGINE.h"



/**
 * Determine if the MidPlnGradPath handle is sane.
 *
 * param MidPlnGradPath
 *     MidPlnGradPath structure in question.
 *
 * return
 *     TRUE if the MidPlnGradPath structure is valid, otherwise FALSE.
 */
Boolean_t MidPlnGradPathIsValid(MidPlnGradPath_pa MidPlnGradPath)
{
    Boolean_t IsValid = FALSE;

    IsValid = (VALID_REF(MidPlnGradPath) &&
               VALID_REF(MidPlnGradPath->BBPIntersectX) && ArrListIsValid(MidPlnGradPath->BBPIntersectX) &&
               VALID_REF(MidPlnGradPath->BBPIntersectY) && ArrListIsValid(MidPlnGradPath->BBPIntersectY) &&
               VALID_REF(MidPlnGradPath->BBPIntersectZ) && ArrListIsValid(MidPlnGradPath->BBPIntersectZ) &&
               VALID_REF(MidPlnGradPath->SeedX) && ArrListIsValid(MidPlnGradPath->SeedX) &&
               VALID_REF(MidPlnGradPath->SeedY) && ArrListIsValid(MidPlnGradPath->SeedY) &&
               VALID_REF(MidPlnGradPath->SeedZ) && ArrListIsValid(MidPlnGradPath->SeedZ) &&
               VALID_REF(MidPlnGradPath->PathLength) && ArrListIsValid(MidPlnGradPath->PathLength) &&
               VALID_REF(MidPlnGradPath->GradPaths)   && ArrListIsValid(MidPlnGradPath->GradPaths));

    /* Require the same count for each coordinate array for the Bounding Path intersections. */
    if (IsValid) IsValid = (ArrListGetCount(MidPlnGradPath->BBPIntersectX) ==
                                ArrListGetCount(MidPlnGradPath->BBPIntersectZ) &&
                                ArrListGetCount(MidPlnGradPath->BBPIntersectY) ==
                                ArrListGetCount(MidPlnGradPath->BBPIntersectZ));
    /* Require the same count for each coordinate array. */
    if (IsValid) IsValid = (ArrListGetCount(MidPlnGradPath->SeedX) ==
                                ArrListGetCount(MidPlnGradPath->GradPaths) &&
                                ArrListGetCount(MidPlnGradPath->SeedY) ==
                                ArrListGetCount(MidPlnGradPath->GradPaths) &&
                                ArrListGetCount(MidPlnGradPath->SeedZ) ==
                                ArrListGetCount(MidPlnGradPath->GradPaths) &&
                                ArrListGetCount(MidPlnGradPath->PathLength) ==
                                ArrListGetCount(MidPlnGradPath->GradPaths));

    ENSURE(VALID_BOOLEAN(IsValid));
    return IsValid;
}





/**
 * Gets the Index'th GradPath handle.
 *
 *
 * param MidPlnGradPath
 *     MidPlnGradPath to get the GradPath from.
 *
 * param Offset
 *     Offset into the GradPath pointer ArrList.
 *
 * return
 *     GradPath handle available, otherewise a handle to NULL.
 */
GradPath_pa MidPlnGradPathGetGP(MidPlnGradPath_pa MidPlnGradPath,
                                LgIndex_t         Offset)
{
    GradPath_pa   Result = NULL;
    ArrListItem_u Item;

    REQUIRE(VALID_REF(MidPlnGradPath));
    REQUIRE(MidPlnGradPathIsValid(MidPlnGradPath));
    REQUIRE(Offset >= 0 && Offset < MidPlnGradPathGetGPCount(MidPlnGradPath));

    Item = ArrListGetItem(MidPlnGradPath->GradPaths, Offset);
    Result = (GradPath_pa)Item.VoidPtr;

    ENSURE(VALID_REF(Result));
    return Result;
}





/**
 * Gets the seed point position for the specified GridPath.
 *
 *
 * param MidPlnGradPath
 *     MidPlnGradPath from which to get the Angle.
 *
 * param Offset
 *     Offset into the Seed(XYZ) ArrLists.
 *
 * return
 *     XYZ position (structure).
 */
XYZ_s MidPlnGradPathGetSeedPos(MidPlnGradPath_pa MidPlnGradPath,
                               LgIndex_t         Offset)
{
    XYZ_s  Result;
    ArrListItem_u Item;

    REQUIRE(VALID_REF(MidPlnGradPath));
    REQUIRE(MidPlnGradPathIsValid(MidPlnGradPath));
    REQUIRE(Offset >= 0 && Offset < MidPlnGradPathGetGPCount(MidPlnGradPath));

    Item = ArrListGetItem(MidPlnGradPath->SeedX, Offset);
    Result.X = Item.Double;
    Item = ArrListGetItem(MidPlnGradPath->SeedY, Offset);
    Result.Y = Item.Double;
    Item = ArrListGetItem(MidPlnGradPath->SeedZ, Offset);
    Result.Z = Item.Double;

    return Result;
}






/**
 * Gets the bundle boundary path intersection with the mid-plane
 * for the specified offset.
 *
 *
 * param MidPlnGradPath
 *     MidPlnGradPath from which to get the Angle.
 *
 * param Offset
 *     Offset into the BBPIntersect(XYZ) ArrLists.
 *
 * return
 *     XYZ position (structure).
 */
XYZ_s MidPlnGradPathGetBBPIPos(MidPlnGradPath_pa MidPlnGradPath,
                               LgIndex_t         Offset)
{
    XYZ_s  Result;
    ArrListItem_u Item;

    REQUIRE(VALID_REF(MidPlnGradPath));
    REQUIRE(MidPlnGradPathIsValid(MidPlnGradPath));
    REQUIRE(Offset >= 0 && Offset < MidPlnGradPathGetBBPICount(MidPlnGradPath));

    Item = ArrListGetItem(MidPlnGradPath->BBPIntersectX, Offset);
    Result.X = Item.Double;
    Item = ArrListGetItem(MidPlnGradPath->BBPIntersectY, Offset);
    Result.Y = Item.Double;
    Item = ArrListGetItem(MidPlnGradPath->BBPIntersectZ, Offset);
    Result.Z = Item.Double;

    return Result;
}






/**
 * Gets the path length for the specified GridPath.
 *
 *
 * param MidPlnGradPath
 *     MidPlnGradPath from which to get the Angle.
 *
 * param Offset
 *     Offset into the Seed(XYZ) ArrLists.
 *
 * return
 *     Path length (double).
 */
double MidPlnGradPathGetPathLen(MidPlnGradPath_pa MidPlnGradPath,
                                LgIndex_t         Offset)
{
    double Result;
    ArrListItem_u Item;

    REQUIRE(VALID_REF(MidPlnGradPath));
    REQUIRE(MidPlnGradPathIsValid(MidPlnGradPath));
    REQUIRE(Offset >= 0 && Offset < MidPlnGradPathGetGPCount(MidPlnGradPath));

    Item = ArrListGetItem(MidPlnGradPath->PathLength, Offset);
    Result = Item.Double;

    return Result;
}







/**
 * Sets the relative start direction vector for the specified GridPath.
 *
 *
 * param MidPlnGradPath
 *     MidPlnGradPath from which to get the Angle.
 *
 * param Offset
 *     Offset into the Seed(XYZ) ArrList.
 *
 * param SeedPos
 *     Seed point position.
 *
 * return
 *     TRUE if it works, FALSE otherwise.
 */
Boolean_t MidPlnGradPathSetSeedPos(MidPlnGradPath_pa MidPlnGradPath,
                                   LgIndex_t         Offset,
                                   XYZ_s             SeedPos)
{
    Boolean_t  IsOk = TRUE;
    ArrListItem_u Item;

    REQUIRE(VALID_REF(MidPlnGradPath));
    REQUIRE(MidPlnGradPathIsValid(MidPlnGradPath));
    REQUIRE(Offset >= 0 && Offset < MidPlnGradPathGetGPCount(MidPlnGradPath));

    Item.Double = SeedPos.X;
    IsOk = ArrListSetItem(MidPlnGradPath->SeedX, Offset, Item);

    if (IsOk)
    {
        Item.Double = SeedPos.Y;
        IsOk = ArrListSetItem(MidPlnGradPath->SeedY, Offset, Item);
    }

    if (IsOk)
    {
        Item.Double = SeedPos.Z;
        IsOk = ArrListSetItem(MidPlnGradPath->SeedZ, Offset, Item);
    }

    ENSURE(VALID_BOOLEAN(IsOk));
    return IsOk;
}









/**
 * Sets the bundle boundary path intersection locations with the
 * mid-plane for the specified offset.
 *
 *
 * param MidPlnGradPath
 *     MidPlnGradPath for which to set the BBP intersection position.
 *
 * param Offset
 *     Offset into the BBPIntersect(XYZ) ArrList.
 *
 * param Position
 *     Position (XYZ) where the bundle boundary path intersects the
 *     mid-plane between the two critical points.
 *
 * return
 *     TRUE if it works, FALSE otherwise.
 */
Boolean_t MidPlnGradPathSetBBPIPos(MidPlnGradPath_pa MidPlnGradPath,
                                   LgIndex_t         Offset,
                                   XYZ_s             Position)
{
    Boolean_t  IsOk = TRUE;
    ArrListItem_u Item;

    REQUIRE(VALID_REF(MidPlnGradPath));
    REQUIRE(MidPlnGradPathIsValid(MidPlnGradPath));
    REQUIRE(Offset >= 0);

    Item.Double = Position.X;
    IsOk = ArrListSetItem(MidPlnGradPath->BBPIntersectX, Offset, Item);

    if (IsOk)
    {
        Item.Double = Position.Y;
        IsOk = ArrListSetItem(MidPlnGradPath->BBPIntersectY, Offset, Item);
    }

    if (IsOk)
    {
        Item.Double = Position.Z;
        IsOk = ArrListSetItem(MidPlnGradPath->BBPIntersectZ, Offset, Item);
    }

    ENSURE(VALID_BOOLEAN(IsOk));
    return IsOk;
}








/**
 * Appends the bundle boundary path intersection locations with the
 * mid-plane.
 *
 * param MidPlnGradPath
 *     MidPlnGradPath for which to append the BBP intersection position.
 *
 * param Position
 *     Position (XYZ) where the bundle boundary path intersects the
 *     mid-plane between the two critical points.
 *
 * return
 *     TRUE if sufficient memory permitted the operation, otherwise FALSE.
 */
Boolean_t MidPlnGradPathAppendBBPIPos(MidPlnGradPath_pa  MidPlnGradPath,
                                      XYZ_s              Position)
{
    Boolean_t IsOk = TRUE;
    LgIndex_t Count;

    REQUIRE(MidPlnGradPathIsValid(MidPlnGradPath));

    Count = MidPlnGradPathGetBBPICount(MidPlnGradPath);

    IsOk = MidPlnGradPathSetBBPIPos(MidPlnGradPath, Count, Position);

    ENSURE(MidPlnGradPathIsValid(MidPlnGradPath));

    ENSURE(VALID_BOOLEAN(IsOk));
    return IsOk;
}










/**
 * Deallocates the MidPlnGradPath handle and set the handle to NULL.
 *
 * note
 *     Item destruction is the responsibility of the caller.
 *
 * param
 *     Reference to a MidPlnGradPath handle.
 */
void MidPlnGradPathDealloc(MidPlnGradPath_pa *MidPlnGradPath)
{

    REQUIRE(VALID_REF(MidPlnGradPath));
    // REQUIRE(MidPlnGradPathIsValid(*MidPlnGradPath) || *MidPlnGradPath == NULL);

    if (*MidPlnGradPath != NULL)
    {
        LgIndex_t ii;
        LgIndex_t Count = ArrListGetCount((*MidPlnGradPath)->GradPaths);

        /* Dealloc the GradPaths */
        for (ii = 0; ii < Count; ii++)
        {
            GradPath_pa GradPath = MidPlnGradPathGetGP((*MidPlnGradPath), ii);
            GradPathDealloc(&GradPath);
        }


        /* release the ArrList's */
        if ((*MidPlnGradPath)->GradPaths != NULL) ArrListDealloc(&((*MidPlnGradPath)->GradPaths));
        if ((*MidPlnGradPath)->SeedX != NULL) ArrListDealloc(&((*MidPlnGradPath)->SeedX));
        if ((*MidPlnGradPath)->SeedY != NULL) ArrListDealloc(&((*MidPlnGradPath)->SeedY));
        if ((*MidPlnGradPath)->SeedZ != NULL) ArrListDealloc(&((*MidPlnGradPath)->SeedZ));
        if ((*MidPlnGradPath)->PathLength != NULL) ArrListDealloc(&((*MidPlnGradPath)->PathLength));
        if ((*MidPlnGradPath)->BBPIntersectX != NULL) ArrListDealloc(&((*MidPlnGradPath)->BBPIntersectX));
        if ((*MidPlnGradPath)->BBPIntersectY != NULL) ArrListDealloc(&((*MidPlnGradPath)->BBPIntersectY));
        if ((*MidPlnGradPath)->BBPIntersectZ != NULL) ArrListDealloc(&((*MidPlnGradPath)->BBPIntersectZ));

        /* release the list structure itself */
        FREE_ITEM(*MidPlnGradPath, "MidPlnGradPath structure");
        *MidPlnGradPath = NULL;
    }

    ENSURE(*MidPlnGradPath == NULL);
}




/**
 * Gets the number of GradPaths/RelNormSeedVector's currently in the MidPlnGradPath
 * (maintained by the MidPlnGradPath array lists).
 *
 * param
 *     MidPlnGradPath structure in question.
 *
 * return
 *     Number of GradPaths/RelNormSeedVector's in the MidPlnGradPath.
 */
LgIndex_t MidPlnGradPathGetGPCount(MidPlnGradPath_pa MidPlnGradPath)
{
    LgIndex_t Result = 0;

    REQUIRE(MidPlnGradPathIsValid(MidPlnGradPath));

    Result = ArrListGetCount(MidPlnGradPath->GradPaths);

    ENSURE(Result >= 0);
    return Result;
}







/**
 * Gets the number of bounding bundle path intersections (BBPIntersect[XYZ])
 * currently in the MidPlnGradPath (maintained by the MidPlnGradPath array lists).
 *
 * param
 *     MidPlnGradPath structure in question.
 *
 * return
 *     Number of BBPIntersect[XYZ] in the MidPlnGradPath.
 */
LgIndex_t MidPlnGradPathGetBBPICount(MidPlnGradPath_pa MidPlnGradPath)
{
    LgIndex_t Result = 0;

    REQUIRE(MidPlnGradPathIsValid(MidPlnGradPath));

    Result = ArrListGetCount(MidPlnGradPath->BBPIntersectX);

    ENSURE(Result >= 0);
    return Result;
}











/**
 * Empties the MidPlnGradPath of all Seeds, GradPaths, and BBPIntersects
 * and resets the other variables.
 *
 *
 * param MidPlnGradPath
 *     MidPlnGradPath to clear.
 */
void MidPlnGradPathClear(MidPlnGradPath_pa MidPlnGradPath)
{
    LgIndex_t ii;
    LgIndex_t Count = ArrListGetCount(MidPlnGradPath->GradPaths);

    REQUIRE(MidPlnGradPathIsValid(MidPlnGradPath));

    for (ii = 0; ii < Count; ii++)
    {
        GradPath_pa GradPath = MidPlnGradPathGetGP(MidPlnGradPath, ii);
        GradPathDealloc(&GradPath);
    }

    ArrListClear(MidPlnGradPath->GradPaths);
    ArrListClear(MidPlnGradPath->SeedX);
    ArrListClear(MidPlnGradPath->SeedY);
    ArrListClear(MidPlnGradPath->SeedZ);
    ArrListClear(MidPlnGradPath->PathLength);
    ArrListClear(MidPlnGradPath->BBPIntersectX);
    ArrListClear(MidPlnGradPath->BBPIntersectZ);
    ArrListClear(MidPlnGradPath->BBPIntersectZ);

    MidPlnGradPath->Position[0]   = 0.0;
    MidPlnGradPath->Position[1]   = 0.0;
    MidPlnGradPath->Position[2]   = 0.0;
    MidPlnGradPath->Normal[0]   = 0.0;
    MidPlnGradPath->Normal[1]   = 0.0;
    MidPlnGradPath->Normal[2]   = 0.0;

    ENSURE(MidPlnGradPathIsValid(MidPlnGradPath) &&
           MidPlnGradPathGetGPCount(MidPlnGradPath) == 0 &&
           MidPlnGradPathGetBBPICount(MidPlnGradPath) == 0);
}






/**
 * Allocates a MidPlnGradPath handle with a suitable default
 * capacity for the contained ArrList's.
 *
 *
 * return
 *     MidPlnGradPath handle if sufficient memory was available,
 *     otherewise a handle to NULL.
 */
MidPlnGradPath_pa MidPlnGradPathAlloc()
{
    MidPlnGradPath_pa Result = NULL;

    Result = ALLOC_ITEM(MidPlnGradPath_s, "MidPlnGradPath structure");
    if (Result != NULL)
    {
        Result->Position[0]   = 0.0;
        Result->Position[1]   = 0.0;
        Result->Position[2]   = 0.0;
        Result->Normal[0]     = 1.0;
        Result->Normal[1]     = 1.0;
        Result->Normal[2]     = 1.0;
        Result->BBPIntersectX = ArrListAlloc(20, ArrListType_Double);
        Result->BBPIntersectY = ArrListAlloc(20, ArrListType_Double);
        Result->BBPIntersectZ = ArrListAlloc(20, ArrListType_Double);
        Result->SeedX = ArrListAlloc(20, ArrListType_Double);
        Result->SeedY = ArrListAlloc(20, ArrListType_Double);
        Result->SeedZ = ArrListAlloc(20, ArrListType_Double);
        Result->GradPaths   = ArrListAlloc(20, ArrListType_VoidPtr);
        Result->PathLength  = ArrListAlloc(20, ArrListType_Double);

        /* If it failed to allocate any of the array lists, clean-up and exit. */
        if (Result->BBPIntersectX == NULL || Result->BBPIntersectY == NULL ||
            Result->BBPIntersectZ == NULL ||
            Result->SeedX == NULL || Result->SeedY == NULL ||
            Result->SeedZ == NULL || Result->GradPaths == NULL)
        {
            if (Result->BBPIntersectX != NULL) ArrListDealloc(&(Result->BBPIntersectX));
            if (Result->BBPIntersectY != NULL) ArrListDealloc(&(Result->BBPIntersectY));
            if (Result->BBPIntersectZ != NULL) ArrListDealloc(&(Result->BBPIntersectZ));
            if (Result->SeedX != NULL) ArrListDealloc(&(Result->SeedX));
            if (Result->SeedY != NULL) ArrListDealloc(&(Result->SeedY));
            if (Result->SeedZ != NULL) ArrListDealloc(&(Result->SeedZ));
            if (Result->PathLength != NULL) ArrListDealloc(&(Result->PathLength));
            if (Result->GradPaths != NULL) ArrListDealloc(&(Result->GradPaths));
            FREE_ITEM(Result, "MidPlnGradPath structure");
            Result = NULL;
        }

    }

    ENSURE(MidPlnGradPathIsValid(Result) || Result == NULL);
    return Result;
}






/**
 * Places GradPath pointer and SeedPosition at the specified offset. If the
 * offset is beyond the end of the array lists, they are sized
 * accordingly and the intervening array values between the last
 * node of the original state and the last node of the new state are
 * assigned 0.
 *
 * note
 *     If SeedVector/GradPaths already exists at the specified location they are
 *     replaced. Therefore, item destruction is the responsibility of
 *     the caller.
 *
 * param MidPlnGradPath
 *     MidPlnGradPath target in which to set the Angle/GradPath.
 * param Offset
 *     Offset into array lists.
 * param SeedVector
 *     Angle to set at the specified offset.
 * param GradPath
 *     GradPath handle to set at the specified offset.
 *
 * return
 *     TRUE if sufficient memory permitted the operation, otherwise FALSE.
 */

Boolean_t   MidPlnGradPathSetGP(MidPlnGradPath_pa MidPlnGradPath,
                                LgIndex_t   Offset,
                                XYZ_s       SeedPos,
                                double      PathLength,
                                GradPath_pa GradPath)
{
    Boolean_t IsOk = TRUE;
    ArrListItem_u Item;


    REQUIRE(MidPlnGradPathIsValid(MidPlnGradPath));
    REQUIRE(Offset >= 0);

    Item.Double = SeedPos.X;
    IsOk = ArrListSetItem(MidPlnGradPath->SeedX, Offset, Item);

    if (IsOk)
    {
        Item.Double = SeedPos.Y;
        IsOk = ArrListSetItem(MidPlnGradPath->SeedY, Offset, Item);
    }

    if (IsOk)
    {
        Item.Double = SeedPos.Z;
        IsOk = ArrListSetItem(MidPlnGradPath->SeedZ, Offset, Item);
    }

    if (IsOk)
    {
        Item.Double = PathLength;
        IsOk = ArrListSetItem(MidPlnGradPath->PathLength, Offset, Item);
    }

    if (IsOk)
    {
        Item.VoidPtr = (void *)GradPath;
        IsOk = ArrListSetItem(MidPlnGradPath->GradPaths, Offset, Item);
    }

    ENSURE(MidPlnGradPathIsValid(MidPlnGradPath));

    ENSURE(VALID_BOOLEAN(IsOk));
    return IsOk;
}






/**
 * Inserts Seed position/GradPath at the specified offset. The arrays
 * will be expanded to accomodate the additional value.
 *
 * note
 *     If SeedVector/GradPath already exists at the specified location they are
 *     replaced. Therefore, item destruction is the responsibility of
 *     the caller.
 *
 * param MidPlnGradPath
 *     MidPlnGradPath target in which to set the coordinates.
 * param Offset
 *     Offset at which to insert the Angle/GradPath.
 * param SeedPos
 *     XYZ seed point position (on mid-plane) to insert at the specified offset.
 * param GradPath
 *     GradPath handle to insert at the specified offset.
 *
 * return
 *     TRUE if sufficient memory permitted the operation, otherwise FALSE.
 */

Boolean_t MidPlnGradPathInsertGP(MidPlnGradPath_pa MidPlnGradPath,
                                 LgIndex_t         Offset,
                                 XYZ_s             SeedPos,
                                 double            PathLength,
                                 GradPath_pa       GradPath)
{
    Boolean_t IsOk = TRUE;
    ArrListItem_u Item;


    REQUIRE(MidPlnGradPathIsValid(MidPlnGradPath));
    REQUIRE(0 <= Offset && Offset <= MidPlnGradPathGetGPCount(MidPlnGradPath));

    Item.Double = SeedPos.X;
    IsOk = ArrListInsertItem(MidPlnGradPath->SeedX, Offset, Item);

    if (IsOk)
    {
        Item.Double = SeedPos.Y;
        IsOk = ArrListInsertItem(MidPlnGradPath->SeedY, Offset, Item);
    }

    if (IsOk)
    {
        Item.Double = SeedPos.Z;
        IsOk = ArrListInsertItem(MidPlnGradPath->SeedZ, Offset, Item);
    }

    if (IsOk)
    {
        Item.Double = PathLength;
        IsOk = ArrListInsertItem(MidPlnGradPath->PathLength, Offset, Item);
    }

    if (IsOk)
    {
        Item.VoidPtr = (void *)GradPath;
        IsOk = ArrListInsertItem(MidPlnGradPath->GradPaths, Offset, Item);
    }

    ENSURE(MidPlnGradPathIsValid(MidPlnGradPath));

    ENSURE(VALID_BOOLEAN(IsOk));
    return IsOk;
}






/**
 * Deletes Seed Position/GradPath at the specified offset.
 *
 *
 * param MidPlnGradPath
 *     MidPlnGradPath target in which to set the coordinates.
 * param Offset
 *     Offset at which to insert the SeedVector/GradPath.
 *
 * return
 *     TRUE if successful operation, otherwise FALSE.
 */

Boolean_t MidPlnGradPathRemoveGP(MidPlnGradPath_pa MidPlnGradPath,
                                 LgIndex_t         Offset)
{
    Boolean_t IsOk = TRUE;
    ArrListItem_u Item;


    REQUIRE(MidPlnGradPathIsValid(MidPlnGradPath));
    REQUIRE(0 <= Offset && Offset <= MidPlnGradPathGetGPCount(MidPlnGradPath));

    Item = ArrListRemoveItem(MidPlnGradPath->SeedX, Offset);

    if (IsOk)
        Item = ArrListRemoveItem(MidPlnGradPath->SeedY, Offset);

    if (IsOk)
        Item = ArrListRemoveItem(MidPlnGradPath->SeedZ, Offset);

    if (IsOk)
        Item = ArrListRemoveItem(MidPlnGradPath->PathLength, Offset);

    if (IsOk)
    {
        GradPath_pa GradPath;
        Item = ArrListRemoveItem(MidPlnGradPath->GradPaths, Offset);
        GradPath = (GradPath_pa)Item.VoidPtr;
        GradPathDealloc(&GradPath);
    }

    ENSURE(MidPlnGradPathIsValid(MidPlnGradPath));

    ENSURE(VALID_BOOLEAN(IsOk));
    return IsOk;
}






/**
 * Appends the Seed Position/GradPath to the array lists in CircleGradPath. The
 * array lists will be expanded to accommodate the additional items.
 *
 * param MidPlnGradPath
 *     MidPlnGradPath target to which the SeedVector/GradPath is to be appended.
 * param SeedPos
 *     Seed point position (on mid-plane) to append to the SeedVector/GradPath.
 * param GradPath
 *     GradPath handle to append to the MidPlnGradPath.
 *
 * return
 *     TRUE if sufficient memory permitted the operation, otherwise FALSE.
 */
Boolean_t MidPlnGradPathAppendGP(MidPlnGradPath_pa  MidPlnGradPath,
                                 XYZ_s              SeedPos,
                                 double             PathLength,
                                 GradPath_pa        GradPath)
{
    Boolean_t IsOk = TRUE;
    LgIndex_t Count;

    REQUIRE(MidPlnGradPathIsValid(MidPlnGradPath));

    Count = MidPlnGradPathGetGPCount(MidPlnGradPath);

    IsOk = MidPlnGradPathInsertGP(MidPlnGradPath, Count, SeedPos, PathLength, GradPath);

    ENSURE(MidPlnGradPathIsValid(MidPlnGradPath));

    ENSURE(VALID_BOOLEAN(IsOk));
    return IsOk;
}










/**
 * Return the ending critical point number for the GradPath at Offset
 * in the array list.
 *
 * param MidPlnGradPath
 *     MidPlnGradPath target to extract the EndCrtPtNum from.
 * param Offset
 *     Offset into array list of GradPath from which to extract EndCrtPtNum.
 *
 * return
 *     EndCrtPtNum for the path.
 */
LgIndex_t MidPlnGradPathGetEndCPNum(MidPlnGradPath_pa  MidPlnGradPath,
                                    LgIndex_t          Offset)
{
    LgIndex_t EndCrtPtNum = -2;

    GradPath_pa GradPath = NULL;
    ArrListItem_u Item;

    REQUIRE(MidPlnGradPathIsValid(MidPlnGradPath));
    REQUIRE(0 <= Offset && Offset < MidPlnGradPathGetGPCount(MidPlnGradPath));

    Item = ArrListGetItem(MidPlnGradPath->GradPaths, Offset);
    GradPath = (GradPath_pa)Item.VoidPtr;

    EndCrtPtNum = GradPath->EndCrtPtNum;

    ENSURE(-1 <= EndCrtPtNum);
    return EndCrtPtNum;
}




/**
 * Return the begining critical point number for the MidPlnGradPath.
 *
 * param MidPlnGradPath
 *     MidPlnGradPath target to extract the EndCrtPtNum from.
 *
 * return
 *     BegCrtPtNum for the path.
 */
LgIndex_t MidPlnGradPathGetBegCPNum(MidPlnGradPath_pa  MidPlnGradPath)
{
    LgIndex_t BegCrtPtNum = -2;
    LgIndex_t Offset = 0;

    GradPath_pa GradPath = NULL;
    ArrListItem_u Item;

    REQUIRE(MidPlnGradPathIsValid(MidPlnGradPath));

    Item = ArrListGetItem(MidPlnGradPath->GradPaths, Offset);
    GradPath = (GradPath_pa)Item.VoidPtr;

    BegCrtPtNum = GradPath->BeginCrtPtNum;

    ENSURE(-1 <= BegCrtPtNum);
    return BegCrtPtNum;
}







/**
 * Return the Offset into MidPlnGradPath->GradPaths of the GradPath with
 * the specified seed point.
 *
 * param MidPlnGradPath
 *     MidPlnGradPath target to extract the EndCrtPtNum from.
 * param Seed
 *     Seed location (compare to MidPlnGradPath->SeedX, etc.).
 *
 * return
 *     Offset of desired GradPath in MidPlnGradPath->GradPaths list, or -1 if fails.
 */
LgIndex_t MidPlnGradPathGetGPBySeed(const MidPlnGradPath_pa  MidPlnGradPath,
                                    const XYZ_s              Seed)
{
    LgIndex_t GradPathOffset = -1;
    LgIndex_t Offset;
    Boolean_t Found = FALSE;
    LgIndex_t NumGradPaths = MidPlnGradPathGetGPCount(MidPlnGradPath);
    double    Tolerance = 1.0e-5;

    REQUIRE(MidPlnGradPathIsValid(MidPlnGradPath));

    for (Offset = 0; !Found && Offset < NumGradPaths; Offset++)
    {
        XYZ_s  GPSeed = MidPlnGradPathGetSeedPos(MidPlnGradPath, Offset);

        if (ABS(GPSeed.X - Seed.X) < Tolerance)
        {
            if (ABS(GPSeed.Y - Seed.Y) < Tolerance)
            {
                if (ABS(GPSeed.Z - Seed.Z) < Tolerance)
                {
                    Found = TRUE;
                    GradPathOffset = Offset;
                }
            }
        }
    }

    ENSURE(-1 <= GradPathOffset);
    return GradPathOffset;
}







/**
 * Return the XYZ location of the intersection of GradPath with MidPlnGradPath.
 *
 * param MidPlnGradPath
 *     MidPlnGradPath plane.
 * param GradPath
 *     GradPath to check for intesection with MidPlnGradPath plane.
 * params X, Y, Z
 *     Pointers to position coordinates for the intersection.
 *
 * return
 *     TRUE if intersection was found, FALSE otherwise.
 */
Boolean_t         MidPlnGradPathGetGPIntersect(const MidPlnGradPath_pa MidPlnGradPath,
                                               const GradPath_pa       GradPath,
                                               double                 *X,
                                               double                 *Y,
                                               double                 *Z)
{
    Boolean_t Found = FALSE;
    Boolean_t IsOk = TRUE;
    LgIndex_t SegNum;
    LgIndex_t GradPathCount, NumGradPathSegs;

    REQUIRE(MidPlnGradPathIsValid(MidPlnGradPath));
    REQUIRE(GradPathIsValid(GradPath));
    REQUIRE(VALID_REF(X));
    REQUIRE(VALID_REF(Y));
    REQUIRE(VALID_REF(Z));

    GradPathCount = GradPathGetCount(GradPath);
    NumGradPathSegs = GradPathCount - 1;

    if (NumGradPathSegs > 0)
    {
        double Xi, Yi, Zi, DXiDotN, junk;
        double XPln = MidPlnGradPath->Position[0];
        double YPln = MidPlnGradPath->Position[1];
        double ZPln = MidPlnGradPath->Position[2];
        double NxPln = MidPlnGradPath->Normal[0];
        double NyPln = MidPlnGradPath->Normal[1];
        double NzPln = MidPlnGradPath->Normal[2];

        IsOk = GradPathGetPoint(GradPath, 0, &Xi, &Yi, &Zi, &junk);

        // Distance from plane, in normal direction, of first point in GradPath
        DXiDotN = (Xi - XPln) * NxPln + (Yi - YPln) * NyPln + (Zi - ZPln) * NzPln;

        for (SegNum = 0; IsOk && !Found && SegNum < NumGradPathSegs; SegNum++)
        {
            double Xip1, Yip1, Zip1, DXip1DotN;

            IsOk = GradPathGetPoint(GradPath, SegNum + 1, &Xip1, &Yip1, &Zip1, &junk);

            // Distance from plane, in normal direction, of i+1 point in GradPath
            DXip1DotN = (Xip1 - XPln) * NxPln + (Yip1 - YPln) * NyPln + (Zip1 - ZPln) * NzPln;

            // Look for first segment where (X-Xplane)*N changes sign for i and ip1
            if (DXiDotN * DXip1DotN < 0.0)
            {
                double Fraction = ABS(DXiDotN) / (ABS(DXiDotN) + ABS(DXip1DotN));
                *X = Xi + Fraction * (Xip1 - Xi);
                *Y = Yi + Fraction * (Yip1 - Yi);
                *Z = Zi + Fraction * (Zip1 - Zi);
                Found = TRUE;
            }
            else
            {
                // Reset "i" values for next segment
                Xi = Xip1;
                Yi = Yip1;
                Zi = Zip1;
                DXiDotN = DXip1DotN;
            }
        }
    }

    ENSURE(VALID_BOOLEAN(Found));
    return Found;
}





/**
 * Search the Long ArrList to see if the specified item is in the list.
 *
 * param LongArrList
 *     Array list (Long type) to be tested.
 * param TestItem
 *     LgIndex_t to test against items in LongArrList.
 *
 * return
 *     TRUE if TestItem is found in LongArrList.
 */
Boolean_t InLongArrList(ArrList_pa LongArrList,
                        LgIndex_t  TestItem)
{
    Boolean_t Found = FALSE;
    LgIndex_t Offset;
    LgIndex_t NumItems = ArrListGetCount(LongArrList);

    REQUIRE(ArrListIsValid(LongArrList));
    REQUIRE(ArrListGetType(LongArrList) == ArrListType_Long);

    for (Offset = 0; !Found && Offset < NumItems; Offset++)
    {

        ArrListItem_u Item = ArrListGetItem(LongArrList, Offset);
        if (TestItem == Item.Long) Found = TRUE;
    }

    ENSURE(VALID_BOOLEAN(Found));
    return Found;
}








/*
 * MidPlnGradPathAdd: Given two critical points, compute the mid-plane position
 *  and orientation. Also compute the intersections of the bundle bounding
 *  gradpaths with the mid-plane. Finally, compute a set of gradpaths that are
 *  within the bundle-bounding gradpaths (therefore should connect the two
 *  critical points), one of which is seeded at the average of all the BBP
 *  intersections. This function sets the initial state of the Mid-Plane Gradient
 *  path object. Iteration to find the minimum-length gradient path connection
 *  is done elsewhere.
 *
 * param ZoneVarInfo
 *     Zone and variable information.
 * param BeginCrtPtNum
 *     Atom or Cage critical point that is (hopefully) the source of all grad
 *     paths in this MidPlnGradPath data structure
 * param EndCrtPtNum
 *     Atom or Cage critical point that is (hopefully) the destination of all
 *     grad paths in this MidPlnGradPath data structure
 * param CritPoints
 *     Critical Points structure
 * param Bundles
 *     Lists containing the Atom, Bond, Ring, Cage combinations that make up the
 *     irreducible bundles.
 * param CPTolerance
 *     Tolerance - distance from critical point at which the pathline is said
 *     to terminate at the critical point.
 * param MidPlnGradPath
 *     MidPlnGradPath data structure to populate for an Atom or Cage critical point.
 *
 * return
 *     TRUE if successful, FALSE if there were errors.
 */
Boolean_t MidPlnGradPathAdd(const ZoneVarInfo_pa  ZoneVarInfo,
                            const LgIndex_t       BeginCrtPtNum,
                            const LgIndex_t       EndCrtPtNum,
                            const CritPoints_pa   CritPoints,
                            const Bundles_pa      Bundles,
                            const float           CPTolerance,
                            MidPlnGradPath_pa     MidPlnGradPath)
{
    Boolean_t IsOk = TRUE;
    double XCrtPtBeg, YCrtPtBeg, ZCrtPtBeg, XCrtPtEnd, YCrtPtEnd, ZCrtPtEnd, dummy;
    char   cdummy;
    EntIndex_t SourceZoneNum = CritPoints->SourceZoneNum;

    REQUIRE(VALID_REF(ZoneVarInfo));
    REQUIRE(CritPointsIsValid(CritPoints));
    REQUIRE(BeginCrtPtNum >= 0 && BeginCrtPtNum < CritPointsGetCount(CritPoints));
    REQUIRE(EndCrtPtNum >= 0 && EndCrtPtNum < CritPointsGetCount(CritPoints));
    REQUIRE(BundlesIsValid(Bundles));
    REQUIRE(CPTolerance > 0.0);
    REQUIRE(MidPlnGradPathIsValid(MidPlnGradPath));
    // Begining and ending critical points must be Cage and Atom
    REQUIRE((CritPointsGetType(CritPoints, BeginCrtPtNum) == -3 &&
             CritPointsGetType(CritPoints, EndCrtPtNum) == 3) ||
            (CritPointsGetType(CritPoints, BeginCrtPtNum) == 3 &&
             CritPointsGetType(CritPoints, EndCrtPtNum) == -3));

    // Find the mid-plane position and orientation (normal vector)
    IsOk = CritPointsGetPoint(CritPoints, BeginCrtPtNum,
                              &XCrtPtBeg, &YCrtPtBeg, &ZCrtPtBeg, &dummy, &cdummy,
                              &dummy, &dummy, &dummy);
    if (IsOk)
        IsOk = CritPointsGetPoint(CritPoints, EndCrtPtNum,
                                  &XCrtPtEnd, &YCrtPtEnd, &ZCrtPtEnd, &dummy, &cdummy,
                                  &dummy, &dummy, &dummy);

    if (IsOk)
    {
        double n0, n1, n2, Length;
        EntIndex_t SourceZoneNum = CritPoints->SourceZoneNum;

        // Position is average of the beginning and ending critical point positions
        MidPlnGradPath->Position[0] = 0.5 * (XCrtPtBeg + XCrtPtEnd);
        MidPlnGradPath->Position[1] = 0.5 * (YCrtPtBeg + YCrtPtEnd);
        MidPlnGradPath->Position[2] = 0.5 * (ZCrtPtBeg + ZCrtPtEnd);

        // Normal is normalized direction from beginning to ending critical points
        n0 = XCrtPtEnd - XCrtPtBeg;
        n1 = YCrtPtEnd - YCrtPtBeg;
        n2 = ZCrtPtEnd - ZCrtPtBeg;
        Length = sqrt(n0 * n0 + n1 * n1 + n2 * n2);
        MidPlnGradPath->Normal[0] = n0 / Length;
        MidPlnGradPath->Normal[1] = n1 / Length;
        MidPlnGradPath->Normal[2] = n2 / Length;
    }

    // Compute GradPath from plane position point. If it connects the Atom & Cage,
    // compute GradPaths for a constellation around that point.
    if (IsOk)
    {
        GradPath_pa GradPath = GradPathAlloc();
        LgIndex_t NumPathPoints = 0;
        double PathLength;
        XYZ_s SeedPos;
        SeedPos.X = MidPlnGradPath->Position[0];
        SeedPos.Y = MidPlnGradPath->Position[1];
        SeedPos.Z = MidPlnGradPath->Position[2];

        // Initial GradPath: Seed at point that straight line between
        // atom/cage intersects the MidPlane.
        GradPath = GradPathAddMidway(ZoneVarInfo, MAXGRADPATHPOINTS, CritPoints,
                                     SeedPos, 1.0, &NumPathPoints);
        if (GradPath == NULL) IsOk = FALSE;

        // Find the arclength of the gradpath
        PathLength = GradPathGetLength(GradPath);

        // Add initial GradPath to MidPlnGradPath
        if (IsOk)
            IsOk = MidPlnGradPathAppendGP(MidPlnGradPath, SeedPos, PathLength, GradPath);
    }



    // If the simple seed point didn't work (went to wrong critical points), find
    // and improved initial seed point using the intersections of the bundle-bounding
    // gradient paths with the mid-plane.
    if (IsOk)
    {
        GradPath_pa GradPath = MidPlnGradPathGetGP(MidPlnGradPath, 0);

        if (GradPath->BeginCrtPtNum != BeginCrtPtNum || GradPath->EndCrtPtNum != EndCrtPtNum)
        {
            char BegCrtPtType = CritPointsGetType(CritPoints, BeginCrtPtNum);
            char EndCrtPtType = CritPointsGetType(CritPoints, EndCrtPtNum);
            LgIndex_t Atom, Cage;

            Bundles_pa BundlesAtomCage = BundlesAlloc();

            if (BegCrtPtType == -3)
            {
                Atom = BeginCrtPtNum;
                Cage = EndCrtPtNum - CritPointsGetBegOffset(CritPoints, 3);
            }
            else
            {
                Atom = EndCrtPtNum;
                Cage = BeginCrtPtNum - CritPointsGetBegOffset(CritPoints, 3);
            }

            IsOk = BundlesGetBondRingFromAtomCage(Bundles, Atom, Cage, BundlesAtomCage);


            if (IsOk)
            {
                LgIndex_t bb;
                LgIndex_t NumBundles = BundlesGetCount(BundlesAtomCage);
                ArrList_pa PathTPZoneList = ArrListAlloc(20, ArrListType_Long);

                for (bb = 0; bb < NumBundles; bb++)
                {
                    LgIndex_t Atom, Bond, Ring, Cage;
                    LgIndex_t BondCrtPtNum, RingCrtPtNum, CageCrtPtNum;

                    IsOk = BundlesGetFromOffset(BundlesAtomCage, bb, &Atom, &Bond, &Ring, &Cage);

                    BondCrtPtNum = Bond + CritPointsGetBegOffset(CritPoints, -1);
                    RingCrtPtNum = Ring + CritPointsGetBegOffset(CritPoints,  1);
                    CageCrtPtNum = Cage + CritPointsGetBegOffset(CritPoints,  3);

                    // Handle Bond (if not already used)
                    if (IsOk)
                    {
                        LgIndex_t PathTPZoneNum = (LgIndex_t)GradPathTPZoneFromBegEndCP(BondCrtPtNum, Atom, SourceZoneNum);

                        if (PathTPZoneNum > 0 && !InLongArrList(PathTPZoneList, PathTPZoneNum))
                        {
                            Boolean_t Found;
                            ArrListItem_u Item;
                            GradPath_pa BondGradPath = NULL;
                            double X, Y, Z;

                            // Add Bond-Atom GradPath zone number to PathTPZoneList
                            Item.Long = PathTPZoneNum;
                            IsOk = ArrListAppendItem(PathTPZoneList, Item);

                            // Compute intersection with plane and add to BBPIntersect[XYZ]
                            BondGradPath = GradPathGetByBegEndCP(BondCrtPtNum, Atom);
                            if (BondGradPath == NULL) IsOk = FALSE;

                            if (IsOk)
                                Found = MidPlnGradPathGetGPIntersect(MidPlnGradPath, BondGradPath, &X, &Y, &Z);

                            if (IsOk && Found)
                            {
                                XYZ_s  Intersection;
                                Intersection.X = X;
                                Intersection.Y = Y;
                                Intersection.Z = Z;
                                IsOk = MidPlnGradPathAppendBBPIPos(MidPlnGradPath, Intersection);
                            }

                            // Cleanup
                            GradPathDealloc(&BondGradPath);
                        }
                    }

                    if (IsOk)
                    {
                        LgIndex_t PathTPZoneNum = (LgIndex_t)GradPathTPZoneFromBegEndCP(BondCrtPtNum, CageCrtPtNum, SourceZoneNum);

                        if (PathTPZoneNum > 0 && !InLongArrList(PathTPZoneList, PathTPZoneNum))
                        {
                            Boolean_t Found;
                            ArrListItem_u Item;
                            GradPath_pa BondGradPath = NULL;
                            double X, Y, Z;

                            // Add Bond-Cage GradPath zone number to PathTPZoneList
                            Item.Long = PathTPZoneNum;
                            IsOk = ArrListAppendItem(PathTPZoneList, Item);

                            // Compute intersection with plane and add to BBPIntersect[XYZ]
                            BondGradPath = GradPathGetByBegEndCP(BondCrtPtNum, CageCrtPtNum);
                            if (BondGradPath == NULL) IsOk = FALSE;

                            if (IsOk)
                                Found = MidPlnGradPathGetGPIntersect(MidPlnGradPath, BondGradPath, &X, &Y, &Z);

                            if (IsOk && Found)
                            {
                                XYZ_s  Intersection;
                                Intersection.X = X;
                                Intersection.Y = Y;
                                Intersection.Z = Z;
                                IsOk = MidPlnGradPathAppendBBPIPos(MidPlnGradPath, Intersection);
                            }

                            // Cleanup
                            GradPathDealloc(&BondGradPath);
                        }
                    }

                    // Handle Ring (if not already used)
                    if (IsOk)
                    {
                        LgIndex_t PathTPZoneNum = (LgIndex_t)GradPathTPZoneFromBegEndCP(RingCrtPtNum, Atom, SourceZoneNum);

                        if (PathTPZoneNum > 0 && !InLongArrList(PathTPZoneList, PathTPZoneNum))
                        {
                            Boolean_t Found;
                            ArrListItem_u Item;
                            GradPath_pa RingGradPath = NULL;
                            double X, Y, Z;

                            // Add Ring-Atom GradPath zone number to PathTPZoneList
                            Item.Long = PathTPZoneNum;
                            IsOk = ArrListAppendItem(PathTPZoneList, Item);

                            // Compute intersection with plane and add to BBPIntersect[XYZ]
                            RingGradPath = GradPathGetByBegEndCP(RingCrtPtNum, Atom);
                            if (RingGradPath == NULL) IsOk = FALSE;

                            if (IsOk)
                                Found = MidPlnGradPathGetGPIntersect(MidPlnGradPath, RingGradPath, &X, &Y, &Z);

                            if (IsOk && Found)
                            {
                                XYZ_s  Intersection;
                                Intersection.X = X;
                                Intersection.Y = Y;
                                Intersection.Z = Z;
                                IsOk = MidPlnGradPathAppendBBPIPos(MidPlnGradPath, Intersection);
                            }

                            // Cleanup
                            GradPathDealloc(&RingGradPath);
                        }
                    }

                    if (IsOk)
                    {
                        LgIndex_t PathTPZoneNum = (LgIndex_t)GradPathTPZoneFromBegEndCP(RingCrtPtNum, CageCrtPtNum, SourceZoneNum);

                        if (PathTPZoneNum > 0 && !InLongArrList(PathTPZoneList, PathTPZoneNum))
                        {
                            Boolean_t Found;
                            ArrListItem_u Item;
                            GradPath_pa RingGradPath = NULL;
                            double X, Y, Z;

                            // Add Ring-Atom GradPath zone number to PathTPZoneList
                            Item.Long = PathTPZoneNum;
                            IsOk = ArrListAppendItem(PathTPZoneList, Item);

                            // Compute intersection with plane and add to BBPIntersect[XYZ]
                            RingGradPath = GradPathGetByBegEndCP(RingCrtPtNum, CageCrtPtNum);
                            if (RingGradPath == NULL) IsOk = FALSE;

                            if (IsOk)
                                Found = MidPlnGradPathGetGPIntersect(MidPlnGradPath, RingGradPath, &X, &Y, &Z);

                            if (IsOk && Found)
                            {
                                XYZ_s  Intersection;
                                Intersection.X = X;
                                Intersection.Y = Y;
                                Intersection.Z = Z;
                                IsOk = MidPlnGradPathAppendBBPIPos(MidPlnGradPath, Intersection);
                            }

                            // Cleanup
                            GradPathDealloc(&RingGradPath);
                        }
                    }
                }

                ArrListDealloc(&PathTPZoneList);
            }

            // Create the new initial seed gradpath from an average of the intersections
            if (IsOk)
            {
                LgIndex_t bbpi;
                GradPath_pa GradPath = GradPathAlloc();
                LgIndex_t NumPathPoints = 0;
                LgIndex_t NumBBPIPos = MidPlnGradPathGetBBPICount(MidPlnGradPath);
                XYZ_s SeedPos;

                SeedPos.X = 0.0;
                SeedPos.Y = 0.0;
                SeedPos.Z = 0.0;
                for (bbpi = 0; IsOk && bbpi < NumBBPIPos; bbpi++)
                {
                    XYZ_s BBPISeedPos = MidPlnGradPathGetBBPIPos(MidPlnGradPath, bbpi);

                    SeedPos.X += BBPISeedPos.X;
                    SeedPos.Y += BBPISeedPos.Y;
                    SeedPos.Z += BBPISeedPos.Z;
                }
                SeedPos.X /= NumBBPIPos;
                SeedPos.Y /= NumBBPIPos;
                SeedPos.Z /= NumBBPIPos;

                // Initial GradPath: Seed at point that straight line between
                // atom/cage intersects the MidPlane.
                GradPath = GradPathAddMidway(ZoneVarInfo, MAXGRADPATHPOINTS, CritPoints,
                                             SeedPos, 1.0, &NumPathPoints);
                if (GradPath == NULL) IsOk = FALSE;

                // Add initial GradPath to MidPlnGradPath
                if (IsOk)
                {
                    double PathLength;
                    GradPath_pa  GradPathOld = MidPlnGradPathGetGP(MidPlnGradPath, 0);
                    GradPathDealloc(&GradPathOld);

                    // Find the arclength of the gradpath
                    PathLength = GradPathGetLength(GradPath);

                    IsOk = MidPlnGradPathSetGP(MidPlnGradPath, 0, SeedPos, PathLength, GradPath);
                }
            }
        }
    }


    ENSURE(VALID_BOOLEAN(IsOk));
    return IsOk;
}





/*
 * MidPlnGradPathMinimizeLen: Find the mid-plane seed point that generates the
 *  minimum length gradient path between the specified BeginCrtPtNum and EndCrtPtNum.
 *  Then generate the minimum length GradPath and return the offset to it.
 *  An existing MidPlnGraPath_pa structure, generated by a successful call to
 *
 * param ZoneVarInfo
 *     Zone and variable information.
 * param BeginCrtPtNum
 *     Atom or Cage critical point that is (hopefully) the source of all grad
 *     paths in this MidPlnGradPath data structure
 * param EndCrtPtNum
 *     Atom or Cage critical point that is (hopefully) the destination of all
 *     grad paths in this MidPlnGradPath data structure
 * param CritPoints
 *     Critical Points structure
 * param Bundles
 *     Lists containing the Atom, Bond, Ring, Cage combinations that make up the
 *     irreducible bundles.
 * param CPTolerance
 *     Tolerance - distance from critical point at which the pathline is said
 *     to terminate at the critical point.
 * param MidPlnGradPath
 *     MidPlnGradPath data structure to populate for an Atom or Cage critical point.
 *
 * return
 *     Offset in MidPlnGradPath of minimum length GradPath number if successful,
 *     -1 if there were errors.
 */
LgIndex_t MidPlnGradPathMinimizeLen(const ZoneVarInfo_pa  ZoneVarInfo,
                                    const LgIndex_t       BeginCrtPtNum,
                                    const LgIndex_t       EndCrtPtNum,
                                    const CritPoints_pa   CritPoints,
                                    const Bundles_pa      Bundles,
                                    const float           CPTolerance,
                                    MidPlnGradPath_pa     MidPlnGradPath)
{
    LgIndex_t MinLenGPNum = -1;
    Boolean_t IsOk = TRUE;
    XYZ_s SeedPos;

    double c1_x, c1_y, c1_z;
    double c2_x, c2_y, c2_z;

    double nx = MidPlnGradPath->Normal[0];
    double ny = MidPlnGradPath->Normal[1];
    double nz = MidPlnGradPath->Normal[2];


    REQUIRE(VALID_REF(ZoneVarInfo));
    REQUIRE(CritPointsIsValid(CritPoints));
    REQUIRE(BeginCrtPtNum >= 0 && BeginCrtPtNum < CritPointsGetCount(CritPoints));
    REQUIRE(EndCrtPtNum >= 0 && EndCrtPtNum < CritPointsGetCount(CritPoints));
    REQUIRE(BundlesIsValid(Bundles));
    REQUIRE(CPTolerance > 0.0);
    REQUIRE(MidPlnGradPathIsValid(MidPlnGradPath));
    // Begining and ending critical points must be Cage and Atom
    REQUIRE((CritPointsGetType(CritPoints, BeginCrtPtNum) == -3 &&
             CritPointsGetType(CritPoints, EndCrtPtNum) == 3) ||
            (CritPointsGetType(CritPoints, BeginCrtPtNum) == 3 &&
             CritPointsGetType(CritPoints, EndCrtPtNum) == -3));

    // Determine the local independent coordinates (c1, c2, Normal).
    if (ABS(ny) < 0.95)
    {
        // (y-axis) x (normal)  gives Coord1-axis
        c1_x =  nz;
        c1_y = 0.0;
        c1_z = -nx;
    }
    else  // (normal) nearly parallel to (y-axis)
    {
        // (normal) x (z-axis) gives Coord1-axis
        c1_x =  ny;
        c1_y = -nx;
        c1_z = 0.0;
    }

    // normalise Coord1-axis components and compute Coord2-axis comps
    if (IsOk)
    {
        double c1_tot = sqrt(c1_x * c1_x + c1_y * c1_y + c1_z * c1_z);
        c1_x = c1_x / c1_tot;
        c1_y = c1_y / c1_tot;
        c1_z = c1_z / c1_tot;

        c2_x = ny * c1_z - nz * c1_y;
        c2_y = nz * c1_x - nx * c1_z;
        c2_z = nx * c1_y - ny * c1_x;
    }

    // Seed additional points near the initial seed point.
    if (IsOk)
    {
        XYZ_s SeedPos0 = MidPlnGradPathGetSeedPos(MidPlnGradPath, 0);
        double PathLen1, PathLen2, PathLen3, PathLen4, PathLen5;
        GradPath_pa GradPath = NULL;
        LgIndex_t NumPathPoints;

        // Get primary GradPath length
        PathLen1 = MidPlnGradPathGetPathLen(MidPlnGradPath, 0);

        // Second GradPath: Primary seed position + dc1
        SeedPos.X = SeedPos0.X + 0.1 * c1_x;
        SeedPos.Y = SeedPos0.Y + 0.1 * c1_y;
        SeedPos.Z = SeedPos0.Z + 0.1 * c1_z;
        // SeedPos.Z = SeedPos0.Z - ( nx * ( SeedPos.X - SeedPos0.X ) +
        //                            ny * ( SeedPos.Y - SeedPos0.Y ) ) / nz;

        // Integrate gradpath from seed position.
        GradPath = GradPathAddMidway(ZoneVarInfo, MAXGRADPATHPOINTS, CritPoints,
                                     SeedPos, 1.0, &NumPathPoints);
        if (GradPath == NULL) IsOk = FALSE;

        // Add second GradPath to MidPlnGradPath
        if (IsOk)
        {
            double PathLength;

            // Find the arclength of the gradpath
            PathLength = GradPathGetLength(GradPath);

            PathLen2 = PathLength;

            IsOk = MidPlnGradPathAppendGP(MidPlnGradPath, SeedPos, PathLength, GradPath);
        }

        // Third GradPath: Primary seed position + dc2
        SeedPos.X = SeedPos0.X + 0.1 * c2_x;
        SeedPos.Y = SeedPos0.Y + 0.1 * c2_y;
        SeedPos.Z = SeedPos0.Z + 0.1 * c2_z;
        // SeedPos.Z = SeedPos0.Z - ( nx * ( SeedPos.X - SeedPos0.X ) +
        //                            ny * ( SeedPos.Y - SeedPos0.Y ) ) / nz;

        // Integrate gradpath from seed position.
        GradPath = GradPathAddMidway(ZoneVarInfo, MAXGRADPATHPOINTS, CritPoints,
                                     SeedPos, 1.0, &NumPathPoints);
        if (GradPath == NULL) IsOk = FALSE;

        // Add third GradPath to MidPlnGradPath
        if (IsOk)
        {
            double PathLength;

            // Find the arclength of the gradpath
            PathLength = GradPathGetLength(GradPath);

            PathLen3 = PathLength;

            IsOk = MidPlnGradPathAppendGP(MidPlnGradPath, SeedPos, PathLength, GradPath);
        }

        // fourth GradPath: Primary seed position - dc1
        SeedPos.X = SeedPos0.X - 0.1 * c1_x;
        SeedPos.Y = SeedPos0.Y - 0.1 * c1_y;
        SeedPos.Z = SeedPos0.Z - 0.1 * c1_z;
        // SeedPos.Z = SeedPos0.Z - ( nx * ( SeedPos.X - SeedPos0.X ) +
        //                            ny * ( SeedPos.Y - SeedPos0.Y ) ) / nz;

        // Integrate gradpath from seed position.
        GradPath = GradPathAddMidway(ZoneVarInfo, MAXGRADPATHPOINTS, CritPoints,
                                     SeedPos, 1.0, &NumPathPoints);
        if (GradPath == NULL) IsOk = FALSE;

        // Add fourth GradPath to MidPlnGradPath
        if (IsOk)
        {
            double PathLength;

            // Find the arclength of the gradpath
            PathLength = GradPathGetLength(GradPath);

            PathLen4 = PathLength;

            IsOk = MidPlnGradPathAppendGP(MidPlnGradPath, SeedPos, PathLength, GradPath);
        }

        // Fifth GradPath: Primary seed position - dc2
        SeedPos.X = SeedPos0.X - 0.1 * c2_x;
        SeedPos.Y = SeedPos0.Y - 0.1 * c2_y;
        SeedPos.Z = SeedPos0.Z - 0.1 * c2_z;
        // SeedPos.Z = SeedPos0.Z - ( nx * ( SeedPos.X - SeedPos0.X ) +
        //                           ny * ( SeedPos.Y - SeedPos0.Y ) ) / nz;

        // Integrate gradpath from seed position.
        GradPath = GradPathAddMidway(ZoneVarInfo, MAXGRADPATHPOINTS, CritPoints,
                                     SeedPos, 1.0, &NumPathPoints);
        if (GradPath == NULL) IsOk = FALSE;

        // Add fifth GradPath to MidPlnGradPath
        if (IsOk)
        {
            double PathLength;

            // Find the arclength of the gradpath
            PathLength = GradPathGetLength(GradPath);

            PathLen5 = PathLength;

            IsOk = MidPlnGradPathAppendGP(MidPlnGradPath, SeedPos, PathLength, GradPath);
        }

        // Sixth GradPath: Primary seed position (+/-dx, +/-dy). +/- depends upon which
        //  pathlength is smaller, +dx,0 or -dx,0 for dx and 0,+dy or 0,-dy for dy.

        if (IsOk)
        {
            double dc1 = 0.1;
            double dc2 = 0.1;
            if (PathLen4 < PathLen2) dc1 = -dc1;
            if (PathLen5 < PathLen3) dc2 = -dc2;
            SeedPos.X = SeedPos0.X + dc1 * c1_x + dc2 * c2_x;
            SeedPos.Y = SeedPos0.Y + dc1 * c1_y + dc2 * c2_y;
            SeedPos.Z = SeedPos0.Z + dc1 * c1_z + dc2 * c2_z;
        }
        // SeedPos.X = SeedPos0.X + 0.1;
        // if (PathLen4 < PathLen2) SeedPos.X = SeedPos0.X - 0.1;
        // SeedPos.Y = SeedPos0.Y + 0.1;
        // if (PathLen5 < PathLen3) SeedPos.Y = SeedPos0.Y - 0.1;
        // SeedPos.Z = SeedPos0.Z - ( nx * ( SeedPos.X - SeedPos0.X ) +
        //                            ny * ( SeedPos.Y - SeedPos0.Y ) ) / nz;

        // Integrate gradpath from seed position.
        GradPath = GradPathAddMidway(ZoneVarInfo, MAXGRADPATHPOINTS, CritPoints,
                                     SeedPos, 1.0, &NumPathPoints);
        if (GradPath == NULL) IsOk = FALSE;

        // Add fifth GradPath to MidPlnGradPath
        if (IsOk)
        {
            double PathLength;

            // Find the arclength of the gradpath
            PathLength = GradPathGetLength(GradPath);

            IsOk = MidPlnGradPathAppendGP(MidPlnGradPath, SeedPos, PathLength, GradPath);
        }
    }


    // Fit a quadratic surface to the PathLengths for the seeded points, and
    // use it to predict the seed location for the minimum length GradPath
    // TODO:
    if (IsOk)
    {
        SurfaceFit_pa SurfaceFit = SurfaceFitAlloc();
        XYZ_s     SeedPos0 = MidPlnGradPathGetSeedPos(MidPlnGradPath, 0);
        LgIndex_t NumPaths = MidPlnGradPathGetGPCount(MidPlnGradPath);
        LgIndex_t ii;
        for (ii = 0; IsOk && ii < NumPaths; ii++)
        {
            double Coord1, Coord2, PathLen;
            XYZ_s SeedPos  = MidPlnGradPathGetSeedPos(MidPlnGradPath, ii);

            // determine independent coordinates
            Coord1 = c1_x * (SeedPos.X - SeedPos0.X) +
                     c1_y * (SeedPos.Y - SeedPos0.Y) +
                     c1_z * (SeedPos.Z - SeedPos0.Z);
            Coord2 = c2_x * (SeedPos.X - SeedPos0.X) +
                     c2_y * (SeedPos.Y - SeedPos0.Y) +
                     c2_z * (SeedPos.Z - SeedPos0.Z);
            PathLen = MidPlnGradPathGetPathLen(MidPlnGradPath, ii);

            IsOk = SurfaceFitAppendPointAtEnd(SurfaceFit, Coord1, Coord2, PathLen);
        }

        // Perform a linear-least-square quadratic surface fit of gradpath lengths (surface tangent
        // orthogonal coordinates being the independent coordinates for the fit.
        IsOk = SurfaceFitCompute(SurfaceFit, SurfFitType_Quadratic);
    }

    // Find the level-point (critical point) of the quadratic surface

    return MinLenGPNum;
}







/*
 * MidPlnGradPathMinLengthGP: Find the minimum length gradient path to the
 *  specified EndCrtPtNum. Use the GradPaths in MidPlnGradPath between GPBegin
 *  and GPEnd.
 *
 * param MidPlnGradPath
 *     MidPlnGradPath data structure to populate for an Atom or Cage critical point.
 * param GPBegin, GPEnd
 *     Range of GradientPaths, in MidPlnGradPath, for testing
 * param EndCrtPtNum
 *     Atom or Cage critical point that is the termination of the desired grad paths.
 * param MinLength
 *     Length of shortest GradPath in range, that terminates at correct critical point.
 * param MinLenGPNum
 *     Number (in list) of shortest GradPath in range that terminates at correct
 *     critical point.
 * param NumGPs
 *     Number of GradPaths, in list, that terminate at the correct critical point
 *
 * return
 *     TRUE if successful (no errors), FALSE if there were errors.
 */
Boolean_t MidPlnGradPathMinLengthGP(const MidPlnGradPath_pa  MidPlnGradPath,
                                    const LgIndex_t          GPBegin,
                                    const LgIndex_t          GPEnd,
                                    const LgIndex_t          EndCrtPtNum,
                                    double                  *MinLength,
                                    LgIndex_t               *MinLenGPNum,
                                    LgIndex_t               *NumGPs)
{
    Boolean_t IsOk = TRUE;
    int ii;

    REQUIRE(MidPlnGradPathIsValid(MidPlnGradPath));
    REQUIRE(GPBegin < GPEnd);
    REQUIRE(EndCrtPtNum >= 0);
    REQUIRE(VALID_REF(MinLength));
    REQUIRE(*MinLength > 0.0);
    REQUIRE(VALID_REF(MinLenGPNum));
    REQUIRE(VALID_REF(NumGPs));
    REQUIRE(*NumGPs >= 0);

    if (GPEnd > GPBegin)
    {
        for (ii = GPBegin; IsOk && ii < GPEnd; ii++)
        {
            GradPath_pa GPTry  = MidPlnGradPathGetGP(MidPlnGradPath, ii);
            if (GPTry != NULL)
            {
                if (GPTry->EndCrtPtNum == EndCrtPtNum)
                {
                    double Length = GradPathGetLength(GPTry);

                    (*NumGPs)++;

                    if (Length < *MinLength)
                    {
                        *MinLength   = Length;
                        *MinLenGPNum = ii;
                    }
                }
            }
            else
                IsOk = FALSE;
        }
    }


    ENSURE(VALID_BOOLEAN(IsOk));
    return IsOk;
}






/*
 * MidPlnGradPathGet2CPPath: Find the minimum length gradient path between the
 *  specified BeginCrtPtNum and EndCrtPtNum. Use the GradPaths in MidPlnGradPath,
 *  and refine the seed point grid on the MidPln (using triangle subdivision) if
 *  necessary.
 *
 * param ZoneVarInfo
 *     Zone and variable information.
 * param BeginCrtPtNum
 *     Atom or Cage critical point that is the source of all grad paths in this
 *     MidPlnGradPath data structure
 * param EndCrtPtNum
 *     Atom or Cage critical point that is the termination of the desired grad paths.
 * param CritPoints
 *     Critical Points structure
 * param PathDir
 *     Direction for pathline integration (Forward or Reverse)
 * param CPTolerance
 *     Tolerance - distance from critical point at which the pathline is said
 *     to terminate at the critical point.
 * param MidPlnGradPath
 *     MidPlnGradPath data structure to populate for an Atom or Cage critical point.
 *
 * return
 *     TRUE if successful, FALSE if there were errors.
 */
/*
GradPath_pa MidPlnGradPath2CPPath(const ZoneVarInfo_pa  ZoneVarInfo,
                                  const LgIndex_t       BeginCrtPtNum,
                                  const LgIndex_t       EndCrtPtNum,
                                  const CritPoints_pa   CritPoints,
                                  const StreamDir_e     PathDir,
                                  const float           CPTolerance,
                                  MidPlnGradPath_pa     MidPlnGradPath)
{
  GradPath_pa Result = NULL;
  Boolean_t   IsOk   = TRUE;
  LgIndex_t   MinLenGPNum = -1;
  LgIndex_t   RootMinLenGPNum;
  double      MinLength = LARGEFLOAT;
  LgIndex_t   NumGPs = 0;
  LgIndex_t   SGPCountOld = 0;
  LgIndex_t   SGPCountNew = 12;

  REQUIRE(VALID_REF(ZoneVarInfo));
  REQUIRE(VALID_REF(CritPoints));
  REQUIRE(BeginCrtPtNum >= 0 && BeginCrtPtNum < CritPointsGetCount(CritPoints) );
  REQUIRE(EndCrtPtNum >= 0 && EndCrtPtNum < CritPointsGetCount(CritPoints) );
  REQUIRE(CritPointsIsValid(CritPoints));
  REQUIRE(VALID_ENUM(PathDir, StreamDir_e));
  REQUIRE(CPTolerance > 0.0);

  //
  // Find the min-length GradPath from the current set that connects BeginCrtPt
  // and EndCrtPtNum
  //
  if (IsOk)
    {
      IsOk = MidPlnGradPathMinLengthGP(MidPlnGradPath, SGPCountOld, SGPCountNew, EndCrtPtNum,
                                       &MinLength, &MinLenGPNum, &NumGPs);
      if (IsOk && MinLenGPNum >= 0)
        Result = MidPlnGradPathGetGP(MidPlnGradPath, MinLenGPNum);

      RootMinLenGPNum = MinLenGPNum;
    }


  //
  // If fewer than 5 GP's connect the two points, subdivide root triangles containing
  // the min-length connecting line.
  //
  if (IsOk)
    {
      LgIndex_t ii;

      for (ii=0; IsOk && ii<20; ii++)
        {
          TriGradPath_pa Triangle = MidPlnGradPathGetRtTri(MidPlnGradPath, ii);
          if (Triangle->GradPathNum[0] == MinLenGPNum ||
              Triangle->GradPathNum[1] == MinLenGPNum ||
              Triangle->GradPathNum[2] == MinLenGPNum)
            {
              IsOk = TriGradPathSubdivide(ZoneVarInfo, CritPoints, Triangle,
                                          CPTolerance, MidPlnGradPath);
            }
        }

      SGPCountOld = SGPCountNew;
      SGPCountNew = MidPlnGradPathGetGPCount(MidPlnGradPath);

      if (SGPCountNew > SGPCountOld)
        {
          IsOk = MidPlnGradPathMinLengthGP(MidPlnGradPath, SGPCountOld, SGPCountNew, EndCrtPtNum,
                                           &MinLength, &MinLenGPNum, &NumGPs);
          if (IsOk && MinLenGPNum >= SGPCountOld)
            Result = MidPlnGradPathGetGP(MidPlnGradPath, MinLenGPNum);
        }
    }

  //
  // If fewer than 5 GP's connect the two points, subdivide the subtriangles containing
  // the min-length connecting line again.
  //
  if (IsOk)
    {
      LgIndex_t ii;

      for (ii=0; IsOk && ii<20; ii++)
        {
          TriGradPath_pa Triangle = MidPlnGradPathGetRtTri(MidPlnGradPath, ii);
          if (Triangle->GradPathNum[0] == RootMinLenGPNum ||
              Triangle->GradPathNum[1] == RootMinLenGPNum ||
              Triangle->GradPathNum[2] == RootMinLenGPNum)
            {
              LgIndex_t jj;
              for (jj=0; IsOk && jj<4; jj++)
                {
                  TriGradPath_pa SubTriangle = TriGradPathGetSubTri(Triangle, jj);
                  if(MinLenGPNum == RootMinLenGPNum ||
                     SubTriangle->GradPathNum[0] == MinLenGPNum ||
                     SubTriangle->GradPathNum[1] == MinLenGPNum ||
                     SubTriangle->GradPathNum[2] == MinLenGPNum)
                    {
                      IsOk = TriGradPathSubdivide(ZoneVarInfo, CritPoints, SubTriangle,
                                                  CPTolerance, MidPlnGradPath);

                      // If all six gradpaths in the subdivided triangle go to the Cage,
                      // assume a quadratic distribution of pathlength and compute the
                      // seed location for the minimum.
// TODO: Haven't tested yet
                      if (IsOk)
                        {
                          int kk;
                          LgIndex_t   GradPathNums[6];
                          GradPath_pa GradPaths[6];
                          TriGradPath_pa SubSubTri0, SubSubTri1, SubSubTri2;
                          SubSubTri0 = TriGradPathGetSubTri(SubTriangle, 0);
                          SubSubTri1 = TriGradPathGetSubTri(SubTriangle, 1);
                          SubSubTri2 = TriGradPathGetSubTri(SubTriangle, 2);
                          GradPathNums[0] = SubTriangle->GradPathNum[0];
                          GradPathNums[1] = SubTriangle->GradPathNum[1];
                          GradPathNums[2] = SubTriangle->GradPathNum[2];
                          GradPathNums[3] = SubSubTri0->GradPathNum[1];
                          GradPathNums[4] = SubSubTri1->GradPathNum[1];
                          GradPathNums[5] = SubSubTri2->GradPathNum[1];
                          for (kk=0; IsOk && kk<6; kk++)
                            {
                              GradPaths[kk] = MidPlnGradPathGetGP(MidPlnGradPath, GradPathNums[kk]);
                            }

                          if (GradPaths[0]->EndCrtPtNum == EndCrtPtNum &&
                              GradPaths[1]->EndCrtPtNum == EndCrtPtNum &&
                              GradPaths[2]->EndCrtPtNum == EndCrtPtNum &&
                              GradPaths[3]->EndCrtPtNum == EndCrtPtNum &&
                              GradPaths[4]->EndCrtPtNum == EndCrtPtNum &&
                              GradPaths[5]->EndCrtPtNum == EndCrtPtNum)
                            {
                              double Len1 = GradPathGetLength(GradPaths[0]);
                              double Len2 = GradPathGetLength(GradPaths[1]);
                              double Len3 = GradPathGetLength(GradPaths[2]);
                              double Len4 = GradPathGetLength(GradPaths[3]);
                              double Len5 = GradPathGetLength(GradPaths[4]);
                              double Len6 = GradPathGetLength(GradPaths[5]);
                              // Based on a quadratic finite-element, solve for the first two area coords
                              double A = 4.0 * ( Len1 + Len3 ) - 8.0 * Len6;
                              double B = 4.0 * ( Len3 + Len4 - Len5 - Len6 );
                              double C = Len1 + 3.0 * Len3 - 4.0 * Len4;
                              double D = B;
                              double E = 4.0 * ( Len2 + Len3 ) - 8.0 * Len5;
                              double F = Len2 + 3.0 * Len3 - 4.0 * Len5;
                              double Denom = B * D - A * E;
                              double AreaCd1, AreaCd2, AreaCd3;
                              double NtCd[6];  // Natural coordinates for min point
                              double XSeed = 0.0;
                              double YSeed = 0.0;
                              double ZSeed = 0.0;
                              XYZ_s XYZNode;
                              if (Denom != 0.0 && A != 0.0)
                                {
                                  AreaCd2 = ( C * D - A * F ) / Denom;
                                  AreaCd1 = ( C - B * AreaCd2 ) / A;
                                  AreaCd3 = 1.0 - AreaCd1 - AreaCd2;
                                  NtCd[0] = AreaCd1 * ( 2.0 * AreaCd1 - 1.0 );
                                  NtCd[1] = AreaCd2 * ( 2.0 * AreaCd2 - 1.0 );
                                  NtCd[2] = AreaCd3 * ( 2.0 * AreaCd3 - 1.0 );
                                  NtCd[3] = 4.0 * AreaCd1 * AreaCd2;
                                  NtCd[4] = 4.0 * AreaCd2 * AreaCd3;
                                  NtCd[5] = 4.0 * AreaCd3 * AreaCd1;
                                }
                              for (kk=0; IsOk && kk<6; kk++)
                                {
                                  XYZNode = MidPlnGradPathGetSeedDir(MidPlnGradPath, GradPathNums[kk]);
                                  XSeed += NtCd[kk] * XYZNode.X;
                                  YSeed += NtCd[kk] * XYZNode.Y;
                                  ZSeed += NtCd[kk] * XYZNode.Z;
                                }
                            }
                        }
                    }
                }
            }
        }

      SGPCountOld = SGPCountNew;
      SGPCountNew = MidPlnGradPathGetGPCount(MidPlnGradPath);

      if (SGPCountNew > SGPCountOld)
        {
          IsOk = MidPlnGradPathMinLengthGP(MidPlnGradPath, SGPCountOld, SGPCountNew, EndCrtPtNum,
                                           &MinLength, &MinLenGPNum, &NumGPs);
          if (IsOk && MinLenGPNum >= SGPCountOld)
            Result = MidPlnGradPathGetGP(MidPlnGradPath, MinLenGPNum);
        }
    }

  //
  // If all points within a triangle and it's subtriangle end at the cage,
  // find the X,Y,Z seed that will give the minimum PathLength assuming that
  // PathLength varies quadratically across the triangle.
  // NOT CURRENTLY USED
  //
  if (0 && IsOk)
    {
      Boolean_t IsDone = FALSE;
      TriGradPath_pa TriPathStack[MAX_TRI_SUBDIVISIONS];
      EntIndex_t     SubTriNumStack[MAX_TRI_SUBDIVISIONS];

      // Start at root and search down branches
      TriGradPath_pa CurrentTri = NULL;
      EntIndex_t     CurSubTriNum = 0;
      EntIndex_t     CurTriNum = 0;
      EntIndex_t     Level = 0;
      while(!IsDone)
        {
          // Special treatment at Root level: Are we done? If not, set CurrentTri
          if (Level == 0)
            {
              if (CurTriNum > 19 || MidPlnGradPath->RootTriangles[CurTriNum] == NULL )
                IsDone = TRUE;
              else
                {
                  CurrentTri = MidPlnGradPath->RootTriangles[CurSubTriNum];
                  CurSubTriNum = 0;
                }
            }

          // Move down a level
          if ( ( Level > 0 && CurSubTriNum < 4 )  && CurrentTri->SubTriangles[CurSubTriNum] != NULL)
            {
              TriPathStack[Level] = CurrentTri;
              SubTriNumStack[Level] = CurTriNum;
              CurrentTri = CurrentTri->SubTriangles[CurSubTriNum];
              CurTriNum = CurSubTriNum;
              CurSubTriNum = 0;
              Level++;
            }
          // CurSubTriNum is a NULL TriGradPath
          else
            {
              // If more SubTri's at level, move laterally to next SubTri at level
              if ( (Level > 0 && CurSubTriNum < 3) || (Level == 0 && CurSubTriNum < 19) )
                {
                  CurSubTriNum++;
                }
              // Move up and laterally
              else
                {
                  // Pop Tri and CurrentTri from stack
                  Level--;
                  CurrentTri = TriPathStack[Level];
                  CurSubTriNum = SubTriNumStack[Level] + 1;
                }
            }
        }

      // Are all six nodes of triangle and subtri's ending at the desired cage?
    }

  if (!IsOk)
    Result = NULL;

  ENSURE(Result == NULL || GradPathIsValid(Result));

  return Result;
}
*/