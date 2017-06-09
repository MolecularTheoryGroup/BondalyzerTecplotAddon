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
// #include "NORMALS.h"
#include "GRADPATH.h"
#include "SPHEREGRADPATH.h"
#include "ENGINE.h"



static XYZ_s Sphere20FacetNodes[] =
{
    {0.000000000E+000, 4.253254042E-001, 2.628655561E-001},
    {0.000000000E+000, -4.253254042E-001, 2.628655561E-001},
    {0.000000000E+000, -4.253254042E-001, -2.628655561E-001},
    {0.000000000E+000, 4.253254042E-001, -2.628655561E-001},
    {2.628655561E-001, 0.000000000E+000, 4.253254042E-001},
    {2.628655561E-001, 0.000000000E+000, -4.253254042E-001},
    { -2.628655561E-001, 0.000000000E+000, -4.253254042E-001},
    { -2.628655561E-001, 0.000000000E+000, 4.253254042E-001},
    {4.253254042E-001, 2.628655561E-001, 0.000000000E+000},
    { -4.253254042E-001, 2.628655561E-001, 0.000000000E+000},
    { -4.253254042E-001, -2.628655561E-001, 0.000000000E+000},
    {4.253254042E-001, -2.628655561E-001, 0.000000000E+000}
};
static LgIndex_t Sphere20FacetTriangulation[] =
{
    0, 4, 8, 11, 5, 2, 6, 10, 9, 7, 0, 4, 4, 4, 4, 7, 1, 10, 2, 2, 5, 5, 5, 6, 3, 9, 0, 0, 5, 5, 5, 3, 8, 0, 0, 2, 2, 11, 1, 4
};




/*
 * Check to see if entry Offset, Offset+1, and Offset+2 of Sphere20FacetTriangulation
 * are unique node numbers that define a triangle.
 *
 * param SphereFacetTriangulation
 *     Strip style specification of node-map.
 *
 * param Offset
 *     Offset into array of first node in triangle.
 *
 * return
 *     TRUE if is a triangle, otherwise FALSE.
 */
Boolean_t IsTriangle(LgIndex_t *SphereFacetTriangulation,
                     int        Offset)
{
    Boolean_t IsTriangle = TRUE;

    REQUIRE(VALID_REF(SphereFacetTriangulation));
    REQUIRE(Offset < 40);

    if (SphereFacetTriangulation[Offset  ] == SphereFacetTriangulation[Offset+1] ||
        SphereFacetTriangulation[Offset+1] == SphereFacetTriangulation[Offset+2] ||
        SphereFacetTriangulation[Offset  ] == SphereFacetTriangulation[Offset+2]) IsTriangle = FALSE;

    ENSURE(VALID_BOOLEAN(IsTriangle));

    return IsTriangle;
}




/**
 * Determine if the TriGradPath handle is sane.
 *
 * param TriGradPath
 *     TriGradPath structure in question.
 *
 * Note: Don't evaluate upper limit of GradPathNum in this test.
 *
 * return
 *     TRUE if the TriGradPath structure is valid, otherwise FALSE.
 */
Boolean_t TriGradPathIsValid(TriGradPath_pa  TriGradPath)
{
    Boolean_t IsValid = FALSE;

    IsValid = (VALID_REF(TriGradPath) &&
               TriGradPath->GradPathNum[0] >= -1 &&
               TriGradPath->GradPathNum[1] >= -1 &&
               TriGradPath->GradPathNum[2] >= -1);

    ENSURE(VALID_BOOLEAN(IsValid));
    return IsValid;
}







/**
 * Gets the Offset, into the SphereGradPath->GradPaths list, of the GradPath at
 * the specified node of the triangle.
 *
 *
 * param TriGradPath
 *     TriGradPath to get the GradPath from.
 *
 * param CornerNumber
 *     CornerNumber of the desired triangle node.
 *
 * return
 *     Offset of GradPath in list, otherewise -1.
 */
LgIndex_t TriGradPathGetGPIndex(TriGradPath_pa TriGradPath,
                                LgIndex_t      CornerNumber)
{
    LgIndex_t Result = -1;

    REQUIRE(VALID_REF(TriGradPath));
    REQUIRE(TriGradPathIsValid(TriGradPath));
    REQUIRE(CornerNumber >= 0 && CornerNumber < 3);

    Result = TriGradPath->GradPathNum[CornerNumber];

    ENSURE(Result >= 0);
    return Result;
}







/**
 * Sets the Offset, into the SphereGradPath->GradPaths list, of the GradPath at
 * the specified node of the triangle.
 *
 *
 * param TriGradPath
 *     TriGradPath to get the GradPath from.
 *
 * param CornerNumber
 *     CornerNumber of the desired triangle node.
 *
 * param GradPathIndex
 *     Offset, into the SphereGradPath->GradPaths list, of the handle to the
 *     GradPath at the CornerNumber triangle node.
 *
 * return
 *     TRUE if successful, otherewise FALSE.
 */
Boolean_t TriGradPathSetGPIndex(TriGradPath_pa TriGradPath,
                                LgIndex_t      CornerNumber,
                                LgIndex_t      GradPathIndex)
{
    Boolean_t IsOk = TRUE;

    REQUIRE(VALID_REF(TriGradPath));
    REQUIRE(TriGradPathIsValid(TriGradPath));
    REQUIRE(CornerNumber >= 0 && CornerNumber < 3);
    REQUIRE(GradPathIndex >= 0);

    TriGradPath->GradPathNum[CornerNumber] = GradPathIndex;

    return IsOk;
}







/**
 * Gets the TriGradPath handle for the specified subtriangle.
 *
 *
 * param TriGradPath
 *     Handle of the TriGradPath to get the GradPath from.
 *
 * param SubTriIndex
 *     Index of the desired sub-triangle. Range: 0-3.
 *
 *           2
 *          /\
 *         /  \
 *        / s2 \
 *     c /______\ b
 *      /\      /\
 *     /  \ s3 /  \
 *    / s0 \  / s1 \
 *   /______\/______\
 *  0       a        1
 *
 * return
 *     Handle of the TriGradPath for the sub-triangle, or NULL.
 */
TriGradPath_pa    TriGradPathGetSubTri(TriGradPath_pa TriGradPath,
                                       LgIndex_t      SubTriIndex)
{
    TriGradPath_pa Result = NULL;

    REQUIRE(VALID_REF(TriGradPath));
    REQUIRE(TriGradPathIsValid(TriGradPath));
    REQUIRE(SubTriIndex >= 0 && SubTriIndex < 4);

    Result = TriGradPath->SubTriangles[SubTriIndex];

    ENSURE(Result == NULL || TriGradPathIsValid(Result));
    return Result;
}








/**
 * Deallocates the TriGradPath handle and set the handle to NULL.
 *
 * note
 *     Item destruction is the responsibility of the caller.
 *
 * param
 *     Reference to a TriGradPath handle.
 */
void TriGradPathDealloc(TriGradPath_pa *TriGradPath)
{

    REQUIRE(VALID_REF(TriGradPath));
    REQUIRE(TriGradPathIsValid(*TriGradPath) || *TriGradPath == NULL);

    if (*TriGradPath != NULL)
    {
        LgIndex_t ii;

        /* Dealloc any sub-triangles.  */
        for (ii = 0; ii < 4; ii++)
        {
            TriGradPath_pa SubTriGradPath = TriGradPathGetSubTri((*TriGradPath), ii);
            if (SubTriGradPath != NULL) TriGradPathDealloc(&SubTriGradPath);
        }

        /* release the list structure itself */
        FREE_ITEM(*TriGradPath, "TriGradPath structure");
        *TriGradPath = NULL;
    }

    ENSURE(*TriGradPath == NULL);
}







/**
 * Empties the TriGradPath of all SubTriGradPaths and resets the
 * other variables.
 *
 *
 * param TriGradPath
 *     TriGradPath to clear.
 */
void TriGradPathClear(TriGradPath_pa TriGradPath)
{
    LgIndex_t ii;

    REQUIRE(TriGradPathIsValid(TriGradPath));

    for (ii = 0; ii < 4; ii++)
    {
        TriGradPath_pa SubTriGradPath = TriGradPathGetSubTri(TriGradPath, ii);
        if (SubTriGradPath != NULL) TriGradPathDealloc(&SubTriGradPath);
        CHECK(TriGradPathGetSubTri(TriGradPath, ii) == NULL);
    }
    TriGradPath->DepthInTree    = 0;
    TriGradPath->GradPathNum[0] = -1;
    TriGradPath->GradPathNum[1] = -1;
    TriGradPath->GradPathNum[2] = -1;
    TriGradPath->NeighborTri[0] = -1;
    TriGradPath->NeighborTri[1] = -1;
    TriGradPath->NeighborTri[2] = -1;

    ENSURE(TriGradPathIsValid(TriGradPath));
}




/**
 * Allocates a TriGradPath handle.
 *
 *
 * return
 *     TriGradPath handle if sufficient memory was available,
 *     otherewise a handle is NULL.
 */
TriGradPath_pa TriGradPathAlloc()
{
    TriGradPath_pa Result = NULL;

    Result = ALLOC_ITEM(TriGradPath_s, "TriGradPath structure");
    if (Result != NULL)
    {
        int ii;
        Result->DepthInTree    = 0;
        for (ii = 0; ii < 3; ii++)
        {
            Result->GradPathNum[ii] = -1;
            Result->NeighborTri[ii] = -1;
        }
        for (ii = 0; ii < 4; ii++)
            Result->SubTriangles[ii] = NULL;

    }

    ENSURE(TriGradPathIsValid(Result) || Result == NULL);
    return Result;
}






/**
 * Sets the TriGradPath handle for the specified subtriangle.
 *
 *
 * param TriGradPath
 *     Handle of the TriGradPath containing the desired sub-triangle.
 *
 * param SubTriIndex
 *     Index of the desired sub-triangle. Range: 0-3.
 *
 * param SubTriGradPath
 *     Handle of the sub-triangle TriGradPath, or NULL if no sub-triangle.
 *
 * return
 *     TRUE if successful, otherewise FALSE.
 */
Boolean_t         TriGradPathSetSubTri(TriGradPath_pa TriGradPath,
                                       LgIndex_t      SubTriIndex,
                                       TriGradPath_pa SubTriGradPath)
{
    Boolean_t IsOk = TRUE;

    REQUIRE(VALID_REF(TriGradPath));
    REQUIRE(TriGradPathIsValid(TriGradPath));
    REQUIRE(SubTriIndex >= 0 && SubTriIndex < 3);
    REQUIRE(SubTriGradPath == NULL || TriGradPathIsValid(SubTriGradPath));

    TriGradPath->SubTriangles[SubTriIndex] = SubTriGradPath;

    return IsOk;
}






/**
 * Determine if the SphereGradPath handle is sane.
 *
 * param SphereGradPath
 *     SphereGradPath structure in question.
 *
 * return
 *     TRUE if the SphereGradPath structure is valid, otherwise FALSE.
 */
Boolean_t SphereGradPathIsValid(SphereGradPath_pa SphereGradPath)
{
    Boolean_t IsValid = FALSE;
    int ii;

    IsValid = (VALID_REF(SphereGradPath) &&
               VALID_REF(SphereGradPath->SeedVectorX) && ArrListIsValid(SphereGradPath->SeedVectorX) &&
               VALID_REF(SphereGradPath->SeedVectorY) && ArrListIsValid(SphereGradPath->SeedVectorY) &&
               VALID_REF(SphereGradPath->SeedVectorZ) && ArrListIsValid(SphereGradPath->SeedVectorZ) &&
               VALID_REF(SphereGradPath->GradPaths)   && ArrListIsValid(SphereGradPath->GradPaths) &&
               VALID_REF(SphereGradPath->RootTriangles));

    /* Require the same count for each coordinate array. */
    if (IsValid) IsValid = (ArrListGetCount(SphereGradPath->SeedVectorX) ==
                                ArrListGetCount(SphereGradPath->GradPaths) &&
                                ArrListGetCount(SphereGradPath->SeedVectorY) ==
                                ArrListGetCount(SphereGradPath->GradPaths) &&
                                ArrListGetCount(SphereGradPath->SeedVectorZ) ==
                                ArrListGetCount(SphereGradPath->GradPaths));

    /* Require that RootTriangles be valid or NULL */
    for (ii = 0; IsValid && ii < 20; ii++)
    {
        TriGradPath_pa RootTriangle = SphereGradPath->RootTriangles[ii];
        if (RootTriangle != NULL)
            IsValid = TriGradPathIsValid(RootTriangle);
    }

    /* Require that Radius be positive */
    if (IsValid) IsValid = (SphereGradPath->Radius > 0.0);

    ENSURE(VALID_BOOLEAN(IsValid));
    return IsValid;
}





/**
 * Gets the Index'th GradPath handle.
 *
 *
 * param SphereGradPath
 *     SphereGradPath to get the GradPath from.
 *
 * param Offset
 *     Offset into the GradPath pointer ArrList.
 *
 * return
 *     GradPath handle available, otherewise a handle to NULL.
 */
GradPath_pa SphereGradPathGetGP(SphereGradPath_pa SphereGradPath,
                                LgIndex_t         Offset)
{
    GradPath_pa   Result = NULL;
    ArrListItem_u Item;

    REQUIRE(VALID_REF(SphereGradPath));
    REQUIRE(SphereGradPathIsValid(SphereGradPath));
    REQUIRE(Offset >= 0 && Offset < SphereGradPathGetGPCount(SphereGradPath));

    Item = ArrListGetItem(SphereGradPath->GradPaths, Offset);
    Result = (GradPath_pa)Item.VoidPtr;

    ENSURE(VALID_REF(Result));
    return Result;
}





/**
 * Gets the relative start direction vector for the specified GridPath.
 *
 *
 * param SphereGradPath
 *     SphereGradPath from which to get the Angle.
 *
 * param Offset
 *     Offset into the Angle ArrList.
 *
 * return
 *     double Angle.
 */
XYZ_s SphereGradPathGetSeedDir(SphereGradPath_pa SphereGradPath,
                               LgIndex_t         Offset)
{
    XYZ_s  Result;
    ArrListItem_u Item;

    REQUIRE(VALID_REF(SphereGradPath));
    REQUIRE(SphereGradPathIsValid(SphereGradPath));
    REQUIRE(Offset >= 0 && Offset < SphereGradPathGetGPCount(SphereGradPath));

    Item = ArrListGetItem(SphereGradPath->SeedVectorX, Offset);
    Result.X = Item.Double;
    Item = ArrListGetItem(SphereGradPath->SeedVectorY, Offset);
    Result.Y = Item.Double;
    Item = ArrListGetItem(SphereGradPath->SeedVectorZ, Offset);
    Result.Z = Item.Double;

    return Result;
}






/**
 * Sets the relative start direction vector for the specified GridPath.
 *
 *
 * param SphereGradPath
 *     SphereGradPath from which to get the Angle.
 *
 * param Offset
 *     Offset into the Angle ArrList.
 *
 * param SeedDir
 *     Normalized vector giving direction from critical point to Node.
 *
 * return
 *     TRUE if it works, FALSE otherwise.
 */
Boolean_t SphereGradPathSetSeedDir(SphereGradPath_pa SphereGradPath,
                                   LgIndex_t         Offset,
                                   XYZ_s             SeedDir)
{
    Boolean_t  IsOk = TRUE;
    ArrListItem_u Item;

    REQUIRE(VALID_REF(SphereGradPath));
    REQUIRE(SphereGradPathIsValid(SphereGradPath));
    REQUIRE(Offset >= 0 && Offset < SphereGradPathGetGPCount(SphereGradPath));

    Item.Double = SeedDir.X;
    IsOk = ArrListSetItem(SphereGradPath->SeedVectorX, Offset, Item);

    if (IsOk)
    {
        Item.Double = SeedDir.Y;
        IsOk = ArrListSetItem(SphereGradPath->SeedVectorY, Offset, Item);
    }

    if (IsOk)
    {
        Item.Double = SeedDir.Z;
        IsOk = ArrListSetItem(SphereGradPath->SeedVectorZ, Offset, Item);
    }

    ENSURE(VALID_BOOLEAN(IsOk));
    return IsOk;
}






/**
 * Deallocates the SphereGradPath handle and set the handle to NULL.
 *
 * note
 *     Item destruction is the responsibility of the caller.
 *
 * param
 *     Reference to a SphereGradPath handle.
 */
void SphereGradPathDealloc(SphereGradPath_pa *SphereGradPath)
{

    REQUIRE(VALID_REF(SphereGradPath));
    // REQUIRE(SphereGradPathIsValid(*SphereGradPath) || *SphereGradPath == NULL);

    if (*SphereGradPath != NULL)
    {
        LgIndex_t ii;
        LgIndex_t Count = ArrListGetCount((*SphereGradPath)->GradPaths);

        /* Dealloc the GradPaths */
        for (ii = 0; ii < Count; ii++)
        {
            GradPath_pa GradPath = SphereGradPathGetGP((*SphereGradPath), ii);
            GradPathDealloc(&GradPath);
        }

        /* Dealloc the root triangles */
        for (ii = 0; ii < 20; ii++)
        {
            TriGradPath_pa TriGradPath = (*SphereGradPath)->RootTriangles[ii];
            if (TriGradPath != NULL) TriGradPathDealloc(&TriGradPath);
        }

        /* release the ArrList's */
        if ((*SphereGradPath)->GradPaths != NULL) ArrListDealloc(&((*SphereGradPath)->GradPaths));
        if ((*SphereGradPath)->SeedVectorX != NULL) ArrListDealloc(&((*SphereGradPath)->SeedVectorX));
        if ((*SphereGradPath)->SeedVectorY != NULL) ArrListDealloc(&((*SphereGradPath)->SeedVectorY));
        if ((*SphereGradPath)->SeedVectorZ != NULL) ArrListDealloc(&((*SphereGradPath)->SeedVectorZ));

        /* release the list structure itself */
        FREE_ITEM(*SphereGradPath, "SphereGradPath structure");
        *SphereGradPath = NULL;
    }

    ENSURE(*SphereGradPath == NULL);
}




/**
 * Gets the number of GradPaths/RelNormSeedVector's currently in the SphereGradPath
 * (maintained by the SphereGradPath array lists).
 *
 * param
 *     SphereGradPath structure in question.
 *
 * return
 *     Number of GradPaths/RelNormSeedVector's in the SphereGradPath.
 */
LgIndex_t SphereGradPathGetGPCount(SphereGradPath_pa SphereGradPath)
{
    LgIndex_t Result = 0;

    REQUIRE(SphereGradPathIsValid(SphereGradPath));

    Result = ArrListGetCount(SphereGradPath->GradPaths);

    ENSURE(Result >= 0);
    return Result;
}










/**
 * Empties the SphereGradPath of all SeedVectors, GradPaths, and RootTriangles
 * and resets the other variables.
 *
 *
 * param SphereGradPath
 *     SphereGradPath to clear.
 */
void SphereGradPathClear(SphereGradPath_pa SphereGradPath)
{
    LgIndex_t ii;
    LgIndex_t Count = ArrListGetCount(SphereGradPath->GradPaths);

    REQUIRE(SphereGradPathIsValid(SphereGradPath));

    for (ii = 0; ii < 20; ii++)
    {
        TriGradPath_pa TriGradPath = SphereGradPathGetRtTri(SphereGradPath, ii);
        TriGradPathClear(TriGradPath);
    }

    for (ii = 0; ii < Count; ii++)
    {
        GradPath_pa GradPath = SphereGradPathGetGP(SphereGradPath, ii);
        GradPathDealloc(&GradPath);
    }

    ArrListClear(SphereGradPath->GradPaths);
    ArrListClear(SphereGradPath->SeedVectorX);
    ArrListClear(SphereGradPath->SeedVectorY);
    ArrListClear(SphereGradPath->SeedVectorZ);

    SphereGradPath->Center[0]   = 0.0;
    SphereGradPath->Center[1]   = 0.0;
    SphereGradPath->Center[2]   = 0.0;
    SphereGradPath->Radius      = 1.0;

    ENSURE(SphereGradPathIsValid(SphereGradPath) &&
           SphereGradPathGetGPCount(SphereGradPath) == 0);
}






/**
 * Allocates a SphereGradPath handle with a suitable default
 * capacity for the contained ArrList's.
 *
 *
 * return
 *     SphereGradPath handle if sufficient memory was available,
 *     otherewise a handle to NULL.
 */
SphereGradPath_pa SphereGradPathAlloc()
{
    SphereGradPath_pa Result = NULL;

    Result = ALLOC_ITEM(SphereGradPath_s, "SphereGradPath structure");
    if (Result != NULL)
    {
        Result->Center[0]   = 0.0;
        Result->Center[1]   = 0.0;
        Result->Center[2]   = 0.0;
        Result->Radius      = 1.0;
        Result->SeedVectorX = ArrListAlloc(60, ArrListType_Double);
        Result->SeedVectorY = ArrListAlloc(60, ArrListType_Double);
        Result->SeedVectorZ = ArrListAlloc(60, ArrListType_Double);
        Result->GradPaths   = ArrListAlloc(60, ArrListType_VoidPtr);

        /* TODO: Set up SeedVectors and RootTriangles. Here or elsewhere? */

        /* If it failed to allocate any of the array lists, clean-up and exit. */
        if (Result->SeedVectorX == NULL || Result->SeedVectorY == NULL ||
            Result->SeedVectorZ == NULL || Result->GradPaths == NULL)
        {
            if (Result->SeedVectorX != NULL) ArrListDealloc(&(Result->SeedVectorX));
            if (Result->SeedVectorY != NULL) ArrListDealloc(&(Result->SeedVectorY));
            if (Result->SeedVectorZ != NULL) ArrListDealloc(&(Result->SeedVectorZ));
            if (Result->GradPaths != NULL) ArrListDealloc(&(Result->GradPaths));
            FREE_ITEM(Result, "SphereGradPath structure");
            Result = NULL;
        }

        /* Allocate TriGradPath structures for the root triangles of the icosohedron. */
        if (Result != NULL)
        {
            int ii;
            for (ii = 0; Result != NULL && ii < 20; ii++)
            {
                Result->RootTriangles[ii] = TriGradPathAlloc();
                if (!TriGradPathIsValid(Result->RootTriangles[ii]))
                {
                    SphereGradPathDealloc(&Result);
                    Result = NULL;
                }
            }
        }
    }

    ENSURE(SphereGradPathIsValid(Result) || Result == NULL);
    return Result;
}






/**
 * Places GradPath pointer and SeedVector at the specified offset. If the
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
 * param SphereGradPath
 *     SphereGradPath target in which to set the Angle/GradPath.
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

Boolean_t   SphereGradPathSetGP(SphereGradPath_pa SphereGradPath,
                                LgIndex_t   Offset,
                                XYZ_s       SeedVector,
                                GradPath_pa GradPath)
{
    Boolean_t IsOk = TRUE;
    ArrListItem_u Item;


    REQUIRE(SphereGradPathIsValid(SphereGradPath));
    REQUIRE(Offset >= 0);

    Item.Double = SeedVector.X;
    IsOk = ArrListSetItem(SphereGradPath->SeedVectorX, Offset, Item);

    if (IsOk)
    {
        Item.Double = SeedVector.Y;
        IsOk = ArrListSetItem(SphereGradPath->SeedVectorY, Offset, Item);
    }

    if (IsOk)
    {
        Item.Double = SeedVector.Z;
        IsOk = ArrListSetItem(SphereGradPath->SeedVectorZ, Offset, Item);
    }

    if (IsOk)
    {
        Item.VoidPtr = (void *)GradPath;
        IsOk = ArrListSetItem(SphereGradPath->GradPaths, Offset, Item);
    }

    ENSURE(SphereGradPathIsValid(SphereGradPath));

    ENSURE(VALID_BOOLEAN(IsOk));
    return IsOk;
}






/**
 * Inserts SeedVector/GradPath at the specified offset. The arrays
 * will be expanded to accomodate the additional value.
 *
 * note
 *     If SeedVector/GradPath already exists at the specified location they are
 *     replaced. Therefore, item destruction is the responsibility of
 *     the caller.
 *
 * param SphereGradPath
 *     SphereGradPath target in which to set the coordinates.
 * param Offset
 *     Offset at which to insert the Angle/GradPath.
 * param Angle
 *     Angle to insert at the specified offset.
 * param GradPath
 *     GradPath handle to insert at the specified offset.
 *
 * return
 *     TRUE if sufficient memory permitted the operation, otherwise FALSE.
 */

Boolean_t SphereGradPathInsertGP(SphereGradPath_pa SphereGradPath,
                                 LgIndex_t         Offset,
                                 XYZ_s             SeedVector,
                                 GradPath_pa       GradPath)
{
    Boolean_t IsOk = TRUE;
    ArrListItem_u Item;


    REQUIRE(SphereGradPathIsValid(SphereGradPath));
    REQUIRE(0 <= Offset && Offset <= SphereGradPathGetGPCount(SphereGradPath));

    Item.Double = SeedVector.X;
    IsOk = ArrListInsertItem(SphereGradPath->SeedVectorX, Offset, Item);

    if (IsOk)
    {
        Item.Double = SeedVector.Y;
        IsOk = ArrListInsertItem(SphereGradPath->SeedVectorY, Offset, Item);
    }

    if (IsOk)
    {
        Item.Double = SeedVector.Z;
        IsOk = ArrListInsertItem(SphereGradPath->SeedVectorZ, Offset, Item);
    }

    if (IsOk)
    {
        Item.VoidPtr = (void *)GradPath;
        IsOk = ArrListInsertItem(SphereGradPath->GradPaths, Offset, Item);
    }

    ENSURE(SphereGradPathIsValid(SphereGradPath));

    ENSURE(VALID_BOOLEAN(IsOk));
    return IsOk;
}






/**
 * Deletes SeedVector/GradPath at the specified offset.
 *
 *
 * param SphereGradPath
 *     SphereGradPath target in which to set the coordinates.
 * param Offset
 *     Offset at which to insert the SeedVector/GradPath.
 *
 * return
 *     TRUE if successful operation, otherwise FALSE.
 */

Boolean_t SphereGradPathRemoveGP(SphereGradPath_pa SphereGradPath,
                                 LgIndex_t         Offset)
{
    Boolean_t IsOk = TRUE;
    ArrListItem_u Item;


    REQUIRE(SphereGradPathIsValid(SphereGradPath));
    REQUIRE(0 <= Offset && Offset <= SphereGradPathGetGPCount(SphereGradPath));

    Item = ArrListRemoveItem(SphereGradPath->SeedVectorX, Offset);

    if (IsOk)
        Item = ArrListRemoveItem(SphereGradPath->SeedVectorY, Offset);

    if (IsOk)
        Item = ArrListRemoveItem(SphereGradPath->SeedVectorZ, Offset);

    if (IsOk)
    {
        GradPath_pa GradPath;
        Item = ArrListRemoveItem(SphereGradPath->GradPaths, Offset);
        GradPath = (GradPath_pa)Item.VoidPtr;
        GradPathDealloc(&GradPath);
    }

    ENSURE(SphereGradPathIsValid(SphereGradPath));

    ENSURE(VALID_BOOLEAN(IsOk));
    return IsOk;
}




/**
 * Appends the SeedVector/GradPath to the array lists in CircleGradPath. The
 * array lists will be expanded to accommodate the additional items.
 *
 * param SphereGradPath
 *     SphereGradPath target to which the SeedVector/GradPath is to be appended.
 * param Angle
 *     Angle to append to the SeedVector/GradPath.
 * param GradPath
 *     GradPath handle to append to the SphereGradPath.
 *
 * return
 *     TRUE if sufficient memory permitted the operation, otherwise FALSE.
 */
Boolean_t SphereGradPathAppendGP(SphereGradPath_pa  SphereGradPath,
                                 XYZ_s              SeedVector,
                                 GradPath_pa        GradPath)
{
    Boolean_t IsOk = TRUE;
    LgIndex_t Count;

    REQUIRE(SphereGradPathIsValid(SphereGradPath));

    Count = SphereGradPathGetGPCount(SphereGradPath);

    IsOk = SphereGradPathInsertGP(SphereGradPath, Count, SeedVector, GradPath);

    ENSURE(SphereGradPathIsValid(SphereGradPath));

    ENSURE(VALID_BOOLEAN(IsOk));
    return IsOk;
}









/**
 * Set the root TriGradPath at the specified offset.
 *
 *
 * param SphereGradPath
 *     SphereGradPath target in which to set the Angle/GradPath.
 * param Offset
 *     Offset into array lists.
 * param TriGradPath
 *     TriGradPath handle to set at the specified offset.
 *
 * return
 *     TRUE if sufficient memory permitted the operation, otherwise FALSE.
 */

Boolean_t   SphereGradPathSetRtTri(SphereGradPath_pa SphereGradPath,
                                   LgIndex_t         Offset,
                                   TriGradPath_pa    TriGradPath)
{
    Boolean_t IsOk = TRUE;

    REQUIRE(SphereGradPathIsValid(SphereGradPath));
    REQUIRE(Offset >= 0 && Offset < 20);

    SphereGradPath->RootTriangles[Offset] = TriGradPath;

    ENSURE(SphereGradPathIsValid(SphereGradPath));

    ENSURE(VALID_BOOLEAN(IsOk));
    return IsOk;
}








/**
 * Get the root TriGradPath at the specified offset.
 *
 *
 * param SphereGradPath
 *     SphereGradPath target in which to set the Angle/GradPath.
 * param Offset
 *     Offset into array lists.
 * param TriGradPath
 *     TriGradPath handle to set at the specified offset.
 *
 * return
 *     TRUE if sufficient memory permitted the operation, otherwise FALSE.
 */

TriGradPath_pa   SphereGradPathGetRtTri(SphereGradPath_pa SphereGradPath,
                                        LgIndex_t         Offset)
{
    TriGradPath_pa Result = NULL;

    REQUIRE(SphereGradPathIsValid(SphereGradPath));
    REQUIRE(Offset >= 0 && Offset < 20);

    if (SphereGradPath->RootTriangles[Offset] != NULL)
        Result = (TriGradPath_pa)SphereGradPath->RootTriangles[Offset];

    ENSURE(TriGradPathIsValid(Result));

    return Result;
}







/**
 * Return the ending critical point number for the GradPath at Offset
 * in the array list.
 *
 * param SphereGradPath
 *     SphereGradPath target to extract the EndCrtPtNum from.
 * param Offset
 *     Offset into array list of GradPath from which to extract EndCrtPtNum.
 *
 * return
 *     EndCrtPtNum for the path.
 */
LgIndex_t SphereGradPathGetEndCPNum(SphereGradPath_pa  SphereGradPath,
                                    LgIndex_t          Offset)
{
    LgIndex_t EndCrtPtNum = -2;

    GradPath_pa GradPath = NULL;
    ArrListItem_u Item;

    REQUIRE(SphereGradPathIsValid(SphereGradPath));
    REQUIRE(0 <= Offset && Offset < SphereGradPathGetGPCount(SphereGradPath));

    Item = ArrListGetItem(SphereGradPath->GradPaths, Offset);
    GradPath = (GradPath_pa)Item.VoidPtr;

    EndCrtPtNum = GradPath->EndCrtPtNum;

    ENSURE(-1 <= EndCrtPtNum);
    return EndCrtPtNum;
}




/**
 * Return the begining critical point number for the SphereGradPath.
 *
 * param SphereGradPath
 *     SphereGradPath target to extract the EndCrtPtNum from.
 *
 * return
 *     BegCrtPtNum for the path.
 */
LgIndex_t SphereGradPathGetBegCPNum(SphereGradPath_pa  SphereGradPath)
{
    LgIndex_t BegCrtPtNum = -2;
    LgIndex_t Offset = 0;

    GradPath_pa GradPath = NULL;
    ArrListItem_u Item;

    REQUIRE(SphereGradPathIsValid(SphereGradPath));

    Item = ArrListGetItem(SphereGradPath->GradPaths, Offset);
    GradPath = (GradPath_pa)Item.VoidPtr;

    BegCrtPtNum = GradPath->BeginCrtPtNum;

    ENSURE(-1 <= BegCrtPtNum);
    return BegCrtPtNum;
}







/**
 * Return the Offset into SphereGradPath->GradPaths of the GradPath with
 * the specified seed point.
 *
 * param SphereGradPath
 *     SphereGradPath target to extract the EndCrtPtNum from.
 * param Seed
 *     Seed location (compare to SphereGradPath->SeedVectorX, etc.).
 *
 * return
 *     Offset of desired GradPath in SphereGradPath->GradPaths list, or -1 if fails.
 */
LgIndex_t SphereGradPathGetGPBySeed(const SphereGradPath_pa  SphereGradPath,
                                    const XYZ_s              Seed)
{
    LgIndex_t GradPathOffset = -1;
    LgIndex_t Offset;
    Boolean_t Found = FALSE;
    LgIndex_t NumGradPaths = SphereGradPathGetGPCount(SphereGradPath);
    double    Tolerance = 1.0e-5;

    REQUIRE(SphereGradPathIsValid(SphereGradPath));

    for (Offset = 0; !Found && Offset < NumGradPaths; Offset++)
    {
        XYZ_s  GPSeed = SphereGradPathGetSeedDir(SphereGradPath, Offset);

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
 * Resample the CircleGradPath so that there are NPointsNew equally
 * spaced grid points along each GradPath length.
 *
 * param CircleGradPath
 *     CircleGradPath structure containing the GradPaths to be
 *     resampled.
 * param NumPointsNew
 *     Pointers to coordinates of the node
 *
 * return
 *     TRUE if it works, FALSE otherwise.
 */
/*
Boolean_t CircleGradPathResample(CircleGradPath_pa CircleGradPath,
                                 LgIndex_t         NumPointsNew)
{
  Boolean_t IsOk = TRUE;

  LgIndex_t   icgp;
  LgIndex_t   NumGradPaths = CircleGradPathGetCount(CircleGradPath);

  REQUIRE(CircleGradPathIsValid(CircleGradPath));
  REQUIRE(NumPointsNew > 1);


  for (icgp=0; IsOk && icgp<NumGradPaths; icgp++)
    {
      GradPath_pa GradPath = CircleGradPathGetGP(CircleGradPath, icgp);
      double      Angle = CircleGradPathGetAngle(CircleGradPath, icgp);

      IsOk = GradPathResample(&GradPath, NumPointsNew);

      IsOk = CircleGradPathSetGP(CircleGradPath, icgp, Angle, GradPath);
    }

  ENSURE(CircleGradPathIsValid(CircleGradPath));
  ENSURE(VALID_BOOLEAN(IsOk));
  return IsOk;
}
*/




/*
 * GradPathSphereAdd: Compute the set of gradient paths (streamlines) resulting
 *  from a seed sphere of radius Radius from point SphereCenter. Terminate each
 *  gradient path at a critical point or the boundary of the domain.
 */
/*
Boolean_t ComputeGradPathAtPtAngle(const ZoneVarInfo_pa  ZoneVarInfo,
                                   const double          XCrtPt,
                                   const double          YCrtPt,
                                   const double          ZCrtPt,
                                   const double          Angle,
                                   const double          n2x,
                                   const double          n2y,
                                   const double          n2z,
                                   const double          n3x,
                                   const double          n3y,
                                   const double          n3z,
                                   const CritPoints_pa   CritPoints,
                                   const StreamDir_e     PathDir,
                                   const float           CPTolerance,
                                   LgIndex_t   *NumPathPoints,
                                   GradPath_pa  GradPath)
{
  Boolean_t IsOk = TRUE;
  double xp, yp, XPos, YPos, ZPos;
  double junk = 0.0;

  REQUIRE(VALID_REF(ZoneVarInfo));
  REQUIRE(CritPointsIsValid(CritPoints));
  REQUIRE(VALID_ENUM(PathDir, StreamDir_e));
  REQUIRE(CPTolerance > 0.0);
  REQUIRE(VALID_REF(NumPathPoints));
  REQUIRE(GradPathIsValid(GradPath));

  *NumPathPoints = 0;

  xp = CPTolerance * cos(Angle);
  yp = CPTolerance * sin(Angle);

  XPos = XCrtPt + xp * n2x + yp * n3x;
  YPos = YCrtPt + xp * n2y + yp * n3y;
  ZPos = ZCrtPt + xp * n2z + yp * n3z;

  // Seed first point
  IsOk = GradPathSetPoint(GradPath, 0, XPos, YPos, ZPos, junk); // corret ChrgDens set later

  // Integrate gradient path line
  IsOk = GradPathAdd(ZoneVarInfo, MAXGRADPATHPOINTS, PathDir,
    CritPoints, CPTolerance, NumPathPoints, GradPath);

  ENSURE(*NumPathPoints > 0);
  ENSURE(GradPathIsValid(GradPath));
  ENSURE(VALID_BOOLEAN(IsOk));

  return IsOk;
}
*/

/*
 * IsValidStripFace: Check to see if the next triade in the FacetTriangles array
 * contains a unique set of nodes. If so, the strip continues across the face connecting
 * triangle Index and Index+1, and the two triangles are face neighbors.
 */
Boolean_t IsValidStripFace(LgIndex_t   Index,
                           LgIndex_t   SizeOfFacetArray,
                           LgIndex_t  *FacetTriangles)
{
    Boolean_t IsValid = TRUE;

    REQUIRE(Index >= 0 && Index < SizeOfFacetArray - 2);

    if (Index >= SizeOfFacetArray - 3 || Index < 0 ||
        FacetTriangles[Index+1] == FacetTriangles[Index+2] ||
        FacetTriangles[Index+2] == FacetTriangles[Index+3] ||
        FacetTriangles[Index+1] == FacetTriangles[Index+3]) IsValid = FALSE;

    ENSURE(VALID_BOOLEAN(IsValid));
    return IsValid;
}




/*
 * SphereGradPathAdd: Compute the set of gradient paths (streamlines) resulting
 *  from a seed sphere of radius CPTolerance from critical point BeginCrtPtNum.
 *  Terminate each gradient path at a critical point or the boundary of the domain.
 *  This function only populates the RootTriangles and the first twelve GradPaths.
 *  Subdivision of the root triangles is done elsewhere.
 *
 * param ZoneVarInfo
 *     Zone and variable information.
 * param BeginCrtPtNum
 *     Atom or Cage critical point that is the source of all grad paths in this
 *     SphereGradPath data structure
 * param CritPoints
 *     Critical Points structure
 * param PathDir
 *     Direction for pathline integration (Forward or Backward)
 * param CPTolerance
 *     Tolerance - distance from critical point at which the pathline is said
 *     to terminate at the critical point.
 * param SphereGradPath
 *     SphereGradPath data structure to populate for an Atom or Cage critical point.
 *
 * return
 *     TRUE if successful, FALSE if there were errors.
 */
Boolean_t SphereGradPathAdd(const ZoneVarInfo_pa  ZoneVarInfo,
                            const LgIndex_t       BeginCrtPtNum,
                            const CritPoints_pa   CritPoints,
                            const StreamDir_e     PathDir,
                            const float           CPTolerance,
                            SphereGradPath_pa     SphereGradPath)
{
    Boolean_t IsOk = TRUE;
    double XCrtPt, YCrtPt, ZCrtPt, dummy;
    char   cdummy;

    REQUIRE(VALID_REF(ZoneVarInfo));
    REQUIRE(VALID_REF(CritPoints));
    REQUIRE(BeginCrtPtNum >= 0 && BeginCrtPtNum < CritPointsGetCount(CritPoints));
    REQUIRE(CritPointsIsValid(CritPoints));
    REQUIRE(VALID_ENUM(PathDir, StreamDir_e));
    REQUIRE(CPTolerance > 0.0);

    // Extract Critical point location and principle direction
    IsOk = CritPointsGetPoint(CritPoints, BeginCrtPtNum,
                              &XCrtPt, &YCrtPt, &ZCrtPt, &dummy, &cdummy,
                              &dummy, &dummy, &dummy);


    //
    // Find seed points on a sphere surrounding the critical point.
    // Sphere is initially an icosohedron centered on the critical
    // point. Later, triangles may be locally refined.
    //
    if (IsOk)
    {
        int ii;

        for (ii = 0; IsOk && ii < 12; ii++)
        {
            double XSeed, YSeed, ZSeed;
            double junk = 0.0;
            LgIndex_t      NumPathPoints = 0;
            GradPath_pa    GradPath = NULL;

            XSeed = XCrtPt + CPTolerance * Sphere20FacetNodes[ii].X;
            YSeed = YCrtPt + CPTolerance * Sphere20FacetNodes[ii].Y;
            ZSeed = ZCrtPt + CPTolerance * Sphere20FacetNodes[ii].Z;

            GradPath = GradPathAlloc();

            GradPathClear(GradPath);

            /* Seed first point */
            GradPath->BeginCrtPtNum = BeginCrtPtNum;
            IsOk = GradPathSetPoint(GradPath, 0, XSeed, YSeed, ZSeed, junk); // corret ChrgDens set later

            /* Integrate gradient path line */
            if (IsOk)
                IsOk = GradPathAdd(ZoneVarInfo, MAXGRADPATHPOINTS, PathDir,
                                   CritPoints, CPTolerance, &NumPathPoints, GradPath);

            /* Add GradPath and Angle to CircleGradPath */
            if (IsOk)
                IsOk = SphereGradPathSetGP(SphereGradPath, ii, Sphere20FacetNodes[ii], GradPath);

        }
    }


    /* Set the root TriangleGradPath structures */
    if (IsOk)
    {
        int ii;
        int jj = 0;
        LgIndex_t FacetTriangArraySize = sizeof(Sphere20FacetTriangulation) / sizeof(LgIndex_t);

        for (ii = 0; IsOk && ii < 20; ii++)
        {
            TriGradPath_pa TriangleGP = SphereGradPath->RootTriangles[ii];
            TriangleGP->DepthInTree = 0;

            TriangleGP->GradPathNum[0] = Sphere20FacetTriangulation[jj];
            TriangleGP->GradPathNum[1] = Sphere20FacetTriangulation[jj+1];
            TriangleGP->GradPathNum[2] = Sphere20FacetTriangulation[jj+2];

            /* For a valid strip sequence, set neigboring cell numbers */
            if (IsValidStripFace(jj, FacetTriangArraySize, Sphere20FacetTriangulation))
            {
                TriangleGP->NeighborTri[0] = ii + 1;
                SphereGradPath->RootTriangles[ii+1]->NeighborTri[2] = ii;
            }

            /* Increment the pointer to the next triangle in the triangle strip array. */
            jj++;
            while ((Sphere20FacetTriangulation[jj  ] == Sphere20FacetTriangulation[jj+1] ||
                    Sphere20FacetTriangulation[jj+1] == Sphere20FacetTriangulation[jj+2]) &&
                   jj < FacetTriangArraySize - 2) jj++;
        }
    }

    /* Find and set the remaining neighbor triangle numbers */
    if (IsOk)
    {
        int Cell;

        for (Cell = 0; IsOk && Cell < 20; Cell++)
        {
            TriGradPath_pa TriangleGP = SphereGradPath->RootTriangles[Cell];
            int Face;

            for (Face = 0; Face < 3; Face++)
            {
                if (TriangleGP->NeighborTri[Face] < 0)
                {
                    Boolean_t Found = FALSE;
                    int       IndxNode1 = Face + 1;
                    int       IndxNode2 = Face + 2;
                    LgIndex_t Node1, Node2, TestNghbr;

                    if (IndxNode1 > 2) IndxNode1 -= 3;
                    if (IndxNode2 > 2) IndxNode2 -= 3;

                    Node1 = TriangleGP->GradPathNum[IndxNode1];
                    Node2 = TriangleGP->GradPathNum[IndxNode2];

                    for (TestNghbr = 0; !Found && TestNghbr < 20; TestNghbr++)
                    {
                        if (TestNghbr != Cell)
                        {
                            TriGradPath_pa TriTstNghbr = SphereGradPath->RootTriangles[TestNghbr];
                            int TstFace;
                            for (TstFace = 0; !Found && TstFace < 3; TstFace++)
                            {
                                int       IndxTstNd1 = TstFace + 1;
                                int       IndxTstNd2 = TstFace + 2;
                                LgIndex_t TstNode1, TstNode2;

                                if (TriTstNghbr->NeighborTri[TstFace] < 0)
                                {
                                    if (IndxTstNd1 > 2) IndxTstNd1 -= 3;
                                    if (IndxTstNd2 > 2) IndxTstNd2 -= 3;

                                    TstNode1 = TriTstNghbr->GradPathNum[IndxTstNd1];
                                    TstNode2 = TriTstNghbr->GradPathNum[IndxTstNd2];

                                    if ((TstNode1 == Node1 && TstNode2 == Node2) ||
                                        (TstNode1 == Node2 && TstNode2 == Node1))
                                    {
                                        Found = TRUE;
                                        TriangleGP->NeighborTri[Face]     = TestNghbr;
                                        TriTstNghbr->NeighborTri[TstFace] = Cell;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    ENSURE(VALID_BOOLEAN(IsOk));
    return IsOk;
}







Boolean_t         TriGradPathSubdivide(const ZoneVarInfo_pa  ZoneVarInfo,
                                       const CritPoints_pa   CritPoints,
                                       TriGradPath_pa        TriGradPath,
                                       const float           CPTolerance,
                                       SphereGradPath_pa     SphereGradPath)
{
    Boolean_t  IsOk  = TRUE;
    LgIndex_t  GradPathNums[3];

    REQUIRE(TriGradPathIsValid(TriGradPath));
    REQUIRE(SphereGradPathIsValid(SphereGradPath));

    /* Create new grad paths (or find them if they already exist) */
    if (IsOk)
    {
        int ii;

        for (ii = 0; IsOk && ii < 3; ii++)
        {
            GradPath_pa GradPath = NULL;
            int iip1 = ii + 1 - ((ii + 1) / 3) * 3;
            LgIndex_t GPNumII   = TriGradPathGetGPIndex(TriGradPath, ii);
            LgIndex_t GPNumIIP1 = TriGradPathGetGPIndex(TriGradPath, iip1);

            XYZ_s SeedII   = SphereGradPathGetSeedDir(SphereGradPath, GPNumII);
            XYZ_s SeedIIP1 = SphereGradPathGetSeedDir(SphereGradPath, GPNumIIP1);
            XYZ_s SeedNew;

            double Radius = sqrt(SeedII.X * SeedII.X + SeedII.Y * SeedII.Y + SeedII.Z * SeedII.Z);
            double RadiusNew;

            SeedNew.X = 0.5 * (SeedII.X + SeedIIP1.X);
            SeedNew.Y = 0.5 * (SeedII.Y + SeedIIP1.Y);
            SeedNew.Z = 0.5 * (SeedII.Z + SeedIIP1.Z);

            RadiusNew = sqrt(SeedNew.X * SeedNew.X + SeedNew.Y * SeedNew.Y + SeedNew.Z * SeedNew.Z);

            SeedNew.X = SeedNew.X * Radius / RadiusNew;
            SeedNew.Y = SeedNew.Y * Radius / RadiusNew;
            SeedNew.Z = SeedNew.Z * Radius / RadiusNew;

            GradPathNums[ii] = SphereGradPathGetGPBySeed(SphereGradPath, SeedNew);

            // If GradPath doesn't exist, create
            if (GradPathNums[ii] < 0)
            {
                GradPath = GradPathAlloc();

                if (GradPath != NULL)
                {
                    double junk = 0.0;
                    char   cjunk;
                    double XSeed, YSeed, ZSeed, XCrtPt, YCrtPt, ZCrtPt;
                    LgIndex_t BegCrtPtNum;
                    LgIndex_t NumPathPoints;
                    StreamDir_e PathDir = StreamDir_Reverse;  // assume atom

                    GradPathClear(GradPath);

                    /* Seed first point */
                    BegCrtPtNum = SphereGradPathGetBegCPNum(SphereGradPath);
                    IsOk = CritPointsGetPoint(CritPoints, BegCrtPtNum, &XCrtPt, &YCrtPt, &ZCrtPt,
                                              &junk, &cjunk, &junk, &junk, &junk);

                    GradPath->BeginCrtPtNum = BegCrtPtNum;
                    XSeed = XCrtPt + CPTolerance * SeedNew.X;
                    YSeed = YCrtPt + CPTolerance * SeedNew.Y;
                    ZSeed = ZCrtPt + CPTolerance * SeedNew.Z;

                    junk = 0.0;
                    IsOk = GradPathSetPoint(GradPath, 0, XSeed, YSeed, ZSeed, junk); // corret ChrgDens set later

                    /* Integrate gradient path line */
                    if (IsOk)
                        IsOk = GradPathAdd(ZoneVarInfo, MAXGRADPATHPOINTS, PathDir,
                                           CritPoints, 1.0, &NumPathPoints, GradPath);

                    /* Add GradPath and Angle to CircleGradPath */
                    if (IsOk)
                        IsOk = SphereGradPathAppendGP(SphereGradPath, SeedNew, GradPath);

                    if (IsOk)
                        GradPathNums[ii] = SphereGradPathGetGPCount(SphereGradPath) - 1;
                }
                else
                    IsOk = FALSE;
            }
        }
    }

    /* Create Sub Triangles */
    if (IsOk)
    {
        int ii;
        for (ii = 0; IsOk && ii < 4; ii++)
        {
            TriGradPath_pa SubTri = TriGradPathAlloc();
            if (TriGradPathIsValid(SubTri))
            {
                TriGradPathClear(SubTri);
                SubTri->DepthInTree = TriGradPath->DepthInTree + 1;
                TriGradPath->SubTriangles[ii] = SubTri;
            }
            else
                IsOk = FALSE;
        }
        for (ii = 0; IsOk && ii < 3; ii++)
        {
            int ThrdGP = ii + 2 - ((ii + 2) / 3) * 3;
            TriGradPath_pa SubTri = TriGradPath->SubTriangles[ii];
            SubTri->GradPathNum[0] = TriGradPath->GradPathNum[ii];
            SubTri->GradPathNum[1] = GradPathNums[ii];
            SubTri->GradPathNum[2] = GradPathNums[ThrdGP];
        }
        if (IsOk)
        {
            TriGradPath_pa SubTri = TriGradPath->SubTriangles[3];
            SubTri->GradPathNum[0] = GradPathNums[0];
            SubTri->GradPathNum[1] = GradPathNums[1];
            SubTri->GradPathNum[2] = GradPathNums[2];
        }
    }

    ENSURE(VALID_BOOLEAN(IsOk));
    return IsOk;
}




/*
 * SphereGradPathMinLengthGP: Find the minimum length gradient path to the
 *  specified EndCrtPtNum. Use the GradPaths in SphereGradPath between SGPBegin
 *  and SGPEnd.
 *
 * param SphereGradPath
 *     SphereGradPath data structure to populate for an Atom or Cage critical point.
 * param SGPBegin, SGPEnd
 *     Range of GradientPaths, in SphereGradPath, for testing
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
Boolean_t SphereGradPathMinLengthGP(const SphereGradPath_pa  SphereGradPath,
                                    const LgIndex_t          SGPBegin,
                                    const LgIndex_t          SGPEnd,
                                    const LgIndex_t          EndCrtPtNum,
                                    double                  *MinLength,
                                    LgIndex_t               *MinLenGPNum,
                                    LgIndex_t               *NumGPs)
{
    Boolean_t IsOk = TRUE;
    int ii;

    REQUIRE(SphereGradPathIsValid(SphereGradPath));
    REQUIRE(SGPBegin < SGPEnd);
    REQUIRE(EndCrtPtNum >= 0);
    REQUIRE(VALID_REF(MinLength));
    REQUIRE(*MinLength > 0.0);
    REQUIRE(VALID_REF(MinLenGPNum));
    REQUIRE(VALID_REF(NumGPs));
    REQUIRE(*NumGPs >= 0);

    if (SGPEnd > SGPBegin)
    {
        for (ii = SGPBegin; IsOk && ii < SGPEnd; ii++)
        {
            GradPath_pa GPTry  = SphereGradPathGetGP(SphereGradPath, ii);
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
 * SphereGradPathGet2CPPath: Find the minimum length gradient path between the
 *  specified BeginCrtPtNum and EndCrtPtNum. Use the GradPaths in SphereGradPath,
 *  and refine the seed point grid on the Sphere (using triangle subdivision) if
 *  necessary.
 *
 * param ZoneVarInfo
 *     Zone and variable information.
 * param BeginCrtPtNum
 *     Atom or Cage critical point that is the source of all grad paths in this
 *     SphereGradPath data structure
 * param EndCrtPtNum
 *     Atom or Cage critical point that is the termination of the desired grad paths.
 * param CritPoints
 *     Critical Points structure
 * param PathDir
 *     Direction for pathline integration (Forward or Reverse)
 * param CPTolerance
 *     Tolerance - distance from critical point at which the pathline is said
 *     to terminate at the critical point.
 * param SphereGradPath
 *     SphereGradPath data structure to populate for an Atom or Cage critical point.
 *
 * return
 *     TRUE if successful, FALSE if there were errors.
 */
GradPath_pa SphereGradPath2CPPath(const ZoneVarInfo_pa  ZoneVarInfo,
                                  const LgIndex_t       BeginCrtPtNum,
                                  const LgIndex_t       EndCrtPtNum,
                                  const CritPoints_pa   CritPoints,
                                  const StreamDir_e     PathDir,
                                  const float           CPTolerance,
                                  SphereGradPath_pa     SphereGradPath)
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
    REQUIRE(BeginCrtPtNum >= 0 && BeginCrtPtNum < CritPointsGetCount(CritPoints));
    REQUIRE(EndCrtPtNum >= 0 && EndCrtPtNum < CritPointsGetCount(CritPoints));
    REQUIRE(CritPointsIsValid(CritPoints));
    REQUIRE(VALID_ENUM(PathDir, StreamDir_e));
    REQUIRE(CPTolerance > 0.0);

    /*
     * Find the min-length GradPath from the current set that connects BeginCrtPt
     * and EndCrtPtNum
     */
    if (IsOk)
    {
        IsOk = SphereGradPathMinLengthGP(SphereGradPath, SGPCountOld, SGPCountNew, EndCrtPtNum,
                                         &MinLength, &MinLenGPNum, &NumGPs);
        if (IsOk && MinLenGPNum >= 0)
            Result = SphereGradPathGetGP(SphereGradPath, MinLenGPNum);

        RootMinLenGPNum = MinLenGPNum;
    }


    /*
     * If fewer than 5 GP's connect the two points, subdivide root triangles containing
     * the min-length connecting line.
     */
    if (IsOk)
    {
        LgIndex_t ii;

        for (ii = 0; IsOk && ii < 20; ii++)
        {
            TriGradPath_pa Triangle = SphereGradPathGetRtTri(SphereGradPath, ii);
            if (Triangle->GradPathNum[0] == MinLenGPNum ||
                Triangle->GradPathNum[1] == MinLenGPNum ||
                Triangle->GradPathNum[2] == MinLenGPNum)
            {
                IsOk = TriGradPathSubdivide(ZoneVarInfo, CritPoints, Triangle,
                                            CPTolerance, SphereGradPath);
            }
        }

        SGPCountOld = SGPCountNew;
        SGPCountNew = SphereGradPathGetGPCount(SphereGradPath);

        if (SGPCountNew > SGPCountOld)
        {
            IsOk = SphereGradPathMinLengthGP(SphereGradPath, SGPCountOld, SGPCountNew, EndCrtPtNum,
                                             &MinLength, &MinLenGPNum, &NumGPs);
            if (IsOk && MinLenGPNum >= SGPCountOld)
                Result = SphereGradPathGetGP(SphereGradPath, MinLenGPNum);
        }
    }

    /*
     * If fewer than 5 GP's connect the two points, subdivide the subtriangles containing
     * the min-length connecting line again.
     */
    if (IsOk)
    {
        LgIndex_t ii;

        for (ii = 0; IsOk && ii < 20; ii++)
        {
            TriGradPath_pa Triangle = SphereGradPathGetRtTri(SphereGradPath, ii);
            if (Triangle->GradPathNum[0] == RootMinLenGPNum ||
                Triangle->GradPathNum[1] == RootMinLenGPNum ||
                Triangle->GradPathNum[2] == RootMinLenGPNum)
            {
                LgIndex_t jj;
                for (jj = 0; IsOk && jj < 4; jj++)
                {
                    TriGradPath_pa SubTriangle = TriGradPathGetSubTri(Triangle, jj);
                    if (MinLenGPNum == RootMinLenGPNum ||
                        SubTriangle->GradPathNum[0] == MinLenGPNum ||
                        SubTriangle->GradPathNum[1] == MinLenGPNum ||
                        SubTriangle->GradPathNum[2] == MinLenGPNum)
                    {
                        IsOk = TriGradPathSubdivide(ZoneVarInfo, CritPoints, SubTriangle,
                                                    CPTolerance, SphereGradPath);

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
                            for (kk = 0; IsOk && kk < 6; kk++)
                            {
                                GradPaths[kk] = SphereGradPathGetGP(SphereGradPath, GradPathNums[kk]);
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
                                double A = 4.0 * (Len1 + Len3) - 8.0 * Len6;
                                double B = 4.0 * (Len3 + Len4 - Len5 - Len6);
                                double C = Len1 + 3.0 * Len3 - 4.0 * Len4;
                                double D = B;
                                double E = 4.0 * (Len2 + Len3) - 8.0 * Len5;
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
                                    AreaCd2 = (C * D - A * F) / Denom;
                                    AreaCd1 = (C - B * AreaCd2) / A;
                                    AreaCd3 = 1.0 - AreaCd1 - AreaCd2;
                                    NtCd[0] = AreaCd1 * (2.0 * AreaCd1 - 1.0);
                                    NtCd[1] = AreaCd2 * (2.0 * AreaCd2 - 1.0);
                                    NtCd[2] = AreaCd3 * (2.0 * AreaCd3 - 1.0);
                                    NtCd[3] = 4.0 * AreaCd1 * AreaCd2;
                                    NtCd[4] = 4.0 * AreaCd2 * AreaCd3;
                                    NtCd[5] = 4.0 * AreaCd3 * AreaCd1;
                                }
                                for (kk = 0; IsOk && kk < 6; kk++)
                                {
                                    XYZNode = SphereGradPathGetSeedDir(SphereGradPath, GradPathNums[kk]);
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
        SGPCountNew = SphereGradPathGetGPCount(SphereGradPath);

        if (SGPCountNew > SGPCountOld)
        {
            IsOk = SphereGradPathMinLengthGP(SphereGradPath, SGPCountOld, SGPCountNew, EndCrtPtNum,
                                             &MinLength, &MinLenGPNum, &NumGPs);
            if (IsOk && MinLenGPNum >= SGPCountOld)
                Result = SphereGradPathGetGP(SphereGradPath, MinLenGPNum);
        }
    }

    /*
     * If all points within a triangle and it's subtriangle end at the cage,
     * find the X,Y,Z seed that will give the minimum PathLength assuming that
     * PathLength varies quadratically across the triangle.
     * NOT CURRENTLY USED
     */
    if (0 && IsOk)
    {
        Boolean_t IsDone = FALSE;
        TriGradPath_pa TriPathStack[MAX_TRI_SUBDIVISIONS];
        EntIndex_t     SubTriNumStack[MAX_TRI_SUBDIVISIONS];

        /* Start at root and search down branches */
        TriGradPath_pa CurrentTri = NULL;
        EntIndex_t     CurSubTriNum = 0;
        EntIndex_t     CurTriNum = 0;
        EntIndex_t     Level = 0;
        while (!IsDone)
        {
            // Special treatment at Root level: Are we done? If not, set CurrentTri
            if (Level == 0)
            {
                if (CurTriNum > 19 || SphereGradPath->RootTriangles[CurTriNum] == NULL)
                    IsDone = TRUE;
                else
                {
                    CurrentTri = SphereGradPath->RootTriangles[CurSubTriNum];
                    CurSubTriNum = 0;
                }
            }

            // Move down a level
            if ((Level > 0 && CurSubTriNum < 4)  && CurrentTri->SubTriangles[CurSubTriNum] != NULL)
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
                if ((Level > 0 && CurSubTriNum < 3) || (Level == 0 && CurSubTriNum < 19))
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

        /* Are all six nodes of triangle and subtri's ending at the desired cage? */
    }

    if (!IsOk)
        Result = NULL;

    ENSURE(Result == NULL || GradPathIsValid(Result));

    return Result;
}