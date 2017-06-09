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
#include "ISOTOPO.h"
#include "SURFELEMMAP.h"
#include <string.h>



/**
 * Determine if the IsoTopo handle is sane.
 *
 * param IsoTopo
 *     IsoTopo structure in question.
 *
 * return
 *     TRUE if the IsoTopo structure is valid, otherwise FALSE.
 */
Boolean_t IsoTopoIsValid(IsoTopo_pa IsoTopo)
{
    Boolean_t IsValid = FALSE;

    IsValid = (VALID_REF(IsoTopo) &&
               VALID_REF(IsoTopo->IsoElemNums) && ArrListIsValid(IsoTopo->IsoElemNums));

    // Unit test some functions
    if (IsValid)
        IsValid = IsoTopoTestPointTriDist();

    ENSURE(VALID_BOOLEAN(IsValid));
    return IsValid;
}




/**
 * Deallocates the IsoTopo handle and set the handle to NULL.
 *
 * note
 *     Item destruction is the responsibility of the caller.
 *
 * param
 *     Reference to a IsoTopo handle.
 */
void IsoTopoDealloc(IsoTopo_pa *IsoTopo)
{
    REQUIRE(VALID_REF(IsoTopo));
    REQUIRE(IsoTopoIsValid(*IsoTopo) || *IsoTopo == NULL);

    if (*IsoTopo != NULL)
    {
        /* release the ArrList's */
        if ((*IsoTopo)->IsoElemNums != NULL) ArrListDealloc(&((*IsoTopo)->IsoElemNums));

        /* release the list structure itself */
        FREE_ITEM(*IsoTopo, "IsoTopo structure");
        *IsoTopo = NULL;
    }

    ENSURE(*IsoTopo == NULL);
}





/**
 * Empties the IsoTopo structure.
 *
 * param IsoTopo
 *     IsoTopo to clear.
 */
void IsoTopoClear(IsoTopo_pa IsoTopo)
{
    REQUIRE(IsoTopoIsValid(IsoTopo));

    ArrListClear(IsoTopo->IsoElemNums);

    ENSURE(IsoTopoIsValid(IsoTopo) && IsoTopoGetElemCount(IsoTopo) == 0);
}






/**
 * Allocates a IsoTopo handle with a suitable default
 * capacity for the contained ArrList's.
 *
 * return
 *     IsoTopo handle if sufficient memory was available,
 *     otherewise a handle to NULL.
 */
IsoTopo_pa IsoTopoAlloc()
{
    IsoTopo_pa Result = NULL;

    Result = ALLOC_ITEM(IsoTopo_s, "IsoTopo structure");
    if (Result != NULL)
    {
        Result->IsoElemNums  = ArrListAlloc(20, ArrListType_Long);

        /* If it failed to allocate any of the array lists, clean-up and exit. */
        if (Result->IsoElemNums == NULL)
        {
            FREE_ITEM(Result, "IsoTopo structure");
            Result = NULL;
        }
    }

    ENSURE(IsoTopoIsValid(Result) || Result == NULL);
    return Result;
}




/**
 * Gets the number of Elements currently in the IsoTopo structure.
 *
 * param
 *     IsoTopo structure in question.
 *
 * return
 *     Number of elements in the IsoTopo structure.
 */
LgIndex_t IsoTopoGetElemCount(const IsoTopo_pa IsoTopo)
{
    LgIndex_t Result = 0;

    REQUIRE(IsoTopoIsValid(IsoTopo));

    Result = ArrListGetCount(IsoTopo->IsoElemNums);

    ENSURE(Result >= 0);
    return Result;
}






/**
 * Places IsoTopo point coordinates and variable value at
 * the specified offset. If the offset is beyond the end of the
 * array lists, they are sized accordingly and the intervening
 * array values between the last node of the original state and
 * the last node of the new state are assigned 0.
 *
 * note
 *     If components already exists at the specified location they are
 *     replaced.
 *
 * param IsoTopo
 *     IsoTopo target in which to set the coordinates and var value.
 * param PointOffset
 *     Offset of the point/node.
 * param Coord1, Coord2, Var
 *     Coordinates and var value to set at the specified offset.
 *
 * return
 *     TRUE if sufficient memory permitted the operation, otherwise FALSE.
 */

Boolean_t IsoTopoSetPointAtOffset(IsoTopo_pa IsoTopo,
                                  LgIndex_t  ElemOffset,
                                  LgIndex_t  IsoElemNum)
{
    Boolean_t IsOk = TRUE;
    ArrListItem_u Item;


    REQUIRE(IsoTopoIsValid(IsoTopo));
    REQUIRE(ElemOffset >= 0);

    Item.Long = IsoElemNum;
    IsOk = ArrListSetItem(IsoTopo->IsoElemNums, ElemOffset, Item);

    ENSURE(IsoTopoIsValid(IsoTopo));

    ENSURE(VALID_BOOLEAN(IsOk));
    return IsOk;
}





/**
 * Inserts IsoTopo element number at the specified
 * offset. The arrays will be expanded to accomodate the additional value.
 *
 *
 * param IsoTopo
 *     IsoTopo target in which to set the coordinates.
 * param ElemOffset
 *     Offset at which to insert the IsoElemNum.
 * param IsoElemNum
 *     Number of an element in the Isosurface zone to set at the specified offset.
 *
 * return
 *     TRUE if sufficient memory permitted the operation, otherwise FALSE.
 */

Boolean_t IsoTopoInsertElemAtOffset(IsoTopo_pa IsoTopo,
                                    LgIndex_t  ElemOffset,
                                    LgIndex_t  IsoElemNum)
{
    Boolean_t IsOk = TRUE;
    ArrListItem_u Item;


    REQUIRE(IsoTopoIsValid(IsoTopo));
    REQUIRE(0 <= ElemOffset && ElemOffset <= IsoTopoGetElemCount(IsoTopo));

    Item.Long = IsoElemNum;
    IsOk = ArrListInsertItem(IsoTopo->IsoElemNums, ElemOffset, Item);

    ENSURE(IsoTopoIsValid(IsoTopo));

    ENSURE(VALID_BOOLEAN(IsOk));
    return IsOk;
}





/**
 * Appends IsoElemNum to IsoTopo structure. The array lists
 * will be expanded to accommodate the additional items.
 *
 * param IsoTopo
 *     IsoTopo target to which the IsoElemNum is to be appended.
 * param IsoElemNum
 *     Element number, in the isosurface zone, to be appended to IsoTopo.
 *
 * return
 *     TRUE if sufficient memory permitted the operation, otherwise FALSE.
 */
Boolean_t IsoTopoAppendElemAtEnd(IsoTopo_pa  IsoTopo,
                                 LgIndex_t   IsoElemNum)
{
    Boolean_t IsOk = TRUE;
    LgIndex_t Count;

    REQUIRE(IsoTopoIsValid(IsoTopo));

    Count = IsoTopoGetElemCount(IsoTopo);

    IsOk = IsoTopoInsertElemAtOffset(IsoTopo, Count, IsoElemNum);

    ENSURE(IsoTopoIsValid(IsoTopo));

    ENSURE(VALID_BOOLEAN(IsOk));
    return IsOk;
}





/**
 * Gets the IsoElemNum from the IsoTopo structure for the
 * element at the specified offset.
 *
 * param IsoTopo
 *     IsoTopo structure containing the desired item.
 * param ElemOffset
 *     Offset to the IsoElemNum of the IsoTopo.
 *
 * return
 *     IsoElemNum if it works, -1 otherwise.
 */

LgIndex_t IsoTopoGetElemNumFromOffset(const IsoTopo_pa  IsoTopo,
                                      const LgIndex_t   ElemOffset)
{
    LgIndex_t IsoElemNum = -1;
    ArrListItem_u Item;

    REQUIRE(IsoTopoIsValid(IsoTopo));
    REQUIRE(ElemOffset >= 0 && ElemOffset < IsoTopoGetElemCount(IsoTopo));

    Item = ArrListGetItem(IsoTopo->IsoElemNums, ElemOffset);
    IsoElemNum = Item.Long;

    ENSURE(IsoElemNum >= -1);
    return IsoElemNum;
}

/**
 * Vector Cross product
 * param V1, V2
 *     Vectors to be crossed
 * return
 *     Vector cross product
 */
static XYZ_s VectorCrossProduct(XYZ_s V1,
                                XYZ_s V2)
{
    XYZ_s Result;
    Result.X = V1.Y * V2.Z - V1.Z * V2.Y;
    Result.Y = V1.Z * V2.X - V1.X * V2.Z;
    Result.Z = V1.X * V2.Y - V1.Y * V2.X;
    return Result;
}






/**
 * Find the distance from a point to a triangular element
 * param P
 *     X, Y, Z coordinates of the point.
 * param Nd1, Nd2, Nd3
 *     X, Y, Z coordinates of nodes 1, 2, and 3 of triangle.
 *
 * return
 *     Smallest of normal, closest-edge, or closest-node distances.
 */
double PointTriangleDistance(XYZ_s  P,
                             XYZ_s  Nd1,
                             XYZ_s  Nd2,
                             XYZ_s  Nd3)
{
    double Dist;

    XYZ_s  E1, E2, E3;   // Edge Vectors
    double E1Length, E2Length, E3Length;
    XYZ_s  Np;
    double NpLength;
    XYZ_s  P0;  // point projected normally to plane of triangle
    XYZ_s  DP10, DP20, DP30;
    XYZ_s  E1CrossDP10, E2CrossDP20, E3CrossDP30;
    double DistNorm;  // normal distance from point to plane of triangle
    double DistLat = LARGEFLOAT;
    Boolean_t IsOutside = FALSE;

    // Compute Edge Vectors
    E1.X = Nd2.X - Nd1.X;
    E1.Y = Nd2.Y - Nd1.Y;
    E1.Z = Nd2.Z - Nd1.Z;
    E1Length = sqrt(E1.X * E1.X + E1.Y * E1.Y + E1.Z * E1.Z);
    E2.X = Nd3.X - Nd2.X;
    E2.Y = Nd3.Y - Nd2.Y;
    E2.Z = Nd3.Z - Nd2.Z;
    E2Length = sqrt(E2.X * E2.X + E2.Y * E2.Y + E2.Z * E2.Z);
    E3.X = Nd1.X - Nd3.X;
    E3.Y = Nd1.Y - Nd3.Y;
    E3.Z = Nd1.Z - Nd3.Z;
    E3Length = sqrt(E3.X * E3.X + E3.Y * E3.Y + E3.Z * E3.Z);

    // Compute the normal vector to the plane of the triangle
    Np = VectorCrossProduct(E1, E2);
    NpLength = sqrt(Np.X * Np.X + Np.Y * Np.Y + Np.Z * Np.Z);
    Np.X /= NpLength;
    Np.Y /= NpLength;
    Np.Z /= NpLength;
    NpLength = 1.0;

    // Project the point normally to the plane of the triangle
    DistNorm = (P.X - Nd1.X) * Np.X + (P.Y - Nd1.Y) * Np.Y + (P.Z - Nd1.Z) * Np.Z;
    P0.X = P.X - DistNorm * Np.X;
    P0.Y = P.Y - DistNorm * Np.Y;
    P0.Z = P.Z - DistNorm * Np.Z;

    // Compute vectors from each of the nodes to P0
    DP10.X = P0.X - Nd1.X;
    DP10.Y = P0.Y - Nd1.Y;
    DP10.Z = P0.Z - Nd1.Z;
    DP20.X = P0.X - Nd2.X;
    DP20.Y = P0.Y - Nd2.Y;
    DP20.Z = P0.Z - Nd2.Z;
    DP30.X = P0.X - Nd3.X;
    DP30.Y = P0.Y - Nd3.Y;
    DP30.Z = P0.Z - Nd3.Z;

    // Lateral distance, if outside edge 1
    E1CrossDP10 = VectorCrossProduct(E1, DP10);
    if ((E1CrossDP10.X * Np.X + E1CrossDP10.Y * Np.Y + E1CrossDP10.Z * Np.Z) < 0.0)    // Outside
    {
        double DistAlongEdge;
        DistAlongEdge = (DP10.X * E1.X + DP10.Y * E1.Y + DP10.Z * E1.Z) / E1Length;

        IsOutside = TRUE;

        // Before beginning of edge: take distance between P and Nd1
        if (DistAlongEdge < 0.0)
            /*
            {
              double dx = P.X - Nd1.X;
              double dy = P.Y - Nd1.Y;
              double dz = P.Z - Nd1.Z;
              DistLat = sqrt(dx * dx + dy * dy + dz * dz);
            } */
            DistLat = sqrt(DP10.X * DP10.X + DP10.Y * DP10.Y + DP10.Z * DP10.Z);
        // After end of edge: take distance between P and Nd2
        else if (DistAlongEdge > 1.0)
            /*
            {
              double dx = P.X - Nd2.X;
              double dy = P.Y - Nd2.Y;
              double dz = P.Z - Nd2.Z;
              DistLat = sqrt(dx * dx + dy * dy + dz * dz);
            } */
            DistLat = sqrt(DP20.X * DP20.X + DP20.Y * DP20.Y + DP20.Z * DP20.Z);
        // Lies along edge - find normal distance from edge
        else
        {
            XYZ_s E1Norm;  // Normal to edge (orthogonal to surface normal)
            E1Norm = VectorCrossProduct(E1, Np);
            E1Norm.X /= E1Length;
            E1Norm.Y /= E1Length;
            E1Norm.Z /= E1Length;
            DistLat = DP10.X * E1Norm.X + DP10.Y * E1Norm.Y + DP10.Z * E1Norm.Z;
        }
    }

    // Lateral distance, if outside edge 2
    E2CrossDP20 = VectorCrossProduct(E2, DP20);
    if ((E2CrossDP20.X * Np.X + E2CrossDP20.Y * Np.Y + E2CrossDP20.Z * Np.Z) < 0.0)    // Outside
    {
        double DistAlongEdge;
        DistAlongEdge = (DP20.X * E2.X + DP20.Y * E2.Y + DP20.Z * E2.Z) / E2Length;

        IsOutside = TRUE;

        // Before beginning of edge: take distance between P and Nd2
        if (DistAlongEdge < 0.0)
            /*
            {
              double dx = P.X - Nd2.X;
              double dy = P.Y - Nd2.Y;
              double dz = P.Z - Nd2.Z;
              DistLat = MIN(DistLat, sqrt(dx * dx + dy * dy + dz * dz));
            } */
            DistLat = MIN(DistLat, sqrt(DP20.X * DP20.X + DP20.Y * DP20.Y + DP20.Z * DP20.Z));
        // After end of edge: take distance between P and Nd3
        else if (DistAlongEdge > 1.0)
            /*
            {
              double dx = P.X - Nd3.X;
              double dy = P.Y - Nd3.Y;
              double dz = P.Z - Nd3.Z;
              DistLat = MIN(DistLat, sqrt(dx * dx + dy * dy + dz * dz));
            } */
            DistLat = MIN(DistLat, sqrt(DP30.X * DP30.X + DP30.Y * DP30.Y + DP30.Z * DP30.Z));
        // Lies along edge - find normal distance from edge
        else
        {
            XYZ_s E2Norm;  // Normal to edge (orthogonal to surface normal)
            E2Norm = VectorCrossProduct(E2, Np);
            E2Norm.X /= E2Length;
            E2Norm.Y /= E2Length;
            E2Norm.Z /= E2Length;
            DistLat = DP20.X * E2Norm.X + DP20.Y * E2Norm.Y + DP20.Z * E2Norm.Z;
        }
    }


    // Lateral distance, if outside edge 3
    E3CrossDP30 = VectorCrossProduct(E3, DP30);
    if ((E3CrossDP30.X * Np.X + E3CrossDP30.Y * Np.Y + E3CrossDP30.Z * Np.Z) < 0.0)    // Outside
    {
        double DistAlongEdge;
        DistAlongEdge = (DP30.X * E3.X + DP30.Y * E3.Y + DP30.Z * E3.Z) / E3Length;

        IsOutside = TRUE;

        // Before beginning of edge: take distance between P and Nd3
        if (DistAlongEdge < 0.0)
            /*
            {
              double dx = P.X - Nd3.X;
              double dy = P.Y - Nd3.Y;
              double dz = P.Z - Nd3.Z;
              DistLat = MIN(DistLat, sqrt(dx * dx + dy * dy + dz * dz));
            } */
            DistLat = MIN(DistLat, sqrt(DP30.X * DP30.X + DP30.Y * DP30.Y + DP30.Z * DP30.Z));
        // After end of edge: take distance between P and Nd1
        else if (DistAlongEdge > 1.0)
            /*
            {
              double dx = P.X - Nd1.X;
              double dy = P.Y - Nd1.Y;
              double dz = P.Z - Nd1.Z;
              DistLat = MIN(DistLat, sqrt(dx * dx + dy * dy + dz * dz));
            } */
            DistLat = MIN(DistLat, sqrt(DP10.X * DP10.X + DP10.Y * DP10.Y + DP10.Z * DP10.Z));
        // Lies along edge - find normal distance from edge
        else
        {
            XYZ_s E3Norm;  // Normal to edge (orthogonal to surface normal)
            E3Norm = VectorCrossProduct(E3, Np);
            E3Norm.X /= E3Length;
            E3Norm.Y /= E3Length;
            E3Norm.Z /= E3Length;
            DistLat = DP30.X * E3Norm.X + DP30.Y * E3Norm.Y + DP30.Z * E3Norm.Z;
        }
    }

    if (IsOutside)
        Dist = sqrt(DistNorm * DistNorm + DistLat * DistLat);
    else
        Dist = ABS(DistNorm);

    ENSURE(Dist >= 0.0);
    return Dist;
}



/**
 * Unit test PointTriangleDistance - the function to calculate
 * distance from a point to a triangular element. Test each of
 * the seven conditions for which PointTriangleDistance computes
 * distance: smallest of normal, closest-edge, or closest-node
 * distances.
 *
 * param - none -
 *
 * return
 *     TRUE if the test passes, FALSE if it fails..
 */
Boolean_t IsoTopoTestPointTriDist()
{
    Boolean_t Result = TRUE;
    XYZ_s     Point, Nd1, Nd2, Nd3;
    double    Dist, Diff;

    Nd1.X = 1.0;
    Nd1.Y = 0.0;
    Nd1.Z = 0.0;

    Nd2.X = 0.0;
    Nd2.Y = 1.0;
    Nd2.Z = 0.0;

    Nd3.X = 0.0;
    Nd3.Y = 0.0;
    Nd3.Z = 1.0;


    // Test Normal distance calculation
    Point.X = 1.0;
    Point.Y = 1.0;
    Point.Z = 1.0;

    Dist = PointTriangleDistance(Point, Nd1, Nd2, Nd3);

    Diff = Dist - 2.0 * sqrt(3.0) / 3.0;
    if (ABS(Diff) > 1.0e-5) Result = FALSE;


    // Test Edge 1 (Nd1 - Nd2) distance calculation
    if (Result)
    {
        Point.X = 1.0;
        Point.Y = 1.0;
        Point.Z = 0.0;

        Dist = PointTriangleDistance(Point, Nd1, Nd2, Nd3);

        Diff = Dist - sqrt(2.0) / 2.0;
        if (ABS(Diff) > 1.0e-5) Result = FALSE;
    }


    // Test Edge 2 (Nd2 - Nd3) distance calculation
    if (Result)
    {
        Point.X = 0.0;
        Point.Y = 1.0;
        Point.Z = 1.0;

        Dist = PointTriangleDistance(Point, Nd1, Nd2, Nd3);

        Diff = Dist - sqrt(2.0) / 2.0;
        if (ABS(Diff) > 1.0e-5) Result = FALSE;
    }


    // Test Edge 3 (Nd3 - Nd1) distance calculation
    if (Result)
    {
        Point.X = 1.0;
        Point.Y = 0.0;
        Point.Z = 1.0;

        Dist = PointTriangleDistance(Point, Nd1, Nd2, Nd3);

        Diff = Dist - sqrt(2.0) / 2.0;
        if (ABS(Diff) > 1.0e-5) Result = FALSE;
    }



    // Test Node 1 distance calculation
    if (Result)
    {
        Point.X = 2.0;
        Point.Y = 0.0;
        Point.Z = 0.0;

        Dist = PointTriangleDistance(Point, Nd1, Nd2, Nd3);

        if (ABS(Dist - 1.0) > 1.0e-5) Result = FALSE;
    }


    // Test Node 2 distance calculation
    if (Result)
    {
        Point.X = 0.0;
        Point.Y = 2.0;
        Point.Z = 0.0;

        Dist = PointTriangleDistance(Point, Nd1, Nd2, Nd3);

        if (ABS(Dist - 1.0) > 1.0e-5) Result = FALSE;
    }



    // Test Node 3 distance calculation
    if (Result)
    {
        Point.X = 0.0;
        Point.Y = 0.0;
        Point.Z = 2.0;

        Dist = PointTriangleDistance(Point, Nd1, Nd2, Nd3);

        if (ABS(Dist - 1.0) > 1.0e-5) Result = FALSE;
    }

    ENSURE(VALID_BOOLEAN(Result));
    return Result;
}



/**
 * Append an Item to the array list only if the item doesn't already
 * exist in the array list.
 *
 * param ArrList
 *     Array list to which the item is to be appended.
 * param Item
 *     Item to be appended (if unique).
 *
 * return
 *     TRUE if it was unique and appended, FALSE otherwise.
 */
Boolean_t ArrListAppendUniqueItem(ArrList_pa    ArrList,
                                  ArrListItem_u Item)
{
    Boolean_t IsUnique = TRUE;
    Boolean_t IsOk = TRUE;
    ArrListItem_u ItemOld;
    int ii;
    LgIndex_t NumItems = ArrListGetCount(ArrList);

    REQUIRE(ArrListIsValid(ArrList));
    REQUIRE(ArrListGetType(ArrList) == ArrListType_Long);  // Limitation for now

    // Check to see if it is unique
    for (ii = 0; IsUnique && ii < NumItems; ii++)
    {
        ItemOld = ArrListGetItem(ArrList, ii);
        if (Item.Long == ItemOld.Long) IsUnique = FALSE;
    }

    if (IsUnique)
        IsOk = ArrListAppendItem(ArrList, Item);

    // Unique but append failed (unlikely)
    CHECK(IsOk);
    if (!IsOk) IsUnique = FALSE;

    ENSURE(VALID_BOOLEAN(IsUnique));
    return IsUnique;
}




/**
 * Extract the topologically independent segment of the specified
 * extracted isosurface zone which is closest to the specified coordinates.
 *
 * param IsoTopo
 *     IsoTopo structure containing the input points for the curve fit.
 * param IsoSurfZone
 *     Number of the extraced isosurface zone.
 * param X, Y, Z
 *     Location of the point defining the desired isosurface segment. Return
 *     the topologically independent segment of the isosurface closest to
 *     this point.
 *
 * return
 *     Zone number of the new isosurface segment (if it worked), zero if it failed.
 */
EntIndex_t IsoTopoCompute(IsoTopo_pa IsoTopo,
                          EntIndex_t IsoSurfZone,
                          double     X,
                          double     Y,
                          double     Z)
{
    EntIndex_t IsoTopoZone = 0;
    EntIndex_t NumZones, NumVars;
    Boolean_t  IsOk = TRUE;
    LgIndex_t  NumNodes, NumElems;
    EntIndex_t ClosestElem;

    REQUIRE(IsoTopoIsValid(IsoTopo));

    /* Get num zones in dataset */
    if (IsOk)
        IsOk = TecUtilDataSetGetInfo(NULL, &NumZones, &NumVars);

    REQUIRE(IsoSurfZone > 0 && IsoSurfZone <= NumZones);


    // Find the closest element to the point
    if (IsOk)
    {
        LgIndex_t    ne;
        EntIndex_t   XVarNum, YVarNum, ZVarNum;
        FieldData_pa XVarFDPtr, YVarFDPtr, ZVarFDPtr;
        double       MinDist = LARGEFLOAT;
        ZoneType_e   ZoneType;
        Boolean_t    IsQuad;
        XYZ_s        Point, Nd1, Nd2, Nd3, Nd4;

        Point.X = X;
        Point.Y = Y;
        Point.Z = Z;

        // Find the number of Elements in the IsoSurfZone
        TecUtilZoneGetInfo(IsoSurfZone, &NumNodes, &NumElems,
                           NULL, NULL, NULL, NULL, NULL, NULL,
                           NULL, NULL, NULL, NULL, NULL);

        // FETriangle or FEQuad?
        ZoneType = TecUtilZoneGetType(IsoSurfZone);
        if (ZoneType = ZoneType_FEQuad)
            IsQuad = TRUE;
        else if (ZoneType = ZoneType_FETriangle)
            IsQuad = FALSE;
        else
            IsOk = FALSE;

        // Get X, Y, Z variable numbers
        XVarNum = TecUtilVarGetNumByAssignment('X');
        YVarNum = TecUtilVarGetNumByAssignment('Y');
        ZVarNum = TecUtilVarGetNumByAssignment('Z');

        if (XVarNum > 0 && XVarNum <= NumVars)
            XVarFDPtr = TecUtilDataValueGetReadableRef(IsoSurfZone, XVarNum);
        else
            IsOk = FALSE;

        if (YVarNum > 0 && YVarNum <= NumVars)
            YVarFDPtr = TecUtilDataValueGetReadableRef(IsoSurfZone, YVarNum);
        else
            IsOk = FALSE;

        if (YVarNum > 0 && YVarNum <= NumVars)
            ZVarFDPtr = TecUtilDataValueGetReadableRef(IsoSurfZone, ZVarNum);
        else
            IsOk = FALSE;

        // Loop over elements of IsoSurfZone to find the one closest to X,Y,Z
        for (ne = 1; IsOk && ne < NumElems; ne++)
        {
            double Dist;
            int ii, ic;
            LgIndex_t NdNums[4], NdNumsTmp[4];
            LgIndex_t Nd1Num, Nd2Num, Nd3Num, Nd4Num;

            // Find unique node numbers
            for (ii = 0; ii < 4; ii++)
                NdNumsTmp[ii] = TecUtilDataNodeGetByZone(IsoSurfZone, ne, ii + 1);

            NdNums[0] = NdNumsTmp[0];
            ic = 1;
            for (ii = 1; ii < 4; ii++)
            {
                int jj;
                Boolean_t IsUnique = TRUE;
                for (jj = 0; jj < ic; jj++)
                    if (NdNumsTmp[ii] == NdNums[jj]) IsUnique = FALSE;

                if (IsUnique)
                {
                    NdNums[ic] = NdNumsTmp[ii];
                    ic++;
                }
            }
            if (ic == 3) NdNums[3] = NdNums[2];

            Nd1Num = NdNums[0];
            Nd1.X = TecUtilDataValueGetByRef(XVarFDPtr, Nd1Num);
            Nd1.Y = TecUtilDataValueGetByRef(YVarFDPtr, Nd1Num);
            Nd1.Z = TecUtilDataValueGetByRef(ZVarFDPtr, Nd1Num);

            Nd2Num = NdNums[1];
            Nd2.X = TecUtilDataValueGetByRef(XVarFDPtr, Nd2Num);
            Nd2.Y = TecUtilDataValueGetByRef(YVarFDPtr, Nd2Num);
            Nd2.Z = TecUtilDataValueGetByRef(ZVarFDPtr, Nd2Num);

            Nd3Num = NdNums[2];
            Nd3.X = TecUtilDataValueGetByRef(XVarFDPtr, Nd3Num);
            Nd3.Y = TecUtilDataValueGetByRef(YVarFDPtr, Nd3Num);
            Nd3.Z = TecUtilDataValueGetByRef(ZVarFDPtr, Nd3Num);

            Dist = PointTriangleDistance(Point, Nd1, Nd2, Nd3);
            if (Dist < MinDist)
            {
                MinDist = Dist;
                ClosestElem = ne;
            }

            // If quad, do the second triangle (nodes 1, 3, and 4)
            Nd4Num = NdNums[3];
            if (IsQuad && Nd1Num != Nd4Num && Nd4Num != Nd3Num)
            {
                Nd4.X = TecUtilDataValueGetByRef(XVarFDPtr, Nd4Num);
                Nd4.Y = TecUtilDataValueGetByRef(YVarFDPtr, Nd4Num);
                Nd4.Z = TecUtilDataValueGetByRef(ZVarFDPtr, Nd4Num);
                Dist = PointTriangleDistance(Point, Nd1, Nd3, Nd4);
                if (Dist < MinDist)
                {
                    MinDist = Dist;
                    ClosestElem = ne;
                }
            }
        }
        // Set the closest element number into the element list.
        IsOk = IsoTopoAppendElemAtEnd(IsoTopo, ClosestElem);
    }

    // Unit test the SurfElemMap module
    CHECK(SurfElemMapTest());


    // Start with the closest cell and add all connected cells
    if (IsOk)
    {
        int nn;
        int NodesPerElem = 3;
        ZoneType_e ZoneType;
        SurfElemMap_pa IsoElemMap = NULL;
        ArrList_pa ElemList = NULL;
        ArrList_pa NodeBoundaryList = NULL;
        ArrList_pa NodeList = NULL;
        ArrListItem_u Item;

        ZoneType = TecUtilZoneGetType(IsoSurfZone);
        if (ZoneType = ZoneType_FEQuad)
            NodesPerElem = 4;
        else if (ZoneType = ZoneType_FETriangle)
            NodesPerElem = 3;
        else
            IsOk = FALSE;

        // Compute the element map of the extracted isosurface
        IsoElemMap = SurfElemMapAlloc(NumNodes);
        IsOk = SurfElemMapCompute(IsoElemMap, IsoSurfZone);

        // Add closest cell and initialize NodeBoundaryList with it's nodes
        if (IsOk)
        {
            ElemList = ArrListAlloc(60, ArrListType_Long);
            NodeBoundaryList = ArrListAlloc(60, ArrListType_Long);
            NodeList = ArrListAlloc(60, ArrListType_Long);
            // Add first element to ElemList
            Item.Long = ClosestElem;
            IsOk = ArrListAppendItem(ElemList, Item);

            // Add unique node number to the NodeBoundaryList
            Item.Long = TecUtilDataNodeGetByZone(IsoSurfZone, ClosestElem, 1);
            IsOk = ArrListAppendItem(NodeBoundaryList, Item);
            IsOk = ArrListAppendItem(NodeList, Item);

            for (nn = 2; nn <= NodesPerElem; nn++)
            {
                Boolean_t IsUnique;
                Item.Long = TecUtilDataNodeGetByZone(IsoSurfZone, ClosestElem, nn);
                IsUnique = ArrListAppendUniqueItem(NodeBoundaryList, Item);
                IsUnique = ArrListAppendUniqueItem(NodeList, Item);
            }
        }

        // Loop over nodes in NodeBoundaryList, each time adding the unique elements
        // that contain the node to the ElemList and adding the unique nodes from
        // those elements to the NodeBoudaryList. Continue until there are no more
        // nodes in the NodeBoundaryList. Only works for 2-manifold surfaces.
        while (IsOk && ArrListGetCount(NodeBoundaryList) > 0)
        {
            int ne;
            LgIndex_t NodeToRemove;
            LgIndex_t NumElemsWithNode;

            // Select first node in list
            Item = ArrListGetItem(NodeBoundaryList, 0);
            NodeToRemove = Item.Long;

            // Cycle through elements containing node
            NumElemsWithNode = SurfElemMapGetElemCountForNode(IsoElemMap, NodeToRemove);

            for (ne = 0; IsOk && ne < NumElemsWithNode; ne++)
            {
                Boolean_t IsUnique = TRUE;
                LgIndex_t NewElemNum = SurfElemMapGetElemNumForNodeOffset(IsoElemMap, NodeToRemove, ne);
                Item.Long = NewElemNum;
                IsUnique = ArrListAppendUniqueItem(ElemList, Item);

                // For unique new element, add the elements unique nodes to NodeBoundaryList
                if (IsUnique)
                {
                    int nn;
                    for (nn = 1; nn <= NodesPerElem; nn++)
                    {
                        Boolean_t IsNodeUnique;
                        Item.Long = TecUtilDataNodeGetByZone(IsoSurfZone, NewElemNum, nn);
                        IsNodeUnique = ArrListAppendUniqueItem(NodeBoundaryList, Item);
                        IsNodeUnique = ArrListAppendUniqueItem(NodeList, Item);
                    }
                }
            }
            // Remove the NodeToRemove from NodeBoundaryList
            Item = ArrListRemoveItem(NodeBoundaryList, 0);
        }

        if (IsOk)
        {
            // Create the new Tecplot zone
            LgIndex_t ITZNumNodes = ArrListGetCount(NodeList);
            LgIndex_t ITZNumElems = ArrListGetCount(ElemList);
            if (TecUtilDataSetAddZone("IsoTopoZone", ITZNumNodes, ITZNumElems, 1,
                                      ZoneType_FEQuad, NULL))
            {
                LgIndex_t ne;
                NodeMap_pa NodeMapIso = NULL;
                NodeMap_pa NodeMapITZ = NULL;
                Set_pa ZonesAdded = TecUtilSetAlloc(TRUE);

                // new zone is always last zone
                TecUtilDataSetGetInfo(NULL, &IsoTopoZone, NULL);

                // Set the node map for the new sub-zone
                NodeMapIso = TecUtilDataNodeGetWritableRef(IsoSurfZone);
                NodeMapITZ = TecUtilDataNodeGetWritableRef(IsoTopoZone);
                for (ne = 1; ne <= ITZNumElems; ne++)
                {
                    LgIndex_t IsoElemNum;

                    // Get the element number in the full surface zone
                    Item = ArrListGetItem(ElemList, ne - 1);
                    IsoElemNum = Item.Long;

                    // Get the node numbers

                    // Convert isosurf node numbers to IsoTopo node numbers

                    // Set the IsoTopo node map for element

                }



                // fill new zone with values for all variables

                // inform Tecplot of new zone
                TecUtilSetAddMember(ZonesAdded, IsoTopoZone, TRUE);
                TecUtilStateChanged(StateChange_ZonesAdded,
                                    (ArbParam_t)ZonesAdded);
                TecUtilSetDealloc(&ZonesAdded);
            }
        }

        // Clean up
        ArrListDealloc(&ElemList);
        ArrListDealloc(&NodeList);
        ArrListDealloc(&NodeBoundaryList);


        // Done - dealloc the element map
        SurfElemMapDealloc(&IsoElemMap);
    }

    REQUIRE(IsoTopoIsValid(IsoTopo));
    ENSURE(VALID_BOOLEAN(IsOk));
    REQUIRE(IsoTopoZone == 0 || IsoTopoZone == NumZones + 1);
    return IsoTopoZone;
}

