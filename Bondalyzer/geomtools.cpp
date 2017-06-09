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
#include "GEOMTOOLS.h"
#include <cmath>




/**
 * Vector Cross product
 * param V1, V2h
 *     Vectors to be crossed
 * return
 *     Vector cross product
 */
inline void calcVectorCrossProduct(XYZ_s & out,
                                   XYZ_s const & V1,
                                   XYZ_s const & V2)
{
    out.X = V1.Y * V2.Z - V1.Z * V2.Y;
    out.Y = V1.Z * V2.X - V1.X * V2.Z;
    out.Z = V1.X * V2.Y - V1.Y * V2.X;
}

double DistanceSquaredXYZ(const XYZ_s A,
						  const XYZ_s B)
{
	double Result = -1;
	double Distance = (B.X - A.X) * (B.X - A.X)
					+ (B.Y - A.Y) * (B.Y - A.Y)
					+ (B.Z - A.Z) * (B.Z - A.Z);
	Result = Distance;
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
double PointTriangleDistance(XYZ_s const &  P,
                             XYZ_s const & Nd1,
                             XYZ_s const & Nd2,
                             XYZ_s const & Nd3,
                             XYZ_pa        PIntersect)
{
    double Dist = LARGEFLOAT;

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

    // Only continue if the nodes are unique (edges have non-zero lengths)
    if (E1Length > SMALLFLOAT && E2Length > SMALLFLOAT && E3Length > SMALLFLOAT)
    {
        // Compute the normal vector to the plane of the triangle
        calcVectorCrossProduct(Np, E1, E2);
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

        // Tentatively save the intersection point
        if (PIntersect != NULL)
        {
            PIntersect->X = P0.X;
            PIntersect->Y = P0.Y;
            PIntersect->Z = P0.Z;
        }

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
        calcVectorCrossProduct(E1CrossDP10, E1, DP10);
        if ((E1CrossDP10.X * Np.X + E1CrossDP10.Y * Np.Y + E1CrossDP10.Z * Np.Z) < 0.0)    // Outside
        {
            double DistAlongEdge;
            DistAlongEdge = (DP10.X * E1.X + DP10.Y * E1.Y + DP10.Z * E1.Z) / E1Length;

            IsOutside = TRUE;

            // Before beginning of edge: take distance between P and Nd1
            if (DistAlongEdge < 0.0)
            {
                DistLat = sqrt(DP10.X * DP10.X + DP10.Y * DP10.Y + DP10.Z * DP10.Z);

                // Modify the intersection point
                if (PIntersect != NULL)
                {
                    PIntersect->X = Nd1.X;
                    PIntersect->Y = Nd1.Y;
                    PIntersect->Z = Nd1.Z;
                }
            }

            // After end of edge: take distance between P and Nd2
            else if (DistAlongEdge > 1.0)
            {
                DistLat = sqrt(DP20.X * DP20.X + DP20.Y * DP20.Y + DP20.Z * DP20.Z);

                // Modify the intersection point
                if (PIntersect != NULL)
                {
                    PIntersect->X = Nd2.X;
                    PIntersect->Y = Nd2.Y;
                    PIntersect->Z = Nd2.Z;
                }
            }

            // Lies along edge - find normal distance from edge
            else
            {
                XYZ_s E1Norm;  // Normal to edge (orthogonal to surface normal)
                calcVectorCrossProduct(E1Norm, E1, Np);
                E1Norm.X /= E1Length;
                E1Norm.Y /= E1Length;
                E1Norm.Z /= E1Length;
                DistLat = DP10.X * E1Norm.X + DP10.Y * E1Norm.Y + DP10.Z * E1Norm.Z;

                // Modify the intersection point
                if (PIntersect != NULL)
                {
                    PIntersect->X = Nd1.X + DistAlongEdge * E1.X;
                    PIntersect->Y = Nd1.Y + DistAlongEdge * E1.Y;
                    PIntersect->Z = Nd1.Z + DistAlongEdge * E1.Z;
                }
            }
        }

        // Lateral distance, if outside edge 2
        calcVectorCrossProduct(E2CrossDP20, E2, DP20);
        if ((E2CrossDP20.X * Np.X + E2CrossDP20.Y * Np.Y + E2CrossDP20.Z * Np.Z) < 0.0)    // Outside
        {
            double DistAlongEdge;
            DistAlongEdge = (DP20.X * E2.X + DP20.Y * E2.Y + DP20.Z * E2.Z) / E2Length;

            IsOutside = TRUE;

            // Before beginning of edge: take distance between P and Nd2
            if (DistAlongEdge < 0.0)
            {
                DistLat = MIN(DistLat, sqrt(DP20.X * DP20.X + DP20.Y * DP20.Y + DP20.Z * DP20.Z));

                // Modify the intersection point
                if (PIntersect != NULL)
                {
                    PIntersect->X = Nd2.X;
                    PIntersect->Y = Nd2.Y;
                    PIntersect->Z = Nd2.Z;
                }
            }

            // After end of edge: take distance between P and Nd3
            else if (DistAlongEdge > 1.0)
            {
                DistLat = MIN(DistLat, sqrt(DP30.X * DP30.X + DP30.Y * DP30.Y + DP30.Z * DP30.Z));

                // Modify the intersection point
                if (PIntersect != NULL)
                {
                    PIntersect->X = Nd3.X;
                    PIntersect->Y = Nd3.Y;
                    PIntersect->Z = Nd3.Z;
                }
            }

            // Lies along edge - find normal distance from edge
            else
            {
                XYZ_s E2Norm;  // Normal to edge (orthogonal to surface normal)
                calcVectorCrossProduct(E2Norm, E2, Np);
                E2Norm.X /= E2Length;
                E2Norm.Y /= E2Length;
                E2Norm.Z /= E2Length;
                DistLat = DP20.X * E2Norm.X + DP20.Y * E2Norm.Y + DP20.Z * E2Norm.Z;

                // Modify the intersection point
                if (PIntersect != NULL)
                {
                    PIntersect->X = Nd2.X + DistAlongEdge * E2.X;
                    PIntersect->Y = Nd2.Y + DistAlongEdge * E2.Y;
                    PIntersect->Z = Nd2.Z + DistAlongEdge * E2.Z;
                }
            }
        }


        // Lateral distance, if outside edge 3
        calcVectorCrossProduct(E3CrossDP30, E3, DP30);
        if ((E3CrossDP30.X * Np.X + E3CrossDP30.Y * Np.Y + E3CrossDP30.Z * Np.Z) < 0.0)    // Outside
        {
            double DistAlongEdge;
            DistAlongEdge = (DP30.X * E3.X + DP30.Y * E3.Y + DP30.Z * E3.Z) / E3Length;

            IsOutside = TRUE;

            // Before beginning of edge: take distance between P and Nd3
            if (DistAlongEdge < 0.0)

            {
                DistLat = MIN(DistLat, sqrt(DP30.X * DP30.X + DP30.Y * DP30.Y + DP30.Z * DP30.Z));

                // Modify the intersection point
                if (PIntersect != NULL)
                {
                    PIntersect->X = Nd3.X;
                    PIntersect->Y = Nd3.Y;
                    PIntersect->Z = Nd3.Z;
                }
            }

            // After end of edge: take distance between P and Nd1
            else if (DistAlongEdge > 1.0)

            {
                DistLat = MIN(DistLat, sqrt(DP10.X * DP10.X + DP10.Y * DP10.Y + DP10.Z * DP10.Z));

                // Modify the intersection point
                if (PIntersect != NULL)
                {
                    PIntersect->X = Nd1.X;
                    PIntersect->Y = Nd1.Y;
                    PIntersect->Z = Nd1.Z;
                }
            }

            // Lies along edge - find normal distance from edge
            else
            {
                XYZ_s E3Norm;  // Normal to edge (orthogonal to surface normal)
                calcVectorCrossProduct(E3Norm, E3, Np);
                E3Norm.X /= E3Length;
                E3Norm.Y /= E3Length;
                E3Norm.Z /= E3Length;
                DistLat = DP30.X * E3Norm.X + DP30.Y * E3Norm.Y + DP30.Z * E3Norm.Z;

                // Modify the intersection point
                if (PIntersect != NULL)
                {
                    PIntersect->X = Nd3.X + DistAlongEdge * E3.X;
                    PIntersect->Y = Nd3.Y + DistAlongEdge * E3.Y;
                    PIntersect->Z = Nd3.Z + DistAlongEdge * E3.Z;
                }
            }
        }

        if (IsOutside)
            Dist = sqrt(DistNorm * DistNorm + DistLat * DistLat);
        else
            Dist = ABS(DistNorm);
    }

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
Boolean_t TestPointTriDist()
{
    Boolean_t Result = TRUE;
    XYZ_s     Point, Nd1, Nd2, Nd3;
    double    Dist, Diff;
    XYZ_s     Position;

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

    Dist = PointTriangleDistance(Point, Nd1, Nd2, Nd3, &Position);

    Diff = Dist - 2.0 * sqrt(3.0) / 3.0;
    if (ABS(Diff) > 1.0e-5) Result = FALSE;


    // Test Edge 1 (Nd1 - Nd2) distance calculation
    if (Result)
    {
        Point.X = 1.0;
        Point.Y = 1.0;
        Point.Z = 0.0;

        Dist = PointTriangleDistance(Point, Nd1, Nd2, Nd3, &Position);

        Diff = Dist - sqrt(2.0) / 2.0;
        if (ABS(Diff) > 1.0e-5) Result = FALSE;
    }


    // Test Edge 2 (Nd2 - Nd3) distance calculation
    if (Result)
    {
        Point.X = 0.0;
        Point.Y = 1.0;
        Point.Z = 1.0;

        Dist = PointTriangleDistance(Point, Nd1, Nd2, Nd3, &Position);

        Diff = Dist - sqrt(2.0) / 2.0;
        if (ABS(Diff) > 1.0e-5) Result = FALSE;
    }


    // Test Edge 3 (Nd3 - Nd1) distance calculation
    if (Result)
    {
        Point.X = 1.0;
        Point.Y = 0.0;
        Point.Z = 1.0;

        Dist = PointTriangleDistance(Point, Nd1, Nd2, Nd3, &Position);

        Diff = Dist - sqrt(2.0) / 2.0;
        if (ABS(Diff) > 1.0e-5) Result = FALSE;
    }



    // Test Node 1 distance calculation
    if (Result)
    {
        Point.X = 2.0;
        Point.Y = 0.0;
        Point.Z = 0.0;

        Dist = PointTriangleDistance(Point, Nd1, Nd2, Nd3, &Position);

        if (ABS(Dist - 1.0) > 1.0e-5) Result = FALSE;
    }


    // Test Node 2 distance calculation
    if (Result)
    {
        Point.X = 0.0;
        Point.Y = 2.0;
        Point.Z = 0.0;

        Dist = PointTriangleDistance(Point, Nd1, Nd2, Nd3, &Position);

        if (ABS(Dist - 1.0) > 1.0e-5) Result = FALSE;
    }



    // Test Node 3 distance calculation
    if (Result)
    {
        Point.X = 0.0;
        Point.Y = 0.0;
        Point.Z = 2.0;

        Dist = PointTriangleDistance(Point, Nd1, Nd2, Nd3, &Position);

        if (ABS(Dist - 1.0) > 1.0e-5) Result = FALSE;
    }

    ENSURE(VALID_BOOLEAN(Result));
    return Result;
}





/**
 * Find the element of the specified zone that is closest to the specified
 * point and return its element number.
 *
 * param SurfZone
 *     Number of the surface zone.
 * param X, Y, Z
 *     Location of the point.
 *
 * return (through pointer) XInter, YInter, ZInter
 *     Location of the closest point
 * return
 *     Element number of the closest element (if it worked), zero if it failed.
 */
LgIndex_t GeomToolsClosestElemNum(EntIndex_t SurfZone, // 1-based
                                  double     X,
                                  double     Y,
                                  double     Z,
                                  double    *XInter,
                                  double    *YInter,
                                  double    *ZInter)
{
    LgIndex_t  NumNodes, NumElems;
    LgIndex_t  ClosestElem = 0;

    /* Inform Tecplot that major data operation is beginning */
    TecUtilDataLoadBegin();

    /* Lock Tecplot */
    TecUtilLockStart(AddOnID);

    REQUIRE(1 <= SurfZone && SurfZone <= TecUtilDataSetGetNumZones());


    // Find the closest element to the point
    LgIndex_t    ne;
    EntIndex_t   XVarNum, YVarNum, ZVarNum;
    FieldData_pa XVarFDPtr, YVarFDPtr, ZVarFDPtr;
    double       MinDist = LARGEFLOAT;
    ZoneType_e   ZoneType;
    Boolean_t    IsQuad;
    XYZ_s        Point, Nd1, Nd2, Nd3, Nd4;
    XYZ_s        IntersectPoint;

    Point.X = X;
    Point.Y = Y;
    Point.Z = Z;

    // Find the number of Elements in the SurfZone
    TecUtilZoneGetInfo(SurfZone, &NumNodes, &NumElems,
                       NULL, NULL, NULL, NULL, NULL, NULL,
                       NULL, NULL, NULL, NULL, NULL);

    // FETriangle or FEQuad?
    ZoneType = TecUtilZoneGetType(SurfZone);
    if (ZoneType = ZoneType_FEQuad)
        IsQuad = TRUE;
    else if (ZoneType = ZoneType_FETriangle)
        IsQuad = FALSE;
    else
        CHECK(FALSE);

    // Get X, Y, Z variable numbers
    XVarNum = TecUtilVarGetNumByAssignment('X');
    YVarNum = TecUtilVarGetNumByAssignment('Y');
    ZVarNum = TecUtilVarGetNumByAssignment('Z');

    XVarFDPtr = TecUtilDataValueGetReadableRef(SurfZone, XVarNum);
    YVarFDPtr = TecUtilDataValueGetReadableRef(SurfZone, YVarNum);
    ZVarFDPtr = TecUtilDataValueGetReadableRef(SurfZone, ZVarNum);

    // Loop over elements of SurfZone to find the one closest to X,Y,Z
    for (ne = 1; ne < NumElems; ne++)
    {
        double Dist;
        int ii, ic;
        LgIndex_t NdNums[4], NdNumsTmp[4];
        LgIndex_t Nd1Num, Nd2Num, Nd3Num, Nd4Num;

        // Find unique node numbers
        for (ii = 0; ii < 4; ii++)
            NdNumsTmp[ii] = TecUtilDataNodeGetByZone(SurfZone, ne, ii + 1);

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

        Dist = PointTriangleDistance(Point, Nd1, Nd2, Nd3, &IntersectPoint);
        if (Dist < MinDist)
        {
            MinDist = Dist;
            ClosestElem = ne;
            *XInter = IntersectPoint.X;
            *YInter = IntersectPoint.Y;
            *ZInter = IntersectPoint.Z;
        }

        // If quad, do the second triangle (nodes 1, 3, and 4)
        Nd4Num = NdNums[3];
        if (IsQuad && Nd1Num != Nd4Num && Nd4Num != Nd3Num)
        {
            Nd4.X = TecUtilDataValueGetByRef(XVarFDPtr, Nd4Num);
            Nd4.Y = TecUtilDataValueGetByRef(YVarFDPtr, Nd4Num);
            Nd4.Z = TecUtilDataValueGetByRef(ZVarFDPtr, Nd4Num);
            Dist = PointTriangleDistance(Point, Nd1, Nd3, Nd4, &IntersectPoint);
            if (Dist < MinDist)
            {
                MinDist = Dist;
                ClosestElem = ne;
                *XInter = IntersectPoint.X;
                *YInter = IntersectPoint.Y;
                *ZInter = IntersectPoint.Z;
            }
        }
    }

    /* Finish lock for this function (may be nested) */
    TecUtilLockFinish(AddOnID);

    /* Inform Tecplot that major data operation is ending */
    TecUtilDataLoadEnd();

    REQUIRE(ClosestElem > 0 && ClosestElem <= NumElems);
    return ClosestElem;
}






/*
 * Find the XYZ intersection of a line segment with a triangle.
 *
 * param triCorner1, triCorner2, triCorner3
 *     XYZ locations of the corners of the triangle
 * param lineSegStart
 *     XYZ location of the start of the line segment
 * param lineSegEnd
 *     XYZ location of the end of the line segment
 *
 * return (through pointer) *naturalTUV
 *     Natural coordinates of intersection on triangle.  Pass NULL to not calculate
 *
 * return (through pointer) *intersectionPt
 *     Pointer to the XYZ values of the intersection point. Only valid if
 *     intersection point within the triangle (i.e. function returns TRUE).
 *
 * return
 *     TRUE if the intersection is within the bounds of the triangle and line segment, false otherwise.
 *
 */
Boolean_t GeomToolsLineSegTriangleIntersection(XYZ_s const & triCorner1,
                                               XYZ_s const & triCorner2,
                                               XYZ_s const & triCorner3,
                                               XYZ_s const & lineSegStart,
                                               XYZ_s const & lineSegEnd,
                                               XYZ_s *       naturalTUV,
                                               XYZ_s *       intersectionPt)
{
    Boolean_t inTriangle = FALSE;

    REQUIRE(VALID_REF_OR_NULL(naturalTUV));
    REQUIRE(VALID_REF(intersectionPt));

    // Edge 1: Vector from triangle corner 1 to triangle corner 2
    XYZ_s edge1;
    edge1.X = triCorner2.X - triCorner1.X;
    edge1.Y = triCorner2.Y - triCorner1.Y;
    edge1.Z = triCorner2.Z - triCorner1.Z;

    // Edge 2: Vector from triangle corner 1 to triangle corner 3
    XYZ_s edge2;
    edge2.X = triCorner3.X - triCorner1.X;
    edge2.Y = triCorner3.Y - triCorner1.Y;
    edge2.Z = triCorner3.Z - triCorner1.Z;

    // Vector from the line segment start to line segment end
    XYZ_s lineSegDir;
    lineSegDir.X = lineSegEnd.X - lineSegStart.X;
    lineSegDir.Y = lineSegEnd.Y - lineSegStart.Y;
    lineSegDir.Z = lineSegEnd.Z - lineSegStart.Z;

    // Save intermediate results
    XYZ_s lineDirXEdge2;
    calcVectorCrossProduct(lineDirXEdge2, lineSegDir, edge2);

    // Compute the determinant
    double det = edge1.X * lineDirXEdge2.X + edge1.Y * lineDirXEdge2.Y + edge1.Z * lineDirXEdge2.Z;

    if (det <= -EPSILON_GT || det >= EPSILON_GT)
    {
        // Vector from triangle corner 1 to line segment start
        XYZ_s corner1ToStart;
        corner1ToStart.X = lineSegStart.X - triCorner1.X;
        corner1ToStart.Y = lineSegStart.Y - triCorner1.Y;
        corner1ToStart.Z = lineSegStart.Z - triCorner1.Z;

        // Compute the first triangle natural coordinate u.
        double invDet = 1.0 / det;
        double uNat = (lineDirXEdge2.X * corner1ToStart.X +
                       lineDirXEdge2.Y * corner1ToStart.Y +
                       lineDirXEdge2.Z * corner1ToStart.Z) * invDet;

        //if (-0.0001 <= uNat && uNat <= 1.0001)
        //if (-0.001 <= uNat && uNat <= 1.001)
        double const natCoordEpsilon = 0.001;
        if (-natCoordEpsilon <= uNat && uNat <= 1.0+natCoordEpsilon)
        {
            // Compute the second triangle natural coordinate v.
            XYZ_s corner1ToStartXEdge1;
            calcVectorCrossProduct(corner1ToStartXEdge1, corner1ToStart, edge1);

            double vNat = (corner1ToStartXEdge1.X * lineSegDir.X +
                           corner1ToStartXEdge1.Y * lineSegDir.Y +
                           corner1ToStartXEdge1.Z * lineSegDir.Z) * invDet;

            //if (-0.0001 <= vNat && uNat+vNat <= 1.0001)
            //if (-0.001 <= vNat && uNat+vNat <= 1.001)
            if (-natCoordEpsilon <= vNat && uNat+vNat <= 1.0+natCoordEpsilon)
            {
                // Calculate XYZ of intersection, and check that it's on the line segment, not just the line
                double oneMinusUMinusV = 1.0 - uNat - vNat;
                double xx = oneMinusUMinusV*triCorner1.X + uNat*triCorner2.X + vNat*triCorner3.X;
                double const epsilon = 0.0;
                if ( lineSegDir.X == 0.0 ||
                     ( lineSegStart.X-epsilon <= xx && xx <= lineSegEnd.X+epsilon ) ||
                     ( lineSegEnd.X-epsilon <= xx && xx <= lineSegStart.X+epsilon ) )
                {
                    double yy = oneMinusUMinusV*triCorner1.Y + uNat*triCorner2.Y + vNat*triCorner3.Y;
                    if ( lineSegDir.Y == 0.0 ||
                         ( lineSegStart.Y-epsilon <= yy && yy <= lineSegEnd.Y+epsilon ) ||
                         ( lineSegEnd.Y-epsilon <= yy && yy <= lineSegStart.Y+epsilon ) )
                    {
                        double zz = oneMinusUMinusV*triCorner1.Z + uNat*triCorner2.Z + vNat*triCorner3.Z;
                        if ( lineSegDir.Z == 0.0 ||
                             ( lineSegStart.Z-epsilon <= zz && zz <= lineSegEnd.Z+epsilon ) ||
                             ( lineSegEnd.Z-epsilon <= zz && zz <= lineSegStart.Z+epsilon ) )
                        {
                            inTriangle = TRUE;
                            intersectionPt->X = xx;
                            intersectionPt->Y = yy;
                            intersectionPt->Z = zz;
                            if ( naturalTUV )
                            {
                                naturalTUV->X = (corner1ToStartXEdge1.X * edge2.X +
                                                 corner1ToStartXEdge1.Y * edge2.Y +
                                                 corner1ToStartXEdge1.Z * edge2.Z) * invDet;
                                naturalTUV->Y = uNat;
                                naturalTUV->Z = vNat;
                            }
                        }
                    }
                }
            }
        }
    }

    ENSURE(VALID_BOOLEAN(inTriangle));
    return inTriangle;
}


/*
 * Find the XYZ intersection(s) of a line segment with a quad.  Note there may be two
 * intersections if the quad is non-planar.
 *
 * param corner1, corner2, corner3, corner4
 *     XYZ locations of the corners of the quad (going around the quad)
 * param lineSegStart
 *     XYZ location of the start of the line segment
 * param lineSegEnd
 *     XYZ location of the end of the line segment
 *
 * return (through pointer) *numIntersections
 *     Pointer to the number of intersections found.
 *
 * return (through pointer) *intersectionPts
 *     Pointer to the XYZ values of the intersection points. Only valid if
 *     intersection point within the quad (i.e. function returns TRUE).
 *
 * return
 *     TRUE if the intersection is within the bounds of the quad and line segment, FALSE otherwise.
 *
 */
Boolean_t GeomToolsLineSegQuadIntersections(XYZ_s const & corner1,
                                            XYZ_s const & corner2,
                                            XYZ_s const & corner3,
                                            XYZ_s const & corner4,
                                            XYZ_s const & lineSegStart,
                                            XYZ_s const & lineSegEnd,
                                            LgIndex_t *   numIntersections,
                                            XYZ_s *       intersectionPts)
{
    Boolean_t result;
    REQUIRE(VALID_REF(numIntersections));
    REQUIRE(VALID_REF(intersectionPts));

    LgIndex_t intersectionCount = 0;
    if ( GeomToolsLineSegTriangleIntersection(corner1, corner2, corner3,
                                              lineSegStart, lineSegEnd,
                                              NULL,
                                              intersectionPts) )
    {
        intersectionCount++;
    }
    if ( GeomToolsLineSegTriangleIntersection(corner1, corner3, corner4,
                                              lineSegStart, lineSegEnd,
                                              NULL,
                                              &intersectionPts[intersectionCount]) )
    {
        intersectionCount++;
    }

    result = (intersectionCount>0);
    *numIntersections = intersectionCount;

    ENSURE(VALID_BOOLEAN(result));
    return result;
}








/*
 * Call the callback function for each XYZ intersection points of a line segment with a surface
 * (Tecplot IJ-ordered zone).
 *
 * param surfaceZoneNum
 *     Tecplot zone number of the surface.  Must be IJ-ordered
 * param lineSegStartPt
 *     XYZ location of start of line segment.
 * param lineSegEndPt
 *     XYZ location of end of line segment
 * param intersectionCallback
 *     callback function
 * param clientData
 *     callback client data
 */
void GeomToolsLineSegSurfaceIntersection(EntIndex_t const        surfaceZoneNum, // 1-based
                                         XYZ_s const &           lineSegStartPt,
                                         XYZ_s const &           lineSegEndPt,
                                         IntersectionCallback_pf intersectionCallback,
                                         void *                  clientData)

{
    REQUIRE(VALID_FN_REF(intersectionCallback));

    /* Inform Tecplot that major data operation is beginning */
    TecUtilDataLoadBegin();

    /* Lock Tecplot */
    TecUtilLockStart(AddOnID);

    /* Check surfaceZoneNum */
    REQUIRE(1 <= surfaceZoneNum && surfaceZoneNum <= TecUtilDataSetGetNumZones());

    // Limited to ordered surface zones for now
    REQUIRE(TecUtilZoneIsOrdered(surfaceZoneNum));

    LgIndex_t iMax, jMax, kMax;
    // Use NULL for values we're not interested in
    TecUtilZoneGetInfo(surfaceZoneNum, &iMax, &jMax, &kMax, NULL, NULL, NULL, NULL, NULL,
                       NULL, NULL, NULL, NULL, NULL);

    // Must be IJ-ordered surface
    REQUIRE(iMax > 2 && jMax > 2 && kMax == 1);

    XYZ_s lineSegMin;
    XYZ_s lineSegMax;
    XYZ_s zoneMin;
    XYZ_s zoneMax;

    // quickly test to see if there is any chance of this working
    lineSegMin.X = MIN(lineSegStartPt.X, lineSegEndPt.X);
    lineSegMax.X = MAX(lineSegStartPt.X, lineSegEndPt.X);
    EntIndex_t xVarNum = TecUtilVarGetNumByAssignment('X');
    FieldData_pa xFDPtr = TecUtilDataValueGetReadableRef(surfaceZoneNum, xVarNum);
    TecUtilDataValueGetMinMaxByRef(xFDPtr, &zoneMin.X, &zoneMax.X);
    if ( lineSegMin.X <= zoneMax.X && lineSegMax.X >= zoneMin.X )
    {
        lineSegMin.Y = MIN(lineSegStartPt.Y, lineSegEndPt.Y);
        lineSegMax.Y = MAX(lineSegStartPt.Y, lineSegEndPt.Y);
        EntIndex_t yVarNum = TecUtilVarGetNumByAssignment('Y');
        FieldData_pa yFDPtr = TecUtilDataValueGetReadableRef(surfaceZoneNum, yVarNum);
        TecUtilDataValueGetMinMaxByRef(yFDPtr, &zoneMin.Y, &zoneMax.Y);
        if ( lineSegMin.Y <= zoneMax.Y && lineSegMax.Y >= zoneMin.Y )
        {
            lineSegMin.Z = MIN(lineSegStartPt.Z, lineSegEndPt.Z);
            lineSegMax.Z = MAX(lineSegStartPt.Z, lineSegEndPt.Z);
            EntIndex_t zVarNum = TecUtilVarGetNumByAssignment('Z');
            FieldData_pa zFDPtr = TecUtilDataValueGetReadableRef(surfaceZoneNum, zVarNum);
            TecUtilDataValueGetMinMaxByRef(zFDPtr, &zoneMin.Z, &zoneMax.Z);
            if ( lineSegMin.Z <= zoneMax.Z && lineSegMax.Z >= zoneMin.Z )
            {
                // Loop over cells of the zone, checking each cell for intersection
                // For now, break each cell into two triangles.
                LgIndex_t iCellMax = iMax-1;
                LgIndex_t jCellMax = jMax-1;
                for (LgIndex_t jj = 0; jj < jCellMax; jj++)
                {
                    for (LgIndex_t ii = 0; ii < iCellMax; ii++)
                    {
                        // Get quad at (i,j):(i+1,j):(i+1,j+1):(i,j+1)
                        LgIndex_t pt1Index = ii + jj*iMax;
                        XYZ_s corner1;
                        corner1.X = TecUtilDataValueGetByRef(xFDPtr, pt1Index+1); // TecUtil is 1-based
                        corner1.Y = TecUtilDataValueGetByRef(yFDPtr, pt1Index+1);
                        corner1.Z = TecUtilDataValueGetByRef(zFDPtr, pt1Index+1);

                        LgIndex_t pt2Index = pt1Index+1;
                        XYZ_s corner2;
                        corner2.X = TecUtilDataValueGetByRef(xFDPtr, pt2Index+1); // TecUtil is 1-based
                        corner2.Y = TecUtilDataValueGetByRef(yFDPtr, pt2Index+1);
                        corner2.Z = TecUtilDataValueGetByRef(zFDPtr, pt2Index+1);

                        XYZ_s corner3;
                        LgIndex_t pt3Index = pt1Index+iMax+1;
                        corner3.X = TecUtilDataValueGetByRef(xFDPtr, pt3Index+1); // TecUtil is 1-based
                        corner3.Y = TecUtilDataValueGetByRef(yFDPtr, pt3Index+1);
                        corner3.Z = TecUtilDataValueGetByRef(zFDPtr, pt3Index+1);

                        XYZ_s corner4;
                        LgIndex_t pt4Index = pt1Index+iMax;
                        corner4.X = TecUtilDataValueGetByRef(xFDPtr, pt4Index+1); // TecUtil is 1-based
                        corner4.Y = TecUtilDataValueGetByRef(yFDPtr, pt4Index+1); // TecUtil is 1-based
                        corner4.Z = TecUtilDataValueGetByRef(zFDPtr, pt4Index+1); // TecUtil is 1-based

                        // Check Quad (i,j):(i+1,j):(i+1,j+1):(i,j+1) 
                        int numIntersections;
                        XYZ_s intersectionPts[2];
                        if ( GeomToolsLineSegQuadIntersections(corner1, corner2, corner3, corner4,
                                                               lineSegStartPt, lineSegEndPt,
                                                               &numIntersections, intersectionPts))
                        {
                            for ( int pp = 0; pp < numIntersections; pp++ )
                                intersectionCallback(intersectionPts[pp], clientData);
                        }
                    }
                }
            }
        }
    }

    /* Finish lock for this function (may be nested) */
    TecUtilLockFinish(AddOnID);

    /* Inform Tecplot that major data operation is ending */
    TecUtilDataLoadEnd();
}




/*
 * Test the GeomTools functions. Currently only tests GeomToolsLineSegTriangleIntersection.
 *
 * return
 *     TRUE if the functions pass the test, FALSE otherwise.
 *
 */
Boolean_t GeomToolsTest()
{
    // Test GeomToolsLineSegTriangleIntersection
    XYZ_s triCorner1;
    triCorner1.X = 0.0;
    triCorner1.Y = 0.0;
    triCorner1.Z = 0.0;

    XYZ_s triCorner2;
    triCorner2.X = 1.0;
    triCorner2.Y = 0.0;
    triCorner2.Z = 0.0;
    
    XYZ_s triCorner3;
    triCorner3.X = 1.0;
    triCorner3.Y = 1.0;
    triCorner3.Z = 0.0;

    XYZ_s lineSegStart;
    lineSegStart.X =  0.25;
    lineSegStart.Y =  0.15;
    lineSegStart.Z = -0.5;

    XYZ_s lineSegEnd;
    lineSegEnd.X = 0.25;
    lineSegEnd.Y = 0.15;
    lineSegEnd.Z = 0.5;

    XYZ_s naturalTUV;
    naturalTUV.X = -99.0; // bogus values to see if they are set correct
    naturalTUV.Y = -99.0;
    naturalTUV.Z = -99.0;

    XYZ_s intersection;
    intersection.X = -99.0; // bogus values to see if they are set correct
    intersection.Y = -99.0;
    intersection.Z = -99.0;

    CHECK(GeomToolsLineSegTriangleIntersection(triCorner1, triCorner2, triCorner3,
                                               lineSegStart, lineSegEnd,
                                               &naturalTUV, &intersection));
    CHECK(0.49999 < naturalTUV.X && naturalTUV.X < 0.50001);
    CHECK(0.24999 < intersection.X && intersection.X < 0.25001);
    CHECK(0.14999 < intersection.Y && intersection.Y < 0.15001);
    CHECK(-0.00001 < intersection.Z && intersection.Z < 0.00001);

    lineSegEnd.Z = -1.5;
    CHECK(!GeomToolsLineSegTriangleIntersection(triCorner1, triCorner2, triCorner3,
                                                lineSegStart, lineSegEnd,
                                                &naturalTUV, &intersection));

    // Rotate 30deg about z-axis and 45deg about x-axis, and translate 0.1, 0.2, 0.3,
    // and try again
    // TODO:

    return TRUE; // will assert if problem
}


/*
	Takes the radius of a circle or sphere and the straight-
	line distance between two points on its surface and returns
	the length of the line on the surface that connects the two points
*/
double	SurfaceLength(const double	Radius,
						const double	PointDistance)
{
	REQUIRE(Radius > 0 && PointDistance > 0);
	
	double Length = 2.0 * Radius * asin(PointDistance / 2.0 / Radius);

	return Length;
}