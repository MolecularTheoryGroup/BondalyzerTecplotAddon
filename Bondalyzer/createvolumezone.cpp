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
#include "GEOMTOOLS.h"
#include "SURFELEMMAP.h"
#include "NORMALS.h"
#include "ZONEVARINFO.h"
#include "CRITPOINTS.h"
#include "GRADPATH.h"
#include "BUNDLES.h"
#include "BONDBUNDLE.h"
#include "CREATEVOLUMEZONE.h"
#include <vector>

typedef signed char BundleInt_t;
static FieldDataType_e bundleIntFieldDataType = FieldDataType_Byte;
static BundleInt_t const noBondBundle = (BundleInt_t)0;
static LgIndex_t const noPt = -1;

static void checkForIntersection(XYZ_s& intersectionPt,
                                 void*  clientData)
{
    Boolean_t *hasIntersection = (Boolean_t*)clientData;
    REQUIRE(VALID_REF(hasIntersection));
    REQUIRE(VALID_BOOLEAN(*hasIntersection));
    *(Boolean_t*)clientData = TRUE;
    ENSURE(VALID_BOOLEAN(*hasIntersection));
}


static Boolean_t checkBondDirection(XYZ_s&         pt,
                                    BundleInt_t    ptBundleNum,
                                    LgIndex_t      otherPtIndex,
                                    FieldData_pa   XD,
                                    FieldData_pa   YD,
                                    FieldData_pa   ZD,
                                    BondBundle_pa  bondBundle,
                                    BundleInt_t*   whichBondBundle)
{
    REQUIRE(ptBundleNum!=noBondBundle);
    REQUIRE(bondBundle->BondNum==ptBundleNum);

    // This assertion keeps firing indicating a given problem, but we can
    // still process other bond bundles, so we make it only show once
    static Boolean_t keepChecking = TRUE;
    if ( keepChecking )
    {
        REQUIRE(whichBondBundle[otherPtIndex]==noBondBundle);
        keepChecking = (whichBondBundle[otherPtIndex]==noBondBundle);
    }

    XYZ_s otherPt;
    otherPt.X = TecUtilDataValueGetByRef(XD, otherPtIndex+1); // tecutil is 1-based
    otherPt.Y = TecUtilDataValueGetByRef(YD, otherPtIndex+1);
    otherPt.Z = TecUtilDataValueGetByRef(ZD, otherPtIndex+1);

    CHECK(pt.X!=otherPt.X || pt.Y!=otherPt.Y || pt.Z!=otherPt.Z); //different points

    Boolean_t hasIntersection = FALSE;
    BondBundleProcessSurfaceLineSegIntersections(bondBundle, pt, otherPt,
                                                 checkForIntersection, (void*)&hasIntersection);

    Boolean_t inSameBundle = !hasIntersection; // if we intersected, we aren't in the same bundle anymore

    ENSURE(VALID_BOOLEAN(inSameBundle));
    return inSameBundle;
}


inline void makeBrick(NodeMap_pa bondZoneNM, // just count if NULL
                      LgIndex_t& bondZoneElem,
                      LgIndex_t  n1,
                      LgIndex_t  n2,
                      LgIndex_t  n3,
                      LgIndex_t  n4,
                      LgIndex_t  n5,
                      LgIndex_t  n6,
                      LgIndex_t  n7,
                      LgIndex_t  n8)
{
    REQUIRE(VALID_REF_OR_NULL(bondZoneNM));
    if ( bondZoneNM != NULL )
    {
        TecUtilDataNodeSetByRef(bondZoneNM, bondZoneElem+1, 1, n1+1); // TecUtil is 1-based
        TecUtilDataNodeSetByRef(bondZoneNM, bondZoneElem+1, 2, n2+1);
        TecUtilDataNodeSetByRef(bondZoneNM, bondZoneElem+1, 3, n3+1);
        TecUtilDataNodeSetByRef(bondZoneNM, bondZoneElem+1, 4, n4+1);
        TecUtilDataNodeSetByRef(bondZoneNM, bondZoneElem+1, 5, n5+1);
        TecUtilDataNodeSetByRef(bondZoneNM, bondZoneElem+1, 6, n6+1);
        TecUtilDataNodeSetByRef(bondZoneNM, bondZoneElem+1, 7, n7+1);
        TecUtilDataNodeSetByRef(bondZoneNM, bondZoneElem+1, 8, n8+1);
    }
    bondZoneElem++;
}


// make a triangular pyramid, aka a tetrahedron
inline void makeTriPyramid(NodeMap_pa bondZoneNM, // just count if NULL
                           LgIndex_t& bondZoneElem,
                           LgIndex_t  n1, // n1-3 are the base
                           LgIndex_t  n2,
                           LgIndex_t  n3,
                           LgIndex_t  top)
{
    REQUIRE(VALID_REF_OR_NULL(bondZoneNM));
    makeBrick(bondZoneNM, bondZoneElem, n1, n2, n3, n3, top, top, top, top);
}

// make a quadrilateral pyramid
inline void makeQuadPyramid(NodeMap_pa bondZoneNM, // just count if NULL
                            LgIndex_t& bondZoneElem,
                            LgIndex_t  n1, // n1-n4 are the base
                            LgIndex_t  n2,
                            LgIndex_t  n3,
                            LgIndex_t  n4,
                            LgIndex_t  top)
{
    REQUIRE(VALID_REF_OR_NULL(bondZoneNM));
    makeBrick(bondZoneNM, bondZoneElem, n1, n2, n3, n4, top, top, top, top);
}

// make a pentagonal pyramid using multiple pyramids
inline void makePentPyramid(NodeMap_pa bondZoneNM, // just count if NULL
                            LgIndex_t& bondZoneElem,
                            LgIndex_t  n1, // n1-n5 are the base
                            LgIndex_t  n2,
                            LgIndex_t  n3,
                            LgIndex_t  n4,
                            LgIndex_t  n5,
                            LgIndex_t  top)
{
    REQUIRE(VALID_REF_OR_NULL(bondZoneNM));
    makeTriPyramid(bondZoneNM, bondZoneElem, n2, n3, n4, top);
    makeQuadPyramid(bondZoneNM, bondZoneElem, n1, n2, n4, n5, top);
}


static void createFacePrisms(LgIndex_t    n1,
                             LgIndex_t    n2,
                             LgIndex_t    n3,
                             LgIndex_t    n4,
                             LgIndex_t    e1_2, // point on edge between n1 & n2
                             LgIndex_t    e2_3,
                             LgIndex_t    e3_4,
                             LgIndex_t    e4_1,
                             LgIndex_t    ff, // ff node
                             LgIndex_t    cc, // cell center
                             NodeMap_pa   bondZoneNM, // just count if NULL
                             LgIndex_t&   bondZoneElem)
{
//     REQUIRE(IMPLICATION((e1_2!=noPt),(n1==noPt)^(n2==noPt)));
//     REQUIRE(IMPLICATION((e2_3!=noPt),(n2==noPt)^(n3==noPt)));
//     REQUIRE(IMPLICATION((e3_4!=noPt),(n3==noPt)^(n4==noPt)));
//     REQUIRE(IMPLICATION((e4_1!=noPt),(n4==noPt)^(n1==noPt)));
    //REQUIRE(cc!=noPt);
    REQUIRE(VALID_REF_OR_NULL(bondZoneNM)); // just count if NULL

    if ( cc == noPt )
    {
        // do nothing
    }
    else if ( n1 == noPt && n2 == noPt && n3 == noPt && n4 == noPt )
    {
        //CHECK(e1_2==noPt && e2_3==noPt && e3_4==noPt && e4_1==noPt );
        // do nothing
    }
    else if ( n1 != noPt && n2 != noPt && n3 != noPt && n4 != noPt )
    {
        makeQuadPyramid(bondZoneNM, bondZoneElem, n1, n2, n3, n4, cc);
    }
    else
    {
        if ( ff == noPt )
        {
            // handle single corners
            if ( e4_1 != noPt && n1 != noPt && e1_2 != noPt )
                makeTriPyramid(bondZoneNM, bondZoneElem, e4_1, n1, e1_2, cc);
            if ( e1_2 != noPt && n2 != noPt && e2_3 != noPt )
                makeTriPyramid(bondZoneNM, bondZoneElem, e1_2, n2, e2_3, cc);
            if ( e2_3 != noPt && n3 != noPt && e3_4 != noPt )
                makeTriPyramid(bondZoneNM, bondZoneElem, e2_3, n3, e3_4, cc);
            if ( e3_4 != noPt && n4 != noPt && e4_1 != noPt )
                makeTriPyramid(bondZoneNM, bondZoneElem, e3_4, n4, e4_1, cc);

            // handle double corners
            if ( e4_1 != noPt && n1 != noPt && n2 != noPt && e2_3 != noPt )
                makeQuadPyramid(bondZoneNM, bondZoneElem, e4_1, n1, n2, e2_3, cc);
            if ( e1_2 != noPt && n2 != noPt && n3 != noPt && e3_4 != noPt )
                makeQuadPyramid(bondZoneNM, bondZoneElem, e1_2, n2, n3, e3_4, cc);
            if ( e2_3 != noPt && n3 != noPt && n4 != noPt && e4_1 != noPt )
                makeQuadPyramid(bondZoneNM, bondZoneElem, e2_3, n3, n4, e4_1, cc);
            if ( e3_4 != noPt && n4 != noPt && n1 != noPt && e1_2 != noPt )
                makeQuadPyramid(bondZoneNM, bondZoneElem, e3_4, n4, n1, e1_2, cc);

            // handle triple corners
            if ( e4_1 != noPt && n1 != noPt && n2 != noPt && n3 != noPt && e3_4 != noPt )
                makePentPyramid(bondZoneNM, bondZoneElem, e4_1, n1, n2, n3, e3_4, cc);
            if ( e1_2 != noPt && n2 != noPt && n3 != noPt && n4 != noPt && e4_1 != noPt )
                makePentPyramid(bondZoneNM, bondZoneElem, e1_2, n2, n3, n4, e4_1, cc);
            if ( e2_3 != noPt && n3 != noPt && n4 != noPt && n1 != noPt && e1_2 != noPt )
                makePentPyramid(bondZoneNM, bondZoneElem, e2_3, n3, n4, n1, e1_2, cc);
            if ( e3_4 != noPt && n4 != noPt && n1 != noPt && n2 != noPt && e2_3 != noPt )
                makePentPyramid(bondZoneNM, bondZoneElem, e3_4, n4, n1, n2, e2_3, cc);
        }
        else // has a face node
        {
            LgIndex_t pts[] = { n1, e1_2, n2, e2_3, n3, e3_4, n4, e4_1 };
            LgIndex_t nPts = sizeof(pts)/sizeof(pts[0]);
            CHECK(nPts==8);
            // if two consecutive points are in, make a tet including face and cell center
            int prevIndex = nPts-1;
            for ( int curIndex = 0; curIndex < nPts; curIndex++ )
            {
                if ( pts[prevIndex] != noPt && pts[curIndex] != noPt )
                    makeTriPyramid(bondZoneNM, bondZoneElem, pts[prevIndex], pts[curIndex], ff, cc);
                prevIndex = curIndex;
            }
            // if two consecutive nodes are in, make a tet out of that edge including face and cell center
            // every other element of pts is a node starting at 0
            prevIndex = nPts-2;
            for ( int curIndex = 0; curIndex < nPts; curIndex += 2 )
            {
                if ( pts[prevIndex] != noPt && pts[curIndex] != noPt )
                    makeTriPyramid(bondZoneNM, bondZoneElem, pts[prevIndex], pts[curIndex], ff, cc);
                prevIndex = curIndex;
            }
        }
    }
}

inline void createElements(LgIndex_t    n1, // point at 1st node
                           LgIndex_t    n2,
                           LgIndex_t    n3,
                           LgIndex_t    n4,
                           LgIndex_t    n5,
                           LgIndex_t    n6,
                           LgIndex_t    n7,
                           LgIndex_t    n8,
                           LgIndex_t    e1_2, // point on edge between n1 & n2
                           LgIndex_t    e2_3,
                           LgIndex_t    e3_4,
                           LgIndex_t    e4_1,
                           LgIndex_t    e1_5,
                           LgIndex_t    e2_6,
                           LgIndex_t    e3_7,
                           LgIndex_t    e4_8,
                           LgIndex_t    e5_6,
                           LgIndex_t    e6_7,
                           LgIndex_t    e7_8,
                           LgIndex_t    e8_5,
                           LgIndex_t    f1_2_3_4, // point on ij-ff @ k-min  (ff with nodes 1,2,3,4)
                           LgIndex_t    f5_8_7_6, // point on ij-ff @ k-max
                           LgIndex_t    f1_4_8_5, // point on jk-ff @ i-min
                           LgIndex_t    f2_6_7_3, // point on jk-ff @ i-max
                           LgIndex_t    f1_5_6_2, // point on ki-ff @ j-min
                           LgIndex_t    f3_7_8_4, // point on ki-ff @ j-max
                           LgIndex_t    cc, // center cell
                           NodeMap_pa   bondZoneNM, // just count if NULL
                           LgIndex_t&   bondZoneElem)
{
//     REQUIRE(IMPLICATION((e1_2!=noPt),(n1==noPt)^(n2==noPt)));
//     REQUIRE(IMPLICATION((e2_3!=noPt),(n2==noPt)^(n3==noPt)));
//     REQUIRE(IMPLICATION((e3_4!=noPt),(n3==noPt)^(n4==noPt)));
//     REQUIRE(IMPLICATION((e4_1!=noPt),(n4==noPt)^(n1==noPt)));
//     REQUIRE(IMPLICATION((e1_5!=noPt),(n1==noPt)^(n5==noPt)));
//     REQUIRE(IMPLICATION((e2_6!=noPt),(n2==noPt)^(n6==noPt)));
//     REQUIRE(IMPLICATION((e3_7!=noPt),(n3==noPt)^(n7==noPt)));
//     REQUIRE(IMPLICATION((e4_8!=noPt),(n4==noPt)^(n8==noPt)));
//     REQUIRE(IMPLICATION((e5_6!=noPt),(n5==noPt)^(n6==noPt)));
//     REQUIRE(IMPLICATION((e6_7!=noPt),(n6==noPt)^(n7==noPt)));
//     REQUIRE(IMPLICATION((e7_8!=noPt),(n7==noPt)^(n8==noPt)));
//     REQUIRE(IMPLICATION((e8_5!=noPt),(n8==noPt)^(n5==noPt)));
    // if there is a cell-center, one of the nodes should be outside the bundle
    // unfortunately not true either
    //REQUIRE(IMPLICATION((cc!=noPt),(n1==noPt)||(n2==noPt)||(n3==noPt)||(n4==noPt)||(n5==noPt)||(n6==noPt)||(n7==noPt)||(n8==noPt)));
    // if there is a cell-center, one of the nodes should be inside the bundle
    // unfortunately, not always case
    //REQUIRE(IMPLICATION((cc!=noPt),(n1!=noPt)||(n2!=noPt)||(n3!=noPt)||(n4!=noPt)||(n5!=noPt)||(n6!=noPt)||(n7!=noPt)||(n8!=noPt)));

    int count = 0;
    if ( n1 >= 0 ) count++;
    if ( n2 >= 0 ) count++;
    if ( n3 >= 0 ) count++;
    if ( n4 >= 0 ) count++;
    if ( n5 >= 0 ) count++;
    if ( n6 >= 0 ) count++;
    if ( n7 >= 0 ) count++;
    if ( n8 >= 0 ) count++;

    if ( count == 0 )
    {
        // easiest case: no element (although theoretically there could be one going through faces without
        // going through nodes, but we skip that possibility for now)
    }
    else if ( count == 8)
    {
        // next easiest case: complete element (again, there could be a part taken out of a ff without
        // taking out any nodes, but we are skipping that case for now)
        makeBrick(bondZoneNM, bondZoneElem, n1, n2, n3, n4, n5, n6, n7, n8);
    }
    else
    {
        // partial cell: idea here it to make prisms with the cell center and the faces
        // make sure to get orientation correct so integration works (no negative volumes)
        createFacePrisms(n1, n2, n3, n4, e1_2, e2_3, e3_4, e4_1, f1_2_3_4, cc, bondZoneNM, bondZoneElem);
        createFacePrisms(n1, n4, n8, n5, e4_1, e4_8, e8_5, e1_5, f1_4_8_5, cc, bondZoneNM, bondZoneElem);
        createFacePrisms(n1, n5, n6, n2, e1_5, e5_6, e2_6, e1_2, f1_5_6_2, cc, bondZoneNM, bondZoneElem);
        createFacePrisms(n2, n6, n7, n3, e2_6, e6_7, e3_7, e2_3, f2_6_7_3, cc, bondZoneNM, bondZoneElem);
        createFacePrisms(n3, n7, n8, n4, e3_7, e7_8, e4_8, e3_4, f3_7_8_4, cc, bondZoneNM, bondZoneElem);
        createFacePrisms(n5, n8, n7, n6, e8_5, e7_8, e6_7, e5_6, f5_8_7_6, cc, bondZoneNM, bondZoneElem);
    }
}


static float const noIntersection = -999.0F;


// add 9 and subtract 10 instead of just subtracting 1 so we can get -.01 to be -1 without using floor()
#define xValueToIIndex(xVal) ( LgIndex_t( double((xVal)) + 9.0 + ijkEpsilon ) - 10 )
#define yValueToJIndex(yVal) ( LgIndex_t( double((yVal)) + 9.0 + ijkEpsilon ) - 10 )
#define zValueToKIndex(zVal) ( LgIndex_t( double((zVal)) + 9.0 + ijkEpsilon ) - 10 )

#define iIndexToXValue(ii) ( double((ii)) + 1.0 )
#define jIndexToYValue(jj) ( double((jj)) + 1.0 )
#define kIndexToZValue(kk) ( double((kk)) + 1.0 )

typedef struct
{
    float rr;
    float ss;
    float tt;
} FloatRST_s;

typedef struct
{
    LgIndex_t  baseZoneIMax;
    LgIndex_t  baseZoneJMax;
    LgIndex_t  baseZoneKMax;
    LgIndex_t  baseZoneIJMax;
    LgIndex_t  baseZoneNumPts;
    float*     iEdgeIntersection;
    float*     jEdgeIntersection;
    float*     kEdgeIntersection;
    LgIndex_t* ijFacePtOffset;
    LgIndex_t* jkFacePtOffset;
    LgIndex_t* kiFacePtOffset;
    LgIndex_t* centerPtOffset;
    FloatRST_s*rstValues;
    LgIndex_t* rstNumValues;
    LgIndex_t  rstMaxValues;
} IntersectionClientData_s;

void recordIEdgeIntersection(XYZ_s& intersectionPt,
                             void*  clientData)
{
    IntersectionClientData_s* intersectionClientData = (IntersectionClientData_s *)clientData;
    REQUIRE(VALID_REF(intersectionClientData));

    LgIndex_t const ii = xValueToIIndex(intersectionPt.X);
    // an i-edge is at a certain jk location and goes from i to i+1
    // that is, an i-edge is (i,j,k) to (i+1,j,k) thus i must be between 0 and imax-1
    if ( 0<=ii && ii<intersectionClientData->baseZoneIMax-1 )
    {
        LgIndex_t const jj = yValueToJIndex(intersectionPt.Y);
        CHECK(0<=jj && jj<intersectionClientData->baseZoneJMax);
        LgIndex_t const kk = zValueToKIndex(intersectionPt.Z);
        CHECK(0<=kk && kk<intersectionClientData->baseZoneKMax);

        LgIndex_t ijk = kk*intersectionClientData->baseZoneIJMax + jj*intersectionClientData->baseZoneIMax + ii;
        CHECK(0.0 <= ijk && ijk < intersectionClientData->baseZoneNumPts);

        float intersectFract = float( intersectionPt.X - iIndexToXValue(ii) );
        CHECK(0.0-ijkEpsilon<=intersectFract && intersectFract<=1.0+ijkEpsilon);

        if ( intersectFract <= 0.0F )
            intersectionClientData->iEdgeIntersection[ijk] = 0.0F;
        else if ( intersectFract >= 1.0F )
            intersectionClientData->iEdgeIntersection[ijk] = 1.0F;
        else
            intersectionClientData->iEdgeIntersection[ijk] = intersectFract;
    }
}


void recordJEdgeIntersection(XYZ_s& intersectionPt,
                             void*  clientData)
{
    IntersectionClientData_s* intersectionClientData = (IntersectionClientData_s *)clientData;
    REQUIRE(VALID_REF(intersectionClientData));

    LgIndex_t const jj = yValueToJIndex(intersectionPt.Y);
    // a j-edge is at a certain ik location and goes from j to j+1
    // that is, j-edge is (i,j,k) to (i,j+1,k) thus j must be between 0 and jmax-1
    if ( 0<=jj && jj<intersectionClientData->baseZoneJMax-1 )
    {
        LgIndex_t const ii = yValueToJIndex(intersectionPt.X);
        CHECK(0<=ii && ii<intersectionClientData->baseZoneIMax);
        LgIndex_t const kk = zValueToKIndex(intersectionPt.Z);
        CHECK(0<=kk && kk<intersectionClientData->baseZoneKMax);

        LgIndex_t ijk = kk*intersectionClientData->baseZoneIJMax + jj*intersectionClientData->baseZoneIMax + ii;
        CHECK(0.0 <= ijk && ijk < intersectionClientData->baseZoneNumPts);

        float intersectFract = float( intersectionPt.Y - jIndexToYValue(jj) );
        CHECK(0.0-ijkEpsilon<=intersectFract && intersectFract<=1.0+ijkEpsilon);

        if ( intersectFract <= 0.0F )
            intersectionClientData->jEdgeIntersection[ijk] = 0.0F;
        else if ( intersectFract >= 1.0F )
            intersectionClientData->jEdgeIntersection[ijk] = 1.0F;
        else
            intersectionClientData->jEdgeIntersection[ijk] = intersectFract;
    }
}


void recordKEdgeIntersection(XYZ_s& intersectionPt,
                             void*  clientData)
{
    IntersectionClientData_s* intersectionClientData = (IntersectionClientData_s *)clientData;
    REQUIRE(VALID_REF(intersectionClientData));

    LgIndex_t const kk = zValueToKIndex(intersectionPt.Z);
    // a k-edge is at a certain ij location and goes from k to k+1
    // that is, k-edge is (i,j,k) to (i,j,k+1) thus k must be between 0 and kmax-1
    if ( 0 <= kk && kk < intersectionClientData->baseZoneKMax-1 )
    {
        LgIndex_t const ii = xValueToIIndex(intersectionPt.X);
        CHECK(0<=ii && ii<intersectionClientData->baseZoneIMax);
        LgIndex_t const jj = yValueToJIndex(intersectionPt.Y);
        CHECK(0<=jj && jj<intersectionClientData->baseZoneJMax);

        LgIndex_t ijk = kk*intersectionClientData->baseZoneIJMax + jj*intersectionClientData->baseZoneIMax + ii;
        CHECK(0.0 <= ijk && ijk < intersectionClientData->baseZoneNumPts);

        float intersectFract = float( intersectionPt.Z - kIndexToZValue(kk) );
        CHECK(0.0-ijkEpsilon<=intersectFract && intersectFract<=1.0+ijkEpsilon);

        if ( intersectFract <= 0.0F )
            intersectionClientData->kEdgeIntersection[ijk] = 0.0F;
        else if ( intersectFract >= 1.0F )
            intersectionClientData->kEdgeIntersection[ijk] = 1.0F;
        else
            intersectionClientData->kEdgeIntersection[ijk] = intersectFract;
    }
}

inline void saveRSTPtInRSTArray(FloatRST_s& ptRST,
                                double      rEpsilon, // how close to edge to allow
                                double      sEpsilon, // how close to edge to allow
                                double      tEpsilon, // how close to edge to allow
                                LgIndex_t   baseZoneI,
                                LgIndex_t   baseZoneJ,
                                LgIndex_t   baseZoneK,
                                LgIndex_t   baseZoneIMax,
                                LgIndex_t   baseZoneJMax,
                                LgIndex_t   baseZoneKMax,
                                LgIndex_t*  offsetArray, // offsets into rstValues
                                FloatRST_s* rstValues,
                                LgIndex_t*  rstNumValues,
                                LgIndex_t   rstMaxValues)
{
    REQUIRE(0<=baseZoneI && baseZoneI<baseZoneIMax);
    REQUIRE(0<=baseZoneJ && baseZoneJ<baseZoneJMax);
    REQUIRE(0<=baseZoneK && baseZoneK<baseZoneKMax);
#if 0
    CHECK(-1.0-2*ijkEpsilon<=rr && rr<=1.0+2*ijkEpsilon);
    CHECK(-1.0-2*ijkEpsilon<=ss && ss<=1.0+2*ijkEpsilon);
    CHECK(-1.0-2*ijkEpsilon<=tt && tt<=1.0+2*ijkEpsilon);

    rr /= float(1.0+2*ijkEpsilon);
    ss /= float(1.0+2*ijkEpsilon);
    tt /= float(1.0+2*ijkEpsilon);

    CHECK(-1.0<=rr && rr<=1.0);
    CHECK(-1.0<=ss && ss<=1.0);
    CHECK(-1.0<=tt && tt<=1.0);
#endif

    // only add points near enough to the center to be good
    if ( -1.0+rEpsilon <= ptRST.rr && ptRST.rr <= 1.0-rEpsilon &&
         -1.0+sEpsilon <= ptRST.ss && ptRST.ss <= 1.0-sEpsilon &&
         -1.0+tEpsilon <= ptRST.tt && ptRST.tt <= 1.0-tEpsilon)
    {
        LgIndex_t baseZoneIJK = (baseZoneK*baseZoneJMax + baseZoneJ)*baseZoneIMax + baseZoneI;
        CHECK(0<=baseZoneIJK && baseZoneIJK<baseZoneIMax*baseZoneJMax*baseZoneKMax);

        LgIndex_t valueOffset = offsetArray[baseZoneIJK];
        Boolean_t updatePt;
        if ( valueOffset==noPt )
        {
            // add a new point to the arrays
            offsetArray[baseZoneIJK] = *rstNumValues;
            (*rstNumValues)++;
            CHECK(*rstNumValues < rstMaxValues);
            valueOffset = offsetArray[baseZoneIJK];
            updatePt = TRUE;
        }
        else
        {
            // if previous point already added, we find the one closest to the center
            double curDist = ptRST.rr*ptRST.rr + ptRST.ss*ptRST.ss + ptRST.tt*ptRST.tt;
            double prevR = rstValues[valueOffset].rr;
            double prevS = rstValues[valueOffset].ss;
            double prevT = rstValues[valueOffset].tt;
            double prevDist = prevR*prevR + prevS*prevS + prevT*prevT;
            updatePt = ( curDist < prevDist );
        }

        if ( updatePt )
        {
            rstValues[valueOffset].rr = ptRST.rr;
            rstValues[valueOffset].ss = ptRST.ss;
            rstValues[valueOffset].tt = ptRST.tt;
        }
    }
}

inline void saveXYZPtInRSTArray(XYZ_s&      ptXYZ,
                                double      rEpsilon, // how close to edge to allow
                                double      sEpsilon, // how close to edge to allow
                                double      tEpsilon, // how close to edge to allow
                                LgIndex_t   baseZoneIMax,
                                LgIndex_t   baseZoneJMax,
                                LgIndex_t   baseZoneKMax,
                                LgIndex_t*  offsetArray, // offsets into rstValues
                                FloatRST_s* rstValues,
                                LgIndex_t*  rstNumValues,
                                LgIndex_t   rstMaxValues)
{
    LgIndex_t baseZoneI = xValueToIIndex(ptXYZ.X);
    LgIndex_t baseZoneJ = yValueToJIndex(ptXYZ.Y);
    LgIndex_t baseZoneK = zValueToKIndex(ptXYZ.Z);

    CHECK(-1<=baseZoneI && baseZoneI<baseZoneIMax); // may reach max due to points on edge
    CHECK(-1<=baseZoneJ && baseZoneJ<baseZoneJMax);
    CHECK(-1<=baseZoneK && baseZoneK<baseZoneKMax);

    if ( 0 <= baseZoneI && baseZoneI < baseZoneIMax && 
         0 <= baseZoneJ && baseZoneJ < baseZoneJMax && 
         0 <= baseZoneK && baseZoneK < baseZoneKMax )
    {
        // create brick natural coordinates (from -1 to 1)
        FloatRST_s ptRST;
        ptRST.rr = float( 2.0*(ptXYZ.X-(double)iIndexToXValue(baseZoneI)) - 1.0 );
        ptRST.ss = float( 2.0*(ptXYZ.Y-(double)jIndexToYValue(baseZoneJ)) - 1.0 );
        ptRST.tt = float( 2.0*(ptXYZ.Z-(double)kIndexToZValue(baseZoneK)) - 1.0 );

        saveRSTPtInRSTArray(ptRST,
            rEpsilon,
            sEpsilon,
            tEpsilon,
            baseZoneI,
            baseZoneJ,
            baseZoneK,
            baseZoneIMax,
            baseZoneJMax,
            baseZoneKMax,
            offsetArray,
            rstValues,
            rstNumValues,
            rstMaxValues);
    }
}

void recordCenterPtIntersection(XYZ_s& intersectionPt,
                                void*  clientData)
{
    IntersectionClientData_s* intersectionClientData = (IntersectionClientData_s *)clientData;
    REQUIRE(VALID_REF(intersectionClientData));

    saveXYZPtInRSTArray(intersectionPt,
                     0.0, // rEpsilon
                     0.0, // sEpsilon
                     0.0, // tEpsilon
                     intersectionClientData->baseZoneIMax,
                     intersectionClientData->baseZoneJMax,
                     intersectionClientData->baseZoneKMax,
                     intersectionClientData->centerPtOffset,
                     intersectionClientData->rstValues,
                     intersectionClientData->rstNumValues,
                     intersectionClientData->rstMaxValues);
}


/*
* BondBundlesCreateVolumeZones: For a given BaseZone, cycle through the
* Bonds and create a volume zone contained within the PlsAtomRingCage
* and MnsAtomRingCage surface.
* 1. Find the Bond-PlsAtom and Bond-MnsAtom directions.
* 2. For each point on the Bond-Ring-Cage surface zones, find the
*    intersection of lines projecting in the Bond-MnsAtom and Bond-PlsAtom
*    directions with the MnsAtomRingCage and PlsAtomRingCage surfaces.
* 3. The K=KAve "plane" of the volume zone is the Bond-Ring-Cage surfaces.
*    Successive planes, on either side of the K=KAve plane, are created by
*    evenly spacing points along the projection lines between the Bond-Ring-Cage
*    surface and the intersction points.
*
* Finally, add the new volume zone numbers to the bondBundle structures.
*
* param baseZoneNum
*     Number of the zone for which the topology is computed.
* param volCritPoints
*     Critical point structure for the volume critical points.
* param BondBundlesList
*     ArrayList of pointers to bondBundle data structures for each bond.
*
* return
*     TRUE if successful, FALSE if there were errors.
*/
Boolean_t BondBundlesCreateVolumeZones(const EntIndex_t    baseZoneNum,
                                       const CritPoints_pa volCritPoints,
                                       BondBundles_pa      bondBundles,
                                       double              xMin,
                                       double              xMax,
                                       double              yMin,
                                       double              yMax,
                                       double              zMin,
                                       double              zMax)
{
    Boolean_t isOk = TRUE;

    REQUIRE(baseZoneNum > 0);
    REQUIRE(CritPointsIsValid(volCritPoints));
    REQUIRE(BondBundlesIsValid(bondBundles));

    // Get num zones in dataset
    EntIndex_t numZones;
    EntIndex_t numVars;
    TecUtilDataSetGetInfo(NULL, &numZones, &numVars);

    LgIndex_t baseZoneIMax;
    LgIndex_t baseZoneJMax;
    LgIndex_t baseZoneKMax;
    TecUtilZoneGetIJK(baseZoneNum, &baseZoneIMax, &baseZoneJMax, &baseZoneKMax);

    // some useful values pre-calculated
    LgIndex_t const baseZoneIMaxM1 = baseZoneIMax-1;
    LgIndex_t const baseZoneJMaxM1 = baseZoneJMax-1;
    LgIndex_t const baseZoneKMaxM1 = baseZoneKMax-1;
    LgIndex_t const baseZoneIMaxM2 = baseZoneIMax-2;
    LgIndex_t const baseZoneJMaxM2 = baseZoneJMax-2;
    LgIndex_t const baseZoneKMaxM2 = baseZoneKMax-2;
    LgIndex_t const baseZoneIJMax  = baseZoneJMax*baseZoneIMax;
    LgIndex_t const baseZoneNumPts = baseZoneIMax*baseZoneJMax*baseZoneKMax;

    // Get the number of bonds
    LgIndex_t numBonds = CritPointsGetEndOffset(volCritPoints, (char)(-1)) -
                         CritPointsGetBegOffset(volCritPoints, (char)(-1));

    LgIndex_t numAtoms = CritPointsGetEndOffset(volCritPoints, (char)(-3)) -
                         CritPointsGetBegOffset(volCritPoints, (char)(-3));

    // Must have already populated the bondBundles surface ArrLists
    REQUIRE(BondBundlesGetCount(bondBundles) == numBonds);

    CHECK(numBonds < 128); // change size of BundleInt_t and value of bundleIntFieldDataType

    // create new variable to eventually hold bundle ids, but set to no-bundle for all current zones
    if (isOk)
    {
        FieldDataType_e* varDataType = ALLOC_ARRAY(numZones, FieldDataType_e, "VarDataType");
        isOk = isOk && (varDataType!=NULL);
        if (isOk)
        {
            for (EntIndex_t iz = 1; iz <= numZones; iz++)
                varDataType[iz-1] = bundleIntFieldDataType;
            isOk = TecUtilDataSetAddVar("bondBundleId", varDataType);
            if (isOk)
                numVars++;
        }
        if (varDataType != NULL)
            FREE_ARRAY(varDataType, "varDataType");
    }

    // Create an array to hold the bundle number for each base-zone grid point
    BundleInt_t* whichBondBundle = ALLOC_ARRAY(baseZoneNumPts, BundleInt_t, "ptInBundle");

    // Map nodes of base zone points into bond zones points
    LgIndex_t* nodeToBondPt = ALLOC_ARRAY(baseZoneNumPts, LgIndex_t, "nodeToBondPt");

    // Map nodes along base zones edges (i,j,&k) to bond zone points (note: for simplicity, the
    // same indexing as nodes is used)
    LgIndex_t* iEdgeToBondPt = ALLOC_ARRAY(baseZoneNumPts, LgIndex_t, "iEdgeToBondPt");
    LgIndex_t* jEdgeToBondPt = ALLOC_ARRAY(baseZoneNumPts, LgIndex_t, "jEdgeToBondPt");
    LgIndex_t* kEdgeToBondPt = ALLOC_ARRAY(baseZoneNumPts, LgIndex_t, "kEdgeToBondPt");
    // Location of intersection across edge (fraction) or -1.0 if none.
    float* iEdgeIntersection = ALLOC_ARRAY(baseZoneNumPts, float, "iEdgeIntersection");
    float* jEdgeIntersection = ALLOC_ARRAY(baseZoneNumPts, float, "jEdgeIntersection");
    float* kEdgeIntersection = ALLOC_ARRAY(baseZoneNumPts, float, "kEdgeIntersection");

    // Map nodes on each ff (ij, jk, and ki) to bond zone points (note: for simplicity, the
    // same indexing as nodes is used)
    LgIndex_t* ijFaceToBondPt = ALLOC_ARRAY(baseZoneNumPts, LgIndex_t, "ijFaceToBondPt");
    LgIndex_t* jkFaceToBondPt = ALLOC_ARRAY(baseZoneNumPts, LgIndex_t, "jkFaceToBondPt");
    LgIndex_t* kiFaceToBondPt = ALLOC_ARRAY(baseZoneNumPts, LgIndex_t, "kiFaceToBondPt");
    // Map nodes on each ff (ij, jk, and ki) to "rstValues"
    LgIndex_t* ijFacePtOffset = ALLOC_ARRAY(baseZoneNumPts, LgIndex_t, "ijFacePtOffset");
    LgIndex_t* jkFacePtOffset = ALLOC_ARRAY(baseZoneNumPts, LgIndex_t, "jkFacePtOffset");
    LgIndex_t* kiFacePtOffset = ALLOC_ARRAY(baseZoneNumPts, LgIndex_t, "kiFacePtOffset");

    // map nodes in center of cell to bond zone points (again for simplicity we use the same
    // indexing as nodes, thus leaving some gaps in the array)
    LgIndex_t* centerToBondPt = ALLOC_ARRAY(baseZoneNumPts, LgIndex_t, "centerToBondPt");
    // indirection into the "rstValues" array for the center of cells
    LgIndex_t* centerPtOffset = ALLOC_ARRAY(baseZoneNumPts, LgIndex_t, "centerPtOffset");

    // natural coordinates for the centers and faces
    LgIndex_t rstNumValues = 0;
    LgIndex_t rstMaxValues = baseZoneNumPts; // TODO: a more dynamic approach would be good
    FloatRST_s* rstValues = ALLOC_ARRAY(rstMaxValues, FloatRST_s, "rstValues");
    if ( whichBondBundle == NULL   ||
         nodeToBondPt == NULL      ||
         iEdgeToBondPt == NULL     ||
         jEdgeToBondPt == NULL     ||
         kEdgeToBondPt == NULL     ||
         iEdgeIntersection == NULL ||
         jEdgeIntersection == NULL ||
         kEdgeIntersection == NULL ||
         ijFaceToBondPt == NULL    ||
         jkFaceToBondPt == NULL    ||
         kiFaceToBondPt == NULL    ||
         ijFacePtOffset == NULL    ||
         jkFacePtOffset == NULL    ||
         kiFacePtOffset == NULL    ||
         centerToBondPt == NULL    ||
         centerPtOffset == NULL    ||
         rstValues == NULL )
        isOk = FALSE;
    if ( isOk )
    {
        for ( int pt = 0; pt < baseZoneNumPts; pt++ )
            whichBondBundle[pt] = noBondBundle;
        // other arrays initialized as used
    }

    // for each bundle calculate which points inside and create a bond bundle zone
    if ( isOk )
    {
        EntIndex_t xVarNum = TecUtilVarGetNumByAssignment('X');
        EntIndex_t yVarNum = TecUtilVarGetNumByAssignment('Y');
        EntIndex_t zVarNum = TecUtilVarGetNumByAssignment('Z');

        for (LgIndex_t bondOffset = 0; isOk && bondOffset < numBonds; bondOffset++)
        {
            TecUtilDataLoadBegin();

            FieldData_pa baseZoneXD = TecUtilDataValueGetReadableNativeRef(baseZoneNum, xVarNum);
            FieldData_pa baseZoneYD = TecUtilDataValueGetReadableNativeRef(baseZoneNum, yVarNum);
            FieldData_pa baseZoneZD = TecUtilDataValueGetReadableNativeRef(baseZoneNum, zVarNum);

            // exploit the simple nature of the IJK base zone to get mins and maxes
            double const baseZoneXMin = TecUtilDataValueGetByRef(baseZoneXD, 1); // 1-based for first point (i=0,j=0,k=0)
            double const baseZoneYMin = TecUtilDataValueGetByRef(baseZoneYD, 1);
            double const baseZoneZMin = TecUtilDataValueGetByRef(baseZoneZD, 1);
            double const baseZoneXMax = TecUtilDataValueGetByRef(baseZoneXD, baseZoneNumPts); // 1-based for last point (i=imax,j=jmax,k=kmax)
            double const baseZoneYMax = TecUtilDataValueGetByRef(baseZoneYD, baseZoneNumPts);
            double const baseZoneZMax = TecUtilDataValueGetByRef(baseZoneZD, baseZoneNumPts);

            IntersectionClientData_s intersectionClientData;
            intersectionClientData.baseZoneIMax = baseZoneIMax;
            intersectionClientData.baseZoneJMax = baseZoneJMax;
            intersectionClientData.baseZoneKMax = baseZoneKMax;
            intersectionClientData.iEdgeIntersection = iEdgeIntersection;
            intersectionClientData.jEdgeIntersection = jEdgeIntersection;
            intersectionClientData.kEdgeIntersection = kEdgeIntersection;
            intersectionClientData.ijFacePtOffset = ijFacePtOffset;
            intersectionClientData.jkFacePtOffset = jkFacePtOffset;
            intersectionClientData.kiFacePtOffset = kiFacePtOffset;
            intersectionClientData.centerPtOffset = centerPtOffset;
            intersectionClientData.rstValues = rstValues;
            intersectionClientData.rstNumValues = &rstNumValues;
            intersectionClientData.rstMaxValues = rstMaxValues;

            // a few precalculated things make things easier in the callback
            intersectionClientData.baseZoneIJMax = baseZoneIJMax;
            intersectionClientData.baseZoneNumPts = baseZoneNumPts;

            BundleInt_t const bondNumber = bondOffset+1;
            BondBundle_pa bondBundle = BondBundlesGetBondBundle(bondBundles, bondOffset);

            EntIndex_t plsAtomZoneNum = bondBundle->BondPlsAtomLineZnNum;

            FieldData_pa plsAtomXD = TecUtilDataValueGetReadableRef(plsAtomZoneNum, xVarNum);
            FieldData_pa plsAtomYD = TecUtilDataValueGetReadableRef(plsAtomZoneNum, yVarNum);
            FieldData_pa plsAtomZD = TecUtilDataValueGetReadableRef(plsAtomZoneNum, zVarNum);

            XYZ_s bondPt;
            bondPt.X = TecUtilDataValueGetByRef(plsAtomXD, 1);
            bondPt.Y = TecUtilDataValueGetByRef(plsAtomYD, 1);
            bondPt.Z = TecUtilDataValueGetByRef(plsAtomZD, 1);

            // find cell that contains bond point
            LgIndex_t bondPtI = xValueToIIndex(bondPt.X);
            CHECK(0<=bondPtI && bondPtI<baseZoneIMaxM1);
            LgIndex_t bondPtJ = yValueToJIndex(bondPt.Y);
            CHECK(0<=bondPtJ && bondPtJ<baseZoneJMaxM1);
            LgIndex_t bondPtK = zValueToKIndex(bondPt.Z);
            CHECK(0<=bondPtK && bondPtK<baseZoneKMaxM1);

            // could try other directions (I+1, J+1, & K+1) but this is working for now.
            LgIndex_t baseZoneIJKOffset = bondPtK*baseZoneIJMax + bondPtJ*baseZoneIMax + bondPtI;
            if ( checkBondDirection(bondPt, bondNumber, baseZoneIJKOffset, baseZoneXD, baseZoneYD, baseZoneZD, bondBundle, whichBondBundle) )
            {
                whichBondBundle[baseZoneIJKOffset] = bondNumber;

                // first find intersections of bond bundle with grid
                {
                    char message[256];
                    sprintf_s(message, "Finding intersections for bond bundle %d", bondNumber);
                    TecUtilStatusStartPercentDone(message, TRUE, TRUE);

                    // initialize intersection arrays (we could loop over edges, but what's the point)
                    for ( int pt = 0; pt < baseZoneNumPts; pt++ )
                    {
                        iEdgeIntersection[pt] = noIntersection;
                        jEdgeIntersection[pt] = noIntersection;
                        kEdgeIntersection[pt] = noIntersection;
                    }

                    // find intersections from jk-plane (i-dir)
                    int count = 0; // for status line operations
                    double const countFactor = 100.0 / (baseZoneKMax*baseZoneJMax + baseZoneJMax*baseZoneIMax + baseZoneKMax*baseZoneIMax);

                    XYZ_s indexDirStartPt;
                    XYZ_s indexDirEndPt;

                    indexDirStartPt.X = iIndexToXValue(-1); // intentionally start outside the box
                    indexDirEndPt.X = iIndexToXValue(baseZoneIMax+1); // intentionally end outside the box

                    for ( int kk = 0; kk < baseZoneKMax; kk++ )
                    {
                        indexDirStartPt.Z = kIndexToZValue(kk);
                        indexDirEndPt.Z = indexDirStartPt.Z;
                        for ( int jj = 0; jj < baseZoneJMax; jj++ )
                        {
                            indexDirStartPt.Y = jIndexToYValue(jj);
                            indexDirEndPt.Y = indexDirStartPt.Y;
                            BondBundleProcessSurfaceLineSegIntersections(bondBundle, indexDirStartPt, indexDirEndPt,
                                                                         recordIEdgeIntersection, &intersectionClientData);
                        }
                        count += baseZoneJMax;
                        TecUtilStatusCheckPercentDone((int)(count*countFactor));
                    }

                    // find intersections from ki-plane (j-dir)
                    indexDirStartPt.Y = jIndexToYValue(-1); // intentionally start outside the box
                    indexDirEndPt.Y = jIndexToYValue(baseZoneJMax+1); // intentionally end outside the box
                    for ( int kk = 0; kk < baseZoneKMax; kk++ )
                    {
                        indexDirStartPt.Z = kIndexToZValue(kk);
                        indexDirEndPt.Z = indexDirStartPt.Z;
                        for ( int ii = 0; ii < baseZoneIMax; ii++ )
                        {
                            indexDirStartPt.X = iIndexToXValue(ii);
                            indexDirEndPt.X = indexDirStartPt.X;
                            BondBundleProcessSurfaceLineSegIntersections(bondBundle, indexDirStartPt, indexDirEndPt,
                                                                         recordJEdgeIntersection, &intersectionClientData);
                        }
                        count += baseZoneIMax;
                        TecUtilStatusCheckPercentDone((int)(count*countFactor));
                    }

                    // find intersections from ij-plane (k-dir)
                    indexDirStartPt.Z = kIndexToZValue(-1); // intentionally start outside the box
                    indexDirEndPt.Z = kIndexToZValue(baseZoneKMax+1); // intentionally end outside the box
                    for ( int jj = 0; jj < baseZoneJMax; jj++ )
                    {
                        indexDirStartPt.Y = jIndexToYValue(jj);
                        indexDirEndPt.Y = indexDirStartPt.Y;
                        for ( int ii = 0; ii < baseZoneIMax; ii++ )
                        {
                            indexDirStartPt.X = iIndexToXValue(ii);
                            indexDirEndPt.X = indexDirStartPt.X;
                            BondBundleProcessSurfaceLineSegIntersections(bondBundle, indexDirStartPt, indexDirEndPt,
                                                                         recordKEdgeIntersection, &intersectionClientData);
                        }
                        count += baseZoneIMax;
                        TecUtilStatusCheckPercentDone((int)(count*countFactor));
                    }

                    TecUtilStatusFinishPercentDone();
                }

                // find all points of the bond bundle's atom-ring-cage surfaces and add them to the
                // center points and use these if we can
                {
                    char message[256];
                    sprintf_s(message, "Adding atom-ring-cage surface points to bond bundle %d", bondNumber);
                    TecUtilStatusStartPercentDone(message, TRUE, TRUE);
                    int count = 0;
                    int const numPasses = 2; // for now
                    double const countFactor = 100.0/double(numPasses);

                    // reset center and ff RST cash & offsets
                    rstNumValues = 0;
                    for ( LgIndex_t pt = 0; pt < baseZoneNumPts; pt++ )
                    {
                        centerPtOffset[pt] = noPt;
                        ijFacePtOffset[pt] = noPt;
                        jkFacePtOffset[pt] = noPt;
                        kiFacePtOffset[pt] = noPt;
                    }

                    for ( int pass = 0; pass < numPasses; pass++ )
                    {
                        SurfType_e const surfType = (pass==0) ? SurfType_PlsAtomRingCage : SurfType_MnsAtomRingCage;
                        LgIndex_t const numSurfs = BondBundleGetSurfCount(bondBundle, surfType);
                        CHECK(numSurfs >= 0);

                        for (LgIndex_t ns = 0; ns < numSurfs; ns++)
                        {
                            EntIndex_t const surfZoneNum = BondBundleGetSurfZnNum(bondBundle, surfType, ns);
                            FieldData_pa surfXD = TecUtilDataValueGetReadableRef(surfZoneNum, xVarNum);
                            FieldData_pa surfYD = TecUtilDataValueGetReadableRef(surfZoneNum, yVarNum);
                            FieldData_pa surfZD = TecUtilDataValueGetReadableRef(surfZoneNum, zVarNum);

                            CHECK(TecUtilZoneIsOrdered(surfZoneNum));
                            LgIndex_t surfIMax, surfJMax, surfKMax;
                            // Use NULL for values we're not interested in
                            TecUtilZoneGetInfo(surfZoneNum, &surfIMax, &surfJMax, &surfKMax,
                                               NULL, NULL, NULL, NULL, NULL,
                                               NULL, NULL, NULL, NULL, NULL);
                            CHECK(surfIMax >= 2 && surfJMax >= 2 && surfKMax == 1);
                            for ( LgIndex_t surfJ = 0; surfJ < surfJMax; surfJ++ )
                            {
                                for ( LgIndex_t surfI = 0; surfI < surfIMax; surfI++ )
                                {
                                    LgIndex_t const surfIJ = surfJ*surfIMax + surfI;
                                    XYZ_s surfPt;
                                    surfPt.X = TecUtilDataValueGetByRef(surfXD, surfIJ+1); // tecutil is 1-based
                                    surfPt.Y = TecUtilDataValueGetByRef(surfYD, surfIJ+1);
                                    surfPt.Z = TecUtilDataValueGetByRef(surfZD, surfIJ+1);
                                    // record the point as the cell center
                                    saveXYZPtInRSTArray(surfPt,
                                                    0.01, // rEpsilon
                                                    0.01, // sEpsilon
                                                    0.01, // tEpsilon
                                                    baseZoneIMax,
                                                    baseZoneJMax,
                                                    baseZoneKMax,
                                                    centerPtOffset,
                                                    rstValues,
                                                    &rstNumValues,
                                                    rstMaxValues);
                                    LgIndex_t const surfPtBaseZoneI = xValueToIIndex(surfPt.X);
                                    LgIndex_t const surfPtBaseZoneJ = yValueToJIndex(surfPt.Y);
                                    LgIndex_t const surfPtBaseZoneK = zValueToKIndex(surfPt.Z);

                                    for ( int surfDir = 0; surfDir < 2; surfDir++ )
                                    {
                                        LgIndex_t surfEndPtIndex;
                                        if ( surfDir == 0 ) // surf i-dir
                                        {
                                            if ( surfI+1 < surfIMax )
                                                surfEndPtIndex = surfIJ+1;
                                            else
                                                continue;
                                        }
                                        else // surf j-dir
                                        {
                                            if ( surfJ+1 < surfJMax )
                                                surfEndPtIndex = surfIJ+surfIMax; // j+1
                                            else
                                                continue;
                                        }
                                        XYZ_s surfEndPt;
                                        surfEndPt.X = TecUtilDataValueGetByRef(surfXD, surfEndPtIndex+1); // tecutil is 1-based
                                        surfEndPt.Y = TecUtilDataValueGetByRef(surfYD, surfEndPtIndex+1);
                                        surfEndPt.Z = TecUtilDataValueGetByRef(surfZD, surfEndPtIndex+1);
                                        LgIndex_t const surfEndPtBaseZoneI = xValueToIIndex(surfEndPt.X);
                                        LgIndex_t const surfEndPtBaseZoneJ = yValueToJIndex(surfEndPt.Y);
                                        LgIndex_t const surfEndPtBaseZoneK = zValueToKIndex(surfEndPt.Z);

                                        LgIndex_t const baseZoneStartI = MIN(surfPtBaseZoneI, surfEndPtBaseZoneI)-1;
                                        LgIndex_t const baseZoneStartJ = MIN(surfPtBaseZoneJ, surfEndPtBaseZoneJ)-1;
                                        LgIndex_t const baseZoneStartK = MIN(surfPtBaseZoneK, surfEndPtBaseZoneK)-1;

                                        LgIndex_t const baseZoneEndI = MAX(surfPtBaseZoneI, surfEndPtBaseZoneI)+1;
                                        LgIndex_t const baseZoneEndJ = MAX(surfPtBaseZoneJ, surfEndPtBaseZoneJ)+1;
                                        LgIndex_t const baseZoneEndK = MAX(surfPtBaseZoneK, surfEndPtBaseZoneK)+1;

                                        XYZ_s corner1;
                                        XYZ_s corner2;
                                        XYZ_s corner3;
                                        XYZ_s corner4;

                                        // ij-plane intersections
                                        corner1.X = baseZoneXMin;
                                        corner1.Y = baseZoneYMin;
                                        corner1.Z = 0; // set in loop below
                                        corner2.X = baseZoneXMin;
                                        corner2.Y = baseZoneYMax;
                                        corner2.Z = 0; // set in loop below
                                        corner3.X = baseZoneXMax;
                                        corner3.Y = baseZoneYMax;
                                        corner3.Z = 0; // set in loop below
                                        corner4.X = baseZoneXMax;
                                        corner4.Y = baseZoneYMin;
                                        corner4.Z = 0; // set in loop below

                                        for ( LgIndex_t kk = baseZoneStartK; kk <= baseZoneEndK; kk++ )
                                        {
                                            double zz = kIndexToZValue(kk);
                                            corner1.Z = zz;
                                            corner2.Z = zz;
                                            corner3.Z = zz;
                                            corner4.Z = zz;

                                            int numIntersections;
                                            XYZ_s intersectionPts[2];
                                            if ( GeomToolsLineSegQuadIntersections(corner1, corner2, corner3, corner4,
                                                                                   surfPt, surfEndPt,
                                                                                   &numIntersections, intersectionPts))
                                            {
                                                for ( int pp = 0; pp < numIntersections; pp++ )
                                                {
                                                    intersectionPts[pp].Z = zz;
                                                    saveXYZPtInRSTArray(intersectionPts[pp],
                                                                        0.01, // rEpsilon
                                                                        0.01, // sEpsilon
                                                                        0.0,  // tEpsilon - allow all the way to edge
                                                                        baseZoneIMax,
                                                                        baseZoneJMax,
                                                                        baseZoneKMax,
                                                                        ijFacePtOffset,
                                                                        rstValues,
                                                                        &rstNumValues,
                                                                        rstMaxValues);
                                                }
                                            }
                                        }
                                        
                                        // jk-plane intersections
                                        corner1.X = 0; // set in loop below
                                        corner1.Y = baseZoneYMin;
                                        corner1.Z = baseZoneZMin;
                                        corner2.X = 0; // set in loop below
                                        corner2.Y = baseZoneYMin;
                                        corner2.Z = baseZoneZMax;
                                        corner3.X = 0; // set in loop below
                                        corner3.Y = baseZoneYMax;
                                        corner3.Z = baseZoneZMax;
                                        corner4.X = 0; // set in loop below
                                        corner4.Y = baseZoneYMax;
                                        corner4.Z = baseZoneZMin;

                                        for ( LgIndex_t ii = baseZoneStartI; ii <= baseZoneEndI; ii++ )
                                        {
                                            double xx = iIndexToXValue(ii);
                                            corner1.X = xx;
                                            corner2.X = xx;
                                            corner3.X = xx;
                                            corner4.X = xx;

                                            int numIntersections;
                                            XYZ_s intersectionPts[2];
                                            if ( GeomToolsLineSegQuadIntersections(corner1, corner2, corner3, corner4,
                                                                                   surfPt, surfEndPt,
                                                                                   &numIntersections, intersectionPts))
                                            {
                                                for ( int pp = 0; pp < numIntersections; pp++ )
                                                {
                                                    intersectionPts[pp].X = xx;
                                                    saveXYZPtInRSTArray(intersectionPts[pp],
                                                                        0.0,  // rEpsilon - allow all the way to edge
                                                                        0.01, // sEpsilon
                                                                        0.01, // tEpsilon
                                                                        baseZoneIMax,
                                                                        baseZoneJMax,
                                                                        baseZoneKMax,
                                                                        jkFacePtOffset,
                                                                        rstValues,
                                                                        &rstNumValues,
                                                                        rstMaxValues);
                                                }
                                            }
                                        }

                                        // ki-plane intersections
                                        corner1.X = baseZoneXMin;
                                        corner1.Y = 0; // set in loop below
                                        corner1.Z = baseZoneZMin;
                                        corner2.X = baseZoneXMin;
                                        corner2.Y = 0; // set in loop below
                                        corner2.Z = baseZoneZMax;
                                        corner3.X = baseZoneXMax;
                                        corner3.Y = 0; // set in loop below
                                        corner3.Z = baseZoneZMax;
                                        corner4.X = baseZoneXMax;
                                        corner4.Y = 0; // set in loop below
                                        corner4.Z = baseZoneZMin;

                                        for ( LgIndex_t jj = baseZoneStartJ; jj <= baseZoneEndJ; jj++ )
                                        {
                                            double yy = jIndexToYValue(jj);
                                            corner1.Y = yy;
                                            corner2.Y = yy;
                                            corner3.Y = yy;
                                            corner4.Y = yy;

                                            int numIntersections;
                                            XYZ_s intersectionPts[2];
                                            if ( GeomToolsLineSegQuadIntersections(corner1, corner2, corner3, corner4,
                                                                                   surfPt, surfEndPt,
                                                                                   &numIntersections, intersectionPts))
                                            {
                                                for ( int pp = 0; pp < numIntersections; pp++ )
                                                {
                                                    intersectionPts[pp].Y = yy;
                                                    saveXYZPtInRSTArray(intersectionPts[pp],
                                                                        0.01, // rEpsilon
                                                                        0.0,  // sEpsilon - allow all the way to edge
                                                                        0.01, // tEpsilon
                                                                        baseZoneIMax,
                                                                        baseZoneJMax,
                                                                        baseZoneKMax,
                                                                        kiFacePtOffset,
                                                                        rstValues,
                                                                        &rstNumValues,
                                                                        rstMaxValues);
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                        count++;
                        TecUtilStatusCheckPercentDone((int)(count*countFactor));
                    }
                    TecUtilStatusFinishPercentDone();
                }

                // now find all points inside the bond bundle
                int pass = 0;
                int numNewPointsAdded = 1; // perform loop at least once
                while ( numNewPointsAdded > 0 )
                {
                    numNewPointsAdded = 0;
                    pass++;

                    char message[256];
                    sprintf_s(message, "Determining points of bond bundle %d, pass %d", bondNumber, pass);
                    TecUtilStatusStartPercentDone(message, TRUE, TRUE);

                    int count = 0;
                    double const countFactor = 100.0/8.0/(double)baseZoneNumPts; // avoid integer overflow

                    // test i-directions
                    for ( LgIndex_t kk = 0; kk < baseZoneKMax; kk++ )
                    {
                        for ( LgIndex_t jj = 0; jj < baseZoneJMax; jj++ )
                        {
                            // sweep i-positive direction
                            for ( LgIndex_t ii = 0; ii < baseZoneIMaxM1; ii++ ) // skip last point
                            {
                                LgIndex_t ijk = kk*baseZoneIJMax + jj*baseZoneIMax + ii;
                                LgIndex_t ip1jk = ijk+1;
                                CHECK(0<=ijk && ijk<baseZoneNumPts);
                                CHECK(0<=ip1jk && ip1jk<baseZoneNumPts);
                                if ( whichBondBundle[ijk] == bondNumber && whichBondBundle[ip1jk] != bondNumber && iEdgeIntersection[ijk] == noIntersection )
                                {
                                    whichBondBundle[ip1jk] = bondNumber;
                                    numNewPointsAdded++;
                                }
                            }
                            count += baseZoneIMax;

                            // sweep i-negative direction
                            for ( LgIndex_t ii = baseZoneIMaxM1; ii >= 1; ii-- ) // skip first point
                            {
                                LgIndex_t ijk = kk*baseZoneIJMax + jj*baseZoneIMax + ii;
                                LgIndex_t im1jk = ijk-1;
                                CHECK(0<=ijk && ijk<baseZoneNumPts);
                                CHECK(0<=im1jk && im1jk<baseZoneNumPts);
                                if ( whichBondBundle[ijk] == bondNumber && whichBondBundle[im1jk] != bondNumber && iEdgeIntersection[im1jk] == noIntersection )
                                {
                                    whichBondBundle[im1jk] = bondNumber;
                                    numNewPointsAdded++;
                                }
                            }
                            count += baseZoneIMax;
                        }
                        TecUtilStatusCheckPercentDone((int)(count*countFactor));
                    }

                    // test j-directions
                    for ( LgIndex_t kk = 0; kk < baseZoneKMax; kk++ )
                    {
                        for ( LgIndex_t ii = 0; ii < baseZoneIMax; ii++ )
                        {
                            // sweep j-positive direction
                            for ( LgIndex_t jj = 0; jj < baseZoneJMaxM1; jj++ ) // skip last point
                            {
                                LgIndex_t ijk = kk*baseZoneIJMax + jj*baseZoneIMax + ii;
                                LgIndex_t ijp1k = ijk+baseZoneIMax;
                                CHECK(0<=ijk && ijk<baseZoneNumPts);
                                CHECK(0<=ijp1k && ijp1k<baseZoneNumPts);
                                if ( whichBondBundle[ijk] == bondNumber && whichBondBundle[ijp1k] != bondNumber && jEdgeIntersection[ijk] == noIntersection )
                                {
                                    whichBondBundle[ijp1k] = bondNumber;
                                    numNewPointsAdded++;
                                }
                            }
                            count += baseZoneJMax;

                            // sweep j-negative direction
                            for ( LgIndex_t jj = baseZoneJMaxM1; jj >= 1; jj-- ) // skip first point
                            {
                                LgIndex_t ijk = kk*baseZoneIJMax + jj*baseZoneIMax + ii;
                                LgIndex_t ijm1k = ijk-baseZoneIMax;
                                CHECK(0<=ijk && ijk<baseZoneNumPts);
                                CHECK(0<=ijm1k && ijm1k<baseZoneNumPts);
                                if ( whichBondBundle[ijk] == bondNumber && whichBondBundle[ijm1k] != bondNumber && jEdgeIntersection[ijm1k] == noIntersection )
                                {
                                    whichBondBundle[ijm1k] = bondNumber;
                                    numNewPointsAdded++;
                                }
                            }
                            count += baseZoneJMax;
                        }
                        TecUtilStatusCheckPercentDone((int)(count*countFactor));
                    }

                    // test k-directions
                    for ( LgIndex_t jj = 0; jj < baseZoneJMax; jj++ )
                    {
                        for ( LgIndex_t ii = 0; ii < baseZoneIMax; ii++ )
                        {
                            // sweep k-positive direction
                            for ( LgIndex_t kk = 0; kk < baseZoneKMaxM1; kk++ ) // skip last point
                            {
                                LgIndex_t ijk = kk*baseZoneIJMax + jj*baseZoneIMax + ii;
                                LgIndex_t ijkp1 = ijk+baseZoneIJMax;
                                CHECK(0<=ijk && ijk<baseZoneNumPts);
                                CHECK(0<=ijkp1 && ijkp1<baseZoneNumPts);
                                if ( whichBondBundle[ijk] == bondNumber && whichBondBundle[ijkp1] != bondNumber && kEdgeIntersection[ijk] == noIntersection )
                                {
                                    whichBondBundle[ijkp1] = bondNumber;
                                    numNewPointsAdded++;
                                }
                            }
                            count += baseZoneKMax;

                            // sweep k-negative direction
                            for ( LgIndex_t kk = baseZoneKMaxM1; kk >= 1; kk-- ) // skip first point
                            {
                                LgIndex_t ijk = kk*baseZoneIJMax + jj*baseZoneIMax + ii;
                                LgIndex_t ijkm1 = ijk-baseZoneIJMax;
                                CHECK(0<=ijk && ijk<baseZoneNumPts);
                                CHECK(0<=ijkm1 && ijkm1<baseZoneNumPts);
                                if ( whichBondBundle[ijk] == bondNumber && whichBondBundle[ijkm1] != bondNumber && kEdgeIntersection[ijkm1] == noIntersection )
                                {
                                    whichBondBundle[ijkm1] = bondNumber;
                                    numNewPointsAdded++;
                                }
                            }
                            count += baseZoneKMax;
                        }
                        TecUtilStatusCheckPercentDone((int)(count*countFactor));
                    }

                    TecUtilStatusFinishPercentDone();
                }

                // make sure we have a cell center point for any cell that isn't entirely inside or outside the bundle
                {
                    char message[256];
                    sprintf_s(message, "Finding cell centers for bond bundle %d", bondNumber);
                    TecUtilStatusStartPercentDone(message, TRUE, TRUE);

                    int count = 0;
                    double const countFactor = 100.0/double(baseZoneKMaxM1);

                    for (LgIndex_t kk = 0; kk < baseZoneKMaxM1; kk++)
                    {
                        LgIndex_t baseZoneKOffset = kk*baseZoneIJMax;
                        for (LgIndex_t jj = 0; jj < baseZoneJMaxM1; jj++)
                        {
                            LgIndex_t baseZoneJKOffset = baseZoneKOffset + jj*baseZoneIMax;
                            for (LgIndex_t ii = 0; ii < baseZoneIMaxM1; ii++)
                            {
                                LgIndex_t const ijk = baseZoneJKOffset+ii;
                                if ( centerPtOffset[ijk] == noPt )
                                {
                                    LgIndex_t const ip1jk     = ijk+1;
                                    LgIndex_t const ip1jp1k   = ip1jk+baseZoneIMax;
                                    LgIndex_t const ijp1k     = ijk+baseZoneIMax;
                                    LgIndex_t const ijkp1     = ijk+baseZoneIJMax;
                                    LgIndex_t const ip1jkp1   = ip1jk+baseZoneIJMax;
                                    LgIndex_t const ip1jp1kp1 = ip1jp1k+baseZoneIJMax;
                                    LgIndex_t const ijp1kp1   = ijp1k+baseZoneIJMax;
                                    LgIndex_t const pp[8] = {ijk, ip1jk, ip1jp1k, ijp1k,
                                                             ijkp1, ip1jkp1, ip1jp1kp1, ijp1kp1};

                                    int cornerCount = 0;
                                    for ( int corner = 0; corner < 8; corner++ )
                                        if ( whichBondBundle[pp[corner]] == bondNumber )
                                            cornerCount++;

                                    if (cornerCount != 0 && cornerCount != 8)
                                    {
                                        // if centerPtOffset isn't set, then there is no bond bundle surface point in the cell,
                                        // so find the intersection of the diagonals of the cell with the bond bundle surface
                                        // this will get cached into centerPtOffset for the other variables
                                        int facePtCount = 0;
                                        FloatRST_s faceRST[6];
                                        if ( ijFacePtOffset[ijk] != noPt )
                                        {
                                            faceRST[facePtCount] = rstValues[ijFacePtOffset[ijk]];
                                            CHECK(faceRST[facePtCount].tt == -1.0F);
                                            facePtCount++;
                                        }
                                        if ( ijFacePtOffset[ijkp1] != noPt )
                                        {
                                            faceRST[facePtCount] = rstValues[ijFacePtOffset[ijkp1]];
                                            CHECK(faceRST[facePtCount].tt == -1.0F);
                                            faceRST[facePtCount].tt = 1.0; // value is for i,j,k+1 cell, set to ijk cell
                                            facePtCount++;
                                        }
                                        if ( jkFacePtOffset[ijk] != noPt )
                                        {
                                            faceRST[facePtCount] = rstValues[jkFacePtOffset[ijk]];
                                            CHECK(faceRST[facePtCount].rr == -1.0F);
                                            facePtCount++;
                                        }
                                        if ( jkFacePtOffset[ip1jk] != noPt )
                                        {
                                            faceRST[facePtCount] = rstValues[jkFacePtOffset[ip1jk]];
                                            CHECK(faceRST[facePtCount].rr == -1.0F);
                                            faceRST[facePtCount].rr = 1.0; // value is for i+1,j,k cell, set to ijk cell
                                            facePtCount++;
                                        }
                                        if ( kiFacePtOffset[ijk] != noPt )
                                        {
                                            faceRST[facePtCount] = rstValues[kiFacePtOffset[ijk]];
                                            CHECK(faceRST[facePtCount].ss == -1.0F);
                                            facePtCount++;
                                        }
                                        if ( kiFacePtOffset[ijp1k] != noPt )
                                        {
                                            faceRST[facePtCount] = rstValues[kiFacePtOffset[ijp1k]];
                                            CHECK(faceRST[facePtCount].ss == -1.0F);
                                            faceRST[facePtCount].ss = 1.0; // value is for i,j+1,k cell, set to ijk cell
                                            facePtCount++;
                                        }
                                        
                                        // if there are two face points, put the center point along the line connecting them
                                        if ( facePtCount == 2 )
                                        {
                                            FloatRST_s centerRST = { 0.0, 0.0, 0.0 };
                                            for ( int ff = 0; ff < facePtCount; ff++ )
                                            {
                                                centerRST.rr += faceRST[ff].rr;
                                                centerRST.ss += faceRST[ff].ss;
                                                centerRST.tt += faceRST[ff].tt;
                                            }
                                            centerRST.rr /= float(facePtCount);
                                            centerRST.ss /= float(facePtCount);
                                            centerRST.tt /= float(facePtCount);
                                            saveRSTPtInRSTArray(centerRST,
                                                                0.001, // rEpsilon
                                                                0.001, // sEpsilon
                                                                0.001, // tEpsilon
                                                                ii,
                                                                jj,
                                                                kk,
                                                                baseZoneIMax,
                                                                baseZoneJMax,
                                                                baseZoneKMax,
                                                                centerPtOffset,
                                                                rstValues,
                                                                &rstNumValues,
                                                                rstMaxValues);
                                            CHECK(centerPtOffset[ijk] != noPt);
                                        }
                                        
                                        if ( centerPtOffset[ijk] == noPt ) // fall back if we still didn't add a point
                                        {
                                            XYZ_s cornerXYZ[8];
                                            for ( int corner = 0; corner < 8; corner++ )
                                            {
                                                cornerXYZ[corner].X = TecUtilDataValueGetByRef(baseZoneXD, pp[corner]+1); // tecUtil is 1-based
                                                cornerXYZ[corner].Y = TecUtilDataValueGetByRef(baseZoneYD, pp[corner]+1); // tecUtil is 1-based
                                                cornerXYZ[corner].Z = TecUtilDataValueGetByRef(baseZoneZD, pp[corner]+1); // tecUtil is 1-based
                                            }

                                            XYZ_s edgeXYZ[12];
                                            for ( int cc = 0; cc < 4; cc++ )
                                            {
                                                int const c2 = (cc+1)%4;
                                                edgeXYZ[cc].X = 0.5 * (cornerXYZ[cc].X + cornerXYZ[c2].X);
                                                edgeXYZ[cc].Y = 0.5 * (cornerXYZ[cc].Y + cornerXYZ[c2].Y);
                                                edgeXYZ[cc].Z = 0.5 * (cornerXYZ[cc].Z + cornerXYZ[c2].Z);

                                                edgeXYZ[cc+4].X = 0.5 * (cornerXYZ[cc].X + cornerXYZ[cc+4].X);
                                                edgeXYZ[cc+4].Y = 0.5 * (cornerXYZ[cc].Y + cornerXYZ[cc+4].Y);
                                                edgeXYZ[cc+4].Z = 0.5 * (cornerXYZ[cc].Z + cornerXYZ[cc+4].Z);

                                                edgeXYZ[cc+8].X = 0.5 * (cornerXYZ[cc+4].X + cornerXYZ[c2+4].X);
                                                edgeXYZ[cc+8].Y = 0.5 * (cornerXYZ[cc+4].Y + cornerXYZ[c2+4].Y);
                                                edgeXYZ[cc+8].Z = 0.5 * (cornerXYZ[cc+4].Z + cornerXYZ[c2+4].Z);
                                            }

                                            BondBundleProcessSurfaceLineSegIntersections(bondBundle, edgeXYZ[0], edgeXYZ[10],
                                                recordCenterPtIntersection, &intersectionClientData);
                                            BondBundleProcessSurfaceLineSegIntersections(bondBundle, edgeXYZ[1], edgeXYZ[11],
                                                recordCenterPtIntersection, &intersectionClientData);
                                            BondBundleProcessSurfaceLineSegIntersections(bondBundle, edgeXYZ[2], edgeXYZ[8],
                                                recordCenterPtIntersection, &intersectionClientData);
                                            BondBundleProcessSurfaceLineSegIntersections(bondBundle, edgeXYZ[3], edgeXYZ[9],
                                                recordCenterPtIntersection, &intersectionClientData);
                                            BondBundleProcessSurfaceLineSegIntersections(bondBundle, edgeXYZ[4], edgeXYZ[6],
                                                recordCenterPtIntersection, &intersectionClientData);
                                            BondBundleProcessSurfaceLineSegIntersections(bondBundle, edgeXYZ[5], edgeXYZ[7],
                                                recordCenterPtIntersection, &intersectionClientData);

                                            BondBundleProcessSurfaceLineSegIntersections(bondBundle, cornerXYZ[0], cornerXYZ[6],
                                                                                         recordCenterPtIntersection, &intersectionClientData);
                                            BondBundleProcessSurfaceLineSegIntersections(bondBundle, cornerXYZ[1], cornerXYZ[7],
                                                                                         recordCenterPtIntersection, &intersectionClientData);
                                            BondBundleProcessSurfaceLineSegIntersections(bondBundle, cornerXYZ[2], cornerXYZ[4],
                                                                                         recordCenterPtIntersection, &intersectionClientData);
                                            BondBundleProcessSurfaceLineSegIntersections(bondBundle, cornerXYZ[3], cornerXYZ[5],
                                                                                         recordCenterPtIntersection, &intersectionClientData);

                                            //CHECK(centerPtOffset[ijk]!=noPt);
                                        }
                                    }
                                }
                            }
                        }
                        count++;
                        TecUtilStatusCheckPercentDone((int)(count*countFactor));
                    }
                    TecUtilStatusFinishPercentDone();
                }

                // create the bond volume zone
                LgIndex_t bondZoneNumPts = 0;
                LgIndex_t bondZoneNumElem = 0;
                if ( isOk )
                {
                    char message[256];
                    sprintf_s(message, "Creating zone for bond bundle %d", bondNumber);
                    TecUtilStatusStartPercentDone(message, TRUE, TRUE);

                    // first count the number of points and create mapping from base zone points to bond zone points
                    // initialize node and edge mapping arrays
                    for ( int pt = 0; pt < baseZoneNumPts; pt++ )
                    {
                        nodeToBondPt[pt] = noPt;
                        iEdgeToBondPt[pt] = noPt;
                        jEdgeToBondPt[pt] = noPt;
                        kEdgeToBondPt[pt] = noPt;
                        centerToBondPt[pt] = noPt; // TODO: with RST offsets, can this be removed
                        ijFaceToBondPt[pt] = noPt;
                        jkFaceToBondPt[pt] = noPt;
                        kiFaceToBondPt[pt] = noPt;
                    }

                    for (LgIndex_t kk = 0; kk < baseZoneKMax; kk++)
                    {
                        LgIndex_t baseZoneKOffset = kk*baseZoneIJMax;
                        for (LgIndex_t jj = 0; jj < baseZoneJMax; jj++)
                        {
                            LgIndex_t baseZoneJKOffset = baseZoneKOffset + jj*baseZoneIMax;
                            for (LgIndex_t ii = 0; ii < baseZoneIMax; ii++)
                            {
                                LgIndex_t ijk = baseZoneJKOffset + ii;

                                // base zone nodes
                                if ( whichBondBundle[ijk] == bondNumber )
                                {
                                    nodeToBondPt[ijk] = bondZoneNumPts;
                                    bondZoneNumPts++;
                                }

                                // i-edge crossings:
                                LgIndex_t ip1jk = ijk+1;
                                if ( iEdgeIntersection[ijk] != noIntersection )
                                {
                                    CHECK(iEdgeIntersection[ijk]!=noIntersection);
                                    iEdgeToBondPt[ijk] = bondZoneNumPts;
                                    bondZoneNumPts++;
                                }

                                // j-edge crossings
                                LgIndex_t ijp1k = ijk+baseZoneIMax;
                                if ( jEdgeIntersection[ijk] != noIntersection )
                                {
                                    CHECK(jEdgeIntersection[ijk]!=noIntersection);
                                    jEdgeToBondPt[ijk] = bondZoneNumPts;
                                    bondZoneNumPts++;
                                }

                                // k-edge crossings
                                LgIndex_t ijkp1 = ijk+baseZoneIJMax;
                                if ( kEdgeIntersection[ijk] != noIntersection )
                                {
                                    CHECK(kEdgeIntersection[ijk]!=noIntersection);
                                    kEdgeToBondPt[ijk] = bondZoneNumPts;
                                    bondZoneNumPts++;
                                }

                                // centers and faces
                                if ( ii < baseZoneIMaxM1 &&
                                     jj < baseZoneJMaxM1 && 
                                     kk < baseZoneKMaxM1 &&
                                     centerPtOffset[ijk] !=noPt )
                                {
                                    centerToBondPt[ijk] = bondZoneNumPts;
                                    bondZoneNumPts++;
                                }

                                if ( ii < baseZoneIMaxM1 &&
                                     jj < baseZoneJMaxM1 &&
                                     ijFacePtOffset[ijk] !=noPt )
                                {
                                    ijFaceToBondPt[ijk] = bondZoneNumPts;
                                    bondZoneNumPts++;
                                }
                                
                                if ( jj < baseZoneJMaxM1 &&
                                     kk < baseZoneKMaxM1 &&
                                     jkFacePtOffset[ijk] !=noPt )
                                {
                                    jkFaceToBondPt[ijk] = bondZoneNumPts;
                                    bondZoneNumPts++;
                                }
                                
                                if ( kk < baseZoneKMaxM1 &&
                                     ii < baseZoneIMaxM1 &&
                                     kiFacePtOffset[ijk] !=noPt )
                                {
                                    kiFaceToBondPt[ijk] = bondZoneNumPts;
                                    bondZoneNumPts++;
                                }
                            }
                        }
                    }
                    CHECK(bondZoneNumPts>0);
                    isOk = isOk && (bondZoneNumPts>0); // just in case

                    // count the number of elements
                    for (LgIndex_t kk = 0; kk < baseZoneKMaxM1; kk++)
                    {
                        LgIndex_t baseZoneKOffset = kk*baseZoneIJMax;
                        for (LgIndex_t jj = 0; jj < baseZoneJMaxM1; jj++)
                        {
                            LgIndex_t baseZoneJKOffset = baseZoneKOffset + jj*baseZoneIMax;
                            for (LgIndex_t ii = 0; ii < baseZoneIMaxM1; ii++)
                            {
                                LgIndex_t const p1 = baseZoneJKOffset + ii;
                                LgIndex_t const p2 = p1+1;
                                LgIndex_t const p3 = p1+baseZoneIMax+1;
                                LgIndex_t const p4 = p1+baseZoneIMax;
                                LgIndex_t const p5 = p1+baseZoneIJMax;
                                LgIndex_t const p6 = p2+baseZoneIJMax;
                                LgIndex_t const p7 = p3+baseZoneIJMax;
                                LgIndex_t const p8 = p4+baseZoneIJMax;
                                createElements(nodeToBondPt[p1], nodeToBondPt[p2], nodeToBondPt[p3], nodeToBondPt[p4],
                                               nodeToBondPt[p5], nodeToBondPt[p6], nodeToBondPt[p7], nodeToBondPt[p8],
                                               iEdgeToBondPt[p1], jEdgeToBondPt[p2], iEdgeToBondPt[p4], jEdgeToBondPt[p1],
                                               kEdgeToBondPt[p1], kEdgeToBondPt[p2], kEdgeToBondPt[p3], kEdgeToBondPt[p4],
                                               iEdgeToBondPt[p5], jEdgeToBondPt[p6], iEdgeToBondPt[p8], jEdgeToBondPt[p5],
                                               ijFaceToBondPt[p1], ijFaceToBondPt[p5],
                                               jkFaceToBondPt[p1], jkFaceToBondPt[p2],
                                               kiFaceToBondPt[p1], kiFaceToBondPt[p4],
                                               centerToBondPt[p1],
                                               NULL/*count only*/,
                                               bondZoneNumElem);
                            }
                        }
                    }
                    CHECK(bondZoneNumElem>0);
                    isOk = isOk && (bondZoneNumElem>0); // just in case

                    // now create the zone
                    FieldDataType_e* varDataType = ALLOC_ARRAY(numVars, FieldDataType_e, "VarDataType");
                    isOk = isOk && (varDataType!=NULL);
                    if (isOk)
                    {
                        for (EntIndex_t iv = 1; iv < numVars; iv++)
                        {
                            if ( TecUtilDataValueGetType(baseZoneNum, iv) == FieldDataType_Double )
                                varDataType[iv-1] = FieldDataType_Double;
                            else
                                varDataType[iv-1] = FieldDataType_Float; // allow for fractional values
                        }
                        varDataType[numVars-1] = bundleIntFieldDataType; // only hold bond bundle values now

                        char zoneNameString[200];
                        sprintf_s(zoneNameString, "Bond %d", bondNumber);

                        isOk = TecUtilDataSetAddZone(zoneNameString, bondZoneNumPts, bondZoneNumElem, 1, ZoneType_FEBrick, varDataType);
                        if (isOk)
                            numZones++;
                    }
                    if (varDataType != NULL)
                        FREE_ARRAY(varDataType, "varDataType");

                    TecUtilStatusFinishPercentDone();
                }

                // set variables
                if ( isOk )
                {
                    char message[256];
                    sprintf_s(message, "Assigning variables for bond bundle zone %d", bondNumber);
                    TecUtilStatusStartPercentDone(message, TRUE, TRUE);
                    int count = 0;
                    double const countFactor = 100.0/double(numVars)/double(baseZoneNumPts); // avoid integer overflow

                    EntIndex_t bondZoneNum = numZones;
                    for (EntIndex_t var = 1; var < numVars; var++) // special handling for var==numVars below
                    {
                        if (TecUtilVarIsEnabled(var))
                        {
                            FieldData_pa baseZoneFD = TecUtilDataValueGetReadableRef(baseZoneNum, var);
                            FieldData_pa bondZoneFD = TecUtilDataValueGetWritableRef(bondZoneNum, var);

                            LgIndex_t bondZoneOffset = 0;
                            // copy values from nodes
                            for (LgIndex_t kk = 0; kk < baseZoneKMax; kk++)
                            {
                                LgIndex_t baseZoneKOffset = kk*baseZoneIJMax;
                                for (LgIndex_t jj = 0; jj < baseZoneJMax; jj++)
                                {
                                    LgIndex_t baseZoneJKOffset = baseZoneKOffset + jj*baseZoneIMax;
                                    for (LgIndex_t ii = 0; ii < baseZoneIMax; ii++)
                                    {
                                        LgIndex_t ijk = baseZoneJKOffset + ii;

                                        // base zone nodes
                                        if ( whichBondBundle[ijk] == bondNumber )
                                        {
                                            CHECK(bondZoneOffset==nodeToBondPt[ijk]);
                                            double val = TecUtilDataValueGetByRef(baseZoneFD, ijk+1); // tecutil is 1-based
                                            TecUtilDataValueSetByRef(bondZoneFD, bondZoneOffset+1, val);
                                            bondZoneOffset++;
                                        }

                                        // i-edge crossings
                                        LgIndex_t ip1jk = ijk+1;
                                        if ( iEdgeToBondPt[ijk] != noPt )
                                        {
                                            CHECK(bondZoneOffset == iEdgeToBondPt[ijk]);
                                            double val1 = TecUtilDataValueGetByRef(baseZoneFD, ijk+1); // tecutil is 1-based
                                            double val2 = TecUtilDataValueGetByRef(baseZoneFD, ip1jk+1);
                                            CHECK(iEdgeIntersection[ijk]!=noIntersection);
                                            double val = (val2-val1)*iEdgeIntersection[ijk] + val1;
                                            TecUtilDataValueSetByRef(bondZoneFD, bondZoneOffset+1, val);
                                            bondZoneOffset++;
                                        }

                                        // j-edge crossings
                                        LgIndex_t ijp1k = ijk+baseZoneIMax;
                                        if ( jEdgeToBondPt[ijk] != noPt )
                                        {
                                            CHECK(bondZoneOffset == jEdgeToBondPt[ijk]);
                                            double val1 = TecUtilDataValueGetByRef(baseZoneFD, ijk+1); // tecutil is 1-based
                                            double val2 = TecUtilDataValueGetByRef(baseZoneFD, ijp1k+1);
                                            CHECK(jEdgeIntersection[ijk]!=noIntersection);
                                            double val = (val2-val1)*jEdgeIntersection[ijk] + val1;
                                            TecUtilDataValueSetByRef(bondZoneFD, bondZoneOffset+1, val);
                                            bondZoneOffset++;
                                        }

                                        // k-edge crossings
                                        LgIndex_t ijkp1 = ijk+baseZoneIJMax;
                                        if ( kEdgeToBondPt[ijk] != noPt )
                                        {
                                            CHECK(bondZoneOffset == kEdgeToBondPt[ijk]);
                                            double val1 = TecUtilDataValueGetByRef(baseZoneFD, ijk+1); // tecutil is 1-based
                                            double val2 = TecUtilDataValueGetByRef(baseZoneFD, ijkp1+1);
                                            CHECK(kEdgeIntersection[ijk]!=noIntersection);
                                            double val = (val2-val1)*kEdgeIntersection[ijk] + val1;
                                            TecUtilDataValueSetByRef(bondZoneFD, bondZoneOffset+1, val);
                                            bondZoneOffset++;
                                        }

                                        // centers and ff nodes
                                        {
                                            for ( int dir = 0; dir < 4; dir++ )
                                            {
                                                LgIndex_t valueOffset;
                                                if ( dir == 0 )
                                                {
                                                    if ( ii >= baseZoneIMaxM1 || 
                                                         jj >= baseZoneJMaxM1 || 
                                                         kk >= baseZoneKMaxM1 )
                                                        continue;
                                                    valueOffset = centerPtOffset[ijk];
                                                }
                                                else if ( dir == 1 )
                                                {
                                                    if ( ii >= baseZoneIMaxM1 || 
                                                         jj >= baseZoneJMaxM1 )
                                                        continue;
                                                    valueOffset = ijFacePtOffset[ijk];
                                                }
                                                else if ( dir == 2 )
                                                {
                                                    if ( jj >= baseZoneJMaxM1 ||
                                                         kk >= baseZoneKMaxM1 )
                                                        continue;
                                                    valueOffset = jkFacePtOffset[ijk];
                                                }
                                                else
                                                {
                                                    if ( kk >= baseZoneKMaxM1 ||
                                                         ii >= baseZoneIMaxM1 )
                                                        continue;
                                                    valueOffset = kiFacePtOffset[ijk];
                                                }
                                                if ( valueOffset != noPt )
                                                {
                                                    double cornerWeight[8] = { 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125 };
                                                    VERIFY(BrickTrilinearWeight(rstValues[valueOffset].rr,
                                                                                rstValues[valueOffset].ss,
                                                                                rstValues[valueOffset].tt,
                                                                                cornerWeight));
                                                    // if at imax, must be on jk-plane, etc.
                                                    CHECK(IMPLICATION(ii==baseZoneIMaxM1, rstValues[valueOffset].rr==-1.0F));
                                                    CHECK(IMPLICATION(jj==baseZoneJMaxM1, rstValues[valueOffset].ss==-1.0F));
                                                    CHECK(IMPLICATION(kk==baseZoneKMaxM1, rstValues[valueOffset].tt==-1.0F));

                                                    LgIndex_t const delI = (ii<baseZoneIMaxM1) ? 1 : 0;
                                                    LgIndex_t const delJ = (jj<baseZoneJMaxM1) ? baseZoneIMax : 0;
                                                    LgIndex_t const delK = (kk<baseZoneKMaxM1) ? baseZoneIJMax : 0;

                                                    LgIndex_t const ip1jk     = ijk+delI;
                                                    LgIndex_t const ijp1k     = ijk+delJ;
                                                    LgIndex_t const ip1jp1k   = ip1jk+delJ;
                                                    LgIndex_t const ijkp1     = ijk+delK;
                                                    LgIndex_t const ip1jkp1   = ip1jk+delK;
                                                    LgIndex_t const ijp1kp1   = ijp1k+delK;
                                                    LgIndex_t const ip1jp1kp1 = ip1jp1k+delK;
                                                    LgIndex_t const pp[8] = {ijk, ip1jk, ip1jp1k, ijp1k,
                                                                             ijkp1, ip1jkp1, ip1jp1kp1, ijp1kp1};

                                                    double val = 0.0;
                                                    for ( int corner = 0; corner < 8; corner++ )
                                                    {
                                                        double cornerVal = TecUtilDataValueGetByRef(baseZoneFD, pp[corner]+1); // tecutil is 1-based
                                                        val += cornerVal*cornerWeight[corner];
                                                    }

                                                    TecUtilDataValueSetByRef(bondZoneFD, bondZoneOffset+1, val);
                                                    bondZoneOffset++;
                                                }
                                            }
                                        }
                                    }
                                }
                                count += baseZoneIJMax;
                                TecUtilStatusCheckPercentDone((int)(count*countFactor));
                            }
                            CHECK(bondZoneOffset==bondZoneNumPts);
                        }
                    }
                    // special handling for var==numVars (the bondBundleID case)
                    CHECK(TecUtilVarIsEnabled(numVars));
                    FieldData_pa bondZoneFD = TecUtilDataValueGetWritableRef(bondZoneNum, numVars);
                    for (LgIndex_t bondZoneOffset = 0; bondZoneOffset < bondZoneNumPts; bondZoneOffset++)
                    {
                        TecUtilDataValueSetByRef(bondZoneFD, bondZoneOffset+1, (double)bondNumber); // tecutil is 1-based
                    }
                    TecUtilStatusFinishPercentDone();


                    // create nodemap
                    if ( isOk )
                    {
                        char message[256];
                        sprintf_s(message, "Assigning node map for bond bundle zone %d", bondNumber);
                        TecUtilStatusStartPercentDone(message, TRUE, TRUE);
                        int count = 0;
                        double const countFactor = 100.0/double(baseZoneNumPts);
                        NodeMap_pa bondZoneNM = TecUtilDataNodeGetWritableRef(bondZoneNum);
                        LgIndex_t bondZoneElem = 0;
                        for (LgIndex_t kk = 0; kk < baseZoneKMaxM1; kk++)
                        {
                            LgIndex_t baseZoneKOffset = kk*baseZoneIJMax;
                            for (LgIndex_t jj = 0; jj < baseZoneJMaxM1; jj++)
                            {
                                LgIndex_t baseZoneJKOffset = baseZoneKOffset + jj*baseZoneIMax;
                                for (LgIndex_t ii = 0; ii < baseZoneIMaxM1; ii++)
                                {
                                    LgIndex_t const p1 = baseZoneJKOffset + ii;
                                    LgIndex_t const p2 = p1+1;
                                    LgIndex_t const p3 = p1+baseZoneIMax+1;
                                    LgIndex_t const p4 = p1+baseZoneIMax;
                                    LgIndex_t const p5 = p1+baseZoneIJMax;
                                    LgIndex_t const p6 = p2+baseZoneIJMax;
                                    LgIndex_t const p7 = p3+baseZoneIJMax;
                                    LgIndex_t const p8 = p4+baseZoneIJMax;
                                    createElements(nodeToBondPt[p1], nodeToBondPt[p2], nodeToBondPt[p3], nodeToBondPt[p4],
                                                   nodeToBondPt[p5], nodeToBondPt[p6], nodeToBondPt[p7], nodeToBondPt[p8],
                                                   iEdgeToBondPt[p1], jEdgeToBondPt[p2], iEdgeToBondPt[p4], jEdgeToBondPt[p1],
                                                   kEdgeToBondPt[p1], kEdgeToBondPt[p2], kEdgeToBondPt[p3], kEdgeToBondPt[p4],
                                                   iEdgeToBondPt[p5], jEdgeToBondPt[p6], iEdgeToBondPt[p8], jEdgeToBondPt[p5],
                                                   ijFaceToBondPt[p1], ijFaceToBondPt[p5],
                                                   jkFaceToBondPt[p1], jkFaceToBondPt[p2],
                                                   kiFaceToBondPt[p1], kiFaceToBondPt[p4],
                                                   centerToBondPt[p1],
                                                   bondZoneNM,
                                                   bondZoneElem);
                                }
                            }
                            count += baseZoneIJMax;
                            TecUtilStatusCheckPercentDone((int)(count*countFactor));
                        }
                        CHECK(bondZoneElem==bondZoneNumElem);
                        TecUtilStatusFinishPercentDone();
                    }
                }
            }
            else
                CHECK(FALSE); // bond point in cell with atom-ring-cage surface, we should probably check other points

            TecUtilDataLoadEnd();
        }
    }

    if (whichBondBundle != NULL)
        FREE_ARRAY(whichBondBundle, "ptInBundle");

    if ( nodeToBondPt != NULL )
        FREE_ARRAY(nodeToBondPt, "nodeToBondPt");

    if ( iEdgeToBondPt != NULL )
        FREE_ARRAY(iEdgeToBondPt, "iEdgeToBondPt");
    if ( jEdgeToBondPt != NULL )
        FREE_ARRAY(jEdgeToBondPt, "jEdgeToBondPt");
    if ( kEdgeToBondPt != NULL )
        FREE_ARRAY(kEdgeToBondPt, "kEdgeToBondPt");
    if ( iEdgeIntersection != NULL )
        FREE_ARRAY(iEdgeIntersection, "iEdgeIntersection");
    if ( jEdgeIntersection != NULL )
        FREE_ARRAY(jEdgeIntersection, "jEdgeIntersection");
    if ( kEdgeIntersection != NULL )
        FREE_ARRAY(kEdgeIntersection, "kEdgeIntersection");

    if ( ijFaceToBondPt != NULL )
        FREE_ARRAY(ijFaceToBondPt, "ijFaceToBondPt");
    if ( jkFaceToBondPt != NULL )
        FREE_ARRAY(jkFaceToBondPt, "jkFaceToBondPt");
    if ( kiFaceToBondPt != NULL )
        FREE_ARRAY(kiFaceToBondPt, "kiFaceToBondPt");
    if ( ijFacePtOffset != NULL )
        FREE_ARRAY(ijFacePtOffset, "ijFacePtOffset");
    if ( jkFacePtOffset != NULL )
        FREE_ARRAY(jkFacePtOffset, "jkFacePtOffset");
    if ( kiFacePtOffset != NULL )
        FREE_ARRAY(kiFacePtOffset, "kiFacePtOffset");

    if ( centerToBondPt != NULL )
        FREE_ARRAY(centerToBondPt, "centerToBondPt");
    if ( centerPtOffset != NULL )
        FREE_ARRAY(centerPtOffset, "centerPtOffset");

    if ( rstValues != NULL )
        FREE_ARRAY(rstValues, "rstValues");

    ENSURE(VALID_BOOLEAN(isOk));
    return isOk;
}

