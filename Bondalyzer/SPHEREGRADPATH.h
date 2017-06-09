/*
*****************************************************************
*****************************************************************
*******                                                  ********
*******    (C) Copyright 1988-2013 by TECPLOT INC.       *******
*******                                                  ********
*****************************************************************
*****************************************************************
*/


#ifndef SPHEREGRADPATH_H_
#define SPHEREGRADPATH_H_

#define MAX_TRI_SUBDIVISIONS 1000


/* private ClientData structure */
/*
typedef struct _GradPath_s
  {
    LgIndex_t  MemUsageReported;
    LgIndex_t  BeginCrtPtNum;
    LgIndex_t  EndCrtPtNum;
    ArrList_pa X;
    ArrList_pa Y;
    ArrList_pa Z;
    ArrList_pa Rho;
  } GradPath_s;

typedef struct _GradPath_s  *GradPath_pa;
*/

/*
 * To refine grid in a region, sub-divide the triangles.
 *
 *    Make new nodes:
 *      a = ( 0 + 1 ) / 2  (mid-point of line 0-1)
 *      b = ( 1 + 2 ) / 2  (mid-point of line 1-2)
 *      c = ( 2 + 0 ) / 2  (mid-point of line 2-0)
 *
 *    Normalize the coordinates of nodes (seed points) a, b, and c
 *
 *    Construct new triangles:
 *       s0: [0, a, c]
 *       s1: [1, b, a]
 *       s2: [2, c, b]
 *       s3: [a, b, c]
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
 */
typedef struct _TriGradPath_s *TriGradPath_pa;

typedef struct _TriGradPath_s
{
    char           DepthInTree;   // Do I need this?
    LgIndex_t      GradPathNum[3]; // Offset of the three corner GradPaths in the SphereGradPath->GradPaths ArrList
    TriGradPath_pa SubTriangles[4];
    char           NeighborTri[3];  // Just at depth=0, Neighbor on face opposite node of same number
} TriGradPath_s;



/* private ClientData structure */
typedef struct _SphereGradPath_s
{
    double         Center[3];
    TriGradPath_pa RootTriangles[20];
    double         Radius;
    ArrList_pa     SeedVectorX; /* List of pointers to X-component of initial vector direction of gradpaths */
    ArrList_pa     SeedVectorY; /* List of pointers to Y-component of initial vector direction of gradpaths */
    ArrList_pa     SeedVectorZ; /* List of pointers to Z-component of initial vector direction of gradpaths */
    ArrList_pa     GradPaths;         /* List of pointers to GradPath structures */
} SphereGradPath_s;

typedef struct _SphereGradPath_s  *SphereGradPath_pa;

Boolean_t         TriGradPathIsValid(TriGradPath_pa TriGradPath);
TriGradPath_pa    TriGradPathAlloc();
void              TriGradPathDealloc(TriGradPath_pa *TriGradPath);
void              TriGradPathClear(TriGradPath_pa TriGradPath);
LgIndex_t         TriGradPathGetGPIndex(TriGradPath_pa TriGradPath,
                                        LgIndex_t      CornerNumber);
Boolean_t         TriGradPathSetGPIndex(TriGradPath_pa TriGradPath,
                                        LgIndex_t      CornerNumber,  // 0-2
                                        LgIndex_t      GradPathIndex);
TriGradPath_pa    TriGradPathGetSubTri(TriGradPath_pa TriGradPath,
                                       LgIndex_t      SubTriIndex);
Boolean_t         TriGradPathSetSubTri(TriGradPath_pa TriGradPath,
                                       LgIndex_t      SubTriIndex,   // 0-3
                                       TriGradPath_pa SubTriGradPath);
Boolean_t         TriGradPathSubdivide(const ZoneVarInfo_pa  ZoneVarInfo,
                                       const CritPoints_pa   CritPoints,
                                       TriGradPath_pa    TriGradPath,
                                       const float           CPTolerance,
                                       SphereGradPath_pa SphereGradPath);


Boolean_t         SphereGradPathIsValid(SphereGradPath_pa SphereGradPath);
SphereGradPath_pa SphereGradPathAlloc();
LgIndex_t         SphereGradPathGetGPCount(SphereGradPath_pa SphereGradPath);
void              SphereGradPathDealloc(SphereGradPath_pa *SphereGradPath);
void              SphereGradPathClear(SphereGradPath_pa SphereGradPath);
GradPath_pa       SphereGradPathGetGP(SphereGradPath_pa SphereGradPath,
                                      LgIndex_t         Offset);
XYZ_s             SphereGradPathGetSeedDir(SphereGradPath_pa SphereGradPath,
                                           LgIndex_t         Offset);
Boolean_t         SphereGradPathSetGP(SphereGradPath_pa SphereGradPath,
                                      LgIndex_t         Offset,
                                      XYZ_s             RelNormSeedVector,
                                      GradPath_pa       GradPath);
Boolean_t         SphereGradPathInsertGP(SphereGradPath_pa SphereGradPath,
                                         LgIndex_t         Offset,
                                         XYZ_s             RelNormSeedVector,
                                         GradPath_pa       GradPath);
Boolean_t         SphereGradPathRemoveGP(SphereGradPath_pa SphereGradPath,
                                         LgIndex_t   Offset);
Boolean_t         SphereGradPathAppendGP(SphereGradPath_pa SphereGradPath,
                                         XYZ_s             RelNormSeedVector,
                                         GradPath_pa       GradPath);
TriGradPath_pa    SphereGradPathGetRtTri(SphereGradPath_pa SphereGradPath,
                                         LgIndex_t         Offset);
Boolean_t         SphereGradPathSetRtTri(SphereGradPath_pa SphereGradPath,
                                         LgIndex_t         Offset,
                                         TriGradPath_pa    TriGradPath);
LgIndex_t         SphereGradPathGetBegCPNum(SphereGradPath_pa SphereGradPath);
// Boolean_t SphereGradPathResample(SphereGradPath_pa SphereGradPath,
//                                 LgIndex_t         NumPointsNew);
Boolean_t SphereGradPathAdd(const ZoneVarInfo_pa  ZoneVarInfo,
                            const LgIndex_t       BeginCrtPtNum,
                            const CritPoints_pa   CritPoints,
                            const StreamDir_e     PathDir,
                            const float           CPTolerance,
                            SphereGradPath_pa     SphereGradPath);
GradPath_pa SphereGradPath2CPPath(const ZoneVarInfo_pa  ZoneVarInfo,
                                  const LgIndex_t       BeginCrtPtNum,
                                  const LgIndex_t       EndCrtPtNum,
                                  const CritPoints_pa   CritPoints,
                                  const StreamDir_e     PathDir,
                                  const float           CPTolerance,
                                  SphereGradPath_pa     SphereGradPath);

#endif /* SPHEREGRADPATH_H_ */