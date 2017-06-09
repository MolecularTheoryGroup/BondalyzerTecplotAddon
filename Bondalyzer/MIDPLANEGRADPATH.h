/*
*****************************************************************
*****************************************************************
*******                                                  ********
*******    (C) Copyright 1988-2013 by TECPLOT INC.       *******
*******                                                  ********
*****************************************************************
*****************************************************************
*/


#ifndef MIDPLANEGRADPATH_H_
#define MIDPLANEGRADPATH_H_


/* private ClientData structure */
typedef struct _MidPlnGradPath_s
{
    double         Position[3];
    double         Normal[3];
    ArrList_pa     BBPIntersectX;  /* List of X-positions for bundle bounding path intersections */
    ArrList_pa     BBPIntersectY;  /* List of Y-positions for bundle bounding path intersections */
    ArrList_pa     BBPIntersectZ;  /* List of Z-positions for bundle bounding path intersections */
    ArrList_pa     SeedX; /* List of X-positions of seed point on mid plane */
    ArrList_pa     SeedY; /* List of Y-positions of seed point on mid plane */
    ArrList_pa     SeedZ; /* List of Z-positions of seed point on mid plane */
    ArrList_pa     GradPaths;         /* List of pointers to GradPath structures */
    ArrList_pa     PathLength;
} MidPlnGradPath_s;

typedef struct _MidPlnGradPath_s  *MidPlnGradPath_pa;

Boolean_t         MidPlnGradPathIsValid(MidPlnGradPath_pa MidPlnGradPath);
MidPlnGradPath_pa MidPlnGradPathAlloc();
LgIndex_t         MidPlnGradPathGetGPCount(MidPlnGradPath_pa MidPlnGradPath);
LgIndex_t         MidPlnGradPathGetBBPICount(MidPlnGradPath_pa MidPlnGradPath);
void              MidPlnGradPathDealloc(MidPlnGradPath_pa *MidPlnGradPath);
void              MidPlnGradPathClear(MidPlnGradPath_pa MidPlnGradPath);
GradPath_pa       MidPlnGradPathGetGP(MidPlnGradPath_pa MidPlnGradPath,
                                      LgIndex_t         Offset);
XYZ_s             MidPlnGradPathGetSeedPos(MidPlnGradPath_pa MidPlnGradPath,
                                           LgIndex_t         Offset);
XYZ_s             MidPlnGradPathGetBBPIPos(MidPlnGradPath_pa MidPlnGradPath,
                                           LgIndex_t         Offset);
double            MidPlnGradPathGetPathLen(MidPlnGradPath_pa MidPlnGradPath,
                                           LgIndex_t         Offset);
Boolean_t         MidPlnGradPathSetBBPIPos(MidPlnGradPath_pa MidPlnGradPath,
                                           LgIndex_t         Offset,
                                           XYZ_s             Position);
Boolean_t         MidPlnGradPathAppendBBPIPos(MidPlnGradPath_pa MidPlnGradPath,
                                              XYZ_s             Position);
Boolean_t         MidPlnGradPathSetGP(MidPlnGradPath_pa MidPlnGradPath,
                                      LgIndex_t         Offset,
                                      XYZ_s             SeedPos,
                                      double            PathLength,
                                      GradPath_pa       GradPath);
Boolean_t         MidPlnGradPathInsertGP(MidPlnGradPath_pa MidPlnGradPath,
                                         LgIndex_t         Offset,
                                         XYZ_s             SeedPos,
                                         double            PathLength,
                                         GradPath_pa       GradPath);
Boolean_t         MidPlnGradPathRemoveGP(MidPlnGradPath_pa MidPlnGradPath,
                                         LgIndex_t   Offset);
Boolean_t         MidPlnGradPathAppendGP(MidPlnGradPath_pa MidPlnGradPath,
                                         XYZ_s             SeedPos,
                                         double            PathLength,
                                         GradPath_pa       GradPath);
LgIndex_t         MidPlnGradPathGetBegCPNum(MidPlnGradPath_pa MidPlnGradPath);
// Boolean_t MidPlnGradPathResample(MidPlnGradPath_pa MidPlnGradPath,
//                                 LgIndex_t         NumPointsNew);
Boolean_t         MidPlnGradPathGetGPIntersect(const MidPlnGradPath_pa MidPlnGradPath,
                                               const GradPath_pa       GradPath,
                                               double                 *X,
                                               double                 *Y,
                                               double                 *Z);
Boolean_t MidPlnGradPathAdd(const ZoneVarInfo_pa  ZoneVarInfo,
                            const LgIndex_t       BeginCrtPtNum,
                            const LgIndex_t       EndCrtPtNum,
                            const CritPoints_pa   CritPoints,
                            const Bundles_pa      Bundles,
                            const float           CPTolerance,
                            MidPlnGradPath_pa     MidPlnGradPath);
LgIndex_t MidPlnGradPathMinimizeLen(const ZoneVarInfo_pa  ZoneVarInfo,
                                    const LgIndex_t       BeginCrtPtNum,
                                    const LgIndex_t       EndCrtPtNum,
                                    const CritPoints_pa   CritPoints,
                                    const Bundles_pa      Bundles,
                                    const float           CPTolerance,
                                    MidPlnGradPath_pa     MidPlnGradPath);
//GradPath_pa MidPlnGradPath2CPPath(const ZoneVarInfo_pa  ZoneVarInfo,
//                                  const LgIndex_t       BeginCrtPtNum,
//                                  const LgIndex_t       EndCrtPtNum,
//                                  const CritPoints_pa   CritPoints,
//                                  const StreamDir_e     PathDir,
//                                  const float           CPTolerance,
//                                  MidPlnGradPath_pa     MidPlnGradPath);

#endif /* MIDPLANEGRADPATH_H_ */