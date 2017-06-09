/*****************************************************************
*****************************************************************
*******                                                  ********
*******  (C) Copyright 2007-2010 by Tecplot, Inc.         ********
*******                                                  ********
*****************************************************************
*****************************************************************
*/

#ifndef CRITPOINT_H_
#define CRITPIONT_H_

#define MAXGRADPATHPOINTS 1000


/* private CritPoint structure */
typedef struct _CritPoints_s
{
    EntIndex_t SourceZoneNum;
    LgIndex_t  MemUsageReported;
    LgIndex_t  NumCrtPts;
    LgIndex_t  NumCrtPtsM3;
    LgIndex_t  NumCrtPtsM1;
    LgIndex_t  NumCrtPtsP1;
    LgIndex_t  NumCrtPtsP3;
    ArrList_pa X;
    ArrList_pa Y;
    ArrList_pa Z;
    ArrList_pa ChrgDens;
    ArrList_pa Type;
    ArrList_pa PrincDirX;
    ArrList_pa PrincDirY;
    ArrList_pa PrincDirZ;
} CritPoints_s;

typedef struct _CritPoints_s  *CritPoints_pa;

Boolean_t BrickTrilinearWeight(double       r,
                               double       s,
                               double       t,
                               double       W[8]);


Boolean_t     CritPointsIsValid(CritPoints_pa CritPoints);
CritPoints_pa CritPointsAlloc();
void          CritPointsClear(CritPoints_pa CritPoints);
Boolean_t     CritPointsRemovePoint(CritPoints_pa CritPoints,
                                    LgIndex_t     PointOffset);
void          CritPointsDealloc(CritPoints_pa *CritPoints);
LgIndex_t     CritPointsGetCount(CritPoints_pa CritPoints);
Boolean_t     CritPointsSetPoint(CritPoints_pa CritPoints,
                                 LgIndex_t     PointOffset,
                                 double        X,
                                 double        Y,
                                 double        Z,
                                 double        ChrgDens,
                                 char          Type,
                                 double        PrincDirX,
                                 double        PrincDirY,
                                 double        PrincDirZ);
Boolean_t   CritPointsInsertPoint(CritPoints_pa CritPoints,
                                  LgIndex_t     PointOffset,
                                  double        X,
                                  double        Y,
                                  double        Z,
                                  double        ChrgDens,
                                  char          Type,
                                  double        PrincDirX,
                                  double        PrincDirY,
                                  double        PrincDirZ);
Boolean_t   CritPointsAppendPoint(CritPoints_pa CritPoints,
                                  double        X,
                                  double        Y,
                                  double        Z,
                                  double        ChrgDens,
                                  char          Type,
                                  double        PrincDirX,
                                  double        PrincDirY,
                                  double        PrincDirZ);
LgIndex_t CritPointsGetBegOffset(const CritPoints_pa CritPoints,
                                 const char          Type);
LgIndex_t CritPointsGetEndOffset(const CritPoints_pa CritPoints,
                                 const char          Type);
char      CritPointsGetType(const CritPoints_pa CritPoints,
                            const LgIndex_t     PointOffset);
Boolean_t CritPointsGetPoint(CritPoints_pa  CritPoints,
                             LgIndex_t      PointOffset,
                             double        *X,
                             double        *Y,
                             double        *Z,
                             double        *ChrgDens,
                             char          *Type,
                             double        *PrincDirX,
                             double        *PrincDirY,
                             double        *PrincDirZ);

Boolean_t ExtractCriticalPoints(const EntIndex_t    ZoneNum,
                                const EntIndex_t    UVarNum,
                                const EntIndex_t    VVarNum,
                                const EntIndex_t    WVarNum,
                                const EntIndex_t    ChrgDensVarNum,
                                const EntIndex_t    TypeVarNum,
                                CritPoints_pa CritPoints);

EntIndex_t CritPointsCreateTPZone(CritPoints_pa CritPoints,
                                  EntIndex_t    ChrgDensVarNum,
                                  EntIndex_t    TypeVarNum,
                                  EntIndex_t    UVarNum,
                                  EntIndex_t    VVarNum,
                                  EntIndex_t    WVarNum);

CritPoints_pa CritPointsGetFromTPZone(EntIndex_t TPZoneNum,
                                      EntIndex_t    ChrgDensVarNum,
                                      EntIndex_t    TypeVarNum,
                                      EntIndex_t    UVarNum,
                                      EntIndex_t    VVarNum,
                                      EntIndex_t    WVarNum);


#endif /* CRITPOINT_H_ */
