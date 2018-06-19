/*
*****************************************************************
*****************************************************************
*******                                                  ********
*******    (C) Copyright 1988-2013 by TECPLOT INC.       *******
*******                                                  ********
*****************************************************************
*****************************************************************
*/


#ifndef CRITPOINT_H_
#define CRITPIONT_H_

#define MAXGRADPATHPOINTS 10000


/* private CritPoint structure */
typedef struct _CritPoints_s
{
    EntIndex_t SourceZoneNum;
    EntIndex_t Dimensions;
    LgIndex_t  MemUsageReported;
    LgIndex_t  NumCrtPts;
    LgIndex_t  NumCrtPtsM3;
    LgIndex_t  NumCrtPtsM1;
    LgIndex_t  NumCrtPtsP1;
    LgIndex_t  NumCrtPtsP3;
    /* The following two are far-field "Ring" and "Cage" CPs */
    LgIndex_t  NumFFCrtPtsP1;
    LgIndex_t  NumFFCrtPtsP3;
	double	   MinCPDistance;
    ArrList_pa X;
    ArrList_pa Y;
    ArrList_pa Z;
    ArrList_pa ChrgDens;
    /*
     * For Volume, Type is  Max     (Atom) = -3,
     *                      Saddle1 (Bond) = -1,
     *                      Saddle2 (Ring) =  1, and
     *                      Min     (Cage) =  3,
     *                      FFRing         = 11,
     *                      FFCage         = 13.
     * For Surface, Type is Max    = -2,
     *                      Saddle =  0,
     *                      Min    =  2.
     */
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

LgIndex_t IndexFromIJK(LgIndex_t    I,
                       LgIndex_t    J,
                       LgIndex_t    K,
                       LgIndex_t IMax,
                       LgIndex_t JMax,
                       LgIndex_t KMax,
                       Boolean_t PeriodicBC);

LgIndex_t PeriodicIndexFromIJK(LgIndex_t    I,
                               LgIndex_t    J,
                               LgIndex_t    K,
                               LgIndex_t IMax,
                               LgIndex_t JMax,
                               LgIndex_t KMax);

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

/*
Boolean_t ExtractCriticalPoints(const EntIndex_t    ZoneNum,
                                const EntIndex_t    UVarNum,
                                const EntIndex_t    VVarNum,
                                const EntIndex_t    WVarNum,
                                const EntIndex_t    ChrgDensVarNum,
                                const EntIndex_t    TypeVarNum,
                                const Boolean_t     PeriodicBC,
                                CritPoints_pa CritPoints);
                                */
Boolean_t ExtractCriticalPoints(const ZoneVarInfo_pa  VolZoneVarInfo,
                                const ZoneVarInfo_pa  SurfZoneVarInfo,
                                const Boolean_t       IsSurf,
                                CritPoints_pa         CritPoints);

EntIndex_t CritPointsCreateTPZone(CritPoints_pa CritPoints,
                                  EntIndex_t    ChrgDensVarNum,
                                  EntIndex_t    TypeVarNum,
                                  EntIndex_t    UVarNum,
                                  EntIndex_t    VVarNum,
								  EntIndex_t    WVarNum,
								  const Boolean_t IsSurfCPs = FALSE);

CritPoints_pa CritPointsGetFromTPZone(EntIndex_t TPZoneNum,
                                      EntIndex_t    ChrgDensVarNum,
                                      EntIndex_t    TypeVarNum,
                                      EntIndex_t    UVarNum,
                                      EntIndex_t    VVarNum,
                                      EntIndex_t    WVarNum);

Boolean_t CritPointsTest();

LgIndex_t VolCPFromAtomIsoTopoSurfCP(const char          VolCPType,
									 const CritPoints_pa CritPoints,
									 const CritPoints_pa SurfCritPoints,
									 const LgIndex_t     AtomNum,
									 const char          SurfCPType,
									 const LgIndex_t     SurfCPOffset,
									 const double        Tolerance);

LgIndex_t VolCPFromAtomIsoTopoSurfCPClosest(const char          VolCPType,
											const CritPoints_pa CritPoints,
											const CritPoints_pa SurfCritPoints,
											const LgIndex_t     AtomNum,
											const char          SurfCPType,
											const LgIndex_t     SurfCPOffset,
											const double        Tolerance);

Boolean_t	SurfCritPtsMoveToVolGPIsoSurfIntersection(
											const CritPoints_pa	SurfCritPoints,
											const CritPoints_pa	VolCritPoints,
											const LgIndex_t		AtomNum,
											const double		GPFoundTolerance);

Boolean_t CritPointGetTypeString(const LgIndex_t CritPointNum,
								 const char      CritPointType,
								 char           *CritPointTypeString);

LgIndex_t BondCPFromAtomIsoTopoMin(CritPoints_pa CritPoints,
								   CritPoints_pa SurfCritPoints,
								   LgIndex_t     AtomNum,
								   LgIndex_t     SurfMinNum,
								   double        Tolerance);

Boolean_t CritPointsMinDistance(CritPoints_pa CritPoints, 
								const Boolean_t TypesSpecified, 
								const char Type1, 
								const char Type2, 
								double *MinDist);

#endif /* CRITPOINT_H_ */
