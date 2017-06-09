/*
*****************************************************************
*****************************************************************
*******                                                  ********
*******    (C) Copyright 1988-2013 by TECPLOT INC.       *******
*******                                                  ********
*****************************************************************
*****************************************************************
*/


#ifndef GRADPATH_H_
#define GRADPATH_H_

// #define MAXGRADPATHPOINTS 1000
#define MAXCGPREFINEMENTS 8

/* For resampling of surface gradpaths, before seeding surface */
#define NUMRESAMPLEDPOINTSSURFACE 50
//#define brSurfaceIMax    7
#define brSurfaceIMax    50
//#define brSurfaceIMax  30
//#define brSurfaceIMax 60

//#define atomISegmentRatio 1.2 // maybe should be power(blah, 1/double(brSurfaceIMax-1) where blah empirically is ~12.9
#define atomISegmentRatio 1.0 //1.2
//#define atomKSegmentRatio 1.3 // maybe should be power(blah, 1/double((brSurfaceKMax-1)/2) where blah empirically is ~12.9
#define atomKSegmentRatio 1.3

//#define numCircularSeeds 9 // +1=Jmax
#define numCircularSeeds 18 // +1=Jmax
//#define numCircularSeeds 36 // +1=Jmax
//#define numCircularSeeds 72 // +1=Jmax
//#define bondBundleKMax 6
//#define bondBundleKMax 11
//#define bondBundleKMax 21
//#define bondBundleKMax 41

/* private ClientData structure */
/*
typedef struct _ZoneVarInfo_s
  {
    EntIndex_t ZoneNum;
    EntIndex_t UVarNum;
    EntIndex_t VVarNum;
    EntIndex_t WVarNum;
    EntIndex_t ChrgDensVarNum;
    EntIndex_t TypeVarNum;
    Boolean_t  PeriodicBC;

    // Temp var used for surface gradpath integration
    LgIndex_t LastElemUsed;
    SurfElemMap_pa SurfElemMap;
    Normals_pa     Normals;
    StreamDir_e    PathDir;
  } ZoneVarInfo_s;

typedef struct _ZoneVarInfo_s  *ZoneVarInfo_pa;
*/

/* private ClientData structure */
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


/* private ClientData structure */
typedef struct _CircleGradPath_s
{
    double      Center[3];
    double      PrincDir[3];
    double      Radius;
    Boolean_t   CompleteCircle;
    ArrList_pa  Angles;     /* List of doubles representing angles around circle */
    ArrList_pa  GradPaths;  /* List of pointers to GradPath structures */
    ArrList_pa  ConnectPathNums; /* List of GradPaths offsets for connectors */
} CircleGradPath_s;

typedef struct _CircleGradPath_s  *CircleGradPath_pa;


Boolean_t   GradPathIsValid(GradPath_pa GradPath);
GradPath_pa GradPathAlloc();
void        GradPathClear(GradPath_pa GradPath);
void        GradPathDealloc(GradPath_pa *GradPath);
LgIndex_t   GradPathGetCount(GradPath_pa GradPath);
Boolean_t   GradPathSetPoint(GradPath_pa GradPath,
                             LgIndex_t   PointOffset,
                             double      X,
                             double      Y,
                             double      Z,
                             double      Rho);
Boolean_t   GradPathInsert(GradPath_pa Target,
                           LgIndex_t   PointOffset,
                           GradPath_pa Source);
Boolean_t   GradPathInsertPoint(GradPath_pa GradPath,
                                LgIndex_t   PointOffset,
                                double      X,
                                double      Y,
                                double      Z,
                                double      Rho);
Boolean_t   GradPathRemovePoint(GradPath_pa GradPath,
                                LgIndex_t   PointOffset);
Boolean_t   GradPathAppendPoint(GradPath_pa GradPath,
                                double      X,
                                double      Y,
                                double      Z,
                                double      Rho);
Boolean_t   GradPathAppend(GradPath_pa Target,
                           GradPath_pa Source);
Boolean_t   GradPathGetPoint(GradPath_pa GradPath,
                             LgIndex_t   PointOffset,
                             double     *X,
                             double     *Y,
                             double     *Z,
                             double     *Rho);
GradPath_pa   GradPathGetFromTPZone(EntIndex_t ZoneNum);
EntIndex_t    GradPathTPZoneFromBegEndCP(LgIndex_t  BegCPNum,
                                         LgIndex_t  EndCPNum,
                                         EntIndex_t SourceZoneNum);
GradPath_pa   GradPathGetByBegEndCP(LgIndex_t     BegCPNum,
                                    LgIndex_t     EndCPNum);

double    GradPathGetLength(GradPath_pa GradPath);
double    GradPathGetDistToPt(GradPath_pa GradPath,
                              XYZ_s       Point);
double    GradPathGetSeparation(GradPath_pa GradPath1, GradPath_pa GradPath2);
Boolean_t GradPathResample(GradPath_pa *GradPath,
                           LgIndex_t    NumPointsNew,
                           double       segmentRatio = 1.0);
Boolean_t GradPathReverse(GradPath_pa *GradPath);

Boolean_t GradPathAdd(const ZoneVarInfo_pa  ZoneVarInfo,
                      const LgIndex_t       MaxNumPathPoints,
                      const StreamDir_e     PathDir,
                      const CritPoints_pa   CritPoints,
                      const double          CPTolerance,
                      LgIndex_t      *NumPathPoints,
                      GradPath_pa     GradPath);
Boolean_t GradPathAddSurf(const ZoneVarInfo_pa ZoneVarInfo, 
						  const LgIndex_t MaxNumPathPoints, 
						  const StreamDir_e PathDir, 
						  const CritPoints_pa CritPoints, 
						  const double CPToleranceOriginal, 
						  LgIndex_t *NumPathPoints,
						  GradPath_pa GradPath);
GradPath_pa GradPathAddMidway(const ZoneVarInfo_pa  ZoneVarInfo,
                              const LgIndex_t       MaxNumPathPoints,
                              const CritPoints_pa   CritPoints,
                              const XYZ_s           SeedPos,
                              const double          CPTolerance,
                              LgIndex_t      *NumPathPoints);
GradPath_pa GradPathSurfAddMidway(const ZoneVarInfo_pa  ZoneVarInfo,
                                  const LgIndex_t       MaxNumPathPoints,
                                  const CritPoints_pa   CritPoints,
                                  const XYZ_s           SeedPos,
                                  const double          CPTolerance,
                                  LgIndex_t            *NumPathPoints);
Boolean_t GradPathAdd2PtMinLen(const ZoneVarInfo_pa  ZoneVarInfo,
                               const LgIndex_t       NumPathPoints,
                               const CritPoints_pa   CritPoints,
                               GradPath_pa     GradPath);
Boolean_t GradPathTest();

Boolean_t         CircleGradPathIsValid(CircleGradPath_pa CircleGradPath);
CircleGradPath_pa CircleGradPathAlloc();
LgIndex_t         CircleGradPathGetCount(CircleGradPath_pa CircleGradPath);
LgIndex_t         CircleGradPathGetCPNCount(CircleGradPath_pa CircleGradPath);
void              CircleGradPathDealloc(CircleGradPath_pa *CircleGradPath);
void              CircleGradPathClear(CircleGradPath_pa CircleGradPath);
GradPath_pa       CircleGradPathGetGP(CircleGradPath_pa CircleGradPath,
                                      LgIndex_t         Index);
Boolean_t         CircleGradPathSetGP(CircleGradPath_pa CircleGradPath,
                                      LgIndex_t   Offset,
                                      double      Angle,
                                      GradPath_pa GradPath);
Boolean_t         CircleGradPathInsertGP(CircleGradPath_pa CircleGradPath,
                                         LgIndex_t   Offset,
                                         double      Angle,
                                         GradPath_pa GradPath);
Boolean_t         CircleGradPathRemoveGP(CircleGradPath_pa CircleGradPath,
                                         LgIndex_t   Offset);
Boolean_t         CircleGradPathAppendGP(CircleGradPath_pa CircleGradPath,
                                         double      Angle,
                                         GradPath_pa GradPath);
Boolean_t         CircleGradPathInsertCPN(CircleGradPath_pa CircleGradPath,
                                          LgIndex_t         Offset,
                                          LgIndex_t         ConnectPathNum);
Boolean_t         CircleGradPathAppendCPN(CircleGradPath_pa CircleGradPath,
                                          LgIndex_t         ConnectPathNum);
double            CircleGradPathGetAngle(CircleGradPath_pa CircleGradPath,
                                         LgIndex_t   Offset);
LgIndex_t         CircleGradPathGetCPN(CircleGradPath_pa CircleGradPath,
                                       LgIndex_t         Offset);
LgIndex_t         CircleGradPathGetEndCPNum(CircleGradPath_pa CircleGradPath,
                                            LgIndex_t         Offset);
Boolean_t CircleGradPathResample(CircleGradPath_pa CircleGradPath,
                                 LgIndex_t         NumPointsNew);
Boolean_t CircleGradPathAdd(const ZoneVarInfo_pa  ZoneVarInfo,
                            const LgIndex_t       BeginCrtPtNum,
                            const CritPoints_pa   CritPoints,
                            const StreamDir_e     PathDir,
                            const double          CPTolerance,
							const double          CircleRadius,
                            CircleGradPath_pa     CircleGradPath);
Boolean_t	GradPathSurfConnectBegEndToVolGPs(
							GradPath_pa			SurfGradPath,
							const CritPoints_pa	SurfCritPoints,
							const CritPoints_pa	VolCritPoints,
							const LgIndex_t		AtomNum);
EntIndex_t CreateConnectorZone(EntIndex_t    BaseZoneNum,
							   EntIndex_t    ChrgDensVarNum,
							   EntIndex_t    TypeVarNum,
							   CritPoints_pa CritPoints,
							   GradPath_pa   GradPath);

Boolean_t	GradPathProjectPointToSurf(const ZoneVarInfo_pa ZoneVarInfo,
									   const GradPath_pa	GradPath,
									   const LgIndex_t		CurrentGPPt,
									   double*				OriginalPtArray,
									   XYZ_s*				NewPt);

Boolean_t	GradPathProjectPointToSurf(const ZoneVarInfo_pa ZoneVarInfo,
									   double				CheckDistanceSqr,
									   XYZ_s				OldPt,
									   XYZ_s*				NewPt);

GradPath_pa GradPathVolGPFromAtomIsoTopoSurfCPClosest(const char		VolCPType1,
													  const char			VolCPType2,
													  const CritPoints_pa CritPoints,
													  const CritPoints_pa SurfCritPoints,
													  const LgIndex_t     AtomNum,
													  const char          SurfCPType,
													  const LgIndex_t     SurfCPOffset,
													  const double        Tolerance);

GradPath_pa CircleGradPathConcatenatePaths(CircleGradPath_pa CircleGradPath,
									 CritPoints_pa     CritPoints,
									 LgIndex_t         CGPNum,
									 LgIndex_t         EndCPNum,
									 LgIndex_t		   CGPWeight,
									 LgIndex_t         IMax);

Boolean_t GradPathSTSetVariables(const ZoneVarInfo_pa VolZoneVarInfo);

Boolean_t GradPathSTSetProperties(const LgIndex_t MaxSteps, 
								  const double StepSize, 
								  const double MinStepSize);

Boolean_t GradPathSTAddSurf(const ZoneVarInfo_pa VolZoneVarInfo, 
							const StreamDir_e PathDir, 
							const CritPoints_pa CritPoints, 
							LgIndex_t *NumPathPoints, 
							const double CPTolerance,
							GradPath_pa GradPath);

Boolean_t CoincidentWithSpecificCP(const XYZ_s TestPt,
								   const CritPoints_pa CritPoints, 
								   const double CPTolerance, 
								   char TestCrtPtType, 
								   LgIndex_t *CritPointNum);

Boolean_t GradPathCheckLenRatio(const GradPath_pa GradPath, 
								const CritPoints_pa CritPoints, 
								const double Radius);
#endif /* GRADPATH_H_ */
