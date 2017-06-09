/*
*****************************************************************
*****************************************************************
*******                                                  ********
*******    (C) Copyright 1988-2013 by TECPLOT INC.       *******
*******                                                  ********
*****************************************************************
*****************************************************************
*/


#ifndef ISOSURFGRADPATH_H_
#define ISOSURFGRADPATH_H_

#define ATOMISOSURFFACTOR 1.05 //1.02

/* private ClientData structure */
typedef struct _IsoSurfGradPath_s
{
	LgIndex_t      AtomNum;
	Bundles_pa     VolBundles;

	EntIndex_t     IsoTopoZone;

	LgIndex_t	   ElementCount;
	double		   AvgElemWidth;
	double		   Radius;
	EntIndex_t     SurfCPZoneNum;
	CritPoints_pa  SurfCritPoints;
	Bundles_pa     SurfBundles;

	ArrList_pa     GradPaths;         /* List of pointers to GradPath structures */
} IsoSurfGradPath_s;

typedef struct _IsoSurfGradPath_s  *IsoSurfGradPath_pa;

Boolean_t          IsoSurfGradPathIsValid(IsoSurfGradPath_pa IsoSurfGradPath);
IsoSurfGradPath_pa IsoSurfGradPathAlloc();
LgIndex_t          IsoSurfGradPathGetGPCount(IsoSurfGradPath_pa IsoSurfGradPath);
void               IsoSurfGradPathDealloc(IsoSurfGradPath_pa *IsoSurfGradPath);
void               IsoSurfGradPathClear(IsoSurfGradPath_pa IsoSurfGradPath);
GradPath_pa        IsoSurfGradPathGetGP(IsoSurfGradPath_pa IsoSurfGradPath,
										LgIndex_t         Offset);
Boolean_t          IsoSurfGradPathSetGP(IsoSurfGradPath_pa IsoSurfGradPath,
										LgIndex_t         Offset,
										GradPath_pa       GradPath);
Boolean_t          IsoSurfGradPathInsertGP(IsoSurfGradPath_pa IsoSurfGradPath,
										   LgIndex_t         Offset,
										   GradPath_pa       GradPath);
Boolean_t          IsoSurfGradPathRemoveGP(IsoSurfGradPath_pa IsoSurfGradPath,
										   LgIndex_t   Offset);
Boolean_t          IsoSurfGradPathAppendGP(IsoSurfGradPath_pa IsoSurfGradPath,
										   GradPath_pa       GradPath);
LgIndex_t          IsoSurfGradPathGetBegCPNum(IsoSurfGradPath_pa IsoSurfGradPath);
EntIndex_t         IsoSurfGradPathGetCPZoneNum(IsoSurfGradPath_pa IsoSurfGradPath);
Bundles_pa         IsoSurfGradPathGetSurfBundles(IsoSurfGradPath_pa IsoSurfGradPath);
CritPoints_pa      IsoSurfGradPathGetSurfCritPoints(IsoSurfGradPath_pa IsoSurfGradPath);
GradPath_pa        IsoSurfGradPathGetGPByBegEndCP(IsoSurfGradPath_pa IsoSurfGradPath,
												  LgIndex_t          BegCPOffset,
												  LgIndex_t          EndCPOffset);
Boolean_t IsoSurfGradPathAdd(const ZoneVarInfo_pa  VolZoneVarInfo,
							 const ZoneVarInfo_pa  SurfZoneVarInfo,
							 const LgIndex_t       AtomCrtPtNum,
							 const CritPoints_pa   VolCritPoints,
							 const Bundles_pa      VolBundles,
							 const float           CPTolerance,
							 IsoSurfGradPath_pa    IsoSurfGradPath);
Boolean_t IsoSurfGradPathSurfConnectSurfSaddleToVolRing(IsoSurfGradPath_pa IsoSurfGradPath, 
														const CritPoints_pa SurfCritPoints, 
														const CritPoints_pa VolCritPoints, 
														const LgIndex_t AtomNum);
#endif /* MIDPLANEGRADPATH_H_ */