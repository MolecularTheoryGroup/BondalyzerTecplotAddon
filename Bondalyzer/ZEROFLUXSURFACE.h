/*****************************************************************
 *****************************************************************
 *******                                                  ********
 *******     (C) Copyright 2013  by Tecplot Inc.          ********
 *******                                                  ********
 *****************************************************************
 *****************************************************************
*/

#ifndef ZEROFLUXSURFACE_H_
#define ZEROFLUXSURFACE_H_

EntIndex_t ZFSurfZoneNumFromAtomRingCage(const EntIndex_t BaseZoneNum,
										 const LgIndex_t  AtomNum,
										 const LgIndex_t  RingNum,
										 const LgIndex_t  CageNum);

EntIndex_t CreateBRSurfaceZone(EntIndex_t        BaseZoneNum,
							   EntIndex_t        ChrgDensVarNum,
							   EntIndex_t        TypeVarNum,
							   CritPoints_pa     CritPoints,
							   CircleGradPath_pa CircleGradPath,
							   LgIndex_t         IMax,
							   LgIndex_t         CGPNumBeg,
							   LgIndex_t         CGPNumEnd);

EntIndex_t CreateZFSurfZoneFrom3CPsGradpathList(EntIndex_t        ChrgDensVarNum,
												EntIndex_t        TypeVarNum,
												CritPoints_pa     CritPoints,
												LgIndex_t			BegCrtPtNum,
												LgIndex_t			FirstCrtPtNum,
												LgIndex_t			LastCrtPtNum,
												ArrList_pa        GradPathList,
												LgIndex_t         IMax);

EntIndex_t CreateZFSurfZoneFromGPList(EntIndex_t        ChrgDensVarNum,
									  EntIndex_t        TypeVarNum,
									  CritPoints_pa     CritPoints,
									  ArrList_pa        GradPathList,
									  LgIndex_t         IMax);

#endif