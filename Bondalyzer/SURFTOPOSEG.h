/*
*****************************************************************
*****************************************************************
*******                                                  ********
*******    (C) Copyright 1988-2013 by TECPLOT INC.       *******
*******                                                  ********
*****************************************************************
*****************************************************************
*/


#ifndef SURFTOPOSEG_H_
#define SURFTOPOSEG_H_

/* private SurfTopoSeg structure */
typedef struct _SurfTopoSeg_s
{
    // SurfTopoSeg Inputs: Coordinates of nearby point used to isolate desired
    // topological component of isosurface.
    // Temporary
    LgIndex_t  FirstElemNumOffset;
    EntIndex_t SurfZone;
    ArrList_pa SurfElemsRemaining;
    ArrList_pa SurfSegZoneNums;
} SurfTopoSeg_s;

typedef struct _SurfTopoSeg_s  *SurfTopoSeg_pa;

Boolean_t       SurfTopoSegIsValid(SurfTopoSeg_pa SurfTopoSeg);
SurfTopoSeg_pa  SurfTopoSegAlloc();
Boolean_t       SurfTopoSegAppendElemAtEnd(SurfTopoSeg_pa  SurfTopoSeg,
                                           LgIndex_t       ElemNum);
LgIndex_t       SurfTopoSegGetElemCount(const SurfTopoSeg_pa SurfTopoSeg);
Boolean_t       SurfTopoSegAppendElem(SurfTopoSeg_pa  SurfTopoSeg,
                                      LgIndex_t       ElemNum);
LgIndex_t       SurfTopoSegRemoveElem(SurfTopoSeg_pa  SurfTopoSeg,
                                      LgIndex_t       Offset);
LgIndex_t       SurfTopoSegGetSegZoneCount(const SurfTopoSeg_pa SurfTopoSeg);
Boolean_t       SurfTopoSegAppendSegZone(SurfTopoSeg_pa  SurfTopoSeg,
                                         EntIndex_t      ZoneNum);

void          SurfTopoSegDealloc(SurfTopoSeg_pa *SurfTopoSeg);
Boolean_t     SurfTopoSegInitialize(SurfTopoSeg_pa SurfTopoSeg,
                                    EntIndex_t     SurfZone);
EntIndex_t    SurfTopoSegComputeSeg(SurfTopoSeg_pa SurfTopoSeg,
                                    LgIndex_t      FirstElemNumOffset);
Boolean_t     SurfTopoSegTest();

Boolean_t SurfTopoDimensions(EntIndex_t SurfZone, 
							 const XYZ_s Center, 
							 LgIndex_t *ElemCount,
							 double *SurfRadius, 
							 double *AvgMaxExtent);

ZoneType_e SurfTopoGetElemNodeXYZByElemNum(EntIndex_t SurfZone, 
										   LgIndex_t ElemNum, 
										   XYZ_s* Nodes);
#endif /* SurfTopoSeg_H_ */
