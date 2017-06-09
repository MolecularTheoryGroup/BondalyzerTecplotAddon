/*****************************************************************
 *****************************************************************
 *******                                                  ********
 *******       (C) Copyright 2007  by Tecplot Inc.        ********
 *******                                                  ********
 *****************************************************************
 *****************************************************************
*/

#ifndef ISOTOPO_H_
#define ISOTOPO_H_

/* private IsoTopo structure */
typedef struct _IsoTopo_s
{
    // IsoTopo Inputs: Coordinates of nearby point used to isolate desired
    // topological component of isosurface.
    double XNear;
    double YNear;
    double ZNear;
    // Temporary
    EntIndex_t IsoZone;
    ArrList_pa IsoElemNums;  // TODO: figure out this confusion
    // IsoTopo Outputs: Fit coefficients.
    EntIndex_t IsoTopoZone;
} IsoTopo_s;

typedef struct _IsoTopo_s  *IsoTopo_pa;

Boolean_t     IsoTopoIsValid(IsoTopo_pa IsoTopo);
IsoTopo_pa    IsoTopoAlloc();
Boolean_t     IsoTopoAppendElemAtEnd(IsoTopo_pa  IsoTopo,
                                     LgIndex_t   IsoElemNum);
LgIndex_t     IsoTopoGetElemCount(const IsoTopo_pa IsoTopo);
void          IsoTopoDealloc(IsoTopo_pa *IsoTopo);
Boolean_t     IsoTopoTestPointTriDist();
//Boolean_t IsoTopoGetCoefficients(IsoTopo_pa  IsoTopo,
//                                    double       **FitCoeffs);
EntIndex_t IsoTopoCompute(IsoTopo_pa IsoTopo,
                          EntIndex_t IsoSurfZone,
                          double X,
                          double Y,
                          double Z);
#endif /* IsoTopo_H_ */