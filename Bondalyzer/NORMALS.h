/*
*****************************************************************
*****************************************************************
*******                                                  ********
*******    (C) Copyright 1988-2013 by TECPLOT INC.       *******
*******                                                  ********
*****************************************************************
*****************************************************************
*/


#ifndef NORMALS_H_
#define NORMALS_H_

#include "SURFELEMMAP.h"

/* private Normals structure */
typedef struct _Normals_s
{
    EntIndex_t SurfZone;
    ArrList_pa Nx;   // X-component of normal vector for node
    ArrList_pa Ny;   // Y-component of normal vector for node
    ArrList_pa Nz;   // Z-component of normal vector for node
} Normals_s;

typedef struct _Normals_s  *Normals_pa;

Boolean_t      NormalsIsValid(Normals_pa Normals);
Normals_pa     NormalsAlloc(const LgIndex_t NumNodes);
void           NormalsClear(Normals_pa Normals);
Boolean_t      NormalsAppendForNode(Normals_pa  Normals,
                                    XYZ_s       Normal);
LgIndex_t      NormalsGetNodeCount(const Normals_pa Normals);
XYZ_s          NormalsGetNormalForNode(const Normals_pa Normals,
                                       const LgIndex_t  NodeNum);
void           NormalsDealloc(Normals_pa *Normals);
Boolean_t      NormalsCompute(Normals_pa Normals,
                              EntIndex_t SurfZone);
Boolean_t      NormalsComputeUsingSEM(Normals_pa     Normals,
                                      const EntIndex_t     SurfZone,
                                      const SurfElemMap_pa SurfElemMap);
// void * SurfElemMap);
Boolean_t NormalsComputeUsingIsoVarGrad(Normals_pa        Normals,
                                        const EntIndex_t  SurfZone,
                                        const EntIndex_t  GradXVar,
                                        const EntIndex_t  GradYVar,
                                        const EntIndex_t  GradZVar);
double    NormalsCompare(Normals_pa   Normals1,
                         Normals_pa   Normals2);

Boolean_t NormalsTest();
#endif /* NORMALS_H_ */