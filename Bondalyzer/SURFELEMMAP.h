/*
*****************************************************************
*****************************************************************
*******                                                  ********
*******    (C) Copyright 1988-2013 by TECPLOT INC.       *******
*******                                                  ********
*****************************************************************
*****************************************************************
*/


#ifndef SURFELEMMAP_H_
#define SURFELEMMAP_H_

/* private SrfElemMap structure */
typedef struct _SurfElemMap_s
{
    EntIndex_t SurfZone;
    ArrList_pa ElemMapListOffset;  // Offset in ElemMapList of beg of ElemList for Node
    ArrList_pa ElemMapList;        // List of connected elements for all nodes in SurfZone
} SurfElemMap_s;

typedef struct _SurfElemMap_s  *SurfElemMap_pa;

Boolean_t      SurfElemMapIsValid(SurfElemMap_pa SurfElemMap);
void           SurfElemMapClear(SurfElemMap_pa SurfElemMap);
SurfElemMap_pa SurfElemMapAlloc(const LgIndex_t NumNodes);
Boolean_t      SurfElemMapAppendElemForNode(SurfElemMap_pa  SurfElemMap,
                                            LgIndex_t       NodeNum,
                                            LgIndex_t       ElemNum);
LgIndex_t      SurfElemMapGetElemCountForNode(const SurfElemMap_pa SurfElemMap,
                                              const LgIndex_t      NodeNum);
LgIndex_t      SurfElemMapGetElemNumForNodeOffset(const SurfElemMap_pa SurfElemMap,
                                                  const LgIndex_t      NodeNum,
                                                  const LgIndex_t      Offset);
void           SurfElemMapDealloc(SurfElemMap_pa *SurfElemMap);
Boolean_t      SurfElemMapCompute(SurfElemMap_pa SurfElemMap,
                                  EntIndex_t     SurfZone);
Boolean_t SurfElemMapTest();
#endif /* SURFELEMMAP_H_ */