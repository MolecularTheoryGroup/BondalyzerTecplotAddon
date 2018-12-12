/*****************************************************************
 *****************************************************************
 *******                                                  ********
 *******  (C) Copyright 2007-2010 by Tecplot, Inc.         ********
 *******                                                  ********
 *****************************************************************
 *****************************************************************
*/

#ifndef BUNDLES_H_
#define BUNDLES_H_

#define MAXGRADPATHPOINTS 1000


/* private Bundles structure */
typedef struct _Bundles_s
{
    EntIndex_t SourceZoneNum;
    LgIndex_t  MemUsageReported;
    // Bundle Atom, Bond, Ring, and Cage. Index is bundle number.
    ArrList_pa Atom;
    ArrList_pa Bond;
    ArrList_pa Ring;
    ArrList_pa CageOrRing2;
    // Bundle may be irreducible or include FarField boundaries
    ArrList_pa CageOrRingCPType;
    ArrList_pa BundleType;
} Bundles_s;

typedef struct _Bundles_s  *Bundles_pa;

Boolean_t     BundlesIsValid(Bundles_pa const Bundles);
Bundles_pa    BundlesAlloc();
void          BundlesClear(Bundles_pa Bundles);
Boolean_t     BundlesRemovePoint(Bundles_pa Bundles,
                                 LgIndex_t  PointOffset);
void          BundlesDealloc(Bundles_pa *Bundles);
LgIndex_t     BundlesGetCount(Bundles_pa Bundles);
Boolean_t     BundlesSetAtOffset(Bundles_pa Bundles,
                                 LgIndex_t  PointOffset,
                                 LgIndex_t  Atom,
                                 LgIndex_t  Bond,
                                 LgIndex_t  Ring,
                                 LgIndex_t  CageOrRing2,
                                 char       CageOrRingCPType,
                                 char       BundleType);
Boolean_t   BundlesInsertAtOffset(Bundles_pa Bundles,
                                  LgIndex_t  PointOffset,
                                  LgIndex_t  Atom,
                                  LgIndex_t  Bond,
                                  LgIndex_t  Ring,
                                  LgIndex_t  CageOrRing2,
                                  char       CageOrRingCPType,
                                  char       BundleType);
Boolean_t   BundlesAppendAtEnd(Bundles_pa Bundles,
                               LgIndex_t  Atom,
                               LgIndex_t  Bond,
                               LgIndex_t  Ring,
                               LgIndex_t  CageOrRing2,
                               char       CageOrRingCPType,
                               char       BundleType);
Boolean_t BundlesGetFromOffset(Bundles_pa  const Bundles,
                               LgIndex_t   const PointOffset,
                               LgIndex_t  *Atom,
                               LgIndex_t  *Bond,
                               LgIndex_t  *Ring,
                               LgIndex_t  *CageOrRing2,
                               char       *CageOrRingCPType,
                               char       *BundleType);
Boolean_t BundlesGetBondRingFromAtomCage(Bundles_pa  const Bundles,
                                         LgIndex_t     const Atom,
                                         LgIndex_t     const Cage,
                                         Bundles_pa          BundlesAtomCage);
Boolean_t BundlesIsEdge(Bundles_pa    const Bundles,
                        LgIndex_t     const BundleNum,
                        char          const BegCrtPtType,
                        LgIndex_t     const BegCrtPtNum,
                        char          const EndCrtPtType,
                        LgIndex_t     const EndCrtPtNum);

Boolean_t BundlesIsSurface(Bundles_pa    const Bundles,
                           LgIndex_t     const BundleNum,
                           char          const BegCrtPtType,
                           LgIndex_t     const BegCrtPtNum,
                           char          const EndCrtPtType,
                           LgIndex_t     const EndCrtPtNum,
                           char          const ThrdCrtPtType,
                           LgIndex_t     const ThrdCrtPtNum);


#endif /* BUNDLES_H_ */
