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

Boolean_t     BundlesIsValid(const Bundles_pa Bundles);
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
Boolean_t BundlesGetFromOffset(const Bundles_pa  Bundles,
                               const LgIndex_t   PointOffset,
                               LgIndex_t  *Atom,
                               LgIndex_t  *Bond,
                               LgIndex_t  *Ring,
                               LgIndex_t  *CageOrRing2,
                               char       *CageOrRingCPType,
                               char       *BundleType);
Boolean_t BundlesGetBondRingFromAtomCage(const Bundles_pa  Bundles,
                                         const LgIndex_t     Atom,
                                         const LgIndex_t     Cage,
                                         Bundles_pa          BundlesAtomCage);
Boolean_t BundlesIsEdge(const Bundles_pa    Bundles,
                        const LgIndex_t     BundleNum,
                        const char          BegCrtPtType,
                        const LgIndex_t     BegCrtPtNum,
                        const char          EndCrtPtType,
                        const LgIndex_t     EndCrtPtNum);

Boolean_t BundlesIsSurface(const Bundles_pa    Bundles,
                           const LgIndex_t     BundleNum,
                           const char          BegCrtPtType,
                           const LgIndex_t     BegCrtPtNum,
                           const char          EndCrtPtType,
                           const LgIndex_t     EndCrtPtNum,
                           const char          ThrdCrtPtType,
                           const LgIndex_t     ThrdCrtPtNum);


#endif /* BUNDLES_H_ */
