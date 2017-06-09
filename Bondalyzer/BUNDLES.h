/*
*****************************************************************
*****************************************************************
*******                                                  ********
*******    (C) Copyright 1988-2013 by TECPLOT INC.       *******
*******                                                  ********
*****************************************************************
*****************************************************************
*/


#ifndef BUNDLES_H_
#define BUNDLES_H_

// #define MAXGRADPATHPOINTS 1000


/* private Bundles structure */
typedef struct _Bundles_s
{
    EntIndex_t SourceZoneNum;
    LgIndex_t  MemUsageReported;
    // Bundle Atom, Bond, Ring, and Cage. Index is bundle number.
    // Use one-based Atom, Bond, Ring, and Cage numbers, so that negative
    // numbers can be used for far-field Rings and Cages
    ArrList_pa Atom;
    ArrList_pa Bond;
    ArrList_pa Ring;
    ArrList_pa Cage;
} Bundles_s;

typedef struct _Bundles_s  *Bundles_pa;

Boolean_t     BundlesIsValid(Bundles_pa Bundles);
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
                                 LgIndex_t  Cage);
Boolean_t   BundlesInsertAtOffset(Bundles_pa Bundles,
                                  LgIndex_t  PointOffset,
                                  LgIndex_t  Atom,
                                  LgIndex_t  Bond,
                                  LgIndex_t  Ring,
                                  LgIndex_t  Cage);
Boolean_t   BundlesAppendAtEnd(Bundles_pa Bundles,
                               LgIndex_t  Atom,
                               LgIndex_t  Bond,
                               LgIndex_t  Ring,
                               LgIndex_t  Cage);
Boolean_t BundlesGetFromOffset(const Bundles_pa  Bundles,
                               const LgIndex_t   PointOffset,
                               LgIndex_t  *Atom,
                               LgIndex_t  *Bond,
                               LgIndex_t  *Ring,
                               LgIndex_t  *Cage);
Boolean_t BundlesGetBondRingFromAtomCage(const Bundles_pa  Bundles,
                                         const LgIndex_t     Atom,
                                         const LgIndex_t     Cage,
                                         Bundles_pa          BundlesAtomCage);

Boolean_t  CreateBundlesForBond(EntIndex_t        BaseZoneNum,
								CritPoints_pa     CritPoints,
								CircleGradPath_pa CircleGradPath,
								LgIndex_t         CGPNumBeg,
								LgIndex_t         CGPNumEnd,
								Bundles_pa        Bundles);


ArrList_pa  CreateAtomCageList(LgIndex_t        Atom,
							   Bundles_pa       Bundles,
							   CritPoints_pa    CritPoints);

ArrList_pa  CreateAtomBondCageList(LgIndex_t        Atom,
								   LgIndex_t        Bond,
								   Bundles_pa       Bundles,
								   CritPoints_pa    CritPoints);

ArrList_pa  CreateRingBondListForAtomCage(LgIndex_t        Atom,
										  LgIndex_t        Cage,
										  Bundles_pa       Bundles,
										  CritPoints_pa    CritPoints);



#endif /* BUNDLES_H_ */
