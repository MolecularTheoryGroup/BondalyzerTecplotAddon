/*
*****************************************************************
*****************************************************************
*******                                                  ********
*******    (C) Copyright 1988-2013 by TECPLOT INC.       *******
*******                                                  ********
*****************************************************************
*****************************************************************
*/


#ifndef BONDBUNDLE_H_
#define BONDBUNDLE_H_

double const ijkEpsilon = 0.00001;

typedef enum
{
    SurfType_BondRingCage,
    SurfType_PlsAtomRingCage,
    SurfType_MnsAtomRingCage,
    END_SurfType_e,
    SurfType_Invalid = BadEnumValue
} SurfType_e;


/* private ClientData structure */
typedef struct _BondBundle_s
{
    LgIndex_t      BondNum;

    LgIndex_t      PlsAtomNum;
    LgIndex_t      MnsAtomNum;

    EntIndex_t     BondPlsAtomLineZnNum;
    EntIndex_t     BondMnsAtomLineZnNum;

    ArrList_pa     BondRingCageSurfZnNums;
    ArrList_pa     PlsAtomRingCageSurfZnNums;
    ArrList_pa     MnsAtomRingCageSurfZnNums;

    EntIndex_t     BaseVolZoneNum;
    EntIndex_t     VolZoneNum;
} BondBundle_s;

typedef struct _BondBundle_s  *BondBundle_pa;

/* private ClientData structure */
typedef struct _BondBundles_s
{
    ArrList_pa     BondBundlesList;
} BondBundles_s;

typedef struct _BondBundles_s  *BondBundles_pa;



Boolean_t          BondBundleIsValid(const BondBundle_pa BondBundle);
BondBundle_pa      BondBundleAlloc();
LgIndex_t          BondBundleGetSurfCount(const BondBundle_pa BondBundle,
                                          const SurfType_e    SurfType);
void               BondBundleDealloc(BondBundle_pa *BondBundle);
void               BondBundleClear(BondBundle_pa BondBundle);
EntIndex_t         BondBundleGetSurfZnNum(const BondBundle_pa BondBundle,
                                          const SurfType_e    SurfType,
                                          const LgIndex_t     Offset);
Boolean_t          BondBundleSetSurfZnNum(BondBundle_pa     BondBundle,
                                          const SurfType_e  SurfType,
                                          const LgIndex_t   Offset,
                                          const EntIndex_t  SurfZnNum);
Boolean_t          BondBundleInsertSurfZnNum(BondBundle_pa     BondBundle,
                                             const SurfType_e  SurfType,
                                             const LgIndex_t   Offset,
                                             const EntIndex_t  SurfZnNum);
Boolean_t          BondBundleRemoveSurfZnNum(BondBundle_pa BondBundle,
                                             const SurfType_e    SurfType,
                                             const LgIndex_t     Offset);
Boolean_t          BondBundleAppendSurfZnNum(BondBundle_pa BondBundle,
                                             const SurfType_e    SurfType,
                                             const EntIndex_t    SurfZnNum);
Boolean_t      BondBundleGetSurfs(const LgIndex_t      BondCrtPtNum,
                                  const Bundles_pa     VolBundles,
                                  BondBundle_pa        BondBundle);
void BondBundleProcessSurfaceLineSegIntersections(BondBundle_pa const     bondBundle,
                                                  XYZ_s const &           lineSegStart,
                                                  XYZ_s const &           lineSegEnd,
                                                  IntersectionCallback_pf intersectionCallback,
                                                  void *                  clientData);
Boolean_t      BondBundleAdd(const ZoneVarInfo_pa  ZoneVarInfo,
                             const LgIndex_t       BondCrtPtNum,
                             const Bundles_pa      VolBundles,
                             BondBundle_pa         BondBundle);

Boolean_t          BondBundlesIsValid(const BondBundles_pa BondBundles);
BondBundles_pa     BondBundlesAlloc(const LgIndex_t ApproxNumBonds);
LgIndex_t          BondBundlesGetCount(const BondBundles_pa BondBundles);
void               BondBundlesDealloc(BondBundles_pa *BondBundles);
void               BondBundlesClear(BondBundles_pa BondBundles);
BondBundle_pa      BondBundlesGetBondBundle(const BondBundles_pa  BondBundles,
                                            const LgIndex_t       Offset);
Boolean_t          BondBundlesSetBondBundle(BondBundles_pa        BondBundles,
                                            const LgIndex_t       Offset,
                                            const BondBundle_pa   BondBundle);

Boolean_t      BondBundlesGetFromZoneAuxData(const EntIndex_t    BaseZoneNum,
                                             const CritPoints_pa VolCritPoints,
                                             BondBundles_pa      BondBundles);
//                                              ArrList_pa          BondBundlesList);
#endif /* BONDBUNDLE_H_ */
