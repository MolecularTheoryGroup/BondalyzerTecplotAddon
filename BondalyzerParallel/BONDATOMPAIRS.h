/*
******************************************************************
******************************************************************
*******                                                   ********
******             (C) 1988-2010 Tecplot, Inc.             *******
*******                                                   ********
******************************************************************
******************************************************************
*/
#ifndef BONDATOMPAIRS_H_
#define BONDATOMPAIRS_H_

/* private BondAtomPairs structure */
typedef struct _BondAtomPairs_s
{
    EntIndex_t SourceZoneNum;
    LgIndex_t  MemUsageReported;
    // Index of following arrays is bond number
    ArrList_pa Atom1;
    ArrList_pa Atom2;
} BondAtomPairs_s;

typedef struct _BondAtomPairs_s  *BondAtomPairs_pa;

Boolean_t        BondAtomPairsIsValid(BondAtomPairs_pa BondAtomPairs);
BondAtomPairs_pa BondAtomPairsAlloc();
void             BondAtomPairsClear(BondAtomPairs_pa BondAtomPairs);
Boolean_t        BondAtomPairsRemovePoint(BondAtomPairs_pa BondAtomPairs,
                                          LgIndex_t  PointOffset);
void             BondAtomPairsDealloc(BondAtomPairs_pa *BondAtomPairs);
LgIndex_t        BondAtomPairsGetCount(BondAtomPairs_pa BondAtomPairs);
Boolean_t        BondAtomPairsSetAtOffset(BondAtomPairs_pa BondAtomPairs,
                                          LgIndex_t  PointOffset,
                                          LgIndex_t  Atom1,
                                          LgIndex_t  Atom2);
Boolean_t        BondAtomPairsInsertAtOffset(BondAtomPairs_pa BondAtomPairs,
                                             LgIndex_t  PointOffset,
                                             LgIndex_t  Atom1,
                                             LgIndex_t  Atom2);
Boolean_t        BondAtomPairsAppendAtEnd(BondAtomPairs_pa BondAtomPairs,
                                          LgIndex_t  Atom1,
                                          LgIndex_t  Atom2);
Boolean_t        BondAtomPairsGetFromOffset(BondAtomPairs_pa  BondAtomPairs,
                                            LgIndex_t   PointOffset,
                                            LgIndex_t  *Atom1,
                                            LgIndex_t  *Atom2);



#endif /* BONDATOMPAIRS_H_ */
