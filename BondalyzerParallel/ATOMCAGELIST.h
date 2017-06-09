/*
******************************************************************
******************************************************************
*******                                                   ********
******             (C) 1988-2010 Tecplot, Inc.             *******
*******                                                   ********
******************************************************************
******************************************************************
*/

#ifndef ATOMCAGELIST_H_
#define ATOMCAGELIST_H_

/* private AtomCageList structure */
typedef struct _AtomCageList_s
{
    EntIndex_t Atom;
    ArrList_pa Cages;
} AtomCageList_s;

typedef struct _AtomCageList_s  *AtomCageList_pa;

Boolean_t        AtomCageListIsValid(AtomCageList_pa AtomCageList);
AtomCageList_pa  AtomCageListAlloc();
void             AtomCageListClear(AtomCageList_pa AtomCageList);
Boolean_t        AtomCageListRemovePoint(AtomCageList_pa AtomCageList,
                                         LgIndex_t  PointOffset);
void             AtomCageListDealloc(AtomCageList_pa *AtomCageList);
LgIndex_t        AtomCageListGetCount(AtomCageList_pa AtomCageList);
Boolean_t        AtomCageListSetAtOffset(AtomCageList_pa AtomCageList,
                                         LgIndex_t  PointOffset,
                                         LgIndex_t  Atom1,
                                         LgIndex_t  Atom2);
Boolean_t        AtomCageListInsertAtOffset(AtomCageList_pa AtomCageList,
                                            LgIndex_t  PointOffset,
                                            LgIndex_t  Atom1,
                                            LgIndex_t  Atom2);
Boolean_t        AtomCageListAppendAtEnd(AtomCageList_pa AtomCageList,
                                         LgIndex_t  Atom1,
                                         LgIndex_t  Atom2);
Boolean_t        AtomCageListGetFromOffset(AtomCageList_pa  AtomCageList,
                                           LgIndex_t   PointOffset,
                                           LgIndex_t  *Atom1,
                                           LgIndex_t  *Atom2);



#endif /* ATOMCAGELISTAt_H_ */
