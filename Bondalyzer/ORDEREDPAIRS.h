/*
*****************************************************************
*****************************************************************
*******                                                  ********
*******    (C) Copyright 1988-2013 by TECPLOT INC.       *******
*******                                                  ********
*****************************************************************
*****************************************************************
*/


#ifndef ORDEREDPAIRS_H_
#define ORDEREDPAIRS_H_

/* private OrderedPairs structure */
typedef struct _OrderedPairs_s
{
    ArrList_pa Item1;
    ArrList_pa Item2;
} OrderedPairs_s;

typedef struct _OrderedPairs_s  *OrderedPairs_pa;

Boolean_t        OrderedPairsIsValid(const OrderedPairs_pa OrderedPairs);
OrderedPairs_pa  OrderedPairsAlloc();
void             OrderedPairsClear(OrderedPairs_pa OrderedPairs);
Boolean_t        OrderedPairsRemovePoint(OrderedPairs_pa OrderedPairs,
                                         const LgIndex_t PointOffset);
void             OrderedPairsDealloc(OrderedPairs_pa *OrderedPairs);
LgIndex_t        OrderedPairsGetCount(const OrderedPairs_pa OrderedPairs);
Boolean_t        OrderedPairsSetAtOffset(OrderedPairs_pa OrderedPairs,
                                         const LgIndex_t PointOffset,
                                         const LgIndex_t Item1,
                                         const LgIndex_t Item2);
Boolean_t        OrderedPairsInsertAtOffset(OrderedPairs_pa OrderedPairs,
                                            const LgIndex_t PointOffset,
                                            const LgIndex_t Item1,
                                            const LgIndex_t Item2);
Boolean_t        OrderedPairsAppendAtEnd(OrderedPairs_pa OrderedPairs,
                                         const LgIndex_t Item1,
                                         const LgIndex_t Item2);
Boolean_t        OrderedPairsGetFromOffset(const OrderedPairs_pa  OrderedPairs,
                                           const LgIndex_t        PointOffset,
                                           LgIndex_t       *Item1,
                                           LgIndex_t       *Item2);
Boolean_t        OrderedPairsIsUniquePair(const OrderedPairs_pa  OrderedPairs,
                                          const LgIndex_t        Item1,
                                          const LgIndex_t        Item2);
LgIndex_t        OrderedPairsGetOffsetFromPair(const OrderedPairs_pa OrderedPairs,
                                               const LgIndex_t       Item1,
                                               const LgIndex_t       Item2);


#endif /* ORDEREDPAIRS_H_ */
