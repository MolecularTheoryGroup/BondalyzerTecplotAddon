/*
*****************************************************************
*****************************************************************
*******                                                  ********
*******    (C) Copyright 1988-2013 by TECPLOT INC.       *******
*******                                                  ********
*****************************************************************
*****************************************************************
*/


#ifndef ARRLISTTOOLS_H_
#define ARRLISTTOOLS_H_

Boolean_t ArrListIsUniqueLongItem(const ArrList_pa    ArrList,
                                  const LgIndex_t     LongArrListItem);

LgIndex_t ArrListAppendUniqueLongItem(ArrList_pa    ArrList,
                                      LgIndex_t     LongArrListItem);

#endif /* ARRLISTTOOLS_H_ */
