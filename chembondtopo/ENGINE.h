#ifndef  ENGINE_H_
#define ENGINE_H_ /* Only include once */
#include "BUNDLES.h"

Boolean_t RecordMacroAndJournal(Strand_t StrandToAve,
                                EntIndex_t VarToAve);


Boolean_t ExtractTopoInZonesAuxData(Bundles_pa *   BundlesPtr,
                                    CritPoints_pa *CritPointsPtr);

ArrList_pa ListOfAtomsForBonds(Boolean_t  IsAll,
                               ArrList_pa BondList,
                               void      *BundlesVoid);
ArrList_pa ListOfCPOutsConnetedToCPIns(char       CPTypeIn,
                                       char       CPTypeOut,
                                       ArrList_pa CPInList,
                                       void      *BundlesVoid);
ArrList_pa ZoneListOfBondLinesForCritPoint(char          CritPointType,
                                           ArrList_pa    CritPointInList,
                                           CritPoints_pa CritPoints,
                                           void         *BundlesVoid);
ArrList_pa ZoneListOfBundleZonesForCritPoint(Boolean_t     GetEdges,
                                             Boolean_t     GetSurfaces,
                                             char          CPType,
                                             ArrList_pa    CPList,
                                             CritPoints_pa CritPoints,
                                             void         *BundlesVoid);
#endif /* ENGINE_H_ */
