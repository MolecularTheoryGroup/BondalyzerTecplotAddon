/*
*****************************************************************
*****************************************************************
*******                                                  ********
*******    (C) Copyright 1988-2013 by TECPLOT INC.       *******
*******                                                  ********
*****************************************************************
*****************************************************************
*/

#ifndef  ENGINE_H_
#define ENGINE_H_ /* Only include once */

#include <fstream>

typedef enum
{
    CompletionLevel_CriticalPoints,
    CompletionLevel_BondAtomLines,
    CompletionLevel_RingCageLines,

    CompletionLevel_BRCSurfaces,
    CompletionLevel_RBASurfaces,
    CompletionLevel_AtomIsoSurfTopo,
    CompletionLevel_ARCSurfaces,
    CompletionLevel_BondBundle,
    CompletionLevel_All,

    CompletionLevel_Any,

    /* BEGINREMOVEFROMADDON */
    END_CompletionLevel_e,
    /* ENDREMOVEFROMADDON */
    CompletionLevel_Invalid = BadEnumValue
} CompletionLevel_e;

Boolean_t RecordMacroAndJournal(Strand_t StrandToAve,
                                EntIndex_t VarToAve);

int PrintTimeDate(std::ofstream & OutFile);


Boolean_t ExtractTopology(EntIndex_t        ZoneNum,
                          Set_pa            refinedGridSet,
                          EntIndex_t        UVarNum,
                          EntIndex_t        VVarNum,
                          EntIndex_t        WVarNum,
                          EntIndex_t        VelMagVarNum,
                          EntIndex_t        ChrgDensVarNum,
                          CompletionLevel_e CompletionLevel,
                          Boolean_t         PeriodicBC);
void CalcGradGradMag(Boolean_t IsPeriodic);
Boolean_t StatusUpdate(unsigned int CurrentNum, unsigned int TotalNum, const char* ProgresssText);


#endif /* ENGINE_H_ */
