#pragma once

#ifdef ENGINE
#include "TecUtil.h"
#else
#include "TECADDON.h"
#endif
#include <string>

Boolean_t STDCALL LoaderCallback(StringList_pa Instructions);

const std::string LoaderName = "Value Load On Demand Loader";

struct ClientDataValues_s
{
    std::string FileName;     /* load/unload/cleanup */
    LgIndex_t   SeekOffset;   /* load/unload/cleanup */
    LgIndex_t   JMax;         /* get/set */
    LgIndex_t   IMax;         /* load/get/set */
    int         NumValues;    /* get/set */
    double      *StagingData; /* load/unload/cleanup/get/set */
};

/**
 * Create a test file of the given specifications, and return the filename
 * of the test file. Returns an empty string if the file is unable to
 * be created.
 */
std::string CreateTestFile(LgIndex_t IMax,
                           LgIndex_t JMax,
                           double    XStart,
                           double    YStart,
                           double    XDelta,
                           double    YDelta);


