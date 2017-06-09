#if !defined ADDON
#define ADDON
#endif /* ADDON */
#include "TECADDON.h"

#define TPINT_EXTERN
#include "zoneinfo.h"

static ZoneInfo_t ZoneInfo;

void SetCurZoneInfo(EntIndex_t Zone)
{
    FrameMode_e FrameMode   = TecUtilFrameGetMode();

    REQUIRE(TecUtilZoneIsEnabled(Zone));
    REQUIRE(FrameMode == Frame_XY || FrameMode == Frame_TwoD || FrameMode == Frame_ThreeD);

    ZoneInfo.Zone          = Zone;
    ZoneInfo.ZoneType      = TecUtilZoneGetType(Zone);

    if (FrameMode == Frame_ThreeD)
    {
        TecUtilZoneGetInfo(Zone, &ZoneInfo.IMax, &ZoneInfo.JMax, &ZoneInfo.KMax,
                           &ZoneInfo.XVar, &ZoneInfo.YVar, &ZoneInfo.ZVar, &ZoneInfo.NMap,
                           &ZoneInfo.UVar, &ZoneInfo.VVar, &ZoneInfo.WVar,
                           NULL, NULL, NULL);
    }
    else if (FrameMode == Frame_TwoD)
    {
        TecUtilZoneGetInfo(Zone, &ZoneInfo.IMax, &ZoneInfo.JMax, &ZoneInfo.KMax,
                           &ZoneInfo.XVar, &ZoneInfo.YVar, NULL, &ZoneInfo.NMap,
                           &ZoneInfo.UVar, &ZoneInfo.VVar, NULL,
                           NULL, NULL, NULL);
        ZoneInfo.ZVar = NULL;
        ZoneInfo.WVar = NULL;
    }
    else
    {
        EntIndex_t           MapZoneNum;
        EntIndex_t           XAxisVar;
        EntIndex_t           YAxisVar;
        SmInteger_t          XAxis;
        SmInteger_t          YAxis;
        FunctionDependency_e FunctionDependency;

        TecUtilZoneGetInfo(Zone, &ZoneInfo.IMax, &ZoneInfo.JMax, &ZoneInfo.KMax,
                           NULL, NULL, NULL, NULL,
                           NULL, NULL, NULL,
                           NULL, NULL, NULL);

        /* Get X for the first map. */
        TecUtilXYMapGetAssignment((EntIndex_t)1,
                                  &MapZoneNum,
                                  &XAxisVar,
                                  &YAxisVar,
                                  &XAxis,
                                  &YAxis,
                                  &FunctionDependency);

        ZoneInfo.XVar = TecUtilDataValueGetRef(Zone, XAxisVar);
        ZoneInfo.YVar = NULL;
        ZoneInfo.ZVar = NULL;
        if (ZoneInfo.ZoneType == ZoneType_FELineSeg)
            ZoneInfo.NMap = TecUtilDataNodeGetWritableRef(Zone);
        else
            ZoneInfo.NMap = NULL;
        ZoneInfo.UVar = NULL;
        ZoneInfo.VVar = NULL;
        ZoneInfo.WVar = NULL;
    }

    CurZoneInfo = &ZoneInfo;
}

