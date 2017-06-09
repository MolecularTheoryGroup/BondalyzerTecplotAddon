#if !defined __ZONEINFO_H
#define __ZONEINFO_H

typedef struct ZoneInfo_s
{
    EntIndex_t   Zone;
    ZoneType_e   ZoneType;
    LgIndex_t    IMax;
    LgIndex_t    JMax;
    LgIndex_t    KMax;
    FieldData_pa XVar;
    FieldData_pa YVar;
    FieldData_pa ZVar;
    NodeMap_pa   NMap;
    FieldData_pa UVar;
    FieldData_pa VVar;
    FieldData_pa WVar;
} ZoneInfo_t;


#ifndef TPINT_EXTERN
extern ZoneInfo_t *CurZoneInfo;
#else
ZoneInfo_t *CurZoneInfo = NULL;
#endif

extern void SetCurZoneInfo(EntIndex_t Zone);

#endif /* __ZONEINFO_H */
