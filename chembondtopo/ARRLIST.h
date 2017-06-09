#if !defined ARRLIST_h
#define ARRLIST_h

typedef enum
{
    ArrListType_UnsignedChar,
    ArrListType_UnsignedShort,
    ArrListType_UnsignedInt,
    ArrListType_UnsignedLong,

    ArrListType_Char,
    ArrListType_Short,
    ArrListType_Int,
    ArrListType_Long,
    ArrListType_Double,

    ArrListType_UnsignedCharPtr,
    ArrListType_UnsignedShortPtr,
    ArrListType_UnsignedIntPtr,
    ArrListType_UnsignedLongPtr,

    ArrListType_CharPtr,
    ArrListType_ShortPtr,
    ArrListType_IntPtr,
    ArrListType_LongPtr,
    ArrListType_DoublePtr,

    ArrListType_VoidPtr,
    ArrListType_FunctionPtr,

    ArrListType_Any,

    /* BEGINREMOVEFROMADDON */
    END_ArrListType_e,
    /* ENDREMOVEFROMADDON */
    ArrListType_Invalid = BadEnumValue
} ArrListType_e;

/* public handles to list and list item */
typedef struct _ArrList_s  *ArrList_pa;


/* union of intrinsic types for sizing */
typedef union
{
    unsigned char  UnsignedChar;
    unsigned short UnsignedShort;
    unsigned int   UnsignedInt;
    unsigned long  UnsignedLong;

    char           Char;
    short          Short;
    int            Int;
    long           Long;
    double         Double;

    unsigned char  *UnsignedCharPtr;
    unsigned short *UnsignedShortPtr;
    unsigned int   *UnsignedIntPtr;
    unsigned long  *UnsignedLongPtr;

    char           *CharPtr;
    short          *ShortPtr;
    int            *IntPtr;
    long           *LongPtr;
    double         *DoublePtr;

    void           *VoidPtr;
    void (*FunctionPtr)(void);
} ArrListItem_u;

extern Boolean_t ArrListIsValid(ArrList_pa ArrList);
extern ArrListType_e ArrListGetType(ArrList_pa ArrList);
extern ArrList_pa ArrListAlloc(LgIndex_t       EstimatedCapacity,
                               ArrListType_e Type);
extern void ArrListDealloc(ArrList_pa *ArrList);
extern LgIndex_t ArrListGetCount(ArrList_pa ArrList);
extern void ArrListClear(ArrList_pa ArrList);
extern void ArrListClearItems(ArrList_pa ArrList,
                              LgIndex_t    ItemOffset,
                              LgIndex_t    Count);
extern void ArrListClearItem(ArrList_pa ArrList,
                             LgIndex_t    ItemOffset);
extern ArrList_pa ArrListRemoveItems(ArrList_pa ArrList,
                                     LgIndex_t    ItemOffset,
                                     LgIndex_t    Count);
extern ArrListItem_u ArrListRemoveItem(ArrList_pa ArrList,
                                       LgIndex_t    ItemOffset);
extern Boolean_t ArrListInsertItem(ArrList_pa    ArrList,
                                   LgIndex_t       ItemOffset,
                                   ArrListItem_u Item);
extern Boolean_t ArrListInsert(ArrList_pa Target,
                               LgIndex_t    ItemOffset,
                               ArrList_pa Source);
extern ArrList_pa ArrListGetItems(ArrList_pa ArrList,
                                  LgIndex_t    ItemOffset,
                                  LgIndex_t    Count);
extern ArrListItem_u ArrListGetItem(ArrList_pa ArrList,
                                    LgIndex_t    ItemOffset);
extern Boolean_t ArrListSetItem(ArrList_pa    ArrList,
                                LgIndex_t   ItemOffset,
                                ArrListItem_u Item);
extern Boolean_t ArrListAppendItem(ArrList_pa    ArrList,
                                   ArrListItem_u Item);
extern Boolean_t ArrListAppend(ArrList_pa Target,
                               ArrList_pa Source);
extern ArrList_pa ArrListCopy(ArrList_pa ArrList);
extern Boolean_t ArrListToNative(ArrList_pa ArrList,
                                 void         **Target,
                                 LgIndex_t    *Count);
extern ArrList_pa ArrListFromNative(void            *Source,
                                    LgIndex_t       Count,
                                    ArrListType_e Type);

#endif /* ARRLIST_h */
