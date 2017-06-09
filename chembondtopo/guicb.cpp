#include "TECADDON.h"
#include "ADDGLBL.h"
#if !defined (MSWIN)
#include <unistd.h>
#endif
#include "GUIDEFS.h"
#include "ADKUTIL.h"
#include "ALLOC.h"
#include "ARRLIST.h"
#include "CRITPOINTS.h"
// #include "BUNDLES.h"
#include "ENGINE.h"
#include "MAIN.h"

//#include "toolbox/tptoolbox.h"
//#include "toolbox/ArgList.h"

CritPoints_pa CritPoints;
Bundles_pa Bundles;

//
// Sidebar version
//

typedef enum
{
    DisplayMode_AllAtomsAndBonds,
    DisplayMode_SingleAtom,
    DisplayMode_SingleBaderAtom,
    DisplayMode_SingleBond,
    DisplayMode_SingleBondBundle,
    END_DisplayMode_e,
    DisplayMode_Invalid = -1
} DisplayMode_e;

static DisplayMode_e curDisplayMode = DisplayMode_Invalid;
static LgIndex_t curAtom = -1;
static LgIndex_t curBond = -1;

static void updateInfoField(void)
{
    char message[256];
    sprintf_s(message, "");
    if ( curAtom >= 1 )
    {
        CHECK(curBond == -1);
        CHECK(curAtom <= CritPoints->NumCrtPtsM3);
        char *atomType = "unknown";
        double charge = 0.0;
        double volume = 0.0;
        sprintf_s(message, "atom %d\n%s\nCharge=%lf\nVolume=%lf", curAtom, atomType, charge, volume);
    }
    else if ( curBond >= 1 )
    {
        CHECK(curAtom == -1);
        CHECK(curBond <= CritPoints->NumCrtPtsM1);
        char *firstAtomType = "??";
        char *secondAtomType = "??";
        double charge = 0.0;
        double volume = 0.0;
        sprintf_s(message, "bond %d\n%s-%s\nCharge=%lf\nVolume=%lf", curBond, firstAtomType, secondAtomType, charge, volume);
    }
    TecGUITextSetString(Info_T_S1, message);
}


DisplayMode_e getDisplayModeFromLgIndex(LgIndex_t Mode)
{
    DisplayMode_e displayMode = DisplayMode_Invalid;
    switch (Mode)
    {
        case 1: displayMode = DisplayMode_AllAtomsAndBonds; break;
        case 2: displayMode = DisplayMode_SingleAtom; break;
        case 3: displayMode = DisplayMode_SingleBond; break;
        case 4: displayMode = DisplayMode_SingleBondBundle; break;
        default: CHECK(FALSE);
    }
    ENSURE(VALID_ENUM(displayMode, DisplayMode_e));
    return displayMode;
}


void updateCurAtomAndBond(LgIndex_t atom,
                          LgIndex_t bond)
{
    LgIndex_t numAtoms = CritPoints->NumCrtPtsM3;
    LgIndex_t numBonds = CritPoints->NumCrtPtsM1;
    REQUIRE(curAtom == -1 || (1<=atom && atom<=numAtoms));
    REQUIRE(curBond == -1 || (1<=bond && bond<=numBonds));
    REQUIRE(IMPLICATION(curAtom != -1, curBond == -1) && IMPLICATION(curBond != -1, curAtom== -1));
    curAtom = atom;
    curBond = bond;
    TecGUIListDeselectAllItems(AtomOrBond_MLST_S1);
    if ( curAtom != -1 )
    {
        TecGUIListSetSelectedItem(AtomOrBond_MLST_S1, curAtom);
    }
    else if ( curBond != -1 )
    {
        if ( curDisplayMode == DisplayMode_AllAtomsAndBonds )
            TecGUIListSetSelectedItem(AtomOrBond_MLST_S1, curBond+CritPoints->NumCrtPtsM3);
        else
            TecGUIListSetSelectedItem(AtomOrBond_MLST_S1, curBond);
    }
}

void fillGuiListWithAtoms(LgIndex_t guiList)
{
    REQUIRE("Valid guiList");
    REQUIRE(VALID_REF(CritPoints));
    LgIndex_t numAtoms = CritPoints->NumCrtPtsM3;
    for ( LgIndex_t atom = 1; atom <= numAtoms; atom++ )
    {
        char string[255];
        sprintf_s(string, "atom %d", atom);
        TecGUIListAppendItem(guiList, string);
    }
}

void fillGuiListWithBonds(LgIndex_t guiList)
{
    REQUIRE("Valid guiList");
    LgIndex_t numBonds = CritPoints->NumCrtPtsM1;
    for ( LgIndex_t bond = 1; bond <= numBonds; bond++ )
    {
        char string[255];
        sprintf_s(string, "bond %d", bond);
        TecGUIListAppendItem(guiList, string);
    }
}

void updateBondOrAtomList(DisplayMode_e displayMode)
{
    REQUIRE(VALID_ENUM(displayMode, DisplayMode_e));

    TecGUIListDeleteAllItems(AtomOrBond_MLST_S1);
    switch (displayMode)
    {
        case DisplayMode_AllAtomsAndBonds :
            {
                TecGUILabelSetText(AtomOrBond_LBL_S1, "Atoms & Bonds");
                fillGuiListWithAtoms(AtomOrBond_MLST_S1);
                fillGuiListWithBonds(AtomOrBond_MLST_S1);
            } break;
        case DisplayMode_SingleAtom :
        case DisplayMode_SingleBaderAtom :
            {
                TecGUILabelSetText(AtomOrBond_LBL_S1, "Atoms");
                fillGuiListWithAtoms(AtomOrBond_MLST_S1);
            } break;
        case DisplayMode_SingleBond :
        case DisplayMode_SingleBondBundle :
            {
                TecGUILabelSetText(AtomOrBond_LBL_S1, "Bonds");
                fillGuiListWithBonds(AtomOrBond_MLST_S1);
            } break;
        default:
            CHECK(FALSE);
    }
    updateCurAtomAndBond(curAtom, curBond);
}

void setDisplayMode(DisplayMode_e displayMode);

Set_pa  orgActiveZoneSet = NULL;
typedef struct
{
    double       meshThickness;
    ColorIndex_t meshColor;
    ColorIndex_t scatterColor;
} OrgZoneProperties_s;

OrgZoneProperties_s *orgZoneProperties = NULL;

EntIndex_t* atomZoneNumbers = NULL;

void atomZoneInfoDealloc(EntIndex_t **atomZoneNumbers)
{
    REQUIRE(VALID_REF(atomZoneNumbers) && VALID_REF_OR_NULL(*atomZoneNumbers));
    if ( *atomZoneNumbers != NULL )
    {
        FREE_ARRAY(*atomZoneNumbers, "atomZoneNumbers");
        *atomZoneNumbers = NULL;
    }
}

/*
 * Creates individual zones for each atom if they do not already exist
 */
void atomZoneInfoAlloc(EntIndex_t **atomZoneNumbers)
{
    REQUIRE(VALID_REF(atomZoneNumbers) && VALID_REF_OR_NULL(*atomZoneNumbers));

    EntIndex_t numAtoms = CritPoints->NumCrtPtsM3;
    if ( *atomZoneNumbers != NULL )
        atomZoneInfoDealloc(atomZoneNumbers);
    *atomZoneNumbers = ALLOC_ARRAY(numAtoms, EntIndex_t, "*atomZoneNumbers");
    CHECK(VALID_REF(*atomZoneNumbers));

    EntIndex_t numAtomZones = 0;
    EntIndex_t critPointZone = -1;
    EntIndex_t numZones = TecUtilDataSetGetNumZones();
    for (EntIndex_t zone = 1; zone <= numZones; zone++)
    {
        /* Set Zone aux data to completely transfer necessary information to Tecplot */
        AuxData_pa zoneAuxDataRef = TecUtilAuxDataZoneGetRef(zone);
        if (zoneAuxDataRef != NULL)
        {
            char* zoneTypeStr = NULL;
            Boolean_t retain; //
            if ( TecUtilAuxDataGetStrItemByName(zoneAuxDataRef,
                                                "CompChem.ZoneType",
                                                &zoneTypeStr,
                                                &retain) )
            {
                if (zoneTypeStr != NULL)
                {
                    if (strcmp(zoneTypeStr, "CriticalPoints") == 0)
                    {
                        // there are critical point zones for iso-surfaces around atoms as well, so check the base zone
                        char* baseZoneStr = NULL;
                        if ( TecUtilAuxDataGetStrItemByName(zoneAuxDataRef,
                                                            "CompChem.BaseZoneNum",
                                                            &baseZoneStr,
                                                            &retain) )
                        {
                            if ( strcmp(baseZoneStr, "1") == 0 )
                            {
                                CHECK(critPointZone==-1);
                                critPointZone = zone;
                            }
                        }
                        TecUtilStringDealloc(&baseZoneStr);
                    }
                    else if (strcmp(zoneTypeStr, "Atom") == 0)
                    {
                        (*atomZoneNumbers)[numAtomZones] = zone;
                        numAtomZones++;
                    }
                    TecUtilStringDealloc(&zoneTypeStr);
                }
            }
        }
    }
    CHECK(critPointZone != -1);
    CHECK(numAtomZones == 0 || numAtomZones == numAtoms); // created or not
    if ( numAtomZones == 0 ) // need to create
    {
        Set_pa activeZoneSet = NULL;
        TecUtilZoneGetActive(&activeZoneSet);
        CHECK(VALID_REF(activeZoneSet));
        Set_pa atomsZoneSet = TecUtilSetAlloc(TRUE);
        CHECK(VALID_REF(atomsZoneSet));
        for (EntIndex_t atom = 1; atom <= numAtoms; atom++)
        {
            TecUtilZoneCopy(critPointZone,
                            atom, atom, 1, // i-range
                            1, 1, 1,       // j-range
                            1, 1, 1);      // k-range
            CHECK(TecUtilDataSetGetNumZones()==numZones+atom);
            (*atomZoneNumbers)[atom-1] = numZones+atom;

            char atomName[500];
            sprintf_s(atomName, "Atom %d", atom);
            TecUtilZoneRename(numZones+atom, atomName);
            AuxData_pa zoneAuxDataRef = TecUtilAuxDataZoneGetRef(numZones+atom);
            if (zoneAuxDataRef != NULL)
            {
                TecUtilAuxDataSetStrItem(zoneAuxDataRef, "CompChem.ZoneType", "Atom", TRUE);
                TecUtilAuxDataSetStrItem(zoneAuxDataRef, "CompChem.AtomType", "Unknown", TRUE);
                TecUtilAuxDataSetStrItem(zoneAuxDataRef, "CompChem.AtomAbbr", "??", TRUE);
                TecUtilAuxDataDeleteItemByName(zoneAuxDataRef, "CompChem.NumCrtPtAtom");
                TecUtilAuxDataDeleteItemByName(zoneAuxDataRef, "CompChem.NumCrtPtBond");
                TecUtilAuxDataDeleteItemByName(zoneAuxDataRef, "CompChem.NumCrtPtRing");
                TecUtilAuxDataDeleteItemByName(zoneAuxDataRef, "CompChem.NumCrtPtCage");
            }
            VERIFY(TecUtilSetAddMember(activeZoneSet, numZones+atom, TRUE));
            VERIFY(TecUtilSetAddMember(atomsZoneSet, numZones+atom, TRUE));
        }
        TecUtilSetRemoveMember(activeZoneSet, critPointZone);
        TecUtilZoneSetActive(activeZoneSet, AssignOp_Equals);
        TecUtilSetDealloc(&activeZoneSet);

        // now copy scatter style from crit-point zone to new atom zones: in an ideal world this would be less code
        ArgList_pa argList = TecUtilArgListAlloc();

        // get symbol shape
        ArbParam_t symbolShape;
        TecUtilArgListClear(argList);
        TecUtilArgListAppendString(argList, SV_P1, SV_FIELDMAP);
        TecUtilArgListAppendString(argList, SV_P2, SV_SCATTER);
        TecUtilArgListAppendString(argList, SV_P3, SV_SYMBOLSHAPE);
        TecUtilArgListAppendString(argList, SV_P4, SV_GEOMSHAPE);
        TecUtilArgListAppendInt(argList, SV_OFFSET1, critPointZone);
        TecUtilArgListAppendArbParamPtr(argList, SV_IVALUE, &symbolShape);
        TecUtilStyleGetLowLevelX(argList);
        // set symbol shape
        TecUtilArgListClear(argList);
        TecUtilArgListAppendString(argList, SV_P1, SV_FIELDMAP);
        TecUtilArgListAppendString(argList, SV_P2, SV_SCATTER);
        TecUtilArgListAppendString(argList, SV_P3, SV_SYMBOLSHAPE);
        TecUtilArgListAppendString(argList, SV_P4, SV_GEOMSHAPE);
        TecUtilArgListAppendSet(argList, SV_OBJECTSET, atomsZoneSet);
        TecUtilArgListAppendArbParam(argList, SV_IVALUE, symbolShape);
        TecUtilStyleSetLowLevelX(argList);

        // get symbol color
        ArbParam_t symbolColor;
        TecUtilArgListClear(argList);
        TecUtilArgListAppendString(argList, SV_P1, SV_FIELDMAP);
        TecUtilArgListAppendString(argList, SV_P2, SV_SCATTER);
        TecUtilArgListAppendString(argList, SV_P3, SV_COLOR);
        TecUtilArgListAppendInt(argList, SV_OFFSET1, critPointZone);
        TecUtilArgListAppendArbParamPtr(argList, SV_IVALUE, &symbolColor);
        TecUtilStyleGetLowLevelX(argList);
        // set symbol color
        TecUtilArgListClear(argList);
        TecUtilArgListAppendString(argList, SV_P1, SV_FIELDMAP);
        TecUtilArgListAppendString(argList, SV_P2, SV_SCATTER);
        TecUtilArgListAppendString(argList, SV_P3, SV_COLOR);
        TecUtilArgListAppendSet(argList, SV_OBJECTSET, atomsZoneSet);
        TecUtilArgListAppendArbParam(argList, SV_IVALUE, symbolColor);
        TecUtilStyleSetLowLevelX(argList);

        // get symbol size
        double symbolSize;
        TecUtilArgListClear(argList);
        TecUtilArgListAppendString(argList, SV_P1, SV_FIELDMAP);
        TecUtilArgListAppendString(argList, SV_P2, SV_SCATTER);
        TecUtilArgListAppendString(argList, SV_P3, SV_FRAMESIZE);
        TecUtilArgListAppendInt(argList, SV_OFFSET1, critPointZone);
        TecUtilArgListAppendDoublePtr(argList, SV_DVALUE, &symbolSize);
        TecUtilStyleGetLowLevelX(argList);
        // set symbol size
        TecUtilArgListClear(argList);
        TecUtilArgListAppendString(argList, SV_P1, SV_FIELDMAP);
        TecUtilArgListAppendString(argList, SV_P2, SV_SCATTER);
        TecUtilArgListAppendString(argList, SV_P3, SV_FRAMESIZE);
        TecUtilArgListAppendSet(argList, SV_OBJECTSET, atomsZoneSet);
        TecUtilArgListAppendDouble(argList, SV_DVALUE, symbolSize);
        TecUtilStyleSetLowLevelX(argList);

        // get symbol size-by-variable state
        ArbParam_t sizeByVariableState;
        TecUtilArgListClear(argList);
        TecUtilArgListAppendString(argList, SV_P1, SV_FIELDMAP);
        TecUtilArgListAppendString(argList, SV_P2, SV_SCATTER);
        TecUtilArgListAppendString(argList, SV_P3, SV_SIZEBYVARIABLE);
        TecUtilArgListAppendInt(argList, SV_OFFSET1, critPointZone);
        TecUtilArgListAppendArbParamPtr(argList, SV_IVALUE, &sizeByVariableState);
        TecUtilStyleGetLowLevelX(argList);
        // set symbol size-by-variable state
        TecUtilArgListClear(argList);
        TecUtilArgListAppendString(argList, SV_P1, SV_FIELDMAP);
        TecUtilArgListAppendString(argList, SV_P2, SV_SCATTER);
        TecUtilArgListAppendString(argList, SV_P3, SV_SIZEBYVARIABLE);
        TecUtilArgListAppendSet(argList, SV_OBJECTSET, atomsZoneSet);
        TecUtilArgListAppendArbParam(argList, SV_IVALUE, sizeByVariableState);
        TecUtilStyleSetLowLevelX(argList);

        TecUtilArgListDealloc(&argList);
        TecUtilSetDealloc(&atomsZoneSet);
    }
}

/**
 */
static void Sidebar1Activate_CB(void)
{
    /*  <<< This function is called when sidebar "Bondalyzer" is activated >>> */

    if ( !TecUtilDataSetIsAvailable() ||
         TecUtilDataSetGetNumZones() < 2 ||
         ( TecUtilFrameGetPlotType() != PlotType_Cartesian2D &&
           TecUtilFrameGetPlotType() != PlotType_Cartesian3D ) )
    {
        // for some reason, doing an error message here doesn't work, so it is in the menu callback
        // TecUtilDialogErrMsg("Current plot does not appear to be from ChemBond");
        TecGUISidebarActivate(TECGUITECPLOTSIDEBAR);
        return;
    }

    ExtractTopoInZonesAuxData(&Bundles, &CritPoints);
    atomZoneInfoAlloc(&atomZoneNumbers);

    // turn status line, toolbar, and menu off
    TecUtilStyleSetLowLevel(NULL, 0.0, (ArbParam_t)FALSE, NULL, AssignOp_Equals,
        SV_INTERFACE, SV_SHOWSTATUSLINE, NULL, NULL, NULL, NULL, FALSE);
    TecUtilToolbarActivate(FALSE);

#if 1
    TecUtilMenuActivate(FALSE);
#else
    TecUtilMenuClearAll();
#endif

    EntIndex_t numZones = TecUtilDataSetGetNumZones();

    // get a copy of the original zones that are active
    if ( orgActiveZoneSet == NULL )
        VERIFY(TecUtilZoneGetActive(&orgActiveZoneSet));
    // get a copy of the original zones mesh thickness
    if ( orgZoneProperties == NULL )
    {
        orgZoneProperties = ALLOC_ARRAY(numZones,OrgZoneProperties_s,"orgZoneProperties");
        CHECK(VALID_REF(orgZoneProperties));
    }

    ArgList_pa argList = TecUtilArgListAlloc();
    CHECK(VALID_REF(argList));
    for ( EntIndex_t zone = 1; zone <= numZones; zone++ )
    {
        ArbParam_t arbParam;
        double doubleVal = 0.0;

        TecUtilArgListClear(argList);
        TecUtilArgListAppendInt(argList,       SV_OFFSET1, zone);
        TecUtilArgListAppendString(argList,    SV_P1,      SV_FIELDMAP);
        TecUtilArgListAppendString(argList,    SV_P2,      SV_MESH);
        TecUtilArgListAppendString(argList,    SV_P3,      SV_LINETHICKNESS);
        TecUtilArgListAppendDoublePtr(argList, SV_DVALUE,  &doubleVal);
        TecUtilStyleGetLowLevelX(argList);
        orgZoneProperties[zone-1].meshThickness = doubleVal;

        TecUtilArgListClear(argList);
        TecUtilArgListAppendInt(argList,       SV_OFFSET1, zone);
        TecUtilArgListAppendString(argList,    SV_P1,      SV_FIELDMAP);
        TecUtilArgListAppendString(argList,    SV_P2,      SV_MESH);
        TecUtilArgListAppendString(argList,    SV_P3,      SV_COLOR);
        TecUtilArgListAppendArbParamPtr(argList, SV_IVALUE, &arbParam);
        TecUtilStyleGetLowLevelX(argList);
        orgZoneProperties[zone-1].meshColor = (ColorIndex_t)arbParam;

        TecUtilArgListClear(argList);
        TecUtilArgListAppendInt(argList,       SV_OFFSET1, zone);
        TecUtilArgListAppendString(argList,    SV_P1,      SV_FIELDMAP);
        TecUtilArgListAppendString(argList,    SV_P2,      SV_SCATTER);
        TecUtilArgListAppendString(argList,    SV_P3,      SV_COLOR);
        TecUtilArgListAppendArbParamPtr(argList, SV_IVALUE, &arbParam);
        TecUtilStyleGetLowLevelX(argList);
        orgZoneProperties[zone-1].scatterColor = (ColorIndex_t)arbParam;
    }
    TecUtilArgListDealloc(&argList);

    setDisplayMode(DisplayMode_AllAtomsAndBonds);

    updateBondOrAtomList(DisplayMode_AllAtomsAndBonds);
    TecUtilMouseSetMode(MouseButtonMode_RotateSpherical);
}


void setMeshThicknessForZone(SetIndex_t zone,
                             double     meshThickness,
                             Set_pa     tempSet, // to avoid allocating/deallocing same set over and over
                             ArgList_pa tempArgList)
{
    REQUIRE(1<=zone && zone<=TecUtilDataSetGetNumZones());
    REQUIRE(meshThickness>0.0);
    REQUIRE(VALID_REF(tempSet) && TecUtilSetIsEmpty(tempSet));
    REQUIRE(VALID_REF(tempArgList));

    TecUtilSetAddMember(tempSet, zone, TRUE);
    TecUtilArgListClear(tempArgList);
    TecUtilArgListAppendSet(tempArgList,    SV_OBJECTSET, tempSet);
    TecUtilArgListAppendString(tempArgList, SV_P1,        SV_FIELDMAP);
    TecUtilArgListAppendString(tempArgList, SV_P2,        SV_MESH);
    TecUtilArgListAppendString(tempArgList, SV_P3,        SV_LINETHICKNESS);
    TecUtilArgListAppendDouble(tempArgList, SV_DVALUE,    meshThickness);
    TecUtilStyleSetLowLevelX(tempArgList);
    TecUtilSetRemoveMember(tempSet, zone);
}

void setMeshColorForZone(SetIndex_t   zone,
                         ColorIndex_t meshColor,
                         Set_pa       tempSet, // to avoid allocating/deallocing same set over and over
                         ArgList_pa   tempArgList)
{
    REQUIRE(1<=zone && zone<=TecUtilDataSetGetNumZones());
    REQUIRE(meshColor!=InvalidColor_C);
    REQUIRE(VALID_REF(tempSet) && TecUtilSetIsEmpty(tempSet));
    REQUIRE(VALID_REF(tempArgList));

    TecUtilSetAddMember(tempSet, zone, TRUE);
    TecUtilArgListClear(tempArgList);
    TecUtilArgListAppendSet(tempArgList,      SV_OBJECTSET, tempSet);
    TecUtilArgListAppendString(tempArgList,   SV_P1,        SV_FIELDMAP);
    TecUtilArgListAppendString(tempArgList,   SV_P2,        SV_MESH);
    TecUtilArgListAppendString(tempArgList,   SV_P3,        SV_COLOR);
    TecUtilArgListAppendArbParam(tempArgList, SV_IVALUE,    meshColor);
    TecUtilStyleSetLowLevelX(tempArgList);
    TecUtilSetRemoveMember(tempSet, zone);
}

void setScatterColorForZone(SetIndex_t   zone,
                            ColorIndex_t scatterColor,
                            Set_pa       tempSet, // to avoid allocating/deallocing same set over and over
                            ArgList_pa   tempArgList)
{
    REQUIRE(1<=zone && zone<=TecUtilDataSetGetNumZones());
    REQUIRE(scatterColor!=InvalidColor_C);
    REQUIRE(VALID_REF(tempSet) && TecUtilSetIsEmpty(tempSet));
    REQUIRE(VALID_REF(tempArgList));

    TecUtilSetAddMember(tempSet, zone, TRUE);
    TecUtilArgListClear(tempArgList);
    TecUtilArgListAppendSet(tempArgList,      SV_OBJECTSET, tempSet);
    TecUtilArgListAppendString(tempArgList,   SV_P1,        SV_FIELDMAP);
    TecUtilArgListAppendString(tempArgList,   SV_P2,        SV_SCATTER);
    TecUtilArgListAppendString(tempArgList,   SV_P3,        SV_COLOR);
    TecUtilArgListAppendArbParam(tempArgList, SV_IVALUE,    scatterColor);
    TecUtilStyleSetLowLevelX(tempArgList);
    TecUtilSetRemoveMember(tempSet, zone);
}

/*
*/
void restoreStyleForZone(SetIndex_t zone,
                         Set_pa     tempSet, // to avoid allocating/deallocing same set over and over
                         ArgList_pa tempArgList)
{
    REQUIRE(1<=zone && zone<=TecUtilDataSetGetNumZones());
    REQUIRE(VALID_REF(tempSet) && TecUtilSetIsEmpty(tempSet));
    REQUIRE(VALID_REF(tempArgList));

    setMeshThicknessForZone(zone, orgZoneProperties[zone-1].meshThickness, tempSet, tempArgList);
    setMeshColorForZone(zone, orgZoneProperties[zone-1].meshColor, tempSet, tempArgList);
    setScatterColorForZone(zone, orgZoneProperties[zone-1].scatterColor, tempSet, tempArgList);
}


/*
*/
void restoreStyleForZoneSet(Set_pa zoneSet)
{
    REQUIRE(VALID_REF(zoneSet));
    ArgList_pa argList = TecUtilArgListAlloc();
    CHECK(VALID_REF(argList));
    Set_pa objectSet = TecUtilSetAlloc(TRUE);
    CHECK(VALID_REF(objectSet));

    SetIndex_t zone = TecUtilSetGetNextMember(zoneSet, TECUTILSETNOTMEMBER);
    while ( zone != TECUTILSETNOTMEMBER )
    {
        restoreStyleForZone(zone, objectSet, argList);
        zone = TecUtilSetGetNextMember(zoneSet, zone);
    }
    TecUtilSetDealloc(&objectSet);
    TecUtilArgListDealloc(&argList);
}


/*
*/
void restoreStyleForAllZones(void)
{
    ArgList_pa argList = TecUtilArgListAlloc();
    CHECK(VALID_REF(argList));
    Set_pa objectSet = TecUtilSetAlloc(TRUE);
    CHECK(VALID_REF(objectSet));

    EntIndex_t numZones = TecUtilDataSetGetNumZones();
    for ( EntIndex_t zone = 1; zone <= numZones; zone++ )
        restoreStyleForZone(zone, objectSet, argList);

    TecUtilSetDealloc(&objectSet);
    TecUtilArgListDealloc(&argList);
}


/**
 */
static void Sidebar1Deactivate_CB(void)
{
    /* <<< This function is called when sidebar "Bondalyzer" is deactivated >>> */

    // turn status line, toolbar, and main menu on
    TecUtilStyleSetLowLevel(NULL, 0.0, (ArbParam_t)TRUE, NULL, AssignOp_Equals,
        SV_INTERFACE, SV_SHOWSTATUSLINE, NULL, NULL, NULL, NULL, FALSE);
    TecUtilToolbarActivate(TRUE);

#if 1
    TecUtilMenuActivate(TRUE);
#else
    // This would be so much easier (and better) if we added TecUtilMenuActivate(Boolean_t)
    TecUtilMenuClearAll();
    Menu_pa main = TecUtilMenuGetMain();
    TecUtilMenuInsertStandard(main, MENU_POSITION_LAST, StandardMenu_File);
    TecUtilMenuInsertStandard(main, MENU_POSITION_LAST, StandardMenu_Edit);
    TecUtilMenuInsertStandard(main, MENU_POSITION_LAST, StandardMenu_View);
    TecUtilMenuInsertStandard(main, MENU_POSITION_LAST, StandardMenu_Plot);
    TecUtilMenuInsertStandard(main, MENU_POSITION_LAST, StandardMenu_Insert);
    TecUtilMenuInsertStandard(main, MENU_POSITION_LAST, StandardMenu_Animate);
    TecUtilMenuInsertStandard(main, MENU_POSITION_LAST, StandardMenu_Data);
    TecUtilMenuInsertStandard(main, MENU_POSITION_LAST, StandardMenu_Frame);
    TecUtilMenuInsertStandard(main, MENU_POSITION_LAST, StandardMenu_Options);
    TecUtilMenuInsertStandard(main, MENU_POSITION_LAST, StandardMenu_Scripting);
    TecUtilMenuInsertStandard(main, MENU_POSITION_LAST, StandardMenu_Tools);
    TecUtilMenuInsertStandard(main, MENU_POSITION_LAST, StandardMenu_Help);
    // can't do anything about other add-ons, but we can re-add ourselves.
    registerMenuCallbacks();
#endif

    // reset style to default and free info used by Bondalyzer
    if ( orgZoneProperties != NULL )
    {
        EntIndex_t numZones = TecUtilDataSetGetNumZones();

        restoreStyleForAllZones();
        FREE_ARRAY(orgZoneProperties, "orgZoneProperties");
    }
    orgZoneProperties = NULL;

    if ( orgActiveZoneSet != NULL )
    {
        TecUtilZoneSetActive(orgActiveZoneSet, AssignOp_Equals);
        TecUtilSetDealloc(&orgActiveZoneSet);
    }
    orgActiveZoneSet = NULL;

    CritPointsDealloc(&CritPoints);
    BundlesDealloc(&Bundles);
    atomZoneInfoDealloc(&atomZoneNumbers);
}


/*
 */
static void updateAtomAndBondZones(LgIndex_t     whichAtom, // -1==none
                                   LgIndex_t     whichBond, // -1==none
                                   DisplayMode_e displayMode)
{
    REQUIRE(whichAtom == -1 || (1<=whichAtom && whichAtom<=CritPoints->NumCrtPtsM3));
    REQUIRE(whichBond == -1 || (1<=whichBond && whichBond<=CritPoints->NumCrtPtsM1));
    REQUIRE(VALID_ENUM(displayMode,DisplayMode_e));
    REQUIRE(VALID_REF(orgActiveZoneSet));

    // create a list of atom zones
    ArrList_pa atomZoneList = NULL;
    if ( whichAtom != -1 )
    {
        if ( displayMode == DisplayMode_SingleBaderAtom )
            atomZoneList = NULL; // not supported yet
        else
        {
            atomZoneList = ArrListAlloc(1, ArrListType_Long);
            CHECK(VALID_REF(atomZoneNumbers));
            EntIndex_t atomZone = atomZoneNumbers[whichAtom-1];
            CHECK(1<=atomZone && atomZone<=TecUtilDataSetGetNumZones());
            ArrListItem_u item;
            item.Long = atomZone;
            VERIFY(ArrListAppendItem(atomZoneList, item));
        }
    }
    
    // create a list of bond zones
    ArrList_pa bondZoneList = NULL;
    if ( whichBond != -1 )
    {
        EntIndex_t numBonds = 1;
        ArrList_pa bondList = ArrListAlloc(numBonds, ArrListType_Long);
        CHECK(VALID_REF(bondList));
        ArrListItem_u item;
        item.Long = whichBond-1;
        VERIFY(ArrListAppendItem(bondList, item));
        // determine the zones correspond to whichBonds and displayMode
        if ( displayMode == DisplayMode_SingleBondBundle )
            bondZoneList = ZoneListOfBundleZonesForCritPoint(FALSE, TRUE, (char)(-1), bondList, CritPoints, Bundles);
        else
            bondZoneList = ZoneListOfBondLinesForCritPoint((char)(-1), bondList, CritPoints, Bundles);
    }

    // create a copy of the original zone set
    Set_pa newActiveZoneSet = TecUtilSetAlloc(TRUE);
    CHECK(VALID_REF(newActiveZoneSet));
    VERIFY(TecUtilSetCopy(newActiveZoneSet, orgActiveZoneSet, TRUE));

    // add the atom zones
    if ( atomZoneList != NULL )
    {
        EntIndex_t numEntries = ArrListGetCount(atomZoneList);
        for (EntIndex_t entry = 0; entry < numEntries; entry++)
        {
            ArrListItem_u item = ArrListGetItem(atomZoneList, entry);
            VERIFY(TecUtilSetAddMember(newActiveZoneSet, item.Long, TRUE));
        }
    }

    // add the bond zones
    if ( bondZoneList != NULL )
    {
        EntIndex_t numEntries = ArrListGetCount(bondZoneList);
        for (EntIndex_t entry = 0; entry < numEntries; entry++)
        {
            ArrListItem_u item = ArrListGetItem(bondZoneList, entry);
            VERIFY(TecUtilSetAddMember(newActiveZoneSet, item.Long, TRUE));
        }
    }

    // activate the zones
    TecUtilZoneSetActive(newActiveZoneSet, AssignOp_Equals);

    // restore style for active zones
    restoreStyleForZoneSet(newActiveZoneSet);

    ArgList_pa argList = TecUtilArgListAlloc();
    Set_pa objectSet = TecUtilSetAlloc(TRUE);
    CHECK(VALID_REF(argList) && VALID_REF(objectSet));

    // modify style for atoms as needed
    #define highlightAtomsForDisplayMode(displayMode) ( (displayMode) == DisplayMode_SingleAtom || (displayMode) == DisplayMode_AllAtomsAndBonds )
    if ( highlightAtomsForDisplayMode(displayMode) && atomZoneList != NULL )
    {
        EntIndex_t numEntries = ArrListGetCount(atomZoneList);
        for (EntIndex_t entry = 0; entry < numEntries; entry++)
        {
            ArrListItem_u item = ArrListGetItem(atomZoneList, entry);
            EntIndex_t zone = item.Long;
            if ( orgZoneProperties[zone-1].scatterColor == White_C ||
                 orgZoneProperties[zone-1].scatterColor == Yellow_C )
                setScatterColorForZone(zone, Custom32_C, objectSet, argList); // dark orange
            else
                setScatterColorForZone(zone, Custom3_C, objectSet, argList); // light orange
        }
    }

    // modify style for bonds as needed
    #define highlightBondsForDisplayMode(displayMode) ( (displayMode) == DisplayMode_SingleBond || (displayMode) == DisplayMode_AllAtomsAndBonds )
    if ( highlightBondsForDisplayMode(displayMode) && bondZoneList != NULL )
    {
        EntIndex_t numEntries = ArrListGetCount(bondZoneList);
        for (EntIndex_t entry = 0; entry < numEntries; entry++)
        {
            ArrListItem_u item = ArrListGetItem(bondZoneList, entry);
            EntIndex_t zone = item.Long;
            setMeshThicknessForZone(zone, 2.0*orgZoneProperties[zone-1].meshThickness, objectSet, argList);
            if ( orgZoneProperties[zone-1].meshColor == Black_C ||
                 orgZoneProperties[zone-1].meshColor == Blue_C ||
                 orgZoneProperties[zone-1].meshColor == Custom1_C )
                setMeshColorForZone(zone, Cyan_C, objectSet, argList);
            else
                setMeshColorForZone(zone, Black_C, objectSet, argList);
        }
    }

    TecUtilSetDealloc(&objectSet);
    TecUtilArgListDealloc(&argList);

    // clean-up
    TecUtilSetDealloc(&newActiveZoneSet);
    ArrListDealloc(&atomZoneList);
    ArrListDealloc(&bondZoneList);

    TecUtilMouseSetMode(MouseButtonMode_RotateRollerBall);
}


/**
 */
static void AtomOrBond_MLST_S1_CB(const LgIndex_t *I)
{
  TecUtilLockStart(AddOnID);
  TRACE1("Multi selection list (AtomOrBond_MLST_S1) item selected,  First Item is: %d\n",*I);
  LgIndex_t selectedItem = TecGUIListGetSelectedItem(AtomOrBond_MLST_S1);
  switch (curDisplayMode)
  {
      case DisplayMode_AllAtomsAndBonds :
          {
              if ( selectedItem == -1 )
              {
                  // nothing Selected
                  curAtom = -1;
                  curBond = -1;
              }
              else if ( selectedItem <= CritPoints->NumCrtPtsM3 )
              {
                  // atom selected
                  curAtom = selectedItem;
                  curBond = -1;
              }
              else if ( selectedItem <= CritPoints->NumCrtPtsM3 + CritPoints->NumCrtPtsM1 )
              {
                  // whichBond selected
                  curAtom = -1;
                  curBond = selectedItem - CritPoints->NumCrtPtsM3;
              }
              else
                  CHECK(FALSE);
          } break;
      case DisplayMode_SingleBaderAtom :
      case DisplayMode_SingleAtom :
          {
              if ( selectedItem == -1 )
              {
                  // nothing Selected
                  curAtom = -1;
                  curBond = -1;
              }
              else
              {
                  // atom selected
                  CHECK(selectedItem <= CritPoints->NumCrtPtsM3);
                  curAtom = selectedItem;
                  curBond = -1;
              }
          } break;
      case DisplayMode_SingleBond :
      case DisplayMode_SingleBondBundle :
          {
              if ( selectedItem == -1 )
              {
                  // nothing Selected
                  curAtom = -1;
                  curBond = -1;
              }
              else
              {
                  // whichBond selected
                  CHECK(selectedItem <= CritPoints->NumCrtPtsM1);
                  curAtom = -1;
                  curBond = selectedItem;
              }
          } break;
      default:
          CHECK(FALSE);
  }
  updateAtomAndBondZones(curAtom, curBond, curDisplayMode);
  updateInfoField();

  TecUtilLockFinish(AddOnID);
}

/**
 */
static void ReturnToTecp_BTN_S1_CB(void)
{
  TecUtilLockStart(AddOnID);
  TRACE("Return to Tecplot Button Pushed\n");

  TecGUISidebarActivate(TECGUITECPLOTSIDEBAR);
  TecUtilLockFinish(AddOnID);
}


char *WhatToShow_OPT_S1_List = "All Atoms & Bonds,Single Atom + Connections,Single Bond + Connections,Single Bond Bundle";



/**
 */
static void WhatToShow_OPT_S1_CB(const LgIndex_t *I)
{
  TecUtilLockStart(AddOnID);
  TRACE1("Option Menu (WhatToShow_OPT_S1) value changed,  New value is: %d\n",*I);
  DisplayMode_e displayMode = getDisplayModeFromLgIndex(*I);
  if ( VALID_ENUM(displayMode, DisplayMode_e) )
      setDisplayMode(displayMode);
  TecUtilLockFinish(AddOnID);
}


//
// Dialog version
//


ArrList_pa ZoneListOld = NULL;

/*
 * Set the Bond, Atom, Ring, and Cage list sensitivity.
 */
void  SetListSensitivity()
{
    LgIndex_t ActiveList = TecGUIOptionMenuGet(CPType_OPT_D1);

    // Set the sensititivity for the choice
    switch (ActiveList)
    {
        case 1:
            TecGUISetSensitivity(BondList_MLST_D1, TRUE);
            TecGUISetSensitivity(AtomList_MLST_D1, FALSE);
            TecGUISetSensitivity(RingList_MLST_D1, FALSE);
            TecGUISetSensitivity(CageList_MLST_D1, FALSE);
            break;
        case 2:
            TecGUISetSensitivity(BondList_MLST_D1, FALSE);
            TecGUISetSensitivity(AtomList_MLST_D1, TRUE);
            TecGUISetSensitivity(RingList_MLST_D1, FALSE);
            TecGUISetSensitivity(CageList_MLST_D1, FALSE);
            break;
        case 3:
            TecGUISetSensitivity(BondList_MLST_D1, FALSE);
            TecGUISetSensitivity(AtomList_MLST_D1, FALSE);
            TecGUISetSensitivity(RingList_MLST_D1, TRUE);
            TecGUISetSensitivity(CageList_MLST_D1, FALSE);
            break;
        case 4:
            TecGUISetSensitivity(BondList_MLST_D1, FALSE);
            TecGUISetSensitivity(AtomList_MLST_D1, FALSE);
            TecGUISetSensitivity(RingList_MLST_D1, FALSE);
            TecGUISetSensitivity(CageList_MLST_D1, TRUE);
            break;
        default:
            TecGUISetSensitivity(BondList_MLST_D1, TRUE);
            TecGUISetSensitivity(AtomList_MLST_D1, TRUE);
            TecGUISetSensitivity(RingList_MLST_D1, TRUE);
            TecGUISetSensitivity(CageList_MLST_D1, TRUE);

    }
}


/*
 * Load the Initial Bonds multi-selection list
 */
Boolean_t InitFullMSLists()
{
    Boolean_t IsOk = TRUE;
    LgIndex_t ActiveList = TecGUIOptionMenuGet(CPType_OPT_D1);
    int ii;

    LgIndex_t NumBonds = CritPoints->NumCrtPtsM1;
    LgIndex_t NumAtoms = CritPoints->NumCrtPtsM3;
    LgIndex_t NumRings = CritPoints->NumCrtPtsP1;
    LgIndex_t NumCages = CritPoints->NumCrtPtsP3;

    TecGUIListDeleteAllItems(BondList_MLST_D1);
    for (ii = 0; ii < NumBonds; ii++)
    {
        char string[10];
        sprintf_s(string, "%d", ii + 1);
        TecGUIListAppendItem(BondList_MLST_D1, string);
    }
    if (ActiveList == 1) TecGUIListSetSelectedItem(BondList_MLST_D1, 1);

    TecGUIListDeleteAllItems(AtomList_MLST_D1);
    for (ii = 0; ii < NumAtoms; ii++)
    {
        char string[10];
        sprintf_s(string, "%d", ii + 1);
        TecGUIListAppendItem(AtomList_MLST_D1, string);
    }
    if (ActiveList == 2) TecGUIListSetSelectedItem(AtomList_MLST_D1, 1);

    TecGUIListDeleteAllItems(RingList_MLST_D1);
    for (ii = 0; ii < NumRings; ii++)
    {
        char string[10];
        sprintf_s(string, "%d", ii + 1);
        TecGUIListAppendItem(RingList_MLST_D1, string);
    }
    if (ActiveList == 3 && NumRings > 0) TecGUIListSetSelectedItem(RingList_MLST_D1, 1);

    TecGUIListDeleteAllItems(CageList_MLST_D1);
    for (ii = 0; ii < NumCages; ii++)
    {
        char string[10];
        sprintf_s(string, "%d", ii + 1);
        TecGUIListAppendItem(CageList_MLST_D1, string);
    }
    if (ActiveList == 4 && NumCages > 0) TecGUIListSetSelectedItem(CageList_MLST_D1, 1);


    return IsOk;
}


/*
 * Load the specified multi-selection list
 *
 * param
 *    ListID: ID of desired list
 * param
 *    CrtPointList: Array list with zero-based values of critical points
 */
void      SetMSList(LgIndex_t  ListID,
                    ArrList_pa CrtPointList)
{
    int ii;
    LgIndex_t NumItems;

    REQUIRE(ArrListIsValid(CrtPointList));

    NumItems = ArrListGetCount(CrtPointList);
    for (ii = 0; ii < NumItems; ii++)
    {
        ArrListItem_u Item = ArrListGetItem(CrtPointList, ii);
        char string[10];
        sprintf_s(string, "%d", Item.Long + 1);
        TecGUIListAppendItem(ListID, string);
    }

}


/*
 * Get the selected item numbers in a multi-selection list and
 * return as an array list.
 *
 * param
 *   ID of the list to get selected values from.
 *
 * return
 *   Array list of zero-based critical point values (converted from one-based)
 */
ArrList_pa  GetSelectedInMSList(LgIndex_t ListID)
{
    ArrList_pa Result = NULL;
    Boolean_t IsOk = TRUE;

    LgIndex_t count;
    LgIndex_t *sel;
    LgIndex_t ii;

    TecGUIListGetSelectedItems(ListID, &sel, &count);
    if (count >= 1)
    {
        ArrListItem_u Item;
        Result = ArrListAlloc(count, ArrListType_Long);
        for (ii = 0; ii < count; ii++)
        {
            Item.Long = sel[ii] - 1;
            IsOk = ArrListAppendItem(Result, Item);
        }
    }
    else // If none selected, assume first item selected
    {
        ArrListItem_u Item;
        Result = ArrListAlloc(count, ArrListType_Long);
        Item.Long = 0;
        IsOk = ArrListAppendItem(Result, Item);
    }


    TecUtilArrayDealloc((void **)&sel); // Clean up when done.

    ENSURE(Result == NULL || ArrListIsValid(Result));
    return Result;
}



/*
 * Based on the CPType option and the independent list (Bond, Atom, Ring,
 * or Cage) list selections, reset the three dependent (Bond, Atom, Ring,
 * or Cage) list items.
 */
ArrList_pa  SetListItems()
{
    ArrList_pa Result = NULL;
    ArrList_pa BondList = NULL;
    ArrList_pa AtomList = NULL;
    ArrList_pa RingList = NULL;
    ArrList_pa CageList = NULL;

    LgIndex_t ActiveList = TecGUIOptionMenuGet(CPType_OPT_D1);

    LgIndex_t NumBonds = CritPoints->NumCrtPtsM1;
    LgIndex_t NumAtoms = CritPoints->NumCrtPtsM3;
    LgIndex_t NumRings = CritPoints->NumCrtPtsP1;
    LgIndex_t NumCages = CritPoints->NumCrtPtsP3;

    // Determine and set the dependent list items for the choice
    switch (ActiveList)
    {
        case 1:  // Bond is independent list
            // Get the selected whichBond list items
            BondList = GetSelectedInMSList(BondList_MLST_D1);

            // Determine and set the connecting Atoms, Rings, and Cages.
            // Atoms
            TecGUIListDeleteAllItems(AtomList_MLST_D1);
            AtomList = ListOfCPOutsConnetedToCPIns((char)(-1), (char)(-3), BondList, Bundles);
            SetMSList(AtomList_MLST_D1, AtomList);
            // Rings
            TecGUIListDeleteAllItems(RingList_MLST_D1);
            RingList = ListOfCPOutsConnetedToCPIns((char)(-1), (char)(1), BondList, Bundles);
            SetMSList(RingList_MLST_D1, RingList);
            // Cages
            TecGUIListDeleteAllItems(CageList_MLST_D1);
            CageList = ListOfCPOutsConnetedToCPIns((char)(-1), (char)(3), BondList, Bundles);
            SetMSList(CageList_MLST_D1, CageList);

            // Set independent list to result and dealloc dependent lists
            ArrListDealloc(&AtomList);
            ArrListDealloc(&RingList);
            ArrListDealloc(&CageList);
            Result = BondList;
            break;

        case 2:  // Atom is independent list
            // Get the selected Atom list items
            AtomList = GetSelectedInMSList(AtomList_MLST_D1);

            // Determine and set the connecting Bonds, Rings, and Cages.
            // Bonds
            TecGUIListDeleteAllItems(BondList_MLST_D1);
            BondList = ListOfCPOutsConnetedToCPIns((char)(-3), (char)(-1), AtomList, Bundles);
            SetMSList(BondList_MLST_D1, BondList);
            // Rings
            TecGUIListDeleteAllItems(RingList_MLST_D1);
            RingList = ListOfCPOutsConnetedToCPIns((char)(-3), (char)(1), AtomList, Bundles);
            SetMSList(RingList_MLST_D1, RingList);
            // Cages
            TecGUIListDeleteAllItems(CageList_MLST_D1);
            CageList = ListOfCPOutsConnetedToCPIns((char)(-3), (char)(3), AtomList, Bundles);
            SetMSList(CageList_MLST_D1, CageList);

            // Set independent list to result and dealloc dependent lists
            ArrListDealloc(&BondList);
            ArrListDealloc(&RingList);
            ArrListDealloc(&CageList);
            Result = AtomList;
            break;
        case 3:  // Ring is independent list
            if (NumRings > 0)
            {
                // Get the selected Ring list items
                RingList = GetSelectedInMSList(RingList_MLST_D1);

                // Determine and set the connecting Bonds, Rings, and Cages.
                // Atoms
                TecGUIListDeleteAllItems(AtomList_MLST_D1);
                AtomList = ListOfCPOutsConnetedToCPIns((char)(1), (char)(-3), RingList, Bundles);
                SetMSList(AtomList_MLST_D1, AtomList);
                // Bonds
                TecGUIListDeleteAllItems(BondList_MLST_D1);
                BondList = ListOfCPOutsConnetedToCPIns((char)(1), (char)(-1), RingList, Bundles);
                SetMSList(BondList_MLST_D1, BondList);
                // Cages
                TecGUIListDeleteAllItems(CageList_MLST_D1);
                CageList = ListOfCPOutsConnetedToCPIns((char)(1), (char)(3), RingList, Bundles);
                SetMSList(CageList_MLST_D1, CageList);

                // Set independent list to result and dealloc dependent lists
                ArrListDealloc(&BondList);
                ArrListDealloc(&AtomList);
                ArrListDealloc(&CageList);
                Result = RingList;
            }
            else
            {
                TecUtilDialogMessageBox("No Rings in topology", MessageBox_Warning);
                Result = ArrListAlloc(10, ArrListType_Long);
            }
            break;
        case 4:  // Cage is independent list
            if (NumCages > 0)
            {
                // Get the selected Cage list items
                CageList = GetSelectedInMSList(CageList_MLST_D1);

                // Determine and set the connecting Bonds, Rings, and Cages.
                // Atoms
                TecGUIListDeleteAllItems(AtomList_MLST_D1);
                AtomList = ListOfCPOutsConnetedToCPIns((char)(3), (char)(-3), CageList, Bundles);
                SetMSList(AtomList_MLST_D1, AtomList);
                // Bonds
                TecGUIListDeleteAllItems(BondList_MLST_D1);
                BondList = ListOfCPOutsConnetedToCPIns((char)(3), (char)(-1), CageList, Bundles);
                SetMSList(BondList_MLST_D1, BondList);
                // Rings
                TecGUIListDeleteAllItems(RingList_MLST_D1);
                RingList = ListOfCPOutsConnetedToCPIns((char)(3), (char)(1), CageList, Bundles);
                SetMSList(RingList_MLST_D1, RingList);

                // Set independent list to result and dealloc dependent lists
                ArrListDealloc(&BondList);
                ArrListDealloc(&AtomList);
                ArrListDealloc(&RingList);
                Result = CageList;
            }
            else
            {
                TecUtilDialogMessageBox("No Cages in topology", MessageBox_Warning);
                Result = ArrListAlloc(10, ArrListType_Long);
            }
            break;
        default:
            CHECK(FALSE);
    }

    ENSURE(ArrListIsValid(Result));
    return Result;
}






/*
 * Based on the ImageType and CPType options and the independent CPList (Bond, Atom, Ring,
 * or Cage) list selections, activate the appropriate Tecplot zones.
 */
void  ActivateZones(ArrList_pa CPList)
{
    Boolean_t GetEdges, GetSurfaces;
    LgIndex_t ActiveCPType = TecGUIOptionMenuGet(CPType_OPT_D1);
    LgIndex_t ImageType = TecGUIOptionMenuGet(ImageType_OPT_D1);

    switch (ImageType)
    {
        case 1:  // Show whichBond lines
        {
            ArrList_pa BondZones = NULL;

            switch (ActiveCPType)
            {
                case 1:  // Bond is independent list
                {
                    // Get the selected whichBond list items
                    ArrList_pa BondList = GetSelectedInMSList(BondList_MLST_D1);

                    // Determine the Bond lines zones corresponding to BondList.
                    BondZones = ZoneListOfBondLinesForCritPoint((char)(-1), BondList, CritPoints, Bundles);

                    // Clean-up
                    ArrListDealloc(&BondList);
                }
                break;

                case 2:  // Atom is independent list
                {
                    // Get the selected Atom list items
                    ArrList_pa AtomList = GetSelectedInMSList(AtomList_MLST_D1);

                    // Determine the Bond lines zones containing Atoms in AtomList.
                    BondZones = ZoneListOfBondLinesForCritPoint((char)(-3), AtomList, CritPoints, Bundles);

                    // Clean-up
                    ArrListDealloc(&AtomList);
                }
                break;

                case 3:  // Ring is independent list
                {
                    // Get the selected Ring list items
                    ArrList_pa RingList = GetSelectedInMSList(RingList_MLST_D1);

                    // Determine the Bond lines zones for bundles containing rings
                    // in RingList.
                    BondZones = ZoneListOfBondLinesForCritPoint((char)(1), RingList, CritPoints, Bundles);

                    // Clean-up
                    ArrListDealloc(&RingList);
                }
                break;

                case 4:  // Cage is independent list
                {
                    // Get the selected Cage list items
                    ArrList_pa CageList = GetSelectedInMSList(CageList_MLST_D1);

                    // Determine the Bond lines zones for bundles containing rings
                    // in RingList.
                    BondZones = ZoneListOfBondLinesForCritPoint((char)(3), CageList, CritPoints, Bundles);

                    // Clean-up
                    ArrListDealloc(&CageList);
                }
                break;
                default:
                    CHECK(FALSE);
            }

            // Deactivate the old whichBond zones
            if (ZoneListOld != NULL)
            {
                Set_pa ZoneSet = TecUtilSetAlloc(FALSE);
                int numEntries = ArrListGetCount(ZoneListOld);
                for (int entry = 0; entry < numEntries; entry++)
                {
                    ArrListItem_u Item = ArrListGetItem(ZoneListOld, entry);
                    TecUtilSetAddMember(ZoneSet, Item.Long, TRUE);
                }
                TecUtilZoneSetActive(ZoneSet, AssignOp_MinusEquals);
                TecUtilSetDealloc(&ZoneSet);

                ArrListDealloc(&ZoneListOld);
            }

            // Activate the whichBond zones
            if (BondZones != NULL)
            {
                Set_pa ZoneSet = TecUtilSetAlloc(FALSE);
                int numEntries = ArrListGetCount(BondZones);
                for (int entry = 0; entry < numEntries; entry++)
                {
                    ArrListItem_u Item = ArrListGetItem(BondZones, entry);
                    TecUtilSetAddMember(ZoneSet, Item.Long, TRUE);
                }
                TecUtilZoneSetActive(ZoneSet, AssignOp_PlusEquals);
                TecUtilSetDealloc(&ZoneSet);

                ZoneListOld = BondZones;
                // ArrListDealloc(&BondZones);
            }
        }
        break;
        case 2:  // Show all bundle critical point path lines
            GetEdges = TRUE;
            GetSurfaces = FALSE;
        case 3:  // Show all bundle surfaces
            if (ImageType == 3)
            {
                GetEdges = FALSE;
                GetSurfaces = TRUE;
            }
        case 4:  // Show all bundle critical point path lines and surfaces
            if (ImageType == 4)
            {
                GetEdges = TRUE;
                GetSurfaces = TRUE;
            }

            {
                ArrList_pa ZoneList = NULL;

                switch (ActiveCPType)
                {
                    case 1:  // Bond is independent list
                    {
                        // Get the selected whichBond list items
                        ArrList_pa BondList = GetSelectedInMSList(BondList_MLST_D1);

                        // Determine the Bond lines zones corresponding to BondList.
                        ZoneList = ZoneListOfBundleZonesForCritPoint(GetEdges, GetSurfaces,
                                                                     (char)(-1), BondList, CritPoints, Bundles);

                        // Clean-up
                        ArrListDealloc(&BondList);
                    }
                    break;

                    case 2:  // Atom is independent list
                    {
                        // Get the selected Atom list items
                        ArrList_pa AtomList = GetSelectedInMSList(AtomList_MLST_D1);

                        // Determine the Bond lines zones containing Atoms in AtomList.
                        ZoneList = ZoneListOfBundleZonesForCritPoint(GetEdges, GetSurfaces,
                                                                     (char)(-3), AtomList, CritPoints, Bundles);

                        // Clean-up
                        ArrListDealloc(&AtomList);
                    }
                    break;

                    case 3:  // Ring is independent list
                    {
                        // Get the selected Ring list items
                        ArrList_pa RingList = GetSelectedInMSList(RingList_MLST_D1);

                        // Determine the Bond lines zones for bundles containing rings
                        // in RingList.
                        ZoneList = ZoneListOfBundleZonesForCritPoint(GetEdges, GetSurfaces,
                                                                     (char)(1), RingList, CritPoints, Bundles);

                        // Clean-up
                        ArrListDealloc(&RingList);
                    }
                    break;

                    case 4:  // Cage is independent list
                    {
                        // Get the selected Cage list items
                        ArrList_pa CageList = GetSelectedInMSList(CageList_MLST_D1);

                        // Determine the Bond lines zones for bundles containing rings
                        // in RingList.
                        ZoneList = ZoneListOfBundleZonesForCritPoint(GetEdges, GetSurfaces,
                                                                     (char)(3), CageList, CritPoints, Bundles);

                        // Clean-up
                        ArrListDealloc(&CageList);
                    }
                    break;
                    default:
                        CHECK(FALSE);
                }

                // Deactivate the old whichBond zones
                if (ZoneListOld != NULL)
                {
                    Set_pa ZoneSet = TecUtilSetAlloc(FALSE);
                    int numEntries = ArrListGetCount(ZoneListOld);
                    for (int entry = 0; entry < numEntries; entry++)
                    {
                        ArrListItem_u Item = ArrListGetItem(ZoneListOld, entry);
                        TecUtilSetAddMember(ZoneSet, Item.Long, TRUE);
                    }
                    TecUtilZoneSetActive(ZoneSet, AssignOp_MinusEquals);
                    TecUtilSetDealloc(&ZoneSet);

                    ArrListDealloc(&ZoneListOld);
                }

                // Activate the whichBond zones
                if (ZoneList != NULL)
                {
                    Set_pa ZoneSet = TecUtilSetAlloc(FALSE);
                    int numEntries = ArrListGetCount(ZoneList);
                    for (int entry = 0; entry < numEntries; entry++)
                    {
                        ArrListItem_u Item = ArrListGetItem(ZoneList, entry);
                        TecUtilSetAddMember(ZoneSet, Item.Long, TRUE);
                    }
                    TecUtilZoneSetActive(ZoneSet, AssignOp_PlusEquals);
                    TecUtilSetDealloc(&ZoneSet);

                    ZoneListOld = ZoneList;
                }
            }
            break;
        default:
            break;
    }


}



/**
 */
static void Dialog1HelpButton_CB(void)
{
    TecUtilLockStart(AddOnID);
    TecUtilDialogMessageBox("On-line Help not available for this dialog.",
                            MessageBox_Information);
    TecUtilLockFinish(AddOnID);
}


/**
 */
static void Dialog1CloseButton_CB(void)
{
    TecUtilLockStart(AddOnID);
    TecGUIDialogDrop(Dialog1Manager);
    TecUtilLockFinish(AddOnID);
}


/**
 */
static void Dialog1Init_CB(void)
{
    Boolean_t IsOk = TRUE;

    TecUtilLockStart(AddOnID);

    // Extract the topology from the dataset and entry auxiliary data
    IsOk = ExtractTopoInZonesAuxData(&Bundles, &CritPoints);

    // Load Multi-Selection lists with defaults (all critical points)
    IsOk = InitFullMSLists();

    if (IsOk)
    {
        ArrList_pa BondList = NULL;

        // Determine and set the connecting Atoms, Rings, and Cages.
        BondList = SetListItems();

        // Activate all appropriate zones
        ActivateZones(BondList);
        ArrListDealloc(&BondList);

        // Set List sensitivities
        SetListSensitivity();
    }

    TecUtilLockFinish(AddOnID);
}


char *SourceZone_OPT_D1_List = "1";



/**
 */
static void SourceZone_OPT_D1_CB(const LgIndex_t *I)
{
    TecUtilLockStart(AddOnID);
    TRACE1("Option Menu (SourceZone_OPT_D1) value changed,  New value is: %d\n", *I);
    TecUtilLockFinish(AddOnID);
}


/**
 */
static void ComputeTopo_BTN_D1_CB(void)
{
    TecUtilLockStart(AddOnID);
    TRACE("Extract Topology Button Pushed\n");
    TecUtilLockFinish(AddOnID);
}


char *CPType_OPT_D1_List = "Bond,Atom,Ring,Cage";



/**
 */
static void CPType_OPT_D1_CB(const LgIndex_t *I)
{
    Boolean_t IsOk = TRUE;

    TecUtilLockStart(AddOnID);
    TRACE1("Option Menu (CPType_OPT_D1) value changed,  New value is: %d\n", *I);

    // Load Multi-Selection lists with defaults (all critical points)
    IsOk = InitFullMSLists();

    SetListSensitivity();

    if (IsOk)
    {
        ArrList_pa BondList = NULL;

        // Determine and set the connecting Atoms, Rings, and Cages.
        BondList = SetListItems();

        // Activate all appropriate zones
        ActivateZones(BondList);
        ArrListDealloc(&BondList);
    }



    TecUtilLockFinish(AddOnID);
}


/**
 */
static void BondList_MLST_D1_CB(const LgIndex_t *I)
{
    ArrList_pa BondList = NULL;

    TecUtilLockStart(AddOnID);
    TRACE1("Multi selection list (BondList_MLST_D1) item selected,  First Item is: %d\n", *I);

    // Determine and set the connecting Atoms, Rings, and Cages.
    BondList = SetListItems();

    // Activate all appropriate zones
    ActivateZones(BondList);
    ArrListDealloc(&BondList);

    // Set List sensitivities
    SetListSensitivity();

    TecUtilLockFinish(AddOnID);
}


/**
 */
static void AtomList_MLST_D1_CB(const LgIndex_t *I)
{
    ArrList_pa AtomList = NULL;

    TecUtilLockStart(AddOnID);
    TRACE1("Multi selection list (AtomList_MLST_D1) item selected,  First Item is: %d\n", *I);

    // Determine and set the connecting Bonds, Rings, and Cages.
    AtomList = SetListItems();

    // Activate all appropriate zones
    ActivateZones(AtomList);
    ArrListDealloc(&AtomList);

    // Set List sensitivities
    SetListSensitivity();

    TecUtilLockFinish(AddOnID);
}


/**
 */
static void RingList_MLST_D1_CB(const LgIndex_t *I)
{
    ArrList_pa RingList = NULL;

    TecUtilLockStart(AddOnID);
    TRACE1("Multi selection list (RingList_MLST_D1) item selected,  First Item is: %d\n", *I);

    // Determine and set the connecting Atoms, Bonds, and Cages.
    RingList = SetListItems();

    // Activate all appropriate
    ActivateZones(RingList);
    ArrListDealloc(&RingList);

    // Set List sensitivities
    SetListSensitivity();

    TecUtilLockFinish(AddOnID);
}


/**
 */
static void CageList_MLST_D1_CB(const LgIndex_t *I)
{
    ArrList_pa CageList = NULL;

    TecUtilLockStart(AddOnID);
    TRACE1("Multi selection list (CageList_MLST_D1) item selected,  First Item is: %d\n", *I);

    // Determine and set the connecting Atoms, Bonds, and Ringss.
    CageList = SetListItems();

    // Activate all appropriate
    ActivateZones(CageList);
    ArrListDealloc(&CageList);

    // Set List sensitivities
    SetListSensitivity();

    TecUtilLockFinish(AddOnID);
}


char *ImageType_OPT_D1_List = "Bonds,All Lines,Surfaces,Lines and Surface";



/**
 */
static void ImageType_OPT_D1_CB(const LgIndex_t *I)
{
    ArrList_pa CPList = NULL;

    TecUtilLockStart(AddOnID);
    TRACE1("Option Menu (ImageType_OPT_D1) value changed,  New value is: %d\n", *I);

    // Determine and set the connecting Atoms, Bonds, and Ringss.
    CPList = SetListItems();

    // Activate all appropriate
    ActivateZones(CPList);
    ArrListDealloc(&CPList);

    // Set List sensitivities
    SetListSensitivity();

    TecUtilLockFinish(AddOnID);
}


/**
 */
static void ShowNeighbor_TOG_D1_CB(const LgIndex_t *I)
{
    TecUtilLockStart(AddOnID);
    TRACE1("Toggle (ShowNeighbor_TOG_D1) value Changed,  New value is: %d\n", *I);
    TecUtilLockFinish(AddOnID);
}

void setDisplayMode(DisplayMode_e displayMode)
{
    REQUIRE(VALID_ENUM(displayMode, DisplayMode_e));
    if ( displayMode != curDisplayMode )
    {
        curDisplayMode = displayMode;
        updateBondOrAtomList(curDisplayMode);
        updateInfoField();
        updateAtomAndBondZones(curAtom, curBond, curDisplayMode);
    }
    ENSURE(VALID_ENUM(curDisplayMode, DisplayMode_e));
}

#include "guibld.cpp"
