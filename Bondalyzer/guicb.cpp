/*
*****************************************************************
*****************************************************************
*******                                                  ********
*******    (C) Copyright 1988-2013 by TECPLOT INC.       *******
*******                                                  ********
*****************************************************************
*****************************************************************
*/

#include <vector>
#include <string>
#include "TECADDON.h"
#include "ADDGLBL.h"
#include "GUIDEFS.h"
#include "ENGINE.h"
#include "SETCHEMBONDSTYLE.h"
#include "ARRLIST.h"
#include "SURFELEMMAP.h"
#include "NORMALS.h"
#include "ZONEVARINFO.h"


Boolean_t selectDefaultRefinedGridsForSourceZone(EntIndex_t sourceZone) 
{
    EntIndex_t numZones;
    Boolean_t  isOk = TecUtilDataSetGetInfo(NULL, &numZones, NULL);
    if (isOk)
    {
        TecGUIListDeselectAllItems(RefinedGrid_MLST_D1);
        EntIndex_t zoneEntry = 0;
        for (EntIndex_t zoneNum = 1; zoneNum <= numZones; zoneNum++)
            if (TecUtilZoneIsEnabled(zoneNum))
            {
                zoneEntry++;
                if ( zoneNum != sourceZone )
                    TecGUIListSetSelectedItem(RefinedGrid_MLST_D1, zoneEntry);
            }
    }
    return isOk;
}

Boolean_t selectZoneInOptionMenu(LgIndex_t  optionMenu,
                                 EntIndex_t zoneToSelect) 
{
    EntIndex_t numZones;
    Boolean_t  isOk = TecUtilDataSetGetInfo(NULL, &numZones, NULL);
    if (isOk)
    {
        isOk = FALSE; // not found yet
        EntIndex_t zoneEntry = 0;
        for (EntIndex_t zoneNum = 1; zoneNum <= numZones; zoneNum++)
        {
            if (TecUtilZoneIsEnabled(zoneNum))
            {
                zoneEntry++;
                if ( zoneNum == zoneToSelect )
                {
                    TecGUIOptionMenuSet(optionMenu, zoneEntry);
                    isOk = TRUE; // found
                    break;
                }
            }
        }
    }
    return isOk;
}

static EntIndex_t getZoneNumFromGuiOffset(LgIndex_t offsetInOptionMenu/*1-based*/)
{
    EntIndex_t result = 0;/*1-based*/
    EntIndex_t numZones;
    if ( TecUtilDataSetGetInfo(NULL, &numZones, NULL) )
    {
        EntIndex_t zoneEntry = 0;
        for (EntIndex_t zoneNum = 1; zoneNum <= numZones; zoneNum++)
        {
            if (TecUtilZoneIsEnabled(zoneNum))
            {
                zoneEntry++;
                if ( zoneEntry == offsetInOptionMenu )
                {
                    result = zoneNum;
                    break;
                }
            }
        }
    }
    return result;
}

Boolean_t selectVarInOptionMenu(LgIndex_t  optionMenu,
                                EntIndex_t varToSelect) 
{
    EntIndex_t numVars;
    Boolean_t  isOk = TecUtilDataSetGetInfo(NULL, NULL, &numVars);
    if (isOk)
    {
        isOk = FALSE; // not found yet
        EntIndex_t varEntry = 0;
        for (EntIndex_t varNum = 1; varNum <= numVars; varNum++)
        {
            if (TecUtilVarIsEnabled(varNum))
            {
                varEntry++;
                if ( varNum == varToSelect )
                {
                    TecGUIOptionMenuSet(optionMenu, varEntry);
                    isOk = TRUE; // found
                    break;
                }
            }
        }
    }
    return isOk;
}

EntIndex_t getVarNumFromGuiOffset(LgIndex_t offsetInOptionMenu/*1-based*/)
{
    EntIndex_t result = 0;/*1-based*/
    EntIndex_t numVars;
    if ( TecUtilDataSetGetInfo(NULL, NULL, &numVars) )
    {
        EntIndex_t varEntry = 0;
        for (EntIndex_t varNum = 1; varNum <= numVars; varNum++)
        {
            if (TecUtilVarIsEnabled(varNum))
            {
                varEntry++;
                if ( varEntry == offsetInOptionMenu )
                {
                    result = varNum;
                    break;
                }
            }
        }
    }
    return result;
}

/**
 */
static void Dialog1Init_CB(void)
{

    TecUtilLockStart(AddOnID);

    EntIndex_t NumZones, NumVars;
    Boolean_t  IsOk = TecUtilDataSetGetInfo(NULL, &NumZones, &NumVars);

    /* Load Tensor component variable option lists. */

    TecGUIOptionMenuDeleteAllItems(SourceZone_OPT_D1);
    TecGUIListDeleteAllItems(RefinedGrid_MLST_D1);
    TecGUIOptionMenuDeleteAllItems(ScalarVarNum_OPT_D1);
    TecGUIOptionMenuDeleteAllItems(GradXVar_OPT_D1);
    TecGUIOptionMenuDeleteAllItems(GradYVar_OPT_D1);
    TecGUIOptionMenuDeleteAllItems(GradZVar_OPT_D1);
    TecGUIOptionMenuDeleteAllItems(GradTot_OPT_D1);
    TecGUIOptionMenuDeleteAllItems(ComputeOptio_OPT_D1);

    EntIndex_t defaultSourceZone = 0;
    for (EntIndex_t zoneNum = 1; zoneNum <= NumZones && IsOk; zoneNum++)
    {
        if (TecUtilZoneIsEnabled(zoneNum))
        {
            char *ZoneName = NULL;
            IsOk = TecUtilZoneGetName(zoneNum, &ZoneName);
            char formattedName[400];
            sprintf_s(formattedName, "%d: %s", zoneNum, ZoneName);
            TecGUIOptionMenuAppendItem(SourceZone_OPT_D1, formattedName);
            TecGUIListAppendItem(RefinedGrid_MLST_D1, formattedName);
            TecUtilStringDealloc(&ZoneName);
            defaultSourceZone = zoneNum;
        }
    }

	defaultSourceZone = ZoneNumByName(std::string("Full Volume"));
	if (defaultSourceZone <= 0)
		defaultSourceZone = 1;

    IsOk = IsOk && selectZoneInOptionMenu(SourceZone_OPT_D1, defaultSourceZone);
    IsOk = IsOk && selectDefaultRefinedGridsForSourceZone(defaultSourceZone);

    for (EntIndex_t VarNum = 1; VarNum <= NumVars; VarNum++)
    {
        char *VarName = NULL;
        if (TecUtilVarIsEnabled(VarNum))
        {
            TecUtilVarGetName(VarNum, &VarName);
            TecGUIOptionMenuAppendItem(ScalarVarNum_OPT_D1, VarName);
            TecGUIOptionMenuAppendItem(GradXVar_OPT_D1, VarName);
            TecGUIOptionMenuAppendItem(GradYVar_OPT_D1, VarName);
            TecGUIOptionMenuAppendItem(GradZVar_OPT_D1, VarName);
            TecGUIOptionMenuAppendItem(GradTot_OPT_D1, VarName);
            TecUtilStringDealloc(&VarName);
        }
    }

	std::vector<std::string> VarNameList;
	EntIndex_t VarNum;

	VarNameList.push_back("Rho");
	VarNameList.push_back("Electron Density");
	VarNum = VarNumByNameList(VarNameList);
	if (VarNum < 0) VarNum = 4;
    selectVarInOptionMenu(ScalarVarNum_OPT_D1, VarNum);

	VarNum = VarNumByName(std::string("X Density Gradient"));
	if (VarNum < 0) VarNum = 5;
	selectVarInOptionMenu(GradXVar_OPT_D1, VarNum);

	VarNum = VarNumByName(std::string("Y Density Gradient"));
	if (VarNum < 0) VarNum = 6;
	selectVarInOptionMenu(GradYVar_OPT_D1, VarNum);

	VarNum = VarNumByName(std::string("Z Density Gradient"));
	if (VarNum < 0) VarNum = 7;
	selectVarInOptionMenu(GradZVar_OPT_D1, VarNum);

	VarNum = VarNumByName(std::string("Density Gradient Magnitude"));
	if (VarNum < 0) VarNum = 5;
	selectVarInOptionMenu(GradTot_OPT_D1, VarNum);
	VarNameList.clear();

// 	if (NumVars >= 4)
// 		selectVarInOptionMenu(ScalarVarNum_OPT_D1, 4);
// 	if (NumVars >= 5)
// 		selectVarInOptionMenu(GradXVar_OPT_D1, 5);
// 	if (NumVars >= 6)
// 		selectVarInOptionMenu(GradYVar_OPT_D1, 6);
// 	if (NumVars >= 7)
// 		selectVarInOptionMenu(GradZVar_OPT_D1, 7);
// 	if (NumVars >= 8)
// 		selectVarInOptionMenu(GradTot_OPT_D1, 8);

    TecGUIOptionMenuAppendItem(ComputeOptio_OPT_D1, "Critical Points");
    TecGUIOptionMenuAppendItem(ComputeOptio_OPT_D1, "Bond Lines");
    TecGUIOptionMenuAppendItem(ComputeOptio_OPT_D1, "Ring-Cage Lines");
    TecGUIOptionMenuAppendItem(ComputeOptio_OPT_D1, "Bond-Ring-Cage Surfaces");
    TecGUIOptionMenuAppendItem(ComputeOptio_OPT_D1, "Ring-Bond-Atom Surfaces");
    TecGUIOptionMenuAppendItem(ComputeOptio_OPT_D1, "Atom Iso-surface Topology");
    TecGUIOptionMenuAppendItem(ComputeOptio_OPT_D1, "Atom-Ring-Cage Surfaces");
    TecGUIOptionMenuAppendItem(ComputeOptio_OPT_D1, "Bond-Bundle Volumes");
    TecGUIOptionMenuAppendItem(ComputeOptio_OPT_D1, "All");
    TecGUIOptionMenuSet(ComputeOptio_OPT_D1, 9);

    

    TecGUIToggleSet(PeriodicBC_TOG_D1, FALSE);

}

/**
 */
static void Dialog1OkButton_CB(void)
{
    Boolean_t isOk = TRUE;

    EntIndex_t sourceZoneOffset = TecGUIOptionMenuGet(SourceZone_OPT_D1);
    EntIndex_t sourceZone = getZoneNumFromGuiOffset(sourceZoneOffset);

    Set_pa refinedGridSet = TecUtilSetAlloc(TRUE);
    isOk = (refinedGridSet != NULL);
    if (isOk)
    {
        LgIndex_t count;
        LgIndex_t *seletedItems;
        TecGUIListGetSelectedItems(RefinedGrid_MLST_D1, &seletedItems, &count);
        for (LgIndex_t ii=0; ii<count && isOk; ii++)
        {
            EntIndex_t zoneNum = getZoneNumFromGuiOffset(seletedItems[ii]);
            if ( zoneNum != sourceZone )
                isOk = TecUtilSetAddMember(refinedGridSet, zoneNum, TRUE);
            else
            {
                TecUtilDialogErrMsg("The source varNum also in the refined grids.");
                isOk = FALSE;
            }
        }
        TecUtilArrayDealloc((void **)&seletedItems);
    }


	// Selecting which scalar property to look at
    EntIndex_t ScalarVarNum = getVarNumFromGuiOffset(TecGUIOptionMenuGet(ScalarVarNum_OPT_D1));

	// Selecting gradient of which property in x,y,z
    EntIndex_t GradXVarNum  = getVarNumFromGuiOffset(TecGUIOptionMenuGet(GradXVar_OPT_D1));
    EntIndex_t GradYVarNum  = getVarNumFromGuiOffset(TecGUIOptionMenuGet(GradYVar_OPT_D1));
    EntIndex_t GradZVarNum  = getVarNumFromGuiOffset(TecGUIOptionMenuGet(GradZVar_OPT_D1));
    EntIndex_t GradTotNum   = getVarNumFromGuiOffset(TecGUIOptionMenuGet(GradTot_OPT_D1));

    CompletionLevel_e CompletionLevel = (CompletionLevel_e)(TecGUIOptionMenuGet(ComputeOptio_OPT_D1) - 1);

    Boolean_t  PeriodicBC = TecGUIToggleGet(PeriodicBC_TOG_D1);

    if (isOk)
    {
        TecGUIDialogDrop(Dialog1Manager);

        // This is a function in 360, and returns the gradient of the system.
		// Here, gradient refers not to the general idea of "the change of a property with spacial displacment",
		// But rather the 1st derivative at each point in the system w.r.t. rho (along x,y,z) by the centered divided difference
		ExtractTopology(sourceZone, refinedGridSet,
                        GradXVarNum, GradYVarNum, GradZVarNum, GradTotNum,
                        ScalarVarNum,
                        CompletionLevel,
                        PeriodicBC);
		// Zone can be 1, 2, 3D, or a collection of sparse points
    }

    if ( refinedGridSet != NULL )
        TecUtilSetDealloc(&refinedGridSet);

	ChemSysView();

    /* Only unlock tecplot here because a modal dialog was launched. */
    if ( isOk )
        TecUtilLockFinish(AddOnID);
}

/**
 */
static void Dialog1CancelButton_CB(void)
{
    /* Only unlock tecplot here because a modal dialog was launched. */
    TecGUIDialogDrop(Dialog1Manager);
    TecUtilLockFinish(AddOnID);
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


char *SourceZone_OPT_D1_List = "Option 1,Option 2,Option 3";



/**
 */
static void SourceZone_OPT_D1_CB(const LgIndex_t *I)
{
    TecUtilLockStart(AddOnID);
    TRACE1("Option Menu (SourceZone_OPT_D1) value changed,  New value is: %d\n", *I);
    EntIndex_t sourceZone = getZoneNumFromGuiOffset(*I);
    if (sourceZone > 0)
        selectDefaultRefinedGridsForSourceZone(sourceZone);
    TecUtilLockFinish(AddOnID);
}


/**
*/
static void RefinedGrid_MLST_D1_CB(const LgIndex_t *I)
{
    TecUtilLockStart(AddOnID);
    TRACE1("Multi selection list (RefinedGrid_MLST_D1) item selected,  First Item is: %d\n",*I);
    TecUtilLockFinish(AddOnID);
}



char *ScalarVarNum_OPT_D1_List = "Option 1,Option 2,Option 3";



/**
 */
static void ScalarVarNum_OPT_D1_CB(const LgIndex_t *I)
{
    TecUtilLockStart(AddOnID);
    TRACE1("Option Menu (ScalarVarNum_OPT_D1) value changed,  New value is: %d\n", *I);
    TecUtilLockFinish(AddOnID);
}


char *GradXVar_OPT_D1_List = "Option 1,Option 2,Option 3";



/**
 */
static void GradXVar_OPT_D1_CB(const LgIndex_t *I)
{
    TecUtilLockStart(AddOnID);
    TRACE1("Option Menu (GradXVar_OPT_D1) value changed,  New value is: %d\n", *I);
    TecUtilLockFinish(AddOnID);
}


char *GradYVar_OPT_D1_List = "Option 1,Option 2,Option 3";



/**
 */
static void GradYVar_OPT_D1_CB(const LgIndex_t *I)
{
    TecUtilLockStart(AddOnID);
    TRACE1("Option Menu (GradYVar_OPT_D1) value changed,  New value is: %d\n", *I);
    TecUtilLockFinish(AddOnID);
}


char *GradZVar_OPT_D1_List = "Option 1,Option 2,Option 3";



/**
 */
static void GradZVar_OPT_D1_CB(const LgIndex_t *I)
{
    TecUtilLockStart(AddOnID);
    TRACE1("Option Menu (GradZVar_OPT_D1) value changed,  New value is: %d\n", *I);
    TecUtilLockFinish(AddOnID);
}


char *GradTot_OPT_D1_List = "Option 1,Option 2,Option 3";



/**
 */
static void GradTot_OPT_D1_CB(const LgIndex_t *I)
{
  TecUtilLockStart(AddOnID);
  TRACE1("Option Menu (GradTot_OPT_D1) value changed,  New value is: %d\n",*I);
  TecUtilLockFinish(AddOnID);
}


char *ComputeOptio_OPT_D1_List = "Option 1,Option 2,Option 3";



/**
 */
static void ComputeOptio_OPT_D1_CB(const LgIndex_t *I)
{
  TecUtilLockStart(AddOnID);
  TRACE1("Option Menu (ComputeOptio_OPT_D1) value changed,  New value is: %d\n",*I);
  TecUtilLockFinish(AddOnID);
}



/**
 */
static void PeriodicBC_TOG_D1_CB(const LgIndex_t *I)
{
    TecUtilLockStart(AddOnID);
    TRACE1("Toggle (PeriodicBC_TOG_D1) Value Changed,  New value is: %d\n", *I);
    TecUtilLockFinish(AddOnID);
}



#include "guibld.cpp"

