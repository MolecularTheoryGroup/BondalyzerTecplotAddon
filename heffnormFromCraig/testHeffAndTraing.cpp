#include "TECADDON.h"
#define _USE_MATH_DEFINES
#include <math.h>
#include <string.h>
#include <string>
#include <sstream>
#include <algorithm> // for std::min
#include <ArgList.h>

#include "ADDONVER.h"
#include "ArgList.h"
#include "StyleValue.h"
#include "Set.h"
#include "StringList.h"
#include "heffNorm.h"
#include "UpdateSphericalTriangulation.h"

using namespace tecplot::toolbox;

AddOn_pa AddOnID;

static char const* UPDATE_SPHERICAL_TRIANGULATION_COMMAND = "Op=UpdateSpherical";

/**
 */
static void STDCALL StateChangeCallback(StateChange_e StateChange)
{
    switch (StateChange)
    {
        case StateChange_VarsAltered:     /* set of altered variables */
        case StateChange_VarsAdded:       /* set of added variables */
        case StateChange_ZonesDeleted:    /* set of deleted zones */
        case StateChange_ZonesAdded:      /* set of added zones */
        case StateChange_NodeMapsAltered: /* set of node maps altered */
        case StateChange_MouseModeUpdate: /* the new mouse mode */
        case StateChange_Style:           /* Style Parameters P1,P2,P3,P4,P5,P6 */
        case StateChange_View:            /* View action (View_e) */
        case StateChange_Streamtrace:     /* Streamtrace action (Streamtrace_e) */
        case StateChange_AuxDataAltered:  /* Name, Auxiliary Location (AuxDataLocation_e), Var/Zone/or Map Num */
        case StateChange_AuxDataAdded:    /* Name, Auxiliary Location (AuxDataLocation_e), Var/Zone/or Map Num */
        case StateChange_AuxDataDeleted:  /* Name, Auxiliary Location (AuxDataLocation_e), Var/Zone/or Map Num */
        case StateChange_VarsDeleted:     /* set of deleted variables (zero based set) */
        case StateChange_VariableLockOn:  /* Locker name, Variable Num, VarLockMode */
        case StateChange_VariableLockOff: /* Unlocker name, Variable Num */
        case StateChange_DataSetLockOn:   /* Locker name */
        case StateChange_DataSetLockOff:  /* Unlocker name */
        case StateChange_TecplotIsInitialized:/* Tecplot is finished initializing */
        case StateChange_FrameDeleted:        /* A frame was delete */
        case StateChange_NewTopFrame:         /* A new frame has become the current frame */
        case StateChange_Text:                /* One or more text elements has changed */
        case StateChange_Geom:                /* One or more geometry elements has changed */
        case StateChange_DataSetReset:        /* A new dataset has been loaded */
        case StateChange_NewLayout:           /* The current layout has been cleared and reset */
        case StateChange_CompleteReset:       /* Anything could have happened */
        case StateChange_LineMapAssignment:   /* A line mapping definition has been altered (includes zone and axis information) */
        case StateChange_ContourLevels:       /* The contour levels have been altered */
        case StateChange_ModalDialogLaunch:   /* A modal dialog has been launched */
        case StateChange_ModalDialogDismiss:  /* A modal dialog has been dismissed */
        case StateChange_QuitTecplot:         /* Tecplot is about to exit */
        case StateChange_ZoneName:            /* The name of a zone has been altered */
        case StateChange_VarName:             /* The name of a variable has been altered */
        case StateChange_LineMapName:           /* The name of an X-Y mapping has been altered */
        case StateChange_LineMapAddDeleteOrReorder: /* The set of existing X-Y mappings has been altered */
        case StateChange_ColorMap:            /* The color mapping has been altered */
        case StateChange_ContourVar:          /* The contour variable has been reassigned */
        case StateChange_NewAxisVariables:    /* The axis variables have been reassigned */
        case StateChange_PickListCleared:     /* All picked objects are unpicked */
        case StateChange_PickListGroupSelect: /* A group of objects has been added to the pick list */
        case StateChange_PickListSingleSelect:/* A single object has been added to or removed from the pick list */
        case StateChange_PickListStyle:       /* An action has been performed on all of the objects in the pick list */
        case StateChange_DataSetFileName:     /* The current data set has been saved to a file */
        case StateChange_DataSetTitle:        /* The current data set title has been changed */
        case StateChange_DrawingInterrupted:  /* The user has interrupted the drawing */
        case StateChange_ImageExported:       /* An image frame was exported */
        case StateChange_PageDeleted:         /* A page was deleted */
        case StateChange_NewTopPage:          /* A different page was made the top page */
        case StateChange_PrintPreviewLaunch:  /* Modeless dialogs should close or disable themselves */
        case StateChange_PrintPreviewDismiss: /* Modeless dialogs can re-launch or enable themselves */
        case StateChange_SuspendInterface:    /* Replaces StateChange_DrawGraphicsOn */
        case StateChange_UnsuspendInterface:  /* Replaces StateChange_DrawGraphicsOff */
        {
            /* TODO: Add code to handle state changes.... */
        } break;
        default: break;
    } /* end switch */
}

/**
*/
char const* performUpdateSphericalTriangulation(
    Index_t& numNewNodes,
    Index_t& numNewTriangles,
    Index_t& numCoincidentNodesIgnored)
{
    using tpcsm::Vec3;

    const char* errMsg = NULL;

    numNewNodes = 0;
    numNewTriangles = 0;
    numCoincidentNodesIgnored = 0;

    if (TecUtilFrameGetPlotType() != PlotType_Cartesian3D)
        errMsg = "Update spherical triangulation test requires a 3-D plot.";
    
    Set_pa activeZones = NULL;
    if ( errMsg == NULL && !TecUtilZoneGetActive(&activeZones))
        errMsg = "Update spherical triangulation encountered memory error allocating active zone set.";

    Set_pa constraintZones = NULL;
    if (errMsg == NULL)
    {
        constraintZones = TecUtilSetAlloc(FALSE);
        if (constraintZones == NULL)
            errMsg = "Update spherical triangulation encountered memory error allocating constraint zone set.";
    }

    EntIndex_t srcSphericalZoneNum = 0;
    if (errMsg == NULL)
    {
        EntIndex_t xVarNum;
        EntIndex_t yVarNum;
        EntIndex_t zVarNum;
        TecUtilAxisGetVarAssignments(&xVarNum, &yVarNum, &zVarNum);

        EntIndex_t numSrcZones = TecUtilDataSetGetNumZones();
        for (EntIndex_t zone = 0; zone < numSrcZones && errMsg == NULL; zone++)
        {
            if (TecUtilZoneIsEnabled(zone+1) && TecUtilZoneIsActive(zone+1))
            {
                // get info about each enabled zone
                LgIndex_t zoneType = TecUtilZoneGetType(zone+1);
                CHECK(VALID_ENUM(zoneType, ZoneType_e));

                LgIndex_t srcI = 0;
                LgIndex_t srcJ = 0;
                LgIndex_t srcK = 0;
                TecUtilZoneGetIJK(zone+1, &srcI, &srcJ, &srcK);
                CHECK(srcI > 0 && srcJ > 0 && srcK > 0);

                if (zoneType == ZoneType_FETriangle ||
                    (zoneType == ZoneType_Ordered && srcI > 1 && srcJ > 1 && srcK == 1))
                {
                    if (srcSphericalZoneNum == 0)
                        srcSphericalZoneNum = zone + 1;
                    else
                        errMsg = "Update spherical triangulation detected multiple spherical surface zones, either FE-triangle or IJK-ordered.";
                }
                else if (zoneType == ZoneType_Ordered && srcI > 1 && srcJ == 1 && srcK == 1)
                {
                    if (!TecUtilSetAddMember(constraintZones, zone + 1, FALSE))
                        errMsg = "Update spherical triangulation encountered memory error filling constraint zone set.";
                }

                if (errMsg != NULL)
                {
                    FieldData_pa XD = TecUtilDataValueGetReadableRef(zone+1, xVarNum);
                    FieldData_pa YD = TecUtilDataValueGetReadableRef(zone+1, yVarNum);
                    FieldData_pa ZD = TecUtilDataValueGetReadableRef(zone+1, zVarNum);

                    if (TecUtilDataValueGetLocationByRef(XD) != ValueLocation_Nodal ||
                        TecUtilDataValueGetLocationByRef(YD) != ValueLocation_Nodal ||
                        TecUtilDataValueGetLocationByRef(ZD) != ValueLocation_Nodal)
                    {
                        errMsg = "Update spherical triangulation requires X, Y and Z variables be nodal for spherical surface and constraint zones.";
                    }
                }
            }
        }
    }

    if (errMsg == NULL)
    {
        if (srcSphericalZoneNum == 0)
            errMsg = "Update spherical triangulation requires a single spherical surface zone, either FE-triangle or IJ-ordered.";
        if (TecUtilSetIsEmpty(constraintZones))
            errMsg = "Update spherical triangulation found no I-ordered zones to use as constraints.";
    }

    if (errMsg == NULL)
    {
        char message[2000];
        _snprintf(message, sizeof(message), "Spherical zone %d found along with %d constraint zone(s) (%d,...,%d).",
                  srcSphericalZoneNum+1,
                  TecUtilSetGetMemberCount(constraintZones),
                  TecUtilSetGetNextMember(constraintZones,TECUTILSETNOTMEMBER),
                  TecUtilSetGetPrevMember(constraintZones, TECUTILSETNOTMEMBER));
        TecUtilDialogMessageBox(message, MessageBoxType_Information);
    }

    if (errMsg == NULL)
    {
        LgIndex_t srcI = 0;
        LgIndex_t srcJ = 0;
        LgIndex_t srcK = 0;
        TecUtilZoneGetIJK(srcSphericalZoneNum, &srcI, &srcJ, &srcK);
        CHECK(srcI > 0 && srcJ > 0 && srcK > 0);

        LgIndex_t zoneType = TecUtilZoneGetType(srcSphericalZoneNum);
        CHECK(VALID_ENUM(zoneType, ZoneType_e));

        EntIndex_t xVarNum;
        EntIndex_t yVarNum;
        EntIndex_t zVarNum;
        TecUtilAxisGetVarAssignments(&xVarNum, &yVarNum, &zVarNum);

        FieldData_pa XD = TecUtilDataValueGetReadableRef(srcSphericalZoneNum, xVarNum);
        FieldData_pa YD = TecUtilDataValueGetReadableRef(srcSphericalZoneNum, yVarNum);
        FieldData_pa ZD = TecUtilDataValueGetReadableRef(srcSphericalZoneNum, zVarNum);

        CHECK(zoneType == ZoneType_FETriangle || (zoneType == ZoneType_Ordered && srcI > 1 && srcJ > 1 && srcK == 1)); // checked above in loop
        CHECK(TecUtilDataValueGetLocationByRef(XD) == ValueLocation_Nodal);
        CHECK(TecUtilDataValueGetLocationByRef(YD) == ValueLocation_Nodal);
        CHECK(TecUtilDataValueGetLocationByRef(ZD) == ValueLocation_Nodal);

        LgIndex_t const numSrcNodes = TecUtilDataValueGetCountByRef(XD);
        CHECK(numSrcNodes == TecUtilDataValueGetCountByRef(YD));

        std::vector<Vec3> startingNodes;
        startingNodes.reserve(numSrcNodes);
        Vec3 minExtents(TecUtilDataValueGetByRef(XD, 1),
                        TecUtilDataValueGetByRef(YD, 1),
                        TecUtilDataValueGetByRef(ZD, 1));
        Vec3 maxExtents = minExtents;
        for (LgIndex_t srcNodeNum = 1; srcNodeNum <= numSrcNodes; srcNodeNum++)
        {
            Vec3 pt(TecUtilDataValueGetByRef(XD, srcNodeNum),
                    TecUtilDataValueGetByRef(YD, srcNodeNum),
                    TecUtilDataValueGetByRef(ZD, srcNodeNum));
            startingNodes.push_back(pt);
            minExtents = minExtents.min(pt);
            maxExtents = maxExtents.max(pt);
        }
        Vec3 const sphereCenter = (minExtents + maxExtents) / 2.0; // better than the average of all nodes for uneven spacing

        Vec3 const sphereExtents = (maxExtents - minExtents) / 2.0;
        double const sphereRadius = sphereExtents.getNorm() / sqrt(3.0);
        // Vec3 scale(1.0, 1.0, 1.0); // scale not needed

        LgIndex_t numStartingTriangles = 0;
        std::vector<TriNodes> startingTriangles;

        if (zoneType == ZoneType_FETriangle)
        {
            numStartingTriangles = srcJ;
            startingTriangles.reserve(numStartingTriangles);
            NodeMap_pa NM = TecUtilDataNodeGetReadableRef(srcSphericalZoneNum);
            for (LgIndex_t elem = 0; elem < numStartingTriangles; elem++)
            {
                Index_t const n1 = TecUtilDataNodeGetByRef(NM, elem + 1, 1) - 1;
                Index_t const n2 = TecUtilDataNodeGetByRef(NM, elem + 1, 2) - 1;
                Index_t const n3 = TecUtilDataNodeGetByRef(NM, elem + 1, 3) - 1;
                CHECK(0 <= n1 && n1 < numSrcNodes);
                CHECK(0 <= n2 && n2 < numSrcNodes);
                CHECK(0 <= n3 && n3 < numSrcNodes);
                startingTriangles.push_back(TriNodes(n1, n2, n3));
            }
        }
        else
        {
            numStartingTriangles = 2 * (srcI - 1)*(srcJ - 1);
            startingTriangles.reserve(numStartingTriangles);
            for (LgIndex_t jj = 0; jj < srcJ - 1; jj++)
                for (LgIndex_t ii = 0; ii < srcI - 1; ii++)
                {
                    TriNodes tri1(jj*srcI + ii, jj*srcI + (ii + 1), (jj + 1)*srcI + ii);
                    TriNodes tri2((jj + 1)*srcI + ii, jj*srcI + (ii + 1), (jj + 1)*srcI + (ii + 1));
                    startingTriangles.push_back(tri1);
                    startingTriangles.push_back(tri2);
                }
        }
        CHECK(numStartingTriangles > 0);
        CHECK(startingTriangles.size() == size_t(numStartingTriangles));

        std::vector<Vec3> constraintNodes; // no constraints for now, TODO: get from zones 2+
        std::vector<Edge> constraintSegments;
        if (errMsg==NULL)
        {
            for (auto zoneIter = TecUtilSetGetNextMember(constraintZones, TECUTILSETNOTMEMBER);
                 zoneIter != TECUTILSETNOTMEMBER;
                 zoneIter = TecUtilSetGetNextMember(constraintZones, zoneIter))
            {
                EntIndex_t zoneNum = EntIndex_t(zoneIter); // zoneNum is 1-based because it comes from a TecUtil function

                LgIndex_t zoneType = TecUtilZoneGetType(zoneNum);
                CHECK(VALID_ENUM(zoneType, ZoneType_e));

                LgIndex_t srcI = 0;
                LgIndex_t srcJ = 0;
                LgIndex_t srcK = 0;
                TecUtilZoneGetIJK(zoneNum, &srcI, &srcJ, &srcK);
                CHECK(srcI > 0 && srcJ > 0 && srcK > 0);

                FieldData_pa XD = TecUtilDataValueGetReadableRef(zoneNum, xVarNum);
                FieldData_pa YD = TecUtilDataValueGetReadableRef(zoneNum, yVarNum);
                FieldData_pa ZD = TecUtilDataValueGetReadableRef(zoneNum, zVarNum);

                CHECK(zoneType == ZoneType_Ordered && srcI > 1 && srcJ == 1 && srcK == 1); // checked above in loop
                CHECK(TecUtilDataValueGetLocationByRef(XD) == ValueLocation_Nodal);
                CHECK(TecUtilDataValueGetLocationByRef(YD) == ValueLocation_Nodal);
                CHECK(TecUtilDataValueGetLocationByRef(ZD) == ValueLocation_Nodal);

                LgIndex_t const numSrcNodes = TecUtilDataValueGetCountByRef(XD);
                CHECK(numSrcNodes == TecUtilDataValueGetCountByRef(YD));

                Index_t const offset = Index_t(constraintNodes.size());
                for (LgIndex_t srcNode = 0; srcNode < numSrcNodes; srcNode++)
                {
                    double const xx = TecUtilDataValueGetByRef(XD, srcNode + 1);
                    double const yy = TecUtilDataValueGetByRef(YD, srcNode + 1);
                    double const zz = TecUtilDataValueGetByRef(ZD, srcNode + 1);

                    constraintNodes.push_back(Vec3(xx,yy,zz));
                    if (srcNode > 0)
                        constraintSegments.push_back(Edge(srcNode-1 + offset, srcNode + offset));
                }
            }
        }

        std::vector<Vec3> newNodes; //out
        std::vector<Index_t> nodeEdge; //out
        std::vector<TriNodes> newTriangles; // out

        if (updateSphericalTriangulation(
                startingNodes,
                startingTriangles,
                sphereCenter,
                sphereRadius,
                constraintNodes,
                constraintSegments,
                false, // includeConstraintSegments
                newNodes,
                nodeEdge,
                newTriangles,
                errMsg))
        {
            errMsg = NULL;

            char const* edgeVarName = "Edge";
            EntIndex_t edgeVarNum = TecUtilVarGetNumByName(edgeVarName);
            if (edgeVarNum == TECUTILSETNOTMEMBER)
            {
                VERIFY(TecUtilDataSetAddVar(edgeVarName, NULL));
                edgeVarNum = TecUtilVarGetNumByName(edgeVarName);
                Set varZoneSet(edgeVarNum);
                CHECK(edgeVarNum != TECUTILSETNOTMEMBER);
                TecUtilStateChanged(StateChange_VarsAdded, ArbParam_t(varZoneSet.getRef()));
            }

            CHECK(newNodes.size() >= 3); // absolute minimum
            CHECK(newTriangles.size() >= 1); // absolute minimum

            std::stringstream ss;
            ss << "Updated triangulation of zone " << srcSphericalZoneNum;

            std::string messsage(ss.str());

            numNewNodes = LgIndex_t(newNodes.size());
            numNewTriangles = LgIndex_t(newTriangles.size());

            VERIFY(TecUtilDataSetAddZone(messsage.c_str(),
                    numNewNodes,
                    numNewTriangles,
                    3,
                    ZoneType_FETriangle,
                    NULL));
            EntIndex_t const dstZoneNum = TecUtilDataSetGetNumZones();
            EntIndex_t const numVars = TecUtilDataSetGetNumVars();

            for (EntIndex_t varNum = 1; varNum <= numVars; varNum++)
            {
                FieldData_pa srcFD = TecUtilDataValueGetReadableRef(srcSphericalZoneNum, varNum);
                CHECK(TecUtilDataValueGetLocationByRef(srcFD) == ValueLocation_Nodal);

                FieldData_pa dstFD = TecUtilDataValueGetReadableRef(dstZoneNum, varNum);
                CHECK(TecUtilDataValueGetLocationByRef(srcFD) == ValueLocation_Nodal);

                for (LgIndex_t dstNode = 0; dstNode < numNewNodes; dstNode++)
                {
                    double val;
                    if (varNum == xVarNum)
                        val = newNodes[dstNode].x();
                    else if (varNum == yVarNum)
                        val = newNodes[dstNode].y();
                    else if (varNum == zVarNum)
                        val = newNodes[dstNode].z();
                    else if (varNum == edgeVarNum)
                        val = nodeEdge[dstNode];
                    else
                        val = double(dstNode); // nodes are rearranged (mainly deleted and re-added
                    TecUtilDataValueSetByRef(dstFD, dstNode + 1, val);
                }
            }

            NodeMap_pa NM = TecUtilDataNodeGetRef(dstZoneNum);
            for (LgIndex_t dstElem = 0; dstElem < numNewTriangles; dstElem++)
            {
                TecUtilDataNodeSetByRef(NM, dstElem + 1, 1, newTriangles[dstElem].v1() + 1);
                TecUtilDataNodeSetByRef(NM, dstElem + 1, 2, newTriangles[dstElem].v2() + 1);
                TecUtilDataNodeSetByRef(NM, dstElem + 1, 3, newTriangles[dstElem].v3() + 1);
            }

            Set dstZoneSet(dstZoneNum);
            TecUtilStateChanged(StateChange_ZonesAdded, ArbParam_t(dstZoneSet.getRef()));

            if ( !TecUtilSetAddMember(activeZones, dstZoneNum, FALSE) )
                errMsg = "Update spherical triangulation encountered memory error updating active zone set."; // hard to see how this would ever happen

            TecUtilSetRemoveMember(activeZones, srcSphericalZoneNum);
            TecUtilZoneSetActive(activeZones, AssignOp_Equals);
        }
    }

    TecUtilSetDealloc(&activeZones);
    TecUtilSetDealloc(&constraintZones);

    return errMsg;
}

/**
*/
static void STDCALL updateSphericalTriangulationMenuCallback(void)
{
    TecUtilLockStart(AddOnID);

    Index_t numNodes = 0;
    Index_t numTriangles = 0;
    Index_t numCoincidentNodesIgnored = 0;

    char const* errMsg = performUpdateSphericalTriangulation(numNodes, numTriangles, numCoincidentNodesIgnored);
    if (errMsg)
        TecUtilDialogMessageBox(errMsg, MessageBox_Error);
    else
    {
        if (TecUtilMacroIsRecordingActive())
            TecUtilMacroRecordAddOnCommand(ADDON_NAME, UPDATE_SPHERICAL_TRIANGULATION_COMMAND);

        std::stringstream ss;
        ss << "Spherical triangulation created. " << numTriangles << " triangles, " << numNodes << " nodes.";
        if (numCoincidentNodesIgnored > 0)
            ss << " (" << numCoincidentNodesIgnored << " coincident points ignored.)";
        TecUtilDialogMessageBox(ss.str().c_str(), MessageBox_Information);
    }

    TecUtilLockFinish(AddOnID);
}

/**
* This function is called when the
* $!EXTENDEDCOMMAND macro command is
* processed.
*/
static Boolean_t STDCALL macroCommandCallback(
    char *macroCommandString,  /* IN */
    char **errMsgPtr)             /* OUT (only if returning FALSE) */
{
    char const* errMsg = NULL;

    TecUtilLockStart(AddOnID);

    if (stricmp(macroCommandString, UPDATE_SPHERICAL_TRIANGULATION_COMMAND) == 0)
    {
        Index_t numNodes = 0;
        Index_t numTriangles = 0;
        Index_t numCoincidentNodesIgnored = 0;
        errMsg = performUpdateSphericalTriangulation(numNodes, numTriangles, numCoincidentNodesIgnored);
    }
    else
    {
        errMsg = "Unknown test triangulation option";
    }

    Boolean_t isOk = TRUE;
    if (errMsg != NULL)
    {
        int const messageSize = int(strlen(errMsg) + 1);
        *errMsgPtr = TecUtilStringAlloc(messageSize, "String for Error Message");
        _snprintf(*errMsgPtr, messageSize, "%s", "Error processing macro command");
        isOk = FALSE;
    }

    TecUtilLockFinish(AddOnID);

    return isOk;
}


/*
*/
namespace {
    EntIndex_t getVarByNames(
        char const* name1,
        char const* name2 = NULL)
    {
        REQUIRE(VALID_NON_ZERO_LEN_STR(name1));
        REQUIRE(name2 == NULL || VALID_NON_ZERO_LEN_STR(name2));

        EntIndex_t result = TecUtilVarGetNumByName(name1);
        if (result == TECUTILSETNOTMEMBER)
        {
            if (name2 == NULL)
                throw std::string("Cannot find variable named '") + name1 + "'";
            else
            {
                result = TecUtilVarGetNumByName(name2);
                if (result == TECUTILSETNOTMEMBER)
                    throw std::string("Cannot find variable named '") + name1 + "' or '" + name2 + "'";
            }
        }

        ENSURE(result >= 0 && result <= TecUtilDataSetGetNumVars());
        return result;
    }
}

/*
* Creates norm field of rotation from he system to ff system which is the
* smallest angle to rotate any of the vectors or their opposites into the other system.
* - All variables are one-based like the TecUtil layer
* - Returns error message on error or NULL if everything is okay
*   The error message is a constant string and should not be freed
*/
char const* calculateHeffNormField(
    EntIndex_t he1xVarNum, EntIndex_t he1yVarNum, EntIndex_t he1zVarNum,
    EntIndex_t he2xVarNum, EntIndex_t he2yVarNum, EntIndex_t he2zVarNum,
    EntIndex_t he3xVarNum, EntIndex_t he3yVarNum, EntIndex_t he3zVarNum,
    EntIndex_t ff1xVarNum, EntIndex_t ff1yVarNum, EntIndex_t ff1zVarNum,
    EntIndex_t ff2xVarNum, EntIndex_t ff2yVarNum, EntIndex_t ff2zVarNum,
    EntIndex_t ff3xVarNum, EntIndex_t ff3yVarNum, EntIndex_t ff3zVarNum)
{
    const char* errMsg = NULL;

    EntIndex_t const numVars = TecUtilDataSetGetNumVars();
    if (he1xVarNum < 1 || he1xVarNum > numVars)
        errMsg = "he1x variable does not exist";
    else if (he1yVarNum < 1 || he1yVarNum > numVars)
        errMsg = "he1y variable does not exist";
    else if (he1zVarNum < 1 || he1zVarNum > numVars)
        errMsg = "he1z variable does not exist";
    else if (he2xVarNum < 1 || he2xVarNum > numVars)
        errMsg = "he2x variable does not exist";
    else if (he2yVarNum < 1 || he2yVarNum > numVars)
        errMsg = "he2y variable does not exist";
    else if (he2zVarNum < 1 || he2zVarNum > numVars)
        errMsg = "he2z variable does not exist";
    else if (he3xVarNum < 1 || he3xVarNum > numVars)
        errMsg = "he3x variable does not exist";
    else if (he3yVarNum < 1 || he3yVarNum > numVars)
        errMsg = "he3y variable does not exist";
    else if (he3zVarNum < 1 || he3zVarNum > numVars)
        errMsg = "he3z variable does not exist";
    else if (ff1xVarNum < 1 || ff1xVarNum > numVars)
        errMsg = "ff1x variable does not exist";
    else if (ff1yVarNum < 1 || ff1yVarNum > numVars)
        errMsg = "ff1y variable does not exist";
    else if (ff1zVarNum < 1 || ff1zVarNum > numVars)
        errMsg = "ff1z variable does not exist";
    else if (ff2xVarNum < 1 || ff2xVarNum > numVars)
        errMsg = "ff2x variable does not exist";
    else if (ff2yVarNum < 1 || ff2yVarNum > numVars)
        errMsg = "ff2y variable does not exist";
    else if (ff2zVarNum < 1 || ff2zVarNum > numVars)
        errMsg = "ff2z variable does not exist";
    else if (ff3xVarNum < 1 || ff3xVarNum > numVars)
        errMsg = "ff3x variable does not exist";
    else if (ff3yVarNum < 1 || ff3yVarNum > numVars)
        errMsg = "ff3y variable does not exist";
    else if (ff3zVarNum < 1 || ff3zVarNum > numVars)
        errMsg = "ff3z variable does not exist";

    EntIndex_t normVarNum = 0;
    if (errMsg == NULL)
    {
        if (!TecUtilDataSetAddVar("HE-FF Norm", NULL))
            errMsg = "Cannot create new variable";
        else
            normVarNum = TecUtilDataSetGetNumVars();
    }

    if (errMsg == NULL)
    {
        // for now, only calculate the norm field for IJK zones
        EntIndex_t numSrcZones = TecUtilDataSetGetNumZones();
        for (EntIndex_t zoneNum = 0; zoneNum <= numSrcZones && errMsg == NULL; zoneNum++)
        {
            if (TecUtilZoneIsEnabled(zoneNum) && TecUtilZoneIsActive(zoneNum))
            {
                // get info about each enabled zone
                LgIndex_t zoneType = TecUtilZoneGetType(zoneNum);
                CHECK(VALID_ENUM(zoneType, ZoneType_e));

                LgIndex_t srcI = 0;
                LgIndex_t srcJ = 0;
                LgIndex_t srcK = 0;
                TecUtilZoneGetIJK(zoneNum, &srcI, &srcJ, &srcK);
                CHECK(srcI > 0 && srcJ > 0 && srcK > 0);
                if (zoneType == ZoneType_Ordered && srcI > 1 && srcJ > 1 && srcK > 1)
                {
                    FieldData_pa he1XD = TecUtilDataValueGetReadableRef(zoneNum, he1xVarNum);
                    FieldData_pa he1YD = TecUtilDataValueGetReadableRef(zoneNum, he1yVarNum);
                    FieldData_pa he1ZD = TecUtilDataValueGetReadableRef(zoneNum, he1zVarNum);
                    FieldData_pa he2XD = TecUtilDataValueGetReadableRef(zoneNum, he2xVarNum);
                    FieldData_pa he2YD = TecUtilDataValueGetReadableRef(zoneNum, he2yVarNum);
                    FieldData_pa he2ZD = TecUtilDataValueGetReadableRef(zoneNum, he2zVarNum);
                    FieldData_pa he3XD = TecUtilDataValueGetReadableRef(zoneNum, he3xVarNum);
                    FieldData_pa he3YD = TecUtilDataValueGetReadableRef(zoneNum, he3yVarNum);
                    FieldData_pa he3ZD = TecUtilDataValueGetReadableRef(zoneNum, he3zVarNum);
                    FieldData_pa ff1XD = TecUtilDataValueGetReadableRef(zoneNum, ff1xVarNum);
                    FieldData_pa ff1YD = TecUtilDataValueGetReadableRef(zoneNum, ff1yVarNum);
                    FieldData_pa ff1ZD = TecUtilDataValueGetReadableRef(zoneNum, ff1zVarNum);
                    FieldData_pa ff2XD = TecUtilDataValueGetReadableRef(zoneNum, ff2xVarNum);
                    FieldData_pa ff2YD = TecUtilDataValueGetReadableRef(zoneNum, ff2yVarNum);
                    FieldData_pa ff2ZD = TecUtilDataValueGetReadableRef(zoneNum, ff2zVarNum);
                    FieldData_pa ff3XD = TecUtilDataValueGetReadableRef(zoneNum, ff3xVarNum);
                    FieldData_pa ff3YD = TecUtilDataValueGetReadableRef(zoneNum, ff3yVarNum);
                    FieldData_pa ff3ZD = TecUtilDataValueGetReadableRef(zoneNum, ff3zVarNum);

                    if (he1XD == NULL)
                        errMsg = "Cannot get field data for he1x";
                    else if (TecUtilDataValueGetLocationByRef(he1XD) != ValueLocation_Nodal)
                        errMsg = "he1x variable not nodal";
                    else if (he1YD == NULL)
                        errMsg = "Cannot get field data for he1y";
                    else if (TecUtilDataValueGetLocationByRef(he1YD) != ValueLocation_Nodal)
                        errMsg = "he1y variable not nodal";
                    else if (he1ZD == NULL)
                        errMsg = "Cannot get field data for he1z";
                    else if (TecUtilDataValueGetLocationByRef(he1ZD) != ValueLocation_Nodal)
                        errMsg = "he1z variable not nodal";
                    else if (he2XD == NULL)
                        errMsg = "Cannot get field data for he2x";
                    else if (TecUtilDataValueGetLocationByRef(he2XD) != ValueLocation_Nodal)
                        errMsg = "he2x variable not nodal";
                    else if (he2YD == NULL)
                        errMsg = "Cannot get field data for he2y";
                    else if (TecUtilDataValueGetLocationByRef(he2YD) != ValueLocation_Nodal)
                        errMsg = "he2y variable not nodal";
                    else if (he2ZD == NULL)
                        errMsg = "Cannot get field data for he2z";
                    else if (TecUtilDataValueGetLocationByRef(he2ZD) != ValueLocation_Nodal)
                        errMsg = "he2z variable not nodal";
                    else if (he3XD == NULL)
                        errMsg = "Cannot get field data for he3x";
                    else if (TecUtilDataValueGetLocationByRef(he3XD) != ValueLocation_Nodal)
                        errMsg = "he3x variable not nodal";
                    else if (he3YD == NULL)
                        errMsg = "Cannot get field data for he3y";
                    else if (TecUtilDataValueGetLocationByRef(he3YD) != ValueLocation_Nodal)
                        errMsg = "he3y variable not nodal";
                    else if (he3ZD == NULL)
                        errMsg = "Cannot get field data for he3z";
                    else if (TecUtilDataValueGetLocationByRef(he3ZD) != ValueLocation_Nodal)
                        errMsg = "he3z variable not nodal";
                    else if (ff1XD == NULL)
                        errMsg = "Cannot get field data for ff1x";
                    else if (TecUtilDataValueGetLocationByRef(ff1XD) != ValueLocation_Nodal)
                        errMsg = "ff1x variable not nodal";
                    else if (ff1YD == NULL)
                        errMsg = "Cannot get field data for ff1y";
                    else if (TecUtilDataValueGetLocationByRef(ff1YD) != ValueLocation_Nodal)
                        errMsg = "ff1y variable not nodal";
                    else if (ff1ZD == NULL)
                        errMsg = "Cannot get field data for ff1z";
                    else if (TecUtilDataValueGetLocationByRef(ff1ZD) != ValueLocation_Nodal)
                        errMsg = "ff1z variable not nodal";
                    else if (ff2XD == NULL)
                        errMsg = "Cannot get field data for ff2x";
                    else if (TecUtilDataValueGetLocationByRef(ff2XD) != ValueLocation_Nodal)
                        errMsg = "ff2x variable not nodal";
                    else if (ff2YD == NULL)
                        errMsg = "Cannot get field data for ff2y";
                    else if (TecUtilDataValueGetLocationByRef(ff2YD) != ValueLocation_Nodal)
                        errMsg = "ff2y variable not nodal";
                    else if (ff2ZD == NULL)
                        errMsg = "Cannot get field data for ff2z";
                    else if (TecUtilDataValueGetLocationByRef(ff2ZD) != ValueLocation_Nodal)
                        errMsg = "ff2z variable not nodal";
                    else if (ff3XD == NULL)
                        errMsg = "Cannot get field data for ff3x";
                    else if (TecUtilDataValueGetLocationByRef(ff3XD) != ValueLocation_Nodal)
                        errMsg = "ff3x variable not nodal";
                    else if (ff3YD == NULL)
                        errMsg = "Cannot get field data for ff3y";
                    else if (TecUtilDataValueGetLocationByRef(ff3YD) != ValueLocation_Nodal)
                        errMsg = "ff3y variable not nodal";
                    else if (ff3ZD == NULL)
                        errMsg = "Cannot get field data for ff3z";
                    else if (TecUtilDataValueGetLocationByRef(ff3ZD) != ValueLocation_Nodal)
                        errMsg = "ff3z variable not nodal";

                    FieldData_pa normVarFD = NULL;
                    if (errMsg == NULL)
                    {
                        normVarFD = TecUtilDataValueGetWritableNativeRef(zoneNum, normVarNum);
                        if (normVarFD == NULL)
                            errMsg = "Cannot get writable field data for norm var";
                        else if (TecUtilDataValueGetLocationByRef(normVarFD) != ValueLocation_Nodal)
                            errMsg = "norm var variable not nodal";
                    }

                    if (errMsg == NULL)
                    {
                        LgIndex_t const numSrcNodes = TecUtilDataValueGetCountByRef(he1XD);
                        CHECK(numSrcNodes == TecUtilDataValueGetCountByRef(ff3ZD)); //spot check

                        for (LgIndex_t index = 1; index <= numSrcNodes; index++)
                        {
                            double const he1x = TecUtilDataValueGetByRef(he1XD, index);
                            double const he1y = TecUtilDataValueGetByRef(he1YD, index);
                            double const he1z = TecUtilDataValueGetByRef(he1ZD, index);
                            double const he2x = TecUtilDataValueGetByRef(he2XD, index);
                            double const he2y = TecUtilDataValueGetByRef(he2YD, index);
                            double const he2z = TecUtilDataValueGetByRef(he2ZD, index);
                            double const he3x = TecUtilDataValueGetByRef(he3XD, index);
                            double const he3y = TecUtilDataValueGetByRef(he3YD, index);
                            double const he3z = TecUtilDataValueGetByRef(he3ZD, index);
                            double const ff1x = TecUtilDataValueGetByRef(ff1XD, index);
                            double const ff1y = TecUtilDataValueGetByRef(ff1YD, index);
                            double const ff1z = TecUtilDataValueGetByRef(ff1ZD, index);
                            double const ff2x = TecUtilDataValueGetByRef(ff2XD, index);
                            double const ff2y = TecUtilDataValueGetByRef(ff2YD, index);
                            double const ff2z = TecUtilDataValueGetByRef(ff2ZD, index);
                            double const ff3x = TecUtilDataValueGetByRef(ff3XD, index);
                            double const ff3y = TecUtilDataValueGetByRef(ff3YD, index);
                            double const ff3z = TecUtilDataValueGetByRef(ff3ZD, index);

                            double const heffNorm = tpcsm::calculateHeffNorm(
                                he1x, he1y, he1z,
                                he2x, he2y, he2z,
                                he3x, he3y, he3z,
                                ff1x, ff1y, ff1z,
                                ff2x, ff2y, ff2z,
                                ff3x, ff3y, ff3z);

                            TecUtilDataValueSetByRef(normVarFD, index, heffNorm);
                        }
                    }
                }
            }
        }
        Set dstVarSet(normVarNum);
        TecUtilStateChanged(StateChange_VarsAdded, ArbParam_t(dstVarSet.getRef()));
    }

    return errMsg;
}

/**
*/
static void STDCALL calculateHeffNormFieldCallback(void)
{
    TecUtilLockStart(AddOnID);

    std::string errMsg;

    try
    {
        EntIndex_t const he1xVarNum = getVarByNames("Eigenvector 1x", "EgnVec11");
        EntIndex_t const he1yVarNum = getVarByNames("Eigenvector 1y", "EgnVec12");
        EntIndex_t const he1zVarNum = getVarByNames("Eigenvector 1z", "EgnVec13");
        EntIndex_t const he2xVarNum = getVarByNames("Eigenvector 2x", "EgnVec21");
        EntIndex_t const he2yVarNum = getVarByNames("Eigenvector 2y", "EgnVec22");
        EntIndex_t const he2zVarNum = getVarByNames("Eigenvector 2z", "EgnVec23");
        EntIndex_t const he3xVarNum = getVarByNames("Eigenvector 3x", "EgnVec31");
        EntIndex_t const he3yVarNum = getVarByNames("Eigenvector 3y", "EgnVec32");
        EntIndex_t const he3zVarNum = getVarByNames("Eigenvector 3z", "EgnVec33");

        EntIndex_t const ff1xVarNum = getVarByNames("t1");
        EntIndex_t const ff1yVarNum = getVarByNames("t2");
        EntIndex_t const ff1zVarNum = getVarByNames("t3");
        EntIndex_t const ff2xVarNum = getVarByNames("n1");
        EntIndex_t const ff2yVarNum = getVarByNames("n2");
        EntIndex_t const ff2zVarNum = getVarByNames("n3");
        EntIndex_t const ff3xVarNum = getVarByNames("b1");
        EntIndex_t const ff3yVarNum = getVarByNames("b2");
        EntIndex_t const ff3zVarNum = getVarByNames("b3");

        char const* result = calculateHeffNormField(
            he1xVarNum, he1yVarNum, he1zVarNum, he2xVarNum, he2yVarNum, he2zVarNum, he3xVarNum, he3yVarNum, he3zVarNum,
            ff1xVarNum, ff1yVarNum, ff1zVarNum, ff2xVarNum, ff2yVarNum, ff2zVarNum, ff3xVarNum, ff3yVarNum, ff3zVarNum);
        if (result != NULL)
            errMsg = result;
    }
    catch (std::string& err)
    {
        errMsg = err;
    }

    if (!errMsg.empty())
        TecUtilDialogMessageBox(errMsg.c_str(), MessageBox_Error);
    else
    {
        //if (TecUtilMacroIsRecordingActive())
        //    TecUtilMacroRecordAddOnCommand(ADDON_NAME, CALCULATE_HEFF_NORM_FIELD_COMMAND);

        std::stringstream ss;
        ss << "HE/FF norm field calculation complete.";
        TecUtilDialogMessageBox(ss.str().c_str(), MessageBox_Information);
    }

    TecUtilLockFinish(AddOnID);
}


/**
 * When Tecplot first loads an add-on, it makes a
 * call to initialize the add-on. This function
 * must be named InitTecAddOn, as shown below.
 */
EXPORTFROMADDON void STDCALL InitTecAddOn(void)
{
    /*
     * NOTE:  TecUtilLockOn MUST be used for InitTecAddOn instead
     *        of TecUtilLockStart because AddonID has yet to be
     *        established.  TecUtilLockOn is in effect an "anonymous"
     *        locking of Tecplot (old style).
     */

    TecUtilLockOn();

    /*
     * The function TecUtilAddOnRegister() is the
     * only function that is REQUIRED to be called from
     * the initialization function.
     *
     * The information you give Tecplot by calling
     * this function will show up in the Help/About Add-ons
     * dialog box.
     */

    /*
     * Note that if your add-on requires a specific version of Tecplot,
     * you would check for that here using TecUtilGetTecplotVersion()
     */

    AddOnID = TecUtilAddOnRegister(110,
                                   ADDON_NAME,
                                   "V" ADDON_VERSION "(" TecVersionId ") " ADDON_DATE,
                                   "Tecplot, Inc.");

    // register state change callback
    {
        ArgList_pa ArgList;
        ArgList = TecUtilArgListAlloc();
        TecUtilArgListAppendFunction(ArgList, SV_CALLBACKFUNCTION, (const void *)StateChangeCallback);
        TecUtilArgListAppendInt(ArgList,      SV_STATECHANGEMODE,        StateChangeMode_v100);
        TecUtilArgListAppendInt(ArgList,      SV_STATECHANGECALLBACKAPI, StateChangeCallbackAPI_ChangeOnly);
        TecUtilStateChangeAddCallbackX(ArgList);
        TecUtilArgListDealloc(&ArgList);
    }

    TecUtilMenuAddOption("Tools",
                         "Update Spherical Triangulation Test",
                         '\0',
                         updateSphericalTriangulationMenuCallback);

    TecUtilMenuAddOption("Tools",
        "Calculate HE/FF Norm Field",
        '\0',
        calculateHeffNormFieldCallback);


    TecUtilMacroAddCommandCallback(ADDON_NAME,
        macroCommandCallback);

    /*
     * See note on TecUtilLockOn at start of this function.
     */
    TecUtilLockOff();
}
