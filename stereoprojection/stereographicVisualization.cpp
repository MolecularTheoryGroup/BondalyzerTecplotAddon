#include "TECADDON.h"
#include "ADDGLBL.h"
#include "GUIDEFS.h"
#if !defined (MSWIN)
#include <unistd.h>
#endif

#include <cstdlib>
#include <cmath>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>

#include "stereographicVisualization.h"
#include "stringformat.h"
#include "ClassMacros.h"

#include "ArgList.h"
#include "StyleValue.h"
#include "Set.h"
#include "StringList.h"
using namespace tecplot::toolbox;

static char const* FRAME_2D_AUX_TAG = "GBA.SP.Frame2dTag";
static char const* TITLE_TEXT_TAG = "GBA.SP.Title";
static double FRAME_POS_EPSLION = 0.01; // in inches on the paper, so 0.01 is rather small
static double FRAME_CASCADE_OFFSET = 0.25;

Boolean_t installStereoVisProbeOverride(StereoVisOp_e stereoVisOp);

namespace {
std::string getFrame2DTag(
    EntIndex_t zoneNum,
    int        subplot) // to support multiple plots per zone
{
    char tag[80];
    snprintf(tag, sizeof(tag), "zone%ldsubplot%ld", long(zoneNum), long(subplot));
    return std::string(tag);
}
}

struct XYZ
{
    XYZ() {}
    XYZ(double xval, double yval, double zval)
    {
        x = xval;
        y = yval;
        z = zval;
    }
    double x;
    double y;
    double z;
    void operator+=(XYZ p)
    {
        x += p.x;
        y += p.y;
        z += p.z;
    }
    void operator-=(XYZ p)
    {
        x -= p.x;
        y -= p.y;
        z -= p.z;
    }
    void operator*=(double val)
    {
        x *= val;
        y *= val;
        z *= val;
    }
    void operator/=(double val)
    {
        double const factor = 1/val;
        x *= factor;
        y *= factor;
        z *= factor;
    }
    XYZ operator+(XYZ const& p) const
    {
        return XYZ(x + p.x, y + p.y, z + p.z);
    }
    XYZ operator-(XYZ const& p) const
    {
        return XYZ(x - p.x, y - p.y, z - p.z);
    }
    // if XYZ is actually a vector, these functions are useful
    double length() const
    {
        return sqrt(x*x + y*y + z*z);
    }
    void normalize()
    {
        *this /= length();
    }
    double dotProduct(XYZ v) const
    {
        return x*v.x + y*v.y + z*v.z;
    }
    XYZ crossProduct(XYZ v) const
    {
        return XYZ(y*v.z - z*v.y, z*v.x - x*v.z, x*v.y - y*v.x);
    }
};

// returns INVALID_UNIQUE_ID if no such frame
UniqueID_t find2DFrameID(
    EntIndex_t zoneNum,
    EntIndex_t subplot) // to support multiple plots per zone
{
    UniqueID_t frame2D = INVALID_UNIQUE_ID;

    std::string frame2dTag = getFrame2DTag(zoneNum, subplot);
    TecUtilFrameLightweightLoopStart();
    do
    {
        AuxData_pa frameAuxData = TecUtilAuxDataFrameGetRef();
        ArbParam_t auxValue = INVALID_UNIQUE_ID;
        AuxDataType_e auxType = AuxDataType_Invalid;
        Boolean_t auxRetain = Boolean_t(-1); // invalid
        if ( TecUtilAuxDataGetItemByName(frameAuxData, FRAME_2D_AUX_TAG, &auxValue, &auxType, &auxRetain) )
        {
            if ( auxType == AuxDataType_String && strcmp((const char*)auxValue,frame2dTag.c_str())==0 )
                frame2D = TecUtilFrameGetUniqueID();
        }
    }
    while (TecUtilFrameLightweightLoopNext());
    TecUtilFrameLightweightLoopEnd();

    return frame2D;
}

class FrameInfo
{
    UNCOPYABLE_CLASS(FrameInfo);
private:
    UniqueID_t m_frameId;
    UniqueID_t m_datasetId;
    double m_xPos;
    double m_yPos;
    double m_width;
    double m_height;
    void getInfoFromCurFrame()
    {
        m_datasetId = TecUtilDataSetGetUniqueID();
        TecUtilFrameGetPosAndSize(&m_xPos, &m_yPos, &m_width, &m_height);
    }
public:
    FrameInfo(UniqueID_t baseFrameId)
        : m_frameId(baseFrameId)
        , m_datasetId(INVALID_UNIQUE_ID)
        , m_xPos(0.0)
        , m_yPos(0.0)
        , m_width(1.0)
        , m_height(1.0)
    {
        if (baseFrameId == TecUtilFrameGetUniqueID())
        {
            getInfoFromCurFrame();
        }
        else
        {
            TecUtilFrameLightweightLoopStart();
            do
            {
                if (baseFrameId == TecUtilFrameGetUniqueID())
                {
                    getInfoFromCurFrame();
                }
            } while (TecUtilFrameLightweightLoopNext());
            TecUtilFrameLightweightLoopEnd();
        }
    }
    double left() const { return m_xPos; }
    double right() const { return m_xPos + m_width; }
    double top() const { return m_yPos; }
    double bottom() const { return m_yPos + m_height; }
    double width() const { return m_width; }
    double height() const { return m_height; }
    UniqueID_t datasetId() const { return m_datasetId; }
};


class DatasetInfo
{
    UNCOPYABLE_CLASS(DatasetInfo);
private:
    UniqueID_t m_datasetId;
    EntIndex_t m_numVars;
    EntIndex_t m_srcZoneNum;
    ZoneType_e m_zoneType;
    LgIndex_t  m_iMax;
    LgIndex_t  m_jMax;
    LgIndex_t  m_kMax;
    StringList m_varNameList;
    Boolean_t m_getStyleInfo;
    char* m_stylesheetFName;
    std::vector<ValueLocation_e> m_valueLocations;
    std::vector<FieldDataType_e> m_fieldDataTypes;
    std::vector<FieldData_pa> m_fieldDataRefs;
    std::vector<void*> m_fieldDataRawPtrs;
    NodeMap_t* m_nodeMapRawPtr;

    void getInfoFromCurDataset()
    {
        m_zoneType = TecUtilZoneGetType(m_srcZoneNum);
        TecUtilZoneGetIJK(m_srcZoneNum, &m_iMax, &m_jMax, &m_kMax);

        m_numVars = TecUtilDataSetGetNumVars();
        m_valueLocations.resize(m_numVars, ValueLocation_Invalid);
        m_fieldDataTypes.resize(m_numVars, FieldDataType_Invalid);
        m_fieldDataRefs.resize(m_numVars, NULL);
        m_fieldDataRawPtrs.resize(m_numVars, NULL);
        for ( EntIndex_t varNum = 1; varNum <= m_numVars; varNum++ )
        {
            char* varName = NULL;
            if ( !TecUtilVarGetName(varNum, &varName) )
                throw std::bad_alloc();
            m_varNameList.append(varName);
            TecUtilStringDealloc(&varName);
            m_valueLocations[varNum-1] = TecUtilDataValueGetLocation(m_srcZoneNum, varNum);
            m_fieldDataRefs[varNum-1] = TecUtilDataValueGetReadableRef(m_srcZoneNum, varNum);
            TecUtilDataValueGetReadableRawPtr(m_srcZoneNum, varNum, &m_fieldDataRawPtrs[varNum-1], &m_fieldDataTypes[varNum-1]);
        }
        if ( m_zoneType != ZoneType_Ordered ) // TODO: Support poly(?)
            TecUtilDataNodeGetReadableRawPtr(m_srcZoneNum, &m_nodeMapRawPtr);

        if ( m_getStyleInfo )
        {
            TecUtilFileGetTempName(&m_stylesheetFName);
            TecUtilWriteStylesheet(m_stylesheetFName, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE);
        }
    }
public:
    ~DatasetInfo()
    {
        if ( m_stylesheetFName != NULL )
        {
            CHECK(m_getStyleInfo);
            remove(m_stylesheetFName);
            TecUtilStringDealloc(&m_stylesheetFName);
        }
    }
    DatasetInfo(
        UniqueID_t baseFrameId,
        EntIndex_t srcZoneNum,
        Boolean_t  getStyleInfo)
        : m_datasetId(INVALID_UNIQUE_ID)
        , m_numVars(-1)
        , m_srcZoneNum(srcZoneNum)
        , m_zoneType(ZoneType_Invalid)
        , m_iMax(-1)
        , m_jMax(-1)
        , m_kMax(-1)
        , m_getStyleInfo(getStyleInfo)
        , m_stylesheetFName(NULL)
        , m_nodeMapRawPtr(NULL)
    {
        if (baseFrameId == TecUtilFrameGetUniqueID())
        {
            getInfoFromCurDataset();
        }
        else
        {
            TecUtilFrameLightweightLoopStart();
            do
            {
                if (baseFrameId == TecUtilFrameGetUniqueID())
                {
                    getInfoFromCurDataset();
                }
            } while (TecUtilFrameLightweightLoopNext());
            TecUtilFrameLightweightLoopEnd();
        }
    }
    UniqueID_t datasetId() const { return m_datasetId; }
    EntIndex_t numVars() const { return m_numVars; }
    StringList const &varNameList() const { return m_varNameList; }
    FieldDataType_e fieldDataType(EntIndex_t varNum) const { return m_fieldDataTypes[varNum-1]; }
    ValueLocation_e valueLocation(EntIndex_t varNum) const { return m_valueLocations[varNum-1]; }
    ZoneType_e zoneType() const { return m_zoneType; }
    LgIndex_t iMax() const { return m_iMax; }
    LgIndex_t jMax() const { return m_jMax; }
    LgIndex_t kMax() const { return m_kMax; }
    char const* stylesheetFileName() const { return m_stylesheetFName; }

    FieldData_pa fieldDataGetRef(EntIndex_t varNum) { return m_fieldDataRefs[varNum - 1]; }
    void* fieldDataGetRawPtr(EntIndex_t varNum) { return m_fieldDataRawPtrs[varNum - 1]; }
    NodeMap_t* nodeMapRawPtr() { return m_nodeMapRawPtr; }

    // convenient access to arrays for sending to TecUtilCreateZoneX
    FieldDataType_e* getFieldDataTypesArray() { return &m_fieldDataTypes[0]; }
    ValueLocation_e* getValueLocationsArray() { return &m_valueLocations[0]; }
};


void createProjectedZone(
    UniqueID_t srcFrameId,
    EntIndex_t srcZoneNum,
    EntIndex_t subplot,
    XYZ        projectionCenter,
    XYZ        projectionTop,
    Boolean_t  copyStyleInfo)
{
    REQUIRE(srcFrameId != INVALID_UNIQUE_ID);
    REQUIRE(srcZoneNum > 0);

    DatasetInfo srcDatasetInfo(srcFrameId, srcZoneNum, copyStyleInfo);

    TecUtilDataSetCreate(getFrame2DTag(srcZoneNum, subplot).c_str(), srcDatasetInfo.varNameList().getRef(), TRUE);
    
    //TecUtilDataSetAddZone(getFrame2DTag(srcZoneNum, offset).c_str(), 2, 2, 1, ZoneType_Ordered, &(srcDatasetInfo.fieldDataTypes()[0]));

    ZoneType_e const dstZoneType = srcDatasetInfo.zoneType();
    LgIndex_t const dstIMax = srcDatasetInfo.iMax();
    LgIndex_t const dstJMax = srcDatasetInfo.jMax();
    LgIndex_t const dstKMax = srcDatasetInfo.kMax();

    ArgList argList;
    argList.appendString(SV_NAME, getFrame2DTag(srcZoneNum, subplot));
    argList.appendInt(SV_ZONETYPE, dstZoneType); // todo: copy from zone
    argList.appendInt(SV_IMAX, dstIMax);
    argList.appendInt(SV_JMAX, dstJMax);
    argList.appendInt(SV_KMAX, dstKMax);
    argList.appendArray(SV_VARDATATYPE, srcDatasetInfo.getFieldDataTypesArray());
    argList.appendArray(SV_VALUELOCATION, srcDatasetInfo.getValueLocationsArray());
    argList.appendInt(SV_DEFERVARCREATION, TRUE);
    // argList.appendInt(SV_DEFERNODEMAPCREATION, TRUE); // todo:implement
    TecUtilDataSetAddZoneX(argList.getRef());

    EntIndex_t const dstZoneNum = TecUtilDataSetGetNumZones();

    for ( EntIndex_t varNum = 1; varNum <= srcDatasetInfo.numVars(); varNum++ )
    {
        LgIndex_t numVals = 0;
        if ( dstZoneType == ZoneType_Ordered )
        {
            if ( srcDatasetInfo.valueLocation(varNum) == ValueLocation_Nodal )
                numVals = dstIMax * dstJMax * dstKMax;
            else
                numVals = std::max(1,dstIMax-1) * std::max(1,dstJMax-1) * std::max(1,dstKMax-1);
        }
        else
        {
            if (srcDatasetInfo.valueLocation(varNum) == ValueLocation_Nodal)
                numVals = dstIMax;
            else
                numVals = dstJMax;
        }

        TecUtilDataValueAlloc(dstZoneNum, varNum);
        FieldData_pa dstFD = TecUtilDataValueGetWritableNativeRef(dstZoneNum, varNum);

        if ( srcDatasetInfo.fieldDataType(varNum) != FieldDataType_Bit )
        {
            void* srcRawPtr = srcDatasetInfo.fieldDataGetRawPtr(varNum);
            TecUtilDataValueArraySetByRef(dstFD, 1, numVals, srcRawPtr);
        }
        else
        {
            FieldData_pa srcFD = srcDatasetInfo.fieldDataGetRef(varNum);
            FieldValueGetFunction_pf getFunc = TecUtilDataValueRefGetGetFunc(srcFD);
            FieldValueSetFunction_pf setFunc = TecUtilDataValueRefGetSetFunc(dstFD);
            for ( EntIndex_t pos = 0; pos < numVals; pos++ )
            {
                double const val = getFunc(srcFD, pos);
                setFunc(dstFD, pos, val);
            }
        }
    }

    // copy connectivity
    if ( dstZoneType != ZoneType_Ordered ) // TODO:Support Poly (?)
    {
        TecUtilDataNodeAlloc(dstZoneNum);
        NodeMap_t* srcNodeRawPtr = srcDatasetInfo.nodeMapRawPtr();
        NodeMap_pa dstNM = TecUtilDataNodeGetWritableRef(dstZoneNum);
        LgIndex_t numValues = dstJMax*dstKMax;

        // TecUtilDataNodeArraySetByRef requires node+1 values, so we have to copy it
        NodeMap_t* copyOfNodeRawPtr = (NodeMap_t*)malloc(sizeof(NodeMap_t)*numValues);
        if ( copyOfNodeRawPtr )
        {
            for ( LgIndex_t count = 0; count < numValues; count++ )
            {
                NodeMap_t const node = srcNodeRawPtr[count]+1;
                copyOfNodeRawPtr[count] = node;
            }

            TecUtilDataNodeArraySetByRef(dstNM, 1, numValues, copyOfNodeRawPtr);
            free(copyOfNodeRawPtr);
        }
    }

    Set dstZoneSet(dstZoneNum);
    TecUtilStateChanged(StateChange_ZonesAdded, ArbParam_t(dstZoneSet.getRef()));

    if ( copyStyleInfo )
    {
        char const* stylesheetName = srcDatasetInfo.stylesheetFileName();
        CHECK(VALID_NON_ZERO_LEN_STR(stylesheetName));
        TecUtilReadStylesheet(stylesheetName, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE);
        TecUtilZoneSetActive(dstZoneSet.getRef(), AssignOp_Equals);
    }
    else
    {
        TecUtilFrameSetPlotType(PlotType_Cartesian2D);
    }

    TecUtilViewFit();
}


EntIndex_t getDataSetVarNum(char const* varName)
{
    REQUIRE(VALID_NON_ZERO_LEN_STR(varName));

    EntIndex_t varNum = TecUtilVarGetNumByName(varName);
    if ( varNum == TECUTILSETNOTMEMBER )
    {
        // var does not exists, add new one
        ArgList argList;
        argList.appendString(SV_NAME, varName);
        //argList.appendArray(SV_VARDATATYPE, (void *)varDataType);
        //argList.appendArray(SV_VALUELOCATION, (void *)valueLocation);
        argList.appendInt(SV_SHAREVARWITHALLZONES, FALSE);
        TecUtilDataSetAddVarX(argList.getRef());
        varNum = TecUtilDataSetGetNumVars();
    }

    ENSURE(1 <= varNum && varNum <= TecUtilDataSetGetNumVars());
    return varNum;
}

void set3DVarAssignments(
    EntIndex_t xVarNum,
    EntIndex_t yVarNum,
    EntIndex_t zVarNum)
{
    REQUIRE(xVarNum >= 1);
    REQUIRE(yVarNum >= 1);
    REQUIRE(zVarNum >= 1);

    StyleValue sv(SV_THREEDAXIS);
    VERIFY(StyleValue::returnCodeOk(sv.set(xVarNum, SV_XDETAIL, SV_VARNUM)));
    VERIFY(StyleValue::returnCodeOk(sv.set(yVarNum, SV_YDETAIL, SV_VARNUM)));
    VERIFY(StyleValue::returnCodeOk(sv.set(zVarNum, SV_ZDETAIL, SV_VARNUM)));
}

void set2DVarAssignments(
    EntIndex_t xVarNum,
    EntIndex_t yVarNum)
{
    REQUIRE(xVarNum >= 1);
    REQUIRE(yVarNum >= 1);

    StyleValue sv(SV_TWODAXIS);
    VERIFY(StyleValue::returnCodeOk(sv.set(xVarNum, SV_XDETAIL, SV_VARNUM)));
    VERIFY(StyleValue::returnCodeOk(sv.set(yVarNum, SV_YDETAIL, SV_VARNUM)));
    
    VERIFY(StyleValue::returnCodeOk(sv.set(-1.0, SV_XDETAIL, SV_RANGEMIN)));
    VERIFY(StyleValue::returnCodeOk(sv.set( 1.0, SV_XDETAIL, SV_RANGEMAX)));
    VERIFY(StyleValue::returnCodeOk(sv.set(-1.0, SV_YDETAIL, SV_RANGEMIN)));
    VERIFY(StyleValue::returnCodeOk(sv.set( 1.0, SV_YDETAIL, SV_RANGEMAX)));
}

void setBlanking(EntIndex_t blankVarNum)
{
    REQUIRE(blankVarNum >= 1);

    StyleValue sv(SV_BLANKING, SV_VALUE);
    VERIFY(StyleValue::returnCodeOk(sv.set(TRUE, SV_INCLUDE)));

    LgIndex_t const constraitGroup = 1;
    VERIFY(StyleValue::returnCodeOk(sv.set(TRUE,          constraitGroup, SV_CONSTRAINT, SV_INCLUDE)));
    VERIFY(StyleValue::returnCodeOk(sv.set(blankVarNum,   constraitGroup, SV_CONSTRAINT, SV_VARA)));
    VERIFY(StyleValue::returnCodeOk(sv.set(RelOp_EqualTo, constraitGroup, SV_CONSTRAINT, SV_RELOP)));
    VERIFY(StyleValue::returnCodeOk(sv.set(1.0,           constraitGroup, SV_CONSTRAINT, SV_VALUECUTOFF)));
}

void setContours(
    UniqueID_t srcFrameID,
    EntIndex_t srcZoneNum,
    UniqueID_t dstFrameID,
    EntIndex_t dstZoneNum)
{
    REQUIRE(srcFrameID != INVALID_UNIQUE_ID);
    REQUIRE(srcZoneNum > 0);
    REQUIRE(dstFrameID != INVALID_UNIQUE_ID);
    REQUIRE(dstZoneNum > 0);

    // Get style from src frame and zone
    #define INVALID_GROUP   EntIndex_t(-1)
    #define INVALID_VARNUM  EntIndex_t(-1)
    #define INVALID_ZONENUM EntIndex_t(-1)

    StyleValue sv;
    ArgList argList;
    bool showMeshLayer = false; // try to set all values to invalid, but for bool that isn't really possible
    MeshType_e meshType = MeshType_Invalid;
    ColorIndex_t meshColor = ColorIndex_t(-1);
    double meshLineThinkness = -1.0;
    LinePattern_e meshLinePattern = LinePattern_Invalid;
    double meshPatternLength = -1.0;

    bool showContourLayer = false;
    ContourType_e srcZoneContourType = ContourType_Invalid;
    ColorIndex_t contourColor = ColorIndex_t(-1);
    double contourLineThinkness = -1.0;
    LinePattern_e contourLinePattern = LinePattern_Invalid;
    double contourPatternLength = -1.0;

    EntIndex_t floodGroup;
    EntIndex_t lineGroup;

    struct 
    {
        EntIndex_t varNum;
        LgIndex_t numLevels;
        double *levels;
        bool showLegend;
    } contourInfo[MaxContourGroups];
    for (EntIndex_t group = 0; group < MaxContourGroups; group++)
    {
        contourInfo[group].varNum = INVALID_VARNUM;
        contourInfo[group].numLevels = -1;
        contourInfo[group].levels = NULL;
        contourInfo[group].showLegend = false;
    }

    TecUtilFrameLightweightLoopStart();
    do
    {
        if (TecUtilFrameGetUniqueID() == srcFrameID)
        {
            VERIFY(StyleValue::returnCodeOk(sv.get(&showMeshLayer, SV_FIELDLAYERS, SV_SHOWMESH, SV_INCLUDE)));
            VERIFY(StyleValue::returnCodeOk(sv.get(&showContourLayer, SV_FIELDLAYERS, SV_SHOWCONTOUR, SV_INCLUDE)));

            VERIFY(StyleValue::returnCodeOk(sv.get(&meshType, srcZoneNum, SV_FIELDMAP, SV_MESH, SV_MESHTYPE)));
            VERIFY(StyleValue::returnCodeOk(sv.get(&meshColor, srcZoneNum, SV_FIELDMAP, SV_MESH, SV_COLOR)));
            VERIFY(StyleValue::returnCodeOk(sv.get(&meshLineThinkness, srcZoneNum, SV_FIELDMAP, SV_MESH, SV_LINETHICKNESS)));
            VERIFY(StyleValue::returnCodeOk(sv.get(&meshLinePattern, srcZoneNum, SV_FIELDMAP, SV_MESH, SV_LINEPATTERN)));
            VERIFY(StyleValue::returnCodeOk(sv.get(&meshPatternLength, srcZoneNum, SV_FIELDMAP, SV_MESH, SV_PATTERNLENGTH)));

            VERIFY(StyleValue::returnCodeOk(sv.get(&srcZoneContourType, srcZoneNum, SV_FIELDMAP, SV_CONTOUR, SV_CONTOURTYPE)));
            VERIFY(StyleValue::returnCodeOk(sv.get(&contourColor, srcZoneNum, SV_FIELDMAP, SV_CONTOUR, SV_COLOR)));
            VERIFY(StyleValue::returnCodeOk(sv.get(&floodGroup, srcZoneNum, SV_FIELDMAP, SV_CONTOUR, SV_FLOODCOLORING)));
            VERIFY(StyleValue::returnCodeOk(sv.get(&lineGroup, srcZoneNum, SV_FIELDMAP, SV_CONTOUR, SV_LINECONTOURGROUP)));
            VERIFY(StyleValue::returnCodeOk(sv.get(&contourLineThinkness, srcZoneNum, SV_FIELDMAP, SV_CONTOUR, SV_LINETHICKNESS)));
            VERIFY(StyleValue::returnCodeOk(sv.get(&contourLinePattern, srcZoneNum, SV_FIELDMAP, SV_CONTOUR, SV_LINEPATTERN)));
            VERIFY(StyleValue::returnCodeOk(sv.get(&contourPatternLength, srcZoneNum, SV_FIELDMAP, SV_CONTOUR, SV_PATTERNLENGTH)));

            for ( EntIndex_t group = 0; group < MaxContourGroups; group++ )
            {
                VERIFY(StyleValue::returnCodeOk(sv.get(&contourInfo[group].varNum, group+1, SV_GLOBALCONTOUR, SV_VAR)));
                if ( contourInfo[group].varNum != 0 )
                {
                    VERIFY(TecUtilContourGetLevels(group+1, &contourInfo[group].numLevels, &contourInfo[group].levels));
                    CHECK(VALID_REF_OR_NULL(contourInfo[group].levels));
                    VERIFY(StyleValue::returnCodeOk(sv.get(&contourInfo[group].showLegend, group+1, SV_GLOBALCONTOUR, SV_LEGEND, SV_SHOW)));
                }
            }
        }
    } while (TecUtilFrameLightweightLoopNext());
    TecUtilFrameLightweightLoopEnd();

    // set layers in dst frame
    VERIFY(StyleValue::returnCodeOk(sv.set(showMeshLayer, SV_FIELDLAYERS, SV_SHOWMESH, SV_INCLUDE)));
    VERIFY(StyleValue::returnCodeOk(sv.set(showContourLayer, SV_FIELDLAYERS, SV_SHOWCONTOUR, SV_INCLUDE)));

    // set contours in dst frame
    CHECK(TecUtilFrameGetUniqueID() == dstFrameID);
    {
        for (EntIndex_t group = 0; group < MaxContourGroups; group++)
        {
            if ( contourInfo[group].varNum != 0 )
            {
                //var
                argList.clear();
                argList.appendInt(SV_CONTOURGROUP, group+1);
                argList.appendInt(SV_VAR, contourInfo[group].varNum);
                TecUtilContourSetVariableX(argList.getRef());

                // levels
                CHECK(TecUtilFrameGetUniqueID() == dstFrameID);
                argList.clear();
                argList.appendInt(SV_CONTOURLEVELACTION, LgIndex_t(ContourLevelAction_New));
                argList.appendInt(SV_CONTOURGROUP, group+1);
                argList.appendInt(SV_NUMVALUES, contourInfo[group].numLevels);
                argList.appendArray(SV_RAWDATA, (void*)contourInfo[group].levels);
                TecUtilContourLevelX(argList.getRef());
                // legend
                VERIFY(StyleValue::returnCodeOk(sv.set(contourInfo[group].showLegend, group+1, SV_GLOBALCONTOUR, SV_LEGEND, SV_SHOW)));
            }
        }
    }

    // set zone style in dst frame
    {
        Set dstZoneSet(dstZoneNum);
        VERIFY(StyleValue::returnCodeOk(sv.set(meshType, dstZoneSet, SV_FIELDMAP, SV_MESH, SV_MESHTYPE)));
        VERIFY(StyleValue::returnCodeOk(sv.set(meshColor, dstZoneSet, SV_FIELDMAP, SV_MESH, SV_COLOR)));
        VERIFY(StyleValue::returnCodeOk(sv.set(meshLineThinkness, dstZoneSet, SV_FIELDMAP, SV_MESH, SV_LINETHICKNESS)));
        VERIFY(StyleValue::returnCodeOk(sv.set(meshLinePattern, dstZoneSet, SV_FIELDMAP, SV_MESH, SV_LINEPATTERN)));
        VERIFY(StyleValue::returnCodeOk(sv.set(meshPatternLength, dstZoneSet, SV_FIELDMAP, SV_MESH, SV_PATTERNLENGTH)));

        VERIFY(StyleValue::returnCodeOk(sv.set(srcZoneContourType, dstZoneSet, SV_FIELDMAP, SV_CONTOUR, SV_CONTOURTYPE)));
        VERIFY(StyleValue::returnCodeOk(sv.set(contourColor, dstZoneSet, SV_FIELDMAP, SV_CONTOUR, SV_COLOR)));
        VERIFY(StyleValue::returnCodeOk(sv.set(contourLineThinkness, dstZoneSet, SV_FIELDMAP, SV_CONTOUR, SV_LINETHICKNESS)));
        VERIFY(StyleValue::returnCodeOk(sv.set(contourLinePattern, dstZoneSet, SV_FIELDMAP, SV_CONTOUR, SV_LINEPATTERN)));
        VERIFY(StyleValue::returnCodeOk(sv.set(contourPatternLength, dstZoneSet, SV_FIELDMAP, SV_CONTOUR, SV_PATTERNLENGTH)));
        VERIFY(StyleValue::returnCodeOk(sv.set(floodGroup, dstZoneSet, SV_FIELDMAP, SV_CONTOUR, SV_FLOODCOLORING)));
        VERIFY(StyleValue::returnCodeOk(sv.set(lineGroup, dstZoneSet, SV_FIELDMAP, SV_CONTOUR, SV_LINECONTOURGROUP)));
    }

    // clean up
    for (EntIndex_t group = 0; group < MaxContourGroups; group++)
        TecUtilArrayDealloc((void **)&contourInfo->levels);
}


void getTransformationCenter(
    EntIndex_t zoneNum,
    double& txorigin,
    double& tyorigin,
    double& tzorigin)
{
    REQUIRE(zoneNum >= 1);
    txorigin = 1.0; // TODO: get actual center based on atom location in aux data
    tyorigin = 1.0;
    tzorigin = 1.0;
}

UniqueID_t updateDependentFrame(
    UniqueID_t    srcFrameId,
    EntIndex_t    srcZoneNum,
    EntIndex_t    subplot, // to support multiple plots per zone
    UniqueID_t    dstFrameId, // if INVALID_UNIQUE_ID, create new frame
    EntIndex_t    dstZoneNum,
    Boolean_t     doProjection,
    XYZ           projectionCenter, // this point on the sphere will be the center of the projection (the "south pole")
    XYZ           projectionTop, // this point on the sphere will end up directly above the center in the projection
    Boolean_t     useCurrentDataset, // true = create new dataset for plotted zone, false = use base frame's dataset for plotted zone
    Boolean_t     copyCurrentStyle) // true = copy frame and style, false = make new frame with new style
{
    REQUIRE(IMPLICATION(srcFrameId == INVALID_UNIQUE_ID, dstFrameId != INVALID_UNIQUE_ID));
    REQUIRE(srcZoneNum > 0);
    REQUIRE(dstZoneNum > 0 || dstFrameId == INVALID_UNIQUE_ID);

    REQUIRE(VALID_BOOLEAN(doProjection));
    REQUIRE(VALID_BOOLEAN(useCurrentDataset));
    REQUIRE(VALID_BOOLEAN(copyCurrentStyle));
    REQUIRE(IMPLICATION(copyCurrentStyle, useCurrentDataset)); // can't copy frame without copying dataset too

    Boolean_t createNewDependantFrame = (dstFrameId == INVALID_UNIQUE_ID);
    if ( createNewDependantFrame ) // create dependent frame
    {
        FrameInfo baseFrameInfo(srcFrameId);

        double newFrameWidth = std::min(baseFrameInfo.width(), baseFrameInfo.height()) / 2.0;
        double newFrameHeight = newFrameWidth;
        double newFrameX = 0.0; // set in loop below
        double newFrameY = 0.0;

        // try frame positions until we find a free one
        Boolean_t continueLoop = TRUE;
        for ( EntIndex_t curPositiontoVerify = 0; continueLoop; curPositiontoVerify++ )
        {
            // set the position
            switch (curPositiontoVerify)
            {
                case 0:  newFrameX = baseFrameInfo.left() - newFrameWidth;  newFrameY = baseFrameInfo.top();                     break; // side upper-left
                case 1:  newFrameX = baseFrameInfo.left() - newFrameWidth;  newFrameY = baseFrameInfo.bottom() - newFrameHeight; break; // side lower-left
                case 2:  newFrameX = baseFrameInfo.right();                 newFrameY = baseFrameInfo.top();                     break; // side upper-right
                case 3:  newFrameX = baseFrameInfo.right();                 newFrameY = baseFrameInfo.bottom() - newFrameHeight; break; // side lower-right
                case 4:  newFrameX = baseFrameInfo.left() - newFrameWidth;  newFrameY = baseFrameInfo.bottom();                  break; // bottom corner left
                case 5:  newFrameX = baseFrameInfo.left();                  newFrameY = baseFrameInfo.bottom();                  break; // bottom left
                case 6:  newFrameX = baseFrameInfo.right() - newFrameWidth; newFrameY = baseFrameInfo.bottom();                  break; // bottom right
                case 7:  newFrameX = baseFrameInfo.right();                 newFrameY = baseFrameInfo.bottom();                  break; // bottom corner right
                case 8:  newFrameX = baseFrameInfo.left() - newFrameWidth;  newFrameY = baseFrameInfo.top() - newFrameHeight;    break; // top corner left
                case 9:  newFrameX = baseFrameInfo.left();                  newFrameY = baseFrameInfo.top() - newFrameHeight;    break; // top left
                case 10: newFrameX = baseFrameInfo.right() - newFrameWidth; newFrameY = baseFrameInfo.top() - newFrameHeight;    break; // top right
                case 11: newFrameX = baseFrameInfo.right();                 newFrameY = baseFrameInfo.top() - newFrameHeight;    break; // top corner right
                case 12: newFrameX = baseFrameInfo.left()+FRAME_CASCADE_OFFSET; newFrameY = baseFrameInfo.top()+FRAME_CASCADE_OFFSET;break; // overlapping top corner (offset a bit)
                default: newFrameX += FRAME_CASCADE_OFFSET; newFrameY += FRAME_CASCADE_OFFSET; // cascade offset from last checked position
            }

            // check to see if frame already there
            continueLoop = FALSE;
            TecUtilFrameLightweightLoopStart();
            do
            {
                double curFrameX;
                double curFrameY;
                double curFrameWidth;
                double curFrameHeight;
                TecUtilFrameGetPosAndSize(&curFrameX, &curFrameY, &curFrameWidth, &curFrameHeight);
                if ( std::abs(curFrameX-newFrameX) < FRAME_POS_EPSLION && std::abs(curFrameY-newFrameY) < FRAME_POS_EPSLION )
                    continueLoop = TRUE; // reverify, but we are in a lightweight frame loop, so just continue
            } while (TecUtilFrameLightweightLoopNext());
            TecUtilFrameLightweightLoopEnd();
        }

        TecUtilDrawGraphics(FALSE);

        dstFrameId = INVALID_UNIQUE_ID;

        if ( copyCurrentStyle )
        {
            TecUtilMouseSetMode(MouseButtonMode_Select);
            TecUtilPickDeselectAll();
            TecUtilPickAddFrameByUniqueID(FALSE/*collecting objects*/, srcFrameId);
            TecUtilPickCopy();
            TecUtilPickPaste();
            TecUtilFrameSetPosAndSize(newFrameX, newFrameY, newFrameWidth, newFrameHeight);
            dstFrameId = TecUtilFrameGetUniqueID();
            dstZoneNum = srcZoneNum;
        }
        else
        {
            TecUtilFrameCreateNew(TRUE/*useSuppliedFrameSize*/, newFrameX, newFrameY, newFrameWidth, newFrameHeight);
            dstFrameId = TecUtilFrameGetUniqueID();
            if ( useCurrentDataset )
            {
                TecUtilFrameSetDataSet(baseFrameInfo.datasetId(), dstFrameId);
                TecUtilFrameSetPlotType(PlotType_Cartesian3D);
                dstZoneNum = srcZoneNum;
            }
            else
            {
                createProjectedZone(srcFrameId, srcZoneNum, subplot, projectionCenter, projectionTop, copyCurrentStyle);
                dstZoneNum = TecUtilDataSetGetNumZones();
                //copyContourStyle(baseFrameId, newFrameId);
            }
        }

        // set up aux data
        AuxData_pa frameAuxData = TecUtilAuxDataFrameGetRef();
        TecUtilAuxDataSetStrItem(frameAuxData, FRAME_2D_AUX_TAG, getFrame2DTag(srcZoneNum, subplot).c_str(), TRUE/*retain*/);

        // set up plot
        Set activeZoneSet(dstZoneNum);
        TecUtilZoneSetActive(activeZoneSet.getRef(), AssignOp_Equals);

        if ( !copyCurrentStyle )
        {
            TecUtilFieldLayerSetIsActive(SV_MESH, TRUE);
            TecUtilFieldLayerSetIsActive(SV_SHOWCONTOUR, TRUE);
            TecUtilFieldLayerSetIsActive(SV_USETRANSLUCENCY, FALSE);
        }

        TecUtilViewFit();
    }

    // create/update transformed variables
    if ( doProjection )
    {
        EntIndex_t const txVarNum = getDataSetVarNum("TransformedX");
        EntIndex_t const tyVarNum = getDataSetVarNum("TransformedY");
        EntIndex_t const tzVarNum = getDataSetVarNum("TransformedZ");
        EntIndex_t const tbVarNum = getDataSetVarNum("BlankVar");

        FieldData_pa XD = TecUtilDataValueGetReadableNativeRef(dstZoneNum, 1); // TODO: Change to the coordinates in use
        FieldData_pa YD = TecUtilDataValueGetReadableNativeRef(dstZoneNum, 2);
        FieldData_pa ZD = TecUtilDataValueGetReadableNativeRef(dstZoneNum, 3);

        FieldData_pa tXD = TecUtilDataValueGetWritableNativeRef(dstZoneNum, txVarNum);
        FieldData_pa tYD = TecUtilDataValueGetWritableNativeRef(dstZoneNum, tyVarNum);
        FieldData_pa tZD = TecUtilDataValueGetWritableNativeRef(dstZoneNum, tzVarNum);
        FieldData_pa tBlankD = TecUtilDataValueGetWritableNativeRef(dstZoneNum, tbVarNum);

        LgIndex_t iMax = 0;
        LgIndex_t jMax = 0;
        LgIndex_t kMax = 0;
        TecUtilZoneGetIJK(dstZoneNum, &iMax, &jMax, &kMax);
        LgIndex_t const numValues = (TecUtilZoneGetType(dstZoneNum) == ZoneType_Ordered) ? iMax*jMax*kMax : iMax;

        // for now calculate the center of the sphere until we can get it from auxdata TODO: get from auxdata
        //getTransformationCenter(zoneNum, txorigin, tyorigin, tzorigin);
        XYZ origin(0.0, 0.0, 0.0);
        for (LgIndex_t index = 1; index <= numValues; index++)
        {
            XYZ p(TecUtilDataValueGetByRef(XD, index),
                  TecUtilDataValueGetByRef(YD, index),
                  TecUtilDataValueGetByRef(ZD, index));
            origin += p;
        }
        origin /= double(numValues);

        // calculate sphere radius TODO: get actual radius from auxdata
        double averageRadius = 0.0;
        for (LgIndex_t index = 1; index <= numValues; index++)
        {
            XYZ p(TecUtilDataValueGetByRef(XD, index),
                  TecUtilDataValueGetByRef(YD, index),
                  TecUtilDataValueGetByRef(ZD, index));
            p -= origin;
            averageRadius += p.length();
        }
        averageRadius /= double(numValues);
        double const scaleFactor = averageRadius > 1.0e-30 ? 1.0/averageRadius : 1.0; // TODO:better handling of small values

        // calculate normalized vector from center of sphere out to the point to be the projection center
        XYZ n = projectionCenter-origin;
        n.normalize();
        // need to rotate so that normalized projection center is at 0,0,-1
        XYZ const targetVector(0, 0, -1);
        XYZ u = n.crossProduct(targetVector); // get vector perpendicular to both n and 0,0,-1
        u.normalize();

        double const cosTheta = n.dotProduct(targetVector);
        double const oneMinusCos2Theta = 1-cosTheta*cosTheta;
        double const sinTheta = oneMinusCos2Theta > 0.0 ? sqrt(oneMinusCos2Theta) : 0.0; // angle must be < 180 so sinTheta is always positive but allow for round off issues

        // now transform the up direction
        double cosTwist;
        double sinTwist;
        {
            XYZ p = projectionTop-origin;
            p.normalize();
            double const rx = (cosTheta + u.x*u.x*(1 - cosTheta))*p.x + (u.x*u.y*(1 - cosTheta) - u.z*sinTheta)*p.y + (u.x*u.z*(1 - cosTheta) + u.y*sinTheta)*p.z;
            double const ry = (u.y*u.x*(1 - cosTheta) + u.z*sinTheta)*p.x + (cosTheta + u.y*u.y*(1 - cosTheta))*p.y + (u.y*u.z*(1 - cosTheta) - u.x*sinTheta)*p.z;
            double const rz = (u.z*u.x*(1 - cosTheta) - u.y*sinTheta)*p.x + (u.z*u.y*(1 - cosTheta) + u.x*sinTheta)*p.y + (cosTheta + u.z*u.z*(1 - cosTheta))*p.z;

            XYZ flat(rx,ry,0);
            flat.normalize();
            cosTwist = flat.dotProduct(XYZ(0,1,0));
            double const oneMinusCos2Twist = 1 - cosTwist*cosTwist;
            sinTwist = oneMinusCos2Twist > 0.0 ? sqrt(oneMinusCos2Twist) : 0.0; // angle must be < 180 so sinTwist is always positive but allow for round off issues

            // check twist by using it on rx and we should get p again
            double const finalx = cosTwist*rx - sinTwist*ry; // ry*rx + sqrt(1-ry^2)*ry
            double const finaly = sinTwist*rx + cosTwist*ry; // -sqrt(1-ry^2)*rx + ry*ry
            double const finalz = rz;

            XYZ final(-finalx, finaly, finalz);
            CHECK(1-1e-6 <= final.length() && final.length() <= 1+1e-6);

            // some kind of check here
        }

        for ( LgIndex_t index = 1; index <= numValues; index++ )
        {
            // get point on unit sphere
            XYZ p(TecUtilDataValueGetByRef(XD, index),
                  TecUtilDataValueGetByRef(YD, index),
                  TecUtilDataValueGetByRef(ZD, index));
            p -= origin;
            p *= scaleFactor;

#if 0
            // rotate so that the direction of inquiry points toward the south pole (the place with no distortion)
            // In other words, rotate so that (nx,ny,nz) maps to (0,0,1) nz having been negated above
            //TODO: Include a "twist" rotation so one can specify the upward direction of the 2D plot
            double rx = n.z*p.x - n.x*p.z - n.y/(1+n.z)*(n.x*p.y - n.y*p.x);
            double ry = n.z*p.y - n.y*p.z + n.x/(1+n.z)*(n.x*p.y - n.y*p.x);
            double rz = n.z*p.z + n.x*p.x + n.y*p.y;
#endif

            double rx = (cosTheta + u.x*u.x*(1-cosTheta))*p.x   + (u.x*u.y*(1-cosTheta)-u.z*sinTheta)*p.y + (u.x*u.z*(1-cosTheta)+u.y*sinTheta)*p.z;
            double ry = (u.y*u.x*(1-cosTheta)+u.z*sinTheta)*p.x + (cosTheta + u.y*u.y*(1-cosTheta))*p.y   + (u.y*u.z*(1-cosTheta)-u.x*sinTheta)*p.z;
            double rz = (u.z*u.x*(1-cosTheta)-u.y*sinTheta)*p.x + (u.z*u.y*(1-cosTheta)+u.x*sinTheta)*p.y + (cosTheta + u.z*u.z*(1-cosTheta))*p.z;

            double blank = 0.0; // FALSE
            static double const epsilon = 5e-3;
            if ( rz > 1.0-epsilon ) // TODO: better processing of singularity (blanking)
            {
                rz = 1.0-epsilon;
                blank = 1.0; // TRUE
            }

            double const tx = ( cosTwist*rx - sinTwist*ry) / (1.0-rz);
            double const ty = ( sinTwist*rx + cosTwist*ry) / (1.0-rz);
            double const tz = rz;

            TecUtilDataValueSetByRef(tXD, index, -tx);
            TecUtilDataValueSetByRef(tYD, index, ty);
            TecUtilDataValueSetByRef(tZD, index, tz);
            TecUtilDataValueSetByRef(tBlankD, index, blank); // TODO: calculate actual blanking value and turn blanking on
        }

        // inform Tecplot of new vars
        {
            Set zoneSet(dstZoneNum);
            Set varSet;
            varSet += txVarNum;
            varSet += tyVarNum;
            varSet += tzVarNum;
            varSet += tbVarNum;

            ArgList argList;
            argList.appendInt(SV_STATECHANGE, StateChange_VarsAltered);
            argList.appendSet(SV_ZONELIST, zoneSet);
            argList.appendSet(SV_VARLIST, varSet);
            TecUtilStateChangedX(argList.getRef());
        }

        if (srcFrameId != INVALID_UNIQUE_ID)
        {
            // assign vars
            set3DVarAssignments(txVarNum, tyVarNum, tzVarNum);
            TecUtilFrameSetPlotType(PlotType_Cartesian2D);
            set2DVarAssignments(txVarNum, tyVarNum);

            // assign style
            setBlanking(tbVarNum);
            setContours(srcFrameId, srcZoneNum, dstFrameId, dstZoneNum);
        }
    }

    TecUtilFrameFitAllToPaper();
    TecUtilWorkViewFitAllFrames();

    TecUtilDrawGraphics(TRUE);
    TecUtilRedrawAll(TRUE);

    return dstFrameId;
}

void createVectorFromCSV(std::vector<std::string>& valueVector, char const* const csv)
{
    char const delimiterChar = ',';
    valueVector.empty();
    for ( char const* csvPos = csv; *csvPos != '\0'; csvPos++ )
    {
        std::string curString;
        for (/*nothing*/; *csvPos != delimiterChar && *csvPos != '\0'; csvPos++)
            curString += *csvPos;
        if (!curString.empty())
            valueVector.push_back(curString);
    }
}


void STDCALL stereoVisProbeCB(
    Boolean_t  wasSuccessful,
    Boolean_t  isNearestPoint,
    ArbParam_t clientData)
{
    REQUIRE(VALID_BOOLEAN(wasSuccessful));
    REQUIRE(VALID_BOOLEAN(isNearestPoint));
    StereoVisOp_e stereoVisOp = StereoVisOp_e(clientData);
    REQUIRE(VALID_ENUM(stereoVisOp, StereoVisOp_e));

    TecUtilLockStart(AddOnID);

    if ( wasSuccessful && ( TecUtilFrameGetPlotType() == PlotType_Cartesian3D || TecUtilFrameGetPlotType() == PlotType_Cartesian2D ) )
    {
        EntIndex_t const zoneNum = TecUtilProbeFieldGetZone();

        if (zoneNum == TECUTILBADZONENUMBER)
        {
            TecUtilDialogErrMsg("Failed to get zone number. Please click on a zone.");
        }
        else
        {
            // get axis vars
            EntIndex_t xVarNum = 1;
            EntIndex_t yVarNum = 2;
            EntIndex_t zVarNum = 3;
            if (TecUtilFrameGetPlotType() == PlotType_Cartesian3D)
            {
                // get actual x, y z,
                TecUtilAxisGetVarAssignments(&xVarNum, &yVarNum, &zVarNum);
            }
            else
            {
                // use the original x,y,z TODO: get from AuxData
            }

            if (stereoVisOp == StereoVisOp_SinglePt)
            {
                EntIndex_t const subplot = 0; // TODO if probing 2D, get from frame

                // get probe location
                double const probeX = TecUtilProbeFieldGetValue(xVarNum);
                double const probeY = TecUtilProbeFieldGetValue(yVarNum);
                double const probeZ = TecUtilProbeFieldGetValue(zVarNum);

                UniqueID_t frame3d;
                UniqueID_t frame2d;
                if ( TecUtilFrameGetPlotType() == PlotType_Cartesian2D )
                {
                    frame3d = INVALID_UNIQUE_ID;
                    frame2d = TecUtilFrameGetUniqueID();
                    // TODO: Check aux data to make sure this is a good frame
                }
                else
                {
                    frame3d = TecUtilFrameGetUniqueID();
                    frame2d = find2DFrameID(zoneNum, subplot);
                }

                if (frame2d == INVALID_UNIQUE_ID)
                {
                    frame2d = updateDependentFrame(
                        frame3d, zoneNum, subplot,
                        INVALID_UNIQUE_ID/*dstFrame:createsNewFrame*/, 0/*dstZone:zoneNumInvalid*/,
                        TRUE, // doProjection
                        XYZ(probeX, probeY, probeZ), // projectionCenter
                        XYZ(probeX, probeY, probeZ+1), // projectionTop
                        FALSE/*useCurrentDataSet*/,
                        FALSE/*copyCurrentStyle*/);
                    if (frame2d == INVALID_UNIQUE_ID)
                        TecUtilDialogErrMsg("Cannot find 2-D frame, and cannot create new frame.");
                    else
                    {
                        char buff[50];
                        snprintf(buff, sizeof(buff), "Stereographic centered at %.2lg,%.2lg,%.2lg", double(probeX),double(probeY),double(probeZ));
                        Text_ID textID = TecUtilTextCreate(CoordSys_Frame, 50.0, 90.0, Units_Frame, 3.6, buff); // same size as default axis titles
                        TecUtilTextSetAnchor(textID, TextAnchor_Center);
                        TecUtilTextSetMacroFunctionCmd(textID, TITLE_TEXT_TAG);
                    }
                }
                else
                {
                    TecUtilFramePopByUniqueID(frame2d);
                    EntIndex_t dstZoneNum = 1; // only one zone in each frame right now
                    UniqueID_t result = updateDependentFrame(
                        frame3d, zoneNum, subplot,
                        frame2d, dstZoneNum,
                        TRUE, // doProjection
                        XYZ(probeX, probeY, probeZ), // projectionCenter
                        XYZ(probeX, probeY, probeZ+1), // projectionTop
                        FALSE/*useCurrentDataSet*/,
                        FALSE/*copyCurrentStyle*/);
                    CHECK(result == frame2d);
                    for ( Text_ID textID = TecUtilTextGetBase(); textID != TECUTILBADID; textID = TecUtilTextGetNext(textID) )
                    {
                        char* macroFunctionCmd = NULL;
                        if ( TecUtilTextGetMacroFunctionCmd(textID, &macroFunctionCmd) )
                        {
                            if ( strcmp(macroFunctionCmd, TITLE_TEXT_TAG) == 0 )
                            {
                                char buff[50];
                                snprintf(buff, sizeof(buff), "Stereographic centered at %.2lg,%.2lg,%.2lg", double(probeX), double(probeY), double(probeZ));
                                TecUtilTextSetString(textID, buff);
                                TecUtilStringDealloc(&macroFunctionCmd);
                                break; // leave loop
                            }
                            TecUtilStringDealloc(&macroFunctionCmd);
                        }
                    }
                }
                TecUtilFramePopByUniqueID(frame3d);
            }
            else if ( stereoVisOp == StereoVisOp_AllCps )
            {
                //TecUtilDialogMessageBox("StereoVisOp_AllCps", MessageBoxType_Information);

                LgIndex_t iMax = 0;
                LgIndex_t jMax = 0;
                LgIndex_t kMax = 0;
                TecUtilZoneGetIJK(zoneNum, &iMax, &jMax, &kMax);
                LgIndex_t const numNodes = ( TecUtilZoneGetType(zoneNum) == ZoneType_Ordered ? iMax*jMax*kMax : iMax );

                AuxData_pa auxDataRef = TecUtilAuxDataZoneGetRef(zoneNum);
                if ( auxDataRef )
                {
                    ArbParam_t    auxDataValue;
                    AuxDataType_e auxDataType;
                    Boolean_t     auxDataRetainFlag;

                    std::string sphereName;
                    if ( TecUtilAuxDataGetItemByName(auxDataRef, "GBA.SphereCP", &auxDataValue, &auxDataType, &auxDataRetainFlag) )
                    {
                        char* name = (char*)auxDataValue;
                        sphereName = name;
                        TecUtilStringDealloc(&name);
                    }

                    std::vector<std::string> numbers;
                    if (TecUtilAuxDataGetItemByName(auxDataRef, "GBA.SphereConstrainedNodeNums", &auxDataValue, &auxDataType, &auxDataRetainFlag))
                    {
                        char* numberCSV = (char*)auxDataValue;
                        createVectorFromCSV(numbers, numberCSV);
                        TecUtilStringDealloc(&numberCSV);
                    }

                    std::vector<std::string> names;
                    if (TecUtilAuxDataGetItemByName(auxDataRef, "GBA.SphereConstrainedNodeIntersectCPNames", &auxDataValue, &auxDataType, &auxDataRetainFlag))
                    {
                        char* nameCSV = (char*)auxDataValue;
                        createVectorFromCSV(names, nameCSV);
                        TecUtilStringDealloc(&nameCSV);
                    }

                    if ( numbers.size() > 0 && names.size() == numbers.size() )
                    {
                        for ( EntIndex_t subplot = 1; subplot <= EntIndex_t(numbers.size()); subplot++ ) // subplot 0 is the probe-at-position projection
                        {
                            char* endOfString = NULL;
                            LgIndex_t const nodeNum = strtol(numbers[subplot-1].c_str(), &endOfString, 10);
                            CHECK(endOfString==NULL || *endOfString=='\0');

                            if ( nodeNum > 0 && nodeNum <= numNodes )
                            {
                                // get probe location
                                double const probeX = TecUtilDataValueGetByZoneVar(zoneNum, xVarNum, nodeNum);
                                double const probeY = TecUtilDataValueGetByZoneVar(zoneNum, yVarNum, nodeNum);
                                double const probeZ = TecUtilDataValueGetByZoneVar(zoneNum, zVarNum, nodeNum);

                                UniqueID_t frame3d;
                                UniqueID_t frame2d;
                                if (TecUtilFrameGetPlotType() == PlotType_Cartesian2D)
                                {
                                    frame3d = INVALID_UNIQUE_ID;
                                    frame2d = TecUtilFrameGetUniqueID();
                                    // TODO: Check aux data to make sure this is a good frame
                                }
                                else
                                {
                                    frame3d = TecUtilFrameGetUniqueID();
                                    frame2d = find2DFrameID(zoneNum, subplot);
                                }

                                if (frame2d == INVALID_UNIQUE_ID)
                                {
                                    frame2d = updateDependentFrame(
                                        frame3d, zoneNum, subplot,
                                        INVALID_UNIQUE_ID/*dstFrame:createsNewFrame*/, 0/*dstZone:zoneNumInvalid*/,
                                        TRUE, // doProjection
                                        XYZ(probeX, probeY, probeZ), // projectionCenter
                                        XYZ(probeX, probeY, probeZ+1), // projectionTop
                                        FALSE/*useCurrentDataSet*/,
                                        FALSE/*copyCurrentStyle*/);
                                    if (frame2d == INVALID_UNIQUE_ID)
                                        TecUtilDialogErrMsg("Cannot find 2-D frame, and cannot create new frame.");
                                    else
                                    {
                                        char buff[200];
                                        snprintf(buff, sizeof(buff), "Stereographic %s on %s (node %ld)", sphereName.c_str(), names[subplot-1].c_str(), long(nodeNum));
                                        Text_ID textID = TecUtilTextCreate(CoordSys_Frame, 50.0, 90.0, Units_Frame, 3.6, buff); // same size as default axis titles
                                        TecUtilTextSetAnchor(textID, TextAnchor_Center);
                                        TecUtilTextSetMacroFunctionCmd(textID, "GBA.SP.Title");
                                    }
                                }
                                else
                                {
                                    // not sure if we need up update the frame here
                                    TecUtilFramePopByUniqueID(frame2d);
                                    EntIndex_t dstZoneNum = 1; // only one zone in each frame right now
                                    UniqueID_t result = updateDependentFrame(
                                        frame3d, zoneNum, subplot,
                                        frame2d, dstZoneNum,
                                        TRUE, // doProjection
                                        XYZ(probeX, probeY, probeZ), // projectionCenter
                                        XYZ(probeX, probeY, probeZ+1), // projectionTop
                                        FALSE/*useCurrentDataSet*/,
                                        FALSE/*copyCurrentStyle*/);
                                    CHECK(result == frame2d);
                                    // update the text too?
                                }
                                    
                                TecUtilFramePopByUniqueID(frame3d);
                            }
                            else
                                TecUtilDialogMessageBox("Zone has invalid values for CP nodes", MessageBoxType_Error);
                        }
                    }
                    else
                        TecUtilDialogMessageBox("Selected zone has invalid aux data for critical paths", MessageBoxType_Error);
                }

            }
            else
            {
                CHECK(stereoVisOp == StereoVisOp_PlotSingleZoneTest);

                EntIndex_t const subplot = 0; // TODO if probing 2D, get from frame

                UniqueID_t frame3d;
                UniqueID_t frame2d;
                if (TecUtilFrameGetPlotType() == PlotType_Cartesian2D)
                {
                    frame3d = INVALID_UNIQUE_ID;
                    frame2d = TecUtilFrameGetUniqueID();
                    // TODO: Check aux data to make sure this is a good frame
                }
                else
                {
                    frame3d = TecUtilFrameGetUniqueID();
                    frame2d = find2DFrameID(zoneNum, subplot);
                }

                if (frame2d == INVALID_UNIQUE_ID)
                {
                    frame2d = updateDependentFrame(
                        frame3d, zoneNum, subplot,
                        INVALID_UNIQUE_ID/*dstFrame:createsNewFrame*/, 0/*dstZone:zoneNumInvalid*/,
                        FALSE, // doProjection
                        XYZ(), // projectionCenter
                        XYZ(), // projectionTop
                        TRUE/*useCurrentDataSet*/,
                        TRUE/*copyCurrentStyle*/);
                    if (frame2d == INVALID_UNIQUE_ID)
                        TecUtilDialogErrMsg("Cannot find 2-D frame, and cannot create new frame.");
                }
                else
                {
                    TecUtilFramePopByUniqueID(frame2d);
                    EntIndex_t dstZoneNum = 1; // only one zone in each frame right now
                    UniqueID_t result = updateDependentFrame(
                        frame3d, zoneNum, subplot,
                        frame2d, dstZoneNum,
                        FALSE, // doProjection
                        XYZ(), // projectionCenter
                        XYZ(), // projectionTop
                        TRUE/*useCurrentDataSet*/,
                        TRUE/*copyCurrentStyle*/);
                    CHECK(result == frame2d);
                }

                TecUtilViewCopy();
                TecUtilFramePopByUniqueID(frame2d);
                TecUtilViewPaste();
                TecUtilViewFit();

                TecUtilFramePopByUniqueID(frame3d);
            }
        }
    }
    // not sure if this line should be inside if or outside
    installStereoVisProbeOverride(stereoVisOp);

    TecUtilLockFinish(AddOnID);
}

Boolean_t installStereoVisProbeOverride(StereoVisOp_e stereoVisOp)
{
    REQUIRE(VALID_ENUM(stereoVisOp,StereoVisOp_e));
    try
    {
        ArgList argList;
        argList.appendFunction(SV_CALLBACKFUNCTION, reinterpret_cast<void*>(stereoVisProbeCB));
        argList.appendString(SV_STATUSLINETEXT, "Click to visualize spheres as separate stereographic plots.");
        argList.appendArbParam(SV_CLIENTDATA, ArbParam_t(stereoVisOp));
        TecUtilProbeInstallCallbackX(argList.getRef());
        return TRUE;
    }
    catch (...)
    {
    	return FALSE;
    }
}


/*
**/
void stereoVis(StereoVisOp_e stereoVisOp)
{
    REQUIRE(VALID_ENUM(stereoVisOp, StereoVisOp_e) || stereoVisOp != StereoVisOp_PlotSingleZoneTest);
    REQUIRE(TecUtilLockIsOn());

    if ( TecUtilDataSetIsAvailable() && TecUtilFrameGetPlotType() == PlotType_Cartesian3D )
    {
        installStereoVisProbeOverride(stereoVisOp);
    }
    else
    {
        TecUtilDialogMessageBox("Stereographic visualizations require a 3-D plot.", MessageBox_Warning);
    }
}

/*
**/
void plotSingleZone()
{
    REQUIRE(TecUtilLockIsOn());

    if (TecUtilDataSetIsAvailable() && (TecUtilFrameGetPlotType() == PlotType_Cartesian3D || TecUtilFrameGetPlotType() == PlotType_Cartesian2D))
    {
        installStereoVisProbeOverride(StereoVisOp_PlotSingleZoneTest);
    }
    else
    {
        TecUtilDialogMessageBox("Select zones in 2-D or 3-D plot.", MessageBox_Warning);
    }
}
