/*
******************************************************************
******************************************************************
*******                                                   ********
******             (C) 1988-2010 Tecplot, Inc.             *******
*******                                                   ********
******************************************************************
******************************************************************
*/

#ifndef ZONEVARINFO_H_
#define ZONEVARINFO_H_

#include <fstream>
#include <vector>
#include <string>

typedef struct _ZoneVarInfo_s
{
    EntIndex_t ZoneNum;
    Set_pa     RefinedZoneNums;
    EntIndex_t UVarNum;
    EntIndex_t VVarNum;
    EntIndex_t WVarNum;
	EntIndex_t DGradXNum;
	EntIndex_t DGradYNum;
	EntIndex_t DGradZNum;
	EntIndex_t XVarNum;
	EntIndex_t YVarNum;
	EntIndex_t ZVarNum;
    EntIndex_t ChrgDensVarNum;
	EntIndex_t GradMagVarNum;
    EntIndex_t TypeVarNum;
	LgIndex_t  NumVars;
    Boolean_t  PeriodicBC;

	// LOG
	std::ofstream*   BondLog;
	Boolean_t  LogBondInfo;

    double     XBeg;
    double     XEnd;
    double     YBeg;
    double     YEnd;
    double     ZBeg;
    double     ZEnd;
 
    ArrList_pa RefinedZoneXBeg;
    ArrList_pa RefinedZoneXEnd;
    ArrList_pa RefinedZoneYBeg;
    ArrList_pa RefinedZoneYEnd;
    ArrList_pa RefinedZoneZBeg;
    ArrList_pa RefinedZoneZEnd;

    // Temp var used for surface gradpath integration
    LgIndex_t LastElemUsed;
    SurfElemMap_pa SurfElemMap;
    Normals_pa     Normals;
    StreamDir_e    PathDir;
	LgIndex_t	   ElementCount;
	double		   AvgElemWidth;
	double		   Radius;
} ZoneVarInfo_s;

typedef struct _ZoneVarInfo_s  *ZoneVarInfo_pa;

SmInteger_t VarNumByNameList(std::vector<std::string> &VarNameList);
SmInteger_t VarNumByName(std::string &VarName);
SmInteger_t ZoneNumByNameList(std::vector<std::string> &ZoneNameList);
SmInteger_t ZoneNumByName(std::string &ZoneName);

Boolean_t ZoneVarInfoIsValid(ZoneVarInfo_pa ZoneVarInfo);

void ZoneVarInfoDealloc(ZoneVarInfo_pa *ZoneVarInfo);

void ZoneVarInfoClear(ZoneVarInfo_pa ZoneVarInfo);

ZoneVarInfo_pa ZoneVarInfoAlloc();

Boolean_t ZoneVarInfoSetMinMax(ZoneVarInfo_pa ZoneVarInfo,
                               Set_pa         RefinedZoneNums);
double ZoneVarInfoGetDXFromZoneNum(ZoneVarInfo_pa ZoneVarInfo,
                                   LgIndex_t      ZoneNum);

#endif /* ZONEVARINFO_H_ */
