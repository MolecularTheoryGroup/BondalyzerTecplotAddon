
#include <string>
#include <vector>

#include "TECADDON.h"
#include "ZONEVARINFO.h"

SmInteger_t VarNumByNameList(std::vector<std::string> &VarNameList)
{
	SmInteger_t VarNum = -1;

	EntIndex_t NumberOfVars = TecUtilDataSetGetNumVars();

	for (EntIndex_t CurrentVar = 1; CurrentVar <= NumberOfVars && VarNum < 0; ++CurrentVar){
		VarName_t CurrentVarName;
		if (TecUtilVarGetName(CurrentVar, &CurrentVarName)){
			for (int CurrentNameCheck = 0; CurrentNameCheck < VarNameList.size() && VarNum < 0; ++CurrentNameCheck){
				if (!strncmp(CurrentVarName, VarNameList[CurrentNameCheck].c_str(), VarNameList[CurrentNameCheck].size())){
					VarNum = CurrentVar;
				}
			}
			TecUtilStringDealloc(&CurrentVarName);
		}
	}

	return VarNum;
}

SmInteger_t VarNumByName(std::string &VarName)
{
	SmInteger_t VarNum = -1;

	EntIndex_t NumberOfVars = TecUtilDataSetGetNumVars();

	for (EntIndex_t CurrentVar = 1; CurrentVar <= NumberOfVars && VarNum < 0; ++CurrentVar){
		VarName_t CurrentVarName;
		if (TecUtilVarGetName(CurrentVar, &CurrentVarName)){
				if (!strncmp(CurrentVarName, VarName.c_str(), VarName.size())){
					VarNum = CurrentVar;
				}
			TecUtilStringDealloc(&CurrentVarName);
		}
	}

	return VarNum;
}

SmInteger_t ZoneNumByNameList(std::vector<std::string> &ZoneNameList)
{
	SmInteger_t ZoneNum = -1;

	EntIndex_t NumberOfZones = TecUtilDataSetGetNumZones();

	for (EntIndex_t CurrentVar = 1; CurrentVar <= NumberOfZones && ZoneNum < 0; ++CurrentVar){
		VarName_t CurrentZoneName;
		if (TecUtilZoneGetName(CurrentVar, &CurrentZoneName)){
			for (int CurrentNameCheck = 0; CurrentNameCheck < ZoneNameList.size() && ZoneNum < 0; ++CurrentNameCheck){
				if (!strncmp(CurrentZoneName, ZoneNameList[CurrentNameCheck].c_str(), ZoneNameList[CurrentNameCheck].size())){
					ZoneNum = CurrentVar;
				}
			}
			TecUtilStringDealloc(&CurrentZoneName);
		}
	}

	return ZoneNum;
}

SmInteger_t ZoneNumByName(std::string &ZoneName)
{
	SmInteger_t ZoneNum = -1;

	EntIndex_t NumberOfZones = TecUtilDataSetGetNumZones();

	for (EntIndex_t CurrentVar = 1; CurrentVar <= NumberOfZones && ZoneNum < 0; ++CurrentVar){
		VarName_t CurrentZoneName;
		if (TecUtilZoneGetName(CurrentVar, &CurrentZoneName)){
				if (!strncmp(CurrentZoneName, ZoneName.c_str(), ZoneName.size())){
					ZoneNum = CurrentVar;
				}
			TecUtilStringDealloc(&CurrentZoneName);
		}
	}

	return ZoneNum;
}