
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

void LogDataValues(std::vector<double> &LevelVector, double Min, double Max, double BaseValue)
{
	if (Max - Min <= 0 && Min != 1e10)
	{
		LevelVector.push_back(0.0);
	}
	else
	{
		int LevelInt = 0, MultInt = 1, ExpInt = 0;
		double TempValue = 0;
		bool BreakLoop = false;
		TempValue = 1.0;
		if (Min != 1e10)
		{
			while (!BreakLoop)
			{
				if (TempValue > Min) TempValue *= 0.1;
				else if (TempValue * 10 <= Min) TempValue *= 10;
				else
				{
					BaseValue = TempValue;
					BreakLoop = true;
				}
			}
		}
		BreakLoop = false;
		while (!BreakLoop)
		{
			TempValue = MultInt * BaseValue * pow(10.0, ExpInt);
			if (TempValue >= Max)
			{
				TempValue = Max;
				BreakLoop = true;
			}
			LevelVector.push_back(TempValue);

			if (MultInt >= 9)
			{
				MultInt = 1;
				ExpInt++;
			}
			else MultInt++;
		}
	}
}

void X2DataValues(std::vector<double> &LevelVector, double Min, double Max)
{
	if (Max - Min <= 0 && Min != 1e10)
	{
		LevelVector.push_back(0.0);
	}
	else
	{
		LevelVector.push_back(Min);
		while (LevelVector.back() < Max){
			LevelVector.push_back(2 * LevelVector.back());
		}
	}
}