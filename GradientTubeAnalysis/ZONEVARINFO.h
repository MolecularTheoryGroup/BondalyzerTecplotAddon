#ifndef ZONEVARINFO_H_
#define ZONEVARINFO_H_

SmInteger_t VarNumByNameList(std::vector<std::string> &VarNameList);
SmInteger_t VarNumByName(std::string &VarName);
SmInteger_t ZoneNumByNameList(std::vector<std::string> &ZoneNameList);
SmInteger_t ZoneNumByName(std::string &ZoneName);

void LogDataValues(std::vector<double> &LevelVector, double Min, double Max, double BaseValue);
void X2DataValues(std::vector<double> &LevelVector, double Min, double Max);

#endif