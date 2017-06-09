#ifndef MACROFUNCTIONS_H_
#define MACROFUNCTIONS_H_




EntIndex_t MacroCreateZoneFromPolylines(std::vector<EntIndex_t> *ZoneNumList, 
										Boolean_t ConnectStartToEnd,
										char* AddOnPath);
Boolean_t MacroIntegrateVarByCellsOverZones(std::vector<ImportType_t> *ResultList,
											EntIndex_t ZoneNumStart, EntIndex_t ZoneNumEnd,
											EntIndex_t IntegratingVarNUm,
											EntIndex_t XYZVarNums[3],
											const char* ResultsFileName);

#endif