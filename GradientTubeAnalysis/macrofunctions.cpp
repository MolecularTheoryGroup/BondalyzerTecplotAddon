
#include <vector>
#include <string>
#include <sstream>
#include <fstream>

#include "TECADDON.h"
#include "ADDGLBL.h"
#include "CSM_DATA_TYPES.h"
#include "MACROFUNCTIONS.h"

#include <armadillo>
using namespace arma;




EntIndex_t MacroCreateZoneFromPolylines(std::vector<EntIndex_t> *ZoneNumList,
											Boolean_t ConnectStartToEnd,
											char* AddOnPath)
{
	EntIndex_t NewZoneNum = -1;

	std::stringstream MacroSS;

	MacroSS << "$!CREATEFESURFACEFROMIORDERED SOURCEZONES = [";

	for (EntIndex_t i = 0; i < ZoneNumList->size(); ++i){
		MacroSS << ZoneNumList->at(i);
		if (i < ZoneNumList->size() - 1) MacroSS << ",";
	}

	MacroSS << "] CONNECTSTARTTOEND = ";
	if (ConnectStartToEnd) 
		MacroSS << "YES";
	else
		MacroSS << "NO";

	Boolean_t MacroSuccess = TecUtilMacroExecuteCommand(MacroSS.str().c_str());

	if (MacroSuccess)
		NewZoneNum = TecUtilDataSetGetNumZones();

	return NewZoneNum;
}
Boolean_t MacroIntegrateVarByCellsOverZones(std::vector<ImportType_t> *ResultList,
												EntIndex_t ZoneNumStart, EntIndex_t ZoneNumEnd,
												EntIndex_t IntegratingVarNum,
												EntIndex_t XYZVarNums[3],
												const char* ResultsFileName)
{
	Boolean_t IsOk = (ZoneNumEnd > ZoneNumStart);

	std::stringstream MacroSS;

	if (IsOk)
		MacroSS << "$!EXTENDEDCOMMAND COMMANDPROCESSORID = 'CFDAnalyzer4' COMMAND = 'Integrate ["
		<< ZoneNumStart << "-" << ZoneNumEnd << "] VariableOption=\\'Scalar\\' XOrigin=0 YOrigin=0 ZOrigin=0 ScalarVar="
		<< IntegratingVarNum << " Absolute=\\'F\\' ExcludeBlanked=\\'F\\' XVariable="
		<< XYZVarNums[0] << " YVariable=" << XYZVarNums[1] << " ZVariable=" << XYZVarNums[2]
		<< " IntegrateOver=\\'Cells\\' IntegrateBy=\\'Zones\\' IRange={MIN = 1 MAX = 0 SKIP = 1} JRange={MIN = 1 MAX = 0 SKIP = 1} KRange={MIN = 1 MAX = 0 SKIP = 1} PlotResults=\\'T\\' PlotAs=\\'Result\\' TimeMin=0 TimeMax=0'";

	if (IsOk)
		IsOk = TecUtilMacroExecuteCommand(MacroSS.str().c_str());

	EntIndex_t FrameNum = 1;
	while (TecUtilFrameGetPlotType() != PlotType_XYLine && IsOk){
		if (FrameNum <= TecUtilFrameGetCount()){
			IsOk = TecUtilFrameActivateByNumber(FrameNum);
			++FrameNum;
		}
	}

	char* FrameName;

	if (IsOk)
		IsOk = TecUtilFrameGetName(&FrameName);
	
	EntIndex_t VarNum = -1;
	if (IsOk)
		VarNum = TecUtilVarGetNumByName("Result");

	EntIndex_t NumZones = TecUtilDataSetGetNumZones();
	for (EntIndex_t i = 1; i <= NumZones && IsOk; ++i){
		FieldData_pa VarDataRef = TecUtilDataValueGetReadableNativeRef(i, VarNum);
		if (VALID_REF(VarDataRef))
			ResultList->push_back(TecUtilDataValueGetByRef(VarDataRef, 1));
	}

	if (IsOk)
		IsOk = (ResultList->size() == NumZones && ResultList->size() == (ZoneNumEnd - ZoneNumStart + 1));

	if (IsOk)
		IsOk = TecUtilFrameDeleteActive();


	/*
	 *	This method saved the integration results to a temporary file, then read
	 *	in the file contents and parsed out the results. it gives low precision 
	 *	and this is a problem becuase the remaining space of a double variable when
	 *	a low precision number is converted from ascii to double is filled with
	 *	random numbers, so it was introducing random positive noise to the results.
	 */
// 	if (IsOk){
// 		MacroSS.clear();
// 		MacroSS.str(std::string());
// 		MacroSS << "$!EXTENDEDCOMMAND COMMANDPROCESSORID = 'CFDAnalyzer4' COMMAND = 'SaveIntegrationResults FileName=\\'"
// 			<< ResultsFileName << "\\''";
// 	}
// 
// 	IsOk = TecUtilMacroExecuteCommand(MacroSS.str().c_str());
// 
// 	std::ifstream ResultsFile;
// 
// 	if (IsOk){
// 		ResultsFile.open(ResultsFileName, std::ios::in);
// 		IsOk = (ResultsFile.is_open());
// 	}
// 
// 	if (IsOk){
// 		if (ResultList->size() > 0) 
// 			ResultList->clear();
// 		ResultsFile.ignore(500, '\n');
// 		std::string LineString;
// 		int ZoneRange = ZoneNumEnd - ZoneNumStart + 1;
// 		for (int i = 0; i < ZoneRange && !ResultsFile.eof(); ++i){
// 			getline(ResultsFile, LineString);
// 			LineString = LineString.substr(LineString.find_last_of(" ") + 1, LineString.length());
// 			ResultList->push_back(stod(LineString));
// 		}
// 		ResultsFile.close();
// 		std::remove(ResultsFileName);
// 	}

	return IsOk;
}
