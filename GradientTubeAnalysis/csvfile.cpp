#include <fstream>
#include <vector>
#include <string>
#include <iomanip>

#include "TECADDON.h"
#include "CSM_DATA_TYPES.h"

#include "CSVFILE.h"

#include <armadillo>
using namespace arma;




Boolean_t CSVCreateFileAndInsertData(char* PathName, std::vector<std::vector<std::vector<ImportType_t> > > *Data, 
	std::vector<EntIndex_t> *IntegrationVarNums)
{
	Boolean_t IsOk = TRUE;
	char* DataSetNameCStr = NULL;
	TecUtilDataSetGetInfo(&DataSetNameCStr, NULL, NULL);
	std::ofstream OutFile;
	std::string DataSetFileNameStr = PathName;
	DataSetFileNameStr += "GTA_Results.csv";
	OutFile.open(DataSetFileNameStr.c_str(), std::ios::out | std::ios::app);
	IsOk = (OutFile.is_open());

	if (IsOk){

		OutFile << "\n\nDataset:," << DataSetNameCStr;
		TecUtilStringDealloc(&DataSetNameCStr);

		unsigned int NumWedgesTot = Data->at(0)[0].size();
		unsigned int NumWedgesHalfSlice = NumWedgesTot / 2;
		double WedgeDegStep = 360. / (double)NumWedgesTot;
		double WedgeRadStep = WedgeDegStep * DEG2RAD;

		unsigned int NumSlices = Data->at(0).size();
		unsigned int NumHalfSlices = 2 * NumSlices;
		double SliceDegStep = 180. / (double)NumSlices;
		double SliceRadStep = SliceDegStep * DEG2RAD;



		for (int VarNum = 0; VarNum < Data->size(); ++VarNum){
			char* VarNameCStr = NULL;
			IsOk = TecUtilVarGetName(IntegrationVarNums->at(VarNum), &VarNameCStr);

			if (VarNameCStr != NULL){
				OutFile << "\n\n\nVariable:" << VarNum + 1 << " of " << Data->size() << "," << VarNameCStr << "\n";
				TecUtilStringDealloc(&VarNameCStr);

				OutFile << "\n,,,Plane #";
				for (unsigned int i = 0; i < NumHalfSlices; ++i)
					OutFile << "," << i + 1;
				OutFile << "\n,,,Degrees";
				for (unsigned int i = 0; i < NumHalfSlices; ++i)
					OutFile << "," << (double)i * SliceDegStep;
				OutFile << "\n,,,Radians";
				for (unsigned int i = 0; i < NumHalfSlices; ++i)
					OutFile << "," << (double)i * SliceRadStep;

				OutFile << ",,Total";

				OutFile << "\nWedge,Degrees,Radians";


				std::vector<ImportType_t> SliceTotals(NumHalfSlices, 0.);

				for (unsigned int WedgeNum = 0; WedgeNum < NumWedgesHalfSlice; ++WedgeNum){
					ImportType_t WedgeTotal = 0;

					OutFile << "\n" << WedgeNum + 1 << "," << WedgeDegStep * WedgeNum << "," << WedgeRadStep * WedgeNum << ",";
					for (unsigned int RunNum = 0; RunNum < 2; ++RunNum){
						for (unsigned int SliceNum = 0; SliceNum < NumSlices; ++SliceNum){
							unsigned int RealSliceNum = RunNum * NumSlices + SliceNum;
							unsigned int RealWedgeNum;
							if (RunNum == 0)
								RealWedgeNum = WedgeNum;
							else
								RealWedgeNum = NumWedgesTot - WedgeNum - 1;

							SliceTotals[RealSliceNum] += Data->at(VarNum)[SliceNum][RealWedgeNum];
							WedgeTotal += Data->at(VarNum)[SliceNum][RealWedgeNum];
							OutFile << "," << std::setprecision(16) << Data->at(VarNum)[SliceNum][RealWedgeNum];
						}
					}
					OutFile << ",," << WedgeTotal;
				}
				OutFile << "\n\n,,Total,";
				for (unsigned int SliceNum = 0; SliceNum < SliceTotals.size(); ++SliceNum)
					OutFile << "," << std::setprecision(16) << SliceTotals[SliceNum];

				ImportType_t Total = 0;
				for (unsigned int SliceNum = 0; SliceNum < SliceTotals.size(); ++SliceNum)
					Total += SliceTotals[SliceNum];
				OutFile << ",," << std::setprecision(16) << Total;
			}
		}
		OutFile.close();
	}

	return IsOk;
}