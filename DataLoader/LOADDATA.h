#ifndef LOADDATA_H_
#define LOADDATA_H_

#include <string>
#include <vector>

#include "CSM_DATA_LOADER_FUNCTIONS.h"

using std::string;
using std::vector;

enum LoadType{
	ADFFile = 0,
	BANDFile,
	VASPFile
};

int LoadADFTape41ASCIIData(char* FileName);
void GetADFASCIITape41FileName();

void NumCellsLoadData();
void SelVarsLoadData();

LoaderStatus_e GetTape41FileNames(Boolean_t IsADFFile);
LoaderStatus_e GetT41Contents(Boolean_t IsADFFile, Boolean_t SaveToList);
void LoadADFTape41Files();
void LoadADFTape41Data();
void LoadBANDTape41Data();
void QuitT41Load();

Boolean_t LoadADFTape21();

Boolean_t GetCHGCARFileName();
Boolean_t GetAECCARFileNames();
Boolean_t LoadVASPData();

void LoadGaussianCubeFiles();

void LoadTurboMoleCubeFiles();
Boolean_t LoadBinaryPLTFileData();

void SpinButtonInt(LgIndex_t DialogIndex, LgIndex_t Difference);
void SpinValueChangedInt(LgIndex_t DialogIndex);

void MakeDensfScriptForZones();
void ImportAdditionalTape41Files(Boolean_t MatchZones = TRUE, Boolean_t MatchDataSet = FALSE);

#endif