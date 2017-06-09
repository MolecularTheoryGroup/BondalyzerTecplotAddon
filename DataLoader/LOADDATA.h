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

const int LoadADFTape41ASCIIData(char* FileName);
void GetADFASCIITape41FileName();

void NumCellsLoadData();
void SelVarsLoadData();

const LoaderStatus_e GetTape41FileNames(const Boolean_t & IsADFFile);
const LoaderStatus_e GetT41Contents(const Boolean_t & IsADFFile, const Boolean_t & SaveToList);
void LoadADFTape41Files();
void LoadADFTape41Data();
void LoadBANDTape41Data();
void QuitT41Load();

const Boolean_t LoadADFTape21();

const Boolean_t GetCHGCARFileName();
const Boolean_t GetAECCARFileNames();
const Boolean_t LoadVASPData();

void LoadGaussianCubeFiles();

void LoadTurboMoleCubeFiles();
Boolean_t LoadBinaryPLTFileData();

void SpinButtonInt(LgIndex_t DialogIndex, LgIndex_t Difference);
void SpinValueChangedInt(LgIndex_t DialogIndex);

void MakeDensfScriptForZones();
void ImportAdditionalTape41Files(const Boolean_t & MatchZones);

#endif