

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <cstdio>
#include <fstream>
#include <vector>
#include <string>
#include <cstring>
#include <sstream>
#include <cctype>
#include <algorithm>

#include <omp.h>

// #include <windows.h>

#include "TECADDON.h"
#include "ADDGLBL.h"
#include "GUIDEFS.h"
#include "CSM_DATA_LOADER_FUNCTIONS.h"
#include "CSM_DATA_SET_INFO.h"
#include "KFc.h"
#include "CSM_DATA_SET_INFO.h"
#include "CSM_GRAD_PATH.h"
#include "CSM_CALC_VARS.h"
#include "LOADDATA.h"

#define BorhToAngstrom 0.529177249
#define AngstromToBohr 1.889725989

using std::vector;
using std::string;
using std::stringstream;
using std::endl;
using std::setprecision;
using std::ifstream;
using std::ofstream;
using std::streampos;
using std::ios;
using std::isdigit;
using std::to_string;
using std::getline;

vector<string> T41LoadUserVarStrings;
vector<string> T41LoadKFVarStrings;
vector<int> T41DataSizes;
KFFile GlobalKFFile;
string FileNameStr;
vector<string> FileNameStrs;

Boolean_t ReplaceDataSet, ImportAtoms;
int GroupIndex;

string WorkingDir;

LoadType FileType;

double MaxPointsForShowVolumeZone = 6e6;

string const DS_ZoneName_Delim = "_-_";
string const T41Ext = ".t41";

string const KFDelim = "%";
vector<FieldDataType_e> const KF_TP_Types = {
	FieldDataType_Int32,
	FieldDataType_Double,
	FieldDataType_Byte,
	FieldDataType_Int32
};
// 	{ KF_T_INTEGER, FieldDataType_Int32, sizeof(int)}, // int, 4 bytes
// 	{ KF_T_DOUBLE, FieldDataType_Double, sizeof(double)}, // double, 8 bytes
// 	{ KF_T_STRING, FieldDataType_Byte, sizeof(char)} // char, 1 byte, in KF files comes in groups of 160 chars
// 	{ KF_T_LOGICAL, FieldDataType_Int32, sizeof(int)}, // 4 bytes
// 	

void GetADFASCIITape41FileName(){
	TecUtilLockStart(AddOnID);

	if (TecUtilDataSetIsAvailable()){
		TecUtilDialogErrMsg("Cannot load file into existing dataset. Start with new file and try again.");
		TecUtilLockFinish(AddOnID);
		return;
	}

	FileNameStr.clear();

	char* FileName = nullptr;
	TecUtilDialogGetFileName(SelectFileOption_ReadSingleFile, &FileName, "ADF Tape41 ASCII", nullptr, "*.txt");
	if (FileName != nullptr)
		LoadADFTape41ASCIIData(FileName);
	TecUtilStringDealloc(&FileName);
	TecUtilLockFinish(AddOnID);
}

LoaderStatus_e GetTape41FileNames(Boolean_t IsADFFile){
	TecUtilLockStart(AddOnID);

// 	if (TecUtilDataSetIsAvailable()){
// 		TecUtilDialogErrMsg("Cannot load file into existing dataset. Start with new file and try again.");
// 		TecUtilLockFinish(AddOnID);
// 		return LoaderFailed;
// 	}

	if (IsADFFile)
		FileType = ADFFile;
	else
		FileType = BANDFile;


	/*
	*	Clear the vectors used to store the selected variables and
	*	section/variable strings
	*/
	T41DataSizes.clear();
	T41LoadKFVarStrings.clear();
	T41LoadUserVarStrings.clear();

	TecGUIListDeleteAllItems(MLT41LoadVa_MLST_D1);


	FileNameStr.clear();

	Boolean_t IsOk = TRUE;
	LoaderStatus_e Status = Loading;
	StringList_pa FileNames = TecUtilStringListAlloc();
	if (IsADFFile){
		IsOk = TecUtilDialogGetFileNames(SelectFileOption_AllowMultiFileRead, &FileNames, "ADF Tape41", nullptr, "*.t41");
	}
	else{
		IsOk = TecUtilDialogGetFileNames(SelectFileOption_AllowMultiFileRead, &FileNames, "BAND Tape41", nullptr, "*.t41");
	}
	if (IsOk && FileNames != nullptr){
		FileNameStrs.resize(TecUtilStringListGetCount(FileNames));
		for (int i = 1; i <= FileNameStrs.size() && Status != LoaderFailed; ++i){
			if (i == 1){
				if (TecUtilDataSetIsAvailable()){
					ReplaceDataSet = TecUtilDialogMessageBox("There is already data loaded. Do you want to replace the current data set with the new data?\n\tYes: Replace current data set\n\tNo: Append to current data set", MessageBoxType_YesNo);
				}
				else{
					ReplaceDataSet = TRUE;
				}
				if (FileNameStrs.size() > 1){
					ImportAtoms = TecUtilDialogMessageBox("Import atomic positions from all files(s)?\n(Click No to import atomic positions from first file only)", MessageBox_YesNo);
				}
				else
					ImportAtoms = TRUE;
			}

			FileNameStr = FileNameStrs[i - 1] = TecUtilStringListGetString(FileNames, i);
			if (i < FileNameStrs.size())
				Status = GetT41Contents(IsADFFile, FALSE);
			else
				Status = GetT41Contents(IsADFFile, TRUE);
		}
	}
	TecUtilStringListDealloc(&FileNames);

	TecUtilLockFinish(AddOnID);
	if (IsOk)
		return Status;
	else 
		return LoaderFailed;
}

Boolean_t IsNumber(string const & s){
	string::const_iterator it = s.begin();
	while (it != s.end() && isdigit(*it)) ++it;
	return (!s.empty() && it == s.end());
}

Boolean_t GetCHGCARFileName(){
	TecUtilLockStart(AddOnID);

	FileType = VASPFile;

	Boolean_t IsOk = TRUE;

	if (TecUtilDataSetIsAvailable()){
		TecUtilDialogErrMsg("Cannot load file into existing dataset. Start with new file and try again.");
		TecUtilLockFinish(AddOnID);
		IsOk = FALSE;
	}

	FileNameStrs.clear();

	if (IsOk){

		char* FileName = nullptr;
		TecUtilDialogGetFileName(SelectFileOption_ReadSingleFile, &FileName, "VASP CHGCAR", nullptr, "*CHGCAR");
		if (FileName != nullptr){
			FileNameStrs.push_back(FileName);
			TecUtilStringDealloc(&FileName);
		}
		else
			IsOk = FALSE;
	}

	TecUtilLockFinish(AddOnID);

	return IsOk;
}

Boolean_t GetAECCARFileNames(){
	TecUtilLockStart(AddOnID);

	vector<string> TmpStrs = {
		"VASP AECCAR0",
		"VASP AECCAR2"
	}; 
	vector<string> TmpStrs2 = {
		"AECCAR0",
		"AECCAR2"
	};

	Boolean_t IsOk = TRUE;

	if (TecUtilDataSetIsAvailable()){
		TecUtilDialogErrMsg("Cannot load file into existing dataset. Start with new file and try again.");
		TecUtilLockFinish(AddOnID);
		IsOk = FALSE;
	}

	FileNameStrs.clear();
	FileNameStr.clear();

	for (int i = 0; i < TmpStrs.size() && IsOk; ++i){

		FileNameStr.clear();

		if (IsOk){

			char* FileName = nullptr;
			TecUtilDialogGetFileName(SelectFileOption_ReadSingleFile, &FileName, TmpStrs[i].c_str(), nullptr, TmpStrs2[i].c_str());
			if (FileName != nullptr){
				FileNameStrs.push_back(FileName);
				TecUtilStringDealloc(&FileName);
			}
			else
				IsOk = FALSE;
		}
	}

	TecUtilLockFinish(AddOnID);

	return IsOk;
}


Boolean_t LoadVASPData(){

	TecUtilLockStart(AddOnID);

	SYSTEM_INFO sysinfo;
	GetSystemInfo(&sysinfo);

	int numCPU = sysinfo.dwNumberOfProcessors;
	omp_set_num_threads(numCPU);

	Boolean_t IsOk = TRUE;

	EntIndex_t VolZoneNum;
	LgIndex_t TotNumPoints = 0;


	double LatticeConstant = -1.0;
// 	mat33 LatticeVector = -ones<vec>(3,3);
	mat LatticeVector(3, 3);
	LatticeVector.fill(-1);
	Boolean_t HasAtomTypes = FALSE;

	for (int FileNum = 0; FileNum < FileNameStrs.size() && IsOk; ++FileNum){
		ifstream InFile(FileNameStrs[FileNum].c_str(), ios::in);
		if (!InFile.is_open()){
			TecUtilDialogErrMsg("Failed to open input file.");
			TecUtilLockFinish(AddOnID);
			return FALSE;
		}

		int const LineLen = 500;

		TecUtilDataLoadBegin();

		InFile.ignore(LineLen, '\n');

		/*
		 *	Read in the lattice constant factor
		 */
		InFile >> LatticeConstant;

		if (LatticeConstant < 0){
			TecUtilDialogErrMsg("Failed to parse lattice constant.");
			TecUtilLockFinish(AddOnID);
			return FALSE;
		}

		/*
		 *	Read in the lattice vector
		 */
		for (int i = 0; i < 3; ++i){
			for (int j = 0; j < 3; ++j){
				InFile >> LatticeVector.at(i, j);
				LatticeVector.at(i, j) *= LatticeConstant;
			}
		}

		int NumAtoms = 0;
		vector<AtomGroup_s> AtomGroupList;

		/*
		 *	Read in the total number of atoms, which is listed as several ints.
		 *	Don't stop until something other than an int is found.
		 */
		
		if (IsOk){
			vector<int> NumAtomList;
			string TmpStr;
			InFile >> TmpStr;

			if (!IsNumber(TmpStr)){
				/*
				*	VASP V5 file, so first get atom types,
				*	then number of each type. This will be
				*	used later to color them appropriately.
				*/
				HasAtomTypes = TRUE;
				vector<AtomColor_s> AtomColorList;
				PopulateAtomColorList(AtomColorList);
				while (true){
					if (!IsNumber(TmpStr)){
						AtomGroupList.emplace_back(TmpStr);
						AtomGroupList.back().AtomColor = GetAtomColor(AtomColorList, TmpStr);
						InFile >> TmpStr;
					}
					else
						break;
				}
			}

			int AtomNum = 1;
			while (true){
				if (IsNumber(TmpStr)){
					NumAtoms += stoi(TmpStr);
					NumAtomList.push_back(stoi(TmpStr));
					if (!HasAtomTypes)
						AtomGroupList.emplace_back(string(string("Atom " + to_string(AtomNum))));
					AtomGroupList[AtomNum - 1].Count = stoi(TmpStr);
					InFile >> TmpStr;
				}
				else
					break;
				AtomNum++;
			}

			if (AtomGroupList.size() <= 0 && NumAtomList.size() <= 0 && AtomGroupList.size() != NumAtomList.size()){
					TecUtilDialogErrMsg("Failed to parse VASP file: Atom types");
					return FALSE;
			}
		}

		/*
		 *	For storing the atomic positions within the cell.
		 *	These will be used to place zones of 1 point for each atom.
		 *	to assist with recognizing the system pre-bondalyzed,
		 *	for identifying non-nuclear maxima, and for easing the
		 *	process of preparing images for papers etc.
		 */

		for (int i = 0; i < AtomGroupList.size(); ++i){
			for (int j = 0; j < AtomGroupList[i].Count; ++j){
				vec3 TmpPoint;
				for (int k = 0; k < 3; ++k)
					InFile >> TmpPoint[k];
				AtomGroupList[i].AddPosition(TmpPoint);
			}
		}

		if (HasAtomTypes){
			/*
			*	There could be more than one group of the same
			*	atom type, so combine groups of the same type.
			*/
			vector<AtomGroup_s> NewAtomGroupList;
			vector<Boolean_t> WasCombined(AtomGroupList.size(), FALSE);

			NewAtomGroupList.push_back(AtomGroupList[0]);
			WasCombined[0] = TRUE;

			for (int i = 1; i < AtomGroupList.size(); ++i){
				Boolean_t IsFound = FALSE;
				for (int j = 0; j < NewAtomGroupList.size() && !IsFound; ++j){
					if (!WasCombined[i] && AtomGroupList[i].Name == NewAtomGroupList[j].Name){
						IsFound = TRUE;
						WasCombined[i] = TRUE;
						NewAtomGroupList[j].Count += AtomGroupList[i].Count;
						for (int k = 0; k < 3; ++k)
							NewAtomGroupList[j].Positions[k].insert(
								NewAtomGroupList[j].Positions[k].begin(),
								AtomGroupList[i].Positions[k].cbegin(),
								AtomGroupList[i].Positions[k].cend()
							);
					}
				}
				if (!IsFound){
					NewAtomGroupList.push_back(AtomGroupList[i]);
				}
			}

			AtomGroupList = NewAtomGroupList;
		}

		InFile.ignore(LineLen, '\n');

		LgIndex_t Mx, My, Mz, NumPts;

		stringstream ss;
		string Title1, Title2;

		InFile.ignore(LineLen, '\n');
		getline(InFile, Title1);

		ss << Title1;
		ss >> Mx >> My >> Mz;

		NumPts = Mx * My * Mz;

		double BlankValue = 1.2345e50;

		int ArrayMemoryKB = sizeof(double) * NumPts / 1024;
		TecUtilMemoryChangeNotify(ArrayMemoryKB);

		vector<vector<vector<double> > > Charge(Mz, vector<vector<double> >(My, vector<double>(Mx, BlankValue)));
		vector<vector<vector<vector<double> > > > SpinFields, TmpCharge;
		vector<vector<vector<vector<vector<double> > > > > TmpSpinFields;

		StatusLaunch("Reading data...", AddOnID, TRUE);

		/*
		 *	Now read in the block of ASCII data that is the
		 *	charge density for the entire system.
		 *	Read in the entire file as a vector of char,
		 *	and then parse through that to get the actual 
		 *	values.
		 */

		/*
		 *	Since we'll read in the whole file, need to
		 *	save the first value in the block so we can
		 *	find the starting point in the char array.
		 *	The first value of the block happens to be the
		 *	next string in the input file
		 */
		string NextValue;
		InFile >> NextValue;

		//	Get file size
		InFile.seekg(0, InFile.end);
		std::streamoff FileSize = InFile.tellg();
		InFile.seekg(0, InFile.beg);

		//	Read in entire file to FileData
		vector<char> FileData(FileSize);
		InFile.read(FileData.data(), FileSize);
		InFile.close();

		//	These are for updating progress bar
// 		int NumPtsTotal = NumPts * static_cast<int>(FileNameStrs.size());
// 		int CurPt = NumPts * FileNum;
		int MzTotal = Mz * static_cast<int>(FileNameStrs.size());
		int CurK = Mz * FileNum;

		/*
		 *	Get the first word from the file, then
		 *	loop over words until first data point
		 *	in the block of data is found again
		 */
		char* TmpCStr = strtok(FileData.data(), " \n\t");
		while (TmpCStr != nullptr && NextValue.compare(TmpCStr))
			TmpCStr = strtok(nullptr, " \n\t");

		/*
		 *	Now simply parse through the char array, saving each
		 *	value to its place in the array of charge density
		 */

// 		for (int i = 0; i < NumPts; ++i){
// 			if (!SetPercent(i + CurPt, NumPtsTotal, "Reading data... Charge Density", AddOnID)){
// 				TecUtilDialogDropPercentDone();
// 				TecUtilMemoryChangeNotify(-ArrayMemoryKB);
// 				TecUtilDataLoadEnd();
// 				TecUtilLockFinish(AddOnID);
// 				return FALSE;
// 			}
// 			Charge[i] = atof(TmpCStr);
// 			TmpCStr = strtok(nullptr, " \n\t");
// 		}

		for (int k = 0; k < Mz; ++k){
			if (!StatusUpdate(k + CurK, MzTotal, "Reading data... Charge Density", AddOnID)){
				StatusDrop(AddOnID);
				TecUtilMemoryChangeNotify(-ArrayMemoryKB);
				TecUtilDataLoadEnd();
				TecUtilLockFinish(AddOnID);
				return FALSE;
			}
			for (int j = 0; j < My; ++j){
				for (int i = 0; i < Mx && TmpCStr != nullptr; ++i){
					Charge[k][j][i] = atof(TmpCStr);
					TmpCStr = strtok(nullptr, " \n\t");
				}
			}
		}

		StatusDrop(AddOnID);

// 		if (Charge[NumPts - 1] == BlankValue){
		if (Charge[Mz - 1][My - 1][Mx - 1] == BlankValue){
			TecUtilDialogErrMsg("CHGCAR file ended sooner than expected. Import failed.");
			TecUtilMemoryChangeNotify(-ArrayMemoryKB);
			TecUtilDataLoadEnd();
			TecUtilLockFinish(AddOnID);
			return FALSE;
		}


		/*
		 *	If another block exists, then it represents
		 *	the difference between the total charge 
		 *	and the charge of each spin density.
		 *	If this is found, parse the char array
		 *	and save the values as with the charge density.
		 */
		int SpinNum = 0;
		while (TmpCStr != nullptr){
			if (!Title1.compare(TmpCStr)){
				TecUtilMemoryChangeNotify(ArrayMemoryKB);
// 				Diff.resize(NumPts, BlankValue);
				SpinFields.emplace_back(Mz, vector<vector<double> >(My, vector<double>(Mx, BlankValue)));
				StatusLaunch("Reading data...", AddOnID, TRUE);
				TmpCStr = strtok(nullptr, " \n\t");

// 				for (int i = 0; i < NumPts; ++i){
// 					if (!SetPercent(i + CurPt, NumPtsTotal, "Reading data... Spin Density", AddOnID)){
// 						TecUtilDialogDropPercentDone();
// 						TecUtilMemoryChangeNotify(-ArrayMemoryKB);
// 						TecUtilDataLoadEnd();
// 						TecUtilLockFinish(AddOnID);
// 						return FALSE;
// 					}
// 					Diff[i] = atof(TmpCStr);
// 					TmpCStr = strtok(nullptr, " \n\t");
// 				}

				for (int k = 0; k < Mz; ++k){
					if (!StatusUpdate(k + CurK, MzTotal, "Reading data... Spin Density " + to_string(SpinFields.size()), AddOnID)){
						StatusDrop(AddOnID);
						TecUtilMemoryChangeNotify(-2 * ArrayMemoryKB);
						TecUtilDataLoadEnd();
						TecUtilLockFinish(AddOnID);
						return FALSE;
					}
					for (int j = 0; j < My; ++j){
						for (int i = 0; i < Mx && TmpCStr != nullptr; ++i){
							SpinFields.back()[k][j][i] = atof(TmpCStr);
							TmpCStr = strtok(nullptr, " \n\t");
						}
					}
				}

				StatusDrop(AddOnID);

// 				if (Diff[NumPts - 1] == BlankValue){
				if (SpinFields.back()[Mz - 1][My - 1][Mx - 1] == BlankValue){
					TecUtilDialogErrMsg("CHGCAR file ended sooner than expected. Import failed.");
					TecUtilMemoryChangeNotify(-2*ArrayMemoryKB);
					TecUtilDataLoadEnd();
					TecUtilLockFinish(AddOnID);
					return FALSE;
				}
			}
			TmpCStr = strtok(nullptr, "\n");
		}

		/*
		 *	Use the lattice vectors to find out the volume of one cell.
		 *	then use that cell to convert charge density for the
		 *	correct volume.
		 */
		double V1, V2, V3, Volume;

		V1 = LatticeVector.at(0, 0) * (LatticeVector.at(1, 1) * LatticeVector.at(2, 2) - LatticeVector.at(1, 2) * LatticeVector.at(2, 1));
		V2 = LatticeVector.at(0, 1) * (LatticeVector.at(1, 2) * LatticeVector.at(2, 0) - LatticeVector.at(1, 0) * LatticeVector.at(2, 2));
		V3 = LatticeVector.at(0, 2) * (LatticeVector.at(1, 0) * LatticeVector.at(2, 1) - LatticeVector.at(1, 1) * LatticeVector.at(2, 0));
		Volume = V1 + V2 + V3;

#pragma omp parallel for
		for (int k = 0; k < Mz; ++k){
			for (int j = 0; j < My; ++j){
				for (int i = 0; i < Mx; ++i){
					Charge[k][j][i] /= Volume;
					for (auto & Spin : SpinFields)
						Spin[k][j][i] /= Volume;
				}
			}
		}

		/*
		 *	If there are multiple files (such as when using
		 *	AECCAR files), then need to save the contents read
		 *	from the first file(s) for later.
		 *	If there's only 1 file, or the current file is
		 *	the last of multiple files, then sum the values
		 *	read from multiple files if necessary and move on
		 *	to load the data into Tec360.
		 */
		if (FileNum < FileNameStrs.size() - 1){
			TmpCharge.push_back(Charge);
			Charge.clear();
			TmpSpinFields.push_back(SpinFields);
			SpinFields.clear();

		}
		else{
			if (FileNameStrs.size() > 1){
				/*
				 *	Add Charge to TmpCharge
				 */


				for (int ChargeNum = 0; ChargeNum < TmpCharge.size(); ++ChargeNum){
					if (TmpCharge[ChargeNum].size() * TmpCharge[ChargeNum][0].size() * TmpCharge[ChargeNum][0][0].size() != NumPts){
						TecUtilDialogErrMsg("Files do not have same system size. Cancelling import.");
						TecUtilMemoryChangeNotify(-ArrayMemoryKB);
						for (auto const & i : TmpSpinFields[ChargeNum])
							TecUtilMemoryChangeNotify(-ArrayMemoryKB);
						TecUtilDataLoadEnd();
						TecUtilLockFinish(AddOnID);
						return FALSE;
					}
					

#pragma omp parallel for
					for (int k = 0; k < Mz; ++k){
						for (int j = 0; j < My; ++j){
							for (int i = 0; i < Mx; ++i){
								Charge[k][j][i] += TmpCharge[ChargeNum][k][j][i];
								for (int iSpin = 0; iSpin < TmpSpinFields[ChargeNum].size(); ++iSpin) {
									SpinFields[iSpin][k][j][i] += TmpSpinFields[ChargeNum][iSpin][k][j][i];
								}
							}
						}
					}

// 					if (!IsPolar && TmpSpinFields.size() > 0 &&
// 						TmpSpinFields[ChargeNum].size() == Mz &&
// 						TmpSpinFields[ChargeNum][0].size() == My &&
// 						TmpSpinFields[ChargeNum][0][0].size() == Mx){
// 						Diff = TmpSpinFields[ChargeNum];
// 						IsPolar = TRUE;
// 					}
				}
			}

			/*
			*	Begin "expand" subroutine.
			*	If the user elected to copy the system 
			*	in x,y,z direction, need to perform those 
			*	calculations.
			*	Actually, this part isn't necessary, as the 
			*	expansion can take place as the data is loaded
			*	into Tec360, so here just getting the number of 
			*	cells to expand in x,y,z direction.
			*/
			int NCells = 1;
			int NCellsX = NCells, NCellsY = NCells, NCellsZ = NCells;

// 			TecGUITextFieldGetLgIndex(XNC_TFS_D2, &NCellsX);
// 			TecGUITextFieldGetLgIndex(YNC_TFS_D2, &NCellsY);
// 			TecGUITextFieldGetLgIndex(ZNC_TFS_D2, &NCellsZ);

			TecUtilDataLoadEnd();

			/*
			 *	begin write_chgcar subroutine equivalent.
			 *	Obtain raw pointers to the x,y,z,rho variables,
			 *	and write the data into tecplot.
			 */

			StringList_pa VarNames = TecUtilStringListAlloc();

			TecUtilStringListAppendString(VarNames, "X");
			TecUtilStringListAppendString(VarNames, "Y");
			TecUtilStringListAppendString(VarNames, "Z");

			int NumVars = 4;
			if (SpinFields.size() == 1){
				NumVars = 6;
				TecUtilStringListAppendString(VarNames, "Electron Density Sum");
				TecUtilStringListAppendString(VarNames, "Electron Density A");
				TecUtilStringListAppendString(VarNames, "Electron Density B");
			}
			else if (SpinFields.size() == 3){
				NumVars = 12;
				TecUtilStringListAppendString(VarNames, "Electron Density total");
				TecUtilStringListAppendString(VarNames, "Electron Density (total + x) / 2");
				TecUtilStringListAppendString(VarNames, "Electron Density (total - x) / 2");
				TecUtilStringListAppendString(VarNames, "Electron Density (total + y) / 2");
				TecUtilStringListAppendString(VarNames, "Electron Density (total - y) / 2");
				TecUtilStringListAppendString(VarNames, "Electron Density (total + z) / 2");
				TecUtilStringListAppendString(VarNames, "Electron Density (total - z) / 2");
				TecUtilStringListAppendString(VarNames, "Electron Density (total + x + y + z) / 6");
				TecUtilStringListAppendString(VarNames, "Electron Density (total - x - y - z) / 6");
			}
			else{
				TecUtilStringListAppendString(VarNames, "Electron Density");
			}

			LgIndex_t MXYZ[3] = { Mx, My, Mz};

			LgIndex_t TMx = Mx * NCellsX;
			LgIndex_t TMy = My * NCellsY;
			LgIndex_t TMz = Mz * NCellsZ;
			TotNumPoints = TMx * TMy * TMz;
			vector<FieldDataType_e> VarDataTypes;
			VarDataTypes.resize(NumVars, FieldDataType_Float);


			if (TecUtilDataSetCreate("CHGCAR", VarNames, TRUE))
				IsOk = TecUtilDataSetAddZone(StringMakeValidZoneName(CSMZoneName.FullVolume + "CHGCAR").c_str(), TMx, TMy, TMz, ZoneType_Ordered, VarDataTypes.data());

			if (IsOk)
				VolZoneNum = TecUtilDataSetGetNumZones();

			SetZoneStyle({ VolZoneNum }, ZoneStyle_Volume);

			vector<FieldDataPointer_c> VarRawPtrs(NumVars);

			TecUtilDataLoadBegin();

			StatusLaunch("Loading data into dataset...", AddOnID, TRUE);

			for (EntIndex_t VarNum = 0; VarNum < NumVars; ++VarNum)
				VarRawPtrs[VarNum].InitializeWritePtr(VolZoneNum, VarNum + 1);

			Boolean_t TaskQuit = FALSE;

			int TotPtsStatus = Mz * NCellsZ / numCPU;

			for (int Zz = 0; Zz < NCellsZ; ++Zz){
#pragma omp parallel for
				for (LgIndex_t k = 1; k <= Mz; ++k){
					if (omp_get_thread_num() == 0 && !StatusUpdate(k + Mz * Zz, TotPtsStatus, "Loading data into dataset...", AddOnID)){
						TaskQuit = TRUE;
#pragma omp flush (TaskQuit)
					}
#pragma omp flush (TaskQuit)
					for (int Zy = 0; Zy < NCellsY && !TaskQuit; ++Zy){
						for (LgIndex_t j = 1; j <= My; ++j){
							for (int Zx = 0; Zx < NCellsX; ++Zx){
								for (LgIndex_t i = 1; i <= Mx; ++i){
									LgIndex_t IJK[3] = { i - 1, j - 1, k - 1 };
									LgIndex_t TIJK[3] = { i + Zx*Mx, j + Zy*My, k + Zz*Mz };
									LgIndex_t TotIndex = IndexFromIJK(i + Zx*Mx, j + Zy*My, k + Zz*Mz, TMx, TMy, TMz, TRUE) - 1;
									for (EntIndex_t VarNum = 0; VarNum < NumVars; ++VarNum){
										double PointValue;
										if (VarNum < 3){
											/*
											 *	Here, for the generation of x,y,z values, using
											 *	the lattice vectors so that non-cartesian systems
											 *	are property represented as well.
											 */
											PointValue = 0.0;
											for (int Dir = 0; Dir < 3; ++Dir)
												PointValue += static_cast<double>(TIJK[Dir]) / static_cast<double>(MXYZ[Dir]) * LatticeVector.at(Dir, VarNum);
										}
										else{
											PointValue = Charge[k - 1][j - 1][i - 1];
											if (VarNum >= 4){
												if (SpinFields.size() == 1) {
													double PtDiff = SpinFields[0][k - 1][j - 1][i - 1];
													switch (VarNum - 3) {
													case 1:
														PointValue = (PointValue + PtDiff) / 2.0;
														break;
													default:
														PointValue = (PointValue - PtDiff) / 2.0;
														break;
													}
												}
												else if (SpinFields.size() == 3) {
													double SpinVal;
													switch (VarNum - 3) {
													case 1: // (total + x) / 2
														SpinVal = SpinFields[0][k - 1][j - 1][i - 1];
														PointValue = (PointValue + SpinVal) / 2.0;
														break;
													case 2: // (total - x) / 2
														SpinVal = SpinFields[0][k - 1][j - 1][i - 1];
														PointValue = (PointValue - SpinVal) / 2.0;
														break;
													case 3: // (total + y) / 2
														SpinVal = SpinFields[1][k - 1][j - 1][i - 1];
														PointValue = (PointValue + SpinVal) / 2.0;
														break;
													case 4: // (total - y) / 2
														SpinVal = SpinFields[1][k - 1][j - 1][i - 1];
														PointValue = (PointValue - SpinVal) / 2.0;
														break;
													case 5: // (total + z) / 2
														SpinVal = SpinFields[2][k - 1][j - 1][i - 1];
														PointValue = (PointValue + SpinVal) / 2.0;
														break;
													case 6: // (total - z) / 2
														SpinVal = SpinFields[2][k - 1][j - 1][i - 1];
														PointValue = (PointValue - SpinVal) / 2.0;
														break;
													case 7: // (total + x + y + z) / 6
														for (auto const & Spin : SpinFields){
															PointValue += Spin[k - 1][j - 1][i - 1];
														}
														PointValue /= 6.0;
														break;
													case 8: // (total - x - y - z) / 6
														for (auto const & Spin : SpinFields) {
															PointValue -= Spin[k - 1][j - 1][i - 1];
														}
														PointValue /= 6.0;
														break;
													}
												}
											}
										}

										VarRawPtrs[VarNum].Write(TotIndex, PointValue);
									}
								}
							}
						}
					}
				}
			}


			StatusDrop(AddOnID);
			if (TaskQuit){
				TecUtilMemoryChangeNotify(-ArrayMemoryKB);
				for (int i = 0; i < SpinFields.size(); ++i)
					TecUtilMemoryChangeNotify(-ArrayMemoryKB);
				TecUtilDataLoadEnd();
				TecUtilLockFinish(AddOnID);
				return FALSE;
			}


			TecUtilDataLoadEnd();

			if (IsOk){
				Set_pa TempSet = TecUtilSetAlloc(FALSE);
				IsOk = TecUtilSetAddMember(TempSet, VolZoneNum, FALSE);
				if (IsOk){
					TecUtilZoneSetActive(TempSet, AssignOp_PlusEquals);
				}
				TecUtilSetDealloc(&TempSet);
			}

			TecUtilFrameSetPlotType(PlotType_Cartesian3D);

			Set_pa ChangedVars = TecUtilSetAlloc(FALSE);
			for (EntIndex_t i = 1; i <= NumVars; ++i)
				TecUtilSetAddMember(ChangedVars, i, FALSE);
			TecUtilStateChanged(StateChange_VarsAltered, (ArbParam_t)ChangedVars);
			TecUtilSetDealloc(&ChangedVars);

			/*
			 *	Now add ordered zones of scatter spheres for each atom
			 *	type. Two considerations: system may be repeated for
			 *	multiple cells, and atoms may lie at the edge of periodic
			 *	boundaries. So this is written to allow both of those things.
			 */

			/*
			 *	Get system boundaries for periodic atoms
			 */
			double XYZMin[3], XYZMax[3];

			for (int i = 0; i < 3; ++i)
				TecUtilVarGetMinMax(i + 1, &XYZMin[i], &XYZMax[i]);

			int NCellsXYZ[3] = { NCellsX, NCellsY, NCellsZ };

			FieldDataPointer_c RhoPtr;
			VolExtentIndexWeights_s VolInfo;
			if (VarNumByName("Electron Density", true) > 0) {
				RhoPtr.InitializeReadPtr(TecUtilDataSetGetNumZones(), VarNumByName("Electron Density", true));
			}
			else if (VarNumByName("rho", true) > 0) {
				RhoPtr.InitializeReadPtr(TecUtilDataSetGetNumZones(), VarNumByName("rho", true));
			}
			else if (VarNumByName("denstiy", true) > 0) {
				RhoPtr.InitializeReadPtr(TecUtilDataSetGetNumZones(), VarNumByName("density", true));
			}

			if (RhoPtr.IsReady()) {
				GetVolInfo(TecUtilDataSetGetNumZones(), { 1,2,3 }, FALSE, VolInfo);
			}

			mat33 FractionalLattice = LatticeVector;
			FractionalLattice.row(0) /= (Mx - 1.);
			FractionalLattice.row(1) /= (My - 1.);
			FractionalLattice.row(2) /= (Mz - 1.);
			mat33 LVTranspose = LatticeVector.t();
			int ColorCount = 0;
			for (int GroupNum = 0; GroupNum < AtomGroupList.size() && IsOk; ++GroupNum) {
				vector<vector<ImportType_t> > XYZ(3, vector<ImportType_t>());
				vector<ImportType_t> Charges;
				for (int AtomNum = 0; AtomNum < AtomGroupList[GroupNum].Count && IsOk; ++AtomNum) {
					vec3 TmpLatPos;
					for (int i = 0; i < 3; ++i)
						TmpLatPos[i] = AtomGroupList[GroupNum].Positions[i][AtomNum];

					TmpLatPos = LatticeVector * TmpLatPos;;

					if (AtomGroupList[GroupNum].Charges.size() <= AtomNum) {
						AtomGroupList[GroupNum].Charges.push_back(16);
					}
					if (RhoPtr.IsReady()) {
						AtomGroupList[GroupNum].Charges[AtomNum] = MAX(16, RhoPtr.At(TmpLatPos, VolInfo));
					}

					vec3 CPTmpXYZ;

					for (int kk = 0; kk < NCellsZ + int(abs(TmpLatPos[2]) < 0.01); ++kk) {
						for (int jj = 0; jj < NCellsY + int(abs(TmpLatPos[1]) < 0.01); ++jj) {
							for (int ii = 0; ii < NCellsX + int(abs(TmpLatPos[0]) < 0.01); ++ii) {
								vec3 IJK;
								IJK << ii * Mx << jj * My << kk * Mz;

								vec3 Pos = (FractionalLattice * IJK) + TmpLatPos;
								for (int dir = 0; dir < 3; ++dir) {
									XYZ[dir].push_back(Pos[dir]);
								}
								Charges.push_back(AtomGroupList[GroupNum].Charges[AtomNum]);
								if (NCellsX == 1) {
									break;
								}
							}
							if (NCellsY == 1) {
								break;
							}
						}
						if (NCellsZ == 1) {
							break;
						}
					}
				}
				// Now copy charge and position to the AtomGroupList
				if (!Charges.empty()) {
					// 			AtomGroupList[GroupNum].Charges.insert(AtomGroupList[GroupNum].Charges.end(), Charges.begin(), Charges.end());
					AtomGroupList[GroupNum].Charges = Charges;
					for (int i = 0; i < 3; ++i) {
						// 				AtomGroupList[GroupNum].Positions[i].insert(AtomGroupList[GroupNum].Positions[i].end(), XYZ[i].begin(), XYZ[i].end());
						AtomGroupList[GroupNum].Positions[i] = XYZ[i];
					}
					AtomGroupList[GroupNum].Count = Charges.size();
				}
			}

			CreateAtomZonesFromAtomGroupList(AtomGroupList, { "X","Y","Z" }, VarDataTypes, 50);

			ChangedVars = TecUtilSetAlloc(FALSE);
			for (EntIndex_t i = 1; i <= NumVars; ++i)
				TecUtilSetAddMember(ChangedVars, i, FALSE);
			TecUtilStateChanged(StateChange_VarsAltered, (ArbParam_t)ChangedVars);
			TecUtilSetDealloc(&ChangedVars);

			/*
			 *	Set style info for the volume zone.
			 *	Specifically, turn off everything you don't need
			 *	to see, and turn on an isosurface so that the
			 *	user can see that the atoms match up with the 
			 *	charge density, basically so that the user can 
			 *	see the basic layout of the system.
			 */
			ArgList_pa argList = TecUtilArgListAlloc();

			TecUtilArgListClear(argList);
			TecUtilArgListAppendString(argList, SV_P1, SV_ISOSURFACEATTRIBUTES);
			TecUtilArgListAppendString(argList, SV_P2, SV_EFFECTS);
			TecUtilArgListAppendString(argList, SV_P3, SV_USETRANSLUCENCY);
			TecUtilArgListAppendInt(argList, SV_OFFSET1, 1);
			TecUtilArgListAppendArbParam(argList, SV_IVALUE, TRUE);
			TecUtilStyleSetLowLevelX(argList);

			TecUtilArgListClear(argList);
			TecUtilArgListAppendString(argList, SV_P1, SV_ISOSURFACEATTRIBUTES);
			TecUtilArgListAppendString(argList, SV_P2, SV_EFFECTS);
			TecUtilArgListAppendString(argList, SV_P3, SV_SURFACETRANSLUCENCY);
			TecUtilArgListAppendInt(argList, SV_OFFSET1, 1);
			TecUtilArgListAppendArbParam(argList, SV_IVALUE, 70);
			TecUtilStyleSetLowLevelX(argList);

			TecUtilArgListClear(argList);
			TecUtilArgListAppendString(argList, SV_P1, SV_FIELDLAYERS);
			TecUtilArgListAppendString(argList, SV_P2, SV_SHOWSCATTER);
			TecUtilArgListAppendArbParam(argList, SV_IVALUE, TRUE);
			TecUtilStyleSetLowLevelX(argList);

			TecUtilArgListClear(argList);
			TecUtilArgListAppendString(argList, SV_P1, SV_ISOSURFACELAYERS);
			TecUtilArgListAppendString(argList, SV_P2, SV_SHOW);
			TecUtilArgListAppendArbParam(argList, SV_IVALUE, TRUE);
			TecUtilStyleSetLowLevelX(argList);

			TecUtilArgListClear(argList);
			TecUtilArgListAppendString(argList, SV_P1, SV_ISOSURFACEATTRIBUTES);
			TecUtilArgListAppendString(argList, SV_P2, SV_DEFINITIONCONTOURGROUP);
			TecUtilArgListAppendInt(argList, SV_OFFSET1, 1);
			TecUtilArgListAppendArbParam(argList, SV_IVALUE, (SmInteger_t)1);
			TecUtilStyleSetLowLevelX(argList);

			TecUtilArgListClear(argList);
			TecUtilArgListAppendString(argList, SV_P1, SV_ISOSURFACEATTRIBUTES);
			TecUtilArgListAppendString(argList, SV_P2, SV_SHOWGROUP);
			TecUtilArgListAppendInt(argList, SV_OFFSET1, 1);
			TecUtilArgListAppendArbParam(argList, SV_IVALUE, TRUE);
			TecUtilStyleSetLowLevelX(argList);

			EntIndex_t IsoSurfaceVarNum = 4;
			double IsoValue = 0.6;

			TecUtilArgListClear(argList);
			TecUtilArgListAppendString(argList, SV_P1, SV_ISOSURFACEATTRIBUTES);
			TecUtilArgListAppendString(argList, SV_P2, SV_ISOVALUE1);
			TecUtilArgListAppendInt(argList, SV_OFFSET1, 1);
			TecUtilArgListAppendDouble(argList, SV_DVALUE, IsoValue);
			TecUtilStyleSetLowLevelX(argList);


			TecUtilArgListClear(argList);
			TecUtilArgListAppendString(argList, SV_P1, SV_GLOBALCONTOUR);
			TecUtilArgListAppendString(argList, SV_P2, SV_VAR);
			TecUtilArgListAppendInt(argList, SV_OFFSET1, 1);
			TecUtilArgListAppendArbParam(argList, SV_IVALUE, IsoSurfaceVarNum);
			TecUtilStyleSetLowLevelX(argList);
			TecUtilArgListClear(argList);
			TecUtilArgListAppendString(argList, SV_P1, SV_ISOSURFACEATTRIBUTES);
			TecUtilArgListAppendString(argList, SV_P2, SV_ISOVALUE1);
			TecUtilArgListAppendInt(argList, SV_OFFSET1, 1);
			TecUtilArgListAppendDouble(argList, SV_DVALUE, IsoValue);
			TecUtilStyleSetLowLevelX(argList);
			TecUtilArgListClear(argList);
			TecUtilArgListAppendString(argList, SV_P1, SV_ISOSURFACEATTRIBUTES);
			TecUtilArgListAppendString(argList, SV_P2, SV_ISOVALUE2);
			TecUtilArgListAppendInt(argList, SV_OFFSET1, 1);
			TecUtilArgListAppendDouble(argList, SV_DVALUE, IsoValue);
			TecUtilStyleSetLowLevelX(argList);
			TecUtilArgListClear(argList);
			TecUtilArgListAppendString(argList, SV_P1, SV_ISOSURFACEATTRIBUTES);
			TecUtilArgListAppendString(argList, SV_P2, SV_ISOVALUE3);
			TecUtilArgListAppendInt(argList, SV_OFFSET1, 1);
			TecUtilArgListAppendDouble(argList, SV_DVALUE, IsoValue);
			TecUtilStyleSetLowLevelX(argList);

			TecUtilArgListDealloc(&argList);

		}
	}


	if (HasAtomTypes)
		AuxDataDataSetSetItem(DLProgramName, "VASP 5");
	else
		AuxDataDataSetSetItem(DLProgramName, "VASP 4 (or earlier)");

	for (int i = 0; i < 3; ++i){
		stringstream ss;
		ss << "{ ";
		for (int j = 0; j < 3; ++j){
			ss << LatticeVector.at(i, j);
			if (j < 2) ss << ", ";
		}
		ss << " }";
		AuxDataDataSetSetItem(DLLatticeVecs[i], ss.str());
	}
	AuxDataDataSetSetItem(DLLatticeConstants, to_string(LatticeConstant));


	/*
	 *	Calculate the gradient and gradient 
	 *	magnitude of charge density.
	 */
// 	CalcGradGradMagForDataset(TRUE, AddOnID);

	if (IsOk && TotNumPoints > MaxPointsForShowVolumeZone){
		Set_pa TempSet = TecUtilSetAlloc(FALSE);
		IsOk = TecUtilSetAddMember(TempSet, VolZoneNum, FALSE);
		if (IsOk){
			TecUtilZoneSetActive(TempSet, AssignOp_MinusEquals);
		}
		TecUtilSetDealloc(&TempSet);
	}

	TecUtilRedrawAll(TRUE);

	TecUtilLockFinish(AddOnID);

	return IsOk;
}

int LoadADFTape41ASCIIData(char* FileNameCStr)
{
	TecUtilLockStart(AddOnID);
	/*
	 * Declaring variables
	 */

//	For file read
	ifstream  InFile;

//	Max number of characters on a single line
	int const LineLen = 500;

//	Strings repeatedly used
	string TempVarType, TempVarName;

//	Names of the input t41 and output .dat files
	string FileNameStr;

//	List of found variables; all possible variables
	vector<T41Var_s> T41VarList, TypeList;

// 	List of atom groups
	vector<AtomGroup_s> AtomGroupList;

//	Total number of found variables; all possible variables
	int NumVars, NumBaseTypes;

//	System IJK and total extent; number of lines per block
	lg_int_t IJK[3] = {0,0,0}, NumPoints, NumLines, RemPoints;

	int BufSize = 1000000;

	/*
	 * Populating list of all possible variables
	 */
	NumBaseTypes = PopulateTypeList(TypeList);

	/*
	 * Parse the input and output filenames from the argument list
	 */

	/*
	 * Now parse the output file argument to define
	 * the working directory, where temporary files
	 * will be stored.
	 */

	FileNameStr = FileNameCStr;
	string DataSetName;


#if defined(WIN32) || defined(_WIN32) || defined(__WIN32) && !defined(__CYGWIN__)
	static string const Slash="\\";
#else
	static string const Slash="/";
#endif

	size_t StrPos = FileNameStr.find_last_of("/\\");
	if (StrPos != string::npos){
		DataSetName = FileNameStr.substr(StrPos + 1, FileNameStr.length());
	}
	StrPos = DataSetName.find_first_of(".");
	if (StrPos != string::npos){
		DataSetName = DataSetName.substr(0, StrPos);
	}


	/*
	 * Open the input file for access
	 */
	InFile.open(FileNameStr.c_str(), ios::in);
	if (!InFile.is_open()){
		TecUtilDialogErrMsg("Failed to open selected ASCII Tape41 file");
	}

	/*
	 * Parse the beginning of the Tape41 txt file to
	 * find the IJK extent of the system
	 */
	TempVarType.resize(LineLen);
	for (int i = 0 ; i < 3 && !InFile.eof() ; ++i){
		while (!InFile.eof()){
			getline(InFile, TempVarType, '\n');
			if (!TempVarType.compare(0, 13, "nr of points "))
				break;
		}
		InFile.ignore(LineLen,'\n');
		InFile >> IJK[i];
	}

//	Calculate and check total number of points
	NumPoints = IJK[0] * IJK[1] * IJK[2];
	if (NumPoints <= 0 || IJK[0] <= 0 || IJK[1] <= 0 || IJK[2] <= 0){
		TecUtilDialogErrMsg("Failed to parse input file or I, J, or K have 0 points!");
		TecUtilLockFinish(AddOnID);
		return 1;
	}

//	Find whether or not there are dangling lines with less than 3 points
	NumLines = NumPoints / 3;
	RemPoints = NumPoints % 3;
	if (RemPoints > 0) ++NumLines;

	NumVars = 0;

	/*
	 * First main loop:
	 * 		Goes line by line through the Tape41 txt file until it finds one
	 * 		containing a variable name (i.e. "Density"), then it knows how
	 * 		many lines to read from the total number of points.
	 * 		This repeats until either the end of the Tape41 file is reached
	 * 		or until the total number of possible variables has been
	 * 		found.
	 *
	 * 		When data for the density gradient is found, it is saved in a
	 * 		temporary binary file for later use.
	 *
	 * 		Also keeps track of which variable appears and when.
	 *
	 * 		The current method does things one line at a time, and is therefore
	 * 		a bit slow. Could speed things up by reading in a set number of lines,
	 * 		processing that buffer, then repeating until the variable is
	 * 		finished.
	 */

	string IsoSurfaceVarName;
	EntIndex_t IsoSurfaceVarNum = 0;
	int CharSize = sizeof(char);
	unsigned int AtomCount = 0;

	StatusLaunch("Indexing variables", AddOnID, FALSE);

	while (!InFile.eof()){
		getline(InFile, TempVarType);
		int TypeNum = IsVarLine(TempVarType);
		if (TypeNum > 0){
			getline(InFile, TempVarName);// , '\n');
			for (int i = 0 ; i < NumBaseTypes ; ++i)
				if (!TempVarName.compare(0, TypeList[i].T41NameStr.length(), TypeList[i].T41NameStr)){
					T41Var_s TempType = TypeList[i];
					if (TypeNum > 3){
						TempType.SetAB(TempVarName);
						TempType.SetType(TempVarType);
					}
					if (!IsInDataTypeList(T41VarList, TempType)){
						if (TempType.VarType != GEOMETRY){

							TempType.ByteOffset = InFile.tellg();

							InFile.ignore(LineLen, '\n');
							InFile >> TempType.FirstValue;

							T41VarList.push_back(TempType);
							++NumVars;

							stringstream ss;
							ss << "Indexing variable " << NumVars << ": " << TempType.NameStr;
							if (!StatusUpdate(0, 1, ss.str(), AddOnID)){
								TecUtilLockFinish(AddOnID);
								return 0;
							}

							if (TempType.UseForIsosurface)
								IsoSurfaceVarName = TempType.NameStr;

							/*
							 *	Skip to the next variable block
							 **/
							//InFile.seekg(TempType.ByteOffset + ((NumLines) * 159 * CharSize));

							for (lg_int_t j = 0; j < NumLines && !InFile.eof(); ++j)
								InFile.ignore(LineLen, '\n');

							//TODO: Need to figure out this seekg problem, would be way faster than ignoring
							//	a bunch of lines.
						}
						else{
							/*
							 *	Now starts the hideous logic of getting the atomic count/type/position/charge information
							 */
							if (AtomCount == 0){
// 								If it's the nuclear count, which we know by the atom count still being 0
								InFile.ignore(LineLen, '\n');
								InFile >> AtomCount;
							}
							else if (!TempVarName.compare(0, 6, string("labels"))){
// 								If it's the atom labels (names)
								InFile.ignore(LineLen, '\n');
								string OldLine, NewLine;
								for (unsigned int j = 0; j < AtomCount && !InFile.eof(); ++j){
									getline(InFile, NewLine);
									InFile.ignore(LineLen, '\n');

									NewLine = NewLine.substr(0, NewLine.find_first_of(" "));
									if (NewLine != OldLine){
										AtomGroupList.push_back(AtomGroup_s(NewLine));
										OldLine = NewLine;
									}
									else{
										size_t AtomNum = AtomGroupList.size();
										if (AtomNum > 0)
											AtomGroupList[AtomNum - 1].Count++;
										else{
											TecUtilDialogErrMsg("Failed Parsing Atom List: Atom Type Names");
											TecUtilLockFinish(AddOnID);
											return 1;
										}
									}
								}
							}
							else if (!TempVarName.compare(0, 6, string("xyznuc"))){
// 								If it's the atom positions
								if (AtomCount == 0){
									TecUtilDialogErrMsg("Failed Parsing Atom List: Count wasn't found");
									TecUtilLockFinish(AddOnID);
									return 1;
								}
								else{
									InFile.ignore(LineLen, '\n');
									for (int GroupNum = 0; GroupNum < AtomGroupList.size() && !InFile.eof(); ++GroupNum){
										for (int AtomNum = 0; AtomNum < AtomGroupList[GroupNum].Count && !InFile.eof(); ++AtomNum){
											vec3 Position;
											for (int j = 0; j < 3 && !InFile.eof(); ++j)
												InFile >> Position[j];
											AtomGroupList[GroupNum].AddPosition(Position);
										}
									}
								}
							}
							else if (!TempVarName.compare(0, 4, string("qtch"))){
// 								If it's the atomic charges
								if (AtomCount == 0){
									TecUtilDialogErrMsg("Failed Parsing Atom List: Count wasn't found");
									TecUtilLockFinish(AddOnID);
									return 1;
								}
								else{
									InFile.ignore(LineLen, '\n');
									for (int GroupNum = 0; GroupNum < AtomGroupList.size() && !InFile.eof(); ++GroupNum){
										for (int AtomNum = 0; AtomNum < AtomGroupList[GroupNum].Count && !InFile.eof(); ++AtomNum){
											double Charge(3);
											InFile >> Charge;
											AtomGroupList[GroupNum].Charges.push_back(Charge);
										}
									}
								}
							}
						}

						break;
					}
				}
		}
	}

	/*
	 *	Now we know all the variables that are in the file, so
	 *	create the data set, all the variables, and then use
	 *	the saved byte offsets to zip through the file and 
	 *	read all the variable blocks into their respective
	 *	tecplot variables.
	 **/
// TODO:	let the user pick which variables they want to load using a separate dialog, then only read in the ones they select.

//	Close Tape41 file
// 	InFile.close();
// 	
// 	Clear the flags of the input file
	InFile.clear();

	TecUtilDataLoadBegin();

	/*
	 *	Make the dataset using the first variable
	 **/
	StringList_pa VarNames = TecUtilStringListAlloc();
	for (int i = 0; i < T41VarList.size(); ++i)
		TecUtilStringListAppendString(VarNames, T41VarList[i].NameStr.c_str());

	vector<FieldDataType_e> VarDataTypes;
	VarDataTypes.resize(T41VarList.size(), FieldDataType_Float);

	if (TecUtilDataSetCreate(StringMakeValidZoneName(DataSetName).c_str(), VarNames, TRUE)){
		TecUtilDataSetAddZone(StringMakeValidZoneName(CSMZoneName.FullVolume + DataSetName).c_str(), IJK[0], IJK[1], IJK[2], ZoneType_Ordered, VarDataTypes.data());
	}
	else{
		TecUtilDialogErrMsg("Failed to create new data set");
		TecUtilLockFinish(AddOnID);
		return 1;
	}
	TecUtilStringListDealloc(&VarNames);
	EntIndex_t ZoneNum, VolZoneNum;
	VolZoneNum = TecUtilDataSetGetNumZones();

	SetZoneStyle({ VolZoneNum }, ZoneStyle_Volume);

	/*
	 *	Now loop over list of variables in the file, reading them in some set number of points at a time
	 *	and giving them to tecplot as chuncks of data
	 **/
	LgIndex_t TotalPoints = NumPoints * static_cast<int>(T41VarList.size());
	LgIndex_t CurrentOverallPoint = 0;

	StatusUpdate(0, 1, "Loading ASCII Tape41 File", AddOnID);

	int BufferMemoryKB = sizeof(ImportType_t) * BufSize / 1024;

	for (int i = 0; i < T41VarList.size(); ++i){
// 		vector<double> DoubleBuffer;
		vector<ImportType_t> DoubleBuffer;
		TecUtilMemoryChangeNotify(BufferMemoryKB);
		DoubleBuffer.reserve(BufSize);
		string TempStr;
		ImportType_t TempDouble = NAN;
		LgIndex_t PointNum = 1;

		FieldData_pa VarRef = TecUtilDataValueGetWritableNativeRef(VolZoneNum, i + 1);
		if (VALID_REF(VarRef)){
			InFile.seekg(T41VarList[i].ByteOffset - 160);
			InFile.ignore(LineLen, '\n');

			while (TempDouble != T41VarList[i].FirstValue){
				getline(InFile, TempStr);
				int LinePoints = 3;
				if (NumPoints <= 3 && RemPoints > 0)
					LinePoints = RemPoints;
				stringstream ss;
				ss << TempStr;
				ss >> TempDouble;
				if (TempDouble == T41VarList[i].FirstValue){
					for (int k = 0; k < LinePoints; ++k){
						if (k > 0)
							ss >> TempDouble;
						DoubleBuffer.push_back(TempDouble);
						CurrentOverallPoint++;
					}
					break;
				}
			}
		}



		for (int LineNum = 1; LineNum < NumLines; ++LineNum){
			getline(InFile, TempStr);
			int LinePoints = 3;
			if (RemPoints > 0 && LineNum >= NumLines - 1)
				LinePoints = RemPoints;
			stringstream ss;
			ss << TempStr;
			for (int k = 0; k < LinePoints; ++k){
				ss >> TempDouble;
				DoubleBuffer.push_back(TempDouble);
				CurrentOverallPoint++;
				if (DoubleBuffer.size() >= BufSize){
					TecUtilDataValueArraySetByRef(VarRef, PointNum, static_cast<int>(DoubleBuffer.size()), DoubleBuffer.data());
					PointNum += static_cast<int>(DoubleBuffer.size());
					DoubleBuffer.clear();
				}
			}

			if (!StatusUpdate(CurrentOverallPoint-1, TotalPoints, "", AddOnID)){
				StatusDrop(AddOnID);
				TecUtilLockFinish(AddOnID);
				return 0;
			}
		}

		if (DoubleBuffer.size() > 0){
			TecUtilDataValueArraySetByRef(VarRef, PointNum, static_cast<int>(DoubleBuffer.size()), DoubleBuffer.data());
			PointNum += static_cast<int>(DoubleBuffer.size());
			DoubleBuffer.clear();
		}
		InFile.clear();

		TecUtilMemoryChangeNotify(-BufferMemoryKB);
	}

	InFile.close();

	StatusDrop(AddOnID);

	/*
	 *	Keep list of created zones for informing Tecplot later
	 */
	Set_pa NewZones = TecUtilSetAlloc(TRUE);
	TecUtilSetAddMember(NewZones, VolZoneNum, TRUE);
	

	/*
	 *	Now create zones for each atom group and set their style settings accordingly
	 */
	vector<AtomColor_s> AtomColorList;
	PopulateAtomColorList(AtomColorList);
	for (int GroupNum = 0; GroupNum < AtomGroupList.size(); ++GroupNum){
		TecUtilDataSetAddZone(StringMakeValidZoneName(AtomGroupList[GroupNum].Name).c_str(), AtomGroupList[GroupNum].Count, 1, 1, ZoneType_Ordered, VarDataTypes.data());

		ZoneNum = TecUtilDataSetGetNumZones();
		FieldData_pa VarRef[3] = {
			TecUtilDataValueGetWritableNativeRef(ZoneNum, TecUtilVarGetNumByName("X")),
			TecUtilDataValueGetWritableNativeRef(ZoneNum, TecUtilVarGetNumByName("Y")),
			TecUtilDataValueGetWritableNativeRef(ZoneNum, TecUtilVarGetNumByName("Z"))
		};
		if (VALID_REF(VarRef[0]) && VALID_REF(VarRef[0]) && VALID_REF(VarRef[0])){
			for (int i = 0; i < 3; ++i){
				TecUtilDataValueArraySetByRef(VarRef[i], 1, AtomGroupList[GroupNum].Count, AtomGroupList[GroupNum].Positions[i].data());
			}
		}

		TecUtilSetAddMember(NewZones, ZoneNum, TRUE);

		Set_pa TempSet = TecUtilSetAlloc(TRUE);
		TecUtilSetAddMember(TempSet, ZoneNum, TRUE);

		TecUtilZoneSetScatter(SV_SHOW, TempSet, 0.0, TRUE);
		TecUtilZoneSetScatter(SV_COLOR, TempSet, 0.0, GetAtomColor(AtomColorList, AtomGroupList[GroupNum].Name));
		TecUtilZoneSetScatter(SV_FRAMESIZE, TempSet, 2*sqrt(AtomGroupList[GroupNum].Charges[0]), 0);
		TecUtilZoneSetScatterSymbolShape(SV_GEOMSHAPE, TempSet, GeomShape_Sphere);

		TecUtilZoneSetMesh(SV_SHOW, TempSet, 0.0, FALSE);
		TecUtilZoneSetContour(SV_SHOW, TempSet, 0.0, FALSE);
		TecUtilZoneSetVector(SV_SHOW, TempSet, 0.0, FALSE);
		TecUtilZoneSetShade(SV_SHOW, TempSet, 0.0, FALSE);
		TecUtilZoneSetEdgeLayer(SV_SHOW, TempSet, 0.0, FALSE);
		

		TecUtilSetDealloc(&TempSet);


	}

	/*
	 *	Inform Tecplot of the new zones
	 */
	TecUtilStateChanged(StateChange_ZonesAdded, (ArbParam_t)NewZones);

	TecUtilFrameSetPlotType(PlotType_Cartesian3D);

	TecUtilZoneSetActive(NewZones, AssignOp_PlusEquals);

	TecUtilSetDealloc(&NewZones);

	ArgList_pa argList = TecUtilArgListAlloc();

	TecUtilArgListClear(argList);
	TecUtilArgListAppendString(argList, SV_P1, SV_ISOSURFACEATTRIBUTES);
	TecUtilArgListAppendString(argList, SV_P2, SV_EFFECTS);
	TecUtilArgListAppendString(argList, SV_P3, SV_USETRANSLUCENCY);
	TecUtilArgListAppendInt(argList, SV_OFFSET1, 1);
	TecUtilArgListAppendArbParam(argList, SV_IVALUE, TRUE);
	TecUtilStyleSetLowLevelX(argList);

	TecUtilArgListClear(argList);
	TecUtilArgListAppendString(argList, SV_P1, SV_ISOSURFACEATTRIBUTES);
	TecUtilArgListAppendString(argList, SV_P2, SV_EFFECTS);
	TecUtilArgListAppendString(argList, SV_P3, SV_SURFACETRANSLUCENCY);
	TecUtilArgListAppendInt(argList, SV_OFFSET1, 1);
	TecUtilArgListAppendArbParam(argList, SV_IVALUE, 70);
	TecUtilStyleSetLowLevelX(argList);

	TecUtilArgListClear(argList);
	TecUtilArgListAppendString(argList, SV_P1, SV_FIELDLAYERS);
	TecUtilArgListAppendString(argList, SV_P2, SV_SHOWSCATTER);
	TecUtilArgListAppendArbParam(argList, SV_IVALUE, TRUE);
	TecUtilStyleSetLowLevelX(argList);

	TecUtilArgListClear(argList);
	TecUtilArgListAppendString(argList, SV_P1, SV_ISOSURFACELAYERS);
	TecUtilArgListAppendString(argList, SV_P2, SV_SHOW);
	TecUtilArgListAppendArbParam(argList, SV_IVALUE, TRUE);
	TecUtilStyleSetLowLevelX(argList);

	TecUtilArgListClear(argList);
	TecUtilArgListAppendString(argList, SV_P1, SV_ISOSURFACEATTRIBUTES);
	TecUtilArgListAppendString(argList, SV_P2, SV_DEFINITIONCONTOURGROUP);
	TecUtilArgListAppendInt(argList, SV_OFFSET1, 1);
	TecUtilArgListAppendArbParam(argList, SV_IVALUE, (SmInteger_t)1);
	TecUtilStyleSetLowLevelX(argList);

	TecUtilArgListClear(argList);
	TecUtilArgListAppendString(argList, SV_P1, SV_ISOSURFACEATTRIBUTES);
	TecUtilArgListAppendString(argList, SV_P2, SV_SHOWGROUP);
	TecUtilArgListAppendInt(argList, SV_OFFSET1, 1);
	TecUtilArgListAppendArbParam(argList, SV_IVALUE, TRUE);
	TecUtilStyleSetLowLevelX(argList);

	IsoSurfaceVarNum = TecUtilVarGetNumByName(IsoSurfaceVarName.c_str());
// 	double ValMax, ValMin;
// 	TecUtilVarGetMinMax(IsoSurfaceVarNum, &ValMin, &ValMax);
// 	double IsoValue = ValMax - 0.05 * (ValMax - ValMin);
	double IsoValue = 0.25;

	TecUtilArgListClear(argList);
	TecUtilArgListAppendString(argList, SV_P1, SV_ISOSURFACEATTRIBUTES);
	TecUtilArgListAppendString(argList, SV_P2, SV_ISOVALUE1);
	TecUtilArgListAppendInt(argList, SV_OFFSET1, 1);
	TecUtilArgListAppendDouble(argList, SV_DVALUE, IsoValue);
	TecUtilStyleSetLowLevelX(argList);


	TecUtilArgListClear(argList);
	TecUtilArgListAppendString(argList, SV_P1, SV_GLOBALCONTOUR);
	TecUtilArgListAppendString(argList, SV_P2, SV_VAR);
	TecUtilArgListAppendInt(argList, SV_OFFSET1, 1);
	TecUtilArgListAppendArbParam(argList, SV_IVALUE, IsoSurfaceVarNum);
	TecUtilStyleSetLowLevelX(argList);
	TecUtilArgListClear(argList);
	TecUtilArgListAppendString(argList, SV_P1, SV_ISOSURFACEATTRIBUTES);
	TecUtilArgListAppendString(argList, SV_P2, SV_ISOVALUE1);
	TecUtilArgListAppendInt(argList, SV_OFFSET1, 1);
	TecUtilArgListAppendDouble(argList, SV_DVALUE, IsoValue);
	TecUtilStyleSetLowLevelX(argList);
	TecUtilArgListClear(argList);
	TecUtilArgListAppendString(argList, SV_P1, SV_ISOSURFACEATTRIBUTES);
	TecUtilArgListAppendString(argList, SV_P2, SV_ISOVALUE2);
	TecUtilArgListAppendInt(argList, SV_OFFSET1, 1);
	TecUtilArgListAppendDouble(argList, SV_DVALUE, IsoValue);
	TecUtilStyleSetLowLevelX(argList);
	TecUtilArgListClear(argList);
	TecUtilArgListAppendString(argList, SV_P1, SV_ISOSURFACEATTRIBUTES);
	TecUtilArgListAppendString(argList, SV_P2, SV_ISOVALUE3);
	TecUtilArgListAppendInt(argList, SV_OFFSET1, 1);
	TecUtilArgListAppendDouble(argList, SV_DVALUE, IsoValue);
	TecUtilStyleSetLowLevelX(argList);

	TecUtilArgListDealloc(&argList);

	TecUtilDataLoadEnd();

	TecUtilViewDataFit();


	/*
	 *	Now we should have all variables read in to the new zone in the new
	 *	dataset. This leaves things the same as they were after using the
	 *	CLI converter program I made. 
	 *	
	 *	Extra variables will be calculated inside Bondalyzer or by specific
	 *	utilities, not by this loader.
	 */

	AuxDataDataSetSetItem(DLProgramName, "SCM ADF");

	CalcGradGradMagForDataset(FALSE, AddOnID);

	if (NumPoints > MaxPointsForShowVolumeZone){
		Set_pa TempSet = TecUtilSetAlloc(FALSE);
		Boolean_t IsOk = TecUtilSetAddMember(TempSet, VolZoneNum, FALSE);
		if (IsOk){
			TecUtilZoneSetActive(TempSet, AssignOp_MinusEquals);
		}
		TecUtilSetDealloc(&TempSet);
	}

	TecUtilLockFinish(AddOnID);
	
	return 0;
}

LoaderStatus_e GetT41Contents(Boolean_t IsADFFile, Boolean_t SaveToList){
	/*
	 *	Create lists of strings that are combined to form all the
	 *	possible combinations of section and variable names
	 */

	Boolean_t IsOk = TRUE;

	vector<string> T41Sec, T41Var, VarStrings, T41VarSuffices, MOTypeStrings;

	if (IsADFFile){
		T41Sec = {
			"SCF",
			"SumFrag",
			"Ortho",
			"Core"
		};
		T41Var = {
			"Density",
			"Fitdensity",
			"Kinetic Energy Density",
			"CoulPot",
			"XCPot",
			"ELF",
			"DensityGradX",
			"DensityGradY",
			"DensityGradZ",
			"DensityHessXX",
			"DensityHessXY",
			"DensityHessXZ",
			"DensityHessYY",
			"DensityHessYZ",
			"DensityHessZZ",
			"DensityLap"
		};
		/*
		*	And the strings that will be used to name the variables
		*/
		VarStrings = {
			"Electron Density",
			"Electron Density Fit",
			"Kinetic Energy Density",
			"Coulomb Potential",
			"Exchange Correlation Potential",
			"Electron Localization Function",
			"X Density Gradient",
			"Y Density Gradient",
			"Z Density Gradient",
			"XX Density Hessian",
			"XY Density Hessian",
			"XZ Density Hessian",
			"YY Density Hessian",
			"YZ Density Hessian",
			"ZZ Density Hessian",
			"Density Laplacian"
		};
		T41VarSuffices = { "_A", "_B", "" };
		MOTypeStrings = { "SCF", "LOC" };
	}
	else{
		T41Sec = {
			"DENSITY",
			"FITDENSITY",
			"COULOMBPOTENTIAL",
			"XCPOTENTIAL",
			"TAU",
			"LDOS",
			"ELF",
			"ATOMIC",
			"FOO"
		};
		T41Var = {
			"total_scf",
			"scf",
			"scf_exact",
			"agrad_total_scf",
			"deformation_scf",
			"density",
			"coulombPot",
			"rho"
		};
		for (int i = 1; i <= 3; ++i) T41Var.push_back("gradRho_" + to_string(i));
		for (int i = 1; i <= 6; ++i) T41Var.push_back("secDerRho_" + to_string(i));

		VarStrings = {
			"Electron Density",
			"Electron Density Fit",
			"Exact Electron Density",
			"Norm of the Gradient of Electron Density",
			"Fitted Deformation Density",
			"Density of Startup Atoms",
			"Coulomb Potential of Startup Atom Density",
			"Electron Density"
		};
		for (string const & s : { "X", "Y", "Z" }) VarStrings.push_back(s + " Density Gradient");
		for (string const & s : { "XX", "XY", "XZ", "YY", "YZ", "ZZ" }) VarStrings.push_back(s + " Density Hessian");
		T41VarSuffices = { "" };
	}

	/*
	 *	Transition density is a little different
	 */
	string T41TransDensSec = "TransDens";
	string TransDensString = "Transition";
	vector<string> T41TransDensSecSuf = { "", "_L1", "_L2" };
	vector<string> T41TransDensVar = { "Fitdensity", "Coulpot" };
	vector<string> TransDensVarStrs = { "Electron Density Fit", "Coulomb Potential" };
	vector<string> T41TransDensVarSuf = { "", "_1", "_2", "_3" };

	/*
	 *	Open Tape41 kf file
	 */

	IsOk = (openKFFile(&GlobalKFFile, const_cast<char*>(FileNameStr.c_str())) > 0);

	if (!IsOk){
		TecUtilDialogErrMsg("Failed to open Tape41 file.");
		return LoaderFailed;
	}


	/*
	*	Get all information about contents of Tape21.
	*	Includes names of sections/variables and length/type of variables
	*/
	vector<string> SectionNames = GetSectionNameList(&GlobalKFFile);

	vector<vector<string> > VariableNames(SectionNames.size());
	vector<vector<int> > TypeList(SectionNames.size()), LengthList(SectionNames.size());

	for (int i = 0; i < SectionNames.size() && IsOk; ++i)
		IsOk = GetVariableInfoForSection(&GlobalKFFile, SectionNames[i].c_str(), VariableNames[i], TypeList[i], LengthList[i]);

	for (string & i : SectionNames)
		i = i.substr(0, i.find_last_not_of(' ') + 1);

	/*
	 *	Get xyz variable info
	 */

	int TotNumPts = -1;

	vector<string> XYZVarStrings = { "x values%x values", "y values%y values", "z values%z values" };
	if (IsOk){
		if (IsADFFile){
			for (int i = 0; i < 3 && IsOk; ++i){
				TotNumPts = getKFVariableLength(&GlobalKFFile, XYZVarStrings[i].c_str());
				IsOk = (TotNumPts > 0);
			}
			if (!IsOk)
				TecUtilDialogErrMsg("Cannot load file without X, Y, and Z grid information. Include the \"save\" keyword in the grid block when calling densf");
		}
		else{
			if (getKFData(&GlobalKFFile, "Grid%total nr of points", &TotNumPts) <= 0){
				IsOk = FALSE;
				TecUtilDialogErrMsg("Failed to read grid information from Tape41");
			}
		}
	}

	/*
	 *	Get list of MOs so they can be put after the main variables
	 */

	vector<string> MOUserVarStrings, MOKFVarStrings;
	vector<double> MOEnergies;
	vector<int> MODataSizes;
	if (IsOk && IsADFFile){
		/*
		*	Then for the MOs.
		*	This one's different. First need to check the symmetries to get
		*	the names of the orbital sections. Then loop over those.
		*/
		int MOLabelSize = getKFVariableLength(&GlobalKFFile, "Grid%labels");
		if (MOLabelSize > 0){
			vector<char> CharBuf(MOLabelSize);
			// 			char *CharBuf = new char[MOLabelSize];
			if (getKFData(&GlobalKFFile, "Grid%labels", reinterpret_cast<void*>(CharBuf.data())) > 0){
				stringstream ss;
				ss << CharBuf.data();
				// 				delete CharBuf;
				vector<string> MOLabels(MOLabelSize / 160);
				for (string & i : MOLabels)
					ss >> i;

				/*
				*	Now we have all the symmetry labels, so we can go through the file
				*	and find out which MOs are present.
				*	To help the user, we'll also save the energy and occupancy of the
				*	orbitals and present them with a sorted list in the order they'd
				*	expect.
				*	Do this for SCM and then for Localized orbitals.
				*/
				for (string const & Type : MOTypeStrings){
					for (string const & Sym : MOLabels){
						for (string const & Spin : T41VarSuffices) {
							/*
							*	Query the MO section and see how many orbitals there are.
							*/
							string MOStr = Type + "_" + Sym + Spin;
							int NumMOSize = getKFVariableLength(&GlobalKFFile, string(MOStr + "%nr of orbitals").c_str());
							if (NumMOSize > 0) {
								/*
								*	Get the MO occupancies and energies
								*/
								int NumMOs;
								if (getKFData(&GlobalKFFile, string(MOStr + "%nr of orbitals").c_str(), reinterpret_cast<void*>(&NumMOs)) > 0) {
									vector<double> Occ(NumMOs), EigVals(NumMOs);
									if (getKFData(&GlobalKFFile, string(MOStr + "%Occupations").c_str(), Occ.data()) > 0
										&& getKFData(&GlobalKFFile, string(MOStr + "%Eigenvalues").c_str(), EigVals.data()) > 0)
									{
										/*
										*	Now actually see if each MO is in the file.
										*	If it is, construct it's string that includes
										*	energy and occupancy and include it in the temporary
										*	MO variable list.
										*/
										for (int i = 0; i < NumMOs; ++i) {
											int TmpSize = getKFVariableLength(&GlobalKFFile, string(MOStr + "%" + to_string(i + 1)).c_str());
											if (TmpSize > 0) {
												/*
												*	MO exists, so add it to temporary lists.
												*/
												stringstream ESS, OccSS;
												string EStr, OccStr;

												ESS.precision(3);
												OccSS.precision(3);

												if (abs(EigVals[i]) >= 1e-3)
													ESS << EigVals[i];
												else
													ESS << std::scientific << EigVals[i];

												if (Occ[i] == static_cast<int>(Occ[i]))
													OccSS << static_cast<int>(Occ[i]);
												else if (Occ[i] >= 1e-3)
													OccSS << Occ[i];
												else
													OccSS << std::scientific << Occ[i];

												ESS >> EStr;
												OccSS >> OccStr;

												MOUserVarStrings.push_back(string("MO E: " + EStr + ", Occ: " + OccStr + ", " + MOStr + "_" + to_string(i + 1)));
												MOKFVarStrings.push_back(string(MOStr + "%" + to_string(i + 1)));
												MOEnergies.push_back(EigVals[i]);
												MODataSizes.push_back(TmpSize);
											}
										}
									}
								}
							}
						}
					}

					/*
					*	Now we have an unsorted list of MO.
					*	Sort the list by energies.
					*	Go bubble sort!
					*/
					if (MOUserVarStrings.size() > 1){
						Boolean_t Swapped = TRUE;
						while (Swapped){
							Swapped = FALSE;
							for (int i = 0; i < MOUserVarStrings.size() - 1; ++i){
								if (MOEnergies[i + 1]<MOEnergies[i]){
									Swapped = TRUE;

									double TmpDbl = MOEnergies[i];
									MOEnergies[i] = MOEnergies[i + 1];
									MOEnergies[i + 1] = TmpDbl;

									int TmpInt = MODataSizes[i];
									MODataSizes[i] = MODataSizes[i + 1];
									MODataSizes[i + 1] = TmpInt;

									string TmpStr = MOUserVarStrings[i];
									MOUserVarStrings[i] = MOUserVarStrings[i + 1];
									MOUserVarStrings[i + 1] = TmpStr;

									TmpStr = MOKFVarStrings[i];
									MOKFVarStrings[i] = MOKFVarStrings[i + 1];
									MOKFVarStrings[i + 1] = TmpStr;
								}
							}
						}
					}

					/*
					 *	MO list has been sorted. It will be added to
					 *	the list for user to select from later.
					 */
				}
			}
		}
	}
	closeKFFile(&GlobalKFFile);

	/*
	 *	Go though the full section%variable contents and try to match expected variables
	 *	from the lists defined above, so that they can be represented using more
	 *	"friendly" names.
	 */

	if (IsOk) for (int SecNum = 0; SecNum < SectionNames.size(); ++SecNum){
		for (int VarNum = 0; VarNum < VariableNames[SecNum].size(); ++VarNum){
			if (TypeList[SecNum][VarNum] = KF_T_DOUBLE && LengthList[SecNum][VarNum] == TotNumPts){
				/*
				 *	Only importing double precision data that has values on the whole grid.
				 */
				string ChkStr = SectionNames[SecNum] + "%" + VariableNames[SecNum][VarNum];
				Boolean_t VarFound = FALSE;
				for (string const & i : XYZVarStrings)
					if (ChkStr == i){
						VarFound = TRUE;
						break;
					}

				if (!VarFound){
					VarFound = (SearchVectorForString(T41LoadKFVarStrings, ChkStr, false) >= 0);
				}
				/*
				 *	First check against "main" variables
				 */
				if (!VarFound) for (string const & Sec : T41Sec){
					// Check if SectionNames[SecNum] starts with Sec
					if (SectionNames[SecNum].substr(0, Sec.length()) == Sec) {
						// Get the suffix (if any) from SectionNames[SecNum]
						string SectionSuffix = "";
						if (SectionNames[SecNum].length() > Sec.length()) {
							SectionSuffix = SectionNames[SecNum].substr(Sec.length());
						}

						unsigned int VarChkNum = 0;
						for (string const & Var : T41Var){
							for (string const & Suffix : T41VarSuffices){
								string Str1 = SectionNames[SecNum] + "%" + Var + Suffix;  // Use actual section name
								if (Str1 == ChkStr){
									string Str2 = VarStrings[VarChkNum] + Suffix + " " + Sec + SectionSuffix;  // Add section suffix
									T41LoadKFVarStrings.push_back(Str1);
									T41LoadUserVarStrings.push_back(Str2);
									T41DataSizes.push_back(LengthList[SecNum][VarNum]);
									VarFound = TRUE;
								}
								if (VarFound) break;
							}
							++VarChkNum;
							if (VarFound) break;
						}
						if (VarFound) break;
					}
				}

				/*
				 *	Then check against transition variables.
				 */
				if (!VarFound && IsADFFile){
					unsigned int VarNum = 0;
					for (string const & Var : T41TransDensVar){
						for (string const & Suffix : T41TransDensVarSuf){
							string Str1 = T41TransDensSec + "%" + Var + Suffix;
							if (Str1 == ChkStr){
								string Str2 = TransDensString + " " + TransDensVarStrs[VarNum] + Suffix;
								T41LoadKFVarStrings.push_back(Str1);
								T41LoadUserVarStrings.push_back(Str2);
								T41DataSizes.push_back(LengthList[SecNum][VarNum]);
								VarFound = TRUE;
							}
							if (VarFound) break;
						}
						++VarNum;
						if (VarFound) break;
					}
				}

				/*
				 *	See if it's an orbital, in which case we skip it until after all the other vars have been added.
				 */
				if (!VarFound && IsADFFile) for (string const & i : MOKFVarStrings){
					if (i == ChkStr){
						VarFound = TRUE;
						break;
					}
				}

				/*
				 *	If no variable matched then just use the variable's name
				 */
				if (!VarFound){
					T41LoadKFVarStrings.push_back(ChkStr);
					T41LoadUserVarStrings.push_back(VariableNames[SecNum][VarNum] + " " + SectionNames[SecNum]);
					T41DataSizes.push_back(LengthList[SecNum][VarNum]);
				}
			}
		}
	}

	if (IsOk && IsADFFile){
		/*
		*	Now append the sorted MO list to the variable list.
		*/

		for (int i = 0; i < MOUserVarStrings.size(); ++i){
			if (SearchVectorForString(T41LoadKFVarStrings, MOKFVarStrings[i], false) < 0){
				T41LoadUserVarStrings.push_back(MOUserVarStrings[i]);
				T41LoadKFVarStrings.push_back(MOKFVarStrings[i]);
				T41DataSizes.push_back(MODataSizes[i]);
			}
		}
	}

	if (IsOk) IsOk = T41LoadUserVarStrings.size() > 0;

	if (!IsOk){
		TecUtilDialogErrMsg("Failed to read Tape41 file.");
		QuitT41Load();
		return LoaderFailed;
	}
	else if (SaveToList){
		for (auto Str = T41LoadUserVarStrings.cbegin(), StrEnd = T41LoadUserVarStrings.cend(); Str != StrEnd; Str++)
			TecGUIListAppendItem(MLT41LoadVa_MLST_D1, Str->c_str());

		TecGUIListSetSelectedItem(MLT41LoadVa_MLST_D1, 1);
		vector<LgIndex_t> TmpInts(T41LoadUserVarStrings.size());
		for (int i = 1; i <= T41LoadUserVarStrings.size(); ++i)
			TmpInts[i - 1] = i;
		TecGUIListSetSelectedItems(MLT41LoadVa_MLST_D1, TmpInts.data(), T41LoadUserVarStrings.size());
		if (T41LoadUserVarStrings.size() == 1){
			if (IsADFFile)
				LoadADFTape41Files();
			else
				SelVarsLoadData();
			
			return Loaded;
		}
		return Loading;
	}

	return Loading;
}

/*
 *	Here's the process of getting function pointers for direct access to field data:
	 GTotVarFDPtr = TecUtilDataValueGetReadableNativeRef(ZoneNum, GradMagVarNum);
	 FieldValueGetFunction_pf GTotVarGetFunction = TecUtilDataValueRefGetGetFunc(GTotVarFDPtr);
	 GTXVarFDPtr = TecUtilDataValueGetWritableNativeRef(ZoneNum, SGradGTXVarNum);
	 FieldValueSetFunction_pf GTXVarSetFunction = TecUtilDataValueRefGetSetFunc(GTXVarFDPtr);
 */

void QuitT41Load(){
	T41DataSizes.clear();
	T41LoadKFVarStrings.clear();
	T41LoadUserVarStrings.clear();
	FileNameStr.clear();
	FileNameStrs.clear();
}

void LoadADFTape41Files(){
	GroupIndex = 50;
	for (string const & i : FileNameStrs){
		FileNameStr = i;
		LoadADFTape41Data();
		ReplaceDataSet = FALSE;
		GroupIndex++;
	}
	QuitT41Load();
}

void LoadADFTape41Data(){

	TecUtilLockStart(AddOnID);

	LgIndex_t *SelectedVarNumsArray = nullptr;
	LgIndex_t NumSelected = -1;
	TecGUIListGetSelectedItems(MLT41LoadVa_MLST_D1, &SelectedVarNumsArray, &NumSelected);

	if (NumSelected < 1){
		TecUtilDialogErrMsg("You must select at least one variable to load from the Tape41 file");
		TecUtilArrayDealloc(reinterpret_cast<void**>(&SelectedVarNumsArray));
		TecGUIDialogLaunch(Dialog1Manager);
		TecUtilLockFinish(AddOnID);
		return;
	}

	string DataSetName;
	vector<LgIndex_t> SelectedVarNums(SelectedVarNumsArray, SelectedVarNumsArray + NumSelected);
	TecUtilArrayDealloc(reinterpret_cast<void**>(&SelectedVarNumsArray));


#if defined(WIN32) || defined(_WIN32) || defined(__WIN32) && !defined(__CYGWIN__)
	string Slash = "\\";
#else
	string Slash = "/";
#endif

	size_t StrPos = FileNameStr.find_last_of("/\\");
	if (StrPos != string::npos){
		DataSetName = FileNameStr.substr(StrPos + 1, FileNameStr.length());
	}
	StrPos = DataSetName.find_last_of(".");
	if (StrPos != string::npos){
		DataSetName = DataSetName.substr(0, StrPos);
	}

	Boolean_t IsOk = TRUE;

	IsOk = (openKFFile(&GlobalKFFile, const_cast<char*>(FileNameStr.c_str())) > 0);

	int DataLength = -1;

	if (IsOk){
		IsOk = FALSE;
		for (int i = 0; i < T41LoadKFVarStrings.size(); ++i){
			T41DataSizes[i] = getKFVariableLength(&GlobalKFFile, T41LoadKFVarStrings[i].c_str());
			if (T41DataSizes[i] > 0){
				if (!IsOk) IsOk = TRUE;
				DataLength = T41DataSizes[i];
			}
			else
				T41DataSizes[i] = 0;
		}
		if (!IsOk){
			closeKFFile(&GlobalKFFile);
			QuitT41Load();
			return;
		}

	}

	int MaxIJK[3];
	int MaxI = -1, MaxJ = -1, MaxK = -1, NumPoints = -1;

	vector<string> TempStrings = { "Grid%nr of points x", "Grid%nr of points y", "Grid%nr of points z" };

	for (int i = 0; i < 3 && IsOk; ++i)
		IsOk = (getKFData(&GlobalKFFile, TempStrings[i].c_str(), &MaxIJK[i]) > 0);

	if (IsOk){
		MaxI = MaxIJK[0];
		MaxJ = MaxIJK[1];
		MaxK = MaxIJK[2];
		IsOk = (getKFData(&GlobalKFFile, "Grid%total nr of points", reinterpret_cast<void*>(&NumPoints)) > 0);
	}

	if (!IsOk || MaxI < 0 || MaxJ < 0 || MaxK < 0 || NumPoints != MaxI*MaxJ*MaxK){
		TecUtilDialogErrMsg("Failed to get system extent");
	}

	TempStrings = { "X", "Y", "Z" };
	for (int i = 2; i >= 0; --i){
		if (SearchVectorForString(T41LoadUserVarStrings, TempStrings[i], false) < 0){
			T41LoadUserVarStrings.insert(T41LoadUserVarStrings.begin(), TempStrings[i]);
			if (i == 2)
				T41DataSizes.insert(T41DataSizes.begin(), 3, DataLength);
		}
	}
// 	T41LoadUserVarStrings.insert(T41LoadUserVarStrings.begin(), TempStrings.cbegin(), TempStrings.cend());

	TempStrings = { "x values%x values", "y values%y values", "z values%z values" };
	for (int i = 2; i >= 0; --i)
		if (SearchVectorForString(T41LoadKFVarStrings, TempStrings[i], false) < 0)
			T41LoadKFVarStrings.insert(T41LoadKFVarStrings.begin(), TempStrings[i]);
// 	T41LoadKFVarStrings.insert(T41LoadKFVarStrings.begin(), TempStrings.cbegin(), TempStrings.cend());
	TempStrings.clear();

	vector<LgIndex_t> TempLgIndex = { 1, 2, 3 };
	for (int & i : SelectedVarNums)
		i += 3;
	SelectedVarNums.insert(SelectedVarNums.begin(), TempLgIndex.cbegin(), TempLgIndex.cend());
	TempLgIndex.clear();

	StringList_pa VarNames = TecUtilStringListAlloc();
	for (int i : SelectedVarNums)
		TecUtilStringListAppendString(VarNames, T41LoadUserVarStrings[i - 1].c_str());

	vector<FieldDataType_e> VarDataTypes;

	/*
	 *	Check that data to be imported is not greater than available memory.
	 *	If it is, have Tacplot unload data between variables.
	 */
	Boolean_t UnloadBetweenVars = FALSE;
	size_t TotalMem = getTotalSystemMemory();
	size_t ImportMem = 0;
	for (int i : SelectedVarNums)
		ImportMem += T41DataSizes[i - 1];
	ImportMem *= 8;

	if (ImportMem > TotalMem)
		UnloadBetweenVars = TRUE;

	if (!UnloadBetweenVars)
		TecUtilDataLoadBegin();

	if (IsOk){
		if (ReplaceDataSet){
			VarDataTypes.reserve(SelectedVarNums.size());
			for (int i : SelectedVarNums){
				if (T41DataSizes[i - 1] > 0)
					VarDataTypes.push_back(FieldDataType_Float);
				else
					VarDataTypes.push_back(FieldDataType_Bit);
			}

			IsOk = TecUtilDataSetCreate(StringMakeValidZoneName(DataSetName).c_str(), VarNames, TRUE);
		}
		else{
			if (!TecUtilDataSetIsAvailable()){
				VarDataTypes.reserve(SelectedVarNums.size());
				for (int i : SelectedVarNums){
					if (T41DataSizes[i - 1] > 0)
						VarDataTypes.push_back(FieldDataType_Float);
					else
						VarDataTypes.push_back(FieldDataType_Bit);
				}

				IsOk = TecUtilDataSetCreate(StringMakeValidZoneName(DataSetName).c_str(), VarNames, TRUE);
			}
			else{
				EntIndex_t NumVars = TecUtilDataSetGetNumVars();
				VarDataTypes.reserve(NumVars);
				for (int i = 1; i <= NumVars; ++i){
					char* TmpCStr;
					TecUtilVarGetName(i, &TmpCStr);
					int Ind = SearchVectorForString(T41LoadUserVarStrings, TmpCStr, false);
					TecUtilStringDealloc(&TmpCStr);
					if (Ind >= 0){
						if (T41DataSizes[Ind] > 0)
							VarDataTypes.push_back(FieldDataType_Float);
						else
							VarDataTypes.push_back(FieldDataType_Bit);
					}
					else
						VarDataTypes.push_back(FieldDataType_Bit);
				}
			}
		}
		if (IsOk)
			IsOk = TecUtilDataSetAddZone(StringMakeValidZoneName(CSMZoneName.FullVolume + DataSetName).c_str(), MaxI, MaxJ, MaxK, ZoneType_Ordered, VarDataTypes.data());
		if (IsOk)
			IsOk = AuxDataZoneSetItem(TecUtilDataSetGetNumZones(), DLZoneType, DLZoneTypeVolumeZone);
	}

	EntIndex_t ZoneNum, VolZoneNum;
	if (IsOk){
		VolZoneNum = TecUtilDataSetGetNumZones();
		IsOk = VolZoneNum > 0;
	}

	vector<FieldDataType_e> ZoneDataTypes(VolZoneNum, FieldDataType_Bit);
	ZoneDataTypes[VolZoneNum - 1] = FieldDataType_Float;

	Set_pa VolZoneSet = TecUtilSetAlloc(TRUE);
	if (IsOk){
		TecUtilSetAddMember(VolZoneSet, VolZoneNum, TRUE);
	}

	double TotNumPoints = MaxI * MaxJ * MaxK;

	StatusLaunch("Loading Tape41 file", AddOnID, TRUE);

	string IsoSurfaceVarName;
	EntIndex_t VarNum;
	for (int i = 0; i < SelectedVarNums.size() && IsOk; ++i){
		if (T41DataSizes[SelectedVarNums[i] - 1] > 0){
			if (UnloadBetweenVars)
				TecUtilDataLoadBegin();
			stringstream ss;
			string CurrentVarName = T41LoadUserVarStrings[SelectedVarNums[i] - 1];

			if (!ReplaceDataSet){
				VarNum = TecUtilVarGetNumByName(CurrentVarName.c_str());
				if (VarNum == TECUTILSETNOTMEMBER){
					IsOk = TecUtilDataSetAddVar(CurrentVarName.c_str(), ZoneDataTypes.data());
					if (IsOk){
						VarNum = TecUtilDataSetGetNumVars();
						VarDataTypes.push_back(FieldDataType_Float);
						IsOk = (VarNum > 0);
					}
				}
			}
			else
				VarNum = i + 1;

			ss << "Loading Tape41 file: Loading variable " << CurrentVarName;
			if (!StatusUpdate(i, SelectedVarNums.size(), ss.str(), AddOnID))
				IsOk = FALSE;
			if (IsoSurfaceVarName.empty() && strncmp("Electron Density", CurrentVarName.c_str(), 16) == 0)
				IsoSurfaceVarName = CurrentVarName;

			FieldData_pa VarRef = TecUtilDataValueGetWritableNativeRef(VolZoneNum, VarNum);
			IsOk = VALID_REF(VarRef);
			if (IsOk){
				int TempBufferMemoryKB = sizeof(double) * T41DataSizes[SelectedVarNums[i] - 1] / 1024;
				int ImportBufferMemoryKB = sizeof(ImportType_t) * T41DataSizes[SelectedVarNums[i] - 1] / 1024;
				TecUtilMemoryChangeNotify(TempBufferMemoryKB + ImportBufferMemoryKB);
				vector<double> TempBuffer(T41DataSizes[SelectedVarNums[i] - 1]);
				IsOk = (getKFData(&GlobalKFFile, T41LoadKFVarStrings[SelectedVarNums[i] - 1].c_str(), TempBuffer.data()) > 0);
				if (IsOk) {
					if (sizeof(ImportType_t) == sizeof(double)) {
						TecUtilDataValueArraySetByRef(VarRef, 1, T41DataSizes[SelectedVarNums[i] - 1], TempBuffer.data());
					}
					else{
						vector<float> TmpFloatData(TempBuffer.begin(), TempBuffer.end());
						TecUtilDataValueArraySetByRef(VarRef, 1, T41DataSizes[SelectedVarNums[i] - 1], TmpFloatData.data());
					}
				}

				TecUtilMemoryChangeNotify(-(TempBufferMemoryKB + ImportBufferMemoryKB));
			}

			if (IsOk){
				Set_pa VarSet = TecUtilSetAlloc(TRUE);
				TecUtilSetAddMember(VarSet, VarNum, TRUE);

				ArgList_pa ArgList = TecUtilArgListAlloc();
				TecUtilArgListAppendInt(ArgList, SV_STATECHANGE, StateChange_VarsAltered);
				TecUtilArgListAppendSet(ArgList, SV_ZONELIST, VolZoneSet);
				TecUtilArgListAppendSet(ArgList, SV_VARLIST, VarSet);
				TecUtilStateChangedX(ArgList);

				TecUtilArgListDealloc(&ArgList);
				TecUtilSetDealloc(&VarSet);
			}
			if (UnloadBetweenVars)
				TecUtilDataLoadEnd();
		}
	}

	if (UnloadBetweenVars)
		TecUtilDataLoadBegin();

	SetZoneStyle({ VolZoneNum }, ZoneStyle_Volume);

	if (!ReplaceDataSet){
		TecUtilZoneSetScatter(SV_SHOW, VolZoneSet, 0.0, FALSE);
		TecUtilZoneSetMesh(SV_SHOW, VolZoneSet, 0.0, FALSE);
		TecUtilZoneSetContour(SV_SHOW, VolZoneSet, 0.0, FALSE);
		TecUtilZoneSetEdgeLayer(SV_SHOW, VolZoneSet, 0.0, TRUE);
	}

	//	Modify group number of volume zone
	if (IsOk){
		ArgList_pa ArgList = TecUtilArgListAlloc();
		TecUtilArgListAppendString(ArgList, SV_P1, SV_FIELDMAP);
		TecUtilArgListAppendString(ArgList, SV_P2, SV_GROUP);
		TecUtilArgListAppendSet(ArgList, SV_OBJECTSET, VolZoneSet);
		TecUtilArgListAppendArbParam(ArgList, SV_IVALUE, (LgIndex_t)(GroupIndex));
		TecUtilStyleSetLowLevelX(ArgList);
		TecUtilArgListDealloc(&ArgList);
	}

	TecUtilSetDealloc(&VolZoneSet);

	if (ReplaceDataSet || ImportAtoms){

		int NumAtoms = -1;
		if (IsOk)
			IsOk = (getKFData(&GlobalKFFile, "Geometry%nnuc", &NumAtoms) > 0);

		vector<AtomGroup_s> AtomGroupList;

		vector<char> AtomLabels(NumAtoms * 160);

		if (IsOk)
			IsOk = (getKFData(&GlobalKFFile, "Geometry%labels", AtomLabels.data()) > 0);

		string OldLine, NewLine;

		vector<char>::iterator LabelPos = AtomLabels.begin();
		for (int i = 0; i < NumAtoms && IsOk; ++i){
			NewLine = string(LabelPos, LabelPos + 160);
			NewLine = NewLine.substr(0, NewLine.find_first_of(" "));
			if (NewLine != OldLine){
				AtomGroupList.push_back(AtomGroup_s(NewLine));
				OldLine = NewLine;
			}
			else{
				size_t AtomNum = AtomGroupList.size();
				if (AtomNum > 0)
					AtomGroupList[AtomNum - 1].Count++;
				else{
					TecUtilDialogErrMsg("Failed Parsing Atom List: Atom Type Names");
					IsOk = FALSE;
				}
			}
			LabelPos += 160;
		}

		AtomLabels.clear();

		vector<double> AtomPositions(3 * NumAtoms);
		if (IsOk)
			IsOk = (getKFData(&GlobalKFFile, "Geometry%xyznuc", AtomPositions.data()) > 0);

		int iNum = 0;
		for (int GroupNum = 0; GroupNum < AtomGroupList.size() && IsOk; ++GroupNum){
			for (int AtomNum = 0; AtomNum < AtomGroupList[GroupNum].Count; ++AtomNum){
				vec3 Position;
				for (int j = 0; j < 3; ++j){
					Position[j] = AtomPositions[iNum];
					++iNum;
				}
				AtomGroupList[GroupNum].AddPosition(Position);
			}
		}

		vector<double> AtomCharges(NumAtoms);
		if (IsOk)
			IsOk = (getKFData(&GlobalKFFile, "Geometry%qtch", AtomCharges.data()) > 0);

		iNum = 0;
		for (int GroupNum = 0; GroupNum < AtomGroupList.size() && IsOk; ++GroupNum){
			for (int AtomNum = 0; AtomNum < AtomGroupList[GroupNum].Count; ++AtomNum){
				AtomGroupList[GroupNum].Charges.push_back(AtomCharges[iNum]);
				++iNum;
			}
		}

		/*
		*	Now create zones for each atom group and set their style settings accordingly
		*/
		Set_pa NewZones = TecUtilSetAlloc(TRUE);
		vector<AtomColor_s> AtomColorList;
		PopulateAtomColorList(AtomColorList);

		CreateAtomZonesFromAtomGroupList(AtomGroupList, { "X","Y","Z" }, VarDataTypes, 50);

//  		for (int GroupNum = 0; GroupNum < AtomGroupList.size() && IsOk; ++GroupNum){
//  			IsOk = TecUtilDataSetAddZone(AtomGroupList[GroupNum].Name.c_str(), AtomGroupList[GroupNum].Count, 1, 1, ZoneType_Ordered, NULL);
//  
//  			if (IsOk){
//  				ZoneNum = TecUtilDataSetGetNumZones();
//  				IsOk = AuxDataZoneSetItem(ZoneNum, DLZoneType, DLZoneTypeNuclearPositions);
//  			}
//  			if (IsOk){
//  				IsOk = AuxDataZoneSetItem(ZoneNum, DLZoneAtomicSpecies, AtomGroupList[GroupNum].Name);
//  			}
//  
//  			if (IsOk){
//  				FieldData_pa VarRef[3] = {
//  					TecUtilDataValueGetWritableNativeRef(ZoneNum, TecUtilVarGetNumByName("X")),
//  					TecUtilDataValueGetWritableNativeRef(ZoneNum, TecUtilVarGetNumByName("Y")),
//  					TecUtilDataValueGetWritableNativeRef(ZoneNum, TecUtilVarGetNumByName("Z"))
//  				};
//  				if (VALID_REF(VarRef[0]) && VALID_REF(VarRef[0]) && VALID_REF(VarRef[0])){
//  					for (int i = 0; i < 3; ++i){
//  						vector<ImportType_t> ImportBuffer(AtomGroupList[GroupNum].Positions[i].begin(), AtomGroupList[GroupNum].Positions[i].end());
//  						// 					TecUtilDataValueArraySetByRef(VarRef[i], 1, AtomGroupList[GroupNum].Count, AtomGroupList[GroupNum].Positions[i].data());
//  						TecUtilDataValueArraySetByRef(VarRef[i], 1, AtomGroupList[GroupNum].Count, ImportBuffer.data());
//  					}
//  				}
//  
//  				TecUtilSetAddMember(NewZones, ZoneNum, FALSE);
//  
//  				Set_pa TempSet = TecUtilSetAlloc(TRUE);
//  				TecUtilSetAddMember(TempSet, ZoneNum, TRUE);
//  
//  				TecUtilZoneSetScatter(SV_SHOW, TempSet, 0.0, TRUE);
//  				TecUtilZoneSetScatter(SV_COLOR, TempSet, 0.0, GetAtomColor(AtomColorList, AtomGroupList[GroupNum].Name));
//  				TecUtilZoneSetScatter(SV_FRAMESIZE, TempSet, MIN(2 * sqrt(AtomGroupList[GroupNum].Charges[0]), 8), 0);
//  				TecUtilZoneSetScatterSymbolShape(SV_GEOMSHAPE, TempSet, GeomShape_Sphere);
//  
//  				TecUtilZoneSetMesh(SV_SHOW, TempSet, 0.0, FALSE);
//  				TecUtilZoneSetContour(SV_SHOW, TempSet, 0.0, FALSE);
//  				TecUtilZoneSetVector(SV_SHOW, TempSet, 0.0, FALSE);
//  				TecUtilZoneSetShade(SV_SHOW, TempSet, 0.0, FALSE);
//  				TecUtilZoneSetEdgeLayer(SV_SHOW, TempSet, 0.0, FALSE);
//  
//  				TecUtilSetDealloc(&TempSet);
//  			}
//  		}

		/*
		*	Inform Tecplot of the new zones
		*/
		if (IsOk){
			//	Modify group number of volume zone
			ArgList_pa ArgList = TecUtilArgListAlloc();
			TecUtilArgListAppendString(ArgList, SV_P1, SV_FIELDMAP);
			TecUtilArgListAppendString(ArgList, SV_P2, SV_GROUP);
			TecUtilArgListAppendSet(ArgList, SV_OBJECTSET, NewZones);
			TecUtilArgListAppendArbParam(ArgList, SV_IVALUE, (LgIndex_t)(GroupIndex));
			TecUtilStyleSetLowLevelX(ArgList);
			TecUtilArgListDealloc(&ArgList);

			TecUtilStateChanged(StateChange_ZonesAdded, reinterpret_cast<ArbParam_t>(NewZones));

			TecUtilFrameSetPlotType(PlotType_Cartesian3D);

			TecUtilZoneSetActive(NewZones, AssignOp_PlusEquals);

			TecUtilSetDealloc(&NewZones);

			ArgList_pa argList = TecUtilArgListAlloc();

			TecUtilArgListClear(argList);
			TecUtilArgListAppendString(argList, SV_P1, SV_ISOSURFACEATTRIBUTES);
			TecUtilArgListAppendString(argList, SV_P2, SV_EFFECTS);
			TecUtilArgListAppendString(argList, SV_P3, SV_USETRANSLUCENCY);
			TecUtilArgListAppendInt(argList, SV_OFFSET1, 1);
			TecUtilArgListAppendArbParam(argList, SV_IVALUE, TRUE);
			TecUtilStyleSetLowLevelX(argList);

			TecUtilArgListClear(argList);
			TecUtilArgListAppendString(argList, SV_P1, SV_ISOSURFACEATTRIBUTES);
			TecUtilArgListAppendString(argList, SV_P2, SV_EFFECTS);
			TecUtilArgListAppendString(argList, SV_P3, SV_SURFACETRANSLUCENCY);
			TecUtilArgListAppendInt(argList, SV_OFFSET1, 1);
			TecUtilArgListAppendArbParam(argList, SV_IVALUE, 70);
			TecUtilStyleSetLowLevelX(argList);

			TecUtilArgListClear(argList);
			TecUtilArgListAppendString(argList, SV_P1, SV_FIELDLAYERS);
			TecUtilArgListAppendString(argList, SV_P2, SV_SHOWSCATTER);
			TecUtilArgListAppendArbParam(argList, SV_IVALUE, TRUE);
			TecUtilStyleSetLowLevelX(argList);

			TecUtilArgListClear(argList);
			TecUtilArgListAppendString(argList, SV_P1, SV_ISOSURFACELAYERS);
			TecUtilArgListAppendString(argList, SV_P2, SV_SHOW);
			TecUtilArgListAppendArbParam(argList, SV_IVALUE, TRUE);
			TecUtilStyleSetLowLevelX(argList);

			TecUtilArgListClear(argList);
			TecUtilArgListAppendString(argList, SV_P1, SV_ISOSURFACEATTRIBUTES);
			TecUtilArgListAppendString(argList, SV_P2, SV_DEFINITIONCONTOURGROUP);
			TecUtilArgListAppendInt(argList, SV_OFFSET1, 1);
			TecUtilArgListAppendArbParam(argList, SV_IVALUE, (SmInteger_t)1);
			TecUtilStyleSetLowLevelX(argList);

			TecUtilArgListClear(argList);
			TecUtilArgListAppendString(argList, SV_P1, SV_ISOSURFACEATTRIBUTES);
			TecUtilArgListAppendString(argList, SV_P2, SV_SHOWGROUP);
			TecUtilArgListAppendInt(argList, SV_OFFSET1, 1);
			TecUtilArgListAppendArbParam(argList, SV_IVALUE, TRUE);
			TecUtilStyleSetLowLevelX(argList);

			EntIndex_t IsoSurfaceVarNum = TecUtilVarGetNumByName(IsoSurfaceVarName.c_str());
			// 	double ValMax, ValMin;
			// 	TecUtilVarGetMinMax(IsoSurfaceVarNum, &ValMin, &ValMax);
			// 	double IsoValue = ValMax - 0.05 * (ValMax - ValMin);
			double IsoValue = 0.25;

			TecUtilArgListClear(argList);
			TecUtilArgListAppendString(argList, SV_P1, SV_ISOSURFACEATTRIBUTES);
			TecUtilArgListAppendString(argList, SV_P2, SV_ISOVALUE1);
			TecUtilArgListAppendInt(argList, SV_OFFSET1, 1);
			TecUtilArgListAppendDouble(argList, SV_DVALUE, IsoValue);
			TecUtilStyleSetLowLevelX(argList);


			TecUtilArgListClear(argList);
			TecUtilArgListAppendString(argList, SV_P1, SV_GLOBALCONTOUR);
			TecUtilArgListAppendString(argList, SV_P2, SV_VAR);
			TecUtilArgListAppendInt(argList, SV_OFFSET1, 1);
			TecUtilArgListAppendArbParam(argList, SV_IVALUE, IsoSurfaceVarNum);
			TecUtilStyleSetLowLevelX(argList);
			TecUtilArgListClear(argList);
			TecUtilArgListAppendString(argList, SV_P1, SV_ISOSURFACEATTRIBUTES);
			TecUtilArgListAppendString(argList, SV_P2, SV_ISOVALUE1);
			TecUtilArgListAppendInt(argList, SV_OFFSET1, 1);
			TecUtilArgListAppendDouble(argList, SV_DVALUE, IsoValue);
			TecUtilStyleSetLowLevelX(argList);
			TecUtilArgListClear(argList);
			TecUtilArgListAppendString(argList, SV_P1, SV_ISOSURFACEATTRIBUTES);
			TecUtilArgListAppendString(argList, SV_P2, SV_ISOVALUE2);
			TecUtilArgListAppendInt(argList, SV_OFFSET1, 1);
			TecUtilArgListAppendDouble(argList, SV_DVALUE, IsoValue);
			TecUtilStyleSetLowLevelX(argList);
			TecUtilArgListClear(argList);
			TecUtilArgListAppendString(argList, SV_P1, SV_ISOSURFACEATTRIBUTES);
			TecUtilArgListAppendString(argList, SV_P2, SV_ISOVALUE3);
			TecUtilArgListAppendInt(argList, SV_OFFSET1, 1);
			TecUtilArgListAppendDouble(argList, SV_DVALUE, IsoValue);
			TecUtilStyleSetLowLevelX(argList);

			TecUtilArgListDealloc(&argList);
		}
	}

	StatusDrop(AddOnID);

	TecUtilDataLoadEnd();

	if (IsOk)
		TecUtilViewDataFit();

	TecUtilStringListDealloc(&VarNames);
	closeKFFile(&GlobalKFFile);

	AuxDataDataSetSetItem(DLProgramName, "SCM ADF");

// 	CalcGradGradMagForDataset(FALSE, AddOnID);

	if (IsOk && TotNumPoints > MaxPointsForShowVolumeZone){
		Set_pa TempSet = TecUtilSetAlloc(FALSE);
		IsOk = TecUtilSetAddMember(TempSet, VolZoneNum, FALSE);
		if (IsOk){
			TecUtilZoneSetActive(TempSet, AssignOp_MinusEquals);
		}
		TecUtilSetDealloc(&TempSet);
	}

	TecUtilLockFinish(AddOnID);


	return;
}

void NumCellsLoadData(){
	if (FileType == BANDFile)
		LoadBANDTape41Data();
	else
		LoadVASPData();
}

void SelVarsLoadData(){
	if (FileType == ADFFile)
		LoadADFTape41Files();
	else
		NumCellsLoadData();
// 		TecGUIDialogLaunch(Dialog2Manager);
}

void LoadBANDTape41Data(){

	TecUtilLockStart(AddOnID);

	SYSTEM_INFO sysinfo;
	GetSystemInfo(&sysinfo);

	int numCPU = sysinfo.dwNumberOfProcessors;
	omp_set_num_threads(numCPU);

	LgIndex_t *SelectedVarNumsArray = nullptr;
	LgIndex_t NumSelected = -1;
	TecGUIListGetSelectedItems(MLT41LoadVa_MLST_D1, &SelectedVarNumsArray, &NumSelected);

	if (NumSelected < 1){
		TecUtilDialogErrMsg("You must select at least one variable to load from the Tape41 file");
		TecUtilArrayDealloc(reinterpret_cast<void**>(&SelectedVarNumsArray));
		TecGUIDialogLaunch(Dialog1Manager);
		TecUtilLockFinish(AddOnID);
		return;
	}

	string DataSetName;
	vector<LgIndex_t> SelectedVarNums(SelectedVarNumsArray, SelectedVarNumsArray + NumSelected);
	TecUtilArrayDealloc(reinterpret_cast<void**>(&SelectedVarNumsArray));


#if defined(WIN32) || defined(_WIN32) || defined(__WIN32) && !defined(__CYGWIN__)
	string Slash = "\\";
#else
	string Slash = "/";
#endif

	size_t StrPos = FileNameStr.find_last_of("/\\");
	if (StrPos != string::npos){
		DataSetName = FileNameStr.substr(StrPos + 1, FileNameStr.length());
	}
	StrPos = DataSetName.find_last_of(".");
	if (StrPos != string::npos){
		DataSetName = DataSetName.substr(0, StrPos);
	}

	Boolean_t IsOk = TRUE;

	/*
	*	Open Tape41 kf file
	*/

	IsOk = (openKFFile(&GlobalKFFile, const_cast<char*>(FileNameStr.c_str())) > 0);

	if (!IsOk){
		TecUtilDialogErrMsg("Failed to open Tape41 file.");
		return;
	}

	int MaxI = -1, MaxJ = -1, MaxK = -1, NumPoints = -1;
	int MaxIJK[3];

	vector<string> TempStrings = { "Grid%nr of points x", "Grid%nr of points y", "Grid%nr of points z" };

	for (int i = 0; i < 3 && IsOk; ++i)
		IsOk = (getKFData(&GlobalKFFile, TempStrings[i].c_str(), &MaxIJK[i]) > 0);

	if (IsOk){
		MaxI = MaxIJK[0];
		MaxJ = MaxIJK[1];
		MaxK = MaxIJK[2];
		IsOk = (getKFData(&GlobalKFFile, "Grid%total nr of points", (void*)&NumPoints) > 0);
	}

	if (!IsOk || MaxI < 0 || MaxJ < 0 || MaxK < 0 || NumPoints != MaxI*MaxJ*MaxK){
		TecUtilDialogErrMsg("Failed to get system extent");
	}

	TempStrings = { "X", "Y", "Z" };
	T41LoadUserVarStrings.insert(T41LoadUserVarStrings.begin(), TempStrings.cbegin(), TempStrings.cend());

	TempStrings = { "x values%x values", "y values%y values", "z values%z values" };
	T41LoadKFVarStrings.insert(T41LoadKFVarStrings.begin(), TempStrings.cbegin(), TempStrings.cend());
	TempStrings.clear();

	mat33 LatticeVector;
	vec3 StartPoint;
// 	double StartPoint[3];
	if (IsOk)
		IsOk = (getKFData(&GlobalKFFile, "Grid%Start_point", StartPoint.memptr()) > 0);
	TempStrings = { "Grid%x-vector", "Grid%y-vector", "Grid%z-vector" };
	for (int i = 0; i < 3 && IsOk; ++i){
		double TmpDbl[3];
		IsOk = (getKFData(&GlobalKFFile, TempStrings[i].c_str(), TmpDbl) > 0);
		if(IsOk){
			for (int j = 0; j < 3; ++j)
				LatticeVector.at(i, j) = TmpDbl[j];
		}
	}
	TempStrings.clear();

	/*
	 *	Multiply elements of lattice vectors by their respective number of points to get
	 *	the full lattice vector (containing direction and magnitude information) back.
	 *	And convert to bohr from angstroms.
	 */
// 	for (int i = 0; i < 3; ++i)
// 		for (int j = 0; j < 3; ++j)
// 			LatticeVector.at(i, j) *= static_cast<double>(MaxIJK[i]) / BorhToAngstrom;

// 	vector<double> Mags(3);
// 	for (int i = 0; i < 3; ++i)
// 		Mags[i] = sqrt(pow(LatticeVector.at(i, 0), 2) + pow(LatticeVector.at(i, 1), 2) + pow(LatticeVector.at(i, 2), 2));

	vector<LgIndex_t> TempLgIndex = { 1, 2, 3 };
	for (int i = 0; i < SelectedVarNums.size(); ++i)
		SelectedVarNums[i] += 3;
	SelectedVarNums.insert(SelectedVarNums.begin(), TempLgIndex.begin(), TempLgIndex.end());
	TempLgIndex.clear();

	T41DataSizes.insert(T41DataSizes.begin(), 3, T41DataSizes[0]);

	int NCells = 1;
	int NCellsX = NCells, NCellsY = NCells, NCellsZ = NCells;

// 	TecGUITextFieldGetLgIndex(XNC_TFS_D2, &NCellsX);
// 	TecGUITextFieldGetLgIndex(YNC_TFS_D2, &NCellsY);
// 	TecGUITextFieldGetLgIndex(ZNC_TFS_D2, &NCellsZ);

	LgIndex_t TotNumPoints = MaxI * MaxJ * MaxK;

	LgIndex_t TMx = MaxI * NCellsX;
	LgIndex_t TMy = MaxJ * NCellsY;
	LgIndex_t TMz = MaxK * NCellsZ;
	TotNumPoints = TMx * TMy * TMz;

	StringList_pa VarNames = TecUtilStringListAlloc();
	for (int i = 0; i < SelectedVarNums.size() && IsOk; ++i)
		TecUtilStringListAppendString(VarNames, T41LoadUserVarStrings[SelectedVarNums[i] - 1].c_str());

	vector<FieldDataType_e> VarDataTypes;
	VarDataTypes.resize(SelectedVarNums.size(), FieldDataType_Float);

	if (IsOk && TecUtilDataSetCreate(StringMakeValidZoneName(DataSetName).c_str(), VarNames, TRUE)){
		IsOk = TecUtilDataSetAddZone(StringMakeValidZoneName(CSMZoneName.FullVolume + DataSetName).c_str(), TMx, TMy, TMz, ZoneType_Ordered, VarDataTypes.data());
	}

// 	EntIndex_t ZoneNum;
	EntIndex_t VolZoneNum;
	if (IsOk)
		VolZoneNum = TecUtilDataSetGetNumZones();


	SetZoneStyle({ VolZoneNum }, ZoneStyle_Volume);

	StatusLaunch("Loading Tape41 file", AddOnID, TRUE);

	string IsoSurfaceVarName;

	vector<vector<double> > VarData(SelectedVarNums.size() - 3);

	for (int i = 3; i < SelectedVarNums.size() && IsOk; ++i){
		stringstream ss;
		ss << "Loading Tape41 file: Loading variable " << T41LoadUserVarStrings[SelectedVarNums[i] - 1];
		if (!StatusUpdate(i, SelectedVarNums.size(), ss.str(), AddOnID))
			IsOk = FALSE;
		string CurrentVarName = T41LoadUserVarStrings[SelectedVarNums[i] - 1];
		if (IsoSurfaceVarName.empty() && strncmp("Electron Density", T41LoadUserVarStrings[SelectedVarNums[i] - 1].c_str(), 16) == 0)
			IsoSurfaceVarName = T41LoadUserVarStrings[SelectedVarNums[i] - 1];
		if (IsOk){
			int TempBufferMemoryKB = sizeof(double) * T41DataSizes[SelectedVarNums[i] - 1] / 1024;
			VarData[i - 3].resize(T41DataSizes[SelectedVarNums[i] - 1]);
			TecUtilMemoryChangeNotify(TempBufferMemoryKB);
			IsOk = (getKFData(&GlobalKFFile, T41LoadKFVarStrings[SelectedVarNums[i] - 1].c_str(), VarData[i - 3].data()) > 0);
		}
	}

	size_t NumVars = SelectedVarNums.size();

	vector<FieldDataPointer_c> VarRawPtrs(NumVars);
	for (EntIndex_t VarNum = 0; VarNum < NumVars; ++VarNum)
		VarRawPtrs[VarNum].InitializeWritePtr(VolZoneNum, VarNum + 1);

	int TotalPoints = MaxK * NCellsZ / numCPU;

	Boolean_t TaskQuit = FALSE;

	TecUtilDataLoadBegin();

	StatusUpdate(0, 1, "Loading data into dataset...", AddOnID);

	for (int Zz = 0; Zz < NCellsZ; ++Zz){
#pragma omp parallel for
		for (LgIndex_t k = 1; k <= MaxK; ++k){
			if (omp_get_thread_num() == 0 && !StatusUpdate(k + MaxK * Zz, TotalPoints, "Loading data into dataset...", AddOnID)){
				TaskQuit = TRUE;
#pragma omp flush (TaskQuit)
			}
#pragma omp flush (TaskQuit)
			for (int Zy = 0; Zy < NCellsY && !TaskQuit; ++Zy){
				for (LgIndex_t j = 1; j <= MaxJ; ++j){
					for (int Zx = 0; Zx < NCellsX; ++Zx){
						for (LgIndex_t i = 1; i <= MaxI; ++i){
							LgIndex_t IJK[3] = { i - 1, j - 1, k - 1 };
							vec TotIJK;
							TotIJK << i + Zx * MaxI - 1 << j + Zy * MaxJ - 1 << k + Zz * MaxK - 1;
// 							LgIndex_t TotIJK[3] = { i + Zx * MaxI - 1, j + Zy * MaxJ - 1, k + Zz * MaxK - 1 };
							LgIndex_t Index = IndexFromIJK(i + Zx*MaxI, j + Zy*MaxJ, k + Zz*MaxK, TMx, TMy, TMz, TRUE) - 1;
							LgIndex_t SmIndex = IndexFromIJK(i, j, k, MaxI, MaxJ, MaxK, TRUE) - 1;

							vec3 XYZ = StartPoint + LatticeVector * TotIJK;
							for (int VarNum = 0; VarNum < 3; ++VarNum){
								VarRawPtrs[VarNum].Write(Index, XYZ[VarNum]);
							}

							for (EntIndex_t VarNum = 3; VarNum < NumVars; ++VarNum){


// 								double Value;
// 								if (VarNum < 3){
// 									/*
// 									*	Here, for the generation of x,y,z values, using
// 									*	the lattice vectors so that non-cartesian systems
// 									*	are property represented as well.
// 									*/
// 									Value = StartPoint[VarNum];
// 									for (int Dir = 0; Dir < 3; ++Dir) {
// 										Value += static_cast<double>(TotIJK[Dir]) / static_cast<double>(MAX(MaxIJK[Dir] - 1, 1)) * LatticeVector.at(Dir, VarNum);
// 									}
// 								}
// 								else {
// 									Value = VarData[VarNum - 3][SmIndex];
// 								}
// 
// 								VarRawPtrs[VarNum].Write(Index, Value);

								VarRawPtrs[VarNum].Write(Index, VarData[VarNum - 3][SmIndex]);
							}
						}
					}
				}
			}
		}
	}

	StatusDrop(AddOnID);
	TecUtilDataLoadEnd();
	if (TaskQuit){
		TecUtilLockFinish(AddOnID);
		return;
	}



	if (IsOk){
		Set_pa TempSet = TecUtilSetAlloc(FALSE);
		IsOk = TecUtilSetAddMember(TempSet, VolZoneNum, FALSE);
		if (IsOk){
			TecUtilZoneSetActive(TempSet, AssignOp_PlusEquals);
		}
		TecUtilSetDealloc(&TempSet);
	}

	TecUtilFrameSetPlotType(PlotType_Cartesian3D);

	Set_pa ChangedVars = TecUtilSetAlloc(FALSE);
	for (EntIndex_t i = 1; i <= NumVars; ++i)
		TecUtilSetAddMember(ChangedVars, i, FALSE);
	TecUtilStateChanged(StateChange_VarsAltered, (ArbParam_t)ChangedVars);
	TecUtilSetDealloc(&ChangedVars);

	int NumAtoms = -1;
	if (IsOk)
		IsOk = (getKFData(&GlobalKFFile, "Geometry%nnuc", &NumAtoms) > 0);

	vector<AtomGroup_s> AtomGroupList;

	vector<char> AtomLabels(NumAtoms * 160);

	if (IsOk)
		IsOk = (getKFData(&GlobalKFFile, "Geometry%labels", AtomLabels.data()) > 0);

	string OldLine, NewLine;

	vector<char>::iterator LabelPos = AtomLabels.begin();
	for (int i = 0; i < NumAtoms && IsOk; ++i){
		NewLine = string(LabelPos, LabelPos + 160);
		NewLine = NewLine.substr(0, NewLine.find_first_of(" "));
		if (NewLine != OldLine){
			AtomGroupList.push_back(AtomGroup_s(NewLine));
			OldLine = NewLine;
		}
		else{
			size_t AtomNum = AtomGroupList.size();
			if (AtomNum > 0)
				AtomGroupList[AtomNum - 1].Count++;
			else{
				TecUtilDialogErrMsg("Failed Parsing Atom List: Atom Type Names");
				IsOk = FALSE;
			}
		}
		LabelPos += 160;
	}

	AtomLabels.clear();

	vector<double> AtomPositions(3 * NumAtoms);
	if (IsOk)
		IsOk = (getKFData(&GlobalKFFile, "Geometry%xyznuc", AtomPositions.data()) > 0);
	
	TecUtilDataLoadBegin();
	FieldDataPointer_c RhoPtr;
	VolExtentIndexWeights_s VolInfo;
	if (VarNumByName("Electron Density", true) > 0) {
		RhoPtr.InitializeReadPtr(TecUtilDataSetGetNumZones(), VarNumByName("Electron Density", true));
	}
	else if (VarNumByName("rho", true) > 0) {
		RhoPtr.InitializeReadPtr(TecUtilDataSetGetNumZones(), VarNumByName("rho", true));
	}
	else if (VarNumByName("denstiy", true) > 0) {
		RhoPtr.InitializeReadPtr(TecUtilDataSetGetNumZones(), VarNumByName("density", true));
	}

	if (RhoPtr.IsReady()) {
		GetVolInfo(TecUtilDataSetGetNumZones(), { 1,2,3 }, FALSE, VolInfo);
	}

	int iNum = 0;
	for (int GroupNum = 0; GroupNum < AtomGroupList.size() && IsOk; ++GroupNum){
		for (int AtomNum = 0; AtomNum < AtomGroupList[GroupNum].Count; ++AtomNum){
			vec3 Position;
			for (int j = 0; j < 3; ++j){
				Position[j] = AtomPositions[iNum];
				++iNum;
			}
			AtomGroupList[GroupNum].AddPosition(Position);

			double CPtmpCharge = 16;
			if (RhoPtr.IsReady()) {
				CPtmpCharge = MAX(16, RhoPtr.At(Position, VolInfo));
			}

			AtomGroupList[GroupNum].Charges.push_back(CPtmpCharge);
		}
	}
	TecUtilDataLoadEnd();


	/*
	*	Now add ordered zones of scatter spheres for each atom
	*	type. Two considerations: system may be repeated for
	*	multiple cells, and atoms may lie at the edge of periodic
	*	boundaries. So this is written to allow both of those things.
	*/

	/*
	*	Get system boundaries for periodic atoms
	*/
	double XYZMin[3], XYZMax[3];

	for (int i = 0; i < 3; ++i)
		TecUtilVarGetMinMax(i + 1, &XYZMin[i], &XYZMax[i]);

	int NCellsXYZ[3] = { NCellsX, NCellsY, NCellsZ };

	ImportType_t StepXYZ[3];

	Set_pa NewZones = TecUtilSetAlloc(FALSE);

	TecUtilDataLoadBegin();

	//IsOk = FALSE;

	mat33 LVTranspose = LatticeVector.t();
	int ColorCount = 0;
	for (int GroupNum = 0; GroupNum < AtomGroupList.size() && IsOk; ++GroupNum){
		vector<vector<ImportType_t> > XYZ(3, vector<ImportType_t>());
		vector<ImportType_t> Charges;
		for (int AtomNum = 0; AtomNum < AtomGroupList[GroupNum].Count && IsOk; ++AtomNum){
			vec3 TmpLatPos;
			for (int i = 0; i < 3; ++i)
				TmpLatPos[i] = AtomGroupList[GroupNum].Positions[i][AtomNum];

			vec3 CPTmpXYZ;

			for (int kk = 0; kk < NCellsZ + int(abs(TmpLatPos[2]) < 0.01); ++kk){
				for (int jj = 0; jj < NCellsY + int(abs(TmpLatPos[1]) < 0.01); ++jj){
					for (int ii = 0; ii < NCellsX + int(abs(TmpLatPos[0]) < 0.01); ++ii){
						vec3 IJK;
						IJK << ii * MaxI << jj * MaxJ << kk * MaxK;

						vec3 Pos = (LatticeVector * IJK) + TmpLatPos;
						for (int dir = 0; dir < 3; ++dir){
							XYZ[dir].push_back(Pos[dir]);
						}
						Charges.push_back(AtomGroupList[GroupNum].Charges[AtomNum]);

						if (NCellsX == 1){
							break;
						}
					}
					if (NCellsY == 1) {
						break;
					}
				}
				if (NCellsZ == 1) {
					break;
				}
			}
		}
		// Now copy charge and position to the AtomGroupList
		if (!Charges.empty()) {
// 			AtomGroupList[GroupNum].Charges.insert(AtomGroupList[GroupNum].Charges.end(), Charges.begin(), Charges.end());
			AtomGroupList[GroupNum].Charges = Charges;
			for (int i = 0; i < 3; ++i) {
// 				AtomGroupList[GroupNum].Positions[i].insert(AtomGroupList[GroupNum].Positions[i].end(), XYZ[i].begin(), XYZ[i].end());
				AtomGroupList[GroupNum].Positions[i] = XYZ[i];
			}
			AtomGroupList[GroupNum].Count = Charges.size();
		}
	}

	CreateAtomZonesFromAtomGroupList(AtomGroupList, { "X","Y","Z" }, VarDataTypes, 50);

	/*
	*	Inform Tecplot of the new zones
	*/
	if (IsOk){
		TecUtilStateChanged(StateChange_ZonesAdded, (ArbParam_t)NewZones);

		TecUtilFrameSetPlotType(PlotType_Cartesian3D);

		TecUtilZoneSetActive(NewZones, AssignOp_PlusEquals);

		TecUtilSetDealloc(&NewZones);

		ArgList_pa argList = TecUtilArgListAlloc();

		TecUtilArgListClear(argList);
		TecUtilArgListAppendString(argList, SV_P1, SV_ISOSURFACEATTRIBUTES);
		TecUtilArgListAppendString(argList, SV_P2, SV_EFFECTS);
		TecUtilArgListAppendString(argList, SV_P3, SV_USETRANSLUCENCY);
		TecUtilArgListAppendInt(argList, SV_OFFSET1, 1);
		TecUtilArgListAppendArbParam(argList, SV_IVALUE, TRUE);
		TecUtilStyleSetLowLevelX(argList);

		TecUtilArgListClear(argList);
		TecUtilArgListAppendString(argList, SV_P1, SV_ISOSURFACEATTRIBUTES);
		TecUtilArgListAppendString(argList, SV_P2, SV_EFFECTS);
		TecUtilArgListAppendString(argList, SV_P3, SV_SURFACETRANSLUCENCY);
		TecUtilArgListAppendInt(argList, SV_OFFSET1, 1);
		TecUtilArgListAppendArbParam(argList, SV_IVALUE, 70);
		TecUtilStyleSetLowLevelX(argList);

		TecUtilArgListClear(argList);
		TecUtilArgListAppendString(argList, SV_P1, SV_FIELDLAYERS);
		TecUtilArgListAppendString(argList, SV_P2, SV_SHOWSCATTER);
		TecUtilArgListAppendArbParam(argList, SV_IVALUE, TRUE);
		TecUtilStyleSetLowLevelX(argList);

		TecUtilArgListClear(argList);
		TecUtilArgListAppendString(argList, SV_P1, SV_ISOSURFACELAYERS);
		TecUtilArgListAppendString(argList, SV_P2, SV_SHOW);
		TecUtilArgListAppendArbParam(argList, SV_IVALUE, TRUE);
		TecUtilStyleSetLowLevelX(argList);

		TecUtilArgListClear(argList);
		TecUtilArgListAppendString(argList, SV_P1, SV_ISOSURFACEATTRIBUTES);
		TecUtilArgListAppendString(argList, SV_P2, SV_DEFINITIONCONTOURGROUP);
		TecUtilArgListAppendInt(argList, SV_OFFSET1, 1);
		TecUtilArgListAppendArbParam(argList, SV_IVALUE, (SmInteger_t)1);
		TecUtilStyleSetLowLevelX(argList);

		TecUtilArgListClear(argList);
		TecUtilArgListAppendString(argList, SV_P1, SV_ISOSURFACEATTRIBUTES);
		TecUtilArgListAppendString(argList, SV_P2, SV_SHOWGROUP);
		TecUtilArgListAppendInt(argList, SV_OFFSET1, 1);
		TecUtilArgListAppendArbParam(argList, SV_IVALUE, TRUE);
		TecUtilStyleSetLowLevelX(argList);

		EntIndex_t IsoSurfaceVarNum = TecUtilVarGetNumByName(IsoSurfaceVarName.c_str());
		// 	double ValMax, ValMin;
		// 	TecUtilVarGetMinMax(IsoSurfaceVarNum, &ValMin, &ValMax);
		// 	double IsoValue = ValMax - 0.05 * (ValMax - ValMin);
		double IsoValue = 0.25;

		TecUtilArgListClear(argList);
		TecUtilArgListAppendString(argList, SV_P1, SV_ISOSURFACEATTRIBUTES);
		TecUtilArgListAppendString(argList, SV_P2, SV_ISOVALUE1);
		TecUtilArgListAppendInt(argList, SV_OFFSET1, 1);
		TecUtilArgListAppendDouble(argList, SV_DVALUE, IsoValue);
		TecUtilStyleSetLowLevelX(argList);


		TecUtilArgListClear(argList);
		TecUtilArgListAppendString(argList, SV_P1, SV_GLOBALCONTOUR);
		TecUtilArgListAppendString(argList, SV_P2, SV_VAR);
		TecUtilArgListAppendInt(argList, SV_OFFSET1, 1);
		TecUtilArgListAppendArbParam(argList, SV_IVALUE, IsoSurfaceVarNum);
		TecUtilStyleSetLowLevelX(argList);
		TecUtilArgListClear(argList);
		TecUtilArgListAppendString(argList, SV_P1, SV_ISOSURFACEATTRIBUTES);
		TecUtilArgListAppendString(argList, SV_P2, SV_ISOVALUE1);
		TecUtilArgListAppendInt(argList, SV_OFFSET1, 1);
		TecUtilArgListAppendDouble(argList, SV_DVALUE, IsoValue);
		TecUtilStyleSetLowLevelX(argList);
		TecUtilArgListClear(argList);
		TecUtilArgListAppendString(argList, SV_P1, SV_ISOSURFACEATTRIBUTES);
		TecUtilArgListAppendString(argList, SV_P2, SV_ISOVALUE2);
		TecUtilArgListAppendInt(argList, SV_OFFSET1, 1);
		TecUtilArgListAppendDouble(argList, SV_DVALUE, IsoValue);
		TecUtilStyleSetLowLevelX(argList);
		TecUtilArgListClear(argList);
		TecUtilArgListAppendString(argList, SV_P1, SV_ISOSURFACEATTRIBUTES);
		TecUtilArgListAppendString(argList, SV_P2, SV_ISOVALUE3);
		TecUtilArgListAppendInt(argList, SV_OFFSET1, 1);
		TecUtilArgListAppendDouble(argList, SV_DVALUE, IsoValue);
		TecUtilStyleSetLowLevelX(argList);

		TecUtilArgListDealloc(&argList);
	}

	TecUtilDataLoadEnd();

	if (IsOk)
		TecUtilViewDataFit();

	TecUtilStringListDealloc(&VarNames);
	closeKFFile(&GlobalKFFile);
	QuitT41Load();

	AuxDataDataSetSetItem(DLProgramName, "SCM BAND");
	for (int i = 0; i < 3; ++i){
		stringstream ss;
		ss << "{ ";
		for (int j = 0; j < 3; ++j){
			ss << LatticeVector.at(i, j);
			if (j < 2) ss << ", ";
		}
		ss << " } [bohr]";
		AuxDataDataSetSetItem(DLLatticeVecs[i], ss.str());
	}
	{
		stringstream ss;
		ss << "{ ";
		for (int i = 0; i < 3; ++i){
			ss << StartPoint[i];
			if (i < 2) ss << ", ";
		}
		ss << " } [bohr]";
		AuxDataDataSetSetItem(DLOrigin, ss.str());
	}

// 	CalcGradGradMagForDataset(TRUE, AddOnID);

	if (IsOk && TotNumPoints > MaxPointsForShowVolumeZone){
		Set_pa TempSet = TecUtilSetAlloc(FALSE);
		IsOk = TecUtilSetAddMember(TempSet, VolZoneNum, FALSE);
		if (IsOk){
			TecUtilZoneSetActive(TempSet, AssignOp_MinusEquals);
		}
		TecUtilSetDealloc(&TempSet);
	}

	TecUtilLockFinish(AddOnID);


	return;
}


void SpinButtonInt(LgIndex_t DialogIndex, LgIndex_t Difference){
	LgIndex_t TheValue;
	if (TecGUITextFieldGetLgIndex(DialogIndex, &TheValue)){
		TheValue += Difference;
		if (TheValue <= 0) TheValue = 1;
		TecGUITextFieldSetLgIndex(DialogIndex, TheValue, FALSE);
	}
}

void SpinValueChangedInt(LgIndex_t DialogIndex){
	LgIndex_t Junk;
	if (!TecGUITextFieldGetLgIndex(DialogIndex, &Junk)){
		TecUtilDialogErrMsg("Only interger values are allowed.");
		TecGUITextFieldSetLgIndex(DialogIndex, 1, FALSE);
	}
}


int dos2unix(string const & FileName)
{
	char ch;
	char temp[MAX_PATH] = "\0";

	//Open the file for reading in binarymode.
	std::ifstream fp_read(FileName, ios::in \
		| ios::binary);
	sprintf(temp, "%s.temp", FileName);
	//Create a temporary file for writing in the binary mode. This
	//file will be created in the same directory as the input file.
	ofstream fp_write(temp, ios::out \
		| ios::trunc \
		| ios::binary);

	while (fp_read.eof() != TRUE)
	{
		fp_read.get(ch);
		//Check for CR (carriage return)
		if ((int)ch == 0x0D)
			continue;
		if (!fp_read.eof())fp_write.put(ch);
	}

	fp_read.close();
	fp_write.close();
	//Delete the existing input file.
	remove(FileName.c_str());
	//Rename the temporary file to the input file.
	rename(temp, FileName.c_str());
	//Delete the temporary file.
	//remove(temp);

	return 0;
}

int dos2unix2(string const & FileName)
{
	char ch;
	char temp[MAX_PATH] = "\0";

	ifstream InFile(FileName, ios::in | ios::binary);

	//	Get file size
// 	InFile.seekg(0, InFile.end);
// 	std::streamoff FileSize = InFile.tellg();
// 	InFile.seekg(0, InFile.beg);
// 
// 	//	Read in entire file to FileData
// 	vector<char> FileData(FileSize);
// 	InFile.read(FileData.data(), FileSize);
// 	InFile.close();

	//string OutStr;

	ofstream OutFile( FileName + ".new.sh", ios::out | ios::trunc | ios::binary);

	while (!InFile.eof()){
		InFile.get(ch);

		if ((int)ch == 0x0D)
			continue;

		if (!InFile.eof()) OutFile.put(ch);
	}

// 	for (char aChar : FileData){
// 		if ((int)aChar == 0x0D)
// 			continue;
// 
// 		OutFile.put(aChar);
// 		//OutStr += aChar;
// 	}


	//OutFile << OutStr;
	InFile.close();
	OutFile.close();

	remove(FileName.c_str());

	return 0;
}

template <typename T>
std::string to_string_with_precision(T const a_value, int const n = 6)
{
	std::ostringstream out;
	out << std::setprecision(n) << a_value;
	return out.str();
}

void MakeDensfScriptForZones(){
	/*
	*	Get the list of selected zones
	*/

	char *TmpName;
	vector<string> ZoneNameList;

	int NumZonesSelected;
	int *SelectedListNums;
	TecGUIListGetSelectedItems(ZoneList_MLST_T1_1, &SelectedListNums, &NumZonesSelected);

	if (NumZonesSelected <= 0){
		TecUtilDialogErrMsg("Select 1 or more zones first.");
		return;
	}

	ZoneNameList.reserve(NumZonesSelected);

	for (int i = 0; i < NumZonesSelected; ++i){
		TmpName = TecGUIListGetString(ZoneList_MLST_T1_1, SelectedListNums[i]);
		string TmpStr = TmpName;
		TecUtilStringDealloc(&TmpName);
		if (TmpStr.find("Iso") > TmpStr.length()){
			ZoneNameList.push_back(TmpStr);
		}
	}
	TecUtilArrayDealloc(reinterpret_cast<void**>(&SelectedListNums));

	NumZonesSelected = static_cast<int>(ZoneNameList.size());

	/*
	* Get dataset info
	*/
	EntIndex_t NumZones;
	TecUtilDataSetGetInfo(&TmpName, &NumZones, nullptr);
	string DataSetName = TmpName;
	TecUtilStringDealloc(&TmpName);

	vector<EntIndex_t> XYZVarNums(3, -1);
	TecUtilAxisGetVarAssignments(&XYZVarNums[0], &XYZVarNums[1], &XYZVarNums[2]);
	if (XYZVarNums[0] <= 0 || XYZVarNums[1] <= 0 || XYZVarNums[2] <= 0){
		vector<vector<string> > TmpStrs = {
			{ "X", "Y", "Z" }, { "I", "J", "K" }
		};
		Boolean_t VarsFound = FALSE;
		for (int i = 0; i < 2 && !VarsFound; ++i){
			for (int j = 0; j < 3; ++j)
				XYZVarNums[j] = TecUtilVarGetNumByName(TmpStrs[i][j].c_str());
			VarsFound = (XYZVarNums[0] > 0 && XYZVarNums[1] > 0 && XYZVarNums[2] > 0);
		}
		if (!VarsFound){
			TecUtilDialogErrMsg("Failed to get x, y, z variable numbers.");
			return;
		}
	}

	/*
	*	Get path for created script file
	*/
	string OutFileName = DataSetName + "_Densf_Script.sh";
	if (!TecUtilDialogGetFileName(SelectFileOption_WriteFile, &TmpName, "Bash shell script", OutFileName.c_str(), "*.sh")){
		return;
	}
	OutFileName = TmpName;
	TecUtilStringDealloc(&TmpName);

	ofstream OutFile(OutFileName, ios::out | ios::trunc);
	if (!OutFile.is_open()){
		TecUtilDialogErrMsg("Failed to open script file for writing.");
		return;
	}

	/*
	* Start the real work
	*/
	LgIndex_t Index;

	/*
	*	Write header for script file. This is where the user can specify the Tape21 file path
	*	or include it as an argument.
	*/

	OutFile << "INPUTFILE=\'\'" << endl << endl

		<< "CALCVARS=\'Density SCF" << endl
		<< "DenGrad" << endl
		<< "DenHess\'" << endl << endl

		<< "if [ -z \"$INPUTFILE\" ] && [ $# -lt 1 ]; then" << endl
		<< "	echo \"Please enter the full path to a Tape21 file either as the only argument to the script or by specifying it in the first line of the script.\"" << endl
		<< "	echo \"Quitting...\"" << endl
		<< "	exit" << endl
		<< "elif [ -z \"$INPUTFILE\" ]; then" << endl
		<< "	INPUTFILE=$1" << endl
		<< "fi" << endl << endl

		<< "INPUTPATH=`dirname \"$INPUTFILE\"`" << endl
		<< "INPUTBASE=`basename \"$INPUTFILE\"`" << endl
		<< "NEWPATH=\"$INPUTPATH/$INPUTBASE.TecplotZoneTape41Files\"" << endl << endl

		<< "echo $INPUTPATH" << endl
		<< "echo $INPUTBASE" << endl
		<< "echo $NEWPATH" << endl

		<< "mkdir \"$NEWPATH\"" << endl << endl;

	/*
	*	Now loop over each selected zone and add the densf call to the script for it.
	*/
	Boolean_t IsOk = TRUE;

	string StatusStr = "Creating densf script";

	StatusLaunch(StatusStr.c_str(), AddOnID, TRUE);

	vector<string> ZoneStrs(NumZonesSelected);
	vector<vector<FieldDataPointer_c> > XYZPtrs(NumZonesSelected, vector<FieldDataPointer_c>(3));

	TecUtilDataLoadBegin();

	for (int SelZoneNum = 0; SelZoneNum < NumZonesSelected; ++SelZoneNum){
		int ZoneNum = ZoneNumByName(ZoneNameList[SelZoneNum]);

		for (int i = 0; i < 3 && IsOk; ++i){
			IsOk = XYZPtrs[SelZoneNum][i].InitializeReadPtr(ZoneNum, XYZVarNums[i]);
		}
		int NumSigFigs = (XYZPtrs[SelZoneNum][0].FDType() == FieldDataType_Double ? 18 : 10);

		ZoneStrs[SelZoneNum].reserve((XYZPtrs[SelZoneNum][0].Size() * NumSigFigs * 3 + 5) + 500);
	}

#ifndef _DEBUG
#pragma omp parallel for schedule(dynamic)
#endif
	for (int SelZoneNum = 0; SelZoneNum < NumZonesSelected; SelZoneNum++){
		EntIndex_t ZoneNum = XYZPtrs[SelZoneNum][0].ZoneNum();
		if (omp_get_thread_num() == 0){
			IsOk = StatusUpdate(SelZoneNum, NumZonesSelected, StatusStr, AddOnID);
#pragma omp flush (IsOk)
		}
#pragma omp flush (IsOk)
		if (ZoneNum >= 0){
			vector<LgIndex_t> MaxIJK = XYZPtrs[SelZoneNum][0].MaxIJK();
			TecUtilZoneGetIJK(ZoneNum, &MaxIJK[0], &MaxIJK[1], &MaxIJK[2]);
			int NumSigFigs = (XYZPtrs[SelZoneNum][0].FDType() == FieldDataType_Double ? 18 : 10);

			if (IsOk){

				ZoneStrs[SelZoneNum] += string("$ADFBIN/densf << eor\n\n")

					+ string("INPUTFILE $INPUTFILE\n")
					+ string("OUTPUTFILE $NEWPATH/") + DataSetName + DS_ZoneName_Delim + ZoneNameList[SelZoneNum] + string(".t41\n\n")

					+ string("GRID Inline\n");

// 				OutFile << "$ADFBIN/densf << eor" << endl << endl
// 
// 					<< "INPUTFILE $INPUTFILE" << endl
// 					<< "OUTPUTFILE $NEWPATH/" << DataSetName << DS_ZoneName_Delim << ZoneNameList[SelZoneNum] << ".t41" << endl << endl
// 
// 					<< "GRID Inline" << endl;

				if (XYZPtrs[SelZoneNum][0].ZoneIsOrdered()){
					for (int kk = 1; kk <= MaxIJK[2] && IsOk; ++kk){
						for (int jj = 1; jj <= MaxIJK[1] && IsOk; ++jj){
							for (int ii = 1; ii <= MaxIJK[0] && IsOk; ++ii){
								Index = IndexFromIJK(ii, jj, kk, MaxIJK[0], MaxIJK[1]) - 1;
								ZoneStrs[SelZoneNum] += "\t";
// 								OutFile << "\t";
								for (int Dir = 0; Dir < 3; ++Dir)
									ZoneStrs[SelZoneNum] += to_string_with_precision(XYZPtrs[SelZoneNum][Dir][Index], NumSigFigs) + " ";
// 								OutFile << endl;
								ZoneStrs[SelZoneNum] += "\n";
							}
						}
					}
				}
				else{
					for (int ii = 1; ii <= MaxIJK[0] && IsOk; ++ii){
						Index = ii - 1;
// 						OutFile << "\t";
						ZoneStrs[SelZoneNum] += "\t";
						for (int Dir = 0; Dir < 3; ++Dir)
							ZoneStrs[SelZoneNum] += to_string_with_precision(XYZPtrs[SelZoneNum][Dir][Index], NumSigFigs) + " ";
// 						OutFile << endl;
						ZoneStrs[SelZoneNum] += "\n";
					}
				}


				ZoneStrs[SelZoneNum] += string("End\n\n")

					+ string("$CALCVARS\n\n")

					+ string("eor\n\n\n");

// 				OutFile << "End" << endl << endl
// 
// 					<< "$CALCVARS" << endl << endl
// 
// 					<< "eor" << endl << endl << endl;

// 				OutFile.flush();
			}
		}
	}

	TecUtilDataLoadEnd();

	StatusStr = "Writing densf script";

	for (int SelZoneNum = 0; SelZoneNum < NumZonesSelected && IsOk; SelZoneNum++){
		IsOk = StatusUpdate(SelZoneNum, NumZonesSelected, StatusStr, AddOnID);
		if (IsOk){
			OutFile << ZoneStrs[SelZoneNum];
		}
	}

	OutFile << "rm ./logfile ./SINFO* ./IINFO*" << endl << endl;

	OutFile << "echo \"Finished\"";

	OutFile.close();

	dos2unix2(OutFileName);
	//dos2unix(OutFileName);

	StatusDrop(AddOnID);

	if (IsOk)
		TecUtilDialogMessageBox("Finished", MessageBoxType_Information);
}


vector<string> SelectDensfTape41Folder(){
	/*
	*	Have user select a folder
	*/
	Boolean_t IsOk = TRUE;

	/*
	*	Get contents of folder and populate file list with it
	*/
	vector<string> T41FileNames;
	StringList_pa StrList;
	if (IsOk){

		IsOk = TecUtilDialogGetFileNames(SelectFileOption_AllowMultiFileRead, &StrList, "Tape41 file(s)", nullptr, "*.t41");
	}

	if (IsOk){
		T41FileNames.resize(TecUtilStringListGetCount(StrList));
		for (int i = 1; i <= T41FileNames.size(); ++i)
			T41FileNames[i-1] = TecUtilStringListGetString(StrList, i);
	}

	return T41FileNames;
}

void ImportAdditionalTape41Files(Boolean_t MatchZones, Boolean_t MatchDataSet){
	Boolean_t IsOk = TRUE;

	/*
	*	Get list of files to import
	*/

	vector<string> T41FileNames = SelectDensfTape41Folder();

	char *TmpName;
	int NumFiles = static_cast<int>(T41FileNames.size());
	IsOk = (NumFiles > 0);

	/*
	*	Get dataset info from Tecplot
	*/
	string DataSetName;
	EntIndex_t NumZones, NumVars;

	if (IsOk){
		IsOk = TecUtilDataSetGetInfo(&TmpName, &NumZones, &NumVars);

		if (IsOk)
			DataSetName = TmpName;

		TecUtilStringDealloc(&TmpName);
	}

#if defined(WIN32) || defined(_WIN32) || defined(__WIN32) && !defined(__CYGWIN__)
	string Slash = "\\";
#else
	string Slash = "/";
#endif

	if (IsOk){
		vector<string> T41Sec = {
			"SCF",
			"SumFrag",
			"Ortho",
			"Core"
		};
		vector<string> T41Var = {
			"Density",
			"Fitdensity",
			"Kinetic Energy Density",
			"CoulPot",
			"XCPot",
			"ELF",
			"DensityGradX",
			"DensityGradY",
			"DensityGradZ",
			"DensityHessXX",
			"DensityHessXY",
			"DensityHessXZ",
			"DensityHessYY",
			"DensityHessYZ",
			"DensityHessZZ",
			"DensityLap"
		};
		/*
		*	And the strings that will be used to name the variables
		*/
		vector<string> VarStrings = {
			"Electron Density",
			"Electron Density Fit",
			"Kinetic Energy Density",
			"Coulomb Potential",
			"Exchange Correlation Potential",
			"Electron Localization Function",
			"X Density Gradient",
			"Y Density Gradient",
			"Z Density Gradient",
			"XX Density Hessian",
			"XY Density Hessian",
			"XZ Density Hessian",
			"YY Density Hessian",
			"YZ Density Hessian",
			"ZZ Density Hessian",
			"Density Laplacian"
		};
		vector<string> T41VarSuffices = { "", "_A", "_B" };

		/*
		*	Transition density is a little different
		*/
		string T41TransDensSec = "TransDens";
		string TransDensString = "Transition";
		vector<string> T41TransDensSecSuf = { "", "_L1", "_L2" };
		vector<string> T41TransDensVar = { "Fitdensity", "Coulpot" };
		vector<string> TransDensVarStrs = { "Electron Density Fit", "Coulomb Potential" };
		vector<string> T41TransDensVarSuf = { "", "_1", "_2", "_3" };


		vector<int> T41DataSizes;
		vector<string> T41LoadKFVarStrings, T41LoadUserVarStrings;
		string FileName, FullPath, ZoneName, InDataSetName;
		EntIndex_t ZoneNum;
		int MaxIJK[3], InNumPoints = 0;
		vector<FieldDataType_e> DataTypes(NumZones, FieldDataType_Float);

		TecUtilDataLoadBegin();

		string const StatusStr = "Importing data";
		StatusLaunch(StatusStr, AddOnID, TRUE);

		for (int FileNum = 0; FileNum < NumFiles; ++FileNum){
			IsOk = StatusUpdate(FileNum, NumFiles, StatusStr.c_str(), AddOnID);

			if (IsOk){
				FullPath = T41FileNames[FileNum];
				FileName = FullPath.substr(FullPath.find_last_of(Slash) + 1, FullPath.length());

				InDataSetName = FileName.substr(0, FileName.find(DS_ZoneName_Delim));
				ZoneName = FileName.substr(FileName.find(DS_ZoneName_Delim) + DS_ZoneName_Delim.length(), FileName.find_last_of(".") - (FileName.find(DS_ZoneName_Delim) + DS_ZoneName_Delim.length()));

				/*
				 *	Zones that have colon ":" in the name have the colon replaced with forward slash "/"
				 *	when the output Tape41 file is made, which is read by Windows as a middle dot, which
				 *	is read by CPP as "" for some reason. These characters are ASCII -17, -128, and -94
				 *	respectively, so replace those in zone name with colon if present.
				 */

				ZoneName = StringReplaceSubString(ZoneName, string({ (char)-17, (char)-128, (char)-94 }), ":");

				if (MatchZones){

					ZoneNum = ZoneNumByName(ZoneName, false, true);

					IsOk = (ZoneNum > 0 && InDataSetName == DataSetName);
					REQUIRE(ZoneNum > 0 && InDataSetName == DataSetName);
				}
				else{
					ZoneName = FileName.substr(0, FileName.find_last_of(".") - 1);
				}
			}

			if (IsOk)
				TecUtilZoneGetIJK(ZoneNum, &MaxIJK[0], &MaxIJK[1], &MaxIJK[2]);

			/*
			*	Open the tape41 file
			*/
			KFFile T41File;
			if (IsOk){
				IsOk = (openKFFile(&T41File, const_cast<char*>(FullPath.c_str())) > 0);
			}

			if (IsOk)
				IsOk = (getKFData(&T41File, "Grid%total nr of points", &InNumPoints) > 0);
			if (IsOk){
				if (TecUtilZoneIsOrdered(ZoneNum))
					IsOk = (InNumPoints == MaxIJK[0] * MaxIJK[1] * MaxIJK[2]);
				else
					IsOk = (InNumPoints == MaxIJK[0]);
			}

			/*
			*	Now get list of all the data from the Tape41.
			*	If a variable is in the Tape41 that doesn't exist in the
			*	tecplot file, then make it (this only applies to the first file,
			*	since they all contain the same variables)
			*/
			/*
			*	Query the Tape41 file for all possible variables (that might
			*	be of interest to a Bondalyzer user)
			*
			*	First for the "main" variables (everything except transition values)
			*/
			T41DataSizes.clear();
			T41LoadKFVarStrings.clear();
			T41LoadUserVarStrings.clear();
			if (IsOk) for (string const & Sec : T41Sec){
				unsigned int VarNum = 0;
				for (string const & Var : T41Var){
					for (string const & Suffix : T41VarSuffices){
						stringstream ss1;
						ss1 << Sec << "%" << Var << Suffix;
						int DataSize = getKFVariableLength(&T41File, ss1.str().c_str());
						if (DataSize > 0){
							stringstream ss2;
							ss2 << VarStrings[VarNum] << Suffix << " " << Sec;
							T41LoadKFVarStrings.push_back(ss1.str());
							T41LoadUserVarStrings.push_back(ss2.str());
							T41DataSizes.push_back(DataSize);
						}
					}
					++VarNum;
				}
			}

			/*
			*	Then for the transition variables
			*/
			unsigned int VarNum = 0;
			for (string const & Var : T41TransDensVar){
				for (string const & Suffix : T41TransDensVarSuf){
					stringstream ss1;
					ss1 << T41TransDensSec << "%" << Var << Suffix;
					int DataSize = getKFVariableLength(&T41File, ss1.str().c_str());
					if (DataSize > 0){
						stringstream ss2;
						ss2 << TransDensString << " " << TransDensVarStrs[VarNum] << Suffix;
						T41LoadKFVarStrings.push_back(ss1.str());
						T41LoadUserVarStrings.push_back(ss2.str());
						T41DataSizes.push_back(DataSize);
					}
				}
				++VarNum;
			}

			/*
			*	Now we know all the variables that are in the file.
			*	Loop over the list, adding each variable.
			*	If the variable already exists, then overwrite the
			*	zone's values for that variable.
			*/
			for (int i = 0; i < T41DataSizes.size() && IsOk; ++i){
				int VarNum = VarNumByName(T41LoadUserVarStrings[i]);
				if (VarNum <= 0){
					TecUtilDataSetAddVar(T41LoadUserVarStrings[i].c_str(), DataTypes.data());
					VarNum = VarNumByName(T41LoadUserVarStrings[i]);
				}
				FieldData_pa VarRef = TecUtilDataValueGetWritableNativeRef(ZoneNum, VarNum);
				IsOk = (VALID_REF(VarRef));
				if (IsOk){
					vector<double> TempBuffer(T41DataSizes[i]);
					IsOk = (getKFData(&T41File, T41LoadKFVarStrings[i].c_str(), TempBuffer.data()) > 0);
					if (IsOk){
						if (TecUtilDataValueGetType(ZoneNum, VarNum) == FieldDataType_Double)
							TecUtilDataValueArraySetByRef(VarRef, 1, T41DataSizes[i], TempBuffer.data());
						else{
// #pragma warning( push )
// #pragma warning( disable : 4244 )
							vector<float> FloatBuffer(TempBuffer.cbegin(), TempBuffer.cend());
// #pragma warning( pop )
							TecUtilDataValueArraySetByRef(VarRef, 1, T41DataSizes[i], FloatBuffer.data());
						}
					}
				}
			}

			if (IsOk)
				closeKFFile(&T41File);
		}

		TecUtilDataLoadEnd();

		StatusDrop(AddOnID);

		if (IsOk)
			TecUtilDialogMessageBox("Finished", MessageBoxType_Information);
	}
}


Boolean_t LoadADFTape21(){
	Boolean_t IsOk = TRUE;

	TecUtilLockStart(AddOnID);

	if (TecUtilDataSetIsAvailable()){
		TecUtilDialogErrMsg("Cannot load file into existing dataset. Start with new file and try again.");
		TecUtilLockFinish(AddOnID);
		return FALSE;
	}

	/*
	 *	Get the Tape21 file path
	 */
	char* FileName = nullptr;
	IsOk = TecUtilDialogGetFileName(SelectFileOption_ReadSingleFile, &FileName, "ADF Tape21", nullptr, "*.t21");

	/*
	 *	Open the Tape21 file
	 */
	if (IsOk && FileName != nullptr){
		FileNameStr = FileName;
		IsOk = (openKFFile(&GlobalKFFile, const_cast<char*>(FileNameStr.c_str())) > 0);
	}

	if (!IsOk){
		TecUtilDialogErrMsg("Failed to open Tape21 file.");
		return LoaderFailed;
	}

	Boolean_t DeleteUselessData = TRUE;

	/*
	 *	Get all information about contents of Tape21.
	 *	Includes names of sections/variables and length/type of variables
	 */
	vector<string> SectionNames = GetSectionNameList(&GlobalKFFile);

	vector<vector<string> > VariableNames(SectionNames.size());
	vector<vector<int> > TypeList(SectionNames.size()), LengthList(SectionNames.size());

	for (int i = 0; i < SectionNames.size() && IsOk; ++i)
		IsOk = GetVariableInfoForSection(&GlobalKFFile, SectionNames[i].c_str(), VariableNames[i], TypeList[i], LengthList[i]);

	for (string & i : SectionNames)
		i = i.substr(0, i.find_last_not_of(' ') + 1);

	/*
	 *	Get the base name of the Tape21 file for nameing the dataset
	 */
#if defined(WIN32) || defined(_WIN32) || defined(__WIN32) && !defined(__CYGWIN__)
	string Slash = "\\";
#else
	string Slash = "/";
#endif

	string DataSetName;

	size_t StrPos = FileNameStr.find_last_of("/\\");
	if (StrPos != string::npos){
		DataSetName = FileNameStr.substr(StrPos + 1, FileNameStr.length());
	}
	StrPos = DataSetName.find_last_of(".");
	if (StrPos != string::npos){
		DataSetName = DataSetName.substr(0, StrPos);
	}


	/*
	 *	Prepare to make the dataset
	 */

	vector<string> XYZVarNames = { "X", "Y", "Z", T21Prefix + "ImportData" };
	StringList_pa VarNames = TecUtilStringListAlloc();
	for (string & i : XYZVarNames) TecUtilStringListAppendString(VarNames, i.c_str());

	if (IsOk)
		IsOk = TecUtilDataSetCreate(StringMakeValidZoneName(DataSetName).c_str(), VarNames, TRUE);

	TecUtilStringListDealloc(&VarNames);

	/*
	*	Need to make XYZ variables in order to place the system.
	*/

	vector<EntIndex_t> XYZVarNums = { 1, 2, 3 };

	/*
	 *	Now loop over all variables.
	 *	Those containing bool/int/real are saved as zones, while
	 *	those containing chars are added as AuxData.
	 */

	EntIndex_t ZoneNum, ImportDataVarNum = TecUtilDataSetGetNumVars();

	Set_pa ZoneSet = TecUtilSetAlloc(TRUE);

	TecUtilDataLoadBegin();

	for (int KFSecNum = 0; KFSecNum < SectionNames.size() && IsOk; ++KFSecNum){
		string SecName = SectionNames[KFSecNum];

		for (int KFVarNum = 0; KFVarNum < VariableNames[KFSecNum].size() && IsOk; ++KFVarNum){
			string VarName = VariableNames[KFSecNum][KFVarNum];
			int VarType = TypeList[KFSecNum][KFVarNum];
			int VarLen = LengthList[KFSecNum][KFVarNum];

			string FullName = SecName + KFDelim + VarName;

#ifdef _DEBUG
			/*
			*	Confirm that the data in the lists we already made agree with
			*	the KF file (still)
			*/
			IsOk = (VarLen == getKFVariableLength(&GlobalKFFile, FullName.c_str()));
			if (IsOk)
				IsOk = (VarType == getKFVariableType(&GlobalKFFile, FullName.c_str()));
#endif // _DEBUG

			if (IsOk){
				if (VarType == KF_T_STRING){
					/*
					 *	Add as aux data to dataset
					 */
					vector<char> TmpChars(VarLen);
					IsOk = (getKFData(&GlobalKFFile, FullName.c_str(), TmpChars.data()) > 0);
					if (IsOk){
						string AuxDataName = T21Prefix + FullName;
						IsOk = AuxDataDataSetSetItem(AuxDataName, TmpChars.data());
					}
				}
				else if (VarLen > 0){
					vector<FieldDataType_e> TypeVec(4, FieldDataType_Bit);
					TypeVec[3] = KF_TP_Types[VarType - 1];

					IsOk = TecUtilDataSetAddZone(StringMakeValidZoneName(T21Prefix + FullName).c_str(), VarLen, 1, 1, ZoneType_Ordered, TypeVec.data());
					ZoneNum = TecUtilDataSetGetNumZones();

					FieldData_pa DataRef = nullptr;
					if (IsOk){
						TecUtilSetAddMember(ZoneSet, ZoneNum, TRUE);
						DataRef = TecUtilDataValueGetWritableNativeRef(ZoneNum, ImportDataVarNum);
						IsOk = VALID_REF(DataRef);
					}

					if (IsOk && VarType == KF_T_DOUBLE){
						vector<double> TmpBuffer(VarLen);
						IsOk = (getKFData(&GlobalKFFile, FullName.c_str(), TmpBuffer.data()) > 0);
						if (IsOk){
							TecUtilDataValueArraySetByRef(DataRef, 1, VarLen, TmpBuffer.data());
						}
					}
					else if (IsOk){
						/*
						 *	the "bool" type in KF files is actually 4 bytes, same as the int type
						 */
						vector<int> TmpBuffer(VarLen);
						IsOk = (getKFData(&GlobalKFFile, FullName.c_str(), TmpBuffer.data()) > 0);
						if (IsOk){
							TecUtilDataValueArraySetByRef(DataRef, 1, VarLen, TmpBuffer.data());
						}
					}
				}
			}
		}
	}

	TecUtilDataLoadEnd();

	TecUtilStateChanged(StateChange_ZonesAdded, reinterpret_cast<ArbParam_t>(ZoneSet));

	TecUtilZoneSetActive(ZoneSet, AssignOp_MinusEquals);
	TecUtilZoneSetContour(SV_SHOW, ZoneSet, 0.0, FALSE);
	TecUtilZoneSetEdgeLayer(SV_SHOW, ZoneSet, 0.0, FALSE);
	TecUtilZoneSetScatter(SV_SHOW, ZoneSet, 0.0, FALSE);
	TecUtilZoneSetShade(SV_SHOW, ZoneSet, 0.0, FALSE);
	TecUtilZoneSetVector(SV_SHOW, ZoneSet, 0.0, FALSE);
	TecUtilZoneSetMesh(SV_SHOW, ZoneSet, 0.0, FALSE);

	TecUtilSetDealloc(&ZoneSet);

	closeKFFile(&GlobalKFFile);



	/*
	 *	Now all Tape21 contents are in the tecplot file.
	 */

	TecUtilFrameSetPlotType(PlotType_Cartesian3D);

	/*
	 *	Load atoms
	 */

	vector<AtomColor_s> AtomColorList;
	PopulateAtomColorList(AtomColorList);
	if (IsOk){
		string AtomTypes = AuxDataDataSetGetItem(T21Prefix + "Geometry%atomtype");
		int NumAtomTypes = AtomTypes.length() / 160;

		vector<AtomGroup_s> AtomGroupList(NumAtomTypes);

		FieldDataPointer_c CumAtomNumPtr, AtomXYZPtr, AtomChgPtr, AtomOrderIndex;

		TecUtilDataLoadBegin();

		int Tmp;

		Tmp = ZoneNumByName(T21Prefix + "Geometry%cum nr of atoms");
		if (Tmp > 0)
			IsOk = CumAtomNumPtr.InitializeReadPtr(Tmp, ImportDataVarNum);
		else IsOk = FALSE;

		Tmp = ZoneNumByName(T21Prefix + "Geometry%xyz");
		if (Tmp > 0)
			IsOk = AtomXYZPtr.InitializeReadPtr(Tmp, ImportDataVarNum);
		else IsOk = FALSE;

		Tmp = ZoneNumByName(T21Prefix + "Geometry%atomtype total charge");
		if (Tmp > 0)
			IsOk = AtomChgPtr.InitializeReadPtr(Tmp, ImportDataVarNum);
		else IsOk = FALSE;

		Tmp = ZoneNumByName(T21Prefix + "Geometry%atom order index");
		if (Tmp > 0)
			AtomOrderIndex.InitializeReadPtr(Tmp, ImportDataVarNum);
		else IsOk = FALSE;

		int NumAtoms;

		Tmp = ZoneNumByName(T21Prefix + "Geometry%nr of atoms");
		if (Tmp > 0)
			NumAtoms = TecUtilDataValueGetByZoneVar(Tmp, ImportDataVarNum, 1);
		else IsOk = FALSE;

		int iNum = 0;
		for (int i = 0; i < NumAtomTypes && IsOk; ++i){
			string TmpStr = AtomTypes.substr(i * 160, (i + 1) * 160);
			TmpStr = TmpStr.substr(0, TmpStr.find_first_of(' '));
			AtomGroupList[i].Name = TmpStr;
			AtomGroupList[i].Count = CumAtomNumPtr[i + 1] - CumAtomNumPtr[i];
			AtomGroupList[i].Charges = vector<ImportType_t>(AtomGroupList[i].Count, AtomChgPtr[i]);

			for (int j = 0; j < AtomGroupList[i].Count; ++j){
				vec3 Pos;
				for (int k = 0; k < 3; ++k){
					Pos[k] = AtomXYZPtr[iNum];
					++iNum;
				}
				AtomGroupList[i].AddPosition(Pos);
			}

			AtomGroupList[i].AtomColor = GetAtomColor(AtomColorList, AtomGroupList[i].Name);
		}

		vector<string> InputAtomLabels(NumAtoms);
		vector<vec3> InputAtomCoords(NumAtoms);

		for (int AtomNum = 0; AtomNum < NumAtoms; ++AtomNum){
			for (int j = 1; j <= NumAtomTypes; ++j){
				if (AtomOrderIndex[AtomNum] > CumAtomNumPtr[j - 1] && AtomOrderIndex[AtomNum] <= CumAtomNumPtr[j]){
					InputAtomLabels[AtomNum] = AtomGroupList[j - 1].Name + to_string(AtomNum + 1);
					for (int Dir = 0; Dir < 3; ++Dir){
						InputAtomCoords[AtomNum][Dir] = AtomXYZPtr[Dir + (AtomNum) * 3];
					}
				}
			}
		}

		

		/*
		 *	Add the atoms as zones for each type
		 */

		vector<FieldDataType_e> TypeVec(4, FieldDataType_Float);
		TypeVec[3] = FieldDataType_Bit;

		ZoneSet = TecUtilSetAlloc(TRUE);

		for (AtomGroup_s & A : AtomGroupList){
			IsOk = TecUtilDataSetAddZone(StringMakeValidZoneName(A.Name).c_str(), A.Count, 1, 1, ZoneType_Ordered, TypeVec.data());
			FieldData_pa DataRef;
			if (IsOk){
				ZoneNum = TecUtilDataSetGetNumZones();
				for (int i = 0; i < 3 && IsOk; ++i){
					DataRef = TecUtilDataValueGetWritableNativeRef(ZoneNum, XYZVarNums[i]);
					if (VALID_REF(DataRef)){
						TecUtilDataValueArraySetByRef(DataRef, 1, A.Count, A.Positions[i].data());
					}
					else IsOk = FALSE;
				}
			}

			IsOk = TecUtilSetAddMember(ZoneSet, ZoneNum, TRUE);

			Set_pa TempSet = TecUtilSetAlloc(TRUE);
			TecUtilSetAddMember(TempSet, ZoneNum, TRUE);

			TecUtilZoneSetScatter(SV_SHOW, TempSet, 0.0, TRUE);
			TecUtilZoneSetScatter(SV_COLOR, TempSet, 0.0, A.AtomColor);
			TecUtilZoneSetScatter(SV_FRAMESIZE, TempSet, MIN(2 * sqrt(A.Charges[0]), 8), 0);
			TecUtilZoneSetScatterSymbolShape(SV_GEOMSHAPE, TempSet, GeomShape_Sphere);

			TecUtilZoneSetMesh(SV_SHOW, TempSet, 0.0, FALSE);
			TecUtilZoneSetContour(SV_SHOW, TempSet, 0.0, FALSE);
			TecUtilZoneSetVector(SV_SHOW, TempSet, 0.0, FALSE);
			TecUtilZoneSetShade(SV_SHOW, TempSet, 0.0, FALSE);
			TecUtilZoneSetEdgeLayer(SV_SHOW, TempSet, 0.0, FALSE);

			TecUtilSetDealloc(&TempSet);
		}

		/*
		 *	That's it for atoms. Now the gui bonds.
		 */

		Set_pa GuiBondZones = TecUtilSetAlloc(TRUE);
		FieldDataPointer_c BondTypePtr, BondAtomNumPtr, AtomXYZInputOrderPtr;
		/*
		 *	The atom coordinates that correspond to the atom indices specified by
		 *	the BondAtomNumPtr are in AtomXYZInputOrderPtr.
		 *	
		 *	Bond 1's atoms' indices are number 1 and number 1 + NumBonds in
		 *	BondAtomNumPtr.
		 *	
		 *	The correspondance between the two different sets
		 *	of atom indices is in AtomOrderIndex.
		 *	Atom 1 is also atom (1 + NumAtoms) in AtomOrderIndex.
		 */

		int NumBonds;
		int NumPtsPerBond = 500;

		Tmp = ZoneNumByName(T21Prefix + "Geometry%guinbonds");
		if (Tmp > 0){
			NumBonds = TecUtilDataValueGetByZoneVar(Tmp, ImportDataVarNum, 1);
			IsOk = (NumBonds > 0);
		}
		else IsOk = FALSE;

		Tmp = ZoneNumByName(T21Prefix + "Geometry%guibondatoms");
		if (Tmp > 0)
			BondAtomNumPtr.InitializeReadPtr(Tmp, ImportDataVarNum);
		else IsOk = FALSE;

		Tmp = ZoneNumByName(T21Prefix + "Geometry%guibondtypes");
		if (Tmp > 0)
			BondTypePtr.InitializeReadPtr(Tmp, ImportDataVarNum);
		else IsOk = FALSE;

		Tmp = ZoneNumByName(T21Prefix + "Geometry%xyz InputOrder");
		if (Tmp > 0)
			AtomXYZInputOrderPtr.InitializeReadPtr(Tmp, ImportDataVarNum);
		else IsOk = FALSE;

		/*
		 *	Make scatter zone for bond midpoints
		 */
		Set_pa GUIBondMidpointZoneSet = TecUtilSetAlloc(FALSE);
		IsOk = TecUtilDataSetAddZone("GUI bond midpoints", NumBonds, 1, 1, ZoneType_Ordered, TypeVec.data());
		ZoneNum = TecUtilDataSetGetNumZones();
		IsOk = TecUtilSetAddMember(GUIBondMidpointZoneSet, ZoneNum, TRUE);
		TecUtilZoneSetMesh(SV_SHOW, GUIBondMidpointZoneSet, 0.0, FALSE);
		TecUtilZoneSetScatterSymbolShape(SV_GEOMSHAPE, GUIBondMidpointZoneSet, GeomShape_Sphere);
		TecUtilZoneSetScatter(SV_FRAMESIZE, GUIBondMidpointZoneSet, 1, 0);
		TecUtilZoneSetScatter(SV_SHOW, GUIBondMidpointZoneSet, 0.0, TRUE);
		TecUtilZoneSetActive(GUIBondMidpointZoneSet, AssignOp_MinusEquals);
		FieldVecPointer_c GUIBondMidpointZoneXYZPtr;
		GUIBondMidpointZoneXYZPtr.InitializeWritePtr(ZoneNum, { 1,2,3 });

		for (int BNum = 0; BNum < NumBonds; ++BNum){
			string BondStr = "GUI bond " + to_string(BNum + 1) + ": ";
			for (int i = 0; i < 2; ++i){
				int AtomNum = BondAtomNumPtr[BNum + (NumBonds * i)];
				for (int j = 1; j <= NumAtomTypes; ++j){
					if (AtomOrderIndex[AtomNum-1] > CumAtomNumPtr[j - 1] && AtomOrderIndex[AtomNum-1] <= CumAtomNumPtr[j]){
						BondStr += AtomGroupList[j - 1].Name + to_string(AtomNum);
						if (i == 0){
							double BondType = BondTypePtr[BNum];
							for (int k = 0; k < BondTypePtr[BNum]; ++k){
								BondStr += "-";
							}
						}
						break;
					}
				}
			}

			IsOk = TecUtilDataSetAddZone(StringMakeValidZoneName(BondStr).c_str(), NumPtsPerBond, 1, 1, ZoneType_Ordered, TypeVec.data());
			if (IsOk){
				ZoneNum = TecUtilDataSetGetNumZones();
				IsOk = TecUtilSetAddMember(ZoneSet, ZoneNum, TRUE);
			}
			if (IsOk)
				IsOk = TecUtilSetAddMember(GuiBondZones, ZoneNum, TRUE);

			vector<vec3> StartEndXYZ(2);

			for (int i = 0; i < 2 && IsOk; ++i){
				for (int j = 0; j < 3 && IsOk; ++j){
					int AtomInd = j + (3 * (BondAtomNumPtr[BNum + i * NumBonds] - 1));
					double TmpVal = AtomXYZInputOrderPtr[AtomInd];
					StartEndXYZ[i][j] = TmpVal;
// 					IsOk = TecUtilDataValueSetByZoneVar(ZoneNum, 
// 						XYZVarNums[j], 
// 						i, 
// 						TmpVal);
				}
			}

			GUIBondMidpointZoneXYZPtr.Write(BNum, (StartEndXYZ[0] + StartEndXYZ[1]) * 0.5);

			vec3 stepVec = (StartEndXYZ[1] - StartEndXYZ[0]) / double(NumPtsPerBond - 1);
			FieldVecPointer_c xyzVecPtr;
			xyzVecPtr.InitializeWritePtr(ZoneNum, { 1,2,3 });
			for (int i = 0; i < NumPtsPerBond; ++i){
				xyzVecPtr.Write(i, StartEndXYZ[0] + stepVec * double(i));
			}

			Set_pa TempSet = TecUtilSetAlloc(TRUE);
			if (IsOk){
				TecUtilSetAddMember(TempSet, ZoneNum, TRUE);

				TecUtilZoneSetMesh(SV_SHOW, TempSet, 0.0, TRUE);
				TecUtilZoneSetMesh(SV_COLOR, TempSet, 0.0, Black_C);
				TecUtilZoneSetScatterSymbolShape(SV_GEOMSHAPE, TempSet, GeomShape_Sphere);
				TecUtilZoneSetScatter(SV_FRAMESIZE, TempSet, 0.5, 0);
				TecUtilZoneSetScatter(SV_SHOW, TempSet, 0.0, FALSE);
				if (BondTypePtr[BNum] == static_cast<int>(BondTypePtr[BNum]))
					TecUtilZoneSetMesh(SV_LINEPATTERN, TempSet, 0.0, LinePattern_Solid);
				else
					TecUtilZoneSetMesh(SV_LINEPATTERN, TempSet, 0.0, LinePattern_LongDash);
				TecUtilZoneSetMesh(SV_LINETHICKNESS, TempSet, BondTypePtr[BNum] * 0.5, 0);

				TecUtilZoneSetScatter(SV_SHOW, TempSet, 0.0, FALSE);
				TecUtilZoneSetContour(SV_SHOW, TempSet, 0.0, FALSE);
				TecUtilZoneSetVector(SV_SHOW, TempSet, 0.0, FALSE);
				TecUtilZoneSetShade(SV_SHOW, TempSet, 0.0, FALSE);
				TecUtilZoneSetEdgeLayer(SV_SHOW, TempSet, 0.0, FALSE);
			}

			TecUtilSetDealloc(&TempSet);
		}

		TecUtilZoneSetActive(ZoneSet, AssignOp_PlusEquals);
		TecUtilStateChanged(StateChange_ZonesAdded, reinterpret_cast<ArbParam_t>(ZoneSet));
		TecUtilSetDealloc(&ZoneSet);


		/*
		 *	gui bonds are created. 
		 *	Now do Bader CPs and bond paths if they exist
		 */

		Tmp = ZoneNumByName(T21Prefix + "Properties%CP number of");
		Set_pa BondPathZones;

		if (IsOk && Tmp > 0){
			/*
			 *	CPs were calculated.
			 */
			ZoneSet = TecUtilSetAlloc(TRUE);

			FieldDataPointer_c CPCoordPtr, CPRankPtr, CPChgPtr, BPNumStepsPtr, BPPropPtr;
			/*
			 *	CP coordinates are stored so the the x coord of all 
			 *	CPs comes first, then the y coord of all...
			 *	The "type" numbers go as {1,2,3,4} = {n,c,b,r}
			 *	
			 *	bond paths are stored in a weird way in BPPropPtr:
			 *	the first value is the x coord for the first point of the
			 *	first bond path. 
			 *	Then follows the x coord for the first point of the 
			 *	remaining paths.
			 *	Then the y coord for the first point of all the paths
			 *	Then the z coord, density, grad XYZ, hess XX XY XZ YY YZ ZZ,
			 *	then the laplacian, I think.
			 *	
			 *	Then this repeats for the second point of all the paths, then
			 *	the third...
			 *	
			 *	So those 13 pieces of information.
			 */

			/*
			 *	First the critical points
			 */

			int NumCPs = TecUtilDataValueGetByZoneVar(Tmp, ImportDataVarNum, 1);
			
			if (NumCPs > 0){
				Tmp = ZoneNumByName(T21Prefix + "Properties%CP coordinates");
				if (IsOk && Tmp > 0)
					CPCoordPtr.InitializeReadPtr(Tmp, ImportDataVarNum);
				else IsOk = FALSE;

				Tmp = ZoneNumByName(T21Prefix + "Properties%CP code number for (Rank,Signatu");
				if (IsOk && Tmp > 0)
					CPRankPtr.InitializeReadPtr(Tmp, ImportDataVarNum);
				else IsOk = FALSE;

				Tmp = ZoneNumByName(T21Prefix + "Properties%CP density at");
				if (IsOk && Tmp > 0)
					CPChgPtr.InitializeReadPtr(Tmp, ImportDataVarNum);
				else IsOk = FALSE;
			}


			if (IsOk && NumCPs > 0){
				vector<string> CPTypeStrs = {
					"Nuclear CPs",
					"Cage CPs",
					"Bond CPs",
					"Ring CPs"
				};

				AtomGroupList.clear();
				AtomGroupList.reserve(4);
				for (string const & i : CPTypeStrs) AtomGroupList.push_back(i);
				for (AtomGroup_s & i : AtomGroupList) i.AtomColor = GetAtomColor(AtomColorList, i.Name);

				for (int i = 0; i < NumCPs; ++i){
					vec3 Pos;
					for (int Dir = 0; Dir < 3; ++Dir)
						Pos[Dir] = CPCoordPtr[i + Dir * NumCPs];

					AtomGroupList[CPRankPtr[i]-1].AddPosition(Pos);
					AtomGroupList[CPRankPtr[i]-1].Charges.push_back(CPChgPtr[i]);
				}

				for (AtomGroup_s & A : AtomGroupList){
					A.Count = A.Positions->size();

					if (A.Count > 0){
						IsOk = TecUtilDataSetAddZone(StringMakeValidZoneName(A.Name).c_str(), A.Count, 1, 1, ZoneType_Ordered, TypeVec.data());
						FieldData_pa DataRef;
						if (IsOk){
							ZoneNum = TecUtilDataSetGetNumZones();
							for (int i = 0; i < 3 && IsOk; ++i){
								DataRef = TecUtilDataValueGetWritableNativeRef(ZoneNum, XYZVarNums[i]);
								if (VALID_REF(DataRef)){
									TecUtilDataValueArraySetByRef(DataRef, 1, A.Count, A.Positions[i].data());
								}
								else IsOk = FALSE;
							}
						}

						double AvgChg = 0.0;
						for (ImportType_t & i : A.Charges) AvgChg += i;
						AvgChg /= A.Count;

						if (A.Name != CPTypeStrs[0])
							IsOk = TecUtilSetAddMember(ZoneSet, ZoneNum, TRUE);

						if (IsOk){
							Set_pa TempSet = TecUtilSetAlloc(TRUE);
							TecUtilSetAddMember(TempSet, ZoneNum, TRUE);

							TecUtilZoneSetScatter(SV_SHOW, TempSet, 0.0, TRUE);
							TecUtilZoneSetScatter(SV_COLOR, TempSet, 0.0, A.AtomColor);
							TecUtilZoneSetScatter(SV_FRAMESIZE, TempSet, MAX(MIN(sqrt(AvgChg), 8), 1), 0);
							TecUtilZoneSetScatterSymbolShape(SV_GEOMSHAPE, TempSet, GeomShape_Sphere);

							TecUtilZoneSetMesh(SV_SHOW, TempSet, 0.0, FALSE);
							TecUtilZoneSetContour(SV_SHOW, TempSet, 0.0, FALSE);
							TecUtilZoneSetVector(SV_SHOW, TempSet, 0.0, FALSE);
							TecUtilZoneSetShade(SV_SHOW, TempSet, 0.0, FALSE);
							TecUtilZoneSetEdgeLayer(SV_SHOW, TempSet, 0.0, FALSE);

							TecUtilSetDealloc(&TempSet);
						}
					}
				}
			}

			/*
			 *	Now the bond paths
			 */

			Tmp = ZoneNumByName(T21Prefix + "Properties%BP number of");
			int NumBPs;
			if (Tmp > 0)
				NumBPs = TecUtilDataValueGetByZoneVar(Tmp, ImportDataVarNum, 1);
			else IsOk = FALSE;

			if (NumBPs > 0){
				Tmp = ZoneNumByName(T21Prefix + "Properties%BP step number");
				if (Tmp > 0)
					BPNumStepsPtr.InitializeReadPtr(Tmp, ImportDataVarNum);
				else IsOk = FALSE;

				Tmp = ZoneNumByName(T21Prefix + "Properties%BPs and their properties");
				if (Tmp > 0)
					BPPropPtr.InitializeReadPtr(Tmp, ImportDataVarNum);
				else IsOk = FALSE;
			}

			if (IsOk && NumBPs > 0){
				/*
				 *	Make a variable for the imported bond paths' density values
				 */
				vector<FieldDataType_e> TmpTypes(TecUtilDataSetGetNumZones(), FieldDataType_Bit);
				IsOk = TecUtilDataSetAddVar(string(T21Prefix + "ImportDensity").c_str(), TmpTypes.data());
			}

			int RhoVarNum = -1;
			if (IsOk){
				RhoVarNum = TecUtilDataSetGetNumVars();
				IsOk = (RhoVarNum > 0);
			}

			if (IsOk && NumBPs > 0){
				BondPathZones = TecUtilSetAlloc(TRUE);
				vector<GradPathBase_c> BPs(NumBPs);

				TypeVec.resize(RhoVarNum, FieldDataType_Float);
				TypeVec[RhoVarNum - 2] = FieldDataType_Bit;

				double MinMaxRho[2] = { 1e100, -1e100 };

				int MaxBPStep = 0;
				for (int i = 0; i < NumBPs; ++i) MaxBPStep = MAX(MaxBPStep, BPNumStepsPtr[i]);
				int GroupLen = MaxBPStep * NumBPs;

				for (int BPNum = 0; BPNum < NumBPs; ++BPNum){
					for (int StepNum = 0; StepNum < BPNumStepsPtr[BPNum]; ++StepNum){
						vec3 Pos;

						for (int Dir = 0; Dir < 3; ++Dir)
							Pos[Dir] = BPPropPtr[BPNum + Dir * GroupLen + StepNum * NumBPs];

						double Rho = BPPropPtr[BPNum + 3 * GroupLen + StepNum * NumBPs];

						BPs[BPNum].PointAppend(Pos, Rho);

						if (Rho > 0.0)
							MinMaxRho[0] = MIN(MinMaxRho[0], Rho);
						MinMaxRho[1] = MAX(MinMaxRho[1], Rho);
					}

					/*
					 *	Get name of bond path by find the Atoms at it's beginning and end.
					 */

					string BondStr = "Bond path " + to_string(BPNum + 1) + ": ";
					for (int i = 0; i < 2; ++i){
						double TmpDistSqr, MinDistSqr = 1e100;
						int MinInd;
						for (int AtomNum = 0; AtomNum < NumAtoms; ++AtomNum){
							TmpDistSqr = DistSqr(BPs[BPNum].XYZAt(i * (BPNumStepsPtr[BPNum] - 1)), InputAtomCoords[AtomNum]);
							if (TmpDistSqr < MinDistSqr){
								MinDistSqr = TmpDistSqr;
								MinInd = AtomNum;
							}
						}
						BondStr += InputAtomLabels[MinInd];
						if (i == 0)
							BondStr += "-";
					}

					ZoneNum = BPs[BPNum].SaveAsOrderedZone(BondStr, TypeVec, XYZVarNums, RhoVarNum);

					IsOk = TecUtilSetAddMember(BondPathZones, ZoneNum, TRUE);
				}

				if (IsOk){
					TecUtilZoneSetActive(GuiBondZones, AssignOp_MinusEquals);
					TecUtilZoneSetActive(BondPathZones, AssignOp_PlusEquals);

					TecUtilZoneSetMesh(SV_LINETHICKNESS, BondPathZones, 0.5, 0);
					TecUtilZoneSetMesh(SV_COLOR, BondPathZones, 0.0, MultiColor8_C);

					ArgList_pa ArgList = TecUtilArgListAlloc();

					TecUtilArgListAppendInt(ArgList, SV_CONTOURGROUP, 8);
					TecUtilArgListAppendInt(ArgList, SV_VAR, RhoVarNum);
					TecUtilContourSetVariableX(ArgList);

					TecUtilArgListClear(ArgList);

					vector<double> ContourLevels = LogLevels(MinMaxRho[0], MinMaxRho[1]);
					TecUtilArgListAppendInt(ArgList, SV_CONTOURLEVELACTION, ContourLevelAction_New);
					TecUtilArgListAppendInt(ArgList, SV_CONTOURGROUP, 8);
					TecUtilArgListAppendInt(ArgList, SV_NUMVALUES, ContourLevels.size());
					TecUtilArgListAppendArray(ArgList, SV_RAWDATA, ContourLevels.data());
					TecUtilContourLevelX(ArgList);

					TecUtilArgListClear(ArgList);

					TecUtilArgListAppendInt(ArgList, SV_CONTOURLABELACTION, ContourLabelAction_DeleteAll);
					TecUtilArgListAppendInt(ArgList, SV_CONTOURGROUP, 8);
					TecUtilContourLabelX(ArgList);

					TecUtilArgListClear(ArgList);

					TecUtilArgListAppendString(ArgList, SV_P1, SV_GLOBALCONTOUR);
					TecUtilArgListAppendString(ArgList, SV_P2, SV_LEGEND);
					TecUtilArgListAppendString(ArgList, SV_P3, SV_SHOW);
					TecUtilArgListAppendInt(ArgList, SV_OFFSET1, 8);
					TecUtilArgListAppendArbParam(ArgList, SV_IVALUE, FALSE);
					TecUtilStyleSetLowLevelX(ArgList);

					TecUtilArgListDealloc(&ArgList);
				}
				TecUtilSetDealloc(&BondPathZones);
			}

			TecUtilZoneSetActive(ZoneSet, AssignOp_PlusEquals);
			TecUtilStateChanged(StateChange_ZonesAdded, reinterpret_cast<ArbParam_t>(ZoneSet));
			TecUtilSetDealloc(&ZoneSet);
		}
		
		TecUtilSetDealloc(&GuiBondZones);

		ArgList_pa argList = TecUtilArgListAlloc();

		TecUtilArgListClear(argList);
		TecUtilArgListAppendString(argList, SV_P1, SV_ISOSURFACEATTRIBUTES);
		TecUtilArgListAppendString(argList, SV_P2, SV_EFFECTS);
		TecUtilArgListAppendString(argList, SV_P3, SV_USETRANSLUCENCY);
		TecUtilArgListAppendInt(argList, SV_OFFSET1, 1);
		TecUtilArgListAppendArbParam(argList, SV_IVALUE, TRUE);
		TecUtilStyleSetLowLevelX(argList);

		TecUtilArgListClear(argList);
		TecUtilArgListAppendString(argList, SV_P1, SV_ISOSURFACEATTRIBUTES);
		TecUtilArgListAppendString(argList, SV_P2, SV_EFFECTS);
		TecUtilArgListAppendString(argList, SV_P3, SV_SURFACETRANSLUCENCY);
		TecUtilArgListAppendInt(argList, SV_OFFSET1, 1);
		TecUtilArgListAppendArbParam(argList, SV_IVALUE, 70);
		TecUtilStyleSetLowLevelX(argList);

		TecUtilArgListClear(argList);
		TecUtilArgListAppendString(argList, SV_P1, SV_FIELDLAYERS);
		TecUtilArgListAppendString(argList, SV_P2, SV_SHOWSCATTER);
		TecUtilArgListAppendArbParam(argList, SV_IVALUE, TRUE);
		TecUtilStyleSetLowLevelX(argList);

		TecUtilArgListClear(argList);
		TecUtilArgListAppendString(argList, SV_P1, SV_FIELDLAYERS);
		TecUtilArgListAppendString(argList, SV_P2, SV_SHOWMESH);
		TecUtilArgListAppendArbParam(argList, SV_IVALUE, TRUE);
		TecUtilStyleSetLowLevelX(argList);

		TecUtilArgListClear(argList);
		TecUtilArgListAppendString(argList, SV_P1, SV_ISOSURFACELAYERS);
		TecUtilArgListAppendString(argList, SV_P2, SV_SHOW);
		TecUtilArgListAppendArbParam(argList, SV_IVALUE, FALSE);
		TecUtilStyleSetLowLevelX(argList);

		TecUtilArgListClear(argList);
		TecUtilArgListAppendString(argList, SV_P1, SV_ISOSURFACEATTRIBUTES);
		TecUtilArgListAppendString(argList, SV_P2, SV_DEFINITIONCONTOURGROUP);
		TecUtilArgListAppendInt(argList, SV_OFFSET1, 1);
		TecUtilArgListAppendArbParam(argList, SV_IVALUE, (SmInteger_t)1);
		TecUtilStyleSetLowLevelX(argList);

		TecUtilArgListClear(argList);
		TecUtilArgListAppendString(argList, SV_P1, SV_ISOSURFACEATTRIBUTES);
		TecUtilArgListAppendString(argList, SV_P2, SV_SHOWGROUP);
		TecUtilArgListAppendInt(argList, SV_OFFSET1, 1);
		TecUtilArgListAppendArbParam(argList, SV_IVALUE, TRUE);
		TecUtilStyleSetLowLevelX(argList);

		TecUtilArgListDealloc(&argList);
		
		TecUtilDataLoadEnd();
	}

	/*
	 *	Delete all the raw Tape21 data
	 */
	Set_pa tmpSet = TecUtilSetAlloc(FALSE);
	for (int z = 1; z <= TecUtilDataSetGetNumZones(); ++z){
		char * zoneName;
		if (TecUtilZoneGetName(z, &zoneName)){
			string nameStr(zoneName);
			if (nameStr.length() >= 4 && nameStr.find(T21Prefix) != string::npos) {
				TecUtilSetAddMember(tmpSet, z, FALSE);
			}
			TecUtilStringDealloc(&zoneName);
		}
	}
	if (!TecUtilSetIsEmpty(tmpSet)) {
		TecUtilDataSetDeleteZone(tmpSet);
	}
	TecUtilSetDealloc(&tmpSet);
	tmpSet = TecUtilSetAlloc(FALSE);
	TecUtilSetAddMember(tmpSet, 4, FALSE);
	TecUtilDataSetDeleteVar(tmpSet);
	TecUtilSetDealloc(&tmpSet);
	

	if (IsOk)
		TecUtilViewDataFit();

	
	TecUtilLockFinish(AddOnID);

	return IsOk;
}

/*
 *	Loads one or more FORMATTED gaussian cube files. Will fail if they're not formatted.
 */
void LoadGaussianCubeFiles()
{
	TecUtilLockStart(AddOnID);

	SYSTEM_INFO sysinfo;
	GetSystemInfo(&sysinfo);

	int numCPU = sysinfo.dwNumberOfProcessors;
	omp_set_num_threads(numCPU);

	StringList_pa FileNames = TecUtilStringListAlloc();
	Boolean_t IsOk = TecUtilDialogGetFileNames(SelectFileOption_AllowMultiFileRead, &FileNames, "Gaussian Formatted Cube File", nullptr, "*.cube");

	Boolean_t UseExistingVolZone = FALSE;
	FileNameStrs.clear();

	if (IsOk && FileNames != nullptr){
		FileNameStrs.resize(TecUtilStringListGetCount(FileNames));
		for (int i = 1; i <= FileNameStrs.size(); ++i){
			if (i == 1){
				if (TecUtilDataSetIsAvailable()){
					ReplaceDataSet = !TecUtilDialogMessageBox("There is already data loaded. Do you want to append to or replace the current data?\n\tYes: Append to current data set\n\tNo: Replace current data set", MessageBoxType_YesNo);
					if (!ReplaceDataSet){
						UseExistingVolZone = TecUtilDialogMessageBox("Import data to existing full volume zone or create new zone?\n\tYes: Import to current zone\n\tNo: Import to new zone", MessageBoxType_YesNo);
						if (!UseExistingVolZone)
							ImportAtoms = TecUtilDialogMessageBox("Import atomic positions from all files(s)?", MessageBox_YesNo);
						else
							ImportAtoms = FALSE;
					}
				}
				else{
					ReplaceDataSet = TRUE;
					ImportAtoms = TRUE;
				}
			}

			FileNameStrs[i - 1] = TecUtilStringListGetString(FileNames, i);
		}

		TecUtilStringListDealloc(&FileNames);

		string HeaderChkStr = " Title Card Required ";

		/*
		 *	Sorry for the rediculous logic below.
		 *	It has to do with the way the order of data in a gaussian cube file,
		 *	so stupid logic can be avoided during the data loading step.
		 *	Blame Gaussian!
		 */
		vector<string> VarChkStrs = {
				"Density",
				"density",
				"Gradient",
				"gradient",
				"Laplacian",
				"laplacian",
				"NormGradient",
				"normGradient",
				"Normgradient",
				"normgradient",
				"Potential"
		};
		vector<vector<string> > VarUserNames = {
			{"Electron Density"},
			{
				"X Density Gradient",
				"Y Density Gradient",
				"Z Density Gradient"
			},
			{"Density Laplacian"},
			{"Density Gradient Magnitude"},
			{"Electrostatic Potential"}
		};
		vector<vector<int> > VarNameIndices = {
			{0},
			{0},
			{0,1},
			{0,1},
			{0,1,2},
			{0,1,2},
			{5,1},
			{5,1},
			{5,1},
			{5,1},
			{6}
		};

		GroupIndex = 50;
		/*
		 *	Loop over each file selected and import the data
		 */
		for (string const & CubeFileName : FileNameStrs){
			//	Open the cube file
			ifstream CubeFile(CubeFileName);
			if (!CubeFile.is_open()){
				TecUtilDialogErrMsg(string("Failed to open " + CubeFileName).c_str());
				return;
			}

			/*
			 *	Parse the file header. Here's the structure:
			 *	
			 *	Header line 1: title and variable requested for the cube file
			 *	Header line 2: 
			 *	NAtoms, X-Origin, Y-Origin, Z-Origin NVal        NVal is the #points/value
			 *	N1, X1, Y1, Z1               # of increments in the slowest running direction
			 *	N2, X2, Y2, Z2
			 *	N3, X3, Y3, Z3               # of increments in the fastest running direction
			 *	IA1, Chg1, X1, Y1, Z1        Atomic number, charge, and coordinates of the first atom
			 *	�
			 *	IAn, Chgn, Xn, Yn, Zn        Atomic number, charge, and coordinates of the last atom
			 *	(N1*N2) records, each of length N3     Values of the density at each point in the grid
			 */

			/*
			 *	Get header lines
			 */

			string H1, H2;

			getline(CubeFile, H1);
			getline(CubeFile, H2);

			/*
			 *	Test that there actually was a header.
			 *	Just need to check that the first word of H1 is not an integer.
			 *	(This assumes an integer will never be the first word in a formatted
			 *	cube file, which may not be a bad assumption now and/or in the future)
			 */

			stringstream ss;
			ss << H1;
			string TmpStr;
			ss >> TmpStr;

			if (IsNumber(TmpStr)){
				TecUtilDialogErrMsg(string("Unformatted Gaussian cube files won't work. Make a formatted cube file and try again.\nFile: " + CubeFileName).c_str());
				CubeFile.close();
				return;
			}

			StatusLaunch("Reading file contents. Please wait.", AddOnID, FALSE);
			/*
			 *	trim off the "Title Card Required" from the front of the header if it's present
			 */
			if (H1.length() > HeaderChkStr.length() && H1.substr(0, HeaderChkStr.length()) == HeaderChkStr){
				H1 = H1.substr(HeaderChkStr.length(), H1.length() - HeaderChkStr.length());
			}
			else
				H1 = H1.substr(H1.find_last_of(" ") + 1, H1.length() - H1.find_last_of(" "));

			/*
			 *	Now get the rest of the header
			 */

			Boolean_t IsMOFile = FALSE;
			int NumAtoms;
			vec3 Origin;
			mat33 LatticeVector;
			vector<int> IJK(3);
			vector<AtomGroup_s> AtomInfo;

			/*
			 *	Line three of the header
			 */
			CubeFile >> NumAtoms;
			/*
			 *	If the number of atoms < 0, then the cube file contains MO information.
			 */
			if (NumAtoms < 0){
				IsMOFile = TRUE;
				NumAtoms = -NumAtoms;
			}
			for (int i = 0; i < 3; ++i)
				CubeFile >> Origin[i];

			/*
			 *	There could be an extra value that I don't care about at the end of this line, so skip it.
			 */
			getline(CubeFile, TmpStr);

			/*
			 *	Lines 4-6 of the header: number of points and lattice directions
			 */

			for (int i = 0; i < 3; ++i){
				CubeFile >> IJK[i];
				for (int j = 0; j < 3; ++j){
					CubeFile >> LatticeVector.at(i, j);
				}
			}

			/*
			 *	Lines 7-(7+NumAtoms) of header: atom info
			 */
			vector<AtomColor_s> AtomColorList;
			PopulateAtomColorList(AtomColorList);

			for (int i = 0; i < NumAtoms; ++i){
				int ElemNum;
				double AtomChg;
				vec3 AtomPos;

				CubeFile >> ElemNum >> AtomChg;
				for (int j = 0; j < 3; ++j)
					CubeFile >> AtomPos[j];

				string AtomSymbol = GetAtomStr(ElemNum);

				/*
				 *	See if we already have an AtomGroup_s for this atom type
				 */
				Boolean_t IsFound = FALSE;
				for (AtomGroup_s & a : AtomInfo){
					if (a.Name == AtomSymbol){
						IsFound = TRUE;

						a.AddPosition(AtomPos);
						a.Charges.push_back(AtomChg);
						a.Count++;

						break;
					}
				}

				if (!IsFound){
					AtomGroup_s NewAtomGroup(AtomSymbol);
					NewAtomGroup.AtomColor = GetAtomColor(AtomColorList, NewAtomGroup.Name);
					NewAtomGroup.AddPosition(AtomPos);
					NewAtomGroup.Charges.push_back(AtomChg);

					AtomInfo.push_back(NewAtomGroup);
				}
			}

			/*
			 *	If MOs in file, read num of MOs which MOs are present
			 */

			int NumMOs = -1;
			vector<int> MOList;
			if (IsMOFile){
				CubeFile >> NumMOs;
				MOList.resize(NumMOs);
				for (int i = 0; i < NumMOs; ++i)
					CubeFile >> MOList[i];
			}



			/*
			*	parse the header lines to get the relevant variable names
			*/
			string VarSuffix = "";
			int GaussVarNum = -1;
			int NumVarBlocks = -1;
			for (int i = 0; i < VarChkStrs.size(); ++i){
				if (H1.length() >= VarChkStrs[i].length() && H1.substr(0, VarChkStrs[i].length()) == VarChkStrs[i]){
					GaussVarNum = i;
					NumVarBlocks = VarNameIndices[i].size();
					break;
				}
			}
			if (GaussVarNum < 0){
				if (IsMOFile)
					NumVarBlocks = NumMOs;
				else
					NumVarBlocks = 1;
			}

			size_t TmpStrPos = H1.find_last_of("=");
			if (TmpStrPos != string::npos){
				VarSuffix = " " + H1.substr(TmpStrPos + 1, H1.length() - TmpStrPos - 1);
			}

			/*
			 *	Header is now read. Next fetch the first value in the file.
			 *	After reading the entire file into memory, we'll use the first
			 *	value to identify when we've arrived at the beginning of the 
			 *	data block.
			 */

			string FirstVal;
			CubeFile >> FirstVal;

			/*
			 *	Read whole file into memory
			 */
			CubeFile.seekg(0, CubeFile.end);
			std::streamoff FileSize = CubeFile.tellg();
			CubeFile.seekg(0, CubeFile.beg);


			vector<char> CubeFileContents(FileSize);
			CubeFile.read(CubeFileContents.data(), FileSize);
			CubeFile.close();

			/*
			 *	Run through the char vector CubeFileContents to find the first data value
			 */

			char* TmpCStr = strtok(CubeFileContents.data(), " \n\t");
			while (TmpCStr != nullptr && FirstVal.compare(TmpCStr))
				TmpCStr = strtok(nullptr, " \n\t");

			/*
			 *	Make dataset and volume zone if necessary
			 */
			StringList_pa VarNames = TecUtilStringListAlloc();
			vector<string> VarNameStrs = { "X", "Y", "Z" };

			for (int i = 0; i < NumVarBlocks; ++i){
				if (IsMOFile)
					VarNameStrs.push_back(H1 + " " + to_string(MOList[i]));
				else{
					for (string const & Str : VarUserNames[VarNameIndices[GaussVarNum][i]])
						VarNameStrs.push_back(Str + VarSuffix);
				}
			}
			for (string const & i : VarNameStrs)
				TecUtilStringListAppendString(VarNames, i.c_str());

#if defined(WIN32) || defined(_WIN32) || defined(__WIN32) && !defined(__CYGWIN__)
			string Slash = "\\";
#else
			string Slash = "/";
#endif

			string DataSetName = CubeFileName.substr(0, CubeFileName.find_last_of("."));
			TmpStrPos = DataSetName.find_last_of(Slash.c_str());
			if (TmpStrPos != string::npos){
				DataSetName = DataSetName.substr(TmpStrPos + 1, DataSetName.length() - TmpStrPos - 1);
			}
			vector<FieldDataType_e> VarDataTypes;
			EntIndex_t ZoneNum, VolZoneNum = -1;
			if (IsOk){
				if (ReplaceDataSet){
					VarDataTypes.resize(VarNameStrs.size(), FieldDataType_Float);
					IsOk = TecUtilDataSetCreate(StringMakeValidZoneName(DataSetName).c_str(), VarNames, TRUE);
				}
				else{
					if (!TecUtilDataSetIsAvailable()){
						VarDataTypes.resize(VarNameStrs.size(), FieldDataType_Float);
						IsOk = TecUtilDataSetCreate(StringMakeValidZoneName(DataSetName).c_str(), VarNames, TRUE);
					}
					else{
						EntIndex_t NumVars = TecUtilDataSetGetNumVars();
						VarDataTypes.reserve(NumVars);
						for (int i = 1; i <= NumVars; ++i){
							char* TmpCStr;
							TecUtilVarGetName(i, &TmpCStr);
							int Ind = SearchVectorForString(VarNameStrs, TmpCStr, true);
							TecUtilStringDealloc(&TmpCStr);
							if (Ind >= 0)
								VarDataTypes.push_back(FieldDataType_Float);
							else
								VarDataTypes.push_back(FieldDataType_Bit);
						}
					}
				}
				if (IsOk){
					if (UseExistingVolZone){
						VolZoneNum = ZoneNumByName("Full Volume");
					}
					if (VolZoneNum <= 0){
						IsOk = TecUtilDataSetAddZone(StringMakeValidZoneName(CSMZoneName.FullVolume + DataSetName).c_str(), IJK[0], IJK[1], IJK[2], ZoneType_Ordered, VarDataTypes.data());
						if (IsOk){
							VolZoneNum = TecUtilDataSetGetNumZones();
							IsOk = VolZoneNum > 0;
						}
					}
				}
			}

			vector<FieldDataType_e> ZoneDataTypes(TecUtilDataSetGetNumZones(), FieldDataType_Bit);
			ZoneDataTypes[VolZoneNum - 1] = FieldDataType_Float;

			StatusDrop(AddOnID);

			TecUtilDataLoadBegin();

			Set_pa VolZoneSet = TecUtilSetAlloc(TRUE);
			if (IsOk){
				SetZoneStyle({ VolZoneNum }, ZoneStyle_Volume);
				TecUtilSetAddMember(VolZoneSet, VolZoneNum, TRUE);
			}

			/*
			 *	Make new variables where necessary and get pointers to variable data
			 */
			int VarNum;
			vector<FieldDataPointer_c> Ptrs(VarNameStrs.size());
			for (int i = 0; i < VarNameStrs.size() && IsOk; ++i){
				if (!ReplaceDataSet){
					VarNum = VarNumByName(VarNameStrs[i]);
					if (VarNum <= 0){
						IsOk = TecUtilDataSetAddVar(VarNameStrs[i].c_str(), ZoneDataTypes.data());
						if (IsOk){
							VarNum = TecUtilDataSetGetNumVars();
							IsOk = (VarNum > 0);
						}
					}
				}
				else
					VarNum = i + 1;

				if (IsOk)
					IsOk = Ptrs[i].InitializeWritePtr(VolZoneNum, VarNum);
			}

			
			/*
			 *	Now just need to parse through the char array and populate the necessary variables
			 */
			int TotalPoints = IJK[0] * IJK[1] * IJK[2] * NumVarBlocks;
			string StatusStr = "Loading Gaussian Cube file: " + DataSetName;
			StatusLaunch(StatusStr, AddOnID, TRUE);

			double Val;
			vec3 TmpIJK = zeros<vec>(3);
			for (int i = 0; i < IJK[0] && IsOk; ++i){
				IsOk = StatusUpdate(i, IJK[0], StatusStr, AddOnID);
				TmpIJK[1] = 0;
				for (int j = 0; j < IJK[1] && IsOk; ++j){
					if (IsMOFile){
						TmpIJK[2] = 0;
						for (int k = 0; k < IJK[2]; ++k){
							int Index = IndexFromIJK(i + 1, j + 1, k + 1, IJK[0], IJK[1]) - 1;
							for (int Dir = 0; Dir < 3; ++Dir){
								Val = Origin[Dir] + dot(TmpIJK, LatticeVector.col(Dir));
								Ptrs[Dir].Write(Index, Val);
							}
							for (VarNum = 3; VarNum < NumVarBlocks + 3 && TmpCStr != nullptr; ++VarNum){
								Val = atof(TmpCStr);
								Ptrs[VarNum].Write(Index, Val);
								TmpCStr = strtok(nullptr, " \n\t");
							}
							++TmpIJK[2];
						}
					}
					else{
						VarNum = 3;
						for (int iVar = 0; iVar < VarNameIndices[GaussVarNum].size(); ++iVar){
							TmpIJK[2] = 0;
							for (int k = 0; k < IJK[2]; ++k){
								int Index = IndexFromIJK(i + 1, j + 1, k + 1, IJK[0], IJK[1]) - 1;
								if (iVar == 0){
									for (int Dir = 0; Dir < 3; ++Dir){
										Val = Origin[Dir] + dot(TmpIJK, LatticeVector.col(Dir));
										Ptrs[Dir].Write(Index, Val);
									}
								}

								for (int jVar = 0; jVar < VarUserNames[VarNameIndices[GaussVarNum][iVar]].size() && TmpCStr != nullptr; ++jVar){
									Val = atof(TmpCStr);
									Ptrs[VarNum + jVar].Write(Index, Val);
									TmpCStr = strtok(nullptr, " \n\t");
								}
								++TmpIJK[2];
							}

							VarNum += VarUserNames[VarNameIndices[GaussVarNum][iVar]].size();
						}
					}
					++TmpIJK[1];
				}
				++TmpIJK[0];
			}

			StatusDrop(AddOnID);

			if (TmpCStr == nullptr 
				&& (TmpIJK[0] < IJK[0] || TmpIJK[1] < IJK[1] || TmpIJK[2] < IJK[2]))
			{
				TecUtilDialogErrMsg(string("Failed to load Gaussian file: " + CubeFileName).c_str());
				TecUtilDataLoadEnd();
				return;
			}


			for (FieldDataPointer_c & i : Ptrs) i.Close();

			TecUtilDataLoadEnd();

			if (!ReplaceDataSet){
				TecUtilZoneSetScatter(SV_SHOW, VolZoneSet, 0.0, FALSE);
				TecUtilZoneSetMesh(SV_SHOW, VolZoneSet, 0.0, FALSE);
				TecUtilZoneSetContour(SV_SHOW, VolZoneSet, 0.0, FALSE);
				TecUtilZoneSetEdgeLayer(SV_SHOW, VolZoneSet, 0.0, TRUE);
			}

			//	Modify group number of volume zone
			if (IsOk){
				ArgList_pa ArgList = TecUtilArgListAlloc();
				TecUtilArgListAppendString(ArgList, SV_P1, SV_FIELDMAP);
				TecUtilArgListAppendString(ArgList, SV_P2, SV_GROUP);
				TecUtilArgListAppendSet(ArgList, SV_OBJECTSET, VolZoneSet);
				TecUtilArgListAppendArbParam(ArgList, SV_IVALUE, (LgIndex_t)(GroupIndex));
				TecUtilStyleSetLowLevelX(ArgList);
				TecUtilArgListDealloc(&ArgList);
			}

			TecUtilSetDealloc(&VolZoneSet);

			/*
			 *	Now create zones for atomss
			 */

			if (IsOk && ReplaceDataSet || ImportAtoms){
				IsOk = CreateAtomZonesFromAtomGroupList(AtomInfo, { "X", "Y", "Z" }, VarDataTypes, GroupIndex);
			}

			/*
			 *	Modify zone styles and such
			 */
			if (IsOk){
				//	Modify group number of volume zone
				ArgList_pa ArgList = TecUtilArgListAlloc();

				TecUtilStateChanged(StateChange_ZonesAdded, reinterpret_cast<ArbParam_t>(VolZoneSet));

				TecUtilFrameSetPlotType(PlotType_Cartesian3D);

				TecUtilZoneSetActive(VolZoneSet, AssignOp_PlusEquals);

				TecUtilSetDealloc(&VolZoneSet);


				TecUtilArgListClear(ArgList);
				TecUtilArgListAppendString(ArgList, SV_P1, SV_ISOSURFACEATTRIBUTES);
				TecUtilArgListAppendString(ArgList, SV_P2, SV_EFFECTS);
				TecUtilArgListAppendString(ArgList, SV_P3, SV_USETRANSLUCENCY);
				TecUtilArgListAppendInt(ArgList, SV_OFFSET1, 1);
				TecUtilArgListAppendArbParam(ArgList, SV_IVALUE, TRUE);
				TecUtilStyleSetLowLevelX(ArgList);

				TecUtilArgListClear(ArgList);
				TecUtilArgListAppendString(ArgList, SV_P1, SV_ISOSURFACEATTRIBUTES);
				TecUtilArgListAppendString(ArgList, SV_P2, SV_EFFECTS);
				TecUtilArgListAppendString(ArgList, SV_P3, SV_SURFACETRANSLUCENCY);
				TecUtilArgListAppendInt(ArgList, SV_OFFSET1, 1);
				TecUtilArgListAppendArbParam(ArgList, SV_IVALUE, 70);
				TecUtilStyleSetLowLevelX(ArgList);

				TecUtilArgListClear(ArgList);
				TecUtilArgListAppendString(ArgList, SV_P1, SV_FIELDLAYERS);
				TecUtilArgListAppendString(ArgList, SV_P2, SV_SHOWSCATTER);
				TecUtilArgListAppendArbParam(ArgList, SV_IVALUE, TRUE);
				TecUtilStyleSetLowLevelX(ArgList);

				TecUtilArgListClear(ArgList);
				TecUtilArgListAppendString(ArgList, SV_P1, SV_ISOSURFACELAYERS);
				TecUtilArgListAppendString(ArgList, SV_P2, SV_SHOW);
				TecUtilArgListAppendArbParam(ArgList, SV_IVALUE, TRUE);
				TecUtilStyleSetLowLevelX(ArgList);

				TecUtilArgListClear(ArgList);
				TecUtilArgListAppendString(ArgList, SV_P1, SV_ISOSURFACEATTRIBUTES);
				TecUtilArgListAppendString(ArgList, SV_P2, SV_DEFINITIONCONTOURGROUP);
				TecUtilArgListAppendInt(ArgList, SV_OFFSET1, 1);
				TecUtilArgListAppendArbParam(ArgList, SV_IVALUE, (SmInteger_t)1);
				TecUtilStyleSetLowLevelX(ArgList);

				TecUtilArgListClear(ArgList);
				TecUtilArgListAppendString(ArgList, SV_P1, SV_ISOSURFACEATTRIBUTES);
				TecUtilArgListAppendString(ArgList, SV_P2, SV_SHOWGROUP);
				TecUtilArgListAppendInt(ArgList, SV_OFFSET1, 1);
				TecUtilArgListAppendArbParam(ArgList, SV_IVALUE, TRUE);
				TecUtilStyleSetLowLevelX(ArgList);

				EntIndex_t IsoSurfaceVarNum = 4;
				// 	double ValMax, ValMin;
				// 	TecUtilVarGetMinMax(IsoSurfaceVarNum, &ValMin, &ValMax);
				// 	double IsoValue = ValMax - 0.05 * (ValMax - ValMin);
				double IsoValue = 0.25;

				TecUtilArgListClear(ArgList);
				TecUtilArgListAppendString(ArgList, SV_P1, SV_ISOSURFACEATTRIBUTES);
				TecUtilArgListAppendString(ArgList, SV_P2, SV_ISOVALUE1);
				TecUtilArgListAppendInt(ArgList, SV_OFFSET1, 1);
				TecUtilArgListAppendDouble(ArgList, SV_DVALUE, IsoValue);
				TecUtilStyleSetLowLevelX(ArgList);


				TecUtilArgListClear(ArgList);
				TecUtilArgListAppendString(ArgList, SV_P1, SV_GLOBALCONTOUR);
				TecUtilArgListAppendString(ArgList, SV_P2, SV_VAR);
				TecUtilArgListAppendInt(ArgList, SV_OFFSET1, 1);
				TecUtilArgListAppendArbParam(ArgList, SV_IVALUE, IsoSurfaceVarNum);
				TecUtilStyleSetLowLevelX(ArgList);
				TecUtilArgListClear(ArgList);
				TecUtilArgListAppendString(ArgList, SV_P1, SV_ISOSURFACEATTRIBUTES);
				TecUtilArgListAppendString(ArgList, SV_P2, SV_ISOVALUE1);
				TecUtilArgListAppendInt(ArgList, SV_OFFSET1, 1);
				TecUtilArgListAppendDouble(ArgList, SV_DVALUE, IsoValue);
				TecUtilStyleSetLowLevelX(ArgList);
				TecUtilArgListClear(ArgList);
				TecUtilArgListAppendString(ArgList, SV_P1, SV_ISOSURFACEATTRIBUTES);
				TecUtilArgListAppendString(ArgList, SV_P2, SV_ISOVALUE2);
				TecUtilArgListAppendInt(ArgList, SV_OFFSET1, 1);
				TecUtilArgListAppendDouble(ArgList, SV_DVALUE, IsoValue);
				TecUtilStyleSetLowLevelX(ArgList);
				TecUtilArgListClear(ArgList);
				TecUtilArgListAppendString(ArgList, SV_P1, SV_ISOSURFACEATTRIBUTES);
				TecUtilArgListAppendString(ArgList, SV_P2, SV_ISOVALUE3);
				TecUtilArgListAppendInt(ArgList, SV_OFFSET1, 1);
				TecUtilArgListAppendDouble(ArgList, SV_DVALUE, IsoValue);
				TecUtilStyleSetLowLevelX(ArgList);

				TecUtilArgListDealloc(&ArgList);

				TecUtilViewDataFit();
			}

			AuxDataDataSetSetItem(DLProgramName, "Gaussian 09");

			ReplaceDataSet = FALSE;
			GroupIndex++;
		}
	}

	TecUtilLockFinish(AddOnID);
}


bool GetFLAPWCHARGEFileHeader(string const & fileName, string & fieldName, int & nx, int & ny, int & ns, float & lx, float & ly, string & firstVal, int *numCharsHeader = nullptr, unsigned long long int *numCharsValBlock = nullptr, unsigned long long int *numCharsBlockTotal = nullptr, int *numPlanes = nullptr) {
	bool IsOk = true;
	bool computeBlockLength = (numCharsHeader != nullptr);
	int headerNumChars = 0;

	ifstream CHARGEFile(fileName);

	if (!CHARGEFile.is_open()) {
		return false;
	}

	string H1, H2;
	getline(CHARGEFile, H1);
	getline(CHARGEFile, H2);

	headerNumChars += H1.length() + H2.length();

	stringstream ss;
	ss << H1;

	fieldName = StringJoin(SplitString(H2, " ", true, true));


	// Split first line and check as valid header (+int +int +int +double +double)
	vector<string> SplitH1 = SplitString(H1, " ", true, true);


	if (SplitH1.size() != 5) {
		return false;
	}

	IsOk = (
		StringIsInt(SplitH1[0]) &&
		StringIsInt(SplitH1[1]) &&
		StringIsInt(SplitH1[2]) &&
		StringIsFloat(SplitH1[3]) &&
		StringIsFloat(SplitH1[4])
		);

	if (!IsOk) {
		return false;
	}

	ss >> nx >> ny >> ns >> lx >> ly; // read in first header line

	/*
	*	Header is now read. Next fetch the first value in the file.
	*	After reading the entire file into memory, we'll use the first
	*	value to identify when we've arrived at the beginning of the
	*	data block.
	*/

	CHARGEFile >> firstVal;

	while (!StringIsFloat(firstVal)){
		CHARGEFile >> firstVal;
	}

	if (computeBlockLength){
		// Determine number of planes using size of single block and total file size
		int nVals = nx * ny; 	// total values for a plane
		int nLinesPerBlock = nVals / 5;
		int nLinesRem = nVals % 5;

		// (12 chars * 5 vals) + (2 spaces * 5) + (1 newline) = 60 + 8 + 1 = 69
		int nCharsPerBlock = 71 * nLinesPerBlock;

		nCharsPerBlock += (12 * nLinesRem) + (2 * (nLinesRem)) + 1;

		*numCharsValBlock = nCharsPerBlock;

		headerNumChars += 1; // a newline is missed doing the two readlines

		nCharsPerBlock += headerNumChars;

		nCharsPerBlock += 1; // newline to next block

		//	Get file size
		CHARGEFile.seekg(0, CHARGEFile.end);
		std::streamoff FileSize = CHARGEFile.tellg();
		CHARGEFile.seekg(0, CHARGEFile.beg);

		*numPlanes = FileSize / nCharsPerBlock;
		if (FileSize % *numPlanes != 0){
			TecUtilDialogErrMsg("Filesize not as expected. Contact support.");
			return false;
		}
		*numCharsHeader = headerNumChars;
		*numCharsBlockTotal = nCharsPerBlock;
	}

	if (!StringIsFloat(firstVal)){
		return false;
	}

	return true;
}



bool GetFLAPWMultiCHARGEFileFieldData(string const & fileName, string const & firstVal, vector<float> & data)
{
	bool IsOk = true;

	ifstream CHARGEFile(fileName);

	if (!CHARGEFile.is_open()) {
		return false;
	}

	/*
	*	Read whole file into memory
	*/
	CHARGEFile.seekg(0, CHARGEFile.end);
	std::streamoff FileSize = CHARGEFile.tellg();
	CHARGEFile.seekg(0, CHARGEFile.beg);


	vector<char> CHARGEFileContents(FileSize);
	CHARGEFile.read(CHARGEFileContents.data(), FileSize);
	CHARGEFile.close();

	/*
	*	Run through the char vector CubeFileContents to find the first data value
	*/

	char* TmpCStr = strtok(CHARGEFileContents.data(), " \n\t");
	while (TmpCStr != nullptr && firstVal.compare(TmpCStr))
		TmpCStr = strtok(nullptr, " \n\t");

	// size of data tells us how many values are expected
	int i = 0;
	while (TmpCStr != nullptr && i < data.size()){
		data[i] = atof(TmpCStr);
		TmpCStr = strtok(nullptr, " \n\t");
		i++;
	}

	if (i < data.size()){
		return false;
	}

	return true;
}

bool GetFLAPWCHARGEFileFieldData(ifstream & CHARGEFile, unsigned long long int pos, vector<char> & CharArray, vector<float> & data)
{
	bool IsOk = true;

	if (!CHARGEFile.is_open()) {
		return false;
	}

	// Position stream to specified position
	CHARGEFile.seekg(pos);
	if (!CHARGEFile.good()) {
		bool isbad = CHARGEFile.bad();
		bool isend = CHARGEFile.eof();
		bool isfail = CHARGEFile.fail();
		return false;
	}
	
	// Find first numeral
	int c = CHARGEFile.peek();
	while (!std::isdigit(c)){
		CHARGEFile.get();
		c = CHARGEFile.peek();
	}

	CHARGEFile.read(CharArray.data(), CharArray.size());
	if (!CHARGEFile.good()){
		bool isbad = CHARGEFile.bad();
		bool isend = CHARGEFile.eof();
		bool isfail = CHARGEFile.fail();
		return false;
	}

	/*
	*	Run through the char vector CubeFileContents to find the first data value
	*/

	char* TmpCStr = strtok(CharArray.data(), " \n\t");
// 	while (TmpCStr != nullptr && firstVal.compare(TmpCStr))
// 		TmpCStr = strtok(nullptr, " \n\t");

	// size of data tells us how many values are expected
	int i = 0;
	while (TmpCStr != nullptr && i < data.size()) {
		data[i] = atof(TmpCStr);
		TmpCStr = strtok(nullptr, " \n\t");
		i++;
	}

	if (i < data.size()) {
		return false;
	}

	return true;
}

// Trump is a super callous fragile racist sexist nazi POTUS


/*
*	Loads one or more 2D CHARGE grids from FLAPW.
*	A 3d grid is composed of multiple 2d grid files.
*	Files must be sorted so that they correspond to consecutive slices of the 3d grid.
*	Header information from first file is used to determine the system extent.
*	Here we assume that the z spacing is the same as the x spacing
*
*	Example CHARGE contents
*
  801 1201    1  20.00000  30.00000 // NumXPts, NumYPts, NumSpins, X extent (bohr), Y extent (bohr)
 total density // name of density
  0.493235E-01  0.490900E-01  0.488877E-01  0.487172E-01  0.485794E-01
  0.484748E-01  0.484044E-01  0.483691E-01  0.483697E-01  0.484074E-01
  0.484833E-01  0.485986E-01  0.487545E-01  0.489524E-01  0.491937E-01
  0.494799E-01  0.498123E-01  0.501927E-01  0.506226E-01  0.511036E-01
  0.516374E-01  0.522257E-01  0.530421E-01  0.537439E-01  0.545290E-01
  0.553822E-01  0.563362E-01  0.573829E-01  0.585322E-01  0.598104E-01
  0.611943E-01  0.627353E-01  0.644193E-01  0.662755E-01  0.683283E-01
*/
void LoadFLAPWMultiCHARGEFiles() {
	StringList_pa FileNameStrList = TecUtilStringListAlloc();
	vector<string> FileNames;
	Boolean_t IsOk = TecUtilDialogGetFileNames(SelectFileOption_AllowMultiFileRead, &FileNameStrList, "FLAPW CHARGE File(s)", nullptr, "*");

	if (IsOk && FileNameStrList != nullptr){
		for (int i = 1; i <= TecUtilStringListGetCount(FileNameStrList); ++i){
			FileNames.emplace_back(TecUtilStringListGetString(FileNameStrList, i));
		}
	}

	if (!IsOk){
		return;
	}

// 	FileNames = {
// 		"Z:\\Tecplot\\StorageWorkspace\\2020_FLAPW-dataLoader\\CHARGE"
// 	};
// 	for (int i = 0; i < 10; ++i) FileNames.push_back(FileNames.front());

	/*
	 * First validata all files and verify they're all for the same system
	 */

	int nx, ny, nz, ns;
	float lx, ly, dx, dy;
	string fieldName;
	vector<string> firstVal(FileNames.size());

	bool filesValid = GetFLAPWCHARGEFileHeader(FileNames.front(), fieldName, nx, ny, ns, lx, ly, firstVal[0]);

	if (!filesValid) {
		TecUtilDialogErrMsg(string("File 1 header invalid (\'" + FileNames[0] + "\')").c_str());
		return;
	}

	for (int i = 1; i < FileNames.size() && filesValid; ++i){
		int nx1, ny1, ns1;
		float lx1, ly1;
		string fieldName1;
		filesValid = GetFLAPWCHARGEFileHeader(FileNames[i], fieldName1, nx1, ny1, ns1, lx1, ly1, firstVal[i])
			&& nx == nx1 && ny == ny1 && ns == ns1 
			&& lx == lx1 && ly == ly1 && fieldName == fieldName1;
		if (!filesValid){
			TecUtilDialogErrMsg(string("File " + to_string(i+1) + " header invalid or mismatch to group (\'" + FileNames[i] + "\')").c_str());
		}
	}

	if (!filesValid){
		return;
	}

	/*
	 * prepare dataset
	 */

	nz = FileNames.size();
	int npts = nx * ny * nz;
	dx = lx / float(nx - 1);
	dy = ly / float(ny - 1);


	vector<int> IJK = { nx,ny,nz };

	/*
	*	Make dataset and volume zone if necessary
	*/
	StringList_pa VarNames = TecUtilStringListAlloc();
	vector<string> VarNameStrs = { "X", "Y", "Z", CSMVarName.Dens};

	for (string const & i : VarNameStrs)
		TecUtilStringListAppendString(VarNames, i.c_str());

#if defined(WIN32) || defined(_WIN32) || defined(__WIN32) && !defined(__CYGWIN__)
	string Slash = "\\";
#else
	string Slash = "/";
#endif

	string DataSetName = FileNames[0].substr(0, FileNames[0].find_last_of("."));
	size_t TmpStrPos = DataSetName.find_last_of(Slash.c_str());
	if (TmpStrPos != string::npos){
		DataSetName = DataSetName.substr(TmpStrPos + 1, DataSetName.length() - TmpStrPos - 1);
	}
	vector<FieldDataType_e> VarDataTypes;
	EntIndex_t ZoneNum, VolZoneNum = -1;
	if (IsOk){
		if (ReplaceDataSet){
			VarDataTypes.resize(VarNameStrs.size(), FieldDataType_Float);
			IsOk = TecUtilDataSetCreate(StringMakeValidZoneName(DataSetName).c_str(), VarNames, TRUE);
		}
		else{
			if (!TecUtilDataSetIsAvailable()){
				VarDataTypes.resize(VarNameStrs.size(), FieldDataType_Float);
				IsOk = TecUtilDataSetCreate(StringMakeValidZoneName(DataSetName).c_str(), VarNames, TRUE);
			}
			else{
				EntIndex_t NumVars = TecUtilDataSetGetNumVars();
				VarDataTypes.reserve(NumVars);
				for (int i = 1; i <= NumVars; ++i){
					char* TmpCStr;
					TecUtilVarGetName(i, &TmpCStr);
					int Ind = SearchVectorForString(VarNameStrs, TmpCStr, true);
					TecUtilStringDealloc(&TmpCStr);
					if (Ind >= 0)
						VarDataTypes.push_back(FieldDataType_Float);
					else
						VarDataTypes.push_back(FieldDataType_Bit);
				}
			}
		}
		//if (VolZoneNum <= 0) {
			IsOk = TecUtilDataSetAddZone(StringMakeValidZoneName(CSMZoneName.FullVolume + DataSetName).c_str(), nx, ny, nz, ZoneType_Ordered, VarDataTypes.data());
			if (IsOk) {
				VolZoneNum = TecUtilDataSetGetNumZones();
				IsOk = VolZoneNum > 0;
			}
		//}
	}

	vector<FieldDataType_e> ZoneDataTypes(TecUtilDataSetGetNumZones(), FieldDataType_Bit);
	ZoneDataTypes[VolZoneNum - 1] = FieldDataType_Float;

	Set_pa VolZoneSet = TecUtilSetAlloc(TRUE);
	if (IsOk){
		SetZoneStyle({ VolZoneNum }, ZoneStyle_Volume);
		TecUtilSetAddMember(VolZoneSet, VolZoneNum, TRUE);
	}

	/*
		*	Make new variables where necessary and get pointers to variable data
		*/
	vector<int> VarNums;
	int VarNum;
// 	vector<FieldDataPointer_c> Ptrs(VarNameStrs.size());
	for (int i = 0; i < VarNameStrs.size() && IsOk; ++i){
		if (!ReplaceDataSet){
			VarNums.push_back(VarNumByName(VarNameStrs[i]));
			if (VarNums.back() <= 0){
				IsOk = TecUtilDataSetAddVar(VarNameStrs[i].c_str(), ZoneDataTypes.data());
				if (IsOk){
					VarNums.back() = TecUtilDataSetGetNumVars();
					IsOk = (VarNum > 0);
				}
			}
		}
		else
			VarNums.push_back(i + 1);

// 		if (IsOk)
// 			IsOk = Ptrs[i].InitializeWritePtr(VolZoneNum, VarNum);
	}

	TecUtilDataLoadBegin(); // keep tecplot from deallocating pointers

	FieldVecPointer_c XYZPtr;
	XYZPtr.InitializeWritePtr(VolZoneNum, vector<int>(VarNums.begin(), VarNums.end() - 1));

// 	FieldDataPointer_c ValPtr;
// 	ValPtr.InitializeWritePtr(VolZoneNum, VarNums.back());
			
	/*
		*	Now just need to parse through the char array and populate the necessary variables
		*/
	int TotalPoints = IJK[0] * IJK[1] * IJK[2];
	mat33 LatticeVector = zeros<mat>(3,3);
	LatticeVector.at(0, 0) = dx;
	LatticeVector.at(1, 1) = dx;
	LatticeVector.at(2, 2) = dx;
	vec3 Origin = zeros<vec>(3);
	string StatusStr = "Loading FLAPW data: " + DataSetName;
	StatusLaunch(StatusStr, AddOnID, TRUE);

	/*
	 *	Loop over files again, reading in data from each.
	 *	rho data can be read into tecplot as one big read, 
	 *	and we'll populate the xyz values as we go.
	 */

	vector<float> SliceData(nx*ny);

	// write rho data as single write for each file
	FieldData_pa VarRef = TecUtilDataValueGetWritableNativeRef(VolZoneNum, VarNums.back());
	if (VALID_REF(VarRef)) {
		vec3 TmpIJK = zeros<vec>(3);
		for (int k = 0; k < IJK[2] && IsOk; ++k) {
// 			IsOk = StatusUpdate(k, IJK[2], StatusStr, AddOnID);
			// Get data from file
			if (GetFLAPWMultiCHARGEFileFieldData(FileNames[k], firstVal[k], SliceData)){
				TecUtilDataValueArraySetByRef(VarRef, IndexFromIJK(1, 1, k + 1, IJK[0], IJK[1]), SliceData.size(), SliceData.data());
			}

			TmpIJK[1] = 0;
			for (int j = 0; j < IJK[1] && IsOk; ++j) {
				IsOk = StatusUpdate(j + IJK[1] * k, IJK[1] * IJK[2], StatusStr, AddOnID);
				TmpIJK[0] = 0;
				for (int i = 0; i < IJK[0]; ++i) {
					int Index = IndexFromIJK(i + 1, j + 1, k + 1, IJK[0], IJK[1]) - 1;
					XYZPtr.Write(Index, Origin + LatticeVector * TmpIJK);

					++TmpIJK[0];
				}
				++TmpIJK[1];
			}
			++TmpIJK[2];
		}
	}

	// write xyz values
// 	StatusStr = "Loading FLAPW data (step 2 of 2): " + DataSetName;
// 	vec3 TmpIJK = zeros<vec>(3);
// 	for (int i = 0; i < IJK[0] && IsOk; ++i) {
// 		IsOk = StatusUpdate(i, IJK[0], StatusStr, AddOnID);
// 		TmpIJK[1] = 0;
// 		for (int j = 0; j < IJK[1] && IsOk; ++j) {
// 			TmpIJK[2] = 0;
// 			for (int k = 0; k < IJK[2]; ++k) {
// 				int Index = IndexFromIJK(i + 1, j + 1, k + 1, IJK[0], IJK[1]) - 1;
// 				XYZPtr.Write(Index, Origin + LatticeVector * TmpIJK);
// 				
// 				++TmpIJK[2];
// 			}
// 			++TmpIJK[1];
// 		}
// 		++TmpIJK[0];
// 	}

	

	StatusDrop(AddOnID);

	XYZPtr.Close();
// 	for (FieldDataPointer_c & i : Ptrs) i.Close();


	if (!ReplaceDataSet){
		TecUtilZoneSetScatter(SV_SHOW, VolZoneSet, 0.0, FALSE);
		TecUtilZoneSetMesh(SV_SHOW, VolZoneSet, 0.0, FALSE);
		TecUtilZoneSetContour(SV_SHOW, VolZoneSet, 0.0, FALSE);
		TecUtilZoneSetEdgeLayer(SV_SHOW, VolZoneSet, 0.0, TRUE);
	}

	//	Modify group number of volume zone
	if (IsOk){
		ArgList_pa ArgList = TecUtilArgListAlloc();
		TecUtilArgListAppendString(ArgList, SV_P1, SV_FIELDMAP);
		TecUtilArgListAppendString(ArgList, SV_P2, SV_GROUP);
		TecUtilArgListAppendSet(ArgList, SV_OBJECTSET, VolZoneSet);
		TecUtilArgListAppendArbParam(ArgList, SV_IVALUE, (LgIndex_t)(GroupIndex));
		TecUtilStyleSetLowLevelX(ArgList);
		TecUtilArgListDealloc(&ArgList);
	}

	if (nz > 1) {
		// 3d, so turn on isosurface and stick a slice at 1/2 z
		TecUtilFrameSetPlotType(PlotType_Cartesian3D);

		double zMin, zMax;
		TecUtilDataValueGetMinMaxByZoneVar(VolZoneNum, VarNums[2], &zMin, &zMax);

		TecUtilMacroExecuteCommand("$!SLICELAYERS SHOW = YES");
		TecUtilMacroExecuteCommand("$!SLICEATTRIBUTES 1  SLICESURFACE = ZPLANES");
		TecUtilMacroExecuteCommand(string("$!SLICEATTRIBUTES 1  PRIMARYPOSITION{ Z = " + to_string((zMax - zMin) * 0.5) + " }").c_str());
		TecUtilMacroExecuteCommand("$!SLICEATTRIBUTES 1  CONTOUR{ CONTOURTYPE = FLOOD }");
		TecUtilMacroExecuteCommand("$!SLICEATTRIBUTES 1  CONTOUR{ FLOODCOLORING = GROUP1 }");

		ArgList_pa argList = TecUtilArgListAlloc();

		TecUtilArgListClear(argList);
		TecUtilArgListAppendString(argList, SV_P1, SV_ISOSURFACEATTRIBUTES);
		TecUtilArgListAppendString(argList, SV_P2, SV_EFFECTS);
		TecUtilArgListAppendString(argList, SV_P3, SV_USETRANSLUCENCY);
		TecUtilArgListAppendInt(argList, SV_OFFSET1, 1);
		TecUtilArgListAppendArbParam(argList, SV_IVALUE, TRUE);
		TecUtilStyleSetLowLevelX(argList);

		TecUtilArgListClear(argList);
		TecUtilArgListAppendString(argList, SV_P1, SV_ISOSURFACEATTRIBUTES);
		TecUtilArgListAppendString(argList, SV_P2, SV_EFFECTS);
		TecUtilArgListAppendString(argList, SV_P3, SV_SURFACETRANSLUCENCY);
		TecUtilArgListAppendInt(argList, SV_OFFSET1, 1);
		TecUtilArgListAppendArbParam(argList, SV_IVALUE, 70);
		TecUtilStyleSetLowLevelX(argList);

		TecUtilArgListClear(argList);
		TecUtilArgListAppendString(argList, SV_P1, SV_FIELDLAYERS);
		TecUtilArgListAppendString(argList, SV_P2, SV_SHOWSCATTER);
		TecUtilArgListAppendArbParam(argList, SV_IVALUE, TRUE);
		TecUtilStyleSetLowLevelX(argList);

		TecUtilArgListClear(argList);
		TecUtilArgListAppendString(argList, SV_P1, SV_ISOSURFACELAYERS);
		TecUtilArgListAppendString(argList, SV_P2, SV_SHOW);
		TecUtilArgListAppendArbParam(argList, SV_IVALUE, TRUE);
		TecUtilStyleSetLowLevelX(argList);

		TecUtilArgListClear(argList);
		TecUtilArgListAppendString(argList, SV_P1, SV_ISOSURFACEATTRIBUTES);
		TecUtilArgListAppendString(argList, SV_P2, SV_DEFINITIONCONTOURGROUP);
		TecUtilArgListAppendInt(argList, SV_OFFSET1, 1);
		TecUtilArgListAppendArbParam(argList, SV_IVALUE, (SmInteger_t)1);
		TecUtilStyleSetLowLevelX(argList);

		TecUtilArgListClear(argList);
		TecUtilArgListAppendString(argList, SV_P1, SV_ISOSURFACEATTRIBUTES);
		TecUtilArgListAppendString(argList, SV_P2, SV_SHOWGROUP);
		TecUtilArgListAppendInt(argList, SV_OFFSET1, 1);
		TecUtilArgListAppendArbParam(argList, SV_IVALUE, TRUE);
		TecUtilStyleSetLowLevelX(argList);

		EntIndex_t IsoSurfaceVarNum = VarNums.back();
		// 	double ValMax, ValMin;
		// 	TecUtilVarGetMinMax(IsoSurfaceVarNum, &ValMin, &ValMax);
		// 	double IsoValue = ValMax - 0.05 * (ValMax - ValMin);
		double IsoValue = 0.25;

		TecUtilArgListClear(argList);
		TecUtilArgListAppendString(argList, SV_P1, SV_ISOSURFACEATTRIBUTES);
		TecUtilArgListAppendString(argList, SV_P2, SV_ISOVALUE1);
		TecUtilArgListAppendInt(argList, SV_OFFSET1, 1);
		TecUtilArgListAppendDouble(argList, SV_DVALUE, IsoValue);
		TecUtilStyleSetLowLevelX(argList);


		TecUtilArgListClear(argList);
		TecUtilArgListAppendString(argList, SV_P1, SV_GLOBALCONTOUR);
		TecUtilArgListAppendString(argList, SV_P2, SV_VAR);
		TecUtilArgListAppendInt(argList, SV_OFFSET1, 1);
		TecUtilArgListAppendArbParam(argList, SV_IVALUE, IsoSurfaceVarNum);
		TecUtilStyleSetLowLevelX(argList);
		TecUtilArgListClear(argList);
		TecUtilArgListAppendString(argList, SV_P1, SV_ISOSURFACEATTRIBUTES);
		TecUtilArgListAppendString(argList, SV_P2, SV_ISOVALUE1);
		TecUtilArgListAppendInt(argList, SV_OFFSET1, 1);
		TecUtilArgListAppendDouble(argList, SV_DVALUE, IsoValue);
		TecUtilStyleSetLowLevelX(argList);
		TecUtilArgListClear(argList);
		TecUtilArgListAppendString(argList, SV_P1, SV_ISOSURFACEATTRIBUTES);
		TecUtilArgListAppendString(argList, SV_P2, SV_ISOVALUE2);
		TecUtilArgListAppendInt(argList, SV_OFFSET1, 1);
		TecUtilArgListAppendDouble(argList, SV_DVALUE, IsoValue);
		TecUtilStyleSetLowLevelX(argList);
		TecUtilArgListClear(argList);
		TecUtilArgListAppendString(argList, SV_P1, SV_ISOSURFACEATTRIBUTES);
		TecUtilArgListAppendString(argList, SV_P2, SV_ISOVALUE3);
		TecUtilArgListAppendInt(argList, SV_OFFSET1, 1);
		TecUtilArgListAppendDouble(argList, SV_DVALUE, IsoValue);
		TecUtilStyleSetLowLevelX(argList);

		TecUtilArgListDealloc(&argList);
	}
	else{
		// only one slice, so turn set to 2d and turn on contours for first zone
		TecUtilFrameSetPlotType(PlotType_Cartesian2D);
		TecUtilMacroExecuteCommand("$!FIELDLAYERS SHOWCONTOUR = YES");
		TecUtilZoneSetContour(SV_SHOW, VolZoneSet, 0.0, TRUE);
	}

	// set contour levels
	double varMin, varMax;
	TecUtilDataValueGetMinMaxByRef(VarRef, &varMin, &varMax);
	vec ContourLevels;
	if (varMin > 0) {
		ContourLevels = LogSpace(varMin, sqrt(varMax), 20);
	}
	else {
		ContourLevels = linspace(varMin, varMax, 20);
	}

	ArgList_pa TempArgList = TecUtilArgListAlloc();
	TecUtilArgListAppendInt(TempArgList, SV_CONTOURGROUP, 1);
	TecUtilArgListAppendInt(TempArgList, SV_VAR, VarNums.back());
	TecUtilContourSetVariableX(TempArgList);
	TecUtilArgListClear(TempArgList);

	TecUtilArgListAppendInt(TempArgList, SV_CONTOURLABELACTION, ContourLabelAction_DeleteAll);
	TecUtilArgListAppendInt(TempArgList, SV_CONTOURGROUP, 1);
	TecUtilContourLabelX(TempArgList);
	TecUtilArgListClear(TempArgList);

	TecUtilArgListAppendInt(TempArgList, SV_CONTOURLEVELACTION, ContourLevelAction_New);
	TecUtilArgListAppendInt(TempArgList, SV_CONTOURGROUP, 1);
	TecUtilArgListAppendInt(TempArgList, SV_NUMVALUES, 20);
	TecUtilArgListAppendArray(TempArgList, SV_RAWDATA, ContourLevels.memptr());
	TecUtilContourLevelX(TempArgList);
	TecUtilArgListClear(TempArgList);

	TecUtilArgListAppendString(TempArgList, SV_P1, SV_GLOBALCONTOUR);
	TecUtilArgListAppendString(TempArgList, SV_P2, SV_LEGEND);
	TecUtilArgListAppendString(TempArgList, SV_P3, SV_SHOW);
	TecUtilArgListAppendInt(TempArgList, SV_OFFSET1, 1);
	TecUtilArgListAppendArbParam(TempArgList, SV_IVALUE, FALSE);
	TecUtilStyleSetLowLevelX(TempArgList);
	TecUtilArgListClear(TempArgList);


	TecUtilDataLoadEnd();

	TecUtilSetDealloc(&VolZoneSet);

	AuxDataDataSetSetItem(DLProgramName, "FLAPW");

	return;
}

/*
*	Loads a 3d CHARGE grids from FLAPW.
*	A 3d grid is composed of multiple 2d grid files, all slapped together into the same ASCII file.
*	Header information from first file section is used to determine the system extent.
*	Here we assume that the z spacing is the same as the x spacing but ask anywyas cause it's not always the case
*
*	Example CHARGE contents
*
  801 1201    1  20.00000  30.00000 // NumXPts, NumYPts, NumSpins, X extent (bohr), Y extent (bohr)
 total density // name of density
  0.493235E-01  0.490900E-01  0.488877E-01  0.487172E-01  0.485794E-01
  0.484748E-01  0.484044E-01  0.483691E-01  0.483697E-01  0.484074E-01
  0.484833E-01  0.485986E-01  0.487545E-01  0.489524E-01  0.491937E-01
  0.494799E-01  0.498123E-01  0.501927E-01  0.506226E-01  0.511036E-01
  0.516374E-01  0.522257E-01  0.530421E-01  0.537439E-01  0.545290E-01
  0.553822E-01  0.563362E-01  0.573829E-01  0.585322E-01  0.598104E-01
  0.611943E-01  0.627353E-01  0.644193E-01  0.662755E-01  0.683283E-01
*/
void LoadFLAPWCHARGEFiles() {
	StringList_pa FileNameStrList = TecUtilStringListAlloc();
	vector<string> FileNames;
	Boolean_t IsOk = TecUtilDialogGetFileNames(SelectFileOption_AllowMultiFileRead, &FileNameStrList, "FLAPW CHARGE File(s)", nullptr, "*");

	if (IsOk && FileNameStrList != nullptr) {
		for (int i = 1; i <= TecUtilStringListGetCount(FileNameStrList); ++i) {
			FileNames.emplace_back(TecUtilStringListGetString(FileNameStrList, i));
		}
	}

	if (!IsOk) {
		return;
	}

	// 	FileNames = {
	// 		"Z:\\Tecplot\\StorageWorkspace\\2020_FLAPW-dataLoader\\CHARGE"
	// 	};
	// 	for (int i = 0; i < 10; ++i) FileNames.push_back(FileNames.front());

		/*
		 * First validata all files and verify they're all for the same system
		 */

	int nx, ny, nz, ns, headerNumChars;
	unsigned long long int numBlockValChars, numBlockCharsTotal;
	float lx, ly, dx, dy, dz;
	string fieldName;
	vector<string> firstVal(FileNames.size());

	bool filesValid = GetFLAPWCHARGEFileHeader(FileNames.front(), fieldName, nx, ny, ns, lx, ly, firstVal[0], &headerNumChars, &numBlockValChars, &numBlockCharsTotal, &nz);

	std::map<string, string> PropNames;
	PropNames["total density"] = "Electron Density";
	PropNames["spin density"] = "Spin density";
	PropNames["spin up"] = "Spin majority density";
	PropNames["spin down"] = "Spin minority density";

	fieldName = PropNames[fieldName];
	
	if (!filesValid) {
		TecUtilDialogErrMsg(string("File header invalid (\'" + FileNames[0] + "\')").c_str());
		return;
	}

	/*
	 * prepare dataset
	 */

	int npts = nx * ny * nz;
	dx = lx / float(nx - 1);
	dy = ly / float(ny - 1);
	dz = dx;

	// Get Z spacing from user if different from the X spacing
	while (true) {
		string zinStr = to_string(dz);
		char * zoutCStr;
		TecUtilDialogGetSimpleText("Enter grid spacing between planes (default is the X spacing)", zinStr.c_str(), &zoutCStr);
		if (StringIsFloat(zoutCStr)){
			dz = atof(zoutCStr);
			break;
		}
		else{
			TecUtilDialogMessageBox("Please enter floating point value for interplane spacing", MessageBoxType_Warning);
		}
		TecUtilStringDealloc(&zoutCStr);
	}


	vector<int> IJK = { nx,ny,nz };


	/*
	*	Make dataset and volume zone if necessary
	*/
	StringList_pa VarNames = TecUtilStringListAlloc();
	vector<string> VarNameStrs = { "X", "Y", "Z", fieldName };

	for (string const & i : VarNameStrs)
		TecUtilStringListAppendString(VarNames, i.c_str());

#if defined(WIN32) || defined(_WIN32) || defined(__WIN32) && !defined(__CYGWIN__)
	string Slash = "\\";
#else
	string Slash = "/";
#endif

	string DataSetName = FileNames[0].substr(0, FileNames[0].find_last_of("."));
	size_t TmpStrPos = DataSetName.find_last_of(Slash.c_str());
	if (TmpStrPos != string::npos) {
		DataSetName = DataSetName.substr(TmpStrPos + 1, DataSetName.length() - TmpStrPos - 1);
	}
	vector<FieldDataType_e> VarDataTypes;
	EntIndex_t ZoneNum, VolZoneNum = -1;
	if (IsOk) {
		if (ReplaceDataSet) {
			VarDataTypes.resize(VarNameStrs.size(), FieldDataType_Float);
			IsOk = TecUtilDataSetCreate(StringMakeValidZoneName(DataSetName).c_str(), VarNames, TRUE);
		}
		else {
			if (!TecUtilDataSetIsAvailable()) {
				VarDataTypes.resize(VarNameStrs.size(), FieldDataType_Float);
				IsOk = TecUtilDataSetCreate(StringMakeValidZoneName(DataSetName).c_str(), VarNames, TRUE);
			}
			else {
				EntIndex_t NumVars = TecUtilDataSetGetNumVars();
				VarDataTypes.reserve(NumVars);
				for (int i = 1; i <= NumVars; ++i) {
					char* TmpCStr;
					TecUtilVarGetName(i, &TmpCStr);
					int Ind = SearchVectorForString(VarNameStrs, TmpCStr, true);
					TecUtilStringDealloc(&TmpCStr);
					if (Ind >= 0)
						VarDataTypes.push_back(FieldDataType_Float);
					else
						VarDataTypes.push_back(FieldDataType_Bit);
				}
			}
		}
		//if (VolZoneNum <= 0) {
		IsOk = TecUtilDataSetAddZone(StringMakeValidZoneName(CSMZoneName.FullVolume + DataSetName).c_str(), nx, ny, nz, ZoneType_Ordered, VarDataTypes.data());
		if (IsOk) {
			VolZoneNum = TecUtilDataSetGetNumZones();
			IsOk = VolZoneNum > 0;
		}
		//}
	}

	vector<FieldDataType_e> ZoneDataTypes(TecUtilDataSetGetNumZones(), FieldDataType_Bit);
	ZoneDataTypes[VolZoneNum - 1] = FieldDataType_Float;

	Set_pa VolZoneSet = TecUtilSetAlloc(TRUE);
	if (IsOk) {
		SetZoneStyle({ VolZoneNum }, ZoneStyle_Volume);
		TecUtilSetAddMember(VolZoneSet, VolZoneNum, TRUE);
	}

	/*
		*	Make new variables where necessary and get pointers to variable data
		*/
	vector<int> VarNums;
	int VarNum;
	// 	vector<FieldDataPointer_c> Ptrs(VarNameStrs.size());
	for (int i = 0; i < VarNameStrs.size() && IsOk; ++i) {
		if (!ReplaceDataSet) {
			VarNums.push_back(VarNumByName(VarNameStrs[i]));
			if (VarNums.back() <= 0) {
				IsOk = TecUtilDataSetAddVar(VarNameStrs[i].c_str(), ZoneDataTypes.data());
				if (IsOk) {
					VarNums.back() = TecUtilDataSetGetNumVars();
					IsOk = (VarNum > 0);
				}
			}
		}
		else
			VarNums.push_back(i + 1);

		// 		if (IsOk)
		// 			IsOk = Ptrs[i].InitializeWritePtr(VolZoneNum, VarNum);
	}


	// 	FieldDataPointer_c ValPtr;
	// 	ValPtr.InitializeWritePtr(VolZoneNum, VarNums.back());

		/*
			*	Now just need to parse through the char array and populate the necessary variables
			*/
	int TotalPoints = IJK[0] * IJK[1] * IJK[2];
	mat33 LatticeVector = zeros<mat>(3, 3);
	LatticeVector.at(0, 0) = dx;
	LatticeVector.at(1, 1) = dy;
	LatticeVector.at(2, 2) = dz;
	vec3 Origin = zeros<vec>(3);
	Origin = (LatticeVector * vec3({ nx - 1., ny - 1., nz - 1. })) * -0.5;
	string StatusStr = "Loading FLAPW data: " + DataSetName;
	StatusLaunch(StatusStr, AddOnID, TRUE);

// //	Generate XYZ Values (didn't speed things up or decrease memory to do it this way)
// 		{
// 			vector<float> X(nx*ny*nz), Y(nx*ny*nz), Z(nx*ny*nz);
// 			vector<float>* XYZ[] = { &X, &Y, &Z };
// 			vec3 TmpIJK = zeros<vec>(3);
// 			vec3 TmpXYZ;
// 	#pragma omp parallel for private(TmpIJK, TmpXYZ)
// 			for (int k = 0; k < IJK[2]; ++k) {
// 				TmpIJK[2] = k;
// 				TmpIJK[1] = 0;
// 				for (int j = 0; j < IJK[1]; ++j) {
// 					TmpIJK[0] = 0;
// 					for (int i = 0; i < IJK[0]; ++i) {
// 						int Index = IndexFromIJK(i + 1, j + 1, k + 1, IJK[0], IJK[1]) - 1;
// 						TmpXYZ = Origin + LatticeVector * TmpIJK;
// 						for (int dir = 0; dir < 3; ++dir) {
// 							XYZ[dir]->at(Index) = TmpXYZ[dir];
// 						}
// 	
// 						++TmpIJK[0];
// 					}
// 					++TmpIJK[1];
// 				}
// 			}
// 	
// 			// Save XYZ values
// 			for (int dir = 0; dir < 3; ++dir) {
// 				TecUtilDataLoadBegin(); // keep tecplot from deallocating pointers
// 				FieldData_pa VarRef = TecUtilDataValueGetWritableNativeRef(VolZoneNum, VarNums[dir]);
// 				if (VALID_REF(VarRef)) {
// 					TecUtilDataValueArraySetByRef(VarRef, 1, XYZ[dir]->size(), XYZ[dir]->data());
// 				}
// 				TecUtilDataLoadEnd();
// 			}
// 		}

	TecUtilDataLoadBegin(); // keep tecplot from deallocating pointers

	FieldVecPointer_c XYZPtr;
	XYZPtr.InitializeWritePtr(VolZoneNum, vector<int>(VarNums.begin(), VarNums.end() - 1));

	vector<float> SliceData(nx*ny);

	// write rho data as single write for each file
	FieldData_pa VarRef = TecUtilDataValueGetWritableNativeRef(VolZoneNum, VarNums.back());
	if (VALID_REF(VarRef)) {
		vec3 TmpIJK = zeros<vec>(3);
		ifstream InFile(FileNames[0]);
		high_resolution_clock::time_point Time1;
		Time1 = high_resolution_clock::now();
		vector<char> CharArray(numBlockValChars);
		for (int k = 0; k < IJK[2] && IsOk; ++k) {
			IsOk = StatusUpdate(k+1, IJK[2], StatusStr, AddOnID, Time1, false);
			// Get data from file
			unsigned long long int pos = k * numBlockCharsTotal + headerNumChars;
			if (GetFLAPWCHARGEFileFieldData(InFile, pos, CharArray, SliceData)) {
				TecUtilDataValueArraySetByRef(VarRef, IndexFromIJK(1, 1, k + 1, IJK[0], IJK[1]), SliceData.size(), SliceData.data());
			}
			else{
				int tmpint = 1;
			}

			TmpIJK[1] = 0;
			for (int j = 0; j < IJK[1] && IsOk; ++j) {
				TmpIJK[0] = 0;
				for (int i = 0; i < IJK[0]; ++i) {
					int Index = IndexFromIJK(i + 1, j + 1, k + 1, IJK[0], IJK[1]) - 1;
					XYZPtr.Write(Index, Origin + LatticeVector * TmpIJK);

					++TmpIJK[0];
				}
				++TmpIJK[1];
			}
			++TmpIJK[2];
		}
		InFile.close();
	}


	StatusDrop(AddOnID);

	XYZPtr.Close();


	if (!ReplaceDataSet) {
		TecUtilZoneSetScatter(SV_SHOW, VolZoneSet, 0.0, FALSE);
		TecUtilZoneSetMesh(SV_SHOW, VolZoneSet, 0.0, FALSE);
		TecUtilZoneSetContour(SV_SHOW, VolZoneSet, 0.0, FALSE);
		TecUtilZoneSetEdgeLayer(SV_SHOW, VolZoneSet, 0.0, TRUE);
	}

	//	Modify group number of volume zone
	if (IsOk) {
		ArgList_pa ArgList = TecUtilArgListAlloc();
		TecUtilArgListAppendString(ArgList, SV_P1, SV_FIELDMAP);
		TecUtilArgListAppendString(ArgList, SV_P2, SV_GROUP);
		TecUtilArgListAppendSet(ArgList, SV_OBJECTSET, VolZoneSet);
		TecUtilArgListAppendArbParam(ArgList, SV_IVALUE, (LgIndex_t)(GroupIndex));
		TecUtilStyleSetLowLevelX(ArgList);
		TecUtilArgListDealloc(&ArgList);
	}

	if (nz > 1) {
		// 3d, so turn on isosurface and stick a slice at 1/2 z
		TecUtilFrameSetPlotType(PlotType_Cartesian3D);

		double zMin, zMax;
		TecUtilDataValueGetMinMaxByZoneVar(VolZoneNum, VarNums[2], &zMin, &zMax);

		TecUtilMacroExecuteCommand("$!SLICELAYERS SHOW = YES");
		TecUtilMacroExecuteCommand("$!SLICEATTRIBUTES 1  SLICESURFACE = ZPLANES");
		TecUtilMacroExecuteCommand(string("$!SLICEATTRIBUTES 1  PRIMARYPOSITION{ Z = " + to_string((zMax - zMin) * 0.5) + " }").c_str());
		TecUtilMacroExecuteCommand("$!SLICEATTRIBUTES 1  CONTOUR{ CONTOURTYPE = FLOOD }");
		TecUtilMacroExecuteCommand("$!SLICEATTRIBUTES 1  CONTOUR{ FLOODCOLORING = GROUP1 }");

		ArgList_pa argList = TecUtilArgListAlloc();

		TecUtilArgListClear(argList);
		TecUtilArgListAppendString(argList, SV_P1, SV_ISOSURFACEATTRIBUTES);
		TecUtilArgListAppendString(argList, SV_P2, SV_EFFECTS);
		TecUtilArgListAppendString(argList, SV_P3, SV_USETRANSLUCENCY);
		TecUtilArgListAppendInt(argList, SV_OFFSET1, 1);
		TecUtilArgListAppendArbParam(argList, SV_IVALUE, TRUE);
		TecUtilStyleSetLowLevelX(argList);

		TecUtilArgListClear(argList);
		TecUtilArgListAppendString(argList, SV_P1, SV_ISOSURFACEATTRIBUTES);
		TecUtilArgListAppendString(argList, SV_P2, SV_EFFECTS);
		TecUtilArgListAppendString(argList, SV_P3, SV_SURFACETRANSLUCENCY);
		TecUtilArgListAppendInt(argList, SV_OFFSET1, 1);
		TecUtilArgListAppendArbParam(argList, SV_IVALUE, 70);
		TecUtilStyleSetLowLevelX(argList);

		TecUtilArgListClear(argList);
		TecUtilArgListAppendString(argList, SV_P1, SV_FIELDLAYERS);
		TecUtilArgListAppendString(argList, SV_P2, SV_SHOWSCATTER);
		TecUtilArgListAppendArbParam(argList, SV_IVALUE, TRUE);
		TecUtilStyleSetLowLevelX(argList);

		TecUtilArgListClear(argList);
		TecUtilArgListAppendString(argList, SV_P1, SV_ISOSURFACELAYERS);
		TecUtilArgListAppendString(argList, SV_P2, SV_SHOW);
		TecUtilArgListAppendArbParam(argList, SV_IVALUE, TRUE);
		TecUtilStyleSetLowLevelX(argList);

		TecUtilArgListClear(argList);
		TecUtilArgListAppendString(argList, SV_P1, SV_ISOSURFACEATTRIBUTES);
		TecUtilArgListAppendString(argList, SV_P2, SV_DEFINITIONCONTOURGROUP);
		TecUtilArgListAppendInt(argList, SV_OFFSET1, 1);
		TecUtilArgListAppendArbParam(argList, SV_IVALUE, (SmInteger_t)1);
		TecUtilStyleSetLowLevelX(argList);

		TecUtilArgListClear(argList);
		TecUtilArgListAppendString(argList, SV_P1, SV_ISOSURFACEATTRIBUTES);
		TecUtilArgListAppendString(argList, SV_P2, SV_SHOWGROUP);
		TecUtilArgListAppendInt(argList, SV_OFFSET1, 1);
		TecUtilArgListAppendArbParam(argList, SV_IVALUE, TRUE);
		TecUtilStyleSetLowLevelX(argList);

		EntIndex_t IsoSurfaceVarNum = VarNums.back();
		// 	double ValMax, ValMin;
		// 	TecUtilVarGetMinMax(IsoSurfaceVarNum, &ValMin, &ValMax);
		// 	double IsoValue = ValMax - 0.05 * (ValMax - ValMin);
		double IsoValue = 0.25;

		TecUtilArgListClear(argList);
		TecUtilArgListAppendString(argList, SV_P1, SV_ISOSURFACEATTRIBUTES);
		TecUtilArgListAppendString(argList, SV_P2, SV_ISOVALUE1);
		TecUtilArgListAppendInt(argList, SV_OFFSET1, 1);
		TecUtilArgListAppendDouble(argList, SV_DVALUE, IsoValue);
		TecUtilStyleSetLowLevelX(argList);


		TecUtilArgListClear(argList);
		TecUtilArgListAppendString(argList, SV_P1, SV_GLOBALCONTOUR);
		TecUtilArgListAppendString(argList, SV_P2, SV_VAR);
		TecUtilArgListAppendInt(argList, SV_OFFSET1, 1);
		TecUtilArgListAppendArbParam(argList, SV_IVALUE, IsoSurfaceVarNum);
		TecUtilStyleSetLowLevelX(argList);
		TecUtilArgListClear(argList);
		TecUtilArgListAppendString(argList, SV_P1, SV_ISOSURFACEATTRIBUTES);
		TecUtilArgListAppendString(argList, SV_P2, SV_ISOVALUE1);
		TecUtilArgListAppendInt(argList, SV_OFFSET1, 1);
		TecUtilArgListAppendDouble(argList, SV_DVALUE, IsoValue);
		TecUtilStyleSetLowLevelX(argList);
		TecUtilArgListClear(argList);
		TecUtilArgListAppendString(argList, SV_P1, SV_ISOSURFACEATTRIBUTES);
		TecUtilArgListAppendString(argList, SV_P2, SV_ISOVALUE2);
		TecUtilArgListAppendInt(argList, SV_OFFSET1, 1);
		TecUtilArgListAppendDouble(argList, SV_DVALUE, IsoValue);
		TecUtilStyleSetLowLevelX(argList);
		TecUtilArgListClear(argList);
		TecUtilArgListAppendString(argList, SV_P1, SV_ISOSURFACEATTRIBUTES);
		TecUtilArgListAppendString(argList, SV_P2, SV_ISOVALUE3);
		TecUtilArgListAppendInt(argList, SV_OFFSET1, 1);
		TecUtilArgListAppendDouble(argList, SV_DVALUE, IsoValue);
		TecUtilStyleSetLowLevelX(argList);

		TecUtilArgListDealloc(&argList);
	}
	else {
		// only one slice, so turn set to 2d and turn on contours for first zone
		TecUtilFrameSetPlotType(PlotType_Cartesian2D);
		TecUtilMacroExecuteCommand("$!FIELDLAYERS SHOWCONTOUR = YES");
		TecUtilZoneSetContour(SV_SHOW, VolZoneSet, 0.0, TRUE);
	}

	// set contour levels
	double varMin, varMax;
	TecUtilDataValueGetMinMaxByRef(VarRef, &varMin, &varMax);
	vec ContourLevels;
	if (varMin > 0) {
		ContourLevels = LogSpace(varMin, sqrt(varMax), 20);
	}
	else {
		ContourLevels = linspace(varMin, varMax, 20);
	}

	ArgList_pa TempArgList = TecUtilArgListAlloc();
	TecUtilArgListAppendInt(TempArgList, SV_CONTOURGROUP, 1);
	TecUtilArgListAppendInt(TempArgList, SV_VAR, VarNums.back());
	TecUtilContourSetVariableX(TempArgList);
	TecUtilArgListClear(TempArgList);

	TecUtilArgListAppendInt(TempArgList, SV_CONTOURLABELACTION, ContourLabelAction_DeleteAll);
	TecUtilArgListAppendInt(TempArgList, SV_CONTOURGROUP, 1);
	TecUtilContourLabelX(TempArgList);
	TecUtilArgListClear(TempArgList);

	TecUtilArgListAppendInt(TempArgList, SV_CONTOURLEVELACTION, ContourLevelAction_New);
	TecUtilArgListAppendInt(TempArgList, SV_CONTOURGROUP, 1);
	TecUtilArgListAppendInt(TempArgList, SV_NUMVALUES, 20);
	TecUtilArgListAppendArray(TempArgList, SV_RAWDATA, ContourLevels.memptr());
	TecUtilContourLevelX(TempArgList);
	TecUtilArgListClear(TempArgList);

	TecUtilArgListAppendString(TempArgList, SV_P1, SV_GLOBALCONTOUR);
	TecUtilArgListAppendString(TempArgList, SV_P2, SV_LEGEND);
	TecUtilArgListAppendString(TempArgList, SV_P3, SV_SHOW);
	TecUtilArgListAppendInt(TempArgList, SV_OFFSET1, 1);
	TecUtilArgListAppendArbParam(TempArgList, SV_IVALUE, FALSE);
	TecUtilStyleSetLowLevelX(TempArgList);
	TecUtilArgListClear(TempArgList);


	TecUtilDataLoadEnd();

	TecUtilSetDealloc(&VolZoneSet);

	AuxDataDataSetSetItem(DLProgramName, "FLAPW");

	return;
}


/*
*	Loads one or more FORMATTED turbomole cube files. Will fail if they're not formatted.
*/
void LoadTurboMoleCubeFiles()
{
	TecUtilLockStart(AddOnID);

	SYSTEM_INFO sysinfo;
	GetSystemInfo(&sysinfo);

	int numCPU = sysinfo.dwNumberOfProcessors;
	omp_set_num_threads(numCPU);

	StringList_pa FileNames = TecUtilStringListAlloc();
	Boolean_t IsOk = TecUtilDialogGetFileNames(SelectFileOption_AllowMultiFileRead, &FileNames, "TurboMole Cube File", nullptr, "*.cub?");

	Boolean_t UseExistingVolZone = FALSE;
	FileNameStrs.clear();

	if (IsOk && FileNames != nullptr){
		FileNameStrs.resize(TecUtilStringListGetCount(FileNames));
		for (int i = 1; i <= FileNameStrs.size(); ++i){
			if (i == 1){
				if (TecUtilDataSetIsAvailable()){
					ReplaceDataSet = !TecUtilDialogMessageBox("There is already data loaded. Do you want to append to or replace the current data?\n\tYes: Append to current data set\n\tNo: Replace current data set", MessageBoxType_YesNo);
					if (!ReplaceDataSet){
						UseExistingVolZone = TecUtilDialogMessageBox("Import data to existing full volume zone or create new zone?\n\tYes: Import to current zone\n\tNo: Import to new zone", MessageBoxType_YesNo);
						if (!UseExistingVolZone)
							ImportAtoms = TecUtilDialogMessageBox("Import atomic positions from all files(s)?", MessageBox_YesNo);
						else
							ImportAtoms = FALSE;
					}
				}
				else{
					ReplaceDataSet = TRUE;
					ImportAtoms = TRUE;
					UseExistingVolZone = TRUE;
				}
			}

			FileNameStrs[i - 1] = TecUtilStringListGetString(FileNames, i);
		}

		TecUtilStringListDealloc(&FileNames);

		/*
		 *	List of cube file title headers and corresponding variable names
		 */

		vector<vector<string> > VarNameList = {
			{
				"density",CSMVarName.Dens
			},
			{
				"gradient of density                                                       (norm)",CSMVarName.DensGradMag
			},
			{
				"2nd derivative of density                                              (trace/3)",CSMVarName.DensAvgCurv
			},
			{
				"kinetic energy density",CSMVarName.Tau
			},
			{
				"laplacian of density",CSMVarName.DensLap
			},
			{
				"electron localiz. fct.",CSMVarName.ELF
			},
			{
				"electrostatic potential",CSMVarName.ElecStatPot
			},
			{
				"electrostatic field                                                       (norm)",CSMVarName.ElecStatField
			},
			{
				"electrostatic fld. grd.                                                (trace/3)",CSMVarName.ElecStatFieldGrad
			},
			{
				"diamagetic shielding                                                   (trace/3)",CSMVarName.DiaMagnShift
			}
		};

		GroupIndex = 50;
		/*
		*	Loop over each file selected and import the data
		*/
		Boolean_t FirstFile = ReplaceDataSet;
		for (string const & CubeFileName : FileNameStrs){
			//	Open the cube file
			ifstream CubeFile(CubeFileName);
			if (!CubeFile.is_open()){
				TecUtilDialogErrMsg(string("Failed to open " + CubeFileName).c_str());
				return;
			}

			/*
			*	Parse the file header. Here's the structure:
			*
			*	Header line 1: variable or property in file
			*	Header line 2: "OUTER/INNER/MIDDLE LOOP: X,Y,Z"
			*	NAtoms, X-Origin, Y-Origin, Z-Origin NVal        NVal is the #points/value
			*	N1, X1, Y1, Z1               # of increments in the slowest running direction
			*	N2, X2, Y2, Z2
			*	N3, X3, Y3, Z3               # of increments in the fastest running direction
			*	IA1, Chg1, X1, Y1, Z1        Atomic number, charge, and coordinates of the first atom
			*	�
			*	IAn, Chgn, Xn, Yn, Zn        Atomic number, charge, and coordinates of the last atom
			*	(N1*N2) records, each of length N3     Values of the density at each point in the grid
			*/

			/*
			*	Get header lines
			*/

			string H1, H2;

			getline(CubeFile, H1);
			getline(CubeFile, H2);

			/*
			*	Test that there actually was a header.
			*	Just need to check that the first word of H1 is not an integer.
			*	(This assumes an integer will never be the first word in a formatted
			*	cube file, which may not be a bad assumption now and/or in the future)
			*/

			stringstream ss;
			ss << H1;
			string TmpStr;
			ss >> TmpStr;

			StatusLaunch("Reading file contents. Please wait.", AddOnID, FALSE);

			/*
			*	Now get the rest of the header
			*/

			string VarName;
			Boolean_t IsMOFile = FALSE;
			int NumAtoms;
			vec3 Origin;
			mat33 LatticeVector;
			vector<int> IJK(3);
			vector<AtomGroup_s> AtomInfo;

			/*
			*	Line three of the header
			*/
			CubeFile >> NumAtoms;
			NumAtoms = NumAtoms >= 0 ? NumAtoms : -NumAtoms;
			/*
			*	If the number of atoms < 0, then the cube file contains MO information.
			*/
			if (H1.substr(0,15) == "amplitude of MO"){
				IsMOFile = TRUE;
			}
			VarName = H1.substr(0,H1.find_last_not_of('  ')+1);
			for (auto const & i : VarNameList){
				if (VarName == i[0]){
					VarName = i[1];
					break;
				}
			}
			for (int i = 0; i < 3; ++i)
				CubeFile >> Origin[i];

			Origin = fliplr(Origin);

			/*
			*	There could be an extra value that I don't care about at the end of this line, so skip it.
			*/
			getline(CubeFile, TmpStr);

			/*
			*	Lines 4-6 of the header: number of points and lattice directions
			*/

			for (int i = 0; i < 3; ++i){
				CubeFile >> IJK[i];
				for (int j = 0; j < 3; ++j){
					CubeFile >> LatticeVector.at(i, j);
				}
			}

			/*
			*	Lines 7-(7+NumAtoms) of header: atom info
			*/
			vector<AtomColor_s> AtomColorList;
			PopulateAtomColorList(AtomColorList);

			for (int i = 0; i < NumAtoms; ++i){
				int ElemNum;
				double AtomChg;
				vec3 AtomPos;

				CubeFile >> ElemNum >> AtomChg;
				AtomChg = ElemNum;
				for (int j = 0; j < 3; ++j)
					CubeFile >> AtomPos[j];

				string AtomSymbol = GetAtomStr(ElemNum);

				/*
				*	See if we already have an AtomGroup_s for this atom type
				*/
				Boolean_t IsFound = FALSE;
				for (AtomGroup_s & a : AtomInfo){
					if (a.Name == AtomSymbol){
						IsFound = TRUE;

						a.AddPosition(AtomPos);
						a.Charges.push_back(AtomChg);
						a.Count++;

						break;
					}
				}

				if (!IsFound){
					AtomGroup_s NewAtomGroup(AtomSymbol);
					NewAtomGroup.AtomColor = GetAtomColor(AtomColorList, NewAtomGroup.Name);
					NewAtomGroup.AddPosition(AtomPos);
					NewAtomGroup.Charges.push_back(AtomChg);

					AtomInfo.push_back(NewAtomGroup);
				}
			}

			int NumVarBlocks = 1;

			/*
			*	Header is now read. Next fetch the first value in the file.
			*	After reading the entire file into memory, we'll use the first
			*	value to identify when we've arrived at the beginning of the
			*	data block.
			*/

			string FirstVal;
			CubeFile >> FirstVal;

			/*
			*	Read whole file into memory
			*/
			CubeFile.seekg(0, CubeFile.end);
			std::streamoff FileSize = CubeFile.tellg();
			CubeFile.seekg(0, CubeFile.beg);


			vector<char> CubeFileContents(FileSize);
			CubeFile.read(CubeFileContents.data(), FileSize);
			CubeFile.close();

			/*
			*	Run through the char vector CubeFileContents to find the first data value
			*/

			char* TmpCStr = strtok(CubeFileContents.data(), " \n\t");
			while (TmpCStr != nullptr && FirstVal.compare(TmpCStr))
				TmpCStr = strtok(nullptr, " \n\t");

			/*
			*	Make dataset and volume zone if necessary
			*/
			StringList_pa VarNames = TecUtilStringListAlloc();
			vector<string> VarNameStrs = { "X", "Y", "Z" , VarName};
			for (string const & i : VarNameStrs)
				TecUtilStringListAppendString(VarNames, i.c_str());

#if defined(WIN32) || defined(_WIN32) || defined(__WIN32) && !defined(__CYGWIN__)
			string Slash = "\\";
#else
			string Slash = "/";
#endif
			string DataSetName = "TurboMole Data";
			EntIndex_t ZoneNum, VolZoneNum = -1;
			if (IsOk){
				if (ReplaceDataSet){
					IsOk = TecUtilDataSetCreate(StringMakeValidZoneName(DataSetName).c_str(), VarNames, TRUE);
				}
				else{
					if (!TecUtilDataSetIsAvailable()){
						IsOk = TecUtilDataSetCreate(StringMakeValidZoneName(DataSetName).c_str(), VarNames, TRUE);
					}
// 					else{
// 						EntIndex_t NumVars = TecUtilDataSetGetNumVars();
// 						VarDataTypes.reserve(NumVars);
// 						for (int i = 1; i <= NumVars; ++i){
// 							char* TmpCStr;
// 							TecUtilVarGetName(i, &TmpCStr);
// 							int Ind = SearchVectorForString(VarNameStrs, TmpCStr, true);
// 							TecUtilStringDealloc(&TmpCStr);
// 							if (Ind >= 0)
// 								VarDataTypes.push_back(FieldDataType_Double);
// 							else
// 								VarDataTypes.push_back(FieldDataType_Bit);
// 						}
// 					}
				}
				if (IsOk){
					if (UseExistingVolZone){
						VolZoneNum = ZoneNumByName("Full Volume");
					}
					if (VolZoneNum <= 0){
						IsOk = TecUtilDataSetAddZone(StringMakeValidZoneName(CSMZoneName.FullVolume + DataSetName).c_str(), IJK[0], IJK[1], IJK[2], ZoneType_Ordered, vector<FieldDataType_e>(TecUtilDataSetGetNumVars(), FieldDataType_Float).data());
						if (IsOk){
							VolZoneNum = TecUtilDataSetGetNumZones();
							IsOk = VolZoneNum > 0;
						}
					}
				}
			}

			vector<FieldDataType_e> ZoneDataTypes(TecUtilDataSetGetNumZones(), FieldDataType_Bit);
			ZoneDataTypes[VolZoneNum - 1] = FieldDataType_Float;

			StatusDrop(AddOnID);

			TecUtilDataLoadBegin();

			Set_pa VolZoneSet = TecUtilSetAlloc(TRUE);
			if (IsOk) {
				SetZoneStyle({ VolZoneNum }, ZoneStyle_Volume);
				TecUtilSetAddMember(VolZoneSet, VolZoneNum, TRUE);
			}

			/*
			*	Make new variables where necessary and get pointers to variable data
			*/
			int VarNum;
			FieldDataPointer_c VarPtr;
			if (!ReplaceDataSet){
				VarNum = VarNumByName(VarName);
				if (VarNum <= 0){
					IsOk = TecUtilDataSetAddVar(VarName.c_str(), nullptr);
					if (IsOk){
						VarNum = TecUtilDataSetGetNumVars();
						IsOk = (VarNum > 0);
					}
				}
			}
			else
				VarNum = 4;

			if (IsOk)
				IsOk = VarPtr.InitializeWritePtr(VolZoneNum, VarNum);


			/*
			*	Now just need to parse through the char array and populate the necessary variables
			*/
			int TotalPoints = IJK[0] * IJK[1] * IJK[2] * NumVarBlocks;
			string StatusStr = "Loading TurboMole Cube file: " + DataSetName;
			StatusLaunch(StatusStr, AddOnID, TRUE);

			// now add the XYZ values (in the order ZYX)
			// x is the innermost loop for how the data is organized
			// 

			FieldVecPointer_c XYZPtr;
			if (FirstFile){
				XYZPtr.InitializeWritePtr(1, { 1, 2, 3 });
			}

			double Val;
			vec3 TmpIJK = zeros<vec>(3);
			for (int i = 0; i < IJK[0] && IsOk; ++i){
				IsOk = StatusUpdate(i, IJK[0], StatusStr, AddOnID);
				TmpIJK[1] = 0;
				for (int j = 0; j < IJK[1] && IsOk; ++j){
					TmpIJK[2] = 0;
					for (int k = 0; k < IJK[2]; ++k){
						int Index = IndexFromIJK(i + 1, j + 1, k + 1, IJK[0], IJK[1]) - 1;
						if (FirstFile){
							XYZPtr.Write(Index, Origin + LatticeVector * TmpIJK);
						}
// 						for (int Dir = 0; Dir < 3 && FirstFile; ++Dir){
// 							Val = Origin[2-Dir] + dot(TmpIJK, LatticeVector.col(Dir));
// 							XYZPtr[Dir].Write(Index, Val);
// 						}
						Val = atof(TmpCStr);
						VarPtr.Write(Index, Val);
						TmpCStr = strtok(nullptr, " \n\t");

						++TmpIJK[2];
					}
					++TmpIJK[1];
				}
				++TmpIJK[0];
			}

			FirstFile = FALSE;

			StatusDrop(AddOnID);

			if (TmpCStr == nullptr
				&& (TmpIJK[0] < IJK[0] || TmpIJK[1] < IJK[1] || TmpIJK[2] < IJK[2]))
			{
				TecUtilDialogErrMsg(string("Failed to load TurboMole file: " + CubeFileName).c_str());
				TecUtilDataLoadEnd();
				return;
			}


			VarPtr.Close();
			XYZPtr.Close();

			TecUtilDataLoadEnd();

			if (!ReplaceDataSet){
				TecUtilZoneSetScatter(SV_SHOW, VolZoneSet, 0.0, FALSE);
				TecUtilZoneSetMesh(SV_SHOW, VolZoneSet, 0.0, FALSE);
				TecUtilZoneSetContour(SV_SHOW, VolZoneSet, 0.0, FALSE);
				TecUtilZoneSetEdgeLayer(SV_SHOW, VolZoneSet, 0.0, TRUE);
			}

			//	Modify group number of volume zone
			if (IsOk){
				ArgList_pa ArgList = TecUtilArgListAlloc();
				TecUtilArgListAppendString(ArgList, SV_P1, SV_FIELDMAP);
				TecUtilArgListAppendString(ArgList, SV_P2, SV_GROUP);
				TecUtilArgListAppendSet(ArgList, SV_OBJECTSET, VolZoneSet);
				TecUtilArgListAppendArbParam(ArgList, SV_IVALUE, (LgIndex_t)(GroupIndex));
				TecUtilStyleSetLowLevelX(ArgList);
				TecUtilArgListDealloc(&ArgList);
			}

			TecUtilSetDealloc(&VolZoneSet);

			/*
			*	Now create zones for atoms
			*/

			if (IsOk && ReplaceDataSet || ImportAtoms){
				IsOk = CreateAtomZonesFromAtomGroupList(AtomInfo, { "X", "Y", "Z" }, vector<FieldDataType_e>(), GroupIndex);
			}

			/*
			*	Modify zone styles and such
			*/
			if (IsOk){
				//	Modify group number of volume zone
				ArgList_pa ArgList = TecUtilArgListAlloc();

				TecUtilStateChanged(StateChange_ZonesAdded, reinterpret_cast<ArbParam_t>(VolZoneSet));

				TecUtilFrameSetPlotType(PlotType_Cartesian3D);

				TecUtilZoneSetActive(VolZoneSet, AssignOp_PlusEquals);

				TecUtilSetDealloc(&VolZoneSet);


				TecUtilArgListClear(ArgList);
				TecUtilArgListAppendString(ArgList, SV_P1, SV_ISOSURFACEATTRIBUTES);
				TecUtilArgListAppendString(ArgList, SV_P2, SV_EFFECTS);
				TecUtilArgListAppendString(ArgList, SV_P3, SV_USETRANSLUCENCY);
				TecUtilArgListAppendInt(ArgList, SV_OFFSET1, 1);
				TecUtilArgListAppendArbParam(ArgList, SV_IVALUE, TRUE);
				TecUtilStyleSetLowLevelX(ArgList);

				TecUtilArgListClear(ArgList);
				TecUtilArgListAppendString(ArgList, SV_P1, SV_ISOSURFACEATTRIBUTES);
				TecUtilArgListAppendString(ArgList, SV_P2, SV_EFFECTS);
				TecUtilArgListAppendString(ArgList, SV_P3, SV_SURFACETRANSLUCENCY);
				TecUtilArgListAppendInt(ArgList, SV_OFFSET1, 1);
				TecUtilArgListAppendArbParam(ArgList, SV_IVALUE, 70);
				TecUtilStyleSetLowLevelX(ArgList);

				TecUtilArgListClear(ArgList);
				TecUtilArgListAppendString(ArgList, SV_P1, SV_FIELDLAYERS);
				TecUtilArgListAppendString(ArgList, SV_P2, SV_SHOWSCATTER);
				TecUtilArgListAppendArbParam(ArgList, SV_IVALUE, TRUE);
				TecUtilStyleSetLowLevelX(ArgList);

				TecUtilArgListClear(ArgList);
				TecUtilArgListAppendString(ArgList, SV_P1, SV_ISOSURFACELAYERS);
				TecUtilArgListAppendString(ArgList, SV_P2, SV_SHOW);
				TecUtilArgListAppendArbParam(ArgList, SV_IVALUE, TRUE);
				TecUtilStyleSetLowLevelX(ArgList);

				TecUtilArgListClear(ArgList);
				TecUtilArgListAppendString(ArgList, SV_P1, SV_ISOSURFACEATTRIBUTES);
				TecUtilArgListAppendString(ArgList, SV_P2, SV_DEFINITIONCONTOURGROUP);
				TecUtilArgListAppendInt(ArgList, SV_OFFSET1, 1);
				TecUtilArgListAppendArbParam(ArgList, SV_IVALUE, (SmInteger_t)1);
				TecUtilStyleSetLowLevelX(ArgList);

				TecUtilArgListClear(ArgList);
				TecUtilArgListAppendString(ArgList, SV_P1, SV_ISOSURFACEATTRIBUTES);
				TecUtilArgListAppendString(ArgList, SV_P2, SV_SHOWGROUP);
				TecUtilArgListAppendInt(ArgList, SV_OFFSET1, 1);
				TecUtilArgListAppendArbParam(ArgList, SV_IVALUE, TRUE);
				TecUtilStyleSetLowLevelX(ArgList);

				EntIndex_t IsoSurfaceVarNum = 4;
				// 	double ValMax, ValMin;
				// 	TecUtilVarGetMinMax(IsoSurfaceVarNum, &ValMin, &ValMax);
				// 	double IsoValue = ValMax - 0.05 * (ValMax - ValMin);
				double IsoValue = 0.25;

				TecUtilArgListClear(ArgList);
				TecUtilArgListAppendString(ArgList, SV_P1, SV_ISOSURFACEATTRIBUTES);
				TecUtilArgListAppendString(ArgList, SV_P2, SV_ISOVALUE1);
				TecUtilArgListAppendInt(ArgList, SV_OFFSET1, 1);
				TecUtilArgListAppendDouble(ArgList, SV_DVALUE, IsoValue);
				TecUtilStyleSetLowLevelX(ArgList);


				TecUtilArgListClear(ArgList);
				TecUtilArgListAppendString(ArgList, SV_P1, SV_GLOBALCONTOUR);
				TecUtilArgListAppendString(ArgList, SV_P2, SV_VAR);
				TecUtilArgListAppendInt(ArgList, SV_OFFSET1, 1);
				TecUtilArgListAppendArbParam(ArgList, SV_IVALUE, IsoSurfaceVarNum);
				TecUtilStyleSetLowLevelX(ArgList);
				TecUtilArgListClear(ArgList);
				TecUtilArgListAppendString(ArgList, SV_P1, SV_ISOSURFACEATTRIBUTES);
				TecUtilArgListAppendString(ArgList, SV_P2, SV_ISOVALUE1);
				TecUtilArgListAppendInt(ArgList, SV_OFFSET1, 1);
				TecUtilArgListAppendDouble(ArgList, SV_DVALUE, IsoValue);
				TecUtilStyleSetLowLevelX(ArgList);
				TecUtilArgListClear(ArgList);
				TecUtilArgListAppendString(ArgList, SV_P1, SV_ISOSURFACEATTRIBUTES);
				TecUtilArgListAppendString(ArgList, SV_P2, SV_ISOVALUE2);
				TecUtilArgListAppendInt(ArgList, SV_OFFSET1, 1);
				TecUtilArgListAppendDouble(ArgList, SV_DVALUE, IsoValue);
				TecUtilStyleSetLowLevelX(ArgList);
				TecUtilArgListClear(ArgList);
				TecUtilArgListAppendString(ArgList, SV_P1, SV_ISOSURFACEATTRIBUTES);
				TecUtilArgListAppendString(ArgList, SV_P2, SV_ISOVALUE3);
				TecUtilArgListAppendInt(ArgList, SV_OFFSET1, 1);
				TecUtilArgListAppendDouble(ArgList, SV_DVALUE, IsoValue);
				TecUtilStyleSetLowLevelX(ArgList);

				TecUtilArgListDealloc(&ArgList);

				TecUtilViewDataFit();
			}

			AuxDataDataSetSetItem(DLProgramName, "Gaussian 09");

			ReplaceDataSet = FALSE;
			GroupIndex++;
		}
	}

	TecUtilLockFinish(AddOnID);
}

Boolean_t LoadBinaryPLTFileData(){
	Boolean_t LoadSuccess = TRUE;

	StringList_pa FileNames;

	LoadSuccess = TecUtilDialogGetFileNames(SelectFileOption_ReadMultiFile, &FileNames, "TurboMole PLT", nullptr, "*.plt");

	if (LoadSuccess){
		FileNameStrs.resize(TecUtilStringListGetCount(FileNames));
		for (int i = 1; i <= FileNameStrs.size(); ++i){
			FileNameStrs[i - 1] = TecUtilStringListGetString(FileNames, i);
		}
		TecUtilStringListDealloc(&FileNames);
	}
	else return LoadSuccess;

	vector<int> XYZVarNums(3,-1);
	int LoadVarNum = -1;

	if (LoadSuccess && TecUtilDataSetIsAvailable()){
		LoadSuccess = TecUtilDialogMessageBox("Current data set will be replace with imported data.\nContinue?", MessageBox_YesNo);
	}
	else ReplaceDataSet = TRUE;

	vector<string> TypeNames = { "Electron Density", "Density Gradient", "Density Hessian" };
	vector<int> TypeNumbers = { 101, 102, 103 };

	for (int FileNum = 0; FileNum < FileNameStrs.size() && LoadSuccess; ++FileNum){
		ifstream File(FileNameStrs[FileNum].c_str(), ios::binary);
		LoadSuccess = File.is_open();

		int IntData[5];
		/*
		* The five integers represent, in order:
		* 	Rank: must be 3 (for 3D)
		* 	Type: will be 101, but don't know why...
		* 	number of points in x dir
		* 	number of points in y dir
		* 	number of points in z dir
		*/
		float MinMaxZYX[6];

		if (LoadSuccess){
			File.read(reinterpret_cast<char*>(IntData), 5 * sizeof(int));
			File.read(reinterpret_cast<char*>(MinMaxZYX), 6 * sizeof(float));

			//		for (int i = 0; i < 5; ++i)
			//			cout << IntData[i] << '\n';
			//
			//		for (int i = 0; i < 6; ++i)
			//					cout << MinMaxZYX[i] << '\n';

			if (IntData[0] != 3){
				std::cerr << "Rank of PLT file must be 3. Quitting.";
				LoadSuccess = FALSE;
			}
		}

		int TypeNum;
		for (TypeNum = 0; TypeNum < TypeNames.size(); ++TypeNum){
			if (TypeNumbers[TypeNum] == IntData[1]) break;
		}
		if (LoadSuccess){
			if (FileNum == 0){
				vector<string> VarNamesVec = { "X", "Y", "Z", TypeNames[TypeNum] };
				StringList_pa VarNames = TecUtilStringListAlloc();
				LoadVarNum = 4;
				for (auto const & i : VarNamesVec)
					TecUtilStringListAppendString(VarNames, i.c_str());
				LoadSuccess = TecUtilDataSetCreate("TurboMole Data", VarNames, TRUE);
				TecUtilStringListDealloc(&VarNames);

				if (LoadSuccess){
					vector<FieldDataType_e> DataTypes(4, FieldDataType_Float);
					LoadSuccess = TecUtilDataSetAddZone("CSMZoneName.FullVolume", IntData[4], IntData[3], IntData[2], ZoneType_Ordered, DataTypes.data());
				}

				// now add the XYZ values (in the order ZYX)
				// x is the innermost loop for how the data is organized
				vector<FieldDataPointer_c> ZYXPtrs(3);
				int Index;
				float ZYXVal[3], DelZYX[3];
				TecUtilDataLoadBegin();
				for (int i = 0; i < 3 && LoadSuccess; ++i){
					LoadSuccess = ZYXPtrs[i].InitializeWritePtr(1, 3 - i);
					DelZYX[i] = (MinMaxZYX[i * 2 + 1] - MinMaxZYX[i * 2]) / (static_cast<float>(IntData[2 + i] - 1));
				}
				ZYXVal[0] = MinMaxZYX[0]; // set z value
				for (int k = 0; k < IntData[2]; ++k){
					ZYXVal[1] = MinMaxZYX[2]; // reset y value
					for (int j = 0; j < IntData[3]; ++j){
						ZYXVal[2] = MinMaxZYX[4]; // reset x value
						for (int i = 0; i < IntData[4]; ++i){
							Index = IndexFromIJK(i + 1, j + 1, k + 1, IntData[4], IntData[3]) - 1;
							for (int Dir = 0; Dir < 3; ++Dir){
								ZYXPtrs[Dir].Write(Index, ZYXVal[Dir]);
							}
							ZYXVal[2] += DelZYX[2];
						}
						ZYXVal[1] += DelZYX[1];
					}
					ZYXVal[0] += DelZYX[0];
				}
				for (auto & i : ZYXPtrs) i.Close();

				TecUtilDataLoadEnd();
			}
			else{
				vector<FieldDataType_e> DataTypes(1, FieldDataType_Float);
				LoadSuccess = TecUtilDataSetAddVar(TypeNames[TypeNum].c_str(), DataTypes.data());
				LoadVarNum = TecUtilDataSetGetNumVars();
			}
		}

		if (LoadSuccess){
			vector<float> Values(IntData[2] * IntData[3] * IntData[4]);
			File.read(reinterpret_cast<char*>(Values.data()), Values.size() * sizeof(float));
			TecUtilDataLoadBegin();
			FieldData_pa LoadVarPtr = TecUtilDataValueGetWritableNativeRef(1, LoadVarNum);
			TecUtilDataValueArraySetByRef(LoadVarPtr, 1, Values.size(), Values.data());
			TecUtilDataLoadEnd();
		}

		if (File.is_open())
			File.close();

	}

	if (LoadSuccess){

		int VolZoneNum = TecUtilDataSetGetNumZones();


		SetZoneStyle({ VolZoneNum }, ZoneStyle_Volume);

		TecUtilFrameSetPlotType(PlotType_Cartesian3D);

		Set_pa VolZoneSet = TecUtilSetAlloc(TRUE);

		TecUtilSetAddMember(VolZoneSet, VolZoneNum, TRUE);

		TecUtilZoneSetScatter(SV_SHOW, VolZoneSet, 0.0, FALSE);
		TecUtilZoneSetMesh(SV_SHOW, VolZoneSet, 0.0, FALSE);
		TecUtilZoneSetContour(SV_SHOW, VolZoneSet, 0.0, FALSE);
		TecUtilZoneSetEdgeLayer(SV_SHOW, VolZoneSet, 0.0, TRUE);
	}


	EntIndex_t IsoSurfaceVarNum = VarNumByName(TypeNames[0]);

	if (IsoSurfaceVarNum > 0){
		ArgList_pa argList = TecUtilArgListAlloc();

		TecUtilArgListClear(argList);
		TecUtilArgListAppendString(argList, SV_P1, SV_ISOSURFACEATTRIBUTES);
		TecUtilArgListAppendString(argList, SV_P2, SV_EFFECTS);
		TecUtilArgListAppendString(argList, SV_P3, SV_USETRANSLUCENCY);
		TecUtilArgListAppendInt(argList, SV_OFFSET1, 1);
		TecUtilArgListAppendArbParam(argList, SV_IVALUE, TRUE);
		TecUtilStyleSetLowLevelX(argList);

		TecUtilArgListClear(argList);
		TecUtilArgListAppendString(argList, SV_P1, SV_ISOSURFACEATTRIBUTES);
		TecUtilArgListAppendString(argList, SV_P2, SV_EFFECTS);
		TecUtilArgListAppendString(argList, SV_P3, SV_SURFACETRANSLUCENCY);
		TecUtilArgListAppendInt(argList, SV_OFFSET1, 1);
		TecUtilArgListAppendArbParam(argList, SV_IVALUE, 70);
		TecUtilStyleSetLowLevelX(argList);

		TecUtilArgListClear(argList);
		TecUtilArgListAppendString(argList, SV_P1, SV_FIELDLAYERS);
		TecUtilArgListAppendString(argList, SV_P2, SV_SHOWSCATTER);
		TecUtilArgListAppendArbParam(argList, SV_IVALUE, TRUE);
		TecUtilStyleSetLowLevelX(argList);

		TecUtilArgListClear(argList);
		TecUtilArgListAppendString(argList, SV_P1, SV_ISOSURFACELAYERS);
		TecUtilArgListAppendString(argList, SV_P2, SV_SHOW);
		TecUtilArgListAppendArbParam(argList, SV_IVALUE, TRUE);
		TecUtilStyleSetLowLevelX(argList);

		TecUtilArgListClear(argList);
		TecUtilArgListAppendString(argList, SV_P1, SV_ISOSURFACEATTRIBUTES);
		TecUtilArgListAppendString(argList, SV_P2, SV_DEFINITIONCONTOURGROUP);
		TecUtilArgListAppendInt(argList, SV_OFFSET1, 1);
		TecUtilArgListAppendArbParam(argList, SV_IVALUE, (SmInteger_t)1);
		TecUtilStyleSetLowLevelX(argList);

		TecUtilArgListClear(argList);
		TecUtilArgListAppendString(argList, SV_P1, SV_ISOSURFACEATTRIBUTES);
		TecUtilArgListAppendString(argList, SV_P2, SV_SHOWGROUP);
		TecUtilArgListAppendInt(argList, SV_OFFSET1, 1);
		TecUtilArgListAppendArbParam(argList, SV_IVALUE, TRUE);
		TecUtilStyleSetLowLevelX(argList);

		// 	double ValMax, ValMin;
		// 	TecUtilVarGetMinMax(IsoSurfaceVarNum, &ValMin, &ValMax);
		// 	double IsoValue = ValMax - 0.05 * (ValMax - ValMin);
		double IsoValue = 0.25;

		TecUtilArgListClear(argList);
		TecUtilArgListAppendString(argList, SV_P1, SV_ISOSURFACEATTRIBUTES);
		TecUtilArgListAppendString(argList, SV_P2, SV_ISOVALUE1);
		TecUtilArgListAppendInt(argList, SV_OFFSET1, 1);
		TecUtilArgListAppendDouble(argList, SV_DVALUE, IsoValue);
		TecUtilStyleSetLowLevelX(argList);


		TecUtilArgListClear(argList);
		TecUtilArgListAppendString(argList, SV_P1, SV_GLOBALCONTOUR);
		TecUtilArgListAppendString(argList, SV_P2, SV_VAR);
		TecUtilArgListAppendInt(argList, SV_OFFSET1, 1);
		TecUtilArgListAppendArbParam(argList, SV_IVALUE, IsoSurfaceVarNum);
		TecUtilStyleSetLowLevelX(argList);
		TecUtilArgListClear(argList);
		TecUtilArgListAppendString(argList, SV_P1, SV_ISOSURFACEATTRIBUTES);
		TecUtilArgListAppendString(argList, SV_P2, SV_ISOVALUE1);
		TecUtilArgListAppendInt(argList, SV_OFFSET1, 1);
		TecUtilArgListAppendDouble(argList, SV_DVALUE, IsoValue);
		TecUtilStyleSetLowLevelX(argList);

		TecUtilArgListDealloc(&argList);
	}

	//	for (auto const & i : Values)
	//		cout << i << '\n';

	return LoadSuccess;
}