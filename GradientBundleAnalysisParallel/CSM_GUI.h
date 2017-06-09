#ifndef CSMGUI_H_
#define CSMGUI_H_

#include <string>
#include <vector>

using std::vector;
using std::string;



/*
*	Default values for GUI elements
*/
const int DefaultCP = 1;
const string DefaultIntVarStr = "Electron Density";
const Boolean_t DefaultIntegrate = TRUE;
const Boolean_t DefaultVolIntegrate = TRUE;
const Boolean_t DefaultSystemIsOpen = TRUE;
// const double DefaultRhoCutoff = 0.001;
const double DefaultRadius = 0.2;
const int DefaultLevel = 3;
const int MinLevel = 0;
const int MaxLevel = 6;
const int DefaultSTPts = 100;
const double DefaultIBDist = 0.05;
const double DefaultIBAng = 20;


const vector<string> ListGetSelectedStrings(const LgIndex_t &  ListID);
const vector<int> ListGetSelectedItemNums(const LgIndex_t & ListID);
void ListPopulateWithVarNames(LgIndex_t ListID);
const vec3 GetCoordsFromListItem(const LgIndex_t & ItemIndex,
	const LgIndex_t & ListIndex,
	string* ItemFullString = NULL,
	string* ItemNameString = NULL,
	int* ItemNumber = NULL,
	int* ZoneNumber = NULL,
	bool* IsCP = NULL,
	vector<int>* NumberOfCPs = NULL);


#endif