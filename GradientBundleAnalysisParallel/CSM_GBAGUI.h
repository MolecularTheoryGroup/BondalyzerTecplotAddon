#ifndef CSMGBAGUI_H_
#define CSMGBAGUI_H_

#include <string>
#include <vector>

using std::vector;
using std::string;



/*
*	Default values for GUI elements
*/
int const DefaultCP = 1;
string const DefaultIntVarStr = "Electron Density";
Boolean_t const DefaultIntegrate = TRUE;
Boolean_t const DefaultVolIntegrate = TRUE;
Boolean_t const DefaultSystemIsOpen = TRUE;
// double const DefaultRhoCutoff = 0.001;
double const DefaultRadius = 0.2;
int const DefaultLevel = 3;
int const MinLevel = 1;
int const MaxLevel = 6;
int const DefaultSTPts = 100;
double const DefaultIBDist = 0.05;
double const DefaultIBAng = 20;


vector<string> ListGetSelectedStrings(LgIndex_t  ListID);
vector<int> ListGetSelectedItemNums(LgIndex_t ListID);
void ListPopulateWithVarNames(LgIndex_t ListID);
vec3 GetCoordsFromListItem(LgIndex_t ItemIndex,
	LgIndex_t ListIndex,
	string* ItemFullString = nullptr,
	string* ItemNameString = nullptr,
	int* ItemNumber = nullptr,
	int* ZoneNumber = nullptr,
	bool* IsCP = nullptr,
	vector<int>* NumberOfCPs = nullptr);


#endif