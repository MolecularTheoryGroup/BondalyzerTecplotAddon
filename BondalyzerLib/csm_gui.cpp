#if defined MSWIN
#include "STDAFX.h"
#endif

#include <vector>
#include <string>

#include "TECADDON.h"

#include "CSM_GUI.h"
#include "CSM_DATA_SET_INFO.h"

using std::string;
using std::vector;
using std::to_string;

const static int LineHeight = 150;
const static int CharWidth = 110;
const static int VertSpacing = 50;
const static int HorzSpacing = 100;

vector<GuiField_c> CSMGuiFields;

int CSMGuiDialogManager;
bool CSMGuiSuccess;
CSMGuiReturnFunc_pf CurrentReturnFunc;

const vector<GuiFieldType_e> 
	IntTypes = {
	Gui_ZoneSelect,
	Gui_ZoneSelectInt,
	Gui_VarSelect,
	Gui_VarSelectInt
},
	BoolTypes = {
		Gui_Toggle,
		Gui_ToggleEnable
	},
	StringTypes = {},
	DoubleTypes = {};

GuiField_c::GuiField_c()
{
}

GuiField_c::GuiField_c(const GuiFieldType_e & Type,
	const string & Label,
	const string & SearchString)
{
	Type_m = Type;
	Label_m = Label;
	SearchString_m = SearchString;
}

GuiField_c::~GuiField_c()
{
}

void GuiField_c::CheckIntType() const{
	bool IsCorrect = false;
	for (const auto & t : IntTypes) 
		IsCorrect = (IsCorrect || Type_m == t);
	if (!IsCorrect) TecUtilDialogErrMsg("Type int requested from non-int CSM gui field");
}

void GuiField_c::CheckDoubleType() const{
	bool IsCorrect = false;
	for (const auto & t : DoubleTypes) IsCorrect = (IsCorrect || Type_m == t);
	if (!IsCorrect) TecUtilDialogErrMsg("Type double requested from non-double CSM gui field");
}

void GuiField_c::CheckBoolType() const{
	bool IsCorrect = false;
	for (const auto & t : BoolTypes) IsCorrect = (IsCorrect || Type_m == t);
	if (!IsCorrect) TecUtilDialogErrMsg("Type bool requested from non-bool CSM gui field");
}

void GuiField_c::CheckStringType() const{
	bool IsCorrect = false;
	for (const auto & t : StringTypes) IsCorrect = (IsCorrect || Type_m == t);
	if (!IsCorrect) TecUtilDialogErrMsg("Type string requested from non-string CSM gui field");
}


const int GuiField_c::GetReturnInt() const{
	CheckIntType();
	return ValueInt_m;
}

const double GuiField_c::GetReturnDouble() const{
	CheckDoubleType();
	return ValueDouble_m;
}

const string GuiField_c::GetReturnString() const{
	CheckStringType();
	return ValueString_m;
}

const bool GuiField_c::GetReturnBool() const{
	CheckBoolType();
	return ValueBool_m;
}


void GuiField_c::SetReturnInt(const int Val){
	CheckIntType();
	ValueInt_m = Val;
}

void GuiField_c::SetReturnDouble(const double Val){
	CheckDoubleType();
	ValueDouble_m = Val;
}

void GuiField_c::SetReturnString(const string Val){
	CheckStringType();
	ValueString_m = Val;
}

void GuiField_c::SetReturnBool(const bool Val){
	CheckBoolType();
	ValueBool_m = Val;
}


void DialogCloseButtonCB(){
	// 	TecUtilLockStart(AddOnID);
	CSMGuiSuccess = false;
	TecGUIDialogDrop(CSMGuiDialogManager);
	// 	TecUtilLockFinish(AddOnID);
}

void DialogOKButtonCB(){

	CSMGuiSuccess = true;

	TecGUIDialogDrop(CSMGuiDialogManager);

	for (auto & f : CSMGuiFields){
		GuiFieldType_e t = f.GetType();

		if (t == Gui_ZoneSelect || t == Gui_VarSelect) f.SetReturnInt(TecGUIOptionMenuGet(f.GetID()));
		else if (t == Gui_ZoneSelectInt || t == Gui_VarSelectInt) f.SetReturnInt(atoi(TecGUITextFieldGetString(f.GetID())));
		else if (t == Gui_Toggle || t == Gui_ToggleEnable) f.SetReturnBool(TecGUIToggleGet(f.GetID()));
	}

	CurrentReturnFunc(CSMGuiSuccess, CSMGuiFields);
}

int TextFieldZoneNumCheck(const char* Text){

	if (!StringIsInt(Text))	TecUtilDialogErrMsg("Please enter an integer zone number");
	else {
		int ZoneNum = atoi(Text);
		if (ZoneNum <= 0 || ZoneNum > TecUtilDataSetGetNumZones()) TecUtilDialogErrMsg("Invalid zone number");
	}

	return 0;
}

int TextFieldVarNumCheck(const char* Text){

	if (!StringIsInt(Text))	TecUtilDialogErrMsg("Please enter an integer variable number");
	else {
		int ZoneNum = atoi(Text);
		if (ZoneNum <= 0 || ZoneNum > TecUtilDataSetGetNumVars()) TecUtilDialogErrMsg("Invalid variable number");
	}

	return 0;
}

static void IntCallbackDoNothing(const int* iVal){
	return;
}


void CSMLaunchGui(const string & Title, const vector<GuiField_c> & Fields)
{
	TecUtilPleaseWait("Loading dialog...", TRUE);

	CSMGuiFields = Fields;

	vector<GuiFieldType_e> FieldType(7);
	for (int i = 0; i < FieldType.size(); ++i) FieldType[i] = GuiFieldType_e(i);
	vector<int> NumFields(FieldType.size(), 0);
	vector<GuiFieldType_e> CheckFields = {
		Gui_ZoneSelect,
		Gui_VarSelect,
		Gui_ToggleEnable
	};

	int MaxLabelWidth = -1,
		MaxZoneVarWidth = -1;

	vector<int> FieldZoneVarNum(CSMGuiFields.size(), -1);

	for (const auto & f : CSMGuiFields) {
		NumFields[int(f.GetType())]++;
		if (f.GetType() != Gui_Toggle && f.GetType() != Gui_ToggleEnable) 
			MaxLabelWidth = MAX(MaxLabelWidth, int(f.GetLabel().length()));
	}

	vector<vector<string> > ZoneVarList(2);
	vector<bool> ZoneVarCommaListValid(2, true);
	vector<string> ZoneVarCommaList(2); // comma-delimited list of zones/vars for option menus
	vector<int> NumZonesVars = { TecUtilDataSetGetNumZones(), TecUtilDataSetGetNumVars() };

	for (int i = 0; i < 2; ++i) if (NumFields[int(CheckFields[i])]){ // Get full list of zone/var names
		char *cStr;
		for (int j = 1; j <= NumZonesVars[i]; ++j){
			if ((i == 0 && TecUtilZoneGetName(j, &cStr)) || TecUtilVarGetName(j, &cStr)){
				ZoneVarList[i].push_back(cStr);
				MaxZoneVarWidth = MAX(MaxZoneVarWidth, int(ZoneVarList[i].back().length()));
				ZoneVarCommaListValid[i] = (ZoneVarCommaListValid[i] && SplitString(ZoneVarList[i].back(), ',').size() == 1);
			}
			TecUtilStringDealloc(&cStr);
		}
		ZoneVarCommaListValid[i] = (ZoneVarCommaListValid[i] && ZoneVarList[i].size() == NumZonesVars[i]);
		if (ZoneVarCommaListValid[i]){
			for (int j = 0; j < NumZonesVars[i] - 1; ++j) ZoneVarCommaList[i] += ZoneVarList[i][j] + ',';
			ZoneVarCommaList[i] += ZoneVarList[i].back();
		}
	}

	for (int f = 0; f < CSMGuiFields.size(); ++f){
		GuiFieldType_e t = CSMGuiFields[f].GetType();
// 		if (t == Gui_ZoneSelect || t == Gui_ZoneSelectInt) FieldZoneVarNum[f] = SearchVectorForString(ZoneVarList[0], CSMGuiFields[f].GetSearchString());
// 		else if (t == Gui_VarSelect || t == Gui_VarSelectInt) FieldZoneVarNum[f] = SearchVectorForString(ZoneVarList[1], CSMGuiFields[f].GetSearchString());
		if (t == Gui_ZoneSelect || t == Gui_ZoneSelectInt) FieldZoneVarNum[f] = ZoneNumByName(CSMGuiFields[f].GetSearchString(), false, true);
		else if (t == Gui_VarSelect || t == Gui_VarSelectInt) FieldZoneVarNum[f] = VarNumByName(CSMGuiFields[f].GetSearchString(), true);
	}

	int X = HorzSpacing, Y = LineHeight;
	int W = CharWidth * (MaxLabelWidth + MaxZoneVarWidth) + HorzSpacing;
	int H = 0;

	int xTmp = X - HorzSpacing + CharWidth * MaxLabelWidth;
	int wTmp = CharWidth * MaxZoneVarWidth;

	for (int t = 0; t < NumFields.size(); ++t){
		if (GuiFieldType_e(t) != Gui_VertSep) H += (LineHeight + VertSpacing) * NumFields[t];
		else H += (2 * VertSpacing) * NumFields[t];
	}
// 	H += LineHeight + VertSpacing;

	CSMGuiDialogManager = TecGUIDialogCreateModeless(MAINDIALOGID, W, H, Title.c_str(), NULL, DialogCloseButtonCB, NULL);

	int fNum = 0;

	for (auto & f : CSMGuiFields){
		if (f.GetType() != Gui_VertSep){
			if (f.GetType() != Gui_Toggle && f.GetType() != Gui_ToggleEnable && f.GetType() != Gui_VertSep){
				TecGUILabelAdd(CSMGuiDialogManager, X, Y + VertSpacing, f.GetLabel().c_str());
				switch (f.GetType()){
					case Gui_ZoneSelectInt:
						f.SetID(TecGUITextFieldAdd(CSMGuiDialogManager, xTmp, Y, wTmp, LineHeight, TextFieldZoneNumCheck));
						TecGUITextFieldSetLgIndex(f.GetID(), MAX(1, FieldZoneVarNum[fNum]), FALSE);
						break;
					case Gui_VarSelectInt:
						f.SetID(TecGUITextFieldAdd(CSMGuiDialogManager, xTmp, Y, wTmp, LineHeight, TextFieldVarNumCheck));
						TecGUITextFieldSetLgIndex(f.GetID(), MAX(1, FieldZoneVarNum[fNum]), FALSE);
						break;
					case Gui_ZoneSelect:
						if (ZoneVarCommaListValid[0])
							f.SetID(TecGUIOptionMenuAdd(CSMGuiDialogManager, xTmp, Y, wTmp, LineHeight, ZoneVarCommaList[0].c_str(), IntCallbackDoNothing));
						else{
							f.SetID(TecGUIOptionMenuAdd(CSMGuiDialogManager, xTmp, Y, wTmp, LineHeight, "", IntCallbackDoNothing));
							for (const auto & s : ZoneVarList[0]) TecGUIOptionMenuAppendItem(f.GetID(), s.c_str());
						}
						TecGUIOptionMenuSet(f.GetID(), MAX(1, FieldZoneVarNum[fNum]));
						break;
					case Gui_VarSelect:
						if (ZoneVarCommaListValid[1])
							f.SetID(TecGUIOptionMenuAdd(CSMGuiDialogManager, xTmp, Y, wTmp, LineHeight, ZoneVarCommaList[1].c_str(), IntCallbackDoNothing));
						else{
							f.SetID(TecGUIOptionMenuAdd(CSMGuiDialogManager, xTmp, Y, wTmp, LineHeight, "", IntCallbackDoNothing));
							for (const auto & s : ZoneVarList[1]) TecGUIOptionMenuAppendItem(f.GetID(), s.c_str());
						}
						TecGUIOptionMenuSet(f.GetID(), MAX(1, FieldZoneVarNum[fNum]));
						break;
				}
				Y += LineHeight + VertSpacing;
			}
			else if (f.GetType() == Gui_Toggle || f.GetType() == Gui_ToggleEnable){
				int w = f.GetLabel().length() * CharWidth;
// 				int x = (W / 2) - (w / 2) + HorzSpacing;
				f.SetID(TecGUIToggleAdd(CSMGuiDialogManager, X, Y + VertSpacing, w, LineHeight, f.GetLabel().c_str(), IntCallbackDoNothing));
				Y += LineHeight + VertSpacing;
			}
			else if (f.GetType() == Gui_VertSep){
				f.SetID(TecGUIVertSeparatorAdd(CSMGuiDialogManager, X, Y, VertSpacing));
				Y += VertSpacing * 2;
			}
		}
		fNum++;
	}

	int OkButton = TecGUIButtonAdd(CSMGuiDialogManager, X + HorzSpacing, Y, W - (4 * HorzSpacing), LineHeight, "OK", DialogOKButtonCB);
	TecGUIButtonSetDefault(CSMGuiDialogManager, OkButton);

	// Now check that Gui_ToggleEnable fields are satisfied
	if (NumFields[int(Gui_ToggleEnable)]){
		for (const auto & f : CSMGuiFields){
			if (f.GetType() == Gui_ToggleEnable){
				vector<string> DependentFields = SplitString(f.GetSearchString(), ',');
				bool IsEnabled = true;
				for (int s = 0; s < DependentFields.size() && IsEnabled; ++s) if (StringIsInt(DependentFields[s])){
					int fi = stoi(DependentFields[s]);
					IsEnabled = ((CSMGuiFields[fi].GetType() == Gui_VarSelect || 
						CSMGuiFields[fi].GetType() == Gui_VarSelectInt || 
						CSMGuiFields[fi].GetType() == Gui_ZoneSelect || 
						CSMGuiFields[fi].GetType() == Gui_ZoneSelectInt) && 
						FieldZoneVarNum[fi] > 0);
				}
				TecGUIToggleSet(f.GetID(), IsEnabled);
			}
		}
	}

	TecUtilPleaseWait("Loading dialog...", FALSE);

	TecGUIDialogLaunch(CSMGuiDialogManager);

	// #ifdef MSWIN
	// 	RECT desktop;
	// 	const HWND hDesktop = GetDesktopWindow();
	// 	GetWindowRect(hDesktop, &desktop);
	// 
	// 	desktop.right; // screen width
	// 	desktop.bottom; // screen height
	// #else
	// 
	// #endif
}

void CSMGui(const string & Title, const vector<GuiField_c> & Fields, CSMGuiReturnFunc_pf ReturnFunc){
	CurrentReturnFunc = ReturnFunc;
	CSMLaunchGui(Title, Fields);
}