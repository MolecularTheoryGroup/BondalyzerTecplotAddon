#if defined MSWIN
#include "STDAFX.h"
#endif

#include <vector>
#include <string>
#include <algorithm>

#include "TECADDON.h"

#include "CSM_CRIT_POINTS.h"
#include "CSM_GUI.h"
#include "CSM_DATA_SET_INFO.h"

using std::string;
using std::vector;
using std::to_string;

const static int LineHeight = 150;
const static int CharWidth = 110;
const static int VertSpacing = 50;
const static int HorzSpacing = 100;

const static int ListNumLines = 10;

vector<GuiField_c> CSMGuiFields;

AddOn_pa CSMGuiAddOnID;

int CSMGuiMultiListID = BADDIALOGID;
bool CSMGuiMultiListIsPointSelect;
int CSMGuiPointSelectZoneNum = -1;
const static string CSMGuiPointSelectDelim = ": ";
vector<Text_ID> CSMGuiPointLabelIDs;

int CSMGuiDialogManager;
bool CSMGuiSuccess;
CSMGuiReturnFunc_pf CurrentReturnFunc;

const vector<GuiFieldType_e> 
	IntTypes = {
	Gui_ZoneSelect,
	Gui_ZoneSelectInt,
	Gui_VarSelect,
	Gui_VarSelectInt,
	Gui_ZonePointSelectMulti
},
	BoolTypes = {
		Gui_Toggle,
		Gui_ToggleEnable
	},
	StringTypes = {},
	DoubleTypes = {
		Gui_Double
	},
	IntVecTypes = {
		Gui_VarSelectMulti,
		Gui_ZoneSelectMulti,
		Gui_ZonePointSelectMulti
	};

GuiField_c::GuiField_c()
{
}

GuiField_c::GuiField_c(const GuiFieldType_e & Type,
	const string & Label,
	const string & Val)
{
	Type_m = Type;
	Label_m = Label;
	InputVal_m = Val;
}

GuiField_c::~GuiField_c()
{
}

void GuiField_c::CheckIntType() const{
#ifdef _DEBUG
	bool IsCorrect = false;
	for (const auto & t : IntTypes) IsCorrect = (IsCorrect || Type_m == t);
	if (!IsCorrect) TecUtilDialogErrMsg("Type int requested from non-int CSM gui field");
#endif
}

void GuiField_c::CheckDoubleType() const{
#ifdef _DEBUG
	bool IsCorrect = false;
	for (const auto & t : DoubleTypes) IsCorrect = (IsCorrect || Type_m == t);
	if (!IsCorrect) TecUtilDialogErrMsg("Type double requested from non-double CSM gui field");
#endif
}

void GuiField_c::CheckBoolType() const{
#ifdef _DEBUG
	bool IsCorrect = false;
	for (const auto & t : BoolTypes) IsCorrect = (IsCorrect || Type_m == t);
	if (!IsCorrect) TecUtilDialogErrMsg("Type bool requested from non-bool CSM gui field");
#endif
}

void GuiField_c::CheckStringType() const{
#ifdef _DEBUG
	bool IsCorrect = false;
	for (const auto & t : StringTypes) IsCorrect = (IsCorrect || Type_m == t);
	if (!IsCorrect) TecUtilDialogErrMsg("Type string requested from non-string CSM gui field");
#endif
}

void GuiField_c::CheckIntVecType() const{
#ifdef _DEBUG
	bool IsCorrect = false;
	for (const auto & t : IntVecTypes) IsCorrect = (IsCorrect || Type_m == t);
	if (!IsCorrect) TecUtilDialogErrMsg("Type integer vector requested from non-IntVec CSM gui field");
#endif
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

const vector<int> GuiField_c::GetReturnIntVec() const{
	CheckIntVecType();
	return ValueIntVec_m;
}


void GuiField_c::SetReturnInt(const int & Val){
	CheckIntType();
	ValueInt_m = Val;
}

void GuiField_c::SetReturnDouble(const double & Val){
	CheckDoubleType();
	ValueDouble_m = Val;
}

void GuiField_c::SetReturnString(const string & Val){
	CheckStringType();
	ValueString_m = Val;
}

void GuiField_c::SetReturnBool(const bool & Val){
	CheckBoolType();
	ValueBool_m = Val;
}

void GuiField_c::SetReturnIntVec(const vector<int> & Val){
	CheckIntVecType();
	ValueIntVec_m = Val;
}

const vec3 GetPointCoordsFromListItem(const int & ItemIndex,
	const int & PointZoneNum){
	vec3 Out;

	int XYZVarNums[3] = { -1, -1, -1 };
	TecUtilAxisGetVarAssignments(&XYZVarNums[0], &XYZVarNums[1], &XYZVarNums[2]);
	if (XYZVarNums[0] < 0 || XYZVarNums[1] < 0 || XYZVarNums[2] < 0){
		for (int i = 1; i <= 3; ++i){
			XYZVarNums[i - 1] = i;
		}
	}

	for (int i = 0; i < 3; ++i) Out[i] = TecUtilDataValueGetByZoneVar(PointZoneNum, XYZVarNums[i], ItemIndex);

	return Out;
}

void CSMGuiLabelSelectedPoints(){
	LgIndex_t * SelectedItems;
	LgIndex_t NumSelected;
	TecGUIListGetSelectedItems(CSMGuiMultiListID, &SelectedItems, &NumSelected);
	if (NumSelected > 0){
		TecUtilDrawGraphics(FALSE);
		TecUtilPickDeselectAll();
		for (int i = 0; i < CSMGuiPointLabelIDs.size(); ++i){
			TecUtilPickText(CSMGuiPointLabelIDs[i]);
		}
		TecUtilPickClear();
		CSMGuiPointLabelIDs.clear();

		double TextOffset = 0.1;
		for (int i = 0; i < NumSelected; ++i){
			char* LabelStr = TecGUIListGetString(CSMGuiMultiListID, SelectedItems[i]);

			vec3 Point = GetPointCoordsFromListItem(SelectedItems[i], CSMGuiPointSelectZoneNum);
			vec3 XYZGrid;
			TecUtilConvert3DPositionToGrid(Point[0], Point[1], Point[2], &XYZGrid[0], &XYZGrid[1], &XYZGrid[2]);

			Text_ID LabelID = TecUtilTextCreate(CoordSys_Grid, XYZGrid[0] + TextOffset, XYZGrid[1] + TextOffset, Units_Frame, 2, LabelStr);
			TecUtilTextBoxSetType(LabelID, TextBox_Filled);
			TecUtilTextBoxSetFillColor(LabelID, White_C);

			Boolean_t IsValid = TecUtilTextIsValid(LabelID);

			CSMGuiPointLabelIDs.push_back(LabelID);

			TecUtilStringDealloc(&LabelStr);
		}
		TecUtilDrawGraphics(TRUE);
	}
	TecUtilArrayDealloc((void **)&SelectedItems);
}

void CSMGuiPointSelectMultiListCallback(const LgIndex_t* Val){
	TecUtilLockStart(CSMGuiAddOnID);
	CSMGuiLabelSelectedPoints();
	TecUtilLockFinish(CSMGuiAddOnID);
}

void CSMGuiMultiListSelectAllButtonCallBack(){
	TecUtilLockStart(CSMGuiAddOnID);
	LgIndex_t* SelectedNums;
	LgIndex_t NumSelected, ItemCount = TecGUIListGetItemCount(CSMGuiMultiListID);
	TecGUIListGetSelectedItems(CSMGuiMultiListID, &SelectedNums, &NumSelected);
	TecUtilArrayDealloc((void**)&SelectedNums);

	if (ItemCount > NumSelected){
		TecGUIListSelectAllItems(CSMGuiMultiListID);
	}
	else if (ItemCount > 0){
		TecGUIListSetSelectedItem(CSMGuiMultiListID, 1);
	}

	CSMGuiPointSelectMultiListCallback(NULL);
	TecUtilLockFinish(CSMGuiAddOnID);
}

void CSMGUIDeleteCPLabels(){
	TecUtilLockStart(CSMGuiAddOnID);
	MouseButtonMode_e MouseMode = TecUtilMouseGetCurrentMode();
	if (CSMGuiPointLabelIDs.size() > 0){
		TecUtilDrawGraphics(FALSE);
		TecUtilPickDeselectAll();
		for (int i = 0; i < CSMGuiPointLabelIDs.size(); ++i){
			TecUtilPickText(CSMGuiPointLabelIDs[i]);
		}
		TecUtilPickClear();
		CSMGuiPointLabelIDs.clear();
		TecUtilDrawGraphics(TRUE);
		TecUtilRedraw(TRUE);
	}
	TecUtilMouseSetMode(MouseMode);
	TecUtilLockFinish(CSMGuiAddOnID);
}

void DialogCloseButtonCB(){
	TecUtilLockStart(CSMGuiAddOnID);
	CSMGuiSuccess = false;
	CSMGUIDeleteCPLabels();

	CSMGuiMultiListID = BADDIALOGID;
	CSMGuiPointSelectZoneNum = -1;

	TecGUIDialogDrop(CSMGuiDialogManager);
	TecUtilLockFinish(CSMGuiAddOnID);
}

void DialogOKButtonCB(){

	CSMGuiSuccess = true;

	CSMGUIDeleteCPLabels();

	CSMGuiMultiListID = BADDIALOGID;
	CSMGuiPointSelectZoneNum = -1;

	TecGUIDialogDrop(CSMGuiDialogManager);

	for (auto & f : CSMGuiFields){
		GuiFieldType_e t = f.GetType();

		if (t <= Gui_VarSelect) f.SetReturnInt(TecGUIOptionMenuGet(f.GetID()));
		else if (t <= Gui_VarSelectInt) f.SetReturnInt(atoi(TecGUITextFieldGetString(f.GetID())));
		else if (t <= Gui_Double) f.SetReturnDouble(atof(TecGUITextFieldGetString(f.GetID())));
		else if (t <= Gui_ToggleEnable) f.SetReturnBool(TecGUIToggleGet(f.GetID()));
		else if (t <= Gui_ZonePointSelectMulti){
			LgIndex_t* SelectedNums;
			LgIndex_t NumSelected;
			vector<int> IntVec;

			TecGUIListGetSelectedItems(f.GetID(), &SelectedNums, &NumSelected);
			for (int i = 0; i < NumSelected; ++i) IntVec.push_back(SelectedNums[i]);
			f.SetReturnIntVec(IntVec);
			TecUtilArrayDealloc((void**)&SelectedNums);

			if (t == Gui_ZonePointSelectMulti) f.SetReturnInt(TecGUIOptionMenuGet(f.GetID(1)));
		}
	}

	CurrentReturnFunc(CSMGuiSuccess, CSMGuiFields);
}

void CSMGuiPointSelectOptionCallback(const int* Val){
	if (CSMGuiMultiListID == BADDIALOGID || CSMGuiPointSelectZoneNum == *Val || !AuxDataZoneItemMatches(*Val, CSMAuxData.CC.ZoneType, CSMAuxData.CC.ZoneTypeCPs)) return;

	TecUtilLockStart(CSMGuiAddOnID);

	CSMGuiPointSelectZoneNum = *Val;
	int IJK[3];
	TecUtilZoneGetIJK(CSMGuiPointSelectZoneNum, &IJK[0], &IJK[1], &IJK[2]);
	int NumPoints = IJK[0];

	bool ZoneIsCP = AuxDataZoneItemMatches(CSMGuiPointSelectZoneNum, CSMAuxData.CC.ZoneType, CSMAuxData.CC.ZoneTypeCPs);
	string LabelBase = "";
	FieldData_pa CPTypeRef;
	if (ZoneIsCP){
		int CPTypeVarNum = VarNumByName(CSMVarName.CritPointType);
		if (CPTypeVarNum > 0){
			CPTypeRef = TecUtilDataValueGetReadableNativeRef(CSMGuiPointSelectZoneNum, CPTypeVarNum);
		}
		else ZoneIsCP = false;
		if (!VALID_REF(CPTypeRef)) ZoneIsCP = false;
	}

	TecGUIListDeleteAllItems(CSMGuiMultiListID);
	for (int i = 0; i < NumPoints; ++i){
		if (ZoneIsCP){
			int CPType = TecUtilDataValueGetByRef(CPTypeRef, i + 1);
			int CPInd = std::find(CPTypeList, CPTypeList + 6, CPType) - CPTypeList;
			LabelBase = CPNameList[CPInd] + " ";
		}
		TecGUIListAppendItem(CSMGuiMultiListID, string(LabelBase + to_string(i + 1)).c_str());
	}

	TecUtilLockFinish(CSMGuiAddOnID);
}

int TextFieldZoneNumCheck(const char* Text){

	if (!StringIsInt(Text)){
		TecUtilDialogErrMsg("Please enter an integer zone number");
		return 0;
	}
	else {
		int ZoneNum = atoi(Text);
		if (ZoneNum <= 0 || ZoneNum > TecUtilDataSetGetNumZones()){
			TecUtilDialogErrMsg("Invalid zone number");
			return 0;
		}
	}

	return 1;
}

int TextFieldVarNumCheck(const char* Text){

	if (!StringIsInt(Text)){
		TecUtilDialogErrMsg("Please enter an integer variable number");
		return 0;
	}
	else {
		int ZoneNum = atoi(Text);
		if (ZoneNum <= 0 || ZoneNum > TecUtilDataSetGetNumVars()){
			TecUtilDialogErrMsg("Invalid variable number");
			return 0;
		}
	}

	return 1;
}

int TextCallbackDoNothing(const char* Text){
	return 0;
}

void IntCallbackDoNothing(const int* Int){
	return;
}

static void VoidCallbackDoNothing(const int* iVal){
	return;
}


void CSMLaunchGui(const string & Title, const vector<GuiField_c> & Fields)
{
	TecUtilPleaseWait("Loading dialog...", TRUE);

	CSMGuiFields = Fields;

	vector<GuiFieldType_e> FieldType(int(Gui_Invalid)+1);
	for (int i = 0; i < FieldType.size(); ++i) FieldType[i] = GuiFieldType_e(i);
	vector<int> NumFields(FieldType.size(), 0);
	vector<vector<GuiFieldType_e> > ZoneVarFields = {
		{
			Gui_ZoneSelect,
			Gui_ZoneSelectInt,
			Gui_ZoneSelectMulti,
			Gui_ZonePointSelectMulti
		},
		{
			Gui_VarSelect,
			Gui_VarSelectInt,
			Gui_VarSelectMulti
		}
	};

	int MaxLabelWidth = -1,
		MaxOnlyLabelWidth = -1,
		MaxZoneVarWidth = -1;

	vector<vector<int> > FieldZoneVarNums(CSMGuiFields.size());

	for (const auto & f : CSMGuiFields) {
		NumFields[int(f.GetType())]++;
		if (f.GetType() < Gui_Toggle) MaxLabelWidth = MAX(MaxLabelWidth, int(f.GetLabel().length()));
		else MaxOnlyLabelWidth = MAX(MaxOnlyLabelWidth, int(f.GetLabel().length()));
	}

	vector<vector<string> > ZoneVarList(2);
	vector<bool> ZoneVarCommaListValid(2, true);
	vector<string> ZoneVarCommaList(2); // comma-delimited list of zones/vars for option menus
	vector<int> NumZonesVars = { TecUtilDataSetGetNumZones(), TecUtilDataSetGetNumVars() };

	for (int i = 0; i < 2; ++i) for (const auto & t : ZoneVarFields[i]) if (NumFields[int(t)]){ // Get full list of zone/var names
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
		break;
	}

	for (int f = 0; f < CSMGuiFields.size(); ++f){
		GuiFieldType_e t = CSMGuiFields[f].GetType();
		if (t == Gui_ZoneSelect || t == Gui_ZoneSelectInt || t == Gui_ZoneSelectMulti || t == Gui_ZonePointSelectMulti){
			for (const auto & s : SplitString(CSMGuiFields[f].GetSearchString(), ',')){
				FieldZoneVarNums[f].push_back(ZoneNumByName(s, false, true));
			}
		}
		else if (t == Gui_VarSelect || t == Gui_VarSelectInt || t == Gui_VarSelectMulti){
			for (const auto & s : SplitString(CSMGuiFields[f].GetSearchString(), ',')){
				FieldZoneVarNums[f].push_back(VarNumByName(s, true));
			}
		}
	}

	int  MultiListNumLines = ListNumLines;
	if (CSMGuiFields.size() < 10) MultiListNumLines *= 2;

	MaxLabelWidth = MAX(MaxLabelWidth, 15);
	int X = HorzSpacing, Y = LineHeight;
	int W = CharWidth * MAX((MaxLabelWidth + MaxZoneVarWidth), MaxOnlyLabelWidth) + HorzSpacing;
	int H = 2 * LineHeight + VertSpacing;

	int xTmp = X - HorzSpacing + CharWidth * MaxLabelWidth;
	int wTmp = CharWidth * MaxZoneVarWidth;

	for (int t = 0; t < NumFields.size(); ++t){
		if (GuiFieldType_e(t) <= Gui_Label) H += (LineHeight + VertSpacing) * NumFields[t];
		else if (GuiFieldType_e(t) <= Gui_ZonePointSelectMulti) H += (LineHeight * (MultiListNumLines + 1.5)) * NumFields[t];
		else H += (2 * VertSpacing) * NumFields[t];
	}

	CSMGuiDialogManager = TecGUIDialogCreateModeless(MAINDIALOGID, W, H, Title.c_str(), NULL, DialogCloseButtonCB, NULL);

	int fNum = 0;

	for (auto & f : CSMGuiFields){
		GuiFieldType_e t = f.GetType();
		if (t < Gui_VertSep){
			TecGUILabelAdd(CSMGuiDialogManager, X, Y + VertSpacing * 0.5, f.GetLabel().c_str());

			if (t <= Gui_VarSelect){ // Gui Field is an option menu for selecting zone/variable
				string TmpStr;
				if (ZoneVarCommaListValid[int(t)]) TmpStr = ZoneVarCommaList[int(t)];
				f.SetID(TecGUIOptionMenuAdd(CSMGuiDialogManager, xTmp, Y, wTmp, LineHeight, ZoneVarCommaList[int(t)].c_str(), VoidCallbackDoNothing));
				if (!ZoneVarCommaListValid[int(t)]) for (const auto & s : ZoneVarList[int(t)]) TecGUIOptionMenuAppendItem(f.GetID(), s.c_str());
				TecGUIOptionMenuSet(f.GetID(), MAX(1, FieldZoneVarNums[fNum].back()));
			}
			else if (t <= Gui_Double){ // Gui field is a text input box (maybe for selecting zone/variable)
				TecGUITextCallback_pf fCB;
				switch (t){
					case Gui_ZoneSelectInt: fCB = TextFieldZoneNumCheck; break;
					case Gui_VarSelectInt: fCB = TextFieldVarNumCheck; break;
					default: fCB = TextCallbackDoNothing; break;
				}
				f.SetID(TecGUITextFieldAdd(CSMGuiDialogManager, xTmp, Y, wTmp, LineHeight, fCB));
				if (t < Gui_Double) TecGUITextFieldSetLgIndex(f.GetID(), MAX(1, FieldZoneVarNums[fNum].back()), FALSE);
				else TecGUITextFieldSetDouble(f.GetID(), stof(f.GetSearchString()), "%f");
			}
			else if (t <= Gui_ToggleEnable){
				int w = f.GetLabel().length() * CharWidth;
				f.SetID(TecGUIToggleAdd(CSMGuiDialogManager, X, Y + VertSpacing * 0.5, w, LineHeight, f.GetLabel().c_str(), VoidCallbackDoNothing));
			}
			else if (t <= Gui_ZonePointSelectMulti){
				vector<string> ListItems;
				TecGUIIntCallback_pf fCB;
				if (t == Gui_ZonePointSelectMulti){
					string TmpStr;
					if (ZoneVarCommaListValid[0]) TmpStr = ZoneVarCommaList[0];
					f.SetID(TecGUIOptionMenuAdd(CSMGuiDialogManager, xTmp, Y, wTmp, LineHeight, ZoneVarCommaList[0].c_str(), CSMGuiPointSelectOptionCallback), 1);
					if (!ZoneVarCommaListValid[0]) for (const auto & s : ZoneVarList[0]) TecGUIOptionMenuAppendItem(f.GetID(1), s.c_str());
					TecGUIOptionMenuSet(f.GetID(1), MAX(1, FieldZoneVarNums[fNum].back()));
					fCB = CSMGuiPointSelectMultiListCallback;
					CSMGuiMultiListIsPointSelect = true;
					Y += LineHeight + VertSpacing;
				}
				else{
					ListItems = ZoneVarList[int(t - Gui_ZoneSelectMulti)]; // This returns the zone list for ZoneSelectMulti and the var list for VarSelectMulti
					CSMGuiMultiListIsPointSelect = false;
					fCB = IntCallbackDoNothing;
				}
				f.SetID(TecGUIListAdd(CSMGuiDialogManager, xTmp, Y, wTmp, MultiListNumLines * LineHeight, TRUE, fCB));
				Y += LineHeight;
				TecGUIButtonAdd(CSMGuiDialogManager, X, Y, CharWidth * 10, LineHeight, "Toggle all", CSMGuiMultiListSelectAllButtonCallBack);
				CSMGuiMultiListID = f.GetID();
				if (t == Gui_ZonePointSelectMulti){
					int Val = MAX(1, FieldZoneVarNums[fNum].back());
					CSMGuiPointSelectOptionCallback(&Val);
				}
				else for (const auto & s : ListItems) TecGUIListAppendItem(f.GetID(), s.c_str());
				if (t != Gui_ZonePointSelectMulti && FieldZoneVarNums[fNum].size() > 0){
					vector<LgIndex_t> SelectNums;
					int MinNum = INT_MAX;
					for (const auto & i : FieldZoneVarNums[fNum]) if (i > 0){
						SelectNums.push_back(i);
						MinNum = MIN(MinNum, i);
					}
					if (SelectNums.size() > 0){
						TecGUIListSetSelectedItems(f.GetID(), SelectNums.data(), SelectNums.size());
						TecGUIListSetTopItemNum(f.GetID(), MinNum);
					}
				}
				Y += LineHeight * (MultiListNumLines - 2);
			}

			Y += LineHeight + VertSpacing;
		}
		else if (t == Gui_VertSep){
			f.SetID(TecGUIVertSeparatorAdd(CSMGuiDialogManager, X, Y, VertSpacing));
			Y += VertSpacing * 2;
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
					GuiFieldType_e t = CSMGuiFields[fi].GetType();
					IsEnabled = ((t == Gui_VarSelect || 
						t == Gui_VarSelectInt || 
						t == Gui_VarSelectMulti ||
						t == Gui_ZoneSelect || 
						t == Gui_ZoneSelectInt ||
						t == Gui_ZoneSelectMulti) && 
						FieldZoneVarNums[fi].size() > 0);
					if (IsEnabled){
						int NumHits = 0;
						for (const auto & i : FieldZoneVarNums[fi]) NumHits += int(i > 0);
						IsEnabled = (NumHits > 0);
					}
				}
				TecGUIToggleSet(f.GetID(), IsEnabled);
			}
		}
	}

	// Set toggle default values based on toggle Val "0" or blank for false and "1" for true
	if (NumFields[int(Gui_Toggle)]){
		for (const auto & f : CSMGuiFields){
			if (f.GetType() == Gui_Toggle){
				TecGUIToggleSet(f.GetID(), f.GetSearchString().length() > 0 && f.GetSearchString() != "0");
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

void CSMGui(const string & Title, const vector<GuiField_c> & Fields, CSMGuiReturnFunc_pf ReturnFunc, const AddOn_pa & InputAddOnID){
	CurrentReturnFunc = ReturnFunc;
	CSMGuiAddOnID = InputAddOnID;
	CSMLaunchGui(Title, Fields);
}