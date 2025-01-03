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

static int const LineHeight = 150;
static int const CharWidth = 110;
static int const VertSpacing = 50;
static int const HorzSpacing = 100;

static int const ListNumLines = 7;

vector<GuiField_c> CSMGuiFields;
vector<GuiField_c> CSMPassthroughFields;

AddOn_pa CSMGuiAddOnID;
bool CSMGuiDialogUp = false;

int CSMGuiMultiListMaxNumMultiListFields = 5;
vector<int> CSMGuiMultiListIDs(5, BADDIALOGID); //remains hardcoded to 5 bacause you have to code extra copies of everything below (see copies of void callback function pointers below)
vector<bool> CSMGuiMultiListIsPointSelects(5);

int CSMGuiMultiListID = BADDIALOGID;
bool CSMGuiMultiListIsPointSelect;
int CSMGuiPointSelectZoneNum = -1;
static string const CSMGuiPointSelectDelim = ": ";
vector<Text_ID> CSMGuiPointLabelIDs;
Set_pa CSMGuiActiveZones = nullptr;
bool CSMGuiProbeFirstInstall = true;
bool CSMGuiDeleteLabelsRunning = false;
MouseButtonMode_e CSMGuiMouseMode = MouseButtonMode_Invalid;
Boolean_t CSMGuiAutoRedrawIsActive = TRUE;

int CSMGuiDialogManager;
bool CSMGuiSuccess;
CSMGuiReturnFunc_pf CurrentReturnFunc;

vector<GuiFieldType_e> const 
	IntTypes = {
		Gui_ZoneSelect,
		Gui_ZoneSelectInt,
		Gui_VarSelect,
		Gui_VarSelectInt,
		Gui_ZonePointSelectMulti,
		Gui_Int,
		Gui_Radio
	},
	BoolTypes = {
		Gui_Toggle,
		Gui_ToggleEnable
	},
	StringTypes = {
		Gui_String
	},
	DoubleTypes = {
		Gui_Double
	},
	IntVecTypes = {
		Gui_VarSelectMulti,
		Gui_ZoneSelectMulti,
		Gui_ZonePointSelectMulti
	},
	StringVecTypes = {
		Gui_VarSelectMulti,
		Gui_ZoneSelectMulti,
		Gui_ZonePointSelectMulti
	};


/*
 *	Begin GuiField_c methods
 */

GuiField_c::GuiField_c()
{
}

GuiField_c::GuiField_c(GuiFieldType_e const & Type,
	string const & Label,
	string const & Val,
	const vector<TecGUIVoidCallback_pf> & CallbackFuntions)
{
	Type_m = Type;
	Label_m = Label;
	InputVal_m = Val;
	CallbackFuntions_m = CallbackFuntions;
}

GuiField_c::~GuiField_c()
{
}

void GuiField_c::CheckIntType() const{
#ifdef _DEBUG
	bool IsCorrect = false;
	for (auto const & t : IntTypes) IsCorrect = (IsCorrect || Type_m == t);
	if (!IsCorrect) TecUtilDialogErrMsg("Type int requested from non-int CSM gui field");
#endif
}

void GuiField_c::CheckDoubleType() const{
#ifdef _DEBUG
	bool IsCorrect = false;
	for (auto const & t : DoubleTypes) IsCorrect = (IsCorrect || Type_m == t);
	if (!IsCorrect) TecUtilDialogErrMsg("Type double requested from non-double CSM gui field");
#endif
}

void GuiField_c::CheckBoolType() const{
#ifdef _DEBUG
	bool IsCorrect = false;
	for (auto const & t : BoolTypes) IsCorrect = (IsCorrect || Type_m == t);
	if (!IsCorrect) TecUtilDialogErrMsg("Type bool requested from non-bool CSM gui field");
#endif
}

void GuiField_c::CheckStringType() const{
#ifdef _DEBUG
	bool IsCorrect = false;
	for (auto const & t : StringTypes) IsCorrect = (IsCorrect || Type_m == t);
	if (!IsCorrect) TecUtilDialogErrMsg("Type string requested from non-string CSM gui field");
#endif
}

void GuiField_c::CheckIntVecType() const{
#ifdef _DEBUG
	bool IsCorrect = false;
	for (auto const & t : IntVecTypes) IsCorrect = (IsCorrect || Type_m == t);
	if (!IsCorrect) TecUtilDialogErrMsg("Type integer vector requested from non-IntVec CSM gui field");
#endif
}

void GuiField_c::CheckStringVecType() const{
#ifdef _DEBUG
	bool IsCorrect = false;
	for (auto const & t : IntVecTypes) IsCorrect = (IsCorrect || Type_m == t);
	if (!IsCorrect) TecUtilDialogErrMsg("Type integer vector requested from non-IntVec CSM gui field");
#endif
}


int GuiField_c::GetReturnInt() const{
	CheckIntType();
	return ValueInt_m;
}

double GuiField_c::GetReturnDouble() const{
	CheckDoubleType();
	return ValueDouble_m;
}

string GuiField_c::GetReturnString() const{
	CheckStringType();
	return ValueString_m;
}

bool GuiField_c::GetReturnBool() const{
	CheckBoolType();
	return ValueBool_m;
}

vector<int> GuiField_c::GetReturnIntVec() const{
	CheckIntVecType();
	return ValueIntVec_m;
}

vector<string> GuiField_c::GetReturnStringVec() const{
	CheckStringVecType();
	return ValueStringVec_m;
}


void GuiField_c::SetReturnInt(int Val){
	CheckIntType();
	ValueInt_m = Val;
}

void GuiField_c::SetReturnDouble(double const & Val){
	CheckDoubleType();
	ValueDouble_m = Val;
}

void GuiField_c::SetReturnString(string const & Val){
	CheckStringType();
	ValueString_m = Val;
}

void GuiField_c::SetReturnBool(bool Val){
	CheckBoolType();
	ValueBool_m = Val;
}

void GuiField_c::SetReturnIntVec(vector<int> const & Val){
	CheckIntVecType();
	ValueIntVec_m = Val;
}

void GuiField_c::SetReturnStringVec(vector<string> const & Val){
	CheckStringVecType();
	ValueStringVec_m = Val;
}

/*
 *	End GuiField_c methods
 */

vec3 const GetPointCoordsFromListItem(int ItemIndex,
	int PointZoneNum){
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

void CSMGUIDeleteCPLabels(AddOn_pa *AddOnID);
void CSMGuiLabelSelectedPoints(AddOn_pa *AddOnID){
	if (AddOnID == nullptr) AddOnID = &CSMGuiAddOnID;
	TecUtilLockStart(*AddOnID);
	CSMGUILock();
	if (CSMGuiMultiListID != BADDIALOGID){
		LgIndex_t * SelectedItems;
		LgIndex_t NumSelected;
		TecGUIListGetSelectedItems(CSMGuiMultiListID, &SelectedItems, &NumSelected);
		CSMGUIDeleteCPLabels();
		if (NumSelected > 0){
			TecUtilDrawGraphics(FALSE);

			Set_pa TmpSet = nullptr;
			if (CSMGuiActiveZones != nullptr && !TecUtilSetIsEmpty(CSMGuiActiveZones) && TecUtilZoneGetActive(&TmpSet)) {
				TecUtilZoneSetActive(CSMGuiActiveZones, AssignOp_Equals);
			}

			double TextOffset = 0.1;
			for (int i = 0; i < NumSelected; ++i){
				char* LabelStr = TecGUIListGetString(CSMGuiMultiListID, SelectedItems[i]);

				vec3 Point = GetPointCoordsFromListItem(SelectedItems[i], CSMGuiPointSelectZoneNum);
				vec3 XYZGrid;
				TecUtilConvert3DPositionToGrid(Point[0], Point[1], Point[2], &XYZGrid[0], &XYZGrid[1], &XYZGrid[2]);

				Text_ID LabelID = TecUtilTextCreate(CoordSys_Grid, XYZGrid[0] + TextOffset, XYZGrid[1] + TextOffset, Units_Frame, 1.5, LabelStr);
				TecUtilTextBoxSetType(LabelID, TextBox_Filled);
				TecUtilTextBoxSetFillColor(LabelID, White_C);

				if (TecUtilTextIsValid(LabelID))
					CSMGuiPointLabelIDs.push_back(LabelID);

				TecUtilStringDealloc(&LabelStr);
			}
			TecUtilRedraw(TRUE);

			if (TmpSet == nullptr) {
				TecUtilDrawGraphics(TRUE);
			}
			else{
				TecUtilZoneSetActive(TmpSet, AssignOp_Equals);
				TecUtilSetDealloc(&TmpSet);
			}


		}
		TecUtilArrayDealloc((void **)&SelectedItems);
	}
	CSMGUIUnlock();
	TecUtilLockFinish(*AddOnID);
}


void CSMGuiSelectPointsWithinRadius(AddOn_pa *AddOnID) {
	if (AddOnID == nullptr) AddOnID = &CSMGuiAddOnID;
	TecUtilLockStart(*AddOnID);
	if (CSMGuiMultiListID != BADDIALOGID) {
		LgIndex_t * SelectedItems;
		LgIndex_t NumSelected;
		TecGUIListGetSelectedItems(CSMGuiMultiListID, &SelectedItems, &NumSelected);
		if (NumSelected > 0) {
			vector<vec3> SelectedPoints, AllPoints;

			SelectedPoints.reserve(NumSelected);
			for (int i = 0; i < NumSelected; ++i) {
				vec3 Point = GetPointCoordsFromListItem(SelectedItems[i], CSMGuiPointSelectZoneNum);
				SelectedPoints.push_back(Point);
			}

			int NumPoints = TecGUIListGetItemCount(CSMGuiMultiListID);
			AllPoints.reserve(NumPoints);
			for (int i = 1; i <= NumPoints; ++i) {
				vec3 Point = GetPointCoordsFromListItem(i, CSMGuiPointSelectZoneNum);
				AllPoints.push_back(Point);
			}

			// Get radius from user
			double Radius = -1.;
			while (Radius < 0.) {
				char * cStr;
				if (TecUtilDialogGetSimpleText("Select CPs around selected CPs within radius:", "4.0", &cStr)) {
					if (StringIsFloat(cStr) || StringIsInt(cStr)) {
						Radius = atof(cStr);
					}
					else {
						TecUtilDialogErrMsg("Please enter positive number");
					}
					TecUtilStringDealloc(&cStr);
				}
				else{
					Radius = -1;
					break;
				}
			}
			double RadiusSqr = Radius * Radius;

			// Now loop over points of associated CP zone to select those with
			// specified distance. (Would like to do all CPs, but this way
			// one can select only a particular type of CP, which is nice.
			// 
			vector<LgIndex_t> NewSelectedItems;
			for (int i = 0; i < AllPoints.size(); ++i){
				auto const & p1 = AllPoints[i];
				for (auto const & p2 : SelectedPoints){
					if (DistSqr(p1,p2) <= RadiusSqr){
						NewSelectedItems.push_back(i + 1);
						break;
					}
				}
			}

			TecGUIListSetSelectedItems(CSMGuiMultiListID, NewSelectedItems.data(), NewSelectedItems.size());
		}
		TecUtilArrayDealloc((void **)&SelectedItems);
	}
	TecUtilMouseSetMode(MouseButtonMode_Select);
	TecUtilLockFinish(*AddOnID);
}

void CSMGuiSelectPointsWithinRadiusButtonCB(){
	CSMGuiSelectPointsWithinRadius(&CSMGuiAddOnID);
	CSMGuiLabelSelectedPoints(&CSMGuiAddOnID);
}



void CSMGuiPointSelectMultiListCallback(const LgIndex_t* Val){
	TecUtilLockStart(CSMGuiAddOnID);
	if (CSMGuiMultiListIsPointSelect){
		CSMGUILock();
		CSMGuiLabelSelectedPoints();
		CSMGUIUnlock();
	}
	TecUtilLockFinish(CSMGuiAddOnID);
}

void CSMGuiMultiListDeselectAll(){
	vector<string> ListItems;
	int NumItems = TecGUIListGetItemCount(CSMGuiMultiListID);
	for (int i = 0; i < NumItems; ++i) ListItems.push_back(TecGUIListGetString(CSMGuiMultiListID, i + 1));
	TecGUIListDeleteAllItems(CSMGuiMultiListID);
	for (auto const & s : ListItems) TecGUIListAppendItem(CSMGuiMultiListID, s.c_str());
}

void CSMGuiPointSelectButtonCB();
void STDCALL CSMGuiPointSelectProbeCB(Boolean_t WasSuccessful,
	Boolean_t isNearestPoint,
	ArbParam_t ClientData){

	if (!isNearestPoint){
		TecUtilDialogErrMsg("You didn't hold the control key when you clicked...");
		return;
	}

	EntIndex_t ZoneNum = TecUtilProbeFieldGetZone();
	if (ZoneNum == TECUTILBADZONENUMBER){
		TecUtilDialogErrMsg("Failed to get zone number!");
		return;
	}

	EntIndex_t CPNum = TecUtilProbeGetPointIndex();

	if (ZoneNum != CSMGuiPointSelectZoneNum){
		bool PointFound = false;
		int CPTypeInd = std::find(CSMAuxData.CC.CPSubTypes.begin(), CSMAuxData.CC.CPSubTypes.end(), AuxDataZoneGetItem(ZoneNum, CSMAuxData.CC.ZoneSubType)) - CSMAuxData.CC.CPSubTypes.begin();
		string SearchString = CPNameList[CPTypeInd] + " " + to_string(CPNum);
		int ItemCount = TecGUIListGetItemCount(CSMGuiMultiListID);
		for (int i = 1; i <= ItemCount && !PointFound; ++i){
			if (SearchString == TecGUIListGetString(CSMGuiMultiListID, i)){
				CPNum = i;
				PointFound = true;
			}
		}
		if (!PointFound){
			TecUtilDialogErrMsg("Select a point from the CP zone");
			return;
		}
	}

	TecUtilLockStart(CSMGuiAddOnID);

	CSMGUILock();

	bool PointSelected = false;
	int NumSelected;
	int * SelectedNums;
	TecGUIListGetSelectedItems(CSMGuiMultiListID, &SelectedNums, &NumSelected);

	vector<int> NewSelectedNums;
	NewSelectedNums.reserve(NumSelected + 1);

	for (int i = 0; i < NumSelected; ++i){
		if (!PointSelected && CPNum == SelectedNums[i]) PointSelected = true;
		if (CPNum != SelectedNums[i]) NewSelectedNums.push_back(SelectedNums[i]);
	}

	if (!PointSelected){
		NewSelectedNums.push_back(CPNum);
// 		TecGUIListSetTopItemNum(CSMGuiMultiListID, CPNum);
	}
	CSMGuiMultiListDeselectAll();
	if (NewSelectedNums.size() > 0){
		TecGUIListSetSelectedItems(CSMGuiMultiListID, NewSelectedNums.data(), NewSelectedNums.size());
	}
	else{
		// Need to just rebuild the list, since Tecplot won't let you deselect everything.
	}
	CSMGuiLabelSelectedPoints();
	CSMGuiPointSelectButtonCB();
// 	TecUtilRedraw(TRUE);

	CSMGUILock();


	TecUtilLockFinish(CSMGuiAddOnID);
}

void CSMGuiPointSelectButtonCB()
{
	TecUtilLockStart(CSMGuiAddOnID);
	if (TecUtilDataSetIsAvailable() && TecUtilFrameGetPlotType() != PlotType_Sketch){
		if (CSMGuiProbeFirstInstall){
// 			CSMGUIDeleteCPLabels();
			CSMGuiProbeFirstInstall = false;
			if (CSMGuiActiveZones != nullptr)
				TecUtilSetDealloc(&CSMGuiActiveZones);
			TecUtilZoneGetActive(&CSMGuiActiveZones);
			CSMGuiMouseMode = TecUtilMouseGetCurrentMode();

			Set_pa ProbeZoneSet = TecUtilSetAlloc(FALSE);
			SetIndex_t Zone;
			TecUtilSetForEachMember(Zone, CSMGuiActiveZones)
			{
				if (AuxDataZoneItemMatches(Zone, CSMAuxData.CC.ZoneType, CSMAuxData.CC.ZoneTypeCPs))
					TecUtilSetAddMember(ProbeZoneSet, Zone, FALSE);
			}
			if (TecUtilSetGetMemberCount(ProbeZoneSet) <= 0)
				TecUtilSetAddMember(ProbeZoneSet, CSMGuiPointSelectZoneNum, FALSE);
			TecUtilZoneSetActive(ProbeZoneSet, AssignOp_Equals);
			TecUtilSetDealloc(&ProbeZoneSet);
		}

		ArgList_pa ProbeArgs = TecUtilArgListAlloc();
		TecUtilArgListAppendFunction(ProbeArgs, SV_CALLBACKFUNCTION, CSMGuiPointSelectProbeCB);
		TecUtilArgListAppendString(ProbeArgs, SV_STATUSLINETEXT, "Select CPs holding ctrl (or Command on Mac)");
		// 		TecUtilArgListAppendArbParam(ProbeArgs, SV_CLIENTDATA, reinterpret_cast<ArbParam_t>(PlaneClientData));
		if (!TecUtilProbeInstallCallbackX(ProbeArgs)){
			TecUtilDialogErrMsg("Failed to install probe callback function!");
		}
		TecUtilArgListDealloc(&ProbeArgs);
	}
	else
		TecUtilDialogErrMsg("To execute this add-on Tecplot must have\n"
		"a data set and the frame must be in XY,\n"
		"2D, or 3D.");
	TecUtilLockFinish(CSMGuiAddOnID);
}

void CSMGuiMultiListInvertSelectionButtonCallBack(){
	TecUtilLockStart(CSMGuiAddOnID);
	LgIndex_t* SelectedNums;
	LgIndex_t NumSelected, ItemCount = TecGUIListGetItemCount(CSMGuiMultiListID);
	TecGUIListGetSelectedItems(CSMGuiMultiListID, &SelectedNums, &NumSelected);

	vector<LgIndex_t> NewSelectedNums;
	NewSelectedNums.reserve(ItemCount - NumSelected);
	for (int i = 1; i <= ItemCount; ++i){
		bool found = false;
		for (int j = 0; j < NumSelected && !found; ++j) found = (SelectedNums[j] == i);
		if (!found) NewSelectedNums.push_back(i);
	}

	CSMGuiMultiListDeselectAll();

	if (NewSelectedNums.size() > 0)
		TecGUIListSetSelectedItems(CSMGuiMultiListID, NewSelectedNums.data(), NewSelectedNums.size());

	TecUtilArrayDealloc((void**)&SelectedNums);

	CSMGuiPointSelectMultiListCallback(nullptr);
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
// 		TecGUIListSetSelectedItem(CSMGuiMultiListID, 1);
		CSMGuiMultiListDeselectAll();
	}

	CSMGuiPointSelectMultiListCallback(nullptr);
	TecUtilLockFinish(CSMGuiAddOnID);
}

void CSMGUIDeleteCPLabels(AddOn_pa *AddOnID){
// 	if (CSMGuiPointLabelIDs.size() > 0){
// 		CSMGuiDeleteLabelsRunning = true;
		if (AddOnID == nullptr) AddOnID = &CSMGuiAddOnID;
		TecUtilLockStart(*AddOnID);
		MouseButtonMode_e MouseMode = TecUtilMouseGetCurrentMode();
		TecUtilDrawGraphics(FALSE);

		Text_ID Text;

		Text = TecUtilTextGetBase();
		while (Text != TECUTILBADID)
		{
			TecUtilTextDelete(Text);
			Text = TecUtilTextGetBase();
		}

// 		for (int i = 0; i < 2; ++i){
// 			TecUtilPickDeselectAll();
// 			for (auto const & l : CSMGuiPointLabelIDs) if (TecUtilTextIsValid(l)) TecUtilPickText(l);
// 
// 			if (TecUtilPickListGetCount() > 0){
// 				TecUtilPickClear();
// 			}
// 		}
// 		if (TecUtilPickListGetCount() > 0){
// 			CSMGuiPointLabelIDs.clear();
// 		}
// 		TecUtilDrawGraphics(TRUE);
// 		TecUtilRedraw(TRUE);
		if (TecUtilMouseIsValidMode(MouseMode))
			TecUtilMouseSetMode(MouseMode);
		TecUtilLockFinish(*AddOnID);
// 		CSMGuiDeleteLabelsRunning = false;
// 	}
}

void DialogCloseButtonCB(){
	TecUtilLockStart(CSMGuiAddOnID);
	CSMGuiDialogUp = false;
	CSMGuiSuccess = false;
	CSMGuiProbeFirstInstall = true;

	if (CSMGuiActiveZones != nullptr){
		if (!TecUtilSetIsEmpty(CSMGuiActiveZones)){
			SetValueReturnCode_e ReturnCode = TecUtilZoneSetActive(CSMGuiActiveZones, AssignOp_Equals);
			TecUtilSetDealloc(&CSMGuiActiveZones);
		}
		CSMGuiActiveZones = nullptr;
	}
	if (TecUtilMouseGetCurrentMode() == MouseButtonMode_Probe && CSMGuiMouseMode != MouseButtonMode_Invalid){
		if (CSMGuiMouseMode == MouseButtonMode_Probe) TecUtilMouseSetMode(MouseButtonMode_RotateRollerBall);
		TecUtilMouseSetMode(CSMGuiMouseMode);
	}

	CSMGuiMultiListID = BADDIALOGID;
	CSMGuiPointSelectZoneNum = -1;

	TecGUIDialogDrop(CSMGuiDialogManager);
	CSMGUIDeleteCPLabels();
	TecUtilDrawGraphics(CSMGuiAutoRedrawIsActive);
	TecUtilLockFinish(CSMGuiAddOnID);
}

void DialogOKButtonCB(){

	CSMGuiSuccess = true;
	CSMGuiDialogUp = false;
	CSMGuiProbeFirstInstall = true;

	TecUtilLockStart(CSMGuiAddOnID);

	if (CSMGuiActiveZones != nullptr){
		if (!TecUtilSetIsEmpty(CSMGuiActiveZones)){
			SetValueReturnCode_e ReturnCode = TecUtilZoneSetActive(CSMGuiActiveZones, AssignOp_Equals);
			TecUtilSetDealloc(&CSMGuiActiveZones);
		}
		CSMGuiActiveZones = nullptr;
	}
	if (TecUtilMouseGetCurrentMode() == MouseButtonMode_Probe && CSMGuiMouseMode != MouseButtonMode_Invalid){
		if (CSMGuiMouseMode == MouseButtonMode_Probe) TecUtilMouseSetMode(MouseButtonMode_RotateRollerBall);
		TecUtilMouseSetMode(CSMGuiMouseMode);
	}

	CSMGuiMultiListID = BADDIALOGID;
	CSMGuiPointSelectZoneNum = -1;

	TecGUIDialogDrop(CSMGuiDialogManager);
	CSMGUIDeleteCPLabels();

	for (auto & f : CSMGuiFields){
		GuiFieldType_e t = f.GetType();

		if (t <= Gui_VarSelect) f.SetReturnInt(TecGUIOptionMenuGet(f.GetID()));
		else if (t <= Gui_Int) f.SetReturnInt(atoi(TecGUITextFieldGetString(f.GetID())));
		else if (t <= Gui_Double) f.SetReturnDouble(atof(TecGUITextFieldGetString(f.GetID())));
		else if (t <= Gui_ToggleEnable) f.SetReturnBool(TecGUIToggleGet(f.GetID()));
		else if (t <= Gui_ZonePointSelectMulti){
			LgIndex_t* SelectedNums;
			LgIndex_t NumSelected;
			vector<int> IntVec;
			vector<string> StringVec;

			TecGUIListGetSelectedItems(f.GetID(), &SelectedNums, &NumSelected);
			for (int i = 0; i < NumSelected; ++i){
				IntVec.push_back(SelectedNums[i]);
				StringVec.push_back(TecGUIListGetString(f.GetID(), SelectedNums[i]));
				/*
				 *	Remove leading number and period if present
				 */
				vector<string> tmpStr = SplitString(StringVec.back(), ". ");
				if (tmpStr.size() > 1 && StringIsInt(tmpStr[0])){
					StringVec.back() = VectorToString(vector<string>(tmpStr.begin() + 1, tmpStr.end()), ". ");
				}
			}
			f.SetReturnStringVec(StringVec);
			f.SetReturnIntVec(IntVec);
			TecUtilArrayDealloc((void**)&SelectedNums);

			if (t == Gui_ZonePointSelectMulti) f.SetReturnInt(TecGUIOptionMenuGet(f.GetID(1)));
		}
		else if (t <= Gui_String) f.SetReturnString(TecGUITextFieldGetString(f.GetID()));
		else if (t <= Gui_Radio) f.SetReturnInt(TecGUIRadioBoxGetToggle(f.GetID()));
	}

	TecUtilDrawGraphics(CSMGuiAutoRedrawIsActive);

	TecUtilLockFinish(CSMGuiAddOnID);

	CurrentReturnFunc(CSMGuiSuccess, CSMGuiFields, CSMPassthroughFields);
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
	vector<int> CPTypeCounts(CPNameList.size(), 1);
	for (int i = 0; i < NumPoints; ++i){
		if (ZoneIsCP){
			int CPType = TecUtilDataValueGetByRef(CPTypeRef, i + 1);
			int CPInd = VectorGetElementNum(CPTypeList, (CPType_e)CPType);
			LabelBase = CPNameList[CPInd] + " " + to_string(CPTypeCounts[CPInd]++);
		}
		TecGUIListAppendItem(CSMGuiMultiListID, ZoneIsCP ? LabelBase.c_str() : to_string(i + 1).c_str());
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

/*
 * In order to distinguish between the buttons of one gui field and another,
 * make multiple (5) copies of void funtion pointers whose sole task is to set
 * the global multilist id for the correct field.
 * Need one such function for each of the 5 possible ids, for each function that
 * uses the global multilist id.
 * How to declare array of function pointers: https://docs.microsoft.com/en-us/troubleshoot/cpp/declare-pointers-to-functions
 */
void CSMGuiMultiListInvertSelectionButtonCallBack0(){
	CSMGuiMultiListID = CSMGuiMultiListIDs[0];
	CSMGuiMultiListIsPointSelect = CSMGuiMultiListIsPointSelects[0];
	CSMGuiMultiListInvertSelectionButtonCallBack();
}
void CSMGuiMultiListInvertSelectionButtonCallBack1() {
	CSMGuiMultiListID = CSMGuiMultiListIDs[1];
	CSMGuiMultiListIsPointSelect = CSMGuiMultiListIsPointSelects[1];
	CSMGuiMultiListInvertSelectionButtonCallBack();
}
void CSMGuiMultiListInvertSelectionButtonCallBack2() {
	CSMGuiMultiListID = CSMGuiMultiListIDs[2];
	CSMGuiMultiListIsPointSelect = CSMGuiMultiListIsPointSelects[2];
	CSMGuiMultiListInvertSelectionButtonCallBack();
}
void CSMGuiMultiListInvertSelectionButtonCallBack3() {
	CSMGuiMultiListID = CSMGuiMultiListIDs[3];
	CSMGuiMultiListIsPointSelect = CSMGuiMultiListIsPointSelects[3];
	CSMGuiMultiListInvertSelectionButtonCallBack();
}
void CSMGuiMultiListInvertSelectionButtonCallBack4() {
	CSMGuiMultiListID = CSMGuiMultiListIDs[4];
	CSMGuiMultiListIsPointSelect = CSMGuiMultiListIsPointSelects[4];
	CSMGuiMultiListInvertSelectionButtonCallBack();
}
void (*CSMGuiMultiListInvertSelectionButtonCallBacks[])() = { 
	CSMGuiMultiListInvertSelectionButtonCallBack0,
	CSMGuiMultiListInvertSelectionButtonCallBack1,
	CSMGuiMultiListInvertSelectionButtonCallBack2,
	CSMGuiMultiListInvertSelectionButtonCallBack3,
	CSMGuiMultiListInvertSelectionButtonCallBack4
};

void CSMGuiMultiListSelectAllButtonCallBack0() {
	CSMGuiMultiListID = CSMGuiMultiListIDs[0];
	CSMGuiMultiListIsPointSelect = CSMGuiMultiListIsPointSelects[0];
	CSMGuiMultiListSelectAllButtonCallBack();
}
void CSMGuiMultiListSelectAllButtonCallBack1() {
	CSMGuiMultiListID = CSMGuiMultiListIDs[1];
	CSMGuiMultiListIsPointSelect = CSMGuiMultiListIsPointSelects[1];
	CSMGuiMultiListSelectAllButtonCallBack();
}
void CSMGuiMultiListSelectAllButtonCallBack2() {
	CSMGuiMultiListID = CSMGuiMultiListIDs[2];
	CSMGuiMultiListIsPointSelect = CSMGuiMultiListIsPointSelects[2];
	CSMGuiMultiListSelectAllButtonCallBack();
}
void CSMGuiMultiListSelectAllButtonCallBack3() {
	CSMGuiMultiListID = CSMGuiMultiListIDs[3];
	CSMGuiMultiListIsPointSelect = CSMGuiMultiListIsPointSelects[3];
	CSMGuiMultiListSelectAllButtonCallBack();
}
void CSMGuiMultiListSelectAllButtonCallBack4() {
	CSMGuiMultiListID = CSMGuiMultiListIDs[4];
	CSMGuiMultiListIsPointSelect = CSMGuiMultiListIsPointSelects[4];
	CSMGuiMultiListSelectAllButtonCallBack();
}
void(*CSMGuiMultiListSelectAllButtonCallBacks[])() = {
	CSMGuiMultiListSelectAllButtonCallBack0,
	CSMGuiMultiListSelectAllButtonCallBack1,
	CSMGuiMultiListSelectAllButtonCallBack2,
	CSMGuiMultiListSelectAllButtonCallBack3,
	CSMGuiMultiListSelectAllButtonCallBack4
};

void CSMGuiPointSelectButtonCB0() {
	CSMGuiMultiListID = CSMGuiMultiListIDs[0];
	CSMGuiMultiListIsPointSelect = CSMGuiMultiListIsPointSelects[0];
	CSMGuiPointSelectButtonCB();
}
void CSMGuiPointSelectButtonCB1() {
	CSMGuiMultiListID = CSMGuiMultiListIDs[1];
	CSMGuiMultiListIsPointSelect = CSMGuiMultiListIsPointSelects[1];
	CSMGuiPointSelectButtonCB();
}
void CSMGuiPointSelectButtonCB2() {
	CSMGuiMultiListID = CSMGuiMultiListIDs[2];
	CSMGuiMultiListIsPointSelect = CSMGuiMultiListIsPointSelects[2];
	CSMGuiPointSelectButtonCB();
}
void CSMGuiPointSelectButtonCB3() {
	CSMGuiMultiListID = CSMGuiMultiListIDs[3];
	CSMGuiMultiListIsPointSelect = CSMGuiMultiListIsPointSelects[3];
	CSMGuiPointSelectButtonCB();
}
void CSMGuiPointSelectButtonCB4() {
	CSMGuiMultiListID = CSMGuiMultiListIDs[4];
	CSMGuiMultiListIsPointSelect = CSMGuiMultiListIsPointSelects[4];
	CSMGuiPointSelectButtonCB();
}
void(*CSMGuiPointSelectButtonCBs[])() = {
	CSMGuiPointSelectButtonCB0,
	CSMGuiPointSelectButtonCB1,
	CSMGuiPointSelectButtonCB2,
	CSMGuiPointSelectButtonCB3,
	CSMGuiPointSelectButtonCB4
};

void CSMGuiPointSelectMultiListCallback0(const LgIndex_t* Val) {
	CSMGuiMultiListID = CSMGuiMultiListIDs[0];
	CSMGuiMultiListIsPointSelect = CSMGuiMultiListIsPointSelects[0];
	CSMGuiPointSelectMultiListCallback(Val);
}
void CSMGuiPointSelectMultiListCallback1(const LgIndex_t* Val) {
	CSMGuiMultiListID = CSMGuiMultiListIDs[1];
	CSMGuiMultiListIsPointSelect = CSMGuiMultiListIsPointSelects[1];
	CSMGuiPointSelectMultiListCallback(Val);
}
void CSMGuiPointSelectMultiListCallback2(const LgIndex_t* Val) {
	CSMGuiMultiListID = CSMGuiMultiListIDs[2];
	CSMGuiMultiListIsPointSelect = CSMGuiMultiListIsPointSelects[2];
	CSMGuiPointSelectMultiListCallback(Val);
}
void CSMGuiPointSelectMultiListCallback3(const LgIndex_t* Val) {
	CSMGuiMultiListID = CSMGuiMultiListIDs[3];
	CSMGuiMultiListIsPointSelect = CSMGuiMultiListIsPointSelects[3];
	CSMGuiPointSelectMultiListCallback(Val);
}
void CSMGuiPointSelectMultiListCallback4(const LgIndex_t* Val) {
	CSMGuiMultiListID = CSMGuiMultiListIDs[4];
	CSMGuiMultiListIsPointSelect = CSMGuiMultiListIsPointSelects[4];
	CSMGuiPointSelectMultiListCallback(Val);
}
void(*CSMGuiPointSelectMultiListCallbacks[])(const LgIndex_t* Val) = {
	CSMGuiPointSelectMultiListCallback0,
	CSMGuiPointSelectMultiListCallback1,
	CSMGuiPointSelectMultiListCallback2,
	CSMGuiPointSelectMultiListCallback3,
	CSMGuiPointSelectMultiListCallback4
};

void CSMGuiPointSelectOptionCallback0(const LgIndex_t* Val) {
	CSMGuiMultiListID = CSMGuiMultiListIDs[0];
	CSMGuiMultiListIsPointSelect = CSMGuiMultiListIsPointSelects[0];
	CSMGuiPointSelectOptionCallback(Val);
}
void CSMGuiPointSelectOptionCallback1(const LgIndex_t* Val) {
	CSMGuiMultiListID = CSMGuiMultiListIDs[1];
	CSMGuiMultiListIsPointSelect = CSMGuiMultiListIsPointSelects[1];
	CSMGuiPointSelectOptionCallback(Val);
}
void CSMGuiPointSelectOptionCallback2(const LgIndex_t* Val) {
	CSMGuiMultiListID = CSMGuiMultiListIDs[2];
	CSMGuiMultiListIsPointSelect = CSMGuiMultiListIsPointSelects[2];
	CSMGuiPointSelectOptionCallback(Val);
}
void CSMGuiPointSelectOptionCallback3(const LgIndex_t* Val) {
	CSMGuiMultiListID = CSMGuiMultiListIDs[3];
	CSMGuiMultiListIsPointSelect = CSMGuiMultiListIsPointSelects[3];
	CSMGuiPointSelectOptionCallback(Val);
}
void CSMGuiPointSelectOptionCallback4(const LgIndex_t* Val) {
	CSMGuiMultiListID = CSMGuiMultiListIDs[4];
	CSMGuiMultiListIsPointSelect = CSMGuiMultiListIsPointSelects[4];
	CSMGuiPointSelectOptionCallback(Val);
}
void(*CSMGuiPointSelectOptionCallbacks[])(const LgIndex_t* Val) = {
	CSMGuiPointSelectOptionCallback0,
	CSMGuiPointSelectOptionCallback1,
	CSMGuiPointSelectOptionCallback2,
	CSMGuiPointSelectOptionCallback3,
	CSMGuiPointSelectOptionCallback4
};

void STDCALL CSMGuiPointSelectProbeCB0(Boolean_t WasSuccessful,
	Boolean_t isNearestPoint,
	ArbParam_t ClientData) {
	CSMGuiMultiListID = CSMGuiMultiListIDs[0];
	CSMGuiMultiListIsPointSelect = CSMGuiMultiListIsPointSelects[0];
	CSMGuiPointSelectProbeCB(WasSuccessful, isNearestPoint, ClientData);
}
void STDCALL CSMGuiPointSelectProbeCB1(Boolean_t WasSuccessful,
	Boolean_t isNearestPoint,
	ArbParam_t ClientData) {
	CSMGuiMultiListID = CSMGuiMultiListIDs[1];
	CSMGuiMultiListIsPointSelect = CSMGuiMultiListIsPointSelects[1];
	CSMGuiPointSelectProbeCB(WasSuccessful, isNearestPoint, ClientData);
}
void STDCALL CSMGuiPointSelectProbeCB2(Boolean_t WasSuccessful,
	Boolean_t isNearestPoint,
	ArbParam_t ClientData) {
	CSMGuiMultiListID = CSMGuiMultiListIDs[2];
	CSMGuiMultiListIsPointSelect = CSMGuiMultiListIsPointSelects[2];
	CSMGuiPointSelectProbeCB(WasSuccessful, isNearestPoint, ClientData);
}
void STDCALL CSMGuiPointSelectProbeCB3(Boolean_t WasSuccessful,
	Boolean_t isNearestPoint,
	ArbParam_t ClientData) {
	CSMGuiMultiListID = CSMGuiMultiListIDs[3];
	CSMGuiMultiListIsPointSelect = CSMGuiMultiListIsPointSelects[3];
	CSMGuiPointSelectProbeCB(WasSuccessful, isNearestPoint, ClientData);
}
void STDCALL CSMGuiPointSelectProbeCB4(Boolean_t WasSuccessful,
	Boolean_t isNearestPoint,
	ArbParam_t ClientData) {
	CSMGuiMultiListID = CSMGuiMultiListIDs[4];
	CSMGuiMultiListIsPointSelect = CSMGuiMultiListIsPointSelects[4];
	CSMGuiPointSelectProbeCB(WasSuccessful, isNearestPoint, ClientData);
}
void(STDCALL *CSMGuiPointSelectProbeCBs[])(Boolean_t WasSuccessful,
	Boolean_t isNearestPoint,
	ArbParam_t ClientData) = {
	CSMGuiPointSelectProbeCB0,
	CSMGuiPointSelectProbeCB1,
	CSMGuiPointSelectProbeCB2,
	CSMGuiPointSelectProbeCB3,
	CSMGuiPointSelectProbeCB4
};

void CSMGuiSelectPointsWithinRadiusButtonCB0() {
	CSMGuiMultiListID = CSMGuiMultiListIDs[0];
	CSMGuiMultiListIsPointSelect = CSMGuiMultiListIsPointSelects[0];
	CSMGuiSelectPointsWithinRadiusButtonCB();
}
void CSMGuiSelectPointsWithinRadiusButtonCB1() {
	CSMGuiMultiListID = CSMGuiMultiListIDs[1];
	CSMGuiMultiListIsPointSelect = CSMGuiMultiListIsPointSelects[1];
	CSMGuiSelectPointsWithinRadiusButtonCB();
}
void CSMGuiSelectPointsWithinRadiusButtonCB2() {
	CSMGuiMultiListID = CSMGuiMultiListIDs[2];
	CSMGuiMultiListIsPointSelect = CSMGuiMultiListIsPointSelects[2];
	CSMGuiSelectPointsWithinRadiusButtonCB();
}
void CSMGuiSelectPointsWithinRadiusButtonCB3() {
	CSMGuiMultiListID = CSMGuiMultiListIDs[3];
	CSMGuiMultiListIsPointSelect = CSMGuiMultiListIsPointSelects[3];
	CSMGuiSelectPointsWithinRadiusButtonCB();
}
void CSMGuiSelectPointsWithinRadiusButtonCB4() {
	CSMGuiMultiListID = CSMGuiMultiListIDs[4];
	CSMGuiMultiListIsPointSelect = CSMGuiMultiListIsPointSelects[4];
	CSMGuiSelectPointsWithinRadiusButtonCB();
}
void(*CSMGuiSelectPointsWithinRadiusButtonCBs[])() = {
	CSMGuiSelectPointsWithinRadiusButtonCB0,
	CSMGuiSelectPointsWithinRadiusButtonCB1,
	CSMGuiSelectPointsWithinRadiusButtonCB2,
	CSMGuiSelectPointsWithinRadiusButtonCB3,
	CSMGuiSelectPointsWithinRadiusButtonCB4
};


void CSMLaunchGui(string const & Title, 
	vector<GuiField_c> const & Fields)
{
	TecUtilPleaseWait("Loading dialog...", TRUE);

	CSMGuiAutoRedrawIsActive = TecUtilAutoRedrawIsActive();

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
		MaxZoneVarWidth = 15;

	vector<vector<int> > FieldZoneVarNums(CSMGuiFields.size());

	for (auto const & f : CSMGuiFields) {
		if (f.GetType() == Gui_Radio){
			NumFields[int(Gui_Radio)] += SplitString(f.GetSearchString(), ",").size();
		}
		else NumFields[int(f.GetType())]++;
		if (f.GetType() < Gui_Toggle) MaxLabelWidth = MAX(MaxLabelWidth, int(f.GetLabel().length()));
		else MaxOnlyLabelWidth = MAX(MaxOnlyLabelWidth, int(f.GetLabel().length()));
	}

	vector<vector<string> > ZoneVarList(2);
	vector<bool> ZoneVarCommaListValid(2, true);
	vector<string> ZoneVarCommaList(2); // comma-delimited list of zones/vars for option menus
	vector<int> NumZonesVars = { TecUtilDataSetGetNumZones(), TecUtilDataSetGetNumVars() };

	for (int i = 0; i < 2; ++i) for (auto const & t : ZoneVarFields[i]) if (NumFields[int(t)]){ // Get full list of zone/var names
		char *cStr;
		for (int j = 1; j <= NumZonesVars[i]; ++j){
			if ((i == 0 && TecUtilZoneGetName(j, &cStr)) || TecUtilVarGetName(j, &cStr)){
				ZoneVarList[i].push_back(to_string(j) + ". " + cStr);
				MaxZoneVarWidth = MAX(MaxZoneVarWidth, int(ZoneVarList[i].back().length()));
				ZoneVarCommaListValid[i] = (ZoneVarCommaListValid[i] && SplitString(ZoneVarList[i].back(), ",").size() == 1);
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
			for (auto const & s : SplitString(CSMGuiFields[f].GetSearchString(), ",")){
				FieldZoneVarNums[f].push_back(ZoneNumByName(s, false, true));
			}
		}
		else if (t == Gui_VarSelect || t == Gui_VarSelectInt || t == Gui_VarSelectMulti){
			for (auto const & s : SplitString(CSMGuiFields[f].GetSearchString(), ",")){
				FieldZoneVarNums[f].push_back(VarNumByName(s, true));
			}
		}
	}

	int  MultiListNumLines = ListNumLines;
	if (CSMGuiFields.size() < 8) MultiListNumLines *= 2;

	if (NumFields[int(Gui_ZonePointSelectMulti)])
		MaxLabelWidth = MAX(MaxLabelWidth, 20);
	else
		MaxLabelWidth = MAX(MaxLabelWidth, 18);

	int X = HorzSpacing, Y = LineHeight;
	int W = CharWidth * MAX((MaxLabelWidth + MaxZoneVarWidth), MaxOnlyLabelWidth) + HorzSpacing;
	int H = 2 * LineHeight + VertSpacing;

	int xTmp = X - HorzSpacing + CharWidth * MaxLabelWidth;
	int wTmp = CharWidth * MaxZoneVarWidth;

	for (int t = 0; t < NumFields.size(); ++t){
		if (GuiFieldType_e(t) <= Gui_ToggleEnable || GuiFieldType_e(t) == Gui_Label) H += (LineHeight + VertSpacing) * NumFields[t];
		else if (GuiFieldType_e(t) <= Gui_ZonePointSelectMulti) H += (LineHeight * (MultiListNumLines + 2.0)) * NumFields[t];
		else if (GuiFieldType_e(t) == Gui_Radio) H += LineHeight * NumFields[t];
			// For radio selection, need to count the lines in each.
// 			for (const auto f : CSMGuiFields)
// 				if (f.GetType() == Gui_Radio) 
// 					H += (LineHeight + VertSpacing) * SplitString(f.GetSearchString(), ",").size();
		else H += (2 * VertSpacing) * NumFields[t];
	}

	CSMGuiDialogManager = TecGUIDialogCreateModeless(MAINDIALOGID, W, H, Title.c_str(), nullptr, DialogCloseButtonCB, nullptr);

	int fNum = 0;
	int ListIDInd = 0;

	for (auto & f : CSMGuiFields){
		GuiFieldType_e t = f.GetType();
		if (t < Gui_VertSep){
			TecGUILabelAdd(CSMGuiDialogManager, X, Y + VertSpacing * 0.5, f.GetLabel().c_str());

			if (t <= Gui_VarSelect){ // Gui Field is an option menu for selecting zone/variable
				if (ZoneVarCommaListValid[int(t)])
					f.SetID(TecGUIOptionMenuAdd(CSMGuiDialogManager, xTmp, Y, wTmp, LineHeight, ZoneVarCommaList[int(t)].c_str(), VoidCallbackDoNothing));
				else{
					f.SetID(TecGUIOptionMenuAdd(CSMGuiDialogManager, xTmp, Y, wTmp, LineHeight, ZoneVarList[int(t)][0].c_str(), VoidCallbackDoNothing));
					for (int i = 1; i < ZoneVarList[int(t)].size(); ++i) TecGUIOptionMenuAppendItem(f.GetID(), ZoneVarList[int(t)][i].c_str());
				}
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
				if (t < Gui_Int) TecGUITextFieldSetLgIndex(f.GetID(), MAX(1, FieldZoneVarNums[fNum].back()), FALSE);
				else if (t < Gui_Double) TecGUITextFieldSetLgIndex(f.GetID(), stoi(f.GetSearchString()), FALSE);
				else TecGUITextFieldSetDouble(f.GetID(), stof(f.GetSearchString()), "%f");
			}
			else if (t <= Gui_ToggleEnable){
				int w = f.GetLabel().length() * CharWidth;
				f.SetID(TecGUIToggleAdd(CSMGuiDialogManager, X, Y + VertSpacing * 0.5, w, LineHeight, f.GetLabel().c_str(), VoidCallbackDoNothing));
			}
			else if (t <= Gui_ZonePointSelectMulti){
				if (ListIDInd >= CSMGuiMultiListMaxNumMultiListFields){
					TecUtilDialogErrMsg("Too many multilist fields used in dialog. Max is 5");
					continue;
				}
				vector<string> ListItems;
				TecGUIIntCallback_pf fCB;
				if (t == Gui_ZonePointSelectMulti){
					TecUtilWorkAreaSuspend(TRUE);
					TecUtilDrawGraphics(FALSE);
					if (ZoneVarCommaListValid[0])
						f.SetID(TecGUIOptionMenuAdd(CSMGuiDialogManager, xTmp, Y, wTmp, LineHeight, ZoneVarCommaList[0].c_str(), CSMGuiPointSelectOptionCallbacks[ListIDInd]), 1);
					else{
						f.SetID(TecGUIOptionMenuAdd(CSMGuiDialogManager, xTmp, Y, wTmp, LineHeight, StringReplaceSubString(ZoneVarList[0][0], ",", ".").c_str(), CSMGuiPointSelectOptionCallbacks[ListIDInd]), 1);
						for (int i = 1; i < ZoneVarList[0].size(); ++i) 
							TecGUIOptionMenuAppendItem(f.GetID(), StringReplaceSubString(ZoneVarList[0][i], ",", ".").c_str());
					}
					TecGUIOptionMenuSet(f.GetID(1), MAX(1, FieldZoneVarNums[fNum].back()));
					fCB = CSMGuiPointSelectMultiListCallbacks[ListIDInd];
					CSMGuiMultiListIsPointSelects[ListIDInd] = true;
					Y += LineHeight + VertSpacing;
				}
				else{
					ListItems = ZoneVarList[int(t - Gui_ZoneSelectMulti)]; // This returns the zone list for ZoneSelectMulti and the var list for VarSelectMulti
					CSMGuiMultiListIsPointSelects[ListIDInd] = false;
					fCB = IntCallbackDoNothing;
				}
				f.SetID(TecGUIListAdd(CSMGuiDialogManager, xTmp, Y, wTmp, MultiListNumLines * LineHeight, TRUE, fCB));
				CSMGuiMultiListIDs[ListIDInd] = f.GetID();
				Y += LineHeight;
				if (t == Gui_ZonePointSelectMulti){
					TecGUIButtonAdd(CSMGuiDialogManager, X, Y - LineHeight * 0.1, CharWidth * 16, LineHeight * 1.2, "Select with mouse", CSMGuiPointSelectButtonCBs[ListIDInd]);
					Y += LineHeight * 1.5;
					TecGUIButtonAdd(CSMGuiDialogManager, X, Y - LineHeight * 0.1, CharWidth * 16, LineHeight * 1.2, "Select within radius", CSMGuiSelectPointsWithinRadiusButtonCBs[ListIDInd]);
					Y += LineHeight * 1.5;
				}
				TecGUIButtonAdd(CSMGuiDialogManager, X, Y - LineHeight * 0.1, CharWidth * 16, LineHeight * 1.2, "Invert selection", CSMGuiMultiListInvertSelectionButtonCallBacks[ListIDInd]);
				Y += LineHeight * 1.5;
				TecGUIButtonAdd(CSMGuiDialogManager, X, Y - LineHeight * 0.1, CharWidth * 16, LineHeight * 1.2, "Toggle all", CSMGuiMultiListSelectAllButtonCallBacks[ListIDInd]);
				if (t == Gui_ZonePointSelectMulti){
					int Val = MAX(1, FieldZoneVarNums[fNum].back());
					CSMGuiPointSelectOptionCallbacks[ListIDInd](&Val);
				}
				else for (auto const & s : ListItems) TecGUIListAppendItem(f.GetID(), s.c_str());
				if (t != Gui_ZonePointSelectMulti && FieldZoneVarNums[fNum].size() > 0){
					vector<LgIndex_t> SelectNums;
					int MinNum = INT_MAX;
					for (auto const & i : FieldZoneVarNums[fNum]) if (i > 0){
						SelectNums.push_back(i);
						MinNum = MIN(MinNum, i);
					}
					if (SelectNums.size() > 0){
						TecGUIListSetSelectedItems(f.GetID(), SelectNums.data(), SelectNums.size());
						TecGUIListSetTopItemNum(f.GetID(), MinNum);
					}
				}
				Y += LineHeight * (MultiListNumLines - 2);
				Y -= LineHeight * 1.5;
				if (t == Gui_ZonePointSelectMulti) Y -= LineHeight * 2.5;

				ListIDInd++;
			}
			else if (t <= Gui_String){
				f.SetID(TecGUITextFieldAdd(CSMGuiDialogManager, xTmp, Y, wTmp, LineHeight, TextCallbackDoNothing));
				TecGUITextFieldSetString(f.GetID(), f.GetSearchString().c_str());
			}
			else if (t <= Gui_Radio){
				vector<string> Choices = SplitString(f.GetSearchString(), ",");
				vector<const char*> ChoicesCStr(5);
				for (int i = 0; i < 5; ++i) ChoicesCStr[i] = (i < Choices.size() ? Choices[i].c_str() : nullptr);

				f.SetID(TecGUIRadioBoxAdd(CSMGuiDialogManager, xTmp, Y, wTmp, (LineHeight + VertSpacing) * Choices.size(),
					ChoicesCStr[0], ChoicesCStr[1], ChoicesCStr[2], ChoicesCStr[3], ChoicesCStr[4], IntCallbackDoNothing));
				Y += (LineHeight) * Choices.size();
			}

			Y += LineHeight + VertSpacing;
		}
		else if (t == Gui_VertSep){
			f.SetID(TecGUIVertSeparatorAdd(CSMGuiDialogManager, X, Y, VertSpacing));
			Y += VertSpacing * 2;
		}
		else if (t == Gui_Label) {
			TecGUILabelAdd(CSMGuiDialogManager, X, Y + VertSpacing * 0.5, f.GetLabel().c_str());
			Y += LineHeight + VertSpacing;
		}
		fNum++;
	}

	int OkButton = TecGUIButtonAdd(CSMGuiDialogManager, X + HorzSpacing, Y - LineHeight * 0.1, W - (4 * HorzSpacing), LineHeight * 1.2, "OK", DialogOKButtonCB);
	TecGUIButtonSetDefault(CSMGuiDialogManager, OkButton);

	// Now check that Gui_ToggleEnable fields are satisfied
	if (NumFields[int(Gui_ToggleEnable)]){
		for (auto const & f : CSMGuiFields){
			if (f.GetType() == Gui_ToggleEnable){
				vector<string> DependentFields = SplitString(f.GetSearchString(), ",");
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
						for (auto const & i : FieldZoneVarNums[fi]) NumHits += int(i > 0);
						IsEnabled = (NumHits > 0);
					}
				}
				TecGUIToggleSet(f.GetID(), IsEnabled);
			}
		}
	}

	// Set toggle default values based on toggle Val "0" or blank for false and "1" for true
	if (NumFields[int(Gui_Toggle)]){
		for (auto const & f : CSMGuiFields){
			if (f.GetType() == Gui_Toggle){
				TecGUIToggleSet(f.GetID(), f.GetSearchString().length() > 0 && f.GetSearchString() != "0");
			}
		}
	}

	TecUtilPleaseWait("Loading dialog...", FALSE);

	CSMGuiDialogUp = true;
	TecGUIDialogLaunch(CSMGuiDialogManager);

	// #ifdef MSWIN
	// 	RECT desktop;
	// 	HWND const hDesktop = GetDesktopWindow();
	// 	GetWindowRect(hDesktop, &desktop);
	// 
	// 	desktop.right; // screen width
	// 	desktop.bottom; // screen height
	// #else
	// 
	// #endif
}

void CSMGui(string const & Title, 
	vector<GuiField_c> const & Fields, 
	CSMGuiReturnFunc_pf ReturnFunc, 
	AddOn_pa const & InputAddOnID,
	vector<GuiField_c> const PassthroughFields)
{
	CurrentReturnFunc = ReturnFunc;
	CSMGuiAddOnID = InputAddOnID;
	CSMPassthroughFields = PassthroughFields;
	CSMLaunchGui(Title, Fields);
}

; void CSMGUILock() {
	TecUtilDrawGraphics(FALSE);
	TecUtilInterfaceSuspend(TRUE);
	TecUtilWorkAreaSuspend(TRUE);
}

void CSMGUIUnlock(){
	TecUtilDrawGraphics(TRUE);
	TecUtilInterfaceSuspend(FALSE);
	TecUtilWorkAreaSuspend(FALSE);
}