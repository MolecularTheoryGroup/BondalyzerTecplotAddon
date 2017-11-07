#pragma once
#ifndef CSMGUI_H_
#define CSMGUI_H_

#ifdef MSWIN
#include <Windows.h>
#endif

using std::vector;
using std::string;

enum GuiFieldType_e{ // comment describes gui field and specifies what should be provided to GuiField_c()'s Val argument
	Gui_ZoneSelect,				// Option menu to select zone: zone name
	Gui_VarSelect,				// Option menu to select var: var name

	Gui_ZoneSelectInt,			// Text input field to select zone: zone name
	Gui_VarSelectInt,			// Text input field to select var: var name
	Gui_Int,					// Text input to specity an integer: default value
	Gui_Double,					// Text input field to specify double value: default value

	Gui_Toggle,					// Toggle (checkbox): "1" for true (checked) or "0" for false (unchecked)
	Gui_ToggleEnable,			// Toggle (checkbox): list of gui field numbers that must have their zone/var number(s) found successfully in order to autoenable the toggle

	Gui_ZoneSelectMulti,		// Multi-select list to select zones: comma-delimited list of zone names
	Gui_VarSelectMulti,			// Multi-select list to select vars: comma-delimited list of var names

	Gui_ZonePointSelectMulti,	// Multi-select list to select points in a zone:
	//	Option menu is placed above multi-select box to specify zone
	//	Val is zone name to load in option menu
	//	** YOU CAN ONLY HAVE ONE IN A SINGLE DIALOG DUE TO LIMITATIONS OF THE LIST CALLBACK FUNCTION!!
	
	Gui_String,					// User-provided string

	Gui_Radio,					// Radio box: comma-delimited list of up to 5 options (only first 5 will be used if > 5 provided)

	Gui_VertSep,				// specify neither Label or Val

	Gui_Label,					// specify Label but no Val

	Gui_Invalid
};

class GuiField_c
{
public:
	GuiField_c();
	GuiField_c(const GuiFieldType_e & Type,
		const string & Label = "",
		const string & Val = "",
		const vector<void*> & CallbackFuntions = vector<void*>());
	~GuiField_c();

	const GuiFieldType_e GetType() const { return Type_m; }
	const string GetLabel() const { return Label_m; }
	const string GetSearchString() const { return InputVal_m; }
	const int GetID(const int & IDNum = 0) const { return FieldID_m[IDNum]; }

	const int GetReturnInt() const;
	const double GetReturnDouble() const;
	const string GetReturnString() const;
	const bool GetReturnBool() const;
	const vector<int> GetReturnIntVec() const;
	const vector<string> GetReturnStringVec() const;

	void SetSearchString(const string & s){ InputVal_m = s; }
	void AppendSearchString(const string & s){ InputVal_m += s; }
	void SetID(const int FieldID, const int & IDNum = 0){ FieldID_m[IDNum] = FieldID; }

	void SetReturnInt(const int & Val);
	void SetReturnDouble(const double & Val);
	void SetReturnString(const string & Val);
	void SetReturnBool(const bool & Val);
	void SetReturnIntVec(const vector<int> & Val);
	void SetReturnStringVec(const vector<string> & Val);

private:
	void CheckIntType() const;
	void CheckDoubleType() const;
	void CheckBoolType() const;
	void CheckStringType() const;
	void CheckIntVecType() const;
	void CheckStringVecType() const;

	GuiFieldType_e Type_m;
	string Label_m;
	string InputVal_m;

	int ValueInt_m = -1;
	double ValueDouble_m = -1.;
	string ValueString_m;
	bool ValueBool_m;
	vector<int> ValueIntVec_m;
	vector<string> ValueStringVec_m;
	int FieldID_m[2];
	vector<void*> CallbackFuntions_m;
};

typedef void(*CSMGuiReturnFunc_pf)(const bool GuiSuccess, 
	const vector<GuiField_c> & Fields,
	const vector<GuiField_c> PassthroughFields);

void CSMGui(const string & Title, 
	const vector<GuiField_c> & Fields, 
	CSMGuiReturnFunc_pf ReturnFunc, 
	const AddOn_pa & InputAddOnID,
	const vector<GuiField_c> PassthroughFields = vector<GuiField_c>());

void CSMGuiLabelSelectedPoints(AddOn_pa *AddOnID = NULL);
void CSMGUIDeleteCPLabels(AddOn_pa *AddOnID = NULL);

void CSMGuiLock();
void CSMGuiUnlock();


#endif