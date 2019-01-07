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
	GuiField_c(GuiFieldType_e const & Type,
		string const & Label = "",
		string const & Val = "",
		const vector<void*> & CallbackFuntions = vector<void*>());
	~GuiField_c();

	GuiFieldType_e GetType() const { return Type_m; }
	string GetLabel() const { return Label_m; }
	string GetSearchString() const { return InputVal_m; }
	int GetID(int IDNum = 0) const { return FieldID_m[IDNum]; }

	int GetReturnInt() const;
	double GetReturnDouble() const;
	string GetReturnString() const;
	bool GetReturnBool() const;
	vector<int> GetReturnIntVec() const;
	vector<string> GetReturnStringVec() const;

	void SetSearchString(string const & s){ InputVal_m = s; }
	void AppendSearchString(string const & s){ InputVal_m += s; }
	void SetID(int const FieldID, int IDNum = 0){ FieldID_m[IDNum] = FieldID; }

	void SetReturnInt(int Val);
	void SetReturnDouble(double const & Val);
	void SetReturnString(string const & Val);
	void SetReturnBool(bool Val);
	void SetReturnIntVec(vector<int> const & Val);
	void SetReturnStringVec(vector<string> const & Val);

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

typedef void(*CSMGuiReturnFunc_pf)(bool const GuiSuccess, 
	vector<GuiField_c> const & Fields,
	vector<GuiField_c> const PassthroughFields);

void CSMGui(string const & Title, 
	vector<GuiField_c> const & Fields, 
	CSMGuiReturnFunc_pf ReturnFunc, 
	AddOn_pa const & InputAddOnID,
	vector<GuiField_c> const PassthroughFields = vector<GuiField_c>());

void CSMGuiLabelSelectedPoints(AddOn_pa *AddOnID = nullptr);
void CSMGUIDeleteCPLabels(AddOn_pa *AddOnID = nullptr);

void CSMGuiLock();
void CSMGuiUnlock();


#endif