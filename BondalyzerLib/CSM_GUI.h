#pragma once
#ifndef CSMGUI_H_
#define CSMGUI_H_

#ifdef MSWIN
#include <Windows.h>
#endif

using std::vector;
using std::string;

enum GuiFieldType_e{
	Gui_ZoneSelect,
	Gui_ZoneSelectInt,
	Gui_VarSelect,
	Gui_VarSelectInt,

	Gui_VertSep,

	Gui_Toggle,
	Gui_ToggleEnable,

	Gui_Invalid = -1
};

class GuiField_c
{
public:
	GuiField_c();
	GuiField_c(const GuiFieldType_e & Type,
		const string & Label = "",
		const string & SearchString = "");
	~GuiField_c();

	const GuiFieldType_e GetType() const { return Type_m; }
	const string GetLabel() const { return Label_m; }
	const string GetSearchString() const { return SearchString_m; }
	const int GetID() const { return FieldID_m; }

	const int GetReturnInt() const;
	const double GetReturnDouble() const;
	const string GetReturnString() const;
	const bool GetReturnBool() const;

	void SetSearchString(const string & s){ SearchString_m = s; }
	void AppendSearchString(const string & s){ SearchString_m += s; }
	void SetID(const int FieldID){ FieldID_m = FieldID; }

	void SetReturnInt(const int Val);
	void SetReturnDouble(const double Val);
	void SetReturnString(const string Val);
	void SetReturnBool(const bool Val);

private:
	void CheckIntType() const;
	void CheckDoubleType() const;
	void CheckBoolType() const;
	void CheckStringType() const;

	GuiFieldType_e Type_m;
	string Label_m;
	string SearchString_m;

	int ValueInt_m = -1;
	double ValueDouble_m = -1.;
	string ValueString_m;
	bool ValueBool_m;
	int FieldID_m;
};

typedef void(*CSMGuiReturnFunc_pf)(const bool GuiSuccess, const vector<GuiField_c> & Fields);

void CSMGui(const string & Title, const vector<GuiField_c> & Fields, CSMGuiReturnFunc_pf ReturnFunc);



#endif