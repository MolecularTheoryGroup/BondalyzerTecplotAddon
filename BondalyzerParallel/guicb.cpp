/*
*****************************************************************
*****************************************************************
*******                                                  ********
*******    (C) Copyright 1988-2013 by TECPLOT INC.       *******
*******                                                  ********
*****************************************************************
*****************************************************************
*/

#include <vector>
#include <string>
#include <sstream>
#include "TECADDON.h"
#include "ADDGLBL.h"
#include "GUIDEFS.h"
#include "ENGINE.h"

#include "CSM_DATA_SET_INFO.h"
#include "CSM_CALC_VARS.h"

using std::stringstream;
using std::string;
using std::to_string;




/*
 *	Calc vars dialog
 */

/**
*/
static void CalcGrad_TOG_T2_1_CB(const LgIndex_t *I)
{
	TecUtilLockStart(AddOnID);
	TRACE1("Toggle (CalcGrad_TOG_T2_1) Value Changed,  New value is: %d\n", *I);
	TecUtilLockFinish(AddOnID);
}


/**
*/
static void CalcMag_TOG_T2_1_CB(const LgIndex_t *I)
{
	TecUtilLockStart(AddOnID);
	TRACE1("Toggle (CalcMag_TOG_T2_1) Value Changed,  New value is: %d\n", *I);
	TecUtilLockFinish(AddOnID);
}


/**
*/
static void CalcHess_TOG_T2_1_CB(const LgIndex_t *I)
{
	TecUtilLockStart(AddOnID);
	TRACE1("Toggle (CalcHess_TOG_T2_1) Value Changed,  New value is: %d\n", *I);
	TecUtilLockFinish(AddOnID);
}


/**
*/
static void CalcLap_TOG_T2_1_CB(const LgIndex_t *I)
{
	TecUtilLockStart(AddOnID);
	TRACE1("Toggle (CalcLap_TOG_T2_1) Value Changed,  New value is: %d\n", *I);
	TecUtilLockFinish(AddOnID);
}


/**
*/
static void CalcGauss_TOG_T2_1_CB(const LgIndex_t *I)
{
	TecUtilLockStart(AddOnID);
	TRACE1("Toggle (CalcGauss_TOG_T2_1) Value Changed,  New value is: %d\n", *I);
	TecUtilLockFinish(AddOnID);
}


/**
*/
static void CalcEVdotGra_TOG_T2_1_CB(const LgIndex_t *I)
{
	TecUtilLockStart(AddOnID);
	TRACE1("Toggle (CalcEVdotGra_TOG_T2_1) Value Changed,  New value is: %d\n", *I);
	TecUtilLockFinish(AddOnID);
}


/**
*/
static void CalcEb1_TOG_T2_1_CB(const LgIndex_t *I)
{
	TecUtilLockStart(AddOnID);
	TRACE1("Toggle (CalcEb1_TOG_T2_1) Value Changed,  New value is: %d\n", *I);
	TecUtilLockFinish(AddOnID);
}


/**
*/
static void CalcEb2_TOG_T2_1_CB(const LgIndex_t *I)
{
	TecUtilLockStart(AddOnID);
	TRACE1("Toggle (CalcEb2_TOG_T2_1) Value Changed,  New value is: %d\n", *I);
	TecUtilLockFinish(AddOnID);
}


/**
*/
static void CalcEigRank_TOG_T2_1_CB(const LgIndex_t *I)
{
	TecUtilLockStart(AddOnID);
	TRACE1("Toggle (CalcEigRank_TOG_T2_1) Value Changed,  New value is: %d\n", *I);
	TecUtilLockFinish(AddOnID);
}


/**
*/
static void CalcBTN_BTN_T2_1_CB(void)
{
	TecUtilLockStart(AddOnID);
	TRACE("Calculate variables Button Pushed\n");

	vector<LgIndex_t> CalcVarVarListDialogs = {
		RhoVarNum_OPT_T1_1,
		GradXVar_OPT_T1_1,
		GradYVar_OPT_T1_1,
		GradZVar_OPT_T1_1,
		HessXX_OPT_T1_1,
		HessXY_OPT_T1_1,
		HessXZ_OPT_T1_1,
		HessYY_OPT_T1_1,
		HessYZ_OPT_T1_1,
		HessZZ_OPT_T1_1
	};
	vector<LgIndex_t> CalcVarsToggleIDs = {
		CalcGrad_TOG_T2_1,
		CalcMag_TOG_T2_1,
		CalcHess_TOG_T2_1,
		CalcES_TOG_T2_1,
		CalcLap_TOG_T2_1,
		CalcGauss_TOG_T2_1,
		CalcEVdotGra_TOG_T2_1,
		CalcEb1_TOG_T2_1,
		CalcEb2_TOG_T2_1,
		CalcEigRank_TOG_T2_1
	};
	vector<CalcVar_e> CalcVarTypes = {
		CalcGradientVectors,
		CalcGradientMagnitude,
		CalcHessian,
		CalcEigenSystem,
		CalcLaplacian,
		CalcGaussianCurvature,
		CalcEigenVectorsDotGradient,
		CalcEberly1Ridge,
		CalcEberly2Ridge,
		CalcEigenRank
	};

	CalcVarsOptions_s Opt;


	Opt.IsPeriodic = TecGUIToggleGet(PeriodicBC_TOG_T1_1);
	Opt.CalcForAllZones = (TecGUIRadioBoxGetToggle(VolOrZone_RADIO_T1_1) == 2);

	if (!Opt.CalcForAllZones){
		Opt.CalcZoneNum = TecGUIOptionMenuGet(SelZone_OPT_T1_1);
	}

	for (int i = 0; i < CalcVarsToggleIDs.size(); ++i)
		if (TecGUIToggleGet(CalcVarsToggleIDs[i]))
			Opt.CalcVarList.push_back(CalcVarTypes[i]);

	Opt.EberlyUseCutoff[0] = TecGUIToggleGet(EB1Cutoff_TOG_T2_1);
	Opt.EberlyUseCutoff[1] = TecGUIToggleGet(EB2Cutoff_TOG_T2_1);
	TecGUITextFieldGetDouble(EB1CutoffVal_TF_T2_1, &Opt.EberlyCutoff[0]);
	TecGUITextFieldGetDouble(EB2CutoffVal_TF_T2_1, &Opt.EberlyCutoff[1]);

	Opt.AddOnID = AddOnID;

	Opt.HasGrad = TecGUIToggleGet(HasGrad_TOG_T1_1);
	Opt.HasHess = TecGUIToggleGet(HasHess_TOG_T1_1);
	Opt.RhoVarNum = TecGUIOptionMenuGet(RhoVarNum_OPT_T1_1);

	if (Opt.HasGrad)
		for (int i = 0; i < 3; ++i)
			Opt.GradVarNums[i] = TecGUIOptionMenuGet(CalcVarVarListDialogs[1 + i]);
	if (Opt.HasHess)
		for (int i = 0; i < 6; ++i)
			Opt.HessVarNums[i] = TecGUIOptionMenuGet(CalcVarVarListDialogs[4 + i]);

	TecGUIDialogDrop(Dialog2Manager);

	CalcVars(Opt);

	TecUtilLockFinish(AddOnID);
}


/**
*/
static void EB2Cutoff_TOG_T2_1_CB(const LgIndex_t *I)
{
	TecUtilLockStart(AddOnID);
	TRACE1("Toggle (EB2Cutoff_TOG_T2_1) Value Changed,  New value is: %d\n", *I);
	TecUtilLockFinish(AddOnID);
}


/**
*/
static LgIndex_t  EB2CutoffVal_TF_T2_1_CB(const char *S)
{
	LgIndex_t IsOk = 1;
	TecUtilLockStart(AddOnID);
	TRACE1("Text field (EB2CutoffVal_TF_T2_1) Value Changed,  New value is: %s\n", S);
	TecUtilLockFinish(AddOnID);
	return (IsOk);
}


/**
*/
static void EB1Cutoff_TOG_T2_1_CB(const LgIndex_t *I)
{
	TecUtilLockStart(AddOnID);
	TRACE1("Toggle (EB1Cutoff_TOG_T2_1) Value Changed,  New value is: %d\n", *I);
	TecUtilLockFinish(AddOnID);
}


/**
*/
static LgIndex_t  EB1CutoffVal_TF_T2_1_CB(const char *S)
{
	LgIndex_t IsOk = 1;
	TecUtilLockStart(AddOnID);
	TRACE1("Text field (EB1CutoffVal_TF_T2_1) Value Changed,  New value is: %s\n", S);
	TecUtilLockFinish(AddOnID);
	return (IsOk);
}


/**
*/
static void CalcES_TOG_T2_1_CB(const LgIndex_t *I)
{
	TecUtilLockStart(AddOnID);
	TRACE1("Toggle (CalcES_TOG_T2_1) Value Changed,  New value is: %d\n", *I);
	TecUtilLockFinish(AddOnID);
}


char *RhoVarNum_OPT_T1_1_List = "Option 1,Option 2,Option 3";



/**
*/
static void RhoVarNum_OPT_T1_1_CB(const LgIndex_t *I)
{
	TecUtilLockStart(AddOnID);
	TRACE1("Option Menu (RhoVarNum_OPT_T1_1) value changed,  New value is: %d\n", *I);
	TecUtilLockFinish(AddOnID);
}


/**
*/
static void HasGrad_TOG_T1_1_CB(const LgIndex_t *I)
{
	TecUtilLockStart(AddOnID);
	TRACE1("Toggle (HasGrad_TOG_T1_1) Value Changed,  New value is: %d\n", *I);
	TecUtilLockFinish(AddOnID);
}


char *GradZVar_OPT_T1_1_List = "Option 1,Option 2,Option 3";



/**
*/
static void GradZVar_OPT_T1_1_CB(const LgIndex_t *I)
{
	TecUtilLockStart(AddOnID);
	TRACE1("Option Menu (GradZVar_OPT_T1_1) value changed,  New value is: %d\n", *I);
	TecUtilLockFinish(AddOnID);
}


char *GradYVar_OPT_T1_1_List = "Option 1,Option 2,Option 3";



/**
*/
static void GradYVar_OPT_T1_1_CB(const LgIndex_t *I)
{
	TecUtilLockStart(AddOnID);
	TRACE1("Option Menu (GradYVar_OPT_T1_1) value changed,  New value is: %d\n", *I);
	TecUtilLockFinish(AddOnID);
}


char *GradXVar_OPT_T1_1_List = "Option 1,Option 2,Option 3";



/**
*/
static void GradXVar_OPT_T1_1_CB(const LgIndex_t *I)
{
	TecUtilLockStart(AddOnID);
	TRACE1("Option Menu (GradXVar_OPT_T1_1) value changed,  New value is: %d\n", *I);
	TecUtilLockFinish(AddOnID);
}


/**
*/
static void HasHess_TOG_T1_1_CB(const LgIndex_t *I)
{
	TecUtilLockStart(AddOnID);
	TRACE1("Toggle (HasHess_TOG_T1_1) Value Changed,  New value is: %d\n", *I);
	TecUtilLockFinish(AddOnID);
}


char *HessXX_OPT_T1_1_List = "Option 1,Option 2,Option 3";



/**
*/
static void HessXX_OPT_T1_1_CB(const LgIndex_t *I)
{
	TecUtilLockStart(AddOnID);
	TRACE1("Option Menu (HessXX_OPT_T1_1) value changed,  New value is: %d\n", *I);
	TecUtilLockFinish(AddOnID);
}


char *HessZZ_OPT_T1_1_List = "Option 1,Option 2,Option 3";



/**
*/
static void HessZZ_OPT_T1_1_CB(const LgIndex_t *I)
{
	TecUtilLockStart(AddOnID);
	TRACE1("Option Menu (HessZZ_OPT_T1_1) value changed,  New value is: %d\n", *I);
	TecUtilLockFinish(AddOnID);
}


char *HessYZ_OPT_T1_1_List = "Option 1,Option 2,Option 3";



/**
*/
static void HessYZ_OPT_T1_1_CB(const LgIndex_t *I)
{
	TecUtilLockStart(AddOnID);
	TRACE1("Option Menu (HessYZ_OPT_T1_1) value changed,  New value is: %d\n", *I);
	TecUtilLockFinish(AddOnID);
}


char *HessYY_OPT_T1_1_List = "Option 1,Option 2,Option 3";



/**
*/
static void HessYY_OPT_T1_1_CB(const LgIndex_t *I)
{
	TecUtilLockStart(AddOnID);
	TRACE1("Option Menu (HessYY_OPT_T1_1) value changed,  New value is: %d\n", *I);
	TecUtilLockFinish(AddOnID);
}


char *HessXZ_OPT_T1_1_List = "Option 1,Option 2,Option 3";



/**
*/
static void HessXZ_OPT_T1_1_CB(const LgIndex_t *I)
{
	TecUtilLockStart(AddOnID);
	TRACE1("Option Menu (HessXZ_OPT_T1_1) value changed,  New value is: %d\n", *I);
	TecUtilLockFinish(AddOnID);
}


char *HessXY_OPT_T1_1_List = "Option 1,Option 2,Option 3";



/**
*/
static void HessXY_OPT_T1_1_CB(const LgIndex_t *I)
{
	TecUtilLockStart(AddOnID);
	TRACE1("Option Menu (HessXY_OPT_T1_1) value changed,  New value is: %d\n", *I);
	TecUtilLockFinish(AddOnID);
}


/**
*/
static void PeriodicBC_TOG_T1_1_CB(const LgIndex_t *I)
{
	TecUtilLockStart(AddOnID);
	TRACE1("Toggle (PeriodicBC_TOG_T1_1) Value Changed,  New value is: %d\n", *I);
	TecUtilLockFinish(AddOnID);
}


/**
*/
static void VolOrZone_RADIO_T1_1_CB(const LgIndex_t *I)
{
	TecUtilLockStart(AddOnID);
	TRACE1("RadioBox (VolOrZone_RADIO_T1_1) Value Changed,  New value is: %d\n", *I);
	TecUtilLockFinish(AddOnID);
}

char *SelZone_OPT_T1_1_List = "Option 1,Option 2,Option 3";



/**
*/
static void SelZone_OPT_T1_1_CB(const LgIndex_t *I)
{
	TecUtilLockStart(AddOnID);
	TRACE1("Option Menu (SelZone_OPT_T1_1) value changed,  New value is: %d\n", *I);
	TecUtilLockFinish(AddOnID);
}

/**
*/
static void Dialog2HelpButton_CB(void)
{
	TecUtilLockStart(AddOnID);
	TecUtilDialogMessageBox("On-line Help not available for this dialog.",
		MessageBox_Information);
	TecUtilLockFinish(AddOnID);
}


/**
*/
static void Dialog2CloseButton_CB(void)
{
	TecUtilLockStart(AddOnID);
	TecGUIDialogDrop(Dialog2Manager);
	TecUtilLockFinish(AddOnID);
}


/**
*/
static void Dialog2Init_CB(void)
{
	TecUtilLockStart(AddOnID);
	
	/*
	 *	Get dataset info
	 */

	EntIndex_t NumZones, NumVars;
	Boolean_t IsOk = TecUtilDataSetGetInfo(NULL, &NumZones, &NumVars);

	/*
	 *	Populate lists of variables
	 */
	vector<LgIndex_t> CalcVarVarListDialogs = {
		RhoVarNum_OPT_T1_1,
		GradXVar_OPT_T1_1,
		GradYVar_OPT_T1_1,
		GradZVar_OPT_T1_1,
		HessXX_OPT_T1_1,
		HessXY_OPT_T1_1,
		HessXZ_OPT_T1_1,
		HessYY_OPT_T1_1,
		HessYZ_OPT_T1_1,
		HessZZ_OPT_T1_1
	};
	vector<LgIndex_t> CalcVarsToggleIDs = {
		CalcGrad_TOG_T2_1,
		CalcMag_TOG_T2_1,
		CalcHess_TOG_T2_1,
		CalcES_TOG_T2_1,
		CalcLap_TOG_T2_1,
		CalcGauss_TOG_T2_1,
		CalcEVdotGra_TOG_T2_1,
		CalcEb1_TOG_T2_1,
		CalcEb2_TOG_T2_1,
		CalcEigRank_TOG_T2_1
	};
	vector<CalcVar_e> CalcVarTypes = {
		CalcGradientVectors,
		CalcGradientMagnitude,
		CalcHessian,
		CalcEigenSystem,
		CalcLaplacian,
		CalcGaussianCurvature,
		CalcEigenVectorsDotGradient,
		CalcEberly1Ridge,
		CalcEberly2Ridge,
		CalcEigenRank
	};
	for (auto i = CalcVarVarListDialogs.cbegin(), End = CalcVarVarListDialogs.cend(); i != End && IsOk; i++){
		TecGUIOptionMenuDeleteAllItems(*i);
		for (EntIndex_t VarNum = 1; VarNum <= NumVars && IsOk; ++VarNum){
			char *TmpName;
			IsOk = TecUtilVarGetName(VarNum, &TmpName);
			if (IsOk)
				TecGUIOptionMenuAppendItem(*i, TmpName);
			TecUtilStringDealloc(&TmpName);
		}
	}
	
	/*
	 *	Populate list of zones
	 */
	TecGUIOptionMenuDeleteAllItems(SelZone_OPT_T1_1);
	for (int i = 1; i <= NumZones && IsOk; ++i){
		char* ZoneName;
		IsOk = TecUtilZoneGetName(i, &ZoneName);
		if (IsOk){
			TecGUIOptionMenuAppendItem(SelZone_OPT_T1_1, string(to_string(i) + ". " + ZoneName).c_str());
		}
		TecUtilStringDealloc(&ZoneName);
	}
	if (NumZones > 0)
		TecGUIOptionMenuSet(SelZone_OPT_T1_1, 1);

	/*
	 *	Attempt to automatically set variables based
	 *	on the names they 'should' have
	 */
	vector<string> RhoVarNames = {
		"Electron Density",
		"rho",
		"Rho"
	},
	GradVarNames = {
		"X Density",
		"Y Density",
		"Z Density"
	},
	HessVarNames = {
		"XX Density",
		"XY Density",
		"XZ Density",
		"YY Density",
		"YZ Density",
		"ZZ Density"
	};

	int NumHits = 0;

	if (IsOk){
		int TmpNum = VarNumByNameList(RhoVarNames);
		if (TmpNum > 0)
			TecGUIOptionMenuSet(RhoVarNum_OPT_T1_1, TmpNum);

		for (int i = 0; i < 3; ++i){
			TmpNum = VarNumByName(GradVarNames[i]);
			if (TmpNum > 0){
				NumHits++;
				TecGUIOptionMenuSet(CalcVarVarListDialogs[1 + i], TmpNum);
			}
		}
		if (NumHits == 3)
			TecGUIToggleSet(HasGrad_TOG_T1_1, TRUE);
		else
			TecGUIToggleSet(HasGrad_TOG_T1_1, FALSE);

		NumHits = 0;
		for (int i = 0; i < 6; ++i){
			TmpNum = VarNumByName(HessVarNames[i]);
			if (TmpNum > 0){
				NumHits++;
				TecGUIOptionMenuSet(CalcVarVarListDialogs[4 + i], TmpNum);
			}
		}
		if (NumHits == 6)
			TecGUIToggleSet(HasHess_TOG_T1_1, TRUE);
		else
			TecGUIToggleSet(HasHess_TOG_T1_1, FALSE);
	}

	/*
	 *	Set other stuff
	 */
	TecGUIToggleSet(PeriodicBC_TOG_T1_1, FALSE);
	TecGUIRadioBoxSetToggle(VolOrZone_RADIO_T1_1, 1);

	for (auto i = CalcVarsToggleIDs.cbegin(), End = CalcVarsToggleIDs.cend(); i != End; ++i)
		TecGUIToggleSet(*i, FALSE);

	TecGUIToggleSet(EB1Cutoff_TOG_T2_1, FALSE);
	TecGUIToggleSet(EB2Cutoff_TOG_T2_1, FALSE);
	TecGUITextFieldSetString(EB1CutoffVal_TF_T2_1, std::to_string(0.001).c_str());
	TecGUITextFieldSetString(EB2CutoffVal_TF_T2_1, std::to_string(0.001).c_str());

	TecGUITabSetCurrentPage(TAB1_TB_D2, 1);

	TecUtilLockFinish(AddOnID);
}


/**
*/
static void TAB1_TBA_D2_CB(const LgIndex_t *I)
{
	TecUtilLockStart(AddOnID);
	TRACE0("Activate callback for tab (TAB1_TBA_D2) called\n");
	TecUtilLockFinish(AddOnID);
}


/**
*/
static void TAB1_TBD_D2_CB(const LgIndex_t *I)
{
	TecUtilLockStart(AddOnID);
	TRACE0("Deactivate callback for tab (TAB1_TBD_D2) called\n");
	TecUtilLockFinish(AddOnID);
}


#include "guibld.cpp"

