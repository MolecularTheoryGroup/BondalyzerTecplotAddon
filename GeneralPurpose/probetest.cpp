
#include "TECADDON.h"
#include "ADDGLBL.h"
#include "GUIDEFS.h"
#if !defined (MSWIN)
#include <unistd.h>
#endif
#include "PROBETEST.h"

#include <sstream>

void ProbeTest_MenuCB()
{
	TecUtilLockOn();
	if (TecUtilDataSetIsAvailable() && TecUtilFrameGetPlotType() != PlotType_Sketch){
		TecUtilProbeInstallCallback(ProbeTest_ProbeCB,
			"Click to run the probe test.");
	}
	else
		TecUtilDialogErrMsg("To execute this add-on Tecplot must have\n"
		"a data set and the frame must be in XY,\n"
		"2D, or 3D.");
	TecUtilLockOff();
}

void STDCALL ProbeTest_ProbeCB(Boolean_t isNearestPoint){
	TecUtilLockStart(AddOnID);

	if (!TecGUIDialogIsUp(Dialog2Manager))
		TecGUIDialogLaunch(Dialog2Manager);

	EntIndex_t VarNum;
	EntIndex_t XYZVarNums[3] = { 0 };

	TecUtilAxisGetVarAssignments(&XYZVarNums[0], &XYZVarNums[1], &XYZVarNums[2]);

	if (XYZVarNums[0] != TECUTILBADVARNUMBER 
		&& XYZVarNums[1] != TECUTILBADVARNUMBER 
		&& XYZVarNums[2] != TECUTILBADVARNUMBER){

		std::stringstream ss;

		if (isNearestPoint)
			ss << "Nearest point probe.\n\n{";
		else
			ss << "Interpolated probe.\n\n{";

		for (VarNum = 0; VarNum < 3; ++VarNum){
			char* VarName;
			if (TecUtilVarGetName(XYZVarNums[VarNum], &VarName)){
				ss << " " << VarName << " ";
				TecUtilStringDealloc(&VarName);
			}
			else
				ss << " ?? ";
		}
		ss << "} = {";

		XYZ_s	PointLoc;
		double* DPtr = &PointLoc.X;

		for (VarNum = 0; VarNum < 3; ++VarNum){
			*DPtr = TecUtilProbeFieldGetValue(XYZVarNums[VarNum]);
			ss << " " << *DPtr << " ";
			DPtr++;
		}

		ss << "}";

		CoordSys_e TheCoordSys = CoordSys_Grid;

		double PaperXY[2] = { TecUtilConvertXPosition(TheCoordSys, CoordSys_Paper, PaperXY[0]),
								TecUtilConvertYPosition(TheCoordSys, CoordSys_Paper, PaperXY[1]) };

		ss << "\nPaper Position: { " << PaperXY[0] << " " << PaperXY[1] << " }";

		TecGUILabelSetText(PointInfoLabel, ss.str().c_str());
	}
	TecUtilLockFinish(AddOnID);
}
