/*
*****************************************************************
*****************************************************************
*******                                                  ********
*******    (C) Copyright 1988-2013 by TECPLOT INC.       *******
*******                                                  ********
*****************************************************************
*****************************************************************
*/

/*
	setchembondstyle.cpp
	Tim Wilson
	Copyright 2013 Tecplot Inc. 

	Macro functions to set frame styles for optimal viewing
	of chemical systems used in the ChemBond add-on.
*/

#include <stdio.h>
#include <string.h>
#include <vector>
#include <cmath>
#include "TECADDON.h"
#include "ADDGLBL.h"
#include "ADDONVER.h"
#include "TASSERT.h"
#include "GUIDEFS.h"
#include "ENGINE.h"
#include "ADKUTIL.h"
#include "ALLOC.h"
#include "ARRLIST.h"
#include "SETCHEMBONDSTYLE.h"

#define countOf(staticArray) ( sizeof((staticArray))/sizeof((staticArray)[0]))


void LogDataValues(std::vector<double> &LevelVector, double Min, double Max, double BaseValue);
void ApplyZoneSettings(std::vector<Set_pa> const &ZoneSets, ZoneSettings const &Zone);
void GroupZones(std::vector <Set_pa> &ZoneSets);
void SetupMultiColorContourLevels(std::vector<Set_pa> const &ZoneSets);
void DeallocSets(std::vector<Set_pa> &ZoneSets);



void ChemSysView()
{
	if (TecUtilDataSetIsAvailable())
	{
		std::vector<Set_pa> ZoneSets;
			
		GroupZones(ZoneSets);

		SetupMultiColorContourLevels(ZoneSets);

		ZoneSettings	Zone;

		int				GroupIndex = 0;



		//***************************************************************
		//===============================================================
		//	Adjust properties for each group.
		//===============================================================
		//***************************************************************

		//===============================================================
		//	InitialDataSetZone Settings
		//===============================================================

		Zone.ZoneShow[GroupIndex] =					FALSE;						//	Boolean_t

		//	Mesh Settings

		Zone.MeshShow[GroupIndex] =					FALSE;						//	Boolean_t
		Zone.MeshType[GroupIndex] =					MeshType_Overlay;			//	MeshType_e
		Zone.MeshColor[GroupIndex] =					Black_C;					//	ColorIndex_t
		/*
			A number of color index constants are \#defined. These include Black_C, Red_C,
			Green_C, Blue_C, Cyan_C, Yellow_C, Purple_C, White_C, Custom1_C through
			Custom56_C, MultiColor_C, NoColor_C, MulitiColor2_C through MulitiColor8_C,
			RGBColor_C, and InvalidColor_C.  
		*/
		Zone.MeshLinePttrn[GroupIndex] =				LinePattern_Solid;			//	LinePattern_e
		Zone.MeshLinePttrnLen[GroupIndex] =			2;		//	Value is in %	//	double
		Zone.MeshLineThickness[GroupIndex] =			0.1;	//	Value is in %	//	double


		//	Contour Settings

		Zone.ContourShow[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.ContourType[GroupIndex] =				ContourType_Flood;			//	ContourType_e
		Zone.FloodColoring[GroupIndex] =				ContourColoring_Group1;		//	ContourColoring_e
		Zone.ContourLineGroup[GroupIndex] =			1;							//	SmInteger_t
		Zone.ContourLineColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.ContourLinePttrn[GroupIndex] =			LinePattern_Solid;			//	LinePattern_e
		Zone.ContourLinePttrnLen[GroupIndex] =		2;		//	Value is in %	//	double
		Zone.ContourLineThickness[GroupIndex] =		0.1;	//	Value is in %	//	double
		Zone.ContourLighting[GroupIndex] =			TRUE;						//	Boolean_t

		//	Vector Settings
		
		Zone.VectorShow[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.VectorType[GroupIndex] =				VectorType_TailAtPoint;		//	VectorType_e
		Zone.VectorArrowType[GroupIndex] =			ArrowheadStyle_Plain;		//	ArrowheadStyle_e
		Zone.VectorLineColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.VectorTangent[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.VectorLinePttrn[GroupIndex] =			LinePattern_Solid;			//	LinePattern_e
		Zone.VectorLinePttrnLen[GroupIndex] =		2;		//	Value is in %	//	double
		Zone.VectorLineThickness[GroupIndex] =		0.1;	//	Value is in %	//	double

		//	Scatter Settings

		Zone.ScatterShow[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.ScatterShape[GroupIndex] =				GeomShape_Sphere;			//	GeomShape_e
		Zone.ScatterOutColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.ScatterInColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.ScatterFillMode[GroupIndex] =			FillMode_UseLineColor;		//	FillMode_e
		Zone.ScatterSize[GroupIndex] =				2.5;	//	value is in %	//	double
		Zone.ScatterLineThickness[GroupIndex] =		0.1;	//	value is in %	//	double

		//	Shade Settings

		Zone.ShadeShow[GroupIndex] =					FALSE;						//	Boolean_t
		Zone.ShadeColor[GroupIndex] =				Black_C;					//	ColorIndex_t
		Zone.ShadeLighting[GroupIndex] =				TRUE;						//	Boolean_t

		//	Edge Settings

		Zone.EdgeShow[GroupIndex] =					TRUE;						//	Boolean_t
		Zone.EdgeType[GroupIndex] =					EdgeType_Borders;	//	EdgeType_e
		Zone.EdgeIBorder[GroupIndex] =				BorderLocation_Both;		//	BorderLocation_e
		Zone.EdgeJBorder[GroupIndex] =				BorderLocation_Both;		//	BorderLocation_e
		Zone.EdgeKBorder[GroupIndex] =				BorderLocation_Both;		//	BorderLocation_e
		Zone.EdgeColor[GroupIndex] =					Black_C;					//	ColorIndex_t
		Zone.EdgeLineThickness[GroupIndex] =			0.1;	//	value is in %	//	double

		//	Points Settings

		Zone.PointsToPlot[GroupIndex] =				PointsToPlot_SurfaceNodes;	//	PointsToPlot_e
		Zone.PointsIndexSkipI[GroupIndex] =			1;							//	LgIndex_t
		Zone.PointsIndexSkipJ[GroupIndex] =			1;							//	LgIndex_t
		Zone.PointsIndexSkipK[GroupIndex] =			1;							//	LgIndex_t
		/*
			1 is the same as "No Skip" in the GUI.
		*/

		//	Surfaces Settings

		Zone.SurfacesToPlot[GroupIndex] =			SurfacesToPlot_ExposedCellFaces;	//	SurfacesToPlot_e
		Zone.RangeBeginI[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeEndI[GroupIndex] =					0;							//	LgIndex_t
		Zone.RangeSkipI[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeBeginJ[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeEndJ[GroupIndex] =					0;							//	LgIndex_t
		Zone.RangeSkipJ[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeBeginK[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeEndK[GroupIndex] =					0;							//	LgIndex_t
		Zone.RangeSkipK[GroupIndex] =				1;							//	LgIndex_t
		/*
			0 is the same as maximum, which sets the end range to the maximum
			index for the dimension. The skip number is to skip every # planes,
			so 1 is to not skip any planes, (show every 1 plane).
		*/

		//	Volume Settings

		Zone.StreamtracesShow[GroupIndex] =			TRUE;						//	Boolean_t
		Zone.IsoSurfacesShow[GroupIndex] =			TRUE;						//	Boolean_t
		Zone.SlicesShow[GroupIndex] =				TRUE;						//	Boolean_t

		//	Effects Settings

		Zone.UseValueBlanking[GroupIndex] =			TRUE;						//	Boolean_t
		//Zone.UseClipPlanes[GroupIndex] =			TRUE;						//	Boolean_t
		Zone.UseSurfaceTranslucency[GroupIndex] =	TRUE;						//	Boolean_t
		Zone.SurfaceTranslucency[GroupIndex] =		50;		//	value is in %	//	SmInteger_t
		Zone.LightingEffect[GroupIndex] =			LightingEffect_Gouraud;		//	LightingEffect_e

		//	Group index and number (increments for next set)

		GroupIndex++;														//	int

		//===============================================================
		//	VolumeCPs Settings
		//===============================================================

		Zone.ZoneShow[GroupIndex] =					TRUE;						//	Boolean_t

		//	Mesh Settings

		Zone.MeshShow[GroupIndex] =					FALSE;						//	Boolean_t
		Zone.MeshType[GroupIndex] =					MeshType_Overlay;			//	MeshType_e
		Zone.MeshColor[GroupIndex] =					Black_C;					//	ColorIndex_t
		/*
			A number of color index constants are \#defined. These include Black_C, Red_C,
			Green_C, Blue_C, Cyan_C, Yellow_C, Purple_C, White_C, Custom1_C through
			Custom56_C, MultiColor_C, NoColor_C, MulitiColor2_C through MulitiColor8_C,
			RGBColor_C, and InvalidColor_C.  
		*/
		Zone.MeshLinePttrn[GroupIndex] =				LinePattern_Solid;			//	LinePattern_e
		Zone.MeshLinePttrnLen[GroupIndex] =			2;		//	Value is in %	//	double
		Zone.MeshLineThickness[GroupIndex] =			0.1;	//	Value is in %	//	double


		//	Contour Settings

		Zone.ContourShow[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.ContourType[GroupIndex] =				ContourType_Flood;			//	ContourType_e
		Zone.FloodColoring[GroupIndex] =				ContourColoring_Group1;		//	ContourColoring_e
		Zone.ContourLineGroup[GroupIndex] =			1;							//	SmInteger_t
		Zone.ContourLineColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.ContourLinePttrn[GroupIndex] =			LinePattern_Solid;			//	LinePattern_e
		Zone.ContourLinePttrnLen[GroupIndex] =		2;		//	Value is in %	//	double
		Zone.ContourLineThickness[GroupIndex] =		0.1;	//	Value is in %	//	double
		Zone.ContourLighting[GroupIndex] =			TRUE;						//	Boolean_t

		//	Vector Settings
		
		Zone.VectorShow[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.VectorType[GroupIndex] =				VectorType_TailAtPoint;		//	VectorType_e
		Zone.VectorArrowType[GroupIndex] =			ArrowheadStyle_Plain;		//	ArrowheadStyle_e
		Zone.VectorLineColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.VectorTangent[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.VectorLinePttrn[GroupIndex] =			LinePattern_Solid;			//	LinePattern_e
		Zone.VectorLinePttrnLen[GroupIndex] =		2;		//	Value is in %	//	double
		Zone.VectorLineThickness[GroupIndex] =		0.1;	//	Value is in %	//	double

		//	Scatter Settings

		Zone.ScatterShow[GroupIndex] =				TRUE;						//	Boolean_t
		Zone.ScatterShape[GroupIndex] =				GeomShape_Sphere;			//	GeomShape_e
		Zone.ScatterOutColor[GroupIndex] =			MultiColor_C;				//	ColorIndex_t
		Zone.ScatterInColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.ScatterFillMode[GroupIndex] =			FillMode_UseLineColor;		//	FillMode_e
		Zone.ScatterSize[GroupIndex] =				2.5;	//	value is in %	//	double
		Zone.ScatterLineThickness[GroupIndex] =		0.1;	//	value is in %	//	double

		//	Shade Settings

		Zone.ShadeShow[GroupIndex] =					FALSE;						//	Boolean_t
		Zone.ShadeColor[GroupIndex] =				Black_C;					//	ColorIndex_t
		Zone.ShadeLighting[GroupIndex] =				TRUE;						//	Boolean_t

		//	Edge Settings

		Zone.EdgeShow[GroupIndex] =					FALSE;						//	Boolean_t
		Zone.EdgeType[GroupIndex] =					EdgeType_BordersAndCreases;	//	EdgeType_e
		Zone.EdgeIBorder[GroupIndex] =				BorderLocation_Both;		//	BorderLocation_e
		Zone.EdgeJBorder[GroupIndex] =				BorderLocation_Both;		//	BorderLocation_e
		Zone.EdgeKBorder[GroupIndex] =				BorderLocation_Both;		//	BorderLocation_e
		Zone.EdgeColor[GroupIndex] =					Black_C;					//	ColorIndex_t
		Zone.EdgeLineThickness[GroupIndex] =			0.1;	//	value is in %	//	double

		//	Points Settings

		Zone.PointsToPlot[GroupIndex] =				PointsToPlot_SurfaceNodes;	//	PointsToPlot_e
		Zone.PointsIndexSkipI[GroupIndex] =			1;							//	LgIndex_t
		Zone.PointsIndexSkipJ[GroupIndex] =			1;							//	LgIndex_t
		Zone.PointsIndexSkipK[GroupIndex] =			1;							//	LgIndex_t
		/*
			1 is the same as "No Skip" in the GUI.
		*/

		//	Surfaces Settings

		Zone.SurfacesToPlot[GroupIndex] =			SurfacesToPlot_ExposedCellFaces;	//	SurfacesToPlot_e
		Zone.RangeBeginI[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeEndI[GroupIndex] =					0;							//	LgIndex_t
		Zone.RangeSkipI[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeBeginJ[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeEndJ[GroupIndex] =					0;							//	LgIndex_t
		Zone.RangeSkipJ[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeBeginK[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeEndK[GroupIndex] =					0;							//	LgIndex_t
		Zone.RangeSkipK[GroupIndex] =				1;							//	LgIndex_t
		/*
			0 is the same as maximum, which sets the end range to the maximum
			index for the dimension. The skip number is to skip every # planes,
			so 1 is to not skip any planes, (show every 1 plane).
		*/

		//	Volume Settings

		Zone.StreamtracesShow[GroupIndex] =			TRUE;						//	Boolean_t
		Zone.IsoSurfacesShow[GroupIndex] =			TRUE;						//	Boolean_t
		Zone.SlicesShow[GroupIndex] =				TRUE;						//	Boolean_t

		//	Effects Settings

		Zone.UseValueBlanking[GroupIndex] =			FALSE;						//	Boolean_t
		//Zone.UseClipPlanes[GroupIndex] =			TRUE;						//	Boolean_t
		Zone.UseSurfaceTranslucency[GroupIndex] =	FALSE;						//	Boolean_t
		Zone.SurfaceTranslucency[GroupIndex] =		50;		//	value is in %	//	SmInteger_t
		Zone.LightingEffect[GroupIndex] =			LightingEffect_Gouraud;		//	LightingEffect_e

		//	Group index and number (increments for next set)

		GroupIndex++;														//	int


		//===============================================================
		//	IsoSurfCPs Settings
		//===============================================================

		Zone.ZoneShow[GroupIndex] =					FALSE;						//	Boolean_t

		//	Mesh Settings

		Zone.MeshShow[GroupIndex] =					FALSE;						//	Boolean_t
		Zone.MeshType[GroupIndex] =					MeshType_Overlay;			//	MeshType_e
		Zone.MeshColor[GroupIndex] =					Black_C;					//	ColorIndex_t
		/*
			A number of color index constants are \#defined. These include Black_C, Red_C,
			Green_C, Blue_C, Cyan_C, Yellow_C, Purple_C, White_C, Custom1_C through
			Custom56_C, MultiColor_C, NoColor_C, MulitiColor2_C through MulitiColor8_C,
			RGBColor_C, and InvalidColor_C.  
		*/
		Zone.MeshLinePttrn[GroupIndex] =				LinePattern_Solid;			//	LinePattern_e
		Zone.MeshLinePttrnLen[GroupIndex] =			2;		//	Value is in %	//	double
		Zone.MeshLineThickness[GroupIndex] =			0.1;	//	Value is in %	//	double


		//	Contour Settings

		Zone.ContourShow[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.ContourType[GroupIndex] =				ContourType_Flood;			//	ContourType_e
		Zone.FloodColoring[GroupIndex] =				ContourColoring_Group3;		//	ContourColoring_e
		Zone.ContourLineGroup[GroupIndex] =			1;							//	SmInteger_t
		Zone.ContourLineColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.ContourLinePttrn[GroupIndex] =			LinePattern_Solid;			//	LinePattern_e
		Zone.ContourLinePttrnLen[GroupIndex] =		2;		//	Value is in %	//	double
		Zone.ContourLineThickness[GroupIndex] =		0.1;	//	Value is in %	//	double
		Zone.ContourLighting[GroupIndex] =			TRUE;						//	Boolean_t

		//	Vector Settings
		
		Zone.VectorShow[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.VectorType[GroupIndex] =				VectorType_TailAtPoint;		//	VectorType_e
		Zone.VectorArrowType[GroupIndex] =			ArrowheadStyle_Plain;		//	ArrowheadStyle_e
		Zone.VectorLineColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.VectorTangent[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.VectorLinePttrn[GroupIndex] =			LinePattern_Solid;			//	LinePattern_e
		Zone.VectorLinePttrnLen[GroupIndex] =		2;		//	Value is in %	//	double
		Zone.VectorLineThickness[GroupIndex] =		0.1;	//	Value is in %	//	double

		//	Scatter Settings

		Zone.ScatterShow[GroupIndex] =				TRUE;						//	Boolean_t
		Zone.ScatterShape[GroupIndex] =				GeomShape_Octahedron;			//	GeomShape_e
		Zone.ScatterOutColor[GroupIndex] =			MultiColor2_C;					//	ColorIndex_t
		Zone.ScatterInColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.ScatterFillMode[GroupIndex] =			FillMode_UseLineColor;		//	FillMode_e
		Zone.ScatterSize[GroupIndex] =				2.5;	//	value is in %	//	double
		Zone.ScatterLineThickness[GroupIndex] =		0.1;	//	value is in %	//	double

		//	Shade Settings

		Zone.ShadeShow[GroupIndex] =					FALSE;						//	Boolean_t
		Zone.ShadeColor[GroupIndex] =				Black_C;					//	ColorIndex_t
		Zone.ShadeLighting[GroupIndex] =				TRUE;						//	Boolean_t

		//	Edge Settings

		Zone.EdgeShow[GroupIndex] =					FALSE;						//	Boolean_t
		Zone.EdgeType[GroupIndex] =					EdgeType_BordersAndCreases;	//	EdgeType_e
		Zone.EdgeIBorder[GroupIndex] =				BorderLocation_Both;		//	BorderLocation_e
		Zone.EdgeJBorder[GroupIndex] =				BorderLocation_Both;		//	BorderLocation_e
		Zone.EdgeKBorder[GroupIndex] =				BorderLocation_Both;		//	BorderLocation_e
		Zone.EdgeColor[GroupIndex] =					Black_C;					//	ColorIndex_t
		Zone.EdgeLineThickness[GroupIndex] =			0.1;	//	value is in %	//	double

		//	Points Settings

		Zone.PointsToPlot[GroupIndex] =				PointsToPlot_SurfaceNodes;	//	PointsToPlot_e
		Zone.PointsIndexSkipI[GroupIndex] =			1;							//	LgIndex_t
		Zone.PointsIndexSkipJ[GroupIndex] =			1;							//	LgIndex_t
		Zone.PointsIndexSkipK[GroupIndex] =			1;							//	LgIndex_t
		/*
			1 is the same as "No Skip" in the GUI.
		*/

		//	Surfaces Settings

		Zone.SurfacesToPlot[GroupIndex] =			SurfacesToPlot_ExposedCellFaces;	//	SurfacesToPlot_e
		Zone.RangeBeginI[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeEndI[GroupIndex] =					0;							//	LgIndex_t
		Zone.RangeSkipI[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeBeginJ[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeEndJ[GroupIndex] =					0;							//	LgIndex_t
		Zone.RangeSkipJ[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeBeginK[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeEndK[GroupIndex] =					0;							//	LgIndex_t
		Zone.RangeSkipK[GroupIndex] =				1;							//	LgIndex_t
		/*
			0 is the same as maximum, which sets the end range to the maximum
			index for the dimension. The skip number is to skip every # planes,
			so 1 is to not skip any planes, (show every 1 plane).
		*/

		//	Volume Settings

		Zone.StreamtracesShow[GroupIndex] =			TRUE;						//	Boolean_t
		Zone.IsoSurfacesShow[GroupIndex] =			TRUE;						//	Boolean_t
		Zone.SlicesShow[GroupIndex] =				TRUE;						//	Boolean_t

		//	Effects Settings

		Zone.UseValueBlanking[GroupIndex] =			TRUE;						//	Boolean_t
		//Zone.UseClipPlanes[GroupIndex] =			TRUE;						//	Boolean_t
		Zone.UseSurfaceTranslucency[GroupIndex] =	FALSE;						//	Boolean_t
		Zone.SurfaceTranslucency[GroupIndex] =		50;		//	value is in %	//	SmInteger_t
		Zone.LightingEffect[GroupIndex] =			LightingEffect_Gouraud;		//	LightingEffect_e

		//	Group index and number (increments for next set)

		GroupIndex++;														//	int


		//===============================================================
		//	ABLines Settings
		//===============================================================

		Zone.ZoneShow[GroupIndex] =					TRUE;						//	Boolean_t

		//	Mesh Settings

		Zone.MeshShow[GroupIndex] =					TRUE;						//	Boolean_t
		Zone.MeshType[GroupIndex] =					MeshType_Overlay;			//	MeshType_e
		Zone.MeshColor[GroupIndex] =					MultiColor4_C;					//	ColorIndex_t
		/*
			A number of color index constants are \#defined. These include Black_C, Red_C,
			Green_C, Blue_C, Cyan_C, Yellow_C, Purple_C, White_C, Custom1_C through
			Custom56_C, MultiColor_C, NoColor_C, MulitiColor2_C through MulitiColor8_C,
			RGBColor_C, and InvalidColor_C.  
		*/
		Zone.MeshLinePttrn[GroupIndex] =				LinePattern_Solid;			//	LinePattern_e
		Zone.MeshLinePttrnLen[GroupIndex] =			2;		//	Value is in %	//	double
		Zone.MeshLineThickness[GroupIndex] =			0.3;		//	Value is in %	//	double


		//	Contour Settings

		Zone.ContourShow[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.ContourType[GroupIndex] =				ContourType_Flood;			//	ContourType_e
		Zone.FloodColoring[GroupIndex] =				ContourColoring_Group1;		//	ContourColoring_e
		Zone.ContourLineGroup[GroupIndex] =			1;							//	SmInteger_t
		Zone.ContourLineColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.ContourLinePttrn[GroupIndex] =			LinePattern_Solid;			//	LinePattern_e
		Zone.ContourLinePttrnLen[GroupIndex] =		2;		//	Value is in %	//	double
		Zone.ContourLineThickness[GroupIndex] =		0.1;	//	Value is in %	//	double
		Zone.ContourLighting[GroupIndex] =			TRUE;						//	Boolean_t

		//	Vector Settings
		
		Zone.VectorShow[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.VectorType[GroupIndex] =				VectorType_TailAtPoint;		//	VectorType_e
		Zone.VectorArrowType[GroupIndex] =			ArrowheadStyle_Plain;		//	ArrowheadStyle_e
		Zone.VectorLineColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.VectorTangent[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.VectorLinePttrn[GroupIndex] =			LinePattern_Solid;			//	LinePattern_e
		Zone.VectorLinePttrnLen[GroupIndex] =		2;		//	Value is in %	//	double
		Zone.VectorLineThickness[GroupIndex] =		0.1;	//	Value is in %	//	double

		//	Scatter Settings

		Zone.ScatterShow[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.ScatterShape[GroupIndex] =				GeomShape_Sphere;			//	GeomShape_e
		Zone.ScatterOutColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.ScatterInColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.ScatterFillMode[GroupIndex] =			FillMode_UseLineColor;		//	FillMode_e
		Zone.ScatterSize[GroupIndex] =				2.5;	//	value is in %	//	double
		Zone.ScatterLineThickness[GroupIndex] =		0.1;	//	value is in %	//	double

		//	Shade Settings

		Zone.ShadeShow[GroupIndex] =					FALSE;						//	Boolean_t
		Zone.ShadeColor[GroupIndex] =				Black_C;					//	ColorIndex_t
		Zone.ShadeLighting[GroupIndex] =				TRUE;						//	Boolean_t

		//	Edge Settings

		Zone.EdgeShow[GroupIndex] =					FALSE;						//	Boolean_t
		Zone.EdgeType[GroupIndex] =					EdgeType_BordersAndCreases;	//	EdgeType_e
		Zone.EdgeIBorder[GroupIndex] =				BorderLocation_Both;		//	BorderLocation_e
		Zone.EdgeJBorder[GroupIndex] =				BorderLocation_Both;		//	BorderLocation_e
		Zone.EdgeKBorder[GroupIndex] =				BorderLocation_Both;		//	BorderLocation_e
		Zone.EdgeColor[GroupIndex] =					Black_C;					//	ColorIndex_t
		Zone.EdgeLineThickness[GroupIndex] =			0.1;	//	value is in %	//	double

		//	Points Settings

		Zone.PointsToPlot[GroupIndex] =				PointsToPlot_SurfaceNodes;	//	PointsToPlot_e
		Zone.PointsIndexSkipI[GroupIndex] =			1;							//	LgIndex_t
		Zone.PointsIndexSkipJ[GroupIndex] =			1;							//	LgIndex_t
		Zone.PointsIndexSkipK[GroupIndex] =			1;							//	LgIndex_t
		/*
			1 is the same as "No Skip" in the GUI.
		*/

		//	Surfaces Settings

		Zone.SurfacesToPlot[GroupIndex] =			SurfacesToPlot_ExposedCellFaces;	//	SurfacesToPlot_e
		Zone.RangeBeginI[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeEndI[GroupIndex] =					0;							//	LgIndex_t
		Zone.RangeSkipI[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeBeginJ[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeEndJ[GroupIndex] =					0;							//	LgIndex_t
		Zone.RangeSkipJ[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeBeginK[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeEndK[GroupIndex] =					0;							//	LgIndex_t
		Zone.RangeSkipK[GroupIndex] =				1;							//	LgIndex_t
		/*
			0 is the same as maximum, which sets the end range to the maximum
			index for the dimension. The skip number is to skip every # planes,
			so 1 is to not skip any planes, (show every 1 plane).
		*/

		//	Volume Settings

		Zone.StreamtracesShow[GroupIndex] =			TRUE;						//	Boolean_t
		Zone.IsoSurfacesShow[GroupIndex] =			TRUE;						//	Boolean_t
		Zone.SlicesShow[GroupIndex] =				TRUE;						//	Boolean_t

		//	Effects Settings

		Zone.UseValueBlanking[GroupIndex] =			FALSE;						//	Boolean_t
		//Zone.UseClipPlanes[GroupIndex] =			TRUE;						//	Boolean_t
		Zone.UseSurfaceTranslucency[GroupIndex] =	FALSE;						//	Boolean_t
		Zone.SurfaceTranslucency[GroupIndex] =		50;		//	value is in %	//	SmInteger_t
		Zone.LightingEffect[GroupIndex] =			LightingEffect_Gouraud;		//	LightingEffect_e

		//	Group index and number (increments for next set)

		GroupIndex++;														//	int


		//===============================================================
		//	ABLinesFF Settings
		//===============================================================

		Zone.ZoneShow[GroupIndex] =					FALSE;						//	Boolean_t

		//	Mesh Settings

		Zone.MeshShow[GroupIndex] =					TRUE;						//	Boolean_t
		Zone.MeshType[GroupIndex] =					MeshType_Overlay;			//	MeshType_e
		Zone.MeshColor[GroupIndex] =					MultiColor4_C;					//	ColorIndex_t
		/*
			A number of color index constants are \#defined. These include Black_C, Red_C,
			Green_C, Blue_C, Cyan_C, Yellow_C, Purple_C, White_C, Custom1_C through
			Custom56_C, MultiColor_C, NoColor_C, MulitiColor2_C through MulitiColor8_C,
			RGBColor_C, and InvalidColor_C.  
		*/
		Zone.MeshLinePttrn[GroupIndex] =				LinePattern_Solid;			//	LinePattern_e
		Zone.MeshLinePttrnLen[GroupIndex] =			2;		//	Value is in %	//	double
		Zone.MeshLineThickness[GroupIndex] =			0.2;	//	Value is in %	//	double


		//	Contour Settings

		Zone.ContourShow[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.ContourType[GroupIndex] =				ContourType_Flood;			//	ContourType_e
		Zone.FloodColoring[GroupIndex] =				ContourColoring_Group1;		//	ContourColoring_e
		Zone.ContourLineGroup[GroupIndex] =			1;							//	SmInteger_t
		Zone.ContourLineColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.ContourLinePttrn[GroupIndex] =			LinePattern_Solid;			//	LinePattern_e
		Zone.ContourLinePttrnLen[GroupIndex] =		2;		//	Value is in %	//	double
		Zone.ContourLineThickness[GroupIndex] =		0.1;	//	Value is in %	//	double
		Zone.ContourLighting[GroupIndex] =			TRUE;						//	Boolean_t

		//	Vector Settings
		
		Zone.VectorShow[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.VectorType[GroupIndex] =				VectorType_TailAtPoint;		//	VectorType_e
		Zone.VectorArrowType[GroupIndex] =			ArrowheadStyle_Plain;		//	ArrowheadStyle_e
		Zone.VectorLineColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.VectorTangent[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.VectorLinePttrn[GroupIndex] =			LinePattern_Solid;			//	LinePattern_e
		Zone.VectorLinePttrnLen[GroupIndex] =		2;		//	Value is in %	//	double
		Zone.VectorLineThickness[GroupIndex] =		0.1;	//	Value is in %	//	double

		//	Scatter Settings

		Zone.ScatterShow[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.ScatterShape[GroupIndex] =				GeomShape_Sphere;			//	GeomShape_e
		Zone.ScatterOutColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.ScatterInColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.ScatterFillMode[GroupIndex] =			FillMode_UseLineColor;		//	FillMode_e
		Zone.ScatterSize[GroupIndex] =				2.5;	//	value is in %	//	double
		Zone.ScatterLineThickness[GroupIndex] =		0.1;	//	value is in %	//	double

		//	Shade Settings

		Zone.ShadeShow[GroupIndex] =					FALSE;						//	Boolean_t
		Zone.ShadeColor[GroupIndex] =				Black_C;					//	ColorIndex_t
		Zone.ShadeLighting[GroupIndex] =				TRUE;						//	Boolean_t

		//	Edge Settings

		Zone.EdgeShow[GroupIndex] =					FALSE;						//	Boolean_t
		Zone.EdgeType[GroupIndex] =					EdgeType_BordersAndCreases;	//	EdgeType_e
		Zone.EdgeIBorder[GroupIndex] =				BorderLocation_Both;		//	BorderLocation_e
		Zone.EdgeJBorder[GroupIndex] =				BorderLocation_Both;		//	BorderLocation_e
		Zone.EdgeKBorder[GroupIndex] =				BorderLocation_Both;		//	BorderLocation_e
		Zone.EdgeColor[GroupIndex] =					Black_C;					//	ColorIndex_t
		Zone.EdgeLineThickness[GroupIndex] =			0.1;	//	value is in %	//	double

		//	Points Settings

		Zone.PointsToPlot[GroupIndex] =				PointsToPlot_SurfaceNodes;	//	PointsToPlot_e
		Zone.PointsIndexSkipI[GroupIndex] =			1;							//	LgIndex_t
		Zone.PointsIndexSkipJ[GroupIndex] =			1;							//	LgIndex_t
		Zone.PointsIndexSkipK[GroupIndex] =			1;							//	LgIndex_t
		/*
			1 is the same as "No Skip" in the GUI.
		*/

		//	Surfaces Settings

		Zone.SurfacesToPlot[GroupIndex] =			SurfacesToPlot_ExposedCellFaces;	//	SurfacesToPlot_e
		Zone.RangeBeginI[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeEndI[GroupIndex] =					0;							//	LgIndex_t
		Zone.RangeSkipI[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeBeginJ[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeEndJ[GroupIndex] =					0;							//	LgIndex_t
		Zone.RangeSkipJ[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeBeginK[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeEndK[GroupIndex] =					0;							//	LgIndex_t
		Zone.RangeSkipK[GroupIndex] =				1;							//	LgIndex_t
		/*
			0 is the same as maximum, which sets the end range to the maximum
			index for the dimension. The skip number is to skip every # planes,
			so 1 is to not skip any planes, (show every 1 plane).
		*/

		//	Volume Settings

		Zone.StreamtracesShow[GroupIndex] =			TRUE;						//	Boolean_t
		Zone.IsoSurfacesShow[GroupIndex] =			TRUE;						//	Boolean_t
		Zone.SlicesShow[GroupIndex] =				TRUE;						//	Boolean_t

		//	Effects Settings

		Zone.UseValueBlanking[GroupIndex] =			TRUE;						//	Boolean_t
		//Zone.UseClipPlanes[GroupIndex] =			TRUE;						//	Boolean_t
		Zone.UseSurfaceTranslucency[GroupIndex] =	FALSE;						//	Boolean_t
		Zone.SurfaceTranslucency[GroupIndex] =		50;		//	value is in %	//	SmInteger_t
		Zone.LightingEffect[GroupIndex] =			LightingEffect_Gouraud;		//	LightingEffect_e

		//	Group index and number (increments for next set)

		GroupIndex++;														//	int


		//===============================================================
		//	ARLines Settings
		//===============================================================

		Zone.ZoneShow[GroupIndex] =					FALSE;						//	Boolean_t

		//	Mesh Settings

		Zone.MeshShow[GroupIndex] =					TRUE;						//	Boolean_t
		Zone.MeshType[GroupIndex] =					MeshType_Overlay;			//	MeshType_e
		Zone.MeshColor[GroupIndex] =					MultiColor5_C;					//	ColorIndex_t
		/*
			A number of color index constants are \#defined. These include Black_C, Red_C,
			Green_C, Blue_C, Cyan_C, Yellow_C, Purple_C, White_C, Custom1_C through
			Custom56_C, MultiColor_C, NoColor_C, MulitiColor2_C through MulitiColor8_C,
			RGBColor_C, and InvalidColor_C.  
		*/
		Zone.MeshLinePttrn[GroupIndex] =				LinePattern_LongDash;			//	LinePattern_e
		Zone.MeshLinePttrnLen[GroupIndex] =			2;		//	Value is in %	//	double
		Zone.MeshLineThickness[GroupIndex] =			0.2;	//	Value is in %	//	double


		//	Contour Settings

		Zone.ContourShow[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.ContourType[GroupIndex] =				ContourType_Flood;			//	ContourType_e
		Zone.FloodColoring[GroupIndex] =				ContourColoring_Group1;		//	ContourColoring_e
		Zone.ContourLineGroup[GroupIndex] =			1;							//	SmInteger_t
		Zone.ContourLineColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.ContourLinePttrn[GroupIndex] =			LinePattern_Solid;			//	LinePattern_e
		Zone.ContourLinePttrnLen[GroupIndex] =		2;		//	Value is in %	//	double
		Zone.ContourLineThickness[GroupIndex] =		0.1;	//	Value is in %	//	double
		Zone.ContourLighting[GroupIndex] =			TRUE;						//	Boolean_t

		//	Vector Settings
		
		Zone.VectorShow[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.VectorType[GroupIndex] =				VectorType_TailAtPoint;		//	VectorType_e
		Zone.VectorArrowType[GroupIndex] =			ArrowheadStyle_Plain;		//	ArrowheadStyle_e
		Zone.VectorLineColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.VectorTangent[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.VectorLinePttrn[GroupIndex] =			LinePattern_Solid;			//	LinePattern_e
		Zone.VectorLinePttrnLen[GroupIndex] =		2;		//	Value is in %	//	double
		Zone.VectorLineThickness[GroupIndex] =		0.1;	//	Value is in %	//	double

		//	Scatter Settings

		Zone.ScatterShow[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.ScatterShape[GroupIndex] =				GeomShape_Sphere;			//	GeomShape_e
		Zone.ScatterOutColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.ScatterInColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.ScatterFillMode[GroupIndex] =			FillMode_UseLineColor;		//	FillMode_e
		Zone.ScatterSize[GroupIndex] =				2.5;	//	value is in %	//	double
		Zone.ScatterLineThickness[GroupIndex] =		0.1;	//	value is in %	//	double

		//	Shade Settings

		Zone.ShadeShow[GroupIndex] =					FALSE;						//	Boolean_t
		Zone.ShadeColor[GroupIndex] =				Black_C;					//	ColorIndex_t
		Zone.ShadeLighting[GroupIndex] =				TRUE;						//	Boolean_t

		//	Edge Settings

		Zone.EdgeShow[GroupIndex] =					FALSE;						//	Boolean_t
		Zone.EdgeType[GroupIndex] =					EdgeType_BordersAndCreases;	//	EdgeType_e
		Zone.EdgeIBorder[GroupIndex] =				BorderLocation_Both;		//	BorderLocation_e
		Zone.EdgeJBorder[GroupIndex] =				BorderLocation_Both;		//	BorderLocation_e
		Zone.EdgeKBorder[GroupIndex] =				BorderLocation_Both;		//	BorderLocation_e
		Zone.EdgeColor[GroupIndex] =					Black_C;					//	ColorIndex_t
		Zone.EdgeLineThickness[GroupIndex] =			0.1;	//	value is in %	//	double

		//	Points Settings

		Zone.PointsToPlot[GroupIndex] =				PointsToPlot_SurfaceNodes;	//	PointsToPlot_e
		Zone.PointsIndexSkipI[GroupIndex] =			1;							//	LgIndex_t
		Zone.PointsIndexSkipJ[GroupIndex] =			1;							//	LgIndex_t
		Zone.PointsIndexSkipK[GroupIndex] =			1;							//	LgIndex_t
		/*
			1 is the same as "No Skip" in the GUI.
		*/

		//	Surfaces Settings

		Zone.SurfacesToPlot[GroupIndex] =			SurfacesToPlot_ExposedCellFaces;	//	SurfacesToPlot_e
		Zone.RangeBeginI[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeEndI[GroupIndex] =					0;							//	LgIndex_t
		Zone.RangeSkipI[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeBeginJ[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeEndJ[GroupIndex] =					0;							//	LgIndex_t
		Zone.RangeSkipJ[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeBeginK[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeEndK[GroupIndex] =					0;							//	LgIndex_t
		Zone.RangeSkipK[GroupIndex] =				1;							//	LgIndex_t
		/*
			0 is the same as maximum, which sets the end range to the maximum
			index for the dimension. The skip number is to skip every # planes,
			so 1 is to not skip any planes, (show every 1 plane).
		*/

		//	Volume Settings

		Zone.StreamtracesShow[GroupIndex] =			TRUE;						//	Boolean_t
		Zone.IsoSurfacesShow[GroupIndex] =			TRUE;						//	Boolean_t
		Zone.SlicesShow[GroupIndex] =				TRUE;						//	Boolean_t

		//	Effects Settings

		Zone.UseValueBlanking[GroupIndex] =			TRUE;						//	Boolean_t
		//Zone.UseClipPlanes[GroupIndex] =			TRUE;						//	Boolean_t
		Zone.UseSurfaceTranslucency[GroupIndex] =	FALSE;						//	Boolean_t
		Zone.SurfaceTranslucency[GroupIndex] =		50;		//	value is in %	//	SmInteger_t
		Zone.LightingEffect[GroupIndex] =			LightingEffect_Gouraud;		//	LightingEffect_e

		//	Group index and number (increments for next set)

		GroupIndex++;														//	int


		//===============================================================
		//	ARLinesFF Settings
		//===============================================================

		Zone.ZoneShow[GroupIndex] =					FALSE;						//	Boolean_t

		//	Mesh Settings

		Zone.MeshShow[GroupIndex] =					TRUE;						//	Boolean_t
		Zone.MeshType[GroupIndex] =					MeshType_Overlay;			//	MeshType_e
		Zone.MeshColor[GroupIndex] =					MultiColor5_C;					//	ColorIndex_t
		/*
			A number of color index constants are \#defined. These include Black_C, Red_C,
			Green_C, Blue_C, Cyan_C, Yellow_C, Purple_C, White_C, Custom1_C through
			Custom56_C, MultiColor_C, NoColor_C, MulitiColor2_C through MulitiColor8_C,
			RGBColor_C, and InvalidColor_C.  
		*/
		Zone.MeshLinePttrn[GroupIndex] =				LinePattern_LongDash;			//	LinePattern_e
		Zone.MeshLinePttrnLen[GroupIndex] =			2;		//	Value is in %	//	double
		Zone.MeshLineThickness[GroupIndex] =			0.2;	//	Value is in %	//	double


		//	Contour Settings

		Zone.ContourShow[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.ContourType[GroupIndex] =				ContourType_Flood;			//	ContourType_e
		Zone.FloodColoring[GroupIndex] =				ContourColoring_Group1;		//	ContourColoring_e
		Zone.ContourLineGroup[GroupIndex] =			1;							//	SmInteger_t
		Zone.ContourLineColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.ContourLinePttrn[GroupIndex] =			LinePattern_Solid;			//	LinePattern_e
		Zone.ContourLinePttrnLen[GroupIndex] =		2;		//	Value is in %	//	double
		Zone.ContourLineThickness[GroupIndex] =		0.1;	//	Value is in %	//	double
		Zone.ContourLighting[GroupIndex] =			TRUE;						//	Boolean_t

		//	Vector Settings
		
		Zone.VectorShow[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.VectorType[GroupIndex] =				VectorType_TailAtPoint;		//	VectorType_e
		Zone.VectorArrowType[GroupIndex] =			ArrowheadStyle_Plain;		//	ArrowheadStyle_e
		Zone.VectorLineColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.VectorTangent[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.VectorLinePttrn[GroupIndex] =			LinePattern_Solid;			//	LinePattern_e
		Zone.VectorLinePttrnLen[GroupIndex] =		2;		//	Value is in %	//	double
		Zone.VectorLineThickness[GroupIndex] =		0.1;	//	Value is in %	//	double

		//	Scatter Settings

		Zone.ScatterShow[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.ScatterShape[GroupIndex] =				GeomShape_Sphere;			//	GeomShape_e
		Zone.ScatterOutColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.ScatterInColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.ScatterFillMode[GroupIndex] =			FillMode_UseLineColor;		//	FillMode_e
		Zone.ScatterSize[GroupIndex] =				2.5;	//	value is in %	//	double
		Zone.ScatterLineThickness[GroupIndex] =		0.1;	//	value is in %	//	double

		//	Shade Settings

		Zone.ShadeShow[GroupIndex] =					FALSE;						//	Boolean_t
		Zone.ShadeColor[GroupIndex] =				Black_C;					//	ColorIndex_t
		Zone.ShadeLighting[GroupIndex] =				TRUE;						//	Boolean_t

		//	Edge Settings

		Zone.EdgeShow[GroupIndex] =					FALSE;						//	Boolean_t
		Zone.EdgeType[GroupIndex] =					EdgeType_BordersAndCreases;	//	EdgeType_e
		Zone.EdgeIBorder[GroupIndex] =				BorderLocation_Both;		//	BorderLocation_e
		Zone.EdgeJBorder[GroupIndex] =				BorderLocation_Both;		//	BorderLocation_e
		Zone.EdgeKBorder[GroupIndex] =				BorderLocation_Both;		//	BorderLocation_e
		Zone.EdgeColor[GroupIndex] =					Black_C;					//	ColorIndex_t
		Zone.EdgeLineThickness[GroupIndex] =			0.1;	//	value is in %	//	double

		//	Points Settings

		Zone.PointsToPlot[GroupIndex] =				PointsToPlot_SurfaceNodes;	//	PointsToPlot_e
		Zone.PointsIndexSkipI[GroupIndex] =			1;							//	LgIndex_t
		Zone.PointsIndexSkipJ[GroupIndex] =			1;							//	LgIndex_t
		Zone.PointsIndexSkipK[GroupIndex] =			1;							//	LgIndex_t
		/*
			1 is the same as "No Skip" in the GUI.
		*/

		//	Surfaces Settings

		Zone.SurfacesToPlot[GroupIndex] =			SurfacesToPlot_ExposedCellFaces;	//	SurfacesToPlot_e
		Zone.RangeBeginI[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeEndI[GroupIndex] =					0;							//	LgIndex_t
		Zone.RangeSkipI[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeBeginJ[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeEndJ[GroupIndex] =					0;							//	LgIndex_t
		Zone.RangeSkipJ[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeBeginK[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeEndK[GroupIndex] =					0;							//	LgIndex_t
		Zone.RangeSkipK[GroupIndex] =				1;							//	LgIndex_t
		/*
			0 is the same as maximum, which sets the end range to the maximum
			index for the dimension. The skip number is to skip every # planes,
			so 1 is to not skip any planes, (show every 1 plane).
		*/

		//	Volume Settings

		Zone.StreamtracesShow[GroupIndex] =			TRUE;						//	Boolean_t
		Zone.IsoSurfacesShow[GroupIndex] =			TRUE;						//	Boolean_t
		Zone.SlicesShow[GroupIndex] =				TRUE;						//	Boolean_t

		//	Effects Settings

		Zone.UseValueBlanking[GroupIndex] =			TRUE;						//	Boolean_t
		//Zone.UseClipPlanes[GroupIndex] =			TRUE;						//	Boolean_t
		Zone.UseSurfaceTranslucency[GroupIndex] =	FALSE;						//	Boolean_t
		Zone.SurfaceTranslucency[GroupIndex] =		50;		//	value is in %	//	SmInteger_t
		Zone.LightingEffect[GroupIndex] =			LightingEffect_Gouraud;		//	LightingEffect_e

		//	Group index and number (increments for next set)

		GroupIndex++;														//	int


		//===============================================================
		//	ACLines Settings
		//===============================================================

		Zone.ZoneShow[GroupIndex] =					FALSE;						//	Boolean_t

		//	Mesh Settings

		Zone.MeshShow[GroupIndex] =					TRUE;						//	Boolean_t
		Zone.MeshType[GroupIndex] =					MeshType_Overlay;			//	MeshType_e
		Zone.MeshColor[GroupIndex] =					MultiColor6_C;					//	ColorIndex_t
		/*
			A number of color index constants are \#defined. These include Black_C, Red_C,
			Green_C, Blue_C, Cyan_C, Yellow_C, Purple_C, White_C, Custom1_C through
			Custom56_C, MultiColor_C, NoColor_C, MulitiColor2_C through MulitiColor8_C,
			RGBColor_C, and InvalidColor_C.  
		*/
		Zone.MeshLinePttrn[GroupIndex] =				LinePattern_Dashed;			//	LinePattern_e
		Zone.MeshLinePttrnLen[GroupIndex] =			2;		//	Value is in %	//	double
		Zone.MeshLineThickness[GroupIndex] =			0.2;	//	Value is in %	//	double


		//	Contour Settings

		Zone.ContourShow[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.ContourType[GroupIndex] =				ContourType_Flood;			//	ContourType_e
		Zone.FloodColoring[GroupIndex] =				ContourColoring_Group1;		//	ContourColoring_e
		Zone.ContourLineGroup[GroupIndex] =			1;							//	SmInteger_t
		Zone.ContourLineColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.ContourLinePttrn[GroupIndex] =			LinePattern_Solid;			//	LinePattern_e
		Zone.ContourLinePttrnLen[GroupIndex] =		2;		//	Value is in %	//	double
		Zone.ContourLineThickness[GroupIndex] =		0.1;	//	Value is in %	//	double
		Zone.ContourLighting[GroupIndex] =			TRUE;						//	Boolean_t

		//	Vector Settings
		
		Zone.VectorShow[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.VectorType[GroupIndex] =				VectorType_TailAtPoint;		//	VectorType_e
		Zone.VectorArrowType[GroupIndex] =			ArrowheadStyle_Plain;		//	ArrowheadStyle_e
		Zone.VectorLineColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.VectorTangent[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.VectorLinePttrn[GroupIndex] =			LinePattern_Solid;			//	LinePattern_e
		Zone.VectorLinePttrnLen[GroupIndex] =		2;		//	Value is in %	//	double
		Zone.VectorLineThickness[GroupIndex] =		0.1;	//	Value is in %	//	double

		//	Scatter Settings

		Zone.ScatterShow[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.ScatterShape[GroupIndex] =				GeomShape_Sphere;			//	GeomShape_e
		Zone.ScatterOutColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.ScatterInColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.ScatterFillMode[GroupIndex] =			FillMode_UseLineColor;		//	FillMode_e
		Zone.ScatterSize[GroupIndex] =				2.5;	//	value is in %	//	double
		Zone.ScatterLineThickness[GroupIndex] =		0.1;	//	value is in %	//	double

		//	Shade Settings

		Zone.ShadeShow[GroupIndex] =					FALSE;						//	Boolean_t
		Zone.ShadeColor[GroupIndex] =				Black_C;					//	ColorIndex_t
		Zone.ShadeLighting[GroupIndex] =				TRUE;						//	Boolean_t

		//	Edge Settings

		Zone.EdgeShow[GroupIndex] =					FALSE;						//	Boolean_t
		Zone.EdgeType[GroupIndex] =					EdgeType_BordersAndCreases;	//	EdgeType_e
		Zone.EdgeIBorder[GroupIndex] =				BorderLocation_Both;		//	BorderLocation_e
		Zone.EdgeJBorder[GroupIndex] =				BorderLocation_Both;		//	BorderLocation_e
		Zone.EdgeKBorder[GroupIndex] =				BorderLocation_Both;		//	BorderLocation_e
		Zone.EdgeColor[GroupIndex] =					Black_C;					//	ColorIndex_t
		Zone.EdgeLineThickness[GroupIndex] =			0.1;	//	value is in %	//	double

		//	Points Settings

		Zone.PointsToPlot[GroupIndex] =				PointsToPlot_SurfaceNodes;	//	PointsToPlot_e
		Zone.PointsIndexSkipI[GroupIndex] =			1;							//	LgIndex_t
		Zone.PointsIndexSkipJ[GroupIndex] =			1;							//	LgIndex_t
		Zone.PointsIndexSkipK[GroupIndex] =			1;							//	LgIndex_t
		/*
			1 is the same as "No Skip" in the GUI.
		*/

		//	Surfaces Settings

		Zone.SurfacesToPlot[GroupIndex] =			SurfacesToPlot_ExposedCellFaces;	//	SurfacesToPlot_e
		Zone.RangeBeginI[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeEndI[GroupIndex] =					0;							//	LgIndex_t
		Zone.RangeSkipI[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeBeginJ[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeEndJ[GroupIndex] =					0;							//	LgIndex_t
		Zone.RangeSkipJ[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeBeginK[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeEndK[GroupIndex] =					0;							//	LgIndex_t
		Zone.RangeSkipK[GroupIndex] =				1;							//	LgIndex_t
		/*
			0 is the same as maximum, which sets the end range to the maximum
			index for the dimension. The skip number is to skip every # planes,
			so 1 is to not skip any planes, (show every 1 plane).
		*/

		//	Volume Settings

		Zone.StreamtracesShow[GroupIndex] =			TRUE;						//	Boolean_t
		Zone.IsoSurfacesShow[GroupIndex] =			TRUE;						//	Boolean_t
		Zone.SlicesShow[GroupIndex] =				TRUE;						//	Boolean_t

		//	Effects Settings

		Zone.UseValueBlanking[GroupIndex] =			TRUE;						//	Boolean_t
		//Zone.UseClipPlanes[GroupIndex] =			TRUE;						//	Boolean_t
		Zone.UseSurfaceTranslucency[GroupIndex] =	FALSE;						//	Boolean_t
		Zone.SurfaceTranslucency[GroupIndex] =		50;		//	value is in %	//	SmInteger_t
		Zone.LightingEffect[GroupIndex] =			LightingEffect_Gouraud;		//	LightingEffect_e

		//	Group index and number (increments for next set)

		GroupIndex++;														//	int


		//===============================================================
		//	ACLinesFF Settings
		//===============================================================

		Zone.ZoneShow[GroupIndex] =					FALSE;						//	Boolean_t

		//	Mesh Settings

		Zone.MeshShow[GroupIndex] =					TRUE;						//	Boolean_t
		Zone.MeshType[GroupIndex] =					MeshType_Overlay;			//	MeshType_e
		Zone.MeshColor[GroupIndex] =					MultiColor6_C;					//	ColorIndex_t
		/*
			A number of color index constants are \#defined. These include Black_C, Red_C,
			Green_C, Blue_C, Cyan_C, Yellow_C, Purple_C, White_C, Custom1_C through
			Custom56_C, MultiColor_C, NoColor_C, MulitiColor2_C through MulitiColor8_C,
			RGBColor_C, and InvalidColor_C.  
		*/
		Zone.MeshLinePttrn[GroupIndex] =				LinePattern_Dashed;			//	LinePattern_e
		Zone.MeshLinePttrnLen[GroupIndex] =			2;		//	Value is in %	//	double
		Zone.MeshLineThickness[GroupIndex] =			0.2;	//	Value is in %	//	double


		//	Contour Settings

		Zone.ContourShow[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.ContourType[GroupIndex] =				ContourType_Flood;			//	ContourType_e
		Zone.FloodColoring[GroupIndex] =				ContourColoring_Group1;		//	ContourColoring_e
		Zone.ContourLineGroup[GroupIndex] =			1;							//	SmInteger_t
		Zone.ContourLineColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.ContourLinePttrn[GroupIndex] =			LinePattern_Solid;			//	LinePattern_e
		Zone.ContourLinePttrnLen[GroupIndex] =		2;		//	Value is in %	//	double
		Zone.ContourLineThickness[GroupIndex] =		0.1;	//	Value is in %	//	double
		Zone.ContourLighting[GroupIndex] =			TRUE;						//	Boolean_t

		//	Vector Settings
		
		Zone.VectorShow[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.VectorType[GroupIndex] =				VectorType_TailAtPoint;		//	VectorType_e
		Zone.VectorArrowType[GroupIndex] =			ArrowheadStyle_Plain;		//	ArrowheadStyle_e
		Zone.VectorLineColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.VectorTangent[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.VectorLinePttrn[GroupIndex] =			LinePattern_Solid;			//	LinePattern_e
		Zone.VectorLinePttrnLen[GroupIndex] =		2;		//	Value is in %	//	double
		Zone.VectorLineThickness[GroupIndex] =		0.1;	//	Value is in %	//	double

		//	Scatter Settings

		Zone.ScatterShow[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.ScatterShape[GroupIndex] =				GeomShape_Sphere;			//	GeomShape_e
		Zone.ScatterOutColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.ScatterInColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.ScatterFillMode[GroupIndex] =			FillMode_UseLineColor;		//	FillMode_e
		Zone.ScatterSize[GroupIndex] =				2.5;	//	value is in %	//	double
		Zone.ScatterLineThickness[GroupIndex] =		0.1;	//	value is in %	//	double

		//	Shade Settings

		Zone.ShadeShow[GroupIndex] =					FALSE;						//	Boolean_t
		Zone.ShadeColor[GroupIndex] =				Black_C;					//	ColorIndex_t
		Zone.ShadeLighting[GroupIndex] =				TRUE;						//	Boolean_t

		//	Edge Settings

		Zone.EdgeShow[GroupIndex] =					FALSE;						//	Boolean_t
		Zone.EdgeType[GroupIndex] =					EdgeType_BordersAndCreases;	//	EdgeType_e
		Zone.EdgeIBorder[GroupIndex] =				BorderLocation_Both;		//	BorderLocation_e
		Zone.EdgeJBorder[GroupIndex] =				BorderLocation_Both;		//	BorderLocation_e
		Zone.EdgeKBorder[GroupIndex] =				BorderLocation_Both;		//	BorderLocation_e
		Zone.EdgeColor[GroupIndex] =					Black_C;					//	ColorIndex_t
		Zone.EdgeLineThickness[GroupIndex] =			0.1;	//	value is in %	//	double

		//	Points Settings

		Zone.PointsToPlot[GroupIndex] =				PointsToPlot_SurfaceNodes;	//	PointsToPlot_e
		Zone.PointsIndexSkipI[GroupIndex] =			1;							//	LgIndex_t
		Zone.PointsIndexSkipJ[GroupIndex] =			1;							//	LgIndex_t
		Zone.PointsIndexSkipK[GroupIndex] =			1;							//	LgIndex_t
		/*
			1 is the same as "No Skip" in the GUI.
		*/

		//	Surfaces Settings

		Zone.SurfacesToPlot[GroupIndex] =			SurfacesToPlot_ExposedCellFaces;	//	SurfacesToPlot_e
		Zone.RangeBeginI[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeEndI[GroupIndex] =					0;							//	LgIndex_t
		Zone.RangeSkipI[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeBeginJ[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeEndJ[GroupIndex] =					0;							//	LgIndex_t
		Zone.RangeSkipJ[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeBeginK[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeEndK[GroupIndex] =					0;							//	LgIndex_t
		Zone.RangeSkipK[GroupIndex] =				1;							//	LgIndex_t
		/*
			0 is the same as maximum, which sets the end range to the maximum
			index for the dimension. The skip number is to skip every # planes,
			so 1 is to not skip any planes, (show every 1 plane).
		*/

		//	Volume Settings

		Zone.StreamtracesShow[GroupIndex] =			TRUE;						//	Boolean_t
		Zone.IsoSurfacesShow[GroupIndex] =			TRUE;						//	Boolean_t
		Zone.SlicesShow[GroupIndex] =				TRUE;						//	Boolean_t

		//	Effects Settings

		Zone.UseValueBlanking[GroupIndex] =			TRUE;						//	Boolean_t
		//Zone.UseClipPlanes[GroupIndex] =			TRUE;						//	Boolean_t
		Zone.UseSurfaceTranslucency[GroupIndex] =	FALSE;						//	Boolean_t
		Zone.SurfaceTranslucency[GroupIndex] =		50;		//	value is in %	//	SmInteger_t
		Zone.LightingEffect[GroupIndex] =			LightingEffect_Gouraud;		//	LightingEffect_e

		//	Group index and number (increments for next set)

		GroupIndex++;														//	int


		//===============================================================
		//	BRLines Settings
		//===============================================================

		Zone.ZoneShow[GroupIndex] =					FALSE;						//	Boolean_t

		//	Mesh Settings

		Zone.MeshShow[GroupIndex] =					TRUE;						//	Boolean_t
		Zone.MeshType[GroupIndex] =					MeshType_Overlay;			//	MeshType_e
		Zone.MeshColor[GroupIndex] =					MultiColor7_C;					//	ColorIndex_t
		/*
			A number of color index constants are \#defined. These include Black_C, Red_C,
			Green_C, Blue_C, Cyan_C, Yellow_C, Purple_C, White_C, Custom1_C through
			Custom56_C, MultiColor_C, NoColor_C, MulitiColor2_C through MulitiColor8_C,
			RGBColor_C, and InvalidColor_C.  
		*/
		Zone.MeshLinePttrn[GroupIndex] =				LinePattern_DashDot;			//	LinePattern_e
		Zone.MeshLinePttrnLen[GroupIndex] =			2;		//	Value is in %	//	double
		Zone.MeshLineThickness[GroupIndex] =			0.2;	//	Value is in %	//	double


		//	Contour Settings

		Zone.ContourShow[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.ContourType[GroupIndex] =				ContourType_Flood;			//	ContourType_e
		Zone.FloodColoring[GroupIndex] =				ContourColoring_Group1;		//	ContourColoring_e
		Zone.ContourLineGroup[GroupIndex] =			1;							//	SmInteger_t
		Zone.ContourLineColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.ContourLinePttrn[GroupIndex] =			LinePattern_Solid;			//	LinePattern_e
		Zone.ContourLinePttrnLen[GroupIndex] =		2;		//	Value is in %	//	double
		Zone.ContourLineThickness[GroupIndex] =		0.1;	//	Value is in %	//	double
		Zone.ContourLighting[GroupIndex] =			TRUE;						//	Boolean_t

		//	Vector Settings
		
		Zone.VectorShow[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.VectorType[GroupIndex] =				VectorType_TailAtPoint;		//	VectorType_e
		Zone.VectorArrowType[GroupIndex] =			ArrowheadStyle_Plain;		//	ArrowheadStyle_e
		Zone.VectorLineColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.VectorTangent[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.VectorLinePttrn[GroupIndex] =			LinePattern_Solid;			//	LinePattern_e
		Zone.VectorLinePttrnLen[GroupIndex] =		2;		//	Value is in %	//	double
		Zone.VectorLineThickness[GroupIndex] =		0.1;	//	Value is in %	//	double

		//	Scatter Settings

		Zone.ScatterShow[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.ScatterShape[GroupIndex] =				GeomShape_Sphere;			//	GeomShape_e
		Zone.ScatterOutColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.ScatterInColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.ScatterFillMode[GroupIndex] =			FillMode_UseLineColor;		//	FillMode_e
		Zone.ScatterSize[GroupIndex] =				2.5;	//	value is in %	//	double
		Zone.ScatterLineThickness[GroupIndex] =		0.1;	//	value is in %	//	double

		//	Shade Settings

		Zone.ShadeShow[GroupIndex] =					FALSE;						//	Boolean_t
		Zone.ShadeColor[GroupIndex] =				Black_C;					//	ColorIndex_t
		Zone.ShadeLighting[GroupIndex] =				TRUE;						//	Boolean_t

		//	Edge Settings

		Zone.EdgeShow[GroupIndex] =					FALSE;						//	Boolean_t
		Zone.EdgeType[GroupIndex] =					EdgeType_BordersAndCreases;	//	EdgeType_e
		Zone.EdgeIBorder[GroupIndex] =				BorderLocation_Both;		//	BorderLocation_e
		Zone.EdgeJBorder[GroupIndex] =				BorderLocation_Both;		//	BorderLocation_e
		Zone.EdgeKBorder[GroupIndex] =				BorderLocation_Both;		//	BorderLocation_e
		Zone.EdgeColor[GroupIndex] =					Black_C;					//	ColorIndex_t
		Zone.EdgeLineThickness[GroupIndex] =			0.1;	//	value is in %	//	double

		//	Points Settings

		Zone.PointsToPlot[GroupIndex] =				PointsToPlot_SurfaceNodes;	//	PointsToPlot_e
		Zone.PointsIndexSkipI[GroupIndex] =			1;							//	LgIndex_t
		Zone.PointsIndexSkipJ[GroupIndex] =			1;							//	LgIndex_t
		Zone.PointsIndexSkipK[GroupIndex] =			1;							//	LgIndex_t
		/*
			1 is the same as "No Skip" in the GUI.
		*/

		//	Surfaces Settings

		Zone.SurfacesToPlot[GroupIndex] =			SurfacesToPlot_ExposedCellFaces;	//	SurfacesToPlot_e
		Zone.RangeBeginI[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeEndI[GroupIndex] =					0;							//	LgIndex_t
		Zone.RangeSkipI[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeBeginJ[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeEndJ[GroupIndex] =					0;							//	LgIndex_t
		Zone.RangeSkipJ[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeBeginK[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeEndK[GroupIndex] =					0;							//	LgIndex_t
		Zone.RangeSkipK[GroupIndex] =				1;							//	LgIndex_t
		/*
			0 is the same as maximum, which sets the end range to the maximum
			index for the dimension. The skip number is to skip every # planes,
			so 1 is to not skip any planes, (show every 1 plane).
		*/

		//	Volume Settings

		Zone.StreamtracesShow[GroupIndex] =			TRUE;						//	Boolean_t
		Zone.IsoSurfacesShow[GroupIndex] =			TRUE;						//	Boolean_t
		Zone.SlicesShow[GroupIndex] =				TRUE;						//	Boolean_t

		//	Effects Settings

		Zone.UseValueBlanking[GroupIndex] =			TRUE;						//	Boolean_t
		//Zone.UseClipPlanes[GroupIndex] =			TRUE;						//	Boolean_t
		Zone.UseSurfaceTranslucency[GroupIndex] =	FALSE;						//	Boolean_t
		Zone.SurfaceTranslucency[GroupIndex] =		50;		//	value is in %	//	SmInteger_t
		Zone.LightingEffect[GroupIndex] =			LightingEffect_Gouraud;		//	LightingEffect_e

		//	Group index and number (increments for next set)

		GroupIndex++;														//	int


		//===============================================================
		//	BRLinesFF Settings
		//===============================================================

		Zone.ZoneShow[GroupIndex] =					FALSE;						//	Boolean_t

		//	Mesh Settings

		Zone.MeshShow[GroupIndex] =					TRUE;						//	Boolean_t
		Zone.MeshType[GroupIndex] =					MeshType_Overlay;			//	MeshType_e
		Zone.MeshColor[GroupIndex] =					MultiColor7_C;					//	ColorIndex_t
		/*
			A number of color index constants are \#defined. These include Black_C, Red_C,
			Green_C, Blue_C, Cyan_C, Yellow_C, Purple_C, White_C, Custom1_C through
			Custom56_C, MultiColor_C, NoColor_C, MulitiColor2_C through MulitiColor8_C,
			RGBColor_C, and InvalidColor_C.  
		*/
		Zone.MeshLinePttrn[GroupIndex] =				LinePattern_DashDot;			//	LinePattern_e
		Zone.MeshLinePttrnLen[GroupIndex] =			2;		//	Value is in %	//	double
		Zone.MeshLineThickness[GroupIndex] =			0.2;	//	Value is in %	//	double


		//	Contour Settings

		Zone.ContourShow[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.ContourType[GroupIndex] =				ContourType_Flood;			//	ContourType_e
		Zone.FloodColoring[GroupIndex] =				ContourColoring_Group1;		//	ContourColoring_e
		Zone.ContourLineGroup[GroupIndex] =			1;							//	SmInteger_t
		Zone.ContourLineColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.ContourLinePttrn[GroupIndex] =			LinePattern_Solid;			//	LinePattern_e
		Zone.ContourLinePttrnLen[GroupIndex] =		2;		//	Value is in %	//	double
		Zone.ContourLineThickness[GroupIndex] =		0.1;	//	Value is in %	//	double
		Zone.ContourLighting[GroupIndex] =			TRUE;						//	Boolean_t

		//	Vector Settings
		
		Zone.VectorShow[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.VectorType[GroupIndex] =				VectorType_TailAtPoint;		//	VectorType_e
		Zone.VectorArrowType[GroupIndex] =			ArrowheadStyle_Plain;		//	ArrowheadStyle_e
		Zone.VectorLineColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.VectorTangent[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.VectorLinePttrn[GroupIndex] =			LinePattern_Solid;			//	LinePattern_e
		Zone.VectorLinePttrnLen[GroupIndex] =		2;		//	Value is in %	//	double
		Zone.VectorLineThickness[GroupIndex] =		0.1;	//	Value is in %	//	double

		//	Scatter Settings

		Zone.ScatterShow[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.ScatterShape[GroupIndex] =				GeomShape_Sphere;			//	GeomShape_e
		Zone.ScatterOutColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.ScatterInColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.ScatterFillMode[GroupIndex] =			FillMode_UseLineColor;		//	FillMode_e
		Zone.ScatterSize[GroupIndex] =				2.5;	//	value is in %	//	double
		Zone.ScatterLineThickness[GroupIndex] =		0.1;	//	value is in %	//	double

		//	Shade Settings

		Zone.ShadeShow[GroupIndex] =					FALSE;						//	Boolean_t
		Zone.ShadeColor[GroupIndex] =				Black_C;					//	ColorIndex_t
		Zone.ShadeLighting[GroupIndex] =				TRUE;						//	Boolean_t

		//	Edge Settings

		Zone.EdgeShow[GroupIndex] =					FALSE;						//	Boolean_t
		Zone.EdgeType[GroupIndex] =					EdgeType_BordersAndCreases;	//	EdgeType_e
		Zone.EdgeIBorder[GroupIndex] =				BorderLocation_Both;		//	BorderLocation_e
		Zone.EdgeJBorder[GroupIndex] =				BorderLocation_Both;		//	BorderLocation_e
		Zone.EdgeKBorder[GroupIndex] =				BorderLocation_Both;		//	BorderLocation_e
		Zone.EdgeColor[GroupIndex] =					Black_C;					//	ColorIndex_t
		Zone.EdgeLineThickness[GroupIndex] =			0.1;	//	value is in %	//	double

		//	Points Settings

		Zone.PointsToPlot[GroupIndex] =				PointsToPlot_SurfaceNodes;	//	PointsToPlot_e
		Zone.PointsIndexSkipI[GroupIndex] =			1;							//	LgIndex_t
		Zone.PointsIndexSkipJ[GroupIndex] =			1;							//	LgIndex_t
		Zone.PointsIndexSkipK[GroupIndex] =			1;							//	LgIndex_t
		/*
			1 is the same as "No Skip" in the GUI.
		*/

		//	Surfaces Settings

		Zone.SurfacesToPlot[GroupIndex] =			SurfacesToPlot_ExposedCellFaces;	//	SurfacesToPlot_e
		Zone.RangeBeginI[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeEndI[GroupIndex] =					0;							//	LgIndex_t
		Zone.RangeSkipI[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeBeginJ[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeEndJ[GroupIndex] =					0;							//	LgIndex_t
		Zone.RangeSkipJ[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeBeginK[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeEndK[GroupIndex] =					0;							//	LgIndex_t
		Zone.RangeSkipK[GroupIndex] =				1;							//	LgIndex_t
		/*
			0 is the same as maximum, which sets the end range to the maximum
			index for the dimension. The skip number is to skip every # planes,
			so 1 is to not skip any planes, (show every 1 plane).
		*/

		//	Volume Settings

		Zone.StreamtracesShow[GroupIndex] =			TRUE;						//	Boolean_t
		Zone.IsoSurfacesShow[GroupIndex] =			TRUE;						//	Boolean_t
		Zone.SlicesShow[GroupIndex] =				TRUE;						//	Boolean_t

		//	Effects Settings

		Zone.UseValueBlanking[GroupIndex] =			TRUE;						//	Boolean_t
		//Zone.UseClipPlanes[GroupIndex] =			TRUE;						//	Boolean_t
		Zone.UseSurfaceTranslucency[GroupIndex] =	FALSE;						//	Boolean_t
		Zone.SurfaceTranslucency[GroupIndex] =		50;		//	value is in %	//	SmInteger_t
		Zone.LightingEffect[GroupIndex] =			LightingEffect_Gouraud;		//	LightingEffect_e

		//	Group index and number (increments for next set)

		GroupIndex++;														//	int


		//===============================================================
		//	BCLines Settings
		//===============================================================

		Zone.ZoneShow[GroupIndex] =					FALSE;						//	Boolean_t

		//	Mesh Settings

		Zone.MeshShow[GroupIndex] =					TRUE;						//	Boolean_t
		Zone.MeshType[GroupIndex] =					MeshType_Overlay;			//	MeshType_e
		Zone.MeshColor[GroupIndex] =					MultiColor8_C;					//	ColorIndex_t
		/*
			A number of color index constants are \#defined. These include Black_C, Red_C,
			Green_C, Blue_C, Cyan_C, Yellow_C, Purple_C, White_C, Custom1_C through
			Custom56_C, MultiColor_C, NoColor_C, MulitiColor2_C through MulitiColor8_C,
			RGBColor_C, and InvalidColor_C.  
		*/
		Zone.MeshLinePttrn[GroupIndex] =				LinePattern_DashDotDot;			//	LinePattern_e
		Zone.MeshLinePttrnLen[GroupIndex] =			2;		//	Value is in %	//	double
		Zone.MeshLineThickness[GroupIndex] =			0.2;	//	Value is in %	//	double


		//	Contour Settings

		Zone.ContourShow[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.ContourType[GroupIndex] =				ContourType_Flood;			//	ContourType_e
		Zone.FloodColoring[GroupIndex] =				ContourColoring_Group1;		//	ContourColoring_e
		Zone.ContourLineGroup[GroupIndex] =			1;							//	SmInteger_t
		Zone.ContourLineColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.ContourLinePttrn[GroupIndex] =			LinePattern_Solid;			//	LinePattern_e
		Zone.ContourLinePttrnLen[GroupIndex] =		2;		//	Value is in %	//	double
		Zone.ContourLineThickness[GroupIndex] =		0.1;	//	Value is in %	//	double
		Zone.ContourLighting[GroupIndex] =			TRUE;						//	Boolean_t

		//	Vector Settings
		
		Zone.VectorShow[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.VectorType[GroupIndex] =				VectorType_TailAtPoint;		//	VectorType_e
		Zone.VectorArrowType[GroupIndex] =			ArrowheadStyle_Plain;		//	ArrowheadStyle_e
		Zone.VectorLineColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.VectorTangent[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.VectorLinePttrn[GroupIndex] =			LinePattern_Solid;			//	LinePattern_e
		Zone.VectorLinePttrnLen[GroupIndex] =		2;		//	Value is in %	//	double
		Zone.VectorLineThickness[GroupIndex] =		0.1;	//	Value is in %	//	double

		//	Scatter Settings

		Zone.ScatterShow[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.ScatterShape[GroupIndex] =				GeomShape_Sphere;			//	GeomShape_e
		Zone.ScatterOutColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.ScatterInColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.ScatterFillMode[GroupIndex] =			FillMode_UseLineColor;		//	FillMode_e
		Zone.ScatterSize[GroupIndex] =				2.5;	//	value is in %	//	double
		Zone.ScatterLineThickness[GroupIndex] =		0.1;	//	value is in %	//	double

		//	Shade Settings

		Zone.ShadeShow[GroupIndex] =					FALSE;						//	Boolean_t
		Zone.ShadeColor[GroupIndex] =				Black_C;					//	ColorIndex_t
		Zone.ShadeLighting[GroupIndex] =				TRUE;						//	Boolean_t

		//	Edge Settings

		Zone.EdgeShow[GroupIndex] =					FALSE;						//	Boolean_t
		Zone.EdgeType[GroupIndex] =					EdgeType_BordersAndCreases;	//	EdgeType_e
		Zone.EdgeIBorder[GroupIndex] =				BorderLocation_Both;		//	BorderLocation_e
		Zone.EdgeJBorder[GroupIndex] =				BorderLocation_Both;		//	BorderLocation_e
		Zone.EdgeKBorder[GroupIndex] =				BorderLocation_Both;		//	BorderLocation_e
		Zone.EdgeColor[GroupIndex] =					Black_C;					//	ColorIndex_t
		Zone.EdgeLineThickness[GroupIndex] =			0.1;	//	value is in %	//	double

		//	Points Settings

		Zone.PointsToPlot[GroupIndex] =				PointsToPlot_SurfaceNodes;	//	PointsToPlot_e
		Zone.PointsIndexSkipI[GroupIndex] =			1;							//	LgIndex_t
		Zone.PointsIndexSkipJ[GroupIndex] =			1;							//	LgIndex_t
		Zone.PointsIndexSkipK[GroupIndex] =			1;							//	LgIndex_t
		/*
			1 is the same as "No Skip" in the GUI.
		*/

		//	Surfaces Settings

		Zone.SurfacesToPlot[GroupIndex] =			SurfacesToPlot_ExposedCellFaces;	//	SurfacesToPlot_e
		Zone.RangeBeginI[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeEndI[GroupIndex] =					0;							//	LgIndex_t
		Zone.RangeSkipI[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeBeginJ[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeEndJ[GroupIndex] =					0;							//	LgIndex_t
		Zone.RangeSkipJ[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeBeginK[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeEndK[GroupIndex] =					0;							//	LgIndex_t
		Zone.RangeSkipK[GroupIndex] =				1;							//	LgIndex_t
		/*
			0 is the same as maximum, which sets the end range to the maximum
			index for the dimension. The skip number is to skip every # planes,
			so 1 is to not skip any planes, (show every 1 plane).
		*/

		//	Volume Settings

		Zone.StreamtracesShow[GroupIndex] =			TRUE;						//	Boolean_t
		Zone.IsoSurfacesShow[GroupIndex] =			TRUE;						//	Boolean_t
		Zone.SlicesShow[GroupIndex] =				TRUE;						//	Boolean_t

		//	Effects Settings

		Zone.UseValueBlanking[GroupIndex] =			TRUE;						//	Boolean_t
		//Zone.UseClipPlanes[GroupIndex] =			TRUE;						//	Boolean_t
		Zone.UseSurfaceTranslucency[GroupIndex] =	FALSE;						//	Boolean_t
		Zone.SurfaceTranslucency[GroupIndex] =		50;		//	value is in %	//	SmInteger_t
		Zone.LightingEffect[GroupIndex] =			LightingEffect_Gouraud;		//	LightingEffect_e

		//	Group index and number (increments for next set)

		GroupIndex++;														//	int


		//===============================================================
		//	BCLinesFF Settings
		//===============================================================

		Zone.ZoneShow[GroupIndex] =					FALSE;						//	Boolean_t

		//	Mesh Settings

		Zone.MeshShow[GroupIndex] =					TRUE;						//	Boolean_t
		Zone.MeshType[GroupIndex] =					MeshType_Overlay;			//	MeshType_e
		Zone.MeshColor[GroupIndex] =					MultiColor8_C;					//	ColorIndex_t
		/*
			A number of color index constants are \#defined. These include Black_C, Red_C,
			Green_C, Blue_C, Cyan_C, Yellow_C, Purple_C, White_C, Custom1_C through
			Custom56_C, MultiColor_C, NoColor_C, MulitiColor2_C through MulitiColor8_C,
			RGBColor_C, and InvalidColor_C.  
		*/
		Zone.MeshLinePttrn[GroupIndex] =				LinePattern_DashDotDot;			//	LinePattern_e
		Zone.MeshLinePttrnLen[GroupIndex] =			2;		//	Value is in %	//	double
		Zone.MeshLineThickness[GroupIndex] =			0.2;	//	Value is in %	//	double


		//	Contour Settings

		Zone.ContourShow[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.ContourType[GroupIndex] =				ContourType_Flood;			//	ContourType_e
		Zone.FloodColoring[GroupIndex] =				ContourColoring_Group1;		//	ContourColoring_e
		Zone.ContourLineGroup[GroupIndex] =			1;							//	SmInteger_t
		Zone.ContourLineColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.ContourLinePttrn[GroupIndex] =			LinePattern_Solid;			//	LinePattern_e
		Zone.ContourLinePttrnLen[GroupIndex] =		2;		//	Value is in %	//	double
		Zone.ContourLineThickness[GroupIndex] =		0.1;	//	Value is in %	//	double
		Zone.ContourLighting[GroupIndex] =			TRUE;						//	Boolean_t

		//	Vector Settings
		
		Zone.VectorShow[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.VectorType[GroupIndex] =				VectorType_TailAtPoint;		//	VectorType_e
		Zone.VectorArrowType[GroupIndex] =			ArrowheadStyle_Plain;		//	ArrowheadStyle_e
		Zone.VectorLineColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.VectorTangent[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.VectorLinePttrn[GroupIndex] =			LinePattern_Solid;			//	LinePattern_e
		Zone.VectorLinePttrnLen[GroupIndex] =		2;		//	Value is in %	//	double
		Zone.VectorLineThickness[GroupIndex] =		0.1;	//	Value is in %	//	double

		//	Scatter Settings

		Zone.ScatterShow[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.ScatterShape[GroupIndex] =				GeomShape_Sphere;			//	GeomShape_e
		Zone.ScatterOutColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.ScatterInColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.ScatterFillMode[GroupIndex] =			FillMode_UseLineColor;		//	FillMode_e
		Zone.ScatterSize[GroupIndex] =				2.5;	//	value is in %	//	double
		Zone.ScatterLineThickness[GroupIndex] =		0.1;	//	value is in %	//	double

		//	Shade Settings

		Zone.ShadeShow[GroupIndex] =					FALSE;						//	Boolean_t
		Zone.ShadeColor[GroupIndex] =				Black_C;					//	ColorIndex_t
		Zone.ShadeLighting[GroupIndex] =				TRUE;						//	Boolean_t

		//	Edge Settings

		Zone.EdgeShow[GroupIndex] =					FALSE;						//	Boolean_t
		Zone.EdgeType[GroupIndex] =					EdgeType_BordersAndCreases;	//	EdgeType_e
		Zone.EdgeIBorder[GroupIndex] =				BorderLocation_Both;		//	BorderLocation_e
		Zone.EdgeJBorder[GroupIndex] =				BorderLocation_Both;		//	BorderLocation_e
		Zone.EdgeKBorder[GroupIndex] =				BorderLocation_Both;		//	BorderLocation_e
		Zone.EdgeColor[GroupIndex] =					Black_C;					//	ColorIndex_t
		Zone.EdgeLineThickness[GroupIndex] =			0.1;	//	value is in %	//	double

		//	Points Settings

		Zone.PointsToPlot[GroupIndex] =				PointsToPlot_SurfaceNodes;	//	PointsToPlot_e
		Zone.PointsIndexSkipI[GroupIndex] =			1;							//	LgIndex_t
		Zone.PointsIndexSkipJ[GroupIndex] =			1;							//	LgIndex_t
		Zone.PointsIndexSkipK[GroupIndex] =			1;							//	LgIndex_t
		/*
			1 is the same as "No Skip" in the GUI.
		*/

		//	Surfaces Settings

		Zone.SurfacesToPlot[GroupIndex] =			SurfacesToPlot_ExposedCellFaces;	//	SurfacesToPlot_e
		Zone.RangeBeginI[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeEndI[GroupIndex] =					0;							//	LgIndex_t
		Zone.RangeSkipI[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeBeginJ[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeEndJ[GroupIndex] =					0;							//	LgIndex_t
		Zone.RangeSkipJ[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeBeginK[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeEndK[GroupIndex] =					0;							//	LgIndex_t
		Zone.RangeSkipK[GroupIndex] =				1;							//	LgIndex_t
		/*
			0 is the same as maximum, which sets the end range to the maximum
			index for the dimension. The skip number is to skip every # planes,
			so 1 is to not skip any planes, (show every 1 plane).
		*/

		//	Volume Settings

		Zone.StreamtracesShow[GroupIndex] =			TRUE;						//	Boolean_t
		Zone.IsoSurfacesShow[GroupIndex] =			TRUE;						//	Boolean_t
		Zone.SlicesShow[GroupIndex] =				TRUE;						//	Boolean_t

		//	Effects Settings

		Zone.UseValueBlanking[GroupIndex] =			TRUE;						//	Boolean_t
		//Zone.UseClipPlanes[GroupIndex] =			TRUE;						//	Boolean_t
		Zone.UseSurfaceTranslucency[GroupIndex] =	FALSE;						//	Boolean_t
		Zone.SurfaceTranslucency[GroupIndex] =		50;		//	value is in %	//	SmInteger_t
		Zone.LightingEffect[GroupIndex] =			LightingEffect_Gouraud;		//	LightingEffect_e

		//	Group index and number (increments for next set)

		GroupIndex++;														//	int


		//===============================================================
		//	RCLines Settings
		//===============================================================

		Zone.ZoneShow[GroupIndex] =					FALSE;						//	Boolean_t

		//	Mesh Settings

		Zone.MeshShow[GroupIndex] =					TRUE;						//	Boolean_t
		Zone.MeshType[GroupIndex] =					MeshType_Overlay;			//	MeshType_e
		Zone.MeshColor[GroupIndex] =					MultiColor8_C;					//	ColorIndex_t
		/*
			A number of color index constants are \#defined. These include Black_C, Red_C,
			Green_C, Blue_C, Cyan_C, Yellow_C, Purple_C, White_C, Custom1_C through
			Custom56_C, MultiColor_C, NoColor_C, MulitiColor2_C through MulitiColor8_C,
			RGBColor_C, and InvalidColor_C.  
		*/
		Zone.MeshLinePttrn[GroupIndex] =				LinePattern_Dotted;			//	LinePattern_e
		Zone.MeshLinePttrnLen[GroupIndex] =			2;		//	Value is in %	//	double
		Zone.MeshLineThickness[GroupIndex] =			0.2;	//	Value is in %	//	double


		//	Contour Settings

		Zone.ContourShow[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.ContourType[GroupIndex] =				ContourType_Flood;			//	ContourType_e
		Zone.FloodColoring[GroupIndex] =				ContourColoring_Group1;		//	ContourColoring_e
		Zone.ContourLineGroup[GroupIndex] =			1;							//	SmInteger_t
		Zone.ContourLineColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.ContourLinePttrn[GroupIndex] =			LinePattern_Solid;			//	LinePattern_e
		Zone.ContourLinePttrnLen[GroupIndex] =		2;		//	Value is in %	//	double
		Zone.ContourLineThickness[GroupIndex] =		0.1;	//	Value is in %	//	double
		Zone.ContourLighting[GroupIndex] =			TRUE;						//	Boolean_t

		//	Vector Settings
		
		Zone.VectorShow[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.VectorType[GroupIndex] =				VectorType_TailAtPoint;		//	VectorType_e
		Zone.VectorArrowType[GroupIndex] =			ArrowheadStyle_Plain;		//	ArrowheadStyle_e
		Zone.VectorLineColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.VectorTangent[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.VectorLinePttrn[GroupIndex] =			LinePattern_Solid;			//	LinePattern_e
		Zone.VectorLinePttrnLen[GroupIndex] =		2;		//	Value is in %	//	double
		Zone.VectorLineThickness[GroupIndex] =		0.1;	//	Value is in %	//	double

		//	Scatter Settings

		Zone.ScatterShow[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.ScatterShape[GroupIndex] =				GeomShape_Sphere;			//	GeomShape_e
		Zone.ScatterOutColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.ScatterInColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.ScatterFillMode[GroupIndex] =			FillMode_UseLineColor;		//	FillMode_e
		Zone.ScatterSize[GroupIndex] =				2.5;	//	value is in %	//	double
		Zone.ScatterLineThickness[GroupIndex] =		0.1;	//	value is in %	//	double

		//	Shade Settings

		Zone.ShadeShow[GroupIndex] =					FALSE;						//	Boolean_t
		Zone.ShadeColor[GroupIndex] =				Black_C;					//	ColorIndex_t
		Zone.ShadeLighting[GroupIndex] =				TRUE;						//	Boolean_t

		//	Edge Settings

		Zone.EdgeShow[GroupIndex] =					FALSE;						//	Boolean_t
		Zone.EdgeType[GroupIndex] =					EdgeType_BordersAndCreases;	//	EdgeType_e
		Zone.EdgeIBorder[GroupIndex] =				BorderLocation_Both;		//	BorderLocation_e
		Zone.EdgeJBorder[GroupIndex] =				BorderLocation_Both;		//	BorderLocation_e
		Zone.EdgeKBorder[GroupIndex] =				BorderLocation_Both;		//	BorderLocation_e
		Zone.EdgeColor[GroupIndex] =					Black_C;					//	ColorIndex_t
		Zone.EdgeLineThickness[GroupIndex] =			0.1;	//	value is in %	//	double

		//	Points Settings

		Zone.PointsToPlot[GroupIndex] =				PointsToPlot_SurfaceNodes;	//	PointsToPlot_e
		Zone.PointsIndexSkipI[GroupIndex] =			1;							//	LgIndex_t
		Zone.PointsIndexSkipJ[GroupIndex] =			1;							//	LgIndex_t
		Zone.PointsIndexSkipK[GroupIndex] =			1;							//	LgIndex_t
		/*
			1 is the same as "No Skip" in the GUI.
		*/

		//	Surfaces Settings

		Zone.SurfacesToPlot[GroupIndex] =			SurfacesToPlot_ExposedCellFaces;	//	SurfacesToPlot_e
		Zone.RangeBeginI[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeEndI[GroupIndex] =					0;							//	LgIndex_t
		Zone.RangeSkipI[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeBeginJ[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeEndJ[GroupIndex] =					0;							//	LgIndex_t
		Zone.RangeSkipJ[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeBeginK[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeEndK[GroupIndex] =					0;							//	LgIndex_t
		Zone.RangeSkipK[GroupIndex] =				1;							//	LgIndex_t
		/*
			0 is the same as maximum, which sets the end range to the maximum
			index for the dimension. The skip number is to skip every # planes,
			so 1 is to not skip any planes, (show every 1 plane).
		*/

		//	Volume Settings

		Zone.StreamtracesShow[GroupIndex] =			TRUE;						//	Boolean_t
		Zone.IsoSurfacesShow[GroupIndex] =			TRUE;						//	Boolean_t
		Zone.SlicesShow[GroupIndex] =				TRUE;						//	Boolean_t

		//	Effects Settings

		Zone.UseValueBlanking[GroupIndex] =			TRUE;						//	Boolean_t
		//Zone.UseClipPlanes[GroupIndex] =			TRUE;						//	Boolean_t
		Zone.UseSurfaceTranslucency[GroupIndex] =	FALSE;						//	Boolean_t
		Zone.SurfaceTranslucency[GroupIndex] =		50;		//	value is in %	//	SmInteger_t
		Zone.LightingEffect[GroupIndex] =			LightingEffect_Gouraud;		//	LightingEffect_e

		//	Group index and number (increments for next set)

		GroupIndex++;														//	int


		//===============================================================
		//	RCLinesFF Settings
		//===============================================================

		Zone.ZoneShow[GroupIndex] =					FALSE;						//	Boolean_t

		//	Mesh Settings

		Zone.MeshShow[GroupIndex] =					TRUE;						//	Boolean_t
		Zone.MeshType[GroupIndex] =					MeshType_Overlay;			//	MeshType_e
		Zone.MeshColor[GroupIndex] =					MultiColor8_C;					//	ColorIndex_t
		/*
			A number of color index constants are \#defined. These include Black_C, Red_C,
			Green_C, Blue_C, Cyan_C, Yellow_C, Purple_C, White_C, Custom1_C through
			Custom56_C, MultiColor_C, NoColor_C, MulitiColor2_C through MulitiColor8_C,
			RGBColor_C, and InvalidColor_C.  
		*/
		Zone.MeshLinePttrn[GroupIndex] =				LinePattern_Dotted;			//	LinePattern_e
		Zone.MeshLinePttrnLen[GroupIndex] =			2;		//	Value is in %	//	double
		Zone.MeshLineThickness[GroupIndex] =			0.2;	//	Value is in %	//	double


		//	Contour Settings

		Zone.ContourShow[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.ContourType[GroupIndex] =				ContourType_Flood;			//	ContourType_e
		Zone.FloodColoring[GroupIndex] =				ContourColoring_Group1;		//	ContourColoring_e
		Zone.ContourLineGroup[GroupIndex] =			1;							//	SmInteger_t
		Zone.ContourLineColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.ContourLinePttrn[GroupIndex] =			LinePattern_Solid;			//	LinePattern_e
		Zone.ContourLinePttrnLen[GroupIndex] =		2;		//	Value is in %	//	double
		Zone.ContourLineThickness[GroupIndex] =		0.1;	//	Value is in %	//	double
		Zone.ContourLighting[GroupIndex] =			TRUE;						//	Boolean_t

		//	Vector Settings
		
		Zone.VectorShow[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.VectorType[GroupIndex] =				VectorType_TailAtPoint;		//	VectorType_e
		Zone.VectorArrowType[GroupIndex] =			ArrowheadStyle_Plain;		//	ArrowheadStyle_e
		Zone.VectorLineColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.VectorTangent[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.VectorLinePttrn[GroupIndex] =			LinePattern_Solid;			//	LinePattern_e
		Zone.VectorLinePttrnLen[GroupIndex] =		2;		//	Value is in %	//	double
		Zone.VectorLineThickness[GroupIndex] =		0.1;	//	Value is in %	//	double

		//	Scatter Settings

		Zone.ScatterShow[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.ScatterShape[GroupIndex] =				GeomShape_Sphere;			//	GeomShape_e
		Zone.ScatterOutColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.ScatterInColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.ScatterFillMode[GroupIndex] =			FillMode_UseLineColor;		//	FillMode_e
		Zone.ScatterSize[GroupIndex] =				2.5;	//	value is in %	//	double
		Zone.ScatterLineThickness[GroupIndex] =		0.1;	//	value is in %	//	double

		//	Shade Settings

		Zone.ShadeShow[GroupIndex] =					FALSE;						//	Boolean_t
		Zone.ShadeColor[GroupIndex] =				Black_C;					//	ColorIndex_t
		Zone.ShadeLighting[GroupIndex] =				TRUE;						//	Boolean_t

		//	Edge Settings

		Zone.EdgeShow[GroupIndex] =					FALSE;						//	Boolean_t
		Zone.EdgeType[GroupIndex] =					EdgeType_BordersAndCreases;	//	EdgeType_e
		Zone.EdgeIBorder[GroupIndex] =				BorderLocation_Both;		//	BorderLocation_e
		Zone.EdgeJBorder[GroupIndex] =				BorderLocation_Both;		//	BorderLocation_e
		Zone.EdgeKBorder[GroupIndex] =				BorderLocation_Both;		//	BorderLocation_e
		Zone.EdgeColor[GroupIndex] =					Black_C;					//	ColorIndex_t
		Zone.EdgeLineThickness[GroupIndex] =			0.1;	//	value is in %	//	double

		//	Points Settings

		Zone.PointsToPlot[GroupIndex] =				PointsToPlot_SurfaceNodes;	//	PointsToPlot_e
		Zone.PointsIndexSkipI[GroupIndex] =			1;							//	LgIndex_t
		Zone.PointsIndexSkipJ[GroupIndex] =			1;							//	LgIndex_t
		Zone.PointsIndexSkipK[GroupIndex] =			1;							//	LgIndex_t
		/*
			1 is the same as "No Skip" in the GUI.
		*/

		//	Surfaces Settings

		Zone.SurfacesToPlot[GroupIndex] =			SurfacesToPlot_ExposedCellFaces;	//	SurfacesToPlot_e
		Zone.RangeBeginI[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeEndI[GroupIndex] =					0;							//	LgIndex_t
		Zone.RangeSkipI[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeBeginJ[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeEndJ[GroupIndex] =					0;							//	LgIndex_t
		Zone.RangeSkipJ[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeBeginK[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeEndK[GroupIndex] =					0;							//	LgIndex_t
		Zone.RangeSkipK[GroupIndex] =				1;							//	LgIndex_t
		/*
			0 is the same as maximum, which sets the end range to the maximum
			index for the dimension. The skip number is to skip every # planes,
			so 1 is to not skip any planes, (show every 1 plane).
		*/

		//	Volume Settings

		Zone.StreamtracesShow[GroupIndex] =			TRUE;						//	Boolean_t
		Zone.IsoSurfacesShow[GroupIndex] =			TRUE;						//	Boolean_t
		Zone.SlicesShow[GroupIndex] =				TRUE;						//	Boolean_t

		//	Effects Settings

		Zone.UseValueBlanking[GroupIndex] =			TRUE;						//	Boolean_t
		//Zone.UseClipPlanes[GroupIndex] =			TRUE;						//	Boolean_t
		Zone.UseSurfaceTranslucency[GroupIndex] =	FALSE;						//	Boolean_t
		Zone.SurfaceTranslucency[GroupIndex] =		50;		//	value is in %	//	SmInteger_t
		Zone.LightingEffect[GroupIndex] =			LightingEffect_Gouraud;		//	LightingEffect_e

		//	Group index and number (increments for next set)

		GroupIndex++;														//	int


		//===============================================================
		//	ABRSurfaces Settings
		//===============================================================

		Zone.ZoneShow[GroupIndex] =					FALSE;						//	Boolean_t

		//	Mesh Settings

		Zone.MeshShow[GroupIndex] =					TRUE;						//	Boolean_t
		Zone.MeshType[GroupIndex] =					MeshType_Overlay;			//	MeshType_e
		Zone.MeshColor[GroupIndex] =					Black_C;					//	ColorIndex_t
		/*
			A number of color index constants are \#defined. These include Black_C, Red_C,
			Green_C, Blue_C, Cyan_C, Yellow_C, Purple_C, White_C, Custom1_C through
			Custom56_C, MultiColor_C, NoColor_C, MulitiColor2_C through MulitiColor8_C,
			RGBColor_C, and InvalidColor_C.  
		*/
		Zone.MeshLinePttrn[GroupIndex] =				LinePattern_Solid;			//	LinePattern_e
		Zone.MeshLinePttrnLen[GroupIndex] =			2;		//	Value is in %	//	double
		Zone.MeshLineThickness[GroupIndex] =			0.1;	//	Value is in %	//	double


		//	Contour Settings

		Zone.ContourShow[GroupIndex] =				TRUE;						//	Boolean_t
		Zone.ContourType[GroupIndex] =				ContourType_Flood;			//	ContourType_e
		Zone.FloodColoring[GroupIndex] =				ContourColoring_Group1;		//	ContourColoring_e
		Zone.ContourLineGroup[GroupIndex] =			1;							//	SmInteger_t
		Zone.ContourLineColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.ContourLinePttrn[GroupIndex] =			LinePattern_Solid;			//	LinePattern_e
		Zone.ContourLinePttrnLen[GroupIndex] =		2;		//	Value is in %	//	double
		Zone.ContourLineThickness[GroupIndex] =		0.1;	//	Value is in %	//	double
		Zone.ContourLighting[GroupIndex] =			TRUE;						//	Boolean_t

		//	Vector Settings
		
		Zone.VectorShow[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.VectorType[GroupIndex] =				VectorType_TailAtPoint;		//	VectorType_e
		Zone.VectorArrowType[GroupIndex] =			ArrowheadStyle_Plain;		//	ArrowheadStyle_e
		Zone.VectorLineColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.VectorTangent[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.VectorLinePttrn[GroupIndex] =			LinePattern_Solid;			//	LinePattern_e
		Zone.VectorLinePttrnLen[GroupIndex] =		2;		//	Value is in %	//	double
		Zone.VectorLineThickness[GroupIndex] =		0.1;	//	Value is in %	//	double

		//	Scatter Settings

		Zone.ScatterShow[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.ScatterShape[GroupIndex] =				GeomShape_Sphere;			//	GeomShape_e
		Zone.ScatterOutColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.ScatterInColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.ScatterFillMode[GroupIndex] =			FillMode_UseLineColor;		//	FillMode_e
		Zone.ScatterSize[GroupIndex] =				2.5;	//	value is in %	//	double
		Zone.ScatterLineThickness[GroupIndex] =		0.1;	//	value is in %	//	double

		//	Shade Settings

		Zone.ShadeShow[GroupIndex] =					TRUE;						//	Boolean_t
		Zone.ShadeColor[GroupIndex] =				Black_C;					//	ColorIndex_t
		Zone.ShadeLighting[GroupIndex] =				TRUE;						//	Boolean_t

		//	Edge Settings

		Zone.EdgeShow[GroupIndex] =					FALSE;						//	Boolean_t
		Zone.EdgeType[GroupIndex] =					EdgeType_BordersAndCreases;	//	EdgeType_e
		Zone.EdgeIBorder[GroupIndex] =				BorderLocation_Both;		//	BorderLocation_e
		Zone.EdgeJBorder[GroupIndex] =				BorderLocation_Both;		//	BorderLocation_e
		Zone.EdgeKBorder[GroupIndex] =				BorderLocation_Both;		//	BorderLocation_e
		Zone.EdgeColor[GroupIndex] =					Black_C;					//	ColorIndex_t
		Zone.EdgeLineThickness[GroupIndex] =			0.1;	//	value is in %	//	double

		//	Points Settings

		Zone.PointsToPlot[GroupIndex] =				PointsToPlot_SurfaceNodes;	//	PointsToPlot_e
		Zone.PointsIndexSkipI[GroupIndex] =			1;							//	LgIndex_t
		Zone.PointsIndexSkipJ[GroupIndex] =			1;							//	LgIndex_t
		Zone.PointsIndexSkipK[GroupIndex] =			1;							//	LgIndex_t
		/*
			1 is the same as "No Skip" in the GUI.
		*/

		//	Surfaces Settings

		Zone.SurfacesToPlot[GroupIndex] =			SurfacesToPlot_ExposedCellFaces;	//	SurfacesToPlot_e
		Zone.RangeBeginI[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeEndI[GroupIndex] =					0;							//	LgIndex_t
		Zone.RangeSkipI[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeBeginJ[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeEndJ[GroupIndex] =					0;							//	LgIndex_t
		Zone.RangeSkipJ[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeBeginK[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeEndK[GroupIndex] =					0;							//	LgIndex_t
		Zone.RangeSkipK[GroupIndex] =				1;							//	LgIndex_t
		/*
			0 is the same as maximum, which sets the end range to the maximum
			index for the dimension. The skip number is to skip every # planes,
			so 1 is to not skip any planes, (show every 1 plane).
		*/

		//	Volume Settings

		Zone.StreamtracesShow[GroupIndex] =			TRUE;						//	Boolean_t
		Zone.IsoSurfacesShow[GroupIndex] =			TRUE;						//	Boolean_t
		Zone.SlicesShow[GroupIndex] =				TRUE;						//	Boolean_t

		//	Effects Settings

		Zone.UseValueBlanking[GroupIndex] =			TRUE;						//	Boolean_t
		//Zone.UseClipPlanes[GroupIndex] =			TRUE;						//	Boolean_t
		Zone.UseSurfaceTranslucency[GroupIndex] =	TRUE;						//	Boolean_t
		Zone.SurfaceTranslucency[GroupIndex] =		50;		//	value is in %	//	SmInteger_t
		Zone.LightingEffect[GroupIndex] =			LightingEffect_Gouraud;		//	LightingEffect_e

		//	Group index and number (increments for next set)

		GroupIndex++;														//	int


		//===============================================================
		//	ABRSurfacesFF Settings
		//===============================================================

		Zone.ZoneShow[GroupIndex] =					FALSE;						//	Boolean_t

		//	Mesh Settings

		Zone.MeshShow[GroupIndex] =					TRUE;						//	Boolean_t
		Zone.MeshType[GroupIndex] =					MeshType_Overlay;			//	MeshType_e
		Zone.MeshColor[GroupIndex] =					Black_C;					//	ColorIndex_t
		/*
			A number of color index constants are \#defined. These include Black_C, Red_C,
			Green_C, Blue_C, Cyan_C, Yellow_C, Purple_C, White_C, Custom1_C through
			Custom56_C, MultiColor_C, NoColor_C, MulitiColor2_C through MulitiColor8_C,
			RGBColor_C, and InvalidColor_C.  
		*/
		Zone.MeshLinePttrn[GroupIndex] =				LinePattern_Solid;			//	LinePattern_e
		Zone.MeshLinePttrnLen[GroupIndex] =			2;		//	Value is in %	//	double
		Zone.MeshLineThickness[GroupIndex] =			0.1;	//	Value is in %	//	double


		//	Contour Settings

		Zone.ContourShow[GroupIndex] =				TRUE;						//	Boolean_t
		Zone.ContourType[GroupIndex] =				ContourType_Flood;			//	ContourType_e
		Zone.FloodColoring[GroupIndex] =				ContourColoring_Group1;		//	ContourColoring_e
		Zone.ContourLineGroup[GroupIndex] =			1;							//	SmInteger_t
		Zone.ContourLineColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.ContourLinePttrn[GroupIndex] =			LinePattern_Solid;			//	LinePattern_e
		Zone.ContourLinePttrnLen[GroupIndex] =		2;		//	Value is in %	//	double
		Zone.ContourLineThickness[GroupIndex] =		0.1;	//	Value is in %	//	double
		Zone.ContourLighting[GroupIndex] =			TRUE;						//	Boolean_t

		//	Vector Settings
		
		Zone.VectorShow[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.VectorType[GroupIndex] =				VectorType_TailAtPoint;		//	VectorType_e
		Zone.VectorArrowType[GroupIndex] =			ArrowheadStyle_Plain;		//	ArrowheadStyle_e
		Zone.VectorLineColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.VectorTangent[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.VectorLinePttrn[GroupIndex] =			LinePattern_Solid;			//	LinePattern_e
		Zone.VectorLinePttrnLen[GroupIndex] =		2;		//	Value is in %	//	double
		Zone.VectorLineThickness[GroupIndex] =		0.1;	//	Value is in %	//	double

		//	Scatter Settings

		Zone.ScatterShow[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.ScatterShape[GroupIndex] =				GeomShape_Sphere;			//	GeomShape_e
		Zone.ScatterOutColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.ScatterInColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.ScatterFillMode[GroupIndex] =			FillMode_UseLineColor;		//	FillMode_e
		Zone.ScatterSize[GroupIndex] =				2.5;	//	value is in %	//	double
		Zone.ScatterLineThickness[GroupIndex] =		0.1;	//	value is in %	//	double

		//	Shade Settings

		Zone.ShadeShow[GroupIndex] =					TRUE;						//	Boolean_t
		Zone.ShadeColor[GroupIndex] =				Black_C;					//	ColorIndex_t
		Zone.ShadeLighting[GroupIndex] =				TRUE;						//	Boolean_t

		//	Edge Settings

		Zone.EdgeShow[GroupIndex] =					FALSE;						//	Boolean_t
		Zone.EdgeType[GroupIndex] =					EdgeType_BordersAndCreases;	//	EdgeType_e
		Zone.EdgeIBorder[GroupIndex] =				BorderLocation_Both;		//	BorderLocation_e
		Zone.EdgeJBorder[GroupIndex] =				BorderLocation_Both;		//	BorderLocation_e
		Zone.EdgeKBorder[GroupIndex] =				BorderLocation_Both;		//	BorderLocation_e
		Zone.EdgeColor[GroupIndex] =					Black_C;					//	ColorIndex_t
		Zone.EdgeLineThickness[GroupIndex] =			0.1;	//	value is in %	//	double

		//	Points Settings

		Zone.PointsToPlot[GroupIndex] =				PointsToPlot_SurfaceNodes;	//	PointsToPlot_e
		Zone.PointsIndexSkipI[GroupIndex] =			1;							//	LgIndex_t
		Zone.PointsIndexSkipJ[GroupIndex] =			1;							//	LgIndex_t
		Zone.PointsIndexSkipK[GroupIndex] =			1;							//	LgIndex_t
		/*
			1 is the same as "No Skip" in the GUI.
		*/

		//	Surfaces Settings

		Zone.SurfacesToPlot[GroupIndex] =			SurfacesToPlot_ExposedCellFaces;	//	SurfacesToPlot_e
		Zone.RangeBeginI[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeEndI[GroupIndex] =					0;							//	LgIndex_t
		Zone.RangeSkipI[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeBeginJ[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeEndJ[GroupIndex] =					0;							//	LgIndex_t
		Zone.RangeSkipJ[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeBeginK[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeEndK[GroupIndex] =					0;							//	LgIndex_t
		Zone.RangeSkipK[GroupIndex] =				1;							//	LgIndex_t
		/*
			0 is the same as maximum, which sets the end range to the maximum
			index for the dimension. The skip number is to skip every # planes,
			so 1 is to not skip any planes, (show every 1 plane).
		*/

		//	Volume Settings

		Zone.StreamtracesShow[GroupIndex] =			TRUE;						//	Boolean_t
		Zone.IsoSurfacesShow[GroupIndex] =			TRUE;						//	Boolean_t
		Zone.SlicesShow[GroupIndex] =				TRUE;						//	Boolean_t

		//	Effects Settings

		Zone.UseValueBlanking[GroupIndex] =			TRUE;						//	Boolean_t
		//Zone.UseClipPlanes[GroupIndex] =			TRUE;						//	Boolean_t
		Zone.UseSurfaceTranslucency[GroupIndex] =	TRUE;						//	Boolean_t
		Zone.SurfaceTranslucency[GroupIndex] =		50;		//	value is in %	//	SmInteger_t
		Zone.LightingEffect[GroupIndex] =			LightingEffect_Gouraud;		//	LightingEffect_e

		//	Group index and number (increments for next set)

		GroupIndex++;														//	int


		//===============================================================
		//	ABCSurfaces Settings
		//===============================================================

		Zone.ZoneShow[GroupIndex] =					FALSE;						//	Boolean_t

		//	Mesh Settings

		Zone.MeshShow[GroupIndex] =					TRUE;						//	Boolean_t
		Zone.MeshType[GroupIndex] =					MeshType_Overlay;			//	MeshType_e
		Zone.MeshColor[GroupIndex] =					Black_C;					//	ColorIndex_t
		/*
			A number of color index constants are \#defined. These include Black_C, Red_C,
			Green_C, Blue_C, Cyan_C, Yellow_C, Purple_C, White_C, Custom1_C through
			Custom56_C, MultiColor_C, NoColor_C, MulitiColor2_C through MulitiColor8_C,
			RGBColor_C, and InvalidColor_C.  
		*/
		Zone.MeshLinePttrn[GroupIndex] =				LinePattern_Solid;			//	LinePattern_e
		Zone.MeshLinePttrnLen[GroupIndex] =			2;		//	Value is in %	//	double
		Zone.MeshLineThickness[GroupIndex] =			0.1;	//	Value is in %	//	double


		//	Contour Settings

		Zone.ContourShow[GroupIndex] =				TRUE;						//	Boolean_t
		Zone.ContourType[GroupIndex] =				ContourType_Flood;			//	ContourType_e
		Zone.FloodColoring[GroupIndex] =				ContourColoring_Group1;		//	ContourColoring_e
		Zone.ContourLineGroup[GroupIndex] =			1;							//	SmInteger_t
		Zone.ContourLineColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.ContourLinePttrn[GroupIndex] =			LinePattern_Solid;			//	LinePattern_e
		Zone.ContourLinePttrnLen[GroupIndex] =		2;		//	Value is in %	//	double
		Zone.ContourLineThickness[GroupIndex] =		0.1;	//	Value is in %	//	double
		Zone.ContourLighting[GroupIndex] =			TRUE;						//	Boolean_t

		//	Vector Settings
		
		Zone.VectorShow[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.VectorType[GroupIndex] =				VectorType_TailAtPoint;		//	VectorType_e
		Zone.VectorArrowType[GroupIndex] =			ArrowheadStyle_Plain;		//	ArrowheadStyle_e
		Zone.VectorLineColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.VectorTangent[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.VectorLinePttrn[GroupIndex] =			LinePattern_Solid;			//	LinePattern_e
		Zone.VectorLinePttrnLen[GroupIndex] =		2;		//	Value is in %	//	double
		Zone.VectorLineThickness[GroupIndex] =		0.1;	//	Value is in %	//	double

		//	Scatter Settings

		Zone.ScatterShow[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.ScatterShape[GroupIndex] =				GeomShape_Sphere;			//	GeomShape_e
		Zone.ScatterOutColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.ScatterInColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.ScatterFillMode[GroupIndex] =			FillMode_UseLineColor;		//	FillMode_e
		Zone.ScatterSize[GroupIndex] =				2.5;	//	value is in %	//	double
		Zone.ScatterLineThickness[GroupIndex] =		0.1;	//	value is in %	//	double

		//	Shade Settings

		Zone.ShadeShow[GroupIndex] =					TRUE;						//	Boolean_t
		Zone.ShadeColor[GroupIndex] =				Black_C;					//	ColorIndex_t
		Zone.ShadeLighting[GroupIndex] =				TRUE;						//	Boolean_t

		//	Edge Settings

		Zone.EdgeShow[GroupIndex] =					FALSE;						//	Boolean_t
		Zone.EdgeType[GroupIndex] =					EdgeType_BordersAndCreases;	//	EdgeType_e
		Zone.EdgeIBorder[GroupIndex] =				BorderLocation_Both;		//	BorderLocation_e
		Zone.EdgeJBorder[GroupIndex] =				BorderLocation_Both;		//	BorderLocation_e
		Zone.EdgeKBorder[GroupIndex] =				BorderLocation_Both;		//	BorderLocation_e
		Zone.EdgeColor[GroupIndex] =					Black_C;					//	ColorIndex_t
		Zone.EdgeLineThickness[GroupIndex] =			0.1;	//	value is in %	//	double

		//	Points Settings

		Zone.PointsToPlot[GroupIndex] =				PointsToPlot_SurfaceNodes;	//	PointsToPlot_e
		Zone.PointsIndexSkipI[GroupIndex] =			1;							//	LgIndex_t
		Zone.PointsIndexSkipJ[GroupIndex] =			1;							//	LgIndex_t
		Zone.PointsIndexSkipK[GroupIndex] =			1;							//	LgIndex_t
		/*
			1 is the same as "No Skip" in the GUI.
		*/

		//	Surfaces Settings

		Zone.SurfacesToPlot[GroupIndex] =			SurfacesToPlot_ExposedCellFaces;	//	SurfacesToPlot_e
		Zone.RangeBeginI[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeEndI[GroupIndex] =					0;							//	LgIndex_t
		Zone.RangeSkipI[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeBeginJ[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeEndJ[GroupIndex] =					0;							//	LgIndex_t
		Zone.RangeSkipJ[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeBeginK[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeEndK[GroupIndex] =					0;							//	LgIndex_t
		Zone.RangeSkipK[GroupIndex] =				1;							//	LgIndex_t
		/*
			0 is the same as maximum, which sets the end range to the maximum
			index for the dimension. The skip number is to skip every # planes,
			so 1 is to not skip any planes, (show every 1 plane).
		*/

		//	Volume Settings

		Zone.StreamtracesShow[GroupIndex] =			TRUE;						//	Boolean_t
		Zone.IsoSurfacesShow[GroupIndex] =			TRUE;						//	Boolean_t
		Zone.SlicesShow[GroupIndex] =				TRUE;						//	Boolean_t

		//	Effects Settings

		Zone.UseValueBlanking[GroupIndex] =			TRUE;						//	Boolean_t
		//Zone.UseClipPlanes[GroupIndex] =			TRUE;						//	Boolean_t
		Zone.UseSurfaceTranslucency[GroupIndex] =	TRUE;						//	Boolean_t
		Zone.SurfaceTranslucency[GroupIndex] =		50;		//	value is in %	//	SmInteger_t
		Zone.LightingEffect[GroupIndex] =			LightingEffect_Gouraud;		//	LightingEffect_e

		//	Group index and number (increments for next set)

		GroupIndex++;														//	int


		//===============================================================
		//	ABCSurfacesFF Settings
		//===============================================================

		Zone.ZoneShow[GroupIndex] =					FALSE;						//	Boolean_t

		//	Mesh Settings

		Zone.MeshShow[GroupIndex] =					TRUE;						//	Boolean_t
		Zone.MeshType[GroupIndex] =					MeshType_Overlay;			//	MeshType_e
		Zone.MeshColor[GroupIndex] =					Black_C;					//	ColorIndex_t
		/*
			A number of color index constants are \#defined. These include Black_C, Red_C,
			Green_C, Blue_C, Cyan_C, Yellow_C, Purple_C, White_C, Custom1_C through
			Custom56_C, MultiColor_C, NoColor_C, MulitiColor2_C through MulitiColor8_C,
			RGBColor_C, and InvalidColor_C.  
		*/
		Zone.MeshLinePttrn[GroupIndex] =				LinePattern_Solid;			//	LinePattern_e
		Zone.MeshLinePttrnLen[GroupIndex] =			2;		//	Value is in %	//	double
		Zone.MeshLineThickness[GroupIndex] =			0.1;	//	Value is in %	//	double


		//	Contour Settings

		Zone.ContourShow[GroupIndex] =				TRUE;						//	Boolean_t
		Zone.ContourType[GroupIndex] =				ContourType_Flood;			//	ContourType_e
		Zone.FloodColoring[GroupIndex] =				ContourColoring_Group1;		//	ContourColoring_e
		Zone.ContourLineGroup[GroupIndex] =			1;							//	SmInteger_t
		Zone.ContourLineColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.ContourLinePttrn[GroupIndex] =			LinePattern_Solid;			//	LinePattern_e
		Zone.ContourLinePttrnLen[GroupIndex] =		2;		//	Value is in %	//	double
		Zone.ContourLineThickness[GroupIndex] =		0.1;	//	Value is in %	//	double
		Zone.ContourLighting[GroupIndex] =			TRUE;						//	Boolean_t

		//	Vector Settings
		
		Zone.VectorShow[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.VectorType[GroupIndex] =				VectorType_TailAtPoint;		//	VectorType_e
		Zone.VectorArrowType[GroupIndex] =			ArrowheadStyle_Plain;		//	ArrowheadStyle_e
		Zone.VectorLineColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.VectorTangent[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.VectorLinePttrn[GroupIndex] =			LinePattern_Solid;			//	LinePattern_e
		Zone.VectorLinePttrnLen[GroupIndex] =		2;		//	Value is in %	//	double
		Zone.VectorLineThickness[GroupIndex] =		0.1;	//	Value is in %	//	double

		//	Scatter Settings

		Zone.ScatterShow[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.ScatterShape[GroupIndex] =				GeomShape_Sphere;			//	GeomShape_e
		Zone.ScatterOutColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.ScatterInColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.ScatterFillMode[GroupIndex] =			FillMode_UseLineColor;		//	FillMode_e
		Zone.ScatterSize[GroupIndex] =				2.5;	//	value is in %	//	double
		Zone.ScatterLineThickness[GroupIndex] =		0.1;	//	value is in %	//	double

		//	Shade Settings

		Zone.ShadeShow[GroupIndex] =					TRUE;						//	Boolean_t
		Zone.ShadeColor[GroupIndex] =				Black_C;					//	ColorIndex_t
		Zone.ShadeLighting[GroupIndex] =				TRUE;						//	Boolean_t

		//	Edge Settings

		Zone.EdgeShow[GroupIndex] =					FALSE;						//	Boolean_t
		Zone.EdgeType[GroupIndex] =					EdgeType_BordersAndCreases;	//	EdgeType_e
		Zone.EdgeIBorder[GroupIndex] =				BorderLocation_Both;		//	BorderLocation_e
		Zone.EdgeJBorder[GroupIndex] =				BorderLocation_Both;		//	BorderLocation_e
		Zone.EdgeKBorder[GroupIndex] =				BorderLocation_Both;		//	BorderLocation_e
		Zone.EdgeColor[GroupIndex] =					Black_C;					//	ColorIndex_t
		Zone.EdgeLineThickness[GroupIndex] =			0.1;	//	value is in %	//	double

		//	Points Settings

		Zone.PointsToPlot[GroupIndex] =				PointsToPlot_SurfaceNodes;	//	PointsToPlot_e
		Zone.PointsIndexSkipI[GroupIndex] =			1;							//	LgIndex_t
		Zone.PointsIndexSkipJ[GroupIndex] =			1;							//	LgIndex_t
		Zone.PointsIndexSkipK[GroupIndex] =			1;							//	LgIndex_t
		/*
			1 is the same as "No Skip" in the GUI.
		*/

		//	Surfaces Settings

		Zone.SurfacesToPlot[GroupIndex] =			SurfacesToPlot_ExposedCellFaces;	//	SurfacesToPlot_e
		Zone.RangeBeginI[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeEndI[GroupIndex] =					0;							//	LgIndex_t
		Zone.RangeSkipI[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeBeginJ[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeEndJ[GroupIndex] =					0;							//	LgIndex_t
		Zone.RangeSkipJ[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeBeginK[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeEndK[GroupIndex] =					0;							//	LgIndex_t
		Zone.RangeSkipK[GroupIndex] =				1;							//	LgIndex_t
		/*
			0 is the same as maximum, which sets the end range to the maximum
			index for the dimension. The skip number is to skip every # planes,
			so 1 is to not skip any planes, (show every 1 plane).
		*/

		//	Volume Settings

		Zone.StreamtracesShow[GroupIndex] =			TRUE;						//	Boolean_t
		Zone.IsoSurfacesShow[GroupIndex] =			TRUE;						//	Boolean_t
		Zone.SlicesShow[GroupIndex] =				TRUE;						//	Boolean_t

		//	Effects Settings

		Zone.UseValueBlanking[GroupIndex] =			TRUE;						//	Boolean_t
		//Zone.UseClipPlanes[GroupIndex] =			TRUE;						//	Boolean_t
		Zone.UseSurfaceTranslucency[GroupIndex] =	TRUE;						//	Boolean_t
		Zone.SurfaceTranslucency[GroupIndex] =		50;		//	value is in %	//	SmInteger_t
		Zone.LightingEffect[GroupIndex] =			LightingEffect_Gouraud;		//	LightingEffect_e

		//	Group index and number (increments for next set)

		GroupIndex++;														//	int


		//===============================================================
		//	ARCSurfaces Settings
		//===============================================================

		Zone.ZoneShow[GroupIndex] =					TRUE;						//	Boolean_t

		//	Mesh Settings

		Zone.MeshShow[GroupIndex] =					TRUE;						//	Boolean_t
		Zone.MeshType[GroupIndex] =					MeshType_Overlay;			//	MeshType_e
		Zone.MeshColor[GroupIndex] =					Black_C;					//	ColorIndex_t
		/*
			A number of color index constants are \#defined. These include Black_C, Red_C,
			Green_C, Blue_C, Cyan_C, Yellow_C, Purple_C, White_C, Custom1_C through
			Custom56_C, MultiColor_C, NoColor_C, MulitiColor2_C through MulitiColor8_C,
			RGBColor_C, and InvalidColor_C.  
		*/
		Zone.MeshLinePttrn[GroupIndex] =				LinePattern_Solid;			//	LinePattern_e
		Zone.MeshLinePttrnLen[GroupIndex] =			2;		//	Value is in %	//	double
		Zone.MeshLineThickness[GroupIndex] =			0.1;	//	Value is in %	//	double


		//	Contour Settings

		Zone.ContourShow[GroupIndex] =				TRUE;						//	Boolean_t
		Zone.ContourType[GroupIndex] =				ContourType_Flood;			//	ContourType_e
		Zone.FloodColoring[GroupIndex] =				ContourColoring_Group1;		//	ContourColoring_e
		Zone.ContourLineGroup[GroupIndex] =			1;							//	SmInteger_t
		Zone.ContourLineColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.ContourLinePttrn[GroupIndex] =			LinePattern_Solid;			//	LinePattern_e
		Zone.ContourLinePttrnLen[GroupIndex] =		2;		//	Value is in %	//	double
		Zone.ContourLineThickness[GroupIndex] =		0.1;	//	Value is in %	//	double
		Zone.ContourLighting[GroupIndex] =			TRUE;						//	Boolean_t

		//	Vector Settings
		
		Zone.VectorShow[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.VectorType[GroupIndex] =				VectorType_TailAtPoint;		//	VectorType_e
		Zone.VectorArrowType[GroupIndex] =			ArrowheadStyle_Plain;		//	ArrowheadStyle_e
		Zone.VectorLineColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.VectorTangent[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.VectorLinePttrn[GroupIndex] =			LinePattern_Solid;			//	LinePattern_e
		Zone.VectorLinePttrnLen[GroupIndex] =		2;		//	Value is in %	//	double
		Zone.VectorLineThickness[GroupIndex] =		0.1;	//	Value is in %	//	double

		//	Scatter Settings

		Zone.ScatterShow[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.ScatterShape[GroupIndex] =				GeomShape_Sphere;			//	GeomShape_e
		Zone.ScatterOutColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.ScatterInColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.ScatterFillMode[GroupIndex] =			FillMode_UseLineColor;		//	FillMode_e
		Zone.ScatterSize[GroupIndex] =				2.5;	//	value is in %	//	double
		Zone.ScatterLineThickness[GroupIndex] =		0.1;	//	value is in %	//	double

		//	Shade Settings

		Zone.ShadeShow[GroupIndex] =					TRUE;						//	Boolean_t
		Zone.ShadeColor[GroupIndex] =				Black_C;					//	ColorIndex_t
		Zone.ShadeLighting[GroupIndex] =				TRUE;						//	Boolean_t

		//	Edge Settings

		Zone.EdgeShow[GroupIndex] =					FALSE;						//	Boolean_t
		Zone.EdgeType[GroupIndex] =					EdgeType_BordersAndCreases;	//	EdgeType_e
		Zone.EdgeIBorder[GroupIndex] =				BorderLocation_Both;		//	BorderLocation_e
		Zone.EdgeJBorder[GroupIndex] =				BorderLocation_Both;		//	BorderLocation_e
		Zone.EdgeKBorder[GroupIndex] =				BorderLocation_Both;		//	BorderLocation_e
		Zone.EdgeColor[GroupIndex] =					Black_C;					//	ColorIndex_t
		Zone.EdgeLineThickness[GroupIndex] =			0.1;	//	value is in %	//	double

		//	Points Settings

		Zone.PointsToPlot[GroupIndex] =				PointsToPlot_SurfaceNodes;	//	PointsToPlot_e
		Zone.PointsIndexSkipI[GroupIndex] =			1;							//	LgIndex_t
		Zone.PointsIndexSkipJ[GroupIndex] =			1;							//	LgIndex_t
		Zone.PointsIndexSkipK[GroupIndex] =			1;							//	LgIndex_t
		/*
			1 is the same as "No Skip" in the GUI.
		*/

		//	Surfaces Settings

		Zone.SurfacesToPlot[GroupIndex] =			SurfacesToPlot_ExposedCellFaces;	//	SurfacesToPlot_e
		Zone.RangeBeginI[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeEndI[GroupIndex] =					0;							//	LgIndex_t
		Zone.RangeSkipI[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeBeginJ[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeEndJ[GroupIndex] =					0;							//	LgIndex_t
		Zone.RangeSkipJ[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeBeginK[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeEndK[GroupIndex] =					0;							//	LgIndex_t
		Zone.RangeSkipK[GroupIndex] =				1;							//	LgIndex_t
		/*
			0 is the same as maximum, which sets the end range to the maximum
			index for the dimension. The skip number is to skip every # planes,
			so 1 is to not skip any planes, (show every 1 plane).
		*/

		//	Volume Settings

		Zone.StreamtracesShow[GroupIndex] =			TRUE;						//	Boolean_t
		Zone.IsoSurfacesShow[GroupIndex] =			TRUE;						//	Boolean_t
		Zone.SlicesShow[GroupIndex] =				TRUE;						//	Boolean_t

		//	Effects Settings

		Zone.UseValueBlanking[GroupIndex] =			TRUE;						//	Boolean_t
		//Zone.UseClipPlanes[GroupIndex] =			TRUE;						//	Boolean_t
		Zone.UseSurfaceTranslucency[GroupIndex] =	TRUE;						//	Boolean_t
		Zone.SurfaceTranslucency[GroupIndex] =		50;		//	value is in %	//	SmInteger_t
		Zone.LightingEffect[GroupIndex] =			LightingEffect_Gouraud;		//	LightingEffect_e

		//	Group index and number (increments for next set)

		GroupIndex++;														//	int


		//===============================================================
		//	ARCSurfacesFF Settings
		//===============================================================

		Zone.ZoneShow[GroupIndex] =					TRUE;						//	Boolean_t

		//	Mesh Settings

		Zone.MeshShow[GroupIndex] =					TRUE;						//	Boolean_t
		Zone.MeshType[GroupIndex] =					MeshType_Overlay;			//	MeshType_e
		Zone.MeshColor[GroupIndex] =					Black_C;					//	ColorIndex_t
		/*
			A number of color index constants are \#defined. These include Black_C, Red_C,
			Green_C, Blue_C, Cyan_C, Yellow_C, Purple_C, White_C, Custom1_C through
			Custom56_C, MultiColor_C, NoColor_C, MulitiColor2_C through MulitiColor8_C,
			RGBColor_C, and InvalidColor_C.  
		*/
		Zone.MeshLinePttrn[GroupIndex] =				LinePattern_Solid;			//	LinePattern_e
		Zone.MeshLinePttrnLen[GroupIndex] =			2;		//	Value is in %	//	double
		Zone.MeshLineThickness[GroupIndex] =			0.1;	//	Value is in %	//	double


		//	Contour Settings

		Zone.ContourShow[GroupIndex] =				TRUE;						//	Boolean_t
		Zone.ContourType[GroupIndex] =				ContourType_Flood;			//	ContourType_e
		Zone.FloodColoring[GroupIndex] =				ContourColoring_Group1;		//	ContourColoring_e
		Zone.ContourLineGroup[GroupIndex] =			1;							//	SmInteger_t
		Zone.ContourLineColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.ContourLinePttrn[GroupIndex] =			LinePattern_Solid;			//	LinePattern_e
		Zone.ContourLinePttrnLen[GroupIndex] =		2;		//	Value is in %	//	double
		Zone.ContourLineThickness[GroupIndex] =		0.1;	//	Value is in %	//	double
		Zone.ContourLighting[GroupIndex] =			TRUE;						//	Boolean_t

		//	Vector Settings
		
		Zone.VectorShow[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.VectorType[GroupIndex] =				VectorType_TailAtPoint;		//	VectorType_e
		Zone.VectorArrowType[GroupIndex] =			ArrowheadStyle_Plain;		//	ArrowheadStyle_e
		Zone.VectorLineColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.VectorTangent[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.VectorLinePttrn[GroupIndex] =			LinePattern_Solid;			//	LinePattern_e
		Zone.VectorLinePttrnLen[GroupIndex] =		2;		//	Value is in %	//	double
		Zone.VectorLineThickness[GroupIndex] =		0.1;	//	Value is in %	//	double

		//	Scatter Settings

		Zone.ScatterShow[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.ScatterShape[GroupIndex] =				GeomShape_Sphere;			//	GeomShape_e
		Zone.ScatterOutColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.ScatterInColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.ScatterFillMode[GroupIndex] =			FillMode_UseLineColor;		//	FillMode_e
		Zone.ScatterSize[GroupIndex] =				2.5;	//	value is in %	//	double
		Zone.ScatterLineThickness[GroupIndex] =		0.1;	//	value is in %	//	double

		//	Shade Settings

		Zone.ShadeShow[GroupIndex] =					TRUE;						//	Boolean_t
		Zone.ShadeColor[GroupIndex] =				Black_C;					//	ColorIndex_t
		Zone.ShadeLighting[GroupIndex] =				TRUE;						//	Boolean_t

		//	Edge Settings

		Zone.EdgeShow[GroupIndex] =					FALSE;						//	Boolean_t
		Zone.EdgeType[GroupIndex] =					EdgeType_BordersAndCreases;	//	EdgeType_e
		Zone.EdgeIBorder[GroupIndex] =				BorderLocation_Both;		//	BorderLocation_e
		Zone.EdgeJBorder[GroupIndex] =				BorderLocation_Both;		//	BorderLocation_e
		Zone.EdgeKBorder[GroupIndex] =				BorderLocation_Both;		//	BorderLocation_e
		Zone.EdgeColor[GroupIndex] =					Black_C;					//	ColorIndex_t
		Zone.EdgeLineThickness[GroupIndex] =			0.1;	//	value is in %	//	double

		//	Points Settings

		Zone.PointsToPlot[GroupIndex] =				PointsToPlot_SurfaceNodes;	//	PointsToPlot_e
		Zone.PointsIndexSkipI[GroupIndex] =			1;							//	LgIndex_t
		Zone.PointsIndexSkipJ[GroupIndex] =			1;							//	LgIndex_t
		Zone.PointsIndexSkipK[GroupIndex] =			1;							//	LgIndex_t
		/*
			1 is the same as "No Skip" in the GUI.
		*/

		//	Surfaces Settings

		Zone.SurfacesToPlot[GroupIndex] =			SurfacesToPlot_ExposedCellFaces;	//	SurfacesToPlot_e
		Zone.RangeBeginI[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeEndI[GroupIndex] =					0;							//	LgIndex_t
		Zone.RangeSkipI[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeBeginJ[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeEndJ[GroupIndex] =					0;							//	LgIndex_t
		Zone.RangeSkipJ[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeBeginK[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeEndK[GroupIndex] =					0;							//	LgIndex_t
		Zone.RangeSkipK[GroupIndex] =				1;							//	LgIndex_t
		/*
			0 is the same as maximum, which sets the end range to the maximum
			index for the dimension. The skip number is to skip every # planes,
			so 1 is to not skip any planes, (show every 1 plane).
		*/

		//	Volume Settings

		Zone.StreamtracesShow[GroupIndex] =			TRUE;						//	Boolean_t
		Zone.IsoSurfacesShow[GroupIndex] =			TRUE;						//	Boolean_t
		Zone.SlicesShow[GroupIndex] =				TRUE;						//	Boolean_t

		//	Effects Settings

		Zone.UseValueBlanking[GroupIndex] =			TRUE;						//	Boolean_t
		//Zone.UseClipPlanes[GroupIndex] =			TRUE;						//	Boolean_t
		Zone.UseSurfaceTranslucency[GroupIndex] =	TRUE;						//	Boolean_t
		Zone.SurfaceTranslucency[GroupIndex] =		50;		//	value is in %	//	SmInteger_t
		Zone.LightingEffect[GroupIndex] =			LightingEffect_Gouraud;		//	LightingEffect_e

		//	Group index and number (increments for next set)

		GroupIndex++;														//	int


		//===============================================================
		//	BRCSurfaces Settings
		//===============================================================

		Zone.ZoneShow[GroupIndex] =					FALSE;						//	Boolean_t

		//	Mesh Settings

		Zone.MeshShow[GroupIndex] =					TRUE;						//	Boolean_t
		Zone.MeshType[GroupIndex] =					MeshType_Overlay;			//	MeshType_e
		Zone.MeshColor[GroupIndex] =					Black_C;					//	ColorIndex_t
		/*
			A number of color index constants are \#defined. These include Black_C, Red_C,
			Green_C, Blue_C, Cyan_C, Yellow_C, Purple_C, White_C, Custom1_C through
			Custom56_C, MultiColor_C, NoColor_C, MulitiColor2_C through MulitiColor8_C,
			RGBColor_C, and InvalidColor_C.  
		*/
		Zone.MeshLinePttrn[GroupIndex] =				LinePattern_Solid;			//	LinePattern_e
		Zone.MeshLinePttrnLen[GroupIndex] =			2;		//	Value is in %	//	double
		Zone.MeshLineThickness[GroupIndex] =			0.1;	//	Value is in %	//	double


		//	Contour Settings

		Zone.ContourShow[GroupIndex] =				TRUE;						//	Boolean_t
		Zone.ContourType[GroupIndex] =				ContourType_Flood;			//	ContourType_e
		Zone.FloodColoring[GroupIndex] =				ContourColoring_Group1;		//	ContourColoring_e
		Zone.ContourLineGroup[GroupIndex] =			1;							//	SmInteger_t
		Zone.ContourLineColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.ContourLinePttrn[GroupIndex] =			LinePattern_Solid;			//	LinePattern_e
		Zone.ContourLinePttrnLen[GroupIndex] =		2;		//	Value is in %	//	double
		Zone.ContourLineThickness[GroupIndex] =		0.1;	//	Value is in %	//	double
		Zone.ContourLighting[GroupIndex] =			TRUE;						//	Boolean_t

		//	Vector Settings
		
		Zone.VectorShow[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.VectorType[GroupIndex] =				VectorType_TailAtPoint;		//	VectorType_e
		Zone.VectorArrowType[GroupIndex] =			ArrowheadStyle_Plain;		//	ArrowheadStyle_e
		Zone.VectorLineColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.VectorTangent[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.VectorLinePttrn[GroupIndex] =			LinePattern_Solid;			//	LinePattern_e
		Zone.VectorLinePttrnLen[GroupIndex] =		2;		//	Value is in %	//	double
		Zone.VectorLineThickness[GroupIndex] =		0.1;	//	Value is in %	//	double

		//	Scatter Settings

		Zone.ScatterShow[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.ScatterShape[GroupIndex] =				GeomShape_Sphere;			//	GeomShape_e
		Zone.ScatterOutColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.ScatterInColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.ScatterFillMode[GroupIndex] =			FillMode_UseLineColor;		//	FillMode_e
		Zone.ScatterSize[GroupIndex] =				2.5;	//	value is in %	//	double
		Zone.ScatterLineThickness[GroupIndex] =		0.1;	//	value is in %	//	double

		//	Shade Settings

		Zone.ShadeShow[GroupIndex] =					TRUE;						//	Boolean_t
		Zone.ShadeColor[GroupIndex] =				Black_C;					//	ColorIndex_t
		Zone.ShadeLighting[GroupIndex] =				TRUE;						//	Boolean_t

		//	Edge Settings

		Zone.EdgeShow[GroupIndex] =					FALSE;						//	Boolean_t
		Zone.EdgeType[GroupIndex] =					EdgeType_BordersAndCreases;	//	EdgeType_e
		Zone.EdgeIBorder[GroupIndex] =				BorderLocation_Both;		//	BorderLocation_e
		Zone.EdgeJBorder[GroupIndex] =				BorderLocation_Both;		//	BorderLocation_e
		Zone.EdgeKBorder[GroupIndex] =				BorderLocation_Both;		//	BorderLocation_e
		Zone.EdgeColor[GroupIndex] =					Black_C;					//	ColorIndex_t
		Zone.EdgeLineThickness[GroupIndex] =			0.1;	//	value is in %	//	double

		//	Points Settings

		Zone.PointsToPlot[GroupIndex] =				PointsToPlot_SurfaceNodes;	//	PointsToPlot_e
		Zone.PointsIndexSkipI[GroupIndex] =			1;							//	LgIndex_t
		Zone.PointsIndexSkipJ[GroupIndex] =			1;							//	LgIndex_t
		Zone.PointsIndexSkipK[GroupIndex] =			1;							//	LgIndex_t
		/*
			1 is the same as "No Skip" in the GUI.
		*/

		//	Surfaces Settings

		Zone.SurfacesToPlot[GroupIndex] =			SurfacesToPlot_ExposedCellFaces;	//	SurfacesToPlot_e
		Zone.RangeBeginI[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeEndI[GroupIndex] =					0;							//	LgIndex_t
		Zone.RangeSkipI[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeBeginJ[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeEndJ[GroupIndex] =					0;							//	LgIndex_t
		Zone.RangeSkipJ[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeBeginK[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeEndK[GroupIndex] =					0;							//	LgIndex_t
		Zone.RangeSkipK[GroupIndex] =				1;							//	LgIndex_t
		/*
			0 is the same as maximum, which sets the end range to the maximum
			index for the dimension. The skip number is to skip every # planes,
			so 1 is to not skip any planes, (show every 1 plane).
		*/

		//	Volume Settings

		Zone.StreamtracesShow[GroupIndex] =			TRUE;						//	Boolean_t
		Zone.IsoSurfacesShow[GroupIndex] =			TRUE;						//	Boolean_t
		Zone.SlicesShow[GroupIndex] =				TRUE;						//	Boolean_t

		//	Effects Settings

		Zone.UseValueBlanking[GroupIndex] =			TRUE;						//	Boolean_t
		//Zone.UseClipPlanes[GroupIndex] =			TRUE;						//	Boolean_t
		Zone.UseSurfaceTranslucency[GroupIndex] =	TRUE;						//	Boolean_t
		Zone.SurfaceTranslucency[GroupIndex] =		50;		//	value is in %	//	SmInteger_t
		Zone.LightingEffect[GroupIndex] =			LightingEffect_Gouraud;		//	LightingEffect_e

		//	Group index and number (increments for next set)

		GroupIndex++;														//	int


		//===============================================================
		//	BRCSurfacesFF Settings
		//===============================================================

		Zone.ZoneShow[GroupIndex] =					FALSE;						//	Boolean_t

		//	Mesh Settings

		Zone.MeshShow[GroupIndex] =					TRUE;						//	Boolean_t
		Zone.MeshType[GroupIndex] =					MeshType_Overlay;			//	MeshType_e
		Zone.MeshColor[GroupIndex] =					Black_C;					//	ColorIndex_t
		/*
			A number of color index constants are \#defined. These include Black_C, Red_C,
			Green_C, Blue_C, Cyan_C, Yellow_C, Purple_C, White_C, Custom1_C through
			Custom56_C, MultiColor_C, NoColor_C, MulitiColor2_C through MulitiColor8_C,
			RGBColor_C, and InvalidColor_C.  
		*/
		Zone.MeshLinePttrn[GroupIndex] =				LinePattern_Solid;			//	LinePattern_e
		Zone.MeshLinePttrnLen[GroupIndex] =			2;		//	Value is in %	//	double
		Zone.MeshLineThickness[GroupIndex] =			0.1;	//	Value is in %	//	double


		//	Contour Settings

		Zone.ContourShow[GroupIndex] =				TRUE;						//	Boolean_t
		Zone.ContourType[GroupIndex] =				ContourType_Flood;			//	ContourType_e
		Zone.FloodColoring[GroupIndex] =				ContourColoring_Group1;		//	ContourColoring_e
		Zone.ContourLineGroup[GroupIndex] =			1;							//	SmInteger_t
		Zone.ContourLineColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.ContourLinePttrn[GroupIndex] =			LinePattern_Solid;			//	LinePattern_e
		Zone.ContourLinePttrnLen[GroupIndex] =		2;		//	Value is in %	//	double
		Zone.ContourLineThickness[GroupIndex] =		0.1;	//	Value is in %	//	double
		Zone.ContourLighting[GroupIndex] =			TRUE;						//	Boolean_t

		//	Vector Settings
		
		Zone.VectorShow[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.VectorType[GroupIndex] =				VectorType_TailAtPoint;		//	VectorType_e
		Zone.VectorArrowType[GroupIndex] =			ArrowheadStyle_Plain;		//	ArrowheadStyle_e
		Zone.VectorLineColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.VectorTangent[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.VectorLinePttrn[GroupIndex] =			LinePattern_Solid;			//	LinePattern_e
		Zone.VectorLinePttrnLen[GroupIndex] =		2;		//	Value is in %	//	double
		Zone.VectorLineThickness[GroupIndex] =		0.1;	//	Value is in %	//	double

		//	Scatter Settings

		Zone.ScatterShow[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.ScatterShape[GroupIndex] =				GeomShape_Sphere;			//	GeomShape_e
		Zone.ScatterOutColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.ScatterInColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.ScatterFillMode[GroupIndex] =			FillMode_UseLineColor;		//	FillMode_e
		Zone.ScatterSize[GroupIndex] =				2.5;	//	value is in %	//	double
		Zone.ScatterLineThickness[GroupIndex] =		0.1;	//	value is in %	//	double

		//	Shade Settings

		Zone.ShadeShow[GroupIndex] =					TRUE;						//	Boolean_t
		Zone.ShadeColor[GroupIndex] =				Black_C;					//	ColorIndex_t
		Zone.ShadeLighting[GroupIndex] =				TRUE;						//	Boolean_t

		//	Edge Settings

		Zone.EdgeShow[GroupIndex] =					FALSE;						//	Boolean_t
		Zone.EdgeType[GroupIndex] =					EdgeType_BordersAndCreases;	//	EdgeType_e
		Zone.EdgeIBorder[GroupIndex] =				BorderLocation_Both;		//	BorderLocation_e
		Zone.EdgeJBorder[GroupIndex] =				BorderLocation_Both;		//	BorderLocation_e
		Zone.EdgeKBorder[GroupIndex] =				BorderLocation_Both;		//	BorderLocation_e
		Zone.EdgeColor[GroupIndex] =					Black_C;					//	ColorIndex_t
		Zone.EdgeLineThickness[GroupIndex] =			0.1;	//	value is in %	//	double

		//	Points Settings

		Zone.PointsToPlot[GroupIndex] =				PointsToPlot_SurfaceNodes;	//	PointsToPlot_e
		Zone.PointsIndexSkipI[GroupIndex] =			1;							//	LgIndex_t
		Zone.PointsIndexSkipJ[GroupIndex] =			1;							//	LgIndex_t
		Zone.PointsIndexSkipK[GroupIndex] =			1;							//	LgIndex_t
		/*
			1 is the same as "No Skip" in the GUI.
		*/

		//	Surfaces Settings

		Zone.SurfacesToPlot[GroupIndex] =			SurfacesToPlot_ExposedCellFaces;	//	SurfacesToPlot_e
		Zone.RangeBeginI[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeEndI[GroupIndex] =					0;							//	LgIndex_t
		Zone.RangeSkipI[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeBeginJ[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeEndJ[GroupIndex] =					0;							//	LgIndex_t
		Zone.RangeSkipJ[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeBeginK[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeEndK[GroupIndex] =					0;							//	LgIndex_t
		Zone.RangeSkipK[GroupIndex] =				1;							//	LgIndex_t
		/*
			0 is the same as maximum, which sets the end range to the maximum
			index for the dimension. The skip number is to skip every # planes,
			so 1 is to not skip any planes, (show every 1 plane).
		*/

		//	Volume Settings

		Zone.StreamtracesShow[GroupIndex] =			TRUE;						//	Boolean_t
		Zone.IsoSurfacesShow[GroupIndex] =			TRUE;						//	Boolean_t
		Zone.SlicesShow[GroupIndex] =				TRUE;						//	Boolean_t

		//	Effects Settings

		Zone.UseValueBlanking[GroupIndex] =			TRUE;						//	Boolean_t
		//Zone.UseClipPlanes[GroupIndex] =			TRUE;						//	Boolean_t
		Zone.UseSurfaceTranslucency[GroupIndex] =	TRUE;						//	Boolean_t
		Zone.SurfaceTranslucency[GroupIndex] =		50;		//	value is in %	//	SmInteger_t
		Zone.LightingEffect[GroupIndex] =			LightingEffect_Gouraud;		//	LightingEffect_e

		//	Group index and number (increments for next set)

		GroupIndex++;														//	int


		//===============================================================
		//	AtomRhoIsoSurfaces Settings
		//===============================================================

		Zone.ZoneShow[GroupIndex] =					FALSE;						//	Boolean_t

		//	Mesh Settings

		Zone.MeshShow[GroupIndex] =					FALSE;						//	Boolean_t
		Zone.MeshType[GroupIndex] =					MeshType_Overlay;			//	MeshType_e
		Zone.MeshColor[GroupIndex] =					MultiColor3_C;					//	ColorIndex_t
		/*
			A number of color index constants are \#defined. These include Black_C, Red_C,
			Green_C, Blue_C, Cyan_C, Yellow_C, Purple_C, White_C, Custom1_C through
			Custom56_C, MultiColor_C, NoColor_C, MulitiColor2_C through MulitiColor8_C,
			RGBColor_C, and InvalidColor_C.  
		*/
		Zone.MeshLinePttrn[GroupIndex] =				LinePattern_Solid;			//	LinePattern_e
		Zone.MeshLinePttrnLen[GroupIndex] =			2;		//	Value is in %	//	double
		Zone.MeshLineThickness[GroupIndex] =			0.1;	//	Value is in %	//	double


		//	Contour Settings

		Zone.ContourShow[GroupIndex] =				TRUE;						//	Boolean_t
		Zone.ContourType[GroupIndex] =				ContourType_Flood;			//	ContourType_e
		Zone.FloodColoring[GroupIndex] =				ContourColoring_Group3;		//	ContourColoring_e
		Zone.ContourLineGroup[GroupIndex] =			3;							//	SmInteger_t
		Zone.ContourLineColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.ContourLinePttrn[GroupIndex] =			LinePattern_Solid;			//	LinePattern_e
		Zone.ContourLinePttrnLen[GroupIndex] =		2;		//	Value is in %	//	double
		Zone.ContourLineThickness[GroupIndex] =		0.1;	//	Value is in %	//	double
		Zone.ContourLighting[GroupIndex] =			TRUE;						//	Boolean_t

		//	Vector Settings
		
		Zone.VectorShow[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.VectorType[GroupIndex] =				VectorType_TailAtPoint;		//	VectorType_e
		Zone.VectorArrowType[GroupIndex] =			ArrowheadStyle_Plain;		//	ArrowheadStyle_e
		Zone.VectorLineColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.VectorTangent[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.VectorLinePttrn[GroupIndex] =			LinePattern_Solid;			//	LinePattern_e
		Zone.VectorLinePttrnLen[GroupIndex] =		2;		//	Value is in %	//	double
		Zone.VectorLineThickness[GroupIndex] =		0.1;	//	Value is in %	//	double

		//	Scatter Settings

		Zone.ScatterShow[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.ScatterShape[GroupIndex] =				GeomShape_Sphere;			//	GeomShape_e
		Zone.ScatterOutColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.ScatterInColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.ScatterFillMode[GroupIndex] =			FillMode_UseLineColor;		//	FillMode_e
		Zone.ScatterSize[GroupIndex] =				2.5;	//	value is in %	//	double
		Zone.ScatterLineThickness[GroupIndex] =		0.1;	//	value is in %	//	double

		//	Shade Settings

		Zone.ShadeShow[GroupIndex] =					FALSE;						//	Boolean_t
		Zone.ShadeColor[GroupIndex] =				Black_C;					//	ColorIndex_t
		Zone.ShadeLighting[GroupIndex] =				TRUE;						//	Boolean_t

		//	Edge Settings

		Zone.EdgeShow[GroupIndex] =					FALSE;						//	Boolean_t
		Zone.EdgeType[GroupIndex] =					EdgeType_BordersAndCreases;	//	EdgeType_e
		Zone.EdgeIBorder[GroupIndex] =				BorderLocation_Both;		//	BorderLocation_e
		Zone.EdgeJBorder[GroupIndex] =				BorderLocation_Both;		//	BorderLocation_e
		Zone.EdgeKBorder[GroupIndex] =				BorderLocation_Both;		//	BorderLocation_e
		Zone.EdgeColor[GroupIndex] =					Black_C;					//	ColorIndex_t
		Zone.EdgeLineThickness[GroupIndex] =			0.1;	//	value is in %	//	double

		//	Points Settings

		Zone.PointsToPlot[GroupIndex] =				PointsToPlot_SurfaceNodes;	//	PointsToPlot_e
		Zone.PointsIndexSkipI[GroupIndex] =			1;							//	LgIndex_t
		Zone.PointsIndexSkipJ[GroupIndex] =			1;							//	LgIndex_t
		Zone.PointsIndexSkipK[GroupIndex] =			1;							//	LgIndex_t
		/*
			1 is the same as "No Skip" in the GUI.
		*/

		//	Surfaces Settings

		Zone.SurfacesToPlot[GroupIndex] =			SurfacesToPlot_ExposedCellFaces;	//	SurfacesToPlot_e
		Zone.RangeBeginI[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeEndI[GroupIndex] =					0;							//	LgIndex_t
		Zone.RangeSkipI[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeBeginJ[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeEndJ[GroupIndex] =					0;							//	LgIndex_t
		Zone.RangeSkipJ[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeBeginK[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeEndK[GroupIndex] =					0;							//	LgIndex_t
		Zone.RangeSkipK[GroupIndex] =				1;							//	LgIndex_t
		/*
			0 is the same as maximum, which sets the end range to the maximum
			index for the dimension. The skip number is to skip every # planes,
			so 1 is to not skip any planes, (show every 1 plane).
		*/

		//	Volume Settings

		Zone.StreamtracesShow[GroupIndex] =			TRUE;						//	Boolean_t
		Zone.IsoSurfacesShow[GroupIndex] =			TRUE;						//	Boolean_t
		Zone.SlicesShow[GroupIndex] =				TRUE;						//	Boolean_t

		//	Effects Settings

		Zone.UseValueBlanking[GroupIndex] =			FALSE;						//	Boolean_t
		//Zone.UseClipPlanes[GroupIndex] =			TRUE;						//	Boolean_t
		Zone.UseSurfaceTranslucency[GroupIndex] =	TRUE;						//	Boolean_t
		Zone.SurfaceTranslucency[GroupIndex] =		50;		//	value is in %	//	SmInteger_t
		Zone.LightingEffect[GroupIndex] =			LightingEffect_Gouraud;		//	LightingEffect_e

		//	Group index and number (increments for next set)

		GroupIndex++;														//	int


		//===============================================================
		//	ExtraIsoSurfaces Settings
		//===============================================================

		Zone.ZoneShow[GroupIndex] =					FALSE;						//	Boolean_t

		//	Mesh Settings

		Zone.MeshShow[GroupIndex] =					TRUE;						//	Boolean_t
		Zone.MeshType[GroupIndex] =					MeshType_Overlay;			//	MeshType_e
		Zone.MeshColor[GroupIndex] =					Black_C;					//	ColorIndex_t
		/*
			A number of color index constants are \#defined. These include Black_C, Red_C,
			Green_C, Blue_C, Cyan_C, Yellow_C, Purple_C, White_C, Custom1_C through
			Custom56_C, MultiColor_C, NoColor_C, MulitiColor2_C through MulitiColor8_C,
			RGBColor_C, and InvalidColor_C.  
		*/
		Zone.MeshLinePttrn[GroupIndex] =				LinePattern_Solid;			//	LinePattern_e
		Zone.MeshLinePttrnLen[GroupIndex] =			2;		//	Value is in %	//	double
		Zone.MeshLineThickness[GroupIndex] =			0.1;	//	Value is in %	//	double


		//	Contour Settings

		Zone.ContourShow[GroupIndex] =				TRUE;						//	Boolean_t
		Zone.ContourType[GroupIndex] =				ContourType_Flood;			//	ContourType_e
		Zone.FloodColoring[GroupIndex] =				ContourColoring_Group1;		//	ContourColoring_e
		Zone.ContourLineGroup[GroupIndex] =			1;							//	SmInteger_t
		Zone.ContourLineColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.ContourLinePttrn[GroupIndex] =			LinePattern_Solid;			//	LinePattern_e
		Zone.ContourLinePttrnLen[GroupIndex] =		2;		//	Value is in %	//	double
		Zone.ContourLineThickness[GroupIndex] =		0.1;	//	Value is in %	//	double
		Zone.ContourLighting[GroupIndex] =			TRUE;						//	Boolean_t

		//	Vector Settings
		
		Zone.VectorShow[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.VectorType[GroupIndex] =				VectorType_TailAtPoint;		//	VectorType_e
		Zone.VectorArrowType[GroupIndex] =			ArrowheadStyle_Plain;		//	ArrowheadStyle_e
		Zone.VectorLineColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.VectorTangent[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.VectorLinePttrn[GroupIndex] =			LinePattern_Solid;			//	LinePattern_e
		Zone.VectorLinePttrnLen[GroupIndex] =		2;		//	Value is in %	//	double
		Zone.VectorLineThickness[GroupIndex] =		0.1;	//	Value is in %	//	double

		//	Scatter Settings

		Zone.ScatterShow[GroupIndex] =				TRUE;						//	Boolean_t
		Zone.ScatterShape[GroupIndex] =				GeomShape_Sphere;			//	GeomShape_e
		Zone.ScatterOutColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.ScatterInColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.ScatterFillMode[GroupIndex] =			FillMode_UseLineColor;		//	FillMode_e
		Zone.ScatterSize[GroupIndex] =				2.5;	//	value is in %	//	double
		Zone.ScatterLineThickness[GroupIndex] =		0.1;	//	value is in %	//	double

		//	Shade Settings

		Zone.ShadeShow[GroupIndex] =					TRUE;						//	Boolean_t
		Zone.ShadeColor[GroupIndex] =				Black_C;					//	ColorIndex_t
		Zone.ShadeLighting[GroupIndex] =				TRUE;						//	Boolean_t

		//	Edge Settings

		Zone.EdgeShow[GroupIndex] =					TRUE;						//	Boolean_t
		Zone.EdgeType[GroupIndex] =					EdgeType_BordersAndCreases;	//	EdgeType_e
		Zone.EdgeIBorder[GroupIndex] =				BorderLocation_Both;		//	BorderLocation_e
		Zone.EdgeJBorder[GroupIndex] =				BorderLocation_Both;		//	BorderLocation_e
		Zone.EdgeKBorder[GroupIndex] =				BorderLocation_Both;		//	BorderLocation_e
		Zone.EdgeColor[GroupIndex] =					Black_C;					//	ColorIndex_t
		Zone.EdgeLineThickness[GroupIndex] =			0.1;	//	value is in %	//	double

		//	Points Settings

		Zone.PointsToPlot[GroupIndex] =				PointsToPlot_SurfaceNodes;	//	PointsToPlot_e
		Zone.PointsIndexSkipI[GroupIndex] =			1;							//	LgIndex_t
		Zone.PointsIndexSkipJ[GroupIndex] =			1;							//	LgIndex_t
		Zone.PointsIndexSkipK[GroupIndex] =			1;							//	LgIndex_t
		/*
			1 is the same as "No Skip" in the GUI.
		*/

		//	Surfaces Settings

		Zone.SurfacesToPlot[GroupIndex] =			SurfacesToPlot_ExposedCellFaces;	//	SurfacesToPlot_e
		Zone.RangeBeginI[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeEndI[GroupIndex] =					0;							//	LgIndex_t
		Zone.RangeSkipI[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeBeginJ[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeEndJ[GroupIndex] =					0;							//	LgIndex_t
		Zone.RangeSkipJ[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeBeginK[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeEndK[GroupIndex] =					0;							//	LgIndex_t
		Zone.RangeSkipK[GroupIndex] =				1;							//	LgIndex_t
		/*
			0 is the same as maximum, which sets the end range to the maximum
			index for the dimension. The skip number is to skip every # planes,
			so 1 is to not skip any planes, (show every 1 plane).
		*/

		//	Volume Settings

		Zone.StreamtracesShow[GroupIndex] =			TRUE;						//	Boolean_t
		Zone.IsoSurfacesShow[GroupIndex] =			TRUE;						//	Boolean_t
		Zone.SlicesShow[GroupIndex] =				TRUE;						//	Boolean_t

		//	Effects Settings

		Zone.UseValueBlanking[GroupIndex] =			TRUE;						//	Boolean_t
		//Zone.UseClipPlanes[GroupIndex] =			TRUE;						//	Boolean_t
		Zone.UseSurfaceTranslucency[GroupIndex] =	TRUE;						//	Boolean_t
		Zone.SurfaceTranslucency[GroupIndex] =		50;		//	value is in %	//	SmInteger_t
		Zone.LightingEffect[GroupIndex] =			LightingEffect_Gouraud;		//	LightingEffect_e

		//	Group index and number (increments for next set)

		GroupIndex++;														//	int


		//===============================================================
		//	BondZones Settings
		//===============================================================

		Zone.ZoneShow[GroupIndex] =					FALSE;						//	Boolean_t

		//	Mesh Settings

		Zone.MeshShow[GroupIndex] =					TRUE;						//	Boolean_t
		Zone.MeshType[GroupIndex] =					MeshType_Overlay;			//	MeshType_e
		Zone.MeshColor[GroupIndex] =					Black_C;					//	ColorIndex_t
		/*
			A number of color index constants are \#defined. These include Black_C, Red_C,
			Green_C, Blue_C, Cyan_C, Yellow_C, Purple_C, White_C, Custom1_C through
			Custom56_C, MultiColor_C, NoColor_C, MulitiColor2_C through MulitiColor8_C,
			RGBColor_C, and InvalidColor_C.  
		*/
		Zone.MeshLinePttrn[GroupIndex] =				LinePattern_Solid;			//	LinePattern_e
		Zone.MeshLinePttrnLen[GroupIndex] =			2;		//	Value is in %	//	double
		Zone.MeshLineThickness[GroupIndex] =			0.1;	//	Value is in %	//	double


		//	Contour Settings

		Zone.ContourShow[GroupIndex] =				TRUE;						//	Boolean_t
		Zone.ContourType[GroupIndex] =				ContourType_Flood;			//	ContourType_e
		Zone.FloodColoring[GroupIndex] =				ContourColoring_Group1;		//	ContourColoring_e
		Zone.ContourLineGroup[GroupIndex] =			1;							//	SmInteger_t
		Zone.ContourLineColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.ContourLinePttrn[GroupIndex] =			LinePattern_Solid;			//	LinePattern_e
		Zone.ContourLinePttrnLen[GroupIndex] =		2;		//	Value is in %	//	double
		Zone.ContourLineThickness[GroupIndex] =		0.1;	//	Value is in %	//	double
		Zone.ContourLighting[GroupIndex] =			TRUE;						//	Boolean_t

		//	Vector Settings
		
		Zone.VectorShow[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.VectorType[GroupIndex] =				VectorType_TailAtPoint;		//	VectorType_e
		Zone.VectorArrowType[GroupIndex] =			ArrowheadStyle_Plain;		//	ArrowheadStyle_e
		Zone.VectorLineColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.VectorTangent[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.VectorLinePttrn[GroupIndex] =			LinePattern_Solid;			//	LinePattern_e
		Zone.VectorLinePttrnLen[GroupIndex] =		2;		//	Value is in %	//	double
		Zone.VectorLineThickness[GroupIndex] =		0.1;	//	Value is in %	//	double

		//	Scatter Settings

		Zone.ScatterShow[GroupIndex] =				TRUE;						//	Boolean_t
		Zone.ScatterShape[GroupIndex] =				GeomShape_Sphere;			//	GeomShape_e
		Zone.ScatterOutColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.ScatterInColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.ScatterFillMode[GroupIndex] =			FillMode_UseLineColor;		//	FillMode_e
		Zone.ScatterSize[GroupIndex] =				2.5;	//	value is in %	//	double
		Zone.ScatterLineThickness[GroupIndex] =		0.1;	//	value is in %	//	double

		//	Shade Settings

		Zone.ShadeShow[GroupIndex] =					TRUE;						//	Boolean_t
		Zone.ShadeColor[GroupIndex] =				Black_C;					//	ColorIndex_t
		Zone.ShadeLighting[GroupIndex] =				TRUE;						//	Boolean_t

		//	Edge Settings

		Zone.EdgeShow[GroupIndex] =					TRUE;						//	Boolean_t
		Zone.EdgeType[GroupIndex] =					EdgeType_BordersAndCreases;	//	EdgeType_e
		Zone.EdgeIBorder[GroupIndex] =				BorderLocation_Both;		//	BorderLocation_e
		Zone.EdgeJBorder[GroupIndex] =				BorderLocation_Both;		//	BorderLocation_e
		Zone.EdgeKBorder[GroupIndex] =				BorderLocation_Both;		//	BorderLocation_e
		Zone.EdgeColor[GroupIndex] =					Black_C;					//	ColorIndex_t
		Zone.EdgeLineThickness[GroupIndex] =			0.1;	//	value is in %	//	double

		//	Points Settings

		Zone.PointsToPlot[GroupIndex] =				PointsToPlot_SurfaceNodes;	//	PointsToPlot_e
		Zone.PointsIndexSkipI[GroupIndex] =			1;							//	LgIndex_t
		Zone.PointsIndexSkipJ[GroupIndex] =			1;							//	LgIndex_t
		Zone.PointsIndexSkipK[GroupIndex] =			1;							//	LgIndex_t
		/*
			1 is the same as "No Skip" in the GUI.
		*/

		//	Surfaces Settings

		Zone.SurfacesToPlot[GroupIndex] =			SurfacesToPlot_ExposedCellFaces;	//	SurfacesToPlot_e
		Zone.RangeBeginI[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeEndI[GroupIndex] =					0;							//	LgIndex_t
		Zone.RangeSkipI[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeBeginJ[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeEndJ[GroupIndex] =					0;							//	LgIndex_t
		Zone.RangeSkipJ[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeBeginK[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeEndK[GroupIndex] =					0;							//	LgIndex_t
		Zone.RangeSkipK[GroupIndex] =				1;							//	LgIndex_t
		/*
			0 is the same as maximum, which sets the end range to the maximum
			index for the dimension. The skip number is to skip every # planes,
			so 1 is to not skip any planes, (show every 1 plane).
		*/

		//	Volume Settings

		Zone.StreamtracesShow[GroupIndex] =			TRUE;						//	Boolean_t
		Zone.IsoSurfacesShow[GroupIndex] =			TRUE;						//	Boolean_t
		Zone.SlicesShow[GroupIndex] =				TRUE;						//	Boolean_t

		//	Effects Settings

		Zone.UseValueBlanking[GroupIndex] =			TRUE;						//	Boolean_t
		//Zone.UseClipPlanes[GroupIndex] =			TRUE;						//	Boolean_t
		Zone.UseSurfaceTranslucency[GroupIndex] =	TRUE;						//	Boolean_t
		Zone.SurfaceTranslucency[GroupIndex] =		50;		//	value is in %	//	SmInteger_t
		Zone.LightingEffect[GroupIndex] =			LightingEffect_Gouraud;		//	LightingEffect_e

		//	Group index and number (increments for next set)

		GroupIndex++;														//	int


		//===============================================================
		//	BondZonesFF Settings
		//===============================================================

		Zone.ZoneShow[GroupIndex] =					FALSE;						//	Boolean_t

		//	Mesh Settings

		Zone.MeshShow[GroupIndex] =					TRUE;						//	Boolean_t
		Zone.MeshType[GroupIndex] =					MeshType_Overlay;			//	MeshType_e
		Zone.MeshColor[GroupIndex] =					Black_C;					//	ColorIndex_t
		/*
			A number of color index constants are \#defined. These include Black_C, Red_C,
			Green_C, Blue_C, Cyan_C, Yellow_C, Purple_C, White_C, Custom1_C through
			Custom56_C, MultiColor_C, NoColor_C, MulitiColor2_C through MulitiColor8_C,
			RGBColor_C, and InvalidColor_C.  
		*/
		Zone.MeshLinePttrn[GroupIndex] =				LinePattern_Solid;			//	LinePattern_e
		Zone.MeshLinePttrnLen[GroupIndex] =			2;		//	Value is in %	//	double
		Zone.MeshLineThickness[GroupIndex] =			0.1;	//	Value is in %	//	double


		//	Contour Settings

		Zone.ContourShow[GroupIndex] =				TRUE;						//	Boolean_t
		Zone.ContourType[GroupIndex] =				ContourType_Flood;			//	ContourType_e
		Zone.FloodColoring[GroupIndex] =				ContourColoring_Group1;		//	ContourColoring_e
		Zone.ContourLineGroup[GroupIndex] =			1;							//	SmInteger_t
		Zone.ContourLineColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.ContourLinePttrn[GroupIndex] =			LinePattern_Solid;			//	LinePattern_e
		Zone.ContourLinePttrnLen[GroupIndex] =		2;		//	Value is in %	//	double
		Zone.ContourLineThickness[GroupIndex] =		0.1;	//	Value is in %	//	double
		Zone.ContourLighting[GroupIndex] =			TRUE;						//	Boolean_t

		//	Vector Settings
		
		Zone.VectorShow[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.VectorType[GroupIndex] =				VectorType_TailAtPoint;		//	VectorType_e
		Zone.VectorArrowType[GroupIndex] =			ArrowheadStyle_Plain;		//	ArrowheadStyle_e
		Zone.VectorLineColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.VectorTangent[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.VectorLinePttrn[GroupIndex] =			LinePattern_Solid;			//	LinePattern_e
		Zone.VectorLinePttrnLen[GroupIndex] =		2;		//	Value is in %	//	double
		Zone.VectorLineThickness[GroupIndex] =		0.1;	//	Value is in %	//	double

		//	Scatter Settings

		Zone.ScatterShow[GroupIndex] =				TRUE;						//	Boolean_t
		Zone.ScatterShape[GroupIndex] =				GeomShape_Sphere;			//	GeomShape_e
		Zone.ScatterOutColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.ScatterInColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.ScatterFillMode[GroupIndex] =			FillMode_UseLineColor;		//	FillMode_e
		Zone.ScatterSize[GroupIndex] =				2.5;	//	value is in %	//	double
		Zone.ScatterLineThickness[GroupIndex] =		0.1;	//	value is in %	//	double

		//	Shade Settings

		Zone.ShadeShow[GroupIndex] =					TRUE;						//	Boolean_t
		Zone.ShadeColor[GroupIndex] =				Black_C;					//	ColorIndex_t
		Zone.ShadeLighting[GroupIndex] =				TRUE;						//	Boolean_t

		//	Edge Settings

		Zone.EdgeShow[GroupIndex] =					TRUE;						//	Boolean_t
		Zone.EdgeType[GroupIndex] =					EdgeType_BordersAndCreases;	//	EdgeType_e
		Zone.EdgeIBorder[GroupIndex] =				BorderLocation_Both;		//	BorderLocation_e
		Zone.EdgeJBorder[GroupIndex] =				BorderLocation_Both;		//	BorderLocation_e
		Zone.EdgeKBorder[GroupIndex] =				BorderLocation_Both;		//	BorderLocation_e
		Zone.EdgeColor[GroupIndex] =					Black_C;					//	ColorIndex_t
		Zone.EdgeLineThickness[GroupIndex] =			0.1;	//	value is in %	//	double

		//	Points Settings

		Zone.PointsToPlot[GroupIndex] =				PointsToPlot_SurfaceNodes;	//	PointsToPlot_e
		Zone.PointsIndexSkipI[GroupIndex] =			1;							//	LgIndex_t
		Zone.PointsIndexSkipJ[GroupIndex] =			1;							//	LgIndex_t
		Zone.PointsIndexSkipK[GroupIndex] =			1;							//	LgIndex_t
		/*
			1 is the same as "No Skip" in the GUI.
		*/

		//	Surfaces Settings

		Zone.SurfacesToPlot[GroupIndex] =			SurfacesToPlot_ExposedCellFaces;	//	SurfacesToPlot_e
		Zone.RangeBeginI[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeEndI[GroupIndex] =					0;							//	LgIndex_t
		Zone.RangeSkipI[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeBeginJ[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeEndJ[GroupIndex] =					0;							//	LgIndex_t
		Zone.RangeSkipJ[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeBeginK[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeEndK[GroupIndex] =					0;							//	LgIndex_t
		Zone.RangeSkipK[GroupIndex] =				1;							//	LgIndex_t
		/*
			0 is the same as maximum, which sets the end range to the maximum
			index for the dimension. The skip number is to skip every # planes,
			so 1 is to not skip any planes, (show every 1 plane).
		*/

		//	Volume Settings

		Zone.StreamtracesShow[GroupIndex] =			TRUE;						//	Boolean_t
		Zone.IsoSurfacesShow[GroupIndex] =			TRUE;						//	Boolean_t
		Zone.SlicesShow[GroupIndex] =				TRUE;						//	Boolean_t

		//	Effects Settings

		Zone.UseValueBlanking[GroupIndex] =			TRUE;						//	Boolean_t
		//Zone.UseClipPlanes[GroupIndex] =			TRUE;						//	Boolean_t
		Zone.UseSurfaceTranslucency[GroupIndex] =	TRUE;						//	Boolean_t
		Zone.SurfaceTranslucency[GroupIndex] =		50;		//	value is in %	//	SmInteger_t
		Zone.LightingEffect[GroupIndex] =			LightingEffect_Gouraud;		//	LightingEffect_e

		//	Group index and number (increments for next set)

		GroupIndex++;														//	int

		//===============================================================
		//	BadZones Settings
		//===============================================================

		Zone.ZoneShow[GroupIndex] =					FALSE;						//	Boolean_t

		//	Mesh Settings

		Zone.MeshShow[GroupIndex] =					FALSE;						//	Boolean_t
		Zone.MeshType[GroupIndex] =					MeshType_Overlay;			//	MeshType_e
		Zone.MeshColor[GroupIndex] =					Black_C;					//	ColorIndex_t
		/*
			A number of color index constants are \#defined. These include Black_C, Red_C,
			Green_C, Blue_C, Cyan_C, Yellow_C, Purple_C, White_C, Custom1_C through
			Custom56_C, MultiColor_C, NoColor_C, MulitiColor2_C through MulitiColor8_C,
			RGBColor_C, and InvalidColor_C.  
		*/
		Zone.MeshLinePttrn[GroupIndex] =				LinePattern_Solid;			//	LinePattern_e
		Zone.MeshLinePttrnLen[GroupIndex] =			2;		//	Value is in %	//	double
		Zone.MeshLineThickness[GroupIndex] =			0.1;	//	Value is in %	//	double


		//	Contour Settings

		Zone.ContourShow[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.ContourType[GroupIndex] =				ContourType_Flood;			//	ContourType_e
		Zone.FloodColoring[GroupIndex] =				ContourColoring_Group1;		//	ContourColoring_e
		Zone.ContourLineGroup[GroupIndex] =			1;							//	SmInteger_t
		Zone.ContourLineColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.ContourLinePttrn[GroupIndex] =			LinePattern_Solid;			//	LinePattern_e
		Zone.ContourLinePttrnLen[GroupIndex] =		2;		//	Value is in %	//	double
		Zone.ContourLineThickness[GroupIndex] =		0.1;	//	Value is in %	//	double
		Zone.ContourLighting[GroupIndex] =			TRUE;						//	Boolean_t

		//	Vector Settings
		
		Zone.VectorShow[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.VectorType[GroupIndex] =				VectorType_TailAtPoint;		//	VectorType_e
		Zone.VectorArrowType[GroupIndex] =			ArrowheadStyle_Plain;		//	ArrowheadStyle_e
		Zone.VectorLineColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.VectorTangent[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.VectorLinePttrn[GroupIndex] =			LinePattern_Solid;			//	LinePattern_e
		Zone.VectorLinePttrnLen[GroupIndex] =		2;		//	Value is in %	//	double
		Zone.VectorLineThickness[GroupIndex] =		0.1;	//	Value is in %	//	double

		//	Scatter Settings

		Zone.ScatterShow[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.ScatterShape[GroupIndex] =				GeomShape_Sphere;			//	GeomShape_e
		Zone.ScatterOutColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.ScatterInColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.ScatterFillMode[GroupIndex] =			FillMode_UseLineColor;		//	FillMode_e
		Zone.ScatterSize[GroupIndex] =				2.5;	//	value is in %	//	double
		Zone.ScatterLineThickness[GroupIndex] =		0.1;	//	value is in %	//	double

		//	Shade Settings

		Zone.ShadeShow[GroupIndex] =					FALSE;						//	Boolean_t
		Zone.ShadeColor[GroupIndex] =				Black_C;					//	ColorIndex_t
		Zone.ShadeLighting[GroupIndex] =				TRUE;						//	Boolean_t

		//	Edge Settings

		Zone.EdgeShow[GroupIndex] =					FALSE;						//	Boolean_t
		Zone.EdgeType[GroupIndex] =					EdgeType_BordersAndCreases;	//	EdgeType_e
		Zone.EdgeIBorder[GroupIndex] =				BorderLocation_Both;		//	BorderLocation_e
		Zone.EdgeJBorder[GroupIndex] =				BorderLocation_Both;		//	BorderLocation_e
		Zone.EdgeKBorder[GroupIndex] =				BorderLocation_Both;		//	BorderLocation_e
		Zone.EdgeColor[GroupIndex] =					Black_C;					//	ColorIndex_t
		Zone.EdgeLineThickness[GroupIndex] =			0.1;	//	value is in %	//	double

		//	Points Settings

		Zone.PointsToPlot[GroupIndex] =				PointsToPlot_SurfaceNodes;	//	PointsToPlot_e
		Zone.PointsIndexSkipI[GroupIndex] =			1;							//	LgIndex_t
		Zone.PointsIndexSkipJ[GroupIndex] =			1;							//	LgIndex_t
		Zone.PointsIndexSkipK[GroupIndex] =			1;							//	LgIndex_t
		/*
			1 is the same as "No Skip" in the GUI.
		*/

		//	Surfaces Settings

		Zone.SurfacesToPlot[GroupIndex] =			SurfacesToPlot_ExposedCellFaces;	//	SurfacesToPlot_e
		Zone.RangeBeginI[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeEndI[GroupIndex] =					0;							//	LgIndex_t
		Zone.RangeSkipI[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeBeginJ[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeEndJ[GroupIndex] =					0;							//	LgIndex_t
		Zone.RangeSkipJ[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeBeginK[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeEndK[GroupIndex] =					0;							//	LgIndex_t
		Zone.RangeSkipK[GroupIndex] =				1;							//	LgIndex_t
		/*
			0 is the same as maximum, which sets the end range to the maximum
			index for the dimension. The skip number is to skip every # planes,
			so 1 is to not skip any planes, (show every 1 plane).
		*/

		//	Volume Settings

		Zone.StreamtracesShow[GroupIndex] =			TRUE;						//	Boolean_t
		Zone.IsoSurfacesShow[GroupIndex] =			TRUE;						//	Boolean_t
		Zone.SlicesShow[GroupIndex] =				TRUE;						//	Boolean_t

		//	Effects Settings

		Zone.UseValueBlanking[GroupIndex] =			TRUE;						//	Boolean_t
		//Zone.UseClipPlanes[GroupIndex] =			TRUE;						//	Boolean_t
		Zone.UseSurfaceTranslucency[GroupIndex] =	TRUE;						//	Boolean_t
		Zone.SurfaceTranslucency[GroupIndex] =		50;		//	value is in %	//	SmInteger_t
		Zone.LightingEffect[GroupIndex] =			LightingEffect_Gouraud;		//	LightingEffect_e

		//	Group index and number (increments for next set)

		GroupIndex++;														//	int


		ApplyZoneSettings(ZoneSets, Zone);

		DeallocSets(ZoneSets);

		//char FinishedMsg[200];
		//sprintf_s(FinishedMsg, "Finished Successfully");
		//TecUtilDialogMessageBox(FinishedMsg, MessageBox_Information);
	}
	else
	{
		TecUtilDialogErrMsg("Data set not available.");
	}

} //	void ChemSysView()



void SGradTopoView()
{
	if (TecUtilDataSetIsAvailable())
	{
		std::vector<Set_pa> ZoneSets;

		GroupZones(ZoneSets);

		SetupMultiColorContourLevels(ZoneSets);

		ZoneSettings	Zone;

		int				GroupIndex = 0;
		
		//***************************************************************
		//===============================================================
		//	Adjust properties for each group.
		//===============================================================
		//***************************************************************

		//===============================================================
		//	InitialDataSetZone Settings
		//===============================================================

		Zone.ZoneShow[GroupIndex] =					FALSE;						//	Boolean_t

		//	Mesh Settings

		Zone.MeshShow[GroupIndex] =					FALSE;						//	Boolean_t
		Zone.MeshType[GroupIndex] =					MeshType_Overlay;			//	MeshType_e
		Zone.MeshColor[GroupIndex] =					Black_C;					//	ColorIndex_t
		/*
			A number of color index constants are \#defined. These include Black_C, Red_C,
			Green_C, Blue_C, Cyan_C, Yellow_C, Purple_C, White_C, Custom1_C through
			Custom56_C, MultiColor_C, NoColor_C, MulitiColor2_C through MulitiColor8_C,
			RGBColor_C, and InvalidColor_C.  
		*/
		Zone.MeshLinePttrn[GroupIndex] =				LinePattern_Solid;			//	LinePattern_e
		Zone.MeshLinePttrnLen[GroupIndex] =			2;		//	Value is in %	//	double
		Zone.MeshLineThickness[GroupIndex] =			0.1;	//	Value is in %	//	double


		//	Contour Settings

		Zone.ContourShow[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.ContourType[GroupIndex] =				ContourType_Flood;			//	ContourType_e
		Zone.FloodColoring[GroupIndex] =				ContourColoring_Group1;		//	ContourColoring_e
		Zone.ContourLineGroup[GroupIndex] =			1;							//	SmInteger_t
		Zone.ContourLineColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.ContourLinePttrn[GroupIndex] =			LinePattern_Solid;			//	LinePattern_e
		Zone.ContourLinePttrnLen[GroupIndex] =		2;		//	Value is in %	//	double
		Zone.ContourLineThickness[GroupIndex] =		0.1;	//	Value is in %	//	double
		Zone.ContourLighting[GroupIndex] =			TRUE;						//	Boolean_t

		//	Vector Settings
		
		Zone.VectorShow[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.VectorType[GroupIndex] =				VectorType_TailAtPoint;		//	VectorType_e
		Zone.VectorArrowType[GroupIndex] =			ArrowheadStyle_Plain;		//	ArrowheadStyle_e
		Zone.VectorLineColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.VectorTangent[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.VectorLinePttrn[GroupIndex] =			LinePattern_Solid;			//	LinePattern_e
		Zone.VectorLinePttrnLen[GroupIndex] =		2;		//	Value is in %	//	double
		Zone.VectorLineThickness[GroupIndex] =		0.1;	//	Value is in %	//	double

		//	Scatter Settings

		Zone.ScatterShow[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.ScatterShape[GroupIndex] =				GeomShape_Sphere;			//	GeomShape_e
		Zone.ScatterOutColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.ScatterInColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.ScatterFillMode[GroupIndex] =			FillMode_UseLineColor;		//	FillMode_e
		Zone.ScatterSize[GroupIndex] =				2.5;	//	value is in %	//	double
		Zone.ScatterLineThickness[GroupIndex] =		0.1;	//	value is in %	//	double

		//	Shade Settings

		Zone.ShadeShow[GroupIndex] =					FALSE;						//	Boolean_t
		Zone.ShadeColor[GroupIndex] =				Black_C;					//	ColorIndex_t
		Zone.ShadeLighting[GroupIndex] =				TRUE;						//	Boolean_t

		//	Edge Settings

		Zone.EdgeShow[GroupIndex] =					FALSE;						//	Boolean_t
		Zone.EdgeType[GroupIndex] =					EdgeType_BordersAndCreases;	//	EdgeType_e
		Zone.EdgeIBorder[GroupIndex] =				BorderLocation_Both;		//	BorderLocation_e
		Zone.EdgeJBorder[GroupIndex] =				BorderLocation_Both;		//	BorderLocation_e
		Zone.EdgeKBorder[GroupIndex] =				BorderLocation_Both;		//	BorderLocation_e
		Zone.EdgeColor[GroupIndex] =					Black_C;					//	ColorIndex_t
		Zone.EdgeLineThickness[GroupIndex] =			0.1;	//	value is in %	//	double

		//	Points Settings

		Zone.PointsToPlot[GroupIndex] =				PointsToPlot_SurfaceNodes;	//	PointsToPlot_e
		Zone.PointsIndexSkipI[GroupIndex] =			1;							//	LgIndex_t
		Zone.PointsIndexSkipJ[GroupIndex] =			1;							//	LgIndex_t
		Zone.PointsIndexSkipK[GroupIndex] =			1;							//	LgIndex_t
		/*
			1 is the same as "No Skip" in the GUI.
		*/

		//	Surfaces Settings

		Zone.SurfacesToPlot[GroupIndex] =			SurfacesToPlot_ExposedCellFaces;	//	SurfacesToPlot_e
		Zone.RangeBeginI[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeEndI[GroupIndex] =					0;							//	LgIndex_t
		Zone.RangeSkipI[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeBeginJ[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeEndJ[GroupIndex] =					0;							//	LgIndex_t
		Zone.RangeSkipJ[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeBeginK[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeEndK[GroupIndex] =					0;							//	LgIndex_t
		Zone.RangeSkipK[GroupIndex] =				1;							//	LgIndex_t
		/*
			0 is the same as maximum, which sets the end range to the maximum
			index for the dimension. The skip number is to skip every # planes,
			so 1 is to not skip any planes, (show every 1 plane).
		*/

		//	Volume Settings

		Zone.StreamtracesShow[GroupIndex] =			TRUE;						//	Boolean_t
		Zone.IsoSurfacesShow[GroupIndex] =			TRUE;						//	Boolean_t
		Zone.SlicesShow[GroupIndex] =				TRUE;						//	Boolean_t

		//	Effects Settings

		Zone.UseValueBlanking[GroupIndex] =			TRUE;						//	Boolean_t
		//Zone.UseClipPlanes[GroupIndex] =			TRUE;						//	Boolean_t
		Zone.UseSurfaceTranslucency[GroupIndex] =	TRUE;						//	Boolean_t
		Zone.SurfaceTranslucency[GroupIndex] =		50;		//	value is in %	//	SmInteger_t
		Zone.LightingEffect[GroupIndex] =			LightingEffect_Gouraud;		//	LightingEffect_e

		//	Group index and number (increments for next set)

		GroupIndex++;														//	int

		//===============================================================
		//	VolumeCPs Settings
		//===============================================================

		Zone.ZoneShow[GroupIndex] =					TRUE;						//	Boolean_t

		//	Mesh Settings

		Zone.MeshShow[GroupIndex] =					FALSE;						//	Boolean_t
		Zone.MeshType[GroupIndex] =					MeshType_Overlay;			//	MeshType_e
		Zone.MeshColor[GroupIndex] =					Black_C;					//	ColorIndex_t
		/*
			A number of color index constants are \#defined. These include Black_C, Red_C,
			Green_C, Blue_C, Cyan_C, Yellow_C, Purple_C, White_C, Custom1_C through
			Custom56_C, MultiColor_C, NoColor_C, MulitiColor2_C through MulitiColor8_C,
			RGBColor_C, and InvalidColor_C.  
		*/
		Zone.MeshLinePttrn[GroupIndex] =				LinePattern_Solid;			//	LinePattern_e
		Zone.MeshLinePttrnLen[GroupIndex] =			2;		//	Value is in %	//	double
		Zone.MeshLineThickness[GroupIndex] =			0.1;	//	Value is in %	//	double


		//	Contour Settings

		Zone.ContourShow[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.ContourType[GroupIndex] =				ContourType_Flood;			//	ContourType_e
		Zone.FloodColoring[GroupIndex] =				ContourColoring_Group1;		//	ContourColoring_e
		Zone.ContourLineGroup[GroupIndex] =			1;							//	SmInteger_t
		Zone.ContourLineColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.ContourLinePttrn[GroupIndex] =			LinePattern_Solid;			//	LinePattern_e
		Zone.ContourLinePttrnLen[GroupIndex] =		2;		//	Value is in %	//	double
		Zone.ContourLineThickness[GroupIndex] =		0.1;	//	Value is in %	//	double
		Zone.ContourLighting[GroupIndex] =			TRUE;						//	Boolean_t

		//	Vector Settings
		
		Zone.VectorShow[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.VectorType[GroupIndex] =				VectorType_TailAtPoint;		//	VectorType_e
		Zone.VectorArrowType[GroupIndex] =			ArrowheadStyle_Plain;		//	ArrowheadStyle_e
		Zone.VectorLineColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.VectorTangent[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.VectorLinePttrn[GroupIndex] =			LinePattern_Solid;			//	LinePattern_e
		Zone.VectorLinePttrnLen[GroupIndex] =		2;		//	Value is in %	//	double
		Zone.VectorLineThickness[GroupIndex] =		0.1;	//	Value is in %	//	double

		//	Scatter Settings

		Zone.ScatterShow[GroupIndex] =				TRUE;						//	Boolean_t
		Zone.ScatterShape[GroupIndex] =				GeomShape_Sphere;			//	GeomShape_e
		Zone.ScatterOutColor[GroupIndex] =			MultiColor_C;				//	ColorIndex_t
		Zone.ScatterInColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.ScatterFillMode[GroupIndex] =			FillMode_UseLineColor;		//	FillMode_e
		Zone.ScatterSize[GroupIndex] =				2.5;	//	value is in %	//	double
		Zone.ScatterLineThickness[GroupIndex] =		0.1;	//	value is in %	//	double

		//	Shade Settings

		Zone.ShadeShow[GroupIndex] =					FALSE;						//	Boolean_t
		Zone.ShadeColor[GroupIndex] =				Black_C;					//	ColorIndex_t
		Zone.ShadeLighting[GroupIndex] =				TRUE;						//	Boolean_t

		//	Edge Settings

		Zone.EdgeShow[GroupIndex] =					FALSE;						//	Boolean_t
		Zone.EdgeType[GroupIndex] =					EdgeType_BordersAndCreases;	//	EdgeType_e
		Zone.EdgeIBorder[GroupIndex] =				BorderLocation_Both;		//	BorderLocation_e
		Zone.EdgeJBorder[GroupIndex] =				BorderLocation_Both;		//	BorderLocation_e
		Zone.EdgeKBorder[GroupIndex] =				BorderLocation_Both;		//	BorderLocation_e
		Zone.EdgeColor[GroupIndex] =					Black_C;					//	ColorIndex_t
		Zone.EdgeLineThickness[GroupIndex] =			0.1;	//	value is in %	//	double

		//	Points Settings

		Zone.PointsToPlot[GroupIndex] =				PointsToPlot_SurfaceNodes;	//	PointsToPlot_e
		Zone.PointsIndexSkipI[GroupIndex] =			1;							//	LgIndex_t
		Zone.PointsIndexSkipJ[GroupIndex] =			1;							//	LgIndex_t
		Zone.PointsIndexSkipK[GroupIndex] =			1;							//	LgIndex_t
		/*
			1 is the same as "No Skip" in the GUI.
		*/

		//	Surfaces Settings

		Zone.SurfacesToPlot[GroupIndex] =			SurfacesToPlot_ExposedCellFaces;	//	SurfacesToPlot_e
		Zone.RangeBeginI[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeEndI[GroupIndex] =					0;							//	LgIndex_t
		Zone.RangeSkipI[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeBeginJ[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeEndJ[GroupIndex] =					0;							//	LgIndex_t
		Zone.RangeSkipJ[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeBeginK[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeEndK[GroupIndex] =					0;							//	LgIndex_t
		Zone.RangeSkipK[GroupIndex] =				1;							//	LgIndex_t
		/*
			0 is the same as maximum, which sets the end range to the maximum
			index for the dimension. The skip number is to skip every # planes,
			so 1 is to not skip any planes, (show every 1 plane).
		*/

		//	Volume Settings

		Zone.StreamtracesShow[GroupIndex] =			TRUE;						//	Boolean_t
		Zone.IsoSurfacesShow[GroupIndex] =			TRUE;						//	Boolean_t
		Zone.SlicesShow[GroupIndex] =				TRUE;						//	Boolean_t

		//	Effects Settings

		Zone.UseValueBlanking[GroupIndex] =			FALSE;						//	Boolean_t
		//Zone.UseClipPlanes[GroupIndex] =			TRUE;						//	Boolean_t
		Zone.UseSurfaceTranslucency[GroupIndex] =	FALSE;						//	Boolean_t
		Zone.SurfaceTranslucency[GroupIndex] =		50;		//	value is in %	//	SmInteger_t
		Zone.LightingEffect[GroupIndex] =			LightingEffect_Gouraud;		//	LightingEffect_e

		//	Group index and number (increments for next set)

		GroupIndex++;														//	int


		//===============================================================
		//	IsoSurfCPs Settings
		//===============================================================

		Zone.ZoneShow[GroupIndex] =					TRUE;						//	Boolean_t

		//	Mesh Settings

		Zone.MeshShow[GroupIndex] =					FALSE;						//	Boolean_t
		Zone.MeshType[GroupIndex] =					MeshType_Overlay;			//	MeshType_e
		Zone.MeshColor[GroupIndex] =					Black_C;					//	ColorIndex_t
		/*
			A number of color index constants are \#defined. These include Black_C, Red_C,
			Green_C, Blue_C, Cyan_C, Yellow_C, Purple_C, White_C, Custom1_C through
			Custom56_C, MultiColor_C, NoColor_C, MulitiColor2_C through MulitiColor8_C,
			RGBColor_C, and InvalidColor_C.  
		*/
		Zone.MeshLinePttrn[GroupIndex] =				LinePattern_Solid;			//	LinePattern_e
		Zone.MeshLinePttrnLen[GroupIndex] =			2;		//	Value is in %	//	double
		Zone.MeshLineThickness[GroupIndex] =			0.1;	//	Value is in %	//	double


		//	Contour Settings

		Zone.ContourShow[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.ContourType[GroupIndex] =				ContourType_Flood;			//	ContourType_e
		Zone.FloodColoring[GroupIndex] =				ContourColoring_Group1;		//	ContourColoring_e
		Zone.ContourLineGroup[GroupIndex] =			1;							//	SmInteger_t
		Zone.ContourLineColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.ContourLinePttrn[GroupIndex] =			LinePattern_Solid;			//	LinePattern_e
		Zone.ContourLinePttrnLen[GroupIndex] =		2;		//	Value is in %	//	double
		Zone.ContourLineThickness[GroupIndex] =		0.1;	//	Value is in %	//	double
		Zone.ContourLighting[GroupIndex] =			TRUE;						//	Boolean_t

		//	Vector Settings
		
		Zone.VectorShow[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.VectorType[GroupIndex] =				VectorType_TailAtPoint;		//	VectorType_e
		Zone.VectorArrowType[GroupIndex] =			ArrowheadStyle_Plain;		//	ArrowheadStyle_e
		Zone.VectorLineColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.VectorTangent[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.VectorLinePttrn[GroupIndex] =			LinePattern_Solid;			//	LinePattern_e
		Zone.VectorLinePttrnLen[GroupIndex] =		2;		//	Value is in %	//	double
		Zone.VectorLineThickness[GroupIndex] =		0.1;	//	Value is in %	//	double

		//	Scatter Settings

		Zone.ScatterShow[GroupIndex] =				TRUE;						//	Boolean_t
		Zone.ScatterShape[GroupIndex] =				GeomShape_Octahedron;			//	GeomShape_e
		Zone.ScatterOutColor[GroupIndex] =			MultiColor2_C;					//	ColorIndex_t
		Zone.ScatterInColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.ScatterFillMode[GroupIndex] =			FillMode_UseLineColor;		//	FillMode_e
		Zone.ScatterSize[GroupIndex] =				2.5;	//	value is in %	//	double
		Zone.ScatterLineThickness[GroupIndex] =		0.1;	//	value is in %	//	double

		//	Shade Settings

		Zone.ShadeShow[GroupIndex] =					FALSE;						//	Boolean_t
		Zone.ShadeColor[GroupIndex] =				Black_C;					//	ColorIndex_t
		Zone.ShadeLighting[GroupIndex] =				TRUE;						//	Boolean_t

		//	Edge Settings

		Zone.EdgeShow[GroupIndex] =					FALSE;						//	Boolean_t
		Zone.EdgeType[GroupIndex] =					EdgeType_BordersAndCreases;	//	EdgeType_e
		Zone.EdgeIBorder[GroupIndex] =				BorderLocation_Both;		//	BorderLocation_e
		Zone.EdgeJBorder[GroupIndex] =				BorderLocation_Both;		//	BorderLocation_e
		Zone.EdgeKBorder[GroupIndex] =				BorderLocation_Both;		//	BorderLocation_e
		Zone.EdgeColor[GroupIndex] =					Black_C;					//	ColorIndex_t
		Zone.EdgeLineThickness[GroupIndex] =			0.1;	//	value is in %	//	double

		//	Points Settings

		Zone.PointsToPlot[GroupIndex] =				PointsToPlot_SurfaceNodes;	//	PointsToPlot_e
		Zone.PointsIndexSkipI[GroupIndex] =			1;							//	LgIndex_t
		Zone.PointsIndexSkipJ[GroupIndex] =			1;							//	LgIndex_t
		Zone.PointsIndexSkipK[GroupIndex] =			1;							//	LgIndex_t
		/*
			1 is the same as "No Skip" in the GUI.
		*/

		//	Surfaces Settings

		Zone.SurfacesToPlot[GroupIndex] =			SurfacesToPlot_ExposedCellFaces;	//	SurfacesToPlot_e
		Zone.RangeBeginI[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeEndI[GroupIndex] =					0;							//	LgIndex_t
		Zone.RangeSkipI[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeBeginJ[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeEndJ[GroupIndex] =					0;							//	LgIndex_t
		Zone.RangeSkipJ[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeBeginK[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeEndK[GroupIndex] =					0;							//	LgIndex_t
		Zone.RangeSkipK[GroupIndex] =				1;							//	LgIndex_t
		/*
			0 is the same as maximum, which sets the end range to the maximum
			index for the dimension. The skip number is to skip every # planes,
			so 1 is to not skip any planes, (show every 1 plane).
		*/

		//	Volume Settings

		Zone.StreamtracesShow[GroupIndex] =			TRUE;						//	Boolean_t
		Zone.IsoSurfacesShow[GroupIndex] =			TRUE;						//	Boolean_t
		Zone.SlicesShow[GroupIndex] =				TRUE;						//	Boolean_t

		//	Effects Settings

		Zone.UseValueBlanking[GroupIndex] =			TRUE;						//	Boolean_t
		//Zone.UseClipPlanes[GroupIndex] =			TRUE;						//	Boolean_t
		Zone.UseSurfaceTranslucency[GroupIndex] =	FALSE;						//	Boolean_t
		Zone.SurfaceTranslucency[GroupIndex] =		50;		//	value is in %	//	SmInteger_t
		Zone.LightingEffect[GroupIndex] =			LightingEffect_Gouraud;		//	LightingEffect_e

		//	Group index and number (increments for next set)

		GroupIndex++;														//	int


		//===============================================================
		//	ABLines Settings
		//===============================================================

		Zone.ZoneShow[GroupIndex] =					TRUE;						//	Boolean_t

		//	Mesh Settings

		Zone.MeshShow[GroupIndex] =					TRUE;						//	Boolean_t
		Zone.MeshType[GroupIndex] =					MeshType_Overlay;			//	MeshType_e
		Zone.MeshColor[GroupIndex] =					MultiColor4_C;					//	ColorIndex_t
		/*
			A number of color index constants are \#defined. These include Black_C, Red_C,
			Green_C, Blue_C, Cyan_C, Yellow_C, Purple_C, White_C, Custom1_C through
			Custom56_C, MultiColor_C, NoColor_C, MulitiColor2_C through MulitiColor8_C,
			RGBColor_C, and InvalidColor_C.  
		*/
		Zone.MeshLinePttrn[GroupIndex] =				LinePattern_Solid;			//	LinePattern_e
		Zone.MeshLinePttrnLen[GroupIndex] =			2;		//	Value is in %	//	double
		Zone.MeshLineThickness[GroupIndex] =			0.3;		//	Value is in %	//	double


		//	Contour Settings

		Zone.ContourShow[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.ContourType[GroupIndex] =				ContourType_Flood;			//	ContourType_e
		Zone.FloodColoring[GroupIndex] =				ContourColoring_Group1;		//	ContourColoring_e
		Zone.ContourLineGroup[GroupIndex] =			1;							//	SmInteger_t
		Zone.ContourLineColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.ContourLinePttrn[GroupIndex] =			LinePattern_Solid;			//	LinePattern_e
		Zone.ContourLinePttrnLen[GroupIndex] =		2;		//	Value is in %	//	double
		Zone.ContourLineThickness[GroupIndex] =		0.1;	//	Value is in %	//	double
		Zone.ContourLighting[GroupIndex] =			TRUE;						//	Boolean_t

		//	Vector Settings
		
		Zone.VectorShow[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.VectorType[GroupIndex] =				VectorType_TailAtPoint;		//	VectorType_e
		Zone.VectorArrowType[GroupIndex] =			ArrowheadStyle_Plain;		//	ArrowheadStyle_e
		Zone.VectorLineColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.VectorTangent[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.VectorLinePttrn[GroupIndex] =			LinePattern_Solid;			//	LinePattern_e
		Zone.VectorLinePttrnLen[GroupIndex] =		2;		//	Value is in %	//	double
		Zone.VectorLineThickness[GroupIndex] =		0.1;	//	Value is in %	//	double

		//	Scatter Settings

		Zone.ScatterShow[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.ScatterShape[GroupIndex] =				GeomShape_Sphere;			//	GeomShape_e
		Zone.ScatterOutColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.ScatterInColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.ScatterFillMode[GroupIndex] =			FillMode_UseLineColor;		//	FillMode_e
		Zone.ScatterSize[GroupIndex] =				2.5;	//	value is in %	//	double
		Zone.ScatterLineThickness[GroupIndex] =		0.1;	//	value is in %	//	double

		//	Shade Settings

		Zone.ShadeShow[GroupIndex] =					FALSE;						//	Boolean_t
		Zone.ShadeColor[GroupIndex] =				Black_C;					//	ColorIndex_t
		Zone.ShadeLighting[GroupIndex] =				TRUE;						//	Boolean_t

		//	Edge Settings

		Zone.EdgeShow[GroupIndex] =					FALSE;						//	Boolean_t
		Zone.EdgeType[GroupIndex] =					EdgeType_BordersAndCreases;	//	EdgeType_e
		Zone.EdgeIBorder[GroupIndex] =				BorderLocation_Both;		//	BorderLocation_e
		Zone.EdgeJBorder[GroupIndex] =				BorderLocation_Both;		//	BorderLocation_e
		Zone.EdgeKBorder[GroupIndex] =				BorderLocation_Both;		//	BorderLocation_e
		Zone.EdgeColor[GroupIndex] =					Black_C;					//	ColorIndex_t
		Zone.EdgeLineThickness[GroupIndex] =			0.1;	//	value is in %	//	double

		//	Points Settings

		Zone.PointsToPlot[GroupIndex] =				PointsToPlot_SurfaceNodes;	//	PointsToPlot_e
		Zone.PointsIndexSkipI[GroupIndex] =			1;							//	LgIndex_t
		Zone.PointsIndexSkipJ[GroupIndex] =			1;							//	LgIndex_t
		Zone.PointsIndexSkipK[GroupIndex] =			1;							//	LgIndex_t
		/*
			1 is the same as "No Skip" in the GUI.
		*/

		//	Surfaces Settings

		Zone.SurfacesToPlot[GroupIndex] =			SurfacesToPlot_ExposedCellFaces;	//	SurfacesToPlot_e
		Zone.RangeBeginI[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeEndI[GroupIndex] =					0;							//	LgIndex_t
		Zone.RangeSkipI[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeBeginJ[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeEndJ[GroupIndex] =					0;							//	LgIndex_t
		Zone.RangeSkipJ[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeBeginK[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeEndK[GroupIndex] =					0;							//	LgIndex_t
		Zone.RangeSkipK[GroupIndex] =				1;							//	LgIndex_t
		/*
			0 is the same as maximum, which sets the end range to the maximum
			index for the dimension. The skip number is to skip every # planes,
			so 1 is to not skip any planes, (show every 1 plane).
		*/

		//	Volume Settings

		Zone.StreamtracesShow[GroupIndex] =			TRUE;						//	Boolean_t
		Zone.IsoSurfacesShow[GroupIndex] =			TRUE;						//	Boolean_t
		Zone.SlicesShow[GroupIndex] =				TRUE;						//	Boolean_t

		//	Effects Settings

		Zone.UseValueBlanking[GroupIndex] =			FALSE;						//	Boolean_t
		//Zone.UseClipPlanes[GroupIndex] =			TRUE;						//	Boolean_t
		Zone.UseSurfaceTranslucency[GroupIndex] =	FALSE;						//	Boolean_t
		Zone.SurfaceTranslucency[GroupIndex] =		50;		//	value is in %	//	SmInteger_t
		Zone.LightingEffect[GroupIndex] =			LightingEffect_Gouraud;		//	LightingEffect_e

		//	Group index and number (increments for next set)

		GroupIndex++;														//	int


		//===============================================================
		//	ABLinesFF Settings
		//===============================================================

		Zone.ZoneShow[GroupIndex] =					TRUE;						//	Boolean_t

		//	Mesh Settings

		Zone.MeshShow[GroupIndex] =					TRUE;						//	Boolean_t
		Zone.MeshType[GroupIndex] =					MeshType_Overlay;			//	MeshType_e
		Zone.MeshColor[GroupIndex] =					MultiColor4_C;					//	ColorIndex_t
		/*
			A number of color index constants are \#defined. These include Black_C, Red_C,
			Green_C, Blue_C, Cyan_C, Yellow_C, Purple_C, White_C, Custom1_C through
			Custom56_C, MultiColor_C, NoColor_C, MulitiColor2_C through MulitiColor8_C,
			RGBColor_C, and InvalidColor_C.  
		*/
		Zone.MeshLinePttrn[GroupIndex] =				LinePattern_Solid;			//	LinePattern_e
		Zone.MeshLinePttrnLen[GroupIndex] =			2;		//	Value is in %	//	double
		Zone.MeshLineThickness[GroupIndex] =			0.2;	//	Value is in %	//	double


		//	Contour Settings

		Zone.ContourShow[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.ContourType[GroupIndex] =				ContourType_Flood;			//	ContourType_e
		Zone.FloodColoring[GroupIndex] =				ContourColoring_Group1;		//	ContourColoring_e
		Zone.ContourLineGroup[GroupIndex] =			1;							//	SmInteger_t
		Zone.ContourLineColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.ContourLinePttrn[GroupIndex] =			LinePattern_Solid;			//	LinePattern_e
		Zone.ContourLinePttrnLen[GroupIndex] =		2;		//	Value is in %	//	double
		Zone.ContourLineThickness[GroupIndex] =		0.1;	//	Value is in %	//	double
		Zone.ContourLighting[GroupIndex] =			TRUE;						//	Boolean_t

		//	Vector Settings
		
		Zone.VectorShow[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.VectorType[GroupIndex] =				VectorType_TailAtPoint;		//	VectorType_e
		Zone.VectorArrowType[GroupIndex] =			ArrowheadStyle_Plain;		//	ArrowheadStyle_e
		Zone.VectorLineColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.VectorTangent[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.VectorLinePttrn[GroupIndex] =			LinePattern_Solid;			//	LinePattern_e
		Zone.VectorLinePttrnLen[GroupIndex] =		2;		//	Value is in %	//	double
		Zone.VectorLineThickness[GroupIndex] =		0.1;	//	Value is in %	//	double

		//	Scatter Settings

		Zone.ScatterShow[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.ScatterShape[GroupIndex] =				GeomShape_Sphere;			//	GeomShape_e
		Zone.ScatterOutColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.ScatterInColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.ScatterFillMode[GroupIndex] =			FillMode_UseLineColor;		//	FillMode_e
		Zone.ScatterSize[GroupIndex] =				2.5;	//	value is in %	//	double
		Zone.ScatterLineThickness[GroupIndex] =		0.1;	//	value is in %	//	double

		//	Shade Settings

		Zone.ShadeShow[GroupIndex] =					FALSE;						//	Boolean_t
		Zone.ShadeColor[GroupIndex] =				Black_C;					//	ColorIndex_t
		Zone.ShadeLighting[GroupIndex] =				TRUE;						//	Boolean_t

		//	Edge Settings

		Zone.EdgeShow[GroupIndex] =					FALSE;						//	Boolean_t
		Zone.EdgeType[GroupIndex] =					EdgeType_BordersAndCreases;	//	EdgeType_e
		Zone.EdgeIBorder[GroupIndex] =				BorderLocation_Both;		//	BorderLocation_e
		Zone.EdgeJBorder[GroupIndex] =				BorderLocation_Both;		//	BorderLocation_e
		Zone.EdgeKBorder[GroupIndex] =				BorderLocation_Both;		//	BorderLocation_e
		Zone.EdgeColor[GroupIndex] =					Black_C;					//	ColorIndex_t
		Zone.EdgeLineThickness[GroupIndex] =			0.1;	//	value is in %	//	double

		//	Points Settings

		Zone.PointsToPlot[GroupIndex] =				PointsToPlot_SurfaceNodes;	//	PointsToPlot_e
		Zone.PointsIndexSkipI[GroupIndex] =			1;							//	LgIndex_t
		Zone.PointsIndexSkipJ[GroupIndex] =			1;							//	LgIndex_t
		Zone.PointsIndexSkipK[GroupIndex] =			1;							//	LgIndex_t
		/*
			1 is the same as "No Skip" in the GUI.
		*/

		//	Surfaces Settings

		Zone.SurfacesToPlot[GroupIndex] =			SurfacesToPlot_ExposedCellFaces;	//	SurfacesToPlot_e
		Zone.RangeBeginI[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeEndI[GroupIndex] =					0;							//	LgIndex_t
		Zone.RangeSkipI[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeBeginJ[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeEndJ[GroupIndex] =					0;							//	LgIndex_t
		Zone.RangeSkipJ[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeBeginK[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeEndK[GroupIndex] =					0;							//	LgIndex_t
		Zone.RangeSkipK[GroupIndex] =				1;							//	LgIndex_t
		/*
			0 is the same as maximum, which sets the end range to the maximum
			index for the dimension. The skip number is to skip every # planes,
			so 1 is to not skip any planes, (show every 1 plane).
		*/

		//	Volume Settings

		Zone.StreamtracesShow[GroupIndex] =			TRUE;						//	Boolean_t
		Zone.IsoSurfacesShow[GroupIndex] =			TRUE;						//	Boolean_t
		Zone.SlicesShow[GroupIndex] =				TRUE;						//	Boolean_t

		//	Effects Settings

		Zone.UseValueBlanking[GroupIndex] =			TRUE;						//	Boolean_t
		//Zone.UseClipPlanes[GroupIndex] =			TRUE;						//	Boolean_t
		Zone.UseSurfaceTranslucency[GroupIndex] =	FALSE;						//	Boolean_t
		Zone.SurfaceTranslucency[GroupIndex] =		50;		//	value is in %	//	SmInteger_t
		Zone.LightingEffect[GroupIndex] =			LightingEffect_Gouraud;		//	LightingEffect_e

		//	Group index and number (increments for next set)

		GroupIndex++;														//	int


		//===============================================================
		//	ARLines Settings
		//===============================================================

		Zone.ZoneShow[GroupIndex] =					TRUE;						//	Boolean_t

		//	Mesh Settings

		Zone.MeshShow[GroupIndex] =					TRUE;						//	Boolean_t
		Zone.MeshType[GroupIndex] =					MeshType_Overlay;			//	MeshType_e
		Zone.MeshColor[GroupIndex] =					MultiColor5_C;					//	ColorIndex_t
		/*
			A number of color index constants are \#defined. These include Black_C, Red_C,
			Green_C, Blue_C, Cyan_C, Yellow_C, Purple_C, White_C, Custom1_C through
			Custom56_C, MultiColor_C, NoColor_C, MulitiColor2_C through MulitiColor8_C,
			RGBColor_C, and InvalidColor_C.  
		*/
		Zone.MeshLinePttrn[GroupIndex] =				LinePattern_LongDash;			//	LinePattern_e
		Zone.MeshLinePttrnLen[GroupIndex] =			2;		//	Value is in %	//	double
		Zone.MeshLineThickness[GroupIndex] =			0.2;	//	Value is in %	//	double


		//	Contour Settings

		Zone.ContourShow[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.ContourType[GroupIndex] =				ContourType_Flood;			//	ContourType_e
		Zone.FloodColoring[GroupIndex] =				ContourColoring_Group1;		//	ContourColoring_e
		Zone.ContourLineGroup[GroupIndex] =			1;							//	SmInteger_t
		Zone.ContourLineColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.ContourLinePttrn[GroupIndex] =			LinePattern_Solid;			//	LinePattern_e
		Zone.ContourLinePttrnLen[GroupIndex] =		2;		//	Value is in %	//	double
		Zone.ContourLineThickness[GroupIndex] =		0.1;	//	Value is in %	//	double
		Zone.ContourLighting[GroupIndex] =			TRUE;						//	Boolean_t

		//	Vector Settings
		
		Zone.VectorShow[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.VectorType[GroupIndex] =				VectorType_TailAtPoint;		//	VectorType_e
		Zone.VectorArrowType[GroupIndex] =			ArrowheadStyle_Plain;		//	ArrowheadStyle_e
		Zone.VectorLineColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.VectorTangent[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.VectorLinePttrn[GroupIndex] =			LinePattern_Solid;			//	LinePattern_e
		Zone.VectorLinePttrnLen[GroupIndex] =		2;		//	Value is in %	//	double
		Zone.VectorLineThickness[GroupIndex] =		0.1;	//	Value is in %	//	double

		//	Scatter Settings

		Zone.ScatterShow[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.ScatterShape[GroupIndex] =				GeomShape_Sphere;			//	GeomShape_e
		Zone.ScatterOutColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.ScatterInColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.ScatterFillMode[GroupIndex] =			FillMode_UseLineColor;		//	FillMode_e
		Zone.ScatterSize[GroupIndex] =				2.5;	//	value is in %	//	double
		Zone.ScatterLineThickness[GroupIndex] =		0.1;	//	value is in %	//	double

		//	Shade Settings

		Zone.ShadeShow[GroupIndex] =					FALSE;						//	Boolean_t
		Zone.ShadeColor[GroupIndex] =				Black_C;					//	ColorIndex_t
		Zone.ShadeLighting[GroupIndex] =				TRUE;						//	Boolean_t

		//	Edge Settings

		Zone.EdgeShow[GroupIndex] =					FALSE;						//	Boolean_t
		Zone.EdgeType[GroupIndex] =					EdgeType_BordersAndCreases;	//	EdgeType_e
		Zone.EdgeIBorder[GroupIndex] =				BorderLocation_Both;		//	BorderLocation_e
		Zone.EdgeJBorder[GroupIndex] =				BorderLocation_Both;		//	BorderLocation_e
		Zone.EdgeKBorder[GroupIndex] =				BorderLocation_Both;		//	BorderLocation_e
		Zone.EdgeColor[GroupIndex] =					Black_C;					//	ColorIndex_t
		Zone.EdgeLineThickness[GroupIndex] =			0.1;	//	value is in %	//	double

		//	Points Settings

		Zone.PointsToPlot[GroupIndex] =				PointsToPlot_SurfaceNodes;	//	PointsToPlot_e
		Zone.PointsIndexSkipI[GroupIndex] =			1;							//	LgIndex_t
		Zone.PointsIndexSkipJ[GroupIndex] =			1;							//	LgIndex_t
		Zone.PointsIndexSkipK[GroupIndex] =			1;							//	LgIndex_t
		/*
			1 is the same as "No Skip" in the GUI.
		*/

		//	Surfaces Settings

		Zone.SurfacesToPlot[GroupIndex] =			SurfacesToPlot_ExposedCellFaces;	//	SurfacesToPlot_e
		Zone.RangeBeginI[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeEndI[GroupIndex] =					0;							//	LgIndex_t
		Zone.RangeSkipI[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeBeginJ[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeEndJ[GroupIndex] =					0;							//	LgIndex_t
		Zone.RangeSkipJ[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeBeginK[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeEndK[GroupIndex] =					0;							//	LgIndex_t
		Zone.RangeSkipK[GroupIndex] =				1;							//	LgIndex_t
		/*
			0 is the same as maximum, which sets the end range to the maximum
			index for the dimension. The skip number is to skip every # planes,
			so 1 is to not skip any planes, (show every 1 plane).
		*/

		//	Volume Settings

		Zone.StreamtracesShow[GroupIndex] =			TRUE;						//	Boolean_t
		Zone.IsoSurfacesShow[GroupIndex] =			TRUE;						//	Boolean_t
		Zone.SlicesShow[GroupIndex] =				TRUE;						//	Boolean_t

		//	Effects Settings

		Zone.UseValueBlanking[GroupIndex] =			TRUE;						//	Boolean_t
		//Zone.UseClipPlanes[GroupIndex] =			TRUE;						//	Boolean_t
		Zone.UseSurfaceTranslucency[GroupIndex] =	FALSE;						//	Boolean_t
		Zone.SurfaceTranslucency[GroupIndex] =		50;		//	value is in %	//	SmInteger_t
		Zone.LightingEffect[GroupIndex] =			LightingEffect_Gouraud;		//	LightingEffect_e

		//	Group index and number (increments for next set)

		GroupIndex++;														//	int


		//===============================================================
		//	ARLinesFF Settings
		//===============================================================

		Zone.ZoneShow[GroupIndex] =					TRUE;						//	Boolean_t

		//	Mesh Settings

		Zone.MeshShow[GroupIndex] =					TRUE;						//	Boolean_t
		Zone.MeshType[GroupIndex] =					MeshType_Overlay;			//	MeshType_e
		Zone.MeshColor[GroupIndex] =					MultiColor5_C;					//	ColorIndex_t
		/*
			A number of color index constants are \#defined. These include Black_C, Red_C,
			Green_C, Blue_C, Cyan_C, Yellow_C, Purple_C, White_C, Custom1_C through
			Custom56_C, MultiColor_C, NoColor_C, MulitiColor2_C through MulitiColor8_C,
			RGBColor_C, and InvalidColor_C.  
		*/
		Zone.MeshLinePttrn[GroupIndex] =				LinePattern_LongDash;			//	LinePattern_e
		Zone.MeshLinePttrnLen[GroupIndex] =			2;		//	Value is in %	//	double
		Zone.MeshLineThickness[GroupIndex] =			0.2;	//	Value is in %	//	double


		//	Contour Settings

		Zone.ContourShow[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.ContourType[GroupIndex] =				ContourType_Flood;			//	ContourType_e
		Zone.FloodColoring[GroupIndex] =				ContourColoring_Group1;		//	ContourColoring_e
		Zone.ContourLineGroup[GroupIndex] =			1;							//	SmInteger_t
		Zone.ContourLineColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.ContourLinePttrn[GroupIndex] =			LinePattern_Solid;			//	LinePattern_e
		Zone.ContourLinePttrnLen[GroupIndex] =		2;		//	Value is in %	//	double
		Zone.ContourLineThickness[GroupIndex] =		0.1;	//	Value is in %	//	double
		Zone.ContourLighting[GroupIndex] =			TRUE;						//	Boolean_t

		//	Vector Settings
		
		Zone.VectorShow[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.VectorType[GroupIndex] =				VectorType_TailAtPoint;		//	VectorType_e
		Zone.VectorArrowType[GroupIndex] =			ArrowheadStyle_Plain;		//	ArrowheadStyle_e
		Zone.VectorLineColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.VectorTangent[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.VectorLinePttrn[GroupIndex] =			LinePattern_Solid;			//	LinePattern_e
		Zone.VectorLinePttrnLen[GroupIndex] =		2;		//	Value is in %	//	double
		Zone.VectorLineThickness[GroupIndex] =		0.1;	//	Value is in %	//	double

		//	Scatter Settings

		Zone.ScatterShow[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.ScatterShape[GroupIndex] =				GeomShape_Sphere;			//	GeomShape_e
		Zone.ScatterOutColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.ScatterInColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.ScatterFillMode[GroupIndex] =			FillMode_UseLineColor;		//	FillMode_e
		Zone.ScatterSize[GroupIndex] =				2.5;	//	value is in %	//	double
		Zone.ScatterLineThickness[GroupIndex] =		0.1;	//	value is in %	//	double

		//	Shade Settings

		Zone.ShadeShow[GroupIndex] =					FALSE;						//	Boolean_t
		Zone.ShadeColor[GroupIndex] =				Black_C;					//	ColorIndex_t
		Zone.ShadeLighting[GroupIndex] =				TRUE;						//	Boolean_t

		//	Edge Settings

		Zone.EdgeShow[GroupIndex] =					FALSE;						//	Boolean_t
		Zone.EdgeType[GroupIndex] =					EdgeType_BordersAndCreases;	//	EdgeType_e
		Zone.EdgeIBorder[GroupIndex] =				BorderLocation_Both;		//	BorderLocation_e
		Zone.EdgeJBorder[GroupIndex] =				BorderLocation_Both;		//	BorderLocation_e
		Zone.EdgeKBorder[GroupIndex] =				BorderLocation_Both;		//	BorderLocation_e
		Zone.EdgeColor[GroupIndex] =					Black_C;					//	ColorIndex_t
		Zone.EdgeLineThickness[GroupIndex] =			0.1;	//	value is in %	//	double

		//	Points Settings

		Zone.PointsToPlot[GroupIndex] =				PointsToPlot_SurfaceNodes;	//	PointsToPlot_e
		Zone.PointsIndexSkipI[GroupIndex] =			1;							//	LgIndex_t
		Zone.PointsIndexSkipJ[GroupIndex] =			1;							//	LgIndex_t
		Zone.PointsIndexSkipK[GroupIndex] =			1;							//	LgIndex_t
		/*
			1 is the same as "No Skip" in the GUI.
		*/

		//	Surfaces Settings

		Zone.SurfacesToPlot[GroupIndex] =			SurfacesToPlot_ExposedCellFaces;	//	SurfacesToPlot_e
		Zone.RangeBeginI[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeEndI[GroupIndex] =					0;							//	LgIndex_t
		Zone.RangeSkipI[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeBeginJ[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeEndJ[GroupIndex] =					0;							//	LgIndex_t
		Zone.RangeSkipJ[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeBeginK[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeEndK[GroupIndex] =					0;							//	LgIndex_t
		Zone.RangeSkipK[GroupIndex] =				1;							//	LgIndex_t
		/*
			0 is the same as maximum, which sets the end range to the maximum
			index for the dimension. The skip number is to skip every # planes,
			so 1 is to not skip any planes, (show every 1 plane).
		*/

		//	Volume Settings

		Zone.StreamtracesShow[GroupIndex] =			TRUE;						//	Boolean_t
		Zone.IsoSurfacesShow[GroupIndex] =			TRUE;						//	Boolean_t
		Zone.SlicesShow[GroupIndex] =				TRUE;						//	Boolean_t

		//	Effects Settings

		Zone.UseValueBlanking[GroupIndex] =			TRUE;						//	Boolean_t
		//Zone.UseClipPlanes[GroupIndex] =			TRUE;						//	Boolean_t
		Zone.UseSurfaceTranslucency[GroupIndex] =	FALSE;						//	Boolean_t
		Zone.SurfaceTranslucency[GroupIndex] =		50;		//	value is in %	//	SmInteger_t
		Zone.LightingEffect[GroupIndex] =			LightingEffect_Gouraud;		//	LightingEffect_e

		//	Group index and number (increments for next set)

		GroupIndex++;														//	int


		//===============================================================
		//	ACLines Settings
		//===============================================================

		Zone.ZoneShow[GroupIndex] =					TRUE;						//	Boolean_t

		//	Mesh Settings

		Zone.MeshShow[GroupIndex] =					TRUE;						//	Boolean_t
		Zone.MeshType[GroupIndex] =					MeshType_Overlay;			//	MeshType_e
		Zone.MeshColor[GroupIndex] =					MultiColor6_C;					//	ColorIndex_t
		/*
			A number of color index constants are \#defined. These include Black_C, Red_C,
			Green_C, Blue_C, Cyan_C, Yellow_C, Purple_C, White_C, Custom1_C through
			Custom56_C, MultiColor_C, NoColor_C, MulitiColor2_C through MulitiColor8_C,
			RGBColor_C, and InvalidColor_C.  
		*/
		Zone.MeshLinePttrn[GroupIndex] =				LinePattern_Dashed;			//	LinePattern_e
		Zone.MeshLinePttrnLen[GroupIndex] =			2;		//	Value is in %	//	double
		Zone.MeshLineThickness[GroupIndex] =			0.2;	//	Value is in %	//	double


		//	Contour Settings

		Zone.ContourShow[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.ContourType[GroupIndex] =				ContourType_Flood;			//	ContourType_e
		Zone.FloodColoring[GroupIndex] =				ContourColoring_Group1;		//	ContourColoring_e
		Zone.ContourLineGroup[GroupIndex] =			1;							//	SmInteger_t
		Zone.ContourLineColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.ContourLinePttrn[GroupIndex] =			LinePattern_Solid;			//	LinePattern_e
		Zone.ContourLinePttrnLen[GroupIndex] =		2;		//	Value is in %	//	double
		Zone.ContourLineThickness[GroupIndex] =		0.1;	//	Value is in %	//	double
		Zone.ContourLighting[GroupIndex] =			TRUE;						//	Boolean_t

		//	Vector Settings
		
		Zone.VectorShow[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.VectorType[GroupIndex] =				VectorType_TailAtPoint;		//	VectorType_e
		Zone.VectorArrowType[GroupIndex] =			ArrowheadStyle_Plain;		//	ArrowheadStyle_e
		Zone.VectorLineColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.VectorTangent[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.VectorLinePttrn[GroupIndex] =			LinePattern_Solid;			//	LinePattern_e
		Zone.VectorLinePttrnLen[GroupIndex] =		2;		//	Value is in %	//	double
		Zone.VectorLineThickness[GroupIndex] =		0.1;	//	Value is in %	//	double

		//	Scatter Settings

		Zone.ScatterShow[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.ScatterShape[GroupIndex] =				GeomShape_Sphere;			//	GeomShape_e
		Zone.ScatterOutColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.ScatterInColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.ScatterFillMode[GroupIndex] =			FillMode_UseLineColor;		//	FillMode_e
		Zone.ScatterSize[GroupIndex] =				2.5;	//	value is in %	//	double
		Zone.ScatterLineThickness[GroupIndex] =		0.1;	//	value is in %	//	double

		//	Shade Settings

		Zone.ShadeShow[GroupIndex] =					FALSE;						//	Boolean_t
		Zone.ShadeColor[GroupIndex] =				Black_C;					//	ColorIndex_t
		Zone.ShadeLighting[GroupIndex] =				TRUE;						//	Boolean_t

		//	Edge Settings

		Zone.EdgeShow[GroupIndex] =					FALSE;						//	Boolean_t
		Zone.EdgeType[GroupIndex] =					EdgeType_BordersAndCreases;	//	EdgeType_e
		Zone.EdgeIBorder[GroupIndex] =				BorderLocation_Both;		//	BorderLocation_e
		Zone.EdgeJBorder[GroupIndex] =				BorderLocation_Both;		//	BorderLocation_e
		Zone.EdgeKBorder[GroupIndex] =				BorderLocation_Both;		//	BorderLocation_e
		Zone.EdgeColor[GroupIndex] =					Black_C;					//	ColorIndex_t
		Zone.EdgeLineThickness[GroupIndex] =			0.1;	//	value is in %	//	double

		//	Points Settings

		Zone.PointsToPlot[GroupIndex] =				PointsToPlot_SurfaceNodes;	//	PointsToPlot_e
		Zone.PointsIndexSkipI[GroupIndex] =			1;							//	LgIndex_t
		Zone.PointsIndexSkipJ[GroupIndex] =			1;							//	LgIndex_t
		Zone.PointsIndexSkipK[GroupIndex] =			1;							//	LgIndex_t
		/*
			1 is the same as "No Skip" in the GUI.
		*/

		//	Surfaces Settings

		Zone.SurfacesToPlot[GroupIndex] =			SurfacesToPlot_ExposedCellFaces;	//	SurfacesToPlot_e
		Zone.RangeBeginI[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeEndI[GroupIndex] =					0;							//	LgIndex_t
		Zone.RangeSkipI[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeBeginJ[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeEndJ[GroupIndex] =					0;							//	LgIndex_t
		Zone.RangeSkipJ[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeBeginK[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeEndK[GroupIndex] =					0;							//	LgIndex_t
		Zone.RangeSkipK[GroupIndex] =				1;							//	LgIndex_t
		/*
			0 is the same as maximum, which sets the end range to the maximum
			index for the dimension. The skip number is to skip every # planes,
			so 1 is to not skip any planes, (show every 1 plane).
		*/

		//	Volume Settings

		Zone.StreamtracesShow[GroupIndex] =			TRUE;						//	Boolean_t
		Zone.IsoSurfacesShow[GroupIndex] =			TRUE;						//	Boolean_t
		Zone.SlicesShow[GroupIndex] =				TRUE;						//	Boolean_t

		//	Effects Settings

		Zone.UseValueBlanking[GroupIndex] =			TRUE;						//	Boolean_t
		//Zone.UseClipPlanes[GroupIndex] =			TRUE;						//	Boolean_t
		Zone.UseSurfaceTranslucency[GroupIndex] =	FALSE;						//	Boolean_t
		Zone.SurfaceTranslucency[GroupIndex] =		50;		//	value is in %	//	SmInteger_t
		Zone.LightingEffect[GroupIndex] =			LightingEffect_Gouraud;		//	LightingEffect_e

		//	Group index and number (increments for next set)

		GroupIndex++;														//	int


		//===============================================================
		//	ACLinesFF Settings
		//===============================================================

		Zone.ZoneShow[GroupIndex] =					TRUE;						//	Boolean_t

		//	Mesh Settings

		Zone.MeshShow[GroupIndex] =					TRUE;						//	Boolean_t
		Zone.MeshType[GroupIndex] =					MeshType_Overlay;			//	MeshType_e
		Zone.MeshColor[GroupIndex] =					MultiColor6_C;					//	ColorIndex_t
		/*
			A number of color index constants are \#defined. These include Black_C, Red_C,
			Green_C, Blue_C, Cyan_C, Yellow_C, Purple_C, White_C, Custom1_C through
			Custom56_C, MultiColor_C, NoColor_C, MulitiColor2_C through MulitiColor8_C,
			RGBColor_C, and InvalidColor_C.  
		*/
		Zone.MeshLinePttrn[GroupIndex] =				LinePattern_Dashed;			//	LinePattern_e
		Zone.MeshLinePttrnLen[GroupIndex] =			2;		//	Value is in %	//	double
		Zone.MeshLineThickness[GroupIndex] =			0.2;	//	Value is in %	//	double


		//	Contour Settings

		Zone.ContourShow[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.ContourType[GroupIndex] =				ContourType_Flood;			//	ContourType_e
		Zone.FloodColoring[GroupIndex] =				ContourColoring_Group1;		//	ContourColoring_e
		Zone.ContourLineGroup[GroupIndex] =			1;							//	SmInteger_t
		Zone.ContourLineColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.ContourLinePttrn[GroupIndex] =			LinePattern_Solid;			//	LinePattern_e
		Zone.ContourLinePttrnLen[GroupIndex] =		2;		//	Value is in %	//	double
		Zone.ContourLineThickness[GroupIndex] =		0.1;	//	Value is in %	//	double
		Zone.ContourLighting[GroupIndex] =			TRUE;						//	Boolean_t

		//	Vector Settings
		
		Zone.VectorShow[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.VectorType[GroupIndex] =				VectorType_TailAtPoint;		//	VectorType_e
		Zone.VectorArrowType[GroupIndex] =			ArrowheadStyle_Plain;		//	ArrowheadStyle_e
		Zone.VectorLineColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.VectorTangent[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.VectorLinePttrn[GroupIndex] =			LinePattern_Solid;			//	LinePattern_e
		Zone.VectorLinePttrnLen[GroupIndex] =		2;		//	Value is in %	//	double
		Zone.VectorLineThickness[GroupIndex] =		0.1;	//	Value is in %	//	double

		//	Scatter Settings

		Zone.ScatterShow[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.ScatterShape[GroupIndex] =				GeomShape_Sphere;			//	GeomShape_e
		Zone.ScatterOutColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.ScatterInColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.ScatterFillMode[GroupIndex] =			FillMode_UseLineColor;		//	FillMode_e
		Zone.ScatterSize[GroupIndex] =				2.5;	//	value is in %	//	double
		Zone.ScatterLineThickness[GroupIndex] =		0.1;	//	value is in %	//	double

		//	Shade Settings

		Zone.ShadeShow[GroupIndex] =					FALSE;						//	Boolean_t
		Zone.ShadeColor[GroupIndex] =				Black_C;					//	ColorIndex_t
		Zone.ShadeLighting[GroupIndex] =				TRUE;						//	Boolean_t

		//	Edge Settings

		Zone.EdgeShow[GroupIndex] =					FALSE;						//	Boolean_t
		Zone.EdgeType[GroupIndex] =					EdgeType_BordersAndCreases;	//	EdgeType_e
		Zone.EdgeIBorder[GroupIndex] =				BorderLocation_Both;		//	BorderLocation_e
		Zone.EdgeJBorder[GroupIndex] =				BorderLocation_Both;		//	BorderLocation_e
		Zone.EdgeKBorder[GroupIndex] =				BorderLocation_Both;		//	BorderLocation_e
		Zone.EdgeColor[GroupIndex] =					Black_C;					//	ColorIndex_t
		Zone.EdgeLineThickness[GroupIndex] =			0.1;	//	value is in %	//	double

		//	Points Settings

		Zone.PointsToPlot[GroupIndex] =				PointsToPlot_SurfaceNodes;	//	PointsToPlot_e
		Zone.PointsIndexSkipI[GroupIndex] =			1;							//	LgIndex_t
		Zone.PointsIndexSkipJ[GroupIndex] =			1;							//	LgIndex_t
		Zone.PointsIndexSkipK[GroupIndex] =			1;							//	LgIndex_t
		/*
			1 is the same as "No Skip" in the GUI.
		*/

		//	Surfaces Settings

		Zone.SurfacesToPlot[GroupIndex] =			SurfacesToPlot_ExposedCellFaces;	//	SurfacesToPlot_e
		Zone.RangeBeginI[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeEndI[GroupIndex] =					0;							//	LgIndex_t
		Zone.RangeSkipI[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeBeginJ[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeEndJ[GroupIndex] =					0;							//	LgIndex_t
		Zone.RangeSkipJ[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeBeginK[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeEndK[GroupIndex] =					0;							//	LgIndex_t
		Zone.RangeSkipK[GroupIndex] =				1;							//	LgIndex_t
		/*
			0 is the same as maximum, which sets the end range to the maximum
			index for the dimension. The skip number is to skip every # planes,
			so 1 is to not skip any planes, (show every 1 plane).
		*/

		//	Volume Settings

		Zone.StreamtracesShow[GroupIndex] =			TRUE;						//	Boolean_t
		Zone.IsoSurfacesShow[GroupIndex] =			TRUE;						//	Boolean_t
		Zone.SlicesShow[GroupIndex] =				TRUE;						//	Boolean_t

		//	Effects Settings

		Zone.UseValueBlanking[GroupIndex] =			TRUE;						//	Boolean_t
		//Zone.UseClipPlanes[GroupIndex] =			TRUE;						//	Boolean_t
		Zone.UseSurfaceTranslucency[GroupIndex] =	FALSE;						//	Boolean_t
		Zone.SurfaceTranslucency[GroupIndex] =		50;		//	value is in %	//	SmInteger_t
		Zone.LightingEffect[GroupIndex] =			LightingEffect_Gouraud;		//	LightingEffect_e

		//	Group index and number (increments for next set)

		GroupIndex++;														//	int


		//===============================================================
		//	BRLines Settings
		//===============================================================

		Zone.ZoneShow[GroupIndex] =					TRUE;						//	Boolean_t

		//	Mesh Settings

		Zone.MeshShow[GroupIndex] =					TRUE;						//	Boolean_t
		Zone.MeshType[GroupIndex] =					MeshType_Overlay;			//	MeshType_e
		Zone.MeshColor[GroupIndex] =					MultiColor7_C;					//	ColorIndex_t
		/*
			A number of color index constants are \#defined. These include Black_C, Red_C,
			Green_C, Blue_C, Cyan_C, Yellow_C, Purple_C, White_C, Custom1_C through
			Custom56_C, MultiColor_C, NoColor_C, MulitiColor2_C through MulitiColor8_C,
			RGBColor_C, and InvalidColor_C.  
		*/
		Zone.MeshLinePttrn[GroupIndex] =				LinePattern_DashDot;			//	LinePattern_e
		Zone.MeshLinePttrnLen[GroupIndex] =			2;		//	Value is in %	//	double
		Zone.MeshLineThickness[GroupIndex] =			0.2;	//	Value is in %	//	double


		//	Contour Settings

		Zone.ContourShow[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.ContourType[GroupIndex] =				ContourType_Flood;			//	ContourType_e
		Zone.FloodColoring[GroupIndex] =				ContourColoring_Group1;		//	ContourColoring_e
		Zone.ContourLineGroup[GroupIndex] =			1;							//	SmInteger_t
		Zone.ContourLineColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.ContourLinePttrn[GroupIndex] =			LinePattern_Solid;			//	LinePattern_e
		Zone.ContourLinePttrnLen[GroupIndex] =		2;		//	Value is in %	//	double
		Zone.ContourLineThickness[GroupIndex] =		0.1;	//	Value is in %	//	double
		Zone.ContourLighting[GroupIndex] =			TRUE;						//	Boolean_t

		//	Vector Settings
		
		Zone.VectorShow[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.VectorType[GroupIndex] =				VectorType_TailAtPoint;		//	VectorType_e
		Zone.VectorArrowType[GroupIndex] =			ArrowheadStyle_Plain;		//	ArrowheadStyle_e
		Zone.VectorLineColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.VectorTangent[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.VectorLinePttrn[GroupIndex] =			LinePattern_Solid;			//	LinePattern_e
		Zone.VectorLinePttrnLen[GroupIndex] =		2;		//	Value is in %	//	double
		Zone.VectorLineThickness[GroupIndex] =		0.1;	//	Value is in %	//	double

		//	Scatter Settings

		Zone.ScatterShow[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.ScatterShape[GroupIndex] =				GeomShape_Sphere;			//	GeomShape_e
		Zone.ScatterOutColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.ScatterInColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.ScatterFillMode[GroupIndex] =			FillMode_UseLineColor;		//	FillMode_e
		Zone.ScatterSize[GroupIndex] =				2.5;	//	value is in %	//	double
		Zone.ScatterLineThickness[GroupIndex] =		0.1;	//	value is in %	//	double

		//	Shade Settings

		Zone.ShadeShow[GroupIndex] =					FALSE;						//	Boolean_t
		Zone.ShadeColor[GroupIndex] =				Black_C;					//	ColorIndex_t
		Zone.ShadeLighting[GroupIndex] =				TRUE;						//	Boolean_t

		//	Edge Settings

		Zone.EdgeShow[GroupIndex] =					FALSE;						//	Boolean_t
		Zone.EdgeType[GroupIndex] =					EdgeType_BordersAndCreases;	//	EdgeType_e
		Zone.EdgeIBorder[GroupIndex] =				BorderLocation_Both;		//	BorderLocation_e
		Zone.EdgeJBorder[GroupIndex] =				BorderLocation_Both;		//	BorderLocation_e
		Zone.EdgeKBorder[GroupIndex] =				BorderLocation_Both;		//	BorderLocation_e
		Zone.EdgeColor[GroupIndex] =					Black_C;					//	ColorIndex_t
		Zone.EdgeLineThickness[GroupIndex] =			0.1;	//	value is in %	//	double

		//	Points Settings

		Zone.PointsToPlot[GroupIndex] =				PointsToPlot_SurfaceNodes;	//	PointsToPlot_e
		Zone.PointsIndexSkipI[GroupIndex] =			1;							//	LgIndex_t
		Zone.PointsIndexSkipJ[GroupIndex] =			1;							//	LgIndex_t
		Zone.PointsIndexSkipK[GroupIndex] =			1;							//	LgIndex_t
		/*
			1 is the same as "No Skip" in the GUI.
		*/

		//	Surfaces Settings

		Zone.SurfacesToPlot[GroupIndex] =			SurfacesToPlot_ExposedCellFaces;	//	SurfacesToPlot_e
		Zone.RangeBeginI[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeEndI[GroupIndex] =					0;							//	LgIndex_t
		Zone.RangeSkipI[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeBeginJ[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeEndJ[GroupIndex] =					0;							//	LgIndex_t
		Zone.RangeSkipJ[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeBeginK[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeEndK[GroupIndex] =					0;							//	LgIndex_t
		Zone.RangeSkipK[GroupIndex] =				1;							//	LgIndex_t
		/*
			0 is the same as maximum, which sets the end range to the maximum
			index for the dimension. The skip number is to skip every # planes,
			so 1 is to not skip any planes, (show every 1 plane).
		*/

		//	Volume Settings

		Zone.StreamtracesShow[GroupIndex] =			TRUE;						//	Boolean_t
		Zone.IsoSurfacesShow[GroupIndex] =			TRUE;						//	Boolean_t
		Zone.SlicesShow[GroupIndex] =				TRUE;						//	Boolean_t

		//	Effects Settings

		Zone.UseValueBlanking[GroupIndex] =			TRUE;						//	Boolean_t
		//Zone.UseClipPlanes[GroupIndex] =			TRUE;						//	Boolean_t
		Zone.UseSurfaceTranslucency[GroupIndex] =	FALSE;						//	Boolean_t
		Zone.SurfaceTranslucency[GroupIndex] =		50;		//	value is in %	//	SmInteger_t
		Zone.LightingEffect[GroupIndex] =			LightingEffect_Gouraud;		//	LightingEffect_e

		//	Group index and number (increments for next set)

		GroupIndex++;														//	int


		//===============================================================
		//	BRLinesFF Settings
		//===============================================================

		Zone.ZoneShow[GroupIndex] =					TRUE;						//	Boolean_t

		//	Mesh Settings

		Zone.MeshShow[GroupIndex] =					TRUE;						//	Boolean_t
		Zone.MeshType[GroupIndex] =					MeshType_Overlay;			//	MeshType_e
		Zone.MeshColor[GroupIndex] =					MultiColor7_C;					//	ColorIndex_t
		/*
			A number of color index constants are \#defined. These include Black_C, Red_C,
			Green_C, Blue_C, Cyan_C, Yellow_C, Purple_C, White_C, Custom1_C through
			Custom56_C, MultiColor_C, NoColor_C, MulitiColor2_C through MulitiColor8_C,
			RGBColor_C, and InvalidColor_C.  
		*/
		Zone.MeshLinePttrn[GroupIndex] =				LinePattern_DashDot;			//	LinePattern_e
		Zone.MeshLinePttrnLen[GroupIndex] =			2;		//	Value is in %	//	double
		Zone.MeshLineThickness[GroupIndex] =			0.2;	//	Value is in %	//	double


		//	Contour Settings

		Zone.ContourShow[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.ContourType[GroupIndex] =				ContourType_Flood;			//	ContourType_e
		Zone.FloodColoring[GroupIndex] =				ContourColoring_Group1;		//	ContourColoring_e
		Zone.ContourLineGroup[GroupIndex] =			1;							//	SmInteger_t
		Zone.ContourLineColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.ContourLinePttrn[GroupIndex] =			LinePattern_Solid;			//	LinePattern_e
		Zone.ContourLinePttrnLen[GroupIndex] =		2;		//	Value is in %	//	double
		Zone.ContourLineThickness[GroupIndex] =		0.1;	//	Value is in %	//	double
		Zone.ContourLighting[GroupIndex] =			TRUE;						//	Boolean_t

		//	Vector Settings
		
		Zone.VectorShow[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.VectorType[GroupIndex] =				VectorType_TailAtPoint;		//	VectorType_e
		Zone.VectorArrowType[GroupIndex] =			ArrowheadStyle_Plain;		//	ArrowheadStyle_e
		Zone.VectorLineColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.VectorTangent[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.VectorLinePttrn[GroupIndex] =			LinePattern_Solid;			//	LinePattern_e
		Zone.VectorLinePttrnLen[GroupIndex] =		2;		//	Value is in %	//	double
		Zone.VectorLineThickness[GroupIndex] =		0.1;	//	Value is in %	//	double

		//	Scatter Settings

		Zone.ScatterShow[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.ScatterShape[GroupIndex] =				GeomShape_Sphere;			//	GeomShape_e
		Zone.ScatterOutColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.ScatterInColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.ScatterFillMode[GroupIndex] =			FillMode_UseLineColor;		//	FillMode_e
		Zone.ScatterSize[GroupIndex] =				2.5;	//	value is in %	//	double
		Zone.ScatterLineThickness[GroupIndex] =		0.1;	//	value is in %	//	double

		//	Shade Settings

		Zone.ShadeShow[GroupIndex] =					FALSE;						//	Boolean_t
		Zone.ShadeColor[GroupIndex] =				Black_C;					//	ColorIndex_t
		Zone.ShadeLighting[GroupIndex] =				TRUE;						//	Boolean_t

		//	Edge Settings

		Zone.EdgeShow[GroupIndex] =					FALSE;						//	Boolean_t
		Zone.EdgeType[GroupIndex] =					EdgeType_BordersAndCreases;	//	EdgeType_e
		Zone.EdgeIBorder[GroupIndex] =				BorderLocation_Both;		//	BorderLocation_e
		Zone.EdgeJBorder[GroupIndex] =				BorderLocation_Both;		//	BorderLocation_e
		Zone.EdgeKBorder[GroupIndex] =				BorderLocation_Both;		//	BorderLocation_e
		Zone.EdgeColor[GroupIndex] =					Black_C;					//	ColorIndex_t
		Zone.EdgeLineThickness[GroupIndex] =			0.1;	//	value is in %	//	double

		//	Points Settings

		Zone.PointsToPlot[GroupIndex] =				PointsToPlot_SurfaceNodes;	//	PointsToPlot_e
		Zone.PointsIndexSkipI[GroupIndex] =			1;							//	LgIndex_t
		Zone.PointsIndexSkipJ[GroupIndex] =			1;							//	LgIndex_t
		Zone.PointsIndexSkipK[GroupIndex] =			1;							//	LgIndex_t
		/*
			1 is the same as "No Skip" in the GUI.
		*/

		//	Surfaces Settings

		Zone.SurfacesToPlot[GroupIndex] =			SurfacesToPlot_ExposedCellFaces;	//	SurfacesToPlot_e
		Zone.RangeBeginI[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeEndI[GroupIndex] =					0;							//	LgIndex_t
		Zone.RangeSkipI[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeBeginJ[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeEndJ[GroupIndex] =					0;							//	LgIndex_t
		Zone.RangeSkipJ[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeBeginK[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeEndK[GroupIndex] =					0;							//	LgIndex_t
		Zone.RangeSkipK[GroupIndex] =				1;							//	LgIndex_t
		/*
			0 is the same as maximum, which sets the end range to the maximum
			index for the dimension. The skip number is to skip every # planes,
			so 1 is to not skip any planes, (show every 1 plane).
		*/

		//	Volume Settings

		Zone.StreamtracesShow[GroupIndex] =			TRUE;						//	Boolean_t
		Zone.IsoSurfacesShow[GroupIndex] =			TRUE;						//	Boolean_t
		Zone.SlicesShow[GroupIndex] =				TRUE;						//	Boolean_t

		//	Effects Settings

		Zone.UseValueBlanking[GroupIndex] =			TRUE;						//	Boolean_t
		//Zone.UseClipPlanes[GroupIndex] =			TRUE;						//	Boolean_t
		Zone.UseSurfaceTranslucency[GroupIndex] =	FALSE;						//	Boolean_t
		Zone.SurfaceTranslucency[GroupIndex] =		50;		//	value is in %	//	SmInteger_t
		Zone.LightingEffect[GroupIndex] =			LightingEffect_Gouraud;		//	LightingEffect_e

		//	Group index and number (increments for next set)

		GroupIndex++;														//	int


		//===============================================================
		//	BCLines Settings
		//===============================================================

		Zone.ZoneShow[GroupIndex] =					TRUE;						//	Boolean_t

		//	Mesh Settings

		Zone.MeshShow[GroupIndex] =					TRUE;						//	Boolean_t
		Zone.MeshType[GroupIndex] =					MeshType_Overlay;			//	MeshType_e
		Zone.MeshColor[GroupIndex] =					MultiColor8_C;					//	ColorIndex_t
		/*
			A number of color index constants are \#defined. These include Black_C, Red_C,
			Green_C, Blue_C, Cyan_C, Yellow_C, Purple_C, White_C, Custom1_C through
			Custom56_C, MultiColor_C, NoColor_C, MulitiColor2_C through MulitiColor8_C,
			RGBColor_C, and InvalidColor_C.  
		*/
		Zone.MeshLinePttrn[GroupIndex] =				LinePattern_DashDotDot;			//	LinePattern_e
		Zone.MeshLinePttrnLen[GroupIndex] =			2;		//	Value is in %	//	double
		Zone.MeshLineThickness[GroupIndex] =			0.2;	//	Value is in %	//	double


		//	Contour Settings

		Zone.ContourShow[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.ContourType[GroupIndex] =				ContourType_Flood;			//	ContourType_e
		Zone.FloodColoring[GroupIndex] =				ContourColoring_Group1;		//	ContourColoring_e
		Zone.ContourLineGroup[GroupIndex] =			1;							//	SmInteger_t
		Zone.ContourLineColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.ContourLinePttrn[GroupIndex] =			LinePattern_Solid;			//	LinePattern_e
		Zone.ContourLinePttrnLen[GroupIndex] =		2;		//	Value is in %	//	double
		Zone.ContourLineThickness[GroupIndex] =		0.1;	//	Value is in %	//	double
		Zone.ContourLighting[GroupIndex] =			TRUE;						//	Boolean_t

		//	Vector Settings
		
		Zone.VectorShow[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.VectorType[GroupIndex] =				VectorType_TailAtPoint;		//	VectorType_e
		Zone.VectorArrowType[GroupIndex] =			ArrowheadStyle_Plain;		//	ArrowheadStyle_e
		Zone.VectorLineColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.VectorTangent[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.VectorLinePttrn[GroupIndex] =			LinePattern_Solid;			//	LinePattern_e
		Zone.VectorLinePttrnLen[GroupIndex] =		2;		//	Value is in %	//	double
		Zone.VectorLineThickness[GroupIndex] =		0.1;	//	Value is in %	//	double

		//	Scatter Settings

		Zone.ScatterShow[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.ScatterShape[GroupIndex] =				GeomShape_Sphere;			//	GeomShape_e
		Zone.ScatterOutColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.ScatterInColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.ScatterFillMode[GroupIndex] =			FillMode_UseLineColor;		//	FillMode_e
		Zone.ScatterSize[GroupIndex] =				2.5;	//	value is in %	//	double
		Zone.ScatterLineThickness[GroupIndex] =		0.1;	//	value is in %	//	double

		//	Shade Settings

		Zone.ShadeShow[GroupIndex] =					FALSE;						//	Boolean_t
		Zone.ShadeColor[GroupIndex] =				Black_C;					//	ColorIndex_t
		Zone.ShadeLighting[GroupIndex] =				TRUE;						//	Boolean_t

		//	Edge Settings

		Zone.EdgeShow[GroupIndex] =					FALSE;						//	Boolean_t
		Zone.EdgeType[GroupIndex] =					EdgeType_BordersAndCreases;	//	EdgeType_e
		Zone.EdgeIBorder[GroupIndex] =				BorderLocation_Both;		//	BorderLocation_e
		Zone.EdgeJBorder[GroupIndex] =				BorderLocation_Both;		//	BorderLocation_e
		Zone.EdgeKBorder[GroupIndex] =				BorderLocation_Both;		//	BorderLocation_e
		Zone.EdgeColor[GroupIndex] =					Black_C;					//	ColorIndex_t
		Zone.EdgeLineThickness[GroupIndex] =			0.1;	//	value is in %	//	double

		//	Points Settings

		Zone.PointsToPlot[GroupIndex] =				PointsToPlot_SurfaceNodes;	//	PointsToPlot_e
		Zone.PointsIndexSkipI[GroupIndex] =			1;							//	LgIndex_t
		Zone.PointsIndexSkipJ[GroupIndex] =			1;							//	LgIndex_t
		Zone.PointsIndexSkipK[GroupIndex] =			1;							//	LgIndex_t
		/*
			1 is the same as "No Skip" in the GUI.
		*/

		//	Surfaces Settings

		Zone.SurfacesToPlot[GroupIndex] =			SurfacesToPlot_ExposedCellFaces;	//	SurfacesToPlot_e
		Zone.RangeBeginI[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeEndI[GroupIndex] =					0;							//	LgIndex_t
		Zone.RangeSkipI[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeBeginJ[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeEndJ[GroupIndex] =					0;							//	LgIndex_t
		Zone.RangeSkipJ[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeBeginK[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeEndK[GroupIndex] =					0;							//	LgIndex_t
		Zone.RangeSkipK[GroupIndex] =				1;							//	LgIndex_t
		/*
			0 is the same as maximum, which sets the end range to the maximum
			index for the dimension. The skip number is to skip every # planes,
			so 1 is to not skip any planes, (show every 1 plane).
		*/

		//	Volume Settings

		Zone.StreamtracesShow[GroupIndex] =			TRUE;						//	Boolean_t
		Zone.IsoSurfacesShow[GroupIndex] =			TRUE;						//	Boolean_t
		Zone.SlicesShow[GroupIndex] =				TRUE;						//	Boolean_t

		//	Effects Settings

		Zone.UseValueBlanking[GroupIndex] =			TRUE;						//	Boolean_t
		//Zone.UseClipPlanes[GroupIndex] =			TRUE;						//	Boolean_t
		Zone.UseSurfaceTranslucency[GroupIndex] =	FALSE;						//	Boolean_t
		Zone.SurfaceTranslucency[GroupIndex] =		50;		//	value is in %	//	SmInteger_t
		Zone.LightingEffect[GroupIndex] =			LightingEffect_Gouraud;		//	LightingEffect_e

		//	Group index and number (increments for next set)

		GroupIndex++;														//	int


		//===============================================================
		//	BCLinesFF Settings
		//===============================================================

		Zone.ZoneShow[GroupIndex] =					TRUE;						//	Boolean_t

		//	Mesh Settings

		Zone.MeshShow[GroupIndex] =					TRUE;						//	Boolean_t
		Zone.MeshType[GroupIndex] =					MeshType_Overlay;			//	MeshType_e
		Zone.MeshColor[GroupIndex] =					MultiColor8_C;					//	ColorIndex_t
		/*
			A number of color index constants are \#defined. These include Black_C, Red_C,
			Green_C, Blue_C, Cyan_C, Yellow_C, Purple_C, White_C, Custom1_C through
			Custom56_C, MultiColor_C, NoColor_C, MulitiColor2_C through MulitiColor8_C,
			RGBColor_C, and InvalidColor_C.  
		*/
		Zone.MeshLinePttrn[GroupIndex] =				LinePattern_DashDotDot;			//	LinePattern_e
		Zone.MeshLinePttrnLen[GroupIndex] =			2;		//	Value is in %	//	double
		Zone.MeshLineThickness[GroupIndex] =			0.2;	//	Value is in %	//	double


		//	Contour Settings

		Zone.ContourShow[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.ContourType[GroupIndex] =				ContourType_Flood;			//	ContourType_e
		Zone.FloodColoring[GroupIndex] =				ContourColoring_Group1;		//	ContourColoring_e
		Zone.ContourLineGroup[GroupIndex] =			1;							//	SmInteger_t
		Zone.ContourLineColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.ContourLinePttrn[GroupIndex] =			LinePattern_Solid;			//	LinePattern_e
		Zone.ContourLinePttrnLen[GroupIndex] =		2;		//	Value is in %	//	double
		Zone.ContourLineThickness[GroupIndex] =		0.1;	//	Value is in %	//	double
		Zone.ContourLighting[GroupIndex] =			TRUE;						//	Boolean_t

		//	Vector Settings
		
		Zone.VectorShow[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.VectorType[GroupIndex] =				VectorType_TailAtPoint;		//	VectorType_e
		Zone.VectorArrowType[GroupIndex] =			ArrowheadStyle_Plain;		//	ArrowheadStyle_e
		Zone.VectorLineColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.VectorTangent[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.VectorLinePttrn[GroupIndex] =			LinePattern_Solid;			//	LinePattern_e
		Zone.VectorLinePttrnLen[GroupIndex] =		2;		//	Value is in %	//	double
		Zone.VectorLineThickness[GroupIndex] =		0.1;	//	Value is in %	//	double

		//	Scatter Settings

		Zone.ScatterShow[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.ScatterShape[GroupIndex] =				GeomShape_Sphere;			//	GeomShape_e
		Zone.ScatterOutColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.ScatterInColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.ScatterFillMode[GroupIndex] =			FillMode_UseLineColor;		//	FillMode_e
		Zone.ScatterSize[GroupIndex] =				2.5;	//	value is in %	//	double
		Zone.ScatterLineThickness[GroupIndex] =		0.1;	//	value is in %	//	double

		//	Shade Settings

		Zone.ShadeShow[GroupIndex] =					FALSE;						//	Boolean_t
		Zone.ShadeColor[GroupIndex] =				Black_C;					//	ColorIndex_t
		Zone.ShadeLighting[GroupIndex] =				TRUE;						//	Boolean_t

		//	Edge Settings

		Zone.EdgeShow[GroupIndex] =					FALSE;						//	Boolean_t
		Zone.EdgeType[GroupIndex] =					EdgeType_BordersAndCreases;	//	EdgeType_e
		Zone.EdgeIBorder[GroupIndex] =				BorderLocation_Both;		//	BorderLocation_e
		Zone.EdgeJBorder[GroupIndex] =				BorderLocation_Both;		//	BorderLocation_e
		Zone.EdgeKBorder[GroupIndex] =				BorderLocation_Both;		//	BorderLocation_e
		Zone.EdgeColor[GroupIndex] =					Black_C;					//	ColorIndex_t
		Zone.EdgeLineThickness[GroupIndex] =			0.1;	//	value is in %	//	double

		//	Points Settings

		Zone.PointsToPlot[GroupIndex] =				PointsToPlot_SurfaceNodes;	//	PointsToPlot_e
		Zone.PointsIndexSkipI[GroupIndex] =			1;							//	LgIndex_t
		Zone.PointsIndexSkipJ[GroupIndex] =			1;							//	LgIndex_t
		Zone.PointsIndexSkipK[GroupIndex] =			1;							//	LgIndex_t
		/*
			1 is the same as "No Skip" in the GUI.
		*/

		//	Surfaces Settings

		Zone.SurfacesToPlot[GroupIndex] =			SurfacesToPlot_ExposedCellFaces;	//	SurfacesToPlot_e
		Zone.RangeBeginI[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeEndI[GroupIndex] =					0;							//	LgIndex_t
		Zone.RangeSkipI[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeBeginJ[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeEndJ[GroupIndex] =					0;							//	LgIndex_t
		Zone.RangeSkipJ[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeBeginK[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeEndK[GroupIndex] =					0;							//	LgIndex_t
		Zone.RangeSkipK[GroupIndex] =				1;							//	LgIndex_t
		/*
			0 is the same as maximum, which sets the end range to the maximum
			index for the dimension. The skip number is to skip every # planes,
			so 1 is to not skip any planes, (show every 1 plane).
		*/

		//	Volume Settings

		Zone.StreamtracesShow[GroupIndex] =			TRUE;						//	Boolean_t
		Zone.IsoSurfacesShow[GroupIndex] =			TRUE;						//	Boolean_t
		Zone.SlicesShow[GroupIndex] =				TRUE;						//	Boolean_t

		//	Effects Settings

		Zone.UseValueBlanking[GroupIndex] =			TRUE;						//	Boolean_t
		//Zone.UseClipPlanes[GroupIndex] =			TRUE;						//	Boolean_t
		Zone.UseSurfaceTranslucency[GroupIndex] =	FALSE;						//	Boolean_t
		Zone.SurfaceTranslucency[GroupIndex] =		50;		//	value is in %	//	SmInteger_t
		Zone.LightingEffect[GroupIndex] =			LightingEffect_Gouraud;		//	LightingEffect_e

		//	Group index and number (increments for next set)

		GroupIndex++;														//	int


		//===============================================================
		//	RCLines Settings
		//===============================================================

		Zone.ZoneShow[GroupIndex] =					TRUE;						//	Boolean_t

		//	Mesh Settings

		Zone.MeshShow[GroupIndex] =					TRUE;						//	Boolean_t
		Zone.MeshType[GroupIndex] =					MeshType_Overlay;			//	MeshType_e
		Zone.MeshColor[GroupIndex] =					MultiColor8_C;					//	ColorIndex_t
		/*
			A number of color index constants are \#defined. These include Black_C, Red_C,
			Green_C, Blue_C, Cyan_C, Yellow_C, Purple_C, White_C, Custom1_C through
			Custom56_C, MultiColor_C, NoColor_C, MulitiColor2_C through MulitiColor8_C,
			RGBColor_C, and InvalidColor_C.  
		*/
		Zone.MeshLinePttrn[GroupIndex] =				LinePattern_Dotted;			//	LinePattern_e
		Zone.MeshLinePttrnLen[GroupIndex] =			2;		//	Value is in %	//	double
		Zone.MeshLineThickness[GroupIndex] =			0.2;	//	Value is in %	//	double


		//	Contour Settings

		Zone.ContourShow[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.ContourType[GroupIndex] =				ContourType_Flood;			//	ContourType_e
		Zone.FloodColoring[GroupIndex] =				ContourColoring_Group1;		//	ContourColoring_e
		Zone.ContourLineGroup[GroupIndex] =			1;							//	SmInteger_t
		Zone.ContourLineColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.ContourLinePttrn[GroupIndex] =			LinePattern_Solid;			//	LinePattern_e
		Zone.ContourLinePttrnLen[GroupIndex] =		2;		//	Value is in %	//	double
		Zone.ContourLineThickness[GroupIndex] =		0.1;	//	Value is in %	//	double
		Zone.ContourLighting[GroupIndex] =			TRUE;						//	Boolean_t

		//	Vector Settings
		
		Zone.VectorShow[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.VectorType[GroupIndex] =				VectorType_TailAtPoint;		//	VectorType_e
		Zone.VectorArrowType[GroupIndex] =			ArrowheadStyle_Plain;		//	ArrowheadStyle_e
		Zone.VectorLineColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.VectorTangent[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.VectorLinePttrn[GroupIndex] =			LinePattern_Solid;			//	LinePattern_e
		Zone.VectorLinePttrnLen[GroupIndex] =		2;		//	Value is in %	//	double
		Zone.VectorLineThickness[GroupIndex] =		0.1;	//	Value is in %	//	double

		//	Scatter Settings

		Zone.ScatterShow[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.ScatterShape[GroupIndex] =				GeomShape_Sphere;			//	GeomShape_e
		Zone.ScatterOutColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.ScatterInColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.ScatterFillMode[GroupIndex] =			FillMode_UseLineColor;		//	FillMode_e
		Zone.ScatterSize[GroupIndex] =				2.5;	//	value is in %	//	double
		Zone.ScatterLineThickness[GroupIndex] =		0.1;	//	value is in %	//	double

		//	Shade Settings

		Zone.ShadeShow[GroupIndex] =					FALSE;						//	Boolean_t
		Zone.ShadeColor[GroupIndex] =				Black_C;					//	ColorIndex_t
		Zone.ShadeLighting[GroupIndex] =				TRUE;						//	Boolean_t

		//	Edge Settings

		Zone.EdgeShow[GroupIndex] =					FALSE;						//	Boolean_t
		Zone.EdgeType[GroupIndex] =					EdgeType_BordersAndCreases;	//	EdgeType_e
		Zone.EdgeIBorder[GroupIndex] =				BorderLocation_Both;		//	BorderLocation_e
		Zone.EdgeJBorder[GroupIndex] =				BorderLocation_Both;		//	BorderLocation_e
		Zone.EdgeKBorder[GroupIndex] =				BorderLocation_Both;		//	BorderLocation_e
		Zone.EdgeColor[GroupIndex] =					Black_C;					//	ColorIndex_t
		Zone.EdgeLineThickness[GroupIndex] =			0.1;	//	value is in %	//	double

		//	Points Settings

		Zone.PointsToPlot[GroupIndex] =				PointsToPlot_SurfaceNodes;	//	PointsToPlot_e
		Zone.PointsIndexSkipI[GroupIndex] =			1;							//	LgIndex_t
		Zone.PointsIndexSkipJ[GroupIndex] =			1;							//	LgIndex_t
		Zone.PointsIndexSkipK[GroupIndex] =			1;							//	LgIndex_t
		/*
			1 is the same as "No Skip" in the GUI.
		*/

		//	Surfaces Settings

		Zone.SurfacesToPlot[GroupIndex] =			SurfacesToPlot_ExposedCellFaces;	//	SurfacesToPlot_e
		Zone.RangeBeginI[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeEndI[GroupIndex] =					0;							//	LgIndex_t
		Zone.RangeSkipI[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeBeginJ[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeEndJ[GroupIndex] =					0;							//	LgIndex_t
		Zone.RangeSkipJ[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeBeginK[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeEndK[GroupIndex] =					0;							//	LgIndex_t
		Zone.RangeSkipK[GroupIndex] =				1;							//	LgIndex_t
		/*
			0 is the same as maximum, which sets the end range to the maximum
			index for the dimension. The skip number is to skip every # planes,
			so 1 is to not skip any planes, (show every 1 plane).
		*/

		//	Volume Settings

		Zone.StreamtracesShow[GroupIndex] =			TRUE;						//	Boolean_t
		Zone.IsoSurfacesShow[GroupIndex] =			TRUE;						//	Boolean_t
		Zone.SlicesShow[GroupIndex] =				TRUE;						//	Boolean_t

		//	Effects Settings

		Zone.UseValueBlanking[GroupIndex] =			TRUE;						//	Boolean_t
		//Zone.UseClipPlanes[GroupIndex] =			TRUE;						//	Boolean_t
		Zone.UseSurfaceTranslucency[GroupIndex] =	FALSE;						//	Boolean_t
		Zone.SurfaceTranslucency[GroupIndex] =		50;		//	value is in %	//	SmInteger_t
		Zone.LightingEffect[GroupIndex] =			LightingEffect_Gouraud;		//	LightingEffect_e

		//	Group index and number (increments for next set)

		GroupIndex++;														//	int


		//===============================================================
		//	RCLinesFF Settings
		//===============================================================

		Zone.ZoneShow[GroupIndex] =					TRUE;						//	Boolean_t

		//	Mesh Settings

		Zone.MeshShow[GroupIndex] =					TRUE;						//	Boolean_t
		Zone.MeshType[GroupIndex] =					MeshType_Overlay;			//	MeshType_e
		Zone.MeshColor[GroupIndex] =					MultiColor8_C;					//	ColorIndex_t
		/*
			A number of color index constants are \#defined. These include Black_C, Red_C,
			Green_C, Blue_C, Cyan_C, Yellow_C, Purple_C, White_C, Custom1_C through
			Custom56_C, MultiColor_C, NoColor_C, MulitiColor2_C through MulitiColor8_C,
			RGBColor_C, and InvalidColor_C.  
		*/
		Zone.MeshLinePttrn[GroupIndex] =				LinePattern_Dotted;			//	LinePattern_e
		Zone.MeshLinePttrnLen[GroupIndex] =			2;		//	Value is in %	//	double
		Zone.MeshLineThickness[GroupIndex] =			0.2;	//	Value is in %	//	double


		//	Contour Settings

		Zone.ContourShow[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.ContourType[GroupIndex] =				ContourType_Flood;			//	ContourType_e
		Zone.FloodColoring[GroupIndex] =				ContourColoring_Group1;		//	ContourColoring_e
		Zone.ContourLineGroup[GroupIndex] =			1;							//	SmInteger_t
		Zone.ContourLineColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.ContourLinePttrn[GroupIndex] =			LinePattern_Solid;			//	LinePattern_e
		Zone.ContourLinePttrnLen[GroupIndex] =		2;		//	Value is in %	//	double
		Zone.ContourLineThickness[GroupIndex] =		0.1;	//	Value is in %	//	double
		Zone.ContourLighting[GroupIndex] =			TRUE;						//	Boolean_t

		//	Vector Settings
		
		Zone.VectorShow[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.VectorType[GroupIndex] =				VectorType_TailAtPoint;		//	VectorType_e
		Zone.VectorArrowType[GroupIndex] =			ArrowheadStyle_Plain;		//	ArrowheadStyle_e
		Zone.VectorLineColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.VectorTangent[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.VectorLinePttrn[GroupIndex] =			LinePattern_Solid;			//	LinePattern_e
		Zone.VectorLinePttrnLen[GroupIndex] =		2;		//	Value is in %	//	double
		Zone.VectorLineThickness[GroupIndex] =		0.1;	//	Value is in %	//	double

		//	Scatter Settings

		Zone.ScatterShow[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.ScatterShape[GroupIndex] =				GeomShape_Sphere;			//	GeomShape_e
		Zone.ScatterOutColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.ScatterInColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.ScatterFillMode[GroupIndex] =			FillMode_UseLineColor;		//	FillMode_e
		Zone.ScatterSize[GroupIndex] =				2.5;	//	value is in %	//	double
		Zone.ScatterLineThickness[GroupIndex] =		0.1;	//	value is in %	//	double

		//	Shade Settings

		Zone.ShadeShow[GroupIndex] =					FALSE;						//	Boolean_t
		Zone.ShadeColor[GroupIndex] =				Black_C;					//	ColorIndex_t
		Zone.ShadeLighting[GroupIndex] =				TRUE;						//	Boolean_t

		//	Edge Settings

		Zone.EdgeShow[GroupIndex] =					FALSE;						//	Boolean_t
		Zone.EdgeType[GroupIndex] =					EdgeType_BordersAndCreases;	//	EdgeType_e
		Zone.EdgeIBorder[GroupIndex] =				BorderLocation_Both;		//	BorderLocation_e
		Zone.EdgeJBorder[GroupIndex] =				BorderLocation_Both;		//	BorderLocation_e
		Zone.EdgeKBorder[GroupIndex] =				BorderLocation_Both;		//	BorderLocation_e
		Zone.EdgeColor[GroupIndex] =					Black_C;					//	ColorIndex_t
		Zone.EdgeLineThickness[GroupIndex] =			0.1;	//	value is in %	//	double

		//	Points Settings

		Zone.PointsToPlot[GroupIndex] =				PointsToPlot_SurfaceNodes;	//	PointsToPlot_e
		Zone.PointsIndexSkipI[GroupIndex] =			1;							//	LgIndex_t
		Zone.PointsIndexSkipJ[GroupIndex] =			1;							//	LgIndex_t
		Zone.PointsIndexSkipK[GroupIndex] =			1;							//	LgIndex_t
		/*
			1 is the same as "No Skip" in the GUI.
		*/

		//	Surfaces Settings

		Zone.SurfacesToPlot[GroupIndex] =			SurfacesToPlot_ExposedCellFaces;	//	SurfacesToPlot_e
		Zone.RangeBeginI[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeEndI[GroupIndex] =					0;							//	LgIndex_t
		Zone.RangeSkipI[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeBeginJ[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeEndJ[GroupIndex] =					0;							//	LgIndex_t
		Zone.RangeSkipJ[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeBeginK[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeEndK[GroupIndex] =					0;							//	LgIndex_t
		Zone.RangeSkipK[GroupIndex] =				1;							//	LgIndex_t
		/*
			0 is the same as maximum, which sets the end range to the maximum
			index for the dimension. The skip number is to skip every # planes,
			so 1 is to not skip any planes, (show every 1 plane).
		*/

		//	Volume Settings

		Zone.StreamtracesShow[GroupIndex] =			TRUE;						//	Boolean_t
		Zone.IsoSurfacesShow[GroupIndex] =			TRUE;						//	Boolean_t
		Zone.SlicesShow[GroupIndex] =				TRUE;						//	Boolean_t

		//	Effects Settings

		Zone.UseValueBlanking[GroupIndex] =			TRUE;						//	Boolean_t
		//Zone.UseClipPlanes[GroupIndex] =			TRUE;						//	Boolean_t
		Zone.UseSurfaceTranslucency[GroupIndex] =	FALSE;						//	Boolean_t
		Zone.SurfaceTranslucency[GroupIndex] =		50;		//	value is in %	//	SmInteger_t
		Zone.LightingEffect[GroupIndex] =			LightingEffect_Gouraud;		//	LightingEffect_e

		//	Group index and number (increments for next set)

		GroupIndex++;														//	int


		//===============================================================
		//	ABRSurfaces Settings
		//===============================================================

		Zone.ZoneShow[GroupIndex] =					FALSE;						//	Boolean_t

		//	Mesh Settings

		Zone.MeshShow[GroupIndex] =					TRUE;						//	Boolean_t
		Zone.MeshType[GroupIndex] =					MeshType_Overlay;			//	MeshType_e
		Zone.MeshColor[GroupIndex] =					Black_C;					//	ColorIndex_t
		/*
			A number of color index constants are \#defined. These include Black_C, Red_C,
			Green_C, Blue_C, Cyan_C, Yellow_C, Purple_C, White_C, Custom1_C through
			Custom56_C, MultiColor_C, NoColor_C, MulitiColor2_C through MulitiColor8_C,
			RGBColor_C, and InvalidColor_C.  
		*/
		Zone.MeshLinePttrn[GroupIndex] =				LinePattern_Solid;			//	LinePattern_e
		Zone.MeshLinePttrnLen[GroupIndex] =			2;		//	Value is in %	//	double
		Zone.MeshLineThickness[GroupIndex] =			0.1;	//	Value is in %	//	double


		//	Contour Settings

		Zone.ContourShow[GroupIndex] =				TRUE;						//	Boolean_t
		Zone.ContourType[GroupIndex] =				ContourType_Flood;			//	ContourType_e
		Zone.FloodColoring[GroupIndex] =				ContourColoring_Group1;		//	ContourColoring_e
		Zone.ContourLineGroup[GroupIndex] =			1;							//	SmInteger_t
		Zone.ContourLineColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.ContourLinePttrn[GroupIndex] =			LinePattern_Solid;			//	LinePattern_e
		Zone.ContourLinePttrnLen[GroupIndex] =		2;		//	Value is in %	//	double
		Zone.ContourLineThickness[GroupIndex] =		0.1;	//	Value is in %	//	double
		Zone.ContourLighting[GroupIndex] =			TRUE;						//	Boolean_t

		//	Vector Settings
		
		Zone.VectorShow[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.VectorType[GroupIndex] =				VectorType_TailAtPoint;		//	VectorType_e
		Zone.VectorArrowType[GroupIndex] =			ArrowheadStyle_Plain;		//	ArrowheadStyle_e
		Zone.VectorLineColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.VectorTangent[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.VectorLinePttrn[GroupIndex] =			LinePattern_Solid;			//	LinePattern_e
		Zone.VectorLinePttrnLen[GroupIndex] =		2;		//	Value is in %	//	double
		Zone.VectorLineThickness[GroupIndex] =		0.1;	//	Value is in %	//	double

		//	Scatter Settings

		Zone.ScatterShow[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.ScatterShape[GroupIndex] =				GeomShape_Sphere;			//	GeomShape_e
		Zone.ScatterOutColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.ScatterInColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.ScatterFillMode[GroupIndex] =			FillMode_UseLineColor;		//	FillMode_e
		Zone.ScatterSize[GroupIndex] =				2.5;	//	value is in %	//	double
		Zone.ScatterLineThickness[GroupIndex] =		0.1;	//	value is in %	//	double

		//	Shade Settings

		Zone.ShadeShow[GroupIndex] =					TRUE;						//	Boolean_t
		Zone.ShadeColor[GroupIndex] =				Black_C;					//	ColorIndex_t
		Zone.ShadeLighting[GroupIndex] =				TRUE;						//	Boolean_t

		//	Edge Settings

		Zone.EdgeShow[GroupIndex] =					FALSE;						//	Boolean_t
		Zone.EdgeType[GroupIndex] =					EdgeType_BordersAndCreases;	//	EdgeType_e
		Zone.EdgeIBorder[GroupIndex] =				BorderLocation_Both;		//	BorderLocation_e
		Zone.EdgeJBorder[GroupIndex] =				BorderLocation_Both;		//	BorderLocation_e
		Zone.EdgeKBorder[GroupIndex] =				BorderLocation_Both;		//	BorderLocation_e
		Zone.EdgeColor[GroupIndex] =					Black_C;					//	ColorIndex_t
		Zone.EdgeLineThickness[GroupIndex] =			0.1;	//	value is in %	//	double

		//	Points Settings

		Zone.PointsToPlot[GroupIndex] =				PointsToPlot_SurfaceNodes;	//	PointsToPlot_e
		Zone.PointsIndexSkipI[GroupIndex] =			1;							//	LgIndex_t
		Zone.PointsIndexSkipJ[GroupIndex] =			1;							//	LgIndex_t
		Zone.PointsIndexSkipK[GroupIndex] =			1;							//	LgIndex_t
		/*
			1 is the same as "No Skip" in the GUI.
		*/

		//	Surfaces Settings

		Zone.SurfacesToPlot[GroupIndex] =			SurfacesToPlot_ExposedCellFaces;	//	SurfacesToPlot_e
		Zone.RangeBeginI[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeEndI[GroupIndex] =					0;							//	LgIndex_t
		Zone.RangeSkipI[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeBeginJ[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeEndJ[GroupIndex] =					0;							//	LgIndex_t
		Zone.RangeSkipJ[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeBeginK[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeEndK[GroupIndex] =					0;							//	LgIndex_t
		Zone.RangeSkipK[GroupIndex] =				1;							//	LgIndex_t
		/*
			0 is the same as maximum, which sets the end range to the maximum
			index for the dimension. The skip number is to skip every # planes,
			so 1 is to not skip any planes, (show every 1 plane).
		*/

		//	Volume Settings

		Zone.StreamtracesShow[GroupIndex] =			TRUE;						//	Boolean_t
		Zone.IsoSurfacesShow[GroupIndex] =			TRUE;						//	Boolean_t
		Zone.SlicesShow[GroupIndex] =				TRUE;						//	Boolean_t

		//	Effects Settings

		Zone.UseValueBlanking[GroupIndex] =			TRUE;						//	Boolean_t
		//Zone.UseClipPlanes[GroupIndex] =			TRUE;						//	Boolean_t
		Zone.UseSurfaceTranslucency[GroupIndex] =	TRUE;						//	Boolean_t
		Zone.SurfaceTranslucency[GroupIndex] =		50;		//	value is in %	//	SmInteger_t
		Zone.LightingEffect[GroupIndex] =			LightingEffect_Gouraud;		//	LightingEffect_e

		//	Group index and number (increments for next set)

		GroupIndex++;														//	int


		//===============================================================
		//	ABRSurfacesFF Settings
		//===============================================================

		Zone.ZoneShow[GroupIndex] =					FALSE;						//	Boolean_t

		//	Mesh Settings

		Zone.MeshShow[GroupIndex] =					TRUE;						//	Boolean_t
		Zone.MeshType[GroupIndex] =					MeshType_Overlay;			//	MeshType_e
		Zone.MeshColor[GroupIndex] =					Black_C;					//	ColorIndex_t
		/*
			A number of color index constants are \#defined. These include Black_C, Red_C,
			Green_C, Blue_C, Cyan_C, Yellow_C, Purple_C, White_C, Custom1_C through
			Custom56_C, MultiColor_C, NoColor_C, MulitiColor2_C through MulitiColor8_C,
			RGBColor_C, and InvalidColor_C.  
		*/
		Zone.MeshLinePttrn[GroupIndex] =				LinePattern_Solid;			//	LinePattern_e
		Zone.MeshLinePttrnLen[GroupIndex] =			2;		//	Value is in %	//	double
		Zone.MeshLineThickness[GroupIndex] =			0.1;	//	Value is in %	//	double


		//	Contour Settings

		Zone.ContourShow[GroupIndex] =				TRUE;						//	Boolean_t
		Zone.ContourType[GroupIndex] =				ContourType_Flood;			//	ContourType_e
		Zone.FloodColoring[GroupIndex] =				ContourColoring_Group1;		//	ContourColoring_e
		Zone.ContourLineGroup[GroupIndex] =			1;							//	SmInteger_t
		Zone.ContourLineColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.ContourLinePttrn[GroupIndex] =			LinePattern_Solid;			//	LinePattern_e
		Zone.ContourLinePttrnLen[GroupIndex] =		2;		//	Value is in %	//	double
		Zone.ContourLineThickness[GroupIndex] =		0.1;	//	Value is in %	//	double
		Zone.ContourLighting[GroupIndex] =			TRUE;						//	Boolean_t

		//	Vector Settings
		
		Zone.VectorShow[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.VectorType[GroupIndex] =				VectorType_TailAtPoint;		//	VectorType_e
		Zone.VectorArrowType[GroupIndex] =			ArrowheadStyle_Plain;		//	ArrowheadStyle_e
		Zone.VectorLineColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.VectorTangent[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.VectorLinePttrn[GroupIndex] =			LinePattern_Solid;			//	LinePattern_e
		Zone.VectorLinePttrnLen[GroupIndex] =		2;		//	Value is in %	//	double
		Zone.VectorLineThickness[GroupIndex] =		0.1;	//	Value is in %	//	double

		//	Scatter Settings

		Zone.ScatterShow[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.ScatterShape[GroupIndex] =				GeomShape_Sphere;			//	GeomShape_e
		Zone.ScatterOutColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.ScatterInColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.ScatterFillMode[GroupIndex] =			FillMode_UseLineColor;		//	FillMode_e
		Zone.ScatterSize[GroupIndex] =				2.5;	//	value is in %	//	double
		Zone.ScatterLineThickness[GroupIndex] =		0.1;	//	value is in %	//	double

		//	Shade Settings

		Zone.ShadeShow[GroupIndex] =					TRUE;						//	Boolean_t
		Zone.ShadeColor[GroupIndex] =				Black_C;					//	ColorIndex_t
		Zone.ShadeLighting[GroupIndex] =				TRUE;						//	Boolean_t

		//	Edge Settings

		Zone.EdgeShow[GroupIndex] =					FALSE;						//	Boolean_t
		Zone.EdgeType[GroupIndex] =					EdgeType_BordersAndCreases;	//	EdgeType_e
		Zone.EdgeIBorder[GroupIndex] =				BorderLocation_Both;		//	BorderLocation_e
		Zone.EdgeJBorder[GroupIndex] =				BorderLocation_Both;		//	BorderLocation_e
		Zone.EdgeKBorder[GroupIndex] =				BorderLocation_Both;		//	BorderLocation_e
		Zone.EdgeColor[GroupIndex] =					Black_C;					//	ColorIndex_t
		Zone.EdgeLineThickness[GroupIndex] =			0.1;	//	value is in %	//	double

		//	Points Settings

		Zone.PointsToPlot[GroupIndex] =				PointsToPlot_SurfaceNodes;	//	PointsToPlot_e
		Zone.PointsIndexSkipI[GroupIndex] =			1;							//	LgIndex_t
		Zone.PointsIndexSkipJ[GroupIndex] =			1;							//	LgIndex_t
		Zone.PointsIndexSkipK[GroupIndex] =			1;							//	LgIndex_t
		/*
			1 is the same as "No Skip" in the GUI.
		*/

		//	Surfaces Settings

		Zone.SurfacesToPlot[GroupIndex] =			SurfacesToPlot_ExposedCellFaces;	//	SurfacesToPlot_e
		Zone.RangeBeginI[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeEndI[GroupIndex] =					0;							//	LgIndex_t
		Zone.RangeSkipI[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeBeginJ[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeEndJ[GroupIndex] =					0;							//	LgIndex_t
		Zone.RangeSkipJ[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeBeginK[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeEndK[GroupIndex] =					0;							//	LgIndex_t
		Zone.RangeSkipK[GroupIndex] =				1;							//	LgIndex_t
		/*
			0 is the same as maximum, which sets the end range to the maximum
			index for the dimension. The skip number is to skip every # planes,
			so 1 is to not skip any planes, (show every 1 plane).
		*/

		//	Volume Settings

		Zone.StreamtracesShow[GroupIndex] =			TRUE;						//	Boolean_t
		Zone.IsoSurfacesShow[GroupIndex] =			TRUE;						//	Boolean_t
		Zone.SlicesShow[GroupIndex] =				TRUE;						//	Boolean_t

		//	Effects Settings

		Zone.UseValueBlanking[GroupIndex] =			TRUE;						//	Boolean_t
		//Zone.UseClipPlanes[GroupIndex] =			TRUE;						//	Boolean_t
		Zone.UseSurfaceTranslucency[GroupIndex] =	TRUE;						//	Boolean_t
		Zone.SurfaceTranslucency[GroupIndex] =		50;		//	value is in %	//	SmInteger_t
		Zone.LightingEffect[GroupIndex] =			LightingEffect_Gouraud;		//	LightingEffect_e

		//	Group index and number (increments for next set)

		GroupIndex++;														//	int


		//===============================================================
		//	ABCSurfaces Settings
		//===============================================================

		Zone.ZoneShow[GroupIndex] =					FALSE;						//	Boolean_t

		//	Mesh Settings

		Zone.MeshShow[GroupIndex] =					TRUE;						//	Boolean_t
		Zone.MeshType[GroupIndex] =					MeshType_Overlay;			//	MeshType_e
		Zone.MeshColor[GroupIndex] =					Black_C;					//	ColorIndex_t
		/*
			A number of color index constants are \#defined. These include Black_C, Red_C,
			Green_C, Blue_C, Cyan_C, Yellow_C, Purple_C, White_C, Custom1_C through
			Custom56_C, MultiColor_C, NoColor_C, MulitiColor2_C through MulitiColor8_C,
			RGBColor_C, and InvalidColor_C.  
		*/
		Zone.MeshLinePttrn[GroupIndex] =				LinePattern_Solid;			//	LinePattern_e
		Zone.MeshLinePttrnLen[GroupIndex] =			2;		//	Value is in %	//	double
		Zone.MeshLineThickness[GroupIndex] =			0.1;	//	Value is in %	//	double


		//	Contour Settings

		Zone.ContourShow[GroupIndex] =				TRUE;						//	Boolean_t
		Zone.ContourType[GroupIndex] =				ContourType_Flood;			//	ContourType_e
		Zone.FloodColoring[GroupIndex] =				ContourColoring_Group1;		//	ContourColoring_e
		Zone.ContourLineGroup[GroupIndex] =			1;							//	SmInteger_t
		Zone.ContourLineColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.ContourLinePttrn[GroupIndex] =			LinePattern_Solid;			//	LinePattern_e
		Zone.ContourLinePttrnLen[GroupIndex] =		2;		//	Value is in %	//	double
		Zone.ContourLineThickness[GroupIndex] =		0.1;	//	Value is in %	//	double
		Zone.ContourLighting[GroupIndex] =			TRUE;						//	Boolean_t

		//	Vector Settings
		
		Zone.VectorShow[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.VectorType[GroupIndex] =				VectorType_TailAtPoint;		//	VectorType_e
		Zone.VectorArrowType[GroupIndex] =			ArrowheadStyle_Plain;		//	ArrowheadStyle_e
		Zone.VectorLineColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.VectorTangent[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.VectorLinePttrn[GroupIndex] =			LinePattern_Solid;			//	LinePattern_e
		Zone.VectorLinePttrnLen[GroupIndex] =		2;		//	Value is in %	//	double
		Zone.VectorLineThickness[GroupIndex] =		0.1;	//	Value is in %	//	double

		//	Scatter Settings

		Zone.ScatterShow[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.ScatterShape[GroupIndex] =				GeomShape_Sphere;			//	GeomShape_e
		Zone.ScatterOutColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.ScatterInColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.ScatterFillMode[GroupIndex] =			FillMode_UseLineColor;		//	FillMode_e
		Zone.ScatterSize[GroupIndex] =				2.5;	//	value is in %	//	double
		Zone.ScatterLineThickness[GroupIndex] =		0.1;	//	value is in %	//	double

		//	Shade Settings

		Zone.ShadeShow[GroupIndex] =					TRUE;						//	Boolean_t
		Zone.ShadeColor[GroupIndex] =				Black_C;					//	ColorIndex_t
		Zone.ShadeLighting[GroupIndex] =				TRUE;						//	Boolean_t

		//	Edge Settings

		Zone.EdgeShow[GroupIndex] =					FALSE;						//	Boolean_t
		Zone.EdgeType[GroupIndex] =					EdgeType_BordersAndCreases;	//	EdgeType_e
		Zone.EdgeIBorder[GroupIndex] =				BorderLocation_Both;		//	BorderLocation_e
		Zone.EdgeJBorder[GroupIndex] =				BorderLocation_Both;		//	BorderLocation_e
		Zone.EdgeKBorder[GroupIndex] =				BorderLocation_Both;		//	BorderLocation_e
		Zone.EdgeColor[GroupIndex] =					Black_C;					//	ColorIndex_t
		Zone.EdgeLineThickness[GroupIndex] =			0.1;	//	value is in %	//	double

		//	Points Settings

		Zone.PointsToPlot[GroupIndex] =				PointsToPlot_SurfaceNodes;	//	PointsToPlot_e
		Zone.PointsIndexSkipI[GroupIndex] =			1;							//	LgIndex_t
		Zone.PointsIndexSkipJ[GroupIndex] =			1;							//	LgIndex_t
		Zone.PointsIndexSkipK[GroupIndex] =			1;							//	LgIndex_t
		/*
			1 is the same as "No Skip" in the GUI.
		*/

		//	Surfaces Settings

		Zone.SurfacesToPlot[GroupIndex] =			SurfacesToPlot_ExposedCellFaces;	//	SurfacesToPlot_e
		Zone.RangeBeginI[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeEndI[GroupIndex] =					0;							//	LgIndex_t
		Zone.RangeSkipI[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeBeginJ[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeEndJ[GroupIndex] =					0;							//	LgIndex_t
		Zone.RangeSkipJ[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeBeginK[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeEndK[GroupIndex] =					0;							//	LgIndex_t
		Zone.RangeSkipK[GroupIndex] =				1;							//	LgIndex_t
		/*
			0 is the same as maximum, which sets the end range to the maximum
			index for the dimension. The skip number is to skip every # planes,
			so 1 is to not skip any planes, (show every 1 plane).
		*/

		//	Volume Settings

		Zone.StreamtracesShow[GroupIndex] =			TRUE;						//	Boolean_t
		Zone.IsoSurfacesShow[GroupIndex] =			TRUE;						//	Boolean_t
		Zone.SlicesShow[GroupIndex] =				TRUE;						//	Boolean_t

		//	Effects Settings

		Zone.UseValueBlanking[GroupIndex] =			TRUE;						//	Boolean_t
		//Zone.UseClipPlanes[GroupIndex] =			TRUE;						//	Boolean_t
		Zone.UseSurfaceTranslucency[GroupIndex] =	TRUE;						//	Boolean_t
		Zone.SurfaceTranslucency[GroupIndex] =		50;		//	value is in %	//	SmInteger_t
		Zone.LightingEffect[GroupIndex] =			LightingEffect_Gouraud;		//	LightingEffect_e

		//	Group index and number (increments for next set)

		GroupIndex++;														//	int


		//===============================================================
		//	ABCSurfacesFF Settings
		//===============================================================

		Zone.ZoneShow[GroupIndex] =					FALSE;						//	Boolean_t

		//	Mesh Settings

		Zone.MeshShow[GroupIndex] =					TRUE;						//	Boolean_t
		Zone.MeshType[GroupIndex] =					MeshType_Overlay;			//	MeshType_e
		Zone.MeshColor[GroupIndex] =					Black_C;					//	ColorIndex_t
		/*
			A number of color index constants are \#defined. These include Black_C, Red_C,
			Green_C, Blue_C, Cyan_C, Yellow_C, Purple_C, White_C, Custom1_C through
			Custom56_C, MultiColor_C, NoColor_C, MulitiColor2_C through MulitiColor8_C,
			RGBColor_C, and InvalidColor_C.  
		*/
		Zone.MeshLinePttrn[GroupIndex] =				LinePattern_Solid;			//	LinePattern_e
		Zone.MeshLinePttrnLen[GroupIndex] =			2;		//	Value is in %	//	double
		Zone.MeshLineThickness[GroupIndex] =			0.1;	//	Value is in %	//	double


		//	Contour Settings

		Zone.ContourShow[GroupIndex] =				TRUE;						//	Boolean_t
		Zone.ContourType[GroupIndex] =				ContourType_Flood;			//	ContourType_e
		Zone.FloodColoring[GroupIndex] =				ContourColoring_Group1;		//	ContourColoring_e
		Zone.ContourLineGroup[GroupIndex] =			1;							//	SmInteger_t
		Zone.ContourLineColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.ContourLinePttrn[GroupIndex] =			LinePattern_Solid;			//	LinePattern_e
		Zone.ContourLinePttrnLen[GroupIndex] =		2;		//	Value is in %	//	double
		Zone.ContourLineThickness[GroupIndex] =		0.1;	//	Value is in %	//	double
		Zone.ContourLighting[GroupIndex] =			TRUE;						//	Boolean_t

		//	Vector Settings
		
		Zone.VectorShow[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.VectorType[GroupIndex] =				VectorType_TailAtPoint;		//	VectorType_e
		Zone.VectorArrowType[GroupIndex] =			ArrowheadStyle_Plain;		//	ArrowheadStyle_e
		Zone.VectorLineColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.VectorTangent[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.VectorLinePttrn[GroupIndex] =			LinePattern_Solid;			//	LinePattern_e
		Zone.VectorLinePttrnLen[GroupIndex] =		2;		//	Value is in %	//	double
		Zone.VectorLineThickness[GroupIndex] =		0.1;	//	Value is in %	//	double

		//	Scatter Settings

		Zone.ScatterShow[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.ScatterShape[GroupIndex] =				GeomShape_Sphere;			//	GeomShape_e
		Zone.ScatterOutColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.ScatterInColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.ScatterFillMode[GroupIndex] =			FillMode_UseLineColor;		//	FillMode_e
		Zone.ScatterSize[GroupIndex] =				2.5;	//	value is in %	//	double
		Zone.ScatterLineThickness[GroupIndex] =		0.1;	//	value is in %	//	double

		//	Shade Settings

		Zone.ShadeShow[GroupIndex] =					TRUE;						//	Boolean_t
		Zone.ShadeColor[GroupIndex] =				Black_C;					//	ColorIndex_t
		Zone.ShadeLighting[GroupIndex] =				TRUE;						//	Boolean_t

		//	Edge Settings

		Zone.EdgeShow[GroupIndex] =					FALSE;						//	Boolean_t
		Zone.EdgeType[GroupIndex] =					EdgeType_BordersAndCreases;	//	EdgeType_e
		Zone.EdgeIBorder[GroupIndex] =				BorderLocation_Both;		//	BorderLocation_e
		Zone.EdgeJBorder[GroupIndex] =				BorderLocation_Both;		//	BorderLocation_e
		Zone.EdgeKBorder[GroupIndex] =				BorderLocation_Both;		//	BorderLocation_e
		Zone.EdgeColor[GroupIndex] =					Black_C;					//	ColorIndex_t
		Zone.EdgeLineThickness[GroupIndex] =			0.1;	//	value is in %	//	double

		//	Points Settings

		Zone.PointsToPlot[GroupIndex] =				PointsToPlot_SurfaceNodes;	//	PointsToPlot_e
		Zone.PointsIndexSkipI[GroupIndex] =			1;							//	LgIndex_t
		Zone.PointsIndexSkipJ[GroupIndex] =			1;							//	LgIndex_t
		Zone.PointsIndexSkipK[GroupIndex] =			1;							//	LgIndex_t
		/*
			1 is the same as "No Skip" in the GUI.
		*/

		//	Surfaces Settings

		Zone.SurfacesToPlot[GroupIndex] =			SurfacesToPlot_ExposedCellFaces;	//	SurfacesToPlot_e
		Zone.RangeBeginI[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeEndI[GroupIndex] =					0;							//	LgIndex_t
		Zone.RangeSkipI[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeBeginJ[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeEndJ[GroupIndex] =					0;							//	LgIndex_t
		Zone.RangeSkipJ[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeBeginK[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeEndK[GroupIndex] =					0;							//	LgIndex_t
		Zone.RangeSkipK[GroupIndex] =				1;							//	LgIndex_t
		/*
			0 is the same as maximum, which sets the end range to the maximum
			index for the dimension. The skip number is to skip every # planes,
			so 1 is to not skip any planes, (show every 1 plane).
		*/

		//	Volume Settings

		Zone.StreamtracesShow[GroupIndex] =			TRUE;						//	Boolean_t
		Zone.IsoSurfacesShow[GroupIndex] =			TRUE;						//	Boolean_t
		Zone.SlicesShow[GroupIndex] =				TRUE;						//	Boolean_t

		//	Effects Settings

		Zone.UseValueBlanking[GroupIndex] =			TRUE;						//	Boolean_t
		//Zone.UseClipPlanes[GroupIndex] =			TRUE;						//	Boolean_t
		Zone.UseSurfaceTranslucency[GroupIndex] =	TRUE;						//	Boolean_t
		Zone.SurfaceTranslucency[GroupIndex] =		50;		//	value is in %	//	SmInteger_t
		Zone.LightingEffect[GroupIndex] =			LightingEffect_Gouraud;		//	LightingEffect_e

		//	Group index and number (increments for next set)

		GroupIndex++;														//	int


		//===============================================================
		//	ARCSurfaces Settings
		//===============================================================

		Zone.ZoneShow[GroupIndex] =					FALSE;						//	Boolean_t

		//	Mesh Settings

		Zone.MeshShow[GroupIndex] =					TRUE;						//	Boolean_t
		Zone.MeshType[GroupIndex] =					MeshType_Overlay;			//	MeshType_e
		Zone.MeshColor[GroupIndex] =					Black_C;					//	ColorIndex_t
		/*
			A number of color index constants are \#defined. These include Black_C, Red_C,
			Green_C, Blue_C, Cyan_C, Yellow_C, Purple_C, White_C, Custom1_C through
			Custom56_C, MultiColor_C, NoColor_C, MulitiColor2_C through MulitiColor8_C,
			RGBColor_C, and InvalidColor_C.  
		*/
		Zone.MeshLinePttrn[GroupIndex] =				LinePattern_Solid;			//	LinePattern_e
		Zone.MeshLinePttrnLen[GroupIndex] =			2;		//	Value is in %	//	double
		Zone.MeshLineThickness[GroupIndex] =			0.1;	//	Value is in %	//	double


		//	Contour Settings

		Zone.ContourShow[GroupIndex] =				TRUE;						//	Boolean_t
		Zone.ContourType[GroupIndex] =				ContourType_Flood;			//	ContourType_e
		Zone.FloodColoring[GroupIndex] =				ContourColoring_Group1;		//	ContourColoring_e
		Zone.ContourLineGroup[GroupIndex] =			1;							//	SmInteger_t
		Zone.ContourLineColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.ContourLinePttrn[GroupIndex] =			LinePattern_Solid;			//	LinePattern_e
		Zone.ContourLinePttrnLen[GroupIndex] =		2;		//	Value is in %	//	double
		Zone.ContourLineThickness[GroupIndex] =		0.1;	//	Value is in %	//	double
		Zone.ContourLighting[GroupIndex] =			TRUE;						//	Boolean_t

		//	Vector Settings
		
		Zone.VectorShow[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.VectorType[GroupIndex] =				VectorType_TailAtPoint;		//	VectorType_e
		Zone.VectorArrowType[GroupIndex] =			ArrowheadStyle_Plain;		//	ArrowheadStyle_e
		Zone.VectorLineColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.VectorTangent[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.VectorLinePttrn[GroupIndex] =			LinePattern_Solid;			//	LinePattern_e
		Zone.VectorLinePttrnLen[GroupIndex] =		2;		//	Value is in %	//	double
		Zone.VectorLineThickness[GroupIndex] =		0.1;	//	Value is in %	//	double

		//	Scatter Settings

		Zone.ScatterShow[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.ScatterShape[GroupIndex] =				GeomShape_Sphere;			//	GeomShape_e
		Zone.ScatterOutColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.ScatterInColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.ScatterFillMode[GroupIndex] =			FillMode_UseLineColor;		//	FillMode_e
		Zone.ScatterSize[GroupIndex] =				2.5;	//	value is in %	//	double
		Zone.ScatterLineThickness[GroupIndex] =		0.1;	//	value is in %	//	double

		//	Shade Settings

		Zone.ShadeShow[GroupIndex] =					TRUE;						//	Boolean_t
		Zone.ShadeColor[GroupIndex] =				Black_C;					//	ColorIndex_t
		Zone.ShadeLighting[GroupIndex] =				TRUE;						//	Boolean_t

		//	Edge Settings

		Zone.EdgeShow[GroupIndex] =					FALSE;						//	Boolean_t
		Zone.EdgeType[GroupIndex] =					EdgeType_BordersAndCreases;	//	EdgeType_e
		Zone.EdgeIBorder[GroupIndex] =				BorderLocation_Both;		//	BorderLocation_e
		Zone.EdgeJBorder[GroupIndex] =				BorderLocation_Both;		//	BorderLocation_e
		Zone.EdgeKBorder[GroupIndex] =				BorderLocation_Both;		//	BorderLocation_e
		Zone.EdgeColor[GroupIndex] =					Black_C;					//	ColorIndex_t
		Zone.EdgeLineThickness[GroupIndex] =			0.1;	//	value is in %	//	double

		//	Points Settings

		Zone.PointsToPlot[GroupIndex] =				PointsToPlot_SurfaceNodes;	//	PointsToPlot_e
		Zone.PointsIndexSkipI[GroupIndex] =			1;							//	LgIndex_t
		Zone.PointsIndexSkipJ[GroupIndex] =			1;							//	LgIndex_t
		Zone.PointsIndexSkipK[GroupIndex] =			1;							//	LgIndex_t
		/*
			1 is the same as "No Skip" in the GUI.
		*/

		//	Surfaces Settings

		Zone.SurfacesToPlot[GroupIndex] =			SurfacesToPlot_ExposedCellFaces;	//	SurfacesToPlot_e
		Zone.RangeBeginI[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeEndI[GroupIndex] =					0;							//	LgIndex_t
		Zone.RangeSkipI[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeBeginJ[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeEndJ[GroupIndex] =					0;							//	LgIndex_t
		Zone.RangeSkipJ[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeBeginK[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeEndK[GroupIndex] =					0;							//	LgIndex_t
		Zone.RangeSkipK[GroupIndex] =				1;							//	LgIndex_t
		/*
			0 is the same as maximum, which sets the end range to the maximum
			index for the dimension. The skip number is to skip every # planes,
			so 1 is to not skip any planes, (show every 1 plane).
		*/

		//	Volume Settings

		Zone.StreamtracesShow[GroupIndex] =			TRUE;						//	Boolean_t
		Zone.IsoSurfacesShow[GroupIndex] =			TRUE;						//	Boolean_t
		Zone.SlicesShow[GroupIndex] =				TRUE;						//	Boolean_t

		//	Effects Settings

		Zone.UseValueBlanking[GroupIndex] =			TRUE;						//	Boolean_t
		//Zone.UseClipPlanes[GroupIndex] =			TRUE;						//	Boolean_t
		Zone.UseSurfaceTranslucency[GroupIndex] =	TRUE;						//	Boolean_t
		Zone.SurfaceTranslucency[GroupIndex] =		50;		//	value is in %	//	SmInteger_t
		Zone.LightingEffect[GroupIndex] =			LightingEffect_Gouraud;		//	LightingEffect_e

		//	Group index and number (increments for next set)

		GroupIndex++;														//	int


		//===============================================================
		//	ARCSurfacesFF Settings
		//===============================================================

		Zone.ZoneShow[GroupIndex] =					FALSE;						//	Boolean_t

		//	Mesh Settings

		Zone.MeshShow[GroupIndex] =					TRUE;						//	Boolean_t
		Zone.MeshType[GroupIndex] =					MeshType_Overlay;			//	MeshType_e
		Zone.MeshColor[GroupIndex] =					Black_C;					//	ColorIndex_t
		/*
			A number of color index constants are \#defined. These include Black_C, Red_C,
			Green_C, Blue_C, Cyan_C, Yellow_C, Purple_C, White_C, Custom1_C through
			Custom56_C, MultiColor_C, NoColor_C, MulitiColor2_C through MulitiColor8_C,
			RGBColor_C, and InvalidColor_C.  
		*/
		Zone.MeshLinePttrn[GroupIndex] =				LinePattern_Solid;			//	LinePattern_e
		Zone.MeshLinePttrnLen[GroupIndex] =			2;		//	Value is in %	//	double
		Zone.MeshLineThickness[GroupIndex] =			0.1;	//	Value is in %	//	double


		//	Contour Settings

		Zone.ContourShow[GroupIndex] =				TRUE;						//	Boolean_t
		Zone.ContourType[GroupIndex] =				ContourType_Flood;			//	ContourType_e
		Zone.FloodColoring[GroupIndex] =				ContourColoring_Group1;		//	ContourColoring_e
		Zone.ContourLineGroup[GroupIndex] =			1;							//	SmInteger_t
		Zone.ContourLineColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.ContourLinePttrn[GroupIndex] =			LinePattern_Solid;			//	LinePattern_e
		Zone.ContourLinePttrnLen[GroupIndex] =		2;		//	Value is in %	//	double
		Zone.ContourLineThickness[GroupIndex] =		0.1;	//	Value is in %	//	double
		Zone.ContourLighting[GroupIndex] =			TRUE;						//	Boolean_t

		//	Vector Settings
		
		Zone.VectorShow[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.VectorType[GroupIndex] =				VectorType_TailAtPoint;		//	VectorType_e
		Zone.VectorArrowType[GroupIndex] =			ArrowheadStyle_Plain;		//	ArrowheadStyle_e
		Zone.VectorLineColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.VectorTangent[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.VectorLinePttrn[GroupIndex] =			LinePattern_Solid;			//	LinePattern_e
		Zone.VectorLinePttrnLen[GroupIndex] =		2;		//	Value is in %	//	double
		Zone.VectorLineThickness[GroupIndex] =		0.1;	//	Value is in %	//	double

		//	Scatter Settings

		Zone.ScatterShow[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.ScatterShape[GroupIndex] =				GeomShape_Sphere;			//	GeomShape_e
		Zone.ScatterOutColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.ScatterInColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.ScatterFillMode[GroupIndex] =			FillMode_UseLineColor;		//	FillMode_e
		Zone.ScatterSize[GroupIndex] =				2.5;	//	value is in %	//	double
		Zone.ScatterLineThickness[GroupIndex] =		0.1;	//	value is in %	//	double

		//	Shade Settings

		Zone.ShadeShow[GroupIndex] =					TRUE;						//	Boolean_t
		Zone.ShadeColor[GroupIndex] =				Black_C;					//	ColorIndex_t
		Zone.ShadeLighting[GroupIndex] =				TRUE;						//	Boolean_t

		//	Edge Settings

		Zone.EdgeShow[GroupIndex] =					FALSE;						//	Boolean_t
		Zone.EdgeType[GroupIndex] =					EdgeType_BordersAndCreases;	//	EdgeType_e
		Zone.EdgeIBorder[GroupIndex] =				BorderLocation_Both;		//	BorderLocation_e
		Zone.EdgeJBorder[GroupIndex] =				BorderLocation_Both;		//	BorderLocation_e
		Zone.EdgeKBorder[GroupIndex] =				BorderLocation_Both;		//	BorderLocation_e
		Zone.EdgeColor[GroupIndex] =					Black_C;					//	ColorIndex_t
		Zone.EdgeLineThickness[GroupIndex] =			0.1;	//	value is in %	//	double

		//	Points Settings

		Zone.PointsToPlot[GroupIndex] =				PointsToPlot_SurfaceNodes;	//	PointsToPlot_e
		Zone.PointsIndexSkipI[GroupIndex] =			1;							//	LgIndex_t
		Zone.PointsIndexSkipJ[GroupIndex] =			1;							//	LgIndex_t
		Zone.PointsIndexSkipK[GroupIndex] =			1;							//	LgIndex_t
		/*
			1 is the same as "No Skip" in the GUI.
		*/

		//	Surfaces Settings

		Zone.SurfacesToPlot[GroupIndex] =			SurfacesToPlot_ExposedCellFaces;	//	SurfacesToPlot_e
		Zone.RangeBeginI[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeEndI[GroupIndex] =					0;							//	LgIndex_t
		Zone.RangeSkipI[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeBeginJ[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeEndJ[GroupIndex] =					0;							//	LgIndex_t
		Zone.RangeSkipJ[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeBeginK[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeEndK[GroupIndex] =					0;							//	LgIndex_t
		Zone.RangeSkipK[GroupIndex] =				1;							//	LgIndex_t
		/*
			0 is the same as maximum, which sets the end range to the maximum
			index for the dimension. The skip number is to skip every # planes,
			so 1 is to not skip any planes, (show every 1 plane).
		*/

		//	Volume Settings

		Zone.StreamtracesShow[GroupIndex] =			TRUE;						//	Boolean_t
		Zone.IsoSurfacesShow[GroupIndex] =			TRUE;						//	Boolean_t
		Zone.SlicesShow[GroupIndex] =				TRUE;						//	Boolean_t

		//	Effects Settings

		Zone.UseValueBlanking[GroupIndex] =			TRUE;						//	Boolean_t
		//Zone.UseClipPlanes[GroupIndex] =			TRUE;						//	Boolean_t
		Zone.UseSurfaceTranslucency[GroupIndex] =	TRUE;						//	Boolean_t
		Zone.SurfaceTranslucency[GroupIndex] =		50;		//	value is in %	//	SmInteger_t
		Zone.LightingEffect[GroupIndex] =			LightingEffect_Gouraud;		//	LightingEffect_e

		//	Group index and number (increments for next set)

		GroupIndex++;														//	int


		//===============================================================
		//	BRCSurfaces Settings
		//===============================================================

		Zone.ZoneShow[GroupIndex] =					FALSE;						//	Boolean_t

		//	Mesh Settings

		Zone.MeshShow[GroupIndex] =					TRUE;						//	Boolean_t
		Zone.MeshType[GroupIndex] =					MeshType_Overlay;			//	MeshType_e
		Zone.MeshColor[GroupIndex] =					Black_C;					//	ColorIndex_t
		/*
			A number of color index constants are \#defined. These include Black_C, Red_C,
			Green_C, Blue_C, Cyan_C, Yellow_C, Purple_C, White_C, Custom1_C through
			Custom56_C, MultiColor_C, NoColor_C, MulitiColor2_C through MulitiColor8_C,
			RGBColor_C, and InvalidColor_C.  
		*/
		Zone.MeshLinePttrn[GroupIndex] =				LinePattern_Solid;			//	LinePattern_e
		Zone.MeshLinePttrnLen[GroupIndex] =			2;		//	Value is in %	//	double
		Zone.MeshLineThickness[GroupIndex] =			0.1;	//	Value is in %	//	double


		//	Contour Settings

		Zone.ContourShow[GroupIndex] =				TRUE;						//	Boolean_t
		Zone.ContourType[GroupIndex] =				ContourType_Flood;			//	ContourType_e
		Zone.FloodColoring[GroupIndex] =				ContourColoring_Group1;		//	ContourColoring_e
		Zone.ContourLineGroup[GroupIndex] =			1;							//	SmInteger_t
		Zone.ContourLineColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.ContourLinePttrn[GroupIndex] =			LinePattern_Solid;			//	LinePattern_e
		Zone.ContourLinePttrnLen[GroupIndex] =		2;		//	Value is in %	//	double
		Zone.ContourLineThickness[GroupIndex] =		0.1;	//	Value is in %	//	double
		Zone.ContourLighting[GroupIndex] =			TRUE;						//	Boolean_t

		//	Vector Settings
		
		Zone.VectorShow[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.VectorType[GroupIndex] =				VectorType_TailAtPoint;		//	VectorType_e
		Zone.VectorArrowType[GroupIndex] =			ArrowheadStyle_Plain;		//	ArrowheadStyle_e
		Zone.VectorLineColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.VectorTangent[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.VectorLinePttrn[GroupIndex] =			LinePattern_Solid;			//	LinePattern_e
		Zone.VectorLinePttrnLen[GroupIndex] =		2;		//	Value is in %	//	double
		Zone.VectorLineThickness[GroupIndex] =		0.1;	//	Value is in %	//	double

		//	Scatter Settings

		Zone.ScatterShow[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.ScatterShape[GroupIndex] =				GeomShape_Sphere;			//	GeomShape_e
		Zone.ScatterOutColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.ScatterInColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.ScatterFillMode[GroupIndex] =			FillMode_UseLineColor;		//	FillMode_e
		Zone.ScatterSize[GroupIndex] =				2.5;	//	value is in %	//	double
		Zone.ScatterLineThickness[GroupIndex] =		0.1;	//	value is in %	//	double

		//	Shade Settings

		Zone.ShadeShow[GroupIndex] =					TRUE;						//	Boolean_t
		Zone.ShadeColor[GroupIndex] =				Black_C;					//	ColorIndex_t
		Zone.ShadeLighting[GroupIndex] =				TRUE;						//	Boolean_t

		//	Edge Settings

		Zone.EdgeShow[GroupIndex] =					FALSE;						//	Boolean_t
		Zone.EdgeType[GroupIndex] =					EdgeType_BordersAndCreases;	//	EdgeType_e
		Zone.EdgeIBorder[GroupIndex] =				BorderLocation_Both;		//	BorderLocation_e
		Zone.EdgeJBorder[GroupIndex] =				BorderLocation_Both;		//	BorderLocation_e
		Zone.EdgeKBorder[GroupIndex] =				BorderLocation_Both;		//	BorderLocation_e
		Zone.EdgeColor[GroupIndex] =					Black_C;					//	ColorIndex_t
		Zone.EdgeLineThickness[GroupIndex] =			0.1;	//	value is in %	//	double

		//	Points Settings

		Zone.PointsToPlot[GroupIndex] =				PointsToPlot_SurfaceNodes;	//	PointsToPlot_e
		Zone.PointsIndexSkipI[GroupIndex] =			1;							//	LgIndex_t
		Zone.PointsIndexSkipJ[GroupIndex] =			1;							//	LgIndex_t
		Zone.PointsIndexSkipK[GroupIndex] =			1;							//	LgIndex_t
		/*
			1 is the same as "No Skip" in the GUI.
		*/

		//	Surfaces Settings

		Zone.SurfacesToPlot[GroupIndex] =			SurfacesToPlot_ExposedCellFaces;	//	SurfacesToPlot_e
		Zone.RangeBeginI[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeEndI[GroupIndex] =					0;							//	LgIndex_t
		Zone.RangeSkipI[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeBeginJ[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeEndJ[GroupIndex] =					0;							//	LgIndex_t
		Zone.RangeSkipJ[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeBeginK[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeEndK[GroupIndex] =					0;							//	LgIndex_t
		Zone.RangeSkipK[GroupIndex] =				1;							//	LgIndex_t
		/*
			0 is the same as maximum, which sets the end range to the maximum
			index for the dimension. The skip number is to skip every # planes,
			so 1 is to not skip any planes, (show every 1 plane).
		*/

		//	Volume Settings

		Zone.StreamtracesShow[GroupIndex] =			TRUE;						//	Boolean_t
		Zone.IsoSurfacesShow[GroupIndex] =			TRUE;						//	Boolean_t
		Zone.SlicesShow[GroupIndex] =				TRUE;						//	Boolean_t

		//	Effects Settings

		Zone.UseValueBlanking[GroupIndex] =			TRUE;						//	Boolean_t
		//Zone.UseClipPlanes[GroupIndex] =			TRUE;						//	Boolean_t
		Zone.UseSurfaceTranslucency[GroupIndex] =	TRUE;						//	Boolean_t
		Zone.SurfaceTranslucency[GroupIndex] =		50;		//	value is in %	//	SmInteger_t
		Zone.LightingEffect[GroupIndex] =			LightingEffect_Gouraud;		//	LightingEffect_e

		//	Group index and number (increments for next set)

		GroupIndex++;														//	int


		//===============================================================
		//	BRCSurfacesFF Settings
		//===============================================================

		Zone.ZoneShow[GroupIndex] =					FALSE;						//	Boolean_t

		//	Mesh Settings

		Zone.MeshShow[GroupIndex] =					TRUE;						//	Boolean_t
		Zone.MeshType[GroupIndex] =					MeshType_Overlay;			//	MeshType_e
		Zone.MeshColor[GroupIndex] =					Black_C;					//	ColorIndex_t
		/*
			A number of color index constants are \#defined. These include Black_C, Red_C,
			Green_C, Blue_C, Cyan_C, Yellow_C, Purple_C, White_C, Custom1_C through
			Custom56_C, MultiColor_C, NoColor_C, MulitiColor2_C through MulitiColor8_C,
			RGBColor_C, and InvalidColor_C.  
		*/
		Zone.MeshLinePttrn[GroupIndex] =				LinePattern_Solid;			//	LinePattern_e
		Zone.MeshLinePttrnLen[GroupIndex] =			2;		//	Value is in %	//	double
		Zone.MeshLineThickness[GroupIndex] =			0.1;	//	Value is in %	//	double


		//	Contour Settings

		Zone.ContourShow[GroupIndex] =				TRUE;						//	Boolean_t
		Zone.ContourType[GroupIndex] =				ContourType_Flood;			//	ContourType_e
		Zone.FloodColoring[GroupIndex] =				ContourColoring_Group1;		//	ContourColoring_e
		Zone.ContourLineGroup[GroupIndex] =			1;							//	SmInteger_t
		Zone.ContourLineColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.ContourLinePttrn[GroupIndex] =			LinePattern_Solid;			//	LinePattern_e
		Zone.ContourLinePttrnLen[GroupIndex] =		2;		//	Value is in %	//	double
		Zone.ContourLineThickness[GroupIndex] =		0.1;	//	Value is in %	//	double
		Zone.ContourLighting[GroupIndex] =			TRUE;						//	Boolean_t

		//	Vector Settings
		
		Zone.VectorShow[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.VectorType[GroupIndex] =				VectorType_TailAtPoint;		//	VectorType_e
		Zone.VectorArrowType[GroupIndex] =			ArrowheadStyle_Plain;		//	ArrowheadStyle_e
		Zone.VectorLineColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.VectorTangent[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.VectorLinePttrn[GroupIndex] =			LinePattern_Solid;			//	LinePattern_e
		Zone.VectorLinePttrnLen[GroupIndex] =		2;		//	Value is in %	//	double
		Zone.VectorLineThickness[GroupIndex] =		0.1;	//	Value is in %	//	double

		//	Scatter Settings

		Zone.ScatterShow[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.ScatterShape[GroupIndex] =				GeomShape_Sphere;			//	GeomShape_e
		Zone.ScatterOutColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.ScatterInColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.ScatterFillMode[GroupIndex] =			FillMode_UseLineColor;		//	FillMode_e
		Zone.ScatterSize[GroupIndex] =				2.5;	//	value is in %	//	double
		Zone.ScatterLineThickness[GroupIndex] =		0.1;	//	value is in %	//	double

		//	Shade Settings

		Zone.ShadeShow[GroupIndex] =					TRUE;						//	Boolean_t
		Zone.ShadeColor[GroupIndex] =				Black_C;					//	ColorIndex_t
		Zone.ShadeLighting[GroupIndex] =				TRUE;						//	Boolean_t

		//	Edge Settings

		Zone.EdgeShow[GroupIndex] =					FALSE;						//	Boolean_t
		Zone.EdgeType[GroupIndex] =					EdgeType_BordersAndCreases;	//	EdgeType_e
		Zone.EdgeIBorder[GroupIndex] =				BorderLocation_Both;		//	BorderLocation_e
		Zone.EdgeJBorder[GroupIndex] =				BorderLocation_Both;		//	BorderLocation_e
		Zone.EdgeKBorder[GroupIndex] =				BorderLocation_Both;		//	BorderLocation_e
		Zone.EdgeColor[GroupIndex] =					Black_C;					//	ColorIndex_t
		Zone.EdgeLineThickness[GroupIndex] =			0.1;	//	value is in %	//	double

		//	Points Settings

		Zone.PointsToPlot[GroupIndex] =				PointsToPlot_SurfaceNodes;	//	PointsToPlot_e
		Zone.PointsIndexSkipI[GroupIndex] =			1;							//	LgIndex_t
		Zone.PointsIndexSkipJ[GroupIndex] =			1;							//	LgIndex_t
		Zone.PointsIndexSkipK[GroupIndex] =			1;							//	LgIndex_t
		/*
			1 is the same as "No Skip" in the GUI.
		*/

		//	Surfaces Settings

		Zone.SurfacesToPlot[GroupIndex] =			SurfacesToPlot_ExposedCellFaces;	//	SurfacesToPlot_e
		Zone.RangeBeginI[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeEndI[GroupIndex] =					0;							//	LgIndex_t
		Zone.RangeSkipI[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeBeginJ[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeEndJ[GroupIndex] =					0;							//	LgIndex_t
		Zone.RangeSkipJ[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeBeginK[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeEndK[GroupIndex] =					0;							//	LgIndex_t
		Zone.RangeSkipK[GroupIndex] =				1;							//	LgIndex_t
		/*
			0 is the same as maximum, which sets the end range to the maximum
			index for the dimension. The skip number is to skip every # planes,
			so 1 is to not skip any planes, (show every 1 plane).
		*/

		//	Volume Settings

		Zone.StreamtracesShow[GroupIndex] =			TRUE;						//	Boolean_t
		Zone.IsoSurfacesShow[GroupIndex] =			TRUE;						//	Boolean_t
		Zone.SlicesShow[GroupIndex] =				TRUE;						//	Boolean_t

		//	Effects Settings

		Zone.UseValueBlanking[GroupIndex] =			TRUE;						//	Boolean_t
		//Zone.UseClipPlanes[GroupIndex] =			TRUE;						//	Boolean_t
		Zone.UseSurfaceTranslucency[GroupIndex] =	TRUE;						//	Boolean_t
		Zone.SurfaceTranslucency[GroupIndex] =		50;		//	value is in %	//	SmInteger_t
		Zone.LightingEffect[GroupIndex] =			LightingEffect_Gouraud;		//	LightingEffect_e

		//	Group index and number (increments for next set)

		GroupIndex++;														//	int


		//===============================================================
		//	AtomRhoIsoSurfaces Settings
		//===============================================================

		Zone.ZoneShow[GroupIndex] =					TRUE;						//	Boolean_t

		//	Mesh Settings

		Zone.MeshShow[GroupIndex] =					FALSE;						//	Boolean_t
		Zone.MeshType[GroupIndex] =					MeshType_Overlay;			//	MeshType_e
		Zone.MeshColor[GroupIndex] =					MultiColor3_C;					//	ColorIndex_t
		/*
			A number of color index constants are \#defined. These include Black_C, Red_C,
			Green_C, Blue_C, Cyan_C, Yellow_C, Purple_C, White_C, Custom1_C through
			Custom56_C, MultiColor_C, NoColor_C, MulitiColor2_C through MulitiColor8_C,
			RGBColor_C, and InvalidColor_C.  
		*/
		Zone.MeshLinePttrn[GroupIndex] =				LinePattern_Solid;			//	LinePattern_e
		Zone.MeshLinePttrnLen[GroupIndex] =			2;		//	Value is in %	//	double
		Zone.MeshLineThickness[GroupIndex] =			0.1;	//	Value is in %	//	double


		//	Contour Settings

		Zone.ContourShow[GroupIndex] =				TRUE;						//	Boolean_t
		Zone.ContourType[GroupIndex] =				ContourType_Flood;			//	ContourType_e
		Zone.FloodColoring[GroupIndex] =				ContourColoring_Group3;		//	ContourColoring_e
		Zone.ContourLineGroup[GroupIndex] =			3;							//	SmInteger_t
		Zone.ContourLineColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.ContourLinePttrn[GroupIndex] =			LinePattern_Solid;			//	LinePattern_e
		Zone.ContourLinePttrnLen[GroupIndex] =		2;		//	Value is in %	//	double
		Zone.ContourLineThickness[GroupIndex] =		0.1;	//	Value is in %	//	double
		Zone.ContourLighting[GroupIndex] =			TRUE;						//	Boolean_t

		//	Vector Settings
		
		Zone.VectorShow[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.VectorType[GroupIndex] =				VectorType_TailAtPoint;		//	VectorType_e
		Zone.VectorArrowType[GroupIndex] =			ArrowheadStyle_Plain;		//	ArrowheadStyle_e
		Zone.VectorLineColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.VectorTangent[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.VectorLinePttrn[GroupIndex] =			LinePattern_Solid;			//	LinePattern_e
		Zone.VectorLinePttrnLen[GroupIndex] =		2;		//	Value is in %	//	double
		Zone.VectorLineThickness[GroupIndex] =		0.1;	//	Value is in %	//	double

		//	Scatter Settings

		Zone.ScatterShow[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.ScatterShape[GroupIndex] =				GeomShape_Sphere;			//	GeomShape_e
		Zone.ScatterOutColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.ScatterInColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.ScatterFillMode[GroupIndex] =			FillMode_UseLineColor;		//	FillMode_e
		Zone.ScatterSize[GroupIndex] =				2.5;	//	value is in %	//	double
		Zone.ScatterLineThickness[GroupIndex] =		0.1;	//	value is in %	//	double

		//	Shade Settings

		Zone.ShadeShow[GroupIndex] =					FALSE;						//	Boolean_t
		Zone.ShadeColor[GroupIndex] =				Black_C;					//	ColorIndex_t
		Zone.ShadeLighting[GroupIndex] =				TRUE;						//	Boolean_t

		//	Edge Settings

		Zone.EdgeShow[GroupIndex] =					FALSE;						//	Boolean_t
		Zone.EdgeType[GroupIndex] =					EdgeType_BordersAndCreases;	//	EdgeType_e
		Zone.EdgeIBorder[GroupIndex] =				BorderLocation_Both;		//	BorderLocation_e
		Zone.EdgeJBorder[GroupIndex] =				BorderLocation_Both;		//	BorderLocation_e
		Zone.EdgeKBorder[GroupIndex] =				BorderLocation_Both;		//	BorderLocation_e
		Zone.EdgeColor[GroupIndex] =					Black_C;					//	ColorIndex_t
		Zone.EdgeLineThickness[GroupIndex] =			0.1;	//	value is in %	//	double

		//	Points Settings

		Zone.PointsToPlot[GroupIndex] =				PointsToPlot_SurfaceNodes;	//	PointsToPlot_e
		Zone.PointsIndexSkipI[GroupIndex] =			1;							//	LgIndex_t
		Zone.PointsIndexSkipJ[GroupIndex] =			1;							//	LgIndex_t
		Zone.PointsIndexSkipK[GroupIndex] =			1;							//	LgIndex_t
		/*
			1 is the same as "No Skip" in the GUI.
		*/

		//	Surfaces Settings

		Zone.SurfacesToPlot[GroupIndex] =			SurfacesToPlot_ExposedCellFaces;	//	SurfacesToPlot_e
		Zone.RangeBeginI[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeEndI[GroupIndex] =					0;							//	LgIndex_t
		Zone.RangeSkipI[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeBeginJ[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeEndJ[GroupIndex] =					0;							//	LgIndex_t
		Zone.RangeSkipJ[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeBeginK[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeEndK[GroupIndex] =					0;							//	LgIndex_t
		Zone.RangeSkipK[GroupIndex] =				1;							//	LgIndex_t
		/*
			0 is the same as maximum, which sets the end range to the maximum
			index for the dimension. The skip number is to skip every # planes,
			so 1 is to not skip any planes, (show every 1 plane).
		*/

		//	Volume Settings

		Zone.StreamtracesShow[GroupIndex] =			TRUE;						//	Boolean_t
		Zone.IsoSurfacesShow[GroupIndex] =			TRUE;						//	Boolean_t
		Zone.SlicesShow[GroupIndex] =				TRUE;						//	Boolean_t

		//	Effects Settings

		Zone.UseValueBlanking[GroupIndex] =			FALSE;						//	Boolean_t
		//Zone.UseClipPlanes[GroupIndex] =			TRUE;						//	Boolean_t
		Zone.UseSurfaceTranslucency[GroupIndex] =	TRUE;						//	Boolean_t
		Zone.SurfaceTranslucency[GroupIndex] =		50;		//	value is in %	//	SmInteger_t
		Zone.LightingEffect[GroupIndex] =			LightingEffect_Gouraud;		//	LightingEffect_e

		//	Group index and number (increments for next set)

		GroupIndex++;														//	int


		//===============================================================
		//	ExtraIsoSurfaces Settings
		//===============================================================

		Zone.ZoneShow[GroupIndex] =					FALSE;						//	Boolean_t

		//	Mesh Settings

		Zone.MeshShow[GroupIndex] =					TRUE;						//	Boolean_t
		Zone.MeshType[GroupIndex] =					MeshType_Overlay;			//	MeshType_e
		Zone.MeshColor[GroupIndex] =					Black_C;					//	ColorIndex_t
		/*
			A number of color index constants are \#defined. These include Black_C, Red_C,
			Green_C, Blue_C, Cyan_C, Yellow_C, Purple_C, White_C, Custom1_C through
			Custom56_C, MultiColor_C, NoColor_C, MulitiColor2_C through MulitiColor8_C,
			RGBColor_C, and InvalidColor_C.  
		*/
		Zone.MeshLinePttrn[GroupIndex] =				LinePattern_Solid;			//	LinePattern_e
		Zone.MeshLinePttrnLen[GroupIndex] =			2;		//	Value is in %	//	double
		Zone.MeshLineThickness[GroupIndex] =			0.1;	//	Value is in %	//	double


		//	Contour Settings

		Zone.ContourShow[GroupIndex] =				TRUE;						//	Boolean_t
		Zone.ContourType[GroupIndex] =				ContourType_Flood;			//	ContourType_e
		Zone.FloodColoring[GroupIndex] =				ContourColoring_Group1;		//	ContourColoring_e
		Zone.ContourLineGroup[GroupIndex] =			1;							//	SmInteger_t
		Zone.ContourLineColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.ContourLinePttrn[GroupIndex] =			LinePattern_Solid;			//	LinePattern_e
		Zone.ContourLinePttrnLen[GroupIndex] =		2;		//	Value is in %	//	double
		Zone.ContourLineThickness[GroupIndex] =		0.1;	//	Value is in %	//	double
		Zone.ContourLighting[GroupIndex] =			TRUE;						//	Boolean_t

		//	Vector Settings
		
		Zone.VectorShow[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.VectorType[GroupIndex] =				VectorType_TailAtPoint;		//	VectorType_e
		Zone.VectorArrowType[GroupIndex] =			ArrowheadStyle_Plain;		//	ArrowheadStyle_e
		Zone.VectorLineColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.VectorTangent[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.VectorLinePttrn[GroupIndex] =			LinePattern_Solid;			//	LinePattern_e
		Zone.VectorLinePttrnLen[GroupIndex] =		2;		//	Value is in %	//	double
		Zone.VectorLineThickness[GroupIndex] =		0.1;	//	Value is in %	//	double

		//	Scatter Settings

		Zone.ScatterShow[GroupIndex] =				TRUE;						//	Boolean_t
		Zone.ScatterShape[GroupIndex] =				GeomShape_Sphere;			//	GeomShape_e
		Zone.ScatterOutColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.ScatterInColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.ScatterFillMode[GroupIndex] =			FillMode_UseLineColor;		//	FillMode_e
		Zone.ScatterSize[GroupIndex] =				2.5;	//	value is in %	//	double
		Zone.ScatterLineThickness[GroupIndex] =		0.1;	//	value is in %	//	double

		//	Shade Settings

		Zone.ShadeShow[GroupIndex] =					TRUE;						//	Boolean_t
		Zone.ShadeColor[GroupIndex] =				Black_C;					//	ColorIndex_t
		Zone.ShadeLighting[GroupIndex] =				TRUE;						//	Boolean_t

		//	Edge Settings

		Zone.EdgeShow[GroupIndex] =					TRUE;						//	Boolean_t
		Zone.EdgeType[GroupIndex] =					EdgeType_BordersAndCreases;	//	EdgeType_e
		Zone.EdgeIBorder[GroupIndex] =				BorderLocation_Both;		//	BorderLocation_e
		Zone.EdgeJBorder[GroupIndex] =				BorderLocation_Both;		//	BorderLocation_e
		Zone.EdgeKBorder[GroupIndex] =				BorderLocation_Both;		//	BorderLocation_e
		Zone.EdgeColor[GroupIndex] =					Black_C;					//	ColorIndex_t
		Zone.EdgeLineThickness[GroupIndex] =			0.1;	//	value is in %	//	double

		//	Points Settings

		Zone.PointsToPlot[GroupIndex] =				PointsToPlot_SurfaceNodes;	//	PointsToPlot_e
		Zone.PointsIndexSkipI[GroupIndex] =			1;							//	LgIndex_t
		Zone.PointsIndexSkipJ[GroupIndex] =			1;							//	LgIndex_t
		Zone.PointsIndexSkipK[GroupIndex] =			1;							//	LgIndex_t
		/*
			1 is the same as "No Skip" in the GUI.
		*/

		//	Surfaces Settings

		Zone.SurfacesToPlot[GroupIndex] =			SurfacesToPlot_ExposedCellFaces;	//	SurfacesToPlot_e
		Zone.RangeBeginI[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeEndI[GroupIndex] =					0;							//	LgIndex_t
		Zone.RangeSkipI[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeBeginJ[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeEndJ[GroupIndex] =					0;							//	LgIndex_t
		Zone.RangeSkipJ[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeBeginK[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeEndK[GroupIndex] =					0;							//	LgIndex_t
		Zone.RangeSkipK[GroupIndex] =				1;							//	LgIndex_t
		/*
			0 is the same as maximum, which sets the end range to the maximum
			index for the dimension. The skip number is to skip every # planes,
			so 1 is to not skip any planes, (show every 1 plane).
		*/

		//	Volume Settings

		Zone.StreamtracesShow[GroupIndex] =			TRUE;						//	Boolean_t
		Zone.IsoSurfacesShow[GroupIndex] =			TRUE;						//	Boolean_t
		Zone.SlicesShow[GroupIndex] =				TRUE;						//	Boolean_t

		//	Effects Settings

		Zone.UseValueBlanking[GroupIndex] =			TRUE;						//	Boolean_t
		//Zone.UseClipPlanes[GroupIndex] =			TRUE;						//	Boolean_t
		Zone.UseSurfaceTranslucency[GroupIndex] =	TRUE;						//	Boolean_t
		Zone.SurfaceTranslucency[GroupIndex] =		50;		//	value is in %	//	SmInteger_t
		Zone.LightingEffect[GroupIndex] =			LightingEffect_Gouraud;		//	LightingEffect_e

		//	Group index and number (increments for next set)

		GroupIndex++;														//	int


		//===============================================================
		//	BondZones Settings
		//===============================================================

		Zone.ZoneShow[GroupIndex] =					FALSE;						//	Boolean_t

		//	Mesh Settings

		Zone.MeshShow[GroupIndex] =					TRUE;						//	Boolean_t
		Zone.MeshType[GroupIndex] =					MeshType_Overlay;			//	MeshType_e
		Zone.MeshColor[GroupIndex] =					Black_C;					//	ColorIndex_t
		/*
			A number of color index constants are \#defined. These include Black_C, Red_C,
			Green_C, Blue_C, Cyan_C, Yellow_C, Purple_C, White_C, Custom1_C through
			Custom56_C, MultiColor_C, NoColor_C, MulitiColor2_C through MulitiColor8_C,
			RGBColor_C, and InvalidColor_C.  
		*/
		Zone.MeshLinePttrn[GroupIndex] =				LinePattern_Solid;			//	LinePattern_e
		Zone.MeshLinePttrnLen[GroupIndex] =			2;		//	Value is in %	//	double
		Zone.MeshLineThickness[GroupIndex] =			0.1;	//	Value is in %	//	double


		//	Contour Settings

		Zone.ContourShow[GroupIndex] =				TRUE;						//	Boolean_t
		Zone.ContourType[GroupIndex] =				ContourType_Flood;			//	ContourType_e
		Zone.FloodColoring[GroupIndex] =				ContourColoring_Group1;		//	ContourColoring_e
		Zone.ContourLineGroup[GroupIndex] =			1;							//	SmInteger_t
		Zone.ContourLineColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.ContourLinePttrn[GroupIndex] =			LinePattern_Solid;			//	LinePattern_e
		Zone.ContourLinePttrnLen[GroupIndex] =		2;		//	Value is in %	//	double
		Zone.ContourLineThickness[GroupIndex] =		0.1;	//	Value is in %	//	double
		Zone.ContourLighting[GroupIndex] =			TRUE;						//	Boolean_t

		//	Vector Settings
		
		Zone.VectorShow[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.VectorType[GroupIndex] =				VectorType_TailAtPoint;		//	VectorType_e
		Zone.VectorArrowType[GroupIndex] =			ArrowheadStyle_Plain;		//	ArrowheadStyle_e
		Zone.VectorLineColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.VectorTangent[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.VectorLinePttrn[GroupIndex] =			LinePattern_Solid;			//	LinePattern_e
		Zone.VectorLinePttrnLen[GroupIndex] =		2;		//	Value is in %	//	double
		Zone.VectorLineThickness[GroupIndex] =		0.1;	//	Value is in %	//	double

		//	Scatter Settings

		Zone.ScatterShow[GroupIndex] =				TRUE;						//	Boolean_t
		Zone.ScatterShape[GroupIndex] =				GeomShape_Sphere;			//	GeomShape_e
		Zone.ScatterOutColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.ScatterInColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.ScatterFillMode[GroupIndex] =			FillMode_UseLineColor;		//	FillMode_e
		Zone.ScatterSize[GroupIndex] =				2.5;	//	value is in %	//	double
		Zone.ScatterLineThickness[GroupIndex] =		0.1;	//	value is in %	//	double

		//	Shade Settings

		Zone.ShadeShow[GroupIndex] =					TRUE;						//	Boolean_t
		Zone.ShadeColor[GroupIndex] =				Black_C;					//	ColorIndex_t
		Zone.ShadeLighting[GroupIndex] =				TRUE;						//	Boolean_t

		//	Edge Settings

		Zone.EdgeShow[GroupIndex] =					TRUE;						//	Boolean_t
		Zone.EdgeType[GroupIndex] =					EdgeType_BordersAndCreases;	//	EdgeType_e
		Zone.EdgeIBorder[GroupIndex] =				BorderLocation_Both;		//	BorderLocation_e
		Zone.EdgeJBorder[GroupIndex] =				BorderLocation_Both;		//	BorderLocation_e
		Zone.EdgeKBorder[GroupIndex] =				BorderLocation_Both;		//	BorderLocation_e
		Zone.EdgeColor[GroupIndex] =					Black_C;					//	ColorIndex_t
		Zone.EdgeLineThickness[GroupIndex] =			0.1;	//	value is in %	//	double

		//	Points Settings

		Zone.PointsToPlot[GroupIndex] =				PointsToPlot_SurfaceNodes;	//	PointsToPlot_e
		Zone.PointsIndexSkipI[GroupIndex] =			1;							//	LgIndex_t
		Zone.PointsIndexSkipJ[GroupIndex] =			1;							//	LgIndex_t
		Zone.PointsIndexSkipK[GroupIndex] =			1;							//	LgIndex_t
		/*
			1 is the same as "No Skip" in the GUI.
		*/

		//	Surfaces Settings

		Zone.SurfacesToPlot[GroupIndex] =			SurfacesToPlot_ExposedCellFaces;	//	SurfacesToPlot_e
		Zone.RangeBeginI[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeEndI[GroupIndex] =					0;							//	LgIndex_t
		Zone.RangeSkipI[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeBeginJ[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeEndJ[GroupIndex] =					0;							//	LgIndex_t
		Zone.RangeSkipJ[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeBeginK[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeEndK[GroupIndex] =					0;							//	LgIndex_t
		Zone.RangeSkipK[GroupIndex] =				1;							//	LgIndex_t
		/*
			0 is the same as maximum, which sets the end range to the maximum
			index for the dimension. The skip number is to skip every # planes,
			so 1 is to not skip any planes, (show every 1 plane).
		*/

		//	Volume Settings

		Zone.StreamtracesShow[GroupIndex] =			TRUE;						//	Boolean_t
		Zone.IsoSurfacesShow[GroupIndex] =			TRUE;						//	Boolean_t
		Zone.SlicesShow[GroupIndex] =				TRUE;						//	Boolean_t

		//	Effects Settings

		Zone.UseValueBlanking[GroupIndex] =			TRUE;						//	Boolean_t
		//Zone.UseClipPlanes[GroupIndex] =			TRUE;						//	Boolean_t
		Zone.UseSurfaceTranslucency[GroupIndex] =	TRUE;						//	Boolean_t
		Zone.SurfaceTranslucency[GroupIndex] =		50;		//	value is in %	//	SmInteger_t
		Zone.LightingEffect[GroupIndex] =			LightingEffect_Gouraud;		//	LightingEffect_e

		//	Group index and number (increments for next set)

		GroupIndex++;														//	int


		//===============================================================
		//	BondZonesFF Settings
		//===============================================================

		Zone.ZoneShow[GroupIndex] =					FALSE;						//	Boolean_t

		//	Mesh Settings

		Zone.MeshShow[GroupIndex] =					TRUE;						//	Boolean_t
		Zone.MeshType[GroupIndex] =					MeshType_Overlay;			//	MeshType_e
		Zone.MeshColor[GroupIndex] =					Black_C;					//	ColorIndex_t
		/*
			A number of color index constants are \#defined. These include Black_C, Red_C,
			Green_C, Blue_C, Cyan_C, Yellow_C, Purple_C, White_C, Custom1_C through
			Custom56_C, MultiColor_C, NoColor_C, MulitiColor2_C through MulitiColor8_C,
			RGBColor_C, and InvalidColor_C.  
		*/
		Zone.MeshLinePttrn[GroupIndex] =				LinePattern_Solid;			//	LinePattern_e
		Zone.MeshLinePttrnLen[GroupIndex] =			2;		//	Value is in %	//	double
		Zone.MeshLineThickness[GroupIndex] =			0.1;	//	Value is in %	//	double


		//	Contour Settings

		Zone.ContourShow[GroupIndex] =				TRUE;						//	Boolean_t
		Zone.ContourType[GroupIndex] =				ContourType_Flood;			//	ContourType_e
		Zone.FloodColoring[GroupIndex] =				ContourColoring_Group1;		//	ContourColoring_e
		Zone.ContourLineGroup[GroupIndex] =			1;							//	SmInteger_t
		Zone.ContourLineColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.ContourLinePttrn[GroupIndex] =			LinePattern_Solid;			//	LinePattern_e
		Zone.ContourLinePttrnLen[GroupIndex] =		2;		//	Value is in %	//	double
		Zone.ContourLineThickness[GroupIndex] =		0.1;	//	Value is in %	//	double
		Zone.ContourLighting[GroupIndex] =			TRUE;						//	Boolean_t

		//	Vector Settings
		
		Zone.VectorShow[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.VectorType[GroupIndex] =				VectorType_TailAtPoint;		//	VectorType_e
		Zone.VectorArrowType[GroupIndex] =			ArrowheadStyle_Plain;		//	ArrowheadStyle_e
		Zone.VectorLineColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.VectorTangent[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.VectorLinePttrn[GroupIndex] =			LinePattern_Solid;			//	LinePattern_e
		Zone.VectorLinePttrnLen[GroupIndex] =		2;		//	Value is in %	//	double
		Zone.VectorLineThickness[GroupIndex] =		0.1;	//	Value is in %	//	double

		//	Scatter Settings

		Zone.ScatterShow[GroupIndex] =				TRUE;						//	Boolean_t
		Zone.ScatterShape[GroupIndex] =				GeomShape_Sphere;			//	GeomShape_e
		Zone.ScatterOutColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.ScatterInColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.ScatterFillMode[GroupIndex] =			FillMode_UseLineColor;		//	FillMode_e
		Zone.ScatterSize[GroupIndex] =				2.5;	//	value is in %	//	double
		Zone.ScatterLineThickness[GroupIndex] =		0.1;	//	value is in %	//	double

		//	Shade Settings

		Zone.ShadeShow[GroupIndex] =					TRUE;						//	Boolean_t
		Zone.ShadeColor[GroupIndex] =				Black_C;					//	ColorIndex_t
		Zone.ShadeLighting[GroupIndex] =				TRUE;						//	Boolean_t

		//	Edge Settings

		Zone.EdgeShow[GroupIndex] =					TRUE;						//	Boolean_t
		Zone.EdgeType[GroupIndex] =					EdgeType_BordersAndCreases;	//	EdgeType_e
		Zone.EdgeIBorder[GroupIndex] =				BorderLocation_Both;		//	BorderLocation_e
		Zone.EdgeJBorder[GroupIndex] =				BorderLocation_Both;		//	BorderLocation_e
		Zone.EdgeKBorder[GroupIndex] =				BorderLocation_Both;		//	BorderLocation_e
		Zone.EdgeColor[GroupIndex] =					Black_C;					//	ColorIndex_t
		Zone.EdgeLineThickness[GroupIndex] =			0.1;	//	value is in %	//	double

		//	Points Settings

		Zone.PointsToPlot[GroupIndex] =				PointsToPlot_SurfaceNodes;	//	PointsToPlot_e
		Zone.PointsIndexSkipI[GroupIndex] =			1;							//	LgIndex_t
		Zone.PointsIndexSkipJ[GroupIndex] =			1;							//	LgIndex_t
		Zone.PointsIndexSkipK[GroupIndex] =			1;							//	LgIndex_t
		/*
			1 is the same as "No Skip" in the GUI.
		*/

		//	Surfaces Settings

		Zone.SurfacesToPlot[GroupIndex] =			SurfacesToPlot_ExposedCellFaces;	//	SurfacesToPlot_e
		Zone.RangeBeginI[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeEndI[GroupIndex] =					0;							//	LgIndex_t
		Zone.RangeSkipI[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeBeginJ[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeEndJ[GroupIndex] =					0;							//	LgIndex_t
		Zone.RangeSkipJ[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeBeginK[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeEndK[GroupIndex] =					0;							//	LgIndex_t
		Zone.RangeSkipK[GroupIndex] =				1;							//	LgIndex_t
		/*
			0 is the same as maximum, which sets the end range to the maximum
			index for the dimension. The skip number is to skip every # planes,
			so 1 is to not skip any planes, (show every 1 plane).
		*/

		//	Volume Settings

		Zone.StreamtracesShow[GroupIndex] =			TRUE;						//	Boolean_t
		Zone.IsoSurfacesShow[GroupIndex] =			TRUE;						//	Boolean_t
		Zone.SlicesShow[GroupIndex] =				TRUE;						//	Boolean_t

		//	Effects Settings

		Zone.UseValueBlanking[GroupIndex] =			TRUE;						//	Boolean_t
		//Zone.UseClipPlanes[GroupIndex] =			TRUE;						//	Boolean_t
		Zone.UseSurfaceTranslucency[GroupIndex] =	TRUE;						//	Boolean_t
		Zone.SurfaceTranslucency[GroupIndex] =		50;		//	value is in %	//	SmInteger_t
		Zone.LightingEffect[GroupIndex] =			LightingEffect_Gouraud;		//	LightingEffect_e

		//	Group index and number (increments for next set)

		GroupIndex++;														//	int

		//===============================================================
		//	BadZones Settings
		//===============================================================

		Zone.ZoneShow[GroupIndex] =					FALSE;						//	Boolean_t

		//	Mesh Settings

		Zone.MeshShow[GroupIndex] =					FALSE;						//	Boolean_t
		Zone.MeshType[GroupIndex] =					MeshType_Overlay;			//	MeshType_e
		Zone.MeshColor[GroupIndex] =					Black_C;					//	ColorIndex_t
		/*
			A number of color index constants are \#defined. These include Black_C, Red_C,
			Green_C, Blue_C, Cyan_C, Yellow_C, Purple_C, White_C, Custom1_C through
			Custom56_C, MultiColor_C, NoColor_C, MulitiColor2_C through MulitiColor8_C,
			RGBColor_C, and InvalidColor_C.  
		*/
		Zone.MeshLinePttrn[GroupIndex] =				LinePattern_Solid;			//	LinePattern_e
		Zone.MeshLinePttrnLen[GroupIndex] =			2;		//	Value is in %	//	double
		Zone.MeshLineThickness[GroupIndex] =			0.1;	//	Value is in %	//	double


		//	Contour Settings

		Zone.ContourShow[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.ContourType[GroupIndex] =				ContourType_Flood;			//	ContourType_e
		Zone.FloodColoring[GroupIndex] =				ContourColoring_Group1;		//	ContourColoring_e
		Zone.ContourLineGroup[GroupIndex] =			1;							//	SmInteger_t
		Zone.ContourLineColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.ContourLinePttrn[GroupIndex] =			LinePattern_Solid;			//	LinePattern_e
		Zone.ContourLinePttrnLen[GroupIndex] =		2;		//	Value is in %	//	double
		Zone.ContourLineThickness[GroupIndex] =		0.1;	//	Value is in %	//	double
		Zone.ContourLighting[GroupIndex] =			TRUE;						//	Boolean_t

		//	Vector Settings
		
		Zone.VectorShow[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.VectorType[GroupIndex] =				VectorType_TailAtPoint;		//	VectorType_e
		Zone.VectorArrowType[GroupIndex] =			ArrowheadStyle_Plain;		//	ArrowheadStyle_e
		Zone.VectorLineColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.VectorTangent[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.VectorLinePttrn[GroupIndex] =			LinePattern_Solid;			//	LinePattern_e
		Zone.VectorLinePttrnLen[GroupIndex] =		2;		//	Value is in %	//	double
		Zone.VectorLineThickness[GroupIndex] =		0.1;	//	Value is in %	//	double

		//	Scatter Settings

		Zone.ScatterShow[GroupIndex] =				FALSE;						//	Boolean_t
		Zone.ScatterShape[GroupIndex] =				GeomShape_Sphere;			//	GeomShape_e
		Zone.ScatterOutColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.ScatterInColor[GroupIndex] =			Black_C;					//	ColorIndex_t
		Zone.ScatterFillMode[GroupIndex] =			FillMode_UseLineColor;		//	FillMode_e
		Zone.ScatterSize[GroupIndex] =				2.5;	//	value is in %	//	double
		Zone.ScatterLineThickness[GroupIndex] =		0.1;	//	value is in %	//	double

		//	Shade Settings

		Zone.ShadeShow[GroupIndex] =					FALSE;						//	Boolean_t
		Zone.ShadeColor[GroupIndex] =				Black_C;					//	ColorIndex_t
		Zone.ShadeLighting[GroupIndex] =				TRUE;						//	Boolean_t

		//	Edge Settings

		Zone.EdgeShow[GroupIndex] =					FALSE;						//	Boolean_t
		Zone.EdgeType[GroupIndex] =					EdgeType_BordersAndCreases;	//	EdgeType_e
		Zone.EdgeIBorder[GroupIndex] =				BorderLocation_Both;		//	BorderLocation_e
		Zone.EdgeJBorder[GroupIndex] =				BorderLocation_Both;		//	BorderLocation_e
		Zone.EdgeKBorder[GroupIndex] =				BorderLocation_Both;		//	BorderLocation_e
		Zone.EdgeColor[GroupIndex] =					Black_C;					//	ColorIndex_t
		Zone.EdgeLineThickness[GroupIndex] =			0.1;	//	value is in %	//	double

		//	Points Settings

		Zone.PointsToPlot[GroupIndex] =				PointsToPlot_SurfaceNodes;	//	PointsToPlot_e
		Zone.PointsIndexSkipI[GroupIndex] =			1;							//	LgIndex_t
		Zone.PointsIndexSkipJ[GroupIndex] =			1;							//	LgIndex_t
		Zone.PointsIndexSkipK[GroupIndex] =			1;							//	LgIndex_t
		/*
			1 is the same as "No Skip" in the GUI.
		*/

		//	Surfaces Settings

		Zone.SurfacesToPlot[GroupIndex] =			SurfacesToPlot_ExposedCellFaces;	//	SurfacesToPlot_e
		Zone.RangeBeginI[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeEndI[GroupIndex] =					0;							//	LgIndex_t
		Zone.RangeSkipI[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeBeginJ[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeEndJ[GroupIndex] =					0;							//	LgIndex_t
		Zone.RangeSkipJ[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeBeginK[GroupIndex] =				1;							//	LgIndex_t
		Zone.RangeEndK[GroupIndex] =					0;							//	LgIndex_t
		Zone.RangeSkipK[GroupIndex] =				1;							//	LgIndex_t
		/*
			0 is the same as maximum, which sets the end range to the maximum
			index for the dimension. The skip number is to skip every # planes,
			so 1 is to not skip any planes, (show every 1 plane).
		*/

		//	Volume Settings

		Zone.StreamtracesShow[GroupIndex] =			TRUE;						//	Boolean_t
		Zone.IsoSurfacesShow[GroupIndex] =			TRUE;						//	Boolean_t
		Zone.SlicesShow[GroupIndex] =				TRUE;						//	Boolean_t

		//	Effects Settings

		Zone.UseValueBlanking[GroupIndex] =			TRUE;						//	Boolean_t
		//Zone.UseClipPlanes[GroupIndex] =			TRUE;						//	Boolean_t
		Zone.UseSurfaceTranslucency[GroupIndex] =	TRUE;						//	Boolean_t
		Zone.SurfaceTranslucency[GroupIndex] =		50;		//	value is in %	//	SmInteger_t
		Zone.LightingEffect[GroupIndex] =			LightingEffect_Gouraud;		//	LightingEffect_e

		//	Group index and number (increments for next set)

		GroupIndex++;														//	int


		ApplyZoneSettings(ZoneSets, Zone);

		DeallocSets(ZoneSets);



		//char FinishedMsg[200];
		//sprintf_s(FinishedMsg, "Finished Successfully");
		//TecUtilDialogMessageBox(FinishedMsg, MessageBox_Information);
	}
	else
	{
		TecUtilDialogErrMsg("Data set not available.");
	}

} //	void SGradTopoView()



void LogDataValues(std::vector<double> &LevelVector, double Min, double Max, double BaseValue)
{
	if (Max - Min <= 0 && Min != 1e10)
	{
		LevelVector.push_back(0.0);
	}
	else
	{
		int LevelInt = 0, MultInt = 1, ExpInt = 0;
		double TempValue = 0;
		bool BreakLoop = false;
		TempValue = 1.0;
		if (Min != 1e10)
		{
			while (!BreakLoop)
			{
				if (TempValue > Min) TempValue *= 0.1;
				else if (TempValue * 10 <= Min) TempValue *= 10;
				else 
				{
					BaseValue = TempValue;
					BreakLoop = true;
				}
			}
		}
		BreakLoop = false;
		while (!BreakLoop)
		{
			TempValue = MultInt * BaseValue * pow(10.0,ExpInt);
			if (TempValue >= Max)
			{
				TempValue = Max;
				BreakLoop = true;
			}
			LevelVector.push_back(TempValue);

			if (MultInt >= 9)
			{
				MultInt = 1;
				ExpInt++;
			}
			else MultInt++;
		}
	}
}	//	void LogDataValues(std::vector<double> &LevelVector, double Min, double Max, double BaseValue)


void ApplyZoneSettings(std::vector<Set_pa> const &ZoneSets, ZoneSettings const &Zone)
{
	//***************************************************************
	//===============================================================
	//	Loop over array of group sets to apply settings
	//===============================================================
	//***************************************************************

	const size_t NumOfSets = ZoneSets.size();
	int			GroupIndex = 0;

	//	Temporary object set and arg list
	Set_pa CurrentObjectSet = TecUtilSetAlloc(TRUE);
	ArgList_pa CurrentArgList = TecUtilArgListAlloc();

	//===============================================================
	//	Activate only layers that are desired
	//===============================================================

	//	Mesh on
	TecUtilArgListClear(CurrentArgList);
	TecUtilArgListAppendString(CurrentArgList, SV_P1, SV_FIELDLAYERS);
	TecUtilArgListAppendString(CurrentArgList, SV_P2, SV_SHOWMESH);
	TecUtilArgListAppendArbParam(CurrentArgList, SV_IVALUE, TRUE);
	TecUtilStyleSetLowLevelX(CurrentArgList);

	//	Contour on
	TecUtilArgListClear(CurrentArgList);
	TecUtilArgListAppendString(CurrentArgList, SV_P1, SV_FIELDLAYERS);
	TecUtilArgListAppendString(CurrentArgList, SV_P2, SV_SHOWCONTOUR);
	TecUtilArgListAppendArbParam(CurrentArgList, SV_IVALUE, TRUE);
	TecUtilStyleSetLowLevelX(CurrentArgList);

	//	Vector off
	TecUtilArgListClear(CurrentArgList);
	TecUtilArgListAppendString(CurrentArgList, SV_P1, SV_FIELDLAYERS);
	TecUtilArgListAppendString(CurrentArgList, SV_P2, SV_SHOWVECTOR);
	TecUtilArgListAppendArbParam(CurrentArgList, SV_IVALUE, FALSE);
	TecUtilStyleSetLowLevelX(CurrentArgList);

	//	Scatter on
	TecUtilArgListClear(CurrentArgList);
	TecUtilArgListAppendString(CurrentArgList, SV_P1, SV_FIELDLAYERS);
	TecUtilArgListAppendString(CurrentArgList, SV_P2, SV_SHOWSCATTER);
	TecUtilArgListAppendArbParam(CurrentArgList, SV_IVALUE, TRUE);
	TecUtilStyleSetLowLevelX(CurrentArgList);

	//	Shade off
	TecUtilArgListClear(CurrentArgList);
	TecUtilArgListAppendString(CurrentArgList, SV_P1, SV_FIELDLAYERS);
	TecUtilArgListAppendString(CurrentArgList, SV_P2, SV_SHOWSHADE);
	TecUtilArgListAppendArbParam(CurrentArgList, SV_IVALUE, FALSE);
	TecUtilStyleSetLowLevelX(CurrentArgList);

	//	Edge off
	TecUtilArgListClear(CurrentArgList);
	TecUtilArgListAppendString(CurrentArgList, SV_P1, SV_FIELDLAYERS);
	TecUtilArgListAppendString(CurrentArgList, SV_P2, SV_SHOWEDGE);
	TecUtilArgListAppendArbParam(CurrentArgList, SV_IVALUE, FALSE);
	TecUtilStyleSetLowLevelX(CurrentArgList);

	//	Mesh on
	TecUtilArgListClear(CurrentArgList);
	TecUtilArgListAppendString(CurrentArgList, SV_P1, SV_FIELDLAYERS);
	TecUtilArgListAppendString(CurrentArgList, SV_P2, SV_SHOWMESH);
	TecUtilArgListAppendArbParam(CurrentArgList, SV_IVALUE, TRUE);
	TecUtilStyleSetLowLevelX(CurrentArgList);

	//	Lighting on
	TecUtilArgListClear(CurrentArgList);
	TecUtilArgListAppendString(CurrentArgList, SV_P1, SV_FIELDLAYERS);
	TecUtilArgListAppendString(CurrentArgList, SV_P2, SV_USELIGHTINGEFFECT);
	TecUtilArgListAppendArbParam(CurrentArgList, SV_IVALUE, TRUE);
	TecUtilStyleSetLowLevelX(CurrentArgList);

	//	Translucency on
	TecUtilArgListClear(CurrentArgList);
	TecUtilArgListAppendString(CurrentArgList, SV_P1, SV_FIELDLAYERS);
	TecUtilArgListAppendString(CurrentArgList, SV_P2, SV_USETRANSLUCENCY);
	TecUtilArgListAppendArbParam(CurrentArgList, SV_IVALUE, TRUE);
	TecUtilStyleSetLowLevelX(CurrentArgList);

	//	Isosurfaces off
	TecUtilArgListClear(CurrentArgList);
	TecUtilArgListAppendString(CurrentArgList, SV_P1, SV_FIELDLAYERS);
	TecUtilArgListAppendString(CurrentArgList, SV_P2, SV_SHOWISOSURFACES);
	TecUtilArgListAppendArbParam(CurrentArgList, SV_IVALUE, TRUE);
	TecUtilStyleSetLowLevelX(CurrentArgList);

	//	Slices off
	TecUtilArgListClear(CurrentArgList);
	TecUtilArgListAppendString(CurrentArgList, SV_P1, SV_FIELDLAYERS);
	TecUtilArgListAppendString(CurrentArgList, SV_P2, SV_SHOWSLICES);
	TecUtilArgListAppendArbParam(CurrentArgList, SV_IVALUE, FALSE);
	TecUtilStyleSetLowLevelX(CurrentArgList);

	//	Streamtraces off
	TecUtilArgListClear(CurrentArgList);
	TecUtilArgListAppendString(CurrentArgList, SV_P1, SV_FIELDLAYERS);
	TecUtilArgListAppendString(CurrentArgList, SV_P2, SV_SHOWSTREAMTRACES);
	TecUtilArgListAppendArbParam(CurrentArgList, SV_IVALUE, TRUE);
	TecUtilStyleSetLowLevelX(CurrentArgList);

	//===============================================================
	//	First loop to apply specific settings
	//===============================================================

	for (GroupIndex = 0 ; GroupIndex < NumOfSets ; GroupIndex++)
	{

		if (!TecUtilSetIsEmpty(ZoneSets[GroupIndex]))
		{
			//	Prepare temporary arg and object set to make changes to
			TecUtilArgListClear(CurrentArgList);
			TecUtilSetClear(CurrentObjectSet);
			TecUtilSetCopy(CurrentObjectSet, ZoneSets[GroupIndex], TRUE);

			//	Modify group number
			TecUtilArgListAppendString(CurrentArgList, SV_P1, SV_FIELDMAP);
			TecUtilArgListAppendString(CurrentArgList, SV_P2, SV_GROUP);
			TecUtilArgListAppendSet(CurrentArgList, SV_OBJECTSET, CurrentObjectSet);
			TecUtilArgListAppendArbParam(CurrentArgList, SV_IVALUE, (LgIndex_t)(GroupIndex + 1));
			TecUtilStyleSetLowLevelX(CurrentArgList);
			TecUtilArgListClear(CurrentArgList);

			//	Modify Mesh Properties

			//	MeshShow
			TecUtilZoneSetMesh(SV_SHOW, CurrentObjectSet, 0.0, Zone.MeshShow[GroupIndex]);

			//	MeshType
			TecUtilZoneSetMesh(SV_MESHTYPE, CurrentObjectSet, 0.0, Zone.MeshType[GroupIndex]);

			//	MeshColor
			TecUtilZoneSetMesh(SV_COLOR, CurrentObjectSet, 0.0, Zone.MeshColor[GroupIndex]);

			//	MeshLinePttrn
			TecUtilZoneSetMesh(SV_LINEPATTERN, CurrentObjectSet, 0.0, Zone.MeshLinePttrn[GroupIndex]);

			//	MeshLinePttrnLen
			TecUtilZoneSetMesh(SV_PATTERNLENGTH, CurrentObjectSet, Zone.MeshLinePttrnLen[GroupIndex], NULL);

			//	MeshLineThickness
			TecUtilZoneSetMesh(SV_LINETHICKNESS, CurrentObjectSet, Zone.MeshLineThickness[GroupIndex], NULL);

			//	Modify Contour Properties

			//	ContourShow
			TecUtilZoneSetContour(SV_SHOW, CurrentObjectSet, 0.0, Zone.ContourShow[GroupIndex]);

			//	ContourType
			TecUtilZoneSetContour(SV_CONTOURTYPE, CurrentObjectSet, 0.0, Zone.ContourType[GroupIndex]);

			//	FloodColoring
			TecUtilZoneSetContour(SV_FLOODCOLORING, CurrentObjectSet, 0.0, Zone.FloodColoring[GroupIndex]);

			//	ContourLineGroup
			TecUtilZoneSetContour(SV_LINECONTOURGROUP, CurrentObjectSet, 0.0, Zone.ContourLineGroup[GroupIndex]);

			//	ContourLineColor
			TecUtilZoneSetContour(SV_COLOR, CurrentObjectSet, 0.0, Zone.ContourLineColor[GroupIndex]);

			//	ContourLinePttrn
			TecUtilZoneSetContour(SV_LINEPATTERN, CurrentObjectSet, 0.0, Zone.ContourLinePttrn[GroupIndex]);

			//	ContourLinePttrnLen
			TecUtilZoneSetContour(SV_PATTERNLENGTH, CurrentObjectSet, Zone.ContourLinePttrnLen[GroupIndex], NULL);

			//	ContourLineThickness
			TecUtilZoneSetContour(SV_LINETHICKNESS, CurrentObjectSet, Zone.ContourLineThickness[GroupIndex], NULL);

			//	ContourLighting
			TecUtilZoneSetContour(SV_USELIGHTINGEFFECT, CurrentObjectSet, 0.0, Zone.ContourLighting[GroupIndex]);

			//	Modify Vector Properties

			//	VectorShow
			TecUtilZoneSetVector(SV_SHOW, CurrentObjectSet, 0.0, Zone.VectorShow[GroupIndex]);

			//	VectorType
			TecUtilZoneSetVector(SV_VECTORTYPE, CurrentObjectSet, 0.0, Zone.VectorType[GroupIndex]);

			//	VectorArrowType
			TecUtilZoneSetVector(SV_ARROWHEADSTYLE, CurrentObjectSet, 0.0, Zone.VectorArrowType[GroupIndex]);

			//	VectorLineColor
			TecUtilZoneSetVector(SV_COLOR, CurrentObjectSet, 0.0, Zone.VectorLineColor[GroupIndex]);

			//	VectorTangent
			TecUtilZoneSetVector(SV_ISTANGENT, CurrentObjectSet, 0.0, Zone.VectorTangent[GroupIndex]);

			//	VectorLinePttrn
			TecUtilZoneSetVector(SV_LINEPATTERN, CurrentObjectSet, 0.0, Zone.VectorLinePttrn[GroupIndex]);

			//	VectorLinePttrnLen
			TecUtilZoneSetVector(SV_PATTERNLENGTH, CurrentObjectSet, Zone.VectorLinePttrnLen[GroupIndex], NULL);

			//	VectorLineThickness
			TecUtilZoneSetVector(SV_LINETHICKNESS, CurrentObjectSet, Zone.VectorLineThickness[GroupIndex], NULL);

			if (GroupIndex != 27){
				//	Modify Scatter Properties

				//	ScatterShow
				TecUtilZoneSetScatter(SV_SHOW, CurrentObjectSet, 0.0, Zone.ScatterShow[GroupIndex]);

				//	ScatterShape
				TecUtilZoneSetScatterSymbolShape(SV_GEOMSHAPE, CurrentObjectSet, Zone.ScatterShape[GroupIndex]);

				//	ScatterOutColor  Outline color
				TecUtilZoneSetScatter(SV_COLOR, CurrentObjectSet, 0.0, Zone.ScatterOutColor[GroupIndex]);

				//	ScatterFillMode
				TecUtilZoneSetScatter(SV_FILLMODE, CurrentObjectSet, 0.0, Zone.ScatterFillMode[GroupIndex]);

				//	ScatterInColor  fill color
				TecUtilZoneSetScatter(SV_FILLCOLOR, CurrentObjectSet, 0.0, Zone.ScatterInColor[GroupIndex]);

				//	ScatterSize
				TecUtilZoneSetScatter(SV_FRAMESIZE, CurrentObjectSet, Zone.ScatterSize[GroupIndex], NULL);

				//	ScatterLineThickness
				TecUtilZoneSetScatter(SV_LINETHICKNESS, CurrentObjectSet, Zone.ScatterLineThickness[GroupIndex], NULL);
			}

			//	Modify Shade Properties

			//	ShadeShow
			TecUtilZoneSetShade(SV_SHOW, CurrentObjectSet, 0.0, Zone.ShadeShow[GroupIndex]);

			//	ShadeColor
			TecUtilZoneSetShade(SV_COLOR, CurrentObjectSet, 0.0, Zone.ShadeColor[GroupIndex]);

			//	ShadeLighting
			TecUtilZoneSetShade(SV_USELIGHTINGEFFECT, CurrentObjectSet, 0.0, Zone.ShadeLighting[GroupIndex]);

			//	Modify Edge Properties

			//	EdgeShow
			TecUtilZoneSetEdgeLayer(SV_SHOW, CurrentObjectSet, 0.0, Zone.EdgeShow[GroupIndex]);

			//	EdgeType
			TecUtilZoneSetEdgeLayer(SV_EDGETYPE, CurrentObjectSet, 0.0, Zone.EdgeType[GroupIndex]);

			//	EdgeIBorder
			TecUtilZoneSetEdgeLayer(SV_IBORDER, CurrentObjectSet, 0.0, Zone.EdgeIBorder[GroupIndex]);

			//	EdgeJBorder
			TecUtilZoneSetEdgeLayer(SV_JBORDER, CurrentObjectSet, 0.0, Zone.EdgeJBorder[GroupIndex]);

			//	EdgeKBorder
			TecUtilZoneSetEdgeLayer(SV_KBORDER, CurrentObjectSet, 0.0, Zone.EdgeKBorder[GroupIndex]);

			//	EdgeColor
			TecUtilZoneSetEdgeLayer(SV_COLOR, CurrentObjectSet, 0.0, Zone.EdgeColor[GroupIndex]);

			//	EdgeLineThickness
			TecUtilZoneSetEdgeLayer(SV_LINETHICKNESS, CurrentObjectSet, Zone.EdgeLineThickness[GroupIndex], NULL);

			//	Modify Points Properties

			//	PointsToPlot
			TecUtilArgListAppendString(CurrentArgList, SV_P1, SV_FIELDMAP);
			TecUtilArgListAppendString(CurrentArgList, SV_P2, SV_POINTS);
			TecUtilArgListAppendString(CurrentArgList, SV_P3, SV_POINTSTOPLOT);
			TecUtilArgListAppendSet(CurrentArgList, SV_OBJECTSET, CurrentObjectSet);
			TecUtilArgListAppendArbParam(CurrentArgList, SV_IVALUE, Zone.PointsToPlot[GroupIndex]);
			TecUtilStyleSetLowLevelX(CurrentArgList);
			TecUtilArgListClear(CurrentArgList);

			//	PointsIndexSkipI
			TecUtilArgListAppendString(CurrentArgList, SV_P1, SV_FIELDMAP);
			TecUtilArgListAppendString(CurrentArgList, SV_P2, SV_POINTS);
			TecUtilArgListAppendString(CurrentArgList, SV_P3, SV_IJKSKIP);
			TecUtilArgListAppendString(CurrentArgList, SV_P4, SV_I);
			TecUtilArgListAppendSet(CurrentArgList, SV_OBJECTSET, CurrentObjectSet);
			TecUtilArgListAppendArbParam(CurrentArgList, SV_IVALUE, Zone.PointsIndexSkipI[GroupIndex]);
			TecUtilStyleSetLowLevelX(CurrentArgList);
			TecUtilArgListClear(CurrentArgList);

			//	PointsIndexSkipJ
			TecUtilArgListAppendString(CurrentArgList, SV_P1, SV_FIELDMAP);
			TecUtilArgListAppendString(CurrentArgList, SV_P2, SV_POINTS);
			TecUtilArgListAppendString(CurrentArgList, SV_P3, SV_IJKSKIP);
			TecUtilArgListAppendString(CurrentArgList, SV_P4, SV_J);
			TecUtilArgListAppendSet(CurrentArgList, SV_OBJECTSET, CurrentObjectSet);
			TecUtilArgListAppendArbParam(CurrentArgList, SV_IVALUE, Zone.PointsIndexSkipJ[GroupIndex]);
			TecUtilStyleSetLowLevelX(CurrentArgList);
			TecUtilArgListClear(CurrentArgList);

			//	PointsIndexSkipK
			TecUtilArgListAppendString(CurrentArgList, SV_P1, SV_FIELDMAP);
			TecUtilArgListAppendString(CurrentArgList, SV_P2, SV_POINTS);
			TecUtilArgListAppendString(CurrentArgList, SV_P3, SV_IJKSKIP);
			TecUtilArgListAppendString(CurrentArgList, SV_P4, SV_K);
			TecUtilArgListAppendSet(CurrentArgList, SV_OBJECTSET, CurrentObjectSet);
			TecUtilArgListAppendArbParam(CurrentArgList, SV_IVALUE, Zone.PointsIndexSkipK[GroupIndex]);
			TecUtilStyleSetLowLevelX(CurrentArgList);
			TecUtilArgListClear(CurrentArgList);

			//	Modify Surfaces Properties

			//	SurfacesToPlot
			TecUtilArgListAppendString(CurrentArgList, SV_P1, SV_FIELDMAP);
			TecUtilArgListAppendString(CurrentArgList, SV_P2, SV_SURFACES);
			TecUtilArgListAppendString(CurrentArgList, SV_P3, SV_SURFACESTOPLOT);
			TecUtilArgListAppendSet(CurrentArgList, SV_OBJECTSET, CurrentObjectSet);
			TecUtilArgListAppendArbParam(CurrentArgList, SV_IVALUE, Zone.SurfacesToPlot[GroupIndex]);
			TecUtilStyleSetLowLevelX(CurrentArgList);
			TecUtilArgListClear(CurrentArgList);

			//	RangeBeginI
			TecUtilArgListAppendString(CurrentArgList, SV_P1, SV_FIELDMAP);
			TecUtilArgListAppendString(CurrentArgList, SV_P2, SV_SURFACES);
			TecUtilArgListAppendString(CurrentArgList, SV_P3, SV_IRANGE);
			TecUtilArgListAppendString(CurrentArgList, SV_P4, SV_MIN);
			TecUtilArgListAppendSet(CurrentArgList, SV_OBJECTSET, CurrentObjectSet);
			TecUtilArgListAppendArbParam(CurrentArgList, SV_IVALUE, Zone.RangeBeginI[GroupIndex]);
			TecUtilStyleSetLowLevelX(CurrentArgList);
			TecUtilArgListClear(CurrentArgList);
			//	RangeEndI
			TecUtilArgListAppendString(CurrentArgList, SV_P1, SV_FIELDMAP);
			TecUtilArgListAppendString(CurrentArgList, SV_P2, SV_SURFACES);
			TecUtilArgListAppendString(CurrentArgList, SV_P3, SV_IRANGE);
			TecUtilArgListAppendString(CurrentArgList, SV_P4, SV_MAX);
			TecUtilArgListAppendSet(CurrentArgList, SV_OBJECTSET, CurrentObjectSet);
			TecUtilArgListAppendArbParam(CurrentArgList, SV_IVALUE, Zone.RangeEndI[GroupIndex]);
			TecUtilStyleSetLowLevelX(CurrentArgList);
			TecUtilArgListClear(CurrentArgList);
			//	RangeSkipI
			TecUtilArgListAppendString(CurrentArgList, SV_P1, SV_FIELDMAP);
			TecUtilArgListAppendString(CurrentArgList, SV_P2, SV_SURFACES);
			TecUtilArgListAppendString(CurrentArgList, SV_P3, SV_IRANGE);
			TecUtilArgListAppendString(CurrentArgList, SV_P4, SV_SKIP);
			TecUtilArgListAppendSet(CurrentArgList, SV_OBJECTSET, CurrentObjectSet);
			TecUtilArgListAppendArbParam(CurrentArgList, SV_IVALUE, Zone.RangeSkipI[GroupIndex]);
			TecUtilStyleSetLowLevelX(CurrentArgList);
			TecUtilArgListClear(CurrentArgList);

			//	RangeBeginJ
			TecUtilArgListAppendString(CurrentArgList, SV_P1, SV_FIELDMAP);
			TecUtilArgListAppendString(CurrentArgList, SV_P2, SV_SURFACES);
			TecUtilArgListAppendString(CurrentArgList, SV_P3, SV_JRANGE);
			TecUtilArgListAppendString(CurrentArgList, SV_P4, SV_MIN);
			TecUtilArgListAppendSet(CurrentArgList, SV_OBJECTSET, CurrentObjectSet);
			TecUtilArgListAppendArbParam(CurrentArgList, SV_IVALUE, Zone.RangeBeginJ[GroupIndex]);
			TecUtilStyleSetLowLevelX(CurrentArgList);
			TecUtilArgListClear(CurrentArgList);
			//	RangeEndJ
			TecUtilArgListAppendString(CurrentArgList, SV_P1, SV_FIELDMAP);
			TecUtilArgListAppendString(CurrentArgList, SV_P2, SV_SURFACES);
			TecUtilArgListAppendString(CurrentArgList, SV_P3, SV_JRANGE);
			TecUtilArgListAppendString(CurrentArgList, SV_P4, SV_MAX);
			TecUtilArgListAppendSet(CurrentArgList, SV_OBJECTSET, CurrentObjectSet);
			TecUtilArgListAppendArbParam(CurrentArgList, SV_IVALUE, Zone.RangeEndJ[GroupIndex]);
			TecUtilStyleSetLowLevelX(CurrentArgList);
			TecUtilArgListClear(CurrentArgList);
			//	RangeSkipJ
			TecUtilArgListAppendString(CurrentArgList, SV_P1, SV_FIELDMAP);
			TecUtilArgListAppendString(CurrentArgList, SV_P2, SV_SURFACES);
			TecUtilArgListAppendString(CurrentArgList, SV_P3, SV_JRANGE);
			TecUtilArgListAppendString(CurrentArgList, SV_P4, SV_SKIP);
			TecUtilArgListAppendSet(CurrentArgList, SV_OBJECTSET, CurrentObjectSet);
			TecUtilArgListAppendArbParam(CurrentArgList, SV_IVALUE, Zone.RangeSkipJ[GroupIndex]);
			TecUtilStyleSetLowLevelX(CurrentArgList);
			TecUtilArgListClear(CurrentArgList);

			//	RangeBeginK
			TecUtilArgListAppendString(CurrentArgList, SV_P1, SV_FIELDMAP);
			TecUtilArgListAppendString(CurrentArgList, SV_P2, SV_SURFACES);
			TecUtilArgListAppendString(CurrentArgList, SV_P3, SV_KRANGE);
			TecUtilArgListAppendString(CurrentArgList, SV_P4, SV_MIN);
			TecUtilArgListAppendSet(CurrentArgList, SV_OBJECTSET, CurrentObjectSet);
			TecUtilArgListAppendArbParam(CurrentArgList, SV_IVALUE, Zone.RangeBeginK[GroupIndex]);
			TecUtilStyleSetLowLevelX(CurrentArgList);
			TecUtilArgListClear(CurrentArgList);
			//	RangeEndK
			TecUtilArgListAppendString(CurrentArgList, SV_P1, SV_FIELDMAP);
			TecUtilArgListAppendString(CurrentArgList, SV_P2, SV_SURFACES);
			TecUtilArgListAppendString(CurrentArgList, SV_P3, SV_KRANGE);
			TecUtilArgListAppendString(CurrentArgList, SV_P4, SV_MAX);
			TecUtilArgListAppendSet(CurrentArgList, SV_OBJECTSET, CurrentObjectSet);
			TecUtilArgListAppendArbParam(CurrentArgList, SV_IVALUE, Zone.RangeEndK[GroupIndex]);
			TecUtilStyleSetLowLevelX(CurrentArgList);
			TecUtilArgListClear(CurrentArgList);
			//	RangeSkipK
			TecUtilArgListAppendString(CurrentArgList, SV_P1, SV_FIELDMAP);
			TecUtilArgListAppendString(CurrentArgList, SV_P2, SV_SURFACES);
			TecUtilArgListAppendString(CurrentArgList, SV_P3, SV_KRANGE);
			TecUtilArgListAppendString(CurrentArgList, SV_P4, SV_SKIP);
			TecUtilArgListAppendSet(CurrentArgList, SV_OBJECTSET, CurrentObjectSet);
			TecUtilArgListAppendArbParam(CurrentArgList, SV_IVALUE, Zone.RangeSkipK[GroupIndex]);
			TecUtilStyleSetLowLevelX(CurrentArgList);
			TecUtilArgListClear(CurrentArgList);

			//	Modify Volume Properties

			//	StreamTracesShow
			TecUtilArgListAppendString(CurrentArgList, SV_P1, SV_FIELDMAP);
			TecUtilArgListAppendString(CurrentArgList, SV_P2, SV_VOLUMEMODE);
			TecUtilArgListAppendString(CurrentArgList, SV_P3, SV_VOLUMEOBJECTSTOPLOT);
			TecUtilArgListAppendString(CurrentArgList, SV_P4, SV_SHOWSTREAMTRACES);
			TecUtilArgListAppendSet(CurrentArgList, SV_OBJECTSET, CurrentObjectSet);
			TecUtilArgListAppendArbParam(CurrentArgList, SV_IVALUE, Zone.StreamtracesShow[GroupIndex]);
			TecUtilStyleSetLowLevelX(CurrentArgList);
			TecUtilArgListClear(CurrentArgList);

			//	IsoSurfacesShow
			TecUtilArgListAppendString(CurrentArgList, SV_P1, SV_FIELDMAP);
			TecUtilArgListAppendString(CurrentArgList, SV_P2, SV_VOLUMEMODE);
			TecUtilArgListAppendString(CurrentArgList, SV_P3, SV_VOLUMEOBJECTSTOPLOT);
			TecUtilArgListAppendString(CurrentArgList, SV_P4, SV_SHOWISOSURFACES);
			TecUtilArgListAppendSet(CurrentArgList, SV_OBJECTSET, CurrentObjectSet);
			TecUtilArgListAppendArbParam(CurrentArgList, SV_IVALUE, Zone.IsoSurfacesShow[GroupIndex]);
			TecUtilStyleSetLowLevelX(CurrentArgList);
			TecUtilArgListClear(CurrentArgList);

			//	SlicesShow
			TecUtilArgListAppendString(CurrentArgList, SV_P1, SV_FIELDMAP);
			TecUtilArgListAppendString(CurrentArgList, SV_P2, SV_VOLUMEMODE);
			TecUtilArgListAppendString(CurrentArgList, SV_P3, SV_VOLUMEOBJECTSTOPLOT);
			TecUtilArgListAppendString(CurrentArgList, SV_P4, SV_SHOWSLICES);
			TecUtilArgListAppendSet(CurrentArgList, SV_OBJECTSET, CurrentObjectSet);
			TecUtilArgListAppendArbParam(CurrentArgList, SV_IVALUE, Zone.SlicesShow[GroupIndex]);
			TecUtilStyleSetLowLevelX(CurrentArgList);
			TecUtilArgListClear(CurrentArgList);

			//	Modify Effects

			//	UseValueBlanking
			TecUtilArgListAppendString(CurrentArgList, SV_P1, SV_FIELDMAP);
			TecUtilArgListAppendString(CurrentArgList, SV_P2, SV_EFFECTS);
			TecUtilArgListAppendString(CurrentArgList, SV_P3, SV_USEVALUEBLANKING);
			TecUtilArgListAppendSet(CurrentArgList, SV_OBJECTSET, CurrentObjectSet);
			TecUtilArgListAppendArbParam(CurrentArgList, SV_IVALUE, Zone.UseValueBlanking[GroupIndex]);
			TecUtilStyleSetLowLevelX(CurrentArgList);
			TecUtilArgListClear(CurrentArgList);

			////	UseClipPlanes
			//TecUtilArgListAppendString(CurrentArgList, SV_P1, SV_FIELDMAP);
			//TecUtilArgListAppendString(CurrentArgList, SV_P2, SV_EFFECTS);
			//TecUtilArgListAppendString(CurrentArgList, SV_P3, SV_USECLIPPLANES);
			//TecUtilArgListAppendSet(CurrentArgList, SV_OBJECTSET, CurrentObjectSet);
			//TecUtilArgListAppendArbParam(CurrentArgList, SV_IVALUE, Zone.UseClipPlanes[GroupIndex]);
			//TecUtilStyleSetLowLevelX(CurrentArgList);
			//TecUtilArgListClear(CurrentArgList);

			//	UseTranslucency
			TecUtilArgListAppendString(CurrentArgList, SV_P1, SV_FIELDMAP);
			TecUtilArgListAppendString(CurrentArgList, SV_P2, SV_EFFECTS);
			TecUtilArgListAppendString(CurrentArgList, SV_P3, SV_USETRANSLUCENCY);
			TecUtilArgListAppendSet(CurrentArgList, SV_OBJECTSET, CurrentObjectSet);
			TecUtilArgListAppendArbParam(CurrentArgList, SV_IVALUE, Zone.UseSurfaceTranslucency[GroupIndex]);
			TecUtilStyleSetLowLevelX(CurrentArgList);
			TecUtilArgListClear(CurrentArgList);

			//	SurfaceTranslucency
			TecUtilArgListAppendString(CurrentArgList, SV_P1, SV_FIELDMAP);
			TecUtilArgListAppendString(CurrentArgList, SV_P2, SV_EFFECTS);
			TecUtilArgListAppendString(CurrentArgList, SV_P3, SV_SURFACETRANSLUCENCY);
			TecUtilArgListAppendSet(CurrentArgList, SV_OBJECTSET, CurrentObjectSet);
			TecUtilArgListAppendArbParam(CurrentArgList, SV_IVALUE, Zone.SurfaceTranslucency[GroupIndex]);
			TecUtilStyleSetLowLevelX(CurrentArgList);
			TecUtilArgListClear(CurrentArgList);

			//	LightingEffect
			TecUtilArgListAppendString(CurrentArgList, SV_P1, SV_FIELDMAP);
			TecUtilArgListAppendString(CurrentArgList, SV_P2, SV_EFFECTS);
			TecUtilArgListAppendString(CurrentArgList, SV_P3, SV_LIGHTINGEFFECT);
			TecUtilArgListAppendSet(CurrentArgList, SV_OBJECTSET, CurrentObjectSet);
			TecUtilArgListAppendArbParam(CurrentArgList, SV_IVALUE, Zone.LightingEffect[GroupIndex]);
			TecUtilStyleSetLowLevelX(CurrentArgList);
			TecUtilArgListClear(CurrentArgList);
		}
	}

	//===============================================================
	//	Second loop to turn zones on/off
	//===============================================================

	//	Prepare temporary arg and object set to make changes to
	TecUtilArgListClear(CurrentArgList);
	TecUtilSetClear(CurrentObjectSet);
	SetIndex_t		TempMember;

	//	Set argument list for active field maps
	TecUtilArgListAppendString(CurrentArgList, SV_P1, SV_ACTIVEFIELDMAPS);

	for (GroupIndex = 0 ; GroupIndex < NumOfSets ; GroupIndex++)
	{
		if (Zone.ZoneShow[GroupIndex] && !TecUtilSetIsEmpty(ZoneSets[GroupIndex]))
		{
			TecUtilSetForEachMember(TempMember,ZoneSets[GroupIndex])
			{
				TecUtilSetAddMember(CurrentObjectSet, TempMember, TRUE);
			}
		}
	}

	TecUtilArgListAppendArbParam(CurrentArgList, SV_IVALUE, (ArbParam_t)CurrentObjectSet);
	TecUtilStyleSetLowLevelX(CurrentArgList);


	//===============================================================
	//	Activate 3d perspective
	//===============================================================

	LgIndex_t	iMax, jMax, kMax, ijkMax;
	EntIndex_t	NumOfZones = TecUtilDataSetGetNumZones();
	//char what[200];
	TecUtilZoneGetIJK(1, &iMax, &jMax, &kMax);
	ijkMax = MAX(iMax, MAX(jMax, kMax));

	//int NumCrtPts = 0;
	//char *TempNum = NULL;
	//AuxData_pa	TempAuxDataRef;
	//char *ZoneName = NULL, *JunkName = NULL;
	////AuxDataType_e JunkDataType = AuxDataType_String;
	//Boolean_t JunkRetain = NULL;

	//char aName[] = "CompChem.NumCrtPtAtom";
	//char bName[] = "CompChem.NumCrtPtBond";
	//char cName[] = "CompChem.NumCrtPtCage";
	//char rName[] = "CompChem.NumCrtPtRing";

	//char *Names[] = {aName, bName, cName, rName};

	//if (NumOfZones > 1)	TempAuxDataRef = TecUtilAuxDataZoneGetRef((EntIndex_t)2);

	//if (TempAuxDataRef != NULL)
	//{
	//	if (TecUtilZoneGetName((EntIndex_t)2, &ZoneName))
	//	{
	//		if (strncmp(ZoneName, "Critical Points Zone ", 21) == 0)
	//		{
	//			for (int i = 0 ; i <= 3 ; i++)
	//			{
	//				TecUtilAuxDataGetStrItemByName(TempAuxDataRef, Names[i], &TempNum, &JunkRetain);
	//				NumCrtPts += (int)TempNum;
	//				sprintf_s(what,200, "tempcrt %d", (int)TempNum);
	//				TecUtilDialogMessageBox(what, MessageBox_Information);
	//			}
	//			TecUtilStringDealloc(&JunkName);
	//			//TecUtilStringDealloc(&TempNum);
	//		}
	//		TecUtilStringDealloc(&ZoneName);
	//	}
	//}
	////sprintf_s(what,200, "ijkmax %d", ijkMax);
	////TecUtilDialogMessageBox(what, MessageBox_Information);
	////sprintf_s(what,200, "crtpoints %d", NumCrtPts);
	////TecUtilDialogMessageBox(what, MessageBox_Information);

	//double	EyeDistance = (double)ijkMax;
	//if (NumCrtPts > 0)
	//{
	//	if ((iMax * jMax * kMax) / NumCrtPts > 5000)
	//	{
	//		EyeDistance = (double)ijkMax * 0.4;
	//		TecUtilDialogMessageBox("yo",MessageBox_Information);
	//	}
	//}

	//TecUtilViewZoom((double)-ijkMax/2, (double)-ijkMax/2, (double)ijkMax/2, (double)ijkMax/2);
	//TecUtilViewSetMagnification((double)1.2);
	
// 	TecUtilArgListClear(CurrentArgList);
// 	TecUtilArgListAppendString(CurrentArgList, SV_P1, SV_THREEDVIEW);
// 	TecUtilArgListAppendString(CurrentArgList, SV_P2, SV_VIEWWIDTH);
// 	TecUtilArgListAppendDouble(CurrentArgList, SV_DVALUE, (double)ijkMax);
// 	TecUtilStyleSetLowLevelX(CurrentArgList);
	
	TecUtilArgListClear(CurrentArgList);
	TecUtilArgListAppendString(CurrentArgList, SV_P1, SV_THREEDVIEW);
	TecUtilArgListAppendString(CurrentArgList, SV_P2, SV_DRAWINPERSPECTIVE);
	TecUtilArgListAppendArbParam(CurrentArgList, SV_IVALUE, (Boolean_t)1);
	TecUtilStyleSetLowLevelX(CurrentArgList);
	
	TecUtilArgListClear(CurrentArgList);
	TecUtilArgListAppendString(CurrentArgList, SV_P1, SV_THREEDVIEW);
	TecUtilArgListAppendString(CurrentArgList, SV_P2, SV_FIELDOFVIEW);
	TecUtilArgListAppendDouble(CurrentArgList, SV_DVALUE, (double)20.0);
	TecUtilStyleSetLowLevelX(CurrentArgList);

	TecUtilSet3DEyeDistance((double)ijkMax*0.4);

	TecUtilViewNiceFit();



	//	deallocate temporary set and arg list
	TecUtilSetDealloc(&CurrentObjectSet);
	TecUtilArgListDealloc(&CurrentArgList);

}	//	void ApplyZoneSettings(std::vector<Set_pa> ZoneSets, ZoneSettings Zone)


void GroupZones(std::vector <Set_pa> &ZoneSets)
{
	//	First create object sets for each zone type.
	//	After each set is created, set it's style.

	//	Get number of zones in the data set
	EntIndex_t NumberOfZones = TecUtilDataSetGetNumZones();

	//	Strings used to identify zones. Replace as required by chembond output.
	//	CompChem.ZoneType for volume and surface CP zones
	char CPZoneType[] =				"CriticalPoints";
	char CPZoneName[] = 			"Critical Points Zone ";

	//	CompChem.BegCrtPtType, CompChem.EndCrtPtType
	//	CompChem.ThrdCrtPtType, CompChem.OppositeCrtPtType
	//	For all lines and surfaces except non-FF bonds
	char AtomCP[] =					"Atom";
	char BondCP[] =					"Bond";
	char RingCP[] =					"Ring";
	char RingCPFF[] =				"RingFF";
	char CageCP[] =					"Cage";
	char CageCPFF[] =				"CageFF";
	char FFCP[] =					"FarField";

	//	For zones lacking aux data so their names are used for
	//	partial string comparisons (luckily they're few)
	char ExtraIsoSurfacesStr[] =	"Iso: C=";
	char AtomRhoIsoStr[] =			"SurfTopoSegZone";
	char BondZoneStr[] =			"Bond ";

	//	Can't predict the name of the data sets that exist before
	//	chembond is run, so just put those in their own group
	//	assuming they dont share the same first characters
	//	as the above 3 strings.
	
	//	allocate sets for each zone type (FF = far field)
	ZoneSets.resize(28);

	Set_pa InitialDataSetZone =	ZoneSets[0] = TecUtilSetAlloc(TRUE);
	Set_pa VolumeCPs =			ZoneSets[1] = TecUtilSetAlloc(TRUE);
	Set_pa IsoSurfCPs =			ZoneSets[2] = TecUtilSetAlloc(TRUE);
	Set_pa ABLines =			ZoneSets[3] = TecUtilSetAlloc(TRUE);
	Set_pa ABLinesFF =			ZoneSets[4] = TecUtilSetAlloc(TRUE);
	Set_pa ARLines =			ZoneSets[5] = TecUtilSetAlloc(TRUE);
	Set_pa ARLinesFF =			ZoneSets[6] = TecUtilSetAlloc(TRUE);
	Set_pa ACLines =			ZoneSets[7] = TecUtilSetAlloc(TRUE);
	Set_pa ACLinesFF =			ZoneSets[8] = TecUtilSetAlloc(TRUE);
	Set_pa BRLines =			ZoneSets[9] = TecUtilSetAlloc(TRUE);
	Set_pa BRLinesFF =			ZoneSets[10] = TecUtilSetAlloc(TRUE);
	Set_pa BCLines =			ZoneSets[11] = TecUtilSetAlloc(TRUE);
	Set_pa BCLinesFF =			ZoneSets[12] = TecUtilSetAlloc(TRUE);
	Set_pa RCLines =			ZoneSets[13] = TecUtilSetAlloc(TRUE);
	Set_pa RCLinesFF =			ZoneSets[14] = TecUtilSetAlloc(TRUE);
	Set_pa ABRSurfaces =		ZoneSets[15] = TecUtilSetAlloc(TRUE);
	Set_pa ABRSurfacesFF =		ZoneSets[16] = TecUtilSetAlloc(TRUE);
	Set_pa ABCSurfaces =		ZoneSets[17] = TecUtilSetAlloc(TRUE);
	Set_pa ABCSurfacesFF =		ZoneSets[18] = TecUtilSetAlloc(TRUE);
	Set_pa ARCSurfaces =		ZoneSets[19] = TecUtilSetAlloc(TRUE);
	Set_pa ARCSurfacesFF =		ZoneSets[20] = TecUtilSetAlloc(TRUE);
	Set_pa BRCSurfaces =		ZoneSets[21] = TecUtilSetAlloc(TRUE);
	Set_pa BRCSurfacesFF =		ZoneSets[22] = TecUtilSetAlloc(TRUE);
	Set_pa AtomRhoIsoSurfaces = ZoneSets[23] = TecUtilSetAlloc(TRUE);
	Set_pa ExtraIsoSurfaces =	ZoneSets[24] = TecUtilSetAlloc(TRUE);
	Set_pa BondZones =			ZoneSets[25] = TecUtilSetAlloc(TRUE);
	Set_pa BondZonesFF =		ZoneSets[26] = TecUtilSetAlloc(TRUE);
	Set_pa BadZones =			ZoneSets[27] = TecUtilSetAlloc(TRUE);

	//	The TecUtilAuxDataGetStrItemByName function takes a string by ref, so
	//	if using complex if statements with that function you need multiple
	//	strings to pass each instance of the function.
	char *ZoneType = NULL;
	char *BegCPType = NULL;
	char *EndCPType = NULL;
	char *ThirdCPType = NULL;
	char *OppCPType = NULL;
	char *FarField = NULL;	//ArbParam_t FarField;
	AuxDataType_e *FarFieldType = NULL;
	char *ZoneName = NULL;

	Boolean_t Retain = FALSE;



	//***************************************************************
	//===============================================================
	//	Group zones based on type into above sets
	//===============================================================
	//***************************************************************

	EntIndex_t CurrentZone = 0;

	

	//	Loop over all zones in the data set
	for (CurrentZone = 1 ; CurrentZone <= NumberOfZones ; CurrentZone++)
	{
		//	Getting reference to current zone so properties can be accessed
		AuxData_pa AuxDataRef = TecUtilAuxDataZoneGetRef(CurrentZone);

		//	Only try to access data if zone reference was retrieved
		if (AuxDataRef != NULL)
		{
			ZoneType = NULL;
			BegCPType = NULL;
			EndCPType = NULL;
			ThirdCPType = NULL;
			OppCPType = NULL;
			FarField = NULL;
			

			//	All lines (since they'll have a beginning CP but no 3rd CP)
			if (TecUtilAuxDataGetStrItemByName(AuxDataRef, "CompChem.BegCrtPtType", &BegCPType, &Retain)
				&& !(TecUtilAuxDataGetStrItemByName(AuxDataRef, "CompChem.ThrdCrtPtType", &ThirdCPType, &Retain))
				&& TecUtilAuxDataGetNumItems(AuxDataRef) < 6)
			{
				//	Need to get the end CP anyways, so use that to be sure it exists
				if (TecUtilAuxDataGetStrItemByName(AuxDataRef, "CompChem.EndCrtPtType", &EndCPType, &Retain))
				{
					//	BA, BC, BR lines
					if (!strcmp(BegCPType,BondCP))
					{
						if (!strcmp(EndCPType, AtomCP)) TecUtilSetAddMember(ABLines, CurrentZone, TRUE);
						if (!strcmp(EndCPType, CageCP)) TecUtilSetAddMember(BCLines, CurrentZone, TRUE);
						if (!strcmp(EndCPType, CageCPFF)) TecUtilSetAddMember(BCLinesFF, CurrentZone, TRUE);
						if (!strcmp(EndCPType, RingCP)) TecUtilSetAddMember(BRLines, CurrentZone, TRUE);
						if (!strcmp(EndCPType, RingCPFF)) TecUtilSetAddMember(BRLinesFF, CurrentZone, TRUE);
						if (!strcmp(EndCPType, BondCP)) TecUtilSetAddMember(BadZones, CurrentZone, TRUE);
					}
					//	AC, ACFF, ARFF lines
					else if (!strcmp(BegCPType, AtomCP))
					{
						if (!strcmp(EndCPType, CageCP)) TecUtilSetAddMember(ACLines, CurrentZone, TRUE);
						if (!strcmp(EndCPType, CageCPFF)) TecUtilSetAddMember(ACLinesFF, CurrentZone, TRUE);
						if (!strcmp(EndCPType, RingCPFF)) TecUtilSetAddMember(ARLinesFF, CurrentZone, TRUE);
						if (!strcmp(EndCPType, AtomCP)) TecUtilSetAddMember(BadZones, CurrentZone, TRUE);
					}
					//	AR, RC and RC far field lines
					else if (!strcmp(BegCPType, RingCP))
					{
						if (!strcmp(EndCPType, CageCP)) TecUtilSetAddMember(RCLines, CurrentZone, TRUE);
						if (!strcmp(EndCPType, AtomCP)) TecUtilSetAddMember(ARLines, CurrentZone, TRUE);
						//	TODO: fix chembond code to name CageFF correctly
						if (!strcmp(EndCPType, FFCP) || !strcmp(EndCPType, CageCPFF)) 
							TecUtilSetAddMember(RCLinesFF, CurrentZone, TRUE);
						if (!strcmp(EndCPType, RingCP)) TecUtilSetAddMember(BadZones, CurrentZone, TRUE);
						if (!strcmp(EndCPType, RingCPFF)) TecUtilSetAddMember(BadZones, CurrentZone, TRUE);
					}
					TecUtilStringDealloc(&BegCPType);
					TecUtilStringDealloc(&EndCPType);
				}
			}
			//	All surfaces not including isosurfaces
			else if (TecUtilAuxDataGetStrItemByName(AuxDataRef, "CompChem.BegCrtPtType", &BegCPType, &Retain)
				&& TecUtilAuxDataGetStrItemByName(AuxDataRef, "CompChem.EndCrtPtType", &EndCPType, &Retain)
				&& TecUtilAuxDataGetStrItemByName(AuxDataRef, "CompChem.ThrdCrtPtType", &ThirdCPType, &Retain))
			{
				//	Far field bond zones
				if (!strcmp(BegCPType, BondCP) 
					&& ((!strcmp(EndCPType, RingCP) && !strcmp(ThirdCPType, RingCP))
					|| (!strcmp(EndCPType, RingCP) && !strcmp(ThirdCPType, FFCP))
					|| (!strcmp(EndCPType, FFCP) && !strcmp(ThirdCPType, RingCP))
					|| (!strcmp(EndCPType, FFCP) && !strcmp(ThirdCPType, FFCP))))
				{
					TecUtilSetAddMember(BondZonesFF, CurrentZone, TRUE);
				}
				//	Well-behaved ABR surfaces
				else if ((!strcmp(BegCPType, RingCP) || !strcmp(BegCPType, RingCPFF))
					&& (!strcmp(EndCPType, AtomCP) || !strcmp(EndCPType, BondCP))
					&& (!strcmp(ThirdCPType, AtomCP) || !strcmp(ThirdCPType, BondCP))
					&& !(!strcmp(BegCPType, EndCPType) && !strcmp(BegCPType, ThirdCPType)
					&& !strcmp(EndCPType, ThirdCPType)))
				{
					//TecUtilAuxDataGetItemByName(AuxDataRef, "CompChem.IsFarField", &FarField, &FarFieldType, &Retain);
					TecUtilAuxDataGetStrItemByName(AuxDataRef, "CompChem.IsFarField", &FarField, &Retain);
					if (!!strcmp(FarField,"TRUE")) TecUtilSetAddMember(ABRSurfaces, CurrentZone, TRUE);
					else TecUtilSetAddMember(ABRSurfacesFF, CurrentZone, TRUE);
					TecUtilStringDealloc(&FarField);
				}
				//	Well-behaved ABC surfaces
				else if ((!strcmp(BegCPType, AtomCP) || !strcmp(BegCPType, BondCP) || !strcmp(BegCPType, CageCP)
						|| !strcmp(BegCPType, CageCPFF))
					&& (!strcmp(EndCPType, AtomCP) || !strcmp(EndCPType, BondCP) || !strcmp(EndCPType, CageCP)
						|| !strcmp(EndCPType, CageCPFF))
					&& (!strcmp(ThirdCPType, AtomCP) || !strcmp(ThirdCPType, BondCP) || !strcmp(ThirdCPType, CageCP)
						|| !strcmp(ThirdCPType, CageCPFF))
					&& (strcmp(BegCPType, EndCPType) && strcmp(BegCPType, ThirdCPType)
						&& strcmp(EndCPType, ThirdCPType))
					&& !((!strcmp(EndCPType,CageCP) && !strcmp(ThirdCPType,CageCPFF))
						&& (!strcmp(EndCPType,CageCPFF) && !strcmp(ThirdCPType,CageCP))
						&& (!strcmp(EndCPType,CageCP) && !strcmp(BegCPType,CageCPFF))
						&& (!strcmp(EndCPType,CageCPFF) && !strcmp(BegCPType,CageCP))
						&& (!strcmp(BegCPType,CageCP) && !strcmp(ThirdCPType,CageCPFF))
						&& (!strcmp(BegCPType,CageCPFF) && !strcmp(ThirdCPType,CageCP))))
				{
					//TecUtilAuxDataGetItemByName(AuxDataRef, "CompChem.IsFarField", &FarField, &FarFieldType, &Retain);
					TecUtilAuxDataGetStrItemByName(AuxDataRef, "CompChem.IsFarField", &FarField, &Retain);
					if (!!strcmp(FarField,"TRUE")) TecUtilSetAddMember(ABCSurfaces, CurrentZone, TRUE);
					else TecUtilSetAddMember(ABCSurfacesFF, CurrentZone, TRUE);
					TecUtilStringDealloc(&FarField);
				}
				//	Well-behaved ARC surfaces
				else if ((!strcmp(BegCPType, AtomCP))
					&& (!strcmp(EndCPType, RingCP) || !strcmp(EndCPType, CageCP)
					|| !strcmp(EndCPType, RingCPFF) || !strcmp(EndCPType, CageCPFF))
					&& (!strcmp(ThirdCPType, RingCP) || !strcmp(ThirdCPType, CageCP)
					|| !strcmp(ThirdCPType, RingCPFF) || !strcmp(ThirdCPType, CageCPFF))
					&& (strcmp(BegCPType, EndCPType) && strcmp(BegCPType, ThirdCPType)
					&& strcmp(EndCPType, ThirdCPType)))
				{
					//TecUtilAuxDataGetItemByName(AuxDataRef, "CompChem.IsFarField", &FarField, &FarFieldType, &Retain);
					TecUtilAuxDataGetStrItemByName(AuxDataRef, "CompChem.IsFarField", &FarField, &Retain);
					if (!!strcmp(FarField,"TRUE")) TecUtilSetAddMember(ARCSurfaces, CurrentZone, TRUE);
					else TecUtilSetAddMember(ARCSurfacesFF, CurrentZone, TRUE);
					TecUtilStringDealloc(&FarField);
				}
				//	ACC surfaces
				else if (!strcmp(BegCPType, AtomCP) 
					&& (strcmp(EndCPType, CageCP) || strcmp(EndCPType, CageCPFF))
					&& (strcmp(ThirdCPType, CageCP) || strcmp(ThirdCPType, CageCPFF)))
						TecUtilSetAddMember(BadZones, CurrentZone, TRUE);
				//	Well-behaved BRC surfaces
				else if (!strcmp(BegCPType, BondCP)
					&& (!strcmp(EndCPType, RingCP) || !strcmp(EndCPType, RingCPFF) || !strcmp(EndCPType, CageCP)
					|| !strcmp(EndCPType, CageCPFF))
					&& (!strcmp(ThirdCPType, RingCP) || !strcmp(ThirdCPType, RingCPFF) || !strcmp(ThirdCPType, CageCP)
					|| !strcmp(ThirdCPType, CageCPFF))
					&& (strcmp(BegCPType, EndCPType) && strcmp(BegCPType, ThirdCPType)
					&& strcmp(EndCPType, ThirdCPType)))
				{
					//TecUtilAuxDataGetItemByName(AuxDataRef, "CompChem.IsFarField", &FarField, &FarFieldType, &Retain);
					TecUtilAuxDataGetStrItemByName(AuxDataRef, "CompChem.IsFarField", &FarField, &Retain);
					if (!!strcmp(FarField,"TRUE")) TecUtilSetAddMember(BRCSurfaces, CurrentZone, TRUE);
					else TecUtilSetAddMember(BRCSurfacesFF, CurrentZone, TRUE);
					TecUtilStringDealloc(&FarField);
				}

				TecUtilStringDealloc(&BegCPType);
				TecUtilStringDealloc(&EndCPType);
				TecUtilStringDealloc(&ThirdCPType);
			}
			//	If zone doesn't contain the aux data used above, then it's either a
			//	bond (nonFF) zone or isosurface.
			else
			{
				if (TecUtilZoneGetName(CurrentZone, &ZoneName))
				{
					if (strncmp(ZoneName, ExtraIsoSurfacesStr, 7) == 0)
					{
						TecUtilSetAddMember(ExtraIsoSurfaces, CurrentZone, TRUE);
					}
					else if (strncmp(ZoneName, AtomRhoIsoStr, 15) == 0)
					{
						TecUtilSetAddMember(AtomRhoIsoSurfaces, CurrentZone, TRUE);
					}
					else if (strncmp(ZoneName, BondZoneStr, 5) == 0)
					{
						TecUtilSetAddMember(BondZones, CurrentZone, TRUE);
					}
					else if (strncmp(ZoneName, CPZoneName, 21) == 0)
					{
						if (TecUtilSetIsEmpty(VolumeCPs)) TecUtilSetAddMember(VolumeCPs, CurrentZone, TRUE);
						else TecUtilSetAddMember(IsoSurfCPs, CurrentZone, TRUE);
					}
					//	If no other zone name comparisons came back true, then 
					//	zone should be one of the initial zones before chembond ran.
					else
					{
						if (TecUtilSetIsEmpty(InitialDataSetZone)) TecUtilSetAddMember(InitialDataSetZone, CurrentZone, TRUE);
						else TecUtilSetAddMember(BadZones, CurrentZone, TRUE);
					}
					TecUtilStringDealloc(&ZoneName);
				}
				else
				{
					char ErrorMsg[200];
					sprintf_s(ErrorMsg, "Properties not available for zone %d", CurrentZone);
					TecUtilDialogErrMsg(ErrorMsg);
				}
			}
		}
		//	Each zone should be reference-able, so print error if not
		else 
		{
			char ErrorMsg[200];
			sprintf_s(ErrorMsg, "Reference not available for zone %d", CurrentZone);
			TecUtilDialogErrMsg(ErrorMsg);
		}	//	if zone reference is set correctly

	}	//	for loop over all zones
}	//	void GroupZones(std::vector <Set_pa> &ZoneSets)


void SetupMultiColorContourLevels(std::vector<Set_pa> const &ZoneSets)
{
	//***************************************************************
	//===============================================================
	//	multicolor settings
	//===============================================================
	//***************************************************************


	//	Get mins and maxes of rho, total gradient to be used in
	//	multicolor schemes.

	char TotGradStr[] = "Density Gradient Magnitude";
	char TotGradStr2[] = "Mag Density Gradient";
	char RhoStr[] = "rho";
	char CPStr[] = "CritPointType";
	char *TempStr =	NULL;
	EntIndex_t RhoVar = 0;
	EntIndex_t GradTotVar = 0;
	EntIndex_t CPVar = 0;
	Boolean_t  NoCPs = FALSE;
	EntIndex_t NumOfVars = TecUtilDataSetGetNumVars();

	for (EntIndex_t VarNum = 1 ; VarNum <= NumOfVars ; VarNum++)
	{
		if (TecUtilVarGetName(VarNum, &TempStr))
		{
			if (!strcmp(TempStr, TotGradStr))
			{
				GradTotVar = VarNum;
			}
			else if (!strcmp(TempStr, TotGradStr2))
			{
				GradTotVar = VarNum;
			}
			else if (!strcmp(TempStr, RhoStr))
			{
				RhoVar = VarNum;
			}
			else if (!strcmp(TempStr, CPStr))
			{
				CPVar = VarNum;
			}
		}
	}
	if (GradTotVar == 0) GradTotVar = 8;
	if (RhoVar == 0) RhoVar = 4;
	if (CPVar == 0) NoCPs = TRUE; //CPVar = 9;

	double ABLineRhoMin = 1e10, ABLineRhoMax= 0,
		ARLineRhoMin = 1e10, ARLineRhoMax = 0,
		ACLineRhoMin = 1e10, ACLineRhoMax = 0,
		BRLineRhoMin = 1e10, BRLineRhoMax = 0,
		BCLineRhoMin = 1e10, BCLineRhoMax = 0,
		AtomIsoGradTotMin = 1e10, AtomIsoGradTotMax = 0;
	double TempRhoMin = 0, TempRhoMax = 0, TempGradTotMin = 0, TempGradTotMax = 0;

	Set_pa	AtomRhoIsoSurfaces = ZoneSets[23];
	Set_pa	ABLines = ZoneSets[3];
	Set_pa	ARLines = ZoneSets[5];
	Set_pa	ACLines = ZoneSets[7];
	Set_pa	BRLines = ZoneSets[9];
	Set_pa	BCLines = ZoneSets[11];
	Set_pa	RCLines = ZoneSets[13];
	


	SetIndex_t	iZone = 0;

	if (!TecUtilSetIsEmpty(AtomRhoIsoSurfaces))
	{
		TecUtilSetForEachMember(iZone, AtomRhoIsoSurfaces)
		{
			if (TecUtilDataValueGetMinMaxByZoneVar((EntIndex_t)iZone, GradTotVar, &TempGradTotMin, &TempGradTotMax))
			{
				if (TempGradTotMin > 0 && TempGradTotMin < AtomIsoGradTotMin) AtomIsoGradTotMin = TempGradTotMin;
				if (TempGradTotMax > AtomIsoGradTotMax) AtomIsoGradTotMax = TempGradTotMax;
			}
		}
	}

	if (!TecUtilSetIsEmpty(ABLines))
	{
		TecUtilSetForEachMember(iZone, ABLines)
		{
			if (TecUtilDataValueGetMinMaxByZoneVar((EntIndex_t)iZone, RhoVar, &TempRhoMin, &TempRhoMax))
			{
				if (TempRhoMin > 0 && TempRhoMin < ABLineRhoMin) ABLineRhoMin = TempRhoMin;
				if (TempRhoMax > ABLineRhoMax) ABLineRhoMax = TempRhoMax;
			}
		}
	}

	if (!TecUtilSetIsEmpty(ARLines))
	{
		TecUtilSetForEachMember(iZone, ARLines)
		{
			if (TecUtilDataValueGetMinMaxByZoneVar((EntIndex_t)iZone, RhoVar, &TempRhoMin, &TempRhoMax))
			{
				if (TempRhoMin > 0 && TempRhoMin < ARLineRhoMin) ARLineRhoMin = TempRhoMin;
				if (TempRhoMax > ARLineRhoMax) ARLineRhoMax = TempRhoMax;
			}
		}
	}

	if (!TecUtilSetIsEmpty(ACLines))
	{
		TecUtilSetForEachMember(iZone, ACLines)
		{
			if (TecUtilDataValueGetMinMaxByZoneVar((EntIndex_t)iZone, RhoVar, &TempRhoMin, &TempRhoMax))
			{
				if (TempRhoMin > 0 && TempRhoMin < ACLineRhoMin) ACLineRhoMin = TempRhoMin;
				if (TempRhoMax > ACLineRhoMax) ACLineRhoMax = TempRhoMax;
			}
		}
	}

	if (!TecUtilSetIsEmpty(BRLines))
	{
		TecUtilSetForEachMember(iZone, BRLines)
		{
			if (TecUtilDataValueGetMinMaxByZoneVar((EntIndex_t)iZone, RhoVar, &TempRhoMin, &TempRhoMax))
			{
				if (TempRhoMin > 0 && TempRhoMin < BRLineRhoMin) BRLineRhoMin = TempRhoMin;
				if (TempRhoMax > BRLineRhoMax) BRLineRhoMax = TempRhoMax;
			}
		}
	}

	if (!TecUtilSetIsEmpty(BCLines))
	{
		TecUtilSetForEachMember(iZone, BCLines)
		{
			if (TecUtilDataValueGetMinMaxByZoneVar((EntIndex_t)iZone, RhoVar, &TempRhoMin, &TempRhoMax))
			{
				if (TempRhoMin > 0 && TempRhoMin < BCLineRhoMin) BCLineRhoMin = TempRhoMin;
				if (TempRhoMax > BCLineRhoMax) BCLineRhoMax = TempRhoMax;
			}
		}
	}
	else
	{
		if (!TecUtilSetIsEmpty(RCLines))
		{
			TecUtilSetForEachMember(iZone, RCLines)
			{
				if (TecUtilDataValueGetMinMaxByZoneVar((EntIndex_t)iZone, RhoVar, &TempRhoMin, &TempRhoMax))
				{
					if (TempRhoMin > 0 && TempRhoMin < BCLineRhoMin) BCLineRhoMin = TempRhoMin;
					if (TempRhoMax > BCLineRhoMax) BCLineRhoMax = TempRhoMax;
				}
			}
		}
	}


	//	Set multicolor (contour color) for rho, gradient total, 
	//	and critical point types
	const int NumOfLevels = 20;
	std::vector<double> ABLineRhoV;
	std::vector<double> ARLineRhoV;
	std::vector<double> ACLineRhoV;
	std::vector<double> BRLineRhoV;
	std::vector<double> BCLineRhoV;
	double AtomIsoSurfaceGradTot[NumOfLevels];
	double VolumeCPRanks[] = {-3, -1, 1, 3};
	double SurfaceCPRanks[] = {-2, 0, 2};

	LogDataValues(ABLineRhoV, ABLineRhoMin, ABLineRhoMax, 1e-1);
	LogDataValues(ARLineRhoV, ARLineRhoMin, ARLineRhoMax, 1e-2);
	LogDataValues(ACLineRhoV, ACLineRhoMin, ACLineRhoMax, 1e-2);
	LogDataValues(BRLineRhoV, BRLineRhoMin, BRLineRhoMax, 1e-2);
	LogDataValues(BCLineRhoV, BCLineRhoMin, BCLineRhoMax, 1e-3);

	double *ABLineRho = &ABLineRhoV[0];
	double *ARLineRho = &ARLineRhoV[0];
	double *ACLineRho = &ACLineRhoV[0];
	double *BRLineRho = &BRLineRhoV[0];
	double *BCLineRho = &BCLineRhoV[0];

	for (int i = 0 ; i < NumOfLevels ; i++)
	{
		if ((AtomIsoGradTotMax - AtomIsoGradTotMin) > 0)
		{
			AtomIsoSurfaceGradTot[i] = AtomIsoGradTotMin + 
				(AtomIsoGradTotMax - AtomIsoGradTotMin) / (NumOfLevels - 1) * i;
		}
		else AtomIsoSurfaceGradTot[i] = 0.0;
	}
	
	double *DataArray[] = {VolumeCPRanks, 
							SurfaceCPRanks, 
							AtomIsoSurfaceGradTot, 
							ABLineRho,
							ARLineRho,
							ACLineRho,
							BRLineRho,
							BCLineRho};

	int DataCounts[] = {4, 3, 
							countOf(AtomIsoSurfaceGradTot),
							(int)ABLineRhoV.size(),
							(int)ARLineRhoV.size(),
							(int)ACLineRhoV.size(),
							(int)BRLineRhoV.size(),
							(int)BCLineRhoV.size()};

	ArgList_pa TempArgList = TecUtilArgListAlloc();

	for (int i = 0; i < countOf(DataCounts) ; i++)
	{
		if (i < 2 && !NoCPs)
		{
			TecUtilArgListAppendInt(TempArgList, SV_CONTOURGROUP, (LgIndex_t)(i+1));
			TecUtilArgListAppendInt(TempArgList, SV_VAR, CPVar);
			TecUtilContourSetVariableX(TempArgList);
			TecUtilArgListClear(TempArgList);
		}
		else if (i == 2)
		{
			TecUtilArgListAppendInt(TempArgList, SV_CONTOURGROUP, (LgIndex_t)(i+1));
			TecUtilArgListAppendInt(TempArgList, SV_VAR, GradTotVar);
			TecUtilContourSetVariableX(TempArgList);
			TecUtilArgListClear(TempArgList);
		}
		else
		{
			TecUtilArgListAppendInt(TempArgList, SV_CONTOURGROUP, (LgIndex_t)(i+1));
			TecUtilArgListAppendInt(TempArgList, SV_VAR, RhoVar);
			TecUtilContourSetVariableX(TempArgList);
			TecUtilArgListClear(TempArgList);
		}
		
		TecUtilArgListAppendInt(TempArgList, SV_CONTOURLEVELACTION, ContourLevelAction_New);
		TecUtilArgListAppendInt(TempArgList, SV_CONTOURGROUP, (LgIndex_t)(i+1));
		TecUtilArgListAppendInt(TempArgList, SV_NUMVALUES, DataCounts[i]);
		TecUtilArgListAppendArray(TempArgList, SV_RAWDATA, DataArray[i]);
		TecUtilContourLevelX(TempArgList);
		TecUtilArgListClear(TempArgList);

		//if (i > 1)
		//{
		//	TecUtilArgListAppendInt(TempArgList, SV_CONTOURLEVELACTION, ContourLevelAction_ResetToNice);
		//	TecUtilArgListAppendInt(TempArgList, SV_CONTOURGROUP, (LgIndex_t)(i+1));
		//	TecUtilArgListAppendInt(TempArgList, SV_APPROXNUMVALUES, (LgIndex_t)200);
		//	TecUtilContourLevelX(TempArgList);
		//	TecUtilArgListClear(TempArgList);
		//}
	}

	for (int i = 1 ; i <= 8 ; i++)
	{
		TecUtilArgListAppendInt(TempArgList, SV_CONTOURLABELACTION, ContourLabelAction_DeleteAll);
		TecUtilArgListAppendInt(TempArgList, SV_CONTOURGROUP, (LgIndex_t)i);
		TecUtilContourLabelX(TempArgList);
		TecUtilArgListClear(TempArgList);
		TecUtilArgListAppendString(TempArgList, SV_P1, SV_GLOBALCONTOUR);
		TecUtilArgListAppendString(TempArgList, SV_P2, SV_LEGEND);
		TecUtilArgListAppendString(TempArgList, SV_P3, SV_SHOW);
		TecUtilArgListAppendInt(TempArgList, SV_OFFSET1, (LgIndex_t)i);
		TecUtilArgListAppendArbParam(TempArgList, SV_IVALUE, FALSE);
		TecUtilStyleSetLowLevelX(TempArgList);
		TecUtilArgListClear(TempArgList);
	}

	TecUtilArgListDealloc(&TempArgList);
}	//	void SetupMultiColorContourLevels(std::vector<Set_pa> ZoneSets)


void DeallocSets(std::vector<Set_pa> &ZoneSets)
{
	const size_t NumOfSets = ZoneSets.size();
	for (int i = 0 ; i < NumOfSets ; i++)
	{
		TecUtilSetDealloc(&ZoneSets[i]);
	}
}	//	void DeallocSets(std::vector<Set_pa> ZoneSets)


/*
	Deletes all zones that may have been created by chembond.
	Returns number of deleted zones (including if 0) if successful, -1 if not.
*/
LgIndex_t	ClearChemBondZones(Set_pa RefinedZoneSet)
{
	Boolean_t	IsOk = TRUE;
	LgIndex_t	NumDeletedZones = 0;
	EntIndex_t	NumberOfZones = TecUtilDataSetGetNumZones();
	if (NumberOfZones > 1)
	{
		Set_pa	DeleteZoneSets = TecUtilSetAlloc(TRUE);
		char	*ZoneName = NULL;
		char	ZoneNameList[][5] =	{"Crit",
									"Atom",
									"Bond",
									"Ring",
									"Cage",
									"Surf",
									"Iso:",
									"Sadd"};

		for (EntIndex_t CurrentZone = 1 ; IsOk && CurrentZone <= NumberOfZones ; CurrentZone++)
		{
			AuxData_pa	AuxDataRef = TecUtilAuxDataZoneGetRef(CurrentZone);

			if (AuxDataRef != NULL)
			{
				if (TecUtilZoneGetName(CurrentZone, &ZoneName))
				{
					Boolean_t	IsFound = FALSE;
					for (int i = 0 ; !IsFound && IsOk && i < 8 ; i++)
					{
						if (strncmp(ZoneName, ZoneNameList[i], 4) == 0)
						{
							IsOk = TecUtilSetAddMember(DeleteZoneSets, CurrentZone, TRUE);
							NumDeletedZones++;
							IsFound = TRUE;
						}
					}
				}
			}
		}
		SetIndex_t RefinedZoneCount = TecUtilSetGetMemberCount(RefinedZoneSet);
		if (RefinedZoneCount > 0 && NumDeletedZones > 0)
			for (SetIndex_t i = RefinedZoneCount ; i > RefinedZoneCount - NumDeletedZones ; i--)
				TecUtilSetRemoveMember(RefinedZoneSet, i);

		if (IsOk) IsOk = TecUtilDataSetDeleteZone(DeleteZoneSets);
		TecUtilSetDealloc(&DeleteZoneSets);
	}
	if (!IsOk) NumDeletedZones = -1;
	return NumDeletedZones;
}	//	Boolean_t	ClearChemBondZones()


Boolean_t	CreateGradVectGradMagVars()
{
	if (TecUtilDataSetGetNumVars() < 4)
	{
		TecUtilDialogErrMsg("The current data set does not have the correct variables (3 spatial {x,y,z} or {i,j,k} and 1 for electron charge density.");
		return 1;
	}

	Boolean_t IsOk = TRUE;

	IsOk = TecUtilMacroExecuteCommand("$!EXTENDEDCOMMAND\nCOMMANDPROCESSORID = 'CFDAnalyzer4'\nCOMMAND = 'SetFieldVariables ConvectionVarsAreMomentum=\'T\' UVar=0 VVar=0 WVar=0 ID1=\'Density\' Variable1=4 ID2=\'NotUsed\' Variable2=0'");

	return IsOk;
}