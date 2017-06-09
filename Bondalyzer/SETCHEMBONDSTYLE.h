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

#ifndef SETCHEMBONDSTYLE_H_
#define SETCHEMBONDSTYLE_H_


void ChemSysView();
//	Sets style for viewing chemical system

void SGradTopoView();

LgIndex_t ClearChemBondZones(Set_pa RefinedZoneSet);
Boolean_t CreateGradVectGradMagVars();
//	Sets style for testing and analysis of SGradTopo add-on


struct ZoneSettings 
{
	//***************************************************************
	//===============================================================
	//	Setup variables for group settings
	//===============================================================
	//***************************************************************

	//	Now all zones should be grouped accordingly into one of the sets allocated above.
	//	Next we apply the desired style settings to each zone.

	//	Declare variables for each style setting that will be changed for each zone.
	//	Defined constants for some variables are listed below the variable declaration.
	//	This makes for easy adjustment of this macro.

	Boolean_t			ZoneShow[50];

	//===============================================================
	//	Mesh Settings
	//===============================================================

	Boolean_t			MeshShow[50];
	MeshType_e			MeshType[50];	
	/*	
		MeshType_Wireframe
		MeshType_HiddenLine	
		MeshType_Overlay	
	*/

	ColorIndex_t		MeshColor[50];
	/*
		A number of color index constants are \#defined. These include Black_C, Red_C,
		Green_C, Blue_C, Cyan_C, Yellow_C, Purple_C, White_C, Custom1_C through
		Custom56_C, MultiColor_C, NoColor_C, MulitiColor2_C through MulitiColor8_C,
		RGBColor_C, and InvalidColor_C.  
		To use an RGB scheme, where the value of R, G, and B are determined
		by a field variable, apply this to the desired group manually.
		This may be added in a future version.
	
		TODO add ability to specify RBG color schemes
			(inter-CP lines should be colored according to rho maybe...
	*/

	LinePattern_e		MeshLinePttrn[50];
	/*	LinePattern_Solid,
		LinePattern_Dashed,
		LinePattern_DashDot,
		LinePattern_Dotted,
		LinePattern_LongDash,
		LinePattern_DashDotDot,
	*/

	double				MeshLinePttrnLen[50];			//	Value is in %
	double				MeshLineThickness[50];			//	Value is in %


	//===============================================================
	//	Contour Settings
	//===============================================================

	Boolean_t			ContourShow[50];
	ContourType_e		ContourType[50];
	/*
		ContourType_Lines,      
		ContourType_Flood,     
		ContourType_Overlay,      
		ContourType_AverageCell, 
		ContourType_PrimaryValue,
	*/

	ContourColoring_e	FloodColoring[50];
	/*	
		ContourColoring_RGB,
		ContourColoring_Group1,
		ContourColoring_Group2,
		ContourColoring_Group3,
		ContourColoring_Group4,
		ContourColoring_Group5,
		ContourColoring_Group6,
		ContourColoring_Group7,
		ContourColoring_Group8,
	*/

	SmInteger_t			ContourLineGroup[50];
	/*
		Similiar to FloodColoring; takes values form the specified variable,
		and so has the same dilemma if files don't have consistent numbering
		for variables.
	*/

	ColorIndex_t		ContourLineColor[50];
	/*
		A number of color index constants are \#defined. These include Black_C, Red_C,
		Green_C, Blue_C, Cyan_C, Yellow_C, Purple_C, White_C, Custom1_C through
		Custom56_C, MultiColor_C, NoColor_C, MulitiColor2_C through MulitiColor8_C,
		RGBColor_C, and InvalidColor_C.  
		To use an RGB scheme, where the value of R, G, and B are determined
		by a field variable, apply this to the desired group manually.
		This may be added in a future version.
	*/

	LinePattern_e		ContourLinePttrn[50];
	/*	LinePattern_Solid,
		LinePattern_Dashed,
		LinePattern_DashDot,
		LinePattern_Dotted,
		LinePattern_LongDash,
		LinePattern_DashDotDot,
	*/

	double				ContourLinePttrnLen[50];			//	Value is in %
	double				ContourLineThickness[50];		//	Value is in %

	Boolean_t			ContourLighting[50];

	//===============================================================
	//	Vector Settings
	//===============================================================
	
	Boolean_t			VectorShow[50];
	VectorType_e		VectorType[50];
	/*
		VectorType_TailAtPoint,
		VectorType_HeadAtPoint, 
		VectorType_MidAtPoint,  
		VectorType_HeadOnly,
	*/

	ArrowheadStyle_e	VectorArrowType[50];
	/*
		ArrowheadStyle_Plain,
		ArrowheadStyle_Filled,
		ArrowheadStyle_Hollow,
	*/
	
	ColorIndex_t		VectorLineColor[50];
	/*
		A number of color index constants are \#defined. These include Black_C, Red_C,
		Green_C, Blue_C, Cyan_C, Yellow_C, Purple_C, White_C, Custom1_C through
		Custom56_C, MultiColor_C, NoColor_C, MulitiColor2_C through MulitiColor8_C,
		RGBColor_C, and InvalidColor_C.  
		To use an RGB scheme, where the value of R, G, and B are determined
		by a field variable, apply this to the desired group manually.
		This may be added in a future version.
	*/

	Boolean_t			VectorTangent[50];
	LinePattern_e		VectorLinePttrn[50];
	/*	LinePattern_Solid,
		LinePattern_Dashed,
		LinePattern_DashDot,
		LinePattern_Dotted,
		LinePattern_LongDash,
		LinePattern_DashDotDot,
	*/

	double				VectorLinePttrnLen[50];			//	Value is in %
	double				VectorLineThickness[50];			//	Value is in %

	//===============================================================
	//	Scatter Settings
	//===============================================================

	Boolean_t			ScatterShow[50];
	GeomShape_e			ScatterShape[50];
	/*
		GeomShape_Square,
		GeomShape_Del,
		GeomShape_Grad,
		GeomShape_RTri,
		GeomShape_LTri,
		GeomShape_Diamond,
		GeomShape_Circle,
		GeomShape_Cube,
		GeomShape_Sphere,
		GeomShape_Octahedron,
		GeomShape_Point,
	*/

	ColorIndex_t		ScatterOutColor[50];
	ColorIndex_t		ScatterInColor[50];
	/*
		A number of color index constants are \#defined. These include Black_C, Red_C,
		Green_C, Blue_C, Cyan_C, Yellow_C, Purple_C, White_C, Custom1_C through
		Custom56_C, MultiColor_C, NoColor_C, MulitiColor2_C through MulitiColor8_C,
		RGBColor_C, and InvalidColor_C.  
		To use an RGB scheme, where the value of R, G, and B are determined
		by a field variable, apply this to the desired group manually.
		This may be added in a future version.
	*/

	FillMode_e			ScatterFillMode[50];
	/*
		FillMode_None,
		FillMode_UseSpecificColor,
		FillMode_UseLineColor,
		FillMode_UseBackgroundColor,
	*/

	double				ScatterSize[50];					//	value is in %
	double				ScatterLineThickness[50];		//	value is in %

	//===============================================================
	//	Shade Settings
	//===============================================================

	Boolean_t			ShadeShow[50];
	ColorIndex_t		ShadeColor[50];
	/*
		A number of color index constants are \#defined. These include Black_C, Red_C,
		Green_C, Blue_C, Cyan_C, Yellow_C, Purple_C, White_C, Custom1_C through
		Custom56_C, MultiColor_C, NoColor_C, MulitiColor2_C through MulitiColor8_C,
		RGBColor_C, and InvalidColor_C.  
		To use an RGB scheme, where the value of R, G, and B are determined
		by a field variable, apply this to the desired group manually.
		This may be added in a future version.
	*/

	Boolean_t			ShadeLighting[50];

	//===============================================================
	//	Edge Settings
	//===============================================================

	Boolean_t			EdgeShow[50];
	EdgeType_e			EdgeType[50];
	/*
		EdgeType_Borders,
		EdgeType_Creases,
		EdgeType_BordersAndCreases,
	*/

	BorderLocation_e	EdgeIBorder[50];
	BorderLocation_e	EdgeJBorder[50];
	BorderLocation_e	EdgeKBorder[50];
	/*
		BorderLocation_None, 
		BorderLocation_Min,  
		BorderLocation_Max, 
		BorderLocation_Both,
	*/

	ColorIndex_t		EdgeColor[50];
	/*
		A number of color index constants are \#defined. These include Black_C, Red_C,
		Green_C, Blue_C, Cyan_C, Yellow_C, Purple_C, White_C, Custom1_C through
		Custom56_C, MultiColor_C, NoColor_C, MulitiColor2_C through MulitiColor8_C,
		RGBColor_C, and InvalidColor_C.  
		To use an RGB scheme, where the value of R, G, and B are determined
		by a field variable, apply this to the desired group manually.
		This may be added in a future version.
	*/

	double				EdgeLineThickness[50];		//	value is in %

	//===============================================================
	//	Points Settings
	//===============================================================

	PointsToPlot_e		PointsToPlot[50];
	/*
		PointsToPlot_SurfaceNodes, 
		PointsToPlot_AllNodes,     
		PointsToPlot_SurfaceCellCenters,
		PointsToPlot_AllCellCenters,
		PointsToPlot_AllConnected,
	*/

	LgIndex_t			PointsIndexSkipI[50];
	LgIndex_t			PointsIndexSkipJ[50];
	LgIndex_t			PointsIndexSkipK[50];
	/*
		1 is the same as "No Skip" in the GUI.
	*/

	//===============================================================
	//	Surfaces Settings
	//===============================================================

	SurfacesToPlot_e	SurfacesToPlot[50];
	/*
		SurfacesToPlot_BoundaryFaces,
		SurfacesToPlot_ExposedCellFaces,
		SurfacesToPlot_IPlanes,
		SurfacesToPlot_JPlanes,
		SurfacesToPlot_KPlanes,
		SurfacesToPlot_IJPlanes,
		SurfacesToPlot_JKPlanes,
		SurfacesToPlot_IKPlanes,
		SurfacesToPlot_IJKPlanes,
		SurfacesToPlot_All,
		SurfacesToPlot_None,
	*/

	LgIndex_t			RangeBeginI[50];
	LgIndex_t			RangeEndI[50];
	LgIndex_t			RangeSkipI[50];

	LgIndex_t			RangeBeginJ[50];
	LgIndex_t			RangeEndJ[50];
	LgIndex_t			RangeSkipJ[50];

	LgIndex_t			RangeBeginK[50];
	LgIndex_t			RangeEndK[50];
	LgIndex_t			RangeSkipK[50];
	/*
		0 is the same as maximum, which sets the end range to the maximum
		index for the dimension. The skip number is to skip every # planes,
		so 1 is to not skip any planes, (show every 1 plane).
	*/

	//===============================================================
	//	Volume Settings
	//===============================================================

	Boolean_t			StreamtracesShow[50];
	Boolean_t			IsoSurfacesShow[50];
	Boolean_t			SlicesShow[50];

	//===============================================================
	//	Effects Settings
	//===============================================================

	Boolean_t			UseValueBlanking[50];
	//Boolean_t			UseClipPlanes[50];
	/*
		I have no idea what's going on here. Tecplot uses a boolean
		for this value, but according to the code generator, the values
		are ~10^8, and when I cycle back and fourth between "All" and
		"None" the values keep changing, so I can't predict what is 
		correct here...
	*/

	Boolean_t			UseSurfaceTranslucency[50];
	SmInteger_t			SurfaceTranslucency[50];		//	value is in %
	LightingEffect_e	LightingEffect[50];
	/*
		LightingEffect_Paneled,
		LightingEffect_Gouraud,
		LightingEffect_None,
	*/
};	//	struct ZoneSettings


#endif