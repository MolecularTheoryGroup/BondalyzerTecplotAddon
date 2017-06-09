#include "MeshData.h"
#include <iostream>
#include <fstream>
#include <string>

#include "TECADDON.h"
#include "CSMDATATYPES.h"

using namespace std;

/* Constructor */
MeshData::MeshData(string filename)
{
	read_data(filename);
}

/* Destructor */
MeshData::~MeshData(){}

/* Reads data from filename into memory */
void MeshData::read_data(string filename)
{
	ifstream dataFile(filename);
	cout << "Here's " << filename << endl;

	/* Read in the number of vertices from the first line in the file */
	string vtxCountString;
	getline(dataFile, vtxCountString);
	vtx_count = atoi(vtxCountString.c_str());

	/* Allocate space for the points in the data array */
	point_coords = new double[vtx_count * 3];

	/* Read in coordinates of vertices */
	for (int i = 0; i < vtx_count * 3; i++)
	{
		string coordString;
		getline(dataFile, coordString, ' ');

		/* This gets rid of strings consisting of just spaces or newlines*/
		if (coordString.length() < 2)
		{
			i--;
			continue;
		}
		
		point_coords[i] = atof(coordString.c_str());
	}

	/* Save coordinates of vertices into Point3d vector */
	for (int i = 0; i < vtx_count * 3; i += 3)
	{
		Point3d p(point_coords[i], point_coords[i + 1], point_coords[i + 2]);
		points.push_back(p);
	}

	/*Deallocate space for points */
	//delete(point_coords);			//TODO : put this line back in


	/*Read in triangles into triangle vector*/
	string triangleString;
	getline(dataFile, triangleString);

	while (triangleString.length() > 0)
	{
		triangleString = triangleString.substr(triangleString.find_first_not_of(" "), triangleString.find_last_not_of(" ") + 1);
		//cout << triangleString << endl;
		
		int secondNumIdx = triangleString.find_first_of(" ");
		int thirdNumIdx = triangleString.find_last_of(" ");
		
		int firstNum = stoi(triangleString.substr(0, secondNumIdx)) - 1;
		int secondNum = stoi(triangleString.substr(secondNumIdx, thirdNumIdx)) - 1;
		int thirdNum = stoi(triangleString.substr(thirdNumIdx, triangleString.length() - 1)) - 1;

		Triangle t(firstNum, secondNum, thirdNum, points);
		triangles.push_back(t);

		getline(dataFile, triangleString);
	}

	/*Convert the triangles into an adjacency list in each vertex*/
	for (Triangle t : triangles)
	{
		int a = t.aIdx;
		int b = t.bIdx;
		int c = t.cIdx;

		points.at(a).addConnection(b);
		points.at(a).addConnection(c);
		points.at(b).addConnection(a);
		points.at(b).addConnection(c);
		points.at(c).addConnection(a);
		points.at(c).addConnection(b);
	}
}

void MeshData::test_build_tets()
{
	int cTri = triangles.size();
	for (int i = 0; i < cTri; i++)
	//for (int i = 0; i < ; i++)
	{
		//get the triangle we're working on
		Triangle t = triangles[i];
		//set the magnitude of it's normal to something useful

		if (i == 261)
		{
			cout << "BREAK HERE";
		}

		t.normal = t.normal.Normalize();
		t.normal = t.normal * ((t.T_P0 - t.T_P1).Magnitude());

		//make a new point along that normal
		CSM_Vec3_s newPointVec = t.centroid + t.normal;
		Point3d newPoint(newPointVec);

		//check if this new point is inside the volume
		//if it isn't, reverse the normal and remake the point
		if (!is_inside(newPoint))
		{
			t.normal = t.normal * -1;
			newPointVec = t.centroid + t.normal;
			newPoint = Point3d(newPointVec);
		}
		//also check the distance of the point to the center of the sphere for debug

		//update the point list/connectivity info
		int newPointIdx = points.size();
		newPoint.addConnection(t.aIdx);
		newPoint.addConnection(t.bIdx);
		newPoint.addConnection(t.cIdx);
		points.push_back(newPoint);

		points[t.aIdx].addConnection(newPointIdx);
		points[t.bIdx].addConnection(newPointIdx);
		points[t.cIdx].addConnection(newPointIdx);


		
		//make and add the new triangles, and remove the old triangle from the front
		Triangle newT1(t.aIdx, t.bIdx, newPointIdx, points);
		Triangle newT2(t.bIdx, t.cIdx, newPointIdx, points);
		Triangle newT3(t.cIdx, t.aIdx, newPointIdx, points);
		
		t.isFront = false;
		newT1.isFront = true;
		newT2.isFront = true;
		newT2.isFront = true;
		
		triangles.push_back(newT1);
		triangles.push_back(newT2);
		triangles.push_back(newT3);
	}
}

bool MeshData::is_inside(Point3d p)
{
	int intersectCount = 0;
	Point3d farPoint(-1000, -1000, -1000); //todo: verify this makes sense as the far point of the intersection ray
	for (const Triangle & t : triangles)
	{
		if (t.isFront && t.intersectWithRay(p, farPoint) == 1)
			intersectCount++;
	}
	if (intersectCount % 2 == 0) //even number of intersections -> outside
		return false;
	else
		return true;
}

void MeshData::debug_print_data()
{
	for (int i = 0; i < vtx_count * 3; i += 3)
	{
		//cout.precision(15);
		std::cout << point_coords[i] << "," << point_coords[i + 1] << ", " << point_coords[i + 2] << endl;
	}
}

void MeshData::debug_print_points_from_vector()
{
	for (const Point3d & point : points)
	{
		cout.precision(15);
		cout << point.x << "," << point.y << "," << point.z << endl;
	}
}

void MeshData::debug_print_triangles()
{
	for (Triangle t : triangles)
	{
		cout << t.aIdx + 1 << "," << t.bIdx + 1 << "," << t.cIdx + 1 << endl;
	}
	cout << "End triangle list" << endl;
}

void MeshData::debug_print_connections()
{
	for (int i = 0; i < points.size(); i++)
	{
		Point3d p = points.at(i);
		cout << i + 1 << " : ";;
		for (int connection : p.adjacentPoints)
		{
			cout << connection + 1 << ",";
		}
		cout << endl;
	}
}

/*
 *	Saves the mesh as a triangle zone.
 *	IterNum < 0 means it's the initial mesh
 *	> 0 means it's an iteration of the volume
 *	meshing algorithm.
 *	Needs to be triangles, and triangles must be sane.
 */
void MeshData::debug_save_tecplot_zone(int IterNum) const
{
	Boolean_t IsOk = TRUE;
	int ZoneNum = -1;
	
	if (IterNum < 0){
		vector<string> XYZNames = { "X", "Y", "Z" };
		StringList_pa XYZVarNames = TecUtilStringListAlloc();
		for (auto i = XYZNames.cbegin(), End = XYZNames.cend(); i != End; i++)
			TecUtilStringListAppendString(XYZVarNames, i->c_str());

		string DataSetName = "VolumeMesherTest";

		vector<FieldDataType_e> DSDataTypes(3, FieldDataType_Double);

		IsOk = TecUtilDataSetCreate(DataSetName.c_str(), XYZVarNames, TRUE);
		if (IsOk){
			TecUtilStringListDealloc(&XYZVarNames);
			IsOk = TecUtilDataSetAddZone("Initial Mesh", points.size(), triangles.size(), 0, ZoneType_FETriangle, DSDataTypes.data());
		}
		if (IsOk){
			TecUtilFrameSetPlotType(PlotType_Cartesian3D);

			ArgList_pa argList = TecUtilArgListAlloc();
			TecUtilArgListClear(argList);
			TecUtilArgListAppendString(argList, SV_P1, SV_FIELDLAYERS);
			TecUtilArgListAppendString(argList, SV_P2, SV_SHOWMESH);
			TecUtilArgListAppendArbParam(argList, SV_IVALUE, TRUE);
			TecUtilStyleSetLowLevelX(argList);

			TecUtilArgListClear(argList);
			TecUtilArgListAppendString(argList, SV_P1, SV_FIELDLAYERS);
			TecUtilArgListAppendString(argList, SV_P2, SV_SHOWEDGE);
			TecUtilArgListAppendArbParam(argList, SV_IVALUE, FALSE);
			TecUtilStyleSetLowLevelX(argList);

			TecUtilArgListDealloc(&argList);
		}
	}
	else{
		int NumZones, NumVars;

		IsOk = TecUtilDataSetGetInfo(NULL, &NumZones, &NumVars);

		vector<FieldDataType_e> DataTypes(NumVars, FieldDataType_Double);


		if (IsOk){
			IsOk = TecUtilDataSetAddZone(string("Mesher " + to_string(IterNum)).c_str(),
				points.size(),
				triangles.size(),
				0,
				ZoneType_FETriangle,
				DataTypes.data());
		}
	}

	Set_pa TmpSet = TecUtilSetAlloc(FALSE);

	if (IsOk){
		ZoneNum = TecUtilDataSetGetNumZones();
		TecUtilSetAddMember(TmpSet, ZoneNum, FALSE);
		if (IterNum < 0)
			TecUtilZoneSetActive(TmpSet, AssignOp_PlusEquals);
		else
			TecUtilZoneSetActive(TmpSet, AssignOp_MinusEquals);
		TecUtilZoneSetScatter(SV_SHOW, TmpSet, 0.0, FALSE);

		ArgList_pa CurrentArgList = TecUtilArgListAlloc();
		TecUtilArgListAppendString(CurrentArgList, SV_P1, SV_FIELDMAP);
		TecUtilArgListAppendString(CurrentArgList, SV_P2, SV_EFFECTS);
		TecUtilArgListAppendString(CurrentArgList, SV_P3, SV_USETRANSLUCENCY);
		TecUtilArgListAppendSet(CurrentArgList, SV_OBJECTSET, TmpSet);
		TecUtilArgListAppendArbParam(CurrentArgList, SV_IVALUE, TRUE);
		TecUtilStyleSetLowLevelX(CurrentArgList);
		TecUtilArgListDealloc(&CurrentArgList);

		TecUtilSetDealloc(&TmpSet);
	}

	vector<EntIndex_t> XYZVarNums = { 1, 2, 3 };
// 	vector<EntIndex_t> XYZVarNums(3);
// 	TecUtilAxisGetVarAssignments(&XYZVarNums[0], &XYZVarNums[1], &XYZVarNums[2]);
// 	for (int i = 0; i < 3 && IsOk; ++i)
// 		IsOk = (XYZVarNums[i] > 0);

	FieldData_pa SurfXYZPtrs[3];

	TecUtilDataLoadBegin();

	for (int i = 0; i < 3 && IsOk; ++i){
		SurfXYZPtrs[i] = TecUtilDataValueGetWritableRef(ZoneNum, XYZVarNums[i]);
		IsOk = VALID_REF(SurfXYZPtrs[i]);
	}

	if (IsOk){
		NodeMap_pa NodeMap;
		NodeMap = TecUtilDataNodeGetWritableRef(ZoneNum);
		IsOk = VALID_REF(NodeMap);


		for (LgIndex_t NodeNum = 0; NodeNum < points.size(); ++NodeNum){
			for (int i = 0; i < 3; ++i){
				TecUtilDataValueSetByRef(SurfXYZPtrs[i], NodeNum + 1, points[NodeNum][i]);
			}
		}
		for (LgIndex_t TriNum = 0; TriNum < triangles.size(); ++TriNum){
			for (int i = 0; i < 3; ++i){
				TecUtilDataNodeSetByRef(NodeMap, TriNum + 1, i + 1, triangles[TriNum][i] + 1);
			}
		}

	}

	TecUtilDataLoadEnd();

	if (IsOk){
		Set_pa TmpSet = TecUtilSetAlloc(TRUE);

		for (int i = 1; i <= 3; ++i)
			TecUtilSetAddMember(TmpSet, i, TRUE);

		TecUtilStateChanged(StateChange_VarsAltered, reinterpret_cast<ArbParam_t>(TmpSet));

		TecUtilSetDealloc(&TmpSet);
	}


	if (IterNum < 0)
		TecUtilReset3DOrigin();

	TecUtilViewNiceFit();

	IsOk = triangles.at(0).debug_save_tecplot_zone("test save 1", 0);
	IsOk = triangles.at(1).debug_save_tecplot_zone("test save 2", 1);
	IsOk = triangles.at(2).debug_save_tecplot_zone("test save 3", 2);
	IsOk = triangles.at(3).debug_save_tecplot_zone("test save 4", 3);
	IsOk = triangles.at(4).debug_save_tecplot_zone("test save 5", 4);
}