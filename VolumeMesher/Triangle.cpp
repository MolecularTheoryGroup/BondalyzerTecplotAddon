#include <string>
#include <vector>

#include "TECADDON.h"

#include "Triangle.h"

using namespace std;

Triangle::Triangle()
{
}
Triangle::~Triangle()
{
	//dealloc centroid?
}
//Triangle::Triangle(Point3d* _a, Point3d* _b, Point3d* _c)
//{
//	a = _a;
//	b = _b;
//	c = _c;
//}
Triangle::Triangle(int _a, int _b, int _c, std::vector<Point3d> pointList)
{
	aIdx = _a;
	bIdx = _b;
	cIdx = _c;

	a = &pointList[_a];
	b = &pointList[_b];
	c = &pointList[_c];

	T_P0 = CSM_Vec3_s(a->x, a->y, a->z);
	T_P1 = CSM_Vec3_s(b->x, b->y, b->z);
	T_P2 = CSM_Vec3_s(c->x, c->y, c->z);

	centroid = CSM_Vec3_s(0, 0, 0);
	findCentroid();
	normal = findNormal();
	isFront = true;
}

CSM_Vec3_s Triangle::findNormal()
{
	CSM_Vec3_s u = T_P1 - T_P0;
	CSM_Vec3_s v = T_P2 - T_P0;
	CSM_Vec3_s n = u.Cross(v);
	return n;
}


void Triangle::findCentroid()
{
	centroid.X = (a->x + b->x + c->x) / 3;
	centroid.Y = (a->y + b->y + c->y) / 3;
	centroid.Z = (a->z + b->z + c->z) / 3;
}
int Triangle::intersectWithRay(Point3d point1, Point3d point2) const
{
	//convert to CSM_Vec datatypes
	CSM_Vec3_s R_P0(point1.x, point1.y, point1.z);
	CSM_Vec3_s R_P1(point2.x, point2.y, point2.z);

	//actual calculations begin here
	CSM_Vec3_s u = T_P1 - T_P0;
	CSM_Vec3_s v = T_P2 - T_P0;
	CSM_Vec3_s n = u.Cross(v);
	CSM_Vec3_s D = R_P1 - R_P0;
	CSM_Vec3_s w0 = R_P0 - T_P0;

	//TODO: check for degeneracy 

	//Check for ray-plane parallelism
	if (abs(n.Dot(D)) < 0.001)
		return -1;

	//Check if ray goes towards triangle plane
	double r = ((-n.Dot(w0))) / (n.Dot(D));
	if (r < 0)
		return -2;

	//Calc intersect or R and the plane
	CSM_Vec3_s P = R_P0 + (D * r);

	CSM_Vec3_s w = P - T_P0;
	double uv = u.Dot(v);
	double wv = w.Dot(v);
	double vv = v.Dot(v);
	double wu = w.Dot(u);
	double uu = u.Dot(u);

	double s = ((uv*wv) - (vv*wu)) / ((uv * uv) - (uu * vv));
	double t = ((uv*wu) - (uu*wv)) / ((uv * uv) - (uu * vv));

	//TODO : boundary cases
	if (s >= 0 && t >= 0 && s + t <= 1)
		return 1;
	else
		return 0;

}


/*
 *	Saves the triangle as a finite element zone whose name is "Tri: " and whatever you give it for the ZoneName.
 *	ColorNum changes the color. There's 65 possible colors and you can get to them all by just iterating
 *	ColorNum.
 */
const Boolean_t Triangle::debug_save_tecplot_zone(const string & ZoneName, const int & ColorNum) const{
	Boolean_t IsOk = (aIdx >= 0 && bIdx >= 0 && cIdx >= 0);
	int ZoneNum = -1;

	int NumZones, NumVars;

	IsOk = TecUtilDataSetGetInfo(NULL, &NumZones, &NumVars);

	vector<FieldDataType_e> DataTypes(NumVars, FieldDataType_Double);


	if (IsOk){
		IsOk = TecUtilDataSetAddZone(string("Tri:  " + ZoneName).c_str(),
			3,
			1,
			0,
			ZoneType_FETriangle,
			DataTypes.data());
	}

	Set_pa TmpSet = TecUtilSetAlloc(FALSE);

	if (IsOk){
		ZoneNum = TecUtilDataSetGetNumZones();
		TecUtilSetAddMember(TmpSet, ZoneNum, FALSE);
		TecUtilZoneSetActive(TmpSet, AssignOp_MinusEquals);
		TecUtilZoneSetScatter(SV_SHOW, TmpSet, 0.0, FALSE);

		TecUtilZoneSetMesh(SV_SHOW, TmpSet, 0.0, TRUE);
		TecUtilZoneSetMesh(SV_COLOR, TmpSet, 0.0, static_cast<ColorIndex_t>(ColorNum % 64));
		TecUtilZoneSetMesh(SV_LINETHICKNESS, TmpSet, 3.0, 0);

		ArgList_pa CurrentArgList = TecUtilArgListAlloc();
		TecUtilArgListAppendString(CurrentArgList, SV_P1, SV_FIELDMAP);
		TecUtilArgListAppendString(CurrentArgList, SV_P2, SV_EFFECTS);
		TecUtilArgListAppendString(CurrentArgList, SV_P3, SV_USETRANSLUCENCY);
		TecUtilArgListAppendSet(CurrentArgList, SV_OBJECTSET, TmpSet);
		TecUtilArgListAppendArbParam(CurrentArgList, SV_IVALUE, FALSE);
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

		vector<CSM_Vec3_s> TmpVecs = { T_P0, T_P1, T_P2 };

		for (int j = 0; j < 3; ++j){
			for (int i = 0; i < 3; ++i){
				TecUtilDataValueSetByRef(SurfXYZPtrs[i], j+1, TmpVecs[j][i]);
			}
		}

		
		for (int i = 0; i < 3; ++i){
			TecUtilDataNodeSetByRef(NodeMap, 1, i + 1, i + 1);
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

	TecUtilViewNiceFit();

	return IsOk;
}