#include "Point3d.h"


Point3d::Point3d(double _x, double _y, double _z)
{
	x = _x;
	y = _y;
	z = _z;
}

Point3d::Point3d(CSM_Vec3_s point)
{
	x = point.X;
	y = point.Y;
	z = point.Z;
}

Point3d::~Point3d()
{
}

void Point3d::addConnection(int connection)
{
	adjacentPoints.push_back(connection);
}
