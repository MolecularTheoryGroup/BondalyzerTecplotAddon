#pragma once
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "Triangle.h"
#include "Point3d.h"

using namespace std;

class MeshData
{
public:
	MeshData(string filename);
	
	void read_data(std::string filename);
	bool is_inside(Point3d p);

	void test_build_tets();

	void debug_print_data();
	void debug_print_points_from_vector();
	void debug_print_triangles();
	void debug_print_connections();
	void debug_save_tecplot_zone(int IterNum) const;


	~MeshData();

	int vtx_count;
	double* point_coords;	//flat array containing coordinates of points (x, y, z)

	vector<Point3d> points;
	vector<Triangle> triangles;
};

