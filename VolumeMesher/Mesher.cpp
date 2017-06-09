// Mesher.cpp : Defines the entry point for the console application.
#include<iostream>
#include<fstream>
#include<string>

#include "MeshData.h"


int main()
{
	MeshData meshData("U:\\Workspace\\small_sphere.dat");
	//MeshData meshData("C:/Users/Tim Wilson/Desktop/sphere.dat");
	//meshData.debug_print_data();
	//meshData.debug_print_points_from_vector();
	//meshData.debug_print_triangles();
	meshData.debug_print_connections();

	meshData.debug_save_tecplot_zone(0);
	meshData.debug_save_tecplot_zone(1);

	//meshData.test_build_tets();

	//meshData.debug_save_tecplot_zone(1);
}


