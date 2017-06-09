//#include "meshdefs.h"
#include <cstdio>
#include <time.h>

#include "meshgen2d_sphere.h"

// ==============================================================================
//
//  Copyright J.A. Rossmanith and the DoGPack Team
//
//  This software is made available for research and instructional use only.
//  You may copy and use this software without charge for these non-commercial
//  purposes, provided that the copyright notice and associated text is
//  reproduced on all copies.  For all other uses (including distribution of
//  modified versions), please contact the author at the address given below.
//
//  *** This software is made available "as is" without any assurance that it
//  *** will work for your purposes.  The software may in fact have defects, so
//  *** use the software at your own risk.
//
// ==============================================================================

int main(int argc, char* argv[])
{
	//
	// NOTE: You should not have to modify this part of the code.
	//       To change parameters, modify the following files:
	//            1. input2D.data
	//            2. SignedDistance.cpp
	//            3. GridSpacing.cpp
	//
  
	// Get current time
	double time1 = time(NULL);

	// Call the mesh generator
//    int meshgen2D_sphere(int argc,char**argv);
	int m = meshgen2D_sphere(argc, argv);
	
	// Get current time
	double time2 = time(NULL);

	// Output elapsed time
	printf(" Total elapsed time in seconds = %11.5e\n\n",
	   time2-time1);

	return m;
}
