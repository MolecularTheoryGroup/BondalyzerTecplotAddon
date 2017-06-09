// ===================================================================
// structs.h
// structs header file
// ===================================================================

#ifndef STRUCTS_H
#define STRUCTS_H

#include <cmath>

// ===================================================================
struct point
{
  double x;
  double y;
  double z;

  point(double xx, double yy, double zz){
	  x = xx;
	  y = yy;
	  z = zz;
  }
  point(){
	  x = y = z = 0.0;
  }

  double magsqr(point op){
	  return pow(x-op.x,2) + pow(y-op.y,2) + pow(z-op.z,2);
  }

  double & operator[](int i){
  	  switch (i){
  	  case 0:
  		  return x;
  		  break;
  	  case 1:
  		  return y;
  		  break;
  	  default:
  		  return z;
  		  break;
  	  }
  }
};
// ===================================================================
struct triangle
{
  int n1;
  int n2;
  int n3;
  int ctr;

  int & operator[](int i){
	  switch (i){
	  case 0:
		  return n1;
		  break;
	  case 1:
		  return n2;
		  break;
	  case 2:
		  return n3;
		  break;
	  default:
		  return ctr;
		  break;
	  }
  }
};
// ===================================================================
struct bar
{
  int n1;
  int n2;
  int bt;
  int bt2;
};
// ===================================================================
struct tetra
{
  int n1;
  int n2;
  int n3;
  int n4;
};
// ===================================================================
struct cart2d
{
  double dx;
  double dy;
  double xlow;
  double ylow;
  double xhigh;
  double yhigh;
  int mx;
  int my;
};
// ===================================================================
struct cart3d
{
  double dx;
  double dy;
  double dz;
  double xlow;
  double ylow;
  double zlow;
  double xhigh;
  double yhigh;
  double zhigh;
  int mx;
  int my;
  int mz;
};
// ===================================================================
#endif
