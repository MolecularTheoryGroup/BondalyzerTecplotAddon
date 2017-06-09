#include "meshdefs.h"
#include "mesh.h"
#include <cstdlib>

#include "screenoutput.h"

//
// Output basic mesh information to screen
//
void ScreenOutput(const mesh& Mesh)
{
  // Compute mesh quality parameters
  double totalarea = Mesh.get_area_prim(1);
  double maxarea = Mesh.get_area_prim(1);
  double minarea = Mesh.get_area_prim(1);
  for (int i=2; i<=Mesh.get_NumPhysElems(); i++)
    {
      double tmp = Mesh.get_area_prim(i);
      totalarea = totalarea + tmp;
      if (tmp < minarea)
	{ minarea = tmp; }
      if (tmp > maxarea)
	{ maxarea = tmp; }
    }

  double minAngle = 180.0;

  for (int i=1; i<=Mesh.get_NumPhysElems(); i++)
    {
      const int i1 = Mesh.get_tnode(i,1);
      const int i2 = Mesh.get_tnode(i,2);
      const int i3 = Mesh.get_tnode(i,3);

      point v12, v23, v31;
      v12.x = Mesh.get_node(i2,1) - Mesh.get_node(i1,1);
      v12.y = Mesh.get_node(i2,2) - Mesh.get_node(i1,2);
      v23.x = Mesh.get_node(i3,1) - Mesh.get_node(i2,1);
      v23.y = Mesh.get_node(i3,2) - Mesh.get_node(i2,2);
      v31.x = Mesh.get_node(i1,1) - Mesh.get_node(i3,1);
      v31.y = Mesh.get_node(i1,2) - Mesh.get_node(i3,2);
      
      double angle1 = acos((v12.x*-v31.x+v12.y*-v31.y)
			   /(sqrt(v12.x*v12.x+v12.y*v12.y)*sqrt(v31.x*v31.x+v31.y*v31.y)));
      double angle2 = acos((v23.x*-v12.x+v23.y*-v12.y)
			   /(sqrt(v23.x*v23.x+v23.y*v23.y)*sqrt(v12.x*v12.x+v12.y*v12.y)));
      double angle3 = acos((v31.x*-v23.x+v31.y*-v23.y)
			   /(sqrt(v31.x*v31.x+v31.y*v31.y)*sqrt(v23.x*v23.x+v23.y*v23.y)));
      if ((angle1*180/pi) < minAngle)
        {  minAngle = angle1*180/pi;  }
      if ((angle2*180/pi) < minAngle)
        {  minAngle = angle2*180/pi;  }
      if ((angle3*180/pi) < minAngle)
        {  minAngle = angle3*180/pi;  }
    }

  // Output summary of results to screen
  printf("\n");
  printf("  SUMMARY OF RESULTS:\n");
  printf("  -------------------\n");
  printf("          Number of Elements:  %8i\n",Mesh.get_NumElems());
  printf(" Number of Physical Elements:  %8i\n",Mesh.get_NumPhysElems());
  printf("    Number of Ghost Elements:  %8i\n",Mesh.get_NumGhostElems());
  printf("             Number of Nodes:  %8i\n",Mesh.get_NumNodes());
  printf("    Number of Physical Nodes:  %8i\n",Mesh.get_NumPhysNodes());
  printf("    Number of Boundary Nodes:  %8i\n",Mesh.get_NumBndNodes());
  printf("             Number of Edges:  %8i\n",Mesh.get_NumEdges());
  printf("    Number of Boundary Edges:  %8i\n",Mesh.get_NumBndEdges());
  printf("\n");
  printf("          Total Area Covered:  %24.16e\n",totalarea);
  printf("     Area Ratio: small/large:  %24.16e\n",minarea/maxarea);
  printf("    Angle Ratio: minAngle/60:  %24.16e\n",minAngle/60.0);
  printf("\n");
}
