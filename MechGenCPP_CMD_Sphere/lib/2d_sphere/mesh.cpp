// --------------------------------------------------------------------------
//  IMPLEMENTATION FILE (mesh.cpp)
//    2d unstructured mesh
// --------------------------------------------------------------------------

#include <cmath>
#include <cstdio>
#include <assert.h>
#include <cstdlib>

#include "constants.h"
#include "unique.h"

#include "mesh.h"

mesh::mesh(int inNumElems, 
	   int inNumPhysElems, 
	   int inNumNodes, 
	   int inNumPhysNodes, 
	   int inNumBndNodes, 
	   int inNumEdges,
	   int inNumBndEdges)
// Constructor
// POST: Creates a mesh
{
  NumElems      = inNumElems;
  NumPhysElems  = inNumPhysElems;
  NumGhostElems = NumElems - NumPhysElems;
  NumNodes      = inNumNodes;
  NumPhysNodes  = inNumPhysNodes; 
  NumBndNodes   = inNumBndNodes;
  NumEdges      = inNumEdges;
  NumBndEdges   = inNumBndEdges;

  if (NumElems<1 || NumPhysElems<1 || NumGhostElems<0 || NumNodes<1 || 
      NumPhysNodes<1 || NumBndNodes<0  || NumEdges<1)
    {
      printf("\n");
      printf(" Error in mesh constructor ... \n");
      printf("         NumElems = %i\n",NumElems);
      printf("     NumPhysElems = %i\n",NumPhysElems);
      printf("    NumGhostElems = %i\n",NumGhostElems);
      printf("         NumNodes = %i\n",NumNodes);
      printf("     NumPhysNodes = %i\n",NumPhysNodes);
      printf("      NumBndNodes = %i\n",NumBndNodes);
      printf("         NumEdges = %i\n",NumEdges);
      printf("      NumBndEdges = %i\n",NumBndEdges);
      printf("\n");
      exit(1);
    }
  
  // node: list of x & y coordinates of nodes (size = NumNodes-by-2)
  node = new dTensor2(NumNodes,2);

  // edge: list of coordinates (x1,y1) and (x2,y2) that make up edge (size = NumEdges-by-4)
  edge = new dTensor2(NumEdges,4);

  // tnode: list of nodes attached to element (size = NumNodes-by-3)
  tnode = new iTensor2(NumElems,3);

  // tedge: list of edges attached to element (size = NumElems-by-3)
  tedge = new iTensor2(NumElems,3);

  // tedge_orientation: list of orientiation (either +1 or -1) for edges attached to element (size = NumElems-by-3)
  tedge_orientation = new iTensor2(NumElems,3);

  // eelem: list of elements attached to edge (size = NumEdges-by-2)
  eelem = new iTensor2(NumEdges,2);

  // enode: list of nodes attached to edge (size = NumEdges-by-2)
  enode = new iTensor2(NumEdges,2);
  
  // list of elements on opposite side of boundary to ghost cell (size = NumGhostElems)
  int tmp;
  if (NumGhostElems==0)
    { tmp = 1; }
  else
    { tmp = NumGhostElems; }
  ghost_link = new iTensor1(tmp);

  // list of nodes on opposite side of boundary to external node
  if ((NumNodes-NumPhysNodes)<1)
    { tmp = 1; }
  else
    { tmp = NumNodes-NumPhysNodes; }
  ext_node_link = new iTensor1(tmp);

  // list of nodes that lie on boundary (size = NumBndNodes)
  if (NumBndNodes==0)
    { tmp = 1; }
  else
    { tmp = NumBndNodes; }
  bnd_node = new iTensor1(tmp);

  // list of edges that lie on boundary (size = NumBndEdges)
  if (NumBndEdges==0)
    { tmp = 1; }
  else
    { tmp = NumBndEdges; }
  bnd_edge = new iTensor1(tmp);

  // element areas (size = NumElems)
  area_prim = new dTensor1(NumElems);

  // dual element areas (size = NumNodes)
  area_dual = new dTensor1(NumNodes);

  // Jacobian matrix (size = NumElems-by-2-by-2)
  Jmat = new dTensor3(NumElems,2,2);
  
  // list of number of elements per node
  NumElemsPerNode = new iTensor1(NumNodes);

  // adjacent: list of elements that are adjacent to current element
  adjacent = new iTensor2(NumElems,3);

  // Skel2Elem: matrix to convert a vector defined on the
  // mesh skeleton to an element-center quantity
  skel2elem_on = false;
  Skel2Elem = new dTensor3(NumElems,20,20);

  // KMI: matrix to convert a vector node information to
  // edge-based quantities
  kmi_on = false;
  KMI = new dTensor3(NumEdges,3,6);

  // Submesh information
  is_submesh = false;
}

mesh::mesh(const mesh& amesh)
// Copy constructor
// POST: New mesh created with size and contents same as amesh
{
  NumElems      = amesh.NumElems;
  NumPhysElems  = amesh.NumPhysElems;
  NumGhostElems = amesh.NumGhostElems;
  NumNodes      = amesh.NumNodes;
  NumPhysNodes  = amesh.NumPhysNodes;
  NumBndNodes   = amesh.NumBndNodes;
  NumEdges      = amesh.NumEdges;
  NumBndEdges   = amesh.NumBndEdges;

  for (int i=1; i<=NumNodes; i++)
    for (int j=1; j<=2; j++)
      {
	node->set(i,j, amesh.node->get(i,j) );
      }

  for (int i=1; i<=NumEdges; i++)
    for (int j=1; j<=4; j++)
      {
	edge->set(i,j, amesh.edge->get(i,j) );
      }

  for (int i=1; i<=NumElems; i++)
    for (int j=1; j<=3; j++)
      {
	tnode->set(i,j, amesh.tnode->get(i,j) );
      }

  for (int i=1; i<=NumElems; i++)
    for (int j=1; j<=3; j++)
      {
	tedge->set(i,j, amesh.tedge->get(i,j) );
      }

  for (int i=1; i<=NumElems; i++)
    for (int j=1; j<=3; j++)
      {
	tedge_orientation->set(i,j, amesh.tedge_orientation->get(i,j) );
      }

  for (int i=1; i<=NumEdges; i++)
    for (int j=1; j<=2; j++)
      {
	eelem->set(i,j, amesh.eelem->get(i,j) );
      }

  for (int i=1; i<=NumEdges; i++)
    for (int j=1; j<=2; j++)
      {
	enode->set(i,j, amesh.enode->get(i,j) );
      }

  for (int i=1; i<=NumGhostElems; i++)
    {
      ghost_link->set(i, amesh.ghost_link->get(i) );
    }
  
  for (int i=1; i<=(NumNodes-NumPhysNodes); i++)
    {
      ext_node_link->set(i, amesh.ext_node_link->get(i) );
    }

  for (int i=1; i<=NumBndNodes; i++)
    {
      bnd_node->set(i, amesh.bnd_node->get(i) );
    }

  for (int i=1; i<=NumBndEdges; i++)
    {
      bnd_edge->set(i, amesh.bnd_edge->get(i) );
    }

  for (int i=1; i<=NumElems; i++)
    {
      area_prim->set(i, amesh.area_prim->get(i) );
    }

  for (int i=1; i<=NumNodes; i++)
    {
      area_dual->set(i, amesh.area_dual->get(i) );
    } 

  for (int i=1; i<=NumElems; i++)
    for (int m1=1; m1<=2; m1++)
      for (int m2=1; m2<=2; m2++)
	{
	  Jmat->set(i,m1,m2, amesh.Jmat->get(i,m1,m2) );
	}
  
  for (int i=1; i<=NumNodes; i++)
    {
      NumElemsPerNode->set(i, amesh.NumElemsPerNode->get(i) );
    }

  for (int i=1; i<=NumElems; i++)
    for (int m=1; m<=3; m++)
      {
	adjacent->set(i,m, amesh.adjacent->get(i,m) );
      } 

  for (int i=1; i<=NumEdges; i++)
    for (int m=1; m<=3; m++)
      for (int n=1; n<=6; n++)
	{
	  KMI->set(i,m,n, amesh.KMI->get(i,m,n) );
	}
  
  // Submesh information
  if (amesh.is_submesh)
    {
      is_submesh      = amesh.is_submesh;
      SubFactor       = amesh.SubFactor;
      SubNumPhysElems = amesh.SubNumPhysElems;
      SubNumPhysNodes = amesh.SubNumPhysNodes;
      SubNumBndNodes  = amesh.SubNumBndNodes;

      for (int i=1; i<=SubNumPhysNodes; i++)
	for (int j=1; j<=2; j++)
	  {
	    sub_node->set(i,j, amesh.sub_node->get(i,j) );
	  }

      for (int i=1; i<=SubNumPhysElems; i++)
	for (int j=1; j<=3; j++)
	  {
	    sub_tnode->set(i,j, amesh.sub_tnode->get(i,j) );
	  }

      for (int i=1; i<=SubNumBndNodes; i++)
        {
	  sub_bnd_node->set(i, amesh.sub_bnd_node->get(i) );
        }

      for (int i=1; i<=SubNumPhysElems; i++)
        {
	  sub_area_prim->set(i, amesh.sub_area_prim->get(i) );
        }

      for (int i=1; i<=NumPhysElems; i++)
	for (int k=1; k<=(SubFactor*SubFactor); k++)
	  {
	    elem_subs->set(i,k, amesh.elem_subs->get(i,k) );
	  }

      for (int i=1; i<=NumPhysElems; i++)
	for (int k=1; k<=(((SubFactor+1)*(SubFactor+2))/2); k++)
	  {
	    node_subs->set(i,k, amesh.node_subs->get(i,k) );
	  }

    }

}

mesh::~mesh()
// Destructor
// POST: mesh no longer exists
{
  delete node;
  delete tnode;
  delete ghost_link;
  delete ext_node_link;
  delete bnd_node;
  delete bnd_edge;
  delete area_prim;
  delete area_dual;
  delete edge;
  delete tedge;
  delete tedge_orientation;
  delete eelem;
  delete enode;
  delete Jmat;
  delete NumElemsPerNode;
  delete adjacent;
  if (skel2elem_on)
    {delete Skel2Elem;}
  if (kmi_on)
    {delete KMI;}
  if (is_submesh)
    {
      delete sub_node;
      delete sub_tnode;
      delete sub_bnd_node;
      delete sub_area_prim;
      delete elem_subs;
      delete node_subs;
    }
}

void mesh::OutputMesh(char* outputdir)
// Ouput all mesh information
{
  // run startscript
  // to create output directory if it does not exist
  // and copy data files to output directory
  char command_str[1024];
  int numchars = snprintf(command_str,1024,
			  "if test -f startscript && test -x startscript;\n"
			  "then ./startscript %s %d\n"
			  "else ${MESHGENCPP}/lib/2d_sphere/startscript %s %d\n"
			  "fi",outputdir,2,outputdir,2);
  assert_lt(numchars,1023);
  assert_gt(numchars,0);
  int exit_status = 0;
  exit_status = system(command_str);
    
  // output: constants
  char fname[1024];
  snprintf(fname,1024,"%s/mesh_params.dat",outputdir);
  FILE* write1 = fopen(fname,"w");
  fprintf(write1,"%8i : NumElems\n",      NumElems);
  fprintf(write1,"%8i : NumPhysElems\n",  NumPhysElems);
  fprintf(write1,"%8i : NumGhostElems\n", NumGhostElems);
  fprintf(write1,"%8i : NumNodes\n",      NumNodes);
  fprintf(write1,"%8i : NumPhysNodes\n",  NumPhysNodes);
  fprintf(write1,"%8i : NumBndNodes\n",   NumBndNodes);
  fprintf(write1,"%8i : NumEdges\n",      NumEdges);
  fprintf(write1,"%8i : NumBndEdges\n",   NumBndEdges);
  fprintf(write1,"%8i : is_submesh\n",    is_submesh);
  fclose(write1);

  // output: NODE
  snprintf(fname,1024,"%s/mesh_node.dat",outputdir);
  FILE* write2 = fopen(fname,"w");
  for (int i=1; i<=NumNodes; i++)
    {
      fprintf(write2,"%24.16e  %24.16e\n",
	      node->get(i,1),
	      node->get(i,2));
    }
  fclose(write2);

  // output: EDGE
  snprintf(fname,1024,"%s/mesh_edge.dat",outputdir);
  FILE* write3 = fopen(fname,"w");
  for (int i=1; i<=NumEdges; i++)
    {
      fprintf(write3,"%24.16e  %24.16e  %24.16e  %24.16e\n",
	      edge->get(i,1),edge->get(i,2),
	      edge->get(i,3),edge->get(i,4));
    }
  fclose(write3);

  // output: TNODE
  snprintf(fname,1024,"%s/mesh_tnode.dat",outputdir);
  FILE* write4 = fopen(fname,"w");
  for (int i=1; i<=NumElems; i++)
    {
      fprintf(write4,"%8i  %8i  %8i\n",
	      tnode->get(i,1),
	      tnode->get(i,2),
	      tnode->get(i,3));
    }
  fclose(write4);

  // output: TEDGE
  snprintf(fname,1024,"%s/mesh_tedge.dat",outputdir);
  FILE* write5a = fopen(fname,"w");
  for (int i=1; i<=NumElems; i++)
    {
      fprintf(write5a,"%8i  %8i  %8i\n",
	      tedge->get(i,1),
	      tedge->get(i,2),
	      tedge->get(i,3));
    }
  fclose(write5a);

  // output: TEDGE_ORIENTATION
  snprintf(fname,1024,"%s/mesh_tedge_orientation.dat",outputdir);
  FILE* write5b = fopen(fname,"w");
  for (int i=1; i<=NumElems; i++)
    {
      fprintf(write5b,"%8i  %8i  %8i\n",
	      tedge_orientation->get(i,1),
	      tedge_orientation->get(i,2),
	      tedge_orientation->get(i,3));
    }
  fclose(write5b);

  // output: EELEM
  snprintf(fname,1024,"%s/mesh_eelem.dat",outputdir);
  FILE* write6 = fopen(fname,"w");
  for (int i=1; i<=NumEdges; i++)
    {
      fprintf(write6,"%8i  %8i\n",
	      eelem->get(i,1),
	      eelem->get(i,2));
    }
  fclose(write6);

  // output: ENODE
  snprintf(fname,1024,"%s/mesh_enode.dat",outputdir);
  FILE* write7 = fopen(fname,"w");
  for (int i=1; i<=NumEdges; i++)
    {
      fprintf(write7,"%8i  %8i\n",
	      enode->get(i,1),
	      enode->get(i,2));
    }
  fclose(write7);

  // output: GHOST_LINK
  snprintf(fname,1024,"%s/mesh_ghost_link.dat",outputdir);
  FILE* write8 = fopen(fname,"w");
  for (int i=1; i<=NumGhostElems; i++)
    {
      fprintf(write8,"%8i\n",ghost_link->get(i));
    }
  fclose(write8);

  // output: EXT_NODE_LINK
  snprintf(fname,1024,"%s/mesh_ext_node_link.dat",outputdir);
  FILE* write8b = fopen(fname,"w");
  for (int i=1; i<=(NumNodes-NumPhysNodes); i++)
    {
      fprintf(write8b,"%8i\n",ext_node_link->get(i));
    }
  fclose(write8b);

  // output: BND_NODE
  snprintf(fname,1024,"%s/mesh_bnd_node.dat",outputdir);
  FILE* write9 = fopen(fname,"w");
  for (int i=1; i<=NumBndNodes; i++)
    {
      fprintf(write9,"%8i\n",bnd_node->get(i));
    }
  fclose(write9);

  // output: BND_EDGE
  snprintf(fname,1024,"%s/mesh_bnd_edge.dat",outputdir);
  FILE* write9b = fopen(fname,"w");
  for (int i=1; i<=NumBndEdges; i++)
    {
      fprintf(write9b,"%8i\n",bnd_edge->get(i));
    }
  fclose(write9b);

  // output: AREA_PRIM
  snprintf(fname,1024,"%s/mesh_area_prim.dat",outputdir);
  FILE* write10 = fopen(fname,"w");
  for (int i=1; i<=NumElems; i++)
    {
      fprintf(write10,"%24.16e\n",area_prim->get(i));
    }
  fclose(write10);

  // output: AREA_DUAL
  snprintf(fname,1024,"%s/mesh_area_dual.dat",outputdir);
  FILE* write11 = fopen(fname,"w");
  for (int i=1; i<=NumNodes; i++)
    {
      fprintf(write11,"%24.16e\n",area_dual->get(i));
    }
  fclose(write11);

  // output: JMAT
  snprintf(fname,1024,"%s/mesh_jmat.dat",outputdir);
  FILE* write12 = fopen(fname,"w");
  for (int i=1; i<=NumElems; i++)
    for (int m1=1; m1<=2; m1++)
      for (int m2=1; m2<=2; m2++)
	{
	  fprintf(write12,"%24.16e\n",Jmat->get(i,m1,m2));
	}
  fclose(write12);

  // output: NUMELEMSPERNODE
  snprintf(fname,1024,"%s/mesh_numelemspernode.dat",outputdir);
  FILE* write13 = fopen(fname,"w");
  for (int i=1; i<=NumNodes; i++)
    {
      fprintf(write13,"%8i\n",NumElemsPerNode->get(i));
    }
  fclose(write13);

  // output: ADJACENT
  snprintf(fname,1024,"%s/mesh_adjacent.dat",outputdir);
  FILE* write14 = fopen(fname,"w");
  for (int i=1; i<=NumElems; i++)
    {
      fprintf(write14,"%8i  %8i  %8i\n",
	      adjacent->get(i,1),
	      adjacent->get(i,2),
	      adjacent->get(i,3));
    }
  fclose(write14);

  // output: SKEL2ELEM
  if (skel2elem_on)
    {
      snprintf(fname,1024,"%s/mesh_skel2elem.dat",outputdir);
      FILE* write15 = fopen(fname,"w");
      for (int i=1; i<=NumElems; i++)
	for (int m=1; m<=20; m++)
	  for (int n=1; n<=20; n++)
	    {
	      fprintf(write15,"%24.16e\n",
		      Skel2Elem->get(i,m,n));
	    }    
      fclose(write15);
    }

  // output: KMI
  if (kmi_on)
    {
      snprintf(fname,1024,"%s/mesh_kmi.dat",outputdir);
      FILE* write16 = fopen(fname,"w");
      for (int i=1; i<=NumEdges; i++)
	for (int m=1; m<=3; m++)
	  for (int n=1; n<=6; n++)
	    {
	      fprintf(write16,"%24.16e\n",
		      KMI->get(i,m,n));
	    }    
      fclose(write16);
    }
    
  // sub-mesh
  if (is_submesh)
    {  OutputSubMesh(outputdir);  }
}

void mesh::InputMesh(char* inputdir)
// Input all mesh information
{
  // input: constants
  int NumElems_in,NumPhysElems_in,NumGhostElems_in; 
  int NumNodes_in,NumPhysNodes_in,NumBndNodes_in;
  int NumEdges_in,NumBndEdges_in;
  
  char fname[1024];
  char buffer[1024];
  snprintf(fname,1024,"%s/mesh_params.dat",inputdir);
  FILE* read1 = fopen(fname,"r");  

  int garbage;
  garbage=fscanf(read1,"%i",&NumElems_in);  
  assert(fgets(buffer, sizeof buffer, read1)!=NULL);
  garbage=fscanf(read1,"%i",&NumPhysElems_in);  
  assert(fgets(buffer, sizeof buffer, read1)!=NULL);
  garbage=fscanf(read1,"%i",&NumGhostElems_in);  
  assert(fgets(buffer, sizeof buffer, read1)!=NULL);
  garbage=fscanf(read1,"%i",&NumNodes_in);  
  assert(fgets(buffer, sizeof buffer, read1)!=NULL);
  garbage=fscanf(read1,"%i",&NumPhysNodes_in);  
  assert(fgets(buffer, sizeof buffer, read1)!=NULL);
  garbage=fscanf(read1,"%i",&NumBndNodes_in);  
  assert(fgets(buffer, sizeof buffer, read1)!=NULL);
  garbage=fscanf(read1,"%i",&NumEdges_in);  
  assert(fgets(buffer, sizeof buffer, read1)!=NULL);
  garbage=fscanf(read1,"%i",&NumBndEdges_in);  
  assert(fgets(buffer, sizeof buffer, read1)!=NULL);
  int tmpint;
  garbage=fscanf(read1,"%i",&tmpint);
  if (tmpint==1)
    { is_submesh = true; }
  else
    { is_submesh = false; }
  assert(fgets(buffer, sizeof buffer, read1)!=NULL);

  fclose(read1);
  
  // Error check
  if (NumElems_in!=NumElems || 
      NumPhysElems_in!=NumPhysElems ||
      NumGhostElems_in!=NumGhostElems ||
      NumNodes_in!=NumNodes ||
      NumPhysNodes_in!=NumPhysNodes ||
      NumBndNodes_in!=NumBndNodes ||
      NumEdges_in!=NumEdges ||
      NumBndEdges_in!=NumBndEdges)
    {
      printf(" ERROR: trying to read-in mesh information \n");
      printf("        with sizes that do not match ... \n");
      printf("\n");
      printf("      NumElems_in = %8i,        NumElems = %8i\n",NumElems_in,NumElems);
      printf("  NumPhysElems_in = %8i,    NumPhysElems = %8i\n",NumPhysElems_in,NumPhysElems);
      printf(" NumGhostElems_in = %8i,   NumGhostElems = %8i\n",NumGhostElems_in,NumGhostElems);
      printf("      NumNodes_in = %8i,        NumNodes = %8i\n",NumNodes_in,NumNodes);
      printf("  NumPhysNodes_in = %8i,    NumPhysNodes = %8i\n",NumPhysNodes_in,NumPhysNodes);
      printf("   NumBndNodes_in = %8i,     NumBndNodes = %8i\n",NumBndNodes_in,NumBndNodes);
      printf("      NumEdges_in = %8i,        NumEdges = %8i\n",NumEdges_in,NumEdges);
      printf("   NumBndEdges_in = %8i,     NumBndEdges = %8i\n",NumBndEdges_in,NumBndEdges);
      printf("\n");
      exit(1);
    }

  // input: NODE
  snprintf(fname,1024,"%s/mesh_node.dat",inputdir);
  FILE* read2 = fopen(fname,"r");  
  for (int i=1; i<=NumNodes; i++)
    {
      double tmp_double;
      garbage=fscanf(read2,"%lf",&tmp_double);
      node->set(i,1, tmp_double );
      garbage=fscanf(read2,"%lf",&tmp_double);
      node->set(i,2, tmp_double );
    }
  fclose(read2);

  // input: EDGE
  snprintf(fname,1024,"%s/mesh_edge.dat",inputdir);
  FILE* read3 = fopen(fname,"r");
  for (int i=1; i<=NumEdges; i++)
    {
      double tmp_double;
      garbage=fscanf(read3,"%lf",&tmp_double);
      edge->set(i,1, tmp_double );
      
      garbage=fscanf(read3,"%lf",&tmp_double);
      edge->set(i,2, tmp_double );
	
      garbage=fscanf(read3,"%lf",&tmp_double);
      edge->set(i,3, tmp_double );

      garbage=fscanf(read3,"%lf",&tmp_double);
      edge->set(i,4, tmp_double );
    }
  fclose(read3);

  // input: TNODE
  snprintf(fname,1024,"%s/mesh_tnode.dat",inputdir);
  FILE* read4 = fopen(fname,"r");
  for (int i=1; i<=NumElems; i++)
    {
      int tmp_int;
      garbage=fscanf(read4,"%i",&tmp_int);
      tnode->set(i,1, tmp_int );
	
      garbage=fscanf(read4,"%i",&tmp_int);
      tnode->set(i,2, tmp_int );

      garbage=fscanf(read4,"%i",&tmp_int);
      tnode->set(i,3, tmp_int );
    }
  fclose(read4);
    
  // input: TEDGE
  snprintf(fname,1024,"%s/mesh_tedge.dat",inputdir);
  FILE* read5a = fopen(fname,"r");
  for (int i=1; i<=NumElems; i++)
    {
      int tmp_int;
      garbage=fscanf(read5a,"%i",&tmp_int);
      tedge->set(i,1, tmp_int );
      
      garbage=fscanf(read5a,"%i",&tmp_int);
      tedge->set(i,2, tmp_int );
      
      garbage=fscanf(read5a,"%i",&tmp_int);
      tedge->set(i,3, tmp_int );
    }
  fclose(read5a);

  // input: TEDGE_ORIENTATION
  snprintf(fname,1024,"%s/mesh_tedge_orientation.dat",inputdir);
  FILE* read5b = fopen(fname,"r");
  for (int i=1; i<=NumElems; i++)
    {
      int tmp_int;
      garbage=fscanf(read5b,"%i",&tmp_int);
      tedge_orientation->set(i,1, tmp_int );

      garbage=fscanf(read5b,"%i",&tmp_int);
      tedge_orientation->set(i,2, tmp_int );
      
      garbage=fscanf(read5b,"%i",&tmp_int);
      tedge_orientation->set(i,3, tmp_int );
    }
  fclose(read5b);

  // input: EELEM
  snprintf(fname,1024,"%s/mesh_eelem.dat",inputdir);
  FILE* read6 = fopen(fname,"r");
  for (int i=1; i<=NumEdges; i++)
    {
      int tmp_int;
      garbage=fscanf(read6,"%i",&tmp_int);
      eelem->set(i,1, tmp_int );
      
      garbage=fscanf(read6,"%i",&tmp_int);
      eelem->set(i,2, tmp_int );
    }
  fclose(read6);

  // input: ENODE
  snprintf(fname,1024,"%s/mesh_enode.dat",inputdir);
  FILE* read7 = fopen(fname,"r");
  for (int i=1; i<=NumEdges; i++)
    {
      int tmp_int;
      garbage=fscanf(read7,"%i",&tmp_int);
      enode->set(i,1, tmp_int );

      garbage=fscanf(read7,"%i",&tmp_int);
      enode->set(i,2, tmp_int );
    }
  fclose(read7);

  // input: GHOST_LINK
  snprintf(fname,1024,"%s/mesh_ghost_link.dat",inputdir);
  FILE* read8 = fopen(fname,"r");
  for (int i=1; i<=NumGhostElems; i++)
    {
      int tmp_int;
      garbage=fscanf(read8,"%i",&tmp_int);
      ghost_link->set(i, tmp_int );
    }
  fclose(read8);

  // input: EXT_NODE_LINK
  snprintf(fname,1024,"%s/mesh_ext_node_link.dat",inputdir);
  FILE* read8b = fopen(fname,"r");
  for (int i=1; i<=(NumNodes-NumPhysNodes); i++)
    {
      int tmp_int;
      garbage=fscanf(read8b,"%i",&tmp_int);
      ext_node_link->set(i, tmp_int );
    }
  fclose(read8b);

  // input: BND_NODE
  snprintf(fname,1024,"%s/mesh_bnd_node.dat",inputdir);
  FILE* read9 = fopen(fname,"r");
  for (int i=1; i<=NumBndNodes; i++)
    {
      int tmp_int;
      garbage=fscanf(read9,"%i",&tmp_int);
      bnd_node->set(i, tmp_int );
    }
  fclose(read9);

  // input: BND_EDGE
  snprintf(fname,1024,"%s/mesh_bnd_edge.dat",inputdir);
  FILE* read9b = fopen(fname,"r");
  for (int i=1; i<=NumBndEdges; i++)
    {
      int tmp_int;
      garbage=fscanf(read9b,"%i",&tmp_int);
      bnd_edge->set(i, tmp_int );
    }
  fclose(read9b);

  // input: AREA_PRIM
  snprintf(fname,1024,"%s/mesh_area_prim.dat",inputdir);
  FILE* read10 = fopen(fname,"r");
  for (int i=1; i<=NumElems; i++)
    {
      double tmp_double;
      garbage=fscanf(read10,"%lf",&tmp_double);
      area_prim->set(i, tmp_double );
    }
  fclose(read10);
  
  // input: AREA_DUAL
  snprintf(fname,1024,"%s/mesh_area_dual.dat",inputdir);
  FILE* read11 = fopen(fname,"r");
  for (int i=1; i<=NumNodes; i++)
    {
      double tmp_double;
      garbage=fscanf(read11,"%lf",&tmp_double);
      area_dual->set(i, tmp_double );
    }
  fclose(read11);

  // input: JMAT
  snprintf(fname,1024,"%s/mesh_jmat.dat",inputdir);
  FILE* read12 = fopen(fname,"r");
  for (int i=1; i<=NumElems; i++)
    for (int m1=1; m1<=2; m1++)
      for (int m2=1; m2<=2; m2++)
	{
	  double tmp_double;
	  garbage=fscanf(read12,"%lf",&tmp_double);
	  Jmat->set(i,m1,m2, tmp_double );
	}
  fclose(read12);

  // input: NUMELEMSPERNODE
  snprintf(fname,1024,"%s/mesh_numelemspernode.dat",inputdir);
  FILE* read13 = fopen(fname,"r");
  for (int i=1; i<=NumNodes; i++)
    {
      int tmp_int;
      garbage=fscanf(read13,"%i",&tmp_int);
      NumElemsPerNode->set(i, tmp_int );
    }
  fclose(read13);

  // input: ADJACENT
  snprintf(fname,1024,"%s/mesh_adjacent.dat",inputdir);
  FILE* read14 = fopen(fname,"r");
  for (int i=1; i<=NumElems; i++)
    {
      int tmp_int;
      garbage=fscanf(read14,"%i",&tmp_int);
      adjacent->set(i,1, tmp_int );
      garbage=fscanf(read14,"%i",&tmp_int);
      adjacent->set(i,2, tmp_int );
      garbage=fscanf(read14,"%i",&tmp_int);
      adjacent->set(i,3, tmp_int );
    }
  fclose(read14);

  // input: SKEL2ELEM
  snprintf(fname,1024,"%s/mesh_skel2elem.dat",inputdir);
  FILE* read15 = fopen(fname,"r");
  skel2elem_on = false;
  if (read15!=NULL)
    {      
      skel2elem_on = true;
      for (int i=1; i<=NumElems; i++)
	for (int m1=1; m1<=20; m1++)
	  for (int m2=1; m2<=20; m2++)
	    {
	      double tmp_double;
	      garbage=fscanf(read15,"%lf",&tmp_double);
	      Skel2Elem->set(i,m1,m2, tmp_double );
	    }
      fclose(read15);
    }
  
  // input: KMI
  snprintf(fname,1024,"%s/mesh_kmi.dat",inputdir);
  FILE* read16 = fopen(fname,"r");
  kmi_on = false;
  if (read16!=NULL)
    {      
      kmi_on = true;
      for (int i=1; i<=NumEdges; i++)
	for (int m1=1; m1<=3; m1++)
	  for (int m2=1; m2<=6; m2++)
	    {
	      double tmp_double;
	      garbage=fscanf(read16,"%lf",&tmp_double);
	      KMI->set(i,m1,m2, tmp_double );
	    }
      fclose(read16);
    }
  
  // sub-mesh
  if (is_submesh)
    {  InputSubMesh(inputdir);  }

}

const int& mesh::get_NumElems() const
// Returns "NumElems"
{
  return NumElems;
}

const int& mesh::get_NumPhysElems() const
// Returns "NumPhysElems"
{
  return NumPhysElems;
}

const int& mesh::get_NumGhostElems() const
// Returns "NumGhostElems"
{
  return NumGhostElems;
}

const int& mesh::get_NumNodes() const
// Returns "NumNodes"
{
  return NumNodes;
}

const int& mesh::get_NumPhysNodes() const
// Returns "NumPhysNodes"
{
  return NumPhysNodes;
}

const int& mesh::get_NumBndNodes() const
// Returns "NumBndNodes"
{
  return NumBndNodes;
}

const int& mesh::get_NumBndEdges() const
// Returns "NumBndEdges"
{
  return NumBndEdges;
}

const int& mesh::get_NumEdges() const
// Returns "NumEdges"
{
  return NumEdges;
}

const double& mesh::get_node(int i,int m) const
// Returns x (m=1) or y (m=2) coordinate of ith node
{
  return node->get(i,m);
}

const double& mesh::get_edge(int i,int m) const
// Returns (x1,y1) and (x2,y2) values of edge:
//     m = 1  --->  x1
//     m = 2  --->  y1
//     m = 3  --->  x2
//     m = 4  --->  y2
{
  return edge->get(i,m);
}

const int& mesh::get_tnode(int i,int m) const
// Returns pointer to node for ith element (m=1,2, or 3)
{
  return tnode->get(i,m);
}

const int& mesh::get_tedge(int i,int m) const
// Returns pointer to edge for ith element (m=1,2, or 3)
{
  return tedge->get(i,m);
}

const int& mesh::get_tedge_orientation(int i,int m) const
// Returns +1 if normal to edge points out of ith element (m=1,2, or 3) -- outward pointing normal
// Returns -1 if normal to edge points into   ith element (m=1,2, or 3) -- inward pointing normal
{
  return tedge_orientation->get(i,m);
}

const int& mesh::get_adjacent(int i, int m) const
// Returns pointer to neighboring elements for ith element (m=1,2, or 3)
{
  return adjacent->get(i,m);
}

const int& mesh::get_eelem(int i,int m) const
// Returns pointer to element for ith edge (m=1 or 2)
{
  return eelem->get(i,m);
}

const int& mesh::get_enode(int i,int m) const
// Returns pointer to node for ith edge (m=1 or 2)
{
  return enode->get(i,m);
}

const int& mesh::get_ghost_link(int i) const
// Returns pointer to element on opposite side of boundary to ghost cell
{
  return ghost_link->get(i);
}

const int& mesh::get_ext_node_link(int i) const
// Returns pointer to node on opposite side of boundary to external node
{
  return ext_node_link->get(i);
}

const int& mesh::get_bnd_node(int i) const
// Returns pointer to node on boundary
{
  return bnd_node->get(i);
}

const double& mesh::get_area_prim(int i) const
// Returns area of element i
{
  return area_prim->get(i);
}

const double& mesh::get_area_dual(int i) const
// Returns area of dual element centered at node i
{
  return area_dual->get(i);
}

const double& mesh::get_longest_edge() const
// Returns area of dual element centered at node i
{
  return longest_edge_len;
}

const double& mesh::get_shortest_edge() const
// Returns area of dual element centered at node i
{
  return shortest_edge_len;
}

const double& mesh::get_jmat(int i, int m1, int m2) const
// Returns (m1,m2) component of Jacobian matrix in element i
{
  return Jmat->get(i,m1,m2);
}

const int& mesh::get_NumElemsPerNode(int i) const
// Returns the number of elements attached to node i
{
  return NumElemsPerNode->get(i);
}

void mesh::set_node(int i,int m, double input)
// Sets x (m=1) or y (m=2) coordinate of ith node
{
  node->set(i,m, input );
}

void mesh::set_edge(int i,int m, double input)
// Sets (x1,y1) and (x2,y2) values of edge:
//     m = 1  --->  x1
//     m = 2  --->  y1
//     m = 3  --->  x2
//     m = 4  --->  y2
{
  edge->set(i,m, input );
}

void mesh::set_tnode(int i,int m, int input)
// Sets pointer to node for ith element (m=1,2, or 3)
{
  tnode->set(i,m, input );
}

void mesh::set_tedge(int i,int m, int input)
// Sets pointer to edge for ith element (m=1,2, or 3)
{
  tedge->set(i,m, input );
}

void mesh::set_tedge_orientation(int i,int m, int input)
// Set value to +1 if normal to edge points out of ith element (m=1,2, or 3) -- outward pointing normal
// Set value to -1 if normal to edge points into   ith element (m=1,2, or 3) -- inward pointing normal
{
  if (input!=1 && input!=-1 && input!=0)
    {
      printf("\n");
      printf(" ERROR in set_tedge_orientation, input must be set to 0, -1 or +1. \n");
      printf(" input = %i\n",input);
      printf("\n");
      exit(1);
    }
  tedge_orientation->set(i,m, input );
}

void mesh::set_eelem(int i,int m, int input)
// Sets pointer to element for ith edge (m=1 or 2)
{
  eelem->set(i,m, input );
}

void mesh::set_enode(int i,int m, int input)
// Sets pointer to node for ith edge (m=1 or 2)
{
  enode->set(i,m, input );
}

void mesh::set_ghost_link(int i, int input)
// Sets pointer to element on opposite side of boundary to ghost cell
{
  ghost_link->set(i, input );
}

void mesh::set_ext_node_link(int i, int input)
// Sets pointer to node on opposite side of boundary to external node
{
  ext_node_link->set(i, input );
}

void mesh::set_bnd_node(int i, int input)
// Sets pointer to node on boundary
{
  bnd_node->set(i, input );
}

void mesh::set_bnd_edge(int i, int input)
// Sets pointer to edge on boundary
{
  bnd_edge->set(i, input );
}

void mesh::set_area_prim(int i, double input)
// Sets area of element i
{
  area_prim->set(i, input );
}

void mesh::set_area_dual(int i, double input)
// Sets area of dual element centered at node i
{
  area_dual->set(i, input );
}

void mesh::set_longest_edge_len( double len )
{
  longest_edge_len = len;
}

void mesh::set_shortest_edge_len( double len )
{
  shortest_edge_len = len;
}

void mesh::Compute_NumElemsPerNode()
// Computes the number of elements attached to each node
{
  for (int i=1; i<=NumNodes; i++)
    {
      NumElemsPerNode->set(i, 0 );
    }
  
  for (int j=1; j<=NumElems; j++)
    for (int k=1; k<=3; k++)
      {
	int i = tnode->get(j,k);
	int num_old = NumElemsPerNode->get(i);
	NumElemsPerNode->set(i,num_old+1);
      }
}

void mesh::ComputeJacobian()
// Compute 2x2 Jacobian transformation matrix for each element i
// The Jacobian matrix tells one how to transform a gradient in
//   in the canonical variables (xi,eta) to a gradient in 
//   physical variables (x,y):
//
//     phi_x = Jmat(1,1)*phi_xi + Jmat(1,2)*phi_eta
//     phi_y = Jmat(2,1)*phi_xi + Jmat(2,2)*phi_eta
{
  // loop over each element
  for (int i=1; i<=NumElems; i++)
    {
      // Find coordinates of current element
      int i1 = tnode->get(i,1);
      int i2 = tnode->get(i,2);
      int i3 = tnode->get(i,3);

      double x1 = node->get(i1,1);
      double y1 = node->get(i1,2);

      double x2 = node->get(i2,1);
      double y2 = node->get(i2,2);
      
      double x3 = node->get(i3,1);
      double y3 = node->get(i3,2);

      double Area = area_prim->get(i);
      
      Jmat->set(i,1,1, (y3-y1)/(2.0*Area) );
      Jmat->set(i,1,2, (y1-y2)/(2.0*Area) );
      Jmat->set(i,2,1, (x1-x3)/(2.0*Area) );
      Jmat->set(i,2,2, (x2-x1)/(2.0*Area) );
    }
}


void mesh::SetAdjacency()
// Find all elements that are adjacent to current element
// and store in "iTensor2* adjacent"
{
  for (int i=1; i<=NumElems; i++)
    {
      int te,e1,e2;

      te = tedge->get(i,1);
      if (te>0)
        {
	  e1 = eelem->get(te,1);
	  e2 = eelem->get(te,2);
	  if (e1!=i)
            { adjacent->set(i,1, e1 ); }
	  else
            { adjacent->set(i,1, e2 ); }
        }
      else
        {
	  adjacent->set(i,1, -1 );
        }

      te = tedge->get(i,2);
      if (te>0)
        {
	  e1 = eelem->get(te,1);
	  e2 = eelem->get(te,2);
	  if (e1!=i)
            { adjacent->set(i,2, e1 ); }
	  else
            { adjacent->set(i,2, e2 ); }
        }
      else
        {
	  adjacent->set(i,2, -1 );
        }
      
      te = tedge->get(i,3);
      if (te>0)
        {
	  e1 = eelem->get(te,1);
	  e2 = eelem->get(te,2);
	  if (e1!=i)
            { adjacent->set(i,3, e1 ); }
	  else
            { adjacent->set(i,3, e2 ); }
        }
      else
        {
	  adjacent->set(i,3, -1 );
        }
    }
}

// get elements of matrix to convert a vector defined on the
// mesh skeleton to an element-center quantity
const double& mesh::get_Skel2Elem(int i, int m, int n) const
{
  return Skel2Elem->get(i,m,n);
}

// set elements of matrix to convert a vector defined on the
// mesh skeleton to an element-center quantity
void mesh::set_Skel2Elem(int i, int m, int n, double val)
{
  Skel2Elem->set(i,m,n, val );
}

const bool& mesh::is_skel2elem_on() const
{ return skel2elem_on; }

void mesh::turn_skel2elem_on()
{ skel2elem_on = true; }

void mesh::turn_skel2elem_off()
{ skel2elem_on = false; }

// get elements of matrix to convert a vector node information to
// edge-based quantities
const double& mesh::get_KMI(int i, int m, int n) const
{
  return KMI->get(i,m,n);
}

// set elements of matrix to convert a vector node information to
// edge-based quantities
void mesh::set_KMI(int i, int m, int n, double val)
{
  KMI->set(i,m,n, val );
}

const bool& mesh::is_kmi_on() const
{ return kmi_on; }

void mesh::turn_kmi_on()
{ kmi_on = true; }

void mesh::turn_kmi_off()
{ kmi_on = false; }


// Create sub-mesh
void mesh::CreateSubMesh(int num_divide, 
			 double xmin, double xmax,
			 double ymin, double ymax,
			 double deps, double (*SignedDistance)(point))
{
  assert(is_submesh==false);

  if (num_divide>1)
    {
      is_submesh = true;
      SubFactor = num_divide;      
    }
  else
    {
      printf("\n");
      printf(" Error in Mesh.CreateSubMesh: Invalid num_divide input parameter.\n");
      printf("         num_divide = %i\n",num_divide);
      printf("\n");
      exit(1);
    }

  SubNumElems = NumElems*SubFactor*SubFactor;
  SubNumPhysElems = NumPhysElems*SubFactor*SubFactor;
  elem_subs = new iTensor2(NumPhysElems,SubFactor*SubFactor);
  node_subs = new iTensor2(NumPhysElems,((SubFactor+1)*(SubFactor+2))/2);
  sub_tnode = new iTensor2(SubNumPhysElems,3);
  sub_area_prim = new dTensor1(SubNumPhysElems);

  const int tmp_num_sub_nodes = (((SubFactor+1)*(SubFactor+2))/2)*NumPhysElems;
  int node_count = 0;
  int elem_count = 0;
  const double dF = 1.0/double(SubFactor);
  const double onethird = 1.0/3.0;
  dTensor2 tmp_sub_node(tmp_num_sub_nodes,3);
  dTensor1 node_val(tmp_num_sub_nodes);
  iTensor1 index(tmp_num_sub_nodes);

  int n = 0;
  iTensor2 nindex(SubFactor+1,SubFactor+1);
  iTensor2 eindex(2*SubFactor-1,SubFactor);

  for (int k=1; k<=(SubFactor+1); k++)
    for (int j=1; j<=(SubFactor+2-k); j++)
      {
	n = n+1;
	nindex.set(j,k, n );
      }

  int e = 0;
  for (int k=1; k<=SubFactor; k++)
    for (int j=1; j<=(2*(SubFactor-k)+1); j++)
      {
	e = e+1;
	eindex.set(j,k, e );
      }

  for (int i=1; i<=NumPhysElems; i++)
    {
      // basic information
      const int i1 = tnode->get(i,1);
      const int i2 = tnode->get(i,2);
      const int i3 = tnode->get(i,3);

      const double x1 = node->get(i1,1);
      const double y1 = node->get(i1,2);

      const double x2 = node->get(i2,1);
      const double y2 = node->get(i2,2);

      const double x3 = node->get(i3,1);
      const double y3 = node->get(i3,2);

      const double xc = onethird*(x1+x2+x3);
      const double yc = onethird*(y1+y2+y3);

      // nodes            
      for (int k=1; k<=(SubFactor+1); k++)
	for (int j=1; j<=(SubFactor+2-k); j++)
	  {
	    node_count = node_count + 1;
	    const double xi  = -onethird + (j-1)*dF;
	    const double eta = -onethird + (k-1)*dF;

	    const double x = xc + xi*(x2-x1) + eta*(x3-x1);
	    const double y = yc + xi*(y2-y1) + eta*(y3-y1);

	    tmp_sub_node.set(node_count,1, x );
	    tmp_sub_node.set(node_count,2, y );
	    
	    node_val.set(node_count, ((x-xmin)/(xmax-xmin)) + 10000.0*(1.0+(y-ymin)/(ymax-ymin)) );
	    index.set(node_count, node_count );
	  }
    }

  QuickSort_double(node_val,index,1,tmp_num_sub_nodes);

  dTensor1 node_val_unique(tmp_num_sub_nodes);
  iTensor1 index_unique(tmp_num_sub_nodes);
  iTensor1 rev_index(tmp_num_sub_nodes);
  for (int k=1; k<=tmp_num_sub_nodes; k++)
    {
      node_val_unique.set(k, 0.0 );
      index_unique.set(k, k );
      rev_index.set(k, k );
    }
  SubNumPhysNodes = Unique_double(node_val,index,node_val_unique,index_unique,rev_index,
				  1,tmp_num_sub_nodes);
  
  sub_node = new dTensor2(SubNumPhysNodes,2);
  for (int i=1; i<=SubNumPhysNodes; i++)
    {
      sub_node->set(i,1, tmp_sub_node.get(index_unique.get(i),1) );
      sub_node->set(i,2, tmp_sub_node.get(index_unique.get(i),2) );
    }

  // elements
  node_count = 0;
  for (int i=1; i<=NumPhysElems; i++)
    {
      int n = 0;      
      for (int k=1; k<=SubFactor; k++)
        {
	  int jtmp = 0;
	  for (int j=1; j<=(2*(SubFactor-k)+1); j++)
            {
	      elem_count = elem_count + 1;
	      n = n+1;
	      
	      elem_subs->set(i,n, elem_count );
	      
	      if (j%2 == 1)
                {
		  jtmp = jtmp + 1;		  
		  sub_tnode->set(elem_count,1, rev_index.get(node_count + nindex.get(jtmp,k  )) );
		  sub_tnode->set(elem_count,2, rev_index.get(node_count + nindex.get(jtmp+1,k)) );
		  sub_tnode->set(elem_count,3, rev_index.get(node_count + nindex.get(jtmp,k+1)) );
                }
	      else
                {
		  sub_tnode->set(elem_count,1, rev_index.get(node_count + nindex.get(jtmp+1,k))   );
		  sub_tnode->set(elem_count,2, rev_index.get(node_count + nindex.get(jtmp+1,k+1)) );
		  sub_tnode->set(elem_count,3, rev_index.get(node_count + nindex.get(jtmp,k+1))   );
                }
            }	 
        }
      node_count = node_count + ((SubFactor+1)*(SubFactor+2))/2;
    }

  // for each super-element store pointers to nodes on this super-element
  for (int i=1; i<=NumPhysElems; i++)
    {      
      iTensor1 sub_node_list(3*SubFactor*SubFactor);
      iTensor1 index(3*SubFactor*SubFactor);
      int ncount = 1;
      for (int k=1; k<=(SubFactor*SubFactor); k++)
        {  	  
	  int t_tmp = elem_subs->get(i,k);
	  
	  sub_node_list.set(ncount,   sub_tnode->get(t_tmp,1) );
	  sub_node_list.set(ncount+1, sub_tnode->get(t_tmp,2) );
	  sub_node_list.set(ncount+2, sub_tnode->get(t_tmp,3) );
	  
	  index.set(ncount,   ncount   );
	  index.set(ncount+1, ncount+1 );
	  index.set(ncount+2, ncount+2 );
	  
	  ncount = ncount+3;
        }
      iTensor1 aout(3*SubFactor*SubFactor);
      iTensor1 iout(3*SubFactor*SubFactor);
      iTensor1 irout(3*SubFactor*SubFactor);
      QuickSort_int(sub_node_list,index,1,3*SubFactor*SubFactor);
      int num_nodes = Unique_int(sub_node_list,index,aout,iout,irout,1,3*SubFactor*SubFactor);
      
      assert( (((SubFactor+1)*(SubFactor+2))/2) == num_nodes );

      // Set pointer to all sub-nodes that live on a given element
      // First organize to make node numbering easy
      const int i1 = tnode->get(i,1);
      const int i2 = tnode->get(i,2);
      const int i3 = tnode->get(i,3);
      const double x1  = node->get(i1,1);
      const double x2  = node->get(i2,1);
      const double x3  = node->get(i3,1);
      const double y1  = node->get(i1,2);
      const double y2  = node->get(i2,2);
      const double y3  = node->get(i3,2);
      const double xc = onethird*(x1+x2+x3);
      const double yc = onethird*(y1+y2+y3);
      const double T  = area_prim->get(i);
      iTensor1 ns_ind(num_nodes);
      dTensor1 ns_val(num_nodes);
      for (int k=1; k<=num_nodes; k++)
        {
	  int ind = aout.get(k);
	  double xx = sub_node->get(ind,1);
	  double yy = sub_node->get(ind,2);
	  double xi  = ((y3-y1)*(xx-xc)+(x1-x3)*(yy-yc))/(2.0*T);
	  double eta = ((y1-y2)*(xx-xc)+(x2-x1)*(yy-yc))/(2.0*T);

	  ns_ind.set(k, k);
	  ns_val.set(k, ((xi+onethird)) + 10000.0*(1.0+(eta+onethird)) ); 	  
        }
      
      QuickSort_double(ns_val,ns_ind,1,num_nodes);

      for (int k=1; k<=num_nodes; k++)
        {  node_subs->set(i,k, aout.get(ns_ind.get(k)) );  }
    }

  // find all sub-nodes that are on a boundary edge
  iTensor2 project_node(NumBndEdges,SubFactor+1);
  for (int i=1; i<=NumBndEdges; i++)
    {
      int e = bnd_edge->get(i);
      int t1 = eelem->get(e,1);
      int t2 = eelem->get(e,2);
      int n1 = enode->get(e,1);
      int n2 = enode->get(e,2);
      
      double x1 = node->get(n1,1);
      double y1 = node->get(n1,2);
      double x2 = node->get(n2,1);
      double y2 = node->get(n2,2);
      
      int t;
      if (t1<=NumPhysElems)
        {  t = t1;  }
      else
        {  t = t2;  }
      
      iTensor1 sub_node_list(3*SubFactor*SubFactor);
      iTensor1 index(3*SubFactor*SubFactor);
      int ncount = 1;
      for (int k=1; k<=(SubFactor*SubFactor); k++)
        {  	  
	  int t_tmp = elem_subs->get(t,k);

	  sub_node_list.set(ncount,   sub_tnode->get(t_tmp,1) );
	  sub_node_list.set(ncount+1, sub_tnode->get(t_tmp,2) );
	  sub_node_list.set(ncount+2, sub_tnode->get(t_tmp,3) );

	  index.set(ncount,   ncount   );
	  index.set(ncount+1, ncount+1 );
	  index.set(ncount+2, ncount+2 );

	  ncount = ncount+3;
        }
      iTensor1 aout(3*SubFactor*SubFactor);
      iTensor1 iout(3*SubFactor*SubFactor);
      iTensor1 irout(3*SubFactor*SubFactor);
      QuickSort_int(sub_node_list,index,1,3*SubFactor*SubFactor);
      int num_nodes = Unique_int(sub_node_list,index,aout,iout,irout,1,3*SubFactor*SubFactor);

      int num_on_edge = 0;      
      for (int k=1; k<=num_nodes; k++)
        {
	  double x = sub_node->get(aout.get(k),1);
	  double y = sub_node->get(aout.get(k),2);

	  double f  = (y2-y1)*(x-x1)-(x2-x1)*(y-y1);

	  if (fabs(f)<=1.0e-13)
            {
	      num_on_edge = num_on_edge+1;
	      project_node.set(i,num_on_edge, aout.get(k) );
            }
        }
      if (num_on_edge!=(SubFactor+1))
        {
	  printf(" e = %i/%i, num_on_edge = %i,    SubFactor+1 = %i\n",e,NumBndEdges,num_on_edge,SubFactor+1);
        }
      assert(num_on_edge==(SubFactor+1));      
    }
  
  // find all boundary nodes
  int tmp_num_bnd_nodes = 0;
  iTensor1 tmp_bnd_nodes(NumBndEdges*(SubFactor+1));
  iTensor1 bnd_index(NumBndEdges*(SubFactor+1));
  iTensor1 bnd_iout(NumBndEdges*(SubFactor+1));
  iTensor1 bnd_riout(NumBndEdges*(SubFactor+1));
  iTensor1 bnd_aout(NumBndEdges*(SubFactor+1));
  for (int i=1; i<=NumBndEdges; i++)
    for (int k=1; k<=(SubFactor+1); k++)
      {
	tmp_num_bnd_nodes = tmp_num_bnd_nodes+1;
	tmp_bnd_nodes.set(tmp_num_bnd_nodes,project_node.get(i,k));
	bnd_index.set(tmp_num_bnd_nodes, tmp_num_bnd_nodes );
      }
  //  SubNumBndNodes = tmp_num_bnd_nodes;
  QuickSort_int(tmp_bnd_nodes,bnd_index,1,tmp_num_bnd_nodes);
  SubNumBndNodes = Unique_int(tmp_bnd_nodes,bnd_index,
			      bnd_aout,bnd_iout,bnd_riout,
			      1,tmp_num_bnd_nodes);
  sub_bnd_node = new iTensor1(SubNumBndNodes);
  for (int i=1; i<=SubNumBndNodes; i++)
    {
      sub_bnd_node->set(i, bnd_aout.get(i) );
    }

  // compute element areas, force CCW orientation
  for (int i=1; i<=SubNumPhysElems; i++)
    {
      const int i1 = sub_tnode->get(i,1);
      const int i2 = sub_tnode->get(i,2);
      const int i3 = sub_tnode->get(i,3);

      double area = 0.5 * ( sub_node->get(i1,1) * (sub_node->get(i2,2)-sub_node->get(i3,2)) + 
			    sub_node->get(i2,1) * (sub_node->get(i3,2)-sub_node->get(i1,2)) +
			    sub_node->get(i3,1) * (sub_node->get(i1,2)-sub_node->get(i2,2)) );
      if (area < 0.0)
        {
	  double temp = sub_tnode->get(i,2);
	  sub_tnode->set(i,2, sub_tnode->get(i,3) );
	  sub_tnode->set(i,3, temp );	  
        }
      sub_area_prim->set(i, fabs(area) );    
    }
  
}

const bool& mesh::get_is_submesh() const
{
  return is_submesh;
}

const int& mesh::get_SubFactor() const
{
  assert(is_submesh==true);
  return SubFactor;
}

const int& mesh::get_SubNumPhysElems() const
{
  assert(is_submesh==true);
  return SubNumPhysElems;
}

const int& mesh::get_SubNumPhysNodes() const
{
  assert(is_submesh==true);
  return SubNumPhysNodes;
}

const int& mesh::get_SubNumBndNodes() const
{
  assert(is_submesh==true);
  return SubNumBndNodes;
}

const double& mesh::get_sub_node(int i, int j) const
{
  assert(is_submesh==true);
  return sub_node->get(i,j);
}

const int& mesh::get_sub_tnode(int i, int j) const
{
  assert(is_submesh==true);
  return sub_tnode->get(i,j);
}

const int& mesh::get_sub_bnd_node(int i) const
{
  assert(is_submesh==true);
  return sub_bnd_node->get(i);
}

const double& mesh::get_sub_area_prim(int i) const
{
  assert(is_submesh==true);
  return sub_area_prim->get(i);
}

const int& mesh::get_elem_subs(int i, int j) const
{
  assert(is_submesh==true);
  return elem_subs->get(i,j);
}

const int& mesh::get_node_subs(int i, int j) const
{
  assert(is_submesh==true);
  return node_subs->get(i,j);
}

void mesh::OutputSubMesh(char* outputdir)
{
  // output: constants
  char fname[1024];
  snprintf(fname,1024,"%s/submesh_params.dat",outputdir);
  FILE* write1 = fopen(fname,"w");
  fprintf(write1,"%8i : SubFactor       \n",SubFactor);
  fprintf(write1,"%8i : SubNumPhysElems \n",SubNumPhysElems);
  fprintf(write1,"%8i : SubNumPhysNodes \n",SubNumPhysNodes);
  fprintf(write1,"%8i : SubNumBndNodes  \n",SubNumBndNodes);
  fclose(write1);

  // output: SUBNODE
  snprintf(fname,1024,"%s/submesh_node.dat",outputdir);
  FILE* write2 = fopen(fname,"w");
  for (int i=1; i<=SubNumPhysNodes; i++)
    {
      fprintf(write2,"%24.16e  %24.16e\n",
	      sub_node->get(i,1),
	      sub_node->get(i,2));
    }
  fclose(write2);

  // output: SUBTNODE
  snprintf(fname,1024,"%s/submesh_tnode.dat",outputdir);
  FILE* write3 = fopen(fname,"w");
  for (int i=1; i<=SubNumPhysElems; i++)
    {
      fprintf(write3,"%8i  %8i  %8i\n",
	      sub_tnode->get(i,1),
	      sub_tnode->get(i,2),
	      sub_tnode->get(i,3));
    }
  fclose(write3);
  
  // output: SUB_BND_NODE
  snprintf(fname,1024,"%s/submesh_bnd_node.dat",outputdir);
  FILE* write4 = fopen(fname,"w");
  for (int i=1; i<=SubNumBndNodes; i++)
    {
      fprintf(write4,"%8i\n",sub_bnd_node->get(i));
    }
  fclose(write4);
  
  // output: SUB_AREA_PRIM
  snprintf(fname,1024,"%s/submesh_area_prim.dat",outputdir);
  FILE* write5 = fopen(fname,"w");
  for (int i=1; i<=SubNumPhysElems; i++)
    {
      fprintf(write5,"%24.16e\n",sub_area_prim->get(i));
    }
  fclose(write5);

  // output: ELEM_SUBS
  snprintf(fname,1024,"%s/submesh_elem_subs.dat",outputdir);
  FILE* write6 = fopen(fname,"w");
  for (int i=1; i<=NumPhysElems; i++)
    {
      for (int j=1; j<=(SubFactor*SubFactor); j++)
        {
	  fprintf(write6,"%8i",elem_subs->get(i,j));
        }
      fprintf(write6,"\n");
    }
  fclose(write6);
  
  // output: NODE_SUBS
  snprintf(fname,1024,"%s/submesh_node_subs.dat",outputdir);
  FILE* write7 = fopen(fname,"w");
  for (int i=1; i<=NumPhysElems; i++)
    {
      for (int j=1; j<=(((SubFactor+1)*(SubFactor+2))/2); j++)
        {
	  fprintf(write7,"%8i",node_subs->get(i,j));
        }
      fprintf(write7,"\n");
    }
  fclose(write7);
}

void mesh::InputSubMesh(char* inputdir)
{
  // input: constants
  char fname[1024];
  char buffer[1024];
  snprintf(fname,1024,"%s/submesh_params.dat",inputdir);
  FILE* read1 = fopen(fname,"r");  
  int garbage;

  garbage=fscanf(read1,"%i",&SubFactor);
  assert(fgets(buffer, sizeof buffer, read1)!=NULL);
  garbage=fscanf(read1,"%i",&SubNumPhysElems);
  assert(fgets(buffer, sizeof buffer, read1)!=NULL);
  garbage=fscanf(read1,"%i",&SubNumPhysNodes);
  assert(fgets(buffer, sizeof buffer, read1)!=NULL);
  garbage=fscanf(read1,"%i",&SubNumBndNodes); 
  assert(fgets(buffer, sizeof buffer, read1)!=NULL);
  fclose(read1);

  // dimension arrays
  sub_node      = new dTensor2(SubNumPhysNodes,2);
  sub_tnode     = new iTensor2(SubNumPhysElems,3);
  sub_bnd_node  = new iTensor1(SubNumBndNodes);
  sub_area_prim = new dTensor1(SubNumPhysElems);
  elem_subs     = new iTensor2(NumPhysElems,SubFactor*SubFactor);
  node_subs     = new iTensor2(NumPhysElems,((SubFactor+1)*(SubFactor+2))/2);
  
  // input: SUBNODE
  snprintf(fname,1024,"%s/submesh_node.dat",inputdir);
  FILE* read2 = fopen(fname,"r");
  for (int i=1; i<=SubNumPhysNodes; i++)
    {
      double tmp_double;
      garbage=fscanf(read2,"%lf",&tmp_double);
      sub_node->set(i,1, tmp_double );
      garbage=fscanf(read2,"%lf",&tmp_double);
      sub_node->set(i,2, tmp_double );
    }
  fclose(read2);

  // input: SUBTNODE
  snprintf(fname,1024,"%s/submesh_tnode.dat",inputdir);
  FILE* read3 = fopen(fname,"r");
  for (int i=1; i<=SubNumPhysElems; i++)
    {
      int tmp_int;
      garbage=fscanf(read3,"%i",&tmp_int);
      sub_tnode->set(i,1, tmp_int );
      garbage=fscanf(read3,"%i",&tmp_int);
      sub_tnode->set(i,2, tmp_int );
      garbage=fscanf(read3,"%i",&tmp_int);
      sub_tnode->set(i,3, tmp_int );
    }
  fclose(read3);

  // input: SUB_BND_NODE
  snprintf(fname,1024,"%s/submesh_bnd_node.dat",inputdir);
  FILE* read4 = fopen(fname,"r");
  for (int i=1; i<=SubNumBndNodes; i++)
    {
      int tmp_int;
      garbage=fscanf(read4,"%i",&tmp_int);
      sub_bnd_node->set(i, tmp_int );
    }
  fclose(read4);
  
  // input: SUB_AREA_PRIM
  snprintf(fname,1024,"%s/submesh_area_prim.dat",inputdir);
  FILE* read5 = fopen(fname,"r");
  for (int i=1; i<=SubNumPhysElems; i++)
    {
      double tmp_double;
      garbage=fscanf(read5,"%lf",&tmp_double);
      sub_area_prim->set(i, tmp_double );
    }
  fclose(read5);

  // input: ELEM_SUBS
  snprintf(fname,1024,"%s/submesh_elem_subs.dat",inputdir);
  FILE* read6 = fopen(fname,"r");
  for (int i=1; i<=NumPhysElems; i++)
    for (int j=1; j<=(SubFactor*SubFactor); j++)
      {
	int tmp_int;
	garbage=fscanf(read6,"%i",&tmp_int);
	elem_subs->set(i,j, tmp_int );
      }
  fclose(read6);

  // input: NODE_SUBS
  snprintf(fname,1024,"%s/submesh_node_subs.dat",inputdir);
  FILE* read7 = fopen(fname,"r");
  for (int i=1; i<=NumPhysElems; i++)
    for (int j=1; j<=(((SubFactor+1)*(SubFactor+2))/2); j++)
      {
	int tmp_int;
	garbage=fscanf(read7,"%i",&tmp_int);
	node_subs->set(i,j, tmp_int );
      }
  fclose(read7);
}


void mesh::QuickSort_int(iTensor1& a, iTensor1& index, int lo, int hi)
{
  //  lo is the lower index, hi is the upper index
  //  hi the region of array a that is to be sorted
  int i = lo;
  int j = hi;
  int x=a.get( (i+j)/2 );
  int h;
  int itmp;

  //  partition
  while(i<=j) 
    {           
      while (a.get(i)<x) 
        {i++;} 

      while (a.get(j)>x) 
        {j--;}

      if (i<=j)
        {
	  h=a.get(i);
	  a.set(i,a.get(j));
	  a.set(j,h);

	  itmp = index.get(i);
	  index.set(i, index.get(j) );
	  index.set(j, itmp );

	  i++; j--;
        }
    }
  
  //  recursion
  if (lo<j) QuickSort_int(a, index, lo, j);
  if (i<hi) QuickSort_int(a, index, i, hi);
}

int mesh::Unique_int(const iTensor1& a_in, 
		     const iTensor1& index_in,
		     iTensor1& a_out, 
		     iTensor1& index_out,
		     iTensor1& rev_index,
		     int lo, 
		     int hi)
{
  int NumUniqueNodes = 1;  
  
  int curr = a_in.get(lo);
  int curr_index = index_in.get(lo);
  a_out.set(NumUniqueNodes,curr);
  index_out.set(NumUniqueNodes, curr_index );
  rev_index.set(index_in.get(lo), 1 );
  
  for (int k=(lo+1); k<=hi; k++)
    {      
      int next = a_in.get(k);
      
      if (curr!=next)
        {
	  NumUniqueNodes = NumUniqueNodes + 1;
	  curr = next;
	  a_out.set(NumUniqueNodes, curr );
	  curr_index = index_in.get(k);
	  index_out.set(NumUniqueNodes, curr_index );
        }
      rev_index.set(index_in.get(k), NumUniqueNodes );
    }

  return NumUniqueNodes;
}


void mesh::QuickSort_double(dTensor1& a, iTensor1& index, int lo, int hi)
{
  //  lo is the lower index, hi is the upper index
  //  hi the region of array a that is to be sorted
  int i = lo;
  int j = hi;
  double x=a.get( (i+j)/2 );
  double h;
  int itmp;
  
  //  partition
  while(i<=j) 
    {           
      while (a.get(i)<x) 
        {i++;} 
      
      while (a.get(j)>x) 
        {j--;}
      
      if (i<=j)
        {
	  h=a.get(i);
	  a.set(i,a.get(j));
	  a.set(j,h);
	  
	  itmp = index.get(i);
	  index.set(i, index.get(j) );
	  index.set(j, itmp );

	  i++; j--;
        }
    }
  
  //  recursion
  if (lo<j) QuickSort_double(a, index, lo, j);
  if (i<hi) QuickSort_double(a, index, i, hi);
}

int mesh::Unique_double(const dTensor1& a_in, 
			const iTensor1& index_in,
			dTensor1& a_out, 
			iTensor1& index_out,
			iTensor1& rev_index,
			int lo, 
			int hi)
{
  int NumUniqueNodes = 1;

  double curr = a_in.get(lo);
  int curr_index = index_in.get(lo);
  a_out.set(NumUniqueNodes,curr);
  index_out.set(NumUniqueNodes, curr_index );
  rev_index.set(index_in.get(lo), 1 );
  double next;

  for (int k=(lo+1); k<=hi; k++)
    {
      curr = a_in.get(k-1);
      next = a_in.get(k);

      if (fabs(next-curr)>1.0e-10)
        {
	  NumUniqueNodes = NumUniqueNodes + 1;
	  a_out.set(NumUniqueNodes, next );
	  curr_index = index_in.get(k);
	  index_out.set(NumUniqueNodes, curr_index );
        }
      rev_index.set(index_in.get(k), NumUniqueNodes );
    }
  
  return NumUniqueNodes;
}


/*
*	Begin MeshInputData.cpp
*/

// Read-in input data
void MeshInputData(double& radius,
	int& level)
{
	FILE* read = fopen("input.data", "r");
	char buffer[256];

	int garbage;
	garbage = fscanf(read, "%lf", &radius);
	assert(fgets(buffer, sizeof buffer, read) != NULL);
	garbage = fscanf(read, "%i", &level);
	assert(fgets(buffer, sizeof buffer, read) != NULL);
	fclose(read);

	return;
}

/*
*	Begin MeshCreateIcosa.cpp
*/

void MeshCreateIcosa(const double radius, int& numpts, int& numtri,
	point*& p, triangle*& t)
{
	p = new point[12];
	t = new triangle[20];
	numpts = 12;
	numtri = 20;

	// Some constants
	const double GR = (1.0 + sqrt(5.0)) / 2.0;
	double scale = 2.0*radius / sqrt(10.0 + 2.0*sqrt(5.0));

	// Points and elements for a regular icosahedron
	p[0].x = p[1].x = p[2].x = p[3].x = p[4].z = p[5].z = 0;
	p[6].z = p[7].z = p[8].y = p[9].y = p[10].y = p[11].y = 0;
	p[0].z = p[2].z = p[4].y = p[6].y = p[8].x = p[9].x = GR;
	p[1].z = p[3].z = p[5].y = p[7].y = p[10].x = p[11].x = -1 * GR;
	p[0].y = p[1].y = p[4].x = p[5].x = p[8].z = p[10].z = 1;
	p[2].y = p[3].y = p[6].x = p[7].x = p[9].z = p[11].z = -1;

	t[0].n1 = t[1].n1 = t[2].n1 = t[3].n1 = t[4].n1 = 0;
	t[7].n1 = t[8].n1 = t[17].n1 = t[18].n1 = t[19].n1 = 1;
	t[0].n2 = t[4].n2 = t[12].n1 = t[13].n1 = t[14].n1 = 2;
	t[5].n1 = t[6].n1 = t[7].n2 = t[8].n2 = t[9].n1 = 3;
	t[1].n2 = t[2].n2 = t[16].n1 = t[17].n2 = t[18].n2 = 4;
	t[6].n2 = t[13].n2 = t[14].n2 = t[15].n1 = t[5].n2 = 5;
	t[2].n3 = t[3].n2 = t[10].n1 = t[18].n3 = t[19].n2 = 6;
	t[5].n3 = t[9].n2 = t[11].n1 = t[12].n2 = t[13].n3 = 7;
	t[0].n3 = t[1].n3 = t[14].n3 = t[15].n2 = t[16].n2 = 8;
	t[6].n3 = t[7].n3 = t[15].n3 = t[16].n3 = t[17].n3 = 9;
	t[3].n3 = t[4].n3 = t[10].n2 = t[11].n2 = t[12].n3 = 10;
	t[8].n3 = t[9].n3 = t[10].n3 = t[11].n3 = t[19].n3 = 11;

	// Scale points so that we have a sphere of radius "radius"
	for (int i = 0; i < 12; i++)
	{
		p[i].x = scale * p[i].x;
		p[i].y = scale * p[i].y;
		p[i].z = scale * p[i].z;
	}

	// Rotate icosahedron so that there is one point at each of the poles
	double ynew;
	double znew;
	double angle = asin(-radius*p[0].y / (pow(p[0].y, 2) + pow(p[0].z, 2)));

	for (int i = 0; i < 12; i++)
	{
		ynew = p[i].y * cos(angle) + p[i].z * sin(angle);
		znew = -p[i].y * sin(angle) + p[i].z * cos(angle);
		p[i].y = ynew;
		p[i].z = znew;
	}

}

/*
*	Begin MeshSubDivide.cpp
*/

void MeshSubdivide(const int level,
	const double radius,
	int& numpts,
	int& numtri,
	point*& p,
	triangle*& t)
{
	//  void QuickSort_bar(bar* b, int* key, int i, int j);

	for (int ell = 1; ell <= level; ell++)
	{
		bar* bars;
		bar* bartemp;
		int* priority;
		bars = new bar[3 * numtri];
		bartemp = new bar[3 * numtri];
		priority = new int[3 * numtri];

		int numbar = 0;
		int temp;

		// find unique listing of bars
		for (int i = 0; i < numtri; i++)
		{
			bartemp[3 * i].n1 = t[i].n1;
			bartemp[3 * i].n2 = t[i].n2;
			bartemp[3 * i + 1].n1 = t[i].n1;
			bartemp[3 * i + 1].n2 = t[i].n3;
			bartemp[3 * i + 2].n1 = t[i].n2;
			bartemp[3 * i + 2].n2 = t[i].n3;
		}

		for (int i = 0; i < 3 * numtri; i++)
		{
			// sort each bar so that n1 < n2
			if (bartemp[i].n1 > bartemp[i].n2)
			{
				temp = bartemp[i].n1;
				bartemp[i].n1 = bartemp[i].n2;
				bartemp[i].n2 = temp;
			}

			// assign priority for sorting list of bars
			priority[i] = numpts*bartemp[i].n1 + bartemp[i].n2;
		}

		// sort the bars according to priority      
		QuickSort_bar(bartemp, priority, 0, 3 * numtri - 1);

		// keep only unique bars
		int value = 0;
		for (int i = 0; i < 3 * numtri; i++)
			if (priority[i] != value)
			{
				bars[numbar].n1 = bartemp[i].n1;
				bars[numbar].n2 = bartemp[i].n2;
				value = priority[i];
				numbar++;
			}

		// store all old points and newly created points
		// (located at bar midpoints) in pnew
		point* pnew;
		pnew = new point[numpts + numbar];
		for (int i = 0; i < numpts; i++)
		{
			pnew[i].x = p[i].x;
			pnew[i].y = p[i].y;
			pnew[i].z = p[i].z;
		}

		for (int i = 0; i < numbar; i++)
		{
			pnew[i + numpts].x = .5 * (p[bars[i].n1].x + p[bars[i].n2].x);
			pnew[i + numpts].y = .5 * (p[bars[i].n1].y + p[bars[i].n2].y);
			pnew[i + numpts].z = .5 * (p[bars[i].n1].z + p[bars[i].n2].z);
			bars[i].bt = i + numpts;
		}
		numpts = numpts + numbar;

		// create new triangle list
		triangle* tnew;
		tnew = new triangle[4 * numtri];

		for (int i = 0; i < numtri; i++)
		{
			int mstop = 0;
			int j = 1;
			int ind1, ind2, ind3;
			int bin1, bin2, tin1, tin2, tin3;

			tin1 = t[i].n1;
			tin2 = t[i].n2;
			tin3 = t[i].n3;

			while (j <= numbar && mstop < 3)
			{
				bin1 = bars[numbar - j].n1;
				bin2 = bars[numbar - j].n2;

				if ((tin1 == bin1) && (tin2 == bin2))
				{
					ind1 = numbar - j;
					mstop = mstop + 1;
				}
				else if ((tin1 == bin1) && (tin3 == bin2))
				{
					ind2 = numbar - j;
					mstop = mstop + 1;
				}
				else if ((tin2 == bin1) && (tin3 == bin2))
				{
					ind3 = numbar - j;
					mstop = mstop + 1;
				}
				j = j + 1;
			}
			tnew[4 * i + 0].n1 = t[i].n1;
			tnew[4 * i + 0].n2 = bars[ind1].bt;
			tnew[4 * i + 0].n3 = bars[ind2].bt;
			tnew[4 * i + 1].n1 = t[i].n2;
			tnew[4 * i + 1].n2 = bars[ind1].bt;
			tnew[4 * i + 1].n3 = bars[ind3].bt;
			tnew[4 * i + 2].n1 = t[i].n3;
			tnew[4 * i + 2].n2 = bars[ind2].bt;
			tnew[4 * i + 2].n3 = bars[ind3].bt;
			tnew[4 * i + 3].n1 = bars[ind1].bt;
			tnew[4 * i + 3].n2 = bars[ind2].bt;
			tnew[4 * i + 3].n3 = bars[ind3].bt;
		}

		delete[]t;
		delete[]p;

		numtri = 4 * numtri;

		t = tnew;
		p = pnew;

		// project new points to sphere surface
		for (int i = 0; i < numpts; i++)
		{
			// Get (x,y,z) Cartesian coordinates
			double x = p[i].x;
			double y = p[i].y;
			double z = p[i].z;

			// Compute spherical coordinates
			double   rad = sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2));
			double theta = acos(z / rad);
			double   phi = atan2(y, x);

			// Project point onto a sphere of radius "radius"
			p[i].x = radius * sin(theta) * cos(phi);
			p[i].y = radius * sin(theta) * sin(phi);
			p[i].z = radius * cos(theta);
		}

		// Report that current level is finished
		printf("  Finished level:  %5i\n", ell);

		delete bars;
		delete bartemp;
		delete priority;
	}

}

/*
*	Begin MeshSphereCoords.cpp
*/

// Store spherical coordinates of p in psphere
void MeshSphereCoords(const int numpts, const point* p, point*& psphere)
{

	for (int i = 0; i < numpts; i++)
	{
		double x = p[i].x;
		double y = p[i].y;
		double z = p[i].z;

		double   rad = sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2));
		double theta = acos(z / rad);
		double   phi = atan2(y, x);

		// Longitude
		psphere[i].x = phi - pi;

		// Latitude
		psphere[i].y = 0.5*pi - theta;
	}

}

/*
*	Begin MeshEdgeData.cpp
*/

//
// Create information about edges
//
void MeshEdgeData(const int& numpts,
	const int& numtri,
	const point p[],
	const point psphere[],
	const triangle t[],
	const double area[],
	int& numedges,
	dTensor2*& edge,
	iTensor2*& tedge,
	iTensor2*& eelem)
{
	// Needed sorting functions
	//  void QuickSort(double*& a, int*& index, int lo, int hi);
	//  void Unique(double*& vec, int*& index, int*& index_reverse, int& size);

	// Create helper variables and arrays
	dTensor2 edge_tmp(3 * numtri, 6);
	double* edge_label;
	edge_label = new double[3 * numtri];

	// Loop over each element
	int k = 0;
	for (int i = 0; i < numtri; i++)
	{
		double xmid, ymid, zmid;
		int tp1 = t[i].n1;
		int tp2 = t[i].n2;
		int tp3 = t[i].n3;

		// look at first edge
		k = k + 1;
		edge_tmp.set(k, 1, p[tp1].x);
		edge_tmp.set(k, 2, p[tp1].y);
		edge_tmp.set(k, 3, p[tp1].z);
		edge_tmp.set(k, 4, p[tp2].x);
		edge_tmp.set(k, 5, p[tp2].y);
		edge_tmp.set(k, 6, p[tp2].z);

		xmid = 0.5*(p[tp1].x + p[tp2].x);
		ymid = 0.5*(p[tp1].y + p[tp2].y);
		zmid = 0.5*(p[tp1].z + p[tp2].z);
		edge_label[k - 1] = xmid + 1e4*ymid + 1e8*zmid;

		// look at second edge
		k = k + 1;
		edge_tmp.set(k, 1, p[tp1].x);
		edge_tmp.set(k, 2, p[tp1].y);
		edge_tmp.set(k, 3, p[tp1].z);
		edge_tmp.set(k, 4, p[tp3].x);
		edge_tmp.set(k, 5, p[tp3].y);
		edge_tmp.set(k, 6, p[tp3].z);

		xmid = 0.5*(p[tp1].x + p[tp3].x);
		ymid = 0.5*(p[tp1].y + p[tp3].y);
		zmid = 0.5*(p[tp1].z + p[tp3].z);
		edge_label[k - 1] = xmid + 1e4*ymid + 1e8*zmid;

		// look at third edge
		k = k + 1;
		edge_tmp.set(k, 1, p[tp2].x);
		edge_tmp.set(k, 2, p[tp2].y);
		edge_tmp.set(k, 3, p[tp2].z);
		edge_tmp.set(k, 4, p[tp3].x);
		edge_tmp.set(k, 5, p[tp3].y);
		edge_tmp.set(k, 6, p[tp3].z);

		xmid = 0.5*(p[tp2].x + p[tp3].x);
		ymid = 0.5*(p[tp2].y + p[tp3].y);
		zmid = 0.5*(p[tp2].z + p[tp3].z);
		edge_label[k - 1] = xmid + 1e4*ymid + 1e8*zmid;
	}

	// Keep only unique edges
	int numedges_tmp = 3 * numtri;
	int* index;
	index = new int[3 * numtri];
	int* index_reverse;
	for (int i = 0; i < (3 * numtri); i++)
	{
		index[i] = i;
	}

	QuickSort(edge_label, index, 0, 3 * numtri - 1);
	Unique(edge_label, index, index_reverse, numedges_tmp);

	// Store information up to this point
	numedges = numedges_tmp;
	edge = new dTensor2(numedges, 4);
	tedge = new iTensor2(numtri, 3);
	eelem = new iTensor2(numedges, 2);

	// edge coordinates
	for (int i = 1; i <= numedges; i++)
	{
		int mu = index[i - 1] + 1;

		double   x1 = edge_tmp.get(mu, 1);
		double   y1 = edge_tmp.get(mu, 2);
		double   z1 = edge_tmp.get(mu, 3);
		double lon1 = atan2(y1, x1) - pi;
		double lat1 = 0.5*pi - acos(z1 / sqrt(x1*x1 + y1*y1 + z1*z1));

		double   x2 = edge_tmp.get(mu, 4);
		double   y2 = edge_tmp.get(mu, 5);
		double   z2 = edge_tmp.get(mu, 6);
		double lon2 = atan2(y2, x2) - pi;
		double lat2 = 0.5*pi - acos(z2 / sqrt(x2*x2 + y2*y2 + z2*z2));

		edge->set(i, 1, lon1);
		edge->set(i, 2, lat1);
		edge->set(i, 3, lon2);
		edge->set(i, 4, lat2);
	}

	// edges attached to element i (physical elements)
	for (int i = 1; i <= numtri; i++)
	{
		tedge->set(i, 1, index_reverse[3 * i - 3] + 1);
		tedge->set(i, 2, index_reverse[3 * i - 2] + 1);
		tedge->set(i, 3, index_reverse[3 * i - 1] + 1);
	}

	// find the elements on either side of each edge
	for (int i = 1; i <= numedges; i++)
	{
		int mfound = 0;
		int k = 1;
		while (mfound < 2 && k <= numtri)
		{
			for (int n = 1; n <= 3; n++)
			{
				int j = tedge->get(k, n);

				if (i == j)
				{
					mfound = mfound + 1;
					eelem->set(i, mfound, k);
				}
			}
			k = k + 1;
		}

		if (mfound != 2)
		{
			printf(" ERROR in MeshEdgeData.cpp:  mfound!=2. Should never get here ... \n");
			printf("      i = %8i\n", i);
			printf("      k = %8i\n", k);
			printf(" mfound = %8i\n", mfound);
			printf("\n");
			exit(1);
		}

	}

	delete edge_label;
	delete index_reverse;
	delete index;

}


/*
*	Begin MeshOrientEdge.cpp
*/

//
// Adjust edge information so that the element to the "left"
// of the edge is always "upwind" of the unit normal to the
// edge and the element to the "right" of the edge is always
// "downwind" of the unit normal to the edge
//
void MeshOrientEdge(mesh& Mesh)
{
	int count1 = 0;
	int count2 = 0;

	// Loop over all Edges
	for (int i = 1; i <= Mesh.get_NumEdges(); i++)
	{
		// Edge coordinates
		double x1 = Mesh.get_edge(i, 1);
		double y1 = Mesh.get_edge(i, 2);
		double x2 = Mesh.get_edge(i, 3);
		double y2 = Mesh.get_edge(i, 4);

		// Check periodicity condition
		double LX = x1 - x2;
		double LY = y2 - y1;
		if (LX >= pi)
		{
			LX = LX - 2 * pi;
		}
		else if (LX <= -pi)
		{
			LX = LX + 2 * pi;
		}

		if (LY >= pi)
		{
			LY = LY - 2 * pi;
		}
		else if (LY <= -pi)
		{
			LY = LY + 2 * pi;
		}

		if (fabs(y1 - 0.5*pi) <= 1.0e-12 || fabs(y1 + 0.5*pi) <= 1.0e-12 ||
			fabs(y2 - 0.5*pi) <= 1.0e-12 || fabs(y2 + 0.5*pi) <= 1.0e-12)
		{
			LX = 0.0e0;
		}

		// Edge length
		double L = sqrt(pow(LX, 2) + pow(LY, 2));

		// Unit normal to edge
		dTensor1 nhat(2);
		nhat.set(1, LY / L);
		nhat.set(2, LX / L);

		// Nodes on either side of edge
		int ileft = Mesh.get_eelem(i, 1);
		int iright = Mesh.get_eelem(i, 2);

		double xc_l = (Mesh.get_node(Mesh.get_tnode(ileft, 1), 1)
			+ Mesh.get_node(Mesh.get_tnode(ileft, 2), 1)
			+ Mesh.get_node(Mesh.get_tnode(ileft, 3), 1)) / 3.0;

		double yc_l = (Mesh.get_node(Mesh.get_tnode(ileft, 1), 2)
			+ Mesh.get_node(Mesh.get_tnode(ileft, 2), 2)
			+ Mesh.get_node(Mesh.get_tnode(ileft, 3), 2)) / 3.0;

		double xc_r = (Mesh.get_node(Mesh.get_tnode(iright, 1), 1)
			+ Mesh.get_node(Mesh.get_tnode(iright, 2), 1)
			+ Mesh.get_node(Mesh.get_tnode(iright, 3), 1)) / 3.0;

		double yc_r = (Mesh.get_node(Mesh.get_tnode(iright, 1), 2)
			+ Mesh.get_node(Mesh.get_tnode(iright, 2), 2)
			+ Mesh.get_node(Mesh.get_tnode(iright, 3), 2)) / 3.0;

		dTensor1 rhat(2);
		double RX = xc_r - xc_l;
		double RY = yc_r - yc_l;

		if (RX >= pi)
		{
			RX = RX - 2 * pi;
		}
		else if (RX <= -pi)
		{
			RX = RX + 2 * pi;
		}

		if (RY >= pi)
		{
			RY = RY - 2 * pi;
		}
		else if (RY <= -pi)
		{
			RY = RY + 2 * pi;
		}

		if (fabs(yc_r - 0.5*pi) <= 1.0e-12 || fabs(yc_r + 0.5*pi) <= 1.0e-12 ||
			fabs(yc_l - 0.5*pi) <= 1.0e-12 || fabs(yc_l + 0.5*pi) <= 1.0e-12)
		{
			RX = 0.0e0;
		}

		rhat.set(1, RX);
		rhat.set(2, RY);

		double ddot = nhat.get(1)*rhat.get(1) + nhat.get(2)*rhat.get(2);

		if (ddot < 0.0) // Then need to interchange nodes
		{
			count2 = count2 + 1;

			// Switch element order
			Mesh.set_eelem(i, 1, iright);
			Mesh.set_eelem(i, 2, ileft);
		}
		else
		{
			count1 = count1 + 1;
		}
	}

	//cout << endl;
	//cout << " NumEdges = " << Mesh.get_NumEdges() << endl;
	//cout << "   count1 = " << count1 << endl;
	//cout << "   count2 = " << count2 << endl;
	//cout << endl;

}

/*
*	Begin MeshComputeAreas.cpp
*/

void MeshComputeAreas(const int numtri, const int numpts,
	const point p[], triangle t[],
	double*& area, double*& cdual)
{
	// **************************************************************************
	// ***************** get triangle areas and dual mesh areas *****************
	// **************************************************************************
	for (int i = 0; i < numpts; i++)
	{
		cdual[i] = 0;
	}

	double A, B, C;
	double xt, yt, zt, dot_prod;
	int tmp_switch;

	for (int i = 0; i < numtri; i++)
	{
		A = p[t[i].n1].z*(p[t[i].n2].y - p[t[i].n3].y) + p[t[i].n2].z*(p[t[i].n3].y - p[t[i].n1].y)
			+ p[t[i].n3].z*(p[t[i].n1].y - p[t[i].n2].y);

		B = p[t[i].n1].x*(p[t[i].n2].z - p[t[i].n3].z) + p[t[i].n2].x*(p[t[i].n3].z - p[t[i].n1].z)
			+ p[t[i].n3].x*(p[t[i].n1].z - p[t[i].n2].z);

		C = p[t[i].n1].y*(p[t[i].n2].x - p[t[i].n3].x) + p[t[i].n2].y*(p[t[i].n3].x - p[t[i].n1].x)
			+ p[t[i].n3].y*(p[t[i].n1].x - p[t[i].n2].x);

		area[i] = 0.5 * sqrt(A*A + B*B + C*C);
		cdual[t[i].n1] += (area[i] / 3.0);
		cdual[t[i].n2] += (area[i] / 3.0);
		cdual[t[i].n3] += (area[i] / 3.0);

		xt = (p[t[i].n1].x + p[t[i].n2].x + p[t[i].n3].x) / 3.0;
		yt = (p[t[i].n1].y + p[t[i].n2].y + p[t[i].n3].y) / 3.0;
		zt = (p[t[i].n1].z + p[t[i].n2].z + p[t[i].n3].z) / 3.0;

		dot_prod = xt*A + yt*B + zt*C;

		if (dot_prod > 0.0)
		{
			tmp_switch = t[i].n2;
			t[i].n2 = t[i].n3;
			t[i].n3 = tmp_switch;
		}
	}

}

/*
*	Begin MeshStore.cpp
*/

void MeshStore(point p[],
	triangle t[],
	int bnd_node[],
	int ghost_link[],
	double area[],
	double cdual[],
	dTensor2*& edge,
	iTensor2*& tedge,
	iTensor2*& eelem,
	mesh& Mesh)
{
	int NumElems = Mesh.get_NumElems();
	int NumGhostElems = Mesh.get_NumGhostElems();
	int NumNodes = Mesh.get_NumNodes();
	int NumBndNodes = Mesh.get_NumBndNodes();
	int NumEdges = Mesh.get_NumEdges();

	// node
	for (int i = 1; i <= NumNodes; i++)
	{
		Mesh.set_node(i, 1, p[i - 1].x);
		Mesh.set_node(i, 2, p[i - 1].y);
	}

	// tnode
	for (int i = 1; i <= NumElems; i++)
	{
		Mesh.set_tnode(i, 1, t[i - 1].n1 + 1);
		Mesh.set_tnode(i, 2, t[i - 1].n2 + 1);
		Mesh.set_tnode(i, 3, t[i - 1].n3 + 1);
	}

	// bnd_node
	for (int i = 1; i <= NumBndNodes; i++)
	{
		Mesh.set_bnd_node(i, bnd_node[i - 1] + 1);
	}

	// ghost_link
	for (int i = 1; i <= NumGhostElems; i++)
	{
		Mesh.set_ghost_link(i, ghost_link[i - 1] + 1);
	}


	// area_prim
	for (int i = 1; i <= NumElems; i++)
	{
		Mesh.set_area_prim(i, area[i - 1]);
	}

	// area_dual
	for (int i = 1; i <= NumNodes; i++)
	{
		Mesh.set_area_dual(i, cdual[i - 1]);
	}

	// edge
	for (int i = 1; i <= NumEdges; i++)
		for (int k = 1; k <= 4; k++)
		{
			Mesh.set_edge(i, k, edge->get(i, k));
		}

	// tedge
	for (int i = 1; i <= NumElems; i++)
		for (int k = 1; k <= 3; k++)
		{
			Mesh.set_tedge(i, k, tedge->get(i, k));
		}

	// eelem
	for (int i = 1; i <= NumEdges; i++)
		for (int k = 1; k <= 2; k++)
		{
			Mesh.set_eelem(i, k, eelem->get(i, k));
		}
}

/*
 *	Begin QuickSort_bar.cpp
 */

// Function to sort a list of "bars"
void QuickSort_bar(bar* b, int* key, int i, int j)
{
	// Swap function
	//  void Swap_bar(bar* b, int* key, int i, int j);

	// Partition function
	//  int Partition_bar(bar* b, int* key, int l, int r, int pivot);

	// Quicksort
	int pivotindex = (i + j) / 2;
	Swap_bar(b, key, j, pivotindex);

	int k = Partition_bar(b, key, i - 1, j, key[j]);
	Swap_bar(b, key, j, k);

	if ((k - i) > 1) QuickSort_bar(b, key, i, k - 1);
	if ((j - k) > 1) QuickSort_bar(b, key, k + 1, j);

	return;
}

// Function to swap bars i and j
void Swap_bar(bar* b, int* key, int i, int j)
{
	int temp;
	temp = key[i];
	key[i] = key[j];
	key[j] = temp;
	temp = b[i].n1;
	b[i].n1 = b[j].n1;
	b[j].n1 = temp;
	temp = b[i].n2;
	b[i].n2 = b[j].n2;
	b[j].n2 = temp;

	return;
}

// Function to partition bar according to key
int Partition_bar(bar* b, int* key, int l, int r, int pivot)
{
	void Swap_bar(bar*, int*, int, int);

	do
	{
		while (key[++l] < pivot);
		while (r && (key[--r]) > pivot);
		Swap_bar(b, key, l, r);
	} while (l < r);
	Swap_bar(b, key, l, r);
	return(l);
}

/*
 *	Begin QuickSort.cpp
 */

void QuickSort(double*& a, int*& index, int lo, int hi)
{
	//  lo is the lower index, hi is the upper index
	//  hi the region of array a that is to be sorted

	int i = lo;
	int j = hi;
	double x = a[(i + j) / 2];
	double h;
	int itmp;

	//  partition
	while (i <= j)
	{
		while (a[i] < x)
		{
			i++;
		}

		while (a[j] > x)
		{
			j--;
		}

		if (i <= j)
		{
			h = a[i];
			a[i] = a[j];
			a[j] = h;

			itmp = index[i];
			index[i] = index[j];
			index[j] = itmp;

			i++; j--;
		}
	}

	//  recursion
	if (lo < j) QuickSort(a, index, lo, j);
	if (i < hi) QuickSort(a, index, i, hi);
}