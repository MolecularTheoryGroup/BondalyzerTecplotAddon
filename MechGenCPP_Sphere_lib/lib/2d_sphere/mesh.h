#ifndef MESH_H
#define MESH_H

#include <assert.h>

#include "../tensors.h"
#include "../structs.h"

// 2D Mesh Object ----------------------------------------
class mesh
{
 public:
  mesh(int inNumElems, 
       int inNumPhysElems, 
       int inNumNodes, 
       int inNumPhysNodes, 
       int inNumBndNodes, 
       int inNumEdges,
       int inNumBndEdges);
  // Constructor
  // POST: Creates a mesh

  mesh(const mesh& another_mesh);
  // Copy constructor
  // POST: New mesh created with size and contents same as another_mesh
  
  ~mesh();
  // Destructor
  // POST: mesh no longer exists

  void OutputMesh(char* outputdir);
  // Ouput all mesh information

  void InputMesh(char* outputdir);
  // Input all mesh information

  const int& get_NumElems() const;
  // Returns "NumElems"

  const int& get_NumPhysElems() const;
  // Returns "NumPhysElems"

  const int& get_NumGhostElems() const;
  // Returns "NumGhostElems"

  const int& get_NumNodes() const;
  // Returns "NumNodes" 

  const int& get_NumPhysNodes() const;
  // Returns "NumPhysNodes

  const int& get_NumBndNodes() const;
  // Returns "NumBndNodes"

  const int& get_NumBndEdges() const;
  // Returns "NumBndEdges"

  const int& get_NumEdges() const;
  // Returns "NumEdges"

  const double& get_node(int i,int m) const;
  // Returns x (m=1) or y (m=2) coordinate of ith node
  
  const double& get_edge(int i,int m) const;
  // Returns (x1,y1) and (x2,y2) values of edge:
  //     m = 1  --->  x1
  //     m = 2  --->  y1
  //     m = 3  --->  x2
  //     m = 4  --->  y2

  const int& get_tnode(int i,int m) const;
  // Returns pointer to node for ith element (m=1,2, or 3)

  const int& get_tedge(int i,int m) const;
  // Returns pointer to edge for ith element (m=1,2, or 3)
  
  const int& get_adjacent(int i, int m) const;
  // Returns pointer to neighboring elements for ith element (m=1,2, or 3)
  
  const int& get_tedge_orientation(int i,int m) const;
  // Returns +1 if normal to edge points out of ith element (m=1,2, or 3) -- outward pointing normal
  // Returns -1 if normal to edge points into   ith element (m=1,2, or 3) -- inward pointing normal

  const int& get_eelem(int i,int m) const;
  // Returns pointer to element for ith edge (m=1 or 2)

  const int& get_enode(int i,int m) const;
  // Returns pointer to node for ith edge (m=1 or 2)

  const int& get_ghost_link(int i) const;
  // Returns pointer to element on opposite side of boundary to ghost cell

  const int& get_ext_node_link(int i) const;
  // Returns pointer to node on opposite side of boundary to external node
  
  const int& get_bnd_node(int i) const;
  // Returns pointer to node on boundary

  const double& get_area_prim(int i) const;
  // Returns area of element i

  const double& get_longest_edge() const;
  // returns the longest edge length

  const double& get_shortest_edge() const;
  // returns the shortest edge length

  const double& get_area_dual(int i) const;
  // Returns area of dual element centered at node i

  const double& get_jmat(int i, int m1, int m2) const;
  // Returns (m1,m2) component of Jacobian matrix in element i
  
  const int& get_NumElemsPerNode(int i) const;
  // Returns the number of elements attached to node i
  
  void set_node(int i,int m, double input);
  // Sets x (m=1) or y (m=2) coordinate of ith node
  
  void set_edge(int i,int m, double input);
  // Sets (x1,y1) and (x2,y2) values of edge:
  //     m = 1  --->  x1
  //     m = 2  --->  y1
  //     m = 3  --->  x2
  //     m = 4  --->  y2

  void set_tnode(int i,int m, int input);
  // Sets pointer to node for ith element (m=1,2, or 3)

  void set_tedge(int i,int m, int input);
  // Sets pointer to edge for ith element (m=1,2, or 3)
  
  void set_tedge_orientation(int i,int m, int input);
  // Set value to +1 if normal to edge points out of ith element (m=1,2, or 3) -- outward pointing normal
  // Set value to -1 if normal to edge points into   ith element (m=1,2, or 3) -- inward pointing normal
  
  void set_eelem(int i,int m, int input);
  // Sets pointer to element for ith edge (m=1 or 2)

  void set_enode(int i,int m, int input);
  // Sets pointer to node for ith edge (m=1 or 2)
  
  void set_ghost_link(int i, int input);
  // Sets pointer to element on opposite side of boundary to ghost cell
  
  void set_ext_node_link(int i, int input);
  // Sets pointer to node on opposite side of boundary to external node
  
  void set_bnd_node(int i, int input);
  // Sets pointer to node on boundary

  void set_bnd_edge(int i, int input);
  // Sets pointer to edge on boundary
  
  void set_area_prim(int i, double input);
  // Sets area of element i
  
  void set_longest_edge_len(double len);
  // Sets the longest edge length

  void set_shortest_edge_len(double len);
  // Sets the shortest edge length

  void set_area_dual(int i, double input);
  // Sets area of dual element centered at node i

  void ComputeJacobian();
  // Compute 2x2 Jacobian transformation matrix for each element i
  // The Jacobian matrix tells one how to transform a gradient in
  //   in the canonical variables (xi,eta) to a gradient in 
  //   physical variables (x,y):
  //
  //     phi_x = Jmat(1,1)*phi_xi + Jmat(1,2)*phi_eta
  //     phi_y = Jmat(2,1)*phi_xi + Jmat(2,2)*phi_eta
  
  void Compute_NumElemsPerNode();
  // Computes the number of elements attached to each node

  void SetAdjacency();
  // Find all elements that are adjacent to current element
  // and store in "iTensor2* adjacent"

  const bool& is_skel2elem_on() const;
  void turn_skel2elem_on();
  void turn_skel2elem_off();
  const double& get_Skel2Elem(int i, int m, int n) const;
  // get elements of matrix to convert a vector defined on the
  // mesh skeleton to an element-center quantity
  
  void set_Skel2Elem(int i, int m, int n, double val);
  // set elements of matrix to convert a vector defined on the
  // mesh skeleton to an element-center quantity
  
  const bool& is_kmi_on() const;
  void turn_kmi_on();
  void turn_kmi_off();
  const double& get_KMI(int i, int m, int n) const;
  // get elements of matrix to convert a vector node information to
  // edge-based quantities
  
  void set_KMI(int i, int m, int n, double val);
  // set elements of matrix to convert a vector node information to
  // edge-based quantities

  // ----------------------------------------------------------
  // Sub-mesh information
  // ----------------------------------------------------------
  void CreateSubMesh(int num_divide, 
		     double xmin, double xmax,
		     double ymin, double ymax,
		     double deps, double (*SignedDistance)(point));
  const bool& get_is_submesh() const;
  const int& get_SubFactor() const;
  const int& get_SubNumPhysElems() const;
  const int& get_SubNumPhysNodes() const;
  const int& get_SubNumBndNodes() const;
  const double& get_sub_node(int i, int j) const;
  const int& get_sub_tnode(int i, int j) const;
  const int& get_sub_bnd_node(int i) const;
  const double& get_sub_area_prim(int i) const;
  const int& get_elem_subs(int i, int j) const;
  const int& get_node_subs(int i, int j) const;
  void OutputSubMesh(char* outputdir);
  void InputSubMesh(char* inputdir);
  // ----------------------------------------------------------
  
 private:
  // Basic constants (NOTE:  NumElems = NumPhysElems + NumGhostElems)
  int NumElems;         // Number of total elements in mesh
  int NumPhysElems;     // Number of physical elements in mesh
  int NumGhostElems;    // Number of ghost elements in mesh
  int NumNodes;         // Number of nodes in mesh
  int NumPhysNodes;     // Number of physical nodes in mesh
  int NumEdges;         // Number of edges in mesh 
  int NumBndEdges;      // Number of edges on the boundary
  int NumBndNodes;      // Number of nodes on boundary
  
  double longest_edge_len;
  double shortest_edge_len;
  
  // node: list of x & y coordinates of nodes (size = NumNodes-by-2)
  dTensor2* node; //(NumNodes,2);

  // tnode: list of nodes attached to element (size = NumNodes-by-3)
  iTensor2* tnode; //(NumElems,3);
  
  // list of elements on opposite side of boundary to ghost cell (size = NumGhostElems)
  iTensor1* ghost_link; //(NumGhostElems);
  
  // list of physical nodes that correspond to given external node
  iTensor1* ext_node_link; //(NumNodes-NumPhysNodes)
  
  // list of nodes that lie on boundary
  iTensor1* bnd_node; //(NumBndNodes)

  // list of edges that lie on boundary
  iTensor1* bnd_edge; //(NumBndEdges)

  // element areas
  dTensor1* area_prim; //(NumElems)
  
  // dual element areas
  dTensor1* area_dual; //(NumNodes)

  // edge: list of coordinates (x1,y1) and (x2,y2) that make up edge (size = NumEdges-by-4)
  dTensor2* edge; //(NumEdges,4);

  // tedge: list of edges attached to element (size = NumElems-by-3)
  iTensor2* tedge; //(NumElems,3);
  
  // tedge: list of orientiation (either +1 or -1) for edges attached to element (size = NumElems-by-3)
  iTensor2* tedge_orientation; //(NumElems,3);

  // eelem: list of elements attached to edge (size = NumEdges-by-2)
  iTensor2* eelem; //(NumEdges,2);

  // enode: list of nodes attached to edge (size = NumEdges-by-2)
  iTensor2* enode; //(NumEdges,2);

  // Jmat: Jacobian transformation matrix for mapping derivatives
  //       in the canonical triangle basis (xi,eta) to derivatives
  //       in the Cartesian coordinates (x,y)
  dTensor3* Jmat; //(NumElems,2,2);

  // A list that stores for each node the number of elements
  // that are attached to this node
  iTensor1* NumElemsPerNode; //(NumNodes)

  // adjacent: list of elements that are adjacent to current element
  iTensor2* adjacent;

  // Skel2Elem: matrix to convert a vector defined on the
  // mesh skeleton to an element-center quantity
  bool skel2elem_on;
  dTensor3* Skel2Elem;

  // KMI: matrix to convert a vector node information to
  // edge-based quantities
  bool kmi_on;
  dTensor3* KMI;
  
  // ----------------------------------------------------------
  // Sub-mesh information
  // ----------------------------------------------------------
  bool is_submesh;
  int SubFactor;
  int SubNumElems;
  int SubNumPhysElems;
  int SubNumPhysNodes;
  int SubNumBndNodes;

  dTensor2* sub_node;       //(SubNumPhysNodes,2)
  iTensor2* sub_tnode;      //(SubNumPhysElems,3)
  iTensor1* sub_bnd_node;   //(SubNumBndNodes)
  dTensor1* sub_area_prim;  //(SubNumPhysElems)
  iTensor2* elem_subs;      //(NumPhysElems,SubFactor^2)
  iTensor2* node_subs;      //(NumPhysElems,((SubFactor+1)*(SubFactor+2))/2)
  
  void QuickSort_int(iTensor1& a, iTensor1& index, int lo, int hi);
  int Unique_int(const iTensor1& a_in,
		 const iTensor1& index_in,
		 iTensor1& a_out, 
		 iTensor1& index_out, 
		 iTensor1& rev_index,
		 int lo, 
		 int hi);
  void QuickSort_double(dTensor1& a, iTensor1& index, int lo, int hi);
  int Unique_double(const dTensor1& a_in,
		    const iTensor1& index_in,
		    dTensor1& a_out, 
		    iTensor1& index_out, 
		    iTensor1& rev_index,
		    int lo, 
		    int hi);
  // ----------------------------------------------------------
};
// ---------------------------------------------------------

/*
*	Begin meshinputdata.h
*/

void MeshInputData(double&, int&);

void MeshCreateIcosa(const double radius, int& numpts, int& numtri,
	point*& p, triangle*& t);

void MeshSubdivide(const int level, const double radius,
	int& numpts, int& numtri, point*& p, triangle*& t);

void MeshSphereCoords(const int, const point* p, point*& psphere);

void MeshEdgeData(const int& numpts,
	const int& numtri,
	const point p[],
	const point psphere[],
	const triangle t[],
	const double area[],
	int& numedges,
	dTensor2*& edge,
	iTensor2*& tedge,
	iTensor2*& eelem);

void MeshOrientEdge(mesh& Mesh);

void MeshComputeAreas(const int numtri,
	const int numpts,
	const point p[],
	triangle t[],
	double*& area,
	double*& cdual);

void MeshStore(point p[],
	triangle t[],
	int bnd_node[],
	int ghost_link[],
	double area[],
	double cdual[],
	dTensor2*& edge,
	iTensor2*& tedge,
	iTensor2*& eelem,
	mesh& Mesh);

void QuickSort_bar(bar* b, int* key, int i, int j);
void Swap_bar(bar* b, int* key, int i, int j);
int Partition_bar(bar* b, int* key, int l, int r, int pivot);

void QuickSort(double*& a, int*& index, int lo, int hi);

#endif
