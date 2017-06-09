
#include "TECADDON.h"
#include "ADDGLBL.h"
#include "GUIDEFS.h"
#if !defined (MSWIN)
#include <unistd.h>
#endif

//#include <cstdlib>
#include <sstream>
#include <string>
#include <vector>

#include "CSMDATATYPES.h"
#include "ZONEVARINFO.h"

#include "MAINFUNCTION.h"
#include "MeshData.h"

using std::vector;
using std::string;
using std::stringstream;

/*
 *	Quick breakdown of the Tecplot type terminology:
 *	 Something_t:		Typedef, typically reduces to a basic datatype like int, double, etc...
 *	 Something_e:		enum
 *	 Something_pa:		pointer to some advanced data type (like a struct) or reference pointer to Tecplot dataset data
 *	 Something_s:		struct
 *	 
 *	 There are others, but not in this file.
 *	 
 *	 I only know which types to use for variables based on return and argument types for TecUtil functions, so RTFM!
 */

void Mesher()
{
	MeshData meshData("U:\\Workspace\\small_sphere.dat");
	//MeshData meshData("D:/Files/Dev/Work/Mesher/sphere.dat");
	//MeshData meshData("C:/Users/Tim Wilson/Desktop/sphere.dat");
	//meshData.debug_print_data();
	//meshData.debug_print_points_from_vector();
	//meshData.debug_print_triangles();
	//meshData.debug_print_connections();
	meshData.debug_save_tecplot_zone(-1);
	meshData.debug_save_tecplot_zone(0);
	meshData.test_build_tets();
	meshData.debug_save_tecplot_zone(1);
}

void MainFunction(){
	TecUtilDialogMessageBox("Hi Jan! This is your life now...", MessageBox_Warning);

	if (!TecUtilDataSetIsAvailable()){
		TecUtilDialogMessageBox("There is no dataset available.", MessageBox_Error);
		return;
	}

	EntIndex_t NumZones = TecUtilDataSetGetNumZones();
	EntIndex_t NumVars = TecUtilDataSetGetNumVars();

	vector<EntIndex_t> XYZRhoVarNums(4, -1);

	/*
	 *	The first call gets the X,Y,Z variable numbers. The second call
	 *	gets the electron density variable number.
	 *	Note that VarNumByName is sort of a wrapper function I wrote to 
	 *	keep from dealing with the normal way of getting a variable
	 *	number according to its name.
	 *  Also note that all four variable numbers are being stored in
	 *  the same vector.
	 */
	TecUtilAxisGetVarAssignments(&XYZRhoVarNums[0], &XYZRhoVarNums[1], &XYZRhoVarNums[2]);
	XYZRhoVarNums[3] = VarNumByName(string("Electron Density"));
	if (XYZRhoVarNums[3] <= 0){
		TecUtilDialogMessageBox("Failed to get electron density var num.", MessageBoxType_Error);
		return;
	}

	/*
	 *	IGNORE MEEEEE!!!!!
	 */
	string AuxDataCheckStr = "GBA.ZoneType";
	string CheckStr = "FEVolumeZone";

	/*
	 *	Loop over every single zone in the system, regardless
	 *	of its type or whether or not it's active.
	 */
	for (EntIndex_t ZoneNum = 1; ZoneNum <= NumZones; ++ZoneNum){
		/*
		 *  Very quick way of skipping over most of the zones,
		 *  which are not finite element or FEQuad.
		 */
		if (TecUtilZoneIsFiniteElement(ZoneNum) && TecUtilZoneGetType(ZoneNum) == ZoneType_FEQuad){
			
			/*
			 * KEEP IGNORING MEEEE!!!!
			 * I'm checking that the FE zone was created by gradient bundle analysis 
			 * (If the TecUtilAuxDataGetItemByName() retuns false, then it wasn't made
			 * by gradient bundle analysis) and that it's the type of volume we're interested in.
			 * More importantly, WHY AREN'T YOU IGNORING ME!?!
			 */
			Boolean_t ZoneOK = TRUE;
			AuxData_pa TempAuxData = TecUtilAuxDataZoneGetRef(ZoneNum);
			ZoneOK = VALID_REF(TempAuxData);
			for (int i = 0; i < AuxDataCheckStr.size() && ZoneOK; ++i){
				char* TempCStr;
				AuxDataType_e ADTJunk;
				Boolean_t BJunk;
				ZoneOK = TecUtilAuxDataGetItemByName(TempAuxData, AuxDataCheckStr.c_str(), (ArbParam_t*)&TempCStr, &ADTJunk, &BJunk);
				if (ZoneOK){
					ZoneOK = (string(TempCStr) == CheckStr);
				}
				TecUtilStringDealloc(&TempCStr);
			}

			/*
			 *	If the zone passed both checks...(you can stop ignoring me now)
			 */
			if (ZoneOK){
				LgIndex_t NumElements = -1;
				LgIndex_t NumNodes = -1;
				LgIndex_t NumNodesPerElem = -1;

				TecUtilZoneGetIJK(ZoneNum, &NumNodes, &NumElements, &NumNodesPerElem);

				ZoneName_t ZoneName;
				/*
				 * Note that this string is deallocated below
				 */
				TecUtilZoneGetName(ZoneNum, &ZoneName);

				stringstream SS;
				SS << "Zone " << ZoneName << " with " << NumElements << " elements (" << NumNodesPerElem
					<< " nodes for each) and " << NumNodes << " nodes total.";

				TecUtilStringDealloc(&ZoneName);

				if (!TecUtilDialogMessageBox(string(SS.str() + string(", Continue?")).c_str(), MessageBoxType_YesNo))
					break;

				/*
				 *	METHOD 1
				 */

				vector<CSM_Vec4_s> FullNodeInfo(NumNodes);

				/*
				 * Looping over all nodes (vertices) of the FE surface according to
				 * node order. This order has no meaning except for the location of the
				 * node's data (position, etc...).
				 */
				for (LgIndex_t NodeNum = 1; NodeNum <= NumNodes; ++NodeNum){
					SS.str(string());
					SS << "Node " << NodeNum << ": {X, Y, Z, Rho} = {";

					/*
					 *	A simplified way to get data from Tecplot:
					 *		Store the variable numbers you need in a vector,
					 *		then loop the GetByZoneVar call, storing the results
					 *		in another vector.
					 */
					for (int i = 0; i < 4; ++i){
						FullNodeInfo[NodeNum - 1][i] = TecUtilDataValueGetByZoneVar(ZoneNum, XYZRhoVarNums[i], NodeNum);

						SS << FullNodeInfo[NodeNum - 1][i];

						if (i < 3) SS << ", ";
					}

					SS << "}";

					if (!TecUtilDialogMessageBox(string(SS.str() + string(", Continue?")).c_str(), MessageBoxType_YesNo))
						break;
				}

				/*
				 *	METHOD 2
				 */

				/*
				 *	Looping over each element of the FE zone according to
				 *	element number. Note that only node index information can
				 *	be retrieved by doing this. The node indices here can be used
				 *	to retrieve actual positions and other values for nodes
				 *	as above.
				 */
				for (LgIndex_t ElemNum = 1; ElemNum <= NumElements; ++ElemNum){
					/*
					 *	A vector to store node indices.
					 *	Note that a NodeMap_t is just an integer, so
					 *	this is just a vector of int to store indices.
					 */
					vector<NodeMap_t> Nodes(NumNodesPerElem);

					SS.str(string());
					SS << "Element " << ElemNum << " has nodes {";

					/*
					 *  Looping over each corner of the element, storing the corner's
					 *  node index.
					 */
					for (LgIndex_t NodeNum = 1; NodeNum <= NumNodesPerElem; ++NodeNum){
						Nodes[NodeNum - 1] = TecUtilDataNodeGetByZone(ZoneNum, ElemNum, NodeNum);
						SS << Nodes[NodeNum - 1];

						if (NodeNum < NumNodesPerElem) SS << ", ";
					}

					SS << "}";

					if (!TecUtilDialogMessageBox(string(SS.str() + string(", Continue?")).c_str(), MessageBoxType_YesNo))
						break;
				}


				/*
				 *	METHOD 3
				 */

				/*
				 *	This is for when you want to navigate the FE surface by node. These
				 *	functions do the work of taking a node index and finding the elements
				 *	to which it belongs.
				 *	This type of navigation won't be necessary thankfully, as there's
				 *	no way to do this without constant calls to tecplot (except by
				 *	getting all the connectivity information yourself and writing the
				 *	functions to navigate it yourself.
				 *
				 *	Note that this is only possible by first getting the NodeToElemMap_pa
				 *	pointer for the zone, and then referencing that when making queries.
				 */
				NodeToElemMap_pa NodeToElemMap = TecUtilDataNodeToElemMapGetReadableRef(ZoneNum);

				for (LgIndex_t NodeNum = 1; NodeNum <= NumNodes; ++NodeNum){
					/*
					 *	Vector for storing the element numbers of a particular node.
					 *	Sized according to the number of elements for the node.
					 */
					vector<LgIndex_t> Elem(TecUtilDataNodeToElemMapGetNumElems(NodeToElemMap, NodeNum));

					SS.str(string());
					SS << "Node " << NodeNum << " is part of " << Elem.size() << " elements: {";

					/*
					 *	Looping over the elements of a node to get their numbers.
					 *	If you then wanted to find the node's neighboring nodes then
					 *	you'd use the second method to get the node's element's node's.
					 */
					for (LgIndex_t ElemNum = 1; ElemNum <= Elem.size(); ++ElemNum){
						Elem[ElemNum - 1] = TecUtilDataNodeToElemMapGetElem(NodeToElemMap, NodeNum, ElemNum);
						SS << Elem[ElemNum - 1];

						if (ElemNum < Elem.size()) SS << ", ";
					}

					SS << "}";

					if (!TecUtilDialogMessageBox(string(SS.str() + string(", Continue?")).c_str(), MessageBoxType_YesNo))
						break;
				}
			}
		}
	}
}