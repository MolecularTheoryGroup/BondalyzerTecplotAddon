/****************************************************************
 * Allows for creation and analysis
 * of a node-edge graph for critical
 * point bookkeeping in Bondalyzer
 *
 * The key functions are AddNode, AddEdge and CheckConnectivity 
 * (look in respective comments for more info on usage)
 *
 * The CritPointGraph class maintains a vector of Node structs, 
 * each with an ID and a list of other Nodes which it connects 
 * to. The nodes are kept in a sorted order (by ID) although 
 * this may not be necessary; the list of edges in each node is 
 * also kept in sorted order.
 *
 * There is a main() function at the bottom of this file for 
 * debug purposes, can be removed in final implementation
 *
 * WARNING : not extensively tested
 * WARNING : EXCEPTIONS CURRENTLY DISABLED!!!!!
 *
 * @author Jan Durakiewicz jdurakie@mymail.mines.edu
 *******************************/
#pragma once
#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <stdlib.h>
#include <time.h>

/*EXCEPTIONS*/
//int NodeAlreadyExists = 10;
//int EdgeAlreadyExists = 11;
//int NoSuchNode = 20;

/****
 * Contains information about a node in a graph : 
 * the node ID, and IDs of nodes it connects to
 ****/
struct Node
{
	unsigned long ID;
	std::vector<unsigned long> edges;
	/*** TODO : payload? ***/
};

/****
 * Manages a list of structs and allows operations on them
 ****/
class CritPointGraph
{
public:
	CritPointGraph();
	~CritPointGraph();

	/****
	* Actual list of nodes;
	* kept sorted by node ID
	****/
	std::vector<Node> graph;

	/****
	* Adds a new node to the graph with given ID
	* @param _ID The ID of the new node to add
	* @throws NodeAlreadyExists if a node with specified ID already exists
	* @author Jan Durakiewicz
	****/
	void AddNode(unsigned long _ID);

	/****
	* Adds a new edge between specified nodes
	* @param _ID1 ID of one edge endpoint node
	* @param _ID2 ID of other edge endpoint node
	* @throws EdgeAlreadyExists if the edge already exists
	* @throws NoSuchNode if one of the endpoint nodes doesn't exist
	* @author Jan Durakiewicz
	****/
	void AddEdge(unsigned long _ID1, unsigned long _ID2);

	/****
	* Checks if two nodes are at a distance of 2 or less 
	* (i.e. if they are directly connected by an edge, or share a neighbor)
	* @param _ID1 first node in pair to test
	* @param _ID2 second node in pair to test
	* @throws NoSuchNode if one of the nodes doesn't exist
	* @return true if the two nodes are connected, false otherwise
	* @author Jan Durakiewicz
	****/
	bool CheckConnectivity(unsigned long _ID1, unsigned long _ID2);

	/****
	* For debugging purposes, prints all nodes and outward edges at each
	* @author Jan Durakiewicz
	****/
	void PrintNodes();

	/****
	* For debugging purposes, dumps all nodes and edges to a CSV file
	* @author Jan Durakiewicz
	****/
	void WriteNodes();

private:
	/****
	* Returns a pointer to a node struct with given ID (if it exists)
	* @param _ID the ID of the node to retrieve
	* @throws NoSuchNode if no node with specified ID exists
	* @return pointer to node struct with given ID, or NULL if it doesn't exist
	* @author Jan Durakiewicz
	****/
	Node* GetNode(unsigned long _ID);

	/****
	* Adds an outgoing edge to a node
	* @param _node pointer to node to add edge to
	* @param _ID the ID of the node that the edge connects this node to
	* @throws EdgeAlreadyExists if the edge is already present from this node
	* @author Jan Durakiewicz
	****/
	void AddEdgeToNode(Node* _node, unsigned long _ID);

	/****
	* Checks if the given node has an edge going to a target node
	* note : just because this returns true does not mean that the target node
	* exists, just that it is listed in the input node's list of edges
	* @param _node node to check for outgoing edge
	* @param _targetID the ID of the node that we're looking for an edge to
	* @return true if the edge exists, false otherwise
	* @author Jan Durakiewicz
	****/
	bool CheckEdgeFromNode(Node* _node, unsigned long _targetID);
};



/*****
 * CODE BODY
 ******/
//CritPointGraph::CritPointGraph()
//{
//}
//
//
//CritPointGraph::~CritPointGraph()
//{
//}
//
//void CritPointGraph::AddNode(unsigned long _ID)
//{
//
//	if (graph.size() >/*!*/ 0)
//	{
//		//the following could be a binary search style insert
//		for (std::vector<Node>::iterator nodeIter = graph.begin(); nodeIter != graph.end(); ++nodeIter)
//		{
//			if (nodeIter->ID == _ID)
//			{
//				throw(NodeAlreadyExists);
//				return;
//			}
//			if (nodeIter->ID > _ID)
//			{
//				Node* newNode = new Node;
//				newNode->ID = _ID;
//				graph.insert(nodeIter, *newNode);
//				std::cout << "Node " << _ID << " added!" << std::endl;
//				return;
//			}
//		}
//	}
//	//happens by default, if the graph is empty or if the new node goes at the end (i.e. doesn't
//	//have to be insterted)
//	Node* newNode = new Node;
//	newNode->ID = _ID;
//	graph.push_back(*newNode);
//}
//
//
//void CritPointGraph::AddEdge(unsigned long _ID1, unsigned long _ID2)
//{
//	Node* node1 = GetNode(_ID1);
//	Node* node2 = GetNode(_ID2);
//
//	AddEdgeToNode(node1, _ID2);
//	AddEdgeToNode(node2, _ID1);
//}
//
//bool CritPointGraph::CheckConnectivity(unsigned long _ID1, unsigned long _ID2)
//{
//	Node* node1 = GetNode(_ID1);
//	Node* node2 = GetNode(_ID2);
//	//first, check if the two nodes are connected to each other
//	if (CheckEdgeFromNode(node1, node2->ID))
//	{
//		return true;
//	}
//
//	//then, check if they have any neighbors in common
//	for (std::vector<unsigned long>::iterator edgeIter = node1->edges.begin(); edgeIter != node1->edges.end(); ++edgeIter)
//	{
//		if (CheckEdgeFromNode(node2, *edgeIter))
//			return true;
//	}
//
//	//if both of the above tests failed, there's no connection
//	return false;
//
//}
//
//void CritPointGraph::PrintNodes()
//{
//	for (std::vector<Node>::iterator nodeIter = graph.begin(); nodeIter != graph.end(); ++nodeIter)
//	{
//		std::cout << (*nodeIter).ID << std::endl;
//		for (std::vector<unsigned long>::iterator edgeIter = nodeIter->edges.begin(); edgeIter != nodeIter->edges.end(); ++edgeIter) //ew
//		{
//			std::cout << "\t " << *edgeIter << std::endl;
//		}
//	}
//	std::cout << std::endl;
//}
//
//void CritPointGraph::WriteNodes()
//{
//	std::ofstream file;
//	file.open("DUMP.csv");
//	for (std::vector<Node>::iterator nodeIter = graph.begin(); nodeIter != graph.end(); ++nodeIter)
//	{
//		//int nodeID = (*nodeIter).ID;
//		file << (*nodeIter).ID + 1;
//
//		for (std::vector<unsigned long>::iterator edgeIter = nodeIter->edges.begin(); edgeIter != nodeIter->edges.end(); ++edgeIter) //ew
//		{
//			file << ", " << (*edgeIter) + 1;
//		}
//		file << std::endl;
//	}
//
//	for (std::vector<Node>::iterator nodeIter = graph.begin(); nodeIter != graph.end(); ++nodeIter)
//	{
//		int Node1 = (*nodeIter).ID;
//		file << Node1 + 1;
//		file << " connected to 5: ";
//		file << std::boolalpha << CheckConnectivity(4, Node1) << std::endl;
//	}
//	file.close();
//}
//
//Node* CritPointGraph::GetNode(unsigned long _ID)
//{
//	//search nodes for one with given ID
//	for (std::vector<Node>::iterator nodeIter = graph.begin(); nodeIter != graph.end(); ++nodeIter)
//	{
//		if (nodeIter->ID == _ID)
//		{
//			return &*nodeIter; // boy that's dumb looking
//		}
//	}
//	throw NoSuchNode;
//	return NULL;
//}
//
//void CritPointGraph::AddEdgeToNode(Node* _node, unsigned long _ID)
//{
//	//works almost identically to AddNode, within the confines of a list of edges
//	for (std::vector<unsigned long>::iterator edgeIter = _node->edges.begin(); edgeIter != _node->edges.end(); ++edgeIter)
//	{
//		if (*edgeIter == _ID)
//		{
//			//throw EdgeAlreadyExists;
//			return;
//		}
//		if (*edgeIter > _ID)
//		{
//			_node->edges.insert(edgeIter, _ID);
//			return;
//		}
//	}
//	_node->edges.push_back(_ID);
//}
//
//bool CritPointGraph::CheckEdgeFromNode(Node* _node, unsigned long _targetID)
//{
//	for (std::vector<unsigned long>::iterator edgeIter = _node->edges.begin(); edgeIter != _node->edges.end(); ++edgeIter)
//	{
//		if (*edgeIter == _targetID)
//		{
//			return true;
//		}
//	}
//
//	return false;
//}



///***
// DEBUG
// ***/
//int main()
//{
//	srand(time(0));
//	CritPointGraph CPGraph;
//	//add some random nodes
//	for (int i = 0; i < 30; i++)
//	{
//		int num = rand() % 10;
//		std::cout << "Adding " << num << std::endl;
//		try{
//			CPGraph.AddNode(num);
//		}
//		catch (int NodeAlreadyExists)
//		{
//			std::cout << "Node with given ID (" << num <<") already exists!" << std::endl;
//		}
//	}
//	//add some random connections between the nodes
//	for (int i = 0; i < 20; i++)
//	{
//		int node1 = rand() % 10;
//		int node2 = rand() % 10;
//		if (node1 == node2)
//			continue;
//		try
//		{
//			CPGraph.AddEdge(node1, node2);
//		}
//		catch (int NoSutchNode)
//		{
//			
//		}
//	}
//
//	CPGraph.PrintNodes();
//	
//	//check if some random nodes are neighbors or not
//	for (int i = 0; i < 10; i++)
//	{
//		int node1 = rand() % 10;
//		int node2 = rand() % 10;
//		if (node1 == node2)
//			continue;
//		try
//		{
//			std::cout << "Checking connectivity between " << node1 << " and " << node2 << " : "
//				<< CPGraph.CheckConnectivity(node1, node2) << std::endl;
//		}
//		catch (int NoSutchNode)
//		{
//		}
//	}
//
//	return 0;
//}