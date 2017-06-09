#include "CRITPOINTGRAPH.h"


CritPointGraph::CritPointGraph()
{
}


CritPointGraph::~CritPointGraph()
{
}

void CritPointGraph::AddNode(unsigned long _ID)
{

	if (graph.size() >/*!*/ 0)
	{
		//the following could be a binary search style insert
		for (std::vector<Node>::iterator nodeIter = graph.begin(); nodeIter != graph.end(); ++nodeIter)
		{
			if (nodeIter->ID == _ID)
			{
				//throw(NodeAlreadyExists);
				return;
			}
			if (nodeIter->ID > _ID)
			{
				Node* newNode = new Node;
				newNode->ID = _ID;
				graph.insert(nodeIter, *newNode);
				std::cout << "Node " << _ID << " added!" << std::endl;
				return;
			}
		}
	}
	//happens by default, if the graph is empty or if the new node goes at the end (i.e. doesn't
	//have to be insterted)
	Node* newNode = new Node;
	newNode->ID = _ID;
	graph.push_back(*newNode);
}


void CritPointGraph::AddEdge(unsigned long _ID1, unsigned long _ID2)
{
	Node* node1 = GetNode(_ID1);
	Node* node2 = GetNode(_ID2);

	AddEdgeToNode(node1, _ID2);
	AddEdgeToNode(node2, _ID1);
}

bool CritPointGraph::CheckConnectivity(unsigned long _ID1, unsigned long _ID2)
{
	Node* node1 = GetNode(_ID1);
	Node* node2 = GetNode(_ID2);
	//first, check if the two nodes are connected to each other
	if (CheckEdgeFromNode(node1, node2->ID))
	{
		return true;
	}

	//then, check if they have any neighbors in common
	for (std::vector<unsigned long>::iterator edgeIter = node1->edges.begin(); edgeIter != node1->edges.end(); ++edgeIter)
	{
		if (CheckEdgeFromNode(node2, *edgeIter))
			return true;
	}

	//if both of the above tests failed, there's no connection
	return false;

}

void CritPointGraph::PrintNodes()
{
	for (std::vector<Node>::iterator nodeIter = graph.begin(); nodeIter != graph.end(); ++nodeIter)
	{
		std::cout << (*nodeIter).ID << std::endl;
		for (std::vector<unsigned long>::iterator edgeIter = nodeIter->edges.begin(); edgeIter != nodeIter->edges.end(); ++edgeIter) //ew
		{
			std::cout << "\t " << *edgeIter << std::endl;
		}
	}
	std::cout << std::endl;
}

void CritPointGraph::WriteNodes()
{
	std::ofstream file;
	file.open("DUMP.csv");
	for (std::vector<Node>::iterator nodeIter = graph.begin(); nodeIter != graph.end(); ++nodeIter)
	{
		//int nodeID = (*nodeIter).ID;
		file << (*nodeIter).ID + 1;

		for (std::vector<unsigned long>::iterator edgeIter = nodeIter->edges.begin(); edgeIter != nodeIter->edges.end(); ++edgeIter) //ew
		{
			file << ", " << (*edgeIter) + 1;
		}
		file << std::endl;
	}

	for (std::vector<Node>::iterator nodeIter = graph.begin(); nodeIter != graph.end(); ++nodeIter)
	{
		int Node1 = (*nodeIter).ID;
		file << Node1 + 1;
		file << " connected to 5: ";
		file << std::boolalpha << CheckConnectivity(4, Node1) << std::endl;
	}
	file.close();
}

Node* CritPointGraph::GetNode(unsigned long _ID)
{
	//search nodes for one with given ID
	for (std::vector<Node>::iterator nodeIter = graph.begin(); nodeIter != graph.end(); ++nodeIter)
	{
		if (nodeIter->ID == _ID)
		{
			return &*nodeIter; // boy that's dumb looking
		}
	}
	//throw NoSuchNode;
	return NULL;
}

void CritPointGraph::AddEdgeToNode(Node* _node, unsigned long _ID)
{
	//works almost identically to AddNode, within the confines of a list of edges
	for (std::vector<unsigned long>::iterator edgeIter = _node->edges.begin(); edgeIter != _node->edges.end(); ++edgeIter)
	{
		if (*edgeIter == _ID)
		{
			//throw EdgeAlreadyExists;
			return;
		}
		if (*edgeIter > _ID)
		{
			_node->edges.insert(edgeIter, _ID);
			return;
		}
	}
	_node->edges.push_back(_ID);
}

bool CritPointGraph::CheckEdgeFromNode(Node* _node, unsigned long _targetID)
{
	for (std::vector<unsigned long>::iterator edgeIter = _node->edges.begin(); edgeIter != _node->edges.end(); ++edgeIter)
	{
		if (*edgeIter == _targetID)
		{
			return true;
		}
	}

	return false;
}