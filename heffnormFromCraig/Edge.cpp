#include <utility>
#include "Edge.h"

Edge MakeEdge(int e1, int e2, bool DoSort) 
{ 
	if (DoSort)
		return (e1 < e2 ? std::make_pair(e1, e2) : std::make_pair(e2, e1));
	else
		return std::make_pair(e1, e2);
}

Edge MakeEdge(Edge e, bool DoSort)
{ 
	return MakeEdge(e.first, e.second, DoSort); 
}