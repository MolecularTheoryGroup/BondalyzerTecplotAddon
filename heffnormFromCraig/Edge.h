#pragma once
#include <utility>
#include "SimpleTypes.h"

typedef std::pair<Index_t, Index_t> Edge;

Edge MakeEdge(int e1, int e2, bool DoSort = true);
Edge MakeEdge(Edge e, bool DoSort = true);