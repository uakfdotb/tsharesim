#ifndef _QUERY_H
#define _QUERY_H

#include "Graph.h"

#include <vector>
#include <queue>

using namespace std;

//#define DEBUG_MODE


class Query{
private:
	Graph& graph;
	
	VertexList nodelist;
	unsigned int nodesize;
public:
	Query(Graph& g);
	virtual ~Query();
	double BiDijkstra(VertexID, VertexID, vector<EdgeID>&);	
	double run(VertexID, VertexID, vector<EdgeID>&);
};

#endif
