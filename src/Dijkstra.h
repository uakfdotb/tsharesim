#ifndef _DIJKSTRA_H
#define _DIJKSTRA_H

#include <vector>
#include <iostream>
#include <queue>

#include "Graph.h"

using namespace std;

//#define DIJKSTRA_DEBUG

class Dijkstra{
private:
	Graph& graph;
	
	VertexList nodelist;
	unsigned int nodesize;
public:
	Dijkstra(Graph& g);
	virtual ~Dijkstra();
	
	double run(VertexID, VertexID, vector<EdgeID>&);
};

#endif
