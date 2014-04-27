#ifndef HD_H
#define HD_H

#include "Graph.h"

#include <queue>

using namespace std;

//#define HD_DEBUG

// highway dimention;
class HD{
private:
	Graph& graph;
	
public:
	HD(Graph&);
	~HD(void);
	
	void findPath(VertexID, VertexID, vector<VertexID>&, double&);
	double run(VertexID, VertexID, VertexID&) const;
};

#endif
