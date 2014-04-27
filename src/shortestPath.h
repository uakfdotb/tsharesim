#ifndef SHORTESTPATH_H
#define SHORTESTPATH_H

#include "lru_cache.h"
#include "Contraction.h"
#include "Query.h"
#include "Dijkstra.h"
#include "PerformanceTimer.h"
#include "GraphUtil.h"
#include "HD.h"
#include "graphlink.h"

struct vertex;

struct distanceElement {
	double distance;
	VertexID meet_node;
};

class ShortestPath {
private:
	Graph *g;
	Query *q;
	Dijkstra *dijk;
	HD *h;
	vector<VertexID> id_map;
	vector<vertex *> &vertices;
	
	typedef plb::LRUCacheH4<uint64_t, distanceElement> dcache;
	typedef plb::LRUCacheH4<uint64_t, queue<vertex *> > pcache;
	
	dcache *distance_cache; //distance cache
	pcache *path_cache; //path cache
	
	PerformanceTimer total_timer;
	PerformanceTimer distance_timer;
	PerformanceTimer path_timer;
	
	long long distance_queries;
	long long path_queries;
	long long distance_misses;
	long long path_misses;

	double distance(vertex *a, vertex *b);

public:
	ShortestPath(vector<vertex *> &cvertices);
	queue<vertex *> shortestPath(vertex *start, vertex *end);
	double pathDistance(vertex *start, queue<vertex *> path); //finds the length of a path returned by shortestPath
	double shortestDistance(vertex *start, vertex *end);
};

#endif
