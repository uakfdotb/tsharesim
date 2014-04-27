#include "Dijkstra.h"

struct QueueComp {
	bool operator()( const pair<double, VertexID>& a, const pair<double, VertexID>& b ) const {
		return ( a.first >= b.first );
	}
};

Dijkstra::Dijkstra(Graph& g):graph(g){
	VertexList(0).swap(nodelist);
	for (int i = 0; i < graph.vertices().size(); i++) {
		nodelist.push_back(graph[i]);
	}
	nodesize = graph.num_vertices();
}
Dijkstra::~Dijkstra(){}
	
double Dijkstra::run(VertexID source, VertexID target, vector<EdgeID>& innerIDs){
	//cout << "Test: source = " << source << ", target = " << target << "." << endl; 
	vector<double> dist (nodesize);
	for ( int vid = 0; vid < nodesize; vid++) {
		dist[vid] = DBL_INFINITY;
	}
	dist[source] = 0.0;
	
	// for debugging;
	//#ifdef DIJKSTRA_DEBUG
		vector<VertexID> father(nodesize, -1);
	//#endif
	
	priority_queue<pair<double, VertexID>, vector<pair<double, VertexID> >, QueueComp> Queue;
	Queue.push(make_pair(0.0, source));
	
	for ( int vid = 0; vid < nodesize; vid++) {
		nodelist[vid].flaga = false;
	}
	
	int num_visited = 0;
	while (Queue.size() > 0) {
		double min_dist = Queue.top().first;
		VertexID vid = Queue.top().second;
		#ifdef DIJKSTRA_DEBUG
			cout << "vid = " << vid << ", min_dist = " << min_dist << endl;
		#endif
		if (vid == target) break;
		Queue.pop();
		if (min_dist > dist[vid]) {/*cout << "Greater!?" << endl;*/ continue;} // lazy update; 
		else nodelist[vid].flaga = true; // settle vertex vid;
		
		num_visited++;
		
		forall_outneighbors(graph, vid, eit){
			if (nodelist[eit->target].flaga) continue; // settled;
			if (dist[eit->target] > min_dist + eit->weight) {
				dist[eit->target] = min_dist + eit->weight;
				Queue.push(make_pair(dist[eit->target], eit->target)); 
				//#ifdef DIJKSTRA_DEBUG
					father[eit->target] = vid;
				//#endif
			}	
		}
		
		// forall_outshortcuts(graph, vid, eit){
			// if (nodelist[eit->target].flaga) continue; // settled;
			// if (dist[eit->target] > min_dist + eit->weight) {
				// dist[eit->target] = min_dist + eit->weight;
				// Queue.push(make_pair(dist[eit->target], eit->target)); 
			// }	
		// }		
	}
	//cout << "num visited dijk = " << num_visited << endl;
	
	
	// print path;
	//#ifdef DIJKSTRA_DEBUG
		vector<int> path;
		path.push_back(target);
		while(true){
			int c = path.back();
			if (father[c] != -1) path.push_back(father[c]);
			else break;
		}
		// cout << "path: ";
		// for ( int i = path.size()-1; i >=0; i-- ) {
			// cout << path[i] << " ";
		// }
		// cout << endl;
	//#endif
	cout << "distance:" << dist[target] << endl;
	
	
	return dist[target];
}
