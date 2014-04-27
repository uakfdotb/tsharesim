#include "includes.h"
#include "vertex.h"
#include "shortestPath.h"

double ShortestPath :: distance(vertex *a, vertex *b) {
	double d1 = a->x - b->x;
	double d2 = a->y - b->y;
	return sqrt(d1 * d1 + d2 * d2);
}

ShortestPath :: ShortestPath(vector<vertex *> &cvertices) : vertices(cvertices) {
	if(EUCLIDEAN_PATH == 0) {
		distance_cache = new dcache(10000000);
		path_cache = new pcache(10000);
	
		//for statistics
		distance_queries = 0;
		path_queries = 0;
		distance_misses = 0;
		path_misses = 0;
		//total_timer.resume();
	
		const rlim_t sSize = 1024 * 1024 * 512;
		struct rlimit rl;

		if (getrlimit(RLIMIT_STACK, &rl) == 0) {
			if (sSize > rl.rlim_cur){
				rl.rlim_cur = sSize;
				if (setrlimit(RLIMIT_STACK, &rl) != 0) {
					cerr << "ShortestPath: could not increase stack size. Program might run into trouble..." << endl;
				}
			}
		}
	
		string filename = inDirectory + "/edges.gra";
		int query_num = 1000;
		bool unpack = true;
		bool save_para = false;
		bool label_generation = false;
		bool doidmap = true;
	
		cout << "ShortestPath: opening " << filename << endl;
		ifstream infile(filename.c_str());
		if (!infile) {
			cout << "ShortestPath: Error: Cannot open " << filename << endl;
			exit(0);
		}
	
		if (doidmap) { // if has id map, we will map the id back to its original id;
			readIDMap(filename, id_map);
		}
	
		string query_file = filename;
	
		// main part;
		// build graph
		cout << "ShortestPath: building graph" << endl;
		g = new Graph(infile);

		// generate queries;
		cout << "ShortestPath: generating queries" << endl;
		vector<VertexID> query;
		checkQuery(*g, filename, query_num, query, save_para);	

		// build contraction hierarchies;
		cout << "ShortestPath: building construction hierarchies" << endl;
		checkHierarchies(*g, filename, true, false);
	
		//q = new Query(*g); // contraction hierarchy method;
		//dijk = new Dijkstra(*g); // baseline method;
		h = new HD(*g); // highway dimension method;
	
		cout << "ShortestPath: finished loading" << endl;
	}
	
	cout << "ShortestPath: finished loading" << endl;
}

queue<vertex *> ShortestPath :: shortestPath(vertex *start, vertex *end) {
	if(EUCLIDEAN_PATH == 0) {
		if(start->id == end->id) {
			//well, if to the same point then we can of course go directly
			queue<vertex *> equeue;
			equeue.push(end);
			return equeue;
		} else {
			path_queries++;
			
			//search cache first
			uint64_t x = ((uint64_t) start->id) * vertices.size() + (uint64_t) end->id;
			pcache::const_iterator it = path_cache->find(x);
	
			if(it != path_cache->end()) {
				return it.value();
			} else {
				//statistic-related operations
				//path_timer.resume();
				path_misses++;
				
				//the path wasn't found in cache, so we have to calculate it
				//the distance might still be there though
				dcache::const_iterator d_it = distance_cache->find(x);
				
				if(d_it == distance_cache->end()) {
					//well, distance wasn't there
					//call shortestDistance so it calculates it, and retrieve from cache
					// (since we need the meet node, not just distance)
					shortestDistance(start, end);
					d_it = distance_cache->find(x);
				}
	
				distanceElement distance_el = d_it.value();
				
				//we execute findPath twice on closer points to minimize execution time
				vector<VertexID> path_from_source, path_to_target;
				h->findPath(start->id, distance_el.meet_node, path_from_source, distance_el.distance);
				h->findPath(end->id, distance_el.meet_node, path_to_target, distance_el.distance);

				//reconstruct the path from the subpaths
				
				queue<vertex *> actualPath;
				for(int i = path_from_source.size() - 2; i >= 0; i--) {
					actualPath.push(vertices[path_from_source[i]]);
				}
	
				for(int i = 1; i < path_to_target.size(); i++) {
					actualPath.push(vertices[path_to_target[i]]);
				}

				//store a copy of the path into the cache
				//must be copy because actualPath might be modified
				queue<vertex *> queuecopy = actualPath;
				(*path_cache)[x] = queuecopy;
		
				//statistic-related operations
				//path_timer.pause();
		
				if(path_misses % 100 == 0) {
					cout << "distance status:   " << (((float) distance_misses) / distance_queries) << "   " << distance_queries << "   " << distance_timer.currRunTime() << endl;
					cout << "path status:   " << (((float) path_misses) / path_queries) << "   " << path_queries << "   " << path_timer.currRunTime() << endl;
					cout << "total time:   " << total_timer.currRunTime() << endl;
				}
		
				if(actualPath.size() == 0) {
					cout << "uh oh, we have a problem; debug: " << path_from_source.size() << "  " << path_to_target.size() << "  " << distance_el.meet_node << endl;
				}
				
				if(abs(pathDistance(start, actualPath) - shortestDistance(start, end)) > 0.00001) {
					cout << "uh oh: sd = " << pathDistance(start, actualPath) << " and " << shortestDistance(start, end) << "   " << start->id << " -> " << end->id << endl;
					cout << "found path:";
					queue<vertex *> pathCopy = actualPath;
					while(!pathCopy.empty()) {
						cout << " " << pathCopy.front()->id;
						pathCopy.pop();
					}
					
					cout << endl;
				}
		
				return actualPath;
			}
		}
	} else {
		//just return direct path to the end point
		//then taxi will follow this possibly fake edge, like Euclidean distance
		
		queue<vertex *> path;
		path.push(end);
		return path;
	}
}

double ShortestPath :: pathDistance(vertex *start, queue<vertex *> path) {
	if(path.size() == 0) return 0;
	
	double totalDistance = distance(start, path.front());
	while(path.size() > 1) {
		vertex *popped = path.front();
		path.pop();
		
		totalDistance += distance(popped, path.front());
	}
	
	return totalDistance;
}

double ShortestPath :: shortestDistance(vertex *start, vertex *end) {
	if(EUCLIDEAN_PATH == 0) {
		distance_queries++;
		
		//find the unique identifier for this distance
		uint64_t x = ((uint64_t) start->id) * vertices.size() + (uint64_t) end->id;
		
		//search the cache
		dcache::const_iterator it = distance_cache->find(x);
	
		if(it != distance_cache->end()) {
			//found in cache! just return then
			return (*distance_cache)[x].distance;
		} else {
			//the distance wasn't in the cache, so we have to calculate it
			
			//statistic-related operations
			//distance_timer.resume();
			distance_misses++;
		
			distanceElement distance_el;
			distance_el.meet_node = numeric_limits<VertexID>::max();
			distance_el.distance = h->run(start->id, end->id, distance_el.meet_node);
			
			//save to cache
			(*distance_cache)[x] = distance_el;
		
			//statistic-related operations
			//distance_timer.pause();
		
			if(distance_misses % 10000 == 0) {
				cout << "distance status:   " << (((float) distance_misses) / distance_queries) << "   " << distance_queries << "   " << distance_timer.currRunTime() << endl;
				cout << "path status:   " << (((float) path_misses) / path_queries) << "   " << path_queries << "   " << path_timer.currRunTime() << endl;
				cout << "total time:   " << total_timer.currRunTime() << endl;
			}
			
			if(distance_el.distance > 100000) cout << "oh no: " << start->id << "," << end->id << endl;
	
			return distance_el.distance;
		}
	} else {
		return distance(start, end);
	}
}
