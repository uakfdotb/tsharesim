#include "Contraction.h"

struct QueueComp {
	bool operator()( const pair<double, VertexID>& a, const pair<double, VertexID>& b ) const {
		return ( a.first >= b.first );
	}
};
struct myQueueComp {
	bool operator()( const pair<int, VertexID>& a, const pair<int, VertexID>& b ) const {
		return ( a.first >= b.first );
	}
};

Contraction::Contraction(Graph & g):graph(g){
	VertexList(0).swap(nodelist);
	nodesize = graph.num_vertices();
	for ( int i = 0; i < nodesize; i++ ) {
		nodelist.push_back(graph[i]);
	}

	//tmpinshortcutlist.resize(nodesize);
	tmpoutshortcutlist.resize(nodesize);

	pathIndex = vector<vector<int> >(nodesize, vector<int>(0));
	pathList = vector<vector<VertexID> >(nodesize, vector<VertexID>(0));

	srand(time(NULL));
}

Contraction::~Contraction(){}


/* this function tries to find alternative shortest paths between u and v with two conditions:
	1. without through nid;
	2. within the distance limit.
*/
bool Contraction::dodgingDijkstra(VertexID source, VertexID target, VertexID nid, double limit, double& real_distance) {
	forall_nodes(graph, vid){
		nodelist[vid].flaga = false;
	}

	//cout << "Test: source = " << source << ", target = " << target << "." << endl;
	vector<float> dist (nodesize);
	forall_nodes(graph, vid){
		dist[vid] = DBL_INFINITY;
	}
	dist[source] = 0.0;

	vector<VertexID> father(nodesize, numeric_limits<VertexID>::max());

	priority_queue<pair<double, VertexID>, vector<pair<double, VertexID> >, QueueComp> Queue;
	Queue.push(make_pair(0.0, source));

	while (Queue.size() > 0) {
		float min_dist = Queue.top().first;
		VertexID vid = Queue.top().second;
		//cout << "vid = " << vid << ", min_dist = " << min_dist << endl;
		if (vid == target) break;
		Queue.pop();
		if (min_dist > dist[vid]) 	continue; 					// lazy update;
		else						nodelist[vid].flaga = true; // settle vertex vid;

		if (min_dist >= limit) break;

		forall_outneighbors(graph, vid, eit){
			if (nodelist[eit->target].flaga) {/*cout << "$ ";*/ continue;} // settled;
			if (dist[eit->target] > min_dist + eit->weight) {
				dist[eit->target] = min_dist + eit->weight;
				Queue.push(make_pair(dist[eit->target], eit->target));
				father[eit->target] = vid; // for test;
			}
		}

		for (int i = 0; i < tmpoutshortcutlist[vid].size(); i++) {
			int tmptarget = tmpoutshortcutlist[vid][i].target;
			if (nodelist[tmptarget].flaga) {/*cout << "$ ";*/ continue;}
			if (dist[tmptarget] > min_dist + tmpoutshortcutlist[vid][i].weight){
				dist[tmptarget] = min_dist + tmpoutshortcutlist[vid][i].weight;
				Queue.push(make_pair(dist[tmptarget], tmptarget));
				father[tmptarget] = vid; // for test;
			}
		}
	}

	// check whether nid has been used or not;
	int c = target;
	bool used = false;
	while (father[c] != -1){
		c = father[c];
		if(c == nid) {used = true; break;}
	}
	real_distance = dist[target];
	if (used) return false; // if the middle node is used;


	forall_nodes(graph, vid){
		nodelist[vid].flaga = false;
	}
	// for test;
	#ifdef CONTRACTION_DEBUG
		vector<int> path;
		path.push_back(target);
		while(true){
			int c = path.back();
			if (father[c] != -1) path.push_back(father[c]);
			else break;
		}
		cout << "path: ";
		for ( int i = path.size()-1; i >=0; i-- ) {
			cout << path[i] << " ";
		}
		cout << endl;
		cout << "dist: " << dist[target] << ", limit: " << limit << endl;
	#endif

	if (dist[target] <= limit) { /*cout << "T" << endl;*/ return true;}
	else { /*cout << "F" << endl;*/ return false;} // the middle node is not used; however their distance is larger than the limit;
}

/*  this function checks whether there is an edge between source and target;
*/
inline bool Contraction::checkDuplicate(VertexID source, VertexID target, Weight weight){
	bool isfound = false;

	forall_outneighbors(graph, source, eit){
		if (eit->target == target) {
			eit->weight = weight;
			isfound = true;
			break;
		}
	}
	for (int i = 0; i < tmpoutshortcutlist[source].size(); i++) {
		if (tmpoutshortcutlist[source][i].target == target){
			tmpoutshortcutlist[source][i].weight = weight;
			isfound = true;
			break;
		}
	}

	if (!isfound) return isfound;

	forall_outneighbors(graph, target, eit){
		if (eit->target == source) {
			eit->weight = weight;
			break;
		}
	}
	for (int i = 0; i < tmpoutshortcutlist[target].size(); i++) {
		if (tmpoutshortcutlist[target][i].target == source){
			tmpoutshortcutlist[target][i].weight = weight;
			break;
		}
	}

	return isfound;
}

/*
	add shortcuts;
*/
void Contraction::packShortcut(int nid,\
				vector<VertexID>& neighbors,\
				vector<double>& distance,\
				vector<vector<VertexID> >& inner_IDs,\
				vector<int>& matched_pair,\
				vector<double>& pair_distance){

	#ifdef CONTRACTION_DEBUG
		cout << "adding shortcuts..." << endl;
	#endif
	assert(matched_pair.size() % 2 == 0);
	for ( int i = 0; i < matched_pair.size(); i+=2 ) {
		int source = neighbors[matched_pair[i]];
		int target = neighbors[matched_pair[i+1]];
		#ifdef CONTRACTION_DEBUG
			cout << "source = " << source << ", target = " << target << endl;
		#endif
		Weight weight = pair_distance[i/2];

		bool isfound = checkDuplicate(source, target, weight); // if found, update its weight;
		if (isfound) {
			#ifdef CONTRACTION_DEBUG
				cout << "duplicate found!" << endl;
			#endif
			continue;
		}else {
			#ifdef CONTRACTION_DEBUG
				cout << "no duplicates!" << endl;
			#endif
			addShortCut(source, target, nid, weight, inner_IDs[matched_pair[i]], inner_IDs[matched_pair[i+1]]);
		}
	}

	removeEdge(nid);
}

/*
	add short cut edges;
*/
void Contraction::addShortCut(VertexID source,\
							  VertexID target,\
							  VertexID nid,\
							  Weight weight,\
							  vector<VertexID>& inner_ids_source,\
							  vector<VertexID>& inner_ids_target){
	if (nodelist[source].rank < UINT_INFINITY) cout << "Wrong-ranking" << endl;
	if (nodelist[target].rank < UINT_INFINITY) cout << "Wrong-ranking" << endl;


	vector<VertexID> source_to_target;
	for (int i = inner_ids_source.size()-1; i >= 0; i--) {
		source_to_target.push_back( inner_ids_source[i] );
	}
	source_to_target.push_back(nid);
	for ( int i = 0; i < inner_ids_target.size(); i++ ) {
		source_to_target.push_back( inner_ids_target[i] );
	}

	vector<VertexID> target_to_source;
	for ( int i = source_to_target.size()-1; i >= 0; i-- ) {
		target_to_source.push_back(source_to_target[i]);
	}

	ShortCutEdge shortcut_edge;
	shortcut_edge.target = source;
	shortcut_edge.flaga = false;
	shortcut_edge.flagb = false;
	shortcut_edge.innerIDs = source_to_target;
	shortcut_edge.weight = weight;

	#ifdef CONTRACTION_DEBUG
		cout<< "print out the edges in the path:";
		cout << source;
		for ( int i = 0; i < source_to_target.size(); i++ ) {
			cout << "->" << source_to_target[i];
		}
		cout << "->" << target << endl;
	#endif

	//tmpinshortcutlist[target].push_back( shortcut_edge );
	shortcut_edge.innerIDs = target_to_source;
	tmpoutshortcutlist[target].push_back( shortcut_edge );

	shortcut_edge.target = target;
	//tmpinshortcutlist[source].push_back( shortcut_edge );
	shortcut_edge.innerIDs = source_to_target;
	tmpoutshortcutlist[source].push_back( shortcut_edge );
}

/* remove edges around the node;
*/
void Contraction::removeEdge(int nid) {
	// remove the edges in original graphs;
	forall_outneighbors(graph, nid, eit){
		eit->flaga = true;
		forall_outneighbors(graph, eit->target, it){
			if (it->target == nid){
				it->flaga = true;
				break;
			}
		}
	}
	// forall_inneighbors(graph, nid, eit){
		// eit->flaga = true;
		// forall_outneighbors(graph, eit->target, it){
			// if ( it->target == nid ) {
				// it->flaga = true;
				// break;
			// }
		// }
	// }

	// remove the shortcuts;
	// for ( int i = 0; i < tmpinshortcutlist[nid].size(); i++ ) {
		// tmpinshortcutlist[nid][i].flaga = true;
		// int vid = tmpinshortcutlist[nid][i].target;
		// for ( int j = 0; j < tmpoutshortcutlist[vid].size(); j++ ){
			// if(tmpoutshortcutlist[vid][j].target == nid) {
				// tmpoutshortcutlist[vid][j].flaga = true;
				// break;
			// }
		// }
	// }
	for (int i = 0; i < tmpoutshortcutlist[nid].size(); i++) {
		tmpoutshortcutlist[nid][i].flaga = true;
		int vid = tmpoutshortcutlist[nid][i].target;
		for (int j = 0; j < tmpoutshortcutlist[vid].size(); j++){
			if (tmpoutshortcutlist[vid][j].target == nid) {
				tmpoutshortcutlist[vid][j].flaga = true;
				break;
			}
		}
	}
}

/*  this function computes the edge difference for vertex nid;
*/
int Contraction::computeEdgeDifference(int nid,\
									vector<VertexID>& neighbors,\
									vector<double>& distance,\
									vector<vector<VertexID> >& inner_IDs,\
									vector<int>& matched_pair,\
									vector<double>& pair_distance){
	vector<VertexID>(0).swap(neighbors);
	vector<double>(0).swap(distance);
	vector<vector<VertexID> >(0).swap(inner_IDs);
	vector<int>(0).swap(matched_pair);
	vector<double>(0).swap(pair_distance);

	// deal with neighbors in original graph;
	forall_outneighbors(graph, nid, eit){
		if (eit->flaga) continue; // removed edges;
		else {
			if (nodelist[eit->target].rank<UINT_INFINITY) cout << "Should not be dealed with!" << endl;
			neighbors.push_back(eit->target);
			distance.push_back(eit->weight);
			// inner_IDs.push_back(vector<EdgeID>(1, graph.get_outedgeindex(nid, eit->target)));
			inner_IDs.push_back(vector<VertexID>(0));
		}
	}

	// deal with neighbors linked with shortcuts;
	// cout << graph.num_vertices() << endl;
	// cout << tmpoutshortcutlist.size() << endl;
	// here we deal with undirected graph;
	for (int i = 0; i < tmpoutshortcutlist[nid].size(); i++){
		if (tmpoutshortcutlist[nid][i].flaga) continue;
		else {
			if (nodelist[tmpoutshortcutlist[nid][i].target].rank<UINT_INFINITY) cout << "Should not be dealed with!" << endl;
			neighbors.push_back(tmpoutshortcutlist[nid][i].target);
			distance.push_back(tmpoutshortcutlist[nid][i].weight);
			inner_IDs.push_back(tmpoutshortcutlist[nid][i].innerIDs);
		}
	}

	// for (int i = 0; i < neighbors.size(); i++) cout << neighbors[i] << " ";
	// cout << endl;

	// calculate edge difference pairwisely;
	// this part could be expensive when average degree is high;
	int edgedifference = 0;
	for (int i = 0; i < neighbors.size(); ++i) {
		for ( int j = i+1; j < neighbors.size(); ++j) {
			if (neighbors[i] == neighbors[j]) continue;
			double real_distance = DBL_INFINITY;
			double limit = distance[i]+distance[j];
			bool isfound = dodgingDijkstra(neighbors[i], neighbors[j], nid, limit, real_distance);

			#ifdef CONTRACTION_DEBUG
				// check the correctness of dodgingDijkstra function;
				cout << "from " << neighbors[i] << " to " << neighbors[j] << " through " << nid << endl;
				double dist = runDijkstra(neighbors[i], neighbors[j]);
				cout << "distance: " << dist << ", limit: " << limit << endl;;
				if (!isfound && dist < min(limit, real_distance)) {
					cout << "Wrong!" << endl;
				}
				cout << "~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
			#endif

			if (!isfound) {
				edgedifference++;
				matched_pair.push_back(i);
				matched_pair.push_back(j);
				pair_distance.push_back(min(real_distance, limit));
			}
		}
	}


	//cout << "[added shortcuts] for " << nid << ": " << edgedifference << endl;
	//cout << "[removed edges]   for " << nid << ": " << neighbors.size() << endl;
	edgedifference = edgedifference - neighbors.size();
	//cout << "[edge difference] for " << nid << ": " << edgedifference << endl;

	return edgedifference;
}


void Contraction::computeShortcuts(bool gp){
	vector<VertexID> neighbors;
	vector<double> distance;
	vector<vector<VertexID> > inner_IDs;
	vector<int> matched_pair;
	vector<double> pair_distance;

	// initial the priority_queue;
	cout << "Calculate the initial queue ...";
	priority_queue<pair<int, VertexID>, vector<pair<int, VertexID> >, myQueueComp> queue;
	forall_nodes(graph, nid) {
		int edgedifference = computeEdgeDifference(nid, neighbors, distance, inner_IDs, matched_pair, pair_distance);
		queue.push(make_pair(edgedifference, nid));
	}
	cout << "Done." << endl;

	// rank the nodes;
	// cout << "Running Queue:~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;

	int current_rank = 0;
	while (queue.size()>1) {
		int edge_difference = queue.top().first;
		VertexID nid = queue.top().second;
		cout << "Queue.size = " << queue.size() << endl;
		//cout << "vid: " << nid << endl;
		queue.pop();
		#ifdef CONTRACTION_DEBUG
			cout << nid << ": edge_difference = " << edge_difference << endl;
		#endif

		int current_edge_difference = computeEdgeDifference(nid, neighbors, distance, inner_IDs, matched_pair, pair_distance);
		if (current_edge_difference <= queue.top().first) { // lazy update;
			#ifdef CONTRACTION_DEBUG
				cout << "[result]: add short cut." << endl;
			#endif
			packShortcut(nid, neighbors, distance, inner_IDs, matched_pair, pair_distance);
			nodelist[nid].rank = current_rank++;
		}else{
			#ifdef CONTRACTION_DEBUG
				cout << "[result]: insert it again" << endl;
			#endif
			queue.push(make_pair(current_edge_difference, nid));
		}
	}
	if (nodelist[queue.top().second].rank == UINT_INFINITY) nodelist[queue.top().second].rank = current_rank;

	#ifdef CONTRACTION_DEBUG
		forall_nodes(graph, nid){
			cout << "rank " << nid << ":" << nodelist[nid].rank << endl;
		}
	#endif


	// insert shortcut into graph;
	#ifdef CONTRACTION_DEBUG
		// cout << "inShortCut:" << endl;
		// for ( int i = 0; i < tmpinshortcutlist.size(); i++ ) {
			// cout << i << ":";
			// for ( int j = 0; j < tmpinshortcutlist[i].size(); j++ ) {
				// cout << tmpinshortcutlist[i][j].target << " ";
			// }
			// cout << endl;
		// }
		cout << "outShortCut:" << endl;
		for ( int i = 0; i < tmpoutshortcutlist.size(); i++ ) {
			cout << i << ":";
			for ( int j = 0; j < tmpoutshortcutlist[i].size(); j++ ) {
				cout << tmpoutshortcutlist[i][j].target << " ";
			}
			cout << endl;
		}
	#endif

	#ifdef CONTRACTION_DEBUG
		checkShortCutCorrectness();
		checkShortCutCompleteness(2000);
	#endif

	// two things are produced: the shortcut edges, and the node rank;
	cout << "insert shortcut ..." << endl;
	graph.insertShortcut(tmpinshortcutlist, tmpoutshortcutlist);
	cout << "done." << endl;
	for ( int i = 0; i < nodesize; i++ ){
		graph[i].rank = nodelist[i].rank;
		graph[i].id = i;
	}

	//cout << "Insert short cut!" << endl;
	//graph.printShortcut();
	// print rank;
	// for ( int i = 0; i < nodelist.size(); i++ ) {
		// cout << i << "::" << nodelist[i].rank << endl;
	// }

	// create labels;
	cout << "insert label ..." << endl;
	generateLabel(gp);
	cout << "done." << endl;

	// check the correctness of labels;
	#ifdef CONTRACTION_DEBUG
		checkLabelCorrectness();
	#endif

	// check the correctness of paths;
	#ifdef CONTRACTION_DEBUG
		checkPathCorrectness();
	#endif
	//exit(0);
	//cout << "Insert label!" << endl;
}

void Contraction::checkShortCutCorrectness(){
	int count = 0;
	int error = 0;
	for (int i = 0; i < tmpoutshortcutlist.size(); i++) {
		for (int j = 0; j < tmpoutshortcutlist[i].size(); j++) {
			count++;
			VertexID u = i;
			VertexID v = tmpoutshortcutlist[i][j].target;
			// print the inner vertex ids;
			// cout << "Deal with vertex pair [" << u << ", " << v << "]:" << endl;
			// cout << "path: ";
			// cout << u << " ";
			// for (int k = 0; k < tmpoutshortcutlist[i][j].innerIDs.size(); k++) {
				// cout << tmpoutshortcutlist[i][j].innerIDs[k] << " ";
			// }
			// cout << v << endl;
			// cout << "distance: " << tmpoutshortcutlist[i][j].weight << endl;

			// the innerIDs should has lower ranks than the two ends;
			for (int k = 0; k < tmpoutshortcutlist[i][j].innerIDs.size(); k++) {
				if (nodelist[tmpoutshortcutlist[i][j].innerIDs[k]].rank > nodelist[u].rank ||\
					nodelist[tmpoutshortcutlist[i][j].innerIDs[k]].rank > nodelist[v].rank) {
					cout << "Wrong: rank!" << endl;
					cout << tmpoutshortcutlist[i][j].innerIDs[k] << ": " << nodelist[tmpoutshortcutlist[i][j].innerIDs[k]].rank << endl;
					cout << u << ": " << nodelist[u].rank << endl;
					cout << v << ": " << nodelist[v].rank << endl;
 				}
			}

			// the shortcut distance between u and v should equal the sum of edge distances from u to v.
			// cout << "Results from dijkstra search:" << endl;
			double dist = runDijkstra(u, v);
			// cout << "distance: " << dist << endl;

			cout << "Comparison distance: ";
			if (fabs(dist - tmpoutshortcutlist[i][j].weight) <= 1e-5) cout << "Right!" << endl;
			else { cout << "Wrong: distance!" << endl; error++;
				cout << "Dijkstra: " << dist << endl;
				cout << "Shortcut: " << tmpoutshortcutlist[i][j].weight << endl;
			}
 		}
	}
	cout << "Total shortcut#: " << count << ", wrong#: " << error << endl;
}

double Contraction::runDijkstra(VertexID source, VertexID target){
	//cout << "Test: source = " << source << ", target = " << target << "." << endl;
	vector<double> dist (nodesize);
	for (int vid = 0; vid < nodesize; vid++) {
		dist[vid] = DBL_INFINITY;
	}
	dist[source] = 0.0;

	vector<VertexID> father(nodesize, -1);

	priority_queue<pair<double, VertexID>, vector<pair<double, VertexID> >, QueueComp> Queue;
	Queue.push(make_pair(0.0, source));

	for ( int vid = 0; vid < nodesize; vid++) {
		nodelist[vid].flaga = false;
	}

	int num_visited = 0;
	while (Queue.size() > 0) {
		double min_dist = Queue.top().first;
		VertexID vid = Queue.top().second;
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
				father[eit->target] = vid;
			}
		}
	}


	// print path;
	vector<int> path;
	path.push_back(target);
	while(true){
		int c = path.back();
		if (father[c] != -1) path.push_back(father[c]);
		else break;
	}
	cout << "path without shortcut: ";
	for ( int i = path.size()-1; i >=0; i-- ) {
		cout << path[i] << " ";
	}
	cout << endl;
	cout << "corresponding distance: ";
	for (int i = path.size()-1; i>=0; i--) {
		cout << dist[path[i]] << " ";
	}
	cout << endl;

	forall_nodes(graph, vid){
		nodelist[vid].flaga = false;
	}

	return dist[target];
}

void Contraction::checkShortCutCompleteness(int sampleSize = 1000){
	int count = 0;
	for (int sample = 0; sample < sampleSize; sample++) {
		VertexID source;
		VertexID target;

		while (true) {
			source = rand() % nodesize;
			target = rand() % nodesize;
			if (source != target) break;
		}

		double result_without_shortcut = runDijkstra(source, target);

		vector<VertexID> fromSourceNode, fromTargetNode;
		vector<Weight> fromSourceDist, fromTargetDist;
		runShortcutDijkstra(source, fromSourceNode, fromSourceDist);
		runShortcutDijkstra(target, fromTargetNode, fromTargetDist);

		double dist = DBL_INFINITY;
		int meet_node = -1;
		int i = 0, j = 0;
		while( i < fromSourceNode.size() && j < fromTargetNode.size() ) {
			if(fromSourceNode[i] == fromTargetNode[j]) {
				if (fromSourceDist[i] + fromTargetDist[j] < dist ) {
					dist = fromSourceDist[i] + fromTargetDist[j];
					meet_node = fromSourceNode[i];
				}
				i++;
				j++;
			}else if (fromSourceNode[i] == target) {
				if (fromSourceDist[i] < dist) {
					dist = fromSourceDist[i];
					meet_node = target;
				}
				i++;
			}else if (fromTargetNode[j] == source) {
				if (fromTargetDist[j] < dist) {
					dist = fromTargetDist[j];
					meet_node = source;
				}
				j++;
			}else if (fromSourceNode[i] < fromTargetNode[j]) i++;
			else j++;
		}

		cout << "Without shortcut: " << result_without_shortcut << endl;
		cout << "With shortcut: " << dist << ", middle point: " << meet_node << endl;
		if (fabs(result_without_shortcut - dist) > 1e-5) {
			cout << "Wrong!" << endl;
			runDijkstra(source, meet_node);
			runDijkstra(target, meet_node);
			count++;
			exit(0);
		}
	}

	cout << "Test size: " << sampleSize << ", wrong: " << count << endl;
}

void Contraction::runShortcutDijkstra(VertexID source, vector<VertexID>& nodeList, vector<Weight>& distanceList){
	//cout << "Test: source = " << source << ", target = " << target << "." << endl;
	vector<double> dist (nodesize);
	for ( int vid = 0; vid < nodesize; vid++) {
		dist[vid] = DBL_INFINITY;
	}
	dist[source] = 0.0;

	vector<VertexID> father(nodesize, -1);

	priority_queue<pair<double, VertexID>, vector<pair<double, VertexID> >, QueueComp> Queue;
	Queue.push(make_pair(0.0, source));

	for ( int vid = 0; vid < nodesize; vid++) {
		nodelist[vid].flaga = false;
	}

	int num_visited = 0;
	while (Queue.size() > 0) {
		double min_dist = Queue.top().first;
		VertexID vid = Queue.top().second;
		//cout << "vid: " << vid << ", min_dist: " << min_dist << endl;
		Queue.pop();
		if (min_dist > dist[vid]) {/*cout << "Greater!?" << endl;*/ continue;} // lazy update;
		else nodelist[vid].flaga = true; // settle vertex vid;

		num_visited++;

		forall_outneighbors(graph, vid, eit){
			if (nodelist[eit->target].flaga) continue; // settled;
			if (nodelist[eit->target].rank < nodelist[vid].rank /*&& vid != source*/) continue; // using rank filter;
			if (dist[eit->target] > min_dist + eit->weight) {
				dist[eit->target] = min_dist + eit->weight;
				Queue.push(make_pair(dist[eit->target], eit->target));
				father[eit->target] = vid;
			}
		}

		for (int i = 0; i < tmpoutshortcutlist[vid].size(); i++) {
			VertexID tmp_t = tmpoutshortcutlist[vid][i].target;
			if (nodelist[tmp_t].flaga) continue;
			if (nodelist[tmp_t].rank < nodelist[vid].rank && vid != source) continue;
			double found = false;
			for (int j = 0; j < tmpoutshortcutlist[vid][i].innerIDs.size(); j++){
				VertexID iv = tmpoutshortcutlist[vid][i].innerIDs[j];
				//if (dist[iv] >= dist[tmp_t]) {found = true; break;}
				if (dist[iv] <= dist[vid]) {found = true; break;}
			}
			if (found) continue;
			if (dist[tmp_t] > min_dist + tmpoutshortcutlist[vid][i].weight) {
				dist[tmp_t] = min_dist + tmpoutshortcutlist[vid][i].weight;
				Queue.push(make_pair(dist[tmp_t], tmp_t));
				father[tmp_t] = vid;
			}
		}
	}


	for (int i = 0; i < nodesize; i++) {
		if (dist[i] < DBL_INFINITY) {
			nodeList.push_back(i);
			distanceList.push_back(dist[i]);
		}
	}

	for (int i = 0; i < nodeList.size(); i++) {
		cout << "[" << nodeList[i] << ", " << distanceList[i] << "] ";
	}
	cout << endl;

	for ( int vid = 0; vid < nodesize; vid++) {
		nodelist[vid].flaga = false;
	}
}

void Contraction::storePathInfo(int source, int target, vector<VertexID>& father, vector<vector<VertexID> >& pathInfo) {
	if (pathInfo[target].size() != 0) return;
	cout << "source: " << source << ", target: " << target << endl;

	if (father[source] != -1 && pathInfo[source].size() == 0) {
		int f = father[source];
		storePathInfo(f, source, father, pathInfo);
	}

	// cout << "source: " << source << ", target: " <<  target << endl;
	// copy father's path;
	for ( int i = 0; i < pathInfo[source].size(); i++ ) {
		// cout << pathInfo[source][i] << " ";
		pathInfo[target].push_back(pathInfo[source][i]);
	}
	// cout << endl;

	// check the original graph;
	// if they are linked by the simple edge, only
	// target node is added to the path;
	bool isfound = false;
	forall_outneighbors(graph, source, eit) {
		if (eit->target == target) {
			// cout << "direct!" << endl;
			isfound = true;
			pathInfo[target].push_back(target);
			break;

			cout << "target (original): " << target << ": ";
			for (int i = 0; i < pathInfo[target].size(); i++) {
				cout << pathInfo[target][i] << " ";
			}
			cout << endl;
		}
	}
	if (isfound) return;
	// check the short cut graph;
	// if they are linked by the shortcut,
	// the inner vertex ids and the target node are added to the path;
	forall_outshortcuts(graph, source, eit){
		if (eit->target == target) {
			// cout << "cut!" << endl;
			for(int i = 0; i < eit->innerIDs.size(); i++){
				pathInfo[target].push_back(eit->innerIDs[i]);
			}
			pathInfo[target].push_back(target);

			cout << "target (shortcut): " << target << ": ";
			for (int i = 0; i < pathInfo[target].size(); i++) {
				cout << pathInfo[target][i] << " ";
			}
			cout << endl;

			return;
		}
	}
}

/* @parameter:
		father: record the predecesor of each node on the short path;
		index_list: record the shortest path end points for each node on the path;
		path_stack: record the current shortest path;
		new_path: when function returns, true; otherwise, false;
*/
void Contraction::buildPathIndex(VertexID source, vector<VertexID>& father,\
								vector<vector<int> >& index_list,\
								vector<VertexID>& path_stack,\
								bool& new_path) {
	// cout << "deal with " << path_stack.back() << endl;
	// cout << "path length " << path_stack.size() << endl;
	// dfs from source to find pathes;
	// a new node is appended;
	VertexID vid = path_stack.back();
	forall_outneighbors(graph, vid, eit) {
		VertexID uid = eit->target;
		if (father[uid] != vid) continue;

		// cout << "child: " << uid << endl;

		if (!new_path) {
			index_list[0][uid] = index_list[0][vid];
			index_list[1][uid] = index_list[1][vid]+1;
			path_stack.push_back(uid);
			pathList[source].push_back(uid);
		}else {
			new_path = false;
			path_stack.push_back(uid);
			index_list[0][uid] = pathList[source].size();
			for (int i = index_list[0][vid]; i < index_list[1][vid]; i++) {
				pathList[source].push_back(pathList[source][i]);
			}
			pathList[source].push_back(uid);
			index_list[1][uid] = pathList[source].size();
		}

		buildPathIndex(source, father, index_list, path_stack, new_path);
		path_stack.pop_back();
	}

	forall_outshortcuts(graph, vid, eit) {
		VertexID uid = eit->target;
		if (father[uid] != vid) continue;

		// cout << "short child: " << uid << endl;

		if (!new_path) {
			// cout << "old?" << endl;
			index_list[0][uid] = index_list[0][vid];
			index_list[1][uid] = index_list[1][vid] + eit->innerIDs.size()+1;
			for (int i = 0; i < eit->innerIDs.size(); i++) {
				pathList[source].push_back(eit->innerIDs[i]);
			}
			path_stack.push_back(uid);
			pathList[source].push_back(uid);
		}else {
			// cout << "b " << index_list[0][vid] << ", e " << index_list[1][vid] << endl;
			index_list[0][uid] = pathList[source].size();
			// cout << "new!" << endl;
			// cout << "vid " << vid << endl;
			new_path = false;
			path_stack.push_back(uid);
			index_list[0][uid] = pathList[source].size();
			for (int i = index_list[0][vid]; i < index_list[1][vid]; i++) {
				pathList[source].push_back(pathList[source][i]);
			}
			for (int i = 0; i < eit->innerIDs.size(); i++) {
				pathList[source].push_back(eit->innerIDs[i]);
			}
			pathList[source].push_back(uid);
			index_list[1][uid] = pathList[source].size();
		}

		buildPathIndex(source, father, index_list, path_stack, new_path);
		path_stack.pop_back();
	}

	new_path = true;
}

// create label for each vertex by Dijkstra algorithm on extended graphs;
void Contraction::createLabel(VertexID source, bool gp){
	//cout << "Test: source = " << source << ", target = " << target << "." << endl;
	vector<double> dist (nodesize);
	for (int vid = 0; vid < nodesize; vid++) dist[vid] = DBL_INFINITY;
	dist[source] = 0.0;

	vector<VertexID> father(nodesize, numeric_limits<VertexID>::max());

	priority_queue<pair<double, VertexID>, vector<pair<double, VertexID> >, QueueComp> Queue;
	Queue.push(make_pair(0.0, source));

	for (int vid = 0; vid < nodesize; vid++) nodelist[vid].flaga = false;

	while (Queue.size() > 0) {
		double min_dist = Queue.top().first;
		VertexID vid = Queue.top().second;
		//cout << "vid = " << vid << ", min_dist = " << min_dist << endl;

		Queue.pop();
		if (min_dist > dist[vid]) {/*cout << "Greater!?" << endl;*/ continue;} // lazy update;
		else nodelist[vid].flaga = true; // settle vertex vid;

		forall_outneighbors(graph, vid, eit){
			if (nodelist[eit->target].flaga) continue; // settled;
			if (/*vid != source &&*/ nodelist[eit->target].rank < nodelist[vid].rank) continue;
			if (dist[eit->target] > min_dist + eit->weight) {
				dist[eit->target] = min_dist + eit->weight;
				Queue.push(make_pair(dist[eit->target], eit->target));
				father[eit->target] = vid;
			}
		}

		forall_outshortcuts(graph, vid, eit){
			if (nodelist[eit->target].flaga) continue; // settled;
			if (/*vid != source &&*/ nodelist[eit->target].rank < nodelist[vid].rank) continue;
			double found = false;
			for (int j = 0; j < eit->innerIDs.size(); j++){
				VertexID iv = eit->innerIDs[j];
				// if (dist[iv] >= dist[eit->target]) {found = true; break;}
				if (dist[iv] <= dist[vid]) {found = true; break;}
			}
			if (found) continue;
			if (dist[eit->target] > min_dist + eit->weight) {
				dist[eit->target] = min_dist + eit->weight;
				Queue.push(make_pair(dist[eit->target], eit->target));
				father[eit->target] = vid;
			}
		}
	}

	if (/*inlabellist.size() < nodesize ||*/ outlabellist.size() < nodesize) {
		//inlabellist.resize(nodesize);
		outlabellist.resize(nodesize);
	}

	// add in inlabellist and outlabellist;
	for (int vid = 0; vid < nodesize; vid++) {
		if (dist[vid] < DBL_INFINITY) {
			Label tmp_label;
			tmp_label.id = vid;
			tmp_label.distance = dist[vid];
			outlabellist[source].push_back(tmp_label);
			//tmp_label.id = source;
			//inlabellist[vid].push_back(tmp_label);
		}
	}

	// for test
	// cout << "For " << source << endl;
	// for (int i = 0; i < outlabellist[source].size(); i++){
		// cout << outlabellist[source][i].id << " " << outlabellist[source][i].distance << endl;
	// }
	// exit(0);

	// for debug,
	/*for (int i = 0; i < outlabellist[source].size(); i++){
		int vid = outlabellist[source][i].id;
		vector<int> path;
		path.push_back(vid);
		while(true){
			int c = path.back();
			if (father[c] != numeric_limits<VertexID>::max()) path.push_back(father[c]);
			else break;
		}
		cout << "path with shortcut: ";
		for ( int i = path.size()-1; i >=0; i-- ) {
			cout << path[i] << " ";
		}
		cout << endl;
		cout << "corresponding distance: ";
		for (int i = path.size()-1; i>=0; i--) {
			cout << dist[path[i]] << " ";
		}
		cout << endl;
	}*/

	/***********************************************************************/
	// for debug, check the correctness of labels;
	/*
	for (int i = 0; i < outlabellist[source].size(); i++){
		int vid = outlabellist[source][i].id;
		double d = outlabellist[source][i].distance;

		cout << "D W: " << source << " vs " << vid << endl;

		double r_dijkstra = runDijkstra(source, vid);

		if (fabs(r_dijkstra - d) >= 1e-5) {
			cerr << "WRONG: " << endl;
			cerr << source << "-->" << vid << ": " << d << " vs " << r_dijkstra << endl;

			// if wrong, compare their paths;
			vector<int> path;
			path.push_back(vid);
			while(true){
				int c = path.back();
				if (father[c] != -1) path.push_back(father[c]);
				else break;
			}
			cout << "path with shortcut: ";
			for ( int i = path.size()-1; i >=0; i-- ) {
				cout << path[i] << " ";
			}
			cout << endl;
			cout << "corresponding distance: ";
			for (int i = path.size()-1; i>=0; i--) {
				cout << dist[path[i]] << " ";
			}
			cout << endl;
		}
	}*/
	/***********************************************************************/

	if (!gp) return;

	// create path information for unpacking path;
	vector<vector<int> > index_list(2, vector<int>(nodesize, 0));
	vector<VertexID> path_stack;
	bool new_path = false;
	path_stack.push_back(source);
	index_list[0][source] = 0; // beginning of the path;
	index_list[1][source] = 1; // end of the path + 1;
	pathList[source].push_back(source);
	buildPathIndex(source, father, index_list, path_stack, new_path);

	// create path index for each node;
	for (int vid = 0; vid < nodesize; vid++) {
		if (dist[vid] < DBL_INFINITY) {
			pathIndex[source].push_back(index_list[0][vid]);
			pathIndex[source].push_back(index_list[1][vid]);
		}
	}

	#ifdef CONTRACTION_DEBUG
		cout << "source: " << source << endl;
		cout << "father: " << endl;
		for (int i = 0; i < father.size(); i++) {
			cout << "[" << i << " " << father[i] << "] ";
		}
		cout << endl;

		cout << "index_list: " << endl;
		for (int i = 0; i < index_list[0].size(); i++) {
			cout << "[" << i << " " << index_list[0][i] << " " << index_list[1][i] << "] ";
		}
		cout << endl;
		cout << "dist: " << endl;
		for (int i = 0; i < dist.size(); i++) {
			cout << "[" << i << " " << dist[i] << "] ";
		}
		cout << endl;
	#endif
}


void Contraction::checkLabelCorrectness(void){
	int count = 0;
	for (int source = 0; source < nodesize; source++) {
		for (int target = source+1; target < nodesize; target++) {

			double result_without_shortcut = runDijkstra(source, target); // using dijkstra algorithm;

			// using label information;
			double dist = DBL_INFINITY;
			int meet_node = -1;
			int i = 0, j = 0;
			while( i < outlabellist[source].size() && j < outlabellist[target].size() ) {
				if(outlabellist[source][i].id == outlabellist[target][j].id) {
					if (outlabellist[source][i].distance + outlabellist[target][j].distance < dist ) {
						dist = outlabellist[source][i].distance + outlabellist[target][j].distance;
						meet_node = outlabellist[source][i].id;
					}
					i++;
					j++;
				}else if (outlabellist[source][i].id == target) {
					if (outlabellist[source][i].distance < dist) {
						dist = outlabellist[source][i].distance;
						meet_node = target;
					}
					i++;
				}else if (outlabellist[target][j].id == source) {
					if (outlabellist[target][j].distance < dist) {
						dist = outlabellist[target][j].distance;
						meet_node = source;
					}
					j++;
				}else if (outlabellist[source][i].id < outlabellist[target][j].id) i++;
				else j++;
			}

			cout << "source: " << source << ", target: " << target << endl;
			cout << "Without shortcut: " << result_without_shortcut << endl;
			cout << "With shortcut: " << dist << ", middle point: " << meet_node << endl;
			if (fabs(result_without_shortcut - dist) > 1e-5) {
				cout << "Wrong!" << endl;
				runDijkstra(source, meet_node);
				runDijkstra(target, meet_node);
				count++;
			}
		}
	}
}

void Contraction::checkPathCorrectness(void){
	int count = 0;
	for (int source = 0; source < nodesize; source++) {
		cout << source << ": #path: " << pathIndex[source].size() << endl;
		for (int j = 0; j < pathIndex[source].size(); j+=2) {

			int begin = pathIndex[source][j];
			int end = pathIndex[source][j+1];

			// print
			cout << "source: " << source << ", target = " << outlabellist[source][j/2].id;
			cout << ", weight = " << outlabellist[source][j/2].distance << endl;
			cout << "begin = " << begin << ", end = " << end << endl;

			cout << "path: ";
			for (int k = begin; k < end; k++) {
				cout << pathList[source][k] << " ";
			}
			cout << endl;

			/* two things to check:
				1. no shortcut;
				2. distance matches the label distance;
			*/
			double dist = 0;
			cout << source << " ";
			for (int k = begin; k < end-1; k++) {
				int uid = pathList[source][k];
				int vid = pathList[source][k+1];

				bool found = false;
				forall_outneighbors(graph, uid, eit) {
					if (eit->target == vid) {
						found = true;
						dist += eit->weight;
						cout << "[" << eit->target << ", " << dist << "] ";
						break;
					}
				}

				if (!found) {
					cout << "the edge can not be found: " << uid << " -> " << vid << endl;
					forall_outshortcuts(graph, uid, eit) {
						if (eit->target == vid) {
							cout << "find shortcut: " << uid << " -> " << vid << endl;
							cout << "inner ids: ";
							for (int i = 0; i < eit->innerIDs.size(); i++){
								cout << eit->innerIDs[i] << " ";
							}
							cout << endl;
						}
					}
				}
			}
			cout << endl;
			if (fabs(dist - outlabellist[source][j/2].distance) > 1e-5) {
				cout << "the distance does not match by using dijkstra: " << dist << endl;
				count++;
			}
		}
	}
	cout << "# wrong paths: " << count << endl;
}

void Contraction::generateLabel(bool gp) {
	for ( int i = 0; i < nodesize; i++ ) {
		cout << "generate label: " << i << endl;
		createLabel(i, gp);
	}

	// createLabel(36719, gp);
	// for (int i=0; i < outlabellist[36719].size(); i++){
		// cout << outlabellist[36719][i].id << " " << outlabellist[36719][i].distance << endl;
	// }

	// exit(0);

	//printLabel();
	//printPath();

	graph.insertShortcut(tmpinshortcutlist, tmpoutshortcutlist);
	graph.insertLabelList(inlabellist, outlabellist);
	//graph.exportPathIndex().swap(pathIndex);
	//graph.exportPathList().swap(pathList);
}

vector<ShortCuts>& Contraction::exportOutShortcut(void) {
	return tmpoutshortcutlist;
}

vector<ShortCuts>& Contraction::exportInShortcut(void) {
	return tmpinshortcutlist;
}

VertexList& Contraction::exportNodeList(void) {
	return nodelist;
}

void Contraction::printLabel(){
	cout << "out label list: " << endl;
	for ( int i = 0; i < outlabellist.size(); i++ ) {
		cout << i << ": ";
		for ( int j = 0; j < outlabellist[i].size(); j++ ) {
			cout << "[" << outlabellist[i][j].id << ", " << outlabellist[i][j].distance << "] ";
		}
		cout << endl;
	}
	// cout << "in label list: " << endl;
	// for ( int i = 0; i < inlabellist.size(); i++ ) {
		// cout << i << ": ";
		// for ( int j = 0; j < inlabellist[i].size(); j++ ) {
			// cout << "[" << inlabellist[i][j].id << ", " << inlabellist[i][j].distance << "] ";
		// }
		// cout << endl;
	// }
}

void Contraction::printPath(){
	cout << "print path information:" << endl;
	for ( int i = 0; i < pathIndex.size(); i++ ) {
		cout << i << ":";
		for ( int j = 0; j < pathIndex[i].size()-1; j++ ) {
			int begin = pathIndex[i][j];
			int end = pathIndex[i][j+1];
			cout << outlabellist[i][j].id << "[";
			for (; begin < end; begin++) {
				cout << pathList[i][begin] << " ";
			}
			cout << "] ";
		}
		cout << endl;
	}
}
