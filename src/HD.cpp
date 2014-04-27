#include "HD.h"

struct QueueComp {
	bool operator()( const pair<double, VertexID>& a, const pair<double, VertexID>& b ) const {
		return ( a.first >= b.first );
	}
};

HD::HD(Graph& g):graph(g){}

HD::~HD(){}

// create label for each vertex by Dijkstra algorithm on extended graphs;
void HD::findPath(VertexID source, VertexID target, vector<VertexID>& path, double& d){
	if (source == target) {
		path.push_back(source);
		return;
	}
	
	//cout << "Test: source = " << source << ", target = " << target << "." << endl; 
	int nodesize = graph.num_vertices();
	
	vector<double> dist (nodesize);
	for (int vid = 0; vid < nodesize; vid++) dist[vid] = DBL_INFINITY;
	dist[source] = 0.0;
	
	vector<VertexID> father(nodesize, numeric_limits<unsigned int>::max());
	
	priority_queue<pair<double, VertexID>, vector<pair<double, VertexID> >, QueueComp> Queue;
	Queue.push(make_pair(0.0, source));
	
	VertexList nodelist(nodesize);	
	for (int vid = 0; vid < nodesize; vid++) {
		nodelist[vid].flaga = false;
		// nodelist[vid].rank = graph[vid].rank;
	}
	
	int i = 0;
	while (Queue.size() > 0) {
		//cout << i++ << endl;
		double min_dist = Queue.top().first;
		VertexID vid = Queue.top().second;
		//cout << "vid = " << vid << ", min_dist = " << min_dist << endl;

		Queue.pop();
		if (min_dist > dist[vid]) {/*cout << "Greater!?" << endl;*/ continue;} // lazy update; 
		else nodelist[vid].flaga = true; // settle vertex vid;
		
		if (vid == target) break; // find the target node;
		
		forall_outneighbors(graph, vid, eit){
			if (nodelist[eit->target].flaga) continue; // settled;
			// if (/*vid != source &&*/ nodelist[eit->target].rank < nodelist[vid].rank) continue;
			if (dist[eit->target] > min_dist + eit->weight) {
				dist[eit->target] = min_dist + eit->weight;
				Queue.push(make_pair(dist[eit->target], eit->target)); 
				father[eit->target] = vid;
			}	
		}
		
		// forall_outshortcuts(graph, vid, eit){
			// if (nodelist[eit->target].flaga) continue; // settled;
			// if (/*vid != source &&*/ nodelist[eit->target].rank < nodelist[vid].rank) continue;			
			// if (dist[eit->target] > min_dist + eit->weight) {
				// dist[eit->target] = min_dist + eit->weight;
				// Queue.push(make_pair(dist[eit->target], eit->target)); 
				// father[eit->target] = vid;
			// }	
		// }		
	}
	
	// create path information for unpacking path;
	vector<VertexID>(0).swap(path);
	path.push_back(target);
	int index = target;	
	while (index != source){
		if (father[index] != numeric_limits<unsigned int>::max()) {
			path.push_back(father[index]);
			index = father[index];
		}else break;
	}
	
	if (path.size()==1) path.push_back(source);
	
	// for test;
	for (int i = path.size()-1; i > 0; i--){
		int a = path[i];
		int b = path[i-1];
		forall_outneighbors(graph, a, eit){
			if (eit->target == b){
				d += eit->weight;
				break;
			}
		}
	}
	// cout << "d = " << d << endl;
	
	// cout << "father: " << endl;
	// for (int i = 0; i < father.size(); i++) {
		// cout << i << ": " << father[i] << endl;
	// }
	
	// cout << "~~~~~" << endl;
	// for (int i = path.size()-1; i >= 0; i--){
		// cout << path[i] << " ";
	// }
	// cout << endl;
}

double HD::run(VertexID source, VertexID target, VertexID& meeting_node) const{
	const LabelList& outlabel = graph.outLabelList[source];
	const LabelList& inlabel = graph.outLabelList[target];
	
	double dist = DBL_INFINITY;
	int meet_node = -1;
	int i = 0, j = 0;
	while(i < outlabel.size() && j < inlabel.size()) {
		if(outlabel[i].id == inlabel[j].id) {
			if (outlabel[i].distance + inlabel[j].distance < dist) {
				dist = outlabel[i].distance + inlabel[j].distance;
				meet_node = inlabel[j].id;
			}
			i++;
			j++;
		}else if (outlabel[i].id == target) {
			if (outlabel[i].distance < dist) {
				dist = outlabel[i].distance;
				meet_node = outlabel[i].id;
			}
			i++;
		}else if (inlabel[j].id == source) {
			if (inlabel[j].distance < dist) {
				dist = inlabel[j].distance;
				meet_node = inlabel[j].id;
			}
			j++;
		}else if (outlabel[i].id < inlabel[j].id) i++;
		else j++;
	}
	
	#ifdef HD_DEBUG
		cout << source << "->outlabellist: ";
		for ( int i = 0; i < outlabel.size(); i++ ) cout << "[" << outlabel[i].id << ", " << outlabel[i].distance << "] ";
		cout << endl;
		cout << target << "->inlabellist: ";	
		for ( int i = 0; i < outlabel.size(); i++ ) cout << "[" << inlabel[i].id << ", " << inlabel[i].distance << "] ";
		cout << endl;	
		cout << endl << "meet at node " << meet_node << endl;
	#endif
	
    /*cout << endl << "meet at node " << meet_node << endl;	
	if (meet_node == source) { cout << 0 << endl; cout << dist << endl; }
	else if (meet_node == target) {cout << dist << endl; cout << 0 << endl;}
	else {
		cout << outlabel[--i].distance << endl;
		cout << inlabel[--j].distance << endl;
	}*/
	meeting_node = meet_node;
	return dist;
}
