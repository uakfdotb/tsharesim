#include "Query.h"

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

Query::Query(Graph& g):graph(g){
	VertexList(0).swap(nodelist);
	for (int i = 0; i < graph.vertices().size(); i++) {
		nodelist.push_back(graph[i]);
	}
	nodesize = graph.num_vertices();
}
Query::~Query(){}

double Query::BiDijkstra(VertexID source, VertexID target, vector<EdgeID>& inneredge){
	priority_queue<pair<double, VertexID>, vector<pair<double, VertexID> >,  QueueComp> forward_queue, backward_queue;
	
	forward_queue.push(make_pair(0, source));
	backward_queue.push(make_pair(0, target));
	
	double min_distance = DBL_INFINITY;
	VertexID meeting_node = -1;
	
	vector<double> dist_forward(nodesize, 0);
	vector<double> dist_backward(nodesize, 0);
	for ( int i = 0; i < nodesize; i++ ) {
		dist_forward[i] = DBL_INFINITY;
		dist_backward[i] = DBL_INFINITY;
	}
	dist_forward[source] = 0;
	dist_backward[target] = 0;
	
	// cout << "Here???" << endl;
	
	// initialize graph;
	for (int vid = 0; vid < nodesize; vid++) {
		nodelist[vid].flaga = false; 
		nodelist[vid].flagb = false;
		nodelist[vid].flagc = false;
		nodelist[vid].flagd = false;
	}
	nodelist[source].flagc = true; // forward visit the node;
	nodelist[target].flagd = true; // backward visit the node;
	
	bool forward_finish = false;
	bool backward_finish = false;
	
	bool forward_insert = true;
	bool backward_insert = true;

	// cout << "there?" << endl;
	int num_visited = 0;
	while(true){
		if (forward_queue.size() == 0) forward_finish = true;
		if (backward_queue.size() == 0) backward_finish = true;
		
		// stop conditions;
		if (forward_finish && backward_finish) break;
		
		// forward one step;
		#ifdef DEBUG_MODE
			cout << "trace forward:" << endl;
		#endif
		while(forward_queue.size() > 0 && !forward_finish){
			double min_dist_forward = forward_queue.top().first;
			VertexID vid_forward = forward_queue.top().second;
			#ifdef DEBUG_MODE
				cout << "top: " << min_dist_forward << "|" << vid_forward << endl;
			#endif		
			
			// case 1: have been settled;
			forward_queue.pop();
			if (nodelist[vid_forward].flaga) continue;
			else nodelist[vid_forward].flaga = true; // settle the node;
			
			// case 2: if its single side distance greater than min dist;
			if (dist_forward[vid_forward] >= min_distance) {
				forward_finish = true;
				break;
			}
			
			/*// case 3: if its both sides distance sum greater than min dist; // something wrong with early termination;
			if(nodelist[vid_forward].flagb) {
				// cout << "mdis = " << min_distance << endl;
				if (min_distance > dist_forward[vid_forward] + dist_backward[vid_forward]) {
					min_distance = dist_forward[vid_forward] + dist_backward[vid_forward];
					meeting_node = vid_forward;
				}
				forward_insert = false;				
			}
			
			// case 4: if one side is fixed, at least the min dist might be reduced;
			else if (nodelist[vid_forward].flagd) {
				if (min_distance > dist_forward[vid_forward] + dist_backward[vid_forward]) {
					min_distance = dist_forward[vid_forward] + dist_backward[vid_forward];
					meeting_node = vid_forward;
				}
			}*/
			
			#ifdef DEBUG_MODE			
				num_visited++;			
			#endif
			
			#ifdef DEBUG_MODE
				cout << "check neighbors: ";
			#endif
			forall_outneighbors(graph, vid_forward, eit){
				// case 0:
				if (eit->target == source) continue;
				// case 1:
				if (eit->target == target) {
					if ( dist_forward[vid_forward] + eit->weight < min_distance ) {
						min_distance = dist_forward[vid_forward] + eit->weight;
						meeting_node = target;
					}
					continue;					
				}
				
				// case 2:
				#ifdef DEBUG_MODE		
					cout << eit->target;
				#endif			
				if (vid_forward != source && nodelist[eit->target].rank < nodelist[vid_forward].rank) { 
					#ifdef DEBUG_MODE	
						cout << "@" << " "; 
					#endif
					continue; 
				}
				#ifdef DEBUG_MODE	
					else{
						cout << " ";
					}
				#endif
				// case 3:
				if (dist_forward[eit->target] > dist_forward[vid_forward] + eit->weight){
					dist_forward[eit->target] = dist_forward[vid_forward] + eit->weight;
					if (dist_forward[eit->target] >= min_distance) continue;
					if (!nodelist[eit->target].flagd && !forward_insert) continue;
					forward_queue.push(make_pair(dist_forward[eit->target], eit->target));
					nodelist[eit->target].flagc = true;
					 // check meet: f->m<-b
					if (nodelist[eit->target].flagd && min_distance > dist_forward[eit->target]\
						+ dist_backward[eit->target]) {
						#ifdef DEBUG_MODE	
							cout << "update distance!" << endl;
						#endif
						min_distance = dist_forward[eit->target] + dist_backward[eit->target];
						meeting_node = eit->target;
					}
				}
			}
			
			forall_outshortcuts(graph, vid_forward, eit){
				if (eit->target == source) continue;			
				if (eit->target == target) {
					if (dist_forward[vid_forward] + eit->weight < min_distance) {
						min_distance = dist_forward[vid_forward] + eit->weight;
						meeting_node = target;
					}
					continue;					
				}			
			
				#ifdef DEBUG_MODE		
					cout << eit->target;
				#endif	
				if (vid_forward != source && nodelist[eit->target].rank < nodelist[vid_forward].rank) { 
					#ifdef DEBUG_MODE	
						cout << "@" << " "; 
					#endif
					continue; 
				}
				#ifdef DEBUG_MODE	
					else{
						cout << " ";
					}
				#endif
					
				if (dist_forward[eit->target] > dist_forward[vid_forward] + eit->weight){
					dist_forward[eit->target] = dist_forward[vid_forward] + eit->weight;
					if (dist_forward[eit->target] >= min_distance) continue;
					if (!nodelist[eit->target].flagd && !forward_insert) continue;
					forward_queue.push(make_pair(dist_forward[eit->target], eit->target));
					nodelist[eit->target].flagc = true;
					 // check meet: f->m<-b
					if (nodelist[eit->target].flagd && min_distance > dist_forward[eit->target]\
						+ dist_backward[eit->target]) {
						#ifdef DEBUG_MODE	
							cout << "update distance!" << endl;
						#endif
						min_distance = dist_forward[eit->target] + dist_backward[eit->target];
						meeting_node = eit->target;
					}					
				}
			}
			#ifdef DEBUG_MODE
				cout << endl;
			#endif	
					
			break;
		}
		
		
		// backward one step;
		#ifdef DEBUG_MODE
			cout << "trace backward:" << endl;
		#endif	
		while(backward_queue.size()>0 && !backward_finish){
			double min_dist_backward = backward_queue.top().first;
			VertexID vid_backward = backward_queue.top().second;
			#ifdef DEBUG_MODE
				cout << "top: " << min_dist_backward << "|" << vid_backward << endl;
			#endif				
			
			backward_queue.pop();
			if (nodelist[vid_backward].flagb) continue;
			else nodelist[vid_backward].flagb = true;
					
			
			if (dist_backward[vid_backward] >= min_distance ) {
				backward_finish = true;
				break;
			}
			/* something wrong with early termination;
			if (nodelist[vid_backward].flaga) {
				if (min_distance > dist_forward[vid_backward] + dist_backward[vid_backward]) {
					min_distance = dist_forward[vid_backward] + dist_backward[vid_backward];
					meeting_node = vid_backward;
				}
				backward_insert = false;
			}else if (nodelist[vid_backward].flagc) {		
				if (min_distance > dist_forward[vid_backward] + dist_backward[vid_backward]) {
					min_distance = dist_forward[vid_backward] + dist_backward[vid_backward];
					meeting_node = vid_backward;
				}
			}*/

			#ifdef DEBUG_MODE			
				num_visited++;			
			#endif				
			
			#ifdef DEBUG_MODE
				cout << "check neighbors: ";
			#endif			
			forall_outneighbors(graph, vid_backward, eit){
				if (eit->target == target) continue;		
				if (eit->target == source) {
					if ( dist_backward[vid_backward] + eit->weight < min_distance ) {
						min_distance = dist_backward[vid_backward] + eit->weight;
						meeting_node = source;
					}
					continue;	
				}			
			
				#ifdef DEBUG_MODE		
					cout << eit->target;
				#endif			
				if (vid_backward != target && nodelist[eit->target].rank < nodelist[vid_backward].rank) { 
					#ifdef DEBUG_MODE	
						cout << "@" << " "; 
					#endif
					continue; 
				}
				#ifdef DEBUG_MODE	
					else{
						cout << " ";
					}
				#endif
				if (dist_backward[eit->target] > dist_backward[vid_backward] + eit->weight){
					dist_backward[eit->target] = dist_backward[vid_backward] + eit->weight;
					if (dist_backward[eit->target] >= min_distance) continue;
					if (!nodelist[eit->target].flagc && !backward_insert) continue;
					backward_queue.push(make_pair(dist_backward[eit->target], eit->target));
					nodelist[eit->target].flagd = true;
					// check meet: f->m<-b
					if (nodelist[eit->target].flagc && min_distance > dist_forward[eit->target]\
						+ dist_backward[eit->target]) {
						#ifdef DEBUG_MODE	
							cout << "update distance!" << endl;
						#endif
						min_distance = dist_forward[eit->target] + dist_backward[eit->target];
						meeting_node = eit->target;
					}					
				}
			}
			
			forall_outshortcuts(graph, vid_backward, eit){
				if (eit->target == target) continue;					
				if (eit->target == source) {
					if ( dist_backward[vid_backward] + eit->weight < min_distance ) {
						min_distance = dist_backward[vid_backward] + eit->weight;
						meeting_node = source;
					}
					continue;					
				}			
				#ifdef DEBUG_MODE		
					cout << eit->target;
				#endif				
				if (vid_backward != target && nodelist[eit->target].rank < nodelist[vid_backward].rank) { 
					#ifdef DEBUG_MODE	
						cout << "@" << " "; 
					#endif
					continue; 
				}
				#ifdef DEBUG_MODE	
					else{
						cout << " ";
					}
				#endif
				if (dist_backward[eit->target] > dist_backward[vid_backward] + eit->weight){
					dist_backward[eit->target] = dist_backward[vid_backward] + eit->weight;
					if (dist_backward[eit->target] >= min_distance) continue;
					if (!nodelist[eit->target].flagc && !backward_insert) continue;
					backward_queue.push(make_pair(dist_backward[eit->target], eit->target));
					nodelist[eit->target].flagd = true;
					 // check meet: f->m<-b
					if (nodelist[eit->target].flagc && min_distance > dist_forward[eit->target]\
						+ dist_backward[eit->target]) {
						#ifdef DEBUG_MODE	
							cout << "update distance!" << endl;
						#endif
						min_distance = dist_forward[eit->target] + dist_backward[eit->target];
						meeting_node = eit->target;
					}					
				}				
			}
			#ifdef DEBUG_MODE
				cout << endl;
			#endif			
			
			break;
		}		
	}
	
	//cout << "num visied cont = " << num_visited << endl;
	
	#ifdef DEBUG_MODE
		cout << "min_distance = " << min_distance << endl;
		cout << "meeting_node = " << meeting_node << endl;
	#endif
	
	return min_distance;
}
	
double Query::run(VertexID source, VertexID target, vector<EdgeID>& inneredge){
	return BiDijkstra(source, target, inneredge);
}
