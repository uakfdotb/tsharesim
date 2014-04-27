#include "GraphUtil.h"

int GraphUtil::BiBFSDist_ptr(Graph& graph, int src, int trg, int* dist, int* que, int& ref, int radius, int& searchspace) {
	searchspace=0;
	if (src==trg) return 0;
	int u, val, tmp, index=0, endindex=0, nid, fradius, bradius, fdist, bdist;
	int findex, fendindex, bindex, bendindex, sdist=MAX_VAL, left, right;
	fradius = radius/2; bradius = radius-fradius;
	// perform forward search from src
	int gsize = graph.num_vertices();
	ref += radius+1;
	que[0]=src; que[gsize-1]=trg;
	dist[src]=ref; dist[trg]=-ref;	// forward distance positive value, backward is negative value
	findex=0; fendindex=1; bindex=gsize-1; bendindex=bindex-1;
	bool forward = true; // forward search
	while (findex<fendindex || bindex>bendindex) {
		if (forward) {
			u = que[findex];
			findex++;
			val = dist[u];
			fdist = val;
		}
		else {
			u = que[bindex];
			bindex--;
			val = dist[u];
			bdist = val;
		}
		if (forward) {
			forall_outneighbors(graph,u,eit) {
				nid=(eit)->target;
				// consider stop condition
				if (dist[nid]<=-ref) {
					if (sdist>val+1-ref-dist[nid]-ref)
						sdist = val+1-ref-dist[nid]-ref;
				}
				// update distance and queue
				if (abs(dist[nid])<ref) {
					dist[nid] = val+1;
					fdist = val+1;
					que[fendindex++] = nid;
				}
			}
		}
		else {
			forall_inneighbors(graph,u,eit) {
				nid=(eit)->target;
				// consider stop condition
				if (dist[nid]>=ref) {
					if (sdist>1-val-ref+dist[nid]-ref)
						sdist = 1-val-ref+dist[nid]-ref;
				}
				// update distance and queue
				if (abs(dist[nid])<ref) {
					dist[nid] = val-1;
					bdist = 1-val;
					que[bendindex--] = nid;
				}
			}
		}
		
		// check stop condition
		if (findex>=fendindex || bindex<=bendindex) {
			break;
		}
		
		tmp = dist[que[findex]]-dist[que[bindex]];
		if (tmp-ref-ref>=sdist) 
			break;
		else if (tmp-ref-ref>radius) {
			break;
		}		

		// change search direction by considering the search space
		/*
		if (bindex-bendindex>fendindex-findex) {
			if (fendindex>findex) forward=true;
			else forward=false;
		}
		else {
			if (bindex>bendindex) forward=false;
			else forward=true;
		}
		*/
			
		if (forward && bindex>bendindex)
			forward=false;
		else if (!forward && findex<fendindex)
			forward=true;
		
	}
	searchspace = fendindex+gsize-bendindex-1-2;
	return sdist>radius?-1:sdist;
}

int GraphUtil::BiBFSDist(Graph& graph, int src, int trg, vector<int>& dist, vector<int>& que,
		int& ref, int radius, int& searchspace) {
	searchspace=0;
	if (src==trg) return 0;
	int u, val, tmp, index=0, endindex=0, nid, fradius, bradius, fdist, bdist;
	int findex, fendindex, bindex, bendindex, sdist=MAX_VAL, left, right;
	fradius = radius/2; bradius = radius-fradius;
	// perform forward search from src
	ref += radius+1;
	int preref = ref;
	que[0]=src; que[que.size()-1]=trg;
	dist[src]=ref; dist[trg]=-ref;	// forward distance positive value, backward is negative value
	findex=0; fendindex=1; bindex=que.size()-1; bendindex=bindex-1;
	bool forward = true; // forward search
	while (findex<fendindex || bindex>bendindex) {
		if (forward) {
			u = que[findex];
			findex++;
			val = dist[u];
			fdist = val;
		}
		else {
			u = que[bindex];
			bindex--;
			val = dist[u];
			bdist = val;
		}
		if (forward) {
			forall_outneighbors(graph,u,eit) {
				nid=(eit)->target;
				// consider stop condition
				if (dist[nid]<=-ref) {
					if (sdist>val+1-ref-dist[nid]-ref) {
						sdist = val+1-ref-dist[nid]-ref;
					}
				}
				// update distance and queue
				if (abs(dist[nid])<ref) {
					dist[nid] = val+1;
					fdist = val+1;
					que[fendindex++] = nid;
				}
			}
		}
		else {
			forall_inneighbors(graph,u,eit) {
				nid=(eit)->target;
				// consider stop condition
				if (dist[nid]>=ref) {
					if (sdist>1-val-ref+dist[nid]-ref) {
						sdist = 1-val-ref+dist[nid]-ref;
					}
				}
				// update distance and queue
				if (abs(dist[nid])<ref) {
					dist[nid] = val-1;
					bdist = 1-val;
					que[bendindex--] = nid;
				}
			}
		}
		
		// check stop condition
		if (findex>=fendindex || bindex<=bendindex) {
			break;
		}
		
		tmp = dist[que[findex]]-dist[que[bindex]];
		if (tmp-ref-ref>=sdist) 
			break;
		else if (tmp-ref-ref>radius) {
			break;
		}		

		// change search direction by considering the search space
		/*
		if (bindex-bendindex>fendindex-findex) {
			if (fendindex>findex) forward=true;
			else forward=false;
		}
		else {
			if (bindex>bendindex) forward=false;
			else forward=true;
		}
		*/
			
		if (forward && bindex>bendindex)
			forward=false;
		else if (!forward && findex<fendindex)
			forward=true;
		
	}
	searchspace = fendindex+que.size()-bendindex-1-2;

	#ifdef GU_DEBUG
	for (int i = 0; i < fendindex; i++)
		cout << que[i] << " ";
	cout << "; ";
	for (int i = que.size()-1; i > bendindex; i--)
		cout << que[i] << " ";
	cout << endl;
	#endif

	return sdist>radius?-1:sdist;
}

vector<int> GraphUtil::BiBFSPath(Graph& graph, int src, int trg, vector<int>& dist, vector<int>& que, vector<int>& prev,
		int& ref, int radius, int& searchspace) {
	searchspace=0;
	vector<int> path;
	if (src==trg) {
		path.push_back(src);
		return path;
	}
	int u, val, tmp, index=0, endindex=0, nid, fradius, bradius, fdist, bdist;
	int findex, fendindex, bindex, bendindex, sdist=MAX_VAL, left, right, fpoint, bpoint;
	fradius = radius/2; bradius = radius-fradius;
	// perform forward search from src
	ref += radius+1;
	int preref = ref;
	que[0]=src; que[que.size()-1]=trg;
	prev[src]=-1; prev[trg]=-1;
	dist[src]=ref; dist[trg]=-ref;	// forward distance positive value, backward is negative value
	findex=0; fendindex=1; bindex=que.size()-1; bendindex=bindex-1;
	bool forward = true; // forward search
	while (findex<fendindex || bindex>bendindex) {
		if (forward) {
			u = que[findex];
			findex++;
			val = dist[u];
			fdist = val;
		}
		else {
			u = que[bindex];
			bindex--;
			val = dist[u];
			bdist = val;
		}
		if (forward) {
			forall_outneighbors(graph,u,eit) {
				nid=(eit)->target;
				// consider stop condition
				if (dist[nid]<=-ref) {
					if (sdist>val+1-ref-dist[nid]-ref) {
						sdist=val+1-ref-dist[nid]-ref;
						fpoint = u;
						bpoint = nid;
					}
				}
				// update distance and queue
				if (abs(dist[nid])<ref) {
					dist[nid] = val+1;
					fdist = val+1;
					que[fendindex++] = nid;
					prev[nid] = u;
				}
			}
		}
		else {
			forall_inneighbors(graph,u,eit) {
				nid=(eit)->target;
				// consider stop condition
				if (dist[nid]>=ref) {
					if (sdist>1-val-ref+dist[nid]-ref) {
						sdist=1-val-ref+dist[nid]-ref;
						bpoint = u;
						fpoint = nid;
					}
				}
				// update distance and queue
				if (abs(dist[nid])<ref) {
					dist[nid] = val-1;
					bdist = 1-val;
					que[bendindex--] = nid;
					prev[nid] = u;
				}
			}
		}
		
		// check stop condition
		if (findex>=fendindex || bindex<=bendindex) {
			break;
		}
		
		tmp = dist[que[findex]]-dist[que[bindex]];
		if (tmp-ref-ref>=sdist) 
			break;
		else if (tmp-ref-ref>radius) {
			break;
		}		

		// change search direction by considering the search space
		/*
		if (bindex-bendindex>fendindex-findex) {
			if (fendindex>findex) forward=true;
			else forward=false;
		}
		else {
			if (bindex>bendindex) forward=false;
			else forward=true;
		}
		*/
			
		if (forward && bindex>bendindex)
			forward=false;
		else if (!forward && findex<fendindex)
			forward=true;
		
	}
	searchspace = fendindex+que.size()-bendindex-1-2;
	if (sdist>radius) return path;

	// recover shortest path
	while (prev[fpoint]!=-1) {
		path.push_back(fpoint);
		fpoint = prev[fpoint];
	}
	path.push_back(src);
	reverse(path.begin(),path.end());
	while (prev[bpoint]!=-1) {
		path.push_back(bpoint);
		bpoint = prev[bpoint];
	}
	path.push_back(trg);
	return path;
}

// bidirectional BFS serves as benchmark
int GraphUtil::BiBFSDist_simple(Graph& graph, int src, int trg, vector<int>& dist, vector<int>& que,
		int& ref, int radius, int& searchspace) {
	int u, val, index=0, endindex=0, nid, fradius, bradius;
	searchspace=0;
	fradius = radius/2; bradius = radius-fradius;
	// perform forward search from src
	ref += radius+1;
	int preref = ref;
	que[0]=src;
	dist[src]=ref;
	endindex=1;
	index=0;
	while (index<endindex) {
		u = que[index];
		index++;
		val = dist[u];
		forall_outneighbors(graph,u,eit) {
			nid=(eit)->target;
			if (dist[nid]<ref) {
				if (nid==trg) {
					searchspace += endindex;
					return val+1-ref;
				}
				dist[nid]=val+1;
				if (val+1-ref<fradius)
					que[endindex++]=nid;
			}
		}
	}	
	searchspace = endindex-1;
	// perform backward search
	ref += radius+1;
	que[0]=trg;
	dist[trg]=ref;
	endindex=1;
	index=0;
	while (index<endindex) {
		u = que[index];
		index++;
		val = dist[u];
		searchspace++;
		forall_inneighbors(graph,u,eit) {
			nid=(eit)->target;
			if (dist[nid]<ref) {
				if (dist[nid]>=preref) {
					searchspace += endindex;
					return val+1-ref+dist[nid]-preref;
				}
				dist[nid]=val+1;
				if (val+1-ref<bradius)
					que[endindex++]=nid;
			} 
		}
	}
	searchspace += endindex-1;
	return -1;
}

// undirected graph
int GraphUtil::BFSDist(Graph& graph, int src, int trg, vector<int>& dist, vector<int>& que, int& ref, 
		int radius, int& searchspace) {
	int u, val, index=0, endindex=0, nid; 
	searchspace=0;
	ref += radius+1;
	que[0]=src;
	dist[src]=ref;
	endindex=1;
	index=0;
	while (index<endindex) {
		u = que[index];
		index++;
		val = dist[u];
		forall_outneighbors(graph,u,eit) {
			nid=(eit)->target;
			if (dist[nid]<ref) {
				dist[nid]=val+1;
				if (trg==nid) {
					searchspace += endindex;
					return val+1-ref;
				}
				if (val+1-ref<radius) {
					que[endindex++]=nid;
				}
			} 
		}
	}
	searchspace = endindex-1;
	return -1;
}

vector<int> GraphUtil::BFSPath(Graph& graph, int src, int trg, vector<int>& dist, vector<int>& que, vector<int>& prev, int& ref, 
		int radius, int& searchspace) {
	vector<int> path;
	if (src==trg) {
		path.push_back(src);
		return path;
	}
	int u, val, index=0, endindex=0, nid; 
	searchspace=0;
	ref += radius+1;
	que[0]=src;
	dist[src]=ref;
	prev[src]=-1;
	endindex=1;
	index=0;
	while (index<endindex) {
		u = que[index];
		index++;
		val = dist[u];
		forall_outneighbors(graph,u,eit) {
			nid=(eit)->target;
			if (dist[nid]<ref) {
				dist[nid]=val+1;
				if (trg==nid) {
					searchspace += endindex;
					// recover shortest path
					prev[nid]=u;
					while (prev[nid]!=-1) {
						path.push_back(nid);
						nid = prev[nid];
					}
					path.push_back(src);
					reverse(path.begin(),path.end());
					return path;
				}
				if (val+1-ref<radius) {
					que[endindex++]=nid;
					prev[nid]=u;
				}
			} 
		}
	}
	searchspace = endindex-1;
	return path;
}

bool GraphUtil::checkUndirectedGraph(Graph& graph) {
	srand48(time(NULL));
	int gsize = graph.num_vertices();
	int s, t;
	for (int i = 0; i < 50; i++) {
		s = lrand48() % gsize;
		forall_outneighbors(graph,s,eit) {
			if (!graph.hasEdge(eit->target,s))
				return false;
		}
	}
	return true;
}

// return maximum degree of non-core nodes
vector<int> GraphUtil::coreNodeSelection_degree(Graph& graph, int corenum, bit_vector* corenodes) {
	int gsize = graph.num_vertices();
	assert(corenum<=gsize);
	int degree = 0;
	multimap<int,int> degmap;
	multimap<int,int>::reverse_iterator mit;
	for (int i = 0; i < gsize; i++) {
		degree = graph.out_degree(i)+graph.in_degree(i);
		degmap.insert(make_pair(degree,i));
	}
	vector<int> results;
	int index = 0;
	mit = degmap.rbegin();
	results.push_back(mit->first);
	while (index<corenum && mit!=degmap.rend()) {
		corenodes->set_one(mit->second);
		index++;
		mit++;
	}
	results.push_back(mit->first);
	return results;
}

vector<int> GraphUtil::coreNodeSelectionWithNodeMap(Graph& graph, int corenum, int type, bit_vector* corenodes, vector<int>& nodemap) {
	int gsize = graph.num_vertices();
	assert(corenum<=gsize);
	int degree=0, outd, ind;
	multimap<int,int> degmap;
	multimap<int,int>::reverse_iterator mit;
	for (int i = 0; i < gsize; i++) {
		if (type==0)
			degree = graph.out_degree(i)+graph.in_degree(i);
		else {
			outd = graph.out_degree(i);
			ind = graph.in_degree(i);
			if (outd==0) outd=1;
			if (ind==0) ind=1;
			degree = outd*ind;
		}
		degmap.insert(make_pair(degree,i));
	}
	nodemap.clear();
	nodemap = vector<int>(gsize,-1);
	vector<int> results;
	int index = 0;
	mit = degmap.rbegin();
	results.push_back(mit->first);
	while (index<corenum && mit!=degmap.rend()) {
		corenodes->set_one(mit->second);
		nodemap[mit->second] = index;
		index++;
		mit++;
	}
	results.push_back(mit->first);
	while (mit!=degmap.rend()) {
		nodemap[mit->second] = index;
		mit++;
		index++;
	}
//	srand (unsigned(time (NULL)));
//	random_shuffle(nodemap.begin()+corenum+1,nodemap.end());
	assert(nodemap.size()==gsize);
	return results;
}

void GraphUtil::computePairwiseDist(Graph& graph, bit_vector* corenodes, int radius, map<int,int>& coreindex, vector<vector<int> >& coredist) {
	int gsize = graph.num_vertices(), ref = 0;
	int src, u, val, index=0, endindex=0, nid; 
	vector<int> que = vector<int>(gsize,0);
	vector<int> dist = vector<int>(gsize,0);
	int counter = 0;
	for (int i = 0; i < gsize; i++) {
		if (corenodes->get(i)) {
			coreindex.insert(make_pair(i,counter));
			counter++;
		}
	}
	coredist = vector<vector<int> >(coreindex.size(),vector<int>());
	map<int,int>::iterator mit;
	for (mit = coreindex.begin(); mit != coreindex.end(); mit++) {
		src = mit->first;
		coredist[mit->second] = vector<int>(mit->second+1,-1);
		ref += radius+1;
		que[0]=src;
		dist[src]=ref;
		endindex=1;
		index=0;
		while (index<endindex) {
			u = que[index];
			index++;
			val = dist[u];
			forall_outneighbors(graph,u,eit) {
				nid=(eit)->target;
				if (dist[nid]<ref) {
					dist[nid]=val+1;
					if (corenodes->get(nid)) {
						if (src>nid)
							coredist[mit->second][coreindex[nid]] = val+1;
					}					
					if (val+1-ref<radius) {
						que[endindex++]=nid;
					}
				} 
			}
		}	
	}
}

int GraphUtil::generateQueriesByNode(Graph& graph, int src, int radius, int num_nodes, vector<int>& dist, int& ref,
		vector<int>& srcvec, vector<int>& trgvec) {
	int u, val, index=0, nid, size=0; 
	ref += radius+1;
	vector<int> que;
	que.push_back(src);
	dist[src]=ref;
	index=0;
	while (index<que.size()) {
		u = que[index];
		index++;
		val = dist[u];
		forall_outneighbors(graph,u,eit) {
			nid=(eit)->target;
			if (dist[nid]<ref) {
				dist[nid]=val+1;
				if (val+1-ref<radius) {
					que.push_back(nid);
				}
				else if (val+1-ref==radius) {
					srcvec.push_back(src);
					trgvec.push_back(nid);
					size++;
					if (size==num_nodes)
						return size;
				}
			} 
		}
	}
	return size;
}

void GraphUtil::generateQueriesByDistanceAndOutput(Graph& graph, int radius, int querynum, ostream& out) {
	int index=0, ref=0, gsize=graph.num_vertices(), size, src, num_nodes;
	vector<int> dist = vector<int>(gsize,0);
	vector<int> srcvec, trgvec;
//	set<int> srcset;
	int timer=0;
	bool stop=false;
	srand (time(NULL));
	while (srcvec.size()<querynum) {
		src = rand()%gsize;
//		if (srcset.find(src)!=srcset.end()) continue;
//		srcset.insert(src);
		num_nodes = 50+rand()%50;
		if (querynum-srcvec.size()<num_nodes)
			num_nodes = querynum-srcvec.size();
		size = generateQueriesByNode(graph, src, radius, num_nodes, dist, ref, srcvec, trgvec);
		if (size==0) {
			if (!stop) {
				timer = 1;
				stop = true;
			}
			else
				timer++;
		}
		else
			stop = false;
		if (timer==20000) {
			cout << "actual query size=" << srcvec.size() << endl;
			break;
		}
	}
	for (int i = 0; i < srcvec.size(); i++) {
		out << srcvec[i] << " " << radius << " " << trgvec[i] << endl;
	}
	out.flush();
}

void GraphUtil::generateQueriesByDistance(Graph& graph, int radius, int querynum, vector<int>& srcvec, vector<int>& trgvec) {
	int index=0, ref=0, gsize=graph.num_vertices(), size, src, num_nodes;
	vector<int> dist = vector<int>(gsize,0);
//	set<int> srcset;
	srand (time(NULL));
	while (srcvec.size()<querynum) {
		src = rand()%gsize;
		num_nodes = 50+rand()%50;
		if (querynum-srcvec.size()<num_nodes)
			num_nodes = querynum-srcvec.size();
		size = generateQueriesByNode(graph, src, radius, num_nodes, dist, ref, srcvec, trgvec);
	}
	assert(srcvec.size()==querynum);
}

int GraphUtil::generateRandQueriesByNode(Graph& graph, int src, int num_nodes, int radius, vector<int>& dist, int& ref,
		vector<int>& srcvec, vector<int>& trgvec, vector<int>& distvec) {
	int u, val, index=0, nid, size=0; 
	ref += radius+1;
	vector<int> que;
	que.push_back(src);
	dist[src]=ref;
	index=0;
	while (index<que.size()) {
		u = que[index];
		index++;
		val = dist[u];
		forall_outneighbors(graph,u,eit) {
			nid=(eit)->target;
			if (dist[nid]<ref) {
				dist[nid]=val+1;
				if (val+1-ref<radius+2) {
					que.push_back(nid);
				}
				const int coin = rand()%10000;
				if (coin<600) {
					srcvec.push_back(src);
					trgvec.push_back(nid);
					distvec.push_back(val+1-ref);
					size++;
					if (size==num_nodes)
						return size;
				}
			} 
		}
	}
	return size;
}

void GraphUtil::generateRandomQueries(Graph& graph, int radius, int querynum, ostream& out) {
	int index=0, ref=0, gsize=graph.num_vertices(), size, src, num_nodes;
	vector<int> dist = vector<int>(gsize,0);
	vector<int> srcvec, trgvec, distvec;
	set<int> srcset;
	srand (time(NULL));
	while (srcvec.size()<querynum) {
		src = rand()%gsize;
		num_nodes = 150+rand()%50;
		if (querynum<num_nodes)
			num_nodes = querynum;
		size = generateRandQueriesByNode(graph, src, num_nodes, radius, dist, ref, srcvec, trgvec, distvec);
	}
	for (int i = 0; i < srcvec.size(); i++) {
	//	out << srcvec[i] << "\t" << trgvec[i] << "\t" << distvec[i] << endl;
		out << srcvec[i] << " " << distvec[i] << " " << trgvec[i] << endl; // to fit sketch query format
	}
	out.flush();
}

int GraphUtil::generateRandQueriesByNode1(Graph& graph, int src, int num_nodes, int radius, vector<int>& dist, int& ref,
		vector<int>& srcvec, vector<int>& trgvec, vector<int>& distvec) {
	int u, val, index=0, nid, size=0; 
	ref += radius+1;
	vector<int> que;
	que.push_back(src);
	dist[src]=ref;
	index=0;
	while (index<que.size()) {
		u = que[index];
		index++;
		val = dist[u];
		forall_outneighbors(graph,u,eit) {
			nid=(eit)->target;
			if (dist[nid]<ref) {
				dist[nid]=val+1;
				if (val+1-ref<radius) {
					que.push_back(nid);
				}
			} 
		}
	}
	if (que.size()<500) return 0;
	
	int left = 0, trg;
	int seg = (int)(que.size()*100/100), start=que.size()-seg;  
	while (left < num_nodes) {
		index = rand()%seg+start;
		if (index<0 || index>=que.size()) continue;
		trg = que[index];
		srcvec.push_back(src);
		trgvec.push_back(trg);
		distvec.push_back(dist[trg]-ref);
		++left;
	}
	
	return size;
}

void GraphUtil::generateRandomQueries1(Graph& graph, int radius, int querynum, ostream& out) {
	int index=0, ref=0, gsize=graph.num_vertices(), size, src, num_nodes;
	vector<int> dist = vector<int>(gsize,0);
	vector<int> srcvec, trgvec, distvec;
	set<int> srcset;
	srand (time(NULL));
	while (srcvec.size()<querynum) {
		src = rand()%gsize;
		num_nodes = 100+rand()%50;
		if (srcvec.size()+num_nodes>querynum)
			num_nodes = querynum-srcvec.size();
		generateRandQueriesByNode1(graph, src, num_nodes, radius+4, dist, ref, srcvec, trgvec, distvec);
	}
	for (int i = 0; i < srcvec.size(); i++) {
	//	out << srcvec[i] << "\t" << trgvec[i] << "\t" << distvec[i] << endl;
		out << srcvec[i] << " " << distvec[i] << " " << trgvec[i] << endl; // to fit sketch query format
	}
	out.flush();
}

void GraphUtil::generateRandomBoundedQueries(Graph& graph, int radius, int querynum, ostream& out) {
	int index=0, ref=0, gsize=graph.num_vertices(), size, src, num_nodes, prev;
	vector<int> dist = vector<int>(gsize,0);
	vector<int> srcvec, trgvec, distvec;
	bool flag=false;
	int unchangecount=0;
	srand (time(NULL));
	while (srcvec.size()<querynum) {
		if (unchangecount>=10000)
			break;
		prev = srcvec.size();
		src = rand()%gsize;
		num_nodes = 100+rand()%50;
		if (srcvec.size()+num_nodes>querynum)
			num_nodes = querynum-srcvec.size();
		generateRandQueriesByNode1(graph, src, num_nodes, radius+1, dist, ref, srcvec, trgvec, distvec);
		if (prev==srcvec.size()) {
			flag = true;
			unchangecount++;
		}
		else {
			flag = false;
			unchangecount = 0;
		}
	}
	for (int i = 0; i < srcvec.size(); i++) {
	//	out << srcvec[i] << "\t" << trgvec[i] << "\t" << distvec[i] << endl;
		out << srcvec[i] << " " << distvec[i] << " " << trgvec[i] << endl; // to fit sketch query format
	}
	out.flush();
}


void GraphUtil::generateCombination(int radius, vector<vector<unsigned int> >& combination) {
	for (int sum = 2; sum <= radius; sum++) {
		for (int i = 1; i < sum; i++) {
			combination.push_back(vector<unsigned int>(3,0));
			combination.back()[0] = i-1;
			combination.back()[1] = sum-i-1;
			combination.back()[2] = sum;
		}
	}
	#ifdef GU_DEBUG
	for (int i = 0; i < combination.size(); i++) {
		cout << "[" << combination[i][0] << " " << combination[i][1] << " " << combination[i][2] << "] ";
	}
	cout << endl;
	#endif
}

void GraphUtil::DFSExpand(Graph& graph, int vid, bit_vector* visited, bit_vector* corenodes, int& size) {
	size++;
	visited->set_one(vid);
	forall_outneighbors(graph,vid,eit) {
		if (!visited->get(eit->target) && !corenodes->get(eit->target)) {
			DFSExpand(graph,eit->target,visited,corenodes,size);
		}
	}
}

// return the size of each connected component after removing high degree nodes
vector<int> GraphUtil::degreeDecomposition(Graph& graph, int coresize) {
	int gsize = graph.num_vertices(), size;
	bit_vector* corenodes = new bit_vector(gsize);
	coreNodeSelection_degree(graph, coresize, corenodes);
	bit_vector* visited = new bit_vector(gsize);
	vector<int> results;
	for (int i = 0; i < gsize; i++) {
		if (corenodes->get(i)||visited->get(i)) continue;
		size = 0;
		DFSExpand(graph,i,visited,corenodes,size);
		results.push_back(size);
	}
	delete visited;
	delete corenodes;
	return results;
}

void GraphUtil::readCpuMemInfo(vector<string>& cpumem) {
	string filepath = string("/proc/meminfo");
	ifstream infile(filepath.c_str(),ios_base::in);
	if (!infile) {
		cerr << "Cannot open " << filepath << endl;
		return;
	}
	string buf;
	while (getline(infile,buf)) {
		if (buf.find("MemTotal:")!=string::npos) {
			cpumem.push_back(buf);
			break;
		}
	}
	infile.close();
	filepath = string("/proc/cpuinfo");
	ifstream cpufile(filepath.c_str(),ios_base::in);
	if (!cpufile) {
		cerr << "Cannot open " << filepath << endl;
		return;
	}
	buf = string("");
	while (getline(cpufile,buf)) {
		if (buf.find("model name")!=string::npos
			|| buf.find("cpu MHz")!=string::npos) {
			cpumem.push_back(buf);
		}
	}
	cpufile.close();
}

void GraphUtil::BFSDistribute(Graph& graph, int src, vector<int>& dist, int& ref, int radius, vector<int>& distribute) {
	int u, val, index=0, endindex=0, nid; 
	ref += radius+1;
	vector<int> que;
	que.push_back(src);
	dist[src]=ref;
	endindex=1;
	index=0;
	while (index<que.size()) {
		u = que[index];
		index++;
		val = dist[u];
		forall_outneighbors(graph,u,eit) {
			nid=(eit)->target;
			if (dist[nid]<ref) {
				dist[nid]=val+1;
				distribute[val+1-ref]++;
				if (val+1-ref<radius) {
					que.push_back(nid);
				}
			} 
		}
	}
}

double GraphUtil::effective90diameter(Graph& graph, int seednum) {
	srand (time(NULL));
	int gsize = graph.num_vertices();
	set<int> seeds;
	int left=0, node;
	while (left<seednum) {
		node = rand()%gsize;
		if (seeds.find(node)==seeds.end()) {
			seeds.insert(node);
			left++;
		}
	}
	
	set<int>::iterator sit;
	vector<int> dist = vector<int>(gsize,0);
	int radius=40, ref = 0;
	vector<int> distribute = vector<int>(radius,0);
	for (sit = seeds.begin(); sit != seeds.end(); sit++) {
		BFSDistribute(graph,*sit,dist,ref,radius+1,distribute);
	}
	
	long sum=0, psum=0, prev=0;
	for (int i = 0; i < distribute.size(); i++) {
//		cout << i << ": " << distribute[i] << endl;
		sum += distribute[i];
	}
	for (int i = 0; i < distribute.size(); i++) {
		prev = psum;
		psum += distribute[i];
		if (psum>=(sum*0.9))
			return ((((sum*0.9)-prev)*1.0)/((psum-prev)*1.0))+i-1;
	}
	return 0;
}

void GraphUtil::generateRandomQuery(int num_vertex, int num_query, vector<VertexID>& query){
	for ( int i = 0; i < num_query; i++ ) {
		while(true){
			int source = rand() % num_vertex;
			int target = rand() % num_vertex;
			if (source != target) {
				query.push_back(source);
				query.push_back(target);
				break;
			}
		}
	}
}

void GraphUtil::generateRandomBoundedQuery(Graph& graph, int num_vertex, int num_query, vector<VertexID>& query, vector<double>& can_dist){
	while (query.size() < 2*num_query) {
		int source = rand() % num_vertex;
		
		// randomly research a node within the distance contraint;
		queue<pair<double, VertexID> > Q;
		Q.push(make_pair(0, source));
		vector<VertexID> candidate;
		vector<double> tmp_dist;
		
		double dist = 0;
		while (Q.size() > 0) {
			VertexID can = Q.front().second;
			double wei = Q.front().first;
			Q.pop();
			
			forall_outneighbors(graph, can, eit){
				if (eit->target == source) continue;
				Q.push(make_pair(eit->weight+wei, eit->target));
				candidate.push_back(eit->target);
				tmp_dist.push_back(eit->weight+wei);
			}
			
			if (candidate.size() > 100 ) break;
		}
		
		// randomly pick one from the candidate;
		if (candidate.size() == 0) continue;
		int rindex = rand() % candidate.size();
		query.push_back(source);
		query.push_back(candidate[rindex]);
		can_dist.push_back(tmp_dist[rindex]);
	}
}

void GraphUtil::generateQuery(Graph& graph, int num_vertex, int num_query, ostream& out, vector<VertexID>& query){
	vector<VertexID>(0).swap(query);
	vector<double> dist;

	// generate local queries;
	// generateRandomBoundedQuery(graph, num_vertex, num_query, query, dist);
	// for ( int i = 0; i < query.size(); i+=2 ){
		// out << query[i] << "\t" << query[i+1] << "\t" << dist[i/2] << endl;
	// }
	
	// generate random queries;
	generateRandomQuery(num_vertex, num_query, query);
	for ( int i = 0; i < query.size(); i+=2 ){
		out << query[i] << "\t" << query[i+1] << "\t" << 1 << endl;
	}
	
	out.flush();
}









