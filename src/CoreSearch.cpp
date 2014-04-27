#include "CoreSearch.h"


CoreSearch::CoreSearch(Graph& graph, int r, int corenum, int ordertype, bool _estflag) 
		: g(graph), radius(r), coresize(corenum), nodeordertype(ordertype), estimate(_estflag) {
	gsize = g.num_vertices();
	corenodes = new bit_vector(gsize);
	unsigned int size = corenum*corenum;
	corepaths = new coresizetype*[size];
	for (unsigned int i = 0; i < size; i++) {
		corepaths[i] = NULL;
	}
	GraphUtil::generateCombination(radius,combination);
	vm_usage = 0; resident_set = 0;
	ref = 0; num_entries = 0;
	prev = NULL;
	dist = NULL;
	mark = NULL;
	pque = NULL;
	#ifdef CS_STAT
	noncorequeries = 0;
	corehitqueries = 0;
	#endif
}

CoreSearch::CoreSearch(Graph& graph, int r, int corenum) : g(graph), radius(r), coresize(corenum) {
	gsize = g.num_vertices();
	corenodes = new bit_vector(gsize);
	unsigned int size = corenum*corenum;
	corepaths = new coresizetype*[size];
	for (unsigned int i = 0; i < size; i++) {
		corepaths[i] = NULL;
	}
	GraphUtil::generateCombination(radius,combination);
	vm_usage = 0; resident_set = 0;
	ref = 0; num_entries = 0;
	nodeordertype = 0;
	estimate = false;
	prev = NULL;
	dist = NULL;
	mark = NULL;
	pque = NULL;
	#ifdef CS_STAT
	noncorequeries = 0;
	corehitqueries = 0;
	#endif
}

CoreSearch::~CoreSearch() {
	if (prev!=NULL) delete[] prev;
	if (dist!=NULL) delete[] dist;
	if (pque!=NULL) delete[] pque;
	delete[] routetables;
	unsigned int size = coresize*coresize;
	for (unsigned int i = 0; i < size; i++) {
		if (corepaths[i]!=NULL)
			delete []corepaths[i];
	}
	delete[] corepaths;
	#ifdef CS_DEBUG
	cout << "runtime1=" << runtime1 << " runtime2=" << runtime2 << " runtime3=" << runtime << endl;
	#endif
	#ifdef CS_STAT
	cout << "num_noncorequery=" << noncorequeries << " num_corehitquery=" << corehitqueries << endl;
	#endif
}

void CoreSearch::setEstimate(bool flag) {
	estimate = flag;
}

void CoreSearch::initHelper() {
	prev = new int[gsize];
	memset(prev,0,gsize*sizeof(int));
	preport = new ushort[gsize];
	memset(preport,0,gsize*sizeof(ushort));
	dist = new int[gsize];
	memset(dist,0,gsize*sizeof(int));
	mark = new int[gsize];
	memset(mark,0,gsize*sizeof(int));
}

void CoreSearch::destroyHelper() {
	delete[] prev;
	delete[] preport;
	delete[] dist;
	delete[] mark;
	delete[] pque;
}

void CoreSearch::corelabeling() {
	double vm, res;
	vector<int> que = vector<int>(gsize,0);
	vector<vector<vector<routeEntry> > > tmp_outtable = vector<vector<vector<routeEntry> > >(gsize);
	vector<vector<vector<routeEntry> > > tmp_intable = vector<vector<vector<routeEntry> > >(gsize);
	for (int i = coresize; i < gsize; i++) {
		tmp_outtable[i] = vector<vector<routeEntry> >(radius,vector<routeEntry>());
		tmp_intable[i] = vector<vector<routeEntry> >(radius,vector<routeEntry>());
	}
	for (int i = 0; i < coresize; i++) {
		coreLabelingNode(i,que,tmp_intable,true);
		coreLabelingNode(i,que,tmp_outtable,false);
		if (i%500==0) {
			cout << i+1 << " corenodes have been processed" << endl;
			Util::process_mem_usage(vm, res);
			if (vm>vm_usage) vm_usage=vm;
			if (res>resident_set) resident_set=res;
		}
	}
	que.clear();
	cout << "Converting tmp_table to concise representation..." << endl;
	convertTables(tmp_outtable,tmp_intable);
	tmp_outtable.clear();
	tmp_intable.clear();
	Util::process_mem_usage(vm, res);
	if (vm>vm_usage) vm_usage=vm;
	if (res>resident_set) resident_set=res;
}

void CoreSearch::convertTables(vector<vector<vector<routeEntry> > >& tmp_outtable, vector<vector<vector<routeEntry> > >& tmp_intable) {
	double vm, res;
	routetables = new nodeRouteTable[gsize];
	int size=0, index=0;
	vector<ushort> tmp_index = vector<ushort>(radius+1,0);
	for (int i = 0; i < coresize; i++) {
		routetables[i].outitemindex = NULL;
		routetables[i].outrouteitems = NULL;
		routetables[i].initemindex = NULL;
		routetables[i].inrouteitems = NULL;
	}
	Util::process_mem_usage(vm, res);
	if (vm>vm_usage) vm_usage=vm;
	if (res>resident_set) resident_set=res;
	for (int i = coresize; i < gsize; i++) {
		size = 0;
		for (int j = 0; j < radius; j++) {
			tmp_index[j] = (ushort)size;
			size += tmp_outtable[i][j].size();
		}
		num_entries += size;
		if (size==0) {
			routetables[i].outrouteitems = NULL;
			routetables[i].outitemindex = NULL;
		}
		else {
			// update index
			routetables[i].outitemindex = new ushort[radius+1];
			for (int j = 0; j < radius; j++) {
				routetables[i].outitemindex[j] = tmp_index[j];
			}
			routetables[i].outitemindex[radius] = (ushort)size;
			// update routeitems
			routetables[i].outrouteitems = new routeEntry[size];
			index = 0;
			for (int j = 0; j < radius; j++) {
				for (int k = 0; k < tmp_outtable[i][j].size(); k++) {
					routetables[i].outrouteitems[index].coreid = tmp_outtable[i][j][k].coreid;
					routetables[i].outrouteitems[index].portid = tmp_outtable[i][j][k].portid;
					index++;
				}
				tmp_outtable[i][j].clear();
			}
		}
		tmp_outtable[i].clear();
		// convert intables
		size = 0;
		for (int j = 0; j < radius; j++) {
			tmp_index[j] = (ushort)size;
			size += tmp_intable[i][j].size();
		}
		num_entries += size;
		if (size==0) {
			routetables[i].inrouteitems = NULL;
			routetables[i].initemindex = NULL;
		}
		else {
			// update index
			routetables[i].initemindex = new ushort[radius+1];
			for (int j = 0; j < radius; j++) {
				routetables[i].initemindex[j] = tmp_index[j];
			}
			routetables[i].initemindex[radius] = (ushort)size;
			// update routeitems
			routetables[i].inrouteitems = new routeEntry[size];
			index = 0;
			for (int j = 0; j < radius; j++) {
				for (int k = 0; k < tmp_intable[i][j].size(); k++) {
					routetables[i].inrouteitems[index].coreid = tmp_intable[i][j][k].coreid;
					routetables[i].inrouteitems[index].portid = tmp_intable[i][j][k].portid;
					index++;
				}
				tmp_intable[i][j].clear();
			}
		}
		tmp_intable[i].clear();
	}
}

void CoreSearch::coreLabelingNode(int vid, vector<int>& que, vector<vector<vector<routeEntry> > >& tmp_table, bool isout) {
	int u, val, index, endindex, nid, mval, pindex, qindex; 
	unsigned int coreindex;
	vector<int> tmp_path;
	vector<int>::reverse_iterator rit;
	ref += radius+1;
	que[0] = vid;
	dist[vid] = ref;
	prev[vid] = -1; // for tracing shortest path
	index = 0;
	endindex = 1;
	#ifdef CS_DEBUG
	cout << "start from " << vid << " isout=" << isout << endl;
	#endif
	while (index<endindex) {
		u=que[index++];
		val=dist[u];
		#ifdef CS_DEBUG
		cout << "touch node " << u << " dist=" << val << " mark=" << mark[u] << " ref=" << ref << endl;
		#endif
		if (u!=vid) {
			// if u is corenode and traversal follows outgoing direction, then we store the shortest path from vid to u (note: vid<u)
			if (u<coresize) {
				mark[u] = ref;
				if (isout) {
					// trace shortest path from current corenode to root vid
					pindex = prev[u];
					while (prev[pindex]!=-1) {
						tmp_path.push_back(pindex);
						pindex = prev[pindex];
					}
					coreindex = computeMatrixIndex(vid,u);
					corepaths[coreindex] = new coresizetype[tmp_path.size()+1];
					// first element is the number of elements in this array
					corepaths[coreindex][0] = (coresizetype)(tmp_path.size());
					// the path should be [vid, x1, x2, x3,..., u]
					qindex = 1;
					for (rit=tmp_path.rbegin(); rit!=tmp_path.rend(); rit++,qindex++)
						corepaths[coreindex][qindex] = (coresizetype)(*rit);
					tmp_path.clear();
					#ifdef CS_DEBUG
					cout << "\tInsert path into index=" << pindex << " based on (" << vid << "," << u << ")" << endl;
					int length = corepaths[coreindex][0];
					cout << "Path(length=" << length+1 << "): " << vid << " ";
					for (int i = 1; i <= length; i++)
						cout << corepaths[coreindex][i] << " ";
					cout << u << endl;
					#endif
				}
			}
		}
		mval=mark[u];
		if (u!=vid && mval!=ref) {
			// if u is non-core node and is NEVER reached by some corenodes, then we store the port and distance to routeTable
			tmp_table[u][val-ref-1].push_back(routeEntry((unsigned int)vid,preport[u]));
			#ifdef CS_DEBUG
			cout << "\tUpdate route table of node[" << u << "]" << endl;
			cout << "Level=" << val-ref << " coreid=" << tmp_table[u][val-ref-1].back().coreid << " portid=" << tmp_table[u][val-ref-1].back().portid << endl;
			#endif
		}
		if (val-ref==radius) continue;
		// traverse all its neighbors with portid
		if (isout) {
			forall_outneighborport(g,u,eit,port_iter) {
				nid=(*eit);
				if (dist[nid]<ref) {
					dist[nid]=val+1;
					prev[nid]=u;
					// revised 11/30/2011 to gaurantee that previous node in the shortest path also record the corenode vid
					if (u==vid || mval!=ref)
						preport[nid]=(*port_iter);
					if (val+1-ref<=radius) {
						que[endindex++]=(nid);
					}
					if (mval==ref && mark[nid]!=ref)
						mark[nid] = ref;
				}
				else if (dist[nid]==val+1 && mval==ref && mark[nid]!=ref)
					mark[nid] = ref;
			}
		}
		else {
			forall_inneighborport(g,u,eit,port_iter) {
				nid=(*eit);
				if (dist[nid]<ref) {
					dist[nid]=val+1;
					prev[nid]=u;
					// revised 11/30/2011 to gaurantee that previous node in the shortest path also record the corenode vid
					if (u==vid || mval!=ref)
						preport[nid]=(*port_iter);
					if (val+1-ref<=radius) {
						que[endindex++]=(nid);
					}
					if (mval==ref && mark[nid]!=ref)
						mark[nid] = ref;
				}
				else if (dist[nid]==val+1 && mval==ref && mark[nid]!=ref)
					mark[nid] = ref;
			}
		}
	}
	
	#ifdef CS_DEBUG
	cout << "visited nodes from " << vid << " : (out=" << isout << ")" << endl;
	for (int i = 1; i < que.size(); i++) {
		if (mark[que[i]]!=ref) {
			cout << i << "\tuncovered node=" << que[i] << "\tdist=" << dist[que[i]]-ref << endl;
		}
		else 
			cout << i << "\tnode=" << que[i] << "\tdist=" << dist[que[i]]-ref << endl;
	}
	#endif
}

void CoreSearch::createLabels() {
	vector<int> fake_srcs, fake_trgs;
	createLabels(fake_srcs,fake_trgs,"");
}

void CoreSearch::createLabels(vector<int>& srcs, vector<int>& trgs, string interResultFilename) {
	double runtime;
	PerformanceTimer timer = PerformanceTimer::start();

	if (interResultFilename!="") {
		cout << "Reading existing label file..." << endl;
		timer.reset();
		ifstream infile(interResultFilename.c_str());
		if (!infile) {
			cerr << "Cannot open label file " << interResultFilename << endl;
			exit(0);
		}
		loadIntermediateResults(infile);
		infile.close();
		runtime = timer.reset();
		cout << "Done. It took " << runtime << " ms" << endl;
		#ifdef CS_DEBUG
		cout << "-------------------------------------------" << endl;
		displayCorepaths(cout);
		cout << "-------------------------------------------" << endl;
		displayRouteTable(cout);
		cout << "-------------------------------------------" << endl;
		#endif
	}
	
	if (nodemap.size()<gsize) {
		cout << "Selecting " << coresize << " nodes with largest degree as core nodes..." << endl;
		timer.reset();
		vector<int> degs = GraphUtil::coreNodeSelectionWithNodeMap(g, coresize, nodeordertype, corenodes, nodemap);
		runtime = timer.reset();
		cout << "Done. It took " << runtime << " ms" << endl;
	}
	
	cout << "Converting queries to new graph queries..." << endl;
	assert(srcs.size()==trgs.size());
	timer.reset();
	for (int i = 0; i < srcs.size(); i++) {
		srcs[i] = nodemap[srcs[i]];
		trgs[i] = nodemap[trgs[i]];
	}
	runtime = timer.reset();
	cout << "Done. It took " << runtime << " ms" << endl;
	
	#ifdef CS_DEBUG
	cout << "------------NodeMap---------------" << endl;
	for (int i = 0; i < nodemap.size(); i++) {
		cout << i << " -> " << nodemap[i] << endl;
	}
	#endif

	cout << "Reorder node ID..." << endl;
	timer.reset();
	g.reorderid(nodemap);
	runtime = timer.reset();
	cout << "Done. It took " << runtime << " ms" << endl;
	
	#ifdef CS_DEBUG
	cout << "New graph based on node degree ordering" << endl;
	g.printGraph();
	#endif
	
	cout << "Building portlist..." << endl;
	timer.reset();
	g.constructPortlist(coresize);
	runtime = timer.reset();
	cout << "Done. It took " << runtime << " ms" << endl;

	cout << "Initializing helper data structures..." << endl;
	timer.reset();
	initHelper();
	runtime = timer.reset();
	cout << "Done. It took " << runtime << " ms" << endl;
	
	if (interResultFilename=="") {
		cout << "Labeling and computing pairwise core-distance..." << endl;
		timer.reset();
		corelabeling();
		runtime = timer.reset();
		cout << "Done. It took " << runtime << " ms" << endl;
	}
	
	cout << "Clearing useless data structures..." << endl;
	timer.reset();
	if (preport!=NULL) delete[] preport; preport = NULL;
	if (mark!=NULL) delete[] mark; mark = NULL;
	if (corenodes!=NULL) delete corenodes; corenodes = NULL;
	pque = new int[gsize]; // for bidirectional core-based search
	runtime = timer.reset();
	cout << "Done. It took " << runtime << " ms" << endl;
	
	#ifdef CS_DEBUG
	cout << "-------------------------------------------" << endl;
	displayCorepaths(cout);
	cout << "-------------------------------------------" << endl;
	displayRouteTable(cout);
	cout << "-------------------------------------------" << endl;
	#endif	
}

vector<int> CoreSearch::shortestpath(int src, int trg, int& distance) {
	#ifdef CS_DEBUG
	cout << "work on shortest path from " << src << " to " << trg << endl;
	#endif

	searchspace=0; comparetimes=0; distance=0;
	vector<int> path; // last element is the length of path
	if (src==trg) {
		path.push_back(trg);
		return path;
	}
	int index=0, endindex=0, sdist=radius+1, u, val, tmp, nid;
	unsigned coreindex;
	if (src<coresize && trg<coresize) {
		#ifdef CS_DEBUG
		cout << "both are corenodes" << endl;
		#endif
		coreindex = computeMatrixIndex(src,trg);
		if (corepaths[coreindex]!=NULL) {
			path.push_back(src);
			endindex = corepaths[coreindex][0];
			for (int i=1; i<=endindex; i++) {
				path.push_back(corepaths[coreindex][i]);
			}
			path.push_back(trg);
		}
		#ifdef CS_STAT
		comparetimes++;
		#endif
		distance = path.size()-1;
		return path;
	}
	else if (src>=coresize && trg<coresize) {
		#ifdef CS_DEBUG
		cout << "only trg is corenode src=" << src << endl;
		#endif
		vector<int> trg_path;
		vector<int>::reverse_iterator rit;
		int tmp_coreid, tmp_portid, coreid, portid;
		unsigned int corepathindex;
		bool inrange = false;
		while (true) {
			if (src==trg) {
				// concatinate two segments: path+trg_path
				path.push_back(trg);
				for (rit = trg_path.rbegin(); rit != trg_path.rend(); rit++)
					path.push_back(*rit);
				distance = path.size()-1;
				return path;
			}
			inrange = false;
			for (int level=0; level<radius; level++) {
				if (level>=sdist) break;
				// check if there is no routeitems in this level
				if (routetables[src].outrouteitems==NULL) continue;
				index = routetables[src].outitemindex[level];
				endindex = routetables[src].outitemindex[level+1];
				for (; index<endindex; index++) {
					#ifdef CS_STAT
					comparetimes++;
					#endif
					tmp_coreid = routetables[src].outrouteitems[index].coreid;
					tmp_portid = routetables[src].outrouteitems[index].portid;
					if (tmp_coreid==trg) {
						inrange = true;
						sdist = level+1;
						portid = tmp_portid;
						break;
					}
					else {
						coreindex = computeMatrixIndex(tmp_coreid,trg);
						if (corepaths[coreindex]==NULL)
							continue;
						else
							val = corepaths[coreindex][0]+1; // path length is the number of intermediate nodes plus one
					}
					if (sdist>val+level+1) {
						sdist = val+level+1;
						coreid = tmp_coreid;
						portid = tmp_portid;
						corepathindex = coreindex;
					}
				}
				if (inrange) break;
			}
			// check existence of shortest path
			if (sdist>radius) {
				distance = sdist;
				return path;
			}
			// update src_path and trg_path from two directions
			if (!inrange) {
				// append path segment from intermediate coreid to trg into trg_path
				index = corepaths[corepathindex][0];
				trg_path.push_back(trg);
				for (int i = 1; i <= index; i++)
					trg_path.push_back(corepaths[corepathindex][i]);	
				trg = coreid;
			}
			// src move forward
			path.push_back(src);
			src = g.get_outneighbor(src,portid);
		}
	}
	else if (src<coresize && trg>=coresize) {
		#ifdef CS_DEBUG
		cout << "only src is corenode trg=" << trg << endl;
		#endif
		vector<int> trg_path;
		vector<int>::reverse_iterator rit;
		int tmp_coreid, tmp_portid, coreid, portid;
		unsigned int corepathindex;
		bool inrange = false;
		while (true) {
			if (src==trg) {
				// concatinate two segments: path+trg_path
				path.push_back(src);
				for (rit = trg_path.rbegin(); rit != trg_path.rend(); rit++)
					path.push_back(*rit);
				distance = path.size()-1;
				return path;
			}
			inrange = false;
			for (int level=0; level<radius; level++) {
				if (level>=sdist) break;
				// check if there is no routeitems in this level
				if (routetables[trg].inrouteitems==NULL) continue;
				index = routetables[trg].initemindex[level];
				endindex = routetables[trg].initemindex[level+1];
				for (; index<endindex; index++) {
					#ifdef CS_STAT
					comparetimes++;
					#endif
					tmp_coreid = routetables[trg].inrouteitems[index].coreid;
					tmp_portid = routetables[trg].inrouteitems[index].portid;
					if (tmp_coreid==src) {
						inrange = true;
						sdist = level+1;
						portid = tmp_portid;
						break;
					}
					else {
						coreindex = computeMatrixIndex(src,tmp_coreid);
						if (corepaths[coreindex]==NULL)
							continue;
						else
							val = corepaths[coreindex][0]+1; // path length is the number of intermediate nodes plus one
					}
					if (sdist>val+level+1) {
						sdist = val+level+1;
						coreid = tmp_coreid;
						portid = tmp_portid;
						corepathindex = coreindex;
					}
				}
				if (inrange) break;
			}
			// check existence of shortest path
			if (sdist>radius) {
				distance = sdist;
				return path;
			}
			// update src_path and trg_path from two directions
			if (!inrange) {
				index = corepaths[corepathindex][0];
				path.push_back(src);
				for (int i = 1; i <= index; i++)
					path.push_back(corepaths[corepathindex][i]);	
				src = coreid;
			}
			// trg move backward
			trg_path.push_back(trg);
			trg = g.get_inneighbor(trg,portid);
		}
	}
	else {
		#ifdef CS_STAT
		noncorequeries++;
		#endif
		#ifdef CS_DEBUG
		cout << "both are NOT corenodes" << endl;
		#endif	
		// both are not core nodes
		// find core-based shortest path and produce upper bound
		int tmp_fcoreid, tmp_bcoreid, fcoreid, fportid, bcoreid, bportid;
		unsigned int corepathindex, meetindex; 
		int findex, fendindex, bindex, bendindex, fcdist=radius+1, bcdist=radius+1;
		#ifdef CS_DEBUG
		gettimeofday(&begin_time, NULL);
		#endif
		if (routetables[src].outrouteitems!=NULL && routetables[trg].inrouteitems!=NULL) {
			for (int i = 0; i < combination.size(); i++) {
				val = combination[i][2];
				if (val>=sdist) break;
				index = combination[i][0];
				endindex = combination[i][1];
				findex = routetables[src].outitemindex[index];
				fendindex = routetables[src].outitemindex[index+1];
				if (findex==fendindex) continue;
				bindex = routetables[trg].initemindex[endindex];
				bendindex = routetables[trg].initemindex[endindex+1];
				if (bindex==bendindex) continue;
				// perform pairwise comparsion
				for (; findex < fendindex; findex++) {
					tmp_fcoreid = routetables[src].outrouteitems[findex].coreid;
					for (int j = bindex; j < bendindex; j++) {
						#ifdef CS_STAT
						comparetimes++;
						#endif
						tmp_bcoreid = routetables[trg].inrouteitems[j].coreid;
						if (tmp_fcoreid==tmp_bcoreid) {
							tmp = 0;
							meetindex = 0; // meet common corenode
						}
						else {
							meetindex = computeMatrixIndex(tmp_fcoreid,tmp_bcoreid);
							if (corepaths[meetindex]==NULL) continue;
							tmp = corepaths[meetindex][0]+1; // path length is the number of intermediate nodes plus one
						}
						if (sdist>val+tmp) {
							sdist = val+tmp; // serve as upper boun in local search without core
							fcdist = index+1;
							bcdist = endindex+1;
							fcoreid = tmp_fcoreid;
							fportid = routetables[src].outrouteitems[findex].portid;
							bcoreid = tmp_bcoreid;
							bportid = routetables[trg].inrouteitems[j].portid;
							corepathindex = meetindex;
						}
					}
				}
			}
		}
		#ifdef CS_DEBUG
		gettimeofday(&end_time, NULL);
		runtime1 += (end_time.tv_sec - begin_time.tv_sec)*1000.0 + (end_time.tv_usec - begin_time.tv_usec)*1.0/1000.0;
		cout << "core-based finding: sdist=" << sdist << " fcoreid=" << fcoreid << " fportid=" << fportid
			<< " bcoreid=" << bcoreid << " bportid=" << bportid << " corepathindex=" << corepathindex << endl;
		#endif
		
		// bidirectional local search
		int fpoint=-1, bpoint=-1; // meeting points on both sides
		if (!estimate) {
			#ifdef CS_DEBUG
			gettimeofday(&begin_time, NULL);
			#endif
			ref += radius+1;
			pque[0]=src; pque[gsize-1]=trg;
			dist[src]=ref; dist[trg]=-ref;	// forward distance positive value, backward is negative value
			prev[src]=-1; prev[trg]=-1;
			findex=0; fendindex=1; bindex=gsize-1; bendindex=bindex-1;
			bool forward = true; // forward search
			while (findex<fendindex || bindex>bendindex) {
				if (forward) {
					u=pque[findex++];
					val=dist[u];
					forall_outneighbors(g,u,eit) {
						nid=(*eit);
						if (dist[nid]<=-ref) {
							if (sdist>val+1-ref-dist[nid]-ref) {
								sdist=val+1-ref-dist[nid]-ref;
								fpoint=u;
								bpoint=nid;
							}
						}
						// update distance and queue
						if (abs(dist[nid])<ref && nid>=coresize) {
							dist[nid]=val+1;
							pque[fendindex++]=nid;
							prev[nid]=u;
						}
					}
				}
				else {
					u=pque[bindex--];
					val=dist[u];
					forall_inneighbors(g,u,eit) {
						nid=(*eit);
						if (dist[nid]>=ref) {
							if (sdist>1-val-ref+dist[nid]-ref) {
								sdist=1-val-ref+dist[nid]-ref;
								bpoint=u;
								fpoint=nid;
							}
						}
						// update distance and queue
						if (abs(dist[nid])<ref && nid>=coresize) {
							dist[nid]=val-1;
							pque[bendindex--]=nid;	
							prev[nid]=u;
						}
					}
				}
				// check stop condition
				if (findex>=fendindex || bindex<=bendindex)
					break;
				tmp = dist[pque[findex]]-dist[pque[bindex]]-ref-ref;
				if (tmp>=sdist || tmp>radius) 
					break;

				// change search direction by considering the search space
				if (bindex-bendindex>fendindex-findex) {
					if (fendindex>findex) forward=true;
					else forward=false;
				}
				else {
					if (bindex>bendindex) forward=false;
					else forward=true;
				}
			}
			searchspace = fendindex+gsize-bendindex-1-2; // not including src and trg
			#ifdef CS_DEBUG
			gettimeofday(&end_time, NULL);
			runtime2 += (end_time.tv_sec - begin_time.tv_sec)*1000.0 + (end_time.tv_usec - begin_time.tv_usec)*1.0/1000.0;	
			cout << "local search finding: fpoint=" << fpoint << " bpoint=" << bpoint << endl;
			#endif
		}
		
		// check the existence of shortest path with length no greater than radius
		if (sdist>radius) {
			distance = sdist;
			return path;
		}
		// determine whether local search or core-based search produces better results
		if (fpoint>=0) {
			// local search is better, recover shortest path based on meetpoints
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
		}
		else {
			#ifdef CS_STAT
			corehitqueries++;
			#endif
			int tmp_coreid, tmp_portid, portid;
			// recover shortest path from src to fcoreid
			path.push_back(src);
			src = g.get_outneighbor(src,fportid);
			fcdist--;
			while (src!=fcoreid) {
				for (int level=0; level<radius; level++) {
					if (level>=fcdist) break;
					// check if there is no routeitems in this level
					if (routetables[src].outrouteitems==NULL) continue;
					index = routetables[src].outitemindex[level];
					endindex = routetables[src].outitemindex[level+1];
					for (; index<endindex; index++) {
						#ifdef CS_STAT
						comparetimes++;
						#endif
						tmp_coreid = routetables[src].outrouteitems[index].coreid;
						tmp_portid = routetables[src].outrouteitems[index].portid;
						if (tmp_coreid==fcoreid) {
							portid = tmp_portid;
							break;
						}
					}
				}
				fcdist--;
				path.push_back(src);
				src = g.get_outneighbor(src,portid);
			}
			path.push_back(fcoreid);
			// recover the core shortest path from fcoreid to bcoreid
			if (corepathindex>0) {
				if (corepaths[corepathindex]!=NULL) {
					index = corepaths[corepathindex][0];
					for (int i = 1; i <= index; i++)
						path.push_back(corepaths[corepathindex][i]);
					path.push_back(bcoreid);
				}
			}
			// recover the shortest path from bcoreid to trg
			int separator = path.size();
			do {
				path.push_back(trg);
				trg = g.get_inneighbor(trg,bportid);
				bcdist--;
				for (int level=0; level<radius; level++) {
					if (level>=bcdist) break;
					// check if there is no routeitems in this level
					if (routetables[trg].inrouteitems==NULL) continue;
					index = routetables[trg].initemindex[level];
					endindex = routetables[trg].initemindex[level+1];
					for (; index<endindex; index++) {
						#ifdef CS_STAT
						comparetimes++;
						#endif
						tmp_coreid = routetables[trg].inrouteitems[index].coreid;
						tmp_portid = routetables[trg].inrouteitems[index].portid;
						if (tmp_coreid==bcoreid) {
							bportid = tmp_portid;
							break;
						}
					}
				}
			} while (trg!=bcoreid);
			reverse(path.begin()+separator,path.end());
		}
	}
	
	#ifdef CS_DEBUG
	cout << "path[" << src << "->" << trg << "]: ";
	for (int i = 0; i < path.size(); i++)
		cout << path[i] << " ";
	cout << endl;
	cout << "length=" << path.size()-1 << " searchspace=" << searchspace << endl;
	#endif
	distance = path.size()-1;
	return path;
}

vector<int> CoreSearch::shortestpath(int src, int trg) {
	int distance;
	vector<int> path = shortestpath(src,trg,distance);
	return path;
}

bool CoreSearch::test_shortestpath(int src, int trg) {
	int newdist;
	vector<int> path = shortestpath(src,trg,newdist);
	int ss;
	int distance = GraphUtil::BiBFSDist_ptr(g,src,trg,dist,pque,ref,radius,ss);
	if (path.size()-1!=distance) {
		cout << "Wrong " << src << "->" << trg << " distance=" << path.size() << " cor_dist=" << distance << endl;
		cout << "path: ";
		for (int i = 0; i < path.size(); i++)
			cout << path[i] << " ";
		cout << endl;
		exit(0);
	}
	return true;
}

int CoreSearch::getSearchspace() {
	return searchspace;
}

int CoreSearch::getComparetimes() {
	return comparetimes;
}

unsigned long CoreSearch::getNumRouteEntries() {
	return num_entries;
}

double CoreSearch::getPeakmemusage() {
	return vm_usage;
}

unsigned long CoreSearch::getCorepathsize() {
	unsigned int begin=0, end=0, len=coresize;
	unsigned long size = 0;
	for (unsigned int i = 0; i < coresize; i++) {
		for (unsigned int j = 0; j < coresize; j++) {
			begin = i*coresize+j;
			if (corepaths[begin]==NULL) continue;
			size += corepaths[begin][0];
		}
	}
	return size;
}

inline unsigned int CoreSearch::computeMatrixIndex(int i, int j) {
	return i*coresize+j;
}

void CoreSearch::loadIntermediateResults(istream& in) {
	double vm, res;
	unsigned long sep = 0l;
	num_entries = 0;
	nodemap = vector<int>(gsize,0);
	for (int i = 0; i < gsize; i++)
		in >> nodemap[i];
	in >> sep;
	assert(sep==0xFFFFFFFFFFFFFFF);
	// read corepaths
	in >> coresize;
	int size = coresize*coresize, index, val;
	corepaths = new coresizetype*[size];
	for (unsigned int i = 0; i < size; i++) {
		in >> index;
		if (index==-1) {
			corepaths[i] = NULL;
			continue;
		}
		corepaths[i] = new coresizetype[index+1];
		corepaths[i][0] = (coresizetype)index;
		for (unsigned int j = 1; j <= index; ++j) {
			in >> val;
			corepaths[i][j] = (coresizetype)val;
		}
	}	
	Util::process_mem_usage(vm, res);
	if (vm>vm_usage) vm_usage=vm;
	if (res>resident_set) resident_set=res;
	in >> sep;
	assert(sep==0xFFFFFFFFFFFFFFF);
	// read routetables
	int begin, end;
	routetables = new nodeRouteTable[gsize];
	for (int i = 0; i < coresize; i++) {
		routetables[i].outitemindex = NULL;
		routetables[i].outrouteitems = NULL;
		routetables[i].initemindex = NULL;
		routetables[i].inrouteitems = NULL;
	}
	Util::process_mem_usage(vm, res);
	if (vm>vm_usage) vm_usage=vm;
	if (res>resident_set) resident_set=res;
	for (int nid = coresize; nid < gsize; nid++) {
		in >> index;
		if (index==-1) {
			routetables[nid].outitemindex = NULL;
			routetables[nid].outrouteitems = NULL;
			continue;
		}
		routetables[nid].outitemindex = new ushort[radius+1];
		routetables[nid].outitemindex[0] = index;
		for (int i = 1; i < radius+1; i++) {
			in >> index;
			routetables[nid].outitemindex[i] = (unsigned int)index;
		}
		size = routetables[nid].outitemindex[radius];
		num_entries += size;
		routetables[nid].outrouteitems = new routeEntry[size];
		for (int i = 0; i < radius; i++) {
			begin = routetables[nid].outitemindex[i];
			end = routetables[nid].outitemindex[i+1];
			if (begin==end) continue;
			for (; begin < end; begin++) {
				in >> index;
				routetables[nid].outrouteitems[begin].coreid = (unsigned int)index;
				in >> index;
				routetables[nid].outrouteitems[begin].portid = (unsigned int)index;
			}
		}
	}
	in >> sep;
	assert(sep==0xFFFFFFFFFFFFFFF);
	for (int nid = coresize; nid < gsize; nid++) {
		in >> index;
		if (index==-1) {
			routetables[nid].initemindex = NULL;
			routetables[nid].inrouteitems = NULL;
			continue;
		}
		routetables[nid].initemindex = new ushort[radius+1];
		routetables[nid].initemindex[0] = index;
		for (int i = 1; i < radius+1; i++) {
			in >> index;
			routetables[nid].initemindex[i] = (unsigned int)index;
		}
		size = routetables[nid].initemindex[radius];
		num_entries += size;
		routetables[nid].inrouteitems = new routeEntry[size];
		for (int i = 0; i < radius; i++) {
			begin = routetables[nid].initemindex[i];
			end = routetables[nid].initemindex[i+1];
			if (begin==end) continue;
			for (; begin < end; begin++) {
				in >> index;
				routetables[nid].inrouteitems[begin].coreid = (unsigned int)index;
				in >> index;
				routetables[nid].inrouteitems[begin].portid = (unsigned int)index;
			}
		}
	}
	Util::process_mem_usage(vm, res);
	if (vm>vm_usage) vm_usage=vm;
	if (res>resident_set) resident_set=res;
}

void CoreSearch::writeIntermediateResults(ostream& out) {
	// store nodemap
	assert(nodemap.size()==gsize);
	for (int i = 0; i < nodemap.size(); i++) {
		out << nodemap[i] << " ";
	}
	out << endl;
	out << 0xFFFFFFFFFFFFFFF << endl; // separator
	out.flush();
	// store corepaths
	out << coresize << endl;
	unsigned int size = coresize*coresize, index;
	for (unsigned int i = 0; i < size; i++) {
		if (corepaths[i]==NULL)
			out << -1 << endl; // header
		else {
			index = corepaths[i][0];
			for (unsigned int j = 0; j <= index; j++)
				out << corepaths[i][j] << " ";
			out << endl;
		}
	}
	out << 0xFFFFFFFFFFFFFFF << endl; // separator
	out.flush();
	// store routetables
	int begin, end;
	for (int nid = coresize; nid < gsize; nid++) {
		if (routetables[nid].outrouteitems==NULL) {
			out << -1 << endl;
			continue;
		}
		for (int i = 0; i < radius+1; i++)
			out << routetables[nid].outitemindex[i] << " ";
		for (int i = 0; i < radius; i++) {
			begin = routetables[nid].outitemindex[i];
			end = routetables[nid].outitemindex[i+1];
			if (begin==end) continue;
			for (; begin < end; begin++) {
				out << 	routetables[nid].outrouteitems[begin].coreid << " "
					<< routetables[nid].outrouteitems[begin].portid << " ";
			}	
		}
		out << endl;
		out.flush();
	}
	out.flush();
	out << 0xFFFFFFFFFFFFFFF << endl; // separator
	for (int nid = coresize; nid < gsize; nid++) {
		if (routetables[nid].inrouteitems==NULL) {
			out << -1 << endl;
			continue;
		}
		for (int i = 0; i < radius+1; i++)
			out << routetables[nid].initemindex[i] << " ";
		for (int i = 0; i < radius; i++) {
			begin = routetables[nid].initemindex[i];
			end = routetables[nid].initemindex[i+1];
			if (begin==end) continue;
			for (; begin < end; begin++) {
				out << 	routetables[nid].inrouteitems[begin].coreid << " "
					<< routetables[nid].inrouteitems[begin].portid << " ";
			}	
		}
		out << endl;
		out.flush();
	}
	out.flush();
}

void CoreSearch::displayCorepaths(ostream& out) {
	out << "----------Core Routing Paths----------" << endl;
	unsigned int begin=0, end=0, len=coresize, size;
	for (unsigned int i = 0; i < coresize; i++) {
		out << "corenode[" << i << "]: ";
		for (unsigned int j = 0; j < coresize; j++) {
			begin = i*coresize+j;
			if (corepaths[begin]==NULL) continue;
			out << "(" << j << ":";
			size = corepaths[begin][0];
			for (int k = 1; k <= size; k++)
				out << " " << corepaths[begin][k];
			out << ") ";
		}
		out << endl;
	}	
}

void CoreSearch::displayRouteTable(ostream& out) {
	out << "----------Core-based Routing Table----------" << endl;
	for (int i = coresize; i < gsize; i++) {
		displayRouteTableByNode(i,out);
	}
}

void CoreSearch::displayRouteTableByNode(int nid, ostream& out) {
	int begin, end;
	out << "OutRouteTable[" << nid << "]" << endl;
	if (routetables[nid].outrouteitems==NULL) {
		out << "Empty" << endl;
	}
	else {
		for (int i = 0; i < radius; i++) {
			begin = routetables[nid].outitemindex[i];
			end = routetables[nid].outitemindex[i+1];
			if (begin==end) continue;
			out << "level " << i+1 << ": ";
			for (; begin < end; begin++) {
				out << 	"[" << routetables[nid].outrouteitems[begin].coreid << ","
					<< routetables[nid].outrouteitems[begin].portid << "] ";
			}
			out << endl;
		}
	}
	out << "InRouteTable[" << nid << "]" << endl;
	if (routetables[nid].inrouteitems==NULL) {
		out << "Empty" << endl;
	}
	else {
		for (int i = 0; i < radius; i++) {
			begin = routetables[nid].initemindex[i];
			end = routetables[nid].initemindex[i+1];
			if (begin==end) continue;
			out << "level " << i+1 << ": ";
			for (; begin < end; begin++) {
				out << 	"[" << routetables[nid].inrouteitems[begin].coreid << ","
					<< routetables[nid].inrouteitems[begin].portid << "] ";
			}
			out << endl;
		}
	}
}

int CoreSearch::distance(int src, int trg) {
	#ifdef CS_DEBUG
	cout << "work on distance from " << src << " to " << trg << endl;
	#endif

	searchspace=0; comparetimes=0;
	if (src==trg) return 0;
	int index=0, endindex=0, sdist=radius+1, u, val, tmp, nid;
	unsigned int coreindex;
	if (src<coresize && trg<coresize) {
		#ifdef CS_DEBUG
		cout << "both are corenodes" << endl;
		#endif
		coreindex = computeMatrixIndex(src,trg);
		if (corepaths[coreindex]!=NULL) {
			endindex = corepaths[coreindex][0];
			return endindex+1;
		}
		#ifdef CS_STAT
		comparetimes++;
		#endif
		return -1;
	}
	else if (src>=coresize && trg<coresize) {
		#ifdef CS_DEBUG
		cout << "only trg is corenode src=" << src << endl;
		#endif
		int tmp_coreid;
		bool inrange = false;
		for (int level=0; level<radius; level++) {
			if (level>=sdist) break;
			// check if there is no routeitems in this level
			if (routetables[src].outrouteitems==NULL) continue;
			index = routetables[src].outitemindex[level];
			endindex = routetables[src].outitemindex[level+1];
			for (; index<endindex; index++) {
				#ifdef CS_STAT
				comparetimes++;
				#endif
				tmp_coreid = routetables[src].outrouteitems[index].coreid;
				if (tmp_coreid==trg) {
					inrange = true;
					sdist = level+1;
					break;
				}
				else {
					coreindex = computeMatrixIndex(tmp_coreid,trg);
					if (corepaths[coreindex]==NULL)
						continue;
					else
						val = corepaths[coreindex][0]+1; // path length is the number of intermediate nodes plus one
				}
				if (sdist>val+level+1) {
					sdist = val+level+1;
				}
			}
			if (inrange) break;
		}
		return sdist>radius?-1:sdist;
	}
	else if (src<coresize&&trg>=coresize) {
		#ifdef CS_DEBUG
		cout << "only src is corenode trg=" << trg << endl;
		#endif
		int tmp_coreid;
		bool inrange = false;
		for (int level=0; level<radius; level++) {
			if (level>=sdist) break;
			// check if there is no routeitems in this level
			if (routetables[trg].inrouteitems==NULL) continue;
			index = routetables[trg].initemindex[level];
			endindex = routetables[trg].initemindex[level+1];
			for (; index<endindex; index++) {
				#ifdef CS_STAT
				comparetimes++;
				#endif
				tmp_coreid = routetables[trg].inrouteitems[index].coreid;
				if (tmp_coreid==src) {
					inrange = true;
					sdist = level+1;
					break;
				}
				else {
					coreindex = computeMatrixIndex(src,tmp_coreid);
					if (corepaths[coreindex]==NULL)
						continue;
					else
						val = corepaths[coreindex][0]+1; // path length is the number of intermediate nodes plus one
				}
				if (sdist>val+level+1) {
					sdist = val+level+1;
				}
			}
			if (inrange) break;
		}
		return sdist>radius?-1:sdist;
	}
	else {
		#ifdef CS_STAT
		noncorequeries++;
		#endif
		#ifdef CS_DEBUG
		cout << "both are NOT corenodes" << endl;
		int fcoreid, fportid, bcoreid, bportid;
		unsigned int corepathindex;
		#endif
		// both are not core nodes
		// find core-based shortest path and produce upper bound
		int tmp_fcoreid, tmp_bcoreid, findex, fendindex, bindex, bendindex, psdist;
		#ifdef CS_DEBUG
		gettimeofday(&begin_time, NULL);
		#endif
		if (routetables[src].outrouteitems!=NULL && routetables[trg].inrouteitems!=NULL) {
			for (int i = 0; i < combination.size(); i++) {
				val = combination[i][2];
				if (val>=sdist) break;
				index = combination[i][0];
				endindex = combination[i][1];
				findex = routetables[src].outitemindex[index];
				fendindex = routetables[src].outitemindex[index+1];
				if (findex==fendindex) continue;
				bindex = routetables[trg].initemindex[endindex];
				bendindex = routetables[trg].initemindex[endindex+1];
				if (bindex==bendindex) continue;
				// perform pairwise comparsion
				for (; findex < fendindex; findex++) {
					tmp_fcoreid = routetables[src].outrouteitems[findex].coreid;
					for (int j = bindex; j < bendindex; j++) {
						#ifdef CS_STAT
						comparetimes++;
						#endif
						tmp_bcoreid = routetables[trg].inrouteitems[j].coreid;
						if (tmp_fcoreid==tmp_bcoreid) {
							tmp = 0;
							index = -1; // meet common corenode
						}
						else {
							coreindex = computeMatrixIndex(tmp_fcoreid,tmp_bcoreid);
							if (corepaths[coreindex]==NULL) continue;
							tmp = corepaths[coreindex][0]+1; // path length is the number of intermediate nodes plus one
						}
						if (sdist>val+tmp) {
							sdist = val+tmp; // serve as upper boun in local search without core
							#ifdef CS_DEBUG
							fcoreid = tmp_fcoreid;
							fportid = routetables[src].outrouteitems[findex].portid;
							bcoreid = tmp_bcoreid;
							bportid = routetables[trg].inrouteitems[j].portid;
							corepathindex = coreindex;
							#endif
						}
					}
				}
			}
		}
		#ifdef CS_STAT
		psdist = sdist;
		#endif
		#ifdef CS_DEBUG
		gettimeofday(&end_time, NULL);
		runtime1 += (end_time.tv_sec - begin_time.tv_sec)*1000.0 + (end_time.tv_usec - begin_time.tv_usec)*1.0/1000.0;
		cout << "core-based finding: sdist=" << sdist << " fcoreid=" << fcoreid << " fportid=" << fportid
			<< " bcoreid=" << bcoreid << " bportid=" << bportid << " corepathindex=" << corepathindex << endl;
		#endif
		
		// bidirectional local search
		if (!estimate) {
			#ifdef CS_DEBUG
			int fpoint=-1, bpoint=-1; // meeting points on both sides
			gettimeofday(&begin_time, NULL);
			#endif
			ref += radius+1;
			pque[0]=src; pque[gsize-1]=trg;
			dist[src]=ref; dist[trg]=-ref;	// forward distance positive value, backward is negative value
			findex=0; fendindex=1; bindex=gsize-1; bendindex=bindex-1;
			bool forward = true; // forward search
			while (findex<fendindex || bindex>bendindex) {
				if (forward) {
					u = pque[findex++];
					val = dist[u];
					forall_outneighbors(g,u,eit) {
						nid=(*eit);
						if (dist[nid]<=-ref) {
							if (sdist>val+1-ref-dist[nid]-ref) {
								sdist=val+1-ref-dist[nid]-ref;
								#ifdef CS_DEBUG
								fpoint = u;
								bpoint = nid;
								#endif
							}
						}
						// update distance and queue
						if (abs(dist[nid])<ref && nid>=coresize) {
							dist[nid] = val+1;
							pque[fendindex++] = nid;
						}
					}
				}
				else {
					u = pque[bindex--];
					val = dist[u];
					forall_inneighbors(g,u,eit) {
						nid=(*eit);
						if (dist[nid]>=ref) {
							if (sdist>1-val-ref+dist[nid]-ref) {
								sdist=1-val-ref+dist[nid]-ref;
								#ifdef CS_DEBUG
								bpoint = u;
								fpoint = nid;
								#endif
							}
						}
						// update distance and queue
						if (abs(dist[nid])<ref && nid>=coresize) {
							dist[nid] = val-1;
							pque[bendindex--] = nid;
						}
					}
				}
				// check stop condition
				if (findex>=fendindex || bindex<=bendindex)
					break;
				tmp = dist[pque[findex]]-dist[pque[bindex]]-ref-ref;
				if (tmp>=sdist-1 || tmp>radius) 
					break;
				// change search direction by considering the search space
				if (bindex-bendindex>fendindex-findex) {
					if (fendindex>findex) forward=true;
					else forward=false;
				}
				else {
					if (bindex>bendindex) forward=false;
					else forward=true;
				}
			}
			searchspace = fendindex+gsize-bendindex-1-2; // not including src and trg
			#ifdef CS_DEBUG
			gettimeofday(&end_time, NULL);
			runtime2 += (end_time.tv_sec - begin_time.tv_sec)*1000.0 + (end_time.tv_usec - begin_time.tv_usec)*1.0/1000.0;
			cout << "local search finding: fpoint=" << fpoint << " bpoint=" << bpoint << endl;
			#endif
		}
		#ifdef CS_STAT
		if (psdist==sdist) corehitqueries++;
		#endif
	}
	return sdist>radius?-1:sdist;
}

bool CoreSearch::test_distance(int src, int trg) {
	int ss;
	int cdist = distance(src,trg);
	int distance = GraphUtil::BiBFSDist_ptr(g,src,trg,dist,pque,ref,radius,ss);
	if (cdist!=distance) {
		cout << "Wrong " << src << "->" << trg << " distance=" << cdist << " correct_dist=" << distance << endl;
		exit(0);
	}
	return true;
}
