#include "Graph.h"

struct Comparison {
	bool operator() ( const Edge& a, const Edge& b ) const {
		return ( a.target < b.target );
	}
}myComparison;

struct LabelComp{
	bool operator() ( const Label& a, const Label& b) const {
		return (a.distance < b.distance);
	}
}myLabelComp;

Graph::Graph() {
	nodesize = 0;
}

Graph::Graph(int size) {
	nodesize = size;
	nodelist = VertexList(size,Vertex());
	for (int index = 0; index < size; index++) {
		nodelist[index].rank = UINT_INFINITY;
		nodelist[index].flaga = false;
		nodelist[index].flagb = false;
	}
	inedgeindex = vector<int>(nodesize+1,0);
	outedgeindex = vector<int>(nodesize+1,0);
	
	inShortCutIndex = vector<int>(0);
	outShortCutIndex = vector<int>(0);
	
	outLabelList = vector<LabelList>(nodesize, LabelList(0));
	inLabelList = vector<LabelList>(nodesize, LabelList(0));	
	
	pathIndex = vector<vector<int> >(nodesize, vector<int>(0));
	pathList = vector<vector<VertexID> >(nodesize, vector<EdgeID>(0));	
}

Graph::Graph(istream& in) {
	readGraph(in);
	sortEdges();
	
	//cout << "Outindex:" << outedgeindex.size() << endl;
	nodesize = outedgeindex.size()-1;
	nodelist = VertexList(nodesize,Vertex());
	for (int index = 0; index < nodesize; index++) {
		nodelist[index].rank = UINT_INFINITY;
		nodelist[index].flaga = false;
		nodelist[index].flagb = false;
	}
	
	inShortCutIndex = vector<int>(0);
	outShortCutIndex = vector<int>(0);	
	
	outLabelList = vector<LabelList>(nodesize, LabelList(0));
	inLabelList = vector<LabelList>(nodesize, LabelList(0));	
	
	pathIndex = vector<vector<int> >(nodesize, vector<int>(0));
	pathList = vector<vector<VertexID> >(nodesize, vector<EdgeID>(0));
} 

// Graph::Graph(istream& in, bool isedgelist) {
	// if (isedgelist) readEdgeList(in);
	// else readGraph(in);
	//sort edges
	// sortEdges();
// }

Graph::~Graph() {}

void Graph::printGraph() {
	//writeGraph(cout);
	cout << "Incoming edgelist:" << endl;
	int begin, end;
	for (int i = 0; i < nodesize; i++) {
		cout << i << ": ";
		begin = inedgeindex[i];
		end = inedgeindex[i+1];
		for (; begin < end; begin++)
			cout << "[" << inedgelist[begin].target << ", " << inedgelist[begin].weight << "] ";
		cout << "#" << endl;
	}
}

void Graph::printShortcut() {
	cout << "Outcoming edgelist:" << endl;
	int begin, end;
	for ( int i = 0; i < nodesize; i++ ){
		cout << i << ": ";
		begin = outShortCutIndex[i];
		end = outShortCutIndex[i+1];
		for (; begin < end; begin++)
			cout << "[" << outShortCutList[begin].target << ", " << outShortCutList[begin].weight << "]";
		cout << "#" << endl;
	}
	
	#ifdef GRAPH_DEBUG	
		for (int i = 0; i < outShortCutIndex.size(); i++) {
			cout << outShortCutIndex[i] << " ";
		}
		cout << endl;
	#endif	
}

void Graph::printRank(){
	cout << "Print node rank:" << endl;
	for ( int i = 0; i < nodelist.size(); i++ ) {
		cout << "node " << i << ": rank " << nodelist[i].rank << endl;
	}
}

void Graph::sortEdges() {
	int begin, end;
	for (int i = 0; i < nodesize; i++) {
		begin = inedgeindex[i];
		end = inedgeindex[i+1];
		sort(inedgelist.begin()+begin,inedgelist.begin()+end, myComparison);
		begin = outedgeindex[i];
		end = outedgeindex[i+1];
		sort(outedgelist.begin()+begin,outedgelist.begin()+end, myComparison);
	}
}

void Graph::clear() {
	nodesize = 0;
	//nodelist.clear();
	inedgelist.clear();
	outedgelist.clear();
	inedgeindex.clear();
	outedgeindex.clear();
}

void Graph::readGraph(istream& fs) {
    string str;
    char* line;
	
	// read the plain graph;
	vector<vector<double> > tmpGraph;
    while (!fs.eof()) {
        getline(fs, str);
		if (str.size()==0) continue; // remove the empty final line;
        line = new char[str.size()+1];
        strcpy(line, str.c_str());
        char* num = strtok(line, " #\t:,;");
        vector<double> tv;
        while (num != NULL) {
            tv.push_back(atof(num));
            num = strtok(NULL, " #\t:,;");
        }
        tmpGraph.push_back(tv);
        delete[] line;
    }
	
	// remove the multiple edges and self loops;
	for (int i = 1; i < tmpGraph.size(); i++) {
		set<int> pool;
		set<int>::iterator it;
		for (int j = 1; j < tmpGraph[i].size();) {
			int v = static_cast<int>(tmpGraph[i][j]);
			it = pool.find(v);
			if (it != pool.end() || v == i-1){
				tmpGraph[i].erase(tmpGraph[i].begin()+j);
				tmpGraph[i].erase(tmpGraph[i].begin()+j);
			}else {
				pool.insert(v);
				j+=2;
			}
		}
	}
	
	
	// for (int i = 0; i < tmpGraph.size(); i++){
		// for(int j = 0; j < tmpGraph[i].size(); j++){
			// cout << tmpGraph[i][j] << " ";
		// }
		// cout << endl;
	// }
	
	// materialize outedgeindex and outedgelist;
	nodesize = static_cast<unsigned int>(tmpGraph[0][0]);
	//cout << "nodesize = " << nodesize << endl;
	for (int i = 1; i < tmpGraph.size(); i++){
		EdgeList tmpOutEdgeList;
		for (int j = 1; j < tmpGraph[i].size(); j=j+2) {
			Edge tmpVertex;
			tmpVertex.source = static_cast<VertexID> (tmpGraph[i][0]);
			tmpVertex.target = static_cast<VertexID> (tmpGraph[i][j]);
			tmpVertex.weight = static_cast<Weight> (tmpGraph[i][j+1]);
			tmpVertex.flaga = false;
			tmpVertex.flagb = false;
			tmpVertex.flagc = false;
			tmpVertex.flagd = false;
			tmpOutEdgeList.push_back(tmpVertex);
		}
		outedgeindex.push_back(outedgelist.size());
		outedgelist.insert(outedgelist.end(), tmpOutEdgeList.begin(), tmpOutEdgeList.end());
	}
	outedgeindex.push_back(outedgelist.size());

	// build inverse graph;
	vector<vector<double> > tmpInverseGraph(tmpGraph.size()-1, vector<double>());
	for( int i = 1; i < tmpGraph.size(); i++ ) {
		for( int j = 1; j < tmpGraph[i].size(); j=j+2 ) {
			VertexID target = static_cast<VertexID>(tmpGraph[i][j]);
			tmpInverseGraph[target].push_back(tmpGraph[i][0]);
			tmpInverseGraph[target].push_back(tmpGraph[i][j+1]);
		}
	}
	vector<vector<double> >(0).swap(tmpGraph);
	
	// for (int i = 0; i < tmpInverseGraph.size(); i++){
		// cout << i << ": ";
		// for(int j = 0; j < tmpInverseGraph[i].size(); j++){
			// cout << tmpInverseGraph[i][j] << " ";
		// }
		// cout << endl;
	// }
	
	// materialize inedgeindex and inedgelist;
	for (int i = 0; i < tmpInverseGraph.size(); i++ ){
		EdgeList tmpInEdgeList;
		for (int j = 0; j < tmpInverseGraph[i].size(); j=j+2 ) {
			Edge tmpVertex;
			tmpVertex.source = i;
			tmpVertex.target = static_cast<VertexID> (tmpInverseGraph[i][j]);
			tmpVertex.weight = static_cast<Weight> ( tmpInverseGraph[i][j+1] );
			tmpVertex.flaga = false;
			tmpVertex.flagb = false;
			tmpVertex.flagc = false;
			tmpVertex.flagd = false;
			tmpInEdgeList.push_back(tmpVertex);
		}
		inedgeindex.push_back(inedgelist.size());
		inedgelist.insert(inedgelist.end(), tmpInEdgeList.begin(), tmpInEdgeList.end());
	}
	inedgeindex.push_back(inedgelist.size());	
}	

// void Graph::readEdgeList(istream& in) {
	// int src, trg;
	// int maxid = -1;
	// vector<vector<int> > tmpoutg = vector<vector<int> >(128,vector<int>());
	// vector<vector<int> > tmping = vector<vector<int> >(128,vector<int>());
	// while (!in.eof()) {
		// in >> src >> trg;
		// if (src>maxid) {
			// maxid = src;
			// if (maxid>=tmpoutg.size()) {
				// tmpoutg.resize(maxid+128,vector<int>());
				// tmping.resize(maxid+128,vector<int>());
			// }
		// }
		// tmpoutg[src].push_back(trg);
		// tmping[trg].push_back(src);
	// }
	// nodesize = maxid+1;
	// vector<int>::iterator vit;
	// for (int i = 0; i < nodesize; i++) {
		// outedgeindex.push_back(outedgelist.size());
		// for (vit = tmpoutg[i].begin(); vit != tmpoutg[i].end(); vit++) 
			// outedgelist.push_back(*vit);
		// tmpoutg[i].clear();
		// inedgeindex.push_back(inedgelist.size());
		// for (vit = tmping[i].begin(); vit != tmping[i].end(); vit++) 
			// inedgelist.push_back(*vit);
		// tmping[i].clear();
	// }
	// outedgeindex.push_back(outedgelist.size());
	// inedgeindex.push_back(inedgelist.size());
	// tmpoutg.clear();
	// nodelist = VertexList(nodesize,Vertex());
	// for (int i = 0; i < nodesize; i++) {
		//nodelist[i].order = i;
		// nodelist[i].flag = false;
	// }
// }

void Graph::writeGraph(ostream& out) {
	cout << "#Vertex=" << nodesize << "\t#Edge="<< outedgelist.size() << endl;
	out << "graph_for_distance" << endl;
	out << nodesize << endl;
	int begin, end;
	for (int i = 0; i < nodesize; i++) {
		out << i << ": ";
		begin = outedgeindex[i];
		end = outedgeindex[i+1];
		for (; begin < end; begin++)
			out << outedgelist[begin].target << " " << outedgelist[begin].weight << " ";
		out << "#" << endl;
	}
}

vector<string>& Graph::split(const string &s, char delim, vector<string> &elems) {
	stringstream ss(s);
	string item;
	while(getline(ss, item, delim)) {
		elems.push_back(item);
	}
	return elems;
}

vector<string> Graph::split(const string &s, char delim) {
	vector<string> elems;
	return split(s, delim, elems);
}

int Graph::num_vertices() {
	return outedgeindex.size()-1;
}

int Graph::num_edges() {
	return outedgelist.size();
}

int Graph::out_degree(int src) {
	return outedgeindex[src+1]-outedgeindex[src];
}

int Graph::in_degree(int src) {
	return inedgeindex[src+1]-inedgeindex[src];
}

// check whether the edge from src to trg is in the graph
bool Graph::hasEdge(int src, int trg) {
	int begin = outedgeindex[src];
	int end = outedgeindex[src+1];
	for (; begin < end; begin++)
		if (outedgelist[begin].target==trg) 
			return true;
	return false;
}

// return vertex list of graph
VertexList& Graph::vertices() {
	return nodelist;
}

Graph& Graph::operator=(const Graph& g) {
	if (this != &g) {
		//nodelist = g.nodelist;
		outedgelist = g.outedgelist;
		outedgeindex = g.outedgeindex;
		inedgelist = g.inedgelist;
		inedgeindex = g.inedgeindex;
		//inportlist = g.inportlist;
		//outportlist = g.outportlist;
		nodesize = g.nodesize;
	}
	return *this;
}

// get a specified vertex property
Vertex& Graph::operator[](const int vid) {
	assert(vid>=0&&vid<nodesize);
	return nodelist[vid];
}

// void Graph::reorderid(const vector<int>& nodemap) {
	// update nodelist
	// VertexList newnodelist = VertexList(nodesize,Vertex());
	// for (int i = 0; i < nodesize; i++) {
		// newnodelist[nodemap[i]] = nodelist[i];
	// }
	// nodelist.clear();
	// nodelist = newnodelist;
	// newnodelist.clear();
	
	// update edgelist and edgeindex
	// vector<vector<int> > newoutedgelist = vector<vector<int> >(nodesize,vector<int>());
	// vector<vector<int> > newinedgelist = vector<vector<int> >(nodesize,vector<int>());
	// int begin, end;
	// for (int i = 0; i < nodesize; i++) {
		// begin = outedgeindex[i];
		// end = outedgeindex[i+1];
		// for (; begin < end; begin++) {
			// newoutedgelist[nodemap[i]].push_back(nodemap[outedgelist[begin]]);
		// }
		// sort(newoutedgelist[nodemap[i]].begin(),newoutedgelist[nodemap[i]].end());
		// begin = inedgeindex[i];
		// end = inedgeindex[i+1];
		// for (; begin < end; begin++) {
			// newinedgelist[nodemap[i]].push_back(nodemap[inedgelist[begin]]);
		// }
		// sort(newinedgelist[nodemap[i]].begin(),newinedgelist[nodemap[i]].end());
	// }
	// outedgelist.clear();
	// outedgeindex.clear();
	// inedgelist.clear();
	// inedgeindex.clear();
	// for (int i = 0; i < nodesize; i++) {
		// outedgeindex.push_back(outedgelist.size());
		// for (int j = 0; j < newoutedgelist[i].size(); j++)
			// outedgelist.push_back(newoutedgelist[i][j]);
		// newoutedgelist[i].clear();
		// inedgeindex.push_back(inedgelist.size());
		// for (int j = 0; j < newinedgelist[i].size(); j++)
			// inedgelist.push_back(newinedgelist[i][j]);
		// newinedgelist[i].clear();
	// }
	// outedgeindex.push_back(outedgelist.size());
	// inedgeindex.push_back(inedgelist.size());
	// newoutedgelist.clear();
	// newinedgelist.clear();
// }

// for core-based routing only
// void Graph::constructPortlist(int coreindex) {
	// inportlist.clear();
	// outportlist.clear();
	// int begin, end, index, port, nid;
	// for (int i = 0; i < nodesize; i++) {
		// build inportlist
		// begin = outedgeindex[i];
		// end = outedgeindex[i+1];
		// for (index = begin; index < end; index++) {
			// nid = outedgelist[index];
			// if (nid>=coreindex) {
				// port = Util::binarysearch_index(inedgelist, inedgeindex[nid], inedgeindex[nid+1], i);
				// port -= inedgeindex[nid];
				// assert(port>-1);
				// inportlist.push_back((degreetype)port);
			// }
			// else 
				// inportlist.push_back(DEGREE_MAX);
		// }
		// build outportlist
		// begin = inedgeindex[i];
		// end = inedgeindex[i+1];
		// for (index = begin; index < end; index++) {
			// nid = inedgelist[index];
			// if (nid>=coreindex) {
				// port = Util::binarysearch_index(outedgelist, outedgeindex[nid], outedgeindex[nid+1], i);
				// port -= outedgeindex[nid];
				// assert(port>-1);
				// outportlist.push_back((degreetype)port);
			// }
			// else 
				// outportlist.push_back(DEGREE_MAX);
		// }
	// }
	
	// #ifdef GRAPH_DEBUG
	// cout << "inportlist coreindex=" << coreindex << endl;
	// int __begin, __end;
	// for (int i = 0; i < outedgeindex.size()-1; i++) {
		// cout << i << ": ";
		// __begin = outedgeindex[i];
		// __end = outedgeindex[i+1];
		// for (int j = __begin; j < __end; j++) {
			// cout << inportlist[j] << " ";
		// }
		// cout << endl;
	// }
	// cout << "outportlist coreindex=" << coreindex << endl;
	// for (int i = 0; i < inedgeindex.size()-1; i++) {
		// cout << i << ": ";
		// __begin = inedgeindex[i];
		// __end = inedgeindex[i+1];
		// for (int j = __begin; j < __end; j++) {
			// cout << outportlist[j] << " ";
		// }
		// cout << endl;
	// }
	// #endif
// }

// void Graph::clearPortlist() {
	// inportlist.clear();
	// outportlist.clear();
// }

// void Graph::writePortlist(ostream& out) {
	// out << "Outportlist_for_search" << endl;
	// int begin, end;
	// for (int i = 0; i < nodesize; i++) {
		// out << i << ": ";
		// begin = outedgeindex[i];
		// end = outedgeindex[i+1];
		// for (; begin < end; begin++)
			// out << inportlist[begin] << " ";
		// out << "#" << endl;
	// }
	// out << "Inportlist_for_search" << endl;
	// for (int i = 0; i < nodesize; i++) {
		// out << i << ": ";
		// begin = inedgeindex[i];
		// end = inedgeindex[i+1];
		// for (; begin < end; begin++)
			// out << outportlist[begin] << " ";
		// out << "#" << endl;
	// }	
// }

// void Graph::printPortlist() {
	// writePortlist(cout);
// }

bool Graph::isDirectedGraph() {
	for (int i = 0; i < inedgeindex.size(); i++) {
		if (inedgeindex[i]!=outedgeindex[i])
			return true;
	}
	return false;
}


void Graph::insertShortcut(vector<ShortCuts>& tmpinshortcut, vector<ShortCuts>& tmpoutshortcut){
	inShortCutIndex.push_back(0);
	for (int i = 0; i < tmpinshortcut.size(); i++){
		for (int j = 0; j < tmpinshortcut[i].size(); j++) {
			inShortCutList.push_back(tmpinshortcut[i][j]);
		}
		//ShortCuts(0).swap(tmpinshortcut[i]);
		inShortCutIndex.push_back(inShortCutList.size());
	}
	
	#ifdef GRAPH_DEBUG 
		cout << "inShortCut:" << endl;
		for ( int i = 0; i < inShortCutIndex.size(); i++ ) {
			cout << inShortCutIndex[i] << " ";
		}
		cout << endl;
	#endif
	
	outShortCutIndex.push_back(0);
	for ( int i = 0; i < tmpoutshortcut.size(); i++ ) {
		for ( int j = 0; j < tmpoutshortcut[i].size(); j++ ) {
			outShortCutList.push_back(tmpoutshortcut[i][j]);
		}
		//ShortCuts(0).swap(tmpoutshortcut[i]);
		outShortCutIndex.push_back(outShortCutList.size());
	}
	
	#ifdef GRAPH_DEBUG 	
		cout << "outShortCut:" << endl;
		for ( int i = 0; i < outShortCutIndex.size(); i++ ) {
			cout << outShortCutIndex[i] << " ";
		}
		cout << endl;
	#endif	
}

void Graph::insertNodeList(VertexList& nodelist) {
	VertexList (0).swap(this->nodelist);
	for ( int i = 0; i < nodelist.size(); i++ ){
		(this->nodelist).push_back(nodelist[i]);	
	}
}

void Graph::insertLabelList(vector<LabelList>& inlabellist, vector<LabelList>& outlabellist){
	vector<LabelList>(0).swap(this->outLabelList);
	vector<LabelList>(0).swap(this->inLabelList);
	
	this->outLabelList.resize(nodesize);
	this->inLabelList.resize(nodesize);
	for ( int i = 0; i < nodesize; i++ ) {
		// for ( int j = 0; j < inlabellist[i].size(); j++ ) {
			// inLabelList[i].push_back( inlabellist[i][j] );
		// }
		for ( int j = 0; j < outlabellist[i].size(); j++ ) {
			outLabelList[i].push_back( outlabellist[i][j] );
		}		
	}
}

LabelList& Graph::exportOutLabel(VertexID source){
	return outLabelList[source];
}
LabelList& Graph::exportInLabel(VertexID target){
	return inLabelList[target];
}

// for each node, sort its label in order of accending weight;
void Graph::sortLabel(){
	for ( int i = 0; i < outLabelList.size(); i++ ) {
		sort(outLabelList[i].begin(), outLabelList[i].end(), myLabelComp);
	}
	for ( int i = 0; i < inLabelList.size(); i++ ) {
		sort(inLabelList[i].begin(), inLabelList[i].end(), myLabelComp);
	}	
}

void Graph::insertInLabel(vector<vector<double> >& data){
	
	// cout << "Data::" << endl;
	// for ( int i = 0; i < data.size(); i++ ) {
		// for ( int j = 0; j < data[i].size(); j++ ) {
			// cout << data[i][j] << " ";
		// }
		// cout << endl;
	// }
	
	
	vector<LabelList>(0).swap(inLabelList);
	
	inLabelList.resize(data.size());
	for ( int i = 0; i < data.size(); i++ ) {
		for ( int j = 1; j < data[i].size(); j+=2 ) {
			Label tmp_label;
			tmp_label.id = static_cast<VertexID>(data[i][j]);
			tmp_label.distance = data[i][j+1];
			inLabelList[i].push_back(tmp_label);
		}
	}
	
	// cout << "in label list: " << endl;
	// for ( int i = 0; i < inLabelList.size(); i++ ) {
		// cout << i << ": ";
		// for ( int j = 0; j < inLabelList[i].size(); j++ ) {
			// cout << "[" << inLabelList[i][j].id << ", " << inLabelList[i][j].distance << "] ";
		// }
		// cout << endl;
	// }	
}

void Graph::insertOutLabel(vector<vector<double> >& data){
	// cout << "Data::" << endl;
	// for ( int i = 0; i < data.size(); i++ ) {
		// for ( int j = 0; j < data[i].size(); j++ ) {
			// cout << data[i][j] << " ";
		// }
		// cout << endl;
	// }


	vector<LabelList>(0).swap(outLabelList);
	
	outLabelList.resize(data.size());
	for ( int i = 0; i < data.size(); i++ ) {
		for ( int j = 1; j < data[i].size(); j+=2 ) {
			Label tmp_label;
			tmp_label.id = static_cast<VertexID>(data[i][j]);
			tmp_label.distance = data[i][j+1];
			outLabelList[i].push_back(tmp_label);
		}
	}
}

vector<vector<int> >& Graph::exportPathIndex(void){
	return pathIndex;
}

vector<vector<VertexID> >& Graph::exportPathList(void){
	return pathList;
}

void Graph::insertPathIndex(vector<vector<double> >& data){
	// cout << "Data::" << endl;
	// for ( int i = 0; i < data.size(); i++ ) {
		// for ( int j = 0; j < data[i].size(); j++ ) {
			// cout << data[i][j] << " ";
		// }
		// cout << endl;
	// }


	vector<vector<int> >(0).swap(pathIndex);
	
	pathIndex.resize(data.size());
	for ( int i = 0; i < data.size(); i++ ) {
		for ( int j = 1; j < data[i].size(); j++ ) {
			pathIndex[i].push_back( static_cast<int>(data[i][j]) );
		}
	}
}

void Graph::insertPathList(vector<vector<double> >& data){
	// cout << "Data::" << endl;
	// for ( int i = 0; i < data.size(); i++ ) {
		// for ( int j = 0; j < data[i].size(); j++ ) {
			// cout << data[i][j] << " ";
		// }
		// cout << endl;
	// }


	vector<vector<VertexID> >(0).swap(pathList);
	
	pathList.resize(data.size());
	for ( int i = 0; i < data.size(); i++ ) {
		for ( int j = 1; j < data[i].size(); j++ ) {
			pathList[i].push_back( static_cast<int>(data[i][j]) );
		}
	}
}

ShortCuts& Graph::exportinShortCutList() {
	return inShortCutList;
}
		
ShortCuts& Graph::exportoutShortCutList(){
	return outShortCutList;
}
		
vector<int>& Graph::exportinShortCutIndex(){
	return inShortCutIndex;
}

vector<int>& Graph::exportoutShortCutIndex(){
	return outShortCutIndex;
}




