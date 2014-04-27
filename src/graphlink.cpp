#include "graphlink.h"

void printUsage() {
	cerr << "Usage of CH:" << endl;
	cerr << "  Note that this program now has two stages: " << endl;
	cerr << "    S1: generate the hierarchical structures for all nodes." << endl;
	cerr << "    	COMMAND: ./run --filename=data.gra --generate_label=yes" << endl;
	cerr << "    S2: query the graph." << endl;
	cerr << "    	COMMAND: ./run --filename=data.gra --generate_label=no" << endl << endl;

	/*
	cerr << "Command line parameters:" << endl << endl;

	cerr << "--help" << endl;
	cerr << "--usage" << endl;
	cerr << "  prints this usage information and exits." << endl << endl;

	cerr << "--query_num=1000" << endl;
	cerr << "  specifies the number of random shortest path queries to be 1000." << endl << endl;

	// cerr << "--resultfile=result.txt" << endl;
	// cerr << "  Specifies the file to output result." << endl << endl;
	
	cerr << "--filename=data.gra" << endl;
	cerr << "  use data.gra as input file." << endl << endl;
	
	cerr << "--unpack=yes[default:no]" << endl;
	cerr << "  print out the path." << endl << endl;	
	
	cerr << "--experiment=1[default:0]" << endl;
	cerr << "  choose the experiment: 1: correctness; 2: performance; 0: both." << endl;
	cerr << "  this paramter is for debugging only." << endl << endl;		
	
	cerr << "--save_para=no[default:yes]" << endl;
	cerr << "  true: use the parameters of last experiment if they are saved before." << endl; 
	cerr << "  false: generate new parameters for this experiment." << endl << endl;

	cerr << "--method=all[default:dijk]" << endl;
	cerr << "  dijk: only run the dijkstra method." << endl; 
	cerr << "  label: only run labeling method." << endl;
	cerr << "  all: run two methods." << endl;
	cerr << "  this parameter is for debugging only. " << endl << endl;
	
	cerr << "--idmap=yes[default:no]" << endl;
	cerr << "  yes: will read the .idmap file from the same directory of graph file." << endl; 
	cerr << "  no: otherwise." << endl << endl;
	
	cerr << "--path=yes[default:no]" << endl;
	cerr << "  yes: path information will be generated when produc labels." << endl; 
	cerr << "  no: otherwise." << endl << endl;	
	
	cerr << "--generate_label=yes[default:no]" << endl;
	cerr << "  yes: generate the hierarchical labels for each nodes." << endl; 
	cerr << "  no: otherwise." << endl << endl;		
	*/
}

void printParameters(const map<string, string>& cmdLineArgs) {
	map<string, string>::const_iterator iter;
	for (iter = cmdLineArgs.begin(); iter != cmdLineArgs.end(); iter++) {
		cout << iter->first << "=";
		if (iter->second=="")
			cout << "EMPTY";
		else 
			cout << iter->second;
		cout << " ";
	}
	cout << endl;
}

void keepResult(const char* resultFileName, vector<string>& results) {
	ofstream out(resultFileName, ios_base::out|ios_base::app);
	for (int i = 0; i < results.size(); i++)
		out << results[i] << "\t";
	out << endl;
	out.close();
}

void readFile(ifstream& in, vector<vector<double> >& data) {
	vector<vector<double> >(0).swap(data);
	string str;
	char* line;
	int count = 0;
	while (!in.eof()) {
		getline(in, str);
		if (str.size() == 0) break;
		line = new char[str.size()+1];
		strcpy(line, str.c_str());
		char* num = strtok(line, " \t#:");
		vector<double> tv;
		while(num != NULL) {
			tv.push_back(atof(num));
			num = strtok(NULL, " \t#:");
		}
		// vector<VertexID> tu;
		// for( int i = 0; i < tv.size(); i++ ) tu.push_back(tv[i]);
		data.push_back(tv);
		
		delete line;
	}
	in.close();
}

// check whether the queries have been generated or not;
// if yes, read the queries; otherwise, create new ones.
void checkQuery(Graph& graph, string file_name, int num_query, vector<VertexID>& query, bool save_para){
	string query_file = file_name;
	size_t found = query_file.find_last_of('.');
	query_file.erase(query_file.begin()+found, query_file.end());
	
	char buff[20];
	sprintf(buff, "%d", num_query);
	query_file += "_";
	query_file += buff;
	query_file += ".query";
	
	ifstream infile(query_file.c_str());
	if(infile && save_para){
		vector<vector<double> > data;
		readFile(infile, data);
		for ( int i = 0; i < data.size(); i++ ) {
			for (int j = 0; j < min(2, data[i].size()); j++) {
				query.push_back(static_cast<VertexID>(data[i][j]));
			}	
		}
	}else{
		cout << "Generate new query file: " + query_file << ".";
		GraphUtil gu;
		ofstream out(query_file.c_str());
		gu.generateQuery(graph, graph.num_vertices(), num_query, out, query);
		cout << " Done!" << endl;
	}
}

void readHierarchies(Graph& graph, string hierarchy_file){
	vector<ShortCuts> outshortcut, inshortcut;
	VertexList nodelist;
	int nodesize = graph.num_vertices();
	
	// outshortcutlist;
	vector<vector<double> > data;
	// cout << "file name:" << hierarchy_file << endl;
	size_t found = hierarchy_file.find_last_of(".");
	hierarchy_file.erase(hierarchy_file.begin()+found, hierarchy_file.end());
	hierarchy_file += ".hierarchy_out";
	ifstream in(hierarchy_file.c_str());
	readFile(in, data);
	
	// for ( int i = 0; i < data.size(); i++ ) {
		// for ( int j = 0;  j < data[i].size(); j++ ) {
			// cout << data[i][j] << " ";
		// }
		// cout << endl;
	// }
	
	// exit(0);
	
	outshortcut.resize(nodesize);
	for ( int i = 0; i < data.size(); i++ ){
		// cout << i << endl;
		int u = static_cast<int>(data[i][0]);
		int v = static_cast<int>(data[i][1]);
		if ( u > nodesize ) continue;
		// cout << i << " " << u << " " << v << endl;		
		ShortCutEdge tmp_sce;
		tmp_sce.target = v;
		tmp_sce.weight = data[i][2];
		for (int j = 3; j < data[i].size(); j++) tmp_sce.innerIDs.push_back(static_cast<int>(data[i][j]));

		outshortcut[u].push_back(tmp_sce);
	}
	
	// inshortcutlist;
	/*found = hierarchy_file.find_last_of("_");
	hierarchy_file.erase(hierarchy_file.begin()+found, hierarchy_file.end());
	hierarchy_file += "_in";
	in.open(hierarchy_file.c_str(), ifstream::in);
	readFile(in, data);*/
	
	
	// for ( int i = 0; i < data.size(); i++ ) {
		// for ( int j = 0;  j < data[i].size(); j++ ) {
			// cout << data[i][j] << " ";
		// }
		// cout << endl;
	// }	
	
	/*inshortcut.resize(nodesize);
	for ( int i = 0; i < data.size(); i++ ){
		int u = static_cast<int>(data[i][0]);
		int v = static_cast<int>(data[i][1]);
		if ( u > nodesize ) continue;
		//cout << i << " " << u << " " << v << endl;			
		ShortCutEdge tmp_sce;
		tmp_sce.target = v;
		tmp_sce.weight = data[i][2];
		for ( int j = 3; j < data[i].size(); j++ ) tmp_sce.innerIDs.push_back(static_cast<int>(data[i][j]));

		inshortcut[u].push_back(tmp_sce);
	}*/
	
	graph.insertShortcut(inshortcut, outshortcut);
	

	// nodelist;
	found = hierarchy_file.find_last_of("_");
	hierarchy_file.erase(hierarchy_file.begin()+found, hierarchy_file.end());
	hierarchy_file += "_node";
	cout << "file name:" << hierarchy_file << endl;	
	in.open(hierarchy_file.c_str(), ifstream::in);
	readFile(in, data);
	
	// for ( int i = 0; i < data.size(); i++ ) {
		// for ( int j = 0;  j < data[i].size(); j++ ) {
			// cout << data[i][j] << " ";
		// }
		// cout << endl;
	// }	
	
	nodelist.resize(nodesize);
	for ( int i = 0; i < data.size(); i++ ) {
		int u = static_cast<int>(data[i][0]);
		int v = static_cast<int>(data[i][1]);
		if (u > nodesize) continue;
		nodelist[u].rank = v;
		nodelist[u].id = u;
	}
	
	graph.insertNodeList(nodelist);
	// for ( int i = 0; i < graph.vertices().size(); i++) {
		// cout << i << " " << graph.vertices()[i].rank << endl;
	// }

	// read labels;
	// nodelist;
	// found = hierarchy_file.find_last_of(".");
	// hierarchy_file.erase(hierarchy_file.begin()+found, hierarchy_file.end());
	// hierarchy_file += ".label_in";
	// cout << "file name:" << hierarchy_file << endl;	
	// in.open(hierarchy_file.c_str(), ifstream::in);
	// readFile(in, data);	
	// graph.insertInLabel(data);
	
	found = hierarchy_file.find_last_of(".");
	hierarchy_file.erase(hierarchy_file.begin()+found, hierarchy_file.end());
	hierarchy_file += ".label_out";
	cout << "file name:" << hierarchy_file << endl;	
	in.open(hierarchy_file.c_str(), ifstream::in);
	readFile(in, data);	
	graph.insertOutLabel(data);	
	
	// read path index and list;
	// found = hierarchy_file.find_last_of(".");
	// hierarchy_file.erase(hierarchy_file.begin()+found, hierarchy_file.end());
	// hierarchy_file += ".path_index";
	// cout << "file name:" << hierarchy_file << endl;	
	// in.open(hierarchy_file.c_str(), ifstream::in);
	// readFile(in, data);	
	// graph.insertPathIndex(data);
	
	// found = hierarchy_file.find_last_of(".");
	// hierarchy_file.erase(hierarchy_file.begin()+found, hierarchy_file.end());
	// hierarchy_file += ".path";
	// cout << "file name:" << hierarchy_file << endl;	
	// in.open(hierarchy_file.c_str(), ifstream::in);
	// readFile(in, data);	
	// graph.insertPathList(data);		
}

void checkWritingCorrectness(Graph& graph, string hierarchy_file) {
	Graph tmp_g(graph.num_vertices());
	cout << "reading" << endl;
	readHierarchies(tmp_g, hierarchy_file);
	cout << "done" << endl;		
			
	// compare outshortcutlist;
	const ShortCuts& outshortcutOrig = graph.exportoutShortCutList();
	const ShortCuts& outshortcutDeri = tmp_g.exportoutShortCutList();
	if (outshortcutOrig.size() != outshortcutDeri.size()) { cout << "Size does not match!" << endl; exit(0); }
	for (int i = 0; i < outshortcutOrig.size(); i++) {
		// check id;
		if (outshortcutOrig[i].target != outshortcutDeri[i].target) {
			cout << "Node id does not match!" << endl;
			exit(0);
		}
		
		// check weight;
		if (fabs(outshortcutOrig[i].weight-outshortcutDeri[i].weight) > 1e-5) {
			cout << "Weight does not match!" << endl;
			exit(0);
		}
		
		// check inner ids;
		if (outshortcutOrig[i].innerIDs.size() != outshortcutDeri[i].innerIDs.size()) {
			cout << "InnerIDs.size does not match!" << endl;
			exit(0);
		}
		for (int j = 0; j < outshortcutOrig[i].innerIDs.size(); j++) {
			if (outshortcutOrig[i].innerIDs[j] != outshortcutDeri[i].innerIDs[j]){
				cout << "InnerIDs content does not match!" << endl;
				exit(0);
			}
		}
	}

	// compare shortcutindex;
	const vector<int>& outshortcutindexOrig = graph.exportoutShortCutIndex();
	const vector<int>& outshortcutindexDeri = tmp_g.exportoutShortCutIndex();
	if (outshortcutindexOrig.size() != outshortcutindexDeri.size()){
		cout << "index size does not match!" << endl;
		exit(0);
	}
	for (int i = 0; i < outshortcutindexOrig.size(); i++) {
		if (outshortcutindexOrig[i] != outshortcutindexDeri[i]) {
			cout << "index content does not match!" << endl;
			exit(0);
		}
	}
	
	
	// compare nodelist;
	int nodesize_o = graph.num_vertices();
	int nodesize_d = tmp_g.num_vertices();
	if (nodesize_o != nodesize_d) {
		cout << "node size does not match!" << endl;
		exit(0);
	}
	for (int i = 0; i < nodesize_o; i++) {
		if (graph[i].rank != tmp_g[i].rank){
			cout << "node rank does not match!" << endl;
			exit(0);
		}
		//cout << i << ": " << graph[i].id << "--" << tmp_g[i].id << endl;
		if (graph[i].id != tmp_g[i].id){
			cout << "node id does not match!" << endl;
			exit(0);
		}
	}		
}

// check whether hierachies have been built or not;
// if yes, read them in; otherwise, build them;
void checkHierarchies(Graph& graph, string file_name, bool save_para, bool gp){
	string hierarchy_file = file_name;
	size_t found = hierarchy_file.find_last_of('.');
	hierarchy_file.erase(hierarchy_file.begin()+found, hierarchy_file.end());
	hierarchy_file += ".hierarchy_out";

	ifstream infile(hierarchy_file.c_str());
	if (infile && save_para){
		cout << "read ... ";
		readHierarchies(graph, hierarchy_file);
		cout << "done ... ";
	}else{
		cout << "Create new hierarchy file: " << hierarchy_file << ".";
		Contraction c(graph);
		cout << "Creating hierarchy ... ";
		c.computeShortcuts(gp);		
		
		const vector<ShortCuts>& outshortcut = c.exportOutShortcut();
		const vector<ShortCuts>& inshortcut  = c.exportInShortcut();
		const VertexList& nodelist = c.exportNodeList();
		
		// write them down;
		// cout << "file name:" << hierarchy_file << endl;
		ofstream out(hierarchy_file.c_str());		
		for ( int i = 0; i < outshortcut.size(); i++ ) {
			for ( int j = 0; j < outshortcut[i].size(); j++ ) {
				out << i << "\t" << outshortcut[i][j].target << "\t" << outshortcut[i][j].weight << ": ";
				for ( int k = 0; k < outshortcut[i][j].innerIDs.size(); k++ ) {
					out << outshortcut[i][j].innerIDs[k] << " ";
				}
				out << endl;
			}
		}
 		out.close();
		
		// found = hierarchy_file.find_last_of("_");
		// hierarchy_file.erase(hierarchy_file.begin()+found, hierarchy_file.end());
		// hierarchy_file += "_in";
		// out.open(hierarchy_file.c_str(), ofstream::out);		
		// for ( int i = 0; i < inshortcut.size(); i++ ) {
			// for ( int j = 0; j < inshortcut[i].size(); j++ ) {
				// out << i << "\t" << inshortcut[i][j].target << "\t" << outshortcut[i][j].weight << ": ";
				// for ( int k = 0; k < inshortcut[i][j].innerIDs.size(); k++ ) {
					// out << inshortcut[i][j].innerIDs[k] << " ";
				// }
				// out << endl;
			// }
		// }
		// out.close();

		found = hierarchy_file.find_last_of("_");
		hierarchy_file.erase(hierarchy_file.begin()+found, hierarchy_file.end());
		hierarchy_file += "_node";
		// cout << "file name:" << hierarchy_file << endl;		
		out.open(hierarchy_file.c_str(), ofstream::out);	
		for ( int i = 0; i < nodelist.size(); i++ ) {
			out << i << ": " << nodelist[i].rank << endl;
		}
		out.close();
		// cout << "!!" << endl;
		
		// found = hierarchy_file.find_last_of(".");
		// hierarchy_file.erase(hierarchy_file.begin()+found, hierarchy_file.end());
		// hierarchy_file += ".label_in";
		// cout << "file name:" << hierarchy_file << endl;		
		// out.open(hierarchy_file.c_str(), ofstream::out);
		// for ( int i = 0; i < nodelist.size(); i++ ) {
			// LabelList& label = graph.exportInLabel(i);
			// out << i << ": ";
			// for ( int j = 0; j < label.size(); j++ ) {
				// out << label[j].id << " " << label[j].distance << " ";
			// }
			// out << endl;
		// }
		// out.close();

		found = hierarchy_file.find_last_of(".");
		hierarchy_file.erase(hierarchy_file.begin()+found, hierarchy_file.end());
		hierarchy_file += ".label_out";
		// cout << "file name:" << hierarchy_file << endl;		
		out.open(hierarchy_file.c_str(), ofstream::out);
		for ( int i = 0; i < nodelist.size(); i++ ) {
			LabelList& label = graph.exportOutLabel(i);
			out << i << ": ";
			for ( int j = 0; j < label.size(); j++ ) {
				out << label[j].id << " " << label[j].distance << " ";
			}
			out << endl;
		}
		out.close();	

		// found = hierarchy_file.find_last_of(".");
		// hierarchy_file.erase(hierarchy_file.begin()+found, hierarchy_file.end());
		// hierarchy_file += ".path_index";
		// cout << "file name:" << hierarchy_file << endl;		
		// out.open(hierarchy_file.c_str(), ofstream::out);
		// vector<vector<int> >& pathIndex = graph.exportPathIndex();
		// for ( int i = 0; i < pathIndex.size(); i++ ) {
			// out << i << ": ";
			// for ( int j = 0; j < pathIndex[i].size(); j++ ) {
				// out << pathIndex[i][j] << " ";
			// }
			// out << endl;
		// }
		// out.close();

		// found = hierarchy_file.find_last_of(".");
		// hierarchy_file.erase(hierarchy_file.begin()+found, hierarchy_file.end());
		// hierarchy_file += ".path";
		// cout << "file name:" << hierarchy_file << endl;		
		// out.open(hierarchy_file.c_str(), ofstream::out);
		// vector<vector<VertexID> >& pathList = graph.exportPathList(); 
		// for ( int i = 0; i < pathList.size(); i++ ) {
			// out << i << ": ";
			// for ( int j = 0; j < pathList[i].size(); j++ ) {
				// out << pathList[i][j] << " ";
			// }
			// out << endl;
		// }
		// out.close();	

		// for testing
		// checkWritingCorrectness(graph, hierarchy_file);
	}
}

void readIDMap(string filename, vector<VertexID>& id_map){
	size_t found = filename.find_last_of(".");
	filename.erase(filename.begin()+found, filename.end());
	filename += ".idmap";
	
	ifstream in(filename.c_str());
	vector<vector<double> > data;
	readFile(in, data);
	
	vector<VertexID>(0).swap(id_map);
	for (int i = 0; i < data.size(); i++) {
		id_map.push_back(static_cast<VertexID>(data[i][1]));
	}
}

bool isNumeric(string input){
	if (input.size() == 0) return false;
	for (int i = 0; i < input.size(); i++) {
		if (input[i] < '0' || input[i] > '9') return false;
	}
	return true;
}
