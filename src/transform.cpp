#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <sys/time.h>
#include <cmath>
#include <stdlib.h>
#include <stdio.h>
#include <set>
#include <map>
#include <algorithm>
#include <cstring>
#include <cassert>

using namespace std;

template<class T>
struct index_comp {
	vector<T>& _x;
	
	index_comp(vector<T>& x): _x(x) {};
	bool operator() (const int a, const int b) const {
		return _x[a]>_x[b];
	}
};

static void strTrimRight(string& str) {
	string whitespaces(" \t\n");
	int index = str.find_last_not_of(whitespaces);
	if (index != string::npos) 
		str.erase(index+1);
	else
		str.clear();
}

static void strTrimLeft(string& str) {
	string whitespaces(" \t");
	int index = str.find_first_not_of(whitespaces);
	if (index != string::npos) 
		str.erase(0,index);
	else
		str.clear();
}

static int str2int(string str) {
	int value = -1;
	istringstream myss(str);
	myss >> value;
	return value;
}

static double str2double(string str) {
	double value = 0.0;
	istringstream myss(str);
	myss >> value;
	return value;
}

vector<string>& split(const string &s, char delim, vector<string> &elems) {
	stringstream ss(s);
	string item;
	while(getline(ss, item, delim)) {
		elems.push_back(item);
	}
	return elems;
}

vector<string> split(const string &s, char delim) {
	vector<string> elems;
	return split(s, delim, elems);
}

static void remapgraph(char* mapfile, char* filename) {
	ifstream dagin(mapfile);
	map<int, int> infomap;
	int code;
	int id;
	string buf;
	int idx;
	while (getline(dagin, buf)) {
		strTrimRight(buf);
		if (buf.length() == 0) break;
		idx = buf.find_first_of("\t");
		id = str2int(buf.substr(0,idx));
		buf.erase(0,idx+1);
		code = str2int(buf);
		cout << id << "\t" << code << endl;
		infomap[id] = code;
	}
	dagin.close();
	
	string filestr(filename);
	filestr = filestr.substr(0,filestr.find_last_of("."));
	filestr = filestr + "_o.bb";
	ofstream out(filestr.c_str());	
	out << "graph_for_backbone" << endl;

	ifstream in(filename);
	getline(in, buf);
	strTrimRight(buf);
	if (buf != "graph_for_backbone") {
		cout << "BAD FILE FORMAT!" << endl;
		exit(0);
	}
	
	int n;
	getline(in, buf);
	istringstream(buf) >> n;	
	out << n << endl;
	string sub;

	int sid = 0;
	int tid = 0;
	while (getline(in, buf)) {
		strTrimRight(buf);
		idx = buf.find(":");
		sid = str2int(buf.substr(0,idx));
		buf.erase(0, idx+2);
		out << infomap[sid] << ": ";
		while (buf.find(" ") != string::npos) {
			sub = buf.substr(0, buf.find(" "));
			istringstream(sub) >> tid;
			buf.erase(0, buf.find(" ")+1);
			out << infomap[tid] << " ";
		}
		out << "#" << endl;
	}
	in.close();
	out.close();
}

static void siftograph(char* siffile) {
	ifstream sifin(siffile);
	int pos, index=0;
	string buf, str;
	map<string,int> labelmap;
	while (getline(sifin,buf)) {
		strTrimRight(buf);
		if (buf.length() == 0) break;
		pos = buf.find_first_of("\t");
		str = buf.substr(0,pos);
		cout << "|" << str << "|\t";
		if (labelmap.find(str)==labelmap.end()) {
			labelmap[str] = index;
			index++;
		}
		buf.erase(0,pos+1);
		buf.erase(0,buf.find_first_of("\t")+1);
		str = buf;
		cout << "|" << str << "|" << endl;
		if (labelmap.find(str)==labelmap.end()) {
			labelmap[str] = index;
			index++;
		}
	}
	sifin.close();
	// write map
	string filestr1(siffile);
	filestr1 = filestr1.substr(0,filestr1.find_last_of("."));
	filestr1 = filestr1 + ".map";
	ofstream out1(filestr1.c_str());
	map<string,int>::iterator mit;
	for (mit = labelmap.begin(); mit != labelmap.end(); mit++)
		out1 << mit->second << "\t" << mit->first << endl;
	out1.close();
	
	sifin.open(siffile);
	string src, trg;
	int sid, tid;
	vector<set<int> > gvec = vector<set<int> >(labelmap.size(),set<int>());
	while (getline(sifin,buf)) {
		strTrimRight(buf);
		if (buf.length() == 0) break;
		pos = buf.find_first_of("\t");
		src = buf.substr(0,pos);
		buf.erase(0,pos+1);
		buf.erase(0,buf.find_first_of("\t")+1);
		trg = buf;
		sid = labelmap[src];
		tid = labelmap[trg];
		gvec[sid].insert(tid);
		gvec[tid].insert(sid);
	}
	sifin.close();
	
	set<int>::iterator sit;
	string filestr(siffile);
	filestr = filestr.substr(0,filestr.find_last_of("."));
	filestr = filestr + ".gra";
	ofstream out(filestr.c_str());
	out << "graph_for_backbone" << endl;
	out << gvec.size() << endl;
	for (int i = 0; i < gvec.size(); i++) {
		out << i << ": ";
		for (sit = gvec[i].begin(); sit != gvec[i].end(); sit++) 
			out << *sit << " ";
		out << "#" << endl;
	}
	out.close();
}

static void mapcountrycode(char* mapfile, char* bbfile) {
	ifstream bbin(bbfile);
	int n;
	string buf;
	getline(bbin, buf);
	istringstream(buf) >> n;
	int idx;
	int sid = 0;
	int tid = 0;
	vector<int> cids;
	while (getline(bbin, buf)) {
		strTrimRight(buf);
		if (buf.length() == 0) break;
		sid = str2int(buf.substr(0,buf.find(":")));
		cids.push_back(sid);
	}	
	cids.erase(cids.begin());
	bbin.close();
	
	
	ifstream dagin(mapfile);
	map<int, string> infomap;
	string code;
	int id;
	/*
	for (int i = 0; i < 196; i++) {
		dagin >> id  >> code;
		infomap[id] = code;
	}
	*/	
	while (getline(dagin, buf)) {
		strTrimRight(buf);
		if (buf.length() == 0) break;
		idx = buf.find_first_of("\t");
		id = str2int(buf.substr(0,idx));
		buf.erase(0,idx+1);
		code = buf;
		cout << id << "\t" << code << endl;
		infomap[id] = code;
	}
	
	dagin.close();
	
		// output to file
	string filestr(bbfile);
	filestr = filestr.substr(0,filestr.find_last_of("."));
	filestr = filestr + ".result";
	ofstream out(filestr.c_str());	
	for (int i = 0; i < cids.size(); i++) {
		out << cids[i] << "\t" << infomap[cids[i]] << endl;
	}
	out.close();
}

static void tranformtradedata(char* filename) {
	ifstream in(filename);
	map<string,int> codemap;
	string buf, code;
	getline(in,buf);
	int index = 0, pos;
	while (getline(in,buf)) {
		strTrimRight(buf);
		if (buf.length() == 0) break;
		pos = buf.find_first_of(" ");
		code = buf.substr(0,pos);
		if (codemap.find(code) == codemap.end()) {
			codemap[code] = index;
			index++;
		}
		buf.erase(0,pos+1);
		buf.erase(0,buf.find_first_of(" ")+1);
		code = buf.substr(0,buf.find_first_of(" "));
		if (codemap.find(code) == codemap.end()) {
			codemap[code] = index;
			index++;
		}
	}
	in.close();
	// write country codemap
	string filestr(filename);
	filestr = filestr.substr(0,filestr.find_last_of("."));
	filestr = filestr + ".map";
	ofstream out(filestr.c_str());
	map<string,int>::iterator mit;
	for (mit = codemap.begin(); mit != codemap.end(); mit++)
		out << mit->second << "\t" << mit->first << endl;
	out.close();
	
	in.open(filename);
	int periods[] = {};
	vector<vector<vector<double> > > data = vector<vector<vector<double> > >(4,vector<vector<double> >(196,vector<double>(196,0)));
	getline(in,buf);
	string codea, codeb;
	int year, inda, indb;
	double expab, impab, expba, impba;
	while (getline(in,buf)) {
		strTrimRight(buf);
		if (buf.length() == 0) break;
		pos = buf.find_first_of(" ");
		codea = buf.substr(0,pos);
		buf.erase(0,pos+1);
		buf.erase(0,buf.find_first_of(" ")+1);
		pos = buf.find_first_of(" ");
		codeb = buf.substr(0,pos);
		buf.erase(0,pos+1);
		buf.erase(0,buf.find_first_of(" ")+1);
		pos = buf.find_first_of(" ");
		year = str2int(buf.substr(0,pos));
		buf.erase(0,pos+1);
		pos = buf.find_first_of(" ");
		expab = str2double(buf.substr(0,pos));
		buf.erase(0,pos+1);
		buf.erase(0,buf.find_first_of(" ")+1);
		pos = buf.find_first_of(" ");
		impab = str2double(buf.substr(0,pos));
		buf.erase(0,pos+1);
		buf.erase(0,buf.find_first_of(" ")+1);
		pos = buf.find_first_of(" ");
		expba = str2double(buf.substr(0,pos));
		buf.erase(0,pos+1);
		buf.erase(0,buf.find_first_of(" ")+1);
		pos = buf.find_first_of(" ");
		impba = str2double(buf.substr(0,pos));
		buf.erase(0,pos+1);
//		cout << codea << "\t" << codeb << "\t" << year << "\t" << expab << "\t" << impab << "\t" << expba << "\t" << impba << endl;
		index = year-1950;
		if (index<0) continue;
		inda = codemap[codea];
		indb = codemap[codeb];
		if (index<=20) { 
			data[0][inda][indb] += (expab+impab+expba+impba)*0.5; 
			data[0][indb][inda] += (expab+impab+expba+impba)*0.5;
		}
		if (index>10 && index<=30) { 
			data[1][inda][indb] += (expab+impab+expba+impba)*0.5; 
			data[1][indb][inda] += (expab+impab+expba+impba)*0.5;
		}
		if (index>20 && index<=40) { 
			data[2][inda][indb] += (expab+impab+expba+impba)*0.5; 
			data[2][indb][inda] += (expab+impab+expba+impba)*0.5;
		}
		if (index>30 && index<=50) { 
			data[3][inda][indb] += (expab+impab+expba+impba)*0.5; 
			data[3][indb][inda] += (expab+impab+expba+impba)*0.5;
		}
	}
	in.close();
	
	vector<vector<set<int> > > gmaps = vector<vector<set<int> > >(4, vector<set<int> >(196,set<int>()));
	string filenames[] = {"trade50_70.gra", "trade60_80.gra", "trade70_90.gra", "trade80_00.gra"};
	vector<int> temp, order;
	for (int i = 0; i < 196; i++)
		temp.push_back(i);
	for (int k = 0; k < 4; k++) {
		for (int i = 0; i < 196; i++) {
			order = temp;
			sort(order.begin(),order.end(),index_comp<double>(data[k][i]));
			/*
			for (int m = 0; m < data[k][i].size(); m++)
				cout << data[k][i][m] << " ";
			cout << endl;
			for (int m = 0; m < order.size(); m++)
				cout << order[m] << " ";
			cout << endl;
			if (true) exit(0);
			*/
			for (int j = 0; j < 15; j++) {
				if (data[k][i][order[j]]>0) {
					gmaps[k][i].insert(order[j]);
			//		gmaps[k][order[j]].insert(i);
				}
			}
		}
	}
	
	set<int>::iterator sit;
	// shrink the graph
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < gmaps[i].size(); j++) {
			for (sit = gmaps[i][j].begin(); sit != gmaps[i][j].end(); ) {
				cout << "check " << j << "->" << *sit << endl;
				if (gmaps[i][*sit].find(j) == gmaps[i][*sit].end())
					gmaps[i][j].erase(sit++);
				else
					sit++;
			}
		}	
	}
	
	for (int i = 0; i < 4; i++) {
		ofstream infile(filenames[i].c_str());
		infile << "graph_for_backbone" << endl;
		infile << 196 << endl;
		for (int j = 0; j < gmaps[i].size(); j++) {
			infile << j << ": ";
			for (sit = gmaps[i][j].begin(); sit != gmaps[i][j].end(); sit++)
				infile << *sit << " ";
			infile << "#" << endl;
		}
		infile.close();
	}
}

static void readpowerlawedges(char* filename) {
	ifstream in(filename);
	int n;
	in >> n;
	cout << "n = " << n << endl;
	string buf;
	int sid, tid;
	set<int>::iterator sit;
	vector<set<int> > gvec = vector<set<int> >(n,set<int>());
	getline(in,buf);
	while (getline(in,buf)) {
		strTrimRight(buf);
		if (buf.length() == 0) break;
		sid = str2int(buf.substr(0,buf.find(" ")));
		tid = str2int(buf.substr(buf.find_last_of(" ")+1));
		gvec[sid].insert(tid);
		gvec[tid].insert(sid);
	}
	in.close();
		// output to file
	string filestr(filename);
	filestr = filestr.substr(0,filestr.find_last_of("."));
	filestr = filestr + ".gra";
	ofstream out(filestr.c_str());
	out << "graph_for_backbone" << endl;
	out << n << endl;
	for (int i = 0; i < gvec.size(); i++) {
		out << i << ": ";
		for (sit = gvec[i].begin(); sit != gvec[i].end(); sit++) 
			out << *sit << " ";
		out << "#" << endl;
	}
	out.close();
}

static void mapauthors(char* dagfile, char* bbfile) {
	ifstream bbin(bbfile);
	int n;
	string buf;
	getline(bbin, buf);
	istringstream(buf) >> n;
	int idx;
	int sid = 0;
	int tid = 0;
	vector<int> authorids;
	while (getline(bbin, buf)) {
		strTrimRight(buf);
		if (buf.length() == 0) break;
		sid = str2int(buf.substr(0,buf.find(":")));
		authorids.push_back(sid);
	}	
	authorids.erase(authorids.begin());
	bbin.close();
	
	ifstream dagin(dagfile);
	getline(dagin, buf);
	strTrimRight(buf);
	buf = buf.substr(buf.find(" ")+1);
	int gsize = str2int(buf);
	map<int, pair<string,int> > infomap;
	for (int i = 0; i < gsize; i++) {
		getline(dagin, buf);
		strTrimRight(buf);
		if (buf.length() == 0) break;
		if (buf[0] == '#') continue;
		idx = buf.find_first_of("\t");
		sid = str2int(buf.substr(0,idx));
		buf.erase(0,idx);
		strTrimLeft(buf);
		idx = buf.find_first_of("\t");
		string name = buf.substr(0,idx);
		buf.erase(0,idx);
		strTrimLeft(buf);
		int pnum = str2int(buf);
		infomap[sid] = pair<string,int>(name,pnum);
	}
	dagin.close();
	
		// output to file
	string filestr(bbfile);
	filestr = filestr.substr(0,filestr.find_last_of("."));
	filestr = filestr + ".result";
	ofstream out(filestr.c_str());	
	for (int i = 0; i < authorids.size(); i++) {
		out << authorids[i] << "\t" << infomap[authorids[i]].first << "\t" << infomap[authorids[i]].second << endl;
	}
	out.close();
}

static void readEdgeList(char *filename) {
	ifstream in(filename);	
	string buf;
	getline(in, buf);
	strTrimRight(buf);
	buf = buf.substr(buf.find(" ")+1);
	int gsize = str2int(buf);
	vector<set<int> > gvec = vector<set<int> >(gsize, set<int>());
	set<int>::iterator sit;
	
	while (getline(in, buf)) 
		if (buf[0] == '*') break;
	
	int src, trg;
	int index, index0;
	while (getline(in, buf)) {
		strTrimRight(buf);
		if (buf.length() == 0) break;
		if (buf[0] == '#') continue;
		index = buf.find_first_of("\t");
		src = str2int(buf.substr(0,index));
		index0 = buf.find_last_of("\t");
		trg = str2int(buf.substr(index+1,index0));
		// undirected graph
		gvec[src].insert(trg);
		gvec[trg].insert(src);
	}
	in.close();	

		// output to file
	string filestr(filename);
	filestr = filestr.substr(0,filestr.find_last_of("."));
	filestr = filestr + ".gra";
	ofstream out(filestr.c_str());
	out << "graph_for_backbone" << endl;
	out << gsize << endl;
	for (int i = 0; i < gvec.size(); i++) {
		out << i << ": ";
		for (sit = gvec[i].begin(); sit != gvec[i].end(); sit++) 
			out << *sit << " ";
		out << "#" << endl;
	}
	out.close();
}

static void edgelistToGraph(char* filename, int gsize) {
	ifstream in(filename);
	vector<set<int> > gvec = vector<set<int> >(gsize, set<int>());
	set<int>::iterator sit;
	
	int src, trg;
	int index = -1;
	string buf;
	while (getline(in, buf)) {
		strTrimRight(buf);
		if (buf.length() == 0) break;
		if (buf[0] == '#') continue;
		index = buf.find_first_of("\t");
		src = str2int(buf.substr(0,index));
		trg = str2int(buf.substr(index+1));
		gvec[src].insert(trg);
		gvec[trg].insert(src);
	}
	in.close();
	
	// output to file
	string filestr(filename);
	filestr = filestr.substr(0,filestr.find_last_of("."));
	filestr = filestr + ".gra";
	ofstream out(filestr.c_str());
	out << "graph_for_backbone" << endl;
	out << gsize << endl;
	for (int i = 0; i < gvec.size(); i++) {
		out << i << ": ";
		for (sit = gvec[i].begin(); sit != gvec[i].end(); sit++) 
			out << *sit << " ";
		out << "#" << endl;
	}
	out.close();
}

static void graphToEdgelist(char* filename) {
	map<int,set<int> > gset;
	ifstream in(filename);
	string buf;
	getline(in, buf);
	strTrimRight(buf);
	if (buf != "graph_for_greach") {
		cout << "BAD FILE FORMAT!" << endl;
	//	exit(0);
	}
	
	int n;
	getline(in, buf);
	istringstream(buf) >> n;

	string filestr(filename);
	filestr = filestr.substr(0,filestr.find_last_of("."));
	filestr = filestr + ".txt";
	ofstream out(filestr.c_str());	
	
	string sub;
	int idx;
	int sid = 0;
	int tid = 0;
	while (getline(in, buf)) {
		strTrimRight(buf);
		idx = buf.find(":");
		buf.erase(0, idx+2);
		while (buf.find(" ") != string::npos) {
			sub = buf.substr(0, buf.find(" "));
			istringstream(sub) >> tid;
			buf.erase(0, buf.find(" ")+1);
			if (gset[sid].find(tid)==gset[sid].end() && gset[tid].find(sid)==gset[tid].end()) {
				gset[sid].insert(tid);
				out << sid << "\t" << tid << endl;
			}
		}
		++sid;
	}
	out.close();
	in.close();
}

static void PowerlawEdgelistToGraph(char* filename) {
	ifstream in(filename);
	string buf;
	getline(in, buf);
	int gsize;
	istringstream(buf) >> gsize;

	vector<set<int> > gvec = vector<set<int> >(gsize, set<int>());
	set<int>::iterator sit;
	
	int src, trg;
	int index = -1;
	while (getline(in, buf)) {
		index = buf.find_first_of(",");
		src = str2int(buf.substr(0,index));
		trg = str2int(buf.substr(index+1));
		gvec[src].insert(trg);
		//gvec[trg].insert(src);
	}
	in.close();
	
	// output to file
	string filestr(filename);
	filestr = filestr.substr(0,filestr.find_last_of("."));
	filestr = filestr + ".gra";
	ofstream out(filestr.c_str());
	out << "graph_for_greach" << endl;
	out << gsize << endl;
	for (int i = 0; i < gvec.size(); i++) {
		out << i << ": ";
		for (sit = gvec[i].begin(); sit != gvec[i].end(); sit++) 
			out << *sit << " ";
		out << "#" << endl;
	}
	out.close();
}

static void PowerlawEdgelistToGraphForGTGraph(char* filename) {
	ifstream in(filename);
	string buf;
	getline(in, buf);
	int gsize;
	istringstream(buf) >> gsize;

	vector<set<int> > gvec = vector<set<int> >(gsize, set<int>());
	set<int>::iterator sit;
	
	int src, trg;
	int index = -1;
	while (getline(in, buf)) {
		index = buf.find_first_of(" ");
		src = str2int(buf.substr(0,index));
		trg = str2int(buf.substr(index+1));
		gvec[src-1].insert(trg-1);
		//gvec[trg].insert(src);
	}
	in.close();
	
	// output to file
	string filestr(filename);
	filestr = filestr.substr(0,filestr.find_last_of("."));
	filestr = filestr + ".gra";
	ofstream out(filestr.c_str());
	out << "graph_for_greach" << endl;
	out << gsize << endl;
	for (int i = 0; i < gvec.size(); i++) {
		out << i << ": ";
		for (sit = gvec[i].begin(); sit != gvec[i].end(); sit++) 
			out << *sit << " ";
		out << "#" << endl;
	}
	out.close();
}

static void PowerlawEdgelistToGraphForBarabasiGraph(char* filename) {
	ifstream in(filename);
	string buf;
	getline(in, buf);
	int gsize;
	istringstream(buf) >> gsize;

	vector<set<int> > gvec = vector<set<int> >(gsize, set<int>());
	set<int>::iterator sit;
	
	int src, trg;
	int index = -1;
	while (getline(in, buf)) {
		index = buf.find_first_of(" ");
		src = str2int(buf.substr(0,index));
		trg = str2int(buf.substr(index+3));
		gvec[src].insert(trg);
		gvec[trg].insert(src);
	}
	in.close();
	
	// output to file
	int edgesize = 0;
	string filestr(filename);
	filestr = filestr.substr(0,filestr.find_last_of("."));
	filestr = filestr + ".gra";
	ofstream out(filestr.c_str());
	out << "graph_for_distance" << endl;
	out << gsize << endl;
	for (int i = 0; i < gvec.size(); i++) {
		out << i << ": ";
		edgesize += gvec[i].size();
		for (sit = gvec[i].begin(); sit != gvec[i].end(); sit++) 
			out << *sit << " ";
		out << "#" << endl;
	}
	out.close();
	cout << "#E=" << edgesize << endl;
}

static void elisttogra(char* filename) {
	ifstream in(filename);
	string buf;
	getline(in, buf);

	int src, trg;
	int index = -1;
	int id = 0;
	map<int,int> idmap;
	map<int,int>::iterator iter;
	map<int,set<int> > graph;
	
	while (getline(in, buf)) {
		index = buf.find_first_of("\t");
		src = str2int(buf.substr(0,index));
		buf = buf.substr(index+1);
		trg = str2int(buf);
		if (idmap.find(src) == idmap.end()) {
			idmap[src] = id;
			id++;
		}
		if (idmap.find(trg) == idmap.end()) {
			idmap[trg] = id;
			id++;
		}
		graph[src].insert(trg);
	}
	
	for (iter = idmap.begin(); iter != idmap.end(); iter++) {
		cout << iter->first << ": " << iter->second << endl;
	}
	
	string filestr(filename);
	filestr = filestr.substr(0,filestr.find_last_of("."));
	string filestr1 = filestr + ".el";
	ofstream out(filestr1.c_str());
	map<int,set<int> >::iterator mit;
	set<int>::iterator sit;
	vector<set<int> > gg = vector<set<int> >(id, set<int>());
	int num=0;
	for (mit = graph.begin(); mit != graph.end(); mit++) {
		for (sit = mit->second.begin(); sit != mit->second.end(); sit++)
			if (mit->first != *sit) {
				out << idmap[mit->first] << "\t" << idmap[*sit] << endl;
				cout << num << "\t" << mit->first<<"->" << *sit << "\t" << idmap[mit->first] << "->" << idmap[*sit] << endl;
				gg[idmap[mit->first]].insert(idmap[*sit]);
				gg[idmap[*sit]].insert(idmap[mit->first]);
				num++;
			}
			else {
				cout << mit->first << "===" << *sit << endl;
			}
	}
	out.close();
	
	cout << "here" << endl;
	string graphfile = filestr+".gra";
	ofstream out1(graphfile.c_str());
	out1 << "graph_for_simple" << endl;
	out1 << id << endl;
	vector<set<int> >::iterator vsit;
	index = 0;
	for (vsit = gg.begin(); vsit != gg.end(); vsit++, index++) {
		out1 << index << ": ";
		for (sit = gg[index].begin(); sit != gg[index].end(); sit++)
			out1 << *sit << " ";
		out1 << "#"  << endl;
	}
	out1.close();
}

static void translate(char* filename) {
	ifstream in(filename);
	string buf;
	getline(in, buf);

	int src, trg;
	int index = -1;
	int id = 0;
	map<int,int> idmap;
	map<int,set<int> > graph;
	
	while (getline(in, buf)) {
		index = buf.find_first_of(",");
		src = str2int(buf.substr(1,index));
		buf = buf.substr(index+1);
		index = buf.find_first_of(",");
		trg = str2int(buf.substr(1,index));
		if (idmap.find(src) == idmap.end()) {
			idmap[src] = id;
			id++;
		}
		if (idmap.find(trg) == idmap.end()) {
			idmap[trg] = id;
			id++;
		}
		graph[src].insert(trg);
	}
	
	string filestr(filename);
	filestr = filestr.substr(0,filestr.find_last_of("."));
	string filestr1 = filestr + ".el";
	ofstream out(filestr1.c_str());
	map<int,set<int> >::iterator mit;
	set<int>::iterator sit;
	for (mit = graph.begin(); mit != graph.end(); mit++) {
		for (sit = mit->second.begin(); sit != mit->second.end(); sit++)
			if (mit->first != *sit)
			out << idmap[mit->first] << "\t" << idmap[*sit] << endl;
	}
	out.close();
	
	string graphfile = filestr+".gra";
	
}

static void degdist(char* filename) {
	map<int,set<int> > gset;
	map<int,set<int> >::iterator msit;
	map<int,int> degs;
	ifstream in(filename);
	string buf;
	getline(in, buf);
	strTrimRight(buf);
	if (buf != "graph_for_greach") {
		cout << "BAD FILE FORMAT!" << endl;
	//	exit(0);
	}
	
	int n;
	getline(in, buf);
	istringstream(buf) >> n;

	string filestr(filename);
	filestr = filestr.substr(0,filestr.find_last_of("."));
	filestr = filestr + ".degs";
	ofstream out(filestr.c_str());	
	
	string sub;
	int idx;
	int sid = 0;
	int tid = 0;
	while (getline(in, buf)) {
		strTrimRight(buf);
		idx = buf.find(":");
		buf.erase(0, idx+2);
		int num = 0;
		while (buf.find(" ") != string::npos) {
			sub = buf.substr(0, buf.find(" "));
			istringstream(sub) >> tid;
			buf.erase(0, buf.find(" ")+1);
			gset[sid].insert(tid);
			num++;
		}
		degs[num]++;
		++sid;
	}
	in.close();
	
	/*
	for (msit = gset.begin(); msit != gset.end(); msit++) {
		degs[msit->second.size()]++;
	}
	*/
	
	// output degree distribution
	map<int,int>::iterator mit;
	for (mit = degs.begin(); mit != degs.end(); mit++) {
		out << mit->first << "\t" << mit->second << endl;
	}
	out.close();
}

static void subgraphToEdgelist(char* filename) {
	map<int,set<int> > gset;
	ifstream in(filename);
	string buf;
	getline(in, buf);
	strTrimRight(buf);
	if (buf != "graph_for_greach") {
		cout << "BAD FILE FORMAT!" << endl;
	//	exit(0);
	}
	
	int n;
	getline(in, buf);
	istringstream(buf) >> n;

	string filestr(filename);
	filestr = filestr.substr(0,filestr.find_last_of("."));
	string ofilestr = filestr + ".txt";
	ofstream out(ofilestr.c_str());	
	
	string sub;
	int idx;
	int sid = 0;
	int tid = 0;
	while (getline(in, buf)) {
		strTrimRight(buf);
		idx = buf.find(":");
		buf.erase(0, idx+2);
		while (buf.find(" ") != string::npos) {
			sub = buf.substr(0, buf.find(" "));
			istringstream(sub) >> tid;
			buf.erase(0, buf.find(" ")+1);
			if (tid>=n) continue;
			gset[sid].insert(tid);
			gset[tid].insert(sid);
		//	cout << sid << " -> " << tid << endl;
		//	if (gset[sid].find(tid)==gset[sid].end() && gset[tid].find(sid)==gset[tid].end()) {
		//		gset[sid].insert(tid);
		//	}
		}
		++sid;
	}

	in.close();
	cout << "here" << endl;
	
	// output edgelist
	string filestr1 = filestr + ".el";
	ofstream out1(filestr1.c_str());
	map<int,set<int> >::iterator mit;
	set<int>::iterator sit;
	out << "graph_for_simple" << endl;
	out << n << endl;
	for (mit = gset.begin(); mit != gset.end(); mit++) {
		out << mit->first << ": ";
		for (sit = mit->second.begin(); sit != mit->second.end(); sit++)
			if (mit->first != *sit) {
				out << *sit << " ";
				if (mit->first < *sit)
					out1<< mit->first << "\t" << *sit << endl;
			}
		out << "#" << endl;
	}
	out1.close();	
	out.close();
}

static void topicmatrixToGraph(char* filename, int gs) {
	vector<set<int> > tg = vector<set<int> >(gs, set<int>());
	ifstream in(filename);
	double input;
	int row = 0, col = 0;
	for (row = 0; row < gs; row++) {
		for (col = 0; col < gs; col++) {
			in >> input;
			if (row == col) continue;
			if (input > 0)
				tg[row].insert(col);
		}
	}
	in.close();
	
	string filestr(filename);
	filestr = filestr.substr(0,filestr.find_last_of("."));
	string filestr1 = filestr + ".gra";
	ofstream out(filestr1.c_str());
	vector<set<int> >::iterator vsit;
	set<int>::iterator sit;
	out << "graph_for_visulization" << endl;
	out << gs << endl;
	row = 0;
	int numedges = 0;
	for (vsit = tg.begin(); vsit != tg.end(); vsit++, row++) {
		out << row << ": ";
		for (sit = vsit->begin(); sit != vsit->end(); sit++) {
			out << *sit << " ";
		}
		numedges += vsit->size();
		out << "#" << endl;
	}
	out.close();
	cout << "#edges=" << numedges/2 << " #E/#V=" << numedges/(2*gs) << endl;
}

static void pywebgraph2gra(char* filename, int gnum, double density) {
	int edgenum = (int)(gnum*density);
	cout << "gnum=" << gnum << " edgenum=" << edgenum << endl;
	ifstream in(filename);
	string buf;
	
	int src, trg, s, t;
	int index = -1;
	int id = 0;
	map<int,int> idmap;
	map<int,int>::iterator iter;
	map<int,set<int> > graph;
	map<int,set<int> >::iterator msit;
	set<int>::iterator sit;
	int maxid = -1;
	set<long> pset;
	long label;
	int ecounter = 0;
	while (getline(in, buf)) {
		index = buf.find_first_of(" ");
		src = str2int(buf.substr(0,index));
		index = buf.find_last_of(" ");
		trg = str2int(buf.substr(index+1));
		if (src<trg) label = (src-1)*gnum+(trg-1);
		else label = (trg-1)*gnum+(src-1);
		if (pset.find(label)!=pset.end()) continue;
		if (src-1>=gnum||trg-1>=gnum) continue;
		graph[src-1].insert(trg-1);
		pset.insert(label);
		if (pset.size()>=edgenum) break;
//		graph[trg-1].insert(src-1);
		if (src-1>maxid) maxid=src-1;
		if (trg-1>maxid) maxid=trg-1;
	}
	cout << "maxid=" << maxid << endl;
	pset.clear();
	multimap<int,int> degmap;
	for (msit = graph.begin(); msit != graph.end(); msit++) {
		degmap.insert(make_pair(msit->second.size(),msit->first));
	}
	cout << "degmap size=" << degmap.size() << endl;
	multimap<int,int>::reverse_iterator mrit;
	int count = 0;
	vector<int> order = vector<int>(gnum,-1);
	for (mrit = degmap.rbegin(); mrit != degmap.rend(); mrit++,count++) {
		if (mrit->second>=gnum) continue;
		order[mrit->second] = count;
	//	cout << mrit->second << " -> " << count << endl;
	}
	vector<set<int> > g = vector<set<int> >(gnum,set<int>());
	for (int i = 0; i < gnum; i++) {
		if (order[i]==-1) {
			order[i] = count++;
		}
	}
	int ne = 0;
	for (mrit = degmap.rbegin(); mrit != degmap.rend(); mrit++,count++) {
		int id = mrit->second;
		for (sit = graph[id].begin(); sit != graph[id].end(); sit++) {
			if (order[id]>order[*sit])
				g[*sit].insert(id);
			else
				g[id].insert(*sit);
			ne++;
		}
	}
	degmap.clear();
	cout << "NE=" << ne << endl;
	/*
	for (int i = 0; i < gnum; i++) {
		cout << i << ": ";
		for (sit = g[i].begin(); sit != g[i].end(); sit++)
			cout << *sit << " ";
		cout << endl;
	}
	*/
	int noaddtime=0;
	srand(time(NULL));
	while (ne<=edgenum && noaddtime<10000000) {
		s = lrand48()%gnum;
		t = lrand48()%gnum;
		if (s==t) continue;
		/*
		cout << s << "->" << t << endl;
		if (g[s].find(t)!=g[s].end()) {
			cout << "find " << t << endl;
			char ch;
			cin >> ch;
		}
		if (g[t].find(s)!=g[t].end()) {
			cout << "find " << s << endl;
			char ch;
			cin >> ch;
		}
		*/
		noaddtime++;
		if (g[s].find(t)!=g[s].end()||g[t].find(s)!=g[t].end()) continue;
		noaddtime = 0;
//		cout << "add " << s << " -> " << t << endl;
		if (order[s]<order[t])
			g[s].insert(t);
		else
			g[t].insert(s);
		ne++;
		/*
		for (int i = 0; i < gnum; i++) {
			cout << i << ": ";
			for (sit = g[i].begin(); sit != g[i].end(); sit++)
				cout << *sit << " ";
			cout << endl;
		}
		*/
	}
	cout << "#E=" << ne << endl;
	cout << "output graph file" << endl;
	string filestr(filename);
	filestr = filestr.substr(0,filestr.find_last_of("."));
	filestr += ".gra";
	ofstream out(filestr.c_str());
	out << "graph_for_reach_scalefree" << endl;
	out << gnum << endl;
	for (int i = 0; i < gnum; i++) {
		out << i << ": ";
		for (sit = g[i].begin(); sit != g[i].end(); sit++)
			out << *sit << " ";
		out << "#" << endl;
	}
	out.close();
	
	// degree distribution
	map<int,int> dist;
	for (int i = 0; i < gnum; i++) {
		dist[g[i].size()]++;
	}
	count = 0;
	map<int,int>::reverse_iterator riter;
	riter = dist.rbegin();
	while (riter!=dist.rend()) { //&&count<1000
		cout << riter->first << "\t" << riter->second << endl;
		count++;
		riter++;
	}
}


// transform power law graph generated by R to graph
// in-degree follow power law distribution
static void relisttogra(char* filename) {
	ifstream in(filename);
	string buf;
	getline(in, buf);

	int src, trg, s, t;
	int index = -1;
	int id = 0;
	map<int,int> idmap;
	map<int,int>::iterator iter;
	map<int,set<int> > graph;
	int maxid = -1;
	
	while (getline(in, buf)) {
		index = buf.find_first_of(" ");
		src = str2int(buf.substr(0,index));
		buf = buf.substr(index+1);
		trg = str2int(buf);
		swap(src,trg);
		graph[src].insert(trg);
		if (src>maxid) maxid=src;
		if (trg>maxid) maxid = trg;
	}
	
	string filestr(filename);
	filestr = filestr.substr(0,filestr.find_last_of("."));
	string filestr1 = filestr + ".gra";
	ofstream out(filestr1.c_str());
	map<int,set<int> >::iterator mit;
	set<int>::iterator sit;
	out << "graph_for_reach" << endl;
	out << maxid << endl;
	index = 0;
	for (mit = graph.begin(); mit != graph.end(); mit++,index++) {
		while (mit->first > index) {
			out << index << ": #" << endl;
			index++;
		}
		out << index << ": ";
		for (sit = mit->second.begin(); sit != mit->second.end(); sit++) {
			out << *sit << " ";
		}
		out << "#" << endl;
	}
	for (; index < maxid; index++)
		out << index << ": #" << endl;
	out.close();
}


// vector<string>& split(const string &s, char delim, vector<string> &elems) {
	// stringstream ss(s);
	// string item;
	// while(getline(ss, item, delim)) {
		// elems.push_back(item);
	// }
	// return elems;
// }

// vector<string> split(const string &s, char delim) {
	// vector<string> elems;
	// return split(s, delim, elems);
// }

// for topology.txt
static void topologytodirectedgra(char* filename) {
	ifstream in(filename);
	string buf;
	int gsize;

	map<int,set<int> > gvec;
	map<int,int> idmap;
	set<int>::iterator sit;
	
	int src, trg, edgenum=0, id=0;
	int index = -1, min=1000000000, max=-100;
	while (getline(in, buf)) {
		if (buf[0]==' ' || buf[0]=='#') continue;
		vector<string> list = split(buf, '\t');
		if (list.size()<2) {
			list = split(buf, ' ');
		}
		assert(list.size()>=2);
		src = atoi(list[0].c_str());
		trg = atoi(list[1].c_str());
		if (src==trg) continue;
		if (idmap.find(src) == idmap.end()) {
			idmap[src] = id;
			id++;
		}
		if (idmap.find(trg) == idmap.end()) {
			idmap[trg] = id;
			id++;
		}
		gvec[idmap[src]].insert(idmap[trg]);
		edgenum++;
		//gvec[trg].insert(src);
	}
	in.close();
	cout << "edgenum=" << edgenum << " #V=" << idmap.size() << endl;
	
	// output to file
	string filestr(filename);
	filestr = filestr.substr(0,filestr.find_last_of("."));
	filestr = filestr + "_d.n3";
	cout << "output: " << filestr << endl;
	ofstream out(filestr.c_str());
	map<int,set<int> >::iterator mit;
	for (int i = 0; i < idmap.size(); i++) {
		for (sit = gvec[i].begin(); sit != gvec[i].end(); sit++)
			out << "<" << i << "> <link> <" << *sit << "> ." << endl;
	}
	
/*	
	filestr = filestr.substr(0,filestr.find_last_of("."));
	filestr = filestr + "_d.gra";
	ofstream out(filestr.c_str());
	out << "graph_for_distance" << endl;
	out << idmap.size() << endl;
	map<int,set<int> >::iterator mit;
	for (int i = 0; i < idmap.size(); i++) {
		out << i << ": ";
		for (sit = gvec[i].begin(); sit != gvec[i].end(); sit++)
			out << *sit << " ";
		out << "#" << endl;
	}
*/
	out.close();
}

static void convert(char* filename) {
	ifstream in(filename);
	string buf;
	getline(in, buf);
	
	int idx, begin, end, sid=0, tid=0, gsize;
	getline(in, buf);
	istringstream(buf) >> gsize;
	string sub;
	vector<set<int> > graph = vector<set<int> >(gsize,set<int>());
	while (getline(in, buf)) {
		begin = buf.find(":");
		end = buf.find_last_of("#");
		if (end-begin-2<1) {
			sid++;
			continue;
		}
		buf = buf.substr(begin+2,end-begin-2);
		vector<string> neighbors = split(buf, ' ');
		for (int i = 0; i < neighbors.size(); i++) {
			if (neighbors[i]=="") continue;
			tid = atoi(neighbors[i].c_str());
			graph[sid].insert(tid);
			graph[tid].insert(sid);
		}
		++sid;
	}
	in.close();
	
	string filestr(filename);
	filestr = filestr.substr(0,filestr.find_last_of("_"));
	string n3filename = filestr + ".n3";
	cout << "output: " << n3filename << endl;
	set<int>::iterator sit;
	ofstream n3out(n3filename.c_str());
	for (int i = 0; i < gsize; i++) {
		for (sit = graph[i].begin(); sit != graph[i].end(); sit++)
			n3out << "<" << i << "> <link> <" << *sit << "> ." << endl;
	}
	n3out.close();
	
	string grafile = filestr + ".gra";
	cout << "output: " << grafile << endl;
	ofstream gout(grafile.c_str());
	gout << "graph_for_distance" << endl;
	gout << gsize << endl;
	for (int i = 0; i < gsize; i++) {
		gout << i << ": ";
		for (sit = graph[i].begin(); sit != graph[i].end(); sit++)
			gout << *sit << " ";
		gout << "#" << endl;
	}
	gout.close();
}


static void networklisttodirectedgra(char* filename) {
	ifstream in(filename);
	string buf;
	int gsize;

	vector<set<int> > gvec = vector<set<int> >(128, set<int>());
	set<int>::iterator sit;
	
	int src, trg, edgenum=0;
	int index = -1, min=1000000000, max=-100;
	while (getline(in, buf)) {
	//	index = buf.find_first_of(" ");
	//	src = str2int(buf.substr(0,index));
	//	trg = str2int(buf.substr(index+1));
		if (buf[0]==' ' || buf[0]=='#') continue;
		vector<string> list = split(buf, ',');
		if (list.size()<2) {
			list = split(buf, ' ');
		}
		assert(list.size()>=2);
		src = atoi(list[0].c_str())-1;
		trg = atoi(list[1].c_str())-1;
//		cout << "src=" << src << " trg=" << trg << endl;
		if (src==trg) continue;
		if (src<min) min=src;
		if (trg<min) min=trg;
		if (src>max) max=src;
		if (trg>max) max=trg;
		
		if (max>=gvec.size()) {
//			cout << "resize" << endl;
			gvec.resize(max+2048,set<int>());
		}
//		cout << "size=" << gvec.size() << endl;
		gvec[src].insert(trg);
		edgenum++;
		//gvec[trg].insert(src);
		if (edgenum%500000==0)
			cout << edgenum << endl;
	}
	in.close();
	gvec.resize(max+1);
	cout << "edgenum=" << edgenum  << " min=" << min << " max=" << max << endl;
	
	// output to file
	string filestr(filename);
	/*
	filestr = filestr.substr(0,filestr.find_last_of("."));
	filestr = filestr + "_d.txt";
	ofstream out(filestr.c_str());
	*/

	
	filestr = filestr + "_d.gra";
	ofstream out(filestr.c_str());
	out << "graph_for_distance" << endl;
	if (min==0)
		out << gvec.size() << endl;
	else
		out << gvec.size()-1 << endl;
	
	if (min==0) {
		for (int i = 0; i < gvec.size(); i++) {
			
			out << i << ": ";
			for (sit = gvec[i].begin(); sit != gvec[i].end(); sit++) 
				out << *sit << " ";
			out << "#" << endl;
			
			/*
			for (sit = gvec[i].begin(); sit != gvec[i].end(); sit++) 
				out << i << "\t" << *sit << endl;
			*/
		}
	}
	else {
		for (int i = 1; i < gvec.size(); i++) {
			/*
			out << i-1 << ": ";
			for (sit = gvec[i].begin(); sit != gvec[i].end(); sit++) 
				out << *sit-1 << " ";
			out << "#" << endl;
			*/
			for (sit = gvec[i].begin(); sit != gvec[i].end(); sit++) 
				out << i-1 << "\t" << *sit-1 << endl;
		}
	}
	out.close();
}

int main(int argc, char* argv[]) {
	int gs = -1;
	char* filename;
	char* dagfilename;
	int pcode = -1;
	int i = 1;
	double density = 5.0;
	while (i < argc) {
		if (strcmp("-g", argv[i]) == 0) {
			i++;
			gs = atoi(argv[i++]);
		}
		else if (strcmp("-density", argv[i]) == 0) {
			i++;
			density = atof(argv[i++]);
		}
		else if (strcmp("-c", argv[i]) == 0) {
			i++;
			pcode = atoi(argv[i++]);
		}
		else if (strcmp("-d", argv[i]) == 0) {
			i++;
			dagfilename = argv[i++];
		}
		else {
			filename = argv[i++];
		}
	}
	if (dagfilename!=NULL)
		cout << "pcode=" << pcode << "\tmapfilename=" << dagfilename << "\tfilename=" << filename << endl;
	
	switch(pcode) {
		case 0: edgelistToGraph(filename, gs); break;
		case 1: readEdgeList(filename); break;
		case 2: mapauthors(dagfilename, filename); break;
		case 3: readpowerlawedges(filename); break;
		case 4: tranformtradedata(filename); break;
		case 5: mapcountrycode(dagfilename,filename); break;
		case 6: siftograph(filename); break;
		case 7: remapgraph(dagfilename, filename); break;
		case 8: graphToEdgelist(filename); break;
		case 9: PowerlawEdgelistToGraphForBarabasiGraph(filename); break;
		case 10: translate(filename); break;
		case 11: elisttogra(filename); break;
		case 12: subgraphToEdgelist(filename); break;
		case 13: degdist(filename); break;
		case 14: topicmatrixToGraph(filename, gs); break;
		case 15: relisttogra(filename); break;
		case 16: pywebgraph2gra(filename, gs, density); break;
		case 17: networklisttodirectedgra(filename); break;
		case 18: topologytodirectedgra(filename); break;
		case 19: convert(filename); break;
	}
}

