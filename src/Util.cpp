#include "Util.h"

int Util::binarysearch_index(const vector<int>& data, int begin, int end, int value) {
	int low = begin, high = end-1, midpoint = 0;
	if (high-low<10) {
		for (; low <= high; low++)
			if (data[low] == value)
				return low;
	}
	else {
		while (low <= high) {
			midpoint = low+(high-low)/2;
			if (value == data[midpoint]) 
				return midpoint;
			else if (value < data[midpoint])
				high = midpoint-1;
			else
				low = midpoint+1;
		}
	}
	return -1;
}

double Util::plogpq(double p, double q) {
	if (isnan(p) || isnan(q)) return MAX_VAL;
	if (p<NN_VAL && q<NN_VAL) return 0;
	if (p<NN_VAL) p=NN_VAL;
	if (q<NN_VAL) q=NN_VAL;
	return p*log2(p/q);
}

double Util::nlogp(int n, double p) {
	if (isnan(p)) return MIN_VAL;
	if (p<NN_VAL) p=NN_VAL;
	return n*1.0*log2(p);
}

double Util::logdist(int n, double p, double q) {
	if (isnan(p) || isnan(q)) return MIN_VAL;
	if (p<NN_VAL) p=NN_VAL;
	if (q<NN_VAL) q=NN_VAL;

	return n*1.0*(log2(p)-log2(q));
}

double Util::longlogdist(long n, long p, long q) {
	/*
	if (isnan(p) || isnan(q)) return MIN_VAL;
	if (p<NN_VAL) p=NN_VAL;
	if (q<NN_VAL) q=NN_VAL;
	*/

	return n*1.0*(log2(p)-log2(q));
}

void Util::printMap(const map<int,int>& data) {
	map<int,int>::const_iterator mit;
	cout << "Map" << endl;
	for (mit=data.begin(); mit != data.end(); mit++)
		cout << "<" << mit->first << "," << mit->second << "> ";
	cout << endl;
}

void Util::printLongMap(const map<long,int>& data) {
	map<long,int>::const_iterator mit;
	cout << "LongMap" << endl;
	for (mit=data.begin(); mit != data.end(); mit++)
		cout << "<" << mit->first << "," << mit->second << "> ";
	cout << endl;
}

void Util::printUTDVec(const UTDVec& tv) {
	vector<unsigned int>::const_iterator vit;
	cout << "UTDVec" << endl;
	for (int i = 0; i < tv.size(); i++) {
		cout << "v " << i << ":";
		for (vit = tv[i].begin(); vit != tv[i].end(); vit++)
			cout << *vit << " ";
		cout << endl;
	}
}

void Util::printTDVec(const TDVec& tv) {
	vector<int>::const_iterator vit;
	cout << "TDVec" << endl;
	for (int i = 0; i < tv.size(); i++) {
		cout << "v " << i << ":";
		for (vit = tv[i].begin(); vit != tv[i].end(); vit++)
			cout << *vit << " ";
		cout << endl;
	}
}

void Util::printSparseVec(const SparseVec& sv) {
	SparseVec::const_iterator sit;
	map<int,int>::const_iterator mit;
	cout << "SparseVec" << endl;
	for (int i = 0; i < sv.size(); i++) {
		cout << "from " << i << ": ";
		for (mit = sv[i].begin(); mit != sv[i].end(); mit++) 
			cout << mit->first << "[" << mit->second << "] ";
		cout << endl;
	}
}

void Util::printSparseVecVL(const SparseVecVL& sv) {
	SparseVecVL::const_iterator sit;
	map<int,long>::const_iterator mit;
	cout << "SparseVecVL" << endl;
	for (int i = 0; i < sv.size(); i++) {
		cout << "from " << i << ": ";
		for (mit = sv[i].begin(); mit != sv[i].end(); mit++) 
			cout << mit->first << "[" << mit->second << "] ";
		cout << endl;
	}
}

void Util::printSparseVecLong(const SparseVecLong& sv) {
	SparseVecVL::const_iterator sit;
	map<long,int>::const_iterator mit;
	cout << "SparseVecLong" << endl;
	for (int i = 0; i < sv.size(); i++) {
		cout << "from " << i << ": ";
		for (mit = sv[i].begin(); mit != sv[i].end(); mit++) 
			cout << mit->first << "[" << mit->second << "] ";
		cout << endl;
	}	
}

void Util::printReducedGraph(const ReducedGraph& rg) {
	ReducedGraph::const_iterator rit;
	vector<int>::const_iterator vit;
	cout << "ReducedGraph #vsize" << rg.size() << endl;
	for (rit = rg.begin(); rit != rg.end(); rit++) {
		cout << "vertex " << rit->first << ": ";
		for (vit = rit->second.begin(); vit != rit->second.end(); vit++)
			cout << *vit << " ";
		cout << "#" << endl;
	}
	cout << endl;
}

void Util::printVecSet(const vector<set<int> >& vs) {
	set<int>::const_iterator sit;
	for (int i = 0; i < vs.size(); i++) {
		cout << i << " : ";
		for (sit = vs[i].begin(); sit != vs[i].end(); sit++) 
			cout << *sit << " ";
		cout << endl;
	}
	cout << endl;
}

void Util::printPairVec(const vector<map<int,bool> >& vm) {
	map<int,bool>::const_iterator mit;
	for (int i = 0; i < vm.size(); i++) {
		cout << "v " << i << " : ";
		for (mit = vm[i].begin(); mit != vm[i].end(); mit++)
			cout << "[" << mit->first << "," << mit->second << "] ";
		cout << endl;
	}
	cout << endl;
}

void Util::printMap(const map<uint,int>& data) {
	map<uint,int>::const_iterator mit;
	for (mit = data.begin(); mit != data.end(); mit++) {
		cout << "{" << mit->first << "," << mit->second << "} ";
	}
	cout << endl;
}

void Util::printMapSet(const map<int,set<int> >& data) {
	map<int,set<int> >::const_iterator msit;
	set<int>::const_iterator sit;
	for (msit = data.begin(); msit != data.end(); msit++) {
		cout << msit->first << ": ";
		for (sit = msit->second.begin(); sit != msit->second.end(); sit++)
			cout << *sit << " ";
		cout << " #" << endl;
	}
}

void Util::printVecPair(const vector<pair<int,int> >& data){
	vector<pair<int,int> >::const_iterator it;
	for (it = data.begin(); it != data.end(); it++)
		cout << "[" << it->first << "," << it->second << "] ";
	cout << endl;
}

void Util::btiming(struct timeval& before_time) {
	gettimeofday(&before_time, NULL);
}

float Util::timing(struct timeval& before_time) {
	struct timeval after_time;
	gettimeofday(&after_time, NULL);
	float run_time = (after_time.tv_sec - before_time.tv_sec)*1000.0 + 
		(after_time.tv_usec - before_time.tv_usec)*1.0/1000.0;
	return run_time;
}

string Util::formatTime(float mstime) {
	int mtime = (int)mstime;
	string timestr = "";
	int unit[] = {60, 60, 60, 1000};
	float quality[] = {0, 0, 0, 0};
	string ustr[] = {"Hours", "Mins", "Secs", "Mss"};
	
	int count = 3;
	while (count>=0) {
		quality[count] = mtime%unit[count];
		mtime = (int)(mtime/unit[count]);
		count--;
	}
	
	timestr = to_string(quality[0])+":"+to_string(quality[1])+":"+to_string(quality[2])
			+" (" +to_string(quality[3])+"ms)";
	return timestr;
}

void Util::uniquevec(vector<int>& vec) {
	sort(vec.begin(), vec.end());
	vector<int>::iterator it = unique(vec.begin(),vec.end());
	vec.resize(it-vec.begin());
}

void Util::uniquevecpair(vector<pair<int,int> >& vec) {
	sort(vec.begin(), vec.end());
	vector<pair<int,int> >::iterator it = unique(vec.begin(),vec.end());
	vec.resize(it-vec.begin());
}

void Util::binaryinsert(vector<pair<int,int> >& vec, pair<int,int> a) {
	vector<pair<int,int> >::iterator it =
		lower_bound(vec.begin(), vec.end(), a, pair_comp());
	if (it == vec.end() || (*it) != a)
		vec.insert(it, a);
}

void Util::complement(ReducedGraph& rg) {
	ReducedGraph newrg;
	ReducedGraph::iterator rit;
	vector<int>::iterator vit;
	for (rit = rg.begin(); rit != rg.end(); rit++) {
		for (vit = rit->second.begin(); vit != rit->second.end(); vit++) {
			newrg[rit->first].push_back(*vit);
			newrg[*vit].push_back(rit->first);
		}
	}
	for (rit = newrg.begin(); rit != newrg.end(); rit++) {
		uniquevec(rit->second);
	}
	ReducedGraph().swap(rg);
	rg = newrg;
}

bool Util::isnan(double x) {
	return (x!=x);
}

void Util::strTrimRight(string& str) {
	string whitespaces(" \t");
	int index = str.find_last_not_of(whitespaces);
	if (index != string::npos) 
		str.erase(index+1);
	else
		str.clear();
}

//////////////////////////////////////////////////////////////////////////////
//
// process_mem_usage(double &, double &) - takes two doubles by reference,
// attempts to read the system-dependent data for a process' virtual memory
// size and resident set size, and return the results in KB.
//
// On failure, returns 0.0, 0.0
void Util::process_mem_usage(double& vm_usage, double& resident_set) {
   using std::ios_base;
   using std::ifstream;
   using std::string;

   vm_usage     = 0.0;
   resident_set = 0.0;

   // 'file' stat seems to give the most reliable results
   pid_t curpid = getpid();
   string filepath = string("/proc/") + to_string(curpid) + string("/stat");
//   cout << "filepath=" << filepath << endl;
   ifstream stat_stream(filepath.c_str(),ios_base::in);
   if (!stat_stream) {
		cerr << "Cannot open " << filepath << endl;
		return;
   }
   // dummy vars for leading entries in stat that we don't care about
   //
   string pid, comm, state, ppid, pgrp, session, tty_nr;
   string tpgid, flags, minflt, cminflt, majflt, cmajflt;
   string utime, stime, cutime, cstime, priority, nice;
   string O, itrealvalue, starttime;

   // the two fields we want
   unsigned long vsize;
   long rss;

   stat_stream >> pid >> comm >> state >> ppid >> pgrp >> session >> tty_nr
               >> tpgid >> flags >> minflt >> cminflt >> majflt >> cmajflt
               >> utime >> stime >> cutime >> cstime >> priority >> nice
               >> O >> itrealvalue >> starttime >> vsize >> rss; // don't care about the rest

   long page_size_kb = sysconf(_SC_PAGE_SIZE) / 1024; // in case x86-64 is configured to use 2MB pages
   vm_usage     = vsize / 1024.0;
   resident_set = rss * page_size_kb;
//   cout << "vm_usage=" << vm_usage << " KB" << endl;
//   cout << "resident_set=" << resident_set << endl;
}

