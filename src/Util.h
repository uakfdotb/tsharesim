#ifndef UTIL_H_
#define UTIL_H_

#include "Config.cpp"
#include "bit_vector.h"
#include "PerformanceTimer.h"

typedef map<int,vector<int> > 	ReducedGraph;
typedef map<int,vector<pair<int,int> > > WeightedGraph;
typedef map<int,map<int,int> >  SparseMatrix;
typedef vector<map<int,double> > SparseDVec;
typedef vector<map<int,int> > 	SparseVec;
typedef vector<map<int,long> >	SparseVecVL;
typedef vector<map<long,int> > 	SparseVecLong;
typedef vector<vector<int> > 	TDVec;
typedef vector<vector<unsigned int> > 	UTDVec;
typedef vector<set<int> >		VecSet;
typedef hash_map<long,pair<int,int> > MapPair;
typedef map<int,vector<vector<double> > > Map2Vec;

struct pair_comp{
	bool operator() (const pair<int,int>& lhs, const pair<int,int>& rhs) const {
		if (lhs.first==rhs.first) return lhs.second<rhs.second;
		return lhs.first<rhs.first;
	}
};

template<class T>
struct val_comp {
	bool operator() (const T& lhs, const T& rhs) const {
		return lhs>rhs;
	}
};

template<class T>
struct index_comp {
	vector<T>& _x;
	
	index_comp(vector<T>& x): _x(x) {};
	bool operator() (const int a, const int b) const {
		return _x[a]>_x[b];
	}
};

struct eindex {
	long eid;
	int inindex;
	int outindex;
};

template<class T>
inline string to_string(const T& t) {
	stringstream ss;
	ss << t;
	return ss.str();
}

class Util {
	public:
		static int binarysearch_index(const vector<int>& data, int begin, int end, int value);
	
		static double plogpq(double p, double q);
		static double nlogp(int n, double p);
		static double logdist(int n, double p, double q);
		static double longlogdist(long n, long p, long q);
		
		static void printMap(const map<int,int>& data);
		static void printLongMap(const map<long,int>& data);
		static void printTDVec(const TDVec& tv);
		static void printUTDVec(const UTDVec& tv);
		static void printReducedGraph(const ReducedGraph& rg);
		static void printSparseVec(const SparseVec& sv);
		static void printSparseVecVL(const SparseVecVL& sv);
		static void printSparseVecLong(const SparseVecLong& svl);
		static void printVecSet(const vector<set<int> >& vs);
		static void printPairVec(const vector<map<int,bool> >& vm);
		static void printMap(const map<uint,int>& data);
		static void printMapSet(const map<int,set<int> >& data);
		static void printVecPair(const vector<pair<int,int> >& data);
		static void btiming(struct timeval& before_time);
		static float timing(struct timeval& before_time);
		static string formatTime(float mstime);
		static bool isnan(double x);
		static void uniquevec(vector<int>& vec);
		static void uniquevecpair(vector<pair<int,int> >& vec);
		static void binaryinsert(vector<pair<int,int> >& vec, pair<int,int> a);
		static void complement(ReducedGraph& rg);
		static void strTrimRight(string& str);
		static void process_mem_usage(double& vm_usage, double& resident_set);
};

#endif
