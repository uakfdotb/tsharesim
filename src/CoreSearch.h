#ifndef _CORE_SEARCH_H_
#define _CORE_SEARCH_H_

#include "GraphUtil.h"

#define CS_STAT
//#define CS_DEBUG

class CoreSearch {
	private:
		Graph& g;
		int radius, coresize, gsize;
		bit_vector* corenodes;
		coresizetype** corepaths; // explicitly record the shortest path, first element is the length of path, the number of corenodes should be no more than 65535 (0xFFFF)
		nodeRouteTable *routetables;
		
		// helper data structures
		bool estimate;
		vector<int> nodemap;
		vector<vector<unsigned int> > combination;
		int ref, searchspace, comparetimes, nodeordertype;
		int *prev, *dist, *mark, *pque;
		ushort* preport;
		double vm_usage, resident_set;
		unsigned long num_entries;
		
		// for performance testing
		int noncorequeries, corehitqueries;
		float runtime1, runtime2, runtime;
		struct timeval after_time, before_time, end_time, begin_time;
		
	public:
		CoreSearch(Graph& graph, int radius, int corenum, int ordertype, bool _estflag);
		CoreSearch(Graph& graph, int radius, int corenum);
		~CoreSearch();
		void corelabeling();
		void coreLabelingNode(int vid, vector<int>& que, vector<vector<vector<routeEntry> > >& tmp_table, bool isout);
		void createLabels();
		void createLabels(vector<int>& srcs, vector<int>& trgs, string interResultFilename);
		void convertTables(vector<vector<vector<routeEntry> > >& tmp_outtable, vector<vector<vector<routeEntry> > >& tmp_intable);
		vector<int> shortestpath(int src, int trg, int& distance);
		vector<int> shortestpath(int src, int trg);
		bool test_shortestpath(int src, int trg);
		int distance(int src, int trg);
		bool test_distance(int src, int trg);
		
		// utility functions
		int getSearchspace();
		int getComparetimes();
		unsigned long getNumRouteEntries();
		unsigned long getCorepathsize();
		double getPeakmemusage();
		void initHelper();
		void destroyHelper();
		void setEstimate(bool flag);
		unsigned int computeMatrixIndex(int i, int j);
		void displayCorepaths(ostream& out);
		void displayRouteTable(ostream& out);
		void displayRouteTableByNode(int nid, ostream& out);
		void loadIntermediateResults(istream& in);
		void writeIntermediateResults(ostream& out);
};

#endif
