#ifndef _GRAPH_UTIL_H_
#define _GRAPH_UTIL_H_

#include "Graph.h"
#include <queue>

//#define GU_DEBUG

class GraphUtil {
	public:
		// benchmarks
		static int  BiBFSDist_ptr(Graph& graph, int src, int trg, int* dist, int* que, int& ref, int radius, int& searchspace);
		static int  BiBFSDist(Graph& graph, int src, int trg, vector<int>& dist, vector<int>& que, int& ref, int radius, int& searchspace);
		static int  BiBFSDist_simple(Graph& graph, int src, int trg, vector<int>& dist, vector<int>& que, int& ref, int radius, int& searchspace);
		static int  BFSDist(Graph& graph, int src, int trg, vector<int>& dist, vector<int>& que, int& ref, int radius, int& searchspace);
		static vector<int>  BiBFSPath(Graph& graph, int src, int trg, vector<int>& dist, vector<int>& que, vector<int>& prev, int& ref, int radius, int& searchspace);
		static vector<int>  BFSPath(Graph& graph, int src, int trg, vector<int>& dist, vector<int>& que, vector<int>& prev, int& ref, int radius, int& searchspace);
		
		// utility functions
		static bool checkUndirectedGraph(Graph& graph);
		static vector<int> coreNodeSelection_degree(Graph& graph, int corenum, bit_vector* corenodes);
		static vector<int> coreNodeSelectionWithNodeMap(Graph& graph, int corenum, int type, bit_vector* corenodes, vector<int>& nodemap);
		static void computePairwiseDist(Graph& graph, bit_vector* corenodes, int radius, map<int,int>& coreindex, vector<vector<int> >& coredist);
		static int generateQueriesByNode(Graph& graph, int src, int radius, int num_nodes, vector<int>& dist, int& ref,
				vector<int>& srcvec, vector<int>& trgvec);
		static void generateQueriesByDistanceAndOutput(Graph& graph, int radius, int querynum, ostream& out);
		static void generateQueriesByDistance(Graph& graph, int radius, int querynum, vector<int>& srcvec, vector<int>& trgvec);
		static void generateCombination(int radius, vector<vector<unsigned int> >& combination);
		static void DFSExpand(Graph& graph, int vid, bit_vector* visited, bit_vector* corenodes, int& size);
		static vector<int> degreeDecomposition(Graph& graph, int coresize);
		static int generateRandQueriesByNode(Graph& graph, int src, int num_nodes, int radius, vector<int>& dist, int& ref,
				vector<int>& srcvec, vector<int>& trgvec, vector<int>& distvec);
		static void generateRandomQueries(Graph& graph, int radius, int querynum, ostream& out);
		static int generateRandQueriesByNode1(Graph& graph, int src, int num_nodes, int radius, vector<int>& dist, int& ref,
				vector<int>& srcvec, vector<int>& trgvec, vector<int>& distvec);
		static void generateRandomQueries1(Graph& graph, int radius, int querynum, ostream& out);
		static void generateRandomBoundedQueries(Graph& graph, int radius, int querynum, ostream& out);
		static void readCpuMemInfo(vector<string>& cpumem);
		static double effective90diameter(Graph& graph, int seednum);
		static void BFSDistribute(Graph& graph, int src, vector<int>& dist, int& ref, int radius, vector<int>& distribute);
		
		// for hub labeling
		static void generateRandomQuery(int, int, vector<VertexID>&);
		static void generateRandomBoundedQuery(Graph&, int, int, vector<VertexID>&, vector<double>&);
		static void generateQuery(Graph&, int, int, ostream&, vector<VertexID>&);
};

#endif
