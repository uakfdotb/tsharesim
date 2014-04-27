#ifndef _STATIC_GRAPH_H
#define _STATIC_GRAPH_H

#include "Util.h"
#include "Definition.h"

//#define GRAPH_DEBUG 

class HD;

// unweighted and directed graph
class Graph {
	private:
		VertexList nodelist;
		EdgeList inedgelist, outedgelist;
		vector<int> inedgeindex, outedgeindex;
		//vector<degreetype> inportlist, outportlist; // entity(x,y) represents the port id of vertex x in the neighbor list of vertex y
		unsigned int nodesize;
		
		ShortCuts inShortCutList, outShortCutList;
		vector<int> inShortCutIndex, outShortCutIndex;
		
		vector<LabelList> outLabelList;
		vector<LabelList> inLabelList;
		
		vector<vector<int> > pathIndex;
		vector<vector<VertexID> > pathList;
		
		friend class HD;
		
	public:
		Graph();
		Graph(int graphsize);
		Graph(istream& gio);
		//Graph(istream& gio, bool isedgelist);
		~Graph();
		void readGraph(istream& gio);
		//void readEdgeList(istream& gio);
		void writeGraph(ostream& gio) ;
		void printGraph();
		void printShortcut();
		void printRank();
		int num_vertices();
		int num_edges();
		VertexList& vertices();
		int in_degree(int nodeid);
		int out_degree(int nodeid);
		bool hasEdge(int src, int trg);	
		void sortEdges();
		//void reorderid(const vector<int>& nodemap);
		Graph& operator=(const Graph& graph);
		Vertex& operator[](const int id);
		bool isDirectedGraph();
		static vector<string>& split(const string &s, char delim, vector<string> &elems);
		static vector<string> split(const string &s, char delim);
		void clear();
		
		void insertShortcut(vector<ShortCuts>&, vector<ShortCuts>&);
		void insertNodeList(VertexList&);
		void insertLabelList(vector<LabelList>&, vector<LabelList>&);
		
		void insertInLabel(vector<vector<double> >&);
		void insertOutLabel(vector<vector<double> >&);
		
		LabelList& exportOutLabel(VertexID);
		LabelList& exportInLabel(VertexID);
		void sortLabel();
		
		vector<vector<int> >& exportPathIndex(void);
		vector<vector<VertexID> >& exportPathList(void);
		void insertPathIndex(vector<vector<double> >&);
		void insertPathList(vector<vector<double> >&);

		ShortCuts& exportinShortCutList();
		ShortCuts& exportoutShortCutList();
		vector<int>& exportinShortCutIndex();
		vector<int>& exportoutShortCutIndex();
		
		// utility functions
		inline int get_inneighbor(int vid, int portid) const {
			return (inedgelist.begin()+inedgeindex[vid]+portid)->target;
		}
		inline int get_outneighbor(int vid, int portid) const {
			return (outedgelist.begin()+outedgeindex[vid]+portid)->target;
		}
		
		inline int get_outedgeindex(int vid, VertexID target) const {
			for (int i = outedgeindex[vid]; i < outedgeindex[vid+1]; i++) {
				if (outedgelist[i].target == target) return i;
			}
		}
		inline Edge& getoutEdge(EdgeID eid) {
			return outedgelist[eid];
		}
		
		inline vector<Edge>::iterator firstoutneighbor(const int vid){
			return outedgelist.begin()+outedgeindex[vid];
		}
		inline vector<Edge>::iterator endofoutneighbor(const int vid){
			return outedgelist.begin()+outedgeindex[vid+1];
		}
		inline vector<Edge>::iterator firstinneighbor(const int vid) {
			return inedgelist.begin()+inedgeindex[vid];
		}
		inline vector<Edge>::iterator endofinneighbor(const int vid){
			return inedgelist.begin()+inedgeindex[vid+1];
		}
		
		inline ShortCuts::iterator firstoutShortCut(const int vid){
			return outShortCutList.begin()+outShortCutIndex[vid];
		}
		inline ShortCuts::iterator endofoutShortCut(const int vid){
			return outShortCutList.begin()+outShortCutIndex[vid+1];
		}
		inline ShortCuts::iterator firstinShortCut(const int vid){
			return inShortCutList.begin()+inShortCutIndex[vid];
		}
		inline ShortCuts::iterator endofinShortCut(const int vid) {
			return inShortCutList.begin()+inShortCutIndex[vid+1];
		}		
};
	

// traverse all nodes
#define forall_nodes(G,nid)  for(unsigned int nid=0; nid<G.num_vertices(); nid++)
// traverse all outgoing neighbor of specific node nid
#define forall_outneighbors(G,nid,outneighbor_ptr)  \
	for (vector<Edge>::iterator outneighbor_ptr=G.firstoutneighbor(nid); outneighbor_ptr!=G.endofoutneighbor(nid); outneighbor_ptr++)
// traverse all incoming neighbor of specific node nid
#define forall_inneighbors(G,nid,inneighbor_ptr)  \
	for (vector<Edge>::iterator inneighbor_ptr=G.firstinneighbor(nid); inneighbor_ptr!=G.endofinneighbor(nid); inneighbor_ptr++)	
// traverse all outgoing shortcuts of specific node nid;
#define forall_outshortcuts(G,nid,outshortcut_ptr)  \
	for (ShortCuts::iterator outshortcut_ptr=G.firstoutShortCut(nid); outshortcut_ptr!=G.endofoutShortCut(nid); outshortcut_ptr++)
// traverse all incoming neighbor of specific node nid
#define forall_inshortcuts(G,nid,inshortcut_ptr)  \
	for (ShortCuts::iterator inshortcut_ptr=G.firstinShortCut(nid); inshortcut_ptr!=G.endofinShortCut(nid); inshortcut_ptr++)
	
	
	
#endif
