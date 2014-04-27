#ifndef CONTRACTION_H_
#define CONTRACTION_H_

#include <iostream>
#include <vector>
#include <queue>

#include "GraphUtil.h"

using namespace std;

//#define CONTRACTION_DEBUG


#define min(a,b) (a<b? a:b)


class Contraction{

private:
	Graph& graph;
	VertexList nodelist;
	unsigned int nodesize;
	
	vector<ShortCuts> tmpinshortcutlist; // for preprocessing;
	vector<ShortCuts> tmpoutshortcutlist;
	
	vector<LabelList> outlabellist;
	vector<LabelList> inlabellist;
	
	vector<vector<int> > pathIndex;
	vector<vector<VertexID> > pathList;	
	
public:
	Contraction(Graph&);
	virtual ~Contraction(void);
	
	bool dodgingDijkstra(VertexID, VertexID, VertexID, double, double&);	
	inline bool checkDuplicate(VertexID, VertexID, Weight);
	void removeEdge(int);
	void packShortcut(int, vector<VertexID>&, vector<double>&, vector<vector<EdgeID> >&, vector<int>&, vector<double>&);
	void addShortCut(VertexID, VertexID, VertexID, Weight, vector<EdgeID>&, vector<EdgeID>&);	
	int  computeEdgeDifference(int, vector<VertexID>&, vector<double>&, vector<vector<EdgeID> >&, vector<int>&, vector<double>&);
	void computeShortcuts(bool);
	
	vector<ShortCuts>& exportOutShortcut(void);	
	vector<ShortCuts>& exportInShortcut(void);	
	VertexList& exportNodeList(void); 

	void createLabel(VertexID, bool);
	void generateLabel(bool);
	void printLabel(void);
	void printPath(void);	
	
	void storePathInfo(int, int, vector<VertexID>&, vector<vector<VertexID> >&);	
	
	void checkShortCutCorrectness(void);
	void checkShortCutCompleteness(int);
	double runDijkstra(VertexID, VertexID);
	void runShortcutDijkstra(VertexID, vector<VertexID>&, vector<Weight>&);
	
	void checkLabelCorrectness(void);	
	void buildPathIndex(VertexID, vector<VertexID>&, vector<vector<int> >&, vector<VertexID>&, bool&);	
	void checkPathCorrectness(void);	
};

#endif
