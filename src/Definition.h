#ifndef DEFINITION_H_
#define DEFINITION_H_

#include <iostream>
#include <limits>


#define UINT_INFINITY 1073741823       //1024*1024*1024-1
#define INT_INFINITY numeric_limits<int>::max()
#define FLT_INFINITY numeric_limits<float>::max()
#define DBL_INFINITY numeric_limits<double>::max()

typedef unsigned int VertexID;
typedef unsigned int EdgeID; 
typedef unsigned int Rank;
typedef float Weight;

struct Edge {
	Weight weight;
	VertexID source : 30;
	bool flaga : 1;    // forward label;
	bool flagb : 1;	
	VertexID target : 30;
	bool flagc : 1;    // forward label;
	bool flagd : 1;	
};
typedef vector<Edge> EdgeList;


struct Vertex{
	unsigned int rank: 30;
	bool flaga: 1;
	bool flagb: 1;
	VertexID id: 30;
	bool flagc: 1;
	bool flagd: 1;
	
	Vertex() {
		flaga = false;
		flagb = false;
		flagc = false;
		flagd = false;
	};
	
	Vertex(int r, VertexID i): rank(r), id(i) {
		flaga = false;
		flagb = false;
		flagc = false;
		flagd = false;		
	};
};
typedef vector<Vertex> VertexList;

struct ShortCutEdge{
	VertexID target: 30;
	bool flaga: 1;
	bool flagb: 1;
	vector<VertexID> innerIDs; // store the intermediate edges;
	Weight weight;
	
	
	ShortCutEdge() {
		flaga = false;
		flagb = false;
	};
};
typedef vector<ShortCutEdge> ShortCuts;

struct Label{
	double distance;
	VertexID id;
};
typedef vector<Label> LabelList;



#endif
