#ifndef __GRAPHLINK_H__
#define __GRAPHLINK_H__

#include "Contraction.h"
#include "Query.h"
#include "Dijkstra.h"
#include "PerformanceTimer.h"
#include "GraphUtil.h"
#include "HD.h"

#include <iostream>
#include <time.h>

using namespace std;

void printUsage();
void printParameters(const map<string, string>& cmdLineArgs);
void keepResult(const char* resultFileName, vector<string>& results);
void readFile(ifstream& in, vector<vector<double> >& data);
void checkQuery(Graph& graph, string file_name, int num_query, vector<VertexID>& query, bool save_para);
void readHierarchies(Graph& graph, string hierarchy_file);
void checkWritingCorrectness(Graph& graph, string hierarchy_file);
void checkHierarchies(Graph& graph, string file_name, bool save_para, bool gp);
void readIDMap(string filename, vector<VertexID>& id_map);
bool isNumeric(string input);

#endif
