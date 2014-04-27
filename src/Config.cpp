#ifndef __CONFIG_H
#define __CONFIG_H

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <set>
#include <map>
#include <list>
#include <deque>
#include <algorithm>
#include <utility>
#include <limits>
#include <cmath>
#include <cstdlib>
#include <cassert>
#include <ext/hash_map>
#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>

using namespace std;
using namespace __gnu_cxx;

typedef unsigned int uint;
typedef unsigned long ulong;
typedef unsigned short ushort;
typedef ushort degreetype;
typedef ushort coresizetype;
typedef map<int,int> HighDegreeMap; // used to explicitly store a small portion of high degree nodes appearing in power-law graph

#define USHORT_NULL	0xFFFF
#define NN_VAL	 1e-20
#define MAX_VAL numeric_limits<int>::max()
#define MIN_VAL numeric_limits<int>::min()
#define USHORT_MAX numeric_limits<unsigned short>::max()
#define DEGREE_MAX numeric_limits<degreetype>::max()

struct routeEntry {
	unsigned int coreid : 16;
	unsigned int portid : 16;
	routeEntry() {
		coreid = 0;
		portid = 0;
	}
	routeEntry(unsigned int _coreid, unsigned int _portid)
		: coreid(_coreid), portid(_portid) {
	}
};

struct nodeRouteTable {
	routeEntry* outrouteitems; // sorted based on their distance to essential core nodes
	ushort* outitemindex; // indexing routeitems with different distance
	routeEntry* inrouteitems;
	ushort* initemindex;
	nodeRouteTable() {
		outrouteitems = NULL;
		outitemindex = NULL;
		inrouteitems = NULL;
		initemindex = NULL;
	}
	~nodeRouteTable() {
		if (outrouteitems!=NULL)
			delete[] outrouteitems;
		if (outitemindex!=NULL)
			delete[] outitemindex;
		if (inrouteitems!=NULL)
			delete[] inrouteitems;
		if (initemindex!=NULL)
			delete[] initemindex;
	}
};

#endif
