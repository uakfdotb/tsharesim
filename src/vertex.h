#ifndef VERTEX_H
#define VERTEX_H

struct vertex;

struct vertex {
	int id;
	vector<vertex *> neighbors;
	
	double x;
	double y;
};

#endif
