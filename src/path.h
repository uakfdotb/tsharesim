#ifndef PATH_H
#define PATH_H

struct PathPoint;
struct vertex;
class PathNode;
class PathHandler;
class TaxiPath;

class PathHandler {
protected:
	double minTour; //length of the best path
	TaxiPath *taxiPath; //used for getting shortest path between vertices
	vector<PathPoint *> &points;
	int N; //number of vertices
	vector<PathPoint *> bestList; //sequence of vertices that represents the best path

public:
	PathHandler(TaxiPath *taxiPath, int M, vector<PathPoint *> &points);
	virtual ~PathHandler();
	
	virtual void reset(); //called after each call to value()
	virtual double length(PathPoint *v1, PathPoint *v2); //shortest distance from v1 to v2
	virtual double lengthVertex(vertex *vert, PathPoint *v2); //shortest distance from vert to v2
	virtual vector<PathPoint *> findPath(vertex *curr) = 0; //finds the shortest path, given that we're currently at curr
	virtual void init( ) = 0; //this function is useless; todo: get rid of it
	
	virtual double getMinTour() { return minTour; }
};

#endif
