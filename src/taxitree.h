#ifndef TAXITREE_H
#define TAXITREE_H

#include "location.h"
#include "taxiPath.h"

class ShortestPath;
struct vertex;

class taxiTree
{
  public:
	int passengers; //number of passengers currently in the taxi
	int hashLong, hashLat; //used to identify our index in the grid
	double distToDest; //distance remaining to nextVertex
	
	vertex *previousVertex; //the vertex that we passed last
	vertex *nextVertex; //the vertex that we're heading to
	queue<vertex *> targetPath; //the path that we're current following
	
	location current; //(x, y) where we are currently at, used for (hashLong, hashLat)
	bool randDest; //true if destination was random
	
	struct timeval startTime;
	struct timeval endTime;
  	
  	TaxiPath *taxiPath;
 
	taxiTree() {}
	~taxiTree() {}
	
	void setInitialPosition(ShortestPath *shortestPath, vertex *vertex); //sets the initial position of the taxi
	double distance(double x1, double y1, double x2, double y2); //returns Euclidean distance between two points
	void setDestination(vertex *newVertex); //sets taxi's destination
	void randomDestination(); //taxi selects a random destination to follow for cruising
	void updateLocation(); //updates current location of taxi by D (speed)
	double value(vertex *source, vertex *dest); //determines the shortest path that includes the request (source, dest), or -1 if no such path exists
	void cancel(); //notification that the last request sent to value() has been assigned to another taxi
	void push(); //notification that the last request sent to value() has been assigned to _this_ taxi
	int getNumberNodes(); //some number for statistical purposes approximating memory used; or simply the number of requests
	void printPoints() { taxiPath->printPoints(); }
	bool dynamicConstraints(); //whether this taxitree supports dynamic constraints
	void setConstraints(double pickup_constraint, double service_constraint);
};

#endif
