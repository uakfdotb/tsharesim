#ifndef TAXIPATH_H
#define TAXIPATH_H

class PathPoint;
class PathBranchAndBound;
class ShortestPath;
struct vertex;

class TaxiPath {
protected:
	ShortestPath *shortestPathC;
	vertex *curr_vert;

public:
	TaxiPath(ShortestPath *shortestPath, vertex *curr); //currNode is index of the vertex that taxi is currently at or heading towards
	virtual ~TaxiPath();

	virtual void moved(double distance) = 0;
	virtual double value(vertex *curr, vertex *source, vertex *dest) = 0; //tests pushing (source, dest)
	virtual void cancel() = 0; //cancels the value() call
	virtual void push() = 0; //completes the pushing after value() is called
	virtual bool step(bool move) = 0; //next() becomes curr(); goes one vertex forward, when curr() is passed
	virtual queue<vertex *> next() = 0; //returns next vertex that the taxi should go to after curr() is passed; -1 if no next
	virtual queue<vertex *> curr(vertex *curr) = 0; //returns the vertex that the taxi is currently at or heading towards

	virtual int size() = 0;
	virtual void printPoints() = 0; //prints current point array
	virtual double shortestPath(vertex *a, vertex *b); //returns shortest path from source to dest
	virtual int getNumberNodes() = 0;

	virtual bool dynamicConstraints() { return false; } //whether this taxipath supports dynamic constraints
	virtual void setConstraints(double pickup_constraint, double service_constraint) {}
};

#endif
