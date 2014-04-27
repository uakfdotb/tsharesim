#ifndef ONETAXIPATH_H
#define ONETAXIPATH_H

class TaxiPath;
class PathHandler;
class ShortestPath;

struct PathPoint {
	int index; //index in points array
	vertex *vert;
	int type; //0 is starting point, 1 is dropoff point, 2 is stand-alone dropoff point (starting point already reached), -1 is root
	int pairIndex; //used to identify pairs, incremental from zero; not set for stand-alone dropoff point
	double remaining; //length remaining; measuerd from root for type=0 and type=2, from pair for type=1
};

class OneTaxiPath: public TaxiPath {
private:
	vector<PathPoint *> points;
	vector<PathPoint *> temp; //the points set when value() is called; points = temp after push() is called
	int nextPair;
	PathHandler *pathb;
	
	bool flag; //1 if value was just called, 0 otherwise

public:
	OneTaxiPath(ShortestPath *shortestPath, vertex *curr); //currNode is index of the vertex that taxi is currently at or heading towards
	~OneTaxiPath();
	
	virtual void moved(double distance); //notification that the taxi has moved a given distance
	virtual double value(vertex *curr, vertex *source, vertex *dest); //tests pushing (source, dest)
	virtual void cancel(); //cancels the value() call
	virtual void push(); //completes the pushing after value() is called
	virtual bool step(bool move); //next() becomes curr(); goes one vertex forward, when curr() is passed
	virtual queue<vertex *> next(); //returns next vertex that the taxi should go to after curr() is passed; -1 if no next
	virtual queue<vertex *> curr(vertex *curr); //returns the vertex that the taxi is currently at or heading towards
	
	virtual int size() { return points.size(); }
	virtual void printPoints(); //prints current point array
	virtual int getNumberNodes() { return points.size(); };
	
	//used with pplog (path point log), to directly feed cases into an algorithm
	virtual void test(vertex *curr, vector<PathPoint *> nPoints);
};

#endif
