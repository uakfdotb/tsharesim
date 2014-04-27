#ifndef LPPATH_H
#define LPPATH_H

struct PathPoint;
class LPPath;
class TaxiPath;

class LPPath : public PathHandler {
public:
	LPPath(TaxiPath *taxiPath, int M, vector<PathPoint *> &points);
	~LPPath();
	
	virtual vector<PathPoint *> findPath(vertex *curr);
	
	//LPPath currently does everything in findPath
	virtual void init() {};
	
private:
	//stores the shortest distance between every pair of points
	//the first index uses zero to represent the current location,
	// and one as the first point
	//the second index uses zero as the first point, since the taxi
	// will never go _to_ the current point, only away from it
	double **costMatrix;
	
	//the earliest time that each point may be scheduled
	//earliest[0] is for the first point
	double *earliest;
	
	//the latest time that each point may be scheduled
	//latest[0] is for the current location
	double *latest;
	
	//this is used to put a non-linear constraint into the solver
	//see the sixth constraint
	//also, note that this matches the costMatrix (first index includes
	// current location)
	double  **m;
	
	//the corresponding pair for each dropoff point
	int *pairIndex;

	//makes sure that d is non-zero
	//iteration determines how much randomness to introduce on d if non-zero
	// this was introduced to prevent the solver from getting stuck
	double getNonZero(double d, int iteration);
};

#endif
