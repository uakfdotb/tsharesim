#ifndef BRANCHPATH_H
#define BRANCHPATH_H

struct PathPoint;
class BranchPathNode;
class BranchPath;
class TaxiPath;

class BranchPathNode {
public:
	int N; //total number of points (pickup and dropoff)
	int level; //number of pickup and dropoff points in this path
	double length, bound; //length to parent and point's bound
	PathPoint *vertex;
	
	BranchPathNode *parent;
	vector<BranchPathNode *> children;
	
	//identifies which points are in the path that this
	// node identifies for fast retrieval
	//also, for pickup points, stores the distance from
	// the root to the point going along the route so
	// that the corresponding dropoff point can use it
	// for feasibility checking
	//for dropoff points, the double is simply -1
	//maps from
	// pairIndex   for pickup points
	// pairIndex+1 for dropoff points
	map<int, double> positions;
	
	BranchPath* branch;
	
	//whether or not this point is feasible
	bool feasible;
	
	BranchPathNode(PathPoint *vert, int N, BranchPath* branch);
	BranchPathNode(BranchPathNode* parent, BranchPath* branch, PathPoint *vert);
	~BranchPathNode();
	bool pathContains(PathPoint *test); //whether the path identified by this node contains the given point
};

class BranchPathComparison
{
public:
	BranchPathComparison() {}
	bool operator() (const BranchPathNode* lhs, const BranchPathNode* rhs) const {
		return lhs->bound > rhs->bound;
	}
};

class BranchPath : public PathHandler {
protected:
	double *minEdge; //each vertex -- used by method bound( )
	priority_queue<BranchPathNode*, vector<BranchPathNode*>, BranchPathComparison> q;

public:
	BranchPath(TaxiPath *taxiPath, int M, vector<PathPoint *> &points);
	~BranchPath();
	
	virtual vector<PathPoint *> findPath(vertex *curr);
	virtual void init( );
	
private:
	double bound(BranchPathNode* n);
	void execute();
};

#endif
