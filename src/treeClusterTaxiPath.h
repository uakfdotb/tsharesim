#ifndef TREECLUSTERTAXIPATH_H
#define TREECLUSTERTAXIPATH_H

class TaxiPath;
class TreeClusterNode;
class ShortestPath;

class TreeClusterNode {
public:
	const static double theta = 0.0005;

	bool start; //true if this is a pickup point, false if this is a dropoff point
	TreeClusterNode *pickup; //NULL if the pickup has already been reached, in which limit should decrease for dropoff point (also NULL if point is pickup)

	TreeClusterNode *parent;
	vector<TreeClusterNode *> children;
	ShortestPath *shortestPath;
	
	long insert_uid; // unique id for the inserted pair
	
	double time; // time relative to parent
	double slackTime;
	double childSlackTime; //maximum of child slack times
	int layer; //layer in the tree
	
	int bestChild; //best child following this one in a route
	
	vertex *vert;
	
	bool chain;
	
	//constructs a root gentree node
	TreeClusterNode(ShortestPath *shortestPath, vertex *vert);
	
	//constructs a normal gentree node
	TreeClusterNode(TreeClusterNode *parent, vertex *vert, bool start, long insert_uid, bool pickupRemoved, double slackTime, double childSlackTime, bool chain, ShortestPath *shortestPath);
	
	//large constructor in case feasibility checking has already been done
	TreeClusterNode(TreeClusterNode *parent, vertex *vert, bool start, long insert_uid, double slackTime, double childSlackTime, double time, TreeClusterNode *pickup, bool chain, ShortestPath *shortestPath);
	
	~TreeClusterNode();
	
	//same as creating a new node with the given parameters, but doesn't allocate memory if it's infeasible
	// instead, if infeasible, this will return null
	static TreeClusterNode *constructCopy(TreeClusterNode *parent, TreeClusterNode *target, TreeClusterNode *branch, double detour, ShortestPath *shortestPath);
	
	//notify the gentree that the taxi has moved a certain distance
	// if pair is valid (-1 is definitely invalid), then it means that pickup point has been reached
	// so corresponding dropoff point will be updated with pickupRemoved = true
	// this way, we can decrease the limit for those dropoff points whose corresponding pickup has been reached
	void step(double distanceTraveled, long pair);
	
	//returns total path time, and stores the best child
	// uses simple DFS-based shortest path algorithm
	// should be called each time path may be affected by a node update (i.e., an insert; no need when we simple do a step)
	double bestTime();

	//returns total number of nodes in the tree
	int getNumberNodes();
	
	//updates the total child slack time
	map<int, double> updateChildSlackTime();
	double getMapSmallest(map<int, double> m);
	
	//returns the best child to follow this node in a route identified by bestTime()
	int getBestChild() { return bestChild; }
	
	// returns whether or not copy was successful; if not successful, parent maybe should delete
	bool copyNodes(vector<TreeClusterNode *> source, TreeClusterNode *branch, double detour);
	
	bool insertNodes(vector<TreeClusterNode *> doInsert, double depth);
	
	//calculates slack time
	double calculateSlackTime();
	
	//creates a copy of this node (identical parent)
	TreeClusterNode *clone();
	
	void print();
};

class TreeClusterTaxiPath : public TaxiPath {
private:
	TreeClusterNode *root;
	TreeClusterNode *rootTemp;
	bool flag; //1 if value was just called, 0 otherwise
	
	long nextPair;

public:
	TreeClusterTaxiPath(ShortestPath *shortestPath, vertex *curr); //currNode is index of the vertex that taxi is currently at or heading towards
	~TreeClusterTaxiPath();
	
	virtual void moved(double distance);
	virtual double value(vertex *curr, vertex *source, vertex *dest); //tests pushing (source, dest)
	virtual void cancel(); //cancels the value() call
	virtual void push(); //completes the pushing after value() is called
	virtual bool step(bool move); //next() becomes curr(); goes one vertex forward, when curr() is passed
	virtual queue<vertex *> next(); //returns next vertex that the taxi should go to after curr() is passed; -1 if no next
	virtual queue<vertex *> curr(vertex *curr); //returns the vertex that the taxi is currently at or heading towards
	
	virtual int size() { return 0; }
	virtual void printPoints(); //prints current point array
	virtual int getNumberNodes() { return root->getNumberNodes(); } //calls getNumberNodes on the root node

private:
	double euclidean(double ax, double ay, double bx, double by); //returns shortest path from source to dest
};

#endif
