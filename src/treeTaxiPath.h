#ifndef TREETAXIPATH_H
#define TREETAXIPATH_H

class TaxiPath;
class TreeNode;
class ShortestPath;

class TreeNode {
public:
	bool start; //true if this is a pickup point, false if this is a dropoff point
	bool pickupRemoved; //true if the pickup has already been reached, in which limit should decrease for dropoff point

	TreeNode *parent;
	vector<TreeNode *> children;
	ShortestPath *shortestPath;
	
	long insert_uid; // unique id for the inserted pair
	
	double time; // time relative to parent
	double absoluteTime;
	double rootTime; // time relative to root
	double pairTime; // if start is false, this is the time between pickup location and this dropoff point
	
	double totalSlackTime; //slack time for this node, considering all of its children
	
	double limit; //tolerance threshold on rootTime or pairTime (i.e., the constraint)
	int bestChild; //best child following this one in a route
	
	vertex *vert;
	
	//constructs a root gentree node
	TreeNode(ShortestPath *shortestPath, vertex *vert);
	
	//constructs a normal gentree node
	TreeNode(TreeNode *parent, vertex *vert, bool start, long insert_uid, double limit, bool pickupRemoved, double totalSlackTime, ShortestPath *shortestPath);
	
	//large constructor in case feasibility checking has already been done
	TreeNode(TreeNode *parent, vertex *vert, bool start, long insert_uid, double limit, bool pickupRemoved, double time, double pairTime, double totalSlackTime, ShortestPath *shortestPath);
	
	~TreeNode();
	
	//same as creating a new node with the given parameters, but doesn't allocate memory if it's infeasible
	// instead, if infeasible, this will return null
	static TreeNode *safeConstructNode(TreeNode *parent, vertex *vert, bool start, long insert_uid, double limit, bool pickupRemoved, double totalSlackTime, ShortestPath *shortestPath);
	
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
	
	//returns the best child to follow this node in a route identified by bestTime()
	int getBestChild() { return bestChild; }
	
	// returns whether or not copy was successful; if not successful, parent maybe should delete
	bool copyNodes(vector<TreeNode *> *source, vector<TreeNode *> *doInsert);
	
	//returns whether or not this node is feasible
	bool feasible();
	
	//this is this node's slack time
	// field totalSlackTime considers children's slack times; this does not
	double slackTime();
	
	//calculates total slack time and stores in totalSlackTime field
	double calculateTotalSlackTime();
	
	//checks that newNode works with at least one of our children
	bool checkSlack(TreeNode *newNode);
	
	//creates a copy of this node with a new parent
	TreeNode *copy(TreeNode *new_parent);
	
	//creates a copy of this node with a new parent
	// or returns null if infeasible
	TreeNode *copySafe(TreeNode *new_parent);
	
	//creates a copy of this node (identical parent)
	TreeNode *clone();
	
	void print();
};

class TreeTaxiPath : public TaxiPath {
private:
	TreeNode *root;
	TreeNode *rootTemp;
	bool flag; //1 if value was just called, 0 otherwise
	
	long nextPair;

public:
	TreeTaxiPath(ShortestPath *shortestPath, vertex *curr); //currNode is index of the vertex that taxi is currently at or heading towards
	~TreeTaxiPath();
	
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
