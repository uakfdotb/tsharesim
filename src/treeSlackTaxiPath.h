#ifndef TREESLACKTAXPATH_H
#define TREESLACKTAXIPATH_H

class TaxiPath;
class TreeSlackNode;
class ShortestPath;

class TreeSlackNode {
public:
	bool start; //true if this is a pickup point, false if this is a dropoff point
	TreeSlackNode *pickup; //NULL if the pickup has already been reached, in which limit should decrease for dropoff point (also NULL if point is pickup)

	TreeSlackNode *parent;
	vector<TreeSlackNode *> children;
	ShortestPath *shortestPath;
	
	long insert_uid; // unique id for the inserted pair
	
	double time; // time from parent to here
	
	//our own slack time: how much detour we can tolerate until
	// we violate our pickup or service constraint
	double slackTime;
	
	//this slack time is the maximum slack time of our branches
	// we only need one branch to work for ourself to be feasible
	double childSlackTime;
	
	int layer; //layer in the tree
	int bestChild; //best child following this one in a route
	vertex *vert; //vertex in graph for this node
	
	//constructs a root gentree node
	TreeSlackNode(ShortestPath* shortestPath, vertex *vert);
	
	//constructs a normal gentree node
	TreeSlackNode(TreeSlackNode *parent, vertex *vert, bool start, long insert_uid, bool pickupRemoved, double slackTime, double childSlackTime, ShortestPath *shortestPath);
	
	//large constructor in case feasibility checking has already been done
	TreeSlackNode(TreeSlackNode *parent, vertex *vert, bool start, long insert_uid, double slackTime, double childSlackTime, double time, TreeSlackNode *pickup, ShortestPath *shortestPath);
	
	~TreeSlackNode();
	
	//same as creating a new node with the given parameters, but doesn't allocate memory if it's infeasible
	// instead, if infeasible, this will return null
	static TreeSlackNode *constructCopy(TreeSlackNode *parent, TreeSlackNode *target, TreeSlackNode *branch, double detour, ShortestPath *shortestPath);
	
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
	bool copyNodes(vector<TreeSlackNode *> source, TreeSlackNode *branch, double detour);
	
	bool insertNodes(vector<TreeSlackNode *> doInsert, double depth);
	
	//calculates slack time
	double calculateSlackTime();
	
	//creates a copy of this node (identical parent)
	TreeSlackNode *clone();
	
	void print();
};

class TreeSlackTaxiPath : public TaxiPath {
private:
	TreeSlackNode *root;
	TreeSlackNode *rootTemp;
	bool flag; //1 if value was just called, 0 otherwise
	
	long nextPair;

public:
	TreeSlackTaxiPath(ShortestPath *shortestPath, vertex *curr); //currNode is index of the vertex that taxi is currently at or heading towards
	~TreeSlackTaxiPath();
	
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
