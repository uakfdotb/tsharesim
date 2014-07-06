#include "includes.h"
#include "shortestPath.h"
#include "vertex.h"
#include "taxiPath.h"
#include "treeTaxiPath.h"

//constructs a root gentree node
TreeNode :: TreeNode::TreeNode(ShortestPath *shortestPath, vertex *vert) {
	parent = NULL;
	
	this->shortestPath = shortestPath;
	this->vert = vert;
	
	insert_uid = -1;
	start = true;

	time = 0;
	absoluteTime = 0;
	rootTime = 0;
	pickupRemoved = false;
	limit = -1;
	totalSlackTime = 0;
}

//constructs a normal gentree node
TreeNode :: TreeNode(TreeNode *parent, vertex *vert, bool start, long insert_uid, double limit, bool pickupRemoved, double totalSlackTime, ShortestPath *shortestPath) {
	this->parent = parent;
	this->vert = vert;
	this->insert_uid = insert_uid;
	this->start = start;
	this->limit = limit;
	this->pickupRemoved = pickupRemoved;
	this->totalSlackTime = totalSlackTime; //our actual slack time can only be less, and it's easier to recalculate later
	this->shortestPath = shortestPath;

	if(parent) {
		time = shortestPath->shortestDistance(parent->vert, vert);
		absoluteTime = parent->absoluteTime + time;
		rootTime = parent->rootTime + time;

		//find our pairTime if this is a dropoff who's pair is still a sub-root tree node
		if(!start && !pickupRemoved) {
			TreeNode *n = parent;
			while(n->insert_uid != insert_uid) {
				n = n->parent;
			}

			pairTime = absoluteTime - n->absoluteTime;
		} else {
			pairTime = 0;
		}
	} else {
		time = 0;
		absoluteTime = 0;
		rootTime = 0;
		pairTime = 0;
	}
}

//large constructor in case feasibility checking has already been done
TreeNode :: TreeNode(TreeNode *parent, vertex *vert, bool start, long insert_uid, double limit, bool pickupRemoved, double time, double pairTime, double totalSlackTime, ShortestPath *shortestPath) {
	this->parent = parent;
	this->vert = vert;
	this->insert_uid = insert_uid;
	this->start = start;
	this->limit = limit;
	this->pickupRemoved = pickupRemoved;
	this->totalSlackTime = totalSlackTime;

	this->time = time;
	absoluteTime = parent->absoluteTime + time;
	rootTime = parent->rootTime + time;
	this->shortestPath = shortestPath;

	this->pairTime = pairTime;
}

TreeNode :: ~TreeNode() {
	for(int i = 0; i < children.size(); i++) {
		delete children[i];
	}
}

//same as creating a new node with the given parameters, but doesn't allocate memory if it's infeasible
// instead, if infeasible, this will return null
TreeNode *TreeNode :: safeConstructNode(TreeNode *parent, vertex *vert, bool start, long insert_uid, double limit, bool pickupRemoved, double totalSlackTime, ShortestPath *shortestPath) {
	double time = shortestPath->shortestDistance(parent->vert, vert);

	if(start || pickupRemoved) {
		if(parent->rootTime + time > limit) {
			return NULL;
		} else {
			return new TreeNode(parent, vert, start, insert_uid, limit, pickupRemoved, time, 0, totalSlackTime, shortestPath);
		}
	} else {
		TreeNode *n = parent;
		while(n->insert_uid != insert_uid) {
			n = n->parent;
		}

		double pairTime = time + parent->absoluteTime - n->absoluteTime;
	
		if(pairTime > limit) {
			return NULL;
		} else {
			return new TreeNode(parent, vert, start, insert_uid, limit, pickupRemoved, time, pairTime, totalSlackTime, shortestPath);
		}
	}
}

//notify the gentree that the taxi has moved a certain distance
// if pair is valid (-1 is definitely invalid), then it means that pickup point has been reached
// so corresponding dropoff point will be updated with pickupRemoved = true
// this way, we can decrease the limit for those dropoff points whose corresponding pickup has been reached
void TreeNode :: step(double distanceTraveled, long pair) {
	if(limit > -1) {
		rootTime -= distanceTraveled;
	
		if(start || pickupRemoved) {
			limit -= distanceTraveled;
		} else if(pair == insert_uid) {
			pickupRemoved = true; //don't decrease pairtime for this one because the current location (root) is now the pair
		}
	}

	for(int i = 0; i < children.size(); i++) {
		children[i]->step(distanceTraveled, pair);
	}
}

//returns total path time, and stores the best child
// uses simple greedy algorithm
// should be called each time path may be affected by a node update (i.e., an insert; no need when we simple do a step)
double TreeNode :: bestTime() {
	if(children.size() > 0) {
		double bestTime = 10000000000; //todo: numeric limits<double>
		bestChild = -1;
		
		for(int i = 0; i < children.size(); i++) {
			double childBestTime = children[i]->bestTime() + children[i]->time;
			
			if(childBestTime < bestTime) {
				bestTime = childBestTime;
				bestChild = i;
			} else if(bestChild == -1) {
				cout << "errrrrr " << children[i]->bestTime() << "," << children[i]->time << "," << bestTime << endl;
			}
		}
		
		return bestTime;
	} else {
		bestChild = -1;
		return 0;
	}
}

int TreeNode :: getNumberNodes() {
	int numNodes = 1;

	for(int i = 0; i < children.size(); i++) {
		numNodes += children[i]->getNumberNodes();
	}
	
	return numNodes;
}

// returns whether or not copy was successful; if not successful, parent maybe should delete
bool TreeNode :: copyNodes(vector<TreeNode *> *sourcePtr, vector<TreeNode *> *doInsertPtr) {
	//recursively copy nodes from other trees if requested
	if(sourcePtr != NULL) {
		vector<TreeNode *> source = *sourcePtr;
		bool fail = false;

		for(int i = 0; i < source.size(); i++) {
			TreeNode *myCopy = source[i]->copy(this);
		
			if(myCopy->feasible()) {
				children.push_back(myCopy);
				
				if(!myCopy->copyNodes(&(source[i]->children), NULL)) {
					//our child has no feasible branches in it's subtree, delete
					children.pop_back();
					delete myCopy;
					fail = true;
				}
			} else {
				//child is not feasible, delete
				delete myCopy;
				fail = true;
			}
		}
		
		//check if we're infeasible
		//make sure to check fail flag because we might not have had to copy any nodes
		if(fail && children.empty()) {
			//we have no feasible branches in our subtree, fail
			return false;
		}
	}
	
	//inserted new nodes into this tree if requested
	if(doInsertPtr != NULL) {
		vector<TreeNode *> doInsert = *doInsertPtr;
		
		if(doInsert.size() > 0) {
			//take the first node to insert and put it directly underneath us
			TreeNode *insertCopy = doInsert[0]->copySafe(this);

			if(insertCopy) {
				bool fail = false;
			
				// first, copy other branches to our inserted copy
				// ignore slack time since this algorithm doesn't use it correctly
				if(children.empty() || checkSlack(insertCopy) || true) {
					if(!(insertCopy->copyNodes(&children, NULL))) {
						fail = true;
					}
				} else {
					fail = true;
				}
			
				if(!fail && doInsert.size() >= 2) {
					// give the insertCopy our doInsert, but with first element that was already copied by it removed
					vector<TreeNode *> doInsertClone;
					for(int i = 1; i < doInsert.size(); i++) {
						doInsertClone.push_back(doInsert[i]);
					}
				
					if(!insertCopy->copyNodes(NULL, &doInsertClone)) {
						fail = true;
					}
				}
			
				// give other children our insertCopy in full
				// note that this must be executed after insertCopy updates because insertCopy copies nodes from these children
				for(int it = 0; it < children.size(); it++) {
					TreeNode *child = children[it];
					
					if(!(child->copyNodes(NULL, doInsertPtr))) {
						// child failed, so delete from children
						delete child;
						children.erase(children.begin() + it);
						it--;
					}
				}
			
				if(!fail) {
					// lastly, insert insertCopy into our children so we don't have infinite loop
					children.push_back(insertCopy);
				} else if(children.empty()) {
					delete insertCopy;
					return false;
				} else {
					delete insertCopy;
				}
			} else {
				// if our inserted copy is not feasible, no other path will be feasible
				//  since we have to insert the same node to our children.. assuming shortest paths are shortest
				return false;
			}
		}
	}
	
	return true;
}

//returns whether or not this node is feasible
bool TreeNode :: feasible() {
	if(start || pickupRemoved) {
		return rootTime <= limit;
	} else {
		return pairTime <= limit;
	}
}

//this is this node's slack time
// field totalSlackTime considers children's slack times; this does not
double TreeNode :: slackTime() {
	if(start || pickupRemoved) {
		return limit - rootTime;
	} else {
		return 100000;
	}
}

//calculates total slack time and stores in totalSlackTime field
double TreeNode :: calculateTotalSlackTime() {
	if(children.empty()) {
		totalSlackTime = slackTime();
	} else {
		totalSlackTime = 0;
	
		for(int i = 0; i < children.size(); i++) {
			totalSlackTime = max(children[i]->calculateTotalSlackTime(), totalSlackTime);
		}
	
		//our slack time always has to be satisfied though
		totalSlackTime = min(totalSlackTime, slackTime());
	}
	
	return totalSlackTime;
}

//checks that newNode works with at least one of our children
bool TreeNode :: checkSlack(TreeNode *newNode) {
	for(int i = 0; i < children.size(); i++) {
		//compare slack time to the increased time due to detour
		// increased time is new tree distance minus old tree distance
		if(children[i]->totalSlackTime > shortestPath->shortestDistance(vert, newNode->vert) + shortestPath->shortestDistance(newNode->vert, children[i]->vert) - children[i]->time) {
			return true;
		}
	}
	
	return false;
}

//creates a copy of this node with a new parent
TreeNode *TreeNode :: copy(TreeNode *new_parent) {
	TreeNode *copy = new TreeNode(new_parent, vert, start, insert_uid, limit, pickupRemoved, totalSlackTime, shortestPath);
	return copy;
}

//creates a copy of this node with a new parent
// or returns null if infeasible
TreeNode *TreeNode :: copySafe(TreeNode *new_parent) {
	TreeNode *copy = safeConstructNode(new_parent, vert, start, insert_uid, limit, pickupRemoved, totalSlackTime, shortestPath);
	return copy;
}

//creates a copy of this node (identical parent)
TreeNode *TreeNode :: clone() {
	TreeNode *clone = new TreeNode(parent, vert, start, insert_uid, limit, pickupRemoved, totalSlackTime, shortestPath);
	
	for(int i = 0; i < children.size(); i++) {
		clone->children.push_back(children[i]->clone());
		clone->children[i]->parent = clone;
	}
	
	return clone;
}

void TreeNode :: print() {
	cout.precision(12);
	cout << vert->id << ": " << rootTime << "/" << (start ? "pickup" : "dropoff") << "/" << limit << " [";
	
	for(int i = 0; i < children.size(); i++) {
		cout << children[i]->vert->id << "  ";
	}
	
	cout << "]" << endl;
	
	for(int i = 0; i < children.size(); i++) {
		children[i]->print();
	}
}

TreeTaxiPath :: TreeTaxiPath(ShortestPath *shortestPath, vertex *curr): TaxiPath(shortestPath, curr) {
	root = new TreeNode(shortestPathC, curr);
	flag = false;
	
	nextPair = 0;
}

TreeTaxiPath :: ~TreeTaxiPath() {
	delete root;
}

void TreeTaxiPath :: moved(double distance) {
	root->step(distance, -1);
}

double TreeTaxiPath :: value(vertex *curr, vertex *source, vertex *dest) {
	flag = true;

	root->vert = curr;

	TreeNode *pick = new TreeNode(NULL, source, true, nextPair++, pickupConstraint, false, 0, shortestPathC);
	TreeNode *drop = new TreeNode(NULL, dest, false, pick->insert_uid, serviceConstraint * shortestPathC->shortestDistance(source, dest), false, 0, shortestPathC);

	vector<TreeNode *> doInsert;
	doInsert.push_back(pick);
	doInsert.push_back(drop);

	rootTemp = root->clone();
	bool copyResult = rootTemp->copyNodes(NULL, &doInsert);

	//clear our insertion vector elements
	delete pick;
	delete drop;

	//update node values to reflect insertion
	double time = rootTemp->bestTime();
	rootTemp->calculateTotalSlackTime();

	if(time != 0 && copyResult) return time;
	else return -1;
}

void TreeTaxiPath :: cancel() {
	if(flag) {
		delete rootTemp;
		
		flag = false;
	}
}

void TreeTaxiPath :: push() {
	delete root;
	root = rootTemp;
	rootTemp = NULL;

	flag = false;
}

bool TreeTaxiPath :: step(bool move) {
	int rootBestChild = root->getBestChild();
	
	for(int i = root->children.size() - 1; i >= 0; i--) {
		if(i != rootBestChild) {
			delete root->children[i];
			root->children.erase(root->children.begin() + i);
		}
	}
	
	TreeNode *removedChild = root->children.back(); //should be only remaining node
	root->children.pop_back();
	
	for(int i = 0; i < removedChild->children.size(); i++) {
	  root->children.push_back(removedChild->children[i]);
	  root->children[i]->parent = root;
	}
	
	root->bestChild = removedChild->getBestChild();
	root->step(0, removedChild->insert_uid); //todo: if move is true we should try to move from current root position to removed node
	
	bool droppedPassenger = !(removedChild->start);
	removedChild->children.clear();
	delete removedChild;
	
	return droppedPassenger; //return true if we dropped a passenger
}

queue<vertex *> TreeTaxiPath :: next() {
	queue<vertex *> ret;
	return ret; //todo
 }

queue<vertex *> TreeTaxiPath :: curr(vertex *curr) {
	root->vert = curr;
	
 	if(root->children.size() > 0) {
 		int rootBestChild = root->getBestChild();
		return shortestPathC->shortestPath(root->vert, root->children[rootBestChild]->vert);
 	} else {
		queue<vertex *> emptyqueue;
		return emptyqueue;
 	}
 }

void TreeTaxiPath :: printPoints() {
	root->print();
	
	cout << endl << "**best**" << endl;
	
	TreeNode *curr = root;
	
	while(curr->children.size() > 0) {
		curr = curr->children[curr->getBestChild()];
		cout << curr->vert->id << endl;
	}
}

double TreeTaxiPath :: euclidean(double ax, double ay, double bx, double by) {
	double d1 = ax - bx;
	double d2 = ay - by;
	return sqrt(d1 * d1 + d2 * d2);
}
