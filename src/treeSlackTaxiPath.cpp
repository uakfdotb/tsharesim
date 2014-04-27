#include "includes.h"
#include "shortestPath.h"
#include "vertex.h"
#include "taxiPath.h"
#include "treeSlackTaxiPath.h"

//constructs a root gentree node
TreeSlackNode :: TreeSlackNode::TreeSlackNode(ShortestPath *shortestPath, vertex *vert) {
	parent = NULL;
	
	this->shortestPath = shortestPath;
	this->vert = vert;
	
	insert_uid = -1;
	start = true;

	layer = 0;
	pickup = NULL;
	time = 0;
	slackTime = 10000;
	childSlackTime = slackTime;
}

//constructs a normal gentree node
TreeSlackNode :: TreeSlackNode(TreeSlackNode *parent, vertex *vert, bool start, long insert_uid, bool pickupRemoved, double slackTime, double childSlackTime, ShortestPath *shortestPath) {
	this->parent = parent;
	this->vert = vert;
	this->insert_uid = insert_uid;
	this->start = start;
	this->slackTime = slackTime;
	this->childSlackTime = childSlackTime;
	this->shortestPath = shortestPath;

	if(parent) {
		time = shortestPath->shortestDistance(parent->vert, this->vert);
		layer = parent->layer + 1;

		//find our pickup if this is a dropoff who's pair is still a sub-root tree node
		if(!start && !pickupRemoved) {
			pickup = parent;
			while(pickup->insert_uid != insert_uid) {
				pickup = pickup->parent;
			}
		} else {
			pickup = NULL;
		}
	} else {
		time = 0;
		layer = 0;
		pickup = NULL;
	}
}

//large constructor in case feasibility checking has already been done
TreeSlackNode :: TreeSlackNode(TreeSlackNode *parent, vertex *vert, bool start, long insert_uid, double slackTime, double childSlackTime, double time, TreeSlackNode *pickup, ShortestPath *shortestPath) {
	this->parent = parent;
	this->vert = vert;
	this->insert_uid = insert_uid;
	this->start = start;
	this->slackTime = slackTime;
	this->childSlackTime = childSlackTime;
	this->time = time;
	this->shortestPath = shortestPath;
	
	if(parent) {
		layer = parent->layer + 1;
	} else {
		layer = 0;
	}
	
	this->pickup = pickup;
}

TreeSlackNode :: ~TreeSlackNode() {
	for(int i = 0; i < children.size(); i++) {
		delete children[i];
	}
}

//same as creating a new node with the given parameters, but doesn't allocate memory if it's infeasible
// instead, if infeasible, this will return null
TreeSlackNode *TreeSlackNode :: constructCopy(TreeSlackNode *parent, TreeSlackNode *target, TreeSlackNode *branch, double detour, ShortestPath *shortestPath) {
	double time = shortestPath->shortestDistance(parent->vert, target->vert);
	
	//we only consider target's child slack time if the node being
	// inserted (branch) is directly above the target node
	//otherwise, we only look at target's slack time, because the
	// children may or may not be affected
	bool direct = branch == parent;

	if(target->start || target->pickup == NULL) {
		if((direct && detour > target->calculateSlackTime()) || detour > target->slackTime) {
			return NULL;
		} else {
			double newChildSlackTime = target->childSlackTime;
			if(direct) newChildSlackTime -= detour;
			
			return new TreeSlackNode(parent, target->vert, target->start, target->insert_uid, target->slackTime - detour, newChildSlackTime, time, NULL, shortestPath);
		}
	} else {
		//check whether the given branch node occurs first, or our pickup node
		// if branch occurs first, then that means that we are extending the distance
		//  from our pickup node and have to subtract from slack times
		// otherwise, though, the service distance is not increasing
		
		bool branchFirst = false;
		TreeSlackNode *n = parent;
		
		while(n->insert_uid != target->insert_uid) {
			if(n == branch) branchFirst = true;
			n = n->parent;
		}
		
		double slackTime = target->slackTime;
		double childSlackTime = target->childSlackTime;
		double totalSlackTime = target->calculateSlackTime();
		
		if(branchFirst) {
			slackTime -= detour;
			
			if(direct) {
				childSlackTime -= detour;
				totalSlackTime -= detour;
			}
		}

		if((direct && totalSlackTime < 0) || slackTime < 0)
			return NULL;
		else
			return new TreeSlackNode(parent, target->vert, target->start, target->insert_uid, slackTime, childSlackTime, time, n, shortestPath);
	}
}

//notify the gentree that the taxi has moved a certain distance
// if pair is valid (-1 is definitely invalid), then it means that pickup point has been reached
// so corresponding dropoff point will be updated with pickupRemoved = true
// this way, we can decrease the limit for those dropoff points whose corresponding pickup has been reached
void TreeSlackNode :: step(double distanceTraveled, long pair) {
	if(!start && pickup != NULL && pair == insert_uid) {
		pickup = NULL;
	}

	for(int i = 0; i < children.size(); i++) {
		children[i]->step(distanceTraveled, pair);
	}
}

//returns total path time, and stores the best child
// uses simple greedy algorithm
// should be called each time path may be affected by a node update (i.e., an insert; no need when we simple do a step)
double TreeSlackNode :: bestTime() {
	if(children.size() > 0) {
		double bestTime = 10000000000; //todo: numeric limits<double>
		
		for(int i = 0; i < children.size(); i++) {
			double childBestTime = children[i]->bestTime() + children[i]->time;
			
			if(childBestTime < bestTime) {
				bestTime = childBestTime;
				bestChild = i;
			}
		}
		
		return bestTime;
	} else {
		return 0;
	}
}

int TreeSlackNode :: getNumberNodes() {
	int numNodes = 1;

	for(int i = 0; i < children.size(); i++) {
		numNodes += children[i]->getNumberNodes();
	}
	
	return numNodes;
}

//updates child slack times in the subtree to accurate values
//returns a map of the slack times
// -2: always should be considered by parents for slack time
// insert id: for dropoff nodes, only consider until we hit corresponding pickup node
map<int, double> TreeSlackNode :: updateChildSlackTime() {
	if(children.size() > 0) {
		//reset our child-based slack time
		this->childSlackTime = -1;
		
		//this will contain the slack times of the passenger
		// pickup/dropoff locations along the child branch
		// that has the maximum overall slack time
		map<int, double> childSlacks;
	
		for(int i = 0; i < children.size(); i++) {
			//get child slack time and update childSlacks
			//we want to take childSlacks for the child with the greatest
			// slack time, because we can use any branch
			//the slack time of a single path in the tree, though, is the
			// minimum slack time of the nodes along the path
			map<int, double> c_slacks = children[i]->updateChildSlackTime();
			
			if(this->start) { //pickup
				//erase the slack time of our corresponding dropoff node
				map<int, double>::iterator it = c_slacks.find(insert_uid);
				if(it != c_slacks.end()) c_slacks.erase(it);
			}
			
			double slack = getMapSmallest(c_slacks);
			
			//if this slack time is more lenient than our current one, use it
			if(this->childSlackTime == -1 || slack > this->childSlackTime) {
				this->childSlackTime = slack;
				childSlacks = c_slacks;
			}
		}
		
		//for our parents, set our own slack time in the map
		//we will always be the first occurence of our insert_uid because 
		if((this->start || this->pickup == NULL) && slackTime < childSlacks[-2]) {
			//our node seems to be the minimum slack time along the
			// branch, so set in the map
			childSlacks[-2] = slackTime;
		} else {
			//if we're a dropoff node, set our slack time in the map
			//this will be done if our slack time is greater than the
			// current minimum, which we don't really want, but I don't
			// want to change the code in case I break something and it
			// shouldn't cause any problems, since it's a more lenient
			// value and won't be considered
			childSlacks[insert_uid] = slackTime;
		}
		
		return childSlacks;
	} else {
		//we have no children, so just use our slack time
		this->childSlackTime = slackTime;
		
		map<int, double> t_slacks;
		
		if(this->start || this->pickup == NULL) {
			//if we're a pickup or stand-alone dropoff node, then all
			// nodes above us should consider our slack time
			t_slacks[-2] = slackTime;
		} else {
			//if we're a dropoff node, our slack time should be ignored
			// once we pass the corresponding pickup node
			t_slacks[-2] = 10000;
			t_slacks[insert_uid] = slackTime;
		}
		
		return t_slacks;
	}
}

//returns the smallest value in a map, or a large
// number if the map is empty
double TreeSlackNode :: getMapSmallest(map<int, double> m) {
	double currMin = -1;

	for(map<int, double>::iterator it = m.begin(); it != m.end(); it++) {
		if(currMin == -1 || (*it).second < currMin) {
			currMin = (*it).second;
		}
	}
	
	if(currMin != -1)
		return currMin;
	else return 10000;
}

// copies source nodes to be underneath this node
// branch is the location where detour is being introduced, so that we can see if slack times should be decreased
// returns whether or not copy was successful; if not successful, parent maybe should delete
bool TreeSlackNode :: copyNodes(vector<TreeSlackNode *> source, TreeSlackNode *branch, double detour) {
	bool fail = false;
	
	for(int i = 0; i < source.size(); i++) {
		TreeSlackNode *myCopy = constructCopy(this, source[i], branch, detour, shortestPath);
		
		if(myCopy) {
			children.push_back(myCopy);
		
			if(!myCopy->copyNodes(source[i]->children, branch, detour)) {
				//something's infeasible in our child's subtree, so delete child
				children.pop_back();
				fail = true;
			}
		} else {
			//child is infeasible, so fail
			fail = true;
		}
	}
	
	//determine whether or not copy was successful
	//we check the fail flag because we may not have had to copy any nodes at all
	if(fail && children.empty()) {
		return false;
	} else {
		return true;
	}
}

//attempt to insert new nodes into the tree
//doInsert will always be either size 2 containing a pickup and
// then a dropoff node, or size 1 and containing a dropoff
//the meaning of depth varies depending on if we're inserting pickup or dropoff
// pickup: distance from the root node
// dropoff: extra distance from the pickup node
bool TreeSlackNode :: insertNodes(vector<TreeSlackNode *> doInsert, double depth) {
	if(doInsert.size() > 0) {
		//try to create the first node to be inserted
		//for pickup nodes, the slack time is initialized to the pickup
		// constraint, so the detour introduced will be the distance from
		// the new node to the root node, or:
		//    depth + distance(here, new node)
		//for dropoff nodes, the slack time is how much extra distance to
		// the pickup node is tolerable, and the depth is set such that the
		// formula above will still work
		TreeSlackNode *insertCopy = constructCopy(this, doInsert[0], this, depth + shortestPath->shortestDistance(this->vert, doInsert[0]->vert), shortestPath);
		
		if(insertCopy) {
			bool fail = false;
			
			// first, copy other branches to our inserted copy
			vector<TreeSlackNode *> copySource;
			for(int i = 0; i < children.size(); i++) {
				copySource.clear();
				copySource.push_back(children[i]);
				
				//the child node will experience a detour between here and the
				// inserted node, and that detour will be equal to the time it
				// takes to go to the new node and back to the child minus the
				// time it takes to go directly from here to the child
				if(!(insertCopy->copyNodes(copySource, insertCopy, insertCopy->time + shortestPath->shortestDistance(insertCopy->vert, children[i]->vert) - children[i]->time))) {
					fail = true;
					break;
				}
			}
			
			if(!fail && doInsert.size() >= 2) {
				//if we haven't failed, we need to give to give the new node
				// the rest of the nodes to be inserted
				insertCopy->updateChildSlackTime();
				
				// give the insertCopy our doInsert, but with first element that was already copied by it removed
				vector<TreeSlackNode *> doInsertClone;
				for(int i = 1; i < doInsert.size(); i++) {
					doInsertClone.push_back(doInsert[i]);
				}
				
				// restart depth counter since this is the corresponding dropoff node
				// it's slack time is remaining service time, and so will only be affected by distance to the pickup node
				// we include the negative number because we only want to include EXTRA time from going directly in the detour
				if(!insertCopy->insertNodes(doInsertClone, -doInsertClone[0]->time)) {
					fail = true;
				}
			}
		
			// give other children our insertCopy in full
			// note that this must be executed after insertCopy updates because insertCopy copies nodes from these children
			for(int it = 0; it < children.size(); it++) {
				TreeSlackNode *child = children[it];
				
				if(!(child->insertNodes(doInsert, depth + child->time))) {
					// child failed, so delete from children
					delete child;
					children.erase(children.begin() + it);
					it--;
				}
			}
		
			if(!fail) {
				// lastly, insert insertCopy into our children (do this last so we don't push insertCopy into itself)
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
	
	return true;
}

//this node's TOTAL slack time
//this depends on both the slack time of just this node
// and slack time of the most lenient branch
double TreeSlackNode :: calculateSlackTime() {
	if(slackTime > childSlackTime)
		return childSlackTime;
	else
		return slackTime;
}

//creates a copy of this node (identical parent)
TreeSlackNode *TreeSlackNode :: clone() {
	TreeSlackNode *clone = new TreeSlackNode(parent, vert, start, insert_uid, pickup == NULL, slackTime, childSlackTime, shortestPath);
	clone->childSlackTime = childSlackTime;
	
	for(int i = 0; i < children.size(); i++) {
		clone->children.push_back(children[i]->clone());
		clone->children[i]->parent = clone;
	}
	
	return clone;
}

void TreeSlackNode :: print() {
	cout.precision(12);
	cout << vert->id << ": " << time << "/" << (start ? "pickup" : "dropoff") << (pickup == NULL ? "" : "-c") << "/" << slackTime << "/" << childSlackTime << " [";
	
	for(int i = 0; i < children.size(); i++) {
		cout << children[i]->vert->id << "  ";
	}
	
	cout << "]" << endl;
	
	for(int i = 0; i < children.size(); i++) {
		children[i]->print();
	}
}

TreeSlackTaxiPath :: TreeSlackTaxiPath(ShortestPath *shortestPath, vertex *curr): TaxiPath(shortestPath, curr) {
	root = new TreeSlackNode(shortestPathC, curr);
	flag = false;
	
	nextPair = 0;
}

TreeSlackTaxiPath :: ~TreeSlackTaxiPath() {
	delete root;
}

void TreeSlackTaxiPath :: moved(double distance) {
	root->step(distance, -1);
}

double TreeSlackTaxiPath :: value(vertex *curr, vertex *source, vertex *dest) {
	flag = true;
	
	root->vert = curr;
	
	//create tree nodes for the pickup and dropoff point
	//for pickup, slack time starts at pickup constraint
	//for dropoff, slack time is extra distance tolerable from pickup node
	TreeSlackNode *pick = new TreeSlackNode(NULL, source, true, nextPair++, PICKUP_CONST, 10000, 0, NULL, shortestPathC);
	double shortestDistance = shortestPathC->shortestDistance(source, dest);
	TreeSlackNode *drop = new TreeSlackNode(NULL, dest, false, pick->insert_uid, (SERVICE_CONST - 1) * shortestDistance, 10000, shortestDistance, pick, shortestPathC);
	
	vector<TreeSlackNode *> doInsert;
	doInsert.push_back(pick);
	doInsert.push_back(drop);
	
	//create a temporary copy of the current tree and insert the nodes
	rootTemp = root->clone();
	bool copyResult = rootTemp->insertNodes(doInsert, 0);
	
	//clear our insertion vector elements
	delete pick;
	delete drop;
	
	//update node values to reflect insertion
	double time = rootTemp->bestTime();
	
	if(time != 0 && copyResult) return time;
	else return -1;
}

void TreeSlackTaxiPath :: cancel() {
	//make sure to check flag because sometimes cancel will
	// be called without a corresponding call to value because
	// of shortest distance filtering on the pickup constraint
	if(flag) {
		delete rootTemp;
		
		flag = false;
	}
}

void TreeSlackTaxiPath :: push() {
	//overwrite our current tree with the temporary copy
	
	if(flag) {
		delete root;
		root = rootTemp;
	
		//update the slack time
		//we don't have to maintain accurate slack times during the insert
		// nodes operation the insertion of a node only affects the slack
		// time of the subtree
		//additionally, the only node that will be further inserted into
		// the subtree is the dropoff node, and that case is handled correctly
		root->updateChildSlackTime();
		
		rootTemp = NULL;

		flag = false;
	}
}

bool TreeSlackTaxiPath :: step(bool move) {
	int rootBestChild = root->getBestChild();
	
	//delete every child except the best option from the root node
	for(int i = root->children.size() - 1; i >= 0; i--) {
		if(i != rootBestChild) {
			delete root->children[i];
			root->children.erase(root->children.begin() + i);
		}
	}
	
	//take the best child and move it's children to be
	// children of the root node
	TreeSlackNode *removedChild = root->children.back(); //should be only remaining node
	root->children.pop_back();
	
	for(int i = 0; i < removedChild->children.size(); i++) {
	  root->children.push_back(removedChild->children[i]);
	  root->children[i]->parent = root;
	}
	
	//update the best child index
	root->bestChild = removedChild->getBestChild();
	
	//notify tree that we've moved
	root->step(0, removedChild->insert_uid); //todo: if move is true we should try to move from current root position to removed node
	
	//update slack times
	root->updateChildSlackTime();
	
	bool droppedPassenger = !(removedChild->start);
	removedChild->children.clear();
	delete removedChild;
	
	return droppedPassenger; //return true if we dropped a passenger
}

queue<vertex *> TreeSlackTaxiPath :: next() {
	queue<vertex *> ret;
	return ret; //todo
}

queue<vertex *> TreeSlackTaxiPath :: curr(vertex *curr) {
	root->vert = curr;

	if(root->children.size() > 0) {
		int rootBestChild = root->getBestChild();
		return shortestPathC->shortestPath(root->vert, root->children[rootBestChild]->vert);
	} else {
		queue<vertex *> emptyqueue;
		return emptyqueue;
	}
}

void TreeSlackTaxiPath :: printPoints() {
	root->print();
	
	cout << endl << "**best**" << endl;
	
	TreeSlackNode *curr = root;
	
	while(curr->children.size() > 0) {
		curr = curr->children[curr->getBestChild()];
		cout << curr->vert->id << endl;
	}
}

double TreeSlackTaxiPath :: euclidean(double ax, double ay, double bx, double by) {
	double d1 = ax - bx;
	double d2 = ay - by;
	return sqrt(d1 * d1 + d2 * d2);
}
