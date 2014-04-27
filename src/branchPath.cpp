#include "includes.h"
#include "path.h"
#include "taxiPath.h"
#include "oneTaxiPath.h"
#include "vertex.h"
#include "branchPath.h"

BranchPathNode :: BranchPathNode(PathPoint *vert, int N, BranchPath* branch) {
	this->N = N;
	this->branch = branch;
	this->vertex = vert;
	level = -1;
	length = 0.0;
	bound = 0.0;
	
	parent = NULL;
	
	if(vertex->type == 0) {
		positions[vertex->pairIndex] = 0;
	} else {
		positions[vertex->pairIndex + 1] = -1;
	}
}
	
BranchPathNode :: BranchPathNode(BranchPathNode* parent, BranchPath* branch, PathPoint *vert) {
	this->parent = parent;
	parent->children.push_back(this);
	
	this->N = parent->N;
	this->branch = branch;
	this->vertex = vert;
	
	level = parent->level + 1;
	
	length = parent->length + branch->length(parent->vertex, vertex);
	bound = 0.0;
	
	//copy positions from parent
	map<int, double>::iterator it;
	for(it = parent->positions.begin(); it != parent->positions.end(); it++) {
		positions[(*it).first] = (*it).second;
	}
	
	//check feasibility
	feasible = true;
	
	if(vertex->type == 0) {
		//this is a pickup point
		if(length > vertex->remaining) {
			feasible = false;
		} else {
			positions[vertex->pairIndex] = length;
		}
	} else if(vertex->type == 2) {
		//this is a stand-alone dropoff point
		if(length > vertex->remaining) {
			feasible = false;
		} else {
			positions[vertex->pairIndex + 1] = -1;
		}
	} else {
		//this is a dropoff point
		map<int, double>::iterator it = positions.find(vertex->pairIndex);
		
		if(it != positions.end()) {
			double pairLength = it->second; //length at corresponding pickup point
	
			if(length - pairLength > vertex->remaining) {
				feasible = false;
			} else {
				positions[vertex->pairIndex + 1] = -1;
			}
		} else {
			feasible = false;
		}
	}
}

BranchPathNode :: ~BranchPathNode() {
	vector<BranchPathNode *>::iterator it;
	
	if(parent) {
		//delete self from parent's child list
		
		for(it = parent->children.begin(); it < parent->children.end(); it++) {
			if((*it) == this) {
				parent->children.erase(it);
				break;
			}
		}
	}
	
	//now delete the subtree, if any
	while(children.size() > 0) {
		delete children.back();
	}
}

bool BranchPathNode :: pathContains(PathPoint *test) {
	int pairIndex = test->pairIndex + 1;
	
	if(test->type == 0) {
		pairIndex--;
	}
	
	map<int, double>::iterator it = positions.find(pairIndex);
	return it != positions.end();
}

BranchPath :: BranchPath(TaxiPath *taxiPath, int M, vector<PathPoint *> &nPoints): PathHandler(taxiPath, M, nPoints) {
	minEdge = new double[M];
}

BranchPath :: ~BranchPath() {
	delete minEdge;
}

vector<PathPoint *> BranchPath :: findPath(vertex *curr) {
	//create a PathPoint so that we can treat the current
	// location as any other point
	PathPoint currentPoint;
	currentPoint.index = -1;
	currentPoint.vert = curr;
	currentPoint.type = -1;
	currentPoint.pairIndex = -5;
	currentPoint.remaining = -1;
	
	//find the minimum outgoing edge from each point
	//this assumes a complete graph (weight of each
	// edge is the shortest distance)
	init();
	
	//add the initial prefix path to the priority queue
	BranchPathNode *root = new BranchPathNode(&currentPoint, N, this);
	root->bound = bound(root);
	q.push( root );
	
	execute();
	
	if(root) {
		delete root; //this clears any remaining BranchPathNode elements
	}
	
	return bestList;
}

void BranchPath :: execute() {
	while (!q.empty( ) ) {
		//remove node with smallest bound from the queue
		BranchPathNode* temp = q.top();
		q.pop();

		if (temp->bound < minTour) {
			//go through every possible vertex from the next node
			for(int i = 0; i < N; i++) {
				if (!temp->pathContains(points[i]) ) {
					//if vertex is not already in the partial tour,
					// form a new partial tour that extends the tour in node
					// temp by appending the vertex
					
					BranchPathNode *u = new BranchPathNode(temp, this, points[i]);
					
					if(!u->feasible) {
						delete u;
						u = NULL;
						continue; //ignore if this is not feasible
					}
					
					if (u->level == N - 2) {
						//if the new partial tour is of length N -1, there is only
						//one possible path that can be formed -- form it now
						
						for(int j = 0; j < N; j++) {
							if (!u->pathContains(points[j]) ) {
								BranchPathNode *uChild = new BranchPathNode(u, this, points[j]);
								
								if (uChild->feasible && uChild->length < minTour) {
									//if this new path is the best so far, save it
									minTour = uChild->length;
									
									bestList.clear();
									BranchPathNode *curr = uChild;
									while(curr->parent) {
										bestList.push_back(curr->vertex);
										curr = curr->parent;
									}
									//don't push the root's vertex because that isn't in the original points
									
									//reverse bestList
									for(int i = 0; i < bestList.size() / 2; i++) {
										PathPoint *tmp = bestList[i];
										bestList[i] = bestList[bestList.size() - i - 1];
										bestList[bestList.size() - i - 1] = tmp;
									}
								}
								
								delete uChild;
								uChild = NULL;
								break;
							}
						}
						
						delete u;
						u = NULL;
					}
					else {
						//if the partial tour is "promising", add node to the priority queue
						u->bound = bound(u);
						
						if (u->bound < minTour)
							q.push(u);
						else {
							delete u;
							u = NULL;
						}
					}
				}
			}
		}
	}
}

double BranchPath :: bound(BranchPathNode* n) {
	//the bound is the current length of the path, plus
	// the minimum outgoing edge from each remaining vertex
	
	double bnd = n->length;
	for (int i = 0; i < N; i++) {
		if (!n->pathContains(points[i]))
			bnd += minEdge[i];
	}
	
	return bnd;
}

void BranchPath :: init( ) {
	//find and record the minimum outgoing edge from each vertex
	for(int i = 0; i < N; i++) {
		double cost = numeric_limits<double>::max();
		
		for(int j = 0; j < N; j++) {
			if(i == j) continue;
			
			//we now use Euclidean distance instead of road network distance for the bound
			// this is because Euclidean distance is faster to compute and still yields a suitable bound
			
			//double len = length(points[i], points[j]);
			double dx = points[i]->vert->x - points[j]->vert->x;
			double dy = points[i]->vert->y - points[j]->vert->y;
			double len = sqrt(dx * dx + dy * dy);
			if (len < cost) cost = len;
		}
		
		minEdge[i] = cost;
	}
}
