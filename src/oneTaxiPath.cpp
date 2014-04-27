#include "includes.h"
#include "path.h"
#include "shortestPath.h"
#include "vertex.h"
#include "taxiPath.h"
#include "branchPath.h"
#include "brutePath.h"
#include "lpPath.h"
#include "oneTaxiPath.h"

OneTaxiPath :: OneTaxiPath(ShortestPath *shortestPath, vertex *curr): TaxiPath(shortestPath, curr) {
	nextPair = 0;
	flag = false;
	
	if(ALGORITHM == 0) {
		pathb = new BranchPath(this, MAX_PASSENGERS * 2, points);
	} else if(ALGORITHM == 1) {
		pathb = new BrutePath(this, MAX_PASSENGERS * 2, points);
	} else if(ALGORITHM == 2) {
		pathb = new LPPath(this, MAX_PASSENGERS * 2, points);
	}
}

OneTaxiPath :: ~OneTaxiPath() {
	points.erase(points.begin(), points.end());
}

void OneTaxiPath :: moved(double distance) {
	//for each pickup point, and also each dropoff point where
	// the passenger was already picked up, we subtract the
	// remaining constraint distance by the distance moved
	for(int i = 0; i < points.size(); i++) {
		if(points[i]->type == 0 || points[i]->type == 2) {
			points[i]->remaining -= distance;
		}
	}
}

double OneTaxiPath :: value(vertex *curr, vertex *source, vertex *dest) {
	this->curr_vert = curr;
	
	//construct two new path points for the pickup
	// and dropoff points of this request
	PathPoint *p1 = new PathPoint;
	p1->index = points.size();
	p1->vert = source;
	p1->type = 0;
	p1->pairIndex = nextPair;
	p1->remaining = PICKUP_CONST;
	points.push_back(p1);
	
	PathPoint *p2 = new PathPoint;
	p2->index = points.size();
	p2->vert = dest;
	p2->type = 1;
	p2->pairIndex = nextPair;
	p2->remaining = shortestPath(source, dest) * SERVICE_CONST;
	points.push_back(p2);
	
	//increment nextPair by two instead of one because some
	// algorithms depend on it o_O
	nextPair += 2;
	
	bool ppLog = false;
	if(LOG_PATHPOINTS) {
		//we want to log the cases where the number of
		// passengers = MAX_PASSENGERS so that we can feed
		// these in directly later to specific algorithms
		// through TEST_MODE=2
		int numPassengers = 0;
		for(int i = 0; i < points.size(); i++) {
			if(points[i]->type == 1) numPassengers++;
		}
		
		if(numPassengers >= MAX_PASSENGERS) {
			ppLog = true;
			
			cout << "pplog:";
			cout << " " << curr->id;
		
			for(int i = 0; i < points.size(); i++) {
				cout << " " << points[i]->vert->id << "/" << points[i]->type << "/" << points[i]->pairIndex << "/" << points[i]->remaining;
			}
		
			cout << endl;
		}
	}
	
	pathb->reset();
	temp = pathb->findPath(curr_vert);
	flag = true;
	
	double result = pathb->getMinTour();
	if(result == numeric_limits<double>::max()) result = -1;
	
	if(ppLog) {
		cout << "pprlog: " << result << endl;
	}
	
	return result;
}

void OneTaxiPath :: cancel() {
	if(flag) {
		//delete the temporary points
		PathPoint *p = points.back();
		points.pop_back();
		delete p;
		
		p = points.back();
		points.pop_back();
		delete p;
		
		temp.clear();
		flag = false;
	}
}

void OneTaxiPath :: push() {
	points = temp;
	
	//re-index the points
	for(int i = 0; i < points.size(); i++) {
		points[i]->index = i;
	}
	
	temp.clear();
	flag = false;
}

bool OneTaxiPath :: step(bool move) {
	int pairIndex = points[0]->pairIndex;
	
	if(move) {
		//deprecated: let taxitree.h handle movement for now
	}
	
	bool droppedPassenger = points[0]->type == 2;
	
	//delete the first point and re-index
	//also, if it was a pickup point, update
	// the corresponding dropoff point
	delete points[0];
	points.erase(points.begin());
	
	for(int i = 0; i < points.size(); i++) {
		points[i]->index = i;
		
		if(points[i]->pairIndex == pairIndex && points[i]->type == 1) {
			points[i]->type = 2;
		}
	}
	
	return droppedPassenger;
}

queue<vertex *> OneTaxiPath :: next() {
	if(points.size() >= 2) {
		return shortestPathC->shortestPath(points[0]->vert, points[1]->vert);
	} else {
		queue<vertex *> emptyqueue;
		return emptyqueue;
	}
}

queue<vertex *> OneTaxiPath :: curr(vertex *curr) {
	curr_vert = curr;
	
	if(points.size() >= 1) {
		return shortestPathC->shortestPath(curr_vert, points[0]->vert);
	} else {
		queue<vertex *> emptyqueue;
		return emptyqueue;
	}
}

void OneTaxiPath :: printPoints() {
	for(int i = 0; i < points.size(); i++) {
		cout << points[i]->vert->id << ":" << points[i]->remaining << endl;
	}
	
	cout << endl;
}

void OneTaxiPath :: test(vertex *curr, vector<PathPoint *> nPoints) {
	this->curr_vert = curr;
	points.clear();
	
	for(int i = 0; i < nPoints.size(); i++) {
		points.push_back(nPoints[i]);
		points[i]->index = i;
	}
	
	pathb->reset();
	temp = pathb->findPath(curr_vert);
	flag = true;
	
	if(pathb->getMinTour() == numeric_limits<double>::max()) cout << "no path found" << endl;
	else {
		double minTour = pathb->getMinTour();
		cout << "minimum tour: " << minTour << endl;
	}
}
