#include "includes.h"
#include "taxiPath.h"
#include "location.h"
#include "vertex.h"
#include "shortestPath.h"
#include "oneTaxiPath.h"
#include "treeTaxiPath.h"
#include "treeSlackTaxiPath.h"
#include "treeClusterTaxiPath.h"
#include "taxitree.h"

void taxiTree :: setInitialPosition(ShortestPath *shortestPath, vertex *vertex) {
	if(ALGORITHM == 0 || ALGORITHM == 1 || ALGORITHM == 2) {
		taxiPath = new OneTaxiPath(shortestPath, vertex);
	} else if(ALGORITHM == 3) {
		taxiPath = new TreeTaxiPath(shortestPath, vertex);
	} else if(ALGORITHM == 4) {
		taxiPath = new TreeSlackTaxiPath(shortestPath, vertex);
	} else if(ALGORITHM == 5) {
		taxiPath = new TreeClusterTaxiPath(shortestPath, vertex);
	} else {
		cout << "Invalid algorithm specified!" << endl;
	}
	
	previousVertex = vertex;
	nextVertex = vertex;
	
	current.longitude = vertex->x;
	current.latitude = vertex->y;
	
	distToDest = -1;
	passengers = 0;
}

double taxiTree :: distance(double x1, double y1, double x2, double y2) {
	double d1 = x2 - x1;
	double d2 = y2 - y1;
	return sqrt(d1 * d1 + d2 * d2);
}

void taxiTree :: setDestination(vertex *newVertex) {
	previousVertex = nextVertex;
	nextVertex = newVertex;
	
	//let the taxiPath know that we're going to be moving from previousVertex to nextVertex
	//this is necessary because value() calls will assume that we have already reached nextVertex
	//this is not a problem because:
	// 1. once we are heading towards nextVertex, we do not turn back along the edge (because this may not be possible)
	// 2. before the first call to setDestination, previousVertex = nextVertex
	taxiPath->moved(distance(previousVertex->x, previousVertex->y, nextVertex->x, nextVertex->y));
	
	distToDest = distance(previousVertex->x, previousVertex->y, nextVertex->x, nextVertex->y);
	randDest = false;
}

void taxiTree :: randomDestination() {
	int randEdge = rand() % nextVertex->neighbors.size();
	setDestination(nextVertex->neighbors[randEdge]);
	randDest = true;
}

void taxiTree :: updateLocation() {
	if(distToDest != -1) { //if we have a destination
		distToDest -= D; //D is the constant taxi speed
		double distanceRemaining = distToDest; //copy distToDest so it's not deleted when we call setDestination
		
		//check if we passed the destination
		if(distanceRemaining <= 0) {
			if(!randDest) {
				//now we update destination
				
				if(!targetPath.empty()) {
					setDestination(targetPath.front());
					targetPath.pop();
				} else {
					bool result = taxiPath->step(false);
			
					if(result) { //we dropped a passenger
						passengers--;
					}
					
					//get path to next pickup/dropoff point
					targetPath = taxiPath->curr(nextVertex);
					
					if(!targetPath.empty()) { //new destination available?
						setDestination(targetPath.front());
						targetPath.pop();
					} else {
						//we have no more passengers to pickup or dropoff
						//we can start cruising now!
						randomDestination();
					}
				}
			} else {
				//we arrived at a random destination, so continue to another random one
				randomDestination();
			}
			
			distToDest += distanceRemaining; //subtract any distance remaining for accuracy, by adding the negative value
		} else {
			//recalculate current location
			double totalDistance = distance(previousVertex->x, previousVertex->y, nextVertex->x, nextVertex->y);
			double remainingPart = (double) distanceRemaining / totalDistance;
		
			if(previousVertex->x > nextVertex->x) {
				current.longitude = (previousVertex->x - nextVertex->x) * remainingPart + nextVertex->x;
			} else {
				current.longitude = nextVertex->x - (nextVertex->x - previousVertex->x) * remainingPart;
			}
		
			if(previousVertex->y > nextVertex->y) {
				current.latitude = (previousVertex->y - nextVertex->y) * remainingPart + previousVertex->y;
			} else {
				current.latitude = nextVertex->y - (nextVertex->y - previousVertex->y) * remainingPart;
			}
		}
	} else {
		//we don't have a destination; assign a random one
		randomDestination();
    }
}

double taxiTree :: value(vertex *source, vertex *dest) {
	if(distance(current.longitude, current.latitude, source->x, source->y) > PICKUP_CONST) {
		//special return code to denote that taxiPath->value was not called, so ->cancel should not be either
		//at the time of this comment, it calls cancel anyway though...
		return -2;
	}
	
	//find the shortest route
	//-1 means no route is feasible
	double val = taxiPath->value(nextVertex, source, dest);
	return val;
}

void taxiTree :: cancel() {
	taxiPath->cancel();
}

void taxiTree :: push() {
	taxiPath->push();
	
	//update destination
	targetPath = taxiPath->curr(nextVertex);
	
	//in case we were on a random destination earlier, we unset the flag
	randDest = false;
	distToDest = 0;
}

int taxiTree :: getNumberNodes() {
	return taxiPath->getNumberNodes();
}

bool taxiTree :: dynamicConstraints() {
	return taxiPath->dynamicConstraints();
}

void taxiTree :: setConstraints(double pickup_constraint, double service_constraint) {
	taxiPath->setConstraints(pickup_constraint, service_constraint);
}
