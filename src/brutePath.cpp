#include "includes.h"
#include "path.h"
#include "taxiPath.h"
#include "oneTaxiPath.h"
#include "brutePath.h"

BrutePath :: BrutePath(TaxiPath *taxiPath, int M, vector<PathPoint *> &nPoints): PathHandler(taxiPath, M, nPoints) {
}

BrutePath :: ~BrutePath() {
}

vector<PathPoint *> BrutePath :: findPath(vertex *curr) {
	//create our own copy of the points because we'll be modifying it
	tempPoints = points;
	
	//iteratively generate all possible permutations of the points and evaluate
	//see http://en.wikipedia.org/wiki/Permutation#Generation_in_lexicographic_order for the algorithm
	while(!executeStep(curr));
	
	return bestList;
}

bool BrutePath :: executeStep(vertex *curr) {
	//test the solution for feasibility
	double currLength = 0;
	positions.clear();
	int it; //use int instead of iterator so that we can compare adjacent elements
	
	for(it = 0; it < tempPoints.size(); it++) {
		if(it > 0) {
			currLength += length(tempPoints[it - 1], tempPoints[it]);
		} else {
			//find distance from current location to first point
			currLength += lengthVertex(curr, tempPoints[it]);
		}
		
		if(tempPoints[it]->type == 0) {
			//this is a pickup point
			if(currLength > tempPoints[it]->remaining) {
				break; //pickup constraint exceeded
			} else {
				positions[tempPoints[it]->pairIndex] = currLength;
			}
		} else if(tempPoints[it]->type == 2) {
			//this is a stand-alone dropoff point
			if(currLength > tempPoints[it]->remaining) {
				break; //service constraint exceeded
			} else {
				positions[tempPoints[it]->pairIndex + 1] = -1;
			}
		} else {
			//this is a dropoff point
			map<int, double>::iterator pos = positions.find(tempPoints[it]->pairIndex);
		
			if(pos != positions.end()) {
				double pairLength = pos->second; //length at corresponding pickup point
	
				if(currLength - pairLength > tempPoints[it]->remaining) {
					break;
				} else {
					positions[tempPoints[it]->pairIndex + 1] = -1;
				}
			} else {
				break;
			}
		}
	}
	
	if(it == tempPoints.size() && currLength < minTour) { //feasible && better than current best
		minTour = currLength;
		bestList = tempPoints;
	}

	int k, l;
	
	//find the largest index k such that a[k] < a[k + 1]
	for(k = tempPoints.size() - 2; k >= 0; k--) {
		if(tempPoints[k]->index < tempPoints[k + 1]->index) {
			break;
		}
	}
	
	if(k < 0) {
		return true; //we're done
	}
	
	//find the largest index l such that a[k] < a[l]
	for(l = tempPoints.size() - 1; l >= 0; l--) {
		if(tempPoints[k]->index < tempPoints[l]->index) {
			break;
		}
	}
	
	if(l < 0) {
		return true; //this should never occur
	}
	
	//swap a[k] with a[l]
	PathPoint* tmp = tempPoints[k];
	tempPoints[k] = tempPoints[l];
	tempPoints[l] = tmp;
	
	//reverse the sequence from a[k + 1] up to and including the final element a[n]
	int numReverse = tempPoints.size() - 1 - k;
	for(int i = 0; i < numReverse / 2; i++) {
		tmp = tempPoints[k + 1 + i];
		tempPoints[k + 1 + i] = tempPoints[tempPoints.size() - i - 1];
		tempPoints[tempPoints.size() - i - 1] = tmp;
	}
	
	return false;
}
