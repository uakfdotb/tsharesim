#include "includes.h"
#include "taxiPath.h"
#include "oneTaxiPath.h"
#include "path.h"

PathHandler :: PathHandler(TaxiPath *taxiPath, int M, vector<PathPoint *> &nPoints): points(nPoints) {
	this->taxiPath = taxiPath;
	minTour = numeric_limits<double>::max();
}

PathHandler :: ~PathHandler() {

}

void PathHandler :: reset() {
	this->N = this->points.size();
	bestList.clear();
	minTour = numeric_limits<double>::max();
}

double PathHandler :: length(PathPoint *v1, PathPoint *v2) {
	return taxiPath->shortestPath(v1->vert, v2->vert);
}

double PathHandler :: lengthVertex(vertex *vert, PathPoint *v2) {
	return taxiPath->shortestPath(vert, v2->vert);
}
