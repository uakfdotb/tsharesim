#include "includes.h"
#include "path.h"
#include "shortestPath.h"
#include "vertex.h"
#include "taxiPath.h"

TaxiPath :: TaxiPath(ShortestPath *shortestPath, vertex *curr) {
	this->shortestPathC = shortestPath;
	this->curr_vert = curr;
}

TaxiPath :: ~TaxiPath() {

}

double TaxiPath :: shortestPath(vertex *a, vertex *b) {
	return shortestPathC->shortestDistance(a, b);
}
