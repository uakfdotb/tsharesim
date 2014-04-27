#ifndef BRUTEPATH_H
#define BRUTEPATH_H

struct PathPoint;
class BrutePath;
class TaxiPath;

class BrutePath : public PathHandler {
private:
	vector<PathPoint *> tempPoints;
	map<int, double> positions;

public:
	BrutePath(TaxiPath *taxiPath, int M, vector<PathPoint *> &points);
	~BrutePath();
	virtual vector<PathPoint *> findPath(vertex *curr);
	virtual void init( ) {};

private:
	bool executeStep(vertex *curr);
};

#endif
