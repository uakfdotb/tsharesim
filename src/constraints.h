#ifndef CONSTRAINTS_H
#define CONSTRAINTS_H

using namespace std;

//this file contains constraints that are used throughout
// the simulation program for various purposes
//for details on the function of each constraint, see
// simulator.cpp, where they are defined

extern int TEST_MODE;
extern int EUCLIDEAN_PATH;
extern bool LOG_PATHPOINTS;

extern string inDirectory;

extern int MAX_PASSENGERS;
extern int MAX_TAXI;
extern const int HASHER;
extern const int HASHINGLONG;
extern const int HASHINGLAT;
extern int SPATIAL_INDEX;
extern int SPATIAL_LIMIT;

extern int ALGORITHM;

extern double SERVICE_CONST;
extern double PICKUP_CONST;

extern const double D;
extern const double LONG_LOW_BOUND;
extern const double LONG_HIGH_BOUND;
extern const double LAT_LOW_BOUND;
extern const double LAT_HIGH_BOUND;
extern const double LONG_RANGE;
extern const double LAT_RANGE;

extern int logicalTime;
extern int logicalTimeLimit;

#endif
