#ifndef INCLUDES_H
#define INCLUDES_H

//this header is included by most of the source files and
// contains general includes that most will use

#include <stdint.h>

//the two includes below are now no longer required because
// the simulator is no longer multithreaded
//#include <pthread.h>
//#include <semaphore.h>

#include <fstream>
#include <algorithm>
#include <deque>
#include <queue>
#include <map>
#include <vector>

#include <unistd.h>
#include <ios>
#include <iostream>
#include <string>
#include <math.h>
#include <sys/time.h>
#include <cstdlib>
#include <limits>
#include <math.h>

#include "constraints.h"

using namespace std;

#endif
