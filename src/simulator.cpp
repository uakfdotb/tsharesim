#include "includes.h"
#include "taxitree.h"
#include "customerqueue.h"
#include "vertex.h"
#include "shortestPath.h"
#include "taxilinklist.h"
#include "oneTaxiPath.h"

#include <csignal>

//BEGIN CONFIGURATION

//whether testing mode is enabled
int TEST_MODE = 0;

//enabling this will make all shortest path calls return a path directly to the destination vertex
int EUCLIDEAN_PATH = 0;

//whether to log PathPoint cases for later analysis with TEST_MODE=2
bool LOG_PATHPOINTS = false;

//the directory to grab input files, including the road network and customer data
string inDirectory = "../first_graph";

//maximum number of passengers assigned to a taxi at a time
int MAX_PASSENGERS = 4;

//number of taxis to simulate
int MAX_TAXIS = 20000;

//number of grid cells per unit distance (in this case, GPS coordinates)
const int HASHER = 200;

//minimim longitude for normalization
const int HASHINGLONG = 120500000;

//minimum latitude for normalization
const int HASHINGLAT = 30500000;

//number of grid cells to search in each direction for taxis
// this is automatically calculated
int SPATIAL_INDEX = 20;

//2 * SPATIAL_INDEX
// this is automatically calculated
int SPATIAL_LIMIT = 40;

//choice of algorithm
//0: branch and bound
//1: brute force
//2: lpsolve
//3: default tree
//4: tree with slack time
//5: tree with slack time, clustering
int ALGORITHM = 0;

//constraint on the service time for a passenger
// service time is the time that the passenger is in the taxi
double SERVICE_CONST = 1.2;

//constraint on the pickup time for a passenger
// pickup time is the time between request submission and passenger pickup
double PICKUP_CONST = 0.08448;

//taxi speed, measured in GPS coordinates per second
const double D = 0.0001408;

//these fields aren't used anymore, apparently
const double LONG_LOW_BOUND = 120.5;
const double LONG_HIGH_BOUND = 122.5;
const double LAT_LOW_BOUND = 30.5;
const double LAT_HIGH_BOUND = 32.0;
const double LONG_RANGE = (LONG_HIGH_BOUND-LONG_LOW_BOUND)/RAND_MAX;
const double LAT_RANGE = (LAT_HIGH_BOUND-LAT_LOW_BOUND)/RAND_MAX;

//start and end time for the simulation
//for the entire dataset, use 0 (midnight) to 63000 (17:30)
int logicalTime = 25200;
int logicalTimeLimit = 28800;

//END CONFIGURATION

//queue of customer requests waiting to be processed
customerqueue customers;

//set of taxis; hopefully 50000 >= MAX_PASSENGERS
taxiTree taxi[50000];

//structure for grid filtering
taxiLinkList spatialHash[HASHER*2][HASHER*2];

//this stores a list of the last several customer response time
// these are time needed to process one customer
deque<long> customerResponseTimes;

//list of last several iteration times
// these are time needed to advance one simulation second
deque<long> iterationTimes;

//contains number of passengers in each taxi
vector<int> passengerSizes;

//list of last several passenger value times
// these are time needed to call value() for specific passenger numbers
deque<long> *passengerValueTimes;

//total number of customers who could not be assigned to any taxi
int unsatisfiedCustomers = 0;

//array of vertices in the graph
vector<vertex *> vertices;

//variables to help print out statistics
double avgTopPassengers = 0.0;
double avgPassengers = 0.0;
int maxPassengers = 0;
int countNoPassengers = 0;

void customerSimulator();
void taxiSimulator();
void updateTaxis();

//used to sort stuff in descending order
bool compdesc (int i,int j) { return (i>j); }


//////////////////////////////////////////////////////////////////////////////
//
// from: http://stackoverflow.com/questions/669438/how-to-get-memory-usage-at-run-time-in-c
// this will only work on Unix-based systems
//
// process_mem_usage(double &, double &) - takes two doubles by reference,
// attempts to read the system-dependent data for a process' virtual memory
// size and resident set size, and return the results in KB.
//
// On failure, returns 0.0, 0.0

void process_mem_usage(double& vm_usage, double& resident_set)
{
	using std::ios_base;
	using std::ifstream;
	using std::string;

	vm_usage		 = 0.0;
	resident_set = 0.0;

	// 'file' stat seems to give the most reliable results
	//
	ifstream stat_stream("/proc/self/stat",ios_base::in);

	// dummy vars for leading entries in stat that we don't care about
	//
	string pid, comm, state, ppid, pgrp, session, tty_nr;
	string tpgid, flags, minflt, cminflt, majflt, cmajflt;
	string utime, stime, cutime, cstime, priority, nice;
	string O, itrealvalue, starttime;

	// the two fields we want
	unsigned long vsize;
	long rss;

	stat_stream >> pid >> comm >> state >> ppid >> pgrp >> session >> tty_nr
				 >> tpgid >> flags >> minflt >> cminflt >> majflt >> cmajflt
				 >> utime >> stime >> cutime >> cstime >> priority >> nice
				 >> O >> itrealvalue >> starttime >> vsize >> rss; // don't care about the rest

	stat_stream.close();

	long page_size_kb = sysconf(_SC_PAGE_SIZE) / 1024; // in case x86-64 is configured to use 2MB pages
	vm_usage = vsize / 1024.0;
	resident_set = rss * page_size_kb;
}

int main(int argc, char* argv[])
{
	//set parameters based on command line arguments
	
	if(argc >= 2) {
		ALGORITHM = atoi(argv[1]);
		
		if(argc >= 3) {
			MAX_PASSENGERS = atoi(argv[2]);
		
			if(argc >= 4) {
				MAX_TAXIS = atoi(argv[3]);
			
				if(argc >= 5) {
					double factor = (double) atoi(argv[4]) / 100;
					SERVICE_CONST = 1 + (SERVICE_CONST - 1) * factor;
					PICKUP_CONST *= factor;
				}
			}
		}
	}
	
	cout << "ALGORITHM=" << ALGORITHM << endl;
	cout << "MAX_PASSENGERS=" << MAX_PASSENGERS << endl;
	cout << "MAX_TAXIS=" << MAX_TAXIS << endl;
	cout << "SERVICE_CONST=" << SERVICE_CONST << endl;
	cout << "PICKUP_CONST=" << PICKUP_CONST << endl;
	
	//calculate the spatial index and spatial limit
	SPATIAL_INDEX = (int) (PICKUP_CONST * HASHER) + 1;
	SPATIAL_LIMIT = SPATIAL_INDEX * 2;
	
	cout << "SPATIAL_INDEX=" << SPATIAL_INDEX << endl;
	cout << "SPATIAL_LIMIT=" << SPATIAL_LIMIT << endl;

	//seed random number generator
	srand(50);
	
	passengerValueTimes = new deque<long>[MAX_PASSENGERS];
	
	//read vertices
	ifstream vertexFile;
	vertexFile.open((inDirectory + "/vertices.out").c_str());
	
	int n;
	vertexFile >> n;
	
	for(int i = 0; i < n; i++) {
		vertex *newVertex = new vertex;
		vertexFile >> newVertex->x;
		vertexFile >> newVertex->y;
		newVertex->id = i;
		
		vertices.push_back(newVertex);
	}
	
	vertexFile.close();
	
	//read edges
	ifstream edgeFile;
	edgeFile.open((inDirectory + "/edges.out").c_str());
	
	edgeFile >> n;
	
	for(int i = 0; i < n; i++) {
		int a, b;
		edgeFile >> a >> b;
		
		vertices[a]->neighbors.push_back(vertices[b]);
		vertices[b]->neighbors.push_back(vertices[a]);
		
		//ignore distance, we recalculate since it's just Euclidean
		double distance;
		edgeFile >> distance;
	}
	
	edgeFile.close();

	//initialize shortest path interface
	ShortestPath *shortestPath = new ShortestPath(vertices);
	
	//some basic testing simulations
	if(TEST_MODE) {
		struct timeval testStartTime;
		struct timeval testEndTime;
		gettimeofday(&testStartTime, NULL);
		
		if(TEST_MODE == 1) {
			//create a taxi and feed it randomly generated requests
			//look mostly at value() calls with MAX_PASSENGERS by moving the taxi
			// once it has reached capacity to the first dropoff point
			
			taxiTree testTaxi;
			testTaxi.setInitialPosition(shortestPath, vertices[0]);
			int count = 0;
		
			for(int k = 0; k < 1000000; k++) {
				int i = rand() % (k/100 + 100);
				int j = rand() % (k/100 + 100);
				if(i == j) continue;
			
				double value = testTaxi.value(vertices[i], vertices[j]);
				if(value != -2) count++;
			
				if(value >= 0) {
					cout << "pushing " << i << ", " << j << endl;
					testTaxi.push();
					testTaxi.passengers++;
					testTaxi.printPoints();
					
					//update location until we drop a passenger off
					while(testTaxi.passengers >= MAX_PASSENGERS) {
						testTaxi.updateLocation();
						//cout << "update location, " << testTaxi.passengers << endl;
					}
				
					cout << "push finished (" << testTaxi.passengers << ")" << endl;
				} else {
					//cout << "canceling with value=" << value << endl;
					testTaxi.cancel();
				}
			}
		
			testTaxi.printPoints();
			
			cout << "count=" << count << endl;
		} else if(TEST_MODE == 2) {
			//from standard input, read oneTaxiPath PathPoints and calculate shortest routes
			//this is also compatible with the basic tree algorithm (ALGORITHM=3)
			
			if(ALGORITHM != 0 && ALGORITHM != 1 && ALGORITHM != 2) {
				cout << "invalid algorithm, not compatible with TEST_MODE=2" << endl;
			}
			
			string line;
			OneTaxiPath ot(shortestPath, vertices[0]);
			
			while(getline(cin, line)) {
				vertex *curr;
				vector<PathPoint *> points;
				
				if(line.length() >= 5 && line[0] == 'p' && line[1] == 'p' && line[2] == 'l' && line[3] == 'o' && line[4] == 'g') {
					istringstream iss(line);
					string pathPointString;
					
					getline(iss, pathPointString, ' ');
					getline(iss, pathPointString, ' ');
					curr = vertices[atoi(pathPointString.c_str())];
					
					while(getline(iss, pathPointString, ' ')) {
						istringstream pathStream(pathPointString);
						string partString;
						PathPoint *point = new PathPoint;
						
						getline(pathStream, partString, '/');
						point->vert = vertices[atoi(partString.c_str())];
						getline(pathStream, partString, '/');
						point->type = atoi(partString.c_str());
						getline(pathStream, partString, '/');
						point->pairIndex = atoi(partString.c_str());
						getline(pathStream, partString, '/');
						point->remaining = atof(partString.c_str());
						
						points.push_back(point);
					}
				}
				
				ot.test(curr, points);
			}
		}
		
		gettimeofday(&testEndTime, NULL);
		long tS = testStartTime.tv_sec*1000000 + (testStartTime.tv_usec);
		long tE = testEndTime.tv_sec*1000000 + (testEndTime.tv_usec);
		long testTime = tE - tS;
		
		cout << "completed in " << testTime << endl;
		return 0;
	}

	//initialize set of taxis to random positions within the city

	for(int i = 0; i < MAX_TAXIS; i++) {
		//select a random vertex to place this taxi
		int randIndex = rand() % vertices.size();
		taxi[i].setInitialPosition(shortestPath, vertices[randIndex]);
		
		//calculate the initial position within the grid for this taxi
		int hashLong = ((((int) (taxi[i].current.longitude*1000000))%HASHINGLONG)/(1000000/HASHER));
		int hashLat = ((((int) (taxi[i].current.latitude*1000000))%HASHINGLAT)/(1000000/HASHER));
		spatialHash[hashLong][hashLat].append(i);
		taxi[i].hashLong = hashLong;
		taxi[i].hashLat = hashLat;
	}

	//run the simulation
	customerSimulator();

	//print any customers remaining in the queue and some statistics
	// the customers only come from last iteration so this output isn't very useful
	customers.display();
	fprintf(stdout, "Logical time = %d\n", logicalTime);
	fprintf(stdout, "Avg Q time=%f\n Serviced=%d\n Queued=%d\n", 
		customers.totalQtime/customers.assignedCount, 
		customers.assignedCount, 
		customers.customerCount);

	return 0;
}

void updateTaxis()
{
	//update the position in the grid of all the taxis
	//also calculate statistics on them
	
	int hashLong, hashLat;
	int passengerTotal;

	passengerTotal = 0;
	maxPassengers = 0;
	countNoPassengers = 0;
	
	passengerSizes.clear();
	
	for(int x = 0; x < MAX_TAXIS; x++) {
		passengerTotal += taxi[x].passengers;
		passengerSizes.push_back(taxi[x].passengers);
		
		maxPassengers = max(maxPassengers, taxi[x].passengers);
		
		if(taxi[x].passengers == 0) {
			countNoPassengers++;
		}
		
		//recalculate the position of this taxi in the grid
		hashLong = ((((int) (taxi[x].current.longitude*1000000))%HASHINGLONG)/(1000000/HASHER));
		hashLat = ((((int) (taxi[x].current.latitude*1000000))%HASHINGLAT)/(1000000/HASHER));
		if(hashLong > HASHER*2-1)
			hashLong = 0;
		if(hashLat > HASHER*2-1)
			hashLat = 0;
		
		bool result = spatialHash[taxi[x].hashLong][taxi[x].hashLat].del(x);
		
		if(!result) {
			cout << x << taxi[x].hashLong << "," << taxi[x].hashLat << "; " << hashLong << "," << hashLat << ";" << taxi[x].current.longitude << "," << taxi[x].current.latitude << endl;
			raise(SIGSEGV);
		}
		
		spatialHash[hashLong][hashLat].append(x);
		taxi[x].hashLong = hashLong;
		taxi[x].hashLat = hashLat;
	}
	
	//find average number of passengers in top 20% taxis
	sort(passengerSizes.begin(), passengerSizes.end(), compdesc);
	
	int passengerTopSize = 0;
	for(int i = 0; i < passengerSizes.size() / 5; i++) {
		passengerTopSize += passengerSizes[i];
	}
	
	avgTopPassengers = (double) passengerTopSize / (passengerSizes.size() / 5);
	avgPassengers = (double) passengerTotal/MAX_TAXIS;
}

//todo: misleading name
void customerSimulator()
{
	//this is the main simulation loop
	//customers are added into the customer queue in this function
	//then, other functions are called to update the taxi and process the requests
	
	int nextID = 1; //the next customer ID to use
	int hashLong, hashLat;
	customer *newCustomer = new customer;
	
	//open the customer request file
	ifstream customerFile;
	customerFile.open((inDirectory + "/newtaxi.dat").c_str());
	customerFile >> newCustomer->time;
	
	while(logicalTime < logicalTimeLimit) {
		//get every customer for this simulated second
		while(!customerFile.eof() && newCustomer->time <= logicalTime) {
			int vertexIndex;
			
			//read the pickup and dropoff points of this request from the file
			customerFile >> vertexIndex;
			newCustomer->pickup = vertices[vertexIndex];
			
			customerFile >> vertexIndex;
			newCustomer->dropoff = vertices[vertexIndex];
			
			newCustomer->hashLong = 
				((((int) (newCustomer->pickup->x*1000000))%HASHINGLONG)/(1000000/HASHER));
			newCustomer->hashLat = 
				((((int) (newCustomer->pickup->y*1000000))%HASHINGLAT)/(1000000/HASHER));
			
			//if logicalTime != 0, then we might not want to
			// store this customer; in that case, we can reuse
			// the allocated newCustomer memory for the next
			// customer
			if(newCustomer->time == logicalTime) {
				customers.insert(newCustomer);
				newCustomer = new customer;
			}
			
			nextID += 2;
			customerFile >> newCustomer->time;
		}
		
		//calculate average customer response time and average iteration time
		long totalCustomerResponseTime = 0;
		long totalIterationTime = 0;
		
		deque<long>::iterator itr;
		for(itr = customerResponseTimes.begin(); itr != customerResponseTimes.end(); itr++) {
			totalCustomerResponseTime += *itr;
		}
		for(itr = iterationTimes.begin(); itr != iterationTimes.end(); itr++) {
			totalIterationTime += *itr;
		}
		
		double averageCustomerResponseTime = (double) totalCustomerResponseTime / customerResponseTimes.size();
		double averageIterationTime = (double) totalIterationTime / iterationTimes.size();
		
		//some other statistics...
		double residentSetSize;
		double vmUsage;
		process_mem_usage(residentSetSize, vmUsage);
		
		fprintf(stdout, "it=%d  %f %1.5f %1.5f %d %d %10.5f %10.5f %1.1f %1.1f %d", logicalTime, customers.totalQtime/customers.assignedCount, avgPassengers, avgTopPassengers, maxPassengers, countNoPassengers, averageCustomerResponseTime, averageIterationTime, residentSetSize, vmUsage, unsatisfiedCustomers);
		
		for(int i = 0; i < MAX_PASSENGERS; i++) {
			long passengerValueTotal = 0;
			for(itr = passengerValueTimes[i].begin(); itr != passengerValueTimes[i].end(); itr++) {
				passengerValueTotal += *itr;
			}
			
			double passengerValueAverage = (double) passengerValueTotal / passengerValueTimes[i].size();
			cout << " " << passengerValueAverage;
		}
		
		cout << endl;
		
		//process the requests and update taxi positions
		taxiSimulator();
		
		//update hash grid and calculate some statistics
		updateTaxis();

		logicalTime++;
	}
	
	customerFile.close();
}

void taxiSimulator()
{
	int i; //the current index in customer queue
	int minTaxi; //the taxi with the shortest route
	
	//these define the box to search for taxis in grid filtering
	int startLong;
	int endLong;
	int startLat;
	int endLat;
	
	node *ptr; //the current node (position in linked list at current grid location)
	customer *newCustomer;
	
	//update the taxi positions
	for(int x = 0; x < MAX_TAXIS; x++) {
		taxi[x].updateLocation();
	}
	
	struct timeval iterationStartTime;
	struct timeval iterationEndTime;
	gettimeofday(&iterationStartTime, NULL);
	
	//todo: get rid of i because it's always zero
	i = 0;
	
	while((newCustomer = customers.peak(i)) != NULL && i < customers.customerCount) {
		struct timeval customerStartTime;
		struct timeval customerEndTime;
		gettimeofday(&customerStartTime, NULL);
		
		//values for finding the taxi with the minimum path for this request
		minTaxi = -1;
		double minValue = -1;
		
		//find taxis close to the request's pickup point using grid filtering
		startLong = newCustomer->hashLong-SPATIAL_INDEX;
		endLong = startLong+SPATIAL_LIMIT;
		startLat = newCustomer->hashLat-SPATIAL_INDEX;
		endLat = startLat+SPATIAL_LIMIT;
		if(startLong < 0)
			startLong = 0;
		else if(startLong > HASHER*2-1)
			startLong = HASHER*2-1;
		if(endLong < 0)
			endLong = 0;
		else if(endLong > HASHER*2-1)
			endLong = HASHER*2-1;
		if(startLat < 0)
			startLat = 0;
		else if(startLat > HASHER*2-1)
			startLat = HASHER*2-1;
		if(endLat < 0)
			endLat = 0;
		else if(endLat > HASHER*2-1)
			endLat = HASHER*2-1;
		
		//loop through the grid
		//each grid cell has a linked list of taxis
		for(int x = startLong; x <= endLong; x++) {
			for(int y = startLat; y <= endLat; y++) {
				//go through the linked list and try each taxi that's not full
				ptr = spatialHash[x][y].returnHead();
				
				while(ptr != NULL) {
					if(taxi[ptr->vehicle].passengers < MAX_PASSENGERS) {
						struct timeval valueStartTime;
						struct timeval valueEndTime;
						gettimeofday(&valueStartTime, NULL);
						
						double value = taxi[ptr->vehicle].value(newCustomer->pickup, newCustomer->dropoff);
						
						//calculate time needed to process this request
						gettimeofday(&valueEndTime, NULL);
						long tS = valueStartTime.tv_sec*1000000 + (valueStartTime.tv_usec);
						long tE = valueEndTime.tv_sec*1000000	+ (valueEndTime.tv_usec);
						long valueTime = tE - tS;
						passengerValueTimes[taxi[ptr->vehicle].passengers].push_back(valueTime);
						
						if((minTaxi == -1 || value < minValue) && value >= 0) { //this taxi is better than our current best
							if(minTaxi != -1) { //we do have a current best
								//let the current best know that it will not be chosen
								taxi[minTaxi].cancel();
							}
							
							//update minimum information
							minTaxi = ptr->vehicle;
							minValue = value;
						} else {
							//make sure this taxi knows it won't be chosen
							taxi[ptr->vehicle].cancel();
						}
					}
					
					ptr = ptr->link;
				}
			}
		}
		
		//check whether we found any taxi
		
		if(minTaxi != -1) {
			//push the request onto the minimum-path taxi
			taxi[minTaxi].passengers++;
			taxi[minTaxi].push();
		} else {
			//increase unsatisfied count for output
			unsatisfiedCustomers++;
		}
		
		//remove this customer from the queue and update queue
		//the customer might have been accepted or rejected, but
		// either way they get an instant response
		newCustomer = customers.remove(i);
		newCustomer->link = NULL;
		delete newCustomer;
		newCustomer = NULL;
		
		gettimeofday(&customerEndTime, NULL);
		
		//calculate time needed to process this request
		long tS = customerStartTime.tv_sec*1000000 + (customerStartTime.tv_usec);
		long tE = customerEndTime.tv_sec*1000000	+ (customerEndTime.tv_usec);
		long customerTime = tE - tS;
		customerResponseTimes.push_back(customerTime);
	}
	
	gettimeofday(&iterationEndTime, NULL);
	
	//calculate time taken for this entire iteration
	long tS = iterationStartTime.tv_sec*1000000 + (iterationStartTime.tv_usec);
	long tE = iterationEndTime.tv_sec*1000000	+ (iterationEndTime.tv_usec);
	long iterationTime = tE - tS;
	iterationTimes.push_back(iterationTime);
	
	//limit sizes of the time history to 1000
	while(iterationTimes.size() > 1000) iterationTimes.pop_front();
	while(customerResponseTimes.size() > 1000) customerResponseTimes.pop_front();
	
	for(int i = 0; i < MAX_PASSENGERS; i++) {
		while(passengerValueTimes[i].size() > 10000) passengerValueTimes[i].pop_front();
	}
}
