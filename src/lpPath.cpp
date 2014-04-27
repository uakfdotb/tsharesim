#include "lpsolve/lp_lib.h"
#include "includes.h"
#include "path.h"
#include "taxiPath.h"
#include "oneTaxiPath.h"
#include "lpPath.h"

long lpTotalTime;
long lpMaxTime;
int lpNum;

LPPath :: LPPath(TaxiPath *taxiPath, int M, vector<PathPoint *> &nPoints): PathHandler(taxiPath, M, nPoints) {
	costMatrix = new double*[M + 1];
	earliest = new double[M];
	latest = new double[M + 1];
	m = new double*[M + 1];
	pairIndex = new int[M];
	
	for(int i = 0; i < M + 1; i++) {
		costMatrix[i] = new double[M];
		m[i] = new double[M];
	}
}

LPPath :: ~LPPath() {
}

vector<PathPoint *> LPPath :: findPath(vertex *curr) {
	lprec *lp;
	
	int numDropoff = 0;
	int numDropped = 0;
	int numPickup = 0;
	
	//find pairs for each dropoff point
	for(int i = 0; i < points.size(); i++) {
		if(points[i]->type == 1) {
			bool foundPair = false;
			
			for(int j = 0; j < points.size(); j++) {
				if(j != i && points[j]->pairIndex == points[i]->pairIndex) {
					pairIndex[i] = j;
					foundPair = true;
					break;
				}
			}
			
			//sometimes, there's an error and the pair cannot be found
			//in that case, print out some debugging information
			if(!foundPair) {
				cout << i << ":" << points[i]->pairIndex << "  ";
				for(int j = 0; j < points.size(); j++) {
					cout << points[j]->type << ":" << points[j]->pairIndex << " ";
				}
				cout << endl;
			}
		}
	}
	
	//occasionally we encounter a model that takes hours or days to solve
	//we set a timeout on the solve function, and then advance to the next iteration
	//as the iteration increases, we introduce more randomness into the model
	// (this is done via the getNonZero function)
	for(int iteration = 0; ; iteration += 10) {
		//calculate cost matrix
		for(int i = 0; i < points.size(); i++) {
			PathPoint *ipoint = points[i];
		
			if(ipoint->type == 0) numPickup++;
			else if(ipoint->type == 1) numDropoff++;
			else if(ipoint->type == 2) numDropped++;
			
			//from this point to itself
			costMatrix[i + 1][i] = getNonZero(0, iteration);
			
			//from this point to every other point
			for(int j = 0; j < points.size(); j++) {
				if(i != j)
					costMatrix[i + 1][j] = getNonZero(length(ipoint, points[j]), iteration);
			}
			
			//from the current location to this point
			costMatrix[0][i] = getNonZero(taxiPath->shortestPath(curr, ipoint->vert), iteration);
		}

	
		//calculate m matrix
		//first, we have to find earliest and latest
		
		//the current location must occur at time zero
		latest[0] = 0;
	
		for(int i = 0; i < points.size(); i++) {
			if(points[i]->type == 0 || points[i]->type == 2) {
				//this is a pickup or stand-alone dropoff point
				//the earliest time occurs when we go directly
				// from the current location to here
				//the latest time is set by the pickup constraint
				
				earliest[i] = costMatrix[0][i];
				latest[i + 1] = points[i]->remaining;
			} else if(points[i]->type == 1) {
				//this is a dropoff point
				//the earliest time occurs when we go directly
				// to the pickup point, then here
				//the latest time occurs when we get to the pickup
				// point the latest, and then here the latest
				// (stretch both pickup and service constraints)
				earliest[i] = costMatrix[0][pairIndex[i]] + costMatrix[pairIndex[i] + 1][i];
				latest[i + 1] = points[pairIndex[i]]->remaining + points[i]->remaining;
			}
		}
		
		//calculate m
		double test;
		for(int i = 0; i < points.size() + 1; i++) {
			for(int j = 0; j < points.size(); j++) {
				test = latest[i] + costMatrix[i][j] - earliest[j];
				if(test > 0) m[i][j] = test;
				else m[i][j] = 0;
			}
		}
		
		//find the number of binary columns
		//each x_ij determines whether or not the taxi will move
		// from i to j
		//in the comments below these movements will be referred
		// to as route segments (_from_ i _to_ j)
		int ncol = (points.size() + 1) * points.size();
		
		//find the total number of columns
		//besides the binary ones, there are ones for the time
		// at which the taxi will reach a point (B_i)
		int ncol_total = ncol + points.size() + 1;
		
		//create the lp instance
		lp = make_lp(0, ncol_total);
		
		//colno and row are used to define the constraints, and
		// later row will store the result from lpsolve
		//colno identifies the variable (column), and row identifies
		// the constants (multiplied by the variable); then, a
		// separate value determines the number of variables
		// that will be read (since we are using a sparse matrix -
		// otherwise we wouldn't need colno)
		//note**: column numbers are labeled starting from 1, not 0
		int *colno = new int[ncol_total];
		REAL *row = new REAL[ncol_total];
		
		//since we're going to be adding constraints equation
		// by equation, we set add row mode to make it faster
		set_add_rowmode(lp, TRUE);
		
		//disable most output from lpsolve
		set_verbose(lp, CRITICAL);
		
		//set timeout of three seconds so we don't spend forever on this model
		set_timeout(lp, 3);
		
		//set up the binary constraints
		for(int i = 0; i < ncol; i++) {
			set_binary(lp, i + 1, TRUE);
		}
		
		//constraints 1 to 3
		//these have one constraint per point
		for(int i = 0; i < points.size(); i++) {
			//1. the total number of route segments to here will
			// be equal to one
			for(int j = 0; j < points.size() + 1; j++) {
				colno[j] = j * points.size() + i + 1;
				row[j] = 1;
			}
			
			add_constraintex(lp, points.size() + 1, row, colno, EQ, 1);
			
			//2. there will be no route segment from here to itself
			colno[0] = (i + 1) * points.size() + i + 1;
			row[0] = 1;
			add_constraintex(lp, 1, row, colno, EQ, 0);
			
			//3. the total number of route segments from here will
			// be less than or equal to one (since the last point in
			// the route will be zero)
			for(int j = 0; j < points.size(); j++) {
				colno[j] = (i + 1) * points.size() + j + 1;
				row[j] = 1;
			}
			
			add_constraintex(lp, points.size(), row, colno, LE, 1);
		}
		
		//4. there will be exactly one route segment from the
		// current location
		for(int i = 0; i < points.size(); i++) {
			colno[i] = i + 1;
			row[i] = 1;
		}
	
		add_constraintex(lp, points.size(), row, colno, EQ, 1);
	
		//5. the relative time that the taxi reaches the current
		// location is zero
		colno[0] = ncol + 1;
		row[0] = 1;
		add_constraintex(lp, 1, row, colno, EQ, 0);
	
		//6. defined for each route segment (i, j)
		//if the segment (i, j) exists, then the time B_j
		// the taxi reaches j will be greater than
		//    B_i + time(i, j)
		// (time is interchangeable with distance)
		//in other words,
		//    B_j >= ( B_i + time(i, j) ) * x_ij
		//
		//**but that's non-linear (since B_i * x_ij)
		//to achieve the if statement, we subtract a large
		// number M from the right and M * x_ij on the left
		//the equation becomes:
		//    B_j - B_i - M*x_ij >= time(i, j) - M
		//
		//m_ij that we found earlier is suitable for M, since
		// if x_ij = 0 the equation reduces to
		//    B_j - B_i >= time(i, j) - M
		// >> M >= B_i + time(i, j) - B_j
		// we used the maximum possible value for B_i (latest[i])
		//  and the minimim for B_j (earliest[j]), so everything
		//  is good :)
		for(int i = 0; i < points.size() + 1; i++) {
			for(int j = 0; j < points.size(); j++) {
				colno[0] = ncol + 1 + i;
				colno[1] = ncol + 1 + j + 1; //make sure to add an extra 1 because we're not including current location
				colno[2] = i * points.size() + j + 1;
			
				double constant = costMatrix[i][j] - m[i][j];
			
				//only use positive constants or it seems to explode
				if(constant >= 0) {
					row[0] = -1;
					row[1] = 1;
					row[2] = -m[i][j];
		
					add_constraintex(lp, 3, row, colno, GE, constant);
				} else {
					row[0] = 1;
					row[1] = -1;
					row[2] = m[i][j];
		
					add_constraintex(lp, 3, row, colno, LE, -constant);
				}
			}
		}
	
		//constraints 7, 8, and 9
		for(int i = 0; i < points.size(); i++) {
			if(points[i]->type == 1) {
				//dropoff point
				
				//make sure to add an extra 1 because we're not including current location
				colno[0] = ncol + 1 + i + 1;
				colno[1] = ncol + 1 + pairIndex[i] + 1;
			
				row[0] = 1;
				row[1] = -1;
			
				//constraints on L_i (= B_i - B_pickup[i])
				
				//7. L_i >= time(pickup[i], i)
				add_constraintex(lp, 2, row, colno, GE, costMatrix[pairIndex[i] + 1][i]);
				
				//8. L_i <= remaining service constraint
				add_constraintex(lp, 2, row, colno, LE, points[i]->remaining);
			} else if(points[i]->type == 0 || points[i]->type == 2) {
				//pickup or stand-alone dropoff point
				colno[0] = ncol + 1 + i + 1;
				row[0] = 1;
				
				//9. B_i <= remaining pickup constraint
				add_constraintex(lp, 1, row, colno, LE, points[i]->remaining);
			}
		}
	
		//10. this used to enforce that all varibles be
		// non-negative, but it seems to be working now
		// (lpsolve makes variables non-negative unless
		// explicitly stated in a constraint)
		for(int i = ncol; i < ncol_total; i++) {
			colno[0] = i + 1;
			row[0] = 1;
			//add_constraintex(lp, 1, row, colno, GE, 0);
		}
		
		//disable rowmode because we're done building model
		set_add_rowmode(lp, FALSE);
		
		//objective function: minimize sum( time(i, j) * x_ij )
		//we maximize the negative though
		// (we could change to set_minim(lp), but it's the same thing)
		for(int i = 0; i < points.size() + 1; i++) {
			for(int j = 0; j < points.size(); j++) {
				colno[i * points.size() + j] = i * points.size() + j + 1;;
				row[i * points.size() + j] = -costMatrix[i][j];
			}
		}
	
		set_obj_fnex(lp, ncol, row, colno);
		set_maxim(lp); //maximize the objective function
		
		struct timeval solveStartTime;
		struct timeval solveEndTime;
		gettimeofday(&solveStartTime, NULL);
		
		int ret = solve(lp);
		
		gettimeofday(&solveEndTime, NULL);
		long tS = solveStartTime.tv_sec*1000000 + (solveStartTime.tv_usec);
		long tE = solveEndTime.tv_sec*1000000 + (solveEndTime.tv_usec);
		long solveTime = tE - tS;
		
		if(iteration == 0 && ret != TIMEOUT) {
			lpTotalTime += solveTime;
			if(solveTime > lpMaxTime) lpMaxTime = solveTime;
			lpNum++;
			
			cout << "lptimestatus: " << lpTotalTime / lpNum << " " << lpMaxTime << " " << lpNum << " " << solveTime << endl;
		}
		
		//if we didn't get the optimal solution, don't continue
		if(ret != OPTIMAL) {
			delete colno;
			delete row;
			delete_lp(lp);
			bestList.clear();
			
			if(ret == TIMEOUT) {
				//if we timed out, then we need to try again
				cout << "timed out on iteration " << iteration << ", advancing..." << endl;
				continue;
			} else {
				return bestList;
			}
		}
	
		get_variables(lp, row); //store variables in our row array
	
		//extract the ordering of the points from the x_ij in the row
		//at the same time, we calculate the route's distance
		
		int previous = 0;
		minTour = 0;
		
		for(int i = 0; i < points.size(); i++) {
			for(int j = 0; j < points.size(); j++) {
				if(row[previous * points.size() + j] == 1) {
					minTour += costMatrix[previous][j];
				
					bestList.push_back(points[j]);
					previous = j + 1;
					break;
				}
			}
		}

		delete colno;
		delete row;
		delete_lp(lp);
		
		//sometimes errors occur, probably because M was
		// too large and double-precision isn't accurate
		// enough
		//in these cases, since they're rare enough, we
		// assume that the model was infeasible
		if(bestList.size() != points.size()) {
			bestList.clear();
			minTour = numeric_limits<double>::max();
			return bestList;
		}
	
		return bestList;
	}
}

double LPPath :: getNonZero(double d, int iteration) {
	if(d < 0.00001 + 0.00002 * iteration)
		return 0.00001;
	else {
		//the code in the if statement won't do anything if iteration = 0
		//nevertheless, we skip for efficiency
		if(iteration != 0) {
			double maxchange = 0.00002 * iteration;
			double change = ((double) (rand() % 1000)) * maxchange / 1000.0;
		
			//either increase or decrease by change
			if(rand() % 2 == 0) {
				d -= change;
			} else {
				d += change;
			}
		}
		
		return d;
	}
}
