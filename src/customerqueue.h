#include <stdlib.h>		 // standard c library
#include <stdio.h>			// standard c io
#include <unistd.h>
#include <sys/time.h>
#include "vertex.h"
using namespace std;

#ifndef CUSTOMERQUEUE_h__
#define CUSTOMERQUEUE_h__

//customer stores the information for one customer request
struct customer {
	//the simulated second at which this request [should be]/was made
	int time;
	
	//the index along the longitudinal access where this request's
	// pickup point falls in the grid
	int hashLong;
	
	//same as above, but for latitude
	int hashLat;
	
	//vertices where the pickup and dropoff points are closest to
	vertex *pickup;
	vertex *dropoff;
	
	//these will be set to the times at which the customer first
	// enters the queue and the time at which it is removed
	//because of how the simulator is programmed, their difference
	// is NOT the customer response time; it has no meaning and is
	// not used anymore
	//this is because all customers for one simulated second are
	// queued before their processing ever begins
	struct timeval start, end;
	
	//the next customer in the queue
	customer *link;
};

class customerqueue {
public:
	//the head of the queue
	customer *head;
	
	//the tail of the queue
	customer *tail;
	
	//the number of customers currently in the queue
	int customerCount;
	
	//the total number of customers that have been
	// removed from the queue
	int assignedCount;
	
	//the total time that customers have spent in the queue
	//sum of x.end - x.start, for each customer x
	double totalQtime;

	customerqueue();
	~customerqueue();
	
	//insert a new customer into the queue
	void insert(customer *toQueue);
	
	//remove a customer from the queue
	customer* remove(int i);
	
	//
	customer* peak(int i);
	
	//display general information about all customers
	// currently residing in the queue
	void display();
};

customerqueue::customerqueue() {
	head = NULL;
	tail = NULL;
	customerCount = 0;
	totalQtime = 0.0;
	assignedCount = 0;
}

customerqueue::~customerqueue() {
	//delete every customer that has not been removed yet
	customer *ptr;

	for(int i = 0; i < customerCount; i++) {
		ptr = remove(i);
		ptr->link = NULL;
		delete ptr;
	}
	
	//null the pointers just in case
	head = NULL;
	tail = NULL;
}

void customerqueue::insert(customer *newCustomer) {
	//make sure that the new customer's link is NULL
	newCustomer->link = NULL;
	
	//record the time at which the customer is entering the queue
	gettimeofday(&newCustomer->start, NULL);

	if(head == NULL) {
		//this customer simply becomes the head
		newCustomer->link = NULL;
		head = newCustomer;
		tail = head;
	}
	
	else {
		tail->link = newCustomer;
		tail = newCustomer;
	}
	
	customerCount++;

	return;
}

customer *customerqueue::remove(int i) {
	customer *tempPtr;

	if(!head) {
		//the queue is empty, there's nothing we can do
		return NULL;
	}
	
	else {
		customer *foundCustomer;
		
		if(i == 0) {
			//the caller wants the head
			//that means we have to set next item as the head
			tempPtr = head;
			head = head->link;
			tempPtr->link = NULL;
			foundCustomer = tempPtr;
		}
		
		else {
			//go down queue for i iterations
			tempPtr = head;
			
			for(int k = 0; k < i - 1 && tempPtr->link != NULL; k++) {
				tempPtr = tempPtr->link;
			}
			
			//now, the customer before foundCustomer must point
			// to the customer after foundCustomer
			foundCustomer = tempPtr->link;
			tempPtr->link = foundCustomer->link;
			foundCustomer->link = NULL;
			
			//also update tail if necessary
			if(foundCustomer == tail) {
				tail = tempPtr;
			}
		}
		
		//calculate the time this customer resided in the queue
		gettimeofday(&foundCustomer->end, NULL);
		long seconds = foundCustomer->end.tv_sec - foundCustomer->start.tv_sec;
		long useconds = foundCustomer->end.tv_usec - foundCustomer->start.tv_usec;
		long mtime = (long) (((seconds) * 1000 + useconds/1000.0) + 0.5);
		
		//update total queue time
		totalQtime += (double) mtime;
		
		//update our counts
		assignedCount++;
		customerCount--;
		
		return foundCustomer;
	}
}

customer *customerqueue::peak(int i) {
	customer *tempPtr = head;
	int k;
	
	//go down queue until we reach the ith item, or the end
	for(k = 0; k < i && tempPtr->link != NULL; k++) {
		tempPtr = tempPtr->link;
	}
	
	if(k != i) {
		//the ith item does not exist
		tempPtr = NULL;
	}

	return tempPtr;
}

void customerqueue::display() {
	// temporary task pointer
	customer *ptr = head;

	if(!head) {
		fprintf(stdout, "Empty Queue\n");
	} else {
		//todo: update code below to print correct statistics
		fprintf(stdout, "\n");
		fprintf(stdout, "Items in queue ******************************\n");
		fprintf(stdout, "---------------------------------------------\n");
		fprintf(stdout, "Status	Time	ID	 D_ID P_Length	toParent latitude longitude\n");
		fprintf(stdout, "---------------------------------------------\n");
		
		while(ptr != NULL) {
			fprintf(stdout, 
				"%d %d %d",
				4, ptr->time,
				4, ptr->hashLong,
				4, ptr->hashLat);
			ptr = ptr->link;
		}
	}
	fprintf(stdout, "\n");

	return;
}

#endif
