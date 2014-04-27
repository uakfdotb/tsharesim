#ifndef TAXILINKLIST_H
#define TAXILINKLIST_H

struct node {
	int vehicle; //the vehicle represented by this node
	node *link; //next node in the linked list, or NULL if this is the tail
};

class taxiLinkList
{
private:
	node *p; //the head of the linked list

public:
	taxiLinkList();
	 ~taxiLinkList();
	node *returnHead();
	void append( int num ); //append a taxi to the linked list
	bool del( int num ); //delete a taxi from the linked list
	void display(); //display the contents
};

taxiLinkList::taxiLinkList()
{
	p = NULL;
}

taxiLinkList::~taxiLinkList()
{
	node *q;
	
	//delete all nodes in the linked list
	while(p != NULL) {
		q = p->link;
		delete p;
		p = q;
	}
}

node * taxiLinkList::returnHead()
{
	return p;
}

void taxiLinkList::append(int num)
{
	node *q,*t;

	if( p == NULL ) {
		//make head
		p = new node;
		p->vehicle = num;
		p->link = NULL;
	}
	else {
		//append
		q = p;
		
		while(q->link != NULL) {
			q = q->link;
		}
		
		t = new node;
		t->vehicle = num;
		t->link = NULL;
		q->link = t;
	}
}

bool taxiLinkList::del( int num )
{
	if(p == NULL) {
		return false;
	}
	
	//store both the current node and previous node
	//when current->vehicle == num, we can update the
	// previous node to point to the current's link
	node *q,*r;
	q = p;
	r = q;
	
	while(q != NULL) {
		if(q->vehicle == num) {
			//if this is the head, we need to update the head
			//otherwise we need to update previous node's link
			if(p == q) {
				p = q->link;
			} else {
				r->link = q->link;
			}
			
			delete q;
			return true;
		}

		r = q;
		q = q->link;
	}
	
	return false;
}

void taxiLinkList::display()
{
	//iterate through the list and display contents
	
	node *q;
	fprintf(stdout,"\n");
	for(q = p ; q != NULL ; q = q->link)
		fprintf(stdout,"vehicle=%d\n", q->vehicle);

}

#endif
