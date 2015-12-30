#ifndef __OPERATION_H
#define __OPERATION_H

#include <iostream>
#include <sstream>

#include "location.h"

using namespace std;


class Operation
{
private:

	unsigned int arrival;	
	unsigned int point;
	double quantity;
	
	
public:
	
	
	Operation();
	
	unsigned int getArrival();
	void setArrival(unsigned int a);
	unsigned int getPoint();
	void setPoint(unsigned int p);
	double getQuantity();
	void setQuantity(double q);	
};

#endif
