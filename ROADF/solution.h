#ifndef __SOLUTION_H
#define __SOLUTION_H

#include <iostream>
#include <sstream>

#include "matrix.h"
#include "driver.h"
#include "trailer.h"
#include "location.h"

#include <iostream>
#include <sstream>

#include "tinyxml.h"
#include "tinystr.h"
#include "matrix.h"
#include "driver.h"
#include "trailer.h"
#include "location.h"

#include <string>
#include <vector>
#include <utility>

#include "operation.h"

using namespace std;


class Shift
{
private:

	unsigned int index;	
	unsigned int driver;
	unsigned int trailer;
	unsigned int start;
    vector<Operation> operations;
    double shiftQuantity;
    double shiftDistance;
    double shiftTime;
    double shiftCost;
	
public:
		
	Shift();
	
	unsigned int getIndex();
	void setIndex(unsigned int i);
	unsigned int getDriver();
	void setDriver(unsigned int d);
	unsigned int getTrailer();
	void setTrailer(unsigned int t);
	unsigned int getStart();
	void setStart(unsigned int s);
	vector<Operation> getOperations();
	void setOperations(vector<Operation> operations);
    void setOperation(Operation o, int index);

    double getShiftQuantity();
    void setShiftQuantity(double sq);
    double getShiftDistance();
    void setShiftDistance(double sd);
    double getShiftTime();
    void setShiftTime(double st);

    double getShiftCost();
    void setShiftCost(double sc);
};


class irpSolution
{
private:

	vector<Shift> shifts;
    vector<unsigned int> representation;

    double cost;
    double totalQuantity;

    vector< vector <double> > trailerQuantities;
    vector< vector <double> > tankQuantities;
	
public:
		
    irpSolution();

	vector<Shift> getShifts();
    void setShifts(vector<Shift> s);
    void setShift(Shift s, int index);
    void setPoint(unsigned int point, int shift, int operation);
    vector<unsigned int> getRepresentation();
    void setRepresentation(vector<unsigned int> representation);


    double getTotalQuantity();
    void setTotalQuantity(double tq);
    double getCost();
    void setCost(double c);

    vector< vector<double> > getTrailerQuantities();
    void setTrailerQuantities(vector< vector<double> > tq);
    vector< vector<double> > getTankQuantities();
    void setTankQuantities(vector< vector<double> > tq);
	
	void loadSolution(const char *pFilename);
    void saveSolution(string pFilename);

    void fromSolutionToRepresentation(irpSolution solution);

};
#endif
