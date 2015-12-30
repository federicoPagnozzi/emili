#ifndef __INSTANCE_H
#define __INSTANCE_H

#include <iostream>
#include <sstream>
#include <fstream>

#include "tinyxml.h"
#include "tinystr.h"
#include "matrix.h"
#include "solution.h"
#include "driver.h"
#include "trailer.h"
#include "location.h"
#include "operation.h"

#include <string>
#include <vector>
#include <float.h>
#include <math.h> 
#include <limits>

using namespace std;


class Instance
{
private:
	string name;
	int unit, horizon;	
	vector< vector<double> > timeMatrices;
	vector<Driver> drivers;
	vector<Trailer> trailers;
	vector<Customer> customers;
	vector< vector<double> > distMatrices;
	
	
public:
	
	
	Instance();
	
	int getUnit();
	int getHorizon();
	vector< vector<double> > getTimeMatrices();
	vector<Driver> getDrivers();
	vector<Trailer> getTrailers();
	vector<Customer> getCustomers();
	vector< vector<double> > getDistMatrices();
	
	void loadInstance(const char *pFilename);
    bool dri01(irpSolution solution);
    bool dri03(irpSolution solution);
    bool dri08(irpSolution solution);
    bool tl01(irpSolution solution);
    bool tl03(irpSolution solution);
    double dyn01(irpSolution solution, bool feasibility);
    bool shi02(irpSolution solution);
    bool shi03(irpSolution solution);
    bool shi05(irpSolution solution);
    bool shi06(irpSolution solution);
    double checkFeasibility(irpSolution solution, double feasibility);
	
    double computeObjective(irpSolution solution);
	
//	void sort (vector< pair<unsigned int,double> > &priorities, vector<unsigned int> urgency);
    void sort (vector<unsigned int> &priorities, vector<double> urgency);
    irpSolution randomSolution2();
	void sortPair(vector< pair< pair<unsigned int,unsigned int>, unsigned int > > &priorities);
    irpSolution randomSolution(double par);
	
    vector<Operation> recursiveRandomShift(irpSolution solution, unsigned int time, unsigned int cumTime, vector< vector<double> > &horizonQuantities, vector<double> &tankQuantities, vector<double> &trailerQuantities, vector<double> &pastQuantities, vector<unsigned int>  &lastOperations, vector< pair<pair<unsigned int,unsigned int>, unsigned int> > &timeWindows, vector<unsigned int> &customerList, vector<Operation> &operations, vector<double> &deliveredQuantities, vector<double> totalForecasts, double timeWeight, double quantityWeight, int ties, vector<bool> refuelFlags);
    irpSolution recursiveRandomSolution(irpSolution solution, unsigned int time, vector< vector<double> > &horizonQuantities, vector<double> &tankQuantities, vector<double> &trailerQuantities, vector<double> &pastQuantities, vector<unsigned int>  &lastOperations, vector< pair<pair<unsigned int,unsigned int>, unsigned int> > &timeWindows, double timeWeight, double quantityWeight, int ties, vector<double> &deliveredQuantities, double totalDistance, vector<double> &totalForecasts);
    irpSolution backTrackingRandomSolution(double timeWeight, double quantityWeight, int ties);

    irpSolution GreedySolution();

//    void twoExchangeNeighborhood(irpSolution solution);
    irpSolution twoExchangeNeighborhood(irpSolution solution, int s1, int o1, int s2, int o2);

    irpSolution rebuildSolution(irpSolution &initialSolution, vector<unsigned int> solutionRepresentation, double refuelRatio, double deliveredQuantityRatio);
    void compute(irpSolution &solution);

};

#endif
