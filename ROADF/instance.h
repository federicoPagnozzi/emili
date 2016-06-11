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
#include <map>
#include <float.h>
#include <math.h> 
#include <limits>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <climits>

using namespace std;


class Instance
{
private:
	string name;
	int unit, horizon;	
    vector< vector<unsigned int> > timeMatrices;
	vector<Driver> drivers;
	vector<Trailer> trailers;
	vector<Customer> customers;
	vector< vector<double> > distMatrices;

//    double maxCapacity;
//    double maxInitialQuantity;

    vector< vector<double> > horizons;
    vector< vector<double> > actualQuantity;

    vector< pair< pair<unsigned int,unsigned int>, unsigned int > > timeWindows;

    map< unsigned int , vector< vector<unsigned int> >   > contributeMatrixes;
    map< unsigned int , vector< vector<double> >   > contributeMatrixesValues;

    vector< vector< vector <unsigned int> > > fastestSubroute;

    unsigned int computeTime(unsigned int succ, unsigned int nextSucc);

	
public:
	
	
	Instance();
	
    string getName();
	int getUnit();
	int getHorizon();
    vector<vector<unsigned int> > getTimeMatrices();
	vector<Driver> getDrivers();
	vector<Trailer> getTrailers();
	vector<Customer> getCustomers();
	vector< vector<double> > getDistMatrices();

    vector< vector<double> > getHorizons();
    vector< vector<double> > getActualQuantity();
    void setHorizons(vector< vector<double> > h);
    void setActualQuantity(vector< vector<double> > aq);

    vector< pair< pair<unsigned int,unsigned int>, unsigned int > > getTimeWindows();
    void setTimeWindows(vector< pair< pair<unsigned int,unsigned int>, unsigned int > > tw);

    map< unsigned int , vector< vector<unsigned int> >   > getContributeMatrixes();
    void setContributeMatrixes(map< unsigned int , vector< vector<unsigned int> >   > cm);

    map< unsigned int , vector< vector<double> >   > getContributeMatrixesValues();
    void setContributeMatrixesValues(map< unsigned int , vector< vector<double> >   > cm);

    vector< vector< vector <unsigned int> > > getFastestSubroute();
    void setFastestSubroute(vector< vector< vector <unsigned int> > > fs);


	void loadInstance(const char *pFilename);
    bool dri01(irpSolution solution);
    bool dri03(irpSolution solution);
    bool dri08(irpSolution solution);
    bool tl01(irpSolution solution);
    bool tl03(irpSolution solution);
    double dyn01(irpSolution solution, bool feasibility);
    double dyn01partial(irpSolution solution, bool feasibility, unsigned int partialHorizon);
    bool shi02(irpSolution solution);
    bool shi03(irpSolution solution);
    bool shi05(irpSolution solution);
    bool shi06(irpSolution solution);
    double checkFeasibility(irpSolution solution, double feasibility);
	
    double computeObjective(irpSolution &solution);
    double computeDeltaObjective(irpSolution oldSolution, irpSolution newSolution, unsigned int shiftIndex);
    double computePartialObjective(irpSolution oldSolution, irpSolution newSolution, unsigned int shiftIndex);
    double computePartialDeltaObjective(irpSolution oldSolution, irpSolution newSolution, unsigned int shiftIndex);

    void computeObjectiveShift(irpSolution &solution);
	
//	void sort (vector< pair<unsigned int,double> > &priorities, vector<unsigned int> urgency);
    void sort (vector<unsigned int> &priorities, vector<double> urgency);
	void sortPair(vector< pair< pair<unsigned int,unsigned int>, unsigned int > > &priorities);
    irpSolution randomSolution(double par);
	
    vector<Operation> recursiveRandomShift(irpSolution solution, unsigned int time, unsigned int cumTime, vector< vector<double> > &horizonQuantities, vector<double> &tankQuantities, vector<double> &trailerQuantities, vector<double> &pastQuantities, vector<unsigned int>  &lastOperations, vector< pair<pair<unsigned int,unsigned int>, unsigned int> > &timeWindows, vector<unsigned int> &customerList, vector<Operation> &operations, vector<double> &deliveredQuantities, vector<double> totalForecasts, double timeWeight, double quantityWeight, int ties, vector<bool> refuelFlags);
    irpSolution recursiveRandomSolution(irpSolution solution, unsigned int time, vector< vector<double> > &horizonQuantities, vector<double> &tankQuantities, vector<double> &trailerQuantities, vector<double> &pastQuantities, vector<unsigned int>  &lastOperations, vector< pair<pair<unsigned int,unsigned int>, unsigned int> > &timeWindows, double timeWeight, double quantityWeight, int ties, vector<double> &deliveredQuantities, double totalDistance, vector<double> &totalForecasts);
    irpSolution backTrackingRandomSolution(double timeWeight, double quantityWeight, double ties);


    irpSolution rebuildSolution(irpSolution &initialSolution, vector<unsigned int> solutionRepresentation, double refuelRatio, double deliveredQuantityRatio, bool originalFlag);

    void computeUrgency(vector<unsigned int> &customerList, vector<double> &urgency,
                         Shift shift, vector<double> deliveredQuantities,
                         double timeWeight, double quantityWeight, int ties,
                         unsigned int prec);
    irpSolution constructSolution(irpSolution solution, double timeWeight, double quantityWeight, int ties,
                                  double servingRatio, double refuelRatio,
                                  unsigned int noCustomer,
                                  unsigned int maxShifts);
    void computeUrgency2(vector<unsigned int> &customerList, vector<double> &urgency,
                         Shift shift, vector<double> deliveredQuantities,
                         double costWeight, double demandWeight,
                         unsigned int prec);
    irpSolution constructSolution2(irpSolution solution, double costWeight, double demandWeight, int ties,
                                  double servingRatio, double refuelRatio,
                                  unsigned int noCustomer,
                                  unsigned int maxShifts);
    irpSolution randomizedConstructSolution(irpSolution solution, double timeWeight, double quantityWeight, int ties,
                                  double servingRatio, double refuelRatio,
                                  unsigned int noCustomer,
                                  unsigned int maxShifts, double randomPick, unsigned int urgencyPolicy);

    irpSolution extendSolution(irpSolution &solution, double servingRatio, double refuelRatio,
                               unsigned int noCustomer,
                               unsigned int maxShift);

    irpSolution exchangeShiftSolution(irpSolution solution, unsigned int shiftIndex, unsigned int operationIndex, double servingRatio, double refuelRatio);
    irpSolution insertShiftSolution(irpSolution solution, unsigned int shiftIndex, unsigned int operationIndex, unsigned int insertedOperation, double servingRatio, double refuelRatio);
    irpSolution removeShiftSolution(irpSolution solution, unsigned int shiftIndex, unsigned int operationIndex, double servingRatio, double refuelRatio);

    irpSolution perturbation(irpSolution solution, unsigned int customer);
};

#endif
