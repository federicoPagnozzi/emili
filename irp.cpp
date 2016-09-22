#include "irp.h"

#include <algorithm>
#include <ctime>
#include <limits.h>
#include <time.h>

#define EPSILON 0.000001
#define FEASIBILITY_PENALTY 1

#define COUT if (0) cout
#define CIN if (0) cin

const void* emili::irp::InventoryRoutingSolution::getRawData()const{

    return (void *) &this->irps;
}

void emili::irp::InventoryRoutingSolution::setRawData(const void* data){
    irpSolution *irps = (irpSolution *)data;
    this->irps = *irps;
}

std::string emili::irp::InventoryRoutingSolution::getSolutionRepresentation(){

    std::ostringstream repr;

    for(int s=0; s<this->irps.getShifts().size(); s++){
        repr<<"SHIFT: "<<this->irps.getShifts()[s].getIndex()<<" "<<this->irps.getShifts()[s].getStart()<<" "<<this->irps.getShifts()[s].getDriver()<<"\n";
        for(int o=0; o<this->irps.getShifts()[s].getOperations().size(); o++)
            repr<<"     "<<this->irps.getShifts()[s].getOperations()[o].getPoint()<<" "
                  <<this->irps.getShifts()[s].getOperations()[o].getArrival()<<" "
                    <<this->irps.getShifts()[s].getOperations()[o].getQuantity()<<"\n";
    }
    return repr.str();
}

emili::irp::InventoryRoutingProblem::InventoryRoutingProblem(char* instance_path){

   irpInstance.loadInstance(instance_path);
   this->time = clock();
}


double emili::irp::InventoryRoutingProblem::evaluateSolution(Solution &solution){

   InventoryRoutingSolution& s = dynamic_cast<InventoryRoutingSolution&> (solution);
   double p = this->irpInstance.computeObjective(s.getIrpSolution());
   double feas = this->irpInstance.checkFeasibility(s.getIrpSolution(), false);
   if(feas){
       double factor = feas/(this->irpInstance.getHorizon()*60);
       p /=  (factor * FEASIBILITY_PENALTY);
       p += 1.0;
   }

   solution.setSolutionValue(p);
   return p;

}


Instance emili::irp::InventoryRoutingProblem::getIrpInstance(){
    return irpInstance;
}

emili::Solution* emili::irp::GreedyInitialSolution::generateSolution(){

    InventoryRoutingProblem& irp = dynamic_cast<InventoryRoutingProblem&> (this->instance);
    InventoryRoutingSolution *irs, *bestIrs;
    irpSolution solution, irps;
    unsigned int feasibleOriginalCounter = 0;

    bool bf = false;
    double bestValue = DBL_MAX;
    for(double tw=0.1; tw<=1.0; tw+=0.1){
        for(double qw=0.0; qw<=1.0; qw+=0.1){
            for(double t=-1; t<=1; t+=2){

                irps = irp.getIrpInstance().randomizedConstructSolution(solution, tw, qw, t, 1.0, 0.0, 0, INT_MAX, this->randomFactor, this->urgencyPolicy);
                irps = irp.getIrpInstance().extendSolution(irps, 1.0, 0.0, 0, INT_MAX);
                irs = new InventoryRoutingSolution(irps);
                instance.evaluateSolution(*irs);

                if(not bf){
                    bestIrs = new InventoryRoutingSolution(irps);
                    instance.evaluateSolution(*bestIrs);
                    bf = true;
                }
                COUT<<"PARAMETERS: "<<tw<<" "<<qw<<" "<<t<<"\n";

                    if(irs->getSolutionValue() < 1.0 /*and bestValue > irs->getSolutionValue()*/)
                    {
                        COUT<<"INITIAL SOLUTION FEASIBLE! \n";

                        delete bestIrs;
                        bestIrs = new InventoryRoutingSolution(irps);

                        instance.evaluateSolution(*bestIrs);
                        bestIrs->getIrpSolution().fromSolutionToRepresentation(bestIrs->getIrpSolution());
                        bestValue = bestIrs->getSolutionValue();
/*
                        string filepath;
                        filepath.append("./Neighborhood/");
                        filepath.append(to_string(feasibleOriginalCounter));
                        filepath.append("OriginalSolution.xml");
                        irs->getIrpSolution().saveSolution(filepath);
                        feasibleOriginalCounter++;
*/
                        COUT<<irs->getSolutionRepresentation();
                        COUT<<"PARAMETERS: "<<tw<<" "<<qw<<" "<<t<<"\n";
                        COUT<<irs->getSolutionValue();
           //             int a;cin>>a;

                        ofstream file;
                        string filepath2;
                        filepath2.append("./Neighborhood/");
                        filepath2.append(irp.getIrpInstance().getName());
                        filepath2.append("/Objective");
                        file.open (filepath2,fstream::app);
                        file.precision(15);
                        file << "-------------------------------------------------" << std::endl;
                        file.close();

                        double time_elapsed = (double)(clock()-irp.getTime())/CLOCKS_PER_SEC;

                        string filepath3;
                        filepath3.append("./Neighborhood/");
                        filepath3.append(irp.getIrpInstance().getName());
                        filepath3.append("/Objective");
                        file.open (filepath3,fstream::app);
                        file.precision(15);
                        file << bestIrs->getSolutionValue() <<" "<< 0 <<" "<< 0 << " "<< time_elapsed<< std::endl;
                        file.close();

                        return bestIrs;
                    
                    }
                    delete irs;
            }
        }
    }

    ofstream file;
    string filepath2;
    filepath2.append("./Neighborhood/");
    filepath2.append(irp.getIrpInstance().getName());
    filepath2.append("/Objective");
    file.open (filepath2,fstream::app);
    file.precision(15);
    file << "-------------------------------------------------" << std::endl;
    file.close();

    double time_elapsed = (double)(clock()-irp.getTime())/CLOCKS_PER_SEC;

    string filepath3;
    filepath3.append("./Neighborhood/");
    filepath3.append(irp.getIrpInstance().getName());
    filepath3.append("/Objective");
    file.open (filepath3,fstream::app);
    file.precision(15);
    file << bestIrs->getSolutionValue() <<" "<< 0 <<" "<< 0 << " "<< time_elapsed<< std::endl;
    file.close();

    return bestIrs;

}

template <typename T>
vector<size_t> sort_indexes(const vector<T> &v) {

  // initialize original index locations
  vector<size_t> idx(v.size());
  for (size_t i = 0; i != idx.size(); ++i) idx[i] = i;

  // sort indexes based on comparing values in v
  sort(idx.begin(), idx.end(),
       [&v](size_t i1, size_t i2) {return v[i1] < v[i2];});

  return idx;
}

emili::Solution* emili::irp::GreedyRandomizedInitialSolution::generateSolution(){

    InventoryRoutingProblem& irp = dynamic_cast<InventoryRoutingProblem&> (this->instance);
    InventoryRoutingSolution *irs, *bestIrs;
    irpSolution solution, irps;
    unsigned int feasibleOriginalCounter = 0;


    bool bf = false;
    vector<InventoryRoutingSolution*> candidateSolutions;
    vector<double> objectiveCandidateSolutions;
    double bestValue = DBL_MAX;
    for(double tw=0.0; tw<=1.0; tw+=0.1){
        for(double qw=0.0; qw<=1.0; qw+=0.1){
            for(double t=-1; t<=1; t+=2){

                irps = irp.getIrpInstance().randomizedConstructSolution(solution, tw, qw, t, 1.0, 0.0, 0, INT_MAX, this->randomFactor, this->urgencyPolicy);
                irps = irp.getIrpInstance().extendSolution(irps, 1.0, 0.0, 0, INT_MAX);
                irs = new InventoryRoutingSolution(irps);
                instance.evaluateSolution(*irs);

                if(not bf){
                    bestIrs = new InventoryRoutingSolution(irps);
                    instance.evaluateSolution(*bestIrs);
                    bestValue = bestIrs->getSolutionValue();
                    bf = true;
                    candidateSolutions.push_back(bestIrs);
                    objectiveCandidateSolutions.push_back(bestValue);
                }
                COUT<<"PARAMETERS: "<<tw<<" "<<qw<<" "<<t<<"\n";

                    if(irs->getSolutionValue() < 1.0)
                    {
                        COUT<<"INITIAL SOLUTION FEASIBLE! \n";

                        bestIrs = new InventoryRoutingSolution(irps);

                        instance.evaluateSolution(*bestIrs);
                        bestValue = bestIrs->getSolutionValue();
/*
                        string filepath;
                        filepath.append("./Neighborhood/");
                        filepath.append(to_string(feasibleOriginalCounter));
                        filepath.append("OriginalSolution.xml");
                        irs->getIrpSolution().saveSolution(filepath);
                        feasibleOriginalCounter++;
*/
                        candidateSolutions.push_back(bestIrs);
                        objectiveCandidateSolutions.push_back(bestValue);



                    }

                    COUT<<irs->getSolutionRepresentation();
                    COUT<<"OBJ VALUE: "<<irs->getSolutionValue();

                    delete irs;
            }
        }
    }
/*
    cout<<"CL:\n";
    for(int i=0; i<candidateSolutions.size(); i++){
        cout<<candidateSolutions[i].second<<"\n";
    }
    vector<long unsigned int> indexes = sort_indexes(candidateSolutions);
    cout<<"CL:\n";
    for(int i=0; i<indexes.size(); i++){
        cout<<candidateSolutions[indexes[i]].second<<"\n";
    }
    */
    double total = 0.0;
    vector<double> cumulatedProbabilities;
/*
    for(int i=0; i<objectiveCandidateSolutions.size(); i++){
        objectiveCandidateSolutions[i] = 1.0/objectiveCandidateSolutions[i];
    }
*/
    for(int i=0; i<objectiveCandidateSolutions.size(); i++){
        total += objectiveCandidateSolutions[i];
    }
    vector<long unsigned int> indexes = sort_indexes(objectiveCandidateSolutions);
    COUT<<"CL:\n";
    for(int i=0; i<indexes.size(); i++){
        COUT<<objectiveCandidateSolutions[indexes[i]]<<" "<<indexes[i]<<" "<<total<< "\n";
        cumulatedProbabilities.push_back((double)objectiveCandidateSolutions[indexes[i]]/(double)total);
    }

    COUT<<"CL:\n";
    for(int i=0; i<indexes.size(); i++){
        if(i>0)
        cumulatedProbabilities[i] += cumulatedProbabilities[i-1];
        COUT<<cumulatedProbabilities[i]<<"\n";
    }

    double randomPick = generateRealRandomNumber();//(double)rand()/RAND_MAX;
    unsigned int pickIndex = 0;
    while(randomPick > cumulatedProbabilities[pickIndex]){
        pickIndex++;
        if(pickIndex >= cumulatedProbabilities.size()-1)
            break;
    }

    bestIrs = candidateSolutions[indexes[pickIndex]];

    COUT<<"BEST PICK: "<<randomPick<<" "<<pickIndex<<" "<<bestIrs->getSolutionValue()<<"\n";
    COUT<<"BEST Value: "<<bestIrs->getSolutionValue()<<"\n";
    COUT<<bestIrs->getSolutionRepresentation();

    ofstream file;
    string filepath2;
    filepath2.append("./Neighborhood/");
    filepath2.append(irp.getIrpInstance().getName());
    filepath2.append("/Objective");
    file.open (filepath2,fstream::app);
    file.precision(15);
    file << "-------------------------------------------------" << std::endl;
    file.close();

    double time_elapsed = (double)(clock()-irp.getTime())/CLOCKS_PER_SEC;

    string filepath3;
    filepath3.append("./Neighborhood/");
    filepath3.append(irp.getIrpInstance().getName());
    filepath3.append("/Objective");
    file.open (filepath3,fstream::app);
    file.precision(15);
    file << bestIrs->getSolutionValue() <<" "<< 0 <<" "<< 0 << " "<< time_elapsed<< std::endl;
    file.close();

    return bestIrs;

}


emili::Solution* emili::irp::GRASP::generateSolution(){

    InventoryRoutingProblem& irp = dynamic_cast<InventoryRoutingProblem&> (this->instance);
    InventoryRoutingSolution *irs, *bestIrs;
    unsigned int feasibleOriginalCounter = 0;

    irpSolution solution, partialSolution;
    unsigned int shiftIndex = 0;
    unsigned int currentShiftIndex = 0;
    double objValue;


    do{
        shiftIndex++;

        vector<InventoryRoutingSolution*> candidateSolutions;
        vector<double> objectiveCandidateSolutions;

        for(double tw=0.0; tw<=1.0; tw+=0.1){
            for(double qw=0.0; qw<=1.0; qw+=0.1){
                for(double t=-1; t<=1; t+=2){
                    partialSolution = irp.getIrpInstance().randomizedConstructSolution(solution, tw, qw, t, 1.0, 0.0, 0, shiftIndex, this->randomFactor, this->urgencyPolicy);
                    irs = new InventoryRoutingSolution(partialSolution);

                    if(shiftIndex > 1){
                        objValue = irp.getIrpInstance().computePartialDeltaObjective(solution, partialSolution, shiftIndex-1);
                   }
                    else
                        objValue = irp.evaluateSolution(*irs);

                    candidateSolutions.push_back(irs);
                    objectiveCandidateSolutions.push_back(objValue);

                }
            }
        }

        double total = 0.0;
        vector<double> cumulatedProbabilities;
/*
        for(int i=0; i<objectiveCandidateSolutions.size(); i++){
            objectiveCandidateSolutions[i] = 1.0/objectiveCandidateSolutions[i];
        }
*/
        for(int i=0; i<objectiveCandidateSolutions.size(); i++){
            total += objectiveCandidateSolutions[i];
        }
        vector<long unsigned int> indexes = sort_indexes(objectiveCandidateSolutions);
        for(int i=0; i<indexes.size(); i++){
            cumulatedProbabilities.push_back((double)objectiveCandidateSolutions[indexes[i]]/(double)total);
            COUT<<objectiveCandidateSolutions[indexes[i]]<<"\n";
        }
        for(int i=0; i<cumulatedProbabilities.size(); i++){
            if(i>0)
                cumulatedProbabilities[i] += cumulatedProbabilities[i-1];
                COUT<<cumulatedProbabilities[i]<<"\n";
        }

        double randomPick = generateRealRandomNumber();/*(double)rand()/RAND_MAX*/;
        unsigned int pickIndex = 0;
        while(randomPick > cumulatedProbabilities[pickIndex]){
            pickIndex++;
            if(pickIndex >= cumulatedProbabilities.size()-1)
                break;
        }


//        solution = candidateSolutions[indexes.front()]->getIrpSolution();
        solution = candidateSolutions[indexes[0/*randomPick*/]]->getIrpSolution();

        bestIrs = new InventoryRoutingSolution(solution);
        COUT<<"PARTIAL OBJ: "<<irp.evaluateSolution(*bestIrs)<<"\n";

        currentShiftIndex = solution.getShifts().size();
        COUT<<shiftIndex<<" "<<currentShiftIndex+1<<" "<<randomPick<<"\n";
        COUT<<bestIrs->getSolutionRepresentation();
    }
    while(shiftIndex < currentShiftIndex+1);


    irpSolution bestSolution = irp.getIrpInstance().extendSolution(bestIrs->getIrpSolution(), 1.0, 0.0, 0, INT_MAX);
    bestIrs = new InventoryRoutingSolution(bestSolution);
    irp.evaluateSolution(*bestIrs);

    COUT<<"PARTIAL OBJ: "<<irp.evaluateSolution(*bestIrs)<<"\n";
    COUT<<bestIrs->getSolutionRepresentation();

    ofstream file;
    string filepath2;
    filepath2.append("./Neighborhood/");
    filepath2.append(irp.getIrpInstance().getName());
    filepath2.append("/Objective");
    file.open (filepath2,fstream::app);
    file.precision(15);
    file << "-------------------------------------------------" << std::endl;
    file.close();

    double time_elapsed = (double)(clock()-irp.getTime())/CLOCKS_PER_SEC;

    string filepath3;
    filepath3.append("./Neighborhood/");
    filepath3.append(irp.getIrpInstance().getName());
    filepath3.append("/Objective");
    file.open (filepath3,fstream::app);
    file.precision(15);
    file << bestIrs->getSolutionValue() <<" "<< 0 <<" "<< 0 << " "<< time_elapsed<< std::endl;
    file.close();

    return bestIrs;

}

emili::Solution* emili::irp::InventoryRoutingSolution::clone(){

    InventoryRoutingSolution *irs = new InventoryRoutingSolution(this->irps);
    irs->setSolutionValue(this->solution_value);
    return irs;
}


irpSolution & emili::irp::InventoryRoutingSolution::getIrpSolution(){
    return this->irps;
}

emili::Solution* emili::irp::GreedyInitialSolution::generateEmptySolution(){

    return new InventoryRoutingSolution(DBL_MAX);
}

emili::Solution* emili::irp::GreedyRandomizedInitialSolution::generateEmptySolution(){

    return new InventoryRoutingSolution(DBL_MAX);
}

emili::Solution* emili::irp::GRASP::generateEmptySolution(){

    return new InventoryRoutingSolution(DBL_MAX);
}



emili::Neighborhood::NeighborhoodIterator emili::irp::irpShiftTwoExchangeNeighborhood::begin(emili::Solution *startSolution)
{

    this->currentNeighboringSolution = startSolution->clone();
    InventoryRoutingSolution *irpStartSolution = dynamic_cast<InventoryRoutingSolution *> (startSolution);

    this->shift = 0;
    this ->operation = -1;
    this->numberOfShifts = irpStartSolution->getIrpSolution().getShifts().size();
    this->numberOfOperations = irpStartSolution->getIrpSolution().getShifts()[shift].getOperations().size();


    if(this->bestValueFound >= DBL_MAX - EPSILON){
        this->bestValueFound = irpStartSolution->getSolutionValue();
 //       this->numberFeasibleSolutions++;
    }
    COUT<<"\nINITIAL BEGIN: \n";

    return emili::Neighborhood::NeighborhoodIterator(this,startSolution);

}



emili::Solution* emili::irp::irpShiftTwoExchangeNeighborhood::step(Solution* currentSolution){

    return this->computeStep(currentSolution);
}

emili::Solution* emili::irp::irpShiftTwoExchangeNeighborhood::computeStep(Solution* currentSolution){

    emili::iteration_increment();

    InventoryRoutingSolution *neighboringSolution = dynamic_cast<InventoryRoutingSolution *> (currentSolution);

    this->numberOfOperations = neighboringSolution->getIrpSolution().getShifts()[this->shift].getOperations().size();
    this->numberOfShifts = neighboringSolution->getIrpSolution().getShifts().size();
    if(this->operation < numberOfOperations - 2 and this->shift < this->numberOfShifts){
        this->operation++;
    }
    else if(this->shift < this->numberOfShifts-1){
        this->shift++;
        this->operation = 0;

    }
    else
        return nullptr;

    COUT<<"\nTWO SHIFT EXCHANGE: "<<this->shift<<" "<<this->operation<<"\n";

    irpSolution irps = neighboringSolution->getIrpSolution();
    irps = irp.getIrpInstance().exchangeShiftSolution(irps, shift, operation, this->deliveredQuantityFactor, this->refuelFactor);

    vector<Shift> shifts = irps.getShifts();
    while(shifts.size() > shift+1)
        shifts.erase(shifts.end());

    irps.setShifts(shifts);
//    irps = irp.getIrpInstance().constructSolution(irps, 0.1, 0.0, 1, 1.0, 0.0, 0, INT_MAX);
    irps = irp.getIrpInstance().randomizedConstructSolution(irps, 0.1, 0.0, 1, this->deliveredQuantityFactor, this->refuelFactor, 0, INT_MAX, this->randomFactor, this->urgencyPolicy);
    irps = irp.getIrpInstance().extendSolution(irps, this->deliveredQuantityFactor, this->refuelFactor, 0, INT_MAX);


    InventoryRoutingSolution irs(irps);
    double newObjValue =  irp.getIrpInstance().computePartialDeltaObjective(neighboringSolution->getIrpSolution(), irs.getIrpSolution(), shift);
    COUT<<"OBJECTIVE: "<<newObjValue<<"\n";
    irs.setSolutionValue(newObjValue);


    if(irs.getSolutionValue() < this->bestValueFound){
        this->numberFeasibleSolutions++;
        this->bestValueFound = irs.getSolutionValue();
        string filepath;

        double time_elapsed = (double)(clock()-irp.getTime())/CLOCKS_PER_SEC;

       /*
       filepath.append("./Neighborhood/");
       filepath.append(this->irp.getIrpInstance().getName());filepath.append("/");
       filepath.append(to_string(this->numberFeasibleSolutions));
       filepath.append("NeighSolution.xml");
       irs.getIrpSolution().saveSolution(filepath);
        */

       ofstream file;
       string filepath2;
       filepath2.append("./Neighborhood/");
       filepath2.append(this->irp.getIrpInstance().getName());
       filepath2.append("/Objective");
       file.open (filepath2,fstream::app);
       file.precision(15);
       file << this->bestValueFound <<" "<< emili::iteration_counter()<<" "<< this->numberFeasibleSolutions << " "<< time_elapsed<< std::endl;
       file.close();


       COUT<<"A BEST FOUND: "<<this->bestValueFound<<"\n";
    }
    else
        this->numberFeasibleSolutions++;

    int a; CIN>>a;
    /*
    * QUI COPIA lo stato interno di irs in currentSolution
    * neighboringSolution e currentSolution puntano allo stesso oggetto
    */
//int b;CIN>>b;
    *neighboringSolution = irs;
    return dynamic_cast<Solution *> (currentSolution);

}


void emili::irp::irpShiftTwoExchangeNeighborhood::reset(){

    this->shift = 0;
    this->operation = -1;
    this->numberOfShifts = 0;
    this->numberOfOperations = 1;
}

emili::Solution* emili::irp::irpShiftTwoExchangeNeighborhood::random(Solution* currentSolution){

    COUT<<"\nSHIFT TWO EXC RANDOM\n";
    InventoryRoutingSolution *neighboringSolution = new InventoryRoutingSolution(*dynamic_cast<InventoryRoutingSolution *> (currentSolution));
    InventoryRoutingSolution *irs;

    this->numberOfShifts = neighboringSolution->getIrpSolution().getShifts().size();
    unsigned int randomShift = generateRandomNumber() % this->numberOfShifts;
    this->numberOfOperations = neighboringSolution->getIrpSolution().getShifts()[randomShift].getOperations().size();
    unsigned int randomOperation = generateRandomNumber() % this->numberOfOperations;

    COUT<<randomShift<<" "<<randomOperation<<" "<<neighboringSolution->getSolutionValue()<<"\n";

    irpSolution irps = neighboringSolution->getIrpSolution();

    irps = irp.getIrpInstance().exchangeShiftSolution(irps, randomShift, randomOperation, this->deliveredQuantityFactor, this->refuelFactor);

    vector<Shift> shifts = irps.getShifts();
    while(shifts.size() > randomShift+1)
        shifts.erase(shifts.end());

    irps.setShifts(shifts);
//    irps = irp.getIrpInstance().constructSolution(irps, 0.1, 0.0, 1, 1.0, 0.0, 0, INT_MAX);
    irps = irp.getIrpInstance().randomizedConstructSolution(irps, 0.1, 0.0, 1, this->deliveredQuantityFactor, this->refuelFactor, 0, INT_MAX, this->randomFactor, this->urgencyPolicy);
    irps = irp.getIrpInstance().extendSolution(irps, this->deliveredQuantityFactor, this->refuelFactor, 0, INT_MAX);

    irs = new InventoryRoutingSolution(irps);
    this->irp.evaluateSolution(*irs);

    /////
//    irs = neighboringSolution;
    /////

    delete neighboringSolution;
    return dynamic_cast<Solution *> (irs);
}

void emili::irp::irpShiftTwoExchangeNeighborhood::reverseLastMove(Solution * step){

    ///fa una copia
    *step = *currentNeighboringSolution;
    this->irp.evaluateSolution(*step);

}

int emili::irp::irpShiftTwoExchangeNeighborhood::size()
{
    return 1;
}



emili::Neighborhood::NeighborhoodIterator emili::irp::irpShiftInsertNeighborhood::begin(emili::Solution *startSolution)
{

    this->currentNeighboringSolution = startSolution->clone();
    InventoryRoutingSolution *irpStartSolution = dynamic_cast<InventoryRoutingSolution *> (startSolution);

    this->shift = 0;
    this ->operation = -1;
    this->numberOfShifts = irpStartSolution->getIrpSolution().getShifts().size();
    this->numberOfOperations = irpStartSolution->getIrpSolution().getShifts()[shift].getOperations().size()+1;
    if(this->bestValueFound >= DBL_MAX - EPSILON){
        this->bestValueFound = irpStartSolution->getSolutionValue();
  //      this->numberFeasibleSolutions++;
    }
    COUT<<"\nINITIAL BEGIN: \n";

    return emili::Neighborhood::NeighborhoodIterator(this,startSolution);

}



emili::Solution* emili::irp::irpShiftInsertNeighborhood::step(Solution* currentSolution){

    return this->computeStep(currentSolution);
}

emili::Solution* emili::irp::irpShiftInsertNeighborhood::computeStep(Solution* currentSolution){

    emili::iteration_increment();

    InventoryRoutingSolution *neighboringSolution = dynamic_cast<InventoryRoutingSolution *> (currentSolution);

    this->numberOfOperations = neighboringSolution->getIrpSolution().getShifts()[this->shift].getOperations().size()+1;
    this->numberOfShifts = neighboringSolution->getIrpSolution().getShifts().size();
    if(this->operation < numberOfOperations - 1 and this->shift < this->numberOfShifts){
        this->operation++;
    }
    else if(this->shift < this->numberOfShifts-1){
        this->shift++;
        this->operation = 0;

    }
    else
        return nullptr;

    COUT<<"\nSHIFT INSERT: "<<this->shift<<" "<<this->operation<<"\n";

    irpSolution irps = neighboringSolution->getIrpSolution();
    irps = irp.getIrpInstance().insertShiftSolution(irps, shift, operation, 2, this->deliveredQuantityFactor, this->refuelFactor);

    vector<Shift> shifts = irps.getShifts();
    while(shifts.size() > shift+1)
        shifts.erase(shifts.end());

    irps.setShifts(shifts);
//    irps = irp.getIrpInstance().constructSolution(irps, 0.1, 0.0, 1, 1.0, 0.0, 0, INT_MAX);
    irps = irp.getIrpInstance().randomizedConstructSolution(irps, 0.1, 0.0, 1, this->deliveredQuantityFactor, this->refuelFactor, 0, INT_MAX, this->randomFactor, this->urgencyPolicy);
    irps = irp.getIrpInstance().extendSolution(irps, this->deliveredQuantityFactor, this->refuelFactor, 0, INT_MAX);

    InventoryRoutingSolution irs(irps);
    double newObjValue =  irp.getIrpInstance().computePartialDeltaObjective(neighboringSolution->getIrpSolution(), irs.getIrpSolution(), shift);
    irs.setSolutionValue(newObjValue);


    if(irs.getSolutionValue() < this->bestValueFound){
        this->numberFeasibleSolutions++;
        this->bestValueFound = irs.getSolutionValue();
        string filepath;

        double time_elapsed = (double)(clock()-irp.getTime())/CLOCKS_PER_SEC;

/*
       filepath.append("./Neighborhood/");
       filepath.append(this->irp.getIrpInstance().getName());filepath.append("/");
       filepath.append(to_string(this->numberFeasibleSolutions));
       filepath.append("NeighSolution.xml");
       irs.getIrpSolution().saveSolution(filepath);
*/

       ofstream file;
       string filepath2;
       filepath2.append("./Neighborhood/");
       filepath2.append(this->irp.getIrpInstance().getName());
       filepath2.append("/Objective");
       file.open (filepath2,fstream::app);
       file.precision(15);
       file << this->bestValueFound <<" "<< emili::iteration_counter()<<" "<< this->numberFeasibleSolutions << " "<< time_elapsed<< std::endl;
       file.close();

       COUT<<"A BEST FOUND: "<<this->bestValueFound<<"\n";
    }
    else
        this->numberFeasibleSolutions++;

    /*
    * QUI COPIA lo stato interno di irs in currentSolution
    * neighboringSolution e currentSolution puntano allo stesso oggetto
    */

    *neighboringSolution = irs;
    return dynamic_cast<Solution *> (currentSolution);

}


void emili::irp::irpShiftInsertNeighborhood::reset(){

    this->shift = 0;
    this->operation = -1;
    this->numberOfShifts = 0;
    this->numberOfOperations = 1;
}

emili::Solution* emili::irp::irpShiftInsertNeighborhood::random(Solution* currentSolution){

    COUT<<"\nSHIFT TWO EXC RANDOM\n";
    InventoryRoutingSolution *neighboringSolution = new InventoryRoutingSolution(*dynamic_cast<InventoryRoutingSolution *> (currentSolution));
    InventoryRoutingSolution *irs;

    this->numberOfShifts = neighboringSolution->getIrpSolution().getShifts().size();
    unsigned int randomShift = generateRandomNumber() % this->numberOfShifts;
    this->numberOfOperations = neighboringSolution->getIrpSolution().getShifts()[randomShift].getOperations().size()+1;
    unsigned int randomOperation = generateRandomNumber() % this->numberOfOperations;
    unsigned int randomInsertedOperation = generateRandomNumber() % irp.getIrpInstance().getCustomers().size();

    COUT<<randomShift<<" "<<randomOperation<<" "<<neighboringSolution->getSolutionValue()<<"\n";

    irpSolution irps = neighboringSolution->getIrpSolution();

    irps = irp.getIrpInstance().insertShiftSolution(irps, randomShift, randomOperation, randomInsertedOperation, this->deliveredQuantityFactor, this->refuelFactor);

    vector<Shift> shifts = irps.getShifts();
    while(shifts.size() > randomShift+1)
        shifts.erase(shifts.end());

    irps.setShifts(shifts);
//    irps = irp.getIrpInstance().constructSolution2(irps, 0.1, 0.0, 1, 1.0, 0.0, 0, INT_MAX);
    irps = irp.getIrpInstance().randomizedConstructSolution(irps, 0.1, 0.0, 1, this->deliveredQuantityFactor, this->refuelFactor, 0, INT_MAX, this->randomFactor, this->urgencyPolicy);
    irps = irp.getIrpInstance().extendSolution(irps, this->deliveredQuantityFactor, this->refuelFactor, 0, INT_MAX);

    irs = new InventoryRoutingSolution(irps);
    this->irp.evaluateSolution(*irs);

    delete neighboringSolution;

    return dynamic_cast<Solution *> (irs);
}

void emili::irp::irpShiftInsertNeighborhood::reverseLastMove(Solution * step){

    ///fa una copia
    *step = *currentNeighboringSolution;
    this->irp.evaluateSolution(*step);

}

int emili::irp::irpShiftInsertNeighborhood::size()
{
    return 1;
}


///////////////////////////////////////////////////////////




emili::Neighborhood::NeighborhoodIterator emili::irp::irpShiftRemoveNeighborhood::begin(emili::Solution *startSolution)
{

    this->currentNeighboringSolution = startSolution->clone();
    InventoryRoutingSolution *irpStartSolution = dynamic_cast<InventoryRoutingSolution *> (startSolution);

    this->shift = 0;
    this ->operation = -1;
    this->numberOfShifts = irpStartSolution->getIrpSolution().getShifts().size();
    this->numberOfOperations = irpStartSolution->getIrpSolution().getShifts()[shift].getOperations().size();

    if(this->bestValueFound >= DBL_MAX - EPSILON){
        this->bestValueFound = irpStartSolution->getSolutionValue();
//        this->numberFeasibleSolutions++;
    }
    COUT<<"\nINITIAL BEGIN: \n";

    return emili::Neighborhood::NeighborhoodIterator(this,startSolution);

}



emili::Solution* emili::irp::irpShiftRemoveNeighborhood::step(Solution* currentSolution){

    return this->computeStep(currentSolution);
}

emili::Solution* emili::irp::irpShiftRemoveNeighborhood::computeStep(Solution* currentSolution){

    emili::iteration_increment();

    InventoryRoutingSolution *neighboringSolution = dynamic_cast<InventoryRoutingSolution *> (currentSolution);

    this->numberOfOperations = neighboringSolution->getIrpSolution().getShifts()[this->shift].getOperations().size();
    this->numberOfShifts = neighboringSolution->getIrpSolution().getShifts().size();
    if(this->operation < numberOfOperations - 1 and this->shift < this->numberOfShifts){
        this->operation++;
    }
    else if(this->shift < this->numberOfShifts-1){
        this->shift++;
        this->operation = 0;

    }
    else
        return nullptr;

    COUT<<"\nSHIFT INSERT: "<<this->shift<<" "<<this->operation<<"\n";

    irpSolution irps = neighboringSolution->getIrpSolution();
    irps = irp.getIrpInstance().removeShiftSolution(irps, shift, operation, this->deliveredQuantityFactor, this->refuelFactor);

    vector<Shift> shifts = irps.getShifts();
    while(shifts.size() > shift+1)
        shifts.erase(shifts.end());

    irps.setShifts(shifts);
//    irps = irp.getIrpInstance().constructSolution(irps, 0.1, 0.0, 1, 1.0, 0.0, 0, INT_MAX);
    irps = irp.getIrpInstance().randomizedConstructSolution(irps, 0.1, 0.0, 1, this->deliveredQuantityFactor, this->refuelFactor, 0, INT_MAX, this->randomFactor, this->urgencyPolicy);
    irps = irp.getIrpInstance().extendSolution(irps, this->deliveredQuantityFactor, this->refuelFactor, 0, INT_MAX);

    InventoryRoutingSolution irs(irps);
    double newObjValue =  irp.getIrpInstance().computePartialDeltaObjective(neighboringSolution->getIrpSolution(), irs.getIrpSolution(), shift);
    irs.setSolutionValue(newObjValue);

    if(irs.getSolutionValue() < this->bestValueFound){
        this->numberFeasibleSolutions++;
        this->bestValueFound = irs.getSolutionValue();
        string filepath;

        double time_elapsed = (double)(clock()-irp.getTime())/CLOCKS_PER_SEC;
/*
       filepath.append("./Neighborhood/");
       filepath.append(this->irp.getIrpInstance().getName());filepath.append("/");
       filepath.append(to_string(this->numberFeasibleSolutions));
       filepath.append("NeighSolution.xml");
       irs.getIrpSolution().saveSolution(filepath);
*/

       ofstream file;
       string filepath2;
       filepath2.append("./Neighborhood/");
       filepath2.append(this->irp.getIrpInstance().getName());
       filepath2.append("/Objective");
       file.open (filepath2,fstream::app);
       file.precision(15);
       file << this->bestValueFound <<" "<< emili::iteration_counter()<<" "<< this->numberFeasibleSolutions << " "<< time_elapsed<< std::endl;
       file.close();

       COUT<<"A BEST FOUND: "<<this->bestValueFound<<"\n";
    }
    else
        this->numberFeasibleSolutions++;

    /*
    * QUI COPIA lo stato interno di irs in currentSolution
    * neighboringSolution e currentSolution puntano allo stesso oggetto
    */
    *neighboringSolution = irs;
    return dynamic_cast<Solution *> (currentSolution);

}


void emili::irp::irpShiftRemoveNeighborhood::reset(){

    this->shift = 0;
    this->operation = -1;
    this->numberOfShifts = 0;
    this->numberOfOperations = 1;
}

emili::Solution* emili::irp::irpShiftRemoveNeighborhood::random(Solution* currentSolution){

    COUT<<"\nSHIFT TWO EXC RANDOM\n";
    InventoryRoutingSolution *neighboringSolution = new InventoryRoutingSolution(*dynamic_cast<InventoryRoutingSolution *> (currentSolution));
    InventoryRoutingSolution *irs;

    this->numberOfShifts = neighboringSolution->getIrpSolution().getShifts().size();
    unsigned int randomShift = generateRandomNumber() % this->numberOfShifts;
    this->numberOfOperations = neighboringSolution->getIrpSolution().getShifts()[randomShift].getOperations().size();
    unsigned int randomOperation = generateRandomNumber() % this->numberOfOperations;

    COUT<<randomShift<<" "<<randomOperation<<" "<<neighboringSolution->getSolutionValue()<<"\n";

    irpSolution irps = neighboringSolution->getIrpSolution();

    irps = irp.getIrpInstance().removeShiftSolution(irps, randomShift, randomOperation, this->deliveredQuantityFactor, this->refuelFactor);

    vector<Shift> shifts = irps.getShifts();
    while(shifts.size() > randomShift+1)
        shifts.erase(shifts.end());

    irps.setShifts(shifts);
//    irps = irp.getIrpInstance().constructSolution(irps, 0.1, 0.0, 1, 1.0, 0.0, 0, INT_MAX);
    irps = irp.getIrpInstance().randomizedConstructSolution(irps, 0.1, 0.0, 1, this->deliveredQuantityFactor, this->refuelFactor, 0, INT_MAX, this->randomFactor, this->urgencyPolicy);
    irps = irp.getIrpInstance().extendSolution(irps, this->deliveredQuantityFactor, this->refuelFactor, 0, INT_MAX);

    irs = new InventoryRoutingSolution(irps);
    this->irp.evaluateSolution(*irs);

    delete neighboringSolution;

    return dynamic_cast<Solution *> (irs);
}

void emili::irp::irpShiftRemoveNeighborhood::reverseLastMove(Solution * step){

    ///fa una copia
    *step = *currentNeighboringSolution;
    this->irp.evaluateSolution(*step);

}

int emili::irp::irpShiftRemoveNeighborhood::size()
{
    return 1;
}

/////////////////////////////////////////






























emili::Neighborhood::NeighborhoodIterator emili::irp::irpRefuelNeighborhood::begin(emili::Solution *startSolution)
{

    /////Salvo la soluzione iniziale per ripristinarla in reverse last move
    this->currentNeighboringSolution = startSolution->clone();
    InventoryRoutingSolution *irpStartSolution = dynamic_cast<InventoryRoutingSolution *> (startSolution);

    this->refuelRatio = 0.0;
    this->deliveredQuantityRatio = -this->deliveredQuantityStep;
    this->shift = 0;
    this->numberOfShifts = irpStartSolution->getIrpSolution().getShifts().size();

    if(this->bestValueFound >= DBL_MAX - EPSILON){
        this->bestValueFound = irpStartSolution->getSolutionValue();
//        this->numberFeasibleSolutions++;
    }


    COUT<<"\nBEGIN INITIAL: "<<this->refuelRatio<<" "<<this->refuelStep<<" \n";

//    startSolution = dynamic_cast<Solution *> (irpStartSolution);

    return emili::Neighborhood::NeighborhoodIterator(this,startSolution);

}

emili::Solution* emili::irp::irpRefuelNeighborhood::step(Solution* currentSolution){

    return this->computeStep(currentSolution);
}

emili::Solution* emili::irp::irpRefuelNeighborhood::computeStep(Solution* currentSolution){

    emili::iteration_increment();

    InventoryRoutingSolution *neighboringSolution = dynamic_cast<InventoryRoutingSolution *> (currentSolution);

    this->numberOfShifts = neighboringSolution->getIrpSolution().getShifts().size();

    if(this->deliveredQuantityRatio <= 1.0-this->deliveredQuantityStep and this->refuelRatio <= 1.0-this->refuelStep + EPSILON and this->shift < this->numberOfShifts){
        this->deliveredQuantityRatio += this->deliveredQuantityStep;
    }
    else if(this->refuelRatio <= 1.0-this->refuelStep and this->shift < this->numberOfShifts){
        this->deliveredQuantityRatio = 0.0;
        this->refuelRatio += this->refuelStep;
    }
    else if(this->shift < this->numberOfShifts){
        this->shift++;
        this->deliveredQuantityRatio = 0.0;
        this->refuelRatio = 0.0;
    }
    else
        return nullptr;

    COUT<<"\nPERTURB RATIO: "<<this->deliveredQuantityRatio<<" "<<this->refuelRatio<<" "<<this->shift;
    COUT<<"     OLD: "<<neighboringSolution->getSolutionValue();


    irpSolution irps = neighboringSolution->getIrpSolution();

    vector<Shift> shifts = irps.getShifts();
    while(shifts.size() > shift)
        shifts.erase(shifts.end());

    irps.setShifts(shifts);
//    irps = irp.getIrpInstance().constructSolution(irps, 0.1, 0.0, 1, 1.0 - this->deliveredQuantityRatio, this->refuelRatio, 0, INT_MAX);
    irps = irp.getIrpInstance().randomizedConstructSolution(irps, 0.1, 0.0, 1, 1.0 - this->deliveredQuantityRatio, this->refuelRatio, 0, INT_MAX, this->randomFactor, this->urgencyPolicy);
    irps = irp.getIrpInstance().extendSolution(irps, 1.0 - this->deliveredQuantityRatio, this->refuelRatio, 0, INT_MAX);

    InventoryRoutingSolution irs(irps);

    double newObjValue =  irp.getIrpInstance().computePartialDeltaObjective(neighboringSolution->getIrpSolution(), irs.getIrpSolution(), shift);
    irs.setSolutionValue(newObjValue);

    if(irs.getSolutionValue() < this->bestValueFound - EPSILON){
      this->bestValueFound = irs.getSolutionValue();
      this->numberFeasibleSolutions++;

        double time_elapsed = (double)(clock()-irp.getTime())/CLOCKS_PER_SEC;

        string filepath;
/*
        filepath.append("./Neighborhood/");
        filepath.append(this->irp.getIrpInstance().getName());filepath.append("/");
        filepath.append(to_string(this->numberFeasibleSolutions));
        filepath.append("NeighSolution.xml");
        irs.getIrpSolution().saveSolution(filepath);
*/

        ofstream file;
        string filepath2;
        filepath2.append("./Neighborhood/");
        filepath2.append(this->irp.getIrpInstance().getName());filepath.append("/");
        filepath2.append("/Objective");
        file.open (filepath2,fstream::app);
        file.precision(15);
        file << this->bestValueFound  <<" "<< emili::iteration_counter()<<" "<< this->numberFeasibleSolutions << " "<< time_elapsed << std::endl;
        file.close();


        COUT<<"A BEST FOUND: "<<this->bestValueFound<<"\n";
   }
    /*
	* QUI COPIA lo stato interno di irs in currentSolution
	* neighboringSolution e currentSolution puntano allo stesso oggetto
	*/
	*neighboringSolution = irs;



    return currentSolution;
}


void emili::irp::irpRefuelNeighborhood::reset(){

    this->refuelRatio = 0.0;
    this->deliveredQuantityRatio = 0.0;
    this->shift = 0;
}

emili::Solution* emili::irp::irpRefuelNeighborhood::random(Solution* currentSolution){

    COUT<<"\nREF RANDOM\n";

    InventoryRoutingSolution *neighboringSolution = new InventoryRoutingSolution(*dynamic_cast<InventoryRoutingSolution *> (currentSolution));
    InventoryRoutingSolution *irs;

    this->numberOfShifts = neighboringSolution->getIrpSolution().getShifts().size();

    double randomRefuel = generateRealRandomNumber();
    double randomDeliveredQuantity = generateRealRandomNumber();
    unsigned int randomShift = generateRandomNumber() % this->numberOfShifts;



    COUT<<"\nPERTURB RATIO: "<<randomDeliveredQuantity<<" "<<randomRefuel<<" "<<randomShift;

    irpSolution irps = neighboringSolution->getIrpSolution();

    vector<Shift> shifts = irps.getShifts();
    while(shifts.size() > randomShift)
        shifts.erase(shifts.end());

    irps.setShifts(shifts);
//    irps = irp.getIrpInstance().constructSolution(irps, 0.1, 0.0, 1, 1.0, 0.0, 0, INT_MAX);
    irps = irp.getIrpInstance().constructSolution(irps, 0.1, 0.0, 1, 1.0 - randomDeliveredQuantity, randomRefuel, 0, INT_MAX);
    irps = irp.getIrpInstance().extendSolution(irps, 1.0 - randomDeliveredQuantity, randomRefuel, 0, INT_MAX);

    irs = new InventoryRoutingSolution(irps);

    this->irp.evaluateSolution(*irs);

    delete neighboringSolution;

    return dynamic_cast<Solution *> (irs);
}

void emili::irp::irpRefuelNeighborhood::reverseLastMove(Solution * step){

    //Ripristino la soluzione iniziale
    *step = *this->currentNeighboringSolution;
    this->irp.evaluateSolution(*step);

}

int emili::irp::irpRefuelNeighborhood::size()
{
    return ( ((1.0)/this->refuelStep) + 1) * ( ((1.0)/this->deliveredQuantityStep) + 1);
}

emili::Solution* emili::irp::irpPerturbation::perturb(Solution* solution){

    InventoryRoutingSolution *perturbedSolution = new InventoryRoutingSolution(*dynamic_cast<InventoryRoutingSolution *> (solution));
    irpSolution irpPerturbedSolution = perturbedSolution->getIrpSolution();
}
