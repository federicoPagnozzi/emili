#include "irp.h"

#define EPSILON 0.000001
#define FEASIBILITY_PENALTY 1
#define POINT_INITIAL_VALUE 0
#define REFUEL_INITIAL_VALUE 0.0
#define REFUEL_STEP 0.1

const void* emili::irp::InventoryRoutingSolution::getRawData()const{

    return (void *) &this->irps;
}

void emili::irp::InventoryRoutingSolution::setRawData(const void* data){
    irpSolution *irps = (irpSolution *)data;
    this->irps = *irps;
}

std::string emili::irp::InventoryRoutingSolution::getSolutionRepresentation(){

    std::ostringstream repr;
    repr<<"\nSOLUTION REPRESENTATION: ";
    for(int s=0; s<this->irps.getShifts().size(); s++)
        for(int o=0; o<this->irps.getShifts()[s].getOperations().size(); o++)
//            if(this->irps.getShifts()[s].getOperations()[o].getPoint()!=1)
            repr <<this->irps.getShifts()[s].getOperations()[o].getPoint()<<" ";
    repr<<"\n";
    return repr.str();
}

emili::irp::InventoryRoutingProblem::InventoryRoutingProblem(char* instance_path){

   irpInstance.loadInstance(instance_path);
}


double emili::irp::InventoryRoutingProblem::evaluateSolution(Solution & solution){

   InventoryRoutingSolution& s = dynamic_cast<InventoryRoutingSolution&> (solution);
   double p = this->irpInstance.computeObjective(s.getIrpSolution());
   double feas = this->irpInstance.checkFeasibility(s.getIrpSolution(), false);
//   cout<<"FEAS: "<<feas<<" "; int a;cin>>a;
   if(feas){
       double factor = feas/(this->irpInstance.getHorizon()*60 * (this->irpInstance.getCustomers().size()-2));
//       cout<<"PENALTY: "<<p<<" "<<feas<<" "<<factor<<" ";
       p /=  (factor * FEASIBILITY_PENALTY);
//       cout<<p<<"\n";
//       int a;cin>>a;
  //     cout<<"NOT FEASIBLE! "<<p<<"\n";
   }

   solution.setSolutionValue(p);
   return p;

}

Instance emili::irp::InventoryRoutingProblem::getIrpInstance(){
    return irpInstance;
}

emili::Solution* emili::irp::GreedyInitialSolution::generateSolution(){

    int a;cin>>a;
    InventoryRoutingProblem& irp = dynamic_cast<InventoryRoutingProblem&> (this->instance);
    InventoryRoutingSolution *irs, *bestIrs;
    irpSolution irps;
    unsigned int feasibleOriginalCounter = 0;

    bool bf = false;
    double bestValue = DBL_MAX;
    for(double tw=0.0; tw<=1.0; tw+=0.1){
        for(double qw=0.0; qw<=1.0; qw+=0.1){
            for(double t=-1; t<=1; t+=1){
                irps = irp.getIrpInstance().backTrackingRandomSolution(tw, qw, t);
                irs = new InventoryRoutingSolution(irps);
                instance.evaluateSolution(*irs);

                cout<<"PARAMETERS: "<<tw<<" "<<qw<<" "<<t;
                if(irp.getIrpInstance().checkFeasibility(irs->getIrpSolution(), false)){
                    cout<<"\nORIGINAL NOT FEASIBLE!\n";
            //        irs->getIrpSolution().saveSolution("OriginalSolution.xml");
                }
                else if(irs->getSolutionValue() < bestValue){
                    cout<<"\nORIGINAL FEASIBLE!\n";
                    bestIrs = new InventoryRoutingSolution(irps);
                    instance.evaluateSolution(*bestIrs);
                    bestIrs->getIrpSolution().fromSolutionToRepresentation(bestIrs->getIrpSolution());
                    bestValue = bestIrs->getSolutionValue();
                    string filepath;
                    filepath.append("./Neighborhood/");
                    filepath.append(to_string(feasibleOriginalCounter));
                    filepath.append("OriginalSolution.xml");
                    irs->getIrpSolution().saveSolution(filepath);
                    feasibleOriginalCounter++;

                        cout<<bestIrs->getSolutionRepresentation();
                        cout<<"    "<<irps.getShifts()[0].getOperations()[0].getPoint();
                    cout<<"\nORIGINAL OBJ VALUE: "<<irs->getSolutionValue()<<"\n";
                    cout<<"\nVALUES: "<<tw<<" "<<qw<<" "<<t<<"\n";
               //     bf = true;break;
                }

                if(bf)break;
            }
            if(bf)break;
        }
        if(bf)break;
    }


    vector<unsigned int> representation;

    InventoryRoutingSolution *rirs;
    irpSolution rebuiltSolution = irp.getIrpInstance().rebuildSolution(bestIrs->getIrpSolution(),bestIrs->getIrpSolution().getRepresentation(), 0.0, 1.0);
    rirs = new InventoryRoutingSolution(rebuiltSolution);
    rirs->getIrpSolution().fromSolutionToRepresentation(rirs->getIrpSolution());
    instance.evaluateSolution(*rirs);
    cout<<rirs->getSolutionRepresentation();
    if(irp.getIrpInstance().checkFeasibility(rirs->getIrpSolution(), false)){
        cout<<"\nRECOSTRUCTION NOT FEASIBLE!\n";
        rirs->getIrpSolution().saveSolution(*new string("RebuiltSolution.xml"));
        cout<<"OBJ VALUE: "<<rirs->getSolutionValue()<<"\n\n";
    }
    else{
        cout<<"\nRECOSTRUCTION FEASIBLE!\n";
        cout<<"OBJ VALUE: "<<rirs->getSolutionValue()<<"\n\n";
    }

    int b;cin>>b;
    return rirs;

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

emili::Neighborhood::NeighborhoodIterator emili::irp::irpTwoExchangeNeighborhood::begin(emili::Solution *startSolution)
{
    this->currentNeighboringSolution = startSolution->clone();
    InventoryRoutingSolution *irpStartSolution = dynamic_cast<InventoryRoutingSolution *> (startSolution);

    this->operation1 = POINT_INITIAL_VALUE;
    this->operation2 = POINT_INITIAL_VALUE;

    irpStartSolution->getIrpSolution().fromSolutionToRepresentation(irpStartSolution->getIrpSolution());
    irp.evaluateSolution(*irpStartSolution);

    this->numberOfOperations1 = irpStartSolution->getIrpSolution().getRepresentation().size();
    this->numberOfOperations2 = irpStartSolution->getIrpSolution().getRepresentation().size();
    cout<<"INITIAL BEGIN: "<<this->numberOfOperations1<<" "<<this->numberOfOperations2<<"\n";
        int a;cin>>a;
    return emili::Neighborhood::NeighborhoodIterator(this,startSolution);

}



emili::Solution* emili::irp::irpTwoExchangeNeighborhood::step(Solution* currentSolution){

    return this->computeStep(currentSolution);
}

emili::Solution* emili::irp::irpTwoExchangeNeighborhood::computeStep(Solution* currentSolution){

    InventoryRoutingSolution *neighboringSolution = dynamic_cast<InventoryRoutingSolution *> (currentSolution);
    neighboringSolution->getIrpSolution().fromSolutionToRepresentation(neighboringSolution->getIrpSolution());

    cout<<"NEIGH: "<<this->irp.evaluateSolution(*neighboringSolution);
    InventoryRoutingSolution *irs;

//    this->currentNeighboringSolution = dynamic_cast<Solution *> (neighboringSolution);

//    cout<<"INITIAL: "<<this->numberOfShifts1<<" "<<this->numberOfOperations1<<" "<<this->numberOfShifts2<<" "<<this->numberOfOperations2<<"\n";

    if(this->operation2 < this->numberOfOperations2 - 1){
        this->operation2++;
    }
    else if(this->operation1 < this->numberOfOperations1 - 1){
        this->operation1++;
        this->operation2 = operation1+1;
        cout<<this->operation1<<" "<<this->operation2<<"\n";
    }
    else
        return nullptr;

    int o1 = this->operation1;
    int o2 = this->operation2;
    unsigned int point1 = neighboringSolution->getIrpSolution().getRepresentation()[o1];
    unsigned int point2 = neighboringSolution->getIrpSolution().getRepresentation()[o2];
//    if(  not(  (o1==o2) or (point1==point2)    )  ){
        cout<<"\nEXCHANGE: "<<o1<<" "<<o2<<" "<<point1<<" "<<point2;

        this->irp.evaluateSolution(*neighboringSolution);
        neighboringSolution->getIrpSolution().fromSolutionToRepresentation(neighboringSolution->getIrpSolution());
        vector<unsigned int> representation = neighboringSolution->getIrpSolution().getRepresentation();
        representation[o1] = point2;
        representation[o2] = point1;
        irpSolution irps = this->irp.getIrpInstance().rebuildSolution(neighboringSolution->getIrpSolution(),representation, 0.0, 0.0);
        irps.fromSolutionToRepresentation(irps);
        irs = new InventoryRoutingSolution(irps);
        this->irp.evaluateSolution(*irs);


        this->numberOfOperations1 = irs->getIrpSolution().getRepresentation().size();
        this->numberOfOperations2 = irs->getIrpSolution().getRepresentation().size();


        cout<<"\nVALUE: "<<irs->getSolutionValue();
//
        if(irp.getIrpInstance().checkFeasibility(irs->getIrpSolution(), false))
            cout<<"\nNEIGH NOT FEASIBLE!\n";
        else if(irs->getSolutionValue() < neighboringSolution->getSolutionValue()){
            cout<<"\nNEIGH FEASIBLE!\n";
            this->numberFeasibleSolutions++;
            string filepath;
            filepath.append("./Neighborhood/");
            filepath.append(this->irp.getIrpInstance().getName());filepath.append("/");
            filepath.append(to_string(this->numberFeasibleSolutions));
            filepath.append("NeighSolution.xml");
            irs->getIrpSolution().saveSolution(filepath);

        }
//    }
    cout<<irs->getSolutionRepresentation();
//    int a; cin>>a;
    return dynamic_cast<Solution *> (irs);

}


void emili::irp::irpTwoExchangeNeighborhood::reset(){

    this->operation1 = POINT_INITIAL_VALUE;
    this->operation2 = POINT_INITIAL_VALUE;
    this->numberOfOperations1 = POINT_INITIAL_VALUE;
    this->numberOfOperations2 = POINT_INITIAL_VALUE;
//    this->numberFeasibleSolutions = 0;
}

emili::Solution* emili::irp::irpTwoExchangeNeighborhood::random(Solution* currentSolution){

        cout<<"\nTOWEXC RANDOM\n";int a;cin>>a;
    InventoryRoutingSolution *neighboringSolution = new InventoryRoutingSolution(*dynamic_cast<InventoryRoutingSolution *> (currentSolution));
    InventoryRoutingSolution *irs;
    irpSolution irps = neighboringSolution->getIrpSolution();
    irps.fromSolutionToRepresentation(irps);
//    neighboringSolution->getIrpSolution().fromSolutionToRepresentation(neighboringSolution->getIrpSolution());
    neighboringSolution = new InventoryRoutingSolution(irps);
    this->irp.evaluateSolution(*neighboringSolution);
    this->operation1 = POINT_INITIAL_VALUE;
    this->operation2 = POINT_INITIAL_VALUE;
    this->numberOfOperations1 = neighboringSolution->getIrpSolution().getRepresentation().size();
    this->numberOfOperations2 = neighboringSolution->getIrpSolution().getRepresentation().size();

    int o1 = generateRandomNumber() % this->numberOfOperations1;
    int o2 = generateRandomNumber() % this->numberOfOperations2;



    unsigned int point1 = neighboringSolution->getIrpSolution().getRepresentation()[o1]/*neighboringSolution->getSolutionRepresentation()[o1]*/;
    unsigned int point2 = neighboringSolution->getIrpSolution().getRepresentation()[o2]/*neighboringSolution->getSolutionRepresentation()[o2]*/;
    cout<<neighboringSolution->getSolutionRepresentation();
    cout<<"\n"<<o1<<" "<<o2<<"    "<<point1<<" "<<point2<<"\n";

    cout<<"BEFORE EXCHANGE VALUE: "<<neighboringSolution->getSolutionValue()<<"\n";
    for(int i=0; i<neighboringSolution->getSolutionRepresentation().size(); i++)
        cout<<neighboringSolution->getSolutionRepresentation()[i]<<" ";

//    if(  not(  (o1==o2) or (point1==point2)  )  ){

//        neighboringSolution->getIrpSolution().fromSolutionToRepresentation(neighboringSolution->getIrpSolution());
        vector<unsigned int> representation = neighboringSolution->getIrpSolution().getRepresentation();
        representation[o1] = point2;
        representation[o2] = point1;

        irps = this->irp.getIrpInstance().rebuildSolution(neighboringSolution->getIrpSolution(),representation, 0.0, 0.0);
        irs = new InventoryRoutingSolution(irps);
        irs->getIrpSolution().fromSolutionToRepresentation(irs->getIrpSolution());
        this->irp.evaluateSolution(*irs);
        cout<<" nEXC VALUE: "<<irs->getSolutionValue();

//    }

    return dynamic_cast<Solution *> (irs);
}

void emili::irp::irpTwoExchangeNeighborhood::reverseLastMove(Solution * step){

    cout<<"////////";
    cout<<"\nTWO EXC RATIO: "<<this->operation1<<" "<<this->operation2;

    /// //assegna un nuovo puntatore
//    step = this->currentNeighboringSolution;
    ///fa una copia
    *step = *currentNeighboringSolution;
    this->irp.evaluateSolution(*step);
    cout<<step->getSolutionRepresentation();
    cout<<"\nREVERSE VALUE: "<<step->getSolutionValue();
    cout<<"///////\n\n\n\n\n";
}

int emili::irp::irpTwoExchangeNeighborhood::size()
{
    return 1;
}


emili::Neighborhood::NeighborhoodIterator emili::irp::irpRefuelNeighborhood::begin(emili::Solution *startSolution)
{
    /////Salvo la soluzione iniziale per ripristinarla in reverse last move
    this->currentNeighboringSolution = startSolution->clone();

    InventoryRoutingSolution *irpStartSolution = dynamic_cast<InventoryRoutingSolution *> (startSolution);
    this->refuelRatio = REFUEL_INITIAL_VALUE;
    this->refuelStep = REFUEL_STEP;
    this->deliveredQuantityRatio = -REFUEL_STEP/*REFUEL_INITIAL_VALUE*/;
    this->deliveredQuantityStep = REFUEL_STEP;
//    this->numberFeasibleSolutions = 0;
//    this->bestValueFound = DBL_MAX;

    irpStartSolution->getIrpSolution().fromSolutionToRepresentation(irpStartSolution->getIrpSolution());
    irp.evaluateSolution(*irpStartSolution);

    cout<<"BEGIN INITIAL: "<<this->refuelRatio<<" "<<this->refuelStep<<" \n";
    startSolution = dynamic_cast<Solution *> (irpStartSolution);
    int a;cin>>a;
    return emili::Neighborhood::NeighborhoodIterator(this,startSolution);

}

emili::Solution* emili::irp::irpRefuelNeighborhood::step(Solution* currentSolution){

    return this->computeStep(currentSolution);
}

emili::Solution* emili::irp::irpRefuelNeighborhood::computeStep(Solution* currentSolution){

    InventoryRoutingSolution *neighboringSolution = dynamic_cast<InventoryRoutingSolution *> (currentSolution);

    cout<<"NEIGH: "<<this->irp.evaluateSolution(*neighboringSolution);

    if(this->deliveredQuantityRatio <= 1.0-this->deliveredQuantityStep and this->refuelRatio <= 1.0){
        this->deliveredQuantityRatio += this->deliveredQuantityStep;
    }
    else if(this->refuelRatio <= 1.0-this->refuelStep){
        this->deliveredQuantityRatio = REFUEL_INITIAL_VALUE;
        this->refuelRatio += this->refuelStep;
    }
    else
        return nullptr;

    cout<<"\nPERTURB RATIO: "<<this->deliveredQuantityRatio<<" "<<this->refuelRatio;

    //Creo la nuova soluzione
    vector<unsigned int> representation = neighboringSolution->getIrpSolution().getRepresentation();
    irpSolution irps = this->irp.getIrpInstance().rebuildSolution(neighboringSolution->getIrpSolution(),representation, this->refuelRatio, this->deliveredQuantityRatio);
    InventoryRoutingSolution *irs = new InventoryRoutingSolution(irps);
    irs->getIrpSolution().fromSolutionToRepresentation(irps);
    irp.evaluateSolution(*irs);

    if(irp.getIrpInstance().checkFeasibility(irs->getIrpSolution(), false))
        cout<<"\nNEIGH NOT FEASIBLE!\n";
    else if(irs->getSolutionValue() < this->bestValueFound){
        cout<<"\nNEIGH FEASIBLE!\n";
        this->numberFeasibleSolutions++;
        this->bestValueFound = irs->getSolutionValue();
        string filepath;
        filepath.append("./Neighborhood/");
        filepath.append(this->irp.getIrpInstance().getName());filepath.append("/");
        filepath.append(to_string(this->numberFeasibleSolutions));
        filepath.append("NeighSolution.xml");
        irs->getIrpSolution().saveSolution(filepath);
//        int a;cin>>a;
    }

   cout<<irs->getSolutionValue()<<"\n";
   cout<<irs->getSolutionRepresentation();
/*
   if(irs->getSolutionValue() < this->bestValueFound - EPSILON){
      this->bestValueFound = irs->getSolutionValue();
      cout<<"\nBEST!: "<<irs->getSolutionValue()<<"\n";
      int a; cin>>a;
   }
*/
//    int a;cin>>a;
   //Ritorno la nuova soluzione
    return dynamic_cast<Solution *> (irs);
}


void emili::irp::irpRefuelNeighborhood::reset(){

    cout<<"\nRESET\n";int a;cin>>a;
    this->refuelRatio = REFUEL_INITIAL_VALUE;
    this->refuelStep = REFUEL_STEP;
    this->deliveredQuantityRatio = REFUEL_INITIAL_VALUE;
    this->deliveredQuantityStep = REFUEL_STEP;
//    this->numberFeasibleSolutions = 0;
//    this->bestValueFound = DBL_MAX;
}

emili::Solution* emili::irp::irpRefuelNeighborhood::random(Solution* currentSolution){

    cout<<"\nREF RANDOM\n";

    InventoryRoutingSolution *neighboringSolution = new InventoryRoutingSolution(*dynamic_cast<InventoryRoutingSolution *> (currentSolution));
    InventoryRoutingSolution *irs;

    double randomRefuel = generateRealRandomNumber();
    double randomDeliveredQuantity = generateRealRandomNumber();

    neighboringSolution->getIrpSolution().fromSolutionToRepresentation(neighboringSolution->getIrpSolution());
    vector<unsigned int> representation = neighboringSolution->getIrpSolution().getRepresentation();
    irpSolution irps = this->irp.getIrpInstance().rebuildSolution(neighboringSolution->getIrpSolution(),representation, randomRefuel, randomDeliveredQuantity);
    irs = new InventoryRoutingSolution(irps);
    irs->getIrpSolution().fromSolutionToRepresentation(irs->getIrpSolution());
    this->irp.evaluateSolution(*irs);

    cout<<"EXC VALUE: "<<irs->getSolutionValue();
    int a;cin>>a;
    return dynamic_cast<Solution *> (irs);
}

void emili::irp::irpRefuelNeighborhood::reverseLastMove(Solution * step){

    cout<<"\n/////";

    //Ripristino la soluzione iniziale
    *step = *this->currentNeighboringSolution;


    cout<<step->getSolutionRepresentation();
    cout<<"\nREVERSE VALUE: "<<step->getSolutionValue();
    cout<<"\n/////";
}

int emili::irp::irpRefuelNeighborhood::size()
{
    return 1/(REFUEL_INITIAL_VALUE * REFUEL_INITIAL_VALUE);
}

emili::Solution* emili::irp::irpPerturbation::perturb(Solution* solution){

    InventoryRoutingSolution *perturbedSolution = new InventoryRoutingSolution(*dynamic_cast<InventoryRoutingSolution *> (solution));
    irpSolution irpPerturbedSolution = perturbedSolution->getIrpSolution();
}

// /home/antoniofisk/Desktop/Uni/MasterThesis/Project/build/Instance_V_1.1.xml IRP stin random
// /home/antoniofisk/Desktop/Uni/MasterThesis/Project/build/Instance_V_1.1.xml IRP best random locmin ref


// /home/antoniofisk/Desktop/Uni/MasterThesis/Project/build/Instance_V_1.1.xml IRP ils best random locmin ref true rndmv twoExc 1 metropolis 3.5 -it 1800
