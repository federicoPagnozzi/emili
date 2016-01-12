#include "irp.h"

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
    repr<<"\nSOLUTION REPRESENTATION: ";
    for(int s=0; s<this->irps.getShifts().size(); s++)
        for(int o=0; o<this->irps.getShifts()[s].getOperations().size(); o++)
            if(this->irps.getShifts()[s].getOperations()[o].getPoint()!=1)
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
   if(feas){
       double factor = feas/(this->irpInstance.getHorizon()*60 /** (this->irpInstance.getCustomers().size()-2)*/);
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

                if(tw == 0.0 and qw == 0.0 and t ==-1){
                    bestIrs = new InventoryRoutingSolution(irps);
                    instance.evaluateSolution(*bestIrs);
                    bestIrs->getIrpSolution().fromSolutionToRepresentation(bestIrs->getIrpSolution());
                }

////                COUT<<"PARAMETERS: "<<tw<<" "<<qw<<" "<<t<<"\n";
           /*     if(irp.getIrpInstance().checkFeasibility(irs->getIrpSolution(), false)){
                    COUT<<"\nORIGINAL NOT FEASIBLE!\n";
                    irs->getIrpSolution().saveSolution("OriginalSolution.xml");
                }
                else*/ if(irs->getSolutionValue() < bestValue){
                    COUT<<"\nORIGINAL FEASIBLE!\n";
                    bestIrs = new InventoryRoutingSolution(irps);
                    instance.evaluateSolution(*bestIrs);
                    bestIrs->getIrpSolution().fromSolutionToRepresentation(bestIrs->getIrpSolution());
                    bestValue = bestIrs->getSolutionValue();


                    COUT<<bestIrs->getSolutionRepresentation();
                    COUT<<"\nORIGINAL OBJ VALUE: "<<irs->getSolutionValue()<<"\n";
                    COUT<<"\nVALUES: "<<tw<<" "<<qw<<" "<<t<<"\n";
                    if(not (irp.getIrpInstance().checkFeasibility(irs->getIrpSolution(), false)))
                    {bf = true;break;
                        /*string filepath;
                        filepath.append("./Neighborhood/");
                        filepath.append(to_string(feasibleOriginalCounter));
                        filepath.append("OriginalSolution.xml");
                        irs->getIrpSolution().saveSolution(filepath);*/
                        feasibleOriginalCounter++;
                    
                    }
                }
    //            bf = true;
                if(bf)break;
            }
            if(bf)break;
        }
        if(bf)break;
    }

/*

    InventoryRoutingSolution *rirs;
    irpSolution rebuiltSolution = irp.getIrpInstance().rebuildSolution(bestIrs->getIrpSolution(),bestIrs->getIrpSolution().getRepresentation(), 0.0, 0.0, true);
    rirs = new InventoryRoutingSolution(rebuiltSolution);
    rirs->getIrpSolution().fromSolutionToRepresentation(rirs->getIrpSolution());
    instance.evaluateSolution(*rirs);
    COUT<<rirs->getSolutionRepresentation();
    if(irp.getIrpInstance().checkFeasibility(rirs->getIrpSolution(), false)){
        COUT<<"\nRECOSTRUCTION NOT FEASIBLE!\n";
    }
    else{
        COUT<<"\nRECOSTRUCTION FEASIBLE!\n";
//        rirs->getIrpSolution().saveSolution(*new string("RebuiltSolution.xml"));
        COUT<<"OBJ VALUE: "<<rirs->getSolutionValue()<<"\n\n";
    }
    COUT<<"OBJ VALUE: "<<rirs->getSolutionValue()<<"\n\n";
*/
    //return rirs;
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

emili::Neighborhood::NeighborhoodIterator emili::irp::irpTwoExchangeNeighborhood::begin(emili::Solution *startSolution)
{

    this->currentNeighboringSolution = startSolution->clone();
    InventoryRoutingSolution *irpStartSolution = dynamic_cast<InventoryRoutingSolution *> (startSolution);

    irpStartSolution->getIrpSolution().fromSolutionToRepresentation(irpStartSolution->getIrpSolution());
    unsigned int representationSize = irpStartSolution->getIrpSolution().getRepresentation().size();

    this->operation1 = 0;//(unsigned int) representationSize * (1.0/(double)this->pointInitialValue);
    this->operation2 = this->operation1;


    irp.evaluateSolution(*irpStartSolution);

    this->numberOfOperations1 = representationSize;
    this->numberOfOperations2 = representationSize;

    if(this->bestValueFound >= DBL_MAX - EPSILON){
        this->bestValueFound = irpStartSolution->getSolutionValue();
        this->numberFeasibleSolutions++;
  /*       string filepath;
         filepath.append("./Neighborhood/");
         filepath.append(this->irp.getIrpInstance().getName());filepath.append("/");
         filepath.append(to_string(this->numberFeasibleSolutions));
         filepath.append("NeighSolution.xml");
         irpStartSolution->getIrpSolution().saveSolution(filepath);*/
    }
//    COUT.clear();
    COUT<<"INITIAL BEGIN: "<<this->numberOfOperations1<<" "<<this->numberOfOperations2<<"\n";

//    std::COUT.setstate(std::ios_base::failbit);
//        int a;CIN>>a;
    return emili::Neighborhood::NeighborhoodIterator(this,startSolution);

}



emili::Solution* emili::irp::irpTwoExchangeNeighborhood::step(Solution* currentSolution){

    return this->computeStep(currentSolution);
}

emili::Solution* emili::irp::irpTwoExchangeNeighborhood::computeStep(Solution* currentSolution){

    InventoryRoutingSolution *neighboringSolution = dynamic_cast<InventoryRoutingSolution *> (currentSolution);

//   COUT<<"NEIGH: "<<this->irp.evaluateSolution(*neighboringSolution);

//    this->currentNeighboringSolution = dynamic_cast<Solution *> (neighboringSolution);


    if(this->operation2 < this->numberOfOperations2 - this->pointStep -1){
        this->operation2+=this->pointStep;
        this->operation1 = this->operation2+1;
    }
  /*  else if(this->operation1 < this->numberOfOperations1 - 1){
        this->operation1++;
        this->operation2 = operation1+1;
        COUT<<this->operation1<<" "<<this->operation2<<"\n";
    }*/
    else
        return nullptr;

    int o1 = this->operation1;
    int o2 = this->operation2;
    unsigned int point1 = neighboringSolution->getIrpSolution().getRepresentation()[o1];
    unsigned int point2 = neighboringSolution->getIrpSolution().getRepresentation()[o2];
       COUT<<"\nEXCHANGE: "<<o1<<" "<<o2<<" "<<point1<<" "<<point2<<"";

        this->irp.evaluateSolution(*neighboringSolution);
//        neighboringSolution->getIrpSolution().fromSolutionToRepresentation(neighboringSolution->getIrpSolution());
        vector<unsigned int> representation = neighboringSolution->getIrpSolution().getRepresentation();
        representation[o1] = point2;
        representation[o2] = point1;
        irpSolution irps = this->irp.getIrpInstance().rebuildSolution(neighboringSolution->getIrpSolution(),representation, 0.0, 0.0, false);
        irps.fromSolutionToRepresentation(irps);
        InventoryRoutingSolution irs(irps);
        this->irp.evaluateSolution(irs);


        this->numberOfOperations1 = irs.getIrpSolution().getRepresentation().size();
        this->numberOfOperations2 = irs.getIrpSolution().getRepresentation().size();

//        COUT<<"\nVALUE: "<<irs.getSolutionValue();
//
        if(irp.getIrpInstance().checkFeasibility(irs.getIrpSolution(), false))
            ;
        else if(irs.getSolutionValue() < this->bestValueFound){
//            COUT<<"\nNEIGH FEASIBLE!\n";
            this->numberFeasibleSolutions++;
            this->bestValueFound = irs.getSolutionValue();
         /*   string filepath;
            filepath.append("./Neighborhood/");
            filepath.append(this->irp.getIrpInstance().getName());filepath.append("/");
            filepath.append(to_string(this->numberFeasibleSolutions));
            filepath.append("NeighSolution.xml");
            irs.getIrpSolution().saveSolution(filepath);*/
            COUT<<"A BEST FOUND: "<<this->bestValueFound<<"\n";

        }
//    COUT<<irs.getSolutionRepresentation();
 //   int a; CIN>>a;

    /*
    * QUI COPIA lo stato interno di irs in currentSolution
    * neighboringSolution e currentSolution puntano allo stesso oggetto
    */
    *neighboringSolution = irs;
    return dynamic_cast<Solution *> (currentSolution);

}


void emili::irp::irpTwoExchangeNeighborhood::reset(){

    this->operation1 = this->pointInitialValue;
    this->operation2 = this->pointInitialValue;
    this->numberOfOperations1 = this->pointInitialValue;
    this->numberOfOperations2 = this->pointInitialValue;
//    this->numberFeasibleSolutions = 0;
}

emili::Solution* emili::irp::irpTwoExchangeNeighborhood::random(Solution* currentSolution){

        COUT<<"\nTOWEXC RANDOM\n";
        //int a;CIN>>a;
    InventoryRoutingSolution *neighboringSolution = new InventoryRoutingSolution(*dynamic_cast<InventoryRoutingSolution *> (currentSolution));
    InventoryRoutingSolution *irs;
    irpSolution irps = neighboringSolution->getIrpSolution();
    irps.fromSolutionToRepresentation(irps);
    neighboringSolution = new InventoryRoutingSolution(irps);
    this->irp.evaluateSolution(*neighboringSolution);
    this->operation1 = this->pointInitialValue;
    this->operation2 = this->pointInitialValue;
    this->numberOfOperations1 = neighboringSolution->getIrpSolution().getRepresentation().size();
    this->numberOfOperations2 = neighboringSolution->getIrpSolution().getRepresentation().size();

    int o1 = generateRandomNumber() % this->numberOfOperations1;
    int o2 = o1 + 1;/*generateRandomNumber() % this->numberOfOperations2;*/

    if(o2 >= this->numberOfOperations2){
        o2--;
        o1--;
    }


    unsigned int point1 = neighboringSolution->getIrpSolution().getRepresentation()[o1]/*neighboringSolution->getSolutionRepresentation()[o1]*/;
    unsigned int point2 = neighboringSolution->getIrpSolution().getRepresentation()[o2]/*neighboringSolution->getSolutionRepresentation()[o2]*/;
//    COUT<<neighboringSolution->getSolutionRepresentation();
    COUT<<"\n"<<o1<<" "<<o2<<"    "<<point1<<" "<<point2<<"\n";

    COUT<<"BEFORE EXCHANGE VALUE: "<<neighboringSolution->getSolutionValue()<<"\n";
//    for(int i=0; i<neighboringSolution->getSolutionRepresentation().size(); i++)
//        COUT<<neighboringSolution->getSolutionRepresentation()[i]<<" ";

//    if(  not(  (o1==o2) or (point1==point2)  )  ){

//        neighboringSolution->getIrpSolution().fromSolutionToRepresentation(neighboringSolution->getIrpSolution());
        vector<unsigned int> representation = neighboringSolution->getIrpSolution().getRepresentation();
        representation[o1] = point2;
        representation[o2] = point1;

        irps = this->irp.getIrpInstance().rebuildSolution(neighboringSolution->getIrpSolution(),representation, 0.0, 0.0, false);
        irs = new InventoryRoutingSolution(irps);
        irs->getIrpSolution().fromSolutionToRepresentation(irs->getIrpSolution());
        this->irp.evaluateSolution(*irs);
        COUT<<" nEXC VALUE: "<<irs->getSolutionValue();

//    }

    return dynamic_cast<Solution *> (irs);
}

void emili::irp::irpTwoExchangeNeighborhood::reverseLastMove(Solution * step){

//    COUT<<"////////";
//    COUT<<"\nTWO EXC RATIO: "<<this->operation1<<" "<<this->operation2;

    /// //assegna un nuovo puntatore
//    step = this->currentNeighboringSolution;
    ///fa una copia
    *step = *currentNeighboringSolution;
    this->irp.evaluateSolution(*step);
//    COUT<<step->getSolutionRepresentation();
//    COUT<<"\nREVERSE VALUE: "<<step->getSolutionValue();
//   COUT<<"///////\n\n\n\n\n";
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
    this->refuelRatio = this->refuelInitialValue;
    this->deliveredQuantityRatio = -this->refuelStep/*this->refuelInitialValue*/;
    this->deliveredQuantityStep = this->refuelStep;
//    this->numberFeasibleSolutions = 0;
//    this->bestValueFound = DBL_MAX;

    irpStartSolution->getIrpSolution().fromSolutionToRepresentation(irpStartSolution->getIrpSolution());
    irp.evaluateSolution(*irpStartSolution);

    if(this->bestValueFound >= DBL_MAX - EPSILON){
        this->bestValueFound = irpStartSolution->getSolutionValue();
        this->numberFeasibleSolutions++;
/*         string filepath;
         filepath.append("./Neighborhood/");
         filepath.append(this->irp.getIrpInstance().getName());filepath.append("/");
         filepath.append(to_string(this->numberFeasibleSolutions));
         filepath.append("NeighSolution.xml");
         irpStartSolution->getIrpSolution().saveSolution(filepath);*/
    }

//    COUT.clear();
    COUT<<"BEGIN INITIAL: "<<this->refuelRatio<<" "<<this->refuelStep<<" \n";
//    std::COUT.setstate(std::ios_base::failbit);
    startSolution = dynamic_cast<Solution *> (irpStartSolution);
//    int a;CIN>>a;
    return emili::Neighborhood::NeighborhoodIterator(this,startSolution);

}

emili::Solution* emili::irp::irpRefuelNeighborhood::step(Solution* currentSolution){

    return this->computeStep(currentSolution);
}

emili::Solution* emili::irp::irpRefuelNeighborhood::computeStep(Solution* currentSolution){

    InventoryRoutingSolution *neighboringSolution = dynamic_cast<InventoryRoutingSolution *> (currentSolution);

//    COUT<<"NEIGH: "<<this->irp.evaluateSolution(*neighboringSolution);

    if(this->deliveredQuantityRatio <= 1.0-this->deliveredQuantityStep and this->refuelRatio <= 1.0){
        this->deliveredQuantityRatio += this->deliveredQuantityStep;
    }
    else if(this->refuelRatio <= 1.0-this->refuelStep){
        this->deliveredQuantityRatio = this->refuelInitialValue;
        this->refuelRatio += this->refuelStep;
    }
    else
        return nullptr;

//    COUT<<"\nPERTURB RATIO: "<<this->deliveredQuantityRatio<<" "<<this->refuelRatio;

    //Creo la nuova soluzione
    vector<unsigned int> representation = neighboringSolution->getIrpSolution().getRepresentation();
    irpSolution irps = this->irp.getIrpInstance().rebuildSolution(neighboringSolution->getIrpSolution(),representation, this->refuelRatio, this->deliveredQuantityRatio, false);
    /*
		Se proprio è necessario creare una nuova soluzione ogni volta almeno creala sullo stack
		così viene cancellata quando questa funzione ritorna.
	*/
	InventoryRoutingSolution irs(irps);
    irs.getIrpSolution().fromSolutionToRepresentation(irps);
    irp.evaluateSolution(irs);

//   COUT<<irs.getSolutionValue()<<"\n";
//   COUT<<irs.getSolutionRepresentation();

    if(irp.getIrpInstance().checkFeasibility(irs.getIrpSolution(), false))
        ;
    else if(irs.getSolutionValue() < this->bestValueFound - EPSILON){
      this->bestValueFound = irs.getSolutionValue();
      this->numberFeasibleSolutions++;
    /*   string filepath;
       filepath.append("./Neighborhood/");
       filepath.append(this->irp.getIrpInstance().getName());filepath.append("/");
       filepath.append(to_string(this->numberFeasibleSolutions));
       filepath.append("NeighSolution.xml");
       irs.getIrpSolution().saveSolution(filepath);*/
//      COUT<<"\nBEST!: "<<irs.getSolutionValue()<<"\n";
//      int a; CIN>>a;
   }

//    int a;CIN>>a;
   //Ritorno la nuova soluzione
	/*
	* QUI COPIA lo stato interno di irs in currentSolution
	* neighboringSolution e currentSolution puntano allo stesso oggetto
	*/
	*neighboringSolution = irs;

    return currentSolution;
}


void emili::irp::irpRefuelNeighborhood::reset(){

//    COUT<<"\nRESET\n";
    //int a;CIN>>a;
    this->refuelRatio = this->refuelInitialValue;
    this->deliveredQuantityRatio = this->refuelInitialValue;
    this->deliveredQuantityStep = this->refuelStep;
//    this->numberFeasibleSolutions = 0;
//    this->bestValueFound = DBL_MAX;
}

emili::Solution* emili::irp::irpRefuelNeighborhood::random(Solution* currentSolution){

    COUT<<"\nREF RANDOM\n";

    InventoryRoutingSolution *neighboringSolution = new InventoryRoutingSolution(*dynamic_cast<InventoryRoutingSolution *> (currentSolution));
    InventoryRoutingSolution *irs;

    double randomRefuel = generateRealRandomNumber();
    double randomDeliveredQuantity = generateRealRandomNumber();

    neighboringSolution->getIrpSolution().fromSolutionToRepresentation(neighboringSolution->getIrpSolution());
    vector<unsigned int> representation = neighboringSolution->getIrpSolution().getRepresentation();
    irpSolution irps = this->irp.getIrpInstance().rebuildSolution(neighboringSolution->getIrpSolution(),representation, randomRefuel, randomDeliveredQuantity, false);
    irs = new InventoryRoutingSolution(irps);
    irs->getIrpSolution().fromSolutionToRepresentation(irs->getIrpSolution());
    this->irp.evaluateSolution(*irs);

    COUT<<"EXC VALUE: "<<irs->getSolutionValue();
//    int a;CIN>>a;
    return dynamic_cast<Solution *> (irs);
}

void emili::irp::irpRefuelNeighborhood::reverseLastMove(Solution * step){

//    COUT<<"\n/////";

    //Ripristino la soluzione iniziale
    *step = *this->currentNeighboringSolution;


//    COUT<<step->getSolutionRepresentation();
//    COUT<<"\nREVERSE VALUE: "<<step->getSolutionValue();


//    COUT<<"\n/////";
}

int emili::irp::irpRefuelNeighborhood::size()
{
    return 1/(this->refuelInitialValue * this->refuelInitialValue);
}

emili::Solution* emili::irp::irpPerturbation::perturb(Solution* solution){

    InventoryRoutingSolution *perturbedSolution = new InventoryRoutingSolution(*dynamic_cast<InventoryRoutingSolution *> (solution));
    irpSolution irpPerturbedSolution = perturbedSolution->getIrpSolution();
}

// /home/antoniofisk/Desktop/Uni/MasterThesis/Project/build/Instance_V_1.1.xml IRP stin random
// /home/antoniofisk/Desktop/Uni/MasterThesis/Project/build/Instance_V_1.1.xml IRP best random locmin ref


// /home/antoniofisk/Desktop/Uni/MasterThesis/Project/build/Instance_V_1.1.xml IRP ils best random locmin ref true rndmv twoExc 1 metropolis 3.5 -it 1800




//1) /home/antoniofisk/Desktop/Uni/MasterThesis/Project/build/Instance_V_1.3.xml IRP ils best random locmin ref true rndmv twoExc 1 metropolis 3.5 -it 1800
//2) /home/antoniofisk/Desktop/Uni/MasterThesis/Project/build/Instance_V_1.3.xml IRP ils best random locmin twoExc true rndmv ref 1 metropolis 3.5 -it 1800
