#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <ctime>
#include <cstring>
#include <algorithm>
#include "generalParser.h"
#include "setup.h"

#include "irp.h"

#include "irpparser.h"

void g2c_info()
{
    std::cout << "usage in grammar2code mode : \n\tEMILI instance_file_path time random_seed" << std::endl;
    exit(0);
}


int main(int argc, char *argv[])
{

    clock_t time = clock();

    float maxTime;
    char *instanceName;
    char *solutionName;
    int teamID = 42;
    int randomSeed = 0;

    bool isDeterministic = false;
    if(argc != 2 and argc != 8 and argc != 10){
        cout<<"\nBad arguments\n";
        exit(0);
    }
    else if(argc==2 and strcmp(argv[1],"-name") == 0){
        cout<<"\nTEAM ID: "<<teamID<<"\n";
        exit(0);
    }
    else{
        for (int i = 1; i < argc; i++) {
            if (i + 1 != argc){
                if (strcmp(argv[i],"-t") == 0) {
                    maxTime = atof(argv[i + 1]);
                    cout<<"\nTIME: "<<maxTime;
                } else if (strcmp(argv[i],"-p") == 0) {
                    instanceName = argv[i + 1];
                    cout<<"\nINSTANCE: "<<instanceName;
                } else if (strcmp(argv[i],"-o") == 0) {
                    solutionName = argv[i + 1];
                    cout<<"\nSOLUTION: "<<solutionName;
                } else if (strcmp(argv[i],"-name") == 0) {
                    cout<<"\nTEAM ID: "<<teamID;
                } else if (strcmp(argv[i],"-s") == 0) {
                    randomSeed = atoi(argv[i + 1]);
                    cout<<"\nSEED: "<<randomSeed;
                    isDeterministic = true;
                }
            }
        }

        cout<<"\n";
    }


    if(isDeterministic)
     emili::initializeRandom(randomSeed);
/*
    emili::irp::InventoryRoutingProblem *instance = new emili::irp::InventoryRoutingProblem(instanceName);
       emili::InitialSolution* in = new emili::irp::GreedyInitialSolution(*instance);
       emili::Termination* te= new emili::LocalMinimaTermination();
       double riv = 0.0;
       double rs = 0.1;
       emili::Neighborhood* ne = new emili::irp::irpRefuelNeighborhood(*instance, riv, rs);
       emili::LocalSearch* ils =  new emili::FirstImprovementSearch(*in,*te,*ne);


       emili::Termination* pft = new emili::TimedTermination(maxTime);

       unsigned int  piv = 5;
       unsigned int ps2 = 3;
       emili::Neighborhood* n = new emili::irp::irpTwoExchangeNeighborhood(*instance, piv, ps2);
       unsigned int num = 3;
       emili::Perturbation* prsp = new emili::RandomMovePertubation(*n,num);

       double start = 1.5613;
       double end = 0.4880;
       double ratio = 0.0488;
       emili::Acceptance* tac = new  emili::Metropolis(start,end,ratio);
       emili::LocalSearch*ls = new emili::IteratedLocalSearch(*ils,*pft,*prsp,*tac);
   //    ls->setSearchTime(getTime(tm,ls->getInitialSolution().getProblem().problemSize()));
       emili::Solution* solution;
       solution = ls->search();
       solution = ls->getBestSoFar();
*/
/*
    emili::irp::InventoryRoutingProblem *instance = new emili::irp::InventoryRoutingProblem(instanceName);
    emili::InitialSolution* in = new emili::irp::GreedyInitialSolution(*instance);
    emili::Termination* te= new emili::LocalMinimaTermination();
    unsigned int  piv = 2;
    unsigned int ps2 = 1;
    emili::Neighborhood* ne = new emili::irp::irpTwoExchangeNeighborhood(*instance, piv, ps2);
    emili::LocalSearch* ils =  new emili::FirstImprovementSearch(*in,*te,*ne);


    emili::Termination* pft = new emili::TimedTermination(maxTime);
    piv = 0;
    ps2 = 3;
    emili::Neighborhood* n = new emili::irp::irpTwoExchangeNeighborhood(*instance, piv, ps2);
    unsigned int num = 3;
    emili::Perturbation* prsp = new emili::RandomMovePertubation(*n,num);
    emili::Acceptance* tac = new emili::ImproveAccept();
    emili::LocalSearch*ls = new emili::IteratedLocalSearch(*ils,*pft,*prsp,*tac);
//    ls->setSearchTime(getTime(tm,ls->getInitialSolution().getProblem().problemSize()));
    emili::Solution* solution;
    solution = ls->search();
    solution = ls->getBestSoFar();
*/

    emili::irp::InventoryRoutingProblem *instance = new emili::irp::InventoryRoutingProblem(instanceName);
    emili::InitialSolution* in = new emili::irp::GreedyInitialSolution(*instance);
    emili::Termination* te= new emili::LocalMinimaTermination();
    unsigned int  piv = -1;
    unsigned int ps2 = 1;
    emili::Neighborhood* ne = new emili::irp::irpTwoExchangeNeighborhood(*instance, piv, ps2);
    emili::LocalSearch* ils =  new emili::FirstImprovementSearch(*in,*te,*ne);


    emili::Termination* pft = new emili::TimedTermination(maxTime);
    double rs = 0.1;
    double riv = -rs;
    emili::Neighborhood* n = new emili::irp::irpRefuelNeighborhood(*instance, riv, rs);
    unsigned int num = 1;
    emili::Perturbation* prsp = new emili::RandomMovePertubation(*n,num);
    emili::Acceptance* tac = new  emili::MetropolisAcceptance(3.5);
    emili::LocalSearch*ls = new emili::IteratedLocalSearch(*ils,*pft,*prsp,*tac);
//    ls->setSearchTime(getTime(tm,ls->getInitialSolution().getProblem().problemSize()));
    emili::Solution* solution;
    solution = ls->search();
    solution = ls->getBestSoFar();

    double time_elapsed = (double)(clock()-time)/CLOCKS_PER_SEC;
    std::cout << "time : " << time_elapsed << std::endl;
    std::cout << "iteration counter : " << emili::iteration_counter()<< std::endl;
    std::cerr << solution->getSolutionValue() << std::endl;
    //cerr << time_elapsed << " ";
    std::cout << "Objective function value: " << solution->getSolutionValue() << std::endl;
    std::cout << "Found solution: ";
    std::cout << solution->getSolutionRepresentation() << std::endl;
//    std::cout << std::endl;

    emili::irp::InventoryRoutingSolution* bestSolution = dynamic_cast<emili::irp::InventoryRoutingSolution*> (solution);


    string filepath;
 //   filepath.append("./Neighborhood/");
    filepath.append("./");
//    filepath.append(this->irp.getIrpInstance().getName());filepath.append("/");
//    filepath.append(to_string(this->numberFeasibleSolutions));
//    filepath.append("./");
    filepath.append(solutionName);filepath.append(".xml");
    bestSolution->getIrpSolution().saveSolution(filepath);


}

int main2(int argc, char *argv[])
{
//std::cout.setstate(std::ios_base::failbit);

prs::emili_header();
    /* initialize random seed: */
    srand ( time(0) );

    /* Create instance object */
    //testTaillardAccel();
    clock_t time = clock();
 //instance.setSilence(true);
#ifdef GRAMMAR2CODE


#else
    //prs::emili_header();
#endif
    /* Read data from file */
    if (argc < 3 )
    {        
#ifndef GRAMMAR2CODE
        prs::info();
#else

        g2c_info();
#endif
      return 1;
    }
   // testNewEvaluationFunction(instance);
    //emili::pfsp::NWPFSP_MS problem(instance);
    //testHeuritstic(problem);
    int pls = 0;
    emili::LocalSearch* ls;
#include "algorithm.h"
#ifndef GRAMMAR2CODE
    std::cout << "searching..." << std::endl;

    prs::GeneralParser ps(argv,argc);

    /*

       HERE you must register your algobuilder to the general parser!!!

     */

    //prs::ParamsParser p;
    //ps.registerBuilder(&p);
    prs::irp::IrpParser p;
    ps.registerBuilder(&p);

    ls = ps.parseParams();
   // testHeuritstic(ps.getInstance());
    if(ls==nullptr)
    {
        return -1;
    }
    pls = ls->getSearchTime();//ps.ils_time;
#else
    pls = atoi(argv[2]);
    int seed = atoi(argv[3]);
    emili::initializeRandom(seed);
    time = clock();
#endif
    emili::Solution* solution;
    if(pls>0)
    {
       solution = ls->timedSearch(pls);
    }
    else
    {
        solution = ls->search();
    }

#ifndef GRAMMAR2CODE

#else
    long int totalWeightedTardiness = problem.computeObjectiveFunction(sol);
    int njobs = problem.getNjobs();
#endif


    solution = ls->getBestSoFar();    
    double time_elapsed = (double)(clock()-time)/CLOCKS_PER_SEC;
    std::cout << "time : " << time_elapsed << std::endl;
    std::cout << "iteration counter : " << emili::iteration_counter()<< std::endl;
    std::cerr << solution->getSolutionValue() << std::endl;
    //cerr << time_elapsed << " ";    
    std::cout << "Objective function value: " << solution->getSolutionValue() << std::endl;
    std::cout << "Found solution: ";
    std::cout << solution->getSolutionRepresentation() << std::endl;
    std::cout << std::endl;

    ofstream file;
    file.open ("./Ciao",fstream::app);
    file.precision(15);
    file<<"time : " << time_elapsed << std::endl;
    file<< "Objective function value: " << solution->getSolutionValue() << std::endl;
    file.close();


    string filepath;
    emili::irp::InventoryRoutingSolution* bestSolution = dynamic_cast<emili::irp::InventoryRoutingSolution*> (solution);
    filepath.append("./BestSolution.xml");
    cout<<"\n"<<bestSolution->getIrpSolution().getShifts().size();
    bestSolution->getIrpSolution().saveSolution(filepath);

}
