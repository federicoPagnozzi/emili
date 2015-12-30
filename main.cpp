#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <ctime>
#include <cstring>
#include <algorithm>
#include "generalParser.h"
#include "setup.h"

#include "irpparser.h"

void g2c_info()
{
    std::cout << "usage in grammar2code mode : \n\tEMILI instance_file_path time random_seed" << std::endl;
    exit(0);
}

void test(){

    emili::irp::InventoryRoutingProblem *irp = new emili::irp::InventoryRoutingProblem("Instance_V_1.1.xml");
    irpSolution s = irp->getIrpInstance().backTrackingRandomSolution(1.0, 1.0, 1);
    emili::irp::InventoryRoutingSolution *irs = new emili::irp::InventoryRoutingSolution(s);
    if(irp->getIrpInstance().checkFeasibility(irs->getIrpSolution(), false))
        cout<<"\nNOT FEASIBLE! \n";
    else
        cout<<"\nFEASIBLE! \n";
    cout<<"CIAOOO";
//    irs->getIrpSolution().saveSolution("Solution_V_1.1.xml");
//    emili::irp::irpPerturbation perturbation = emili::irp::irpPerturbation();
//    perturbation.perturb(irs);

}

int main(int argc, char *argv[])
{

 //   test();

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
}
