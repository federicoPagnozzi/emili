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

 /*   ofstream file;
    file.open ("./Ciao");
    file.precision(15);
    file<<"time : " << time_elapsed << std::endl;
    file<< "Objective function value: " << solution->getSolutionValue() << std::endl;
    file.close();*/
/*
    emili::irp::InventoryRoutingSolution* bestSolution = dynamic_cast<emili::irp::InventoryRoutingSolution*> (solution);
    string filepath;
    filepath.append("./BestSolution.xml");
    cout<<"\n"<<bestSolution->getIrpSolution().getShifts().size();
    bestSolution->getIrpSolution().saveSolution(filepath);
    */
}
