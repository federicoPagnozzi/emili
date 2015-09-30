#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <ctime>
#include <cstring>
#include <algorithm>
#include "generalParser.h"
// #include "paramsparser.h"
#include "setup.h"

#include "SA/sa_pfsp_parser.h"
#include "SA/sa_qap_parser.h"
#include "QAP/qapinitialsolution.h"
#include "QAP/qapneighborhood.h"
#include "QAP/qap.h"


void g2c_info()
{
    std::cout << "usage in grammar2code mode : \n\tEMILI instance_file_path time random_seed" << std::endl;
    exit(0);
}
int main(int argc, char *argv[])
{

    /** /emili::initializeRandom(atoi(argv[2]));

    string file(argv[1]);
    QAPInstance* inst = new QAPInstance(file);

    std::cout << inst->toString() << std::endl;

    qap::QAP* prob = new qap::QAP(*inst);

    QAPRandomInitialSolution* initsol = new QAPRandomInitialSolution(*prob);
    emili::Solution* sol = initsol->generateSolution();

    std::cout << sol->getSolutionRepresentation() << std::endl;
    std::cout << sol->getSolutionValue() << std::endl;


    QAPFirst2optNeighborhood* neigh2 = new QAPFirst2optNeighborhood(*prob);

    emili::Solution* solf = neigh2->step(sol);

    std::cout << solf->getSolutionRepresentation() << std::endl;
    std::cout << solf->getSolutionValue() << std::endl;

    QAPBest2optNeighborhood* neigh = new QAPBest2optNeighborhood(*prob);

    emili::Solution* sol2 = neigh->step(sol);

    std::cout << sol2->getSolutionRepresentation() << std::endl;
    std::cout << sol2->getSolutionValue() << std::endl;


    QAPExchangeNeighborhood* neighe = new QAPExchangeNeighborhood(*prob);

    emili::Solution* incumbent = neighe->random(sol2);

    std::cout << incumbent->getSolutionRepresentation() << std::endl;
    std::cout << incumbent->getSolutionValue() << std::endl;

    return 0;

    emili::Solution* ithSolution = nullptr;
    emili::Solution* bestOfTheIteration = sol2;
    for(emili::Neighborhood::NeighborhoodIterator iter = neighe->begin(sol2);iter!=neighe->end();++iter)
    {
        ithSolution = *iter;

        std::cout << ithSolution->getSolutionRepresentation() << std::endl;
        std::cout << ithSolution->getSolutionValue() << std::endl;

        /*if(bestOfTheIteration->operator >( *ithSolution)){
            if(bestOfTheIteration!=sol2)
            delete bestOfTheIteration;

            bestOfTheIteration = ithSolution;

        }
        else
        {
            delete ithSolution;
        }* /

    }


    return 0;/ **/

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
    //SAPFSPParser p;
    SAQAPParser p;
    prs::GeneralParser ps(argv,argc);
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
    double time_elapsed = (double)(clock()-time)/CLOCKS_PER_SEC;
    cout << "time : " << time_elapsed << std::endl;
    cout << "iteration counter : " << emili::iteration_counter()<< std::endl;
    //cerr << time_elapsed << " ";    
    cout << "Objective function value: " << solution->getSolutionValue() << endl;
    cout << "Found solution: ";
    cout << solution->getSolutionRepresentation() << std::endl;
    cout << endl;

    cerr << solution->getSolutionValue() << endl;
}
