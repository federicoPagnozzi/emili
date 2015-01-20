#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <ctime>
#include <cstring>
#include "paramsparser.h"
#include "setup.h"


void g2c_info()
{
    std::cout << "usage in grammar2code mode : \n\tEMILI instance_file_path time random_seed" << std::endl;
    exit(0);
}


int main(int argc, char *argv[])
{


    /* initialize random seed: */
    srand ( time(0) );

    /* Create instance object */
    PfspInstance instance;
    clock_t time = clock();
#ifdef GRAMMAR2CODE
    instance.setSilence(true);
#else
    prs::emili();
#endif
    /* Read data from file */
    if (argc < 3 || !instance.readDataFromFile(argv[1]) )
    {        
#ifndef GRAMMAR2CODE
        prs::info();
#else

        g2c_info();
#endif
      return 1;
    }
    emili::pfsp::PermutationFlowShop problem(instance);    
    int pls = 0;
    emili::LocalSearch* ls;
#include "algorithm.h"
#ifndef GRAMMAR2CODE
    std::cout << "searching..." << std::endl;
    prs::ParamsParser ps(argv,argc,problem);
    ls = ps.parseParams();
    pls = ps.ils_time;
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
    std::vector < int > *sol = (std::vector < int > *)solution->getRawData();
    long int totalWeightedTardiness = instance.computeWT(*sol);
    double time_elapsed = (double)(clock()-time)/CLOCKS_PER_SEC;
    cout << "time : " << time_elapsed << std::endl;
    cout << "iteration counter " << emili::iteration_counter()<< std::endl;
    //cerr << time_elapsed << " ";
    cout << "Found solution: ";
    for (int i = 1; i <= instance.getNbJob(); ++i)
      cout << (*sol)[i] << " " ;

    cout << endl;
    cout << "Total weighted tardiness: " << totalWeightedTardiness << endl;
    cerr << totalWeightedTardiness << endl;
    exit(0);
}
