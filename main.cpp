#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <ctime>
#include <cstring>
#include "paramsparser.h"



int main(int argc, char *argv[])
{
    prs::emili();

    /* initialize random seed: */
    srand ( time(0) );

    /* Create instance object */
    PfspInstance instance;

    /* Read data from file */
    if (argc < 2 || !instance.readDataFromFile(argv[1]) )
    {
      prs::info();
      return 1;
    }
    emili::pfsp::PermutationFlowShop problem(instance);    
    std::cout << "searching..." << std::endl;
    int pls = 0;
    emili::LocalSearch* ls;
#include "algorithm.h"
#ifndef GRAMMAR2CODE
    prs::ParamsParser ps(argv,argc,problem);
    ls = ps.parseParams();
    pls = ps.ils_time;
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
    clock_t time = clock();
    double time_elapsed = (double)(clock()-time)/CLOCKS_PER_SEC;
    cout << "time : " << time_elapsed << std::endl;
    std::vector < int > *sol = (std::vector < int > *)solution->getRawData();
    cout << "iteration counter " << emili::iteration_counter()<< std::endl;
    cerr << time_elapsed << " ";
    cout << "Found solution: ";
    for (int i = 1; i <= instance.getNbJob(); ++i)
      cout << (*sol)[i] << " " ;
    cout << endl;
    long int totalWeightedTardiness = instance.computeWT(*sol);
    cout << "Total weighted tardiness: " << totalWeightedTardiness << endl;
    cerr << totalWeightedTardiness << endl;
    exit(0);
}


