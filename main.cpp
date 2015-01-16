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

    emili::pfsp::PfspRandomInitialSolution r(problem);

    emili::Solution* j = r.generateSolution();

    std::vector< int >* test = (std::vector< int >*)j->getRawData();

    std::vector<std::vector <int> > pt(instance.getNbMac()+1,std::vector<int>(problem.getNjobs()+1,0));
    long int cazz;
    clock_t asd = clock();
    for(int i=0;i<1000;i++){
    cazz = problem.computeWT(*test);
    }
    double time_e = (double)(clock()-asd)/CLOCKS_PER_SEC;
    std::cout << "normale -> " << cazz  << " tempo-> " << time_e << std::endl;
    asd = clock();
    for(int i=0;i<1000;i++){
    cazz = instance.computeWT(*test,pt,1,101);
    }
    time_e = (double)(clock()-asd)/CLOCKS_PER_SEC;
    std::cout << "nuovo -> " << cazz  << " tempo-> " << time_e << std::endl;



    std::swap((*test)[70],(*test)[79]);
    asd = clock();
    for(int i=0;i<1000;i++){

    cazz = problem.computeWT(*test);
    }
    time_e = (double)(clock()-asd)/CLOCKS_PER_SEC;
    std::cout << "normale -> " << cazz  << " tempo-> " << time_e << std::endl;

    asd = clock();
    for(int i=0;i<1000;i++){

    cazz = instance.computeWT(*test,pt,70,79);
    }
    time_e = (double)(clock()-asd)/CLOCKS_PER_SEC;
    std::cout << "nuovo -> " << cazz  << " tempo-> " << time_e << std::endl;
    exit(0);

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
#ifndef GRAMMAR2CODE    
    double time_elapsed = (double)(clock()-time)/CLOCKS_PER_SEC;
    cout << "time : " << time_elapsed << std::endl;
    cout << "iteration counter " << emili::iteration_counter()<< std::endl;
    cerr << time_elapsed << " ";
    cout << "Found solution: ";
    for (int i = 1; i <= instance.getNbJob(); ++i)
      cout << (*sol)[i] << " " ;

    cout << endl;
    cout << "Total weighted tardiness: " << totalWeightedTardiness << endl;
    cerr << totalWeightedTardiness << endl;
#else
    std::cout << totalWeightedTardiness << std::endl;
#endif
    exit(0);
}
