#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <ctime>
#include <cstring>
#include <algorithm>
#include "paramsparser.h"
#include "setup.h"


void g2c_info()
{
    std::cout << "usage in grammar2code mode : \n\tEMILI instance_file_path time random_seed" << std::endl;
    exit(0);
}


void testNewEvaluationFunction(PfspInstance& instance)
{

    std::vector<int > prevJob(instance.getNbMac(),0);
    std::vector<int > previousMachineEndTime(instance.getNbJob()+1,0);
    emili::pfsp::PFSP_WT pro(instance);
    emili::pfsp::PfspNEHwslackInitialSolution p(pro);
    emili::pfsp::PermutationFlowShopSolution* s = (emili::pfsp::PermutationFlowShopSolution*) p.generateSolution();
    std::vector<int > sol(s->getJobSchedule());
    std::vector<int > test(sol);
    std::vector<int > test1(sol);

    clock_t start = clock();
    long k = 0;

    for (int var = 1; var < instance.getNbJob(); ++var) {
       std::swap(test[var],test[var+1]);
        k = instance.computeWT(test,prevJob,var,previousMachineEndTime);

    }
    std::cout << " New time -> " << (clock()-start)/(float)CLOCKS_PER_SEC << std::endl;
    std::cout << " Value -> " << k << std::endl;

    start = clock();
    k = 0;



    for (int var = 1; var < instance.getNbJob(); ++var) {
        std::swap(test1[var],test1[var+1]);
        k = instance.computeWT(test1);

    }
    std::cout << " Normal time -> " << (clock()-start)/(float)CLOCKS_PER_SEC << std::endl;
    std::cout << " Value -> " << k << std::endl;
    exit(0);
}

void testHeuritstics2(emili::pfsp::PermutationFlowShop& problem){
    emili::pfsp::PfspRandomInitialSolution rnd(problem);
    emili::pfsp::LITSolution tests(problem);
    emili::pfsp::PfspNEHwslackInitialSolution nwslack(problem);
    emili::pfsp::PfspSlackInitialSolution slack(problem);
    emili::pfsp::RZSolution rz(problem);
    emili::pfsp::NeRZ2Solution nrz2(problem);
    emili::pfsp::NeRZSolution nrz(problem);
    emili::pfsp::LRSolution lr(problem,10);
    emili::pfsp::NLRSolution nlr(problem,10);

    clock_t start = clock();
    emili::pfsp::PermutationFlowShopSolution* sorl = (emili::pfsp::PermutationFlowShopSolution*)tests.generateSolution();
    std::cerr <<(clock()-start)/(float)CLOCKS_PER_SEC;
    std::cout << sorl->getSolutionValue();

    start = clock();
    emili::pfsp::PermutationFlowShopSolution* rns = (emili::pfsp::PermutationFlowShopSolution*)rnd.generateSolution();
    std::cerr << "\t" <<(clock()-start)/(float)CLOCKS_PER_SEC;
    std::cout << "\t " << rns->getSolutionValue();

    start = clock();
    emili::pfsp::PermutationFlowShopSolution* sss =(emili::pfsp::PermutationFlowShopSolution*) slack.generateSolution();
    std::cerr << "\t" <<(clock()-start)/(float)CLOCKS_PER_SEC;
    std::cout << "\t" << sss->getSolutionValue();

    start = clock();
    emili::pfsp::PermutationFlowShopSolution* nws = (emili::pfsp::PermutationFlowShopSolution*)nwslack.generateSolution();
    std::cerr << "\t" << (clock()-start)/(float)CLOCKS_PER_SEC;
    std::cout << "\t" << nws->getSolutionValue();

    start = clock();
    emili::pfsp::PermutationFlowShopSolution* rzs = (emili::pfsp::PermutationFlowShopSolution*)rz.generateSolution();
    std::cerr << "\t" <<(clock()-start)/(float)CLOCKS_PER_SEC;
    std::cout << "\t" << rzs->getSolutionValue();

    start = clock();
    emili::pfsp::PermutationFlowShopSolution* nrzs = (emili::pfsp::PermutationFlowShopSolution*)nrz.generateSolution();
    std::cerr << "\t" << (clock()-start)/(float)CLOCKS_PER_SEC;
    std::cout << "\t" << nrzs->getSolutionValue();

    start = clock();
    emili::pfsp::PermutationFlowShopSolution* nrz2s = (emili::pfsp::PermutationFlowShopSolution*)nrz2.generateSolution();
    std::cerr << "\t" << (clock()-start)/(float)CLOCKS_PER_SEC;
    std::cout << "\t" << nrz2s->getSolutionValue();

    start = clock();
    emili::pfsp::PermutationFlowShopSolution* lrs = (emili::pfsp::PermutationFlowShopSolution*)lr.generateSolution();
    std::cerr << "\t" << (clock()-start)/(float)CLOCKS_PER_SEC;
    std::cout << "\t" << lrs->getSolutionValue();

    start = clock();
    emili::pfsp::PermutationFlowShopSolution* nlrs = (emili::pfsp::PermutationFlowShopSolution*)nlr.generateSolution();
    std::cerr << "\t" << (clock()-start)/(float)CLOCKS_PER_SEC << std::endl;
    std::cout << "\t" << nlrs->getSolutionValue() << std::endl;
    exit(0);
}

void testHeuritstic(emili::pfsp::PermutationFlowShop& problem){
    emili::pfsp::PfspRandomInitialSolution rnd(problem);
    emili::pfsp::LITSolution tests(problem);
    emili::pfsp::PfspNEHwslackInitialSolution nwslack(problem);
    emili::pfsp::PfspSlackInitialSolution slack(problem);
    emili::pfsp::RZSolution rz(problem);
    emili::pfsp::NeRZ2Solution nrz2(problem);
    emili::pfsp::NeRZSolution nrz(problem);
    emili::pfsp::LRSolution lr(problem,10);
    emili::pfsp::NLRSolution nlr(problem,10);
    emili::pfsp::MNEH mnr(problem);

    clock_t start = clock();
    emili::pfsp::PermutationFlowShopSolution* sorl = (emili::pfsp::PermutationFlowShopSolution*)tests.generateSolution();
    std::cout << " LIT time -> " << (clock()-start)/(float)CLOCKS_PER_SEC << std::endl;
    std::vector< int > sol =  sorl->getJobSchedule();
    cout << "Found solution: ";
    for (int i = 1; i <= problem.getNjobs(); ++i)
      cout << sol[i] << " " ;
    cout << endl;
    std::cout << "LIT -> " << sorl->getSolutionValue() << std::endl;

    start = clock();
    emili::pfsp::PermutationFlowShopSolution* rns = (emili::pfsp::PermutationFlowShopSolution*)rnd.generateSolution();
    std::cout << " RANDOM time -> " << (clock()-start)/(float)CLOCKS_PER_SEC << std::endl;
    sol =  rns->getJobSchedule();
    cout << "Found solution: ";
    for (int i = 1; i <= problem.getNjobs(); ++i)
      cout << sol[i] << " " ;
    cout << endl;
    std::cout << "rns -> " << rns->getSolutionValue() << std::endl;

    start = clock();
    emili::pfsp::PermutationFlowShopSolution* sss = (emili::pfsp::PermutationFlowShopSolution*)slack.generateSolution();
    std::cout << " SLACK time -> " << (clock()-start)/(float)CLOCKS_PER_SEC << std::endl;
    sol =  sss->getJobSchedule();
    cout << "Found solution: ";
    for (int i = 1; i <= problem.getNjobs(); ++i)
      cout << sol[i] << " " ;
    cout << endl;
    std::cout << "sss -> " << sss->getSolutionValue() << std::endl;

    start = clock();
    emili::pfsp::PermutationFlowShopSolution* nws = (emili::pfsp::PermutationFlowShopSolution*)nwslack.generateSolution();
    std::cout << " NWSLACK time -> " << (clock()-start)/(float)CLOCKS_PER_SEC << std::endl;
    sol =  nws->getJobSchedule();
    cout << "Found solution: ";
    for (int i = 1; i <= problem.getNjobs(); ++i)
      cout << sol[i] << " " ;
    cout << endl;
    std::cout << "nws -> " << nws->getSolutionValue() << std::endl;

    start = clock();
    emili::pfsp::PermutationFlowShopSolution* rzs = (emili::pfsp::PermutationFlowShopSolution*)rz.generateSolution();
    std::cout << " RZ time -> " << (clock()-start)/(float)CLOCKS_PER_SEC << std::endl;
    sol =  rzs->getJobSchedule();
    cout << "Found solution: ";
    for (int i = 1; i <= problem.getNjobs(); ++i)
      cout << sol[i] << " " ;
    cout << endl;
    std::cout << "rz -> " << rzs->getSolutionValue() << std::endl;

    start = clock();
    emili::pfsp::PermutationFlowShopSolution* nrzs = (emili::pfsp::PermutationFlowShopSolution*)nrz.generateSolution();
    std::cout << " NRZ time -> " << (clock()-start)/(float)CLOCKS_PER_SEC << std::endl;
    sol =  nrzs->getJobSchedule();
    cout << "Found solution: ";
    for (int i = 1; i <= problem.getNjobs(); ++i)
      cout << sol[i] << " " ;
    cout << endl;
    std::cout << "nrz -> " << nrzs->getSolutionValue() << std::endl;

    start = clock();
    emili::pfsp::PermutationFlowShopSolution* nrz2s =(emili::pfsp::PermutationFlowShopSolution*) nrz2.generateSolution();
    std::cout << " NRZ2 time -> " << (clock()-start)/(float)CLOCKS_PER_SEC << std::endl;
    sol =  nrz2s->getJobSchedule();
    cout << "Found solution: ";
    for (int i = 1; i <= problem.getNjobs(); ++i)
      cout << sol[i] << " " ;
    cout << endl;
    std::cout << "nrz2 -> " << nrz2s->getSolutionValue() << std::endl;

    start = clock();
    emili::pfsp::PermutationFlowShopSolution* lrs = (emili::pfsp::PermutationFlowShopSolution*)lr.generateSolution();
    std::cout << " LR time -> " << (clock()-start)/(float)CLOCKS_PER_SEC << std::endl;
    sol =  lrs->getJobSchedule();
    cout << "Found solution: ";
    for (int i = 1; i <= problem.getNjobs(); ++i)
      cout << sol[i] << " " ;
    cout << endl;
    std::cout << "LR -> " << lrs->getSolutionValue() << std::endl;

    start = clock();
    emili::pfsp::PermutationFlowShopSolution* nlrs = (emili::pfsp::PermutationFlowShopSolution*)nlr.generateSolution();
    std::cout << " NLR time -> " << (clock()-start)/(float)CLOCKS_PER_SEC << std::endl;
    sol =  nlrs->getJobSchedule();
    cout << "Found solution: ";
    for (int i = 1; i <= problem.getNjobs(); ++i)
      cout << sol[i] << " " ;
    cout << endl;
    std::cout << "NLR -> " << nlrs->getSolutionValue() << std::endl;


    start = clock();
    emili::pfsp::PermutationFlowShopSolution* mnehs = (emili::pfsp::PermutationFlowShopSolution*)mnr.generateSolution();
    std::cout << " NLR time -> " << (clock()-start)/(float)CLOCKS_PER_SEC << std::endl;
    sol =  mnehs->getJobSchedule();
    cout << "Found solution: ";
    for (int i = 1; i <= problem.getNjobs(); ++i)
      cout << sol[i] << " " ;
    cout << endl;
    std::cout << "MNEH -> " << mnehs->getSolutionValue() << std::endl;
    exit(0);
}


    //char* file = "/Users/federicopagnozzi/sviluppo/QTCREAT/EMILI-build/instances/testInstance";
    //Users/federicopagnozzi/Desktop/phd/PFSWTinstances/Taillard_DueDates/DD_Ta086.txt
    //char* file = "/Users/federicopagnozzi/sviluppo/QTCREAT/EMILI-build/instances/testruiz.txt";




int main(int argc, char *argv[])
{

    /* initialize random seed: */
    srand ( time(0) );

    /* Create instance object */
    //testTaillardAccel();
    clock_t time = clock();
 //instance.setSilence(true);
#ifdef GRAMMAR2CODE


#else
    prs::emili_header();
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
    prs::ParamsParser ps(argv,argc);
    ls = ps.parseParams();
   // testHeuritstic(ps.getInstance());
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
    std::vector < int >& sol = ((emili::pfsp::PermutationFlowShopSolution*)solution)->getJobSchedule();
#ifndef GRAMMAR2CODE
    emili::pfsp::PermutationFlowShop& prob = ps.getInstance(); 
    long int totalWeightedTardiness = prob.computeObjectiveFunction(sol);
    int njobs = prob.getNjobs();
#else
    long int totalWeightedTardiness = problem.computeObjectiveFunction(sol);
    int njobs = problem.getNjobs();
#endif    
    double time_elapsed = (double)(clock()-time)/CLOCKS_PER_SEC;
    cout << "time : " << time_elapsed << std::endl;
    cout << "iteration counter " << emili::iteration_counter()<< std::endl;
    //cerr << time_elapsed << " ";
    cout << "Found solution: ";
    for (int i = 1; i <= njobs; ++i)
      cout << sol[i] << " " ;

    cout << endl;
    cout << "Objective function value: " << totalWeightedTardiness << endl;
    cerr << totalWeightedTardiness << endl;
    exit(0);
}
