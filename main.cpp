//
//  Created by Federico Pagnozzi on 28/11/14.
//  Copyright (c) 2014 Federico Pagnozzi. All rights reserved.
//  This file is distributed under the BSD 2-Clause License. See LICENSE.TXT
//  for details.

#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <ctime>
#include <cstring>
#include <algorithm>
#include "generalParser.h"
#include "examtt/examttparser.h"
#include "setup.h"

#include "examtt/examtt.h"
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


int main(int argc, const char *argv[])
{
    // prs::emili_header();

    /* initialize random seed: */
    srand ( time(0) );

    /* Create instance object */
    //testTaillardAccel();
    clock_t time = clock();

    // examTT test

    if(0){
        emili::ExamTT::test();
        double time_elapsed = (double)(clock() - time) / CLOCKS_PER_SEC;
        std::cout << "Time " << time_elapsed << std::endl;
        return 0;
    }


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
    
    prs::ExamTT::ExamTTParser p;
    // SAPFSPParser p;
    // SAQAPParser p;
    
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
    solution = ls->getBestSoFar();    
    // std::cout << "Number of delete a:" << emili::ExamTT::ExamTTSolution::numberOfDeletes << std::endl;
    double time_elapsed = (double)(clock()-time)/CLOCKS_PER_SEC;
    std::cout << "time : " << time_elapsed << std::endl;
    std::cout << "iteration counter : " << std::fixed << emili::iteration_counter()<< std::endl;
    std::cerr << std::fixed << solution->getSolutionValue() << std::endl;
    //cerr << time_elapsed << " ";    
    std::cout << "Objective function value: " << std::fixed << solution->getSolutionValue() << std::endl;
    std::cout << "Found solution: ";
    std::cout << std::fixed << solution->getSolutionRepresentation() << std::endl;
    std::cout << std::endl;
    std::cout << "Number of clones: " << emili::ExamTT::ExamTTSolution::numberOfClones << std::endl;
    // std::cout << "Number of delete: " << emili::ExamTT::ExamTTSolution::numberOfDeletes << std::endl;

    // std::cerr << std::fixed << solution->getSolutionValue() << endl;

    return 0;
}
