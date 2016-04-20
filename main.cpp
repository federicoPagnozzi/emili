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

#ifndef GRAMMAR2CODE

int main(int argc, const char *argv[])
{
    // prs::emili_header();

    srand(time(0)); // will probably changed by the parser

    clock_t time = clock();

    if(argc < 3)
    {
        prs::info();
        return 1;
    }

    prs::ExamTT::ExamTTParser p;

    prs::GeneralParser ps(argv,argc);
    ps.registerBuilder(&p);

    emili::LocalSearch* ls = nullptr;

    try {
        ls = ps.parseParams();
    } catch(prs::NoSearch) {
        std::cout << "No Search" << std::endl;
    }

    if(ls == nullptr) {
        return -1;
    }

    int pls = ls->getSearchTime(); // ps.ils_time;

    emili::Solution* returnedSolution = pls > 0 ? ls->timedSearch(pls) : ls->search();

    emili::Solution* solution = ls->getBestSoFar();

    using std::cout;
    using std::cerr;
    using std::endl;

    if(returnedSolution != solution) {
        cout << "Warning: " << "LocalSearch::search did not returned the same as BestSoFar" << endl;
    }

    cout
        << "time : " << ((double)(clock() - time)/CLOCKS_PER_SEC) << endl
        << "iteration counter : " << std::fixed << emili::iteration_counter()<< endl
        << "Objective function value: " << std::fixed << solution->getSolutionValue() << endl
        << "Found solution: " << std::fixed << solution->getSolutionRepresentation() << endl
        << "numberOfClones: " << emili::ExamTT::ExamTTSolution::numberOfClones << endl
        << "numberOfTotalCompute: " << emili::ExamTT::ExamTTSolution::numberOfTotalCompute << endl
    ;

    cerr << std::fixed << solution->getSolutionValue() << endl;

    return 0;
}

#else


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

    /* Read data from file */
    if (argc < 3 )
    {
        g2c_info();
        return 1;
    }
    // testNewEvaluationFunction(instance);
    // emili::pfsp::NWPFSP_MS problem(instance);
    // testHeuritstic(problem);


    pls = atoi(argv[2]);
    int seed = atoi(argv[3]);
    emili::initializeRandom(seed);
    time = clock();

    emili::Solution* returnedSolution = pls > 0 ? ls->timedSearch(pls) : ls->search();

    long int totalWeightedTardiness = problem.computeObjectiveFunction(sol);
    int njobs = problem.getNjobs();

    emili::Solution* solution = ls->getBestSoFar();

    using std::cout;
    using std::cerr;
    using std::endl;

    if(returnedSolution != solution)
        cout << "Warning: " << "LocalSearch::search did not returned the same as BestSoFar" << endl;

    cout
        << "time : " << ((double)(clock() - time)/CLOCKS_PER_SEC) << endl
        << "iteration counter : " << std::fixed << emili::iteration_counter()<< endl
        << "Objective function value: " << std::fixed << solution->getSolutionValue() << endl
        << "Found solution: " << std::fixed << solution->getSolutionRepresentation() << endl
        << "numberOfClones: " << emili::ExamTT::ExamTTSolution::numberOfClones << endl
        << "numberOfTotalCompute: " << emili::ExamTT::ExamTTSolution::numberOfTotalCompute << endl
    ;

    cerr << std::fixed << solution->getSolutionValue() << endl;

    return 0;
}

#endif
