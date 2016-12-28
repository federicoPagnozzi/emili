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

#include <execinfo.h>
#include <signal.h>
#include <stdlib.h>
#include <unistd.h>

void handler(int sig) {
    void *array[10];
    size_t size;

    // get void*'s for all entries on the stack
    size = backtrace(array, 10);

    // print out all the frames to stderr
    if(sig == SIGSEGV)
        std::cerr << "Error: SEGFAULT" << std::endl;
    else
        std::cerr << "Error: signal " << sig << std::endl;
    backtrace_symbols_fd(array, size, STDERR_FILENO);
    exit(1);
}

int main(int argc, const char *argv[])
{
    // prs::emili_header();

    signal(SIGSEGV, handler);

    srand(time(0)); // will probably be changed by the parser

    clock_t time = clock();

    if(argc < 3)
    {
        prs::info();
        return 1;
    }

    prs::ExamTT::ExamTTParser examParser;

    prs::GeneralParser ps(argv,argc);
    ps.registerBuilder(&examParser);

    emili::LocalSearch* ls = nullptr;

    try {
        ls = ps.parseParams();
    } catch(prs::NoSearch&) {
        std::cout << "No Search" << std::endl;
    } catch(prs::ErrorExpected& e) {
        std::cerr << e.what() << std::endl;
    } catch(prs::ParsingError& e) {
        std::cerr << "PARSING ERROR: " << e.what() << std::endl;
    }

    if(ls == nullptr) {
        return -1;
    }

    if(ps.noSearch) {
        std::cout << "No search" << std::endl;
        return 0;
    }

    int pls = ls->getSearchTime(); // ps.ils_time;

    ls->theInstance = &examParser.instance; // only used in finalise
    emili::Solution* solution = pls > 0 ? ls->timedSearch(pls) : ls->search();
    // if pls > 0 and solution times out, look at finalise function in emilibase.cpp

    emili::Solution* solutionBestSoFar = ls->getBestSoFar();

    using std::cout;
    using std::cerr;
    using std::endl;

    if(solutionBestSoFar != solution) {
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

    examParser.instance.finaliseSolution(solution); // will apply final hardweight
    cerr << std::fixed << solution->getSolutionValue() << endl; // to inform hook-run about the value found

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
