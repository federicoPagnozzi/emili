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
//#define MAIN_NEW
#include "pfsp/pfspBuilder.h"
//#include "template/problem_builder.h"
//#include "setup.h"
#include <sys/types.h>

#ifdef EM_LIB
#include <dirent.h>
#include <dlfcn.h>
#include <sstream>
#define SO_LIB ".so"
#define A_LIB ".a"
#endif

//#include "SA/sa_pfsp_parser.h"
//#include "SA/sa_qap_parser.h"
#include "SA/SABuilder.h"
#include "QAP/qapBuilder.h"
#include "stats/stats_builder.h"
#include "TSP/tspBuilder.h"

void g2c_info()
{
    std::cout << "usage in grammar2code mode : \n\tEMILI instance_file_path time random_seed" << std::endl;
    exit(0);
}

#ifndef MAIN_NEW
int main(int argc, char *argv[])
{
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
    float pls = 0;
    emili::LocalSearch* ls;
#include "algorithm.h"
#ifndef GRAMMAR2CODE
    prs::GeneralParser ps(argv,argc);
    ps.registerBuilder(&p);
    ls = ps.parseParams();
   // testHeuritstic(ps.getInstance());
    if(ls==nullptr)
    {
       // std::cout << "EXITING" << std::endl;
        exit(-1);
     //   return -1;
    }
    pls = ls->getSearchTime();//ps.ils_time;
#else
    pls = atoi(argv[2]);
    int seed = atoi(argv[3]);
    emili::initializeRandom(seed);
    time = clock();
#endif
    emili::Solution* solution;
    std::cout << "searching..." << std::endl;
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
    std::cout << "iteration counter : " << std::fixed << emili::iteration_counter()<< std::endl;
    std::cerr << std::fixed << solution->getSolutionValue() << std::endl;
    //cerr << time_elapsed << " ";    
    std::cout << "Objective function value: " << std::fixed << solution->getSolutionValue() << std::endl;
    std::cout << "Found solution: ";
    std::cout << std::fixed << solution->getSolutionRepresentation() << std::endl;
    std::cout << std::endl;

    // std::cerr << std::fixed << solution->getSolutionValue() << endl;
    
    return 0;
}
#else
#ifdef EM_LIB

typedef prs::Builder* (*getBuilderFcn)(prs::GeneralParserE* ge);

void loadBuilders(prs::GeneralParserE& ps)
{
   std::string so_ext(SO_LIB);
   std::string a_ext(A_LIB);
   const char* lib_dir = std::getenv("EMILI_LIBS");
   if(!lib_dir)
   {
       lib_dir = "./";
   }

    DIR* dp = opendir(lib_dir);
    dirent* den;
    if (dp != NULL){
       while (den = readdir(dp)){
        std::string file(den->d_name);
        bool load=false;
        if(file.size() > so_ext.size())
        {
          if(std::equal(file.begin() + file.size() - so_ext.size(), file.end(), so_ext.begin()))
          {
              load = true;
          }
          else if(std::equal(file.begin() + file.size() - a_ext.size(), file.end(), a_ext.begin()))
          {
              load = true;
          }
          if(load)
          {
             std::ostringstream oss;
             oss << lib_dir << "/" << file;
             void* lib = dlopen(oss.str().c_str(),RTLD_LAZY);
             getBuilderFcn *build = (getBuilderFcn*) dlsym(lib,"getBuilder");
             prs::Builder* bld = (*build)(&ps);
             ps.addBuilder(bld);
          }

        }
       }
    }

   /*else
   {
      std::cerr << "the EMILI_LIBS environmental variable is not set!" << std::endl;
      exit(-1);
   }*/
}
#endif

int main(int argc, char *argv[])
{
    prs::emili_header();
    srand ( time(0) );
    clock_t time = clock();
    if (argc < 3 )
    {
        prs::info();
        return 1;
    }
    float pls = 0;
    emili::LocalSearch* ls;

    prs::GeneralParserE  ps(argv,argc);
    prs::EmBaseBuilder emb(ps,ps.getTokenManager());
    prs::PfspBuilder pfspb(ps,ps.getTokenManager());
    prs::QAPBuilder qap(ps,ps.getTokenManager());
    prs::SABuilder sab(ps,ps.getTokenManager());
    prs::MABuilder mab(ps,ps.getTokenManager());
    prs::stats::StatsBuilder sb(ps,ps.getTokenManager());
    prs::TSPBuilder tsp(ps,ps.getTokenManager());
    //prs::problemX::ProblemXBuilder px(ps,ps.getTokenManager());
    ps.addBuilder(&emb);
    //ps.addBuilder(&px);
#ifdef EM_LIB
    loadBuilders(ps);
#else
    ps.addBuilder(&qap);
    ps.addBuilder(&pfspb);
    ps.addBuilder(&sab);
    ps.addBuilder(&mab);
    ps.addBuilder(&sb);
    ps.addBuilder(&tsp);
#endif
    ls = ps.parseParams();
    if(ls!=nullptr)
    {
        pls = ls->getSearchTime();//ps.ils_time;
        emili::Solution* solution;
        std::cout << "searching..." << std::endl;
        if(pls>0)
        {
            solution = ls->timedSearch(pls);
        }
        else
        {
            solution = ls->search();
        }
        if(!emili::get_print())
        {
            solution = ls->getBestSoFar();
            double time_elapsed = (double)(clock()-time)/CLOCKS_PER_SEC;
            double solval = solution->getSolutionValue();
            std::cout << "time : " << time_elapsed << std::endl;
            std::cout << "iteration counter : " << emili::iteration_counter()<< std::endl;
          //  std::cerr << solution->getSolutionValue() << std::endl;            
            std::cout << "Objective function value: "<< std::fixed << solval << std::endl;
            std::cerr << std::fixed << solval << std::endl;
          /*  std::cout << "Found solution: ";
            std::cout << solution->getSolutionRepresentation() << std::endl;
            std::cout << std::endl;*/
        }
        delete ls;
    }
}
#endif

