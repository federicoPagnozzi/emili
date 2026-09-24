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
    clock_t time = clock();
    if (argc < 3 )
    {
        prs::info();
        return 1;
    }
    float pls = 0;
    emili::LocalSearch* ls = nullptr;

    prs::GeneralParserE  ps(argv,argc);
    prs::EmBaseBuilder emb(ps,ps.getTokenManager());
    prs::PfspBuilder pfspb(ps,ps.getTokenManager());
    //prs::problemX::ProblemXBuilder px(ps,ps.getTokenManager());
    ps.addBuilder(&emb);
    //ps.addBuilder(&px);
#ifdef EM_LIB
    loadBuilders(ps);
#else
    ps.addBuilder(&pfspb);
#endif
    try
    {
        ls = ps.parseParams();
    }
    catch(std::exception& e)
    {
        std::cerr << e.what() << std::endl;
        return 255; // same process exit status exit(-1) produced
    }
    if(ls!=nullptr)
    {
        pls = ls->getSearchTime();//ps.ils_time;
        emili::Solution* searchResult;   // caller-owned clone returned by search()
        std::cout << "searching..." << std::endl;
        if(pls>0)
        {
            searchResult = ls->timedSearch(pls);
            emili::printFinalReport();
        }
        else
        {
            searchResult = ls->search();
        }
        if(!emili::get_print())
        {
            emili::Solution* solution = ls->getBestSoFar(); // component-owned, do NOT delete
            if(solution == nullptr)
            {   // e.g. a root that keeps no best (nols): report the caller-owned result
                solution = searchResult;
            }
            double time_elapsed = (double)(clock()-time)/CLOCKS_PER_SEC;
            double solval = solution->getSolutionValue();
            std::cout << "time : " << time_elapsed << std::endl;
            std::cout << "iteration counter : " << emili::iteration_counter()<< std::endl;
          //  std::cerr << solution->getSolutionValue() << std::endl;            
            std::cout << "Objective function value: "<< std::fixed << solval << std::endl;
            std::cerr << std::fixed << solval << std::endl;
            std::cout << "Found solution: ";
            std::cout << solution->getSolutionRepresentation() << std::endl;
            std::cout << std::endl;
        }
        delete searchResult; // caller owns the clone returned by search()
        // ls and every other component are owned by ps's ComponentRegistry
        // and destroyed, in reverse construction order, when ps goes out of
        // scope below (improvement plan P2.3). Do not delete ls here.
    }
}
