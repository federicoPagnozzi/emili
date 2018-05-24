//
//  Created by Federico Pagnozzi on 28/11/14.
//  Copyright (c) 2014 Federico Pagnozzi. All rights reserved.
//  This file is distributed under the BSD 2-Clause License. See LICENSE.TXT
//  for details.

#ifndef PARAMSPARSER_H
#define PARAMSPARSER_H

#include "../SA/sa.h"

#include "../SA/sa_constants.h"
#include "../SA/sa_init_temp.h"
#include "../SA/sa_acceptance_criteria.h"
#include "../SA/sa_cooling.h"
#include "../SA/sa_termination_criteria.h"
#include "../SA/sa_templength.h"
#include "../SA/sa_exploration.h"
#include "../SA/sa_temperature_restart.h"

#include "../generalParser.h"
#include "permutationflowshop.h"

using namespace emili::sa;
namespace prs
{
//void emili_header();
void info();

class ParamsParser: public AlgoBuilder
{
protected:    
    emili::pfsp::PermutationFlowShop* instance;
    std::vector< emili::pfsp::PermutationFlowShop* > istances;
    /**  ALGOS */
    emili::LocalSearch* eparams(prs::TokenManager& tm);
    emili::LocalSearch* search(prs::TokenManager& tm);
    emili::LocalSearch* ils(prs::TokenManager& tm);
    emili::LocalSearch* gvns(prs::TokenManager& tm);
    emili::LocalSearch* ch6_params(prs::TokenManager& tm);
    emili::BestTabuSearch* tparams(prs::TokenManager& tm);
    emili::LocalSearch* vparams(prs::TokenManager& tm);
    void params(prs::TokenManager& tm);
    /** INITIAL SOLUTION*/
    emili::InitialSolution* init(prs::TokenManager& tm);
    /** TERMINATION*/
    emili::Termination* term(prs::TokenManager& tm);
    /** NEIGHBORHOOD*/
    emili::Neighborhood* neigh(prs::TokenManager& tm,bool checkExist);
    /** PERTURBATION*/
    emili::Perturbation* per(prs::TokenManager& tm);
    /** ACCEPTANCE*/
    emili::Acceptance* acc(prs::TokenManager& tm);
    /** TABU TENURE */
    emili::TabuMemory* tmemory(emili::pfsp::PfspNeighborhood* n,prs::TokenManager& tm);
    /** NEIGHBORHOOD UTILS*/
    void neighs(prs::TokenManager& tm);
    void neighs1(prs::TokenManager& tm);
    /** Problem load*/
    void problem(prs::TokenManager& tm);
    virtual std::string availableProblems() const;
public:
    virtual bool isParsable(std::string& problem) ;
    virtual emili::LocalSearch* buildAlgo(prs::TokenManager& tm);
    virtual std::string info();
    ParamsParser() { }
    ~ParamsParser() { delete instance; for(int i=0;i<istances.size();i++)delete istances[i];}
/*
        / **
     * identify cooling scheme
     * @param  tm TokenManager
     * @return    SACooling object
     * /
    SACooling*       COOL(prs::TokenManager& tm,
                              SAInitTemp *it,
                              emili::Neighborhood *nei,
                              emili::pfsp::PermutationFlowShop *instance);

    / **
     * identify acceptance criterion
     * @param  tm TokenManager
     * @return    SAAcceptance object
     * /
    SAAcceptance*    ACCEPTANCE(prs::TokenManager& tm,
                                SAInitTemp *inittemp,
                                 emili::Neighborhood *nei,
                                 emili::Problem* instance);

    / **
     * identify termination criterion
     * @param  tm TokenManager
     * @return    SATermination object
     * /
    SATermination*   TERMINATION(prs::TokenManager& tm,
                                 SAInitTemp* inittemp,
                                 emili::Neighborhood *nei);

    / **
     * identify Neighborhood
     * @param  tm TokenManager
     * @return    Neighborhood object
     * /
    emili::Neighborhood*  NEIGH(prs::TokenManager& tm);

    / **
     * identify initial temperature
     * @param  tm      TokenManager
     * @param  initsol initial solution
     * @return         InitTemp object
     * /
    SAInitTemp*      INITTEMP(prs::TokenManager&      tm,
                              emili::InitialSolution* initsol,
                              emili::Neighborhood *nei,
                              emili::pfsp::PermutationFlowShop *instance);

    /  **
     * identify initial solution builder
     * @param  tm TokenManager
     * @return    InitialSolution object
     * /
    emili::InitialSolution* INITSOL(prs::TokenManager& tm);


    SATempLength* TEMPLENGTH(prs::TokenManager& tm,
                             emili::Neighborhood* neigh,
                             emili::Problem* instance);

    SAExploration* EXPLORATION(prs::TokenManager& tm,
                                        emili::Neighborhood* neigh,
                                        SAAcceptance *acc,
                                        SACooling *cool,
                                        SATermination *term);

    SATempRestart *TEMPRESTART(prs::TokenManager& tm,
                               SAInitTemp *it,
                               emili::Neighborhood* neigh);*/
};
}
#endif // PARAMSPARSER_H
