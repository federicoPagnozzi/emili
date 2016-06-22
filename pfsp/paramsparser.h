//
//  Created by Federico Pagnozzi on 28/11/14.
//  Copyright (c) 2014 Federico Pagnozzi. All rights reserved.
//  This file is distributed under the BSD 2-Clause License. See LICENSE.TXT
//  for details.

#ifndef PARAMSPARSER_H
#define PARAMSPARSER_H
#include "../generalParser.h"
#include "permutationflowshop.h"
namespace prs
{
//void emili_header();
void info();

class ParamsParser: public AlgoBuilder
{
protected:    
    emili::pfsp::PermutationFlowShop* istance;
    std::vector< emili::pfsp::PermutationFlowShop* > istances;
    /* ALGOS */
    emili::LocalSearch* eparams(prs::TokenManager& tm);
    emili::LocalSearch* search(prs::TokenManager& tm);
    emili::LocalSearch* ils(prs::TokenManager& tm);
    emili::LocalSearch* gvns(prs::TokenManager& tm);
    emili::LocalSearch* ch6_params(prs::TokenManager& tm);
    emili::BestTabuSearch* tparams(prs::TokenManager& tm);
    emili::LocalSearch* vparams(prs::TokenManager& tm);
    void params(prs::TokenManager& tm);
    /*INITIAL SOLUTION*/
    emili::InitialSolution* init(prs::TokenManager& tm);
    /*TERMINATION*/
    emili::Termination* term(prs::TokenManager& tm);
    /*NEIGHBORHOOD*/
    emili::pfsp::PfspNeighborhood* neigh(prs::TokenManager& tm,bool checkExist);
    /*PERTURBATION*/
    emili::Perturbation* per(prs::TokenManager& tm);
    /*ACCEPTANCE*/
    emili::Acceptance* acc(prs::TokenManager& tm);
    /*TABU TENURE */
    emili::TabuMemory* tmemory(emili::pfsp::PfspNeighborhood* n,prs::TokenManager& tm);
    /*NEIGHBORHOOD UTILS*/
    void neighs(prs::TokenManager& tm);
    void neighs1(prs::TokenManager& tm);
    /*Problem load*/
    void problem(prs::TokenManager& tm);
    virtual std::string availableProblems() const;
public:
    virtual bool isParsable(std::string& problem) ;
    virtual emili::LocalSearch* buildAlgo(prs::TokenManager& tm);
    virtual std::string info();
    ParamsParser() { }
    ~ParamsParser() { delete istance; for(int i=0;i<istances.size();i++)delete istances[i];}
};
}
#endif // PARAMSPARSER_H
