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
    emili::LocalSearch* eparams(prs::TokenManager& tm);
    emili::LocalSearch* search(prs::TokenManager& tm);
    emili::LocalSearch* ils(prs::TokenManager& tm);
    emili::LocalSearch* gvns(prs::TokenManager& tm);
    emili::BestTabuSearch* tparams(prs::TokenManager& tm);
    emili::TabuMemory* tmemory(emili::pfsp::PfspNeighborhood* n,prs::TokenManager& tm);
    void params(prs::TokenManager& tm);
    emili::LocalSearch* vparams(prs::TokenManager& tm);
    emili::InitialSolution* init(prs::TokenManager& tm);
    emili::Termination* term(prs::TokenManager& tm);
    emili::Acceptance* acc(prs::TokenManager& tm);
    emili::Perturbation* per(prs::TokenManager& tm);
    emili::pfsp::PfspNeighborhood* neigh(prs::TokenManager& tm,bool checkExist);
    void neighs(prs::TokenManager& tm);
    void neighs1(prs::TokenManager& tm);
    void problem(prs::TokenManager& tm);

    virtual std::string availableProblems() const;
public:
    virtual bool isParsable(std::string& problem) ;
    virtual emili::LocalSearch* buildAlgo(prs::TokenManager& tm);
    virtual std::string info();
    ParamsParser() { }
};
}
#endif // PARAMSPARSER_H
