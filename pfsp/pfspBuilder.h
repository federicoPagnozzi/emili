//
//  Created by Federico Pagnozzi on 28/11/14.
//  Copyright (c) 2014 Federico Pagnozzi. All rights reserved.
//  This file is distributed under the BSD 2-Clause License. See LICENSE.TXT
//  for details.

#ifndef PFSPBUILDER_H
#define  PFSPBUILDER_H
#include "../generalParser.h"
#include "permutationflowshop.h"
namespace prs
{
void info();

class PfspBuilder: public Builder
{
public:
    PfspBuilder(GeneralParserE& generalParser,TokenManager& tokenManager):
        Builder(generalParser,tokenManager) { }
    virtual bool isCompatibleWith(char* problem_definition);
    virtual bool canOpenInstance(char* problem_definition);
    virtual emili::Problem* openInstance();
    virtual emili::LocalSearch* buildAlgo();
    virtual emili::InitialSolution* buildInitialSolution();
    virtual emili::Neighborhood* buildNeighborhood();
    virtual emili::Termination* buildTermination();
    virtual emili::Perturbation* buildPerturbation();
    virtual emili::Acceptance* buildAcceptance();
    virtual emili::TabuMemory* buildTabuTenure();
    virtual emili::Problem* buildProblem();    
};
}
#endif //  PFSPBUILDER_H
