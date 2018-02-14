//
//  Created by Federico Pagnozzi on 12/12/17.
//  Copyright (c) 2014 Federico Pagnozzi. All rights reserved.
//  This file is distributed under the BSD 2-Clause License. See LICENSE.TXT
//  for details.

#ifndef ProblemXBUILDER_H
#define  ProblemXBUILDER_H
#include "../generalParser.h"
namespace prs
{
namespace problemX
{

class ProblemXBuilder: public Builder
{
public:
    ProblemXBuilder(GeneralParserE& generalParser,TokenManager& tokenManager):Builder(generalParser,tokenManager) { }
    virtual bool isCompatibleWith(char* problem_definition);
    virtual bool canOpenInstance(char* problem_definition);
    virtual emili::Problem* openInstance();    
    virtual emili::InitialSolution* buildInitialSolution();
    virtual emili::Neighborhood* buildNeighborhood();    
    virtual emili::Perturbation* buildPerturbation();
};
}
}
#endif //  PFSPBUILDER_H
