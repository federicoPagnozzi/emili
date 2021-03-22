#ifndef TSPBUILDER_H
#define TSPBUILDER_H

#include "../generalParser.h"
#include "tsp.h"
#include "tspinitialsolution.h"
#include "tspneighborhood.h"

#include <cstring>

#define TSPPROBLEMNAME "TSP"
namespace prs {

class TSPBuilder: public prs::Builder{
public:
    TSPBuilder(GeneralParserE& generalParser,TokenManager& tokenManager):Builder(generalParser,tokenManager) { }
    virtual bool isCompatibleWith(char* problem_definition);
    virtual bool canOpenInstance(char* problem_definition);
    virtual emili::Problem* openInstance();
    virtual emili::InitialSolution* buildInitialSolution();
    virtual emili::Neighborhood* buildNeighborhood();
    virtual emili::Problem* buildProblem();
};


}
#endif
