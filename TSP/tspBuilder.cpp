#include "tspBuilder.h"
#include "tspconstants.h"
#include <cstring>
#define TSPPROBLEMNAME "TSP"

bool prs::TSPBuilder::isCompatibleWith(char* problem_definition)
{
    if(strcmp(problem_definition,TSPPROBLEMNAME)==0)
    {
        return true;
    }
    return false;
}

bool prs::TSPBuilder::canOpenInstance(char* problem_definition)
{
    if(strcmp(problem_definition,TSPPROBLEMNAME)==0)
    {
        return true;
    }
    return false;
}

emili::Problem* prs::TSPBuilder::openInstance()
{
   tm.next();
   printTab("Travelling Salesman Problem");
   return new emili::tsp::TSP(tm.tokenAt(1));
}

emili::InitialSolution* prs::TSPBuilder::buildInitialSolution()
{
    prs::incrementTabLevel();
    emili::InitialSolution* init = nullptr;
    if(tm.checkToken(INITIAL_RANDOM))
    {
        emili::tsp::TSP* instance =(emili::tsp::TSP*) gp.getInstance();
        prs::printTab("Random initial solution");
        init = new emili::tsp::TSPRandomInitialSolution(*instance);
    }
    prs::decrementTabLevel();
    return init;
}

emili::Neighborhood* prs::TSPBuilder::buildNeighborhood()
{
    prs::incrementTabLevel();
    emili::Neighborhood* neigh=nullptr;
    emili::tsp::TSP* instance =(emili::tsp::TSP*) gp.getInstance();
    if(tm.checkToken(NEIGHBORHOOD_EXCHANGE)) {
        prs::printTab( "Exchange neighborhood");
        neigh = new emili::tsp::TSPExchangeNeighborhood(*instance);
    } else if (tm.checkToken(NEIGHBORHOOD_INSERT)) {
        prs::printTab( "Insert neighborhood");
        neigh = new emili::tsp::TSPInsertNeighborhood(*instance);
    }
    prs::decrementTabLevel();
    return neigh;
}

emili::Problem* prs::TSPBuilder::buildProblem(){
    return new emili::tsp::TSP(tm.tokenAt(1));
}
