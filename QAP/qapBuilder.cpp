#include "qapBuilder.h"
#include <cstring>
#define QAPPROBLEMNAME "QAP"

/* QAP Initial Solution*/
#define INITIAL_RANDOM "random"

/* QAP neighborhoods*/
#define NEIGHBORHOOD_INSERT "insert"
#define NEIGHBORHOOD_EXCHANGE "exchange"
#define BEST2OPT_EXCHANGE "best2opt"
#define FIRST2OPT_EXCHANGE "first2opt"

bool prs::QAPBuilder::isCompatibleWith(char* problem_definition)
{
    if(strcmp(problem_definition,QAPPROBLEMNAME)==0)
    {
        return true;
    }
    return false;
}

bool prs::QAPBuilder::canOpenInstance(char* problem_definition)
{
    if(strcmp(problem_definition,QAPPROBLEMNAME)==0)
    {
        return true;
    }
    return false;
}

emili::Problem* prs::QAPBuilder::openInstance()
{
   tm.next();
   printTab("Quadratic Assignment Problem");
   return new emili::qap::QAP(tm.tokenAt(1));
}

emili::InitialSolution* prs::QAPBuilder::buildInitialSolution()
{
    prs::incrementTabLevel();
    emili::InitialSolution* init = nullptr;
    if(tm.checkToken(INITIAL_RANDOM))
    {
        emili::qap::QAP* instance =(emili::qap::QAP*) gp.getInstance();
        prs::printTab("Random initial solution");
        init = new emili::qap::QAPRandomInitialSolution(*instance);
    }
    prs::decrementTabLevel();
    return init;
}

emili::Neighborhood* prs::QAPBuilder::buildNeighborhood()
{
    prs::incrementTabLevel();
    emili::Neighborhood* neigh=nullptr;
    emili::qap::QAP* instance =(emili::qap::QAP*) gp.getInstance();
    if(tm.checkToken(NEIGHBORHOOD_EXCHANGE)) {
        prs::printTab( "Exchange neighborhood");
        neigh = new emili::qap::QAPExchangeNeighborhood(*instance);
    } else if (tm.checkToken(NEIGHBORHOOD_INSERT)) {
        prs::printTab( "Insert neighborhood");
        neigh = new emili::qap::QAPInsertNeighborhood(*instance);
    } else if (tm.checkToken(BEST2OPT_EXCHANGE)) {
        prs::printTab( "Best improvement 2-opt neighborhood");
        neigh = new emili::qap::QAPBest2optNeighborhood(*instance);
    } else if (tm.checkToken(FIRST2OPT_EXCHANGE)) {
        prs::printTab( "First improvement 2-opt neighborhood");
        neigh = new emili::qap::QAPFirst2optNeighborhood(*instance);
    }
    prs::decrementTabLevel();
    return neigh;
}

emili::Problem* prs::QAPBuilder::buildProblem(){
    return new emili::qap::QAP(tm.tokenAt(1));
}
