#ifndef PARAMSPARSER_H
#define PARAMSPARSER_H
#include "emilibase.h"
#include "permutationflowshop.h"

namespace prs
{
void emili_header();
void info();

class ParamsParser
{
protected:
    emili::pfsp::PermutationFlowShop* istance;
    char** tokens;
    int numberOfTokens;
    int currentToken;
    char* nextToken();
    emili::LocalSearch* eparams();
    emili::LocalSearch* search();
    emili::LocalSearch* ils();
    emili::LocalSearch* ig();
    emili::LocalSearch* gvns();
    int ilstime();
    emili::TabuSearch* tparams();
    emili::TabuMemory* tmemory(emili::pfsp::PfspNeighborhood* n);
    std::pair<int,int> tsettings();
    int ttsize();
    int ttiter();
    void params();
    int getSeed();
    emili::LocalSearch* vparams();
    emili::InitialSolution* init();
    emili::Termination* term();
    emili::Acceptance* acc();
    emili::Perturbation* per();
    emili::pfsp::PfspNeighborhood* neigh();

    emili::pfsp::PfspNeighborhood* neighV();
    void neighs();
    void neighs1();
    void problem();
    int number();

    float decimal();
public:
    int ils_time;
    ParamsParser(char** tokens,int numberOfTokens,emili::pfsp::PermutationFlowShop* is):tokens(tokens),numberOfTokens(numberOfTokens),currentToken(2),istance(is),ils_time(-123) { }
    ParamsParser(char** tokens,int numberOfTokens):tokens(tokens),numberOfTokens(numberOfTokens),currentToken(2),ils_time(-123) { }
    emili::LocalSearch* parseParams();
    emili::pfsp::PermutationFlowShop& getInstance() { return *istance;}

};
}
#endif // PARAMSPARSER_H
