#ifndef EXAMTT_PARAMSPARSER_H
#define EXAMTT_PARAMSPARSER_H
#include "../generalParser.h"
#include "examtt.h"

#include <utility>

#include "../SA/sa_qap_parser.h"

namespace prs
{
namespace  ExamTT
{

/*
    This class models a parser to instantiate an algorithm from the arguments.
    Every method that loads a component will(should) terminate the execution of
    EMILI with errors if one of the needed parameters is missing.
*/
class ExamTTParser: public AlgoBuilder
{
protected:
    void errorExpected(prs::TokenManager& tm, std::string name, std::vector<std::string> const& tokens);
protected:    
    //insert a variable to hold the problem instance
    emili::ExamTT::ExamTT instance;
    emili::LocalSearch* eparams(prs::TokenManager& tm);
    emili::LocalSearch* search(prs::TokenManager& tm);
    emili::LocalSearch* ils(prs::TokenManager& tm);
    emili::BestTabuSearch* tparams(prs::TokenManager& tm);
    emili::TabuMemory* tmemory(emili::Neighborhood* n,prs::TokenManager& tm);
    std::tuple<emili::InitialSolution*, emili::Termination*, emili::Neighborhood*> params(prs::TokenManager& tm);
    emili::LocalSearch* vparams(prs::TokenManager& tm);
    emili::InitialSolution* init(prs::TokenManager& tm);
    emili::Termination* term(prs::TokenManager& tm);
    emili::Acceptance* acc(prs::TokenManager& tm);
    emili::Perturbation* per(prs::TokenManager& tm);

    emili::Neighborhood* neigh(prs::TokenManager& tm);
    emili::Neighborhood* neighV(prs::TokenManager& tm);
    std::vector<emili::Neighborhood *> neighs(prs::TokenManager& tm);

    void neighs1(prs::TokenManager& tm, std::vector<emili::Neighborhood *> &nes);
    void problem(prs::TokenManager& tm);

    struct SAParser : SAQAPParser {
        SimulatedAnnealing* buildSA(prs::TokenManager& tm, emili::InitialSolution* initsol, emili::Neighborhood* nei);
    } sa;

    virtual std::string availableProblems() const;
public:
    virtual bool isParsable(std::string &problem) ;
    virtual emili::LocalSearch* buildAlgo(prs::TokenManager& tm);
    virtual std::string info();
    ExamTTParser() { }
};
}
}
#endif // EXAMTT_PARAMSPARSER_H
