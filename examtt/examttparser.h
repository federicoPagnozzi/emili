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
    std::ostringstream generic;
    void genericError(std::ostream&);
    void genericError(std::string name);
    void genericError(std::ostringstream&);
    void errorExpected(prs::TokenManager& tm, std::string name, std::vector<std::string> const& tokens);
protected:    
    //insert a variable to hold the problem instance
    emili::ExamTT::ExamTT instance;
    emili::LocalSearch* search(prs::TokenManager& tm, bool mustHaveInit, string prefix = "search");
    emili::BestTabuSearch* tparams(prs::TokenManager& tm);
    emili::TabuMemory* tmemory(emili::Neighborhood* n,prs::TokenManager& tm);
    emili::LocalSearch* vnd(prs::TokenManager& tm, bool mustHaveInit, string prefix = "vnd");

    emili::InitialSolution* initializer(prs::TokenManager& tm, string prefix = "initializer");
    emili::Termination* termination(prs::TokenManager& tm, string prefix = "termination");
    emili::Acceptance* acceptance(prs::TokenManager& tm, string prefix = "acceptance");
    emili::Perturbation* perturbation(prs::TokenManager& tm, string prefix = "perturbation");
    emili::Constructor* constructor(prs::TokenManager&, string prefix = "constructor");
    emili::Destructor* destructor(prs::TokenManager&, string prefix = "destructor");
    emili::ExamTT::InsertHeuristic* insertHeuristic(prs::TokenManager& tm, string prefix = "insertHeuristic");

    emili::Neighborhood* neigh(prs::TokenManager& tm, bool errorIfNotFound = true);
    std::vector<emili::Neighborhood *> neighs(prs::TokenManager& tm);

    void problem(prs::TokenManager& tm);

    struct SAParser : SAQAPParser {
        SimulatedAnnealing* buildSA(prs::TokenManager& tm, emili::InitialSolution* initsol, emili::Neighborhood* nei);
    } sa;

    virtual std::string availableProblems() const;
    emili::Problem* getInstance() override { return &instance; }
public:
    virtual bool isParsable(std::string &problem) ;
    virtual emili::LocalSearch* buildAlgo(prs::TokenManager& tm);
    virtual std::string info();
    ExamTTParser() { }
};
}
}
#endif // EXAMTT_PARAMSPARSER_H
