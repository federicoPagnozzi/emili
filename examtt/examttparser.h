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
public:
    emili::ExamTT::ExamTT instance;
protected:    
    //insert a variable to hold the problem instance
    std::string instanceFilename;
    emili::LocalSearch* search(prs::TokenManager& tm, bool mustHaveInit, string prefix = "search");
    emili::BestTabuSearch* tabu(prs::TokenManager& tm);
    emili::TabuMemory* tabuMemory(emili::Neighborhood* n,prs::TokenManager& tm);
    emili::LocalSearch* vnd(prs::TokenManager& tm, bool mustHaveInit, string prefix = "vnd");
    void globalParams(prs::TokenManager& tm);

    emili::InitialSolution* initializer(prs::TokenManager& tm, string prefix = "initializer");
    emili::Termination* termination(prs::TokenManager& tm, string prefix = "termination");
    emili::Acceptance* acceptance(prs::TokenManager& tm, string prefix = "acceptance");
    emili::Perturbation* perturbation(prs::TokenManager& tm, string prefix = "perturbation");
    emili::Constructor* constructor(prs::TokenManager&, string prefix = "constructor");
    emili::Destructor* destructor(prs::TokenManager&, string prefix = "destructor");
    emili::ExamTT::InsertHeuristic* insertHeuristic(prs::TokenManager& tm, string prefix = "insertHeuristic");

    emili::Neighborhood* neigh(prs::TokenManager& tm, bool errorIfNotFound = true, string prefix = "neighborhood");
    std::vector<emili::Neighborhood *> neighs(prs::TokenManager& tm);

    void problem(prs::TokenManager& tm);

    double getHardweightFromFeatures(bool USE_FORMULA_1 = true);

    struct SAParser : SAQAPParser {
        SimulatedAnnealing* buildSA(prs::TokenManager& tm, emili::InitialSolution* initsol, emili::Neighborhood* nei);
    } sa;

    std::string availableProblems() const override;
    std::vector<std::string> availableProblemsList() const override;
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
