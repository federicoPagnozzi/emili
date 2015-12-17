#ifndef EXAMTT_PARAMSPARSER_H
#define EXAMTT_PARAMSPARSER_H
#include "../generalParser.h"
#include "examtt.h"

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
    /* Called by buildAlgo it calls all the other methods.
     * Add here the call to your metaheuristcs if you don't want
     * to let the user include it in the ILS.
    */
    emili::LocalSearch* eparams(prs::TokenManager& tm);
    /* This method load all the metaheuristics that
     * can be included in a ILS.
     */
    emili::LocalSearch* search(prs::TokenManager& tm);
    /*This method load a ILS and it's components.
     * It calls search() to build it's inner metaheuristic
    */
    emili::LocalSearch* ils(prs::TokenManager& tm);
    /* Even though this method returns a BestTabuSearch it's able to
       load best and first tabu search.
     */
    emili::BestTabuSearch* tparams(prs::TokenManager& tm);
    /* This method should be implemented to load the tabutenure objects
     * for the problem that is supported by this parser
     */
    emili::TabuMemory* tmemory(emili::Neighborhood* n,prs::TokenManager& tm);
    /* This method loads the parameters for the localsearch*/
    void params(prs::TokenManager& tm);
    /* This method loads the parameters for the VND*/
    emili::LocalSearch* vparams(prs::TokenManager& tm);
    /*This method should be implemented to load the initial Solution objects
     * for the problem that is supported by this parser
     * */
    emili::InitialSolution* init(prs::TokenManager& tm);
    /*This method loads the termination criteria*/
    emili::Termination* term(prs::TokenManager& tm);
    /*This method loads the acceptance criteria*/
    emili::Acceptance* acc(prs::TokenManager& tm);
    /*This method should be implemented to load the perturbations
     * for the problem that is supported by this parser
     * */
    emili::Perturbation* per(prs::TokenManager& tm);

    emili::Neighborhood* neigh(prs::TokenManager& tm);
    emili::Neighborhood* neighV(prs::TokenManager& tm);
    void neighs(prs::TokenManager& tm);
    void neighs1(prs::TokenManager& tm);
   /*This method loads the problem instance*/
    void problem(prs::TokenManager& tm);

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
