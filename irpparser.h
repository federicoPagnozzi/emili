#ifndef IRPPARSER_H
#define IRPPARSER_H
#include "generalParser.h"
#include "irp.h"

#define PROBLEM_IRP "IRP"
namespace prs
{
namespace irp
{

/* This class should do the parsing of the arguments to build an algorithm for the IRP*/
class IrpParser: public AlgoBuilder
{
protected:
    /*Should return all the different problem for which this class can build an algorithm*/
    virtual std::string availableProblems() const;
public:
    /*This method tells the general parser if this class can parse the arguments for problem*/
     virtual bool isParsable(std::string& problem);
    /*This method should return human readable informations to let the user know how to instantiate
      an algorithm for this problem : which kind of metaheuristics are supported (ILS, TABU, SA etc...),
      which kind of neighborhood, termination criteria, acceptance , perturbations and so on...
    */
    virtual std::string info();

    /*
        This method is called by the general parser to build the algorithm.
        It must use the TokenManager to parse the arguments and build the algorithm accordigly.
    */

    virtual emili::LocalSearch* buildAlgo(prs::TokenManager& tm);
};


}
}

#endif // IRPPARSER_H
