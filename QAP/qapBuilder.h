#ifndef QAPBUILDER_H
#define QAPBUILDER_H

#include "../generalParser.h"
#include "qap.h"
#include "qapinitialsolution.h"
#include "qapneighborhood.h"

#include <cstring>

#define QAPPROBLEMNAME "QAP"
namespace prs {

class QAPBuilder: public prs::Builder{
public:
    QAPBuilder(GeneralParserE& generalParser,TokenManager& tokenManager):Builder(generalParser,tokenManager) { }
    virtual bool isCompatibleWith(char* problem_definition);
    virtual bool canOpenInstance(char* problem_definition);
    virtual emili::Problem* openInstance();
    virtual emili::InitialSolution* buildInitialSolution();
    virtual emili::Neighborhood* buildNeighborhood();
    virtual emili::Problem* buildProblem();
};


}
#endif
