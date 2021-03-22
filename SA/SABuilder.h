#ifndef SABUILDER_H
#define SABUILDER_H

#include "sa.h"
#include "metropolis.h"

#include "sa_constants.h"
#include "sa_init_temp.h"
#include "sa_acceptance_criteria.h"
#include "sa_cooling.h"
#include "sa_termination_criteria.h"
#include "sa_templength.h"
#include "sa_exploration.h"
#include "sa_temperature_restart.h"

#include "../generalParser.h"


#define COMPONENT_COOLING      0xD1
#define COMPONENT_TEMP_RESTART 0xD2
#define COMPONENT_TEMP_LENGTH  0xD3
#define COMPONENT_EXPLORATION  0xD4
#define COMPONENT_INIT_TEMP  0xD5
#define COMPONENT_SA_ACCEPTANCE 0xD6

namespace prs {


class SABuilder: public prs::Builder
{
public:
    SABuilder(GeneralParserE& generalParser,TokenManager& tokenManager):Builder(generalParser,tokenManager) { }
    virtual bool isCompatibleWith(char* problem_definition);
    virtual prs::Component buildComponent(int type);
    virtual emili::LocalSearch* buildAlgo();
    virtual emili::Termination* buildTermination();
    virtual emili::Acceptance* buildAcceptance();
    virtual emili::sa::SACooling* buildCooling();
    virtual emili::sa::SATempRestart* buildTempRestart();
    virtual emili::sa::SATempLength* buildTempLength();
    virtual emili::sa::SAExploration* buildExploration();
    virtual emili::sa::SAInitTemp* buildInitTemp();
};

class MABuilder: public prs::Builder
{
public:
    MABuilder(GeneralParserE& generalParser,TokenManager& tokenManager):Builder(generalParser,tokenManager) { }
    virtual bool isCompatibleWith(char* problem_definition) { return true; }
    virtual emili::LocalSearch* buildAlgo();
};


}
#endif
