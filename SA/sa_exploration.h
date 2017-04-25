#ifndef SA_EXPLORATION_H
#define SA_EXPLORATION_H


#include <string>
#include <iostream>
#include <vector>

#include "../emilibase.h"
#include "sa_common.h"
#include "sa_acceptance_criteria.h"
#include "sa_cooling.h"
#include "sa_termination_criteria.h"



class SAExploration {

protected:
    std::string type;
    emili::Neighborhood* neigh;
    SAAcceptance* acceptance;
    SACooling* cooling;
    SATermination* term;
    std::string tc_type;

public:
    SAExploration(emili::Neighborhood* _neigh,
                  SAAcceptance* _acceptance,
                  SACooling* _cooling,
                  SATermination* _term,
                  std::string _type):
        neigh(_neigh),
        acceptance(_acceptance),
        cooling(_cooling),
        term(_term),
        tc_type(_term->getType()),
        type(_type) { }

    std::string getType(void) {
        return type;
    }

    virtual emili::Solution* nextSolution(emili::Solution *startingSolution,
                                          SAStatus& status)=0;

}; // SAExploration


class SARandomExploration: public SAExploration {

public:
    SARandomExploration(emili::Neighborhood* _neigh,
                        SAAcceptance* _acceptance,
                        SACooling* _cooling,
                        SATermination* _term):
        SAExploration(_neigh,
                      _acceptance,
                      _cooling,
                      _term,
                      SARANDOMEXPLORATION) { }

    virtual emili::Solution* nextSolution(emili::Solution *startingSolution,
                                          SAStatus& status);

}; // SARandomExploration


class SASequentialExploration: public SAExploration {

public:
    SASequentialExploration(emili::Neighborhood* _neigh,
                            SAAcceptance* _acceptance,
                            SACooling* _cooling,
                            SATermination* _term):
        SAExploration(_neigh,
                      _acceptance,
                      _cooling,
                      _term,
                      SASEQUENTIALEXPLORATION) { }

    virtual emili::Solution* nextSolution(emili::Solution *startingSolution,
                                          SAStatus& status);

}; // SASequentialExploration


/**
 * Modified simulated annealing algorithms for the flow shop sequencing problem
 * Ishibuchi, Hisao and Misaki, Shinta and Tanaka, Hideo
 */
class SABestOfKExploration: public SAExploration {

long k;

public:
    SABestOfKExploration(emili::Neighborhood* _neigh,
                         SAAcceptance* _acceptance,
                         SACooling* _cooling,
                         SATermination* _term,
                         long _k):
        k(_k),
        SAExploration(_neigh,
                      _acceptance,
                      _cooling,
                      _term,
                      SABESTOFKEXPLORATION) { }

    virtual emili::Solution* nextSolution(emili::Solution *startingSolution,
                                          SAStatus& status);

}; // SABestOfKExploration

/**
 * Modified simulated annealing algorithms for the flow shop sequencing problem
 * Ishibuchi, Hisao and Misaki, Shinta and Tanaka, Hideo
 */
class SAFirstBestOfKExploration: public SAExploration {

long k;

public:
    SAFirstBestOfKExploration(emili::Neighborhood* _neigh,
                              SAAcceptance* _acceptance,
                              SACooling* _cooling,
                              SATermination* _term,
                              long _k):
        k(_k),
        SAExploration(_neigh,
                      _acceptance,
                      _cooling,
                      _term,
                      SAFIRSTBESTOFKEXPLORATION) { }

    virtual emili::Solution* nextSolution(emili::Solution *startingSolution,
                                          SAStatus& status);

}; // SAFirstBestOfKExploration


class SAFirstImprovementExploration: public SAExploration {

public:
    SAFirstImprovementExploration(emili::Neighborhood* _neigh,
                                  SAAcceptance* _acceptance,
                            SACooling* _cooling,
                                  SATermination* _term):
        SAExploration(_neigh,
                      _acceptance,
                      _cooling,
                      _term,
                      SAFIRSTIMPROVEMENTEXPLORATION) { }

    virtual emili::Solution* nextSolution(emili::Solution *startingSolution,
                                          SAStatus& status);

}; // SAFirstImprovementExploration


class SABEstImprovementExploration: public SAExploration {

public:
    SABEstImprovementExploration(emili::Neighborhood* _neigh,
                                 SAAcceptance* _acceptance,
                            SACooling* _cooling,
                                 SATermination* _term):
        SAExploration(_neigh,
                      _acceptance,
                      _cooling,
                      _term,
                      SABESTIMPROVEMENTEXPLORATION) { }

    virtual emili::Solution* nextSolution(emili::Solution *startingSolution,
                                          SAStatus& status);

}; // SABEstImprovementExploration

#endif
