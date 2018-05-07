#ifndef SA_EXPLORATION_H
#define SA_EXPLORATION_H


#include <string>
#include <iostream>
#include <vector>
#include <cmath>

#include "../emilibase.h"
#include "sa_common.h"
#include "sa_acceptance_criteria.h"
#include "sa_cooling.h"
#include "sa_termination_criteria.h"

namespace emili {
namespace sa {


class SAExploration {

protected:
    std::string type;
    emili::Neighborhood* neigh;
    SAAcceptance* acceptance;
    SACooling* cooling;
    SATermination* term;
    std::string tc_type;

public:
    SAExploration(std::string _type):
        type(_type) { }

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

    virtual void setNeighborhood(emili::Neighborhood* neighborhood)
    {
        neigh = neighborhood;
    }

    virtual void setAcceptance(emili::sa::SAAcceptance* accept)
    {
        acceptance = accept;
    }

    virtual void setCooling(emili::sa::SACooling* cool)
    {
        cooling = cool;
    }

    virtual void setTermination(emili::sa::SATermination* termination)
    {
        term = termination;
        tc_type = term->getType();
    }

}; // SAExploration


class SARandomExploration: public SAExploration {

public:
    SARandomExploration():
        SAExploration(SARANDOMEXPLORATION) { }

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
    SASequentialExploration():
        SAExploration(SASEQUENTIALEXPLORATION) { }

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
    SABestOfKExploration(long _k):
        k(_k),
        SAExploration(SABESTOFKEXPLORATION) { }

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
    SAFirstBestOfKExploration(long _k):
            k(_k),
            SAExploration(SAFIRSTBESTOFKEXPLORATION) { }

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


/**
 * Modified simulated annealing algorithms for the flow shop sequencing problem
 * Ishibuchi, Hisao and Misaki, Shinta and Tanaka, Hideo
 */
class SANSBestOfKExploration: public SAExploration {

long k;

public:
    SANSBestOfKExploration(long _k):
        k(_k),
        SAExploration(SANSBESTOFKEXPLORATION) { }

    SANSBestOfKExploration(emili::Neighborhood* _neigh,
                         SAAcceptance* _acceptance,
                         SACooling* _cooling,
                         SATermination* _term,
                         double _k):
        k(round(_k * _neigh->size())),
        SAExploration(_neigh,
                      _acceptance,
                      _cooling,
                      _term,
                      SANSBESTOFKEXPLORATION) { }

    virtual emili::Solution* nextSolution(emili::Solution *startingSolution,
                                          SAStatus& status);

}; // SANSBestOfKExploration

/**
 * Modified simulated annealing algorithms for the flow shop sequencing problem
 * Ishibuchi, Hisao and Misaki, Shinta and Tanaka, Hideo
 */
class SANSFirstBestOfKExploration: public SAExploration {

long k;

public:
    SANSFirstBestOfKExploration(long _k):
            k(_k),
            SAExploration(SANSFIRSTBESTOFKEXPLORATION) { }

    SANSFirstBestOfKExploration(emili::Neighborhood* _neigh,
                              SAAcceptance* _acceptance,
                              SACooling* _cooling,
                              SATermination* _term,
                              double _k):
        k(round(_k * _neigh->size())),
        SAExploration(_neigh,
                      _acceptance,
                      _cooling,
                      _term,
                      SANSFIRSTBESTOFKEXPLORATION) { }

    virtual emili::Solution* nextSolution(emili::Solution *startingSolution,
                                          SAStatus& status);

}; // SANSFirstBestOfKExploration

class SAFirstImprovementExploration: public SAExploration {

public:
    SAFirstImprovementExploration():
        SAExploration(SAFIRSTIMPROVEMENTEXPLORATION) { }

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
    SABEstImprovementExploration():
        SAExploration(SABESTIMPROVEMENTEXPLORATION) { }

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

}
}
#endif
