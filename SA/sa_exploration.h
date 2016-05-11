#ifndef SA_EXPLORATION_H
#define SA_EXPLORATION_H


#include <string>
#include <iostream>
#include <vector>

#include "../emilibase.h"
#include "sa_acceptance_criteria.h"
#include "sa_termination_criteria.h"
#include "sa_common.h"



class SAExploration {

protected:
    std::string type;
    emili::Neighborhood* neigh;
    SAAcceptance* acceptance;
    SATermination* term;
    std::string tc_type;

public:
    SAExploration(emili::Neighborhood* _neigh,
                  SAAcceptance* _acceptance,
                  SATermination* _term,
                  std::string _type):
        neigh(_neigh),
        acceptance(_acceptance),
        term(_term),
        tc_type(_term->getType()),
        type(_type) { }

    virtual ~SAExploration() {}

    std::string getType(void) {
        return type;
    }

    /**
     * @brief nextSolution
     * @param startingSolution
     * @param status
     * @return next solution, may or may not be startingSolution
     */
    virtual emili::Solution* nextSolution(emili::Solution *startingSolution,
                                          SAStatus& status)=0;

}; // SAExploration


class SARandomExploration: public SAExploration {

public:
    SARandomExploration(emili::Neighborhood* _neigh,
                        SAAcceptance* _acceptance,
                        SATermination* _term):
        SAExploration(_neigh,
                      _acceptance,
                      _term,
                      SARANDOMEXPLORATION) { }

    virtual emili::Solution* nextSolution(emili::Solution *startingSolution,
                                          SAStatus& status);

}; // SARandomExploration


class SARandomExplorationNoCopy : public SAExploration {
public:
    SARandomExplorationNoCopy(emili::Neighborhood *_neigh, SAAcceptance *_acceptance, SATermination *_term)
        : SAExploration(_neigh, _acceptance, _term, "SARANDOMEXPLORATIONNOCOPY") {

    }

    // int nacc = 0, nnacc = 0;

    virtual emili::Solution* nextSolution(emili::Solution *startingSolution, SAStatus &status);
};

class SARandomExplorationNoCopyDebug : public SARandomExplorationNoCopy {
public:
    SARandomExplorationNoCopyDebug(emili::Neighborhood *_neigh, SAAcceptance *_acceptance, SATermination *_term) : SARandomExplorationNoCopy(_neigh, _acceptance, _term) {}

    int nacc = 0, nnacc = 0;

    virtual emili::Solution* nextSolution(emili::Solution *startingSolution, SAStatus &status);
};


class SASequentialExploration: public SAExploration {

public:
    SASequentialExploration(emili::Neighborhood* _neigh,
                            SAAcceptance* _acceptance,
                            SATermination* _term):
        SAExploration(_neigh,
                      _acceptance,
                      _term,
                      SASEQUENTIALEXPLORATION) { }

    virtual emili::Solution* nextSolution(emili::Solution *startingSolution,
                                          SAStatus& status);

}; // SASequentialExploration


class SAFirstImprovementExploration: public SAExploration {

public:
    SAFirstImprovementExploration(emili::Neighborhood* _neigh,
                                  SAAcceptance* _acceptance,
                                  SATermination* _term):
        SAExploration(_neigh,
                      _acceptance,
                      _term,
                      SAFIRSTIMPROVEMENTEXPLORATION) { }

    virtual emili::Solution* nextSolution(emili::Solution *startingSolution,
                                          SAStatus& status);

}; // SAFirstImprovementExploration


class SABEstImprovementExploration: public SAExploration {

public:
    SABEstImprovementExploration(emili::Neighborhood* _neigh,
                                 SAAcceptance* _acceptance,
                                 SATermination* _term):
        SAExploration(_neigh,
                      _acceptance,
                      _term,
                      SABESTIMPROVEMENTEXPLORATION) { }

    virtual emili::Solution* nextSolution(emili::Solution *startingSolution,
                                          SAStatus& status);

}; // SABEstImprovementExploration

#endif
