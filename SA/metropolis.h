#ifndef METROPOLIS_H
#define METROPOLIS_H

#include <string>
#include <iostream>
#include <vector>
#include <functional>

#include "sa.h"

#include "../emilibase.h"
namespace emili {
namespace metropolis {

class MetropolisAlgorithm: public AlternateLocalSearch {

private:
  emili::sa::SAAcceptance* acceptance;
  emili::sa::SATermination* terminationCondition;

public:
  MetropolisAlgorithm(InitialSolution& is,
                      emili::LocalSearch* localsearch1,
                      emili::LocalSearch* localsearch2,
                      emili::sa::SAAcceptance* acceptance,
                      emili::sa::SATermination* terminationCondition):
     acceptance(acceptance),
     terminationCondition(terminationCondition),
     emili::AlternateLocalSearch(is, localsearch1, localsearch2) { }

    emili::Solution* search();
    emili::Solution* search(emili::Solution* );

    virtual emili::Solution* getBestSoFar() { return status->getBestSolution();}

};

}

}
#endif
