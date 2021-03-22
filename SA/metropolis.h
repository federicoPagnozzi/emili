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

    virtual emili::Solution* getBestSoFar() {
      emili::Solution*  sls1 = this->ls1->getBestSoFar();
      //std::cout << "in metropolis, SA1: " << std::fixed << sls1->getSolutionValue() << std::endl;

      emili::Solution*  sls2 = this->ls2->getBestSoFar();
      //std::cout << "in metropolis, SA2: " << std::fixed << sls2->getSolutionValue() << std::endl;

      emili::Solution*  sls = status->getBestSolution();
      //std::cout << "in metropolis: " << std::fixed << sls->getSolutionValue() << std::endl;

      // all these checks should not be necessary, but who knows...
      if (sls != nullptr) {
        if (!(sls1 == nullptr) && !(sls2 == nullptr)) {
          if (sls1->_vptr == sls2->_vptr && sls1->_vptr == sls->_vptr) {
            if (*sls1 <= *sls2 && *sls1 <= *sls) return sls1;
            if (*sls2 <= *sls1 && *sls2 <= *sls) return sls2;
            if (*sls <= *sls1 && *sls <= *sls2) return sls;
          } else {
            if (sls1->_vptr != sls2->_vptr) {
              if (sls1->_vptr == sls->_vptr) {
                if (*sls1 <= *sls) return sls1;
                else return sls;
              } else if (sls2->_vptr == sls->_vptr) {
                if (*sls2 <= *sls) return sls2;
                else return sls;
              } else { 
                return sls;
              }
            }
          }
        }
      } else {
        if (!(sls1 == nullptr) && (sls2 == nullptr)) {
          return sls1;
        }
        if ((sls1 == nullptr) && !(sls2 == nullptr)) {
          return sls2;
        }
        return nullptr;
      }
      return status->getBestSolution();
    }

};

}

}
#endif
