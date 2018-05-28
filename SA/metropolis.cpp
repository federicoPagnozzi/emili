#include "metropolis.h"

using namespace emili::sa;
using namespace emili::metropolis;

emili::Solution* MetropolisAlgorithm::search() {
    emili::Solution* current = this->init->generateSolution();
    emili::Solution* sol = MetropolisAlgorithm::search(current);
    return sol;
} // end search


emili::Solution* MetropolisAlgorithm::search(emili::Solution* initial) {
  emili::Solution* incumbent = this->init->generateEmptySolution();
  emili::Solution* tmpSol1;
  emili::Solution* tmpSol2;
  emili::Solution* accepted;
  *incumbent = *initial;
  bestSoFar  = incumbent;

  do {
    tmpSol1 = this->ls1->search(bestSoFar);
    tmpSol2 = this->ls2->search(tmpSol1);
    accepted = acceptance->accept(bestSoFar, tmpSol2);

    if (accepted == bestSoFar) {
      delete tmpSol2;
    } else {
      delete bestSoFar;
      bestSoFar = accepted;
    }

    delete tmpSol1;
  } while(!terminationCondition->terminate((SAStatus *)status));

  return this->getBestSoFar();
} // end search

