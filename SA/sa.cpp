#include "sa.h"

emili::Solution* SimulatedAnnealing::search() {
    emili::Solution* current = init->generateSolution();
    emili::Solution* sol = search(current);
    if(current!=sol)
    delete current;

    return sol;
} // end search

emili::Solution* SimulatedAnnealing::search(emili::Solution* initial) {
    emili::Solution* incumbent = init->generateEmptySolution();
    // emili::Solution* accepted;
    *incumbent = *initial;
    bestSoFar  = incumbent;

    status->best = incumbent->clone();
    status->best_cost = incumbent->getSolutionValue();

    acceptanceCriterion->setCurrentTemp(temp);

    neigh->reset();

    do {

        bestSoFar = exploration->nextSolution(bestSoFar, *status);
        temp = coolingScheme->update_cooling(temp);
        acceptanceCriterion->setCurrentTemp(temp);

    } while(!terminationCriterion->terminate(*status));

    delete bestSoFar;

    return status->best;
} // end search

void SimulatedAnnealing::reset(void) {
    temp = init_temp;
}
