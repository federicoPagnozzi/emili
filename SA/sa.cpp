#include "sa.h"
using namespace emili::sa;
emili::Solution* SimulatedAnnealing::search() {
    emili::Solution* current = init->generateSolution();
    emili::Solution* sol = SimulatedAnnealing::search(current);
    if(current!=sol)
    delete current;

    return sol;
} // end search



emili::Solution* SimulatedAnnealing::search(emili::Solution* initial) {
    emili::Solution* incumbent = init->generateEmptySolution();
    // emili::Solution* accepted;
    *incumbent = *initial;
    bestSoFar  = incumbent;


    sastatus->best_cost = incumbent->getSolutionValue();
    sastatus->new_best_solution(incumbent->clone(),sastatus->best_cost,sastatus->temp);
    acceptanceCriterion->setCurrentTemp(sastatus->temp);

    neigh->reset();

    do {

        bestSoFar = exploration->nextSolution(bestSoFar, *sastatus);
        //status->temp = coolingScheme->update_cooling(status->temp);
        //acceptanceCriterion->setCurrentTemp(status->temp);

    } while(!terminationCriterion->terminate(*sastatus));

    delete bestSoFar;

    return sastatus->getBestSolution();
} // end search

void SimulatedAnnealing::reset(void) {
    temp = sastatus->init_temp;
    sastatus->temp = sastatus->init_temp;
}


void SimulatedAnnealing::setSearchTime(int _time) {
    if (sastatus->tc_type == NEVERTERM)
        coolingScheme->set_search_time(_time);

    this->seconds = _time;
}

