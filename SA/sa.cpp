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

    status->best = incumbent->clone();
    status->best_cost = incumbent->getSolutionValue();

    acceptanceCriterion->setCurrentTemp(status->temp);

    neigh->reset();

    do {

        bestSoFar = exploration->nextSolution(bestSoFar, *status);
        //status->temp = coolingScheme->update_cooling(status->temp);
        //acceptanceCriterion->setCurrentTemp(status->temp);

    } while(!terminationCriterion->terminate(*status));

    delete bestSoFar;

    return status->best;
} // end search

void SimulatedAnnealing::reset(void) {
    temp = status->init_temp;
    status->temp = status->init_temp;
}


void SimulatedAnnealing::setSearchTime(int _time) {
    if (status->tc_type == NEVERTERM)
        coolingScheme->set_search_time(_time);

    this->seconds = _time;
}

void set_status(SAStatus* _status) {
    // status = _status;

    // status->set_types(terminationCriterion->getType(),
    //                   acceptanceCriterion->getType(),
    //                   tempLength->getType(),
    //                   temprestart->getType());
    // status->init_temp = init_temp;
    // status->final_temp = initialTemperature->getMinTemp();
    // status->temp = init_temp;
    // acceptanceCriterion->setCurrentTemp(status->temp);
    

    // /**
    //  * initialization of attribute depends on termination criteria
    //  * but in sa_termination_criteria.h I have to include sa_common.h
    //  * therefore it sucks a bit but I have to initialize this here.
    //  */
    // status->tenure = 1;

    // if (status->tc_type == LASTACCRATETERM) {
    //   status->tenure = terminationCriterion->getTenure();
    // } else if (status->tr_type == SALASTRATERESTART        ||
    //            status->tr_type == SALASTRATEREHEAT         ||
    //            status->tr_type == SALOCALMINENHANCEDREHEAT   ) {
    //   status->tenure = temprestart->getTenure();
    // }

    // status->last_accepted = (short *)
    //     malloc(status->tenure * sizeof(short));

    // // try at least status->tenure solutions
    // // otherwise it will terminate immediately
    // for (int i = 0 ; i < status->tenure ; i++) {
    //   status->last_accepted[i] = 1;
    // }

    // status->final_temp = initialTemperature->getMinTemp();
    // status->init_prob = initialTemperature->getInit_prob();
    // status->neigh_size = neighborhood->size();

    // acceptanceCriterion->set_status(status);
    // temprestart->set_status(status);
    // tempLength->set_status(status);
    // coolingScheme->set_status(status);
    // initialTemperature->set_status(status);
}
