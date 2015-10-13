/*
 * =====================================================================================
 *
 *       Filename:  sa.cpp
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  06/07/2015 17:54:05
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *        Company:  
 *
 * =====================================================================================
 */

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

    std::string tc_type = terminationCriterion->getType();
    std::string ac_type = acceptanceCriterion->getType();
    std::string tl_type = tempLength->getType();

    acceptanceCriterion->setCurrentTemp(temp);

    neigh->reset();

    coolingScheme->setMaxIterations(tempLength->getLength());

    do {

        bestSoFar = exploration->nextSolution(bestSoFar, status);
        temp = coolingScheme->update_cooling(temp);
        acceptanceCriterion->setCurrentTemp(temp);

    } while(!terminationCriterion->terminate(status));

    return bestSoFar;
} // end search

void SimulatedAnnealing::reset(void) {
    temp = init_temp;
}
