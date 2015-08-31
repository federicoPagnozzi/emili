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
    *incumbent = *initial;
    bestSoFar  = incumbent;

    std::string type = terminationCriterion->getType();

    do {

        counter++;

        if(bestSoFar->operator > (*incumbent)) {
            delete bestSoFar;
            bestSoFar = incumbent;
            if (type == MAXBADITERS) {
                counter = 0;
            }
        } else if (incumbent != bestSoFar) {
            delete incumbent;
        }

        // std::cout << counter << std::endl;

        incumbent = neighbh->random(bestSoFar);

        incumbent = acceptanceCriterion->accept(bestSoFar,
                                                incumbent);

        temp = coolingScheme->update_cooling(temp);
        acceptanceCriterion->setCurrentTemp(temp);

    } while(!terminationCriterion->terminate(counter));

    return bestSoFar;
} // end search

void SimulatedAnnealing::reset(void) {
    temp = init_temp;
}
