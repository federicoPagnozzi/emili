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

    int counter = 0;

    do {

        counter++;
        std::cout << "counter " << counter << std::endl;

        if(bestSoFar->operator > (*incumbent)) {
            delete bestSoFar;
            bestSoFar = incumbent;
        } else if (incumbent != bestSoFar) {
            delete incumbent;
        }

        std::cout << "generate random solution" << std::endl;
        incumbent = neighbh->random(bestSoFar);
        std::cout << "generated" << std::endl;

        incumbent = acceptanceCriterion->accept(bestSoFar, incumbent);
        std::cout << "accepted" << std::endl;

        temp = coolingScheme->update_cooling(temp);
        std::cout << "cooling updated" << std::endl;
        acceptanceCriterion->setCurrentTemp(temp);
        std::cout << "end of cycle" << std::endl;

    } while(!terminationCriterion->terminate(counter));

    return bestSoFar;
} // end search

void SimulatedAnnealing::reset(void) {
    temp = init_temp;
}