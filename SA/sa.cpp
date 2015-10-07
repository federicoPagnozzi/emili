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

    coolingScheme->setMaxIterations(tempLength->getLength());

    do {

        /*incumbent = neighbh->random(bestSoFar);
        accepted = acceptanceCriterion->accept(bestSoFar,
                                               incumbent);


        if (accepted == bestSoFar) {
            delete incumbent;
        } else {
            delete bestSoFar;
            if (tc_type == MAXBADITERS) {
                counter = 0;
            }
        }*/

        bestSoFar = exploration->nextSolution(bestSoFar, &counter);
        //std::cout << bestSoFar->getSolutionRepresentation() << std::endl;

        // bestSoFar = accepted;
         
        //std::cout << counter << std::endl;

        temp = coolingScheme->update_cooling(temp);
        acceptanceCriterion->setCurrentTemp(temp);

    } while(!terminationCriterion->terminate(counter));

    return bestSoFar;
} // end search

void SimulatedAnnealing::reset(void) {
    temp = init_temp;
}
