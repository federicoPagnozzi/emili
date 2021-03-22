#include "sa.h"
using namespace emili::sa;
emili::Solution* SimulatedAnnealing::search() {
    emili::Solution* current = init->generateSolution();
    emili::Solution* sol = SimulatedAnnealing::search(current);
    //if(current!=sol)
    //delete current;

    return sol;
} // end search



emili::Solution* SimulatedAnnealing::search(emili::Solution* initial) {
    reset();
    emili::Solution* incumbent = init->generateEmptySolution();
    //emili::Solution* accepted;
    incumbent = initial;
    bestSoFar = incumbent;
    sastatus->init_cost = initial->getSolutionValue();

    //printf("\nSA::search INITIAL SOLUTION %f\n", initial->getSolutionValue());
    //sastatus->best_cost = incumbent->getSolutionValue();
    /**/if (1 || sastatus->best == nullptr ||
       initial->getSolutionValue() < sastatus->best_cost) {
      //printf("\n\n\nFIRST time %f\n", incumbent->getSolutionValue());
      sastatus->new_best_solution_silent(initial, initial->getSolutionValue(), sastatus->temp);
    } else {
      //printf("best is %f, incumbent is %f\n", sastatus->best_cost,  incumbent->getSolutionValue());
    }/**/
    acceptanceCriterion->setCurrentTemp(sastatus->temp);

    neigh->reset();

    do {

        emili::Solution* tempSolSA = exploration->nextSolution(bestSoFar, *sastatus);
        if (tempSolSA != bestSoFar) {
          /** /
          //delete bestSoFar;
          *bestSoFar = *tempSolSA;
          / **/
          /**/
          delete bestSoFar;
          bestSoFar = tempSolSA;
          /**/
          /*std::cout << bestSoFar << std::endl;
          std::cout << std::fixed << bestSoFar->getSolutionValue() << std::endl;*/
        } else {
        }
      /*  printf("%s\n", sastatus->best->getSolutionRepresentation().c_str());
        printf("%f\n", sastatus->best_cost);
        printf("%f\n", bestSoFar->getSolutionValue());*/
        //printf("%f \n", sastatus->getBestSolution()->getSolutionValue());
        //printf("%f \n", getSearchStatus().getBestSolution()->getSolutionValue());
        //status->temp = coolingScheme->update_cooling(status->temp);
        //acceptanceCriterion->setCurrentTemp(status->temp);
        //std::cout << std::fixed << "TEMPSTATS " << sastatus->total_counter << " " << sastatus->temp << " " << bestSoFar->getSolutionValue() << std::endl;
    } while(!terminationCriterion->terminate(*sastatus));
      /*  printf("%f\n", bestSoFar->getSolutionValue());
        printf("%s\n", sastatus->best->getSolutionRepresentation().c_str());
        printf("%f\n", sastatus->best_cost);*/

    //delete bestSoFar;
    //sastatus->print();
    //printf("returning from SA %f\n",  sastatus->best_cost);
    //printf("returning from SA %f\n",  sastatus->getBestSolution()->getSolutionValue());
    return sastatus->getBestSolution();
} // end search

void SimulatedAnnealing::reset(void) {
    sastatus->resetCounters();
    temp = sastatus->init_temp;
    sastatus->temp = sastatus->init_temp;
}


void SimulatedAnnealing::setSearchTime(int _time) {
    if (sastatus->tc_type == NEVERTERM)
        coolingScheme->set_search_time(_time);

    this->seconds = _time;
}

emili::Solution* SimulatedAnnealingIncumbent::search(emili::Solution* initial) {
    reset();
    emili::Solution* incumbent = init->generateEmptySolution();
    incumbent = initial;
    bestSoFar = incumbent;

    //std::cout << initial->getSolutionRepresentation() << std::endl;
    if (sastatus->best == nullptr ||
       initial->getSolutionValue() < sastatus->best_cost) {
      //sastatus->new_best_solution(initial, initial->getSolutionValue(), sastatus->temp);
      sastatus->new_best_solution_silent(initial, initial->getSolutionValue(), sastatus->temp);
    }
    //std::cout << "after set" << std::endl;
    //sastatus->new_best_solution_silent(initial, initial->getSolutionValue(), sastatus->temp);
    acceptanceCriterion->setCurrentTemp(sastatus->temp);
    //std::cout << initial->getSolutionValue() << std::endl;
    neigh->reset();

    do {
        emili::Solution* tempSolSA = exploration->nextSolution(bestSoFar, *sastatus);
        if (tempSolSA != bestSoFar) {
          *bestSoFar = *tempSolSA;
        }
    } while(!terminationCriterion->terminate(*sastatus));
    return bestSoFar;
} // end search

