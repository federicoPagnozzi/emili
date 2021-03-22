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

  if (status->best == nullptr ||
       incumbent->getSolutionValue() < status->best->getSolutionValue()) {
    printf("FIRST time metropolis %f\n", incumbent->getSolutionValue());
    status->newBestSolution(incumbent->clone());
  }

  ((emili::sa::SimulatedAnnealing *)this->ls1)->sastatus->resetCounters();
  ((emili::sa::SimulatedAnnealing *)this->ls1)->sastatus->resetCounters();

  do {

   // printf("bestsofar %f\n\n", bestSoFar->getSolutionValue());

    //printf("\n\n\nVAFFANCULOOOOOOO before ls1\n");
/*    ((SAStatus &)this->ls1->getSearchStatus()).new_best_solution(status->getBestSolution(),
                                                                status->getBestSolution()->getSolutionValue(),
                                                                0);*/
    ((emili::sa::SimulatedAnnealing *)this->ls1)->sastatus->new_best_solution_silent(bestSoFar, bestSoFar->getSolutionValue(), 0);
    //((emili::sa::SimulatedAnnealing *)this->ls1)->sastatus->best_cost = bestSoFar->getSolutionValue();
    tmpSol1 = this->ls1->search(bestSoFar);
    status->counter += ((emili::sa::SimulatedAnnealing *)this->ls1)->sastatus->counter;
    status->total_counter +=  ((emili::sa::SimulatedAnnealing *)this->ls1)->sastatus->local_counter;
    status->not_improved +=  ((emili::sa::SimulatedAnnealing *)this->ls1)->sastatus->not_improved;

    //printf("this status %f\n", status->best->getSolutionValue());
    //printf("ls1 status %f\n",  ((emili::sa::SimulatedAnnealing *)this->ls1)->sastatus->best_cost);

    /*if (((emili::sa::SimulatedAnnealing *)this->ls1)->sastatus->best_cost < status->best->getSolutionValue()) {
      status->newBestSolution(((emili::sa::SimulatedAnnealing *)this->ls1)->sastatus->best);
      printf("update metropolis after 1 %f\n", status->best->getSolutionValue());
    }*/

    ((emili::sa::SimulatedAnnealing *)this->ls1)->sastatus->resetCounters(status->total_counter);
    ((emili::sa::SimulatedAnnealing *)this->ls2)->sastatus->resetCounters(status->total_counter);

    //printf("\n\n\n ls1 to ls2\n");

    ((emili::sa::SimulatedAnnealing *)this->ls2)->sastatus->best_cost = ((emili::sa::SimulatedAnnealing *)this->ls1)->sastatus->best_cost;
    status->best = ((emili::sa::SimulatedAnnealing *)this->ls1)->sastatus->best->clone();

    //((emili::sa::SimulatedAnnealing *)this->ls2)->sastatus->best = tmpSol1;
    //((emili::sa::SimulatedAnnealing *)this->ls2)->sastatus->best_cost = tmpSol1->getSolutionValue();
    ((emili::sa::SimulatedAnnealing *)this->ls2)->sastatus->new_best_solution_silent(status->best, status->best->getSolutionValue(), 0);

    bestSoFar = this->ls2->search(tmpSol1);

    status->counter += ((emili::sa::SimulatedAnnealing *)this->ls2)->sastatus->counter;
    status->total_counter +=  ((emili::sa::SimulatedAnnealing *)this->ls2)->sastatus->local_counter;
    status->not_improved +=  ((emili::sa::SimulatedAnnealing *)this->ls2)->sastatus->not_improved;

    //printf("this status %f\n", status->best->getSolutionValue());
    //printf("ls2 status %f\n",  ((emili::sa::SimulatedAnnealing *)this->ls2)->sastatus->best_cost);

    /*if (((emili::sa::SimulatedAnnealing *)this->ls2)->sastatus->best_cost < status->best->getSolutionValue()) {
      status->newBestSolution(((emili::sa::SimulatedAnnealing *)this->ls2)->sastatus->best);
      printf("update metropolis after 2 %f\n", status->best->getSolutionValue());
    }*/
    ((emili::sa::SimulatedAnnealing *)this->ls2)->sastatus->resetCounters(status->total_counter);

    ((emili::sa::SimulatedAnnealing *)this->ls1)->sastatus->best_cost = ((emili::sa::SimulatedAnnealing *)this->ls2)->sastatus->best_cost;
    status->best = ((emili::sa::SimulatedAnnealing *)this->ls2)->sastatus->best->clone();

    //printf("\n\n\nls2 done\n");
    //((SAStatus &)this->ls2->getSearchStatus()).print();
    //accepted = acceptance->accept(bestSoFar, tmpSol2);
    //((SAStatus &)this->ls2->getSearchStatus()).print();

    /*if (accepted == bestSoFar) {
      //delete tmpSol2;
    } else {
      //delete bestSoFar;
      bestSoFar = accepted;
    }*/
    //bestSoFar = tmpSol2;
    //this->changeTurn();

    /*if (tmpSol1 == nullptr)
      printf("oh fuck\n");
    if (tmpSol1 != bestSoFar)*/
    //  delete tmpSol1;
    //delete bestSoFar;
    //
    //printf("bestsofar at the end %f\n\n", bestSoFar->getSolutionValue());

    //printf("\n\nthis status at the end: %f\n\n\n", status->getBestSolution()->getSolutionValue());

  } while(!terminationCondition->terminate(*((SAStatus*)status)));

  return this->getBestSoFar();
} // end search

