// #include "metropolis.h"
// using namespace emili::sa;
// emili::Solution* Metropolis::search() {
//     emili::Solution* current = init->generateSolution();
//     emili::Solution* sol = Metropolis::search(current);
//     if(current!=sol)
//     delete current;

//     return sol;
// } // end search



// emili::Solution* Metropolis::search(emili::Solution* initial) {
//     emili::Solution* incumbent = init->generateEmptySolution();
//     emili::Solution* tmpSol;
//     *incumbent = *initial;
//     bestSoFar  = incumbent;

//     do {
//         tmpSol = ftsa->search(bestSoFar);
//         bestSoFar = dsa->search(tmpSol);
//     } while(!terminationCriterion->terminate(*status));

//     delete bestSoFar;
//     delete tmpSol;

//     return status->best;
// } // end search

// void Metropolis::reset(void) {
//     status->temp = status->init_temp;
// }


// void Metropolis::setSearchTime(int _time) {
//     if (status->tc_type == NEVERTERM)
//         coolingScheme->set_search_time(_time);

//     this->seconds = _time;

// }
