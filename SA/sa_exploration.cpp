#include "sa_exploration.h"

emili::Solution* SARandomExploration::nextSolution(emili::Solution *startingSolution,
                                                   SAStatus& status) {

    status.increment_counters();

    emili::Solution* incumbent = neigh->random(startingSolution);
    emili::Solution* accepted = acceptance->accept(startingSolution,
                                                   incumbent);

    if (accepted == startingSolution) {
        delete incumbent;
        status.not_accepted_sol();
    } else {
        delete startingSolution;
        status.accepted_sol(accepted->getSolutionValue());
    }
    startingSolution = accepted;

    return startingSolution;

}

emili::Solution* SARandomExplorationNoCopy::nextSolution(emili::Solution *startingSolution, SAStatus &status) {
   status.increment_counters();

   double costBefore = startingSolution->getSolutionValue();
   neigh->randomStep(startingSolution);
   double delta = startingSolution->getSolutionValue() - costBefore;

   if(! acceptance->acceptViaDelta(startingSolution, delta)) {
       neigh->reverseLastRandomStep(startingSolution);
       status.not_accepted_sol();
   } else {
       status.accepted_sol(startingSolution->getSolutionValue());
   }

   return startingSolution;
}


emili::Solution* SASequentialExploration::nextSolution(emili::Solution *startingSolution,
                                                       SAStatus& status) {


    emili::Solution* incumbent = startingSolution->clone();
    emili::Solution* accepted;
    emili::Solution* ithSolution;

    emili::Neighborhood::NeighborhoodIterator iter = neigh->begin(incumbent);
    ithSolution = *iter; // isn't there a mistake ? ithSolution does not vary

    for(;
        iter!=neigh->end();
        ++iter) {

        status.increment_counters();
 
        accepted = acceptance->accept(incumbent,
                                      ithSolution);

        if (accepted == incumbent) {
            status.not_accepted_sol();
        } else {
            delete startingSolution;
            *accepted = *ithSolution;
            status.accepted_sol(accepted->getSolutionValue());
            break;
        }

    } 

    delete incumbent;
    
    return accepted;
}
