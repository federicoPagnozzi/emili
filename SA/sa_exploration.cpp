#include "sa_exploration.h"

#include <cassert>

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
       // assert(costBefore == startingSolution->getSolutionValue());
       status.not_accepted_sol();
       // nnacc++;
   } else {
       status.accepted_sol(startingSolution->getSolutionValue());
       // nacc++;
   }

   return startingSolution;
}


emili::Solution* SASequentialExploration::nextSolution(emili::Solution *startingSolution,
                                                       SAStatus& status) {


    emili::Solution* incumbent = startingSolution->clone();
    emili::Solution* accepted;
    emili::Solution* ithSolution;

    emili::Neighborhood::NeighborhoodIterator iter = neigh->begin(incumbent);
    ithSolution = *iter;

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
