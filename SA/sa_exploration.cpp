#include "sa_exploration.h"

emili::Solution* SARandomExploration::nextSolution(emili::Solution *startingSolution,
                                                   SAStatus& status) {


    status.counter += 1;
    status.total_counter += 1;

    emili::Solution* incumbent = neigh->random(startingSolution);
    emili::Solution* accepted = acceptance->accept(startingSolution,
                                                   incumbent);

    if (accepted == startingSolution) {
        delete incumbent;
        status.not_accepted_sol();
    } else {
        delete startingSolution;
        status.accepted_sol();
    }
    startingSolution = accepted;

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

        status.counter += 1;
        status.total_counter += 1;
 
        accepted = acceptance->accept(incumbent,
                                      ithSolution);

        if (accepted == incumbent) {
            status.not_accepted_sol();
        } else {
            delete startingSolution;
            *accepted = *ithSolution;
            status.accepted_sol();
            break;
        }

    } 

    delete incumbent;
    
    return accepted;
}
