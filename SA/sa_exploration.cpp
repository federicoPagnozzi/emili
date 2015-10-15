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
        if (tc_type == LASTACCRATETERM) {
            status.last_accepted[status.index] = 0;
            status.index = (status.index + 1) % status.tenure;
        }
    } else {
        delete startingSolution;
        status.accepted += 1;
        status.counter = 0;
        
        if (tc_type == LASTACCRATETERM) {
            status.last_accepted[status.index] = 1;
            status.index = (status.index + 1) % status.tenure;
        }
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
            // delete ithSolution;
            if (tc_type == LASTACCRATETERM) {
                status.last_accepted[status.index] = 0;
                status.index = (status.index + 1) % status.tenure;
            }
        } else {
            delete startingSolution;
            *accepted = *ithSolution;
            status.accepted += 1;
            status.counter = 0;

            if (tc_type == LASTACCRATETERM) {
                status.last_accepted[status.index] = 1;
                status.index = (status.index + 1) % status.tenure;
            }
            break;
        }

    } 

    delete incumbent;
    
    return accepted;
}
