#include "sa_exploration.h"

emili::Solution* SARandomExploration::nextSolution(emili::Solution *startingSolution,
                                                   sa_status* status) {


    status->counter += 1;
    status->total_counter += 1;

    emili::Solution* incumbent = neigh->random(startingSolution);
    //std::cout << startingSolution->getSolutionRepresentation();
    //std::cout << incumbent->getSolutionRepresentation();
    emili::Solution* accepted = acceptance->accept(startingSolution,
                                                   incumbent);

    //std::cout << startingSolution << " " << incumbent << " " << accepted << std::endl;
    //std::cout << startingSolution->getSolutionValue() << " " << incumbent->getSolutionValue() << " " << accepted->getSolutionValue() << std::endl;

    if (accepted == startingSolution) {
        delete incumbent;
        if (tc_type == LASTACCRATETERM) {
            status->last_accepted[status->index] = 0;
            status->index = (status->index + 1) % status->tenure;
        }
    } else {
        delete startingSolution;
        status->accepted += 1;
        if (tc_type == MAXBADITERS) {
            status->counter = 0;
        } else if (tc_type == LASTACCRATETERM) {
            status->last_accepted[status->index] = 1;
            status->index = (status->index + 1) % status->tenure;
        }
    }
    startingSolution = accepted;

    return startingSolution;

}


emili::Solution* SASequentialExploration::nextSolution(emili::Solution *startingSolution,
                                                       sa_status* status) {


    emili::Solution* incumbent;
    emili::Solution* accepted;

    for(emili::Neighborhood::NeighborhoodIterator iter = neigh->begin(startingSolution);
        iter!=neigh->end();
        ++iter) {

        incumbent = *iter;

        status->counter += 1;
        status->total_counter += 1;
    
        //std::cout << startingSolution->getSolutionRepresentation();
        //std::cout << incumbent->getSolutionRepresentation();
    
        accepted = acceptance->accept(startingSolution,
                                      incumbent);

        //std::cout << startingSolution << " " << incumbent << " " << accepted << std::endl;
        //std::cout << startingSolution->getSolutionValue() << " " << incumbent->getSolutionValue() << " " << accepted->getSolutionValue() << std::endl;

        if (accepted == startingSolution) {
            delete incumbent;
            if (tc_type == LASTACCRATETERM) {
                status->last_accepted[status->index] = 0;
                status->index = (status->index + 1) % status->tenure;
            }
        } else {
            delete startingSolution;
            status->accepted += 1;
            if (tc_type == MAXBADITERS) {
                status->counter = 0;
            } else if (tc_type == LASTACCRATETERM) {
                status->last_accepted[status->index] = 1;
                status->index = (status->index + 1) % status->tenure;
            }
            break;
        }

    } 

    startingSolution = accepted;

    return startingSolution;

}
