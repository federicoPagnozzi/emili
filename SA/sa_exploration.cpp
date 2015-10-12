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


    neigh->reset();

    emili::Solution* incumbent = startingSolution->clone();
    emili::Solution* accepted;
    emili::Solution* ithSolution;

    std::cout << startingSolution->getSolutionRepresentation() << std::endl;
    std::cout << "PPP" << status->total_counter << std::endl;
    std::cout << "YYY" << std::endl;

    emili::Neighborhood::NeighborhoodIterator iter = neigh->begin(incumbent);
    ithSolution = *iter;
        
    std::cout << "MMMMMMM" << std::endl;

    std::cout << ithSolution << std::endl;
    //std::cout << ithSolution->getSolutionRepresentation() << std::endl;
    //    getchar();
    std::cout << "NNN" << std::endl;


    for(;
        iter!=neigh->end();
        ++iter) {

        std::cout << "YJYFHGFHJGJHG" << std::endl;

        status->counter += 1;
        status->total_counter += 1;

        std::cout << status->total_counter << " " << status->accepted << std::endl;
    
        std::cout << ithSolution->getSolutionRepresentation();
        std::cout << ithSolution->getSolutionValue() << std::endl;
    
        accepted = acceptance->accept(incumbent,
                                      ithSolution);

        std::cout << "AAA" << std::endl;
        std::cout << incumbent->getSolutionRepresentation();
        std::cout << ithSolution->getSolutionRepresentation();
        std::cout << accepted->getSolutionRepresentation();

        //std::cout << startingSolution << " " << incumbent << " " << accepted << std::endl;
        //std::cout << startingSolution->getSolutionValue() << " " << incumbent->getSolutionValue() << " " << accepted->getSolutionValue() << std::endl;

        if (accepted == incumbent) {
            // delete ithSolution;
            if (tc_type == LASTACCRATETERM) {
                status->last_accepted[status->index] = 0;
                status->index = (status->index + 1) % status->tenure;
            }
            std::cout << "A" << std::endl;
        } else {
            // delete startingSolution;
            status->accepted += 1;
            if (tc_type == MAXBADITERS) {
                status->counter = 0;
            } else if (tc_type == LASTACCRATETERM) {
                status->last_accepted[status->index] = 1;
                status->index = (status->index + 1) % status->tenure;
            }
            std::cout << "B" << std::endl;
            break;
        }

    } 

    getchar();

    // delete ithSolution;
    // delete incumbent;
    
    startingSolution = accepted->clone();

    std::cout << startingSolution->getSolutionRepresentation();

    return startingSolution;
    
    //return accepted;

}
