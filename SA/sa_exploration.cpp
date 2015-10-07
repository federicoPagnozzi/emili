#include "sa_exploration.h"

emili::Solution* SARandomExploration::nextSolution(emili::Solution *startingSolution,
                                                           int* counter) {


    counter[0]++;

    emili::Solution* incumbent = neigh->random(startingSolution);
    //std::cout << startingSolution->getSolutionRepresentation();
    //std::cout << incumbent->getSolutionRepresentation();
    emili::Solution* accepted = acceptance->accept(startingSolution,
                                                   incumbent);
    std::string tc_type = term->getType();

    //std::cout << startingSolution << " " << incumbent << " " << accepted << std::endl;
    //std::cout << startingSolution->getSolutionValue() << " " << incumbent->getSolutionValue() << " " << accepted->getSolutionValue() << std::endl;

    if (accepted == startingSolution) {
        delete incumbent;
    } else {
        delete startingSolution;
        if (tc_type == MAXBADITERS) {
            counter[0] = 0;
        }
    }
    startingSolution = accepted;

    return startingSolution;

}


emili::Solution* SASequentialExploration::nextSolution(emili::Solution *startingSolution,
                                                           int* counter) {


    std::string tc_type = term->getType();

    emili::Solution* incumbent;
    emili::Solution* accepted;

    for(emili::Neighborhood::NeighborhoodIterator iter = neigh->begin(startingSolution);
        iter!=neigh->end();
        ++iter) {

        incumbent = *iter;

        counter[0]++;
    
        //std::cout << startingSolution->getSolutionRepresentation();
        //std::cout << incumbent->getSolutionRepresentation();
    
        accepted = acceptance->accept(startingSolution,
                                      incumbent);

        //std::cout << startingSolution << " " << incumbent << " " << accepted << std::endl;
        //std::cout << startingSolution->getSolutionValue() << " " << incumbent->getSolutionValue() << " " << accepted->getSolutionValue() << std::endl;

        if (accepted == startingSolution) {
            delete incumbent;
        } else {
            delete startingSolution;
            if (tc_type == MAXBADITERS) {
                counter[0] = 0;
            }
            break;
        }

    } 

    startingSolution = accepted;

    return startingSolution;

}
