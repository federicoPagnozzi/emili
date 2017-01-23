#include "sa_exploration.h"

emili::Solution* SARandomExploration::nextSolution(emili::Solution *startingSolution,
                                                   SAStatus& status) {

    status.increment_counters();

    emili::Solution* incumbent = neigh->random(startingSolution);
    emili::Solution* accepted = acceptance->accept(startingSolution,
                                                   incumbent);
    status.temp = cooling->update_cooling(status.temp);
    acceptance->setCurrentTemp(status.temp);

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


emili::Solution* SASequentialExploration::nextSolution(emili::Solution *startingSolution,
                                                       SAStatus& status) {


    emili::Solution* incumbent = startingSolution->clone();
    emili::Solution* accepted;
    emili::Solution* ithSolution;

    emili::Neighborhood::NeighborhoodIterator iter = neigh->begin(incumbent);
    ithSolution = *iter;

    bool noneaccepted = true;

    for(;
        iter!=neigh->end();
        ++iter) {

        status.increment_counters();
        accepted = acceptance->accept(incumbent,
                                      ithSolution);
        status.temp = cooling->update_cooling(status.temp);
        acceptance->setCurrentTemp(status.temp);


        if (accepted == incumbent) {
            status.not_accepted_sol();
        } else {
            delete startingSolution;
            *accepted = *ithSolution;
            status.accepted_sol(accepted->getSolutionValue());
            noneaccepted = false;
            break;
        }
        

    } 

    delete incumbent;
    
    if (noneaccepted) {
        return startingSolution;
    }
    
    return accepted;
}
