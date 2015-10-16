#include "sa_acceptance_criteria.h"

emili::Solution* SAMetropolisAcceptance::accept(emili::Solution *current_solution,
	                                            emili::Solution *new_solution) {

    double cs = current_solution->getSolutionValue();
    double ns = new_solution->getSolutionValue();
    
    if (ns > cs) {
        double prob = std::exp((cs-ns) / temperature);

        if (prob < 1.0 && emili::generateRealRandomNumber() > prob) {
            return current_solution;
        }
    } else if (ns < status->best_cost) {
        status->new_best_solution(new_solution, ns);
    }

    return new_solution;

}


emili::Solution* SABasicAcceptance::accept(emili::Solution *current_solution,
	                                       emili::Solution *new_solution) {

    double cs = current_solution->getSolutionValue();
    double ns = new_solution->getSolutionValue();
    
    if (ns > cs) {
        return current_solution;
    } else if (ns < status->best_cost) {
        status->new_best_solution(new_solution, ns);
    }

    return new_solution;

}


emili::Solution* SAGeometricAcceptance::accept(emili::Solution *current_solution,
                                               emili::Solution *new_solution) {

    double cs = current_solution->getSolutionValue();
    double ns = new_solution->getSolutionValue();
    
    if (ns > cs) {
        double prob = start_temp * std::pow(rate, step - 1);

        if (prob < 1.0 && emili::generateRealRandomNumber() > prob) {
            return current_solution;
        }
    } else if (ns < status->best_cost) {
        status->new_best_solution(new_solution, ns);
    }

    return new_solution;

}


emili::Solution* SADeterministicAcceptance::accept(emili::Solution *current_solution,
                                                   emili::Solution *new_solution) {

    double cs = current_solution->getSolutionValue();
    double ns = new_solution->getSolutionValue();
    
    if (ns > cs * (1 + delta)) {
        return current_solution;
    } else if (ns < status->best_cost) {
        status->new_best_solution(new_solution, ns);
    }

    return new_solution;

}
