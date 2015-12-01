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
        status->new_best_solution(new_solution, ns, temperature);
    }

    return new_solution;

}


emili::Solution* SAMetropolisWithForcedAcceptance::accept(emili::Solution *current_solution,
                                                          emili::Solution *new_solution) {

    double cs = current_solution->getSolutionValue();
    double ns = new_solution->getSolutionValue();
    
    if (!status->force_accept && ns > cs) {
        double prob = std::exp((cs-ns) / temperature);

        if (prob < 1.0 && emili::generateRealRandomNumber() > prob) {
            return current_solution;
        }
    } else if (ns < status->best_cost) {
        status->new_best_solution(new_solution, ns, temperature);
    }

    return new_solution;

}


emili::Solution* SAApproxExpAcceptance::accept(emili::Solution *current_solution,
                                                emili::Solution *new_solution) {

    double cs = current_solution->getSolutionValue();
    double ns = new_solution->getSolutionValue();
    
    if (ns > cs) {
        double prob = 1 - ((cs-ns) / temperature);

        if (prob < 1.0 && emili::generateRealRandomNumber() > prob) {
            return current_solution;
        }
    } else if (ns < status->best_cost) {
        status->new_best_solution(new_solution, ns, temperature);
    }

    return new_solution;

}


emili::Solution* GeneralizedSAAcceptance::accept(emili::Solution *current_solution,
                                                 emili::Solution *new_solution) {

    double cs = current_solution->getSolutionValue();
    double ns = new_solution->getSolutionValue();
    
    if (ns > cs) {
        double prob = std::exp(std::pow(std::abs(cs), g) * (cs-ns) / temperature);

        if (prob < 1.0 && emili::generateRealRandomNumber() > prob) {
            return current_solution;
        }
    } else if (ns < status->best_cost) {
        status->new_best_solution(new_solution, ns, temperature);
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
        status->new_best_solution(new_solution, ns, temperature);
    }

    return new_solution;

}


emili::Solution* SAGeometricAcceptance::accept(emili::Solution *current_solution,
                                               emili::Solution *new_solution) {

    double cs = current_solution->getSolutionValue();
    double ns = new_solution->getSolutionValue();
    
    if (ns > cs) {
        double prob = start_temp * std::pow(rate, status->step - 1);

        if (prob < 1.0 && emili::generateRealRandomNumber() > prob) {
            return current_solution;
        }
    } else if (ns < status->best_cost) {
        status->new_best_solution(new_solution, ns, temperature);
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
        status->new_best_solution(new_solution, ns, temperature);
    }

    return new_solution;

}


emili::Solution* GreatDelugeAcceptance::accept(emili::Solution *current_solution,
                                               emili::Solution *new_solution) {

    double cs = current_solution->getSolutionValue();
    double ns = new_solution->getSolutionValue();
    
    if (ns > temperature) {
        return current_solution;
    } else if (ns < status->best_cost) {
        status->new_best_solution(new_solution, ns, temperature);
    }

    return new_solution;

}


emili::Solution* RecordToRecordAcceptance::accept(emili::Solution *current_solution,
                                                  emili::Solution *new_solution) {

    double cs = current_solution->getSolutionValue();
    double ns = new_solution->getSolutionValue();
    
    if (ns > status->best_cost * (1 + deviation/100.0)) {
        return current_solution;
    } else if (ns < status->best_cost) {
        status->new_best_solution(new_solution, ns, temperature);
    }

    return new_solution;

}


emili::Solution* LAHCAcceptance::accept(emili::Solution *current_solution,
                                        emili::Solution *new_solution) {

    double cs = current_solution->getSolutionValue();
    double ns = new_solution->getSolutionValue();

    int v = status->total_counter % tenure;

    if (ns > cs && ns > cost_list[v]) {
        return current_solution;
    } else if (ns < status->best_cost) {
        status->new_best_solution(new_solution, ns, temperature);
    }

    cost_list[v] = ns;

    return new_solution;

}
