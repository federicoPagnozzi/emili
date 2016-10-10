#include "sa_acceptance_criteria.h"

emili::Solution* SAMetropolisAcceptance::accept(emili::Solution *current_solution,
                                                emili::Solution *new_solution) {

    double cs = current_solution->getSolutionValue();
    double ns = new_solution->getSolutionValue();

    
    if (ns > cs) {
        double prob = std::exp((cs-ns) / temperature); // (1.3806503e-23 * temperature) ?

        if (prob < 1.0 && emili::generateRealRandomNumber() > prob) {
            return current_solution;
        }
    } else if (ns < status->best_cost) {
        status->new_best_solution(new_solution, ns, temperature);
    }

    return new_solution;

}

bool SAMetropolisAcceptance::acceptViaDelta(emili::Solution *new_solution, double delta)
{
    // double cs = new_solution->getSolutionValue() - delta;
    double ns = new_solution->getSolutionValue();
    // delta == ns - cs
    // cs = ns - delta

    if(delta > 0) {
        // not improving
        double prob = std::exp(- delta / temperature); // (1.3806503e-23 * temperature) ?

        if(prob < 1.0 && emili::generateRealRandomNumber() > prob) {
            return false;
        }
    } else if(ns < status->best_cost) {
        // ns is better, We keep a copy
        status->new_best_solution(new_solution, ns, temperature);
    }

    return true;
}

bool SAMetropolisAcceptanceDebug::acceptViaDelta(emili::Solution *new_solution, double delta)
{
    // double cs = new_solution->getSolutionValue() - delta;
    double ns = new_solution->getSolutionValue();
    // delta == ns - cs
    // cs = ns - delta

    if(delta > 0) {
        // not improving
        double prob = std::exp(- delta / temperature); // (1.3806503e-23 * temperature) ?
        nnotimproving++;
        sumprob += prob;

        if(prob < 1.0 && emili::generateRealRandomNumber() > prob) {
            return false;
        }
    } else if(ns < status->best_cost) {
        // ns is better, We keep a copy
        status->new_best_solution(new_solution, ns, temperature);
    }

    return true;
}


emili::Solution* SAMetropolisWithForcedAcceptance::accept(emili::Solution *current_solution,
                                                          emili::Solution *new_solution) {

    double cs = current_solution->getSolutionValue();
    double ns = new_solution->getSolutionValue();
    
    // std::cout << "cost=" << cs << " newcost=" << ns << " " << " e^(Delta/T)=" << std::exp((cs-ns) / temperature) << std::endl;
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

bool SAMetropolisWithForcedAcceptance::acceptViaDelta(emili::Solution *new_solution, double delta) {
    auto ns = new_solution->getSolutionValue();
    // delta = ns - cs
    if(!status->force_accept && delta > 0) {
        double prob = std::exp(- delta / temperature);
        if(prob < 1 && emili::generateRandomNumber() > prob)
            return false;
    } else if(ns < status->best_cost) {
        status->new_best_solution(new_solution, ns, temperature);
    }

    return true;
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

bool SAApproxExpAcceptance::acceptViaDelta(emili::Solution *new_solution, double delta) {
    double ns = new_solution->getSolutionValue();

    if(delta > 0) {
        double prob = 1 - (- delta / temperature);

        if (prob < 1.0 && emili::generateRealRandomNumber() > prob)
            return false;
    } else if (ns < status->best_cost) {
        status->new_best_solution(new_solution, ns, temperature);
    }

    return new_solution;
}


emili::Solution* GeneralizedSAAcceptance::accept(emili::Solution *current_solution,
                                                 emili::Solution *new_solution) {

    double ns = new_solution->getSolutionValue();
    double cs = current_solution->getSolutionValue();

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

bool GeneralizedSAAcceptance::acceptViaDelta(emili::Solution *new_solution, double delta) {
    double ns = new_solution->getSolutionValue();
    double cs = ns - delta;

    if(delta > 0) {
        double prob = std::exp(std::pow(std::abs(cs), g) * (-delta) / temperature);
        if (prob < 1.0 && emili::generateRealRandomNumber() > prob)
            return false;
    } else if(ns < status->best_cost) {
        status->new_best_solution(new_solution, ns, temperature);
    }

    return true;
}


emili::Solution* SABasicAcceptance::accept(emili::Solution *current_solution,
                                           emili::Solution *new_solution) {

    double ns = new_solution->getSolutionValue();
    double cs = current_solution->getSolutionValue();

    if (ns > cs) {
        return current_solution;
    } else if (ns < status->best_cost) {
        status->new_best_solution(new_solution, ns, temperature);
    }

    return new_solution;

}

bool SABasicAcceptance::acceptViaDelta(emili::Solution *new_solution, double delta) {
    double ns = new_solution->getSolutionValue();

    if(delta > 0)
        return false;
    else if(ns < status->best_cost)
        status->new_best_solution(new_solution, ns, temperature);
    return true;
}


emili::Solution* SAGeometricAcceptance::accept(emili::Solution *current_solution,
                                               emili::Solution *new_solution) {

    double ns = new_solution->getSolutionValue();
    double cs = current_solution->getSolutionValue();

    if (ns > cs) {
        double prob = status->init_prob * std::pow(rate, status->step - 1);

        if (prob < 1.0 && emili::generateRealRandomNumber() > prob) {
            return current_solution;
        }
    } else if (ns < status->best_cost) {
        status->new_best_solution(new_solution, ns, temperature);
    }

    return new_solution;

}

bool SAGeometricAcceptance::acceptViaDelta(emili::Solution *new_solution, double delta) {
    double ns = new_solution->getSolutionValue();

    if(delta > 0) {
        double prob = status->init_prob * std::pow(rate, status->step - 1);
        if (prob < 1.0 && emili::generateRealRandomNumber() > prob)
            return false;
    } else if (ns < status->best_cost) {
        status->new_best_solution(new_solution, ns, temperature);
    }

    return true;
}


emili::Solution* SADeterministicAcceptance::accept(emili::Solution *current_solution,
                                                   emili::Solution *new_solution) {

    double ns = new_solution->getSolutionValue();
    double cs = current_solution->getSolutionValue();

    if (ns > cs * (1 + deltaParam)) {
        return current_solution;
    } else if (ns < status->best_cost) {
        status->new_best_solution(new_solution, ns, temperature);
    }

    return new_solution;
}

bool SADeterministicAcceptance::acceptViaDelta(emili::Solution *new_solution, double delta) {
    double ns = new_solution->getSolutionValue();
    double cs = ns - delta;

    if (ns > cs * (1 + deltaParam))
        return false;
    else if (ns < status->best_cost)
        status->new_best_solution(new_solution, ns, temperature);

    return true;
}


emili::Solution* GreatDelugeAcceptance::accept(emili::Solution *current_solution,
                                               emili::Solution *new_solution) {

    double cs = current_solution->getSolutionValue();
    double ns = new_solution->getSolutionValue();
    
    if (ns > temperature)
        return current_solution;
    else if (ns < status->best_cost)
        status->new_best_solution(new_solution, ns, temperature);

    return new_solution;
}

bool GreatDelugeAcceptance::acceptViaDelta(emili::Solution *new_solution, double delta) {
    double ns = new_solution->getSolutionValue();
    double cs = ns - delta;

    if(ns > temperature)
        return false;
    else if (ns < status->best_cost)
        status->new_best_solution(new_solution, ns, temperature);

    return true;
}


emili::Solution* RecordToRecordAcceptance::accept(emili::Solution *current_solution,
                                                  emili::Solution *new_solution) {

    double cs = current_solution->getSolutionValue();
    double ns = new_solution->getSolutionValue();
    
    if (ns > status->best_cost * (1 + deviation/100.0))
        return current_solution;
    else if (ns < status->best_cost)
        status->new_best_solution(new_solution, ns, temperature);

    return new_solution;

}

bool RecordToRecordAcceptance::acceptViaDelta(emili::Solution *new_solution, double delta) {
    double ns = new_solution->getSolutionValue();
    double cs = ns - delta;

    if(ns > status->best_cost * (1 + deviation/100.0))
        return false;
    else if (ns < status->best_cost)
        status->new_best_solution(new_solution, ns, temperature);

    return true;
}


emili::Solution* LAHCAcceptance::accept(emili::Solution *current_solution,
                                        emili::Solution *new_solution) {

    double cs = current_solution->getSolutionValue();
    double ns = new_solution->getSolutionValue();

    int v = status->total_counter % tenure;

    if (ns > cs && ns > cost_list[v])
        return current_solution;
    else if (ns < status->best_cost)
        status->new_best_solution(new_solution, ns, temperature);

    cost_list[v] = ns;

    return new_solution;

}

bool LAHCAcceptance::acceptViaDelta(emili::Solution *new_solution, double delta) {
    double ns = new_solution->getSolutionValue();

    int v = status->total_counter % tenure;

    if (delta > 0 && ns > cost_list[v])
        return false;
    else if (ns < status->best_cost)
        status->new_best_solution(new_solution, ns, temperature);

    cost_list[v] = ns;

    return true;
}
