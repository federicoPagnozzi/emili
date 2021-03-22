#include "sa_acceptance_criteria.h"
using namespace emili::sa;
emili::Solution* SAMetropolisAcceptance::accept(emili::Solution *current_solution,
                                                emili::Solution *new_solution) {

    double cs = current_solution->getSolutionValue();
    double ns = new_solution->getSolutionValue();
  
    //std::cout << std::fixed << cs << " " << ns << " " << status->total_counter << " " << temperature <<  " " << emili::getCurrentExecutionTime() << std::endl;
    
    //if (new_solution > current_solution) {
    if (ns > cs) {
        double prob = std::exp(-std::abs(cs-ns) / temperature); // (1.3806503e-23 * temperature) ?
        
        //std::cout << std::fixed << ns << " " << status->total_counter << " " << temperature <<  " " << emili::getCurrentExecutionTime() << " " << prob << std::endl;

        if (prob < 1.0 && emili::generateRealRandomNumber() > prob) {
            /** /acc_tracker[acc_pointer] = 0;
            acc_pointer = (acc_pointer+1) % acc_tsize;
            if (acc_pointer == 0) {
                double ratio = 0.0;
                for (int i = 0 ; i < acc_tsize ; i++) ratio = ratio + acc_tracker[i];
                std::cout << std::fixed << "ratio: " << 100 * ratio / acc_tsize << std::endl;
            }/ **/
  //          std::cout << "R";
            return current_solution;
        }


    } else if (ns < status->best_cost) {
    //} else if (new_solution < status->best) {
        status->new_best_solution(new_solution, ns, temperature);
    }
    
    //std::cout << std::fixed << ns << " " << status->total_counter << " " << temperature <<  " " << emili::getCurrentExecutionTime() << std::endl;

    /** /acc_tracker[acc_pointer] = 1;
    acc_pointer = (acc_pointer+1) % acc_tsize;
    if (acc_pointer == 0) {
       double ratio = 0.0;
       for (int i = 0 ; i < acc_tsize ; i++) ratio = ratio + acc_tracker[i];
       std::cout << std::fixed << "ratio: " << 100 * ratio / acc_tsize << std::endl;
    }/ **/
   // std::cout << "A";
    return new_solution;

}

emili::Solution* SAPrecomputedMetropolisAcceptance::accept(emili::Solution *current_solution,
                                                           emili::Solution *new_solution) {

    double cs = current_solution->getSolutionValue();
    double ns = new_solution->getSolutionValue();


    double ratio = -std::abs((cs - ns) / temperature);

    if (new_solution > current_solution) {
        if (ratio < -5.3) {
            return current_solution;
        } else if (ratio < 0) {
            double prob = probs[(int)(ratio / delta)];

            if (prob < 1.0 && emili::generateRealRandomNumber() > prob) {
                return current_solution;
            }
        }
    } else if (new_solution < status->best) {
        status->new_best_solution(new_solution, ns, temperature);
    }

    return new_solution;

}


emili::Solution* SAPrecomputedMetropolisWithForcedAcceptance::accept(emili::Solution *current_solution,
                                                                     emili::Solution *new_solution) {

    double cs = current_solution->getSolutionValue();
    double ns = new_solution->getSolutionValue();


    double ratio = -1.0 * std::abs((cs - ns) / temperature);

    if (new_solution > current_solution) {
        if (!status->force_accept && ratio < -5.3) {
            return current_solution;
        } else if (!status->force_accept && ratio < 0) {
            double prob = probs[(int)( -ratio / delta)];

            if (prob < 1.0 && emili::generateRealRandomNumber() > prob) {
                return current_solution;
            }

    // std::cout << std::fixed << ns << " " << status->total_counter << " " << temperature <<  " " << emili::getCurrentExecutionTime() << " " << prob << std::endl;
      }
    } else if (new_solution < status->best) {
        status->new_best_solution(new_solution, ns, temperature);
    }

    return new_solution;

}


emili::Solution* SAMetropolisWithForcedAcceptance::accept(emili::Solution *current_solution,
                                                          emili::Solution *new_solution) {

    double cs = current_solution->getSolutionValue();
    double ns = new_solution->getSolutionValue();
    
    // std::cout << cs << " " << ns << " " << " " << std::exp((cs-ns) / temperature) << std::endl;
    if (!status->force_accept && new_solution > current_solution) {
        double prob = std::exp(-std::abs(cs-ns) / temperature);

        if (prob < 1.0 && emili::generateRealRandomNumber() > prob) {
            /** /acc_tracker[acc_pointer] = 0;
            acc_pointer = (acc_pointer+1) % acc_tsize;
            if (acc_pointer == 0) {
                double ratio = 0.0;
                for (int i = 0 ; i < acc_tsize ; i++) ratio = ratio + acc_tracker[i];
                std::cout << std::fixed << "ratio: " << 100 * ratio / acc_tsize << std::endl;
            }/ **/
            return current_solution;
        }
    } else if (new_solution < status->best) {
        status->new_best_solution(new_solution, ns, temperature);
    }
    
    //std::cout << std::fixed << ns << " " << status->total_counter << " " << temperature <<  " " << emili::getCurrentExecutionTime() << std::endl;

    /** /acc_tracker[acc_pointer] = 1;
    acc_pointer = (acc_pointer+1) % acc_tsize;
    if (acc_pointer == 0) {
       double ratio = 0.0;
       for (int i = 0 ; i < acc_tsize ; i++) ratio = ratio + acc_tracker[i];
       std::cout << std::fixed << "ratio: " << 100 * ratio / acc_tsize << std::endl;
    }/ **/

    return new_solution;

}


emili::Solution* SAApproxExpAcceptance::accept(emili::Solution *current_solution,
                                                emili::Solution *new_solution) {

    double cs = current_solution->getSolutionValue();
    double ns = new_solution->getSolutionValue();
    
    if (new_solution > current_solution) {
        double prob = 1 - (-std::abs(cs-ns) / temperature);

        if (prob < 1.0 && emili::generateRealRandomNumber() > prob) {
            return current_solution;
        }
    } else if (new_solution < status->best) {
        status->new_best_solution(new_solution, ns, temperature);
    }

    return new_solution;

}


emili::Solution* GeneralizedSAAcceptance::accept(emili::Solution *current_solution,
                                                 emili::Solution *new_solution) {

    double cs = current_solution->getSolutionValue();
    double ns = new_solution->getSolutionValue();
    
    if (new_solution > current_solution) {
        double prob = std::exp(beta * std::pow(std::abs(cs), g) * (-1)*std::abs(cs-ns) / temperature);

        if (prob < 1.0 && emili::generateRealRandomNumber() > prob) {
            return current_solution;
        }
    } else if (new_solution < status->best) {
        status->new_best_solution(new_solution, ns, temperature);
    }

    return new_solution;

}


emili::Solution* SABasicAcceptance::accept(emili::Solution *current_solution,
                                           emili::Solution *new_solution) {

    double cs = current_solution->getSolutionValue();
    double ns = new_solution->getSolutionValue();
    
    if (new_solution > current_solution) {
        return current_solution;
    } else if (new_solution < status->best) {
        status->new_best_solution(new_solution, ns, temperature);
    }

    return new_solution;

}


emili::Solution* SAGeometricAcceptance::accept(emili::Solution *current_solution,
                                               emili::Solution *new_solution) {

    double cs = current_solution->getSolutionValue();
    double ns = new_solution->getSolutionValue();
    
    if (new_solution > current_solution) {
        double prob = status->init_prob * std::pow(rate, status->step - 1);

        if (prob < 1.0 && emili::generateRealRandomNumber() > prob) {
            return current_solution;
        }
    } else if (new_solution < status->best) {
        status->new_best_solution(new_solution, ns, temperature);
    }

    return new_solution;

}


emili::Solution* SADeterministicAcceptance::accept(emili::Solution *current_solution,
                                                   emili::Solution *new_solution) {

    double cs = current_solution->getSolutionValue();
    double ns = new_solution->getSolutionValue();
    
/*<<<<<<< HEAD
    if (ns > cs + temperature) {
  //      std::cout << "R";
        return current_solution;
    } else if (ns < status->best_cost) {
    //    std::cout << "B";
=======*/
    bool reject = true;
    if (new_solution > current_solution && ns > cs) {
      // minimization problem, worsening move
      reject = ns > cs + temperature;
    } else if (new_solution > current_solution && ns < cs) {
      // maximization problem, worsening move
      reject = ns < cs - temperature;
    } else {
      // whatever problem, it's an improving move
      reject = false;
    }

    if (reject) {
        return current_solution;
    } else if (new_solution < status->best) {
//>>>>>>> 18fad98daa8ee678663774f2e0cae7fb0b1df89f
        status->new_best_solution(new_solution, ns, temperature);
    }
//std::cout << "A";
    return new_solution;

}


emili::Solution* GreatDelugeAcceptance::accept(emili::Solution *current_solution,
                                               emili::Solution *new_solution) {

    double cs = current_solution->getSolutionValue();
    double ns = new_solution->getSolutionValue();
   
    bool reject = true;
    if (new_solution > current_solution && ns > cs) {
      // minimization problem, worsening move
      reject = ns > temperature;
    } else if (new_solution > current_solution && ns < cs) {
      // maximization problem, worsening move
      // pretty sure it will be impossible to make this work in practice...
      reject = 1.0/ns < 1.0/temperature;
    } else {
      // whatever problem, it's an improving move
      reject = false;
    }

 
    if (reject) {
        return current_solution;
    } else if (new_solution < status->best) {
        status->new_best_solution(new_solution, ns, temperature);
    }

    return new_solution;

}


emili::Solution* RecordToRecordAcceptance::accept(emili::Solution *current_solution,
                                                  emili::Solution *new_solution) {

    double cs = current_solution->getSolutionValue();
    double ns = new_solution->getSolutionValue();
    
    bool reject = true;
    if (new_solution > current_solution && ns > cs) {
      // minimization problem, worsening move
      reject = ns > status->best_cost * (1 + deviation/100.0);
    } else if (new_solution > current_solution && ns < cs) {
      // maximization problem, worsening move
      reject = ns < status->best_cost * (1 - deviation/100.0);
    } else {
      // whatever problem, it's an improving move
      reject = false;
    }

    if (reject) {
        return current_solution;
    } else if (new_solution < status->best) {
        status->new_best_solution(new_solution, ns, temperature);
    }

    return new_solution;

}


emili::Solution* LAHCAcceptance::accept(emili::Solution *current_solution,
                                        emili::Solution *new_solution) {

    double cs = current_solution->getSolutionValue();
    double ns = new_solution->getSolutionValue();

    int v = status->total_counter % tenure;

/*<<<<<<< HEAD
    if (ns > cs && ns > cost_list[v]) {
//        std::cout << "R";
        return current_solution;
    } else if (ns < status->best_cost) {
//std::cout << "B";

=======*/
    bool reject = true;
    if (new_solution > current_solution && ns > cs) {
      // minimization problem, worsening move
      reject = ns > cs && ns > cost_list[v];
    } else if (new_solution > current_solution && ns < cs) {
      // maximization problem, worsening move
      reject = ns < cs && ns < cost_list[v];
    } else {
      // whatever problem, it's an improving move
      reject = false;
    }

    if (reject) {
        return current_solution;
    } else if (new_solution < status->best) {
//>>>>>>> 18fad98daa8ee678663774f2e0cae7fb0b1df89f
        status->new_best_solution(new_solution, ns, temperature);
    }

    cost_list[v] = ns;
//std::cout << "A";

    return new_solution;

}

emili::Solution* LAHCNSAcceptance::accept(emili::Solution *current_solution,
                                        emili::Solution *new_solution) {

    double cs = current_solution->getSolutionValue();
    double ns = new_solution->getSolutionValue();

    int v = status->total_counter % tenure;

    bool reject = true;
    if (new_solution > current_solution && ns > cs) { 
      // minimization problem, worsening move
      reject = ns > cs && ns > cost_list[v];
    } else if (new_solution > current_solution && ns < cs) {
      // maximization problem, worsening move
      reject = ns < cs && ns < cost_list[v];
    } else {
      // whatever problem, it's an improving move
      reject = false;
    }

    if (reject) {
        return current_solution;
    } else if (new_solution < status->best) {
        status->new_best_solution(new_solution, ns, temperature);
    }

    cost_list[v] = ns;

    return new_solution;

}

emili::Solution* LAHCPSAcceptance::accept(emili::Solution *current_solution,
                                        emili::Solution *new_solution) {

    double cs = current_solution->getSolutionValue();
    double ns = new_solution->getSolutionValue();

    int v = status->total_counter % tenure;

    bool reject = true;
    if (new_solution > current_solution && ns > cs) { 
      // minimization problem, worsening move
      reject = ns > cs && ns > cost_list[v];
    } else if (new_solution > current_solution && ns < cs) {
      // maximization problem, worsening move
      reject = ns < cs && ns < cost_list[v];
    } else {
      // whatever problem, it's an improving move
      reject = false;
    }

    if (reject) {
        return current_solution;
    } else if (new_solution < status->best) {
        status->new_best_solution(new_solution, ns, temperature);
    }

    cost_list[v] = ns;

    return new_solution;

}


emili::Solution* SABoundedMetropolisAcceptance::accept(emili::Solution *current_solution,
                                                       emili::Solution *new_solution) {

    double cs = current_solution->getSolutionValue();
    double ns = new_solution->getSolutionValue();

    if (new_solution > current_solution) {
        double prob = std::exp(-std::abs(cs-ns) / temperature); // (1.3806503e-23 * temperature) ?
        double rd = 0.0;
        if (ns > cs) {
            // minimization problem
            rd = (ns - cs) / cs - 1;
        } else {
            // maximization problem
            rd = (cs - ns) / ns - 1;
        }

        if (prob < 1.0 && emili::generateRealRandomNumber() > prob && rd < reldelta) {
            return current_solution;
        }


    } else if (new_solution < status->best) {
        status->new_best_solution(new_solution, ns, temperature);
    }

    return new_solution;

}


emili::Solution* SAAcceptanceAll::accept(emili::Solution *current_solution,
                                         emili::Solution *new_solution) {

    double ns = new_solution->getSolutionValue();
        //printf(" %f status->best_cost: %f\n", ns, status->best_cost);

    if (new_solution < status->best) {
      status->new_best_solution(new_solution,
                                new_solution->getSolutionValue(),
                                temperature);
      //status->print();
    }

    return new_solution;

}

