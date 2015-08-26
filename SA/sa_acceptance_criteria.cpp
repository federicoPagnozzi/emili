#include "sa_acceptance_criteria.h"

emili::Solution* SAMetropolisAcceptance::accept(emili::Solution *current_solution,
	                                     emili::Solution *new_solution) {

    if (counter == interval && temperature > end_temp) {
        temperature = (alpha * temperature) - rate;        
        counter=0;
    }

    counter++;
    double cs = current_solution->getSolutionValue();
    double ns = new_solution->getSolutionValue();
    
    if (ns > cs) {
        double prob = std::exp((cs-ns) / temperature);

        if (prob < 1.0 && emili::generateRealRandomNumber() > prob) {
            return current_solution;
        }
    }

    return new_solution;

}


emili::Solution* SABasicAcceptance::accept(emili::Solution *current_solution,
	                                emili::Solution *new_solution) {

    if (counter == interval && temperature > end_temp) {
        temperature = (alpha * temperature) - rate;        
        counter=0;
    }

    counter++;
    double cs = current_solution->getSolutionValue();
    double ns = new_solution->getSolutionValue();
    
    if (ns > cs) {
        return current_solution;
    }

    return new_solution;

}
