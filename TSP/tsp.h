#ifndef TSPPROBLEM_H
#define TSPPROBLEM_H


#include <iostream>
#include <cstdlib>

#include "tspinstance.h"
#include "tspsolution.h"

#include "../emilibase.h"

namespace emili{

namespace tsp {

class TSP: public emili::Problem {

protected:
    TSPInstance instance;

public:
    TSP(TSPInstance& problemInstance):
        instance(problemInstance) { }
    TSP(char* instance_path):
        instance(TSPInstance(instance_path)) { }

    TSPInstance* getInstance(void) {
        return (&instance);
    }

    virtual int problemSize() {
        return instance.getn();
    } 

    /**
     computes the objective function value of solution.
     */
    virtual double calcObjectiveFunctionValue(emili::Solution& solution);
    virtual double evaluateSolution(emili::Solution& solution);
    virtual double computeObjectiveFunction(std::vector< int > & partial_solution);
    virtual double computeObjectiveFunction(std::vector< int > & partial_solution, int size);

}; // TSP

}
} // namespace tsp

#endif
