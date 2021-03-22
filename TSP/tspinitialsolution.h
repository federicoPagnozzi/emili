#ifndef TSPINITIALSOLUTION_H
#define TSPINITIALSOLUTION_H


#include <algorithm>
#include <vector>

#include "../emilibase.h"

#include "tspinstance.h"
#include "tspsolution.h"
#include "tsp.h"

namespace emili {
namespace tsp {


class TSPInitialSolution: public emili::InitialSolution {

protected:
    TSP& problem;

public:
    TSPInitialSolution(TSP& problem_instance):
        emili::InitialSolution(problem_instance),
        problem(problem_instance) { }

    virtual emili::Solution* generateSolution(void)=0;

    virtual emili::Solution* generateEmptySolution(void) {
        std::vector< int > empty(problem.getInstance()->getn());
        TSPSolution* sol = new TSPSolution(empty);
        sol->setSolutionValue(std::numeric_limits<double>::max());
        return sol;
    }

    std::vector< int > generateEmptyIntVector(void) {
        std::vector< int > empty(problem.getInstance()->getn());
        return empty;
    }
    
}; // TSPInitialSolution


class TSPRandomInitialSolution: public TSPInitialSolution {

public:
    TSPRandomInitialSolution(tsp::TSP& problem_instance):
        TSPInitialSolution(problem_instance) { }

    virtual emili::Solution* generateSolution(void) {
        std::vector< int > rnd(problem.getInstance()->getn());
        std::iota(rnd.begin(), rnd.end(), 0);
        std::shuffle(rnd.begin(), rnd.end(), emili::getRandomGenerator());
        TSPSolution* sol = new TSPSolution(rnd);
        double value = problem.evaluateSolution(*sol);
        return sol;
    }

    std::vector< int > generateRandomPermutation(void) {
        std::vector< int > rnd(problem.getInstance()->getn());
        std::iota(rnd.begin(), rnd.end(), 0);
        std::shuffle(rnd.begin(), rnd.end(), emili::getRandomGenerator());
        return rnd;
    }

}; // TSPRandomInitialSolution


}
}
#endif
