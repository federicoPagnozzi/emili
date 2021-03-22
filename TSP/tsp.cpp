#include "tsp.h"

using namespace emili::tsp;

double TSP::evaluateSolution(emili::Solution& solution) {
    TSPSolution& s = dynamic_cast<TSPSolution&> (solution);
    double p = this->getInstance()->computeObjectiveFunction(&s);
    solution.setSolutionValue(p);
    return p;
}

double TSP::computeObjectiveFunction(std::vector< int > & partial_solution) {
    int i, j, n = this->getInstance()->getn();
    double value = 0.0;

    std::vector<std::vector< long > > A = this->getInstance()->getDistanceMatrix(); 

    std::cout << "\n\n\nOBJ  : ";
    for (i = 0 ; i < n-1 ; i++) {
            value += A[partial_solution[i]][partial_solution[i+1]];
            std::cout << A[partial_solution[i]][partial_solution[i+1]] << " ";
    }
    value += A[partial_solution[n-1]][partial_solution[0]];
    std::cout << A[partial_solution[n-1]][partial_solution[0]] << std::endl;

    return value;
}

double TSP::computeObjectiveFunction(std::vector< int > & partial_solution, int size) {
    return this->computeObjectiveFunction(partial_solution);
}

double TSP::calcObjectiveFunctionValue(emili::Solution& solution) {
    TSPSolution& s = dynamic_cast<TSPSolution&> (solution);
    return this->getInstance()->computeObjectiveFunction(&s);
}
