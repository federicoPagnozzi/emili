#include "qap.h"


double qap::QAP::evaluateSolution(emili::Solution& solution) {
    QAPSolution& s = dynamic_cast<QAPSolution&> (solution);
    double p = this->getInstance().computeObjectiveFunction(&s);
    solution.setSolutionValue(p);
    return p;
}

int qap::QAP::computeObjectiveFunction(std::vector< int > & partial_solution) {
    return 0;
}

int qap::QAP::computeObjectiveFunction(std::vector< int > & partial_solution, int size) {
    return 0;
}
