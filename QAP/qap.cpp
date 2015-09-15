#include "qap.h"


double qap::QAP::evaluateSolution(emili::Solution& solution) {
    QAPSolution& s = dynamic_cast<QAPSolution&> (solution);
    double p = this->getInstance()->computeObjectiveFunction(&s);
    solution.setSolutionValue(p);
    return p;
}

int qap::QAP::computeObjectiveFunction(std::vector< int > & partial_solution) {
    int i, j, n = this->getInstance()->getn();
    float value = 0.0;

    std::vector<std::vector< int > > A = this->getInstance()->getA(); 
    std::vector<std::vector< int > > B = this->getInstance()->getB();

    for (i = 0 ; i < n ; i++) {
        for (j = 0 ; j < n ; j++) {
            value += A[i][j] * B[partial_solution[i]][partial_solution[j]];
        }
    }

    if (this->getInstance()->is_made_symmetric())
        return value / 2;

    return value;
}

int qap::QAP::computeObjectiveFunction(std::vector< int > & partial_solution, int size) {
    return this->computeObjectiveFunction(partial_solution);
}
