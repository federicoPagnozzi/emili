#ifndef QAPPROBLEM_H
#define QAPPROBLEM_H


#include <iostream>
#include <cstdlib>

#include "qapinstance.h"
#include "qapsolution.h"

#include "../emilibase.h"

namespace emili{

namespace qap {

class QAP: public emili::Problem {

protected:
    QAPInstance instance;

public:
    QAP(QAPInstance& problemInstance):
        instance(problemInstance) { }
    QAP(char* instance_path):
        instance(QAPInstance(instance_path)) { }

    QAPInstance* getInstance(void) {
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

}; // QAP

}
} // namespace qap

#endif
