#ifndef QAPINITIALSOLUTION_H
#define QAPINITIALSOLUTION_H


#include <algorithm>
#include <vector>

#include "../emilibase.h"

#include "qapinstance.h"
#include "qapsolution.h"
#include "qap.h"


class QAPInitialSolution: public emili::InitialSolution {

protected:
    qap::QAP& problem;

public:
    QAPInitialSolution(qap::QAP& problem_instance):
        emili::InitialSolution(problem_instance),
        problem(problem_instance) { }

    virtual emili::Solution* generateSolution(void)=0;

    virtual emili::Solution* generateEmptySolution(void) {
        std::vector< int > empty(problem.getInstance()->getn());
        QAPSolution* sol = new QAPSolution(empty);
        sol->setSolutionValue(std::numeric_limits<double>::max());
        return sol;
    }
    
}; // QAPInitialSolution


class QAPRandomInitialSolution: public QAPInitialSolution {

public:
    QAPRandomInitialSolution(qap::QAP& problem_instance):
        QAPInitialSolution(problem_instance) { }

    virtual emili::Solution* generateSolution(void) {
        std::vector< int > rnd(problem.getInstance()->getn());
        std::iota(rnd.begin(), rnd.end(), 0);
        std::shuffle(rnd.begin(), rnd.end(), emili::getRandomGenerator());
        QAPSolution* sol = new QAPSolution(rnd);
        double value = problem.evaluateSolution(*sol);
        return sol;
    }

}; // QAPRandomInitialSolution


#endif
