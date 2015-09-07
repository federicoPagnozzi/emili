#ifndef QAPINSTANCE_H
#define QAPINSTANCE_H


#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "qapsolution.h"


using namespace std;


typedef int matrixEl;

/**
 * QAPInstance
 */
class QAPInstance {

private:
    int n;
    vector< vector< matrixEl > > A;
    vector< vector< matrixEl > > B;

    vector< int > optPerm;

    QAPSolution *solution;

    float optValue;

    /* bestKnownValue is the value read in Thomas' QAP instances */
    float bestKnownValue;

public:
    QAPInstance(void);
    /**
     * initialize an instance from a QAPLib instance file
     *
     * Please provide an existing && valid file
     */
    QAPInstance(string QAPLibFile);
    ~QAPInstance(void);

    /* Read/Write private attributes */
    int getn(void) {
        return n;
    }

    void setn(int _n) {
        n = _n;
    }

    vector< vector< matrixEl > > getA(void) {
        return A;
    }
    vector< vector< matrixEl > > getB(void) {
        return B;
    }

    void setA(vector< vector< matrixEl > > _A) {
        A = _A;
    }
    void setB(vector< vector< matrixEl > > _B) {
        B = _B;
    }

    QAPSolution* getSolution(void) {
        return solution;
    }

    void setSolution(QAPSolution* _solution) {
        solution = _solution;
    }

    string toString(void);

    float computeObjectiveFunction(QAPSolution* _solution);

    QAPSolution* getRandomSolution(void) {

    }

}; // QAPInstance

#endif
