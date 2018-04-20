#ifndef QAPINSTANCE_H
#define QAPINSTANCE_H


#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "qapsolution.h"

namespace emili {
namespace qap {

using namespace std;

#define XOR(x,y) ((x && !y) || (!x && y))


typedef long matrixEl;

/**
 * QAPInstance
 */
class QAPInstance {

private:
    int n;
    /* matrices: they might not be the original ones
       after symmetry manipulation */
    // distance matrix
    vector< vector< matrixEl > > A;
    // flow matrix
    vector< vector< matrixEl > > B;

    // original distance matrix
    vector< vector< matrixEl > > A_orig;
    // original flow matrix
    vector< vector< matrixEl > > B_orig;

    vector< int > optPerm;

    QAPSolution *solution;

    float optValue;

    /* bestKnownValue is the value read in Thomas' QAP instances */
    float bestKnownValue;

    /* symmetry flags */
    bool null_diagonal_flag  = false;  /* at least one matrix has zero diagonal: TRUE */
    bool d_symmetric_flag    = false;  /* if first (d) matrix is symmetric: TRUE */
    bool f_symmetric_flag    = false;  /* if second (f) matrix is symmetric: TRUE */
    bool make_symmetric_flag = false;  /* convert asymmetric instance into symmetric 
                                          instance: TRUE
                                        */
    /* symmetry checking and handling */
    bool check_null_diagonal(vector< vector< matrixEl > > matrix);
    bool check_symmetry(vector< vector< matrixEl > > matrix);
    vector< vector< matrixEl > > make_matrix_symmetric(vector< vector< matrixEl > > matrix);

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

    double computeObjectiveFunction(QAPSolution* _solution);

    bool is_null_diagonal(void) {
        return null_diagonal_flag;
    }

    bool is_d_symmetric(void) {
        return d_symmetric_flag;
    }

    bool is_f_symmetric(void) {
        return f_symmetric_flag;
    }

    bool is_made_symmetric(void) {
        return make_symmetric_flag;
    }

    void read_line(ifstream& s) {
        char tmp[10000];
        s.getline(tmp, sizeof(tmp));
    }

}; // QAPInstance

}
}
#endif
