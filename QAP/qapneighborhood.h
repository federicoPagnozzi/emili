#ifndef QAPNEIGHBORHOOD_H
#define QAPNEIGHBORHOOD_H


#include <cassert>
#include <limits>

#include "../emilibase.h"

#include "qapsolution.h"
#include "qapinstance.h"
#include "qapinitialsolution.h"
#include "qap.h"
namespace emili {
namespace qap{

/**
 * basic neighborhood for QAP
 */
class QAPNeighborhood: public emili::Neighborhood {

protected:
    qap::QAP& problem_instance;
    virtual emili::Solution* computeStep(emili::Solution* step)=0;
    virtual void reverseLastMove(emili::Solution* step)=0;

public:
    QAPNeighborhood(qap::QAP& problem_instance):
        problem_instance(problem_instance) { }
    ~QAPNeighborhood(void);

    virtual emili::Solution* step(emili::Solution* currentSolution);
    virtual void reset()=0;
    virtual std::pair<int,int> lastMove() { return std::pair<int,int>(0,0); }
    virtual int size()=0;

    // i don't think a setter method is useful.
    qap::QAP& getProblemInstance(void) {
        return problem_instance;
    }

}; // QAPNeighborhood


/**
 * 
 * INSERTION NEIGHBORHOOD
 * 
 */


/**
 * insertion neighborhood for QAP
 */
class QAPInsertNeighborhood: public QAPNeighborhood {

protected:
    int start_position;
    int end_position;
    int sp_iterations;
    int ep_iterations;
    std::vector < int > current;
    int current_value;
    int n;
    bool first_iter;
    bool symmetric, make_symmetric;
    QAPInstance* instance;
    virtual emili::Solution* computeStep(emili::Solution* value);
    virtual void reverseLastMove(emili::Solution* step) { }

public:
    QAPInsertNeighborhood(qap::QAP& problem_instance):
        start_position(0),
        end_position(0),
        sp_iterations(1),
        ep_iterations(1),
        QAPNeighborhood(problem_instance) {
            QAPInstance* instance = problem_instance.getInstance();
            n = instance->getn();
            //d = instance->getA();
            //f = instance->getB();

            //vector< matrixEl > vtemp(n, 0);
            //move_values = vector< vector< matrixEl > >(n, vtemp);

            /*if ( instance->is_made_symmetric() ||
                (instance->is_d_symmetric() && instance->is_f_symmetric() && instance->is_null_diagonal())
               ) {
                symmetric = true;
            } else {
                symmetric = false;
            }

            make_symmetric = instance->is_made_symmetric();*/
        }
    virtual ~QAPInsertNeighborhood(void);

    virtual void reset();
    virtual emili::Solution* random(emili::Solution *currentSolution);//=0
    virtual std::pair<int,int> lastMove() { return std::pair<int,int>(end_position,start_position); }
    virtual NeighborhoodIterator begin(emili::Solution *base);//=0

    virtual int size();

    virtual emili::Solution* step(emili::Solution* currentSolution);

}; // QAPInsertNeighborhood


/**
 * 
 * EXCHANGE NEIGHBORHOOD
 * 
 */


class QAPExchangeNeighborhood: public QAPNeighborhood {

protected:
    int start_position;
    int end_position;
    int sp_iterations;
    int ep_iterations;
    std::vector < int > current;
    double current_value;

    QAPInstance* instance;

    // additional info for speedup
    vector< vector< matrixEl > > move_values;
    vector< vector< matrixEl > > d;
    vector< vector< matrixEl > > f;
    //vector< matrixEl >* x;
    int n;
    bool first_iter;
    bool symmetric, make_symmetric;

    virtual emili::Solution* computeStep(emili::Solution* value);
    virtual void reverseLastMove(emili::Solution* step);

    double computeDelta(int u, int v, vector< int >& x);

public:
    QAPExchangeNeighborhood(qap::QAP& problem_instance):
        start_position(0),
        end_position(0),
        sp_iterations(1),
        ep_iterations(1),
        first_iter(true),
        QAPNeighborhood(problem_instance) {
            QAPInstance* instance = problem_instance.getInstance();
            n = instance->getn();
            d = instance->getA();
            f = instance->getB();
            
            vector< matrixEl > vtemp(n, 0);
            move_values = vector< vector< matrixEl > >(n, vtemp);

            if ( instance->is_made_symmetric() ||
                (instance->is_d_symmetric() && instance->is_f_symmetric() && instance->is_null_diagonal())
               ) {
                symmetric = true;
            } else {
                symmetric = false;
            }

            make_symmetric = instance->is_made_symmetric();
        }

    virtual ~QAPExchangeNeighborhood(void) { }

    virtual void reset();
    /**
     * just return a permutation
     */
    virtual emili::Solution* random(emili::Solution *currentSolution);
    virtual std::pair<int,int> lastMove() { return std::pair<int,int>(end_position,start_position); }
    virtual NeighborhoodIterator begin(emili::Solution *base);

    virtual int size();

    //virtual emili::Solution* step(emili::Solution* currentSolution)=0;

}; // QAPExchangeNeighborhood


class QAPFirst2optNeighborhood: public QAPExchangeNeighborhood {

protected:
    int start_position;
    int end_position;
    int sp_iterations;
    int ep_iterations;
    std::vector < int > current;
    int current_value;
    emili::Solution* first2opt_symmetric(emili::Solution* _value, bool make_symmetric_flag);
    emili::Solution* first2opt_asymmetric(emili::Solution* _value, bool make_symmetric_flag);
    virtual emili::Solution* computeStep(emili::Solution* value);

public:
    QAPFirst2optNeighborhood(qap::QAP& problem_instance):
        start_position(0),
        end_position(0),
        sp_iterations(1),
        ep_iterations(1),
        QAPExchangeNeighborhood(problem_instance) { }
    ~QAPFirst2optNeighborhood(void);

    virtual void reset();
    virtual std::pair<int,int> lastMove() { return std::pair<int,int>(end_position,start_position); }
    virtual NeighborhoodIterator begin(emili::Solution *base);

    /**
     * just return a permutation
     */
    virtual emili::Solution* random(emili::Solution *currentSolution);

    virtual int size();

    virtual emili::Solution* step(emili::Solution* currentSolution);

}; // QAPFirst2optNeighborhood


class QAPBest2optNeighborhood: public QAPExchangeNeighborhood {

protected:
    int start_position;
    int end_position;
    int sp_iterations;
    int ep_iterations;
    std::vector < int > current;
    int current_value;
    emili::Solution* best2opt_symmetric(emili::Solution* value, bool make_symmetric_flag);
    emili::Solution* best2opt_asymmetric(emili::Solution* value, bool make_symmetric_flag);
    virtual emili::Solution* computeStep(emili::Solution* value);

public:
    QAPBest2optNeighborhood(qap::QAP& problem_instance):
        start_position(0),
        end_position(0),
        sp_iterations(1),
        ep_iterations(1),
        QAPExchangeNeighborhood(problem_instance) { }
    ~QAPBest2optNeighborhood(void);

    virtual void reset();
    virtual std::pair<int,int> lastMove() { return std::pair<int,int>(end_position,start_position); }
    virtual NeighborhoodIterator begin(emili::Solution *base);

    /**
     * just return a permutation
     */
    virtual emili::Solution* random(emili::Solution *currentSolution);

    virtual int size();

    virtual emili::Solution* step(emili::Solution* currentSolution);

}; // QAPBest2optNeighborhood


}
}

#endif
