#ifndef TSPNEIGHBORHOOD_H
#define TSPNEIGHBORHOOD_H


#include <cassert>
#include <limits>

#include "../emilibase.h"

#include "tspsolution.h"
#include "tspinstance.h"
#include "tspinitialsolution.h"
#include "tsp.h"
namespace emili {
namespace tsp{

/**
 * basic neighborhood for TSP
 */
class TSPNeighborhood: public emili::Neighborhood {

protected:
    tsp::TSP& problem_instance;
    virtual emili::Solution* computeStep(emili::Solution* step)=0;
    virtual void reverseLastMove(emili::Solution* step)=0;

public:
    TSPNeighborhood(tsp::TSP& problem_instance):
        problem_instance(problem_instance) { }
    ~TSPNeighborhood(void);

    virtual emili::Solution* step(emili::Solution* currentSolution);
    virtual void reset()=0;
    virtual std::pair<int,int> lastMove() { return std::pair<int,int>(0,0); }
    virtual int size()=0;

    // i don't think a setter method is useful.
    tsp::TSP& getProblemInstance(void) {
        return problem_instance;
    }

}; // TSPNeighborhood


/**
 * 
 * INSERTION NEIGHBORHOOD
 * 
 */


/**
 * insertion neighborhood for TSP
 */
class TSPInsertNeighborhood: public TSPNeighborhood {

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
    TSPInstance* instance;
    virtual emili::Solution* computeStep(emili::Solution* value);
    virtual void reverseLastMove(emili::Solution* step) { }

public:
    TSPInsertNeighborhood(tsp::TSP& problem_instance):
        start_position(0),
        end_position(0),
        sp_iterations(1),
        ep_iterations(1),
        TSPNeighborhood(problem_instance) {
            TSPInstance* instance = problem_instance.getInstance();
            n = instance->getn();
        }
    virtual ~TSPInsertNeighborhood(void);

    virtual void reset();
    virtual emili::Solution* random(emili::Solution *currentSolution);//=0
    virtual std::pair<int,int> lastMove() { return std::pair<int,int>(end_position,start_position); }
    virtual NeighborhoodIterator begin(emili::Solution *base);//=0

    virtual int size();

    virtual emili::Solution* step(emili::Solution* currentSolution);

}; // TSPInsertNeighborhood


/**
 * 
 * EXCHANGE NEIGHBORHOOD
 * 
 */


class TSPExchangeNeighborhood: public TSPNeighborhood {

protected:
    int start_position;
    int end_position;
    int sp_iterations;
    int ep_iterations;
    std::vector < int > current;
    double current_value;

    TSPInstance* instance;

    // additional info for speedup
    vector< vector< matrixEl > > move_values;
    vector< vector< matrixEl > > d;
    //vector< matrixEl >* x;
    int n;
    bool first_iter;
    bool symmetric, make_symmetric;

    virtual emili::Solution* computeStep(emili::Solution* value);
    virtual void reverseLastMove(emili::Solution* step);

    double computeDelta(int u, int v, vector< int >& x);

public:
    TSPExchangeNeighborhood(tsp::TSP& problem_instance):
        start_position(0),
        end_position(0),
        sp_iterations(1),
        ep_iterations(1),
        first_iter(true),
        TSPNeighborhood(problem_instance) {
            TSPInstance* instance = problem_instance.getInstance();
            n = instance->getn();
            d = instance->getDistanceMatrix();
            
            vector< matrixEl > vtemp(n, 0);
            move_values = vector< vector< matrixEl > >(n, vtemp);
        }

    virtual ~TSPExchangeNeighborhood(void) { }

    virtual void reset();
    /**
     * just return a permutation
     */
    virtual emili::Solution* random(emili::Solution *currentSolution);
    virtual std::pair<int,int> lastMove() { return std::pair<int,int>(end_position,start_position); }
    virtual NeighborhoodIterator begin(emili::Solution *base);

    virtual int size();

    //virtual emili::Solution* step(emili::Solution* currentSolution)=0;

}; // TSPExchangeNeighborhood


}
}

#endif
