#ifndef SA_INIT_TEMP_H
#define SA_INIT_TEMP_H

#include <string>
#include <iostream>
#include <vector>
#include <functional>
#include <cmath>

#include "../emilibase.h"


/**
 * Generic initial temperature setup for SA.
 *
 * The actual scheme has to be implemented with derived classes.
 */
class SAInitTemp {

protected:
    emili::Solution  *solution;
    double            init_temp;

public:

    /**
     * Empty constructor
     */
    SAInitTemp(void):
        solution(nullptr) { }
    /**
     * Constructor: set initial temperature starting from a valid solution.
     */
    SAInitTemp(emili::Solution *solution):
        solution(solution) { }


    /**
     * set initial temperature.
     *
     * If the initial temperature is set by the user,
     * value is the starting value for the temperature.
     *
     * If the initial temperature is set starting from a solution,
     * then value is a coefficient that rescales the cost of the solution.
     * 
     * @param  value initial temperature or coefficient.
     * @return       initial temperature
     */
    virtual double set(double value)=0;

}; // class SAInitTemp


/************************************************
 *
 * Implementation of subclasses
 * 
 ************************************************/


/**
 * Set the initial temperature at a fixed value.
 *
 * It just uses the value provieded by the user.
 */
class FixedInitTemp: public SAInitTemp {
public:

    FixedInitTemp(void):
        SAInitTemp() { }

    virtual double set(double value) {
        init_temp = value;
        return init_temp;
    }

}; // FixedInitTemp


/**
 * Set the initial temperature starting from an initial solution.
 *
 * The initial temperature is computed according to the formula:
 * init_temp = coeff * cost(solution)
 */
class InitTempFromSolution: public SAInitTemp {
public:

    InitTempFromSolution(emili::Solution *solution):
        SAInitTemp(solution) { }

    virtual double set(double value) {
        init_temp = value * solution->getSolutionValue();;
        return init_temp;
    }

}; // InitTempFromSolution


/**
 * Do a random walk and take the (absolute value of the) highest gap
 * as initial temperature
 */
class RandomWalkInitTemp: public SAInitTemp {
protected:
    emili::InitialSolution *is;
    int length;

public:

    RandomWalkInitTemp(emili::InitialSolution* _is,
                       int _length):
        is(_is),
        length(_length) { }

    virtual double set(double value) {
        int i;

        emili::Solution *s1;
        emili::Solution *s2 = is->generateSolution();
        
        double c1,
               c2 = s2->getSolutionValue();
        double maxdelta,
               bestcost = c2;

        solution = s2;

        for (i = 0 ; i < length ; i++) {
            s1 = s2;
            s2 = is->generateSolution();
            c1 = c2;
            c2 = s2->getSolutionValue();

            if (abs(c2 - c1) > maxdelta) {
                maxdelta = abs(c2 - c1);
            }

            if (c2 < bestcost) {
                delete solution;
                solution = s2;
                bestcost = c2;
            } else {
                delete s1;
            }
        }

        if (solution != s2)
            delete s2;

        return value * maxdelta;
    }

}; // RandomWalkInitTemp


#endif

