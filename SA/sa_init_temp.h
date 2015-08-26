#ifndef SA_INIT_TEMP_H
#define SA_INIT_TEMP_H

#include <string>
#include <iostream>
#include <vector>
#include <functional>

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

#endif
