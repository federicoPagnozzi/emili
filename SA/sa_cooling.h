#ifndef SA_COOLING_H
#define SA_COOLING_H


#include "../emilibase.h"


/**
 * Generic cooling scheme for SA.
 *
 * The actual scheme has to be implemented with derived classes.
 */
class SACooling {

protected:
    double a;
    double b;
    int maxIterations;
    int counter;
    int step;

public:

    SACooling(double a, double b):
        a(a),
        b(b),
        maxIterations(0),
        counter(0),
        step(1) { }

    /**
     * SA cooling scheme
     * @param  temp initial temperature
     * @return      updated temperature
     */
    virtual double update_cooling(double temp)=0;


    /**
     * set the maximum number of iterations to be performed
     * at the same temperature
     * @param maxIt maximum number of iterations
     */
    void setMaxIterations(int maxIt) {
        maxIterations = maxIt;
    }

    int getStep(void) {
        return step;
    }

}; // SACooling


/************************************************
 *
 * Implementation of subclasses
 * 
 ************************************************/


/**
 * Geometric cooling scheme:
 * 
 * t_{i+1} = a * b^{t_{i}}
 */
class GeomCooling: public SACooling {
public:

    GeomCooling(double a, double b):
        SACooling(a, b) { }

    virtual double update_cooling(double temp) {
        counter++;

        if (counter >= maxIterations) {
            counter = 0;
            step++;
            return (a * std::pow(b, temp));
        }

        return(temp);
    }
}; // GeomCooling


/**
 * Markov-based cooling, guaranteed to converge to the optimal solution
 * in infinite time.
 */
class MarkovCooling: public SACooling{
public:

    MarkovCooling(double a, double b):
        SACooling(a, b) { }

    virtual double update_cooling(double temp) {
        counter++;

        if (counter >= maxIterations) {
            counter = 0;
            step++;
            return (a / (b + std::log(temp)));
        }

        return(temp);
    }
}; // MarkovCooling


/**
 * Log-based cooling.
 */
class LogCooling: public SACooling {

public:
    LogCooling(double a, double b):
        SACooling(a, b) { }

    virtual double update_cooling(double temp) {
        counter++;

        if (counter >= maxIterations) {
            counter = 0;
            step++;
            return (a / std::log(temp + b));
        }

        return(temp);
    }

}; // LogCooling


/**
 * http://www.sciencedirect.com/science/article/pii/037722179090301Q#
 */
class ConstantCooling: public SACooling {

public:
    ConstantCooling(double a, double b):
        SACooling(a, b) { }

    virtual double update_cooling(double temp) {
        counter++;

        if (counter >= maxIterations) {
            counter = 0;
            step++;
            return (a / (1 + b*temp));
        }

        return(temp);
    }


}; // ConstantCooling


/**
 * LundyMeesCooling
 * a = 1
 * b << T_0
 *
 * http://link.springer.com/article/10.1007/BF01582166
 * 
 */
class LundyMeesCooling: public SACooling {

public:
    LundyMeesCooling(double a, double b):
        SACooling(a, b) { }

    virtual double update_cooling(double temp) {
        counter++;

        if (counter >= maxIterations) {
            counter = 0;
            step++;
            return (temp / (a + b*temp));
        }

        return(temp);
    }

}; // LundyMeesCooling


/**
 * Linear cooling
 * 0 < a < 1
 *
 * equal to b*temp, where b = 1 - a
 */
class LinearCooling: public SACooling {

public:
    LinearCooling(double a):
        SACooling(a, 0) { }

    virtual double update_cooling(double temp) {
        counter++;

        if (counter >= maxIterations) {
            counter = 0;
            step++;
            return a*temp;
        }

        return(temp);
    }

}; // LinearCooling


/**
 * No Cooling - contant temperature
 */
class NoCooling: public SACooling {

public:
    NoCooling(void):
        SACooling(0, 0) { }

    virtual double update_cooling(double temp) {
        return temp;
    }

}; // NoCooling

#endif
