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

public:

    SACooling(double a, double b):
        a(a),
        b(b) { }

    /**
     * SA cooling scheme
     * @param  temp initial temperature
     * @return      updated temperature
     */
    virtual double update_cooling(double temp)=0;

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
        return (a * std::pow(b, temp));
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
        return (a / (b + std::log(temp)));
    }
}; // MarkovCooling

#endif
