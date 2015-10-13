#ifndef SA_COOLING_H
#define SA_COOLING_H


#include "sa_init_temp.h"
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
    double reset_threshold;
    double inittemp;

public:

    SACooling(double a, double b, SAInitTemp *it):
        a(a),
        b(b),
        maxIterations(0),
        counter(0),
        step(1),
        inittemp(it->get()) { }

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

    void setResetThreshold(float value) {
        reset_threshold = value;
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

    GeomCooling(double a, double b, SAInitTemp *it):
        SACooling(a, b, it) { }

    virtual double update_cooling(double temp) {
        counter++;

        if (counter >= maxIterations) {
            counter = 0;
            step++;
            float tmp = a * std::pow(b, temp);

            if (tmp <= reset_threshold)
                return inittemp;

            return tmp;
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

    MarkovCooling(double a, double b, SAInitTemp *it):
        SACooling(a, b, it) { }

    virtual double update_cooling(double temp) {
        counter++;

        if (counter >= maxIterations) {
            counter = 0;
            step++;
            float tmp = (a / (b + std::log(temp)));

            if (tmp <= reset_threshold)
                return inittemp;

            return tmp;

        }

        return(temp);
    }
}; // MarkovCooling


/**
 * Log-based cooling.
 */
class LogCooling: public SACooling {

public:
    LogCooling(double a, double b, SAInitTemp *it):
        SACooling(a, b, it) { }

    virtual double update_cooling(double temp) {
        counter++;

        if (counter >= maxIterations) {
            counter = 0;
            step++;
            float tmp = (a / std::log(temp + b));

            if (tmp <= reset_threshold)
                return inittemp;

            return tmp;

        }

        return(temp);
    }

}; // LogCooling


/**
 * http://www.sciencedirect.com/science/article/pii/037722179090301Q#
 */
class ConstantCooling: public SACooling {

public:
    ConstantCooling(double a, double b, SAInitTemp *it):
        SACooling(a, b, it) { }

    virtual double update_cooling(double temp) {
        counter++;

        if (counter >= maxIterations) {
            counter = 0;
            step++;
            float tmp = (a / (1 + b*temp));

            if (tmp <= reset_threshold)
                return inittemp;

            return tmp;

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
    LundyMeesCooling(double a, double b, SAInitTemp *it):
        SACooling(a, b, it) { }

    virtual double update_cooling(double temp) {
        counter++;

        if (counter >= maxIterations) {
            counter = 0;
            step++;
            float tmp = (temp / (a + b*temp));

            if (tmp <= reset_threshold)
                return inittemp;

            return tmp;

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
    LinearCooling(double a, SAInitTemp *it):
        SACooling(a, 0, it) { }

    virtual double update_cooling(double temp) {
        counter++;

        if (counter >= maxIterations) {
            counter = 0;
            step++;
            float tmp = a*temp;

            if (tmp <= reset_threshold)
                return inittemp;

            return tmp;

        }

        return(temp);
    }

}; // LinearCooling


/**
 * No Cooling - contant temperature
 */
class NoCooling: public SACooling {

public:
    NoCooling(SAInitTemp *it):
        SACooling(0, 0, it) { }

    virtual double update_cooling(double temp) {
        return temp;
    }

}; // NoCooling

#endif
