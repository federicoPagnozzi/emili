#ifndef SA_COOLING_H
#define SA_COOLING_H

#include "../emilibase.h"

#include "sa_common.h"
#include "sa_init_temp.h"
#include "sa_templength.h"
#include "sa_temperature_restart.h"


/**
 * Generic cooling scheme for SA.
 *
 * The actual scheme has to be implemented with derived classes.
 */
class SACooling {

protected:
    double a;
    double b;
    int counter;
    double inittemp;

    SAStatus* status;
    SATempRestart* tempRestart;
    SATempLength*  tempLength;

public:

    SACooling(double a, double b, SAInitTemp *it):
        a(a),
        b(b),
        counter(0),
        inittemp(it->get()) {
            if (a < b) {
                std::swap(a, b);
            }
        }

    /**
     * SA cooling scheme
     * @param  temp initial temperature
     * @return      updated temperature
     */
    virtual double update_cooling(double temp)=0;

    void set_status(SAStatus* _status) {
        status = _status;
    }

    float get_init_temp(void) {
        return inittemp;
    }

    void setTempRestart(SATempRestart* _tempRestart) {
        tempRestart = _tempRestart;
    }

    void setTempLength(SATempLength* _tempLength) {
        tempLength = _tempLength;
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

        if (tempLength->isCoolingTime(counter)) {
            counter = 0;
            status->step = status->step + 1;
            float tmp = a * std::pow(b, temp);

            return tempRestart->adjust(tmp);
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

        if (tempLength->isCoolingTime(counter)) {
            counter = 0;
            status->step = status->step + 1;
            float tmp = (a / (b + std::log(temp)));

            return tempRestart->adjust(tmp);

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

        if (tempLength->isCoolingTime(counter)) {
            counter = 0;
            status->step = status->step + 1;
            float tmp = (a / std::log(temp + b));

            return tempRestart->adjust(tmp);

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

        if (tempLength->isCoolingTime(counter)) {
            counter = 0;
            status->step = status->step + 1;
            float tmp = (a / (1 + b*temp));

            return tempRestart->adjust(tmp);

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

        if (tempLength->isCoolingTime(counter)) {
            counter = 0;
            status->step = status->step + 1;
            float tmp = (temp / (a + b*temp));

            return tempRestart->adjust(tmp);

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

        if (tempLength->isCoolingTime(counter)) {
            counter = 0;
            status->step = status->step + 1;
            float tmp = a*temp;

            return tempRestart->adjust(tmp);

        }

        return(temp);
    }

}; // LinearCooling


class SATemperatureBandCooling: public SACooling {

public:

    SATemperatureBandCooling(double a, SAInitTemp *it):
        SACooling(a * it->get(), it->get(), it) { }

    virtual double update_cooling(double temp) {
        counter++;

        if (tempLength->isCoolingTime(counter)) {
            counter = 0;
            status->step = status->step + 1;
            float tmp = b + emili::generateRealRandomNumber() * std::abs(b - a);

            // return tempRestart->adjust(tmp);
            return tmp;
        }

        return(temp);
    }
    
}; // SATemperatureBandCooling


/**
 * No Cooling - constant temperature
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
