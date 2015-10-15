#ifndef SA_TEMP_RESTART_H
#define SA_TEMP_RESTART_H


#include "sa_init_temp.h"


class SATempRestart {

protected:
    float value;

public:
    SATempRestart(void) { }

    virtual float adjust(float temp)=0;

}; // SaTempRestart



class SANoRestart: public SATempRestart {
public:
    SANoRestart(void):
        SATempRestart( ) { }

    virtual float adjust(float temp) {
        return temp;
    }

}; // SANoRestart


/**
 * temperature restart when reaches a minimum
 */
class SAMinRestart: public SATempRestart {

protected:
    float reset_threshold;
    float init_temp;

public:
    SAMinRestart(SAInitTemp *it,
                 float _value):
        reset_threshold(_value),
        init_temp(it->get()),
        SATempRestart( ) { }

    virtual float adjust(float temp) {
        if (temp <= reset_threshold)
            return init_temp;
        return temp;
    }

}; // SADeltaRestart


/**
 * temperature percentage
 */
class SAPercRestart: public SATempRestart {

protected:
    float reset_threshold;
    float init_temp;

public:
    SAPercRestart(SAInitTemp *it,
                  float _value):
        reset_threshold(_value * it->get() / 100.0),
        init_temp(it->get()),
        SATempRestart() { }

    virtual float adjust(float temp) {
        return value;
    }

}; // SAPercRestart


#endif
