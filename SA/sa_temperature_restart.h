#ifndef SA_TEMP_RESTART_H
#define SA_TEMP_RESTART_H


#include "sa_init_temp.h"


class SATempRestart {

protected:
    float value;

public:
    SATempRestart(float _value):
        value(_value) { }

    virtual float restartAt(void) {
        return value;
    }

    void setValue(float _value) {
        value = _value;
    }

}; // SaTempRestart



class SANoRestart: public SATempRestart {
public:
    SANoRestart(void):
        SATempRestart(-1) { }

}; // SANoRestart


/**
 * temperature delta
 */
class SADeltaRestart: public SATempRestart {
public:
    SADeltaRestart(float _value):
        SATempRestart(_value) { }

}; // SADeltaRestart


/**
 * temperature percentage
 */
class SAPercRestart: public SATempRestart {
public:
    SAPercRestart(SAInitTemp *it, float _value):
        SATempRestart(_value * it->get() / 100.0) { }

}; // SAPercRestart

#endif
