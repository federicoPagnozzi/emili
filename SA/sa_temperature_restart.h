#ifndef SA_TEMP_RESTART_H
#define SA_TEMP_RESTART_H


#include "sa_init_temp.h"
#include "sa_cooling.h"


class SATempRestart {

protected:
    float value;

public:
    SATempRestart(float _value,
                  SACooling& cooling):
        value(_value) {
            cooling.setResetThreshold(value);
        }

    virtual float restartAt(void) {
        return value;
    }

    void setValue(float _value) {
        value = _value;
    }

}; // SaTempRestart



class SANoRestart: public SATempRestart {
public:
    SANoRestart(SACooling& cooling):
        SATempRestart(-1, cooling) { }

}; // SANoRestart


/**
 * temperature restart when reaches a minimum
 */
class SAMinRestart: public SATempRestart {
public:
    SAMinRestart(float _value,
                 SACooling& cooling):
        SATempRestart(_value, cooling) { }

}; // SADeltaRestart


/**
 * temperature percentage
 */
class SAPercRestart: public SATempRestart {
public:
    SAPercRestart(SAInitTemp *it,
                  float _value,
                  SACooling& cooling):
        SATempRestart(_value * it->get() / 100.0, cooling) { }

}; // SAPercRestart


class SAAcceptanceRestart: public SATempRestart {

}; // SAAcceptanceRestart

#endif
