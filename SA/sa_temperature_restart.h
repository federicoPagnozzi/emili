#ifndef SA_TEMP_RESTART_H
#define SA_TEMP_RESTART_H


#include "sa_init_temp.h"
#include "sa_common.h"


class SATempRestart {

protected:
    SAStatus* status;
    std::string type;

public:
    SATempRestart(std::string type):
        type(type) { }

    virtual float adjust(float temp)=0;

    void set_status(SAStatus* _status) {
        status = _status;
    }

    std::string getType() {
        return type;
    }

    virtual int getTenure(void) {
        return 0;
    }

}; // SaTempRestart



class SANoRestart: public SATempRestart {
public:
    SANoRestart(void):
        SATempRestart(SANOTEMPRESTART) { }

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
        SATempRestart(SAMINTEMPRESTART) { }

    virtual float adjust(float temp) {
        if (temp <= reset_threshold) {
            status->temp_restarts += 1;
            return init_temp;
        }
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
        SATempRestart(SAPERCTEMPRESTART) { }

    virtual float adjust(float temp) {
        if (temp <= reset_threshold) {
            status->temp_restarts += 1;
            return init_temp;
        }
        return temp;
    }

}; // SAPercRestart



class SALowRateRestart: public SATempRestart {

protected:
    float rate_threshold;
    float init_temp;

public:
    SALowRateRestart(SAInitTemp *it,
                     float _value):
        rate_threshold(_value),
        init_temp(it->get()),
        SATempRestart(SALOWRATERESTART) { }

    virtual float adjust(float temp) {
        if ((1.0 * status->accepted / status->total_counter) < rate_threshold) {
            status->temp_restarts += 1;
            return init_temp;
        }
        return temp;
    }

}; // SALowRateRestart


class SALowRateRestartToBest: public SATempRestart {

protected:
    float rate_threshold;
    float init_temp;

public:
    SALowRateRestartToBest(SAInitTemp *it,
                           float _value):
        rate_threshold(_value),
        init_temp(it->get()),
        SATempRestart(SALOWRATERESTARTBEST) { }

    virtual float adjust(float temp) {
        if ((1.0 * status->accepted / status->total_counter) < rate_threshold) {
            status->temp_restarts += 1;
            return status->best_temp;
        }
        return temp;
    }

}; // SALowRateRestartToBest


class SALastRateRestart: public SATempRestart {

protected:
    float rate_threshold;
    int tenure;
    float init_temp;


    int total_accepted(void) {
        int i, tot = 0;
        for (i = 0 ; i < status->tenure ; i++) {
            tot += status->last_accepted[i];
        }
        return(tot);
    }

public:
    SALastRateRestart(SAInitTemp *it,
                      int _tenure,
                      float _value):
        tenure(_tenure),
        rate_threshold(_value),
        init_temp(it->get()),
        SATempRestart(SALASTRATERESTART) { }

    virtual float adjust(float temp) {
        if ((1.0 * total_accepted() / status->tenure) < rate_threshold) {
            status->temp_restarts += 1;
            return init_temp;
        }
        return temp;
    }

    virtual int getTenure(void) {
        return tenure;
    }

}; // SALastRateRestart


/**
 * 0 < value < 1
 * Abramson, Krishnamoorthy, Dand - SA Cooling Schedules for the School Timetabling Problem
 */
class SALowRateReheat: public SATempRestart {

protected:
    float rate_threshold;
    float value;

public:
    SALowRateReheat(SAInitTemp *it,
                    float threshold,
                    float _value):
        rate_threshold(threshold),
        value(_value),
        SATempRestart(SALOWRATEREHEAT) { }

    virtual float adjust(float temp) {
        if ((1.0 * status->accepted / status->total_counter) < rate_threshold) {
            status->temp_restarts += 1;
            return temp / value;
        }
        return temp;
    }

}; // SALowRateReheat


/**
 * 0 < value < 1
 * Abramson, Krishnamoorthy, Dand - SA Cooling Schedules for the School Timetabling Problem
 */
class SALastRateReheat: public SATempRestart {

protected:
    int tenure;
    float rate_threshold;
    float value;


    int total_accepted(void) {
        int i, tot = 0;
        for (i = 0 ; i < status->tenure ; i++) {
            tot += status->last_accepted[i];
        }
        return(tot);
    }

public:
    SALastRateReheat(SAInitTemp *it,
                     int _tenure,
                     float threshold,
                     float _value):
        tenure(_tenure),
        rate_threshold(threshold),
        value(_value),
        SATempRestart(SALASTRATEREHEAT) { }

    virtual float adjust(float temp) {
        if ((1.0 * total_accepted() / status->tenure) < rate_threshold) {
            status->temp_restarts += 1;
            return temp / value;
        }
        return temp;
    }

    virtual int getTenure(void) {
        return tenure;
    }

}; // SALastRateReheat


/**
 * 0 < value < 1
 * Abramson, Krishnamoorthy, Dand - SA Cooling Schedules for the School Timetabling Problem
 */
class SALocalMinReheat: public SATempRestart {

protected:
    int tenure;
    float value;

public:
    SALocalMinReheat(SAInitTemp *it,
                     int _tenure,
                     float _value):
        tenure(_tenure),
        value(_value),
        SATempRestart(SALOCALMINREHEAT) { }

    virtual float adjust(float temp) {
        if (status->not_improved > tenure) {
            status->temp_restarts += 1;
            return temp / value;
        }
        return temp;
    }

    virtual int getTenure(void) {
        return tenure;
    }

}; // SALocalMinReheat


class SALocalMinRestartToBest: public SATempRestart {

protected:
    int tenure;

public:
    SALocalMinRestartToBest(SAInitTemp *it,
                            int _tenure):
        tenure(_tenure),
        SATempRestart(SALOCALMINRESTARTBEST) { }

    virtual float adjust(float temp) {
        if (status->not_improved > tenure) {
            status->temp_restarts += 1;
            return status->best_temp;
        }
        return temp;
    }

    virtual int getTenure(void) {
        return tenure;
    }

}; // SALocalMinRestartToBest


/**
 * 0 < value < 1
 * Abramson, Krishnamoorthy, Dand - SA Cooling Schedules for the School Timetabling Problem
 *
 * modified to avoid value to take non-positive values
 */
class SALocalMinEnhancedReheat: public SATempRestart {

protected:
    int tenure;
    float value;
    float epsilon;

public:
    SALocalMinEnhancedReheat(SAInitTemp *it,
                             int _tenure,
                             float _value,
                             float _epsilon):
        tenure(_tenure),
        value(_value),
        epsilon(_epsilon),
        SATempRestart(SALOCALMINENHANCEDREHEAT) { }

    virtual float adjust(float temp) {
        if (status->not_improved > tenure) {
            value = std::max(epsilon, value - epsilon);
            status->temp_restarts += 1;
            return temp / value;
        }
        return temp;
    }

    virtual int getTenure(void) {
        return tenure;
    }

}; // SALocalMinEnhancedReheat


#endif
