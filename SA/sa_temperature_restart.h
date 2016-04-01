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
            status->temp_counter = 0;
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
            status->temp_counter = 0;
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
            status->temp_counter = 0;
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
            status->temp_counter = 0;
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
            status->temp_counter = 0;
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
            status->temp_counter = 0;
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
            status->temp_counter = 0;
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
            status->temp_counter = 0;
            return temp / value;
        }
        return temp;
    }

    virtual int getTenure(void) {
        return tenure;
    }

}; // SALocalMinReheat


// e.g. meller-bozer
class SALocalMinTempRestart: public SATempRestart {

protected:
    int tenure;

public:
    SALocalMinTempRestart(SAInitTemp *it,
                            int _tenure):
        tenure(_tenure),
        SATempRestart(SALOCALMINRESTARTBEST) { }

    virtual float adjust(float temp) {
        if (status->not_improved > tenure) {
            status->temp_restarts += 1;
            status->temp_counter = 0;
            return status->best_temp;
        }
        return temp;
    }

    virtual int getTenure(void) {
        return tenure;
    }

}; // SALocalMinTempRestart


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
            status->temp_counter = 0;
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
            status->temp_counter = 0;
            return temp / value;
        }
        return temp;
    }

    virtual int getTenure(void) {
        return tenure;
    }

}; // SALocalMinEnhancedReheat



class SAMaxItersTempRestart: public SATempRestart {

protected:
    int tenure;
    float init_temp;

public:
    SAMaxItersTempRestart(SAInitTemp *it,
                          int _tenure):
        tenure(_tenure),
        init_temp(it->get()),
        SATempRestart(SAMAXITERSTEMPRESTART) { }

    virtual float adjust(float temp) {
        if (status->temp_counter > tenure) {
            status->temp_counter = 0;
            status->temp_restarts += 1;
            return init_temp;
        }
        return temp;
    }

    virtual int getTenure(void) {
        return tenure;
    }

}; // SAMaxItersTempRestart


class SANeighSizeMaxItersTempRestart: public SATempRestart {

protected:
    int tenure;
    float init_temp;

public:
    SANeighSizeMaxItersTempRestart(SAInitTemp *it,
                                   emili::Neighborhood* neigh,
                                   float _coeff):
        tenure(_coeff * neigh->size()),
        init_temp(it->get()),
        SATempRestart(SANEIGHSIZEMAXITERSTEMPRESTART) { }

    virtual float adjust(float temp) {
        if (status->temp_counter > tenure) {
            status->temp_counter = 0;
            status->temp_restarts += 1;
            return init_temp;
        }
        return temp;
    }

    virtual int getTenure(void) {
        return tenure;
    }

}; // SANeighSizeMaxItersTempRestart


class SASquaredNeighSizeMaxItersTempRestart: public SATempRestart {

protected:
    int tenure;
    float init_temp;

public:
    SASquaredNeighSizeMaxItersTempRestart(SAInitTemp *it,
                                          emili::Neighborhood* neigh,
                                          float _coeff):
        tenure(_coeff * neigh->size()),
        init_temp(it->get()),
        SATempRestart(SASQUAREDNEIGHSIZEMAXITERSTEMPRESTART) {
            double n = ceil(sqrt(2 * neigh->size()));
            tenure = (int)_coeff * n * n;
        }

    virtual float adjust(float temp) {
        if (status->temp_counter > tenure) {
            status->temp_counter = 0;
            status->temp_restarts += 1;
            // std::cout << status->temp_counter << " " << status->temp_restarts << " " << init_temp << "  *****" << std::endl;
            return init_temp;
        }
            // std::cout << status->temp_counter << " " << status->temp_restarts << " " << temp << std::endl;
        return temp;
    }

    virtual int getTenure(void) {
        return tenure;
    }

}; // SASquaredNeighSizeMaxItersTempRestart


class SAMaxItersReheat: public SATempRestart {

protected:
    int tenure;
    float alpha;

public:
    SAMaxItersReheat(SAInitTemp *it,
                          int _tenure,
                          float _alpha):
        tenure(_tenure),
        alpha(_alpha),
        SATempRestart(SAMAXITERSREHEAT) { }

    virtual float adjust(float temp) {
        if (status->temp_counter > tenure) {
            status->temp_counter = 0;
            status->temp_restarts += 1;
            return temp / alpha;
        }
        return temp;
    }

    virtual int getTenure(void) {
        return tenure;
    }

}; // SAMaxItersReheat


class SANeighSizeMaxItersReheat: public SATempRestart {

protected:
    int tenure;
    float alpha;

public:
    SANeighSizeMaxItersReheat(SAInitTemp *it,
                              emili::Neighborhood* neigh,
                              float _coeff,
                              float _alpha):
        tenure(_coeff * neigh->size()),
        alpha(_alpha),
        SATempRestart(SANEIGHSIZEMAXITERSREHEAT) { }

    virtual float adjust(float temp) {
        if (status->temp_counter > tenure) {
            status->temp_counter = 0;
            status->temp_restarts += 1;
            return temp / alpha;
        }
        return temp;
    }

    virtual int getTenure(void) {
        return tenure;
    }

}; // SANeighSizeMaxItersReheat


class SAMaxStepsTempRestart: public SATempRestart {

protected:
    int tenure;
    float init_temp;

public:
    SAMaxStepsTempRestart(SAInitTemp *it,
                          int _tenure):
        tenure(_tenure),
        init_temp(it->get()),
        SATempRestart(SAMAXSTEPSTEMPRESTART) { }

    virtual float adjust(float temp) {
        if (status->step > tenure) {
            status->step = 0;
            status->temp_restarts += 1;
            return init_temp;
        }
        return temp;
    }

    virtual int getTenure(void) {
        return tenure;
    }

}; // SAMaxStepsTempRestart


class SANeighSizeMaxStepsTempRestart: public SATempRestart {

protected:
    int tenure;
    float init_temp;

public:
    SANeighSizeMaxStepsTempRestart(SAInitTemp *it,
                                   emili::Neighborhood* neigh,
                                   float _coeff):
        tenure(_coeff * neigh->size()),
        init_temp(it->get()),
        SATempRestart(SANEIGHSIZEMAXSTEPSTEMPRESTART) { }

    virtual float adjust(float temp) {
        if (status->step > tenure) {
            status->step = 0;
            status->temp_restarts += 1;
            return init_temp;
        }
        return temp;
    }

    virtual int getTenure(void) {
        return tenure;
    }

}; // SANeighSizeMaxStepsTempRestart


class SAMaxStepsReheat: public SATempRestart {

protected:
    int tenure;
    float alpha;

public:
    SAMaxStepsReheat(SAInitTemp *it,
                          int _tenure,
                          float _alpha):
        tenure(_tenure),
        alpha(_alpha),
        SATempRestart(SAMAXSTEPSREHEAT) { }

    virtual float adjust(float temp) {
        if (status->step > tenure) {
            status->step = 0;
            status->temp_restarts += 1;
            return temp / alpha;
        }
        return temp;
    }

    virtual int getTenure(void) {
        return tenure;
    }

}; // SAMaxStepsReheat


class SANeighSizeMaxStepsReheat: public SATempRestart {

protected:
    int tenure;
    float alpha;

public:
    SANeighSizeMaxStepsReheat(SAInitTemp *it,
                              emili::Neighborhood* neigh,
                              float _coeff,
                              float _alpha):
        tenure(_coeff * neigh->size()),
        alpha(_alpha),
        SATempRestart(SANEIGHSIZEMAXSTEPSREHEAT) { }

    virtual float adjust(float temp) {
        if (status->step > tenure) {
            status->step = 0;
            status->temp_restarts += 1;
            return temp / alpha;
        }
        return temp;
    }

    virtual int getTenure(void) {
        return tenure;
    }

}; // SANeighSizeMaxStepsReheat


#endif
