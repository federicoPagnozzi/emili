#ifndef SA_COOLING_H
#define SA_COOLING_H

#include "../emilibase.h"

#include "sa_common.h"
#include "sa_init_temp.h"
#include "sa_templength.h"
#include "sa_temperature_restart.h"

#include "../pfsp/pfspinstance.h"


/**
 * Generic cooling scheme for SA.
 *
 * The actual scheme has to be implemented with derived classes.
 */
class SACooling {

protected:
    double a;
    double b;
    long counter;
    double inittemp;

    int search_time_seconds;

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

    virtual void set_status(SAStatus* _status) {
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

    void setA(double _a) {
        a = _a;
    }

    void setB(double _b) {
        b = _b;
    }

    virtual void set_search_time(int ts) {
        search_time_seconds = ts;
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

    MarkovCooling(double b, SAInitTemp *it):
        SACooling(it->get(), b, it) { }

    virtual double update_cooling(double temp) {
        counter++;

        if (tempLength->isCoolingTime(counter)) {
            counter = 0;
            status->step = status->step + 1;
            float tmp = (a / (b + std::log(status->step)));// instead of log(temp) + b

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
    LogCooling(double b, SAInitTemp *it):
        SACooling(it->get(), b, it) { }

    virtual double update_cooling(double temp) {
        counter++;

        if (tempLength->isCoolingTime(counter)) {
            counter = 0;
            status->step = status->step + 1;
            float tmp = (a / std::log(status->step + b)); // instead of temp + b

            return tempRestart->adjust(tmp);

        }

        return(temp);
    }

}; // LogCooling


/**
 * NO... http://www.sciencedirect.com/science/article/pii/037722179090301Q#
 *
 * Szu, Hartley - fast simulated annealing
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
 * for testing Connolly settings
 * 
 */
class LundyMeesConnollyCooling: public SACooling {


public:
    LundyMeesConnollyCooling(SAInitTemp *it):
        SACooling(0, 0, it) { }

    void set_search_time(int ts) {
        search_time_seconds = ts;

        double estimated_no_of_moves = search_time_seconds / status->move_time;
        b = (status->init_temp - status->final_temp) /
                (estimated_no_of_moves * status->init_temp * status->final_temp);
    }

    virtual double update_cooling(double temp) {
        counter++;

        b = (status->init_temp - status->final_temp) /
                (status->neigh_size * 50 * status->init_temp * status->final_temp);

        if (tempLength->isCoolingTime(counter)) {
            counter = 0;
            status->step = status->step + 1;
            float tmp = (temp / (1 + b*temp));

            return tempRestart->adjust(tmp);

        }

        return(temp);
    }

}; // LundyMeesConnollyCooling


/**
 * Linear cooling
 * 0 < a < 1
 *
 * temp' = a*temp
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
            double tmp = a*temp;

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
 * Andersen, Vidal, Iversen
 * Design of a teleprocessing communication network using SA
 *
 * T_f = 0
 */
class SAQuadraticCooling: public SACooling {
protected:
    float c;
public:
    SAQuadraticCooling(SAInitTemp *it):
        c(it->get()),
        SACooling(0, 0, it) {
            a = -inittemp; // T_i - T_f
            b = 2 * inittemp; // 2 * (T_f - T_i)
        }

    void setTempLength(SATempLength* _tempLength) {
        tempLength = _tempLength;
        a = a / (tempLength->getLength() * tempLength->getLength());
        b = b / tempLength->getLength();
    }

    virtual double update_cooling(double temp) {
        counter++;

        if (tempLength->isCoolingTime(counter)) {
            counter = 0;
            status->step = status->step + 1;
            float tmp = a * (status->step * status->step) + b * status->step + c;

            return tempRestart->adjust(tmp);
        }

        return(temp);
    }

}; // SAQuadraticCooling


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


// connolly paper
class Q87Cooling: public SACooling {

protected:
    int max_reject;
    bool in_fixed_temp_state;

public:
    Q87Cooling(double a, double b, int _max_reject, SAInitTemp *it):
        max_reject(_max_reject),
        in_fixed_temp_state(false),
        SACooling(a, b, it) { }

    virtual double update_cooling(double temp) {
        counter++;

        if (status->force_accept) {
            status->force_accept = false;
        }

        if (!in_fixed_temp_state              &&
            tempLength->isCoolingTime(counter)  ) {

            if (status->counter > max_reject) {
                // stop cooling
                status->force_accept = true;
                setA(1);
                setB(0);
                in_fixed_temp_state = true;
                // std::cout << "stopped cooling at: " << status->best_temp << std::endl;
                return status->best_temp;
            }

            counter = 0;
            status->step = status->step + 1;
            float tmp = (temp / (a + b*temp));

            tmp = tempRestart->adjust(tmp);
            // std::cout << "cooling at: " << tmp << std::endl;
            return tmp;

        }

        return(temp);
    }

}; // Q87Cooling


// connolly paper
class ConnollyQ87Cooling: public SACooling {

protected:
    int max_reject;
    bool in_fixed_temp_state;

public:
    ConnollyQ87Cooling(SAInitTemp *it,
                       emili::Neighborhood *nei):
        max_reject(nei->size()),
        in_fixed_temp_state(false),
        SACooling(1, 0, it) { }

    void set_search_time(int ts) {
        search_time_seconds = ts;

        double estimated_no_of_moves = search_time_seconds / status->move_time;
        b = (status->init_temp - status->final_temp) /
                (estimated_no_of_moves * status->init_temp * status->final_temp);
                
    }

    void set_status(SAStatus* _status) {
        status = _status;

        b = (status->init_temp - status->final_temp) /
                (status->neigh_size * 50 * status->init_temp * status->final_temp);
    }

    virtual double update_cooling(double temp) {
        counter++;
        if (status->force_accept) {
            status->force_accept = false;
        }

        if (!in_fixed_temp_state) {

            if (status->counter > max_reject) {
                // stop cooling
                status->force_accept = true;
                setA(1);
                setB(0);
                in_fixed_temp_state = true;
                // std::cout << "stopped cooling at: " << status->best_temp << std::endl;
                return status->best_temp;
            }

            counter = 0;
            status->step = status->step + 1;
            
            float tmp = (temp / (1 + b*temp));

            tmp = tempRestart->adjust(tmp);
            // std::cout << "cooling at: " << tmp << std::endl;
            return tmp;

        }

        return(temp);
    }

}; // ConnollyQ87Cooling


// osman potts
class OsmanPottsPFSPCooling: public SACooling {
protected:
emili::pfsp::PermutationFlowShop *instance;
double K;

public:
    OsmanPottsPFSPCooling(SAInitTemp *_it,
                          emili::pfsp::PermutationFlowShop *_instance):
                        instance(_instance),
                        SACooling(1, 0, _it) {
                            K = std::max(3300 * std::log(_instance->getNjobs()) + 7500 * std::log(_instance->getNmachines()) - 18250, 2000.0);
                        }

    void set_status(SAStatus* _status) {
        status = _status;

        b = (status->init_temp - status->final_temp) /
                ((K-1) * status->init_temp * status->final_temp);
    }

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

}; // OsmanPottsPFSPCooling

/**
 * Arithmetic Cooling for GDA
 */
class ArithmeticCooling: public SACooling {

public:
    ArithmeticCooling(double a, SAInitTemp *it):
        SACooling(a, 0, it) { }

    virtual double update_cooling(double temp) {
        counter++;
        
        if (tempLength->isCoolingTime(counter)) {
            counter = 0;
            status->step = status->step + 1;
            double tmp = std::max(temp - a, 0.0);

            return tempRestart->adjust(tmp);

        }

        return(temp);
    }

}; // ArithmeticCooling



/**
 * Old Bachelor Acceptance - it's actually a nonmonotonic cooling scheme
 * Old bachelor acceptance: A new class of non-monotone threshold accepting methods
 * Hu, Te C and Kahng, Andrew B and Tsao, Chung-Wen Albert
 */
/**
 * version 1
 */
class OldBachelor1: public SACooling {

protected:
    long M;
    double delta;
    double a;
    double b;
    double c;
    emili::Problem *prob;

public:
    OldBachelor1 (long _M,
                  double _delta,
                  double _a,
                  double _b,
                  double _c,
                  SAInitTemp *it,
                  emili::Problem* _prob):
                  M(_M),
                  delta(_delta),
                  a(_a),
                  b(_b),
                  c(_c),
                  prob(_prob),
                  SACooling(0, 0, it) {
                    a = std::lrint(a * prob->problemSize());
                  }

    virtual double update_cooling(double temp) {
        counter++;

        if (tempLength->isCoolingTime(counter)) {
            counter = 0;
            status->step = status->step + 1;

            // status->counter == age variable
            double tmp = (std::pow(status->counter / a, b) - 1) * delta;
            tmp = tmp * std::pow(std::fmax(0, 1 - (status->total_counter / M)), c); 

            return tempRestart->adjust(tmp);

        }

        return(temp);
    }

}; // OldBachelor1

/**
 * version 2
 */
class OldBachelor2: public SACooling {

protected:
    long M;
    double delta;
    double d;
    emili::Neighborhood* neigh;
    long age_count;

public:
    OldBachelor2 (long _M,
                  double _delta,
                  double _d,
                  SAInitTemp *it,
                  emili::Neighborhood* _neigh):
                  M(_M),
                  delta(_delta),
                  neigh(_neigh),
                  SACooling(0, 0, it) {
                    d = std::sqrt(_d * neigh->size());
                    age_count = 1;
                  }

    virtual double update_cooling(double temp) {
        counter++;

        if (tempLength->isCoolingTime(counter)) {
            counter = 0;
            status->step = status->step + 1;

            // status->counter == age variable
            double tmp = temp;
            if (status->counter == 0) {
                if (status->counter < d) {
                    age_count++;
                } else {
                    age_count = 1;
                }
                tmp = tmp - age_count * delta * (1 - status->total_counter / M);
            } else {
                age_count++;
                tmp = tmp + (delta / d) * (1 - status->total_counter / M);
            } 

            return tempRestart->adjust(tmp);

        }

        return(temp);
    }

}; // OldBachelor2

#endif
