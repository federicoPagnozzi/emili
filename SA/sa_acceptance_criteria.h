#ifndef SA_ACCEPTANCE_H
#define SA_ACCEPTANCE_H

#include <string>
#include <iostream>
#include <vector>
#include <functional>
#include <cfloat>
#include <stdexcept>

#include "sa_constants.h"
#include "sa_common.h"

#include "../emilibase.h"

/**
 * Acceptance criteria for Simulated Annealing
 */
class SAAcceptance: public emili::Acceptance
{

protected:

    float temperature;
    float start_temp;

    std::string type;

    SAStatus* status;

public:

    /**
     * constructors
     *
     * @param  initial_temperature initial temperature
     * @param  final_temperature   final_temperature
     */
    SAAcceptance(std::string type,
                 float initial_temperature):
                type(type),
                temperature(initial_temperature),
                start_temp(initial_temperature),
                status(nullptr) { }

    /**
     * Acceptance method
     * 
     * @param  current_emili::Solution emili::Solution 1
     * @param  new_emili::Solution     emili::Solution 2
     * @return                  accepted solution (must be one of the two ?)
     */
    virtual emili::Solution* accept(emili::Solution *current_solution,
                                    emili::Solution *new_solution)=0;

    /**
     * @brief acceptViaDelta
     * @param new_solution the solution modified
     * @param delta = cost(new_solution) - cost(before)
     * @return true iff the change is accepted
     */
    virtual bool acceptViaDelta(emili::Solution* new_solution, double delta) {
        throw std::invalid_argument("not implemented");
    }


    /**
     * update temperature
     * @param  temp new temperature
     */
    virtual void setCurrentTemp(float temp) {
        temperature = temp;
    }

    std::string getType() {
        return type;
    }

    void set_status(SAStatus* _status) {
        status = _status;
    }

}; // class SAAcceptance


class SAMetropolisAcceptance: public SAAcceptance {
public:
    SAMetropolisAcceptance(float initial_temperature):
                SAAcceptance(METROPOLIS,
                             initial_temperature) { }

    virtual emili::Solution* accept(emili::Solution *current_solution,
                                    emili::Solution *new_solution);

    virtual bool acceptViaDelta(emili::Solution *new_solution, double delta) override;
}; // SAMetropolisAcceptance


// connolly paper
class SAMetropolisWithForcedAcceptance: public SAAcceptance {
public:
    SAMetropolisWithForcedAcceptance(float initial_temperature):
                SAAcceptance(METROPOLISWFORCED,
                             initial_temperature) { }

    virtual emili::Solution* accept(emili::Solution *current_solution,
                                    emili::Solution *new_solution);

}; // SAMetropolisWithForcedAcceptance


class SAApproxExpAcceptance: public SAAcceptance {
public:
    SAApproxExpAcceptance(float initial_temperature):
                SAAcceptance(APPROXEXPACC,
                             initial_temperature) { }

    virtual emili::Solution* accept(emili::Solution *current_solution,
                                    emili::Solution *new_solution);

}; // SAApproxExpAcceptance


class SABasicAcceptance: public SAAcceptance {
public:
    SABasicAcceptance(void):
                SAAcceptance(BASICACC,
                             0) { }

    virtual emili::Solution* accept(emili::Solution *current_solution,
                                    emili::Solution *new_solution);
    
}; // SABasicAcceptance


/**
 * Bohachevsky-Johnson-Stein, Generalized Simulated Annealing for function optimization
 */
class GeneralizedSAAcceptance: public SAAcceptance {

protected:
    float g;

public:
    GeneralizedSAAcceptance(float initial_temperature,
                            float g):
                SAAcceptance(GENSAACC,
                             initial_temperature) { }

    virtual emili::Solution* accept(emili::Solution *current_solution,
                                    emili::Solution *new_solution);

}; // GeneralizedSAAcceptance


/**
 * http://www.sciencedirect.com/science/article/pii/030505489090001N#
 *
 * Ogbu - Smith
 * 
 * 0 < initial_acceptance < 1
 */
class SAGeometricAcceptance: public SAAcceptance {

protected:
    float rate;
    int step;

public:
    SAGeometricAcceptance(float initial_acceptance,
                          float reducing_factor):
                rate(reducing_factor),
                SAAcceptance(GEOMACC,
                             initial_acceptance) { }

    virtual emili::Solution* accept(emili::Solution *current_solution,
                                    emili::Solution *new_solution);

    void setStep(int newStep) {
        step = newStep;
    }

    int getStep(void) {
        return step;
    }

}; // SAGeometricAcceptance



class SADeterministicAcceptance: public SAAcceptance {

protected:
    float delta;

public:
    SADeterministicAcceptance(float _delta):
                delta(_delta),
                SAAcceptance(DETERMINISTICACC,
                             0) { }

    virtual emili::Solution* accept(emili::Solution *current_solution,
                                    emili::Solution *new_solution);

}; // SADeterministicAcceptance


/**
 * Dueck - Great deluge
 */
class GreatDelugeAcceptance: public SAAcceptance {

public:
    GreatDelugeAcceptance(void):
        SAAcceptance(GDAACC,
                     0) { }

    virtual emili::Solution* accept(emili::Solution *current_solution,
                                    emili::Solution *new_solution);

}; // GreatDelugeAcceptance


/**
 * Dueck - Record-to-Record travel
 *
 * acceptance is based NOT on record + deviation
 * but on temperature * (1 + deviation/100)
 * 0 < deviation  <= 100
 */
class RecordToRecordAcceptance: public SAAcceptance {

protected:
    float deviation;

public:
    RecordToRecordAcceptance(float _deviation):
        deviation(_deviation),
        SAAcceptance(RTRACC,
                     0) { }

    virtual emili::Solution* accept(emili::Solution *current_solution,
                                    emili::Solution *new_solution);

}; // RecordToRecordAcceptance


/**
 * Burke-Bykov, late acceptance Hill climbing
 */
class LAHCAcceptance: public SAAcceptance {

protected:
    int    tenure;
    float *cost_list;

public:
    LAHCAcceptance(int _tenure):
        tenure(_tenure),
        SAAcceptance(LAHCACC,
                     0) {
            cost_list = (float *)malloc(sizeof(float) * tenure);
            for (int i = 0 ; i < tenure ; i++) {
                cost_list[i] = FLT_MAX;
            }
        }

    virtual emili::Solution* accept(emili::Solution *current_solution,
                                    emili::Solution *new_solution);

    ~LAHCAcceptance(void) {
        free(cost_list);
    }

};  // LAHCAcceptance

#endif
