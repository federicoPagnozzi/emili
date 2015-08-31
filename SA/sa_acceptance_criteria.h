#ifndef SA_ACCEPTANCE_H
#define SA_ACCEPTANCE_H

#include <string>
#include <iostream>
#include <vector>
#include <functional>

#include "sa_constants.h"
#include "../emilibase.h"

/**
 * Acceptance criteria for Simulated Annealing
 */
class SAAcceptance: public emili::Acceptance
{

protected:

    float temperature;
    float start_temp;
    float end_temp;

    std::string type;

public:

    /**
     * constructors
     *
     * @param  initial_temperature initial temperature
     * @param  final_temperature   final_temperature
     */
    SAAcceptance(std::string type,
                 float initial_temperature,
                 float final_temperature):
                type(type),
                temperature(initial_temperature),
                start_temp(initial_temperature),
                end_temp(final_temperature) { }

    /**
     * Acceptance method
     * 
     * @param  current_emili::Solution emili::Solution 1
     * @param  new_emili::Solution     emili::Solution 2
     * @return                  accepted solution
     */
    virtual emili::Solution* accept(emili::Solution *current_solution,
                                    emili::Solution *new_solution)=0;

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

}; // class SAAcceptance


class SAMetropolisAcceptance: public SAAcceptance {
public:
    SAMetropolisAcceptance(float initial_temperature,
                           float final_temperature):
                SAAcceptance(METROPOLIS,
                             initial_temperature,
                             final_temperature) { }

    virtual emili::Solution* accept(emili::Solution *current_solution,
                                    emili::Solution *new_solution);

}; // SAMetropolisAcceptance


class SABasicAcceptance: public SAAcceptance {
public:
    SABasicAcceptance(float initial_temperature,
                      float final_temperature):
                SAAcceptance(BASICACC,
                             initial_temperature,
                             final_temperature) { }

    virtual emili::Solution* accept(emili::Solution *current_solution,
                                    emili::Solution *new_solution);
    
}; // SABasicAcceptance



/**
 * http://www.sciencedirect.com/science/article/pii/030505489090001N#
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
                             initial_acceptance,
                             0) { }

    virtual emili::Solution* accept(emili::Solution *current_solution,
                                    emili::Solution *new_solution);

    void setStep(int newStep) {
        step = newStep;
    }

    int getStep(void) {
        return step;
    }

}; // SAGeometricAcceptance

#endif
