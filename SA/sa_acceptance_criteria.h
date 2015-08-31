#ifndef SA_ACCEPTANCE_H
#define SA_ACCEPTANCE_H

#include <string>
#include <iostream>
#include <vector>
#include <functional>

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

public:

    /**
     * constructors
     *
     * @param  initial_temperature initial temperature
     * @param  final_temperature   final_temperature
     */
    SAAcceptance(float initial_temperature,
                 float final_temperature):
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

}; // class SAAcceptance


class SAMetropolisAcceptance: public SAAcceptance {
public:
    SAMetropolisAcceptance(float initial_temperature,
                           float final_temperature):
                SAAcceptance(initial_temperature,
                             final_temperature) { }

    virtual emili::Solution* accept(emili::Solution *current_solution,
                                    emili::Solution *new_solution);

}; // SAMetropolisAcceptance


class SABasicAcceptance: public SAAcceptance {
public:
    SABasicAcceptance(float initial_temperature,
                      float final_temperature):
                SAAcceptance(initial_temperature,
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
    float: rate;

public:
    SAGeometricAcceptance(float initial_acceptance,
                          float reducing_factor):
                rate(reducing_factor),
                SAAcceptance(initial_acceptance,
                             0) { }

    virtual emili::Solution* accept(emili::Solution *current_solution,
                                    emili::Solution *new_solution);

}; // SAGeometricAcceptance

#endif
