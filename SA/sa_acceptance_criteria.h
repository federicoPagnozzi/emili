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
    int   interval;
    int   counter;
    float rate;
    float alpha;

public:

    /**
     * constructors
     *
     * @param  initial_temperature initial temperature
     * @param  final_temperature   final_temperature
     * @param  descending_ratio    ratio of descent
     * @param  iterations          number of iterations (optional)
     * @param  alpha               alpha (optional)
     */
    SAAcceptance(float initial_temperature,
                 float final_temperature,
                 float descending_ratio):
                temperature(initial_temperature),
                start_temp(initial_temperature),
                end_temp(final_temperature),
                rate(descending_ratio),
                interval(1),
                counter(0),
                alpha(1) {}
    SAAcceptance(float initial_temperature,
                 float final_temperature,
                 float descending_ratio,
                 int   iterations):
                temperature(initial_temperature),
                start_temp(initial_temperature),
                end_temp(final_temperature),
                rate(descending_ratio),
                interval(iterations),
                counter(0),
                alpha(1) { }
    SAAcceptance(float initial_temperature,
                 float final_temperature,
                 float descending_ratio,
                 int   iterations,
                 float alpha):
                temperature(initial_temperature),
                start_temp(initial_temperature),
                end_temp(final_temperature),
                rate(descending_ratio),
                interval(iterations),
                counter(0),
                alpha(alpha) { }

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
                           float final_temperature,
                           float descending_ratio):
                SAAcceptance(initial_temperature,
                             final_temperature,
                             descending_ratio) { }
    SAMetropolisAcceptance(float initial_temperature,
                           float final_temperature,
                           float descending_ratio,
                           int   iterations):
                SAAcceptance(initial_temperature,
                             final_temperature,
                             descending_ratio,
                             iterations) { }
    SAMetropolisAcceptance(float initial_temperature,
                           float final_temperature,
                           float descending_ratio,
                           int   iterations,
                           float alpha):
                SAAcceptance(initial_temperature,
                             final_temperature,
                             descending_ratio,
                             iterations,
                             alpha) { }

    virtual emili::Solution* accept(emili::Solution *current_solution,
                                    emili::Solution *new_solution);

}; // SAMetropolisAcceptance


class SABasicAcceptance: public SAAcceptance {
public:
    SABasicAcceptance(float initial_temperature,
                      float final_temperature,
                      float descending_ratio):
                SAAcceptance(initial_temperature,
                             final_temperature,
                             descending_ratio) { }
    SABasicAcceptance(float initial_temperature,
                      float final_temperature,
                      float descending_ratio,
                      int   iterations):
                SAAcceptance(initial_temperature,
                             final_temperature,
                             descending_ratio,
                             iterations) { }
    SABasicAcceptance(float initial_temperature,
                      float final_temperature,
                      float descending_ratio,
                      int   iterations,
                      float alpha):
                SAAcceptance(initial_temperature,
                             final_temperature,
                             descending_ratio,
                             iterations,
                             alpha) { }

    virtual emili::Solution* accept(emili::Solution *current_solution,
                                    emili::Solution *new_solution);
    
}; // SABasicAcceptance

#endif
