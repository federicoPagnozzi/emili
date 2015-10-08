#ifndef SA_TERMINATION_CRITERIA_H
#define SA_TERMINATION_CRITERIA_H

#include <string>
#include <iostream>
#include <vector>
#include <functional>

#include <string>

#include "sa_common.h"

#include "../emilibase.h"


class SATermination: public emili::Termination {

protected:
     std::string type;

public:
    SATermination(std::string type):
        type(type) { }

    virtual bool terminate(emili::Solution* currentSolution,
                           emili::Solution* newSolution)=0;

    virtual bool terminate(sa_status *status)=0;

    virtual void reset()=0;

    std::string getType() {
        return type;
    }

    virtual int getTenure(void) {
        return 0;
    }

}; // SATermination


class SAMaxBadIterTermination: public SATermination {

protected:
    int maxBadIterations;

public:
    SAMaxBadIterTermination(int maxBadIterations):
        maxBadIterations(maxBadIterations),
        SATermination(MAXBADITERS) { }

    virtual bool terminate(emili::Solution* currentSolution,
                           emili::Solution* newSolution) {
        return true;
    }

    virtual bool terminate(sa_status *status) {
        if (status->counter >= maxBadIterations) return true;
        return false;
    }

    virtual void reset() {
        // counter = 0;
    }

}; // SAMaxBadIterTermination


class SAMaxIterTermination: public SATermination {

protected:
    int maxIterations;

public:
    SAMaxIterTermination(int maxIterations):
        maxIterations(maxIterations),
        SATermination(MAXITERS) { }

    virtual bool terminate(emili::Solution* currentSolution,
                           emili::Solution* newSolution) {
        return true;
    }

    virtual bool terminate(sa_status *status) {
        if (status->counter >= maxIterations) return true;
        return false;
    }

    // does nothing
    virtual void reset() {
        // counter = 0;
    }

}; // SAMaxIterTermination


class SAWhileTrueTermination: public SATermination {

public:
    SAWhileTrueTermination(void):
        SATermination(NEVERTERM) { }

    virtual bool terminate(emili::Solution* currentSolution,
                           emili::Solution* newSolution) {
        return false;
    }

    virtual bool terminate(sa_status *status) {
        return false;
    }

    // does nothing
    virtual void reset() {
        // counter = 0;
    }

}; // SAWhileTrueTermination



class SAAcceptanceRateTermination: public SATermination {

protected:
    float rate;

public:
    SAAcceptanceRateTermination(float _rate):
        rate(_rate),
        SATermination(ACCRATETERM) { }

    virtual bool terminate(emili::Solution* currentSolution,
                           emili::Solution* newSolution) {
        return false;
    }

    virtual bool terminate(sa_status *status) {
        if ((1.0 * status->accepted / status->total_counter) < rate)
            return true;
        
        return false;
    }

    // does nothing
    virtual void reset() {
        // counter = 0;
    }

}; // SAAcceptanceRateTermination



class SALastAcceptanceRateTermination: public SATermination {

protected:
    float rate;
    int   tenure;

    int total_accepted(sa_status *status) {
        int i, tot = 0;
        for (i = 0 ; i < status->tenure ; i++) {
            tot += status->last_accepted[i];
        }
        return(tot);
    }

public:
    SALastAcceptanceRateTermination(int   _tenure,
                                    float _rate):
        tenure(_tenure),
        rate(_rate),
        SATermination(LASTACCRATETERM) { }

    virtual bool terminate(emili::Solution* currentSolution,
                           emili::Solution* newSolution) {
        return false;
    }

    virtual bool terminate(sa_status *status) {
        if ((1.0 * total_accepted(status) / status->tenure) < rate)
            return true;
        
        return false;
    }

    // does nothing
    virtual void reset() {
        // counter = 0;
    }

    virtual int getTenure(void) {
        return tenure;
    }

}; // SALastAcceptanceRateTermination


#endif
