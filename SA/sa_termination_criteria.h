#ifndef SA_TERMINATION_CRITERIA_H
#define SA_TERMINATION_CRITERIA_H

#include <string>
#include <iostream>
#include <vector>
#include <functional>

#include <string>

#include "sa_common.h"

#include "../emilibase.h"
namespace emili {
namespace sa{

class SATermination: public emili::Termination {

protected:
     std::string type;

public:
    SATermination(std::string type):
        type(type) { }

    virtual bool terminate(emili::Solution* currentSolution,
                           emili::Solution* newSolution)=0;

    virtual bool terminate(SAStatus& status)=0;

    virtual void reset()=0;

    std::string getType() {
        return type;
    }

    virtual void setMinTemp(double value) {}

    virtual void setNeighborhoodsize(int neighborhood_size) {}

    virtual int getTenure(void) {
        return 0;
    }

}; // SATermination


class SAMaxBadIterTermination: public SATermination {

protected:
    long maxBadIterations;

public:
    SAMaxBadIterTermination(int maxBadIterations):
        maxBadIterations(maxBadIterations),
        SATermination(MAXBADITERS) { }

    virtual bool terminate(emili::Solution* currentSolution,
                           emili::Solution* newSolution) {
        return true;
    }

    virtual bool terminate(SAStatus& status) {
        if (status.counter >= maxBadIterations) return true;
        return false;
    }

    virtual void reset() {
        // counter = 0;
    }

}; // SAMaxBadIterTermination



class SAMaxIterTermination: public SATermination {

protected:
    long maxIterations;

public:
    SAMaxIterTermination(int maxIterations):
        maxIterations(maxIterations),
        SATermination(MAXITERS) { }

    virtual bool terminate(emili::Solution* currentSolution,
                           emili::Solution* newSolution) {
        return true;
    }

    virtual bool terminate(SAStatus& status) {
        if (status.total_counter >= maxIterations) return true;
        return false;
    }

    // does nothing
    virtual void reset() {
        // counter = 0;
    }

}; // SAMaxIterTermination


class SAMaxLocalIterTermination: public SATermination {

protected:
    long maxIterations;

public:
    SAMaxLocalIterTermination(int maxIterations):
        maxIterations(maxIterations),
        SATermination(MAXLOCALITERS) { }

    virtual bool terminate(emili::Solution* currentSolution,
                           emili::Solution* newSolution) {
        return true;
    }

    virtual bool terminate(SAStatus& status) {
        if (status.local_counter >= maxIterations) return true;
        return false;
    }

    // does nothing
    virtual void reset() {
        // counter = 0;
    }

}; // SAMaxIterTermination




class SANeighSizeIterTermination: public SATermination {

protected:
    double coeff;
    emili::Neighborhood *neigh;
    long maxIterations;

public:
    SANeighSizeIterTermination(double _coeff):
        neigh(nullptr),
        coeff(_coeff),
        maxIterations(coeff),
        SATermination(NEIGHSIZEITERTERM) { }

    SANeighSizeIterTermination(emili::Neighborhood *nei, double _coeff):
        neigh(nei),
        coeff(_coeff),
        maxIterations(neigh->size() * coeff),
        SATermination(NEIGHSIZEITERTERM) { }

    virtual bool terminate(emili::Solution* currentSolution,
                           emili::Solution* newSolution) {
        return true;
    }

    virtual bool terminate(SAStatus& status) {
        if (status.total_counter >= maxIterations) return true;
        return false;
    }

    // does nothing
    virtual void reset() {
        // counter = 0;
    }

    virtual void setNeighborhoodsize(int neighborhood_size)
    {
        maxIterations = coeff*neighborhood_size;
    }

}; // SANeighSizeIterTermination


// Burkard-Rendl
class SASquaredNeighSizeIterTermination: public SATermination {

protected:
    double coeff;
    emili::Neighborhood *neigh;
    long maxIterations;

public:
    SASquaredNeighSizeIterTermination(double _coeff):
        coeff(_coeff),
        SATermination(SQUAREDNSITERTERM) {
            double n = ceil(sqrt(2 * neigh->size()));
            maxIterations = _coeff * n * n;
        }

    SASquaredNeighSizeIterTermination(emili::Neighborhood *nei, double _coeff):
        neigh(nei),
        coeff(_coeff),
        SATermination(SQUAREDNSITERTERM) {
            double n = ceil(sqrt(2 * neigh->size()));
            maxIterations = _coeff * n * n;
        }

    virtual bool terminate(emili::Solution* currentSolution,
                           emili::Solution* newSolution) {
        return true;
    }

    virtual bool terminate(SAStatus& status) {
        if (status.total_counter >= maxIterations) return true;
        return false;
    }

    // does nothing
    virtual void reset() {
        // counter = 0;
    }

    virtual void setNeighborhoodsize(int neighborhood_size)
    {
        double n = ceil(sqrt(2 * neighborhood_size));
        maxIterations = coeff * n * n;
    }

}; // SASquaredNeighSizeIterTermination


class SAWhileTrueTermination: public SATermination {

public:
    SAWhileTrueTermination(void):
        SATermination(NEVERTERM) { }

    virtual bool terminate(emili::Solution* currentSolution,
                           emili::Solution* newSolution) {
        if(emili::getKeepGoing()) {
            return false;
        } else {
            return true;
        }
    }

    virtual bool terminate(SAStatus& status) {
        if(emili::getKeepGoing()) {
            return false;
        } else {
            return true;
        }
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

    virtual bool terminate(SAStatus& status) {
        if ((1.0 * status.accepted / status.total_counter) < rate)
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

    int total_accepted(SAStatus& status) {
        int i, tot = 0;
        for (i = 0 ; i < status.tenure ; i++) {
            tot += status.last_accepted[i];
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

    virtual bool terminate(SAStatus& status) {
        if ((1.0 * total_accepted(status) / status.tenure) < rate)
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


/**
 * meller-bozer, new simulated annealing algorithm for facility layout
 */
class SaMaxTempRestartsTermination: public SATermination {

protected:
    int mtr;

public:
    SaMaxTempRestartsTermination(int _mtr):
        mtr(_mtr),
        SATermination(MAXTEMPRESTARTSTERM) { }

    virtual bool terminate(emili::Solution* currentSolution,
                           emili::Solution* newSolution) {
        return false;
    }

    virtual bool terminate(SAStatus& status) {
        if (status.temp_restarts > mtr)
            return true;
        
        return false;
    }

    virtual void reset() {
        // counter = 0;
    }

}; // SaMaxTempRestartsTermination


/**
 * ker yan tam - a simulated annealing for allocating to manifacturing cells
 */
class SALocalMinTermination: public SATermination {

protected:
    int tenure;

public:
    SALocalMinTermination(int _tenure):
        tenure(_tenure),
        SATermination(LOCALMINTERM) { }

    virtual bool terminate(emili::Solution* currentSolution,
                           emili::Solution* newSolution) {
        return false;
    }

    virtual bool terminate(SAStatus& status) {
        if (status.not_improved > tenure) {
            return true;
        }
        return false;
    }

    virtual void reset() {
        // counter = 0;
    }

}; // SALocalMinTermination


/**
 * ker yan tam - a simulated annealing for allocating to manifacturing cells
 */
class SANeighSizeLocalMinTermination: public SATermination {

protected:
    int tenure;
    float coeff;
public:
    SANeighSizeLocalMinTermination(float _coeff):
        tenure(1),
        coeff(_coeff),
        SATermination(NEIGHSIZELOCALMINTERM) { }

    SANeighSizeLocalMinTermination(emili::Neighborhood *neigh, float _coeff):
        tenure((int)(_coeff * neigh->size())),
        SATermination(NEIGHSIZELOCALMINTERM) { }

    virtual bool terminate(emili::Solution* currentSolution,
                           emili::Solution* newSolution) {
        return false;
    }

    virtual bool terminate(SAStatus& status) {
        if (status.not_improved > tenure) {
            return true;
        }
        return false;
    }

    virtual void reset() {
        // counter = 0;
    }

    virtual void setNeighborhoodsize(int neighborhood_size)
    {
        tenure = coeff*neighborhood_size;
    }

}; // SANeighSizeLocalMinTermination


/**
 * jajodia et al - class, computerized layout solution using simulated annealing
 */
class SAMaxStepsTermination: public SATermination {

protected:
    int maxsteps;

public:
    SAMaxStepsTermination(long _ms):
        maxsteps(_ms),
        SATermination(MAXSTEPSTERM) { }

    virtual bool terminate(emili::Solution* currentSolution,
                           emili::Solution* newSolution) {
        return false;
    }

    virtual bool terminate(SAStatus& status) {
        if (status.step > maxsteps) {
            return true;
        }
        return false;
    }

    virtual void reset() {
        // counter = 0;
    }

}; // SAMaxStepsTermination


/**
 */
class SAMinTempTermination: public SATermination {

protected:
    double mintemp;

public:
    SAMinTempTermination():
        mintemp(0),
        SATermination(MINTEMPTERM) { }

    SAMinTempTermination(double _mt):
        mintemp(_mt),
        SATermination(MINTEMPTERM) { }

    virtual bool terminate(emili::Solution* currentSolution,
                           emili::Solution* newSolution) {
        return false;
    }

    virtual bool terminate(SAStatus& status) {
        if (status.temp <= mintemp) {
            return true;
        }
        return false;
    }

    virtual void reset() {
        // counter = 0;
    }

    virtual void setMinTemp(double value)
    {
        mintemp = value;
    }

}; // SAMinTempTermination

}
}

#endif
