#ifndef SA_TERMINATION_CRITERIA_H
#define SA_TERMINATION_CRITERIA_H

#include <string>
#include <iostream>
#include <vector>
#include <functional>

#include <string>

#include "../emilibase.h"


class SATermination: public emili::Termination {

protected:
     std::string type;

public:
    SATermination(std::string type):
        type(type) { }

    virtual bool terminate(emili::Solution* currentSolution,
                           emili::Solution* newSolution)=0;

    virtual bool terminate(int counter)=0;

    virtual void reset()=0;

    std::string getType() {
        return type;
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

    virtual bool terminate(int counter) {
        if (counter >= maxBadIterations) return true;
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

    virtual bool terminate(int counter) {
        if (counter >= maxIterations) return true;
        return false;
    }

    // does nothing
    virtual void reset() {
        // counter = 0;
    }

}; // SAMaxIterTermination


class SAWhileTrueTermination: public SATermination {

protected:
    int maxIterations;

public:
    SAWhileTrueTermination(void):
        SATermination(NEVERTERM) { }

    virtual bool terminate(emili::Solution* currentSolution,
                           emili::Solution* newSolution) {
        return false;
    }

    virtual bool terminate(int counter) {
        return false;
    }

    // does nothing
    virtual void reset() {
        // counter = 0;
    }

}; // SAWhileTrueTermination


#endif
