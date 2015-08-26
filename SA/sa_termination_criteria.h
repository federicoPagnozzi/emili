#ifndef SA_TERMINATION_CRITERIA_H
#define SA_TERMINATION_CRITERIA_H

#include <string>
#include <iostream>
#include <vector>
#include <functional>

#include "../emilibase.h"


class SATermination: public emili::Termination {

public:
	virtual bool terminate(emili::Solution* currentSolution,
		                   emili::Solution* newSolution)=0;

	virtual bool terminate(int counter)=0;

	virtual void reset()=0;

}; // SATermination


class SAMaxBadIterTermination: public SATermination {

protected:
	int maxBadIterations;

public:
    SAMaxBadIterTermination(int maxBadIterations):
        maxBadIterations(maxBadIterations) { }

    virtual bool terminate(emili::Solution* currentSolution,
		                   emili::Solution* newSolution) {
    	return true;
    }

	virtual bool terminate(int counter) {
		if (counter >= maxBadIterations) return true;
		return false;
	}

	virtual void reset() {
		maxBadIterations = 0;
	}

}; // SAMaxBadIterTermination


#endif
