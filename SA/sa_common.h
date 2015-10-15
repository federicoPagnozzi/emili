#ifndef SA_COMMON_H
#define SA_COMMON_H


#include "../emilibase.h"

#include "sa_constants.h"


#define nullptr NULL


class SAStatus {

public:
    int    counter;
    int    total_counter;
    int    accepted;
    float  rate;
    short *last_accepted;
    int    tenure;
    int    index;

    emili::Solution *best;
    float best_cost;

    SAStatus(void) {
        counter = 0;
        total_counter = 0;
        accepted = 0;
        index = 0;
        rate = 1.;
        tenure = 0;
    }

}; // SAStatus;

#endif
