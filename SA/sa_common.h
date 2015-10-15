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

    bool keep_last;

    std::string tc_type; // termination criterion
    std::string ac_type; // acceptance criterion
    std::string tl_type; // temperature length
    std::string tr_type; // temperature restart

    SAStatus(void) {
        counter = 0;
        total_counter = 0;
        accepted = 0;
        index = 0;
        rate = 1.;
        tenure = 0;

        keep_last = false;
    }


    void set_types(std::string _tc_type,
                   std::string _ac_type,
                   std::string _tl_type,
                   std::string _tr_type) {

        tc_type = _tc_type;
        ac_type = _ac_type;
        tl_type = _tl_type;
        tr_type = _tr_type;

        if (tc_type == LASTACCRATETERM   ||
            tr_type == SALASTRATERESTART ||
            tr_type == SALASTRATEREHEAT    ) {
            keep_last = true;
        }
    }


    void accepted_sol(void) {
        accepted += 1;
        counter = 0;

        if (keep_last) {
            last_accepted[index] = 1;
            index = (index + 1) % tenure;
        }
    }

    void not_accepted_sol(void) {
        if (keep_last) {
            last_accepted[index] = 0;
            index = (index + 1) % tenure;
        }
    }

}; // SAStatus;

#endif
