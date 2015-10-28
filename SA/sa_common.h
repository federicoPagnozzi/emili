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
    int    curr_accepted;
    float  rate;
    short *last_accepted;
    int    tenure;
    int    index;
    int    step;

    int    not_improved;

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
        curr_accepted = 0;
        index = 0;
        rate = 1.;
        tenure = 0;
        not_improved = 0;
        step = 0;

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

        if (tc_type == LASTACCRATETERM          ||
            tr_type == SALASTRATERESTART        ||
            tr_type == SALASTRATEREHEAT         ||
            tr_type == SALOCALMINENHANCEDREHEAT   ) {
            keep_last = true;
        }
    }


    void increment_counters(void) {
        counter += 1;
        total_counter += 1;
        not_improved += 1;
    }


    void new_best_solution(emili::Solution* sol, float cost) {
        delete best;
        best = sol->clone();
        best_cost = cost;
        not_improved = 0;
    }


    void accepted_sol(float cost) {
        accepted += 1;
        curr_accepted += 1;
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
