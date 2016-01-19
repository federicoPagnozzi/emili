#ifndef SA_COMMON_H
#define SA_COMMON_H


#include "../emilibase.h"

#include "sa_constants.h"


#define nullptr NULL


class SAStatus {

public:
    long   counter;
    long   total_counter;
    long   temp_counter;
    long   accepted;
    long   curr_accepted;
    double  rate;
    short *last_accepted;
    long   tenure;
    long   index;
    long   step;
    int    temp_restarts;

    long   neigh_size;

    long   not_improved;

    emili::Solution *best;
    double best_cost;
    double best_temp;

    double init_temp;
    double final_temp;
    double temp;

    double init_prob;

    bool keep_last;
    bool force_accept;

    std::string tc_type; // termination criterion
    std::string ac_type; // acceptance criterion
    std::string tl_type; // temperature length
    std::string tr_type; // temperature restart

    SAStatus(void) {
        counter = 0;
        total_counter = 0;
        temp_counter = 0;
        accepted = 0;
        curr_accepted = 0;
        index = 0;
        rate = 1.;
        tenure = 0;
        not_improved = 0;
        step = 0;
        best_temp = 0;
        neigh_size = 0;
        temp_restarts = 0;
        init_prob = 1.0;

        keep_last = false;
        force_accept = false;
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
        temp_counter += 1;
        not_improved += 1;
    }


    void new_best_solution(emili::Solution* sol, double cost, double temp) {
        delete best;
        best = sol->clone();
        best_cost = cost;
        best_temp = temp;
        //std::cout << std::fixed << "New best solution found: " << best->getSolutionRepresentation();
        //std::cout << std::fixed << "of cost " << cost << " at iteration " << total_counter << std::endl;
        std::cout << std::fixed << cost << " " << total_counter << " ";// << temp << " ";
        std::cout << std::fixed << emili::getCurrentExecutionTime() << std::endl;
        not_improved = 0;
        counter = 0;
    }


    void accepted_sol(double cost) {
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
