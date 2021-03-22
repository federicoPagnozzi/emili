#ifndef SA_COMMON_H
#define SA_COMMON_H


#include "../emilibase.h"

#include "sa_constants.h"


#define nullptr NULL

// a reference to SearchStatus sstatus can be changed to SAStatus like this
// SAStatus& ss = dynamic_cast<SAStatus&> (sstatus);
class SAStatus : public emili::SearchStatus{

public:
    emili::Solution* best;
    /*long counter;
    long total_counter;
    long not_improved;*/
    long   local_counter;
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

    double move_time;

    double best_cost;
    double best_temp;

    double init_cost;

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

    SAStatus(void) :emili::SearchStatus() {
        local_counter = 0;
        counter = 0;
        total_counter = 0;
        not_improved = 0;
        temp_counter = 0;
        accepted = 0;
        curr_accepted = 0;
        index = 0;
        rate = 1.;
        tenure = 0;        
        step = 0;
        best_temp = 0;
        neigh_size = 0;
        temp_restarts = 0;
        init_prob = 1.0;
        move_time = 0.0;
        keep_last = false;
        force_accept = false;
        init_cost = 0.0;
        best = nullptr;
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


    virtual void increment_counters(void) {
        counter += 1;
        total_counter += 1;
        temp_counter += 1;
        not_improved += 1;
        local_counter += 1;
    }

    void new_best_solution_silent(emili::Solution* sol, double cost, double temp) {
        //delete best;
        best = sol->clone();
        newBestSolution(best);
        best_cost = cost;
        best_temp = temp;
        //std::cout << std::fixed << "New best solution found: " << best->getSolutionRepresentation();
        //std::cout << std::fixed << "of cost " << cost << " at iteration " << total_counter << std::endl;
        //fprintf(stdout, "NEWBEST %f %ld %f %f\n", cost, total_counter, temp, emili::getCurrentExecutionTime());
        //fflush(stdout);
        /*std::cout << std::fixed;
        std::cout << cost;
        std::cout << " ";
        std::cout << total_counter;
        std::cout << " ";
        std::cout << temp;
        std::cout << " ";
        std::cout << emili::getCurrentExecutionTime();
        std::cout << std::endl;
        std::cout << std::flush;*/
        counter = 0;
    }

    void new_best_solution(emili::Solution* sol, double cost, double temp) {
        //delete best;
        best = sol->clone();
        newBestSolution(best);
        best_cost = cost;
        best_temp = temp;
        //std::cout << std::fixed << "New best solution found: " << best->getSolutionRepresentation();
        //std::cout << std::fixed << "of cost " << cost << " at iteration " << total_counter << std::endl;
        fprintf(stdout, "NEWBEST %f %ld %f %f\n", cost, total_counter, temp, emili::getCurrentExecutionTime());
        fflush(stdout);
        /*std::cout << std::fixed;
        std::cout << cost;
        std::cout << " ";
        std::cout << total_counter;
        std::cout << " ";
        std::cout << temp;
        std::cout << " ";
        std::cout << emili::getCurrentExecutionTime();
        std::cout << std::endl;
        std::cout << std::flush;*/
        counter = 0;
    }


    void accepted_sol(double cost) {
        accepted += 1;
        curr_accepted += 1;
        counter = 0;
        //std::cout << "accepted " << counter << std::endl;

        if (keep_last) {
            last_accepted[index] = 1;
            index = (index + 1) % tenure;
        }
    }

    void not_accepted_sol(void) {
        //std::cout << "NOT" << std::endl;
        if (keep_last) {
            last_accepted[index] = 0;
            index = (index + 1) % tenure;
        }
    }

    void print(void) {
        std::cout << "\nSA Status: \n"
                  << "  initial temperature : " << init_temp << "\n"
                  << "  final temperature :   " << final_temp << "\n"
                  << "  current temperature : " << temp << "\n"
                  << "  initial probability : " << init_prob << "\n\n"
                  << "  total # of moves :    " << total_counter << "\n"
                  << "  # of accepted moves : " << accepted << "\n"
                  << "  # of steps :          " << step << "\n"
                  << "  # of temp restarts :  " << temp_restarts << "\n"
                  << "  neighbourhood size :  " << neigh_size << "\n\n"
                  << "  best cost :           " << best_cost << "\n"
                  << "  best solution temp :  " << best_temp << "\n" << std::endl;
    }

    virtual void resetCounters() {
      accepted = 0;
      curr_accepted = 0;
      counter = 0;
      total_counter = 0;
      local_counter = 0;
      not_improved = 0;
      temp_counter = 0;
      temp_restarts = 0;
      if (last_accepted != nullptr) {
        free(last_accepted);
      }
      step = 0;
      index = 0;
      last_accepted = (short *)malloc(tenure * sizeof(short));      
      for (int i = 0 ; i < tenure ; i++) {
        last_accepted[i] = 1;
      }
    }

    void resetCounters(long _total_c) {
      this->resetCounters();
      total_counter = _total_c;
    }

    /*emili::Solution* getBestSolution(void) {
       return getBestSolution();
    }*/

}; // SAStatus;

#endif
