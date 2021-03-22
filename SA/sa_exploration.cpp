#include "sa_exploration.h"
using namespace emili::sa;
emili::Solution* SARandomExploration::nextSolution(emili::Solution *startingSolution,
                                                   SAStatus& status) {

    bool noneaccepted = true;

    status.increment_counters();

    //emili::Solution* startingSolutionBackup = startingSolution->clone();
    emili::Solution* incumbent = neigh->random(startingSolution);
    emili::Solution* accepted = acceptance->accept(startingSolution,
                                                   incumbent);
    if (accepted == startingSolution) {
        //std::cout << "randomexpl, delete incumbent line 19 " << std::endl;
        //std::cout << incumbent << std::endl;
        delete incumbent;
        status.not_accepted_sol();
    } else {
        //std::cout << "randomexpl, delete startingSolution line 14 " << std::endl;
        /*delete startingSolution;
        startingSolution = nullptr;*/
        status.accepted_sol(accepted->getSolutionValue());
        noneaccepted = false;
        //std::cout << "startingSolution " << status.counter << std::endl;
    }
    //startingSolution = accepted;

    status.temp = cooling->update_cooling(status.temp);
    acceptance->setCurrentTemp(status.temp);

    if (noneaccepted) {
        return startingSolution;
    }

    return accepted;

    //return accepted;
}


emili::Solution* SASequentialExploration::nextSolution(emili::Solution *startingSolution,
                                                       SAStatus& status) {


    emili::Solution* incumbent = startingSolution->clone();
    emili::Solution* accepted;
    emili::Solution* ithSolution;

    //neigh->reset();
    emili::Neighborhood::NeighborhoodIterator iter = neigh->begin(incumbent);
    // std::cout << incumbent->getSolutionRepresentation() << " " << incumbent->getSolutionValue() << std::endl;
    ithSolution = *iter;

    bool noneaccepted = true;

    for(;
        iter!=neigh->end();
        ++iter) {

        status.increment_counters();
        accepted = acceptance->accept(incumbent,
                                      ithSolution);

        //std::cout << accepted->getSolutionRepresentation() << " " << accepted->getSolutionValue() << std::endl;
        if (accepted == incumbent) {
            status.not_accepted_sol();
        } else {
            //delete startingSolution;
            *accepted = *ithSolution;
            status.accepted_sol(accepted->getSolutionValue());
            noneaccepted = false;
            break;
        }
        
        status.temp = cooling->update_cooling(status.temp);
        acceptance->setCurrentTemp(status.temp);


    } 

    delete incumbent;
    
    if (noneaccepted) {
        return startingSolution;
    }
    
    return accepted;
}


emili::Solution* SABestOfKExploration::nextSolution(emili::Solution *startingSolution,
                                                    SAStatus& status) {

    double ci, cg;

    emili::Solution* incumbent = neigh->random(startingSolution);
    ci = incumbent->getSolutionValue();

    for (long i = 1 ; i < k ; i++) {
        status.increment_counters();

        emili::Solution* generated = neigh->random(startingSolution);
        cg = generated->getSolutionValue();

        if (ci < cg) {
            delete generated;
        } else {
            delete incumbent;
            incumbent = generated;
            ci = cg;
        }
    }

    emili::Solution* accepted = acceptance->accept(startingSolution,
                                                   incumbent);
    status.temp = cooling->update_cooling(status.temp);
    acceptance->setCurrentTemp(status.temp);

    if (accepted == startingSolution) {
        delete incumbent;
        status.not_accepted_sol();
    } else {
        //delete startingSolution;
        status.accepted_sol(accepted->getSolutionValue());
    }
    //startingSolution = accepted;

    return accepted;//startingSolution;

}

emili::Solution* SABestOfKSequentialExploration::nextSolution(emili::Solution *startingSolution,
                                                    SAStatus& status) {

    double ci, cg;

    long i = 1;

    emili::Solution* incumbent = startingSolution->clone();
    emili::Solution* accepted;
    emili::Solution* generated;
    emili::Solution* candidate = incumbent;

    ci = incumbent->getSolutionValue();

    emili::Neighborhood::NeighborhoodIterator iter = neigh->begin(incumbent);
    generated = *iter;

    bool noneaccepted = true;

    for(;
        iter!=neigh->end() && i < k;
        ++iter) {

        status.increment_counters();

        cg = generated->getSolutionValue();

        if (cg < ci) {
            candidate = generated;
            ci = cg;
        }
        
        i++;

    } 

    accepted = acceptance->accept(incumbent,
                                  candidate);
    status.temp = cooling->update_cooling(status.temp);
    acceptance->setCurrentTemp(status.temp);


    if (accepted == incumbent) {
        status.not_accepted_sol();
    } else {
        //delete startingSolution;
        *accepted = *candidate;
        status.accepted_sol(accepted->getSolutionValue());
        noneaccepted = false;
    }

    delete incumbent;
    
    if (noneaccepted) {
        return startingSolution;
    }
    
    return accepted;

}

emili::Solution* SANSBestOfKSequentialExploration::nextSolution(emili::Solution *startingSolution,
                                                    SAStatus& status) {

    double ci, cg;
    long i = 0;

    // k CANNOT EXCEED THE NEIGHBOURHOOD SIZE!!!
    k = std::min(k, neighsize);
    
    // stats
    double orig_ci; // original cost of initial solution 
    // number of improving, worsening, neutral moves
    long num_better = 0,
         num_worse = 0,
         num_equal = 0;
    double gap;
    double gap_sum = 0.0; // gap sum
    double gaps[k]; // all the gaps
    double maxgap, mingap; // max gap, min gap, ça va sans dire


    emili::Solution* incumbent = startingSolution->clone();
    emili::Solution* accepted;
    emili::Solution* generated;
    emili::Solution* candidate = startingSolution->clone();

    //printf("COST OF STARTING: %f\n", startingSolution->getSolutionValue());

    ci = incumbent->getSolutionValue();

    orig_ci = ci;
    maxgap = 0.0;
    mingap = ci;

    emili::Neighborhood::NeighborhoodIterator iter = neigh->begin(incumbent);
    generated = *iter;

    bool noneaccepted = true;

    for(;
        iter!=neigh->end() &&
        i < k;
        ++iter) {

        status.increment_counters();

        cg = generated->getSolutionValue();

        if (cg < ci) {
            delete candidate;
            candidate = generated->clone();
            ci = cg;
        }

        //stats
        if (cg < orig_ci) {
            num_better++;
        } else if (cg == orig_ci) {
            num_equal++;
        } else {
            num_worse++;
        }
        gap = abs(cg - orig_ci);
        gaps[i] = gap;
        if (gap < mingap && gap > 0) mingap = gap;
        if (gap > maxgap) maxgap = gap;
        gap_sum += gap;
        
        //printf("%ld %f\n", i, cg);
        
        i++;

    } 

    if (1 || status.total_counter % 100 == 1) {
        // print:
        // orig_ci, final_ci, %better, %equal, %worse, mingap, avggap, maxgap, stddevgap, 
        double avggap = gap_sum / k;
        double stddevgap = 0.0, tmpstd;
        for (long j = 0 ; j < k ; j++) {
            tmpstd = gaps[j] - avggap;
            stddevgap = stddevgap + tmpstd * tmpstd;
        }
        stddevgap = sqrt(stddevgap / (k - 1));
        fprintf(stdout, "RUNTIMESTATS %f %f %f %f %f %f %f %f %f\n",
            orig_ci, ci,
            (num_better * 1.0) / k,
            (num_equal * 1.0) / k,
            (num_worse * 1.0) / k,
            mingap, avggap, maxgap, stddevgap
            );
        fflush(stdout);
    }

    //for (long j = 0 ; j < 100 ; j++) {
        //printf("COST OF WTF: %f\n", accepted->getSolutionValue());
        //accepted = neigh->random(accepted);
        /*accepted = acceptance->accept(candidate,
                                      neigh->random(incumbent));*/
    /**/accepted = acceptance->accept(incumbent,
                                  candidate);/**/
    //}
    status.temp = cooling->update_cooling(status.temp);
    acceptance->setCurrentTemp(status.temp);


    if (accepted == incumbent) {
        status.not_accepted_sol();
    } else {
        //delete startingSolution;
        //*accepted = *candidate;
        status.accepted_sol(accepted->getSolutionValue());
        noneaccepted = false;
    }

    delete incumbent;
    
    if (noneaccepted) {
        return startingSolution;
    }

    //printf("COST OF ACCEPTED: %f\n", accepted->getSolutionValue());    

    return accepted;

}

emili::Solution* SANSBestOfKRandomExploration::nextSolution(emili::Solution *startingSolution,
                                                    SAStatus& status) {
    double ci, cg;
    long i = 0;

    // k CANNOT EXCEED THE NEIGHBOURHOOD SIZE!!!
    k = std::min(k, neighsize);
    
    // stats
    double orig_ci; // original cost of initial solution 
    // number of improving, worsening, neutral moves
    long num_better = 0,
         num_worse = 0,
         num_equal = 0;
    double gap;
    double gap_sum = 0.0; // gap sum
    double gaps[k]; // all the gaps
    double maxgap, mingap; // max gap, min gap, ça va sans dire

    emili::Solution* incumbent = startingSolution->clone();
    emili::Solution* accepted;
    emili::Solution* generated;
    emili::Solution* candidate = startingSolution->clone();

    //printf("COST OF STARTING: %f\n", startingSolution->getSolutionValue());

    ci = incumbent->getSolutionValue();

    orig_ci = ci;
    maxgap = 0.0;
    mingap = ci;

    emili::Neighborhood::NeighborhoodIterator iter = neigh->begin(incumbent);
    generated = *iter;

    bool noneaccepted = true;

    for(;
        iter!=neigh->end() &&
        i < k;
        ++iter) {

        status.increment_counters();

        cg = generated->getSolutionValue();

        if (cg < ci) {
            delete candidate;
            candidate = generated->clone();
            ci = cg;
        }

        //stats
        if (cg < orig_ci) {
            num_better++;
        } else if (cg == orig_ci) {
            num_equal++;
        } else {
            num_worse++;
        }
        gap = abs(cg - orig_ci);
        gaps[i] = gap;
        if (gap < mingap && gap > 0) mingap = gap;
        if (gap > maxgap) maxgap = gap;
        gap_sum += gap;
        
        i++;

    } 

    if (1 || status.total_counter % 100 == 1) {
        // print:
        // orig_ci, final_ci, %better, %equal, %worse, mingap, avggap, maxgap, stddevgap, 
        double avggap = gap_sum / k;
        double stddevgap = 0.0, tmpstd;
        for (long j = 0 ; j < k ; j++) {
            tmpstd = gaps[j] - avggap;
            stddevgap = stddevgap + tmpstd * tmpstd;
        }
        stddevgap = sqrt(stddevgap / (k - 1));
        fprintf(stdout, "RUNTIMESTATS %f %f %f %f %f %f %f %f %f\n",
            orig_ci, ci,
            (num_better * 1.0) / k,
            (num_equal * 1.0) / k,
            (num_worse * 1.0) / k,
            mingap, avggap, maxgap, stddevgap
            );
        fflush(stdout);
    }

    //for (long j = 0 ; j < 100 ; j++) {
        //printf("COST OF WTF: %f\n", accepted->getSolutionValue());
        //accepted = neigh->random(accepted);
        /**/accepted = acceptance->accept(startingSolution,
                                      neigh->random(startingSolution));/**/
    /** /accepted = acceptance->accept(incumbent,
                                  candidate);/ **/
    //}
    status.temp = cooling->update_cooling(status.temp);
    acceptance->setCurrentTemp(status.temp);

    if (accepted == startingSolution) {
        status.not_accepted_sol();
    } else {
        //delete startingSolution;
        //*accepted = *candidate;
        status.accepted_sol(accepted->getSolutionValue());
        noneaccepted = false;
    }

    delete incumbent;

    if (noneaccepted) {
        return startingSolution;
    }

    //printf("COST OF ACCEPTED: %f\n", accepted->getSolutionValue());    

    return accepted;

}


emili::Solution* SAFirstBestOfKExploration::nextSolution(emili::Solution *startingSolution,
                                                         SAStatus& status) {

    double ci, cg, cs = startingSolution->getSolutionValue();

    emili::Solution* incumbent = neigh->random(startingSolution);
    ci = incumbent->getSolutionValue();

    for (long i = 1 ; i < k ; i++) {
        status.increment_counters();

        if (ci < cs) {
            //delete startingSolution;
            status.accepted_sol(incumbent->getSolutionValue());
            return incumbent;
        }

        emili::Solution* generated = neigh->random(startingSolution);
        cg = generated->getSolutionValue();

        if (ci < cg) {
            delete generated;
        } else {
            delete incumbent;
            incumbent = generated;
            ci = cg;
        }
    }

    emili::Solution* accepted = acceptance->accept(startingSolution,
                                                   incumbent);
    status.temp = cooling->update_cooling(status.temp);
    acceptance->setCurrentTemp(status.temp);

    if (accepted == startingSolution) {
        delete incumbent;
        status.not_accepted_sol();
    } else {
        //delete startingSolution;
        status.accepted_sol(accepted->getSolutionValue());
    }
    //startingSolution = accepted;

    return accepted;//startingSolution;

}


emili::Solution* SANSBestOfKExploration::nextSolution(emili::Solution *startingSolution,
                                                    SAStatus& status) {

    double ci, cg;

    emili::Solution* incumbent = neigh->random(startingSolution);
    ci = incumbent->getSolutionValue();

    for (long i = 1 ; i < k ; i++) {
        status.increment_counters();

        emili::Solution* generated = neigh->random(startingSolution);
        cg = generated->getSolutionValue();

        if (ci < cg) {
            delete generated;
        } else {
            delete incumbent;
            incumbent = generated;
            ci = cg;
        }
    }

    emili::Solution* accepted = acceptance->accept(startingSolution,
                                                   incumbent);
    status.temp = cooling->update_cooling(status.temp);
    acceptance->setCurrentTemp(status.temp);

    if (accepted == startingSolution) {
        delete incumbent;
        status.not_accepted_sol();
    } else {
        //delete startingSolution;
        status.accepted_sol(accepted->getSolutionValue());
    }
    //startingSolution = accepted;

    return accepted;//startingSolution;

}


emili::Solution* SANSFirstBestOfKExploration::nextSolution(emili::Solution *startingSolution,
                                                         SAStatus& status) {

    double ci, cg, cs = startingSolution->getSolutionValue();

    emili::Solution* incumbent = neigh->random(startingSolution);
    ci = incumbent->getSolutionValue();

    for (long i = 1 ; i < k ; i++) {
        status.increment_counters();

        if (ci < cs) {
            delete startingSolution;
            status.accepted_sol(incumbent->getSolutionValue());
            return incumbent;
        }

        emili::Solution* generated = neigh->random(startingSolution);
        cg = generated->getSolutionValue();

        if (ci < cg) {
            delete generated;
        } else {
            delete incumbent;
            incumbent = generated;
            ci = cg;
        }
    }

    emili::Solution* accepted = acceptance->accept(startingSolution,
                                                   incumbent);
    status.temp = cooling->update_cooling(status.temp);
    acceptance->setCurrentTemp(status.temp);

    if (accepted == startingSolution) {
        delete incumbent;
        status.not_accepted_sol();
    } else {
        //delete startingSolution;
        status.accepted_sol(accepted->getSolutionValue());
    }
    //startingSolution = accepted;

    return accepted;//startingSolution;

}
