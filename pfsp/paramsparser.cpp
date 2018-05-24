//
//  Created by Federico Pagnozzi on 28/11/14.
//  Copyright (c) 2014 Federico Pagnozzi. All rights reserved.
//  This file is distributed under the BSD 2-Clause License. See LICENSE.TXT
//  for details.
#include "paramsparser.h"
#include <cstdlib>
#include <cstdio>
#include <ctime>
#include <cstring>
#include <iostream>
#include <sstream>
#include <algorithm>
#include "pfspinstance.h"

/* Algos */
#define IG "ig"
#define ILS "ils"
#define TABU "tabu"
#define FIRST "first"
#define BEST "best"
#define TB_FIRST "tfirst"
#define TB_BEST "tbest"
#define VND "vnd"
#define GVNS_ILS "gvns"
#define CH6_LS "ch6"
#define TEST_INIT "stin"
#define EMPTY_LOCAL "nols"

#define ALBERTOSA "SA"

/* tabu tenure types */
#define TABU_MEMORY_MOVES "move"
#define TABU_MEMORY_HASHES "hash"
#define TABU_MEMORY_SOLUTIONS "solution"
#define TABU_MEMORY_TSAB "tsabm"
#define TABU_MEMORY_TSAB_TEST "tsabmt"
#define TABU_MEMORY_VALUE "value"

/* modifiers */
#define RO "-ro"
#define IT "-it"
#define TS "-ts"
#define TI "-ti"
#define IN "-in"
#define TE "-te"
#define NE "-ne"
#define NS "-ns"
#define RNDSEED "rnds"
#define DEFAULT_TS 10
#define DEFAULT_TI 10
#define DEFAULT_IT 0


/* Permutation flowshop*/
#define PROBLEM_PFS_WT "PFSP_WT"
#define PROBLEM_PFS_WE "PFSP_WE"
#define PROBLEM_PFS_TCT "PFSP_TCT"
#define PROBLEM_PFS_T "PFSP_T"
#define PROBLEM_PFS_E "PFSP_E"
#define PROBLEM_PFS_WCT "PFSP_WCT"
#define PROBLEM_PFS_MS "PFSP_MS"

/* no wait permutation flowshop*/
#define PROBLEM_NWPFS_MS "NWPFSP_MS"
#define PROBLEM_NWPFS_WT "NWPFSP_WT"
#define PROBLEM_NWPFS_WE "NWPFSP_WE"
#define PROBLEM_NWPFS_WCT "NWPFSP_WCT"
#define PROBLEM_NWPFS_TCT "NWPFSP_TCT"
#define PROBLEM_NWPFS_T "NWPFSP_T"
#define PROBLEM_NWPFS_E "NWPFSP_E"

/* no idle permutation flowshop*/
#define PROBLEM_NIPFS_MS "NIPFSP_MS"
#define PROBLEM_NIPFS_WT "NIPFSP_WT"
#define PROBLEM_NIPFS_WE "NIPFSP_WE"
#define PROBLEM_NIPFS_WCT "NIPFSP_WCT"
#define PROBLEM_NIPFS_TCT "NIPFSP_TCT"
#define PROBLEM_NIPFS_T "NIPFSP_T"
#define PROBLEM_NIPFS_E "NIPFSP_E"

/* Sequence dependent setup times */
#define PROBLEM_SDSTPFS_MS "SDSTPFS_MS"
#define PROBLEM_SDSTPFS_WT "SDSTPFS_WT"
#define PROBLEM_SDSTPFS_WE "SDSTPFS_WE"
#define PROBLEM_SDSTPFS_T "SDSTPFS_T"
#define PROBLEM_SDSTPFS_E "SDSTPFS_E"
#define PROBLEM_SDSTPFS_TCT "SDSTPFS_TCT"
#define PROBLEM_SDSTPFS_WCT "SDSTPFS_WCT"


/* initial solution heuristics */
#define INITIAL_NEH "neh"
#define INITIAL_NEHRS "nehrs"
#define INITIAL_NEHEDD "nehedd"
#define INITIAL_NEHFF "nehff"
#define INITIAL_NEHLS "nehls"
#define INITIAL_NEHEDDLS "neheddls"
#define INITIAL_NEHFFLS "nehffls"
#define INITIAL_RANDOM "random"
#define INITIAL_RANDOM_ITERATED "irandom"
#define INITIAL_SLACK "slack"
#define INITIAL_SRZ "srz"
#define INITIAL_LIT "lit"
#define INITIAL_RZ "rz"
#define INITIAL_NRZ "nrz"
#define INITIAL_NRZ2 "nrz2"
#define INITIAL_NRZ2FF "nrz2ff"
#define INITIAL_LR "lr"
#define INITIAL_LR_NM "lrnm"
#define INITIAL_NLR "nlr"
#define INITIAL_MNEH "mneh"
#define INITIAL_WNSLACK "nwslack"
#define INITIAL_FRB5 "frb5"

/* Termination criteria*/
#define TERMINATION_MAXSTEPS "maxstep"
#define TERMINATION_MAXSTEPS_OR_LOCMIN "msorlocmin"
#define TERMINATION_MAXSTEPS_OR_LOCMIN "msorlocmin"
#define TERMINATION_TIME "time"
#define TERMINATION_LOCMIN "locmin"
#define TERMINATION_ITERA "iteration"
#define TERMINATION_WTRUE "true"
#define TERMINATION_SOA "soater"

/*
 *  permutation flowshop neighborhoods
 *
 */

/* Generic */
#define NEIGHBORHOOD_INSERT "insert"
#define NEIGHBORHOOD_BACK_INSERT "binsert"
#define NEIGHBORHOOD_FORW_INSERT "finsert"
#define NEIGHBORHOOD_TWO_INSERT "tinsert"
#define NEIGHBORHOOD_TRANSPOSE "transpose"
#define NEIGHBORHOOD_XTRANSPOSE "xtranspose"
#define NEIGHBORHOOD_EXCHANGE "exchange"
#define NEIGHBORHOOD_ADAPTIVE_INSERT "adpinsert"
#define NEIGHBORHOOD_RANDCONHE "rch"

/* Weighted Tardiness*/
#define NEIGHBORHOOD_ATX_EXCHANGE "atxexchange"
#define NEIGHBORHOOD_HATX_EXCHANGE "hatxexchange"
#define NEIGHBORHOOD_EATX_EXCHANGE "eatxexchange"
#define NEIGHBORHOOD_OPT_EXCHANGE "oexchange"
#define NEIGHBORHOOD_TAx_INSERT "txinsert"
#define NEIGHBORHOOD_ATAx_INSERT "atxinsert"
#define NEIGHBORHOOD_OPT_INSERT "oinsert"
#define NEIGHBORHOOD_HATAx_INSERT "hatxinsert"
#define NEIGHBORHOOD_NATAx_INSERT "natxinsert"
#define NEIGHBORHOOD_NATA2x_INSERT "natx2insert"
#define NEIGHBORHOOD_THATAx_INSERT "thatxinsert"
#define NEIGHBORHOOD_FATAx_INSERT "fatxinsert"
#define NEIGHBORHOOD_PATAx_INSERT "patxinsert"
#define NEIGHBORHOOD_SATAx_INSERT "satxinsert"
#define NEIGHBORHOOD_EATAx_INSERT "eatxinsert"
#define NEIGHBORHOOD_TATAx_INSERT "tatxinsert"

/* Total Completion Time*/
#define NEIGHBORHOOD_NATA_TCT_INSERT "ntctinsert"
#define NEIGHBORHOOD_RZ_TCT_INSERT "nrztctinsert"

/* Total Tardiness*/
#define NEIGHBORHOOD_NATA_TT_INSERT "nttinsert"

/* Makespan */
#define NEIGHBORHOOD_TA_INSERT "tainsert"
#define NEIGHBORHOOD_FSTA_INSERT "fstainsert"
#define NEIGHBORHOOD_CSTA_INSERT "cstainsert"

/* No idle makespan*/
#define NEIGHBORHOOD_NITA_INSERT "ntainsert"


/* Sequence Dependent Setup times makespan*/
#define NEIGHBORHOOD_SDSTTA_INSERT "sdsttainsert"

/*
 * END Neighborhoods
 */


/* permutation flowshop solution perturbations */
#define PERTURBATION_RANDOM_MOVE "rndmv"
#define PERTURBATION_VNRANDOM_MOVE "vrndmv"
#define PERTURBATION_NOPER "noper"
#define PERTURBATION_RND "randpert"
#define PERTURBATION_NRZ "nrzper"
#define PERTURBATION_TMIIG "tmiigper"
#define PERTURBATION_SOA "igper"
#define PERTURBATION_SOA_LEGACY "soaper"
#define PERTURBATION_TEST "testper"
#define PERTURBATION_IGLS "igls"
#define PERTURBATION_RSLS "rsls"
#define PERTURBATION_RSffLS "rsffls"
#define PERTURBATION_RS "rsper"
#define PERTURBATION_RSFF "rsff"
#define PERTURBATION_IGIO "igio"
#define PERTURBATION_RSIO "rsio"
#define PERTURBATION_CP3 "cp3"

/* acceptance criteria*/
#define ACCEPTANCE_PROB "prob"
#define ACCEPTANCE_METRO "metropolis"
#define ACCEPTANCE_RS "rsacc"
#define ACCEPTANCE_PMETRO "pmetro"
#define ACCEPTANCE_TMIIG "tmiigacc"
#define ACCEPTANCE_IMPROVE_PLATEAU "implat"
#define ACCEPTANCE_TEST "testacc"
#define ACCEPTANCE_SOA "soaacc"
#define ACCEPTANCE_ALWAYS "always"
#define ACCEPTANCE_INTENSIFY "intensify"
#define ACCEPTANCE_DIVERSIFY "diversify"
#define ACCEPTANCE_IMPROVE "improve"
#define ACCEPTANCE_SA_METRO "sa_metropolis"
#define ACCEPTANCE_SA "saacc"
#define ACCEPTANCE_KAR "karacc"
char* problem_type;

emili::pfsp::PermutationFlowShop* instantiateProblem(char* t, PfspInstance i)
{
    emili::pfsp::PermutationFlowShop* prob;
    if(strcmp(t,PROBLEM_PFS_WT)==0)
    {
#ifdef ENABLE_SSE
        prs::printTab("SSE run path enabled");
#endif
        prs::printTab("Permutation Flow Shop Weighted Tardiness");
        prob = new emili::pfsp::PFSP_WT(i);
    }else if(strcmp(t,PROBLEM_PFS_E)==0)
    {
#ifdef ENABLE_SSE
        prs::printTab("SSE run path enabled");
#endif
        prs::printTab("Permutation Flow Shop Earliness");
        prob = new emili::pfsp::PFSP_E(i);
    }else if(strcmp(t,PROBLEM_PFS_WE)==0)
    {
#ifdef ENABLE_SSE
        prs::printTab("SSE run path enabled");
#endif
        prs::printTab("Permutation Flow Shop Weighted Earliness");
        prob = new emili::pfsp::PFSP_WE(i);
    }else if(strcmp(t,PROBLEM_PFS_T)==0)
    {
#ifdef ENABLE_SSE
        prs::printTab("SSE run path enabled");
#endif

        prs::printTab("Permutation Flow Shop Tardiness");
        prob = new emili::pfsp::PFSP_T(i);
    }else if(strcmp(t,PROBLEM_PFS_MS)==0)
    {
#ifdef ENABLE_SSE
        prs::printTab("SSE run path enabled");
#endif
        prs::printTab("Permutation Flow Shop Make Span");
        prob = new emili::pfsp::PFSP_MS(i);
    }
    else if(strcmp(t,PROBLEM_PFS_TCT)==0)
    {
#ifdef ENABLE_SSE
        prs::printTab("SSE run path enabled");
#endif
        prs::printTab("Permutation Flow Shop TCT");
        prob = new emili::pfsp::PFSP_TCT(i);
    }
    else if(strcmp(t,PROBLEM_PFS_WCT)==0)
    {
#ifdef ENABLE_SSE
        prs::printTab("SSE run path enabled");
#endif

        prs::printTab("Permutation Flow Shop WCT");
        prob = new emili::pfsp::PFSP_WCT(i);
    }
    else if(strcmp(t,PROBLEM_NWPFS_MS)==0)
    {
        prs::printTab("No Wait Permutation Flow Shop Make Span");
        prob = new emili::pfsp::NWPFSP_MS(i);
    }
    else if(strcmp(t,PROBLEM_NWPFS_WT)==0)
    {
        prs::printTab("No Wait Permutation Flow Shop Weighted Tardiness");
        prob = new emili::pfsp::NWPFSP_WT(i);
    }
    else if(strcmp(t,PROBLEM_NWPFS_WE)==0)
    {
        prs::printTab("No Wait Permutation Flow Shop Weighted Earliness");
        prob = new emili::pfsp::NWPFSP_WE(i);
    }
    else if(strcmp(t,PROBLEM_NWPFS_T)==0)
    {
        prs::printTab("No Wait Permutation Flow Shop Tardiness");
        prob = new emili::pfsp::NWPFSP_T(i);
    }
    else if(strcmp(t,PROBLEM_NWPFS_E)==0)
    {
        prs::printTab("No Wait Permutation Flow Shop Earliness");
        prob = new emili::pfsp::NWPFSP_E(i);
    }
    else if(strcmp(t,PROBLEM_NWPFS_TCT)==0)
    {
        prs::printTab("No Wait Permutation Flow Shop TCT");
        prob = new emili::pfsp::NWPFSP_TCT(i);
    }
    else if(strcmp(t,PROBLEM_NWPFS_WCT)==0)
    {
        prs::printTab("No Wait Permutation Flow Shop WCT");
        prob = new emili::pfsp::NWPFSP_WCT(i);
    }
    else if(strcmp(t,PROBLEM_NIPFS_MS)==0)
    {
        prs::printTab("No Idle Permutation Flow Shop Make Span" );
        prob = new emili::pfsp::NI_A_PFSP_MS(i);
    }
    else if(strcmp(t,PROBLEM_NIPFS_WT)==0)
    {
        prs::printTab("No Idle Permutation Flow Shop Weighted Tardiness" );
        prob = new emili::pfsp::NIPFSP_WT(i);
    }
    else if(strcmp(t,PROBLEM_NIPFS_WE)==0)
    {
        prs::printTab("No Idle Permutation Flow Shop Weighted Earliness" );
        prob = new emili::pfsp::NIPFSP_WE(i);
    }
    else if(strcmp(t,PROBLEM_NIPFS_T)==0)
    {
        prs::printTab("No Idle Permutation Flow Shop Tardiness" );
        prob = new emili::pfsp::NIPFSP_T(i);
    }
    else if(strcmp(t,PROBLEM_NIPFS_E)==0)
    {
        prs::printTab("No Idle Permutation Flow Shop Earliness" );
        prob = new emili::pfsp::NIPFSP_E(i);
    }
    else if(strcmp(t,PROBLEM_NIPFS_TCT)==0)
    {
        prs::printTab("No Idle Permutation Flow Shop total completion time" );
        prob = new emili::pfsp::NIPFSP_TCT(i);
    }
    else if(strcmp(t,PROBLEM_NIPFS_WCT)==0)
    {
        prs::printTab("No Idle Permutation Flow Shop Weighted completion time" );
        prob = new emili::pfsp::NIPFSP_WCT(i);
    }
    else if(strcmp(t,PROBLEM_SDSTPFS_MS)==0)
    {
        prs::printTab("Sequence dependent setup times Make Span");
        prob = new emili::pfsp::SDSTFSP_MS(i);
    }
    else if(strcmp(t,PROBLEM_SDSTPFS_WT)==0)
    {
        prs::printTab("Sequence dependent setup times weighted tardiness");
        prob = new emili::pfsp::SDSTFSP_WT(i);
    }
    else if(strcmp(t,PROBLEM_SDSTPFS_WE)==0)
    {
        prs::printTab("Sequence dependent setup times weighted earliness");
        prob = new emili::pfsp::SDSTFSP_WE(i);
    }
    else if(strcmp(t,PROBLEM_SDSTPFS_T)==0)
    {
        prs::printTab("Sequence dependent setup times Tardiness");
        prob = new emili::pfsp::SDSTFSP_T(i);
    }
    else if(strcmp(t,PROBLEM_SDSTPFS_E)==0)
    {
        prs::printTab("Sequence dependent setup times Earliness");
        prob = new emili::pfsp::SDSTFSP_E(i);
    }
    else if(strcmp(t,PROBLEM_SDSTPFS_TCT)==0)
    {
        prs::printTab("Sequence dependent setup times Total completion time");
        prob = new emili::pfsp::SDSTFSP_TCT(i);
    }
    else if(strcmp(t,PROBLEM_SDSTPFS_WCT)==0)
    {
        prs::printTab("Sequence dependent setup times weighted completion time");
        prob = new emili::pfsp::SDSTFSP_WCT(i);
    }
    else
    {
        std::cerr<< "'" << t << "' -> ERROR a problem was expected! " << std::endl;
        prs::info();
        exit(-1);
    }
    return prob;
}


/*void prs::emili_header()
{
    std::cout << "\t ______ __  __ _____ _      _____ " << std::endl;
    std::cout << "\t|  ____|  \\/  |_   _| |    |_   _|" << std::endl;
    std::cout << "\t| |__  | \\  / | | | | |      | |  " << std::endl;
    std::cout << "\t|  __| | |\\/| | | | | |      | |  " << std::endl;
    std::cout << "\t| |____| |  | |_| |_| |____ _| |_ " << std::endl;
    std::cout << "\t|______|_|  |_|_____|______|_____|" << std::endl;
    std::cout << std::endl;
}*/

std::string prs::ParamsParser::info()
{
    std::ostringstream oss;
    oss << "Usage:\n\n";
    oss << "EMILI INSTANCE_FILE_PATH PFS_PROBLEM <LOCAL_SEARCH | ITERATED_LOCAL_SEARCH | TABU_SEARCH | VND_SEARCH> [rnds seed]\n\n";
    oss << "PROBLEM               = "<<PROBLEM_PFS_WT<< " " <<PROBLEM_PFS_WE<< " " <<PROBLEM_PFS_TCT<< " " <<PROBLEM_PFS_MS
              << " " <<PROBLEM_PFS_WCT<< " " <<PROBLEM_PFS_T<< " " <<PROBLEM_PFS_E<< " "<<PROBLEM_NWPFS_WT<< " " <<PROBLEM_NWPFS_WE
              << " " <<PROBLEM_NWPFS_TCT<< " " <<PROBLEM_NWPFS_MS<< " " <<PROBLEM_NWPFS_T<< " " <<PROBLEM_NWPFS_E
              << PROBLEM_NIPFS_MS <<" "<<PROBLEM_NIPFS_WT<< " " <<PROBLEM_NIPFS_WE<< " " <<PROBLEM_NIPFS_TCT<< " " <<PROBLEM_NIPFS_MS<< " "
              << " " <<PROBLEM_NIPFS_T<< " " <<PROBLEM_NIPFS_E <<"\n";
    oss << "LOCAL_SEARCH          = SEARCH_TYPE INITIAL_SOLUTION TERMINATION NEIGHBORHOOD" <<"\n";
    oss << "ITERATED_LOCAL_SEARCH = ils LOCAL_SEARCH TERMINATION PERTURBATION ACCEPTANCE -it seconds" << "\n";
    oss << "TABU_SEARCH           = tabu < first | best > INITIAL_SOLUTION TERMINATION NEIGHBORHOOD TABU_MEMORY " << "\n";
    oss << "VND_SEARCH            = vnd < first | best > INITIAL_SOLUTION TERMINATION NEIGHBORHOOD1 NEIGHBORHOOD2 ... NEIGHBORHOODn" << "\n";
    oss << "GVNS_SEARCH           = gvns INITIAL_SOLUTION PERTURBATION1 PERTURBATION2 -it seconds" << "\n";
    oss << "SEARCH_TYPE           = first | best | tabu | vnd | ils" << "\n";
    oss << "INITIAL_SOLUTION      = random | slack | nwslack | lit | rz | nrz | nrz2 | lr size(int)| nlr size(int) | mneh" <<"\n";
    oss << "TERMINATION           = true | time float | locmin | soater | iteration int | maxstep int" << "\n";
    oss << "NEIGHBORHOOD          = transpose | exchange | insert | binsert | finsert | tinsert | "<< NEIGHBORHOOD_TA_INSERT << " | " << NEIGHBORHOOD_NITA_INSERT<<"\n";
    oss << "PERTURBATION           = igper int | testper | rndmv NEIGHBORHOOD #moves(int) | noper (int) | nrzper (int) | tmiigper (int) (int) | igls (int) LOCAL_SEARCH | rsls (int) LOCAL_SEARCH" << "\n";
    oss << "ACCEPTANCE            = soaacc float | testacc #swaps(int) | metropolis start_temperature(float) | always (intensify | diversify) | improve | sa_metropolis start_temp end_temp ratio | pmetro start_temp end_temp ratio frequence(int) | saacc start_temp end_temp ratio frequence(int) alpha ]0,1] | tmiigacc start_temperature(float) | implat number_of_non_improving_steps_accepted plateau_threshold" << std::endl;
    oss << "TABU_MEMORY           = move size(int) | hash size(int) | solution size(int) | tsabm size(int)" << "\n";
   // std::cout << " syntax->EMILI instancefile search_type intial_solution termination neighborhood" << std::endl;
    return oss.str();
}

emili::Neighborhood* ne = nullptr;
emili::InitialSolution* in= nullptr;
emili::Termination* te= nullptr;
emili::TabuMemory* tmem= nullptr;
emili::Termination* ilt= nullptr;

std::vector< emili::Neighborhood*> nes;

emili::LocalSearch* prs::ParamsParser::eparams(prs::TokenManager& tm)
{
    prs::incrementTabLevel();
    emili::LocalSearch* ls;
    if(tm.checkToken(ILS))
    {
        printTab("ILS");
        ls = ils(tm);
    }  
    else if(tm.checkToken(GVNS_ILS))
    {
        printTab("GVNS...");
        ls = gvns(tm);
    }
    /*else if(tm.checkToken(ALBERTOSA))
    {
        / **
spostare          * /
            emili::InitialSolution* initsol    = init(tm);
    emili::Neighborhood*    nei        = neigh(tm, true);
    SAInitTemp*      inittemp   = INITTEMP(tm, initsol, nei, instance);
    SAAcceptance*    acceptance = ACCEPTANCE(tm, inittemp, nei, instance);
    SACooling*       cooling    = COOL(tm, inittemp, nei, instance);
    SATempRestart*   temprestart = TEMPRESTART(tm, inittemp, nei);
    cooling->setTempRestart(temprestart);
    SATermination*     term       = TERMINATION(tm, inittemp, nei); // termin(tm);
    SATempLength*    templ      = TEMPLENGTH(tm, nei, instance);
    cooling->setTempLength(templ);
    SAExploration* explo = EXPLORATION(tm, nei, acceptance, cooling, term);

    return new SimulatedAnnealing(initsol,
                                  inittemp,
                                  acceptance,
                                  cooling,
                                  temprestart,
                                  term,
                                  templ,
                                  explo,
                                  nei,
                                  NULL);
    / **
     * 
     * /
    }*/
    else
    {             
        ls = search(tm);
    }

    prs::decrementTabLevel();
    return ls;
}



emili::LocalSearch* prs::ParamsParser::search(prs::TokenManager& tm)
{
    prs::incrementTabLevel();
    emili::LocalSearch* ls;
    if(tm.checkToken(ILS))
    {
        printTab("ILS ");
        ls = ils(tm);

    }else if(tm.checkToken(TABU))
    {
        printTab("TABU SEARCH");
        ls = tparams(tm);
    }
    else if(tm.checkToken(FIRST))
    {
        printTab("FIRST IMPROVEMENT");
        params(tm);
        ls =  new emili::FirstImprovementSearch(*in,*te,*ne);
    }
    else if(tm.checkToken(BEST))
    {
        printTab("BEST IMPROVEMENT");
        params(tm);
        ls =  new emili::BestImprovementSearch(*in,*te,*ne);
    }
    else if(tm.checkToken(TB_FIRST))
    {
        printTab("FIRST IMPROVEMENT");
        params(tm);
        emili::Problem* p = (emili::Problem*)instantiateProblem(tm.nextToken(),instance->getInstance());
        ls =  new emili::TieBrakingFirstImprovementSearch(*in,*te,*ne,*p);
    }
    else if(tm.checkToken(TB_BEST))
    {
        printTab("BEST IMPROVEMENT");
        params(tm);
        emili::Problem* p = (emili::Problem*)instantiateProblem(tm.nextToken(),instance->getInstance());
        ls =  new emili::TieBrakingBestImprovementSearch(*in,*te,*ne,*p);
    }
    else if(tm.checkToken(CH6_LS))
    {
        printTab("CH6");
        ls = ch6_params(tm);
    }
    else if(tm.checkToken(VND))
    {
        printTab("VND SEARCH");
        ls = vparams(tm);
    }
    else if(tm.checkToken(TEST_INIT))
    {       
        emili::InitialSolution* ini = init(tm);        
        /*ls = new emili::EmptyLocalSearch(*ini);*/
        clock_t time = clock();
        emili::Solution* s = ini->generateSolution();
        double time_elapsed = (double)(clock()-time)/CLOCKS_PER_SEC;
        std::cout << "time : " << time_elapsed << std::endl;
        std::cout << "iteration counter : " << emili::iteration_counter()<< std::endl;
        std::cout << "Objective function value: " << s->getSolutionValue() << std::endl;
        std::cout << "Found solution: ";
        std::cout << s->getSolutionRepresentation() << std::endl;
        std::cout << std::endl;
        std::cerr << s-> getSolutionValue() << std::endl;
        exit(123);
    }
    else if(tm.checkToken(EMPTY_LOCAL))
    {
        printTab("NO LOCAL SEARCH");
        emili::InitialSolution* ini = init(tm);
        ls = new emili::EmptyLocalSearch(*ini);
    }
    else
    {
        std::cerr<< "'" << tm.peek() << "' -> ERROR a search definition was expected! " << std::endl;
        std::cout << info() << std::endl;
        exit(-1);
    }
    prs::decrementTabLevel();
    return ls;

}

emili::LocalSearch* prs::ParamsParser::ils(prs::TokenManager& tm)
{

    emili::LocalSearch* ls = search(tm);
    //ils_time = ilstime();
    emili::Termination* pft = term(tm);
    //emili::pfsp::PfspRandomSwapPertub* prsp = new emili::pfsp::PfspRandomSwapPertub(istance);
    int rpc = 5;
    emili::Perturbation* prsp = per(tm);
    //emili::AcceptanceCriteria* tac = new emili::pfsp::PfspTestAcceptance(istance);
    //emili::AcceptanceCriteria* tac = new emili::MetropolisAcceptance(1);
    emili::Acceptance* tac = acc(tm);//new emili::pfsp::SOAacceptance(1.2f);
    emili::LocalSearch* iils = new emili::IteratedLocalSearch(*ls,*pft,*prsp,*tac);   
    return iils;
}

emili::Perturbation* prs::ParamsParser::per(prs::TokenManager& tm)
{
    prs::incrementTabLevel();
    std::ostringstream oss;
    emili::Perturbation* per;
    if(tm.checkToken(PERTURBATION_SOA) || tm.checkToken(PERTURBATION_SOA_LEGACY))
    {
        int nj = instance->getNjobs();
        int n = tm.getInteger();
        n = n<nj?n:nj-1;
        oss << "NEH destruct/construct perturbation which use objective function. number of job erased: "<<n;
        printTab(oss.str().c_str());
        per = new emili::pfsp::IGPerturbation(n,*instance);
    }else if(tm.checkToken(PERTURBATION_RS))
    {
        int nj = instance->getNjobs();
        int n = tm.getInteger();
        n = n<nj?n:nj-1;
        oss << "NEH destruct/construct perturbation. number of job erased: "<<n;
        printTab(oss.str().c_str());
        per = new emili::pfsp::RSPerturbation(n,*instance);
    }
    else if(tm.checkToken(PERTURBATION_RSFF))
        {
        int nj = instance->getNjobs();
        int n = tm.getInteger();
        n = n<nj?n:nj-1;
            oss << "NEH destruct/construct perturbation with tbff tie breaking. number of job erased: "<<n;
            printTab(oss.str().c_str());
            per = new emili::pfsp::RSffPerturbation(n,*instance);
        }
    else if(tm.checkToken(PERTURBATION_IGLS))
    {
        int nj = instance->getNjobs()-2;
        int k = tm.getInteger();
        int n = k<nj?k:nj-1;
        oss.str(""); oss  << "IG perturbation with local search applied on the partial solution. d = "<<n;
        printTab(oss.str().c_str());
        if(n > 0)
        {
            PfspInstance pfs = this->instance->getInstance();
            pfs.setNbJob(pfs.getNbJob()-n);
            emili::pfsp::PermutationFlowShop * pfse = instantiateProblem(problem_type,pfs);
            emili::pfsp::PermutationFlowShop* is = this->instance;
            this->instance = pfse;
            emili::LocalSearch* ll = search(tm);
            this->instance = is;
            istances.push_back(pfse);
            per = new emili::pfsp::IgLsPerturbation(n,*instance,ll);
        } else {
             per = new emili::pfsp::IGPerturbation(1,*instance);
        }
    }
    else if(tm.checkToken(PERTURBATION_RSLS))
    {
        int nj = instance->getNjobs()-2;
        int n = tm.getInteger();
        n = n<nj?n:nj-1;
        oss.str(""); oss  << "IG perturbation with local search applied on the partial solution. d = "<<n;
        printTab(oss.str().c_str());
        if(n > 0)
        {
            PfspInstance pfs = this->instance->getInstance();
            pfs.setNbJob(pfs.getNbJob()-n);
            emili::pfsp::PermutationFlowShop * pfse = instantiateProblem(problem_type,pfs);
            emili::pfsp::PermutationFlowShop* is = this->instance;
            this->instance = pfse;
            emili::LocalSearch* ll = search(tm);
            this->instance = is;
            istances.push_back(pfse);
            per = new emili::pfsp::RSLSPerturbation(n,*instance,ll);
        } else {
            per = new emili::pfsp::RSPerturbation(n,*instance);
        }
    }
    else if(tm.checkToken(PERTURBATION_RSffLS))
    {
        int nj = instance->getNjobs()-2;
        int n = tm.getInteger();
        n = n<nj?n:nj-1;
        oss.str(""); oss  << "IG perturbation with tbff tie breaking and local search applied on the partial solution. d = "<<n;
        printTab(oss.str().c_str());
        if(n > 0)
        {
            PfspInstance pfs = this->instance->getInstance();
            pfs.setNbJob(pfs.getNbJob()-n);
            emili::pfsp::PermutationFlowShop * pfse = instantiateProblem(problem_type,pfs);
            emili::pfsp::PermutationFlowShop* is = this->instance;
            this->instance = pfse;
            emili::LocalSearch* ll = search(tm);
            this->instance = is;
            istances.push_back(pfse);
            per = new emili::pfsp::RSffLSPerturbation(n,*instance,ll);
        } else {
            per = new emili::pfsp::RSffPerturbation(n,*instance);
        }
    }
    else if(tm.checkToken(PERTURBATION_TEST))
    {
        oss.str(""); oss<< "Random swap test perturbation.";
        printTab(oss.str().c_str());
        per = new emili::pfsp::PfspRandomSwapPertub(*instance);
    }else if(tm.checkToken(PERTURBATION_RANDOM_MOVE))
    {
        printTab("Random move perturbation.");
        emili::Neighborhood* n = neigh(tm,true);
        int num = tm.getInteger();
        prs::incrementTabLevel();
        oss.str(""); oss  << "number of moves per perturbation step " << num;
        printTab(oss.str().c_str());
        prs::decrementTabLevel();
        per = new emili::RandomMovePerturbation(*n,num);
    }
    else if(tm.checkToken(PERTURBATION_NOPER))
    {
        printTab("No PERTURBATION.");
        per = new emili::NoPerturbation();
    }
    else if(tm.checkToken(PERTURBATION_NRZ))
    {
        int nj = instance->getNjobs();
        int n = tm.getInteger();
        n = n<nj?n:nj-1;
        oss.str(""); oss  << "neh rz destruct/construct PERTURBATION. number of job erased: "<<n;
        printTab(oss.str().c_str());
        per = new emili::pfsp::NRZPerturbation(n,*instance);
    }else if(tm.checkToken(PERTURBATION_VNRANDOM_MOVE))
    {
        printTab("Random move perturbation." );
        prs::incrementTabLevel();
        int num = tm.getInteger();
        oss.str(""); oss  << "number of moves per perturbation step " << num << ".\n\t";
        printTab(oss.str().c_str());
        int iter = tm.getInteger();
        oss.str(""); oss  << "number of iteration before changing the neighborhood " << iter << ".\n\t";
        printTab(oss.str().c_str());
        nes.clear();
        prs::decrementTabLevel();
        neighs(tm);
        per = new emili::VNRandomMovePerturbation(nes,num,iter);
    }
    else if(tm.checkToken(PERTURBATION_TMIIG))
    {
        int nj = instance->getNjobs();
        int n = tm.getInteger();
        n = n<nj?n:nj-1;
        int tsize = tm.getInteger();
        oss.str(""); oss  << "TMIIG PERTURBATION. Number of job erased " << n << ". tabu list size " << tsize <<".\n\t";
        printTab(oss.str().c_str());
        per = new emili::pfsp::TMIIGPerturbation(n,*instance,tsize);
    }
    else if(tm.checkToken(PERTURBATION_IGIO))
    {
        int nj = instance->getNjobs();
        int n = tm.getInteger();
        n = n<nj?n:nj-1;

        oss.str(""); oss  << "IG perturbation that inserts first the removed job with max sum of processing times. d= " << n <<".\n\t";
        printTab(oss.str().c_str());
        per = new emili::pfsp::IGIOPerturbation(n,*instance);
    }
    else if(tm.checkToken(PERTURBATION_RSIO))
    {
        int nj = instance->getNjobs();
        int n = tm.getInteger();
        n = n<nj?n:nj-1;

        oss.str(""); oss  << "IG perturbation that inserts first the removed job with max sum of processing times using taillard acceleration. d= " << n <<".\n\t";
        printTab(oss.str().c_str());
        per = new emili::pfsp::RSIOPerturbation(n,*instance);
    }
    else if(tm.checkToken(PERTURBATION_CP3))
    {
        int d = tm.getInteger();
        int omega = tm.getInteger();
        float pc = tm.getDecimal();
        oss.str(""); oss  << "Compound perturbation :  d= " << d << ", omega= " << omega << ",pc= "<< pc;
        printTab(oss.str().c_str());
        per = new emili::pfsp::CompoundPerturbation(*instance,omega,d,pc);
    }
    else
    {
        std::cerr<< "'" << tm.peek() << "' -> ERROR a perturbation criteria specification was expected! " << std::endl;
        std::cout << info() << std::endl;
        exit(-1);
    }
    prs::decrementTabLevel();
    return per;
}

emili::Acceptance* prs::ParamsParser::acc(prs::TokenManager& tm)
{
    prs::incrementTabLevel();
    emili::Acceptance* acc;
    std::ostringstream oss;
    if(tm.checkToken(ACCEPTANCE_SOA))
    {
        float n = tm.getDecimal();
        oss.str(""); oss  << "soa metropolis like acceptance. temperature : "<<n;
        printTab(oss.str().c_str());
        acc = new  emili::pfsp::SOAacceptance(n);
    }
    else if(tm.checkToken(ACCEPTANCE_TEST))
    {
        int n = tm.getInteger();
        oss.str(""); oss  << "Probabilistic Acceptance. improving solution accepted"<<n<<" % of the time";
        printTab(oss.str().c_str());
        acc = new  emili::pfsp::PfspTestAcceptance(*instance,n);
    }
    else  if(tm.checkToken(ACCEPTANCE_METRO))
    {
        float n = tm.getDecimal();
        oss.str(""); oss  << "metropolis acceptance. temperature : "<<n;
        printTab(oss.str().c_str());
        acc = new  emili::MetropolisAcceptance(n);
    }
    else  if(tm.checkToken(ACCEPTANCE_RS))
    {
        float n = tm.getDecimal();
        const std::vector < std::vector < long int > >& pm = instance->getProcessingTimesMatrix();
        int nj = instance->getNjobs();
        int nm = instance->getNmachines();

        float temp = 0;
        for(int i = 1; i<= nj; i++ )
        {
            for(int j=1; j<=nm; j++)
            {
                temp += pm[i][j];
            }
        }

        temp = n*(temp/(nj*nm))/10;

        oss.str(""); oss  << "metropolis like Ruiz Stuetzle 2006 acceptance. temperature : "<<temp;
        printTab(oss.str().c_str());
        acc = new  emili::MetropolisAcceptance(temp);
    }
    else  if(tm.checkToken(ACCEPTANCE_KAR))
    {
        float n = tm.getDecimal();
        int nj = instance->getNjobs();
        int lb = instance->getInstance().computeMSLB();
        std::vector<long>& dd = instance->getDueDates();
        float temp = 0;
        for(int i = 1; i<= nj; i++ )
        {
            temp += lb - dd[i];
        }

        temp = n*(temp/nj)/10;

        oss.str(""); oss  << "metropolis like Kar2016 acceptance. temperature : "<<temp;
        printTab(oss.str().c_str());
        acc = new  emili::MetropolisAcceptance(temp);
    }
    else  if(tm.checkToken(ACCEPTANCE_ALWAYS))
    {


        emili::accept_candidates accc;
        std::string t1;
        if(tm.checkToken(ACCEPTANCE_INTENSIFY))
        {
            accc = emili::ACC_INTENSIFICATION;
            t1 = ACCEPTANCE_INTENSIFY;
        }
        else if(tm.checkToken(ACCEPTANCE_DIVERSIFY))
        {
            t1 = ACCEPTANCE_DIVERSIFY;
            accc = emili::ACC_DIVERSIFICATION;
        }
        else
        {
            std::cerr<< "'" << *tm << "' -> ERROR " << ACCEPTANCE_INTENSIFY << " or " << ACCEPTANCE_DIVERSIFY <<" was expected! " << std::endl;
            std::cout << info() << std::endl;
        exit(-1);
        }
        oss.str(""); oss  << "Acceptance always "<< t1;
        printTab(oss.str().c_str());
        acc = new  emili::AlwaysAccept(accc);
    }
    else  if(tm.checkToken(ACCEPTANCE_IMPROVE))
    {

        printTab( "improve acceptance");

        acc = new  emili::ImproveAccept();
    }
    else  if(tm.checkToken(ACCEPTANCE_SA_METRO))
    {
        float start =tm.getDecimal();
        float end =tm.getDecimal();
        float ratio =tm.getDecimal();
        oss.str(""); oss  << "metropolis acceptance. start ,end , ratio : "<< start << ", "<< end << "," << ratio;
        printTab(oss.str().c_str());
        acc = new  emili::Metropolis(start,end,ratio);
    }
    else  if(tm.checkToken(ACCEPTANCE_PMETRO))
    {
        float start =tm.getDecimal();
        float end =tm.getDecimal();
        float ratio =tm.getDecimal();
        int iterations = tm.getInteger();
        oss.str(""); oss  << "metropolis acceptance. start ,end , ratio, frequence : "<< start << ", "<< end << "," << ratio <<","<< iterations;
        printTab(oss.str().c_str());
        acc = new  emili::Metropolis(start,end,ratio,iterations);
    }
    else if(tm.checkToken(ACCEPTANCE_SA))
    {
        float start =tm.getDecimal();
        float end =tm.getDecimal();
        float ratio =tm.getDecimal();
        int iterations = tm.getInteger();
        float alpha =tm.getDecimal();
        oss.str(""); oss  << "metropolis acceptance. start ,end , ratio, frequence, alpha : "<< start << ", "<< end << "," << ratio <<","<< iterations << "," << alpha;
        printTab(oss.str().c_str());
        acc = new  emili::Metropolis(start,end,ratio,iterations,alpha);
    }
    else if(tm.checkToken(ACCEPTANCE_TMIIG))
    {
        float t0 =tm.getDecimal();
        float t=0.0f;
        int nj = instance->getNjobs();
        int nm = instance->getNmachines();
        const std::vector< std::vector < long > >& p = instance->getProcessingTimesMatrix();

        for(int i=1; i<=nj ; i++)
        {
            for(int j=1;j<=nm; j++)
            {
                t += p[i][j];
            }
        }
        t = (t*t0)/(10.0f*nj*nm);
        oss.str(""); oss  << "TMIIG metropolis like acceptance criterion. temperature " << t;
        printTab(oss.str().c_str());
        acc = new  emili::MetropolisAcceptance(t);
    }else if(tm.checkToken(ACCEPTANCE_IMPROVE_PLATEAU))
    {
        int plateau_steps = tm.getInteger();
        int threshold = tm.getInteger();
        oss.str(""); oss  << "Accept a diversification solution if it improves on the intensification otherwise it will accept "<< plateau_steps << " non improving steps once it reaches the threshold of " << threshold;
        printTab(oss.str().c_str());
        acc = new  emili::AcceptPlateau(plateau_steps,threshold);
    }
    else
    {
        std::cerr<< "'" << *tm << "' -> ERROR an acceptance criteria specification was expected! " << std::endl;
        std::cout << info() << std::endl;
        exit(-1);
    }
    prs::decrementTabLevel();
    return acc;
}


emili::LocalSearch* prs::ParamsParser::gvns(prs::TokenManager& tm)
{
  return nullptr;
}

emili::BestTabuSearch* prs::ParamsParser::tparams(prs::TokenManager& tm)
{
    emili::pfsp::PfspNeighborhood* nei = dynamic_cast<emili::pfsp::PfspNeighborhood*>(ne);
    if(tm.checkToken(BEST))
    {
        params(tm);
        tmem = tmemory(nei,tm);
        return new emili::BestTabuSearch(*in,*te,*ne,*tmem);
    }
    else if(tm.checkToken(FIRST))
    {
        params(tm);
        tmem = tmemory(nei,tm);
        return new emili::FirstTabuSearch(*in,*te,*ne,*tmem);
    }
    else
    {
        std::cerr<< "'" << *tm << "' -> ERROR a pivotal rule (best or first) for the tabu search was expected! \n" << std::endl;
        std::cout << info() << std::endl;
        exit(-1);
    }
}

emili::TabuMemory* prs::ParamsParser::tmemory(emili::pfsp::PfspNeighborhood* n,prs::TokenManager& tm)
{
    prs::incrementTabLevel();

    std::ostringstream oss;
    emili::TabuMemory* tmem;


    printTab(oss.str().c_str());

    if(tm.checkToken(TABU_MEMORY_MOVES))
    {
        oss.str(""); oss << "USING MOVES\n\t";
        printTab(oss.str().c_str());
        int ts = tm.getInteger();
        oss << "Tabu tenure size " << ts;
        tmem = new  emili::pfsp::PfspMovesMemory(ts , n);
    }
    else if(tm.checkToken(TABU_MEMORY_HASHES))
    {
        oss.str(""); oss << "USING HASHES\n\t";
        printTab(oss.str().c_str());
        int ts = tm.getInteger();
        oss << "Tabu tenure size " << ts;
        tmem = new  emili::pfsp::PfspTabuHashMemory(ts);
    }
    else if(tm.checkToken(TABU_MEMORY_SOLUTIONS))
    {
        oss.str(""); oss << "USING FULL SOLUtiON\n\t";
        printTab(oss.str().c_str());
        int ts = tm.getInteger();
        oss << "Tabu tenure size " << ts;
        tmem = new  emili::pfsp::PfspFullSolutionMemory(ts);
    }
    else if(tm.checkToken(TABU_MEMORY_TSAB))
    {
        oss.str(""); oss << "USING TSAB\n\t";
        printTab(oss.str().c_str());
        int ts = tm.getInteger();
        oss << "Tabu tenure size " << ts;
        tmem = new  emili::pfsp::TSABMemory(ts , n);
    }
    else if(tm.checkToken(TABU_MEMORY_TSAB_TEST))
    {
        oss.str(""); oss << "USING TSAB\n\t";
        printTab(oss.str().c_str());
        int ts = tm.getInteger();
        oss << "Tabu tenure size " << ts;
        tmem = new  emili::pfsp::TSABtestMemory(ts , n);
    }
    else if(tm.checkToken(TABU_MEMORY_VALUE))
    {
        oss.str(""); oss << "USING VALUE\n\t" ;
        printTab(oss.str().c_str());
        int ts = tm.getInteger();
        oss << "Tabu tenure size " << ts;
        tmem = new  emili::pfsp::PfspTabuValueMemory(ts);
    }
    else
    {
        std::cerr<< "'" << *tm << "' -> ERROR a memory specification for the tabu search was expected! " << std::endl;
        std::cout << info() << std::endl;
        exit(-1);
    }
    prs::decrementTabLevel();
    return tmem;
}

void prs::ParamsParser::params(prs::TokenManager& tm)
{
    in = init(tm);
    te = term(tm);
    ne = neigh(tm,true);
}

emili::LocalSearch* prs::ParamsParser::ch6_params(prs::TokenManager& tm)
{
    in = init(tm);
    te = term(tm);
    ne = neigh(tm,true);
    emili::Neighborhood* ne2 = neigh(tm,true);
    return new emili::pfsp::CH6(*in,*te,*ne,*ne2);

}

emili::LocalSearch* prs::ParamsParser::vparams(prs::TokenManager& tm)
{
    prs::incrementTabLevel();
    emili::LocalSearch* ls;
    if(tm.checkToken(FIRST))
    {
        printTab("FIRST IMPROVEMENT VND");
        in = init(tm);
        te = term(tm);
        neighs(tm);
        ls =  new emili::VNDSearch<emili::FirstImprovementSearch>(*in,*te,nes);
    }
    else if(tm.checkToken(BEST))
    {
       printTab("BEST IMPROVEMENT VND");
       in = init(tm);
       te = term(tm);
       neighs(tm);
        ls =  new emili::VNDSearch<emili::BestImprovementSearch>(*in,*te,nes);
    }
    else
    {
        std::cerr<< "'" << *tm << "' -> ERROR a valid type of search must be specified (first,best) " << std::endl;

        std::cout << info() << std::endl;
        exit(-1);
    }
    prs::decrementTabLevel();
    return ls;
}

emili::InitialSolution* prs::ParamsParser::init(prs::TokenManager& tm)
{
    prs::incrementTabLevel();
    std::ostringstream oss;
    emili::InitialSolution* init;
    if(tm.checkToken(INITIAL_RANDOM))
    {
        printTab("Random initial solution");
        init = new emili::pfsp::PfspRandomInitialSolution(*instance);
    }else if(tm.checkToken(INITIAL_RANDOM_ITERATED))
    {
        printTab("Random initial solution");
        int n = tm.getInteger();
        init = new emili::pfsp::RandomInitialSolution(*instance,n);
    }
    else if(tm.checkToken(INITIAL_SLACK))
    {
        printTab("SLACK initial solution");
        init = new emili::pfsp::PfspSlackInitialSolution(*instance);
    }else if(tm.checkToken(INITIAL_WNSLACK))
    {
        printTab( "NEH WSLACK initial solution");
        //init = new testIS(istance);
        init = new emili::pfsp::PfspNEHwslackInitialSolution(*instance);
    }
    else if(tm.checkToken(INITIAL_LIT))
        {
            printTab( "Less idle times initial solution");
            //return new testIS(istance);
            init = new emili::pfsp::LITSolution(*instance);
        }
    else if(tm.checkToken(INITIAL_RZ))
        {
            printTab( "rz initial solution");
            //return new testIS(istance);
            init = new emili::pfsp::RZSolution(*instance);
        }
    else if(tm.checkToken(INITIAL_NRZ))
        {
            printTab( "neh rz initial solution");
            //return new testIS(istance);
            init = new emili::pfsp::NeRZSolution(*instance);
        }
    else if(tm.checkToken(INITIAL_NRZ2))
        {
            printTab( "neh rz initial solution without improvement phase");
            //return new testIS(*istance);
            init = new emili::pfsp::NeRZ2Solution(*instance);
        }
    else if(tm.checkToken(INITIAL_SRZ))
        {
            printTab( "srz intial solution generator");
            //return new testIS(*istance);
            init = new emili::pfsp::SRZSolution(*instance);
        }
    else if(tm.checkToken(INITIAL_NRZ2FF))
        {
            printTab( "nehff rz initial solution without improvement phase");
            //return new testIS(*istance);
            init = new emili::pfsp::NfRZ2Solution(*instance);
        }
    else if(tm.checkToken(INITIAL_LR))
        {
            int n = tm.getInteger();
            oss.str(""); oss << "LR initial solution with "<<n<<" starting sequences";
            printTab(oss.str().c_str());
            // testIS(*istance);
            init = new emili::pfsp::LRSolution(*instance,n);
        }
    else if(tm.checkToken(INITIAL_LR_NM))
        {
            int n = instance->getNjobs()/instance->getNmachines();
            oss.str(""); oss << "LR initial solution with "<<n<<" starting sequences";
            printTab(oss.str().c_str());
            // testIS(*istance);
            init = new emili::pfsp::LRSolution(*instance,n);
        }
    else if(tm.checkToken(INITIAL_NLR))
        {
        int n = tm.getInteger();
        oss.str(""); oss << "NLR initial solution with "<<n<<" starting sequences";
        //return new testIS(*istance);printTab(oss.str().c_str());
        printTab(oss.str().c_str());
        init = new emili::pfsp::NLRSolution(*instance,n);
        }
    else if(tm.checkToken(INITIAL_MNEH))
        {
            printTab( "mneh initial solution");
            //return new testIS(istance);
            init = new emili::pfsp::MNEH(*instance);
        }
    else if(tm.checkToken(INITIAL_NEH))
    {
        printTab( "NEH initial solution");
        //return new testIS(istance);
        init = new emili::pfsp::NEH(*instance);
    }
    else if(tm.checkToken(INITIAL_NEHRS))
    {
        printTab( "NEHRS (random restarts) initial solution");
        //return new testIS(istance);
        int iterations = tm.getInteger();
        oss.str("");oss<<"number of restarts: " << iterations;
        printTab(oss.str().c_str());
        init = new emili::pfsp::NEHRS(*instance,iterations);
    }
    else if(tm.checkToken(INITIAL_NEHEDD))
    {
        printTab( "NEHedd initial solution");
        //return new testIS(istance);
        init = new emili::pfsp::NEHedd(*instance);
    }
    else if(tm.checkToken(INITIAL_NEHFF))
    {
        printTab( "NEHFF initial solution");
        //return new testIS(istance);
        init = new emili::pfsp::NEHff(*instance);
    }
    else if(tm.checkToken(INITIAL_NEHLS))
    {
        printTab( "NEHls initial solution");
        PfspInstance pfs = this->instance->getInstance();
        emili::pfsp::PermutationFlowShop * pfse = instantiateProblem(problem_type,pfs);
        emili::pfsp::PermutationFlowShop* is = this->instance;
        this->instance = pfse;
        emili::LocalSearch* ll = search(tm);
        this->instance = is;
        istances.push_back(pfse);
        init = new emili::pfsp::NEHls(*instance,ll);
    }
    else if(tm.checkToken(INITIAL_NEHEDDLS))
    {
        printTab( "NEHls initial solution");
        PfspInstance pfs = instance->getInstance();
        emili::pfsp::PermutationFlowShop * pfse = instantiateProblem(problem_type,pfs);
        emili::pfsp::PermutationFlowShop* is = this->instance;
        this->instance = pfse;
        emili::LocalSearch* ll = search(tm);
        this->instance = is;
        init = new emili::pfsp::NEHeddLS(*instance,ll);
    }
    else if(tm.checkToken(INITIAL_FRB5))
    {
        printTab( "FRB5 initial solution");
        PfspInstance pfs = this->instance->getInstance();
        emili::pfsp::PermutationFlowShop * pfse = instantiateProblem(problem_type,pfs);
        emili::InitialSolution* in = new emili::pfsp::PfspRandomInitialSolution(*pfse);
        emili::Termination* term = new emili::LocalMinimaTermination();
        emili::Neighborhood* nei = new emili::pfsp::TaillardAcceleratedInsertNeighborhood(*pfse);
        emili::LocalSearch* ll = new emili::FirstImprovementSearch(*in,*term,*nei);
        istances.push_back(pfse);
        init = new emili::pfsp::NEHls(*instance,ll);
    }
    else if(tm.checkToken(INITIAL_NEHFFLS))
    {
        printTab( "NEHffls initial solution");
        PfspInstance pfs = this->instance->getInstance();
        emili::pfsp::PermutationFlowShop * pfse = instantiateProblem(problem_type,pfs);
        emili::pfsp::PermutationFlowShop* is = this->instance;
        this->instance = pfse;
        emili::LocalSearch* ll = search(tm);
        this->instance = is;
        istances.push_back(pfse);
        init = new emili::pfsp::NEHffls(*instance,ll);
    }
    else
    {
        std::cerr<< "'" << *tm << "' -> ERROR a initial solution generator specification was expected! (random,slack)" << std::endl;

        std::cout << info() << std::endl;
        exit(-1);
    }
    prs::decrementTabLevel();
    return init;
}

emili::Termination* prs::ParamsParser::term(prs::TokenManager& tm)
{
    prs::incrementTabLevel();
    emili::Termination* term;
    if(tm.checkToken(TERMINATION_LOCMIN))
    {
        printTab("Local minima termination");
        term = new emili::LocalMinimaTermination();
    }
    else if(tm.checkToken(TERMINATION_WTRUE))
    {
        printTab("While true termination");
        term = new emili::WhileTrueTermination();
    }
    else if(tm.checkToken(TERMINATION_ITERA))
    {

        int ti = tm.getInteger();
        std::ostringstream oss;
        oss << "Relaxed local minima termination. number of max iterations "<< ti;
        printTab(oss.str().c_str());
        term =  new emili::pfsp::PfspTerminationIterations(ti);
    }
    else if(tm.checkToken(TERMINATION_SOA))
    {
        printTab("Max iteration number termination");
        int ti = instance->getNjobs();
         ti = 2*(ti-1);
        term =  new emili::pfsp::SOAtermination(ti);
    }
    else if(tm.checkToken(TERMINATION_TIME))
    {

        float time =tm.getDecimal();
        if(time==0){
            time = 1;
        }
        std::ostringstream oss;
        oss << "Timed termination. ratio: " << time;
        printTab(oss.str().c_str());
        term =  new emili::TimedTermination(time);
    }
    else if(tm.checkToken(TERMINATION_MAXSTEPS))
    {
        int steps = tm.getInteger();
        std::ostringstream oss;
        oss << "Max Steps termination. # steps: "<< steps;
        printTab(oss.str().c_str());
        term = new emili::MaxStepsTermination(steps);
    }
    else if(tm.checkToken(TERMINATION_MAXSTEPS_OR_LOCMIN))
    {
        int steps = tm.getInteger();
        std::ostringstream oss;
        oss << "Max Steps termination or when reaching locmin. # steps: "<< steps;
        printTab(oss.str().c_str());
        term = new emili::MaxStepsOrLocmin(steps);
    }
    else
    {
        std::cerr<< "'" << *tm << "' -> ERROR a termination criteria specification was expected! " << std::endl;
        std::cout << info() << std::endl;
        exit(-1);
    }
    prs::decrementTabLevel();
    return term;
}

emili::Neighborhood* prs::ParamsParser::neigh(prs::TokenManager& tm,bool checkExist)
{
    prs::incrementTabLevel();
    emili::Neighborhood* neigh = nullptr;
    if(tm.checkToken(NEIGHBORHOOD_INSERT))
    {
        printTab( "Insert Neighborhood");
        neigh = new emili::pfsp::PfspInsertNeighborhood(*instance);
    }
    else if(tm.checkToken(NEIGHBORHOOD_ADAPTIVE_INSERT))
    {
        printTab( "Insert Neighborhood");
        neigh = new emili::pfsp::PfspInsertNeighborhood(*instance);
    }
    else  if(tm.checkToken(NEIGHBORHOOD_FORW_INSERT))
    {
        printTab( "Forward insert Neighborhood");
        neigh = new emili::pfsp::PfspForwardInsertNeighborhood(*instance);
    }
    else  if(tm.checkToken(NEIGHBORHOOD_BACK_INSERT))
    {
        printTab( "Backward Insert Neighborhood");
        neigh = new emili::pfsp::PfspBackwardInsertNeighborhood(*instance);
    }
    else if(tm.checkToken(NEIGHBORHOOD_EXCHANGE))
    {
        printTab( "Exchange neighborhood");
        neigh = new emili::pfsp::PfspExchangeNeighborhood(*instance);
    }
    else if(tm.checkToken(NEIGHBORHOOD_ATX_EXCHANGE))
    {
        printTab( "Exchange neighborhood");
        neigh = new emili::pfsp::AxtExchange(*instance);
    }
    else if(tm.checkToken(NEIGHBORHOOD_OPT_EXCHANGE))
    {
        printTab( "Optimized Exchange neighborhood");
        neigh = new emili::pfsp::OptExchange(*instance);
    }
    else if(tm.checkToken(NEIGHBORHOOD_HATX_EXCHANGE))
       {
           printTab( "Exchange neighborhood with speedup");
           neigh = new emili::pfsp::HaxtExchange(*instance);
       }
       else if(tm.checkToken(NEIGHBORHOOD_EATX_EXCHANGE))
       {
           printTab( "Exchange neighborhood with speedup");
           neigh = new emili::pfsp::EaxtExchange(*instance);
       }
    else if(tm.checkToken(NEIGHBORHOOD_TRANSPOSE))
    {
        printTab( "Transpose neighborhood");
        neigh = new emili::pfsp::PfspTransposeNeighborhood(*instance);
    }
    else if(tm.checkToken(NEIGHBORHOOD_TWO_INSERT))
    {
        printTab( "Two insert neighborhood");
        neigh = new emili::pfsp::PfspTwoInsertNeighborhood(*instance);
    }
    else if(tm.checkToken(NEIGHBORHOOD_XTRANSPOSE))
    {
        printTab( "XTranspose neighborhood");
        neigh = new emili::pfsp::XTransposeNeighborhood(*instance);
    }
    else if(tm.checkToken(NEIGHBORHOOD_TA_INSERT))
    {
        printTab( "Insert with Taillard Acceleration");
        neigh = new emili::pfsp::TaillardAcceleratedInsertNeighborhood(*instance);
    }
    else if(tm.checkToken(NEIGHBORHOOD_FSTA_INSERT))
    {
        printTab( "Insert with Taillard Acceleration that updates the base solution after each improvement");
        neigh = new emili::pfsp::FSTaillardAcceleratedInsertNeighborhood(*instance);
    }
    else if(tm.checkToken(NEIGHBORHOOD_CSTA_INSERT))
    {
        printTab( "Insert with Taillard Acceleration that evaluates all the possible insertion points");
        neigh = new emili::pfsp::CSTaillardAcceleratedInsertNeighborhood(*instance);
    }
    else if(tm.checkToken(NEIGHBORHOOD_TAx_INSERT))
    {
        printTab( "Insert with Taillard Acceleration(Experimental)");
        neigh = new emili::pfsp::TAxInsertNeighborhood(*instance);
    }
    else if(tm.checkToken(NEIGHBORHOOD_OPT_INSERT))
    {
        printTab( "Delta Evaluation Insert for Weighted Tardiness with tail improvement");
        neigh = new emili::pfsp::OptInsert(*instance);
    }
    else if(tm.checkToken(NEIGHBORHOOD_ATAx_INSERT))
    {
        printTab( "Atx Delta Evaluation Insert for Weighted Tardiness");
        neigh = new emili::pfsp::AtxNeighborhood(*instance);
    }
    else if(tm.checkToken(NEIGHBORHOOD_HATAx_INSERT))
    {
        printTab( "Approximated Insert with Taillard Acceleration for Weighted Tardiness with no threshold");
        neigh = new emili::pfsp::HeavilyApproximatedTaillardAcceleratedInsertNeighborhood(*instance);
    }
    else if(tm.checkToken(NEIGHBORHOOD_NATAx_INSERT))
    {
        printTab( "Approximated Insert for Weighted Tardiness with 1 level approximation");
        neigh = new emili::pfsp::NatxNeighborhood(*instance);
    }
    else if(tm.checkToken(NEIGHBORHOOD_NATA2x_INSERT))
    {
        printTab( "Improved Approximated Insert for Weighted Tardiness with 1 level approximation and online tuned threshold");
        neigh = new emili::pfsp::Natx2Neighborhood(*instance);
    }
    else if(tm.checkToken(NEIGHBORHOOD_EATAx_INSERT))
    {
        printTab( "Approximated Insert for Weighted Tardiness with 2 levels of approximation");
        neigh = new emili::pfsp::EatxNeighborhood(*instance);
    }
    else if(tm.checkToken(NEIGHBORHOOD_THATAx_INSERT))
    {
        printTab( "Approximated Insert for Weighted Tardiness with 3 levels of approximation");
        neigh = new emili::pfsp::ThatxNeighborhood(*instance);
    }
    else if(tm.checkToken(NEIGHBORHOOD_PATAx_INSERT))
    {
        printTab( "Approximated Insert for Weighted Tardiness with 5 levels of approximation");
        neigh = new emili::pfsp::PatxNeighborhood(*instance);
    }
    else if(tm.checkToken(NEIGHBORHOOD_SATAx_INSERT))
    {
        printTab( "Approximated Insert for Weighted Tardiness with 6 levels of approximation");
        neigh = new emili::pfsp::SatxNeighborhood(*instance);
    }
    else if(tm.checkToken(NEIGHBORHOOD_FATAx_INSERT))
    {
        printTab( "Approximated Insert for Weighted Tardiness with 4 levels of approximation");
        neigh = new emili::pfsp::FatxNeighborhood(*instance);
    }
    else if(tm.checkToken(NEIGHBORHOOD_TATAx_INSERT))
    {
        printTab( "Approximated Insert for Weighted Tardiness with settable threshold");
        float start_level = tm.getDecimal();       
        neigh = new emili::pfsp::TatxNeighborhood(start_level,*instance);
    }
    else if(tm.checkToken(NEIGHBORHOOD_NITA_INSERT))
    {
        printTab( "Insert with Taillard Acceleration for no idle make span ");
        neigh = new emili::pfsp::NoIdleAcceleratedInsertNeighborhood(*instance);
    }
    else if(tm.checkToken(NEIGHBORHOOD_NATA_TCT_INSERT))
    {
        printTab( "Improved Approximated Insert for Total Completion Times with 1 level approximation and online tuned threshold");
        neigh = new emili::pfsp::NatxTCTNeighborhood(*instance);
    }
    else if(tm.checkToken(NEIGHBORHOOD_RZ_TCT_INSERT))
    {
        printTab( "iRZ neighborhood see PanRui2012");
        neigh = new emili::pfsp::NrzTCTNeighborhood(*instance);
    }
    else if(tm.checkToken(NEIGHBORHOOD_NATA_TT_INSERT))
    {
        printTab( "Improved Approximated Insert for Total Tardiness with 1 level approximation and online tuned threshold");
        neigh = new emili::pfsp::NatxTTNeighborhood(*instance);
    }
    else if(tm.checkToken(NEIGHBORHOOD_RANDCONHE))
    {
         printTab("Random Constructive Heuristic Neighborhood ");
         emili::InitialSolution* in = init(tm);
         neigh = new emili::RandomConstructiveHeuristicNeighborhood(*in);
    }
    else if(tm.checkToken(NEIGHBORHOOD_SDSTTA_INSERT))
    {
        printTab( "Taillard acceleration for Sequence dependent setup times");
        neigh = new emili::pfsp::SDSTTaillardAcceleratedInsertNeighborhood(*instance);
    }
    else
    {
        if(checkExist)
        {
            std::cerr<< "'" << *tm << "' -> ERROR a neighborhood specification was expected! " << std::endl;
            std::cout << info() << std::endl;
            exit(-1);
        }
    }
    prs::decrementTabLevel();
    return neigh;
}

void prs::ParamsParser::neighs(prs::TokenManager& tm)
{
    std::vector<emili::Neighborhood*> vnds;
    vnds.push_back(neigh(tm,true));
    nes = vnds;
    neighs1(tm);
}

void prs::ParamsParser::neighs1(prs::TokenManager& tm)
{
    emili::Neighborhood* n = neigh(tm,false);
	if(n!=nullptr)
	{
           nes.push_back(n);
           neighs1(tm);
	}

}


void prs::ParamsParser::problem(prs::TokenManager& tm)
{
    PfspInstance i;    
    problem_type = tm.nextToken();
    bool ok;
    std::string pro(problem_type);
    std::string sdst("SDSTPFS");

    if(pro.find(sdst) != std::string::npos)
    {
        ok = i.readSeqDepDataFromFile(tm.tokenAt(1));
    }
    else
    {
        ok = i.readDataFromFile(tm.tokenAt(1));

    }

    if(ok)
     {
         instance = instantiateProblem(problem_type, i);
        return;
     }

        std::cout << info() << std::endl;
        exit(-1);
}
#include "pfspBuilder.h"

emili::LocalSearch* prs::ParamsParser::buildAlgo(prs::TokenManager& tm)
{
    /*tm.move(0);
    prs::GeneralParserE  ps(tm);
    prs::EmBaseBuilder emb(ps,ps.getTokenManager());
    prs::PfspBuilder pfspb(ps,ps.getTokenManager());
    ps.addBuilder(&emb);
    ps.addBuilder(&pfspb);
    emili::LocalSearch* local = ps.parseParams();*/
    problem(tm);  
    emili::LocalSearch* local = eparams(tm);
    std::cout << "------" << std::endl;
    return local;
}

bool prs::ParamsParser::isParsable(std::string &problem)
{
    if(strcmp(problem.c_str(),PROBLEM_PFS_WT)==0)
    {
        return true;
    }else if(strcmp(problem.c_str(),PROBLEM_PFS_E)==0)
    {
        return true;
    }else if(strcmp(problem.c_str(),PROBLEM_PFS_WE)==0)
    {
        return true;
    }else if(strcmp(problem.c_str(),PROBLEM_PFS_T)==0)
    {
        return true;
    }else if(strcmp(problem.c_str(),PROBLEM_PFS_MS)==0)
    {
        return true;
    }
    else if(strcmp(problem.c_str(),PROBLEM_PFS_TCT)==0)
    {
        return true;
    }else if(strcmp(problem.c_str(),PROBLEM_PFS_WCT)==0)
    {
        return true;
    }
    else if(strcmp(problem.c_str(),PROBLEM_NIPFS_MS)==0)
    {
        return true;
    }
    else if(strcmp(problem.c_str(),PROBLEM_NIPFS_E)==0)
    {
        return true;
    }else if(strcmp(problem.c_str(),PROBLEM_NIPFS_T)==0)
    {
        return true;
    }
    else if(strcmp(problem.c_str(),PROBLEM_NIPFS_WT)==0)
    {
        return true;
    }
    else if(strcmp(problem.c_str(),PROBLEM_NIPFS_WE)==0)
    {
        return true;
    }
    else if(strcmp(problem.c_str(),PROBLEM_NIPFS_TCT)==0)
    {
        return true;
    }else if(strcmp(problem.c_str(),PROBLEM_NIPFS_WCT)==0)
    {
        return true;
    }
    else if(strcmp(problem.c_str(),PROBLEM_NWPFS_MS)==0)
    {
        return true;
    }
    else if(strcmp(problem.c_str(),PROBLEM_NWPFS_E)==0)
    {
        return true;
    }else if(strcmp(problem.c_str(),PROBLEM_NWPFS_T)==0)
    {
        return true;
    }
    else if(strcmp(problem.c_str(),PROBLEM_NWPFS_WT)==0)
    {
        return true;
    }
    else if(strcmp(problem.c_str(),PROBLEM_NWPFS_WE)==0)
    {
        return true;
    }
    else if(strcmp(problem.c_str(),PROBLEM_NWPFS_TCT)==0)
    {
        return true;
    }else if(strcmp(problem.c_str(),PROBLEM_NWPFS_WCT)==0)
    {
        return true;
    }
    else if(strcmp(problem.c_str(),PROBLEM_SDSTPFS_MS)==0)
    {
        return true;
    }
    else if(strcmp(problem.c_str(),PROBLEM_SDSTPFS_E)==0)
    {
        return true;
    }else if(strcmp(problem.c_str(),PROBLEM_SDSTPFS_T)==0)
    {
        return true;
    }
    else if(strcmp(problem.c_str(),PROBLEM_SDSTPFS_WT)==0)
    {
        return true;
    }
    else if(strcmp(problem.c_str(),PROBLEM_SDSTPFS_WE)==0)
    {
        return true;
    }
    else if(strcmp(problem.c_str(),PROBLEM_SDSTPFS_TCT)==0)
    {
        return true;
    }else if(strcmp(problem.c_str(),PROBLEM_SDSTPFS_WCT)==0)
    {
        return true;
    }
    else
    {
        return false;
    }
}

std::string prs::ParamsParser::availableProblems() const
{
    std::ostringstream oss;
    oss <<PROBLEM_PFS_WT<< " " <<PROBLEM_PFS_WE<< " " <<PROBLEM_PFS_TCT<< " " <<PROBLEM_PFS_MS
       << " " <<PROBLEM_PFS_WCT<< " " <<PROBLEM_PFS_T<< " " <<PROBLEM_PFS_E<< " "<<PROBLEM_NWPFS_WT<< " " <<PROBLEM_NWPFS_WE
       << " " <<PROBLEM_NWPFS_TCT<< " " <<PROBLEM_NWPFS_MS<< " " <<PROBLEM_NWPFS_T<< " " <<PROBLEM_NWPFS_E
       << PROBLEM_NIPFS_MS <<" "<<PROBLEM_NIPFS_WT<< " " <<PROBLEM_NIPFS_WE<< " " <<PROBLEM_NIPFS_TCT<< " " <<PROBLEM_NIPFS_MS<< " "
       << " " <<PROBLEM_NIPFS_T<< " " <<PROBLEM_NIPFS_E << " " << PROBLEM_SDSTPFS_MS;

    return oss.str();
}

/*SAInitTemp* prs::ParamsParser::INITTEMP(prs::TokenManager& tm,
                                   emili::InitialSolution *initsol,
                                   emili::Neighborhood *nei,
                                   emili::pfsp::PermutationFlowShop *instance) {

    if (tm.checkToken(FIXEDINITTEMP)) {
        double             value     = tm.getDecimal();
        SAInitTemp* init_temp = new FixedInitTemp();
        init_temp->set(value);
        return init_temp;
    } else if (tm.checkToken(INITFROMSOL)) {
        double                  value    = tm.getDecimal();
        emili::Solution* is = initsol->generateSolution();
        SAInitTemp* init_temp     = new InitTempFromSolution(is);
        init_temp->set(value);
        return init_temp;
    } else if (tm.checkToken(RANDOMWALKINITTEMP)) {
        int length = tm.getInteger();
        double value = tm.getDecimal();
        SAInitTemp* init_temp = new RandomWalkInitTemp(initsol, nei, length);
        init_temp->set(value);
        return init_temp;
    } else if (tm.checkToken(CONNOLLYRWIT)) {
        int length = tm.getInteger();
        double value = tm.getDecimal();
        SAInitTemp* init_temp = new ConnollyRandomWalkInitTemp(initsol, nei, length);
        init_temp->set(value);
        return init_temp;
    } else if (tm.checkToken(RANDOMWALKAVGINITTEMP)) {
        int length = tm.getInteger();
        double value = tm.getDecimal();
        SAInitTemp* init_temp = new RandomWalkAvgInitTemp(initsol, nei, length);
        init_temp->set(value);
        return init_temp;
    } else if (tm.checkToken(RANDOMWALKINITPROB)) {
        float ip = tm.getDecimal();
        int length = tm.getInteger();
        double value = tm.getDecimal();
        SAInitTemp* init_temp = new RandomWalkInitProb(initsol, nei, ip, length);
        init_temp->set(value);
        return init_temp;
    } else if (tm.checkToken(MISEVICIUSINITTEMP)) {
        int length = tm.getInteger();
        float l11 = tm.getDecimal();
        float l12 = tm.getDecimal();
        float l21 = tm.getDecimal();
        float l22 = tm.getDecimal();
        SAInitTemp* init_temp = new MiseviciusInitTemp(initsol, nei, length, l11, l12, l21, l22);
        init_temp->set(1);
        return init_temp;
    } else if (tm.checkToken(SIMPLEMISEVICIUSINITTEMP)) {
        int length = tm.getInteger();
        float l1 = tm.getDecimal();
        float l2 = tm.getDecimal();
        SAInitTemp* init_temp = new SimplifiedMiseviciusInitTemp(initsol, nei, length, l1, l2);
        init_temp->set(1);
        return init_temp;
    } / *else if (tm.checkToken(OSMANPOTTSINITTEMP)) {
        float dc = tm.getDecimal();
        float tf = tm.getDecimal();
        float coeff = tm.getDecimal();
        SAInitTemp* init_temp = new OsmanPottsInitTemp(initsol, nei, instance, dc, tf);
        init_temp->set(coeff);
        return init_temp;
    }  else if (tm.checkToken(RANDOMWALKSTATSINITTEMP)) {
        int length = tm.getInteger();
        double value = tm.getDecimal();
        SAInitTemp* init_temp = new RandomWalkStatsInitTemp(initsol, nei, length);
        init_temp->set(value);
        return init_temp;
    } * / else {
        std::cerr << "SAInitTemp expected, not found : " << std::endl;
        std::cerr << tm.peek() << std::endl;
        exit(1);
    }

}


SAAcceptance* prs::ParamsParser::ACCEPTANCE(prs::TokenManager& tm,
                                      SAInitTemp *inittemp,
                                 emili::Neighborhood *nei,
                                 emili::Problem* instance) {

    if (tm.checkToken(METROPOLIS)) {
        return new SAMetropolisAcceptance(inittemp->get());
    } else if (tm.checkToken(METROPOLISWFORCED)) {
        return new SAMetropolisWithForcedAcceptance(inittemp->get());
    } else if (tm.checkToken(BASICACC)) {
        return new SABasicAcceptance();
    } else if (tm.checkToken(APPROXEXPACC)) {
        return new SAApproxExpAcceptance(inittemp->get());
    } else if (tm.checkToken(GEOMACC)) {
        double rf = tm.getDecimal();
        return new SAGeometricAcceptance(inittemp->get(), rf);
    } else if (tm.checkToken(GENSAACC)) {
        double beta = tm.getDecimal();
        double g = tm.getDecimal();
        return new GeneralizedSAAcceptance(inittemp->get(), beta, g);
    } else if (tm.checkToken(DETERMINISTICACC)) {
        return new SADeterministicAcceptance(inittemp->get());
    } else if (tm.checkToken(GDAACC)) {
        return new GreatDelugeAcceptance();
    } else if (tm.checkToken(RTRACC)) {
        double de = tm.getDecimal();
        return new RecordToRecordAcceptance(de);
    } else if (tm.checkToken(LAHCACC)) {
        double te = tm.getInteger();
        return new LAHCAcceptance(te);
    } else if (tm.checkToken(LAHCNSACC)) {
        double te = tm.getDecimal();
        return new LAHCNSAcceptance(te, nei);
    } else if (tm.checkToken(PRECOMPUTEDMETROPOLIS)) {
        int te = tm.getInteger();
        return new SAPrecomputedMetropolisAcceptance(inittemp->get(), te);
    } else if (tm.checkToken(PRECOMPUTEDMETROPOLISWFORCED)) {
        int te = tm.getInteger();
        return new SAPrecomputedMetropolisWithForcedAcceptance(inittemp->get(), te);
    } else if (tm.checkToken(BOUNDEDMETROPOLIS)) {
        double rd = tm.getDecimal();
        return new SABoundedMetropolisAcceptance(inittemp->get(), rd);
    } else if (tm.checkToken(ALLACC)) {
        return new SAAcceptanceAll();
    } else {
        std::cerr << "SAAcceptance expected, not found : " << std::endl;
        std::cerr << tm.peek() << std::endl;
        exit(1);
    }

}


SATermination* prs::ParamsParser::TERMINATION(prs::TokenManager& tm,
                                              SAInitTemp* inittemp,
                                              emili::Neighborhood *nei) {

    if (tm.checkToken(MAXBADITERS)) {
        int    mb = tm.getInteger();
        return new SAMaxBadIterTermination(mb);
    } else if (tm.checkToken(MAXITERS)) {
        int    mi = tm.getInteger();
        return new SAMaxIterTermination(mi);
    } else if (tm.checkToken(NEVERTERM)) {
        return new SAWhileTrueTermination();
    } else if (tm.checkToken(ACCRATETERM)) {
        float rate = tm.getDecimal();
        return new SAAcceptanceRateTermination(rate);
    } else if (tm.checkToken(LASTACCRATETERM)) {
        int te = tm.getInteger();
        float rate = tm.getDecimal();
        return new SALastAcceptanceRateTermination(te, rate);
    } else if (tm.checkToken(MAXTEMPRESTARTSTERM)) {
        int tr = tm.getInteger();
        return new SaMaxTempRestartsTermination(tr);
    } else if (tm.checkToken(NEIGHSIZEITERTERM)) {
        float co = tm.getDecimal();
        return new SANeighSizeIterTermination(nei, co);
    } else if (tm.checkToken(SQUAREDNSITERTERM)) {
        float co = tm.getDecimal();
        return new SASquaredNeighSizeIterTermination(nei, co);
    } else if (tm.checkToken(LOCALMINTERM)) {
        int te = tm.getInteger();
        return new SALocalMinTermination(te);
    } else if (tm.checkToken(NEIGHSIZELOCALMINTERM)) {
        float co = tm.getDecimal();
        return new SANeighSizeLocalMinTermination(nei, co);
    } else if (tm.checkToken(MAXSTEPSTERM)) {
        int ms = tm.getInteger();
        return new SAMaxStepsTermination(ms);
    } else if (tm.checkToken(MINTEMPTERM)) {
        char* p = tm.peek();
        double v;
        try {
            v = std::stod(p);
        } catch(...) {
            return new SAMinTempTermination(inittemp->getMinTemp());
        }
        v = tm.getDecimal();
        return new SAMinTempTermination(v);
    } else {
        std::cerr << "SATermination expected, not found : " << std::endl;
        std::cerr << tm.peek() << std::endl;
        exit(1);
    }


}


SACooling* prs::ParamsParser::COOL(prs::TokenManager& tm,
                              SAInitTemp *it,
                              emili::Neighborhood *nei,
                              emili::pfsp::PermutationFlowShop *instance) {

    if (tm.checkToken(GEOM)) {
        float a = tm.getDecimal();
        float b = tm.getDecimal();
        return new GeomCooling(a,b, it);
    } else if (tm.checkToken(MARKOV)) {
        float b = tm.getDecimal();
        return new MarkovCooling(b, it);
    } else if (tm.checkToken(LOGCOOLING)) {
        float b = tm.getDecimal();
        return new LogCooling(b, it);
    } else if (tm.checkToken(CONSTCOOLING)) {
        float a = tm.getDecimal();
        float b = tm.getDecimal();
        return new ConstantCooling(a,b, it);
    } else if (tm.checkToken(LUNDYMEES)) {
        float a = tm.getDecimal();
        float b = tm.getDecimal();
        return new LundyMeesCooling(a,b, it);
    } else if (tm.checkToken(LUNDYMEESCONNOLLY)) {
        return new LundyMeesConnollyCooling(it);
    } else if (tm.checkToken(Q87COOLING)) {
        float a = tm.getDecimal();
        float b = tm.getDecimal();
        int   m = tm.getInteger();
        return new Q87Cooling(a, b, m, it);
    } else if (tm.checkToken(LINEARCOOLING)) {
        float a = tm.getDecimal();
        return new LinearCooling(a, it);
    } else if (tm.checkToken(NOCOOLING)) {
        return new NoCooling(it);
    } else if (tm.checkToken(TEMPBANDCOOLING)) {
        float a = tm.getDecimal();
        return new SATemperatureBandCooling(a, it);
    } else if (tm.checkToken(QUADRATICCOOLING)) {
        return new SAQuadraticCooling(it);
    } else if (tm.checkToken(OSMANPOTTSPFSPCOOLING)) {
        return new OsmanPottsPFSPCooling(it, instance);
    } else if (tm.checkToken(ARITHMETICCOOLING)) {
        double a = tm.getDecimal();
        return new ArithmeticCooling(a, it);
    } else if (tm.checkToken(OBA1)) {
        long   M = tm.getInteger();
        long   delta = tm.getDecimal();
        float a = tm.getDecimal();
        float b = tm.getDecimal();
        float c = tm.getDecimal();
        return new OldBachelor1(M, delta, a, b, c, it, instance);
    } else if (tm.checkToken(OBA2)) {
        long   M = tm.getInteger();
        long   delta = tm.getDecimal();
        float  d = tm.getDecimal();
        return new OldBachelor2(M, delta, d, it, nei);
    } else {
        std::cerr << "SACooling expected, not found : " << std::endl;
        std::cerr << tm.peek() << std::endl;
        exit(1);
    }

}


SATempLength* prs::ParamsParser::TEMPLENGTH(prs::TokenManager& tm,
                                       emili::Neighborhood* neigh,
                                       emili::Problem* instance) {

    if (tm.checkToken(CONSTTEMPLEN)) {
        int l = tm.getInteger();
        return new ConstantTempLength(l);
    } else if (tm.checkToken(NEIGHSIZETEMPLEN)) {
        float a = tm.getDecimal();
        return new NeighSizeTempLength(neigh, a);
    } else if (tm.checkToken(PROBSIZETEMPLEN)) {
        float a = tm.getDecimal();
        return new ProblemSizeTempLength(instance, a);
    } else if (tm.checkToken(SQUAREDPROBSIZETEMPLEN)) {
        float a = tm.getDecimal();
        return new SquaredProblemSizeTempLength(instance, a);
    } else if (tm.checkToken(BRNEIGHSIZETEMPLEN)) {
        float a = tm.getDecimal();
        return new BurkardRendlNeighSizeTempLength(neigh, a);
    } else if (tm.checkToken(MAXACCEPTEDTEMPLEN)) {
        int a = tm.getInteger();
        return new MaxAcceptedTempLength(a);
    } else if (tm.checkToken(CAPPEDMAXACCEPTEDTEMPLEN)) {
        int a = tm.getInteger();
        int c = tm.getInteger();
        return new CappedMaxAcceptedTempLength(a, c);
    } else if (tm.checkToken(NEIGHCAPPEDMAXACCEPTEDTEMPLEN)) {
        float a = tm.getDecimal();
        int c = tm.getInteger();
        return new NeighSizeCappedMaxAcceptedTempLength(neigh, a, c);
    } else if (tm.checkToken(ARITMTEMPLEN)) {
        int a = tm.getInteger();
        int b = tm.getInteger();
        return new ArithmeticTempLength(a, b);
    } else if (tm.checkToken(GEOMTEMPLEN)) {
        int a = tm.getInteger();
        float b = tm.getDecimal();
        return new GeomTempLength(a, b);
    } else if (tm.checkToken(LOGTEMPLEN)) {
        int a = tm.getInteger();
        int b = tm.getInteger();
        return new LogTempLength(a, b);
    } else if (tm.checkToken(EXPTEMPLEN)) {
        int a = tm.getInteger();
        float b = tm.getDecimal();
        return new ExpTempLength(a, b);
    } else if (tm.checkToken(NOTEMPLEN)) {
        return new NoTempLength();
    } else if (tm.checkToken(BRGEOMTEMPLEN)) {
        float b = tm.getDecimal();
        return new BurkardRendlGeomTempLength(neigh, b);
    } else {
        std::cerr << "SATempLength expected, not found : " << std::endl;
        std::cerr << tm.peek() << std::endl;
        exit(1);
    }
}

SAExploration* prs::ParamsParser::EXPLORATION(prs::TokenManager& tm,
                                        emili::Neighborhood* neigh,
                                        SAAcceptance *acc,
                                        SACooling *cool,
                                        SATermination *term) {

    if (tm.checkToken(SARANDOMEXPLORATION)) {
        return new SARandomExploration(neigh, acc, cool, term);
    } else if (tm.checkToken(SASEQUENTIALEXPLORATION)) {
        return new SASequentialExploration(neigh, acc, cool, term);
    } else if (tm.checkToken(SABESTOFKEXPLORATION)) {
        long k = tm.getInteger();
        return new SABestOfKExploration(neigh, acc, cool, term, k);
    } else if (tm.checkToken(SABESTOFKSEQUENTIALEXPLORATION)) {
        long k = tm.getInteger();
        return new SABestOfKSequentialExploration(neigh, acc, cool, term, k);
    } else if (tm.checkToken(SANSBESTOFKSEQUENTIALEXPLORATION)) {
        double k = tm.getDecimal();
        return new SANSBestOfKSequentialExploration(neigh, acc, cool, term, k);
    } else if (tm.checkToken(SANSBESTOFKRANDOMEXPLORATION)) {
        double k = tm.getDecimal();
        return new SANSBestOfKRandomExploration(neigh, acc, cool, term, k);
    } else if (tm.checkToken(SAFIRSTBESTOFKEXPLORATION)) {
        long k = tm.getInteger();
        return new SAFirstBestOfKExploration(neigh, acc, cool, term, k);
    } else if (tm.checkToken(SANSBESTOFKEXPLORATION)) {
        double k = tm.getDecimal();
        return new SANSBestOfKExploration(neigh, acc, cool, term, k);
    } else if (tm.checkToken(SANSFIRSTBESTOFKEXPLORATION)) {
        double k = tm.getDecimal();
        return new SANSFirstBestOfKExploration(neigh, acc, cool, term, k);
    } else {
        std::cerr << "SAExploration expected, not found : " << std::endl;
        std::cerr << tm.peek() << std::endl;
        exit(1);
    }

}

SATempRestart* prs::ParamsParser::TEMPRESTART(prs::TokenManager& tm,
                                         SAInitTemp *it,
                                         emili::Neighborhood* neigh) {

    if (tm.checkToken(SANOTEMPRESTART)) {
        return new SANoRestart();
    } else if (tm.checkToken(SAMINTEMPRESTART)) {
        float va = tm.getDecimal();
        return new SAMinRestart(it, va);
    } else if (tm.checkToken(SAPERCTEMPRESTART)) {
        float va = tm.getDecimal();
        return new SAPercRestart(it, va);
    } else if (tm.checkToken(SALOWRATERESTART)) {
        float va = tm.getDecimal();
        return new SALowRateRestart(it, va);
    } else if (tm.checkToken(SALASTRATERESTART)) {
        int   te = tm.getInteger();
        float va = tm.getDecimal();
        return new SALastRateRestart(it, te, va);
    } else if (tm.checkToken(SALOWRATEREHEAT)) {
        float th = tm.getDecimal();
        float va = tm.getDecimal();
        return new SALowRateReheat(it, th, va);
    } else if (tm.checkToken(SALASTRATEREHEAT)) {
        int   te = tm.getInteger();
        float th = tm.getDecimal();
        float va = tm.getDecimal();
        return new SALastRateReheat(it, te, th, va);
    } else if (tm.checkToken(SALOCALMINREHEAT)) {
        int   te = tm.getInteger();
        float va = tm.getDecimal();
        return new SALocalMinReheat(it, te, va);
    } else if (tm.checkToken(SALOWRATERESTARTBEST)) {
        float va = tm.getDecimal();
        return new SALowRateRestartToBest(it, va);
    } else if (tm.checkToken(SALOCALMINRESTARTBEST)) {
        int   te = tm.getInteger();
        return new SALocalMinRestartToBest(it, te);
    } else if (tm.checkToken(SALOCALMINTEMPRESTART)) {
        int   te = tm.getInteger();
        return new SALocalMinTempRestart(it, te);
    } else if (tm.checkToken(SALOCALMINENHANCEDREHEAT)) {
        int   te = tm.getInteger();
        float va = tm.getDecimal();
        float ep = tm.getDecimal();
        return new SALocalMinEnhancedReheat(it, te, va, ep);
    } else if (tm.checkToken(SAMAXITERSTEMPRESTART)) {
        int   te = tm.getInteger();
        return new SAMaxItersTempRestart(it, te);
    } else if (tm.checkToken(SAMAXITERSREHEAT)) {
        int   te = tm.getInteger();
        float va = tm.getDecimal();
        return new SAMaxItersReheat(it, te, va);
    } else if (tm.checkToken(SANEIGHSIZEMAXITERSTEMPRESTART)) {
        float   te = tm.getDecimal();
        return new SANeighSizeMaxItersTempRestart(it, neigh, te);
    } else if (tm.checkToken(SASQUAREDNEIGHSIZEMAXITERSTEMPRESTART)) {
        float   te = tm.getDecimal();
        return new SASquaredNeighSizeMaxItersTempRestart(it, neigh, te);
    } else if (tm.checkToken(SANEIGHSIZEMAXITERSREHEAT)) {
        float   te = tm.getDecimal();
        float va = tm.getDecimal();
        return new SANeighSizeMaxItersReheat(it, neigh, te, va);
    } else if (tm.checkToken(SAMAXSTEPSTEMPRESTART)) {
        int   te = tm.getInteger();
        return new SAMaxStepsTempRestart(it, te);
    } else if (tm.checkToken(SAMAXSTEPSREHEAT)) {
        int   te = tm.getInteger();
        float va = tm.getDecimal();
        return new SAMaxStepsReheat(it, te, va);
    } else if (tm.checkToken(SANEIGHSIZEMAXSTEPSTEMPRESTART)) {
        float   te = tm.getDecimal();
        return new SANeighSizeMaxStepsTempRestart(it, neigh, te);
    } else if (tm.checkToken(SANEIGHSIZEMAXSTEPSREHEAT)) {
        float   te = tm.getDecimal();
        float va = tm.getDecimal();
        return new SANeighSizeMaxStepsReheat(it, neigh, te, va);
    } else {
        std::cerr << "SATempRestart expected, not found : " << std::endl;
        std::cerr << tm.peek() << std::endl;
        exit(1);
    }
}*/
