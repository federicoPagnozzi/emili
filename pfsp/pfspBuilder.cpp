//
//  Created by Federico Pagnozzi on 28/11/14.
//  Copyright (c) 2014 Federico Pagnozzi. All rights reserved.
//  This file is distributed under the BSD 2-Clause License. See LICENSE.TXT
//  for details.
#include "pfspBuilder.h"
#include <cstdlib>
#include <cstdio>
#include <ctime>
#include <cstring>
#include <iostream>
#include <sstream>
#include <algorithm>
#include "pfspinstance.h"

/* Algos */
#define CH6_LS "ch6"
#define RIS_LS "ris"
#define NI_RIS_LS "niris"
#define NW_RIS_LS "nwris"
#define RNW_RIS_LS "rnwris"
#define SWP_INC_LS "swpinc"
#define STH_LS "sthp"
#define STHF_LS "sth"

/* tabu tenure types */
#define TABU_MEMORY_MOVES "move"
#define TABU_MEMORY_MOVES2 "move2"
#define TABU_MEMORY_HASHES "hash"
#define TABU_MEMORY_SOLUTIONS "solution"
#define TABU_MEMORY_TSAB "tsabm"
#define TABU_MEMORY_TSAB_TEST "tsabmt"
#define TABU_MEMORY_VALUE "value"

/* modifiers */
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
#define INITIAL_LIT "lit"
#define INITIAL_RZ "rz"
#define INITIAL_NRZ "nrz"
#define INITIAL_NRZ2 "nrz2"
#define INITIAL_NRZ2FF "nrz2ff"
#define INITIAL_SRZ "srz"
#define INITIAL_LR "lr"
#define INITIAL_LR_NM "lrnm"
#define INITIAL_NLR "nlr"
#define INITIAL_MNEH "mneh"
#define INITIAL_RMNEH "rmneh"
#define INITIAL_WNSLACK "nwslack"
#define INITIAL_FRB5 "frb5"
#define INITIAL_CSFRB5 "csfrb5"
#define INITIAL_FRB5_GENERAL "gfrb5"
#define INITIAL_BS "bs"
#define INITIAL_BS2 "bs2"
#define INITIAL_BS2N "bs2n"
#define INITIAL_BSNN "bsnn"
#define INITIAL_BS2NF "bs2nf"
#define INITIAL_FF "ff"
#define INITIAL_FFN "ffn"
#define INITIAL_BSCHO "bscho"
#define INITIAL_BSCHR "bschr"
#define INITIAL_BSCH "bsch"


/* Termination criteria*/
#define TERMINATION_ITERA "iteration"
#define TERMINATION_SOA "soater"
#define TERMINATION_KAR "karter"
#define TERMINATION_MAXSTEPS_WITHNOIMPROV "msnoimprov"

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
#define NEIGHBORHOOD_KAR "karnghb"
#define NEIGHBORHOOD_CS_INSERT "csinsert"

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

/* No wait makespan*/
#define NEIGHBORHOOD_NW_INSERT "nwinsert"
#define NEIGHBORHOOD_NW_TWO_INSERT "nwtinsert"
#define NEIGHBORHOOD_NW_EXCHANGE "nwexchange"
#define NEIGHBORHOOD_NW_TRANSPOSE "nwtranspose"

/* Sequence Dependent Setup times makespan*/
#define NEIGHBORHOOD_SDSTTA_INSERT "sdsttainsert"
#define NEIGHBORHOOD_SDST_CS_INSERT "sdstcsinsert"

/*
 * END Neighborhoods
 */


/* permutation flowshop solution perturbations */
#define PERTURBATION_RND "randpert"
#define PERTURBATION_NRZ "nrzper"
#define PERTURBATION_TMIIG "tmiigper"
#define PERTURBATION_NWTMIIG "nwtmiigper"
#define PERTURBATION_SOA "igper"
#define PERTURBATION_SOA_LEGACY "soaper"
#define PERTURBATION_TEST "testper"
#define PERTURBATION_IGLS "igls"
#define PERTURBATION_NWIGLS "nwigls"
#define PERTURBATION_RSLS "rsls"
#define PERTURBATION_RSffLS "rsffls"
#define PERTURBATION_RS "rsper"
#define PERTURBATION_RSFF "rsff"
#define PERTURBATION_IGIO "igio"
#define PERTURBATION_RSIO "rsio"
#define PERTURBATION_CP3 "cp3"
#define PERTURBATION_IG_OPTIMIZED "igoper"
#define PERTURBATION_IGLS_OPTIMIZED "igols"
#define PERTURBATION_NWIG "nwig"
#define PERTURBATION_NIIG "niig"
#define PERTURBATION_MPTLM "mptlm"
#define PERTURBATION_RESTART "restart"
#define PERTURBATION_RESTART_LS "restartls"
#define PERTURBATION_IG_SDST "sdstigo"
#define PERTURBATION_IGLS_SDST "sdstigols"




/* acceptance criteria*/
#define ACCEPTANCE_PROB "prob"
#define ACCEPTANCE_RS "rsacc"
#define ACCEPTANCE_KAR "karacc"
#define ACCEPTANCE_TMIIG "tmiigacc"
#define ACCEPTANCE_TEST "testacc"
#define ACCEPTANCE_SOA "soaacc"

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
char* problem_string;
emili::pfsp::PermutationFlowShop* loadProblem(char* t, PfspInstance i);

std::string info_pfsp()
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


emili::LocalSearch* prs::PfspBuilder::buildAlgo()
{
    prs::incrementTabLevel();
    emili::LocalSearch* ls = nullptr;
    if(tm.checkToken(CH6_LS))
    {
        printTab("CH6");
        emili::InitialSolution* in = retrieveComponent(COMPONENT_INITIAL_SOLUTION_GENERATOR).get<emili::InitialSolution>();
        emili::Termination* te = retrieveComponent(COMPONENT_TERMINATION_CRITERION).get<emili::Termination>();
        emili::Neighborhood* ne = retrieveComponent(COMPONENT_NEIGHBORHOOD).get<emili::Neighborhood>();
        emili::Neighborhood* ne2 = retrieveComponent(COMPONENT_NEIGHBORHOOD).get<emili::Neighborhood>();
        ls = new emili::pfsp::CH6(*in,*te,*ne,*ne2);
    }
    else if(tm.checkToken(RIS_LS))
    {
        printTab("RIS");
        emili::pfsp::PermutationFlowShop* instance =(emili::pfsp::PermutationFlowShop*) gp.getInstance();
        emili::InitialSolution* in = retrieveComponent(COMPONENT_INITIAL_SOLUTION_GENERATOR).get<emili::InitialSolution>();
        ls = new emili::pfsp::RIS(*instance,*in);
    }
    else if(tm.checkToken(NI_RIS_LS))
    {
        printTab("NoIdle RIS");
        emili::pfsp::PermutationFlowShop* instance =(emili::pfsp::PermutationFlowShop*) gp.getInstance();
        emili::InitialSolution* in = retrieveComponent(COMPONENT_INITIAL_SOLUTION_GENERATOR).get<emili::InitialSolution>();
        ls = new emili::pfsp::NoIdle_RIS(*instance,*in);
    }
    else if(tm.checkToken(NW_RIS_LS))
    {
        printTab("NoWait RIS");
        emili::pfsp::NWPFSP_MS* instance =(emili::pfsp::NWPFSP_MS*) gp.getInstance();
        emili::InitialSolution* in = retrieveComponent(COMPONENT_INITIAL_SOLUTION_GENERATOR).get<emili::InitialSolution>();
        ls = new emili::pfsp::NoWait_RIS(*instance,*in);
    }
    else if(tm.checkToken(RNW_RIS_LS))
    {
        printTab("Random NoWait RIS");
        emili::pfsp::NWPFSP_MS* instance =(emili::pfsp::NWPFSP_MS*) gp.getInstance();
        emili::InitialSolution* in = retrieveComponent(COMPONENT_INITIAL_SOLUTION_GENERATOR).get<emili::InitialSolution>();
        ls = new emili::pfsp::RandomNoWait_RIS(*instance,*in);
    }
    else if(tm.checkToken(SWP_INC_LS))
    {
        printTab("SwapInc local search");
        emili::InitialSolution* in = retrieveComponent(COMPONENT_INITIAL_SOLUTION_GENERATOR).get<emili::InitialSolution>();
        int r = 3;
        ls = new emili::pfsp::SwapIncLocalSearch(r,*in);
    }
    else if(tm.checkToken(STH_LS))
    {
        printTab("STH");
        emili::pfsp::SDSTFSP_MS* instance =(emili::pfsp::SDSTFSP_MS*) gp.getInstance();
        emili::InitialSolution* in = retrieveComponent(COMPONENT_INITIAL_SOLUTION_GENERATOR).get<emili::InitialSolution>();
        int b = tm.getInteger();
        printTabPlusOne("b",b);
        ls = new emili::pfsp::STH(b,*instance,*in);
    }
    else if(tm.checkToken(STHF_LS))
    {
        printTab("STH");
        emili::pfsp::SDSTFSP_MS* instance =(emili::pfsp::SDSTFSP_MS*) gp.getInstance();
        emili::InitialSolution* in = retrieveComponent(COMPONENT_INITIAL_SOLUTION_GENERATOR).get<emili::InitialSolution>();
        int n4 = instance->getNjobs()/4;
        n4 = n4==0?1:n4;
        int b  = n4 + emili::generateRandomNumber()%n4;
        //int b = tm.getInteger();
        printTabPlusOne("b",b);
        ls = new emili::pfsp::STH(b,*instance,*in);
    }

    prs::decrementTabLevel();
    return ls;

}

emili::Perturbation* prs::PfspBuilder:: buildPerturbation()
{
    prs::incrementTabLevel();
    emili::pfsp::PermutationFlowShop* instance =(emili::pfsp::PermutationFlowShop*) gp.getInstance();
    std::ostringstream oss;
    emili::Perturbation* per=nullptr;
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
    }else if(tm.checkToken(PERTURBATION_NWIG))
    {
        int nj = instance->getNjobs();
        int n = tm.getInteger();
        n = n<nj?n:nj-1;
        oss << "No wait optimized NEH destruct/construct perturbation. number of job erased: "<<n;
        printTab(oss.str().c_str());
        per = new emili::pfsp::NWIGPerturbation(n,*((emili::pfsp::NWPFSP_MS*)instance));
    }
    else if(tm.checkToken(PERTURBATION_NIIG))
        {
            int nj = instance->getNjobs();
            int n = tm.getInteger();
            n = n<nj?n:nj-1;
            oss << "No idle optimized NEH destruct/construct perturbation. number of job erased: "<<n;
            printTab(oss.str().c_str());
            per = new emili::pfsp::NoIdleIGper(n,*((emili::pfsp::NWPFSP_MS*)instance));
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
            PfspInstance pfs = instance->getInstance();
            pfs.setNbJob(pfs.getNbJob()-n);
            emili::pfsp::PermutationFlowShop * pfse = loadProblem(problem_string,pfs);
            gp.setInstance(pfse);
            emili::LocalSearch* ll = retrieveComponent(COMPONENT_ALGORITHM).get<emili::LocalSearch>();
            gp.setInstance(instance);
            //instances.push_back(pfse);
            per = new emili::pfsp::IgLsPerturbation(n,*instance,ll);
        } else {
             per = new emili::pfsp::IGPerturbation(1,*instance);
        }
    }
    else if(tm.checkToken(PERTURBATION_NWIGLS))
    {
        int nj = instance->getNjobs()-2;
        int k = tm.getInteger();
        int n = k<nj?k:nj-1;
        oss.str(""); oss  << "IG perturbation with local search applied on the partial solution. d = "<<n;
        printTab(oss.str().c_str());
        if(n > 0)
        {
            PfspInstance pfs = instance->getInstance();
            pfs.setNbJob(pfs.getNbJob()-n);
            emili::pfsp::PermutationFlowShop * pfse = loadProblem(problem_string,pfs);
            gp.setInstance(pfse);
            emili::LocalSearch* ll = retrieveComponent(COMPONENT_ALGORITHM).get<emili::LocalSearch>();
            gp.setInstance(instance);
            //instances.push_back(pfse);
            per = new emili::pfsp::NwIgLsPerturbation(n,*((emili::pfsp::NWPFSP_MS*)instance),ll);
        } else {
             per = new emili::pfsp::NWIGPerturbation(1,*((emili::pfsp::NWPFSP_MS*)instance));
        }
    }
    else if(tm.checkToken(PERTURBATION_IGLS_OPTIMIZED))
    {
        int nj = instance->getNjobs()-2;
        int k = tm.getInteger();
        int n = k<nj?k:nj-1;
        oss.str(""); oss  << "IG perturbation with local search applied on the partial solution. d = "<<n;
        printTab(oss.str().c_str());
        if(n > 0)
        {
            PfspInstance pfs = instance->getInstance();
            pfs.setNbJob(pfs.getNbJob()-n);
            emili::pfsp::PermutationFlowShop * pfse = loadProblem(problem_string,pfs);
            gp.setInstance(pfse);
            emili::LocalSearch* ll = retrieveComponent(COMPONENT_ALGORITHM).get<emili::LocalSearch>();
            gp.setInstance(instance);
            //instances.push_back(pfse);
            per = new emili::pfsp::IGOLsPerturbation(n,*instance,ll);
        } else {
             per = new emili::pfsp::IGOPerturbation(1,*instance);
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
            PfspInstance pfs = instance->getInstance();
            pfs.setNbJob(pfs.getNbJob()-n);
            emili::pfsp::PermutationFlowShop * pfse = loadProblem(problem_string,pfs);
            gp.setInstance(pfse);
            emili::LocalSearch* ll = retrieveComponent(COMPONENT_ALGORITHM).get<emili::LocalSearch>();
            gp.setInstance(instance);
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
            PfspInstance pfs = instance->getInstance();
            pfs.setNbJob(pfs.getNbJob()-n);
            emili::pfsp::PermutationFlowShop * pfse = loadProblem(problem_string,pfs);
            gp.setInstance(pfse);
            emili::LocalSearch* ll = retrieveComponent(COMPONENT_ALGORITHM).get<emili::LocalSearch>();
            gp.setInstance(instance);
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
    }
    else if(tm.checkToken(PERTURBATION_NRZ))
    {
        int nj = instance->getNjobs();
        int n = tm.getInteger();
        n = n<nj?n:nj-1;
        oss.str(""); oss  << "neh rz destruct/construct PERTURBATION. number of job erased: "<<n;
        printTab(oss.str().c_str());
        per = new emili::pfsp::NRZPerturbation(n,*instance);
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
    else if(tm.checkToken(PERTURBATION_NWTMIIG))
    {
        int nj = instance->getNjobs();
        int n = tm.getInteger();
        n = n<nj?n:nj-1;
        int tsize = tm.getInteger();
        oss.str(""); oss  << "No wait optimized TMIIG PERTURBATION. Number of job erased " << n << ". tabu list size " << tsize <<".\n\t";
        printTab(oss.str().c_str());
        per = new emili::pfsp::NWTMIIGPerturbation(n,*((emili::pfsp::NWPFSP_MS*)instance),tsize);
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
    else if(tm.checkToken(PERTURBATION_IG_OPTIMIZED))
    {
        int nj = instance->getNjobs();
        int n = tm.getInteger();
        n = n<nj?n:nj-1;

        oss.str(""); oss  << "IG perturbation with general optimization. d= " << n <<".\n\t";
        printTab(oss.str().c_str());
        per = new emili::pfsp::IGOPerturbation(n,*instance);
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
    else if(tm.checkToken(PERTURBATION_IG_SDST))
    {
        int nj = instance->getNjobs();
        int n = tm.getInteger();
        n = n<nj?n:nj-1;

        oss.str(""); oss  << "SDST IG perturbation with general optimization. d= " << n <<".\n\t";
        printTab(oss.str().c_str());
        per = new emili::pfsp::SDSTIGOPerturbation(n,*instance);
    }
    else if(tm.checkToken(PERTURBATION_IGLS_SDST))
    {
        int nj = instance->getNjobs()-2;
        int k = tm.getInteger();
        int n = k<nj?k:nj-1;
        oss.str(""); oss  << "SDST IG perturbation with local search applied on the partial solution. d = "<<n;
        printTab(oss.str().c_str());
        if(n > 0)
        {
            PfspInstance pfs = instance->getInstance();
            pfs.setNbJob(pfs.getNbJob()-n);
            emili::pfsp::PermutationFlowShop * pfse = loadProblem(problem_string,pfs);
            gp.setInstance(pfse);
            emili::LocalSearch* ll = retrieveComponent(COMPONENT_ALGORITHM).get<emili::LocalSearch>();
            gp.setInstance(instance);
            //instances.push_back(pfse);
            per = new emili::pfsp::SDSTIGOLsPerturbation(n,*instance,ll);
        } else {
             per = new emili::pfsp::SDSTIGOPerturbation(1,*instance);
        }
    }
    else if(tm.checkToken(PERTURBATION_RESTART))
    {
        int n = tm.getInteger();
        oss.str(""); oss  << "Restart perturbation, n =" << n << "";
        printTab(oss.str().c_str());
        emili::InitialSolution* init = retrieveComponent(COMPONENT_INITIAL_SOLUTION_GENERATOR).get<emili::InitialSolution>();
        per = new emili::pfsp::RestartPerturbation(n,init);
    }
    else if(tm.checkToken(PERTURBATION_RESTART_LS))
    {
        int n = tm.getInteger();
        oss.str(""); oss  << "Restart perturbation, n =" << n << "";
        printTab(oss.str().c_str());
        emili::InitialSolution* init = retrieveComponent(COMPONENT_INITIAL_SOLUTION_GENERATOR).get<emili::InitialSolution>();
        emili::LocalSearch* ll = retrieveComponent(COMPONENT_ALGORITHM).get<emili::LocalSearch>();
        per = new emili::pfsp::RestartPerturbation(n,init,ll);
    }
    else if(tm.checkToken(PERTURBATION_MPTLM))
        {

            int n = tm.getInteger();
            oss.str(""); oss  << "mPTLM inspired perturbation (1-alpha)*np = " << n << "";
            printTab(oss.str().c_str());
            per = new emili::pfsp::MPTLMPerturbation(n,*((emili::pfsp::NWPFSP_MS*)instance));
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

    prs::decrementTabLevel();
    return per;
}

emili::Acceptance* prs::PfspBuilder::buildAcceptance()
{
    prs::incrementTabLevel();
    emili::pfsp::PermutationFlowShop* instance =(emili::pfsp::PermutationFlowShop*) gp.getInstance();
    emili::Acceptance* acc = nullptr;
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

        temp = n*(temp/(nj*nm))/10.0;

        oss.str(""); oss  << "metropolis like Ruiz Stuetzle 2006 acceptance. temperature : " << temp;
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
    }

    prs::decrementTabLevel();
    return acc;
}


emili::TabuMemory* prs::PfspBuilder::buildTabuTenure()
{
    prs::incrementTabLevel();
    std::ostringstream oss;
    emili::TabuMemory* tmem = nullptr;
    //printTab(oss.str().c_str());
    if(tm.checkToken(TABU_MEMORY_MOVES))
    {
        oss.str(""); oss << "USING MOVES\n\t";
        int ts = tm.getInteger();
        oss << "Tabu tenure size " << ts;
        printTab(oss.str().c_str());
        tmem = new  emili::pfsp::PfspMovesMemory(ts);
    }
    else if(tm.checkToken(TABU_MEMORY_MOVES2))
    {
        oss.str(""); oss << "USING MOVES2\n\t";
        int ts = tm.getInteger();
        oss << "Tabu tenure size " << ts;
        printTab(oss.str().c_str());
        tmem = new  emili::pfsp::PfspMovesMemory2(ts);
    }
    else if(tm.checkToken(TABU_MEMORY_HASHES))
    {
        oss.str(""); oss << "USING HASHES\n\t";
        int ts = tm.getInteger();
        oss << "Tabu tenure size " << ts;
        printTab(oss.str().c_str());
        tmem = new  emili::pfsp::PfspTabuHashMemory(ts);
    }
    else if(tm.checkToken(TABU_MEMORY_SOLUTIONS))
    {
        oss.str(""); oss << "USING FULL SOLUtiON\n\t";        
        int ts = tm.getInteger();
        oss << "Tabu tenure size " << ts;
        printTab(oss.str().c_str());
        tmem = new  emili::pfsp::PfspFullSolutionMemory(ts);
    }
    else if(tm.checkToken(TABU_MEMORY_TSAB))
    {
        oss.str(""); oss << "USING TSAB\n\t";        
        int ts = tm.getInteger();
        oss << "Tabu tenure size " << ts;
        printTab(oss.str().c_str());
        tmem = new  emili::pfsp::TSABMemory(ts);
    }
    else if(tm.checkToken(TABU_MEMORY_TSAB_TEST))
    {
        oss.str(""); oss << "USING TSAB\n\t";        
        int ts = tm.getInteger();
        oss << "Tabu tenure size " << ts;
        printTab(oss.str().c_str());
        tmem = new  emili::pfsp::TSABtestMemory(ts);
    }
    else if(tm.checkToken(TABU_MEMORY_VALUE))
    {
        oss.str(""); oss << "USING VALUE\n\t" ;        
        int ts = tm.getInteger();
        oss << "Tabu tenure size " << ts;
        printTab(oss.str().c_str());
        tmem = new  emili::pfsp::PfspTabuValueMemory(ts);
    }

    prs::decrementTabLevel();
    return tmem;
}

emili::InitialSolution* prs::PfspBuilder::buildInitialSolution()
{
    prs::incrementTabLevel();
    std::ostringstream oss;
    emili::InitialSolution* init = nullptr;
    emili::pfsp::PermutationFlowShop* instance =(emili::pfsp::PermutationFlowShop*) gp.getInstance();

    if(tm.checkToken(INITIAL_RANDOM))
    {
        printTab("Random initial solution");
        init = new emili::pfsp::PfspRandomInitialSolution(*instance);
    }else if(tm.checkToken(INITIAL_RANDOM_ITERATED))
    {
        printTab("Iterated random initial solution");
        int n = tm.getInteger();
        oss.str("");oss << "number of random solutions generated: "<< n;
        printTabPlusOne(oss.str().c_str());
        init = new emili::pfsp::RandomInitialSolution(*instance,n);
    }
    else if(tm.checkToken(INITIAL_SLACK))
    {
        printTab("SLACK initial solution");
        init = new emili::pfsp::PfspSlackInitialSolution(*instance);
    }else if(tm.checkToken(INITIAL_WNSLACK))
    {
        printTab( "NEH WSLACK initial solution");
        //init = new testIS(instance);
        init = new emili::pfsp::PfspNEHwslackInitialSolution(*instance);
    }
    else if(tm.checkToken(INITIAL_LIT))
        {
            printTab( "Less idle times initial solution");
            //return new testIS(instance);
            init = new emili::pfsp::LITSolution(*instance);
        }
    else if(tm.checkToken(INITIAL_RZ))
        {
            printTab( "rz initial solution");
            //return new testIS(instance);
            init = new emili::pfsp::RZSolution(*instance);
        }
    else if(tm.checkToken(INITIAL_NRZ))
        {
            printTab( "neh rz initial solution");
            //return new testIS(instance);
            init = new emili::pfsp::NeRZSolution(*instance);
        }
    else if(tm.checkToken(INITIAL_NRZ2))
        {
            printTab( "neh rz initial solution without improvement phase");
            //return new testIS(*instance);
            init = new emili::pfsp::NeRZ2Solution(*instance);
        }
    else if(tm.checkToken(INITIAL_NRZ2FF))
        {
            printTab( "nehff rz initial solution without improvement phase");
            //return new testIS(*instance);
            init = new emili::pfsp::NfRZ2Solution(*instance);
        }
    else if(tm.checkToken(INITIAL_SRZ))
        {
            printTab( "srz initial solution generator");
            //return new testIS(*instance);
            init = new emili::pfsp::SRZSolution(*instance);
        }
    else if(tm.checkToken(INITIAL_LR))
        {
            int n = tm.getInteger();
            printTab("LR initial solution");
            oss.str(""); oss << "starting sequences "<<n;
            printTabPlusOne(oss.str().c_str());
            // testIS(*instance);
            init = new emili::pfsp::LRSolution(*instance,n);
        }
    else if(tm.checkToken(INITIAL_LR_NM))
        {            
            int n = instance->getNjobs()/instance->getNmachines();
            printTab("LR initial solution");
            oss.str(""); oss << "starting sequences "<<n;
            printTabPlusOne(oss.str().c_str());
            // testIS(*instance);
            init = new emili::pfsp::LRSolution(*instance,n);
        }
    else if(tm.checkToken(INITIAL_NLR))
        {
        int n = tm.getInteger();
        printTab("NLR initial solution");
        oss.str(""); oss << "starting sequences "<<n;
        printTabPlusOne(oss.str().c_str());
        init = new emili::pfsp::NLRSolution(*instance,n);
        }
    else if(tm.checkToken(INITIAL_MNEH))
        {
            printTab( "mneh initial solution");
            //return new testIS(instance);
            init = new emili::pfsp::MNEH(*instance);
        }
    else if(tm.checkToken(INITIAL_RMNEH))
        {
            printTab( "mneh initial solution");
            //return new testIS(instance);
            init = new emili::pfsp::RMNEH(*instance);
        }
    else if(tm.checkToken(INITIAL_NEH))
    {
        printTab( "NEH initial solution");
        //return new testIS(instance);
        init = new emili::pfsp::NEH(*instance);
    }
    else if(tm.checkToken(INITIAL_NEHRS))
    {
        printTab( "NEHRS (random restarts) initial solution");
        //return new testIS(instance);
        int iterations = tm.getInteger();
        oss.str("");oss<<"number of restarts: " << iterations;
        printTabPlusOne(oss.str().c_str());
        init = new emili::pfsp::NEHRS(*instance,iterations);
    }
    else if(tm.checkToken(INITIAL_NEHEDD))
    {
        printTab( "NEHedd initial solution");
        //return new testIS(instance);
        init = new emili::pfsp::NEHedd(*instance);
    }
    else if(tm.checkToken(INITIAL_NEHFF))
    {
        printTab( "NEHFF initial solution");
        //return new testIS(instance);
        init = new emili::pfsp::NEHff(*instance);
    }
    else if(tm.checkToken(INITIAL_NEHLS))
    {
        printTab( "NEHls initial solution");
        PfspInstance pfs = instance->getInstance();
        emili::pfsp::PermutationFlowShop * pfse = loadProblem(problem_string,pfs);
        gp.setInstance(pfse);
        emili::LocalSearch* ll = retrieveComponent(COMPONENT_ALGORITHM).get<emili::LocalSearch>();
        gp.setInstance(instance);
        init = new emili::pfsp::NEHls(*instance,ll);
    }
    else if(tm.checkToken(INITIAL_NEHEDDLS))
    {
        printTab( "NEHls initial solution");
        PfspInstance pfs = instance->getInstance();
        emili::pfsp::PermutationFlowShop * pfse = loadProblem(problem_string,pfs);
        gp.setInstance(pfse);
        emili::LocalSearch* ll = retrieveComponent(COMPONENT_ALGORITHM).get<emili::LocalSearch>();
        gp.setInstance(instance);
        init = new emili::pfsp::NEHeddLS(*instance,ll);
    }
    else if(tm.checkToken(INITIAL_FRB5))
    {
        printTab( "FRB5 initial solution");
        PfspInstance pfs = instance->getInstance();
        emili::pfsp::PermutationFlowShop * pfse = loadProblem(problem_string,pfs);
        emili::InitialSolution* in = new emili::pfsp::PfspRandomInitialSolution(*pfse);
        emili::Termination* term = new emili::LocalMinimaTermination();
        emili::Neighborhood* nei = new emili::pfsp::TaillardAcceleratedInsertNeighborhood(*pfse);
        emili::LocalSearch* ll = new emili::FirstImprovementSearch(*in,*term,*nei);        
        init = new emili::pfsp::NEHls(*instance,ll);
    }
    else if(tm.checkToken(INITIAL_CSFRB5))
    {
        printTab( "CSFRB5 initial solution");
        PfspInstance pfs = instance->getInstance();
        emili::pfsp::PermutationFlowShop * pfse = loadProblem(problem_string,pfs);
        emili::InitialSolution* in = new emili::pfsp::PfspRandomInitialSolution(*pfse);
        emili::Termination* term = new emili::LocalMinimaTermination();
        emili::Neighborhood* nei = new emili::pfsp::CSTaillardAcceleratedInsertNeighborhood(*pfse);
        emili::LocalSearch* ll = new emili::FirstImprovementSearch(*in,*term,*nei);
        init = new emili::pfsp::NEHls(*instance,ll);
    }
    else if(tm.checkToken(INITIAL_FRB5_GENERAL))
    {
        printTab( "FRB5 initial solution");
        PfspInstance pfs = instance->getInstance();
        emili::pfsp::PermutationFlowShop * pfse = loadProblem(problem_string,pfs);
        emili::InitialSolution* in = new emili::pfsp::PfspRandomInitialSolution(*pfse);
        emili::Termination* term = new emili::LocalMinimaTermination();
        emili::Neighborhood* nei = new emili::pfsp::PfspInsertNeighborhood(*pfse);
        emili::LocalSearch* ll = new emili::FirstImprovementSearch(*in,*term,*nei);
        init = new emili::pfsp::NEHls(*instance,ll);
    }
    else if(tm.checkToken(INITIAL_NEHFFLS))
    {
        printTab( "NEHffls initial solution");
        PfspInstance pfs =instance->getInstance();
        emili::pfsp::PermutationFlowShop * pfse = loadProblem(problem_string,pfs);
        gp.setInstance(pfse);
        emili::LocalSearch* ll = retrieveComponent(COMPONENT_ALGORITHM).get<emili::LocalSearch>();
        gp.setInstance(instance);
        init = new emili::pfsp::NEHffls(*instance,ll);
    }
    else if(tm.checkToken(INITIAL_BS))
    {
        printTab(" BS based initial solution");
        double a = tm.getDecimal();
        printTabPlusOne("a",a);
        double b = tm.getDecimal();
        printTabPlusOne("b",b);
        double c = tm.getDecimal();
        printTabPlusOne("c",c);
        double e = tm.getDecimal();
        printTabPlusOne("e",e);
        int gamma = tm.getInteger();
        printTabPlusOne("gamma",gamma);
        init = new emili::pfsp::BeamSearchHeuristic(*instance,gamma,a,b,c,e);
    }
    else if(tm.checkToken(INITIAL_BSNN))
    {
        printTab(" BS based initial solution");
        double a = tm.getDecimal();
        printTabPlusOne("a",a);
        double b = tm.getDecimal();
        printTabPlusOne("b",b);
        double c = tm.getDecimal();
        printTabPlusOne("c",c);
        double e = tm.getDecimal();
        printTabPlusOne("e",e);
        int gamma = instance->getNjobs()/10;
        if(gamma==0)
            gamma = instance->getNjobs();
        printTabPlusOne("gamma",gamma);        
        init = new emili::pfsp::BeamSearchHeuristic(*instance,gamma,a,b,c,e);
    }
    else if(tm.checkToken(INITIAL_BS2))
    {
        printTab(" BS based initial solution");
        double a = tm.getDecimal();
        printTabPlusOne("a",a);
        double b = tm.getDecimal();
        printTabPlusOne("b",b);
        double c = tm.getDecimal();
        printTabPlusOne("c",c);
        double e = tm.getDecimal();
        printTabPlusOne("e",e);
        int gamma = tm.getInteger();
        if(gamma==0 || gamma> instance->getNjobs())
            gamma = instance->getNjobs();
        printTabPlusOne("gamma",gamma);
        init = new emili::pfsp::BSheuristic(*instance,gamma,a,b,c,e);
    }
    else if(tm.checkToken(INITIAL_BS2N))
    {
        printTab(" BS based initial solution");
        double a = tm.getDecimal();
        printTabPlusOne("a",a);
        double b = tm.getDecimal();
        printTabPlusOne("b",b);
        double c = tm.getDecimal();
        printTabPlusOne("c",c);
        double e = tm.getDecimal();
        printTabPlusOne("e",e);
        int gamma = instance->getNjobs()*tm.getDecimal();
        if(gamma==0 || gamma> instance->getNjobs())
            gamma = instance->getNjobs();
        printTabPlusOne("gamma",gamma);
        init = new emili::pfsp::BSheuristic(*instance,gamma,a,b,c,e);
    }
    else if(tm.checkToken(INITIAL_BS2NF))
    {
        printTab(" BS based initial solution tuned by Irace");
        double a = 0.5662;
        printTabPlusOne("a",a);
        double b = 2.4326;
        printTabPlusOne("b",b);
        double c = 52.8383;
        printTabPlusOne("c",c);
        double e = 46.9122;
        printTabPlusOne("e",e);
        int gamma = instance->getNjobs()*0.1;
        if(gamma==0 || gamma> instance->getNjobs())
            gamma = instance->getNjobs();
        printTabPlusOne("gamma",gamma);
        init = new emili::pfsp::BSheuristic(*instance,gamma,a,b,c,e);
    }
    else if(tm.checkToken(INITIAL_FF))
    {
        printTab(" FF initial solution");
        double a = tm.getDecimal();
        printTabPlusOne("a",a);
        double b = tm.getDecimal();
        printTabPlusOne("b",b);
        int x = tm.getInteger();
        printTabPlusOne("x",x);
        init = new emili::pfsp::FFheuristic(*instance,x,a,b);
    }
    else if(tm.checkToken(INITIAL_FFN))
    {
        printTab(" FF initial solution");
        double a = tm.getDecimal();
        printTabPlusOne("a",a);
        double b = tm.getDecimal();
        printTabPlusOne("b",b);
        int x = instance->getNjobs();
        printTabPlusOne("x",x);
        init = new emili::pfsp::FFheuristic(*instance,x,a,b);
    }
    else if(tm.checkToken(INITIAL_BSCHO))
    {
        printTab(" BS based initial solution");
        double a = tm.getDecimal();
        printTabPlusOne("a",a);
        double b = tm.getDecimal();
        printTabPlusOne("b",b);
        double c = tm.getDecimal();
        printTabPlusOne("c",c);
        int gamma = tm.getInteger();
        if(gamma==0 || gamma> instance->getNjobs())
            gamma = instance->getNjobs();
        printTabPlusOne("x",gamma);
        init = new emili::pfsp::BSCH(*instance,gamma,a,b,c);
    }
    else if(tm.checkToken(INITIAL_BSCHR))
    {
        printTab(" BS based initial solution");
        double a = tm.getDecimal();
        printTabPlusOne("a",a);
        double b = tm.getDecimal();
        printTabPlusOne("b",b);
        double c = tm.getDecimal();
        printTabPlusOne("c",c);
        int gamma = instance->getNjobs()*tm.getDecimal();
        if(gamma==0 || gamma> instance->getNjobs())
            gamma = instance->getNjobs();
        printTabPlusOne("x",gamma);
        init = new emili::pfsp::BSCH(*instance,gamma,a,b,c);
    }
    else if(tm.checkToken(INITIAL_BSCH))
    {
        printTab(" BS based initial solution");
        double a = 9;
        printTabPlusOne("a",a);
        double b = 3;
        printTabPlusOne("b",b);
        double c = 7;
        printTabPlusOne("c",c);
        float r = tm.getDecimal();
        int gamma = instance->getNjobs()*r;
        if(gamma==0 || gamma > instance->getNjobs())
            gamma = instance->getNjobs();
        printTabPlusOne("x",gamma);
        init = new emili::pfsp::BSCH(*instance,gamma,a,b,c);
    }


    prs::decrementTabLevel();
    return init;
}

emili::Termination* prs::PfspBuilder::buildTermination()
{
    prs::incrementTabLevel();

    emili::Termination* term=nullptr;
    emili::pfsp::PermutationFlowShop* instance =(emili::pfsp::PermutationFlowShop*) gp.getInstance();
    if(tm.checkToken(TERMINATION_ITERA))
    {

        int ti = tm.getInteger();        
        printTab("Relaxed local minima termination");
        std::ostringstream oss;
        oss << "number of max iterations "<< ti;
        printTabPlusOne(oss.str().c_str());
        term =  new emili::pfsp::PfspTerminationIterations(ti);
    }
    else if(tm.checkToken(TERMINATION_SOA))
    {
        printTab("Max iteration number termination");
        int ti = instance->getNjobs();
         ti = 2*(ti-1);
         std::ostringstream oss;
         oss << "number of max iterations "<< ti;
         printTabPlusOne(oss.str().c_str());
        term =  new emili::pfsp::SOAtermination(ti);
    }
    if(tm.checkToken(TERMINATION_KAR))
    {
        int ti = instance->getNjobs();
        printTab("Kar termination");
        term =  new emili::pfsp::KarTermination(ti);
    }
    if(tm.checkToken(TERMINATION_MAXSTEPS_WITHNOIMPROV))
    {
        int n = instance->getNjobs();
        printTab("Termination that stops after n not improving steps");
        term = new emili::MaxStepsNoImprov(n);
    }


    prs::decrementTabLevel();
    return term;
}

emili::Neighborhood* prs::PfspBuilder::buildNeighborhood()
{
    prs::incrementTabLevel();
    emili::pfsp::PfspNeighborhood* neigh = nullptr;
    emili::pfsp::PermutationFlowShop* instance =(emili::pfsp::PermutationFlowShop*) gp.getInstance();
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
    else if(tm.checkToken(NEIGHBORHOOD_SDSTTA_INSERT))
    {
        printTab( "Taillard acceleration for Sequence dependent setup times");
        neigh = new emili::pfsp::SDSTTaillardAcceleratedInsertNeighborhood(*instance);
    }
    else if(tm.checkToken(NEIGHBORHOOD_KAR))
    {
        printTab( "KAR2016 Neighborhood");
        neigh = new emili::pfsp::KarNeighborhood(*instance);
    }
    else if(tm.checkToken(NEIGHBORHOOD_NW_INSERT))
    {
        printTab("No wait delta evaluation insert");
        neigh = new emili::pfsp::NoWaitAcceleratedInsertNeighborhood(*((emili::pfsp::NWPFSP_MS*)instance));
    }
    else if(tm.checkToken(NEIGHBORHOOD_NW_TWO_INSERT))
    {
        printTab("No wait delta evaluation insert");
        neigh = new emili::pfsp::NoWaitAcceleratedTwoInsertNeighborhood(*((emili::pfsp::NWPFSP_MS*)instance));
    }
    else if(tm.checkToken(NEIGHBORHOOD_NW_EXCHANGE))
    {
        printTab("No wait delta evaluation exchange");
        neigh = new emili::pfsp::NoWaitAcceleratedExchangeNeighborhood(*((emili::pfsp::NWPFSP_MS*)instance));
    }
    else if(tm.checkToken(NEIGHBORHOOD_NW_TRANSPOSE))
    {
        printTab("No wait delta evaluation transpose");
        neigh = new emili::pfsp::NoWaitAcceleratedTransposeNeighborhood(*((emili::pfsp::NWPFSP_MS*)instance));
    }
    else if(tm.checkToken(NEIGHBORHOOD_CS_INSERT))
    {
        printTab("Insert Neighborhood that returns only the best insertion");
        neigh = new emili::pfsp::CSInsertNeighborhood(*instance);
    }
    else if(tm.checkToken(NEIGHBORHOOD_SDST_CS_INSERT))
    {
        printTab("SDST Insert Neighborhood that returns only the best insertion");
        neigh = new emili::pfsp::SDSTCSInsertNeighborhood(*instance);
    }

    prs::decrementTabLevel();
    return neigh;
}
emili::Problem* prs::PfspBuilder::buildProblem()
{
    emili::pfsp::PermutationFlowShop* instance =(emili::pfsp::PermutationFlowShop*) gp.getInstance();
    return loadProblem(tm.nextToken(),instance->getInstance());
}

emili::pfsp::PermutationFlowShop* loadProblem(char* t, PfspInstance i)
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

emili::Problem* prs::PfspBuilder::openInstance()
{
    PfspInstance i;    
    problem_string = tm.nextToken();
    bool ok;
    std::string pro(problem_string);
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
         emili::pfsp::PermutationFlowShop* instance = loadProblem(problem_string, i);
         return instance;
     }

        std::cout << info_pfsp() << std::endl;        
        exit(-1);
}

bool isParsable(std::string &problem)
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

bool prs::PfspBuilder::isCompatibleWith(char *problem_definition)
{
    std::string s(problem_definition);
    return isParsable(s);
}

bool prs::PfspBuilder::canOpenInstance(char *problem_definition)
{
    std::string s(problem_definition);
    return isParsable(s);
}

extern "C" {
    prs::Builder* getBuilder(prs::GeneralParserE* ge)
    {
        return new prs::PfspBuilder(*ge,ge->getTokenManager());
    }
}
