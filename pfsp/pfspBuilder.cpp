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

/* tabu tenure types */
#define TABU_MEMORY_MOVES "move"
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
#define INITIAL_NEHFFLS "nehffls"
#define INITIAL_RANDOM "random"
#define INITIAL_RANDOM_ITERATED "irandom"
#define INITIAL_SLACK "slack"
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
#define INITIAL_FRB5_GENERAL "gfrb5"

/* Termination criteria*/
#define TERMINATION_ITERA "iteration"
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

/*
 * END Neighborhoods
 */


/* permutation flowshop solution perturbations */
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
#define ACCEPTANCE_RS "rsacc"
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

    prs::decrementTabLevel();
    return ls;

}

emili::Perturbation* prs::PfspBuilder::buildPerturbation()
{
    prs::incrementTabLevel();
    emili::pfsp::PermutationFlowShop* istance =(emili::pfsp::PermutationFlowShop*) gp.getInstance();
    std::ostringstream oss;
    emili::Perturbation* per=nullptr;
    if(tm.checkToken(PERTURBATION_SOA) || tm.checkToken(PERTURBATION_SOA_LEGACY))
    {
        int nj = istance->getNjobs();
        int n = tm.getInteger();
        n = n<nj?n:nj-1;
        oss << "NEH destruct/construct perturbation which use objective function. number of job erased: "<<n;
        printTab(oss.str().c_str());
        per = new emili::pfsp::IGPerturbation(n,*istance);
    }else if(tm.checkToken(PERTURBATION_RS))
    {
        int nj = istance->getNjobs();
        int n = tm.getInteger();
        n = n<nj?n:nj-1;
        oss << "NEH destruct/construct perturbation. number of job erased: "<<n;
        printTab(oss.str().c_str());
        per = new emili::pfsp::RSPerturbation(n,*istance);
    }
    else if(tm.checkToken(PERTURBATION_RSFF))
        {
        int nj = istance->getNjobs();
        int n = tm.getInteger();
        n = n<nj?n:nj-1;
            oss << "NEH destruct/construct perturbation with tbff tie breaking. number of job erased: "<<n;
            printTab(oss.str().c_str());
            per = new emili::pfsp::RSffPerturbation(n,*istance);
        }
    else if(tm.checkToken(PERTURBATION_IGLS))
    {
        int nj = istance->getNjobs()-2;
        int k = tm.getInteger();
        int n = k<nj?k:nj-1;
        oss.str(""); oss  << "IG perturbation with local search applied on the partial solution. d = "<<n;
        printTab(oss.str().c_str());
        if(n > 0)
        {
            PfspInstance pfs = istance->getInstance();
            pfs.setNbJob(pfs.getNbJob()-n);
            emili::pfsp::PermutationFlowShop * pfse = loadProblem(problem_string,pfs);
            gp.setInstance(pfse);
            emili::LocalSearch* ll = retrieveComponent(COMPONENT_ALGORITHM).get<emili::LocalSearch>();
            gp.setInstance(istance);
            //istances.push_back(pfse);
            per = new emili::pfsp::IgLsPerturbation(n,*istance,ll);
        } else {
             per = new emili::pfsp::IGPerturbation(1,*istance);
        }
    }
    else if(tm.checkToken(PERTURBATION_RSLS))
    {
        int nj = istance->getNjobs()-2;
        int n = tm.getInteger();
        n = n<nj?n:nj-1;
        oss.str(""); oss  << "IG perturbation with local search applied on the partial solution. d = "<<n;
        printTab(oss.str().c_str());
        if(n > 0)
        {
            PfspInstance pfs = istance->getInstance();
            pfs.setNbJob(pfs.getNbJob()-n);
            emili::pfsp::PermutationFlowShop * pfse = loadProblem(problem_string,pfs);
            gp.setInstance(pfse);
            emili::LocalSearch* ll = retrieveComponent(COMPONENT_ALGORITHM).get<emili::LocalSearch>();
            gp.setInstance(istance);
            per = new emili::pfsp::RSLSPerturbation(n,*istance,ll);
        } else {
            per = new emili::pfsp::RSPerturbation(n,*istance);
        }
    }
    else if(tm.checkToken(PERTURBATION_RSffLS))
    {
        int nj = istance->getNjobs()-2;
        int n = tm.getInteger();
        n = n<nj?n:nj-1;
        oss.str(""); oss  << "IG perturbation with tbff tie breaking and local search applied on the partial solution. d = "<<n;
        printTab(oss.str().c_str());
        if(n > 0)
        {
            PfspInstance pfs = istance->getInstance();
            pfs.setNbJob(pfs.getNbJob()-n);
            emili::pfsp::PermutationFlowShop * pfse = loadProblem(problem_string,pfs);
            gp.setInstance(pfse);
            emili::LocalSearch* ll = retrieveComponent(COMPONENT_ALGORITHM).get<emili::LocalSearch>();
            gp.setInstance(istance);
            per = new emili::pfsp::RSffLSPerturbation(n,*istance,ll);
        } else {
            per = new emili::pfsp::RSffPerturbation(n,*istance);
        }
    }
    else if(tm.checkToken(PERTURBATION_TEST))
    {
        oss.str(""); oss<< "Random swap test perturbation.";
        printTab(oss.str().c_str());
        per = new emili::pfsp::PfspRandomSwapPertub(*istance);
    }
    else if(tm.checkToken(PERTURBATION_NRZ))
    {
        int nj = istance->getNjobs();
        int n = tm.getInteger();
        n = n<nj?n:nj-1;
        oss.str(""); oss  << "neh rz destruct/construct PERTURBATION. number of job erased: "<<n;
        printTab(oss.str().c_str());
        per = new emili::pfsp::NRZPerturbation(n,*istance);
    }
    else if(tm.checkToken(PERTURBATION_TMIIG))
    {
        int nj = istance->getNjobs();
        int n = tm.getInteger();
        n = n<nj?n:nj-1;
        int tsize = tm.getInteger();
        oss.str(""); oss  << "TMIIG PERTURBATION. Number of job erased " << n << ". tabu list size " << tsize <<".\n\t";
        printTab(oss.str().c_str());
        per = new emili::pfsp::TMIIGPerturbation(n,*istance,tsize);
    }
    else if(tm.checkToken(PERTURBATION_IGIO))
    {
        int nj = istance->getNjobs();
        int n = tm.getInteger();
        n = n<nj?n:nj-1;

        oss.str(""); oss  << "IG perturbation that inserts first the removed job with max sum of processing times. d= " << n <<".\n\t";
        printTab(oss.str().c_str());
        per = new emili::pfsp::IGIOPerturbation(n,*istance);
    }
    else if(tm.checkToken(PERTURBATION_RSIO))
    {
        int nj = istance->getNjobs();
        int n = tm.getInteger();
        n = n<nj?n:nj-1;

        oss.str(""); oss  << "IG perturbation that inserts first the removed job with max sum of processing times using taillard acceleration. d= " << n <<".\n\t";
        printTab(oss.str().c_str());
        per = new emili::pfsp::RSIOPerturbation(n,*istance);
    }
    else if(tm.checkToken(PERTURBATION_CP3))
    {
        int d = tm.getInteger();
        int omega = tm.getInteger();
        float pc = tm.getDecimal();
        oss.str(""); oss  << "Compound perturbation :  d= " << d << ", omega= " << omega << ",pc= "<< pc;
        printTab(oss.str().c_str());
        per = new emili::pfsp::CompoundPerturbation(*istance,omega,d,pc);
    }

    prs::decrementTabLevel();
    return per;
}

emili::Acceptance* prs::PfspBuilder::buildAcceptance()
{
    prs::incrementTabLevel();
    emili::pfsp::PermutationFlowShop* istance =(emili::pfsp::PermutationFlowShop*) gp.getInstance();
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
        acc = new  emili::pfsp::PfspTestAcceptance(*istance,n);
    }
    else  if(tm.checkToken(ACCEPTANCE_RS))
    {
        float n = tm.getDecimal();
        const std::vector < std::vector < long int > >& pm = istance->getProcessingTimesMatrix();
        int nj = istance->getNjobs();
        int nm = istance->getNmachines();

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
    else if(tm.checkToken(ACCEPTANCE_TMIIG))
    {
        float t0 =tm.getDecimal();
        float t=0.0f;
        int nj = istance->getNjobs();
        int nm = istance->getNmachines();
        const std::vector< std::vector < long > >& p = istance->getProcessingTimesMatrix();

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
    printTab(oss.str().c_str());

    if(tm.checkToken(TABU_MEMORY_MOVES))
    {
        oss.str(""); oss << "USING MOVES\n\t";
        printTab(oss.str().c_str());
        int ts = tm.getInteger();
        oss << "Tabu tenure size " << ts;
        tmem = new  emili::pfsp::PfspMovesMemory(ts);
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
        tmem = new  emili::pfsp::TSABMemory(ts);
    }
    else if(tm.checkToken(TABU_MEMORY_TSAB_TEST))
    {
        oss.str(""); oss << "USING TSAB\n\t";
        printTab(oss.str().c_str());
        int ts = tm.getInteger();
        oss << "Tabu tenure size " << ts;
        tmem = new  emili::pfsp::TSABtestMemory(ts);
    }
    else if(tm.checkToken(TABU_MEMORY_VALUE))
    {
        oss.str(""); oss << "USING VALUE\n\t" ;
        printTab(oss.str().c_str());
        int ts = tm.getInteger();
        oss << "Tabu tenure size " << ts;
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
        PfspInstance pfs = instance->getInstance();
        emili::pfsp::PermutationFlowShop * pfse = loadProblem(problem_string,pfs);
        gp.setInstance(pfse);
        emili::LocalSearch* ll = retrieveComponent(COMPONENT_ALGORITHM).get<emili::LocalSearch>();
        gp.setInstance(instance);
        init = new emili::pfsp::NEHls(*instance,ll);
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


    prs::decrementTabLevel();
    return init;
}

emili::Termination* prs::PfspBuilder::buildTermination()
{
    prs::incrementTabLevel();

    emili::Termination* term=nullptr;
    emili::pfsp::PermutationFlowShop* istance =(emili::pfsp::PermutationFlowShop*) gp.getInstance();
    if(tm.checkToken(TERMINATION_ITERA))
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
        int ti = istance->getNjobs();
         ti = 2*(ti-1);
        term =  new emili::pfsp::SOAtermination(ti);
    }

    prs::decrementTabLevel();
    return term;
}

emili::Neighborhood* prs::PfspBuilder::buildNeighborhood()
{
    prs::incrementTabLevel();
    emili::pfsp::PfspNeighborhood* neigh = nullptr;
    emili::pfsp::PermutationFlowShop* istance =(emili::pfsp::PermutationFlowShop*) gp.getInstance();
    if(tm.checkToken(NEIGHBORHOOD_INSERT))
    {
        printTab( "Insert Neighborhood");
        neigh = new emili::pfsp::PfspInsertNeighborhood(*istance);
    }
    else if(tm.checkToken(NEIGHBORHOOD_ADAPTIVE_INSERT))
    {
        printTab( "Insert Neighborhood");
        neigh = new emili::pfsp::PfspInsertNeighborhood(*istance);
    }
    else  if(tm.checkToken(NEIGHBORHOOD_FORW_INSERT))
    {
        printTab( "Forward insert Neighborhood");
        neigh = new emili::pfsp::PfspForwardInsertNeighborhood(*istance);
    }
    else  if(tm.checkToken(NEIGHBORHOOD_BACK_INSERT))
    {
        printTab( "Backward Insert Neighborhood");
        neigh = new emili::pfsp::PfspBackwardInsertNeighborhood(*istance);
    }
    else if(tm.checkToken(NEIGHBORHOOD_EXCHANGE))
    {
        printTab( "Exchange neighborhood");
        neigh = new emili::pfsp::PfspExchangeNeighborhood(*istance);
    }
    else if(tm.checkToken(NEIGHBORHOOD_ATX_EXCHANGE))
    {
        printTab( "Exchange neighborhood");
        neigh = new emili::pfsp::AxtExchange(*istance);
    }
    else if(tm.checkToken(NEIGHBORHOOD_OPT_EXCHANGE))
    {
        printTab( "Optimized Exchange neighborhood");
        neigh = new emili::pfsp::OptExchange(*istance);
    }
    else if(tm.checkToken(NEIGHBORHOOD_HATX_EXCHANGE))
       {
           printTab( "Exchange neighborhood with speedup");
           neigh = new emili::pfsp::HaxtExchange(*istance);
       }
       else if(tm.checkToken(NEIGHBORHOOD_EATX_EXCHANGE))
       {
           printTab( "Exchange neighborhood with speedup");
           neigh = new emili::pfsp::EaxtExchange(*istance);
       }
    else if(tm.checkToken(NEIGHBORHOOD_TRANSPOSE))
    {
        printTab( "Transpose neighborhood");
        neigh = new emili::pfsp::PfspTransposeNeighborhood(*istance);
    }
    else if(tm.checkToken(NEIGHBORHOOD_TWO_INSERT))
    {
        printTab( "Two insert neighborhood");
        neigh = new emili::pfsp::PfspTwoInsertNeighborhood(*istance);
    }
    else if(tm.checkToken(NEIGHBORHOOD_XTRANSPOSE))
    {
        printTab( "XTranspose neighborhood");
        neigh = new emili::pfsp::XTransposeNeighborhood(*istance);
    }
    else if(tm.checkToken(NEIGHBORHOOD_TA_INSERT))
    {
        printTab( "Insert with Taillard Acceleration");
        neigh = new emili::pfsp::TaillardAcceleratedInsertNeighborhood(*istance);
    }
    else if(tm.checkToken(NEIGHBORHOOD_FSTA_INSERT))
    {
        printTab( "Insert with Taillard Acceleration that updates the base solution after each improvement");
        neigh = new emili::pfsp::FSTaillardAcceleratedInsertNeighborhood(*istance);
    }
    else if(tm.checkToken(NEIGHBORHOOD_CSTA_INSERT))
    {
        printTab( "Insert with Taillard Acceleration that evaluates all the possible insertion points");
        neigh = new emili::pfsp::CSTaillardAcceleratedInsertNeighborhood(*istance);
    }
    else if(tm.checkToken(NEIGHBORHOOD_TAx_INSERT))
    {
        printTab( "Insert with Taillard Acceleration(Experimental)");
        neigh = new emili::pfsp::TAxInsertNeighborhood(*istance);
    }
    else if(tm.checkToken(NEIGHBORHOOD_OPT_INSERT))
    {
        printTab( "Delta Evaluation Insert for Weighted Tardiness with tail improvement");
        neigh = new emili::pfsp::OptInsert(*istance);
    }
    else if(tm.checkToken(NEIGHBORHOOD_ATAx_INSERT))
    {
        printTab( "Atx Delta Evaluation Insert for Weighted Tardiness");
        neigh = new emili::pfsp::AtxNeighborhood(*istance);
    }
    else if(tm.checkToken(NEIGHBORHOOD_HATAx_INSERT))
    {
        printTab( "Approximated Insert with Taillard Acceleration for Weighted Tardiness with no threshold");
        neigh = new emili::pfsp::HeavilyApproximatedTaillardAcceleratedInsertNeighborhood(*istance);
    }
    else if(tm.checkToken(NEIGHBORHOOD_NATAx_INSERT))
    {
        printTab( "Approximated Insert for Weighted Tardiness with 1 level approximation");
        neigh = new emili::pfsp::NatxNeighborhood(*istance);
    }
    else if(tm.checkToken(NEIGHBORHOOD_NATA2x_INSERT))
    {
        printTab( "Improved Approximated Insert for Weighted Tardiness with 1 level approximation and online tuned threshold");
        neigh = new emili::pfsp::Natx2Neighborhood(*istance);
    }
    else if(tm.checkToken(NEIGHBORHOOD_EATAx_INSERT))
    {
        printTab( "Approximated Insert for Weighted Tardiness with 2 levels of approximation");
        neigh = new emili::pfsp::EatxNeighborhood(*istance);
    }
    else if(tm.checkToken(NEIGHBORHOOD_THATAx_INSERT))
    {
        printTab( "Approximated Insert for Weighted Tardiness with 3 levels of approximation");
        neigh = new emili::pfsp::ThatxNeighborhood(*istance);
    }
    else if(tm.checkToken(NEIGHBORHOOD_PATAx_INSERT))
    {
        printTab( "Approximated Insert for Weighted Tardiness with 5 levels of approximation");
        neigh = new emili::pfsp::PatxNeighborhood(*istance);
    }
    else if(tm.checkToken(NEIGHBORHOOD_SATAx_INSERT))
    {
        printTab( "Approximated Insert for Weighted Tardiness with 6 levels of approximation");
        neigh = new emili::pfsp::SatxNeighborhood(*istance);
    }
    else if(tm.checkToken(NEIGHBORHOOD_FATAx_INSERT))
    {
        printTab( "Approximated Insert for Weighted Tardiness with 4 levels of approximation");
        neigh = new emili::pfsp::FatxNeighborhood(*istance);
    }
    else if(tm.checkToken(NEIGHBORHOOD_TATAx_INSERT))
    {
        printTab( "Approximated Insert for Weighted Tardiness with settable threshold");
        float start_level = tm.getDecimal();       
        neigh = new emili::pfsp::TatxNeighborhood(start_level,*istance);
    }
    else if(tm.checkToken(NEIGHBORHOOD_NITA_INSERT))
    {
        printTab( "Insert with Taillard Acceleration for no idle make span ");
        neigh = new emili::pfsp::NoIdleAcceleratedInsertNeighborhood(*istance);
    }
    else if(tm.checkToken(NEIGHBORHOOD_NATA_TCT_INSERT))
    {
        printTab( "Improved Approximated Insert for Total Completion Times with 1 level approximation and online tuned threshold");
        neigh = new emili::pfsp::NatxTCTNeighborhood(*istance);
    }
    else if(tm.checkToken(NEIGHBORHOOD_RZ_TCT_INSERT))
    {
        printTab( "iRZ neighborhood see PanRui2012");
        neigh = new emili::pfsp::NrzTCTNeighborhood(*istance);
    }
    else if(tm.checkToken(NEIGHBORHOOD_NATA_TT_INSERT))
    {
        printTab( "Improved Approximated Insert for Total Tardiness with 1 level approximation and online tuned threshold");
        neigh = new emili::pfsp::NatxTTNeighborhood(*istance);
    }
    prs::decrementTabLevel();
    return neigh;
}
emili::Problem* prs::PfspBuilder::buildProblem()
{
    emili::pfsp::PermutationFlowShop* istance =(emili::pfsp::PermutationFlowShop*) gp.getInstance();
    return loadProblem(tm.nextToken(),istance->getInstance());
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
         emili::pfsp::PermutationFlowShop* istance = loadProblem(problem_string, i);
         return istance;
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