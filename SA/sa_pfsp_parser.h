/*#ifndef SAPFSPPARSER_H
#define SAPFSPPARSER_H

#include "sa.h"

#include "sa_constants.h"
#include "sa_init_temp.h"
#include "sa_acceptance_criteria.h"
#include "sa_cooling.h"
#include "sa_termination_criteria.h"
#include "sa_templength.h"
#include "sa_exploration.h"
#include "sa_temperature_restart.h"


#include "../emilibase.h"
#include "../generalParser.h"
#include "../pfsp/permutationflowshop.h"

/* Permutation flowshop * /
#define PROBLEM_PFS_WT "PFSP_WT"
#define PROBLEM_PFS_WE "PFSP_WE"
#define PROBLEM_PFS_TCT "PFSP_TCT"
#define PROBLEM_PFS_T "PFSP_T"
#define PROBLEM_PFS_E "PFSP_E"
#define PROBLEM_PFS_WCT "PFSP_WCT"
#define PROBLEM_PFS_MS "PFSP_MS"

/* no wait permutation flowshop * /
#define PROBLEM_NWPFS_MS "NWPFSP_MS"
#define PROBLEM_NWPFS_WT "NWPFSP_WT"
#define PROBLEM_NWPFS_WE "NWPFSP_WE"
#define PROBLEM_NWPFS_TCT "NWPFSP_TCT"
#define PROBLEM_NWPFS_T "NWPFSP_T"
#define PROBLEM_NWPFS_E "NWPFSP_E"

/* no idle permutation flowshop* /
#define PROBLEM_NIPFS_MS "NIPFSP_MS"
#define PROBLEM_NIPFS_WT "NIPFSP_WT"
#define PROBLEM_NIPFS_WE "NIPFSP_WE"
#define PROBLEM_NIPFS_TCT "NIPFSP_TCT"
#define PROBLEM_NIPFS_T "NIPFSP_T"
#define PROBLEM_NIPFS_E "NIPFSP_E"

/* Sequence dependent setup times * /
#define PROBLEM_SDSTPFS_MS "SDSTPFS_MS"


/* initial solution heuristics * /
#define INITIAL_RANDOM "random"
#define INITIAL_SLACK "slack"
#define INITIAL_LIT "lit"
#define INITIAL_RZ "rz"
#define INITIAL_NRZ "nrz"
#define INITIAL_NRZ2 "nrz2"
#define INITIAL_LR "lr"
#define INITIAL_NLR "nlr"
#define INITIAL_MNEH "mneh"
#define INITIAL_WNSLACK "nwslack"

/* Termination criteria* /
#define TERMINATION_MAXSTEPS "maxstep"
#define TERMINATION_TIME "time"
#define TERMINATION_LOCMIN "locmin"
#define TERMINATION_ITERA "iteration"
#define TERMINATION_WTRUE "true"
#define TERMINATION_SOA "soater"

/* permutation flowshop neighborhoods* /
#define NEIGHBORHOOD_INSERT "insert"
#define NEIGHBORHOOD_BACK_INSERT "binsert"
#define NEIGHBORHOOD_FORW_INSERT "finsert"
#define NEIGHBORHOOD_TWO_INSERT "tinsert"
#define NEIGHBORHOOD_TRANSPOSE "transpose"
#define NEIGHBORHOOD_XTRANSPOSE "xtranspose"
#define NEIGHBORHOOD_EXCHANGE "exchange"
#define NEIGHBORHOOD_TA_INSERT "tainsert"
#define NEIGHBORHOOD_TAx_INSERT "txinsert"
#define NEIGHBORHOOD_ATAx_INSERT "atxinsert"
#define NEIGHBORHOOD_HATAx_INSERT "hatxinsert"
#define NEIGHBORHOOD_NITA_INSERT "ntainsert"

using namespace emili::sa;

class SAPFSPParser: public prs::AlgoBuilder {

protected:

    emili::pfsp::PermutationFlowShop* instantiateSAPFSPProblem(char* t, PfspInstance i);

    /**
     * an instance of PFS problem
     * /
    emili::pfsp::PermutationFlowShop* instance;

    /**
     * identify cooling scheme
     * @param  tm TokenManager
     * @return    SACooling object
     * /
    SACooling*       COOL(prs::TokenManager& tm,
                              SAInitTemp *it,
                              emili::Neighborhood *nei,
                              emili::pfsp::PermutationFlowShop *instance);

    /**
     * identify acceptance criterion
     * @param  tm TokenManager
     * @return    SAAcceptance object
     * /
    SAAcceptance*    ACCEPTANCE(prs::TokenManager& tm,
                                SAInitTemp *inittemp);

    /**
     * identify termination criterion
     * @param  tm TokenManager
     * @return    SATermination object
     * /
    SATermination*   TERMINATION(prs::TokenManager& tm,
                                 emili::Neighborhood *nei);

    /**
     * identify Neighborhood
     * @param  tm TokenManager
     * @return    Neighborhood object
     * /
    emili::Neighborhood*  NEIGH(prs::TokenManager& tm);

    /**
     * identify initial temperature
     * @param  tm      TokenManager
     * @param  initsol initial solution
     * @return         InitTemp object
     * /
    SAInitTemp*      INITTEMP(prs::TokenManager&      tm,
                              emili::InitialSolution* initsol,
                              emili::Neighborhood *nei,
                              emili::pfsp::PermutationFlowShop *instance);

    /**
     * identify initial solution builder
     * @param  tm TokenManager
     * @return    InitialSolution object
     * /
    emili::InitialSolution* INITSOL(prs::TokenManager& tm);


    SATempLength* TEMPLENGTH(prs::TokenManager& tm,
                             emili::Neighborhood* neigh,
                             emili::Problem* instance);

    emili::InitialSolution* init(prs::TokenManager& tm);
    emili::Termination* termin(prs::TokenManager& tm);
    emili::pfsp::PfspNeighborhood* neigh(prs::TokenManager& tm);

    SAExploration* EXPLORATION(prs::TokenManager& tm,
                                        emili::Neighborhood* neigh,
                                        SAAcceptance *acc,
                                        SACooling *cool,
                                        SATermination *term);

    SATempRestart *TEMPRESTART(prs::TokenManager& tm,
                               SAInitTemp *it,
                               emili::Neighborhood* neigh);

    /**
     * load the instance
     * @param tm token manager
     * /
    void problem(prs::TokenManager& tm);


public:

    /**
     * algorithm builder, according to grammar
     * @param  tm TokenManager
     * @return    assembled algorithm
     * /
    virtual emili::LocalSearch* buildAlgo(prs::TokenManager& tm);

    virtual bool isParsable(std::string& problem) ;

    virtual std::string info();

}; // class SAPFSPParser

#endif
*/