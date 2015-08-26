#ifndef SAPFSPPARSER_H
#define SAPFSPPARSER_H

#include "sa.h"

#include "sa_constants.h"
#include "sa_init_temp.h"
#include "sa_acceptance_criteria.h"
#include "sa_cooling.h"
#include "sa_termination_criteria.h"
#include "sa_neighborhood.h"

#include "../emilibase.h"
#include "../generalParser.h"
#include "../permutationflowshop.h"


class SAPFSPParser: public prs::AlgoBuilder {

protected:

    /**
     * an instance of PFS problem
     */
    emili::pfsp::PermutationFlowShop* instance;

    /**
     * identify cooling scheme
     * @param  tm TokenManager
     * @return    SACooling object
     */
    SACooling*       COOL(prs::TokenManager& tm);

    /**
     * identify acceptance criterion
     * @param  tm TokenManager
     * @return    SAAcceptance object
     */
    SAAcceptance*    ACCEPTANCE(prs::TokenManager& tm);

    /**
     * identify termination criterion
     * @param  tm TokenManager
     * @return    SATermination object
     */
    SATermination*   TERMINATION(prs::TokenManager& tm);

    /**
     * identify Neighborhood
     * @param  tm TokenManager
     * @return    Neighborhood object
     */
    emili::Neighborhood*  NEIGH(prs::TokenManager& tm);

    /**
     * identify initial temperature
     * @param  tm      TokenManager
     * @param  initsol initial solution
     * @return         InitTemp object
     */
    SAInitTemp*      INITTEMP(prs::TokenManager&      tm,
    	                      emili::InitialSolution* initsol);

    /**
     * identify initial solution builder
     * @param  tm TokenManager
     * @return    InitialSolution object
     */
    emili::InitialSolution* INITSOL(prs::TokenManager& tm);

    emili::InitialSolution* init(prs::TokenManager& tm);
    emili::Termination* termin(prs::TokenManager& tm);
    emili::pfsp::PfspNeighborhood* neigh(prs::TokenManager& tm);

    /**
     * load the instance
     * @param tm token manager
     */
    void problem(prs::TokenManager& tm);


public:

    /**
     * algorithm builder, according to grammar
     * @param  tm TokenManager
     * @return    assembled algorithm
     */
    virtual emili::LocalSearch* buildAlgo(prs::TokenManager& tm);

    virtual bool isParsable(std::string& problem) ;

    virtual std::string info();

}; // class SAPFSPParser

#endif
