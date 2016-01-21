#ifndef SABUILDER_H
#define SABUILDER_H

#include "../emilibase.h"
#include "../generalParser.h"

#include "sa_constants.h"
#include "sa_acceptance_criteria.h"
#include "sa_cooling.h"
#include "sa_termination_criteria.h"
#include "sa_neighborhood.h"
#include "sa_init_temp.h"

#include "sa.h"

class SABuilder: public prs::AlgoBuilder {

protected:

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
     * @return    SANeighborhood object
     */
    emili::Neighborhood*  NEIGH(prs::TokenManager& tm);

    /**
     * identify initial temperature
     * @param  tm TokenManager
     * @return    SAInitTemp object
     */
    SAInitTemp*      INITTEMP(prs::TokenManager& tm);

    /**
     * identify initial solution builder
     * @param  tm TokenManager
     * @return    InitialSolution object
     */
    emili::InitialSolution* INITSOL(prs::TokenManager& tm);

public:

    /**
     * algorithm builder, according to grammar
     * @param  tm TokenManager
     * @return    assembled algorithm
     */
    virtual emili::LocalSearch* buildAlgo(prs::TokenManager& tm);

}; // class SABuilder

#endif
