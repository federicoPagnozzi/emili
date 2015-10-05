#ifndef SAQAPPARSER_H
#define SAQAPPARSER_H

#include "sa.h"

#include "sa_constants.h"
#include "sa_init_temp.h"
#include "sa_acceptance_criteria.h"
#include "sa_cooling.h"
#include "sa_termination_criteria.h"
#include "sa_templength.h"
#include "sa_exploration.h"


#include "../emilibase.h"
#include "../generalParser.h"
#include "../QAP/qapneighborhood.h"
#include "../QAP/qap.h"


class SAQAPParser: public prs::AlgoBuilder {

protected:

    /**
     * an instance of QAP problem
     */
    qap::QAP* instance;

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


    SATempLength* TEMPLENGTH(prs::TokenManager& tm,
                                       emili::Neighborhood* neigh);

    emili::InitialSolution* init(prs::TokenManager& tm);
    QAPNeighborhood* neigh(prs::TokenManager& tm);

    SAExploration* EXPLORATION(prs::TokenManager& tm,
                                        emili::Neighborhood* neigh,
                                        SAAcceptance *acc,
                                        SATermination *term);

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

}; // class SAQAPParser

#endif
