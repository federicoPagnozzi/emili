#include "SABuilder.h"


SAInitTemp* SABuilder::INITTEMP(prs::TokenManager& tm) {

    if (tm.checkToken(FIXEDINITTEMP)) {
        double      value     = tm.getDecimal();
        SAInitTemp* init_temp = new FixedInitTemp();
        init_temp->set(value);
        return init_temp;
    } else if (tm.checkToken(INITFROMSOL)) {
        emili::Solution* init_sol        = nullptr; // tm.getInitialSolution(); // ?
        double                  value    = tm.getDecimal();
        SAInitTemp* init_temp            = new InitTempFromSolution(init_sol);
        init_temp->set(value);
        return init_temp;
    } else {
        std::cerr << "SAInitTemp expected, not found : " << std::endl;
        std::cerr << tm.peek() << std::endl;
        exit(1);
    }

}


SAAcceptance* SABuilder::ACCEPTANCE(prs::TokenManager& tm) {

    if (tm.checkToken(METROPOLIS)) {
        double it = tm.getDecimal();
        double ft = tm.getDecimal();
        double dr = tm.getDecimal();
        int    ir = tm.getInteger();
        double al = tm.getDecimal();
        return new SAMetropolisAcceptance(it, ft, dr, ir, al);
    } else if (tm.checkToken(BASICACC)) {
        double it = tm.getDecimal();
        double ft = tm.getDecimal();
        double dr = tm.getDecimal();
        int    ir = tm.getInteger();
        double al = tm.getDecimal();
        return new SABasicAcceptance(it, ft, dr, ir, al);
    } else {
        std::cerr << "SAAcceptance expected, not found : " << std::endl;
        std::cerr << tm.peek() << std::endl;
        exit(1);
    }

}


SATermination* SABuilder::TERMINATION(prs::TokenManager& tm) {

    if (tm.checkToken(MAXBADITERS)) {
        int    mb = tm.getInteger();
        return new SAMaxBadIterTermination(mb);
    } else {
        std::cerr << "SATermination expected, not found : " << std::endl;
        std::cerr << tm.peek() << std::endl;
        exit(1);
    }

}


SACooling* SABuilder::COOL(prs::TokenManager& tm) {

    if (tm.checkToken(GEOM)) {
        float a = tm.getDecimal();
        float b = tm.getDecimal();
        return new GeomCooling(a,b);
    } else if (tm.checkToken(MARKOV)) {
        float a = tm.getDecimal();
        float b = tm.getDecimal();
        return new MarkovCooling(a,b);
    } else {
        std::cerr << "SACooling expected, not found : " << std::endl;
        std::cerr << tm.peek() << std::endl;
        exit(1);
    }

}


emili::Neighborhood*  SABuilder::NEIGH(prs::TokenManager& tm) {return nullptr;}

emili::InitialSolution*  SABuilder::INITSOL(prs::TokenManager& tm) {return nullptr;}


emili::LocalSearch* SABuilder::buildAlgo(prs::TokenManager& tm) {

    emili::InitialSolution*     initsol    = nullptr; // INITSOL(tm);
    emili::Neighborhood* neigh      = NEIGH(tm);
    SAInitTemp*          inittemp   = INITTEMP(tm);
    SAAcceptance*        acceptance = ACCEPTANCE(tm);
    SACooling*           cooling    = COOL(tm);
    SATermination*       term       = TERMINATION(tm);

    return new SimulatedAnnealing(initsol,
                                  inittemp,
                                  acceptance,
                                  cooling,
                                  term,
                                  neigh);

}
    