#include "sa_qap_parser.h"

#include <cstdlib>
#include <cstdio>
#include <ctime>
#include <cstring>
#include <iostream>
#include <sstream>
#include <algorithm>
#include "../QAP/qapneighborhood.h"
#include "../QAP/qap.h"

/* QAP */
// #define QAP "QAP"

/* initial solution heuristics */
#define INITIAL_RANDOM "random"

/* permutation flowshop neighborhoods*/
/*#define NEIGHBORHOOD_INSERT "insert"
#define NEIGHBORHOOD_BACK_INSERT "binsert"
#define NEIGHBORHOOD_FORW_INSERT "finsert"
#define NEIGHBORHOOD_TWO_INSERT "tinsert"
#define NEIGHBORHOOD_TRANSPOSE "transpose"
#define NEIGHBORHOOD_XTRANSPOSE "xtranspose"*/
#define NEIGHBORHOOD_EXCHANGE "exchange"
#define BEST2OPT_EXCHANGE "best2opt"
#define FIRST2OPT_EXCHANGE "first2opt"
/*#define NEIGHBORHOOD_TA_INSERT "tainsert"
#define NEIGHBORHOOD_TAx_INSERT "txinsert"
#define NEIGHBORHOOD_ATAx_INSERT "atxinsert"
#define NEIGHBORHOOD_HATAx_INSERT "hatxinsert"
#define NEIGHBORHOOD_NITA_INSERT "ntainsert"*/


// char* problem_type;

//using namespace prs;

SAInitTemp* SAQAPParser::INITTEMP(prs::TokenManager& tm,
                                  emili::InitialSolution *initsol) {

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
        SAInitTemp* init_temp = new RandomWalkInitTemp(initsol, length);
        init_temp->set(value);
        return init_temp;
    } else if (tm.checkToken(RANDOMWALKAVGINITTEMP)) {
        int length = tm.getInteger();
        double value = tm.getDecimal();
        SAInitTemp* init_temp = new RandomWalkAvgInitTemp(initsol, length);
        init_temp->set(value);
        return init_temp;
    } else {
        std::cerr << "SAInitTemp expected, not found : " << std::endl;
        std::cerr << tm.peek() << std::endl;
        exit(1);
    }

}


SAAcceptance* SAQAPParser::ACCEPTANCE(prs::TokenManager& tm) {

    if (tm.checkToken(METROPOLIS)) {
        double it = tm.getDecimal();
        double ft = tm.getDecimal();
        return new SAMetropolisAcceptance(it, ft);
    } else if (tm.checkToken(BASICACC)) {
        return new SABasicAcceptance();
    } else if (tm.checkToken(GEOMACC)) {
        double ia = tm.getDecimal();
        double rf = tm.getDecimal();
        return new SAGeometricAcceptance(ia, rf);
    } else if (tm.checkToken(DETERMINISTICACC)) {
        double de = tm.getDecimal();
        return new SADeterministicAcceptance(de);
    } else if (tm.checkToken(LAHCACC)) {
        double te = tm.getInteger();
        return new LAHCAcceptance(te);
    } else {
        std::cerr << "SAAcceptance expected, not found : " << std::endl;
        std::cerr << tm.peek() << std::endl;
        exit(1);
    }

}


SATermination* SAQAPParser::TERMINATION(prs::TokenManager& tm) {

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
    } else {
        std::cerr << "SATermination expected, not found : " << std::endl;
        std::cerr << tm.peek() << std::endl;
        exit(1);
    }

}


SACooling* SAQAPParser::COOL(prs::TokenManager& tm,
                             SAInitTemp *it) {

    if (tm.checkToken(GEOM)) {
        float a = tm.getDecimal();
        float b = tm.getDecimal();
        return new GeomCooling(a,b, it);
    } else if (tm.checkToken(MARKOV)) {
        float a = tm.getDecimal();
        float b = tm.getDecimal();
        return new MarkovCooling(a,b, it);
    } else if (tm.checkToken(LOGCOOLING)) {
        float a = tm.getDecimal();
        float b = tm.getDecimal();
        return new LogCooling(a,b, it);
    } else if (tm.checkToken(CONSTCOOLING)) {
        float a = tm.getDecimal();
        float b = tm.getDecimal();
        return new ConstantCooling(a,b, it);
    } else if (tm.checkToken(LUNDYMEES)) {
        float a = tm.getDecimal();
        float b = tm.getDecimal();
        return new LundyMeesCooling(a,b, it);
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
    } else {
        std::cerr << "SACooling expected, not found : " << std::endl;
        std::cerr << tm.peek() << std::endl;
        exit(1);
    }

}


SATempLength* SAQAPParser::TEMPLENGTH(prs::TokenManager& tm,
                                       emili::Neighborhood* neigh) {

    if (tm.checkToken(CONSTTEMPLEN)) {
        int l = tm.getInteger();
        return new ConstantTempLength(l);
    } else if (tm.checkToken(NEIGHSIZETEMPLEN)) {
        float a = tm.getDecimal();
        return new NeighSizeTempLength(neigh, a);
    } else if (tm.checkToken(MAXACCEPTEDTEMPLEN)) {
        int a = tm.getInteger();
        return new MaxAcceptedTempLength(a);
    } else if (tm.checkToken(NOTEMPLEN)) {
        return new NoTempLength();
    } else {
        std::cerr << "SATempLength expected, not found : " << std::endl;
        std::cerr << tm.peek() << std::endl;
        exit(1);
    }

}


emili::InitialSolution* SAQAPParser::init(prs::TokenManager& tm)
{
    prs::incrementTabLevel();
    std::ostringstream oss;
    emili::InitialSolution* init;
    if(tm.checkToken(INITIAL_RANDOM))
    {
        prs::printTab("Random initial solution");
        init = new QAPRandomInitialSolution(*instance);
    }
    else
    {
        std::cerr<< "'" << *tm << "' -> ERROR an initial solution generator specification was expected! (random,slack)" << std::endl;

        std::cout << SAQAPParser::info() << std::endl;
        exit(-1);
    }
    prs::decrementTabLevel();
    return init;
}


QAPNeighborhood* SAQAPParser::neigh(prs::TokenManager& tm)
{
    prs::incrementTabLevel();
    QAPNeighborhood* neigh;
    if(tm.checkToken(NEIGHBORHOOD_EXCHANGE)) {
        prs::printTab( "Exchange neighborhood");
        neigh = new QAPExchangeNeighborhood(*instance);
    } else if (tm.checkToken(BEST2OPT_EXCHANGE)) {
        prs::printTab( "Best improvement 2-opt neighborhood");
        neigh = new QAPBest2optNeighborhood(*instance);
    } else if (tm.checkToken(FIRST2OPT_EXCHANGE)) {
        prs::printTab( "First improvement 2-opt neighborhood");
        neigh = new QAPFirst2optNeighborhood(*instance);
    } else {
        std::cerr<< "'" << *tm << "' -> ERROR a neighborhood specification was expected! " << std::endl;
        std::cout << SAQAPParser::info() << std::endl;
        exit(-1);
    }
    prs::decrementTabLevel();
    return neigh;
}


void SAQAPParser::problem(prs::TokenManager& tm) {
    instance = new qap::QAP(tm.tokenAt(1));
    return;
}

SAExploration* SAQAPParser::EXPLORATION(prs::TokenManager& tm,
                                        emili::Neighborhood* neigh,
                                        SAAcceptance *acc,
                                        SATermination *term) {

    if (tm.checkToken(SARANDOMEXPLORATION)) {
        return new SARandomExploration(neigh, acc, term);
    } else if (tm.checkToken(SASEQUENTIALEXPLORATION)) {
        return new SASequentialExploration(neigh, acc, term);
    } else {
        std::cerr << "SAExploration expected, not found : " << std::endl;
        std::cerr << tm.peek() << std::endl;
        exit(1);
    }

}

SATempRestart* SAQAPParser::TEMPRESTART(prs::TokenManager& tm,
                                        SAInitTemp *it) {

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
    } else if (tm.checkToken(SALOCALMINENHANCEDREHEAT)) {
        int   te = tm.getInteger();
        float va = tm.getDecimal();
        float ep = tm.getDecimal();
        return new SALocalMinEnhancedReheat(it, te, va, ep);
    } else {
        std::cerr << "SATempRestart expected, not found : " << std::endl;
        std::cerr << tm.peek() << std::endl;
        exit(1);
    }


}


emili::LocalSearch* SAQAPParser::buildAlgo(prs::TokenManager& tm) {

    problem(tm);
    emili::InitialSolution* initsol    = init(tm);
    emili::Neighborhood*    nei        = neigh(tm);
    SAInitTemp*      inittemp   = INITTEMP(tm, initsol);
    SAAcceptance*    acceptance = ACCEPTANCE(tm);
    SACooling*       cooling    = COOL(tm, inittemp);
    SATempRestart*   temprestart = TEMPRESTART(tm, inittemp);
    cooling->setTempRestart(temprestart);
    SATermination*     term       = TERMINATION(tm); // termin(tm);
    SATempLength*    templ      = TEMPLENGTH(tm, nei);
    cooling->setTempLength(templ);
    SAExploration* explo = EXPLORATION(tm, nei, acceptance, term);

    return new SimulatedAnnealing(initsol,
                                  inittemp,
                                  acceptance,
                                  cooling,
                                  temprestart,
                                  term,
                                  templ,
                                  explo,
                                  nei);

}

bool SAQAPParser::isParsable(string &problem)
{
    /*if(strcmp(problem.c_str(),PROBLEM_SDSTPFS_MS)==0)
    {*/
        return true;
    /*}
    else
    {
        return false;
    }*/
}

std::string SAQAPParser::info()
{
    ostringstream oss;
    oss << "Usage:\n\n";
    oss << "EMILI INSTANCE_FILE_PATH PFS_PROBLEM <LOCAL_SEARCH | ITERATED_LOCAL_SEARCH | TABU_SEARCH | VND_SEARCH> [rnds seed]\n\n";
    /*oss << "PROBLEM               = "<<PROBLEM_PFS_WT<< " " <<PROBLEM_PFS_WE<< " " <<PROBLEM_PFS_TCT<< " " <<PROBLEM_PFS_MS
              << " " <<PROBLEM_PFS_WCT<< " " <<PROBLEM_PFS_T<< " " <<PROBLEM_PFS_E<< " "<<PROBLEM_NWPFS_WT<< " " <<PROBLEM_NWPFS_WE
              << " " <<PROBLEM_NWPFS_TCT<< " " <<PROBLEM_NWPFS_MS<< " " <<PROBLEM_NWPFS_T<< " " <<PROBLEM_NWPFS_E
              << PROBLEM_NIPFS_MS <<" "<<PROBLEM_NIPFS_WT<< " " <<PROBLEM_NIPFS_WE<< " " <<PROBLEM_NIPFS_TCT<< " " <<PROBLEM_NIPFS_MS<< " "
              << " " <<PROBLEM_NIPFS_T<< " " <<PROBLEM_NIPFS_E <<"\n";
              */
    oss << "LOCAL_SEARCH          = SEARCH_TYPE INITIAL_SOLUTION TERMINATION NEIGHBORHOOD" <<"\n";
    oss << "ITERATED_LOCAL_SEARCH = ils LOCAL_SEARCH TERMINATION PERTUBATION ACCEPTANCE -it seconds" << "\n";
    oss << "TABU_SEARCH           = tabu INITIAL_SOLUTION TERMINATION NEIGHBORHOOD TABU_MEMORY" << "\n";
    oss << "VND_SEARCH            = vnd < first | best > INITIAL_SOLUTION TERMINATION NEIGHBORHOOD1 NEIGHBORHOOD2 ... NEIGHBORHOODn" << "\n";
    oss << "GVNS_SEARCH           = gvns INITIAL_SOLUTION PERTUBATION1 PERTUBATION2 -it seconds" << "\n";
    oss << "SEARCH_TYPE           = first | best | tabu | vnd | ils" << "\n";
    oss << "INITIAL_SOLUTION      = random | slack | nwslack | lit | rz | nrz | nrz2 | lr size(int)| nlr size(int) | mneh" <<"\n";
    oss << "TERMINATION           = true | time float | locmin | soater | iteration int | maxstep int" << "\n";
    //oss << "NEIGHBORHOOD          = transpose | exchange | insert | binsert | finsert | tinsert | "<< NEIGHBORHOOD_TA_INSERT << " | " << NEIGHBORHOOD_NITA_INSERT<<"\n";
    oss << "PERTUBATION           = igper int | testper | rndmv NEIGHBORHOOD #moves(int) | noper (int) | nrzper (int) | tmiigper (int) (int) | igls (int) LOCAL_SEARCH | rsls (int) LOCAL_SEARCH" << "\n";
    oss << "ACCEPTANCE            = soaacc float | testacc #swaps(int) | metropolis start_temperature(float) | always (intensify | diversify) | improve | sa_metropolis start_temp end_temp ratio | pmetro start_temp end_temp ratio frequence(int) | saacc start_temp end_temp ratio frequence(int) alpha ]0,1] | tmiigacc start_temperature(float) | implat number_of_non_improving_steps_accepted plateau_threshold" << std::endl;
    oss << "TABU_MEMORY           = move size(int) | hash size(int) | solution size(int) | tsabm size(int)" << "\n";
   // std::cout << " syntax->EMILI instancefile search_type intial_solution termination neighborhood" << std::endl;
    return oss.str();
}