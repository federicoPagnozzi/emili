/*#include "sa_pfsp_parser.h"

#include <cstdlib>
#include <cstdio>
#include <ctime>
#include <cstring>
#include <iostream>
#include <sstream>
#include <algorithm>
#include "../pfsp/pfspinstance.h"
#include "../pfsp/permutationflowshop.h"


//using namespace prs;

emili::pfsp::PermutationFlowShop* SAPFSPParser::instantiateSAPFSPProblem(char* t, PfspInstance i)
{
    emili::pfsp::PermutationFlowShop* prob;
    if(strcmp(t,PROBLEM_PFS_WT)==0)
    {

        prs::printTab("Permutation Flow Shop Weighted Tardiness");
        prob = new emili::pfsp::PFSP_WT(i);
    }else if(strcmp(t,PROBLEM_NWPFS_MS)==0)
    {
        prs::printTab("No Wait Permutation Flow Shop Make Span");
        prob = new emili::pfsp::NWPFSP_MS(i);
    }else if(strcmp(t,PROBLEM_PFS_E)==0)
    {
        prs::printTab("Permutation Flow Shop Earliness");
        prob = new emili::pfsp::PFSP_E(i);
    }else if(strcmp(t,PROBLEM_PFS_WE)==0)
    {
        prs::printTab("Permutation Flow Shop Weighted Earliness");
        prob = new emili::pfsp::PFSP_WE(i);
    }else if(strcmp(t,PROBLEM_PFS_T)==0)
    {
        prs::printTab("Permutation Flow Shop Tardiness");
        prob = new emili::pfsp::PFSP_T(i);
    }else if(strcmp(t,PROBLEM_PFS_MS)==0)
    {
        prs::printTab("Permutation Flow Shop Make Span");
        prob = new emili::pfsp::PFSP_MS(i);
    }
    else if(strcmp(t,PROBLEM_NIPFS_MS)==0)
            {
                prs::printTab("No Idle Permutation Flow Shop Make Span" );
                prob = new emili::pfsp::NI_A_PFSP_MS(i);
            }
    else if(strcmp(t,PROBLEM_SDSTPFS_MS)==0)
    {
        prs::printTab("Sequence dependent setup times Make Span");
        prob = new emili::pfsp::SDSTFSP_MS(i);
    }
    else
    {
        std::cerr<< "'" << t << "' -> ERROR a problem was expected! " << std::endl;
        prs::info();
    exit(-1);
    }
    return prob;
}

SAInitTemp* SAPFSPParser::INITTEMP(prs::TokenManager& tm,
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
    } else if (tm.checkToken(OSMANPOTTSINITTEMP)) {
        float dc = tm.getDecimal();
        float tf = tm.getDecimal();
        float coeff = tm.getDecimal();
        SAInitTemp* init_temp = new OsmanPottsInitTemp(initsol, nei, instance, dc, tf);
        init_temp->set(coeff);
        return init_temp;
    } else {
        std::cerr << "SAInitTemp expected, not found : " << std::endl;
        std::cerr << tm.peek() << std::endl;
        exit(1);
    }

}


SAAcceptance* SAPFSPParser::ACCEPTANCE(prs::TokenManager& tm,
                                      SAInitTemp *inittemp) {

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


SATermination* SAPFSPParser::TERMINATION(prs::TokenManager& tm,
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
    } else {
        std::cerr << "SATermination expected, not found : " << std::endl;
        std::cerr << tm.peek() << std::endl;
        exit(1);
    }


}


SACooling* SAPFSPParser::COOL(prs::TokenManager& tm,
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
    } else  {
        std::cerr << "SACooling expected, not found : " << std::endl;
        std::cerr << tm.peek() << std::endl;
        exit(1);
    }

}


SATempLength* SAPFSPParser::TEMPLENGTH(prs::TokenManager& tm,
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


emili::InitialSolution* SAPFSPParser::init(prs::TokenManager& tm)
{
    prs::incrementTabLevel();
    std::ostringstream oss;
    emili::InitialSolution* init;
    if(tm.checkToken(INITIAL_RANDOM))
    {
        prs::printTab("Random initial solution");
        init = new emili::pfsp::PfspRandomInitialSolution(*instance);
    }
    else if(tm.checkToken(INITIAL_SLACK))
    {
        prs::printTab("SLACK initial solution");
        init = new emili::pfsp::PfspSlackInitialSolution(*instance);
    }else if(tm.checkToken(INITIAL_WNSLACK))
    {
        prs::printTab( "NEH WSLACK initial solution");
        //init = new testIS(instance);
        init = new emili::pfsp::PfspNEHwslackInitialSolution(*instance);
    }
    else if(tm.checkToken(INITIAL_LIT))
        {
            prs::printTab( "Less idle times initial solution");
            //return new testIS(instance);
            init = new emili::pfsp::LITSolution(*instance);
        }
    else if(tm.checkToken(INITIAL_RZ))
        {
            prs::printTab( "rz initial solution");
            //return new testIS(instance);
            init = new emili::pfsp::RZSolution(*instance);
        }
    else if(tm.checkToken(INITIAL_NRZ))
        {
            prs::printTab( "neh rz initial solution");
            //return new testIS(instance);
            init = new emili::pfsp::NeRZSolution(*instance);
        }
    else if(tm.checkToken(INITIAL_NRZ2))
        {
            prs::printTab( "neh rz initial solution without improvement phase");
            //return new testIS(*instance);
            init = new emili::pfsp::NeRZ2Solution(*instance);
        }
    else if(tm.checkToken(INITIAL_LR))
        {
            int n = tm.getInteger();
            oss.str(""); oss << "LR initial solution with "<<n<<" starting sequences";
            prs::printTab(oss.str().c_str());
            // testIS(*instance);
            init = new emili::pfsp::LRSolution(*instance,n);
        }
    else if(tm.checkToken(INITIAL_NLR))
        {
        int n = tm.getInteger();
        oss.str(""); oss << "NLR initial solution with "<<n<<" starting sequences";
        //return new testIS(*instance);printTab(oss.str().c_str());
        prs::printTab(oss.str().c_str());
        init = new emili::pfsp::NLRSolution(*instance,n);
        }
    else if(tm.checkToken(INITIAL_MNEH))
        {
            prs::printTab( "mneh initial solution");
            //return new testIS(instance);
            init = new emili::pfsp::MNEH(*instance);
        }
    else
    {
        std::cerr<< "'" << *tm << "' -> ERROR a initial solution generator specification was expected! (random,slack)" << std::endl;

        std::cout << SAPFSPParser::info() << std::endl;
        exit(-1);
    }
    prs::decrementTabLevel();
    return init;
}



emili::Termination* SAPFSPParser::termin(prs::TokenManager& tm)
{
    prs::incrementTabLevel();
    emili::Termination* term;
    if(tm.checkToken(TERMINATION_LOCMIN))
    {
        prs::printTab("Local minima termination");
        term = new emili::LocalMinimaTermination();
    }
    else if(tm.checkToken(TERMINATION_WTRUE))
    {
        prs::printTab("While true termination");
        term = new emili::WhileTrueTermination();
    }
    else if(tm.checkToken(TERMINATION_ITERA))
    {

        int ti = tm.getInteger();
        std::ostringstream oss;
        oss << "Relaxed local minima termination. number of max iterations "<< ti;
        prs::printTab(oss.str().c_str());
        term =  new emili::pfsp::PfspTerminationIterations(ti);
    }
    else if(tm.checkToken(TERMINATION_SOA))
    {
        prs::printTab("Max iteration number termination");
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
        prs::printTab(oss.str().c_str());
        term =  new emili::TimedTermination(time);
    }
    else if(tm.checkToken(TERMINATION_MAXSTEPS))
    {
        int steps = tm.getInteger();
        std::ostringstream oss;
        oss << "Max Steps termination. # steps: "<< steps;
        prs::printTab(oss.str().c_str());
        term = new emili::MaxStepsTermination(steps);
    }
    else
    {
        std::cerr<< "'" << *tm << "' -> ERROR a termination criteria specification was expected! " << std::endl;
        std::cout << SAPFSPParser::info() << std::endl;
        exit(-1);
    }
    prs::decrementTabLevel();
    return term;
}

emili::pfsp::PfspNeighborhood* SAPFSPParser::neigh(prs::TokenManager& tm)
{
    prs::incrementTabLevel();
    emili::pfsp::PfspNeighborhood* neigh;
    if(tm.checkToken(NEIGHBORHOOD_INSERT))
    {
        prs::printTab( "Insert Neighborhood");
        neigh = new emili::pfsp::PfspInsertNeighborhood(*instance);
    }
    else  if(tm.checkToken(NEIGHBORHOOD_FORW_INSERT))
    {
        prs::printTab( "Forward insert Neighborhood");
        neigh = new emili::pfsp::PfspForwardInsertNeighborhood(*instance);
    }
    else  if(tm.checkToken(NEIGHBORHOOD_BACK_INSERT))
    {
        prs::printTab( "Backward Insert Neighborhood");
        neigh = new emili::pfsp::PfspBackwardInsertNeighborhood(*instance);
    }
    else if(tm.checkToken(NEIGHBORHOOD_EXCHANGE))
    {
        prs::printTab( "Exchange neighborhood");
        neigh = new emili::pfsp::PfspExchangeNeighborhood(*instance);
    }
    else if(tm.checkToken(NEIGHBORHOOD_TRANSPOSE))
    {
        prs::printTab( "Transpose neighborhood");
        neigh = new emili::pfsp::PfspTransposeNeighborhood(*instance);
    }
    else if(tm.checkToken(NEIGHBORHOOD_TWO_INSERT))
    {
        prs::printTab( "Two insert neighborhood");
        neigh = new emili::pfsp::PfspTwoInsertNeighborhood(*instance);
    }
    else if(tm.checkToken(NEIGHBORHOOD_XTRANSPOSE))
    {
        prs::printTab( "XTranspose neighborhood");
        neigh = new emili::pfsp::XTransposeNeighborhood(*instance);
    }
    else if(tm.checkToken(NEIGHBORHOOD_TA_INSERT))
    {
        prs::printTab( "Insert with Taillard Acceleration");
        neigh = new emili::pfsp::TaillardAcceleratedInsertNeighborhood(*instance);
    }
    else if(tm.checkToken(NEIGHBORHOOD_TAx_INSERT))
    {
        prs::printTab( "Insert with Taillard Acceleration(Experimental)");
        neigh = new emili::pfsp::TAxInsertNeighborhood(*instance);
    }
    else if(tm.checkToken(NEIGHBORHOOD_ATAx_INSERT))
    {
        prs::printTab( "Approximated Insert with Taillard Acceleration(Experimental) for Weighted Tardiness");
        neigh = new emili::pfsp::OptInsert(*instance);
    }
    else if(tm.checkToken(NEIGHBORHOOD_HATAx_INSERT))
    {
        prs::printTab( "Heavily Approximated Insert with Taillard Acceleration(Experimental) for Weighted Tardiness");
        neigh = new emili::pfsp::HeavilyApproximatedTaillardAcceleratedInsertNeighborhood(*instance);
    }
    else if(tm.checkToken(NEIGHBORHOOD_NITA_INSERT))
    {
        prs::printTab( "Insert with Taillard Acceleration for no idle make span ");
        neigh = new emili::pfsp::NoIdleAcceleratedInsertNeighborhood(*instance);
    }
    else
    {
        std::cerr<< "'" << *tm << "' -> ERROR a neighborhood specification was expected! " << std::endl;
        std::cout << SAPFSPParser::info() << std::endl;
        exit(-1);
    }
    prs::decrementTabLevel();
    return neigh;
}


void SAPFSPParser::problem(prs::TokenManager& tm)
{
    PfspInstance i;    
    char *problem_type = tm.nextToken();
    bool ok;

    if(tm.checkToken(PROBLEM_SDSTPFS_MS))
    {
        ok = i.readSeqDepDataFromFile(tm.tokenAt(1));
    }
    else
    {
        ok = i.readDataFromFile(tm.tokenAt(1));

    }

    if(ok)
     {
         instance = instantiateSAPFSPProblem(problem_type, i);
        return;
     }

        std::cout << SAPFSPParser::info() << std::endl;
        exit(-1);
}


SAExploration* SAPFSPParser::EXPLORATION(prs::TokenManager& tm,
                                        emili::Neighborhood* neigh,
                                        SAAcceptance *acc,
                                        SACooling *cool,
                                        SATermination *term) {

    if (tm.checkToken(SARANDOMEXPLORATION)) {
        return new SARandomExploration(neigh, acc, cool, term);
    } else if (tm.checkToken(SASEQUENTIALEXPLORATION)) {
        return new SASequentialExploration(neigh, acc, cool, term);
    } else {
        std::cerr << "SAExploration expected, not found : " << std::endl;
        std::cerr << tm.peek() << std::endl;
        exit(1);
    }

}

SATempRestart* SAPFSPParser::TEMPRESTART(prs::TokenManager& tm,
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
}

emili::LocalSearch* SAPFSPParser::buildAlgo(prs::TokenManager& tm) {

    problem(tm);
    emili::InitialSolution* initsol    = init(tm);
    emili::Neighborhood*    nei        = neigh(tm);
    SAInitTemp*      inittemp   = INITTEMP(tm, initsol, nei, instance);
    SAAcceptance*    acceptance = ACCEPTANCE(tm, inittemp);
    SACooling*       cooling    = COOL(tm, inittemp, nei, instance);
    SATempRestart*   temprestart = TEMPRESTART(tm, inittemp, nei);
    cooling->setTempRestart(temprestart);
    SATermination*     term       = TERMINATION(tm, nei); // termin(tm);
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

}


bool SAPFSPParser::isParsable(std::string &problem)
{
    if(strcmp(problem.c_str(),PROBLEM_PFS_WT)==0)
    {
        return true;
    }else if(strcmp(problem.c_str(),PROBLEM_NWPFS_MS)==0)
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
    else if(strcmp(problem.c_str(),PROBLEM_NIPFS_MS)==0)
            {
        return true;
            }
    else if(strcmp(problem.c_str(),PROBLEM_SDSTPFS_MS)==0)
    {
        return true;
    }
    else
    {
        return false;
    }
}

std::string SAPFSPParser::info()
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
    oss << "ITERATED_LOCAL_SEARCH = ils LOCAL_SEARCH TERMINATION PERTUBATION ACCEPTANCE -it seconds" << "\n";
    oss << "TABU_SEARCH           = tabu INITIAL_SOLUTION TERMINATION NEIGHBORHOOD TABU_MEMORY" << "\n";
    oss << "VND_SEARCH            = vnd < first | best > INITIAL_SOLUTION TERMINATION NEIGHBORHOOD1 NEIGHBORHOOD2 ... NEIGHBORHOODn" << "\n";
    oss << "GVNS_SEARCH           = gvns INITIAL_SOLUTION PERTUBATION1 PERTUBATION2 -it seconds" << "\n";
    oss << "SEARCH_TYPE           = first | best | tabu | vnd | ils" << "\n";
    oss << "INITIAL_SOLUTION      = random | slack | nwslack | lit | rz | nrz | nrz2 | lr size(int)| nlr size(int) | mneh" <<"\n";
    oss << "TERMINATION           = true | time float | locmin | soater | iteration int | maxstep int" << "\n";
    oss << "NEIGHBORHOOD          = transpose | exchange | insert | binsert | finsert | tinsert | "<< NEIGHBORHOOD_TA_INSERT << " | " << NEIGHBORHOOD_NITA_INSERT<<"\n";
    oss << "PERTUBATION           = igper int | testper | rndmv NEIGHBORHOOD #moves(int) | noper (int) | nrzper (int) | tmiigper (int) (int) | igls (int) LOCAL_SEARCH | rsls (int) LOCAL_SEARCH" << "\n";
    oss << "ACCEPTANCE            = soaacc float | testacc #swaps(int) | metropolis start_temperature(float) | always (intensify | diversify) | improve | sa_metropolis start_temp end_temp ratio | pmetro start_temp end_temp ratio frequence(int) | saacc start_temp end_temp ratio frequence(int) alpha ]0,1] | tmiigacc start_temperature(float) | implat number_of_non_improving_steps_accepted plateau_threshold" << std::endl;
    oss << "TABU_MEMORY           = move size(int) | hash size(int) | solution size(int) | tsabm size(int)" << "\n";
   // std::cout << " syntax->EMILI instancefile search_type intial_solution termination neighborhood" << std::endl;
    return oss.str();
}
*/