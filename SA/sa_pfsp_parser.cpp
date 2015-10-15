#include "sa_pfsp_parser.h"

#include <cstdlib>
#include <cstdio>
#include <ctime>
#include <cstring>
#include <iostream>
#include <sstream>
#include <algorithm>
#include "../pfspinstance.h"
#include "../permutationflowshop.h"

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
#define PROBLEM_NWPFS_TCT "NWPFSP_TCT"
#define PROBLEM_NWPFS_T "NWPFSP_T"
#define PROBLEM_NWPFS_E "NWPFSP_E"

/* no idle permutation flowshop*/
#define PROBLEM_NIPFS_MS "NIPFSP_MS"
#define PROBLEM_NIPFS_WT "NIPFSP_WT"
#define PROBLEM_NIPFS_WE "NIPFSP_WE"
#define PROBLEM_NIPFS_TCT "NIPFSP_TCT"
#define PROBLEM_NIPFS_T "NIPFSP_T"
#define PROBLEM_NIPFS_E "NIPFSP_E"

/* Sequence dependent setup times */
#define PROBLEM_SDSTPFS_MS "SDSTPFS_MS"


/* initial solution heuristics */
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

/* Termination criteria*/
#define TERMINATION_MAXSTEPS "maxstep"
#define TERMINATION_TIME "time"
#define TERMINATION_LOCMIN "locmin"
#define TERMINATION_ITERA "iteration"
#define TERMINATION_WTRUE "true"
#define TERMINATION_SOA "soater"

/* permutation flowshop neighborhoods*/
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


char* problem_type;

//using namespace prs;

emili::pfsp::PermutationFlowShop* instantiateProblem(char* t, PfspInstance i)
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
    } else {
        std::cerr << "SAInitTemp expected, not found : " << std::endl;
        std::cerr << tm.peek() << std::endl;
        exit(1);
    }

}


SAAcceptance* SAPFSPParser::ACCEPTANCE(prs::TokenManager& tm) {

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
    } else {
        std::cerr << "SAAcceptance expected, not found : " << std::endl;
        std::cerr << tm.peek() << std::endl;
        exit(1);
    }

}


SATermination* SAPFSPParser::TERMINATION(prs::TokenManager& tm) {

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


SACooling* SAPFSPParser::COOL(prs::TokenManager& tm,
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
    } else {
        std::cerr << "SACooling expected, not found : " << std::endl;
        std::cerr << tm.peek() << std::endl;
        exit(1);
    }

}


SATempLength* SAPFSPParser::TEMPLENGTH(prs::TokenManager& tm,
                                       emili::Neighborhood* neigh) {

    if (tm.checkToken(CONSTTEMPLEN)) {
        int l = tm.getInteger();
        return new ConstantTempLength(l);
    } else if (tm.checkToken(NEIGHSIZETEMPLEN)) {
        float a = tm.getDecimal();
        return new NeighSizeTempLength(neigh, a);
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
        neigh = new emili::pfsp::ApproximatedTaillardAcceleratedInsertNeighborhood(*instance);
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
    problem_type = tm.nextToken();
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
         instance = instantiateProblem(problem_type, i);
        return;
     }

        std::cout << SAPFSPParser::info() << std::endl;
        exit(-1);
}


SAExploration* SAPFSPParser::EXPLORATION(prs::TokenManager& tm,
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

SATempRestart* SAPFSPParser::TEMPRESTART(prs::TokenManager& tm,
                                         SAInitTemp *it) {

    if (tm.checkToken(SANOTEMPRESTART)) {
        return new SANoRestart();
    } else if (tm.checkToken(SAMINTEMPRESTART)) {
        float va = tm.getDecimal();
        return new SAMinRestart(it, va);
    } else if (tm.checkToken(SAPERCTEMPRESTART)) {
        float va = tm.getDecimal();
        return new SAPercRestart(it, va);
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
    SAInitTemp*      inittemp   = INITTEMP(tm, initsol);
    SAAcceptance*    acceptance = ACCEPTANCE(tm);
    SACooling*       cooling    = COOL(tm, inittemp);
    SATempRestart*   temprestart = TEMPRESTART(tm, inittemp);
    cooling->setTempRestart(temprestart);
    SATermination*     term       = TERMINATION(tm); // termin(tm);
    SATempLength*    templ      = TEMPLENGTH(tm, nei);
    SAExploration* explo = EXPLORATION(tm, nei, acceptance, term);

    return new SimulatedAnnealing(initsol,
                                  inittemp,
                                  acceptance,
                                  cooling,
                                  term,
                                  templ,
                                  explo,
                                  nei);

}


bool SAPFSPParser::isParsable(string &problem)
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
    ostringstream oss;
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