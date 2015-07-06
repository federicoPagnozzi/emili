#include "paramsparser.h"
#include <cstdlib>
#include <cstdio>
#include <ctime>
#include <cstring>
#include <iostream>
#include <sstream>
#include <algorithm>
#include "pfspinstance.h"

/* Algos */
#define IG "ig"
#define ILS "ils"
#define TABU "tabu"
#define FIRST "first"
#define BEST "best"
#define VND "vnd"
#define GVNS_ILS "gvns"

/* tabu tenure types */
#define TABU_MEMORY_MOVES "move"
#define TABU_MEMORY_HASHES "hash"
#define TABU_MEMORY_SOLUTIONS "solution"
#define TABU_MEMORY_TSAB "tsabm"
#define TABU_MEMORY_TSAB_TEST "tsabmt"
#define TABU_MEMORY_VALUE "value"

/* modifiers */
#define RO "-ro"
#define IT "-it"
#define TS "-ts"
#define TI "-ti"
#define IN "-in"
#define TE "-te"
#define NE "-ne"
#define NS "-ns"
#define RNDSEED "rnds"
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

/* permutation flowshop solution pertubations */
#define PERTUBATION_RANDOM_MOVE "rndmv"
#define PERTUBATION_VNRANDOM_MOVE "vnrmv"
#define PERTUBATION_NOPER "noper"
#define PERTUBATION_RND "randpert"
#define PERTUBATION_NRZ "nrzper"
#define PERTUBATION_TMIIG "tmiigper"
#define PERTUBATION_SOA "igper"
#define PERTUBATION_SOA_LEGACY "soaper"
#define PERTUBATION_TEST "testper"
#define PERTUBATION_IGLS "igls"
#define PERTUBATION_RSLS "rsls"
#define PERTUBATION_RS "rsper"

/* acceptance criteria*/
#define ACCEPTANCE_PROB "prob"
#define ACCEPTANCE_METRO "metropolis"
#define ACCEPTANCE_PMETRO "pmetro"
#define ACCEPTANCE_TMIIG "tmiigacc"
#define ACCEPTANCE_IMPROVE_PLATEAU "implat"
#define ACCEPTANCE_TEST "testacc"
#define ACCEPTANCE_SOA "soaacc"
#define ACCEPTANCE_ALWAYS "always"
#define INTENSIFY "intensify"
#define DIVERSIFY "diversify"
#define ACCEPTANCE_IMPROVE "improve"
#define ACCEPTANCE_SA_METRO "sa_metropolis"
#define ACCEPTANCE_SA "saacc"

char* problem_type;

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

std::string prs::ParamsParser::info()
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

emili::pfsp::PfspNeighborhood* ne = nullptr;
emili::InitialSolution* in= nullptr;
emili::Termination* te= nullptr;
emili::TabuMemory* tmem= nullptr;
emili::Termination* ilt= nullptr;

std::vector< emili::Neighborhood*> nes;

emili::LocalSearch* prs::ParamsParser::eparams(prs::TokenManager& tm)
{
    prs::incrementTabLevel();
    emili::LocalSearch* ls;
    if(tm.checkToken(ILS))
    {
        printTab("ILS");
        ls = ils(tm);
    }  
    else if(tm.checkToken(GVNS_ILS))
    {
        printTab("GVNS...");
        ls = gvns(tm);
    }
    else
    {             
        ls = search(tm);
    }

    prs::decrementTabLevel();
    return ls;
}



emili::LocalSearch* prs::ParamsParser::search(prs::TokenManager& tm)
{
    prs::incrementTabLevel();
    emili::LocalSearch* ls;
    if(tm.checkToken(ILS))
    {
        printTab("ILS ");
        ls = ils(tm);

    }else if(tm.checkToken(TABU))
    {
        printTab("TABU SEARCH");
        ls = tparams(tm);
    }
    else if(tm.checkToken(FIRST))
    {
        printTab("FIRST IMPROVEMENT");
        params(tm);
        ls =  new emili::FirstImprovementSearch(*in,*te,*ne);
    }
    else if(tm.checkToken(BEST))
    {
        printTab("BEST IMPROVEMENT");
        params(tm);
        ls =  new emili::BestImprovementSearch(*in,*te,*ne);
    }
    else if(tm.checkToken(VND))
    {
        printTab("VND SEARCH");
        ls = vparams(tm);
    }
    else
    {
        std::cerr<< "'" << tm.peek() << "' -> ERROR a search definition was expected! " << std::endl;
        std::cout << info() << std::endl;
        exit(-1);
    }
    prs::decrementTabLevel();
    return ls;

}

emili::LocalSearch* prs::ParamsParser::ils(prs::TokenManager& tm)
{

    emili::LocalSearch* ls = search(tm);

    //ils_time = ilstime();
    emili::Termination* pft = term(tm);
    //emili::pfsp::PfspRandomSwapPertub* prsp = new emili::pfsp::PfspRandomSwapPertub(istance);
    int rpc = 5;
    emili::Perturbation* prsp = per(tm);
    //emili::AcceptanceCriteria* tac = new emili::pfsp::PfspTestAcceptance(istance);
    //emili::AcceptanceCriteria* tac = new emili::MetropolisAcceptance(1);
    emili::Acceptance* tac = acc(tm);//new emili::pfsp::SOAacceptance(1.2f);
    emili::LocalSearch* iils = new emili::IteratedLocalSearch(*ls,*pft,*prsp,*tac);   
    return iils;
}

emili::Perturbation* prs::ParamsParser::per(prs::TokenManager& tm)
{
    prs::incrementTabLevel();
    std::ostringstream oss;
    emili::Perturbation* per;
    if(tm.checkToken(PERTUBATION_SOA) || tm.checkToken(PERTUBATION_SOA_LEGACY))
    {
        int n = tm.getInteger();

        oss << "NEH destruct/construct pertubation which use objective function. number of job erased: "<<n;
        printTab(oss.str().c_str());
        per = new emili::pfsp::IGPerturbation(n,*istance);
    }else if(tm.checkToken(PERTUBATION_RS))
    {
        int n = tm.getInteger();

        oss << "NEH destruct/construct pertubation. number of job erased: "<<n;
        printTab(oss.str().c_str());
        per = new emili::pfsp::RSPertubation(n,*istance);
    }
    else if(tm.checkToken(PERTUBATION_IGLS))
    {
        int n = tm.getInteger();
        oss.str(""); oss  << "IG pertubation with local search applied on the partial solution. d = "<<n;
        printTab(oss.str().c_str());
        PfspInstance pfs = this->istance->getInstance();
        pfs.setNbJob(pfs.getNbJob()-n);
        emili::pfsp::PermutationFlowShop * pfse = instantiateProblem(problem_type,pfs);
        emili::pfsp::PermutationFlowShop* is = this->istance;
        this->istance = pfse;
        emili::LocalSearch* ll = search(tm);
        this->istance = is;
        per = new emili::pfsp::IgLsPertubation(n,*istance,ll);
    }
    else if(tm.checkToken(PERTUBATION_RSLS))
    {
        int n = tm.getInteger();
        oss.str(""); oss  << "IG pertubation with local search applied on the partial solution. d = "<<n;
        printTab(oss.str().c_str());
        PfspInstance pfs = this->istance->getInstance();
        pfs.setNbJob(pfs.getNbJob()-n);
        emili::pfsp::PermutationFlowShop * pfse = instantiateProblem(problem_type,pfs);
        emili::pfsp::PermutationFlowShop* is = this->istance;
        this->istance = pfse;
        emili::LocalSearch* ll = search(tm);
        this->istance = is;
        per = new emili::pfsp::RSLSPertubation(n,*istance,ll);
    }
    else if(tm.checkToken(PERTUBATION_TEST))
    {
        oss.str(""); oss<< "Random swap test pertubation.";
        printTab(oss.str().c_str());
        per = new emili::pfsp::PfspRandomSwapPertub(*istance);
    }else if(tm.checkToken(PERTUBATION_RANDOM_MOVE))
    {
        printTab("Random move perturbation.");
        emili::Neighborhood* n = neigh(tm);
        int num = tm.getInteger();
        prs::incrementTabLevel();
        oss.str(""); oss  << "number of moves per pertubation step " << num;
        printTab(oss.str().c_str());
        prs::decrementTabLevel();
        per = new emili::RandomMovePertubation(*n,num);
    }
    else if(tm.checkToken(PERTUBATION_NOPER))
    {
        printTab("No pertubation.");
        per = new emili::NoPertubation();
    }
    else if(tm.checkToken(PERTUBATION_NRZ))
    {
        int n = tm.getInteger();
        oss.str(""); oss  << "neh rz destruct/construct pertubation. number of job erased: "<<n;
        printTab(oss.str().c_str());
        per = new emili::pfsp::NRZPertubation(n,*istance);
    }else if(tm.checkToken(PERTUBATION_VNRANDOM_MOVE))
    {
        printTab("Random move perturbation." );
        prs::incrementTabLevel();
        int num = tm.getInteger();
        oss.str(""); oss  << "number of moves per pertubation step " << num << ".\n\t";
        printTab(oss.str().c_str());
        int iter = tm.getInteger();
        oss.str(""); oss  << "number of iteration before changing the neighborhood " << iter << ".\n\t";
        printTab(oss.str().c_str());
        nes.clear();
        prs::decrementTabLevel();
        neighs(tm);
        per = new emili::VNRandomMovePertubation(nes,num,iter);
    }
    else if(tm.checkToken(PERTUBATION_TMIIG))
    {
        int num = tm.getInteger();
        int tsize = tm.getInteger();
        oss.str(""); oss  << "TMIIG pertubation. Number of job erased " << num << ". tabu list size " << tsize <<".\n\t";
        printTab(oss.str().c_str());
        per = new emili::pfsp::TMIIGPertubation(num,*istance,tsize);
    }
    else
    {
        std::cerr<< "'" << tm.peek() << "' -> ERROR a pertubation criteria specification was expected! " << std::endl;
        std::cout << info() << std::endl;
        exit(-1);
    }
    prs::decrementTabLevel();
    return per;
}

emili::Acceptance* prs::ParamsParser::acc(prs::TokenManager& tm)
{
    prs::incrementTabLevel();
    emili::Acceptance* acc;
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
    else  if(tm.checkToken(ACCEPTANCE_METRO))
    {
        float n = tm.getDecimal();
        oss.str(""); oss  << "metropolis acceptance. temperature : "<<n;
        printTab(oss.str().c_str());
        acc = new  emili::MetropolisAcceptance(n);
    }
    else  if(tm.checkToken(ACCEPTANCE_ALWAYS))
    {


        emili::accept_candidates accc;
        char* t1;
        if(tm.checkToken(INTENSIFY))
        {
            accc = emili::ACC_INTENSIFICATION;
            t1 = INTENSIFY;
        }
        else if(tm.checkToken(DIVERSIFY))
        {
            t1 = DIVERSIFY;
            accc = emili::ACC_DIVERSIFICATION;
        }
        else
        {
            std::cerr<< "'" << *tm << "' -> ERROR " << INTENSIFY << " or " << DIVERSIFY <<" was expected! " << std::endl;
            std::cout << info() << std::endl;
        exit(-1);
        }
        oss.str(""); oss  << "Acceptance always "<< t1;
        printTab(oss.str().c_str());
        acc = new  emili::AlwaysAccept(accc);
    }
    else  if(tm.checkToken(ACCEPTANCE_IMPROVE))
    {

        printTab( "improve acceptance");

        acc = new  emili::ImproveAccept();
    }
    else  if(tm.checkToken(ACCEPTANCE_SA_METRO))
    {
        float start =tm.getDecimal();
        float end =tm.getDecimal();
        float ratio =tm.getDecimal();
        oss.str(""); oss  << "metropolis acceptance. start ,end , ratio : "<< start << ", "<< end << "," << ratio;
        printTab(oss.str().c_str());
        acc = new  emili::Metropolis(start,end,ratio);
    }
    else  if(tm.checkToken(ACCEPTANCE_PMETRO))
    {
        float start =tm.getDecimal();
        float end =tm.getDecimal();
        float ratio =tm.getDecimal();
        int iterations = tm.getInteger();
        oss.str(""); oss  << "metropolis acceptance. start ,end , ratio, frequence : "<< start << ", "<< end << "," << ratio <<","<< iterations;
        printTab(oss.str().c_str());
        acc = new  emili::Metropolis(start,end,ratio,iterations);
    }
    else if(tm.checkToken(ACCEPTANCE_SA))
    {
        float start =tm.getDecimal();
        float end =tm.getDecimal();
        float ratio =tm.getDecimal();
        int iterations = tm.getInteger();
        float alpha =tm.getDecimal();
        oss.str(""); oss  << "metropolis acceptance. start ,end , ratio, frequence, alpha : "<< start << ", "<< end << "," << ratio <<","<< iterations << "," << alpha;
        printTab(oss.str().c_str());
        acc = new  emili::Metropolis(start,end,ratio,iterations,alpha);
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
    }else if(tm.checkToken(ACCEPTANCE_IMPROVE_PLATEAU))
    {
        int plateau_steps = tm.getInteger();
        int threshold = tm.getInteger();
        oss.str(""); oss  << "Accept a diversification solution if it improves on the intensification otherwise it will accept "<< plateau_steps << " non improving steps once it reaches the threshold of " << threshold;
        printTab(oss.str().c_str());
        acc = new  emili::AcceptPlateau(plateau_steps,threshold);
    }
    else
    {
        std::cerr<< "'" << *tm << "' -> ERROR an acceptance criteria specification was expected! " << std::endl;
        std::cout << info() << std::endl;
        exit(-1);
    }
    prs::decrementTabLevel();
    return acc;
}


emili::LocalSearch* prs::ParamsParser::gvns(prs::TokenManager& tm)
{
    emili::InitialSolution* is = init(tm);
    emili::Perturbation* p1 = per(tm);
    emili::Perturbation* p2 = per(tm);
    emili::pfsp::GVNS_innerloop* gvi = new emili::pfsp::GVNS_innerloop(*is);
    emili::Perturbation* p3 = per(tm);
    emili::Perturbation* p4 = per(tm);

    std::vector< emili::Perturbation* > p;
    p.push_back(p1);
    p.push_back(p2);

    std::vector< emili::Perturbation* > pl;
    pl.push_back(p3);
    pl.push_back(p4);

    emili::GVNS* gg = new emili::GVNS(*gvi,p);
    return new emili::GVNS(*gvi,pl);
}

emili::TabuSearch* prs::ParamsParser::tparams(prs::TokenManager& tm)
{
    params(tm);
    tmem = tmemory(ne,tm);
    //std::pair<int,int > tset = tsettings();
   // tm->setTabuTenure(tset.first);
    //emili::pfsp::PfspTerminationIterations* ptc = new emili::pfsp::PfspTerminationIterations(tset.second);
    return new emili::TabuSearch(*in,*te,*ne,*tmem);
}

emili::TabuMemory* prs::ParamsParser::tmemory(emili::pfsp::PfspNeighborhood* n,prs::TokenManager& tm)
{
    prs::incrementTabLevel();

    std::ostringstream oss;
    emili::TabuMemory* tmem;


    printTab(oss.str().c_str());

    if(tm.checkToken(TABU_MEMORY_MOVES))
    {
        oss.str(""); oss << "USING MOVES\n\t";
        printTab(oss.str().c_str());
        int ts = tm.getInteger();
        oss << "Tabu tenure size " << ts;
        tmem = new  emili::pfsp::PfspMovesMemory(ts , *n);
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
        tmem = new  emili::pfsp::TSABMemory(ts , *n);
    }
    else if(tm.checkToken(TABU_MEMORY_TSAB_TEST))
    {
        oss.str(""); oss << "USING TSAB\n\t";
        printTab(oss.str().c_str());
        int ts = tm.getInteger();
        oss << "Tabu tenure size " << ts;
        tmem = new  emili::pfsp::TSABtestMemory(ts , *n);
    }
    else if(tm.checkToken(TABU_MEMORY_VALUE))
    {
        oss.str(""); oss << "USING VALUE\n\t" ;
        printTab(oss.str().c_str());
        int ts = tm.getInteger();
        oss << "Tabu tenure size " << ts;
        tmem = new  emili::pfsp::PfspTabuValueMemory(ts);
    }
    else
    {
        std::cerr<< "'" << *tm << "' -> ERROR a memory specification for the tabu search was expected! " << std::endl;
        std::cout << info() << std::endl;
        exit(-1);
    }
    prs::decrementTabLevel();
    return tmem;
}

void prs::ParamsParser::params(prs::TokenManager& tm)
{
    in = init(tm);
    te = term(tm);
    ne = neigh(tm);
}

emili::LocalSearch* prs::ParamsParser::vparams(prs::TokenManager& tm)
{
    prs::incrementTabLevel();
    emili::LocalSearch* ls;
    if(tm.checkToken(FIRST))
    {
        printTab("FIRST IMPROVEMENT VND");
        in = init(tm);
        te = term(tm);
        neighs(tm);
        ls =  new emili::VNDSearch<emili::FirstImprovementSearch>(*in,*te,nes);
    }
    else if(tm.checkToken(BEST))
    {
       printTab("BEST IMPROVEMENT VND");
       in = init(tm);
       te = term(tm);
       neighs(tm);
        ls =  new emili::VNDSearch<emili::BestImprovementSearch>(*in,*te,nes);
    }
    else
    {
        std::cerr<< "'" << *tm << "' -> ERROR a valid type of search must be specified (first,best) " << std::endl;

        std::cout << info() << std::endl;
        exit(-1);
    }
    prs::decrementTabLevel();
    return ls;
}

emili::InitialSolution* prs::ParamsParser::init(prs::TokenManager& tm)
{
    prs::incrementTabLevel();
    std::ostringstream oss;
    emili::InitialSolution* init;
    if(tm.checkToken(INITIAL_RANDOM))
    {
        printTab("Random initial solution");
        init = new emili::pfsp::PfspRandomInitialSolution(*istance);
    }
    else if(tm.checkToken(INITIAL_SLACK))
    {
        printTab("SLACK initial solution");
        init = new emili::pfsp::PfspSlackInitialSolution(*istance);
    }else if(tm.checkToken(INITIAL_WNSLACK))
    {
        printTab( "NEH WSLACK initial solution");
        //init = new testIS(istance);
        init = new emili::pfsp::PfspNEHwslackInitialSolution(*istance);
    }
    else if(tm.checkToken(INITIAL_LIT))
        {
            printTab( "Less idle times initial solution");
            //return new testIS(istance);
            init = new emili::pfsp::LITSolution(*istance);
        }
    else if(tm.checkToken(INITIAL_RZ))
        {
            printTab( "rz initial solution");
            //return new testIS(istance);
            init = new emili::pfsp::RZSolution(*istance);
        }
    else if(tm.checkToken(INITIAL_NRZ))
        {
            printTab( "neh rz initial solution");
            //return new testIS(istance);
            init = new emili::pfsp::NeRZSolution(*istance);
        }
    else if(tm.checkToken(INITIAL_NRZ2))
        {
            printTab( "neh rz initial solution without improvement phase");
            //return new testIS(*istance);
            init = new emili::pfsp::NeRZ2Solution(*istance);
        }
    else if(tm.checkToken(INITIAL_LR))
        {
            int n = tm.getInteger();
            oss.str(""); oss << "LR initial solution with "<<n<<" starting sequences";
            printTab(oss.str().c_str());
            // testIS(*istance);
            init = new emili::pfsp::LRSolution(*istance,n);
        }
    else if(tm.checkToken(INITIAL_NLR))
        {
        int n = tm.getInteger();
        oss.str(""); oss << "NLR initial solution with "<<n<<" starting sequences";
        //return new testIS(*istance);printTab(oss.str().c_str());
        printTab(oss.str().c_str());
        init = new emili::pfsp::NLRSolution(*istance,n);
        }
    else if(tm.checkToken(INITIAL_MNEH))
        {
            printTab( "mneh initial solution");
            //return new testIS(istance);
            init = new emili::pfsp::MNEH(*istance);
        }
    else
    {
        std::cerr<< "'" << *tm << "' -> ERROR a initial solution generator specification was expected! (random,slack)" << std::endl;

        std::cout << info() << std::endl;
        exit(-1);
    }
    prs::decrementTabLevel();
    return init;
}

emili::Termination* prs::ParamsParser::term(prs::TokenManager& tm)
{
    prs::incrementTabLevel();
    emili::Termination* term;
    if(tm.checkToken(TERMINATION_LOCMIN))
    {
        printTab("Local minima termination");
        term = new emili::LocalMinimaTermination();
    }
    else if(tm.checkToken(TERMINATION_WTRUE))
    {
        printTab("While true termination");
        term = new emili::WhileTrueTermination();
    }
    else if(tm.checkToken(TERMINATION_ITERA))
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
    else if(tm.checkToken(TERMINATION_TIME))
    {

        float time =tm.getDecimal();
        if(time==0){
            time = 1;
        }
        std::ostringstream oss;
        oss << "Timed termination. ratio: " << time;
        printTab(oss.str().c_str());
        term =  new emili::TimedTermination(time);
    }
    else if(tm.checkToken(TERMINATION_MAXSTEPS))
    {
        int steps = tm.getInteger();
        std::ostringstream oss;
        oss << "Max Steps termination. # steps: "<< steps;
        printTab(oss.str().c_str());
        term = new emili::MaxStepsTermination(steps);
    }
    else
    {
        std::cerr<< "'" << *tm << "' -> ERROR a termination criteria specification was expected! " << std::endl;
        std::cout << info() << std::endl;
        exit(-1);
    }
    prs::decrementTabLevel();
    return term;
}

emili::pfsp::PfspNeighborhood* prs::ParamsParser::neigh(prs::TokenManager& tm)
{
    prs::incrementTabLevel();
    emili::pfsp::PfspNeighborhood* neigh;
    if(tm.checkToken(NEIGHBORHOOD_INSERT))
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
    else if(tm.checkToken(NEIGHBORHOOD_TAx_INSERT))
    {
        printTab( "Insert with Taillard Acceleration(Experimental)");
        neigh = new emili::pfsp::TAxInsertNeighborhood(*istance);
    }
    else if(tm.checkToken(NEIGHBORHOOD_ATAx_INSERT))
    {
        printTab( "Approximated Insert with Taillard Acceleration(Experimental) for Weighted Tardiness");
        neigh = new emili::pfsp::ApproximatedTaillardAcceleratedInsertNeighborhood(*istance);
    }
    else if(tm.checkToken(NEIGHBORHOOD_HATAx_INSERT))
    {
        printTab( "Heavily Approximated Insert with Taillard Acceleration(Experimental) for Weighted Tardiness");
        neigh = new emili::pfsp::HeavilyApproximatedTaillardAcceleratedInsertNeighborhood(*istance);
    }
    else if(tm.checkToken(NEIGHBORHOOD_NITA_INSERT))
    {
        printTab( "Insert with Taillard Acceleration for no idle make span ");
        neigh = new emili::pfsp::NoIdleAcceleratedInsertNeighborhood(*istance);
    }
    else
    {
        std::cerr<< "'" << *tm << "' -> ERROR a neighborhood specification was expected! " << std::endl;
        std::cout << info() << std::endl;
        exit(-1);
    }
    prs::decrementTabLevel();
    return neigh;
}

emili::pfsp::PfspNeighborhood* prs::ParamsParser::neighV(prs::TokenManager& tm)
{
    prs::incrementTabLevel();
    emili::pfsp::PfspNeighborhood* neigh;
    if(tm.checkToken(NEIGHBORHOOD_INSERT))
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
    else if(tm.checkToken(NEIGHBORHOOD_TAx_INSERT))
    {
        printTab( "Insert with Taillard Acceleration");
        neigh = new emili::pfsp::TAxInsertNeighborhood(*istance);
    }    
    else if(tm.checkToken(NEIGHBORHOOD_NITA_INSERT))
    {
        printTab( "Insert with Taillard Acceleration for no idle make span ");
        neigh = new emili::pfsp::NoIdleAcceleratedInsertNeighborhood(*istance);
    }
    else
    {	
        neigh =  nullptr;
    }
    prs::decrementTabLevel();
    return neigh;
}

void prs::ParamsParser::neighs(prs::TokenManager& tm)
{
    std::vector<emili::Neighborhood*> vnds;
    vnds.push_back(neigh(tm));
    nes = vnds;
    neighs1(tm);
}

void prs::ParamsParser::neighs1(prs::TokenManager& tm)
{
    emili::Neighborhood* n = neighV(tm);
	if(n!=nullptr)
	{
           nes.push_back(n);
           neighs1(tm);
	}

}


void prs::ParamsParser::problem(prs::TokenManager& tm)
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
         istance = instantiateProblem(problem_type, i);
        return;
     }

        std::cout << info() << std::endl;
        exit(-1);
}


emili::LocalSearch* prs::ParamsParser::buildAlgo(prs::TokenManager& tm)
{
    problem(tm);
    emili::LocalSearch* local = eparams(tm);
    std::cout << "------" << std::endl;
    return local;
}

bool prs::ParamsParser::isParsable(string &problem)
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

std::string prs::ParamsParser::availableProblems() const
{
    ostringstream oss;
    oss <<PROBLEM_PFS_WT<< " " <<PROBLEM_PFS_WE<< " " <<PROBLEM_PFS_TCT<< " " <<PROBLEM_PFS_MS
       << " " <<PROBLEM_PFS_WCT<< " " <<PROBLEM_PFS_T<< " " <<PROBLEM_PFS_E<< " "<<PROBLEM_NWPFS_WT<< " " <<PROBLEM_NWPFS_WE
       << " " <<PROBLEM_NWPFS_TCT<< " " <<PROBLEM_NWPFS_MS<< " " <<PROBLEM_NWPFS_T<< " " <<PROBLEM_NWPFS_E
       << PROBLEM_NIPFS_MS <<" "<<PROBLEM_NIPFS_WT<< " " <<PROBLEM_NIPFS_WE<< " " <<PROBLEM_NIPFS_TCT<< " " <<PROBLEM_NIPFS_MS<< " "
       << " " <<PROBLEM_NIPFS_T<< " " <<PROBLEM_NIPFS_E << PROBLEM_SDSTPFS_MS;

    return oss.str();
}
