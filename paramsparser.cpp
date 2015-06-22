#include "paramsparser.h"
#include <cstdlib>
#include <cstdio>
#include <ctime>
#include <cstring>
#include <iostream>
#include <sstream>
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
#define DEFAULT_IT -10


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
#define NEIGHBORHOOD_NITA_INSERT "ntainsert"

/* permutation flowshop solution pertubations */
#define PERTUBATION_RANDOM_MOVE "rndmv"
#define PERTUBATION_VNRANDOM_MOVE "vnrmv"
#define PERTUBATION_NOPER "noper"
#define PERTUBATION_RND "randpert"
#define PERTUBATION_NRZ "nrzper"
#define PERTUBATION_TMIIG "tmiigper"
#define PERTUBATION_SOA "soaper"
#define PERTUBATION_TEST "testper"
#define PERTUBATION_IGLS "igls"

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


int tab_level = 0;

void printTab(const char* string)
{
    for(int i=0;i<tab_level; i++)
    {
        std::cout << "  ";
    }

    std::cout << string << std::endl;
}

char* problem_type;

emili::pfsp::PermutationFlowShop* instantiateProblem(char* t, PfspInstance i)
{
    emili::pfsp::PermutationFlowShop* prob;
    if(strcmp(t,PROBLEM_PFS_WT)==0)
    {
        printTab("Permutation Flow Shop Weighted Tardiness");
        prob = new emili::pfsp::PFSP_WT(i);
    }else if(strcmp(t,PROBLEM_NWPFS_MS)==0)
    {
        printTab("No Wait Permutation Flow Shop Make Span");
        prob = new emili::pfsp::NWPFSP_MS(i);
    }else if(strcmp(t,PROBLEM_PFS_E)==0)
    {
        printTab("Permutation Flow Shop Earliness");
        prob = new emili::pfsp::PFSP_E(i);
    }else if(strcmp(t,PROBLEM_PFS_WE)==0)
    {
        printTab("Permutation Flow Shop Weighted Earliness");
        prob = new emili::pfsp::PFSP_WE(i);
    }else if(strcmp(t,PROBLEM_PFS_T)==0)
    {
        printTab("Permutation Flow Shop Tardiness");
        prob = new emili::pfsp::PFSP_T(i);
    }else if(strcmp(t,PROBLEM_PFS_MS)==0)
    {
        printTab("Permutation Flow Shop Make Span");
        prob = new emili::pfsp::PFSP_MS(i);
    }
    else if(strcmp(t,PROBLEM_NIPFS_MS)==0)
            {
                printTab("No Idle Permutation Flow Shop Make Span" );
                prob = new emili::pfsp::NI_A_PFSP_MS(i);
            }
    else
    {
        std::cerr<< "'" << t << "' -> ERROR a problem was expected! " << std::endl;
        prs::info();
    exit(-1);
    }
    return prob;
}


void prs::emili_header()
{
    std::cout << "\t ______ __  __ _____ _      _____ " << std::endl;
    std::cout << "\t|  ____|  \\/  |_   _| |    |_   _|" << std::endl;
    std::cout << "\t| |__  | \\  / | | | | |      | |  " << std::endl;
    std::cout << "\t|  __| | |\\/| | | | | |      | |  " << std::endl;
    std::cout << "\t| |____| |  | |_| |_| |____ _| |_ " << std::endl;
    std::cout << "\t|______|_|  |_|_____|______|_____|" << std::endl;
    std::cout << std::endl;
}

void prs::info()
{
    std::cout << "Usage:" << std::endl;
    std::cout << std::endl;
    std::cout << "EMILI INSTANCE_FILE_PATH PROBLEM <LOCAL_SEARCH | ITERATED_LOCAL_SEARCH | TABU_SEARCH | VND_SEARCH> [rnds seed]" << std::endl;
    std::cout << std::endl;
    std::cout << "PROBLEM               = "<<PROBLEM_PFS_WT<< " " <<PROBLEM_PFS_WE<< " " <<PROBLEM_PFS_TCT<< " " <<PROBLEM_PFS_MS
              << " " <<PROBLEM_PFS_WCT<< " " <<PROBLEM_PFS_T<< " " <<PROBLEM_PFS_E<< " "<<PROBLEM_NWPFS_WT<< " " <<PROBLEM_NWPFS_WE
              << " " <<PROBLEM_NWPFS_TCT<< " " <<PROBLEM_NWPFS_MS<< " " <<PROBLEM_NWPFS_T<< " " <<PROBLEM_NWPFS_E
              << PROBLEM_NIPFS_MS <<" "<<PROBLEM_NIPFS_WT<< " " <<PROBLEM_NIPFS_WE<< " " <<PROBLEM_NIPFS_TCT<< " " <<PROBLEM_NIPFS_MS<< " "
              << " " <<PROBLEM_NIPFS_T<< " " <<PROBLEM_NIPFS_E << std::endl;
    std::cout << "LOCAL_SEARCH          = SEARCH_TYPE INITIAL_SOLUTION TERMINATION NEIGHBORHOOD" << std::endl;
    std::cout << "ITERATED_LOCAL_SEARCH = ils LOCAL_SEARCH TERMINATION PERTUBATION ACCEPTANCE -it seconds" << std::endl;
    std::cout << "TABU_SEARCH           = tabu INITIAL_SOLUTION TERMINATION NEIGHBORHOOD TABU_MEMORY" << std::endl;
    std::cout << "VND_SEARCH            = vnd < first | best > INITIAL_SOLUTION TERMINATION NEIGHBORHOOD1 NEIGHBORHOOD2 ... NEIGHBORHOODn" << std::endl;
    std::cout << "GVNS_SEARCH           = gvns INITIAL_SOLUTION PERTUBATION1 PERTUBATION2 -it seconds" << std::endl;
    std::cout << "SEARCH_TYPE           = first | best | tabu | vnd | ils" << std::endl;
    std::cout << "INITIAL_SOLUTION      = random | slack | nwslack | lit | rz | nrz | nrz2 | lr size(int)| nlr size(int) | mneh" << std::endl;
    std::cout << "TERMINATION           = true | time float | locmin | soater | iteration int | maxstep int" << std::endl;
    std::cout << "NEIGHBORHOOD          = transpose | exchange | insert | binsert | finsert | tinsert | "<< NEIGHBORHOOD_TA_INSERT << " | " << NEIGHBORHOOD_NITA_INSERT<< std::endl;
    std::cout << "PERTUBATION           = soaper int | testper | rndmv NEIGHBORHOOD #moves(int) | noper (int) | nrzper (int) | tmiigper (int) (int)" << std::endl;
    std::cout << "ACCEPTANCE            = soaacc float | testacc #swaps(int) | metropolis start_temperature(float) | always (intensify | diversify) | improve | sa_metropolis start_temp end_temp ratio | pmetro start_temp end_temp ratio frequence(int) | saacc start_temp end_temp ratio frequence(int) alpha ]0,1] | tmiigacc start_temperature(float) | implat number_of_non_improving_steps_accepted plateau_threshold" << std::endl;
    std::cout << "TABU_MEMORY           = move size(int) | hash size(int) | solution size(int) | tsabm size(int)" << std::endl;
   // std::cout << " syntax->EMILI instancefile search_type intial_solution termination neighborhood" << std::endl;
}

void check(char* t,const char* message)
{
    if(t==nullptr)
    {
        std::cerr <<"PARSING ERROR"<< message << std::endl;
        prs::info();
        exit(-1);
    }
}

emili::pfsp::PfspNeighborhood* ne = nullptr;
emili::InitialSolution* in= nullptr;
emili::Termination* te= nullptr;
emili::TabuMemory* tmem= nullptr;
emili::Termination* ilt= nullptr;

std::vector< emili::Neighborhood*> nes;

emili::LocalSearch* prs::ParamsParser::eparams()
{

    char* t = nextToken();
    tab_level++;
    check(t,"SEARCH TYPE AND PARAMETERS MISSING!!!");
    emili::LocalSearch* ls;
    if(strcmp(t,ILS)==0)
    {
        printTab("ILS");
        ls = ils();

    }
    else if(strcmp(t,IG)==0)
    {
        printTab("Iterated Greedy");
        ls = ig();

    }
    else if(strcmp(t,GVNS_ILS)==0)
    {
        printTab("GVNS...");
        ls = gvns();
        ils_time = ilstime();
        ils_time = (istance->getNjobs()*(istance->getNmachines()/2)*ils_time)/1000;
        ls->setSearchTime(ils_time);

    }
    else
    {
        currentToken--;
        //std::cout << "LocalSearch \n\t";
        ls = search();

        ils_time = ilstime();
        ls->setSearchTime(ils_time);
    }
    int seed = getSeed();
    std::ostringstream oss;
    oss << "\tRANDOM SEED " << seed;

    printTab(oss.str().c_str());
    emili::initializeRandom(seed);
    tab_level--;

    /*std::cout << "test random " << emili::generateRandomNumber() << " " << emili::generateRandomNumber()<< " " << emili::generateRandomNumber()<< " " << emili::generateRandomNumber()<< " " << emili::generateRandomNumber()<< " " << emili::generateRandomNumber()<< " " << emili::generateRandomNumber() <<std::endl;
    std::cout << "test random real" << emili::generateRealRandomNumber() << " " << emili::generateRealRandomNumber() <<std::endl;*/
    return ls;
}

int prs::ParamsParser::getSeed()
{
    char* t = nextToken();
    if(t != nullptr && strcmp(t,RNDSEED)==0)
    {
        //currentToken--;

        return number();
    }
    return 0;
}

emili::LocalSearch* prs::ParamsParser::search()
{
    tab_level++;
    char* t = nextToken();
    check(t,"SEARCH PARAMETERS MISSING!!!");
    emili::LocalSearch* ls;
    if(strcmp(t,ILS)==0)
    {
        printTab("ILS ");
        ls = ils();

    }else if(strcmp(t,TABU)==0)
    {
        printTab("TABU SEARCH");
        ls = tparams();
    }
    else if(strcmp(t,FIRST)==0)
    {
        printTab("FIRST IMPROVEMENT");
        params();
        ls =  new emili::FirstImprovementSearch(*in,*te,*ne);
    }
    else if(strcmp(t,BEST)==0)
    {
        printTab("BEST IMPROVEMENT");
        params();
        ls =  new emili::BestImprovementSearch(*in,*te,*ne);
    }
    else if(strcmp(t,VND)==0)
    {
        printTab("VND SEARCH");
        ls = vparams();
    }
    else
    {
        std::cerr<< "'" << t << "' -> ERROR a search definition was expected! " << std::endl;
        prs::info();
        exit(-1);
    }
    tab_level--;
    return ls;

}

emili::LocalSearch* prs::ParamsParser::ils()
{

    emili::LocalSearch* ls = search();

    //ils_time = ilstime();
    emili::Termination* pft = term();
    //emili::pfsp::PfspRandomSwapPertub* prsp = new emili::pfsp::PfspRandomSwapPertub(istance);
    int rpc = 5;
    emili::Perturbation* prsp = per();
    //emili::AcceptanceCriteria* tac = new emili::pfsp::PfspTestAcceptance(istance);
    //emili::AcceptanceCriteria* tac = new emili::MetropolisAcceptance(1);
    emili::Acceptance* tac = acc();//new emili::pfsp::SOAacceptance(1.2f);
    emili::LocalSearch* iils = new emili::IteratedLocalSearch(*ls,*pft,*prsp,*tac);
    ils_time = ilstime();
    if(ils_time>0)
    {
        iils->setSearchTime(ils_time);
        //std::cerr <<"ERROR for ils a time has to be provided"<< std::endl;
        //exit(-1);
    }

    return iils;
}

emili::Perturbation* prs::ParamsParser::per()
{
    tab_level++;
    char* t = nextToken();
    check(t,"PERTUBATION CRITERIA EXPECTED!");
    std::ostringstream oss;
    emili::Perturbation* per;
    if(strcmp(t,PERTUBATION_SOA)==0)
    {
        int n = number();

        oss << "wslack destruct/construct pertubation. number of job erased: "<<n;
        printTab(oss.str().c_str());
        per = new emili::pfsp::SOAPerturbation(n,*istance);
    }
    else if(strcmp(t,PERTUBATION_IGLS)==0)
    {
        int n = number();
        oss.str(""); oss  << "IG pertubation with local search applied on the partial solution. d = "<<n;
        printTab(oss.str().c_str());
        PfspInstance pfs = this->istance->getInstance();
        pfs.setNbJob(pfs.getNbJob()-n);
        emili::pfsp::PermutationFlowShop * pfse = instantiateProblem(problem_type,pfs);
        emili::pfsp::PermutationFlowShop* is = this->istance;
        this->istance = pfse;
        emili::LocalSearch* ll = search();
        this->istance = is;
        per = new emili::pfsp::IgLsPertubation(n,*istance,ll);
    }
    else if(strcmp(t,PERTUBATION_TEST)==0)
    {
        oss.str(""); oss<< "Random swap test pertubation.";
        printTab(oss.str().c_str());
        per = new emili::pfsp::PfspRandomSwapPertub(*istance);
    }else if(strcmp(t,PERTUBATION_RANDOM_MOVE)==0)
    {
        printTab("Random move perturbation.");
        emili::Neighborhood* n = neigh();
        int num = number();
        tab_level++;
        oss.str(""); oss  << "number of moves per pertubation step " << num;
        printTab(oss.str().c_str());
        tab_level--;
        per = new emili::RandomMovePertubation(*n,num);
    }
    else if(strcmp(t,PERTUBATION_NOPER)==0)
    {
        printTab("No pertubation.");
        per = new emili::NoPertubation();
    }
    else if(strcmp(t,PERTUBATION_NRZ)==0)
    {
        int n = number();
        oss.str(""); oss  << "neh rz destruct/construct pertubation. number of job erased: "<<n;
        printTab(oss.str().c_str());
        per = new emili::pfsp::NRZPertubation(n,*istance);
    }else if(strcmp(t,PERTUBATION_VNRANDOM_MOVE)==0)
    {
        printTab("Random move perturbation." );
        tab_level++;
        int num = number();
        oss.str(""); oss  << "number of moves per pertubation step " << num << ".\n\t";
        printTab(oss.str().c_str());
        int iter = number();
        oss.str(""); oss  << "number of iteration before changing the neighborhood " << iter << ".\n\t";
        printTab(oss.str().c_str());
        nes.clear();
        tab_level--;
        neighs();
        per = new emili::VNRandomMovePertubation(nes,num,iter);
    }
    else if(strcmp(t,PERTUBATION_TMIIG)==0)
    {
        int num = number();
        int tsize = number();
        oss.str(""); oss  << "TMIIG pertubation. Number of job erased " << num << ". tabu list size " << tsize <<".\n\t";
        printTab(oss.str().c_str());
        per = new emili::pfsp::TMIIGPertubation(num,*istance,tsize);
    }
    else
    {
        std::cerr<< "'" << t << "' -> ERROR a pertubation criteria specification was expected! " << std::endl;
        prs::info();
        exit(-1);
    }
    tab_level --;
    return per;
}

emili::Acceptance* prs::ParamsParser::acc()
{
    tab_level++;
    char* t = nextToken();
    check(t,"ACCEPTANCE CRITERIA EXPECTED!");
    emili::Acceptance* acc;
    std::ostringstream oss;
    if(strcmp(t,ACCEPTANCE_SOA)==0)
    {
        float n = decimal();
        oss.str(""); oss  << "soa metropolis like acceptance. temperature : "<<n;
        printTab(oss.str().c_str());
        acc = new  emili::pfsp::SOAacceptance(n);
    }
    else if(strcmp(t,ACCEPTANCE_TEST)==0)
    {
        int n = number();
        oss.str(""); oss  << "Probabilistic Acceptance. improving solution accepted"<<n<<" % of the time";
        printTab(oss.str().c_str());
        acc = new  emili::pfsp::PfspTestAcceptance(*istance,n);
    }
    else  if(strcmp(t,ACCEPTANCE_METRO)==0)
    {
        float n = decimal();
        oss.str(""); oss  << "metropolis acceptance. temperature : "<<n;
        printTab(oss.str().c_str());
        acc = new  emili::MetropolisAcceptance(n);
    }
    else  if(strcmp(t,ACCEPTANCE_ALWAYS)==0)
    {
        char* t1 = nextToken();

        emili::accept_candidates accc;
        if(strcmp(t1,INTENSIFY)==0)
        {
            accc = emili::ACC_INTENSIFICATION;
        }
        else if(strcmp(t1,DIVERSIFY)==0)
        {
            accc = emili::ACC_DIVERSIFICATION;
        }
        else
        {
            std::cerr<< "'" << t1 << "' -> ERROR " << INTENSIFY << " or " << DIVERSIFY <<" was expected! " << std::endl;
            prs::info();
        exit(-1);
        }
        oss.str(""); oss  << "Acceptance always "<< t1;
        printTab(oss.str().c_str());
        acc = new  emili::AlwaysAccept(accc);
    }
    else  if(strcmp(t,ACCEPTANCE_IMPROVE)==0)
    {

        printTab( "improve acceptance");

        acc = new  emili::ImproveAccept();
    }
    else  if(strcmp(t,ACCEPTANCE_SA_METRO)==0)
    {
        float start = decimal();
        float end = decimal();
        float ratio = decimal();
        oss.str(""); oss  << "metropolis acceptance. start ,end , ratio : "<< start << ", "<< end << "," << ratio;
        printTab(oss.str().c_str());
        acc = new  emili::Metropolis(start,end,ratio);
    }
    else  if(strcmp(t,ACCEPTANCE_PMETRO)==0)
    {
        float start = decimal();
        float end = decimal();
        float ratio = decimal();
        int iterations = number();
        oss.str(""); oss  << "metropolis acceptance. start ,end , ratio, frequence : "<< start << ", "<< end << "," << ratio <<","<< iterations;
        printTab(oss.str().c_str());
        acc = new  emili::Metropolis(start,end,ratio,iterations);
    }
    else if(strcmp(t,ACCEPTANCE_SA)==0)
    {
        float start = decimal();
        float end = decimal();
        float ratio = decimal();
        int iterations = number();
        float alpha = decimal();
        oss.str(""); oss  << "metropolis acceptance. start ,end , ratio, frequence, alpha : "<< start << ", "<< end << "," << ratio <<","<< iterations << "," << alpha;
        printTab(oss.str().c_str());
        acc = new  emili::Metropolis(start,end,ratio,iterations,alpha);
    }
    else if(strcmp(t,ACCEPTANCE_TMIIG)==0)
    {
        float t0 = decimal();
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
    }else if(strcmp(t,ACCEPTANCE_IMPROVE_PLATEAU)==0)
    {
        int plateau_steps = number();
        int threshold = number();
        oss.str(""); oss  << "Accept a diversification solution if it improves on the intensification otherwise it will accept "<< plateau_steps << " non improving steps once it reaches the threshold of " << threshold;
        printTab(oss.str().c_str());
        acc = new  emili::AcceptPlateau(plateau_steps,threshold);
    }
    else
    {
        std::cerr<< "'" << t << "' -> ERROR an acceptance criteria specification was expected! " << std::endl;
        prs::info();
        exit(-1);
    }
    tab_level--;
    return acc;
}



emili::LocalSearch* prs::ParamsParser::ig()
{

    emili::Constructor* ls = new emili::pfsp::NEHSlackConstructor(*istance);//search();
    ils_time = ilstime();
    if(ils_time<=0)
    {
        std::cerr <<"ERROR for ils a time has to be provided"<< std::endl;
        prs::info();
        exit(-1);
    }
    //ils_time = ilstime();
    emili::WhileTrueTermination* pft = new emili::WhileTrueTermination;
    emili::Destructor* prsp = new emili::pfsp::PfspDestructorTest(*istance);
    //emili::AcceptanceCriteria* tac = new emili::pfsp::PfspTestAcceptance(istance);
    emili::Acceptance* tac = new emili::MetropolisAcceptance(10000);
    emili::LocalSearch* iils = new emili::IteratedGreedy(*ls,*pft,*prsp,*tac);
    iils->setSearchTime(ils_time);
    return iils;
}

emili::LocalSearch* prs::ParamsParser::gvns()
{
    emili::InitialSolution* is = init();
    emili::Perturbation* p1 = per();
    emili::Perturbation* p2 = per();
    emili::pfsp::GVNS_innerloop* gvi = new emili::pfsp::GVNS_innerloop(*is);
    emili::Perturbation* p3 = per();
    emili::Perturbation* p4 = per();

    std::vector< emili::Perturbation* > p;
    p.push_back(p1);
    p.push_back(p2);

    std::vector< emili::Perturbation* > pl;
    pl.push_back(p3);
    pl.push_back(p4);

    emili::GVNS* gg = new emili::GVNS(*gvi,p);
    return new emili::GVNS(*gvi,pl);
}


int prs::ParamsParser::ilstime()
{

    char* t = nextToken();
    if(t!=nullptr)
    {
        if(strcmp(t,IT)==0)
        {
            int n = number();
            std::ostringstream oss;
            oss << "ILS time secs : " << n;
            printTab(oss.str().c_str());

            return n;
        }
        else
        {
            currentToken--;
        }

    }

    return DEFAULT_IT;
}

emili::TabuSearch* prs::ParamsParser::tparams()
{
    params();
    tmem = tmemory(ne);
    //std::pair<int,int > tset = tsettings();
   // tm->setTabuTenure(tset.first);
    //emili::pfsp::PfspTerminationIterations* ptc = new emili::pfsp::PfspTerminationIterations(tset.second);
    return new emili::TabuSearch(*in,*te,*ne,*tmem);
}

emili::TabuMemory* prs::ParamsParser::tmemory(emili::pfsp::PfspNeighborhood* n)
{
    tab_level++;
    char* t = nextToken();
    check(t,"TABU MEMORY PARAMETERS EXPECTED!");
    int ts = number();
    std::ostringstream oss;
    emili::TabuMemory* tmem;

    oss << "Tabu tenure size " << ts;
    printTab(oss.str().c_str());

    if(strcmp(t,TABU_MEMORY_MOVES)==0)
    {
        oss.str(""); oss << "USING MOVES\n\t";
        printTab(oss.str().c_str());
        tmem = new  emili::pfsp::PfspMovesMemory(ts , *n);
    }
    else if(strcmp(t,TABU_MEMORY_HASHES)==0)
    {
        oss.str(""); oss << "USING HASHES\n\t";
        printTab(oss.str().c_str());
        tmem = new  emili::pfsp::PfspTabuHashMemory(ts);
    }
    else if(strcmp(t,TABU_MEMORY_SOLUTIONS)==0)
    {
        oss.str(""); oss << "USING FULL SOLUtiON\n\t";
        printTab(oss.str().c_str());
        tmem = new  emili::pfsp::PfspFullSolutionMemory(ts);
    }
    else if(strcmp(t,TABU_MEMORY_TSAB)==0)
    {
        oss.str(""); oss << "USING TSAB\n\t";
        printTab(oss.str().c_str());
        tmem = new  emili::pfsp::TSABMemory(ts , *n);
    }
    else if(strcmp(t,TABU_MEMORY_TSAB_TEST)==0)
    {
        oss.str(""); oss << "USING TSAB\n\t";
        printTab(oss.str().c_str());
        tmem = new  emili::pfsp::TSABtestMemory(ts , *n);
    }
    else if(strcmp(t,TABU_MEMORY_VALUE)==0)
    {
        oss.str(""); oss << "USING VALUE\n\t" ;
        printTab(oss.str().c_str());
        tmem = new  emili::pfsp::PfspTabuValueMemory(ts);
    }
    else
    {
        std::cerr<< "'" << t << "' -> ERROR a memory specification for the tabu search was expected! " << std::endl;
        prs::info();
        exit(-1);
    }
    tab_level--;
    return tmem;
}

std::pair<int,int> prs::ParamsParser::tsettings()
{

    int ts = ttsize();
    int ti = ttiter();
    return std::pair<int,int>(ts,ti);
}

int prs::ParamsParser::ttsize()
{
    char* t = nextToken();
    if(t!=nullptr)
    {
        if(strcmp(t,TS)==0)
        {
            std::cout << "Table tenure size : ";
            return number();
        }else
        {
            currentToken--;
        }

    }
    return DEFAULT_TS;
}

int prs::ParamsParser::ttiter()
{
    char* t = nextToken();
    if(t!=nullptr)
    {
        if(strcmp(t,TI)==0)
        {
            std::cout << "Number of tabu iterations : ";
            return number();
        }else
        {
            currentToken--;
        }

    }
    return DEFAULT_TI;
}

void prs::ParamsParser::params()
{
    in = init();
    te = term();
    ne = neigh();
}

emili::LocalSearch* prs::ParamsParser::vparams()
{
    tab_level++;
    char* t = nextToken();
    check(t,"TYPE OF SEARCH FOR VND MISSING!!!");
    in = init();
    te = term();
    neighs();
    emili::LocalSearch* ls;
    if(strcmp(t,FIRST)==0)
    {
        printTab("FIRST IMPROVEMENT VND");
        ls =  new emili::VNDSearch<emili::FirstImprovementSearch>(*in,*te,nes);
    }
    else if(strcmp(t,BEST)==0)
    {
       printTab("BEST IMPROVEMENT VND");
        ls =  new emili::VNDSearch<emili::BestImprovementSearch>(*in,*te,nes);
    }
    else
    {
        std::cerr<< "'" << t << "' -> ERROR a valid type of search must be specified (first,best) " << std::endl;

        prs::info();
        exit(-1);
    }
    tab_level--;
    return ls;
}

emili::InitialSolution* prs::ParamsParser::init()
{
    tab_level++;
    char* t = nextToken();
    check(t,"INITIAL SOLUTION GENERATOR EXPECTED!");
    std::ostringstream oss;
    emili::InitialSolution* init;
    if(strcmp(t,INITIAL_RANDOM)==0)
    {
        printTab("Random initial solution");
        init = new emili::pfsp::PfspRandomInitialSolution(*istance);
    }
    else if(strcmp(t,INITIAL_SLACK)==0)
    {
        printTab("SLACK initial solution");
        init = new emili::pfsp::PfspSlackInitialSolution(*istance);
    }else if(strcmp(t,INITIAL_WNSLACK)==0)
    {
        printTab( "NEH WSLACK initial solution");
        //init = new testIS(istance);
        init = new emili::pfsp::PfspNEHwslackInitialSolution(*istance);
    }
    else if(strcmp(t,INITIAL_LIT)==0)
        {
            printTab( "Less idle times initial solution");
            //return new testIS(istance);
            init = new emili::pfsp::LITSolution(*istance);
        }
    else if(strcmp(t,INITIAL_RZ)==0)
        {
            printTab( "rz initial solution");
            //return new testIS(istance);
            init = new emili::pfsp::RZSolution(*istance);
        }
    else if(strcmp(t,INITIAL_NRZ)==0)
        {
            printTab( "neh rz initial solution");
            //return new testIS(istance);
            init = new emili::pfsp::NeRZSolution(*istance);
        }
    else if(strcmp(t,INITIAL_NRZ2)==0)
        {
            printTab( "neh rz initial solution without improvement phase");
            //return new testIS(*istance);
            init = new emili::pfsp::NeRZ2Solution(*istance);
        }
    else if(strcmp(t,INITIAL_LR)==0)
        {
            int n = number();
            oss.str(""); oss << "LR initial solution with "<<n<<" starting sequences";
            printTab(oss.str().c_str());
            // testIS(*istance);
            init = new emili::pfsp::LRSolution(*istance,n);
        }
    else if(strcmp(t,INITIAL_NLR)==0)
        {
        int n = number();
        oss.str(""); oss << "NLR initial solution with "<<n<<" starting sequences";
        //return new testIS(*istance);printTab(oss.str().c_str());
        printTab(oss.str().c_str());
        init = new emili::pfsp::NLRSolution(*istance,n);
        }
    else if(strcmp(t,INITIAL_MNEH)==0)
        {
            printTab( "mneh initial solution");
            //return new testIS(istance);
            init = new emili::pfsp::MNEH(*istance);
        }
    else
    {
        std::cerr<< "'" << t << "' -> ERROR a initial solution generator specification was expected! (random,slack)" << std::endl;

        prs::info();
        exit(-1);
    }
    tab_level--;
    return init;
}

emili::Termination* prs::ParamsParser::term()
{
    tab_level++;
    char* t = nextToken();
    check(t,"TERMINATION CRITERIA EXPECTED!");
    emili::Termination* term;
    if(strcmp(t,TERMINATION_LOCMIN)==0)
    {
        printTab("Local minima termination");
        term = new emili::LocalMinimaTermination();
    }
    else if(strcmp(t,TERMINATION_WTRUE)==0)
    {
        printTab("While true termination");
        term = new emili::WhileTrueTermination();
    }
    else if(strcmp(t,TERMINATION_ITERA)==0)
    {

        int ti = number();
        std::ostringstream oss;
        oss << "Relaxed local minima termination. number of max iterations "<< ti;
        printTab(oss.str().c_str());
        term =  new emili::pfsp::PfspTerminationIterations(ti);
    }
    else if(strcmp(t,TERMINATION_SOA)==0)
    {
        printTab("Max iteration number termination");
        int ti = istance->getNjobs();
         ti = 2*(ti-1);
        term =  new emili::pfsp::SOAtermination(ti);
    }
    else if(strcmp(t,TERMINATION_TIME)==0)
    {

        float time = decimal();
        if(time==0){
            time = 1;
        }
        std::ostringstream oss;
        oss << "Timed termination. ratio: " << time;
        printTab(oss.str().c_str());
        term =  new emili::TimedTermination(time);
    }
    else if(strcmp(t,TERMINATION_MAXSTEPS)==0)
    {
        int steps = number();
        std::ostringstream oss;
        oss << "Max Steps termination. # steps: "<< steps;
        printTab(oss.str().c_str());
        term = new emili::MaxStepsTermination(steps);
    }
    else
    {
        std::cerr<< "'" << t << "' -> ERROR a termination criteria specification was expected! " << std::endl;
        prs::info();
        exit(-1);
    }
    tab_level--;
    return term;
}

emili::pfsp::PfspNeighborhood* prs::ParamsParser::neigh()
{
    tab_level++;
    char* t = nextToken();
    check(t,"NEIGHBORHOOD EXPECTED!");

    emili::pfsp::PfspNeighborhood* neigh;
    if(strcmp(t,NEIGHBORHOOD_INSERT)==0)
    {
        printTab( "Insert Neighborhood");
        neigh = new emili::pfsp::PfspInsertNeighborhood(*istance);
    }
    else  if(strcmp(t,NEIGHBORHOOD_FORW_INSERT)==0)
    {
        printTab( "Forward insert Neighborhood");
        neigh = new emili::pfsp::PfspForwardInsertNeighborhood(*istance);
    }
    else  if(strcmp(t,NEIGHBORHOOD_BACK_INSERT)==0)
    {
        printTab( "Backward Insert Neighborhood");
        neigh = new emili::pfsp::PfspBackwardInsertNeighborhood(*istance);
    }
    else if(strcmp(t,NEIGHBORHOOD_EXCHANGE)==0)
    {
        printTab( "Exchange neighborhood");
        neigh = new emili::pfsp::PfspExchangeNeighborhood(*istance);
    }
    else if(strcmp(t,NEIGHBORHOOD_TRANSPOSE)==0)
    {
        printTab( "Transpose neighborhood");
        neigh = new emili::pfsp::PfspTransposeNeighborhood(*istance);
    }
    else if(strcmp(t,NEIGHBORHOOD_TWO_INSERT)==0)
    {
        printTab( "Two insert neighborhood");
        neigh = new emili::pfsp::PfspTwoInsertNeighborhood(*istance);
    }
    else if(strcmp(t,NEIGHBORHOOD_XTRANSPOSE)==0)
    {
        printTab( "XTranspose neighborhood");
        neigh = new emili::pfsp::XTransposeNeighborhood(*istance);
    }
    else if(strcmp(t,NEIGHBORHOOD_TA_INSERT)==0)
    {
        printTab( "Insert with Taillard Acceleration");
        neigh = new emili::pfsp::TaillardAcceleratedInsertNeighborhood(*istance);
    }
    else if(strcmp(t,NEIGHBORHOOD_TAx_INSERT)==0)
    {
        printTab( "Insert with Taillard Acceleration");
        neigh = new emili::pfsp::TAxInsertNeighborhood(*istance);
    }
    else if(strcmp(t,NEIGHBORHOOD_NITA_INSERT)==0)
    {
        printTab( "Insert with Taillard Acceleration for no idle make span ");
        neigh = new emili::pfsp::NoIdleAcceleratedInsertNeighborhood(*istance);
    }
    else
    {
        std::cerr<< "'" << t << "' -> ERROR a neighborhood specification was expected! " << std::endl;
        prs::info();
        exit(-1);
    }
    tab_level--;
    return neigh;
}

emili::pfsp::PfspNeighborhood* prs::ParamsParser::neighV()
{
    tab_level++;
    char* t = nextToken();
    check(t,"NEIGHBORHOOD EXPECTED!");
    emili::pfsp::PfspNeighborhood* neigh;
    if(strcmp(t,NEIGHBORHOOD_INSERT)==0)
    {
        printTab( "Insert Neighborhood");
        neigh = new emili::pfsp::PfspInsertNeighborhood(*istance);
    }
    else  if(strcmp(t,NEIGHBORHOOD_FORW_INSERT)==0)
    {
        printTab( "Forward insert Neighborhood");
        neigh = new emili::pfsp::PfspForwardInsertNeighborhood(*istance);
    }
    else  if(strcmp(t,NEIGHBORHOOD_BACK_INSERT)==0)
    {
        printTab( "Backward Insert Neighborhood");
        neigh = new emili::pfsp::PfspBackwardInsertNeighborhood(*istance);
    }
    else if(strcmp(t,NEIGHBORHOOD_EXCHANGE)==0)
    {
        printTab( "Exchange neighborhood");
        neigh = new emili::pfsp::PfspExchangeNeighborhood(*istance);
    }
    else if(strcmp(t,NEIGHBORHOOD_TRANSPOSE)==0)
    {
        printTab( "Transpose neighborhood");
        neigh = new emili::pfsp::PfspTransposeNeighborhood(*istance);
    }
    else if(strcmp(t,NEIGHBORHOOD_TWO_INSERT)==0)
    {
        printTab( "Two insert neighborhood");
        neigh = new emili::pfsp::PfspTwoInsertNeighborhood(*istance);
    }
    else if(strcmp(t,NEIGHBORHOOD_XTRANSPOSE)==0)
    {
        printTab( "XTranspose neighborhood");
        neigh = new emili::pfsp::XTransposeNeighborhood(*istance);
    }
    else if(strcmp(t,NEIGHBORHOOD_TA_INSERT)==0)
    {
        printTab( "Insert with Taillard Acceleration");
        neigh = new emili::pfsp::TaillardAcceleratedInsertNeighborhood(*istance);
    }
    else if(strcmp(t,NEIGHBORHOOD_TAx_INSERT)==0)
    {
        printTab( "Insert with Taillard Acceleration");
        neigh = new emili::pfsp::TAxInsertNeighborhood(*istance);
    }
    else if(strcmp(t,NEIGHBORHOOD_NITA_INSERT)==0)
    {
        printTab( "Insert with Taillard Acceleration for no idle make span ");
        neigh = new emili::pfsp::NoIdleAcceleratedInsertNeighborhood(*istance);
    }
    else
    {
	currentToken--;
        neigh =  nullptr;
    }
    tab_level--;
    return neigh;
}
void prs::ParamsParser::neighs()
{
    std::vector<emili::Neighborhood*> vnds;
    vnds.push_back(neigh());
    nes = vnds;
    neighs1();
}

void prs::ParamsParser::neighs1()
{
    char* t = nextToken();
    if(t != nullptr )
    {
        currentToken--;
	emili::Neighborhood* n = neighV();
	if(n!=nullptr)
	{
           nes.push_back(n);
           neighs1();
	}
    }
}

int prs::ParamsParser::number()
{
    char* t = nextToken();
    check(t,"A NUMBER WAS EXPECTED!");
    int k = atoi(t);
   // std::cout << k << "\n\t";
    return k;
}

float prs::ParamsParser::decimal()
{
    char* t = nextToken();
    check(t,"A DECIMAL NUMBER WAS EXPECTED!");
    float k = atof(t);
    //std::cout << k << "\n\t";
    return k;
}


void prs::ParamsParser::problem()
{
    PfspInstance i;

    if(i.readDataFromFile(tokens[1]))
    {
        char* t = nextToken();
        problem_type = t;
        istance = instantiateProblem(t, i);
    return;
    }
        info();
        exit(-1);
}

char* prs::ParamsParser::nextToken()
{
    if(currentToken < numberOfTokens)
    {
        char* token = tokens[currentToken];
        currentToken++;
        return token;
    }
    else
    {
        return nullptr;
    }
}

emili::LocalSearch* prs::ParamsParser::parseParams()
{
    problem();
    emili::LocalSearch* local = eparams();
    std::cout << "------" << std::endl;
 return local;
}

