#include "paramsparser.h"
#include <cstdlib>
#include <cstdio>
#include <ctime>
#include <cstring>

#define IG "ig"
#define ILS "ils"
#define TABU "tabu"
#define FIRST "first"
#define BEST "best"
#define VND "vnd"
#define IT "-it"
#define MOVES "move"
#define HASHES "hash"
#define SOLUTIONS "solution"
#define TS "-ts"
#define TI "-ti"
#define IN "-in"
#define TE "-te"
#define NE "-ne"
#define NS "-ns"
#define RANDOM "random"
#define SLACK "slack"
#define LOCMIN "locmin"
#define ITERA "iteration"
#define MAXSTEPS "maxstep"
#define TIME "time"
#define WNSLACK "nwslack"
#define WTRUE "true"
#define INSERT "insert"
#define BACK_INSERT "binsert"
#define FORW_INSERT "finsert"
#define TRANSPOSE "transpose"
#define EXCHANGE "exchange"
#define RNDSEED "rnds"
#define PRT_RND "randpert"
#define ACC_PROB "prob"
#define ACC_METRO "metropolis"
#define SOA_PER "soaper"
#define TEST_PER "testper"
#define TEST_ACC "testacc"
#define SOA_ACC "soaacc"
#define SOA_TER "soater"
#define DEFAULT_TS 10
#define DEFAULT_TI 10
#define DEFAULT_IT -10


void prs::emili()
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
    std::cout << "EMILI INSTANCE_FILE_PATH <LOCAL_SEARCH | ITERATED_LOCAL_SEARCH | TABU_SEARCH | VND_SEARCH> [rnds seed]" << std::endl;
    std::cout << std::endl;
    std::cout << "LOCAL_SEARCH          = SEARCH_TYPE INITIAL_SOLUTION TERMINATION NEIGHBORHOOD" << std::endl;
    std::cout << "ITERATED_LOCAL_SEARCH = ils LOCAL_SEARCH TERMINATION PERTUBATION ACCEPTANCE -it seconds" << std::endl;
    std::cout << "TABU_SEARCH           = tabu INITIAL_SOLUTION TERMINATION NEIGHBORHOOD TABU_MEMORY" << std::endl;
    std::cout << "VND_SEARCH            = vnd < first | best > INITIAL_SOLUTION TERMINATION NEIGHBORHOOD1 NEIGHBORHOOD2 ... NEIGHBORHOODn" << std::endl;
    std::cout << "SEARCH_TYPE           = first | best | tabu | vnd | ils" << std::endl;
    std::cout << "INITIAL_SOLUTION      = random | slack | nwslack " << std::endl;
    std::cout << "TERMINATION           = true | time int | locmin | soater | iteration int | maxsteps int" << std::endl;
    std::cout << "NEIGHBORHOOD          = transpose | exchange | insert | binsert | finsert" << std::endl;
    std::cout << "PERTUBATION           = soaper int | testper " << std::endl;
    std::cout << "ACCEPTANCE            = soaacc float | testacc int | metropolis float" << std::endl;
    std::cout << "TABU_MEMORY           = move int | hash int | solution int" << std::endl;
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

emili::pfsp::PfspNeighborhood* ne;
emili::InitialSolution* in;
emili::Termination* te;
emili::TabuMemory* tmem;
emili::Termination* ilt;

std::vector< emili::Neighborhood*> nes;

emili::LocalSearch* prs::ParamsParser::eparams()
{

    char* t = nextToken();

    check(t,"SEARCH TYPE AND PARAMETERS MISSING!!!");
    emili::LocalSearch* ls;
    if(strcmp(t,ILS)==0)
    {
        std::cout << "ILS \n\t";
        ls = ils();

    }
    else if(strcmp(t,IG)==0)
    {
        std::cout << "Iterated Greedy \n\t";
        ls = ig();

    }
    else
    {
        currentToken--;
        std::cout << "LocalSearch \n\t";
        ls = search();

    ils_time = ilstime();
    ls->setSearchTime(ils_time);
    }
    int seed = getSeed();
    std::cout << "\tRANDOM SEED " << seed<< "\n\t" ;
    emili::initializeRandom(seed);
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
    char* t = nextToken();
    check(t,"SEARCH PARAMETERS MISSING!!!");
    emili::LocalSearch* ls;
    if(strcmp(t,TABU)==0)
    {
        std::cout << "TABU SEARCH\n\t";
        ls = tparams();
    }
    else if(strcmp(t,FIRST)==0)
    {
        std::cout << "FIRST IMPROVEMENT \n\t";
        params();
        ls =  new emili::FirstImprovementSearch(*in,*te,*ne);
    }
    else if(strcmp(t,BEST)==0)
    {
        std::cout << "BEST IMPROVEMENT\n\t";
        params();
        ls =  new emili::BestImprovementSearch(*in,*te,*ne);
    }
    else if(strcmp(t,VND)==0)
    {
        std::cout << "VND SEARCH\n\t";
        ls = vparams();
    }
    else
    {
        std::cerr<< "'" << t << "' -> ERROR a search definition was expected! " << std::endl;
        exit(-1);
    }
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
    emili::AcceptanceCriteria* tac = acc();//new emili::pfsp::SOAacceptance(1.2f);
    emili::LocalSearch* iils = new emili::IteratedLocalSearch(*ls,*pft,*prsp,*tac);
    ils_time = ilstime();
    if(ils_time<=0)
    {
        std::cerr <<"ERROR for ils a time has to be provided"<< std::endl;
        exit(-1);
    }
    iils->setSearchTime(ils_time);
    return iils;
}

emili::Perturbation* prs::ParamsParser::per()
{
    char* t = nextToken();
    check(t,"PERTUBATION CRITERIA EXPECTED!");
    if(strcmp(t,SOA_PER)==0)
    {
        int n = number();
        std::cout << "wslack destruct/construct pertubation. number of job erased: "<<n<<"\n\t";

        return new emili::pfsp::SOAPerturbation(n,istance);
    }
    else if(strcmp(t,TEST_PER)==0)
    {
        std::cout << "Random swap test pertubation. \n\t";
        return new emili::pfsp::PfspRandomSwapPertub(istance);
    }
    else
    {
        std::cerr<< "'" << t << "' -> ERROR a pertubation criteria specification was expected! " << std::endl;
        exit(-1);
    }
}

emili::AcceptanceCriteria* prs::ParamsParser::acc()
{
    char* t = nextToken();
    check(t,"ACCEPTANCE CRITERIA EXPECTED!");
    if(strcmp(t,SOA_ACC)==0)
    {
        float n = decimal();
        std::cout << "soa metropolis like acceptance. temperature : "<<n<<"\n\t";

        return new emili::pfsp::SOAacceptance(n);
    }
    else if(strcmp(t,TEST_ACC)==0)
    {
        int n = number();
        std::cout << "Random swap test pertubation. improving solution accepted"<<n<<"% of the time.\n\t";
        return new emili::pfsp::PfspTestAcceptance(istance,n);
    }
    else  if(strcmp(t,ACC_METRO)==0)
    {
        float n = decimal();
        std::cout << "metropolis acceptance. temperature : "<<n<<"\n\t";

        return new emili::MetropolisAcceptance(n);
    }
    else
    {
        std::cerr<< "'" << t << "' -> ERROR an acceptance criteria specification was expected! " << std::endl;
        exit(-1);
    }
}



emili::LocalSearch* prs::ParamsParser::ig()
{

    emili::Constructor* ls = new emili::pfsp::NEHSlackConstructor(istance);//search();
    ils_time = ilstime();
    if(ils_time<=0)
    {
        std::cerr <<"ERROR for ils a time has to be provided"<< std::endl;
        exit(-1);
    }
    //ils_time = ilstime();
    emili::WhileTrueTermination* pft = new emili::WhileTrueTermination;
    emili::Destructor* prsp = new emili::pfsp::PfspDestructorTest(istance);
    //emili::AcceptanceCriteria* tac = new emili::pfsp::PfspTestAcceptance(istance);
    emili::AcceptanceCriteria* tac = new emili::MetropolisAcceptance(10000);
    emili::LocalSearch* iils = new emili::IteratedGreedy(*ls,*pft,*prsp,*tac);
    iils->setSearchTime(ils_time);
    return iils;
}


int prs::ParamsParser::ilstime()
{
    char* t = nextToken();
    if(t!=nullptr)
    {
        if(strcmp(t,IT)==0)
        {
            int n = number();
            std::cout << "ILS time secs : " << n << std::endl;
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
    char* t = nextToken();
    check(t,"TABU MEMORY PARAMETERS EXPECTED!");
    int ts = ttsize();
    if(strcmp(t,MOVES)==0)
    {
        std::cout << "USING MOVES\n\t";
        return new emili::pfsp::PfspMovesMemory(ts , *n);
    }
    else if(strcmp(t,HASHES)==0)
    {
        std::cout << "USING HASHES\n\t";
        return new emili::pfsp::PfspTabuHashMemory(ts);
    }
    else if(strcmp(t,SOLUTIONS)==0)
    {
        std::cout << "USING FULL SOLUtiON\n\t";
        return new emili::pfsp::PfspFullSolutionMemory(ts);
    }
    else
    {
        std::cerr<< "'" << t << "' -> ERROR a memory specification for the tabu search was expected! " << std::endl;
        exit(-1);
    }
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
    char* t = nextToken();
    check(t,"TYPE OF SEARCH FOR VND MISSING!!!");
    in = init();
    te = term();
    neighs();
    emili::LocalSearch* ls;
    if(strcmp(t,FIRST)==0)
    {
        std::cout << "FIRST IMPROVEMENT VND \n\t";
        ls =  new emili::VNDSearch<emili::FirstImprovementSearch>(*in,*te,nes);
    }
    else if(strcmp(t,BEST)==0)
    {
        std::cout << "BEST IMPROVEMENT VND \n\t";
        ls =  new emili::VNDSearch<emili::BestImprovementSearch>(*in,*te,nes);
    }
    else
    {
        std::cerr<< "'" << t << "' -> ERROR a valid type of search must be specified (first,best) " << std::endl;

        exit(-1);
    }
    return ls;
}

emili::InitialSolution* prs::ParamsParser::init()
{
    char* t = nextToken();
    check(t,"INITIAL SOLUTION GENERATOR EXPECTED!");
    if(strcmp(t,RANDOM)==0)
    {
        std::cout << "Random initial solution\n\t";
        return new emili::pfsp::PfspRandomInitialSolution(istance);
    }
    else if(strcmp(t,SLACK)==0)
    {
        std::cout << "SLACK initial solution\n\t";
        return new emili::pfsp::PfspSlackInitialSolution(istance);
    }else if(strcmp(t,WNSLACK)==0)
    {
        std::cout << "NEH WSLACK initial solution\n\t";
        //return new testIS(istance);
        return new emili::pfsp::PfspNEHwslackInitialSolution(istance);
    }
    else
    {
        std::cerr<< "'" << t << "' -> ERROR a initial solution generator specification was expected! (random,slack)" << std::endl;

        exit(-1);
    }
}

emili::Termination* prs::ParamsParser::term()
{
    char* t = nextToken();
    check(t,"TERMINATION CRITERIA EXPECTED!");
    if(strcmp(t,LOCMIN)==0)
    {
        std::cout << "Local minima termination\n\t";
        return new emili::LocalMinimaTermination();
    }
    else if(strcmp(t,WTRUE)==0)
    {
        std::cout << "While true termination\n\t";

        return new emili::WhileTrueTermination();
    }
    else if(strcmp(t,ITERA)==0)
    {

        int ti = number();
        std::cout << "Relaxed local minima termination. number of max iterations "<< ti <<"\n\t";
        return new emili::pfsp::PfspTerminationIterations(ti);
    }
    else if(strcmp(t,SOA_TER)==0)
    {
        std::cout << "Max iteration number termination\n\t";
        int ti = istance.getNjobs();
         ti = 2*(ti-1);
        return new emili::pfsp::SOAtermination(ti);
    }
    else if(strcmp(t,TIME)==0)
    {

        int time = number();
        std::cout << "Timed termination. secs: " << time << "\n\t";
        return new emili::TimedTermination(time);
    }
    else if(strcmp(t,MAXSTEPS)==0)
    {
        int steps = number();
        std::cout << "Max Steps termination. # steps: "<< steps << "\n\t";
        return new emili::MaxStepsTermination(steps);
    }
    else
    {
        std::cerr<< "'" << t << "' -> ERROR a termination criteria specification was expected! " << std::endl;
        exit(-1);
    }
}

emili::pfsp::PfspNeighborhood* prs::ParamsParser::neigh()
{
    char* t = nextToken();
    check(t,"NEIGHBORHOOD EXPECTED!");
    if(strcmp(t,INSERT)==0)
    {
        std::cout << "Insert Neighborhood\n\t";
        return new emili::pfsp::PfspInsertNeighborhood(istance);
    }
    else  if(strcmp(t,FORW_INSERT)==0)
    {
        std::cout << "Forward insert Neighborhood\n\t";
        return new emili::pfsp::PfspForwardInsertNeighborhood(istance);
    }
    else  if(strcmp(t,BACK_INSERT)==0)
    {
        std::cout << "Backward Insert Neighborhood\n\t";
        return new emili::pfsp::PfspBackwardInsertNeighborhood(istance);
    }
    else if(strcmp(t,EXCHANGE)==0)
    {
        std::cout << "Exchange neighborhood\n\t";
        return new emili::pfsp::PfspExchangeNeighborhood(istance);
    }
    else if(strcmp(t,TRANSPOSE)==0)
    {
        std::cout << "Transpose neighborhood\n\t";
        return new emili::pfsp::PfspTransposeNeighborhood(istance);
    }
    else
    {
        std::cerr<< "'" << t << "' -> ERROR a neighborhood specification was expected! " << std::endl;
        exit(-1);
    }
}

emili::pfsp::PfspNeighborhood* prs::ParamsParser::neighV()
{
    char* t = nextToken();
    check(t,"NEIGHBORHOOD EXPECTED!");
    if(strcmp(t,INSERT)==0)
    {
        std::cout << "Insert Neighborhood\n\t";
        return new emili::pfsp::PfspInsertNeighborhood(istance);
    }
    else if(strcmp(t,EXCHANGE)==0)
    {
        std::cout << "Exchange neighborhood\n\t";
        return new emili::pfsp::PfspExchangeNeighborhood(istance);
    }
    else if(strcmp(t,TRANSPOSE)==0)
    {
        std::cout << "Transpose neighborhood\n\t";
        return new emili::pfsp::PfspTransposeNeighborhood(istance);
    }
    else
    {
	currentToken--;
        return nullptr;
    }
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
    emili::LocalSearch* local = eparams();
    std::cout << "------" << std::endl;
 return local;
}

