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
#define TIME "time"
#define WNSLACK "nwslack"
#define WTRUE "true"
#define INSERT "insert"
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

class testIS: public emili::InitialSolution
{
public:
    testIS(emili::Problem& problem_instance):emili::InitialSolution(problem_instance){}
    /*
        The generated solution must be a valid solution for the problem with
        the appropriate data structures used for the implemented solution.
    */
    virtual emili::Solution* generateSolution()
    {
        //int arr[] = {0,39,83,45,12,86,4,43,60,31,87,24,53,18,90,37,70,93,41,82,73,38,79,9,6,1,5,15,44,61,2,27,97,77,98,47,59,58,26,63,80,46,65,89,68,78,51,25,55,75,66,17,84,28,36,33,99,48,72,54,16,3,30,96,64,95,20,50,21,92,8,91,34,14,11,42,74,56,29,76,94,22,7,100,23,32,13,62,52,85,67,19,69,35,71,49,10,57,88,40,81};
        int arr[] = {0,83,39,45,12,86,4,43,60,31,87,24,53,18,90,37,70,93,41,82,73,38,79,9,6,1,5,15,44,61,2,27,97,77,98,47,59,58,26,63,80,46,65,89,68,78,51,25,55,75,66,17,84,28,36,33,99,48,72,54,16,3,30,96,64,95,20,50,21,92,8,91,34,14,11,42,74,56,29,76,94,22,7,100,23,32,13,62,52,85,67,19,69,35,71,49,10,57,88,40,81};
        std::vector< int > v(arr, arr + sizeof(arr) / sizeof(int) );
        /*
        v.push_back(0);
        v.push_back(82);
        v.push_back(38);
        v.push_back(44);
        v.push_back(11);
        v.push_back(85);
        v.push_back(3);
        v.push_back(42);
        v.push_back(59);
        v.push_back(30);
        v.push_back(86);
        v.push_back(23);
        v.push_back(52);
        v.push_back(17);
        v.push_back(89);
        v.push_back(36);
        v.push_back(69);
        v.push_back(92);
        v.push_back(40);
        v.push_back(81);
        v.push_back(72);
        v.push_back(37);
        v.push_back(78);
        v.push_back(8 );
        v.push_back(5);
        v.push_back(0 );
        v.push_back(4);
        v.push_back(14);
        v.push_back(43);
        v.push_back(60);
        v.push_back(1);
        v.push_back(26);
        v.push_back(96);
        v.push_back(76);
        v.push_back(97);
        v.push_back(46);
        v.push_back(58);
        v.push_back(57);
        v.push_back(25);
        v.push_back(62);
        v.push_back(79);
        v.push_back(45);
        v.push_back(64);
        v.push_back(88);
        v.push_back(67);
        v.push_back(77);
        v.push_back(50);
        v.push_back(24);
        v.push_back(54);
        v.push_back(74);
        v.push_back(65);
        v.push_back(16);
        v.push_back(83);
        v.push_back(27);
        v.push_back(35);
        v.push_back(32);
        v.push_back(98);
        v.push_back(47);
        v.push_back(71);
        v.push_back(53);
        v.push_back(15);
        v.push_back(2);
        v.push_back(29);
        v.push_back(95);
        v.push_back(63);
        v.push_back(94);
        v.push_back(19);
        v.push_back(49);
        v.push_back(20);
        v.push_back(91);
        v.push_back(7);
        v.push_back(90);
        v.push_back(33);
        v.push_back(13);
        v.push_back(10);
        v.push_back(41);
        v.push_back(73);
        v.push_back(55);
        v.push_back(28);
        v.push_back(75);
        v.push_back(93);
        v.push_back(21);
        v.push_back(6);
        v.push_back(99);
        v.push_back(22);
        v.push_back(31);
        v.push_back(12);
        v.push_back(61);
        v.push_back(51);
        v.push_back(84);
        v.push_back(66);
        v.push_back(18);
        v.push_back(68);
        v.push_back(34);
        v.push_back(70);
        v.push_back(48);
        v.push_back(9);
        v.push_back(56);
        v.push_back(87);
        v.push_back(39);
        v.push_back(80);
        for(std::vector< int > ::const_iterator it = v.begin()+1;it!=v.end();++it)
        {
            *it++;
        }*/
        emili::Solution* test = new emili::pfsp::PermutationFlowShopSolution(v);

        instance.evaluateSolution(*test);
        return test;
    }

    virtual emili::Solution* generateEmptySolution()
    {
        return generateSolution();
    }
};


void info()
{
    std::cout << " syntax for local search -> EMILI instancefile search_type intial_solution termination neighborhood [parameters]" << std::endl;
    std::cout << " syntax for iterated local search -> EMILI instancefile ils search_type intial_solution termination neighborhood ilstermination perturbation acceptance -it seconds" << std::endl;
    std::cout << "search_type = first | best | tabu | vnd | ils" << std::endl;
    std::cout << "intial_solution = random | slack | nwslack " << std::endl;
    std::cout << "termination = true | time | locmin | soater | iteration [-it number of iterations]" << std::endl;
    std::cout << "neighborhood = transpose | exchange | insert " << std::endl;
    std::cout << "pertubaiton = soaper int | testper " << std::endl;
    std::cout << "acceptance = soaacc float | testacc int | metropolis float" << std::endl;
   // std::cout << " syntax->EMILI instancefile search_type intial_solution termination neighborhood" << std::endl;
}

void check(char* t,const char* message)
{
    if(t==nullptr)
    {
        std::cerr <<"PARSING ERROR"<< message << std::endl;
        info();
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
    std::cout << "RANDOM SEED " << seed<< "\n\t" ;
    emili::initializeRandom(seed);
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
            std::cout << "ILS time secs : ";
            return number();
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
        std::cout << "Relaxed local minima termination\n\t";
        int ti = ttiter();
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
        std::cout << "Timed termination\n\t";
        int time = number();
        return new emili::TimedTermination(time);
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

