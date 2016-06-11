#include "irpparser.h"
#include <cstdlib>
#include <cstdio>
#include <ctime>
#include <cstring>
#include <iostream>
#include <sstream>
#include <algorithm>

/* Algos */
#define IG "ig"
#define ILS "ils"
#define TABU "tabu"
#define FIRST "first"
#define BEST "best"
#define VND "vnd"
#define TEST_INIT "stin"

#define IRP_TWOEXCHANGE "twoExc"
#define IRP_SHIFTTWOEXCHANGE "shiftTwoExc"
#define IRP_SHIFTINSERT "shiftIns"
#define IRP_SHIFTREMOVE "shiftRem"
#define IRP_REFUEL "ref"

/* initial solution heuristics */
#define INITIAL_GREEDY "greedy"
#define INITIAL_GREEDY_RANDOMIZED "greedyRandomized"
#define INITIAL_GRASP "grasp"

/* Termination criteria*/
#define TERMINATION_MAXSTEPS "maxstep"
#define TERMINATION_TIME "time"
#define TERMINATION_LOCMIN "locmin"
#define TERMINATION_ITERA "iteration"
#define TERMINATION_WTRUE "true"
#define TERMINATION_SOA "soater"

/* acceptance criteria*/
#define ACCEPTANCE_PROB "prob"
#define ACCEPTANCE_METRO "metropolis"
#define ACCEPTANCE_RS "rsacc"
#define ACCEPTANCE_PMETRO "pmetro"
#define ACCEPTANCE_TMIIG "tmiigacc"
#define ACCEPTANCE_IMPROVE_PLATEAU "implat"
#define ACCEPTANCE_TEST "testacc"
#define ACCEPTANCE_SOA "soaacc"
#define ACCEPTANCE_ALWAYS "always"
#define ACCEPTANCE_INTENSIFY "intensify"
#define ACCEPTANCE_DIVERSIFY "diversify"
#define ACCEPTANCE_IMPROVE "improve"
#define ACCEPTANCE_SA_METRO "sa_metropolis"
#define ACCEPTANCE_SA "saacc"

/* perturbations*/
#define PERTURBATION_RANDOM_MOVE "rndmv"



std::string prs::irp::IrpParser::info()
{
    ostringstream oss;
    oss << "Usage:\n\n";
    oss << "EMILI INSTANCE_FILE_PATH IRP_PROBLEM <LOCAL_SEARCH | ITERATED_LOCAL_SEARCH | TABU_SEARCH | VND_SEARCH> [rnds seed]\n\n";
    /*TODO INSERT HOW TO CALL YOUR PARSER HERE*/
    return oss.str();
}


emili::InitialSolution* in= nullptr;
emili::Termination* te= nullptr;
emili::TabuMemory* tmem= nullptr;
emili::Termination* ilt= nullptr;
emili::Neighborhood* ne = nullptr;
std::vector< emili::Neighborhood*> nes;

emili::LocalSearch* prs::irp::IrpParser::eparams(prs::TokenManager& tm)
{
    prs::incrementTabLevel();
    emili::LocalSearch* ls;
    if(tm.checkToken(ILS))
    {
        printTab("ILS");
        ls = ils(tm);
    }  
    else
    {             
        ls = search(tm);
    }

    prs::decrementTabLevel();
    return ls;
}



emili::LocalSearch* prs::irp::IrpParser::search(prs::TokenManager& tm)
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
    else if(tm.checkToken(TEST_INIT))
    {
        emili::InitialSolution* ini = init(tm);
        emili::Solution* s = ini->generateSolution();
        std::cout << s->getSolutionRepresentation() << std::endl;
        std::cout << s-> getSolutionValue() << std::endl;
        std::cerr << s-> getSolutionValue() << std::endl;
        exit(-1);
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

emili::LocalSearch* prs::irp::IrpParser::ils(prs::TokenManager& tm)
{
    /*Parse the inner local search*/
    emili::LocalSearch* ls = search(tm);
    /*Parse the termination*/
    emili::Termination* pft = term(tm);
    /*Parse the perturbation*/
    emili::Perturbation* prsp = per(tm);
    /*Parse the Acceptance*/
    emili::Acceptance* tac = acc(tm);
    /*Build the ils and returns it*/
    emili::LocalSearch* iils = new emili::IteratedLocalSearch(*ls,*pft,*prsp,*tac);   
    return iils;
}

emili::Perturbation* prs::irp::IrpParser::per(prs::TokenManager& tm)
{
    /*increment tab for nicer output*/
    prs::incrementTabLevel();
    emili::Perturbation* per;
    std::ostringstream oss;


    /*TODO HERE GOES THE CODE TO PARSE THE PERTURBATION*/
    if(tm.checkToken(PERTURBATION_RANDOM_MOVE))
        {
            printTab("Random move perturbation.");
            emili::Neighborhood* n = neigh(tm);
            int num = tm.getInteger();
            prs::incrementTabLevel();
            oss.str(""); oss  << "number of moves per perturbation step " << num;
            printTab(oss.str().c_str());
            prs::decrementTabLevel();
            per = new emili::RandomMovePertubation(*n,num);
        }
    else
        {
            std::cerr<< "'" << tm.peek() << "' -> ERROR a perturbation criteria specification was expected! " << std::endl;
            std::cout << info() << std::endl;
            exit(-1);
        }

    prs::decrementTabLevel();
    return per;
}

emili::Acceptance* prs::irp::IrpParser::acc(prs::TokenManager& tm)
{
    prs::incrementTabLevel();
    emili::Acceptance* acc;
    std::ostringstream oss;
    if(tm.checkToken(ACCEPTANCE_METRO))
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
        if(tm.checkToken(ACCEPTANCE_INTENSIFY))
        {
            accc = emili::ACC_INTENSIFICATION;
            t1 = ACCEPTANCE_INTENSIFY;
        }
        else if(tm.checkToken(ACCEPTANCE_DIVERSIFY))
        {
            t1 = ACCEPTANCE_DIVERSIFY;
            accc = emili::ACC_DIVERSIFICATION;
        }
        else
        {
            std::cerr<< "'" << *tm << "' -> ERROR " << ACCEPTANCE_INTENSIFY << " or " << ACCEPTANCE_DIVERSIFY <<" was expected! " << std::endl;
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

emili::BestTabuSearch* prs::irp::IrpParser::tparams(prs::TokenManager& tm)
{
    if(tm.checkToken(BEST))
    {
        params(tm);
        tmem = tmemory(ne,tm);
        return new emili::BestTabuSearch(*in,*te,*ne,*tmem);
    }
    else if(tm.checkToken(FIRST))
    {
        params(tm);
        tmem = tmemory(ne,tm);
        return new emili::FirstTabuSearch(*in,*te,*ne,*tmem);
    }
    else
    {
        std::cerr<< "'" << *tm << "' -> ERROR a pivotal rule (best or first) for the tabu search was expected! \n" << std::endl;
        std::cout << info() << std::endl;
        exit(-1);
    }
}

emili::TabuMemory* prs::irp::IrpParser::tmemory(emili::Neighborhood* n,prs::TokenManager& tm)
{
    prs::incrementTabLevel();


    emili::TabuMemory* tmem;

    /*TODO HERE GOES THE CODE TO INSTATIATE A TABUTENURE IMPLEMENTATION*/

    prs::decrementTabLevel();
    return tmem;
}

void prs::irp::IrpParser::params(prs::TokenManager& tm)
{
    in = init(tm);
    te = term(tm);
    ne = neigh(tm);
}

emili::LocalSearch* prs::irp::IrpParser::vparams(prs::TokenManager& tm)
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

emili::InitialSolution* prs::irp::IrpParser::init(prs::TokenManager& tm)
{
    prs::incrementTabLevel();    
    emili::InitialSolution* init;

    if(tm.checkToken(INITIAL_GREEDY))
    {
     /*TODO HERE THE CODE TO INSTATIATE A RANDOM INITIAL SOLUTION*/

        double randomFactor = tm.getDecimal();
        unsigned int urgencyPolicy = tm.getInteger();
        init = new emili::irp::GreedyInitialSolution(*this->istance, randomFactor, urgencyPolicy);

    }
    else if(tm.checkToken(INITIAL_GREEDY_RANDOMIZED))
    {
     /*TODO HERE THE CODE TO INSTATIATE A RANDOM INITIAL SOLUTION*/

        double randomFactor = tm.getDecimal();
        unsigned int urgencyPolicy = tm.getInteger();
        init = new emili::irp::GreedyRandomizedInitialSolution(*this->istance, randomFactor, urgencyPolicy);

    }
    else if(tm.checkToken(INITIAL_GRASP))
    {
     /*TODO HERE THE CODE TO INSTATIATE A RANDOM INITIAL SOLUTION*/

        double randomFactor = tm.getDecimal();
        unsigned int urgencyPolicy = tm.getInteger();
        init = new emili::irp::GRASP(*this->istance, randomFactor, urgencyPolicy);

    }
    /*TODO EXTEND THIS IF/ELSE STATEMENT TO INSTATIATE MORE INTIAL SOLUTION HEURISTICS*/
    else
    {
        std::cerr<< "'" << *tm << "' -> ERROR a initial solution generator specification was expected! (random,slack)" << std::endl;

        std::cout << info() << std::endl;
        exit(-1);
    }
    prs::decrementTabLevel();
    return init;
}

emili::Termination* prs::irp::IrpParser::term(prs::TokenManager& tm)
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

emili::Neighborhood* prs::irp::IrpParser::neigh(prs::TokenManager& tm)
{
    prs::incrementTabLevel();
    emili::Neighborhood* neigh;
/*    if(tm.checkToken(IRP_TWOEXCHANGE)){
        int piv = tm.getInteger();
        unsigned int ps = tm.getInteger();
        return new emili::irp::irpTwoExchangeNeighborhood(*istance, piv, ps);
    }
    else */if(tm.checkToken(IRP_SHIFTTWOEXCHANGE)){
        /*TODO HERE INSERT THE CODE FOR PARSING THE NEIGHBORHOODS*/
        double randomFactor = tm.getDecimal();
        unsigned int urgencyPolicy = tm.getInteger();
        double refuelFactor = tm.getDecimal();
        double deliveredQuantityFactor = tm.getDecimal();
        return new emili::irp::irpShiftTwoExchangeNeighborhood(*istance, randomFactor, urgencyPolicy, refuelFactor, deliveredQuantityFactor);
    }
    else if(tm.checkToken(IRP_SHIFTINSERT)){
        /*TODO HERE INSERT THE CODE FOR PARSING THE NEIGHBORHOODS*/
        double randomFactor = tm.getDecimal();
        unsigned int urgencyPolicy = tm.getInteger();
        double refuelFactor = tm.getDecimal();
        double deliveredQuantityFactor = tm.getDecimal();
        return new emili::irp::irpShiftInsertNeighborhood(*istance, randomFactor, urgencyPolicy, refuelFactor, deliveredQuantityFactor);
    }
    else if(tm.checkToken(IRP_SHIFTREMOVE)){
        /*TODO HERE INSERT THE CODE FOR PARSING THE NEIGHBORHOODS*/
        double randomFactor = tm.getDecimal();
        unsigned int urgencyPolicy = tm.getInteger();
        double refuelFactor = tm.getDecimal();
        double deliveredQuantityFactor = tm.getDecimal();
        return new emili::irp::irpShiftRemoveNeighborhood(*istance, randomFactor, urgencyPolicy, refuelFactor, deliveredQuantityFactor);
    }
    else if(tm.checkToken(IRP_REFUEL)){
        /*TODO HERE INSERT THE CODE FOR PARSING THE NEIGHBORHOODS*/
        double refuelStep = tm.getDecimal();
        double deliveredQuantityStep = tm.getDecimal();
        double randomFactor = tm.getDecimal();
        unsigned int urgencyPolicy = tm.getInteger();
        return new emili::irp::irpRefuelNeighborhood(*istance, refuelStep, deliveredQuantityStep,randomFactor, urgencyPolicy);
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

emili::Neighborhood* prs::irp::IrpParser::neighV(prs::TokenManager& tm)
{
    prs::incrementTabLevel();
    emili::Neighborhood* neigh;
/*    if(tm.checkToken(IRP_TWOEXCHANGE)){
        int piv = tm.getInteger();
        unsigned int ps = tm.getInteger();
        return new emili::irp::irpTwoExchangeNeighborhood(*istance, piv, ps);
    }
    else */if(tm.checkToken(IRP_SHIFTTWOEXCHANGE)){
        /*TODO HERE INSERT THE CODE FOR PARSING THE NEIGHBORHOODS*/
        double randomFactor = tm.getDecimal();
        unsigned int urgencyPolicy = tm.getInteger();
        double refuelFactor = tm.getDecimal();
        double deliveredQuantityFactor = tm.getDecimal();
        return new emili::irp::irpShiftTwoExchangeNeighborhood(*istance, randomFactor, urgencyPolicy, refuelFactor, deliveredQuantityFactor);
    }
    else if(tm.checkToken(IRP_SHIFTINSERT)){
        /*TODO HERE INSERT THE CODE FOR PARSING THE NEIGHBORHOODS*/
        double randomFactor = tm.getDecimal();
        unsigned int urgencyPolicy = tm.getInteger();
        double refuelFactor = tm.getDecimal();
        double deliveredQuantityFactor = tm.getDecimal();
        return new emili::irp::irpShiftInsertNeighborhood(*istance, randomFactor, urgencyPolicy, refuelFactor, deliveredQuantityFactor);
    }
    else if(tm.checkToken(IRP_SHIFTREMOVE)){
        /*TODO HERE INSERT THE CODE FOR PARSING THE NEIGHBORHOODS*/
        double randomFactor = tm.getDecimal();
        unsigned int urgencyPolicy = tm.getInteger();
        double refuelFactor = tm.getDecimal();
        double deliveredQuantityFactor = tm.getDecimal();
        return new emili::irp::irpShiftRemoveNeighborhood(*istance, randomFactor, urgencyPolicy, refuelFactor, deliveredQuantityFactor);
    }
    else if(tm.checkToken(IRP_REFUEL)){
        /*TODO HERE INSERT THE CODE FOR PARSING THE NEIGHBORHOODS*/
        double refuelStep = tm.getDecimal();
        double deliveredQuantityStep = tm.getDecimal();
        double randomFactor = tm.getDecimal();
        unsigned int urgencyPolicy = tm.getInteger();
        return new emili::irp::irpRefuelNeighborhood(*istance, refuelStep, deliveredQuantityStep, randomFactor, urgencyPolicy);
    }
    else
    {	
        neigh =  nullptr;
    }
    prs::decrementTabLevel();
    return neigh;
}

void prs::irp::IrpParser::neighs(prs::TokenManager& tm)
{
    std::vector<emili::Neighborhood*> vnds;
    vnds.push_back(neigh(tm));
    nes = vnds;
    neighs1(tm);
}

void prs::irp::IrpParser::neighs1(prs::TokenManager& tm)
{
    emili::Neighborhood* n = neighV(tm);
	if(n!=nullptr)
	{
           nes.push_back(n);
           neighs1(tm);
	}

}


void prs::irp::IrpParser::problem(prs::TokenManager& tm)
{

    char *problem_type = tm.nextToken();
    /*TODO HERE THE CODE TO INSTANTIATE THE PROBLEM*/
    /*the problem should be a global variable accessible to all the function of the parser
      so that all the components that need it can access it.
    */

    this->istance = new emili::irp::InventoryRoutingProblem(tm.tokenAt(1));
}


emili::LocalSearch* prs::irp::IrpParser::buildAlgo(prs::TokenManager& tm)
{
    problem(tm);  
    emili::LocalSearch* local = eparams(tm);
    std::cout << "------" << std::endl;
    return local;
}

bool prs::irp::IrpParser::isParsable(string &problem)
{
    /*your parser can support more than one problem definition, extend this if/else to insert
      them*/
    if(strcmp(problem.c_str(),PROBLEM_IRP/*TODO INSERT YOUR PROBLEM DEFINITION HERE*/)==0)
    {

        return true;
    }
    else
    {
        return false;
    }
}

std::string prs::irp::IrpParser::availableProblems() const
{
    ostringstream oss;
    oss << "IRP";

    return oss.str();
}


