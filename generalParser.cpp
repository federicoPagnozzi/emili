//
//  Created by Federico Pagnozzi on 28/11/14.
//  Copyright (c) 2014 Federico Pagnozzi. All rights reserved.
//  This file is distributed under the BSD 2-Clause License. See LICENSE.TXT
//  for details.
#include "generalParser.h"
#include <cstdlib>
#include <cstdio>
#include <ctime>
#include <cstring>
#include <iostream>
#include <sstream>
#include <algorithm>

/* modifiers */
#define RO "-ro"
#define IT "-it"
#define RNDSEED "rnds"
#define PRINT_SOLUTION "ps"
#define DEFAULT_TS 10
#define DEFAULT_TI 10
#define DEFAULT_IT 0
#define GIT_COMMIT_NUMBER "5d88679782f6b293057e3c98f04b056de6c05b35"
/*Base Algos */
#define IG "ig"
#define ILS "ils"
#define FEASIBLE_ILS "fils"
#define TABU "tabu"
#define FIRST "first"
#define FEASIBLE_FIRST "ffirst"
#define TB_FIRST "tfirst"
#define BEST "best"
#define FEASIBLE_BEST "fbest"
#define TB_BEST "tbest"
#define VND "vnd"
#define GVNS_ILS "gvns"
#define TEST_INIT "stin"
#define EMPTY_LOCAL "nols"
/*Base Termination criteria*/
#define TERMINATION_MAXSTEPS "maxstep"
#define TERMINATION_MAXSTEPS_OR_LOCMIN "msorlocmin"
#define TERMINATION_MAXSTEPS_OR_LOCMIN "msorlocmin"
#define TERMINATION_TIME "time"
#define TERMINATION_TIMERO "timero"
#define TERMINATION_LOCMIN "locmin"
#define TERMINATION_ITERA "iteration"
#define TERMINATION_WTRUE "true"
#define TERMINATION_SOA "soater"
/*base permutation flowshop solution perturbations */
#define PERTURBATION_RANDOM_MOVE "rndmv"
#define PERTURBATION_VNRANDOM_MOVE "vnrmv"
#define PERTURBATION_NOPER "noper"
/*base acceptance criteria*/
#define ACCEPTANCE_PROB "prob"
#define ACCEPTANCE_METRO "metropolis"
#define ACCEPTANCE_PMETRO "pmetro"
#define ACCEPTANCE_IMPROVE_PLATEAU "implat"
#define ACCEPTANCE_ALWAYS "always"
#define ACCEPTANCE_INTENSIFY "intensify"
#define ACCEPTANCE_DIVERSIFY "diversify"
#define ACCEPTANCE_IMPROVE "improve"
#define ACCEPTANCE_SA_METRO "sa_metropolis"
#define ACCEPTANCE_SA "saacc"
/*base neighborhoods*/
#define NEIGHBORHOOD_RANDCONHE "rch"

int tab_level = 0;

void prs::printTab(const char* string)
{
    for(int i=0;i<tab_level; i++)
    {
        std::cout << "  ";
    }

    std::cout << string << std::endl;
}

void prs::incrementTabLevel()
{
    tab_level++;
}

void prs::decrementTabLevel()
{
    tab_level--;
    tab_level= tab_level<0?0:tab_level;
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
    std::cout << "commit : " << GIT_COMMIT_NUMBER << std::endl;
}

void prs::info()
{
    std::cout << "Usage:" << std::endl;
    std::cout << std::endl;
    std::cout << "EMILI INSTANCE_FILE_PATH PROBLEM <PROBLEM OPTIONS> [-it|-ro time] [rnds seed] [ps]" << std::endl;
    std::cout << std::endl;
}


void prs::check(char* t,const char* message)
{
    if(t==nullptr)
    {
        prs::info();
        std::cerr <<"PARSING ERROR "<< message << std::endl;
        exit(-1);
    }
}

char* prs::TokenManager::nextToken()
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

char* prs::TokenManager::peek()
{
    if(currentToken < numberOfTokens)
    {
        char* token = tokens[currentToken];
     //   currentToken++;
        return token;
    }
    else
    {
        return " ";
    }
}

char* prs::TokenManager::operator *()
{
    return peek();
}

void prs::TokenManager::next()
{
    if(currentToken<numberOfTokens)
    {
        currentToken++;
    }
}

void prs::TokenManager::operator ++(int k)
{
    next();
}

char* prs::TokenManager::tokenAt(int i)
{
    if(i < numberOfTokens)
    {
        return tokens[i];
    }
    return nullptr;
}

bool prs::TokenManager::checkToken(std::string &token)
{
    if(currentToken < numberOfTokens)
    {
        if(strcmp(tokens[currentToken],token.c_str())==0)
        {
            currentToken++;
            return true;
        }
    }
        return false;
}

bool prs::TokenManager::checkToken(const char* token)
{
    if(currentToken < numberOfTokens)
    {
        if(strcmp(tokens[currentToken],token)==0)
        {
            currentToken++;
            return true;
        }
    }
        return false;
}

bool prs::TokenManager::hasMoreTokens()
{
    return currentToken < numberOfTokens;
}

bool prs::TokenManager::move(int i)
{
    if(i>=0 && i<numberOfTokens)
    {
            previousCurrentToken = currentToken;
            currentToken = i;
            return true;
    }
    return false;
}

void prs::TokenManager::restore()
{
    currentToken = previousCurrentToken;
}

int prs::TokenManager::seek(const char * token)
{
    for(int i=0;i<numberOfTokens;i++)
    {
        if(strcmp(tokens[i],token)==0)
        {
            return i;
        }
    }
    return -1;
}


int prs::TokenManager::getInteger()
{
    char* t = peek();
    check(t,"A NUMBER WAS EXPECTED!");
    std::string num(t);
    if( !std::all_of(num.begin(),num.end(),::isdigit))
       check(nullptr,"A NUMBER WAS EXPECTED!");

    int k = atoi(t);
   // std::cout << k << "\n\t";
    next();
    return k;
}

float prs::TokenManager::getDecimal()
{
    char* t = peek();
    check(t,"A DECIMAL NUMBER WAS EXPECTED!");
    std::string num(t);

    if(std::none_of(num.begin(),num.end(),::isdigit))
        check(nullptr,"A DECIMAL NUMBER WAS EXPECTED!");

    float k = atof(t);
    //std::cout << k << "\n\t";
    next();
    return k;
}

bool prs::AlgoBuilder::operator ==(const AlgoBuilder& b)
{
    if(this->availableProblems()==b.availableProblems())
    {
        return true;
    }
    return false;
}


float getTime(prs::TokenManager& tm,int problemSize)
{
        if(tm.checkToken(IT))
        {
            int n = tm.getDecimal();
            std::ostringstream oss;
            oss << "Run time secs : " << n;
            //printTab(oss.str().c_str());
            std::cout << oss.str() << std::endl;
            return n;
        }
        else if(tm.checkToken(RO))
        {
            float d = tm.getDecimal();
            float time = d*problemSize;
            //int n = floorf(time);
            std::ostringstream oss;
            oss << "Rho = "<< d << " Run time secs : " << time;
            //printTab(oss.str().c_str());
            std::cout << oss.str() << std::endl;
            return time;
        }


    return DEFAULT_IT;
}

int getSeed(prs::TokenManager& tm)
{
    int rnds = 0;
    if(tm.move(tm.seek(RNDSEED)))
    {
        tm++;
        rnds = tm.getInteger();
        tm.restore();
    }
    std::cout << "Random seed : " << rnds << std::endl;
    return rnds;
}


emili::LocalSearch* prs::GeneralParser::parseParams()
{
    tm++;
    tm++;
    char* p = tm.peek();
    if(p != nullptr)
    {
    std::string prob(p);
    for(std::vector< AlgoBuilder*> ::iterator iter= builders.begin(); iter!=builders.end(); ++iter)
    {
        AlgoBuilder* bld = *iter;
        if(bld->isParsable(prob))
        {
            emili::initializeRandom(getSeed(tm));
            emili::LocalSearch* ls = bld->buildAlgo(tm);                        
            ls->setSearchTime(getTime(tm,ls->getInitialSolution().getProblem().problemSize()));            
            if(tm.checkToken(PRINT_SOLUTION))
            {
                emili::set_print(true);
                if(ls->getSearchTime() > 0){
                    std::cout << "The solution will be print at the end of the execution" << std::endl;
                }
            }
            else
            {
                emili::set_print(false);
            }
            return ls;
        }
    }
    std::cout << "------" << std::endl;
    std::cout << "No parser for "<< p << " available" << std::endl;
   }
    std::cerr << "A problem was expected!" << std::endl;
    // info
    prs::info();
 return nullptr;
}

void prs::GeneralParser::registerBuilder(AlgoBuilder* builder)
{
    builders.push_back(builder);
}

void prs::GeneralParser::removeBuilder(AlgoBuilder* builder)
{
    int ind = 0;
    for(;ind<builders.size();ind++)
    {
        if(builder == builders[ind])
        {
            builders.erase(builders.begin()+ind);
            return;
        }
    }
}
/*NEW PARSER STUFF*/
prs::Component& prs::Component::operator=(const Component& a)
{
    this->type = a.type;
    this->rawComponent = a.rawComponent;
    this->token = a.token;
}


prs::Component prs::Builder::retrieveComponent(int type,bool allowEmpty)
{

        Component c = gp.buildComponent(type);
        if( c.is(COMPONENT_NULL) || !c.is(type))
        {
            if(!(c.is(COMPONENT_EMPTY) && allowEmpty))
            gp.fatalError(c.getType(),type);
        }
        return c;
}

prs::Component prs::Builder::buildComponent(int type)
{
    Component c;
    void* raw_pointer=nullptr;
    switch(type)
    {    
    case COMPONENT_ALGORITHM: raw_pointer = buildAlgo();break;
    case COMPONENT_INITIAL_SOLUTION_GENERATOR: raw_pointer = buildInitialSolution();break;
    case COMPONENT_TERMINATION_CRITERION: raw_pointer =  buildTermination();break;    
    case COMPONENT_NEIGHBORHOOD: raw_pointer = buildNeighborhood();c.setType(COMPONENT_EMPTY);break;
    case COMPONENT_PERTURBATION: raw_pointer =  buildPerturbation();break;
    case COMPONENT_ACCEPTANCE: raw_pointer =  buildAcceptance();break;
    case COMPONENT_TABU_TENURE: raw_pointer =  buildTabuTenure();break;
    case COMPONENT_PROBLEM: raw_pointer = buildProblem();break;
    }

    if(raw_pointer != nullptr)
    {
        c.setType(type);
        c.setRawComponent(raw_pointer);
    }
    return c;
}

std::vector<emili::Neighborhood*> prs::Builder::buildNeighborhoodVector()
{
    std::vector<emili::Neighborhood*> nes;
    prs::Component c = retrieveComponent(COMPONENT_NEIGHBORHOOD);
    while(!c.is(COMPONENT_EMPTY))
    {
        nes.push_back(c.get<emili::Neighborhood>());
        c = retrieveComponent(COMPONENT_NEIGHBORHOOD,true);
    }
    return nes;
}


void prs::GeneralParserE::addBuilder(Builder *b)
{
    allbuilders.push_back(b);
}

prs::Component prs::GeneralParserE::buildComponent(int type)
{
    Component last;
    for(std::vector<Builder*>::iterator iter = activeBuilders.begin(); iter != activeBuilders.end();iter++)
    {
        Builder* b = *iter;
        Component c = b->buildComponent(type);
        if(!c.is(COMPONENT_NULL))
        {
            if(c.is(COMPONENT_EMPTY))
            {
                last.setType(COMPONENT_EMPTY);
            }
            else
            {
                return c;
            }
        }
    }
    return last;
}

emili::LocalSearch* prs::GeneralParserE::parseParams()
{
    tm++;
    tm++;
    char* p = tm.peek();
    if(p != nullptr)
    {
    std::string prob(p);
    Builder* probBuilder = nullptr;
    for(std::vector< Builder*> ::iterator iter= allbuilders.begin(); iter!=allbuilders.end(); ++iter)
    {
        Builder* bld = *iter;

        if(bld->isCompatibleWith(prob))
        {
           activeBuilders.push_back(bld);
           if(bld->canOpenInstance(prob))
           {
               probBuilder = bld;
           }
        }
    }
    if(probBuilder != nullptr)
    {
        emili::initializeRandom(getSeed(tm));
        this->instance = probBuilder->openInstance();
        emili::LocalSearch* ls = buildComponent(COMPONENT_ALGORITHM).get<emili::LocalSearch>();
        ls->setSearchTime(getTime(tm,ls->getInitialSolution().getProblem().problemSize()));        
        if(tm.seek(PRINT_SOLUTION)>0)
        {
            emili::set_print(true);
            if(ls->getSearchTime() > 0){
                std::cout << "The solution will be print at the end of the execution" << std::endl;
            }
        }
        else
        {
            std::cout << "The print option is not enabled" << std::endl;
            emili::set_print(false);
        }
        return ls;
    }
    else
    {
        std::cout << "------" << std::endl;
        std::cout << "No parser for "<< p << " available" << std::endl;
    }
   }
    std::cerr << "A problem was expected!" << std::endl;
    // info
    prs::info();
 return nullptr;
}

void prs::GeneralParserE::fatalError(int received_type,int expected_type)
{
    std::cerr << "FATAL ERROR!\n";
    std::cerr << "I was expecting a " << typeName(expected_type) << " for token '" << tm.peek() << "' but I received " << typeName(received_type)<< std::endl;
    exit(-1);
}

std::string prs::GeneralParserE::typeName(int type)
{
    switch(type)
    {
    case COMPONENT_EMPTY: return std::string("EPSILON");
    case COMPONENT_NULL: return std::string("NULL");
    case COMPONENT_ALGORITHM: return std::string("Algorithm");
    case COMPONENT_INITIAL_SOLUTION_GENERATOR: return std::string("Initial solution");
    case COMPONENT_TERMINATION_CRITERION: return std::string("Termination");
    case COMPONENT_NEIGHBORHOOD: return std::string("Neighborhood");
    case COMPONENT_TABU_TENURE: return std::string("Tabu tenure");
    case COMPONENT_NEIGHBORHOOD_OR_EMPTY: return std::string("Neighborhood or empty");
    case COMPONENT_PERTURBATION: return std::string("Perturbation");
    case COMPONENT_ACCEPTANCE: return std::string("Acceptance");
    }
    /*TODO - if Gp does not have a name for the component it should ask the builders*/
    for(std::vector<Builder*>::iterator iter = activeBuilders.begin(); iter != activeBuilders.end();iter++)
    {
        Builder* b = *iter;
        std::string s = b->typeName(type);
        if(!s.empty())
        {
            return s;
        }
    }
    std::ostringstream oss;
    oss << type;
    return oss.str();
}



emili::LocalSearch* prs::EmBaseBuilder::buildAlgo()
{
    prs::incrementTabLevel();
    emili::LocalSearch* ls = nullptr;
    if(tm.checkToken(ILS))
    {
        printTab("ILS ");
        emili::LocalSearch* lls = retrieveComponent(COMPONENT_ALGORITHM).get<emili::LocalSearch>();
        emili::Termination* pft = retrieveComponent(COMPONENT_TERMINATION_CRITERION).get<emili::Termination>();
        emili::Perturbation* prsp = retrieveComponent(COMPONENT_PERTURBATION).get<emili::Perturbation>();        
        emili::Acceptance* tac = retrieveComponent(COMPONENT_ACCEPTANCE).get<emili::Acceptance>();
        ls = new emili::IteratedLocalSearch(*lls,*pft,*prsp,*tac);
    }else if(tm.checkToken(FEASIBLE_ILS))
    {
        printTab("ILS that returns a feasible solution if it generates one");
        emili::LocalSearch* lls = retrieveComponent(COMPONENT_ALGORITHM).get<emili::LocalSearch>();
        emili::Termination* pft = retrieveComponent(COMPONENT_TERMINATION_CRITERION).get<emili::Termination>();
        emili::Perturbation* prsp = retrieveComponent(COMPONENT_PERTURBATION).get<emili::Perturbation>();
        emili::Acceptance* tac = retrieveComponent(COMPONENT_ACCEPTANCE).get<emili::Acceptance>();
        ls = new emili::FeasibleIteratedLocalSearch(*lls,*pft,*prsp,*tac);
    }else if(tm.checkToken(TABU))
    {
        printTab("TABU SEARCH");
        bool best = -1;
        if(tm.checkToken(BEST))
        {
            best = 1;
        }else if(tm.checkToken(FIRST))
        {
           best = 0;
        }
        if(best >= 0)
        {
            emili::InitialSolution* in = retrieveComponent(COMPONENT_INITIAL_SOLUTION_GENERATOR).get<emili::InitialSolution>();
            emili::Termination* te = retrieveComponent(COMPONENT_TERMINATION_CRITERION).get<emili::Termination>();
            emili::Neighborhood* ne = retrieveComponent(COMPONENT_NEIGHBORHOOD).get<emili::Neighborhood>();
            emili::TabuMemory* tmem = retrieveComponent(COMPONENT_TABU_TENURE).get<emili::TabuMemory>();
            if(best)
            {
                ls =  new emili::BestTabuSearch(*in,*te,*ne,*tmem);
            }
            else
            {
                ls = new emili::FirstTabuSearch(*in,*te,*ne,*tmem);
            }
        }
    }
    else if(tm.checkToken(FIRST))
    {
        printTab("FIRST IMPROVEMENT");
        emili::InitialSolution* in = retrieveComponent(COMPONENT_INITIAL_SOLUTION_GENERATOR).get<emili::InitialSolution>();
        emili::Termination* te = retrieveComponent(COMPONENT_TERMINATION_CRITERION).get<emili::Termination>();
        emili::Neighborhood* ne = retrieveComponent(COMPONENT_NEIGHBORHOOD).get<emili::Neighborhood>();
        ls =  new emili::FirstImprovementSearch(*in,*te,*ne);
    }
    else if(tm.checkToken(BEST))
    {
        printTab("BEST IMPROVEMENT");
        emili::InitialSolution* in = retrieveComponent(COMPONENT_INITIAL_SOLUTION_GENERATOR).get<emili::InitialSolution>();
        emili::Termination* te = retrieveComponent(COMPONENT_TERMINATION_CRITERION).get<emili::Termination>();
        emili::Neighborhood* ne = retrieveComponent(COMPONENT_NEIGHBORHOOD).get<emili::Neighborhood>();
        ls =  new emili::BestImprovementSearch(*in,*te,*ne);
    }
    else if(tm.checkToken(FEASIBLE_FIRST))
    {
        printTab("FIRST IMPROVEMENT that returns a feasible solution if it generates one");
        emili::InitialSolution* in = retrieveComponent(COMPONENT_INITIAL_SOLUTION_GENERATOR).get<emili::InitialSolution>();
        emili::Termination* te = retrieveComponent(COMPONENT_TERMINATION_CRITERION).get<emili::Termination>();
        emili::Neighborhood* ne = retrieveComponent(COMPONENT_NEIGHBORHOOD).get<emili::Neighborhood>();
        ls =  new emili::FeasibleFirstImprovementSearch(*in,*te,*ne);
    }
    else if(tm.checkToken(FEASIBLE_BEST))
    {
        printTab("BEST IMPROVEMENT that returns a feasible solution if it generates one");
        emili::InitialSolution* in = retrieveComponent(COMPONENT_INITIAL_SOLUTION_GENERATOR).get<emili::InitialSolution>();
        emili::Termination* te = retrieveComponent(COMPONENT_TERMINATION_CRITERION).get<emili::Termination>();
        emili::Neighborhood* ne = retrieveComponent(COMPONENT_NEIGHBORHOOD).get<emili::Neighborhood>();
        ls =  new emili::FeasibleBestImprovementSearch(*in,*te,*ne);
    }
    else if(tm.checkToken(TB_FIRST))
    {
        printTab("BEST IMPROVEMENT");
        emili::InitialSolution* in = retrieveComponent(COMPONENT_INITIAL_SOLUTION_GENERATOR).get<emili::InitialSolution>();
        emili::Termination* te = retrieveComponent(COMPONENT_TERMINATION_CRITERION).get<emili::Termination>();
        emili::Neighborhood* ne = retrieveComponent(COMPONENT_NEIGHBORHOOD).get<emili::Neighborhood>();
        emili::Problem* p = retrieveComponent(COMPONENT_PROBLEM).get<emili::Problem>();
        ls =  new emili::TieBrakingFirstImprovementSearch(*in,*te,*ne,*p);
    }
    else if(tm.checkToken(TB_BEST))
    {
        printTab("BEST IMPROVEMENT");
        emili::InitialSolution* in = retrieveComponent(COMPONENT_INITIAL_SOLUTION_GENERATOR).get<emili::InitialSolution>();
        emili::Termination* te = retrieveComponent(COMPONENT_TERMINATION_CRITERION).get<emili::Termination>();
        emili::Neighborhood* ne = retrieveComponent(COMPONENT_NEIGHBORHOOD).get<emili::Neighborhood>();
        emili::Problem* p = retrieveComponent(COMPONENT_PROBLEM).get<emili::Problem>();
        ls =  new emili::TieBrakingBestImprovementSearch(*in,*te,*ne,*p);
    }
    else if(tm.checkToken(VND))
    {
        //printTab("VND SEARCH");
        bool best = -1;
        if(tm.checkToken(BEST))
        {
            printTab("BEST IMPROVEMENT VND");
            best = 1;
        }else if(tm.checkToken(FIRST))            
        {
           printTab("FIRST IMPROVEMENT VND");
           best = 0;
        }
        if(best >= 0)
        {
            emili::InitialSolution* in = retrieveComponent(COMPONENT_INITIAL_SOLUTION_GENERATOR).get<emili::InitialSolution>();
            emili::Termination* te = retrieveComponent(COMPONENT_TERMINATION_CRITERION).get<emili::Termination>();
            std::vector<emili::Neighborhood*> nes = buildNeighborhoodVector();
            if(best)
            {
                ls =  new emili::VNDSearch<emili::BestImprovementSearch>(*in,*te,nes);
            }
            else
            {
                 ls =  new emili::VNDSearch<emili::FirstImprovementSearch>(*in,*te,nes);
            }
        }
    }
    else if(tm.checkToken(TEST_INIT))
    {
        emili::InitialSolution* ini = retrieveComponent(COMPONENT_INITIAL_SOLUTION_GENERATOR).get<emili::InitialSolution>();
        /*ls = new emili::EmptyLocalSearch(*ini);*/
        clock_t time = clock();
        emili::Solution* s = ini->generateSolution();
        double time_elapsed = (double)(clock()-time)/CLOCKS_PER_SEC;
        std::cout << "time : " << time_elapsed << std::endl;
        std::cout << s->getSolutionRepresentation() << std::endl;
        std::cout << s-> getSolutionValue() << std::endl;
        std::cerr << s-> getSolutionValue() << std::endl;
        exit(123);
    }
    else if(tm.checkToken(EMPTY_LOCAL))
    {
        printTab("NO LOCAL SEARCH");
        emili::InitialSolution* ini = retrieveComponent(COMPONENT_INITIAL_SOLUTION_GENERATOR).get<emili::InitialSolution>();
        ls = new emili::EmptyLocalSearch(*ini);
    }

    prs::decrementTabLevel();
    return ls;
}

emili::Termination* prs::EmBaseBuilder::buildTermination()
{
    prs::incrementTabLevel();
    emili::Termination* term=nullptr;
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
    else if(tm.checkToken(TERMINATION_TIMERO))
    {

        float time =tm.getDecimal();
        if(time==0){
            time = 1;
        }
        std::ostringstream oss;
        int rtime = time * gp.getInstance()->problemSize();
        oss << "Timed termination";
        printTab(oss.str().c_str());
        oss.str("");
        prs::incrementTabLevel();
        oss << "ro : " << time;
        printTab(oss.str().c_str());
        oss.str("");
        oss << "running time : " << rtime;
        printTab(oss.str().c_str());
        prs::decrementTabLevel();
        term =  new emili::TimedTermination(rtime);
    }
    else if(tm.checkToken(TERMINATION_MAXSTEPS))
    {
        int steps = tm.getInteger();
        std::ostringstream oss;
        oss << "Max Steps termination. # steps: "<< steps;
        printTab(oss.str().c_str());
        term = new emili::MaxStepsTermination(steps);
    }
    else if(tm.checkToken(TERMINATION_MAXSTEPS_OR_LOCMIN))
    {
        int steps = tm.getInteger();
        std::ostringstream oss;
        oss << "Max Steps termination or when reaching locmin. # steps: "<< steps;
        printTab(oss.str().c_str());
        term = new emili::MaxStepsOrLocmin(steps);
    }


    prs::decrementTabLevel();
    return term;
}

emili::Perturbation* prs::EmBaseBuilder::buildPerturbation()
{
    prs::incrementTabLevel();    
    std::ostringstream oss;
    emili::Perturbation* per=nullptr;
    if(tm.checkToken(PERTURBATION_RANDOM_MOVE))
    {
        printTab("Random move perturbation.");
        emili::Neighborhood* n = retrieveComponent(COMPONENT_NEIGHBORHOOD).get<emili::Neighborhood>();
        int num = tm.getInteger();
        prs::incrementTabLevel();
        oss.str(""); oss  << "number of moves per perturbation step " << num;
        printTab(oss.str().c_str());
        prs::decrementTabLevel();
        per = new emili::RandomMovePerturbation(*n,num);
    }
    else if(tm.checkToken(PERTURBATION_NOPER))
    {
        printTab("No PERTURBATION.");
        per = new emili::NoPerturbation();
    }
    else if(tm.checkToken(PERTURBATION_VNRANDOM_MOVE))
    {
        printTab("Random move perturbation." );
        prs::incrementTabLevel();
        int num = tm.getInteger();
        oss.str(""); oss  << "number of moves per perturbation step " << num << ".\n\t";
        printTab(oss.str().c_str());
        int iter = tm.getInteger();
        oss.str(""); oss  << "number of iteration before changing the neighborhood " << iter << ".\n\t";
        printTab(oss.str().c_str());
        prs::decrementTabLevel();
        std::vector<emili::Neighborhood*> nes = this->buildNeighborhoodVector();
        per = new emili::VNRandomMovePerturbation(nes,num,iter);
    }

    prs::decrementTabLevel();
    return per;
}

emili::Neighborhood* prs::EmBaseBuilder::buildNeighborhood()
{
    prs::incrementTabLevel();
    emili::Neighborhood* neigh = nullptr;
    if(tm.checkToken(NEIGHBORHOOD_RANDCONHE))
    {
        printTab("Random Constructive Heuristic Neighborhood ");
        emili::InitialSolution* heuristic = retrieveComponent(COMPONENT_INITIAL_SOLUTION_GENERATOR).get<emili::InitialSolution>();
        neigh = new emili::RandomConstructiveHeuristicNeighborhood(*heuristic);
    }
    prs::decrementTabLevel();
    return neigh;
}

emili::Acceptance* prs::EmBaseBuilder::buildAcceptance()
{
    prs::incrementTabLevel();   
    emili::Acceptance* acc = nullptr;
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
            //std::cout << info() << std::endl;
            exit(-123);
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
        oss.str(""); oss  << "metropolis acceptance.";
        printTab(oss.str().c_str());
        prs::incrementTabLevel();
        oss.str(""); oss << "start: "<< start;
        printTab(oss.str().c_str());
        oss.str(""); oss << "end: "<< end;
        printTab(oss.str().c_str());
        oss.str(""); oss << "ratio: "<< ratio;
        printTab(oss.str().c_str());
        prs::decrementTabLevel();
        acc = new  emili::Metropolis(start,end,ratio);
    }
    else  if(tm.checkToken(ACCEPTANCE_PMETRO))
    {
        float start =tm.getDecimal();
        float end =tm.getDecimal();
        float ratio =tm.getDecimal();
        int iterations = tm.getInteger();
        oss.str(""); oss  << "metropolis acceptance.";
        printTab(oss.str().c_str());
        prs::incrementTabLevel();
        oss.str(""); oss << "start: "<< start;
        printTab(oss.str().c_str());
        oss.str(""); oss << "end: "<< end;
        printTab(oss.str().c_str());
        oss.str(""); oss << "ratio: "<< ratio;
        printTab(oss.str().c_str());
        oss.str(""); oss << "frequence: "<< iterations;
        printTab(oss.str().c_str());
        prs::decrementTabLevel();
        acc = new  emili::Metropolis(start,end,ratio,iterations);
    }
    else if(tm.checkToken(ACCEPTANCE_SA))
    {
        float start =tm.getDecimal();
        float end =tm.getDecimal();
        float ratio =tm.getDecimal();
        int iterations = tm.getInteger();
        float alpha =tm.getDecimal();
        oss.str(""); oss  << "metropolis acceptance.";
        printTab(oss.str().c_str());
        prs::incrementTabLevel();
        oss.str(""); oss << "start: "<< start;
        printTab(oss.str().c_str());
        oss.str(""); oss << "end: "<< end;
        printTab(oss.str().c_str());
        oss.str(""); oss << "ratio: "<< ratio;
        printTab(oss.str().c_str());
        oss.str(""); oss << "frequence: "<< iterations;
        printTab(oss.str().c_str());
        oss.str(""); oss << "alpha: "<< alpha;
        printTab(oss.str().c_str());
        prs::decrementTabLevel();
        acc = new  emili::Metropolis(start,end,ratio,iterations,alpha);
    }
    else if(tm.checkToken(ACCEPTANCE_IMPROVE_PLATEAU))
    {
        int plateau_steps = tm.getInteger();
        int threshold = tm.getInteger();
        oss.str(""); oss  << "Accept a diversification solution if it improves on the intensification otherwise it will accept "<< plateau_steps << " non improving steps once it reaches the threshold of " << threshold;
        printTab(oss.str().c_str());
        acc = new  emili::AcceptPlateau(plateau_steps,threshold);
    }


    prs::decrementTabLevel();
    return acc;
}


