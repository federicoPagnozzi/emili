//
//  Created by Federico Pagnozzi on 28/11/14.
//  Copyright (c) 2014 Federico Pagnozzi. All rights reserved.
//  This file is distributed under the BSD 2-Clause License. See LICENSE.TXT
//  for details.

#include "emilibase.h"
#include <cstdlib>
#include <cstdio>
#include <signal.h>
#include <ctime>
#include <sstream>

#if defined(_WIN32) || defined(_WIN64)
//no signals compilation path for windows.
#define NOSIG 1
#else

#include <sys/time.h>

#endif

#include <iostream>
#include <assert.h>
/**
 * WARNING!!!
 * Adding data structures to a solution subclass could broken this method

emili::Solution& emili::Solution::operator=(const emili::Solution& a)
{
    this->instance = a.instance;
    this->setRawData(a.getRawData());
    return *this;
}
 */
/**
 * RANDOM NUMBER GENERATOR
 */
#ifndef NOC11

std::mt19937 generator;
std::uniform_int_distribution<int> distribution;
std::uniform_real_distribution<float> realdistr;
void emili::initializeRandom(int seed)
{
    generator = std::mt19937(seed);
    //rand = std::bind(distribution,generator);
}

std::mt19937& emili::getRandomGenerator()
{
    return generator;
}

#else
//Random generation compilation path for compilers that don't support c++11
std::tr1::mt19937 generator;
std::tr1::uniform_int<int> distribution;
std::tr1::uniform_real<float> realdistr;
void emili::initializeRandom(int seed)
{
    generator = std::tr1::mt19937(seed);
    //rand = std::bind(distribution,generator);
}

std::tr1::mt19937& emili::getRandomGenerator()
{
    return generator;
}

#endif


int emili::generateRandomNumber()
{
   // auto rand = std::bind(distribution,generator);
   return distribution(generator);
}

float emili::generateRealRandomNumber()
{
    return realdistr(generator);
}


/**
 * TIMED SEARCH CODE
 */
bool print;
bool keep_going;
bool timer_keep_going;
clock_t endTime;
clock_t beginTime;
clock_t s_time;
long final_iteration_counter = 0;
emili::LocalSearch* localsearch = nullptr;
emili::Solution* s_cap = nullptr;
double sol_val = -1;

void emili::initializeTimerBaseSolution(emili::Solution* base_solution)
{
    s_cap = base_solution;
}

bool emili::getKeepGoing()
{
    return keep_going;
}

double emili::getCurrentExecutionTime()
{
    return (double)((clock()-beginTime)/ (double)CLOCKS_PER_SEC);
}

emili::LocalSearch* emili::getAlgo()
{
    return localsearch;
}

void emili::setRootAlgorithm(emili::LocalSearch *ls)
{
    if(ls != nullptr)
    {
        localsearch = ls;
    }
}

static void finalise (int _)
{
    keep_going = false;
    endTime = clock();
    *s_cap = *localsearch->getBestSoFar();
    if(s_cap != nullptr)
    {
        sol_val = s_cap->getSolutionValue();
        if(print)
        {
            //messages << "CPU time: " << (endTime - beginTime) / (float)CLOCKS_PER_SEC << std::endl;
            //messages << "iteration counter : " << emili::iteration_counter()<< std::endl;
            //messages << "objective function value : "<< std::fixed << sol_val << std::endl;
            final_iteration_counter = emili::iteration_counter();
            //messages << "solution : " << s_cap->getSolutionRepresentation() << std::endl;
            //std::cout << "Reached at time: " << (s_time - beginTime) / (float)CLOCKS_PER_SEC << std::endl;
            //std::cerr << (endTime - beginTime) / (float)CLOCKS_PER_SEC << " ";
        }
        else
        {
            /*std::cout << "CPU time: " << (endTime - beginTime) / (float)CLOCKS_PER_SEC << std::endl;
            std::cout << "iteration counter : " << emili::iteration_counter()<< std::endl;
            std::cerr << std::fixed << sol_val << std::endl;
            std::cerr << std::flush;*/
        }
    }
    else
    {
        if(!print)
        {         
            std::cout  << "No valid solution found!" << std::endl;
        }
    }
    //std::cout << std::flush;
    /*if(!print)
    {
       exit(0);
    }*/
}


#ifndef NOSIG
struct itimerval timer;
struct itimerval termination_timer;

void timeUp(int _)
{
    timer_keep_going = false;
}

void setTerminationTimer(int time)
{

    termination_timer.it_value.tv_sec = time;
    termination_timer.it_value.tv_usec = 0;
    termination_timer.it_interval.tv_sec = 0;
    termination_timer.it_interval.tv_usec = 0;
    signal(SIGPROF, timeUp);
    signal(SIGINT, timeUp);
    if (setitimer (ITIMER_PROF, &termination_timer, NULL) != 0) {
        printf("error in setitimer\n");
        exit(10);
    }else{
       timer_keep_going = true;
    }
}

static inline bool isTimerUp()
{

      itimerval current_timer;
       getitimer(ITIMER_PROF, &current_timer);
      return (current_timer.it_value.tv_sec != 0 ||
              current_timer.it_value.tv_usec != 0);

}

void lastPrint()
{
    if(s_cap != nullptr)
    {
        double sol_val2 = s_cap->getSolutionValue();
        std::cout << "CPU time: " << (endTime - beginTime) / (float)CLOCKS_PER_SEC << std::endl;
        std::cout << "iteration counter : " << emili::iteration_counter()<< std::endl;
        std::cout << "objective function value : "<< std::fixed << sol_val << std::endl;
        std::cout << "solution : " << s_cap->getSolutionRepresentation() << std::endl;
        if(sol_val2 != sol_val)
        {
            std::cout << "value of printed solution : " << sol_val2 << std::endl;
        }
    }
    else {
        std::cout  << "No valid solution found!" << std::endl;
    }
}

int max_time = -1 ;
static inline void setTimer(float maxTime)
{
    keep_going = true;
    int secs = floorf(maxTime);
    float usecs = (maxTime-secs)*1000*1000;
    timer.it_value.tv_sec = secs;
    timer.it_value.tv_usec = usecs;
    timer.it_interval.tv_sec = 0;
    timer.it_interval.tv_usec = 0;
    emili::iteration_counter_zero();
    signal(SIGPROF, finalise);
    signal(SIGINT, finalise);
    if(print)
        atexit(lastPrint);
    if (setitimer (ITIMER_PROF, &timer, NULL) != 0) {
        printf("error in setitimer\n");
        exit(10);
    }else{
        std::cout << "timer set " << maxTime << " seconds " << std::endl;       
        max_time = maxTime;
    }
}

static inline void stopTimer()
{
    std::cout << "timer stopped" << std::endl;

    struct itimerval zero_timer = { 0 };
   setitimer(ITIMER_PROF, &zero_timer, &timer);

}
#else


static inline void stopTimer()
{
    std::cout << "timer stopped" << std::endl;
}

static inline void setTimer(float maxTime)
{
    //std::cout << "timer set"
}
#endif

/**
 * Timer HOOK
 */

/**
 * Iteration counter
 */
static unsigned long iteration_counter_ ;

void emili::iteration_counter_zero()
{
    iteration_counter_ = 0;
}


unsigned long emili::iteration_counter(){
    return iteration_counter_;
}

void emili::iteration_increment(){
    iteration_counter_++;
}

void emili::iteration_decrement(){
    iteration_counter_--;
}

/**
 * Print Solution info
 */
double cbest = std::numeric_limits<double>::max();
inline void emili::printSolstats(emili::Solution* sol)
{
#ifdef WITH_STATS
    double nbest = sol->getSolutionValue();
    if(print && cbest > nbest)
    {
      cbest = nbest;
      std::cout << (clock() - beginTime) / (float)CLOCKS_PER_SEC << " , " << sol->getSolutionValue() << " , " << iteration_counter_ << "\n";
    }
#endif
}

inline void emili::printSearchstats(emili::SearchStatus* status)
{
#ifdef WITH_STATS
    emili::Solution* sol = status->getBestSolution();
    if(print && cbest > sol->getSolutionValue())
    {
      cbest = sol->getSolutionValue();
      std::cout << (clock() - beginTime) / (float)CLOCKS_PER_SEC << " , "
                << sol->getSolutionValue() << " , "
                << iteration_counter_      << " , "
                << status->total_counter   << " , "
                << status->not_improved    << "\n";
    }
#endif
}

void emili::set_print(bool p)
{
    print = p;
}

bool emili::get_print()
{
    return print;
}

/**
  * Problem Implementation
  */

double emili::Problem::solutionDistance(Solution& solution1, Solution& solution2)
{
    return evaluateSolution(solution1)-evaluateSolution(solution2);
}

int emili::Problem::problemSize()
{
    return 1;
}

/**
 * Solution implementation
 */
emili::Solution& emili::Solution::operator=(const emili::Solution& a)
{
    this->solution_value = a.solution_value;
    this->setRawData(a.getRawData());
    return *this;
}

bool emili::Solution::operator<(emili::Solution& a)
{

    return solution_value < a.solution_value;
}

bool emili::Solution::operator<=(emili::Solution& a)
{    
    return solution_value <= a.solution_value;
}

bool emili::Solution::operator>=(emili::Solution& a)
{
    return solution_value >= a.solution_value;
}

bool emili::Solution::operator>(emili::Solution& a)
{
    return solution_value > a.solution_value;
}

bool emili::Solution::operator==(emili::Solution& a)
{
    return solution_value == a.solution_value;
}

double emili::Solution::getSolutionValue()
{    
    return solution_value;
}

std::string emili::Solution::getSolutionRepresentation()
{
    std::string s("");
    return s;
}

void emili::Solution::setSolutionValue(double value)
{
    this->solution_value = value;
}

/**
 * Base implementation of Neighborhood class
 */
emili::Neighborhood::NeighborhoodIterator& emili::Neighborhood::NeighborhoodIterator::operator =(const emili::Neighborhood::NeighborhoodIterator& iter)
{
    line_ = iter.line_;
    n = iter.n;
    return *this;
}

bool emili::Neighborhood::NeighborhoodIterator::operator ==(const emili::Neighborhood::NeighborhoodIterator& iter)
{
    return (this->n == iter.n) && (this->line_ == iter.line_);
}

bool emili::Neighborhood::NeighborhoodIterator::operator !=(const emili::Neighborhood::NeighborhoodIterator& iter)
{
    return (this->n != iter.n) || (this->line_ != iter.line_);
}

emili::Neighborhood::NeighborhoodIterator& emili::Neighborhood::NeighborhoodIterator::operator++()
{    
    line_->setSolutionValue(base_value);
    n->reverseLastMove(line_);
    this->line_ = n->computeStep(this->line_);
    return *this;
}

emili::Solution* emili::Neighborhood::NeighborhoodIterator::operator *()
{
    return line_;
}

emili::Neighborhood::NeighborhoodIterator emili::Neighborhood::begin(emili::Solution *startSolution)
{
    return emili::Neighborhood::NeighborhoodIterator(this,startSolution);
}

emili::Neighborhood::NeighborhoodIterator emili::Neighborhood::end()
{
    return emili::Neighborhood::NeighborhoodIterator(this,nullptr);
}

emili::Solution* emili::RandomConstructiveHeuristicNeighborhood::computeStep(Solution *step)
{
     if(state==0)
     {
         state=1;
         emili::Solution* nsol = heuristic->generateSolution();
         *step = *nsol;
         delete nsol;
         return step;
     }
     else
     {
         state=0;
         return nullptr;
     }
}
/**
 * SearchStatus
 */

void emili::SearchStatus::incrementCounters()
{
    counter += 1;
    total_counter += 1;
    not_improved += 1;
}

void emili::SearchStatus::newBestSolution(Solution* new_best)
{
    not_improved = 0;
    //delete best;
    /*    if(best==nullptr)
        {
            best = new_best->clone();
        }
        else
        {
            *best = *new_best;
        }
    */
    best = new_best;
    if(new_best->isFeasible())
    {
        if(feasible_best==nullptr)
        {
            feasible_best = new_best->clone();
        }
        else
        {
            *feasible_best = *new_best;
        }
    }
    printSearchstats(this);
}

emili::Solution* emili::SearchStatus::getBestSolution()
{
    return best;
}

emili::Solution* emili::SearchStatus::getFeasibleBestSolution()
{    
    return feasible_best;
}

void emili::SearchStatus::resetCounters()
{
    counter = 0;
    total_counter = 0;
    not_improved = 0;
}

void emili::SearchAlgorithm::setBest(Solution* newbest)
{
    status->newBestSolution(newbest);
}


/**
 * LocalSearch base class ( Old neighborhood concept)
 */

emili::Solution* emili::LocalSearch::search()
{
    neighbh->reset();
    emili::Solution* current = init->generateSolution();
    printSolstats(current);
    emili::Solution* sol = search(current);
    //if(current!=sol)
    //    delete current;

    return sol;
}

emili::Solution* emili::LocalSearch::timedSearch(float time_seconds)
{
    neighbh->reset();
    setTimer(time_seconds);
    beginTime = clock();    
    emili::Solution* current = init->generateSolution();
    emili::Solution* sol = search(current);
    //if(current!=sol)
    //    delete current;
    stopTimer();
    return sol;
}

emili::Solution* emili::LocalSearch::search(emili::Solution* initial)
{
        termcriterion->reset();
        neighbh->reset();        
        emili::Solution* newSolution = init->generateEmptySolution();
        *newSolution = *initial;
        do
        { 
            newSolution = neighbh->step(bestSoFar);
            if(bestSoFar->operator >(*newSolution))
            {

                *bestSoFar = *newSolution;
            }
            else
            {
                delete newSolution;
            }

        }while(!termcriterion->terminate(bestSoFar,newSolution));

        return bestSoFar->clone();
}

emili::Solution* emili::LocalSearch::timedSearch(float time_seconds, Solution *initial)
{   
    setTimer(time_seconds);
    beginTime = clock();    
    emili::Solution* s = search(initial);
    stopTimer();
    return s;
}

emili::Solution* emili::LocalSearch::timedSearch()
{
    setTimer(seconds);
    beginTime = clock();
    emili::Solution* current = init->generateSolution();
    emili::Solution* sol = search(current);
   // if(current!=sol)
   //     delete current;
    stopTimer();
    return sol;
}

emili::Solution* emili::LocalSearch::timedSearch(Solution *initial)
{
    setTimer(seconds);
    beginTime = clock();    
    emili::Solution* sol = search(initial);
    stopTimer();
    return sol;
}

float emili::LocalSearch::getSearchTime()
{
    return this->seconds;
}

void emili::LocalSearch::setSearchTime(float time)
{
#ifdef NOSIG
    if(time > 0)
    {
    emili::TimedTermination* tt = new emili::TimedTermination(time);
    delete termcriterion;
    termcriterion = tt;
    std::cout << "timer set " << time << " seconds " << std::endl;
    }
#endif
    this->seconds = time;
}

emili::LocalSearch::~LocalSearch()
 {
    if(init != nullptr)
    {
        delete init;
    }
    if(termcriterion != nullptr)
    {
       delete termcriterion;
    }
    if(neighbh != nullptr)
    {
       delete neighbh;
    }
      delete bestSoFar;
}

emili::Termination& emili::LocalSearch::getTermination()
{
    return *this->termcriterion;
}

emili::Neighborhood& emili::LocalSearch::getNeighborhood()
{
    return *this->neighbh;
}

emili::InitialSolution& emili::LocalSearch::getInitialSolution()
{
    return *this->init;
}

/**
 * Empty local search
 */
emili::Solution* emili::EmptyLocalSearch::search()
{
    bestSoFar = init->generateSolution();
    return bestSoFar;
}

emili::Solution* emili::EmptyLocalSearch::timedSearch(int seconds)
{
    setTimer(seconds);
    beginTime = clock();    
    bestSoFar = init->generateSolution();
    stopTimer();
    return bestSoFar;
}

emili::Solution* emili::EmptyLocalSearch::timedSearch()
{
    return timedSearch(this->getSearchTime());
}

/**
 * Feasible Local Search stuff
 */

emili::Solution* emili::LocalSearch::getBestSoFar()
{   
    return status->getBestSolution();
}

void emili::LocalSearch::setBestSoFar(Solution *newBest)
{
    setBest(newBest);
}

/**
 * Best improvement local search
 */
emili::Solution* emili::BestImprovementSearch::search(emili::Solution* initial)
{
        termcriterion->reset();
        neighbh->reset();
        emili::Solution* incumbent = initial->clone();
        emili::Solution* ithSolution;
        Neighborhood::NeighborhoodIterator end = neighbh->end();
        do
        {                       
            *bestSoFar = *incumbent;            
            Neighborhood::NeighborhoodIterator iter = neighbh->begin(bestSoFar);
            ithSolution = *iter;
            for(;iter!=end;++iter)
            {
                if(incumbent->operator >( *ithSolution)){                    
                    *incumbent = *ithSolution;
                    printSolstats(incumbent);
                }                
            }
            delete ithSolution;
        }while(!termcriterion->terminate(bestSoFar,incumbent));
        if(*bestSoFar > *incumbent)
        {
            *bestSoFar = *incumbent;
        }
        delete incumbent;
        return bestSoFar->clone();
}

/**
 * TieBraking Best Improvement local search
 */
emili::Solution* emili::TieBrakingBestImprovementSearch::search(emili::Solution* initial)
{
        termcriterion->reset();
        neighbh->reset();
        emili::Solution* incumbent = initial->clone();
        emili::Solution* ithSolution;
        Neighborhood::NeighborhoodIterator end = neighbh->end();
        do
        {
            *bestSoFar = *incumbent;
            Neighborhood::NeighborhoodIterator iter = neighbh->begin(bestSoFar);
            ithSolution = *iter;
            for(;iter!=end;++iter)
            {
                if(incumbent->operator >( *ithSolution)){
                    *incumbent = *ithSolution;
                    printSolstats(incumbent);
                }else if(incumbent->operator ==( *ithSolution))// if the two solution have the same cost
                {
                    //Compare the two solution using the cost for another problem
                   if(tiebraker.calcObjectiveFunctionValue(*incumbent) >
                           tiebraker.calcObjectiveFunctionValue((*ithSolution)))
                   {

                       *incumbent = *ithSolution;
                       printSolstats(incumbent);
                   }
                }
            }
            delete ithSolution;
        }while(!termcriterion->terminate(bestSoFar,incumbent));
        if(*bestSoFar > *incumbent)
        {
            *bestSoFar = *incumbent;
        }
        delete incumbent;
        return bestSoFar->clone();
}

/**
 * Feasible Best improvement local search
 */
emili::Solution* emili::FeasibleBestImprovementSearch::search(emili::Solution* initial)
{
        termcriterion->reset();
        neighbh->reset();
        emili::Solution* incumbent = initial->clone();
        emili::Solution* ithSolution;
        Neighborhood::NeighborhoodIterator end = neighbh->end();
        do
        {
            //*bestSoFar = *incumbent;
            setBest(incumbent);
            Neighborhood::NeighborhoodIterator iter = neighbh->begin(bestSoFar);
            ithSolution = *iter;
            for(;iter!=end;++iter)
            {
                if(incumbent->operator >( *ithSolution)){
                    *incumbent = *ithSolution;
                    printSolstats(incumbent);
                }
            }
            delete ithSolution;
        }while(!termcriterion->terminate(bestSoFar,incumbent));
        if(*bestSoFar > *incumbent)
        {
            *bestSoFar = *incumbent;
        }
        delete incumbent;
        return bestSoFar->clone();
}

/**
 * First improvement local search
 */
emili::Solution* emili::FirstImprovementSearch::search(emili::Solution* initial)
{
        termcriterion->reset();
        neighbh->reset();
        emili::Solution* incumbent = initial->clone();
        emili::Solution* ithSolution;

        Neighborhood::NeighborhoodIterator end = neighbh->end();
        do{

            *bestSoFar = *incumbent;            
            Neighborhood::NeighborhoodIterator iter = neighbh->begin(incumbent);
            ithSolution = *iter;            

            for(;iter!=end;++iter)
            {               
                if(incumbent->operator >(*ithSolution)){
                    *incumbent=*ithSolution;
                    printSolstats(incumbent);
                    break;
                }
            }
            delete ithSolution;
        }while(!termcriterion->terminate(bestSoFar,incumbent));
        if(*bestSoFar > *incumbent)
        {
            *bestSoFar = *incumbent;
        }
        delete incumbent;
        return bestSoFar->clone();
}

/**
 * Tie Braking First improvement local search
 */
emili::Solution* emili::TieBrakingFirstImprovementSearch::search(emili::Solution* initial)
{
        termcriterion->reset();
        neighbh->reset();
        emili::Solution* incumbent = initial->clone();
        emili::Solution* ithSolution;

        //bestSoFar->setSolutionValue(bestSoFar->getSolutionValue()+1);
        Neighborhood::NeighborhoodIterator end = neighbh->end();
        do{

            *bestSoFar = *incumbent;
            Neighborhood::NeighborhoodIterator iter = neighbh->begin(incumbent);
            ithSolution = *iter;

            for(;iter!=end;++iter)
            {
                if(incumbent->operator >(*ithSolution)){
                    *incumbent=*ithSolution;
                    printSolstats(incumbent);
                    break;
                }else if(incumbent->operator ==( *ithSolution))// if the two solution have the same cost
                {
                    //Compare the two solution using the cost for another problem
                   if(tiebraker.calcObjectiveFunctionValue(*incumbent) >
                           tiebraker.calcObjectiveFunctionValue((*ithSolution)))
                   {

                       *incumbent = *ithSolution;
                       printSolstats(incumbent);
                       break;
                   }
                }
            }
            delete ithSolution;
        }while(!termcriterion->terminate(bestSoFar,incumbent));
        if(*bestSoFar > *incumbent)
        {
            *bestSoFar = *incumbent;
        }
        delete incumbent;
        return bestSoFar->clone();
}

/**
 * Feasible First improvement local search
 */
emili::Solution* emili::FeasibleFirstImprovementSearch::search(emili::Solution* initial)
{
        termcriterion->reset();
        neighbh->reset();
        emili::Solution* incumbent = initial->clone();
        emili::Solution* ithSolution;

        //bestSoFar->setSolutionValue(bestSoFar->getSolutionValue()+1);
        Neighborhood::NeighborhoodIterator end = neighbh->end();
        do{

            setBest(incumbent);
            Neighborhood::NeighborhoodIterator iter = neighbh->begin(incumbent);
            ithSolution = *iter;

            for(;iter!=end;++iter)
            {
                if(incumbent->operator >(*ithSolution)){
                    *incumbent=*ithSolution;
                    printSolstats(incumbent);
                    break;
                }
            }
            delete ithSolution;
        }while(!termcriterion->terminate(bestSoFar,incumbent));
        if(*bestSoFar > *incumbent)
        {
            *bestSoFar = *incumbent;
        }
        delete incumbent;
        return bestSoFar->clone();
}

/** OLD

emili::Solution* emili::FirstImprovementSearch::search(emili::Solution* initial)
{
        termcriterion->reset();
        neighbh->reset();
        bestSoFar = init->generateEmptySolution();
        emili::Solution* incumbent = initial->clone();
        emili::Solution* ithSolution;

        //bestSoFar->setSolutionValue(bestSoFar->getSolutionValue()+1);

        do{

            *bestSoFar = *incumbent;
            for(Neighborhood::NeighborhoodIterator iter = neighbh->begin(incumbent);iter!=neighbh->end();++iter)
            {
                 ithSolution = *iter;
                if(incumbent->operator >(*ithSolution)){
                    *incumbent = *ithSolution;
                    break;
                }
            }
            delete ithSolution;
        }while(!termcriterion->terminate(bestSoFar,incumbent));
        delete incumbent;
        return bestSoFar;
}
*/



/**
 * TABU SEARCH
 */
emili::Solution* emili::BestTabuSearch::search()
{
    tabuMemory.reset();
     return LocalSearch::search();
}


emili::Solution* emili::BestTabuSearch::search(emili::Solution *initial)
{
    termcriterion->reset();
    neighbh->reset();
    emili::Solution* incumbent = initial->clone();    
    emili::Solution* ithSolution = nullptr;
    do
    {
        if(*bestSoFar > *incumbent){
            *bestSoFar = *incumbent;            
        }

        Neighborhood::NeighborhoodIterator iter = neighbh->begin(bestSoFar);
        if(iter!=neighbh->end())
        {
           ithSolution = *iter;
           *incumbent = *ithSolution;


        for(;iter!=neighbh->end();++iter)
        {
            ithSolution = *iter;            
            if(*incumbent > *ithSolution && (tabuMemory.tabu_check(ithSolution) || *ithSolution < *bestSoFar))                
            {
                tabuMemory.registerMove(incumbent,ithSolution);// make the tabu memory record the move used on incumbent to generate ithSolution
                *incumbent = *ithSolution;
                printSolstats(incumbent);
            }
        }
        delete ithSolution;
        tabuMemory.forbid(incumbent);
        }
    }while(!termcriterion->terminate(bestSoFar,incumbent));
    if(*bestSoFar > *incumbent)
    {
        *bestSoFar = *incumbent;
    }
    delete incumbent;
    return bestSoFar->clone();

}


emili::Solution* emili::FirstTabuSearch::search(emili::Solution *initial)
{
    termcriterion->reset();
    neighbh->reset();
    emili::Solution* incumbent = initial->clone();
    emili::Solution* ithSolution;
        do{
            if(*bestSoFar > *incumbent)
                *bestSoFar = *incumbent;

            Neighborhood::NeighborhoodIterator iter = neighbh->begin(bestSoFar);
            if(iter!=neighbh->end())
            {
               ithSolution = *iter;
               *incumbent = *ithSolution;


            for(;iter!=neighbh->end();++iter)
            {
                ithSolution = *iter;                
                if(incumbent->operator >(*ithSolution)&& (tabuMemory.tabu_check(ithSolution) || *ithSolution < *bestSoFar)){                    
                    tabuMemory.registerMove(incumbent,ithSolution);
                    *incumbent = *ithSolution;
                    printSolstats(incumbent);
                    break;
                }
            }
         delete ithSolution;         
        tabuMemory.forbid(incumbent);
        }
    }while(!termcriterion->terminate(bestSoFar,incumbent));
    delete incumbent;
    return bestSoFar->clone();
}
/**
emili::Solution* emili::TabuSearch::search(emili::Solution *initial)
{
    termcriterion->reset();
    neighbh->reset();
    emili::Solution* current = init->generateEmptySolution();
    emili::Solution* newSolution = init->generateEmptySolution();
    *newSolution = *initial;
    emili::Solution* newS = nullptr;
    do
    {
       if(newSolution!=current)
       {
        delete current;
       }
        current = newSolution;
        neighbh->reset();
        Solution* best = neighbh->step(current);
        double bestSol = best->getSolutionValue();

        for(Neighborhood::NeighborhoodIterator iter = neighbh->begin(current);iter!=neighbh->end();++iter)
        {
            newS = *iter;
            double newSol = newS->getSolutionValue();
            if(bestSol > newSol){
                if(tabuMemory.tabu_check(newS))//<- Aspiration goes here.
                {
                    best = newS;
                    bestSol = newSol;
                }
            }

        }
        newSolution = best;
        tabuMemory.forbid(newSolution);
    }while(!termcriterion->terminate(current,newSolution));
    return current;
}


*/

emili::Solution* emili::RandomMovePerturbation::perturb(Solution *solution)
{
    Solution* ret = explorer.random(solution);

    for (int var = 1; var < numberOfSteps; ++var) {
        Solution* temp = ret;
        ret = explorer.random(temp);
        delete temp;
    }

    return ret;
}

emili::Solution* emili::VNRandomMovePerturbation::perturb(Solution *solution)
{

    Solution* ret = explorers[currentExplorer]->random(solution);
    for (int var = 1; var < numberOfSteps; ++var) {
        Solution* temp = ret;
        ret = explorers[currentExplorer]->random(temp);
        delete temp;
    }

    if(!(currentIteration <= numberOfIterations))
    {
        currentIteration=0;
        currentExplorer = (currentExplorer+1)%explorers.size();
    }
    else
    {
        currentIteration++;
    }

    return ret;
}

emili::Solution* emili::RandomPerturbationSet::perturb(Solution *solution)
{
    int p = emili::generateRandomNumber()%size;
    return perturbations[p]->perturb(solution);
}

emili::Solution* emili::ComplexPerturbation::perturb(Solution *solution)
{
    Solution* div = p->perturb(solution);
    Solution* inte = ls->search(div);
    delete div;
    return inte;
}

emili::Solution* emili::MRSILSPerturbation::perturb(Solution *solution)
{
    emili::Solution* toper = solution;
    int psize = solution_pool.size();
    bool toinsert = true;
    //Check if solution is already present in solution_pool
    for(int i = 0; i< psize; i++)
    {
        if(*solution_pool[i] == *solution)
        {
            toinsert = false;
            //if already in pool stop loop
            break;
        }
    }
    if(toinsert)
    {
        //add solution to solution_poll
        solution_pool.push_back(solution->clone());
        //update worst solution
        if(*solution_pool[worst] < *solution)
        {
            worst = psize;
        }
        //check if pool is full
        if(psize == pool_size)
        {
            //Delete the worst
            std::swap(solution_pool[worst],solution_pool[psize]);
            //Free Solution memory
            delete solution_pool[psize];
            //Clean solution_pool
            solution_pool.pop_back();
            //Update worst
            for(int i = 0; i< psize; i++)
            {
                if(*solution_pool[i] > *solution_pool[worst])
                {
                    worst = i;
                }
            }
        }
        else
        {
            psize++;
        }
    }
    if(psize >= pool_size)
    {
        //Take a random solution from the pool for the perturbation
        int k = emili::generateRandomNumber()%psize;
        toper = solution_pool[k];
    }

    return p->perturb(toper);
}

emili::Solution* emili::AlwaysAccept::accept(Solution *intensification_solution, Solution *diversification_solution)
{
    if(acc==ACC_DIVERSIFICATION)
    {
        return diversification_solution;
    }
    else
    {
        return intensification_solution;
    }
}

emili::Solution* emili::AcceptImproveEqual::accept(Solution *intensification_solution, Solution *diversification_solution)
{
    if(*diversification_solution <= *intensification_solution)
    {
        return diversification_solution;
    }
    return intensification_solution;
}

emili::Solution* emili::ImproveAccept::accept(Solution *intensification_solution, Solution *diversification_solution)
{
    Solution* k = intensification_solution;
    if(intensification_solution->operator >(*diversification_solution)){
        k = diversification_solution;
    }
    return k;
}

emili::Solution* emili::AcceptPlateau::accept(Solution *intensification_solution, Solution *diversification_solution)
{
    if(*diversification_solution <= *intensification_solution)
    {
        return diversification_solution;
        this->current_step=0;
        this->threshold_status=0;
    }
    else
    {
        threshold_status++;
        if(threshold_status >= this->plateau_threshold)
        {
            if(current_step <= max_plateau_steps)
            {
                current_step++;
                return diversification_solution;
            }
            else
            {
                threshold_status = 0;
                current_step = 0;
            }

        }
    }
    return intensification_solution;
}

emili::Solution* emili::AcceptExplore::accept(Solution *intensification_solution, Solution *diversification_solution)
{
    if(*diversification_solution < *intensification_solution)
    {
        iteration = 0;
    }
    else if(iteration < k)
    {
        iteration++;
    }
    else
    {
        iteration=0;
        Solution* bestSoFar = emili::getAlgo()->getBestSoFar();
        *diversification_solution = *bestSoFar;
    }
    return diversification_solution;
}


/**
 * Iterated Local Search
 */

emili::Solution* emili::IteratedLocalSearch::search(){
    termcriterion->reset();
    acc.reset();
    Solution* current = init->generateSolution();
    printSolstats(current);
    Solution* ret = search(current);
    delete current;
    return ret;
}

emili::Solution* emili::IteratedLocalSearch::search(emili::Solution* initial){
    termcriterion->reset();
    acc.reset();        
    emili::Solution* s = ls.search(initial);
    *bestSoFar = *s;
    emili::Solution* s_s = nullptr;
    emili::Solution* s_p = nullptr;
    //initialization done
    do{
        if(s_p != s && s_p != nullptr)
            delete s_p;
        //Perturbation step
        s_p = pert.perturb(s);
        //local search on s_p
        if(s!=s_s && s_s != nullptr)
            delete s_s;
        s_s = ls.search(s_p);
        delete s_p;
        //best solution
        if(*s_s < *bestSoFar)
        {
            *bestSoFar = *s_s;
         printSolstats(bestSoFar);
            //s_time = clock();
        }
        //acceptance step
        s_p = s;
        s = acc.accept(s_p,s_s);
    }while(!termcriterion->terminate(s_p,s));
    delete s_p;
    delete s_s;
    return bestSoFar->clone();
}

emili::Solution* emili::IteratedLocalSearch::timedSearch(float maxTime)
{
        termcriterion->reset();
        acc.reset();
        setTimer(maxTime);
        /**
            search start
        */
        beginTime = clock();        
        emili::Solution*  s = ls.search();
        *bestSoFar = *s ;
        emili::Solution* s_s = nullptr;
        emili::Solution* s_p = nullptr;
        //initialization done
        do{
            if(s_p != s && s_p != nullptr)
                delete s_p;
            //Perturbation step
            s_p = pert.perturb(s);
            //local search on s_p
            if(s!=s_s && s_s != nullptr)
                delete s_s;
            s_s = ls.search(s_p);
            delete s_p;
            //best solution
            if(*s_s < *bestSoFar)
            {
                *bestSoFar = *s_s;
                printSolstats(bestSoFar);
                //s_time = clock();
            }
            //acceptance step
            s_p = s;
            s = acc.accept(s_p,s_s);
        }while(!termcriterion->terminate(s_p,s) && keep_going);
        delete s_p;
        delete s_s;
        stopTimer();
        return bestSoFar->clone();
}

emili::Solution* emili::IteratedLocalSearch::timedSearch(float maxTime,emili::Solution* initial)
{
        termcriterion->reset();
        acc.reset();
        setTimer(maxTime);
        /**
            search start
        */
        beginTime = clock();                
        emili::Solution* s = ls.search(initial);
        *bestSoFar = *s;
        emili::Solution* s_s = nullptr;
        emili::Solution* s_p = nullptr;
        //initialization done
        do{
            if(s_p != s && s_p != nullptr)
                delete s_p;
            //Perturbation step
            s_p = pert.perturb(s);
            //local search on s_p
            if(s!=s_s && s_s != nullptr)
                delete s_s;
            s_s = ls.search(s_p);
            delete s_p;
            //best solution
            if(*s_s < *bestSoFar)
            {
                *bestSoFar = *s_s;             
                printSolstats(bestSoFar);
                //s_time = clock();
            }
            //acceptance step
            s_p = s;
            s = acc.accept(s_p,s_s);
        }while(!termcriterion->terminate(s_p,s) && keep_going);
        delete s_p;
        delete s_s;
        stopTimer();
        return bestSoFar->clone();
}

emili::Solution* emili::IteratedLocalSearch::getBestSoFar()
{
    emili::Solution* bestOfInnerLocal = this->ls.getBestSoFar();

    if(bestOfInnerLocal != nullptr &&  bestOfInnerLocal->operator <(*bestSoFar))
    {
        return bestOfInnerLocal;
    }
    return bestSoFar;
}

/**
 * Feasible Iterated Local Search
 */

emili::Solution* emili::FeasibleIteratedLocalSearch::search(){
    termcriterion->reset();
    acc.reset();
    Solution* current = init->generateSolution();
    printSolstats(current);
    Solution* ret = search(current);
    delete current;
    return ret;
}

emili::Solution* emili::FeasibleIteratedLocalSearch::search(emili::Solution* initial){
    termcriterion->reset();
    acc.reset();
    emili::Solution* s = ls.search(initial);
    //*bestSoFar = *s;
    setBest(s);
    emili::Solution* s_s = nullptr;
    emili::Solution* s_p = nullptr;
    //initialization done
    do{
        if(s_p != s && s_p != nullptr)
            delete s_p;
        //Perturbation step
        s_p = pert.perturb(s);
        //local search on s_p
        if(s!=s_s && s_s != nullptr)
            delete s_s;
        s_s = ls.search(s_p);
        delete s_p;
        //best solution
        if(*s_s < *bestSoFar)
        {
            //*bestSoFar = *s_s;
            setBest(s_s);
         printSolstats(bestSoFar);
            //s_time = clock();
        }
        //acceptance step
        s_p = s;
        s = acc.accept(s_p,s_s);
    }while(!termcriterion->terminate(s_p,s));
    delete s_p;
    delete s_s;
    return bestSoFar->clone();
}

emili::Solution* emili::FeasibleIteratedLocalSearch::timedSearch(float maxTime)
{
        termcriterion->reset();
        acc.reset();
        setTimer(maxTime);
        /**
            search start
        */
        beginTime = clock();
        emili::Solution*  s = ls.search();
        //*bestSoFar = *s;
        setBest(s);
        emili::Solution* s_s = nullptr;
        emili::Solution* s_p = nullptr;
        //initialization done
        do{
            if(s_p != s && s_p != nullptr)
                delete s_p;
            //Perturbation step
            s_p = pert.perturb(s);
            //local search on s_p
            if(s!=s_s && s_s != nullptr)
                delete s_s;
            s_s = ls.search(s_p);
            delete s_p;
            //best solution
            if(*s_s < *bestSoFar)
            {
                //*bestSoFar = *s_s;
                setBest(s_s);
                printSolstats(bestSoFar);
                //s_time = clock();
            }
            //acceptance step
            s_p = s;
            s = acc.accept(s_p,s_s);
        }while(!termcriterion->terminate(s_p,s) && keep_going);
        delete s_p;
        delete s_s;
        stopTimer();
        return bestSoFar->clone();
}

emili::Solution* emili::FeasibleIteratedLocalSearch::timedSearch(float maxTime,emili::Solution* initial)
{
        termcriterion->reset();
        acc.reset();
        setTimer(maxTime);
        /**
            search start
        */
        beginTime = clock();
        emili::Solution* s = ls.search(initial);
        //*bestSoFar = *s;
        setBest(s);
        emili::Solution* s_s = nullptr;
        emili::Solution* s_p = nullptr;
        //initialization done
        do{
            if(s_p != s && s_p != nullptr)
                delete s_p;
            //Perturbation step
            s_p = pert.perturb(s);
            //local search on s_p
            if(s!=s_s && s_s != nullptr)
                delete s_s;
            s_s = ls.search(s_p);
            delete s_p;
            //best solution
            if(*s_s < *bestSoFar)
            {
                //*bestSoFar = *s_s;
                setBest(s_s);
                printSolstats(bestSoFar);
                //s_time = clock();
            }
            //acceptance step
            s_p = s;
            s = acc.accept(s_p,s_s);
        }while(!termcriterion->terminate(s_p,s) && keep_going);
        delete s_p;
        delete s_s;
        stopTimer();
        return bestSoFar->clone();
}

emili::Solution* emili::FeasibleIteratedLocalSearch::getBestSoFar()
{
    emili::Solution* bestOfInnerLocal = IteratedLocalSearch::getBestSoFar();
    emili::Solution* feasibleBest = status->getFeasibleBestSolution();
    if(feasibleBest != nullptr )
    {
        if(bestOfInnerLocal->isFeasible() && bestOfInnerLocal->operator <(*feasibleBest))
        {
            return bestOfInnerLocal;
        }
        return feasibleBest;
    }
    return bestOfInnerLocal;

}

/**
 * Never stopping termination
 */

bool emili::WhileTrueTermination::terminate(Solution* currentSolution, Solution* newSolution)
{
    if(keep_going)
    {
        return false;
    }
    else
    {
        return true;
    }
}

/**
 * Timed termination
 */

bool emili::TimedTermination::terminate(Solution *currentSolution, Solution *newSolution)
{
    clock_t test = clock();
    float time = (test-start)/ (float)CLOCKS_PER_SEC;
#ifdef NOSIG
    if(time > secs)
    {
        finalise(2);
    }
#endif
    return !(time < secs);
}

void emili::TimedTermination::reset()
{

    //setTerminationTimer(this->secs);
   /**  if(secs<0)
    {
    if(max_time > 0)
        secs = (int)max_time*_ratio;
    else
        secs = (int)30*_ratio;
     }*/
    secs = _ratio;
    start = clock();
#ifdef NOSIG
    beginTime = start;
#endif
}

/** emili::Solution* emili::VNDSearch::searchOneNeigh(Solution *initial, emili::Neighborhood* n)
{
    termcriterion->reset();
    n->reset();
    emili::Solution* bestSoFar = init->generateEmptySolution();
    emili::Solution* incumbent = init->generateEmptySolution();
    *incumbent = *initial;
    emili::Solution* ithSolution = nullptr;

    do
    {
        delete bestSoFar;
        bestSoFar = incumbent;
        n->reset();
        Solution* bestOfTheIteration = incumbent;

        for(Neighborhood::NeighborhoodIterator iter = n->begin(bestSoFar);iter!=n->end();++iter)
        {
            ithSolution = *iter;

            if(bestOfTheIteration->operator >( *ithSolution)){
                bestOfTheIteration = ithSolution;

            }

        }
        incumbent = bestOfTheIteration;

    }while(!termcriterion->terminate(bestSoFar,incumbent));
    return bestSoFar;
}*/

/**
 * LocalMinima terminations
 */
bool emili::LocalMinimaTermination::terminate(Solution* currentSolution,Solution* newSolution)
{
   //std::cout << currentSolution.getSolutionValue() << " <= " << newSolution.getSolutionValue()<<std::endl;
    if(newSolution == nullptr)
    {
        return true;
    }
    else
    {
        return currentSolution->operator <=(*newSolution);
    }
}

/**
 * MaxSteps Termination
 */
bool emili::MaxStepsTermination::terminate(Solution *currentSolution, Solution *newSolution)

{
    if(current_step >= max_steps_){
        if(*currentSolution > *newSolution)
        {
            *currentSolution = *newSolution;
        }
        return true;
    }
    else
    {
        current_step++;
        return false;
    }
}

/**
 * MaxSteps or Local minimum termination
 */

bool emili::MaxStepsOrLocmin::terminate(Solution *currentSolution, Solution *newSolution)
{
        if(current_step >= max_steps_){
            if(*currentSolution > *newSolution)
            {
                *currentSolution = *newSolution;
            }
            return true;
        }
        else
        {
            current_step++;
            return currentSolution->operator <=(*newSolution);
        }
}

void emili::MaxStepsTermination::reset()
{
    current_step = 0;
}

/**
 * Piped Local Search
*/

emili::Solution* emili::PipeSearch::search(Solution *initial)
{
    Solution* current = init->generateEmptySolution();
    *bestSoFar = *initial;
    *current  = *initial;
    for(std::vector< emili::LocalSearch*>::iterator iter = lss.begin();iter!=lss.end();++iter)
    {
        Solution* ithSolution = (*iter)->search(current);
        *current = *ithSolution;
        if(current->operator <(*bestSoFar))
        {            
            *bestSoFar = *ithSolution;
        }
        delete ithSolution;
    }
    delete current;
    return bestSoFar;
}

/**
 *  METROPOLIS ACCEPTANCE
 */

float emili::MetropolisAcceptance::getTemp()
{
    return temperature;
}

void emili::MetropolisAcceptance::setTemp(float temp)
{
    this->temperature = temp;
}

emili::Solution* emili::MetropolisAcceptance::accept(Solution *intensification_solution, Solution *diversification_solution)
{
    float intens = intensification_solution->getSolutionValue();
    float divers = diversification_solution->getSolutionValue();
    if(diversification_solution->operator >(*intensification_solution))
    {
        float prob = std::exp((intens-divers)/temperature);
        if(prob < 1.0 && generateRealRandomNumber()>prob)
        {
            return intensification_solution;
        }
    }
    return diversification_solution;
}

emili::Solution* emili::Metropolis::accept(Solution *intensification_solution, Solution *diversification_solution)
{
    if(counter == interval && temperature > end_temp)
    {     
        temperature = (alpha * temperature) - beta;
        counter=0;
    }    
    counter++;
    float intens = intensification_solution->getSolutionValue();
    float divers = diversification_solution->getSolutionValue();
    if(diversification_solution->operator >(*intensification_solution))
    { 
        float prob = std::exp((intens-divers)/temperature);
        if(prob < 1.0 && generateRealRandomNumber()>prob)
        {
            return intensification_solution;
        }
    }

    return diversification_solution;
}

void emili::Metropolis::reset()
{
    temperature = start_temp;
    counter = 0;
}

/**  GVNS */

emili::Solution* emili::LS_VND::search(emili::Solution *initial)
{
    int i = 0;
    int k = neigh.size();
    Solution* incumbent = neigh[i]->search(initial);
    *bestSoFar = *incumbent;
    do{

        Solution* new_s = neigh[i]->search(incumbent);
        if(*new_s < *incumbent)
        {
            delete incumbent;
            incumbent = new_s;
            if(incumbent->operator <(*bestSoFar))
            {
                *bestSoFar = *incumbent;
            }
            i = 0;
        }
        else
        {
            i = i+1;
            delete new_s;
        }
    }while(i < k);
    return bestSoFar->clone();
}


emili::Solution* emili::Shake::perturb(Solution *solution)
{
        return shake(solution,emili::generateRandomNumber()%Kmax);
}

emili::Solution* emili::PerShake::shake(Solution *s, int n)
{
    Perturbation* p = shakes[n];
    return p->perturb(s);
}

emili::Solution* emili::NeighborhoodShake::shake(Solution *s, int n)
{
    int k = n%n_num;
    int size = n/n_num + 1 ;
    //std::cout << "N " << n << " K " << k << " S " << size << std::endl;
    Neighborhood* p = shakes[k];
    return p->random(s,size);
}

emili::Solution* emili::NeighborhoodChange::accept(Solution *intensification_solution, Solution *diversification_solution)
{
        int n = 1;
        return neighborhoodChange(intensification_solution,diversification_solution,n);
}

emili::Solution* emili::AccNeighborhoodChange::neighborhoodChange(emili::Solution *intensification_solution,emili::Solution *diversification_solution, int &n)
{
    Solution* accepted = acc->accept(intensification_solution,diversification_solution);
    if(accepted == intensification_solution)
    {
        n++;
    }
    else
    {
        n=0;
    }
    return accepted;
}


emili::Solution* emili::GVNS::search(Solution* initial)
{
        int k = 0;
        int k_max = shaker.getKmax();
        termcriterion->reset();
        changer.reset();
        changer.setKmax(k_max);
        emili::Solution* s = this->init->generateEmptySolution();
        *s = *initial;
        *bestSoFar = *s;
        emili::Solution* s_s = nullptr;
        emili::Solution* s_p = nullptr;
        //initialization done
        do{
            k = 0;
            do{              
                if(s_p != s && s_p != nullptr)
                    delete s_p;
                //Shake step
               // std::cout << "pre shake"<< k << std::endl;
                s_p = shaker.shake(s,k);
                //std::cout << "post shake" << std::endl;
                //local search on s_p
                if(s!=s_s && s_s != nullptr)
                    delete s_s;
               // std::cout << "pre ls" << std::endl;
                s_s = ls.search(s_p);
               // std::cout << "post ls" << std::endl;
                delete s_p;
                //best solution
                if(*s_s < *bestSoFar)
                {
                    *bestSoFar = *s_s;
                    printSolstats(bestSoFar);
                    //s_time = clock();
                }
                //neighborhood change step
                s_p = s;
               // std::cout << "pre nc " << std::endl;
                s = changer.neighborhoodChange(s_p,s_s,k);                
                //std::cout << "post nc" << std::endl;
            }while (k < k_max);
        }while(!termcriterion->terminate(s_p,s));
       // std::cout << "returning" << std::endl;

        delete s_p;
        delete s_s;
        return bestSoFar->clone();

}

emili::Solution* emili::GVNS::getBestSoFar()
{
    emili::Solution* bestOfInnerLocal = this->ls.getBestSoFar();

    if(bestOfInnerLocal != nullptr &&  bestOfInnerLocal->operator <(*bestSoFar))
    {
        return bestOfInnerLocal;
    }
    return bestSoFar;
}

emili::Solution* emili::ComposedInitialSolution::generateEmptySolution()
{
    return is.generateEmptySolution();
}

emili::Solution* emili::ComposedInitialSolution::generateSolution()
{
    Solution* s = is.generateSolution();
    Solution* ss = ls.search(s);
    if(s!=ss)
        delete s;
    return ss;
}

emili::Solution* emili::AlternateLocalSearch::search(emili::Solution * solution)
{
    emili::Solution* s = solution->clone();
    std::cout << "START" << std::endl;
    while (!termcriterion->terminate(nullptr, nullptr)) {
      if(turn)
      {
          //std::cout << "run first " << std::endl;
          //Solution* sol1
          s = ls1->search(s);
          //std::cout << "done first " << s->getSolutionValue() << std::endl;
          /*if (sol1 != s) {
             delete s;
             s = sol1;
          }*/
      }
      else
      {
          //std::cout << "run second " << std::endl;
          //Solution* sol2 
          s = ls2->search(s);
          /*if (sol2 != s) {
               delete s;
               s = sol2;
          }*/
          //std::cout << "done second " << s->getSolutionValue() << std::endl;
      }
      turn = !turn;
    }
    return s;
}

emili::Solution* emili::AlternateLocalSearch::getBestSoFar()
{
    emili::Solution* best = ls1->getBestSoFar();
    emili::Solution* best2 = ls2->getBestSoFar();


    /*if ( best == nullptr)
      printf("oh fuck\n");
    printf("best: %f\n", best->getSolutionValue());
    if ( best2 == nullptr)
      printf("oh fuck\n");
    printf("best2: %f\n", best2->getSolutionValue());*/

    if (best == nullptr && best2 == nullptr) {
      //std::cout << "WTF .. both null" << std::endl;
      return nullptr;
    } else if (best != nullptr && best2 == nullptr) {
      //std::cout << "WTF .. ls1 not null" << std::endl; 
      return best;
    } else if (best == nullptr && best2 != nullptr) {
      //std::cout << "WTF .. ls2 not null" << std::endl;
      return best2; // should never be here
    }
    
    /*std::cout << "WTF both not null" << std::endl;
    double xxx = best->getSolutionValue();
    std::cout << best->getSolutionRepresentation() << std::endl;
    
    std::cout << best2->getSolutionValue() << std::endl;
    std::cout << best2->getSolutionRepresentation() << std::endl;*/
    if(*best < *best2)
    {
        return best;
    }
    else
    {
        return best2;
    }

}

bool emili::MaxStepsNoImprov::terminate(Solution *currentSolution, Solution *newSolution)
{
    if(*newSolution < *currentSolution)
    {
        h = 0;
    }
    else
    {
        h++;
        if(h <= max_h)
        {
            return true;
        }
    }
    return false;
}

void emili::MaxStepsNoImprov::reset()
{
    h= 0;
}
