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
#include "pfsp/permutationflowshop.h"
/*
 * WARNING!!!
 * Adding data structures to a solution subclass could broken this method

emili::Solution& emili::Solution::operator=(const emili::Solution& a)
{
    this->instance = a.instance;
    this->setRawData(a.getRawData());
    return *this;
}
 */
/*
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


/*
 * TIMED SEARCH CODE
 */
bool print;
bool keep_going;
bool timer_keep_going;
std::ostringstream messages;
std::string lastMessage;
clock_t endTime;
clock_t beginTime;
clock_t s_time;
emili::LocalSearch* localsearch;


double emili::getCurrentExecutionTime()
{
    return (double)((clock()-beginTime)/ (double)CLOCKS_PER_SEC);
}

emili::LocalSearch* emili::getAlgo()
{
    return localsearch;
}

static void finalise (int _)
{
    keep_going = false;
    endTime = clock();
    emili::Solution* s_cap = localsearch->getBestSoFar();
    if(s_cap != nullptr)
    {
        double sol_val = s_cap->getSolutionValue();
        if(print)
        {
            messages << "CPU time: " << (endTime - beginTime) / (float)CLOCKS_PER_SEC << std::endl;
            messages << "iteration counter : " << emili::iteration_counter()<< std::endl;
            messages << "objective function value : "<< sol_val << std::endl;
            messages << "solution : " << s_cap->getSolutionRepresentation() << std::endl;
            //std::cout << "Reached at time: " << (s_time - beginTime) / (float)CLOCKS_PER_SEC << std::endl;
            //std::cerr << (endTime - beginTime) / (float)CLOCKS_PER_SEC << " ";
        }
        else
        {
            std::cout << "CPU time: " << (endTime - beginTime) / (float)CLOCKS_PER_SEC << std::endl;
            std::cout << "iteration counter : " << emili::iteration_counter()<< std::endl;
            std::cerr << sol_val << std::endl;
            std::cerr << std::flush;
        }
    }
    else
    {
        if(print)
        {
            messages << "No valid solution found!" << std::endl;
        }
        else
        {
            std::cout  << "No valid solution found!" << std::endl;

        }
    }
    //std::cout << std::flush;
    if(print)
    {
        lastMessage = messages.str();
    }
    else
    {
        exit(0);
    }
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
    std::cout << lastMessage << std::endl;
}

int max_time = -1 ;
static inline void setTimer(int maxTime)
{
    keep_going = true;
    timer.it_value.tv_sec = maxTime;
    timer.it_value.tv_usec = 0;
    timer.it_interval.tv_sec = 0;
    timer.it_interval.tv_usec = 0;
    emili::iteration_counter_zero();
    signal(SIGPROF, finalise);
    signal(SIGINT, finalise);
    if (setitimer (ITIMER_PROF, &timer, NULL) != 0) {
        printf("error in setitimer\n");
        exit(10);
    }else{
        std::cout << "timer set " << maxTime << " seconds " << std::endl;
        if(print)
            atexit(lastPrint);
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

static inline void setTimer(int maxTime)
{
    //std::cout << "timer set"
}
#endif

/*
 * Timer HOOK
 */

/*
 * Iteration counter
 */
static unsigned long iteration_counter_ ;

void emili::iteration_counter_zero()
{
    iteration_counter_ = 0;
}


int emili::iteration_counter(){
    return iteration_counter_;
}

void emili::iteration_increment(){
    iteration_counter_++;
}

void emili::iteration_decrement(){
    iteration_counter_--;
}

/*
 * Print Solution info
 */

inline void emili::printSolstats(emili::Solution* sol)
{
#ifdef WITH_STATS
    if(print)
    {
      std::cout << (clock() - beginTime) / (float)CLOCKS_PER_SEC << " , " << sol->getSolutionValue() << " , " << iteration_counter_ << "\n";
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


/*
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

/*
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

/*
 * LocalSearch base class ( Old neighborhood concept)
 */

emili::Solution* emili::LocalSearch::search()
{
    neighbh->reset();
    emili::Solution* current = init->generateSolution();
    emili::Solution* sol = search(current);
    if(current!=sol)
        delete current;

    return sol;
}

emili::Solution* emili::LocalSearch::timedSearch(int time_seconds)
{
    neighbh->reset();
    emili::Solution* current = init->generateSolution();
    emili::Solution* sol = timedSearch(time_seconds,current);
    if(current!=sol)
        delete current;

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

emili::Solution* emili::LocalSearch::timedSearch(int time_seconds, Solution *initial)
{

    setTimer(time_seconds);
    beginTime = clock();
    localsearch = this;
    emili::Solution* s = search(initial);
    stopTimer();
    return s;
}

emili::Solution* emili::LocalSearch::timedSearch()
{
    neighbh->reset();
    emili::Solution* current = init->generateSolution();
    emili::Solution* sol = timedSearch(seconds,current);
    if(current!=sol)
        delete current;

    return sol;
}

emili::Solution* emili::LocalSearch::timedSearch(Solution *initial)
{
    neighbh->reset();
    emili::Solution* sol = timedSearch(seconds,initial);
    return sol;
}

int emili::LocalSearch::getSearchTime()
{
    return this->seconds;
}

void emili::LocalSearch::setSearchTime(int time)
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




/*
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
                }                
            }
            delete ithSolution;
        }while(!termcriterion->terminate(bestSoFar,incumbent));
        delete incumbent;
        return bestSoFar->clone();
}


/*
 * First improvement local search
 */
emili::Solution* emili::FirstImprovementSearch::search(emili::Solution* initial)
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
                    break;
                }                                
            }
            delete ithSolution;
        }while(!termcriterion->terminate(bestSoFar,incumbent));
        delete incumbent;
        return bestSoFar->clone();
}

/*OLD

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



/*
 * TABU SEARCH
 */
emili::Solution* emili::BestTabuSearch::search()
{
    tabuMemory.reset();
    emili::Solution* current = init->generateSolution();
    emili::Solution* sol = search(current);
    if(current!=sol)
    delete current;

    return sol;
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
            tabuMemory.registerMove(incumbent,ithSolution);// make the tabu memory record the move used on incumbent to generate ithSolution
            if(*incumbent > *ithSolution && (tabuMemory.tabu_check(ithSolution) || *ithSolution < *bestSoFar))
            {
                *incumbent = *ithSolution;
            }
        }

        delete ithSolution;
        tabuMemory.forbid(incumbent);
        }
    }while(!termcriterion->terminate(bestSoFar,incumbent));
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
                tabuMemory.registerMove(incumbent,ithSolution);
                if(incumbent->operator >(*ithSolution)&& (tabuMemory.tabu_check(ithSolution) || *ithSolution < *bestSoFar)){                    
                    *incumbent = *ithSolution;
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
/*
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


/*
 * Iterated Local Search
 */

emili::Solution* emili::IteratedLocalSearch::search(){
    termcriterion->reset();
    acc.reset();
    Solution* current = init->generateSolution();
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

emili::Solution* emili::IteratedLocalSearch::timedSearch(int maxTime)
{
        termcriterion->reset();
        acc.reset();
        localsearch = this;
        setTimer(maxTime);
        /*
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

emili::Solution* emili::IteratedLocalSearch::timedSearch(int maxTime,emili::Solution* initial)
{
        termcriterion->reset();
        acc.reset();
        setTimer(maxTime);
        localsearch = this;
        /*
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

/*
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

/*
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
   /* if(secs<0)
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

/*emili::Solution* emili::VNDSearch::searchOneNeigh(Solution *initial, emili::Neighborhood* n)
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

/*
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

/*
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

/*
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

/*
 * Piped Local Search
*/

emili::Solution* emili::PipeSearch::search(Solution *initial)
{
    Solution* bestSoFar = init->generateEmptySolution();
    bestSoFar->operator =(*initial);
    Solution* ithSolution = bestSoFar;
    for(std::vector< emili::LocalSearch*>::iterator iter = lss.begin();iter!=lss.end();++iter)
    {
        ithSolution = (*iter)->search(ithSolution);
        if(ithSolution->operator <(*bestSoFar))
        {
            delete bestSoFar;
            bestSoFar = ithSolution;
        }
    }
    return bestSoFar;
}

/*
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
        temperature = (alpha * temperature) - rate;        
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

/* GVNS */

emili::Solution* emili::GVNS::search(Solution* initial)
{
        int k = 0;
        int k_max = perturbations.size();
        bestSoFar = initial;
        ls.setBestSoFar(initial);
        emili::Solution*  s = ls.search(bestSoFar);
        *bestSoFar = *s ;

        emili::Solution* s_s = nullptr;
        //initialization done
        do{

            //Perturbation step
            emili::Solution* s_p = perturbations[k]->perturb(s);

            //local search on s_p
            if(s!=s_s && s_s != nullptr)
                delete s_s;

            //ls.setBestSoFar(bestSoFar);
            s_s = ls.search(s_p);

            //best solution
            if(*s_s < *bestSoFar)
            {
                *bestSoFar = *s_s;
                //ls.setBestSoFar(bestSoFar);
            }

            if(*s_s < *s_p)
            {
                s = s_s;
                k = 0 ;
                delete s_p;
            }
            else
            {
                delete s_p;
                k++;
            }

        }while(k < k_max && keep_going);

        return bestSoFar;
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

