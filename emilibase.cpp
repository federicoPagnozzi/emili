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

void emili::initializeRandomFromRandom() {
    initializeRandom( getRandomSeedFromRandom() );
}

int emili::getRandomSeedFromRandom() {
#ifndef NOC11
    std::random_device rd;
    return rd();
#else
    return time(0);
#endif
}

#ifndef NOC11

std::mt19937 generator;
std::uniform_int_distribution<int> distribution; // (0,maxint)
std::uniform_real_distribution<float> realdistr; // (0,1)
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
std::tr1::uniform_int<int> distribution; // (0,maxint)
std::tr1::uniform_real<float> realdistr; // (0,1)
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

int emili::generateRandRange(int n)
{
    return generateRandInt(0, n - 1);
}

int emili::generateRandRange(int from, int to)
{
    return from + generateRandRange(to - from);
}

int emili::generateRandInt(int from, int to)
{
    return std::uniform_int_distribution<>(from, to)(generator);
}

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

bool keep_going;
bool timer_keep_going;
clock_t endTime;
clock_t beginTime;
clock_t s_time;
emili::LocalSearch* localsearch;

double emili::getCurrentExecutionTime()
{
    return (double)((clock()-beginTime)/ (double)CLOCKS_PER_SEC);
}

emili::Problem* emili::LocalSearch::theInstance = nullptr;

static void finalise (int _)
{
    keep_going = false;
    endTime = clock();
    std::cout << "CPU time: " << (endTime - beginTime) / (float)CLOCKS_PER_SEC << std::endl;
    emili::Solution* bestSoFar = localsearch->getBestSoFar();

    if(emili::LocalSearch::theInstance)
        emili::LocalSearch::theInstance->finaliseSolution(bestSoFar); // will apply "final" weight

    if(bestSoFar != nullptr)
    {
        double sol_val = bestSoFar->getSolutionValue();
        std::cout << "iteration counter : " << emili::iteration_counter()<< std::endl;
        std::cout << sol_val << std::endl;
        // std::cout << "Reached at time: " << (s_time - beginTime) / (float)CLOCKS_PER_SEC << std::endl;
        // std::cerr << (endTime - beginTime) / (float)CLOCKS_PER_SEC << " ";
        std::cerr << std::fixed << sol_val << std::endl; // to give hook-run the best so far solution
    }
    else
    {
        std::cout << "No valid solution found!" << std::endl;
        std::cerr << "timeout and no bestSoFar" << std::endl; // to inform hook-run
    }
#ifndef NOSIG
    _Exit(EXIT_SUCCESS);
#else
    exit(0);
#endif
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

emili::Solution::Value emili::Solution::getSolutionValue()
{    
    return solution_value;
}

std::string emili::Solution::getSolutionRepresentation()
{
    return "";
}

void emili::Solution::setSolutionValue(emili::Solution::Value value)
{
    this->solution_value = value;
}

void emili::Solution::swap(emili::Solution * other) {
    // very bad implementation, should swap data
    auto s = other->clone();
    *other = *this;
    *this = *s;
    delete s;
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
    n->reverseLastMove(line_);
    line_->setSolutionValue(this->base_->getSolutionValue());
    this->line_ = n->computeStep(this->line_);
    return *this;
}

emili::Solution* emili::Neighborhood::NeighborhoodIterator::operator *()
{
    return line_;
}

emili::Neighborhood::NeighborhoodIterator emili::Neighborhood::begin(emili::Solution *startSolution)
{
    // do not reset here because of FirstImprovement behaviour !
    if(needToResetWhenInstanceChanged) // unless problem (like kempe neigh)
        reset();
    return emili::Neighborhood::NeighborhoodIterator(this,startSolution);
}

emili::Neighborhood::NeighborhoodIterator emili::Neighborhood::end()
{
    return emili::Neighborhood::NeighborhoodIterator(this,nullptr);
}

void emili::Neighborhood::iterate(emili::Solution *base, std::function<void ()> yield) {
    auto it = begin(base);
    while(it != end()) {
        yield();
        ++it;
    }
}

/*
 * LocalSearch base class ( Old neighborhood concept)
 */

emili::Solution* emili::LocalSearch::search()
{
    if(init == nullptr)
        throw std::invalid_argument("local search was not constructed with initializer");
    neighbh->reset();
    emili::Solution* current = init->generateSolution();
    emili::Solution* sol = search(current);
    if(current!=sol)
        delete current;

    return sol;
}

emili::Solution* emili::LocalSearch::search(emili::Solution* initial)
{
    termcriterion->reset();
    neighbh->reset();
    bestSoFar = init->generateEmptySolution();

    emili::Solution* newSolution = init->generateEmptySolution();
    *newSolution = *initial;

    do
    {
        newSolution = neighbh->step(bestSoFar);
        if(*newSolution < *bestSoFar) {
            delete bestSoFar;
            bestSoFar = newSolution;
        } else {
            delete newSolution;
        }
    } while(!termcriterion->terminate(bestSoFar,newSolution)); // what happens when newSolution was just deleted ?

    return bestSoFar;
}

void emili::LocalSearch::searchInPlace(Solution *initial) {
    // default, dummy implementation
    auto sol = search(initial);
    if(sol != initial) {
        *initial = *sol;
        delete sol; // move semantics... (or swap)
    }
    bestSoFar = initial; // if the caller did keep the reference
}

emili::Solution* emili::LocalSearch::timedSearch(int time_seconds)
{
    setTimer(time_seconds);
    beginTime = clock();
    localsearch = this;
    emili::Solution* s = search();
    stopTimer();
    setBestSoFar(s);
    return s;
}

emili::Solution* emili::LocalSearch::timedSearch(int time_seconds, Solution *initial)
{
    setTimer(time_seconds);
    beginTime = clock();
    localsearch = this;
    emili::Solution* s = search(initial);
    stopTimer();
    setBestSoFar(s);
    return s;
}

emili::Solution* emili::LocalSearch::timedSearch()
{
    setTimer(seconds);
    beginTime = clock();
    localsearch = this;
    emili::Solution* s = search();
    stopTimer();
    setBestSoFar(s);
    return s;
}

emili::Solution* emili::LocalSearch::timedSearch(Solution *initial)
{
    setTimer(seconds);
    beginTime = clock();
    localsearch = this;
    emili::Solution* s = search(initial);
    stopTimer();
    setBestSoFar(s);
    return s;
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

emili::BestImprovementSearch::BestImprovementSearch(emili::Termination &terminationcriterion, emili::Neighborhood &neighborh)
    : emili::LocalSearch()
{
    this->init = nullptr;
    this->termcriterion = &terminationcriterion;
    this->neighbh = &neighborh;
    const_cast<Behaviour&>(behaviour) = Behaviour::FUNC;
}

emili::BestImprovementSearch::BestImprovementSearch(emili::InitialSolution &initialSolutionGenerator, emili::Termination &terminationcriterion, emili::Neighborhood &neighborh)
    : emili::LocalSearch(initialSolutionGenerator,terminationcriterion,neighborh)
{
    const_cast<Behaviour&>(behaviour) = Behaviour::FUNC;
}

/*
 * Best improvement local search
 */
emili::Solution* emili::BestImprovementSearch::search(emili::Solution* initial)
{
    termcriterion->reset();
    neighbh->reset();
    emili::Solution* incumbent = initial->clone();
    bestSoFar = init ? init->generateEmptySolution() : incumbent->clone();
    emili::Solution* ithSolution;
    Neighborhood::NeighborhoodIterator end = neighbh->end();
    do
    {
        *bestSoFar = *incumbent;
        Neighborhood::NeighborhoodIterator iter = neighbh->begin(bestSoFar);
        ithSolution = *iter;
        for(;iter!=end;++iter)
        {
            if(*ithSolution < *incumbent) {
                *incumbent = *ithSolution;
            }
        }
        delete ithSolution;
    } while(!termcriterion->terminate(bestSoFar, incumbent));
    delete incumbent;
    return bestSoFar;
}

/*
 * First improvement local search
 */

emili::FirstImprovementSearch::FirstImprovementSearch(emili::Termination &terminationcriterion, emili::Neighborhood &neighborh)
    : emili::LocalSearch()
{
    this->init = nullptr;
    this->termcriterion = &terminationcriterion;
    this->neighbh = &neighborh;
    const_cast<Behaviour&>(behaviour) = Behaviour::FUNC;
}


emili::FirstImprovementSearch::FirstImprovementSearch(emili::InitialSolution &initialSolutionGenerator, emili::Termination &terminationcriterion, emili::Neighborhood &neighborh)
    : emili::LocalSearch(initialSolutionGenerator,terminationcriterion,neighborh)
{
    const_cast<Behaviour&>(behaviour) = Behaviour::FUNC;
}

emili::Solution* emili::FirstImprovementSearch::search(emili::Solution* initial)
{
    termcriterion->reset();
    neighbh->reset();
    bestSoFar = init ? init->generateEmptySolution() : initial->clone();
    emili::Solution* incumbent = initial->clone();
    emili::Solution* ithSolution;

    //bestSoFar->setSolutionValue(bestSoFar->getSolutionValue()+1);
    Neighborhood::NeighborhoodIterator end = neighbh->end();
    do {
        *bestSoFar = *incumbent;
        // std::cout << "Best = Incumbent = " << bestSoFar << ", " << bestSoFar->getSolutionValue() << std::endl;
        Neighborhood::NeighborhoodIterator iter = neighbh->begin(incumbent); // incumbent may be different from initial ? This causes problem in Kempe Neighborhood
        ithSolution = *iter;
        // std::cout << "  Ith = " << ithSolution << ", " << ithSolution->getSolutionValue() << std::endl;
        int nIter = 0;
        for(;iter!=end;++iter, ++nIter)
        {
            // std::cout << "  Ith after " << nIter << " = " << ithSolution << ", " << ithSolution->getSolutionValue() << " : " << ithSolution->getSolutionRepresentation() <<  std::endl;
            if(*ithSolution < *incumbent) {
                *incumbent = *ithSolution;
                break;
            }
        }
        delete ithSolution;
        // std::cout << "  Incumbent in the end after " << nIter << " = " << incumbent->getSolutionValue() << std::endl;
    }while(!termcriterion->terminate(bestSoFar,incumbent));
    delete incumbent;
    return bestSoFar;
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
    if(current != sol)
        delete current;

    return sol;
}


emili::Solution* emili::BestTabuSearch::search(emili::Solution *initial)
{
    termcriterion->reset();
    neighbh->reset();
    emili::Solution* incumbent = initial->clone();
    bestSoFar = initial->clone();
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
    return bestSoFar;

}


emili::Solution* emili::FirstTabuSearch::search(emili::Solution *initial)
{
    termcriterion->reset();
    neighbh->reset();
    emili::Solution* incumbent = initial->clone();
    bestSoFar = initial->clone();
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
    return bestSoFar;
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

emili::Solution* emili::RandomMovePerturbationInPlace::perturb(Solution *solution) {
    for (int var = 0; var < numberOfSteps; ++var)
        explorer.randomStep(solution);
    return solution;
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

emili::Solution* emili::VNRandomMovePerturbationInPlace::perturb(Solution *solution)
{

    for (int var = 0; var < numberOfSteps; ++var)
        explorers[currentExplorer]->randomStep(solution);

    if(!(currentIteration <= numberOfIterations))
    {
        currentIteration = 0;
        currentExplorer = (currentExplorer + 1) % explorers.size();
    }
    else
    {
        currentIteration++;
    }

    return solution;
}

emili::Solution* emili::AlwaysAccept::accept(Solution *intensification_solution, Solution *diversification_solution)
{
    if(acc == ACC_DIVERSIFICATION)
        return diversification_solution;
    else
        return intensification_solution;
}

bool emili::AlwaysAccept::acceptViaDelta(Solution *diversification_solution, double delta)
{
    return acc == ACC_DIVERSIFICATION;
}

emili::Solution* emili::AcceptImproveEqual::accept(Solution *intensification_solution, Solution *diversification_solution)
{
    if(*diversification_solution <= *intensification_solution)
        return diversification_solution;
    else
        return intensification_solution;
}

bool emili::AcceptImproveEqual::acceptViaDelta(Solution *diversification_solution, double delta)
{
    return delta <= 0;
}

emili::Solution* emili::ImproveAccept::accept(Solution *intensification_solution, Solution *diversification_solution)
{
    if(*diversification_solution < *intensification_solution)
        return diversification_solution;
    else
        return intensification_solution;
}

bool emili::ImproveAccept::acceptViaDelta(Solution *newSolution, double delta)
{
    return delta < 0;
}

emili::Solution* emili::AcceptPlateau::accept(Solution *intensification_solution, Solution *diversification_solution)
{
    if(*diversification_solution <= *intensification_solution)
    {
        this->current_step = 0;
        this->threshold_status = 0;
        return diversification_solution;
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

bool emili::AcceptPlateau::acceptViaDelta(Solution *newSolution, double delta) {
    if(delta <= 0) {
        this->current_step = 0;
        this->threshold_status = 0;
        return true;
    } else {
        threshold_status++;
        if(threshold_status >= this->plateau_threshold)
        {
            if(current_step <= max_plateau_steps)
            {
                current_step++;
                return true;
            }
            else
            {
                threshold_status = 0;
                current_step = 0;
            }
        }
    }
    return false;
}


/*
 * Iterated Local Search
 */

emili::Solution* emili::IteratedLocalSearch::search(){
    termcriterion->reset();
    acc.reset();
    bestSoFar = ls.search();
    return search(bestSoFar);
}

emili::Solution *emili::MyIteratedLocalSearch::search() {
    termcriterion->reset();
    acc->reset();
    bestSoFar = ls->search();
    auto r = search(bestSoFar);
    assert(r == bestSoFar);
    return r;
}

emili::Solution* emili::IteratedLocalSearch::search(emili::Solution* initial){
    using std::cout;
    using std::endl;

    termcriterion->reset();
    acc.reset();
    bestSoFar = ls.search(initial);
    emili::Solution* s = bestSoFar->clone(); // what is the diff between clone() and x = init->generateEmptySolution(); *s = x;
    emili::Solution* s_s = nullptr;

    do {
        emili::Solution* s_p = pert.perturb(s);

        // local search on s_p
        if(s != s_s)
            delete s_s; // may be null

        s_s = ls.search(s_p);

        // best solution
        if(*s_s < *bestSoFar) {
            *bestSoFar = *s_s;
            //s_time = clock();
        }

        if(s != s_p)
            delete s_p;

        s_p = s;
        s = acc.accept(s_p, s_s);
        if(s != s_p)
            delete s_p;
    } while(! termcriterion->terminate(s,s_s));

    return bestSoFar;
}

emili::Solution *emili::MyIteratedLocalSearch::search(emili::Solution * initial) {
    termcriterion->reset();
    acc->reset();
    bestSoFar = ls->searchBehave(initial, Behaviour::FUNC);

    Solution* sol = bestSoFar->clone();
    for(;;) {
        Solution* newSol = per->perturbBehave(sol, Behaviour::FUNC);
        newSol = ls->searchBehave(newSol, Behaviour::VOID);

        if(*newSol < *bestSoFar)
            *bestSoFar = *newSol;

        bool done;
        if(acc->accept(sol, newSol) == newSol) {
            done = termcriterion->terminate(sol, newSol);
            delete sol;
            sol = newSol;
        } else {
            delete newSol;
            done = termcriterion->terminate(sol, sol);
        }
        if(done)
            break;
    }

    delete sol;

    return bestSoFar;
}

/*
emili::Solution* emili::IteratedLocalSearch::timedSearch(int maxTime)
{
        termcriterion->reset();
        acc.reset();
        localsearch = this;
        setTimer(maxTime);
        // search start
        beginTime = clock();
        bestSoFar = init->generateSolution();
        emili::Solution*  s = ls.search(bestSoFar);
        *bestSoFar = *s ;

        emili::Solution* s_s = nullptr;
        //initialization done
        do{
            //iteration_increment();
            //Perturbation step
            emili::Solution* s_p = pert.perturb(s);
           // std::cout << s_p->getSolutionValue() << std::endl;
            //local search on s_p
            if(s!=s_s && s_s != nullptr)
                delete s_s;

            s_s = ls.search(s_p);

            //best solution
            if(*s_s < *bestSoFar)
            {

                *bestSoFar = *s_s;
                //s_time = clock();
            //    std::cout << bestSoFar->getSolutionValue() << std::endl;
            }
            if(s != s_p)
                delete s_p;
            //acceptance step
            s_p = s;
            s = acc.accept(s_p,s_s);
            if(s != s_p)
                delete s_p;            
            //std::cout << "accepted fitness -> " << s->getSolutionValue() << std::endl;
            //end loop
        }while(!termcriterion->terminate(s,s_s) && keep_going);
        stopTimer();
        return bestSoFar;
}

emili::Solution* emili::IteratedLocalSearch::timedSearch(int maxTime,emili::Solution* initial)
{
        termcriterion->reset();
        acc.reset();
        setTimer(maxTime);
        localsearch = this;
        // search start
        beginTime = clock();
        bestSoFar = initial;
        bestSoFar = ls.search(initial);
        emili::Solution* s = bestSoFar;
        emili::Solution* s_s;
        //initialization done
        do{

            //Perturbation step
            emili::Solution* s_p = pert.perturb(s);
            //local search on s_p
            s_s = ls.search(s_p);
            delete s_p;
            //best solution
            if(*s_s < *bestSoFar)
            {
                bestSoFar = s_s;
                //s_time = clock();
            }
            //acceptance step
            s_p = s;
            s = acc.accept(s_p,s_s);
            if(s == s_p)
            {
                if(bestSoFar != s_s)
                    delete s_s;
            }
            else
            {
                if(bestSoFar != s_p)
                delete s_p;
            }
            //std::cout << "accepted fitness -> " << s->getSolutionValue() << std::endl;
            //end loop
        }while(!termcriterion->terminate(s,s_s) && keep_going);
        stopTimer();
        return bestSoFar;
}
*/

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
    return false;
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
bool emili::LocalMinimaTermination::terminate(Solution* currentSolution, Solution* newSolution)
{
    return newSolution == nullptr || currentSolution->operator <=(*newSolution);
}

/*
 * MaxSteps Termination
 */
bool emili::MaxStepsTermination::terminate(Solution *currentSolution, Solution *newSolution)
{
    if(current_step > max_steps_) {
        return true;
    }
    else {
        current_step++;
        return false;
    }
}

void emili::MaxStepsTermination::reset()
{
    current_step = 0;
}

bool emili::MaxStepsTerminationDebug::terminate(emili::Solution *currentSolution, emili::Solution *newSolution) {
    bool r = MaxStepsTermination::terminate(currentSolution, newSolution);

    double percent = 100 * current_step / max_steps_;
    while(percent >= currentPercent) {
        currentPercent += stepPercent;
        if(! prefix.empty())
            std::cout << prefix << ": ";
        std::cout << (int)(percent) << "%" << " " << currentSolution->getSolutionValue() << " vs " << newSolution->getSolutionValue() << std::endl;
    }

    return r;
}

void emili::MaxStepsTerminationDebug::reset()
{
    MaxStepsTermination::reset();
    currentPercent = 0;
}


/*
 * Piped Local Search
*/

emili::Solution* emili::PipeSearch::search(Solution *initial)
{
    Solution* bestSoFar = init->generateEmptySolution();
    *bestSoFar = *initial;
    Solution* ithSolution = bestSoFar;
    for(std::vector< emili::LocalSearch*>::iterator iter = lss.begin();iter!=lss.end();++iter)
    {
        ithSolution = (*iter)->search(ithSolution);
        if(*ithSolution < *bestSoFar)
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

bool emili::MetropolisAcceptance::acceptViaDelta(Solution* diversification_solution, double delta)
{
    // delta = divers - intens
    if(delta < 0) {
        float prob = std::exp(-delta/temperature);
        if(prob < 1.0 && generateRealRandomNumber() > prob)
            return false;
    }
    return true;
}

emili::Solution* emili::Metropolis::accept(Solution *intensification_solution, Solution *diversification_solution)
{
    if(counter == interval && temperature > end_temp)
    {     
        temperature = (alpha * temperature) - rate;
        counter = 0;
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

bool emili::Metropolis::acceptViaDelta(Solution *diversification_solution, double delta)
{
    if(counter == interval && temperature > end_temp)
    {
        temperature = (alpha * temperature) - rate;
        counter = 0;
    }
    counter++;
    if(delta < 0) {
        float prob = std::exp(-delta/temperature);
        if(prob < 1.0 && generateRealRandomNumber() > prob)
            return false;
    }

    return true;
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

