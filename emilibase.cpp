#include "emilibase.h"
#include <cstdlib>
#include <cstdio>
#include <signal.h>
#include <ctime>
#include <sys/time.h>
#include <iostream>
#include <assert.h>
#include "permutationflowshop.h"
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

bool keep_going;
bool timer_keep_going;
clock_t endTime;
clock_t beginTime;
clock_t s_time;
emili::Solution* s_cap;
struct itimerval timer;
struct itimerval termination_timer;

static void finalise (int _)
{

    keep_going = false;
    endTime = clock();
    std::cout << "CPU time: " << (endTime - beginTime) / (float)CLOCKS_PER_SEC << std::endl;
    if(s_cap)
    {
        cout << "iteration counter " << emili::iteration_counter()<< std::endl;
        std::cout << s_cap->getSolutionValue() << std::endl;

       //std::cout << "Reached at time: " << (s_time - beginTime) / (float)CLOCKS_PER_SEC << std::endl;
        //std::cerr << (endTime - beginTime) / (float)CLOCKS_PER_SEC << " ";
        std::cerr << s_cap->getSolutionValue() << std::endl;
        std::cerr << std::flush;
    }
    else
    {
        std::cout << "No valid solution found!" << std::endl;
    }
    std::cout << std::flush;
    _Exit(EXIT_SUCCESS);
}

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

/*
 * Timer HOOK
 */

/*
 * Iteration counter
 */
static int iteration_counter_ ;

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

bool emili::Solution::operator>(emili::Solution& a)
{
    return solution_value > a.solution_value;
}

double emili::Solution::getSolutionValue()
{
    return solution_value;
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
    this->line_ = n->computeStep(this->base_);
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
    delete current;
    return sol;
}

emili::Solution* emili::LocalSearch::timedSearch(int time_seconds)
{
    neighbh->reset();
    emili::Solution* current = init->generateSolution();
    emili::Solution* sol = timedSearch(time_seconds,current);
    delete current;
    return sol;
}

emili::Solution* emili::LocalSearch::search(emili::Solution* initial)
{
        termcriterion->reset();
        neighbh->reset();
        emili::Solution* current = init->generateEmptySolution();
        emili::Solution* newSolution = init->generateEmptySolution();

        *newSolution = *initial;

        do
        { 

            delete current;
            current = newSolution;
            newSolution = neighbh->step(current);

        }while(!termcriterion->terminate(current,newSolution));

        return current;
}

emili::Solution* emili::LocalSearch::timedSearch(int time_seconds, Solution *initial)
{
    termcriterion->reset();
    neighbh->reset();
    setTimer(time_seconds);
    beginTime = clock();
    emili::Solution* current = init->generateEmptySolution();
    emili::Solution* newSolution = init->generateEmptySolution();
    s_cap = initial;
    *newSolution = *initial;

    do
    {

        delete current;
        current = newSolution;
        newSolution = neighbh->step(current);
        if(newSolution->operator <(*s_cap)){
            s_cap = newSolution;
        }
    }while(!termcriterion->terminate(current,newSolution)&&keep_going&& isTimerUp());
    stopTimer();
    return current;
}

emili::Solution* emili::LocalSearch::timedSearch()
{
    neighbh->reset();
    emili::Solution* current = init->generateSolution();
    emili::Solution* sol = timedSearch(seconds,current);
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
        emili::Solution* bestSoFar = init->generateEmptySolution();
        emili::Solution* incumbent = init->generateEmptySolution();
        *incumbent = *initial;        
        emili::Solution* ithSolution = nullptr;

        do
        {

            if(bestSoFar != incumbent)
            {
            delete bestSoFar;
            }
            bestSoFar = incumbent;            
            neighbh->reset();
            Solution* bestOfTheIteration = incumbent;

            for(Neighborhood::NeighborhoodIterator iter = neighbh->begin(bestSoFar);iter!=neighbh->end();++iter)
            {
                ithSolution = *iter;

                if(bestOfTheIteration->operator >( *ithSolution)){
                    if(bestOfTheIteration!=incumbent)
                    delete bestOfTheIteration;

                    bestOfTheIteration = ithSolution;

                }
                else
                {
                    delete ithSolution;
                }

            }
            incumbent = bestOfTheIteration;

        }while(!termcriterion->terminate(bestSoFar,incumbent));

        return bestSoFar;
}

emili::Solution* emili::BestImprovementSearch::timedSearch(int seconds, Solution *initial)
{
    termcriterion->reset();
    setTimer(seconds);
    neighbh->reset();
    emili::Solution* bestSoFar = init->generateEmptySolution();
    emili::Solution* incumbent = init->generateEmptySolution();
    *incumbent = *initial;
    s_cap = initial;
    emili::Solution* ithSolution = nullptr;

    do
    {
        s_cap = incumbent;
        if(bestSoFar != incumbent)
        {
        delete bestSoFar;
        }
        bestSoFar = incumbent;
        neighbh->reset();
        Solution* bestOfTheIteration = incumbent;

        for(Neighborhood::NeighborhoodIterator iter = neighbh->begin(bestSoFar);iter!=neighbh->end();++iter)
        {
            ithSolution = *iter;

            if(bestOfTheIteration->operator >( *ithSolution)){
                if(bestOfTheIteration!=incumbent)
                delete bestOfTheIteration;

                bestOfTheIteration = ithSolution;

            }
            else
            {
                delete ithSolution;
            }

        }
        incumbent = bestOfTheIteration;

    }while(!termcriterion->terminate(bestSoFar,incumbent)&&keep_going&& isTimerUp());
    stopTimer();
    return bestSoFar;
}

/*
 * First improvement local search
 */

emili::Solution* emili::FirstImprovementSearch::search(emili::Solution* initial)
{
        termcriterion->reset();
        neighbh->reset();
        emili::Solution* bestSoFar = init->generateEmptySolution();
        emili::Solution* incumbent = bestSoFar;
        *incumbent = *initial;
        //bestSoFar->setSolutionValue(bestSoFar->getSolutionValue()+1);

        do{
            if(bestSoFar != incumbent)
            {
                delete bestSoFar;
            }
            bestSoFar = incumbent;
            for(Neighborhood::NeighborhoodIterator iter = neighbh->begin(incumbent);iter!=neighbh->end();++iter)
            {
                emili::Solution* ithSolution = *iter;
                if(incumbent->operator >(*ithSolution)){
                    if(incumbent!=bestSoFar)
                    delete incumbent;

                    incumbent = ithSolution;
                    break;
                }
                else
                {
                    delete ithSolution;
                }

            }

        }while(!termcriterion->terminate(bestSoFar,incumbent));
        return bestSoFar;
}

emili::Solution* emili::FirstImprovementSearch::timedSearch(int seconds, Solution *initial)
{
    termcriterion->reset();
    setTimer(seconds);
    neighbh->reset();
    emili::Solution* bestSoFar = init->generateEmptySolution();
    emili::Solution* incumbent = bestSoFar;
    *incumbent = *initial;
    s_cap = initial;
    bestSoFar->setSolutionValue(bestSoFar->getSolutionValue()+1);

    do{
        s_cap = incumbent;
        if(bestSoFar != incumbent)
        {
        delete bestSoFar;
        }
        bestSoFar = incumbent;
        for(Neighborhood::NeighborhoodIterator iter = neighbh->begin(incumbent);iter!=neighbh->end();++iter)
        {
            emili::Solution* ithSolution = *iter;
            if(incumbent->operator >(*ithSolution)){

                if(incumbent!=bestSoFar)
                    delete incumbent;

                incumbent = ithSolution;
                break;
            }
            else
            {
              delete ithSolution;
            }
        }



    }while(!termcriterion->terminate(bestSoFar,incumbent)&&keep_going&& isTimerUp());
    stopTimer();
    return bestSoFar;
}

/*
 * TABU SEARCH
 */
emili::Solution* emili::TabuSearch::search()
{
    tabuMemory.reset();
    emili::Solution* current = init->generateSolution();
    emili::Solution* sol = search(current);
    delete current;
    return sol;
}

emili::Solution* emili::TabuSearch::timedSearch(int seconds)
{
    tabuMemory.reset();
    emili::Solution* current = init->generateSolution();
    emili::Solution* sol = timedSearch(seconds,current);
    delete current;
    return sol;
}


emili::Solution* emili::TabuSearch::search(emili::Solution *initial)
{
    termcriterion->reset();
    neighbh->reset();    
    emili::Solution* incumbent = init->generateEmptySolution();
    *incumbent = *initial;
    emili::Solution* bestSoFar = incumbent;
    emili::Solution* ithSolution = nullptr;
    do
    {
        if(bestSoFar->operator >(*incumbent)){
            delete bestSoFar;            
            bestSoFar = incumbent;
        }

        neighbh->reset();
        Solution* bestOfTheIteration = incumbent;

        for(Neighborhood::NeighborhoodIterator iter = neighbh->begin(bestSoFar);iter!=neighbh->end();++iter)
        {
            ithSolution = *iter;            

            if(bestOfTheIteration->operator >( *ithSolution)){
                tabuMemory.registerMove(incumbent,ithSolution);

                if(tabuMemory.tabu_check(ithSolution))//<- Aspiration goes here.
                {
                    if(bestOfTheIteration!=incumbent)
                        delete bestOfTheIteration;

                    bestOfTheIteration = ithSolution;
                }
                else
                {
                    delete ithSolution;
                }

            }
            else
            {
                delete ithSolution;
            }


        }
        if(incumbent!=bestSoFar)
            delete incumbent;

        incumbent = bestOfTheIteration;
        tabuMemory.forbid(incumbent);
    }while(!termcriterion->terminate(bestSoFar,incumbent));
    return bestSoFar;
}

emili::Solution* emili::TabuSearch::timedSearch(int seconds, Solution *initial)
{
    termcriterion->reset();
    setTimer(seconds);
    neighbh->reset();
    emili::Solution* incumbent = init->generateEmptySolution();
    *incumbent = *initial;
    s_cap = incumbent;
    emili::Solution* ithSolution = nullptr;
    do
    {
        if(s_cap->operator >(*incumbent)){
            delete s_cap;
            s_cap = incumbent;
        }

        neighbh->reset();
        Solution* bestOfTheIteration = incumbent;

        for(Neighborhood::NeighborhoodIterator iter = neighbh->begin(s_cap);iter!=neighbh->end();++iter)
        {
            ithSolution = *iter;



            if(bestOfTheIteration->operator >( *ithSolution)){
                tabuMemory.registerMove(incumbent,ithSolution);

                if(tabuMemory.tabu_check(ithSolution))//<- Aspiration goes here.
                {
                    if(bestOfTheIteration!=incumbent)
                        delete bestOfTheIteration;

                    bestOfTheIteration = ithSolution;
                }
                else
                {
                    delete ithSolution;
                }

            }
            else
            {
                delete ithSolution;
            }



        }
        if(incumbent!=s_cap)
            delete incumbent;
        incumbent = bestOfTheIteration;
        tabuMemory.forbid(incumbent);
    }while(!termcriterion->terminate(s_cap,incumbent)&&keep_going&& isTimerUp());
    stopTimer();
    return s_cap;
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

emili::Solution* emili::RandomMovePertubation::perturb(Solution *solution)
{
    Solution* ret = explorer.random(solution);

    for (int var = 1; var < numberOfSteps; ++var) {
        Solution* temp = ret;
        ret = explorer.random(temp);
        delete temp;
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

emili::Solution* emili::ImproveAccept::accept(Solution *intensification_solution, Solution *diversification_solution)
{
    Solution* k = intensification_solution;
    if(intensification_solution->operator >(*diversification_solution)){
        k = diversification_solution;
    }
    return k;
}


/*
 * Iterated Local Search
 */

emili::Solution* emili::IteratedLocalSearch::search(){
    termcriterion->reset();
    emili::Solution* s_cap = ls.search();    
    emili::Solution* s = init->generateEmptySolution();
    s_cap = ls.search(s);
     *s = *s_cap;
    emili::Solution* s_s = nullptr;
    //initialization done
    do{

        //Pertubation step
        emili::Solution* s_p = pert.perturb(s);
        //local search on s_p
        if(s!=s_s && s_s != nullptr)
            delete s_s;
        s_s = ls.search(s_p);
        //best solution
        if(*s_s < *s_cap)
        {

            *s_cap = *s_s;
            //s_time = clock();
        }
        delete s_p;
        //acceptance step
        s_p = s;
        s = acc.accept(s_p,s_s);
        if(s != s_p)
            delete s_p;
        //std::cout << "accepted fitness -> " << s->getSolutionValue() << std::endl;
        //end loop
    }while(!termcriterion->terminate(s,s_s));
    return s_cap;
}

emili::Solution* emili::IteratedLocalSearch::search(emili::Solution* initial){
    termcriterion->reset();
    emili::Solution* scap = ls.search(initial);
    emili::Solution* s = init->generateEmptySolution();
     *s = *scap;
    emili::Solution* s_s = nullptr;
    //initialization done
    do{

        //Pertubation step
        emili::Solution* s_p = pert.perturb(s);
        //local search on s_p
        if(s!=s_s && s_s != nullptr)
            delete s_s;
        s_s = ls.search(s_p);
        //best solution
        if(*s_s < *scap)
        {

            *scap = *s_s;
            //s_time = clock();
        }
        delete s_p;
        //acceptance step
        s_p = s;
        s = acc.accept(s_p,s_s);
        if(s != s_p)
            delete s_p;
    }while(!termcriterion->terminate(s,s_s));
    return scap;
}

emili::Solution* emili::IteratedLocalSearch::timedSearch(int maxTime)
{
        termcriterion->reset();
        setTimer(maxTime);
        /*
            search start
        */
        beginTime = clock();
        s_cap = init->generateSolution();
        emili::Solution*  s = ls.search(s_cap);
        *s_cap = *s ;
        emili::Solution* s_s = nullptr;
        //initialization done
        do{

            //Pertubation step
            emili::Solution* s_p = pert.perturb(s);            
            //local search on s_p
            if(s!=s_s && s_s != nullptr)
                delete s_s;
            s_s = ls.search(s_p);
            //best solution
            if(*s_s < *s_cap)
            {

                *s_cap = *s_s;
                //s_time = clock();
            }
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
        return s_cap;
}

emili::Solution* emili::IteratedLocalSearch::timedSearch(int maxTime,emili::Solution* initial)
{
        termcriterion->reset();
        setTimer(maxTime);
        /*
            search start
        */
        beginTime = clock();
        s_cap = initial;
        s_cap = ls.search(initial);
        emili::Solution* s = s_cap;
        emili::Solution* s_s;
        //initialization done
        do{

            //Pertubation step
            emili::Solution* s_p = pert.perturb(s);
            //local search on s_p
            s_s = ls.search(s_p);
            delete s_p;
            //best solution
            if(*s_s < *s_cap)
            {
                s_cap = s_s;
                //s_time = clock();
            }
            //acceptance step
            s_p = s;
            s = acc.accept(s_p,s_s);
            if(s == s_p)
            {
                if(s_cap != s_s)
                    delete s_s;
            }
            else
            {
                if(s_cap != s_p)
                delete s_p;
            }
            //std::cout << "accepted fitness -> " << s->getSolutionValue() << std::endl;
            //end loop
        }while(!termcriterion->terminate(s,s_s) && keep_going);
        stopTimer();
        return s_cap;
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
    return !(time < secs);
}

void emili::TimedTermination::reset()
{

    //setTerminationTimer(this->secs);
    start = clock();
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
    if(current_step > max_steps_){
        return true;
    }
    else
    {
        current_step++;
        return false;
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
    if(temperature > end_temp)
    {
        temperature = temperature - rate;
    }
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

