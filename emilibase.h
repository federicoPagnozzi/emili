#ifndef EMILIBASE_H
#define EMILIBASE_H

/*

                                 ______ __  __ _____ _      _____
                                |  ____|  \/  |_   _| |    |_   _|
                                | |__  | \  / | | | | |      | |
                                |  __| | |\/| | | | | |      | |
                                | |____| |  | |_| |_| |____ _| |_
                                |______|_|  |_|_____|______|_____|


E M I L I stands for Easily Modifiable (or Moddable) Iterated Local search Implementation.

P.S.
It's a very bad acronym but it's the best I could came up with in 5 minutes...

In order to use this framework for running your things you should implement
at least a problem specific extension of the class Problem, Solution, InitialSolution,
Neighborhood, Pertubation and AcceptanceCriteria.
*/


#include <string>
#include <iostream>
#include <vector>
#ifndef NOC11
#include <random>
#else
//if the compiler does not support c++11 this compile path will be selected.
#include <tr1/random>
#define nullptr NULL
#endif
#include <functional>



namespace emili{

void iteration_counter_zero();
int iteration_counter();
void iteration_increment();

class Solution;
/*
 * The istance of the problem to solve.
*/
class Problem{
public:

    virtual double evaluateSolution(Solution & solution)=0;

};

/*
This class models a solution to an optimization problem.
I'm not sure if the fact that the solution must contain an instance
of the base class problem
is a good thing...
*/

class Solution
{
protected:    
    double solution_value;
public:
    Solution(double value):solution_value(value)    {    }
    virtual Solution& operator=(const Solution& a);
    virtual bool operator<(Solution& a);
    virtual bool operator<=(Solution& a);
    virtual bool operator>(Solution& a);
    /*
     * It's ugly I know, but every problem has its own data structures.
     * The next version will have an object designed to be a data carrier.
    */
    virtual const void* getRawData()const=0;
    virtual void setRawData(const void* data)=0;
    //virtual Solution clone()=0;
    virtual double getSolutionValue();
    virtual void setSolutionValue(double value);
    virtual ~Solution() {}
};

/*
    The initial solution generator
*/
class InitialSolution
{
protected:
    Problem& instance;
public:
    InitialSolution(Problem& problem_instance):instance(problem_instance){}
    /*
        The generated solution must be a valid solution for the problem with
        the appropriate data structures used for the implemented solution.
    */
    virtual Solution* generateSolution()=0;
    virtual Solution* generateEmptySolution()=0;
    virtual ~InitialSolution() {}
};

/*
    Extends this to implement your termination criteria
*/
class Termination
{
public:    
    /*
     * The name of the parameters are purely indicative
     * this method shall return true if the termination condition has been reached (otherwise false)
    */
    virtual bool terminate(Solution* currentSolution, Solution* newSolution)=0;
    /*
     *it should be possible to reset the state ( e.g. counters) of a Termination rule.
     */
    virtual void reset()=0;
};

/*
 * This termination criteria checks if newSolution improves on currentSolution,
 * if it does not improve the termination condition is verified and terminate returns true.
 */
class LocalMinimaTermination: public emili::Termination
{
public:
    LocalMinimaTermination() {}
    virtual bool terminate(Solution* currentSolution, Solution* newSolution);
    virtual void reset() { }
};
/*
 *  The equivalent of a while(true) in the search loop...
*/
class WhileTrueTermination: public Termination
{
public:
    virtual bool terminate(Solution* currentSolution, Solution* newSolution);
    virtual void reset() { }
};

/*
 * This starts a system timer whe reset it's called and
 * checks if the timer is expired every time terminate it's called.
 * This way there is not need anymore for the timedsearch method.
 */
class TimedTermination: public Termination
{
protected:
    int secs;
public:
    TimedTermination(int seconds):secs(seconds) { }
    TimedTermination():secs(1) { }
    virtual bool terminate(Solution *currentSolution, Solution *newSolution);
    virtual void reset();
};

class MaxStepsTermination : public Termination
{
protected:
    int max_steps_;
    int current_step;
public:
    MaxStepsTermination(int max_steps):max_steps_(max_steps), current_step(0){ }
    virtual bool terminate(Solution *currentSolution, Solution *newSolution);
    virtual void reset();
};

/*
    The class models a neighborhood of a solution

*/
class Neighborhood
{
protected:
    virtual Solution* computeStep(Solution* step)=0;
public:
    //unsigned long num();
       class NeighborhoodIterator : public std::iterator<std::forward_iterator_tag, emili::Solution> {
       public:
           NeighborhoodIterator(emili::Neighborhood* n,emili::Solution* startSolution):base_(startSolution),n(n)
           {
               if(startSolution != nullptr )
               {
                  line_ = n->computeStep(base_);
               }
               else
               {
                   line_=nullptr;
               }
           }
           NeighborhoodIterator& operator=(const NeighborhoodIterator& iter);
           bool operator==(const NeighborhoodIterator& iter);
           bool operator != (const NeighborhoodIterator& iter);
           NeighborhoodIterator& operator++();
           //NeighborhoodIterator& operator++(int);
           emili::Solution* operator*();
       private:
           emili::Solution* base_;
           emili::Solution* line_;
           emili::Neighborhood* n;
       };
       virtual NeighborhoodIterator begin(emili::Solution* base);
       virtual NeighborhoodIterator end();
    /*this method returns a solution in the decided neighborhood
     * of the currentSolution */
    virtual Solution* step(Solution* currentSolution)=0;
    /*
     * The state of the Neighborhood object may need to be restored to
     * initial conditions between local search calls
     * (e.g. first improvement strategies for permutation flow shop).
     */
    virtual void reset()=0;
    /*
     * A method that returns a random solution in the neighborhood has to be provided
     */
    virtual Solution* random(Solution* currentSolution) = 0;
    virtual ~Neighborhood() {}
};

/*
  This class models a very general local search.
*/
class LocalSearch
{
protected:
InitialSolution* init;
Termination* termcriteria;
Neighborhood* neighbh;
int seconds;
    LocalSearch() { }
public:
    LocalSearch(InitialSolution& initialSolutionGenerator ,Termination& terminationCriteria, Neighborhood& neighborh):
    init(&initialSolutionGenerator),termcriteria(&terminationCriteria),neighbh(&neighborh),seconds(0)    {    }

    LocalSearch(InitialSolution& initialSolutionGenerator ,Termination& terminationCriteria, Neighborhood& neighborh, int time):
    init(&initialSolutionGenerator),termcriteria(&terminationCriteria),neighbh(&neighborh),seconds(time)    {    }
    /*
     * search use the InitialSolutionGenerator instance
     * to generate the first solution for the local search
    */
    virtual Solution* search();
    /*
     * a starting solution can be also provided
     */
    virtual Solution* search(Solution* initial);
    /*
     * this method ends the execution of the algorithm when the termination criteria is true or
     * after the amount of seconds provided as argument (regardless of the value of the termination).
     */
    virtual Solution* timedSearch(int seconds);
    virtual Solution* timedSearch(int seconds, Solution* initial);
    /*
     * In order to make easier the creation of batchs of LocalSearch objects the time of execution
     * can be inserted as an instance variable in the constructor of the object so these methos below
     * rely on the value of that variable for the execution time.
     */
    virtual Solution* timedSearch();
    virtual Solution* timedSearch(Solution* initial);
    virtual void setSearchTime(int time);
    virtual int getSearchTime();
    emili::Termination& getTermination();
    emili::Neighborhood& getNeighborhood();
    emili::InitialSolution& getInitialSolution();
    virtual ~LocalSearch() { delete init; delete termcriteria; delete neighbh;}

};

/*
 this class models a best improvement local search using the dumb neighboor
 with the iterator interface.
 */
class BestImprovementSearch : public emili::LocalSearch
{
public:
    BestImprovementSearch(InitialSolution& initialSolutionGenerator ,Termination& terminationCriteria, Neighborhood& neighborh):emili::LocalSearch(initialSolutionGenerator,terminationCriteria,neighborh) {}
    virtual Solution* search(emili::Solution* initial);
    virtual Solution* search()
    {
        return emili::LocalSearch::search();
    }
    virtual Solution* timedSearch(int seconds, Solution* initial);
};

/*
 this class models a first improvement local search using the dumb neighboor
 with the iterator interface.
 */
class FirstImprovementSearch : public emili::LocalSearch
{
public:
    FirstImprovementSearch(InitialSolution& initialSolutionGenerator ,Termination& terminationCriteria, Neighborhood& neighborh):emili::LocalSearch(initialSolutionGenerator,terminationCriteria,neighborh) {}
    virtual Solution* search(emili::Solution* intial);
    virtual Solution* timedSearch(int seconds, Solution* initial);
};

/*
* The pertubation phase of the ils.
*/
class Perturbation
{
  public:
    virtual Solution* perturb(Solution* solution)=0;
};

class AcceptanceCriteria
{
public:
    /*
     *  the accept method decides the direction of the search by searching between intensification and diversification,
     *  the IteratedLocalSearch class calls this method putting as first paramenter the solution used for pertubation
     *  in the last iteration and as second parameter the result of the local search around the pertubed solution.
    */
    virtual Solution* accept(Solution* intensification_solution,Solution* diversification_solution)=0;
};
/*
 * IteratedLocalSearch it's a general implementation of iterated local search.
*/
class IteratedLocalSearch: public emili::LocalSearch
{
protected:
    LocalSearch& ls;
    Perturbation& pert;
    AcceptanceCriteria& acc;
public:
    IteratedLocalSearch(LocalSearch& localsearch,Termination& terminationCriteria,Perturbation& perturb,AcceptanceCriteria& accept):emili::LocalSearch(localsearch.getInitialSolution(),terminationCriteria,localsearch.getNeighborhood()),ls(localsearch),pert(perturb),acc(accept){}

    virtual Solution* search();
    virtual Solution* search(emili::Solution* initial);
    virtual Solution* timedSearch(int seconds);
    virtual Solution* timedSearch(int seconds,emili::Solution* initial);
};

/*
    This class models the memory of a tabu search.
*/
class TabuMemory
{
protected:
    int tabutenure;
public:
    TabuMemory(int tenureSize):tabutenure(tenureSize) { }
    TabuMemory():tabutenure(1) { }
    /*
     * tabu_check determines if the input it's a forbidden solution.
     * this method should return true if the solution is not tabu and false in the other case,
     */
    virtual bool tabu_check(emili::Solution* solution) = 0;
    /*
     * this method should mark the input as a forbidden solution.
     */
    virtual void forbid(emili::Solution* solution) = 0;
    /*
     * TabuSearch calls this method to let the tabumemory elaborate and store
     * a step (or move).
     * It should be overwritten only if the kind of memory
     * to be implemented needs to know this information.
     */
    virtual void registerMove(emili::Solution* base,emili::Solution* solution) {}
    /*
     * resets the state of the memory.
     * To be used to restore the object state between two distinct search.
     */
    virtual void reset()=0;
    virtual void setTabuTenure(int tt) { this->tabutenure = tt; }

};

class TabuNeighborhood: public emili::Neighborhood
{
protected:
    int tabutenure;
public:
    TabuNeighborhood(int tt_size):tabutenure(tt_size) { }
};

class TabuSearch: public emili::LocalSearch
{
protected:
    emili::TabuMemory& tabuMemory;
public:
    TabuSearch(InitialSolution& initialSolutionGenerator ,Termination& terminationCriteria, Neighborhood& neighborh,TabuMemory& tabut):
    emili::LocalSearch(initialSolutionGenerator,terminationCriteria,neighborh),tabuMemory(tabut) {    }
    virtual emili::Solution* search(emili::Solution* initial);
    virtual emili::Solution* search();
    virtual Solution* timedSearch(int seconds, Solution* initial);
    virtual Solution* timedSearch(int seconds);

};

/*
 * Variable Neighborhood Descent implementation
 * accept as template parameters a LocalSearch and uses it with the various kinds of neighborhoods
 */

template <class T>
class VNDSearch: public T
{
protected:
    std::vector < emili::Neighborhood* > neigh;
public:
    VNDSearch(emili::InitialSolution& is, emili::Termination& tc, std::vector< emili::Neighborhood* > n):T(is,tc,*n[0]),neigh(n) { }
    VNDSearch(T& ls, std::vector<emili::Neighborhood*> n):T(ls),neigh(n) { }
    virtual emili::Solution* search(emili::Solution *initial)
    {
        int i = 0;
        Solution* bestSoFar = this->init->generateEmptySolution();
        Solution* incumbent = this->init->generateEmptySolution();
        *incumbent  = *initial;
        do{
            this->neighbh = neigh[i];
            Solution* new_s = T::search(incumbent);
            if(new_s->operator <(*incumbent))
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
            }
        }while(i < neigh.size());
        return bestSoFar;
    }

};

class EmptyNeighBorHood: public emili::Neighborhood
{
protected:
    virtual Solution* computeStep(Solution *step);
public:
    virtual Solution* step(Solution *currentSolution);
    virtual void reset();
    virtual Solution* random(Solution *currentSolution);
};


class PipeSearch: public emili::LocalSearch
{
protected:
    std::vector <emili::LocalSearch*> lss;
public:
    PipeSearch(InitialSolution& is,std::vector< emili::LocalSearch*> lss):emili::LocalSearch(is,lss[0]->getTermination(),lss[0]->getNeighborhood()),lss(lss) { }
    virtual Solution* search(Solution* initial);
};

/*
   a consistent and centralized implementation of the random function
   that uses the marsenne twister.
 */


void initializeRandom(int seed);
#ifdef NOC11
std::tr1::mt19937& getRandomGenerator();
#else
std::mt19937& getRandomGenerator();
#endif
int generateRandomNumber();
float generateRealRandomNumber();

/*
 * Metropolis acceptance creteria implementation (fixed temperature)
 */
class MetropolisAcceptance: public emili::AcceptanceCriteria
{
protected:
    float temperature;
public:
    MetropolisAcceptance(float startTemp):temperature(startTemp) { }
    void setTemp(float temp);
    float getTemp();
    virtual Solution* accept(Solution* intensification_solution,Solution* diversification_solution);
};

/*
 *  This class models the pertubation algorithm that destruct the solution
 *  used in the iterated greedy algorithms
 */
class Destructor: public emili::Perturbation
{
public:
    virtual emili::Solution* destruct(Solution* solutioon)=0;
    virtual emili::Solution* perturb(Solution *solution) { return destruct(solution); }
};

/*
 * this class models the algorithm that returns a solution
 * starting from a partial one.
 */
class Constructor: public emili::LocalSearch
{
public:
  Constructor():emili::LocalSearch()
  {
      this->neighbh = nullptr;//new emili::EmptyNeighBorHood();
      this->init = nullptr;
      this->termcriteria = nullptr;
  }
 virtual emili::Solution* construct(emili::Solution* partial) = 0;
 virtual emili::Solution* constructFull() = 0;
 virtual emili::Solution* search() {return constructFull();}
 virtual emili::Solution* search(emili::Solution* initial) { return construct(initial);}
 virtual emili::Solution* timedSearch(int seconds, Solution *initial) { return construct(initial);}
};

/*
 * This class models an iterated greedy heuristics
 */
class IteratedGreedy : public emili::IteratedLocalSearch
{
public:
    IteratedGreedy(Constructor& c,Termination& t,Destructor& d,AcceptanceCriteria& ac):emili::IteratedLocalSearch(c,t,d,ac) { }
};

}
#endif // EMILIBASE_H
