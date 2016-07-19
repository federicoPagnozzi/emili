//
//  Created by Federico Pagnozzi on 28/11/14.
//  Copyright (c) 2014 Federico Pagnozzi. All rights reserved.
//  This file is distributed under the BSD 2-Clause License. See LICENSE.TXT
//  for details.

#ifndef EMILIBASE_H
#define EMILIBASE_H
/*

                                 ______ __  __ _____ _      _____
                                |  ____|  \/  |_   _| |    |_   _|
                                | |__  | \  / | | | | |      | |
                                |  __| | |\/| | | | | |      | |
                                | |____| |  | |_| |_| |____ _| |_
                                |______|_|  |_|_____|______|_____|


E.M.I.L.I. stands for Easily Modifiable (or Moddable) Iterated Local search Implementation.

P.S.
It's a very bad acronym but it's the best I could came up with in 5 minutes...

In order to use this framework for running your things you should implement
at least a problem specific extension of the class Problem, Solution, InitialSolution,
Neighborhood and Perturbation.
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
#include <stdexcept>


namespace emili{

void iteration_counter_zero();
int iteration_counter();
void iteration_increment();
void iteration_decrement();
class Solution;
/*
 * The instance of the problem to solve.
*/
class Problem{
public:

    virtual double evaluateSolution(Solution & solution)=0;
    virtual int problemSize() {return 1;}

};

/*
This class models a solution to an optimization problem.
I'm not sure if the fact that the solution must contain an instance
of the base class problem
is a good thing...
*/

class Solution
{
public:
    typedef double Value;
protected:    
    double solution_value;
    /*
     * It's ugly I know, but every problem has its own data structures.
     * The next version will have an object designed to be a data carrier.
    */
    virtual const void* getRawData()const=0;
    virtual void setRawData(const void* data)=0;    
public:
    Solution(double value):solution_value(value) {}
    virtual Solution& operator=(const Solution& a);    
    virtual bool operator<(Solution& a);
    virtual bool operator<=(Solution& a);
    virtual bool operator>=(Solution& a);
    virtual bool operator>(Solution& a);
    //virtual Solution clone()=0;
    virtual std::string getSolutionRepresentation();
    virtual Value getSolutionValue();
    virtual void setSolutionValue(Value value);
    virtual Solution* clone()=0;
    virtual void swap(Solution*);
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
    virtual Problem& getProblem() {return instance;}
    virtual ~InitialSolution() {}
};


/*
    Extends this to implement your termination criterion
*/
class Termination
{
public:    
    /*
     * The name of the parameters are purely indicative
     * this method shall return true if the termination condition has been reached (otherwise false)
    */
    virtual bool terminate(Solution* currentSolution, Solution* newSolution)=0;

    /**
     * Given a new solution and a "delta". Return true if the delta is accepted or not.
     * if ONLY_NEED_DELTA == true, this method is prefered to terminate
     */
    virtual bool terminateViaDelta(Solution *newSolution, double delta) {
        throw std::invalid_argument("not implemented");
    }

    /**
     * @brief ONLY_NEED_DELTA
     * if true, will call terminateViaDelta
     */
    const bool ONLY_NEED_DELTA = false;

    /*
     *it should be possible to reset the state ( e.g. counters) of a Termination rule.
     */
    virtual void reset()=0;

    virtual ~Termination() {}
};

/*
 * This termination criterion checks if newSolution improves on currentSolution,
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
 *
 * It checks if the timer is expired every time terminate it's called.
 *
 */
class TimedTermination: public Termination
{
protected:
    int secs;
    float _ratio;
    clock_t start;
public:
    TimedTermination(float ratio):secs(-1),_ratio(ratio) { }
    TimedTermination():secs(-1),_ratio(1) { }
    virtual bool terminate(Solution *currentSolution, Solution *newSolution);
    virtual void reset();
};
/*
 * As the name suggests after max_steps this Termination returns true.
 */
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
 * The class models a neighborhood of a solution
 * - A basic neighborhood is ITERABLE
 * - A random neighborhood can produce one RANDOM neighbor (SimulatedAnnealing)
 * - Note that "random exclusive" neighborhhods are not iterable (p probability insert vs 1-p probability swap for instance)
 *
 * reset() is called before begin() calls
 * Ways of iterating (depending on implementation)
 *
 * 1) Iterator
 * for(auto it = neigh.base(sol); it != neigh.end(); ++it)
 *     Solution* a = *it; // may need cast to your problem
 *     ...
 *
 * for(Solution* a : it)
 *     ...
 *
 * 2) iterate function (only for in place iterations)
 * TSPSolution* sol;
 * neigh.iterate(sol, [sol]{
 *     sol->graph...
 * });
*/
class Neighborhood
{
protected:
    /*
     * Takes step wich is a pointer to the base solution of the neighborhood
     * and applies the move.
     */
    virtual Solution* computeStep(Solution* step)=0;

    /**
     * First call (in begin()) to computeStep will be computeFirstStep
     * Useful when the neighborhood wants to setup init variables
     *
     * defaults to computeStep
     */
    virtual Solution* computeFirstStep(Solution* step) {
        return computeStep(step);
    }

    /*
     * Takes a solution and undoes the last move
     * Asserts computeStep or computeFirstStep was last non const method called
     * Asserts solution has not been changed since
     */
    virtual void reverseLastMove(Solution* step)=0;
public:
    //unsigned long num();

    class NeighborhoodIterator : public std::iterator<std::forward_iterator_tag, emili::Solution> {
    public:
        NeighborhoodIterator(emili::Neighborhood* n,emili::Solution* startSolution) : base_(startSolution), n(n)
        {
            if(startSolution != nullptr){
                line_ = base_->clone();
                line_ = n->computeFirstStep(line_);
            } else {
                line_ = nullptr;
            }
        }

        NeighborhoodIterator& operator=(const NeighborhoodIterator& iter);
        bool operator ==(const NeighborhoodIterator& iter);
        bool operator !=(const NeighborhoodIterator& iter);
        NeighborhoodIterator& operator++(); // only ++x semantics (iterator not copiable)
        emili::Solution* operator*();
    private:
        emili::Solution* base_;
        emili::Solution* line_;
        emili::Neighborhood* n;
    };

    const bool needToResetWhenInstanceChanged = false;

    /**
     * Used in iterator based iteration :
     * for(auto it = begin(base); it != end(); ++it)
     *    yield(*it);
     */
    virtual NeighborhoodIterator begin(emili::Solution* base);

    /**
     * Used in iterator based iteration
     */
    virtual NeighborhoodIterator end();

    /**
     * @brief see stdIterate
     * will reset the neighborhod in the beginning
     */
    struct StdReadyToIterate {
        Neighborhood* self;
        emili::Solution* sol;

        NeighborhoodIterator begin() {
            self->reset();
            return self->begin(sol);
        }

        NeighborhoodIterator end() {
            return self->end();
        }
    };

    /**
     * the neighborhood it reset in the beginning.
     * for(Solution* s : neigh->stdIterate(base))
     *      yield(s)
     */
    StdReadyToIterate stdIterate(emili::Solution* base) { return {this, base}; }


    /** this method returns a solution in the decided neighborhood
     * of the currentSolution */
    virtual Solution* step(Solution* currentSolution)=0;

    /*
     * The state of the Neighborhood object may need to be restored to
     * initial conditions between local search calls
     * (e.g. first improvement strategies for permutation flow shop).
     *
     * After reset, the state should the same as the one on creation.
     * In the constructor, reset() should be called
     * neighborhood.begin() will call reset()
     */
    virtual void reset()=0;

   /**
     * Alternative to computeStep and reverseLastStep
     * Iterate the neighborhood by modifying the state of <base>
     * yield will be called when the solution is been switched in a neighbor
     * yield may throw an exception if the iteration wants to stop
     *
     * generally an iterate function is way more easy to code than computeStep/reverseLastStep
     * reset() is not called
     */
    virtual void iterate(Solution* base, std::function<void()> yield);

    /*
     * @brief returns a new random solution in the neighborhood has to be provided
     * Used in simulated annealing.
     * computeStep in more preferred because it's generally fast and does not copy.
     */
    virtual Solution* random(Solution* currentSolution) = 0;

    /* // I suggest this default implementation
    {
        Solution* sol = currentSolution->clone();
        randomStep(sol);
        return sol;
    }
    */

    /**
     * @brief Do a random step on the solution.
     * May be called multiple times on multiple solutions.
     */
    virtual void randomStep(Solution* currentSolution) {
        throw std::invalid_argument("not implemented");
    }

    /**
     * @brief reverse last random step.
     * Asserts the last non const call was randomStep()
     * Asserts the solution has not been changed since
     */
    virtual void reverseLastRandomStep(Solution *currentSolution) {
        throw std::invalid_argument("not implemented");
    }

    /*
     * This method returns the size of the neighborhood
    */
    virtual int size()=0;
    virtual ~Neighborhood() {}
};

typedef Neighborhood::NeighborhoodIterator NeighborhoodIterator;

class EmptyNeighBorHood: public emili::Neighborhood
{
protected:
    virtual Solution* computeStep(Solution *step) {return nullptr;}
    virtual void reverseLastMove(Solution *step) { }
public:
    EmptyNeighBorHood() {}
    virtual Solution* step(Solution *currentSolution) {return currentSolution;}
    virtual void reset() { }
    virtual Solution* random(Solution *currentSolution) { return currentSolution;}
    virtual int size() { return 0;}
};

/**
 * @brief Behaviour of function y = f(x)
 * Imagine sol is a list
 *
 * FUNC | RETURN_NEW_NOT_MODIFY
 *    def f(x):
 *       y = list(x);
 *       y[0] = 5;
 *       return y
 * assert x unchanged
 * assert y is not x
 * assert x and y not deleted
 *
 * VOID | MODIFY_AND_RETURN_SAME
 *    def f(x):
 *       x[0] = 5
 *       return x
 *
 * assert x is y or x is deleted
 *
 * CAUTION, when calling applyWithBehaviour(... VOID), use :
 * x = search(x, VOID)
 * and NOT : search(x, VOID)
 * @see applyWithBehaviour
 *
 * MIX | RETURN_NEW_AND_MODIFY
 *    def f(x):
 *       x[0] = 5
 *       return list(x)
 *
 * assert x unchanged or x changed
 * assert y is not x
 * assert x not deleted and y not deleted
 *
 * UNKNOWN
 * assert x unchanged or x changed
 * assert y is x or y is not x
 * This shouldn't happen. Implementation must know if they return a new solution or not.
 */
enum struct Behaviour {
    FUNC,
    VOID,
    MIX
};

/**
 * Make clone and delete depending on needs
 *
 * y = f(x)
 * if target is VOID
 *    either x == y and x is modified
 *    either x != y and x is deleted
 * => when calling VOID use x = search(x, VOID)
 */
template <typename T>
Solution* applyWithBehaviour(Behaviour me, Behaviour target, Solution* initial, T search) {
    if(me == target)
        return search(initial);

    if(target == Behaviour::VOID) {
        Solution* y = search(initial);
        if(y != initial) // Normal func or Normal mix should have y != initial. In case it is lying about the "new", we don't do anything.
            delete initial;
        return y; // x = search(x, VOID)
    } else if(target == Behaviour::MIX) {
        if(me == Behaviour::VOID) {
            search(initial); // only semantics y is x
            return initial->clone();
        } else { // FUNC
            return search(initial); // lost semantics "NOT_MODIFY"
        }
    } else { // if(target == Behaviour::FUNC) {
        if(me == Behaviour::VOID) {
            Solution* x = initial->clone();
            search(x);
            return x;
        } else { // MIX
            Solution* x = initial->clone();
            Solution* y = search(x);
            delete x;
            return y;
        }
    }
}

/*
  This class models a very general local search.

  - Either that can start from nothing
    - see Solution* search())
    - init != nullptr
  - Either than can improve a solution
    - Solution* search(Solution*)
    - void searchInPlace(Solution*)
    - init == nullptr || init != nullptr
*/
class LocalSearch
{
protected:
    InitialSolution* init = nullptr;
    Termination* termcriterion = nullptr;
    Neighborhood* neighbh = nullptr;

    Solution* bestSoFar = nullptr;
    int seconds = 0;

    LocalSearch() { }

public:
    LocalSearch(InitialSolution& initialSolutionGenerator, Termination& terminationcriterion, Neighborhood& neighborh):
        init(&initialSolutionGenerator), termcriterion(&terminationcriterion), neighbh(&neighborh) {}

    LocalSearch(InitialSolution* initialSolutionGenerator, Termination* terminationcriterion, Neighborhood* neighborh):
        init(initialSolutionGenerator), termcriterion(terminationcriterion), neighbh(neighborh) {}

    LocalSearch(InitialSolution& initialSolutionGenerator ,Termination& terminationcriterion, Neighborhood& neighborh, int time):
        init(&initialSolutionGenerator),termcriterion(&terminationcriterion),neighbh(&neighborh), seconds(time) {}

    LocalSearch(InitialSolution* initialSolutionGenerator, Termination* terminationcriterion, Neighborhood* neighborh, int time):
        init(initialSolutionGenerator),termcriterion(terminationcriterion),neighbh(neighborh),seconds(time) {}

    /**
     * @brief Search use the InitialSolutionGenerator instance
     * to generate the first solution for the local search
     * @return new solution
     */
    virtual Solution* search();

    /**
     * @brief Search starting from a starting solution
     * @return new solution or initial modified (@see behaviour())
     */
    virtual Solution* search(Solution* initial);

    /**
     * @brief impose a behaviour
     */
    Solution* searchBehave(Solution* initial, Behaviour target) {
        return applyWithBehaviour(behaviour, target, initial, [this](Solution* sol){
            return search(sol);
        });
    }

    /**
     * @brief behaviour of search(Solution) -> Solution
     */
    const Behaviour behaviour = Behaviour::MIX; // MIX is the more general

    /**
     * @brief impose a VOID behavior
     * @deprecated
     */
    virtual void searchInPlace(Solution* initial);

    /*
     * this method ends the execution of the algorithm when the termination criterion is true or
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
    virtual Solution* getBestSoFar() { return bestSoFar;}
    virtual void setBestSoFar(Solution* newBest) {this->bestSoFar=newBest;}
    virtual ~LocalSearch() { delete init; delete termcriterion; delete neighbh;}

};


class EmptyLocalSearch: public emili::LocalSearch
{
public:
    EmptyLocalSearch(InitialSolution& in) : emili::LocalSearch() {
        this->init = &in;
        this->neighbh = new emili::EmptyNeighBorHood();
        this->termcriterion = new emili::MaxStepsTermination(0);
        const_cast<Behaviour&>(behaviour) = Behaviour::FUNC;
    }

    virtual Solution* search(Solution* initial) { return initial->clone(); }
    virtual Solution* timedSearch(int seconds, Solution *initial) { return initial->clone(); }
    virtual Solution* timedSearch(Solution* initial) { return initial->clone(); }
};

class IdentityLocalSearch: public emili::LocalSearch
{
public:
    IdentityLocalSearch():emili::LocalSearch() {
        this->init = nullptr;
        this->neighbh = new emili::EmptyNeighBorHood();
        this->termcriterion = new emili::MaxStepsTermination(0);
        const_cast<Behaviour&>(behaviour) = Behaviour::VOID;
    }

    IdentityLocalSearch(InitialSolution& in):emili::LocalSearch() {
        this->init = &in;
        this->neighbh = new emili::EmptyNeighBorHood();
        this->termcriterion = new emili::MaxStepsTermination(0);
        const_cast<Behaviour&>(behaviour) = Behaviour::VOID;
    }

    virtual Solution* search(Solution* initial) { return initial; }
    virtual Solution* timedSearch(int seconds, Solution *initial) { return initial; }
    virtual Solution* timedSearch(Solution* initial) {return initial; }
};

/*
 this class models a best improvement local search using the dumb neighborhood
 with the iterator interface.
 */
class BestImprovementSearch : public emili::LocalSearch
{
public:
    BestImprovementSearch(Termination& terminationcriterion, Neighborhood& neighborh);
    BestImprovementSearch(InitialSolution& initialSolutionGenerator, Termination& terminationcriterion, Neighborhood& neighborh);

    virtual Solution* search(emili::Solution* initial);

    virtual Solution* search()
    {
        return emili::LocalSearch::search();
    }
};

/*
 this class models a first improvement local search using the dumb neighborhood
 with the iterator interface.
 */
class FirstImprovementSearch : public emili::LocalSearch
{
public:
    FirstImprovementSearch(Termination& terminationcriterion, Neighborhood& neighborh);
    FirstImprovementSearch(InitialSolution& initialSolutionGenerator ,Termination& terminationcriterion, Neighborhood& neighborh);
    virtual Solution* search(emili::Solution* intial);
};

/*
* The pertubation phase of the ils.
*/
class Perturbation
{
  public:
    /**
     * @brief perturb a solution
     * @return new solution or modified solution
     */
    virtual Solution* perturb(Solution* solution)=0;
    virtual ~Perturbation() {}

    Solution* perturbBehave(Solution* solution, Behaviour target) {
        return applyWithBehaviour(behaviour, target, solution, [this](Solution* sol){
            return perturb(sol);
        });
    }

    Behaviour behaviour = Behaviour::MIX;
};

/*
 * NO pertubation
 * in place implementation
 */
class NoPerturbation: public emili::Perturbation
{
public:
    NoPerturbation() {
        behaviour = Behaviour::VOID;
    }

    virtual Solution* perturb(Solution *solution) { return solution; }
};

/*
    Performs a series of random steps in the given neighborhood.
*/
class RandomMovePerturbation : public emili::Perturbation
{
protected:
    Neighborhood& explorer;
    int numberOfSteps;
public:
    RandomMovePerturbation(Neighborhood& neighboorhod, int number_of_steps):explorer(neighboorhod),numberOfSteps(number_of_steps) {
        behaviour = Behaviour::FUNC;
    }

    virtual Solution* perturb(Solution* solution);
};

class RandomMovePerturbationInPlace : public RandomMovePerturbation
{
public:
    RandomMovePerturbationInPlace(Neighborhood& neighboorhod, int number_of_steps)
        : RandomMovePerturbation(neighboorhod, number_of_steps)
    {
        behaviour = Behaviour::VOID;
    }

    virtual Solution* perturb(Solution* solution);
};

class VNRandomMovePerturbation : public emili::Perturbation
{
protected:
  std::vector< Neighborhood* > explorers;
  int numberOfSteps;
  int numberOfIterations;
  int currentIteration;
  int currentExplorer;
public:
  VNRandomMovePerturbation(std::vector< Neighborhood* > neighborhoods, int number_of_steps, int number_of_iterations):explorers(neighborhoods),numberOfSteps(number_of_steps),numberOfIterations(number_of_iterations),currentIteration(0),currentExplorer(0) { }
  virtual Solution* perturb(Solution *solution);
};

class VNRandomMovePerturbationInPlace : public VNRandomMovePerturbation
{
public:
    VNRandomMovePerturbationInPlace(std::vector< Neighborhood* > neighborhoods, int number_of_steps, int number_of_iterations)
        : VNRandomMovePerturbation(neighborhoods, number_of_steps, number_of_iterations)
    {
        behaviour = Behaviour::VOID;
    }

    virtual Solution* perturb(Solution *solution);
};

class Acceptance
{
public:
    /*
     *  the accept method decides the direction of the search by choosing between intensification and diversification,
     *  the IteratedLocalSearch class calls this method putting as first paramenter the solution used for pertubation
     *  in the last iteration and as second parameter the result of the local search around the pertubed solution.
    */
    virtual Solution* accept(Solution* intensification_solution,Solution* diversification_solution)=0;

    /**
     * Given a new solution (diversification) and a "delta". Return true if the delta is accepted or not.
     * if ONLY_NEED_DELTA == true, this method is prefered to accept(Solution, Solution)
     */
    virtual bool acceptViaDelta(Solution *newSolution, double delta) {
        throw std::invalid_argument("not implemented");
    }

    /**
     * @brief ONLY_NEED_DELTA
     * if true, will call acceptViaDelta
     */
    const bool ONLY_NEED_DELTA = false;

    /**
     * @brief reset will be used if the algo is used repeatedly
     */
    virtual void reset() { }
    virtual ~Acceptance() {}
};

enum accept_candidates {ACC_INTENSIFICATION,ACC_DIVERSIFICATION};

class AlwaysAccept : public emili::Acceptance
{
protected:
    accept_candidates acc;
public:
    AlwaysAccept(accept_candidates choice):acc(choice) { }
    virtual Solution* accept(Solution *intensification_solution, Solution *diversification_solution);
};

class ImproveAccept : public emili::Acceptance
{
public:
    virtual Solution* accept(Solution *intensification_solution, Solution *diversification_solution);
};

class AcceptImproveEqual : public emili::Acceptance
{
public:
    virtual Solution* accept(Solution *intensification_solution, Solution *diversification_solution);
};

class AcceptPlateau : public emili::Acceptance
{
protected:
    int max_plateau_steps;
    int current_step;
    int threshold_status;
    int plateau_threshold;
public:
    AcceptPlateau(int maxNonImprovingSteps,int threshold):max_plateau_steps(maxNonImprovingSteps),plateau_threshold(threshold),current_step(0),threshold_status(0) { }
    virtual Solution* accept(Solution *intensification_solution, Solution *diversification_solution);
};

/*
 * IteratedLocalSearch it's a general implementation of iterated local search.
*/
class IteratedLocalSearch: public emili::LocalSearch
{
protected:
    LocalSearch& ls;
    Perturbation& pert;
    Acceptance& acc;

    // this.init must always be ls.init

public:
    IteratedLocalSearch(LocalSearch& localsearch,Termination& terminationcriterion,Perturbation& perturb,Acceptance& accept):emili::LocalSearch(localsearch.getInitialSolution(),terminationcriterion,localsearch.getNeighborhood()),ls(localsearch),pert(perturb),acc(accept){}

    virtual Solution* search();
    virtual Solution* search(emili::Solution* initial);
    virtual Solution* timedSearch(int seconds);
    virtual Solution* timedSearch(int seconds,emili::Solution* initial);
    virtual Solution* getBestSoFar();
};

class MyIteratedLocalSearch : public emili::LocalSearch
{
protected:
    LocalSearch* ls;
    Perturbation* per;
    Acceptance* acc;
public:
    MyIteratedLocalSearch(LocalSearch* ls, Termination* term, Perturbation* per, Acceptance* acc)
        : emili::LocalSearch()
    {
        this->termcriterion = term;
        this->ls = ls;
        this->per = per;
        this->acc = acc;

        const_cast<Behaviour&>(behaviour) = Behaviour::FUNC;
    }

    ~MyIteratedLocalSearch() {
        delete ls;
        delete per;
        delete acc;
    }

    Solution* search() override;
    Solution* search(Solution*) override;
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

class BestTabuSearch: public emili::LocalSearch
{
protected:
    emili::TabuMemory& tabuMemory;
public:
    BestTabuSearch(InitialSolution& initialSolutionGenerator ,Termination& terminationcriterion, Neighborhood& neighborh,TabuMemory& tabut):
    emili::LocalSearch(initialSolutionGenerator,terminationcriterion,neighborh),tabuMemory(tabut) {    }
    virtual emili::Solution* search(emili::Solution* initial);
    virtual emili::Solution* search();
};

class FirstTabuSearch: public emili::BestTabuSearch
{
public:
    FirstTabuSearch(InitialSolution& initialSolutionGenerator ,Termination& terminationcriterion, Neighborhood& neighborh,TabuMemory& tabut):
    emili::BestTabuSearch(initialSolutionGenerator,terminationcriterion,neighborh,tabut) {    }
    virtual emili::Solution* search(emili::Solution* initial);
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

        this->neighbh = neigh[0];
        Solution* incumbent = T::search(initial);
        int i = 0;
        do{
            this->neighbh = neigh[i];
            Solution* new_s = T::search(incumbent);
            if(*new_s < *incumbent)
            {
                delete incumbent;
                incumbent = new_s;
                i = 0;
            }
            else
            {
                i = i+1;
            }
        }while(i < neigh.size());
        return incumbent;
    }
};


class GVNS: public emili::LocalSearch
{
protected:
    std::vector < emili::Perturbation* > perturbations;
    emili::LocalSearch& ls;
public:
    GVNS(emili::LocalSearch& localsearch,std::vector< emili::Perturbation* > perturbs):emili::LocalSearch(),perturbations(perturbs),ls(localsearch) {this->init = &ls.getInitialSolution();this->neighbh = new emili::EmptyNeighBorHood(); }
    virtual Solution* search(emili::Solution* initial);
    virtual Solution* getBestSoFar();
    virtual ~GVNS() { delete neighbh;}
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
void initializeRandomFromRandom();
int getRandomSeedFromRandom();
#ifdef NOC11
std::tr1::mt19937& getRandomGenerator();
#else
std::mt19937& getRandomGenerator();
#endif

int generateRandomNumber();

/**
 * @return uniform variable between 0 and 1
 */
float generateRealRandomNumber();

/**
 * return x uniformly in 0 ... n-1
 */
int generateRandRange(int n);

/**
 * return x uniformly in from ... to-1
 */
int generateRandRange(int from, int to);

/**
 * return x uniformly in from ... to
 */
int generateRandInt(int from, int to);

/*TIME RELATED STUFF */
// this function returns the time from the beginning of the execution in seconds
double getCurrentExecutionTime();

/*
 * Metropolis acceptance criterion implementation (fixed temperature)
 */
class MetropolisAcceptance: public emili::Acceptance
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
 * Proper Metropolis acceptance criterion
 */
class Metropolis: public emili::Acceptance
{
protected:
    float temperature;
    float start_temp;
    float end_temp;
    int interval;
    int counter;
    float rate;
    float alpha;
public:
    Metropolis(float initial_temperature,float final_temperature,float descending_ratio):temperature(initial_temperature),start_temp(initial_temperature),end_temp(final_temperature),rate(descending_ratio),interval(1),counter(0),alpha(1) { }
    Metropolis(float initial_temperature,float final_temperature,float descending_ratio,int iterations):temperature(initial_temperature),start_temp(initial_temperature),end_temp(final_temperature),rate(descending_ratio),interval(iterations),counter(0),alpha(1) { }
    Metropolis(float initial_temperature, float final_temperature, float descending_ratio, int iterations, float alpha):temperature(initial_temperature),start_temp(initial_temperature),end_temp(final_temperature),rate(descending_ratio),interval(iterations),counter(0),alpha(alpha) { }
    virtual Solution* accept(Solution *intensification_solution, Solution *diversification_solution);
    virtual void reset();
};

/*
 *  This class models the pertubation algorithm that destruct the solution
 *  used in the iterated greedy algorithms
 */
class Destructor: public emili::Perturbation
{
public:
    virtual emili::Solution* destruct(Solution* solution)=0;
    virtual emili::Solution* perturb(Solution *solution) { return destruct(solution); }

    emili::Solution* destructBehave(Solution* initial, Behaviour target) {
        return applyWithBehaviour(behaviour, target, initial, [this](Solution* sol){
            return destruct(sol);
        });
    }
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
      this->neighbh = nullptr;//new emili::EmptyNeighborHood();
      this->init = nullptr;
      this->termcriterion = nullptr;
  }

  void setInitialSolution(InitialSolution* in) {
      init = in;
  }

  /**
   * @brief Construct from partial solution
   * @return new solution or partial modified (@see behaviour)
   */
  virtual emili::Solution* construct(emili::Solution* partial) = 0;

  /**
  * @brief Construct solution from nothing
  * @return new solution fully constructed (not partial)
  */
  virtual emili::Solution* constructFull() = 0;
  virtual emili::Solution* search() {return constructFull(); }
  virtual emili::Solution* search(emili::Solution* initial) { return construct(initial); }
  virtual emili::Solution* timedSearch(int seconds, Solution *initial) { return construct(initial); }

  emili::Solution* constructBehave(Solution* initial, Behaviour target) {
      return applyWithBehaviour(behaviour, target, initial, [this](Solution* sol){
          return construct(sol);
      });
  }
};

class ConstructDestructPertub : public emili::Perturbation {
    Constructor* cons;
    Destructor* des;
public:
    ConstructDestructPertub(Constructor* cons_, Destructor* des_)
        : cons(cons_), des(des_)
    {
        behaviour = Behaviour::VOID;
    }

    Solution* perturb(Solution* p) override {
        p = des->destructBehave(p, Behaviour::VOID);
        p = cons->constructBehave(p, Behaviour::VOID);
        return p;
    }
};

/*
 * This class models an iterated greedy heuristics
 */
class IteratedGreedy : public emili::IteratedLocalSearch
{
public:
    IteratedGreedy(Constructor& c,Termination& t,Destructor& d,Acceptance& ac):emili::IteratedLocalSearch(c,t,d,ac) { }
};

class InitAndSearch : public emili::LocalSearch
{
protected:
    LocalSearch* ls = nullptr;
public:
    // only search->search(sol) will be called, not search->search()
    InitAndSearch(InitialSolution* init, LocalSearch* ls) : LocalSearch() {
        this->init = init;
        this->ls = ls;
        const_cast<Behaviour&>(behaviour) = ls->behaviour;
    }

    Solution* search() override {
        return searchBehave(init->generateSolution(), Behaviour::VOID);
    }

    Solution* search(Solution* sol) override {
        return ls->search(sol);
    }
};

/**
 * @brief MyIteratedGreedy
 * search() uses cons->constructFull()
 * search(sol) uses ls(sol), destr(sol), cons(sol)
 */
class MyIteratedGreedy : public emili::LocalSearch {
    LocalSearch* ls;
    Destructor* dest;
    Constructor* cons;
    Acceptance* acc;
public:
    MyIteratedGreedy(Constructor* c, Termination* t, Destructor* d, Acceptance* ac, LocalSearch* ls)
        : emili::LocalSearch()
    {
        this->ls = ls;
        this->termcriterion = t;
        this->dest = d;
        this->cons = c;
        this->acc = ac;
        // init and neigh are forgotten

        const_cast<Behaviour&>(behaviour) = Behaviour::FUNC;
    }

    ~MyIteratedGreedy() {
        delete ls;
        delete dest;
        delete cons;
        delete acc;
    }

    /**
     * @return a new solution
     */
    Solution* search() override {
        return searchBehave(cons->constructFull(), Behaviour::VOID);
    }

    Solution* search(Solution* solInit) {
        bestSoFar = solInit->clone();

        // idea : if acc->ONLY_NEED_DELTA && term->ONLY_NEED_DELTA, then try VOID|VOID|VOID behaviour, doing nothing if accepted, using REVERT if not accepted
        // determines a way to quickly know if 3 REVERTs are better than 1 CLONE.

        bool terminate = false;
        Solution* sol = solInit->clone(); // FUNC behaviour
        while(! terminate) {
            // std::cout << "Starting from " << sol << ", " << sol->getSolutionValue() << " ";
            // std::cout.flush();

            Solution* newSol = dest->destructBehave(sol, Behaviour::FUNC);
            newSol = cons->searchBehave(newSol, Behaviour::VOID);
            newSol = ls->searchBehave(newSol, Behaviour::VOID);

            // std::cout << "Found:" << newSol->getSolutionValue() << std::endl;
            if(*newSol < *bestSoFar)
                *bestSoFar = *newSol;

            terminate = termcriterion->terminate(sol, newSol);

            if(acc->accept(sol, newSol) == newSol) {
                delete sol;
                sol = newSol;
            } else {
                delete newSol;
                // sol = sol;
            }
        }

        delete sol;

        return bestSoFar;
    }
};

class MyIteratedGreedyInit : public MyIteratedGreedy {
public:
    MyIteratedGreedyInit(InitialSolution* i, Constructor* c, Termination* t, Destructor* d, Acceptance* ac, LocalSearch* ls)
        : MyIteratedGreedy(c,t,d,ac,ls)
    {
        this->init = i;
    }

    Solution* search() override {
        return searchBehave(init->generateSolution(), Behaviour::VOID);
    }

};

/*
class SimulatedAnnealing : public emili::LocalSearch
{
protected:
emili::Acceptance* acceptance;
public:
SimulatedAnnealing(InitialSolution* initial,Neighborhood* neigh,Termination* term,Acceptance* acc):emili::LocalSearch(*initial,*term,*neigh),acceptance(acc) { }
virtual Solution* getBestSoFar();
virtual Solution* search(Solution *initial);
virtual ~SimulatedAnnealing() { delete acceptance;}
};
*/
}
#endif // EMILIBASE_H
