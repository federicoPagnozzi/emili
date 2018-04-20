//
//  Created by Federico Pagnozzi on 12/12/17.
//  Copyright (c) 2014 Federico Pagnozzi. All rights reserved.
//  This file is distributed under the BSD 2-Clause License. See LICENSE.TXT
//  for details.

#ifndef PROBLEMX_H
#define PROBLEMX_H
#include "../emilibase.h"

void notyet();



namespace emili{
namespace problemX{

class ProblemX: public emili::Problem{
public:     
    /**
     * @brief calcObjectiveFunctionValue
     * This function calculates the objective function value of solution and returns
     * it as a double.
     * @param solution
     * the solution to use in the calculation of the objective function value.
     * @return
     * a double representing the cost of solution
     */
    virtual double calcObjectiveFunctionValue(Solution& solution);
    /**
     * @brief evaluateSolution
     *  This function calculates the objective function value and, using setSolutionValue,
     *  update the solution value in the Solution class.
     * @param solution
     * @return
     * the return value should be the objective function value of the solution
     */
    virtual double evaluateSolution(Solution & solution);
    /**
     * @brief problemSize
     * @return
     *  if overloaded this function should return the problem size as an integer number
     *  This value is used by some Termination criterion and by the -ro running option
     */
    virtual int problemSize() {return 1;}
};

/**
 * @brief The Solution class
 * This class models a solution to an optimization problem.
 * I'm not sure if the fact that the solution must contain an instance
 * of the base class problem
 * is a good thing...
 */
class SolutionProblemX: public emili::Solution
{
protected:
    /**
     * @brief getRawData
     * It's ugly I know, but every problem has its own data structures.
     * The next version will have an object designed to be a data carrier.
     * @return
     *  a pointer to the raw data which is the problem dependent data structures of the solution
     */
    virtual const void* getRawData()const;
    /**
     * @brief setRawData
     * changes the rawdata of the solution to data.
     *  The definition of the corresponding instance variable is not here because this variable
     *  can be defined, in the child class, as a pointer to the actual type of raw data
     *  so you can avoid the use of a pointer to void.
     * @param data
     */
    virtual void setRawData(const void* data);
public:
    /**
     * @brief getSolutionRepresentation
     * @return
     * returns an empty string. A child class should return a representation of
     * the solution. ( e.g. for permutation flowshop a permutaion ).
     */
    virtual std::string getSolutionRepresentation();
    /**
     * @brief clone
     * this method should return a pointer to a clone of the Solution.
     * The clone is expected to be allocated on the heap (so created with 'new').
     * @return
     * return a pointer to a clone of the Solution.
     */
    virtual Solution* clone();
    /**
     *@brief isFeasible
     * This method should be overwritten if the problem has
     * unfeasible solutions.
     * this method should return true if the Solution is feasible and
     * false in the other case.
     */
    virtual bool isFeasible() {return true;}

    /**
     * @brief ~Solution
     * Destructor of Solution.
     * This method should always be overloaded and when called delete the raw data
     * (if allocated on the heap).
     */
    virtual ~SolutionProblemX() {}
};

/**
 * @brief The InitialSolution class
 * The initial solution generator
 */
class InitialSolutionProblemX: public emili::InitialSolution
{
public:
    /**
     * @brief InitialSolution
     *  The constructor only needs the problem instance to build an initial solution
     * @param problem_instance
     */
    InitialSolutionProblemX(Problem& problem_instance):emili::InitialSolution(problem_instance){}
    /**
     * @brief generateSolution
     * The generated solution must be a valid solution for the problem with
     * the appropriate data structures used for the implemented solution.
     * @return
     *  A new solution for instance
     */
    virtual Solution* generateSolution();
    /**
     * @brief generateEmptySolution
     * The method should generate an empty solution correctly allocated in
     * memory. This method is necessary because the LocalSearch class
     * and his extensions do not have any clou on how to allocate a solution.
     * @return
     * A solution object with all his data structures initialiazed
     */
    virtual Solution* generateEmptySolution();
    /**
     * @brief ~InitialSolution
     * This method should be overloaded by the child class to delete
     *  anything it allocates on the heap
     */
    virtual ~InitialSolutionProblemX() {}
};

/**
 * @brief The Neighborhood class
 *     The class models the neighborhood of a solution
 *     This class should return the neighbors of a base solution
 *     given a specific neighborhood relation.
 */
class NeighborhoodProblemX: public emili::Neighborhood
{
protected:
    /**
     * @brief computeStep
     * Takes step wich is a pointer to the base solution of the neighborhood
     * and applies the move.
     *
     * @return
     *  A solution representing the next neighbor
     */
    virtual Solution* computeStep(Solution* step);
    /**
     * @brief reverseLastMove
     * Takes a solution and undoes the last move
     * @param step
     * the solution to change back
     */
    virtual void reverseLastMove(Solution* step);
public:
       /**
        * @brief begin
        * Returns iterator pointing to the first neighbor
        * of base emili::Solution
        * @param base
        * base solution
        * @return
        * NeighborhoodIterator
        */
       virtual NeighborhoodIterator begin(emili::Solution* base);
       /** @brief reset
     * The state of the Neighborhood object may need to be restored to
     * initial conditions between local search calls
     * (e.g. first improvement strategies for permutation flow shop).
     */
    virtual void reset();
    /** @brief random
     * A method that returns a random solution in the neighborhood has to be provided.
     * This method shoudl return a new solution.
     */
    virtual Solution* random(Solution* currentSolution);
    /** @brief size
     * This method returns the size of the neighborhood
    */
    virtual int size(){return 1;}
    virtual ~NeighborhoodProblemX() {}
};

/** @brief The Perturbation class
* The pertubation phase of the ils.
*/
class PerturbationProblemX: public emili::Perturbation
{
  public:
    /**
     * @brief perturb
     * The method returns a perturbed solution generated starting from solution
     * @param solution
     * The solution to perturb
     * @return
     * perturbed solution
     */
    virtual Solution* perturb(Solution* solution);

    virtual ~PerturbationProblemX() { }
};
}
}
#endif // EMILIBASE_H
