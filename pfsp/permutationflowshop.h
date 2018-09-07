//
//  Created by Federico Pagnozzi on 28/11/14.
//  Copyright (c) 2014 Federico Pagnozzi. All rights reserved.
//  This file is distributed under the BSD 2-Clause License. See LICENSE.TXT
//  for details.
#ifndef PERMUTATIONFLOWSHOP_H
#define PERMUTATIONFLOWSHOP_H
#include "../emilibase.h"
#include "pfspinstance.h"
#include <iostream>
#include <cstdlib>
#include <limits>

/**
 *
 *  Permutation Flow shop components for EMILI
 *
*/
namespace emili
{
namespace pfsp
{
/** Permutation Flowshop problem implementation:
  The class uses code based on the code written by Jeremie ( see pfspinstance.h)
*/
class PermutationFlowShop: public emili::Problem
{
protected:
PfspInstance instance;
public:
    // Constructor that uses a PfspInstance implementation
    PermutationFlowShop(PfspInstance& problemInstance):instance(problemInstance) { }
    //Constructor that loads the instance from file path
    PermutationFlowShop(char* instance_path):instance()
    {
        /**  Read data from file */
        if (! instance.readDataFromFile(instance_path) ){
            exit(-1);
        }
    }
    /**
     computes the objective function value of solution.
     */
    virtual double calcObjectiveFunctionValue(Solution &solution);
    //implementation of evaluate solution
    virtual double evaluateSolution(emili::Solution& solution);
    /**  This method returns the number of jobs*/
    int getNjobs();
    /**  This method returns the number of machines*/
    int getNmachines();
    /**  This method returns the due date given the job*/
    int getDueDate(int job);
    /**  This method returns the priority given the job*/
    int getPriority(int job);
    /**  This method returns the due dates for all the jobs*/
    std::vector< long int >& getDueDates();
    /**  This method returns the priorities for all the jobs*/
    std::vector< long int >& getPriorities();
    /**  This method returns the processing time matrix so that it can be used by particular implementations of neighborhoods*/
    const std::vector< std::vector < long int > > & getProcessingTimesMatrix();
    /**  this method returns the problem size (used by some timed termination criteria)*/
    virtual int problemSize(){ return instance.getNbMac()*instance.getNbJob();}
    /**  This method returns the pfspinstance object that incapsulate the actual computation of the objective functions*/
    PfspInstance& getInstance();
    /** this method computes the makespan for the given solution.
    * The method is here because there are some heuristics that use it.
    */
    int computeMS(std::vector< int > & partial_solution);
    int computeMS(std::vector< int >& partial, int size);
    /**
     * The classes that extends this class to implement a PFSP objective
     * should implement these methods to calculate the objective functions
    */
    virtual int computeObjectiveFunction(std::vector< int > & partial_solution)=0;
    virtual int computeObjectiveFunction(std::vector< int > & partial_solution, int size)=0;
    virtual int computeObjectiveFunction(std::vector<int> &solution, std::vector<int>& makespans, int size);

    virtual long int computeObjectiveFunctionFromHead(std::vector<int> &solution, int starting_point, std::vector < std::vector < int > >& head,int njobs);
    virtual long int computeObjectiveFunctionFromHead(std::vector<int> &solution, int starting_point, std::vector < std::vector < int > >& head);
    /** This methods compute the matrices to implement Taillard's acceleration*/
    void computeTAmatrices(std::vector<int> &sol,std::vector< std::vector < int > >& head, std::vector< std::vector< int > >& tail);
    void computeTAmatrices(std::vector<int> &sol,std::vector< std::vector < int > >& head, std::vector< std::vector< int > >& tail,int size);
    void computeNoIdleTAmatrices(std::vector<int> &sol,std::vector< std::vector < int > >& head, std::vector< std::vector< int > >& tail);
    void computeHead(std::vector<int> &sol,std::vector< std::vector< int > >& head, int njobs) ;
    /** Old methods used by some particular and exceptional speed-ups*/
    int computeObjectiveFunction(std::vector<int> &sol,std::vector<int>& prevJob,int job,std::vector<int>& previousMachineEndTime);
    int computeObjectiveFunction(std::vector< int > & sol,std::vector<std::vector<int > >& previousMachineEndTimeMatrix, int start_i, int end_i);
    void computeWTs(std::vector<int> &sol,std::vector<int>& prevJob,int job,std::vector<int>& previousMachineEndTime);
    void computeTails(std::vector<int> &sol, std::vector< std::vector< std::vector< int > > > & tails);

    virtual ~PermutationFlowShop() { }

};
/**  CLASSIC PERMUTATION FLOW SHOP*/
/** Weighted Tardiness*/
class PFSP_WT: public PermutationFlowShop
{
public:
    PFSP_WT(PfspInstance& problemInstance):PermutationFlowShop(problemInstance) { }
    PFSP_WT(char* instance_path):PermutationFlowShop(instance_path) { }
    virtual int computeObjectiveFunction(std::vector< int > & partial_solution);
    virtual int computeObjectiveFunction(std::vector< int > & partial_solution, int size);
};
/** Weighted completion time*/
class PFSP_WCT: public PermutationFlowShop
{
public:
    PFSP_WCT(PfspInstance& problemInstance):PermutationFlowShop(problemInstance) { }
    PFSP_WCT(char* instance_path):PermutationFlowShop(instance_path) { }
    virtual int computeObjectiveFunction(std::vector< int > & partial_solution);
    virtual int computeObjectiveFunction(std::vector< int > & partial_solution, int size);
};
/** Total completion time*/
class PFSP_TCT: public PermutationFlowShop
{
public:
    PFSP_TCT(PfspInstance& problemInstance):PermutationFlowShop(problemInstance) { }
    PFSP_TCT(char* instance_path):PermutationFlowShop(instance_path) { }
    virtual int computeObjectiveFunction(std::vector< int > & partial_solution);
    virtual int computeObjectiveFunction(std::vector< int > & partial_solution, int size);
    virtual int computeObjectiveFunction(std::vector<int> &solution, std::vector<int>& makespans, int size);
};
/** Weighted Earliness*/
class PFSP_WE: public PermutationFlowShop
{
public:
    PFSP_WE(PfspInstance& problemInstance):PermutationFlowShop(problemInstance) { }
    PFSP_WE(char* instance_path):PermutationFlowShop(instance_path) { }
    virtual int computeObjectiveFunction(std::vector< int > & partial_solution);
    virtual int computeObjectiveFunction(std::vector< int > & partial_solution, int size);
};
/** Tardiness*/
class PFSP_T: public PermutationFlowShop
{
public:
    PFSP_T(PfspInstance& problem_instance):PermutationFlowShop(problem_instance) { }
    PFSP_T(char* instance_path):PermutationFlowShop(instance_path) { }
    virtual int computeObjectiveFunction(std::vector<int> &partial_solution);
    virtual int computeObjectiveFunction(std::vector<int> &partial_solution,int size);
    virtual int computeObjectiveFunction(std::vector<int> &solution, std::vector<int>& makespans, int size);
};
/** Earliness*/
class PFSP_E: public PermutationFlowShop
{
public:
    PFSP_E(PfspInstance& problem_instance):PermutationFlowShop(problem_instance) { }
    PFSP_E(char* instance_path):PermutationFlowShop(instance_path) { }
    virtual int computeObjectiveFunction(std::vector<int> &partial_solution);
    virtual int computeObjectiveFunction(std::vector<int> &partial_solution,int size);
};
/** Make span*/
class PFSP_MS: public PermutationFlowShop
{
public:
    PFSP_MS(PfspInstance& problem_instance):PermutationFlowShop(problem_instance) { }
    PFSP_MS(char* instance_path):PermutationFlowShop(instance_path) { }
    virtual int computeObjectiveFunction(std::vector<int> &partial_solution,int size);
    virtual int computeObjectiveFunction(std::vector<int> &partial_solution);
};


/**  NO WAIT PERMUTATION FLOW SHOP*/

/** Make span*/
class NWPFSP_MS: public PermutationFlowShop
{
protected:
    std::vector< std::vector < int > > distances;
    void computeNoWaitTimeDistances();
public:
    NWPFSP_MS(PfspInstance& problem_instance):
         PermutationFlowShop(problem_instance),
         distances(problem_instance.getProcessingTimesMatrix().size(),std::vector< int > (problem_instance.getProcessingTimesMatrix().size(),0))
         {
             computeNoWaitTimeDistances();
         }
    NWPFSP_MS(char* instance_path):PermutationFlowShop(instance_path) { }
    virtual int computeObjectiveFunction(std::vector<int> &partial_solution,int size);
    virtual int computeObjectiveFunction(std::vector<int> &partial_solution);
    const std::vector< std::vector < int > >& getDistances();
};
/** Weighted Tardiness*/
class NWPFSP_WT: public PermutationFlowShop
{
public:
    NWPFSP_WT(PfspInstance& problem_instance):PermutationFlowShop(problem_instance) { }
    NWPFSP_WT(char* instance_path):PermutationFlowShop(instance_path) { }
    virtual int computeObjectiveFunction(std::vector<int> &partial_solution,int size);
    virtual int computeObjectiveFunction(std::vector<int> &partial_solution);
};
/** Weighted Earliness*/
class NWPFSP_WE: public PermutationFlowShop
{
public:
    NWPFSP_WE(PfspInstance& problem_instance):PermutationFlowShop(problem_instance) { }
    NWPFSP_WE(char* instance_path):PermutationFlowShop(instance_path) { }
    virtual int computeObjectiveFunction(std::vector<int> &partial_solution,int size);
    virtual int computeObjectiveFunction(std::vector<int> &partial_solution);
};
/** Tardiness*/
class NWPFSP_T: public PermutationFlowShop
{
public:
    NWPFSP_T(PfspInstance& problem_instance):PermutationFlowShop(problem_instance) { }
    NWPFSP_T(char* instance_path):PermutationFlowShop(instance_path) { }
    virtual int computeObjectiveFunction(std::vector<int> &partial_solution,int size);
    virtual int computeObjectiveFunction(std::vector<int> &partial_solution);
};
/** Earliness*/
class NWPFSP_E: public PermutationFlowShop
{
public:
    NWPFSP_E(PfspInstance& problem_instance):PermutationFlowShop(problem_instance) { }
    NWPFSP_E(char* instance_path):PermutationFlowShop(instance_path) { }
    virtual int computeObjectiveFunction(std::vector<int> &partial_solution,int size);
    virtual int computeObjectiveFunction(std::vector<int> &partial_solution);
};
/** Weighted completion time*/
class NWPFSP_WCT: public PermutationFlowShop
{
public:
    NWPFSP_WCT(PfspInstance& problem_instance):PermutationFlowShop(problem_instance) { }
    NWPFSP_WCT(char* instance_path):PermutationFlowShop(instance_path) { }
    virtual int computeObjectiveFunction(std::vector<int> &partial_solution,int size);
    virtual int computeObjectiveFunction(std::vector<int> &partial_solution);
};
/** Total completion time*/
class NWPFSP_TCT: public PermutationFlowShop
{
public:
    NWPFSP_TCT(PfspInstance& problem_instance):PermutationFlowShop(problem_instance) { }
    NWPFSP_TCT(char* instance_path):PermutationFlowShop(instance_path) { }
    virtual int computeObjectiveFunction(std::vector<int> &partial_solution,int size);
    virtual int computeObjectiveFunction(std::vector<int> &partial_solution);
};

/**  NO IDLE PERMUTATION FLOW SHOP*/

/** Make span*/
class NIPFSP_MS: public PermutationFlowShop
{
public:
    NIPFSP_MS(PfspInstance& problem_instance):PermutationFlowShop(problem_instance) { }
    NIPFSP_MS(char* instance_path):PermutationFlowShop(instance_path) { }
    virtual int computeObjectiveFunction(std::vector<int> &partial_solution);
    virtual int computeObjectiveFunction(std::vector<int> &partial_solution,int size);
};
/** Make span with accelerations*/
class NI_A_PFSP_MS: public PermutationFlowShop
{
protected:
    long int nims_base;
    void calc_nims_base();
public:
    NI_A_PFSP_MS(PfspInstance& problem_instance):PermutationFlowShop(problem_instance),nims_base(0) { calc_nims_base();}
    NI_A_PFSP_MS(char* instance_path):PermutationFlowShop(instance_path),nims_base(0) { calc_nims_base();}
    virtual int computeObjectiveFunction(std::vector<int> &partial_solution);
    virtual int computeObjectiveFunction(std::vector<int> &partial_solution,int size);
};
/** Weighted Tardiness*/
class NIPFSP_WT: public PermutationFlowShop
{
public:
    NIPFSP_WT(PfspInstance& problem_instance):PermutationFlowShop(problem_instance) { }
    NIPFSP_WT(char* instance_path):PermutationFlowShop(instance_path) { }
    virtual int computeObjectiveFunction(std::vector<int> &partial_solution);
    virtual int computeObjectiveFunction(std::vector<int> &partial_solution,int size);
};
/** Weighted Earliness*/
class NIPFSP_WE: public PermutationFlowShop
{
public:
    NIPFSP_WE(PfspInstance& problem_instance):PermutationFlowShop(problem_instance) { }
    NIPFSP_WE(char* instance_path):PermutationFlowShop(instance_path) { }
    virtual int computeObjectiveFunction(std::vector<int> &partial_solution);
    virtual int computeObjectiveFunction(std::vector<int> &partial_solution,int size);
};
/** Tardiness*/
class NIPFSP_T: public PermutationFlowShop
{
public:
    NIPFSP_T(PfspInstance& problem_instance):PermutationFlowShop(problem_instance) { }
    NIPFSP_T(char* instance_path):PermutationFlowShop(instance_path) { }
    virtual int computeObjectiveFunction(std::vector<int> &partial_solution);
    virtual int computeObjectiveFunction(std::vector<int> &partial_solution,int size);
};
/** Earliness*/
class NIPFSP_E: public PermutationFlowShop
{
public:
    NIPFSP_E(PfspInstance& problem_instance):PermutationFlowShop(problem_instance) { }
    NIPFSP_E(char* instance_path):PermutationFlowShop(instance_path) { }
    virtual int computeObjectiveFunction(std::vector<int> &partial_solution);
    virtual int computeObjectiveFunction(std::vector<int> &partial_solution,int size);
};
/** Weighted completion time*/
class NIPFSP_WCT: public PermutationFlowShop
{
public:
    NIPFSP_WCT(PfspInstance& problem_instance):PermutationFlowShop(problem_instance) { }
    NIPFSP_WCT(char* instance_path):PermutationFlowShop(instance_path) { }
    virtual int computeObjectiveFunction(std::vector<int> &partial_solution);
    virtual int computeObjectiveFunction(std::vector<int> &partial_solution,int size);
};
/** Total completion time*/
class NIPFSP_TCT: public PermutationFlowShop
{
public:
    NIPFSP_TCT(PfspInstance& problem_instance):PermutationFlowShop(problem_instance) { }
    NIPFSP_TCT(char* instance_path):PermutationFlowShop(instance_path) { }
    virtual int computeObjectiveFunction(std::vector<int> &partial_solution);
    virtual int computeObjectiveFunction(std::vector<int> &partial_solution,int size);
};

/** Sequence dependent setup times Permutation flowshop*/

/** Make span*/
class SDSTFSP_MS: public PermutationFlowShop
{
public:
    SDSTFSP_MS(PfspInstance& problem_instance):PermutationFlowShop(problem_instance) { }
    SDSTFSP_MS(char* instance_path):PermutationFlowShop(instance_path) { }
    virtual int computeObjectiveFunction(std::vector<int> &partial_solution);
    virtual int computeObjectiveFunction(std::vector<int> &partial_solution,int size);
    virtual int computeObjectiveFunction(std::vector<int> &solution, std::vector<int>& makespans, int size);
    virtual long int computeObjectiveFunctionFromHead(std::vector<int> &solution, int starting_point, std::vector<std::vector<int> > &head, int njobs);
};
/** Weighted Tardiness*/
class SDSTFSP_WT: public SDSTFSP_MS
{
public:
    SDSTFSP_WT(PfspInstance& problem_instance):SDSTFSP_MS(problem_instance) { }
    SDSTFSP_WT(char* instance_path):SDSTFSP_MS(instance_path) { }
    virtual int computeObjectiveFunction(std::vector<int> &partial_solution);
    virtual int computeObjectiveFunction(std::vector<int> &partial_solution,int size);
    virtual int computeObjectiveFunction(std::vector<int> &solution, std::vector<int>& makespans, int size);
};
/** Weighted Earliness*/
class SDSTFSP_WE: public SDSTFSP_MS
{
public:
    SDSTFSP_WE(PfspInstance& problem_instance):SDSTFSP_MS(problem_instance) { }
    SDSTFSP_WE(char* instance_path):SDSTFSP_MS(instance_path) { }
    virtual int computeObjectiveFunction(std::vector<int> &partial_solution);
    virtual int computeObjectiveFunction(std::vector<int> &partial_solution,int size);
    virtual int computeObjectiveFunction(std::vector<int> &solution, std::vector<int>& makespans, int size);
};
/** Tardiness*/
class SDSTFSP_T: public SDSTFSP_MS
{
public:
    SDSTFSP_T(PfspInstance& problem_instance):SDSTFSP_MS(problem_instance) { }
    SDSTFSP_T(char* instance_path):SDSTFSP_MS(instance_path) { }
    virtual int computeObjectiveFunction(std::vector<int> &partial_solution);
    virtual int computeObjectiveFunction(std::vector<int> &partial_solution,int size);
    virtual int computeObjectiveFunction(std::vector<int> &solution, std::vector<int>& makespans, int size);
};
/** Earliness*/
class SDSTFSP_E: public SDSTFSP_MS
{
public:
    SDSTFSP_E(PfspInstance& problem_instance):SDSTFSP_MS(problem_instance) { }
    SDSTFSP_E(char* instance_path):SDSTFSP_MS(instance_path) { }
    virtual int computeObjectiveFunction(std::vector<int> &partial_solution);
    virtual int computeObjectiveFunction(std::vector<int> &partial_solution,int size);
    virtual int computeObjectiveFunction(std::vector<int> &solution, std::vector<int>& makespans, int size);
};
/** Total completion time*/
class SDSTFSP_TCT: public SDSTFSP_MS
{
public:
    SDSTFSP_TCT(PfspInstance& problem_instance):SDSTFSP_MS(problem_instance) { }
    SDSTFSP_TCT(char* instance_path):SDSTFSP_MS(instance_path) { }
    virtual int computeObjectiveFunction(std::vector<int> &partial_solution);
    virtual int computeObjectiveFunction(std::vector<int> &partial_solution,int size);
    virtual int computeObjectiveFunction(std::vector<int> &solution, std::vector<int>& makespans, int size);
};
/** Weighted completion time*/
class SDSTFSP_WCT: public SDSTFSP_MS
{
public:
    SDSTFSP_WCT(PfspInstance& problem_instance):SDSTFSP_MS(problem_instance) { }
    SDSTFSP_WCT(char* instance_path):SDSTFSP_MS(instance_path) { }
    virtual int computeObjectiveFunction(std::vector<int> &partial_solution);
    virtual int computeObjectiveFunction(std::vector<int> &partial_solution,int size);
    virtual int computeObjectiveFunction(std::vector<int> &solution, std::vector<int>& makespans, int size);
};


/** This class implements the Solution for the Permutation FlowShop problem
  It uses a vector of ints for storing the job sequence.
*/
class PermutationFlowShopSolution: public emili::Solution
{
protected:
    std::vector< int > solution;

    virtual const void* getRawData()const;
    virtual void setRawData(const void* data);
public:
    PermutationFlowShopSolution(double p_value):emili::Solution(p_value),solution()
    {}

    PermutationFlowShopSolution(std::vector< int >& solution):emili::Solution(1e9),solution(solution)
    {}

    PermutationFlowShopSolution(double p_value,std::vector< int >& solution):emili::Solution(p_value),solution(solution)
    {}
    /** Returns the job sequence that rapresents this solution*/
    virtual std::vector< int >& getJobSchedule();
    /** Set the job sequence*/
    virtual void setJobSchedule(std::vector<int>& newSeq);
    /** Returns a printable version of the job sequence*/
    virtual std::string getSolutionRepresentation();
    /** Implements the clone method of emili::Solution*/
    virtual emili::Solution* clone();
    /** Overrides the default == operator of Solution class*/
    virtual bool operator==(Solution& a);
    /** Destructor of the class*/
    virtual ~PermutationFlowShopSolution();
};

class PfspInitialSolution: public emili::InitialSolution
{
protected:
    PermutationFlowShop& pis;
    virtual Solution* generate() = 0 ;
public:
    PfspInitialSolution(PermutationFlowShop& problem_instance):emili::InitialSolution(problem_instance),pis(problem_instance) { }
    virtual Solution* generateSolution();
    /** This method generates a new empty solution by instantiating an empty vector of int of the correct size
    and setting the solution value to the biggest double number*/
    virtual Solution* generateEmptySolution();

};

class PfspRandomInitialSolution: public emili::pfsp::PfspInitialSolution
{
protected:    
    virtual Solution* generate();
public:
    PfspRandomInitialSolution(PermutationFlowShop& problem_instance):emili::pfsp::PfspInitialSolution(problem_instance){ }
};

class RandomInitialSolution: public emili::pfsp::PfspRandomInitialSolution
{
protected:
    virtual Solution* generate();
    int numOfSols;
public:
    RandomInitialSolution(PermutationFlowShop& problem_instance, int number_of_solutions):emili::pfsp::PfspRandomInitialSolution(problem_instance),numOfSols(number_of_solutions) { }
};

class PfspSlackInitialSolution: public emili::pfsp::PfspInitialSolution
{
protected:
    virtual Solution* generate();
public:
    PfspSlackInitialSolution(PermutationFlowShop& problem_instance):emili::pfsp::PfspInitialSolution(problem_instance){}    
};

class PfspNEHwslackInitialSolution: public emili::pfsp::PfspInitialSolution
{
protected:
    virtual Solution* generate();
public:
    PfspNEHwslackInitialSolution(PermutationFlowShop& problem_instance):emili::pfsp::PfspInitialSolution(problem_instance){}
};

class NEH: public emili::pfsp::PfspInitialSolution
{
protected:
    virtual Solution* generate();
public:
    NEH(PermutationFlowShop& problem_instance):emili::pfsp::PfspInitialSolution(problem_instance) {}
};

class NEHRS: public emili::pfsp::PfspInitialSolution
{
protected:
    int iterations;
    virtual Solution* generate();
public:
    NEHRS(PermutationFlowShop& problem_instance,int number_of_iterations):emili::pfsp::PfspInitialSolution(problem_instance),iterations(number_of_iterations) {}
};

class NEHedd: public emili::pfsp::NEH
{
protected:
    virtual Solution* generate();
public:
    NEHedd(PermutationFlowShop& problem_instance):emili::pfsp::NEH(problem_instance) { }
};

class NEHls: public NEH
{
protected:
    emili::LocalSearch* _ls;
    virtual Solution* generate();
public:
    NEHls(PermutationFlowShop& problem_instance,emili::LocalSearch* ls):emili::pfsp::NEH(problem_instance),_ls(ls) {}
};

class NEHeddLS: public NEHls
{
protected:
    virtual Solution* generate();
public:
    NEHeddLS(PermutationFlowShop& problem_instance,emili::LocalSearch* ls):emili::pfsp::NEHls(problem_instance, ls) {}
};

class NEHffls: public NEH
{
protected:
    emili::LocalSearch* _ls;
    virtual Solution* generate();
public:
    NEHffls(PermutationFlowShop& problem_instance,emili::LocalSearch* ls):emili::pfsp::NEH(problem_instance),_ls(ls) {}
};

class NEHff: public emili::pfsp::PfspInitialSolution
{
protected:
    virtual Solution* generate();
public:
    NEHff(PermutationFlowShop &problem_instance):emili::pfsp::PfspInitialSolution(problem_instance) {}
};

/** Less idle times construction heuristic from
        Wang CG, Chu CB, Proth JM. Heuristic approaches for n/m/F/SCi, scheduling
        problems. European Journal of Operational Research 1997;96(3):636–44.

 */
class LITSolution: public emili::pfsp::PfspInitialSolution
{
protected:
    virtual Solution* generate();
public:
    LITSolution(PermutationFlowShop& problem_instance):emili::pfsp::PfspInitialSolution(problem_instance) { }
};

class RZSolution: public emili::pfsp::PfspInitialSolution
{
protected:
    virtual Solution* generate();
public:
    RZSolution(PermutationFlowShop& problem):emili::pfsp::PfspInitialSolution(problem) { }
};

class MNEH: public emili::pfsp::PfspInitialSolution
{
protected:
    virtual Solution* generate();
public:
    MNEH(PermutationFlowShop& problem):emili::pfsp::PfspInitialSolution(problem) { }
};

class RMNEH: public emili::pfsp::PfspInitialSolution
{
protected:    
    virtual Solution* generate();
public:
    RMNEH(PermutationFlowShop& problem):emili::pfsp::PfspInitialSolution(problem) { }
};

class NeRZSolution: public emili::pfsp::PfspInitialSolution
{
protected:
    virtual Solution* generate();
public:
    NeRZSolution(PermutationFlowShop& problem):emili::pfsp::PfspInitialSolution(problem) { }
};

class NeRZ2Solution: public emili::pfsp::PfspInitialSolution
{
protected:
    virtual Solution* generate();
public:
    NeRZ2Solution(PermutationFlowShop& problem):emili::pfsp::PfspInitialSolution(problem) { }
};

class SRZSolution: public emili::pfsp::PfspInitialSolution
{
protected:
    virtual Solution* generate();
public:
    SRZSolution(PermutationFlowShop& problem):emili::pfsp::PfspInitialSolution(problem) { }
};

class NfRZ2Solution: public emili::pfsp::PfspInitialSolution
{
protected:
    virtual Solution* generate();
public:
    NfRZ2Solution(PermutationFlowShop& problem):emili::pfsp::PfspInitialSolution(problem) { }
};

class LRSolution: public emili::pfsp::PfspInitialSolution
{
protected:
    int number_of_sequences;
    virtual Solution* generate();
public:
    LRSolution(PermutationFlowShop& problem):emili::pfsp::PfspInitialSolution(problem),number_of_sequences(1) { }
    LRSolution(PermutationFlowShop& problem,int number_of_sequences):emili::pfsp::PfspInitialSolution(problem),number_of_sequences(number_of_sequences) { }
};

class NLRSolution: public emili::pfsp::LRSolution
{
protected:
    virtual Solution* generate();
public:
    NLRSolution(PermutationFlowShop& problem):emili::pfsp::LRSolution(problem) { }
    NLRSolution(PermutationFlowShop& problem,int number_of_sequences):emili::pfsp::LRSolution(problem,number_of_sequences) { }
};

class NRZPerturbation: public emili::Perturbation
{
protected:
    int d;
    emili::pfsp::PermutationFlowShop& prob;
public:
    NRZPerturbation(int d_parameter, emili::pfsp::PermutationFlowShop& problem):d(d_parameter),prob(problem) { }
    virtual emili::Solution* perturb(Solution *solution);
};

class TMIIGPerturbation: public emili::Perturbation
{
protected:
    int d;
    int tbsize;
    emili::pfsp::PermutationFlowShop& instance;
    std::vector< std::vector < int > > tblist;
public:
    TMIIGPerturbation(int d_parameter, emili::pfsp::PermutationFlowShop& problem,int tabu_list_size):d(d_parameter),instance(problem),tbsize(tabu_list_size),tblist(problem.getNjobs()+1 ) { }
    virtual emili::Solution* perturb(Solution *solution);
};

class NWTMIIGPerturbation: public emili::pfsp::TMIIGPerturbation
{
protected:
    const std::vector < std::vector < int > >& distances;
public:
    NWTMIIGPerturbation(int d_parameter, emili::pfsp::NWPFSP_MS& problem,int tabu_list_size):TMIIGPerturbation(d_parameter,problem,tabu_list_size),distances(problem.getDistances()){ }
    virtual emili::Solution* perturb(Solution *solution);
};

class IGPerturbation: public emili::Perturbation
{
protected:
    int d;
    emili::pfsp::PermutationFlowShop& instance;
public:
    IGPerturbation(int d_parameter, emili::pfsp::PermutationFlowShop& problem):d(d_parameter),instance(problem) { }
    virtual emili::Solution* perturb(Solution *solution);
};

class NWIGPerturbation: public emili::pfsp::IGPerturbation
{
protected:
    const std::vector < std::vector < int > >& distances;
public:
    NWIGPerturbation(int d_parameter, emili::pfsp::NWPFSP_MS& problem):emili::pfsp::IGPerturbation(d_parameter,problem),distances(problem.getDistances()) { }
    virtual emili::Solution* perturb(Solution *solution);
};

class IGOPerturbation: public emili::Perturbation
{
protected:
    int d;
    emili::pfsp::PermutationFlowShop& instance;
    std::vector < std::vector < int > > head;
    const std::vector< std::vector < long int > >& pmatrix;
    int nmac;
public:
    IGOPerturbation(int d_parameter, emili::pfsp::PermutationFlowShop& problem):d(d_parameter),instance(problem),pmatrix(problem.getProcessingTimesMatrix()),nmac(problem.getNmachines()),head(problem.getNmachines()+1,std::vector< int > (problem.getNjobs()+1,0)) { }
    virtual emili::Solution* perturb(Solution *solution);
};

class SDSTIGOPerturbation: public emili::pfsp::IGOPerturbation
{
public:
    SDSTIGOPerturbation(int dparameter, emili::pfsp::PermutationFlowShop& problem):
        IGOPerturbation(dparameter,problem) { }

    virtual emili::Solution* perturb(Solution *solution);
};

class IGIOPerturbation: public emili::Perturbation
{
protected:
    int d;
    std::vector< int > weights;
    emili::pfsp::PermutationFlowShop& instance;
    void updateWeights();
public:
    IGIOPerturbation(int d_parameter, emili::pfsp::PermutationFlowShop& problem):d(d_parameter),instance(problem),weights(problem.getNjobs()+1,0) { updateWeights();}
    virtual emili::Solution* perturb(Solution *solution);
};

class RSIOPerturbation: public emili::pfsp::IGIOPerturbation
{
protected:
    std::vector < std::vector < int > > head;
    std::vector < std::vector < int > > tail;
    const std::vector < std::vector < long int > >& pmatrix;
public:
    RSIOPerturbation(int d_parameter, emili::pfsp::PermutationFlowShop& problem):emili::pfsp::IGIOPerturbation(d_parameter,problem),head(problem.getNmachines()+1,std::vector< int > (problem.getNjobs()+1,0)),tail(problem.getNmachines()+1,std::vector< int >(problem.getNjobs()+1,0)),pmatrix(problem.getProcessingTimesMatrix()) {}
    virtual emili::Solution* perturb(Solution *solution);
};

class RSPerturbation: public emili::Perturbation
{
protected:
    int d;
    emili::pfsp::PermutationFlowShop& instance;
    std::vector < std::vector < int > > head;
    std::vector < std::vector < int > > tail;
    const std::vector < std::vector < long int > >& pmatrix;

public:
    RSPerturbation(int d_param, emili::pfsp::PermutationFlowShop& problem):d(d_param),instance(problem),head(problem.getNmachines()+1,std::vector< int > (problem.getNjobs()+1,0)),tail(problem.getNmachines()+1,std::vector< int >(problem.getNjobs()+1,0)),pmatrix(problem.getProcessingTimesMatrix()) { }
    virtual emili::Solution* perturb(Solution *solution);
};

class RSffPerturbation: public emili::Perturbation
{
protected:
    int d;
    emili::pfsp::PermutationFlowShop& instance;
    std::vector < std::vector < int > > head;
    std::vector < std::vector < int > > tail;
    const std::vector < std::vector < long int > >& pmatrix;

public:
    RSffPerturbation(int d_param, emili::pfsp::PermutationFlowShop& problem):d(d_param),instance(problem),head(problem.getNmachines()+1,std::vector< int > (problem.getNjobs()+1,0)),tail(problem.getNmachines()+1,std::vector< int >(problem.getNjobs()+1,0)),pmatrix(problem.getProcessingTimesMatrix()) { }
    virtual emili::Solution* perturb(Solution *solution);
};

class IgLsPerturbation: public emili::pfsp::IGPerturbation
{
protected:
    emili::LocalSearch* ls;
public:
    IgLsPerturbation(int d_parameter, emili::pfsp::PermutationFlowShop& problem, emili::LocalSearch* ls): emili::pfsp::IGPerturbation(d_parameter,problem),ls(ls) {/**    */}
    virtual emili::Solution* perturb(Solution *solution);
    ~IgLsPerturbation() { delete ls;}
};

class NwIgLsPerturbation: public emili::pfsp::IgLsPerturbation
{
protected:
    const std::vector< std::vector < int > >& distances;
public:
    NwIgLsPerturbation(int d_parameter, emili::pfsp::NWPFSP_MS& problem, emili::LocalSearch* ls): emili::pfsp::IgLsPerturbation(d_parameter,problem,ls),distances(problem.getDistances()) {/**    */}
    virtual emili::Solution* perturb(Solution *solution);
};

class IGOLsPerturbation: public emili::pfsp::IGOPerturbation
{
protected:
    emili::LocalSearch* ls;
public:
    IGOLsPerturbation(int d_parameter, emili::pfsp::PermutationFlowShop& problem, emili::LocalSearch* ls):emili::pfsp::IGOPerturbation(d_parameter, problem), ls(ls) { }
    virtual emili::Solution* perturb(Solution *solution);
    ~IGOLsPerturbation() { delete ls;}
};

class SDSTIGOLsPerturbation: public emili::pfsp::IGOLsPerturbation
{    
public:
    SDSTIGOLsPerturbation(int d_parameter, emili::pfsp::PermutationFlowShop& problem, emili::LocalSearch* ls):emili::pfsp::IGOLsPerturbation(d_parameter, problem,ls) { }
    virtual emili::Solution* perturb(Solution *solution);
};

class RSLSPerturbation: public emili::Perturbation
{
protected:
    int d;
    emili::pfsp::PermutationFlowShop& instance;
    std::vector < std::vector < int > > head;
    std::vector < std::vector < int > > tail;
    const std::vector < std::vector < long int > >& pmatrix;
    emili::LocalSearch* ls;
public:
    RSLSPerturbation(int d_param, emili::pfsp::PermutationFlowShop& problem, emili::LocalSearch* ls):d(d_param),instance(problem),head(problem.getNmachines()+1,std::vector< int > (problem.getNjobs()+1,0)),tail(problem.getNmachines()+1,std::vector< int >(problem.getNjobs()+1,0)),pmatrix(problem.getProcessingTimesMatrix()),ls(ls) { }
    virtual emili::Solution* perturb(Solution *solution);
    ~RSLSPerturbation() { delete ls;}
};

class RSffLSPerturbation: public emili::pfsp::RSLSPerturbation
{
public:
    RSffLSPerturbation(int d_param, emili::pfsp::PermutationFlowShop& problem, emili::LocalSearch* ls):emili::pfsp::RSLSPerturbation(d_param,problem,ls) { }
    virtual emili::Solution* perturb(Solution *solution);
};


class RestartPerturbation: public emili::Perturbation
{
protected:
    int num_of_solutions;
    bool locser;
    emili::InitialSolution* initial;
    emili::LocalSearch* ls;
public:
    RestartPerturbation(int np, emili::InitialSolution* init, emili::LocalSearch* ll):num_of_solutions(np),initial(init),ls(ll),locser(true) {}
    RestartPerturbation(int np, emili::InitialSolution* init):num_of_solutions(np),initial(init),ls(nullptr),locser(false) {}

    emili::Solution* perturb(Solution *solution);
    ~RestartPerturbation();
};

class MPTLMPerturbation: public emili::Perturbation
{
protected:
    const std::vector < std::vector < int > >& distances;
    int num_of_solutions;
    emili::pfsp::NWPFSP_MS& pis;
    int njobs;
public:
    MPTLMPerturbation(int np, emili::pfsp::NWPFSP_MS& init):num_of_solutions(np),pis(init),njobs(init.getNjobs()),distances(init.getDistances()) {}

    emili::Solution* perturb(Solution *solution);

};

class PfspNeighborhood: public emili::Neighborhood
{
protected:
    PermutationFlowShop& pis;    
    int njobs;
   virtual Solution* computeStep(Solution* step) =0;
public:
    PfspNeighborhood(PermutationFlowShop& problem):pis(problem),njobs(problem.getNjobs()){}
    virtual Solution* step(Solution* currentSolution);
    virtual PermutationFlowShop& getProblem() { return pis; }
    virtual void setNjobs(int num_of_jobs) {njobs=num_of_jobs;}
    virtual void reset();
    virtual std::pair<int,int> lastMove() { return std::pair<int,int>(0,0); }
    virtual Solution* random(Solution *currentSolution,int size) = 0;
    virtual int size();
};

/**
 * Basic insert neighborhood
 */
class PfspInsertNeighborhood: public emili::pfsp::PfspNeighborhood
{
protected:
    int start_position;
    int end_position;
    int sp_iterations;
    int ep_iterations;
    std::vector < int > current;
    int current_value;
    virtual Solution* computeStep(Solution* value);
    virtual void reverseLastMove(Solution *step);
public:
    PfspInsertNeighborhood(PermutationFlowShop& problem):PfspNeighborhood(problem),start_position(0),end_position(0),sp_iterations(1),ep_iterations(1){}
    virtual void reset();
    virtual Solution* random(Solution *currentSolution);
    virtual Solution* random(Solution *currentSolution,int size);
    virtual std::pair<int,int> lastMove() { return std::pair<int,int>(end_position,start_position); }
    virtual NeighborhoodIterator begin(Solution *base);    
};

/**
 * Insert neighborhood with Taillard's acceleration
 */
class TaillardAcceleratedInsertNeighborhood: public emili::pfsp::PfspInsertNeighborhood
{
protected:
    std::vector < std::vector < int > > head;
    std::vector < std::vector < int > > tail;
    const std::vector < std::vector < long int > >& pmatrix;
    const int nmac;
    void computeTAmatrices(std::vector<int>& sol);
    virtual Solution* computeStep(Solution *value);
public:
    TaillardAcceleratedInsertNeighborhood(PermutationFlowShop& problem):emili::pfsp::PfspInsertNeighborhood(problem),head(problem.getNmachines()+1,std::vector< int > (problem.getNjobs()+1,0)),tail(problem.getNmachines()+1,std::vector< int >(problem.getNjobs()+1,0)),pmatrix(problem.getProcessingTimesMatrix()),nmac(problem.getNmachines()) { }
    virtual NeighborhoodIterator begin(Solution *base);
};

/**
 * Insert neighborhood with Taillard's acceleration
 * that does a full scan each iteration
 */
class CSTaillardAcceleratedInsertNeighborhood: public emili::pfsp::TaillardAcceleratedInsertNeighborhood
{
protected:
    virtual Solution* computeStep(Solution *value);
public:
    CSTaillardAcceleratedInsertNeighborhood(PermutationFlowShop& problem):emili::pfsp::TaillardAcceleratedInsertNeighborhood(problem){ }
};

/**
 * Insert neighborhood
 * that does a full scan each iteration for plain Permutation Flowshop
 */
class CSInsertNeighborhood: public emili::pfsp::TaillardAcceleratedInsertNeighborhood
{
protected:
    virtual Solution* computeStep(Solution *value);
public:
    CSInsertNeighborhood(PermutationFlowShop& problem):emili::pfsp::TaillardAcceleratedInsertNeighborhood(problem){ }
};

/**
 * Insert neighborhood
 * that does a full scan each iteration for Flowshop with Sequence Dependent Setup Times
 */
class SDSTCSInsertNeighborhood: public emili::pfsp::CSInsertNeighborhood
{
protected:
    virtual Solution* computeStep(Solution *value);
public:
    SDSTCSInsertNeighborhood(PermutationFlowShop& problem):emili::pfsp::CSInsertNeighborhood(problem){ }
     virtual NeighborhoodIterator begin(Solution *base);
};

/**
 * Insert neighborhood with Taillard's acceleration
 * that changes the base solution after each improvement
 */
class FSTaillardAcceleratedInsertNeighborhood: public emili::pfsp::TaillardAcceleratedInsertNeighborhood
{
protected:
    bool improved;    
    virtual Solution* computeStep(Solution *value);
    virtual void reverseLastMove(Solution *step);
    int current_value;
public:
    FSTaillardAcceleratedInsertNeighborhood(PermutationFlowShop& problem):emili::pfsp::TaillardAcceleratedInsertNeighborhood(problem),improved(false),current_value(0){ }
    virtual NeighborhoodIterator begin(Solution *base);
};

/**
 * One level approximation no threshold for Weigthed Tardiness
 */
class HeavilyApproximatedTaillardAcceleratedInsertNeighborhood: public emili::pfsp::TaillardAcceleratedInsertNeighborhood
{
protected:
    virtual void computeHead(std::vector<int>& sol);
    virtual Solution* computeStep(Solution *value);
    std::vector< long int >& duedates;
    std::vector< long int >& priorities;
public:
    HeavilyApproximatedTaillardAcceleratedInsertNeighborhood(PermutationFlowShop& problem):emili::pfsp::TaillardAcceleratedInsertNeighborhood(problem),priorities(problem.getPriorities()),duedates(problem.getDueDates()) { }
    virtual NeighborhoodIterator begin(Solution *base);
};

/**
 * Kar2016 random neighborhood
 */

class KarNeighborhood: public emili::pfsp::HeavilyApproximatedTaillardAcceleratedInsertNeighborhood
{
protected:
    int lastMoveType;
    virtual Solution* computeStep(Solution *value);
    virtual void reverseLastMove(Solution *step);
public:
    KarNeighborhood(PermutationFlowShop& problem):HeavilyApproximatedTaillardAcceleratedInsertNeighborhood(problem),lastMoveType(0) { }
    virtual Solution* random(Solution *currentSolution);
    virtual NeighborhoodIterator begin(Solution *base);
};

/**  This Insert recomputes the objective function value only for the modified parts of the solution
 * for Weigthed Tardiness
 * */
class OptInsert: public emili::pfsp::TaillardAcceleratedInsertNeighborhood
{
protected:
    virtual Solution* computeStep(Solution *value);
public:
    OptInsert(PermutationFlowShop& problem):emili::pfsp::TaillardAcceleratedInsertNeighborhood(problem) { }
};

/**
 * One level approximation
 */
class NatxNeighborhood: public emili::pfsp::HeavilyApproximatedTaillardAcceleratedInsertNeighborhood
{
protected:    
    virtual Solution* computeStep(Solution *value);
public:
    NatxNeighborhood(PermutationFlowShop& problem):emili::pfsp::HeavilyApproximatedTaillardAcceleratedInsertNeighborhood(problem) { }
};
/**
 * One level approximation with experimental performance improvement tricks
 */
class Natx2Neighborhood: public emili::pfsp::HeavilyApproximatedTaillardAcceleratedInsertNeighborhood
{
protected:
    int thresh;
    int value_wt;
    virtual Solution* computeStep(Solution *value);
public:
    Natx2Neighborhood(PermutationFlowShop& problem):emili::pfsp::HeavilyApproximatedTaillardAcceleratedInsertNeighborhood(problem),thresh(problem.getNjobs()/2) { }
    Natx2Neighborhood(PermutationFlowShop& problem, int starting_threshold):emili::pfsp::HeavilyApproximatedTaillardAcceleratedInsertNeighborhood(problem),thresh(starting_threshold) { }
    virtual NeighborhoodIterator begin(Solution *base);
};
/**
 * zero level approximation
 **/
class AtxNeighborhood: public emili::pfsp::HeavilyApproximatedTaillardAcceleratedInsertNeighborhood
{
protected:
    virtual Solution* computeStep(Solution *value);
public:
    AtxNeighborhood(PermutationFlowShop& problem):emili::pfsp::HeavilyApproximatedTaillardAcceleratedInsertNeighborhood(problem) { }
};

/**
 * Two level approximation
 */
class EatxNeighborhood: public emili::pfsp::HeavilyApproximatedTaillardAcceleratedInsertNeighborhood
{
protected:
    virtual Solution* computeStep(Solution *value);
public:
    EatxNeighborhood(PermutationFlowShop& problem):emili::pfsp::HeavilyApproximatedTaillardAcceleratedInsertNeighborhood(problem) { }
};

/**
 * Three level approximation
 */
class ThatxNeighborhood: public emili::pfsp::HeavilyApproximatedTaillardAcceleratedInsertNeighborhood
{
protected:
    virtual Solution* computeStep(Solution *value);
public:
    ThatxNeighborhood(PermutationFlowShop& problem):emili::pfsp::HeavilyApproximatedTaillardAcceleratedInsertNeighborhood(problem) { }
};

/**
 * Four level approximation
 */
class FatxNeighborhood: public emili::pfsp::HeavilyApproximatedTaillardAcceleratedInsertNeighborhood
{
protected:
    virtual Solution* computeStep(Solution *value);
public:
    FatxNeighborhood(PermutationFlowShop& problem):emili::pfsp::HeavilyApproximatedTaillardAcceleratedInsertNeighborhood(problem) { }
};
/**
 * Five level approximation
 */
class PatxNeighborhood: public emili::pfsp::HeavilyApproximatedTaillardAcceleratedInsertNeighborhood
{
protected:
    virtual Solution* computeStep(Solution *value);
public:
    PatxNeighborhood(PermutationFlowShop& problem):emili::pfsp::HeavilyApproximatedTaillardAcceleratedInsertNeighborhood(problem) { }
};
/**
 * Six level approximation
 */
class SatxNeighborhood: public emili::pfsp::HeavilyApproximatedTaillardAcceleratedInsertNeighborhood
{
protected:
    virtual Solution* computeStep(Solution *value);
public:
    SatxNeighborhood(PermutationFlowShop& problem):emili::pfsp::HeavilyApproximatedTaillardAcceleratedInsertNeighborhood(problem) { }
};
/**
 * One level approximation Threshold testbed
 */
class TatxNeighborhood: public emili::pfsp::HeavilyApproximatedTaillardAcceleratedInsertNeighborhood
{
protected:
    int aptre;
    virtual Solution* computeStep(Solution *value);
public:
    TatxNeighborhood(float approximation_start_threshold, PermutationFlowShop& problem):emili::pfsp::HeavilyApproximatedTaillardAcceleratedInsertNeighborhood(problem),aptre(approximation_start_threshold*problem.getNjobs()) { }
};

/**  Neighborhoods based on approximation speed-up for other objectives
 *
 * */

/**
 * Total Completion Time
 * One level approximation with experimental performance improvement tricks
 */
class NatxTCTNeighborhood: public emili::pfsp::Natx2Neighborhood
{
protected:
    virtual Solution* computeStep(Solution *value);
public:
    NatxTCTNeighborhood(PermutationFlowShop& problem):emili::pfsp::Natx2Neighborhood(problem){ }
    NatxTCTNeighborhood(PermutationFlowShop& problem, int starting_threshold):emili::pfsp::Natx2Neighborhood(problem,starting_threshold) { }
};

/**
 * Total Completion Time
 * RZ kind of Neighborhood
 *
 */
class NrzTCTNeighborhood: public emili::pfsp::Natx2Neighborhood
{
protected:
    virtual Solution* computeStep(Solution *value);
    std::vector< int > seed_seq;
public:
    NrzTCTNeighborhood(PermutationFlowShop& problem):emili::pfsp::Natx2Neighborhood(problem){ }
    NrzTCTNeighborhood(PermutationFlowShop& problem, int starting_threshold):emili::pfsp::Natx2Neighborhood(problem,starting_threshold) { }
    virtual NeighborhoodIterator begin(Solution *base);
};

/**
 * Total Tardiness
 * One level approximation with experimental performance improvement tricks
 */
class NatxTTNeighborhood: public emili::pfsp::Natx2Neighborhood
{
protected:
    virtual Solution* computeStep(Solution *value);
public:
    NatxTTNeighborhood(PermutationFlowShop& problem):emili::pfsp::Natx2Neighborhood(problem){ }
    NatxTTNeighborhood(PermutationFlowShop& problem, int starting_threshold):emili::pfsp::Natx2Neighborhood(problem,starting_threshold) { }
};

class NoIdleAcceleratedInsertNeighborhood: public TaillardAcceleratedInsertNeighborhood
{
protected:
    virtual Solution* computeStep(Solution *value);    
public:
    NoIdleAcceleratedInsertNeighborhood(PermutationFlowShop& problem):TaillardAcceleratedInsertNeighborhood(problem) { }
    virtual NeighborhoodIterator begin(Solution *base);
};

class NoWaitAcceleratedNeighborhood: public PfspInsertNeighborhood
{
protected:
    const std::vector<std::vector < int > >& distance;
    const std::vector< std::vector< long int > >& pmatrix;
    const int nmac;
//    void computeNoWaitTimeDistances();
//    virtual Solution* computeStep(Solution *value);
public:
    NoWaitAcceleratedNeighborhood(NWPFSP_MS& problem):PfspInsertNeighborhood(problem),pmatrix(problem.getProcessingTimesMatrix()),distance(problem.getDistances()),nmac(problem.getNmachines()){ }
//    virtual NeighborhoodIterator begin(Solution *base);
    const std::vector<std::vector < int > >& getDistance() { return distance;}
};

class NoWaitAcceleratedInsertNeighborhood: public NoWaitAcceleratedNeighborhood
{
protected:
    virtual Solution* computeStep(Solution *value);
public:
    NoWaitAcceleratedInsertNeighborhood(NWPFSP_MS& problem):NoWaitAcceleratedNeighborhood(problem) { }
    virtual Solution* random(Solution *currentSolution);
};

class NoWaitAcceleratedTwoInsertNeighborhood: public NoWaitAcceleratedNeighborhood
{
protected:
    virtual Solution* computeStep(Solution *value);
    virtual void reverseLastMove(Solution *step);
public:
    NoWaitAcceleratedTwoInsertNeighborhood(NWPFSP_MS& problem):NoWaitAcceleratedNeighborhood(problem) { }
    virtual Solution* random(Solution *currentSolution);
};

class NoWaitAcceleratedExchangeNeighborhood: public NoWaitAcceleratedNeighborhood
{
protected:
    virtual Solution* computeStep(Solution *value);
    virtual void reverseLastMove(Solution *step);
public:
    NoWaitAcceleratedExchangeNeighborhood(NWPFSP_MS& problem):NoWaitAcceleratedNeighborhood(problem) { }
    virtual Solution* random(Solution *currentSolution);
};

class NoWaitAcceleratedTransposeNeighborhood: public NoWaitAcceleratedNeighborhood
{
protected:
    virtual Solution* computeStep(Solution *value);
    virtual void reverseLastMove(Solution *step);
public:
    NoWaitAcceleratedTransposeNeighborhood(NWPFSP_MS& problem):NoWaitAcceleratedNeighborhood(problem) { }
    virtual Solution* random(Solution *currentSolution);
};

/**
 * Insert neighborhood with Taillard's acceleration
 * that does a full scan each iteration
 */
class SDSTTaillardAcceleratedInsertNeighborhood: public emili::pfsp::TaillardAcceleratedInsertNeighborhood
{
protected:
    virtual Solution* computeStep(Solution *value);
    std::vector< std::vector < std::vector< int > > >& setUpTimes;
public:
    SDSTTaillardAcceleratedInsertNeighborhood(PermutationFlowShop& problem):emili::pfsp::TaillardAcceleratedInsertNeighborhood(problem),setUpTimes(problem.getInstance().getSetUpTimes()){ }
    virtual NeighborhoodIterator begin(Solution *base);
};

class TAxInsertNeighborhood: public emili::pfsp::PfspInsertNeighborhood
{
protected:
    std::vector < std::vector < int > > head;
    std::vector < std::vector < std::vector < int > > > tails;
    const std::vector < std::vector < long int > >& pmatrix;
    const int nmac;

    virtual Solution* computeStep(Solution *value);
public:
    TAxInsertNeighborhood(PermutationFlowShop& problem):emili::pfsp::PfspInsertNeighborhood(problem),head(problem.getNmachines()+1,std::vector< int > (problem.getNjobs()+1,0)),tails(problem.getNjobs()+1,std::vector< std::vector< int > >(problem.getNmachines()+1,std::vector< int >(problem.getNjobs()+1,0))),pmatrix(problem.getProcessingTimesMatrix()),nmac(pis.getNmachines()) { }
    virtual NeighborhoodIterator begin(Solution *base);
};

class PfspBackwardInsertNeighborhood: public PfspInsertNeighborhood
{
protected:
    virtual Solution* computeStep(Solution* value);
public:
    PfspBackwardInsertNeighborhood(PermutationFlowShop& problem):PfspInsertNeighborhood(problem) { }
    virtual Solution* random(Solution *currentSolution);
};

class PfspForwardInsertNeighborhood: public PfspInsertNeighborhood
{
protected:
    virtual Solution* computeStep(Solution* value);
public:    
    PfspForwardInsertNeighborhood(PermutationFlowShop& problem):PfspInsertNeighborhood(problem) { }
    virtual Solution* random(Solution *currentSolution);
};

class PfspTwoInsertNeighborhood: public PfspInsertNeighborhood
{
protected:
    virtual Solution* computeStep(Solution *value);
    virtual void reverseLastMove(Solution *step);
public:
    PfspTwoInsertNeighborhood(PermutationFlowShop& problem):PfspInsertNeighborhood(problem) { }
    virtual Solution* random(Solution *currentSolution);
};

class PfspExchangeNeighborhood: public emili::pfsp::PfspNeighborhood
{
protected:
    int start_position;
    int end_position;
    int sp_iterations;
    int ep_iterations;
    virtual Solution* computeStep(Solution* value);
    virtual void reverseLastMove(Solution *step);
public:
    PfspExchangeNeighborhood(PermutationFlowShop& problem):PfspNeighborhood(problem),start_position(0),end_position(0),sp_iterations(1),ep_iterations(1){}
    virtual void reset();
    virtual Solution* random(Solution *currentSolution);
    virtual Solution* random(Solution *currentSolution,int size);
    virtual std::pair<int,int> lastMove() { return std::pair<int,int>(end_position,start_position); }    
    virtual NeighborhoodIterator begin(Solution *base);
};

class AxtExchange: public emili::pfsp::PfspExchangeNeighborhood
{
protected:
    std::vector < std::vector < int > > head;
    const std::vector < std::vector < long int > >& pmatrix;
    const int nmac;
    int thresh;
    virtual Solution* computeStep(Solution *value);
    virtual void computeHead(std::vector<int>& sol);
public:
    AxtExchange(PermutationFlowShop& problem):emili::pfsp::PfspExchangeNeighborhood(problem),head(problem.getNmachines()+1,std::vector< int > (problem.getNjobs()+1,0)),pmatrix(problem.getProcessingTimesMatrix()),nmac(problem.getNmachines()),thresh(problem.getNjobs()/2) { }
    virtual NeighborhoodIterator begin(Solution *base);
};

class HaxtExchange: public emili::pfsp::AxtExchange
{
protected:
    int threshold;
   virtual Solution* computeStep(Solution *value);
public:
HaxtExchange(PermutationFlowShop& problem):emili::pfsp::AxtExchange(problem),threshold(problem.getNjobs()/2) {}
};

class EaxtExchange: public emili::pfsp::AxtExchange
{
protected:
   virtual Solution* computeStep(Solution *value);
public:
EaxtExchange(PermutationFlowShop& problem):emili::pfsp::AxtExchange(problem) {}
};


class OptExchange: public emili::pfsp::AxtExchange
{
protected:
    virtual Solution* computeStep(Solution *value);
public:
    OptExchange(PermutationFlowShop& problem):emili::pfsp::AxtExchange(problem) { }
};

class PfspTransposeNeighborhood: public emili::pfsp::PfspNeighborhood
{
protected:
    int start_position;
    int sp_iterations;    
    virtual Solution* computeStep(Solution* value);
    virtual void reverseLastMove(Solution *step);
public:
    PfspTransposeNeighborhood(PermutationFlowShop& problem):PfspNeighborhood(problem),start_position(0),sp_iterations(1){}
    virtual void reset();
    virtual Solution* random(Solution *currentSolution);
    virtual Solution* random(Solution *currentSolution,int size);
    virtual std::pair<int,int> lastMove() { return std::pair<int,int>(start_position+1,start_position); }
    virtual NeighborhoodIterator begin(Solution *base);
    virtual int size();
};

class XTransposeNeighborhood: public emili::pfsp::PfspTransposeNeighborhood
{
    //TODO FINISH THIS THING!!!!
protected:
    std::vector<int> prevJob;
    std::vector<int> previousMachineEndTime;
    int last_saved_position;
    virtual Solution* computeStep(Solution *value);
public:
    XTransposeNeighborhood(PermutationFlowShop& problem):emili::pfsp::PfspTransposeNeighborhood(problem),prevJob(problem.getNmachines()+1,0),previousMachineEndTime(problem.getNjobs()+1,0),last_saved_position(-1) { }
    virtual NeighborhoodIterator begin(Solution *base);
};
/**
 * @brief The InversionNeighborhood class
 * Neighborhood used in SDST
 * See Sioud Gagné 2018
 */
class InversionNeighborhood: public emili::pfsp::PfspTransposeNeighborhood
{
protected:
    virtual Solution* computeStep(Solution* value);
    virtual void reverseLastMove(Solution *step);
    int gsize;
public:
    InversionNeighborhood(PermutationFlowShop& problem, int group_size):
        emili::pfsp::PfspTransposeNeighborhood(problem),
        gsize(group_size)
        {}

    virtual Solution* random(Solution *currentSolution);
    virtual Solution* random(Solution *currentSolution,int size);
    virtual std::pair<int,int> lastMove()
            {return std::pair<int,int>(start_position,start_position+gsize);}
    virtual int size(){return pis.getNjobs()/gsize;}
};

/**
 * @brief The ScrambleNeighborhood class
 * Neighborhood used in SDST, scrambles randomly
 * a given group of jobs. Works similarly to InversionNeighborhood
 */
class ScrambleNeighborhood: public emili::pfsp::InversionNeighborhood
{
public:
    ScrambleNeighborhood(PermutationFlowShop& problem, int group_size):
        emili::pfsp::InversionNeighborhood(problem,group_size) {}
    virtual Solution* random(Solution *currentSolution);
    virtual Solution* random(Solution *currentSolution,int size);
};


/**
 * @brief The PfspTerminationClassic class
 */
class PfspTerminationClassic: public emili::Termination
{
public:
    PfspTerminationClassic() {}
    virtual bool terminate(Solution* currentSolution, Solution* newSolution);
    virtual void reset() { }
};

class PfspRandomSwapPertub: public emili::Perturbation
{
protected:
    PermutationFlowShop& pfs;
public:
    PfspRandomSwapPertub(PermutationFlowShop& problem_instance):pfs(problem_instance) { }
    virtual Solution* perturb(Solution* solution);

};

class PfspTestAcceptance: public emili::Acceptance
{
protected:
    PermutationFlowShop& pfs;
    int percentage;
public:
    PfspTestAcceptance(PermutationFlowShop& problem_instance):pfs(problem_instance),percentage(70) { }
    PfspTestAcceptance(PermutationFlowShop& problem_instance,int perc):pfs(problem_instance),percentage(perc) { }
    virtual Solution* accept(Solution* candidate1, Solution* candidate2);
};

class SOAacceptance: public emili::MetropolisAcceptance
{
public:
    SOAacceptance(float start_temp):emili::MetropolisAcceptance(start_temp) { }
    virtual Solution* accept(Solution* candidate1, Solution* candidate2);
};

class SOAtermination: public emili::Termination
{
protected:
    int numberOfSteps;
    int currentStep;
public:
    SOAtermination(int number_of_steps):numberOfSteps(number_of_steps),currentStep(0) { }
    virtual bool terminate(Solution* currentSolution, Solution* newSolution);
    void reset();
};

class PfspTerminationIterations: public emili::Termination
{
protected:
    int iterations;
    int maxIterations;
public:
    PfspTerminationIterations(int max_badIterations):maxIterations(max_badIterations),iterations(0) { }
    virtual bool terminate(Solution* currentSolution, Solution* newSolution);
    void reset();
};

class KarTermination: public emili::Termination
{
protected:
    int iterations;
    int maxIterations;
public:
    KarTermination(int max_iterations):maxIterations(max_iterations),iterations(0) { }
    virtual bool terminate(Solution *currentSolution, Solution *newSolution);
    virtual void reset();
};

class PfspTabuHashMemory: public emili::TabuMemory
{
protected:
    std::vector < size_t > tabuVector;
    int tt_index;
    size_t calc_hash(std::vector< int >& vec_sol);
public:
    PfspTabuHashMemory(int tabuTenure):emili::TabuMemory(tabuTenure),tabuVector(),tt_index(0) { }
    PfspTabuHashMemory():emili::TabuMemory(),tabuVector(),tt_index(0) { }
    /**
     * this method should return true if the solution is not tabu and false in the other case,
     */
    virtual bool tabu_check(Solution *solution);
    virtual void forbid(Solution* solution);
    virtual void reset();
};

class PfspTabuValueMemory: public emili::TabuMemory
{
    std::vector < double> tabuVector;
    int tt_index;
public:
    PfspTabuValueMemory(int tabuTenure):emili::TabuMemory(tabuTenure),tabuVector(),tt_index(0) { }
    PfspTabuValueMemory():emili::TabuMemory(),tabuVector(),tt_index(0) { }
    /**
     * this method should return true if the solution is not tabu and false in the other case,
     */
    virtual bool tabu_check(Solution *solution);
    virtual void forbid(Solution* solution);
    virtual void reset();

};

class PfspFullSolutionMemory: public emili::TabuMemory
{
  std::vector < std::vector < int > > tabuVector;

  int tt_index;
public:
    PfspFullSolutionMemory(int tabtenure):emili::TabuMemory(tabtenure),tt_index(0) { }
    PfspFullSolutionMemory():emili::TabuMemory(),tt_index(0) { }
    /**
     * this method should return true if the solution is not tabu and false in the other case,
     */
    virtual bool tabu_check(Solution *solution);
    virtual void forbid(Solution *solution);
    virtual void reset();
};

class PfspMovesMemory: public emili::TabuMemory
{
protected:
    std::vector < std::pair < int,int > > tabuVector;
    emili::pfsp::PfspNeighborhood* neigh;
    int tt_index;
    std::pair <int,int> lastMove;
    virtual bool tabu_check(std::pair< int,int > value);
  public:
      PfspMovesMemory(int tabtenure,emili::pfsp::PfspNeighborhood* n):emili::TabuMemory(tabtenure),tt_index(0),neigh(n),lastMove(0,0) { }
      PfspMovesMemory(emili::pfsp::PfspNeighborhood* n):emili::TabuMemory(),tt_index(0),neigh(n),lastMove(0,0) { }
      PfspMovesMemory(int tabtenure):emili::TabuMemory(tabtenure),tt_index(0),neigh(nullptr),lastMove(0,0) { }
      /**
       * this method should return true if the solution is not tabu and false in the other case,
       */
      virtual bool tabu_check(Solution *solution);
      virtual void forbid(Solution *solution);
      virtual void registerMove(emili::Solution* base,emili::Solution* solution);
      virtual void reset();
      virtual void setNeighborhood(Neighborhood *neighborhood) {neigh = (emili::pfsp::PfspNeighborhood*)neighborhood;}
};

class PfspMovesMemory2: public PfspMovesMemory
{
  public:
      PfspMovesMemory2(int tabtenure,emili::pfsp::PfspNeighborhood* n):PfspMovesMemory(tabtenure,n) { }
      PfspMovesMemory2(emili::pfsp::PfspNeighborhood* n):PfspMovesMemory(n) { }
      PfspMovesMemory2(int tabtenure):PfspMovesMemory(tabtenure) { }
      /**
       * this method should return true if the solution is not tabu and false in the other case,
       */
      virtual void registerMove(emili::Solution* base,emili::Solution* solution);
};


class TSABtestMemory: public emili::pfsp::PfspMovesMemory
{

protected:
    virtual bool tabu_check(std::pair<int, int> value);

public:
    TSABtestMemory(int tabtenure,emili::pfsp::PfspNeighborhood* n):emili::pfsp::PfspMovesMemory(tabtenure,n) { }
    TSABtestMemory(int tabtenure):emili::pfsp::PfspMovesMemory(tabtenure) { }
    TSABtestMemory(emili::pfsp::PfspNeighborhood* n):emili::pfsp::PfspMovesMemory(n) { }
    virtual void forbid(Solution *solution);
};

class TSABMemory: public emili::pfsp::PfspMovesMemory
{

protected:
    virtual bool tabu_check(std::pair<int, int> value,std::vector< int >& solution);

public:
    TSABMemory(int tabtenure,emili::pfsp::PfspNeighborhood* n):emili::pfsp::PfspMovesMemory(tabtenure,n) { }
    TSABMemory(emili::pfsp::PfspNeighborhood* n):emili::pfsp::PfspMovesMemory(n) { }
    TSABMemory(int tabtenure):emili::pfsp::PfspMovesMemory(tabtenure) { }
    virtual bool tabu_check(Solution *solution);
    virtual void forbid(Solution *solution);
};

class VNDBestSearch: public emili::BestImprovementSearch
{
protected:
    emili::BestImprovementSearch bs1;
    emili::BestImprovementSearch bs2;
    emili::BestImprovementSearch bs3;

public:
    VNDBestSearch(emili::InitialSolution& is,emili::Termination& tc,emili::pfsp::PfspNeighborhood& n1,emili::pfsp::PfspNeighborhood& n2,emili::pfsp::PfspNeighborhood& n3):emili::BestImprovementSearch(is,tc,n1),bs1(is,tc,n1),bs2(is,tc,n2),bs3(is,tc,n3) { }
    virtual emili::Solution* search(Solution* initial);
};

class VNDFirstSearch: public emili::FirstImprovementSearch
{
protected:
    emili::FirstImprovementSearch bs1;
    emili::FirstImprovementSearch bs2;
    emili::FirstImprovementSearch bs3;

public:
    VNDFirstSearch(emili::InitialSolution& is,emili::Termination& tc,emili::pfsp::PfspNeighborhood& n1,emili::pfsp::PfspNeighborhood& n2,emili::pfsp::PfspNeighborhood& n3):emili::FirstImprovementSearch(is,tc,n1),
        bs1(is,tc,n1),bs2(is,tc,n2),bs3(is,tc,n3) { }
    virtual emili::Solution* search(Solution* initial);
};


class GVNS_RIS_Neighborhood: public emili::pfsp::PfspNeighborhood
{
protected:
    emili::Solution* reference;
    int njobs;
    int index;
    virtual Solution* computeStep(Solution *step);
    virtual void reverseLastMove(Solution *step) { }
public:
    GVNS_RIS_Neighborhood(emili::pfsp::PermutationFlowShop& problem):emili::pfsp::PfspNeighborhood(problem),njobs(problem.getNjobs()),index(1) { }
    void setReference(emili::Solution* ref) {this->reference = ref;}
    virtual void reset();
    virtual Solution* random(Solution *currentSolution);
    virtual Solution* random(Solution *currentSolution,int size);
    virtual std::pair<int,int> lastMove() { return std::pair<int,int>(index,index); }
    virtual NeighborhoodIterator begin(Solution *base);
};

class GVNS_innerloop: public emili::LocalSearch
{
protected:
    emili::pfsp::GVNS_RIS_Neighborhood* rneigh;
public:
    GVNS_innerloop(InitialSolution& initialSolutionGenerator);
    virtual Solution* search(emili::Solution* initial);
};
/**
class IILS_neighborhood: public emili::pfsp::PfspInsertNeighborhood
{
protected:
    int ns_appl;
    int ns;
    int insmo;
    int ins;
    int CN;
    virtual Solution* computeStep(Solution *step);
    virtual void reverseLastMove(Solution *step);
public:
    IILS_neighborhood(emili::pfsp::PermutationFlowShop& prob):emili::pfsp::PfspInsertNeighborhood(prob),ns(0),insmo(0),ins(0),CN(0),ns_appl(0) { }
    virtual void reset();
    virtual NeighborhoodIterator begin(Solution *base);
};

class IILS_perturbation: public emili::Perturbation
{

};
*/
/**
 * Implementation of TSM algorithm from:
 * i, X., Chen, L., Xu, H., & Gupta, J. N. D. (2015). Trajectory scheduling methods
 *  for minimizing total tardiness in a flowshop.
 *  Operations Research Perspectives, 2, 13-23. doi:10.1016/j.orp.2014.12.001
 */

class CompoundPerturbation : public emili::Perturbation
{
protected:
    emili::pfsp::PermutationFlowShop& pis;
    int omega;
    float pc;
    int nbj;
    int d;
    int calc_distance(std::vector< int >& x, std::vector< int >& y);
    emili::pfsp::PfspInsertNeighborhood ins;
    emili::pfsp::PfspTransposeNeighborhood tra;
public:
    CompoundPerturbation(emili::pfsp::PermutationFlowShop& problem):pis(problem),ins(problem),tra(problem),nbj(problem.getNjobs()),omega(30),d(3),pc(0.2) {}
    CompoundPerturbation(emili::pfsp::PermutationFlowShop& problem,int phy_size,int number_of_perturbations, float perturbation_probability):pis(problem),ins(problem),tra(problem),nbj(problem.getNjobs()),omega(phy_size),d(number_of_perturbations),pc(perturbation_probability) {}
    virtual Solution* perturb(Solution *solution);


};


class CH6 : public emili::FirstImprovementSearch
{
protected:
    emili::Neighborhood* neigh1;
    emili::Neighborhood* neigh2;
public:
    CH6(emili::InitialSolution& is, emili::Termination& tc, emili::Neighborhood&  n1,emili::Neighborhood&  n2):emili::FirstImprovementSearch(is,tc,n1),neigh1(&n1),neigh2(&n2) { }
    //CH6(emili::BestImprovementSearch& ls, std::vector<emili::Neighborhood*> n):emili::BestImprovementSearch(ls),neigh(n) { }
    virtual emili::Solution* search(emili::Solution *initial)
    {
        emili::Solution* incumbent = initial->clone();
        emili::Solution* new_s;
        emili::Solution* new_s2 = incumbent->clone();
        int i = 1;
        while(i){
            //** bestSoFar = *incumbent;
            this->neighbh = neigh1;
            new_s = emili::FirstImprovementSearch::search(new_s2);
            this->neighbh = neigh2;
            delete new_s2;
            new_s2 = emili::FirstImprovementSearch::search(new_s);
            if(*new_s2 < *incumbent)
            {
                //delete incumbent;
                *incumbent = *new_s2;
                //std::cout << "xxxx ," << incumbent->getSolutionValue() << " " << emili::iteration_counter() << std::endl;
            }
            else
            {
                //delete new_s2;
                i = 0;
            }
            delete new_s;
        }
        return incumbent;
    }
};

class RIS : public emili::LocalSearch
{
protected:
    int ni_position;
    int njob;
    virtual int neh_ig(std::vector<int>& pi, int k);
    void invertPerturbation(std::vector<int>& pi, std::vector<int>& pi_i);
    emili::pfsp::PermutationFlowShop& instance;
public:
    RIS(emili::pfsp::PermutationFlowShop& problem, emili::InitialSolution& is):emili::LocalSearch(),instance(problem),njob(problem.getNjobs())
    {
        this->init = &is;
        this->neighbh = new emili::EmptyNeighBorHood();
        this->termcriterion = new emili::LocalMinimaTermination();
        this->bestSoFar = is.generateEmptySolution();
    }
    virtual emili::Solution* search(Solution *initial);
};

class NoIdle_RIS : public RIS
{
protected:
    std::vector < std::vector < int > > head;
    std::vector < std::vector < int > > tail;
    const std::vector < std::vector < long int > >& pmatrix;
    const int nmac;
    PfspInstance& pis;
    virtual int neh_ig(std::vector<int> &pi, int k);
public:
    NoIdle_RIS(emili::pfsp::PermutationFlowShop& problem, emili::InitialSolution& is):emili::pfsp::RIS(problem,is),head(problem.getNmachines()+1,std::vector< int > (problem.getNjobs()+1,0)),tail(problem.getNmachines()+1,std::vector< int >(problem.getNjobs()+1,0)),pmatrix(problem.getProcessingTimesMatrix()),nmac(problem.getNmachines()),pis(problem.getInstance()) { }
};

class NoIdleIGper : public RSPerturbation
{
public:
    NoIdleIGper(int d, emili::pfsp::PermutationFlowShop& problem):RSPerturbation(d,problem) {}
    virtual emili::Solution* perturb(Solution *solution);
};

class NoWait_RIS : public RIS
{
protected:
    const std::vector < std::vector < long int > >& pmatrix;
    const std::vector<std::vector < int > >& distance;
    const int nmac;
    PfspInstance& pis;
    virtual int neh_ig(std::vector<int> &pi, int k);
public:
    NoWait_RIS(NWPFSP_MS& problem, emili::InitialSolution& is):
        emili::pfsp::RIS(problem,is),
        distance(problem.getDistances()),
        nmac(problem.getNmachines()),
        pmatrix(problem.getProcessingTimesMatrix()),
        pis(problem.getInstance())
    { }
    const std::vector<std::vector < int > >& getDistance() { return distance;}
};

class RandomNoWait_RIS : public NoWait_RIS
{
protected:
    long cmax;
    emili::pfsp::PfspRandomInitialSolution rand;
public:
    RandomNoWait_RIS(NWPFSP_MS& problem, emili::InitialSolution& is):
        emili::pfsp::NoWait_RIS(problem,is),
        cmax(std::numeric_limits<long>::max()),
        rand(problem)
        { }
    virtual emili::Solution* search(Solution *initial);
};

class BeamSearchHeuristic : public PfspInitialSolution
{
protected:
    int _gamma;
    double _a;
    double _b;
    double _c;
    double _e;
    std::vector<double> xi;
    std::vector<int> xi_order;
    const std::vector< std::vector< long > >& pi;
    const std::vector<long>& dueDates;
    void buildXi();
    int njobs;
    int nmacs;
    virtual emili::Solution* generate();
public:
    BeamSearchHeuristic(PermutationFlowShop& problem, int gamma, double a, double b,double c,double e):
        emili::pfsp::PfspInitialSolution(problem),
        _gamma(gamma),
        _a(a),
        _b(b),
        _c(c),
        _e(e),
        pi(problem.getProcessingTimesMatrix()),
        xi(problem.getNjobs()+1,0),
        njobs(problem.getNjobs()),
        nmacs(problem.getNmachines()),
        dueDates(problem.getDueDates())
    {
        buildXi();        
    }

    class bs_node{
    protected:
        std::vector<int> scheduled;
        std::vector<int> unscheduled;
        std::vector< int > completionTimes;
       //std::vector< float > tpj;
       //double tpd;
        bs_node* father;

        std::vector< bs_node* > children;
        BeamSearchHeuristic& init;
        int m;
        int n;
        void evaluateNode();
        double calcW();
    public:
        int k;
        int kjob;
        double g_value;
        double TE;
        double TI;
        double TT;
        double T;
        double L;        
        bs_node(BeamSearchHeuristic& bs):
            father(nullptr),
            init(bs),
            completionTimes(bs.nmacs+1,0),            
            m(bs.nmacs),
            n(bs.njobs),
            g_value(0),
            k(0),
            L(0.0),
            T(0)
        {
            scheduled.push_back(0);
            std::vector<int>& xi = init.getXi_order();
            kjob = xi[0];
            scheduled.push_back(kjob);
            unscheduled = xi;
            unscheduled.erase(unscheduled.begin());
            for(int i=1;i<=m;i++)
                completionTimes[i] = completionTimes[i-1]+bs.pi[kjob][i];

            TE = std::max(((int)bs.dueDates[kjob])-completionTimes[m],0);
            TT = std::max(completionTimes[m]-((int)bs.dueDates[kjob]),0);
            TI = 0;
        }

        bs_node(BeamSearchHeuristic& bs,int start):
            father(nullptr),
            init(bs),
            completionTimes(bs.nmacs+1,0),
            m(bs.nmacs),
            n(bs.njobs),
            g_value(0),
            k(0),
            L(0.0),
            T(0)
        {
            scheduled.push_back(0);
            std::vector<int>& xi = init.getXi_order();

            kjob = xi[start];
            scheduled.push_back(kjob);
            unscheduled = xi;
            unscheduled.erase(unscheduled.begin()+start);
            for(int i=1;i<=m;i++)
                completionTimes[i] = completionTimes[i-1]+bs.pi[kjob][i];

            TE = std::max(((int)bs.dueDates[kjob])-completionTimes[m],0);
            TT = std::max(completionTimes[m]-((int)bs.dueDates[kjob]),0);
            TI = 0;
        }

        bs_node(bs_node& fat,int i):
            father(&fat),
            init(fat.init),
            m(fat.m),
            n(fat.n),
            completionTimes(fat.completionTimes),
            k(fat.k+1),
            TE(fat.TE),
            TT(fat.TT),
            TI(fat.TI),
            L(0.0),
            T(0)
        {
            scheduled = fat.scheduled;
            std::vector<int>& xi = fat.unscheduled;
            kjob = xi[i];           
            scheduled.push_back(kjob);
            unscheduled = xi;
            unscheduled.erase(unscheduled.begin()+i);
            evaluateNode();
        }

        bs_node(bs_node& node):
            father(node.father),
            init(node.init),
            m(node.m),
            n(node.n),
            completionTimes(node.completionTimes),           
            k(node.k),
            kjob(node.kjob),
            TE(node.TE),
            TT(node.TT),
            TI(node.TI),
            scheduled(node.scheduled),
            unscheduled(node.unscheduled),
            g_value(node.g_value),
            children(node.children)
            {         }
       void buildChildren();
      // bs_node* generateSequence();

       virtual bool operator < (bs_node& a)
       {
           return g_value < a.g_value;
       }

       virtual std::vector<int>& getPermutation()
       {
           return scheduled;
       }
       std::vector< bs_node* >& getChildren()
       {
           return this->children;
       }

      double calcG(double W);

       ~bs_node()
       {
          /*  std::vector<bs_node*>::iterator iter = children.begin();
            for(;iter!=children.end();++iter)
            {
                delete *iter;
            }*/
       }

    };

    virtual std::vector<double>& getXi();
    virtual std::vector<int>& getXi_order();
};

class FFheuristic : public PfspInitialSolution{
protected:
    int x;
    float a;
    float b;
    virtual Solution* generate();
public:
    FFheuristic(PermutationFlowShop& problem_instance, int num_sols, float a_opt, float b_opt):
        emili::pfsp::PfspInitialSolution(problem_instance),
        a(a_opt),
        b(b_opt),
        x(num_sols) { }


};

class BSCH : public PfspInitialSolution{
protected:
    int x;
    float ap;
    float bp;
    float cp;
    virtual Solution* generate();
public:
    BSCH(PermutationFlowShop& instance,int l,float a,float b,float c):
        emili::pfsp::PfspInitialSolution(instance),
        x(l),
        ap(a),
        bp(b),
        cp(c){ }
};

class BSheuristic : public PfspInitialSolution{
protected:
    int x;
    float ap;
    float bp;
    float cp;
    float ep;
    virtual Solution* generate();
public:
    BSheuristic(PermutationFlowShop& instance,int l,float a,float b,float c,float e):
        emili::pfsp::PfspInitialSolution(instance),
        x(l),
        ap(a),
        bp(b),
        cp(c),
        ep(e){ }
};

class SwapIncLocalSearch: public emili::LocalSearch
{
protected:
    int _r;
    emili::pfsp::PermutationFlowShop& p;
    int njobs;
public:
    SwapIncLocalSearch(int r,InitialSolution& is):
      emili::LocalSearch(),
      _r(r),
      p((PermutationFlowShop&)is.getProblem())
    {
        this->init = &is;
        this->neighbh = new emili::EmptyNeighBorHood();
        this->termcriterion = new emili::LocalMinimaTermination();
        this->bestSoFar = is.generateEmptySolution();
        this->njobs = p.getNjobs();
    }
    emili::Solution* search(emili::Solution* initial);
};

class STH : public emili::LocalSearch
{
protected:
    int b;
    const std::vector< std::vector <long int> >& processingTimesMatrix;
    std::vector< std::vector <long int> > sumOfST;
    emili::pfsp::SDSTFSP_MS& prob;
    void initialize_sumOfST();
public:
    STH(int _b, emili::pfsp::SDSTFSP_MS& problem,emili::InitialSolution& is):
        b(_b),
        prob(problem),
        processingTimesMatrix(problem.getProcessingTimesMatrix()),
        sumOfST(problem.getNjobs()+1,std::vector<long int>(problem.getNjobs()+1,0))
    {
        this->init = &is;
        this->neighbh = new emili::EmptyNeighBorHood();
        this->termcriterion = new emili::LocalMinimaTermination();
        this->bestSoFar = is.generateEmptySolution();
        initialize_sumOfST();
    }

    virtual emili::Solution* search(emili::Solution* initial);
};

}
}

#endif // PERMUTATIONFLOWSHOP_H
