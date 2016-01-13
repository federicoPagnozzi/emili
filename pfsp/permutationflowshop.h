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

/*
 *
 *  Permutation Flow shop components for EMILI
 *
*/
namespace emili
{
namespace pfsp
{
/*Permutation Flowshop problem implementation:
  The class uses the implementation code originally from Jeremie ( see pfspinstance.h)
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
        /* Read data from file */
        if (! instance.readDataFromFile(instance_path) ){
            exit(-1);
        }
    }
    //implementation of evaluate solution
    virtual double evaluateSolution(emili::Solution& solution);
    /* This method returns the number of jobs*/
    int getNjobs();
    /* This method returns the number of machines*/
    int getNmachines();
    /* This method returns the due date given the job*/
    int getDueDate(int job);
    /* This method returns the priority given the job*/
    int getPriority(int job);
    /* This method returns the due dates for all the jobs*/
    std::vector< long int >& getDueDates();
    /* This method returns the priorities for all the jobs*/
    std::vector< long int >& getPriorities();
    /* This method returns the processing time matrix so that it can be used by particular implementations of neighborhoods*/
    const std::vector< std::vector < long int > > & getProcessingTimesMatrix();
    /* this method returns the problem size (used by some timed termination criteria)*/
    virtual int problemSize(){ return instance.getNbMac()*instance.getNbJob();}
    /* This method returns the pfspinstance object that incapsulate the actual computation of the objective functions*/
    PfspInstance& getInstance();
    /*this method computes the makespan for the given solution.
    * The method is here because there are some heuristics that use it.
    */
    int computeMS(std::vector< int > & partial_solution);
    int computeMS(std::vector< int >& partial, int size);
    /*
     * The classes that extends this class to implement a PFSP objective
     * should implement these methods to calculate the objective functions
    */
    virtual int computeObjectiveFunction(std::vector< int > & partial_solution)=0;
    virtual int computeObjectiveFunction(std::vector< int > & partial_solution, int size)=0;

    /*This methods compute the matrices to implement Taillard's acceleration*/
    void computeTAmatrices(std::vector<int> &sol,std::vector< std::vector < int > >& head, std::vector< std::vector< int > >& tail);
    void computeTAmatrices(std::vector<int> &sol,std::vector< std::vector < int > >& head, std::vector< std::vector< int > >& tail,int size);
    void computeNoIdleTAmatrices(std::vector<int> &sol,std::vector< std::vector < int > >& head, std::vector< std::vector< int > >& tail);

    /*Old methods used by some particular and exceptional speed-ups*/
    int computeObjectiveFunction(std::vector<int> &sol,std::vector<int>& prevJob,int job,std::vector<int>& previousMachineEndTime);
    int computeObjectiveFunction(std::vector< int > & sol,std::vector<std::vector<int > >& previousMachineEndTimeMatrix, int start_i, int end_i);
    void computeWTs(std::vector<int> &sol,std::vector<int>& prevJob,int job,std::vector<int>& previousMachineEndTime);
    void computeTails(std::vector<int> &sol, std::vector< std::vector< std::vector< int > > > & tails);
};
/* CLASSIC PERMUTATION FLOW SHOP*/
/*Weighted Tardiness*/
class PFSP_WT: public PermutationFlowShop
{
public:
    PFSP_WT(PfspInstance& problemInstance):PermutationFlowShop(problemInstance) { }
    PFSP_WT(char* instance_path):PermutationFlowShop(instance_path) { }
    virtual int computeObjectiveFunction(std::vector< int > & partial_solution);
    virtual int computeObjectiveFunction(std::vector< int > & partial_solution, int size);
};
/*Weighted completion time*/
class PFSP_WCT: public PermutationFlowShop
{
public:
    PFSP_WCT(PfspInstance& problemInstance):PermutationFlowShop(problemInstance) { }
    PFSP_WCT(char* instance_path):PermutationFlowShop(instance_path) { }
    virtual int computeObjectiveFunction(std::vector< int > & partial_solution);
    virtual int computeObjectiveFunction(std::vector< int > & partial_solution, int size);
};
/*Total completion time*/
class PFSP_TCT: public PermutationFlowShop
{
public:
    PFSP_TCT(PfspInstance& problemInstance):PermutationFlowShop(problemInstance) { }
    PFSP_TCT(char* instance_path):PermutationFlowShop(instance_path) { }
    virtual int computeObjectiveFunction(std::vector< int > & partial_solution);
    virtual int computeObjectiveFunction(std::vector< int > & partial_solution, int size);
};
/*Weighted Earliness*/
class PFSP_WE: public PermutationFlowShop
{
public:
    PFSP_WE(PfspInstance& problemInstance):PermutationFlowShop(problemInstance) { }
    PFSP_WE(char* instance_path):PermutationFlowShop(instance_path) { }
    virtual int computeObjectiveFunction(std::vector< int > & partial_solution);
    virtual int computeObjectiveFunction(std::vector< int > & partial_solution, int size);
};
/*Tardiness*/
class PFSP_T: public PermutationFlowShop
{
public:
    PFSP_T(PfspInstance& problem_instance):PermutationFlowShop(problem_instance) { }
    PFSP_T(char* instance_path):PermutationFlowShop(instance_path) { }
    virtual int computeObjectiveFunction(std::vector<int> &partial_solution);
    virtual int computeObjectiveFunction(std::vector<int> &partial_solution,int size);
};
/*Earliness*/
class PFSP_E: public PermutationFlowShop
{
public:
    PFSP_E(PfspInstance& problem_instance):PermutationFlowShop(problem_instance) { }
    PFSP_E(char* instance_path):PermutationFlowShop(instance_path) { }
    virtual int computeObjectiveFunction(std::vector<int> &partial_solution);
    virtual int computeObjectiveFunction(std::vector<int> &partial_solution,int size);
};
/*Make span*/
class PFSP_MS: public PermutationFlowShop
{
public:
    PFSP_MS(PfspInstance& problem_instance):PermutationFlowShop(problem_instance) { }
    PFSP_MS(char* instance_path):PermutationFlowShop(instance_path) { }
    virtual int computeObjectiveFunction(std::vector<int> &partial_solution,int size);
    virtual int computeObjectiveFunction(std::vector<int> &partial_solution);
};


/* NO WAIT PERMUTATION FLOW SHOP*/

/*Make span*/
class NWPFSP_MS: public PermutationFlowShop
{
public:
    NWPFSP_MS(PfspInstance& problem_instance):PermutationFlowShop(problem_instance) { }
    NWPFSP_MS(char* instance_path):PermutationFlowShop(instance_path) { }
    virtual int computeObjectiveFunction(std::vector<int> &partial_solution,int size);
    virtual int computeObjectiveFunction(std::vector<int> &partial_solution);
};
/*Weighted Tardiness*/
class NWPFSP_WT: public PermutationFlowShop
{
public:
    NWPFSP_WT(PfspInstance& problem_instance):PermutationFlowShop(problem_instance) { }
    NWPFSP_WT(char* instance_path):PermutationFlowShop(instance_path) { }
    virtual int computeObjectiveFunction(std::vector<int> &partial_solution,int size);
    virtual int computeObjectiveFunction(std::vector<int> &partial_solution);
};
/*Weighted Earliness*/
class NWPFSP_WE: public PermutationFlowShop
{
public:
    NWPFSP_WE(PfspInstance& problem_instance):PermutationFlowShop(problem_instance) { }
    NWPFSP_WE(char* instance_path):PermutationFlowShop(instance_path) { }
    virtual int computeObjectiveFunction(std::vector<int> &partial_solution,int size);
    virtual int computeObjectiveFunction(std::vector<int> &partial_solution);
};
/*Tardiness*/
class NWPFSP_T: public PermutationFlowShop
{
public:
    NWPFSP_T(PfspInstance& problem_instance):PermutationFlowShop(problem_instance) { }
    NWPFSP_T(char* instance_path):PermutationFlowShop(instance_path) { }
    virtual int computeObjectiveFunction(std::vector<int> &partial_solution,int size);
    virtual int computeObjectiveFunction(std::vector<int> &partial_solution);
};
/*Earliness*/
class NWPFSP_E: public PermutationFlowShop
{
public:
    NWPFSP_E(PfspInstance& problem_instance):PermutationFlowShop(problem_instance) { }
    NWPFSP_E(char* instance_path):PermutationFlowShop(instance_path) { }
    virtual int computeObjectiveFunction(std::vector<int> &partial_solution,int size);
    virtual int computeObjectiveFunction(std::vector<int> &partial_solution);
};
/*Weighted completion time*/
class NWPFSP_WCT: public PermutationFlowShop
{
public:
    NWPFSP_WCT(PfspInstance& problem_instance):PermutationFlowShop(problem_instance) { }
    NWPFSP_WCT(char* instance_path):PermutationFlowShop(instance_path) { }
    virtual int computeObjectiveFunction(std::vector<int> &partial_solution,int size);
    virtual int computeObjectiveFunction(std::vector<int> &partial_solution);
};
/*Total completion time*/
class NWPFSP_TCT: public PermutationFlowShop
{
public:
    NWPFSP_TCT(PfspInstance& problem_instance):PermutationFlowShop(problem_instance) { }
    NWPFSP_TCT(char* instance_path):PermutationFlowShop(instance_path) { }
    virtual int computeObjectiveFunction(std::vector<int> &partial_solution,int size);
    virtual int computeObjectiveFunction(std::vector<int> &partial_solution);
};

/* NO IDLE PERMUTATION FLOW SHOP*/

/*Make span*/
class NIPFSP_MS: public PermutationFlowShop
{
public:
    NIPFSP_MS(PfspInstance& problem_instance):PermutationFlowShop(problem_instance) { }
    NIPFSP_MS(char* instance_path):PermutationFlowShop(instance_path) { }
    virtual int computeObjectiveFunction(std::vector<int> &partial_solution);
    virtual int computeObjectiveFunction(std::vector<int> &partial_solution,int size);
};
/*Make span with accelerations*/
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
/*Weighted Tardiness*/
class NIPFSP_WT: public PermutationFlowShop
{
public:
    NIPFSP_WT(PfspInstance& problem_instance):PermutationFlowShop(problem_instance) { }
    NIPFSP_WT(char* instance_path):PermutationFlowShop(instance_path) { }
    virtual int computeObjectiveFunction(std::vector<int> &partial_solution);
    virtual int computeObjectiveFunction(std::vector<int> &partial_solution,int size);
};
/*Weighted Earliness*/
class NIPFSP_WE: public PermutationFlowShop
{
public:
    NIPFSP_WE(PfspInstance& problem_instance):PermutationFlowShop(problem_instance) { }
    NIPFSP_WE(char* instance_path):PermutationFlowShop(instance_path) { }
    virtual int computeObjectiveFunction(std::vector<int> &partial_solution);
    virtual int computeObjectiveFunction(std::vector<int> &partial_solution,int size);
};
/*Tardiness*/
class NIPFSP_T: public PermutationFlowShop
{
public:
    NIPFSP_T(PfspInstance& problem_instance):PermutationFlowShop(problem_instance) { }
    NIPFSP_T(char* instance_path):PermutationFlowShop(instance_path) { }
    virtual int computeObjectiveFunction(std::vector<int> &partial_solution);
    virtual int computeObjectiveFunction(std::vector<int> &partial_solution,int size);
};
/*Earliness*/
class NIPFSP_E: public PermutationFlowShop
{
public:
    NIPFSP_E(PfspInstance& problem_instance):PermutationFlowShop(problem_instance) { }
    NIPFSP_E(char* instance_path):PermutationFlowShop(instance_path) { }
    virtual int computeObjectiveFunction(std::vector<int> &partial_solution);
    virtual int computeObjectiveFunction(std::vector<int> &partial_solution,int size);
};
/*Weighted completion time*/
class NIPFSP_WCT: public PermutationFlowShop
{
public:
    NIPFSP_WCT(PfspInstance& problem_instance):PermutationFlowShop(problem_instance) { }
    NIPFSP_WCT(char* instance_path):PermutationFlowShop(instance_path) { }
    virtual int computeObjectiveFunction(std::vector<int> &partial_solution);
    virtual int computeObjectiveFunction(std::vector<int> &partial_solution,int size);
};
/*Total completion time*/
class NIPFSP_TCT: public PermutationFlowShop
{
public:
    NIPFSP_TCT(PfspInstance& problem_instance):PermutationFlowShop(problem_instance) { }
    NIPFSP_TCT(char* instance_path):PermutationFlowShop(instance_path) { }
    virtual int computeObjectiveFunction(std::vector<int> &partial_solution);
    virtual int computeObjectiveFunction(std::vector<int> &partial_solution,int size);
};

/*Sequence dependent setup times Permutation flowshop*/

/*Make span*/
class SDSTFSP_MS: public PermutationFlowShop
{
public:
    SDSTFSP_MS(PfspInstance& problem_instance):PermutationFlowShop(problem_instance) { }
    SDSTFSP_MS(char* instance_path):PermutationFlowShop(instance_path) { }
    virtual int computeObjectiveFunction(std::vector<int> &partial_solution);
    virtual int computeObjectiveFunction(std::vector<int> &partial_solution,int size);
};
/*Weighted Tardiness*/
class SDSTFSP_WT: public PermutationFlowShop
{
public:
    SDSTFSP_WT(PfspInstance& problem_instance):PermutationFlowShop(problem_instance) { }
    SDSTFSP_WT(char* instance_path):PermutationFlowShop(instance_path) { }
    virtual int computeObjectiveFunction(std::vector<int> &partial_solution);
    virtual int computeObjectiveFunction(std::vector<int> &partial_solution,int size);
};
/*Weighted Earliness*/
class SDSTFSP_WE: public PermutationFlowShop
{
public:
    SDSTFSP_WE(PfspInstance& problem_instance):PermutationFlowShop(problem_instance) { }
    SDSTFSP_WE(char* instance_path):PermutationFlowShop(instance_path) { }
    virtual int computeObjectiveFunction(std::vector<int> &partial_solution);
    virtual int computeObjectiveFunction(std::vector<int> &partial_solution,int size);
};
/*Tardiness*/
class SDSTFSP_T: public PermutationFlowShop
{
public:
    SDSTFSP_T(PfspInstance& problem_instance):PermutationFlowShop(problem_instance) { }
    SDSTFSP_T(char* instance_path):PermutationFlowShop(instance_path) { }
    virtual int computeObjectiveFunction(std::vector<int> &partial_solution);
    virtual int computeObjectiveFunction(std::vector<int> &partial_solution,int size);
};
/*Earliness*/
class SDSTFSP_E: public PermutationFlowShop
{
public:
    SDSTFSP_E(PfspInstance& problem_instance):PermutationFlowShop(problem_instance) { }
    SDSTFSP_E(char* instance_path):PermutationFlowShop(instance_path) { }
    virtual int computeObjectiveFunction(std::vector<int> &partial_solution);
    virtual int computeObjectiveFunction(std::vector<int> &partial_solution,int size);
};
/*Total completion time*/
class SDSTFSP_TCT: public PermutationFlowShop
{
public:
    SDSTFSP_TCT(PfspInstance& problem_instance):PermutationFlowShop(problem_instance) { }
    SDSTFSP_TCT(char* instance_path):PermutationFlowShop(instance_path) { }
    virtual int computeObjectiveFunction(std::vector<int> &partial_solution);
    virtual int computeObjectiveFunction(std::vector<int> &partial_solution,int size);
};
/*Weighted completion time*/
class SDSTFSP_WCT: public PermutationFlowShop
{
public:
    SDSTFSP_WCT(PfspInstance& problem_instance):PermutationFlowShop(problem_instance) { }
    SDSTFSP_WCT(char* instance_path):PermutationFlowShop(instance_path) { }
    virtual int computeObjectiveFunction(std::vector<int> &partial_solution);
    virtual int computeObjectiveFunction(std::vector<int> &partial_solution,int size);
};


/*This class implements the Solution for the Permutation FlowShop problem
  It uses astd::vector of ints for storing the job sequence.
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
    /*Returns the job sequence that rapresents this solution*/
    virtual std::vector< int >& getJobSchedule();
    /*Returns a printable version of the job sequence*/
    virtual std::string getSolutionRepresentation();
    /*Implements the clone method of emili::Solution*/
    virtual emili::Solution* clone();
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
    /*This method generates a new empty solution by instantiating an emptystd::vector of int of the correct size
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

class NEHff: public emili::pfsp::PfspInitialSolution
{
protected:
    virtual Solution* generate();
public:
    NEHff(PermutationFlowShop &problem_instance):emili::pfsp::PfspInitialSolution(problem_instance) {}
};

/*Less idle times construction heuristic from
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

class SlackConstructor: public emili::Constructor
{
protected:
   PermutationFlowShop& pis;
public:
   SlackConstructor(PermutationFlowShop& problem):emili::Constructor(),pis(problem) { }
   virtual emili::Solution* construct(Solution *partial);
   virtual emili::Solution* constructFull();
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

class NEHSlackConstructor: public emili::Constructor
{
protected:
   PermutationFlowShop& pis;
public:
   NEHSlackConstructor(PermutationFlowShop& problem):emili::Constructor(),pis(problem) { }
   virtual emili::Solution* construct(Solution *partial);
   virtual emili::Solution* constructFull();
};

class PfspDestructor: public emili::Destructor
{
protected:
    emili::pfsp::PermutationFlowShop& instance;
public:
    PfspDestructor(emili::pfsp::PermutationFlowShop& ist):instance(ist) { }
    virtual emili::Solution* destruct(Solution* solutioon);
};

class SOADestructor: public emili::Destructor
{
protected:
    int d;
    emili::pfsp::PermutationFlowShop& instance;
public:
    SOADestructor(int d_parameter, emili::pfsp::PermutationFlowShop& inst):d(d_parameter),instance(inst) {}
    virtual emili::Solution* destruct(Solution *solutioon);
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

class IGPerturbation: public emili::Perturbation
{
protected:
    int d;
    emili::pfsp::PermutationFlowShop& instance;
public:
    IGPerturbation(int d_parameter, emili::pfsp::PermutationFlowShop& problem):d(d_parameter),instance(problem) { }
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
    IgLsPerturbation(int d_parameter, emili::pfsp::PermutationFlowShop& problem, emili::LocalSearch* ls): emili::pfsp::IGPerturbation(d_parameter,problem),ls(ls) {/*   */}
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
};

class PfspDestructorTest: public emili::Destructor
{
protected:
    PermutationFlowShop& istance;
public:
    PfspDestructorTest(PermutationFlowShop& instance):istance(instance) { }
    virtual emili::Solution* destruct(Solution* solutioon);
};


class PfspNeighborhood: public emili::Neighborhood
{
protected:
    PermutationFlowShop& pis;    
   virtual Solution* computeStep(Solution* step) =0;
public:
    PfspNeighborhood(PermutationFlowShop& problem):pis(problem){}
    virtual Solution* step(Solution* currentSolution);
    virtual void reset();
    virtual std::pair<int,int> lastMove() { return std::pair<int,int>(0,0); }
    virtual int size();
};


class PfspInsertNeighborhood: public emili::pfsp::PfspNeighborhood
{
protected:
    int start_position;
    int end_position;
    int sp_iterations;
    int ep_iterations;
    int njobs;
    std::vector < int > current;
    int current_value;
    virtual Solution* computeStep(Solution* value);
    virtual void reverseLastMove(Solution *step);
public:
    PfspInsertNeighborhood(PermutationFlowShop& problem):PfspNeighborhood(problem),start_position(0),end_position(0),njobs(problem.getNjobs()),sp_iterations(1),ep_iterations(1){}
    virtual void reset();
    virtual Solution* random(Solution *currentSolution);
    virtual std::pair<int,int> lastMove() { return std::pair<int,int>(end_position,start_position); }
    virtual NeighborhoodIterator begin(Solution *base);
};

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


class HeavilyApproximatedTaillardAcceleratedInsertNeighborhood: public emili::pfsp::TaillardAcceleratedInsertNeighborhood
{
protected:
    virtual void computeHead(std::vector<int>& sol);
    virtual Solution* computeStep(Solution *value);
public:
    HeavilyApproximatedTaillardAcceleratedInsertNeighborhood(PermutationFlowShop& problem):emili::pfsp::TaillardAcceleratedInsertNeighborhood(problem) { }
    virtual NeighborhoodIterator begin(Solution *base);
};

/* This Insert recomputes the objective function value only for the modified parts of the solution
 *
 * */
class OptInsert: public emili::pfsp::TaillardAcceleratedInsertNeighborhood
{
protected:
    virtual Solution* computeStep(Solution *value);
public:
    OptInsert(PermutationFlowShop& problem):emili::pfsp::TaillardAcceleratedInsertNeighborhood(problem) { }
};

class NatxNeighborhood: public emili::pfsp::HeavilyApproximatedTaillardAcceleratedInsertNeighborhood
{
protected:    
    virtual Solution* computeStep(Solution *value);
public:
    NatxNeighborhood(PermutationFlowShop& problem):emili::pfsp::HeavilyApproximatedTaillardAcceleratedInsertNeighborhood(problem) { }
};

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

class AtxNeighborhood: public emili::pfsp::HeavilyApproximatedTaillardAcceleratedInsertNeighborhood
{
protected:
    virtual Solution* computeStep(Solution *value);
public:
    AtxNeighborhood(PermutationFlowShop& problem):emili::pfsp::HeavilyApproximatedTaillardAcceleratedInsertNeighborhood(problem) { }
};

class EatxNeighborhood: public emili::pfsp::HeavilyApproximatedTaillardAcceleratedInsertNeighborhood
{
protected:
    virtual Solution* computeStep(Solution *value);
public:
    EatxNeighborhood(PermutationFlowShop& problem):emili::pfsp::HeavilyApproximatedTaillardAcceleratedInsertNeighborhood(problem) { }
};

class TatxNeighborhood: public emili::pfsp::HeavilyApproximatedTaillardAcceleratedInsertNeighborhood
{
protected:
    float aptre;
    int aplev;
    virtual Solution* computeStep(Solution *value);
public:
    TatxNeighborhood(float approximation_start_threshold, int approximation_level, PermutationFlowShop& problem):emili::pfsp::HeavilyApproximatedTaillardAcceleratedInsertNeighborhood(problem),aplev(problem.getNmachines()-approximation_level),aptre(approximation_start_threshold*problem.getNjobs()) { }
};

class NoIdleAcceleratedInsertNeighborhood: public TaillardAcceleratedInsertNeighborhood
{
protected:
    virtual Solution* computeStep(Solution *value);    
public:
    NoIdleAcceleratedInsertNeighborhood(PermutationFlowShop& problem):TaillardAcceleratedInsertNeighborhood(problem) { }
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
    int njobs;
    virtual Solution* computeStep(Solution* value);
    virtual void reverseLastMove(Solution *step);
public:
    PfspExchangeNeighborhood(PermutationFlowShop& problem):PfspNeighborhood(problem),start_position(0),end_position(0),njobs(problem.getNjobs()),sp_iterations(1),ep_iterations(1){}
    virtual void reset();
    virtual Solution* random(Solution *currentSolution);
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
    int njobs;    
    virtual Solution* computeStep(Solution* value);
    virtual void reverseLastMove(Solution *step);
public:
    PfspTransposeNeighborhood(PermutationFlowShop& problem):PfspNeighborhood(problem),start_position(0),njobs(problem.getNjobs()),sp_iterations(1){}
    virtual void reset();
    virtual Solution* random(Solution *currentSolution);
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

class PfspTabuHashMemory: public emili::TabuMemory
{
protected:
    std::vector < size_t > tabuVector;
    int tt_index;
    size_t calc_hash(std::vector< int >& vec_sol);
public:
    PfspTabuHashMemory(int tabuTenure):emili::TabuMemory(tabuTenure),tabuVector(),tt_index(0) { }
    PfspTabuHashMemory():emili::TabuMemory(),tabuVector(),tt_index(0) { }
    /*
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
    /*
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
    /*
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
    emili::pfsp::PfspNeighborhood& neigh;
    int tt_index;
    std::pair <int,int> lastMove;
    virtual bool tabu_check(std::pair< int,int > value);
  public:
      PfspMovesMemory(int tabtenure,emili::pfsp::PfspNeighborhood& n):emili::TabuMemory(tabtenure),tt_index(0),neigh(n),lastMove(0,0) { }
      PfspMovesMemory(emili::pfsp::PfspNeighborhood& n):emili::TabuMemory(),tt_index(0),neigh(n),lastMove(0,0) { }
      /*
       * this method should return true if the solution is not tabu and false in the other case,
       */
      virtual bool tabu_check(Solution *solution);
      virtual void forbid(Solution *solution);
      virtual void registerMove(emili::Solution* base,emili::Solution* solution);
      virtual void reset();
};

class TSABtestMemory: public emili::pfsp::PfspMovesMemory
{

protected:
    virtual bool tabu_check(std::pair<int, int> value);

public:
    TSABtestMemory(int tabtenure,emili::pfsp::PfspNeighborhood& n):emili::pfsp::PfspMovesMemory(tabtenure,n) { }
    TSABtestMemory(emili::pfsp::PfspNeighborhood& n):emili::pfsp::PfspMovesMemory(n) { }
    virtual void forbid(Solution *solution);
};

class TSABMemory: public emili::pfsp::PfspMovesMemory
{

protected:
    virtual bool tabu_check(std::pair<int, int> value,std::vector< int >& solution);

public:
    TSABMemory(int tabtenure,emili::pfsp::PfspNeighborhood& n):emili::pfsp::PfspMovesMemory(tabtenure,n) { }
    TSABMemory(emili::pfsp::PfspNeighborhood& n):emili::pfsp::PfspMovesMemory(n) { }
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

}
}

#endif // PERMUTATIONFLOWSHOP_H