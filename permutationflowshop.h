#ifndef PERMUTATIONFLOWSHOP_H
#define PERMUTATIONFLOWSHOP_H
#include "emilibase.h"
#include "pfspinstance.h"
#include <iostream>
#include <cstdlib>


namespace emili
{
namespace pfsp
{


class PermutationFlowShop: public emili::Problem
{
protected:
PfspInstance instance;
public:
    PermutationFlowShop(PfspInstance& problemInstance):instance(problemInstance)
    {

    }

    PermutationFlowShop(char* instance_path):instance()
    {
        /* Read data from file */
        if (! instance.readDataFromFile(instance_path) ){
            exit(-1);
        }
    }

    virtual double evaluateSolution(emili::Solution& solution);
    int getNjobs();
    int getDueDate(int job);
    int getPriority(int job);
    int computeMS(std::vector< int > & partial_solution);
    int computeWT(std::vector< int > & partial_solution);
    PfspInstance& getInstance();
};

class PermutationFlowShopSolution: public emili::Solution
{
protected:
    std::vector< int > solution;
public:
    PermutationFlowShopSolution(double p_value):emili::Solution(p_value),solution()
    {}

    PermutationFlowShopSolution(std::vector< int >& solution):emili::Solution(1e9),solution(solution)
    {}

    PermutationFlowShopSolution(double p_value,std::vector< int >& solution):emili::Solution(p_value),solution(solution)
    {}

    virtual const void* getRawData()const;
    virtual void setRawData(const void* data);
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
    /*The em*/
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

class SlackConstructor: public emili::Constructor
{
protected:
   PermutationFlowShop& pis;
public:
   SlackConstructor(PermutationFlowShop& problem):emili::Constructor(),pis(problem) { }
   virtual emili::Solution* construct(Solution *partial);
   virtual emili::Solution* constructFull();
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
    emili::pfsp::PermutationFlowShop instance;
public:
    PfspDestructor(emili::pfsp::PermutationFlowShop ist):instance(ist) { }
    virtual emili::Solution* destruct(Solution* solutioon);
};

class SOADestructor: public emili::Destructor
{
protected:
    int d;
    emili::pfsp::PermutationFlowShop instance;
public:
    SOADestructor(int d_parameter, emili::pfsp::PermutationFlowShop inst):d(d_parameter),instance(inst) {}
    virtual emili::Solution* destruct(Solution *solutioon);
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
    virtual PermutationFlowShopSolution* computeStep(std::vector< int > & solution,double value)=0;
    virtual Solution* computeStep(Solution* step)
    {
        return this->step(step);
    }

public:
    PfspNeighborhood(PermutationFlowShop& problem):pis(problem){}
    virtual Solution* step(Solution* currentSolution);
    virtual void reset();
    virtual std::pair<int,int> lastMove() { return std::pair<int,int>(0,0); }
};

class PfspBestImprovExchangeNeighborhood: public emili::pfsp::PfspNeighborhood
{
protected:
    int start_position;
    int end_position;
    int njobs;
    PfspInstance& instance;
    virtual PermutationFlowShopSolution* computeStep(std::vector< int > & solution,double value);
public:
    PfspBestImprovExchangeNeighborhood(PermutationFlowShop& problem):PfspNeighborhood(problem),start_position(1),end_position(1),njobs(problem.getNjobs()),instance(problem.getInstance()){}
    virtual Solution* random(Solution *currentSolution);

};

class PfspBestImprovInsertNeighborhood: public emili::pfsp::PfspNeighborhood
{
protected:
    int start_position;
    int end_position;
    int njobs;
    PfspInstance& instance;
    virtual PermutationFlowShopSolution* computeStep(std::vector< int > & solution,double value);
public:
    PfspBestImprovInsertNeighborhood(PermutationFlowShop& problem):PfspNeighborhood(problem),start_position(0),end_position(0),njobs(problem.getNjobs()),instance(problem.getInstance()){}
    virtual Solution* random(Solution *currentSolution);

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
    PfspInstance& instance;
    bool start;
    virtual PermutationFlowShopSolution* computeStep(std::vector<int> &solution,double value);
public:
    PfspInsertNeighborhood(PermutationFlowShop& problem):PfspNeighborhood(problem),start_position(0),end_position(0),njobs(problem.getNjobs()),instance(problem.getInstance()),sp_iterations(1),ep_iterations(1){}
    virtual void reset();
    virtual Solution* random(Solution *currentSolution);
    virtual std::pair<int,int> lastMove() { return std::pair<int,int>(end_position,start_position); }
    virtual NeighborhoodIterator begin(Solution *base);
};

class PfspExchangeNeighborhood: public emili::pfsp::PfspNeighborhood
{
protected:
    int start_position;
    int end_position;
    int sp_iterations;
    int ep_iterations;
    int njobs;
    PfspInstance& instance;
    virtual PermutationFlowShopSolution* computeStep(std::vector<int> &solution,double value);
public:
    PfspExchangeNeighborhood(PermutationFlowShop& problem):PfspNeighborhood(problem),start_position(0),end_position(0),njobs(problem.getNjobs()),instance(problem.getInstance()),sp_iterations(1),ep_iterations(1){}
    virtual void reset();
    virtual Solution* random(Solution *currentSolution);
    virtual std::pair<int,int> lastMove() { return std::pair<int,int>(end_position,start_position); }
    virtual NeighborhoodIterator begin(Solution *base);
};

class PfspTransposeNeighborhood: public emili::pfsp::PfspNeighborhood
{
protected:
    int start_position;
    int sp_iterations;
    int njobs;
    PfspInstance& instance;
    virtual PermutationFlowShopSolution* computeStep(std::vector< int > & solution,double value);
public:
    PfspTransposeNeighborhood(PermutationFlowShop& problem):PfspNeighborhood(problem),start_position(0),njobs(problem.getNjobs()),instance(problem.getInstance()),sp_iterations(1){}
    virtual void reset();
    virtual Solution* random(Solution *currentSolution);
    virtual std::pair<int,int> lastMove() { return std::pair<int,int>(start_position+1,start_position); }
    virtual NeighborhoodIterator begin(Solution *base);
};


class PfspBestImprovTransposeNeighborhood: public emili::pfsp::PfspNeighborhood
{
protected:
    int start_position;
    int end_position;
    int njobs;
    PfspInstance& instance;
    virtual PermutationFlowShopSolution* computeStep(std::vector< int > & solution,double value);
public:
    PfspBestImprovTransposeNeighborhood(PermutationFlowShop& problem):PfspNeighborhood(problem),start_position(0),end_position(0),njobs(problem.getNjobs()),instance(problem.getInstance()){}
    virtual Solution* random(Solution *currentSolution);
};

class PfspFirstImprovExchangeNeighborhood: public emili::pfsp::PfspNeighborhood
{
protected:
    int start_position;
    int end_position;
    virtual PermutationFlowShopSolution* computeStep(std::vector< int > & solution,double value);
public:
    PfspFirstImprovExchangeNeighborhood(PermutationFlowShop& problem):PfspNeighborhood(problem),start_position(1),end_position(1) { }
    virtual void reset();
    virtual Solution* random(Solution *currentSolution);
};

class PfspFirstImprovInsertNeighborhood: public emili::pfsp::PfspNeighborhood
{
protected:
    int start_position;
    int end_position;
    virtual PermutationFlowShopSolution* computeStep(std::vector< int > & solution,double value);
public:
    PfspFirstImprovInsertNeighborhood(PermutationFlowShop& problem):PfspNeighborhood(problem),start_position(1),end_position(1) { }
    virtual void reset();
    virtual Solution* random(Solution *currentSolution);
};

class PfspFirstImprovTransposeNeighborhood: public emili::pfsp::PfspNeighborhood
{
protected:
    int start_position;
    int end_position;
    virtual PermutationFlowShopSolution* computeStep(std::vector< int > & solution,double value);
public:
    PfspFirstImprovTransposeNeighborhood(PermutationFlowShop& problem):PfspNeighborhood(problem),start_position(1),end_position(1) { }
    virtual void reset();
    virtual Solution* random(Solution *currentSolution);
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
    PermutationFlowShop pfs;
public:
    PfspRandomSwapPertub(PermutationFlowShop problem_instance):pfs(problem_instance) { }
    virtual Solution* perturb(Solution* solution);

};

class PfspTestAcceptance: public emili::AcceptanceCriteria
{
protected:
    PermutationFlowShop pfs;
    int percentage;
public:
    PfspTestAcceptance(PermutationFlowShop problem_instance):pfs(problem_instance),percentage(70) { }
    PfspTestAcceptance(PermutationFlowShop problem_instance,int perc):pfs(problem_instance),percentage(perc) { }
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

class PfspTabuNeighborhood: public PfspNeighborhood, public TabuNeighborhood
{
protected:
    emili::Solution* computeStep(Solution *step)
    {
        return emili::pfsp::PfspNeighborhood::computeStep(step);
    }

public:
    PfspTabuNeighborhood(PermutationFlowShop& problem,int tabuTenureSize):emili::pfsp::PfspNeighborhood(problem),emili::TabuNeighborhood(tabuTenureSize) { }
    virtual Solution* step(Solution* currentSolution)
    {
        return emili::pfsp::PfspNeighborhood::step(currentSolution);
    }

    virtual void reset(){

    }
};

class PfspTabuInsertNeighborhood: public PfspTabuNeighborhood
{
protected:
    class NeighborhoodMove
    {
        int i;
        int j;
    public:
        NeighborhoodMove(int i_p,int j_p):i(i_p),j(j_p) { }
        bool operator ==( NeighborhoodMove& b)
        {
            return ((i == b.i) && (j==b.j)) || ((i == b.j) && (j==b.i)) ;
            //return (i+j)-(b.i+b.j) == 0;
            //what's faster?
        }
    };
    std::vector < NeighborhoodMove > tabuTable;
    int njobs;
    PfspInstance& instance;
    virtual emili::pfsp::PermutationFlowShopSolution* computeStep(std::vector< int > & solution,double value);
    virtual bool notTabu(NeighborhoodMove a);
    virtual void updateTabuTable(NeighborhoodMove a);
public:
    PfspTabuInsertNeighborhood(PermutationFlowShop& problem,int tabuTenureSize):emili::pfsp::PfspTabuNeighborhood(problem,tabuTenureSize),tabuTable(),njobs(problem.getNjobs()),instance(problem.getInstance()){}

};

class PfspTabuHashMemory: public emili::TabuMemory
{
protected:
    std::vector < size_t > tabuVector;
    int tt_index;
    size_t calc_hash(std::vector< int > * vec_sol);
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
    bool tabu_check(std::pair< int,int > value);
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
}
}

#endif // PERMUTATIONFLOWSHOP_H
