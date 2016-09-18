#ifndef IRP_H
#define IRP_H
#include "emilibase.h"
#include "ROADF/instance.h"

namespace emili
{
namespace irp
{

class InventoryRoutingProblem: public emili::Problem
{
protected:
    /*Internal structure of the invetory routing problem*/
    Instance irpInstance;
    clock_t time;
public:
    InventoryRoutingProblem(char* instance_path);
    /*implement this method to get the objective Function value*/
    virtual double evaluateSolution(Solution & solution);
//    virtual double evaluatePartialSolution(Solution &oldSolution, Solution &newSolution, unsigned int shiftIndex);

    /*Methods that are needed  for IRP */
    Instance getIrpInstance();
    double getTime(){return this->time;}
};

class InventoryRoutingSolution: public emili::Solution
{
protected:
    irpSolution irps;
    /*Implement these methods to make the copy constructor works*/
    virtual const void* getRawData()const;
    virtual void setRawData(const void* data);


    /*Internal structure of the Solution*/

public:
    /*Base constructor*/
    InventoryRoutingSolution(double solution_value):emili::Solution(solution_value) { }
    InventoryRoutingSolution(irpSolution irps):emili::Solution(DBL_MAX) {this->irps = irps;}
    /*This method must build a clone of the Solution on the heap and returns the pointer to it*/
    virtual Solution* clone();

    /*Methods needed for the Internal use of the Solution
     * ( how the solution is seen by Neighborhoods, Perturbations, Initial Solution generators... etc)*/
    irpSolution &getIrpSolution();

    virtual std::string getSolutionRepresentation();
};

class GreedyInitialSolution : public emili::InitialSolution
{
private:
    double randomFactor;
    unsigned int urgencyPolicy;
public:
    GreedyInitialSolution(Problem& instance, double rf, unsigned int up):emili::InitialSolution(instance), randomFactor(rf), urgencyPolicy(up){}
    virtual Solution* generateSolution();
    virtual Solution* generateEmptySolution();
};

class GreedyRandomizedInitialSolution : public emili::InitialSolution
{
private:
    double randomFactor;
    unsigned int urgencyPolicy;
public:
 GreedyRandomizedInitialSolution(Problem& instance, double rf, unsigned int up):emili::InitialSolution(instance), randomFactor(rf), urgencyPolicy(up){}    virtual Solution* generateSolution();
    virtual Solution* generateEmptySolution();
};

class GRASP : public emili::InitialSolution
{
private:
    double randomFactor;
    unsigned int urgencyPolicy;
public:
 GRASP(Problem& instance, double rf, unsigned int up):emili::InitialSolution(instance), randomFactor(rf), urgencyPolicy(up){}    virtual Solution* generateSolution();
    virtual Solution* generateEmptySolution();
};



class irpShiftTwoExchangeNeighborhood : public emili::Neighborhood
{
protected:
    InventoryRoutingProblem irp;
    Solution * currentNeighboringSolution;
    int shift;
    int operation;
    int numberOfOperations;
    int numberOfShifts;
    double randomFactor;
    unsigned int urgencyPolicy;
    double refuelFactor;
    double deliveredQuantityFactor;
    int numberFeasibleSolutions; //in a neighborhood
    double bestValueFound;
    virtual Solution* computeStep(Solution* step);
    virtual void reverseLastMove(Solution* step);
public:

     /*
      * Initialize shift and operation size, 4 indexes
      * this method needs to be overidden if there are things that a neighborhood has to do or reset when begin is called*/
     virtual NeighborhoodIterator begin(emili::Solution* base);

    /*this method returns a solution in the decided neighborhood
     * of the currentSolution */
    virtual Solution* step(Solution* currentSolution);
    /*
     * The state of the Neighborhood object may need to be restored to
     * initial conditions between local search calls
     * (e.g. first improvement strategies for permutation flow shop).
     */
    virtual void reset();
    /*
     * A method that returns a random solution in the neighborhood has to be provided
     */
    virtual Solution* random(Solution* currentSolution);
    /*
     * This method returns the size of the neighborhood
    */
    virtual int size();

    irpShiftTwoExchangeNeighborhood(InventoryRoutingProblem &i, double rf, unsigned int up, double ref, double dqf):
        irp(i),numberFeasibleSolutions(0),bestValueFound(DBL_MAX),
        randomFactor(rf),urgencyPolicy(up), refuelFactor(ref), deliveredQuantityFactor(dqf) {}
};

class irpShiftInsertNeighborhood : public emili::Neighborhood
{
protected:
    InventoryRoutingProblem irp;
    Solution * currentNeighboringSolution;
    int shift;
    int operation;
    int numberOfOperations;
    int numberOfShifts;
    double randomFactor;
    unsigned int urgencyPolicy;
    double refuelFactor;
    double deliveredQuantityFactor;
    int numberFeasibleSolutions; //in a neighborhood
    double bestValueFound;
    virtual Solution* computeStep(Solution* step);
    virtual void reverseLastMove(Solution* step);
public:

     /*
      * Initialize shift and operation size, 4 indexes
      * this method needs to be overidden if there are things that a neighborhood has to do or reset when begin is called*/
     virtual NeighborhoodIterator begin(emili::Solution* base);

    /*this method returns a solution in the decided neighborhood
     * of the currentSolution */
    virtual Solution* step(Solution* currentSolution);
    /*
     * The state of the Neighborhood object may need to be restored to
     * initial conditions between local search calls
     * (e.g. first improvement strategies for permutation flow shop).
     */
    virtual void reset();
    /*
     * A method that returns a random solution in the neighborhood has to be provided
     */
    virtual Solution* random(Solution* currentSolution);
    /*
     * This method returns the size of the neighborhood
    */
    virtual int size();

    irpShiftInsertNeighborhood(InventoryRoutingProblem &i, double rf, unsigned int up, double ref, double dqf):
        irp(i),numberFeasibleSolutions(0),bestValueFound(DBL_MAX),
        randomFactor(rf),urgencyPolicy(up),refuelFactor(ref), deliveredQuantityFactor(dqf){}
};



class irpShiftRemoveNeighborhood : public emili::Neighborhood
{
protected:
    InventoryRoutingProblem irp;
    Solution * currentNeighboringSolution;
    int shift;
    int operation;
    int numberOfOperations;
    int numberOfShifts;
    double randomFactor;
    unsigned int urgencyPolicy;
    double refuelFactor;
    double deliveredQuantityFactor;
    int numberFeasibleSolutions; //in a neighborhood
    double bestValueFound;
    virtual Solution* computeStep(Solution* step);
    virtual void reverseLastMove(Solution* step);
public:

     /*
      * Initialize shift and operation size, 4 indexes
      * this method needs to be overidden if there are things that a neighborhood has to do or reset when begin is called*/
     virtual NeighborhoodIterator begin(emili::Solution* base);

    /*this method returns a solution in the decided neighborhood
     * of the currentSolution */
    virtual Solution* step(Solution* currentSolution);
    /*
     * The state of the Neighborhood object may need to be restored to
     * initial conditions between local search calls
     * (e.g. first improvement strategies for permutation flow shop).
     */
    virtual void reset();
    /*
     * A method that returns a random solution in the neighborhood has to be provided
     */
    virtual Solution* random(Solution* currentSolution);
    /*
     * This method returns the size of the neighborhood
    */
    virtual int size();

    irpShiftRemoveNeighborhood(InventoryRoutingProblem &i, double rf, unsigned int up, double ref, double dqf):
        irp(i),numberFeasibleSolutions(0),bestValueFound(DBL_MAX),
        randomFactor(rf),urgencyPolicy(up),refuelFactor(ref), deliveredQuantityFactor(dqf){}
};



class irpRefuelNeighborhood : public emili::Neighborhood
{
protected:
    InventoryRoutingProblem irp;
    Solution * currentNeighboringSolution;
    int shift;
    int numberOfShifts;
    double refuelRatio;
    double deliveredQuantityRatio;
    double refuelStep;
    double deliveredQuantityStep;
    double randomFactor;
    unsigned int urgencyPolicy;
    int numberFeasibleSolutions; //in a neighborhood
    double bestValueFound;
    virtual Solution* computeStep(Solution* step);
    virtual void reverseLastMove(Solution* step);

public:
     virtual NeighborhoodIterator begin(emili::Solution* base);

    /*this method returns a solution in the decided neighborhood
     * of the currentSolution */
    virtual Solution* step(Solution* currentSolution);
    /*
     * The state of the Neighborhood object may need to be restored to
     * initial conditions between local search calls
     * (e.g. first improvement strategies for permutation flow shop).
     */
    virtual void reset();
    /*
     * A method that returns a random solution in the neighborhood has to be provided
     */
    virtual Solution* random(Solution* currentSolution);
    /*
     * This method returns the size of the neighborhood
    */
    virtual int size();

    irpRefuelNeighborhood(InventoryRoutingProblem &i, double rs, double dqs, double rf, unsigned int up):irp(i),numberFeasibleSolutions(0),bestValueFound(DBL_MAX),
        refuelStep(rs),deliveredQuantityStep(dqs),randomFactor(rf), urgencyPolicy(up){}
};

class irpPerturbation : public emili::Perturbation
{
public:
  virtual Solution* perturb(Solution* solution);

};

}
}
#endif // IRP_H
