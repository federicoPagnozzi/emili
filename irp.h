#ifndef IRP_H
#include "emilibase.h"
#define IRP_H

namespace emili
{
namespace irp
{

class InvetoryRoutingProblem: public emili::Problem
{
protected:
    /*Internal structure of the invetory routing problem*/
public:
    InvetoryRoutingProblem(char* instance_path);
    /*implement this method to get the objective Function value*/
    virtual double evaluateSolution(Solution & solution);

    /*Methods that are needed  for IRP */
};

class InvetoryRoutingSolution: public emili::Solution
{
protected:
    /*Implement these methods to make the copy constructor works*/
    virtual const void* getRawData()const=0;
    virtual void setRawData(const void* data)=0;

    /*Internal structure of the Solution*/

public:
    /*Base constructor*/
    InvetoryRoutingSolution(double solution_value):emili::Solution(solution_value) { }
    /*This method must build a clone of the Solution on the heap and returns the pointer to it*/
    virtual Solution* clone();

    /*Methods needed for the Internal use of the Solution
     * ( how the solution is seen by Neighborhoods, Perturbations, Initial Solution generators... etc)*/

};




}
}
#endif // IRP_H
