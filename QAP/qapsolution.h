#ifndef QAPSOLUTION_H
#define QAPSOLUTION_H


#include <string>
#include <vector>

#include "../emilibase.h"


class QAPSolution: public emili::Solution {

protected:
    double solutionValue;
    std::vector< int > solution;

    virtual const void* getRawData()const;
    virtual void setRawData(const void* data);

public:
    QAPSolution(std::vector< int >& solution):
        emili::Solution(1e9),
        solutionValue(0.0),
        solution(solution) { }

    virtual double getSolutionValue(void) {
        return solutionValue;
    }

    std::vector<int> getSolution(void) {
        return solution;
    }

    void setSolutionValue(double value) {
        solutionValue = value;
    }

    /**
     * does not set solutionValue too
     * @param _solution 
     */
    void setSolution(std::vector< int > _solution) {
        solution = _solution;
    }
   
    virtual std::string getSolutionRepresentation(void);

}; // QAPSolution


#endif
