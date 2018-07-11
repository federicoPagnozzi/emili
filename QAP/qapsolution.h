#ifndef QAPSOLUTION_H
#define QAPSOLUTION_H


#include <string>
#include <vector>

#include "../emilibase.h"

namespace emili {
namespace qap{

class QAPSolution: public emili::Solution {

protected:
    //double solutionValue;
    std::vector< int > solution;

    virtual const void* getRawData()const;
    virtual void setRawData(const void* data);

public:
    QAPSolution(std::vector< int >& solution):
        emili::Solution(1e19),
        solution(solution) { }

    virtual double getSolutionValue(void) {
        return solution_value;
    }

    std::vector<int>& getSolution(void) {
        return solution;
    }

    void setSolutionValue(double value) {
        solution_value = value;
    }

    /**
     * does not set solutionValue too
     * @param _solution 
     */
    void setSolution(std::vector< int > _solution) {
        solution = _solution;
    }
   
    virtual std::string getSolutionRepresentation(void);
    
    
    virtual emili::Solution* clone() {
        emili::Solution* newsol = new QAPSolution(solution);
        newsol->setSolutionValue(solution_value);
        return newsol;
    }

}; // QAPSolution


}
}
#endif
