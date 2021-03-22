#include "tspsolution.h"

using namespace emili::tsp;
const void* TSPSolution::getRawData()const {
    return (void*)&solution;
}

void TSPSolution::setRawData(const void* data) {
    std::vector< int >* sol = (std::vector< int >*) data;
    solution = *sol;
}



std::string TSPSolution::getSolutionRepresentation(void) {

    std::string solrep = "[ ";

    for (int i = 0; i < this->solution.size(); i++) {
        solrep.append(std::to_string(this->solution[i]));
        solrep.append(" ");
    }

    solrep.append("]");
    
    return solrep;

}

bool TSPSolution::operator<(TSPSolution& a)
{

    return solution_value < a.solution_value;
}

bool TSPSolution::operator<=(TSPSolution& a)
{
    return solution_value <= a.solution_value;
}

bool TSPSolution::operator>=(TSPSolution& a)
{
    return solution_value >= a.solution_value;
}

bool TSPSolution::operator>(TSPSolution& a)
{
    return solution_value > a.solution_value;
}

bool TSPSolution::operator==(TSPSolution& a)
{
    return solution_value == a.solution_value;
}

