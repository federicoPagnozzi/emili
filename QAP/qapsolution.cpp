#include "qapsolution.h"


const void* QAPSolution::getRawData()const {
    return (void*)&solution;
}

void QAPSolution::setRawData(const void* data) {
    std::vector< int >* sol = (std::vector< int >*) data;
    solution = *sol;
}



std::string QAPSolution::getSolutionRepresentation(void) {

    std::string solrep = "\n[ ";

    for (int i = 0; i < this->solution.size(); i++) {
        solrep.append(std::to_string(this->solution[i]));
        solrep.append(" ");
    }

    solrep.append("]\n");
    
    return solrep;

}