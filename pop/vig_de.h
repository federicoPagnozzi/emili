#ifndef VIG_DE_H
#define VIG_DE_H
#include "../pfsp/permutationflowshop.h"

namespace emili {
namespace pop {

class vIG_DE : public emili::LocalSearch
{
protected:

std::vector< emili::Solution* > pop;
std::vector< std::pair<int,float> > xij;

emili::InitialSolution* secondary;
emili::pfsp::IGPerturbation* igp;
emili::pfsp::RIS* ris;
emili::pfsp::PermutationFlowShop* pfs;
int pop_size;
float Fr;
float Cr;

void mutatePopulation(int& a, int& b, int& c,int i);


public:
vIG_DE(Termination& terminationcriterion,InitialSolution& first_pop,
       InitialSolution& pop_initializer, int popsize, float mscale_factor, float cross_factor):
       emili::LocalSearch(first_pop,terminationcriterion,*(new emili::EmptyNeighBorHood())),
       pop_size(popsize),Fr(mscale_factor),Cr(cross_factor),secondary(&pop_initializer),
       xij(std::vector< std::pair<int,float> >(popsize,std::pair<int,float>(0,0.0f))),
       pop(popsize,nullptr)
    {
        pfs = (emili::pfsp::PermutationFlowShop*) &init->getProblem();
        igp = new emili::pfsp::IGPerturbation(4,*pfs);
        ris = new emili::pfsp::RIS(*pfs,*init);
    }

virtual emili::Solution* search(Solution* initial);
virtual emili::Solution* getBestSoFar();
};

}
}
#endif // VIG_DE_H
