#ifndef MANEH_H
#define MANEH_H

#include "../pfsp/permutationflowshop.h"
namespace emili {
namespace pop{

class MANEH : public emili::LocalSearch
{
protected:
    std::vector< emili::Solution* > pop;
    emili::LocalSearch* ls1;
    emili::LocalSearch* ls2;
    emili::Perturbation* mutation;
    int PS; // population size
    int sPS;
    int njobs;
    float Pm; //prob mutation
    float Pc; //prob crossover
    float Bratio; // bratio ?? 0.2
    float ls; // prob localsearch
    float alpha; // PS*alpha gives the sPS
    emili::pfsp::PermutationFlowShop& pfs;
    void check_and_replace(emili::Solution* s);
    void compute_nhem(std::vector< std::vector< float >>& nehm, std::vector< emili::Solution* >& spop);
    emili::Solution* RSC(emili::Solution* tmpl, std::vector< std::vector< float > >& nehm);
public:
    MANEH(Termination& termination, InitialSolution& popinit,
          int popsize, emili::LocalSearch* firstLS, emili::LocalSearch* secondLS,
          emili::Perturbation* mut, float Pmp, float Pcp, float Bratiop,
          float lsp,float alphap):
        emili::LocalSearch(popinit,termination,*(new emili::EmptyNeighBorHood())),
        ls1(firstLS),ls2(secondLS),mutation(mut),PS(popsize),Pm(Pmp),Pc(Pcp),Bratio(Bratiop),
        ls(lsp),alpha(alphap),pfs((emili::pfsp::PermutationFlowShop&)popinit.getProblem()),
        sPS(alphap*popsize)
    {
        njobs = pfs.getNjobs();
    }

    virtual emili::Solution* search(Solution* initial);
};

}
}
#endif // MANEH_H
