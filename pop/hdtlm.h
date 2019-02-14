#ifndef HDTLM_H
#define HDTLM_H
#include "../pfsp/permutationflowshop.h"

namespace emili{
namespace pop {

class InsertPathRelink : public emili::TwoSolutionPerturbation
{
protected:
    emili::Perturbation* p;
    emili::pfsp::PermutationFlowShop& prob;
    int njobs;
public:
    InsertPathRelink(emili::Perturbation* pr, emili::pfsp::PermutationFlowShop& problem):
        emili::TwoSolutionPerturbation(),p(pr),prob(problem),njobs(problem.getNjobs()){ }
    virtual emili::Solution* perturb(emili::Solution* s1, emili::Solution* s2);
};

class HDTLM : public emili::LocalSearch
{
protected:
    std::vector<emili::Solution* > pop;
    int PS;
    int ePS;
    int njobs;
    emili::Perturbation* per;
    emili::pop::InsertPathRelink* isp;
    emili::InitialSolution* popinit; //second initial solution
    emili::pfsp::PermutationFlowShop& pfs;
    emili::LocalSearch* p3; //insert local search
    emili::LocalSearch* p5; // RIS
    emili::LocalSearch* p8; // ILS
    float lambda;
    float alpha;
    float beta;
    float ls;
    void check_and_replace(emili::Solution* s);
    void procedure4();
    void procedure7();
public:
    HDTLM(Termination& termination, InitialSolution& popinit1, InitialSolution& popinit2,
          emili::Perturbation& perturbation, emili::LocalSearch* ls1, emili::LocalSearch* ls2,
          emili::LocalSearch* ls3, int popsize, float _lambda, float _alpha, float _beta, float _ls):
          emili::LocalSearch(popinit1,termination,*(new emili::EmptyNeighBorHood())),
          pfs((emili::pfsp::PermutationFlowShop&)popinit1.getProblem()), popinit(&popinit2),
          per(&perturbation),p3(ls1),p5(ls2),p8(ls3),lambda(_lambda),alpha(_alpha),beta(_beta),
          ls(_ls),isp(new emili::pop::InsertPathRelink(&perturbation,(emili::pfsp::PermutationFlowShop&)popinit1.getProblem()))
    {
        njobs = pfs.getNjobs();
    }

    virtual emili::Solution* search(Solution* initial);
};


}
}
#endif // HDTLM_H
