#include "../emilibase.h"
namespace emili {
namespace pop {


/**
 * @brief The MBO class
 * Migrating Birds Optimization algorithm of Sioud Gagné 2017
 * that for some reasons ( Unknown to me ) is so
 * popular...
 */
class EMBO : public emili::LocalSearch
{
protected:    
    emili::Solution* leader;
    emili::LocalSearch* ls;
    emili::Perturbation* pert;
    std::vector< emili::Solution* > leftwing;
    std::vector< emili::Solution* > rightwing;
    int pop_size;
    int k; // k solutions evaluated for the leader
    int x; // k-x solutions evaluated for the other solutions
    int m;
    int age;
    std::vector< int > leftwing_age;
    std::vector< int > rightwing_age;
    int leader_age;
    float q0;
public:
    EMBO(InitialSolution& initialSolutionGenerator ,
        Termination& terminationcriterion,
        Neighborhood& neighborh,
        Perturbation* pert_,
        LocalSearch* ls_,
        int popsize, int kp, int xp, int mp, int _age, float _q0):
        emili::LocalSearch(initialSolutionGenerator,terminationcriterion,neighborh),
        pert(pert_),
        ls(ls_),
        pop_size(popsize),
        k(kp),
        x(xp),
        m(mp),
        age(_age),
        leftwing_age(popsize,0),
        rightwing_age(popsize,0),
        leader_age(0),
        q0(_q0)
    { }

    virtual Solution* search(Solution* initial);
    virtual Solution* getBestSoFar();
};

}
}
