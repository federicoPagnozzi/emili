#include "maneh.h"
bool check_probability(float probability)
{
    return emili::generateRealRandomNumber() < probability ;
}

emili::Solution* tournament_selection(std::vector< emili::Solution* > pop, int size)
{
    int psize = pop.size();
    int first = emili::generateRandomNumber()%psize;
    emili::Solution* selected = pop[first];

    for(int i=1; i < size ; i++)
    {
        int job = emili::generateRandomNumber()%psize;
        if(*selected > *pop[job])
            selected = pop[job];
    }

    return selected;
}

void erase_job(std::vector< int >& permutation, int job)
{

    for(int i=0; i< permutation.size();i++)
    {
        if(permutation[i] == job)
        {
            permutation.erase(permutation.begin()+i);
            break;
        }
    }

}

void random_permutation(std::vector<int>& ranper, int n)
{
        std::vector<int> already_taken(n+1,0);
        already_taken[0] = 1;
        int j = 1;
        while(j<=n)
        {
            int job = (rand()%n)+1;
            if(already_taken[job] == 0)
            {
            ranper[j] = job;
            already_taken[job] = 1;
            j++;
            }
        }
}

void calc_nhm(std::vector< emili::pfsp::PermutationFlowShopSolution* >& spop, std::vector< std::vector< int > >& nhm, int n, int pop)
{
    for(int i=1;i<=n ;i++)
    for(int j=1;j<=n ;j++)
       for(int k=0; k<pop; k++)
       {
           emili::pfsp::PermutationFlowShopSolution& kpop = *spop[k];
        if(kpop[i] == j)
        {
           nhm[j][i]++;
        }
       }
}

void calc_ehm(std::vector< emili::pfsp::PermutationFlowShopSolution* >& spop, std::vector< std::vector< int > >& ehm, int n, int pop)
{
    for(int i=1;i<=n ;i++)
    for(int j=1;j<=n ;j++)
       for(int k=0; k<pop; k++)
       {
          emili::pfsp::PermutationFlowShopSolution& kpop = *spop[k];
        for(int h=1; h<=n ; h++)
        {
            int hp = h+1;
            if(h==n)
            {
                hp = 1;
            }
            if(kpop[h] == i && kpop[hp] == j)
            {
                ehm[j][i]++;
            }
        }
       }
}

void calc_nehm(std::vector< std::vector< float > >& nehm,std::vector< std::vector< int > >& ehm,std::vector< std::vector< int > >& nhm, int poposize, float Bratio, int n )
{
    for(int i=1;i<=n ;i++)
        for(int j=1;j<=n ;j++)
        {
            nehm[j][i] = nhm[j][i] + ehm[j][i] + (1/(float)n + 1/(float)(n-1))*poposize*Bratio;
        }
}


void emili::pop::MANEH::compute_nhem(std::vector< std::vector< float >>& nehm, std::vector< emili::Solution* >& spop)
{
    emili::pfsp::PermutationFlowShopSolution** ps = (emili::pfsp::PermutationFlowShopSolution**)spop.data();
    std::vector< emili::pfsp::PermutationFlowShopSolution* > stuff(ps,ps+sPS);
    std::vector< std::vector< int > > ehm(njobs+1,std::vector< int >(njobs+1,0));
    std::vector< std::vector< int > > nhm(njobs+1,std::vector< int >(njobs+1,0));
    calc_ehm(stuff,ehm,njobs,sPS);
    calc_nhm(stuff,nhm,njobs,sPS);
    calc_nehm(nehm,ehm,nhm,PS,Bratio,njobs);
}


void random_sample_crossover(std::vector< std::vector< float > >& nehm, std::vector< int>& templt, std::vector< int>& res,int n)
{
   // std::cout << "\n RSC " << std::endl;
    std::vector<int> ranper(n+1,0);
    random_permutation(ranper, n);

    int p=1;
    int beforejob = (rand()%n)+1;
    int rsize = n;
    while(p <= n)
    {
       int flag = emili::generateRandomNumber();
       if(flag)
          {
        res[p] = templt[0];
        templt.erase(templt.begin());
        erase_job(ranper, res[p]);
          }
       else
          {
        //roulette wheel from NEHM and ranper
        int pir = ranper.size();
        float nejp = 0.0;
        std::vector< float > weights(n+1,1.0);
        for(int i = 1; i < pir; i++ )
        {
            nejp += nehm[ranper[i]][p];
        }
        float pw = 0.0;
        for(int i=1 ; i< pir; i++)
        {
            pw += nehm[ranper[i]][p] / nejp;
            weights[ranper[i]] = pw;
        }
        //sample node x with probability ....
        float prob = emili::generateRandomNumber();
        int x = ranper[0];
        for(int i=1; i<pir ; i++)
        {
            if(prob <= weights[i])
            {
                x = ranper[i];
                break;
            }
        }
        // add x to res
        res[p] = x;
        // delete x from templt and ranper
        erase_job(templt, x);
        erase_job(ranper, x);
          }
       beforejob = res[p];
       p++;
    }
   /* for(int i=1; i<=n; i++)
             std::cout << " " << res[i];
        std::cout << std::endl;
*/
}

emili::Solution* emili::pop::MANEH::RSC(emili::Solution* tmpl, std::vector< std::vector< float > >& nehm)
{
    emili::Solution* res = nullptr;
    emili::pfsp::PermutationFlowShopSolution* t = (emili::pfsp::PermutationFlowShopSolution*)tmpl;
    if(check_probability(Pc))
    {
        std::vector<int> templ(t->getJobSchedule());
        templ.erase(templ.begin());
        std::vector<int> vres(njobs+1,0);
        random_sample_crossover(nehm,templ,vres,njobs);
        res = new emili::pfsp::PermutationFlowShopSolution(vres);
        pfs.evaluateSolution(*res);
        return res;
    }
    else
    {
        return tmpl->clone();
    }
}

void emili::pop::MANEH::check_and_replace(emili::Solution* s)
{
    if(*s < *pop[PS-1])
    {
        bool unique = true;
        for(auto iter=pop.begin(); iter!=pop.end(); ++iter)
        {
            if((*iter)->operator<(*s))
            {
                break;
            }

            if((*iter)->getSolutionValue() == s->getSolutionValue())
            {
                unique = false;
            }
        }
        if(unique)
        {
            emili::Solution* t = pop[PS-1];
            pop[PS-1] = s;
            delete t;
            int idx = PS-1;
            for(int i=PS-2;i>=0;i--)
            {
                if(*pop[i] > *s)
                {
                    idx = i;
                }
                else
                {
                    break;
                }
            }
            if(idx < PS-1)
            {
                std::swap(pop[idx],pop[PS-1]);
            }
            if(idx == 0)
            {
                if(*bestSoFar > *pop[0])
                {
                    *bestSoFar = *pop[0];
                }
            }
            return;
        }
    }
    delete s;
}


emili::Solution* emili::pop::MANEH::search(emili::Solution* initial)
{
    //Initialize pop
    *bestSoFar = *initial;
    pop.push_back(initial->clone());
    for(int i = 0; i < PS ; i++)
    {
        pop.push_back(init->generateSolution());
    }
    // order pop by ms
    std::sort(pop.begin(),pop.end(),[](emili::Solution* i1,emili::Solution* i2){return *i1 < *i2;});
    int popcount = 0;
    // Selects alpha*Ps individuals with truncation selection and build nehm
    std::vector< emili::Solution* > alphaps(pop.begin(),pop.begin()+sPS);
    std::vector< std::vector< float > > nehm(njobs+1,std::vector< float >(njobs+1,0.0));
    do{
        emili::Solution* parent1 = tournament_selection(pop,2);
        emili::Solution* parent2 = tournament_selection(pop,2);

        emili::Solution* child1 = RSC(parent1,nehm);
        emili::Solution* child2 = RSC(parent2,nehm);

        // Child1 = mutate(Child1,mc)
        if(check_probability(Pc))
        {
            emili::Solution* m = mutation->perturb(child1);
            delete child1;
            child1 = m;
        }
        // Child2 = mutate(Child2,mc)
        if(check_probability(Pc))
        {
            emili::Solution* m = mutation->perturb(child2);
            delete child2;
            child2 = m;
        }
        // Child1 = ls1(Child1,ls)
        if(check_probability(ls))
        {
            emili::Solution* m = ls1->search(child1);
            delete child1;
            child1 = m;
        }
        // Child2 = ls1(Child2,ls)
        if(check_probability(ls))
        {
            emili::Solution* m = ls1->search(child2);
            delete child2;
            child2 = m;
        }

        // if child1/2 < worst or child1/2 is unique
        // replace worst with child1/2
        check_and_replace(child1);
        check_and_replace(child2);
        popcount +=2;
        if(popcount > (PS/2))
        {
            //Update best
            emili::Solution* ns = ls2->search(pop[0]);
            if(*ns < *pop[0])
            {
                emili::Solution* t = pop[0];
                pop[0] = ns;
                delete t;
            }
            else
            {
                delete ns;
            }
            // if best < bestsofar
            // update bestsofar
            if(*pop[0] < *bestSoFar)
            {
                *bestSoFar = *pop[0];
            }
            // Selects alpha*Ps individuals with truncation selection and build nehm
            std::vector< emili::Solution* > spops(pop.begin(),pop.begin()+sPS);
            compute_nhem(nehm,spops);
            popcount = 0;
        }


    }
    while(!termcriterion->terminate(pop[0],pop[1]));
    return bestSoFar;
}
