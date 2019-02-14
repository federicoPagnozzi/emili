#include "hdtlm.h"
#include <algorithm>

void job_at_index(std::vector<int>& permutation, std::vector<int>& indexes, int njobs)
{
        for(int i=1; i<=njobs; i++ )
                indexes[permutation[i]] = i;
}

void calc_Imat(std::vector< emili::pfsp::PermutationFlowShopSolution* > & stuf, std::vector< std::vector < int >>& Im, int njobs, int eps)
{
        for(int i = 1; i <= njobs; i++)
                {

                        for(int j=1; j<= njobs; j++)
                        {
                                for(int k=0; k < eps ; k++)
                                {
                                        int appear = 0;
                                        for(int l = 1; l<=j; l++)
                                                if((*stuf[k])[l] == i )
                                                {
                                                   appear = 1;
                                                   break;
                                                }
                                        Im[i][j] += appear;

                                }
                        }
                }
}

void calc_pmat(std::vector < std::vector< float >>& p_m,std::vector< std::vector < int >> & Im,std::vector< std::vector < int >> & Imb,float alpha, float beta, int njobs, int eps)
{
        std::vector<float> p_col(njobs+1,0.0);
        for(int i = 1; i<=njobs; i++ )
        {
                for(int j = 1; j<= njobs; j++)
                {
                        //strange formula
                        float Iij = Im[i][j];
                        p_m[i][j] = (1-beta) * ((1-alpha)*p_m[i][j]+ (Iij * (alpha/(j*eps)))  ) + (beta/j)*Imb[i][j];
                        p_col[i] += p_m[i][j];
                }
        }


        for(int i = 1; i<=njobs; i++ )
        {
                for(int j = 1; j<= njobs; j++)
                {
                        p_m[i][j] = p_m[i][j]/p_col[i];
                }
        }
}

void calc_consensu_permutation(std::vector< std::vector < int >>& pop, int ps, std::vector< int >& pcon, int njobs )
{
        std::vector< float > sigma(njobs+1, 0);
        for(int i = 1; i <= njobs; i++)
        {
                for(int k=0; k<ps; k++)
                {
                        sigma[i] += (pop[k][i]-1);
                }
                sigma[i] = sigma[i]/(float)ps;
        }

        std::vector<int> visited(njobs+1,0);
        for(int j = 1; j <= njobs ; j++)
        {
                float lmin=ps;
                int l = 1;
                for(int k=1; k<= njobs;k++)
                {
                        if(!visited[k] && sigma[k] < lmin)
                        {
                                l = k;
                                lmin = sigma[k];
                        }
                }
                //std::cout << "| p[" << l << "] = " << j << " |" << sigma[l] << std::endl;

                pcon[l] = j;
                visited[l] = 1;
        }

}

void sample_permutation( std::vector< std::vector < float >>& pm, std::vector<int>& pcon, std::vector<int>& newp, int njobs)
{
        std::vector<int> vp(njobs+1,0);
        for(int i = 1; i <= njobs; i++ )
        {
                int pos = 1;
                float prob = emili::generateRealRandomNumber();//fakeRand(i);
                float cas = 0.0;

                for(int j=1; j <= njobs; j++)
                {
                        cas += pm[i][j];
                        if(prob < cas)
                        {
                                pos = j;
                                break;
                        }
                }
                vp[pcon[i]] = pos;
        }

        for(int i = 1; i<= njobs; i++)
        {
                int lp = vp[i] ;
                if(newp[lp] == 0)
                {
                        newp[lp] = pcon[i];
                }
                else
                {
                        int h = lp+1 ;
                        for(;h<=njobs;h++)
                        {
                                if(newp[h]==0)
                                   break;
                        }
                        int d = lp-1;
                        for(;d>0;d--)
                        {
                                if(newp[d]==0)
                                   break;
                        }

                        if((h-lp) < (lp-d))
                        {
                                newp[h] = pcon[i];
                        }
                        else
                        {
                                newp[d] = pcon[i];
                        }
                }
        }
}

void pmx(std::vector<int>& parent1, std::vector<int>& parent2, std::vector<int>& child, int njobs)
{
        std::vector<int> j2i(njobs+1,0);
        job_at_index(parent2, j2i, njobs);

        int lb = emili::generateRandomNumber()%(njobs);
        int up = emili::generateRandomNumber()%(njobs);
        if(lb>up)
        {
                std::swap(lb,up);
        }
        std::vector<int> to_insert;
        std::vector<int> to_check(njobs+1,1);

        for(int i=lb; i <= up;i++)
        {

                child[i] = parent1[i];
                to_check[parent1[i]] = 0;
        }

        for(int i=lb; i<=up; i++)
        {
                if(to_check[parent2[i]])
                {
                        //to_insert.push_back(parent2[i]);
                        int idx = j2i[parent1[i]];
                        while(lb<=idx && idx<=up)
                        {
                                idx = j2i[parent1[idx]];
                        }
                        child[idx] = parent2[i];
                }
        }

        for(int i=1; i<=njobs; i++ )
        {
                if(child[i] != 0)
                     j2i[child[i]] = 0;
        }

        int k = 1;
        for(int i=1; i<= njobs; i++)
        {
                if(j2i[parent2[i]] != 0)
                {
                        while(child[k] != 0)
                                k++;

                        child[k] = parent2[i];
                }
        }
}

emili::Solution* emili::pop::InsertPathRelink::perturb(emili::Solution* s1, emili::Solution* s2)
{
    emili::pfsp::PermutationFlowShopSolution* ss1= (emili::pfsp::PermutationFlowShopSolution*) s1;
    emili::pfsp::PermutationFlowShopSolution* ss2= (emili::pfsp::PermutationFlowShopSolution*) s2;

    std::vector<int>& spoint = ss1->getJobSchedule();
    std::vector<int>& epoint = ss2->getJobSchedule();
    std::vector<int> best(njobs+1,0);
    int best_value = s1->getSolutionValue();

    int i = 1;
    int j = 1;
    std::vector<int> temp(spoint);
    while(i != njobs)
    {
        if(epoint[i] != spoint[i])
        {
            int job = epoint[i];
            for(j = i; j<=njobs; j++)
            {
                if(temp[j] == job)
                    break;
            }
            temp.erase(temp.begin()+j);
            temp.insert(temp.begin()+i,job);
            int value = prob.computeObjectiveFunction(temp);
            if(value < best_value)
            {
                best = temp;
                best_value = value;
            }
        }
        i++;
    }

    if(best_value == s1->getSolutionValue())
    {
        return p->perturb(s1);
    }
    else
    {
        emili::Solution* ret = new emili::pfsp::PermutationFlowShopSolution(best_value,best);
        return ret;
    }
}

void emili::pop::HDTLM::check_and_replace(emili::Solution* s)
{
    if(*s < *pop[PS-1])
    {
        bool unique = true;
        for(auto iter=pop.begin(); iter!=pop.end(); ++iter)
        {
            if((*iter)->getSolutionValue() == s->getSolutionValue())
            {
                unique = false;
                break;
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
                pop.erase(pop.begin()+PS-1);
                pop.insert(pop.begin()+idx,s);
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

void emili::pop::HDTLM::procedure4()
{
    //Calculate consensus ranking permutation
    std::vector < std::vector < float > > p_m(njobs+1, std::vector< float>(njobs+1, 1.0/njobs) );
    std::vector < std::vector< int > > Im (njobs+1, std::vector< int > (njobs+1,0));
    std::vector < std::vector< int > > Imb (njobs+1, std::vector< int > (njobs+1,0));
    emili::pfsp::PermutationFlowShopSolution** ps = (emili::pfsp::PermutationFlowShopSolution**)pop.data();
    std::vector< emili::pfsp::PermutationFlowShopSolution* > stuff(ps,ps+PS);
    //Order solutions ?
    // Select lamba * PS
    // Calculate the probabilistic model ( procedure 1)
    calc_Imat(stuff,Im,njobs,ePS);
    calc_Imat(stuff,Imb,njobs,1);
    calc_pmat(p_m,Im,Imb,alpha,beta,njobs,ePS);
    // for each learner PI_i
    for(int i = 0; i < PS ; i++)
    {
    //   PI = sample P (procedure 2)
        std::vector< int>& PI_I = stuff[i]->getJobSchedule();
        std::vector<int > PI(njobs+1,0);
        sample_permutation(p_m,PI_I,PI,njobs);
    //   PI_p = pmx(PI, PI_i)
        std::vector<int> PI_p(njobs+1,0);
        pmx(PI,PI_I,PI_p,njobs);
        emili::Solution* newp = new emili::pfsp::PermutationFlowShopSolution(pfs.computeObjectiveFunction(PI_p),PI_p);
    //   if rand < ls p3(pi_p)
        double r= emili::generateRealRandomNumber();
        if(r < ls)
        {
            emili::Solution* t = p3->search(newp);
            if(*t < *newp)
            {
                emili::Solution* toDel = newp;
                newp = t;
                delete toDel;
            }
            else
            {
                delete t;
            }
        }
    //   if pi_p < pi_i and pi_p is unique -> pi_i = pi_p
        if(*newp < *stuff[i])
            check_and_replace(newp);
    }
    // update best
    if(*bestSoFar > *pop[0])
    {
        *bestSoFar = *pop[0];
    }
}

void emili::pop::HDTLM::procedure7()
{
    // Super Elite pop
    // perturb pop[0]
    emili::Solution* pert = per->perturb(pop[0]);
    // nel = ris(pop[0])
    emili::Solution* rissed = p5->search(pert);
    delete pert;
    // if nel < pop[0] -> pop[0] = nel
    if(*rissed < *pop[0])
    {
        emili::Solution* s = pop[0];
        pop[0] = rissed;
        delete s;
    }
    else
    {
        delete rissed;
    }
    // middle layer
    // for each elite pop ep_i (pop[1,ePS])
    //   pi_new = PMX(ep_i,pop[0])
    //   if pi_new < ep_i -> ep_i = pi_new

    // bottom layer
    // for each of the remaining opi
    // rep = random ep_i
    // pi_new = path_relink_shift(opi,epi)
    // if pi_new < opi -> opi = pi_new

    // update best
}

emili::Solution* emili::pop::HDTLM::search(emili::Solution* initial)
{
    //Initialize pop
    *bestSoFar = *initial;
    int pophalf = PS/2;
    pop.push_back(initial->clone());
    for(int i=1; i<pophalf;i++)
        pop.push_back(init->generateSolution());
    for(int i=0; i<pophalf;i++)
        pop.push_back(popinit->generateSolution());

    //Order solutions and select best as teacher ( the teacher is the bestsofar)
    std::sort(pop.begin(),pop.end(),[](emili::Solution* i1, emili::Solution* i2){return *i1 < *i2;});
    if(*pop[0] < *bestSoFar)
    {
        *bestSoFar = *pop[0];
    }

    do
    {
    //while(!termination)
    //procedure 4
     procedure4();
    //procedure 7
     procedure7();
    emili::Solution* pnbest = p8->search(pop[0]);
    if(*pnbest < *pop[0])
    {
        emili::Solution* t = pop[0];
        pop[0] = pnbest;
        delete t;
    }
    else
    {
        delete pnbest;
    }
    if(*pop[0] < *bestSoFar)
    {
        *bestSoFar = *pop[0];
    }
    }while(!termcriterion->terminate(pop[0],pop[1]));

    return bestSoFar;
}
