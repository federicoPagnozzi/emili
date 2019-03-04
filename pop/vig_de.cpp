#include "vig_de.h"
#include <cmath>

void emili::pop::vIG_DE::mutatePopulation(int& a, int& b, int& c,int i)
{
    int a1 = emili::generateRandomNumber()%pop_size;
    int b1 = emili::generateRandomNumber()%pop_size;
    int c1 = emili::generateRandomNumber()%pop_size;

    while(a==i)
        a = emili::generateRandomNumber()%pop_size;
    while(a1==i && a1==a)
        a1 = emili::generateRandomNumber()%pop_size;
    if(*pop[a1] < *pop[a])
        a = a1;

    while(b==i && b==a)
        b = emili::generateRandomNumber()%pop_size;
    while(b1==i && b1==a && b1==b)
        b1 = emili::generateRandomNumber()%pop_size;
    if(*pop[b1] < *pop[b])
        b = b1;

    while(c==i && c==a && c==b)
        c = emili::generateRandomNumber()%pop_size;
    while(c1==i && c1==a && c1==b && c1==c)
        c1 = emili::generateRandomNumber()%pop_size;
    if(*pop[c1] < *pop[c])
        c = c1;

}

emili::Solution* emili::pop::vIG_DE:: search(Solution* initial)
{
    // Initialize pop
    int njobs = pfs->getNjobs();
    xij[0].first = emili::generateRandomNumber()%njobs+1;
    emili::Solution* temp = init->generateSolution();
    igp->setDestruction(xij[0].first);
    pop[0] = igp->perturb(temp);
    *bestSoFar = *pop[0];
    delete temp;
    long int sumfpi = pop[0]->getSolutionValue();
    for( int i = 1; i < pop_size; i++)
    {
        xij[i].first = emili::generateRandomNumber()%njobs+1;
        temp = secondary->generateSolution();
        igp->setDestruction(xij[0].first);
        pop[i] = igp->perturb(temp);
        delete temp;
        sumfpi += pop[i]->getSolutionValue();
        if(*pop[i] < *bestSoFar)
        {
            *bestSoFar = *pop[i];
        }
    }
    xij[0].second = 1 - (pop[0]->getSolutionValue()/sumfpi);    
    for( int i = 1; i < pop_size; i++)
    {
        xij[i].second = 1 - (pop[i]->getSolutionValue()/sumfpi);
    }
    while(!termcriterion->terminate(pop[0],pop[1]))
    {
        for(int i=0; i<pop_size; i++)
        {
            emili::iteration_increment();
            // mutate
            int a = emili::generateRandomNumber()%pop_size;
            int b = emili::generateRandomNumber()%pop_size;
            int c = emili::generateRandomNumber()%pop_size;
            mutatePopulation(a,b,c,i);            
            // Create Mutation
            int V_d = xij[a].first + Fr * std::abs(xij[b].first - xij[c].first);
            float V_p = xij[a].second + Fr * (xij[b].second - xij[c].second);       
            // Create trial
            int U_d = Cr*V_d + (1-Cr)*xij[i].first;
            float U_p = Cr*V_p + (1-Cr)*xij[i].second;
            // Verify that trial is not out of range
            if(U_d >= njobs)
            {
                U_d = 1+(njobs-1)*emili::generateRealRandomNumber();
            }
            if(U_p >= 1)
            {
                U_p = emili::generateRealRandomNumber();
            }
            // evaluate trial
            emili::Solution* trial=nullptr;
            float thr = emili::generateRealRandomNumber();
            if(thr < U_p)
            {                
                igp->setDestruction(U_d);
                emili::Solution* temp = igp->perturb(pop[i]);                
                trial = ris->search(temp);                
                delete temp;
            //}

            // update pop
            //if(trial != nullptr)
            //{
                if(*trial < *pop[i])
                {
                    //delete pop[i];
                    *pop[i] = *trial;
                    if(*pop[i] < *bestSoFar)
                    {
                        //update best so far
                        *bestSoFar = *pop[i];
                    }
                }
                delete trial;
            }

        }
    }


    return pop[0];
}

emili::Solution* emili::pop::vIG_DE::getBestSoFar()
{
    emili::Solution* b = ris->getBestSoFar();
    if(*b < *bestSoFar)
    {
        return b;
    }
    return bestSoFar;
}

