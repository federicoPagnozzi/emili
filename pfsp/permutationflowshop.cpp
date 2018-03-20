//
//  Created by Federico Pagnozzi on 28/11/14.
//  Copyright (c) 2014 Federico Pagnozzi. All rights reserved.
//  This file is distributed under the BSD 2-Clause License. See LICENSE.TXT
//  for details.
#include "permutationflowshop.h"
#include <cstdlib>
#include <climits>
#include <string>
#include <sstream>
#include <assert.h>
#include <algorithm>
#include <limits>
#include "sse_functions.h"

#ifdef NOC11
struct rzcomp{
    std::vector< float >* tas;
    bool operator()(int i1,int i2){
        if((*tas)[i1]==(*tas)[i2] && i2!=0)
            return i1>i2;
        else
           return (*tas)[i1] < (*tas)[i2];
}
}rzco;

struct stdstartComp
{
    std::vector< float >* stds;
    bool operator()(int i1,int i2){return (*stds)[i1] < (*stds)[i2];}
}stdc;

struct igioComp
{
    std::vector< int >* stds;
    bool operator()(int i1,int i2){return (*stds)[i1] > (*stds)[i2];}
}igioc;
#endif

std::vector< float > inline bs_start_sequence(emili::pfsp::PermutationFlowShop& prob)
{
    std::vector< float > xi;
    xi.push_back(0);
    // pmat[ machine ][ job ]
    const std::vector< std::vector < long int > > & pmat = prob.getProcessingTimesMatrix();
    int njob = prob.getNjobs();
    int nmac = prob.getNmachines();
    float k = (njob-2)/(float)4;
    for(int j = 1 ; j <= njob; j++)
    {
        float pi = pmat[1][j];
        float wj = 0;

        for( int i = 2; i <= nmac ; i++)
        {
            wj = (nmac*pi)/(i-1);
            pi += pmat[i][j];

        }
        xi[j] = pi + k * wj;
    }
    return xi;
}



std::vector< int > inline std_start_sequence(emili::pfsp::PermutationFlowShop& prob)
{
    int njobs = prob.getNjobs();
    int nmac = prob.getNmachines();
    const std::vector < std::vector < long int > >& priorities = prob.getProcessingTimesMatrix();
    std::vector< float > stds(njobs+1,0);
    std::vector< int > seq;
    seq.push_back(0);
    int i=1;
    for(;i<=njobs;i++)
    {
        seq.push_back(i);
        float avg = 0.0;
           for(int m=1;m<=nmac;m++)
           {
               avg += priorities[i][m];
           }
           avg = avg/nmac;
           stds[i] = avg;
    }

    float mm1 = 1/(nmac-1);
    for(i=1;i<=njobs;i++)
    {
        float sTd=0.0;
        for(int m=1;m <= nmac;m++)
        {
            sTd += powf((priorities[i][m]-stds[i]),2.0f);
        }
        sTd = sTd * mm1;
        stds[i] = sqrtf(sTd);
    }
#ifndef NOC11
    std::sort(seq.begin()+1,seq.end(),[stds](int i1,int i2){return stds[i1] < stds[i2];});
#else
    stdc.stds = &stds;
    std::sort(seq.begin()+1,seq.end(),stdc);
#endif

    return seq;

}

std::vector< int > inline rz_seed_sequence_simplified(emili::pfsp::PermutationFlowShop& prob)
{
    int machines = prob.getNmachines();
    int jobs = prob.getNjobs();
    const std::vector < std::vector < long int > >& priorities = prob.getProcessingTimesMatrix();
        std::vector< int > temp;
        std::vector< float > tas(jobs+1,0);
        temp.push_back(0);
        for(int i=1;i<=jobs;i++ )
        {
            float tai= 0;

            for (int j=1;j<=machines;j++)
            {
                tai += ( machines - j + 1 ) * priorities[i][j];
            }

            int pri = prob.getPriority(i);
            if(pri>0)
                tai = tai/pri;

            tas[i] = tai;
            temp.push_back(i);
        }

#ifndef NOC11
        std::sort(temp.begin(),temp.end(),[tas](int i1,int i2){
                                                                if(tas[i1]==tas[i2] && i2!=0)
                                                                    return i1>i2;
                                                                else
                                                                   return tas[i1] < tas[i2];
        });
#else
    rzco.tas = &tas;
    std::sort(temp.begin(),temp.end(),rzco);
#endif


        prob.computeObjectiveFunction(temp);





    return temp;
}

std::vector< int > inline rz_seed_sequence(emili::pfsp::PermutationFlowShop& prob)
{
    std::vector< int > best;
    int wbest = 0;
    int machines = prob.getNmachines();
    int jobs = prob.getNjobs();
    const std::vector < std::vector < long int > >& priorities = prob.getProcessingTimesMatrix();

    for(int k=1; k<=machines; k++)
    {
        std::vector< int > temp;
        std::vector< float > tas(jobs+1,0);
        temp.push_back(0);
        for(int i=1;i<=jobs;i++ )
        {
            float tai= 0;

            for (int j=k;j<=machines;j++)
            {
                tai += ( machines - j + 1 ) * priorities[i][j];
            }

            int pri = prob.getPriority(i);
            if(pri>0)
                tai = tai/pri;

            tas[i] = tai;
            temp.push_back(i);
        }

#ifndef NOC11
        std::sort(temp.begin(),temp.end(),[tas](int i1,int i2){
                                                                if(tas[i1]==tas[i2] && i2!=0)
                                                                    return i1>i2;
                                                                else
                                                                   return tas[i1] < tas[i2];
        });
#else
    rzco.tas = &tas;
    std::sort(temp.begin(),temp.end(),rzco);
#endif


        int w = prob.computeObjectiveFunction(temp);

        if(wbest == 0 || wbest > w)
        {

            wbest = w;
            best = temp;
        }
    }


    return best;
}

std::vector< int > inline rz_seed_sequence(std::vector< int > partial, std::vector< int > removed,emili::pfsp::PermutationFlowShop& prob)
{
    std::vector< int > best(removed);
    std::vector< int > best_sol(partial);
    best_sol.insert(best_sol.end(),removed.begin(),removed.end());
    int wbest = prob.computeObjectiveFunction(best_sol);
    int machines = prob.getNmachines();
    int jobs = removed.size();
    const std::vector < std::vector < long int > >& priorities = prob.getProcessingTimesMatrix();

    for(int k=1; k<=machines; k++)
    {
        std::vector< int > temp;
        std::vector< float > tas(prob.getNjobs()+1,0);
        //temp.push_back(0);
        for(int i=0;i<jobs;i++ )
        {
            float tai= 0;

            for (int j=k;j<=machines;j++)
            {
                tai += ( machines - j + 1 ) * priorities[removed[i]][j];
            }
            int pri = prob.getPriority(i);
            if(pri>0)
                tai = tai/pri;
            tas[removed[i]] = tai;
            temp.push_back(removed[i]);
        }
#ifndef NOC11
    std::sort(temp.begin(),temp.end(),[tas](int i1,int i2){
                                                            if(tas[i1]==tas[i2])
                                                                return i1>i2;
                                                            else
                                                               return tas[i1] < tas[i2];
    });
#else
        rzco.tas = &tas;
        std::sort(temp.begin(),temp.end(),rzco);
#endif

        std::vector< int > test_sol(partial);
        test_sol.insert(test_sol.end(),temp.begin(),temp.end());
        int w = prob.computeObjectiveFunction(test_sol);

        if( wbest > w)
        {

            wbest = w;
            best = temp;
            best_sol = test_sol;
        }
    }


    return best_sol;
}


std::vector< int > inline rz_improvement_phase(std::vector<int>& start_seq, emili::pfsp::PermutationFlowShop& prob)
{
    int jobs = prob.getNjobs();
    std::vector < int > res(start_seq);
    int wres = prob.computeObjectiveFunction(start_seq);
    for(int i = 1; i <= jobs; i++ )
    {
        int spos = start_seq[i];
        for(int j = 1; j<= jobs; j++ )
        {
            if( j != i )
            {
                std::vector< int > temp(start_seq);
                temp.erase(temp.begin()+i);
                temp.insert(temp.begin()+j,spos);
                int wtemp = prob.computeObjectiveFunction(temp);
                if( wres > wtemp)
                {
                    wres = wtemp;
                    res = temp;
                }
            }
        }
    }

    return res;
}


std::vector< int > inline neh(std::vector< int >& partial,int nbJobs,emili::pfsp::PermutationFlowShop& pis)
{
    std::vector< int >  sol(nbJobs+1,0);
    sol[1] = partial[1];
    sol[2] = partial[2];
    int wt = pis.computeObjectiveFunction(sol);
    sol[1] = partial[2];
    sol[2] = partial[1];
    int wt2 = pis.computeObjectiveFunction(sol);
    if(wt2>wt){
        sol[1] = partial[1];
        sol[2] = partial[2];
    }
    for(int i=3; i <= nbJobs; i++)
    {
        int candidate = partial[i];
        int bestpos = 1;
        sol.insert(sol.begin()+1,candidate);
        int wt_min = pis.computeObjectiveFunction(sol);
        for(int j = 2;j<i;j++)
        {
            sol.erase(sol.begin()+(j-1));
            sol.insert(sol.begin()+j,candidate);
            int wt = pis.computeObjectiveFunction(sol);
            if(wt<wt_min)
            {
                wt_min = wt;
                bestpos = j;
            }
        }
    }
    sol.erase(sol.begin()+nbJobs+1,sol.end());
    return sol;
}

std::vector< int > inline neh2(std::vector<int >& _fsp, int N, emili::pfsp::PermutationFlowShop& pis)
{
            int min;
            int tmp,ind;
    std::vector< int >  solTMP(N+1,0);
            solTMP[1]=_fsp[2];
            solTMP[2]=_fsp[1];

            int mS=pis.computeObjectiveFunction(_fsp,2);//compute_total_wt(_fsp,2);
            if(pis.computeObjectiveFunction(solTMP,2)<mS){//compute_total_wt(solTMP,2)<mS){
                _fsp[1]=solTMP[1];
                _fsp[2]=solTMP[2];
            }

            for(int k=3;k<=N;k++){
                    min = std::numeric_limits<int>::max();//min=10000000;
                for(int r=1; r<=k; r++){

                    for(int h=1; h<r; h++)
                        solTMP[h]=_fsp[h];
                    solTMP[r]=_fsp[k];
                    for(int h=r+1; h<=k; h++)
                        solTMP[h]=_fsp[h-1];

                    tmp=pis.computeObjectiveFunction(solTMP,k);//compute_total_wt(solTMP,k+1);
                    if(tmp<min){
                        min=tmp;
                        ind=r;
                    }

                }

                for(int h=0; h<ind; h++)
                    solTMP[h]=_fsp[h];
                solTMP[ind]=_fsp[k];
                for(int h=ind+1; h<=k; h++)
                    solTMP[h]=_fsp[h-1];

                for(int h=0; h<=k; ++h)
                    _fsp[h]=solTMP[h];
            }

         return _fsp;
}

std::vector< int > inline nehrs(std::vector<int >& _fsp, int N, emili::pfsp::PermutationFlowShop& pis)
{
            int min;
            int tmp,ind;
    std::vector< int >  solTMP(N+1,0);
 /*          solTMP[1]=_fsp[2];
               solTMP[2]=_fsp[1];
            int mS=pis.computeObjectiveFunction(_fsp,2);//compute_total_wt(_fsp,2);
            if(pis.computeObjectiveFunction(solTMP,2)<mS){//compute_total_wt(solTMP,2)<mS){
                _fsp[1]=solTMP[1];
                _fsp[2]=solTMP[2];
            }
*/
            int ind1 = emili::generateRandomNumber()%N +1;
            int ind2 = emili::generateRandomNumber()%N +1;

            int swp = _fsp[1];
            _fsp[1] = _fsp[ind1];
            _fsp[ind1] = swp;

            swp = _fsp[2];
            _fsp[2] = _fsp[ind2];
            _fsp[ind2] = swp;

            /*solTMP[1]=_fsp[2];
            solTMP[2]=_fsp[1];
            int mS=pis.computeObjectiveFunction(_fsp,2);//compute_total_wt(_fsp,2);
            if(pis.computeObjectiveFunction(solTMP,2)<mS){//compute_total_wt(solTMP,2)<mS){
                _fsp[1]=solTMP[1];
                _fsp[2]=solTMP[2];
            }*/
            for(int k=3;k<=N;k++){
                    min = std::numeric_limits<int>::max();//min=10000000;
                for(int r=1; r<=k; r++){

                    for(int h=1; h<r; h++)
                        solTMP[h]=_fsp[h];
                    solTMP[r]=_fsp[k];
                    for(int h=r+1; h<=k; h++)
                        solTMP[h]=_fsp[h-1];

                    tmp=pis.computeObjectiveFunction(solTMP,k);//compute_total_wt(solTMP,k+1);
                    if(tmp<min){
                        min=tmp;
                        ind=r;
                    }

                }

                for(int h=0; h<ind; h++)
                    solTMP[h]=_fsp[h];
                solTMP[ind]=_fsp[k];
                for(int h=ind+1; h<=k; h++)
                    solTMP[h]=_fsp[h-1];

                for(int h=0; h<=k; ++h)
                    _fsp[h]=solTMP[h];
            }

         return _fsp;
}

int inline nwnehrs(std::vector<int >& _fsp, int N, emili::pfsp::PermutationFlowShop& pis, const std::vector < std::vector < int > >& distances)
{
            int min;
            int tmp,ind;
            std::vector< int >  solTMP(N+1,0);
 /*          solTMP[1]=_fsp[2];
               solTMP[2]=_fsp[1];
            int mS=pis.computeObjectiveFunction(_fsp,2);//compute_total_wt(_fsp,2);
            if(pis.computeObjectiveFunction(solTMP,2)<mS){//compute_total_wt(solTMP,2)<mS){
                _fsp[1]=solTMP[1];
                _fsp[2]=solTMP[2];
            }
*/
            int ind1 = emili::generateRandomNumber()%N +1;
            int ind2 = emili::generateRandomNumber()%N +1;

            int swp = _fsp[1];
            _fsp[1] = _fsp[ind1];
            _fsp[ind1] = swp;

            swp = _fsp[2];
            _fsp[2] = _fsp[ind2];
            _fsp[ind2] = swp;
            int makespan = distances[0][_fsp[1]]+distances[_fsp[1]][_fsp[2]];
            for(int k=3;k<=N;k++){
                    min = std::numeric_limits<int>::max();//min=10000000;
                for(int r=1; r<=k; r++){

                    for(int h=1; h<r; h++)
                        solTMP[h]=_fsp[h];
                    solTMP[r]=_fsp[k];
                    for(int h=r+1; h<=k; h++)
                        solTMP[h]=_fsp[h-1];

                 //   int tmp2=pis.computeObjectiveFunction(solTMP,k);//compute_total_wt(solTMP,k+1);
                    int kk = solTMP[r];
                    int kmo = solTMP[r-1];
                    int kpo = solTMP[(r+1)%(N+1)];
                    tmp = makespan + distances[kmo][kk] + distances[kk][kpo] - distances[kmo][kpo];
                   // assert(tmp2 == tmp);
                    if(tmp<min){
                        min=tmp;
                        ind=r;
                    }

                }

                for(int h=0; h<ind; h++)
                    solTMP[h]=_fsp[h];
                solTMP[ind]=_fsp[k];
                for(int h=ind+1; h<=k; h++)
                    solTMP[h]=_fsp[h-1];

                for(int h=0; h<=k; ++h)
                    _fsp[h]=solTMP[h];

                 makespan = min;
            }

         return makespan;
}

std::vector< int > inline nehls(std::vector<int >& _fsp, int N, emili::pfsp::PermutationFlowShop& pis,emili::LocalSearch* ls)
{
            emili::pfsp::PermutationFlowShopSolution sol(std::numeric_limits<int>::max());
            PfspInstance& pfinstance = ((emili::pfsp::PermutationFlowShop&)ls->getInitialSolution().getProblem()).getInstance();
            emili::pfsp::PfspNeighborhood& pneigh =((emili::pfsp::PfspNeighborhood&)ls->getNeighborhood());
            int min;
            int tmp,ind;
    std::vector< int >  solTMP(N+1,0);
            solTMP[1]=_fsp[2];
            solTMP[2]=_fsp[1];

            int mS=pis.computeObjectiveFunction(_fsp,2);//compute_total_wt(_fsp,2);
            if(pis.computeObjectiveFunction(solTMP,2)<mS){//compute_total_wt(solTMP,2)<mS){
                _fsp[1]=solTMP[1];
                _fsp[2]=solTMP[2];
            }

            for(int k=3;k<=N;k++){
                    min = std::numeric_limits<int>::max();//min=10000000;
                for(int r=1; r<=k; r++){

                    for(int h=1; h<r; h++)
                        solTMP[h]=_fsp[h];
                    solTMP[r]=_fsp[k];
                    for(int h=r+1; h<=k; h++)
                        solTMP[h]=_fsp[h-1];

                    tmp=pis.computeObjectiveFunction(solTMP,k);//compute_total_wt(solTMP,k+1);
                    if(tmp<min){
                        min=tmp;
                        ind=r;
                    }

                }

                for(int h=0; h<ind; h++)
                    solTMP[h]=_fsp[h];
                solTMP[ind]=_fsp[k];
                for(int h=ind+1; h<=k; h++)
                    solTMP[h]=_fsp[h-1];

                for(int h=0; h<=k; ++h)
                    _fsp[h]=solTMP[h];
                //sol.setSolutionValue(min);
                //will it work??
        /**         for(int i =1; i<=k;i++)
                {
                    std::cout << " " << _fsp[i];
                }
                std::cout << std::endl;
                std::cout << min << std::endl;*/
                sol.setJobSchedule(_fsp);
                pfinstance.setNbJob(k);
                pneigh.setNjobs(k);
                emili::pfsp::PermutationFlowShopSolution* s2 = (emili::pfsp::PermutationFlowShopSolution*)ls->search(&sol);
                _fsp = s2->getJobSchedule();
                delete s2;
            }

         return _fsp;
}

std::vector< int > inline nehffls(std::vector<int >& _fsp,
                                int N,
                                emili::pfsp::PermutationFlowShop& pis,emili::LocalSearch* ls
                                )
{
    emili::pfsp::PermutationFlowShopSolution sol(std::numeric_limits<int>::max());
    PfspInstance& pfinstance = ((emili::pfsp::PermutationFlowShop&)ls->getInitialSolution().getProblem()).getInstance();
    emili::pfsp::PfspNeighborhood& pneigh =((emili::pfsp::PfspNeighborhood&)ls->getNeighborhood());
    //tmp=compute_total_wt(solTMP,sizePartial+1);
    //                  std::cout << "start perturb" << std::endl;
    //check why plus 1
    const std::vector < std::vector < long int > >& pmatrix = pis.getProcessingTimesMatrix();

    int min;
    int tmp,ind;
    std::vector< int >  solTMP(N+1,0);
    int m = pis.getNmachines();
    std::vector < std::vector < int > > head(m+1,std::vector< int > (pis.getNjobs()+1,0));
    std::vector < std::vector < int > > tail(m+1,std::vector< int > (pis.getNjobs()+1,0));


    solTMP[1]=_fsp[2];
    solTMP[2]=_fsp[1];

    int mS=pis.computeMS(_fsp,2);//compute_total_wt(_fsp,2);
    if(pis.computeMS(solTMP,2)<mS){//compute_total_wt(solTMP,2)<mS){
        _fsp[1]=solTMP[1];
        _fsp[2]=solTMP[2];
    }

    for(int k=3;k<=N;k++){
            min = std::numeric_limits<int>::max();//min=10000000;
            //pis.computeTAmatrices(solTMP,head,tail,solTMP.size());
#ifdef ENABLE_SSE
        computeHEADandTAIL(solTMP,head,tail,pmatrix,solTMP.size()-1,m);
#else
        pis.computeTAmatrices(solTMP,head,tail,solTMP.size());
#endif
            std::vector< int >  ptb;
            int kk = _fsp[k];
        for(int r=1; r<=k; r++){

            for(int h=1; h<r; h++)
                solTMP[h]=_fsp[h];
            solTMP[r]=_fsp[k];
            for(int h=r+1; h<=k; h++)
                solTMP[h]=_fsp[h-1];

            //tmp=pis.computeObjectiveFunction(solTMP,k);//compute_total_wt(solTMP,k+1);

            long int c_cur = head[1][r-1]+pmatrix[kk][1];
            long int c_max = c_cur+tail[1][r];
            for (int i = 2; i <= m; ++i) {
                int c_pm = head[i][r-1];
                if(c_pm < c_cur)
                {
                    c_cur = c_cur + pmatrix[kk][i];
                }
                else
                {
                    c_cur = c_pm + pmatrix[kk][i];
                }
                long int c_can = (c_cur+tail[i][r]);
                c_max = c_max>c_can?c_max:c_can;
            }
            tmp = c_max;

            if(tmp<min){
                min=tmp;
                ind=r;
                ptb.clear();
                ptb.push_back(ind);
            }else if(tmp==min)
            {
                ptb.push_back(r);
            }

        }
    int tb = ptb.size();
        if(tb > 1 && k<N)// if there are ties to break...
        {
            //tie breaker!
            int bp = ptb[0]; // insert position for the minimum idle time
            int itbp = std::numeric_limits<int> ::max(); // minimum idle time estimation
            for(int l=0; l < tb; l++) // for each tie location
            {
                int ptbl = ptb[l];
                int itdd = 0; // idle time estimation
                int fil = 0; // make span of kk when inserted in ptbl
                if(ptbl == k) // if insert position is the last position in the partial solution
                {
                    fil =  head[1][ptbl-1]+pmatrix[kk][1];
                    for(int i=2 ; i <= m; i++)
                    {
                        fil = std::max(head[i][ptbl-1],fil) + pmatrix[kk][i];
                        itdd += fil-head[i][ptbl-1]-pmatrix[kk][i]; //

                    }
                }
                else
                {
                    int sptbl = solTMP[ptbl];
                    fil =  head[1][ptbl-1]+pmatrix[kk][1];
                    int filp =  fil+ pmatrix[sptbl][1];
                    for(int i=2 ; i <= m; i++)
                    {
                        fil = std::max(head[i][ptbl-1],fil) + pmatrix[kk][i];
                        itdd += fil-head[i][ptbl-1]-pmatrix[kk][i]+pmatrix[sptbl][i]+std::max(filp-fil,0);
                        //itdd += fil-head[i][ptbl-1]+pmatrix[sptbl][i]+std::max(filp-fil,0);
                        filp = std::max(fil,filp) + pmatrix[sptbl][i];
                    }
                }
                if(itbp > itdd)
                {
                    itbp = itdd;
                    bp = ptbl;
                }
            }
            ind = bp;
        }
        //end tie braking
        for(int h=0; h<ind; h++)
            solTMP[h]=_fsp[h];
        solTMP[ind]=_fsp[k];
        for(int h=ind+1; h<=k; h++)
            solTMP[h]=_fsp[h-1];

        for(int h=0; h<=k; ++h)
            _fsp[h]=solTMP[h];

        sol.setJobSchedule(_fsp);
        pfinstance.setNbJob(k);
        pneigh.setNjobs(k);
        emili::pfsp::PermutationFlowShopSolution* s2 = (emili::pfsp::PermutationFlowShopSolution*)ls->search(&sol);
        _fsp = s2->getJobSchedule();
        delete s2;
    }

 return _fsp;
}

std::vector< int > inline nehff(std::vector<int >& _fsp,
                                int N,
                                emili::pfsp::PermutationFlowShop& pis
                                )
{
    //tmp=compute_total_wt(solTMP,sizePartial+1);
    //                  std::cout << "start perturb" << std::endl;
    //check why plus 1
    const std::vector < std::vector < long int > >& pmatrix = pis.getProcessingTimesMatrix();

    int min;
    int tmp,ind;
    std::vector< int >  solTMP(N+1,0);
    int m = pis.getNmachines();
    std::vector < std::vector < int > > head(m+1,std::vector< int > (pis.getNjobs()+1,0));
    std::vector < std::vector < int > > tail(m+1,std::vector< int > (pis.getNjobs()+1,0));


    solTMP[1]=_fsp[2];
    solTMP[2]=_fsp[1];

    int mS=pis.computeMS(_fsp,2);//compute_total_wt(_fsp,2);
    if(pis.computeMS(solTMP,2)<mS){//compute_total_wt(solTMP,2)<mS){
        _fsp[1]=solTMP[1];
        _fsp[2]=solTMP[2];
    }

    for(int k=3;k<=N;k++){
            min = std::numeric_limits<int>::max();//min=10000000;
            //pis.computeTAmatrices(solTMP,head,tail,solTMP.size());
#ifdef ENABLE_SSE
        computeHEADandTAIL(solTMP,head,tail,pmatrix,solTMP.size()-1,m);
#else
        pis.computeTAmatrices(solTMP,head,tail,solTMP.size());
#endif
            std::vector< int >  ptb;
            int kk = _fsp[k];
        for(int r=1; r<=k; r++){

            for(int h=1; h<r; h++)
                solTMP[h]=_fsp[h];
            solTMP[r]=_fsp[k];
            for(int h=r+1; h<=k; h++)
                solTMP[h]=_fsp[h-1];

            //tmp=pis.computeObjectiveFunction(solTMP,k);//compute_total_wt(solTMP,k+1);

            long int c_cur = head[1][r-1]+pmatrix[kk][1];
            long int c_max = c_cur+tail[1][r];
            for (int i = 2; i <= m; ++i) {
                int c_pm = head[i][r-1];
                if(c_pm > c_cur)
                {
                    c_cur = c_pm;
                }                
                c_cur = c_cur + pmatrix[kk][i];
                long int c_can = (c_cur+tail[i][r]);
                c_max = c_max>c_can?c_max:c_can;
            }
            tmp = c_max;

            if(tmp<min){
                min=tmp;
                ind=r;
                ptb.clear();
                ptb.push_back(ind);
            }else if(tmp==min)
            {
                ptb.push_back(r);
            }

        }
    int tb = ptb.size();
        if(tb > 1 && k<N)// if there are ties to break...
        {
            //tie breaker!
            int bp = ptb[0]; // insert position for the minimum idle time
            int itbp = std::numeric_limits<int> ::max(); // minimum idle time estimation
            for(int l=0; l < tb; l++) // for each tie location
            {
                int ptbl = ptb[l];
                int itdd = 0; // idle time estimation
                int fil = 0; // make span of kk when inserted in ptbl
                if(ptbl == k) // if insert position is the last position in the partial solution
                {
                    fil =  head[1][ptbl-1]+pmatrix[kk][1];
                    for(int i=2 ; i <= m; i++)
                    {
                        fil = std::max(head[i][ptbl-1],fil) + pmatrix[kk][i];
                        itdd += fil-head[i][ptbl-1]-pmatrix[kk][i]; //

                    }
                }
                else
                {
                    int sptbl = solTMP[ptbl];
                    fil =  head[1][ptbl-1]+pmatrix[kk][1];
                    int filp =  fil+ pmatrix[sptbl][1];
                    for(int i=2 ; i <= m; i++)
                    {
                        fil = std::max(head[i][ptbl-1],fil) + pmatrix[kk][i];
                        itdd += fil-head[i][ptbl-1]-pmatrix[kk][i]+pmatrix[sptbl][i]+std::max(filp-fil,0);
                        //itdd += fil-head[i][ptbl-1]+pmatrix[sptbl][i]+std::max(filp-fil,0);
                        filp = std::max(fil,filp) + pmatrix[sptbl][i];
                    }
                }
                if(itbp > itdd)
                {
                    itbp = itdd;
                    bp = ptbl;
                }
            }
            ind = bp;
        }
        //end tie braking
        for(int h=0; h<ind; h++)
            solTMP[h]=_fsp[h];
        solTMP[ind]=_fsp[k];
        for(int h=ind+1; h<=k; h++)
            solTMP[h]=_fsp[h-1];

        for(int h=0; h<=k; ++h)
            _fsp[h]=solTMP[h];
    }

 return _fsp;
}

std::vector< int > inline slack_construct(std::vector< int >& partial, int nbJobs,emili::pfsp::PermutationFlowShop& pis)
{
    int Ci = pis.computeMS(partial);
    std::vector<bool> assigned(nbJobs+1, false);
    for(std::vector< int >::const_iterator iter = partial.begin();iter !=  partial.end();++iter)
    {
        assigned[*iter] = true;
    }
//    assigned[0] = false;
    for (int var = 0; var < nbJobs; ++var) {
        int kish = partial[var+1];
        if(kish==0)
        {          
            int minJ = 0 ;
            int minE = INT_MAX;
            for (int jb = 1; jb <= nbJobs; ++jb) {
                if(!assigned[jb]){
                    int dd = pis.getDueDate(jb);
                    int pr = pis.getPriority(jb);
                    partial[var+1] = jb;
                    Ci = pis.computeMS(partial,var+1);
                    int wej = pr * (dd-Ci);
                    partial[var+1] = 0;
                    if(minE > wej){
                        minJ = jb;
                        minE = wej;
                    }
                }
            }
            partial[var+1] = minJ;
            assigned[minJ] = true;
            //Ci = pis.computeMS(partial);
        }

        //cout << "partial makespan: " << Ci << std::endl;
    }

    return partial;
}

std::vector< int > inline lit_construct(std::vector<int>& partial,std::vector<int> toFill, emili::pfsp::PermutationFlowShop& problem)
{
    PfspInstance instance = problem.getInstance();
    int njobs = instance.getNbJob();
    int nmac = instance.getNbMac();
    std::vector<int> res(partial);
    for(int i = partial.size(); i<= njobs; i++)
    {
        int index = 0;
        int nextJob = toFill[0];
        std::vector< int > prevJobC(nmac+1,0);
        std::vector< int > pmacend(njobs+1,0);
        instance.computeWTs(res,prevJobC,i-1,pmacend);
        int nja = instance.computeIdleTimeCoeff(prevJobC,nextJob);
        for(int j = 1;j<toFill.size();j++)
        {
            int tja = instance.computeIdleTimeCoeff(prevJobC,toFill[j]);
            if(tja > nja)
            {
                nja = tja;
                nextJob = toFill[j];
                index = j;
            }
        }
        res.push_back(nextJob);
        toFill.erase(toFill.begin()+index);
    }

    return res;
}

std::vector< float > inline lr_index(std::vector< int >& s, std::vector<int>& u, emili::pfsp::PermutationFlowShop& prob)
{
    PfspInstance instance = prob.getInstance();
    int njobs = instance.getNbJob();
    int nmac = instance.getNbMac();
    std::vector< int > lastJobCompletionTimes(nmac+1,0);
    std::vector< int > pmacend(njobs+1,0);
    if(s.size() > 1)
    {
        instance.computeWTs(s,lastJobCompletionTimes,s.size()-1,pmacend);
    }

    if(u.size() == 1)
    {
        return std::vector< float > (1,0.0f);
    }

    const std::vector< std::vector< long > >& ctimesMatrix = prob.getProcessingTimesMatrix();

    std::vector< float > findex(njobs+1,std::numeric_limits<float>::max());
    std::vector< float > tpj = std::vector< float> (nmac+1,0);


        for(int j=1;j<=nmac;j++)
        {
            for(int i=0;i<u.size();i++)
            {
                tpj[j] += ctimesMatrix[u[i]][j];
            }
            //tpj[j] = tpj[j];
        }

    int k = s.size();
    for(int i=0;i < u.size();i++)
    {
        int i_job = u[i];
        float ITik = 0;
        float i_comp = lastJobCompletionTimes[1]+ctimesMatrix[i_job][1];
        float p_comp = i_comp + (tpj[1]-ctimesMatrix[i_job][1])/(float)(u.size()-1);
        for(int j=2;j<=nmac;j++)
        {
            /** ITK calculation*/
            float wjk = nmac/(j+k*((nmac-j)/(float)(njobs-2)));

            float tcomp = std::max(i_comp-lastJobCompletionTimes[j],0.0f);


            ITik += wjk*tcomp;

            i_comp = tcomp>0?i_comp:lastJobCompletionTimes[j];
            i_comp += ctimesMatrix[i_job][j];

            /**  P job completion time calculation*/
            float tpcj = (tpj[j]-ctimesMatrix[i_job][j])/(float)(u.size()-1);

            if(p_comp + tpcj > i_comp)
            {
                p_comp = p_comp + tpcj;
            }
            else
            {
                p_comp = i_comp + tpcj;
            }
        }
        float ATik = i_comp+p_comp;

        findex[i_job] = ((u.size()-1)*ITik)+ATik;
    }

    return findex;
}


int generateRndPos(int min, int max)
{
  return (  emili::generateRandomNumber()%max + min );
}

double emili::pfsp::PermutationFlowShop::calcObjectiveFunctionValue(Solution &solution)
{
    emili::pfsp::PermutationFlowShopSolution& s = dynamic_cast<emili::pfsp::PermutationFlowShopSolution&> (solution);
    return computeObjectiveFunction(s.getJobSchedule());
}

double emili::pfsp::PermutationFlowShop::evaluateSolution(emili::Solution& solution)
{    
    emili::pfsp::PermutationFlowShopSolution& s = dynamic_cast<emili::pfsp::PermutationFlowShopSolution&> (solution);
    double p = computeObjectiveFunction(s.getJobSchedule());
    solution.setSolutionValue(p);
    return p;
}

int emili::pfsp::PermutationFlowShop::getNjobs(){
    return instance.getNbJob();
}

int emili::pfsp::PermutationFlowShop::getNmachines()
{
    return instance.getNbMac();
}

int emili::pfsp::PermutationFlowShop::getDueDate(int job)
{
    return instance.getDueDate(job);
}

PfspInstance& emili::pfsp::PermutationFlowShop::getInstance()
{
    return instance;
}

int emili::pfsp::PermutationFlowShop::getPriority(int job)
{
    return instance.getPriority(job);
}

int emili::pfsp::PermutationFlowShop::computeMS(std::vector< int > & partial_solution)
{
    return instance.computeMS(partial_solution);
}

int emili::pfsp::PermutationFlowShop::computeMS(std::vector< int > & partial_solution,int size)
{
    return instance.computeMS(partial_solution,size);
}



int emili::pfsp::PermutationFlowShop::computeObjectiveFunction(std::vector<int> &sol,std::vector<int>& prevJob,int job,std::vector<int>& previousMachineEndTime)
{
    return instance.computeWT(sol,prevJob,job,previousMachineEndTime);
}

int emili::pfsp::PermutationFlowShop::computeObjectiveFunction(std::vector<int> &sol, std::vector<std::vector<int> > &previousMachineEndTimeMatrix, int start_i, int end_i)
{
    return instance.computeWT(sol,previousMachineEndTimeMatrix,start_i,end_i);
}

int emili::pfsp::PermutationFlowShop::computeObjectiveFunction(std::vector<int> &solution, std::vector<int>& makespans, int size)
{
    return this->computeObjectiveFunction(solution);
}

long int emili::pfsp::PermutationFlowShop::computeObjectiveFunctionFromHead(std::vector<int> &solution, int starting_point, std::vector < std::vector < int > >& head, int njobs)
{
    int nmac = this->getNmachines();
    //int njobs = this->getNjobs();
    std::vector< int > makespans = head[nmac];
//    instance.updateHead(solution,starting_point,head,makespans);
//    return computeObjectiveFunction(solution,makespans,njobs);
    const std::vector < std::vector < long int > >& pmatrix = this->getProcessingTimesMatrix();
    int sol_i = solution[starting_point];
    int ins_pos[nmac+1];
    long int c_cur = head[1][starting_point-1]+pmatrix[sol_i][1];
    ins_pos[1] = c_cur;
   // std::vector< int > makespans = head[nmac];
    //(njobs+1,0);
    //for( int i=0; i < starting_point; i++)
     //   makespans[i] = head[nmac][i];

    for (int i = 2; i <= nmac; ++i) {
        int c_pm = head[i][starting_point-1];
        if(c_pm < c_cur)
        {
            c_cur = c_cur + pmatrix[sol_i][i];
        }
        else
        {
            c_cur = c_pm + pmatrix[sol_i][i];
        }
        ins_pos[i] =  c_cur;

    }

    //int wt = (std::max(c_cur - this.getDueDate(sol_i), 0L) * this.getPriority(sol_i));
    makespans[starting_point] = c_cur;
//    long int pre_c_cur = c_cur;

 //   for (int j = 1; j< starting_point; ++j )
  //  {
  //      wt += (std::max((long int)head[nmac][j] - pis.getDueDate(solution[j]), 0L) * pis.getPriority(solution[j]));
  //  }

    for(int k=starting_point+1; k<= njobs; k++)
    {
        int job = solution[k];
        c_cur = ins_pos[1] + pmatrix[job][1];
        ins_pos[1] = c_cur;
        for(int m=2; m <= nmac ; m++)
        {
            int c_pm = ins_pos[m];
            if(c_pm > c_cur)
            {
                c_cur = c_pm;
            }
            c_cur += pmatrix[job][m];
            ins_pos[m] = c_cur;
        }
        makespans[k] = c_cur;
     //   wt += (std::max(pre_c_cur - pis.getDueDate(solution[k]), 0L) * pis.getPriority(solution[k]));
    }
    //value->setSolutionValue(wt);

    return computeObjectiveFunction(solution,makespans,njobs);
}

long int emili::pfsp::PermutationFlowShop::computeObjectiveFunctionFromHead(std::vector<int> &solution, int starting_point, std::vector < std::vector < int > >& head)
{
    int njobs = this->getNjobs();
    return computeObjectiveFunctionFromHead(solution,starting_point,head,njobs);
}


void emili::pfsp::PermutationFlowShop::computeWTs(std::vector<int> &sol,std::vector<int>& prevJob,int job,std::vector<int>& previousMachineEndTime)
{
   return instance.computeWTs(sol,prevJob,job,previousMachineEndTime);
}

void emili::pfsp::PermutationFlowShop::computeTAmatrices(std::vector<int> &sol,std::vector< std::vector < int > >& head, std::vector< std::vector< int > >& tail)
{
    instance.computeTAmatrices(sol,head,tail);
}

void emili::pfsp::PermutationFlowShop::computeTAmatrices(std::vector<int> &sol,std::vector< std::vector < int > >& head, std::vector< std::vector< int > >& tail, int size)
{
    instance.computeTAmatrices(sol,head,tail,size);
}

void emili::pfsp::PermutationFlowShop::computeNoIdleTAmatrices(std::vector<int> &sol, std::vector<std::vector<int> > &head, std::vector<std::vector<int> > &tail)
{
    instance.computeNoIdleTAmatrices(sol,head,tail);
}

void emili::pfsp::PermutationFlowShop::computeTails(std::vector<int> &sol, std::vector< std::vector< std::vector< int > > > & tails)
{
    instance.computeTails(sol,tails);
}

std::vector< long int >& emili::pfsp::PermutationFlowShop::getDueDates()
{
    return instance.getDueDates();
}

std::vector< long int >& emili::pfsp::PermutationFlowShop::getPriorities()
{
    return instance.getPriorities();
}

const void* emili::pfsp::PermutationFlowShopSolution::getRawData()const
{

    return &solution;

}

void emili::pfsp::PermutationFlowShopSolution::setRawData(const void *data)
{
    //
    std::vector < int >* data_vector = (std::vector< int >*) data;

    this->solution = *data_vector;


}


emili::Solution* emili::pfsp::PermutationFlowShopSolution::clone()
{
    return new emili::pfsp::PermutationFlowShopSolution(this->solution_value,this->solution);
}

std::vector< int > & emili::pfsp::PermutationFlowShopSolution::getJobSchedule()
{
    return this->solution;
}

void emili::pfsp::PermutationFlowShopSolution::setJobSchedule(std::vector<int>& newSeq)
{
    this->solution = newSeq;
}


const std::vector< std::vector < long int > > & emili::pfsp::PermutationFlowShop::getProcessingTimesMatrix()
{
    return instance.getProcessingTimesMatrix();
}

emili::pfsp::PermutationFlowShopSolution::~PermutationFlowShopSolution()
{
 /** nothing to delete*/
}

int emili::pfsp::PFSP_WT::computeObjectiveFunction(std::vector< int > & partial_solution)
{
    return instance.computeWT(partial_solution);
}

int emili::pfsp::PFSP_WT::computeObjectiveFunction(std::vector< int > & partial_solution,int size)
{
    return instance.computeWT(partial_solution,size);
}

int emili::pfsp::PFSP_WCT::computeObjectiveFunction(std::vector< int > & partial_solution)
{
    return instance.computeWCT(partial_solution);
}

int emili::pfsp::PFSP_WCT::computeObjectiveFunction(std::vector< int > & partial_solution,int size)
{
    return instance.computeWCT(partial_solution,size);
}

int emili::pfsp::PFSP_TCT::computeObjectiveFunction(std::vector< int > & partial_solution)
{
    return instance.computeTCT(partial_solution);
}

int emili::pfsp::PFSP_TCT::computeObjectiveFunction(std::vector< int > & partial_solution,int size)
{
    return instance.computeTCT(partial_solution,size);
}

int emili::pfsp::PFSP_TCT::computeObjectiveFunction(std::vector<int> &solution, std::vector<int>& makespans, int size)
{
    return instance.computeT(solution,makespans,size);
}

int emili::pfsp::PFSP_WE::computeObjectiveFunction(std::vector< int > & partial_solution)
{
    return -instance.computeWE(partial_solution);
}

int emili::pfsp::PFSP_WE::computeObjectiveFunction(std::vector< int > & partial_solution,int size)
{
    return -instance.computeWE(partial_solution,size);
}

int emili::pfsp::PFSP_T::computeObjectiveFunction(std::vector<int> &partial_solution)
{
    return instance.computeT(partial_solution);
}

int emili::pfsp::PFSP_T::computeObjectiveFunction(std::vector<int> &partial_solution, int size)
{
    return instance.computeT(partial_solution,size);
}

int emili::pfsp::PFSP_T::computeObjectiveFunction(std::vector<int> &solution, std::vector<int>& makespans, int size)
{
    return instance.computeT(solution,makespans,size);
}

int emili::pfsp::PFSP_E::computeObjectiveFunction(std::vector<int> &partial_solution)
{
    return instance.computeE(partial_solution);
}

int emili::pfsp::PFSP_E::computeObjectiveFunction(std::vector<int> &partial_solution, int size)
{
    return instance.computeE(partial_solution,size);
}

int emili::pfsp::PFSP_MS::computeObjectiveFunction(std::vector<int> &partial_solution)
{
    return instance.computeMS(partial_solution);
}

int emili::pfsp::PFSP_MS::computeObjectiveFunction(std::vector<int> &partial_solution, int size)
{
    return instance.computeMS(partial_solution,size);
}

void emili::pfsp::NWPFSP_MS::computeNoWaitTimeDistances()
{
    const std::vector< std::vector < long int > >& pmatrix = this->instance.getProcessingTimesMatrix();
    int nmac = this->instance.getNbMac();
    int nJobs = pmatrix.size()-1;

    for(int i = 1 ; i <= nJobs ; i++)
    {
        int p = 0;
        for(int m=1; m<=nmac; m++)
        {
            p+=pmatrix[i][m];
        }

        distances[0][i] = p;

        for(int j=1; j <= nJobs; j++)
        {
            if(i!=j)
            {
                int max_ctd = 0;
                for (int k = 1; k<= nmac ; k++)
                {
                    int temp_ctd = pmatrix[i][k];
                    for(int h=k; h <= nmac ; h++ )
                    {
                        temp_ctd += pmatrix[j][h]-pmatrix[i][h];
                    }
                    if(temp_ctd > max_ctd){
                        max_ctd = temp_ctd;
                    }
                }
                distances[i][j] = max_ctd;
            }
        }
    }
}

const std::vector< std::vector < int > >& emili::pfsp::NWPFSP_MS::getDistances()
{
    return distances;
}

int emili::pfsp::NWPFSP_MS::computeObjectiveFunction(std::vector<int> &partial_solution)
{
    return instance.computeNWMS(partial_solution);
}

int emili::pfsp::NWPFSP_MS::computeObjectiveFunction(std::vector<int> &partial_solution, int size)
{
    return instance.computeNWMS(partial_solution,size);
}

int emili::pfsp::NWPFSP_WT::computeObjectiveFunction(std::vector<int> &partial_solution)
{
    return instance.computeNWWT(partial_solution);
}

int emili::pfsp::NWPFSP_WT::computeObjectiveFunction(std::vector<int> &partial_solution, int size)
{
    return instance.computeNWWT(partial_solution,size);
}

int emili::pfsp::NWPFSP_WE::computeObjectiveFunction(std::vector<int> &partial_solution)
{
    return instance.computeNWWE(partial_solution);
}

int emili::pfsp::NWPFSP_WE::computeObjectiveFunction(std::vector<int> &partial_solution, int size)
{
    return instance.computeNWWE(partial_solution,size);
}

int emili::pfsp::NWPFSP_E::computeObjectiveFunction(std::vector<int> &partial_solution)
{
    return instance.computeNWE(partial_solution);
}

int emili::pfsp::NWPFSP_E::computeObjectiveFunction(std::vector<int> &partial_solution, int size)
{
    return instance.computeNWE(partial_solution,size);
}

int emili::pfsp::NWPFSP_T::computeObjectiveFunction(std::vector<int> &partial_solution)
{
    return instance.computeNWT(partial_solution);
}

int emili::pfsp::NWPFSP_T::computeObjectiveFunction(std::vector<int> &partial_solution, int size)
{
    return instance.computeNWT(partial_solution,size);
}

int emili::pfsp::NWPFSP_TCT::computeObjectiveFunction(std::vector<int> &partial_solution)
{
    return instance.computeNWTCT(partial_solution);
}

int emili::pfsp::NWPFSP_TCT::computeObjectiveFunction(std::vector<int> &partial_solution, int size)
{
    return instance.computeNWTCT(partial_solution,size);
}

int emili::pfsp::NWPFSP_WCT::computeObjectiveFunction(std::vector<int> &partial_solution)
{
    return instance.computeNWWCT(partial_solution);
}

int emili::pfsp::NWPFSP_WCT::computeObjectiveFunction(std::vector<int> &partial_solution, int size)
{
    return instance.computeNWWCT(partial_solution,size);
}

int emili::pfsp::NIPFSP_MS::computeObjectiveFunction(std::vector<int> &partial_solution)
{
    return instance.computeNIMS(partial_solution);
}

int emili::pfsp::NIPFSP_MS::computeObjectiveFunction(std::vector<int> &partial_solution, int size)
{
    return instance.computeNIMS(partial_solution,size);
}

int emili::pfsp::NIPFSP_WT::computeObjectiveFunction(std::vector<int> &partial_solution)
{
    return instance.computeNIWT(partial_solution);
}

int emili::pfsp::NIPFSP_WT::computeObjectiveFunction(std::vector<int> &partial_solution, int size)
{
    return instance.computeNIWT(partial_solution,size);
}

int emili::pfsp::NIPFSP_WE::computeObjectiveFunction(std::vector<int> &partial_solution)
{
    return instance.computeNIWE(partial_solution);
}

int emili::pfsp::NIPFSP_WE::computeObjectiveFunction(std::vector<int> &partial_solution, int size)
{
    return instance.computeNIWE(partial_solution,size);
}

int emili::pfsp::NIPFSP_E::computeObjectiveFunction(std::vector<int> &partial_solution)
{
    return instance.computeNIE(partial_solution);
}

int emili::pfsp::NIPFSP_E::computeObjectiveFunction(std::vector<int> &partial_solution, int size)
{
    return instance.computeNIE(partial_solution,size);
}

int emili::pfsp::NIPFSP_T::computeObjectiveFunction(std::vector<int> &partial_solution)
{
    return instance.computeNIT(partial_solution);
}

int emili::pfsp::NIPFSP_T::computeObjectiveFunction(std::vector<int> &partial_solution, int size)
{
    return instance.computeNIT(partial_solution,size);
}

int emili::pfsp::NIPFSP_WCT::computeObjectiveFunction(std::vector<int> &partial_solution)
{
    return instance.computeNIWCT(partial_solution);
}

int emili::pfsp::NIPFSP_WCT::computeObjectiveFunction(std::vector<int> &partial_solution, int size)
{
    return instance.computeNIWCT(partial_solution,size);
}

int emili::pfsp::NIPFSP_TCT::computeObjectiveFunction(std::vector<int> &partial_solution)
{
    return instance.computeNITCT(partial_solution);
}

int emili::pfsp::NIPFSP_TCT::computeObjectiveFunction(std::vector<int> &partial_solution, int size)
{
    return instance.computeNITCT(partial_solution,size);
}

void emili::pfsp::NI_A_PFSP_MS::calc_nims_base()
{
    const std::vector < std::vector < long > > & ptm = instance.getProcessingTimesMatrix();
    for (int i = 1; i <= instance.getNbJob(); ++i) {
        nims_base += ptm[i][1];
    }
}

int emili::pfsp::NI_A_PFSP_MS::computeObjectiveFunction(std::vector<int> &partial_solution)
{
    return instance.computeNIMS(partial_solution,nims_base);
}

int emili::pfsp::NI_A_PFSP_MS::computeObjectiveFunction(std::vector<int> &partial_solution, int size)
{
    return instance.computeNIMS(partial_solution,size);
}

int emili::pfsp::SDSTFSP_MS::computeObjectiveFunction(std::vector<int> &partial_solution)
{
    return instance.computeSDSTMS(partial_solution);
}

int emili::pfsp::SDSTFSP_MS::computeObjectiveFunction(std::vector<int> &partial_solution, int size)
{
    return instance.computeSDSTMS(partial_solution,size);
}

int emili::pfsp::SDSTFSP_WT::computeObjectiveFunction(std::vector<int> &partial_solution)
{
    return instance.computeSDSTWT(partial_solution);
}

int emili::pfsp::SDSTFSP_WT::computeObjectiveFunction(std::vector<int> &partial_solution, int size)
{
    return instance.computeSDSTWT(partial_solution,size);
}

int emili::pfsp::SDSTFSP_WE::computeObjectiveFunction(std::vector<int> &partial_solution)
{
    return instance.computeSDSTWE(partial_solution);
}

int emili::pfsp::SDSTFSP_WE::computeObjectiveFunction(std::vector<int> &partial_solution, int size)
{
    return instance.computeSDSTWE(partial_solution,size);
}

int emili::pfsp::SDSTFSP_T::computeObjectiveFunction(std::vector<int> &partial_solution)
{
    return instance.computeSDSTT(partial_solution);
}

int emili::pfsp::SDSTFSP_T::computeObjectiveFunction(std::vector<int> &partial_solution, int size)
{
    return instance.computeSDSTT(partial_solution,size);
}

int emili::pfsp::SDSTFSP_E::computeObjectiveFunction(std::vector<int> &partial_solution)
{
    return instance.computeSDSTE(partial_solution);
}

int emili::pfsp::SDSTFSP_E::computeObjectiveFunction(std::vector<int> &partial_solution, int size)
{
    return instance.computeSDSTE(partial_solution,size);
}

int emili::pfsp::SDSTFSP_TCT::computeObjectiveFunction(std::vector<int> &partial_solution)
{
    return instance.computeSDSTTCT(partial_solution);
}

int emili::pfsp::SDSTFSP_TCT::computeObjectiveFunction(std::vector<int> &partial_solution, int size)
{
    return instance.computeSDSTTCT(partial_solution,size);
}

int emili::pfsp::SDSTFSP_WCT::computeObjectiveFunction(std::vector<int> &partial_solution)
{
    return instance.computeSDSTWCT(partial_solution);
}

int emili::pfsp::SDSTFSP_WCT::computeObjectiveFunction(std::vector<int> &partial_solution, int size)
{
    return instance.computeSDSTWCT(partial_solution,size);
}

emili::Solution* emili::pfsp::PfspInitialSolution::generateEmptySolution()
{
    std::vector< int > empty(pis.getNjobs()+1);
    emili::pfsp::PermutationFlowShopSolution*  em =  new emili::pfsp::PermutationFlowShopSolution(empty);
    em->setSolutionValue(std::numeric_limits<double>::max());
    return em;
}

emili::Solution* emili::pfsp::PfspInitialSolution::generateSolution()
{
    return generate();
}

emili::Solution* emili::pfsp::PfspRandomInitialSolution::generate()
{
    int nbJobs = pis.getNjobs();
    std::vector< int >   sol(nbJobs+1, 0);
  std::vector<bool> alreadyTaken(nbJobs+1, false); // nbJobs elements with value false
  std::vector<int > choosenNumber(nbJobs+1, 0);

  int nbj;
  int rnd, i, j, nbFalse;

  nbj = 0;
  for (i = nbJobs; i >= 1; --i)
  {
    rnd = generateRndPos(1, i);
    nbFalse = 0;

    /**  find the rndth cell with value = false : */
    for (j = 1; nbFalse < rnd; ++j)
      if ( ! alreadyTaken[j] )
        ++nbFalse;
    --j;

    sol[j] = i;

    ++nbj;
    choosenNumber[nbj] = j;

    alreadyTaken[j] = true;
  }
  PermutationFlowShopSolution* s = new PermutationFlowShopSolution(sol);
  pis.evaluateSolution(*s);
  return s;
}

emili::Solution* emili::pfsp::RandomInitialSolution::generate()
{
    emili::Solution* base = emili::pfsp::RandomInitialSolution::generate();
    for (int i=1;i<numOfSols;i++)
    {
        emili::Solution* news = emili::pfsp::RandomInitialSolution::generate();
        if(base->getSolutionValue() > news->getSolutionValue())
        {
            delete base;
            base = news;
        }
        else
        {
            delete news;
        }
    }
    return base;
}

emili::Solution* emili::pfsp::PfspSlackInitialSolution::generate()
{
    int nbJobs = pis.getNjobs();
    std::vector< int >  sol(nbJobs+1, 0);
    int Ci  = 0;
    std::vector<bool> assigned(nbJobs+1, false);
    std::vector<int> partial(nbJobs+1,0);
    for (int var = 0; var < nbJobs; ++var) {
        int minJ = 0 ;
        int minE = INT_MAX;
        for (int jb = 1; jb <= nbJobs; ++jb) {
            if(!assigned[jb]){
                int dd = pis.getDueDate(jb);
                int pr = pis.getPriority(jb);
                int wej = pr * (dd-Ci);
                if(minE > wej){
                    minJ = jb;
                    minE = wej;
                }
            }
        }
        partial[var+1] = minJ;
        assigned[minJ] = true;
        Ci = pis.computeMS(partial);
        //cout << "partial makespan: " << Ci << std::endl;
    }
    sol = partial;
    PermutationFlowShopSolution* s = new PermutationFlowShopSolution(sol);
    pis.evaluateSolution(*s);
    return s;
}

emili::Solution* emili::pfsp::PfspNEHwslackInitialSolution::generate()
{

    //SLACK: At time t, the job with the minimum value of dj − Cj (s) is selected.
    int nbJobs = pis.getNjobs(); // number of jobs in the problem
    std::vector< int >  sol(nbJobs+1, 0); //vector that contains the solution
    int Ci  = 0;
    std::vector<bool> assigned(nbJobs+1, false); // std::vector that indicates if a job it's already assigned to the partial solution
    std::vector<int> partial(nbJobs+1,0); // partial solution
    for (int var = 0; var < nbJobs; ++var) {
        int minJ = 0 ;
        int minE = INT_MAX;
        for (int jb = 1; jb <= nbJobs; ++jb) { //this loop finds the job that added to the partial solution achieves the minimum earliness
            if(!assigned[jb]){ // if jb it's not already assigned to the partial solution
                int dd = pis.getDueDate(jb); // get the due date of jb
                int pr = pis.getPriority(jb); // get the weight of jb
                partial[var+1] = jb; // insert jb at the end of the partial solution
                Ci = pis.computeMS(partial,var+1); // compute partial solution make span
                int wej = pr * (dd-Ci); // weighted earliness of jb when put in var+1
                partial[var+1] = 0; // deletes jb from the partial solution
                if(minE > wej){
                    minJ = jb;
                    minE = wej; //updates the current minimum earliness
                }
            }
        }
        partial[var+1] = minJ; // updates the partial solution inserting at the end the job with the minimum earliness
        assigned[minJ] = true; // updates the assigned vector
    }
    //int partial_w = pis.computeWT(partial);
    sol = neh2(partial,nbJobs,pis);
    PermutationFlowShopSolution* s = new PermutationFlowShopSolution(sol);
    pis.evaluateSolution(*s);
    //std::vector< std::vector< int > > etm = std::vector< std::vector< int > >(pis.getNmachines()+1,std::vector<int>(pis.getNjobs()+1,0));
    //clock_t start = clock();
    //int new_value =  pis.computeObjectiveFunction(sol,etm,1,pis.getNjobs()+1);

    //s->setSolutionValue(new_value);
    return s;
}

emili::Solution* emili::pfsp::NEH::generate()
{
    // NEH initial solution
    int njobs = pis.getNjobs();
    int nmac = pis.getNmachines();
    std::vector< int > tpt(njobs+1,0);
    std::vector< int > order;
    const std::vector< std::vector < long > >& ptm = pis.getProcessingTimesMatrix();
  //  order.push_back(0);
    for (int i = 1; i <= njobs; ++i) {
        int tpti = 0;
        for (int k = 1; k <= nmac; ++k) {
            tpti += ptm[i][k];
        }
        tpt[i] = tpti;
        order.push_back(i);
    }
   // std::sort(order.begin(),order.end(),[tpt](int i1,int i2){return tpt[i1] > tpt[i2];});
#ifndef NOC11
    std::sort(order.begin(),order.end(),[tpt](int i1,int i2){
        if(tpt[i1] == tpt[i2])
        {
            return i1>i2;
        }
        else
        {
            return tpt[i1] > tpt[i2];
        }
        });
#else
    igioc.stds=&tpt;
    std::sort(order.begin(),order.end(),igioc);
#endif
//    order.erase(order.begin()+njobs);
    order.insert(order.begin(),0);
    order = neh2(order,njobs,pis);
    PermutationFlowShopSolution* s = new PermutationFlowShopSolution(order);
    pis.evaluateSolution(*s);
    return s;
}

emili::Solution* emili::pfsp::NEHRS::generate()
{
    // NEH initial solution
    int njobs = pis.getNjobs();
    int nmac = pis.getNmachines();
    std::vector< int > tpt(njobs+1,0);
    std::vector< int > order;
    const std::vector< std::vector < long > >& ptm = pis.getProcessingTimesMatrix();
  //  order.push_back(0);
    for (int i = 1; i <= njobs; ++i) {
        int tpti = 0;
        for (int k = 1; k <= nmac; ++k) {
            tpti += ptm[i][k];
        }
        tpt[i] = tpti;
        order.push_back(i);
    }
   // std::sort(order.begin(),order.end(),[tpt](int i1,int i2){return tpt[i1] > tpt[i2];});
#ifndef NOC11
    std::sort(order.begin(),order.end(),[tpt](int i1,int i2){
        if(tpt[i1] == tpt[i2])
        {
            return i1>i2;
        }
        else
        {
            return tpt[i1] > tpt[i2];
        }
        });
#else
    igioc.stds=&tpt;
    std::sort(order.begin(),order.end(),igioc);
#endif
//    order.erase(order.begin()+njobs);
    order.insert(order.begin(),0);
    std::vector<int> original_order = order;
    std::vector<int> best_permutation;
    int value = std::numeric_limits<int>::max();
    int index = this->iterations;
    do
    {
        order = original_order;
        std::vector<int> ord = nehrs(order,njobs,pis);
        int v = pis.computeObjectiveFunction(order);
        if(v < value)
        {
           best_permutation = ord;
           value = v;
        }
        index--;
    }while(index > 0);
    PermutationFlowShopSolution* s = new PermutationFlowShopSolution(best_permutation);
    s->setSolutionValue(value);
    return s;
}

/**
 * Orders the jobs by due date to produce the seed sequence for neh
 */
emili::Solution* emili::pfsp::NEHedd::generate()
{
    // NEH initial solution
    int njobs = pis.getNjobs();
    std::vector< int > tpt(njobs+1,0);
    std::vector< int > order;
 //   order.push_back(0);
    for (int i = 1; i <= njobs; ++i) {
        tpt[i] = pis.getDueDate(i);
        order.push_back(i);
    }
#ifndef NOC11
    std::sort(order.begin(),order.end(),[tpt](int i1,int i2){return tpt[i1] < tpt[i2];});
#else
    igioc.stds=&tpt;
    std::sort(order.begin(),order.end(),igioc);
#endif
//    order.erase(order.begin()+njobs);
    order.insert(order.begin(),0);
    order = neh2(order,njobs,pis);
    PermutationFlowShopSolution* s = new PermutationFlowShopSolution(order);
    pis.evaluateSolution(*s);
    return s;
}

emili::Solution* emili::pfsp::NEHls::generate()
{
    // NEH initial solution
    int njobs = pis.getNjobs();
    int nmac = pis.getNmachines();
    std::vector< int > tpt(njobs+1,0);
    std::vector< int > order;
    const std::vector< std::vector < long > >& ptm = pis.getProcessingTimesMatrix();
//    order.push_back(0);
    for (int i = 1; i <= njobs; ++i) {
        int tpti = 0;
        for (int k = 1; k <= nmac; ++k) {
            tpti += ptm[i][k];
        }
        tpt[i] = tpti;
        order.push_back(i);
    }
#ifndef NOC11
    std::sort(order.begin(),order.end(),[tpt](int i1,int i2){
        if(tpt[i1] == tpt[i2])
        {
            return i1>i2;
        }
        else
        {
            return tpt[i1] > tpt[i2];
        }
        });
#else
    igioc.stds=&tpt;
    std::sort(order.begin(),order.end(),igioc);
#endif
//    order.erase(order.begin()+njobs);
    order.insert(order.begin(),0);
#ifdef WITH_STATS
    bool print_state = get_print();
    set_print(false);
#endif
    order = nehls(order,njobs,pis,_ls);
#ifdef WITH_STATS
    set_print(print_state);
#endif
    PermutationFlowShopSolution* s = new PermutationFlowShopSolution(order);
    pis.evaluateSolution(*s);
    return s;
}

emili::Solution* emili::pfsp::NEHeddLS::generate()
{
    // NEH initial solution
    int njobs = pis.getNjobs();
    std::vector< int > tpt(njobs+1,0);
    std::vector< int > order;
 //   order.push_back(0);
    for (int i = 1; i <= njobs; ++i) {
        tpt[i] = pis.getDueDate(i);
        order.push_back(i);
    }
#ifndef NOC11
    std::sort(order.begin(),order.end(),[tpt](int i1,int i2){return tpt[i1] < tpt[i2];});
#else
    igioc.stds=&tpt;
    std::sort(order.begin(),order.end(),igioc);
#endif
//    order.erase(order.begin()+njobs);
    order.insert(order.begin(),0);
#ifdef WITH_STATS
    bool print_state = get_print();
    set_print(false);
#endif
    order = nehls(order,njobs,pis,_ls);
#ifdef WITH_STATS
    set_print(print_state);
#endif
    PermutationFlowShopSolution* s = new PermutationFlowShopSolution(order);
    pis.evaluateSolution(*s);
    return s;
}

emili::Solution* emili::pfsp::NEHffls::generate()
{
    // NEH initial solution
    int njobs = pis.getNjobs();
    int nmac = pis.getNmachines();
    std::vector< int > tpt(njobs+1,0);
    std::vector< int > order;
    const std::vector< std::vector < long > >& ptm = pis.getProcessingTimesMatrix();
 //   order.push_back(0);
    for (int i = 1; i <= njobs; ++i) {
        int tpti = 0;
        for (int k = 1; k <= nmac; ++k) {
            tpti += ptm[i][k];
        }
        tpt[i] = tpti;
        order.push_back(i);
    }
#ifndef NOC11
    std::sort(order.begin(),order.end(),[tpt](int i1,int i2){return tpt[i1] > tpt[i2];});
#else
    igioc.stds=&tpt;
    std::sort(order.begin(),order.end(),igioc);
#endif
  //  order.erase(order.begin()+njobs);
    order.insert(order.begin(),0);
#ifdef WITH_STATS
    bool print_state = get_print();
    set_print(false);
#endif
    order = nehffls(order,njobs,pis,_ls);
#ifdef WITH_STATS
    set_print(print_state);
#endif
    PermutationFlowShopSolution* s = new PermutationFlowShopSolution(order);
    pis.evaluateSolution(*s);
    return s;
}

emili::Solution* emili::pfsp::NEHff::generate()
{
    // NEH initial solution
    int njobs = pis.getNjobs();
    int nmac = pis.getNmachines();
    std::vector< int > tpt(njobs+1,0);
    std::vector< int > order;
    const std::vector< std::vector < long > >& ptm = pis.getProcessingTimesMatrix();
//    order.push_back(0);
    for (int i = 1; i <= njobs; ++i) {
        int tpti = 0;
        for (int k = 1; k <= nmac; ++k) {
            tpti += ptm[i][k];
        }
        tpt[i] = tpti;
        order.push_back(i);
    }
#ifndef NOC11
    std::sort(order.begin(),order.end(),[tpt](int i1,int i2){return tpt[i1] > tpt[i2];});
#else
    igioc.stds=&tpt;
    std::sort(order.begin(),order.end(),igioc);
#endif
//    order.erase(order.begin()+njobs);
    order.insert(order.begin(),0);
    order = nehff(order,njobs,pis);
    PermutationFlowShopSolution* s = new PermutationFlowShopSolution(order);
    pis.evaluateSolution(*s);
    return s;
}

emili::Solution* emili::pfsp::LITSolution::generate()
{
    int nbJobs = pis.getNjobs();
    std::vector< int > res;
    res.push_back(0);
    res.push_back(1);
    std::vector< int > toFill;
    for(int i=2;i<=nbJobs;i++)
    {
        toFill.push_back(i);
    }
    std::vector< int > sol = lit_construct(res,toFill,pis);
    int best = pis.computeObjectiveFunction(sol);
    for(int i=2; i<=nbJobs; i++)
    {
        int j = res[1];
        res.erase(res.begin()+1);
        res.push_back(toFill[0]);
        toFill.erase(toFill.begin());
        toFill.push_back(j);
        std::vector<int> temp = lit_construct(res,toFill,pis);
        int bt = pis.computeObjectiveFunction(temp);
        if(bt<best)
        {
            sol = temp;
            best = bt;
        }
    }
   // sol = neh2(sol,nbJobs,pis);
   // best = pis.computeWT(sol);
    return new PermutationFlowShopSolution(best,sol);
}

emili::Solution* emili::pfsp::RZSolution::generate()
{
    std::vector< int > initial = rz_seed_sequence(pis);
    initial = rz_improvement_phase(initial,pis);
    //initial = neh2(initial,pis.getNjobs(),pis);
    PermutationFlowShopSolution* s = new PermutationFlowShopSolution(initial);
    pis.evaluateSolution(*s);
    return s;
}

emili::Solution* emili::pfsp::NeRZSolution::generate()
{
    std::vector< int > initial = rz_seed_sequence(pis);
    initial = rz_improvement_phase(initial,pis);
    initial = neh2(initial,pis.getNjobs(),pis);
    PermutationFlowShopSolution* s = new PermutationFlowShopSolution(initial);
    pis.evaluateSolution(*s);
    return s;
}

emili::Solution* emili::pfsp::NeRZ2Solution::generate()
{
    std::vector< int > initial = rz_seed_sequence(pis);
    //initial = rz_improvement_phase(initial,pis);
    initial = neh2(initial,pis.getNjobs(),pis);

    PermutationFlowShopSolution* s = new PermutationFlowShopSolution(initial);
    pis.evaluateSolution(*s);
    return s;
}

emili::Solution* emili::pfsp::NfRZ2Solution::generate()
{
    std::vector< int > initial = rz_seed_sequence(pis);
    //initial = rz_improvement_phase(initial,pis);
    initial = nehff(initial,pis.getNjobs(),pis);

    PermutationFlowShopSolution* s = new PermutationFlowShopSolution(initial);
    pis.evaluateSolution(*s);
    return s;
}

emili::Solution* emili::pfsp::SRZSolution::generate()
{
    std::vector< int > initial = rz_seed_sequence_simplified(pis);
    initial = neh2(initial,pis.getNjobs(),pis);
    PermutationFlowShopSolution* s = new PermutationFlowShopSolution(initial);
    pis.evaluateSolution(*s);
    return s;
}

emili::Solution* emili::pfsp::MNEH::generate()
{
    std::vector< int > initial = std_start_sequence(pis);
    //initial = rz_improvement_phase(initial,pis);
    initial = neh2(initial,pis.getNjobs(),pis);

    PermutationFlowShopSolution* s = new PermutationFlowShopSolution(initial);
    pis.evaluateSolution(*s);
    return s;
}

emili::Solution* emili::pfsp::RMNEH::generate()
{
    std::vector< int > initial = std_start_sequence(pis);
    //initial = rz_improvement_phase(initial,pis);
    initial = nehrs(initial,pis.getNjobs(),pis);

    PermutationFlowShopSolution* s = new PermutationFlowShopSolution(initial);
    pis.evaluateSolution(*s);
    return s;
}

std::vector< int > inline lr_solution_sequence(int start,std::vector< int > u,std::vector< int > initial,emili::pfsp::PermutationFlowShop& pis)
{
    std::vector< float > fndx = lr_index(initial,u,pis);
#ifndef NOC11
    std::sort(u.begin(),u.end(),[fndx](int i1,int i2){return fndx[i1] < fndx[i2];});
#else
    stdc.stds = &fndx;
    std::sort(u.begin(),u.end(),stdc);
#endif
    initial.push_back(u[start]);
    u.erase(u.begin()+start);

    int usize = u.size()-1;

    for(int i=0; i< usize;i++)
    {
        fndx = lr_index(initial,u,pis);
#ifndef NOC11
    std::sort(u.begin(),u.end(),[fndx](int i1,int i2){return fndx[i1] < fndx[i2];});
#else
    stdc.stds = &fndx;
    std::sort(u.begin(),u.end(),stdc);
#endif
        initial.push_back(u[0]);
        u.erase(u.begin());
    }

    initial.push_back(u[0]);
/*    for(int i = 0; i< initial.size() ; i++)
        std::cout << " " << initial[i];
    std::cout << " -> " << initial.size() <<  "\n";*/
    return initial;
}

emili::Solution* emili::pfsp::LRSolution::generate()
{
    std::vector< int > initial;
    initial.push_back(0);
    std::vector< int > u;
    for(int i=1; i<= pis.getNjobs() ; i++)
    {
        u.push_back(i);
    }

    //initial = lr_solution_sequence(0,u,initial,pis);
    std::vector< int > best = lr_solution_sequence(0,u,initial,pis);
    int best_wt = pis.computeObjectiveFunction(best);
    for(int i=1; i < number_of_sequences; i++)
    {
        std::vector< int > temp = lr_solution_sequence(i,u,initial,pis);
        int temp_wt = pis.computeObjectiveFunction(temp);
        if(best_wt>temp_wt)
        {
            best = temp;
            best_wt = temp_wt;
        }
    }
    PermutationFlowShopSolution* s = new PermutationFlowShopSolution(best_wt,best);

    return s;
}

emili::Solution* emili::pfsp::NLRSolution::generate()
{
    std::vector< int > initial;
    initial.push_back(0);
    std::vector< int > u;
    for(int i=1; i<= pis.getNjobs() ; i++)
    {
        u.push_back(i);
    }

    //initial = lr_solution_sequence(0,u,initial,pis);
    std::vector< int > best = lr_solution_sequence(0,u,initial,pis);
    int best_wt = pis.computeObjectiveFunction(best);
    for(int i=1; i < number_of_sequences; i++)
    {
        std::vector< int > temp = lr_solution_sequence(i,u,initial,pis);
        int temp_wt = pis.computeObjectiveFunction(temp);
        if(best_wt>temp_wt)
        {
            best = temp;
            best_wt = temp_wt;
        }
    }

    best = neh2(best,pis.getNjobs(),pis);
    best_wt = pis.computeObjectiveFunction(best);
    PermutationFlowShopSolution* s = new PermutationFlowShopSolution(best_wt,best);

    return s;
}

emili::Solution* emili::pfsp::NRZPerturbation::perturb(Solution *solution)
{

    std::vector< int > removed;
    std::vector< int > solPartial(((emili::pfsp::PermutationFlowShopSolution*)solution)->getJobSchedule());

    int sops = solPartial.size()-1;
    for(int k = 0; k < d; k++) {
        int index = (emili::generateRandomNumber()%sops)+1;
        removed.push_back(solPartial[index]);
        solPartial.erase(solPartial.begin() + index);
        sops--;
    }

    std::vector< int > nrz_sol = rz_seed_sequence(solPartial,removed,prob);
    nrz_sol = neh2(nrz_sol,prob.getNjobs(),prob);

    emili::Solution* s = new emili::pfsp::PermutationFlowShopSolution(nrz_sol);
    prob.evaluateSolution(*s);
    return s;
}

emili::Solution* emili::pfsp::TMIIGPerturbation::perturb(Solution *solution)
{
    int index;
    int min;
    int k,tmp=0,ind=1;
    std::vector< int > removed;
    std::vector< int > solPartial(((emili::pfsp::PermutationFlowShopSolution*)solution)->getJobSchedule());
    //std::cout << "partial size " << solPartial.size() << std::endl;
    int size = solPartial.size();
    std::vector< int > solTMP(size,0); 
    int sops = size-1;
    for(int k = 0; k < d; k++) {
        index = (emili::generateRandomNumber()%sops)+1;
        //std::cout << index << " " ;//<< std::endl;
        removed.push_back(solPartial[index]);
        tblist[solPartial[index]].push_back(solPartial[index-1]);
        if(tblist[solPartial[index]].size() > tbsize)
        {
            tblist[solPartial[index]].erase(tblist[solPartial[index]].begin());
        }
        solPartial.erase(solPartial.begin() + index);
        sops--;
    }

    for(int l=0;l<removed.size();l++){
        sops++;
        k=removed[l];
        min = std::numeric_limits<int>::max();

        for(int r=1; r<sops; r++){
            if(std::find(tblist[k].begin(),tblist[k].end(),solPartial[r]) == tblist[k].end())
            {
            for(int h=1; h<r; h++)
                solTMP[h]=solPartial[h];
            solTMP[r]=k;
            for(int h=r+1; h<=sops; h++)
                solTMP[h]=solPartial[h-1];


            //tmp=compute_total_wt(solTMP,sizePartial+1);
            //                  std::cout << "start perturb" << std::endl;
            //check why plus 1
            tmp = instance.computeObjectiveFunction(solTMP,sops);

            if(tmp<min){
                min=tmp;
                ind=r;
            }
            }
        }

        solPartial.insert(solPartial.begin()+ind,k);

        //std::cout << "end insert " << solPartial.size() << std::endl;
    }
   // assert(min == instance.computeObjectiveFunction(solPartial));
    emili::Solution* s = new emili::pfsp::PermutationFlowShopSolution(min,solPartial);
    //instance.evaluateSolution(*s);
    return s;
}

emili::Solution* emili::pfsp::NWTMIIGPerturbation::perturb(Solution *solution)
{
    int index;
    int min;
    int k,partial_ms=0,ind=1;
    std::vector< int > removed;
    std::vector< int > solPartial(((emili::pfsp::PermutationFlowShopSolution*)solution)->getJobSchedule());
    int makespan = solution->getSolutionValue();    
    int size = solPartial.size();
    std::vector< int > solTMP(size,0);
    int sops = size-1;
    for(int k = 0; k < d; k++) {
        index = (emili::generateRandomNumber()%sops)+1;
        int indexj = solPartial[index];
        int indexmo = solPartial[index-1];
        int indexpo = solPartial[(index+1)%(sops+1)];
        removed.push_back(solPartial[index]);
        tblist[solPartial[index]].push_back(solPartial[index-1]);
        if(tblist[solPartial[index]].size() > tbsize)
        {
            tblist[solPartial[index]].erase(tblist[solPartial[index]].begin());
        }
        makespan = makespan - distances[indexmo][indexj]
                 - distances[indexj][indexpo] + distances[indexmo][indexpo];        
        solPartial.erase(solPartial.begin() + index);
        sops--;
    }    

    for(int l=0;l<removed.size();l++){
        sops++;
        k=removed[l];
        min = std::numeric_limits<int>::max();
        int ok = 0;
        for(int r=1; r<sops+1; r++){
            if(std::find(tblist[k].begin(),tblist[k].end(),solPartial[r]) == tblist[k].end())
            {
                ok++;
            for(int h=1; h<r; h++)
                solTMP[h]=solPartial[h];
            solTMP[r]=k;
            for(int h=r+1; h<=sops; h++)
                solTMP[h]=solPartial[h-1];
            int kmo = solTMP[r-1];
            int kpo = solTMP[(r+1)%(sops+1)];            
            partial_ms = makespan + distances[kmo][k] + distances[k][kpo] - distances[kmo][kpo];                      
            if(partial_ms<min)
            {                
                min = partial_ms;
                ind = r;
            }
            }
        }
        if(ok==0)
        {
            removed.push_back(k);
            sops--;
        }
        else
        {
            solPartial.insert(solPartial.begin()+ind,k);
            makespan = min;
        }
    }
   // assert(min == instance.computeObjectiveFunction(solPartial));
    emili::Solution* s = new emili::pfsp::PermutationFlowShopSolution(min,solPartial);
    //instance.evaluateSolution(*s);
    return s;
}

emili::Solution* emili::pfsp::IGPerturbation::perturb(Solution *solution)
{
    //emili::iteration_increment();

    int index;
    int min;
    int k,tmp=0,ind=1;

    std::vector< int > removed;
    std::vector< int > solPartial(((emili::pfsp::PermutationFlowShopSolution*)solution)->getJobSchedule());
    //std::cout << "partial size " << solPartial.size() << std::endl;
    int size = solPartial.size();
    std::vector< int > solTMP(size,0);
    int sops = size-1;
    for(int k = 0; k < d; k++) {
        index = (emili::generateRandomNumber()%sops)+1;
        //std::cout << index << " " ;//<< std::endl;
        removed.push_back(solPartial[index]);
        solPartial.erase(solPartial.begin() + index);
        sops--;
    }


    //
    // Local search on partial
    //
    for(int l=0;l<removed.size();l++){
        sops++;
        k=removed[l];
        min = std::numeric_limits<int>::max();

        for(int r=1; r<sops; r++){

            for(int h=1; h<r; h++)
                solTMP[h]=solPartial[h];
            solTMP[r]=k;
            for(int h=r+1; h<=sops; h++)
                solTMP[h]=solPartial[h-1];


            //tmp=compute_total_wt(solTMP,sizePartial+1);
            //                  std::cout << "start perturb" << std::endl;
            //check why plus 1
            tmp = instance.computeObjectiveFunction(solTMP,sops);

            if(tmp<min){
                min=tmp;
                ind=r;
            }

        }
        solPartial.insert(solPartial.begin()+ind,k);

        //std::cout << "end insert " << solPartial.size() << std::endl;
    }

    //assert(min == instance.computeObjectiveFunction(solPartial));
    emili::pfsp::PermutationFlowShopSolution* s = new emili::pfsp::PermutationFlowShopSolution(min,solPartial);
    //instance.evaluateSolution(*s);
    //std::vector< std::vector< int > >  etm = std::vector< std::vector< int > >(instance.getNmachines()+1,std::vector<int>(instance.getNjobs()+1,0));
    //clock_t start = clock();
    //int new_value =  instance.computeObjectiveFunction(solPartial,etm,1,instance.getNjobs()+1);

    //s->setSolutionValue(new_value);
    return s;
}

emili::Solution* emili::pfsp::NWIGPerturbation::perturb(Solution *solution)
{
    int index;
    int min;
    int k,partial_ms=0,ind=1;
    std::vector< int > removed;
    std::vector< int > solPartial(((emili::pfsp::PermutationFlowShopSolution*)solution)->getJobSchedule());
    int makespan = solution->getSolutionValue();
    int size = solPartial.size();
    std::vector< int > solTMP(size,0);
    int sops = size-1;
    for(int k = 0; k < d; k++) {
        index = (emili::generateRandomNumber()%sops)+1;
        int indexj = solPartial[index];
        int indexmo = solPartial[index-1];
        int indexpo = solPartial[(index+1)%(sops+1)];
        removed.push_back(indexj);
        makespan = makespan - distances[indexmo][indexj]
                 - distances[indexj][indexpo] + distances[indexmo][indexpo];
        solPartial.erase(solPartial.begin() + index);
        sops--;
    }

    for(int l=0;l<removed.size();l++){
        sops++;
        k=removed[l];
        min = std::numeric_limits<int>::max();
        for(int r=1; r<sops; r++){

            for(int h=1; h<r; h++)
                solTMP[h]=solPartial[h];
            solTMP[r]=k;
            for(int h=r+1; h<=sops; h++)
                solTMP[h]=solPartial[h-1];

            int kmo = solTMP[r-1];
            int kpo = solTMP[(r+1)%(sops+1)];
            partial_ms = makespan + distances[kmo][k] + distances[k][kpo] - distances[kmo][kpo];
            if(partial_ms<min)
            {
                min = partial_ms;
                ind = r;
            }
        }
        solPartial.insert(solPartial.begin()+ind,k);
        makespan = min;
    }
    emili::pfsp::PermutationFlowShopSolution* s = new emili::pfsp::PermutationFlowShopSolution(min,solPartial);
    return s;
}

emili::Solution* emili::pfsp::IGOPerturbation::perturb(Solution *solution)
{
    int index;
    int min;
    int ind=1;
    std::vector< int > removed;
    std::vector< int > solPartial(((emili::pfsp::PermutationFlowShopSolution*)solution)->getJobSchedule());
    int size = solPartial.size();
    std::vector< int > solTMP(size,0);
    int sops = size-1;
    for(int k = 0; k < d; k++) {
        index = (emili::generateRandomNumber()%sops)+1;
        removed.push_back(solPartial[index]);
        solPartial.erase(solPartial.begin() + index);
        sops--;
    }

    for(int l=0;l<removed.size();l++){
        sops++;
        int k=removed[l];
        min = std::numeric_limits<int>::max();
#ifdef ENABLE_SSE
        //std::cout << "pre head" << "\n";
    computeHEAD(solPartial,head,pmatrix,sops-1,nmac);
        //std::cout << "post head" << "\n";
#else
    instance.getInstance().computeHead(solPartial,head,sops);
#endif
        for(int r=1; r<sops; r++){

            for(int h=1; h<r; h++)
                solTMP[h]=solPartial[h];
            solTMP[r]=k;
            for(int h=r+1; h<=sops; h++)
                solTMP[h]=solPartial[h-1];

//            int tmp = instance.computeObjectiveFunction(solTMP,sops);
            int tmp = instance.computeObjectiveFunctionFromHead(solTMP,r,head,sops);
//            std::cout << sops << " <> " << tmp << " - " << tmp1 << std::endl;
//            assert(tmp==tmp1);
            if(tmp<min){
                min=tmp;
                ind=r;
            }

        }
        solPartial.insert(solPartial.begin()+ind,k);
    }

    //assert(min == instance.computeObjectiveFunction(solPartial));
    emili::pfsp::PermutationFlowShopSolution* s = new emili::pfsp::PermutationFlowShopSolution(min,solPartial);
    return s;
}

void emili::pfsp::IGIOPerturbation::updateWeights()
{
    const std::vector < std::vector < long > >& pmat = instance.getProcessingTimesMatrix();
    int nmac = instance.getNmachines();
    int nj = instance.getNjobs();
    for(int i=1; i<= nj; i++)
    {
        //
        weights[i] = pmat[i][1];
        for(int m=2;m<=nmac;m++)
        {
            weights[i] += pmat[i][m];
        }
    }
}

emili::Solution* emili::pfsp::IGIOPerturbation::perturb(Solution *solution)
{
    int index;
    int min;
    int k,tmp=0,ind=1;
    std::vector< int > removed;
    std::vector< int > solPartial(((emili::pfsp::PermutationFlowShopSolution*)solution)->getJobSchedule());
    int size = solPartial.size();
    std::vector< int > solTMP(size,0);
    int sops = size-1;
    for(int k = 0; k < d; k++) {
        index = (emili::generateRandomNumber()%sops)+1;
        removed.push_back(solPartial[index]);
        solPartial.erase(solPartial.begin() + index);
        sops--;
    }

    std::vector < int >& w = this->weights;
#ifndef NOC11
    std::sort(removed.begin(),removed.end(),[w](int i1,int i2){ return w[i1] > w[i2];});
#else
    igioc.stds=&w;
    std::sort(removed.begin(),removed.end(),igioc);
#endif

    for(int l=0;l<removed.size();l++){
        sops++;
        k=removed[l];
        min = std::numeric_limits<int>::max();
        for(int r=1; r<sops; r++){
            for(int h=1; h<r; h++)
                solTMP[h]=solPartial[h];
            solTMP[r]=k;
            for(int h=r+1; h<=sops; h++)
                solTMP[h]=solPartial[h-1];
            tmp = instance.computeObjectiveFunction(solTMP,sops);
            if(tmp<min){
                min=tmp;
                ind=r;
            }

        }
        solPartial.insert(solPartial.begin()+ind,k);
    }

    emili::pfsp::PermutationFlowShopSolution* s = new emili::pfsp::PermutationFlowShopSolution(min,solPartial);
    return s;
}

emili::Solution* emili::pfsp::RSIOPerturbation::perturb(Solution *solution)
{
    int index;
    int min;
    int k,tmp=0,ind=1;
    int nmac = instance.getNmachines();
    std::vector< int > removed;
    std::vector< int > solPartial(((emili::pfsp::PermutationFlowShopSolution*)solution)->getJobSchedule());
    //std::cout << "partial size " << solPartial.size() << std::endl;
    int size = solPartial.size();
    std::vector< int > solTMP(size,0);
    int sops = size-1;
    for(int k = 0; k < d; k++) {
        index = (emili::generateRandomNumber()%sops)+1;
        //std::cout << index << " " ;//<< std::endl;
        removed.push_back(solPartial[index]);
        solPartial.erase(solPartial.begin() + index);
        sops--;
    }

    std::vector < int >& w = this->weights;
#ifndef NOC11
    std::sort(removed.begin(),removed.end(),[w](int i1,int i2){ return w[i1] > w[i2];});
#else
    igioc.stds=&w;
    std::sort(removed.begin(),removed.end(),igioc);
#endif

    for(int l=0;l<removed.size();l++){
        k=removed[l];
        sops++;
        min = std::numeric_limits<int>::max();

#ifdef ENABLE_SSE
        computeHEADandTAIL(solPartial,head,tail,pmatrix,sops-1,nmac);
#else
        instance.computeTAmatrices(solPartial,head,tail,sops);
#endif
        for(int r=1; r<sops; r++){

            for(int h=1; h<r; h++)
                solTMP[h]=solPartial[h];
            solTMP[r]=k;
            for(int h=r+1; h<=sops; h++)
                solTMP[h]=solPartial[h-1];


            //tmp=compute_total_wt(solTMP,sizePartial+1);
            //                  std::cout << "start perturb" << std::endl;
            //check why plus 1
            long int c_cur = head[1][r-1]+pmatrix[k][1];
            long int c_max = c_cur+tail[1][r];
            for (int i = 2; i <= nmac; ++i) {
                int c_pm = head[i][r-1];
                if(c_pm > c_cur)
                {
                    c_cur = c_pm;
                }
                c_cur = c_cur + pmatrix[k][i];
                long int c_can = (c_cur+tail[i][r]);
                c_max = c_max>c_can?c_max:c_can;
            }

            tmp = c_max;
            if(tmp<min){
                min=tmp;
                ind=r;
            }

        }
        solPartial.insert(solPartial.begin()+ind,k);

        //std::cout << "end insert " << solPartial.size() << std::endl;
    }

    //assert(min == instance.computeObjectiveFunction(solPartial));
    emili::pfsp::PermutationFlowShopSolution* s = new emili::pfsp::PermutationFlowShopSolution(min,solPartial);
    //instance.evaluateSolution(*s);
    return s;
}

emili::Solution* emili::pfsp::RSPerturbation::perturb(Solution *solution)
{
    int index;
    int min;
    int k,tmp=0,ind=1;
    int nmac = instance.getNmachines();
    std::vector< int > removed;
    std::vector< int > solPartial(((emili::pfsp::PermutationFlowShopSolution*)solution)->getJobSchedule());
    //std::cout << "partial size " << solPartial.size() << std::endl;
    int size = solPartial.size();
    std::vector< int > solTMP(size,0);
    int sops = size-1;
    for(int k = 0; k < d; k++) {
        index = (emili::generateRandomNumber()%sops)+1;
        //std::cout << index << " " ;//<< std::endl;
        removed.push_back(solPartial[index]);
        solPartial.erase(solPartial.begin() + index);
        sops--;
    }

    for(int l=0;l<removed.size();l++){
        sops++;
        k=removed[l];
        min = std::numeric_limits<int>::max();
#ifdef ENABLE_SSE
        computeHEADandTAIL(solPartial,head,tail,pmatrix,sops-1,nmac);
#else
        instance.computeTAmatrices(solPartial,head,tail,sops);
#endif
        for(int r=1; r<sops; r++){

            for(int h=1; h<r; h++)
                solTMP[h]=solPartial[h];
            solTMP[r]=k;
            for(int h=r+1; h<=sops; h++)
                solTMP[h]=solPartial[h-1];


            //tmp=compute_total_wt(solTMP,sizePartial+1);
            //                  std::cout << "start perturb" << std::endl;
            //check why plus 1
            long int c_cur = head[1][r-1]+pmatrix[k][1];
            long int c_max = c_cur+tail[1][r];
            for (int i = 2; i <= nmac; ++i) {
                int c_pm = head[i][r-1];
                if(c_pm > c_cur)
                {
                    c_cur = c_pm;
                }
                c_cur = c_cur + pmatrix[k][i];
                long int c_can = (c_cur+tail[i][r]);
                c_max = c_max>c_can?c_max:c_can;
            }
            tmp = c_max;

            if(tmp<min){
                min=tmp;
                ind=r;
            }

        }
        solPartial.insert(solPartial.begin()+ind,k);
        //std::cout << "end insert " << solPartial.size() << std::endl;
    }

    //assert(min == instance.computeObjectiveFunction(solPartial));
    emili::pfsp::PermutationFlowShopSolution* s = new emili::pfsp::PermutationFlowShopSolution(min,solPartial);
    //instance.evaluateSolution(*s);
    return s;
}

emili::Solution* emili::pfsp::RSffPerturbation::perturb(Solution *solution)
{
    int index;
    int min;
    int k,tmp=0,ind=1;
    int nmac = instance.getNmachines();
    std::vector< int > removed;
    std::vector< int > solPartial(((emili::pfsp::PermutationFlowShopSolution*)solution)->getJobSchedule());
    //std::cout << "partial size " << solPartial.size() << std::endl;
    int size = solPartial.size();
    std::vector< int > solTMP(size,0);
    int sops = size-1;
    for(int k = 0; k < d; k++) {
        index = (emili::generateRandomNumber()%sops)+1;
        //std::cout << index << " " ;//<< std::endl;
        removed.push_back(solPartial[index]);
        solPartial.erase(solPartial.begin() + index);
        sops--;
    }

    for(int l=0;l<removed.size();l++){
        k=removed[l];
        sops++;
        min = std::numeric_limits<int>::max();
#ifdef ENABLE_SSE
        computeHEADandTAIL(solPartial,head,tail,pmatrix,sops-1,nmac);
#else
        instance.computeTAmatrices(solPartial,head,tail,sops);
#endif
         std::vector< int >  ptb;
        for(int r=1; r<sops; r++){

            for(int h=1; h<r; h++)
                solTMP[h]=solPartial[h];
            solTMP[r]=k;
            for(int h=r+1; h<=sops; h++)
                solTMP[h]=solPartial[h-1];


            //tmp=compute_total_wt(solTMP,sizePartial+1);
            //                  std::cout << "start perturb" << std::endl;
            //check why plus 1
            long int c_cur = head[1][r-1]+pmatrix[k][1];
            long int c_max = c_cur+tail[1][r];
            for (int i = 2; i <= nmac; ++i) {
                int c_pm = head[i][r-1];
                if(c_pm > c_cur)
                {
                    c_cur = c_pm;
                }
                c_cur = c_cur + pmatrix[k][i];
                long int c_can = (c_cur+tail[i][r]);
                c_max = c_max>c_can?c_max:c_can;
            }
            tmp = c_max;//instance.computeObjectiveFunction(solTMP,sizePartial);

            //assert(c_max == tmp);
            if(tmp<min){
                min=tmp;
                ind=r;
                ptb.clear();
                ptb.push_back(ind);
            }else if(tmp==min)
            {
                ptb.push_back(r);
            }

        }
        /** TIE BREAKING FF*/
        int tb = ptb.size();

            if(tb > 1 && l<sops)
            {
                //tie breaker!
                int bp = ptb[0];
                int itbp = std::numeric_limits<int> ::max();
                for(int l2=0; l2 < tb; l2++)
                {
                    int ptbl = ptb[l2];
                    int itdd = 0;
                    int fil = 0;
                    if(ptbl == l)
                    {
                        fil =  head[1][ptbl-1]+pmatrix[k][1]; //???????????
                        for(int i=2 ; i <= nmac; i++)
                        {
                            itdd += fil-head[i][ptbl-1]-pmatrix[k][i];
                            fil = std::max(head[i][ptbl-1],fil) + pmatrix[k][i];
                        }
                    }
                    else
                    {
                        int sptbl = solTMP[ptbl];
                        int filp =  fil+ pmatrix[sptbl][1];
                        for(int i=2 ; i <= nmac; i++)
                        {
                            fil = std::max(head[i][ptbl-1],fil) + pmatrix[k][i];
                            itdd += fil-head[i][ptbl-1]-pmatrix[k][i]+pmatrix[sptbl][i]+std::max(filp-fil,0);
                            filp = std::max(fil,filp) + pmatrix[sptbl][i];
                        }
                    }
                    if(itbp > itdd)
                    {
                        itbp = itdd;
                        bp = ptbl;
                    }
                }
                ind = bp;
            }
        /** END TIE BREAKING*/
        solPartial.insert(solPartial.begin()+ind,k);        
        //std::cout << "end insert " << solPartial.size() << std::endl;
    }
    emili::pfsp::PermutationFlowShopSolution* s = new emili::pfsp::PermutationFlowShopSolution(min,solPartial);
    //instance.evaluateSolution(*s);
    return s;
}

emili::Solution* emili::pfsp::IgLsPerturbation::perturb(Solution *solution)
{
    //emili::iteration_increment();
   // assert(solution->getSolutionValue() >= 0);
    int index;
    int min;
    int k,tmp=0,ind=1;

    std::vector< int > removed;
    std::vector< int > solPartial(((emili::pfsp::PermutationFlowShopSolution*)solution)->getJobSchedule());
    //std::cout << "partial size " << solPartial.size() << std::endl;
    int size = solPartial.size();
    std::vector< int > solTMP(size,0);
    int sops = size-1;
    for(int k = 0; k < d; k++) {
        index = (emili::generateRandomNumber()%sops)+1;
        //std::cout << index << " " ;//<< std::endl;
        removed.push_back(solPartial[index]);
        solPartial.erase(solPartial.begin() + index);
        sops--;
    }


    //
    // Local search on partial
    //
    emili::pfsp::PermutationFlowShopSolution s(solPartial);
    s.setSolutionValue(instance.computeObjectiveFunction(solPartial,sops));
#ifdef WITH_STATS
    bool print_state = get_print();
    set_print(false);
#endif
    emili::pfsp::PermutationFlowShopSolution* s_n =(emili::pfsp::PermutationFlowShopSolution*) ls->search(&s);
#ifdef WITH_STATS
    set_print(print_state);
#endif


    solPartial = s_n->getJobSchedule();
    //delete s_n;
    for(int l=0;l<removed.size();l++){
        sops++;
        k=removed[l];
        min = std::numeric_limits<int>::max();
        for(int r=1; r<sops; r++){

            for(int h=1; h<r; h++)
                solTMP[h]=solPartial[h];
            solTMP[r]=k;
            for(int h=r+1; h<=sops; h++)
                solTMP[h]=solPartial[h-1];


            //tmp=compute_total_wt(solTMP,sizePartial+1);
            //                  std::cout << "start perturb" << std::endl;
            //check why plus 1
            tmp = instance.computeObjectiveFunction(solTMP,sops);

            if(tmp<min){
                min=tmp;
                ind=r;
            }

        }
        solPartial.insert(solPartial.begin()+ind,k);
        //sops++;
        //std::cout << "end insert " << solPartial.size() << std::endl;
    }

    //assert(min == instance.computeObjectiveFunction(solPartial));
    s_n->setJobSchedule(solPartial);

   // instance.evaluateSolution(*s_n);
    //assert(min==s_n->getSolutionValue());
    s_n->setSolutionValue(min);
    return s_n;
}

emili::Solution* emili::pfsp::NwIgLsPerturbation::perturb(Solution *solution)
{
    //emili::iteration_increment();
   // assert(solution->getSolutionValue() >= 0);
    int index;
    int min;
    int k,partial_ms=0,ind=1;

    std::vector< int > removed;
    std::vector< int > solPartial(((emili::pfsp::PermutationFlowShopSolution*)solution)->getJobSchedule());
    //std::cout << "partial size " << solPartial.size() << std::endl;
    int size = solPartial.size();
    std::vector< int > solTMP(size,0);
    int sops = size-1;
    for(int k = 0; k < d; k++) {
        index = (emili::generateRandomNumber()%sops)+1;
        //std::cout << index << " " ;//<< std::endl;
        removed.push_back(solPartial[index]);
        solPartial.erase(solPartial.begin() + index);
        sops--;
    }


    //
    // Local search on partial
    //
    emili::pfsp::PermutationFlowShopSolution s(solPartial);
    s.setSolutionValue(instance.computeObjectiveFunction(solPartial,sops));
#ifdef WITH_STATS
    bool print_state = get_print();
    set_print(false);
#endif
    emili::pfsp::PermutationFlowShopSolution* s_n =(emili::pfsp::PermutationFlowShopSolution*) ls->search(&s);
#ifdef WITH_STATS
    set_print(print_state);
#endif


    solPartial = s_n->getJobSchedule();
    int makespan = s_n->getSolutionValue();
    //delete s_n;
    for(int l=0;l<removed.size();l++){
        sops++;
        k=removed[l];
        min = std::numeric_limits<int>::max();
        for(int r=1; r<sops; r++){

            for(int h=1; h<r; h++)
                solTMP[h]=solPartial[h];
            solTMP[r]=k;
            for(int h=r+1; h<=sops; h++)
                solTMP[h]=solPartial[h-1];


            int kmo = solTMP[r-1];
            int kpo = solTMP[(r+1)%(sops+1)];
            partial_ms = makespan + distances[kmo][k] + distances[k][kpo] - distances[kmo][kpo];
            if(partial_ms<min)
            {
                min = partial_ms;
                ind = r;
            }

        }
        solPartial.insert(solPartial.begin()+ind,k);
        makespan = min;
        //sops++;
        //std::cout << "end insert " << solPartial.size() << std::endl;
    }

    //assert(min == instance.computeObjectiveFunction(solPartial));
    s_n->setJobSchedule(solPartial);
    //instance.evaluateSolution(*s_n);
    //assert(min==this->instance.computeObjectiveFunction(solPartial));
    s_n->setSolutionValue(min);
    return s_n;
}

emili::Solution* emili::pfsp::IGOLsPerturbation::perturb(Solution *solution)
{
    int index;
    int min;
    int ind=1;
    std::vector< int > removed;
    std::vector< int > solPartial(((emili::pfsp::PermutationFlowShopSolution*)solution)->getJobSchedule());
    int size = solPartial.size();
    std::vector< int > solTMP(size,0);
    int sops = size-1;
    for(int k = 0; k < d; k++) {
        index = (emili::generateRandomNumber()%sops)+1;
        removed.push_back(solPartial[index]);
        solPartial.erase(solPartial.begin() + index);
        sops--;
    }

    //
    // Local search on partial
    //
    emili::pfsp::PermutationFlowShopSolution s(solPartial);
    s.setSolutionValue(instance.computeObjectiveFunction(solPartial,sops));
#ifdef WITH_STATS
    bool print_state = get_print();
    set_print(false);
#endif
    emili::pfsp::PermutationFlowShopSolution* s_n =(emili::pfsp::PermutationFlowShopSolution*) ls->search(&s);
#ifdef WITH_STATS
    set_print(print_state);
#endif


    solPartial = s_n->getJobSchedule();

    for(int l=0;l<removed.size();l++){
        sops++;
        int k=removed[l];
        min = std::numeric_limits<int>::max();
#ifdef ENABLE_SSE
        //std::cout << "pre head" << "\n";
    computeHEAD(solPartial,head,pmatrix,sops-1,nmac);
        //std::cout << "post head" << "\n";
#else
    instance.getInstance().computeHead(solPartial,head,sops);
#endif
        for(int r=1; r<sops; r++){

            for(int h=1; h<r; h++)
                solTMP[h]=solPartial[h];
            solTMP[r]=k;
            for(int h=r+1; h<=sops; h++)
                solTMP[h]=solPartial[h-1];

//            int tmp = instance.computeObjectiveFunction(solTMP,sops);
            int tmp = instance.computeObjectiveFunctionFromHead(solTMP,r,head,sops);
//            std::cout << sops << " <> " << tmp << " - " << tmp1 << std::endl;
//            assert(tmp==tmp1);
            if(tmp<min){
                min=tmp;
                ind=r;
            }

        }
        solPartial.insert(solPartial.begin()+ind,k);
    }

    //assert(min == instance.computeObjectiveFunction(solPartial));
    //emili::pfsp::PermutationFlowShopSolution* ns = new emili::pfsp::PermutationFlowShopSolution(min,solPartial);
    //return ns;
    // delete s;
     s_n->setJobSchedule(solPartial);

    // instance.evaluateSolution(*s_n);
     //assert(min==s_n->getSolutionValue());
     s_n->setSolutionValue(min);
     return s_n;
}

emili::Solution* emili::pfsp::RSLSPerturbation::perturb(Solution *solution)
{
    //emili::iteration_increment();

    int index;
    int min;
    int k,tmp=0,ind=1;

    std::vector< int > removed;
    std::vector< int > solPartial(((emili::pfsp::PermutationFlowShopSolution*)solution)->getJobSchedule());
    //std::cout << "partial size " << solPartial.size() << std::endl;
    int size = solPartial.size();
    std::vector< int > solTMP(size,0);

    int sops = size-1;
    for(int k = 0; k < d; k++) {
        index = (emili::generateRandomNumber()%sops)+1;
        //std::cout << index << " " ;//<< std::endl;
        removed.push_back(solPartial[index]);
        solPartial.erase(solPartial.begin() + index);
        sops--;
    }

    //
    // Local search on partial
    //
    emili::pfsp::PermutationFlowShopSolution s(solPartial);
    s.setSolutionValue(instance.computeObjectiveFunction(solPartial,sops));
#ifdef WITH_STATS
    bool print_state = get_print();
    set_print(false);
#endif
    emili::pfsp::PermutationFlowShopSolution* s_n =(emili::pfsp::PermutationFlowShopSolution*) ls->search(&s);
#ifdef WITH_STATS
    set_print(print_state);
#endif
    solPartial = s_n->getJobSchedule();
    int mac = instance.getNmachines();
    for(int l=0;l<removed.size();l++){
        sops++;
        k=removed[l];
        min = std::numeric_limits<int>::max();
#ifdef ENABLE_SSE
        computeHEADandTAIL(solPartial,head,tail,pmatrix,sops-1,mac);
#else
        instance.computeTAmatrices(solPartial,head,tail,sops);
#endif
        for(int r=1; r<sops; r++){

            for(int h=1; h<r; h++)
                solTMP[h]=solPartial[h];
            solTMP[r]=k;
            for(int h=r+1; h<=sops; h++)
                solTMP[h]=solPartial[h-1];


            //tmp=compute_total_wt(solTMP,sizePartial+1);
            //                  std::cout << "start perturb" << std::endl;
            //check why plus 1
            //tmp = instance.computeObjectiveFunction(solTMP,sizePartial);
            long int c_cur = head[1][r-1]+pmatrix[k][1];
            long int c_max = c_cur+tail[1][r];
            for (int i = 2; i <= mac; ++i) {
                int c_pm = head[i][r-1];
                if(c_pm > c_cur)
                {
                    c_cur = c_pm;
                }
                c_cur = c_cur + pmatrix[k][i];
                long int c_can = (c_cur+tail[i][r]);
                c_max = c_max>c_can?c_max:c_can;
            }
            tmp = c_max;
            //assert(c_max == instance.computeMS(solTMP,sops));
            if(tmp<min){
                min=tmp;
                ind=r;
            }

        }
        solPartial.insert(solPartial.begin()+ind,k);
        //sops++;
        //std::cout << "end insert " << solPartial.size() << std::endl;

    }

   // delete s;
    s_n->setJobSchedule(solPartial);

   // instance.evaluateSolution(*s_n);
    //assert(min==s_n->getSolutionValue());
    s_n->setSolutionValue(min);
    return s_n;
}

emili::Solution* emili::pfsp::RSffLSPerturbation::perturb(Solution *solution)
{
    int index;
    int min;
    int k,tmp=0,ind=1;
    int nmac = instance.getNmachines();
    std::vector< int > removed;
    std::vector< int > solPartial(((emili::pfsp::PermutationFlowShopSolution*)solution)->getJobSchedule());
    //std::cout << "partial size " << solPartial.size() << std::endl;
    int size = solPartial.size();
    std::vector< int > solTMP(size,0);

    int sops = size-1;
    for(int k = 0; k < d; k++) {
        index = (emili::generateRandomNumber()%sops)+1;
        //std::cout << index << " " ;//<< std::endl;
        removed.push_back(solPartial[index]);
        solPartial.erase(solPartial.begin() + index);
        sops--;
    }

    //
    // Local search on partial
    //
    emili::pfsp::PermutationFlowShopSolution s(solPartial);
    s.setSolutionValue(instance.computeObjectiveFunction(solPartial,sops));
#ifdef WITH_STATS
    bool print_state = get_print();
    set_print(false);
#endif
    emili::pfsp::PermutationFlowShopSolution* s_n =(emili::pfsp::PermutationFlowShopSolution*) ls->search(&s);
#ifdef WITH_STATS
    set_print(print_state);
#endif
    solPartial = s_n->getJobSchedule();

    for(int l=0;l<removed.size();l++){
        k=removed[l];
        sops++;
        min = std::numeric_limits<int>::max();

#ifdef ENABLE_SSE
        computeHEADandTAIL(solPartial,head,tail,pmatrix,sops-1,nmac);
#else
        instance.computeTAmatrices(solPartial,head,tail,sops);
#endif
         std::vector< int >  ptb;
        for(int r=1; r<sops; r++){

            for(int h=1; h<r; h++)
                solTMP[h]=solPartial[h];
            solTMP[r]=k;
            for(int h=r+1; h<=sops; h++)
                solTMP[h]=solPartial[h-1];


            //tmp=compute_total_wt(solTMP,sizePartial+1);
            //                  std::cout << "start perturb" << std::endl;
            //check why plus 1
            long int c_cur = head[1][r-1]+pmatrix[k][1];
            long int c_max = c_cur+tail[1][r];
            for (int i = 2; i <= nmac; ++i) {
                int c_pm = head[i][r-1];
                if(c_pm < c_cur)
                {
                    c_cur = c_cur + pmatrix[k][i];
                }
                else
                {
                    c_cur = c_pm + pmatrix[k][i];
                }
                long int c_can = (c_cur+tail[i][r]);
                c_max = c_max>c_can?c_max:c_can;
            }
            tmp = c_max;//instance.computeObjectiveFunction(solTMP,sizePartial);

            //assert(c_max == tmp);
            if(tmp<min){
                min=tmp;
                ind=r;
                ptb.clear();
                ptb.push_back(ind);
            }else if(tmp==min)
            {
                ptb.push_back(r);
            }

        }
        /** TIE BREAKING FF*/
        int tb = ptb.size();

            if(tb > 1 && l<sops)
            {
                //tie breaker!
                int bp = ptb[0];
                int itbp = std::numeric_limits<int> ::max();
                for(int l2=0; l2 < tb; l2++)
                {
                    int ptbl = ptb[l2];
                    int itdd = 0;
                    int fil = 0;
                    if(ptbl == l)
                    {
                        fil =  head[1][ptbl-1]+pmatrix[k][1]; //???????????
                        for(int i=2 ; i <= nmac; i++)
                        {
                            itdd += fil-head[i][ptbl-1]-pmatrix[k][i];
                            fil = std::max(head[i][ptbl-1],fil) + pmatrix[k][i];
                        }
                    }
                    else
                    {
                        int sptbl = solTMP[ptbl];
                        int filp =  fil+ pmatrix[sptbl][1];
                        for(int i=2 ; i <= nmac; i++)
                        {
                            fil = std::max(head[i][ptbl-1],fil) + pmatrix[k][i];
                            itdd += fil-head[i][ptbl-1]-pmatrix[k][i]+pmatrix[sptbl][i]+std::max(filp-fil,0);
                            filp = std::max(fil,filp) + pmatrix[sptbl][i];
                        }
                    }
                    if(itbp > itdd)
                    {
                        itbp = itdd;
                        bp = ptbl;
                    }
                }
                ind = bp;
            }
        /** END TIE BREAKING*/
        solPartial.insert(solPartial.begin()+ind,k);
        //std::cout << "end insert " << solPartial.size() << std::endl;
    }

    s_n->setJobSchedule(solPartial);

   // instance.evaluateSolution(*s_n);
    //assert(min==s_n->getSolutionValue());
    s_n->setSolutionValue(min);
    return s_n;
}

emili::Solution* emili::pfsp::RestartPerturbation::perturb(Solution *solution)
{
    emili::Solution* best;
    emili::Solution* current;
    best = initial->generateSolution();

    for(int i=0;i<num_of_solutions;i++)
    {
        current = initial->generateSolution();
        if( *current < *best)
        {
            delete best;
            best = current;
        }
        else
        {
            delete current;
        }

    }

    if(locser)
    {
        current = ls->search(best);
        if(current!=best)
            delete best;

        return current;
    }

    return best;
}

emili::Solution* emili::pfsp::MPTLMPerturbation::perturb(Solution *solution)
{
    std::vector<int> sol = ((PermutationFlowShopSolution*)solution)->getJobSchedule();
    std::vector<int> best = sol;
    std::vector<int> current = sol;
    int bms = nwnehrs(best,njobs,pis,distances);
    int cms = 0;
    for(int i=0;i<num_of_solutions;i++)
    {
        //current = sol;
        cms = nwnehrs(current,njobs,pis,distances);
        //cms = pis.computeObjectiveFunction(current);
       // std::cout << i+1 << "] current " << cms << " best " << bms << std::endl;
        if( cms < bms)
        {
            best = current;
            bms = cms;
        }

    }
    //std::cout << "final " << bms << "\n";
    //assert(bms = pis.computeObjectiveFunction(best));
    return new PermutationFlowShopSolution(bms,best);
}

emili::pfsp::RestartPerturbation::~RestartPerturbation()
{
    if(ls!=nullptr)
        delete ls;
    delete initial;
}


/**
Construct the solution inserting one job at a time, by always selecting
the one that minimizes the weighted earlyness.

The weighted earlyness of job Ji is computed as wi · (di − Ci).

Note: the solution is constructed incrementally, and at each iteration
Ci corresponds to the makespan of the partial solution

void emili::pfsp::PfspSlackInitialSolution::slackInitial(std::vector<int> & sol)
{
    int noj = pis.getNjobs();
    int Ci  = 0;
    std::vector<bool> assigned(noj+1, false);
    std::vector<int> partial(noj+1,0);
    for (int var = 0; var < noj; ++var) {
        int minJ = 0 ;
        int minE = 4*10e9;
        for (int jb = 1; jb <= noj; ++jb) {
            if(!assigned[jb]){
                int dd = pis.getDueDate(jb);
                int pr = pis.getPriority(jb);
                int wej = pr * (dd-Ci);
                if(minE > wej){
                    minJ = jb;
                    minE = wej;
                }
            }
        }
        partial[var+1] = minJ;
        assigned[minJ] = true;
        //Ci = pis.computeMS(partial);
        //cout << "partial makespan: " << Ci << std::endl;
    }
    sol = partial;

}

*/

emili::Solution* emili::pfsp::PfspNeighborhood::step(emili::Solution *currentSolution)
{
    return computeStep(currentSolution);
}

int emili::pfsp::PfspNeighborhood::size()
{
    return (pis.getNjobs()*(pis.getNjobs()-1))/2;
}

void emili::pfsp::PfspNeighborhood::reset()
{
    /** No counters to reset*/
}

emili::Neighborhood::NeighborhoodIterator emili::pfsp::PfspInsertNeighborhood::begin(emili::Solution *base)
{
    ep_iterations = 1;
    sp_iterations = 1;
    return emili::Neighborhood::NeighborhoodIterator(this,base);
}

emili::Neighborhood::NeighborhoodIterator emili::pfsp::TaillardAcceleratedInsertNeighborhood::begin(emili::Solution *base)
{
    ep_iterations = 1;
    sp_iterations = 1;
    std::vector< int > sol(((emili::pfsp::PermutationFlowShopSolution*)base)->getJobSchedule());
    sol.erase(sol.begin()+start_position);
#ifdef ENABLE_SSE
            computeHEADandTAIL(sol,head,tail,pmatrix,njobs-1,nmac);
#else
            computeTAmatrices(sol);
#endif
    return emili::Neighborhood::NeighborhoodIterator(this,base);
}

emili::Neighborhood::NeighborhoodIterator emili::pfsp::FSTaillardAcceleratedInsertNeighborhood::begin(emili::Solution *base)
{
    ep_iterations = 1;
    sp_iterations = 1;
    std::vector< int > sol(((emili::pfsp::PermutationFlowShopSolution*)base)->getJobSchedule());
    sol.erase(sol.begin()+start_position);
#ifdef ENABLE_SSE
            computeHEADandTAIL(sol,head,tail,pmatrix,njobs-1,nmac);
#else
            computeTAmatrices(sol);
#endif
    current_value = base->getSolutionValue();
    improved = false;
    return emili::Neighborhood::NeighborhoodIterator(this,base);
}

emili::Neighborhood::NeighborhoodIterator emili::pfsp::HeavilyApproximatedTaillardAcceleratedInsertNeighborhood::begin(Solution *base)
{
    ep_iterations = 1;
    sp_iterations = 1;
    std::vector< int > sol(((emili::pfsp::PermutationFlowShopSolution*)base)->getJobSchedule());
    sol.erase(sol.begin()+start_position);
#ifdef ENABLE_SSE
    computeHEAD(sol,head,pmatrix,njobs-1,nmac);
#else
    computeHead(sol);    
#endif
    return emili::Neighborhood::NeighborhoodIterator(this,base);
}

emili::Neighborhood::NeighborhoodIterator emili::pfsp::Natx2Neighborhood::begin(Solution *base)
{
    ep_iterations = 1;
    sp_iterations = 1;
    std::vector< int > sol(((emili::pfsp::PermutationFlowShopSolution*)base)->getJobSchedule());
    sol.erase(sol.begin()+start_position);
#ifdef ENABLE_SSE
    computeHEAD(sol,head,pmatrix,njobs-1,nmac);
#else
    computeHead(sol);
#endif
    value_wt = base->getSolutionValue();
    return emili::Neighborhood::NeighborhoodIterator(this,base);
}

emili::Neighborhood::NeighborhoodIterator emili::pfsp::NrzTCTNeighborhood::begin(Solution *base)
{
    ep_iterations = 1;
    sp_iterations = 1;
    std::vector< int > sol(((emili::pfsp::PermutationFlowShopSolution*)base)->getJobSchedule());
    this->seed_seq = sol;
    return emili::Neighborhood::NeighborhoodIterator(this,base);
}

emili::Neighborhood::NeighborhoodIterator emili::pfsp::NoIdleAcceleratedInsertNeighborhood::begin(Solution *base)
{
    ep_iterations = 1;
    sp_iterations = 1;
    std::vector< int > sol(((emili::pfsp::PermutationFlowShopSolution*)base)->getJobSchedule());
    sol.erase(sol.begin()+start_position);
    pis.computeNoIdleTAmatrices(sol,head,tail);
    return emili::Neighborhood::NeighborhoodIterator(this,base);
}

/*emili::Neighborhood::NeighborhoodIterator emili::pfsp::NoWaitAcceleratedNeighborhood::begin(Solution *base)
{
    ep_iterations = 1;
    sp_iterations = 1;
    computeNoWaitTimeDistances(base);
    return emili::Neighborhood::NeighborhoodIterator(this,base);
}*/

emili::Neighborhood::NeighborhoodIterator emili::pfsp::TAxInsertNeighborhood::begin(emili::Solution *base)
{
    ep_iterations = 1;
    sp_iterations = 1;
    std::vector< int > sol(((emili::pfsp::PermutationFlowShopSolution*)base)->getJobSchedule());
    sol.erase(sol.begin()+start_position);
    pis.computeTAmatrices(sol,head,tails[njobs]);
    pis.computeTails(sol,tails);

    return emili::Neighborhood::NeighborhoodIterator(this,base);
}

emili::Neighborhood::NeighborhoodIterator emili::pfsp::PfspExchangeNeighborhood::begin(emili::Solution *base)
{
    ep_iterations = 1;
    sp_iterations = 1;
    return emili::Neighborhood::NeighborhoodIterator(this,base);
}

emili::Neighborhood::NeighborhoodIterator emili::pfsp::AxtExchange::begin(Solution *base)
{
    ep_iterations = 1;
    sp_iterations = 1;
    std::vector< int > sol(((emili::pfsp::PermutationFlowShopSolution*)base)->getJobSchedule());
    sol.erase(sol.begin()+start_position);
#ifdef ENABLE_SSE
    computeHEAD(sol,head,pmatrix,njobs-1,nmac);
#else
    computeHead(sol);
#endif
    return emili::Neighborhood::NeighborhoodIterator(this,base);
}

emili::Neighborhood::NeighborhoodIterator emili::pfsp::PfspTransposeNeighborhood::begin(emili::Solution *base)
{    
    sp_iterations = 1;
    return emili::Neighborhood::NeighborhoodIterator(this,base);
}

emili::Neighborhood::NeighborhoodIterator emili::pfsp::XTransposeNeighborhood::begin(emili::Solution *base)
{
    last_saved_position = -1;
    return emili::Neighborhood::NeighborhoodIterator(this,base);
}
emili::Solution* emili::pfsp::PfspBackwardInsertNeighborhood::computeStep(emili::Solution* value)
{
    emili::iteration_increment();
    if(sp_iterations > njobs)
    {
        return nullptr;
    }
    else{
        end_position = ((end_position)%njobs)+1;
        if(ep_iterations < njobs){
            ep_iterations++;
           /**  if(ep_iterations == sp_iterations){
                ep_iterations++;
                end_position++;
                std::cout << "BOOM!" << std::endl;
            }*/
            if(end_position == start_position-1)
            {
                end_position=((end_position+1)%njobs)+1;
                ep_iterations+=2;
                if(ep_iterations > njobs && sp_iterations+1 > njobs)
                    return nullptr;
            }

        }
        else
        {
            sp_iterations++;
            ep_iterations = 1;
            start_position = ((start_position)%njobs)+1;
            if(end_position == start_position-1)
            {
                end_position=((end_position+1)%njobs)+1;
                ep_iterations+=2;
                if(ep_iterations > njobs && sp_iterations+1 > njobs)
                    return nullptr;
            }

        }

        std::vector < int >& newsol = ((emili::pfsp::PermutationFlowShopSolution*)value)->getJobSchedule();
        int sol_i = newsol[start_position];
        newsol.erase(newsol.begin()+start_position);
        newsol.insert(newsol.begin()+end_position,sol_i);
        long int new_value = pis.computeObjectiveFunction(newsol);
        value->setSolutionValue(new_value);
        return value;
    }
}

emili::Solution* emili::pfsp::PfspForwardInsertNeighborhood::computeStep(emili::Solution* value)
{
    emili::iteration_increment();
    if(sp_iterations > njobs)
    {
        return nullptr;
    }
    else{
        end_position = ((end_position)%njobs)+1;
        if(ep_iterations < njobs){
            ep_iterations++;
           /**  if(ep_iterations == sp_iterations){
                ep_iterations++;
                end_position++;
                std::cout << "BOOM!" << std::endl;
            }*/
            if(end_position == start_position-1)
            {
                end_position=((end_position+1)%njobs)+1;
                ep_iterations+=2;
                if(ep_iterations > njobs && sp_iterations+1 > njobs)
                    return nullptr;
            }

        }
        else
        {
            sp_iterations++;
            ep_iterations = 1;
            start_position = ((start_position)%njobs)+1;
            if(end_position == start_position-1)
            {
                end_position=((end_position+1)%njobs)+1;
                ep_iterations+=2;
                if(ep_iterations > njobs && sp_iterations+1 > njobs)
                    return nullptr;
            }

        }
        std::vector < int >& newsol = ((emili::pfsp::PermutationFlowShopSolution*)value)->getJobSchedule();
        int sol_i = newsol[start_position];
        newsol.erase(newsol.begin()+start_position);
        newsol.insert(newsol.begin()+end_position,sol_i);
        long int new_value = pis.computeObjectiveFunction(newsol);
        value->setSolutionValue(new_value);
        return value;
    }

}

emili::Solution* emili::pfsp::PfspInsertNeighborhood::computeStep(emili::Solution* value)
{
    emili::iteration_increment();
    if(sp_iterations > njobs)
    {
        return nullptr;
    }
    else{
        end_position = ((end_position)%njobs)+1;
        if(ep_iterations < njobs){
            ep_iterations++;

            if(end_position == start_position-1)
            {
                end_position=((end_position+1)%njobs)+1;
                ep_iterations+=2;
                if(ep_iterations > njobs && sp_iterations+1 > njobs)
                    return nullptr;
            }

        }
        else
        {
            sp_iterations++;
            ep_iterations = 1;
            start_position = ((start_position)%njobs)+1;
            if(end_position == start_position-1)
            {
                end_position=((end_position+1)%njobs)+1;
                ep_iterations+=2;
                if(ep_iterations > njobs && sp_iterations+1 > njobs)
                    return nullptr;
            }

        }
        std::vector < int >& newsol = ((emili::pfsp::PermutationFlowShopSolution*)value)->getJobSchedule();
        int sol_i = newsol[start_position];
        newsol.erase(newsol.begin()+start_position);
        newsol.insert(newsol.begin()+end_position,sol_i);
        long int new_value = pis.computeObjectiveFunction(newsol);    
        value->setSolutionValue(new_value);
        return value;
    }
}

void emili::pfsp::PfspInsertNeighborhood::reverseLastMove(Solution *step)
{
    std::vector < int >& newsol = ((emili::pfsp::PermutationFlowShopSolution*)step)->getJobSchedule();
    int sol_i = newsol[end_position];
    newsol.erase(newsol.begin()+end_position);
    newsol.insert(newsol.begin()+start_position,sol_i);
}

emili::Neighborhood::NeighborhoodIterator emili::pfsp::KarNeighborhood::begin(Solution *base)
{
    ep_iterations = 1;
    sp_iterations = 1;
    std::vector< int > sol(((emili::pfsp::PermutationFlowShopSolution*)base)->getJobSchedule());
    //sol.erase(sol.begin()+start_position);
#ifdef ENABLE_SSE
    computeHEAD(sol,head,pmatrix,njobs,nmac);
#else
    njobs++;
    computeHead(sol);
    njobs--;
#endif
    return emili::Neighborhood::NeighborhoodIterator(this,base);
}

emili::Solution* emili::pfsp::KarNeighborhood::computeStep(emili::Solution* value)
{
    emili::iteration_increment();
    if(sp_iterations > njobs)
    {
        return nullptr;
    }
    else{
       sp_iterations++;
        std::vector < int >& newsol = ((emili::pfsp::PermutationFlowShopSolution*)value)->getJobSchedule();
        float p = emili::generateRealRandomNumber();
        start_position = (emili::generateRandomNumber()%njobs)+1;
        end_position = (emili::generateRandomNumber()%njobs)+1;
        if(p<0.5)
        {
            int sol_i = newsol[start_position];
            newsol.erase(newsol.begin()+start_position);
            newsol.insert(newsol.begin()+end_position,sol_i);
            lastMoveType=0;
        }
        else
        {
            std::swap(newsol[start_position],newsol[end_position]);
            lastMoveType=1;
        }
        //long int new_value = pis.computeObjectiveFunction(newsol);
        int starting = start_position>end_position?end_position:start_position;
        long int new_value = pis.computeObjectiveFunctionFromHead(newsol,starting,this->head,njobs);
        //std::cout << new_v << " " << new_value << std::endl;
      //  assert(new_v==new_value);
        value->setSolutionValue(new_value);
        return value;
    }
}

void emili::pfsp::KarNeighborhood::reverseLastMove(Solution *step)
{
    std::vector < int >& newsol = ((emili::pfsp::PermutationFlowShopSolution*)step)->getJobSchedule();
    if(lastMoveType)
    {
       std::swap(newsol[start_position],newsol[end_position]);
    }
    else
    {
        int sol_i = newsol[end_position];
        newsol.erase(newsol.begin()+end_position);
        newsol.insert(newsol.begin()+start_position,sol_i);
    }
}

void emili::pfsp::TaillardAcceleratedInsertNeighborhood::computeTAmatrices(std::vector<int> &sol)
{     
    int j,m;
    int jobNumber;
    int end_i = njobs-1;
   // std::vector< std::vector < int >> head(previousMachineEndTimeMatri);
    int prevj = 0;
    int postj = 0;
    int k;
    for(j=1;j<njobs;j++)
    {
        k = njobs-j;
        jobNumber = sol[j];
        prevj = prevj + pmatrix[jobNumber][1];
        postj = postj + pmatrix[sol[k]][nmac];
        head[1][j] = prevj;        
        tail[nmac][k] = postj;        
    }

      for ( j = 1; j <= end_i; ++j )
        {
            k = njobs-j;
            long int previousJobEndTime = head[1][j];
            long int postJobEndTime = tail[nmac][k];

            jobNumber = sol[j];

            for ( m = 2; m <= nmac; ++m )
            {
                int n = nmac-m+1;

            if ( head[m][j-1] > previousJobEndTime )
            {
                head[m][j] = head[m][j-1] + pmatrix[jobNumber][m];

            }
            else
            {
                head[m][j] = previousJobEndTime + pmatrix[jobNumber][m];
            }

            if ( tail[n][k+1] > postJobEndTime )
            {
                tail[n][k] = tail[n][k+1] + pmatrix[sol[k]][n];

            }
            else
            {
                tail[n][k] = postJobEndTime + pmatrix[sol[k]][n];
            }            
            previousJobEndTime = head[m][j];
            postJobEndTime = tail[n][k];
        }
    }
}

emili::Solution* emili::pfsp::TaillardAcceleratedInsertNeighborhood::computeStep(emili::Solution *value)
{
    emili::iteration_increment();
    if(sp_iterations > njobs)
    {        
        return nullptr;
    }
    else
    {
        end_position = ((end_position)%njobs)+1;
        std::vector < int >& newsol = ((emili::pfsp::PermutationFlowShopSolution*)value)->getJobSchedule();
        int sol_i;
        if(ep_iterations < njobs){
            ep_iterations++;
            /** if(ep_iterations == sp_iterations){
                ep_iterations++;
                end_position++;
            }*/
            if(end_position == start_position-1)
            {
                end_position=((end_position+1)%njobs)+1;
                ep_iterations+=2;
                if(ep_iterations > njobs && sp_iterations+1 > njobs)
                    return nullptr;
            }
            sol_i = newsol[start_position];
            newsol.erase(newsol.begin()+start_position);
        }
        else
        {
            sp_iterations++;
            ep_iterations = 1;
            start_position = ((start_position)%njobs)+1;
            if(end_position == start_position-1)
            {
                end_position=((end_position+1)%njobs)+1;
                ep_iterations+=2;
                if(ep_iterations > njobs && sp_iterations+1 > njobs)
                    return nullptr;
            }
            sol_i = newsol[start_position];
            newsol.erase(newsol.begin()+start_position);
#ifdef ENABLE_SSE
            computeHEADandTAIL(newsol,head,tail,pmatrix,njobs-1,nmac);
#else
            computeTAmatrices(newsol);

#endif
        }
        newsol.insert(newsol.begin()+end_position,sol_i);
        long int c_cur = head[1][end_position-1]+pmatrix[sol_i][1];
        long int c_max = c_cur+tail[1][end_position];
        for (int i = 2; i <= nmac; ++i) {
            int c_pm = head[i][end_position-1];

            if(c_pm > c_cur)
            {
                c_cur = c_pm;
            }

            c_cur = c_cur + pmatrix[sol_i][i];

            long int c_can = (c_cur+tail[i][end_position]);

            if(c_can>value->getSolutionValue())
            {
                value->setSolutionValue(c_can);
                return value;
            }

            c_max = c_max>c_can?c_max:c_can;
        }
        //long int old_v  = pis.computeObjectiveFunction(newsol);
        //std::cout << c_max << " - " << old_v << std::endl;
        //assert(c_max == old_v);
        value->setSolutionValue(c_max);
        return value;
    }
}

emili::Solution* emili::pfsp::CSTaillardAcceleratedInsertNeighborhood::computeStep(emili::Solution *value)
{
    emili::iteration_increment();
    if(sp_iterations > njobs)
    {
        return nullptr;
    }
    else
    {
        std::vector < int >& newsol = ((emili::pfsp::PermutationFlowShopSolution*)value)->getJobSchedule();

        int best_inspos = 1;
        int best_cmax = std::numeric_limits<int>::max();
        end_position = 1;
        start_position = ((start_position)%njobs)+1;
        int sol_i = newsol[start_position];
        newsol.erase(newsol.begin()+start_position);
#ifdef ENABLE_SSE
        computeHEADandTAIL(newsol,head,tail,pmatrix,njobs-1,nmac);
#else
        computeTAmatrices(newsol);

#endif

        for(int k=1; k<=njobs; k++)
        {

            if(k != start_position)
            {
                long int c_cur = head[1][k-1]+pmatrix[sol_i][1];
                long int c_max = c_cur+tail[1][k];
                for (int i = 2; i <= nmac; ++i) {
                    int c_pm = head[i][k-1];

                    if(c_pm > c_cur)
                    {
                        c_cur = c_pm;
                    }

                    c_cur = c_cur + pmatrix[sol_i][i];
                    long int c_can = (c_cur+tail[i][k]);
                    c_max = c_max>c_can?c_max:c_can;
                }

                if(c_max < best_cmax)
                {
                    best_cmax = c_max;
                    best_inspos = k;

                }

            }

        }
        //long int old_v  = pis.computeObjectiveFunction(newsol);
        //std::cout << c_max << " - " << old_v << std::endl;
        //assert(c_max == old_v);
        end_position = best_inspos;
        newsol.insert(newsol.begin()+best_inspos,sol_i);
        value->setSolutionValue(best_cmax);
        sp_iterations++;
        return value;
    }
}

emili::Solution* emili::pfsp::FSTaillardAcceleratedInsertNeighborhood::computeStep(emili::Solution *value)
{
    emili::iteration_increment();
    if(sp_iterations > njobs)
    {
        return nullptr;
    }
    else
    {
        end_position = ((end_position)%njobs)+1;
         std::vector < int >& newsol = ((emili::pfsp::PermutationFlowShopSolution*)value)->getJobSchedule();
       // int ktest = pis.computeMS(newsol);
       // assert(pis.computeMS(newsol)==current_value);   
        int sol_i;
        if(ep_iterations < njobs){
            ep_iterations++;
            /** if(ep_iterations == sp_iterations){
                ep_iterations++;
                end_position++;
            }*/
            if(end_position == start_position-1)
            {
                end_position=((end_position+1)%njobs)+1;
                ep_iterations+=2;
                if(ep_iterations > njobs && sp_iterations+1 > njobs)
                    return nullptr;
            }
            sol_i = newsol[start_position];
            newsol.erase(newsol.begin()+start_position);
        }
        else
        {
            sp_iterations++;
            ep_iterations = 1;
            start_position = ((start_position)%njobs)+1;            
            if(end_position == start_position-1)
            {
                end_position=((end_position+1)%njobs)+1;
                ep_iterations+=2;
                if(ep_iterations > njobs && sp_iterations+1 > njobs)
                    return nullptr;
            }
            sol_i = newsol[start_position];
            newsol.erase(newsol.begin()+start_position);
#ifdef ENABLE_SSE
            computeHEADandTAIL(newsol,head,tail,pmatrix,njobs-1,nmac);
#else
            computeTAmatrices(newsol);

#endif
        }

        newsol.insert(newsol.begin()+end_position,sol_i);
        long int c_cur = head[1][end_position-1]+pmatrix[sol_i][1];
        long int c_max = c_cur+tail[1][end_position];        
        for (int i = 2; i <= nmac; ++i) {
            int c_pm = head[i][end_position-1];

            if(c_pm > c_cur)
            {
                c_cur = c_pm;
            }

            c_cur = c_cur + pmatrix[sol_i][i];

            long int c_can = (c_cur+tail[i][end_position]);

            c_max = c_max>c_can?c_max:c_can;
        }
        //long int old_vi  = pis.computeObjectiveFunction(newsol);
        //std::cout << c_max << " - " << old_v << std::endl;
        //assert(c_max == old_vi);

        if(c_max < current_value)
        {
           current_value = c_max;
           improved = true;
        }
        value->setSolutionValue(c_max);
        return value;
    }
}

void emili::pfsp::FSTaillardAcceleratedInsertNeighborhood::reverseLastMove(Solution *step)
{
    if(!improved)
    {
        std::vector < int >& newsol = ((emili::pfsp::PermutationFlowShopSolution*)step)->getJobSchedule();
        int sol_i = newsol[end_position];
        newsol.erase(newsol.begin()+end_position);
        newsol.insert(newsol.begin()+start_position,sol_i);
    }
    else
    {
        ep_iterations = 1;
        std::vector< int > sol(((emili::pfsp::PermutationFlowShopSolution*)step)->getJobSchedule());        
        sol.erase(sol.begin()+start_position);
    #ifdef ENABLE_SSE
                computeHEADandTAIL(sol,head,tail,pmatrix,njobs-1,nmac);
    #else
                computeTAmatrices(sol);
    #endif
        improved = false;
    }
}

emili::Solution* emili::pfsp::OptInsert::computeStep(emili::Solution *value)
{
    emili::iteration_increment();
    if(sp_iterations > njobs)
    {
        return nullptr;
    }
    else
    {
        end_position = ((end_position)%njobs)+1;
        std::vector < int >& newsol = ((emili::pfsp::PermutationFlowShopSolution*)value)->getJobSchedule();
        int sol_i;
        if(ep_iterations < njobs){
            ep_iterations++;
            /** if(ep_iterations == sp_iterations){
                ep_iterations++;
                end_position++;
            }*/
            if(end_position == start_position-1)
            {
                end_position=((end_position+1)%njobs)+1;
                ep_iterations+=2;
                if(ep_iterations > njobs && sp_iterations+1 > njobs)
                    return nullptr;
            }
            sol_i = newsol[start_position];
            newsol.erase(newsol.begin()+start_position);
        }
        else
        {
            sp_iterations++;
            ep_iterations = 1;
            start_position = ((start_position)%njobs)+1;
            if(end_position == start_position-1)
            {
                end_position=((end_position+1)%njobs)+1;
                ep_iterations+=2;
                if(ep_iterations > njobs && sp_iterations+1 > njobs)
                    return nullptr;
            }
            sol_i = newsol[start_position];
            newsol.erase(newsol.begin()+start_position);
#ifdef ENABLE_SSE
            computeHEADandTAIL(newsol,head,tail,pmatrix,njobs-1,nmac);
#else
            computeTAmatrices(newsol);

#endif
        }
        newsol.insert(newsol.begin()+end_position,sol_i);
        int ins_pos[nmac+1];
        long int c_cur = head[1][end_position-1]+pmatrix[sol_i][1];
        ins_pos[1] = c_cur;
        long int c_max = c_cur+tail[1][end_position];
        for (int i = 2; i <= nmac; ++i) {
            int c_pm = head[i][end_position-1];
            if(c_pm > c_cur)
            {
                c_cur = c_pm;
            }
            c_cur += pmatrix[sol_i][i];
            ins_pos[i] =  c_cur;
            long int c_can = (c_cur+tail[i][end_position]);
            c_max = c_max>c_can?c_max:c_can;
        }

        int wt = (std::max(c_cur - pis.getDueDate(sol_i), 0L) * pis.getPriority(sol_i));

        if(end_position < njobs)
            wt += (std::max(c_max - pis.getDueDate(newsol[njobs]), 0L) * pis.getPriority(newsol[njobs]));

        long int pre_c_cur = c_cur;

        for (int j = 1; j< end_position; ++j )
        {
            wt += (std::max((long int)head[nmac][j] - pis.getDueDate(newsol[j]), 0L) * pis.getPriority(newsol[j]));
        }

        int pre_wt = wt;
        for(int k=end_position+1; k< njobs; k++)
        {
            int job = newsol[k];
            pre_c_cur = ins_pos[1] + pmatrix[job][1];
            ins_pos[1] = pre_c_cur;
            for(int m=2; m <= nmac ; m++)
            {
                int c_pm = ins_pos[m];
                if(c_pm > pre_c_cur)
                {
                    pre_c_cur = c_pm;
                }
                pre_c_cur += pmatrix[job][m];
                ins_pos[m] = pre_c_cur;
            }
            pre_wt += (std::max(pre_c_cur - pis.getDueDate(newsol[k]), 0L) * pis.getPriority(newsol[k]));
        }
        value->setSolutionValue(pre_wt);
        return value;
    }
}

void emili::pfsp::HeavilyApproximatedTaillardAcceleratedInsertNeighborhood::computeHead(std::vector<int>& sol)
{
    int j,m;
    int jobNumber;
    int prevj = 0;
    for(j=1;j<njobs;j++)
    {
        jobNumber = sol[j];
        prevj = prevj + pmatrix[jobNumber][1];
        head[1][j] = prevj;
    }

      for ( j = 1; j < njobs; ++j )
        {
            long int previousJobEndTime = head[1][j];
            jobNumber = sol[j];

            for ( m = 2; m <= nmac; ++m )
            {
                if ( head[m][j-1] > previousJobEndTime )
                {
                    head[m][j] = head[m][j-1] + pmatrix[jobNumber][m];
                }
                else
                {
                    head[m][j] = previousJobEndTime + pmatrix[jobNumber][m];
                }
                previousJobEndTime = head[m][j];
            }
    }
}

emili::Solution* emili::pfsp::HeavilyApproximatedTaillardAcceleratedInsertNeighborhood::computeStep(emili::Solution *value)
{
    emili::iteration_increment();
    if(sp_iterations > njobs)
    {
        return nullptr;
    }
    else
    {
        end_position = ((end_position)%njobs)+1;
        std::vector < int >& newsol = ((emili::pfsp::PermutationFlowShopSolution*)value)->getJobSchedule();
        int sol_i;
        if(ep_iterations < njobs){
            ep_iterations++;
            /** if(ep_iterations == sp_iterations){
                ep_iterations++;
                end_position++;
            }*/
            if(end_position == start_position-1)
            {
                end_position=((end_position+1)%njobs)+1;
                ep_iterations+=2;
                if(ep_iterations > njobs && sp_iterations+1 > njobs)
                    return nullptr;
            }
            sol_i = newsol[start_position];
            newsol.erase(newsol.begin()+start_position);
        }
        else
        {
            sp_iterations++;
            ep_iterations = 1;
            start_position = ((start_position)%njobs)+1;
            if(end_position == start_position-1)
            {
                end_position=((end_position+1)%njobs)+1;
                ep_iterations+=2;
                if(ep_iterations > njobs && sp_iterations+1 > njobs)
                    return nullptr;
            }
            sol_i = newsol[start_position];
            newsol.erase(newsol.begin()+start_position);
#ifdef ENABLE_SSE
            computeHEAD(newsol,head,pmatrix,njobs-1,nmac);
#else
            computeHead(newsol);
#endif
        }

        newsol.insert(newsol.begin()+end_position,sol_i);
        //std::vector< int > ins_pos(nmac+1,0);
        int ins_pos[nmac+1];
        long int c_cur = head[1][end_position-1]+pmatrix[sol_i][1];
        ins_pos[1] = c_cur;

        for (int i = 2; i <= nmac; ++i) {
           /**  int c_pm = head[i][end_position-1];
            if(c_pm < c_cur)
            {
                c_cur = c_cur + pmatrix[sol_i][i];
            }
            else
            {
                c_cur = c_pm + pmatrix[sol_i][i];
            }*/
            int c_pm = head[i][end_position-1];
            if(c_pm > c_cur)
            {
                c_cur = c_pm;
            }
            c_cur += pmatrix[sol_i][i];
           ins_pos[i] =  c_cur;

        }

        int wt = (std::max(c_cur - pis.getDueDate(sol_i), 0L) * pis.getPriority(sol_i));

        long int pre_c_cur = c_cur;


        for (int j = 1; j< end_position; ++j )
        {
            wt += (std::max((long int)head[nmac][j] - pis.getDueDate(newsol[j]), 0L) * pis.getPriority(newsol[j]));
        }

        int pre_wt = wt;

        for(int k=end_position+1; k<= njobs; k++)
        {
            pre_c_cur = pre_c_cur + pmatrix[newsol[k]][nmac];
            wt += (std::max(pre_c_cur - pis.getDueDate(newsol[k]), 0L) * pis.getPriority(newsol[k]));
        }


        int value_wt = value->getSolutionValue();

        if(wt < value_wt)
        {
            for(int k=end_position+1; k<= njobs; k++)
            {
                int job = newsol[k];
                pre_c_cur = ins_pos[1] + pmatrix[job][1];
                ins_pos[1] = pre_c_cur;
                for(int m=2; m <= nmac ; m++)
                {
                    int c_pm = ins_pos[m];
                    if(c_pm > pre_c_cur)
                    {
                        pre_c_cur = c_pm;
                    }
                    pre_c_cur += pmatrix[job][m];
                    ins_pos[m] = pre_c_cur;
                }
               pre_wt += (std::max(pre_c_cur - pis.getDueDate(newsol[k]), 0L) * pis.getPriority(newsol[k]));
              if(pre_wt > value_wt)
               {
                  // cout << "exited at " << k << std::endl;

                   value->setSolutionValue(pre_wt);
                   return value;
                }

            }
             //   std::cout << pre_wt << " "<< pis.computeObjectiveFunction(newsol);
            //assert( pre_wt == pis.computeObjectiveFunction(newsol));

            value->setSolutionValue(pre_wt);
           //pis.evaluateSolution(*news);
        /**    if(news->getSolutionValue() > value->getSolutionValue())
           {
               emili::iteration_increment();
           }*/
        }
        else
        {
            value->setSolutionValue(wt);
        }

        //long int old_v  = pis.computeObjectiveFunction(newsol);
        //std::cout << c_max << " - " << old_v << std::endl;
        //assert(c_max == old_v);
        return value;
    }
}

emili::Solution* emili::pfsp::NatxNeighborhood::computeStep(emili::Solution *value)
{
    emili::iteration_increment();
    if(sp_iterations > njobs)
    {
        return nullptr;
    }
    else
    {
        end_position = ((end_position)%njobs)+1;
        std::vector < int >& newsol = ((emili::pfsp::PermutationFlowShopSolution*)value)->getJobSchedule();
        int sol_i;
        if(ep_iterations < njobs){
            ep_iterations++;
            /** if(ep_iterations == sp_iterations){
                ep_iterations++;
                end_position++;
            }*/
            if(end_position == start_position-1)
            {
                end_position=((end_position+1)%njobs)+1;
                ep_iterations+=2;
                if(ep_iterations > njobs && sp_iterations+1 > njobs)
                    return nullptr;
            }
            sol_i = newsol[start_position];
            newsol.erase(newsol.begin()+start_position);
        }
        else
        {
            sp_iterations++;
            ep_iterations = 1;
            start_position = ((start_position)%njobs)+1;
            if(end_position == start_position-1)
            {
                end_position=((end_position+1)%njobs)+1;
                ep_iterations+=2;
                if(ep_iterations > njobs && sp_iterations+1 > njobs)
                    return nullptr;
            }
            sol_i = newsol[start_position];
            newsol.erase(newsol.begin()+start_position);
#ifdef ENABLE_SSE
            computeHEAD(newsol,head,pmatrix,njobs-1,nmac);
#else
            computeHead(newsol);
#endif
        }


        newsol.insert(newsol.begin()+end_position,sol_i);
        //std::vector< int > ins_pos(nmac+1,0);
        int ins_pos[nmac+1];
        long int c_cur = head[1][end_position-1]+pmatrix[sol_i][1];
        ins_pos[1] = c_cur;

        for (int i = 2; i <= nmac; ++i) {
            int c_pm = head[i][end_position-1];
            if(c_pm > c_cur)
            {
                c_cur = c_pm;
            }
            c_cur += pmatrix[sol_i][i];
           ins_pos[i] =  c_cur;

        }
        long int pre_c_cur = c_cur;
        int wt = (std::max(c_cur - pis.getDueDate(sol_i), 0L) * pis.getPriority(sol_i));

        for (int j = 1; j< end_position; ++j )
        {
            wt += (std::max((long int)head[nmac][j] - pis.getDueDate(newsol[j]), 0L) * pis.getPriority(newsol[j]));
        }

        int pre_wt = wt;
        if(end_position > njobs/2)
        for(int k=end_position+1; k<= njobs; k++)
        {
            pre_c_cur = pre_c_cur + pmatrix[newsol[k]][nmac];
            wt += (std::max(pre_c_cur - pis.getDueDate(newsol[k]), 0L) * pis.getPriority(newsol[k]));
        }


        int value_wt = value->getSolutionValue();

        if(wt < value_wt)
        {
            //emili::iteration_decrement();
            for(int k=end_position+1; k<= njobs; k++)
            {
                int job = newsol[k];
                pre_c_cur = ins_pos[1] + pmatrix[job][1];
                ins_pos[1] = pre_c_cur;
                for(int m=2; m <= nmac ; m++)
                {
                    int c_pm = ins_pos[m];
                    if(c_pm > pre_c_cur)
                    {
                        pre_c_cur = c_pm;
                    }
                    pre_c_cur += pmatrix[job][m];
                    ins_pos[m] = pre_c_cur;
                }
               pre_wt += (std::max(pre_c_cur - pis.getDueDate(newsol[k]), 0L) * pis.getPriority(newsol[k]));
              if(pre_wt > value_wt)
                {
                   value->setSolutionValue(pre_wt);
                   return value;
                }

            }

            value->setSolutionValue(pre_wt);

        }
        else
        {
            value->setSolutionValue(wt);
        }

        return value;
    }
}

emili::Solution* emili::pfsp::AtxNeighborhood::computeStep(emili::Solution *value)
{
    emili::iteration_increment();
    if(sp_iterations > njobs)
    {
        return nullptr;
    }
    else
    {
        end_position = ((end_position)%njobs)+1;
        std::vector < int >& newsol = ((emili::pfsp::PermutationFlowShopSolution*)value)->getJobSchedule();
        int sol_i;
        if(ep_iterations < njobs){
            ep_iterations++;
            /** if(ep_iterations == sp_iterations){
                ep_iterations++;
                end_position++;
            }*/
            if(end_position == start_position-1)
            {
                end_position=((end_position+1)%njobs)+1;
                ep_iterations+=2;
                if(ep_iterations > njobs && sp_iterations+1 > njobs)
                    return nullptr;
            }
            sol_i = newsol[start_position];
            newsol.erase(newsol.begin()+start_position);
        }
        else
        {
            sp_iterations++;
            ep_iterations = 1;
            start_position = ((start_position)%njobs)+1;
            if(end_position == start_position-1)
            {
                end_position=((end_position+1)%njobs)+1;
                ep_iterations+=2;
                if(ep_iterations > njobs && sp_iterations+1 > njobs)
                    return nullptr;
            }
            sol_i = newsol[start_position];
            newsol.erase(newsol.begin()+start_position);
#ifdef ENABLE_SSE
            computeHEAD(newsol,head,pmatrix,njobs-1,nmac);
#else
            computeHead(newsol);
#endif
        }


        newsol.insert(newsol.begin()+end_position,sol_i);
        //std::vector< int > ins_pos(nmac+1,0);
        int ins_pos[nmac+1];
        long int c_cur = head[1][end_position-1]+pmatrix[sol_i][1];
        ins_pos[1] = c_cur;

        for (int i = 2; i <= nmac; ++i) {
            int c_pm = head[i][end_position-1];
            if(c_pm > c_cur)
            {
                c_cur = c_pm;
            }
            c_cur += pmatrix[sol_i][i];
           ins_pos[i] =  c_cur;

        }
        long int pre_c_cur = c_cur;
        int wt = (std::max(c_cur - pis.getDueDate(sol_i), 0L) * pis.getPriority(sol_i));

        for (int j = 1; j< end_position; ++j )
        {
            wt += (std::max((long int)head[nmac][j] - pis.getDueDate(newsol[j]), 0L) * pis.getPriority(newsol[j]));
        }

        for(int k=end_position+1; k<= njobs; k++)
            {
                int job = newsol[k];
                pre_c_cur = ins_pos[1] + pmatrix[job][1];
                ins_pos[1] = pre_c_cur;
                for(int m=2; m <= nmac ; m++)
                {
                    int c_pm = ins_pos[m];
                    if(c_pm > pre_c_cur)
                    {
                        pre_c_cur = c_pm;
                    }
                    pre_c_cur += pmatrix[job][m];
                    ins_pos[m] = pre_c_cur;
                }
               wt += (std::max(pre_c_cur - pis.getDueDate(newsol[k]), 0L) * pis.getPriority(newsol[k]));
            }

            value->setSolutionValue(wt);

        }
        return value;

}

emili::Solution* emili::pfsp::Natx2Neighborhood::computeStep(emili::Solution *value)
{

    emili::iteration_increment();
    if(sp_iterations > njobs)
    {
        return nullptr;
    }
    else
    {
        end_position = ((end_position)%njobs)+1;
        std::vector < int >& newsol = ((emili::pfsp::PermutationFlowShopSolution*)value)->getJobSchedule();
        int sol_i;
        if(ep_iterations < njobs){
            ep_iterations++;
            /** if(ep_iterations == sp_iterations){
                ep_iterations++;
                end_position++;
            }*/
            if(end_position == start_position-1)
            {
                end_position=((end_position+1)%njobs)+1;
                ep_iterations+=2;
                if(ep_iterations > njobs && sp_iterations+1 > njobs)
                    return nullptr;
            }
            sol_i = newsol[start_position];
            newsol.erase(newsol.begin()+start_position);
        }
        else
        {
            sp_iterations++;
            ep_iterations = 1;
            start_position = ((start_position)%njobs)+1;
            if(end_position == start_position-1)
            {
                end_position=((end_position+1)%njobs)+1;
                ep_iterations+=2;
                if(ep_iterations > njobs && sp_iterations+1 > njobs)
                    return nullptr;
            }
            sol_i = newsol[start_position];
            newsol.erase(newsol.begin()+start_position);
#ifdef ENABLE_SSE
            computeHEAD(newsol,head,pmatrix,njobs-1,nmac);
#else
            computeHead(newsol);
#endif
        }

        newsol.insert(newsol.begin()+end_position,sol_i);
/** #ifdef ENABLE_SSE
        float ins_pos[nmac+1];
#else*/
        int ins_pos[nmac+1];
//#endif
        long int c_cur = head[1][end_position-1]+pmatrix[sol_i][1];
        ins_pos[1] = c_cur;

        for (int i = 2; i <= nmac; ++i) {
            int c_pm = head[i][end_position-1];
            if(c_pm > c_cur)
            {
                c_cur = c_pm;
            }
            c_cur += pmatrix[sol_i][i];
           ins_pos[i] =  c_cur;
        }
        long int pre_c_cur = c_cur;
        int wt = (std::max(c_cur - duedates[sol_i], 0L) * priorities[sol_i]);

        for (int j = 1; j< end_position; ++j )
        {
            wt += (std::max((long int)head[nmac][j] - duedates[newsol[j]], 0L) * priorities[newsol[j]]);
        }

        int pre_wt = wt;
        if(end_position > thresh)
        for(int k=end_position+1; k<= njobs; k++)
        {
            int nk = newsol[k];
            pre_c_cur = pre_c_cur + pmatrix[nk][nmac];
            wt += (std::max(pre_c_cur - duedates[nk], 0L) * priorities[nk]);
        }

        if(wt < value_wt)
        {
            if(end_position > thresh)
            {
                thresh++;
            }
            /**
#ifdef ENABLE_SSE
            std::vector< long int > pmet(njobs+1,0);
            computePMakespans(newsol,pmet,pis.getProcessingTimesMatrix(),njobs+1,nmac,end_position+1,ins_pos);
           for(int k=end_position+1; k<= njobs; k++)
            {
               int job = newsol[k];
               pre_wt += (std::max(pmet[k] - pis.getDueDate(job), 0L) * pis.getPriority(job));
               if(pre_wt > value_wt)
                 {
                   break;
                 }
           }

#else
*/
            for(int k=end_position+1; k<= njobs; k++)
            {
                int job = newsol[k];
                pre_c_cur = ins_pos[1] + pmatrix[job][1];
                ins_pos[1] = pre_c_cur;
                //std::cout << "start critico " << std::endl;
                for(int m=2; m <= nmac ; m++)
                {                  
                    int c_pm = ins_pos[m];
                    if(c_pm > pre_c_cur)
                    {
                        pre_c_cur = c_pm;
                    }
                    pre_c_cur += pmatrix[job][m];
                    ins_pos[m] = pre_c_cur;
                }
            //    std::cout << "end critico " << std::endl;

               pre_wt += (std::max(pre_c_cur - duedates[job], 0L) * priorities[job]);
              if(pre_wt > value_wt)
                {
                  // value->setSolutionValue(pre_wt);
                 //  double t = ((double)clock()-s)/CLOCKS_PER_SEC;
                 // std::cout << "natx2 : " << t << std::endl;
                  // return value;
                  break;
                }

            }
//#endif

            value->setSolutionValue(pre_wt);

        }
        else
        {

           thresh--;
            value->setSolutionValue(wt);
        }

        return value;
    }
}

emili::Solution* emili::pfsp::EatxNeighborhood::computeStep(emili::Solution *value)
{
    emili::iteration_increment();
    if(sp_iterations > njobs)
    {
        return nullptr;
    }
    else
    {
        end_position = ((end_position)%njobs)+1;
        std::vector < int >& newsol = ((emili::pfsp::PermutationFlowShopSolution*)value)->getJobSchedule();
        int sol_i;
        if(ep_iterations < njobs){
            ep_iterations++;
            /** if(ep_iterations == sp_iterations){
                ep_iterations++;
                end_position++;
            }*/
            if(end_position == start_position-1)
            {
                end_position=((end_position+1)%njobs)+1;
                ep_iterations+=2;
                if(ep_iterations > njobs && sp_iterations+1 > njobs)
                    return nullptr;
            }
            sol_i = newsol[start_position];
            newsol.erase(newsol.begin()+start_position);
        }
        else
        {
            sp_iterations++;
            ep_iterations = 1;
            start_position = ((start_position)%njobs)+1;
            if(end_position == start_position-1)
            {
                end_position=((end_position+1)%njobs)+1;
                ep_iterations+=2;
                if(ep_iterations > njobs && sp_iterations+1 > njobs)
                    return nullptr;
            }
            sol_i = newsol[start_position];
            newsol.erase(newsol.begin()+start_position);
#ifdef ENABLE_SSE
            computeHEAD(newsol,head,pmatrix,njobs-1,nmac);
#else
            computeHead(newsol);
#endif
        }


        newsol.insert(newsol.begin()+end_position,sol_i);
        //std::vector< int > ins_pos(nmac+1,0);
        int ins_pos[nmac+1];
        long int c_cur = head[1][end_position-1]+pmatrix[sol_i][1];
        ins_pos[1] = c_cur;

        for (int i = 2; i <= nmac; ++i) {
            int c_pm = head[i][end_position-1];
            if(c_pm > c_cur)
            {
                c_cur = c_pm;
            }
            c_cur += pmatrix[sol_i][i];
           ins_pos[i] =  c_cur;

        }
        long int ppre_c_cur = ins_pos[nmac-1];
        long int pre_c_cur = c_cur;

        int wt = (std::max(c_cur - pis.getDueDate(sol_i), 0L) * pis.getPriority(sol_i));
        for (int j = 1; j< end_position; ++j )
        {
            wt += (std::max((long int)head[nmac][j] - pis.getDueDate(newsol[j]), 0L) * pis.getPriority(newsol[j]));
        }

        int pre_wt = wt;

        if(end_position >njobs/2)
        for(int k=end_position+1; k<= njobs; k++)
        {
            int job = newsol[k];
            ppre_c_cur = ppre_c_cur+pmatrix[job][nmac-1];
            pre_c_cur = std::max(pre_c_cur,ppre_c_cur) + pmatrix[job][nmac];
            wt += (std::max(pre_c_cur - pis.getDueDate(newsol[k]), 0L) * pis.getPriority(newsol[k]));
        }


        int value_wt = value->getSolutionValue();

        if(wt < value_wt)
        {
            for(int k=end_position+1; k<= njobs; k++)
            {
                int job = newsol[k];
                pre_c_cur = ins_pos[1] + pmatrix[job][1];
                ins_pos[1] = pre_c_cur;
                for(int m=2; m <= nmac ; m++)
                {
                    int c_pm = ins_pos[m];
                    if(c_pm > pre_c_cur)
                    {
                        pre_c_cur = c_pm;
                    }
                    pre_c_cur += pmatrix[job][m];
                    ins_pos[m] = pre_c_cur;
                }
               pre_wt += (std::max(pre_c_cur - pis.getDueDate(newsol[k]), 0L) * pis.getPriority(newsol[k]));
              if(pre_wt > value_wt)
                {
                   value->setSolutionValue(pre_wt);
                   return value;
                }

            }

            value->setSolutionValue(pre_wt);

        }
        else
        {
            value->setSolutionValue(wt);
        }
        return value;
    }
}

emili::Solution* emili::pfsp::ThatxNeighborhood::computeStep(emili::Solution *value)
{
    emili::iteration_increment();
    if(sp_iterations > njobs)
    {
        return nullptr;
    }
    else
    {
        end_position = ((end_position)%njobs)+1;
        std::vector < int >& newsol = ((emili::pfsp::PermutationFlowShopSolution*)value)->getJobSchedule();
        int sol_i;
        if(ep_iterations < njobs){
            ep_iterations++;
            /** if(ep_iterations == sp_iterations){
                ep_iterations++;
                end_position++;
            }*/
            if(end_position == start_position-1)
            {
                end_position=((end_position+1)%njobs)+1;
                ep_iterations+=2;
                if(ep_iterations > njobs && sp_iterations+1 > njobs)
                    return nullptr;
            }
            sol_i = newsol[start_position];
            newsol.erase(newsol.begin()+start_position);
        }
        else
        {
            sp_iterations++;
            ep_iterations = 1;
            start_position = ((start_position)%njobs)+1;
            if(end_position == start_position-1)
            {
                end_position=((end_position+1)%njobs)+1;
                ep_iterations+=2;
                if(ep_iterations > njobs && sp_iterations+1 > njobs)
                    return nullptr;
            }
            sol_i = newsol[start_position];
            newsol.erase(newsol.begin()+start_position);
#ifdef ENABLE_SSE
            computeHEAD(newsol,head,pmatrix,njobs-1,nmac);
#else
            computeHead(newsol);
#endif
        }


        newsol.insert(newsol.begin()+end_position,sol_i);
        //std::vector< int > ins_pos(nmac+1,0);
        int ins_pos[nmac+1];
        long int c_cur = head[1][end_position-1]+pmatrix[sol_i][1];
        ins_pos[1] = c_cur;

        for (int i = 2; i <= nmac; ++i) {
            int c_pm = head[i][end_position-1];
            if(c_pm > c_cur)
            {
                c_cur = c_pm;
            }
            c_cur += pmatrix[sol_i][i];
           ins_pos[i] =  c_cur;

        }
        long int pppre_c_cur = ins_pos[nmac-2];
        long int ppre_c_cur = ins_pos[nmac-1];
        long int pre_c_cur = c_cur;

        int wt = (std::max(c_cur - pis.getDueDate(sol_i), 0L) * pis.getPriority(sol_i));
        for (int j = 1; j< end_position; ++j )
        {
            wt += (std::max((long int)head[nmac][j] - pis.getDueDate(newsol[j]), 0L) * pis.getPriority(newsol[j]));
        }

        int pre_wt = wt;

        if(end_position >njobs/2)
        for(int k=end_position+1; k<= njobs; k++)
        {
            int job = newsol[k];
            pppre_c_cur = pppre_c_cur + pmatrix[job][nmac-2];
            ppre_c_cur = std::max(ppre_c_cur,pppre_c_cur)+pmatrix[job][nmac-1];
            pre_c_cur = std::max(pre_c_cur,ppre_c_cur) + pmatrix[job][nmac];
            wt += (std::max(pre_c_cur - pis.getDueDate(newsol[k]), 0L) * pis.getPriority(newsol[k]));
        }


        int value_wt = value->getSolutionValue();

        if(wt < value_wt)
        {
            for(int k=end_position+1; k<= njobs; k++)
            {
                int job = newsol[k];
                pre_c_cur = ins_pos[1] + pmatrix[job][1];
                ins_pos[1] = pre_c_cur;
                for(int m=2; m <= nmac ; m++)
                {
                    int c_pm = ins_pos[m];
                    if(c_pm > pre_c_cur)
                    {
                        pre_c_cur = c_pm;
                    }
                    pre_c_cur += pmatrix[job][m];
                    ins_pos[m] = pre_c_cur;
                }
               pre_wt += (std::max(pre_c_cur - pis.getDueDate(newsol[k]), 0L) * pis.getPriority(newsol[k]));
              if(pre_wt > value_wt)
                {
                   value->setSolutionValue(pre_wt);
                   return value;
                }

            }

            value->setSolutionValue(pre_wt);

        }
        else
        {
            value->setSolutionValue(wt);
        }
        return value;
    }
}

emili::Solution* emili::pfsp::FatxNeighborhood::computeStep(emili::Solution *value)
{
    emili::iteration_increment();
    if(sp_iterations > njobs)
    {
        return nullptr;
    }
    else
    {
        end_position = ((end_position)%njobs)+1;
        std::vector < int >& newsol = ((emili::pfsp::PermutationFlowShopSolution*)value)->getJobSchedule();
        int sol_i;
        if(ep_iterations < njobs){
            ep_iterations++;
            /** if(ep_iterations == sp_iterations){
                ep_iterations++;
                end_position++;
            }*/
            if(end_position == start_position-1)
            {
                end_position=((end_position+1)%njobs)+1;
                ep_iterations+=2;
                if(ep_iterations > njobs && sp_iterations+1 > njobs)
                    return nullptr;
            }
            sol_i = newsol[start_position];
            newsol.erase(newsol.begin()+start_position);
        }
        else
        {
            sp_iterations++;
            ep_iterations = 1;
            start_position = ((start_position)%njobs)+1;
            if(end_position == start_position-1)
            {
                end_position=((end_position+1)%njobs)+1;
                ep_iterations+=2;
                if(ep_iterations > njobs && sp_iterations+1 > njobs)
                    return nullptr;
            }
            sol_i = newsol[start_position];
            newsol.erase(newsol.begin()+start_position);
#ifdef ENABLE_SSE
            computeHEAD(newsol,head,pmatrix,njobs-1,nmac);
#else
            computeHead(newsol);
#endif
        }


        newsol.insert(newsol.begin()+end_position,sol_i);
        //std::vector< int > ins_pos(nmac+1,0);
        int ins_pos[nmac+1];
        long int c_cur = head[1][end_position-1]+pmatrix[sol_i][1];
        ins_pos[1] = c_cur;

        for (int i = 2; i <= nmac; ++i) {
            int c_pm = head[i][end_position-1];
            if(c_pm > c_cur)
            {
                c_cur = c_pm;
            }
            c_cur += pmatrix[sol_i][i];
           ins_pos[i] =  c_cur;

        }
        long int ppppre_c_cur = ins_pos[nmac-3];
        long int pppre_c_cur = ins_pos[nmac-2];
        long int ppre_c_cur = ins_pos[nmac-1];
        long int pre_c_cur = c_cur;

        int wt = (std::max(c_cur - pis.getDueDate(sol_i), 0L) * pis.getPriority(sol_i));
        for (int j = 1; j< end_position; ++j )
        {
            wt += (std::max((long int)head[nmac][j] - pis.getDueDate(newsol[j]), 0L) * pis.getPriority(newsol[j]));
        }

        int pre_wt = wt;

        if(end_position >njobs/2)
        for(int k=end_position+1; k<= njobs; k++)
        {
            int job = newsol[k];
            ppppre_c_cur = ppppre_c_cur + pmatrix[job][nmac-3];
            pppre_c_cur = std::max(pppre_c_cur,ppppre_c_cur)+pmatrix[job][nmac-2];
            ppre_c_cur = std::max(ppre_c_cur,pppre_c_cur)+pmatrix[job][nmac-1];
            pre_c_cur = std::max(pre_c_cur,ppre_c_cur) + pmatrix[job][nmac];
            wt += (std::max(pre_c_cur - pis.getDueDate(newsol[k]), 0L) * pis.getPriority(newsol[k]));
        }


        int value_wt = value->getSolutionValue();

        if(wt < value_wt)
        {
            for(int k=end_position+1; k<= njobs; k++)
            {
                int job = newsol[k];
                pre_c_cur = ins_pos[1] + pmatrix[job][1];
                ins_pos[1] = pre_c_cur;
                for(int m=2; m <= nmac ; m++)
                {
                    int c_pm = ins_pos[m];
                    if(c_pm > pre_c_cur)
                    {
                        pre_c_cur = c_pm;
                    }
                    pre_c_cur += pmatrix[job][m];
                    ins_pos[m] = pre_c_cur;
                }
               pre_wt += (std::max(pre_c_cur - pis.getDueDate(newsol[k]), 0L) * pis.getPriority(newsol[k]));
              if(pre_wt > value_wt)
                {
                   value->setSolutionValue(pre_wt);
                   return value;
                }

            }

            value->setSolutionValue(pre_wt);

        }
        else
        {
            value->setSolutionValue(wt);
        }
        return value;
    }
}

emili::Solution* emili::pfsp::PatxNeighborhood::computeStep(emili::Solution *value)
{
    emili::iteration_increment();
    if(sp_iterations > njobs)
    {
        return nullptr;
    }
    else
    {
        end_position = ((end_position)%njobs)+1;
        std::vector < int >& newsol = ((emili::pfsp::PermutationFlowShopSolution*)value)->getJobSchedule();
        int sol_i;
        if(ep_iterations < njobs){
            ep_iterations++;
            /** if(ep_iterations == sp_iterations){
                ep_iterations++;
                end_position++;
            }*/
            if(end_position == start_position-1)
            {
                end_position=((end_position+1)%njobs)+1;
                ep_iterations+=2;
                if(ep_iterations > njobs && sp_iterations+1 > njobs)
                    return nullptr;
            }
            sol_i = newsol[start_position];
            newsol.erase(newsol.begin()+start_position);
        }
        else
        {
            sp_iterations++;
            ep_iterations = 1;
            start_position = ((start_position)%njobs)+1;
            if(end_position == start_position-1)
            {
                end_position=((end_position+1)%njobs)+1;
                ep_iterations+=2;
                if(ep_iterations > njobs && sp_iterations+1 > njobs)
                    return nullptr;
            }
            sol_i = newsol[start_position];
            newsol.erase(newsol.begin()+start_position);
#ifdef ENABLE_SSE
            computeHEAD(newsol,head,pmatrix,njobs-1,nmac);
#else
            computeHead(newsol);
#endif
        }

        newsol.insert(newsol.begin()+end_position,sol_i);
        //std::vector< int > ins_pos(nmac+1,0);
        int ins_pos[nmac+1];
        long int c_cur = head[1][end_position-1]+pmatrix[sol_i][1];
        ins_pos[1] = c_cur;

        for (int i = 2; i <= nmac; ++i) {
            int c_pm = head[i][end_position-1];
            if(c_pm > c_cur)
            {
                c_cur = c_pm;
            }
            c_cur += pmatrix[sol_i][i];
           ins_pos[i] =  c_cur;

        }
        long int pppppre_c_cur = ins_pos[nmac-4];
        long int ppppre_c_cur = ins_pos[nmac-3];
        long int pppre_c_cur = ins_pos[nmac-2];
        long int ppre_c_cur = ins_pos[nmac-1];
        long int pre_c_cur = c_cur;

        int wt = (std::max(c_cur - pis.getDueDate(sol_i), 0L) * pis.getPriority(sol_i));
        for (int j = 1; j< end_position; ++j )
        {
            wt += (std::max((long int)head[nmac][j] - pis.getDueDate(newsol[j]), 0L) * pis.getPriority(newsol[j]));
        }

        int pre_wt = wt;

        if(end_position >njobs/2)
        for(int k=end_position+1; k<= njobs; k++)
        {
            int job = newsol[k];
            pppppre_c_cur = pppppre_c_cur + pmatrix[job][nmac-4];
            ppppre_c_cur = std::max(ppppre_c_cur,pppppre_c_cur)+pmatrix[job][nmac-3];
            pppre_c_cur = std::max(pppre_c_cur,ppppre_c_cur)+pmatrix[job][nmac-2];
            ppre_c_cur = std::max(ppre_c_cur,pppre_c_cur)+pmatrix[job][nmac-1];
            pre_c_cur = std::max(pre_c_cur,ppre_c_cur) + pmatrix[job][nmac];
            wt += (std::max(pre_c_cur - pis.getDueDate(newsol[k]), 0L) * pis.getPriority(newsol[k]));
        }


        int value_wt = value->getSolutionValue();

        if(wt < value_wt)
        {
            for(int k=end_position+1; k<= njobs; k++)
            {
                int job = newsol[k];
                pre_c_cur = ins_pos[1] + pmatrix[job][1];
                ins_pos[1] = pre_c_cur;
                for(int m=2; m <= nmac ; m++)
                {
                    int c_pm = ins_pos[m];
                    if(c_pm > pre_c_cur)
                    {
                        pre_c_cur = c_pm;
                    }
                    pre_c_cur += pmatrix[job][m];
                    ins_pos[m] = pre_c_cur;
                }
               pre_wt += (std::max(pre_c_cur - pis.getDueDate(newsol[k]), 0L) * pis.getPriority(newsol[k]));
              if(pre_wt > value_wt)
                {
                   value->setSolutionValue(pre_wt);
                   return value;
                }

            }

            value->setSolutionValue(pre_wt);

        }
        else
        {
            value->setSolutionValue(wt);
        }
        return value;
    }
}

emili::Solution* emili::pfsp::SatxNeighborhood::computeStep(emili::Solution *value)
{
    emili::iteration_increment();
    if(sp_iterations > njobs)
    {
        return nullptr;
    }
    else
    {
        end_position = ((end_position)%njobs)+1;
        std::vector < int >& newsol = ((emili::pfsp::PermutationFlowShopSolution*)value)->getJobSchedule();
        int sol_i;
        if(ep_iterations < njobs){
            ep_iterations++;
            /** if(ep_iterations == sp_iterations){
                ep_iterations++;
                end_position++;
            }*/
            if(end_position == start_position-1)
            {
                end_position=((end_position+1)%njobs)+1;
                ep_iterations+=2;
                if(ep_iterations > njobs && sp_iterations+1 > njobs)
                    return nullptr;
            }
            sol_i = newsol[start_position];
            newsol.erase(newsol.begin()+start_position);
        }
        else
        {
            sp_iterations++;
            ep_iterations = 1;
            start_position = ((start_position)%njobs)+1;
            if(end_position == start_position-1)
            {
                end_position=((end_position+1)%njobs)+1;
                ep_iterations+=2;
                if(ep_iterations > njobs && sp_iterations+1 > njobs)
                    return nullptr;
            }
            sol_i = newsol[start_position];
            newsol.erase(newsol.begin()+start_position);
#ifdef ENABLE_SSE
            computeHEAD(newsol,head,pmatrix,njobs-1,nmac);
#else
            computeHead(newsol);
#endif
        }


        newsol.insert(newsol.begin()+end_position,sol_i);
        //std::vector< int > ins_pos(nmac+1,0);
        int ins_pos[nmac+1];
        long int c_cur = head[1][end_position-1]+pmatrix[sol_i][1];
        ins_pos[1] = c_cur;

        for (int i = 2; i <= nmac; ++i) {
            int c_pm = head[i][end_position-1];
            if(c_pm > c_cur)
            {
                c_cur = c_pm;
            }
            c_cur += pmatrix[sol_i][i];
           ins_pos[i] =  c_cur;

        }
        long int ppppppre_c_cur = ins_pos[nmac-5];
        long int pppppre_c_cur = ins_pos[nmac-4];
        long int ppppre_c_cur = ins_pos[nmac-3];
        long int pppre_c_cur = ins_pos[nmac-2];
        long int ppre_c_cur = ins_pos[nmac-1];
        long int pre_c_cur = c_cur;

        int wt = (std::max(c_cur - pis.getDueDate(sol_i), 0L) * pis.getPriority(sol_i));
        for (int j = 1; j< end_position; ++j )
        {
            wt += (std::max((long int)head[nmac][j] - pis.getDueDate(newsol[j]), 0L) * pis.getPriority(newsol[j]));
        }

        int pre_wt = wt;

        if(end_position >njobs/2)
        for(int k=end_position+1; k<= njobs; k++)
        {
            int job = newsol[k];
            ppppppre_c_cur = ppppppre_c_cur + pmatrix[job][nmac-5];
            pppppre_c_cur = std::max(pppppre_c_cur,ppppppre_c_cur)+pmatrix[job][nmac-4];
            ppppre_c_cur = std::max(ppppre_c_cur,pppppre_c_cur)+pmatrix[job][nmac-3];
            pppre_c_cur = std::max(pppre_c_cur,ppppre_c_cur)+pmatrix[job][nmac-2];
            ppre_c_cur = std::max(ppre_c_cur,pppre_c_cur)+pmatrix[job][nmac-1];
            pre_c_cur = std::max(pre_c_cur,ppre_c_cur) + pmatrix[job][nmac];
            wt += (std::max(pre_c_cur - pis.getDueDate(newsol[k]), 0L) * pis.getPriority(newsol[k]));
        }


        int value_wt = value->getSolutionValue();

        if(wt < value_wt)
        {
            for(int k=end_position+1; k<= njobs; k++)
            {
                int job = newsol[k];
                pre_c_cur = ins_pos[1] + pmatrix[job][1];
                ins_pos[1] = pre_c_cur;
                for(int m=2; m <= nmac ; m++)
                {
                    int c_pm = ins_pos[m];
                    if(c_pm > pre_c_cur)
                    {
                        pre_c_cur = c_pm;
                    }
                    pre_c_cur += pmatrix[job][m];
                    ins_pos[m] = pre_c_cur;
                }
               pre_wt += (std::max(pre_c_cur - pis.getDueDate(newsol[k]), 0L) * pis.getPriority(newsol[k]));
              if(pre_wt > value_wt)
                {
                   value->setSolutionValue(pre_wt);
                   return value;
                }

            }

            value->setSolutionValue(pre_wt);

        }
        else
        {
            value->setSolutionValue(wt);
        }
        return value;
    }
}


emili::Solution* emili::pfsp::TatxNeighborhood::computeStep(emili::Solution *value)
{
    emili::iteration_increment();
    if(sp_iterations > njobs)
    {
        return nullptr;
    }
    else
    {
        end_position = ((end_position)%njobs)+1;
        std::vector < int >& newsol = ((emili::pfsp::PermutationFlowShopSolution*)value)->getJobSchedule();
        int sol_i;
        if(ep_iterations < njobs){
            ep_iterations++;
            /** if(ep_iterations == sp_iterations){
                ep_iterations++;
                end_position++;
            }*/
            if(end_position == start_position-1)
            {
                end_position=((end_position+1)%njobs)+1;
                ep_iterations+=2;
                if(ep_iterations > njobs && sp_iterations+1 > njobs)
                    return nullptr;
            }
            sol_i = newsol[start_position];
            newsol.erase(newsol.begin()+start_position);
        }
        else
        {
            sp_iterations++;
            ep_iterations = 1;
            start_position = ((start_position)%njobs)+1;
            if(end_position == start_position-1)
            {
                end_position=((end_position+1)%njobs)+1;
                ep_iterations+=2;
                if(ep_iterations > njobs && sp_iterations+1 > njobs)
                    return nullptr;
            }
            sol_i = newsol[start_position];
            newsol.erase(newsol.begin()+start_position);
#ifdef ENABLE_SSE
            computeHEAD(newsol,head,pmatrix,njobs-1,nmac);
#else
            computeHead(newsol);
#endif
        }

         newsol.insert(newsol.begin()+end_position,sol_i);

        //std::vector< int > ins_pos(nmac+1,0);
        int ins_pos[nmac+1];
        long int c_cur = head[1][end_position-1]+pmatrix[sol_i][1];
        ins_pos[1] = c_cur;

        for (int i = 2; i <= nmac; ++i) {
            int c_pm = head[i][end_position-1];
            if(c_pm > c_cur)
            {
                c_cur = c_pm;
            }
            c_cur += pmatrix[sol_i][i];
           ins_pos[i] =  c_cur;

        }
        long int pre_c_cur = c_cur;
        int wt = (std::max(c_cur - pis.getDueDate(sol_i), 0L) * pis.getPriority(sol_i));

        for (int j = 1; j< end_position; ++j )
        {
            wt += (std::max((long int)head[nmac][j] - pis.getDueDate(newsol[j]), 0L) * pis.getPriority(newsol[j]));
        }

        int pre_wt = wt;

        if(end_position > aptre)
        for(int k=end_position+1; k<= njobs; k++)
        {
            pre_c_cur = pre_c_cur + pmatrix[newsol[k]][nmac];
            wt += (std::max(pre_c_cur - pis.getDueDate(newsol[k]), 0L) * pis.getPriority(newsol[k]));
        }


        int value_wt = value->getSolutionValue();

        if(wt < value_wt)
        {
            //emili::iteration_decrement();
            for(int k=end_position+1; k<= njobs; k++)
            {
                int job = newsol[k];
                pre_c_cur = ins_pos[1] + pmatrix[job][1];
                ins_pos[1] = pre_c_cur;
                for(int m=2; m <= nmac ; m++)
                {
                    int c_pm = ins_pos[m];
                    if(c_pm > pre_c_cur)
                    {
                        pre_c_cur = c_pm;
                    }
                    pre_c_cur += pmatrix[job][m];
                    ins_pos[m] = pre_c_cur;
                }
               pre_wt += (std::max(pre_c_cur - pis.getDueDate(newsol[k]), 0L) * pis.getPriority(newsol[k]));
              if(pre_wt > value_wt)
                {
                   value->setSolutionValue(pre_wt);
                   return value;
                }

            }

            value->setSolutionValue(pre_wt);

        }
        else
        {
            value->setSolutionValue(wt);
        }

        return value;
    }
}

/**  Approximation based neighborhoods for other objectives ( no Weighted Tardiness)
 *
 * */

emili::Solution* emili::pfsp::NatxTCTNeighborhood::computeStep(emili::Solution *value)
{
    emili::iteration_increment();
    if(sp_iterations > njobs)
    {
        return nullptr;
    }
    else
    {
        end_position = ((end_position)%njobs)+1;
        std::vector < int >& newsol = ((emili::pfsp::PermutationFlowShopSolution*)value)->getJobSchedule();
        int sol_i;
        if(ep_iterations < njobs){
            ep_iterations++;
            /** if(ep_iterations == sp_iterations){
                ep_iterations++;
                end_position++;
            }*/
            if(end_position == start_position-1)
            {
                end_position=((end_position+1)%njobs)+1;
                ep_iterations+=2;
                if(ep_iterations > njobs && sp_iterations+1 > njobs)
                    return nullptr;
            }
            sol_i = newsol[start_position];
            newsol.erase(newsol.begin()+start_position);
        }
        else
        {
            sp_iterations++;
            ep_iterations = 1;
            start_position = ((start_position)%njobs)+1;
            if(end_position == start_position-1)
            {
                end_position=((end_position+1)%njobs)+1;
                ep_iterations+=2;
                if(ep_iterations > njobs && sp_iterations+1 > njobs)
                    return nullptr;
            }
            sol_i = newsol[start_position];
            newsol.erase(newsol.begin()+start_position);
#ifdef ENABLE_SSE
            computeHEAD(newsol,head,pmatrix,njobs-1,nmac);
#else
            computeHead(newsol);
#endif
        }

        newsol.insert(newsol.begin()+end_position,sol_i);
/** #ifdef ENABLE_SSE
        float ins_pos[nmac+1];
#else*/
        int ins_pos[nmac+1];
//#endif
        long int c_cur = head[1][end_position-1]+pmatrix[sol_i][1];
        ins_pos[1] = c_cur;

        for (int i = 2; i <= nmac; ++i) {
            int c_pm = head[i][end_position-1];
            if(c_pm > c_cur)
            {
                c_cur = c_pm;
            }
            c_cur += pmatrix[sol_i][i];
           ins_pos[i] =  c_cur;
        }
        long int pre_c_cur = c_cur;
        int wt = c_cur;

        for (int j = 1; j< end_position; ++j )
        {
            wt += head[nmac][j];
        }

        int pre_wt = wt;
        if(end_position > thresh)
        for(int k=end_position+1; k<= njobs; k++)
        {
            int nk = newsol[k];
            pre_c_cur = pre_c_cur + pmatrix[nk][nmac];
            wt += pre_c_cur;
        }

        if(wt < value_wt)
        {
            if(end_position > thresh)
            {
                thresh++;
            }
            /**
#ifdef ENABLE_SSE
            std::vector< long int > pmet(njobs+1,0);
            computePMakespans(newsol,pmet,pis.getProcessingTimesMatrix(),njobs+1,nmac,end_position+1,ins_pos);
           for(int k=end_position+1; k<= njobs; k++)
            {
               int job = newsol[k];
               pre_wt += (std::max(pmet[k] - pis.getDueDate(job), 0L) * pis.getPriority(job));
               if(pre_wt > value_wt)
                 {
                   break;
                 }
           }

#else
*/
            for(int k=end_position+1; k<= njobs; k++)
            {
                int job = newsol[k];
                pre_c_cur = ins_pos[1] + pmatrix[job][1];
                ins_pos[1] = pre_c_cur;
                //std::cout << "start critico " << std::endl;
                for(int m=2; m <= nmac ; m++)
                {
                    int c_pm = ins_pos[m];
                    if(c_pm > pre_c_cur)
                    {
                        pre_c_cur = c_pm;
                    }
                    pre_c_cur += pmatrix[job][m];
                    ins_pos[m] = pre_c_cur;
                }
            //    std::cout << "end critico " << std::endl;

               pre_wt += pre_c_cur;
              if(pre_wt > value_wt)
                {
                  // value->setSolutionValue(pre_wt);
                 //  double t = ((double)clock()-s)/CLOCKS_PER_SEC;
                 // std::cout << "natx2 : " << t << std::endl;
                  // return value;
                  break;
                }
            }
//#endif

            value->setSolutionValue(pre_wt);

        }
        else
        {
           thresh--;
            value->setSolutionValue(wt);
        }
        return value;
    }
}

/**  Approximation based neighborhoods for other objectives ( no Weighted Tardiness)
 *
 * */

emili::Solution* emili::pfsp::NrzTCTNeighborhood::computeStep(emili::Solution *value)
{
    emili::iteration_increment();
    if(sp_iterations > njobs)
    {
        return nullptr;
    }
    else
    {

    std::vector < int >& newsol = ((emili::pfsp::PermutationFlowShopSolution*)value)->getJobSchedule();
    int sol_i;

    sp_iterations++;
    //  ep_iterations = 1;
    start_position = ((start_position)%njobs)+1;
    int seed = seed_seq[start_position];
    int actual_pos = 1;
    for(int i=1;i<=njobs;i++)
        if(seed == newsol[i])
            actual_pos = i;

    sol_i = seed;
    newsol.erase(newsol.begin()+actual_pos);
#ifdef ENABLE_SSE
    computeHEAD(newsol,head,pmatrix,njobs-1,nmac);
#else
    computeHead(newsol);
#endif

    int best_wt =  std::numeric_limits< int >::max();
    int best_ins = 1;
    for(best_ins=1; best_ins <= njobs;best_ins++)
    {
        //  newsol.insert(newsol.begin()+end_position,sol_i);
        // compute Ct for sol_i
        int ins_pos[nmac+1];
        long int c_cur = head[1][best_ins-1]+pmatrix[sol_i][1];
        ins_pos[1] = c_cur;

        for (int i = 2; i <= nmac; ++i) {
            int c_pm = head[i][best_ins-1];
            if(c_pm > c_cur)
            {
                c_cur = c_pm;
            }
            c_cur += pmatrix[sol_i][i];
            ins_pos[i] =  c_cur;
        }

        int wt = c_cur;

        // add Ct for jobs before end_position
        for (int j = 1; j< best_ins; ++j )
        {
            wt += head[nmac][j];
        }

        //compute Ct for jobs after end_position
        for(int k=best_ins; k< njobs; k++)
        {
            int job = newsol[k];
            c_cur = ins_pos[1] + pmatrix[job][1];
            ins_pos[1] = c_cur;
            for(int m=2; m <= nmac ; m++)
            {
                int c_pm = ins_pos[m];
                if(c_pm > c_cur)
                {
                    c_cur = c_pm;
                }
                c_cur += pmatrix[job][m];
                ins_pos[m] = c_cur;
            }
            wt += c_cur;
        }

        if(wt < best_wt)
        {
            best_wt = wt;
            end_position = best_ins;
        }

    }
    newsol.insert(newsol.begin()+end_position,sol_i);
    value->setSolutionValue(best_wt);
    return value;
    }

}



emili::Solution* emili::pfsp::NatxTTNeighborhood::computeStep(emili::Solution *value)
{

    emili::iteration_increment();
    if(sp_iterations > njobs)
    {
        return nullptr;
    }
    else
    {
        end_position = ((end_position)%njobs)+1;
        std::vector < int >& newsol = ((emili::pfsp::PermutationFlowShopSolution*)value)->getJobSchedule();
        int sol_i;
        if(ep_iterations < njobs){
            ep_iterations++;
            /** if(ep_iterations == sp_iterations){
                ep_iterations++;
                end_position++;
            }*/
            if(end_position == start_position-1)
            {
                end_position=((end_position+1)%njobs)+1;
                ep_iterations+=2;
                if(ep_iterations > njobs && sp_iterations+1 > njobs)
                    return nullptr;
            }
            sol_i = newsol[start_position];
            newsol.erase(newsol.begin()+start_position);
        }
        else
        {
            sp_iterations++;
            ep_iterations = 1;
            start_position = ((start_position)%njobs)+1;
            if(end_position == start_position-1)
            {
                end_position=((end_position+1)%njobs)+1;
                ep_iterations+=2;
                if(ep_iterations > njobs && sp_iterations+1 > njobs)
                    return nullptr;
            }
            sol_i = newsol[start_position];
            newsol.erase(newsol.begin()+start_position);
#ifdef ENABLE_SSE
            computeHEAD(newsol,head,pmatrix,njobs-1,nmac);
#else
            computeHead(newsol);
#endif
        }

        newsol.insert(newsol.begin()+end_position,sol_i);
/** #ifdef ENABLE_SSE
        float ins_pos[nmac+1];
#else*/
        int ins_pos[nmac+1];
//#endif
        long int c_cur = head[1][end_position-1]+pmatrix[sol_i][1];
        ins_pos[1] = c_cur;

        for (int i = 2; i <= nmac; ++i) {
            int c_pm = head[i][end_position-1];
            if(c_pm > c_cur)
            {
                c_cur = c_pm;
            }
            c_cur += pmatrix[sol_i][i];
           ins_pos[i] =  c_cur;
        }
        long int pre_c_cur = c_cur;
        int wt = (std::max(c_cur - duedates[sol_i], 0L));

        for (int j = 1; j< end_position; ++j )
        {
            wt += (std::max((long int)head[nmac][j] - duedates[newsol[j]], 0L));
        }

        int pre_wt = wt;
        if(end_position > thresh)
        for(int k=end_position+1; k<= njobs; k++)
        {
            int nk = newsol[k];
            pre_c_cur = pre_c_cur + pmatrix[nk][nmac];
            wt += (std::max(pre_c_cur - duedates[nk], 0L));
        }

        if(wt < value_wt)
        {
            if(end_position > thresh)
            {
                thresh++;
            }

            for(int k=end_position+1; k<= njobs; k++)
            {
                int job = newsol[k];
                pre_c_cur = ins_pos[1] + pmatrix[job][1];
                ins_pos[1] = pre_c_cur;
                //std::cout << "start critico " << std::endl;
                for(int m=2; m <= nmac ; m++)
                {
                    int c_pm = ins_pos[m];
                    if(c_pm > pre_c_cur)
                    {
                        pre_c_cur = c_pm;
                    }
                    pre_c_cur += pmatrix[job][m];
                    ins_pos[m] = pre_c_cur;
                }
            //    std::cout << "end critico " << std::endl;

               pre_wt += (std::max(pre_c_cur - duedates[job], 0L));
              if(pre_wt > value_wt)
                {
                  // value->setSolutionValue(pre_wt);
                 //  double t = ((double)clock()-s)/CLOCKS_PER_SEC;
                 // std::cout << "natx2 : " << t << std::endl;
                  // return value;
                  break;
                }

            }
//#endif

            value->setSolutionValue(pre_wt);

        }
        else
        {

           thresh--;
            value->setSolutionValue(wt);
        }

        return value;
    }
}


/**
 * No Indle accelerated neighborhood
 **/
emili::Solution* emili::pfsp::NoIdleAcceleratedInsertNeighborhood::computeStep(emili::Solution *value)
{
    emili::iteration_increment();
    if(sp_iterations > njobs)
    {
        return nullptr;
    }
    else
    {
        end_position = ((end_position)%njobs)+1;
        std::vector < int >& newsol = ((emili::pfsp::PermutationFlowShopSolution*)value)->getJobSchedule();
        int sol_i;
        if(ep_iterations < njobs){
            ep_iterations++;
            if(end_position == start_position-1)
            {
                end_position=((end_position+1)%njobs)+1;
                ep_iterations+=2;
                if(ep_iterations > njobs && sp_iterations+1 > njobs)
                    return nullptr;
            }
            sol_i = newsol[start_position];
            newsol.erase(newsol.begin()+start_position);
        }
        else
        {
            sp_iterations++;
            ep_iterations = 1;
            start_position = ((start_position)%njobs)+1;
            if(end_position == start_position-1)
            {
                end_position=((end_position+1)%njobs)+1;
                ep_iterations+=2;
                if(ep_iterations > njobs && sp_iterations+1 > njobs)
                    return nullptr;
            }
            sol_i = newsol[start_position];
            newsol.erase(newsol.begin()+start_position);
            pis.computeNoIdleTAmatrices(newsol,head,tail);
        }




        newsol.insert(newsol.begin()+end_position,sol_i);
        long int c_cur = head[1][end_position-1]+pmatrix[sol_i][1];
        long int c_max = c_cur+tail[1][end_position];
        long int a = 0;
        long int aa = 0;
        for (int i = 2; i <= pis.getNmachines(); ++i) {
            int c_pm = head[i][end_position-1];
            if(c_pm + aa < c_cur)
            {
                aa += c_cur - (c_pm+aa);
                c_cur = c_cur + pmatrix[sol_i][i];
            }
            else
            {
                c_cur = c_pm + aa + pmatrix[sol_i][i];
            }
            long int c_can = (c_cur+a+tail[i][end_position]);
            c_max = c_max>c_can?c_max:c_can;
            a = a+c_max-c_can;
        }
        //long int old_v  = pis.computeMS(newsol);
       // long int old_v2 = pis.computeObjectiveFunction(newsol);
        //std::cout << c_max << " - " << old_v << std::endl;

        //assert(c_max == old_v2);
        value->setSolutionValue(c_max);
        return value;
    }
}
/*
  Sequence Dependent Setup Times Taillard Acceleration
*/

emili::Neighborhood::NeighborhoodIterator emili::pfsp::SDSTTaillardAcceleratedInsertNeighborhood::begin(Solution *base)
{
    ep_iterations = 1;
    sp_iterations = 1;
//    std::vector< int > sol(((emili::pfsp::PermutationFlowShopSolution*)base)->getJobSchedule());
//    sol.erase(sol.begin()+start_position);
//#ifdef ENABLE_SSE
//    pis.getInstance().computeSDSTTAmatrices(sol,head,tail,njobs);
    return emili::Neighborhood::NeighborhoodIterator(this,base);
}

emili::Solution* emili::pfsp::SDSTTaillardAcceleratedInsertNeighborhood::computeStep(emili::Solution *value)
{
    emili::iteration_increment();
    if(sp_iterations > njobs)
    {
        return nullptr;
    }
    else
    {
        std::vector < int >& newsol = ((emili::pfsp::PermutationFlowShopSolution*)value)->getJobSchedule();

        int best_inspos = 1;
        int best_cmax = std::numeric_limits<int>::max();
        end_position = 1;
        start_position = ((start_position)%njobs)+1;
        int sol_i = newsol[start_position];
        newsol.erase(newsol.begin()+start_position);
//#ifdef ENABLE_SSE
        pis.getInstance().computeSDSTTAmatrices(newsol,head,tail,njobs);

        for(int k=1; k<=njobs; k++)
        {
            int kp1 = newsol[k-1];
            if(k != start_position)
            {
                long int c_cur = head[1][k-1]+pmatrix[sol_i][1]+setUpTimes[1][kp1][sol_i];
                long int c_max = c_cur+tail[1][k]+setUpTimes[1][sol_i][newsol[k]];
                for (int i = 2; i <= nmac; ++i) {
                    int c_pm = setUpTimes[i][kp1][sol_i] + head[i][k-1];

                    if(c_pm > c_cur)
                    {
                        c_cur = c_pm;
                    }

                    c_cur = c_cur + pmatrix[sol_i][i];
                    long int c_can = (c_cur+tail[i][k]+setUpTimes[i][sol_i][newsol[k]]);
                    c_max = c_max>c_can?c_max:c_can;
                }

                if(c_max < best_cmax)
                {
                    best_cmax = c_max;
                    best_inspos = k;

                }

            }

        }
        end_position = best_inspos;
        newsol.insert(newsol.begin()+best_inspos,sol_i);
      //  long int old_v  = pis.computeObjectiveFunction(newsol);
        //std::cout << best_cmax << " - " << old_v << std::endl;
        //assert(best_cmax == old_v);
        value->setSolutionValue(best_cmax);
        sp_iterations++;
        return value;
    }
}
/*
 * No Wait Delta evaluation
 *
 */

emili::Solution* emili::pfsp::NoWaitAcceleratedInsertNeighborhood::computeStep(Solution *value)
{
    emili::iteration_increment();
    if(sp_iterations > njobs)
    {
        return nullptr;
    }
    else{
        end_position = ((end_position)%njobs)+1;
        if(ep_iterations < njobs){
            ep_iterations++;

            if(end_position == start_position-1)
            {
                end_position=((end_position+1)%njobs)+1;
                ep_iterations+=2;
                if(ep_iterations > njobs && sp_iterations+1 > njobs)
                    return nullptr;
            }

        }
        else
        {
            sp_iterations++;
            ep_iterations = 1;
            start_position = ((start_position)%njobs)+1;
            if(end_position == start_position-1)
            {
                end_position=((end_position+1)%njobs)+1;
                ep_iterations+=2;
                if(ep_iterations > njobs && sp_iterations+1 > njobs)
                    return nullptr;
            }

        }
        std::vector < int >& newsol = ((emili::pfsp::PermutationFlowShopSolution*)value)->getJobSchedule();
/*
        std::ostringstream oss;
        oss << "*********\n";
        oss << "sp " << start_position<<"\n";
        oss << "ep " << end_position<<"\n";
        for(int i=0; i <= njobs ; i++)
            oss << " " << newsol[i];
        oss << "\n";
*/
        int j = newsol[start_position];
        int jpo = newsol[(start_position+1)%(njobs+1)];
        int jmo = newsol[start_position-1];
        newsol.erase(newsol.begin()+start_position);
        int k,kmo;
        //end_position != njobs ? k = newsol[end_position] : k = 0;
        k = newsol[end_position%njobs];
        kmo = newsol[end_position-1];
        newsol.insert(newsol.begin()+end_position,j);        
        long int new_value = value->getSolutionValue();
        int delta = distance[kmo][j] + distance[j][k]
                  - distance[kmo][k] + distance[jmo][jpo]
                  - distance[jmo][j] - distance[j][jpo];
        new_value += delta;
        value->setSolutionValue(new_value);
/*
       long int old_value = pis.computeObjectiveFunction(newsol);
       if(new_value != old_value)
       {
        std::cout << oss.str();
        std::cout << "j " << j << "\n";
        std::cout << "jpo " << jpo << "\n";
        std::cout << "jmo " << jmo << "\n";
        std::cout << "k " << k << "\n";
        std::cout << "kmo " << kmo << "\n";
        std::cout << "new " << new_value << " = " << old_value << "\n";
        for(int i=0; i <= njobs ; i++)
            std::cout << " " << newsol[i];
        std::cout << "\n";
        std::cout << "+++++++++++\n";
        assert(new_value == old_value);
        }
*/
        return value;
    }
}

emili::Solution* emili::pfsp::NoWaitAcceleratedInsertNeighborhood::random(Solution *currentSolution)
{    
    std::vector < int > newsol(((emili::pfsp::PermutationFlowShopSolution*)currentSolution)->getJobSchedule());
    int sposition = (emili::generateRandomNumber()%njobs) +1;
    int eposition = (emili::generateRandomNumber()%njobs) +1;
    int j = newsol[sposition];
    int jpo = newsol[(sposition+1)%(njobs+1)];
    int jmo = newsol[sposition-1];
    newsol.erase(newsol.begin()+sposition);
    int k,kmo;
    //end_position != njobs ? k = newsol[end_position] : k = 0;
    k = newsol[eposition%njobs];
    kmo = newsol[eposition-1];
    newsol.insert(newsol.begin()+eposition,j);
    long int new_value = currentSolution->getSolutionValue();
    int delta = distance[kmo][j] + distance[j][k]
              - distance[kmo][k] + distance[jmo][jpo]
              - distance[jmo][j] - distance[j][jpo];
    new_value += delta;
    return new emili::pfsp::PermutationFlowShopSolution(new_value,newsol);
}

emili::Solution* emili::pfsp::NoWaitAcceleratedTwoInsertNeighborhood::computeStep(emili::Solution* value)
{

    emili::iteration_increment();
    if(sp_iterations >= njobs)
    {
        return nullptr;
    }
    else{
        end_position = ((end_position+1)%njobs);
        end_position==0?end_position=1:1+1;
        if(ep_iterations < njobs){
            ep_iterations++;
           /**  if(ep_iterations == sp_iterations){
                ep_iterations++;
                end_position++;
                std::cout << "BOOM!" << std::endl;
            }*/
            if(end_position == start_position-1)
            {
                end_position=((end_position+2)%njobs);
                end_position==0?end_position=1:1+1;
                ep_iterations+=2;
                if(ep_iterations > njobs && sp_iterations+1 > njobs)
                    return nullptr;
            }

        }
        else
        {
            sp_iterations++;
            ep_iterations = 1;
            start_position = ((start_position+1)%njobs);
            start_position==0?start_position=1:1+1;
            if(end_position == start_position-1)
            {
                end_position=((end_position+2)%njobs);
                end_position==0?end_position=1:1+1;
                ep_iterations+=2;
                if(ep_iterations > njobs && sp_iterations+1 > njobs)
                    return nullptr;
            }

        }

       std::vector < int >& newsol = ((emili::pfsp::PermutationFlowShopSolution*)value)->getJobSchedule();
/*
       std::cout << "++++++++++++\n";
       for(int i = 0 ; i <= njobs ; i++)
           std::cout << " " << newsol[i];
       std::cout << "\n";
*/
       int j = newsol[start_position];
       int j2 = newsol[start_position+1];
       int jpo = newsol[(start_position+2)%(njobs+1)];
       int jmo = newsol[start_position-1];
/*
       std::cout << "sol_i " << j << "\nsol_i2 " << j2 << "\n";
       std::cout << "startpos " << start_position << "\n";
       std::cout << "endpos " << end_position << "\n";
*/
       newsol.erase(newsol.begin()+start_position+1);
       newsol.erase(newsol.begin()+start_position);
       //std::cout << "epos " << epos << "\n";
       int k = newsol[end_position%(njobs-1)];
       int kmo = newsol[end_position-1];
       newsol.insert(newsol.begin()+end_position,j2);
       newsol.insert(newsol.begin()+end_position,j);
       long int new_value = value->getSolutionValue();
       int delta = distance[kmo][j] + distance[j2][k]
                 - distance[kmo][k] + distance[jmo][jpo]
                 - distance[jmo][j] - distance[j2][jpo];
       new_value += delta;
/*
       long int old_value = pis.computeObjectiveFunction(newsol);
       std::cout << "j " << j << "\n";
       std::cout << "j2 " << j2 << "\n";
       std::cout << "jpo " << jpo << "\n";
       std::cout << "jmo " << jmo << "\n";
       std::cout << "k " << k << "\n";
       std::cout << "kmo " << kmo << "\n";
       std::cout << "new " << new_value << " = " << old_value << "\n";
       for(int i = 0 ; i <= njobs ; i++)
           std::cout << " " << newsol[i];
       std::cout << "\n";
       std::cout << "*************\n";
       assert(old_value == new_value);
*/
       value->setSolutionValue(new_value);
       return value;
    }
}

emili::Solution* emili::pfsp::NoWaitAcceleratedTwoInsertNeighborhood::random(Solution *currentSolution)
{    
    std::vector < int > newsol(((emili::pfsp::PermutationFlowShopSolution*)currentSolution)->getJobSchedule());
    int sposition = (emili::generateRandomNumber()%(njobs));
    if(sposition==0)
        sposition=1;
    int eposition = (emili::generateRandomNumber()%(njobs-1));
    if(eposition==0)
        eposition=1;
    int j = newsol[sposition];
    int j2 = newsol[sposition+1];
    int jpo = newsol[(sposition+2)%(njobs+1)];
    int jmo = newsol[sposition-1];
    newsol.erase(newsol.begin()+sposition+1);
    newsol.erase(newsol.begin()+sposition);
    int k = newsol[eposition%(njobs-1)];
    int kmo = newsol[eposition-1];
    newsol.insert(newsol.begin()+eposition,j2);
    newsol.insert(newsol.begin()+eposition,j);
    long int new_value = currentSolution->getSolutionValue();
    int delta = distance[kmo][j] + distance[j2][k]
              - distance[kmo][k] + distance[jmo][jpo]
              - distance[jmo][j] - distance[j2][jpo];
    new_value += delta;
    return new PermutationFlowShopSolution(new_value,newsol);
}

void emili::pfsp::NoWaitAcceleratedTwoInsertNeighborhood::reverseLastMove(Solution *step)
{
    std::vector < int >& newsol = ((emili::pfsp::PermutationFlowShopSolution*)step)->getJobSchedule();
    int sol_i = newsol[end_position];
    int sol_i2 = newsol[end_position+1];
    newsol.erase(newsol.begin()+end_position+1);
    newsol.erase(newsol.begin()+end_position);
    newsol.insert(newsol.begin()+start_position,sol_i2);
    newsol.insert(newsol.begin()+start_position,sol_i);
}

void emili::pfsp::NoWaitAcceleratedExchangeNeighborhood::reverseLastMove(Solution *step)
{
    std::vector < int >& newsol = ((emili::pfsp::PermutationFlowShopSolution*)step)->getJobSchedule();
    std::swap(newsol[start_position],newsol[end_position]);
}

emili::Solution* emili::pfsp::NoWaitAcceleratedExchangeNeighborhood::computeStep(Solution* value)
{
    emili::iteration_increment();
     if(sp_iterations >= (njobs-1))
     {
         return nullptr;
     }
     else
     {
         if(ep_iterations < njobs){
             ep_iterations++;
         }
         else
         {
             sp_iterations++;
             ep_iterations = sp_iterations+1;
             start_position = (start_position%njobs)+1;
             end_position = start_position;
         }
         end_position = (end_position%njobs)+1;
         std::vector < int >& newsol = ((emili::pfsp::PermutationFlowShopSolution*)value)->getJobSchedule();
         int j = newsol[start_position];
         int jpo = newsol[(start_position+1)%(njobs+1)];
         int jmo = newsol[start_position-1];
         int k = newsol[end_position];
         int kpo = newsol[(end_position+1)%(njobs+1)];
         int kmo = newsol[end_position-1];
         long int new_value = value->getSolutionValue();
         int delta = distance[kmo][j] + distance[j][kpo]
                   + distance[jmo][k] + distance[k][jpo]
                   - distance[kmo][k] - distance[k][kpo]
                   - distance[jmo][j] - distance[j][jpo];
         if(kmo == j || kpo == j)
         {
             delta += distance[j][k]+distance[k][j];
         }
  /*        std::ostringstream oss;
        oss << "Cmax " << new_value << " delta " << delta << "\n";
         new_value += delta;
         for(int i=0; i <= njobs ; i++)
             oss << " " << newsol[i];
        oss << "\n";
*/
         std::swap(newsol[start_position],newsol[end_position]);
/*
         long int old_value = pis.computeObjectiveFunction(newsol);
         if(new_value != old_value)
         {
               std::cout << oss.str();
               std::cout << "j " << j << "\n";
               std::cout << "jpo " << jpo << "\n";
               std::cout << "jmo " << jmo << "\n";
               std::cout << "k " << k << "\n";
               std::cout << "kmo " << kmo << "\n";
               std::cout << "kpo " << kpo << "\n";

               std::cout << "new " << new_value << " = " << old_value << "\n";
               for(int i=0; i <= njobs ; i++)
                   std::cout << " " << newsol[i];
               std::cout << "\n";
               std::cout << "+++++++++++\n";
               assert(old_value==new_value);
         }
*/
         value->setSolutionValue(new_value);
         return value;
     }
}

emili::Solution* emili::pfsp::NoWaitAcceleratedExchangeNeighborhood::random(Solution *currentSolution)
{
    std::vector < int > newsol(((emili::pfsp::PermutationFlowShopSolution*)currentSolution)->getJobSchedule());
    int sposition = (emili::generateRandomNumber()%njobs)+1;
    int eposition = (emili::generateRandomNumber()%njobs)+1;

    int j = newsol[sposition];
    int jpo = newsol[(sposition+1)%(njobs+1)];
    int jmo = newsol[sposition-1];
    int k = newsol[eposition];
    int kpo = newsol[(eposition+1)%(njobs+1)];
    int kmo = newsol[eposition-1];
    long int new_value = currentSolution->getSolutionValue();
    int delta = distance[kmo][j] + distance[j][kpo]
              + distance[jmo][k] + distance[k][jpo]
              - distance[kmo][k] - distance[k][kpo]
              - distance[jmo][j] - distance[j][jpo];
    if(kmo == j || kpo == j)
    {
        delta += distance[j][k]+distance[k][j];
    }
    std::swap(newsol[sposition],newsol[eposition]);
    new_value += delta;
   return new PermutationFlowShopSolution(new_value,newsol);
}

void emili::pfsp::NoWaitAcceleratedTransposeNeighborhood::reverseLastMove(Solution *step)
{
    std::vector < int >& newsol = ((emili::pfsp::PermutationFlowShopSolution*)step)->getJobSchedule();
    int endpos = start_position<njobs?start_position+1:1;
    std::swap(newsol[start_position],newsol[endpos]);
}

emili::Solution* emili::pfsp::NoWaitAcceleratedTransposeNeighborhood::computeStep(Solution *value)
{
    emili::iteration_increment();
    if(sp_iterations >= njobs)
    {
        return nullptr;
    }
    else
    {
        sp_iterations++;
        start_position = (start_position%njobs)+1;
        std::vector < int >& newsol = ((emili::pfsp::PermutationFlowShopSolution*)value)->getJobSchedule();
        int end_position = start_position<njobs?start_position+1:1;
        int j = newsol[start_position];
        int jpo = newsol[(start_position+1)%(njobs+1)];
        int jmo = newsol[start_position-1];
        int k = newsol[end_position];
        int kpo = newsol[(end_position+1)%(njobs+1)];
        int kmo = newsol[end_position-1];
        long int new_value = value->getSolutionValue();
        int delta = distance[kmo][j] + distance[j][kpo]
                  + distance[jmo][k] + distance[k][jpo]
                  - distance[kmo][k] - distance[k][kpo]
                  - distance[jmo][j] - distance[j][jpo];
        if(kmo == j)
        {
            delta += distance[j][k]+distance[k][j];
        }
        //int delta2 = distance[jmo][k]+distance[k][j]+distance[j][kpo]
        //          - distance[jmo][j]-distance[j][k]-distance[k][kpo];
        new_value += delta;
        std::swap(newsol[start_position],newsol[end_position]);        
        value->setSolutionValue(new_value);
        return value;
    }
}

emili::Solution* emili::pfsp::NoWaitAcceleratedTransposeNeighborhood::random(Solution *currentSolution)
{
    std::vector < int > newsol(((emili::pfsp::PermutationFlowShopSolution*)currentSolution)->getJobSchedule());
    int sposition = (emili::generateRandomNumber()%njobs)+1;
    int eposition = sposition<njobs?sposition+1:1;
    int j = newsol[sposition];
    int jpo = newsol[(sposition+1)%(njobs+1)];
    int jmo = newsol[sposition-1];
    int k = newsol[eposition];
    int kpo = newsol[(eposition+1)%(njobs+1)];
    int kmo = newsol[eposition-1];
    long int new_value = currentSolution->getSolutionValue();
    int delta = distance[kmo][j] + distance[j][kpo]
              + distance[jmo][k] + distance[k][jpo]
              - distance[kmo][k] - distance[k][kpo]
              - distance[jmo][j] - distance[j][jpo];
    if(kmo == j || kpo == j)
    {
        delta += distance[j][k]+distance[k][j];
    }
    std::swap(newsol[sposition],newsol[eposition]);
    new_value += delta;
    return new PermutationFlowShopSolution(new_value,newsol);
}



emili::Solution* emili::pfsp::TAxInsertNeighborhood::computeStep(emili::Solution *value)
{
    emili::iteration_increment();
    if(sp_iterations >= njobs)
    {
        return nullptr;
    }
    else
    {
        std::vector < int >& newsol = ((emili::pfsp::PermutationFlowShopSolution*)value)->getJobSchedule();
        int sol_i;
        if(ep_iterations < njobs){
            ep_iterations++;
            if(ep_iterations == sp_iterations){
                ep_iterations++;
                end_position++;
            }
            sol_i = newsol[start_position];
            newsol.erase(newsol.begin()+start_position);
        }
        else
        {
            sp_iterations++;
            ep_iterations = 1;
            start_position = ((start_position)%njobs)+1;
            sol_i = newsol[start_position];
            newsol.erase(newsol.begin()+start_position);
            pis.computeTAmatrices(newsol,head,tails[njobs]);
            pis.computeTails(newsol,tails);


        }
        end_position = ((end_position)%njobs)+1;

        newsol.insert(newsol.begin()+end_position,sol_i);

        //std::vector< int > ins_pos(nmac+1,0);
        int ins_pos[nmac+1];
        long int c_cur = head[1][end_position-1]+pmatrix[sol_i][1];
        ins_pos[1] = c_cur;
        for (int i = 2; i <= nmac; ++i) {
            int c_pm = head[i][end_position-1];
            if(c_pm < c_cur)
            {
                c_cur = c_cur + pmatrix[sol_i][i];
            }
            else
            {
                c_cur = c_pm + pmatrix[sol_i][i];
            }
            ins_pos[i] =  c_cur;
        }



        int wt = (std::max(c_cur - pis.getDueDate(newsol[end_position]), 0L) * pis.getPriority(newsol[end_position]));
        for(int k=end_position+1; k<= njobs; k++)
        {
            std::vector < std::vector < int > >& tail = tails[k];

            long int c_max = ins_pos[1]+tail[1][end_position];
            for (int i = 2; i <= pis.getNmachines(); ++i) {
                long int c_can = (ins_pos[i]+tail[i][end_position]);
                c_max = c_max>c_can?c_max:c_can;
            }

            wt += (std::max(c_max - pis.getDueDate(newsol[k]), 0L) * pis.getPriority(newsol[k]));
        }



        for (int j = 1; j< end_position; ++j )
            wt += (std::max((long int)head[nmac][j] - pis.getDueDate(newsol[j]), 0L) * pis.getPriority(newsol[j]));


        //long int old_v  = pis.computeObjectiveFunction(newsol);
        //std::cout << c_max << " - " << old_v << std::endl;
      //  assert(wt == old_v);
     //   int robj = pis.computeObjectiveFunction(newsol);
        //std::cout << "wt " << wt << " <-> "<< robj << std::endl;

       // assert(wt==robj);
        value->setSolutionValue(wt);
        return value;
    }
}

emili::Solution* emili::pfsp::PfspTwoInsertNeighborhood::computeStep(emili::Solution* value)
{

    emili::iteration_increment();
    if(sp_iterations >= njobs)
    {
        return nullptr;
    }
    else{
        end_position = ((end_position+1)%njobs);
        end_position==0?end_position=1:1+1;
        if(ep_iterations < njobs){
            ep_iterations++;
           /**  if(ep_iterations == sp_iterations){
                ep_iterations++;
                end_position++;
                std::cout << "BOOM!" << std::endl;
            }*/
            if(end_position == start_position-1)
            {
                end_position=((end_position+2)%njobs);
                end_position==0?end_position=1:1+1;
                ep_iterations+=2;
                if(ep_iterations > njobs && sp_iterations+1 > njobs)
                    return nullptr;
            }

        }
        else
        {
            sp_iterations++;
            ep_iterations = 1;
            start_position = ((start_position+1)%njobs);
            start_position==0?start_position=1:1+1;
            if(end_position == start_position-1)
            {
                end_position=((end_position+2)%njobs);
                end_position==0?end_position=1:1+1;
                ep_iterations+=2;
                if(ep_iterations > njobs && sp_iterations+1 > njobs)
                    return nullptr;
            }

        }

       std::vector < int >& newsol = ((emili::pfsp::PermutationFlowShopSolution*)value)->getJobSchedule();
  /*     std::cout << "++++++++++++\n";
       for(int i = 0 ; i <= njobs ; i++)
           std::cout << " " << newsol[i];
       std::cout << "\n";*/
       // first insert
       int sol_i = newsol[start_position];
       int sol_i2 = newsol[start_position+1];
 /*      std::cout << "sol_i " << sol_i << "\nsol_i2 " << sol_i2 << "\n";
       std::cout << "startpos" << start_position << "\n";
       std::cout << "endpos " << end_position << "\n";
 */    newsol.erase(newsol.begin()+start_position+1);
       newsol.erase(newsol.begin()+start_position);
       newsol.insert(newsol.begin()+end_position,sol_i2);
       newsol.insert(newsol.begin()+end_position,sol_i);
/*       for(int i = 0 ; i <= njobs ; i++)
           std::cout << " " << newsol[i];
       std::cout << "\n";
       std::cout << "*************\n";*/
       long int new_value = pis.computeObjectiveFunction(newsol);
       value->setSolutionValue(new_value);
       return value;
    }
}

void emili::pfsp::PfspTwoInsertNeighborhood::reverseLastMove(Solution *step)
{
    std::vector < int >& newsol = ((emili::pfsp::PermutationFlowShopSolution*)step)->getJobSchedule();
    int sol_i = newsol[end_position];
    int sol_i2 = newsol[end_position+1];
    newsol.erase(newsol.begin()+end_position+1);
    newsol.erase(newsol.begin()+end_position);
    newsol.insert(newsol.begin()+start_position,sol_i2);
    newsol.insert(newsol.begin()+start_position,sol_i);
}

emili::Solution* emili::pfsp::PfspInsertNeighborhood::random(Solution *currentSolution,int size)
{
    if(size > 0)
    {
    std::vector < int > newsol(((emili::pfsp::PermutationFlowShopSolution*)currentSolution)->getJobSchedule());    
    int njobs = pis.getNjobs();    
    std::vector < int > indexes;
    std::vector < int > kn(njobs+1,0);
    if(size >= njobs)
    {
        size = njobs-1;
    }

    int k=0;
    while(k <= size)
    {
        int ind = (emili::generateRandomNumber()%njobs)+1;
        if(kn[ind] == 0)
        {
            kn[ind] = 1;
            indexes.push_back(ind);
            k++;
        }
    }
    int i = 0;
    for(; i < size; i+=2)
    {
        int sol_i = newsol[indexes[i]];
        newsol.erase(newsol.begin()+indexes[i]);
        newsol.insert(newsol.begin()+indexes[i+1],sol_i);
    }
    i++;
    if(i < size)
    {
        int sol_i = newsol[indexes[i]];
        newsol.erase(newsol.begin()+indexes[i]);
        newsol.insert(newsol.begin()+indexes[0],sol_i);
    }
        long int value = pis.computeObjectiveFunction(newsol);
        return new emili::pfsp::PermutationFlowShopSolution(value,newsol);
    }
    return currentSolution->clone();
}

emili::Solution* emili::pfsp::PfspInsertNeighborhood::random(Solution *currentSolution)
{

    std::vector < int > newsol(((emili::pfsp::PermutationFlowShopSolution*)currentSolution)->getJobSchedule());
    int njobs = pis.getNjobs();

    int best_i = (emili::generateRandomNumber()%njobs)+1;
    int best_j = (emili::generateRandomNumber()%njobs)+1;
    int sol_i = newsol[best_i];
    newsol.erase(newsol.begin()+best_i);
    newsol.insert(newsol.begin()+best_j,sol_i);
    long int value = pis.computeObjectiveFunction(newsol);
    return new emili::pfsp::PermutationFlowShopSolution(value,newsol);
}

emili::Solution* emili::pfsp::KarNeighborhood::random(Solution *currentSolution)
{
    std::vector < int > newsol(((emili::pfsp::PermutationFlowShopSolution*)currentSolution)->getJobSchedule());
    int njobs = pis.getNjobs();
    float p = emili::generateRealRandomNumber();
    int best_i = (emili::generateRandomNumber()%njobs)+1;
    int best_j = (emili::generateRandomNumber()%njobs)+1;
    if(p < 0.5)
    {
        int sol_i = newsol[best_i];
        newsol.erase(newsol.begin()+best_i);
        newsol.insert(newsol.begin()+best_j,sol_i);
    }
    else
    {
        std::swap(newsol[best_i],newsol[best_j]);
    }
    long int value = pis.computeObjectiveFunction(newsol);
    return new emili::pfsp::PermutationFlowShopSolution(value,newsol);
}

emili::Solution* emili::pfsp::PfspForwardInsertNeighborhood::random(Solution *currentSolution)
{
    int njobs = pis.getNjobs();
    std::vector < int > newsol(((emili::pfsp::PermutationFlowShopSolution*)currentSolution)->getJobSchedule());
    int best_i = (emili::generateRandomNumber()%njobs)+1;
    int best_j = (emili::generateRandomNumber()%(njobs-best_i+1))+best_i;
    int sol_i = newsol[best_i];
    newsol.erase(newsol.begin()+best_i);
    newsol.insert(newsol.begin()+best_j,sol_i);
    long int value = pis.computeObjectiveFunction(newsol);
    return new emili::pfsp::PermutationFlowShopSolution(value,newsol);
}

emili::Solution* emili::pfsp::PfspBackwardInsertNeighborhood::random(Solution *currentSolution)
{
    int njobs = pis.getNjobs();
    std::vector < int > newsol(((emili::pfsp::PermutationFlowShopSolution*)currentSolution)->getJobSchedule());
    int best_i = (emili::generateRandomNumber()%njobs)+1;
    int best_j = (emili::generateRandomNumber()%best_i)+1;
    int sol_i = newsol[best_i];
    newsol.erase(newsol.begin()+best_i);
    newsol.insert(newsol.begin()+best_j,sol_i);
    long int value = pis.computeObjectiveFunction(newsol);
    return new emili::pfsp::PermutationFlowShopSolution(value,newsol);
}

emili::Solution* emili::pfsp::PfspTwoInsertNeighborhood::random(Solution *currentSolution)
{
    std::vector < int > newsol(((emili::pfsp::PermutationFlowShopSolution*)currentSolution)->getJobSchedule());
    int njobs = pis.getNjobs();
    int best_i = (emili::generateRandomNumber()%njobs)+1;
    int best_j = (emili::generateRandomNumber()%njobs)+1;

    //first insert
    int sol_i = newsol[best_i];
    newsol.erase(newsol.begin()+best_i);
    newsol.insert(newsol.begin()+best_j,sol_i);
    //second insert
    sol_i = newsol[best_i];
    newsol.erase(newsol.begin()+best_i);
    newsol.insert(newsol.begin()+best_j,sol_i);

    long int value = pis.computeObjectiveFunction(newsol);
    return new emili::pfsp::PermutationFlowShopSolution(value,newsol);
}

emili::Solution* emili::pfsp::PfspExchangeNeighborhood::computeStep(emili::Solution* value)
{
   emili::iteration_increment();
    if(sp_iterations >= (njobs-1))
    {        
        return nullptr;
    }
    else
    {
        if(ep_iterations < njobs){
            ep_iterations++;
        }
        else
        {
            sp_iterations++;
            ep_iterations = sp_iterations+1;
            start_position = (start_position%njobs)+1;
            end_position = start_position;
        }
        end_position = (end_position%njobs)+1;
        std::vector < int >& newsol = ((emili::pfsp::PermutationFlowShopSolution*)value)->getJobSchedule();
        std::swap(newsol[start_position],newsol[end_position]);
        long int new_value = pis.computeObjectiveFunction(newsol);
        value->setSolutionValue(new_value);
        return value;
    }
}

void emili::pfsp::PfspExchangeNeighborhood::reverseLastMove(Solution *step)
{
    std::vector < int >& newsol = ((emili::pfsp::PermutationFlowShopSolution*)step)->getJobSchedule();
    std::swap(newsol[start_position],newsol[end_position]);
}

void emili::pfsp::AxtExchange::computeHead(std::vector<int>& sol)
{
    int j,m;
    int jobNumber;
    int prevj = 0;
    for(j=1;j<njobs;j++)
    {
        jobNumber = sol[j];
        prevj = prevj + pmatrix[jobNumber][1];
        head[1][j] = prevj;
    }

      for ( j = 1; j < njobs; ++j )
        {
            long int previousJobEndTime = head[1][j];
            jobNumber = sol[j];

            for ( m = 2; m <= nmac; ++m )
            {
                if ( head[m][j-1] > previousJobEndTime )
                {
                    head[m][j] = head[m][j-1] + pmatrix[jobNumber][m];
                }
                else
                {
                    head[m][j] = previousJobEndTime + pmatrix[jobNumber][m];
                }
                previousJobEndTime = head[m][j];
            }
    }
}

emili::Solution* emili::pfsp::AxtExchange::computeStep(emili::Solution *value)
{
    emili::iteration_increment();
    if(sp_iterations >= (njobs-1))
    {
        return nullptr;
    }
    else
    {
        std::vector < int >& newsol = ((emili::pfsp::PermutationFlowShopSolution*)value)->getJobSchedule();
        int sol_i;
        if(ep_iterations < njobs){
            ep_iterations++;
        }
        else
        {
            sp_iterations++;
            ep_iterations = sp_iterations+1;
            start_position = (start_position%njobs)+1;
            end_position = start_position;
            sol_i = newsol[start_position];
            newsol.erase(newsol.begin()+start_position);
#ifdef ENABLE_SSE
            computeHEAD(newsol,head,pmatrix,njobs-1,nmac);
#else
            computeHead(newsol);
#endif
            newsol.insert(newsol.begin()+start_position,sol_i);
        }
        end_position = ((end_position)%njobs)+1;

        std::swap(newsol[start_position],newsol[end_position]);
        //std::vector< int > ins_pos(nmac+1,0);

        int ms_pos = start_position>end_position?end_position:start_position;
        sol_i = newsol[ms_pos];

        int ins_pos[nmac+1];
        long int c_cur = head[1][ms_pos-1]+pmatrix[sol_i][1];
        ins_pos[1] = c_cur;

        for (int i = 2; i <= nmac; ++i) {
            int c_pm = head[i][ms_pos-1];
            if(c_pm > c_cur)
            {
                c_cur = c_pm;
            }
            c_cur += pmatrix[sol_i][i];
            ins_pos[i] =  c_cur;

        }

        int wt = (std::max(c_cur - pis.getDueDate(sol_i), 0L) * pis.getPriority(sol_i));


        for (int j = 1; j< ms_pos; ++j )
        {
            wt += (std::max((long int)head[nmac][j] - pis.getDueDate(newsol[j]), 0L) * pis.getPriority(newsol[j]));
        }

        int pre_wt = wt;

        if(ms_pos > thresh)
        {
            long int pppre_c_cur = ins_pos[nmac-2];
            long int ppre_c_cur = ins_pos[nmac-1];
            long int pre_c_cur = c_cur;
            for(int k=ms_pos+1; k<= njobs; k++)
            {
                int nk = newsol[k];
                pppre_c_cur = pppre_c_cur + pmatrix[nk][nmac-2];
                ppre_c_cur = std::max(ppre_c_cur,pppre_c_cur) + pmatrix[nk][nmac-1];
                pre_c_cur = std::max(pre_c_cur,ppre_c_cur) + pmatrix[nk][nmac];
                wt += (std::max(pre_c_cur - pis.getDueDate(nk), 0L) * pis.getPriority(nk));
            }
        }

        int value_wt = value->getSolutionValue();


        if(wt < value_wt)
        {

             //emili::iteration_increment();
            for(int k=ms_pos+1; k<= njobs; k++)
            {
                int job = newsol[k];
                c_cur = ins_pos[1] + pmatrix[job][1];
                ins_pos[1] = c_cur;
                for(int m=2; m <= nmac ; m++)
                {
                    int c_pm = ins_pos[m];
                    if(c_pm > c_cur)
                    {
                        c_cur = c_pm;
                    }
                    c_cur += pmatrix[job][m];
                    ins_pos[m] = c_cur;
                }
                pre_wt += (std::max(c_cur - pis.getDueDate(job), 0L) * pis.getPriority(job));
               if(pre_wt > value_wt)
                {
                    value->setSolutionValue(pre_wt);
                    return value;
                }

            }
            value->setSolutionValue(pre_wt);

        }
        else
        {


            value->setSolutionValue(wt);
        }


        //long int old_v  = pis.computeObjectiveFunction(newsol);
        //std::cout << c_max << " - " << old_v << std::endl;
        //assert(c_max == old_v);
        return value;
    }
}


emili::Solution* emili::pfsp::HaxtExchange::computeStep(emili::Solution *value)
{
    emili::iteration_increment();
    if(sp_iterations >= (njobs-1))
    {
        return nullptr;
    }
    else
    {
        std::vector < int >& newsol = ((emili::pfsp::PermutationFlowShopSolution*)value)->getJobSchedule();
        int sol_i;
        if(ep_iterations < njobs){
            ep_iterations++;
        }
        else
        {
            sp_iterations++;
            ep_iterations = sp_iterations+1;
            start_position = (start_position%njobs)+1;
            end_position = start_position;
            sol_i = newsol[start_position];
            newsol.erase(newsol.begin()+start_position);
#ifdef ENABLE_SSE
            computeHEAD(newsol,head,pmatrix,njobs-1,nmac);
#else
            computeHead(newsol);
#endif
            newsol.insert(newsol.begin()+start_position,sol_i);
        }
        end_position = ((end_position)%njobs)+1;

        std::swap(newsol[start_position],newsol[end_position]);
        //std::vector< int > ins_pos(nmac+1,0);

        int ms_pos = start_position>end_position?end_position:start_position;
        sol_i = newsol[ms_pos];

        int ins_pos[nmac+1];
        long int c_cur = head[1][ms_pos-1]+pmatrix[sol_i][1];
        ins_pos[1] = c_cur;

        for (int i = 2; i <= nmac; ++i) {
            int c_pm = head[i][ms_pos-1];
            if(c_pm > c_cur)
            {
                c_cur = c_pm;
            }
            c_cur += pmatrix[sol_i][i];
            ins_pos[i] =  c_cur;

        }

        int wt = (std::max(c_cur - pis.getDueDate(sol_i), 0L) * pis.getPriority(sol_i));
        long int pre_c_cur = c_cur;

        for (int j = 1; j< ms_pos; ++j )
        {
            wt += (std::max((long int)head[nmac][j] - pis.getDueDate(newsol[j]), 0L) * pis.getPriority(newsol[j]));
        }

        int pre_wt = wt;

        if(ms_pos > threshold)
            for(int k=ms_pos+1; k<= njobs; k++)
            {
                pre_c_cur = pre_c_cur + pmatrix[newsol[k]][nmac];
                wt += (std::max(pre_c_cur - pis.getDueDate(newsol[k]), 0L) * pis.getPriority(newsol[k]));
            }

        int value_wt = value->getSolutionValue();


        if(wt < value_wt)
        {
           // if(ms_pos > threshold)
          //      threshold++;
            for(int k=ms_pos+1; k<= njobs; k++)
            {
                int job = newsol[k];
                c_cur = ins_pos[1] + pmatrix[job][1];
                ins_pos[1] = c_cur;
                for(int m=2; m <= nmac ; m++)
                {
                    int c_pm = ins_pos[m];
                    if(c_pm > c_cur)
                    {
                        c_cur = c_pm;
                    }
                    c_cur += pmatrix[job][m];
                    ins_pos[m] = c_cur;
                }
                pre_wt += (std::max(c_cur - pis.getDueDate(job), 0L) * pis.getPriority(job));
               if(pre_wt > value_wt)
                {
                    value->setSolutionValue(pre_wt);
                    return value;
                }

            }
                     value->setSolutionValue(pre_wt);

        }
        else
        {
        //    threshold--;
            value->setSolutionValue(wt);
        }
        return value;
    }
}

emili::Solution* emili::pfsp::EaxtExchange::computeStep(emili::Solution *value)
{
    emili::iteration_increment();
    if(sp_iterations >= (njobs-1))
    {
        return nullptr;
    }
    else
    {
        std::vector < int >& newsol = ((emili::pfsp::PermutationFlowShopSolution*)value)->getJobSchedule();
        int sol_i;
        if(ep_iterations < njobs){
            ep_iterations++;
        }
        else
        {
            sp_iterations++;
            ep_iterations = sp_iterations+1;
            start_position = (start_position%njobs)+1;
            end_position = start_position;
            sol_i = newsol[start_position];
            newsol.erase(newsol.begin()+start_position);
#ifdef ENABLE_SSE
            computeHEAD(newsol,head,pmatrix,njobs-1,nmac);
#else
            computeHead(newsol);
#endif
            newsol.insert(newsol.begin()+start_position,sol_i);
        }
        end_position = ((end_position)%njobs)+1;

        std::swap(newsol[start_position],newsol[end_position]);
        //std::vector< int > ins_pos(nmac+1,0);

        int ms_pos = start_position>end_position?end_position:start_position;
        sol_i = newsol[ms_pos];

        int ins_pos[nmac+1];
        long int c_cur = head[1][ms_pos-1]+pmatrix[sol_i][1];
        ins_pos[1] = c_cur;

        for (int i = 2; i <= nmac; ++i) {
            int c_pm = head[i][ms_pos-1];
            if(c_pm > c_cur)
            {
                c_cur = c_pm;
            }
            c_cur += pmatrix[sol_i][i];
            ins_pos[i] =  c_cur;

        }

        int wt = (std::max(c_cur - pis.getDueDate(sol_i), 0L) * pis.getPriority(sol_i));
        //long int pppre_c_cur = ins_pos[nmac-2];
        long int ppre_c_cur = ins_pos[nmac-1];
        long int pre_c_cur = c_cur;

        for (int j = 1; j< ms_pos; ++j )
        {
            wt += (std::max((long int)head[nmac][j] - pis.getDueDate(newsol[j]), 0L) * pis.getPriority(newsol[j]));
        }

        int pre_wt = wt;

        if(ms_pos > (njobs/2+10))
            for(int k=ms_pos+1; k<= njobs; k++)
            {
           //     pppre_c_cur = pppre_c_cur + pmatrix[newsol[k]][nmac-2];
            //    ppre_c_cur = std::max(ppre_c_cur,pppre_c_cur) + pmatrix[newsol[k]][nmac-1];
                    ppre_c_cur = ppre_c_cur + pmatrix[newsol[k]][nmac-1];
                pre_c_cur = std::max(pre_c_cur,ppre_c_cur) + pmatrix[newsol[k]][nmac];
                wt += (std::max(pre_c_cur - pis.getDueDate(newsol[k]), 0L) * pis.getPriority(newsol[k]));
            }

        int value_wt = value->getSolutionValue();


        if(wt < value_wt)
        {
             //emili::iteration_decrement();
            for(int k=ms_pos+1; k<= njobs; k++)
            {
                int job = newsol[k];
                pre_c_cur = ins_pos[1] + pmatrix[job][1];
                ins_pos[1] = pre_c_cur;
                for(int m=2; m <= nmac ; m++)
                {
                    int c_pm = ins_pos[m];
                    if(c_pm > pre_c_cur)
                    {
                        pre_c_cur = c_pm;
                    }
                    pre_c_cur += pmatrix[job][m];
                    ins_pos[m] = pre_c_cur;
                }
                pre_wt += (std::max(pre_c_cur - pis.getDueDate(newsol[k]), 0L) * pis.getPriority(newsol[k]));
               if(pre_wt > value_wt)
                {
                    value->setSolutionValue(pre_wt);
                    return value;
                }

            }
            value->setSolutionValue(pre_wt);

        }
        else
        {
            value->setSolutionValue(wt);
        }


        //long int old_v  = pis.computeObjectiveFunction(newsol);
        //std::cout << c_max << " - " << old_v << std::endl;
        //assert(c_max == old_v);
        return value;
    }
}



emili::Solution* emili::pfsp::OptExchange::computeStep(emili::Solution *value)
{
    emili::iteration_increment();
    if(sp_iterations >= (njobs-1))
    {
        return nullptr;
    }
    else
    {
       std::vector < int >& newsol = ((emili::pfsp::PermutationFlowShopSolution*)value)->getJobSchedule();
        int sol_i;
        if(ep_iterations < njobs){
            ep_iterations++;
        }
        else
        {
            sp_iterations++;
            ep_iterations = sp_iterations+1;
            start_position = (start_position%njobs)+1;
            end_position = start_position;
            sol_i = newsol[start_position];
            newsol.erase(newsol.begin()+start_position);
#ifdef ENABLE_SSE
            computeHEAD(newsol,head,pmatrix,njobs-1,nmac);
#else
            computeHead(newsol);
#endif
            newsol.insert(newsol.begin()+start_position,sol_i);
        }
        end_position = ((end_position)%njobs)+1;

        std::swap(newsol[start_position],newsol[end_position]);

        int ms_pos = start_position>end_position?end_position:start_position;
        sol_i = newsol[ms_pos];
        int ins_pos[nmac+1];
        long int c_cur = head[1][ms_pos-1]+pmatrix[sol_i][1];
        ins_pos[1] = c_cur;

        for (int i = 2; i <= nmac; ++i) {
            int c_pm = head[i][ms_pos-1];
            if(c_pm < c_cur)
            {
                c_cur = c_cur + pmatrix[sol_i][i];
            }
            else
            {
                c_cur = c_pm + pmatrix[sol_i][i];
            }
            ins_pos[i] =  c_cur;

        }

        int wt = (std::max(c_cur - pis.getDueDate(sol_i), 0L) * pis.getPriority(sol_i));

        long int pre_c_cur = c_cur;

        for (int j = 1; j< ms_pos; ++j )
        {
            wt += (std::max((long int)head[nmac][j] - pis.getDueDate(newsol[j]), 0L) * pis.getPriority(newsol[j]));
        }

        for(int k=ms_pos+1; k<= njobs; k++)
        {
            int job = newsol[k];
            pre_c_cur = ins_pos[1] + pmatrix[job][1];
            ins_pos[1] = pre_c_cur;
            for(int m=2; m <= nmac ; m++)
            {
                int c_pm = ins_pos[m];
                if(c_pm > pre_c_cur)
                {
                    pre_c_cur = c_pm;
                }
                pre_c_cur += pmatrix[job][m];
                ins_pos[m] = pre_c_cur;
            }
            wt += (std::max(pre_c_cur - pis.getDueDate(newsol[k]), 0L) * pis.getPriority(newsol[k]));
        }
        value->setSolutionValue(wt);


        return value;
    }
}


emili::Solution* emili::pfsp::PfspExchangeNeighborhood::random(Solution *currentSolution)
{

    std::vector < int > newsol(((emili::pfsp::PermutationFlowShopSolution*)currentSolution)->getJobSchedule());
    int njobs = pis.getNjobs();
    int best_i = (emili::generateRandomNumber()%njobs)+1;
    int best_j = (emili::generateRandomNumber()%njobs)+1;
    std::swap(newsol[best_i],newsol[best_j]);
    long int value = pis.computeObjectiveFunction(newsol);
    return new emili::pfsp::PermutationFlowShopSolution(value,newsol);
}

emili::Solution* emili::pfsp::PfspExchangeNeighborhood::random(Solution *currentSolution,int size)
{
    if(size > 0)
    {
    std::vector < int > newsol(((emili::pfsp::PermutationFlowShopSolution*)currentSolution)->getJobSchedule());
    int njobs = pis.getNjobs();
    std::vector < int > indexes;
    std::vector < int > kn(njobs+1,0);
    if(size >= njobs)
    {
        size = njobs-1;
    }

    int k=0;
    while(k <= size)
    {
        int ind = (emili::generateRandomNumber()%njobs)+1;
        if(kn[ind] == 0)
        {
            kn[ind] = 1;
            indexes.push_back(ind);
            k++;
        }
    }
    int i = 0;
    for(; i < size; i+=2)
    {
        std::swap(newsol[indexes[i]],newsol[indexes[i+1]]);
    }
    i++;
    if(i < size)
    {
       std::swap(newsol[indexes[i]],newsol[indexes[0]]);
    }
        long int value = pis.computeObjectiveFunction(newsol);
        return new emili::pfsp::PermutationFlowShopSolution(value,newsol);
    }
    return currentSolution->clone();
}


void emili::pfsp::PfspExchangeNeighborhood::reset()
{    
    start_position = 1;
    end_position = 1;
}

emili::Solution* emili::pfsp::PfspTransposeNeighborhood::computeStep(emili::Solution* value)
{

    emili::iteration_increment();
    if(sp_iterations >= njobs)
    {
        return nullptr;
    }
    else
    {
        sp_iterations++;
        start_position = (start_position%njobs)+1;
        std::vector < int >& newsol = ((emili::pfsp::PermutationFlowShopSolution*)value)->getJobSchedule();
        int endpos = start_position<njobs?start_position+1:1;     
        std::swap(newsol[start_position],newsol[endpos]);        
        long int new_value = pis.computeObjectiveFunction(newsol);        
        value->setSolutionValue(new_value);
        return value;
    }
}

void emili::pfsp::PfspTransposeNeighborhood::reverseLastMove(Solution *step)
{
    std::vector < int >& newsol = ((emili::pfsp::PermutationFlowShopSolution*)step)->getJobSchedule();
    int endpos = start_position<njobs?start_position+1:1;
    std::swap(newsol[start_position],newsol[endpos]);
}

int emili::pfsp::PfspTransposeNeighborhood::size()
{
    return pis.getNjobs();
}

emili::Solution* emili::pfsp::XTransposeNeighborhood::computeStep(emili::Solution* value)
{

    emili::iteration_increment();
    if(sp_iterations >= njobs)
    {
        return nullptr;
    }
    else
    {
        sp_iterations++;
        start_position = (start_position%njobs)+1;
        //start_position = njobs-sp_iterations+1;
        emili::pfsp::PermutationFlowShopSolution* p = (emili::pfsp::PermutationFlowShopSolution*) value;
        std::vector < int >& newsol = ((emili::pfsp::PermutationFlowShopSolution*)value)->getJobSchedule();
        int endpos =  start_position<njobs?start_position+1:1;
        std::swap(newsol[start_position],newsol[endpos]);

        long int new_value = 0;
        //clock_t start = clock();

            //clock_t start = clock();


           //new_value =  pis.computeObjectiveFunction(newsol,etm,start_position,endpos);
        value->setSolutionValue(new_value);
         //std::cout <<  (clock()-start) << ", " ;//<< std::endl;
            //assert(new_value == pis.computeObjectiveFunction(newsol));

        return value;
    }
}


emili::Solution* emili::pfsp::PfspTransposeNeighborhood::random(Solution *currentSolution)
{
    int njobs = pis.getNjobs()-1;
    std::vector < int > newsol(((emili::pfsp::PermutationFlowShopSolution*)currentSolution)->getJobSchedule());
    int best_i = (emili::generateRandomNumber()%njobs)+1;
    std::swap(newsol[best_i],newsol[best_i+1]);
    long int value = pis.computeObjectiveFunction(newsol);
    return new emili::pfsp::PermutationFlowShopSolution(value,newsol);
}

emili::Solution* emili::pfsp::PfspTransposeNeighborhood::random(Solution *currentSolution,int size)
{
    if(size > 0)
    {
        std::vector < int > newsol(((emili::pfsp::PermutationFlowShopSolution*)currentSolution)->getJobSchedule());
        int njobs = pis.getNjobs()-1;
        std::vector < int > indexes;
        std::vector < int > kn(njobs+1,0);
        if(size >= njobs)
        {
            size = njobs-2;
        }

        int k=0;
        while(k < size)
        {
            int ind = (emili::generateRandomNumber()%njobs)+1;
            if(kn[ind] == 0)
            {
                kn[ind] = 1;
                indexes.push_back(ind);
                k++;
            }
        }
        int i = 0;
        for(; i < size; i++)
        {
            std::swap(newsol[indexes[i]],newsol[indexes[i]+1]);
        }
        long int value = pis.computeObjectiveFunction(newsol);
        return new emili::pfsp::PermutationFlowShopSolution(value,newsol);
    }
    return currentSolution->clone();
}

void emili::pfsp::PfspTransposeNeighborhood::reset()
{
    sp_iterations = 1;
    start_position = 0;
}




void emili::pfsp::PfspInsertNeighborhood::reset()
{
    start_position = 1;
    end_position = 1;

}

bool emili::pfsp::PfspTerminationClassic::terminate(Solution* currentSolution,Solution* newSolution)
{
   //std::cout << currentSolution.getSolutionValue() << " <= " << newSolution.getSolutionValue()<<std::endl;
    if(newSolution == nullptr)
    {
        return true;
    }
    else
    {
        return currentSolution->operator <=(*newSolution);
    }
}

emili::Solution* emili::pfsp::PfspRandomSwapPertub::perturb(emili::Solution* solution)
{

    std::vector < int > perturbed(((emili::pfsp::PermutationFlowShopSolution*)solution)->getJobSchedule());
    int n = pfs.getNjobs()-1;
    int pos1 = emili::generateRandomNumber()%n +1;
    int pos2 = emili::generateRandomNumber()%n +1;
    int pos3 = emili::generateRandomNumber()%n +1;
    int pos4 = emili::generateRandomNumber()%n +1;
    int swap = perturbed[pos2];
    perturbed[pos2] = perturbed[pos1];
    perturbed[pos1] = swap;
    swap = perturbed[pos4];
    perturbed[pos4] = perturbed[pos3];
    perturbed[pos3] = swap;
    emili::pfsp::PermutationFlowShopSolution* nsol = new emili::pfsp::PermutationFlowShopSolution(perturbed);
    pfs.evaluateSolution(*nsol);
    return nsol;
}

emili::Solution* emili::pfsp::PfspTestAcceptance::accept(Solution *candidate1, Solution *candidate2)
{
    int chance = emili::generateRandomNumber()%100;
    emili::Solution* c = candidate1;
    emili::Solution* c2 = candidate2;
    if(candidate1->operator >(*candidate2)){
        c = candidate2;
        c2 = candidate1;
    }

    if(chance <= percentage)
    {
        return c;
    }
    else
    {
        return c2;
    }

}

emili::Solution* emili::pfsp::SOAacceptance::accept(Solution *intensification_solution, Solution *diversification_solution)
{
    float intens = intensification_solution->getSolutionValue();
    float divers = diversification_solution->getSolutionValue();
    if(diversification_solution->operator >(*intensification_solution))
    {
        float prob = std::exp(100.0f*((intens-divers)/intens)/temperature);
        if(prob < 1.0 && emili::generateRealRandomNumber()>prob)
        {
            return intensification_solution;
        }
    }
    return diversification_solution;
}

bool emili::pfsp::SOAtermination::terminate(Solution *currentSolution, Solution *newSolution)
{
    if(currentStep < numberOfSteps){
        currentStep++;
        return false;
    }
    else
    {
        return true;
    }
}

void emili::pfsp::SOAtermination::reset()
{
    currentStep=0;
}

bool emili::pfsp::PfspTerminationIterations::terminate(Solution* currentSolution, Solution* newSolution)
{

    if(iterations < maxIterations)
    {
        if(newSolution != nullptr && currentSolution->operator <=(*newSolution))
        {
            iterations++;
        }

            return false;
    }
    return true;
}

void emili::pfsp::PfspTerminationIterations::reset()
{
    iterations=0;
}

bool emili::pfsp::KarTermination::terminate(Solution *currentSolution, Solution *newSolution)
{
    if(iterations < maxIterations)
    {
        if(*newSolution < *currentSolution)
        {
            //std::cout << "iterations " << iterations << "\n news <- " << newSolution->getSolutionValue() << "\n current <- " <<currentSolution->getSolutionValue() << "\n";
            iterations = 0;
        }
        else
        {
           iterations++;
        }
         return false;
    }
    return true;
}

void emili::pfsp::KarTermination::reset()
{
    iterations=0;
}



size_t emili::pfsp::PfspTabuHashMemory::calc_hash(std::vector< int >& v )
{
    std::hash< std::string> hasher;
    std::ostringstream oss;
    for(std::vector< int > ::iterator iter = v.begin();iter!= v.end();++iter )
    {
        oss << *iter;
    }
    return hasher(oss.str());
}

bool emili::pfsp::PfspTabuHashMemory::tabu_check(emili::Solution* toCheck)
{  
    std::vector< int > & v = ((emili::pfsp::PermutationFlowShopSolution*)toCheck)->getJobSchedule();
    size_t c_hash = calc_hash(v);
    for(std::vector< size_t>::iterator iter = tabuVector.begin();iter!=tabuVector.end();++iter)
    {
        if(c_hash == *iter){
            return false;
        }
    }
    return true;
}

void emili::pfsp::PfspTabuHashMemory::forbid(Solution *solution)
{
    std::vector < int >& v = ((emili::pfsp::PermutationFlowShopSolution*)solution)->getJobSchedule();
    if(tabu_check(solution))
    {
        if(tt_index < this->tabutenure){
            tabuVector.push_back(calc_hash(v));
            tt_index++;
        }
        else
        {
            tabuVector.erase(tabuVector.begin());
            tabuVector.push_back(calc_hash(v));
        }
    }
}

std::string emili::pfsp::PermutationFlowShopSolution::getSolutionRepresentation()
{
    std::ostringstream oss;    
    int size = solution.size();
    oss << "{ ";

    for (int i = 1; i < size-1; ++i)
      oss << solution[i] << ", " ;

    oss << solution[size-1];
    oss << " }";
    return oss.str();
}

void emili::pfsp::PfspTabuHashMemory::reset()
{
    tt_index = 0;
    tabuVector.clear();
}


bool emili::pfsp::PfspTabuValueMemory::tabu_check(emili::Solution* toCheck)
{

    double value = toCheck->getSolutionValue();
    for(std::vector< double>::iterator iter = tabuVector.begin();iter!=tabuVector.end();++iter)
    {
        if( value == *iter){
            return false;
        }
    }
    return true;
}

void emili::pfsp::PfspTabuValueMemory::forbid(Solution *solution)
{

    if(tabu_check(solution))
    {
        if(tt_index < this->tabutenure){
            tabuVector.push_back(solution->getSolutionValue());
            tt_index++;
        }
        else
        {
            tabuVector.erase(tabuVector.begin());
            tabuVector.push_back(solution->getSolutionValue());
        }
    }
}

void emili::pfsp::PfspTabuValueMemory::reset()
{
    tt_index = 0;
    tabuVector.clear();
}

bool emili::pfsp::PfspFullSolutionMemory::tabu_check(emili::Solution* toCheck)
{

    std::vector< int >& value = ((emili::pfsp::PermutationFlowShopSolution*)toCheck)->getJobSchedule();

    for(std::vector< std::vector<int > >::iterator iter = tabuVector.begin();iter!=tabuVector.end();++iter)
    {
        std::vector< int > t = *iter ;
        bool ret = false;
        for(int i=0; i< t.size() ; i++)
        {
            if(t[i]!=value[i])
            {
                ret = true;
                break;
            }
        }
        if(!ret)
        {
            return ret;
        }
    }
    return true;
}

void emili::pfsp::PfspFullSolutionMemory::forbid(Solution *solution)
{

    if(tabu_check(solution))
    {
        std::vector< int >& solution_vector = ((emili::pfsp::PermutationFlowShopSolution*)solution)->getJobSchedule();
        if(tt_index < this->tabutenure){
            tabuVector.push_back(solution_vector);
            tt_index++;
        }
        else
        {
            tabuVector.erase(tabuVector.begin());            
            tabuVector.push_back(solution_vector);

        }
    }
}

void emili::pfsp::PfspFullSolutionMemory::reset()
{
    tt_index = 0;
    tabuVector.clear();
}

bool emili::pfsp::PfspMovesMemory::tabu_check(emili::Solution* toCheck)
{

    return tabu_check(lastMove);
}

bool emili::pfsp::PfspMovesMemory::tabu_check(std::pair< int,int > value)
{
    for(std::vector< std::pair<int,int > >::iterator iter = tabuVector.begin();iter!=tabuVector.end();++iter)
    {
        std::pair< int ,int> t = *iter ;
        if(value==t)
        {
            return false;
        }
    }
    return true;
}

void emili::pfsp::PfspMovesMemory::forbid(Solution *solution)
{
    //std::cout << lastMove.first << " <-> " << lastMove.second << std::endl;
    if(tabu_check(lastMove))
    {
        if(tt_index < this->tabutenure){
            tabuVector.push_back(lastMove);
            tt_index++;
        }
        else
        {
            tabuVector.erase(tabuVector.begin());
            tabuVector.push_back(lastMove);
        }
    }
}

void emili::pfsp::PfspMovesMemory::registerMove(emili::Solution* base,emili::Solution* solution)
{
    lastMove = neigh->lastMove();
}

void emili::pfsp::PfspMovesMemory2::registerMove(emili::Solution* base,emili::Solution* solution)
{
    std::pair<int,int> lm = neigh->lastMove();

    int job = ((emili::pfsp::PermutationFlowShopSolution*) base)->getJobSchedule()[lm.second];
    lastMove = std::pair<int,int>(job,lm.first);
}

void emili::pfsp::PfspMovesMemory::reset()
{
    tt_index = 0;
    tabuVector.clear();
}

bool emili::pfsp::TSABtestMemory::tabu_check(std::pair< int,int > value)
{
    int k = value.first;
    int l = value.second;
    if(k<l)
    {
        for(std::vector< std::pair<int,int > >::iterator iter = tabuVector.begin();iter!=tabuVector.end();++iter)
        {
            std::pair< int ,int> t = *iter ;
            if(t.first==k && (t.second > k && t.second <= l) )
            {
                return false;
            }
        }
    }
    else
    {
        for(std::vector< std::pair<int,int > >::iterator iter = tabuVector.begin();iter!=tabuVector.end();++iter)
        {
            std::pair< int ,int> t = *iter ;
            if(t.second==k && (t.first >= l && t.first < k) )
            {
                return false;
            }
        }
    }
    return true;
}

void emili::pfsp::TSABtestMemory::forbid(Solution *solution)
{
    std::pair< int , int > fmove;
    int k = lastMove.first;
    int l = lastMove.second;    
    if(k<l)
    {
        fmove = std::pair< int, int > (k,k+1);
    }
    else
    {
        fmove = std::pair< int, int > (k-1,k);
    }

    if(tabu_check(lastMove))
    {
        if(tt_index < this->tabutenure){

            tabuVector.push_back(fmove);
            tt_index++;
        }
        else
        {
            tabuVector.erase(tabuVector.begin());
            tabuVector.push_back(fmove);
        }
    }
}

bool emili::pfsp::TSABMemory::tabu_check(emili::Solution* toCheck)
{
    std::vector< int > & pi = ((emili::pfsp::PermutationFlowShopSolution*)toCheck)->getJobSchedule();
    return tabu_check(lastMove,pi);
}

bool emili::pfsp::TSABMemory::tabu_check(std::pair< int,int > value,std::vector< int>& pi)
{
    int k = lastMove.first;
    int l = lastMove.second;

    int size = k>l?k-l:l-k;

    for(int i=0;i<size;i++)
    {
        std::pair< int , int > toTest;
        if(k<l)
        {
            toTest = std::pair<int , int >(pi[k+1+i],pi[k]);
        }
        else
        {
            toTest = std::pair<int, int >( pi[k] , pi[l+i] );
        }

        for(std::vector< std::pair<int,int > >::iterator iter = tabuVector.begin();iter!=tabuVector.end();++iter)
        {
            std::pair< int ,int> t = *iter ;
            if(toTest==t)
            {               
                return false;
            }
        }
    }

    return true;
}

void emili::pfsp::TSABMemory::forbid(Solution *solution)
{
    std::pair< int , int > fmove;
    int k = lastMove.first;
    int l = lastMove.second;
    std::vector< int > & pi =((emili::pfsp::PermutationFlowShopSolution*)solution)->getJobSchedule();
    if(k<l)
    {
        fmove = std::pair< int, int > (pi[l],pi[k]);
    }
    else
    {
        fmove = std::pair< int, int > (pi[k],pi[l]);
    }

    if(tabu_check(lastMove,pi))
    {
        if(tt_index < this->tabutenure){

            tabuVector.push_back(fmove);
            tt_index++;
        }
        else
        {
            tabuVector.erase(tabuVector.begin());
            tabuVector.push_back(fmove);
        }
    }
}



emili::Solution* emili::pfsp::VNDBestSearch::search(emili::Solution* initial)
{

   emili::Solution* temp = bs1.search(initial);

   temp = bs2.search(temp);

   temp = bs3.search(temp);
   return temp;
}

emili::Solution* emili::pfsp::VNDFirstSearch::search(emili::Solution* initial)
{

   emili::Solution* temp = bs1.search(initial);

   temp = bs2.search(temp);

   temp = bs3.search(temp);

   return temp;
}

emili::pfsp::GVNS_innerloop::GVNS_innerloop(InitialSolution& initialSolutionGenerator)
{
    this->init = &initialSolutionGenerator;
    this->termcriterion = new emili::LocalMinimaTermination();
    this->rneigh = new emili::pfsp::GVNS_RIS_Neighborhood((emili::pfsp::PermutationFlowShop&)initialSolutionGenerator.getProblem());
    this->neighbh = this->rneigh;
}
//Users/federicopagnozzi/Desktop/phd/PFSP_NO_IDLE/PFSP_NOIDLE/PFSP_NoIdle/I_7_500_50_4.txt NIPFSP_MS gvns nwslack soaper 8 rndmv insert 5 rndmv insert 1 rndmv transpose 1 -it 30 rnds 30
emili::Solution* emili::pfsp::GVNS_innerloop::search(emili::Solution *initial)
{
    termcriterion->reset();
    neighbh->reset();

    bestSoFar = init->generateEmptySolution();
    emili::Solution* incumbent = bestSoFar;

    *incumbent = *initial;
    //bestSoFar->setSolutionValue(bestSoFar->getSolutionValue()+1);

    do{
        if(bestSoFar != incumbent)
        {
            delete bestSoFar;
            bestSoFar = incumbent;
        }
        rneigh->setReference(bestSoFar);
        for(Neighborhood::NeighborhoodIterator iter = neighbh->begin(incumbent);iter!=neighbh->end();++iter)
        {
            emili::Solution* ithSolution = *iter;
            if(incumbent->operator >(*ithSolution)){
                if(incumbent!=bestSoFar)
                delete incumbent;

                incumbent = ithSolution;
                break;
            }
            else
            {
                delete ithSolution;
            }

        }

    }while(!termcriterion->terminate(bestSoFar,incumbent));
    return bestSoFar;
}

void emili::pfsp::GVNS_RIS_Neighborhood::reset()
{
    index = 0;
}

emili::Solution* emili::pfsp::GVNS_RIS_Neighborhood::random(emili::Solution* currentSolution)
{
    return currentSolution;
}

emili::Solution* emili::pfsp::GVNS_RIS_Neighborhood::random(emili::Solution* currentSolution,int size)
{
    return currentSolution;
}

emili::Neighborhood::NeighborhoodIterator emili::pfsp::GVNS_RIS_Neighborhood::begin(Solution *base)
{
    return emili::Neighborhood::NeighborhoodIterator(this,base);
}

emili::Solution* emili::pfsp::GVNS_RIS_Neighborhood::computeStep(Solution *step)
{
    if(index > njobs)
    {
        return nullptr;
    }
    else
    {
        std::vector< int >& bsf = ((emili::pfsp::PermutationFlowShopSolution*)reference)->getJobSchedule();
        std::vector< int >& base = ((emili::pfsp::PermutationFlowShopSolution*)step)->getJobSchedule();

        int k = bsf[index];
        int k_index = 1;

        for(int i=1; i<= njobs; ++i)
        {
            if(base[i] == k)
            {
                k_index = i;
                break;
            }
        }

        std::vector< int > bestCombination(base);
        std::vector< int > workCombination(base);
        long int best_res = step->getSolutionValue();

        for(int i=1; i<=njobs;++i )
        {
            if(i!=k_index)
            {
                workCombination.erase(workCombination.begin()+k_index);
                workCombination.insert(workCombination.begin()+i,k);
                long int nvalue = pis.computeObjectiveFunction(workCombination);
                if(nvalue < best_res)
                {
                    best_res = nvalue;
                    bestCombination = workCombination;
                }
                workCombination = base;
            }
        }
        index++;

           //std::cout<< best_res << std::endl;

        return new emili::pfsp::PermutationFlowShopSolution(best_res,bestCombination);
    }
}

int emili::pfsp::CompoundPerturbation::calc_distance(std::vector<int> &x, std::vector<int> &y)
{
            int arrX[nbj];
            int arrY[nbj];
            int E[nbj][nbj];
            int F[nbj][nbj];


            int i,j;
            for(i=0; i<nbj; i++)
            {

                arrX[x[i+1]] = i;
                arrY[y[i+1]] = i;
            }

            for(i=0; i<nbj; i++)
                for(j=0; j<nbj; j++)
                {
                    if( arrX[i] > arrX[j] )
                        E[i][j] = 1;
                    else
                        E[i][j] = 0;

                    if( arrY[i] < arrY[j] )
                        F[i][j] = 1;
                    else
                        F[i][j] = 0;
                }

            int count = 0;
            for(i=0; i<nbj; i++)
                for(j=0; j<nbj; j++)
                    if( E[i][j]==1&&F[i][j]==1 )
                        count++;
            return count;
}

emili::Solution* emili::pfsp::CompoundPerturbation::perturb(emili::Solution *solution)
{

    int i=0;
    std::vector < emili::pfsp::PermutationFlowShopSolution* > phy(omega);
    std::vector< int > distance_vector(omega);
    emili::Solution* best_sf = getAlgo()->getBestSoFar();
    std::vector< int >& sol_schedule = ((emili::pfsp::PermutationFlowShopSolution*)best_sf)->getJobSchedule();
    emili::pfsp::PermutationFlowShopSolution* best_per=nullptr;
    int min_value = std::numeric_limits<int>::max();
    do
    {
        emili::pfsp::PermutationFlowShopSolution* candidate = (emili::pfsp::PermutationFlowShopSolution*)solution->clone();
        for(int j=0;j<d;j++)
        {
            emili::pfsp::PermutationFlowShopSolution* cand = candidate;
            float p = emili::generateRealRandomNumber();
            if(p <= pc)
            {
                candidate = (emili::pfsp::PermutationFlowShopSolution*)ins.random(cand);
                delete cand;
            }
            else
            {
                candidate = (emili::pfsp::PermutationFlowShopSolution*)tra.random(cand);
                delete cand;
            }
        }
        int dd = calc_distance(candidate->getJobSchedule(),sol_schedule);
        if(dd > 0)
        {
            int val = candidate->getSolutionValue();
            if(min_value > val)
            {
                best_per = candidate;
                min_value = val;
            }

            phy[i] = candidate;
            distance_vector[i] = dd;
            i++;
        }
        else
        {
            delete candidate;
        }
    }while(i < omega);

    emili::pfsp::PermutationFlowShopSolution* toRet = nullptr;
    int min_d = nbj*nbj;
    int min_d_index=0;
    if(min_value < best_sf->getSolutionValue())
    {
        toRet = best_per;
    }

    if(toRet==nullptr)
    {
        for(i=0;i<omega;i++)
        {
            if(distance_vector[i]<min_d)
            {
                min_d = distance_vector[i];
                min_d_index = i;
            }
        }
        toRet = phy[min_d_index];
    }

    return toRet->clone();

}


int emili::pfsp::RIS::neh_ig(std::vector<int>& solPartial, int x)
{
   int k = solPartial[x];
   solPartial.erase(solPartial.begin()+x);
   int tmp=0,ind=1;
   int min = std::numeric_limits<int>::max();
   int sops = njob+1;
   std::vector< int > solTMP(sops,0);
   for(int r=1; r<sops; r++){
    emili::iteration_increment();
       for(int h=1; h<r; h++)
           solTMP[h]=solPartial[h];
       solTMP[r]=k;
       for(int h=r+1; h<=sops; h++)
           solTMP[h]=solPartial[h-1];

       tmp = instance.computeObjectiveFunction(solTMP);

       if(tmp<min){
           min=tmp;
           ind=r;
       }

   }
   ni_position = ind;
   solPartial.insert(solPartial.begin()+ind,k);
   return min;
}

void emili::pfsp::RIS::invertPerturbation(std::vector<int>& pi, std::vector<int>& pi_i)
{
    for(int i = 1; i<=njob;i++)
    {
        pi_i[pi[i]] = i;
    }
}

emili::Solution* emili::pfsp::RIS::search(Solution *initial)
{
    *bestSoFar = *initial;
    int h = 1;
    int i = 1;
    std::vector<int> pi = ((emili::pfsp::PermutationFlowShopSolution*)initial)->getJobSchedule();
    emili::LocalSearch* ls = emili::getAlgo();
    emili::Solution* best = bestSoFar;
    if(ls!=nullptr)
    {
        best = ls->getBestSoFar();
    }
    std::vector<int> pi_r = ((emili::pfsp::PermutationFlowShopSolution*)best)->getJobSchedule();
    int pi_value = initial->getSolutionValue();
    int initial_value = pi_value;
    std::vector<int> pi_i(njob+1,0);
    invertPerturbation(pi,pi_i);
    while(i < njob+1)
    {
        h = h%njob+1;
        int job = pi_r[h];
        int k = pi_i[job];
        int value = neh_ig(pi,k);
        if(value < pi_value)
        {
            pi_value = value;
            i = 1;
            invertPerturbation(pi,pi_i);
        }
        else
        {
            h++;
            i++;
            pi.erase(pi.begin()+ni_position);
            pi.insert(pi.begin()+k,job);
        }

    }
    if(initial_value==pi_value)
    {
        return bestSoFar->clone();
    }
    else
    {
        emili::Solution* last = new PermutationFlowShopSolution(pi_value,pi);
        *bestSoFar = *last;
        return last;
    }
}

int emili::pfsp::NoIdle_RIS::neh_ig(std::vector<int>& solPartial, int x)
{
   int k = solPartial[x];
   solPartial.erase(solPartial.begin()+x);
   int tmp=0,ind=1;
   int min = std::numeric_limits<int>::max();
   int sops = njob+1;
   std::vector< int > solTMP(sops,0);
   pis.computeNoIdleTAmatrices(solPartial,head,tail);
   for(int r=1; r<sops; r++){
    emili::iteration_increment();
       for(int h=1; h<r; h++)
           solTMP[h]=solPartial[h];
       solTMP[r]=k;
       for(int h=r+1; h<=sops; h++)
           solTMP[h]=solPartial[h-1];

       //tmp = instance.computeObjectiveFunction(solTMP,njob);

       long int c_cur = head[1][r-1]+pmatrix[k][1];
       long int c_max = c_cur+tail[1][r];
       long int a = 0;
       long int aa = 0;
       for (int i = 2; i <= nmac; ++i) {
           int c_pm = head[i][r-1];
           if(c_pm + aa < c_cur)
           {
               aa += c_cur - (c_pm+aa);
               c_cur = c_cur + pmatrix[k][i];
           }
           else
           {
               c_cur = c_pm + aa + pmatrix[k][i];
           }
           long int c_can = (c_cur+a+tail[i][r]);
           c_max = c_max>c_can?c_max:c_can;
           a = a+c_max-c_can;
       }

       //assert(c_max == tmp);
       tmp = c_max;
       if(tmp<min){
           min=tmp;
           ind=r;
       }

   }
   ni_position = ind;
   solPartial.insert(solPartial.begin()+ind,k);
   return min;
}

emili::Solution* emili::pfsp::NoIdleIGper::perturb(Solution *solution)
{
    int index;
    int min;
    int k,tmp=0,ind=1;
    int nmac = instance.getNmachines();
    PfspInstance& pis = instance.getInstance();
    std::vector< int > removed;
    std::vector< int > solPartial(((emili::pfsp::PermutationFlowShopSolution*)solution)->getJobSchedule());
    //std::cout << "partial size " << solPartial.size() << std::endl;
    int size = solPartial.size();
    std::vector< int > solTMP(size,0);
    int sops = size-1;
    for(int k = 0; k < d; k++) {
        index = (emili::generateRandomNumber()%sops)+1;
        //std::cout << index << " " ;//<< std::endl;
        removed.push_back(solPartial[index]);
        solPartial.erase(solPartial.begin() + index);
        sops--;
    }

    for(int l=0;l<removed.size();l++){
        sops++;
        k=removed[l];
        min = std::numeric_limits<int>::max();
        pis.computeNoIdleTAmatrices(solPartial,head,tail,sops);
        for(int r=1; r<sops; r++){

            for(int h=1; h<r; h++)
                solTMP[h]=solPartial[h];
            solTMP[r]=k;
            for(int h=r+1; h<=sops; h++)
                solTMP[h]=solPartial[h-1];

            long int c_cur = head[1][r-1]+pmatrix[k][1];
            long int c_max = c_cur+tail[1][r];
            long int a = 0;
            long int aa = 0;
            for (int i = 2; i <= nmac; ++i) {
                int c_pm = head[i][r-1];
                if(c_pm + aa < c_cur)
                {
                    aa += c_cur - (c_pm+aa);
                    c_cur = c_cur + pmatrix[k][i];
                }
                else
                {
                    c_cur = c_pm + aa + pmatrix[k][i];
                }
                long int c_can = (c_cur+a+tail[i][r]);

                c_max = c_max>c_can?c_max:c_can;
                a = a+c_max-c_can;
            }
            //tmp = instance.computeObjectiveFunction(solTMP,sops);

            //assert(c_max == tmp);
            tmp = c_max;
            if(tmp<min){
                min=tmp;
                ind=r;
            }

        }
        solPartial.insert(solPartial.begin()+ind,k);
        //std::cout << "end insert " << solPartial.size() << std::endl;
    }

    //assert(min == instance.computeObjectiveFunction(solPartial));
    emili::pfsp::PermutationFlowShopSolution* s = new emili::pfsp::PermutationFlowShopSolution(min,solPartial);
    //instance.evaluateSolution(*s);
    return s;
}

int emili::pfsp::NoWait_RIS::neh_ig(std::vector<int>& solPartial, int x)
{
   /* the jobs around the one to insert
      [1,2,....,jmo,j,jpo,...njob]*/

   int j = solPartial[x];
   int jpo = solPartial[(x+1)%(njob+1)];
   int jmo = solPartial[x-1];
   // The value of the makespan without j
   long int base_value = pis.computeNWMS(solPartial)
                        + distance[jmo][jpo]
                        - distance[jmo][j]
                        - distance[j][jpo];

   solPartial.erase(solPartial.begin()+x);
   int tmp=0,ind=1;
   int min = std::numeric_limits<int>::max();
   int sops = njob+1;
   std::vector< int > solTMP(sops,0);

   for(int r=1; r<sops; r++){
    emili::iteration_increment();
       for(int h=1; h<r; h++)
           solTMP[h]=solPartial[h];
       solTMP[r]=j;
       for(int h=r+1; h<=sops; h++)
           solTMP[h]=solPartial[h-1];

      // tmp = instance.computeObjectiveFunction(solTMP,njob);
       /*Here I take the job before and after j*/
       int k = solPartial[r%njob];
       int kmo = solPartial[r-1];
       int delta = distance[kmo][j] + distance[j][k] - distance[kmo][k] ;
       int c_max = base_value + delta;

       //assert(c_max == tmp);
       tmp = c_max;
       if(tmp<min){
           min=tmp;
           ind=r;
       }

   }
   ni_position = ind;
   solPartial.insert(solPartial.begin()+ind,j);
   return min;
}

emili::Solution* emili::pfsp::RandomNoWait_RIS::search(Solution *initial)
{
    *bestSoFar = *initial;
    Solution* random = nullptr;
    int h = 1;
    int i = 1;
    std::vector<int> pi = ((emili::pfsp::PermutationFlowShopSolution*)initial)->getJobSchedule();
    emili::LocalSearch* ls = emili::getAlgo();
    emili::Solution* best = bestSoFar;
    if(ls!=nullptr)
    {
        best = ls->getBestSoFar();
    }
    std::vector<int> pi_r = ((emili::pfsp::PermutationFlowShopSolution*)best)->getJobSchedule();
    if(cmax == initial->getSolutionValue())
    {
        random = rand.generateSolution();
        pi_r = ((emili::pfsp::PermutationFlowShopSolution*)random)->getJobSchedule();
    }

    int pi_value = initial->getSolutionValue();
    int initial_value = pi_value;
    std::vector<int> pi_i(njob+1,0);
    invertPerturbation(pi,pi_i);
    while(i < njob+1)
    {
        h = h%njob+1;
        int job = pi_r[h];
        int k = pi_i[job];
        int value = neh_ig(pi,k);
        if(value < pi_value)
        {
            pi_value = value;
            i = 1;
            invertPerturbation(pi,pi_i);
        }
        else
        {
            h++;
            i++;
            pi.erase(pi.begin()+ni_position);
            pi.insert(pi.begin()+k,job);
        }

    }
    if(random !=nullptr)
    {
        delete random;
    }
    if(initial_value==pi_value)
    {
        cmax = initial_value;
        return bestSoFar->clone();
    }
    else
    {
        emili::Solution* last = new PermutationFlowShopSolution(pi_value,pi);
        cmax = pi_value;
        *bestSoFar = *last;
        return last;
    }
}

void emili::pfsp::BeamSearchHeuristic::buildXi()
{
    double wprod = (njobs-2.0)/4.0;
    std::vector<double> wj(njobs+1,0.0);
    for(int j=1; j<=njobs ; j++ )
    {

        double pij = pi[j][1];
        for(int i = 2; i<=nmacs ; i++)
        {
            pij = pij + pi[j][i];
            double wjprime = 0;
            for(int ip = 1; ip < i ; ip++ )
            {
                wjprime += pi[j][ip];
            }
            wjprime = (wjprime*nmacs)/(i-1);
            wj[j] = wj[j] + wjprime;
        }

        wj[j] = wj[j]*wprod;
        xi[j] = pij+wj[j];
        xi_order.push_back(j);
    }
    std::vector<double>& xil = xi;

    std::sort(xi_order.begin(),xi_order.end(),[xil,wj](int i1,int i2){
                                                            if(xil[i1]==xil[i2] && i2!=0)
                                                                return wj[i1]<wj[i2];
                                                            else
                                                               return xil[i1] < xil[i2];
    });
/*
    for(int i = 0; i<=njobs; i++)
        std::cout << " " << xi[i];
    std::cout << "\n";

    for(int i = 0; i<=njobs; i++)
        std::cout << " " << xi_order[i];
    std::cout << "\n";

    for(int i=1; i<=njobs; i++)
    {
        std::cout << i << ": ";
        for(int k = 1; k<=nmacs; k++)
             std::cout << pi[i][k] << " ";

        std::cout << std::endl;
    }
*/
}
std::vector<double>& emili::pfsp::BeamSearchHeuristic::getXi()
{
    return xi;
}

std::vector<int>& emili::pfsp::BeamSearchHeuristic::getXi_order()
{
    return xi_order;
}

emili::Solution* emili::pfsp::BeamSearchHeuristic::generate()
{


    //start.buildChildren();

    int k = 1;

    std::vector<bs_node*> current_iteration;
    int g = _gamma<njobs?_gamma:njobs;
   // for(int x=0;x<g;x++)
    //{
       bs_node* start = new bs_node(*this,0);
      current_iteration.push_back(start);
    //}

    //std::sort(current_iteration.begin(),current_iteration.end(),[](bs_node* i1,bs_node* i2){return *i1<*i2;});

    do{
       // std::cout << "#child " << current_iteration.size() << std::endl;
        //current_iteration[0]->buildChildren();
        //std::vector<bs_node*> next_iteration = current_iteration[0]->getChildren();
        std::vector<bs_node*> next_iteration;
        int size = current_iteration.size();

        size = _gamma<size?_gamma:size;

        //for(int i=1;i < size; i++)
        //std::cout << k << " " << k << std::endl;
        for(int i=0;i < size; i++)
        {
           // std::cout << " ITERATION " << i <<"/" << size << std::endl;
         //   std::cout << "Node value " << current_iteration[i]->g_value << std::endl;
            //current_iteration[i]->buildChildren();
            //std::vector<bs_node*> ccc = current_iteration[i]->getChildren();
            //next_iteration.insert(next_iteration.end(),ccc.begin(),ccc.end());
            int csize = njobs-k;
            double W = 0.0;
            for(int c=0;c < csize; c++)
            {
                bs_node* b = new bs_node(*current_iteration[i],c);
                W += b->TT;
                next_iteration.push_back(b);
            }

            for(int c=0;c < csize; c++)
            {
                int node = i*csize+c;
                next_iteration[node]->calcG(W);
            }

        }

  //      std::vector< bs_node* >::iterator iter = next_iteration.begin();
    //    for(;iter!=next_iteration.end();++iter)
      //      (*iter)->calcG(W);
        std::vector< bs_node* >::iterator iter = current_iteration.begin();
            for(;iter!=current_iteration.end();++iter)
            {
                delete *iter;
            }
         current_iteration = next_iteration;
         std::sort(current_iteration.begin(),current_iteration.end(),[](bs_node* i1,bs_node* i2){return *i1<*i2;});

       /*  std::cout << k << ": ";
         for(;iter!=current_iteration.end();++iter)
         {
             //std::cout << "K " << (*iter)->k << " J " << (*iter)->kjob << " TE " << (*iter)->TE << " TT " << (*iter)->TT << " TI " << (*iter)->TI << " G " <<  (*iter)->g_value << std::endl;
             std::cout << "( " << (*iter)->kjob << "," << (*iter)->g_value << ") ";
         }
         std::cout << std::endl;*/
         k++;
    }while(k<(njobs));
    /*std::vector< bs_node* >::iterator iter = current_iteration.begin();
      for(;iter!=current_iteration.end();++iter)
    {
    std::vector<int> seq = (*iter)->getPermutation();
    double res = pis.computeObjectiveFunction(seq);
    std::cout << "K " << (*iter)->k << " G "<< (*iter)->g_value << " V " << res << std::endl;
     }*/
    std::vector<int> seq = current_iteration[0]->getPermutation();
    double res = current_iteration[0]->g_value;
    std::vector< bs_node* >::iterator iter = current_iteration.begin();
        for(;iter!=current_iteration.end();++iter)
        {
            delete *iter;
        }
    //pis.computeObjectiveFunction(seq);
//    std::cout << "sol " << seq.size() << " children " <<  std::endl;

    return new emili::pfsp::PermutationFlowShopSolution(res,seq);
}

/*
void emili::pfsp::BeamSearchHeuristic::bs_node::evaluateNode()
{

    completionTimes[1] = father->completionTimes[1]+init.pi[kjob][1];
    double I = 0.0;
    //tpj[1] = tpj[1] - init.pi[kjob][1];
    double W = father->completionTimes[1]+tpj[1];
    for(int i=2;i<=m;i++)
    {
        completionTimes[i] = std::max(father->completionTimes[i],completionTimes[i-1])+init.pi[kjob][i];
      //  tpj[i] = tpj[i] - init.pi[kjob][i];
        W = std::max((double)father->completionTimes[i],W)+tpj[i];
        double iden = i-1+((k-1)*(m-i+1)/(n-2));
        I = I + (m * std::max(completionTimes[i-1]-father->completionTimes[i],0))/iden;
    }
    //tpd = tpd - init.dueDates[kjob];

    //W =  W - tpd;
    //W = std::max(W,0.0);

    double E = std::max((int)init.dueDates[kjob]-completionTimes[m],0);
    double T = std::max(completionTimes[m]-(int)init.dueDates[kjob],0);
    TE = TE + E;
    TT = TT + T;
    TI = TI + I;

    //std::cout << "K " << k << " TE " << TE << " TT " << TT << " TI " << TI << " W " << W << " TPD " <<  tpd << std::endl;
    double F = TI * (n-k-1)/n
             + init._a * TE * (2*n-k-1)/(2*n)
             + init._b * TT * (k-1+n)/(2*n);

    double L = I * (n - k - 1)
             + init._c*E
             + (init._e)/(n-k+1)*W;

    g_value = F + L;
}*/

double emili::pfsp::BeamSearchHeuristic::bs_node::calcW()
{
    std::vector<int>::iterator iter = unscheduled.begin();
    std::vector<int> ctimes(completionTimes);
    double W = TT;
    for(;iter!=unscheduled.end();++iter)
    {
        int cj = *iter;
        ctimes[1] += init.pi[cj][1];
        for(int i=2;i<=m;i++)
        {
            ctimes[i] = std::max(ctimes[i],ctimes[i-1])+init.pi[cj][i];
        }
        //ctimes[m] = ctimes[m]+init.pi[cj][m];
        W = W + std::max(ctimes[m]-(int)init.dueDates[cj],0);
    }
    //std::vector< int > asd = scheduled;
    //asd.insert(asd.end(),unscheduled.begin(),unscheduled.end());
    //assert(W == init.pis.computeObjectiveFunction(asd));
    return W;
}

void emili::pfsp::BeamSearchHeuristic::bs_node::evaluateNode()
{

    completionTimes[1] += init.pi[kjob][1];

    double I = 0.0;
    double km = k;
    double nd = n;
    double md = m;
    for(int i=2;i<=m;i++)
    {
        double iden = md / (i-1+(km-1)*(md-i+1)/(nd-2));
        double idletime = std::max(completionTimes[i-1]-completionTimes[i],0);
        I = I + idletime*iden;
        completionTimes[i] = std::max(completionTimes[i],completionTimes[i-1])+init.pi[kjob][i];
    }    
    T = std::max(completionTimes[m]-(int)init.dueDates[kjob],0);
    TT = TT + T;

    if(k == (n-1))
    {
        g_value = TT;

    }
    else
    {
        double E = std::max((int)init.dueDates[kjob]-completionTimes[m],0);
        TE = TE + E;
        TI = TI + I;
       // std::cout << "K " << k << " F " << father->kjob << " J " << kjob << " I " << I << " E " << E << std::endl;
        double F = TI * (nd-km-1)/nd
                + init._a * TE * (2*nd-km-1)/(2*nd)
                + init._b * TT * (km-1+nd)/(2*nd);

        //double W = calcW();
        L = I * (nd - km - 1) + E*init._c;// + ((init._e)/(n-k+1))*W;
        g_value = F ;//+ L ;
        //if(k < 4)

    }
}

double emili::pfsp::BeamSearchHeuristic::bs_node::calcG(double W)
{
   // double F = g_value;
    if(k< (n-1))
    {
        //L = L + (init._e)/(n-k+1)*W;
        double nd = n-k;
        L = L + (init._e)/(nd+1)*W;
        g_value = g_value+L;
    }

    //std::cout << "K " << k << " F "<< father->kjob << " J " << kjob << " TE " << TE << " TT " << TT << " TI " << TI
    //        << " W " << W << " L " <<  L << " G " << g_value << " F "<< F << std::endl;

    return g_value;
}

/*
void emili::pfsp::BeamSearchHeuristic::bs_node::buildTree()
{
    std::vector< bs_node* > offspring;

    int size = n-k;
    for(int i=0;i < size; i++)
    {
        bs_node* b = new bs_node(*this,i);
        offspring.push_back(b);
    }

    int childnum = std::min(size,init._gamma);
    for(int i=0; i < childnum ; i++)
    {
        double best = offspring[0]->g_value;
        int b = 0;
        for(int j=0; j < size; j++)
        {
            if(best > offspring[j]->g_value)
            {
                b = j;
                best = offspring[j]->g_value;
            }
        }
        children.push_back(offspring[b]);
        offspring.erase(offspring.begin()+b);
        size--;
    }
    int remaining_nodes = offspring.size();
    for(int i=0;i < remaining_nodes;i++)
    {
        delete offspring[i];
    }
}*/

void emili::pfsp::BeamSearchHeuristic::bs_node::buildChildren()
{
    int size = n-k;
    for(int i=0;i < size; i++)
    {
        bs_node* b = new bs_node(*this,i);
        children.push_back(b);
    }
}


/*emili::pfsp::BeamSearchHeuristic::bs_node* emili::pfsp::BeamSearchHeuristic::bs_node::generateSequence()
{
    count++;
    emili::pfsp::BeamSearchHeuristic::bs_node* best = this;
    if((n-k-1) == 0)
    {
/
        std::cout << "K: " << k << " L: " << kjob << " G: "
                  << g_value << " un: "<< unscheduled.size()
                  << " sc: " << scheduled.size() << " this " << (long)this  << std::endl;
        what++;
        if(what>1)
        std::cout << "what! " << what << std::endl;
        /
        scheduled.push_back(unscheduled[0]);
        this->g_value = init.pis.computeObjectiveFunction(scheduled);
    }
    else
    {
        buildTree();
        if(k==n-2)
        {
            kc++;
        }
        best = (children[0])->generateSequence();
        int size = children.size();
        for(int i = 1 ;i<size;i++)
        {
            emili::pfsp::BeamSearchHeuristic::bs_node* current = (children[i])->generateSequence();
            if(*current < *best)
            {
                delete best;
                best = current;
            }
            else
            {
               if(best!=current)
                 delete current;
            }
        }
    }
    return best;
}
*/
