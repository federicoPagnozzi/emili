#include "permutationflowshop.h"
#include <cstdlib>
#include <climits>
#include <string>
#include <sstream>
#include <assert.h>
#include <algorithm>
#include <limits>

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

    std::sort(seq.begin()+1,seq.end(),[stds](int i1,int i2){return stds[i1] < stds[i2];});

    return seq;

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
        std::vector< int > tas(jobs+1,0);
        temp.push_back(0);
        for(int i=1;i<=jobs;i++ )
        {
            int tai= 0;

            for (int j=k;j<=machines;j++)
            {
                tai += ( machines - j + 1 ) * priorities[i][j];
            }
            tai = tai/prob.getPriority(i);
            tas[i] = tai;
            temp.push_back(i);
        }


        std::sort(temp.begin(),temp.end(),[tas](int i1,int i2){
                                                                if(tas[i1]==tas[i2] && i2!=0)
                                                                    return i1>i2;
                                                                else
                                                                   return tas[i1] < tas[i2];
        });

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
    std:vector< int > best_sol(partial);
    best_sol.insert(best_sol.end(),removed.begin(),removed.end());
    int wbest = prob.computeObjectiveFunction(best_sol);
    int machines = prob.getNmachines();
    int jobs = removed.size();
    const std::vector < std::vector < long int > >& priorities = prob.getProcessingTimesMatrix();

    for(int k=1; k<=machines; k++)
    {
        std::vector< int > temp;
        std::vector< int > tas(prob.getNjobs()+1,0);
        //temp.push_back(0);
        for(int i=0;i<jobs;i++ )
        {
            int tai= 0;

            for (int j=k;j<=machines;j++)
            {
                tai += ( machines - j + 1 ) * priorities[removed[i]][j];
            }
            tai = tai/prob.getPriority(removed[i]);
            tas[removed[i]] = tai;
            temp.push_back(removed[i]);
        }

    std::sort(temp.begin(),temp.end(),[tas](int i1,int i2){
                                                            if(tas[i1]==tas[i2])
                                                                return i1>i2;
                                                            else
                                                               return tas[i1] < tas[i2];
    });

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
    std::vector < int > res;
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
            pis.computeTAmatrices(solTMP,head,tail,solTMP.size());
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
                        itdd += fil-head[i][ptbl-1]-pmatrix[kk][i]+pmatrix[sptbl][i]+max(filp-fil,0);
                        //itdd += fil-head[i][ptbl-1]+pmatrix[sptbl][i]+max(filp-fil,0);
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
    vector<bool> assigned(nbJobs+1, false);
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
            /*ITK calculation*/
            float wjk = nmac/(j+k*((nmac-j)/(float)(njobs-2)));

            float tcomp = std::max(i_comp-lastJobCompletionTimes[j],0.0f);


            ITik += wjk*tcomp;

            i_comp = tcomp>0?i_comp:lastJobCompletionTimes[j];
            i_comp += ctimesMatrix[i_job][j];

            /* P job completion time calculation*/
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



int emili::pfsp::PermutationFlowShop::computeObjectiveFunction(vector<int> &sol,vector<int>& prevJob,int job,vector<int>& previousMachineEndTime)
{
    return instance.computeWT(sol,prevJob,job,previousMachineEndTime);
}

int emili::pfsp::PermutationFlowShop::computeObjectiveFunction(vector<int> &sol, vector<vector<int> > &previousMachineEndTimeMatrix, int start_i, int end_i)
{
    return instance.computeWT(sol,previousMachineEndTimeMatrix,start_i,end_i);
}

void emili::pfsp::PermutationFlowShop::computeWTs(vector<int> &sol,vector<int>& prevJob,int job,vector<int>& previousMachineEndTime)
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


const std::vector< std::vector < long int > > & emili::pfsp::PermutationFlowShop::getProcessingTimesMatrix()
{
    return instance.getProcessingTimesMatrix();
}

emili::pfsp::PermutationFlowShopSolution::~PermutationFlowShopSolution()
{
 /*nothing to delete*/
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

int emili::pfsp::PFSP_E::computeObjectiveFunction(std::vector<int> &partial_solution)
{
    return instance.computeE(partial_solution);
}

int emili::pfsp::PFSP_E::computeObjectiveFunction(std::vector<int> &partial_solution, int size)
{
    return instance.computeE(partial_solution,size);
}

#pragma region NEWCODE

// Due to the behaviour of normal Flowshop modules, we override getNmachines to get 
// the number of stages. That in their internal behaviour will be what is expected. 
int emili::pfsp::HybridFlowShop::getNmachines()
{
	return instance.getNbStages();
}

std::vector<int> emili::pfsp::HybridFlowShop::getStages()
{
	return instance.getStages();
}

int emili::pfsp::HybridFlowShop::computeMS(std::vector< int > & partial_solution)
{
	return instance.computeHMS(partial_solution);
}

int emili::pfsp::HybridFlowShop::computeMS(std::vector< int > & partial_solution, int size)
{
	return instance.computeHMS(partial_solution, size);
}

int emili::pfsp::HFSP_MS::computeObjectiveFunction(std::vector<int> &partial_solution)
{
	return instance.computeHMS(partial_solution);
}

int emili::pfsp::HFSP_MS::computeObjectiveFunction(std::vector<int> &partial_solution, int size)
{
	return instance.computeHMS(partial_solution, size);
}

int emili::pfsp::HFSP_TCT::computeObjectiveFunction(std::vector<int> &partial_solution)
{
	return instance.computeHTCT(partial_solution);
}

int emili::pfsp::HFSP_TCT::computeObjectiveFunction(std::vector<int> &partial_solution, int size)
{
	return instance.computeHTCT(partial_solution, size);
}

int emili::pfsp::HFSP_WT::computeObjectiveFunction(std::vector<int> &partial_solution)
{
	return instance.computeHWT(partial_solution);
}

int emili::pfsp::HFSP_WT::computeObjectiveFunction(std::vector<int> &partial_solution, int size)
{
	return instance.computeHWT(partial_solution, size);
}

int emili::pfsp::HFSP_WE::computeObjectiveFunction(std::vector<int> &partial_solution)
{
	return instance.computeHWE(partial_solution);
}

int emili::pfsp::HFSP_WE::computeObjectiveFunction(std::vector<int> &partial_solution, int size)
{
	return instance.computeHWE(partial_solution, size);
}

int emili::pfsp::HFSP_WET::computeObjectiveFunction(std::vector<int> &partial_solution)
{
	return instance.computeHWET(partial_solution);
}

int emili::pfsp::HFSP_WET::computeObjectiveFunction(std::vector<int> &partial_solution, int size)
{
	return instance.computeHWET(partial_solution, size);
}
/*
#pragma region SSD

int emili::pfsp::HFSP_SDST_MS::computeObjectiveFunction(std::vector<int> &partial_solution)
{
	return instance.computeHSDSTMS(partial_solution);
}

int emili::pfsp::HFSP_SDST_MS::computeObjectiveFunction(std::vector<int> &partial_solution, int size)
{
	return instance.computeHSDSTMS(partial_solution, size);
}

int emili::pfsp::HFSP_SDST_TCT::computeObjectiveFunction(std::vector<int> &partial_solution)
{
	return instance.computeHSDSTTCT(partial_solution);
}

int emili::pfsp::HFSP_SDST_TCT::computeObjectiveFunction(std::vector<int> &partial_solution, int size)
{
	return instance.computeHSDSTTCT(partial_solution, size);
}

int emili::pfsp::HFSP_SDST_WT::computeMS(std::vector<int> &partial_solution)
{
	return instance.computeHSDSTMS(partial_solution);
}

int emili::pfsp::HFSP_SDST_WT::computeMS(std::vector<int> &partial_solution, int size)
{
	return instance.computeHSDSTMS(partial_solution, size);
}

int emili::pfsp::HFSP_SDST_WT::computeObjectiveFunction(std::vector<int> &partial_solution)
{
	return instance.computeHSDSTWT(partial_solution);
}

int emili::pfsp::HFSP_SDST_WT::computeObjectiveFunction(std::vector<int> &partial_solution, int size)
{
	return instance.computeHSDSTWT(partial_solution, size);
}

int emili::pfsp::HFSP_SDST_WE::computeObjectiveFunction(std::vector<int> &partial_solution)
{
	return instance.computeHSDSTWE(partial_solution);
}

int emili::pfsp::HFSP_SDST_WE::computeObjectiveFunction(std::vector<int> &partial_solution, int size)
{
	return instance.computeHSDSTWE(partial_solution, size);
}

int emili::pfsp::HFSP_SDST_WET::computeObjectiveFunction(std::vector<int> &partial_solution)
{
	return instance.computeHSDSTWET(partial_solution);
}

int emili::pfsp::HFSP_SDST_WET::computeObjectiveFunction(std::vector<int> &partial_solution, int size)
{
	return instance.computeHSDSTWET(partial_solution, size);
}
#pragma endregion SSD
*/
#pragma endregion NEWCODE

int emili::pfsp::PFSP_MS::computeObjectiveFunction(std::vector<int> &partial_solution)
{
    return instance.computeMS(partial_solution);
}

int emili::pfsp::PFSP_MS::computeObjectiveFunction(std::vector<int> &partial_solution, int size)
{
    return instance.computeMS(partial_solution,size);
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
  vector<bool> alreadyTaken(nbJobs+1, false); // nbJobs elements with value false
  vector<int > choosenNumber(nbJobs+1, 0);

  int nbj;
  int rnd, i, j, nbFalse;

  nbj = 0;
  for (i = nbJobs; i >= 1; --i)
  {
    rnd = generateRndPos(1, i);
    nbFalse = 0;

    /* find the rndth cell with value = false : */
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

emili::Solution* emili::pfsp::PfspSlackInitialSolution::generate()
{
    int nbJobs = pis.getNjobs();
    std::vector< int >  sol(nbJobs+1, 0);
    int Ci  = 0;
    vector<bool> assigned(nbJobs+1, false);
    vector<int> partial(nbJobs+1,0);
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
    vector<bool> assigned(nbJobs+1, false); // vector that indicates if a job it's already assigned to the partial solution
    vector<int> partial(nbJobs+1,0); // partial solution
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
    order.push_back(0);
    for (int i = 1; i <= njobs; ++i) {
        int tpti = 0;
        for (int k = 1; k <= nmac; ++k) {
            tpti += ptm[i][k];
        }
        tpt[i] = tpti;
        order.push_back(i);
    }
    std::sort(order.begin(),order.end(),[tpt](int i1,int i2){return tpt[i1] > tpt[i2];});
    order.erase(order.begin()+njobs);
    order.insert(order.begin(),0);
    order = neh2(order,njobs,pis);
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
    order.push_back(0);
    for (int i = 1; i <= njobs; ++i) {
        int tpti = 0;
        for (int k = 1; k <= nmac; ++k) {
            tpti += ptm[i][k];
        }
        tpt[i] = tpti;
        order.push_back(i);
    }
    std::sort(order.begin(),order.end(),[tpt](int i1,int i2){return tpt[i1] > tpt[i2];});
    order.erase(order.begin()+njobs);
    order.insert(order.begin(),0);
    order = nehff(order,njobs,pis);
    PermutationFlowShopSolution* s = new PermutationFlowShopSolution(order);
    pis.evaluateSolution(*s);
    return s;
}

emili::Solution* emili::pfsp::SlackConstructor::construct(Solution *partial)
{
    int nbJobs = pis.getNjobs();
    std::vector< int >  sol(nbJobs+1, 0);
    std::vector< int >&  p = ((emili::pfsp::PermutationFlowShopSolution*)partial)->getJobSchedule();
    int Ci = pis.computeMS(p);
    vector<bool> assigned(nbJobs+1, false);
    for(std::vector< int >::const_iterator iter = p.begin();iter !=  p.end();++iter)
    {
        assigned[*iter] = true;
    }
//    assigned[0] = false;
    for (int var = 0; var < nbJobs; ++var) {
        int kish = p[var+1];
        if(kish==0)
        {
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
             p[var+1] = minJ;
            assigned[minJ] = true;
            Ci = pis.computeMS(p);
        }
        //cout << "partial makespan: " << Ci << std::endl;
    }
    sol = p;
    PermutationFlowShopSolution* s = new PermutationFlowShopSolution(sol);
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

emili::Solution* emili::pfsp::SlackConstructor::constructFull()
{
    int nbJobs = pis.getNjobs();
    std::vector< int >  sol(nbJobs+1, 0);
    int Ci  = 0;
    vector<bool> assigned(nbJobs+1, false);
    vector<int> partial(nbJobs+1,0);
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

emili::Solution* emili::pfsp::MNEH::generate()
{
    std::vector< int > initial = std_start_sequence(pis);
    //initial = rz_improvement_phase(initial,pis);
    initial = neh2(initial,pis.getNjobs(),pis);

    PermutationFlowShopSolution* s = new PermutationFlowShopSolution(initial);
    pis.evaluateSolution(*s);
    return s;
}

std::vector< int > inline lr_solution_sequence(int start,std::vector< int > u,std::vector< int > initial,emili::pfsp::PermutationFlowShop& pis)
{
    std::vector< float > fndx = lr_index(initial,u,pis);
    std::sort(u.begin(),u.end(),[fndx](int i1,int i2){return fndx[i1] < fndx[i2];});
    initial.push_back(u[start]);
    u.erase(u.begin());

    int usize = u.size()-1;

    for(int i=0; i< usize;i++)
    {
        fndx = lr_index(initial,u,pis);
        std::sort(u.begin(),u.end(),[fndx](int i1,int i2){return fndx[i1] < fndx[i2];});
        initial.push_back(u[0]);
        u.erase(u.begin());
    }

    initial.push_back(u[0]);

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


emili::Solution* emili::pfsp::NEHSlackConstructor::construct(Solution *partial)
{
    int nbJobs = pis.getNjobs();
    std::vector< int >  sol(nbJobs+1, 0);

    std::vector<int > part = slack_construct(((emili::pfsp::PermutationFlowShopSolution*)partial)->getJobSchedule(),nbJobs,pis);
    sol = neh(part,nbJobs,pis);
    PermutationFlowShopSolution* s = new PermutationFlowShopSolution(sol);
    pis.evaluateSolution(*s);
    return s;

}



emili::Solution* emili::pfsp::NEHSlackConstructor::constructFull()
{
    int nbJobs = pis.getNjobs();
    std::vector< int >  sol(nbJobs+1, 0);
    int Ci  = 0;
    vector<bool> assigned(nbJobs+1, false);
    vector<int> partial(nbJobs+1,0);
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
    sol = neh(partial,nbJobs,pis);
    PermutationFlowShopSolution* s = new PermutationFlowShopSolution(sol);
    pis.evaluateSolution(*s);
    return s;
}

emili::Solution* emili::pfsp::PfspDestructor::destruct(Solution *solutioon)
{

    std::vector< int > des(((emili::pfsp::PermutationFlowShopSolution*)solutioon)->getJobSchedule());
    int size = des.size();
    int start_position = emili::generateRandomNumber()%(size-1);
    int num_postion = emili::generateRandomNumber()%(size-start_position-1);
    des.erase(des.begin()+start_position,des.begin()+start_position+num_postion);
    des.insert(des.begin()+start_position,num_postion,0);

    int nbJobs = instance.getNjobs();
    std::vector<int> res = slack_construct(des,nbJobs,instance);
    res = neh2(res,nbJobs,instance);
    emili::Solution* s = new emili::pfsp::PermutationFlowShopSolution(res);
    instance.evaluateSolution(*s);
    return s;
}

emili::Solution* emili::pfsp::SOADestructor::destruct(Solution *solutioon)
{

    std::vector< int > des (((emili::pfsp::PermutationFlowShopSolution*)solutioon)->getJobSchedule());
    int size = des.size();
    for (int var = 0; var < d; ++var) {
        int num = emili::generateRandomNumber()%(size-1)+1;
        if(des[num]!=0){
            des[num]=0;
        }else{
            var--;
        }
    }
    int nbJobs = instance.getNjobs();
    std::vector<int> res = slack_construct(des,nbJobs,instance);
    res = neh2(res,nbJobs,instance);
    emili::Solution* s = new emili::pfsp::PermutationFlowShopSolution(res);
    instance.evaluateSolution(*s);
    return s;
}

emili::Solution* emili::pfsp::NRZPertubation::perturb(Solution *solution)
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

emili::Solution* emili::pfsp::TMIIGPertubation::perturb(Solution *solution)
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
        tblist[solPartial[index]].push_back(solPartial[index-1]);
        if(tblist[solPartial[index]].size() > tbsize)
        {
            tblist[solPartial[index]].erase(tblist[solPartial[index]].begin());
        }
        solPartial.erase(solPartial.begin() + index);
        sops--;
    }

    for(int l=0;l<removed.size();l++){
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
        sops++;
        //std::cout << "end insert " << solPartial.size() << std::endl;
    }


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
        sops++;
        //std::cout << "end insert " << solPartial.size() << std::endl;
    }


    emili::pfsp::PermutationFlowShopSolution* s = new emili::pfsp::PermutationFlowShopSolution(min,solPartial);
    //instance.evaluateSolution(*s);
    //std::vector< std::vector< int > >  etm = std::vector< std::vector< int > >(instance.getNmachines()+1,std::vector<int>(instance.getNjobs()+1,0));
    //clock_t start = clock();
    //int new_value =  instance.computeObjectiveFunction(solPartial,etm,1,instance.getNjobs()+1);

    //s->setSolutionValue(new_value);
    return s;
}

emili::Solution* emili::pfsp::RSPertubation::perturb(Solution *solution)
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
        solPartial.erase(solPartial.begin() + index);
        sops--;
    }



    for(int l=0;l<removed.size();l++){
        k=removed[l];

        min = std::numeric_limits<int>::max();

        instance.computeTAmatrices(solPartial,head,tail,sops+1);

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
            for (int i = 2; i <= instance.getNmachines(); ++i) {
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
            }

        }
        solPartial.insert(solPartial.begin()+ind,k);
        sops++;
        //std::cout << "end insert " << solPartial.size() << std::endl;
    }


    emili::pfsp::PermutationFlowShopSolution* s = new emili::pfsp::PermutationFlowShopSolution(min,solPartial);
    instance.evaluateSolution(*s);
    return s;
}

emili::Solution* emili::pfsp::RSffPertubation::perturb(Solution *solution)
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

        min = std::numeric_limits<int>::max();

        instance.computeTAmatrices(solPartial,head,tail,sops+1);
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
        /*TIE BREAKING FF*/
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
                            itdd += fil-head[i][ptbl-1]-pmatrix[k][i]+pmatrix[sptbl][i]+max(filp-fil,0);
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
        /*END TIE BREAKING*/
        solPartial.insert(solPartial.begin()+ind,k);
        sops++;
        //std::cout << "end insert " << solPartial.size() << std::endl;
    }


    emili::pfsp::PermutationFlowShopSolution* s = new emili::pfsp::PermutationFlowShopSolution(min,solPartial);
    instance.evaluateSolution(*s);
    return s;
}

emili::Solution* emili::pfsp::IgLsPertubation::perturb(Solution *solution)
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
      emili::pfsp::PermutationFlowShopSolution* s = new emili::pfsp::PermutationFlowShopSolution(solPartial);
      s->setSolutionValue(instance.computeObjectiveFunction(solPartial,sops));
      s =(emili::pfsp::PermutationFlowShopSolution*) ls->search(s);
      solPartial = s->getJobSchedule();

    for(int l=0;l<removed.size();l++){
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
        sops++;
        //std::cout << "end insert " << solPartial.size() << std::endl;
    }

    delete s;
    s = new emili::pfsp::PermutationFlowShopSolution(solPartial);
    instance.evaluateSolution(*s);
    return s;
}

emili::Solution* emili::pfsp::RSLSPertubation::perturb(Solution *solution)
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
      emili::pfsp::PermutationFlowShopSolution* s = new emili::pfsp::PermutationFlowShopSolution(solPartial);
      s->setSolutionValue(instance.computeObjectiveFunction(solPartial,sops));
      s =(emili::pfsp::PermutationFlowShopSolution*) ls->search(s);
      solPartial = s->getJobSchedule();

    for(int l=0;l<removed.size();l++){
        k=removed[l];
        min = std::numeric_limits<int>::max();
        instance.computeTAmatrices(solPartial,head,tail,sops+1);
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
            for (int i = 2; i <= instance.getNmachines(); ++i) {
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
            tmp = c_max;

            if(tmp<min){
                min=tmp;
                ind=r;
            }

        }
        solPartial.insert(solPartial.begin()+ind,k);
        sops++;
        //std::cout << "end insert " << solPartial.size() << std::endl;
    }

    delete s;
    s = new emili::pfsp::PermutationFlowShopSolution(solPartial);
    instance.evaluateSolution(*s);
    return s;
}


emili::Solution* emili::pfsp::PfspDestructorTest::destruct(Solution *solutioon)
{

    std::vector< int > des(((emili::pfsp::PermutationFlowShopSolution*)solutioon)->getJobSchedule());
    int size = des.size();    
    int hsize = size/2;
    int start_position =  emili::generateRandomNumber()%(hsize-1);
    int start_position_2 = emili::generateRandomNumber()%(hsize-1)+hsize;
    int num_postion = size/10;
    int end_pos = (start_position+num_postion)<hsize?(start_position+num_postion):hsize;
    int end_pos_2 = (start_position_2+num_postion)<size?(start_position_2+num_postion):size;
    //std::cout << "start_position , num_position " << start_position << " , " << end_pos << std::endl;
    //std::cout << "start_position2 , num_position2 " << start_position_2 << " , " << end_pos_2 << std::endl;
    des.erase(des.begin()+start_position,des.begin()+end_pos);
    des.insert(des.begin()+start_position,(end_pos-start_position),0);
    des.erase(des.begin()+start_position_2,des.begin()+end_pos_2);
    des.insert(des.begin()+start_position_2,(end_pos_2-start_position_2),0);
    int nbJobs = istance.getNjobs();
    std::vector<int> res = slack_construct(des,nbJobs,istance);
    res = neh(res,nbJobs,istance);
    emili::Solution* s = new emili::pfsp::PermutationFlowShopSolution(res);
    istance.evaluateSolution(*s);
    return s;
}

/*
Construct the solution inserting one job at a time, by always selecting
the one that minimizes the weighted earlyness.

The weighted earlyness of job Ji is computed as wi · (di − Ci).

Note: the solution is constructed incrementally, and at each iteration
Ci corresponds to the makespan of the partial solution

void emili::pfsp::PfspSlackInitialSolution::slackInitial(std::vector<int> & sol)
{
    int noj = pis.getNjobs();
    int Ci  = 0;
    vector<bool> assigned(noj+1, false);
    vector<int> partial(noj+1,0);
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
    /*No counters to reset*/
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
    computeTAmatrices(sol);
    return emili::Neighborhood::NeighborhoodIterator(this,base);
}

emili::Neighborhood::NeighborhoodIterator emili::pfsp::HeavilyApproximatedTaillardAcceleratedInsertNeighborhood::begin(Solution *base)
{
    ep_iterations = 1;
    sp_iterations = 1;
    std::vector< int > sol(((emili::pfsp::PermutationFlowShopSolution*)base)->getJobSchedule());
    sol.erase(sol.begin()+start_position);
    computeHead(sol);    
    return emili::Neighborhood::NeighborhoodIterator(this,base);
}

emili::Neighborhood::NeighborhoodIterator emili::pfsp::Natx2Neighborhood::begin(Solution *base)
{
    ep_iterations = 1;
    sp_iterations = 1;
    std::vector< int > sol(((emili::pfsp::PermutationFlowShopSolution*)base)->getJobSchedule());
    sol.erase(sol.begin()+start_position);
    computeHead(sol);
    value_wt = base->getSolutionValue();
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
    computeHead(sol);
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

    if(sp_iterations >= njobs)
    {
        return nullptr;
    }
    else
    {
        if(end_position < (start_position-1)){
            ep_iterations++;
            /*if(ep_iterations == sp_iterations){
                ep_iterations++;
                end_position++;
            }*/
        }
        else
        {
            sp_iterations++;
            ep_iterations = 1;
            start_position = ((start_position)%njobs)+1;
            end_position = 0;

        }
        end_position = ((end_position)%njobs)+1;
        std::vector < int > newsol(((emili::pfsp::PermutationFlowShopSolution*)value)->getJobSchedule());
        int sol_i = newsol[start_position];

        newsol.erase(newsol.begin()+start_position);
        newsol.insert(newsol.begin()+end_position,sol_i);
        long int new_value = pis.computeObjectiveFunction(newsol);
        return new emili::pfsp::PermutationFlowShopSolution(new_value,newsol);
    }
}

emili::Solution* emili::pfsp::PfspForwardInsertNeighborhood::computeStep(emili::Solution* value)
{
    if(sp_iterations >= njobs)
    {
        return nullptr;
    }
    else
    {
        if(end_position < njobs){
            ep_iterations++;
            if(ep_iterations == sp_iterations){
                ep_iterations++;
                end_position++;
            }
        }
        else
        {
            sp_iterations++;
            ep_iterations = sp_iterations;
            start_position = ((start_position)%njobs)+1;
            end_position = start_position;

        }
        end_position = ((end_position)%njobs)+1;     
        std::vector < int > newsol(((emili::pfsp::PermutationFlowShopSolution*)value)->getJobSchedule());
        int sol_i = newsol[start_position];
        newsol.erase(newsol.begin()+start_position);
        newsol.insert(newsol.begin()+end_position,sol_i);
        long int new_value = pis.computeObjectiveFunction(newsol);


        return new emili::pfsp::PermutationFlowShopSolution(new_value,newsol);
    }

}

emili::Solution* emili::pfsp::PfspInsertNeighborhood::computeStep(emili::Solution* value)
{
    emili::iteration_increment();
    if(sp_iterations >= njobs)
    {
        return nullptr;
    }
    else
    {
        if(ep_iterations < njobs){
            ep_iterations++;
            if(ep_iterations == sp_iterations){
                ep_iterations++;
                end_position++;
            }
        }
        else
        {
            sp_iterations++;
            ep_iterations = 1;
            start_position = ((start_position)%njobs)+1;
            //end_position = start_position;

        }        
        end_position = ((end_position)%njobs)+1;        
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
            computeTAmatrices(newsol);

        }
        end_position = ((end_position)%njobs)+1;



        newsol.insert(newsol.begin()+end_position,sol_i);
        long int c_cur = head[1][end_position-1]+pmatrix[sol_i][1];
        long int c_max = c_cur+tail[1][end_position];
        for (int i = 2; i <= pis.getNmachines(); ++i) {
            int c_pm = head[i][end_position-1];
            if(c_pm < c_cur)
            {
                c_cur = c_cur + pmatrix[sol_i][i];
            }
            else
            {
                c_cur = c_pm + pmatrix[sol_i][i];
            }
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

emili::Solution* emili::pfsp::ApproximatedTaillardAcceleratedInsertNeighborhood::computeStep(emili::Solution *value)
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
            computeTAmatrices(newsol);

        }
        end_position = ((end_position)%njobs)+1;



        newsol.insert(newsol.begin()+end_position,sol_i);
        //std::vector< int > ins_pos(nmac+1,0);
        //int ins_pos[nmac+1];
		std::vector< int > ins_pos(nmac);
        long int c_cur = head[1][end_position-1]+pmatrix[sol_i][1];
        ins_pos[1] = c_cur;
        long int c_max = c_cur+tail[1][end_position];
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
       /* for(int k=end_position+1; k< njobs; k++)
        {
            pre_c_cur = pre_c_cur + pmatrix[newsol[k]][nmac];
            wt += (std::max(pre_c_cur - pis.getDueDate(newsol[k]), 0L) * pis.getPriority(newsol[k]));
        }


*/

      //  int value_wt = value->getSolutionValue();
  //      if(wt < value_wt)
   //     {

            for(int k=end_position+1; k< njobs; k++)
            {
                int job = newsol[k];
                pre_c_cur = ins_pos[1] + pmatrix[job][1];
                ins_pos[1] = pre_c_cur;
                for(int m=2; m <= nmac ; m++)
                {
                    int c_pm = ins_pos[m];
                    if(c_pm < pre_c_cur)
                    {
                        pre_c_cur += pmatrix[job][m];
                    }
                    else
                    {
                        pre_c_cur = c_pm + pmatrix[job][m];
                    }
                    ins_pos[m] = pre_c_cur;
                }
                pre_wt += (std::max(pre_c_cur - pis.getDueDate(newsol[k]), 0L) * pis.getPriority(newsol[k]));
     /*           if(pre_wt > value_wt)
                 {
                    // cout << "exited at " << k << std::endl;

                     news->setSolutionValue(pre_wt);
                     return news;
                  } */
            }
             //   std::cout << pre_wt << " "<< pis.computeObjectiveFunction(newsol);
            //assert( pre_wt == pis.computeObjectiveFunction(newsol));

            value->setSolutionValue(pre_wt);
           //pis.evaluateSolution(*news);
        /*   if(news->getSolutionValue() > value->getSolutionValue())
           {
               emili::iteration_increment();
           }*/
        //}


        //long int old_v  = pis.computeObjectiveFunction(newsol);
        //std::cout << c_max << " - " << old_v << std::endl;
        //assert(c_max == old_v);
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
            computeHead(newsol);

        }
        end_position = ((end_position)%njobs)+1;



        newsol.insert(newsol.begin()+end_position,sol_i);
        //std::vector< int > ins_pos(nmac+1,0);
        //int ins_pos[nmac+1];
		std::vector< int > ins_pos(nmac + 1);
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
                    if(c_pm < pre_c_cur)
                    {
                        pre_c_cur += pmatrix[job][m];
                    }
                    else
                    {
                        pre_c_cur = c_pm + pmatrix[job][m];
                    }
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
        /*   if(news->getSolutionValue() > value->getSolutionValue())
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
            computeHead(newsol);

        }
        end_position = ((end_position)%njobs)+1;



        newsol.insert(newsol.begin()+end_position,sol_i);
        //std::vector< int > ins_pos(nmac+1,0);
        //int ins_pos[nmac+1];
		std::vector< int > ins_pos(nmac + 1);
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
                    if(c_pm < pre_c_cur)
                    {
                        pre_c_cur += pmatrix[job][m];
                    }
                    else
                    {
                        pre_c_cur = c_pm + pmatrix[job][m];
                    }
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

emili::Solution* emili::pfsp::Natx2Neighborhood::computeStep(emili::Solution *value)
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
            computeHead(newsol);
        }
        end_position = ((end_position)%njobs)+1;
        newsol.insert(newsol.begin()+end_position,sol_i);

        //int ins_pos[nmac+1];
		std::vector< int > ins_pos(nmac + 1);
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
                    if(c_pm < pre_c_cur)
                    {
                        pre_c_cur += pmatrix[job][m];
                    }
                    else
                    {
                        pre_c_cur = c_pm + pmatrix[job][m];
                    }
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

emili::Solution* emili::pfsp::EatxNeighborhood::computeStep(emili::Solution *value)
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
            computeHead(newsol);

        }
        end_position = ((end_position)%njobs)+1;



        newsol.insert(newsol.begin()+end_position,sol_i);
        //std::vector< int > ins_pos(nmac+1,0);
        //int ins_pos[nmac+1];
		std::vector< int > ins_pos(nmac + 1);
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
                    if(c_pm < pre_c_cur)
                    {
                        pre_c_cur += pmatrix[job][m];
                    }
                    else
                    {
                        pre_c_cur = c_pm + pmatrix[job][m];
                    }
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
            computeHead(newsol);

        }
        end_position = ((end_position)%njobs)+1;



        newsol.insert(newsol.begin()+end_position,sol_i);
        //std::vector< int > ins_pos(nmac+1,0);
        //int ins_pos[nmac+1];
        //int ins_pos_t[nmac+1];
		std::vector< int > ins_pos(nmac + 1);
		std::vector< int > ins_pos_t(nmac + 1);
        long int c_cur = head[1][end_position-1]+pmatrix[sol_i][1];
        ins_pos[1] = c_cur;
        ins_pos_t[1] = c_cur;
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
           ins_pos_t[i] = c_cur;
        }



        int wt = (std::max(c_cur - pis.getDueDate(sol_i), 0L) * pis.getPriority(sol_i));
        for (int j = 1; j< end_position; ++j )
        {
            wt += (std::max((long int)head[nmac][j] - pis.getDueDate(newsol[j]), 0L) * pis.getPriority(newsol[j]));
        }

        int pre_wt = wt;
        long int pre_c_cur=0;
        if(end_position > aptre)
        for(int k=end_position+1; k<= njobs; k++)
        {
            int job = newsol[k];
            pre_c_cur = ins_pos_t[aplev] + pmatrix[job][aplev];
            ins_pos_t[aplev] = pre_c_cur;
            for(int m=aplev+1; m <= nmac ; m++)
            {
                int c_pm = ins_pos_t[m];
                if(c_pm < pre_c_cur)
                {
                    pre_c_cur += pmatrix[job][m];
                }
                else
                {
                    pre_c_cur = c_pm + pmatrix[job][m];
                }
                ins_pos_t[m] = pre_c_cur;
            }
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
                    if(c_pm < pre_c_cur)
                    {
                        pre_c_cur += pmatrix[job][m];
                    }
                    else
                    {
                        pre_c_cur = c_pm + pmatrix[job][m];
                    }
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


emili::Solution* emili::pfsp::NoIdleAcceleratedInsertNeighborhood::computeStep(emili::Solution *value)
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
            pis.computeNoIdleTAmatrices(newsol,head,tail);
        }
        end_position = ((end_position)%njobs)+1;



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
        //int ins_pos[nmac+1];
		std::vector< int > ins_pos(nmac + 1);
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

    if(sp_iterations >= njobs)
    {
        return nullptr;
    }
    else
    {
        if(ep_iterations < njobs){
            ep_iterations++;
            if(ep_iterations == sp_iterations){
                ep_iterations++;
                end_position++;
            }
        }
        else
        {
            sp_iterations++;
            ep_iterations = 1;
            start_position = ((start_position)%njobs)+1;
            //end_position = start_position;

        }
        end_position = ((end_position)%njobs)+1;

       std::vector < int >& newsol = ((emili::pfsp::PermutationFlowShopSolution*)value)->getJobSchedule();

        // first insert
        int sol_i = newsol[start_position];
        newsol.erase(newsol.begin()+start_position);
        newsol.insert(newsol.begin()+end_position,sol_i);
        // second insert
        sol_i = newsol[start_position];
        newsol.erase(newsol.begin()+start_position);
        newsol.insert(newsol.begin()+end_position,sol_i);

        long int new_value = pis.computeObjectiveFunction(newsol);
        value->setSolutionValue(new_value);
        return value;
    }
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
            computeHead(newsol);
            newsol.insert(newsol.begin()+start_position,sol_i);
        }
        end_position = ((end_position)%njobs)+1;

        std::swap(newsol[start_position],newsol[end_position]);
        //std::vector< int > ins_pos(nmac+1,0);

        int ms_pos = start_position>end_position?end_position:start_position;
        sol_i = newsol[ms_pos];

        //int ins_pos[nmac+1];
		std::vector< int > ins_pos(nmac + 1);
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
        long int pppre_c_cur = ins_pos[nmac-2];
        long int ppre_c_cur = ins_pos[nmac-1];
        long int pre_c_cur = c_cur;

        for (int j = 1; j< ms_pos; ++j )
        {
            wt += (std::max((long int)head[nmac][j] - pis.getDueDate(newsol[j]), 0L) * pis.getPriority(newsol[j]));
        }

        int pre_wt = wt;

        if(ms_pos > (njobs/2))
            for(int k=ms_pos+1; k<= njobs; k++)
            {
                pppre_c_cur = pppre_c_cur + pmatrix[newsol[k]][nmac-2];
                ppre_c_cur = std::max(ppre_c_cur,pppre_c_cur) + pmatrix[newsol[k]][nmac-1];
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
                    if(c_pm < pre_c_cur)
                    {
                        pre_c_cur += pmatrix[job][m];
                    }
                    else
                    {
                        pre_c_cur = c_pm + pmatrix[job][m];
                    }
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
            computeHead(newsol);
            newsol.insert(newsol.begin()+start_position,sol_i);
        }
        end_position = ((end_position)%njobs)+1;

        std::swap(newsol[start_position],newsol[end_position]);

        int ms_pos = start_position>end_position?end_position:start_position;
        sol_i = newsol[ms_pos];

        //int ins_pos[nmac+1];
		std::vector< int > ins_pos(nmac + 1);
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
                if(c_pm < pre_c_cur)
                {
                    pre_c_cur += pmatrix[job][m];
                }
                else
                {
                    pre_c_cur = c_pm + pmatrix[job][m];
                }
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


void emili::pfsp::PfspExchangeNeighborhood::reset()
{    
    start_position = 1;
    end_position = 2;
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
    emili::iteration_increment();
    return true;
}

void emili::pfsp::PfspTerminationIterations::reset()
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
    oss << "[ ";

    for (int i = 1; i < solution.size(); ++i)
      oss << solution[i] << " " ;

    oss << " ]";
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
    lastMove = neigh.lastMove();
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
