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

            int mS=pis.computeObjectiveFunction(_fsp);//compute_total_wt(_fsp,2);
            if(pis.computeObjectiveFunction(solTMP,solTMP.size())<mS){//compute_total_wt(solTMP,2)<mS){
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

                    tmp=pis.computeObjectiveFunction(solTMP,solTMP.size());//compute_total_wt(solTMP,k+1);
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

int emili::pfsp::NIPFSP_MS::computeObjectiveFunction(std::vector<int> &partial_solution)
{
    return instance.computeNIMS(partial_solution);
}

int emili::pfsp::NIPFSP_MS::computeObjectiveFunction(std::vector<int> &partial_solution, int size)
{
    return instance.computeNIMS(partial_solution,size);
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
        //cout << "partial makespan: " << Ci << std::endl;
    }
    //int partial_w = pis.computeWT(partial);
    sol = neh2(partial,nbJobs,pis);
    PermutationFlowShopSolution* s = new PermutationFlowShopSolution(sol);
    //pis.evaluateSolution(*s);
    std::vector< std::vector< int > > etm = std::vector< std::vector< int > >(pis.getNmachines()+1,std::vector<int>(pis.getNjobs()+1,0));
    //clock_t start = clock();
    int new_value =  pis.computeObjectiveFunction(sol,etm,1,pis.getNjobs()+1);

    s->setSolutionValue(new_value);
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
    int k,tmp,ind;
    std::vector< int > removed;
    std::vector< int > solPartial(((emili::pfsp::PermutationFlowShopSolution*)solution)->getJobSchedule());
    //std::cout << "partial size " << solPartial.size() << std::endl;
    int size = solPartial.size();
    std::vector< int > solTMP(size,0);
    int sizePartial;
    int sops = solPartial.size()-1;
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

    sizePartial = solPartial.size();
    for(int l=0;l<removed.size();l++){
        k=removed[l];
        min = std::numeric_limits<int>::max();

        for(int r=1; r<sizePartial; r++){
            if(std::find(tblist[k].begin(),tblist[k].end(),solPartial[r]) == tblist[k].end())
            {
            for(int h=1; h<r; h++)
                solTMP[h]=solPartial[h];
            solTMP[r]=k;
            for(int h=r+1; h<=sizePartial; h++)
                solTMP[h]=solPartial[h-1];


            //tmp=compute_total_wt(solTMP,sizePartial+1);
            //                  std::cout << "start perturb" << std::endl;
            //check why plus 1
            tmp = instance.computeObjectiveFunction(solTMP,sizePartial);

            if(tmp<min){
                min=tmp;
                ind=r;
            }
            }
        }
        solPartial.insert(solPartial.begin()+ind,k);
        sizePartial++;
        //std::cout << "end insert " << solPartial.size() << std::endl;
    }


    emili::Solution* s = new emili::pfsp::PermutationFlowShopSolution(solPartial);
    instance.evaluateSolution(*s);
    return s;
}

emili::Solution* emili::pfsp::SOAPerturbation::perturb(Solution *solution)
{
    //emili::iteration_increment();

    int index;
    int min;
    int k,tmp,ind;

    std::vector< int > removed;
    std::vector< int > solPartial(((emili::pfsp::PermutationFlowShopSolution*)solution)->getJobSchedule());
    //std::cout << "partial size " << solPartial.size() << std::endl;
    int size = solPartial.size();
    std::vector< int > solTMP(size,0);
    int sizePartial;
    int sops = size-1;
    for(int k = 0; k < d; k++) {
        index = (emili::generateRandomNumber()%sops)+1;
        //std::cout << index << " " ;//<< std::endl;
        removed.push_back(solPartial[index]);
        solPartial.erase(solPartial.begin() + index);
        sops--;
    }

    sizePartial = solPartial.size();
    for(int l=0;l<removed.size();l++){
        k=removed[l];
        min = std::numeric_limits<int>::max();

        for(int r=1; r<sizePartial; r++){

            for(int h=1; h<r; h++)
                solTMP[h]=solPartial[h];
            solTMP[r]=k;
            for(int h=r+1; h<=sizePartial; h++)
                solTMP[h]=solPartial[h-1];


            //tmp=compute_total_wt(solTMP,sizePartial+1);
            //                  std::cout << "start perturb" << std::endl;
            //check why plus 1
            tmp = instance.computeObjectiveFunction(solTMP,sizePartial);

            if(tmp<min){
                min=tmp;
                ind=r;
            }

        }
        solPartial.insert(solPartial.begin()+ind,k);
        sizePartial++;
        //std::cout << "end insert " << solPartial.size() << std::endl;
    }


    emili::pfsp::PermutationFlowShopSolution* s = new emili::pfsp::PermutationFlowShopSolution(solPartial);
    //instance.evaluateSolution(*s);
    std::vector< std::vector< int > >  etm = std::vector< std::vector< int > >(instance.getNmachines()+1,std::vector<int>(instance.getNjobs()+1,0));
    //clock_t start = clock();
    int new_value =  instance.computeObjectiveFunction(solPartial,etm,1,instance.getNjobs()+1);

    s->setSolutionValue(new_value);
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
    pis.computeTAmatrices(sol,head,tail);
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

emili::Neighborhood::NeighborhoodIterator emili::pfsp::PfspTransposeNeighborhood::begin(emili::Solution *base)
{    
    //sp_iterations = 1;
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
        std::vector < int > newsol(((emili::pfsp::PermutationFlowShopSolution*)value)->getJobSchedule());
        int sol_i = newsol[start_position];
        newsol.erase(newsol.begin()+start_position);
        newsol.insert(newsol.begin()+end_position,sol_i);
        long int new_value = pis.computeObjectiveFunction(newsol);

        return new emili::pfsp::PermutationFlowShopSolution(new_value,newsol);
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
        std::vector < int > newsol(((emili::pfsp::PermutationFlowShopSolution*)value)->getJobSchedule());
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
            pis.computeTAmatrices(newsol,head,tail);

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
            c_max = c_max>c_can?c_max:c_can;
        }
        //long int old_v  = pis.computeObjectiveFunction(newsol);
        //std::cout << c_max << " - " << old_v << std::endl;
        //assert(c_max == old_v);
        return new emili::pfsp::PermutationFlowShopSolution(c_max,newsol);
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
        std::vector < int > newsol(((emili::pfsp::PermutationFlowShopSolution*)value)->getJobSchedule());
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

        return new emili::pfsp::PermutationFlowShopSolution(c_max,newsol);
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
        std::vector < int > newsol(((emili::pfsp::PermutationFlowShopSolution*)value)->getJobSchedule());
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
        int nmac = pis.getNmachines();
        std:vector< long int > partials(njobs+1,0);
        for(int i = 1; i<end_position; i++)
        {
            partials[i] = head[nmac][i];
        }

        std::vector< int > ins_pos(nmac+1,0);

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

        partials[end_position] = c_cur;

        for(int k=end_position+1; k<= njobs; k++)
        {
            std::vector < std::vector < int > >& tail = tails[k];

            long int c_max = ins_pos[1]+tail[1][end_position];
            for (int i = 2; i <= pis.getNmachines(); ++i) {
                long int c_can = (ins_pos[i]+tail[i][end_position]);
                c_max = c_max>c_can?c_max:c_can;
            }
            partials[k] = c_max;
        }


        int wt = 0;
        for (int j = 1; j<= njobs; ++j )
            wt += (std::max(partials[j] - pis.getDueDate(newsol[j]), 0L) * pis.getPriority(newsol[j]));


        //long int old_v  = pis.computeObjectiveFunction(newsol);
        //std::cout << c_max << " - " << old_v << std::endl;
        //assert(c_max == old_v);
     //   int robj = pis.computeObjectiveFunction(newsol);
        //std::cout << "wt " << wt << " <-> "<< robj << std::endl;

       // assert(wt==robj);
        return new emili::pfsp::PermutationFlowShopSolution(wt,newsol);
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

        std::vector < int > newsol(((emili::pfsp::PermutationFlowShopSolution*)value)->getJobSchedule());

        // first insert
        int sol_i = newsol[start_position];
        newsol.erase(newsol.begin()+start_position);
        newsol.insert(newsol.begin()+end_position,sol_i);
        // second insert
        sol_i = newsol[start_position];
        newsol.erase(newsol.begin()+start_position);
        newsol.insert(newsol.begin()+end_position,sol_i);

        long int new_value = pis.computeObjectiveFunction(newsol);
        return new emili::pfsp::PermutationFlowShopSolution(new_value,newsol);
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
        std::vector < int > newsol(((emili::pfsp::PermutationFlowShopSolution*)value)->getJobSchedule());
        std::swap(newsol[start_position],newsol[end_position]);
        long int new_value = pis.computeObjectiveFunction(newsol);
        return new emili::pfsp::PermutationFlowShopSolution(new_value,newsol);
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
        std::vector < int > newsol(((emili::pfsp::PermutationFlowShopSolution*)value)->getJobSchedule());
        int endpos = start_position<njobs?start_position+1:1;
        //int endpos = (start_position%njobs)+1;
        std::swap(newsol[start_position],newsol[endpos]);
        //clock_t start = clock();
        long int new_value = pis.computeObjectiveFunction(newsol);
        //std::cout <<  (clock()-start) << ", " ;//<< std::endl;
        return new emili::pfsp::PermutationFlowShopSolution(new_value,newsol);
    }
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
        std::vector < int > newsol(p->getJobSchedule());
        int endpos =  start_position<njobs?start_position+1:1;
        std::swap(newsol[start_position],newsol[endpos]);

        long int new_value = 0;
        //clock_t start = clock();

            //clock_t start = clock();


            //new_value =  pis.computeObjectiveFunction(newsol,etm,start_position,endpos);

         //std::cout <<  (clock()-start) << ", " ;//<< std::endl;
            //assert(new_value == pis.computeObjectiveFunction(newsol));

        return new emili::pfsp::PermutationFlowShopSolution(new_value,newsol);
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
