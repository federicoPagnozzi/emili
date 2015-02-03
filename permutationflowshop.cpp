#include "permutationflowshop.h"
#include <cstdlib>
#include <climits>
#include <string>
#include <sstream>
#include <assert.h>


std::vector< int > inline neh(std::vector< int >& partial,int nbJobs,emili::pfsp::PermutationFlowShop& pis)
{
    std::vector< int >  sol(nbJobs+1,0);
    sol[1] = partial[1];
    sol[2] = partial[2];
    int wt = pis.computeWT(sol);
    sol[1] = partial[2];
    sol[2] = partial[1];
    int wt2 = pis.computeWT(sol);
    if(wt2>wt){
        sol[1] = partial[1];
        sol[2] = partial[2];
    }
    for(int i=3; i <= nbJobs; i++)
    {
        int candidate = partial[i];
        int bestpos = 1;
        sol.insert(sol.begin()+1,candidate);
        int wt_min = pis.computeWT(sol);
        for(int j = 2;j<i;j++)
        {
            sol.erase(sol.begin()+(j-1));
            sol.insert(sol.begin()+j,candidate);
            int wt = pis.computeWT(sol);
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

            int mS=pis.computeWT(_fsp);//compute_total_wt(_fsp,2);
            if(pis.computeWT(solTMP)<mS){//compute_total_wt(solTMP,2)<mS){
                _fsp[1]=solTMP[1];
                _fsp[2]=solTMP[2];
            }

            for(int k=3;k<=N;k++){
                    min=10000000;
                for(int r=1; r<=k; r++){

                    for(int h=1; h<r; h++)
                        solTMP[h]=_fsp[h];
                    solTMP[r]=_fsp[k];
                    for(int h=r+1; h<=k; h++)
                        solTMP[h]=_fsp[h-1];

                    tmp=pis.computeWT(solTMP);//compute_total_wt(solTMP,k+1);
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

int generateRndPos(int min, int max)
{
  return (  emili::generateRandomNumber()%max + min );
}

double emili::pfsp::PermutationFlowShop::evaluateSolution(emili::Solution& solution)
{
    std::vector< int >* current_solution = (std::vector< int >*) solution.getRawData();
    double p = instance.computeWT(*current_solution);
    solution.setSolutionValue(p);
    return p;
}

int emili::pfsp::PermutationFlowShop::getNjobs(){
    return instance.getNbJob();
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

int emili::pfsp::PermutationFlowShop::computeWT(std::vector< int > & partial_solution)
{
    return instance.computeWT(partial_solution);
}

int emili::pfsp::PermutationFlowShop::computeWT(std::vector< int > & partial_solution,int size)
{
    return instance.computeWT(partial_solution,size);
}

int emili::pfsp::PermutationFlowShop::computeWT(vector<int> &sol,vector<int>& prevJob,int job,vector<int>& previousMachineEndTime)
{
    return instance.computeWT(sol,prevJob,job,previousMachineEndTime);
}

const void* emili::pfsp::PermutationFlowShopSolution::getRawData()const
{
    return &solution;
}

void emili::pfsp::PermutationFlowShopSolution::setRawData(const void *data)
{
    std::vector < int >* data_vector = (std::vector< int >*) data;
    this->solution = *data_vector;
}

emili::pfsp::PermutationFlowShopSolution::~PermutationFlowShopSolution()
{
 /*nothing to delete*/
}

emili::Solution* emili::pfsp::PfspInitialSolution::generateEmptySolution()
{
    std::vector< int > empty(pis.getNjobs()+1);
    return new emili::pfsp::PermutationFlowShopSolution(empty);

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
    pis.evaluateSolution(*s);
    return s;
}

emili::Solution* emili::pfsp::SlackConstructor::construct(Solution *partial)
{
    int nbJobs = pis.getNjobs();
    std::vector< int >  sol(nbJobs+1, 0);
    std::vector< int > * p = (std::vector<int >*) partial->getRawData();
    int Ci = pis.computeMS(*p);
    vector<bool> assigned(nbJobs+1, false);
    for(std::vector< int >::const_iterator iter = p->begin();iter !=  p->end();++iter)
    {
        assigned[*iter] = true;
    }
//    assigned[0] = false;
    for (int var = 0; var < nbJobs; ++var) {
        int kish = p->operator [](var+1);
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
            p->operator [](var+1) = minJ;
            assigned[minJ] = true;
            Ci = pis.computeMS(*p);
        }
        //cout << "partial makespan: " << Ci << std::endl;
    }
    sol = *p;
    PermutationFlowShopSolution* s = new PermutationFlowShopSolution(sol);
    pis.evaluateSolution(*s);
    return s;

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

emili::Solution* emili::pfsp::NEHSlackConstructor::construct(Solution *partial)
{
    int nbJobs = pis.getNjobs();
    std::vector< int >  sol(nbJobs+1, 0);
    std::vector< int > * p = (std::vector<int >*) partial->getRawData();
    std::vector<int > part = slack_construct(*p,nbJobs,pis);
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
     std::vector< int > * p = (std::vector< int > *) solutioon->getRawData();
    std::vector< int > des(*p);
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
    std::vector< int > * p = (std::vector< int > *) solutioon->getRawData();
    std::vector< int > des (*p);
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



emili::Solution* emili::pfsp::SOAPerturbation::perturb(Solution *solution)
{
    //emili::iteration_increment();

    int index;
    int min;
    int k,tmp,ind;
    std::vector< int > * p = (std::vector< int > *) solution->getRawData();
    std::vector< int > removed;
    std::vector< int > solPartial(*p);
    //std::cout << "partial size " << solPartial.size() << std::endl;
    int size = p->size();
    std::vector< int > solTMP(size,0);
    int sizePartial;
    int sops = solPartial.size()-1;
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
            tmp = instance.computeWT(solTMP,sizePartial);

            if(tmp<min){
                min=tmp;
                ind=r;
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

emili::Solution* emili::pfsp::PfspDestructorTest::destruct(Solution *solutioon)
{
     std::vector< int > * p = (std::vector< int > *) solutioon->getRawData();
    std::vector< int > des(*p);
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
        std::vector< int > * solution = (std::vector< int > *) value->getRawData();
        std::vector < int > newsol(*solution);
        int sol_i = newsol[start_position];

        newsol.erase(newsol.begin()+start_position);
        newsol.insert(newsol.begin()+end_position,sol_i);
        long int new_value = instance.computeWT(newsol);

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
        std::vector< int > * solution = (std::vector< int > *) value->getRawData();
        std::vector < int > newsol(*solution);
        int sol_i = newsol[start_position];
        newsol.erase(newsol.begin()+start_position);
        newsol.insert(newsol.begin()+end_position,sol_i);
        long int new_value = instance.computeWT(newsol);

        return new emili::pfsp::PermutationFlowShopSolution(new_value,newsol);
    }

}

emili::Solution* emili::pfsp::PfspInsertNeighborhood::computeStep(emili::Solution* value)
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
        std::vector< int > * solution = (std::vector< int > *) value->getRawData();
        std::vector < int > newsol(*solution);
        int sol_i = newsol[start_position];
        newsol.erase(newsol.begin()+start_position);
        newsol.insert(newsol.begin()+end_position,sol_i);
        long int new_value = instance.computeWT(newsol);

        return new emili::pfsp::PermutationFlowShopSolution(new_value,newsol);
    }
}

emili::Solution* emili::pfsp::PfspInsertNeighborhood::random(Solution *currentSolution)
{

    std::vector < int > newsol = *((std::vector<int>*)currentSolution->getRawData());
    int best_i = (emili::generateRandomNumber()%100)+1;
    int best_j = (emili::generateRandomNumber()%100)+1;
    int sol_i = newsol[best_i];
    newsol.erase(newsol.begin()+best_i);
    newsol.insert(newsol.begin()+best_j,sol_i);
    long int value = instance.computeWT(newsol);
    return new emili::pfsp::PermutationFlowShopSolution(value,newsol);
}

emili::Solution* emili::pfsp::PfspForwardInsertNeighborhood::random(Solution *currentSolution)
{

    std::vector < int > newsol = *((std::vector<int>*)currentSolution->getRawData());
    int best_i = (emili::generateRandomNumber()%100)+1;    
    int best_j = (emili::generateRandomNumber()%(100-best_i+1))+best_i;
    int sol_i = newsol[best_i];
    newsol.erase(newsol.begin()+best_i);
    newsol.insert(newsol.begin()+best_j,sol_i);
    long int value = instance.computeWT(newsol);
    return new emili::pfsp::PermutationFlowShopSolution(value,newsol);
}

emili::Solution* emili::pfsp::PfspBackwardInsertNeighborhood::random(Solution *currentSolution)
{

    std::vector < int > newsol = *((std::vector<int>*)currentSolution->getRawData());
    int best_i = (emili::generateRandomNumber()%100)+1;
    int best_j = (emili::generateRandomNumber()%best_i)+1;
    int sol_i = newsol[best_i];
    newsol.erase(newsol.begin()+best_i);
    newsol.insert(newsol.begin()+best_j,sol_i);
    long int value = instance.computeWT(newsol);
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
        std::vector< int > * solution = (std::vector< int > *) value->getRawData();
        std::vector < int > newsol(*solution);
        std::swap(newsol[start_position],newsol[end_position]);
        long int new_value = instance.computeWT(newsol);
        return new emili::pfsp::PermutationFlowShopSolution(new_value,newsol);
    }
}

emili::Solution* emili::pfsp::PfspExchangeNeighborhood::random(Solution *currentSolution)
{

    std::vector < int > newsol = *((std::vector<int>*)currentSolution->getRawData());
    int best_i = (emili::generateRandomNumber()%100)+1;
    int best_j = (emili::generateRandomNumber()%100)+1;
    std::swap(newsol[best_i],newsol[best_j]);
    long int value = instance.computeWT(newsol);
    return new emili::pfsp::PermutationFlowShopSolution(value,newsol);
}


void emili::pfsp::PfspExchangeNeighborhood::reset()
{    
    start_position = 1;
    end_position = 2;
}

emili::Solution* emili::pfsp::PfspTransposeNeighborhood::computeStep(emili::Solution* value)
{
    if(sp_iterations >= njobs)
    {
        return nullptr;
    }
    else
    {
        sp_iterations++;
        start_position = (start_position%njobs)+1;
        std::vector< int > * solution = (std::vector< int > *) value->getRawData();
        std::vector < int > newsol(*solution);
        int endpos = start_position<njobs?start_position+1:1;        
        std::swap(newsol[start_position],newsol[endpos]);
        long int new_value = instance.computeWT(newsol);
        return new emili::pfsp::PermutationFlowShopSolution(new_value,newsol);
    }
}

emili::Solution* emili::pfsp::PfspTransposeNeighborhood::random(Solution *currentSolution)
{

    std::vector < int > newsol = *((std::vector<int>*)currentSolution->getRawData());
    int best_i = (emili::generateRandomNumber()%100);
    std::swap(newsol[best_i],newsol[best_i+1]);
    long int value = instance.computeWT(newsol);
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
    std::vector < int >* sol_data = (std::vector < int >*)solution->getRawData();
    std::vector < int > perturbed(*sol_data);
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

size_t emili::pfsp::PfspTabuHashMemory::calc_hash(std::vector< int > * v )
{
    std::hash< std::string> hasher;
    std::ostringstream oss;
    for(std::vector< int > ::iterator iter = v->begin();iter!= v->end();++iter )
    {
        oss << *iter;
    }
    return hasher(oss.str());
}

bool emili::pfsp::PfspTabuHashMemory::tabu_check(emili::Solution* toCheck)
{  
    std::vector< int > * v = (std::vector< int > *)toCheck->getRawData();
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
    std::vector < int > * v = (std::vector < int >* ) solution->getRawData();
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

    std::vector< int > *value = (std::vector< int >*) toCheck->getRawData();

    for(std::vector< std::vector<int > >::iterator iter = tabuVector.begin();iter!=tabuVector.end();++iter)
    {
        std::vector< int > t = *iter ;
        bool ret = false;
        for(int i=0; i< t.size() ; i++)
        {
            if(t[i]!=(*value)[i])
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
        std::vector< int >* solution_vector = (std::vector< int >*) solution->getRawData();
        if(tt_index < this->tabutenure){
            tabuVector.push_back(*solution_vector);
            tt_index++;
        }
        else
        {
            tabuVector.erase(tabuVector.begin());            
            tabuVector.push_back(*solution_vector);

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
