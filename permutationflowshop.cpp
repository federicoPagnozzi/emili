#include "permutationflowshop.h"
#include <cstdlib>
#include <string>
#include <sstream>

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
        int minE = 4*10e9;
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

    std::vector< int >* csol = (std::vector< int >*)currentSolution->getRawData();

    emili::pfsp::PermutationFlowShopSolution* ret = computeStep(*csol,currentSolution->getSolutionValue());

    return ret;
}

void emili::pfsp::PfspNeighborhood::reset()
{
    /*No counters to reset*/
}

emili::pfsp::PermutationFlowShopSolution* emili::pfsp::PfspBestImprovExchangeNeighborhood::computeStep(std::vector< int > & solution,double value_c)
{
    std::vector < int > newsol = solution;
    long int value = value_c;
    int best_i = 1;
    int best_j = 1;
    for(int i = 1; i< njobs; i++)
    {        
        int posb = newsol[i];
        for(int j = i+1; j<=njobs; j++)
        {
            newsol[i] = newsol[j];
            newsol[j] = posb;
            long int new_value = instance.computeWT(newsol);
            if(new_value < value){
                best_i = i;
                best_j = j;
                value = new_value;
            }
            newsol[j] = newsol[i];
            newsol[i] = posb;
        }
    }
    int posb = newsol[best_i];
    newsol[best_i] = newsol[best_j];
    newsol[best_j] = posb;

    if(value < value_c)
    {
      return new emili::pfsp::PermutationFlowShopSolution(value,newsol);
    }
    else
    {
        return nullptr;
    }
}
emili::Solution* emili::pfsp::PfspBestImprovExchangeNeighborhood::random(Solution *currentSolution)
{

    std::vector < int > newsol = *((std::vector<int>*)currentSolution->getRawData());
    int best_i = (emili::generateRandomNumber()%100)+1;
    int best_j = (emili::generateRandomNumber()%100)+1;
    int posb = newsol[best_i];
    newsol[best_i] = newsol[best_j];
    newsol[best_j] = posb;
    long int value = instance.computeWT(newsol);
    return new emili::pfsp::PermutationFlowShopSolution(value,newsol);
}

/*
emili::Solution* emili::pfsp::PfspBestImprovExchangeNeighborhood::computeStep(std::vector< int > & solution,double value_c)
{
     clock_t time = clock();
    std::vector < int > newsol;
    long int value = value_c;
    std::vector < int > best = solution;
    for(int i = 1; i< njobs; i++)
    {
        for(int j = 1; j<=njobs; j++)
        {
            newsol = solution;
            int posb = newsol[i];
            newsol[i] = newsol[j];
            newsol[j] = posb;
            long int new_value = instance.computeWT(newsol);
            if(new_value < value){
                best = newsol;
                value = new_value;
            }

        }
    }
    double time_elapsed = (double)(clock()-time)/CLOCKS_PER_SEC;
    cout << ".>stepTime : " << time_elapsed << std::endl;
    return new emili::pfsp::PermutationFlowShopSolution(value,best);
}

emili::Solution* emili::pfsp::PfspInsertNeighborhood::computeStep(std::vector< int > & solution,double value)
{
    if(start)
    {
        current = solution;
        current_value = current[start_position];
        current.erase(current.begin()+start_position);
        start = false;
    }

    if(start_position >= njobs)
    {
        return nullptr;
    }
    else
    {
        if(end_position < njobs){
            end_position++;
        }
        else
        {
            start_position++;
            end_position = 1;
            current = solution;
            current_value = current[start_position];
            current.erase(current.begin()+start_position);
        }
        std::vector < int > newsol;
        newsol = current;
        newsol.insert(newsol.begin()+end_position,current_value);
        long int new_value = instance.computeWT(newsol);
        return new emili::pfsp::PermutationFlowShopSolution(new_value,newsol);
    }
}
*/
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
    sp_iterations = 1;
    return emili::Neighborhood::NeighborhoodIterator(this,base);
}

emili::pfsp::PermutationFlowShopSolution* emili::pfsp::PfspInsertNeighborhood::computeStep(std::vector< int > & solution,double value)
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

        }        
        end_position = ((end_position)%njobs)+1;

        std::vector < int > newsol;
        int sol_i = solution[start_position];
        newsol = solution;
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





emili::pfsp::PermutationFlowShopSolution* emili::pfsp::PfspExchangeNeighborhood::computeStep(std::vector<int> &solution, double value)
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
        std::vector < int > newsol;
        newsol = solution;
        int posb = newsol[start_position];
        newsol[start_position] = newsol[end_position];
        newsol[end_position] = posb;
        long int new_value = instance.computeWT(newsol);
        return new emili::pfsp::PermutationFlowShopSolution(new_value,newsol);
    }
}

emili::Solution* emili::pfsp::PfspExchangeNeighborhood::random(Solution *currentSolution)
{

    std::vector < int > newsol = *((std::vector<int>*)currentSolution->getRawData());
    int best_i = (emili::generateRandomNumber()%100)+1;
    int best_j = (emili::generateRandomNumber()%100)+1;
    int posb = newsol[best_i];
    newsol[best_i] = newsol[best_j];
    newsol[best_j] = posb;
    long int value = instance.computeWT(newsol);
    return new emili::pfsp::PermutationFlowShopSolution(value,newsol);
}


void emili::pfsp::PfspExchangeNeighborhood::reset()
{    
    start_position = 1;
    end_position = 2;
}

emili::pfsp::PermutationFlowShopSolution* emili::pfsp::PfspTransposeNeighborhood::computeStep(std::vector<int> &solution, double value)
{
    if(sp_iterations >= njobs)
    {
        return nullptr;
    }
    else
    {
        sp_iterations++;
        start_position = (start_position%njobs)+1;
        std::vector < int > newsol;
        newsol = solution;
        int posb = newsol[start_position];
        int endpos = start_position<njobs?start_position+1:1;
        newsol[start_position] = newsol[endpos];
        newsol[endpos] = posb;
        long int new_value = instance.computeWT(newsol);
        return new emili::pfsp::PermutationFlowShopSolution(new_value,newsol);
    }
}

emili::Solution* emili::pfsp::PfspTransposeNeighborhood::random(Solution *currentSolution)
{

    std::vector < int > newsol = *((std::vector<int>*)currentSolution->getRawData());
    int best_i = (emili::generateRandomNumber()%100);

    int posb = newsol[best_i];
    newsol[best_i] = newsol[best_i+1];
    newsol[best_i+1] = posb;
    long int value = instance.computeWT(newsol);
    return new emili::pfsp::PermutationFlowShopSolution(value,newsol);
}

void emili::pfsp::PfspTransposeNeighborhood::reset()
{
    //sp_iterations = 1;
    start_position = 1;
}


emili::pfsp::PermutationFlowShopSolution* emili::pfsp::PfspBestImprovInsertNeighborhood::computeStep(std::vector< int > & solution,double valuec)
{
    std::vector < int > newsol;
    long int value = valuec;
    std::vector < int > best = solution;
    for(int i = 1; i< njobs; i++)
    {
        int sol_i = solution[i];
        for(int j = 1; j<=njobs; j++)
        {
            newsol = solution;
            newsol.erase(newsol.begin()+i);
            newsol.insert(newsol.begin()+j,sol_i);
            long int new_value = instance.computeWT(newsol);
            if(new_value < value){
                best = newsol;
                value = new_value;
            }

        }
    }
    if(value < valuec)
    {
      return new emili::pfsp::PermutationFlowShopSolution(value,best);
    }
    else
    {
        return nullptr;
    }
}

emili::Solution* emili::pfsp::PfspBestImprovInsertNeighborhood::random(Solution *currentSolution)
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


 /* std::vector< int > emili::pfsp::PfspBestImprovInsertNeighborhood::computeStep(std::vector< int > & solution)
{
    std::vector < int > newsol;
    long int value = instance.computeWT(solution);
    std::vector < int > best = solution;
    for(int i = 1; i< njobs; i++)
    {
        int sol_i = solution[i];
        for(int j = 1; j<=njobs; j++)
        {
            newsol = solution;
            newsol.erase(newsol.begin()+i);
            newsol.insert(newsol.begin()+j,sol_i);
            long int new_value = instance.computeWT(newsol);
            if(new_value < value){
                best = newsol;
                value = new_value;
            }

        }
    }
    return best;
}*/


emili::pfsp::PermutationFlowShopSolution* emili::pfsp::PfspBestImprovTransposeNeighborhood::computeStep(std::vector<int> &solution,double valuec)
{
    int njobs = pis.getNjobs()-1;
    std::vector < int > newsol;
    std::vector < int > currentSol(solution);
    PfspInstance& instance = pis.getInstance();
    long int value = valuec;
    std::vector < int > best = currentSol;
    for(int i = 1; i< njobs; i++)
        {
            newsol = currentSol;
            int posb = newsol[i];
            newsol[i] = newsol[i+1];
            newsol[i+1] = posb;
            long int new_value = instance.computeWT(newsol);
            if(new_value < value){
                best = newsol;
                value = new_value;
            }

        }
    if(value < valuec)
    {
      return new emili::pfsp::PermutationFlowShopSolution(value,best);
    }
    else
    {
        return nullptr;
    }
}

emili::Solution* emili::pfsp::PfspBestImprovTransposeNeighborhood::random(Solution *currentSolution)
{

    std::vector < int > newsol = *((std::vector<int>*)currentSolution->getRawData());
    int best_i = (emili::generateRandomNumber()%100);

    int posb = newsol[best_i];
    newsol[best_i] = newsol[best_i+1];
    newsol[best_i+1] = posb;
    long int value = instance.computeWT(newsol);
    return new emili::pfsp::PermutationFlowShopSolution(value,newsol);
}


emili::pfsp::PermutationFlowShopSolution* emili::pfsp::PfspFirstImprovExchangeNeighborhood::computeStep(std::vector< int > & solution,double valuec)
{
    int njobs = pis.getNjobs();
    std::vector < int > newsol;
    std::vector < int > currentSol(solution);
    PfspInstance& instance = pis.getInstance();
    long int value = valuec;
    std::vector < int > best = currentSol;
    for(int i = start_position; i< njobs; i++)
        for(int j = end_position; j<=njobs; j++)
        {
            newsol = currentSol;
            int posb = newsol[i];
            newsol[i] = newsol[j];
            newsol[j] = posb;
            long int new_value = instance.computeWT(newsol);
            if(new_value < value){
                best = newsol;
                value = new_value;
                start_position = i;
                end_position = j;
                return new emili::pfsp::PermutationFlowShopSolution(value,best);
            }

        }

    if(value < valuec)
    {
      return new emili::pfsp::PermutationFlowShopSolution(value,best);
    }
    else
    {
        return nullptr;
    }

}

emili::Solution* emili::pfsp::PfspFirstImprovExchangeNeighborhood::random(Solution *currentSolution)
{

    std::vector < int > newsol = *((std::vector<int>*)currentSolution->getRawData());
    int best_i = (emili::generateRandomNumber()%100)+1;
    int best_j = (emili::generateRandomNumber()%100)+1;
    int posb = newsol[best_i];
    newsol[best_i] = newsol[best_j];
    newsol[best_j] = posb;
    long int value = pis.getInstance().computeWT(newsol);
    return new emili::pfsp::PermutationFlowShopSolution(value,newsol);
}

emili::pfsp::PermutationFlowShopSolution* emili::pfsp::PfspFirstImprovInsertNeighborhood::computeStep(std::vector< int > & solution,double valuec)
{
    int njobs = pis.getNjobs();
    std::vector < int > newsol;
    std::vector < int > currentSol(solution);
    PfspInstance& instance = pis.getInstance();
    long int value = valuec;
    std::vector < int > best = currentSol;
    for(int i = start_position; i< njobs; i++)
        for(int j = end_position; j<=njobs; j++)
        {
            newsol = currentSol;
            newsol.erase(newsol.begin()+i);
            newsol.insert(newsol.begin()+j,solution[i]);
            long int new_value = instance.computeWT(newsol);
            if(new_value < value){
                best = newsol;
                value = new_value;
                start_position = i;
                end_position = j;
                return new emili::pfsp::PermutationFlowShopSolution(value,best);
            }

        }

    if(value < valuec)
    {
      return new emili::pfsp::PermutationFlowShopSolution(value,best);
    }
    else
    {
        return nullptr;
    }
}

emili::Solution* emili::pfsp::PfspFirstImprovInsertNeighborhood::random(Solution *currentSolution)
{

    std::vector < int > newsol = *((std::vector<int>*)currentSolution->getRawData());
    int best_i = (emili::generateRandomNumber()%100)+1;
    int best_j = (emili::generateRandomNumber()%100)+1;
    int sol_i = newsol[best_i];
    newsol.erase(newsol.begin()+best_i);
    newsol.insert(newsol.begin()+best_j,sol_i);
    long int value = pis.getInstance().computeWT(newsol);
    return new emili::pfsp::PermutationFlowShopSolution(value,newsol);
}


emili::pfsp::PermutationFlowShopSolution* emili::pfsp::PfspFirstImprovTransposeNeighborhood::computeStep(std::vector<int> &solution,double valuec)
{
    int njobs = pis.getNjobs()-1;
    std::vector < int > newsol;
    std::vector < int > currentSol(solution);
    PfspInstance& instance = pis.getInstance();
    long int value = valuec;
    std::vector < int > best = currentSol;
    for(int i = start_position; i< njobs; i++)
        {
            newsol = currentSol;
            int posb = newsol[i];
            newsol[i] = newsol[i+1];
            newsol[i+1] = posb;
            long int new_value = instance.computeWT(newsol);
            if(new_value < value){
                best = newsol;
                value = new_value;
                start_position = i;
                return new emili::pfsp::PermutationFlowShopSolution(value,best);
            }

        }

    if(value < valuec)
    {
      return new emili::pfsp::PermutationFlowShopSolution(value,best);
    }
    else
    {
        return nullptr;
    }
}

emili::Solution* emili::pfsp::PfspFirstImprovTransposeNeighborhood::random(Solution *currentSolution)
{

    std::vector < int > newsol = *((std::vector<int>*)currentSolution->getRawData());
    int best_i = (emili::generateRandomNumber()%100);

    int posb = newsol[best_i];
    newsol[best_i] = newsol[best_i+1];
    newsol[best_i+1] = posb;

    long int value = pis.getInstance().computeWT(newsol);
    return new emili::pfsp::PermutationFlowShopSolution(value,newsol);
}


void emili::pfsp::PfspInsertNeighborhood::reset()
{
    start_position = 1;
    end_position = 1;
    start = true;
}

void emili::pfsp::PfspFirstImprovExchangeNeighborhood::reset()
{

    start_position = 1;
    end_position = 1;
}

void emili::pfsp::PfspFirstImprovInsertNeighborhood::reset()
{
    start_position = 1;
    end_position = 1;
}

void emili::pfsp::PfspFirstImprovTransposeNeighborhood::reset()
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

    if(chance <= 70)
    {
        return c;
    }
    else
    {
        return c2;
    }

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

bool emili::pfsp::PfspTabuInsertNeighborhood::notTabu(NeighborhoodMove a)
{
    for(std::vector<NeighborhoodMove>::iterator it = tabuTable.begin(); it != tabuTable.end(); ++it) {
        if(*it == a){
            return false;
        }
    }
    return true;
}

void emili::pfsp::PfspTabuInsertNeighborhood::updateTabuTable(NeighborhoodMove a){
    tabuTable.insert(tabuTable.begin(),a);
    int size = tabuTable.size();
    if(size > tabutenure){
        tabuTable.erase(tabuTable.end()-1);
    }
}

emili::pfsp::PermutationFlowShopSolution* emili::pfsp::PfspTabuInsertNeighborhood::computeStep(std::vector<int> &solution,double valuec)
{
    std::vector < int > newsol;
    long int value = valuec;
    std::vector < int > best = solution;
    for(int i = 1; i< njobs; i++)
    {
        int sol_i = solution[i];
        for(int j = 1; j<=njobs; j++)
        {
            newsol = solution;
            newsol.erase(newsol.begin()+i);
            newsol.insert(newsol.begin()+j,sol_i);
            long int new_value = instance.computeWT(newsol);
            if(new_value < value){
                NeighborhoodMove move(i,j);
                if(notTabu(move)){
                    best = newsol;
                    value = new_value;
                    updateTabuTable(move);
                }
            }

        }
    }
    return new emili::pfsp::PermutationFlowShopSolution(value,best);
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
