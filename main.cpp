#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <ctime>
#include <cstring>
#include "permutationflowshop.h"
#include "paramsparser.h"
#define FIRST_IMPROVEMENT 1
#define BEST_IMPROVEMENT 2
#define RANDOM_IS 1
#define SLACK_IS 2
#define EXCHANGE 3
#define TRANSPOSE 5
#define INSERT 7
#define VND1 11
#define VND2 13
#define TABU 17


using namespace std;

int mainss()
{

    srand(time(0));
    std::vector< int > vv;
    char path[] = "/Users/federicopagnozzi/sviluppo/QTCREAT/exercise1-build/instances/70x20_1";
    emili::pfsp::PermutationFlowShop problem(path);  
    emili::InitialSolution* rand_solution = new emili::pfsp::PfspSlackInitialSolution(problem);
    emili::pfsp::PfspNeighborhood* een = new emili::pfsp::PfspInsertNeighborhood(problem);

    //emili::Termination* ptc = new emili::pfsp::PfspTerminationClassic(problem);

    emili::Termination* ptc = new emili::pfsp::PfspTerminationIterations(50);
    emili::pfsp::PfspTabuValueMemory prh(10);
    //emili::pfsp::PfspMovesMemory prh(10,*een);
    //emili::BestImprovementSearch(rand_solution,ptc,een);
    emili::LocalSearch* ls = new emili::TabuSearch(*rand_solution,*ptc,*een,prh);

    clock_t time = clock();

    emili::Solution* sol_d = ls->search();
    clock_t end = clock();
    std::cout << "Local search required -> " << ((end-time)/(float)CLOCKS_PER_SEC) << std::endl;
    vv = *((std::vector<int>*)(*sol_d).getRawData());
    for (int i = 1; i <= problem.getNjobs(); ++i)
    cout << vv[i] << " " ;
    cout << endl;
    cout << "Local Search Slack BestImprov Insert -> ";
    cout << problem.evaluateSolution(*sol_d)<<std::endl;
    een->reset();

    //emili::pfsp::PfspTerminationIterations pft(problem,10);
    emili::WhileTrueTermination pft;
    emili::pfsp::PfspRandomSwapPertub perturb(problem);
    emili::pfsp::PfspTestAcceptance acc(problem);
    emili::IteratedLocalSearch ils(*ls,pft,perturb,acc);
    //emili::Solution* sol_ils = ils.search();
    emili::Solution* sol_ils = ils.timedSearch(15);

    vv = *((std::vector<int>*)(*sol_ils).getRawData());
    for (int i = 1; i <= problem.getNjobs(); ++i)
    cout << vv[i] << " " ;
    cout << endl;
    cout << "Iterated Local search results -> ";
    cout << problem.evaluateSolution(*sol_ils)<<std::endl;

    //emili::pfsp::PermutationFlowShopSolution p(vv);
    cout << "Hello World!" << endl;
    return 0;
}

/***********************************************************************/

void usage(char* progName)
{
    cout << "Usage : " << progName << " <instance_file> <first | best | tabuv | tabuh | tabum> <random | slack> <exchange | transpose | insert > [ils] [time] [tabu_tenure] [tabu_iterations]" << endl;
   // cout << "        vnd1: transpose -> exchange -> insert" <<endl;
    //cout << "        vnd2: transpose -> insert   -> exchange" <<endl;
    exit(0);
}



int main(int argc, char *argv[])
{
    /* initialize random seed: */
    srand ( time(0) );

    /* Create instance object */
    PfspInstance instance;

    /* Read data from file */
    if (! instance.readDataFromFile(argv[1]) )
      return 1;

    emili::pfsp::PermutationFlowShop problem(instance);
    prs::ParamsParser ps(argv,argc,problem);
    emili::LocalSearch* ls = ps.parseParams();    
    std::cout << "searching..." << std::endl;
    clock_t time = clock();
    emili::Solution* solution;
    if(ps.ils_time>0)
    {
       solution = ls->timedSearch(ps.ils_time);
    }
    else
    {

        solution = ls->search();
    }
    double time_elapsed = (double)(clock()-time)/CLOCKS_PER_SEC;
    cout << "time : " << time_elapsed << std::endl;
    std::vector < int > *sol = (std::vector < int > *)solution->getRawData();
    cerr << time_elapsed << " ";
    cout << "Found solution: ";
    for (int i = 1; i <= instance.getNbJob(); ++i)
      cout << (*sol)[i] << " " ;
    cout << endl;
    long int totalWeightedTardiness = instance.computeWT(*sol);
    cout << "Total weighted tardiness: " << totalWeightedTardiness << endl;
    cerr << totalWeightedTardiness << endl;
    exit(0);
}

int mainold(int argc, char *argv[])
{

    emili::pfsp::PfspNeighborhood* n;
    emili::InitialSolution* is;
    emili::LocalSearch* ls;
    emili::Termination* ptc;
    int timespan;
    cout << "Parameters: \n";
    int init = 0;
    int search = 0;
    int neigh = 0;
    int vnd = 0;
    if (argc < 5)
    {
      usage(argv[0]);
    }
    else
    {
        /* initialize random seed: */
        srand ( time(0) );

        /* Create instance object */
        PfspInstance instance;

        /* Read data from file */
        if (! instance.readDataFromFile(argv[1]) )
          return 1;

        emili::pfsp::PermutationFlowShop problem(instance);

        if(strcmp(argv[3],"random") == 0)
        {
            init = RANDOM_IS;
             is = new emili::pfsp::PfspRandomInitialSolution(problem);

            cout << "\tINITIAL SOLUTION RANDOM "<<endl;
        }
        else if( strcmp(argv[3],"slack") == 0)
        {
            init = SLACK_IS;
            is = new emili::pfsp::PfspSlackInitialSolution(problem);

            cout << "\tINITIAL SOLUTION SLACK "<<endl;
        }
        else
        {
            usage(argv[0]);
        }

        if(strcmp(argv[4],"exchange") == 0)
        {
            neigh = EXCHANGE;

                 n = new emili::pfsp::PfspExchangeNeighborhood(problem);


            cout << "\tEXCHANGE NEIGHBORHOOD "<<endl;
        }
        else if(strcmp(argv[4],"transpose") == 0)
        {
            neigh = TRANSPOSE;

                n = new emili::pfsp::PfspTransposeNeighborhood(problem);


            cout << "\tTRANSPOSE NEIGHBORHOOD "<<endl;
        }
        else if(strcmp(argv[4],"insert") == 0 )
        {
            neigh = INSERT;

                n = new emili::pfsp::PfspInsertNeighborhood(problem);


            cout << "\tINSERT NEIGHBORHOOD "<<endl;
        }

        int iterations = 50;
        int tabutenure = 10;

        if(argc >= 8){
            tabutenure = atoi(argv[7]);
        }

        if(argc >= 9){
            iterations = atoi(argv[8]);
        }

        if(strcmp(argv[2],"first") == 0)
        {
            search = FIRST_IMPROVEMENT;
            emili::pfsp::PfspTerminationClassic ptc;
            ls = new emili::FirstImprovementSearch(*is,ptc,*n);
            cout << "\tFIRST IMPROVEMENT "<<endl;
        }
        else if(strcmp(argv[2],"best") == 0)
        {
            search = BEST_IMPROVEMENT;
            emili::pfsp::PfspTerminationClassic ptc;
            ls = new emili::BestImprovementSearch(*is,ptc,*n);
            cout << "\tBEST IMPROVEMENT "<< search<<endl;
        }
        else if(strcmp(argv[2],"tabum") == 0)
        {
            search = TABU;
            emili::pfsp::PfspTerminationIterations ptc(iterations);
            emili::pfsp::PfspMovesMemory mm(tabutenure,*n);
            ls = new emili::TabuSearch(*is,ptc,*n,mm);
            cout << "\tTabu moves memory\n\ttabuTenure: "<< tabutenure<< "\n\titerations: "<< iterations<<endl;
        }
        else if(strcmp(argv[2],"tabuh") == 0)
        {
            search = TABU;
            emili::pfsp::PfspTerminationIterations ptc(iterations);
            emili::pfsp::PfspTabuHashMemory mm(tabutenure);
            ls = new emili::TabuSearch(*is,ptc,*n,mm);
            cout << "\tTabu hash memory\n\ttabuTenure: "<< tabutenure<< "\n\titerations: "<< iterations<<endl;
        }
        else if(strcmp(argv[2],"tabuv") == 0)
        {
            search = TABU;
            emili::pfsp::PfspTerminationIterations ptc(iterations);
            emili::pfsp::PfspTabuValueMemory mm(tabutenure);
            ls = new emili::TabuSearch(*is,ptc,*n,mm);
            cout << "\tTabu fullsolution memory\n\ttabuTenure: "<< tabutenure<< "\n\titerations: "<< iterations<<endl;
        }
        else if(strcmp(argv[2],"vnd1") == 0 )
                {
                    neigh = VND1;
                    emili::pfsp::PfspTerminationClassic ptc;
                    n = new emili::pfsp::PfspExchangeNeighborhood(problem);
                    emili::pfsp::PfspTransposeNeighborhood n2(problem);
                    emili::pfsp::PfspInsertNeighborhood n3(problem);
                    std::vector< emili::Neighborhood*> nn(3);
                    nn[0] = &n2;
                    nn[1] = &n3;
                    nn[2] = n;
                    ls = new emili::VNDSearch< emili::FirstImprovementSearch >(*is,ptc,nn);
                    cout << "\tVND1 NEIGHBORHOOD "<<endl;
                }
                else if(strcmp(argv[2],"vnd2") == 0 )
                {
                    neigh = VND2;
                    emili::pfsp::PfspTerminationClassic ptc;
                    n = new emili::pfsp::PfspExchangeNeighborhood(problem);
                    emili::pfsp::PfspTransposeNeighborhood n2(problem);
                    emili::pfsp::PfspInsertNeighborhood n3(problem);
                    std::vector< emili::Neighborhood*> nn(3);
                    nn[0] = &n2;
                    nn[1] = n;
                    nn[2] = &n3;
                    ls = new emili::VNDSearch< emili::BestImprovementSearch >(*is,ptc,nn);
                    //ls = new emili::pfsp::VNDFirstSearch(*is,ptc,n2,n3,*n);
                    cout << "\t neighborhood vector size : " << nn.size() << std::endl;
                    cout << "\tVND2 NEIGHBORHOOD "<<endl;
                }

        else
        {
            usage(argv[0]);
        }

        /*else if(strcmp(argv[4],"vnd1") == 0 )
        {
            neigh = VND1;
            vnd = 1;
            cout << "\tVND1 NEIGHBORHOOD "<<endl;
        }
        else if(strcmp(argv[4],"vnd2") == 0 )
        {
            neigh = VND2;
            vnd = 2;
            cout << "\tVND2 NEIGHBORHOOD "<<endl;
        }
        else
        {
            usage(argv[0]);
        }
        */






    std::cout << "searching..." << std::endl;
    clock_t time = clock();    
    emili::Solution* solution = ls->search();
    double time_elapsed = (double)(clock()-time)/CLOCKS_PER_SEC;
    cout << "time : " << time_elapsed << std::endl;
    std::vector < int > *sol = (std::vector < int > *)solution->getRawData();
    cerr << time_elapsed << " ";
    cout << "Found solution: ";
    for (int i = 1; i <= instance.getNbJob(); ++i)
      cout << (*sol)[i] << " " ;
    cout << endl;
    long int totalWeightedTardiness = instance.computeWT(*sol);
    cout << "Total weighted tardiness: " << totalWeightedTardiness << endl;
    cerr << totalWeightedTardiness << endl;
    if(argc >= 6 && strcmp(argv[5],"ils")== 0)
    {
        int time = 10;
        if(argc >= 7){
            time = atoi(argv[6]);
        }
    emili::WhileTrueTermination pft;
    emili::pfsp::PfspRandomSwapPertub perturb(problem);
    emili::pfsp::PfspTestAcceptance acc(problem);
    emili::IteratedLocalSearch ils(*ls,pft,perturb,acc);
    //emili::Solution* sol_ils = ils.search();
    emili::Solution* sol_ilsa = ils.timedSearch(time);
    std::vector < int >* sol_ils = (std::vector < int > *)sol_ilsa->getRawData();
    cout << "Found solution: ";
    for (int i = 1; i <= instance.getNbJob(); ++i)
      cout << (*sol_ils)[i] << " " ;
    cout << endl;
    totalWeightedTardiness = instance.computeWT(*sol_ils);
    cout << "Total weighted tardiness: " << totalWeightedTardiness << endl;
    cerr << totalWeightedTardiness << endl;
    }
    delete n;

    delete is;

    delete solution;

     }
    return 0;
}

void testP()
{

}


