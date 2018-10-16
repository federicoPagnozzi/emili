#include "mbo.h"

void inline restart_bird(emili::Solution* bird,emili::Neighborhood* restarter)
{
    emili::Solution* res = restarter->random(bird);
    *bird = *res;
    delete res;
}

inline void update_wing(int x,int k,int i,
                        std::vector<emili::Solution*>& nghbs,
                        emili::Solution* lefb,
                        emili::Perturbation* pert,
                        emili::Neighborhood* neigh,
                        std::vector<int>& wing_age,
                        int age)
{
    for(int o=0;o<x;o++)
    {
        emili::Solution* s = pert->perturb(lefb);
        nghbs.push_back(s);
    }
    // sort by fitness
    std::sort(nghbs.begin(),nghbs.end(),[](emili::Solution* s1, emili::Solution* s2)
              {
              return *s1 < *s2;
              });
    // check if bird has to be updated
    if(*lefb > *nghbs[0])
    {
        *lefb = *nghbs[0];
        wing_age[i] = 0;
        delete nghbs[0];
        nghbs.erase(nghbs.begin());
    }
    else
    {
        //update age
        //TODO check if age is over limit
       wing_age[i] ++;
       if(wing_age[i] >= age)
       {
           restart_bird(lefb,neigh);
           wing_age[i] = 0;
       }
    }
    // delete last x neighbors
    for(int o=(k-x);o<(k-1);o++)
    {
        delete nghbs[o];
        nghbs.erase(nghbs.begin()+o);
    }
}

inline void change_leader(emili::Solution* leader,
                          int& leader_age,
                          std::vector< emili::Solution* >& wing,
                          std::vector< int >& wing_age,int i =0)
{
    //swap leader pointer
    emili::Solution* l = leader;
    leader = wing[i];
    wing.erase(wing.begin()+i);
    wing.push_back(l);
    //swap solutions age
    int ola = leader_age;
    leader_age = wing_age[i];
    wing_age.erase(wing_age.begin()+i);
    wing_age.push_back(ola);
}

inline void prand_leader_update(emili::Solution* leader,
                                int& leader_age,
                                std::vector< emili::Solution* >& wing,
                                std::vector< int >& wing_age,
                                int size)
{
    //calculate sum of inverse of solution's age
    float den = 0;
    std::vector< int > psol;
    for(int o = 0 ; o < size; o++)
    {
        den += 1.0/(float)wing_age[o];
    }
    // select new leader solution based on
    // probability and age
    float prob = emili::generateRealRandomNumber();
    float level = 0;
    for(int o=0 ; o < size; o++)
    {
        level += (1.0/(float)wing_age[o])/den;
        if(prob < level)
        {
            //Update leader and exit
            change_leader(leader,leader_age,wing,wing_age,o);
            break;
        }
     }
}

emili::Solution* emili::pop::EMBO::search(emili::Solution* initial)
{
    //intialize pop    
    //initialize leftwing and rightwing
    for(int i =0; i< pop_size ; i++)
    {
        emili::Solution* s = init->generateSolution();
        if(i%2 == 0)
            leftwing.push_back(s);
        else
            rightwing.push_back(s);
    }
    int lwsize = leftwing.size();
    int rwsize = rightwing.size();
    this->leader = init->generateSolution();

    do
    {
       for(int j=0; j<m; j++)
       {
           //improve leader
           emili::Solution* nl = ls->search(leader);
           *leader = *nl;
           delete nl;
           //generate k neighbors
           std::vector< emili::Solution* > nghbs;
           emili::Solution* bestN=pert->perturb(leader);
           nghbs.push_back(bestN);
           for(int o=0;o<k;o++)
           {
               //Check tabulist
                emili::Solution* s = pert->perturb(leader);
                nghbs.push_back(s);
           }
           //order neighbors
           std::sort(nghbs.begin(),nghbs.end(),[](emili::Solution* s1, emili::Solution* s2)
                     {
                     return *s1 < *s2;
                     });
           //update leader or leader age
           if(*leader > *nghbs[0])
           {
                *leader = *nghbs[0];
                leader_age = 0;
                delete nghbs[0];
                nghbs.erase(nghbs.begin());
           }
           else
           {

               leader_age ++;
               //check if leader age is too big
               if(leader_age >= age)
               {
                   restart_bird(leader,neighbh);
                   leader_age = 0;
               }

           }
           std::vector<Solution*> lnghbs;
           std::vector<Solution*> rnghbs;
           for(int o=0;o<(k-x);o++)
           {
                lnghbs[o] = nghbs[o];
                rnghbs[o] = nghbs[o];
           }

           for(int i=0;i<lwsize;i++)
           {
               //integrate k-x best neighbors               
               //improve solutions by generating x neighbors
               emili::Solution* lefb = leftwing[i];
               update_wing(x,k,i,lnghbs,lefb,pert,neighbh,leftwing_age,age);

           }
           //clean lnghbs
           for(int i=0; i < lnghbs.size(); i++)
               delete lnghbs[i];

           for(int i=0;i<rwsize;i++)
           {
               //integrate k-x best neighbors
               //improve solutions by generating x neighbors
               emili::Solution* lefb = rightwing[i];
               update_wing(x,k,i,rnghbs,lefb,pert,neighbh,rightwing_age,age);
           }
           //clean rnghbs
           for(int i=0; i< rnghbs.size(); i++)
               delete rnghbs[i];
       }
       //Update best solution
       if(*leader < *bestSoFar)
       {
           *bestSoFar = *leader;
       }
       //Select new leader and move current leader on the back
       //
       float q = emili::generateRealRandomNumber();

           int t = emili::generateRandomNumber()%2;
           if(t)
           {
              if(q <= q0)
              {
                change_leader(leader,leader_age,leftwing,leftwing_age);
              }
              else
              {
                prand_leader_update(leader,leader_age,leftwing,leftwing_age,lwsize);
              }
           }
           else
           {
              if(q <= q0)
              {
               change_leader(leader,leader_age,rightwing,rightwing_age);
              }
              else
              {
                prand_leader_update(leader,leader_age,leftwing,leftwing_age,rwsize);
              }
           }

    }while(!termcriterion->terminate(initial,leader));
    return bestSoFar;
}
