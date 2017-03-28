//
//  main.cpp
//  OrCode
//
//  Created by garazi on 25/01/2017.
//  Copyright © 2017 Garazi. All rights reserved.
//

#include<iostream>
#include<fstream>
#include <sstream>
#include<string>
#include<cfloat>
#include<cmath>
#include <stdio.h>
#include <stdlib.h>
#include <cstring>
#include <ctime>
#include <limits.h>

#include<time.h>
#include <iomanip>
#include <cstdlib>

#include "openfiles.hpp"
#include "Instance.hpp"
#include "actualizar.hpp"
#include "matrix.hpp"
#include "InitialSolution.hpp"

#include "Neighborhood.hpp"
#include "Veicoli.hpp"
#include "RichiesteServizio.hpp"
//#define seed 6611
#define CONT 0
#define DEBUG 0
#define WRITE 0
//THINGS TO CHANGE

//eightstep evaluation scheme
#define EIGHT 1
//adjustment procedures
#define ADJ 1





int vrp() {
    
    //cout << "VARIABLE NEIGHBORHOOD SEARCH" << endl;
    string names[24];
    int zee;
  
    
    string instanc, instance;
    
    
   
    float T_max, T_red, n_imp;
    int Max_It;
    T_max=1.2;
    T_red=300;
    n_imp=400;
    
    string b;
    string c="/Users/garazi/Documents/Garazi/Data/new_instances/file.txt";
    ifstream ar11(c.c_str());
    
	   while (!ar11.eof()) {
           for(zee=0; zee<24; zee++){
               ar11>> names[zee];
               
           }
       }
	   ar11.close();
    
    
    //cout << "We are using the file " << b << endl;
    
    
    for(zee=1; zee<4; zee++){
        b=names[zee];
        std::cout << "File " << b << std::endl;
        bool solamtrovata;
        
        int ze;
        
        
        Instance* inst=NULL;
        
        
        inst=new Instance();
        
        inst->read_instance(b, zee);
        
      
        
        int numVeicoli=inst->numVeicoli0;
        
        int numRichieste=inst->numRichieste0;
        
        int numLocation=inst->numLocation0;
        
        
        string d="/Users/garazi/Documents/Garazi/Data/new_instances/seeds.txt";
        
        int SEEDS[5];
        ifstream ar12(d.c_str());
        
        while (!ar12.eof()) {
            for(ze=0; ze<5; ze++){
                ar12>> SEEDS[ze];
            }
        }
        ar12.close();
        
        
        
        inst->TimeWindowTightening();
       
        
        
        std::vector<std::vector<double>>& D=inst->Dist;
        std::vector<std::vector<double>>& MatTemp=inst->T;
        std::vector<RichiesteServizio*>& ric = inst->rc;
        std::vector<Veicoli*>& veic = inst->vec;
        std::vector<std::vector<int>>& MatCompVei=inst->MCV;
    
        int i;
        for(i=0; i<numRichieste; i++){
            ric[i]->id=i;
        }
        for(i=0; i<numRichieste; i++){
            ric[i]->display_ric();
        }
        
        for(i=0; i<numRichieste; i++){
            ric[i]->Barellato_SediaRotelle=0;
        }
        for(i=0; i<numVeicoli+1; i++){
            veic[i]->id=i;
        }

  
        
        
        for(ze=0; ze<5; ze++){
            solamtrovata=false;
            clock_t start, end, end1, st1, st2, en2;
            start=clock();
            
            
            int I;

            
            int seed=SEEDS[ze];
            srand(seed);
            
            
            ostringstream os;
            string name;


            st1=clock();
            SolutionVRP* Sol=NULL;
            
            Sol=new SolutionVRP();

            
            Sol->add_num_empty_route(numVeicoli);
            Sol->Initialize(veic, numVeicoli, MatTemp);
            Sol->setSolutionValue(DBL_MAX);
            Sol->numAddRoutes=INT_MAX;
            
            
            
            st2=clock();
            Sol=InitialSolutionBraekersF2(Sol, numVeicoli,  numRichieste,  ric,  veic, D, MatCompVei,  MatTemp);
            std:: cout << "Best Initial Solution is :\n";
            
            std::cout << "costF: " << Sol->getSolutionValue() << std::endl;

           
            en2=clock();
            double temp1=0;
            temp1=((double)(en2)-(st2))/CLOCKS_PER_SEC;
            std::cout << "Tempo di calcolo : " << temp1 << std::endl;
            Sol->DisplaySolution();
  
           // Sol->DisplaySolution();
            
            
           // double tempo=0;
            //tempo=((double)(end1)-(start))/CLOCKS_PER_SEC;
           // std::cout << "Tempo di calcolo : " << tempo << std::endl;
            if(Sol->numAddRoutes==0){
                
                std::cout << "Feasible" << std::endl;
            }
           // Sol->DisplaySolution();
            
            SolutionVRP* Sol_Curr=NULL;
            
            Sol_Curr=new SolutionVRP();
            
            Sol_Curr->CopySolution(Sol);
            
            //return 0;
            
            
            SolutionVRP* Sol_Best=NULL;
            Sol_Best=new SolutionVRP();
            Sol_Best->CopySolution(Sol);
            
           
            //Parameter initialization.
            float T=0;
            float r;
            
            int MaxNumIterations=350000;
            int imp=0;
            int numofoperator=0;
            
            
            for(i=0;i<MaxNumIterations; i++){
                
                //std::cout << i << std::endl;
                if(i%10000==0)
                {
                    std::cout << i << std::endl;
                }
                imp+=1;
                if(Sol_Curr->numAddRoutes>0){
                    numofoperator=5;
                }else{
                    numofoperator=4;
                }
                //Apply each operator in the order choosen.
                int* OR;
                OR=new int[numofoperator];
                for(int j=0;j<numofoperator;j++){
                    //apply the neighborhood solution.
                    OR[j]=j;
                    //acceptance criteria is met
                    //update imp
                }
                OR=random_order2(numofoperator,OR);
                if(i%10000==0)
                {
                   // for(int j=0; j<numofoperator;j++){
                     //   std::cout << OR[j] << " " ;
                    //}std::cout << "\n";
                }
                for(int j=0;j<numofoperator;j++){
                    //apply the neighborhood solution.
                    if(j==0){
                        //RELOCATE
                       // std::cout << "relI" << std::endl;
                       Sol_Curr=Relocate_NeighborhoodF(Sol_Curr, numVeicoli, MatCompVei, ric, MatTemp,  veic, D, numRichieste);
                        //Sol_Curr->DisplaySolution();
                      //  std::cout << "relF" << std::endl;
                        if((Sol->numRoutes+Sol->numAddRoutes>numVeicoli && Sol_Curr->numRoutes+Sol_Curr->numAddRoutes<Sol->numRoutes+Sol->numAddRoutes) || Sol_Curr->getSolutionValue()<Sol->getSolutionValue()+T){
                            Sol->CopySolution(Sol_Curr);
                            if((Sol_Best->numRoutes+Sol_Best->numAddRoutes>numVeicoli && Sol->numRoutes+Sol->numAddRoutes<Sol_Best->numRoutes+Sol_Best->numAddRoutes) || Sol->getSolutionValue()<Sol_Best->getSolutionValue()){
                                
                                Sol_Best->CopySolution(Sol);
                                imp=0;
                                //std::cout <<"changed" <<std::endl;
                                std::cout << "New Cost: " << Sol_Best->getSolutionValue() << " Relocate " << "It: " << i << std::endl;
                            }
                            
                        }else{
                            Sol_Curr->CopySolution(Sol);
                        }
                    }else{
                        if(j==1){
                            //EXCHANGE
                            
                        }else{if(j==2){
                            //2-OPT
                           // std::cout << "2optI" << std::endl;
                            Sol_Curr=two_opt( Sol_Curr,  numVeicoli,  D, veic, ric, numRichieste);
                          //  std::cout << "2optF" << std::endl;
                            if((Sol->numRoutes+Sol->numAddRoutes>numVeicoli && Sol_Curr->numRoutes+Sol_Curr->numAddRoutes<Sol->numRoutes+Sol->numAddRoutes) || Sol_Curr->getSolutionValue()<Sol->getSolutionValue()+T){
                                Sol->CopySolution(Sol_Curr);
                                if((Sol_Best->numRoutes+Sol_Best->numAddRoutes>numVeicoli && Sol->numRoutes+Sol->numAddRoutes<Sol_Best->numRoutes+Sol_Best->numAddRoutes) || Sol->getSolutionValue()<Sol_Best->getSolutionValue()){
                                    
                                    Sol_Best->CopySolution(Sol);
                                    imp=0;
                                    //std::cout <<"changed" <<std::endl;
                                    std::cout << "New Cost: " << Sol_Best->getSolutionValue() <<  " 2-opt " << " It: " << i << std::endl;
                                }
                                
                            }else{
                                Sol_Curr->CopySolution(Sol);
                            }
                            
                        }else{if(j==3){
                            //R-4-OPT
                            //std::cout << "4optI" << std::endl;
                            Sol_Curr=r_4_opt(Sol_Curr,  numVeicoli, D, ric,  veic);
                          //  std::cout << "4optF" << std::endl;
                            if((Sol->numRoutes+Sol->numAddRoutes>numVeicoli && Sol_Curr->numRoutes+Sol_Curr->numAddRoutes<Sol->numRoutes+Sol->numAddRoutes) || Sol_Curr->getSolutionValue()<Sol->getSolutionValue()+T){
                                Sol->CopySolution(Sol_Curr);
                                if((Sol_Best->numRoutes+Sol_Best->numAddRoutes>numVeicoli && Sol->numRoutes+Sol->numAddRoutes<Sol_Best->numRoutes+Sol_Best->numAddRoutes) || Sol->getSolutionValue()<Sol_Best->getSolutionValue()){
                                    
                                    Sol_Best->CopySolution(Sol);
                                    imp=0;
                                    //std::cout <<"changed" <<std::endl;
                                    std::cout << "New Cost: " << Sol_Best->getSolutionValue() << " r-4-opt " << "It: " << i << std::endl;
                                }
                                
                            }else{
                                Sol_Curr->CopySolution(Sol);
                            }
                            
                            
                            
                        }else{
                            //ELIMINATE
                            Sol_Curr=Eliminate_NeighborhoodF(Sol_Curr,  numVeicoli,  MatCompVei,  ric,  MatTemp,  veic,  D, numRichieste);
                            if((Sol->numRoutes+Sol->numAddRoutes>numVeicoli && Sol_Curr->numRoutes+Sol_Curr->numAddRoutes<Sol->numRoutes+Sol->numAddRoutes) || Sol_Curr->getSolutionValue()<Sol->getSolutionValue()+T){
                                Sol->CopySolution(Sol_Curr);
                                if((Sol_Best->numRoutes+Sol_Best->numAddRoutes>numVeicoli && Sol->numRoutes+Sol->numAddRoutes<Sol_Best->numRoutes+Sol_Best->numAddRoutes) || Sol->getSolutionValue()<Sol_Best->getSolutionValue()){
                                    
                                    Sol_Best->CopySolution(Sol);
                                    imp=0;
                                    //std::cout <<"changed" <<std::endl;
                                    std::cout << "New Cost: " << Sol_Best->getSolutionValue() << " Eliminate " << "It: " << i << std::endl;
                                }
                                
                            }else{
                                Sol_Curr->CopySolution(Sol);
                            }
                            
                        }}}}
                    
                }
                
                
                //acceptance criteria is met
                
                //update imp
                if(imp>0){
                    T=T-T_max/T_red;
                    if(T<0){
                        r=((float) rand())/RAND_MAX;
                        T=r*T_max;
                        
                        if(imp>n_imp*Sol_Best->get_num_routes()){
                            
                            if(Sol_Best->getSolutionValue()!=Sol_Curr->getSolutionValue()){
                                Sol->CopySolution(Sol_Best);
                                Sol_Curr->CopySolution(Sol_Best);
                            }
                            imp=0;
                        }
                        
                    }
                }
                
                
                
                
                
                
                
                delete[] OR;
                //un parametro que lea en un file... tambien tengo que pasar por un file 360 o 420 o asi.
                
            }//END OF FOR
            
            
            end=clock();
            double tempo1=0;
            tempo1=((double)(end)-(start))/CLOCKS_PER_SEC;
            std::cout<<"Solution: " << Sol_Best->getSolutionValue() << std::endl;
            
           // Sol_Best->DisplaySolution();
            if(Sol_Best->numAddRoutes==0){
                
                std::cout << "Feasible" << std::endl;
            }
            std::cout << "Tempo di calcolo : " << tempo1 << std::endl;
            
            delete Sol;
            
            delete Sol_Best;
            
            
            
            
            delete Sol_Curr;
            
        }
        
        
    }
    //system("pause");
    
    //std::cout << "fin de verdad" << endl;
    
    return 0;
    
}

