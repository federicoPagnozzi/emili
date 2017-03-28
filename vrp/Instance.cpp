//
//  Instance.cpp
//  C-B
//
//  Created by garazi on 10/03/2017.
//  Copyright © 2017 Garazi. All rights reserved.
//

#include "Instance.hpp"
#include<fstream>
#include <sstream>

Instance::Instance() {

   
   
    
    
}

Instance::~Instance() {

    
    
}


void Instance::read_Distance(std::string a){

    std::string c="/Users/garazi/Documents/Garazi/Data/new_instances/"+a+"/MatriceDistanze.txt";
    std::ifstream ar1(c.c_str());
    int i,j;
    std::vector <double> d(numLocation0+1);

    while (!ar1.eof()){
        for(i=0; i<numLocation0+1; i++){
            Dist.push_back(d);
            for(j=0; j<numLocation0+1; j++){
                ar1 >> Dist[i][j];
            }
        }
    }
    ar1.close();
}

void Instance::read_Time(std::string a){
    
    std::string c="/Users/garazi/Documents/Garazi/Data/new_instances/"+a+"/MatriceTempi.txt";
    std::ifstream ar1(c.c_str());
    int i,j;
    std::vector <double> d(numLocation0+1);

    while (!ar1.eof()){
        for(i=0; i<numLocation0+1; i++){
             T.push_back(d);
            for(j=0; j<numLocation0+1; j++){
                ar1 >> T[i][j];
            }
        }
    }
    ar1.close();
}

void Instance::read_MatCompVei(std::string a){
    
    std::string c="/Users/garazi/Documents/Garazi/Data/new_instances/"+a+"/MatriceCompVeicoli.txt";
    std::ifstream ar1(c.c_str());
    int i,j;
    std::vector <int> d(numVeicoli0);

    while (!ar1.eof()){
        for(i=0; i<numRichieste0; i++){
            MCV.push_back(d);
            for(j=0; j<numVeicoli0; j++){
                ar1 >> MCV[i][j];
            }
        }
    }
    ar1.close();
    
}

void Instance::read_ric(std::string a){
    std::string c="/Users/garazi/Documents/Garazi/Data/new_instances/"+a+"/RichiesteServizio.txt";
    std::ifstream ar1(c.c_str());
    int i;
    //int T;
   // double b;
    while (!ar1.eof()) {
        for(i=0; i<numRichieste0; i++){
            
            rc.push_back(new RichiesteServizio);
            
            ar1 >> rc[i]->id;
            ar1 >> rc[i]->tipo;
            ar1 >> rc[i]->oc;
            ar1 >> rc[i]->timewinPmin;
            ar1 >> rc[i]->timewinPmax;
            ar1 >> rc[i]->timewinDmin;
            ar1 >> rc[i]->timewinDmax;
            ar1 >> rc[i]->locationP;
            ar1 >> rc[i]->locationD;
            ar1 >> rc[i]->paziente;
            ar1 >> rc[i]->idSosta;
            ar1 >> rc[i]->Barellato_SediaRotelle;
            ar1 >> rc[i]->obbligo_associazione;
            ar1 >> rc[i]->staff;
            ar1 >> rc[i]->seated;
            ar1 >> rc[i]->stretcher;
            ar1 >> rc[i]->wheelchair;
            ar1 >> rc[i]->RideTime;
            ar1 >> rc[i]->st;
            
        }
    }
    ar1.close();
}
void Instance::read_veic(std::string a){
    std::string c="/Users/garazi/Documents/Garazi/Data/new_instances/"+a+"/Veicoli.txt";
    std::ifstream ar1(c.c_str());
    int i;

    while (!ar1.eof()) {
        for(i=0; i<numVeicoli0+1; i++){
            vec.push_back(new Veicoli);
            ar1 >> vec[i]->id;
            ar1 >> vec[i]->location;
            ar1 >> vec[i]->cap;
            ar1 >> vec[i]->disponibili;
            ar1 >> vec[i]->tipo;
            ar1 >> vec[i]->distanza;
            ar1 >> vec[i]->costo1;
            ar1 >> vec[i]->costo2;
            ar1 >> vec[i]->costoNumeroPazientiTrasportati;
            ar1 >> vec[i]->staff;
            ar1 >> vec[i]->seated;
            ar1 >> vec[i]->stretcher;
            ar1 >> vec[i]->wheelchair;
            ar1>> vec[i]->MaxRoute;

        }
    }
    ar1.close();
}

void Instance::read_num(std::string a){
    std::string c="/Users/garazi/Documents/Garazi/Data/new_instances/"+a+"/DimensioniIstanza.txt";
    std::ifstream ar1(c.c_str());
    
    ar1 >> numRichieste0;
    ar1 >> numVeicoli0;
    ar1 >> numLocation0;
    
    ar1.close();

}

void Instance::read_instance(std::string a, int i){
    read_num(a);
    read_Distance(a);
    read_Time(a);
    read_MatCompVei(a);
    read_ric(a);
    read_veic(a);
    
    id=i;
    name=a;
    
    
}

void Instance::TimeWindowTightening(){
    
    int i;
    for(i=0; i<numRichieste0;i++){

        if(rc[i]->tipo==1){
            if(rc[i]->timewinDmin-rc[i]->RideTime -rc[i]->st<0){
                rc[i]->timewinPmin=0;
            }else{
                rc[i]->timewinPmin=rc[i]->timewinDmin-rc[i]->RideTime -rc[i]->st;
            }
            
            
            
            if(rc[i]->timewinDmax-rc[i]->st-T[rc[i]->locationP][rc[i]->locationD]>1440){
                rc[i]->timewinPmax=1440;
            }else{
                rc[i]->timewinPmax=rc[i]->timewinDmax-rc[i]->st-T[rc[i]->locationP][rc[i]->locationD];
            }
        }
        if(rc[i]->tipo==2){
            if(rc[i]->timewinPmin+rc[i]->st+T[rc[i]->locationP][rc[i]->locationD]<0){
                rc[i]->timewinDmin=0;
            }else{
                rc[i]->timewinDmin=rc[i]->timewinPmin+rc[i]->st+T[rc[i]->locationP][rc[i]->locationD];
            }
            if(rc[i]->timewinPmax+rc[i]->st+rc[i]->RideTime>1440){
                rc[i]->timewinDmax=1440;
            }else{
                
                rc[i]->timewinDmax=rc[i]->timewinPmax+rc[i]->st+rc[i]->RideTime;
            }
            
            
            
        }
        
        
        
    }
    
    
    
}

double Instance::calcObjectiveFunctionValue(emili::Solution& solution)
{
     SolutionVRP& sol= (SolutionVRP&) solution;
     return sol.calculate_dist(Dist);
}
double Instance::evaluateSolution(emili::Solution & solution)
{
    double cost = calcObjectiveFunctionValue(solution);
    solution.setSolutionValue(cost);
    return cost;
}

int Instance::problemSize()
{
    return 1;
}
