//
//  InitialSolution.cpp
//  OrCode
//
//  Created by garazi on 25/01/2017.
//  Copyright © 2017 Garazi. All rights reserved.
//

#include "InitialSolution.hpp"
#include<iostream>
#include<fstream>
#include<string>
#include<cmath>
#include <cstring>
#include <ctime>
//#include <sys/time.h>
#include "openfiles.hpp"
#include "matrix.hpp"

#include "actualizar.hpp"

#include "Veicoli.hpp"
#include "RichiesteServizio.hpp"
#include "Route.hpp"



#include <float.h>

/*

Solution* InitialSolutionBraekers(int numVeicoli, int numRichieste, RichiesteServizio* ric, Veicoli* veic, double** D, int** MatCompVei, double** MatTemp, double MinAt, int MaxTimeRoute, double MaxRideTime){
    
    
    int i, j, l , g, uu;
    
    double bestdist;
    int bestpl, bestdl, bestv;
    bool inserted;
    int bestaddpl, bestadddl,bestaddv;
    
    Solution* InitialSol=NULL;
    
    InitialSol=new Solution();
    
    InitialSol->cost=DBL_MAX;
    Route* routad=NULL;
    routad=new Route();
    
    //0
    InitialSol->numAddRoutes=INT_MAX;
    
    Solution* Solaux=NULL;
    
    Route* rout=NULL;
    rout=new Route();
    
    //Ordering vector initialized
    int* E;
    E=new int[numRichieste];
    Solaux=new Solution();
    //The initial solution is calculated 1000 times and the best solution is chosen
    for(uu=0;uu<1000; uu++){
        
        Solaux->ClearSol();
        
        
        
        
        //Initialize the solution one route per vehicle
        Solaux->add_num_empty_route(numVeicoli);
        
        Solaux->Initialize(veic, numVeicoli, MatTemp, MaxTimeRoute);
        //std::cout << "New order" << std::endl;
        //E is a vector containing the order in which the requests will be inserted
        if(uu<1){
            E=timewinminP_order(ric, numRichieste,E);
        }else{
            E=random_order( ric,  numRichieste, E);
        }
        
               //Print order
        //std::cout << "ordering number " << uu << endl;
        //for(i=0;i<numRichieste;i++){
        //std::cout<< E[i] << " " ;
        
        //}std::cout << "\n";
        //Each request will be inserted in all the vehicles. Cost will be kept
        for(i=0;i<numRichieste;i++){
            //	std::cout<< "request-> "<< E[i] << "\n";
            inserted=false;
            //Create a Route to insert everything on it.
            
            bool feas1, feas2, feas3,feasride, feas4, feas5, feas6;
            //double Bestdisteffect=DBL_MAX;
            bestdist=DBL_MAX;
            for(j=0;j<numVeicoli;j++){
                //std::cout<< "vehicle-> "<< j << "\n";
                // See if the vehicle is compatible with the request to be inserted.
                
                
                rout->clear_route();
                
                
                
                Solaux->copyvehicleinroute(rout, j);
                //
                int originallength=rout->length;
                for(l=0; l<originallength-1; l++)
                {
                    //copy current Solution in route.
                    //Solaux->copyvehicleinroute(rout, j);
                    //Insert pickup in position l
                    //rout->display_route();
                    
                    rout->insert_pickup(E[i], l, ric);
                    
                    rout->calculate_earliest_latest(ric, MatTemp, MaxTimeRoute);
                    
                    rout->calculate_capacity(ric, veic);
                    
                    feas2=rout->check_feasibility_P_cap(E[i],l, veic,ric);
                    
                    //	std:: cout << "inserting pickup on " << l << "\n";
                    feas1=rout->check_feasibility_P_tw(E[i],l, ric, MatTemp);
                    
                    if(feas1>0 && feas2>0){
                        //		std::cout << "pickup of request " << E[i] << "compatible in position " << l << endl;
                        //
                        
                        for (g=l; g<originallength-1; g++){
                            
                            rout->insert_delivery(E[i], g, ric);
                            
                            
                            rout->calculate_capacity(ric, veic);
                            
                            rout->calculate_earliest_latest(ric, MatTemp, MaxTimeRoute);
                            
                            
                            feasride=rout->ridetime_check( g,  ric,  MaxRideTime);
                            
                            if(feasride>0){//std::cout <<"true" << endl;
                                double disteffect;
                                bool gg;
                                //	std:: cout << "length: " << rout->length << "l : " << l << "g : " << g << endl;
                                if(l==g){
                                    
                                    disteffect=0;
                                    disteffect=rout->effect_ondistance(E[i], l, g, D,  ric);
                                    //std:: cout << "length: " << rout->length << "l : " << l << "g : " << g << endl;
                                    //std::cout <<  disteffect << " < " << bestdist << "\n";
                                    if(disteffect<bestdist){
                                        
                                        feas3=rout->check_feasibillity_D_tw1(E[i],l,ric,MatTemp);
                                        //std::cout << feas3 << "  feas3 " << endl;
                                        
                                        if(feas3>0){
                                            //calculate 8 step evaluation scheme
                                            //	std::cout << "check delivery tw feasibility OK " << endl;
                                            
                                            gg=rout->Eightstepevaluationscheme(MatTemp, ric,  MaxRideTime,  MaxTimeRoute);
                                            //	rout->display_route();
                                            if(gg>0){
                                                bestdist=disteffect;
                                                bestpl=l;
                                                bestdl=g;
                                                bestv=j;
                                                inserted=true;
                                                
                                            }
                                        }
                                        
                                    }
                                }else{
                                    feas4=rout->check_feasibility_D_tw2( l,  g, ric, MatTemp);
                                    //	std::cout << feas3 << "  1feas3 " << endl;
                                    feas5=rout->check_feasibility_capacity(l,g,veic);
                                    //	std::cout << feas3 << "  2feas3 " << endl;
                                    //Insert E{i} in route j in position l pickup g delivery
                                    //rout.insert_req_pos(E[i], l, g, ric);
                                    if(feas4>0 && feas5>0){
                                        disteffect=rout->effect_ondistance(E[i], l, g, D,  ric);
                                        //std::cout <<  disteffect << " < " << MINI[0] << "\n";
                                        //	std:: cout << "length: " << rout->length << "l : " << l << "g : " << g << endl;
                                        //	std::cout <<  disteffect << " < " << bestdist << "\n";
                                        if(disteffect<bestdist){
                                            //calculate 8 step evaluation scheme
                                            feas6=rout->check_feasibillity_D_tw1(E[i],l,ric,MatTemp);
                                            //	std::cout << feas3 << "  feas3 " << endl;
                                            if(feas6>0){
                                                //calculate 8 step evaluation scheme
                                                //	std::cout << "check delivery tw feasibility OK " << endl;
                                                
                                                gg=rout->Eightstepevaluationscheme(MatTemp, ric,  MaxRideTime,  MaxTimeRoute);
                                                //	rout->display_route();
                                                if(gg>0){
                                                    bestdist=disteffect;
                                                    bestpl=l;
                                                    bestdl=g;
                                                    bestv=j;
                                                    inserted=true;
                                                    
                                                }
                                            }
                                            
                                        }else{
                                            //Continue with the following insertion, This position is not considered
                                        }
                                    }
                                }
                                
                            }
                            
                            rout->delete_delivery(g);
                            
                        }
                        
                        
                    }
                    rout->delete_pickup(l, MatTemp);
                    //std:: cout << "pickup deleted" << endl;
                    
                }
                
                
                
            }
            
            
            if(inserted>0){
                Solaux->route[bestv]->insert_req_pos(E[i],bestpl,bestdl,ric, MatTemp);
                Solaux->route[bestv]->calculate_capacity(ric,veic);
                
                //Solaux->route[MINI[3]]->display_route();
                Solaux->route[bestv]->calculate_earliest_latest(ric, MatTemp, MaxTimeRoute);
                Solaux->route[bestv]->Eightstepevaluationscheme(MatTemp,ric,MaxRideTime,MaxTimeRoute);
                //Solaux->copyrouteinsol(rout);
                //	std::cout << "request" << E[i] << "inserted in vehicle, " << bestv<< "\n";
                //Solaux->route[bestv]->display_route();
                Solaux->updatecost(numVeicoli);
                Solaux->updateusedvehicles(numVeicoli);
                
                
            }else{
                //std::cout<< "false" << std::endl;
                if(Solaux->numAddRoutes==0){
                    //Add a new one an insert it
                    //std::cout<< "false" << std::endl;
                    Solaux->addAdditionalRoute(numVeicoli, MatTemp, veic, MaxTimeRoute);
                    //std::cout<< numVeicoli+Solaux->numAddRoutes << std::endl;
                    Solaux->route[(numVeicoli+Solaux->numAddRoutes)-1]->insert_req_pos(E[i],0,0,ric, MatTemp);
                    //	std::cout << "request" <<  E[i] << "inserted in additional vehicle , " << (numVeicoli+Solaux->numAddRoutes)-1 << "\n";
                    //std::cout << "Add route added" << std::endl;
                    Solaux->route[(numVeicoli+Solaux->numAddRoutes)-1]->calculate_earliest_latest(ric, MatTemp, MaxTimeRoute);
                    Solaux->route[(numVeicoli+Solaux->numAddRoutes)-1]->Eightstepevaluationscheme(MatTemp,ric,MaxRideTime,MaxTimeRoute);
                    Solaux->updatecost(numVeicoli);
                    Solaux->updateusedvehicles(numVeicoli);
                    //	Solaux->route[(numVeicoli+Solaux->numAddRoutes)-1]->display_route();
                    
                }else{
                    
                    
                    double bestadddist=DBL_MAX;
                    bool addinsertion;
                    //Try to insert it in the ones that there are....
                    for(int t=0;t<Solaux->numAddRoutes; t++){
                        //try on each additional vehicle that exists
                        routad->clear_route();
                        Solaux->copyvehicleinroute(routad, t+Solaux->numRoutes);
                        //std::cout << t << std::endl;
                        addinsertion=false;
                        int originallength=routad->length;
                        for(l=0; l<originallength-1; l++)
                        {
                            
                            
                            routad->insert_pickup(E[i], l, ric);
                            routad->calculate_earliest_latest(ric, MatTemp, MaxTimeRoute);
                            
                            routad->calculate_capacity(ric, veic);
                            
                            feas2=routad->check_feasibility_P_cap(E[i],l, veic,ric);
                            
                            
                            
                            //std:: cout << "inserting pickup on " << l << "\n";
                            feas1=routad->check_feasibility_P_tw(E[i],l, ric, MatTemp);
                            
                            if(feas1>0 && feas2>0){
                                //std::cout << "pickup of request " << E[i] << "compatible in position " << l << endl;
                                //
                                
                                for (g=l; g<originallength-1; g++){
                                    
                                    routad->insert_delivery(E[i], g, ric);
                                    //		std::cout << "LENGTH " <<  rout->length << endl;
                                    //		rout->display_route();
                                    //routad->calculate_earliest_latest(ric, MatTemp, MaxTimeRoute);
                                    
                                    routad->calculate_capacity(ric, veic);
                                    //		rout->display_route();
                                    routad->calculate_earliest_latest(ric, MatTemp, MaxTimeRoute);
                                    //		std::cout << "EARLIEST and LATEST + length"  << rout->length << endl;
                                    //		rout->display_route();
                                    
                                    //checkridetime....
                                    
                                    feasride=routad->ridetime_check( g,  ric,  MaxRideTime);
                                    //	std::cout << "feasride " <<  feasride << endl;
                                    if(feasride>0){//std::cout <<"true" << endl;
                                        double distaddeffect;
                                        bool gg;
                                        //std:: cout << "length: " << rout->length << "l : " << l << "g : " << g << endl;
                                        if(l==g){
                                            
                                            distaddeffect=0;
                                            distaddeffect=routad->effect_ondistance(E[i], l, g, D,  ric);
                                            //std:: cout << "length: " << rout->length << "l : " << l << "g : " << g << endl;
                                            //std::cout <<  disteffect << " < " << bestdist << "\n";
                                            if(distaddeffect<bestadddist){
                                                
                                                feas3=routad->check_feasibillity_D_tw1(E[i],l,ric,MatTemp);
                                                //std::cout << feas3 << "  feas3 " << endl;
                                                //	routad->display_route();
                                                //std:: cout << ric[E[i]].timewinPmin << " " << ric[E[i]].timewinPmax << " " << ric[E[i]].timewinDmin << " " << ric[E[i]].timewinDmax << std::endl;
                                                if(feas3>0){
                                                    //calculate 8 step evaluation scheme
                                                    //	std::cout << "check delivery tw feasibility OK " << endl;
                                                    
                                                    gg=routad->Eightstepevaluationscheme(MatTemp, ric,  MaxRideTime,  MaxTimeRoute);
                                                    //	rout->display_route();
                                                    if(gg>0){
                                                        bestadddist=distaddeffect;
                                                        bestaddpl=l;
                                                        bestadddl=g;
                                                        bestaddv=t+Solaux->numRoutes;
                                                        addinsertion=true;
                                                        
                                                    }
                                                }
                                                
                                            }
                                        }else{
                                            feas4=routad->check_feasibility_D_tw2( l,  g, ric, MatTemp);
                                            //	std::cout << feas4 << "  1feas3 " << endl;
                                            feas5=routad->check_feasibility_capacity(l,g,veic);
                                            //	std::cout << feas5 << "  2feas3 " << endl;
                                            //Insert E{i} in route j in position l pickup g delivery
                                            //rout.insert_req_pos(E[i], l, g, ric);
                                            if(feas4>0 && feas5>0){
                                                distaddeffect=routad->effect_ondistance(E[i], l, g, D,  ric);
                                                //std::cout <<  disteffect << " < " << MINI[0] << "\n";
                                                //	std:: cout << "length: " << rout->length << "l : " << l << "g : " << g << endl;
                                                //	std::cout <<  disteffect << " < " << bestdist << "\n";
                                                if(distaddeffect<bestadddist){
                                                    //calculate 8 step evaluation scheme
                                                    feas6=routad->check_feasibillity_D_tw1(E[i],l,ric,MatTemp);
                                                    //std::cout << feas3 << "  feas6 " << endl;
                                                    if(feas6>0){
                                                        //calculate 8 step evaluation scheme
                                                        //	std::cout << "check delivery tw feasibility OK " << endl;
                                                        
                                                        gg=routad->Eightstepevaluationscheme(MatTemp, ric,  MaxRideTime,  MaxTimeRoute);
                                                        //	rout->display_route();
                                                        if(gg>0){
                                                            bestdist=distaddeffect;
                                                            bestaddpl=l;
                                                            bestadddl=g;
                                                            bestaddv=t+Solaux->numRoutes;
                                                            addinsertion=true;
                                                            
                                                        }
                                                    }
                                                    
                                                }else{
                                                    //Continue with the following insertion, This position is not considered
                                                }
                                            }
                                        }
                                        
                                    }
                                    
                                    routad->delete_delivery(g);
                                    //std:: cout << "delivery deleted" << endl;
                                }
                                
                                
                            }
                            routad->delete_pickup(l, MatTemp);
                            //std:: cout << "pickup deleted" << endl;
                            
                        }
                        
                        
                        
                        
                        
                        
                        
                    }
                    //bool addinsertion=true;
                    //if not possible, create a new one, and insert it on it.
                    if(addinsertion>0){
                        Solaux->route[bestaddv]->insert_req_pos(E[i],bestaddpl,bestadddl,ric, MatTemp);
                        Solaux->route[bestaddv]->calculate_capacity(ric,veic);
                        
                        //Solaux->route[MINI[3]]->display_route();
                        Solaux->route[bestaddv]->calculate_earliest_latest(ric, MatTemp, MaxTimeRoute);
                        Solaux->route[bestaddv]->Eightstepevaluationscheme(MatTemp,ric,MaxRideTime,MaxTimeRoute);
                        //Solaux->copyrouteinsol(rout);
                        //std::cout << "request" <<  E[i] << "inserted in vehicle additional, " << bestv<< "\n";
                        //Solaux->route[bestaddv]->display_route();
                        Solaux->updatecost(numVeicoli);
                        Solaux->updateusedvehicles(numVeicoli);
                        
                        
                    }else{
                        Solaux->addAdditionalRoute(numVeicoli, MatTemp, veic, MaxTimeRoute);
                        Solaux->route[(numVeicoli+Solaux->numAddRoutes)-1]->insert_req_pos(E[i],0,0,ric, MatTemp);
                        //std::cout << "request" <<  E[i] << "inserted in additional vehicle , " << (numVeicoli+Solaux->numAddRoutes)-1 << "\n";
                        Solaux->route[(numVeicoli+Solaux->numAddRoutes)-1]->calculate_earliest_latest(ric, MatTemp, MaxTimeRoute);
                        Solaux->route[(numVeicoli+Solaux->numAddRoutes)-1]->Eightstepevaluationscheme(MatTemp,ric,MaxRideTime,MaxTimeRoute);
                        //Solaux->copyrouteinsol(rout);
                        //std::cout << "request" <<  E[i] << "inserted in vehicle additional, " << bestv<< "\n";
                        //Solaux->route[(numVeicoli+Solaux->numAddRoutes)-1]->display_route();
                        Solaux->updatecost(numVeicoli);
                        Solaux->updateusedvehicles(numVeicoli);
                    }
                    
                }
                //INSERT IN A NEW VEHICLE
                //std:: cout << "A NEW VEHICLE HAS TO BE CREATED" << endl;
                
                //create an additional route, where the route will be ienserted.
                
                
                
            }
            
            
            
        }
        
        
        if(Solaux->numAddRoutes<=InitialSol->numAddRoutes){
            if(Solaux->cost<InitialSol->cost){
                //std::cout <<"SOLAUX" << std::endl;
                if(uu>0){
                    InitialSol->delete_routes();}else{
                    }
                InitialSol->CopySolution(Solaux);
                
            }
            
        }
        
        
    }
    
    delete routad;
    
    delete rout;
    
    delete[] E;
    //delete Solaux;
    delete Solaux;
    //std:: cout << "FIN" <<endl;
    return InitialSol;
}


Solution* InitialSolutionBraekers2(int numVeicoli, int numRichieste, RichiesteServizio* ric, Veicoli* veic, double** D, int** MatCompVei, double** MatTemp, double MinAt, int MaxTimeRoute, double MaxRideTime){
    
    
    int i, j, l , g, uu;
    
    double bestdist;
    int bestpl, bestdl, bestv;
    bool inserted;
    int bestaddpl, bestadddl,bestaddv;
    
    Solution* InitialSol=NULL;
    
    InitialSol=new Solution();
    
    InitialSol->cost=DBL_MAX;
    Route* routad=NULL;
    routad=new Route();
    
    //0
    InitialSol->numAddRoutes=INT_MAX;
    
    Solution* Solaux=NULL;
    
    Route* rout=NULL;
    rout=new Route();
    
    //Ordering vector initialized
    int* E;
    E=new int[numRichieste];
    Solaux=new Solution();
    //The initial solution is calculated 1000 times and the best solution is chosen
    for(uu=0;uu<1000; uu++){
        
        Solaux->ClearSol();
        
        
        
        
        //Initialize the solution one route per vehicle
        Solaux->add_num_empty_route(numVeicoli);
        
        Solaux->Initialize(veic, numVeicoli, MatTemp, MaxTimeRoute);
        //std::cout << "New order" << std::endl;
        //E is a vector containing the order in which the requests will be inserted
        if(uu<1){
            E=timewinminP_order(ric, numRichieste,E);
        }else{
            E=random_order( ric,  numRichieste, E);
        }
        //Print order
        //std::cout << "ordering number " << uu << endl;
        //for(i=0;i<numRichieste;i++){
        //std::cout<< E[i] << " " ;
        
        //}std::cout << "\n";
        //Each request will be inserted in all the vehicles. Cost will be kept
        for(i=0;i<numRichieste;i++){
            //	std::cout<< "request-> "<< E[i] << "\n";
            inserted=false;
            //Create a Route to insert everything on it.
            
            bool feas1, feas2, feas3,feasride, feas4, feas5, feas6;
            //double Bestdisteffect=DBL_MAX;
            bestdist=DBL_MAX;
            for(j=0;j<numVeicoli;j++){
                //std::cout<< "vehicle-> "<< j << "\n";
                // See if the vehicle is compatible with the request to be inserted.
                
                
                rout->clear_route();
                Solaux->copyvehicleinroute(rout, j);
                //
                int originallength=rout->length;
                for(l=0; l<originallength-1; l++)
                {
                    //copy current Solution in route.
                    //Solaux->copyvehicleinroute(rout, j);
                    //Insert pickup in position l
                    //rout->display_route();
                    
                    rout->insert_pickup(E[i], l, ric);
                    
                    rout->calculate_earliest_latest(ric, MatTemp, MaxTimeRoute);
                    
                    rout->calculate_capacity(ric, veic);
                    
                    feas2=rout->check_feasibility_P_cap(E[i],l, veic,ric);
                    
                    
                    
                    //	std:: cout << "inserting pickup on " << l << "\n";
                    feas1=rout->check_feasibility_P_tw(E[i],l, ric, MatTemp);
                    
                    if(feas1>0 && feas2>0){
                        //		std::cout << "pickup of request " << E[i] << "compatible in position " << l << endl;
                        //
                        
                        for (g=l; g<originallength-1; g++){
                            
                            rout->insert_delivery(E[i], g, ric);
                            
                            //rout->calculate_earliest_latest(ric, MatTemp, MaxTimeRoute);
                            
                            rout->calculate_capacity(ric, veic);
                            
                            rout->calculate_earliest_latest(ric, MatTemp, MaxTimeRoute);
                            
                            
                            feasride=rout->ridetime_check( g,  ric,  MaxRideTime);
                            
                            if(feasride>0){//std::cout <<"true" << endl;
                                double disteffect;
                                bool gg;
                                //	std:: cout << "length: " << rout->length << "l : " << l << "g : " << g << endl;
                                if(l==g){
                                    
                                    disteffect=0;
                                    disteffect=rout->effect_ondistance(E[i], l, g, D,  ric);
                                    //std:: cout << "length: " << rout->length << "l : " << l << "g : " << g << endl;
                                    //std::cout <<  disteffect << " < " << bestdist << "\n";
                                    if(disteffect<bestdist){
                                        
                                        feas3=rout->check_feasibillity_D_tw1(E[i],l,ric,MatTemp);
                                        //std::cout << feas3 << "  feas3 " << endl;
                                        
                                        if(feas3>0){
                                            //calculate 8 step evaluation scheme
                                            //	std::cout << "check delivery tw feasibility OK " << endl;
                                            
                                            gg=rout->Eightstepevaluationscheme(MatTemp, ric,  MaxRideTime,  MaxTimeRoute);
                                            //	rout->display_route();
                                            if(gg>0){
                                                bestdist=disteffect;
                                                bestpl=l;
                                                bestdl=g;
                                                bestv=j;
                                                inserted=true;
                                                
                                            }
                                        }
                                        
                                    }
                                }else{
                                    feas4=rout->check_feasibility_D_tw2( l,  g, ric, MatTemp);
                                    //	std::cout << feas3 << "  1feas3 " << endl;
                                    feas5=rout->check_feasibility_capacity(l,g,veic);
                                    //	std::cout << feas3 << "  2feas3 " << endl;
                                    //Insert E{i} in route j in position l pickup g delivery
                                    //rout.insert_req_pos(E[i], l, g, ric);
                                    if(feas4>0 && feas5>0){
                                        disteffect=rout->effect_ondistance(E[i], l, g, D,  ric);
                                        //std::cout <<  disteffect << " < " << MINI[0] << "\n";
                                        //	std:: cout << "length: " << rout->length << "l : " << l << "g : " << g << endl;
                                        //	std::cout <<  disteffect << " < " << bestdist << "\n";
                                        if(disteffect<bestdist){
                                            //calculate 8 step evaluation scheme
                                            feas6=rout->check_feasibillity_D_tw1(E[i],l,ric,MatTemp);
                                            //	std::cout << feas3 << "  feas3 " << endl;
                                            if(feas6>0){
                                                //calculate 8 step evaluation scheme
                                                //	std::cout << "check delivery tw feasibility OK " << endl;
                                                
                                                gg=rout->Eightstepevaluationscheme(MatTemp, ric,  MaxRideTime,  MaxTimeRoute);
                                                //	rout->display_route();
                                                if(gg>0){
                                                    bestdist=disteffect;
                                                    bestpl=l;
                                                    bestdl=g;
                                                    bestv=j;
                                                    inserted=true;
                                                    
                                                }
                                            }
                                            
                                        }else{
                                            //Continue with the following insertion, This position is not considered
                                        }
                                    }
                                }
                                
                            }
                            
                            rout->delete_delivery(g);
                            
                        }
                        
                        
                    }
                    rout->delete_pickup(l, MatTemp);
                    //std:: cout << "pickup deleted" << endl;
                    
                }
                
                
                
            }
            
            
            if(inserted>0){
                Solaux->route[bestv]->insert_req_pos(E[i],bestpl,bestdl,ric, MatTemp);
                Solaux->route[bestv]->calculate_capacity(ric,veic);
                
                //Solaux->route[MINI[3]]->display_route();
                Solaux->route[bestv]->calculate_earliest_latest(ric, MatTemp, MaxTimeRoute);
                Solaux->route[bestv]->Eightstepevaluationscheme(MatTemp,ric,MaxRideTime,MaxTimeRoute);
                //Solaux->copyrouteinsol(rout);
                //	std::cout << "request" << E[i] << "inserted in vehicle, " << bestv<< "\n";
                //Solaux->route[bestv]->display_route();
                Solaux->updatecost(numVeicoli);
                Solaux->updateusedvehicles(numVeicoli);
                
                
            }else{
                //std::cout<< "false" << std::endl;
                if(Solaux->numAddRoutes==0){
                    //Add a new one an insert it
                    //std::cout<< "false" << std::endl;
                    Solaux->addAdditionalRoute(numVeicoli, MatTemp, veic, MaxTimeRoute);
                    //std::cout<< numVeicoli+Solaux->numAddRoutes << std::endl;
                    Solaux->route[(numVeicoli+Solaux->numAddRoutes)-1]->insert_req_pos(E[i],0,0,ric, MatTemp);
                    //	std::cout << "request" <<  E[i] << "inserted in additional vehicle , " << (numVeicoli+Solaux->numAddRoutes)-1 << "\n";
                    //std::cout << "Add route added" << std::endl;
                    Solaux->route[(numVeicoli+Solaux->numAddRoutes)-1]->calculate_earliest_latest(ric, MatTemp, MaxTimeRoute);
                    Solaux->route[(numVeicoli+Solaux->numAddRoutes)-1]->Eightstepevaluationscheme(MatTemp,ric,MaxRideTime,MaxTimeRoute);
                    Solaux->updatecost(numVeicoli);
                    Solaux->updateusedvehicles(numVeicoli);
                    //	Solaux->route[(numVeicoli+Solaux->numAddRoutes)-1]->display_route();
                    
                }else{
                    
                    
                    double bestadddist=DBL_MAX;
                    bool addinsertion;
                    //Try to insert it in the ones that there are....
                    for(int t=0;t<Solaux->numAddRoutes; t++){
                        //try on each additional vehicle that exists
                        routad->clear_route();
                        Solaux->copyvehicleinroute(routad, t+Solaux->numRoutes);
                        //std::cout << t << std::endl;
                        addinsertion=false;
                        int originallength=routad->length;
                        for(l=0; l<originallength-1; l++)
                        {
                            
                            
                            routad->insert_pickup(E[i], l, ric);
                            routad->calculate_earliest_latest(ric, MatTemp, MaxTimeRoute);
                            
                            routad->calculate_capacity(ric, veic);
                            
                            feas2=routad->check_feasibility_P_cap(E[i],l, veic,ric);
                            
                            
                            
                            //std:: cout << "inserting pickup on " << l << "\n";
                            feas1=routad->check_feasibility_P_tw(E[i],l, ric, MatTemp);
                            
                            if(feas1>0 && feas2>0){
                                //std::cout << "pickup of request " << E[i] << "compatible in position " << l << endl;
                                //
                                
                                for (g=l; g<originallength-1; g++){
                                    
                                    routad->insert_delivery(E[i], g, ric);
                                    //		std::cout << "LENGTH " <<  rout->length << endl;
                                    //		rout->display_route();
                                    //routad->calculate_earliest_latest(ric, MatTemp, MaxTimeRoute);
                                    
                                    routad->calculate_capacity(ric, veic);
                                    //		rout->display_route();
                                    routad->calculate_earliest_latest(ric, MatTemp, MaxTimeRoute);
                                    //		std::cout << "EARLIEST and LATEST + length"  << rout->length << endl;
                                    //		rout->display_route();
                                    
                                    //checkridetime....
                                    
                                    feasride=routad->ridetime_check( g,  ric,  MaxRideTime);
                                    //	std::cout << "feasride " <<  feasride << endl;
                                    if(feasride>0){//std::cout <<"true" << endl;
                                        double distaddeffect;
                                        bool gg;
                                        //std:: cout << "length: " << rout->length << "l : " << l << "g : " << g << endl;
                                        if(l==g){
                                            
                                            distaddeffect=0;
                                            distaddeffect=routad->effect_ondistance(E[i], l, g, D,  ric);
                                            //std:: cout << "length: " << rout->length << "l : " << l << "g : " << g << endl;
                                            //std::cout <<  disteffect << " < " << bestdist << "\n";
                                            if(distaddeffect<bestadddist){
                                                
                                                feas3=routad->check_feasibillity_D_tw1(E[i],l,ric,MatTemp);
                                                //std::cout << feas3 << "  feas3 " << endl;
                                                //	routad->display_route();
                                                //std:: cout << ric[E[i]].timewinPmin << " " << ric[E[i]].timewinPmax << " " << ric[E[i]].timewinDmin << " " << ric[E[i]].timewinDmax << std::endl;
                                                if(feas3>0){
                                                    //calculate 8 step evaluation scheme
                                                    //	std::cout << "check delivery tw feasibility OK " << endl;
                                                    
                                                    gg=routad->Eightstepevaluationscheme(MatTemp, ric,  MaxRideTime,  MaxTimeRoute);
                                                    //	rout->display_route();
                                                    if(gg>0){
                                                        bestadddist=distaddeffect;
                                                        bestaddpl=l;
                                                        bestadddl=g;
                                                        bestaddv=t+Solaux->numRoutes;
                                                        addinsertion=true;
                                                        
                                                    }
                                                }
                                                
                                            }
                                        }else{
                                            feas4=routad->check_feasibility_D_tw2( l,  g, ric, MatTemp);
                                            //	std::cout << feas4 << "  1feas3 " << endl;
                                            feas5=routad->check_feasibility_capacity(l,g,veic);
                                            //	std::cout << feas5 << "  2feas3 " << endl;
                                            //Insert E{i} in route j in position l pickup g delivery
                                            //rout.insert_req_pos(E[i], l, g, ric);
                                            if(feas4>0 && feas5>0){
                                                distaddeffect=routad->effect_ondistance(E[i], l, g, D,  ric);
                                                //std::cout <<  disteffect << " < " << MINI[0] << "\n";
                                                //	std:: cout << "length: " << rout->length << "l : " << l << "g : " << g << endl;
                                                //	std::cout <<  disteffect << " < " << bestdist << "\n";
                                                if(distaddeffect<bestadddist){
                                                    //calculate 8 step evaluation scheme
                                                    feas6=routad->check_feasibillity_D_tw1(E[i],l,ric,MatTemp);
                                                    //std::cout << feas3 << "  feas6 " << endl;
                                                    if(feas6>0){
                                                        //calculate 8 step evaluation scheme
                                                        //	std::cout << "check delivery tw feasibility OK " << endl;
                                                        
                                                        gg=routad->Eightstepevaluationscheme(MatTemp, ric,  MaxRideTime,  MaxTimeRoute);
                                                        //	rout->display_route();
                                                        if(gg>0){
                                                            bestdist=distaddeffect;
                                                            bestaddpl=l;
                                                            bestadddl=g;
                                                            bestaddv=t+Solaux->numRoutes;
                                                            addinsertion=true;
                                                            
                                                        }
                                                    }
                                                    
                                                }else{
                                                    //Continue with the following insertion, This position is not considered
                                                }
                                            }
                                        }
                                        
                                    }
                                    
                                    routad->delete_delivery(g);
                                    //std:: cout << "delivery deleted" << endl;
                                }
                                
                                
                            }
                            routad->delete_pickup(l, MatTemp);
                            //std:: cout << "pickup deleted" << endl;
                            
                        }
                        
                        
                        
                        
                        
                        
                        
                    }
                    //bool addinsertion=true;
                    //if not possible, create a new one, and insert it on it.
                    if(addinsertion>0){
                        Solaux->route[bestaddv]->insert_req_pos(E[i],bestaddpl,bestadddl,ric, MatTemp);
                        Solaux->route[bestaddv]->calculate_capacity(ric,veic);
                        
                        //Solaux->route[MINI[3]]->display_route();
                        Solaux->route[bestaddv]->calculate_earliest_latest(ric, MatTemp, MaxTimeRoute);
                        Solaux->route[bestaddv]->Eightstepevaluationscheme(MatTemp,ric,MaxRideTime,MaxTimeRoute);
                        //Solaux->copyrouteinsol(rout);
                        //std::cout << "request" <<  E[i] << "inserted in vehicle additional, " << bestv<< "\n";
                        //Solaux->route[bestaddv]->display_route();
                        Solaux->updatecost(numVeicoli);
                        Solaux->updateusedvehicles(numVeicoli);
                        
                        
                    }else{
                        Solaux->addAdditionalRoute(numVeicoli, MatTemp, veic, MaxTimeRoute);
                        Solaux->route[(numVeicoli+Solaux->numAddRoutes)-1]->insert_req_pos(E[i],0,0,ric, MatTemp);
                        //std::cout << "request" <<  E[i] << "inserted in additional vehicle , " << (numVeicoli+Solaux->numAddRoutes)-1 << "\n";
                        Solaux->route[(numVeicoli+Solaux->numAddRoutes)-1]->calculate_earliest_latest(ric, MatTemp, MaxTimeRoute);
                        Solaux->route[(numVeicoli+Solaux->numAddRoutes)-1]->Eightstepevaluationscheme(MatTemp,ric,MaxRideTime,MaxTimeRoute);
                        //Solaux->copyrouteinsol(rout);
                        //std::cout << "request" <<  E[i] << "inserted in vehicle additional, " << bestv<< "\n";
                        //Solaux->route[(numVeicoli+Solaux->numAddRoutes)-1]->display_route();
                        Solaux->updatecost(numVeicoli);
                        Solaux->updateusedvehicles(numVeicoli);
                    }
                    
                }
                //INSERT IN A NEW VEHICLE
                //std:: cout << "A NEW VEHICLE HAS TO BE CREATED" << endl;
                
                //create an additional route, where the route will be ienserted.
                
                
                
            }
            
            
            
        }
        
        
        if(Solaux->numAddRoutes<=InitialSol->numAddRoutes){
            if(Solaux->cost<InitialSol->cost){
                //std::cout <<"SOLAUX" << std::endl;
                if(uu>0){
                    InitialSol->delete_routes();}else{
                    }
                InitialSol->CopySolution(Solaux);
                
            }
            
        }
        
        
    }
    
    delete routad;
    
    delete rout;
    
    delete[] E;
    //delete Solaux;
    delete Solaux;
    //std:: cout << "FIN" <<endl;
    return InitialSol;
}



Solution* InitialSolutionBraekersF(Solution* InitialSol, int numVeicoli, int numRichieste, RichiesteServizio* ric, Veicoli* veic, double** D, int** MatCompVei, double** MatTemp, double MinAt, int MaxTimeRoute, double MaxRideTime){
    
    
    int i, j, l , g, uu;
    
    double bestdist;
    int bestpl, bestdl, bestv;
    bool inserted;
    int bestaddpl, bestadddl,bestaddv;
    
    
    Solution* Solaux=NULL;
    
    //Ordering vector initialized
    int* E;
    E=new int[numRichieste];
    Solaux=new Solution();
    //The initial solution is calculated 1000 times and the best solution is chosen
    for(uu=0;uu<1000; uu++){
        
        Solaux->ClearSol();
        
        
        //Initialize the solution one route per vehicle
        Solaux->add_num_empty_route(numVeicoli);
        
        Solaux->Initialize(veic, numVeicoli, MatTemp, MaxTimeRoute);
        
        //E is a vector containing the order in which the requests will be inserted
        if(uu<1){ //First ordered according to PickupEarliest TW
            E=timewinminP_order(ric, numRichieste,E);
        }else{ //Randomly ordered
            E=random_order( ric,  numRichieste, E);
        }
     
        //Each request will be inserted in all the vehicles.
        for(i=0;i<numRichieste;i++){
           
            inserted=false;
            //Create a Route to insert everything on it.
            bool feas1, feas2, feas3,feasride, feas4, feas5, feas6;
            //double Bestdisteffect=DBL_MAX;
            bestdist=DBL_MAX;
            //Find the best vehicle for the request to be inserted
            for(j=0;j<numVeicoli;j++){
                if(MatCompVei[E[i]][j]>0){
                int originallength=Solaux->route[j]->length;
                
                for(l=0; l<originallength-1; l++)
                {
                    
                    Solaux->route[j]->insert_pickup(E[i], l, ric);
                    
                    Solaux->route[j]->calculate_earliest_latest(ric, MatTemp, MaxTimeRoute);
                    
                    Solaux->route[j]->calculate_capacity(ric, veic);
                    //Check the feasibility of the insertion of the pickup node on position l.
                   // feas1=Solaux->route[j]->check_feasibility_P(E[i],l, ric, MatTemp, veic);
                    feas2=Solaux->route[j]->check_feasibility_P_cap(E[i],l, veic,ric);
                    
                    
                    
                    //	std:: cout << "inserting pickup on " << l << "\n";
                    feas1=Solaux->route[j]->check_feasibility_P_tw(E[i],l, ric, MatTemp);
                    
                   // std::cout << "feas1&2" << feas1 << " " << feas2 << std::endl;
                    if(feas1>0 && feas2>0){
                        //If the it is feasible, insert the delivery node.
                        for (g=l; g<originallength-1; g++){
                            
                            Solaux->route[j]->insert_delivery(E[i], g, ric);
                            Solaux->route[j]->calculate_earliest_latest(ric, MatTemp, MaxTimeRoute);
                            feasride=Solaux->route[j]->ridetime_check( g,  ric,  MaxRideTime);
                           // std::cout << "ride" << feasride << std::endl;
                            if(feasride>0){//std::cout <<"true" << endl;
                                double disteffect;
                                bool gg;
                                disteffect=0;
                                feas4=Solaux->route[j]->check_feasibility_D_tw2( l,  g, ric, MatTemp);
                                Solaux->route[j]->calculate_capacity(ric, veic);
                                feas5=Solaux->route[j]->check_feasibility_capacity(l,g,veic);
                               // std::cout<< "4&5" << feas4 << " " << feas5 << std::endl;
                                if(feas4>0 && feas5>0){
                                    disteffect=Solaux->route[j]->effect_ondistance(E[i], l, g, D,  ric);
                                    std::cout<<disteffect << std::endl;
                                    if(disteffect<bestdist){
                                            //calculate 8 step evaluation scheme
                                        feas6=Solaux->route[j]->check_feasibillity_D_tw1(E[i],g,ric,MatTemp);
                                            //	std::cout << "6" << feas6 << endl;
                                        if(feas6>0){
                                                //calculate 8 step evaluation scheme
                                                //	std::cout << "check delivery tw feasibility OK " << endl;
                                                
                                            gg=Solaux->route[j]->Eightstepevaluationscheme(MatTemp, ric,  MaxRideTime,  MaxTimeRoute);
                                                //	rout->display_route();
                                           // std::cout << "gg   " << gg << std::endl;
                                            if(gg>0){
                                                bestdist=disteffect;
                                                bestpl=l;
                                                bestdl=g;
                                                bestv=j;
                                                inserted=true;
                                                    
                                            }
                                        }
                                            
                                    }else{
                                            //Continue with the following insertion, This position is not considered
                                    }
                                }
                                
                            }
                            
                            Solaux->route[j]->delete_delivery(g);
                            
                        }
                        
                        
                    }
                    Solaux->route[j]->delete_pickup(l, MatTemp);
                    //std:: cout << "pickup deleted" << endl;
                
                }
                
                }
                
            }
            
            
            if(inserted>0){
                Solaux->route[bestv]->insert_req_pos(E[i],bestpl,bestdl,ric, MatTemp);
                //Solaux->route[bestv]->calculate_capacity(ric,veic);
                //Solaux->route[bestv]->calculate_earliest_latest(ric, MatTemp, MaxTimeRoute);
                //Solaux->route[bestv]->Eightstepevaluationscheme(MatTemp,ric,MaxRideTime,MaxTimeRoute);
                //Solaux->copyrouteinsol(rout);
                	//std::cout << "request" << E[i] << "inserted in vehicle, " << bestv<< "  in pos : " << bestpl << "  " << bestdl << "\n";
                //Solaux->route[bestv]->display_route();
                Solaux->updatecost(numVeicoli);
                Solaux->updateusedvehicles(numVeicoli);
               
                
            }else{
                //std::cout<< "false" << std::endl;
                if(Solaux->numAddRoutes==0){
                    //Add a new one an insert it
                    //std::cout<< "false" << std::endl;
                    Solaux->addAdditionalRoute(numVeicoli, MatTemp, veic, MaxTimeRoute);
                    //std::cout<< numVeicoli+Solaux->numAddRoutes << std::endl;
                    Solaux->route[(numVeicoli+Solaux->numAddRoutes)-1]->insert_req_pos(E[i],0,0,ric, MatTemp);
                    	//std::cout << "request" <<  E[i] << "inserted in additional vehicle , " << (numVeicoli+Solaux->numAddRoutes)-1 << "   in pos: 0 0 \n";
                    //std::cout << "Add route added" << std::endl;
                  //  Solaux->route[(numVeicoli+Solaux->numAddRoutes)-1]->calculate_earliest_latest(ric, MatTemp, MaxTimeRoute);
                   // Solaux->route[(numVeicoli+Solaux->numAddRoutes)-1]->Eightstepevaluationscheme(MatTemp,ric,MaxRideTime,MaxTimeRoute);
                    Solaux->updatecost(numVeicoli);
                    Solaux->updateusedvehicles(numVeicoli);
                    //	Solaux->route[(numVeicoli+Solaux->numAddRoutes)-1]->display_route();
                    
                }else{
                    
                    
                    double bestadddist=DBL_MAX;
                    bool addinsertion;
                    //Try to insert it in the ones that there are....
                    for(int t=0;t<Solaux->numAddRoutes; t++){
                        //try on each additional vehicle that exists
                        //routad->clear_route();
                       // Solaux->copyvehicleinroute(routad, t+Solaux->numRoutes);
                        //std::cout << t << std::endl;
                        addinsertion=false;
                        int originallength=Solaux->route[t+Solaux->numRoutes]->length;
                        for(l=0; l<originallength-1; l++)
                        {
                            
                            
                            Solaux->route[t+Solaux->numRoutes]->insert_pickup(E[i], l, ric);
                            Solaux->route[t+Solaux->numRoutes]->calculate_earliest_latest(ric, MatTemp, MaxTimeRoute);
                            
                            Solaux->route[t+Solaux->numRoutes]->calculate_capacity(ric, veic);
                            
                            feas2=Solaux->route[t+Solaux->numRoutes]->check_feasibility_P_cap(E[i],l, veic,ric);
                            
                            
                            
                            //std:: cout << "inserting pickup on " << l << "\n";
                            feas1=Solaux->route[t+Solaux->numRoutes]->check_feasibility_P_tw(E[i],l, ric, MatTemp);
                            //feas1=Solaux->route[t+Solaux->numRoutes]->check_feasibility_P(E[i],l, ric, MatTemp, veic);
                            if(feas1>0 && feas2>0){
                                //std::cout << "pickup of request " << E[i] << "compatible in position " << l << endl;
                                //
                                
                                for (g=l; g<originallength-1; g++){
                                    
                                    Solaux->route[t+Solaux->numRoutes]->insert_delivery(E[i], g, ric);
                                    //		std::cout << "LENGTH " <<  rout->length << endl;
                                    //		rout->display_route();
                                    //routad->calculate_earliest_latest(ric, MatTemp, MaxTimeRoute);
                                    
                                    
                                    //		rout->display_route();
                                    Solaux->route[t+Solaux->numRoutes]->calculate_earliest_latest(ric, MatTemp, MaxTimeRoute);
                                    //		std::cout << "EARLIEST and LATEST + length"  << rout->length << endl;
                                    //		rout->display_route();
                                    
                                    //checkridetime....
                                    
                                    feasride=Solaux->route[t+Solaux->numRoutes]->ridetime_check( g,  ric,  MaxRideTime);
                                    //	std::cout << "feasride " <<  feasride << endl;
                                    if(feasride>0){//std::cout <<"true" << endl;
                                        double distaddeffect;
                                        bool gg;
                                        distaddeffect=0;
                                        
                                        feas4=Solaux->route[t+Solaux->numRoutes]->check_feasibility_D_tw2( l,  g, ric, MatTemp);
                                        Solaux->route[t+Solaux->numRoutes]->calculate_capacity(ric, veic);
                                        feas5=Solaux->route[t+Solaux->numRoutes]->check_feasibility_capacity(l,g,veic);
                                            //	std::cout << feas5 << "  2feas3 " << endl;
                                            //Insert E{i} in route j in position l pickup g delivery
                                            //rout.insert_req_pos(E[i], l, g, ric);
                                        if(feas4>0 && feas5>0){
                                            distaddeffect=Solaux->route[t+Solaux->numRoutes]->effect_ondistance(E[i], l, g, D,  ric);
                                                //std::cout <<  distaddeffect << "\n";
                                                //	std:: cout << "length: " << rout->length << "l : " << l << "g : " << g << endl;
                                                //	std::cout <<  disteffect << " < " << bestdist << "\n";
                                            if(distaddeffect<bestadddist){
                                                    //calculate 8 step evaluation scheme
                                                feas6=Solaux->route[t+Solaux->numRoutes]->check_feasibillity_D_tw1(E[i],g,ric,MatTemp);
                                                    //std::cout << feas3 << "  feas6 " << endl;
                                                if(feas6>0){
                                                        //calculate 8 step evaluation scheme
                                                        //	std::cout << "check delivery tw feasibility OK " << endl;
                                                        
                                                    gg=Solaux->route[t+Solaux->numRoutes]->Eightstepevaluationscheme(MatTemp, ric,  MaxRideTime,  MaxTimeRoute);
                                                        //	rout->display_route();
                                                    if(gg>0){
                                                        bestdist=distaddeffect;
                                                        bestaddpl=l;
                                                        bestadddl=g;
                                                        bestaddv=t+Solaux->numRoutes;
                                                        addinsertion=true;
                                                            
                                                    }
                                                }
                                                    
                                            }else{
                                                    //Continue with the following insertion, This position is not considered
                                            }
                                        }
                                     
                                        
                                    }
                                    
                                    Solaux->route[t+Solaux->numRoutes]->delete_delivery(g);
                                    //std:: cout << "delivery deleted" << endl;
                                }
                                
                                
                            }
                            Solaux->route[t+Solaux->numRoutes]->delete_pickup(l, MatTemp);
                            //std:: cout << "pickup deleted" << endl;
                            
                        }
                        
                        
                        
                        
                        
                        
                        
                    }
                    //bool addinsertion=true;
                    //if not possible, create a new one, and insert it on it.
                    if(addinsertion>0){
                        Solaux->route[bestaddv]->insert_req_pos(E[i],bestaddpl,bestadddl,ric, MatTemp);
                        //Solaux->route[bestaddv]->calculate_capacity(ric,veic);
                        
                        //Solaux->route[MINI[3]]->display_route();
                       // Solaux->route[bestaddv]->calculate_earliest_latest(ric, MatTemp, MaxTimeRoute);
                        //Solaux->route[bestaddv]->Eightstepevaluationscheme(MatTemp,ric,MaxRideTime,MaxTimeRoute);
                        //Solaux->copyrouteinsol(rout);
                       // std::cout << "request" <<  E[i] << "inserted in vehicle additional, " << bestaddv<< "  in pos : " << bestaddpl << "  " << bestadddl <<  "\n";
                        //Solaux->route[bestaddv]->display_route();
                        Solaux->updatecost(numVeicoli);
                        Solaux->updateusedvehicles(numVeicoli);
                        
                        
                    }else{
                        Solaux->addAdditionalRoute(numVeicoli, MatTemp, veic, MaxTimeRoute);
                        Solaux->route[(numVeicoli+Solaux->numAddRoutes)-1]->insert_req_pos(E[i],0,0,ric, MatTemp);
                        //std::cout << "request" <<  E[i] << "inserted in additional vehicle , " << (numVeicoli+Solaux->numAddRoutes)-1 << "\n";
                       // Solaux->route[(numVeicoli+Solaux->numAddRoutes)-1]->calculate_earliest_latest(ric, MatTemp, MaxTimeRoute);
                       // Solaux->route[(numVeicoli+Solaux->numAddRoutes)-1]->Eightstepevaluationscheme(MatTemp,ric,MaxRideTime,MaxTimeRoute);
                        //Solaux->copyrouteinsol(rout);
                       // std::cout << "request" <<  E[i] << "inserted in vehicle additional, " << (numVeicoli+Solaux->numAddRoutes)-1<< "   in pos: 0 0 \n";
                        //Solaux->route[(numVeicoli+Solaux->numAddRoutes)-1]->display_route();
                        Solaux->updatecost(numVeicoli);
                        Solaux->updateusedvehicles(numVeicoli);
                    }
                    
                }
                //INSERT IN A NEW VEHICLE
                //std:: cout << "A NEW VEHICLE HAS TO BE CREATED" << endl;
                
                //create an additional route, where the route will be ienserted.
                
                
                
            }
            
            
            
        }
        
        
        if(Solaux->numAddRoutes<=InitialSol->numAddRoutes){
            if(Solaux->cost<InitialSol->cost){
                //std::cout <<"SOLAUX" << std::endl;
               // if(uu>0){
                //    InitialSol->delete_routes();}else{
                  //  }
                InitialSol->CopySolution(Solaux);
              //  std::cout << "it:  "<<  uu ;
                
                //std::cout<<"   COST: " << InitialSol->cost;

                if(InitialSol->numAddRoutes>0){
                    //std::cout << "    Not feasible" <<std::endl;}else{std::cout << "    Feasible" <<std::endl;
                }
                
                for(unsigned y=0; y<InitialSol->numRoutes+InitialSol->numAddRoutes; y++){
                    InitialSol->route[y]->calculate_capacity(ric,veic);
                    InitialSol->route[y]->Eightstepevaluationscheme(MatTemp,ric,MaxRideTime,MaxTimeRoute);
                    InitialSol->route[y]->calculate_earliest_latest(ric, MatTemp, MaxTimeRoute);}
                
               // std::cout<<"COST: " << InitialSol->cost << std::endl;
                
                
            }
            
        }
        
        
    }
    
    //delete routad;
    
  //  delete rout;
    
    delete[] E;
    //delete Solaux;
    delete Solaux;
    
    
    
    //std:: cout << "FIN" <<endl;
    return InitialSol;
}

*/

SolutionVRP* InitialSolutionBraekersF2(SolutionVRP* InitialSol, int numVeicoli, int numRichieste, std::vector<RichiesteServizio*> &ric , std::vector<Veicoli*> &veic, std::vector<std::vector<double>> &D, std::vector<std::vector<int>> &MatCompVei, std::vector<std::vector<double>> &MatTemp){
    
    
    int i, j, l , g, uu;
    
    double bestdist;
    int bestpl, bestdl, bestv;
    bool inserted;
    int bestaddpl, bestadddl,bestaddv;
    
    
    SolutionVRP* Solaux=NULL;
    
    //Ordering vector initialized
    int* E;
    E=new int[numRichieste];
    Solaux=new SolutionVRP();
    int* rid;
    rid=new int[2*numRichieste+2];
    int* typ;
    typ=new int[2*numRichieste+2];
    int* loc;
    loc=new int[2*numRichieste+2];
    double* earl;
    earl=new double[2*numRichieste+2];
    //The initial solution is calculated 1000 times and the best solution is chosen
    for(uu=0;uu<1000; uu++){
        //std::cout << uu << std::endl;
        Solaux->ClearSol();
        
        //Initialize the solution one route per vehicle
        Solaux->add_num_empty_route(numVeicoli);
        
        Solaux->Initialize(veic, numVeicoli, MatTemp);
        
        //E is a vector containing the order in which the requests will be inserted
        if(uu<1){ //First ordered according to PickupEarliest TW
            E=timewinminP_order(ric, numRichieste,E);
        }else{ //Randomly ordered
            E=random_order( ric,  numRichieste, E);
        }
        
        //Each request will be inserted in all the vehicles.
        for(i=0;i<numRichieste;i++){
            
            inserted=false;
            //Create a Route to insert everything on it.
            bool feas1, feas2, feas3,feasride, feas4, feas5, feas6;
            //double Bestdisteffect=DBL_MAX;
            bestdist=DBL_MAX;
            //Find the best vehicle for the request to be inserted
            for(j=0;j<numVeicoli;j++){
                if(MatCompVei[E[i]][j]>0){
                    int originallength=Solaux->route[j]->length;
                    //std::cout<<originallength<<std::endl;
                    
                    for(l=0; l<originallength-1; l++)
                    {
                        
                        //Solaux->route[j]->insert_pickup(E[i], l, ric);
                        
                      //  Solaux->route[j]->calculate_earliest_latest(ric, MatTemp, MaxTimeRoute);
                        
                       // Solaux->route[j]->calculate_capacity(ric, veic);
                        //Check the feasibility of the insertion of the pickup node on position l.
                        // feas1=Solaux->route[j]->check_feasibility_P(E[i],l, ric, MatTemp, veic);
                        //feas2=Solaux->route[j]->check_feasibility_P_cap(E[i],l, veic,ric);
                        feas2=Solaux->route[j]->capacity_P_feasibility( l, E[i],  veic, ric);
                        
                        
                        //	std:: cout << "inserting pickup on " << l << "\n";
                        //feas1=Solaux->route[j]->check_feasibility_P_tw(E[i],l, ric, MatTemp);
                        feas1=Solaux->route[j]->tw_P_feasibility(ric, MatTemp,  l, E[i]);
                        //std::cout<<feas1 << " " << feas2 << std::endl;
                       // std::cout << "feas1&2" << feas1 << " " << feas2 << std::endl;
                        if(feas1>0 && feas2>0){
                            //If the it is feasible, insert the delivery node.
                            //for (g=l; g<originallength-1; g++){
                            g=l;
                            while(g<originallength-1){
                                double disteffect;
                                bool gg;
                                
                                if(l==g){
                                    
                                    disteffect=0;
                                    disteffect=Solaux->route[j]->effect_of_inserting_req_on_pos(E[i], l, g, D, ric);
                                    //std:: cout << "length: " << rout->length << "l : " << l << "g : " << g << endl;
                                    //std::cout <<  disteffect << " < " << bestdist << "\n";
                                    if(disteffect<bestdist){
                                        //calculate 8 step evaluation scheme
                                        //feas6=Solaux->route[j]->check_feasibillity_D_tw1(E[i],l,ric,MatTemp);
                                    
                                        
                                        
                                        rid=Solaux->route[j]->calcrid(rid, l, g, E[i]);
                                        typ=Solaux->route[j]->calctyp(typ, l, g);
                                        loc=Solaux->route[j]->calcloc(loc, E[i],ric, l, g);
                                        earl=Solaux->route[j]->calcearl(earl, E[i], ric, MatTemp, l, g, loc, rid, typ);
                                        
                                        
                                        feas6=Solaux->route[j]->check_feas_D_tw1(earl, E[i], ric, g, l, MatTemp, veic);
                                        //std::cout << "6" << feas6 << endl;
                                        if(feas6>0){
                                            //calculate 8 step evaluation scheme
                                            //	std::cout << "check delivery tw feasibility OK " << endl;
                                            
                                            //gg=Solaux->route[j]->Eightstepevaluationscheme(MatTemp, ric,  MaxRideTime,  g);
                                            gg=Solaux->route[j]->EightStepEvaluation(MatTemp, ric, l, g, E[i], veic);
                                            //	rout->display_route();
                                            //  std::cout << "gg   " << gg << std::endl;
                                            
                                            
                                            if(gg>0){
                                                bestdist=disteffect;
                                                bestpl=l;
                                                bestdl=g;
                                                bestv=j;
                                                inserted=true;
                                                
                                            }
                                        }
                                      
                                    }else{
                                        //Continue with the following insertion, This position is not considered
                                    }
                                }else{
                                
                                //std::cout << l << " " << g << std::endl;
                                    
                                rid=Solaux->route[j]->calcrid(rid, l, g, E[i]);
                                typ=Solaux->route[j]->calctyp(typ, l, g);
                                loc=Solaux->route[j]->calcloc(loc, E[i],ric, l, g);
                                earl=Solaux->route[j]->calcearl(earl, E[i], ric, MatTemp, l, g, loc, rid, typ);
                                //Solaux->route[j]->insert_delivery(E[i], g, ric);
                                //Solaux->route[j]->calculate_earliest_latest(ric, MatTemp, MaxTimeRoute);
                                //feasride=Solaux->route[j]->ridetime_check( g,  ric,  MaxRideTime);
                                feasride=Solaux->route[j]->ridetime_feas_D( g,  l, ric, E[i], MatTemp, earl);
                              //  std::cout << "ride" << feasride << std::endl;
                               // for(int yyu=0; yyu<g-l+2; yyu++){
                                //    std::cout<< rid[yyu] << "\t"<< typ[yyu] << "\t"<< loc[yyu] << "\t"<< earl[yyu] << std::endl;
                                
                               // }
                                
                                
                                
                                if(feasride>0){//std::cout <<"true" << endl;
                                    
                                    disteffect=0;
                                    //feas4=Solaux->route[j]->check_feasibility_D_tw2( l,  g, ric, MatTemp);
                                    
                                    feas4=check_feas_D_tw2(earl,  typ, rid,  l,  g, ric);
                                    
                                    //Solaux->route[j]->calculate_capacity(ric, veic);
                                    //feas5=Solaux->route[j]->check_feasibility_capacity(l,g,veic);
                                    
                                    feas5=Solaux->route[j]->check_cap_from(l, g, veic, E[i], ric);
                                    
                                    //std::cout<< "4&5" << feas4 << " " << feas5 << std::endl;
                                    
                                    if(feas4>0 && feas5>0){
                                        //disteffect=Solaux->route[j]->effect_ondistance(E[i], l, g, D,  ric);
                                        disteffect=Solaux->route[j]->effect_of_inserting_req_on_pos(E[i], l, g, D, ric);
                                        
                                       // std::cout << disteffect << std::endl;
                                        
                                        if(disteffect<bestdist){
                                            //calculate 8 step evaluation scheme
                                            //feas6=Solaux->route[j]->check_feasibillity_D_tw1(E[i],l,ric,MatTemp);
                                            feas6=Solaux->route[j]->check_feas_D_tw1(earl, E[i], ric, g, l, MatTemp, veic);
                                            	//std::cout << "6" << feas6 << endl;
                                            if(feas6>0){
                                                //calculate 8 step evaluation scheme
                                                //	std::cout << "check delivery tw feasibility OK " << endl;
                                                
                                                //gg=Solaux->route[j]->Eightstepevaluationscheme(MatTemp, ric,  MaxRideTime,  g);
                                                gg=Solaux->route[j]->EightStepEvaluation(MatTemp, ric, l, g, E[i], veic);
                                                //	rout->display_route();
                                              //  std::cout << "gg   " << gg << std::endl;
                                                
                                                
                                                if(gg>0){
                                                    bestdist=disteffect;
                                                    bestpl=l;
                                                    bestdl=g;
                                                    bestv=j;
                                                    inserted=true;
                                                    
                                                }
                                            }
                                            
                                        }else{
                                            //Continue with the following insertion, This position is not considered
                                        }
                                    }else{
                                        
                                        g=originallength-1;
                                    }
                                    
                                }else{
                                    g=originallength-1;
                                    
                                }
                                
                               // Solaux->route[j]->delete_delivery(g);
                                
    
                            }
                                g+=1;

                            }//end while
                            
                            
                        }
                       // Solaux->route[j]->delete_pickup(l, MatTemp);
                        //std:: cout << "pickup deleted" << endl;
                        
                    }
                    
                }
                
            }
            
            
            if(inserted>0){
                Solaux->route[bestv]->insert_req_pos(E[i],bestpl,bestdl,ric, MatTemp);
                Solaux->route[bestv]->calculate_capacity(ric,veic);
                Solaux->route[bestv]->calculate_earliest_latest(ric, MatTemp, veic);
                Solaux->route[bestv]->Eightstepevaluationscheme(MatTemp,ric,veic);
                //Solaux->copyrouteinsol(rout);
                	//std::cout << "request" << E[i] << "inserted in vehicle, " << bestv<< "  in pos : " << bestpl << "  " << bestdl << "\n";
                //Solaux->route[bestv]->display_route();
               // std::cout<< "VEH:  " << bestv << " " << bestpl << " " << bestdl << " " << std::endl;
                
                Solaux->updatecost(numVeicoli);
                Solaux->updateusedvehicles(numVeicoli);
                
                
            }else{
               // std::cout<< "false" << std::endl;
                if(Solaux->numAddRoutes==0){
                    //Add a new one an insert it
                    //std::cout<< "false" << std::endl;
                    Solaux->addAdditionalRoute(numVeicoli, MatTemp, veic);
                    //std::cout<< numVeicoli+Solaux->numAddRoutes << std::endl;
                    Solaux->route[(numVeicoli+Solaux->numAddRoutes)-1]->insert_req_pos(E[i],0,0,ric, MatTemp);
                    	//std::cout << "request" <<  E[i] << "inserted in additional vehicle , " << (numVeicoli+Solaux->numAddRoutes)-1 << "  in pos : 0 0 \n";
                    //std::cout << "Add route added" << std::endl;
                      Solaux->route[(numVeicoli+Solaux->numAddRoutes)-1]->calculate_earliest_latest(ric, MatTemp, veic);
                     Solaux->route[(numVeicoli+Solaux->numAddRoutes)-1]->Eightstepevaluationscheme(MatTemp,ric,veic);
                    Solaux->updatecost(numVeicoli);
                    Solaux->updateusedvehicles(numVeicoli);
                    //	Solaux->route[(numVeicoli+Solaux->numAddRoutes)-1]->display_route();
                    
                }else{
                    
                    
                    double bestadddist=DBL_MAX;
                    bool addinsertion;
                    //Try to insert it in the ones that there are....
                    for(int t=0;t<Solaux->numAddRoutes; t++){
                        //try on each additional vehicle that exists
                        //routad->clear_route();
                        // Solaux->copyvehicleinroute(routad, t+Solaux->numRoutes);
                        //std::cout << t << std::endl;
                        addinsertion=false;
                        int originallength=Solaux->route[t+Solaux->numRoutes]->length;
                        for(l=0; l<originallength-1; l++)
                        {
                            feas2=Solaux->route[t+Solaux->numRoutes]->capacity_P_feasibility( l, E[i],  veic, ric);
                            
                            
                            //	std:: cout << "inserting pickup on " << l << "\n";
                            //feas1=Solaux->route[j]->check_feasibility_P_tw(E[i],l, ric, MatTemp);
                            feas1=Solaux->route[t+Solaux->numRoutes]->tw_P_feasibility(ric, MatTemp,  l, E[i]);
                            
                          
                            if(feas1>0 && feas2>0){
                                //std::cout << "pickup of request " << E[i] << "compatible in position " << l << endl;
                                //
                                g=l;
                               // for (g=l; g<originallength-1; g++){
                                    while(g<originallength-1){
                                    double distaddeffect;
                                    bool gg;
                                    if(l==g){
                                        distaddeffect=0;
                                        distaddeffect=Solaux->route[t+Solaux->numRoutes]->effect_of_inserting_req_on_pos(E[i], l, g, D, ric);
                                        //std:: cout << "length: " << rout->length << "l : " << l << "g : " << g << endl;
                                        //std::cout <<  disteffect << " < " << bestdist << "\n";
                                        if(distaddeffect<bestadddist){
                                            //calculate 8 step evaluation scheme
                                            //feas6=Solaux->route[j]->check_feasibillity_D_tw1(E[i],l,ric,MatTemp);

                                            
                                            
                                            rid=Solaux->route[t+Solaux->numRoutes]->calcrid(rid, l, g, E[i]);
                                            typ=Solaux->route[t+Solaux->numRoutes]->calctyp(typ, l, g);
                                            loc=Solaux->route[t+Solaux->numRoutes]->calcloc(loc, E[i],ric, l, g);
                                            earl=Solaux->route[t+Solaux->numRoutes]->calcearl(earl, E[i], ric, MatTemp, l, g, loc, rid, typ);
                                            
                                            
                                            feas6=Solaux->route[t+Solaux->numRoutes]->check_feas_D_tw1(earl, E[i], ric, g, l, MatTemp,veic);
                                            //std::cout << "6" << feas6 << endl;
                                            if(feas6>0){
                                                //calculate 8 step evaluation scheme
                                                //	std::cout << "check delivery tw feasibility OK " << endl;
                                                
                                                //gg=Solaux->route[j]->Eightstepevaluationscheme(MatTemp, ric,  MaxRideTime,  g);
                                                gg=Solaux->route[t+Solaux->numRoutes]->EightStepEvaluation(MatTemp, ric, l, g, E[i],veic);
                                                //	rout->display_route();
                                                //  std::cout << "gg   " << gg << std::endl;
                                                
                                                
                                                if(gg>0){
                                                    bestdist=distaddeffect;
                                                    bestaddpl=l;
                                                    bestadddl=g;
                                                    bestaddv=t+Solaux->numRoutes;
                                                    addinsertion=true;
                                                    
                                                }
                                            }
                                        
     
                                        }
                                    
                                    }else{
                                
                                    rid=Solaux->route[t+Solaux->numRoutes]->calcrid(rid, l, g, E[i]);
            
                                    
                                    typ=Solaux->route[t+Solaux->numRoutes]->calctyp(typ, l, g);
                                    loc=Solaux->route[t+Solaux->numRoutes]->calcloc(loc, E[i],ric, l, g);
                                    earl=Solaux->route[t+Solaux->numRoutes]->calcearl(earl, E[i], ric, MatTemp, l, g, loc, rid, typ);
                                    //Solaux->route[j]->insert_delivery(E[i], g, ric);
                                    //Solaux->route[j]->calculate_earliest_latest(ric, MatTemp, MaxTimeRoute);
                                    //feasride=Solaux->route[j]->ridetime_check( g,  ric,  MaxRideTime);
                                    feasride=Solaux->route[t+Solaux->numRoutes]->ridetime_feas_D( g,  l, ric, E[i], MatTemp, earl);

                                    
                                    
                                    if(feasride>0){//std::cout <<"true" << endl;
                                        
                                        distaddeffect=0;
                                        feas4=check_feas_D_tw2(earl,  typ, rid,  l,  g, ric);
                                        
                                        //Solaux->route[j]->calculate_capacity(ric, veic);
                                        //feas5=Solaux->route[j]->check_feasibility_capacity(l,g,veic);
                                        
                                        feas5=Solaux->route[t+Solaux->numRoutes]->check_cap_from(l, g, veic, E[i], ric);
                                        if(feas4>0 && feas5>0){
                                           // distaddeffect=Solaux->route[t+Solaux->numRoutes]->effect_ondistance(E[i], l, g, D,  ric);
                                            distaddeffect=Solaux->route[t+Solaux->numRoutes]->effect_of_inserting_req_on_pos(E[i], l, g, D, ric);
                                            //std::cout <<  disteffect << " < " << MINI[0] << "\n";
                                            //	std:: cout << "length: " << rout->length << "l : " << l << "g : " << g << endl;
                                            	//std::cout <<  distaddeffect  << " < " << bestdist << "\n";
                                            if(distaddeffect<bestadddist){
                                                //calculate 8 step evaluation scheme
                                                //feas6=Solaux->route[t+Solaux->numRoutes]->check_feasibillity_D_tw1(E[i],l,ric,MatTemp);
                                                feas6=Solaux->route[t+Solaux->numRoutes]->check_feas_D_tw1(earl, E[i], ric, g, l, MatTemp, veic);
                                                //std::cout << feas3 << "  feas6 " << endl;
                                                if(feas6>0){
                                                    //calculate 8 step evaluation scheme
                                                    //	std::cout << "check delivery tw feasibility OK " << endl;
                                                    gg=Solaux->route[t+Solaux->numRoutes]->EightStepEvaluation(MatTemp, ric, l, g, E[i], veic);
                                                    //gg=Solaux->route[t+Solaux->numRoutes]->Eightstepevaluationscheme(MatTemp, ric,  MaxRideTime,  MaxTimeRoute);
                                                    //	rout->display_route();
                                                    if(gg>0){
                                                        bestdist=distaddeffect;
                                                        bestaddpl=l;
                                                        bestadddl=g;
                                                        bestaddv=t+Solaux->numRoutes;
                                                        addinsertion=true;
                                                        
                                                    }
                                                }
                                                
                                            }else{
                                                //Continue with the following insertion, This position is not considered
                                            }
                                        }else{
                                        
                                            g=originallength-1;
                                        }
    
                                    }else{
                                        g=originallength-1;
                                    
                                    }
   
                                }
                                        
                                        g+=1;
                                }// end of the while
                                
                            }
                           // Solaux->route[t+Solaux->numRoutes]->delete_pickup(l, MatTemp);
                            //std:: cout << "pickup deleted" << endl;
                            
                        }
                        
                        
                        
                        
                        
                        
                        
                    }
                    //bool addinsertion=true;
                    //if not possible, create a new one, and insert it on it.
                    if(addinsertion>0){
                        Solaux->route[bestaddv]->insert_req_pos(E[i],bestaddpl,bestadddl,ric, MatTemp);
                        Solaux->route[bestaddv]->calculate_capacity(ric,veic);
                        
                        //Solaux->route[MINI[3]]->display_route();
                         Solaux->route[bestaddv]->calculate_earliest_latest(ric, MatTemp, veic);
                        Solaux->route[bestaddv]->Eightstepevaluationscheme(MatTemp,ric, veic);
                        //Solaux->copyrouteinsol(rout);
                      //  std::cout << "request" <<  E[i] << "inserted in vehicle additional, " << bestaddv<< "  in pos : " << bestaddpl << "  " << bestadddl <<  "\n";
                        //Solaux->route[bestaddv]->display_route();
                        Solaux->updatecost(numVeicoli);
                        Solaux->updateusedvehicles(numVeicoli);
                        
                        
                    }else{
                        Solaux->addAdditionalRoute(numVeicoli, MatTemp, veic);
                        Solaux->route[(numVeicoli+Solaux->numAddRoutes)-1]->insert_req_pos(E[i],0,0,ric, MatTemp);
                        //std::cout << "request" <<  E[i] << "inserted in additional vehicle , " << (numVeicoli+Solaux->numAddRoutes)-1 << "\n";
                         Solaux->route[(numVeicoli+Solaux->numAddRoutes)-1]->calculate_earliest_latest(ric, MatTemp, veic);
                         Solaux->route[(numVeicoli+Solaux->numAddRoutes)-1]->Eightstepevaluationscheme(MatTemp,ric,veic);
                        //Solaux->copyrouteinsol(rout);
                        //std::cout << "request" <<  E[i] << "inserted in vehicle additional, " << (numVeicoli+Solaux->numAddRoutes)-1<< "   in pos: 0 0 \n";
                        //Solaux->route[(numVeicoli+Solaux->numAddRoutes)-1]->display_route();
                        Solaux->updatecost(numVeicoli);
                        Solaux->updateusedvehicles(numVeicoli);
                    }
                    
                }
                //INSERT IN A NEW VEHICLE
                //std:: cout << "A NEW VEHICLE HAS TO BE CREATED" << endl;
                
                //create an additional route, where the route will be ienserted.
                
                
                
            }
            
            
            
        }
        
        
        if(Solaux->numAddRoutes<=InitialSol->numAddRoutes){
            if(*Solaux<*InitialSol){
                
                
                //std::cout <<"SOLAUX" << std::endl;
                // if(uu>0){
                //    InitialSol->delete_routes();}else{
                //  }
                InitialSol->CopySolution(Solaux);
               // std::cout << "it:  "<<  uu ;
                
                //std::cout<<"   COST: " << InitialSol->cost;
                
                if(InitialSol->numAddRoutes>0){
                   // std::cout << "    Not feasible" <<std::endl;}else{std::cout << "    Feasible" <<std::endl;
                }
                //InitialSol->DisplaySolution();
                
               // for(unsigned y=0; y<InitialSol->numRoutes+InitialSol->numAddRoutes; y++){
                 //   InitialSol->route[y]->calculate_capacity(ric,veic);
                  //  InitialSol->route[y]->Eightstepevaluationscheme(MatTemp,ric,veic);
                  //  InitialSol->route[y]->calculate_earliest_latest(ric, MatTemp, veic);}
            }
            
        }
        
        
    }
    
    //delete routad;
    
    //  delete rout;
    delete[] rid;
    delete[] typ;
    delete[] loc;
    delete[] earl;
    delete[] E;
    //delete Solaux;
    delete Solaux;
    
    std::cout <<"INIT: " << InitialSol->getSolutionValue() << std::endl;
    //std::cout << InitialSol->numAddRoutes << std::endl;

    return InitialSol;
}


emili::Solution* IS::generateSolution()
{
    SolutionVRP* sol = new SolutionVRP();
    return InitialSolutionBraekersF2(sol,inst.numVeicoli0,inst.numRichieste0,inst.rc,inst.vec,inst.Dist,inst.MCV,inst.T);
}



emili::Solution* IS::generateEmptySolution()
{
    return new SolutionVRP();
}
