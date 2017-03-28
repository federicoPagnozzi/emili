//
//  Neighborhood.cpp
//  OrCode
//
//  Created by garazi on 25/01/2017.
//  Copyright © 2017 Garazi. All rights reserved.
//

#include "Neighborhood.hpp"
#include<iostream>
#include<fstream>
#include<string>
#include<cmath>
#include <cstring>
#include <ctime>
#include <float.h>
#include "SolutionVRP.hpp"
#include "matrix.hpp"
#include "Route.hpp"
/*
Solution* Relocate_Neighborhood2(Solution* Sol, int numVeicoli, int** MatCompVei, RichiesteServizio* ric, double** MatTemp, int MaxTimeRoute, Veicoli* veic, int MaxRideTime, double** D){
    
    int vei;
    int i, j, l, g;
    
    do{
        vei=rand()%(numVeicoli+Sol->numAddRoutes);
        //take out from this route all the requests and insert them
    }while(Sol->route[vei]->numRicRoute==0);
    //std::cout<< vei << std:: endl;
    
    //take out all the requests from the route
    int* E;
    
    int numofreq=0;
    numofreq=Sol->route[vei]->numRicRoute;
    E=new int[numofreq];
    
    E=Sol->route[vei]->count_request(E);
    
    //std::cout<< numofreq << std:: endl;
    //for(i=0; i<numofreq; i++){
    //std::cout<< E[i] << std:: endl;
    //}
    
    Sol->delete_route(vei,numVeicoli);
    //Sol->DisplaySolution();
    
    int bestpl, bestdl, bestv;
    int bestaddpl, bestadddl, bestaddv;
    bool inserted,changed;
    changed=true;
    double bestdist;
    
    
    for(i=0;i<numofreq;i++){
        
        //std::cout<< "request-> "<< E[i] << "\n";
        inserted=false;
        //Create a Route to insert everything on it.
        
        bool feas1, feas2, feas3, compatibility,feasride, feas4, feas5, feas6;
        //double Bestdisteffect=DBL_MAX;
        bestdist=DBL_MAX;
        for(j=0;j<numVeicoli;j++){
            if(j!=vei){
                //std::cout<< "vehicle-> "<< j << "\n";
                // See if the vehicle is compatible with the request to be inserted.
                if(MatCompVei[E[i]][j]==1){
                    
                    compatibility=true;
                    //	std::cout << "request " << E[i] << "compatible in vehicle " << j << std::endl;
                    
                }else{
                    //	std::cout << "request " << E[i] << "NO compatible in vehicle " << j << std::endl;
                    compatibility=false;
                    
                }
                if(compatibility==true){
                    Route* rout=NULL;
                    rout=new Route();
                    rout=Sol->copyvehicleinroute(rout, j);
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
                        
                        if(feas1==true && feas2==true){
                            //		std::cout << "pickup of request " << E[i] << "compatible in position " << l << endl;
                            //
                            
                            for (g=l; g<originallength-1; g++){
                                
                                rout->insert_delivery(E[i], g, ric);
                                
                                rout->calculate_earliest_latest(ric, MatTemp, MaxTimeRoute);
                                
                                rout->calculate_capacity(ric, veic);
                                
                                rout->calculate_earliest_latest(ric, MatTemp, MaxTimeRoute);
                                
                                
                                feasride=rout->ridetime_check( g,  ric,  MaxRideTime);
                                
                                if(feasride==true){//std::cout <<"true" << endl;
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
                                            
                                            if(feas3==true){
                                                //calculate 8 step evaluation scheme
                                                //	std::cout << "check delivery tw feasibility OK " << endl;
                                                
                                                gg=rout->Eightstepevaluationscheme(MatTemp, ric,  MaxRideTime,  MaxTimeRoute);
                                                //	rout->display_route();
                                                if(gg==true){
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
                                        if(feas4==true && feas5==true){
                                            disteffect=rout->effect_ondistance(E[i], l, g, D,  ric);
                                            //std::cout <<  disteffect << " < " << MINI[0] << "\n";
                                            //	std:: cout << "length: " << rout->length << "l : " << l << "g : " << g << endl;
                                            //	std::cout <<  disteffect << " < " << bestdist << "\n";
                                            if(disteffect<bestdist){
                                                //calculate 8 step evaluation scheme
                                                feas6=rout->check_feasibillity_D_tw1(E[i],l,ric,MatTemp);
                                                //	std::cout << feas3 << "  feas3 " << endl;
                                                if(feas6==true){
                                                    //calculate 8 step evaluation scheme
                                                    //	std::cout << "check delivery tw feasibility OK " << endl;
                                                    
                                                    gg=rout->Eightstepevaluationscheme(MatTemp, ric,  MaxRideTime,  MaxTimeRoute);
                                                    //	rout->display_route();
                                                    if(gg==true){
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
                    
                    delete rout;
                    
                }
                
            }
        }
        if(inserted==true){
            //	std::cout<< "true" << std::endl;
            
            Sol->route[bestv]->insert_req_pos(E[i],bestpl,bestdl,ric, MatTemp);
            Sol->route[bestv]->calculate_capacity(ric,veic);
            
            
            Sol->route[bestv]->calculate_earliest_latest(ric, MatTemp, MaxTimeRoute);
            Sol->route[bestv]->Eightstepevaluationscheme(MatTemp,ric,MaxRideTime,MaxTimeRoute);
            //Solaux->copyrouteinsol(rout);
            //	std::cout << "request" << E[i] << "inserted in vehicle, " << bestv<< "\n";
            //Solaux->route[bestv]->display_route();
            Sol->updatecost(numVeicoli);
            Sol->updateusedvehicles(numVeicoli);
            //Sol->route[bestv]->display_route();
            //std::cin >> j;
        }else{
            //	std::cout<< "false" << std::endl;
            
            if(Sol->numAddRoutes==0){
                //Add a new one an insert it
                //std::cout<< "false" << std::endl;
                Sol->addAdditionalRoute(numVeicoli, MatTemp, veic, MaxTimeRoute);
                //	std::cout<< numVeicoli+Sol->numAddRoutes << std::endl;
                
                //	Sol->route[(numVeicoli+Sol->numAddRoutes)-1]->display_route();
                //	std::cout << "req" << E[i] << std::endl;
                //	std::cout << Sol->route[(numVeicoli+Sol->numAddRoutes)-1]->length << std::endl;
                Sol->route[(numVeicoli+Sol->numAddRoutes)-1]->insert_req_pos(E[i],0,0,ric, MatTemp);
                //Sol->route[(numVeicoli+Sol->numAddRoutes)-1]->display_route();
                //	std::cout << Sol->route[(numVeicoli+Sol->numAddRoutes)-1]->length << std::endl;
                //std::cin >> j;
                //	std::cout << "request" <<  E[i] << "inserted in additional vehicle , " << (numVeicoli+Solaux->numAddRoutes)-1 << "\n";
                //std::cout << "Add route added" << std::endl;
                Sol->route[(numVeicoli+Sol->numAddRoutes)-1]->calculate_earliest_latest(ric, MatTemp, MaxTimeRoute);
                Sol->route[(numVeicoli+Sol->numAddRoutes)-1]->Eightstepevaluationscheme(MatTemp,ric,MaxRideTime,MaxTimeRoute);
                Sol->updatecost(numVeicoli);
                Sol->updateusedvehicles(numVeicoli);
                //Sol->route[(numVeicoli+Sol->numAddRoutes)-1]->display_route();
                
            }else{
                
                
                double bestadddist=DBL_MAX;
                bool addinsertion;
                //Try to insert it in the ones that there are....
                for(int t=0;t<Sol->numAddRoutes; t++){
                    //try on each additional vehicle that exists
                    Route* routad=NULL;
                    routad=new Route();
                    routad=Sol->copyvehicleinroute(routad, t+Sol->numRoutes);
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
                        
                        if(feas1==true && feas2==true){
                            //std::cout << "pickup of request " << E[i] << "compatible in position " << l << endl;
                            //
                            
                            for (g=l; g<originallength-1; g++){
                                
                                routad->insert_delivery(E[i], g, ric);
                                //		std::cout << "LENGTH " <<  rout->length << endl;
                                //		rout->display_route();
                                routad->calculate_earliest_latest(ric, MatTemp, MaxTimeRoute);
                                
                                routad->calculate_capacity(ric, veic);
                                //		rout->display_route();
                                routad->calculate_earliest_latest(ric, MatTemp, MaxTimeRoute);
                                //		std::cout << "EARLIEST and LATEST + length"  << rout->length << endl;
                                //		rout->display_route();
                                
                                //checkridetime....
                                
                                feasride=routad->ridetime_check( g,  ric,  MaxRideTime);
                                //	std::cout << "feasride " <<  feasride << endl;
                                if(feasride==true){//std::cout <<"true" << endl;
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
                                            if(feas3==true){
                                                //calculate 8 step evaluation scheme
                                                //	std::cout << "check delivery tw feasibility OK " << endl;
                                                
                                                gg=routad->Eightstepevaluationscheme(MatTemp, ric,  MaxRideTime,  MaxTimeRoute);
                                                //	rout->display_route();
                                                if(gg==true){
                                                    bestadddist=distaddeffect;
                                                    bestaddpl=l;
                                                    bestadddl=g;
                                                    bestaddv=t+Sol->numRoutes;
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
                                        if(feas4==true && feas5==true){
                                            distaddeffect=routad->effect_ondistance(E[i], l, g, D,  ric);
                                            //std::cout <<  disteffect << " < " << MINI[0] << "\n";
                                            //	std:: cout << "length: " << rout->length << "l : " << l << "g : " << g << endl;
                                            //	std::cout <<  disteffect << " < " << bestdist << "\n";
                                            if(distaddeffect<bestadddist){
                                                //calculate 8 step evaluation scheme
                                                feas6=routad->check_feasibillity_D_tw1(E[i],l,ric,MatTemp);
                                                //std::cout << feas3 << "  feas6 " << endl;
                                                if(feas6==true){
                                                    //calculate 8 step evaluation scheme
                                                    //	std::cout << "check delivery tw feasibility OK " << endl;
                                                    
                                                    gg=routad->Eightstepevaluationscheme(MatTemp, ric,  MaxRideTime,  MaxTimeRoute);
                                                    //	rout->display_route();
                                                    if(gg==true){
                                                        bestdist=distaddeffect;
                                                        bestaddpl=l;
                                                        bestadddl=g;
                                                        bestaddv=t+Sol->numRoutes;
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
                    
                    delete routad;
                    
                    
                    
                    
                    
                }
                //bool addinsertion=true;
                //if not possible, create a new one, and insert it on it.
                if(addinsertion==true){
                    Sol->route[bestaddv]->insert_req_pos(E[i],bestaddpl,bestadddl,ric, MatTemp);
                    Sol->route[bestaddv]->calculate_capacity(ric,veic);
                    
                    //Solaux->route[MINI[3]]->display_route();
                    Sol->route[bestaddv]->calculate_earliest_latest(ric, MatTemp, MaxTimeRoute);
                    Sol->route[bestaddv]->Eightstepevaluationscheme(MatTemp,ric,MaxRideTime,MaxTimeRoute);
                    //Solaux->copyrouteinsol(rout);
                    //std::cout << "request" <<  E[i] << "inserted in vehicle additional, " << bestv<< "\n";
                    //Solaux->route[bestaddv]->display_route();
                    Sol->updatecost(numVeicoli);
                    Sol->updateusedvehicles(numVeicoli);
                    
                    
                }else{
                    Sol->addAdditionalRoute(numVeicoli, MatTemp, veic, MaxTimeRoute);
                    Sol->route[(numVeicoli+Sol->numAddRoutes)-1]->insert_req_pos(E[i],0,0,ric, MatTemp);
                    //std::cout << "request" <<  E[i] << "inserted in additional vehicle , " << (numVeicoli+Solaux->numAddRoutes)-1 << "\n";
                    Sol->route[(numVeicoli+Sol->numAddRoutes)-1]->calculate_earliest_latest(ric, MatTemp, MaxTimeRoute);
                    Sol->route[(numVeicoli+Sol->numAddRoutes)-1]->Eightstepevaluationscheme(MatTemp,ric,MaxRideTime,MaxTimeRoute);
                    //Solaux->copyrouteinsol(rout);
                    //std::cout << "request" <<  E[i] << "inserted in vehicle additional, " << bestv<< "\n";
                    //Solaux->route[(numVeicoli+Solaux->numAddRoutes)-1]->display_route();
                    Sol->updatecost(numVeicoli);
                    Sol->updateusedvehicles(numVeicoli);
                }
                
            }
            //INSERT IN A NEW VEHICLE
            //std:: cout << "A NEW VEHICLE HAS TO BE CREATED" << endl;
            
            //create an additional route, where the route will be ienserted.
            
            
            
        }
        
        
        
    }
    
    
    
    
    //
    delete[] E;
    return Sol;
}
*/

SolutionVRP* Exchange_F(SolutionVRP* Sol, int numVeicoli, std::vector<RichiesteServizio*> & ric, std::vector<std::vector<double>> & MatTemp,std::vector<Veicoli*> &veic, int numRichieste, std::vector<std::vector<int>> &MatCompVei){
    
    int vei;
    //randomly chose a route...
    do{
        vei=rand()%(numVeicoli+Sol->numAddRoutes);
        //take out from this route all the requests and insert them
    }while(Sol->route[vei]->numRicRoute==0);
    
    int* rid;
    rid=new int[2*numRichieste+2];
    int* typ;
    typ=new int[2*numRichieste+2];
    int* loc;
    loc=new int[2*numRichieste+2];
    double* earl;
    earl=new double[2*numRichieste+2];
    
    
    int* E;
    
    int numofreq=0;
    numofreq=Sol->route[vei]->numRicRoute;
    E=new int[numofreq];
    int* F=new int[numRichieste];
    E=Sol->route[vei]->count_request(E);
    bool feas2, feas1, feas6;
    //First Swap
    
    for(int i=0;i<numofreq;i++){
        int p,d, p1, d1;
        p1=Sol->route[vei]->find_pickup(E[i]);
        d1=Sol->route[vei]->find_delivery(E[i]);
        
        for(int j=0; j<numVeicoli; j++){
            if(j!=vei){
                
                if(MatCompVei[E[i]][j]>0){
                    int numofreq2=0;
                    numofreq2=Sol->route[j]->numRicRoute;
                
                    F=Sol->route[j]->count_request(F);
                
                    for(int k=0; k<numofreq2; k++){
                    
                    //first thing, check feasibility of inserting request F[k] in positions p1 and d1, if it is feasible. //replace it in the end...

                        if(MatCompVei[F[k]][vei]>0){
                        
                    //if it is feasible, check where you can insert in the best position.
                            feas2=Sol->route[vei]->capacity_feasibility2(p1, F[k], veic, ric);
                        
                            feas1=Sol->route[vei]->tw_P_feasibility_2(ric, MatTemp,  p1, F[k]);
                        
                            if(feas1>0 && feas2>0){
                                
                                if(p1==d1+1){
                                    
                                    //feas6=Sol->route[vei]->checkfeasD_tw1();
                                    
                                    
                                    
                                
                                }else{
                                
                                
                                
                                
                                
                                }
                                
                                
                                
                                
                                
                        
                            }
                        
                        
                    
                    
                        }
                    
                    }
                
                
                
                
                }
            }
            
        }
        
    }
    
    //
    delete[] F;
    delete[] E;

    
    
    
    
    return Sol;
}





SolutionVRP* two_opt(SolutionVRP* Sol, int numVeicoli, std::vector<std::vector<double>> & D, std::vector<Veicoli*> &veic, std::vector<RichiesteServizio*> &ric, int numRichieste){
    int i, j, k, l;
    bool found=false;
    int v1, v2, pos1, pos2;
    double disteffect;
    double d;
    bool f1, f2, f3, f4;
    int numveh=Sol->numAddRoutes+Sol->numRoutes;
    disteffect=DBL_MAX;

    
    for(i=0;i<numveh;i++){
        for(j=1;j<Sol->route[i]->length-2;j++){
            if(Sol->route[i]->cap1[j]+Sol->route[i]->cap2[j]+Sol->route[i]->cap3[j]+Sol->route[i]->cap4[j]<1){
                //std::cout << "after position " << j << std::endl;
                for(k=0;k<numveh;k++){
                    //disteffect=Sol->route[i]->totaldist+Sol->route[k]->totaldist;
                    if(k!=i){
                        for(l=1;l<Sol->route[k]->length-2;l++){
                            if(Sol->route[k]->cap1[l]+Sol->route[k]->cap2[l]+Sol->route[k]->cap3[l]+Sol->route[k]->cap4[l]<1){
                                
                                //no esta bien->
                                d=Sol->calculate_effect_2(i,j,k,l,D);
                                
                                if(d<disteffect){
                                    f1=Sol->cap_feas(i,j,k,l, veic);
                                    //std::cout << f1 << std::endl;
                                    if(f1==true){
                                        f2=Sol->tw_feas(i,j,k,l,ric, D);
                                        if(f2==true){
                 
                                            f3=Sol->Eightstepevaluationscheme_newr(D, ric,  i, j, k, l,  numRichieste, veic);
                                            f4=Sol->Eightstepevaluationscheme_newr(D, ric, k,l,i, j, numRichieste, veic);
                                            if(f3==true && f4==true){
                                                disteffect=d;
                                                found=true;
                                                v1=i;
                                                v2=k;
                                                pos1=j;
                                                pos2=l;
                                            }
 
                                        }
                                    }
                                }
                            }
                        }
                    }
                    
                }
               
            }
        }
        
    }
    if(found==true){
        Sol->changetworoutes(v1,pos1, v2, pos2, ric, D, veic);
        
        
    }
    //if inserted is OK then redo all the routes and try
    
    return Sol;
}


SolutionVRP* r_4_opt(SolutionVRP* Sol, int numVeicoli, std::vector<std::vector<double>> & D, std::vector<RichiesteServizio*> &ric, std::vector<Veicoli*> &veic){
    int vei;
    
    double ordist, newdist;
    
    do{
            vei=rand()%(numVeicoli+Sol->numAddRoutes);
            //take out from this route all the requests and insert them
        }while(Sol->route[vei]->numRicRoute==0);

       // Route* Rout=NULL;
        //Rout= new Route();
        int best_i, best_j, best_z, best_u;
        bool changed=false;
        ordist=Sol->route[vei]->totaldist;

        for(int u=1; u<Sol->route[vei]->length-3; u++){

           // Rout=Sol->copyvehicleinroute(Rout, vei);
            //u=starting point in position i the first in position j the second in position z the third
            for(int i=0; i<3; i++){
                for(int j=0; j<3; j++){
                    for(int z=0;z<3;z++){
                        if(i!=j && j!=z && z!=i){
                            //std::cout << i<< j << z << std::endl;
                            //create new route
                            //Rout->change_order(u, i , j, z);
                            bool f=Sol->route[vei]->delivery_before_pickup_2(u, i, j, z);

                           // std::cout<< "U  " << u <<  "   i -- > " << i <<  "  Type " << Sol->route[vei]->type[u+i] << " " << Sol->route[vei]->Ricid[u+i] <<"  j -- > " << j <<  "  Type " << Sol->route[vei]->type[u+j] << " " << Sol->route[vei]->Ricid[u+j] << "z -- > " << z <<  "  Type " << Sol->route[vei]->type[u+z] << " " << Sol->route[vei]->Ricid[u+z] << "  RESULT: " << f << std::endl;



                            if(f==true){
                                newdist=Sol->route[vei]->calculatedist_2(D,i,j,z,u);
                                //std::cout<< newdist << "<" << ordist << std::endl;
                                if(newdist<ordist){
                                    //Rout->calculate_earliest_latest(ric, D, MaxTimeRoute);
                                    bool f1=Sol->route[vei]->check_feas1(u, ric, D,i,j,z);


                                    if(f1==true){
                                        //Rout->calculate_capacity(ric, veic);
                                        bool f2=Sol->route[vei]->check_cap(u, veic, i, j, z, ric);


                                        if(f2==true){
                                            bool f3=Sol->route[vei]->EightStepEvaluation_2(D, ric, u, i, j , z, veic);


                                            if(f3==true){
                                                changed=true;
                                                //save i j z and u
                                                best_i=i;
                                                best_j=j;
                                                best_z=z;
                                                best_u=u;
                                                ordist=newdist;

                                                //std::cout<< "change" << i << j << z << u << std::endl;
                                            }


                                        }

                                    }

                                }

                            }

                        }

                    }

                }
            }

        }

        if(changed==true){
            //std::cout<< "Change start point in " << vei << " start " << best_u << " " << best_i << " " << best_j << " " << best_z << std::endl;

            Sol->route[vei]->change_order(best_u, best_i , best_j, best_z);
            Sol->route[vei]->calculate_earliest_latest(ric, D, veic);
            Sol->route[vei]->calculate_capacity(ric, veic);

            Sol->route[vei]->Eightstepevaluationscheme(D, ric, veic);

            Sol->route[vei]->totaldist=Sol->route[vei]->calculatedist(D);

            Sol->updatecost(numVeicoli);
            Sol->route[vei]->display_route();

           // return 0;

        }
    
    
    
    return Sol;
}

/*
Solution* Eliminate_Neighborhood(Solution* Sol, int numVeicoli, int** MatCompVei, RichiesteServizio* ric, double** MatTemp, int MaxTimeRoute, Veicoli* veic, int MaxRideTime, double** D){
    
    
    int vei;
    int i,j, l, g;
    //randomly chose a route...
    
    
    Solution* Solaux=NULL;
    
    Solaux=new Solution();
    
    Solaux->CopySolution(Sol);
    
    do{
        vei=rand()%(numVeicoli+Solaux->numAddRoutes);
        //take out from this route all the requests and insert them
    }while(Solaux->route[vei]->numRicRoute==0);
    
    //take out all the requests from the route
    int* E;
    
    E=new int[Solaux->route[vei]->numRicRoute];
    
    int numofreq=0;
    numofreq=Solaux->route[vei]->numRicRoute;
    E=Solaux->route[vei]->count_request(E);
    
    
    E=random_order2(numofreq, E);
    
    
    Solaux->delete_route(vei,numVeicoli);
    int bestpl, bestdl, bestv;
    //take all the request
    bool inserted,changed;
    changed=true;
    double bestdist;
    i=0;
    while(i<numofreq){
        //	std::cout<< "request-> "<< E[i] << "\n";
        inserted=false;
        //Create a Route to insert everything on it.
        
        bool feas1, feas2, feas3, compatibility,feasride, feas4, feas5, feas6;
        //double Bestdisteffect=DBL_MAX;
        bestdist=DBL_MAX;
        for(j=0;j<numVeicoli;j++){
            if(j!=vei){
                //std::cout<< "vehicle-> "<< j << "\n";
                // See if the vehicle is compatible with the request to be inserted.
                if(MatCompVei[E[i]][j]==1){
                    
                    compatibility=true;
                    //std::cout << "request " << E[i] << "compatible in vehicle " << j << endl;
                }else{
                    compatibility=false;
                }
                if(compatibility==true){
                    Route* rout=NULL;
                    rout=new Route();
                    rout=Solaux->copyvehicleinroute(rout, j);
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
                        
                        if(feas1==true && feas2==true){
                            //		std::cout << "pickup of request " << E[i] << "compatible in position " << l << endl;
                            //
                            
                            for (g=l; g<originallength-1; g++){
                                
                                rout->insert_delivery(E[i], g, ric);
                                
                                rout->calculate_earliest_latest(ric, MatTemp, MaxTimeRoute);
                                
                                rout->calculate_capacity(ric, veic);
                                
                                rout->calculate_earliest_latest(ric, MatTemp, MaxTimeRoute);
                                
                                
                                feasride=rout->ridetime_check( g,  ric,  MaxRideTime);
                                
                                if(feasride==true){//std::cout <<"true" << endl;
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
                                            
                                            if(feas3==true){
                                                //calculate 8 step evaluation scheme
                                                //	std::cout << "check delivery tw feasibility OK " << endl;
                                                
                                                gg=rout->Eightstepevaluationscheme(MatTemp, ric,  MaxRideTime,  MaxTimeRoute);
                                                //	rout->display_route();
                                                if(gg==true){
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
                                        if(feas4==true && feas5==true){
                                            disteffect=rout->effect_ondistance(E[i], l, g, D,  ric);
                                            //std::cout <<  disteffect << " < " << MINI[0] << "\n";
                                            //	std:: cout << "length: " << rout->length << "l : " << l << "g : " << g << endl;
                                            //	std::cout <<  disteffect << " < " << bestdist << "\n";
                                            if(disteffect<bestdist){
                                                //calculate 8 step evaluation scheme
                                                feas6=rout->check_feasibillity_D_tw1(E[i],l,ric,MatTemp);
                                                //	std::cout << feas3 << "  feas3 " << endl;
                                                if(feas6==true){
                                                    //calculate 8 step evaluation scheme
                                                    //	std::cout << "check delivery tw feasibility OK " << endl;
                                                    
                                                    gg=rout->Eightstepevaluationscheme(MatTemp, ric,  MaxRideTime,  MaxTimeRoute);
                                                    //	rout->display_route();
                                                    if(gg==true){
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
                    
                    delete rout;
                    
                }
            }
        }
        if(inserted==true){
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
            i=numofreq;
            changed=false;
            
            
        }
        
        
        i+=1;
        
    }
    if(changed==true){
        
        Sol->copyexistingsol(Solaux);
    }
    
    
    delete Solaux;
    delete[] E;
    
    
    return Sol;
}


*/

/*
Solution* Relocate_Neighborhood(Solution* Sol, int numVeicoli, int** MatCompVei, RichiesteServizio* ric, double** MatTemp, int MaxTimeRoute, Veicoli* veic, int MaxRideTime, double** D){
    
    int vei;
    int i, j, l, g;
    
    do{
        vei=rand()%(numVeicoli+Sol->numAddRoutes);
        //take out from this route all the requests and insert them
    }while(Sol->route[vei]->numRicRoute==0);
    
    //  std::cout << Sol->route[vei]->length << std::endl;
    //Sol->route[vei]->numRicRoute=((Sol->route[vei]->length)-2)/2;
    //take out all the requests from the route
    int* E;
    
    int numofreq=0;
    numofreq=Sol->route[vei]->numRicRoute;
    E=new int[numofreq];
    
    E=Sol->route[vei]->count_request(E);
    //std::cout<< numofreq << std:: endl;
    //for(i=0; i<numofreq; i++){
    //    std::cout<< E[i] << " ";
    //}std::cout<< std:: endl;
    //   std::cout<< std:: endl;
    
    //Sol->delete_route(vei,numVeicoli);
    
    
    int bestpl=0, bestdl=0, bestv=0;
    
    bool inserted,changed;
    changed=true;
    inserted=false;
    double bestdist, bestreq=0;
    double or_dist, new_dist, dif;
    //Sol->route[vei]->display_route();
    Route* rout1=NULL;
    rout1=new Route();
    
    or_dist=Sol->route[vei]->totaldist;
    Route* rout=NULL;
    rout=new Route();
    
    //original distance normal no crear la ruta
    
    for(i=0;i<numofreq;i++){
        // std::cout << " attempt no: " ;
        //std::cout << i << std::endl;
        
        rout1->clear_route();
        
        
        
        rout1=Sol->copyvehicleinroute(rout1, vei);
        // rout1->display_route();
        
        // or_dist=rout1->calculatedist(MatTemp);
        // std::cout<< E[i] << std::endl;
        int req=E[i];
        rout1->delete_req(req, MatTemp);
        new_dist=rout1->totaldist;
        dif=new_dist-or_dist;
        //calcular distancia, sin hacer quita y pon. quitando la request req, calcular la distancia. que sera new_dist. Asi me ahorro route y un copy y un todo de eso
        double new_dist2=Sol->route[vei]->new_dist_without_req(req,  MatTemp);
        
        //Aqui calcular la distancia que tendria si quitasemos req
        std::cout << new_dist << " "  << new_dist2 << std::endl;
        
        //calc dist withou
               
        
        //Create a Route to insert everything on it.
        
        bool feas1, feas2, feas3, compatibility,feasride, feas4, feas5, feas6;
        //double Bestdisteffect=DBL_MAX;
        bestdist=DBL_MAX;
        
        for(j=0;j<numVeicoli;j++){
            if(j!=vei){
                //std::cout<< "vehicle-> "<< j << "\n";
                // See if the vehicle is compatible with the request to be inserted.
                
                if(MatCompVei[req][j]>0){
                
                //lo que seria si insertasemos
                
                rout->clear_route();
                rout=Sol->copyvehicleinroute(rout, j);
                //
                int originallength=Sol->route[j]->length;
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
                    
                     //calcular si es feasible pero sin calcular las cosas....
                    
                    //	std:: cout << "inserting pickup on " << l << "\n";
                    feas1=rout->check_feasibility_P_tw(E[i],l, ric, MatTemp);
                    
                    if(feas1>0 && feas2>0){
                        //		std::cout << "pickup of request " << E[i] << "compatible in position " << l << endl;
                        //
                        
                        for (g=l; g<originallength-1; g++){
                             //calcular  los cachos de vectores desde l a g;.......
                            
                            rout->insert_delivery(E[i], g, ric);
                        
                            rout->calculate_earliest_latest(ric, MatTemp, MaxTimeRoute);
                            
                            feasride=rout->ridetime_check( g,  ric,  MaxRideTime);
                            //calcular si es feasible pero sin calcular las cosas....
                            
                            if(feasride>0){//std::cout <<"true" << endl;
                                double disteffect=0;
                                bool gg;
                                //	std:: cout << "length: " << rout->length << "l : " << l << "g : " << g << endl;
                              
                                    feas4=rout->check_feasibility_D_tw2( l,  g, ric, MatTemp);
                                    //	std::cout << feas3 << "  1feas3 " << endl;
                                    rout->calculate_capacity(ric, veic);
                                    feas5=rout->check_feasibility_capacity(l,g,veic);
                                    //	std::cout << feas3 << "  2feas3 " << endl;
                                    //Insert E{i} in route j in position l pickup g delivery
                                    //rout.insert_req_pos(E[i], l, g, ric);
                                    if(feas4>0 && feas5>0){
                                        disteffect=rout->effect_ondistance(E[i], l, g, D,  ric)+dif;
                                        //disteffect=disteffect+dif;
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
                                                    bestreq=E[i];
                                                    inserted=true;
                                                    
                                                }
                                            }
                                            
                                        }else{
                                            //Continue with the following insertion, This position is not considered
                                        }
                                    }
                                //}
                                
                            }
                            
                            rout->delete_delivery(g);
                            
                        }
                        
                        
                    }
                    rout->delete_pickup(l, MatTemp);
                    //std:: cout << "pickup deleted" << endl;
                    
                }
                
                
                }//close MatCompVei
            }
        }
        //std::cout<< numofreq << std:: endl;
        //std::cout<<E[i] << std:: endl;
        
        
    }
    
    //std::cout<<vei << std:: endl;
    if(inserted>0){
        //std::cout<< "true" << std::endl;
        // Sol->route[bestv]->display_route();
        Sol->route[bestv]->insert_req_pos(bestreq,bestpl,bestdl,ric, MatTemp);
        Sol->route[bestv]->calculate_capacity(ric,veic);
        // Sol->route[bestv]->display_route();
        Sol->route[vei]->delete_req(bestreq, MatTemp);
        
        Sol->route[bestv]->calculate_earliest_latest(ric, MatTemp, MaxTimeRoute);
        Sol->route[bestv]->Eightstepevaluationscheme(MatTemp,ric,MaxRideTime,MaxTimeRoute);
        Sol->route[vei]->calculate_capacity(ric,veic);
        Sol->route[vei]->calculate_earliest_latest(ric, MatTemp, MaxTimeRoute);
        Sol->route[vei]->Eightstepevaluationscheme(MatTemp,ric,MaxRideTime,MaxTimeRoute);
        Sol->route[bestv]->totaldist=Sol->route[bestv]->calculatedist(MatTemp);
        Sol->route[vei]->totaldist=Sol->route[vei]->calculatedist(MatTemp);
        //Solaux->copyrouteinsol(rout);
        //	std::cout << "request" << E[i] << "inserted in vehicle, " << bestv<< "\n";
        //Solaux->route[bestv]->display_route();
        Sol->updatecost(numVeicoli);
        Sol->updateusedvehicles(numVeicoli);
        //Sol->route[bestv]->display_route();
        //std::cin >> j;
        
    }else{
        
        
        
        
        
        
        std::cout<< "false" << std::endl;
    }
    delete rout;
    
    delete rout1;
    //
    delete[] E;
    return Sol;
}
*/




SolutionVRP* Relocate_NeighborhoodF(SolutionVRP* Sol, int numVeicoli, std::vector<std::vector<int>> &MatCompVei, std::vector<RichiesteServizio*> &ric, std::vector<std::vector<double>> & MatTemp, std::vector<Veicoli*> &veic, std::vector<std::vector<double>> & D, int numRichieste){
    
    int vei;
    int i, j, l, g;
    
    int* rid;
    rid=new int[2*numRichieste+2];
    int* typ;
    typ=new int[2*numRichieste+2];
    int* loc;
    loc=new int[2*numRichieste+2];
    double* earl;
    earl=new double[2*numRichieste+2];
    
    do{
        vei=rand()%(numVeicoli+Sol->numAddRoutes);
        //take out from this route all the requests and insert them
    }while(Sol->route[vei]->numRicRoute==0);
    
    //  std::cout << Sol->route[vei]->length << std::endl;
    //Sol->route[vei]->numRicRoute=((Sol->route[vei]->length)-2)/2;
    //take out all the requests from the route
    int* E;
    
    int numofreq=0;
    numofreq=Sol->route[vei]->numRicRoute;
    E=new int[numofreq];
    
    E=Sol->route[vei]->count_request(E);
    //std::cout<< numofreq << std:: endl;
    //for(i=0; i<numofreq; i++){
    //    std::cout<< E[i] << " ";
    //}std::cout<< std:: endl;
    //   std::cout<< std:: endl;
    
    //Sol->delete_route(vei,numVeicoli);
    
    
    int bestpl=0, bestdl=0, bestv=0;
    
    bool inserted,changed;
    changed=true;
    
    double bestdist, bestreq=0;
    double or_dist, new_dist, dif;
    //Sol->route[vei]->display_route();
    //Route* rout1=NULL;
    //rout1=new Route();
    
    or_dist=Sol->route[vei]->totaldist;
    //Route* rout=NULL;
    //rout=new Route();
    inserted=false;
    //original distance normal no crear la ruta
    
    for(i=0;i<numofreq;i++){
        // std::cout << " attempt no: " ;
        //std::cout << i << std::endl;
        
        //rout1->clear_route();
        
        
        
       // rout1=Sol->copyvehicleinroute(rout1, vei);
        // rout1->display_route();
        
        // or_dist=rout1->calculatedist(MatTemp);
        // std::cout<< E[i] << std::endl;
        int req=E[i];
       // rout1->delete_req(req, MatTemp);
        //new_dist=Sol->route[vei]->calc_cost_without_req();
        
       // new_dist=rout1->totaldist;
        
        //calcular distancia, sin hacer quita y pon. quitando la request req, calcular la distancia. que sera new_dist. Asi me ahorro route y un copy y un todo de eso
        
        
        
        new_dist=Sol->route[vei]->new_dist_without_req(req,  MatTemp);
        
        //Aqui calcular la distancia que tendria si quitasemos req
        //std::cout << new_dist << " "  << new_dist2 << std::endl;
        
        //calc dist withou
        dif=new_dist-or_dist;
        
        //Create a Route to insert everything on it.
        
        bool feas1, feas2,feasride, feas4, feas5, feas6;
        //double Bestdisteffect=DBL_MAX;
        bestdist=DBL_MAX;
        
        for(j=0;j<numVeicoli;j++){
            if(j!=vei){
                
                if(MatCompVei[req][j]>0){
                                      int originallength=Sol->route[j]->length;
                    for(l=0; l<originallength-1; l++){
                        
                        feas2=Sol->route[j]->capacity_P_feasibility( l, E[i],  veic, ric);
                        
                        feas1=Sol->route[j]->tw_P_feasibility(ric, MatTemp,  l, E[i]);
                        
                        if(feas1>0 && feas2>0){
                            //		std::cout << "pickup of request " << E[i] << "compatible in position " << l << endl;
                            //
                            
                            g=l;
                            // for (g=l; g<originallength-1; g++){
                            while(g<originallength-1){                                //calcular  los cachos de vectores desde l a g;.......
                                double disteffect=0;
                                bool gg;
                                
                                if(g==l){
                                    disteffect=0;
                                    disteffect=Sol->route[j]->effect_of_inserting_req_on_pos(E[i], l, g, D, ric)+dif;
                                    if(disteffect<bestdist){
                                        
                                        rid=Sol->route[j]->calcrid(rid, l, g, E[i]);
                                        typ=Sol->route[j]->calctyp(typ, l, g);
                                        loc=Sol->route[j]->calcloc(loc, E[i],ric, l, g);
                                        earl=Sol->route[j]->calcearl(earl, E[i], ric, MatTemp, l, g, loc, rid, typ);
                                        
                                        
                                        feas6=Sol->route[j]->check_feas_D_tw1(earl, E[i], ric, g, l, MatTemp, veic);
                                        //std::cout << "6" << feas6 << endl;
                                        if(feas6>0){
                                            //calculate 8 step evaluation scheme
                                            //	std::cout << "check delivery tw feasibility OK " << endl;
                                            
                                            //gg=Solaux->route[j]->Eightstepevaluationscheme(MatTemp, ric,  MaxRideTime,  g);
                                            gg=Sol->route[j]->EightStepEvaluation(MatTemp, ric, l, g, E[i], veic);
                                            //	rout->display_route();
                                            //  std::cout << "gg   " << gg << std::endl;
                                            
                                            
                                            if(gg>0){
                                                bestdist=disteffect;
                                                bestpl=l;
                                                bestdl=g;
                                                bestv=j;
                                                bestreq=E[i];
                                                inserted=true;
                                                
                                            }
                                        }

                                    }else{
                                        //Continue with the following insertion, This position is not considered
                                    }
                                    
                                    
                                
                                }else{
                                    
                                    
                                    rid=Sol->route[j]->calcrid(rid, l, g, E[i]);
                                    typ=Sol->route[j]->calctyp(typ, l, g);
                                    loc=Sol->route[j]->calcloc(loc, E[i],ric, l, g);
                                    earl=Sol->route[j]->calcearl(earl, E[i], ric, MatTemp, l, g, loc, rid, typ);
                                    //Solaux->route[j]->insert_delivery(E[i], g, ric);
                                    //Solaux->route[j]->calculate_earliest_latest(ric, MatTemp, MaxTimeRoute);
                                    //feasride=Solaux->route[j]->ridetime_check( g,  ric,  MaxRideTime);
                                    feasride=Sol->route[j]->ridetime_feas_D( g,  l, ric, E[i], MatTemp, earl);
                                
                               // rout->insert_delivery(E[i], g, ric);
                                
                                //rout->calculate_earliest_latest(ric, MatTemp, MaxTimeRoute);
                                
                                //feasride=rout->ridetime_check( g,  ric,  MaxRideTime);
                                //calcular si es feasible pero sin calcular las cosas....
                                
                                    if(feasride>0){//std::cout <<"true" << endl;
                                        disteffect=0;
                                        //feas4=Solaux->route[j]->check_feasibility_D_tw2( l,  g, ric, MatTemp);
                                        
                                        feas4=check_feas_D_tw2(earl,  typ, rid,  l,  g, ric);
                                    //	std:: cout << "length: " << rout->length << "l : " << l << "g : " << g << endl;
                                    
                                        //feas4=rout->check_feasibility_D_tw2( l,  g, ric, MatTemp);
                                        //	std::cout << feas3 << "  1feas3 " << endl;
                                        //rout->calculate_capacity(ric, veic);
                                        //feas5=rout->check_feasibility_capacity(l,g,veic);
                                        //	std::cout << feas3 << "  2feas3 " << endl;
                                        //Insert E{i} in route j in position l pickup g delivery
                                        //rout.insert_req_pos(E[i], l, g, ric);
                                        
                                        feas5=Sol->route[j]->check_cap_from(l, g, veic, E[i], ric);
                                        
                                        if(feas4>0 && feas5>0){
                                            disteffect=Sol->route[j]->effect_of_inserting_req_on_pos(E[i], l, g, D, ric)+dif;
                                            //disteffect=rout->effect_ondistance(E[i], l, g, D,  ric)+dif;
                                            //disteffect=disteffect+dif;
                                            //std::cout <<  disteffect << " < " << MINI[0] << "\n";
                                            //	std:: cout << "length: " << rout->length << "l : " << l << "g : " << g << endl;
                                            //	std::cout <<  disteffect << " < " << bestdist << "\n";
                                            if(disteffect<bestdist){
                                                //calculate 8 step evaluation scheme
                                                feas6=Sol->route[j]->check_feas_D_tw1(earl, E[i], ric, g, l, MatTemp, veic);
                                                //	std::cout << feas3 << "  feas3 " << endl;
                                                if(feas6>0){
                                                    //calculate 8 step evaluation scheme
                                                    //	std::cout << "check delivery tw feasibility OK " << endl;
                                                    
                                                    gg=Sol->route[j]->EightStepEvaluation(MatTemp, ric, l, g, E[i], veic);
                                                    //	rout->display_route();
                                                    if(gg>0){
                                                        bestdist=disteffect;
                                                        bestpl=l;
                                                        bestdl=g;
                                                        bestv=j;
                                                        bestreq=E[i];
                                                        inserted=true;
                                                        
                                                    }
                                                }
                                                
                                            }else{
                                                //Continue with the following insertion, This position is not considered
                                            }
                                        }else{
                                            
                                            g=originallength-1;
                                        }
                                    //}
                                    
                                    }else{
                                        g=originallength-1;
                                        
                                    }
                                
                                    
                                    
                                    
                               // rout->delete_delivery(g);
                                }
                                g+=1;
                            }//end while
                            
                            
                        }
                        //rout->delete_pickup(l, MatTemp);
                        //std:: cout << "pickup deleted" << endl;
                        
                    }
                    
                    
                }//close MatCompVei
            }
        }
        //std::cout<< numofreq << std:: endl;
        //std::cout<<E[i] << std:: endl;
        
        
    }
    
    
    
    //std::cout<<vei << std:: endl;
    if(inserted>0){
        //std::cout<< "true" << std::endl;
        // Sol->route[bestv]->display_route();
        Sol->route[bestv]->insert_req_pos(bestreq,bestpl,bestdl,ric, MatTemp);
        Sol->route[bestv]->calculate_capacity(ric,veic);
        // Sol->route[bestv]->display_route();
        Sol->route[vei]->delete_req(bestreq, MatTemp);
        
        Sol->route[bestv]->calculate_earliest_latest(ric, MatTemp, veic);
        Sol->route[bestv]->Eightstepevaluationscheme(MatTemp,ric,veic);
        Sol->route[vei]->calculate_capacity(ric,veic);
        Sol->route[vei]->calculate_earliest_latest(ric, MatTemp, veic);
        Sol->route[vei]->Eightstepevaluationscheme(MatTemp,ric,veic);
        Sol->route[bestv]->totaldist=Sol->route[bestv]->calculatedist(MatTemp);
        Sol->route[vei]->totaldist=Sol->route[vei]->calculatedist(MatTemp);
        //Solaux->copyrouteinsol(rout);
        //	std::cout << "request" << E[i] << "inserted in vehicle, " << bestv<< "\n";
        //Solaux->route[bestv]->display_route();
        Sol->updatecost(numVeicoli);
        Sol->updateusedvehicles(numVeicoli);
        //Sol->route[bestv]->display_route();
        //std::cin >> j;
        
    }else{
        
        
        //std::cout<< "false" << std::endl;
    }
    //delete rout;
    
    //delete rout1;
    //
    
    delete[] rid;
    delete[] typ;
    delete[] loc;
    delete[] earl;
    delete[] E;
    return Sol;
}

void Relocate_Neighborhood_2(SolutionVRP* Sol,int vei,int& pickup,int& delivery,int& bestRequest, int& bestVei, int numVeicoli, std::vector<std::vector<int>> &MatCompVei, std::vector<RichiesteServizio*> &ric, std::vector<std::vector<double>> & MatTemp, std::vector<Veicoli*> &veic, std::vector<std::vector<double>> & D, int numRichieste){

    //int vei;
    int i, j, l, g;

    int* rid;
    rid=new int[2*numRichieste+2];
    int* typ;
    typ=new int[2*numRichieste+2];
    int* loc;
    loc=new int[2*numRichieste+2];
    double* earl;
    earl=new double[2*numRichieste+2];

    //do{
      //  vei=rand()%(numVeicoli+Sol->numAddRoutes);
        //take out from this route all the requests and insert them
    //}while(Sol->route[vei]->numRicRoute==0);

    //  std::cout << Sol->route[vei]->length << std::endl;
    //Sol->route[vei]->numRicRoute=((Sol->route[vei]->length)-2)/2;
    //take out all the requests from the route
    int* E;

    int numofreq=0;
    numofreq=Sol->route[vei]->numRicRoute;
    E=new int[numofreq];

    E=Sol->route[vei]->count_request(E);
    //std::cout<< numofreq << std:: endl;
    //for(i=0; i<numofreq; i++){
    //    std::cout<< E[i] << " ";
    //}std::cout<< std:: endl;
    //   std::cout<< std:: endl;

    //Sol->delete_route(vei,numVeicoli);


    int bestpl=0, bestdl=0, bestv=0;

    bool inserted,changed;
    changed=true;

    double bestdist, bestreq=0;
    double or_dist, new_dist, dif;
    //Sol->route[vei]->display_route();
    //Route* rout1=NULL;
    //rout1=new Route();

    or_dist=Sol->route[vei]->totaldist;
    //Route* rout=NULL;
    //rout=new Route();
    inserted=false;
    //original distance normal no crear la ruta

    for(i=0;i<numofreq;i++){
        // std::cout << " attempt no: " ;
        //std::cout << i << std::endl;

        //rout1->clear_route();



       // rout1=Sol->copyvehicleinroute(rout1, vei);
        // rout1->display_route();

        // or_dist=rout1->calculatedist(MatTemp);
        // std::cout<< E[i] << std::endl;
        int req=E[i];
       // rout1->delete_req(req, MatTemp);
        //new_dist=Sol->route[vei]->calc_cost_without_req();

       // new_dist=rout1->totaldist;

        //calcular distancia, sin hacer quita y pon. quitando la request req, calcular la distancia. que sera new_dist. Asi me ahorro route y un copy y un todo de eso



        new_dist=Sol->route[vei]->new_dist_without_req(req,  MatTemp);

        //Aqui calcular la distancia que tendria si quitasemos req
        //std::cout << new_dist << " "  << new_dist2 << std::endl;

        //calc dist withou
        dif=new_dist-or_dist;

        //Create a Route to insert everything on it.

        bool feas1, feas2,feasride, feas4, feas5, feas6;
        //double Bestdisteffect=DBL_MAX;
        bestdist=DBL_MAX;

        for(j=0;j<numVeicoli;j++){
            if(j!=vei){

                if(MatCompVei[req][j]>0){
                                      int originallength=Sol->route[j]->length;
                    for(l=0; l<originallength-1; l++){

                        feas2=Sol->route[j]->capacity_P_feasibility( l, E[i],  veic, ric);

                        feas1=Sol->route[j]->tw_P_feasibility(ric, MatTemp,  l, E[i]);

                        if(feas1>0 && feas2>0){
                            //		std::cout << "pickup of request " << E[i] << "compatible in position " << l << endl;
                            //

                            g=l;
                            // for (g=l; g<originallength-1; g++){
                            while(g<originallength-1){                                //calcular  los cachos de vectores desde l a g;.......
                                double disteffect=0;
                                bool gg;

                                if(g==l){
                                    disteffect=0;
                                    disteffect=Sol->route[j]->effect_of_inserting_req_on_pos(E[i], l, g, D, ric)+dif;
                                    if(disteffect<bestdist){

                                        rid=Sol->route[j]->calcrid(rid, l, g, E[i]);
                                        typ=Sol->route[j]->calctyp(typ, l, g);
                                        loc=Sol->route[j]->calcloc(loc, E[i],ric, l, g);
                                        earl=Sol->route[j]->calcearl(earl, E[i], ric, MatTemp, l, g, loc, rid, typ);


                                        feas6=Sol->route[j]->check_feas_D_tw1(earl, E[i], ric, g, l, MatTemp, veic);
                                        //std::cout << "6" << feas6 << endl;
                                        if(feas6>0){
                                            //calculate 8 step evaluation scheme
                                            //	std::cout << "check delivery tw feasibility OK " << endl;

                                            //gg=Solaux->route[j]->Eightstepevaluationscheme(MatTemp, ric,  MaxRideTime,  g);
                                            gg=Sol->route[j]->EightStepEvaluation(MatTemp, ric, l, g, E[i], veic);
                                            //	rout->display_route();
                                            //  std::cout << "gg   " << gg << std::endl;


                                            if(gg>0){
                                                bestdist=disteffect;
                                                bestpl=l;
                                                bestdl=g;
                                                bestv=j;
                                                bestreq=E[i];
                                                inserted=true;

                                            }
                                        }

                                    }else{
                                        //Continue with the following insertion, This position is not considered
                                    }



                                }else{


                                    rid=Sol->route[j]->calcrid(rid, l, g, E[i]);
                                    typ=Sol->route[j]->calctyp(typ, l, g);
                                    loc=Sol->route[j]->calcloc(loc, E[i],ric, l, g);
                                    earl=Sol->route[j]->calcearl(earl, E[i], ric, MatTemp, l, g, loc, rid, typ);
                                    //Solaux->route[j]->insert_delivery(E[i], g, ric);
                                    //Solaux->route[j]->calculate_earliest_latest(ric, MatTemp, MaxTimeRoute);
                                    //feasride=Solaux->route[j]->ridetime_check( g,  ric,  MaxRideTime);
                                    feasride=Sol->route[j]->ridetime_feas_D( g,  l, ric, E[i], MatTemp, earl);

                               // rout->insert_delivery(E[i], g, ric);

                                //rout->calculate_earliest_latest(ric, MatTemp, MaxTimeRoute);

                                //feasride=rout->ridetime_check( g,  ric,  MaxRideTime);
                                //calcular si es feasible pero sin calcular las cosas....

                                    if(feasride>0){//std::cout <<"true" << endl;
                                        disteffect=0;
                                        //feas4=Solaux->route[j]->check_feasibility_D_tw2( l,  g, ric, MatTemp);

                                        feas4=check_feas_D_tw2(earl,  typ, rid,  l,  g, ric);
                                    //	std:: cout << "length: " << rout->length << "l : " << l << "g : " << g << endl;

                                        //feas4=rout->check_feasibility_D_tw2( l,  g, ric, MatTemp);
                                        //	std::cout << feas3 << "  1feas3 " << endl;
                                        //rout->calculate_capacity(ric, veic);
                                        //feas5=rout->check_feasibility_capacity(l,g,veic);
                                        //	std::cout << feas3 << "  2feas3 " << endl;
                                        //Insert E{i} in route j in position l pickup g delivery
                                        //rout.insert_req_pos(E[i], l, g, ric);

                                        feas5=Sol->route[j]->check_cap_from(l, g, veic, E[i], ric);

                                        if(feas4>0 && feas5>0){
                                            disteffect=Sol->route[j]->effect_of_inserting_req_on_pos(E[i], l, g, D, ric)+dif;
                                            //disteffect=rout->effect_ondistance(E[i], l, g, D,  ric)+dif;
                                            //disteffect=disteffect+dif;
                                            //std::cout <<  disteffect << " < " << MINI[0] << "\n";
                                            //	std:: cout << "length: " << rout->length << "l : " << l << "g : " << g << endl;
                                            //	std::cout <<  disteffect << " < " << bestdist << "\n";
                                            if(disteffect<bestdist){
                                                //calculate 8 step evaluation scheme
                                                feas6=Sol->route[j]->check_feas_D_tw1(earl, E[i], ric, g, l, MatTemp, veic);
                                                //	std::cout << feas3 << "  feas3 " << endl;
                                                if(feas6>0){
                                                    //calculate 8 step evaluation scheme
                                                    //	std::cout << "check delivery tw feasibility OK " << endl;

                                                    gg=Sol->route[j]->EightStepEvaluation(MatTemp, ric, l, g, E[i], veic);
                                                    //	rout->display_route();
                                                    if(gg>0){
                                                        bestdist=disteffect;
                                                        bestpl=l;
                                                        bestdl=g;
                                                        bestv=j;
                                                        bestreq=E[i];
                                                        inserted=true;

                                                    }
                                                }

                                            }else{
                                                //Continue with the following insertion, This position is not considered
                                            }
                                        }else{

                                            g=originallength-1;
                                        }
                                    //}

                                    }else{
                                        g=originallength-1;

                                    }




                               // rout->delete_delivery(g);
                                }
                                g+=1;
                            }//end while


                        }
                        //rout->delete_pickup(l, MatTemp);
                        //std:: cout << "pickup deleted" << endl;

                    }


                }//close MatCompVei
            }
        }
        //std::cout<< numofreq << std:: endl;
        //std::cout<<E[i] << std:: endl;


    }



    //std::cout<<vei << std:: endl;
    if(inserted>0){
        //std::cout<< "true" << std::endl;
        // Sol->route[bestv]->display_route();
        pickup = Sol->route[vei]->find_pickup(bestreq);
        delivery = Sol->route[vei]->find_delivery(bestreq);
        bestRequest = bestreq;
        bestVei = bestv;
        Sol->route[bestv]->insert_req_pos(bestreq,bestpl,bestdl,ric, MatTemp);
        Sol->route[bestv]->calculate_capacity(ric,veic);
        // Sol->route[bestv]->display_route();
        Sol->route[vei]->delete_req(bestreq, MatTemp);

        Sol->route[bestv]->calculate_earliest_latest(ric, MatTemp, veic);
        Sol->route[bestv]->Eightstepevaluationscheme(MatTemp,ric,veic);
        Sol->route[vei]->calculate_capacity(ric,veic);
        Sol->route[vei]->calculate_earliest_latest(ric, MatTemp, veic);
        Sol->route[vei]->Eightstepevaluationscheme(MatTemp,ric,veic);
        Sol->route[bestv]->totaldist=Sol->route[bestv]->calculatedist(MatTemp);
        Sol->route[vei]->totaldist=Sol->route[vei]->calculatedist(MatTemp);
        //Solaux->copyrouteinsol(rout);
        //	std::cout << "request" << E[i] << "inserted in vehicle, " << bestv<< "\n";
        //Solaux->route[bestv]->display_route();
        Sol->updatecost(numVeicoli);
        Sol->updateusedvehicles(numVeicoli);
        //Sol->route[bestv]->display_route();
        //std::cin >> j;

    }else{


        //std::cout<< "false" << std::endl;
    }
    //delete rout;

    //delete rout1;
    //

    delete[] rid;
    delete[] typ;
    delete[] loc;
    delete[] earl;
    delete[] E;
   // return Sol;
}


SolutionVRP* Eliminate_NeighborhoodF(SolutionVRP* Sol, int numVeicoli, std::vector<std::vector<int>> & MatCompVei, std::vector<RichiesteServizio*> &ric, std::vector<std::vector<double>> & MatTemp, std::vector<Veicoli*> & veic, std::vector<std::vector<double>> & D, int numRichieste){
    
    
    int vei;
    int i,j, l, g;
    //randomly chose a route...
    
    int* rid;
    rid=new int[2*numRichieste+2];
    int* typ;
    typ=new int[2*numRichieste+2];
    int* loc;
    loc=new int[2*numRichieste+2];
    double* earl;
    earl=new double[2*numRichieste+2];
    SolutionVRP* Solaux=NULL;
    
    Solaux=new SolutionVRP();
    
    Solaux->CopySolution(Sol);
    
    do{
        vei=rand()%(numVeicoli+Solaux->numAddRoutes);
        //take out from this route all the requests and insert them
    }while(Solaux->route[vei]->numRicRoute==0);
    
    //take out all the requests from the route
    int* E;
    
    E=new int[Solaux->route[vei]->numRicRoute];
    
    int numofreq=0;
    numofreq=Solaux->route[vei]->numRicRoute;
    E=Solaux->route[vei]->count_request(E);
    
    
    E=random_order2(numofreq, E);
    
    
    
    int bestpl, bestdl, bestv;
    //take all the request
    bool inserted,changed;
    changed=true;
    double bestdist;
    i=0;
    while(i<numofreq){
       
        inserted=false;
        
        
        bool feas1, feas2, feas3, compatibility,feasride, feas4, feas5, feas6;
    
        bestdist=DBL_MAX;
        for(j=0;j<numVeicoli;j++){
            if(j!=vei){
                
                if(MatCompVei[E[i]][j]==1){
                    
                    compatibility=true;
                    
                }else{
                    compatibility=false;
                }
                if(compatibility==true){
                    int originallength=Solaux->route[j]->length;
                    for(l=0; l<originallength-1; l++){
                        feas2=Solaux->route[j]->capacity_P_feasibility( l, E[i],  veic, ric);
                        feas1=Solaux->route[j]->tw_P_feasibility(ric, MatTemp,  l, E[i]);
                  
                        if(feas1>0 && feas2>0){
                            
                            
                            //for (g=l; g<originallength-1; g++){
                            g=l;
                            while(g<originallength-1){
                                double disteffect=0;
                                bool gg;
                                if(g==l){
                                    disteffect=0;
                                    disteffect=Solaux->route[j]->effect_of_inserting_req_on_pos(E[i], l, g, D, ric);
                                    
                                    if(disteffect<bestdist){
                                        rid=Solaux->route[j]->calcrid(rid, l, g, E[i]);
                                        typ=Solaux->route[j]->calctyp(typ, l, g);
                                        loc=Solaux->route[j]->calcloc(loc, E[i],ric, l, g);
                                        earl=Solaux->route[j]->calcearl(earl, E[i], ric, MatTemp, l, g, loc, rid, typ);
                                        
                                        
                                        feas6=Sol->route[j]->check_feas_D_tw1(earl, E[i], ric, g, l, MatTemp, veic);
                                        
                                        if(feas6>0){
                                            gg=Solaux->route[j]->EightStepEvaluation(MatTemp, ric, l, g, E[i], veic);
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
                                    rid=Solaux->route[j]->calcrid(rid, l, g, E[i]);
                                    typ=Solaux->route[j]->calctyp(typ, l, g);
                                    loc=Solaux->route[j]->calcloc(loc, E[i],ric, l, g);
                                    earl=Solaux->route[j]->calcearl(earl, E[i], ric, MatTemp, l, g, loc, rid, typ);

                                    feasride=Solaux->route[j]->ridetime_feas_D( g,  l, ric, E[i],MatTemp, earl);
                                    if(feasride>0){
                                        disteffect=0;
                                        feas4=check_feas_D_tw2(earl,  typ, rid,  l,  g, ric);
                                        feas5=Solaux->route[j]->check_cap_from(l, g, veic, E[i], ric);
                                        if(feas4>0 && feas5>0){
                                            disteffect=Solaux->route[j]->effect_of_inserting_req_on_pos(E[i], l, g, D, ric);
                                            if(disteffect<bestdist){
                                                feas6=Solaux->route[j]->check_feas_D_tw1(earl, E[i], ric, g, l, MatTemp, veic);
                                                if(feas6>0){
                                                    gg=Solaux->route[j]->EightStepEvaluation(MatTemp, ric, l, g, E[i], veic);
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
                                            g=originallength-1;

                                        }
                                    }else{
                                        g=originallength-1;
                                    }
                                }
                                g=+1;
                            }//end while
                        }
                    }
                }
            }
        }
        
        if(inserted==true){
            Solaux->route[bestv]->insert_req_pos(E[i],bestpl,bestdl,ric, MatTemp);
            Solaux->route[bestv]->calculate_capacity(ric,veic);
            
            //Solaux->route[MINI[3]]->display_route();
            Solaux->route[bestv]->calculate_earliest_latest(ric, MatTemp, veic);
            Solaux->route[bestv]->Eightstepevaluationscheme(MatTemp,ric,veic);
            //Solaux->copyrouteinsol(rout);
            //	std::cout << "request" << E[i] << "inserted in vehicle, " << bestv<< "\n";
            //Solaux->route[bestv]->display_route();
            Solaux->updatecost(numVeicoli);
            //Solaux->updateusedvehicles(numVeicoli);
            
            
        }else{
            //std::cout<< "false" << std::endl;
            i=numofreq;
            changed=false;
            
            
        }
        
        
        i+=1;
        
        
        
    }
    if(changed==true){
        Solaux->delete_route(vei,numVeicoli);
        
        Sol->copyexistingsol(Solaux);
    }
                                                    
    
    delete[] rid;
    delete[] typ;
    delete[] loc;
    delete[] earl;
    delete Solaux;
    delete[] E;
    

    return Sol;
}

emili::Solution* RelocateNeighborhood::computeStep(emili::Solution* step)
{
    if(vei==(num_r-1))
    {
        return nullptr;
    }
    else
    {
        SolutionVRP* sol = (SolutionVRP*) step;
        pickup = -1;
        Relocate_Neighborhood_2(sol,vei,pickup,delivery,bestRequest,bestVei,inst.numVeicoli0,inst.MCV,inst.rc,inst.T,inst.vec,inst.Dist,inst.numRichieste0);
        vei++;
    }
    return step;
}

void RelocateNeighborhood::reverseLastMove(emili::Solution* step)
{
    if(pickup>0)
    {
        //reverse move
        SolutionVRP* Sol = (SolutionVRP*) step;
        Sol->route[vei]->insert_req_pos(bestRequest,pickup,delivery,inst.rc, inst.T);
        Sol->route[vei]->calculate_capacity(inst.rc,inst.vec);
        // Sol->route[bestv]->display_route();
        Sol->route[bestVei]->delete_req(bestRequest, inst.T);

        Sol->route[vei]->calculate_earliest_latest(inst.rc, inst.T, inst.vec);
        Sol->route[vei]->Eightstepevaluationscheme(inst.T,inst.rc,inst.vec);
        Sol->route[bestVei]->calculate_capacity(inst.rc,inst.vec);
        Sol->route[bestVei]->calculate_earliest_latest(inst.rc, inst.T, inst.vec);
        Sol->route[bestVei]->Eightstepevaluationscheme(inst.T,inst.rc,inst.vec);
        Sol->route[vei]->totaldist=Sol->route[vei]->calculatedist(inst.Dist);
        Sol->route[bestVei]->totaldist=Sol->route[bestVei]->calculatedist(inst.Dist);
        //Solaux->copyrouteinsol(rout);
        //	std::cout << "request" << E[i] << "inserted in vehicle, " << bestv<< "\n";
        //Solaux->route[bestv]->display_route();
        Sol->updatecost(inst.numVeicoli0);
        Sol->updateusedvehicles(inst.numVeicoli0);

    }

}

RelocateNeighborhood::NeighborhoodIterator RelocateNeighborhood::begin(emili::Solution* base)
{
    SolutionVRP* sol = (SolutionVRP*) base;
    vei = 0;
    num_r = sol->numRoutes+sol->numAddRoutes;
    pickup = -1;
    delivery = -1;
    bestRequest = -1;
    bestVei = -1;
    return NeighborhoodIterator(this,base);
}

void RelocateNeighborhood::reset()
{

}

emili::Solution* RelocateNeighborhood::random(emili::Solution* currentSolution)
{
    SolutionVRP* sol = (SolutionVRP*) currentSolution;
    SolutionVRP* nsol = new SolutionVRP();
    nsol->CopySolution(sol);
    int random_vei = emili::generateRandomNumber()%(sol->numRoutes+sol->numAddRoutes);
    Relocate_Neighborhood_2(nsol,random_vei,pickup,delivery,bestRequest,bestVei,inst.numVeicoli0,inst.MCV,inst.rc,inst.T,inst.vec,inst.Dist,inst.numRichieste0);
    return nsol;
}

