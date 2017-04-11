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
                                        
                                        Sol->route[j]->calcrid(rid, l, g, E[i]);
                                        Sol->route[j]->calctyp(typ, l, g);
                                        Sol->route[j]->calcloc(loc, E[i],ric, l, g);
                                        Sol->route[j]->calcearl(earl, E[i], ric, MatTemp, l, g, loc, rid, typ);
                                        
                                        
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
                                    
                                    
                                    Sol->route[j]->calcrid(rid, l, g, E[i]);
                                    Sol->route[j]->calctyp(typ, l, g);
                                    Sol->route[j]->calcloc(loc, E[i],ric, l, g);
                                    Sol->route[j]->calcearl(earl, E[i], ric, MatTemp, l, g, loc, rid, typ);
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
                                        Solaux->route[j]->calcrid(rid, l, g, E[i]);
                                        Solaux->route[j]->calctyp(typ, l, g);
                                        Solaux->route[j]->calcloc(loc, E[i],ric, l, g);
                                        Solaux->route[j]->calcearl(earl, E[i], ric, MatTemp, l, g, loc, rid, typ);
                                        
                                        
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
                                    Solaux->route[j]->calcrid(rid, l, g, E[i]);
                                    Solaux->route[j]->calctyp(typ, l, g);
                                    Solaux->route[j]->calcloc(loc, E[i],ric, l, g);
                                    Solaux->route[j]->calcearl(earl, E[i], ric, MatTemp, l, g, loc, rid, typ);

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
        
        Sol->copyexistingsol(Solaux, vei);
    }
                                                    
    
    delete[] rid;
    delete[] typ;
    delete[] loc;
    delete[] earl;
    delete Solaux;
    delete[] E;
    

    return Sol;
}

//Relocate Neighborhood for EMILI
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
   // std::cout << "VEI ::  "<< vei << std::endl;
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
                    if(Sol->route[j] == nullptr || Sol->route.size() < j )
                    {
                        std::cout << "J : " << j << "\n";
                        std::cout << numVeicoli << " NumADDR: " << Sol->numAddRoutes << " NumR: " << Sol->numRoutes << " Size: " << Sol->route.size() << "\n";

                    }
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

                                        Sol->route[j]->calcrid(rid, l, g, E[i]);
                                        Sol->route[j]->calctyp(typ, l, g);
                                        Sol->route[j]->calcloc(loc, E[i],ric, l, g);
                                        Sol->route[j]->calcearl(earl, E[i], ric, MatTemp, l, g, loc, rid, typ);


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


                                    Sol->route[j]->calcrid(rid, l, g, E[i]);
                                    Sol->route[j]->calctyp(typ, l, g);
                                    Sol->route[j]->calcloc(loc, E[i],ric, l, g);
                                    Sol->route[j]->calcearl(earl, E[i], ric, MatTemp, l, g, loc, rid, typ);
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
        if(Sol->route[vei]->numRicRoute>0){
        Sol->route[vei]->calculate_capacity(ric,veic);
        Sol->route[vei]->calculate_earliest_latest(ric, MatTemp, veic);
        Sol->route[vei]->Eightstepevaluationscheme(MatTemp,ric,veic);}
        Sol->route[bestv]->totaldist=Sol->route[bestv]->calculatedist(MatTemp);
       // std::cout<<"BEFORE " << Sol->route[vei]->totaldist << std::endl;
        Sol->route[vei]->totaldist=Sol->route[vei]->calculatedist(MatTemp);


        //std::cout<<Sol->route[vei]->totaldist << std::endl; //Solaux->copyrouteinsol(rout);
        //	std::cout << "request" << E[i] << "inserted in vehicle, " << bestv<< "\n";
        //Solaux->route[bestv]->display_route();
        //std::cout <<"REL BEFORE " << sol->getSolutionValue() << std::endl;
        Sol->updatecost(numVeicoli);
        //std::cout << Sol->route[0]->totaldist << Sol->route[0]->calculatedist(MatTemp) << std::endl;
        //std::cout << Sol->route[1]->totaldist << Sol->route[1]->calculatedist(MatTemp) << std::endl;
        //std::cout << Sol->route[2]->totaldist << Sol->route[2]->calculatedist(MatTemp) << std::endl;

        //std::cout <<"REL etwwtew " << Sol->getSolutionValue() << std::endl;
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

emili::Solution* RelocateNeighborhood::computeStep(emili::Solution* step)
{
    if(vei>=(num_r-1))
    {
        return nullptr;
    }
    else
    {
        SolutionVRP* sol = (SolutionVRP*) step;
        pickup = -1;
        while(sol->route[vei]->numRicRoute==0 && vei<num_r)
            vei++;
        //std::cout << sol->numAddRoutes << std::endl;
        //std::cout <<"REL BEFORE " << sol->getSolutionValue() << std::endl;
        //sol->DisplaySolution();

        Relocate_Neighborhood_2(sol,vei,pickup,delivery,bestRequest,bestVei,inst.numVeicoli0,inst.MCV,inst.rc,inst.T,inst.vec,inst.Dist,inst.numRichieste0);

        //std::cout << vei << " " << pickup << " " << delivery << " " << bestRequest << " " << bestVei << std::endl;
        //std::cout << sol->numAddRoutes << std::endl;
        //std::cout << sol->route[0]->calculatedist(inst.Dist) << std::endl;
        //std::cout << sol->route[1]->calculatedist(inst.Dist) << std::endl;
        //std::cout << sol->route[2]->calculatedist(inst.Dist) << std::endl;
        //std::cout <<"REL " << sol->getSolutionValue() << std::endl;
        //sol->DisplaySolution();
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



        Sol->route[vei-1]->insert_request_on(bestRequest,pickup,delivery,inst.rc, inst.T);

        Sol->route[vei-1]->calculate_capacity(inst.rc,inst.vec);
        // Sol->route[bestv]->display_route();
        Sol->route[bestVei]->delete_req(bestRequest, inst.T);

        Sol->route[vei-1]->calculate_earliest_latest(inst.rc, inst.T, inst.vec);
        Sol->route[vei-1]->Eightstepevaluationscheme(inst.T,inst.rc,inst.vec);
        Sol->route[bestVei]->calculate_capacity(inst.rc,inst.vec);
        Sol->route[bestVei]->calculate_earliest_latest(inst.rc, inst.T, inst.vec);
        Sol->route[bestVei]->Eightstepevaluationscheme(inst.T,inst.rc,inst.vec);
        Sol->route[vei-1]->totaldist=Sol->route[vei-1]->calculatedist(inst.Dist);
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
    while(sol->route[random_vei]->numRicRoute==0){
        random_vei=emili::generateRandomNumber()%(sol->numRoutes+sol->numAddRoutes);
        //take out from this route all the requests and insert them
    }
    //if(random_vei==1){
      //  sol->route[random_vei]->display_route();

   // }
    //std::cout << "random vei " << random_vei << "\n";

    //int random_vei = emili::generateRandomNumber()%(sol->numRoutes+sol->numAddRoutes);
    Relocate_Neighborhood_2(nsol,random_vei,pickup,delivery,bestRequest,bestVei,inst.numVeicoli0,inst.MCV,inst.rc,inst.T,inst.vec,inst.Dist,inst.numRichieste0);
    return nsol;
}

//Two Opt Neighborhood for EMILI
void two_opt_2(SolutionVRP* Sol, int vei, int& bestpos1, int& bestv2, int& bestpos2, int numVeicoli, std::vector<std::vector<double>> & D, std::vector<Veicoli*> &veic, std::vector<RichiesteServizio*> &ric, int numRichieste){
    int i, j, k, l;
    bool found=false;
    int v1, v2, pos1, pos2;
    double disteffect;
    double d;
    bool f1, f2, f3, f4;
    int numveh=Sol->numAddRoutes+Sol->numRoutes;
    disteffect=DBL_MAX;

    i=vei;
    //for(i=0;i<numveh;i++){
        for(j=1;j<Sol->route[i]->length-2;j++){
            if(Sol->route[i]->cap1[j]+Sol->route[i]->cap2[j]+Sol->route[i]->cap3[j]+Sol->route[i]->cap4[j]<1){
                //std::cout << "after position " << j << std::endl;
                for(k=i+1;k<numveh;k++){
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

    //}
    if(found==true){

        bestpos1=pos1;
        bestpos2=pos2;
        bestv2=v2;
       // std::cout << "asfsfa" << std::endl;
        //std::cout << v1 << " " << pos1 << " " << bestv2 << " " << pos2 << std::endl;
        //std::cout << v1 << std::endl;
       // Sol->route[v1]->display_route();
        //std::cout << v2 << std::endl;
       // Sol->route[v2]->display_route();
        Sol->changetworoutes(v1,pos1, v2, pos2, ric, D, veic);
        //std::cout << "asfsfa" << std::endl;

    }
    //if inserted is OK then redo all the routes and try

    //return Sol;
}


emili::Solution* TwoOptNeighborhood::computeStep(emili::Solution* step)
{
    if(vei==(num_r-1))
    {
        return nullptr;
    }
    else
    {
        SolutionVRP* sol = (SolutionVRP*) step;
        bestv2 = -1;
        while(sol->route[vei]->numRicRoute==0 && vei<num_r)
            vei++;

        two_opt_2(sol, vei, bestpos1,  bestv2, bestpos2,inst.numVeicoli0, inst.Dist, inst.vec, inst.rc, inst.numRichieste0);
        //std::cout << sol->numAddRoutes << std::endl;
        //std::cout <<"TOPT " << sol->getSolutionValue() << std::endl;
        vei++;
    }
    return step;
}


void TwoOptNeighborhood::reverseLastMove(emili::Solution* step)
{
    if(bestv2>=0)
    {
        //reverse move
        SolutionVRP* Sol = (SolutionVRP*) step;
        //std::cout << "asfsfa" << std::endl;
        Sol->changetworoutes(vei-1, bestpos1, bestv2, bestpos2, inst.rc, inst.Dist, inst.vec);
        //std::cout << "asfsfa" << std::endl;
    }

}

TwoOptNeighborhood::NeighborhoodIterator TwoOptNeighborhood::begin(emili::Solution* base)
{
    SolutionVRP* sol = (SolutionVRP*) base;
    vei = 0;
    num_r = sol->numRoutes+sol->numAddRoutes;
    bestpos1 = -1;
    bestv2 = -1;
    bestpos2 = -1;

    return NeighborhoodIterator(this,base);
}

void TwoOptNeighborhood::reset()
{

}

emili::Solution* TwoOptNeighborhood::random(emili::Solution* currentSolution)
{
    SolutionVRP* sol = (SolutionVRP*) currentSolution;
    SolutionVRP* nsol = new SolutionVRP();
    nsol->CopySolution(sol);
    int random_vei = emili::generateRandomNumber()%(sol->numRoutes+sol->numAddRoutes);
    two_opt_2(nsol, random_vei, bestpos1,  bestv2, bestpos2,inst.numVeicoli0, inst.Dist, inst.vec, inst.rc, inst.numRichieste0);

    return nsol;
}

//Eliminate Neighborhood for EMILI

void Eliminate_NeighborhoodF_2(SolutionVRP* Sol, int vei, int numVeicoli, std::vector<std::vector<int>> & MatCompVei, std::vector<RichiesteServizio*> &ric, std::vector<std::vector<double>> & MatTemp, std::vector<Veicoli*> & veic, std::vector<std::vector<double>> & D, int numRichieste){


   // int vei;
    int i,j, l, g;
    //randomly chose a route...

    //int* rid;
    //rid=new int[2*numRichieste+2];
    //int* typ;
    //typ=new int[2*numRichieste+2];
    //int* loc;
    //loc=new int[2*numRichieste+2];
    //double* earl;
    //earl=new double[2*numRichieste+2];
    int loc[2*numRichieste+2];
    double earl[2*numRichieste+2];
    int typ[2*numRichieste+2];
    int rid[2*numRichieste+2];
    SolutionVRP* Solaux=NULL;

    Solaux=new SolutionVRP();

    Solaux->CopySolution(Sol);

   // do{
     //   vei=rand()%(numVeicoli+Solaux->numAddRoutes);
        //take out from this route all the requests and insert them
    //}while(Solaux->route[vei]->numRicRoute==0);

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
    //std::cout<<numofreq<<std::endl;
    while(i<numofreq){

        inserted=false;


        bool feas1, feas2, feas3, compatibility,feasride, feas4, feas5, feas6;

        bestdist=DBL_MAX;
        for(j=0;j<numVeicoli+Solaux->numAddRoutes;j++){

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
                                        Solaux->route[j]->calcrid(rid, l, g, E[i]);
                                        Solaux->route[j]->calctyp(typ, l, g);
                                        Solaux->route[j]->calcloc(loc, E[i],ric, l, g);
                                        Solaux->route[j]->calcearl(earl, E[i], ric, MatTemp, l, g, loc, rid, typ);


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
                                    Solaux->route[j]->calcrid(rid, l, g, E[i]);
                                    Solaux->route[j]->calctyp(typ, l, g);
                                    Solaux->route[j]->calcloc(loc, E[i],ric, l, g);
                                    Solaux->route[j]->calcearl(earl, E[i], ric, MatTemp, l, g, loc, rid, typ);

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
                                g++;
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



        Sol->copyexistingsol(Solaux, vei);
    }


   // delete[] rid;
    //delete[] typ;
    //delete[] loc;
    //delete[] earl;


    delete Solaux;
    delete[] E;


    //return Sol;
}



emili::Solution* EliminateNeighborhood::computeStep(emili::Solution* step)
{
    if(vei>=(num_r-1))
    {
        return nullptr;
    }
    else
    {
        SolutionVRP* sol = (SolutionVRP*) step;
        while(sol->route[vei]->numRicRoute==0 && vei<num_r)
            vei++;


        Eliminate_NeighborhoodF_2(sol, vei, inst.numVeicoli0,inst.MCV, inst.rc, inst.T, inst.vec, inst.Dist, inst.numRichieste0);

        //r_4_opt_2(sol,  vei,  best_start, best_first, best_second, best_third, inst.numVeicoli0, inst.Dist, inst.rc, inst.vec);

        vei++;
    }
    return step;
}


void EliminateNeighborhood::reverseLastMove(emili::Solution* step)
{
    SolutionVRP* sol = (SolutionVRP*) step;
    sol->CopySolution(base_solution);
}

EliminateNeighborhood::NeighborhoodIterator EliminateNeighborhood::begin(emili::Solution* base)
{
    SolutionVRP* sol = (SolutionVRP*) base;
    vei = 0;
    num_r = sol->numRoutes+sol->numAddRoutes;
    base_solution->CopySolution(sol);

    return NeighborhoodIterator(this,base);
}

void EliminateNeighborhood::reset()
{

}

emili::Solution* EliminateNeighborhood::random(emili::Solution* currentSolution)
{
    SolutionVRP* sol = (SolutionVRP*) currentSolution;
    SolutionVRP* nsol = new SolutionVRP();
    nsol->CopySolution(sol);
    int random_vei = emili::generateRandomNumber()%(sol->numRoutes+sol->numAddRoutes);
    while(sol->route[random_vei]->numRicRoute==0){
        random_vei=emili::generateRandomNumber()%(sol->numRoutes+sol->numAddRoutes);
        //take out from this route all the requests and insert them
    }
    //int random_vei = emili::generateRandomNumber()%(sol->numRoutes+sol->numAddRoutes);
    //r_4_opt_2(nsol,  random_vei,  best_start, best_first, best_second, best_third, inst.numVeicoli0, inst.Dist, inst.rc, inst.vec);
    Eliminate_NeighborhoodF_2(nsol, random_vei, inst.numVeicoli0,inst.MCV, inst.rc, inst.T, inst.vec, inst.Dist, inst.numRichieste0);

    return nsol;
}


//4OPT Neighborhood for EMILI


void r_4_opt_2(SolutionVRP* Sol, int vei, int& best_start, int& best_first, int& best_second, int& best_third, int numVeicoli, std::vector<std::vector<double>> & D, std::vector<RichiesteServizio*> &ric, std::vector<Veicoli*> &veic){


    double ordist, newdist;

    //do{
        //    vei=rand()%(numVeicoli+Sol->numAddRoutes);
            //take out from this route all the requests and insert them
     //   }while(Sol->route[vei]->numRicRoute==0);

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
            best_start=best_u;
            best_first=best_i;
            best_second=best_j;
            best_third=best_z;
            Sol->route[vei]->change_order(best_u, best_i , best_j, best_z);
            Sol->route[vei]->calculate_earliest_latest(ric, D, veic);
            Sol->route[vei]->calculate_capacity(ric, veic);

            Sol->route[vei]->Eightstepevaluationscheme(D, ric, veic);

            Sol->route[vei]->totaldist=Sol->route[vei]->calculatedist(D);

            Sol->updatecost(numVeicoli);
            //Sol->route[vei]->display_route();

           // return 0;

        }



    //return Sol;
}

emili::Solution* FourOptNeighborhood::computeStep(emili::Solution* step)
{
    if(vei>=(num_r-1))
    {
        return nullptr;
    }
    else
    {
        SolutionVRP* sol = (SolutionVRP*) step;
        best_start = -1;
        while(sol->route[vei]->numRicRoute==0 && vei<num_r)
            vei++;
        r_4_opt_2(sol,  vei,  best_start, best_first, best_second, best_third, inst.numVeicoli0, inst.Dist, inst.rc, inst.vec);
        //std::cout << sol->numAddRoutes << std::endl;
        //std::cout <<"FOPT " << sol->getSolutionValue() << std::endl;
        vei++;
    }
    return step;
}


void FourOptNeighborhood::reverseLastMove(emili::Solution* step)
{
    if(best_start>0)
    {
        //reverse move
        SolutionVRP* Sol = (SolutionVRP*) step;
        Sol->route[vei-1]->reverse_order(best_start, best_first, best_second, best_third);

        Sol->route[vei-1]->calculate_earliest_latest(inst.rc, inst.Dist, inst.vec);
        Sol->route[vei-1]->calculate_capacity(inst.rc, inst.vec);

        Sol->route[vei-1]->Eightstepevaluationscheme(inst.Dist, inst.rc, inst.vec);

        Sol->route[vei-1]->totaldist=Sol->route[vei-1]->calculatedist(inst.Dist);

        Sol->updatecost(inst.numVeicoli0);


    }

}

FourOptNeighborhood::NeighborhoodIterator FourOptNeighborhood::begin(emili::Solution* base)
{
    SolutionVRP* sol = (SolutionVRP*) base;
    vei = 0;
    num_r = sol->numRoutes+sol->numAddRoutes;
    best_start = -1;
    best_first = -1;
    best_second = -1;
    best_third = -1;

    return NeighborhoodIterator(this,base);
}

void FourOptNeighborhood::reset()
{

}

emili::Solution* FourOptNeighborhood::random(emili::Solution* currentSolution)
{
    SolutionVRP* sol = (SolutionVRP*) currentSolution;
    SolutionVRP* nsol = new SolutionVRP();
    nsol->CopySolution(sol);
    int random_vei = emili::generateRandomNumber()%(sol->numRoutes+sol->numAddRoutes);
    r_4_opt_2(nsol,  random_vei,  best_start, best_first, best_second, best_third, inst.numVeicoli0, inst.Dist, inst.rc, inst.vec);

    return nsol;
}


//EXCHANGE FOR EMILI

void Exchange_Neighborhood(SolutionVRP* Sol, int vei, int &best_v, int &best_r1, int& best_r2, int &or_pickup_pos_1, int &or_pickup_pos_2, int &or_delivery_pos_1, int &or_delivery_pos_2, int numVeicoli, std::vector<RichiesteServizio*> & ric, std::vector<std::vector<double>> & MatTemp,std::vector<Veicoli*> &veic, int numRichieste, std::vector<std::vector<int>> &MatCompVei, std::vector<std::vector<double>> & D){


int rid [2*numRichieste+2];
int typ [2*numRichieste+2];
int loc [2*numRichieste+2];
double earl [2*numRichieste+2];

   // int* rid;
   // rid=new int[2*numRichieste+2];
   // int* typ;
   // typ=new int[2*numRichieste+2];
  //  int* loc;
  //  loc=new int[2*numRichieste+2];
   // double* earl;
   // earl=new double[2*numRichieste+2];
    bool inserted=false;
    int p1best, d1best, p2best, d2best, bestr1, bestr2, besv, gbest, lbest, bestv;
    int* E;
    int g;
    int numofreq=0;
    numofreq=Sol->route[vei]->numRicRoute;
    E=new int[numofreq];


    int* F=new int[numRichieste];
    E=Sol->route[vei]->count_request(E);
    int rr1, rr2;

    bool feas2, feas1, feas6, gg;
    bool feas22, feas11, feas66, gg1, feasride2, feasride1, feas41, feas51;
    //First Swap
    double bestdist=DBL_MAX;
    for(int i=0;i<numofreq;i++){
        E=Sol->route[vei]->count_request(E);
        rr1=E[i];
        //std::cout << " I " << i << std::endl;
        int p2,d2, p1, d1;
        p1=Sol->route[vei]->find_pickup(rr1);
        d1=Sol->route[vei]->find_delivery(rr1);

        for(int j=0; j<numVeicoli; j++){
           // std::cout << "Veh " << j << std::endl;
            if(j!=vei){
                //std::cout << i << std::endl;
                //std::cout << E[i] << " " << rr1 << std::endl;
                //std::cout << " E" << std::endl;
                //for(int iii=0; iii<numofreq; iii++){
                  //  std::cout << E[iii] << " " ;
                //}std::cout << std::endl;
                if(MatCompVei[rr1][j]>0){
                    int numofreq2=0;
                    numofreq2=Sol->route[j]->numRicRoute;

                    F=Sol->route[j]->count_request(F);



                    for(int k=0; k<numofreq2; k++){
                        F=Sol->route[j]->count_request(F);
                        rr2=F[k];
                  //      std::cout << " F" << std::endl;
                       // std::cout << rr2 << std::endl;
                    //    for(int kk=0; kk<numofreq2; kk++){
                      //      std::cout << F[kk] << " " ;
                        //}std::cout << std::endl;
                        //std::cout << k << std::endl;

                        //std::cout << F[k] << " " << vei << std::endl;
                        if(MatCompVei[rr2][vei]>0){
                            p2=Sol->route[j]->find_pickup(rr2);
                            d2=Sol->route[j]->find_delivery(rr2);
                            //check the feasibility of E[i] in position d2 and p2
                            //if it is feasible, check where you can insert in the best position.
                            feas2=Sol->route[j]->capacity_feasibility2(p2, rr1, veic, ric);

                            feas1=Sol->route[j]->tw_P_feasibility_2(ric, MatTemp,  p2, rr1);

                            if(feas1>0 && feas2>0){

                                if(p2+1==d2){
                                    feas6=Sol->route[j]->check_feas_FirstRequest(p2, d2, MatTemp, ric, rr1);
                                    if(feas6>0){
                                        gg=Sol->route[j]->EightStep_2(p2,d2, rr1, MatTemp, ric ,veic);

                                        if(gg>0){//it is feasible to insert the first request on the position of the second one.
                                            //Now I am going to try and insert the second request in all the possible positions in the first one.
                                            for(int l=1; l<Sol->route[vei]->length-2; l++){

                                                feas22=Sol->route[vei]->capacity_P_feasibility_3(l, rr2,  veic, ric, p1, d1, rr1);

                                                feas11=Sol->route[vei]->tw_P_feasibility_3(ric, MatTemp,  l, rr2, p1, d1, loc, typ, rid);
                                                if(feas11>0 && feas22>0){

                                                    g=l+1;
                                                    while(g<Sol->route[vei]->length-1){
                                                        double disteffect=0;


                                                        Sol->route[vei]->calc_loc3(rr2, ric, p1, d1, l , g, loc);
                                                        //std::cout << "first vei: " << vei << " " << p1 << " " << d1 << " " << E[i] << "second vei: " << j << " " << p2 << " " << d2 << " " << F[k] << " IN " << l << " " << g << std::endl;
                                                       // for(int iii=0; iii<Sol->route[vei]->length; iii++){
                                                       //     std::cout << loc[iii] << std::endl;

                                                       // }

                                                        if(g==l+1){

                                                            disteffect=0;
                                                            disteffect=Sol->route[vei]->effect_of_inserting_req_on_pos_3(rr2, l, g, D, ric, p1, d1, loc);
                                                            if(disteffect<bestdist){
                                                                Sol->route[vei]->calc_typ3(rr2, p1, d1, l , g, typ);

                                                                 Sol->route[vei]->calc_rid3(rr2,  p1, d1, l , g, rid);
                                                              //  std::cout << "RID" << std::endl;
                                                              //  for(int iii=0; iii<Sol->route[vei]->length; iii++){
                                                              //      std::cout << rid[iii] << std::endl;

                                                              //  }

                                                             //   return 0;
                                                              //  earl=Sol->route[vei]->calc_earl3(F[k], ric, p1, d1, l , g);

                                                                feas66=Sol->route[vei]->check_feas_D_tw1_3(rr2,ric,l,g,MatTemp,veic, loc, rid, typ);
                                                                if(feas66>0){
                                                                    gg1=Sol->route[vei]->EightStep_3(MatTemp, ric, veic, loc, rid, typ);

                                                                    if(gg1>0){
                                                                        bestdist=disteffect;
                                                                        bestr1=rr1;
                                                                        bestr2=rr2;
                                                                        p1best=p1;
                                                                        d1best=d1;
                                                                        p2best=p2;
                                                                        d2best=d2;
                                                                        gbest=g;
                                                                        lbest=l;
                                                                        bestv=j;

                                                                        inserted=true;
                                                                    }
                                                                }
                                                            }
                                                        }else{
                                                            Sol->route[vei]->calc_typ3(rr2,  p1, d1, l , g, typ);
                                                            Sol->route[vei]->calc_rid3(rr2, p1, d1, l , g, rid);

                                                            Sol->route[vei]->calc_earl3( ric, typ, loc, rid, earl, MatTemp);

                                                            feasride2=Sol->route[vei]->calc_feas_rtime(l,g,ric,rr2,MatTemp, earl);
                                                            if(feasride2>0){
                                                                disteffect=0;
                                                                bool feas44, feas55;
                                                                feas44=Sol->route[vei]->check_feas_D_tw22(earl, typ, rid, l, g, ric);

                                                                feas55=Sol->route[vei]->check_cap_from22(l, g, veic, rr2, ric, rid, typ);

                                                                if(feas44>0 && feas55>0){
                                                                    disteffect=Sol->route[vei]->effect_of_inserting_req_on_pos_3(rr2, l, g, D, ric, p1, d1, loc);

                                                                    if(disteffect<bestdist){

                                                                        feas66=Sol->route[vei]->check_feas_D_tw1_3(rr2,ric,l,g,MatTemp,veic, loc, rid, typ);
                                                                        if(feas66>0){
                                                                            gg1=Sol->route[vei]->EightStep_3(MatTemp, ric, veic, loc, rid, typ);

                                                                            if(gg1>0){
                                                                                bestdist=disteffect;
                                                                                bestr1=rr1;
                                                                                bestr2=rr2;
                                                                                p1best=p1;
                                                                                d1best=d1;
                                                                                p2best=p2;
                                                                                d2best=d2;
                                                                                gbest=g;
                                                                                lbest=l;
                                                                                bestv=j;
                                                                                inserted=true;


                                                                            }


                                                                        }

                                                                    }

                                                                }else{g=Sol->route[vei]->length-1;}
                                                            }else{g=Sol->route[vei]->length-1;}

                                                        }

                                                        g++;
                                                    }
                                                }

                                            }



                                        }
                                    }//


                                }else{//esto es de la primera parte! De lo de pickup
                                    Sol->route[j]->calc_earl1(p2, d2, rr1, earl, MatTemp, ric);
                                    feasride1=Sol->route[j]->ridetime_feas_D22(d2, ric, rr1, MatTemp, earl);
                                    if(feasride1>0){
                                        feas41=Sol->route[j]->check_feas_D_tw221(earl, p2, d2, ric);
                                        feas51=Sol->route[j]->check_cap_from23(p2,d2,veic,  rr1, ric);
                                        if(feas41>0 && feas51>0){

                                            feas6=Sol->route[j]->check_feas_FirstRequest(p2, d2, MatTemp, ric, rr1);
                                            if(feas6>0){
                                                gg=Sol->route[j]->EightStep_2(p2,d2, rr1, MatTemp, ric ,veic);

                                                if(gg>0){//it is feasible to insert the first request on the position of the second one.
                                                    //Now I am going to try and insert the second request in all the possible positions in the first one.
                                                    for(int l=1; l<Sol->route[vei]->length-2; l++){

                                                        feas22=Sol->route[vei]->capacity_P_feasibility_3(l, rr2,  veic, ric, p1, d1, rr1);

                                                        feas11=Sol->route[vei]->tw_P_feasibility_3(ric, MatTemp,  l, rr2, p1, d1, loc, typ, rid);
                                                        if(feas11>0 && feas22>0){

                                                            g=l+1;
                                                            while(g<Sol->route[vei]->length-1){
                                                                double disteffect=0;


                                                                Sol->route[vei]->calc_loc3(rr2, ric, p1, d1, l , g, loc);



                                                                if(g==l+1){

                                                                    disteffect=0;
                                                                    disteffect=Sol->route[vei]->effect_of_inserting_req_on_pos_3(rr2, l, g, D, ric, p1, d1, loc);
                                                                    if(disteffect<bestdist){
                                                                        Sol->route[vei]->calc_typ3(rr2, p1, d1, l , g, typ);
                                                                        Sol->route[vei]->calc_rid3(rr2, p1, d1, l , g, rid);
                                                                        //  earl=Sol->route[vei]->calc_earl3(F[k], ric, p1, d1, l , g);

                                                                        feas66=Sol->route[vei]->check_feas_D_tw1_3(rr2,ric,l,g,MatTemp,veic, loc, rid, typ);
                                                                        if(feas66>0){
                                                                            gg1=Sol->route[vei]->EightStep_3(MatTemp, ric, veic, loc, rid, typ);

                                                                            if(gg1>0){
                                                                                bestdist=disteffect;
                                                                                bestr1=rr1;
                                                                                bestr2=rr2;
                                                                                p1best=p1;
                                                                                d1best=d1;
                                                                                p2best=p2;
                                                                                d2best=d2;
                                                                                gbest=g;
                                                                                lbest=l;
                                                                                bestv=j;
                                                                                inserted=true;
                                                                            }
                                                                        }
                                                                    }
                                                                }else{
                                                                    Sol->route[vei]->calc_typ3(rr2, p1, d1, l , g, typ);
                                                                    Sol->route[vei]->calc_rid3(rr2, p1, d1, l , g, rid);
                                                                    Sol->route[vei]->calc_earl3( ric, typ, loc, rid, earl, MatTemp);

                                                                    feasride2=Sol->route[vei]->calc_feas_rtime(l,g,ric,F[k],MatTemp, earl);
                                                                    if(feasride2>0){
                                                                        bool feas44, feas55;
                                                                        disteffect=0;
                                                                        feas44=Sol->route[vei]->check_feas_D_tw22(earl, typ, rid, l, g, ric);

                                                                        feas55=Sol->route[vei]->check_cap_from22(l, g, veic, rr2, ric, rid, typ);

                                                                        if(feas44>0 && feas55>0){
                                                                            disteffect=Sol->route[vei]->effect_of_inserting_req_on_pos_3(rr2, l, g, D, ric, p1, d1, loc);

                                                                            if(disteffect<bestdist){

                                                                                feas66=Sol->route[vei]->check_feas_D_tw1_3(rr2,ric,l,g,MatTemp,veic, loc, rid, typ);
                                                                                if(feas66>0){
                                                                                    gg1=Sol->route[vei]->EightStep_3(MatTemp, ric, veic, loc, rid, typ);

                                                                                    if(gg1>0){
                                                                                        bestdist=disteffect;
                                                                                        bestr1=rr1;
                                                                                        bestr2=rr2;
                                                                                        p1best=p1;
                                                                                        d1best=d1;
                                                                                        p2best=p2;
                                                                                        d2best=d2;
                                                                                        gbest=g;
                                                                                        lbest=l;
                                                                                        bestv=j;
                                                                                        inserted=true;


                                                                                    }


                                                                                }

                                                                            }

                                                                        }else{g=Sol->route[vei]->length-1;}
                                                                    }else{g=Sol->route[vei]->length-1;}

                                                                }

                                                                g++;
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
                       // std::cout << " LAST F " << std::endl;
                        //for(int kk=0; kk<numofreq2; kk++){
                         //   std::cout << F[kk] << " " ;
                        //}std::cout << std::endl;
                    }




                }
            }

        }
        //std::cout << " LAST E" << std::endl;
        //for(int iii=0; iii<numofreq; iii++){
         //   std::cout << E[iii] << " " ;
        //}std::cout << std::endl;
    }


    if(inserted>0){
        best_v=bestv;
       // std::cout << "1" << std::endl;

      //  std::cout << bestr1 << " " << p1best << " " << d1best << std::endl;
        //primero takeout request
       // std::cout << "antes" << std::endl;
        //Sol->route[vei]->display_route();
        best_r1=bestr1;
        best_r2=bestr2;
        or_pickup_pos_1=p1best;
        or_delivery_pos_1=d1best;
        Sol->route[vei]->delete_req_2(bestr1,  MatTemp, p1best, d1best);
        or_pickup_pos_2=p2best;
        or_delivery_pos_2=d2best;
        //std::cout << "despues" << std::endl;
        //Sol->route[vei]->display_route();
        //std::cout << "antes" << std::endl;
        //Sol->route[bestv]->display_route();
        //std::cout << bestr2 << " " << p2best << " " << d2best << std::endl;
        Sol->route[bestv]->delete_req_2(bestr2,  MatTemp, p2best, d2best);
        //std::cout << "despues" << std::endl;
        //Sol->route[bestv]->display_route();
        //then insert request
        //std::cout << bestr2 << " " << lbest << " " << gbest << std::endl;
        //std::cout << "antes" << std::endl;
        //Sol->route[vei]->display_route();

        Sol->route[vei]->insert_request_on(bestr2, lbest, gbest, ric, MatTemp);
        //std::cout << "despues" << std::endl;
        //Sol->route[vei]->display_route();
        //std::cout << bestr1 << " " << p2best << " " << d2best << std::endl;
        //std::cout << "antes" << std::endl;
        //Sol->route[bestv]->display_route();
        Sol->route[bestv]->insert_request_on(bestr1, p2best, d2best, ric, MatTemp);
        //std::cout << "despues" << std::endl;
       // Sol->route[bestv]->display_route();
    //insert whatever in whatever


        Sol->route[bestv]->calculate_capacity(ric,veic);

        Sol->route[bestv]->calculate_earliest_latest(ric, MatTemp, veic);
        Sol->route[bestv]->Eightstepevaluationscheme(MatTemp,ric,veic);
        Sol->route[vei]->calculate_capacity(ric,veic);
        Sol->route[vei]->calculate_earliest_latest(ric, MatTemp, veic);
        Sol->route[vei]->Eightstepevaluationscheme(MatTemp,ric,veic);
        Sol->route[bestv]->totaldist=Sol->route[bestv]->calculatedist(MatTemp);
        Sol->route[vei]->totaldist=Sol->route[vei]->calculatedist(MatTemp);



        Sol->updatecost(numVeicoli);
        Sol->updateusedvehicles(numVeicoli);
       // std::cout<< "adf: " << std::endl;
       // Sol->DisplaySolution();



    }


    //
    delete[] F;
    delete[] E;

    //delete[] loc;
   // delete[] typ;
   // delete[] rid;
    //delete[] earl;




}

emili::Solution* ExchangeNeighborhood::computeStep(emili::Solution* step)
{
    if(vei>=(num_r-1))
    {
        return nullptr;
    }
    else
    {
        SolutionVRP* sol = (SolutionVRP*) step;
        best_v = -1;
        while(sol->route[vei]->numRicRoute==0 && vei<num_r)
            vei++;
        Exchange_Neighborhood(sol, vei, best_v, best_r1, best_r2, or_pickup_pos_1, or_pickup_pos_2, or_delivery_pos_1, or_delivery_pos_2, inst.numVeicoli0, inst.rc, inst.T,inst.vec, inst.numRichieste0, inst.MCV, inst.Dist);

        //std::cout << sol->numAddRoutes << std::endl;
        //std::cout <<"FOPT " << sol->getSolutionValue() << std::endl;
        vei++;
    }
    return step;
}


void ExchangeNeighborhood::reverseLastMove(emili::Solution* step)
{
    if(best_v>-1)
    {
        //reverse move
        SolutionVRP* Sol = (SolutionVRP*) step;

        Sol->route[vei-1]->delete_req(best_r2, inst.T);
        Sol->route[vei-1]->insert_request_on(best_r1, or_pickup_pos_1, or_delivery_pos_1, inst.rc, inst.T);
        Sol->route[best_v]->delete_req(best_r1, inst.T);
        Sol->route[best_v]->insert_request_on(best_r2, or_pickup_pos_2, or_delivery_pos_2, inst.rc, inst.T);
        Sol->route[vei-1]->calculate_capacity(inst.rc,inst.vec);
        Sol->route[best_v]->calculate_capacity(inst.rc,inst.vec);
        Sol->route[vei-1]->calculate_earliest_latest(inst.rc, inst.T, inst.vec);
        Sol->route[vei-1]->Eightstepevaluationscheme(inst.T,inst.rc,inst.vec);
        Sol->route[best_v]->calculate_earliest_latest(inst.rc, inst.T, inst.vec);
        Sol->route[best_v]->Eightstepevaluationscheme(inst.T,inst.rc,inst.vec);
        Sol->route[vei-1]->totaldist=Sol->route[vei-1]->calculatedist(inst.Dist);
        Sol->route[best_v]->totaldist=Sol->route[best_v]->calculatedist(inst.Dist);
        //Solaux->copyrouteinsol(rout);
        //	std::cout << "request" << E[i] << "inserted in vehicle, " << bestv<< "\n";
        //Solaux->route[bestv]->display_route();
        Sol->updatecost(inst.numVeicoli0);
        Sol->updateusedvehicles(inst.numVeicoli0);

    }

}

ExchangeNeighborhood::NeighborhoodIterator ExchangeNeighborhood::begin(emili::Solution* base)
{
    SolutionVRP* sol = (SolutionVRP*) base;
    vei = 0;
    num_r = sol->numRoutes+sol->numAddRoutes;
    best_v = -1;
    best_r1 = -1;
    best_r2 = -1;
    or_pickup_pos_1 = -1;
    or_pickup_pos_2 = -1;
    or_delivery_pos_1 = -1;
    or_delivery_pos_2 = -1;
    return NeighborhoodIterator(this,base);
}

void ExchangeNeighborhood::reset()
{

}

emili::Solution* ExchangeNeighborhood::random(emili::Solution* currentSolution)
{
    SolutionVRP* sol = (SolutionVRP*) currentSolution;
    SolutionVRP* nsol = new SolutionVRP();
    nsol->CopySolution(sol);
    int random_vei = emili::generateRandomNumber()%(sol->numRoutes+sol->numAddRoutes);
    Exchange_Neighborhood(sol, random_vei, best_v, best_r1, best_r2, or_pickup_pos_1, or_pickup_pos_2, or_delivery_pos_1, or_delivery_pos_2, inst.numVeicoli0, inst.rc, inst.T,inst.vec, inst.numRichieste0, inst.MCV, inst.Dist);

    return nsol;
}


//Exchange vehicle for EMILI


void Exchange_Vehicles(SolutionVRP* Sol, int vei, int& vei1, int numVeicoli, int numRichieste, std::vector<std::vector<int>> &MatCompVei, std::vector<std::vector<double>> &Dist, std::vector<std::vector<double>> &Time, std::vector<RichiesteServizio*> &ric, std::vector<Veicoli*> &veic){

    int numreq1=Sol->route[vei]->numRicRoute;
    int numreq2;
    double dist, bestdist=FLT_MAX;

    int best_v;
    int* E=new int[numreq1];
    bool feas, feas1, inserted=false;
    E=Sol->route[vei]->count_request(E);
    int* F=new int[numRichieste];
    for(int i=0; i<numVeicoli+Sol->numAddRoutes; i++){

        if(i!=vei){
            //Are all the requests from vei compatible with i, and are all requests from i compatible with vei
            feas=feasible_with_vehicle(numreq1, E, i, MatCompVei);


            if(feas>0){


                if(Sol->route[i]->numRicRoute>0){
                    numreq2=Sol->route[i]->numRicRoute;
                    F=Sol->route[i]->count_request(F);

                    feas1=feasible_with_vehicle(numreq2, F, vei, MatCompVei);

                    if(feas1>0){

                            //calculate the effect of the distance choose the best and then
                        dist=Sol->calculate_distance_effect(i, vei, Dist);
                        if(bestdist<dist){

                            best_v=i;

                            inserted=true;
                        }
                    }

                }else{

                    //calculate the effect of the distance choose the best and then
                    dist=Sol->calculate_distance_effect(i, vei, Dist);
                    if(bestdist<dist){

                        best_v=i;

                        inserted=true;
                    }



                }



                }
        }

    }


    if(inserted>0){

        Sol->ExchangeRoutes(vei, best_v, veic, Time, ric, Dist);
    //change both routes
        vei1=best_v;

    }


}




emili::Solution* ExchangeVehicleNeighborhood::computeStep(emili::Solution* step)
{
    if(vei>=(num_r-1))
    {
        return nullptr;
    }
    else
    {
        SolutionVRP* sol = (SolutionVRP*) step;
        vei1 = -1;
        while(sol->route[vei]->numRicRoute==0 && vei<num_r)
            vei++;

        Exchange_Vehicles(sol, vei,vei1, inst.numVeicoli0, inst.numRichieste0, inst.MCV, inst.Dist, inst.T, inst.rc, inst.vec);

        //std::cout << sol->numAddRoutes << std::endl;
        //std::cout <<"FOPT " << sol->getSolutionValue() << std::endl;
        vei++;
    }
    return step;
}


void ExchangeVehicleNeighborhood::reverseLastMove(emili::Solution* step)
{
    if(vei1>-1)
    {
        //reverse move
        SolutionVRP* Sol = (SolutionVRP*) step;

        Sol->ExchangeRoutes(vei, vei1, inst.vec, inst.T, inst.rc, inst.Dist);


    }

}

ExchangeVehicleNeighborhood::NeighborhoodIterator ExchangeVehicleNeighborhood::begin(emili::Solution* base)
{
    SolutionVRP* sol = (SolutionVRP*) base;
    vei = 0;
    num_r = sol->numRoutes+sol->numAddRoutes;
    vei1 = -1;

    return NeighborhoodIterator(this,base);
}

void ExchangeVehicleNeighborhood::reset()
{

}

emili::Solution* ExchangeVehicleNeighborhood::random(emili::Solution* currentSolution)
{
    SolutionVRP* sol = (SolutionVRP*) currentSolution;
    SolutionVRP* nsol = new SolutionVRP();
    nsol->CopySolution(sol);
    int random_vei = emili::generateRandomNumber()%(sol->numRoutes+sol->numAddRoutes);
    Exchange_Vehicles(sol, random_vei, vei1, inst.numVeicoli0, inst.numRichieste0, inst.MCV, inst.Dist, inst.T, inst.rc, inst.vec);

    return nsol;
}


