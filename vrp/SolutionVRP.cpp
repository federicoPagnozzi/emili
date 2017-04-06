//
//  Solution.cpp
//  OrCode
//
//  Created by garazi on 25/01/2017.
//  Sol->setSolutionValue(right © 2017 Garazi. All rights reserved.
//

#include "SolutionVRP.hpp"
#include<iostream>

#include "Veicoli.hpp"
#include "Route.hpp"


SolutionVRP::SolutionVRP():emili::Solution(std::numeric_limits<double>::max()) {
    
    numRoutes=0;
    numAddRoutes=0;
    
}

SolutionVRP::~SolutionVRP() {

    for(int i=0; i<route.size(); i++)
    {
        delete route[i];
    }

}


void SolutionVRP::add_empty_route(){
    
    route.push_back(new Route());
    
    
    
}

void SolutionVRP::add_num_empty_route(int n){
    
    //route.reserve(n);
    for(int i=0; i<n; i++){
        add_empty_route();
        
    }
    
}
int SolutionVRP::get_num_routes(){
    
    return route.size();
    
    
}
void SolutionVRP::ClearSol(){
    
    deleteallroutes(numRoutes+numAddRoutes);
    
    
    numRoutes=0;
    numAddRoutes=0;
    solution_value=0;
    
    
    
    
}
void SolutionVRP::Initialize(std::vector<Veicoli*> &veic, int numVeicoli, std::vector<std::vector<double>> &Time){
    
    
    
    if(route.size()>numVeicoli){
        for(int i=numVeicoli;i<route.size();i++){
            for(int a=0;a<route[i]->locations.size(); a++){
                route[i]->remove_pos(a);
            }
            delete route[i];
            route.erase(route.begin()+i);
        }
    }
    
    for(int i=0; i<numVeicoli; i++){
        
        route[i]->InitializeSol(veic,  numVeicoli, Time, i);
        
    }
    
    
    solution_value=0;
    numRoutes=numVeicoli;
    numAddRoutes=0;
    
    
}

void SolutionVRP::copyrouteinsol(Route* r){
    
    route[r->id]->id=r->id;
    route[r->id]->length=r->length;
    route[r->id]->numRicRoute=r->numRicRoute;
    route[r->id]->totaldist=r->totaldist;
    route[r->id]->locations=r->locations;
    route[r->id]->Ricid=r->Ricid;
    route[r->id]->type=r->type;
    route[r->id]->arrival=r->arrival;
    route[r->id]->departure=r->departure;
    route[r->id]->waiting=r->waiting;
    route[r->id]->earliest=r->earliest;
    route[r->id]->latest=r->latest;
    route[r->id]->cap1=r->cap1;
    route[r->id]->maxcap1=r->maxcap1;
    route[r->id]->cap2=r->cap2;
    route[r->id]->maxcap2=r->maxcap2;
    route[r->id]->cap3=r->cap3;
    route[r->id]->maxcap3=r->maxcap3;
    route[r->id]->cap4=r->cap4;
    route[r->id]->maxcap4=r->maxcap4;
    
    
}

Route* SolutionVRP::copyvehicleinroute(Route* r, int j){
    
    r->id=route[j]->id;
    
    r->length=route[j]->length;
    
    r->numRicRoute=route[j]->numRicRoute;
    r->totaldist=route[j]->totaldist;
    r->locations=route[j]->locations;
    r->Ricid=route[j]->Ricid;
    r->type=route[j]->type;
    r->arrival=route[j]->arrival;
    r->departure=route[j]->departure;
    r->waiting=route[j]->waiting;
    r->earliest=route[j]->earliest;
    r->latest=route[j]->latest;
    r->cap1=route[j]->cap1;
    r->maxcap1=route[j]->maxcap1;
    r->cap2=route[j]->cap2;
    r->maxcap2=route[j]->maxcap2;
    r->cap3=route[j]->cap3;
    r->maxcap3=route[j]->maxcap3;
    r->cap4=route[j]->cap4;
    r->maxcap4=route[j]->maxcap4;
    
    
    return r;
}

void SolutionVRP::updateusedvehicles(int numVeicoli){
    int cont=0;
    for(int i=0;i<numVeicoli; i++){
        if(route[i]->length>2){
            cont+=1;
        }
    }
    numRoutes=cont;
}
void SolutionVRP::updatecost(int numVeicoli){
    double c=0;
    for(int i=0; i<numVeicoli+numAddRoutes; i++){
        if(route[i]->length>2){
            
            c+=route[i]->totaldist;
            
        }
        
    }
    solution_value=c;
}

void SolutionVRP::addAdditionalRoute(int numVeicoli, std::vector<std::vector<double>>  &Time, std::vector<Veicoli*> &veic){
    int a=numVeicoli+numAddRoutes;
   // std::cout << a << std::endl;
    add_empty_route();
    
    route[a]->Initializeaddroute(a, Time, veic, numVeicoli);
    
    numAddRoutes+=1;
    
    
    
}

double SolutionVRP::calculate_dist(std::vector<std::vector<double>> & D)
{
    double normale=0;
    for (int i = 0; i < numRoutes+numAddRoutes; ++i) {
        normale += route[i]->calculatedist(D);
    }
    return normale;
}

void SolutionVRP::copyexistingsol(SolutionVRP* Sol){
    
    numRoutes=Sol->numRoutes;
    numAddRoutes=Sol->numAddRoutes;
    solution_value=Sol->solution_value;
    for(int i=0;i<numRoutes+numAddRoutes;i++){
        
        route[i]->CopyRoute(Sol->route[i]);
    }
    
}

void SolutionVRP::CopySolution(SolutionVRP* Sol){
    
    
    
    deleteallroutes(numRoutes+numAddRoutes);
    
    
    //route.clear();
    numRoutes=Sol->numRoutes;
    numAddRoutes=Sol->numAddRoutes;
    solution_value=Sol->solution_value;
    
    
    add_num_empty_route(Sol->numRoutes+Sol->numAddRoutes);
    for(int i=0;i<numRoutes+numAddRoutes;i++){
        
        route[i]->CopyRoute(Sol->route[i]);
    }
    
}


void SolutionVRP::DisplaySolution(){
    std::cout<<"cost: " << solution_value << std::endl;
    for(int i=0;i<numRoutes+numAddRoutes;i++){
        std::cout<<"route " << i << std::endl;
        route[i]->display_route();
        
    }
    
    
    
}

void SolutionVRP::deleteallroutes(int n){
 //   int l;
    for(int i=n-1;i>-1;i--){
//        l=route[i]->locations.size();
//        for(int j=l-1;j>-1;j--){
//            route[i]->remove_pos(j);
//        }
        route[i]->locations.clear();
        route[i]->locations.clear();
        route[i]->Ricid.clear();
        route[i]->type.clear();
        route[i]->arrival.clear();
        route[i]->departure.clear();
        route[i]->waiting.clear();
        route[i]->earliest.clear();
        route[i]->latest.clear();
        route[i]->cap1.clear();
        route[i]->maxcap1.clear();
        route[i]->cap2.clear();
        route[i]->maxcap2.clear();
        route[i]->cap3.clear();
        route[i]->maxcap3.clear();
        route[i]->cap4.clear();
        route[i]->maxcap4.clear();
        delete route[i];
        route.erase(route.begin()+i);
//        
//        
//        
    }
    //route.clear();
}


void SolutionVRP::delete_route(int a, int numVeicoli){
    int size_route =route[a]->length;
    if(a>=numVeicoli){
        
        
        
        for (int i =size_route-1; i > -1; i--) {
            route[a]->remove_pos(i);
        }
        
        delete route[a];
        
        route.erase(route.begin() + a);
        numAddRoutes-=1;
        
        for(int i=a; i<numAddRoutes+numRoutes; i++){
            
            route[a]->rewriteid(a);
            
        }
        
    }else{
        for (int i =size_route-2; i > 0; i--) {
            route[a]->remove_pos(i);
        }
        route[a]->length=2;
        route[a]->numRicRoute=0;
        route[a]->totaldist=0;
        
    }
    updatecost(numVeicoli);
    updateusedvehicles(numVeicoli);
    
}

bool SolutionVRP::accepted(int T, SolutionVRP* Sol_Curr, int numVeicoli){
    
    bool accep;
    
    if((numRoutes+numAddRoutes>numVeicoli && numAddRoutes<Sol_Curr->numAddRoutes) || solution_value<Sol_Curr->solution_value+T){
        accep=true;
    }else{
        accep=false;
        
    }
    
    return accep;
}

bool SolutionVRP::bestaccepted(SolutionVRP* Sol_Curr, int numVeicoli){
    bool accep;
    
    
    return accep;
}

void SolutionVRP::delete_routes(){
    
    
    for(int i=0;i<numRoutes+numAddRoutes;i++){
        
        
        delete route[i];
        
        route.clear();
        
        
    }
    
}

double SolutionVRP::calculate_effect_2(int i, int j, int k, int l, std::vector<std::vector<double>> & D){
    double d1=0, d2=0, d=0;
    //calcular toda...
    for(int ii=0; ii<j;ii++){
        d1+=D[route[i]->locations[ii]][route[i]->locations[ii+1]];
    }
    for(int ii=0; ii<l; ii++){
        d2+=D[route[k]->locations[ii]][route[k]->locations[ii+1]];
    }
    for(int ii=l+1;ii<route[k]->length-1;ii++){
        d1+=D[route[k]->locations[ii]][route[k]->locations[ii+1]];
    }
    for(int ii=j+1;ii<route[i]->length-1;ii++){
        d2+=D[route[i]->locations[ii]][route[i]->locations[ii+1]];
    }
    d1+=D[route[i]->locations[j]][route[k]->locations[l+1]];
    
    d2=D[route[k]->locations[l]][route[i]->locations[j+1]];
    
    d=(d1-route[i]->totaldist)+(d2-route[k]->totaldist);
    
    return d;
}
double SolutionVRP::calculate_effect(int i, int j, int k, int l, std::vector<std::vector<double>> & D){
    double d1=0, d2=0, d=0;
    //no tiene sentido.........
    
    d1=D[route[i]->locations[j]][route[i]->locations[j+1]]-D[route[i]->locations[j]][route[k]->locations[l+1]];
    d2=D[route[k]->locations[l]][route[k]->locations[l+1]]-D[route[k]->locations[l]][route[i]->locations[j+1]];
    d=d1+d2;
    
    return d;
}

bool SolutionVRP::cap_feas(int i, int j,int k,int l,std::vector<Veicoli*> &veic){
    
    int veh1, veh2;
    veh1=route[i]->veh;
    veh2=route[k]->veh;
    bool f1, f2, f;
    if(route[veh1]->maxcap1[j+1]<=veic[veh2]->staff && route[veh1]->maxcap2[j+1]<=veic[veh2]->seated && route[veh1]->maxcap3[j+1]<=veic[veh2]->stretcher && route[veh1]->maxcap4[j+1]<=veic[veh2]->wheelchair){
        f1=true;}else{
            f1=false;
        }
    if(route[veh2]->maxcap1[l+1]<=veic[veh1]->staff && route[veh2]->maxcap2[l+1]<=veic[veh1]->seated && route[veh2]->maxcap3[l+1]<=veic[veh1]->stretcher && route[veh2]->maxcap4[l+1]<=veic[veh1]->wheelchair){
        f2=true;}else{
            f2=false;
        }
    
    if(f1==true && f2==true){
        f=true;
    }else{f=false;}
    return f;
}


bool SolutionVRP::tw_feas(int i,int j,int k,int l, std::vector<RichiesteServizio*> &ric,std::vector<std::vector<double>> &D){
    
    bool f, f1, f2;
    
    if(route[i]->earliest[j]+ric[route[i]->Ricid[j]]->st+D[route[i]->locations[j]][route[k]->locations[l+1]]<=route[k]->latest[l+1]+D[route[k]->locations[(route[k]->length)-2]][route[i]->locations[0]]){
        
        f1=true;
        
    }else{
        f1=false;
    }
    
    if(route[k]->earliest[l]+ric[route[k]->Ricid[l]]->st+D[route[k]->locations[l]][route[i]->locations[j+1]]<=route[i]->latest[j+1]+D[route[i]->locations[(route[i]->length)-2]][route[k]->locations[0]]){
        
        f2=true;
        
    }else{
        f2=false;
    }
    
    if(f1==true && f2==true){
        f=true;
    }else{f=false;}
    return f;
}
 //CAMBIA ESTO !!!! PORQUE SI ES MAS LARGA UNA QUE LA OTRA HAY Q CAMBIAR ALGO ! PIENSA
 
void SolutionVRP::createnewroute(int i, int j, int k, int l, std::vector<RichiesteServizio*> &ric, std::vector<std::vector<double>> &D, std::vector<Veicoli*> &veic){
    
  
    int loc_aux[route[i]->length-(j+2)];
    int rd_aux[route[i]->length-(j+2)];
    int ty_aux[route[i]->length-(j+2)];
    int orlength=route[i]->length;
    int o, p;

    for(o=0;o<orlength-(j+2); o++)
    {
        loc_aux[o]=route[i]->locations[o+j+1];
        rd_aux[o]=route[i]->Ricid[o+j+1];
        ty_aux[o]=route[i]->type[o+j+1];
    }
   // for(o=0;o<orlength-(j+2); o++)
   // {
   //     std::cout << loc_aux[o]  << " " ;
   //     std::cout << rd_aux[o]  << " " ;
   //     std::cout << ty_aux[o]  << " " ;
   //     std::cout << std::endl;
   // }
    int nl=(j+1)+(route[k]->length-(l+1));
    //std::cout << nl << std::endl;
    route[i]->resize_from(nl);

    
    for(o=0; o<route[k]->length-(l+2); o++){
        route[i]->locations[o+j+1]=route[k]->locations[o+l+1];
        route[i]->Ricid[o+j+1]=route[k]->Ricid[o+l+1];
        route[i]->type[o+j+1]=route[k]->type[o+l+1];
    }
    
    //route[i]->display_route();

    int nl1=(l+1)+(orlength-(j+1));

    //std::cout << nl1 << std::endl;
    route[k]->resize_from(nl1);
    //std::cout << nl << std::endl;
    for(o=0; o<orlength-(j+2); o++){
        route[k]->locations[o+l+1]=loc_aux[o];
        route[k]->Ricid[o+l+1]=rd_aux[o];
        route[k]->type[o+l+1]=ty_aux[o];
    }

    
    //CALCULATE 8 STEP EVALUATION SCHEME BUT WITHOUT THE BOOL, Try to change it, try to make the otherone faster... ????????????
    
    
    
}

void SolutionVRP::changetworoutes(int i, int j, int k, int l, std::vector<RichiesteServizio*> & ric, std::vector<std::vector<double>> &D ,std::vector<Veicoli*> & veic){
    
    //createnewroute(i,j,k, l, ric, D, route[i], veic, MaxTimeRoute);

    createnewroute(i,j,k, l, ric, D, veic);

    
    route[i]->numRicRoute=(route[i]->length-2)/2;
    route[k]->numRicRoute=(route[k]->length-2)/2;
    
    route[i]->totaldist=route[i]->calculatedist(D);
    route[k]->totaldist=route[k]->calculatedist(D);
    
    
    //calculate from j and from l
    route[i]->calculate_capacity(ric, veic);
    route[k]->calculate_capacity(ric,veic);
    
    route[i]->calculate_earliest_latest(ric, D, veic);
    route[k]->calculate_earliest_latest(ric, D, veic);
    route[i]->Eightstepevaluationscheme(D, ric, veic);
    route[k]->Eightstepevaluationscheme(D, ric, veic);
    
}


bool SolutionVRP::Eightstepevaluationscheme_newr(std::vector<std::vector<double>> & Time, std::vector<RichiesteServizio*> & ric, int a, int b, int c, int d, int numRichieste, std::vector<Veicoli*> & veic){
    int i, req;
    bool feas=true;
    int len;
    int p;
    
   // int len=length+2;
    double* arr, *dep, * wait;
    arr=new double[2*numRichieste+2];
    dep=new double[2*numRichieste+2];
    wait=new double[2*numRichieste+2];
    int* loc, *rid, *typ;
    loc=new int[2*numRichieste+2];
    rid=new int[2*numRichieste+2];
    typ=new int[2*numRichieste+2];
    
    for(i=0;i<b+1;i++){
        loc[i]=route[a]->locations[i];
        rid[i]=route[a]->Ricid[i];
        typ[i]=route[a]->type[i];
    }
    i=d+1;
    p=b+1;
    while(i<route[c]->length){
        loc[p]=route[c]->locations[i];
        rid[p]=route[c]->Ricid[i];
        typ[p]=route[c]->type[i];
        i+=1;
        p+=1;
    }
    
    len=p;
    
    
    
    double F, v1=0, v2=0, v3=0;
    
    
    for(i=0; i<len; i++){
        arr[i]=0;
        dep[i]=0;
        wait[i]=0;
    }
    for(i=1; i<len-1; i++){
        
        arr[i]=dep[i-1]+Time[loc[i-1]][loc[i]];
        
        if(typ[i]==1){
            if(arr[i]-ric[rid[i]]->timewinPmin>=0 && arr[i]-ric[rid[i]]->timewinPmax<=0){
                dep[i]=arr[i]+ric[rid[i]]->st;
                wait[i]=0;
                
            }else{if(arr[i]-ric[rid[i]]->timewinPmin<0){
                
                dep[i]=ric[rid[i]]->timewinPmin+ric[rid[i]]->st;
                wait[i]=dep[i]-ric[rid[i]]->st-arr[i];
            }else{
                dep[i]=arr[i]+ric[rid[i]]->st;
                
                v1+=arr[i]-ric[rid[i]]->timewinPmax;
                //Viol+=arr[j]-B[int(Sol1[i][j][1])][4];
            }
            }
        }else{//if it is a delivery node
            if(arr[i]-ric[rid[i]]->timewinDmin>=0 && arr[i]-ric[rid[i]]->timewinDmax<=0){
                dep[i]=arr[i]+ric[rid[i]]->st;
                wait[i]=0;
            }else{ if(arr[i]-ric[rid[i]]->timewinDmin<0){
                dep[i]=ric[rid[i]]->timewinDmin+ric[rid[i]]->st;
                wait[i]=dep[i]-ric[rid[i]]->st-arr[i];
            }else{
                dep[i]=arr[i]+ric[rid[i]]->st;
                
                v1+=arr[i]-ric[rid[i]]->timewinDmax;
                //Viol+=arr[j]-B[int(Sol1[i][j][1])][6];
            }}
        }
        
    }
    arr[len-1]=dep[len-2]+Time[loc[len-1]][loc[len-2]];
    dep[len-1]=0;
    wait[len-1]=0;
    
    //for(i=0;i<length;i++){
    //	std::cout << arr[i] << " " << dep[i] << " " << wait[i] << "\n" ;
    //}
    
    
    if(v1>0){
        //go to step 8
        //std:: cout << "false \n";
    }else{
        //std:: cout << "true \n";
        //calculate the forward slack
        //F=calculateF(1, ric, MaxRideTime);
        F=calcF(rid, wait,  arr,  dep, 1, ric, typ,  len);
        double G;
        G=0;
        
        for(i=1; i<len; i++){
            G+=wait[i];
        }
        if(F<=G){
            G=F;
        }
        //std:: cout << "G is" << G << "\n" ;
        arr[1]=(dep[1]-ric[rid[1]]->st)+G;
        dep[0]=arr[1]-Time[loc[0]][loc[1]];
        v1=0;
        for(i=1; i<len-1; i++){//hasta menos uno porque not iene niungun sentido que se calcule lo del deposito
            arr[i]=dep[i-1]+Time[loc[i]][loc[i-1]];
            if(typ[i]==1){ //if it is a pickup node
                if(arr[i]-ric[rid[i]]->timewinPmin>=0 && arr[i]-ric[rid[i]]->timewinPmax<=0){
                    dep[i]=arr[i]+ric[rid[i]]->st;
                    wait[i]=0;
                }else{if(arr[i]-ric[rid[i]]->timewinPmin<0){
                    
                    dep[i]=ric[rid[i]]->timewinPmin+ric[rid[i]]->st;
                    wait[i]=dep[i]-ric[rid[i]]->st-arr[i];
                }else{
                    dep[i]=arr[i]+ric[rid[i]]->st;
                    
                    v1+=arr[i]-ric[rid[i]]->timewinPmax;
                    //Viol+=arr[j]-B[int(Sol1[i][j][1])][4];
                }}
            }else{//if it is a delivery node
                if(arr[i]-ric[rid[i]]->timewinDmin>=0 && arr[i]-ric[rid[i]]->timewinDmax<=0){
                    dep[i]=arr[i]+ric[rid[i]]->st;
                    wait[i]=0;
                }else{ if(arr[i]-ric[rid[i]]->timewinDmin<0){
                    dep[i]=ric[rid[i]]->timewinDmin+ric[rid[i]]->st;
                    wait[i]=dep[i]-ric[rid[i]]->st-arr[i];
                }else{
                    dep[i]=arr[i]+ric[rid[i]]->st;
                    
                    v1+=arr[i]-ric[rid[i]]->timewinDmax;
                    //Viol+=arr[j]-B[int(Sol1[i][j][1])][6];
                }}
            }
        }
        arr[len-1]=dep[len-2]+Time[loc[len-1]][loc[len-2]];
        dep[len-1]=0;
        wait[len-1]=0;
        
        //for(i=0;i<length;i++){
        //	std::cout << arr[i] << " " << dep[i] << " " << wait[i] << "\n" ;
        //}
        
        
        
        v3=0;
        for(int ii=0; ii<len; ii++){
            
            if(typ[ii]==1){
                
                req=rid[ii];
                for(int jj=ii+1;jj<len;jj++){
                    
                    if(typ[jj]==-1){ //es nodo de delivery
                        if(rid[jj]==req){
                            
                            if(dep[jj]-ric[req]->st-dep[ii]-ric[req]->RideTime>0){
                                v3+=(dep[jj]-ric[req]->st-dep[ii]-ric[req]->RideTime);
                                
                            }else{
                                
                            }
                        }
                    }
                    
                }
            }
            
        }
        if(v3==0)
        {}else{
            v1=0;
            //Viol=0;
            for(i=1; i<len-1; i++){
                if(typ[i]==1){
                    // F=calculateF(i, ric, MaxRideTime);
                    F=calcF(rid, wait,  arr,  dep, i, ric, typ,  len);
                    //set waits
                    G=0;
                    for(int p=i+1; p<len; p++){
                        G+=wait[p];
                    }
                    if(F<G){
                        G=F;
                    }else{
                        
                    }
                    //calculate the others
                    dep[i]=dep[i]+G;
                    wait[i]=dep[i]-ric[rid[i]]->st-arr[i];
                    
                    
                    for(int kk=i+1;kk<len-1;kk++){
                        arr[kk]=dep[kk-1]+Time[loc[kk-1]][loc[kk]];
                        if(typ[kk]==1){ //if it is a pickup node
                            if(arr[kk]-ric[rid[kk]]->timewinPmin>=0 && arr[kk]-ric[rid[kk]]->timewinPmax<=0){
                                dep[kk]=arr[kk]+ric[rid[kk]]->st;
                                wait[kk]=0;
                            }else{ if(arr[kk]-ric[rid[kk]]->timewinPmin<0){
                                
                                dep[kk]=ric[rid[kk]]->timewinPmin+ric[rid[kk]]->st;
                                wait[kk]=dep[kk]-ric[rid[kk]]->st-arr[kk];
                            }else{
                                
                                dep[kk]=arr[kk]+ric[rid[kk]]->st;
                                v1+=arr[kk]-ric[rid[kk]]->timewinPmax;
                                //Viol+=arr[kk]-B[int(Sol1[i][kk][1])][4];
                            }}
                        }else{ //if it is a delivery node
                            //delivery home
                            
                            //if it is a delivery node
                            if(arr[kk]-ric[rid[kk]]->timewinDmin>=-0 && arr[kk]-ric[rid[kk]]->timewinDmax<=0  ){
                                dep[kk]=arr[kk]+ric[rid[kk]]->st;
                                wait[kk]=0;
                            }else{ if(arr[kk]-ric[rid[kk]]->timewinDmin<-0){
                                dep[kk]=ric[rid[kk]]->timewinDmin+ric[rid[kk]]->st;
                                wait[kk]=dep[kk]-ric[rid[kk]]->st-arr[kk];
                            }else{
                                
                                dep[kk]=arr[kk]+ric[rid[kk]]->st;
                                v1+=arr[kk]-ric[rid[kk]]->timewinDmax;
                                //Viol+=arr[kk]-B[int(Sol1[i][kk][1])][6];
                            }}
                            
                            
                        }
                        
                        
                    }
                    
                    arr[len-1]=dep[len-2]+Time[loc[len-2]][loc[len-1]];
                    
                }else{
                    
                }
                
                
                
                
                
            }
            
            
            arr[1]=dep[1]-ric[rid[1]]->st;
            dep[0]=arr[1]-Time[loc[0]][loc[1]];
            
            
            
            
        }
        
        
    }
    
    v3=0;
    for(int ii=0; ii<len; ii++){
        
        if(typ[ii]==1){
            
            req=rid[ii];
            for(int jj=ii+1;jj<len;jj++){
                
                if(typ[jj]==-1){ //es nodo de delivery
                    if(rid[jj]==req){
                        if(dep[jj]-ric[req]->st-dep[ii]-ric[req]->RideTime>0){
                            
                            
                            v3+=(dep[jj]-ric[req]->st-dep[ii]-ric[req]->RideTime);
                        }else{
                        }
                    }
                }
                
            }
        }
        
    }
    
    
    if(arr[len-1]-dep[0]-veic[route[a]->veh]->MaxRoute>0){
        
        //V[i][1]=arr[length-1]-dep[0]-MaxTimeRoute;
        v2=arr[len-1]-dep[0]-veic[route[a]->veh]->MaxRoute;
        
        
    }else{v2=0;}
    
    
    if(v1>0 || v2>0 || v3>0){
        feas=false;
    }
    
    delete[] loc;
    delete[] rid;
    delete[] typ;
    delete[] arr;
    delete[] dep;
    delete[] wait;
    
    return feas;
}


emili::Solution* SolutionVRP::clone()
{
    SolutionVRP* sol = new SolutionVRP();
    sol->CopySolution(this);
    return sol;
}


const void* SolutionVRP::getRawData() const
{
    return (void*) this;
}

void SolutionVRP::setRawData(const void* data)
{
    SolutionVRP* sol = (SolutionVRP*) data;
    this->CopySolution(sol);
}

bool SolutionVRP::isFeasible(){

    int feasible;

    if(numAddRoutes>0){
        feasible=false;
    }
    else{
        feasible=true;

    }


}
