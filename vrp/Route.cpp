//
//  Route.cpp
//  OrCode
//
//  Created by garazi on 25/01/2017.
//  Copyright © 2017 Garazi. All rights reserved.
//



#include"Veicoli.hpp"
#include"Route.hpp"
#include"RichiesteServizio.hpp"
#include<cfloat>
#include<iostream>


Route::Route() {
    id=0;
    length=0;
    numRicRoute=0;
    totaldist=0;
}


void Route::Initializeaddroute(int a, std::vector<std::vector<double>> &Time,  std::vector<Veicoli*> &veic , int numVeicoli){
    
    id=a;
    length=2;
    veh=numVeicoli;
    numRicRoute=0;
    
    
    locations.clear();
    Ricid.clear();
    type.clear();
    arrival.clear();
    departure.clear();
    waiting.clear();
    earliest.clear();
    latest.clear();
    cap1.clear();
    maxcap1.clear();
    cap2.clear();
    maxcap2.clear();
    cap3.clear();
    maxcap3.clear();
    cap4.clear();
    maxcap4.clear();
    locations.insert(locations.begin(), veic[veh]->location);
    
    locations.insert(locations.begin()+1, veic[veh]->location);
    Ricid.insert(Ricid.begin(), -2);
    
    Ricid.insert(Ricid.begin()+1, -2);
    type.insert(type.begin(), 0);
    
    type.insert(type.begin()+1, 0);
    arrival.insert(arrival.begin(), 0);
    
    arrival.insert(arrival.begin()+1, 0);
    
    departure.insert(departure.begin(), 0);
    
    departure.insert(departure.begin()+1, -0);
    waiting.insert(waiting.begin(), 0);
    	
    waiting.insert(waiting.begin()+1, 0);
    
    earliest.insert(earliest.begin(),0);
    
    
    
    earliest.insert(earliest.begin()+1,Time[locations[0]][locations[1]]);
    
    
    latest.insert(latest.begin(),veic[veh]->MaxRoute-Time[locations[0]][locations[1]]);
    
    latest.insert(latest.begin(),veic[veh]->MaxRoute);
    cap1.insert(cap1.begin(), 0);
    cap1.insert(cap1.begin()+1, 0);
    cap2.insert(cap2.begin(), 0);
    cap2.insert(cap2.begin()+1, 0);
    cap3.insert(cap3.begin(), 0);
    cap3.insert(cap3.begin()+1, 0);
    cap4.insert(cap4.begin(), 0);
    cap4.insert(cap4.begin()+1, 0);
    maxcap1.insert(maxcap1.begin(), veic[veh]->staff);
    maxcap1.insert(maxcap1.begin()+1, 0);
    maxcap2.insert(maxcap2.begin(), veic[veh]->seated);
    maxcap2.insert(maxcap2.begin()+1, 0);
    maxcap3.insert(maxcap3.begin(), veic[veh]->stretcher);
    maxcap3.insert(maxcap3.begin()+1, 0);
    maxcap4.insert(maxcap4.begin(), veic[veh]->wheelchair);
    maxcap4.insert(maxcap4.begin()+1, 0);
    
}
Route::~Route() {
//    
//    
//    
//    
//    locations.clear();
//    Ricid.clear();
//    type.clear();
//    arrival.clear();
//    departure.clear();
//    waiting.clear();
//    earliest.clear();
//    latest.clear();
//    cap1.clear();
//    maxcap1.clear();
//    cap2.clear();
//    maxcap2.clear();
//    cap3.clear();
//    maxcap3.clear();
//    cap4.clear();
//    maxcap4.clear();
    
}

void Route::InitializeSol(std::vector<Veicoli*> &veic, int numVeicoli, std::vector<std::vector<double>> &Time, int i){
    
    locations.clear();
    Ricid.clear();
    type.clear();
    arrival.clear();
    departure.clear();
    waiting.clear();
    earliest.clear();
    latest.clear();
    cap1.clear();
    maxcap1.clear();
    cap2.clear();
    maxcap2.clear();
    cap3.clear();
    maxcap3.clear();
    cap4.clear();
    maxcap4.clear();
    
    id=i;
    length=2;
    
    veh=veic[i]->id;
    
    numRicRoute=0;
    locations.insert(locations.begin(), veic[i]->location);
    
    locations.insert(locations.begin()+1, veic[i]->location);
    Ricid.insert(Ricid.begin(), -2);
    
    Ricid.insert(Ricid.begin()+1, -2);
    type.insert(type.begin(), 0);
    
    type.insert(type.begin()+1, 0);
    arrival.insert(arrival.begin(), 0);
    
    arrival.insert(arrival.begin()+1, 0);
    
    departure.insert(departure.begin(), 0);
    
    departure.insert(departure.begin()+1, -0);
    waiting.insert(waiting.begin(), 0);
    
    waiting.insert(waiting.begin()+1, 0);
    
    earliest.insert(earliest.begin(),0);
    
    
    
    earliest.insert(earliest.begin()+1,Time[veic[i]->location][veic[i]->location]);
    
    
    latest.insert(latest.begin(),veic[i]->MaxRoute-Time[veic[i]->location][veic[i]->location]);
    
    latest.insert(latest.begin(),veic[i]->MaxRoute);
    cap1.insert(cap1.begin(), 0);
    cap1.insert(cap1.begin()+1, 0);
    cap2.insert(cap2.begin(), 0);
    cap2.insert(cap2.begin()+1, 0);
    cap3.insert(cap3.begin(), 0);
    cap3.insert(cap3.begin()+1, 0);
    cap4.insert(cap4.begin(), 0);
    cap4.insert(cap4.begin()+1, 0);
    maxcap1.insert(maxcap1.begin(), 0);
    maxcap1.insert(maxcap1.begin()+1, 0);
    maxcap2.insert(maxcap2.begin(), 0);
    maxcap2.insert(maxcap2.begin()+1, 0);
    maxcap3.insert(maxcap3.begin(), 0);
    maxcap3.insert(maxcap3.begin()+1, 0);
    maxcap4.insert(maxcap4.begin(), 0);
    maxcap4.insert(maxcap4.begin()+1, 0);
    
    
    
    
    
}

void Route::clear_route(){
    
    
    
    locations.clear();
    Ricid.clear();
    type.clear();
    arrival.clear();
    departure.clear();
    waiting.clear();
    earliest.clear();
    latest.clear();
    cap1.clear();
    maxcap1.clear();
    cap2.clear();
    maxcap2.clear();
    cap3.clear();
    maxcap3.clear();
    cap4.clear();
    maxcap4.clear();
    
    
}
void Route::insert_pickup(int req, int l, RichiesteServizio* ric){
    locations.insert(locations.begin()+(l+1),ric[req].locationP);
    Ricid.insert(Ricid.begin()+(l+1), req);
    type.insert(type.begin()+(l+1), 1);
    arrival.insert(arrival.begin()+(l+1),0);
    departure.insert(departure.begin()+(l+1), 0);
    waiting.insert(waiting.begin()+(l+1), 0);
    earliest.insert(earliest.begin()+(l+1),0);
    latest.insert(latest.begin()+(l+1), 0);
    cap1.insert(cap1.begin()+(l+1), 0);
    maxcap1.insert(maxcap1.begin()+(l+1), 0);
    cap2.insert(cap2.begin()+(l+1), 0);
    maxcap2.insert(maxcap2.begin()+(l+1), 0);
    cap3.insert(cap3.begin()+(l+1), 0);
    maxcap3.insert(maxcap3.begin()+(l+1), 0);
    cap4.insert(cap4.begin()+(l+1), 0);
    maxcap4.insert(maxcap4.begin()+(l+1), 0);
    length=locations.size();
}
void Route::insert_delivery(int req, int g, RichiesteServizio* ric){
    
    locations.insert(locations.begin()+(g+2),ric[req].locationD);
    
    Ricid.insert(Ricid.begin()+(g+2), req);
    type.insert(type.begin()+(g+2), -1);
    
    arrival.insert(arrival.begin()+(g+2),0);
    departure.insert(departure.begin()+(g+2), 0);
    
    waiting.insert(waiting.begin()+(g+2), 0);
    earliest.insert(earliest.begin()+(g+2),0);
    latest.insert(latest.begin()+(g+2), 0);
    cap1.insert(cap1.begin()+(g+2), 0);
    maxcap1.insert(maxcap1.begin()+(g+2), 0);
    cap2.insert(cap2.begin()+(g+2), 0);
    maxcap2.insert(maxcap2.begin()+(g+2), 0);
    cap3.insert(cap3.begin()+(g+2), 0);
    maxcap3.insert(maxcap3.begin()+(g+2), 0);
    cap4.insert(cap4.begin()+(g+2), 0);
    maxcap4.insert(maxcap4.begin()+(g+2), 0);
    
    length=locations.size();
    numRicRoute=(length-2)/2;
    
    
}

void Route::delete_delivery(int g){
    
    locations.erase(locations.begin()+g+2);
    Ricid.erase(Ricid.begin()+g+2);
    type.erase(type.begin()+g+2);
    
    //numRicRoute=(length-2)/2;
    arrival.erase(arrival.begin()+g+2);
    departure.erase(departure.begin()+g+2);
    waiting.erase(waiting.begin()+g+2);
    earliest.erase(earliest.begin()+g+2);
    latest.erase(latest.begin()+g+2);
    cap1.erase(cap1.begin()+g+2);
    maxcap1.erase(maxcap1.begin()+g+2);
    cap2.erase(cap2.begin()+g+2);
    maxcap2.erase(maxcap2.begin()+g+2);
    cap3.erase(cap3.begin()+g+2);
    maxcap3.erase(maxcap3.begin()+g+2);
    cap4.erase(cap4.begin()+g+2);
    maxcap4.erase(maxcap4.begin()+g+2);
    length=locations.size();
    
}
void Route::display_route(){
    std::cout << "loc \t Ric \t type \t ar \t de \t wa \t ear \t lat \t c1 \t mc1 \t c2 \t mc2 \t c3 \t mc3 \t c4 \t mc4 \n";
    for(int i=0; i<length; i++){
        
        std::cout << locations[i] << "\t" << Ricid[i] << "\t" << type[i]  << "\t" << arrival[i] << "\t" << departure[i] << "\t" << waiting[i] << "\t" << earliest[i] << "\t";
        std::cout << latest[i] << "\t" << cap1[i] << "\t" << maxcap1[i] << "\t" << cap2[i] << "\t" << maxcap2[i] << "\t" << cap3[i] << "\t" << maxcap3[i] << "\t" << cap4[i] << "\t" << maxcap4[i] << "\n";
        
    }
    
    
    
}
/*
void Route::delete_pickup(int l, double** Time){
    
    
    
    locations.erase(locations.begin()+l+1);
    Ricid.erase(Ricid.begin()+l+1);
    type.erase(type.begin()+l+1);
    length=locations.size();
    numRicRoute=(length-2)/2;
    arrival.erase(arrival.begin()+l+1);
    departure.erase(departure.begin()+l+1);
    waiting.erase(waiting.begin()+l+1);
    earliest.erase(earliest.begin()+l+1);
    latest.erase(latest.begin()+l+1);
    cap1.erase(cap1.begin()+l+1);
    maxcap1.erase(maxcap1.begin()+l+1);
    cap2.erase(cap2.begin()+l+1);
    maxcap2.erase(maxcap2.begin()+l+1);
    cap3.erase(cap3.begin()+l+1);
    maxcap3.erase(maxcap3.begin()+l+1);
    cap4.erase(cap4.begin()+l+1);
    maxcap4.erase(maxcap4.begin()+l+1);
    totaldist=calculatedist(Time);
    
}

*/
/*
void Route::insert_req(int req, int l, int g, RichiesteServizio* ric, double** Time){
    locations.insert(locations.begin()+(l),ric[req].locationP);
    locations.insert(locations.begin()+(g),ric[req].locationD);
    Ricid.insert(Ricid.begin()+(l), req);
    Ricid.insert(Ricid.begin()+(g), req);
    type.insert(type.begin()+(l), 1);
    type.insert(type.begin()+(g), -1);
    arrival.insert(arrival.begin()+(l),0);
    departure.insert(departure.begin()+(l), 0);
    waiting.insert(waiting.begin()+(l), 0);
    earliest.insert(earliest.begin()+(l),0);
    latest.insert(latest.begin()+(l), 0);
    cap1.insert(cap1.begin()+(l), 0);
    maxcap1.insert(maxcap1.begin()+(l), 0);
    cap2.insert(cap2.begin()+(l), 0);
    maxcap2.insert(maxcap2.begin()+(l), 0);
    cap3.insert(cap3.begin()+(l), 0);
    maxcap3.insert(maxcap3.begin()+(l), 0);
    cap4.insert(cap4.begin()+(l), 0);
    maxcap4.insert(maxcap4.begin()+(l), 0);
    arrival.insert(arrival.begin()+(g),0);
    departure.insert(departure.begin()+(g), 0);
    
    waiting.insert(waiting.begin()+(g), 0);
    earliest.insert(earliest.begin()+(g),0);
    latest.insert(latest.begin()+(g), 0);
    cap1.insert(cap1.begin()+(g), 0);
    maxcap1.insert(maxcap1.begin()+(g), 0);
    cap2.insert(cap2.begin()+(g), 0);
    maxcap2.insert(maxcap2.begin()+(g), 0);
    cap3.insert(cap3.begin()+(g), 0);
    maxcap3.insert(maxcap3.begin()+(g), 0);
    cap4.insert(cap4.begin()+(g), 0);
    maxcap4.insert(maxcap4.begin()+(g), 0);
    
    length=locations.size();
    numRicRoute=(length-2)/2;
    
    totaldist=calculatedist(Time);
    
    
    
}
 
 
 */
void Route::insert_req_pos(int req, int l, int g, std::vector<RichiesteServizio*> &ric, std::vector<std::vector<double>> &Time){
    
    
    locations.insert(locations.begin()+(l+1),ric[req]->locationP);
    locations.insert(locations.begin()+(g+2),ric[req]->locationD);
    
    Ricid.insert(Ricid.begin()+(l+1), req);
    Ricid.insert(Ricid.begin()+(g+2), req);
    type.insert(type.begin()+(l+1), 1);
    type.insert(type.begin()+(g+2), -1);
    arrival.insert(arrival.begin()+(l+1),0);
    departure.insert(departure.begin()+(l+1), 0);
    waiting.insert(waiting.begin()+(l+1), 0);
    earliest.insert(earliest.begin()+(l+1),0);
    latest.insert(latest.begin()+(l+1), 0);
    cap1.insert(cap1.begin()+(l+1), 0);
    maxcap1.insert(maxcap1.begin()+(l+1), 0);
    cap2.insert(cap2.begin()+(l+1), 0);
    maxcap2.insert(maxcap2.begin()+(l+1), 0);
    cap3.insert(cap3.begin()+(l+1), 0);
    maxcap3.insert(maxcap3.begin()+(l+1), 0);
    cap4.insert(cap4.begin()+(l+1), 0);
    maxcap4.insert(maxcap4.begin()+(l+1), 0);
    arrival.insert(arrival.begin()+(g+2),0);
    departure.insert(departure.begin()+(g+2), 0);
    
    waiting.insert(waiting.begin()+(g+2), 0);
    earliest.insert(earliest.begin()+(g+2),0);
    latest.insert(latest.begin()+(g+2), 0);
    cap1.insert(cap1.begin()+(g+2), 0);
    maxcap1.insert(maxcap1.begin()+(g+2), 0);
    cap2.insert(cap2.begin()+(g+2), 0);
    maxcap2.insert(maxcap2.begin()+(g+2), 0);
    cap3.insert(cap3.begin()+(g+2), 0);
    maxcap3.insert(maxcap3.begin()+(g+2), 0);
    cap4.insert(cap4.begin()+(g+2), 0);
    maxcap4.insert(maxcap4.begin()+(g+2), 0);
    
    length=locations.size();
    numRicRoute=(length-2)/2;
    
    totaldist=calculatedist(Time);
    
}

void Route::insert_request_on(int req, int p, int d, std::vector<RichiesteServizio*> &ric, std::vector<std::vector<double>> &Time){




    locations.insert(locations.begin()+p,ric[req]->locationP);
    locations.insert(locations.begin()+(d),ric[req]->locationD);
    Ricid.insert(Ricid.begin()+(p), req);
    Ricid.insert(Ricid.begin()+(d), req);
    type.insert(type.begin()+(p), 1);
    type.insert(type.begin()+(d), -1);

    arrival.insert(arrival.begin()+(p),0);
    departure.insert(departure.begin()+(p), 0);
    waiting.insert(waiting.begin()+(p), 0);
    earliest.insert(earliest.begin()+(p),0);
    latest.insert(latest.begin()+(p), 0);
    cap1.insert(cap1.begin()+(p), 0);
    maxcap1.insert(maxcap1.begin()+(p), 0);
    cap2.insert(cap2.begin()+(p), 0);
    maxcap2.insert(maxcap2.begin()+(p), 0);
    cap3.insert(cap3.begin()+(p), 0);
    maxcap3.insert(maxcap3.begin()+(p), 0);
    cap4.insert(cap4.begin()+(p), 0);
    maxcap4.insert(maxcap4.begin()+(p), 0);


    arrival.insert(arrival.begin()+(d),0);
    departure.insert(departure.begin()+(d), 0);
    waiting.insert(waiting.begin()+(d), 0);
    earliest.insert(earliest.begin()+(d),0);
    latest.insert(latest.begin()+(d), 0);
    cap1.insert(cap1.begin()+(d), 0);
    maxcap1.insert(maxcap1.begin()+(d), 0);
    cap2.insert(cap2.begin()+(d), 0);
    maxcap2.insert(maxcap2.begin()+(d), 0);
    cap3.insert(cap3.begin()+(d), 0);
    maxcap3.insert(maxcap3.begin()+(d), 0);
    cap4.insert(cap4.begin()+(d), 0);
    maxcap4.insert(maxcap4.begin()+(d), 0);

    length=locations.size();
    numRicRoute=(length-2)/2;

    totaldist=calculatedist(Time);

}



void Route::calculate_times(RichiesteServizio* ric, double** Time){
    
    int i;
    for (i=2; i<length; i++){
        
        
    }
    
}
bool Route::check_feasibility_P(int req,int pre,RichiesteServizio* ric, double** Time, Veicoli* veic){
    bool feas1;
    
    
    if(earliest[pre+1]<=ric[req].timewinPmax){
        
        feas1=true;
    }else{
        feas1=false;
    }
    
    if(feas1==true){
    
        if(cap3[pre+1]>veic[id].stretcher){
            if(cap2[pre+1]>veic[id].seated){
                if(cap1[pre+1]>veic[id].staff){
                    feas1=false;
                    
                }else{
                    
                }
            }else{
                if(cap2[pre+1]+cap1[pre+1]>veic[id].staff+veic[id].seated){
                    feas1=false;
                }else{
                    
                    
                }
                
            }
        }else{
            if(cap2[pre+1]>veic[id].seated){
                if(cap3[pre+1]+cap1[pre+1]>veic[id].staff+veic[id].stretcher){
                    
                    feas1=false;
                }else{}
            }else{
                if(cap2[pre+1]+cap3[pre+1]+cap1[pre+1]>veic[id].staff+veic[id].stretcher+veic[id].seated){
                    feas1=false;
                    
                }else{
                    
                }
                
                
                
            }
            
        }
        
        
        
        if(cap3[pre+1]>veic[id].stretcher){
            if(cap2[pre+1]>veic[id].seated){
                feas1=false;
            }else{
                
            }
        }else{
            if(cap2[pre+1]+cap3[pre+1]>veic[id].stretcher+veic[id].seated){
                feas1=false;
            }else{
                
            }
            
        }
        
        
        
        
        if(cap4[pre+1]>veic[id].wheelchair){
            
            feas1=false;
            
        }else{}
        
        if(cap3[pre+1]>veic[id].stretcher){
            feas1=false;
        }else{
            
        }

    
    }else{
    
        feas1=false;
    }
    
    
    
    return feas1;
}
bool Route::check_feasibility_P_tw(int req, int pre, RichiesteServizio* ric, double** Time){
    
    bool feas;
    
    
    if(earliest[pre+1]<=ric[req].timewinPmax){
        
        feas=true;
    }else{
        feas=false;
    }
    
    return feas;
}
bool Route::check_feasibility_P_cap(int req, int pre, Veicoli* veic, RichiesteServizio* ric){
    bool feas=true;
    if(cap3[pre+1]>veic[id].stretcher){
        if(cap2[pre+1]>veic[id].seated){
            if(cap1[pre+1]>veic[id].staff){
                feas=false;
                
            }else{
                
            }
        }else{
            if(cap2[pre+1]+cap1[pre+1]>veic[id].staff+veic[id].seated){
                feas=false;
            }else{
                
                
            }
            
        }
    }else{
        if(cap2[pre+1]>veic[id].seated){
            if(cap3[pre+1]+cap1[pre+1]>veic[id].staff+veic[id].stretcher){
                
                feas=false;
            }else{}
        }else{
            if(cap2[pre+1]+cap3[pre+1]+cap1[pre+1]>veic[id].staff+veic[id].stretcher+veic[id].seated){
                feas=false;
                
            }else{
                
            }
            
            
            
        }
        
    }
    
    
    
    if(cap3[pre+1]>veic[id].stretcher){
        if(cap2[pre+1]>veic[id].seated){
            feas=false;
        }else{
            
        }
    }else{
        if(cap2[pre+1]+cap3[pre+1]>veic[id].stretcher+veic[id].seated){
            feas=false;
        }else{
            
        }
        
    }
    
    
    
    
    if(cap4[pre+1]>veic[id].wheelchair){
        
        feas=false;
        
    }else{}
    
    if(cap3[pre+1]>veic[id].stretcher){
        feas=false;
    }else{
        
    }
    
    
    
    
    return feas;
}

double Route::effect_ondistance(int req, int p, int d, double** Dist, RichiesteServizio* ric){
    double dif;
    double newdist=0;
    int i;
    for(i=1;i<locations.size();i++){
        
        newdist+=Dist[locations[i-1]][locations[i]];
    }
    
    
    
    dif=newdist-totaldist;
    return dif;
}


bool Route::check_feasibillity_D_tw1(int req,int pre, RichiesteServizio* ric, double** Time){
    
    bool feas;
    
    //for(int i=0; i<length; i++){
    
     //   std::cout << latest[i] << " " ;}
   // std::cout << std::endl;
   // std::cout << earliest[pre+2] << "<=" << ric[req].timewinDmax << " &&" <<  earliest[pre+3] << "<=" << latest[pre+3] << std::endl;
    
    if(earliest[pre+2]<=ric[req].timewinDmax && earliest[pre+3]<=latest[pre+3]){
        feas=true;
    }else{
        feas=false;
    }
    
    
    return feas;
}

bool Route::check_feasibility_D_tw2( int pre, int suc, RichiesteServizio* ric, double** Time){
    
    bool feas=true;
    
    int i;
    //earliest of the pickup insertion node i
    i=pre;
    while(i<suc){
        //	std::cout << i << "\n";
        if(type[i+2]==1){
            if(earliest[i+2]<=ric[Ricid[i+2]].timewinPmax){
                feas=true;
            }else{
                feas=false;
                i=suc;
                //std::cout << i << "\n";
            }
        }else{
            if(earliest[i+2]<=ric[Ricid[i+2]].timewinDmax){
                feas=true;
            }else{
                feas=false;
                i=suc;
            }
            
        }
        i=i+1;
        //std::cout << i << "\n";
        
        
    }
    if(earliest[suc+2]>ric[Ricid[suc+2]].timewinDmax){
        feas=false;
    }
    return feas;
}

void Route::calculate_capacity(std::vector<RichiesteServizio*> &ric, std::vector<Veicoli*> &veic){
    
    cap1[0]=0;
    cap2[0]=0;
    cap3[0]=0;
    cap4[0]=0;
    
    
    for(int i=1; i<length-1; i++){
        if(type[i]==1){
            
            cap1[i]=cap1[i-1]+ric[Ricid[i]]->staff;
            cap2[i]=cap2[i-1]+ric[Ricid[i]]->seated;
            cap3[i]=cap3[i-1]+ric[Ricid[i]]->stretcher;
            cap4[i]=cap4[i-1]+ric[Ricid[i]]->wheelchair;
        }else{
            cap1[i]=cap1[i-1]-ric[Ricid[i]]->staff;
            cap2[i]=cap2[i-1]-ric[Ricid[i]]->seated;
            cap3[i]=cap3[i-1]-ric[Ricid[i]]->stretcher;
            cap4[i]=cap4[i-1]-ric[Ricid[i]]->wheelchair;
        }
    }
    cap1[length-1]=0;
    cap2[length-1]=0;
    cap3[length-1]=0;
    cap4[length-1]=0;
    
    maxcap1[length-1]=0;
    maxcap2[length-1]=0;
    maxcap3[length-1]=0;
    maxcap4[length-1]=0;
    
    for(int i=length-2; i>-1;i--){
        if(type[i]==1){
            maxcap1[i]=maxcap1[i+1]+ric[Ricid[i]]->staff;
            maxcap2[i]=maxcap2[i+1]+ric[Ricid[i]]->seated;
            maxcap3[i]=maxcap3[i+1]+ric[Ricid[i]]->stretcher;
            maxcap4[i]=maxcap4[i+1]+ric[Ricid[i]]->wheelchair;
        }else{
            if(type[i]==-1){
                maxcap1[i]=maxcap1[i+1]-ric[Ricid[i]]->staff;
                maxcap2[i]=maxcap2[i+1]-ric[Ricid[i]]->seated;
                maxcap3[i]=maxcap3[i+1]-ric[Ricid[i]]->stretcher;
                maxcap4[i]=maxcap4[i+1]-ric[Ricid[i]]->wheelchair;
            }else{
                maxcap1[i]=maxcap1[i+1];
                maxcap2[i]=maxcap2[i+1];
                maxcap3[i]=maxcap3[i+1];
                maxcap4[i]=maxcap4[i+1];
                
            
            }
        }
        if(maxcap1[i]<0){
            maxcap1[i]=0;
        }
        if(maxcap2[i]<0){
            maxcap2[i]=0;
        }
        if(maxcap3[i]<0){
            maxcap3[i]=0;
        }
        if(maxcap4[i]<0){
            maxcap4[i]=0;
        }
    }
    
}

void Route::calculate_earliest_latest(std::vector<RichiesteServizio*> &ric,  std::vector<std::vector<double>> &Time, std::vector<Veicoli*> &veic){
    
    earliest[0]=0;
    earliest[1]=Time[locations[0]][locations[1]];
    if(earliest[1]<ric[Ricid[1]]->timewinPmin){
        earliest[1]=ric[Ricid[1]]->timewinPmin;
    }

    for(int i=2; i<length-1;i++){
        RichiesteServizio& rich =*ric[Ricid[i]];
        
        earliest[i]=earliest[i-1]+ric[Ricid[i-1]]->st+Time[locations[i-1]][locations[i]];
        if(type[i]==1){
            
            if(earliest[i]<rich.timewinPmin){
                earliest[i]=rich.timewinPmin;
            }
        }else{
            if(type[i]==-1){
                
                if(earliest[i]<rich.timewinDmin){
                    earliest[i]=rich.timewinDmin;
                }}
            
        }
        
    }
    earliest[length-1]=earliest[length-2]+ric[Ricid[length-2]]->st+Time[locations[length-2]][locations[length-1]];
    
    
    latest[length-1]=veic[veh]->MaxRoute;
    for(int i=length-2; i>0;i--){
        RichiesteServizio& rich =*ric[Ricid[i]];
        latest[i]=latest[i+1]-Time[locations[i]][locations[i+1]]-rich.st;
        
        if(type[i]==1){
            if(latest[i]>rich.timewinPmax){
                latest[i]=rich.timewinPmax;
            }
        }else{if(type[i]==-1){
            if(latest[i]>rich.timewinDmax){
                latest[i]=rich.timewinDmax;
            }
        }
        }
    }
    
    latest[0]=latest[1]-Time[locations[0]][locations[1]];
    
    
   
}


bool Route::check_feasibility_capacity(int pre,int suc, Veicoli* veic){
    
    Veicoli& ve=veic[id];
    
    bool feas=true;
    int i;
    i=pre;
    while(i<suc){
        if(cap3[i+2]>ve.stretcher){
            if(cap2[i+2]>ve.seated){
                if(cap1[i+2]>ve.staff){
                    feas=false;
                    i=suc;
                    
                }else{
                    
                }
            }else{
                if(cap2[i+2]+cap1[i+2]>ve.staff+ve.seated){
                    feas=false;
                    i=suc;
                }else{
                    
                    
                }
                
            }
        }else{
            if(cap2[i+2]>ve.seated){
                if(cap3[i+2]+cap1[i+2]>ve.staff+ve.stretcher){
                    
                    feas=false;
                    i=suc;
                }else{}
            }else{
                if(cap2[i+2]+cap3[i+2]+cap1[i+2]>ve.staff+ve.stretcher+ve.seated){
                    feas=false;
                    i=suc;
                    
                }else{
                    
                }
                
                
                
            }
            
        }
        
        
        
        if(cap3[i+2]>ve.stretcher){
            if(cap2[i+2]>ve.seated){
                feas=false;
                i=suc;
            }else{
                
            }
        }else{
            if(cap2[i+2]+cap3[i+2]>ve.stretcher+ve.seated){
                feas=false;
                i=suc;
            }else{
                
            }
            
        }
        
        
        
        
        if(cap4[i+2]>ve.wheelchair){
            
            feas=false;
            i=suc;
            
        }else{}
        
        if(cap3[i+2]>ve.stretcher){
            feas=false;
            i=suc;
        }else{
            
        }
        
        i=i+1;
    }
    if(cap3[suc+2]>ve.stretcher){
        if(cap2[suc+2]>ve.seated){
            if(cap1[suc+2]>ve.staff){
                feas=false;
                
                
            }else{
                
            }
        }else{
            if(cap2[suc+2]+cap1[suc+2]>ve.staff+ve.seated){
                feas=false;
                
            }else{
                
                
            }
            
        }
    }else{
        if(cap2[suc+2]>ve.seated){
            if(cap3[suc+2]+cap1[suc+2]>ve.staff+ve.stretcher){
                
                feas=false;
                
            }else{}
        }else{
            if(cap2[suc+2]+cap3[suc+2]+cap1[suc+2]>ve.staff+ve.stretcher+ve.seated){
                feas=false;
                
                
            }else{
                
            }
            
            
            
        }
        
    }
    
    
    
    if(cap3[suc+2]>ve.stretcher){
        if(cap2[suc+2]>ve.seated){
            feas=false;
            
        }else{
            
        }
    }else{
        if(cap2[suc+2]+cap3[suc+2]>ve.stretcher+ve.seated){
            feas=false;
            
        }else{
            
        }
        
    }
    
    
    
    
    if(cap4[suc+2]>ve.wheelchair){
        
        feas=false;
        
        
    }else{}
    
				if(cap3[suc+2]>ve.stretcher){
                    feas=false;
                    
                }else{
                    
                }
    
    
    return feas;
}

bool Route::ridetime_check(int suc, RichiesteServizio* ric, int MaxRideTime){
    bool feas=true;
    if(earliest[suc+2]-ric[Ricid[suc+2]].timewinPmax-ric[Ricid[suc+2]].st>MaxRideTime){
        
        feas=false;
        
    }
   // for(int i=0; i<length; i++){
    //    std::cout<< earliest[i] << " " ;
   // }std::cout<< std::endl;
    
    return feas;
}

bool Route::Eightstepevaluationscheme(std::vector<std::vector<double>> &Time, std::vector<RichiesteServizio*> &ric, std::vector<Veicoli*> &veic){
    int i, req;
    double F;
    bool feas=true;
    double v1=0, v2=0, v3=0;
    
    for(i=0; i<length; i++){
        arrival[i]=0;
        departure[i]=0;
        waiting[i]=0;
    }
    
    
    
    for(i=1; i<length-1; i++){
        
        arrival[i]=departure[i-1]+Time[locations[i-1]][locations[i]];
        
        if(type[i]==1){
            if(arrival[i]-ric[Ricid[i]]->timewinPmin>=0 && arrival[i]-ric[Ricid[i]]->timewinPmax<=0){
                departure[i]=arrival[i]+ric[Ricid[i]]->st;
                waiting[i]=0;
                
            }else{if(arrival[i]-ric[Ricid[i]]->timewinPmin<0){
                
                departure[i]=ric[Ricid[i]]->timewinPmin+ric[Ricid[i]]->st;
                waiting[i]=departure[i]-ric[Ricid[i]]->st-arrival[i];
            }else{
                departure[i]=arrival[i]+ric[Ricid[i]]->st;
                
                v1+=arrival[i]-ric[Ricid[i]]->timewinPmax;
                //Viol+=Arrival[j]-B[int(Sol1[i][j][1])][4];
            }
            }
        }else{//if it is a delivery node
            if(arrival[i]-ric[Ricid[i]]->timewinDmin>=0 && arrival[i]-ric[Ricid[i]]->timewinDmax<=0){
                departure[i]=arrival[i]+ric[Ricid[i]]->st;
                waiting[i]=0;
            }else{ if(arrival[i]-ric[Ricid[i]]->timewinDmin<0){
                departure[i]=ric[Ricid[i]]->timewinDmin+ric[Ricid[i]]->st;
                waiting[i]=departure[i]-ric[Ricid[i]]->st-arrival[i];
            }else{
                departure[i]=arrival[i]+ric[Ricid[i]]->st;
                
                v1+=arrival[i]-ric[Ricid[i]]->timewinDmax;
                //Viol+=Arrival[j]-B[int(Sol1[i][j][1])][6];
            }}
        }
        
    }
    arrival[length-1]=departure[length-2]+Time[locations[length-1]][locations[length-2]];
    departure[length-1]=0;
    waiting[length-1]=0;
    
    //for(i=0;i<length;i++){
    //	std::cout << arrival[i] << " " << departure[i] << " " << waiting[i] << "\n" ;
    //}
    
    
    if(v1>0){
        //go to step 8
        //std:: cout << "false \n";
    }else{
        //std:: cout << "true \n";
        //calculate the forward slack
        F=calculateF(1, ric);
        
        double G;
        G=0;
        
        for(i=1; i<length; i++){
            G+=waiting[i];
        }
        if(F<=G){
            G=F;
        }
        //std:: cout << "G is" << G << "\n" ;
        arrival[1]=(departure[1]-ric[Ricid[1]]->st)+G;
        departure[0]=arrival[1]-Time[locations[0]][locations[1]];
        v1=0;
        for(i=1; i<length-1; i++){//hasta menos uno porque not iene niungun sentido que se calcule lo del deposito
            arrival[i]=departure[i-1]+Time[locations[i]][locations[i-1]];
            if(type[i]==1){ //if it is a pickup node
                if(arrival[i]-ric[Ricid[i]]->timewinPmin>=0 && arrival[i]-ric[Ricid[i]]->timewinPmax<=0){
                    departure[i]=arrival[i]+ric[Ricid[i]]->st;
                    waiting[i]=0;
                }else{if(arrival[i]-ric[Ricid[i]]->timewinPmin<0){
                    
                    departure[i]=ric[Ricid[i]]->timewinPmin+ric[Ricid[i]]->st;
                    waiting[i]=departure[i]-ric[Ricid[i]]->st-arrival[i];
                }else{
                    departure[i]=arrival[i]+ric[Ricid[i]]->st;
                    
                    v1+=arrival[i]-ric[Ricid[i]]->timewinPmax;
                    //Viol+=Arrival[j]-B[int(Sol1[i][j][1])][4];
                }}
            }else{//if it is a delivery node
                if(arrival[i]-ric[Ricid[i]]->timewinDmin>=0 && arrival[i]-ric[Ricid[i]]->timewinDmax<=0){
                    departure[i]=arrival[i]+ric[Ricid[i]]->st;
                    waiting[i]=0;
                }else{ if(arrival[i]-ric[Ricid[i]]->timewinDmin<0){
                    departure[i]=ric[Ricid[i]]->timewinDmin+ric[Ricid[i]]->st;
                    waiting[i]=departure[i]-ric[Ricid[i]]->st-arrival[i];
                }else{
                    departure[i]=arrival[i]+ric[Ricid[i]]->st;
                    
                    v1+=arrival[i]-ric[Ricid[i]]->timewinDmax;
                    //Viol+=Arrival[j]-B[int(Sol1[i][j][1])][6];
                }}
            }
        }
        arrival[length-1]=departure[length-2]+Time[locations[length-1]][locations[length-2]];
        departure[length-1]=0;
        waiting[length-1]=0;
        
        //for(i=0;i<length;i++){
        //	std::cout << arrival[i] << " " << departure[i] << " " << waiting[i] << "\n" ;
        //}
        
        
        
        v3=0;
        for(int ii=0; ii<length; ii++){
            
            if(type[ii]==1){
                
                req=Ricid[ii];
                for(int jj=ii+1;jj<length;jj++){
                    
                    if(type[jj]==-1){ //es nodo de delivery
                        if(Ricid[jj]==req){
                            
                            if(departure[jj]-ric[req]->st-departure[ii]-ric[req]->RideTime>0){
                                v3+=(departure[jj]-ric[req]->st-departure[ii]-ric[req]->RideTime);
                                
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
            for(i=1; i<length-1; i++){
                if(type[i]==1){
                    F=calculateF(i, ric);
                    //set waitings
                    G=0;
                    for(int p=i+1; p<length; p++){
                        G+=waiting[p];
                    }
                    if(F<G){
                        G=F;
                    }else{
                        
                    }
                    //calculate the others
                    departure[i]=departure[i]+G;
                    waiting[i]=departure[i]-ric[Ricid[i]]->st-arrival[i];
                    
                    
                    for(int kk=i+1;kk<length-1;kk++){
                        arrival[kk]=departure[kk-1]+Time[locations[kk-1]][locations[kk]];
                        if(type[kk]==1){ //if it is a pickup node
                            if(arrival[kk]-ric[Ricid[kk]]->timewinPmin>=0 && arrival[kk]-ric[Ricid[kk]]->timewinPmax<=0){
                                departure[kk]=arrival[kk]+ric[Ricid[kk]]->st;
                                waiting[kk]=0;
                            }else{ if(arrival[kk]-ric[Ricid[kk]]->timewinPmin<0){
                                
                                departure[kk]=ric[Ricid[kk]]->timewinPmin+ric[Ricid[kk]]->st;
                                waiting[kk]=departure[kk]-ric[Ricid[kk]]->st-arrival[kk];
                            }else{
                                
                                departure[kk]=arrival[kk]+ric[Ricid[kk]]->st;
                                v1+=arrival[kk]-ric[Ricid[kk]]->timewinPmax;
                                //Viol+=Arrival[kk]-B[int(Sol1[i][kk][1])][4];
                            }}
                        }else{ //if it is a delivery node
                            //delivery home
                            
                            //if it is a delivery node
                            if(arrival[kk]-ric[Ricid[kk]]->timewinDmin>=-0 && arrival[kk]-ric[Ricid[kk]]->timewinDmax<=0  ){
                                departure[kk]=arrival[kk]+ric[Ricid[kk]]->st;
                                waiting[kk]=0;
                            }else{ if(arrival[kk]-ric[Ricid[kk]]->timewinDmin<-0){
                                departure[kk]=ric[Ricid[kk]]->timewinDmin+ric[Ricid[kk]]->st;
                                waiting[kk]=departure[kk]-ric[Ricid[kk]]->st-arrival[kk];
                            }else{
                                
                                departure[kk]=arrival[kk]+ric[Ricid[kk]]->st;
                                v1+=arrival[kk]-ric[Ricid[kk]]->timewinDmax;
                                //Viol+=Arrival[kk]-B[int(Sol1[i][kk][1])][6];
                            }}
                            
                            
                        }
                        
                        
                    }
                    
                    arrival[length-1]=departure[length-2]+Time[locations[length-2]][locations[length-1]];
                    
                }else{
                    
                }
                
                
                
                
                
            }
            
            
            arrival[1]=departure[1]-ric[Ricid[1]]->st;
            departure[0]=arrival[1]-Time[locations[0]][locations[1]];
            
            
            
            
        }
        
        
    }
    
    v3=0;
    for(int ii=0; ii<length; ii++){
        
        if(type[ii]==1){
            
            req=Ricid[ii];
            for(int jj=ii+1;jj<length;jj++){
                
                if(type[jj]==-1){ //es nodo de delivery
                    if(Ricid[jj]==req){
                        if(departure[jj]-ric[req]->st-departure[ii]-ric[req]->RideTime>0){
                            
                            
                            v3+=(departure[jj]-ric[req]->st-departure[ii]-ric[req]->RideTime);
                        }else{
                        }
                    }
                }
                
            }
        }
        
    }
    
    
    if(arrival[length-1]-departure[0]-veic[veh]->MaxRoute>0){
        
        //V[i][1]=arrival[length-1]-departure[0]-MaxTimeRoute;
        v2=arrival[length-1]-departure[0]-veic[veh]->MaxRoute;
        
        
    }else{v2=0;}
    
    
    if(v1>0 || v2>0 || v3>0){
        feas=false;
    }
    
    
    return feas;
}

double Route::calculateF(int a, std::vector<RichiesteServizio*>  ric){
    double F, F1;
    F=DBL_MAX;
    //double* F1;
    int i, j, l;
    //F1=new double[length-a];
    double sum, M1, M2, timetravelled, M;
    for(i=a; i<length; i++){
        sum=0;
        
        for(j=a+1; j<i+1; j++){
            
            sum+=waiting[j];
            
        }
        
        
        if(type[i]==0){
            
            M1=1440-(departure[i]-ric[Ricid[i-1]]->st);
            
            M2=1440-0;
            
        }else{
            
            if(type[i]==1){//es de pickup
                
                M1=ric[Ricid[i]]->timewinPmax-(departure[i]-ric[Ricid[i]]->st);
                M2=DBL_MAX;
            }else{
                M1=ric[Ricid[i]]->timewinDmax-(departure[i]-ric[Ricid[i]]->st);
                bool yes=false;
                timetravelled=0;
                for(l=i-1; l>0; l--){
                    
                    if(Ricid[i]==Ricid[l]){
                        timetravelled=arrival[i]-departure[l];
                        if(l<a){yes=true;}
                    }
                }
                if(yes==true){
                    M2=ric[Ricid[i]]->RideTime-timetravelled;
                }else{M2=FLT_MAX;}
            }}
        if(M1<=0){M1=0;}
        if(M2<=0){M2=0;}
        
        if(M1<=M2){
            M=M1;
        }else{M=M2;}
        F1=sum+M;
        
        if(F>F1){
            
            F=F1;
            
        }
        
        
        
    }
    
    
    
    
    return F;
}
double Route::calculatedist_2(std::vector<std::vector<double>> &Time, int i, int j, int z, int u){
    
    double dist=0;
    int l, loc_p;
    for(int ii=0;ii<u-1; ii++){
        
        dist+=Time[locations[ii]][locations[ii+1]];
        
    }
    
    int l_p=locations[u-1];
    for(int ii=0;ii<3;ii++){
    
        if(ii==i){
            dist+=Time[l_p][locations[u]];
            l_p=locations[u];
            
        }else{
            if(ii==j){
                dist+=Time[l_p][locations[u+1]];
                l_p=locations[u+1];
                
            }else{
                dist+=Time[l_p][locations[u+2]];
                l_p=locations[u+2];

            }
        }
    }
    
    for(int ii=u+2; ii<length-1; ii++){
        dist+=Time[l_p][locations[ii+1]];
        l_p=locations[ii+1];
    
    }
    
    return dist;
}
double Route::calculatedist(std::vector<std::vector<double>> &Time){
    
    double dist=0;
    
    for(int i=0;i<length-1; i++){
        
        dist+=Time[locations[i]][locations[i+1]];
        
    }
    
    
    return dist;
}

void Route::CopyRoute(Route* route){
    
    id=route->id;
    veh=route->veh;
    length=route->length;
    numRicRoute=route->numRicRoute;
    totaldist=route->totaldist;

    locations=route->locations;
    Ricid=route->Ricid;
    type=route->type;
    arrival=route->arrival;
    departure=route->departure;
    waiting=route->waiting;
    earliest=route->earliest;
    latest=route->latest;
    cap1=route->cap1;
    maxcap1=route->maxcap1;
    cap2=route->cap2;
    maxcap2=route->maxcap2;
    cap3=route->cap3;
    maxcap3=route->maxcap3;
    cap4=route->cap4;
    maxcap4=route->maxcap4;
    
}


void Route::remove_pos(int p){
    
    locations.erase(locations.begin()+p);
    Ricid.erase(Ricid.begin()+p);
    type.erase(type.begin()+p);
    //length=locations.size();
    arrival.erase(arrival.begin()+p);
    departure.erase(departure.begin()+p);
    waiting.erase(waiting.begin()+p);
    earliest.erase(earliest.begin()+p);
    latest.erase(latest.begin()+p);
    cap1.erase(cap1.begin()+p);
    maxcap1.erase(maxcap1.begin()+p);
    cap2.erase(cap2.begin()+p);
    maxcap2.erase(maxcap2.begin()+p);
    cap3.erase(cap3.begin()+p);
    maxcap3.erase(maxcap3.begin()+p);
    cap4.erase(cap4.begin()+p);
    maxcap4.erase(maxcap4.begin()+p);
    
    
    
    
}


void Route::rewriteid(int a){
    
    id=a;
    
    
    
}


void Route::count_request(int* E){
    
    int req;
    int k=0;
    for(int i=0; i<length; i++){
        
        if(type[i]==1){
            req=Ricid[i];
            for(int j=i+1; j<length; j++){
                if(type[j]==-1){
                    if(Ricid[j]==req){
                        
                        E[k]=req;
                        k+=1;
                    }
                }
            }
        }
    }
    
    
    
    

}


int Route::find_pickup(int m){
    
    int p=0;
    for(int i=0; i<length; i++){
        if(Ricid[i]==m){
            if(type[i]==1){
                p=i;
            }
        }
        
        
    }
    return p;
}

int Route::find_delivery(int m){
    
    int p=0;
    for(int i=0; i<length; i++){
        if(Ricid[i]==m){
            if(type[i]==-1){
                p=i;
            }
        }
        
        
    }
    return p;
}

void Route::delete_req(int m,  std::vector<std::vector<double>>  &Time){
    int a, b;
    a=find_pickup(m);
    
    b=find_delivery(m);
    // std::cout << " cheking " << std::endl;
    //std::cout << a << " " << b << std::endl;
    remove_pos(b);
    remove_pos(a);
    length=locations.size();
    numRicRoute=(length-2)/2;
    totaldist=calculatedist(Time);
}
int Route::count_req3(int i){
    int a=0;
    
    for(int j=i; j<i+3; j++){
        if(type[j]==1){
            a+=1;
            
        }
        
    }
    
    return a;
}
void::Route::reverse_order(int u, int i, int j, int z){

    //the one in pos i goes in 0, the one in pos j goes in 1, the one in pos z goes in 2

    int a1=locations[u];
    int a2=locations[u+1];
    int a3=locations[u+2];

    int b1=Ricid[u];
    int b2=Ricid[u+1];
    int b3=Ricid[u+2];

    int c1=type[u];
    int c2=type[u+1];
    int c3=type[u+2];


    if(i==0){

    }else{
        if(i==1){
            locations[u]=a2;
            Ricid[u]=b2;
            type[u]=c2;

        }else{
            locations[u]=a3;
            Ricid[u]=b3;
            type[u]=c3;

        }

    }

    if(j==0){
        locations[u+1]=a1;
        Ricid[u+1]=b1;
        type[u+1]=c1;

    }else{
        if(j==1){

        }else{
            locations[u+1]=a3;
            Ricid[u+1]=b3;
            type[u+1]=c3;
        }

    }
    if(z==0){
        locations[u+2]=a1;
        Ricid[u+2]=b1;
        type[u+2]=c1;

    }else{if(z==1){
            locations[u+2]=a2;
            Ricid[u+2]=b2;
            type[u+2]=c2;

        }else{}}


}


void Route::change_order(int u, int i , int j, int z){



    int a1=locations[u];
    int a2=locations[u+1];
    int a3=locations[u+2];

    int b1=Ricid[u];
    int b2=Ricid[u+1];
    int b3=Ricid[u+2];

    int c1=type[u];
    int c2=type[u+1];
    int c3=type[u+2];




    for(int ii=0; ii<3; ii++){
        if(ii==i){
            locations[u+ii]=a1;
            Ricid[u+ii]=b1;
            type[u+ii]=c1;

        }else{
            if(j==ii){
                locations[u+ii]=a2;
                Ricid[u+ii]=b2;
                type[u+ii]=c2;


            }else{
                locations[u+ii]=a3;
                Ricid[u+ii]=b3;
                type[u+ii]=c3;

            }


        }

    }

    
}
bool Route::delivery_before_pickup_2(int u, int i, int j, int z){
    bool feas=true;
    int r, r1;
    
    int a;
    
    for(int ii=0; ii<3;ii++){
        if(feas==true){
        if(i==ii){
            r=Ricid[u];
            a=type[u];
        }else{if(j==ii){
            r=Ricid[u+1];
            a=type[u+1];
            }else{
            
            r=Ricid[u+2];
            a=type[u+2];
            }}
        
        if(a==-1){
            for(int jj=ii+1; jj<3; jj++){
                if(i==jj){
                    r1=Ricid[u];
                    
                }else{if(j==jj){
                    r1=Ricid[u+1];
                    
                }else{
                    r1=Ricid[u+2];
                    
                }}
                
                
                if(r1==r){
                    feas=false;
                }
                
            }
        }
    }
    }
    return feas;
}

bool Route::delivery_before_pickup(int u){
    bool feas=true;
    int r;
    for(int i=0; i<3;i++){
        r=Ricid[u+i];
        if(type[u+i]==-1){
            for(int j=i+1; j<3; j++){
                if(Ricid[u+j]==r){
                    feas=false;
                }
                
            }
        }
    }
    
    return feas;
}



bool Route::check_feas1(int u, std::vector<RichiesteServizio*> &ric, std::vector<std::vector<double>> & Time, int i, int j, int z){
    bool f;
    
    //calculate the earliest
    //calculate the latest
    double ear_0=0;
    double ear;
    
    //only calculate desde earl
    
    ear_0=earliest[u-1];

    int ii_0=u-1;
    for(int ii=0; ii<3; ii++){
        
        if(ii==i){
            
            if(type[ii_0]==0){
                ear=+Time[locations[ii_0]][locations[u]];
                
            }else{ear=ear_0+ric[Ricid[ii_0]]->st+Time[locations[ii_0]][locations[u]];}
            if(type[u]==1){
                
                if(ear<ric[Ricid[u]]->timewinPmin){
                    ear=ric[Ricid[u]]->timewinPmin;
                }
            }else{
                if(type[u]==-1){
                    
                    if(ear<ric[Ricid[u]]->timewinDmin){
                        ear=ric[Ricid[u]]->timewinDmin;
                    }}
                
            }
            ii_0=u;
            
            ear_0=ear;
            
        }else{if(ii==j){
            if(type[ii_0]==0){
                ear=+Time[locations[ii_0]][locations[u+1]];
                
            }else{ear=ear_0+ric[Ricid[ii_0]]->st+Time[locations[ii_0]][locations[u+1]];}
            if(type[u+1]==1){
                
                if(ear<ric[Ricid[u+1]]->timewinPmin){
                    ear=ric[Ricid[u+1]]->timewinPmin;
                }
            }else{
                if(type[u+1]==-1){
                    
                    if(ear<ric[Ricid[u+1]]->timewinDmin){
                        ear=ric[Ricid[u+1]]->timewinDmin;
                    }}
                
            }
            
            
            ear_0=ear;
            
            ii_0=u+1;
        
        }else{
            if(type[ii_0]==0){
                ear=+Time[locations[ii_0]][locations[u+2]];
                
            }else{ear=ear_0+ric[Ricid[ii_0]]->st+Time[locations[ii_0]][locations[u+2]];}
            if(type[u+2]==1){
                
                if(ear<ric[Ricid[u+2]]->timewinPmin){
                    ear=ric[Ricid[u+2]]->timewinPmin;
                }
            }else{
                if(type[u+2]==-1){
                    
                    if(ear<ric[Ricid[u+2]]->timewinDmin){
                        ear=ric[Ricid[u+2]]->timewinDmin;
                    }}
                
            }
            
            
            ear_0=ear;
            
            ii_0=u+2;
        
        
        }}
    
    
    
    }
    
    ear=ear_0+ric[Ricid[ii_0]]->st+Time[locations[ii_0]][locations[u+3]];
        if(type[u+3]==1){
            if(ear<ric[Ricid[u+3]]->timewinPmin){
                ear=ric[Ricid[u+3]]->timewinPmin;
            }

        }else{

            if(type[u+3]==-1){
                if(ear<ric[Ricid[u+3]]->timewinDmin){
                    ear=ric[Ricid[u+3]]->timewinDmin;
                }
            }
        }

        if(ear<=latest[u+3]){
            f=true;
        }else{
            f=false;
        }


        /*if(type[ii_0]==1){
            if(ear_0<=ric[Ricid[ii_0]]->timewinPmax){
                f=true;

            }else{
                f=false;

            }
        }else{
            if(ear_0<=ric[Ricid[ii_0]]->timewinDmax){
                f=true;

            }else{
                f=false;

            }


        }*/
    return f;
    
}



bool Route::check_cap(int u, std::vector<Veicoli*> &veic, int i, int j, int z, std::vector<RichiesteServizio*> &ric){
    

    bool f=true;

    int ii;

    int ii_c;
    int c1, c2, c3, c4;
    int c10, c20, c30, c40;
    c10=cap1[u-1];
    c20=cap2[u-1];
    c30=cap3[u-1];
    c40=cap4[u-1];
    for(ii=0;ii<3;ii++){
        if(f==true){
            if(ii==i){
                ii_c=u;
            
            }else{
                if(ii==j){
                    ii_c=u+1;
                    
                }else{
                    ii_c=u+2;
                
                }
            }
            if(type[ii_c]==1){
                c1=c10+ric[Ricid[ii_c]]->staff;
                c2=c20+ric[Ricid[ii_c]]->seated;
                c3=c30+ric[Ricid[ii_c]]->stretcher;
                c4=c40+ric[Ricid[ii_c]]->wheelchair;
                }else{
                    c1=c10-ric[Ricid[ii_c]]->staff;
                    c2=c20-ric[Ricid[ii_c]]->seated;
                    c3=c30-ric[Ricid[ii_c]]->stretcher;
                    c4=c40-ric[Ricid[ii_c]]->wheelchair;
                }

            c10=c1;
            c20=c2;
            c30=c3;
            c40=c4;

        if(c3>veic[veh]->stretcher){
            if(c2>veic[veh]->seated){
                if(c1>veic[veh]->staff){
                    f=false;
                    
                }
            }else{
                if(c2+c1>veic[veh]->staff+veic[veh]->seated){
                    f=false;
                    
                }
                
            }
        }else{
            if(c2>veic[veh]->seated){
                if(c3+c1>veic[veh]->staff+veic[veh]->stretcher){
                    
                    f=false;
                    
                }
            }else{
                if(c2+c3+c1>veic[veh]->staff+veic[veh]->stretcher+veic[veh]->seated){
                    f=false;
                    
                    
                }
                
                
                
            }
            
        }
        
        
        
        if(c3>veic[veh]->stretcher){
            if(c2>veic[veh]->seated){
                f=false;
                
            }
        }else{
            if(c2+c3>veic[veh]->stretcher+veic[veh]->seated){
                f=false;
                
            }else{
                
            }
            
        }
        
        
        
        
        if(c4>veic[veh]->wheelchair){
            
            f=false;
            
            
        }
        
        if(c3>veic[veh]->stretcher){
            f=false;
            
        }
        
        }
    }
    
    
    
    
    
    
    
    return f;
}


void Route::change_vehicle(int a, Veicoli* veic){
    
    
    locations[0]=veic[a].location;
    locations[length-1]=veic[a].location;
    
    id=a;
    
    return;

}

double Route::new_dist_without_req(int req, std::vector<std::vector<double>> & MatTemp){
    int a, b;
    a=find_pickup(req);
    b=find_delivery(req);
    double cost=0;
    
    for(int i=0;i<a-1; i++){
        cost+=MatTemp[locations[i]][locations[i+1]];
    }
    for(int i=a+1; i<b-1; i++ ){
        cost+=MatTemp[locations[i]][locations[i+1]];
    }
    
    for(int i=b+1;i<length-1;i++){
        cost+=MatTemp[locations[i]][locations[i+1]];
    
    }
    if(a+1==b){
        cost+=MatTemp[locations[a-1]][locations[b+1]];
    }else{
        cost+=MatTemp[locations[a-1]][locations[a+1]];
        cost+=MatTemp[locations[b-1]][locations[b+1]];
    }
    return cost;
}

bool Route::capacity_feasibility2(int p1,int req, std::vector<Veicoli*> &veic, std::vector<RichiesteServizio*>  &ric){

    bool feas=true;
    
    int c1=0, c2=0, c3=0, c4=0;
    
    c1=cap1[p1-1]+ric[req]->staff;
    c2=cap2[p1-1]+ric[req]->seated;
    c3=cap3[p1-1]+ric[req]->stretcher;
    c4=cap4[p1-1]+ric[req]->wheelchair;
    
    
    if(c3>veic[veh]->stretcher){
        if(c2>veic[veh]->seated){
            if(c1>veic[veh]->staff){
                feas=false;
                
            }
        }else{
            if(c2+c1>veic[veh]->staff+veic[veh]->seated){
                feas=false;
            }
        }
    }else{
        if(c2>veic[veh]->seated){
            if(c3+c1>veic[veh]->staff+veic[veh]->stretcher){
                
                feas=false;
            }
        }else{
            if(c2+c3+c1>veic[veh]->staff+veic[veh]->stretcher+veic[veh]->seated){
                feas=false;
                
            }
            
            
            
        }
        
    }
    
    
    
    if(c3>veic[veh]->stretcher){
        if(c2>veic[veh]->seated){
            feas=false;
        }
    }else{
        if(c2+c3>veic[veh]->stretcher+veic[veh]->seated){
            feas=false;
        }
        
    }
    
    
    
    
    if(c4>veic[veh]->wheelchair){
        
        feas=false;
        
    }
    
    if(c3>veic[veh]->stretcher){
        feas=false;
    }
    
    
    
    return feas;
}
bool Route::capacity_P_feasibility(int l, int req, std::vector<Veicoli*> &veic, std::vector<RichiesteServizio*>  &ric){
    bool feas=true;
    int c1=0, c2=0, c3=0, c4=0;
    
    c1=cap1[l]+ric[req]->staff;
    c2=cap2[l]+ric[req]->seated;
    c3=cap3[l]+ric[req]->stretcher;
    c4=cap4[l]+ric[req]->wheelchair;
    
    
    
    if(c3>veic[veh]->stretcher){
        if(c2>veic[veh]->seated){
            if(c1>veic[veh]->staff){
                feas=false;
                
            }else{
                
            }
        }else{
            if(c2+c1>veic[veh]->staff+veic[veh]->seated){
                feas=false;
            }else{
                
                
            }
            
        }
    }else{
        if(c2>veic[veh]->seated){
            if(c3+c1>veic[veh]->staff+veic[veh]->stretcher){
                
                feas=false;
            }else{}
        }else{
            if(c2+c3+c1>veic[veh]->staff+veic[veh]->stretcher+veic[veh]->seated){
                feas=false;
                
            }else{
                
            }
            
            
            
        }
        
    }
    
    
    
    if(c3>veic[veh]->stretcher){
        if(c2>veic[veh]->seated){
            feas=false;
        }
    }else{
        if(c2+c3>veic[veh]->stretcher+veic[veh]->seated){
            feas=false;
        }
        
    }
    
    
    
    
    if(c4>veic[veh]->wheelchair){
        
        feas=false;
        
    }
    
    if(c3>veic[veh]->stretcher){
        feas=false;
    }

    
    return feas;
}
bool Route::tw_P_feasibility_2(std::vector<RichiesteServizio*> &ric, std::vector<std::vector<double>>&Time, int p1, int req){

    bool feas;
    RichiesteServizio& rich=*ric[req];
    
    double earl;
    if(p1==1){
        earl=Time[locations[0]][rich.locationP];
    }else{
        earl=earliest[p1-1]+ric[Ricid[p1-1]]->st+Time[locations[p1-1]][rich.locationP];
    }
    
    if(earl<rich.timewinPmin){
        earl=rich.timewinPmin;
    }
    
    if(earl<=rich.timewinPmax){
        feas=true;
    }else{
        feas=false;
    }
    
    return feas;
}
bool Route::tw_P_feasibility(std::vector<RichiesteServizio*> &ric, std::vector<std::vector<double>>&Time, int l, int req){
    
    bool feas;
    
    RichiesteServizio& rich=*ric[req];
    
    double earl;
    if(l==0){
        earl=Time[locations[0]][rich.locationP];
    }else{
        earl=earliest[l]+ric[Ricid[l]]->st+Time[locations[l]][rich.locationP];
    }
    
    if(earl<rich.timewinPmin){
        earl=rich.timewinPmin;
    }
    
    if(earl<=rich.timewinPmax){
        feas=true;
    }else{
        feas=false;
    }
    
    return feas;
}
void Route::calcrid(int* rid,int l, int g, int req ){

    rid[0]=req;
    for(int i=0; i<g-l; i++){
        rid[i+1]=Ricid[l+i+1];
    }
    rid[g-l+1]=req;
  //  std::cout<<" Rid  " ;
  //  for(int i=0;i<g-l+2; i++){
   //     std::cout<<rid[i]<<" ";
        
   // }std::cout<<std::endl;
    //return rid;
}

void  Route::calctyp(int* typ,int l, int g){
    typ[0]=1;
    for(int i=0; i<g-l; i++){
        typ[i+1]=type[l+i+1];
    }
    typ[g-l+1]=-1;
    //std::cout<<" Typ  " ;
    //for(int i=0;i<g-l+2; i++){
     //   std::cout<<typ[i]<<" ";
        
   // }std::cout<<std::endl;
    //return typ;
}

void  Route::calcloc(int * loc, int req, std::vector<RichiesteServizio*> &ric, int l, int g){
    loc[0]=ric[req]->locationP;
    for(int i=0; i<g-l; i++){
        loc[i+1]=locations[l+i+1];
    }
    loc[g-l+1]=ric[req]->locationD;
   // std::cout<<" LOC  " ;
  //  for(int i=0;i<g-l+2; i++){
   //     std::cout<<loc[i]<<" ";
        
   // }std::cout<<std::endl;
    // return loc;
}

void  Route::calcearl(double* earl, int req, std::vector<RichiesteServizio*> &ric, std::vector<std::vector<double>> &Time, int l, int g, int* loc, int* rid, int* typ ){
    

    double earl0;
    RichiesteServizio& rich=*ric[req];
    if(l==0){
        earl0=Time[locations[0]][rich.locationP];
    }else{
        earl0=earliest[l]+ric[Ricid[l]]->st+Time[locations[l]][rich.locationP];
    }
    
    if(earl0<rich.timewinPmin){
        earl0=rich.timewinPmin;
    }
    
    earl[0]=earl0;
    for(int i=1; i<g-l+1; i++){
        
        RichiesteServizio& rich1=*ric[rid[i]];
        if(typ[i-1]==0){
            earl[i]=earl[i-1]+Time[loc[i-1]][loc[i]];
        }else{
            earl[i]=earl[i-1]+ric[rid[i-1]]->st+Time[loc[i-1]][loc[i]];
        }
        if(typ[i]==1){
            if(earl[i]<rich1.timewinPmin){
                earl[i]=rich1.timewinPmin;
            }
            
        }else{
            
            if(typ[i]==-1){
                if(earl[i]<rich1.timewinDmin){
                    earl[i]=rich1.timewinDmin;
                }
            }
        }
        
    }
    earl[g-l+1]=earl[g-l]+ric[rid[g-l]]->st+Time[loc[g-l]][loc[g-l+1]];
    if(earl[g-l+1]<rich.timewinDmin){
        earl[g-l+1]=rich.timewinDmin;
    }


    
    
   // return earl;
}

bool Route::ridetime_feas_D(int g, int l, std::vector<RichiesteServizio*> &ric, int req, std::vector<std::vector<double>> &Time, double* earl){
    
    
    bool feas;
    //puedo calcular desde earl l
       feas=true;
    if(earl[g-l+1]-ric[req]->timewinPmax-ric[req]->st>ric[req]->RideTime){
        feas=false;
    }
    
   // std::cout<<" earl  " ;
    //  for(int i=0;i<g-l+2; i++){
    //     std::cout<<earl[i]<<" ";
    
   //  }std::cout<<std::endl;
    
    return feas;
}


bool check_feas_D_tw2(double* earl, int* typ, int* rid, int l, int g, std::vector<RichiesteServizio*> &ric){
    int i=0;
    bool feas=true;
    int a=g-l+2;
    while(i<a-2){
    
        if(typ[i+1]==1){
            if(earl[i+1]<=ric[rid[i+1]]->timewinPmax){
                feas=true;
            
            }else{
                feas=false;
                i=a-2;
            }
        
        
        }else{
            if(earl[i+1]<=ric[rid[i+1]]->timewinDmax){
                feas=true;
                
            }else{
                feas=false;
                i=a-2;
            }

        
        
        
        }
    i=i+1;
    }
    
    if(earl[a-1]>ric[rid[a-1]]->timewinDmax){
        feas=false;
    }
    

    return feas;

}

bool Route::check_cap_from(int l, int g,std::vector<Veicoli*> &veic, int req, std::vector<RichiesteServizio*> &ric){
    bool feas=true;
    
    int i;
    
    Veicoli& ve=*veic[veh];
   
    
    int c10, c20, c30, c40;
    int c1, c2, c3, c4;
    c10=cap1[l]+ric[req]->staff;
    c20=cap2[l]+ric[req]->seated;
    c30=cap3[l]+ric[req]->stretcher;
    c40=cap4[l]+ric[req]->wheelchair;
    i=1;
    //std::cout<< g << l << std::endl;
    
    while(i<g-l+1){
        if(feas==true){
        int r=Ricid[i+l];
        //std::cout << r << std::endl;
        
        if(type[i+l]==1){
            //std::cout <<type[i+l]<< std::endl;
            c1=c10+ric[r]->staff;
            c2=c20+ric[r]->seated;
            c3=c30+ric[r]->stretcher;
            c4=c40+ric[r]->wheelchair;
        }else{
            c1=c10-ric[r]->staff;
            c2=c20-ric[r]->seated;
            c3=c30-ric[r]->stretcher;
            c4=c40-ric[r]->wheelchair;
        
        }
            
        if(c3>ve.stretcher){
            if(c2>ve.seated){
                if(c1>ve.staff){
                    feas=false;
                    i=g-l+1;
                }
            }else{
                if(c2+c1>ve.staff+ve.seated){
                    feas=false;
                    i=g-l+1;
                }
            }
        }else{
            if(c2>ve.seated){
                if(c3+c1>ve.staff+ve.stretcher){
                    
                    feas=false;
                    i=g-l+1;
                }
            }else{
                if(c2+c3+c1>ve.staff+ve.stretcher+ve.seated){
                    feas=false;
                    i=g-l+1;
                    
                }
            }
            
        }
        if(c3>ve.stretcher){
            if(c2>ve.seated){
                feas=false;
                i=g-l+1;
            }
        }else{
            if(c2+c3>ve.stretcher+ve.seated){
                feas=false;
                i=g-l+1;
            }
        }
        if(c4>ve.wheelchair){
            
            feas=false;
            i=g-l+1;
            
        }
        
        if(c3>ve.stretcher){
            feas=false;
            i=g-l+1;
        }
        
        
    
        i=i+1;
        c10=c1;
        c20=c2;
        c30=c3;
        c40=c4;
        
        
        }else{
            i=g-l+1;
        }

    }
    
    if(feas==true){
        c1=c10-ric[req]->staff;
        c2=c20-ric[req]->seated;
        c3=c30-ric[req]->stretcher;
        c4=c40-ric[req]->wheelchair;
        if(c3>ve.stretcher){
            if(c2>ve.seated){
                if(c1>ve.staff){
                    feas=false;
                    i=g-l+1;
                }
            }else{
                if(c2+c1>ve.staff+ve.seated){
                    feas=false;
                    i=g-l+1;
                }
            }
        }else{
            if(c2>ve.seated){
                if(c3+c1>ve.staff+ve.stretcher){
                    
                    feas=false;
                    i=g-l+1;
                }
            }else{
                if(c2+c3+c1>ve.staff+ve.stretcher+ve.seated){
                    feas=false;
                    i=g-l+1;
                    
                }
            }
            
        }
        if(c3>ve.stretcher){
            if(c2>ve.seated){
                feas=false;
                i=g-l+1;
            }
        }else{
            if(c2+c3>ve.stretcher+ve.seated){
                feas=false;
                i=g-l+1;
            }
        }
        if(c4>ve.wheelchair){
            
            feas=false;
            i=g-l+2;
            
        }
        
        if(c3>ve.stretcher){
            feas=false;
            i=g-l+1;
        }


    
    }
    
    
    
    return feas;

}

double Route::effect_of_inserting_req_on_pos(int req, int l, int g, std::vector<std::vector<double>> &Dist, std::vector<RichiesteServizio*> &ric){

    double ef=0, ef1=0;
    

    for(int i=1; i<l+1; i++){
    
        ef+=Dist[locations[i-1]][locations[i]];
    }
    ef+=Dist[locations[l]][ric[req]->locationP];
    if(l==g){
        ef+=Dist[ric[req]->locationP][ric[req]->locationD];
    }else{
        ef+=Dist[ric[req]->locationP][locations[l+1]];
        for(int i=l+2; i<g+1; i++){
            ef+=Dist[locations[i-1]][locations[i]];
        
        }
        
        ef+=Dist[locations[g]][ric[req]->locationD];
    }
    
    ef+=Dist[ric[req]->locationD][locations[g+1]];
    for(int i=g+2; i<length; i++){
    
        ef+=Dist[locations[i-1]][locations[i]];
    }
    
    ef1=ef-totaldist;
                                 
    
    return ef1;
}


bool Route::check_feas_D_tw1(double* earl, int req, std::vector<RichiesteServizio*> &ric, int g, int l, std::vector<std::vector<double>> &Time, std::vector<Veicoli*> &veic){

    bool feas;

        //double* lat;
        double lat0=0;
        int p=length+2;
        //lat=new double[p];
        //for(int i=0; i<p; i++){
        //    lat[i]=0;
       // }
       // lat[p-1]=veic[veh]->MaxRoute;
       lat0=veic[veh]->MaxRoute;
        for(int i=p-2; i>g+2; i--){

            RichiesteServizio& rich=*ric[Ricid[i-2]];
            //lat[i]=lat[i+1]-Time[locations[i-2]][locations[i-1]]-rich.st;
            lat0=lat0-Time[locations[i-2]][locations[i-1]]-rich.st;
            //assert(lat[i]==lat0);
            if(type[i-2]==1){
                if(lat0>rich.timewinPmax){
                    lat0=rich.timewinPmax;
                }
            }else{if(type[i-2]==-1){
                if(lat0>rich.timewinDmax){
                    lat0=rich.timewinDmax;
                }
            }
            }

        }

        lat0=lat0-Time[ric[req]->locationD][locations[g+1]]-ric[req]->st;
        if(lat0>ric[req]->timewinDmax){
            lat0=ric[req]->timewinDmax;
        }

        //hasta aqui esta bien

        if(g==l){
        /*
            lat[g+1]=lat[g+2]-Time[ric[req]->locationP][ric[req]->locationD]-ric[req]->st;
            if(lat[g+1]>ric[req]->timewinPmax){
                lat[g+1]=ric[req]->timewinPmax;
            }*/

        }else{

            lat0=lat0-Time[locations[g]][ric[req]->locationD]-ric[Ricid[g]]->st;
            if(type[g]==1){
                if(lat0>ric[Ricid[g]]->timewinPmax){
                    lat0=ric[Ricid[g]]->timewinPmax;
                }
            }else{if(type[g]==-1){
                if(lat0>ric[Ricid[g]]->timewinDmax){
                    lat0=ric[Ricid[g]]->timewinDmax;
                }
            }
            }


            for(int i=g; i>l+1; i--){
                RichiesteServizio& rich=*ric[Ricid[i-1]];
                lat0=lat0-Time[locations[i-1]][locations[i]]-rich.st;
                if(type[i-1]==1){
                    if(lat0>rich.timewinPmax){
                        lat0=rich.timewinPmax;
                    }
                }else{if(type[i-1]==-1){
                    if(lat0>rich.timewinDmax){
                        lat0=rich.timewinDmax;
                    }
                }
                }


            }






        }



        //???

        //???

      /*  double earll;
        earll=earl[g-l+1]+ric[req]->st+Time[ric[req]->locationD][locations[g+1]];

        if(type[g+1]==1){
            if(earll<ric[Ricid[g+1]]->timewinPmin){
                earll=ric[Ricid[g+1]]->timewinPmin;
            }

        }else{

            if(type[g+1]==-1){
                if(earll<ric[Ricid[g+1]]->timewinDmin){
                    earll=ric[Ricid[g+1]]->timewinDmin;
                }
            }

        }*/
      //  std::cout << earll << " = " << earl[g-l+1]<< "+" <<ric[req].st<< "+" << Time[ric[req].locationD][locations[g+1]]<<std::endl;
        if(earl[g-l+1]<=ric[req]->timewinDmax && earl[1]<=lat0){
            feas=true;
        }else{
            feas=false;
        }
       // for(int i=0; i<p; i++){
       //     std::cout<< lat[i] << " " ;

       // }std::cout<< std::endl;

        //std::cout<< earl[g-l+1] << "<=" << ric[req].timewinDmax << " && " << earll<<"<=" <<lat[g+3] <<std::endl;

        //delete[] lat;
        return feas;
}


bool Route::EightStepEvaluation(std::vector<std::vector<double>> &Time, std::vector<RichiesteServizio*>&ric, int l, int g, int req, std::vector<Veicoli*>&veic){
    int i;
    bool feas=true;
    int len=length+2;

    /*double* arr, *dep, * wait;
    arr=new double[len];
    dep=new double[len];
    wait=new double[len];


    int* loc, *rid, *typ;
    loc=new int[len];
    rid=new int[len];
    typ=new int[len];*/

    double arr[len];
    double dep[len];
    double wait[len];
    int loc[len];
    int rid[len];
    int typ[len];

    for(i=0;i<l+1;i++){
        loc[i]=locations[i];
        rid[i]=Ricid[i];
        typ[i]=type[i];
    }
    loc[l+1]=ric[req]->locationP;
    rid[l+1]=req;
    typ[l+1]=1;
    for(i=l+1;i<g+1;i++){
        loc[i+1]=locations[i];
        rid[i+1]=Ricid[i];
        typ[i+1]=type[i];
    }
    loc[g+2]=ric[req]->locationD;
    rid[g+2]=req;
    typ[g+2]=-1;
    for(i=g+1;i<length;i++){
        loc[i+2]=locations[i];
        rid[i+2]=Ricid[i];
        typ[i+2]=type[i];
    }
    
    
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
                            
                            if(dep[jj]-ric[rid[jj]]->st-dep[ii]-ric[rid[ii]]->RideTime>0){
                                v3+=(dep[jj]-ric[rid[jj]]->st-dep[ii]-ric[rid[ii]]->RideTime);
                                
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
                        if(dep[jj]-ric[rid[jj]]->st-dep[ii]-ric[rid[ii]]->RideTime>0){
                            
                            
                            v3+=(dep[jj]-ric[rid[jj]]->st-dep[ii]-ric[rid[ii]]->RideTime);
                        }else{
                        }
                    }
                }
                
            }
        }
        
    }
    
    
    if(arr[len-1]-dep[0]-veic[veh]->MaxRoute>0){
        
        //V[i][1]=arr[length-1]-dep[0]-MaxTimeRoute;
        v2=arr[len-1]-dep[0]-veic[veh]->MaxRoute;
        
        
    }else{v2=0;}
    
    
    if(v1>0 || v2>0 || v3>0){
        feas=false;
    }
    /*
    delete[] loc;
    delete[] rid;
    delete[] typ;
    delete[] arr;
    delete[] dep;
    delete[] wait;
*/
    return feas;
}

double calcF(int* rid, double* wait, double* arr, double* dep, int a, std::vector<RichiesteServizio*> &ric, int* typ, int len){

    double F;
    double F1;
    int i, j, l;
    F=DBL_MAX;
    //F1=new double[len-a];
    double sum, M1, M2, timetravelled, M;
    
    for(i=a; i<len; i++){
        sum=0;
        
        for(j=a+1; j<i+1; j++){
            
            sum+=wait[j];
            
        }
        
        
        if(typ[i]==0){
            
            M1=1440-(dep[i]-ric[rid[i-1]]->st);
            
            M2=1440-0;
            
        }else{
            
            if(typ[i]==1){//es de pickup
                
                M1=ric[rid[i]]->timewinPmax-(dep[i]-ric[rid[i]]->st);
                M2=DBL_MAX;
            }else{
                M1=ric[rid[i]]->timewinDmax-(dep[i]-ric[rid[i]]->st);
                bool yes=false;
                timetravelled=0;
                for(l=i-1; l>0; l--){
                    
                    if(rid[i]==rid[l]){
                        timetravelled=arr[i]-dep[l];
                        if(l<a){yes=true;}
                    }
                }
                if(yes==true){
                    M2=ric[rid[i]]->RideTime-timetravelled;
                }else{M2=FLT_MAX;}
            }}
        if(M1<=0){M1=0;}
        if(M2<=0){M2=0;}
        
        if(M1<=M2){
            M=M1;
        }else{M=M2;}
        
        
        
        F1=sum+M;
        if(F>F1){
            
            F=F1;
            
        }
   
        
    }
    
    
    
   // delete[] F1;

    return F;
}


void Route::resize_from(int nl){
    
    //std::cout << y << std::endl;



       if(length<nl){
           //insert the difference

           int y1=nl-length;
           locations.insert(locations.end()-2,y1,0);
           Ricid.insert(Ricid.end()-2,y1,0);
           type.insert(type.end()-2,y1,0);
           arrival.insert(arrival.end()-2,y1,0);
           departure.insert(departure.end()-2,y1,0);
           waiting.insert(waiting.end()-2,y1,0);
           earliest.insert(earliest.end()-2,y1,0);
           latest.insert(latest.end()-2,y1,0);
           cap1.insert(cap1.end()-2,y1,0);
           maxcap1.insert(maxcap1.end()-2,y1,0);
           cap2.insert(cap2.end()-2,y1,0);
           maxcap2.insert(maxcap2.end()-2,y1,0);
           cap3.insert(cap3.end()-2,y1,0);
           maxcap3.insert(maxcap3.end()-2,y1,0);
           cap4.insert(cap4.end()-2,y1,0);
           maxcap4.insert(maxcap4.end()-2,y1,0);
           length=locations.size();

       }else{
           if(length>nl){
               //delete the difference
               int y=length-nl;
               locations.erase(locations.end()-1-(y),locations.end()-1);
               Ricid.erase(Ricid.end()-1-(y),Ricid.end()-1);
               type.erase(type.end()-1-y,type.end()-1);
               arrival.erase(arrival.end()-1-y,arrival.end()-1);
               departure.erase(departure.end()-1-y,departure.end()-1);
               waiting.erase(waiting.end()-1-y,waiting.end()-1);
               earliest.erase(earliest.end()-1-y,earliest.end()-1);
               latest.erase(latest.end()-1-y,latest.end()-1);
               cap1.erase(cap1.end()-1-y,cap1.end()-1);
               maxcap1.erase(maxcap1.end()-1-y,maxcap1.end()-1);
               cap2.erase(cap2.end()-1-y,cap2.end()-1);
               maxcap2.erase(maxcap2.end()-1-y,maxcap2.end()-1);
               cap3.erase(cap3.end()-1-y,cap3.end()-1);
               maxcap3.erase(maxcap3.end()-1-y,maxcap3.end()-1);
               cap4.erase(cap4.end()-1-y,cap4.end()-1);
               maxcap4.erase(maxcap4.end()-1-y,maxcap4.end()-1);
               length=locations.size();



           }
       }



}



bool Route::EightStepEvaluation_2(std::vector<std::vector<double>> &Time, std::vector<RichiesteServizio*>&ric, int u, int ii, int j, int z, std::vector<Veicoli*>&veic){
    int i;
    bool feas=true;
    int len=length;
    double* arr, *dep, * wait;
    arr=new double[len];
    dep=new double[len];
    wait=new double[len];
    int* loc, *rid, *typ;
    loc=new int[len];
    rid=new int[len];
    typ=new int[len];
    int req;
    for(i=0;i<u;i++){
        loc[i]=locations[i];
        rid[i]=Ricid[i];
        typ[i]=type[i];
    }
    int ii_0;
    for(i=0; i<3; i++){
        if(i==ii){
        
            ii_0=u;
        }else{if(j==i){
        
            ii_0=u+1;
        }else{
            ii_0=u+2;
        
        }}
    
        loc[u+i]=locations[ii_0];
        rid[u+i]=Ricid[ii_0];
        typ[u+i]=type[ii_0];
    }
    
    
    for(i=u+3;i<length; i++){
        loc[i]=locations[i];
        rid[i]=Ricid[i];
        typ[i]=type[i];
    
    }
    //std::cout << "u -> " << u << "i -> " << ii  << "j -> " << j  << "z -> " << z << "\n";

    //for(i=0;i<len;i++){
    //    std::cout<<loc[i] << "\t" << rid[i] << "\t" <<  typ[i] << "\n";
    
   // }
    
    
    
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
                            
                            if(dep[jj]-ric[rid[jj]]->st-dep[ii]-ric[rid[ii]]->RideTime>0){
                                v3+=(dep[jj]-ric[rid[jj]]->st-dep[ii]-ric[rid[ii]]->RideTime);
                                
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
                        if(dep[jj]-ric[rid[jj]]->st-dep[ii]-ric[rid[ii]]->RideTime>0){
                            
                            
                            v3+=(dep[jj]-ric[rid[jj]]->st-dep[ii]-ric[rid[ii]]->RideTime);
                        }else{
                        }
                    }
                }
                
            }
        }
        
    }
    
    
    if(arr[len-1]-dep[0]-veic[veh]->MaxRoute>0){
        
        //V[i][1]=arr[length-1]-dep[0]-MaxTimeRoute;
        v2=arr[len-1]-dep[0]-veic[veh]->MaxRoute;
        
        
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



bool Route::check_feas_FirstRequest(int p2, int d2, std::vector<std::vector<double>> &Time, std::vector<RichiesteServizio*> &ric, int req){
    //Checks wheather all vertices between p2 and d2 can be reached in time (earliest(d2)<=TWdeliveryreq && earliest(p2+1)<=latest(p2+1)
    bool f;
    double earl=earliest[p2-1];
    double earl0, lat, lat0;
    int loc0=locations[p2-1];
    int rr0;


    //P2 earliest calculation
    if(p2-1==0){
        earl=earl+Time[loc0][ric[req]->locationP];


    }else{
        earl=earl+ric[req]->st+Time[loc0][ric[req]->locationP];

    }
    if(earl<ric[req]->timewinPmin){
        earl=ric[req]->timewinPmin;
    }

    loc0=ric[req]->locationP;
    rr0=req;


    for(int i=p2+1; i<d2; i++){
        RichiesteServizio& rich =*ric[Ricid[i]];

        earl=earl+ric[rr0]->st+Time[loc0][locations[i]];

        if(type[i]==1){
            if(earl<rich.timewinPmin){
                earl=rich.timewinPmin;
            }
        }else{
            if(type[i]==-1){
                if(earl<rich.timewinDmin){
                    earl=rich.timewinDmin;
                }

            }

        }
        if(i==p2+1){earl0=earl;}
        rr0=Ricid[i];
        loc0=locations[i];

    }

    earl=earl+ric[rr0]->st+Time[loc0][ric[req]->locationD];

    if(earl<ric[req]->timewinDmin){
        earl=ric[req]->timewinDmin;
    }


    loc0=ric[req]->locationD;
    rr0=req;

    if(d2==p2+1){
        earl0=earl;
    }


    lat0=latest[d2+1];

    lat=lat0-Time[locations[d2+1]][ric[req]->locationD]-ric[req]->st;

    if(lat>ric[req]->timewinDmax){
        lat=ric[req]->timewinDmax;
    }

    lat0=lat;
    loc0=ric[req]->locationD;
    rr0=req;

    for(int i=d2-1; i<p2;i++){

        lat=lat0-Time[loc0][locations[i]]-ric[Ricid[i]]->st;
        if(type[i]==1){
            if(lat>ric[req]->timewinPmax){
                lat=ric[req]->timewinPmax;
            }

        }else{
            if(type[i]==-1){
                if(lat>ric[req]->timewinDmax){
                    lat=ric[req]->timewinDmax;
                }
            }

        }

        lat0=lat;
        loc0=locations[i];
        rr0=Ricid[i];

    }


    if(earl<ric[req]->timewinDmax && earl0<lat0){
        f=true;
    }else{

        f=false;
    }




    return f;
}


bool Route::EightStep_2(int p2, int d2, int r, std::vector<std::vector<double>> &Time, std::vector<RichiesteServizio*>&ric, std::vector<Veicoli*>&veic){

    int i;
    bool feas=true;
    int len=length;
    double* arr, *dep, * wait;
    arr=new double[len];
    dep=new double[len];
    wait=new double[len];
    int* loc, *rid, *typ;
    loc=new int[len];
    rid=new int[len];
    typ=new int[len];
    int req;
    for(i=0;i<len;i++){
        loc[i]=locations[i];
        rid[i]=Ricid[i];
        typ[i]=type[i];
    }
    loc[p2]=ric[r]->locationP;
    loc[d2]=ric[r]->locationD;
    typ[p2]=1;
    typ[d2]=-1;
    rid[p2]=r;
    rid[d2]=r;

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

                            if(dep[jj]-ric[rid[jj]]->st-dep[ii]-ric[rid[ii]]->RideTime>0){
                                v3+=(dep[jj]-ric[rid[jj]]->st-dep[ii]-ric[rid[ii]]->RideTime);

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
                        if(dep[jj]-ric[rid[jj]]->st-dep[ii]-ric[rid[ii]]->RideTime>0){


                            v3+=(dep[jj]-ric[rid[jj]]->st-dep[ii]-ric[rid[ii]]->RideTime);
                        }else{
                        }
                    }
                }

            }
        }

    }


    if(arr[len-1]-dep[0]-veic[veh]->MaxRoute>0){

        //V[i][1]=arr[length-1]-dep[0]-MaxTimeRoute;
        v2=arr[len-1]-dep[0]-veic[veh]->MaxRoute;


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

bool Route::capacity_P_feasibility_3( int l, int req, std::vector<Veicoli*> &veic, std::vector<RichiesteServizio*>  &ric, int p2, int d2, int req2){


    bool feas=true;
    int c1=0, c2=0, c3=0, c4=0;

    if(p2>=l){

        c1=cap1[l-1]+ric[req]->staff;
        c2=cap1[l-1]+ric[req]->seated;
        c3=cap1[l-1]+ric[req]->stretcher;
        c4=cap1[l-1]+ric[req]->wheelchair;


    }else{
        if(d2>l){
            c1=cap1[l+1]+ric[req]->staff;
            c2=cap1[l+1]+ric[req]->seated;
            c3=cap1[l+1]+ric[req]->stretcher;
            c4=cap1[l+1]+ric[req]->wheelchair;


        }else{
        if(cap1[l]>0){
            c1=cap1[l]-ric[req2]->staff+ric[req]->staff;
        }else{
            c1=cap1[l]+ric[req]->staff;
        }
        if(cap2[l]>0){
            c2=cap2[l]-ric[req2]->seated+ric[req]->seated;
        }else{
            c2=cap2[l]+ric[req]->seated;
        }
        if(cap3[l]>0){
            c3=cap3[l]-ric[req2]->stretcher+ric[req]->stretcher;
        }else{
            c3=cap3[l]+ric[req]->stretcher;
        }
        if(cap4[l]>0){
            c4=cap4[l]-ric[req2]->wheelchair+ric[req]->wheelchair;
        }else{
            c4=cap4[l]+ric[req]->wheelchair;
        }

    }
    }

    if(c3>veic[veh]->stretcher){
        if(c2>veic[veh]->seated){
            if(c1>veic[veh]->staff){
                feas=false;

            }else{

            }
        }else{
            if(c2+c1>veic[veh]->staff+veic[veh]->seated){
                feas=false;
            }else{


            }

        }
    }else{
        if(c2>veic[veh]->seated){
            if(c3+c1>veic[veh]->staff+veic[veh]->stretcher){

                feas=false;
            }else{}
        }else{
            if(c2+c3+c1>veic[veh]->staff+veic[veh]->stretcher+veic[veh]->seated){
                feas=false;

            }else{

            }



        }

    }



    if(c3>veic[veh]->stretcher){
        if(c2>veic[veh]->seated){
            feas=false;
        }else{

        }
    }else{
        if(c2+c3>veic[veh]->stretcher+veic[veh]->seated){
            feas=false;
        }else{

        }

    }




    if(c4>veic[veh]->wheelchair){

        feas=false;

    }else{}

    if(c3>veic[veh]->stretcher){
        feas=false;
    }else{

    }

    return feas;
}


bool Route::tw_P_feasibility_3(std::vector<RichiesteServizio*> &ric, std::vector<std::vector<double>>&Time, int l, int req, int p2, int d2, int* loc, int*typ, int*rq){
    int j=0;
    int i=0;
    while(i<length){
        if(i==p2 || i==d2){
            i++;
        }else{
            if(j==l){
                j++;
            }else{
                loc[j]=locations[i];
                rq[j]=Ricid[i];
                typ[j]=type[i];
                j++;
                i++;
            }
        }
    }
    rq[l]=req;
    loc[l]=ric[req]->locationP;
    typ[l]=1;
   // std::cout << "ST" << std::endl;
    //for(i=0; i<length; i++){
   //     std::cout<< rq[i] << " " << typ[i] << " " << loc[j] << std::endl;


   // }

    bool feas;

    RichiesteServizio& rich=*ric[req];
    double  earl;

    earl=Time[locations[0]][loc[1]];
    if(earl<ric[rq[1]]->timewinPmin){
        earl=ric[rq[1]]->timewinPmin;

    }

    for(int i=2; i<l+1; i++){
        earl=earl+ric[rq[i-1]]->st+Time[locations[i-1]][locations[i]];
        if(typ[i]==1){
            if(earl<ric[rq[1]]->timewinPmin){
                earl=ric[rq[1]]->timewinPmin;

            }

        }else{
            if(typ[i]==-1){
                if(earl<ric[rq[1]]->timewinDmin){
                    earl=ric[rq[1]]->timewinDmin;

                }

            }

        }



    }





    if(earl<=rich.timewinPmax){
        feas=true;
    }else{
        feas=false;
    }

    return feas;
}

//  earl=Sol->route[vei]->calc_earl3(F[k], ric, p1, d1, l , g);
void Route::calc_rid3(int req, int p, int d, int l, int g, int* rid ){

    int j=0;
    int i=0;
    while(i<length){
        if(i==p || i==d){
            i++;
        }else{
            if(j==l ||  j==g){
                j++;
            }else{
                rid[j]=Ricid[i];
                j++;
                i++;
            }
        }
    }
    rid[l]=req;
    rid[g]=req;



}
void Route::calc_typ3(int req, int p, int d, int l, int g, int* typ ){

    int j=0;
    int i=0;
    while(i<length){
        if(i==p || i==d){
            i++;
        }else{
            if(j==l ||  j==g){
                j++;
            }else{
                typ[j]=type[i];
                j++;
                i++;
            }
        }
    }
    typ[l]=1;
    typ[g]=-1;




}

void Route::calc_earl3(std::vector<RichiesteServizio*> &ric, int*typ, int*loc, int*rid, double*earl,std::vector<std::vector<double>> &Time ){

    earl[1]=Time[loc[0]][loc[1]];
    if(earl[1]<ric[rid[1]]->timewinPmin){
        earl[1]=ric[rid[1]]->timewinPmin;
    }


    for(int i=2; i<length; i++){
        earl[i]=earl[i-1]+ric[rid[i-1]]->st+Time[loc[i-1]][loc[i]];
        if(typ[i]==1){
            if(earl[i]<ric[rid[i]]->timewinDmin){
                earl[i]=ric[rid[i]]->timewinDmin;
            }
        }else{
            if(typ[i]==-1){
                if(earl[i]<ric[rid[i]]->timewinDmin){
                    earl[i]=ric[rid[i]]->timewinDmin;
                }
            }

        }


    }



}
void Route::calc_loc3(int req, std::vector<RichiesteServizio*> &ric, int p, int d, int l, int g, int* loc ){

    int j=0;
    int i=0;
    while(i<length){
        if(i==p || i==d){
            i++;
        }else{
            if(j==l ||  j==g){
                j++;
            }else{
                loc[j]=locations[i];
                j++;
                i++;
            }
        }
    }
    loc[l]=ric[req]->locationP;
    loc[g]=ric[req]->locationD;


}
double Route::effect_of_inserting_req_on_pos_3(int req, int l, int g, std::vector<std::vector<double>> &Dist, std::vector<RichiesteServizio*> &ric, int p2, int d2, int* loc){

    double ef=0, ef1=0;





    for(int i=1; i<length; i++){

        ef+=Dist[loc[i-1]][loc[i]];

    }
    ef1=ef-totaldist;

    return ef1;
}

bool Route::check_feas_D_tw1_3(int req,std::vector<RichiesteServizio*> &ric,int l,int g,std::vector<std::vector<double>> &Time,std::vector<Veicoli*> &veic, int* loc, int* rid, int* typ){
    //Checks wheather all vertices between l and g can be reached in time (earliest(g)<=TWdeliveryreq && earliest(l+1)<=latest(l+1)
    bool f;

    double earl=0;
    double earl0=0;
    double lat;
    earl=Time[loc[0]][loc[1]];
    if(earl<ric[rid[1]]->timewinPmin){
        earl=ric[rid[1]]->timewinPmin;
    }

    for(int i=2; i<g+1; i++){
        earl=earl+ric[rid[i-1]]->st+Time[loc[i-1]][loc[i]];
        if(typ[i]==1){
            if(earl<ric[rid[i]]->timewinDmin){
                earl=ric[rid[i]]->timewinDmin;
            }
        }else{
            if(typ[i]==-1){
                if(earl<ric[rid[i]]->timewinDmin){
                    earl=ric[rid[i]]->timewinDmin;
                }
            }

        }
        if(i==l+1){
            earl0=earl;
        }

    }

    lat=veic[veh]->MaxRoute;



    for(int i=length-2; i>l; i--){
        RichiesteServizio& rich =*ric[Ricid[i]];

        lat=lat-Time[loc[i]][loc[i+1]]-rich.st;

        if(typ[i]==1){
            if(lat>rich.timewinPmax){
                lat=rich.timewinPmax;
            }


        }else{
            if(lat>rich.timewinDmax){
                lat=rich.timewinDmax;
            }

        }
    }

    if(earl<ric[rid[g]]->timewinDmax && earl0<lat){
        f=true;
    }else{
        f=false;
    }



    return f;


}
bool Route::check_feas_D_tw22(double* earl, int* typ, int *rid, int l, int g, std::vector<RichiesteServizio*> &ric){
    int i=l+1;
    bool feas=true;
    int a=g;
    while(i<a){

        if(typ[i]==1){
            if(earl[i]<=ric[rid[i]]->timewinPmax){
                feas=true;

            }else{
                feas=false;
                i=a;
            }


        }else{
            if(earl[i]<=ric[rid[i]]->timewinDmax){
                feas=true;

            }else{
                feas=false;
                i=a;
            }




        }
        i=i+1;
    }

    if(earl[l+1]>ric[rid[l+1]]->timewinDmax){
        feas=false;
    }


    return feas;




}


bool Route::check_feas_D_tw221(double* earl, int l, int g, std::vector<RichiesteServizio*> &ric){
    int i=l+1;
    bool feas=true;
    int a=g;
    while(i<a){

        if(type[i]==1){
            if(earl[i]<=ric[Ricid[i]]->timewinPmax){
                feas=true;

            }else{
                feas=false;
                i=a;
            }


        }else{
            if(earl[i]<=ric[Ricid[i]]->timewinDmax){
                feas=true;

            }else{
                feas=false;
                i=a;
            }




        }
        i=i+1;
    }

    if(earl[l+1]>ric[Ricid[l+1]]->timewinDmax){
        feas=false;
    }


    return feas;




}


bool Route::ridetime_feas_D22(int g, std::vector<RichiesteServizio*> &ric, int req, std::vector<std::vector<double>> &Time, double* earl){


    bool feas;
    //puedo calcular desde earl l
    feas=true;
    if(earl[g]-ric[req]->timewinPmax-ric[req]->st>ric[req]->RideTime){
        feas=false;
    }

    // std::cout<<" earl  " ;
    //  for(int i=0;i<g-l+2; i++){
    //     std::cout<<earl[i]<<" ";

    //  }std::cout<<std::endl;

    return feas;
}

void Route::calc_earl1(int p2, int d2, int req, double * earl, std::vector<std::vector<double>> &Time, std::vector<RichiesteServizio*> &ric){

    int loc0, rq0;
    if(p2==1){

        earl[1]=Time[locations[0]][ric[req]->locationP];
        if(earl[1]<ric[req]->timewinPmin){

            earl[1]=ric[req]->timewinPmin;
        }


        loc0=ric[req]->locationP;
        rq0=req;
    }else{
        earl[1]=Time[locations[0]][locations[1]];
        if(earl[1]<ric[Ricid[1]]->timewinPmin){

            earl[1]=ric[Ricid[1]]->timewinPmin;
        }

        loc0=locations[1];
        rq0=Ricid[1];
    }

   // std::cout<< p2 << " " << d2 <<  " " << length << std::endl;

    for(int i=2; i<p2; i++){
        //std::cout << i << std::endl;
        earl[i]=ric[rq0]->st+Time[loc0][locations[i]];
        if(type[i]==1){
            if(earl[i]<ric[Ricid[i]]->timewinPmin){

                earl[i]=ric[Ricid[i]]->timewinPmin;
            }


        }else{
            if(earl[i]<ric[Ricid[i]]->timewinDmin){

                earl[i]=ric[Ricid[i]]->timewinDmin;
            }



        }


        loc0=locations[i];
        rq0=Ricid[i];
    }
    if(p2>1){
    earl[p2]=ric[req]->st+Time[loc0][ric[req]->locationP];

    if(earl[p2]<ric[req]->timewinPmin){

        earl[p2]=ric[req]->timewinPmin;
    }

        loc0=ric[req]->locationP;
        rq0=req;
    }


    for(int i=p2+1; i<d2; i++){
        earl[i]=ric[rq0]->st+Time[loc0][locations[i]];
        if(type[i]==1){
            if(earl[i]<ric[Ricid[i]]->timewinPmin){

                earl[i]=ric[Ricid[i]]->timewinPmin;
            }


        }else{
            if(earl[i]<ric[Ricid[i]]->timewinDmin){

                earl[i]=ric[Ricid[i]]->timewinDmin;
            }



        }
        loc0=locations[i];
        rq0=Ricid[i];

    }
    earl[d2]=ric[rq0]->st+Time[loc0][ric[req]->locationD];



}
void Route::delete_req_2(int m,  std::vector<std::vector<double>>  &Time, int a, int b){


    // std::cout << " cheking " << std::endl;
    //std::cout << a << " " << b << std::endl;
    remove_pos(b);
    remove_pos(a);
    length=locations.size();
    numRicRoute=(length-2)/2;
    totaldist=calculatedist(Time);
}




bool Route::check_cap_from23(int l, int g,std::vector<Veicoli*> &veic, int req, std::vector<RichiesteServizio*> &ric){
    bool feas=true;

    int i;

    Veicoli& ve=*veic[veh];


    int c10, c20, c30, c40;
    int c1, c2, c3, c4;
    c10=cap1[l-1]+ric[req]->staff;
    c20=cap2[l-1]+ric[req]->seated;
    c30=cap3[l-1]+ric[req]->stretcher;
    c40=cap4[l-1]+ric[req]->wheelchair;
    i=l+1;
    //std::cout<< g << l << std::endl;

    while(i<g){
        if(feas==true){
            int r=Ricid[i];
            //std::cout << r << std::endl;

            if(type[i]==1){
                //std::cout <<type[i+l]<< std::endl;
                c1=c10+ric[r]->staff;
                c2=c20+ric[r]->seated;
                c3=c30+ric[r]->stretcher;
                c4=c40+ric[r]->wheelchair;
            }else{
                c1=c10-ric[r]->staff;
                c2=c20-ric[r]->seated;
                c3=c30-ric[r]->stretcher;
                c4=c40-ric[r]->wheelchair;

            }

            if(c3>ve.stretcher){
                if(c2>ve.seated){
                    if(c1>ve.staff){
                        feas=false;
                        i=g;
                    }
                }else{
                    if(c2+c1>ve.staff+ve.seated){
                        feas=false;
                        i=g;
                    }
                }
            }else{
                if(c2>ve.seated){
                    if(c3+c1>ve.staff+ve.stretcher){

                        feas=false;
                        i=g;
                    }
                }else{
                    if(c2+c3+c1>ve.staff+ve.stretcher+ve.seated){
                        feas=false;
                        i=g;

                    }
                }

            }
            if(c3>ve.stretcher){
                if(c2>ve.seated){
                    feas=false;
                    i=g;
                }
            }else{
                if(c2+c3>ve.stretcher+ve.seated){
                    feas=false;
                    i=g;
                }
            }
            if(c4>ve.wheelchair){

                feas=false;
                i=g;

            }

            if(c3>ve.stretcher){
                feas=false;
                i=g;
            }



            i=i+1;
            c10=c1;
            c20=c2;
            c30=c3;
            c40=c4;


        }else{
            i=g;
        }

    }

    if(feas==true){
        c1=c10-ric[req]->staff;
        c2=c20-ric[req]->seated;
        c3=c30-ric[req]->stretcher;
        c4=c40-ric[req]->wheelchair;
        if(c3>ve.stretcher){
            if(c2>ve.seated){
                if(c1>ve.staff){
                    feas=false;

                }
            }else{
                if(c2+c1>ve.staff+ve.seated){
                    feas=false;

                }
            }
        }else{
            if(c2>ve.seated){
                if(c3+c1>ve.staff+ve.stretcher){

                    feas=false;

                }
            }else{
                if(c2+c3+c1>ve.staff+ve.stretcher+ve.seated){
                    feas=false;


                }
            }

        }
        if(c3>ve.stretcher){
            if(c2>ve.seated){
                feas=false;

            }
        }else{
            if(c2+c3>ve.stretcher+ve.seated){
                feas=false;

            }
        }
        if(c4>ve.wheelchair){

            feas=false;


        }

        if(c3>ve.stretcher){
            feas=false;

        }

    }



    return feas;

}

bool Route::check_cap_from22(int l, int g, std::vector<Veicoli*> &veic, int req, std::vector<RichiesteServizio*> &ric, int* rid, int*typ ){

    int i;

    Veicoli& ve=*veic[veh];

    bool feas=true;
    int c10, c20, c30, c40;
    int c1, c2, c3, c4;
    c10=0;
    c20=0;
    c30=0;
    c40=0;

    i=1;
    //std::cout<< g << l << std::endl;

    while(i<g){
        if(feas==true){
            int r=rid[i];
            //std::cout << r << std::endl;

            if(typ[i]==1){
                //std::cout <<type[i+l]<< std::endl;
                c1=c10+ric[r]->staff;
                c2=c20+ric[r]->seated;
                c3=c30+ric[r]->stretcher;
                c4=c40+ric[r]->wheelchair;
            }else{
                c1=c10-ric[r]->staff;
                c2=c20-ric[r]->seated;
                c3=c30-ric[r]->stretcher;
                c4=c40-ric[r]->wheelchair;

            }

            if(c3>ve.stretcher){
                if(c2>ve.seated){
                    if(c1>ve.staff){
                        feas=false;
                        i=g;
                    }
                }else{
                    if(c2+c1>ve.staff+ve.seated){
                        feas=false;
                        i=g;
                    }
                }
            }else{
                if(c2>ve.seated){
                    if(c3+c1>ve.staff+ve.stretcher){

                        feas=false;
                        i=g;
                    }
                }else{
                    if(c2+c3+c1>ve.staff+ve.stretcher+ve.seated){
                        feas=false;
                        i=g;

                    }
                }

            }
            if(c3>ve.stretcher){
                if(c2>ve.seated){
                    feas=false;
                    i=g;
                }
            }else{
                if(c2+c3>ve.stretcher+ve.seated){
                    feas=false;
                    i=g;
                }
            }
            if(c4>ve.wheelchair){

                feas=false;
                i=g;

            }

            if(c3>ve.stretcher){
                feas=false;
                i=g;
            }



            i=i+1;
            c10=c1;
            c20=c2;
            c30=c3;
            c40=c4;


        }else{
            i=g;
        }

    }

    if(feas==true){
        c1=c10-ric[req]->staff;
        c2=c20-ric[req]->seated;
        c3=c30-ric[req]->stretcher;
        c4=c40-ric[req]->wheelchair;
        if(c3>ve.stretcher){
            if(c2>ve.seated){
                if(c1>ve.staff){
                    feas=false;

                }
            }else{
                if(c2+c1>ve.staff+ve.seated){
                    feas=false;

                }
            }
        }else{
            if(c2>ve.seated){
                if(c3+c1>ve.staff+ve.stretcher){

                    feas=false;
                }
            }else{
                if(c2+c3+c1>ve.staff+ve.stretcher+ve.seated){
                    feas=false;


                }
            }

        }
        if(c3>ve.stretcher){
            if(c2>ve.seated){
                feas=false;

            }
        }else{
            if(c2+c3>ve.stretcher+ve.seated){
                feas=false;

            }
        }
        if(c4>ve.wheelchair){

            feas=false;


        }

        if(c3>ve.stretcher){
            feas=false;

        }

    }


    return feas;
}
bool Route::calc_feas_rtime(int l,int g, std::vector<RichiesteServizio*> &ric,int req ,std::vector<std::vector<double>> &Time, double* earl){
    bool f;

    f=true;
    if(earl[g]-ric[req]->timewinPmax-ric[req]->st>ric[req]->RideTime){
        f=false;
    }

    return f;

}

bool Route::EightStep_3(std::vector<std::vector<double>> &Time, std::vector<RichiesteServizio*>&ric, std::vector<Veicoli*>&veic, int*loc, int*rid, int*typ){

    int i;
    bool feas=true;
    int len;
    len=length;


    double *arr, *dep, * wait;


    arr=new double[length];

    dep=new double[length];
    wait=new double[length];

    int req;


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

                            if(dep[jj]-ric[rid[jj]]->st-dep[ii]-ric[rid[ii]]->RideTime>0){
                                v3+=(dep[jj]-ric[rid[jj]]->st-dep[ii]-ric[rid[ii]]->RideTime);

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
                        if(dep[jj]-ric[rid[jj]]->st-dep[ii]-ric[rid[ii]]->RideTime>0){


                            v3+=(dep[jj]-ric[rid[jj]]->st-dep[ii]-ric[rid[ii]]->RideTime);
                        }else{
                        }
                    }
                }

            }
        }

    }


    if(arr[len-1]-dep[0]-veic[veh]->MaxRoute>0){

        //V[i][1]=arr[length-1]-dep[0]-MaxTimeRoute;
        v2=arr[len-1]-dep[0]-veic[veh]->MaxRoute;


    }else{v2=0;}


    if(v1>0 || v2>0 || v3>0){
        feas=false;
    }


    delete[] arr;
    delete[] dep;
    delete[] wait;

    return feas;

}

bool feasible_with_vehicle(int m, int* E, int vei, std::vector<std::vector<int>> &MatCompVei){

    bool feas=true;
    int i=0;
    while(i<m){
        if(MatCompVei[E[i]][vei]<1){
            feas=false;
            i=m;
        }
        i++;
    }



    return feas;
}

