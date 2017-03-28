//
//  Route.hpp
//  OrCode
//
//  Created by garazi on 25/01/2017.
//  Copyright © 2017 Garazi. All rights reserved.
//

#ifndef Route_hpp
#define Route_hpp

#include <stdio.h>
#include <vector>
#include"Veicoli.hpp"

#include"RichiesteServizio.hpp"
class Route{
    
public:
    
    int id;
    int veh;
    int length;
    int numRicRoute;
    double totaldist;
    
    
    std::vector<int> locations;
    std::vector<int> Ricid;
    std::vector<int> type;
    std::vector<double> arrival;
    std::vector<double> departure;
    std::vector<double> waiting;
    std::vector<double> earliest;
    std::vector<double> latest;
    std::vector <int> cap1;
    std::vector<int> maxcap1;
    std::vector <int> cap2;
    std::vector<int> maxcap2;
    std::vector <int> cap3;
    std::vector<int> maxcap3;
    std::vector <int> cap4;
    std::vector<int> maxcap4;
    
    Route();
    
    
    ~Route();
public:
    void resize_from(int nl);
    void InitializeSol(std::vector<Veicoli*> &veic, int numVeicoli, std::vector<std::vector<double>> &Time, int i);
    void clear_route();
    void insert_pickup(int req, int l, RichiesteServizio* ric);
    void insert_delivery(int req, int g, RichiesteServizio* ric);
    void delete_delivery(int g);
    //void delete_pickup(int g, double** Time);
    void insert_req_pos(int req, int l, int g, std::vector<RichiesteServizio*> &ric, std::vector<std::vector<double>> &Time);
    //void insert_req(int req, int l, int g, RichiesteServizio* ric, double** Time);
    void calculate_times(RichiesteServizio* ric, double** Time);
    //Check feasibility of inserting pickup vertex i, in position pre+1. Time Window
    bool check_feasibility_P_tw(int req, int pre, RichiesteServizio* ric, double** Time);
    //Check feasibility of inserting pickup vertex i,  in position pre+1. Capacity
    bool check_feasibility_P_cap(int req, int pre, Veicoli* veic, RichiesteServizio* ric);
    
    double effect_ondistance(int req, int p, int d, double** Dist, RichiesteServizio* ric);
    void display_route();
    bool check_feasibillity_D_tw1(int req,int pre, RichiesteServizio* ric, double** MatTemp);
    bool check_feasibility_D_tw2(int pre, int suc, RichiesteServizio* ric, double** Time);
    
    void calculate_earliest_latest(std::vector<RichiesteServizio*> &ric,  std::vector<std::vector<double>> &Time, std::vector<Veicoli*> &veic);
    void calculate_capacity(std::vector<RichiesteServizio*> &ric, std::vector<Veicoli*> &veic);
    bool check_feasibility_capacity(int pre,int suc, Veicoli* veic);
    bool ridetime_check(int suc, RichiesteServizio* ric, int MaxRideTime);
    bool Eightstepevaluationscheme(std::vector<std::vector<double>> &Time, std::vector<RichiesteServizio*> &ric, std::vector<Veicoli*> &veic);
    double calculateF(int a, std::vector<RichiesteServizio*> ric);
    bool check_feasibility_P(int req,int pre,RichiesteServizio* ric, double** Time, Veicoli* veic);
    double calculatedist(std::vector<std::vector<double>> &Time);
    
    void Initializeaddroute(int a, std::vector<std::vector<double>> &Time,  std::vector<Veicoli*>  &veic, int numVeicoli);
    void CopyRoute(Route* route);
    
    void remove_pos(int p);
    void rewriteid(int a);
    
    int* count_request(int* E);
    int find_pickup(int m);
    int find_delivery(int m);
    void delete_req(int m,  std::vector<std::vector<double>> &Time);
    int count_req3(int i);
    void change_order(int u, int i , int j, int z);
    bool delivery_before_pickup(int u);
    bool check_feas1(int u, std::vector<RichiesteServizio*> &ric, std::vector<std::vector<double>> & Time, int i, int j, int z);
    bool check_cap(int u, std::vector<Veicoli*> &veic, int i, int j, int z, std::vector<RichiesteServizio*> &ric);
    void change_vehicle(int a, Veicoli* veic);
    double new_dist_without_req(int req,std::vector<std::vector<double>> & MatTemp);
    bool capacity_P_feasibility(int l, int req, std::vector<Veicoli*> &veic, std::vector<RichiesteServizio*>  &ric);
    bool tw_P_feasibility(std::vector<RichiesteServizio*> &ric, std::vector<std::vector<double>> &Time, int l, int req);
    
    int* calcrid(int* rid,int l, int g, int req );
    int* calctyp(int* typ,int l, int g);
    int* calcloc(int * loc, int req, std::vector<RichiesteServizio*> &ric, int l, int g);
    double* calcearl(double* earl, int req, std::vector<RichiesteServizio*> &ric, std::vector<std::vector<double>> &Time, int l, int g, int* loc, int* rid, int* typ );
    bool ridetime_feas_D(int g, int l,std::vector<RichiesteServizio*> &ric, int req, std::vector<std::vector<double>> &Time, double* earl);
    
    bool check_cap_from(int l, int g, std::vector<Veicoli*> &veic, int req,std::vector<RichiesteServizio*> &ric);
    double effect_of_inserting_req_on_pos(int req, int l, int g, std::vector<std::vector<double>> &Dist, std::vector<RichiesteServizio*> &ric);
    bool check_feas_D_tw1(double* earl, int req, std::vector<RichiesteServizio*>&ric, int g, int l, std::vector<std::vector<double>>  &Time, std::vector<Veicoli*> &veic);
    bool EightStepEvaluation(std::vector<std::vector<double>> &Time, std::vector<RichiesteServizio*> &ric, int l, int g, int req, std::vector<Veicoli*>&veic);
    
    bool delivery_before_pickup_2(int u, int i, int j, int z);
    double calculatedist_2(std::vector<std::vector<double>> &Time, int i, int j, int z, int u);
    bool EightStepEvaluation_2(std::vector<std::vector<double>> &Time, std::vector<RichiesteServizio*>&ric, int u, int ii, int j, int z, std::vector<Veicoli*>&veic);
    
    bool capacity_feasibility2(int p1,int req, std::vector<Veicoli*> &veic, std::vector<RichiesteServizio*>  &ric);
    bool tw_P_feasibility_2(std::vector<RichiesteServizio*> &ric, std::vector<std::vector<double>>&Time, int p1, int req);
    bool check_feas_D_tw1_2(double* earl, int req, std::vector<RichiesteServizio*> &ric, int g, int l, std::vector<std::vector<double>> &Time, std::vector<Veicoli*> &veic);
};
bool check_feas_D_tw2(double* earl, int* typ, int* rid, int l, int g,std::vector<RichiesteServizio*> &ric);
double calcF(int* rid, double* wait, double* earl, double* dep, int a, std::vector<RichiesteServizio*> &ric, int* typ, int len);


#endif /* Route_hpp */
