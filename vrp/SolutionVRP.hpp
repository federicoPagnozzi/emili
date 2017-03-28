//
//  Solution.hpp
//  OrCode
//
//  Created by garazi on 25/01/2017.
//  Copyright © 2017 Garazi. All rights reserved.
//

#ifndef SolutionVRP_hpp
#define SolutionVRP_hpp

#include <stdio.h>

#include "../emilibase.h"
#include "Veicoli.hpp"
#include <vector>
#include "Route.hpp"


class SolutionVRP : public emili::Solution{
protected:
    virtual const void* getRawData() const;


    virtual void setRawData(const void* data);

    
public:
    
    int numRoutes;
    int numAddRoutes;
    std::vector<Route*> route;
    
    SolutionVRP();
    ~SolutionVRP();
public:
    
    
    void add_empty_route();
    void add_num_empty_route(int n);
    void delete_routes();
    int get_num_routes();
    void deleteallroutes(int n);
    void  Initialize(std::vector<Veicoli*> &veic, int numVeicoli, std::vector<std::vector<double>> &Time);
    Route* copyvehicleinroute(Route* r, int j);
    void copyrouteinsol(Route* r);
    void updatecost(int numVeicoli);
    void updateusedvehicles(int numVeicoli);
    
    void addAdditionalRoute(int numVeicoli, std::vector<std::vector<double>>&Time, std::vector<Veicoli*> &veic);
    
    void CopySolution(SolutionVRP* Sol);
    void DisplaySolution();
    void delete_route(int a, int numVeicoli);
    void copyexistingsol(SolutionVRP* Sol);
    
    bool accepted(int T, SolutionVRP* Sol_Curr, int numVeicoli);
    bool bestaccepted(SolutionVRP* Sol_Curr, int numVeicoli);
    double calculate_effect(int i, int j, int k, int l, std::vector<std::vector<double>> & D);
    
    bool cap_feas(int i, int j,int k,int l,std::vector<Veicoli*> &veic);
    
    bool tw_feas(int i,int j,int k,int l, std::vector<RichiesteServizio*> &ric,std::vector<std::vector<double>> &D);
    void createnewroute(int i, int j, int k, int l, std::vector<RichiesteServizio*> &ric, std::vector<std::vector<double>>& D, std::vector<Veicoli*> &veic);
    void changetworoutes(int i, int j, int k, int l, std::vector<RichiesteServizio*> & ric , std::vector<std::vector<double>> &D,std::vector<Veicoli*> & veic);
    void ClearSol();
    bool Eightstepevaluationscheme_newr(std::vector<std::vector<double>> & Time, std::vector<RichiesteServizio*> & ric, int a, int b, int c, int d, int numRichieste, std::vector<Veicoli*> & veic);
    double calculate_effect_2(int i, int j, int k, int l, std::vector<std::vector<double>> & D);
    double calculate_dist( std::vector<std::vector<double>> & D);

    virtual Solution* clone();

};



#endif /* Solution_hpp */
