//
//  Neighborhood.hpp
//  OrCode
//
//  Created by garazi on 25/01/2017.
//  Copyright © 2017 Garazi. All rights reserved.
//

#ifndef Neighborhood_hpp
#define Neighborhood_hpp

#include <stdio.h>
#include "Instance.hpp"
#include "SolutionVRP.hpp"
#include "Veicoli.hpp"
#include "RichiesteServizio.hpp"

SolutionVRP* Relocate_Neighborhood(SolutionVRP* Sol, int numVeicoli, int** MatCompVei, RichiesteServizio* ric, double** MatTemp, int MaxTimeRoute, Veicoli* veic, int MaxRideTime, double** D);
SolutionVRP* Relocate_Neighborhood2(SolutionVRP* Sol, int numVeicoli, int** MatCompVei, RichiesteServizio* ric, double** MatTemp, int MaxTimeRoute, Veicoli* veic, int MaxRideTime, double** D);

SolutionVRP* Exchange_F(SolutionVRP* Sol, int numVeicoli, std::vector<RichiesteServizio*> & ric, std::vector<std::vector<double>> & MatTemp,std::vector<Veicoli*> &veic, int numRichieste);
SolutionVRP* two_opt(SolutionVRP* Sol, int numVeicoli, std::vector<std::vector<double>> & D, std::vector<Veicoli*> &veic, std::vector<RichiesteServizio*> &ric, int numRichieste);

SolutionVRP* r_4_opt(SolutionVRP* Sol, int numVeicoli, std::vector<std::vector<double>> & D, std::vector<RichiesteServizio*> &ric, std::vector<Veicoli*> &veic);

SolutionVRP* Eliminate_Neighborhood(SolutionVRP* Sol, int numVeicoli, int** MatCompVei, RichiesteServizio* ric, double** MatTemp, int MaxTimeRoute, Veicoli* veic, int MaxRideTime, double** D);
SolutionVRP* Relocate_NeighborhoodF(SolutionVRP* Sol, int numVeicoli, std::vector<std::vector<int>> &MatCompVei, std::vector<RichiesteServizio*> &ric, std::vector<std::vector<double>> & MatTemp, std::vector<Veicoli*> &veic, std::vector<std::vector<double>> & D, int numRichieste);

SolutionVRP* Eliminate_NeighborhoodF(SolutionVRP* Sol, int numVeicoli, std::vector<std::vector<int>> &MatCompVei, std::vector<RichiesteServizio*> &ric, std::vector<std::vector<double>> & MatTemp, std::vector<Veicoli*> &veic, std::vector<std::vector<double>> & D, int numRichieste);

void Relocate_Neighborhood_2(SolutionVRP* Sol,int vei,int& pickup,int& delivery,int& bestRequest, int& bestVei, int numVeicoli, std::vector<std::vector<int>> &MatCompVei, std::vector<RichiesteServizio*> &ric, std::vector<std::vector<double>> & MatTemp, std::vector<Veicoli*> &veic, std::vector<std::vector<double>> & D, int numRichieste);
void two_opt_2(SolutionVRP* Sol,int vei, int& bestpos1, int& bestv2, int& bestpos2,int numVeicoli, std::vector<std::vector<double>> & D, std::vector<Veicoli*> &veic, std::vector<RichiesteServizio*> &ric, int numRichieste);


void Eliminate_NeighborhoodF_2(SolutionVRP* Sol,int vei, int numVeicoli, std::vector<std::vector<int>> &MatCompVei, std::vector<RichiesteServizio*> &ric, std::vector<std::vector<double>> & MatTemp, std::vector<Veicoli*> &veic, std::vector<std::vector<double>> & D, int numRichieste);

void r_4_opt_2(SolutionVRP* Sol, int vei, int& best_start, int& best_first, int& best_second, int& best_third, int numVeicoli, std::vector<std::vector<double>> & D, std::vector<RichiesteServizio*> &ric, std::vector<Veicoli*> &veic);
void Exchange_Neighborhood(SolutionVRP* Sol, int vei, int &best_v, int &best_r1, int& best_r2, int &or_pickup_pos_1, int &or_pickup_pos_2, int &or_delivery_pos_1, int &or_delivery_pos_2, int numVeicoli, std::vector<RichiesteServizio*> & ric, std::vector<std::vector<double>> & MatTemp,std::vector<Veicoli*> &veic, int numRichieste, std::vector<std::vector<int>> &MatCompVei, std::vector<std::vector<double>> & D);
void Exchange_Vehicles(SolutionVRP* Sol, int vei, int& vei1, int numVeicoli, int numRichieste, std::vector<std::vector<int>> &MatCompVei, std::vector<std::vector<double>> &Dist, std::vector<std::vector<double>> &Time, std::vector<RichiesteServizio*> &ric, std::vector<Veicoli*> &veic);
class RelocateNeighborhood : public emili::Neighborhood
{
protected:
    int vei; //vehicle to modify
    int num_r; // number of routes
    // reverse last move variables
    int pickup; // pickup of last relocation
    int delivery; // delivery of last relocation
    int bestRequest; // best request of last relocation
    int bestVei; // best vehicle of last relocation
    Instance& inst;
    virtual emili::Solution* computeStep(emili::Solution* step);
    virtual void reverseLastMove(emili::Solution* step);
public:
    RelocateNeighborhood(Instance& instance):inst(instance) {}
    virtual NeighborhoodIterator begin(emili::Solution* base);
    virtual void reset();
    virtual emili::Solution* random(emili::Solution* currentSolution);
    virtual int size() { return num_r;}
};

class TwoOptNeighborhood: public emili::Neighborhood
{
 protected:
    int vei;//vehicle to modify route1
    int num_r; //number of routes
    // reverse last move variables
    int bestpos1; //best position of the route1 of last change
    int bestv2; //best route2 to exchange of last change
    int bestpos2;//best position of the route2  of last change
    Instance& inst;
    virtual emili::Solution* computeStep(emili::Solution* step);
    virtual void reverseLastMove(emili::Solution* step);


 public:
    TwoOptNeighborhood(Instance& instance):inst(instance) {}
    virtual NeighborhoodIterator begin(emili::Solution* base);
    virtual void reset();
    virtual emili::Solution* random(emili::Solution* currentSolution);
    virtual int size() { return num_r;}

};

class FourOptNeighborhood: public emili::Neighborhood
{
protected:
    int vei; //vehicle to modify
    int num_r; // number of routes
    // reverse last move variables
    int best_start; // where the change has started of last change
    int best_first; // first node is in location best_first at last change
    int best_second; // second node is in location best_second at last change
    int best_third; //third node is in location best_third at last change
    Instance& inst;
    virtual emili::Solution* computeStep(emili::Solution* step);
    virtual void reverseLastMove(emili::Solution* step);
public:
    FourOptNeighborhood(Instance& instance):inst(instance) {}
    virtual NeighborhoodIterator begin(emili::Solution* base);
    virtual void reset();
    virtual emili::Solution* random(emili::Solution* currentSolution);
    virtual int size() { return num_r;}
};

class EliminateNeighborhood : public emili::Neighborhood
{
protected:
    int vei; //vehicle to modify
    int num_r; // number of routes
    SolutionVRP* base_solution;
    Instance& inst;
    virtual emili::Solution* computeStep(emili::Solution* step);
    virtual void reverseLastMove(emili::Solution* step);
public:
    EliminateNeighborhood(Instance& instance):inst(instance),base_solution(new SolutionVRP()) {}
    virtual NeighborhoodIterator begin(emili::Solution* base);
    virtual void reset();
    virtual emili::Solution* random(emili::Solution* currentSolution);
    virtual int size() { return num_r;}
};

class ExchangeNeighborhood : public emili::Neighborhood
{
protected:
    int vei; //vehicle to modify
    int num_r; // number of routes

    int best_v; //best vehicle in where the request to exchange is
    int best_r1; //request to exchange from first vehicle (vei)
    int best_r2; //request to exchange from second vehicle (best_v)
    int or_pickup_pos_1;//where in first vehicle is the second request inserted pickup
    int or_pickup_pos_2;//where in first vehicle is the second request inserted delivery
    int or_delivery_pos_1;//where in first vehicle is the second request inserted pickup
    int or_delivery_pos_2;//where in first vehicle is the second request inserted delivery
    SolutionVRP* base_solution;
    Instance& inst;
    virtual emili::Solution* computeStep(emili::Solution* step);
    virtual void reverseLastMove(emili::Solution* step);
public:
    ExchangeNeighborhood(Instance& instance):inst(instance),base_solution(new SolutionVRP()) {}
    virtual NeighborhoodIterator begin(emili::Solution* base);
    virtual void reset();
    virtual emili::Solution* random(emili::Solution* currentSolution);
    virtual int size() { return num_r;}
};

class ExchangeVehicleNeighborhood : public emili::Neighborhood
{
protected:
    int vei; //vehicle to modify
    int num_r; // number of routes

    int vei1;
    SolutionVRP* base_solution;
    Instance& inst;
    virtual emili::Solution* computeStep(emili::Solution* step);
    virtual void reverseLastMove(emili::Solution* step);
public:
    ExchangeVehicleNeighborhood(Instance& instance):inst(instance),base_solution(new SolutionVRP()) {}
    virtual NeighborhoodIterator begin(emili::Solution* base);
    virtual void reset();
    virtual emili::Solution* random(emili::Solution* currentSolution);
    virtual int size() { return num_r;}
};

#endif /* Neighborhood_hpp */

