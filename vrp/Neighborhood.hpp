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

#endif /* Neighborhood_hpp */

