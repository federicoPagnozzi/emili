//
//  InitialSolution.hpp
//  OrCode
//
//  Created by garazi on 25/01/2017.
//  Copyright © 2017 Garazi. All rights reserved.
//

#ifndef InitialSolution_hpp
#define InitialSolution_hpp

#include <stdio.h>

#include "Veicoli.hpp"
#include "RichiesteServizio.hpp"
#include "SolutionVRP.hpp"
#include "Instance.hpp"

SolutionVRP* InitialSolutionBraekersF2(SolutionVRP* InitialSol, int numVeicoli, int numRichieste, std::vector<RichiesteServizio*> &ric, std::vector<Veicoli*> &veic, std::vector<std::vector<double>> &D, std::vector<std::vector<int>> &MatCompVei, std::vector<std::vector<double>> &MatTemp);

class IS :public emili::InitialSolution
{
protected:
    Instance& inst;
public:
    IS(Instance& instance):emili::InitialSolution(instance),inst(instance) {}

    virtual emili::Solution* generateSolution();


    virtual emili::Solution* generateEmptySolution();

};


/*
Solution* InitialSolutionBraekers(int numVeicoli, int numRichieste, RichiesteServizio* ric, Veicoli* veic, double** D, int** MatCompVei, double** MatTemp, double MinAt, int MaxTimeRoute, double MaxRideTime);
Solution* InitialSolutionBraekers2(int numVeicoli, int numRichieste, RichiesteServizio* ric, Veicoli* veic, double** D, int** MatCompVei, double** MatTemp, double MinAt, int MaxTimeRoute, double MaxRideTime);
Solution* InitialSolutionBraekersF(Solution* InitialSol, int numVeicoli, int numRichieste, RichiesteServizio* ric, Veicoli* veic, double** D, int** MatCompVei, double** MatTemp, double MinAt, int MaxTimeRoute, double MaxRideTime);*/

#endif /* InitialSolution_hpp */
