//
//  Instance.hpp
//  C-B
//
//  Created by garazi on 10/03/2017.
//  Copyright © 2017 Garazi. All rights reserved.
//

#ifndef Instance_hpp
#define Instance_hpp
#include <vector>
#include <stdio.h>
#include "RichiesteServizio.hpp"
#include "Veicoli.hpp"
#include <string>
#include "../emilibase.h"
#include "SolutionVRP.hpp"

class Instance : public emili::Problem
{
public:
    int id;
    std::string name;
    int numRichieste0;
    int numVeicoli0;
    int numLocation0;
    
    std::vector<RichiesteServizio*> rc;
    std::vector<Veicoli*> vec;
    
    std::vector< std::vector<double> > T;
    std::vector<std::vector<double>> Dist;
    std::vector<std::vector<int>> MCV;

    Instance();
    ~Instance();
public:
    void read_Distance(std::string a);
    void read_Time(std::string a);
    void read_MatCompVei(std::string a);
    void read_ric(std::string a);
    void read_veic(std::string a);
    void read_instance(std::string a, int i);
    
    void read_num(std::string);
    
    void TimeWindowTightening();

    virtual double calcObjectiveFunctionValue(emili::Solution& solution);
    virtual double evaluateSolution(emili::Solution & solution);
    virtual int problemSize();

};
#endif /* Instance_hpp */
