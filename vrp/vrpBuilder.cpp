//
//  Created by Federico Pagnozzi on 28/11/14.
//  Copyright (c) 2014 Federico Pagnozzi. All rights reserved.
//  This file is distributed under the BSD 2-Clause License. See LICENSE.TXT
//  for details.
#include "vrpBuilder.h"
#include <cstdlib>
#include <cstdio>
#include <ctime>
#include <cstring>
#include <iostream>
#include <sstream>
#include <algorithm>
#include "Instance.hpp"
#include "Neighborhood.hpp"

/* Algos */

/* tabu tenure types */

/* VRP PROBLEM*/
#define PROBLEM_VRP "VRP"

/* initial solution heuristics */
#define INITIAL_IS "is"


/* Termination criteria*/


/*
 *  VRP neighborhoods
 *
 */
#define NEIGHBORHOOD_RELOCATE "relocate"
#define NEIGHBORHOOD_ELIMINATE "eliminate"
#define NEIGHBORHOOD_TWO_OPT "topt"
#define NEIGHBORHOOD_FOUR_OPT "fopt"
#define NEIGHBORHOOD_EXCHANGE "exchange"
#define NEIGHBORHOOD_EXCHANGE_VEHICLE "vehicle_exchange"



/*
 * END Neighborhoods
 */


/* permutation flowshop solution perturbations */

/* acceptance criteria*/

/*void prs::emili_header()
{
    std::cout << "\t ______ __  __ _____ _      _____ " << std::endl;
    std::cout << "\t|  ____|  \\/  |_   _| |    |_   _|" << std::endl;
    std::cout << "\t| |__  | \\  / | | | | |      | |  " << std::endl;
    std::cout << "\t|  __| | |\\/| | | | | |      | |  " << std::endl;
    std::cout << "\t| |____| |  | |_| |_| |____ _| |_ " << std::endl;
    std::cout << "\t|______|_|  |_|_____|______|_____|" << std::endl;
    std::cout << std::endl;
}*/

emili::InitialSolution* prs::VrpBuilder::buildInitialSolution()
{
    prs::incrementTabLevel();
    std::ostringstream oss;
    emili::InitialSolution* init = nullptr;
    Instance* instance =(Instance*) gp.getInstance();
    if(tm.checkToken(INITIAL_IS))
    {
        printTab("IS initial solution");
        init = new IS(*instance);
    }
    prs::decrementTabLevel();
    return init;
}

emili::Neighborhood* prs::VrpBuilder::buildNeighborhood()
{
    prs::incrementTabLevel();
    emili::Neighborhood* neigh = nullptr;
    Instance* instance =(Instance*) gp.getInstance();
    if(tm.checkToken(NEIGHBORHOOD_RELOCATE))
    {
        printTab( "Relocate Neighborhood");
        neigh = new RelocateNeighborhood(*instance);
    }
    else if(tm.checkToken(NEIGHBORHOOD_TWO_OPT))
    {
        printTab(" Two-Opt Neighborhood");
        neigh= new TwoOptNeighborhood(*instance);
    }
    else if(tm.checkToken(NEIGHBORHOOD_FOUR_OPT))
    {
        printTab("Four-Opt Neighborhood");
        neigh=new FourOptNeighborhood(*instance);
    }else if(tm.checkToken(NEIGHBORHOOD_ELIMINATE))
    {
        printTab("Eliminate Neighborhood");
        neigh=new EliminateNeighborhood(*instance);
    }else if(tm.checkToken(NEIGHBORHOOD_EXCHANGE))
    {
        printTab("Eliminate Neighborhood");
        neigh=new ExchangeNeighborhood(*instance);
    }else if(tm.checkToken(NEIGHBORHOOD_EXCHANGE_VEHICLE))
    {
        printTab("Eliminate Neighborhood");
        neigh=new ExchangeVehicleNeighborhood(*instance);
    }
    prs::decrementTabLevel();
    return neigh;
}
emili::Problem* prs::VrpBuilder::buildProblem()
{

    return gp.getInstance();
}

emili::Problem* prs::VrpBuilder::openInstance()
{
    //a2-16hetIUY
    std::string instance_string = tm.tokenAt(1);
    Instance* inst=nullptr;
    std::cout << instance_string << "\n";
    inst = new Instance();
    inst->read_instance(instance_string, 1);
    inst->TimeWindowTightening();

    tm.next();
    return inst;
}

bool prs::VrpBuilder::isParsable(std::string &problem)
{
    if(strcmp(problem.c_str(),PROBLEM_VRP)==0)
    {
        return true;
    }
    return false;
}

bool prs::VrpBuilder::isCompatibleWith(char *problem_definition)
{
    std::string s(problem_definition);
    return isParsable(s);
}

bool prs::VrpBuilder::canOpenInstance(char *problem_definition)
{
    std::string s(problem_definition);
    return isParsable(s);
}


