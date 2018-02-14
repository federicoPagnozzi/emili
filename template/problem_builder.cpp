//
//  Created by Federico Pagnozzi on 12/12/17.
//  Copyright (c) 2014 Federico Pagnozzi. All rights reserved.
//  This file is distributed under the BSD 2-Clause License. See LICENSE.TXT
//  for details.

#include "problem_builder.h"
#include "problem_template.h"

#define PROBLEMX "px"
#define INITIAL_PROBLEMX "ipx"
#define NEIGHBORHOOD_PROBLEMX "npx"
#define PERTURBATION_PROBLEMX "ppx"

    
    bool prs::problemX::ProblemXBuilder::isCompatibleWith(char* problem_definition)
    {
        if(strcmp(problem_definition,PROBLEMX)==0)
        {
            return true;
        }
        return false;
    }
    
    bool prs::problemX::ProblemXBuilder::canOpenInstance(char* problem_definition)
    {
        if(strcmp(problem_definition,PROBLEMX)==0)
        {
            return true;
        }
        return false;
    }
    
    emili::Problem* prs::problemX::ProblemXBuilder::openInstance()
    {
        //instance file path ( can be relative )
        char* problem_string = tm.nextToken();
        notyet();
        return nullptr;
    }
    
    
    emili::InitialSolution* prs::problemX::ProblemXBuilder::buildInitialSolution()
    {
        prs::incrementTabLevel();           
        emili::InitialSolution* init = nullptr;
        emili::problemX::ProblemX* instance =(emili::problemX::ProblemX*) gp.getInstance();
        if(tm.checkToken(INITIAL_PROBLEMX))
        {
            printTab( "PROBLEMX initial solution");
          //  int index = tm.getInteger();
            init = new emili::problemX::InitialSolutionProblemX(*instance);
        }
         prs::decrementTabLevel();

        return init;

    }
    emili::Neighborhood* prs::problemX::ProblemXBuilder::buildNeighborhood()
    {           
            prs::incrementTabLevel();           
            emili::Neighborhood* neigh = nullptr;
            //emili::problemX::ProblemX* instance =(emili::problemX::ProblemX*) gp.getInstance();
            if(tm.checkToken(NEIGHBORHOOD_PROBLEMX))
            {
                printTab( "PROBLEMX Neighborhood");
                neigh = new emili::problemX::NeighborhoodProblemX();
            }
             prs::decrementTabLevel();

            return neigh;
    }
    emili::Perturbation* prs::problemX::ProblemXBuilder::buildPerturbation()
    {
         prs::incrementTabLevel();           
            emili::Perturbation* per = nullptr;
            if(tm.checkToken(PERTURBATION_PROBLEMX))
            {
                printTab( "PROBLEMX Perturbation");
                //emili::Neighborhood* n = gp.buildComponent(COMPONENT_NEIGHBORHOOD).get<emili::Neighborhood>();
                per = new emili::problemX::PerturbationProblemX();
            }
             prs::decrementTabLevel();

            return per;
    }  

/*extern "C" {
    prs::Builder* getBuilder(prs::GeneralParserE* ge)
    {
        return new prs::problemX::ProblemXBuilder(*ge,ge->getTokenManager());
    }
}
*/
