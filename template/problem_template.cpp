//
//  Created by Federico Pagnozzi on 12/12/17.
//  Copyright (c) 2014 Federico Pagnozzi. All rights reserved.
//  This file is distributed under the BSD 2-Clause License. See LICENSE.TXT
//  for details.

#include "problem_template.h"
void notyet()
{
    std::cout << "Not implemented yet \n";
    //std::exit(-5);
}
   
      double emili::problemX::ProblemX::calcObjectiveFunctionValue(Solution& solution)
      {
            notyet();
            return 0;
      }
 
      double emili::problemX::ProblemX::evaluateSolution(Solution & solution)
       {
        notyet();
        return 0;
    }
 

  
      const void* emili::problemX::SolutionProblemX::getRawData()const
       {
        notyet();
        return nullptr;
    }
 
      void emili::problemX::SolutionProblemX::setRawData(const void* data)
      {
        notyet();
 
    }
 
      std::string emili::problemX::SolutionProblemX::getSolutionRepresentation()
      {

        return std::string("TO implement");
    }
 
      emili::Solution* emili::problemX::SolutionProblemX::clone()
    {
        notyet();
        return nullptr;
    }
 
 
      emili::Solution* emili::problemX::InitialSolutionProblemX::generateSolution()
    {
        notyet();
        return nullptr;
    }
 
 
      emili::Solution* emili::problemX::InitialSolutionProblemX::generateEmptySolution()
    {
        notyet();
        return nullptr;
    }
 

      emili::Solution* emili::problemX::NeighborhoodProblemX::computeStep(Solution* step)
    {
        notyet();
        return nullptr;
    }
 

      void emili::problemX::NeighborhoodProblemX::reverseLastMove(Solution* step)
    {
        notyet();

    }
 

      emili::Neighborhood::NeighborhoodIterator emili::problemX::NeighborhoodProblemX::begin(emili::Solution* base)
    {
        notyet();
        return emili::Neighborhood::NeighborhoodIterator(this,nullptr);
    }
 

      void emili::problemX::NeighborhoodProblemX::reset()
      {

      }

      emili::Solution* emili::problemX::NeighborhoodProblemX::random(Solution* currentSolution)
      {
          notyet();
          return nullptr;
      }

      /*int emili::problemX::NeighborhoodProblemX::size()
      {
          notyet();
          return 0;
      }*/
     
      emili::Solution* emili::problemX::PerturbationProblemX::perturb(Solution* solution)
      {
          notyet();
          return nullptr;
      }
