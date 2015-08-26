/*
 * =====================================================================================
 *
 *       Filename:  sa.h
 *
 *    Description:  Simulated Annealing
 *
 *        Version:  1.0
 *        Created:  06/07/2015 16:16:42
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *        Company:  
 *
 * =====================================================================================
 */

#ifndef SA_H
#define SA_H

#include <string>
#include <iostream>
#include <vector>
#include <functional>

#include "sa_constants.h"
#include "sa_acceptance_criteria.h"
#include "sa_cooling.h"
#include "sa_termination_criteria.h"
#include "sa_neighborhood.h"
#include "sa_init_temp.h"

 #include "../emilibase.h"

#define nullptr NULL

class SimulatedAnnealing: public emili::LocalSearch
{

protected:
    double            temp;
    double            init_temp;
    SAInitTemp       *initialTemperature;
    SAAcceptance     *acceptanceCriterion;
    SACooling        *coolingScheme;
    SATermination    *terminationCriterion;


public:
    SimulatedAnnealing(emili::InitialSolution  *initialSolutionGenerator,
                       SAInitTemp       *initialTemperature,
                       SAAcceptance     *acceptanceCriterion,
                       SACooling        *coolingScheme,
                       SATermination      *terminationCriterion,
                       emili::Neighborhood     *neighborhood):
                      initialTemperature(initialTemperature),
                      acceptanceCriterion(acceptanceCriterion),
                      coolingScheme(coolingScheme),
                      terminationCriterion(terminationCriterion),
                      emili::LocalSearch(*initialSolutionGenerator,
                                         *terminationCriterion,
                                         *neighborhood) { }

    virtual emili::Solution* search(emili::Solution* initial);
    virtual emili::Solution* search();

    virtual void reset();

    virtual ~SimulatedAnnealing() {
      delete initialTemperature;
      delete acceptanceCriterion;
      delete coolingScheme;
    };

}; // class SimulatedAnnealing

#endif
