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

#include "sa_common.h"
#include "sa_exploration.h"
#include "sa_constants.h"
#include "sa_acceptance_criteria.h"
#include "sa_cooling.h"
#include "sa_termination_criteria.h"
#include "sa_neighborhood.h"
#include "sa_init_temp.h"
#include "sa_templength.h"

#include "../emilibase.h"


class SimulatedAnnealing: public emili::LocalSearch
{

protected:
    double            temp;
    double            init_temp;
    SAInitTemp       *initialTemperature;
    SAAcceptance     *acceptanceCriterion;
    SACooling        *coolingScheme;
    SATermination    *terminationCriterion;
    SATempLength     *tempLength;
    SAExploration    *exploration;
    sa_status        *status;


public:
    SimulatedAnnealing(emili::InitialSolution  *initialSolutionGenerator,
                       SAInitTemp       *initialTemperature,
                       SAAcceptance     *acceptanceCriterion,
                       SACooling        *coolingScheme,
                       SATermination    *terminationCriterion,
                       SATempLength     *tempLength,
                       SAExploration    *exploration,
                       emili::Neighborhood     *neighborhood):
                      initialTemperature(initialTemperature),
                      acceptanceCriterion(acceptanceCriterion),
                      coolingScheme(coolingScheme),
                      terminationCriterion(terminationCriterion),
                      exploration(exploration),
                      tempLength(tempLength),
                      emili::LocalSearch(*initialSolutionGenerator,
                                         *terminationCriterion,
                                         *neighborhood) {
                        status = (sa_status *)malloc(sizeof(sa_status));
                        status->counter = 0;
                        status->total_counter = 0;
                        status->accepted = 0;

                        if (terminationCriterion->getType() == LASTACCRATETERM) {
                          status->tenure = terminationCriterion->getTenure();
                          status->last_accepted = (short *)
                              malloc(status->tenure * sizeof(short));

                          // try at least status->tenure solutions
                          // otherwise it will terminate immediately
                          for (int i = 0 ; i < status->tenure ; i++) {
                            status->last_accepted[i] = 1;
                          }

                          status->index = 0;
                        }
                      }

    virtual emili::Solution* search(emili::Solution* initial);
    virtual emili::Solution* search();

    virtual void reset();

    virtual ~SimulatedAnnealing() {
      delete initialTemperature;
      delete acceptanceCriterion;
      delete coolingScheme;
      if (terminationCriterion->getType() == LASTACCRATETERM) {
        free(status->last_accepted);
      }
      free (status);
    };

}; // class SimulatedAnnealing

#endif
