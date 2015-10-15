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
    emili::Neighborhood *neigh;
    SAStatus         *status;


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

                        neigh = neighborhood;

                        init_temp = initialTemperature->get();
                        temp = init_temp;

                        //status = (sa_status *)malloc(sizeof(sa_status));
                        status = new SAStatus();

                        /**
                         * initialization of attribute depends on termination criteria
                         * but in sa_termination_criteria.h I have to include sa_common.h
                         * therefore I have to initialize this here.
                         */
                        if (terminationCriterion->getType() == LASTACCRATETERM) {
                          status->tenure = terminationCriterion->getTenure();
                          status->last_accepted = (short *)
                              malloc(status->tenure * sizeof(short));

                          // try at least status->tenure solutions
                          // otherwise it will terminate immediately
                          for (int i = 0 ; i < status->tenure ; i++) {
                            status->last_accepted[i] = 1;
                          }

                        }

                        std::string tc_type = terminationCriterion->getType();
                        std::string ac_type = acceptanceCriterion->getType();
                        std::string tl_type = tempLength->getType();

                        acceptanceCriterion->set_status(status);
                      }

    virtual emili::Solution* search(emili::Solution* initial);
    virtual emili::Solution* search();

    virtual void reset();

    virtual ~SimulatedAnnealing() {
      delete initialTemperature;
      delete acceptanceCriterion;
      delete coolingScheme;
      delete exploration;
      delete tempLength;
      if (terminationCriterion->getType() == LASTACCRATETERM) {
        free(status->last_accepted);
      }
      delete (status);
    };

}; // class SimulatedAnnealing

#endif
