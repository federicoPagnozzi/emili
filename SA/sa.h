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
#include "sa_temperature_restart.h"

#include "../emilibase.h"


class SimulatedAnnealing: public emili::LocalSearch
{

protected:
    double            temp;
    double            init_temp;
    SAInitTemp       *initialTemperature;
    SAAcceptance     *acceptanceCriterion;
    SACooling        *coolingScheme;
    SATempRestart    *temprestart;
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
                       SATempRestart    *temprestart,
                       SATermination    *terminationCriterion,
                       SATempLength     *tempLength,
                       SAExploration    *exploration,
                       emili::Neighborhood     *neighborhood):
                      initialTemperature(initialTemperature),
                      acceptanceCriterion(acceptanceCriterion),
                      coolingScheme(coolingScheme),
                      temprestart(temprestart),
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

                        status->set_types(terminationCriterion->getType(),
                                          acceptanceCriterion->getType(),
                                          tempLength->getType(),
                                          temprestart->getType());

                        /**
                         * initialization of attribute depends on termination criteria
                         * but in sa_termination_criteria.h I have to include sa_common.h
                         * therefore it sucks a bit but I have to initialize this here.
                         */
                        /**
                         * dummy values that will never be touched unless needed,
                         * or unless I fuck up something
                         */
                        status->tenure = 1;

                        if (status->tc_type == LASTACCRATETERM) {
                          status->tenure = terminationCriterion->getTenure();
                        } else if (status->tr_type == SALASTRATERESTART        ||
                                   status->tr_type == SALASTRATEREHEAT         ||
                                   status->tr_type == SALOCALMINENHANCEDREHEAT   ) {
                          status->tenure = temprestart->getTenure();
                        }

                        status->last_accepted = (short *)
                            malloc(status->tenure * sizeof(short));

                        // try at least status->tenure solutions
                        // otherwise it will terminate immediately
                        for (int i = 0 ; i < status->tenure ; i++) {
                          status->last_accepted[i] = 1;
                        }

                        acceptanceCriterion->set_status(status);
                        temprestart->set_status(status);
                        tempLength->set_status(status);
                        coolingScheme->set_status(status);
                      }

    virtual emili::Solution* search(emili::Solution* initial);
    virtual emili::Solution* search();

    virtual void reset();

    virtual ~SimulatedAnnealing() {
      delete initialTemperature;
      delete acceptanceCriterion;
      delete coolingScheme;
      delete temprestart;
      delete exploration;
      delete tempLength;
      free(status->last_accepted);
      delete (status);
    };

}; // class SimulatedAnnealing

#endif
