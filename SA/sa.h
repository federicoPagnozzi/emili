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
namespace emili {
namespace sa {


class SimulatedAnnealing: public emili::LocalSearch
{

protected:
    double            temp;
    double            init_temp;
    double            final_temp;
    SAInitTemp       *initialTemperature;
    SAAcceptance     *acceptanceCriterion;
    SACooling        *coolingScheme;
    SATempRestart    *temprestart;
    SATermination    *terminationCriterion;
    SATempLength     *tempLength;
    SAExploration    *exploration;
    emili::Neighborhood *neigh;


public:
    SAStatus         *sastatus;

    SimulatedAnnealing(emili::InitialSolution  *initialSolutionGenerator,
                       SAInitTemp       *initialTemperature,
                       SAAcceptance     *acceptanceCriterion,
                       SACooling        *coolingScheme,
                       SATempRestart    *temprestart,
                       SATermination    *terminationCriterion,
                       SATempLength     *tempLength,
                       SAExploration    *exploration,
                       emili::Neighborhood     *neighborhood,
                       SAStatus*         _status):
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

        status = _status;
        if (status == NULL) {
            status = new SAStatus();
          }
        sastatus = (SAStatus*) status;
        neigh = neighborhood;

        init_temp = initialTemperature->get();
        temp = init_temp;



        sastatus->set_types(terminationCriterion->getType(),
                            acceptanceCriterion->getType(),
                            tempLength->getType(),
                            temprestart->getType());
        sastatus->init_temp = init_temp;
        sastatus->final_temp = initialTemperature->getMinTemp();
        sastatus->temp = init_temp;
        acceptanceCriterion->setCurrentTemp(sastatus->temp);
        /**
         * initialization of attribute depends on termination criteria
         * but in sa_termination_criteria.h I have to include sa_common.h
         * therefore it sucks a bit but I have to initialize this here.
         */
        sastatus->tenure = 1;

        if (sastatus->tc_type == LASTACCRATETERM) {
            sastatus->tenure = terminationCriterion->getTenure();
        } else if (sastatus->tr_type == SALASTRATERESTART        ||
                   sastatus->tr_type == SALASTRATEREHEAT         ||
                   sastatus->tr_type == SALOCALMINENHANCEDREHEAT   ) {
            sastatus->tenure = temprestart->getTenure();
        }

        sastatus->last_accepted = (short *)
                malloc(sastatus->tenure * sizeof(short));

        // try at least status->tenure solutions
        // otherwise it will terminate immediately
        for (int i = 0 ; i < sastatus->tenure ; i++) {
            sastatus->last_accepted[i] = 1;
        }

        sastatus->final_temp = initialTemperature->getMinTemp();
        sastatus->init_prob = initialTemperature->getInit_prob();
        sastatus->neigh_size = neighborhood->size();

        acceptanceCriterion->set_status(sastatus);
        temprestart->set_status(sastatus);
        tempLength->set_status(sastatus);
        coolingScheme->set_status(sastatus);
        initialTemperature->set_status(sastatus);
        exploration->set_status(sastatus);

        sastatus->print();

    }

    virtual emili::Solution* search(emili::Solution* initial);
    virtual emili::Solution* search();
    
    virtual void setSearchTime(int time);

    virtual void reset();

    virtual emili::Solution* getBestSoFar() { 
      //std::cout << "in SA: " << std::fixed << sastatus->getBestSolution()->getSolutionValue() << std::endl;
      //return sastatus->best->clone();//getBestSolution();
      return sastatus->getBestSolution();
    }

    virtual void setSearchStatus(emili::SearchStatus* _status) {
      sastatus = (SAStatus *)_status;
      status = _status;
/**/      acceptanceCriterion->set_status(sastatus);
      temprestart->set_status(sastatus);
      tempLength->set_status(sastatus);
      coolingScheme->set_status(sastatus);
      initialTemperature->set_status(sastatus);
      exploration->set_status(sastatus);/**/
    }
    virtual SAStatus* getStatus(void) {
      return sastatus;
    }
    virtual SearchStatus* getSearchStatus(void) {
      return sastatus;
    }

    virtual ~SimulatedAnnealing() {
      delete initialTemperature;
      delete acceptanceCriterion;
      delete coolingScheme;
      delete temprestart;
      delete exploration;
      delete tempLength;
      /*if (status != NULL) {
        free(((SAStatus*)status)->last_accepted);
        delete (status);
      }*/
    }

}; // class SimulatedAnnealing

class SimulatedAnnealingIncumbent: public SimulatedAnnealing
{

public:

    SimulatedAnnealingIncumbent(emili::InitialSolution  *initialSolutionGenerator,
                       SAInitTemp       *initialTemperature,
                       SAAcceptance     *acceptanceCriterion,
                       SACooling        *coolingScheme,
                       SATempRestart    *temprestart,
                       SATermination    *terminationCriterion,
                       SATempLength     *tempLength,
                       SAExploration    *exploration,
                       emili::Neighborhood     *neighborhood,
                       SAStatus*         _status):
        SimulatedAnnealing(
            initialSolutionGenerator,
            initialTemperature,
            acceptanceCriterion,
            coolingScheme,
            temprestart,
            terminationCriterion,
            tempLength,
            exploration,
            neighborhood,
            _status) { }

    virtual emili::Solution* search(emili::Solution* initial);


    virtual ~SimulatedAnnealingIncumbent() {
      delete initialTemperature;
      delete acceptanceCriterion;
      delete coolingScheme;
      delete temprestart;
      delete exploration;
      delete tempLength;
      /*if (status != NULL) {
        free(((SAStatus*)status)->last_accepted);
        delete (status);
      }*/
    }

}; // class SimulatedAnnealingIncumbent

}

}
#endif
