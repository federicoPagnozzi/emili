// #ifndef METROPOLIS_H
// #define METROPOLIS_H

// #include <string>
// #include <iostream>
// #include <vector>
// #include <functional>

// #include "sa.h"

// #include "../emilibase.h"
// namespace emili {
// namespace metropolis {


// class Metropolis: public emili::LocalSearch
// {

// protected:
//     SimulatedAnnealing*  ftsa; // fixed temperature SA
//     SimulatedAnnealing*  dsa;  // diversification SA (accept all)
//     SATermination*       terminationCriterion
//     SAStatus*            status;


// public:
//     Metropolis(SimulatedAnnealing* _ftsa,
//                SimulatedAnnealing* _dsa,
//                SATermination*      _terminationCriterion):
//                       ftsa(_ftsa),
//                       dsa(_dsa),
//                       terminationCriterion(_terminationCriterion) {

//                         neigh = neighborhood;

//                         init_temp = initialTemperature->get();
//                         temp = init_temp;

//                         status = new SAStatus();

//                         dsa->set_status(status);
//                         ftsa->set_status(status);

//                         status->print();

//                         status->set_types(terminationCriterion->getType(),
//                                           acceptanceCriterion->getType(),
//                                           tempLength->getType(),
//                                           temprestart->getType());
//                         status->init_temp = init_temp;
//                         status->final_temp = initialTemperature->getMinTemp();
//                         status->temp = init_temp;
//                         acceptanceCriterion->setCurrentTemp(status->temp);
                        

//                         /**
//                          * initialization of attribute depends on termination criteria
//                          * but in sa_termination_criteria.h I have to include sa_common.h
//                          * therefore it sucks a bit but I have to initialize this here.
//                          */
//                         status->tenure = 1;

//                         if (status->tc_type == LASTACCRATETERM) {
//                           status->tenure = terminationCriterion->getTenure();
//                         } else if (status->tr_type == SALASTRATERESTART        ||
//                                    status->tr_type == SALASTRATEREHEAT         ||
//                                    status->tr_type == SALOCALMINENHANCEDREHEAT   ) {
//                           status->tenure = temprestart->getTenure();
//                         }

//                         status->last_accepted = (short *)
//                             malloc(status->tenure * sizeof(short));

//                         // try at least status->tenure solutions
//                         // otherwise it will terminate immediately
//                         for (int i = 0 ; i < status->tenure ; i++) {
//                           status->last_accepted[i] = 1;
//                         }

//                         status->final_temp = initialTemperature->getMinTemp();
//                         status->init_prob = initialTemperature->getInit_prob();
//                         status->neigh_size = neighborhood->size();

//                         acceptanceCriterion->set_status(status);
//                         //temprestart->set_status(status);
//                         tempLength->set_status(status);
//                         coolingScheme->set_status(status);
//                         initialTemperature->set_status(status);

//                         status->print();
//                       }

//     virtual emili::Solution* search(emili::Solution* initial);
//     virtual emili::Solution* search();
    
//     virtual void setSearchTime(int time);

//     virtual void reset();

//     virtual emili::Solution* getBestSoFar() { return status->best;}

//     virtual ~Metropolis() {
//       delete ftsa;
//       delete dsa;
//       if (status != NULL) {
//         free(status->last_accepted);
//         delete (status);
//       }
//     }

// }; // class Metropolis

// }

// }
// #endif
