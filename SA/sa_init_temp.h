#ifndef SA_INIT_TEMP_H
#define SA_INIT_TEMP_H

#include <string>
#include <iostream>
#include <vector>
#include <functional>
#include <cmath>

#include "sa_common.h"

#include "../emilibase.h"


/**
 * Generic initial temperature setup for SA.
 *
 * The actual scheme has to be implemented with derived classes.
 */
class SAInitTemp {

protected:
    emili::Solution  *solution;
    double            init_temp;
    SAStatus         *status;
    double            maxdelta,
                      mindelta;

public:

    /**
     * Empty constructor
     */
    SAInitTemp(void):
        solution(nullptr) { }
    /**
     * Constructor: set initial temperature starting from a valid solution.
     */
    SAInitTemp(emili::Solution *solution):
        solution(solution) { }


    /**
     * set initial temperature.
     *
     * If the initial temperature is set by the user,
     * value is the starting value for the temperature.
     *
     * If the initial temperature is set starting from a solution,
     * then value is a coefficient that rescales the cost of the solution.
     * 
     * @param  value initial temperature or coefficient.
     * @return       initial temperature
     */
    virtual double set(double value)=0;
    virtual double getMinTemp(void)=0;

    virtual double get(void) {
        return init_temp;
    }

    void set_status(SAStatus* _status) {
        status = _status;
    }

    virtual double getInit_prob(void) {
        return 1;
    }

}; // class SAInitTemp


/************************************************
 *
 * Implementation of subclasses
 * 
 ************************************************/


/**
 * Set the initial temperature at a fixed value.
 *
 * It just uses the value provieded by the user.
 */
class FixedInitTemp: public SAInitTemp {
public:

    FixedInitTemp(void):
        SAInitTemp() { }

    virtual double set(double value) {
        init_temp = value;
        return init_temp;
    }

    virtual double getMinTemp(void) {
        return 0.0001;
    }

}; // FixedInitTemp


/**
 * Set the initial temperature starting from an initial solution.
 *
 * The initial temperature is computed according to the formula:
 * init_temp = coeff * cost(solution)
 */
class InitTempFromSolution: public SAInitTemp {
public:

    InitTempFromSolution(emili::Solution *solution):
        SAInitTemp(solution) { }

    virtual double set(double value) {
        init_temp = value * solution->getSolutionValue();;
        return init_temp;
    }

    virtual double getMinTemp(void) {
        return 0.0001;
    }

}; // InitTempFromSolution


/**
 * Do a random walk and take the (absolute value of the) highest gap
 * as initial temperature
 * Minimum temperature is the smallest non-zero gap.
 */
class RandomWalkInitTemp: public SAInitTemp {
protected:
    emili::InitialSolution *is;
    emili::Neighborhood *nei;
    int length;

public:

    RandomWalkInitTemp(emili::InitialSolution* _is,
                       emili::Neighborhood *_nei,
                       int _length):
        is(_is),
        nei(_nei),
        length(_length) { }

    virtual double set(double value) {
        int i;

        emili::Solution *s1;
        emili::Solution *s2 = is->generateSolution();
        
        double c1,
               c2 = s2->getSolutionValue();
        double bestcost = c2;

        maxdelta = 0;
        mindelta = c2;

        for (i = 0 ; i < length ; i++) {
            s1 = s2;
            s2 = nei->random(s1); //is->generateSolution();
            c1 = c2;
            c2 = s2->getSolutionValue();

            if (abs(c2 - c1) > maxdelta) {
                maxdelta = abs(c2 - c1);
            } else if (abs(c2 - c1) > 0 && abs(c2 - c1) < mindelta) {
                mindelta = abs(c2 - c1);
            }

            delete s1;
        }

        delete s2;

        init_temp = value * maxdelta;

        return init_temp;
    }

    virtual double getMinTemp(void) {
        return mindelta;
    }

    virtual double getInit_prob(void) {
        return std::exp(maxdelta / init_temp);
    }

}; // RandomWalkInitTemp



class ConnollyRandomWalkInitTemp: public SAInitTemp {
protected:
    emili::InitialSolution *is;
    emili::Neighborhood *nei;
    int length;

public:

    ConnollyRandomWalkInitTemp(emili::InitialSolution* _is,
                               emili::Neighborhood *_nei,
                               int _length):
        is(_is),
        nei(_nei),
        length(_length) { }

    virtual double set(double value) {
        int i;

        emili::Solution *s1;
        emili::Solution *s2 = is->generateSolution();
        
        double c1,
               c2 = s2->getSolutionValue();
        double bestcost = c2;

        maxdelta = 0;
        mindelta = c2;

        for (i = 0 ; i < length ; i++) {
            s1 = s2;
            s2 = nei->random(s1); //is->generateSolution();
            c1 = c2;
            c2 = s2->getSolutionValue();

            if (abs(c2 - c1) > maxdelta) {
                maxdelta = abs(c2 - c1);
            } else if (abs(c2 - c1) > 0 && abs(c2 - c1) < mindelta) {
                mindelta = abs(c2 - c1);
            }

            delete s1;
        }

        delete s2;

        init_temp = value * (mindelta + (maxdelta - mindelta) / 10);
        
        return init_temp;
    }

    virtual double getMinTemp(void) {
        return mindelta;
    }

    virtual double getInit_prob(void) {
        return std::exp(maxdelta / init_temp);
    }

}; // ConnollyRandomWalkInitTemp


/**
 * Do a random walk and take the average of the (absolute values of the) gaps
 * as initial temperature
 */
class RandomWalkAvgInitTemp: public SAInitTemp {
protected:
    emili::InitialSolution *is;
    emili::Neighborhood *nei;
    int length;

public:

    RandomWalkAvgInitTemp(emili::InitialSolution* _is,
                       emili::Neighborhood *_nei,
                       int _length):
        is(_is),
        nei(_nei),
        length(_length) { }

    virtual double set(double value) {
        int i;

        emili::Solution *s1;
        emili::Solution *s2 = is->generateSolution();
        
        double c1,
               c2 = s2->getSolutionValue();
        double bestcost = c2,
               costsum = 0;

        for (i = 0 ; i < length ; i++) {
            s1 = s2;
            s2 = nei->random(s1);; //is->generateSolution();
            c1 = c2;
            c2 = s2->getSolutionValue();
            costsum += abs(c2 - c1);

            delete s1;
        }

        delete s2;

        init_temp = value * costsum / length;

        return init_temp;
    }

    virtual double getMinTemp(void) {
        return 0.0001;
    }

    virtual double getInit_prob(void) {
        return std::exp(maxdelta / init_temp);
    }

}; // RandomWalkAvgInitTemp



/**
 * t = avg(delta) / ln(desired init prob9)
 */
class RandomWalkInitProb: public SAInitTemp {
protected:
    emili::InitialSolution *is;
    emili::Neighborhood *nei;
    float init_prob;
    int length;

public:

    RandomWalkInitProb(emili::InitialSolution* _is,
                       emili::Neighborhood *_nei,
                       float _init_prob,
                       int _length):
        is(_is),
        nei(_nei),
        init_prob(_init_prob),
        length(_length) { }

    virtual double set(double value) {
        int i;

        emili::Solution *s1;
        emili::Solution *s2 = is->generateSolution();
        
        double c1,
               c2 = s2->getSolutionValue();
        double bestcost = c2,
               costsum = 0;

        for (i = 0 ; i < length ; i++) {
            s1 = s2;
            s2 = nei->random(s1);; //is->generateSolution();
            c1 = c2;
            c2 = s2->getSolutionValue();
            costsum += abs(c2 - c1);

            delete s1;
        }

        delete s2;

        init_temp = std::abs(value * (costsum / length) / std::log(init_prob));

        return init_temp;
    }

    virtual double getMinTemp(void) {
        return 0.0001;
    }

    virtual double getInit_prob(void) {
        return init_prob;
    }

}; // RandomWalkInitProb



/**
 * misevicius - new sa
 */
class MiseviciusInitTemp: public SAInitTemp {
protected:
    emili::InitialSolution *is;
    emili::Neighborhood *nei;
    int length;
    float l11, l12, l21, l22;
    double avgdelta, ti, tf;

public:

    MiseviciusInitTemp(emili::InitialSolution* _is,
                       emili::Neighborhood *_nei,
                       int _length,
                       float _l11,
                       float _l12,
                       float _l21,
                       float _l22):
        is(_is),
        nei(_nei),
        length(_length),
        l11(_l11),
        l12(_l12),
        l21(_l21),
        l22(_l22) {
            if (l11 + l12 > 1) {
                l11 = l11 / (l11 + l12);
                l12 = l12 / (l11 + l12);
            }
            if (l21 + l22 > 1) {
                l21 = l21 / (l21 + l22);
                l22 = l22 / (l21 + l22);
            }

        }

    virtual double set(double value) {
        int i;

        emili::Solution *s1;
        emili::Solution *s2 = is->generateSolution();
        
        double c1,
               c2 = s2->getSolutionValue();
        double bestcost = c2,
               costsum = 0;

        for (i = 0 ; i < length ; i++) {
            s1 = s2;
            s2 = nei->random(s1);; //is->generateSolution();
            c1 = c2;
            c2 = s2->getSolutionValue();
            costsum += abs(c2 - c1);
            if (abs(c2 - c1) > maxdelta) {
                maxdelta = abs(c2 - c1);
            } else if (abs(c2 - c1) > 0 && abs(c2 - c1) < mindelta) {
                mindelta = abs(c2 - c1);
            }

            delete s1;
        }

        delete s2;

        avgdelta = costsum / length;

        ti = (1 - l11 - l12) * mindelta + l11 * avgdelta + l12 * maxdelta;
        tf = (1 - l21 - l22) * mindelta + l12 * avgdelta + l22 * maxdelta;

        return ti;
    }

    virtual double getMinTemp(void) {
        return tf;
    }

    virtual double getInit_prob(void) {
        return std::exp(maxdelta / ti);
    }

}; // MiseviciusInitTemp

/**
 * misevicius - new sa
 */
class SimplifiedMiseviciusInitTemp: public SAInitTemp {
protected:
    emili::InitialSolution *is;
    emili::Neighborhood *nei;
    int length;
    float l1, l2;
    double mindelta, avgdelta, ti, tf;

public:

    SimplifiedMiseviciusInitTemp(emili::InitialSolution* _is,
                       emili::Neighborhood *_nei,
                       int _length,
                       float _l1,
                       float _l2):
        is(_is),
        nei(_nei),
        length(_length),
        l1(_l1),
        l2(_l2) {
            if (l2 >= l1) {
                std::swap(l1, l2);
            }
        }

    virtual double set(double value) {
        int i;

        emili::Solution *s1;
        emili::Solution *s2 = is->generateSolution();
        
        double c1,
               c2 = s2->getSolutionValue();
        double bestcost = c2,
               costsum = 0;

        for (i = 0 ; i < length ; i++) {
            s1 = s2;
            s2 = nei->random(s1);; //is->generateSolution();
            c1 = c2;
            c2 = s2->getSolutionValue();
            costsum += abs(c2 - c1);
            if (abs(c2 - c1) > 0 && abs(c2 - c1) < mindelta) {
                mindelta = abs(c2 - c1);
            }

            delete s1;
        }

        delete s2;

        avgdelta = costsum / length;

        ti = (1 - l1) * mindelta + l1 * avgdelta;
        tf = (1 - l2) * mindelta + l2 * avgdelta;

        return ti;
    }

    virtual double getMinTemp(void) {
        return tf;
    }

    virtual double getInit_prob(void) {
        return std::exp(avgdelta / ti);
    }

}; // SimplifiedMiseviciusInitTemp



/**
 * see Moscato-Fontanari, Stochastic vs. deterministic update in SA
 *
 * to be considered later...
 */
class BestRatioInitTemp: public SAInitTemp {
protected:
    emili::InitialSolution* is;
    int num_trials;
    float max_temp,
          min_temp,
          target_ratio;

    /*
    emili::InitialSolution* initsol;
    emili::Neighborhood*    nei;
    SAInitTemp*      inittemp;
    SAAcceptance*    acceptance;
    SACooling*       cooling;
    SATempRestart*   temprestart;
    SATermination*   term;
    SATempLength*    templ;
    SAExploration*   explo;
    */

public:
    BestRatioInitTemp(emili::InitialSolution* _is,
                      int   _num_trials,
                      float _max_temp,
                      float _min_temp,
                      float _target_ratio):
        is(_is),
        num_trials(_num_trials),
        max_temp(_max_temp),
        min_temp(_min_temp),
        target_ratio(_target_ratio) { }



}; // BestRatioInitTemp


#endif
