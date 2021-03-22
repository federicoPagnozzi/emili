#ifndef SA_TEMPLENGTH_H
#define SA_TEMPLENGTH_H

#include "sa_common.h"

#include "sa_constants.h"
#include "../emilibase.h"
namespace emili {
namespace sa{

class SATempLength {

protected:
    long length;
    std::string type;
    SAStatus* status;

public:
    SATempLength(std::string type):
        type(type){}

    SATempLength(std::string type, int length):
        type(type),
        length(length) { }

    void set_status(SAStatus* _status) {
        status = _status;
    }

    void setLength(int len) {
        length = len;
    }

    int getLength(void) {
        return(length);
    }

    std::string getType() {
        return type;
    }

    virtual bool isCoolingTime(int counter)=0;

    virtual void setNeighborhoodSize(int neighborhood_size) {}

}; // SATempLength


/**
 * constant number of iterations at the same temperature
 */
class ConstantTempLength: public SATempLength {

public:
    ConstantTempLength(int length):
        SATempLength(CONSTTEMPLEN, length) { }

    bool isCoolingTime(int counter) {
        if (counter >= length)
            return true;
        return false;
    }

}; // ConstantTempLength


/**
 ** random number of iterations at the same temperature
 **/
class RandomTempLength: public SATempLength {
protected:
  int lblength, ublength;

public:
    RandomTempLength(int lblength, int ublength):
        lblength(lblength),
        ublength(ublength),
        SATempLength(RANDOMTEMPLEN, 0) {
           if (lblength > ublength) {
               int tmp = lblength;
               lblength = ublength;
               ublength = tmp;
           }
           length = lblength + int(emili::generateRealRandomNumber() * std::abs(ublength - lblength));
        }

    bool isCoolingTime(int counter) {
        if (counter >= length) {
            length = lblength + int(emili::generateRealRandomNumber() * std::abs(ublength - lblength));
            return true;
        }
        return false;
    }

}; // RandomTempLength


/**
 * Length dependant on the size of the neighborhood:
 * alpha * |neigh|
 */
class NeighSizeTempLength: public SATempLength {

protected:
    emili::Neighborhood* neigh;
    float alpha;

public:
    NeighSizeTempLength(float alpha):
        alpha(alpha),
        SATempLength(NEIGHSIZETEMPLEN) { }

    NeighSizeTempLength(emili::Neighborhood* neigh,
                       float alpha):
        neigh(neigh),
        alpha(alpha),
        SATempLength(NEIGHSIZETEMPLEN,
                     std::lrint(alpha * neigh->size())) { }

    bool isCoolingTime(int counter) {
        if (counter >= length)
            return true;
        return false;
    }

    virtual void setNeighborhoodSize(int neighborhood_size)
    {
        this->length = std::lrint(alpha * neighborhood_size);
    }

}; // NeighSizeTempLength


/**
 * Length dependant on the size of the problem:
 * alpha * |n|
 * PROBSIZETEMPLEN
 */
class ProblemSizeTempLength: public SATempLength {

protected:
    emili::Problem *prob;
    float alpha;

public:
    ProblemSizeTempLength(emili::Problem* prob,
                       float alpha):
        prob(prob),
        alpha(alpha),
        SATempLength(PROBSIZETEMPLEN,
                     std::lrint(alpha * prob->problemSize())) {
            // std::cout << "initial temp length: " << length << std::endl << std::endl;
        }

    bool isCoolingTime(int counter) {
        if (counter >= length) {
            return true;
        }
        return false;
    }

}; // ProblemSizeTempLength

/**
 * Length dependant on the size of the problem:
 * alpha * |n|^2
 * SQUAREDPROBSIZETEMPLEN
 */
class SquaredProblemSizeTempLength: public SATempLength {

protected:
    emili::Problem* prob;
    float alpha;

public:
    SquaredProblemSizeTempLength(emili::Problem* prob,
                       float alpha):
        prob(prob),
        alpha(alpha),
        SATempLength(SQUAREDPROBSIZETEMPLEN,
                     std::lrint(alpha * prob->problemSize())) { }

    bool isCoolingTime(int counter) {
        if (counter >= length) {
            return true;
        }
        return false;
    }

}; // SquaredProblemSizeTempLength


/**
 * Length dependant on the size of the instance:
 * Burkard-Rendl (From COnnolly Paper)
 * 1/2 * |n|^2
 */
class BurkardRendlNeighSizeTempLength: public SATempLength {

protected:
    emili::Neighborhood* neigh;
    float alpha;

public:
    BurkardRendlNeighSizeTempLength(float alpha):
        alpha(alpha),
        SATempLength(BRNEIGHSIZETEMPLEN) {
            // std::cout << "initial temp length: " << length << std::endl << std::endl;
        }

    BurkardRendlNeighSizeTempLength(emili::Neighborhood* neigh,
                       float alpha):
        neigh(neigh),
        alpha(alpha),
        SATempLength(BRNEIGHSIZETEMPLEN,
                     std::lrint(alpha * std::ceil(std::sqrt(neigh->size())) * std::ceil(std::sqrt(neigh->size())))) {
            // std::cout << "initial temp length: " << length << std::endl << std::endl;
        }

    bool isCoolingTime(int counter) {
        if (counter >= length)
            return true;
        return false;
    }

    virtual void setNeighborhoodSize(int neighborhood_size)
    {
        length = std::lrint(alpha * std::ceil(std::sqrt(neighborhood_size)) * std::ceil(std::sqrt(neighborhood_size)));
    }

}; // BurkardRendlNeighSizeTempLength


class BurkardRendlGeomTempLength: public SATempLength {

protected:
    float c;

public:
    BurkardRendlGeomTempLength(float _c):
        c(_c),
        SATempLength(BRGEOMTEMPLEN) { }

    BurkardRendlGeomTempLength(emili::Neighborhood* neigh, float _c):
        c(_c),
        SATempLength(BRGEOMTEMPLEN,
            std::lrint(0.5 * std::ceil(std::sqrt(neigh->size())) * std::ceil(std::sqrt(neigh->size())))) {
            // std::cout << std::ceil(std::sqrt(neigh->size())) << " " << 0.5 * std::ceil(std::sqrt(neigh->size())) * std::ceil(std::sqrt(neigh->size())) << " " << std::lrint(0.5 * std::ceil(std::sqrt(neigh->size())) * std::ceil(std::sqrt(neigh->size()))) << std::endl;
            // std::cout << "initial temp length: " << length << std::endl << std::endl;
        }

    bool isCoolingTime(int counter) {
        if (counter >= length) {
            length  = length * c;
            return true;
        }
        return false;
    }

    virtual void setNeighborhoodSize(int neighborhood_size)
    {
        length = std::lrint(0.5 * std::ceil(std::sqrt(neighborhood_size)) * std::ceil(std::sqrt(neighborhood_size)));
    }


}; // BurkardRendlGeomTempLength


/**
 * Length based on a maximum number of accepted solutions. (Schaerf? Jajodia?)
 */
class MaxAcceptedTempLength: public SATempLength {

public:
    MaxAcceptedTempLength(int length):
        SATempLength(MAXACCEPTEDTEMPLEN, length) { }

    bool isCoolingTime(int counter) {
        if (status->curr_accepted >= length) {
            status->curr_accepted = 0;
            return true;
        }
        return false;
    }

}; // MaxAcceptedTempLength


/**
 * Length based on a maximum number of accepted solutions OR max length (Jajodia - class)
 */
class CappedMaxAcceptedTempLength: public SATempLength {

protected:
    long cap;

public:
    CappedMaxAcceptedTempLength(int length, long _cap):
        cap(_cap),
        SATempLength(CAPPEDMAXACCEPTEDTEMPLEN, length) { }

    bool isCoolingTime(int counter) {
        if (status->curr_accepted >= length || counter > cap) {
            status->curr_accepted = 0;
            return true;
        }
        return false;
    }

}; // CappedMaxAcceptedTempLength


/**
 * Length dependant on the size of the neighborhood:
 * alpha * |neigh|
 * (Jajodia - class)
 */
class NeighSizeCappedMaxAcceptedTempLength: public SATempLength {

protected:
    emili::Neighborhood* neigh;
    float alpha;
    long cap;

public:
    NeighSizeCappedMaxAcceptedTempLength(float alpha,long _cap):
        alpha(alpha),
        cap(_cap),
        SATempLength(NEIGHCAPPEDMAXACCEPTEDTEMPLEN) { }

    NeighSizeCappedMaxAcceptedTempLength(emili::Neighborhood* neigh,
                       float alpha,
                       long _cap):
        neigh(neigh),
        alpha(alpha),
        cap(_cap * std::lrint(alpha * neigh->size())),
        SATempLength(NEIGHCAPPEDMAXACCEPTEDTEMPLEN,
                     std::lrint(alpha * neigh->size())) { }

    bool isCoolingTime(int counter) {
        if (status->curr_accepted >= length || counter >= cap) {
            status->curr_accepted = 0;
            return true;
        }
        return false;
    }

    virtual void setNeighborhoodSize(int neighborhood_size)
    {
        length = std::lrint(alpha * neighborhood_size);
        cap = cap * length;
    }

}; // NeighSizeCappedMaxAcceptedTempLength


/**
 * @article{rose1990temperature,
  title={Temperature measurement and equilibrium dynamics of simulated annealing placements},
  author={Rose, Jonathan and Klebsch, Wolfgang and Wolf, Jurgen},
  journal={Computer-Aided Design of Integrated Circuits and Systems, IEEE Transactions on},
  volume={9},
  number={3},
  pages={253--259},
  year={1990},
  publisher={IEEE}
}

Simulated annealing for manufacturing systems layout design

    Abdelghani Souilah

 */
class ArithmeticTempLength: public SATempLength {

protected:
    int inc;

public:
    ArithmeticTempLength(int length, int _inc):
        inc(_inc),
        SATempLength(ARITMTEMPLEN, length) { }

    bool isCoolingTime(int counter) {
        if (counter >= length) {
            length += inc;
            return true;
        }
        return false;
    }

}; // ArithmeticTempLength


class GeomTempLength: public SATempLength {

protected:
    float c;

public:
    GeomTempLength(int length, float _c):
        c(_c),
        SATempLength(GEOMTEMPLEN, length) { }

    bool isCoolingTime(int counter) {
        if (counter >= length) {
            length  = length * c;
            return true;
        }
        return false;
    }

}; // GeomTempLength


class LogTempLength: public SATempLength {

protected:
    int c;

public:
    LogTempLength(int length, int _c):
        c(_c),
        SATempLength(LOGTEMPLEN, length) { }

    bool isCoolingTime(int counter) {
        if (counter >= length) {
            length  = c / log(status->temp);
            return true;
        }
        return false;
    }

}; // LogTempLength


class ExpTempLength: public SATempLength {

protected:
    float c;

public:
    ExpTempLength(int length, float _c):
        c(_c),
        SATempLength(EXPTEMPLEN, length) { }

    bool isCoolingTime(int counter) {
        if (counter >= length) {
            length  = std::pow(length, 1/c);
            return true;
        }
        return false;
    }

}; // ExpTempLength

/**
 * No Temp Length option.
 */
class NoTempLength: public SATempLength {

public:
    NoTempLength(void):
        SATempLength(NOTEMPLEN, 0) { }

    bool isCoolingTime(int counter) {
        return false;
    }

}; // NoTempLength

}
}


#endif
