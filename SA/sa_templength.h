#ifndef SA_TEMPLENGTH_H
#define SA_TEMPLENGTH_H

#include "sa_common.h"

#include "sa_constants.h"
#include "../emilibase.h"


class SATempLength {

protected:
    int length;
    std::string type;
    SAStatus* status;

public:
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
 * Length dependant on the size of the neighborhood:
 * alpha * |neigh|
 */
class NeighSizeTempLength: public SATempLength {

protected:
    emili::Neighborhood* neigh;
    float alpha;

public:
    NeighSizeTempLength(emili::Neighborhood* neigh,
                       float alpha):
        neigh(neigh),
        alpha(alpha),
        SATempLength(NEIGHSIZETEMPLEN,
                     (int)std::lrint(alpha * neigh->size())) { }

    bool isCoolingTime(int counter) {
        if (counter >= length)
            return true;
        return false;
    }

}; // NeighSizeTempLength


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
    BurkardRendlNeighSizeTempLength(emili::Neighborhood* neigh,
                       float alpha):
        neigh(neigh),
        alpha(alpha),
        SATempLength(BRNEIGHSIZETEMPLEN,
                     (int)alpha * std::ceil(std::sqrt(neigh->size())) * std::ceil(std::sqrt(neigh->size()))) { }

    bool isCoolingTime(int counter) {
        if (counter >= length)
            return true;
        return false;
    }

}; // BurkardRendlNeighSizeTempLength


class BurkardRendlGeomTempLength: public SATempLength {

protected:
    float c;

public:
    BurkardRendlGeomTempLength(emili::Neighborhood* neigh, float _c):
        c(_c),
        SATempLength(BRGEOMTEMPLEN,
            (int)0.5 * std::ceil(std::sqrt(neigh->size())) * std::ceil(std::sqrt(neigh->size()))) { }

    bool isCoolingTime(int counter) {
        if (counter >= length) {
            length  = length * c;
            return true;
        }
        return false;
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
    NeighSizeCappedMaxAcceptedTempLength(emili::Neighborhood* neigh,
                       float alpha,
                       long _cap):
        neigh(neigh),
        alpha(alpha),
        cap(_cap),
        SATempLength(NEIGHCAPPEDMAXACCEPTEDTEMPLEN,
                     (int)std::lrint(alpha * neigh->size())) { }

    bool isCoolingTime(int counter) {
        if (status->curr_accepted >= length || counter >= cap) {
            status->curr_accepted = 0;
            return true;
        }
        return false;
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
            length  = length / c;
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
            length  = c / log(status->step);
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


#endif
