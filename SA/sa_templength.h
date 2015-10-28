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
 * Length based on a maximum number of accepted solutions.
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


#endif
