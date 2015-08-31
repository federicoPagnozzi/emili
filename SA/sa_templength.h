#ifndef SA_TEMPLENGTH_H
#define SA_TEMPLENGTH_H


#include "../emilibase.h"


class SATempLength {

protected:
	int length;
	std::string type;

public:
	SATempLength(std::string type, int length):
	    type(type),
	    length(length) { }

	void setLength(int len) {
		length = len;
	}

	int getLength(void) {
		return(length);
	}

	std::string getType() {
	    return type;
	}

}; // SATempLength


/**
 * constant number of iterations at the same temperature
 */
class ConstantTempLength: public SATempLength {

public:
	ConstantTempLength(int length):
	    SATempLength(CONSTTEMPLEN, length) { }

}; // ConstantTempLength



class NeighSizeTempLength: public SATempLength {

protected:
	emili::Neighborhood* neigh;
	float alpha;

public:
	ConstantTempLength(emili::Neighborhood* neigh,
		               float alpha):
	    neigh(neigh),
	    alpha(alpha),
	    SATempLength(NEIGHSIZETEMPLEN,
	    	         std::lrint(alpha * neigh->size())) { }


}; // NeighSizeTempLength


#endif
