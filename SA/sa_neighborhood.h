#ifndef SA_NEIGHBORHOOD_H
#define SA_NEIGHBORHOOD_H

#include "../emilibase.h"
namespace emili {
namespace sa {


class SANeighborhood: public emili::Neighborhood {

public:
	virtual emili::Solution* random(emili::Solution* currentSolution);

}; // SANeighborhood

}
}
#endif
