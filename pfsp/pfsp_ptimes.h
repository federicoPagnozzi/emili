#ifndef PFSP_PTIMES_H
#define PFSP_PTIMES_H

#include <cstdint>

// Read-only view over the flat int32 processing-time matrix of a PfspInstance.
// A view stays valid until allowMatrixMemory (or an instance loader) is called
// again on that instance, the instance is assigned to, or it is destroyed;
// setTime writes in place and keeps views valid.
struct PfspTimesView {
    const int32_t* data;   // (nbJob+1)*(nbMac+1) entries, row-major by job, row 0 and column 0 are zero
    int stride;            // nbMac+1
    int32_t at(int job, int machine) const { return data[job * stride + machine]; }
};

#endif
