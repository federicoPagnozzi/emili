#ifndef PFSP_PTIMES_H
#define PFSP_PTIMES_H

#include <cstdint>

// Read-only view over the flat int32 processing-time matrix of a PfspInstance.
// Valid for the lifetime of the instance (the buffer is never reallocated
// after allowMatrixMemory).
struct PfspTimesView {
    const int32_t* data;   // (nbJob+1)*(nbMac+1) entries, row-major by job, row 0 and column 0 are zero
    int stride;            // nbMac+1
    int32_t at(int job, int machine) const { return data[job * stride + machine]; }
};

#endif
