#ifndef PFSP_SIMD_H
#define PFSP_SIMD_H

// Compile-time SIMD dispatch for the PFSP fast path.
// Include this header, never sse_functions.h / neon_functions.h directly.
// Inert when ENABLE_SIMD is not defined (scalar build).
//
// The kernels (same signatures on SSE and NEON; pt is the int32 PfspTimesView):
//   computePartialMakespans(sol, previousMachineEndTime, pt, nbJob, nbMac)
//   computeHEAD(sol, head, pt, nbJob, nbMac)
//   computeTAIL(sol, tail, pt, nbJob, nbMac)
//   computeHEADandTAIL(sol, head, tail, pt, nbJob, nbMac)
// plus the makespan kernel behind PfspInstance::computeMS: on NEON
// pfsp_neon_makespan (neon_functions.h); on SSE computeMS is defined in
// pfspinstance.cpp and still reads the long matrix.
// Limits and contracts are documented at the top of each ISA header.
#ifdef ENABLE_SIMD
#  if defined(__SSE__)
#    define PFSP_SIMD_NAME "SSE"
#    include "sse_functions.h"
#  elif defined(__ARM_NEON)
#    define PFSP_SIMD_NAME "NEON"
#    include "neon_functions.h"
#  else
#    error "ENABLE_SIMD set but neither __SSE__ nor __ARM_NEON is defined"
#  endif
#endif

#endif
