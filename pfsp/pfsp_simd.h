#ifndef PFSP_SIMD_H
#define PFSP_SIMD_H

// Compile-time SIMD dispatch for the PFSP fast path.
// Include this header, never sse_functions.h / neon_functions.h directly.
// Inert when ENABLE_SIMD is not defined (scalar build).
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
