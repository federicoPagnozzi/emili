//
//  ARM NEON kernels of the PFSP fast path (the counterpart of sse_functions.h).
//  Include pfsp_simd.h, never this header directly.
//
//  All kernels read the flat int32 processing-time view (row 0 / column 0 are
//  zero, so job id 0 is an all-zero row) and follow the scalar contracts:
//    C[j][m] = max(C[j-1][m], C[j][m-1]) + p(sol[j], m),  C[0][*] = C[*][0] = 0
//    T[j][m] = max(T[j][m+1], T[j+1][m]) + p(sol[j], m),  T[*][nbMac+1] = T[nbJob+1][*] = 0
//  where nbJob is the number of leading jobs sol[1..nbJob] that count. Unlike
//  the SSE kernels they never resize or modify sol / previousMachineEndTime and
//  never write past index nbJob. Arithmetic is int32, widened to long on store.
//
#ifndef PFSP_NEON_FUNCTIONS_H
#define PFSP_NEON_FUNCTIONS_H

#include <arm_neon.h>
#include <algorithm>
#include <cstdint>
#include <vector>

#include "pfsp_ptimes.h"

// [0, v0, v1, v2]: lane k receives lane k-1, lane 0 is zero
static inline int32x4_t pfsp_shift_up1(int32x4_t v)
{
    return vextq_s32(vdupq_n_s32(0), v, 3);
}

// [v0+k, v0+v1+k, v0+v1+v2+k, v0+..+v3+k]
static inline int32x4_t pfsp_prefix_sum(int32x4_t v, int32_t k)
{
    v = vaddq_s32(v, pfsp_shift_up1(v));                     // v0, v0+v1, v1+v2, v2+v3
    v = vaddq_s32(v, vextq_s32(vdupq_n_s32(0), v, 2));       // v0, v0+v1, v0+..+v2, v0+..+v3
    return vaddq_s32(v, vdupq_n_s32(k));
}

// [pt.at(j0,m0), pt.at(j1,m1), pt.at(j2,m2), pt.at(j3,m3)]
static inline int32x4_t pfsp_gather(const PfspTimesView& pt, int j0, int j1, int j2, int j3,
                                    int m0, int m1, int m2, int m3)
{
    int32x4_t v = vdupq_n_s32(pt.at(j0, m0));
    v = vsetq_lane_s32(pt.at(j1, m1), v, 1);
    v = vsetq_lane_s32(pt.at(j2, m2), v, 2);
    v = vsetq_lane_s32(pt.at(j3, m3), v, 3);
    return v;
}

// p(job, m) for any m, 0 outside the machines 1..nbMac (fill / drain lanes)
static inline int32_t pfsp_time_or_zero(const PfspTimesView& pt, int job, int m, int nbMac)
{
    return (m >= 1 && m <= nbMac) ? pt.at(job, m) : 0;
}

/*  Anti-diagonal wavefront over groups of 4 jobs J0..J3 = sol[j..j+3].
 *
 *  At step s (s = 1 .. nbMac+3) lane k computes job Jk on machine m = s-k:
 *
 *        v (after step s)           up = shift_up1(v_prev), lane 0 <- L[s]
 *    [ C0,s  C1,s-1  C2,s-2  C3,s-3 ]   [ L[s]  C0,s-1  C1,s-2  C2,s-3 ]
 *
 *    v = max(v_prev, up) + [p(J0,s), p(J1,s-1), p(J2,s-2), p(J3,s-3)]
 *
 *  v_prev lane k is Ck,m-1 and up lane k is Ck-1,m (lane k-1 one step
 *  earlier), lane 0 takes the carry L[m]: the previous group's last job on
 *  machine m. Lane k only reads lanes < k, so lanes past nbJob (job id 0 in
 *  the remainder group) never reach a stored value.
 *
 *  fill  s = 1..3      : lanes k >= s are before machine 1; they stay 0 by
 *                        induction (v_prev = up = 0, p = 0).
 *  steady s = 4..nbMac : all lanes on a real machine, unchecked gather.
 *  drain s = nbMac+1.. : lanes k < s-nbMac are done; lane k's final value
 *                        Ck,nbMac is captured at step nbMac+k, before later
 *                        steps can touch it.
 *  Lane 3 at step s >= 4 produces C3,s-3, the carry for the next group.
 *  nbMac < 4 needs no special path: the steady loop is empty and fill / drain
 *  overlap, every step using the checked gather.
 *
 *  Store == true writes out[j+k] = Ck,nbMac for the real lanes. Returns the
 *  completion time of sol[nbJob] on machine nbMac. Requires nbJob >= 1.
 */
template <bool Store>
static inline int32_t pfsp_neon_wavefront(const std::vector<int>& sol, long* out,
                                          const PfspTimesView& pt, int nbJob, int nbMac)
{
    std::vector<int32_t> L(nbMac + 1, 0);
    int32_t fin[4] = {0, 0, 0, 0};
    int lanes = 0;

    for (int j = 1; j <= nbJob; j += 4) {
        lanes = std::min(4, nbJob - j + 1);
        const int J0 = sol[j];
        const int J1 = lanes > 1 ? sol[j+1] : 0;
        const int J2 = lanes > 2 ? sol[j+2] : 0;
        const int J3 = lanes > 3 ? sol[j+3] : 0;

        int32x4_t v = vdupq_n_s32(0);

        // one checked step (fill / drain): any s in 1..nbMac+3
        auto checked_step = [&](int s) {
            int32x4_t up = pfsp_shift_up1(v);
            if (s <= nbMac)
                up = vsetq_lane_s32(L[s], up, 0);
            int32x4_t p = vdupq_n_s32(pfsp_time_or_zero(pt, J0, s, nbMac));
            p = vsetq_lane_s32(pfsp_time_or_zero(pt, J1, s - 1, nbMac), p, 1);
            p = vsetq_lane_s32(pfsp_time_or_zero(pt, J2, s - 2, nbMac), p, 2);
            p = vsetq_lane_s32(pfsp_time_or_zero(pt, J3, s - 3, nbMac), p, 3);
            v = vaddq_s32(vmaxq_s32(v, up), p);
            if (s >= 4)                       // s-3 <= nbMac since s <= nbMac+3
                L[s - 3] = vgetq_lane_s32(v, 3);
            if (s >= nbMac) {                 // lane s-nbMac just finished machine nbMac
                int32_t tmp[4];
                vst1q_s32(tmp, v);
                fin[s - nbMac] = tmp[s - nbMac];
            }
        };

        int s = 1;
        for (; s <= 3; ++s)
            checked_step(s);
        for (; s <= nbMac; ++s) {
            // v -> [C0,s C1,s-1 C2,s-2 C3,s-3], all lanes on machines 1..nbMac
            const int32x4_t up = vsetq_lane_s32(L[s], pfsp_shift_up1(v), 0);
            v = vaddq_s32(vmaxq_s32(v, up), pfsp_gather(pt, J0, J1, J2, J3, s, s - 1, s - 2, s - 3));
            L[s - 3] = vgetq_lane_s32(v, 3);
        }
        if (nbMac >= 4)                       // step nbMac ran in the steady loop
            fin[0] = vgetq_lane_s32(v, 0);
        for (; s <= nbMac + 3; ++s)
            checked_step(s);

        if (Store)
            for (int k = 0; k < lanes; ++k)
                out[j + k] = fin[k];
    }
    return fin[lanes - 1];
}

/** previousMachineEndTime[j] = C[j][nbMac] for j = 1..nbJob, and
 *  previousMachineEndTime[0] = 0 as in the scalar kernel. */
inline void computePartialMakespans(std::vector<int>& sol, std::vector< long >& previousMachineEndTime, const PfspTimesView& pt, int nbJob, int nbMac)
{
    previousMachineEndTime[0] = 0;
    if (nbJob >= 1)
        pfsp_neon_wavefront<true>(sol, previousMachineEndTime.data(), pt, nbJob, nbMac);
}

/** C[size][nbMac], 0 for size == 0 (PfspInstance::computeMS on NEON). */
inline long pfsp_neon_makespan(std::vector<int>& sol, const PfspTimesView& pt, int size, int nbMac)
{
    if (size <= 0)
        return 0;
    return pfsp_neon_wavefront<false>(sol, nullptr, pt, size, nbMac);
}

// placeholder: replaced in Task 6
inline void computeHEAD(std::vector<int> &sol, std::vector< std::vector < int > >& head, const PfspTimesView& pt, int nbJob, int nbMac)
{
    for (int j = 1; j <= nbJob; ++j)
        for (int m = 1; m <= nbMac; ++m) {
            const int prevJob = j > 1 ? head[m][j-1] : 0;
            const int prevMac = m > 1 ? head[m-1][j] : 0;
            head[m][j] = std::max(prevJob, prevMac) + pt.at(sol[j], m);
        }
}

// placeholder: replaced in Task 7
inline void computeTAIL(std::vector<int> &sol, std::vector< std::vector < int > >& tail, const PfspTimesView& pt, int nbJob, int nbMac)
{
    for (int j = nbJob; j >= 1; --j)
        for (int m = nbMac; m >= 1; --m) {
            const int nextJob = j < nbJob ? tail[m][j+1] : 0;
            const int nextMac = m < nbMac ? tail[m+1][j] : 0;
            tail[m][j] = std::max(nextJob, nextMac) + pt.at(sol[j], m);
        }
}

// placeholder: replaced in Task 7
inline void computeHEADandTAIL(std::vector<int> &sol, std::vector< std::vector < int > >& head, std::vector< std::vector < int > >& tail, const PfspTimesView& pt, int nbJob, int nbMac)
{
    computeHEAD(sol, head, pt, nbJob, nbMac);
    computeTAIL(sol, tail, pt, nbJob, nbMac);
}

#endif // PFSP_NEON_FUNCTIONS_H
