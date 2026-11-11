// simd_parity: checks the SIMD data structures / kernels against the scalar code.
// usage: simd_parity <instance> [seed] [rounds]
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <memory>
#include <numeric>
#include <random>
#include <vector>

#include "../../pfsp/pfspinstance.h"
#include "../../pfsp/pfsp_ptimes.h"
#include "../../pfsp/pfsp_simd.h"

// Strict bound checks (sentinel cells never written, vector sizes unchanged,
// previousMachineEndTime[0] == 0) are compiled for NEON only. The SSE
// computePartialMakespans pads sol / previousMachineEndTime with push_back and
// pop_back and stores whole 4-job groups, so for nbJob % 4 != 0 it writes
// previousMachineEndTime[nbJob+1 .. nbJob+3] and leaves index 0 alone; the SSE
// head/tail kernels have not been audited against the sentinel contract and
// are assumed to read row/column 0 of head/tail as zeros. Value comparisons on
// the contract cells run on every SIMD path.
#if defined(ENABLE_SIMD) && defined(__ARM_NEON)
#define PFSP_PARITY_STRICT 1
#else
#define PFSP_PARITY_STRICT 0
#endif

// The SSE computeMS does not support size 0; scalar and NEON return 0.
#if defined(ENABLE_SIMD) && defined(__SSE__)
#define PFSP_PARITY_MS_EMPTY 0
#else
#define PFSP_PARITY_MS_EMPTY 1
#endif

// Failure messages print matrix cells as cell=[machine][job] (head / tail /
// view layout: row = machine, column = job), vectors as name[index].

static int check_view(PfspInstance& inst)
{
    const PfspTimesView view = inst.getTimesView();
    const int nj = inst.getNbJob();
    const int nm = inst.getNbMac();

    for (int m = 0; m <= nm; ++m) {
        if (view.at(0, m) != 0) {
            std::printf("kernel=view size=%d cell=[%d][0] expected=0 got=%d\n", nj, m, (int)view.at(0, m));
            return 1;
        }
    }
    for (int j = 1; j <= nj; ++j) {
        if (view.at(j, 0) != 0) {
            std::printf("kernel=view size=%d cell=[0][%d] expected=0 got=%d\n", nj, j, (int)view.at(j, 0));
            return 1;
        }
        for (int m = 1; m <= nm; ++m) {
            const int expected = (int)inst.getTime(j, m);
            if (view.at(j, m) != expected) {
                std::printf("kernel=view size=%d cell=[%d][%d] expected=%d got=%d\n", nj, m, j, expected, (int)view.at(j, m));
                return 1;
            }
        }
    }

    PfspInstance copy(inst);
    const PfspTimesView cview = copy.getTimesView();
    if (cview.data == view.data) {
        std::printf("kernel=view-copy size=%d cell=[0][0] expected=distinct-buffer got=shared-buffer\n", nj);
        return 1;
    }
    for (int j = 0; j <= nj; ++j) {
        for (int m = 0; m <= nm; ++m) {
            if (cview.at(j, m) != view.at(j, m)) {
                std::printf("kernel=view-copy size=%d cell=[%d][%d] expected=%d got=%d\n", nj, m, j, (int)view.at(j, m), (int)cview.at(j, m));
                return 1;
            }
        }
    }

    // setTime writes in place: a view taken earlier must observe it
    const long saved = inst.getTime(1, 1);
    inst.setTime(1, 1, saved + 7);
    if (view.at(1, 1) != (int32_t)(saved + 7)) {
        std::printf("kernel=view-settime size=%d cell=[1][1] expected=%d got=%d\n", nj, (int)(saved + 7), (int)view.at(1, 1));
        return 1;
    }
    inst.setTime(1, 1, saved);

    std::printf("view: OK\n");
    return 0;
}

// Reference completion times: C[j][m] = max(C[j-1][m], C[j][m-1]) + p(sol[j], m),
// with C[0][*] = C[*][0] = 0.
static std::vector<std::vector<long>> ref_completion(PfspInstance& inst, const std::vector<int>& sol, int n)
{
    const PfspTimesView pt = inst.getTimesView();
    const int nm = inst.getNbMac();
    std::vector<std::vector<long>> C(n + 1, std::vector<long>(nm + 1, 0));
    for (int j = 1; j <= n; ++j)
        for (int m = 1; m <= nm; ++m)
            C[j][m] = std::max(C[j-1][m], C[j][m-1]) + pt.at(sol[j], m);
    return C;
}

// Identity, the reversed permutation (n, n-1, ..., 1) and `rounds` seeded
// shuffles; sol[0] = 0 is the unused slot. NEON only: identity and reversed
// again with sol[0] = -7. The NEON kernels never read sol[0] by construction
// (Dir::job maps 1..nbJob onto 1..nbJob), so a read would index the view with
// a negative job id; the SSE kernels were not audited for this.
// Note: std::shuffle with mt19937 is not guaranteed to give the same
// permutation across standard libraries, so a failing seed may not reproduce
// on another platform.
static std::vector<std::vector<int>> test_permutations(int nj, unsigned seed, int rounds)
{
    std::vector<std::vector<int>> perms;
    std::vector<int> sol(nj + 1);
    std::iota(sol.begin(), sol.end(), 0);
    perms.push_back(sol);
    std::vector<int> rev = sol;
    std::reverse(rev.begin() + 1, rev.end());
    perms.push_back(rev);
#if defined(ENABLE_SIMD) && defined(__ARM_NEON)
    sol[0] = -7;
    rev[0] = -7;
    perms.push_back(sol);
    perms.push_back(rev);
    sol[0] = 0;
#endif
    std::mt19937 rng(seed);
    for (int r = 0; r < rounds; ++r) {
        std::shuffle(sol.begin() + 1, sol.end(), rng);
        perms.push_back(sol);
    }
    return perms;
}

static int check_sol_unchanged(const char* kernel, const std::vector<int>& sol, const std::vector<int>& orig, int n)
{
    if (sol.size() != orig.size()) {
        std::printf("kernel=%s size=%d cell=sol.size expected=%zu got=%zu\n", kernel, n, orig.size(), sol.size());
        return 1;
    }
    for (size_t i = 0; i < sol.size(); ++i) {
        if (sol[i] != orig[i]) {
            std::printf("kernel=%s size=%d cell=sol[%zu] expected=%d got=%d\n", kernel, n, i, orig[i], sol[i]);
            return 1;
        }
    }
    return 0;
}

// computeMS is a PfspInstance member on every build path (scalar, SSE, NEON).
static int check_makespan(PfspInstance& inst, const std::vector<std::vector<int>>& perms)
{
    const int nj = inst.getNbJob();
    const int nm = inst.getNbMac();
    for (const std::vector<int>& perm : perms) {
        const std::vector<std::vector<long>> C = ref_completion(inst, perm, nj);
        std::vector<int> sol = perm;
        const long full = inst.computeMS(sol);
        if (full != C[nj][nm]) {
            std::printf("kernel=makespan(sol) size=%d nbMac=%d expected=%ld got=%ld\n", nj, nm, C[nj][nm], full);
            return 1;
        }
        if (check_sol_unchanged("makespan(sol)", sol, perm, nj))
            return 1;
        for (int n = 1; n <= nj; ++n) {
            const long got = inst.computeMS(sol, n);
            if (got != C[n][nm]) {
                std::printf("kernel=makespan(sol,n) size=%d nbMac=%d expected=%ld got=%ld\n", n, nm, C[n][nm], got);
                return 1;
            }
            if (check_sol_unchanged("makespan(sol,n)", sol, perm, n))
                return 1;
        }
#if PFSP_PARITY_MS_EMPTY
        const long empty = inst.computeMS(sol, 0);
        if (empty != 0) {
            std::printf("kernel=makespan(sol,0) size=0 nbMac=%d expected=0 got=%ld\n", nm, empty);
            return 1;
        }
#endif
    }
    return 0;
}

#ifdef ENABLE_SIMD
// Reference tails: T[j][m] = max(T[j][m+1], T[j+1][m]) + p(sol[j], m), with
// T[*][nbMac+1] = T[n+1][*] = 0. Indexed T[j][m].
static std::vector<std::vector<long>> ref_tail(PfspInstance& inst, const std::vector<int>& sol, int n)
{
    const PfspTimesView pt = inst.getTimesView();
    const int nm = inst.getNbMac();
    std::vector<std::vector<long>> T(n + 2, std::vector<long>(nm + 2, 0));
    for (int j = n; j >= 1; --j)
        for (int m = nm; m >= 1; --m)
            T[j][m] = std::max(T[j][m+1], T[j+1][m]) + pt.at(sol[j], m);
    return T;
}

static int check_partial_makespans(PfspInstance& inst, const std::vector<std::vector<int>>& perms)
{
    const int nj = inst.getNbJob();
    const int nm = inst.getNbMac();
    for (const std::vector<int>& perm : perms) {
        const std::vector<std::vector<long>> C = ref_completion(inst, perm, nj);
        for (int n = 1; n <= nj; ++n) {
            // the full permutation is in sol, the kernel is told only the first n jobs count
            std::vector<int> sol = perm;
            std::vector<long> pmet(nj + 1, -1);
            computePartialMakespans(sol, pmet, inst.getTimesView(), n, nm);
            for (int j = 1; j <= n; ++j) {
                if (pmet[j] != C[j][nm]) {
                    std::printf("kernel=partial_makespans size=%d nbMac=%d cell=[%d][%d] expected=%ld got=%ld\n",
                                n, nm, nm, j, C[j][nm], pmet[j]);   // [machine nbMac][job j] = pmet[j]
                    return 1;
                }
            }
            if (check_sol_unchanged("partial_makespans", sol, perm, n))
                return 1;
#if PFSP_PARITY_STRICT
            if (pmet.size() != (size_t)(nj + 1)) {
                std::printf("kernel=partial_makespans size=%d cell=pmet.size expected=%d got=%zu\n", n, nj + 1, pmet.size());
                return 1;
            }
            // the scalar kernel sets previousMachineEndTime[0] = 0
            if (pmet[0] != 0) {
                std::printf("kernel=partial_makespans size=%d cell=pmet[0] expected=0 got=%ld\n", n, pmet[0]);
                return 1;
            }
            for (int j = n + 1; j <= nj; ++j) {
                if (pmet[j] != -1) {
                    std::printf("kernel=partial_makespans size=%d cell=pmet[%d] expected=-1(unwritten) got=%ld\n", n, j, pmet[j]);
                    return 1;
                }
            }
#endif
        }
    }
    return 0;
}

// head[m][j] == C[j][m] (want_head) / tail[m][j] == T[j][m] (want_tail) for
// m = 1..nbMac, j = 1..n; with strict checks every other cell keeps the fill.
static int check_matrix(const char* kernel, const std::vector<std::vector<int>>& M,
                        const std::vector<std::vector<long>>& ref, int n, int nm, int fill)
{
    for (size_t m = 0; m < M.size(); ++m) {
        for (size_t j = 0; j < M[m].size(); ++j) {
            const bool contract = m >= 1 && (int)m <= nm && j >= 1 && (int)j <= n;
            if (contract) {
                const long expected = ref[j][m];   // ref is indexed [job][machine]
                if (M[m][j] != expected) {
                    std::printf("kernel=%s size=%d nbMac=%d cell=[%zu][%zu] expected=%ld got=%d\n",
                                kernel, n, nm, m, j, expected, M[m][j]);
                    return 1;
                }
            } else if (PFSP_PARITY_STRICT && M[m][j] != fill) {
                std::printf("kernel=%s size=%d nbMac=%d cell=[%zu][%zu] expected=%d(unwritten) got=%d\n",
                            kernel, n, nm, m, j, fill, M[m][j]);
                return 1;
            }
        }
    }
    return 0;
}

enum HeadTailKind { HT_HEAD, HT_TAIL, HT_BOTH };

static int check_head_tail(PfspInstance& inst, const std::vector<std::vector<int>>& perms, HeadTailKind kind)
{
    const int nj = inst.getNbJob();
    const int nm = inst.getNbMac();
    const int fill = PFSP_PARITY_STRICT ? -1 : 0;
    const char* name = kind == HT_HEAD ? "head" : kind == HT_TAIL ? "tail" : "head_and_tail";
    for (const std::vector<int>& perm : perms) {
        for (int n = 1; n <= nj; ++n) {
            std::vector<int> sol = perm;
            std::vector<std::vector<int>> head(nm + 2, std::vector<int>(nj + 2, fill));
            std::vector<std::vector<int>> tail(nm + 2, std::vector<int>(nj + 2, fill));
            if (kind == HT_HEAD)
                computeHEAD(sol, head, inst.getTimesView(), n, nm);
            else if (kind == HT_TAIL)
                computeTAIL(sol, tail, inst.getTimesView(), n, nm);
            else
                computeHEADandTAIL(sol, head, tail, inst.getTimesView(), n, nm);
            if (kind != HT_TAIL && check_matrix(name, head, ref_completion(inst, perm, n), n, nm, fill))
                return 1;
            if (kind != HT_HEAD && check_matrix(name, tail, ref_tail(inst, perm, n), n, nm, fill))
                return 1;
            if (check_sol_unchanged(name, sol, perm, n))
                return 1;
        }
    }
    return 0;
}

// Taillard identity, a cross-check of head and tail together: the longest path
// through cell (j, m) is head[m][j] + tail[m][j] - p(sol[j], m), so for every
// position j = 1..n its maximum over m is the makespan C[n][nbMac]. Checked on
// computeHEADandTAIL and on computeHEAD + computeTAIL. Today these are the same
// two calls on NEON; both variants are kept to guard a future fused kernel.
// The matrices are zero-filled here; stray-write detection lives in
// check_head_tail (sentinel fill).
static int check_head_tail_identity(PfspInstance& inst, const std::vector<std::vector<int>>& perms)
{
    const int nj = inst.getNbJob();
    const int nm = inst.getNbMac();
    const PfspTimesView pt = inst.getTimesView();
    for (const std::vector<int>& perm : perms) {
        const std::vector<std::vector<long>> C = ref_completion(inst, perm, nj);
        for (int n = 1; n <= nj; ++n) {
            for (int fused = 0; fused <= 1; ++fused) {
                std::vector<int> sol = perm;
                std::vector<std::vector<int>> head(nm + 2, std::vector<int>(nj + 2, 0));
                std::vector<std::vector<int>> tail(nm + 2, std::vector<int>(nj + 2, 0));
                if (fused) {
                    computeHEADandTAIL(sol, head, tail, pt, n, nm);
                } else {
                    computeHEAD(sol, head, pt, n, nm);
                    computeTAIL(sol, tail, pt, n, nm);
                }
                for (int j = 1; j <= n; ++j) {
                    long longest = 0;
                    for (int m = 1; m <= nm; ++m)
                        longest = std::max(longest, (long)head[m][j] + tail[m][j] - pt.at(perm[j], m));
                    if (longest != C[n][nm]) {
                        std::printf("kernel=head_tail_identity(%s) size=%d nbMac=%d job=%d expected=%ld got=%ld\n",
                                    fused ? "head_and_tail" : "head+tail", n, nm, j, C[n][nm], longest);
                        return 1;
                    }
                }
            }
        }
    }
    return 0;
}

#if defined(__ARM_NEON)
static int check_neon_helpers()
{
    int32_t out[4];

    // 3 jobs x 4 machines, stride 5
    std::vector<int32_t> buf(4 * 5);
    for (size_t i = 0; i < buf.size(); ++i)
        buf[i] = (int32_t)(i * 10);
    const PfspTimesView pt{buf.data(), 5};
    vst1q_s32(out, pfsp_gather(pt, 1, 2, 3, 0, 4, 3, 2, 1));
    const int32_t ga[4] = {pt.at(1, 4), pt.at(2, 3), pt.at(3, 2), pt.at(0, 1)};
    for (int k = 0; k < 4; ++k)
        if (out[k] != ga[k]) {
            std::printf("kernel=pfsp_gather lane=%d expected=%d got=%d\n", k, ga[k], out[k]);
            return 1;
        }
    std::printf("neon_helpers: OK\n");
    return 0;
}
#endif
#endif

static std::unique_ptr<PfspInstance> make_const_instance(int nj, int nm, long p)
{
    std::unique_ptr<PfspInstance> s(new PfspInstance());
    s->setSilence(true);
    s->setNbJob(nj);
    s->setNbMac(nm);
    s->allowMatrixMemory(nj, nm);
    for (int j = 1; j <= nj; ++j)
        for (int m = 1; m <= nm; ++m)
            s->setTime(j, m, p);
    return s;
}

int main(int argc, char** argv)
{
    if (argc < 2) {
        std::fprintf(stderr, "usage: %s <instance> [seed] [rounds]\n", argv[0]);
        return 2;
    }
    PfspInstance inst;
    inst.setSilence(true);
    // char* loader: handles the Taillard-style files (e.g. DD_Ta055.txt) and
    // delegates to the std::string loader for the 12345 format
    if (!inst.readDataFromFile(argv[1])) {
        std::fprintf(stderr, "cannot read %s\n", argv[1]);
        return 2;
    }
    if (inst.getNbJob() < 1 || inst.getNbMac() < 1) {
        std::fprintf(stderr, "empty instance %s\n", argv[1]);
        return 2;
    }
    const unsigned seed = argc > 2 ? (unsigned)std::strtoul(argv[2], nullptr, 10) : 42u;
    const int rounds = argc > 3 ? std::atoi(argv[3]) : 20;

    int fail = check_view(inst);

    // The file instance plus synthetic 9-job instances with few machines and a
    // 13 x 8 one (times 1..99 from the seeded RNG). nbMac 1..3 is NEON / scalar only: the
    // SSE kernels index machines up to 4 unconditionally.
    std::vector<std::unique_ptr<PfspInstance>> synth;
#if defined(ENABLE_SIMD) && defined(__SSE__)
    const int min_synth_mac = 4;
#else
    const int min_synth_mac = 1;
#endif
    {
        std::mt19937 rng(seed);
        std::uniform_int_distribution<int> ptime(1, 99);
        for (int nm = min_synth_mac; nm <= 5; ++nm) {
            synth.emplace_back(new PfspInstance());
            PfspInstance& s = *synth.back();
            s.setSilence(true);
            s.setNbJob(9);
            s.setNbMac(nm);
            s.allowMatrixMemory(9, nm);
            for (int j = 1; j <= 9; ++j)
                for (int m = 1; m <= nm; ++m)
                    s.setTime(j, m, ptime(rng));
        }
        // 13 jobs x 8 machines: steady wavefront steps (nbMac >= 4) with every
        // remainder-group size (n = 1..13 covers n % 4 = 0..3, incl. 1,2,3,4,5,8,13)
        synth.emplace_back(new PfspInstance());
        PfspInstance& s = *synth.back();
        s.setSilence(true);
        s.setNbJob(13);
        s.setNbMac(8);
        s.allowMatrixMemory(13, 8);
        for (int j = 1; j <= 13; ++j)
            for (int m = 1; m <= 8; ++m)
                s.setTime(j, m, ptime(rng));
#if defined(ENABLE_SIMD) && defined(__ARM_NEON)
        // 6 jobs x (PFSP_NEON_STACK_MACHINES + 3) machines: the NEON carry's
        // heap fallback (the stack buffer covers nbMac <= PFSP_NEON_STACK_MACHINES)
        const int big_mac = PFSP_NEON_STACK_MACHINES + 3;
        synth.emplace_back(new PfspInstance());
        PfspInstance& b = *synth.back();
        b.setSilence(true);
        b.setNbJob(6);
        b.setNbMac(big_mac);
        b.allowMatrixMemory(6, big_mac);
        for (int j = 1; j <= 6; ++j)
            for (int m = 1; m <= big_mac; ++m)
                b.setTime(j, m, ptime(rng));
        // 5 jobs at the stack-buffer boundary (nbMac == PFSP_NEON_STACK_MACHINES,
        // last stack case) and one past it (first heap-fallback case); appended
        // after the instances above so their random times are unchanged
        for (int edge_mac : {PFSP_NEON_STACK_MACHINES, PFSP_NEON_STACK_MACHINES + 1}) {
            synth.emplace_back(new PfspInstance());
            PfspInstance& e = *synth.back();
            e.setSilence(true);
            e.setNbJob(5);
            e.setNbMac(edge_mac);
            e.allowMatrixMemory(5, edge_mac);
            for (int j = 1; j <= 5; ++j)
                for (int m = 1; m <= edge_mac; ++m)
                    e.setTime(j, m, ptime(rng));
        }
#endif
    }
    // Magnitude cases, 7 jobs x 20 machines, constant times. Largest completion
    // time is (nbJob + nbMac - 1) * p = 26 * p.
    //   p = 99:     2,574
    //   p = 100000: 2,600,000 < 2^24 = 16,777,216, so SSE float lanes are exact too
    synth.push_back(make_const_instance(7, 20, 99));
    synth.push_back(make_const_instance(7, 20, 100000));
#if defined(ENABLE_SIMD) && defined(__ARM_NEON)
    // NEON only (the SSE lanes are float): p = 20,000,001 gives a largest
    // completion time of 26 * 20,000,001 = 520,000,026, above 2^24 (float
    // lanes would round; p is odd so it is not even representable) but below
    // 2^31 (int32 lanes are exact).
    synth.push_back(make_const_instance(7, 20, 20000001));
#endif
    // in-process instances (setNbJob / setNbMac / allowMatrixMemory / setTime)
    // must get a correct view too: zero row 0 and column 0
    fail |= check_view(*synth.front());
    std::vector<PfspInstance*> all{&inst};
    for (const std::unique_ptr<PfspInstance>& s : synth)
        all.push_back(s.get());

    int f_ms = 0;
#ifdef ENABLE_SIMD
    int f_pm = 0, f_head = 0, f_ht = 0, f_tail = 0, f_id = 0;
#if defined(__ARM_NEON)
    fail |= check_neon_helpers();
#endif
#endif
    for (PfspInstance* in : all) {
        const std::vector<std::vector<int>> perms = test_permutations(in->getNbJob(), seed, rounds);
#ifdef ENABLE_SIMD
        f_pm |= check_partial_makespans(*in, perms);
        f_head |= check_head_tail(*in, perms, HT_HEAD);
        f_ht |= check_head_tail(*in, perms, HT_BOTH);
        f_tail |= check_head_tail(*in, perms, HT_TAIL);
        f_id |= check_head_tail_identity(*in, perms);
#endif
        f_ms |= check_makespan(*in, perms);
    }
#ifdef ENABLE_SIMD
    if (!f_pm) std::printf("partial_makespans: OK\n");
    if (!f_head) std::printf("head: OK\n");
    if (!f_ht) std::printf("head_and_tail: OK\n");
    if (!f_tail) std::printf("tail: OK\n");
    if (!f_id) std::printf("head_tail_identity: OK\n");
    fail |= f_pm | f_head | f_ht | f_tail | f_id;
#endif
    if (!f_ms) std::printf("makespan: OK\n");
    fail |= f_ms;
    return fail;
}
