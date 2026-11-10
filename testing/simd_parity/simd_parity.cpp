// simd_parity: checks the SIMD data structures / kernels against the scalar code.
// usage: simd_parity <instance> [seed] [rounds]
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <numeric>
#include <random>
#include <vector>

#include "../../pfsp/pfspinstance.h"
#include "../../pfsp/pfsp_ptimes.h"
#include "../../pfsp/pfsp_simd.h"

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

#ifdef ENABLE_SIMD
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

static int check_partial_makespans_one(PfspInstance& inst, std::vector<int>& sol, const char* label)
{
    const int n = inst.getNbJob();
    const int nm = inst.getNbMac();
    const std::vector<std::vector<long>> C = ref_completion(inst, sol, n);
    std::vector<long> pmet(n + 1, 0);
    computePartialMakespans(sol, pmet, inst.getTimesView(), n, nm);
    for (int j = 1; j <= n; ++j) {
        if (pmet[j] != C[j][nm]) {
            std::printf("kernel=partial_makespans(%s) size=%d cell=[%d][%d] expected=%ld got=%ld\n",
                        label, n, nm, j, C[j][nm], pmet[j]);
            return 1;
        }
    }
    return 0;
}

static int check_partial_makespans(PfspInstance& inst, std::vector<int>& sol, unsigned seed, int rounds)
{
    if (check_partial_makespans_one(inst, sol, "identity"))
        return 1;
    std::mt19937 rng(seed);
    for (int r = 0; r < rounds; ++r) {
        std::shuffle(sol.begin() + 1, sol.end(), rng);
        if (check_partial_makespans_one(inst, sol, "shuffle"))
            return 1;
    }
    std::printf("partial_makespans: OK\n");
    return 0;
}
#endif

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
#ifdef ENABLE_SIMD
    // sol[0] = 0 is the unused slot, jobs are 1..nbJob
    std::vector<int> sol(inst.getNbJob() + 1);
    std::iota(sol.begin(), sol.end(), 0);
    fail |= check_partial_makespans(inst, sol, seed, rounds);
#else
    (void)seed;
    (void)rounds;
#endif
    return fail;
}
