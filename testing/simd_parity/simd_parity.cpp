// simd_parity: checks the SIMD data structures / kernels against the scalar code.
// usage: simd_parity <instance> [seed] [rounds]
#include <cstdio>
#include <cstdlib>

#include "../../pfsp/pfspinstance.h"
#include "../../pfsp/pfsp_ptimes.h"

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
    return check_view(inst);
}
