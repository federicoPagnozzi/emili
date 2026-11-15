#!/bin/bash
# The Linux CI gate (improvement plan P2 Task 7).
#
# Runs the same steps in two places so they cannot drift:
#   - .github/workflows/ci.yml on ubuntu-latest (x86-64) and ubuntu-24.04-arm
#     (AArch64), every push / pull request
#   - testing/ci_docker.sh inside an Ubuntu container on a developer machine
#
# Steps: scalar release build -> simd_parity + goldens + negative tests; SIMD
# release build (-DUSE_SIMD=ON: SSE on x86-64, NEON on AArch64) -> simd_parity
# + goldens + negative tests; then scalar and SIMD sanitizer builds (ASan +
# UBSan, and LeakSanitizer which is only available on Linux) with
# -Wall -Wextra -Werror -> simd_parity + goldens (+ negative tests for scalar).
#
# Environment: JOBS (parallel build jobs), BUILD_REL / BUILD_ASAN /
# BUILD_SIMD / BUILD_SIMD_ASAN (build dirs, kept distinct from the macOS ones
# so both can coexist on a bind mount).
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${ROOT}"
JOBS="${JOBS:-$(nproc 2>/dev/null || echo 4)}"
BUILD_REL="${BUILD_REL:-build-ci}"
BUILD_ASAN="${BUILD_ASAN:-build-ci-asan}"
BUILD_SIMD="${BUILD_SIMD:-build-ci-simd}"
BUILD_SIMD_ASAN="${BUILD_SIMD_ASAN:-build-ci-simd-asan}"
PARITY_INSTANCE="testing/DD_Ta055.txt"

step() { echo; echo "=== $* ==="; }

# check_path <build dir> <scalar|simd>: the golden harness strips the startup
# banner, so run one cheap instance and assert the build really took (simd) or
# avoided (scalar) the SIMD code path.
check_path() {
    local build="$1" expect="$2" out
    out="$("./${build}/emili" testing/DD_Ta055.txt PFSP_MS first neh locmin insert rnds 42)"
    if [ "${expect}" = "simd" ]; then
        if ! printf '%s\n' "${out}" | grep -q 'SIMD run path enabled:'; then
            echo "FAIL: ${build} is a SIMD leg but did not report 'SIMD run path enabled:'"; exit 1
        fi
    else
        if printf '%s\n' "${out}" | grep -q 'SIMD run path enabled:'; then
            echo "FAIL: ${build} is a scalar leg but reported 'SIMD run path enabled:'"; exit 1
        fi
    fi
    echo "path check OK: ${build} is ${expect}"
}

step "toolchain"
cmake --version | head -1
"${CXX:-c++}" --version | head -1
echo "architecture: $(uname -m)"

step "normalize check"
bash testing/normalize_check.sh

step "release build (${BUILD_REL})"
cmake -S . -B "${BUILD_REL}" > /dev/null
cmake --build "${BUILD_REL}" -j"${JOBS}"
check_path "${BUILD_REL}" scalar

step "simd_parity (scalar release)"
"./${BUILD_REL}/simd_parity" "${PARITY_INSTANCE}"

step "goldens (release)"
./testing/golden_tests.sh check "./${BUILD_REL}/emili"

step "negative tests (release)"
./testing/negative_tests.sh "./${BUILD_REL}/emili"

step "simd release build (${BUILD_SIMD})"
cmake -S . -B "${BUILD_SIMD}" -DUSE_SIMD=ON > /dev/null
cmake --build "${BUILD_SIMD}" -j"${JOBS}"
check_path "${BUILD_SIMD}" simd

step "simd_parity (simd release)"
"./${BUILD_SIMD}/simd_parity" "${PARITY_INSTANCE}"

step "goldens (simd release)"
./testing/golden_tests.sh check "./${BUILD_SIMD}/emili"

step "negative tests (simd release)"
./testing/negative_tests.sh "./${BUILD_SIMD}/emili"

step "sanitizer build (${BUILD_ASAN}): ASan + UBSan + LSan, -Wall -Wextra -Werror"
cmake -S . -B "${BUILD_ASAN}" -DO3_FLAGS=OFF -DDEBUG_FLAGS=ON > /dev/null
cmake --build "${BUILD_ASAN}" -j"${JOBS}"
check_path "${BUILD_ASAN}" scalar

# halt_on_error: the first sanitizer report aborts the run, which the golden
# script records as a nonzero exit -> [DIFF] -> gate failure.
export ASAN_OPTIONS="detect_leaks=1:halt_on_error=1"
export UBSAN_OPTIONS="halt_on_error=1:print_stacktrace=1"
# Deliberate process-lifetime objects (P2.3 audit table) are listed in
# testing/lsan.supp; anything else that leaks fails the gate.
export LSAN_OPTIONS="suppressions=${ROOT}/testing/lsan.supp:print_suppressions=0"

step "simd_parity (sanitizers)"
"./${BUILD_ASAN}/simd_parity" "${PARITY_INSTANCE}"

step "goldens (sanitizers)"
./testing/golden_tests.sh check "./${BUILD_ASAN}/emili"

step "negative tests (sanitizers)"
./testing/negative_tests.sh "./${BUILD_ASAN}/emili"

step "simd sanitizer build (${BUILD_SIMD_ASAN}): ASan + UBSan + LSan, -Wall -Wextra -Werror"
cmake -S . -B "${BUILD_SIMD_ASAN}" -DO3_FLAGS=OFF -DDEBUG_FLAGS=ON -DUSE_SIMD=ON > /dev/null
cmake --build "${BUILD_SIMD_ASAN}" -j"${JOBS}"
check_path "${BUILD_SIMD_ASAN}" simd

step "simd_parity (simd sanitizers)"
"./${BUILD_SIMD_ASAN}/simd_parity" "${PARITY_INSTANCE}"

step "goldens (simd sanitizers)"
./testing/golden_tests.sh check "./${BUILD_SIMD_ASAN}/emili"

echo
echo "CI GATE PASSED"
