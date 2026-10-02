#!/bin/bash
# The Linux CI gate (improvement plan P2 Task 7).
#
# Runs the same steps in two places so they cannot drift:
#   - .github/workflows/ci.yml on ubuntu-latest (every push / pull request)
#   - testing/ci_docker.sh inside an Ubuntu container on a developer machine
#
# Steps: release build -> goldens + negative tests; then a sanitizer build
# (ASan + UBSan, and LeakSanitizer which is only available on Linux) with
# -Wall -Wextra -Werror -> goldens + negative tests again.
#
# Environment: JOBS (parallel build jobs), BUILD_REL / BUILD_ASAN (build dirs,
# kept distinct from the macOS ones so both can coexist on a bind mount).
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${ROOT}"
JOBS="${JOBS:-$(nproc 2>/dev/null || echo 4)}"
BUILD_REL="${BUILD_REL:-build-ci}"
BUILD_ASAN="${BUILD_ASAN:-build-ci-asan}"

step() { echo; echo "=== $* ==="; }

step "toolchain"
cmake --version | head -1
"${CXX:-c++}" --version | head -1

step "release build (${BUILD_REL})"
cmake -S . -B "${BUILD_REL}" > /dev/null
cmake --build "${BUILD_REL}" -j"${JOBS}"

step "goldens (release)"
./testing/golden_tests.sh check "./${BUILD_REL}/emili"

step "negative tests (release)"
./testing/negative_tests.sh "./${BUILD_REL}/emili"

step "sanitizer build (${BUILD_ASAN}): ASan + UBSan + LSan, -Wall -Wextra -Werror"
cmake -S . -B "${BUILD_ASAN}" -DO3_FLAGS=OFF -DDEBUG_FLAGS=ON > /dev/null
cmake --build "${BUILD_ASAN}" -j"${JOBS}"

# halt_on_error: the first sanitizer report aborts the run, which the golden
# script records as a nonzero exit -> [DIFF] -> gate failure.
export ASAN_OPTIONS="detect_leaks=1:halt_on_error=1"
export UBSAN_OPTIONS="halt_on_error=1:print_stacktrace=1"
# Deliberate process-lifetime objects (P2.3 audit table) are listed in
# testing/lsan.supp; anything else that leaks fails the gate.
export LSAN_OPTIONS="suppressions=${ROOT}/testing/lsan.supp:print_suppressions=0"

step "goldens (sanitizers)"
./testing/golden_tests.sh check "./${BUILD_ASAN}/emili"

step "negative tests (sanitizers)"
./testing/negative_tests.sh "./${BUILD_ASAN}/emili"

echo
echo "CI GATE PASSED"
