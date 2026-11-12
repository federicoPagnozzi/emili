#!/bin/bash
# Run the Linux CI gate locally, inside Docker, without pushing.
#
# Executes testing/ci_linux.sh - the exact script .github/workflows/ci.yml
# runs on ubuntu-latest and ubuntu-24.04-arm - in an Ubuntu 24.04 container
# with the repository bind-mounted. Build trees land in build-ci*/ (gitignored)
# next to the macOS ones, so both can coexist.
#
# Usage: testing/ci_docker.sh            # full gate, host architecture
#        JOBS=8 testing/ci_docker.sh     # more parallel build jobs
#        PLATFORM=linux/amd64 testing/ci_docker.sh   # force x86-64 (emulated
#                                        # on Apple silicon: slow) or
#        PLATFORM=linux/arm64 ...        # AArch64
#        testing/ci_docker.sh bash       # drop into a shell in the image instead
#
# With PLATFORM set the image is tagged emili-ci:ubuntu24-<arch> and the build
# directories get an -<arch> suffix, so amd64 and arm64 containers sharing the
# same bind-mounted checkout never reuse each other's CMake caches.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PLATFORM="${PLATFORM:-}"

PLATFORM_ARGS=()
BUILD_ARGS=()
if [ -n "${PLATFORM}" ]; then
    ARCH="${PLATFORM#linux/}"
    IMAGE="${IMAGE:-emili-ci:ubuntu24-${ARCH}}"
    PLATFORM_ARGS=(--platform "${PLATFORM}")
    BUILD_ARGS=(
        -e BUILD_REL="build-ci-${ARCH}"
        -e BUILD_ASAN="build-ci-asan-${ARCH}"
        -e BUILD_SIMD="build-ci-simd-${ARCH}"
        -e BUILD_SIMD_ASAN="build-ci-simd-asan-${ARCH}"
    )
else
    IMAGE="${IMAGE:-emili-ci:ubuntu24}"
fi

echo "=== building image ${IMAGE} (cached after the first run) ==="
docker build -q ${PLATFORM_ARGS[@]+"${PLATFORM_ARGS[@]}"} -t "${IMAGE}" "${ROOT}/testing/docker"

# --cap-add SYS_PTRACE: LeakSanitizer stops the process with ptrace at exit
# to scan for unreachable memory; the default seccomp profile forbids it.
if [ "$#" -gt 0 ]; then
    exec docker run --rm -it ${PLATFORM_ARGS[@]+"${PLATFORM_ARGS[@]}"} --cap-add SYS_PTRACE \
        -v "${ROOT}:/work" -w /work \
        -e JOBS="${JOBS:-}" ${BUILD_ARGS[@]+"${BUILD_ARGS[@]}"} "${IMAGE}" "$@"
fi
exec docker run --rm -t ${PLATFORM_ARGS[@]+"${PLATFORM_ARGS[@]}"} --cap-add SYS_PTRACE \
    -v "${ROOT}:/work" -w /work \
    -e JOBS="${JOBS:-}" ${BUILD_ARGS[@]+"${BUILD_ARGS[@]}"} "${IMAGE}" bash testing/ci_linux.sh
