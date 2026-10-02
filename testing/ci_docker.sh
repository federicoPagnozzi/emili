#!/bin/bash
# Run the Linux CI gate locally, inside Docker, without pushing.
#
# Executes testing/ci_linux.sh - the exact script .github/workflows/ci.yml
# runs on ubuntu-latest - in an Ubuntu 24.04 container with the repository
# bind-mounted. Build trees land in build-ci/ and build-ci-asan/ (gitignored)
# next to the macOS ones, so both can coexist.
#
# Usage: testing/ci_docker.sh            # full gate
#        JOBS=8 testing/ci_docker.sh     # more parallel build jobs
#        testing/ci_docker.sh bash       # drop into a shell in the image instead
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
IMAGE="${IMAGE:-emili-ci:ubuntu24}"

echo "=== building image ${IMAGE} (cached after the first run) ==="
docker build -q -t "${IMAGE}" "${ROOT}/testing/docker"

# --cap-add SYS_PTRACE: LeakSanitizer stops the process with ptrace at exit
# to scan for unreachable memory; the default seccomp profile forbids it.
if [ "$#" -gt 0 ]; then
    exec docker run --rm -it --cap-add SYS_PTRACE -v "${ROOT}:/work" -w /work \
        -e JOBS="${JOBS:-}" "${IMAGE}" "$@"
fi
exec docker run --rm -t --cap-add SYS_PTRACE -v "${ROOT}:/work" -w /work \
    -e JOBS="${JOBS:-}" "${IMAGE}" bash testing/ci_linux.sh
