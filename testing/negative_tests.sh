#!/bin/bash
# Parser error-path tests (P1.2): bad input => nonzero exit + message on stderr.
set -u
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
EMILI_BIN="${1:-${SCRIPT_DIR}/../build/emili}"
INSTANCE="${SCRIPT_DIR}/DD_Ta055.txt"
FAILED=0

expect_error() {
    local desc="$1"; shift
    local pattern="$1"; shift
    err="$("${EMILI_BIN}" "$@" 2>&1 >/dev/null)"
    status=$?
    if [ "${status}" -eq 0 ]; then
        echo "[FAIL] ${desc}: expected nonzero exit"; FAILED=1; return
    fi
    if ! printf '%s' "${err}" | grep -q "${pattern}"; then
        echo "[FAIL] ${desc}: stderr missing '${pattern}'"; echo "${err}" | head -5; FAILED=1; return
    fi
    echo "[OK] ${desc} (exit ${status})"
}

expect_error "unknown neighborhood token" "FATAL ERROR" \
    "${INSTANCE}" PFSP_MS first neh locmin bogus_neighborhood rnds 42
expect_error "truncated algorithm description" "FATAL ERROR" \
    "${INSTANCE}" PFSP_MS ils first neh locmin insert rnds 42
expect_error "non-numeric where number expected" "EXPECTED" \
    "${INSTANCE}" PFSP_MS ils first neh locmin insert maxstep notanumber rndmv insert 3 improve rnds 42

[ "${FAILED}" -eq 0 ] && echo "All negative tests passed." && exit 0
exit 1
