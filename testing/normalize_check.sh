#!/bin/bash
# Checks that golden_tests.sh normalize() drops the SIMD banner line (any ISA).
set -u
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
eval "$(sed -n '/^normalize()/,/^}/p' "${SCRIPT_DIR}/golden_tests.sh")"
out="$(printf 'x\n\tSIMD run path enabled: NEON\ny\n' | normalize)"
[ "${out}" = $'x\ny' ] && echo "normalize_check: OK" || { echo "normalize_check: FAIL"; exit 1; }
