#!/bin/bash
# Golden-trajectory regression harness (improvement plan P0.1).
#
# For each line of golden_configs.txt (format: name|instance|args) runs
#   emili <testing/instance> <args>
# captures stdout and stderr separately, normalizes lines that legitimately
# vary between runs (wall time, commit id), and either records the result as
# a golden or diffs it against the stored golden.
#
# Usage:
#   golden_tests.sh record [emili_binary]   # (re)record goldens - needs explicit sign-off
#   golden_tests.sh check  [emili_binary]   # regression gate; exit 1 on any diff
set -u
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MODE="${1:-check}"
EMILI_BIN="${2:-${SCRIPT_DIR}/../build/emili}"
GOLDEN_DIR="${SCRIPT_DIR}/goldens"
CONFIGS="${SCRIPT_DIR}/golden_configs.txt"
TMP_DIR="$(mktemp -d)"
trap 'rm -rf "${TMP_DIR}"' EXIT
FAILED=0

if [ ! -x "${EMILI_BIN}" ]; then
    echo "emili binary not found/executable: ${EMILI_BIN}" >&2
    exit 2
fi
mkdir -p "${GOLDEN_DIR}"

normalize() {
    # The instance path is printed by the loaders; replace the absolute
    # testing directory so goldens are checkout-independent (CI runs them
    # from /work, developers from wherever the repo lives).
    sed -E -e 's/^time : .*/time : <NORMALIZED>/' \
           -e 's/^CPU time: .*/CPU time: <NORMALIZED>/' \
           -e 's/^commit : .*/commit : <NORMALIZED>/' \
           -e "s#${SCRIPT_DIR}#<TESTING_DIR>#g"
}

while IFS='|' read -r name instance args; do
    case "${name}" in ''|'#'*) continue;; esac
    stdout_f="${TMP_DIR}/${name}.out"
    stderr_f="${TMP_DIR}/${name}.err"
    "${EMILI_BIN}" "${SCRIPT_DIR}/${instance}" ${args} > "${stdout_f}" 2> "${stderr_f}"
    status=$?
    combined="${TMP_DIR}/${name}.combined"
    {
        echo "=== exit ${status} ==="
        echo "=== stdout ==="
        normalize < "${stdout_f}"
        echo "=== stderr ==="
        normalize < "${stderr_f}"
    } > "${combined}"
    golden="${GOLDEN_DIR}/${name}.golden"
    if [ "${MODE}" = "record" ]; then
        cp "${combined}" "${golden}"
        echo "[RECORDED] ${name} (exit ${status})"
        if [ "${status}" -ne 0 ]; then
            echo "  WARNING: nonzero exit while recording - config probably invalid" >&2
            FAILED=1
        fi
    else
        if [ ! -f "${golden}" ]; then
            echo "[MISSING] ${name} - run 'golden_tests.sh record' first" >&2
            FAILED=1
            continue
        fi
        if diff -q "${golden}" "${combined}" > /dev/null; then
            echo "[OK] ${name}"
        else
            echo "[DIFF] ${name}"
            diff -u "${golden}" "${combined}" | head -60
            FAILED=1
        fi
    fi
done < "${CONFIGS}"

if [ "${FAILED}" -ne 0 ]; then
    echo "GOLDEN CHECK FAILED" >&2
    exit 1
fi
echo "All goldens match."
exit 0
