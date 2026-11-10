#!/bin/bash
# Benchmark emili binaries against each other (e.g. master scalar vs branch
# scalar vs branch NEON) and print a GitHub-flavoured Markdown table.
#
# usage: testing/bench_simd.sh <label>=<binary> [<label>=<binary> ...]
# env:   RUNS (default 5)  timed runs per row and binary, median reported
#
# Rows: every line of testing/golden_configs.txt (on its own instance) plus the
# generated 200x20 / 500x20 rows below.  Binaries are interleaved per run, not
# per binary.  Rows whose median is under 0.05 s are timed as batches of 20
# consecutive runs (per-run mean of the median batch) and marked with '*'.
# The objective and the iteration counter must be identical across binaries for
# every row; otherwise the script exits 1.  Progress goes to stderr.
set -u

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
RUNS="${RUNS:-5}"
BATCH=20
FAST=0.05
SPREAD=0.15

if [ $# -lt 1 ]; then
    echo "usage: $0 <label>=<binary> ..." >&2
    exit 2
fi
LABELS=(); BINS=()
for a in "$@"; do
    case "$a" in
        *=*) LABELS+=("${a%%=*}"); BINS+=("${a#*=}") ;;
        *) echo "bad argument '$a' (want label=binary)" >&2; exit 2 ;;
    esac
done
NB=${#BINS[@]}
for b in "${BINS[@]}"; do
    [ -x "$b" ] || { echo "binary not executable: $b" >&2; exit 2; }
done

TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT

# now: seconds as a float, portable (python3 or perl).
if command -v python3 >/dev/null 2>&1; then
    now() { python3 -c 'import time;print("%.6f"%time.time())'; }
else
    now() { perl -MTime::HiRes=time -e 'printf "%.6f\n", time'; }
fi

# timed <binary> <instance> <count> <args...>: wall seconds for <count>
# consecutive runs (output discarded).  Timestamps are taken by a single
# python3/perl process that spawns the runs, so interpreter start-up is not
# part of the measured interval.
timed() {
    local bin="$1" inst="$2" count="$3"; shift 3
    if command -v python3 >/dev/null 2>&1; then
        python3 - "$bin" "$inst" "$count" "$@" <<'PY'
import subprocess, sys, time
bin_, inst, count = sys.argv[1], sys.argv[2], int(sys.argv[3])
cmd = [bin_, inst] + sys.argv[4:]
t = time.perf_counter()
for _ in range(count):
    subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
print("%.6f" % (time.perf_counter() - t))
PY
    else
        perl -MTime::HiRes=time -e '
            my ($bin,$inst,$count,@a)=@ARGV; my $t=time;
            for (1..$count) { system("$bin $inst @a >/dev/null 2>&1") }
            printf "%.6f\n", time-$t;' "$bin" "$inst" "$count" "$@"
    fi
}

# median / min / max of whitespace separated numbers -> "med min max"
stats() {
    printf '%s\n' "$@" | sort -g | awk '{v[NR]=$1} END{
        m=(NR%2)?v[(NR+1)/2]:(v[NR/2]+v[NR/2+1])/2; printf "%.6f %.6f %.6f", m, v[1], v[NR]}'
}

ROWS_NAME=(); ROWS_INST=(); ROWS_ARGS=()
add_row() { ROWS_NAME+=("$1"); ROWS_INST+=("$2"); ROWS_ARGS+=("$3"); }
while IFS='|' read -r name inst args; do
    case "$name" in ''|'#'*) continue ;; esac
    add_row "$name" "testing/$inst" "$args"
done < "${ROOT}/testing/golden_configs.txt"
for sz in 200x20 500x20; do
    I="testing/bench/rand_${sz}.txt"
    add_row "${sz}_tainsert" "$I" "PFSP_MS first neh locmin tainsert rnds 42"
    add_row "${sz}_insert" "$I" "PFSP_MS first neh locmin insert rnds 42"
    add_row "${sz}_csinsert" "$I" "PFSP_MS first neh locmin csinsert rnds 42"
    add_row "${sz}_nehff" "$I" "PFSP_MS first nehff locmin insert rnds 42"
    add_row "${sz}_ils_igoper" "$I" "PFSP_MS ils first neh locmin tainsert maxstep 30 igoper 4 improve rnds 42"
    add_row "${sz}_ils_rsper" "$I" "PFSP_MS ils first neh locmin tainsert maxstep 30 rsper 4 improve rnds 42"
    add_row "${sz}_wt_tainsert" "$I" "PFSP_WT first neh locmin tainsert rnds 7"
done

cd "$ROOT"
host_cpu=""
if command -v sysctl >/dev/null 2>&1; then host_cpu="$(sysctl -n machdep.cpu.brand_string 2>/dev/null || true)"; fi
if [ -z "$host_cpu" ] && command -v lscpu >/dev/null 2>&1; then
    host_cpu="$(lscpu | sed -n 's/^Model name:[[:space:]]*//p' | head -1)"
fi
echo "# SIMD benchmark"
echo
echo "- machine: $(uname -m), ${host_cpu:-unknown CPU}"
echo "- date: $(date '+%Y-%m-%d %H:%M:%S %Z')"
echo "- runs per cell: ${RUNS} (median); rows marked * are per-run means of batches of ${BATCH}"
echo
printf '| row | instance |'
for l in "${LABELS[@]}"; do printf ' %s (s) |' "$l"; done
echo ' speed-up |'
printf '|---|---|'
for l in "${LABELS[@]}"; do printf -- '---:|'; done
echo '---:|'

fail=0
NR_=${#ROWS_NAME[@]}
r=0
while [ $r -lt $NR_ ]; do
    name="${ROWS_NAME[$r]}"; inst="${ROWS_INST[$r]}"; args="${ROWS_ARGS[$r]}"
    r=$((r+1))
    echo "[$r/$NR_] $name" >&2

    # correctness: objective and iteration counter must agree across binaries
    ref=""
    for ((i=0;i<NB;i++)); do
        out="$("${BINS[$i]}" "$inst" $args 2>/dev/null)"
        sig="$(printf '%s\n' "$out" | grep -E '^(Objective function value:|iteration counter :)' | tr '\n' ';')"
        if [ -z "$sig" ]; then
            echo "FAIL: row $name: no objective/iteration output from ${LABELS[$i]}" >&2; exit 1
        fi
        if [ -z "$ref" ]; then ref="$sig"; reflabel="${LABELS[$i]}"
        elif [ "$sig" != "$ref" ]; then
            echo "FAIL: row $name: trajectories diverged: ${reflabel} gave [$ref] but ${LABELS[$i]} gave [$sig]" >&2
            exit 1
        fi
    done

    measure() { # sets MED[i] SPR[i] and BATCHED
        local n per i k t first
        BATCHED=0
        # probe single-run median on the first binary to decide batching
        first=()
        for ((k=0;k<3;k++)); do first+=("$(timed "${BINS[0]}" "$inst" 1 $args)"); done
        set -- $(stats "${first[@]}")
        if awk -v m="$1" -v f="$FAST" 'BEGIN{exit !(m<f)}'; then BATCHED=1; fi
        if [ $BATCHED = 1 ]; then n=$BATCH; else n=1; fi
        for ((i=0;i<NB;i++)); do : > "$TMP/t$i"; done
        for ((k=0;k<RUNS;k++)); do
            for ((i=0;i<NB;i++)); do
                t="$(timed "${BINS[$i]}" "$inst" $n $args)"
                awk -v t="$t" -v n="$n" 'BEGIN{printf "%.6f\n", t/n}' >> "$TMP/t$i"
            done
        done
        for ((i=0;i<NB;i++)); do
            set -- $(stats $(cat "$TMP/t$i"))
            MED[$i]=$1
            SPR[$i]=$(awk -v a="$2" -v b="$3" -v m="$1" 'BEGIN{printf "%.4f", (m>0)?(b-a)/m:0}')
        done
    }
    MED=(); SPR=()
    measure
    noisy=0
    for ((i=0;i<NB;i++)); do
        awk -v s="${SPR[$i]}" -v lim="$SPREAD" 'BEGIN{exit !(s>lim)}' && noisy=1
    done
    note=""
    if [ $noisy = 1 ]; then
        echo "  spread > 15% for $name, rerunning once" >&2
        measure
        note=" (rerun)"
    fi
    mark=""; [ $BATCHED = 1 ] && mark="*"
    printf '| %s%s%s | %s |' "$name" "$mark" "$note" "$(basename "$inst")"
    for ((i=0;i<NB;i++)); do printf ' %.3f |' "${MED[$i]}"; done
    awk -v a="${MED[0]}" -v b="${MED[$((NB-1))]}" 'BEGIN{printf " %.2f |\n", (b>0)?a/b:0}'
done
exit $fail
