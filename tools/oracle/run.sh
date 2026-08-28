#!/usr/bin/env bash
# Build the IRBEM oracle and regenerate the measurements in space/irbem/docs/ERROR_BUDGET.md.
#
# DEV-ONLY, and never part of the QA gate. IRBEM is LGPL-3.0 and cheatah-space is MIT, so the
# library is built here, run as a BLACK BOX through its documented C entry points, and never
# linked into anything we ship. The checkout lives in the git-ignored space/irbem/vendor/.
#
# It is built TWICE on purpose:
#   as-shipped  — IRBEM's own FFLAGS, which contain no -O at all
#   -O2         — what we benchmark against, because a speedup quoted against an unoptimized
#                 reference is inflated before any of our own work is measured
#
#   tools/oracle/run.sh            # build both, run the study
#   tools/oracle/run.sh --no-build # re-run the study against existing builds
set -uo pipefail
cd "$(git rev-parse --show-toplevel)"

VENDOR="space/irbem/vendor/IRBEM"
OUT="${IRBEM_ORACLE_OUT:-/tmp/irbem-builds}"
AS_SHIPPED_FFLAGS="-fpic -fno-second-underscore -std=legacy -ffixed-line-length-none"
O2_FFLAGS="$AS_SHIPPED_FFLAGS -O2 -ffp-contract=off -fno-fast-math"

fail() { printf '\033[31m[oracle] FAILED: %s\033[0m\n' "$*"; exit 1; }
note() { printf '\n\033[1m[oracle] %s\033[0m\n' "$*"; }

[ -d "$VENDOR" ] || fail "no oracle at $VENDOR — clone https://github.com/PRBEM/IRBEM there (git-ignored)"
command -v gfortran >/dev/null 2>&1 || fail "gfortran not found — 'sudo apt install gfortran'"
mkdir -p "$OUT"

if [ "${1:-}" != "--no-build" ]; then
    for lane in asshipped O2; do
        case "$lane" in
            asshipped) flags="$AS_SHIPPED_FFLAGS" ;;
            O2)        flags="$O2_FFLAGS" ;;
        esac
        note "building the oracle ($lane): FFLAGS=$flags"
        ( cd "$VENDOR" && make clean >/dev/null 2>&1
          make OS=linux64 ENV=gfortran64 -j"$(nproc)" FFLAGS="$flags" ) >"$OUT/build_$lane.log" 2>&1 \
            || { tail -20 "$OUT/build_$lane.log"; fail "oracle build ($lane)"; }
        cp "$VENDOR/bin/libirbem.linux64.gfortran64.so" "$OUT/libirbem-$lane.so" \
            || fail "oracle build ($lane) produced no library"
    done
fi

note "building the convergence driver"
g++ -O2 -std=c++20 -Wall -Wextra tools/oracle/convergence.cpp -ldl -o "$OUT/convergence" \
    || fail "convergence driver"

note "convergence study (-O2 oracle) -> $OUT/convergence-O2.txt"
"$OUT/convergence" "$OUT/libirbem-O2.so" >"$OUT/convergence-O2.txt" 2>&1 || fail "convergence run"
sed -n '1,12p' "$OUT/convergence-O2.txt"

note "discretization at IRBEM's recommended resolution (options(3)=options(4)=0)"
awk '$3=="0" && $4=="0" && NF>=7 {print $6}' "$OUT/convergence-O2.txt" | sort -g \
  | awk '{a[NR]=$1} END {printf "  L* relative move vs (9,9): min %s  max %s  (n=%d)\n", a[1], a[NR], NR}'

note "performance baseline — mean wall time per L* evaluation at default resolution"
for lane in asshipped O2; do
    [ -f "$OUT/libirbem-$lane.so" ] || continue
    "$OUT/convergence" "$OUT/libirbem-$lane.so" 2>/dev/null \
      | awk -v L="$lane" '$3=="0" && $4=="0" {s+=$NF; n++}
             END {if (n>0) printf "  %-10s %6.1f ms/point  (n=%d)\n", L, s/n*1000, n}'
done

note "done — fold these into space/irbem/docs/ERROR_BUDGET.md (§2 and §5)"
