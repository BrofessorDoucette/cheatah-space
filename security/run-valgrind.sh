#!/usr/bin/env bash
# Run cheatah-space's unit tests under Valgrind memcheck — a second memory checker alongside ASan
# (each catches things the other can miss). Uses the plain `debug` build, because Valgrind cannot run
# an ASan-instrumented binary.
#
#   security/run-valgrind.sh
set -uo pipefail
cd "$(git rev-parse --show-toplevel)"

SUPP="security/valgrind.supp"

echo "[valgrind] building (debug)…"
cmake --preset debug          >/tmp/cheatah_space_vg_cfg.log   2>&1 || { tail -20 /tmp/cheatah_space_vg_cfg.log;   exit 1; }
cmake --build --preset debug  >/tmp/cheatah_space_vg_build.log 2>&1 || { tail -30 /tmp/cheatah_space_vg_build.log; exit 1; }

VG=(valgrind --tool=memcheck --leak-check=full
    --show-leak-kinds=definite,indirect
    --errors-for-leak-kinds=definite,indirect
    --error-exitcode=1 --suppressions="$SUPP")

UNIT_BINS=(cheatah_space_tests)

status=0
total_ran=0
for t in "${UNIT_BINS[@]}"; do
    bin="./build/debug/bin/$t"
    [ -x "$bin" ] || { echo "[valgrind] missing $bin — cannot cover all unit tests"; status=1; continue; }
    log="/tmp/cheatah_space_vg_$t.log"
    echo "[valgrind] memcheck: $t"
    if ! "${VG[@]}" "$bin" >"$log" 2>&1; then
        echo "[valgrind] ERRORS/LEAKS in $t:"; tail -50 "$log"; status=1
    fi
    # Coverage assertion: confirm every test actually executed under Valgrind (a crash/abort would run
    # fewer than registered, silently shrinking coverage).
    ran=$(grep -oE '\[=+\] [0-9]+ test' "$log" | tail -1 | grep -oE '[0-9]+' || true)
    passed=$(grep -oE '\[ *PASSED *\] [0-9]+ test' "$log" | tail -1 | grep -oE '[0-9]+' || true)
    failed=$(grep -cE '\[ *FAILED *\]' "$log" || true)
    : "${ran:=0}"; : "${passed:=0}"; : "${failed:=0}"
    if [ "$ran" -eq 0 ]; then
        echo "[valgrind] $t executed 0 tests under Valgrind — coverage incomplete"; status=1
    elif [ "$passed" -ne "$ran" ] || [ "$failed" -ne 0 ]; then
        echo "[valgrind] $t: only $passed/$ran tests passed under Valgrind (failures: $failed)"; status=1
    else
        echo "[valgrind] $t: $ran/$ran unit tests executed clean under Valgrind"
        total_ran=$((total_ran + ran))
    fi
done

if [ "$status" -eq 0 ]; then
    echo "[valgrind] 100% unit-test coverage: $total_ran tests, no errors or definite/indirect leaks"
fi
exit "$status"
