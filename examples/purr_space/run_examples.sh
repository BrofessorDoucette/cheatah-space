#!/usr/bin/env bash
# Compile and run every example, writing PNGs into out/.
#
# Not part of the QA gate, and deliberately skippable: the examples need the cheatah toolchain,
# the sibling cheatah-plot checkout, and a corpus file that is fetched rather than committed.
# Any of those being absent is a SKIP, not a failure — a gate that fails on a fresh clone is a
# gate people learn to bypass.
#
#   examples/purr_space/run_examples.sh
#   CHEATAH_PLOT_DIR=../cheatah-plot examples/purr_space/run_examples.sh
set -uo pipefail
cd "$(git rev-parse --show-toplevel)"
REPO="$PWD"

CHEATAH_DIR="${CHEATAH_DIR:-$(cd "$REPO/../cheatah" 2>/dev/null && pwd || true)}"
find_tool() { local n="$1"; for c in release debug asan; do
    [ -x "$CHEATAH_DIR/build/$c/bin/$n" ] && { echo "$CHEATAH_DIR/build/$c/bin/$n"; return 0; }; done; return 1; }
PURRC="$(find_tool purrc)"; CHEATAH="$(find_tool cheatah)"
[ -n "$PURRC" ] && [ -n "$CHEATAH" ] || { echo "[examples] no cheatah toolchain — skipping."; exit 0; }

PLOT_DIR="${CHEATAH_PLOT_DIR:-$(cd "$REPO/../cheatah-plot" 2>/dev/null && pwd || true)}"
ROOTS=(--import-root "$REPO")
if [ -n "$PLOT_DIR" ] && [ -f "$PLOT_DIR/plot/plot.hpp" ]; then
    ROOTS+=(--import-root "$PLOT_DIR")
else
    echo "[examples] cheatah-plot not found — the plotting example will be skipped."
fi

CORPUS="$REPO/space/cdf/vendor/corpus/tier1/omni_hro2_1min_20150101_v01.cdf"
[ -f "$CORPUS" ] || echo "[examples] corpus absent — run: scripts/cdf-corpus.sh fetch --tier 1"

W="$(mktemp -d)"; trap 'rm -rf "$W"' EXIT
mkdir -p examples/purr_space/out
fails=0
for ex in examples/purr_space/[0-9]*.purr; do
    name="$(basename "$ex" .purr)"
    if grep -q '^import plot' "$ex" && [ ${#ROOTS[@]} -lt 4 ]; then
        echo "[examples] SKIP $name (needs cheatah-plot)"; continue
    fi
    if ! "$PURRC" "${ROOTS[@]}" "$ex" -o "$W/$name.so" 2>"$W/$name.err"; then
        echo "[examples] COMPILE FAILED $name"; sed 's/^/    /' "$W/$name.err" | head -10; fails=$((fails + 1)); continue
    fi
    # Run from the example's own directory: the programs use relative paths for both the corpus
    # and their output, so that a reader can copy a line straight out of one.
    if ! ( cd examples/purr_space && "$CHEATAH" "$W/$name.so" ); then
        echo "[examples] RUN FAILED $name"; fails=$((fails + 1))
    fi
done
[ "$fails" -eq 0 ] || { echo "[examples] $fails example(s) failed."; exit 1; }
echo "[examples] all examples ran."
