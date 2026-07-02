#!/usr/bin/env bash
# qa_gate.sh — the gate that guards every push (see .githooks/pre-push).
#
# It compiles each cheatah test in tests/ against the in-repo `space` package and runs it
# through the cheatah runtime, failing if ANY test fails to compile, prints a FAIL line, or
# does not finish with `RESULT: PASS`. Because compiling a test instantiates the C++ module
# templates (time/time.hpp) through purrc, this also gates that the C++ extension compiles
# cleanly for both the scalar and the ndarray (vectorized) concept paths.
#
# Cross-platform note: this is the POSIX/bash gate (Linux, macOS, WSL, Git-Bash). A native
# Windows PowerShell gate can mirror it; the pre-push hook prefers this when bash is present.
#
# Override the toolchain location with CHEATAH_DIR (default: the sibling ../cheatah checkout).
set -u

# --- locate the repo and the toolchain ---------------------------------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
CHEATAH_DIR="${CHEATAH_DIR:-$(cd "$REPO_ROOT/../cheatah" 2>/dev/null && pwd || true)}"

red()   { printf '\033[31m%s\033[0m\n' "$*"; }
green() { printf '\033[32m%s\033[0m\n' "$*"; }
bold()  { printf '\033[1m%s\033[0m\n' "$*"; }

if [ -z "${CHEATAH_DIR:-}" ] || [ ! -d "$CHEATAH_DIR" ]; then
    red "qa_gate: cannot find the cheatah toolchain (set CHEATAH_DIR or place it at ../cheatah)."
    exit 2
fi

# Find a built purrc + cheatah runtime; prefer release, then debug, then asan.
find_tool() {
    local name="$1"
    for cfg in release debug asan; do
        if [ -x "$CHEATAH_DIR/build/$cfg/bin/$name" ]; then
            echo "$CHEATAH_DIR/build/$cfg/bin/$name"; return 0
        fi
    done
    return 1
}
PURRC="$(find_tool purrc || true)"
CHEATAH="$(find_tool cheatah || true)"

if [ -z "$PURRC" ] || [ -z "$CHEATAH" ]; then
    bold "qa_gate: building the cheatah toolchain (no prebuilt purrc/cheatah found)…"
    cmake -S "$CHEATAH_DIR" -B "$CHEATAH_DIR/build/release" -DCMAKE_BUILD_TYPE=Release >/dev/null \
        && cmake --build "$CHEATAH_DIR/build/release" --target purrc cheatah >/dev/null || {
            red "qa_gate: failed to build the cheatah toolchain."; exit 2; }
    PURRC="$(find_tool purrc)"; CHEATAH="$(find_tool cheatah)"
fi

bold "qa_gate: purrc   = $PURRC"
bold "qa_gate: runtime = $CHEATAH"

# --- run every cheatah test --------------------------------------------------------------
WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT
fails=0
ran=0

# The shared test helper is pure cheatah; emit it to an importable header so `import check`
# resolves (a sibling .purr is not auto-imported). It is dev-only, so it lives in WORK and is
# added as a second import root alongside the repo (which provides the `space` package).
if ! "$PURRC" --emit-library --transparent "$REPO_ROOT/tests/check.purr" -o "$WORK/check.hpp" 2>"$WORK/check.err"; then
    red "qa_gate: failed to build the test helper (tests/check.purr):"
    sed 's/^/    /' "$WORK/check.err"
    exit 1
fi

shopt -s nullglob
tests=("$REPO_ROOT"/tests/test_*.purr)
if [ ${#tests[@]} -eq 0 ]; then
    red "qa_gate: no tests found in tests/ — refusing to pass a gate with nothing to check."
    exit 1
fi

for t in "${tests[@]}"; do
    name="$(basename "$t" .purr)"
    ran=$((ran + 1))
    mod="$WORK/$name.so"
    bold "── $name ──────────────────────────────────────────"
    if ! "$PURRC" --import-root "$REPO_ROOT" --import-root "$WORK" "$t" -o "$mod" 2>"$WORK/$name.err"; then
        red "  COMPILE FAILED"; sed 's/^/    /' "$WORK/$name.err"; fails=$((fails + 1)); continue
    fi
    out="$("$CHEATAH" "$mod" 2>&1)"
    echo "$out" | sed 's/^/  /'
    if echo "$out" | grep -q "FAIL" || ! echo "$out" | grep -q "RESULT: PASS"; then
        red "  TEST FAILED"; fails=$((fails + 1))
    fi
done

echo
if [ "$fails" -eq 0 ]; then
    green "qa_gate: PASS — $ran test file(s) green."
    exit 0
fi
red "qa_gate: FAIL — $fails of $ran test file(s) failed."
exit 1
