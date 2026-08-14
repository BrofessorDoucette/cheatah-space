#!/usr/bin/env bash
# QA gate for cheatah-space — the quality checks that must pass before a push. Invoked by the git
# pre-push hook (.githooks/pre-push) and runnable by hand. Exits non-zero to BLOCK the push.
#
#   1. Unit coverage (hard gate): clang source-based coverage of cheatah_space_tests must be 100% lines
#      AND functions over space/**/*.hpp; the README coverage table must be committed in sync.
#   2. Doc coverage (hard gate): 100% Javadoc on the hand-authored space C++ API
#      (scripts/doc_coverage.sh), and every `@par Example` @code{.purr} block COMPILES
#      (scripts/check_doc_examples.sh).
#   3. Module sidecars (hard gate): every module header verifies against its .sha512 sidecar — the
#      checksums that make `biome add cheatah-space` resolvable on the extension path. (Sibling
#      extensions also check generated headers against their `.purr` sources here; cheatah-space's
#      headers are hand-authored C++ — they ARE the source — so sync means sidecars only.)
#   4. Configure + build (debug).
#   5. System tests (hard gate): every systests/test_*.purr compiles against the space package via
#      purrc and runs under the cheatah runtime, finishing with `RESULT: PASS`.
#   5b. Biome-install sandbox (hard gate): a fresh `biome add cheatah-space` user can compile + run it.
#   6. Unit tests (hard gate): ctest.
#   7. ASan + UBSan (hard gate): build + run the suite under sanitizers.
#   8. Valgrind memcheck (hard gate): run every unit test under Valgrind.
#   9. cppcheck (hard gate): performance + security static analysis.
#   10. Private-reference scan (hard gate): no non-public sibling-project names in the tree.
#
# The toolchain (purrc + cheatah runtime + the stdlib sources the C++ tests compile) is the sibling
# ../cheatah checkout; override with CHEATAH_DIR. Skips (discouraged; for fast local iteration):
# QA_GATE_SKIP_COVERAGE, QA_GATE_SKIP_DOCS, QA_GATE_SKIP_ASAN, QA_GATE_SKIP_VALGRIND. QA_GATE_SKIP=1
# bypasses everything.
set -uo pipefail
cd "$(git rev-parse --show-toplevel)"
REPO_ROOT="$(pwd)"
CHEATAH_DIR="${CHEATAH_DIR:-$(cd "$REPO_ROOT/../cheatah" 2>/dev/null && pwd || true)}"
export CHEATAH_DIR

if [ "${QA_GATE_SKIP:-0}" = "1" ]; then
    printf '\n[qa-gate] QA_GATE_SKIP=1 — skipping the QA gate (NOT recommended).\n'; exit 0
fi

bold()  { printf '\n\033[1m[qa-gate] %s\033[0m\n' "$*"; }
green() { printf '\033[32m%s\033[0m\n' "$*"; }
red()   { printf '\033[31m%s\033[0m\n' "$*"; }
fail()  { printf '\n\033[31m[qa-gate] FAILED: %s\033[0m\n' "$*"; exit 1; }

# 1. Unit coverage: refresh the README table, fail if it drifted, hard-fail unless 100% ----------
if [ "${QA_GATE_SKIP_COVERAGE:-0}" = "1" ]; then
    bold "Skipping coverage stage (QA_GATE_SKIP_COVERAGE=1)."
else
    bold "Measuring unit-test coverage (clang source-based) + refreshing the README table…"
    bash scripts/coverage.sh update-readme >/tmp/cheatah_space_coverage.log 2>&1 || { tail -30 /tmp/cheatah_space_coverage.log; fail "coverage report"; }
    if git rev-parse HEAD >/dev/null 2>&1 && ! git diff --quiet -- README.md; then
        printf '\n[qa-gate] The README coverage table is out of date. Updated it to:\n\n'
        git --no-pager diff -- README.md | sed -n '/coverage:start/,/coverage:end/p'
        fail "README coverage table changed — 'git add README.md && git commit', then push again"
    fi
    cat /tmp/cheatah_space_coverage.log
    covnums=$(sed -n 's/.*lines [0-9.]*% (\([0-9]*\)\/\([0-9]*\)), functions [0-9.]*% (\([0-9]*\)\/\([0-9]*\)).*/\1 \2 \3 \4/p' /tmp/cheatah_space_coverage.log)
    [ -n "$covnums" ] || fail "could not parse the coverage summary (coverage.sh output changed?)"
    read -r lcov_n lcov_d fcov_n fcov_d <<<"$covnums"
    if [ "$lcov_n" != "$lcov_d" ] || [ "$fcov_n" != "$fcov_d" ]; then
        fail "unit-test coverage below 100% — lines $lcov_n/$lcov_d, functions $fcov_n/$fcov_d (find gaps with: scripts/coverage.sh show <file>)"
    fi
    bold "Unit-test coverage: 100% lines ($lcov_n/$lcov_d) + functions ($fcov_n/$fcov_d)."
fi

# 2. Doc coverage: 100% Javadoc + every doc example compiles (hard gate) --------------------------
if [ "${QA_GATE_SKIP_DOCS:-0}" = "1" ]; then
    bold "Skipping documentation-coverage stage (QA_GATE_SKIP_DOCS=1)."
else
    bold "Checking documentation coverage (100% Javadoc)…"
    bash scripts/doc_coverage.sh || fail "documentation coverage below 100% — document the entities listed above"
    bold "Compiling the documentation examples (@par Example blocks)…"
    bash scripts/check_doc_examples.sh || fail "a documentation example does not compile — fix the @code{.purr} block"
fi

# 3. Module sidecars: every module header must verify against its .sha512 (what `biome add
#    cheatah-space` resolves on the extension path). The headers are hand-authored — they are the
#    source of truth — so there is no generated-header drift to check, only the signatures. --------
bold "Checking module sidecars are in sync…"
bash scripts/sign-modules.sh >/dev/null || fail "sign-modules.sh"
# Robust whether or not the tree is committed yet: every sidecar must verify against its header bytes.
while IFS= read -r sc; do
    d="$(dirname "$sc")"
    ( cd "$d" && sha512sum -c "$(basename "$sc")" >/dev/null 2>&1 ) || fail "module sidecar mismatch: $sc (run scripts/sign-modules.sh)"
done < <(find space -name '*.hpp.sha512')
# If the tree IS committed, also refuse a stale (uncommitted) sidecar so it can't be pushed out of sync.
if git rev-parse HEAD >/dev/null 2>&1; then
    [ -z "$(git status --porcelain -- 'space/*.hpp.sha512')" ] || fail "a module sidecar drifted — run scripts/sign-modules.sh, commit it, push again"
fi

# 4. Configure + build (debug) -------------------------------------------------------------------
bold "Configuring + building (debug)…"
cmake --preset debug         >/tmp/cheatah_space_cfg_debug.log   2>&1 || { tail -20 /tmp/cheatah_space_cfg_debug.log;   fail "configure (debug)"; }
cmake --build --preset debug >/tmp/cheatah_space_build_debug.log 2>&1 || { tail -30 /tmp/cheatah_space_build_debug.log; fail "debug build"; }

# 5. System tests: compile each systests/test_*.purr against space/ and run it --------------------
[ -n "${CHEATAH_DIR:-}" ] && [ -d "$CHEATAH_DIR" ] || fail "cannot find the cheatah toolchain (set CHEATAH_DIR or place it at ../cheatah)"
find_tool() { local n="$1"; for c in release debug asan; do
    [ -x "$CHEATAH_DIR/build/$c/bin/$n" ] && { echo "$CHEATAH_DIR/build/$c/bin/$n"; return 0; }; done; return 1; }
PURRC="$(find_tool purrc || true)"; CHEATAH="$(find_tool cheatah || true)"
if [ -z "$PURRC" ] || [ -z "$CHEATAH" ]; then
    bold "building the cheatah toolchain (no prebuilt purrc/cheatah found)…"
    cmake -S "$CHEATAH_DIR" -B "$CHEATAH_DIR/build/release" -DCMAKE_BUILD_TYPE=Release >/dev/null \
      && cmake --build "$CHEATAH_DIR/build/release" >/dev/null \
      || fail "failed to build the cheatah toolchain"
    PURRC="$(find_tool purrc)"; CHEATAH="$(find_tool cheatah)"
fi
bold "Running cheatah (.purr) system tests…  purrc=$PURRC"
WORK="$(mktemp -d)"; trap 'rm -rf "$WORK"' EXIT
# The shared test helper is pure cheatah; emit it to an importable header so `import check` resolves.
"$PURRC" --emit-library --transparent "$REPO_ROOT/systests/check.purr" -o "$WORK/check.hpp" 2>"$WORK/check.err" \
  || { red "failed to build the test helper (systests/check.purr):"; sed 's/^/    /' "$WORK/check.err"; exit 1; }
shopt -s nullglob
tests=("$REPO_ROOT"/systests/test_*.purr)
[ ${#tests[@]} -gt 0 ] || fail "no systests/test_*.purr — refusing to pass a gate with nothing to check"
# space.time is pure cheatah (no GPU, no SDK), so the systests compile against the headers alone —
# compiling each one also gates that the C++ templates instantiate cleanly (scalar + ndarray paths).
sfails=0; sran=0
for t in "${tests[@]}"; do
    name="$(basename "$t" .purr)"; sran=$((sran + 1)); mod="$WORK/$name.so"
    bold "── $name ──"
    if ! "$PURRC" --import-root "$REPO_ROOT" --import-root "$WORK" "$t" -o "$mod" 2>"$WORK/$name.err"; then
        red "  COMPILE FAILED"; sed 's/^/    /' "$WORK/$name.err"; sfails=$((sfails + 1)); continue
    fi
    out="$("$CHEATAH" "$mod" 2>&1)"; echo "$out" | sed 's/^/  /'
    if echo "$out" | grep -q "FAIL" || ! echo "$out" | grep -q "RESULT: PASS"; then
        red "  TEST FAILED"; sfails=$((sfails + 1))
    fi
done
[ "$sfails" -eq 0 ] || fail "$sfails of $sran system test file(s) failed"
green "[qa-gate] system tests: $sran/$sran green."

# 5b. Biome-install sandbox: prove a standard cheatah install can `biome add cheatah-space` and use
#     it (extension path + sidecar; no git, no --import-root) — the first-real-consumer contract. ----
bold "Sandboxing the 'biome add cheatah-space' user experience…"
bash scripts/test-biome-install.sh || fail "biome-install sandbox"

# 6. Unit tests (hard gate) ----------------------------------------------------------------------
# Exclude the `qa_gate` ctest entry itself — it shells back into THIS script (infinite recursion).
bold "Running unit test suite…"
ctest --preset debug --output-on-failure --exclude-regex '^qa_gate$' || fail "unit tests"

# 7. Sanitizers: ASan + UBSan (hard gate) --------------------------------------------------------
if [ "${QA_GATE_SKIP_ASAN:-0}" = "1" ]; then
    bold "Skipping sanitizer stage (QA_GATE_SKIP_ASAN=1)."
else
    bold "Configuring + building (ASan + UBSan)…"
    cmake --preset asan         >/tmp/cheatah_space_cfg_asan.log   2>&1 || { tail -20 /tmp/cheatah_space_cfg_asan.log;   fail "configure (asan)"; }
    cmake --build --preset asan >/tmp/cheatah_space_build_asan.log 2>&1 || { tail -30 /tmp/cheatah_space_build_asan.log; fail "asan build"; }
    bold "Running unit test suite under ASan + UBSan…"
    UBSAN_OPTIONS="print_stacktrace=1:halt_on_error=1" ASAN_OPTIONS="detect_leaks=1:abort_on_error=1" \
        ctest --preset asan --output-on-failure --exclude-regex '^qa_gate$' || fail "sanitizer (ASan/UBSan) tests"
fi

# 8. Valgrind memcheck (hard gate) ---------------------------------------------------------------
if [ "${QA_GATE_SKIP_VALGRIND:-0}" = "1" ]; then
    bold "Skipping Valgrind stage (QA_GATE_SKIP_VALGRIND=1)."
elif ! command -v valgrind >/dev/null 2>&1; then
    fail "valgrind not installed (install it, or set QA_GATE_SKIP_VALGRIND=1)"
else
    bold "Running unit tests under Valgrind memcheck…"
    bash security/run-valgrind.sh >/tmp/cheatah_space_valgrind.log 2>&1 || { tail -50 /tmp/cheatah_space_valgrind.log; fail "valgrind memcheck"; }
    tail -1 /tmp/cheatah_space_valgrind.log
fi

# 9. Static analysis: cppcheck (hard gate) -------------------------------------------------------
bold "Running cppcheck (performance + security)…"
bash scripts/cppcheck.sh || fail "cppcheck (performance/security findings)"

# 10. Private-reference scan (hard gate): this repo is public — no sibling-project names --------
bold "Scanning for private-project references…"
bash scripts/check_no_private_refs.sh || fail "a private-project reference is in the tree"

bold "QA gate PASSED — push may proceed."
exit 0
