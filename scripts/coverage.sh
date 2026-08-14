#!/usr/bin/env bash
# Measure cheatah-space test coverage with clang source-based coverage — the in-process unit tests
# (cheatah_space_tests) that exercise the hand-authored space/ headers.
#
#   scripts/coverage.sh               # per-file summary report for space/
#   scripts/coverage.sh show <file>   # uncovered lines of one file, e.g. space/time/time.hpp
#   scripts/coverage.sh funcs <file>  # per-function coverage of one file
#   scripts/coverage.sh update-readme # rewrite the coverage table in README.md
#
# The summary comes straight from `llvm-cov report`'s TOTAL row. space.time's templates are
# single-expression functions with no internal branches, so per-instantiation line accounting is
# exact — every instantiation the test TU creates is also executed by the suite (the suite is built
# that way on purpose), which is precisely what 100% lines + functions asserts. No merged-view
# post-processing is needed (and the QA tooling stays bash + LLVM, per the house rule).
set -uo pipefail
cd "$(git rev-parse --show-toplevel)"

B=build/cov
# clang source-based coverage (-fprofile-instr-generate -fcoverage-mapping is clang syntax), so pin
# clang for this build tree regardless of the machine's default compiler. Override via CLANGXX/CLANG.
CLANGXX="${CLANGXX:-clang++}"; CLANG="${CLANG:-clang}"
command -v "$CLANGXX" >/dev/null 2>&1 || { echo "coverage: $CLANGXX not found — clang is required for source-based coverage"; exit 1; }
cmake -S . -B "$B" -G Ninja -DCMAKE_BUILD_TYPE=Debug -DCHEATAH_SPACE_BUILD_TESTS=ON \
  -DCMAKE_C_COMPILER="$CLANG" -DCMAKE_CXX_COMPILER="$CLANGXX" \
  -DCMAKE_CXX_FLAGS="-fprofile-instr-generate -fcoverage-mapping" \
  -DCMAKE_EXE_LINKER_FLAGS="-fprofile-instr-generate -fcoverage-mapping" >/tmp/cheatah_space_cov_cfg.log 2>&1 \
  || { tail -20 /tmp/cheatah_space_cov_cfg.log; exit 1; }
cmake --build "$B" --target cheatah_space_tests >/tmp/cheatah_space_cov_build.log 2>&1 \
  || { tail -25 /tmp/cheatah_space_cov_build.log; exit 1; }

( cd "$B"
  LLVM_PROFILE_FILE=t1.profraw ./bin/cheatah_space_tests >/dev/null 2>&1
  llvm-profdata merge -sparse t1.profraw -o merged.profdata )

OBJS=(./"$B"/bin/cheatah_space_tests)
PROF="-instr-profile=$B/merged.profdata"
# The hand-authored, host-testable space surface. git's '*' spans '/', so this matches
# space/space.hpp and the submodule headers alike. `--cached --others --exclude-standard` also lists
# not-yet-committed headers while still honouring .gitignore.
SRCS=$(git ls-files --cached --others --exclude-standard 'space/*.hpp')

case "${1:-report}" in
    show)  llvm-cov show   "${OBJS[@]}" $PROF "${2:?usage: coverage.sh show <file>}" 2>/dev/null \
             | grep -nE '\|[[:space:]]*0\|' || echo "all lines covered in ${2}" ;;
    funcs) llvm-cov report "${OBJS[@]}" $PROF -show-functions "${2:?usage: coverage.sh funcs <file>}" 2>/dev/null ;;
    update-readme)
        # llvm-cov report TOTAL columns (after the filename): Regions MissedRegions RCover
        # Functions MissedFuncs FExec Lines MissedLines LCover Branches MissedBranches BCover.
        total=$(llvm-cov report "${OBJS[@]}" $PROF $SRCS 2>/dev/null | awk '$1=="TOTAL"{$1="";print}')
        [ -n "$total" ] || { echo "coverage: llvm-cov produced no TOTAL row"; exit 1; }
        read -r regions mreg rcov funcs mfun fexec lines mlin lcov branches mbr bcov <<<"$total"
        lcn=$((lines - mlin)); fcn=$((funcs - mfun))
        # Rewrite the marker-bounded coverage table in README.md (awk, in place via a temp file).
        awk -v lcov="$lcov" -v lcn="$lcn" -v lt="$lines" \
            -v fexec="$fexec" -v fcn="$fcn" -v ft="$funcs" \
            -v rcov="$rcov" -v bcov="$bcov" '
            /<!-- coverage:start -->/ {
                print "<!-- coverage:start -->"
                print "| Metric | space package |"
                print "|--------|---------------|"
                printf "| **Lines** | %s (%d/%d) |\n", lcov, lcn, lt
                printf "| **Functions** | %s (%d/%d) |\n", fexec, fcn, ft
                printf "| Regions | %s |\n", rcov
                printf "| Branches | %s |\n", bcov
                print "<!-- coverage:end -->"
                inblock = 1; seen = 1; next
            }
            /<!-- coverage:end -->/ && inblock { inblock = 0; next }
            !inblock { print }
            END { exit seen ? 0 : 2 }
        ' README.md > README.md.tmp || { rm -f README.md.tmp; echo "coverage markers not found in README.md"; exit 1; }
        mv README.md.tmp README.md
        echo "README coverage table: lines $lcov ($lcn/$lines), functions $fexec ($fcn/$funcs)"
        ;;
    *)     llvm-cov report "${OBJS[@]}" $PROF $SRCS 2>/dev/null ;;
esac
