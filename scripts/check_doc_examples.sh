#!/usr/bin/env bash
# check_doc_examples.sh — every documentation example COMPILES, or the gate fails.
#
# The house convention: every public space function's docstring carries an `@par Example` with a
# `@code{.purr} … @endcode` block, and each block is a COMPLETE program (its own imports, ready
# to copy-paste). This script extracts every block from the hand-authored module headers
# (space/space.hpp + space/*/*.hpp — cheatah-space's docstrings live in the C++, there are no
# `.purr` sources) and compiles each with purrc against this checkout — so an example that drifts
# from the API breaks the push instead of lying on the docs site. (Blocks are compiled, not run:
# compilation is the truth bar the docs promise.)
#
#   scripts/check_doc_examples.sh          # check every block
#   scripts/check_doc_examples.sh -v       # also print each block's source on failure
set -uo pipefail
cd "$(git rev-parse --show-toplevel)"
REPO="$PWD"
CHEATAH_DIR="${CHEATAH_DIR:-$(cd "$PWD/../cheatah" 2>/dev/null && pwd || true)}"
VERBOSE=0; [ "${1:-}" = "-v" ] && VERBOSE=1

find_tool() {
    local n="$1"
    for c in release debug asan; do
        [ -x "$CHEATAH_DIR/build/$c/bin/$n" ] && { echo "$CHEATAH_DIR/build/$c/bin/$n"; return 0; }
    done
    command -v "$n" 2>/dev/null
}
PURRC="$(find_tool purrc)"
[ -n "$PURRC" ] && [ -x "$PURRC" ] || {
    echo "doc-examples: no purrc toolchain (set CHEATAH_DIR or place cheatah at ../cheatah)"; exit 2; }

W="$(mktemp -d)"; trap 'rm -rf "$W"' EXIT
total=0; failed=0

for src in space/space.hpp space/*/*.hpp space/*/*/*.hpp; do
    [ -f "$src" ] || continue
    # Extract every @code{.purr} … @endcode block with its starting line number, stripping the
    # Javadoc block-comment's leading " * " from each code line.
    awk -v src="$src" '
        /@code\{\.purr\}/ { inblock = 1; start = NR; n = 0; next }
        /@endcode/ && inblock {
            file = sprintf("BLOCK\t%s\t%d", src, start); print file
            for (i = 1; i <= n; i++) print lines[i]
            print "ENDBLOCK"
            inblock = 0; next
        }
        inblock {
            line = $0
            sub(/^[[:space:]]*\*[[:space:]]?/, "", line)
            lines[++n] = line
        }
    ' "$src"
done > "$W/blocks.txt"

block_file=""; block_line=""; body="$W/body.purr"
while IFS= read -r line; do
    case "$line" in
        BLOCK*)
            block_file="$(printf '%s' "$line" | cut -f2)"
            block_line="$(printf '%s' "$line" | cut -f3)"
            : > "$body"
            ;;
        ENDBLOCK)
            total=$((total + 1))
            if ! out="$("$PURRC" --import-root "$REPO" "$body" -o "$W/ex.so" 2>&1)"; then
                failed=$((failed + 1))
                echo "doc-examples: FAILED — $block_file:$block_line does not compile:"
                echo "$out" | sed 's/^/    /' | head -12
                [ "$VERBOSE" = "1" ] && { echo "    ── block ──"; sed 's/^/    /' "$body"; }
            fi
            ;;
        *) printf '%s\n' "$line" >> "$body" ;;
    esac
done < "$W/blocks.txt"

if [ "$total" -eq 0 ]; then
    echo "doc-examples: no @code{.purr} blocks found — the convention requires them on every public fn."
    exit 1
fi
[ "$failed" -eq 0 ] || { echo "doc-examples: $failed of $total example blocks failed."; exit 1; }
echo "doc-examples: all $total example blocks compile."
