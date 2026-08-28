#!/usr/bin/env bash
# @test TAG gate for space.irbem — every public function names a test, and every test it names EXISTS.
#
# scripts/doc_coverage.sh proves a function is DOCUMENTED. It says nothing about whether the
# documented behaviour is CHECKED anywhere, and `@complexity`/`@alloc` are claims a reader has no
# way to audit: "O(1), no allocation" is prose until something measures it. The house doc contract
# closes that with a third tag — `@test <Suite.TestName>` — naming the GoogleTest case that
# exercises the entity, so the documentation carries a pointer to its own evidence.
#
# A tag like that rots in one specific way: the test gets renamed or deleted and the tag keeps
# claiming coverage that is no longer there. A DEAD TAG IS WORSE THAN A MISSING ONE — a missing tag
# is visibly missing, a dead tag reads as proof. So this gate checks both directions:
#
#   1. every public function in space/irbem/**.hpp carries at least one `@test`;
#   2. every name a `@test` gives is a real test case, as reported by the test binaries themselves
#      (`--gtest_list_tests`) — not by a grep of the sources, which would still pass for a test that
#      no longer compiles into anything;
#   3. one tag per line, per the house doc contract.
#
# The entity list comes from Doxygen's XML — the same parser scripts/doc_coverage.sh runs, so the
# two gates agree by construction about what "a public entity" is (`detail` namespaces are
# implementation and excluded by the Doxyfile; private members are excluded here). Writing a second,
# hand-rolled C++ parser to answer a question Doxygen already answers is how the two definitions
# drift apart. The XML is regenerated into a temp directory rather than reused from docs/xml,
# because a stale tree silently checks yesterday's headers — and it is a temp directory rather than
# docs/xml so two gate runs cannot race each other.
#
# The "test exists" side has two states between live and dead, and they are not the same defect:
#
#   • the test is declared in tests/*.cpp and its source file is NOT named in tests/CMakeLists.txt.
#     Then it compiles into nothing, runs nowhere, and the tag is as empty as a dead one. Hard fail.
#   • the test is declared in a source file tests/CMakeLists.txt DOES name, but this build tree did
#     not produce that binary — the GPU suite is its own target behind CHEATAH_SPACE_GPU, so a
#     machine with no device legitimately has no such binary to list. The test exists; this tree
#     cannot see it. Reported as a warning, not a failure, because a gate that cannot pass without
#     a GPU would be turned off on the machines that need it most.
#
#   scripts/check_doc_tags.sh          # check; exit 1 on any missing, dead or orphaned tag
#   scripts/check_doc_tags.sh report   # same, plus the full per-entity tag table
#
# Environment:
#   DOXYGEN               doxygen binary (default: `doxygen`, else ~/Tools/doxygen-1.16.1/bin)
#   SPACE_TEST_BINS       space-separated test binaries to list tests from
#                         (default: every executable matching build/*/bin/*_tests)
set -uo pipefail
cd "$(git rev-parse --show-toplevel)"

MODE="${1:-check}"

# The headers under the gate. space/time is a different module with a different contract; scoping
# this to space/irbem is deliberate, and widening it is a one-line change here plus the tags.
SCOPE="space/irbem/"

DOXYGEN="${DOXYGEN:-doxygen}"
command -v "$DOXYGEN" >/dev/null 2>&1 || DOXYGEN="$HOME/Tools/doxygen-1.16.1/bin/doxygen"
command -v "$DOXYGEN" >/dev/null 2>&1 || { echo "doc-tags: doxygen not found (install it, or set \$DOXYGEN)"; exit 1; }

work="$(mktemp -d)"
trap 'rm -rf "$work"' EXIT

# ---------------------------------------------------------------------------------------------
# 1. The entity list, from Doxygen's XML.
# ---------------------------------------------------------------------------------------------
( cat Doxyfile
  echo "OUTPUT_DIRECTORY=$work"
  echo "XML_OUTPUT=xml"
  echo "GENERATE_XML=YES"
  echo "GENERATE_HTML=NO"
  echo "GENERATE_LATEX=NO"
  echo "HAVE_DOT=NO"
  echo "QUIET=YES"
) | "$DOXYGEN" - 2>"$work/doxygen.log" >/dev/null

[ -d "$work/xml" ] || { echo "doc-tags: doxygen produced no XML"; sed -n '1,20p' "$work/doxygen.log"; exit 1; }

# Doxygen wraps each documented entity in <memberdef …>…</memberdef>, pretty-printed across many
# lines. Fold each one onto a single line first, so the rest is ordinary line-at-a-time text
# processing rather than an XML parser written in awk.
#
# Emitted, tab-separated: id, file, line, name, argsstring, then every @test name.
# `@test` is one of Doxygen's own commands, so it lands in the XML as an <xrefsect> titled "Test"
# with one <para> per tag — no ALIASES entry is needed in the Doxyfile for it.
cat "$work"/xml/*.xml \
  | tr '\n' ' ' \
  | sed -e 's|<memberdef |\n<memberdef |g' -e 's|</memberdef>|</memberdef>\n|g' \
  | grep '^<memberdef ' \
  | awk -v scope="$SCOPE" '
    function field(tag,   s, e) {              # <tag>…</tag>, first occurrence
        s = index($0, "<" tag ">");
        if (s == 0) return "";
        e = index(substr($0, s), "</" tag ">");
        if (e == 0) return "";
        return substr($0, s + length(tag) + 2, e - length(tag) - 3);
    }
    function attr(name,   s, rest, e) {        # name="…", first occurrence
        s = index($0, name "=\"");
        if (s == 0) return "";
        rest = substr($0, s + length(name) + 2);
        e = index(rest, "\"");
        return substr(rest, 1, e - 1);
    }
    {
        kind = attr("kind"); prot = attr("prot");
        # Free and member functions, plus the friend operators the frame types declare inline.
        if (kind != "function" && kind != "friend") next;
        if (prot != "public") next;

        # <location file="…" line="…"> — the declaration site, which is what a reader greps for.
        loc = index($0, "<location ");
        if (loc == 0) next;
        rest = $0; $0 = substr($0, loc); file = attr("file"); line = attr("line"); $0 = rest;
        if (index(file, scope) != 1) next;

        id = attr("id"); name = field("name"); args = field("argsstring");
        gsub(/&amp;/, "\\&", args); gsub(/&lt;/, "<", args); gsub(/&gt;/, ">", args);
        qname = field("qualifiedname"); if (qname != "") name = qname;

        # The @test names: every <para> inside the xrefsect Doxygen titles "Test".
        tags = ""; sect = $0;
        s = index(sect, "<xrefsect id=\"test_");
        if (s > 0) {
            sect = substr(sect, s);
            e = index(sect, "</xrefsect>");
            if (e > 0) sect = substr(sect, 1, e);
            while ((s = index(sect, "<para>")) > 0) {
                sect = substr(sect, s + 6);
                e = index(sect, "</para>");
                if (e == 0) break;
                t = substr(sect, 1, e - 1);
                sect = substr(sect, e + 7);
                gsub(/^[ \t]+/, "", t); gsub(/[ \t]+$/, "", t);
                if (t != "") tags = tags "\t" t;
            }
        }
        print id "\t" file "\t" line "\t" name "\t" args tags;
    }' \
  | sort -u -t'	' -k1,1 > "$work/entities.tsv"

entities=$(wc -l < "$work/entities.tsv")
[ "$entities" -gt 0 ] || { echo "doc-tags: parsed no public functions under $SCOPE — the XML shape changed"; exit 1; }

# ---------------------------------------------------------------------------------------------
# 2. The tests that actually exist, from the binaries and (for diagnosis only) from the sources.
# ---------------------------------------------------------------------------------------------
bins="${SPACE_TEST_BINS:-}"
if [ -z "$bins" ]; then
    for b in build/*/bin/*_tests; do [ -x "$b" ] && bins="$bins $b"; done
fi
[ -n "${bins// /}" ] || { echo "doc-tags: no test binary found (build one, or set \$SPACE_TEST_BINS)"; exit 1; }

: > "$work/live.txt"
for b in $bins; do
    [ -x "$b" ] || { echo "doc-tags: \$SPACE_TEST_BINS names '$b', which is not executable"; exit 1; }
    # gtest prints "Suite." at column 0 and "  TestName" indented, with an optional
    # "  # TypeParam = …" trailer on either.
    "$b" --gtest_list_tests 2>/dev/null | sed 's/[ \t]*#.*$//' | awk '
        /^[A-Za-z_]/  { suite = $1; next }
        /^[ \t]+[A-Za-z_]/ { gsub(/^[ \t]+|[ \t]+$/, ""); if ($0 != "") print suite $0 }' >> "$work/live.txt"
done
sort -u "$work/live.txt" -o "$work/live.txt"
live=$(wc -l < "$work/live.txt")
[ "$live" -gt 0 ] || { echo "doc-tags: the test binaries listed no tests — is '$bins' a GoogleTest binary?"; exit 1; }

# Every TEST/TEST_F/TEST_P in the suite sources, tagged with the file that declares it. Used ONLY
# to tell a stale tag apart from a test this build tree happens not to have built. Recorded as
# "Suite.Name<TAB>file" so the diagnosis can name the source and check it against the CMake list.
grep -HoE '\bTEST(_F|_P)?\([ \t]*[A-Za-z_][A-Za-z0-9_]*[ \t]*,[ \t]*[A-Za-z_][A-Za-z0-9_]*' tests/*.cpp 2>/dev/null \
  | awk '{ i = index($0, ":"); file = substr($0, 1, i - 1); rest = substr($0, i + 1);
           s = substr(rest, index(rest, "(") + 1);
           gsub(/[ \t]/, "", s); sub(/,/, ".", s);
           print s "\t" file }' \
  | sort -u > "$work/declared.txt"

# ---------------------------------------------------------------------------------------------
# 3. Compare.
# ---------------------------------------------------------------------------------------------
: > "$work/missing.txt"; : > "$work/dead.txt"; : > "$work/orphan.txt"
: > "$work/malformed.txt"; : > "$work/unbuilt.txt"
tagged=0
while IFS=$'\t' read -r id file line name args rest; do
    if [ -z "$rest" ]; then
        printf '%s:%s  %s%s\n' "$file" "$line" "$name" "$args" >> "$work/missing.txt"
        continue
    fi
    tagged=$((tagged + 1))
    printf '%s' "$rest" | tr '\t' '\n' | while read -r t; do
        [ -n "$t" ] || continue
        case "$t" in
            *.*) ;;
            *) printf '%s:%s  %s — @test "%s" is not a Suite.TestName\n' \
                   "$file" "$line" "$name" "$t" >> "$work/malformed.txt"; continue ;;
        esac
        grep -qxF "$t" "$work/live.txt" && continue

        src="$(awk -F'\t' -v want="$t" '$1 == want { print $2; exit }' "$work/declared.txt")"
        if [ -z "$src" ]; then
            printf '%s:%s  %s — @test %s\n' "$file" "$line" "$name" "$t" >> "$work/dead.txt"
        elif ! grep -qF "$(basename "$src")" tests/CMakeLists.txt; then
            printf '%s:%s  %s — @test %s (declared in %s, which no target compiles)\n' \
                "$file" "$line" "$name" "$t" "$src" >> "$work/orphan.txt"
        else
            printf '%s:%s  %s — @test %s (in %s; that target was not built in this tree)\n' \
                "$file" "$line" "$name" "$t" "$src" >> "$work/unbuilt.txt"
        fi
    done
done < "$work/entities.tsv"

# The house doc contract is one tag per line: a second tag after a `@test` on the same line is
# invisible to a reader skimming for evidence, and Doxygen folds it into the tag's own text.
: > "$work/twotags.txt"
while IFS= read -r hit; do
    [ -n "$hit" ] && printf '%s\n' "$hit" >> "$work/twotags.txt"
done < <(grep -nE '@test\b' -r --include='*.hpp' "$SCOPE" \
         | grep -E '@test\b.*@[a-z]+\b|@[a-z]+\b[^*]*@test\b' || true)

# ---------------------------------------------------------------------------------------------
# 4. Report.
# ---------------------------------------------------------------------------------------------
if [ "$MODE" = "report" ]; then
    echo "--- entity → @test ---"
    while IFS=$'\t' read -r id file line name args rest; do
        printf '%-52s %s\n' "$name" "$(printf '%s' "$rest" | sed 's/^\t//; s/\t/, /g')"
    done < "$work/entities.tsv" | sort
    echo "---"
fi

fail=0
report_block() {   # $1 = file, $2 = headline, $3 = remedy
    local n; n=$(wc -l < "$1")
    [ "$n" -eq 0 ] && return 0
    echo "doc-tags: FAIL — $n $2:"
    sed 's/^/    /' "$1"
    echo "    $3"
    echo
    fail=1
}

report_block "$work/missing.txt" "public function(s) under $SCOPE carry no @test tag" \
    "Add '@test <Suite.TestName>' naming the test that genuinely exercises each one."
report_block "$work/dead.txt" "@test tag(s) name a test that does not exist" \
    "The test was renamed or deleted; a tag claiming coverage that is not there is a defect."
report_block "$work/orphan.txt" "@test tag(s) name a test no target compiles" \
    "Its source file is not in tests/CMakeLists.txt, so the test never runs anywhere — add it there."
report_block "$work/malformed.txt" "@test tag(s) are not of the form Suite.TestName" \
    "One tag, one test case, spelled exactly as --gtest_list_tests reports it."
report_block "$work/twotags.txt" "line(s) put a second doc tag on a @test line" \
    "House contract: one tag per line."

unbuilt=$(wc -l < "$work/unbuilt.txt")
if [ "$unbuilt" -ne 0 ]; then
    echo "doc-tags: NOTE — $unbuilt @test tag(s) name a test this build tree did not build:"
    sed 's/^/    /' "$work/unbuilt.txt"
    echo "    The test is declared and its file IS in tests/CMakeLists.txt, so it is not a dead tag."
    echo "    To put it under the gate here, configure that target (the GPU suite needs"
    echo "    -DCHEATAH_SPACE_GPU=ON) and re-run, or point \$SPACE_TEST_BINS at its binary."
    echo
fi

if [ "$fail" -ne 0 ]; then
    printf 'doc-tags: %s public function(s) under %s, %s tagged; %s test case(s) live in %s\n' \
        "$entities" "$SCOPE" "$tagged" "$live" "$(echo $bins | tr ' ' ',')"
    exit 1
fi

printf 'doc-tags: 100%% — all %s public functions under %s carry a @test, and every one of them names a live test case (%s cases listed from%s).\n' \
    "$entities" "$SCOPE" "$live" "$(echo " $bins")"
exit 0
