#!/usr/bin/env bash
# Hard gate: no third-party reference implementation is ever committed to this repo.
#
# cheatah-space builds space.cdf and space.irbem from published specifications and papers, and
# keeps the reference implementations (NASA's CDF library, PRBEM/IRBEM) strictly as OPTIONAL,
# DEV-ONLY oracles under space/<mod>/vendor/, which .gitignore excludes.
#
# That arrangement is load-bearing for two separate reasons, and a convention alone will not hold
# it: someone eventually runs `git add -f`, or edits .gitignore, and nothing notices.
#
#   LICENSING. NASA's CDF distribution is NOT public domain. Per CDF_copyright.txt it "may be
#   copied or redistributed as long as it is not sold for profit" and carries modification-notice
#   requirements. cheatah-space ships under its own LICENSE; committing NASA's source would put
#   an incompatible redistribution term inside a repo that does not carry one.
#
#   PROVENANCE. "We implemented this from the spec" stops being checkable the moment their source
#   is sitting in the tree. Keeping it out is what keeps the clean-room claim true.
#
# Offline, no network, runs in milliseconds. Checks the INDEX and the working tree, not history.
set -uo pipefail
cd "$(git rev-parse --show-toplevel)" || { echo "check-no-vendored-nasa: not in a git checkout"; exit 1; }

fail() { printf '\n\033[31m[no-vendored-nasa] FAILED: %s\033[0m\n' "$*" >&2; exit 1; }

# 1. Nothing under any vendor/ directory may be tracked.
tracked_vendor="$(git ls-files -- 'space/*/vendor/*' 'space/*/vendor' 2>/dev/null)"
if [ -n "$tracked_vendor" ]; then
    printf '\n[no-vendored-nasa] These vendored files are TRACKED but must never be committed:\n\n'
    printf '%s\n' "$tracked_vendor" | sed 's/^/    /'
    printf '\nRemove them from the index (keeping them on disk):\n\n'
    printf '    git rm -r --cached space/cdf/vendor space/irbem/vendor\n\n'
    fail "vendored reference implementation is staged for commit"
fi

# 2. .gitignore must still exclude them. A deleted line here is how (1) silently becomes possible.
for d in space/cdf/vendor/ space/irbem/vendor/; do
    grep -qxF "$d" .gitignore || fail ".gitignore no longer lists '$d' — restore it (see the '# --- cheatah-space ---' block)"
done

# 3. No tracked file may carry NASA's copyright notice. This catches the case (1) cannot: a file
#    copied OUT of vendor/ into the source tree, where .gitignore does not apply. Checked against
#    tracked files only, and this script excludes itself so the sentinel below is not self-matching.
needle='This software may be copied or redistributed as long as it is not sold'
hits="$(git grep -l -F "$needle" -- . ':(exclude)scripts/check_no_vendored_nasa.sh' 2>/dev/null)"
if [ -n "$hits" ]; then
    printf '\n[no-vendored-nasa] These tracked files carry the NASA CDF copyright notice:\n\n'
    printf '%s\n' "$hits" | sed 's/^/    /'
    printf '\nspace.cdf is implemented from the published Internal Format Description. If you need\n'
    printf 'to consult the reference implementation, read it in space/cdf/vendor/ — do not copy it in.\n\n'
    fail "third-party source copied into the tracked tree"
fi

# 4. No tracked file may cite IRBEM's FORTRAN SOURCE as where a value came from.
#
#    This is the space.irbem counterpart of (3), and it catches a failure (3) cannot: not source
#    copied in, but a NUMBER read out of source/*.f and transcribed with a citation. That is a
#    provenance breach even though not one line of their code is present, because the model is
#    then derived from LGPL-3 source rather than from the papers — and the derivation trail says
#    so in writing.
#
#    It is not hypothetical. wave4-wip's ext_t01.hpp carries three such citations naming
#    `source/Tsyganenko01.f` and specific line numbers. They were written honestly, by an agent
#    following a rule that permitted "a single number", and they are exactly why that rule is now
#    "read vendor/IRBEM/docs/source/**.rst, never source/*.f". Nothing noticed for a month.
#
#    The oracle is a BLACK BOX: dlopen it, call its documented C entry points, measure what comes
#    back. A constant that cannot be reached that way or from a paper does not get transcribed —
#    the model documents the gap instead. See space/irbem/docs/VERIFICATION.md.
fortran_cite='source/[A-Za-z0-9_]+\.f\b'
hits="$(git grep -l -E "$fortran_cite" -- . ':(exclude)scripts/check_no_vendored_nasa.sh' 2>/dev/null)"
if [ -n "$hits" ]; then
    printf '\n[no-vendored-nasa] These tracked files cite IRBEM Fortran source as a provenance:\n\n'
    git grep -n -E "$fortran_cite" -- . ':(exclude)scripts/check_no_vendored_nasa.sh' | sed 's/^/    /'
    printf '\nThe oracle is read as a BLACK BOX, never as source. Re-derive the value from the\n'
    printf 'published paper or by probing the library through its C entry points, and cite THAT.\n'
    printf 'If neither reaches it, the model documents a measured gap instead of transcribing.\n\n'
    fail "a value is attributed to IRBEM's Fortran source"
fi

echo "[no-vendored-nasa] clean — no vendored reference implementation is tracked, and no value cites their source."
