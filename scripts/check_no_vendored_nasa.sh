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

echo "[no-vendored-nasa] clean — no vendored reference implementation is tracked."
