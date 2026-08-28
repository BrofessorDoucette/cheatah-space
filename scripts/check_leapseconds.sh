#!/usr/bin/env bash
# Is our embedded leap-second table still what NASA publishes?
#
#   scripts/check_leapseconds.sh          # compare; WARN (exit 0) when offline
#   scripts/check_leapseconds.sh --strict # a fetch failure is also a failure
#
# space/cdf/leapseconds.hpp carries NASA's table compiled in, because a header-only numeric
# library that opens a file or reads an environment variable mid-conversion has global state, and
# the same CDF would then decode to different timestamps on two machines. The cost of that choice
# is that the table can go stale, and this is what pays it: a cheap out-of-band check that the
# snapshot still matches upstream.
#
# NETWORK-OPTIONAL BY DESIGN. With no network this WARNS and exits 0. That is deliberate: the QA
# gate must pass on a laptop on a plane, and a gate that fails when the wifi is down is a gate
# people learn to bypass with --no-verify — which loses every other check too, not just this one.
# What actually keeps the table current is the scheduled CI job running this with a network.
#
# The comparison is against a NORMALIZED form — comments and blank lines dropped, whitespace
# collapsed — so upstream reformatting or an edited comment raises no false alarm, while any
# change to the DATA does.
set -uo pipefail
cd "$(git rev-parse --show-toplevel)" || { echo "check-leapseconds: not in a git checkout"; exit 1; }

HEADER="space/cdf/leapseconds.hpp"
STRICT=0
[ "${1:-}" = "--strict" ] && STRICT=1

warn() { printf '\033[33m[leapseconds] %s\033[0m\n' "$*"; }
fail() { printf '\n\033[31m[leapseconds] FAILED: %s\033[0m\n' "$*" >&2; exit 1; }

[ -f "$HEADER" ] || fail "missing $HEADER"

# The URL and the expected digest are single-sourced from the header itself, so they cannot drift
# apart from the table they describe.
url=$(sed -n 's/.*return "\(https:\/\/cdf\.gsfc\.nasa\.gov[^"]*\)".*/\1/p' "$HEADER" | head -1)
want=$(sed -n 's/^    return "\([0-9a-f]\{64\}\)";$/\1/p' "$HEADER" | head -1)
[ -n "$url" ]  || fail "could not read leap_seconds_url() out of $HEADER"
[ -n "$want" ] || fail "could not read leap_seconds_sha256() out of $HEADER"

tmp=$(mktemp); trap 'rm -f "$tmp" "$tmp.norm"' EXIT
if ! curl -fsSL --max-time 60 -o "$tmp" "$url" 2>/dev/null; then
    if [ "$STRICT" = "1" ]; then fail "could not fetch $url (--strict)"; fi
    warn "could not reach $url — skipping the currency check (offline?)."
    warn "The embedded table is used as-is. Run this with a network before a release."
    exit 0
fi

# Normalize exactly as leap_seconds_sha256() documents.
awk '!/^;/ && NF { $1 = $1; print }' "$tmp" > "$tmp.norm"
got=$(sha256sum "$tmp.norm" | cut -d' ' -f1)
rows=$(wc -l < "$tmp.norm")

if [ "$got" = "$want" ]; then
    echo "[leapseconds] current — upstream still matches the embedded table ($rows rows)."
    exit 0
fi

printf '\n\033[31m[leapseconds] UPSTREAM HAS CHANGED\033[0m\n\n'
printf '  expected %s\n  got      %s\n  rows     %s\n\n' "$want" "$got" "$rows"
printf '  A leap second may have been added, or NASA may have corrected an entry. Either way this\n'
printf '  changes the meaning of every TT2000 timestamp we convert, so it is a REVIEWED change:\n\n'
printf '    1. scripts/update_leapseconds.sh          # prints the new table for review\n'
printf '    2. paste it into %s, update\n' "$HEADER"
printf '       leap_seconds_sha256(), leap_seconds_verified_on() and\n'
printf '       leap_seconds_known_good_through()\n'
printf '    3. update the expected row count in tests/space_cdf_leapseconds_test.cpp\n\n'
printf '  New rows upstream that we do not have:\n\n'
awk '!/^;/ && NF { $1 = $1; print "      " $0 }' "$tmp" | tail -5
printf '\n'
fail "the embedded leap-second table is out of date"
