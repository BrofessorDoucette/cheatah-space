#!/usr/bin/env bash
# Fetch NASA's leap-second table and PRINT it as the C++ initializer, for a human to review.
#
#   scripts/update_leapseconds.sh
#
# It deliberately does NOT edit space/cdf/leapseconds.hpp. Adding a leap second changes the
# meaning of every TT2000 timestamp in the estate, and the house rule is that these headers are
# hand-authored — they ARE the source, which scripts/sign-modules.sh and the QA gate both state
# in as many words. A script that rewrote the header would make the table a generated artifact
# with different review rules, and would let a change to every timestamp we produce land without
# anyone reading it. Verification is automated (scripts/check_leapseconds.sh); authorship is not.
#
# What to do with the output:
#   1. Diff it against the kLeapSeconds initializer in space/cdf/leapseconds.hpp.
#   2. Paste in the new rows.
#   3. Update leap_seconds_sha256() (printed below), leap_seconds_upstream_updated(),
#      leap_seconds_verified_on() (today), and leap_seconds_known_good_through()
#      (today + 6 months — IERS Bulletin C's minimum notice).
#   4. Update the expected counts in tests/space_cdf_leapseconds_test.cpp; the pre-1972
#      divergence count against NASA's arithmetic must be RE-DERIVED, not adjusted to fit.
set -uo pipefail
cd "$(git rev-parse --show-toplevel)" || { echo "update-leapseconds: not in a git checkout"; exit 1; }

HEADER="space/cdf/leapseconds.hpp"
url=$(sed -n 's/.*return "\(https:\/\/cdf\.gsfc\.nasa\.gov[^"]*\)".*/\1/p' "$HEADER" | head -1)
[ -n "$url" ] || { echo "could not read leap_seconds_url() out of $HEADER" >&2; exit 1; }

tmp=$(mktemp); trap 'rm -f "$tmp" "$tmp.norm"' EXIT
curl -fsSL --max-time 60 -o "$tmp" "$url" || { echo "could not fetch $url" >&2; exit 1; }
awk '!/^;/ && NF { $1 = $1; print }' "$tmp" > "$tmp.norm"

rows=$(wc -l < "$tmp.norm")
drift=$(awk '$5 != 0 || $6 != 0' "$tmp.norm" | wc -l)
stamp=$(sed -n 's/^;[[:space:]]*Updated:[[:space:]]*\([0-9]\{8\}\).*/\1/p' "$tmp" | head -1)

printf '// --- paste into detail::kLeapSeconds in %s ---\n' "$HEADER"
# Columns: Year Month Day LeapSeconds Drift1(MJD) Drift2(s/day). The last three are scaled to
# exact integers: x1e7 for the two second-valued columns, so no double touches the table.
# NOTE the parentheses around the concatenations. In awk `+` binds TIGHTER than concatenation,
# so `a[1] frac + 0` means `a[1] (frac + 0)` — which silently turned "10" plus a "0000000" pad
# into "100" instead of 100000000, i.e. 10 microseconds where 10 seconds belongs. Also note
# `n = split(...)` rather than `length(array)`: only gawk implements the latter.
awk '{
    na = split($4, a, "."); frac = (na > 1) ? a[2] : ""
    while (length(frac) < 7) { frac = frac "0" }
    leap = (a[1] frac) + 0
    nb = split($6, b, "."); bfrac = (nb > 1) ? b[2] : ""
    while (length(bfrac) < 7) { bfrac = bfrac "0" }
    rate = (b[1] bfrac) + 0
    split($5, c, "."); mjd = c[1] + 0
    printf "    {%4d, %2d, %2d, %12d, %6d, %8d},\n", $1, $2, $3, leap, mjd, rate
}' "$tmp.norm"

printf '\n// leap_second_count       -> %s\n' "$rows"
printf '// leap_second_drift_rows  -> %s\n' "$drift"
printf '// leap_seconds_upstream_updated() -> "%s"\n' "${stamp:-UNKNOWN}"
printf '// leap_seconds_sha256()           -> "%s"\n' "$(sha256sum "$tmp.norm" | cut -d' ' -f1)"
printf '// leap_seconds_verified_on()      -> {%s}\n' "$(date -u +'%Y, %-m, %-d')"
printf '// leap_seconds_known_good_through -> {%s}  (today + 6 months)\n' \
       "$(date -u -d '+6 months' +'%Y, %-m, %-d' 2>/dev/null || echo 'compute by hand')"
