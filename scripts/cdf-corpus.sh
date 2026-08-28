#!/usr/bin/env bash
# The CDF differential corpus: real files from NASA's archive, plus the generated matrix that
# covers what the archive does not.
#
#   scripts/cdf-corpus.sh fetch [--tier N]   # download + verify (default: tiers 0 and 1)
#   scripts/cdf-corpus.sh verify             # offline: hash what is on disk against the manifest
#   scripts/cdf-corpus.sh probe              # re-derive the manifest's tag column from the files
#   scripts/cdf-corpus.sh derive             # GENERATE tier3 with NASA's cdfconvert
#   scripts/cdf-corpus.sh status
#   scripts/cdf-corpus.sh clean [--all]
#
# The committed part is space/cdf/corpus/manifest.tsv. The files are not: they live in the
# .gitignore'd space/cdf/vendor/corpus/ and are fetched on demand. Tests that need an absent file
# SKIP with the command to fetch it — never fail for absence, so the offline QA gate stays green.
#
# WHY TIER3 EXISTS. CDF defines four compression modes. Not one file in NASA's public archive
# uses three of them: every compressed variable in all 16 manifest rows is GZIP, and none is
# whole-file (CCR) compressed. RLE, Huffman, adaptive Huffman and CCR would therefore ship
# completely untested against real bytes. `derive` closes that by asking NASA's own cdfconvert to
# re-encode a known file into every mode, which also makes those files an oracle: we know exactly
# what values must come back out, because we know what went in.
#
# Tier3 is NOT in the committed manifest. Its bytes depend on the oracle build, so pinning a hash
# would pin our own toolchain rather than an upstream artifact; it gets a generated manifest
# instead, written next to the files.
set -uo pipefail
cd "$(git rev-parse --show-toplevel)" || { echo "cdf-corpus: not in a git checkout"; exit 1; }
REPO_ROOT="$(pwd)"

MANIFEST="$REPO_ROOT/space/cdf/corpus/manifest.tsv"
CORPUS="$REPO_ROOT/space/cdf/vendor/corpus"
ORACLE="$REPO_ROOT/space/cdf/vendor/install-ng"

bold() { printf '\n\033[1m[cdf-corpus] %s\033[0m\n' "$*"; }
info() { printf '[cdf-corpus] %s\n' "$*"; }
warn() { printf '\033[33m[cdf-corpus] %s\033[0m\n' "$*"; }
die()  { printf '\n\033[31m[cdf-corpus] FAILED: %s\033[0m\n' "$*" >&2; exit 1; }

[ -f "$MANIFEST" ] || die "missing $MANIFEST"

# Every data row of the manifest: tier, file, url, sha256, bytes, tags.
rows() { grep -v '^#' "$MANIFEST" | awk -F'\t' 'NF>=6 && $1 != "tier"'; }

# --- fetch ------------------------------------------------------------------------------------
do_fetch() {
    local want_tier="" ok=0 skipped=0
    [ "${1:-}" = "--tier" ] && { want_tier="tier${2:?--tier needs a number}"; }
    command -v curl >/dev/null || die "curl is required"

    bold "Fetching corpus${want_tier:+ (up to $want_tier)} …"
    while IFS=$'\t' read -r tier file url sha bytes _tags; do
        # Default to tiers 0+1: tier2 is 450 MB and only benchmarks need it.
        if [ -n "$want_tier" ]; then [ "$tier" \> "$want_tier" ] && continue
        else [ "$tier" = "tier2" ] && { skipped=$((skipped+1)); continue; }; fi

        local dest="$CORPUS/$tier/$file"
        mkdir -p "$CORPUS/$tier"
        if [ -f "$dest" ] && [ "$(sha256sum "$dest" | cut -d' ' -f1)" = "$sha" ]; then
            ok=$((ok+1)); continue
        fi
        info "$tier/$file ($(numfmt --to=iec "$bytes" 2>/dev/null || echo "$bytes B"))"
        curl -fsSL --max-time 3000 -o "$dest.part" "$url" || die "download failed: $url"
        mv -f "$dest.part" "$dest"

        local got; got="$(sha256sum "$dest" | cut -d' ' -f1)"
        if [ "$got" != "$sha" ]; then
            rm -f "$dest"
            printf '\n  expected %s\n  got      %s\n\n' "$sha" "$got"
            printf '  NASA reprocesses and re-versions archive files. If this file legitimately\n'
            printf '  changed upstream, re-run `scripts/cdf-corpus.sh probe` and commit the new row\n'
            printf '  as a REVIEWED change — a corpus that silently follows upstream cannot detect\n'
            printf '  a regression, because it has no fixed point to regress from.\n\n'
            die "$file does not match its pinned digest (deleted)"
        fi
        ok=$((ok+1))
    done < <(rows)
    info "$ok file(s) present and verified${skipped:+, $skipped tier2 file(s) skipped (use --tier 2)}"
}

# --- verify (offline) --------------------------------------------------------------------------
do_verify() {
    local ok=0 bad=0 absent=0
    while IFS=$'\t' read -r tier file _url sha _bytes _tags; do
        local dest="$CORPUS/$tier/$file"
        if [ ! -f "$dest" ]; then absent=$((absent+1)); continue; fi
        if [ "$(sha256sum "$dest" | cut -d' ' -f1)" = "$sha" ]; then ok=$((ok+1))
        else bad=$((bad+1)); warn "MISMATCH $tier/$file"; fi
    done < <(rows)
    info "verified: $ok ok, $bad mismatched, $absent not fetched"
    [ "$bad" -eq 0 ] || die "$bad corpus file(s) do not match the manifest"
}

# --- probe ------------------------------------------------------------------------------------
# Re-derives the tags column from the files themselves, so the manifest states measurements
# rather than claims. Prints the rows; diffing against the committed manifest is the point.
do_probe() {
    [ -x "$ORACLE/bin/cdfdump" ] || die "needs the oracle — run scripts/cdf-oracle.sh fetch && build"
    export CDF_LIB="$ORACLE/lib"
    # -dump nodata: metadata only. The full dump of a 250 MB file takes minutes and tells us
    # nothing extra — every tag below is metadata.
    local D="$ORACLE/bin/cdfdump -dump nodata"

    # Everything reads $dump from a here-string, never `printf ... | grep`. That is not style.
    # `grep -q` exits at the first match, which SIGPIPEs the writer, and with `set -o pipefail`
    # (line 25) the pipeline then reports 141 even though the match SUCCEEDED — so `&& tags+=`
    # silently does not run. It is a race on how much the writer flushed first, so it only bites
    # on large dumps: it silently dropped tags from exactly the four files whose metadata exceeds
    # ~50 KB. A here-string is not a pipeline and cannot do this.
    local dump=""
    has()  { grep -qE "$1" <<<"$dump"; }                       # regex present?
    field() { awk -F': *' "/^$1:/{print \$2; exit}" <<<"$dump"; }  # one header field

    while IFS=$'\t' read -r tier file url sha bytes recorded_tags; do
        local f="$CORPUS/$tier/$file"
        if [ ! -f "$f" ]; then
            warn "not fetched, keeping recorded tags: $tier/$file"
            printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$tier" "$file" "$url" "$sha" "$bytes" "$recorded_tags"
            continue
        fi

        local magic tags=()
        magic=$(xxd -p -l4 "$f")
        dump=$($D "${f%.cdf}" 2>/dev/null)
        [ -n "$dump" ] || die "cdfdump produced nothing for $tier/$file — the oracle is broken, and a silent empty dump would look like a file with no features"

        case "$magic" in
            cdf30001) tags+=(v3) ;;
            cdf26002) tags+=(v2.6) ;;
            0000ffff) tags+=(pre-v2.6) ;;
            *)        tags+=("magic-$magic") ;;
        esac
        tags+=("$(field Encoding | tr 'A-Z' 'a-z')")
        tags+=("$(field Majority | tr 'A-Z' 'a-z')")
        [ "$(field Format)" = "MULTI" ] && tags+=(multifile)
        [ "$(field Checksum)" = "MD5" ] && tags+=(checksum)

        local fc; fc=$(field Compression)   # the FILE-level (CCR) compression
        [ -n "$fc" ] && [ "$fc" != "None" ] && tags+=("ccr-$(tr 'A-Z.' 'a-z-' <<<"$fc")")

        local c
        for c in RLE GZIP HUFF AHUFF; do
            has "\(Compression: $c" && tags+=("var-$(tr 'A-Z' 'a-z' <<<"$c")")
        done

        has 'CDF_EPOCH16'   && tags+=(epoch16)
        has 'CDF_EPOCH[^1]' && tags+=(epoch)
        has 'CDF_TT2000'    && tags+=(tt2000)
        has 'CDF_U?CHAR'    && tags+=(char)
        has 'sRecords\.PAD'  && tags+=(sparse-pad)
        has 'sRecords\.PREV' && tags+=(sparse-prev)

        [ "$(field NumrVars)" -gt 0 ] 2>/dev/null && tags+=(rvars)
        [ "$(field NumzVars)" -gt 0 ] 2>/dev/null && tags+=(zvars)

        printf '%s\t%s\t%s\t%s\t%s\t%s\n' \
            "$tier" "$file" "$url" "$sha" "$(stat -c%s "$f")" "$(IFS=,; printf '%s' "${tags[*]}")"
    done < <(rows)
}

# --- derive -------------------------------------------------------------------------------------
# The covering matrix. Each cell is one cdfconvert re-encode of a known-good source, so what must
# come back out is known exactly. This is the ONLY coverage RLE/HUFF/AHUFF/CCR will ever get.
do_derive() {
    [ -x "$ORACLE/bin/cdfconvert" ] || die "needs the oracle — run scripts/cdf-oracle.sh fetch && build"
    export CDF_LIB="$ORACLE/lib"
    local CVT="$ORACLE/bin/cdfconvert" SRC="$CORPUS/tier0/test_alltypes"
    [ -f "$SRC.cdf" ] || die "tier0 not fetched — run scripts/cdf-corpus.sh fetch"

    local OUT="$CORPUS/tier3"
    mkdir -p "$OUT"; rm -f "$OUT"/*.cdf "$OUT/manifest.tsv"
    bold "Generating tier3 …"

    local made=0 failed=0
    emit() { # emit <name> <cdfconvert args...>
        local name="$1"; shift
        if "$CVT" "$SRC" "$OUT/$name" "$@" -delete >/dev/null 2>&1 && [ -f "$OUT/$name.cdf" ]; then
            printf '%s\t%s\t%s\n' "$name.cdf" "$(sha256sum "$OUT/$name.cdf" | cut -d' ' -f1)" "$*" >> "$OUT/manifest.tsv"
            made=$((made+1))
        else
            warn "could not generate $name ($*)"; failed=$((failed+1))
        fi
    }

    # Per-variable compression x blocking factor — the four codecs, which is the whole point.
    local c bf
    for c in none rle.0 huff.0 ahuff.0 gzip.1 gzip.6 gzip.9; do
        for bf in 1 8; do emit "vc_${c%%.*}${c#*.}_bf$bf" -compression "vars:$c:$bf"; done
    done
    # Whole-file (CCR) compression, with and without per-variable compression underneath.
    for c in rle.0 huff.0 ahuff.0 gzip.6; do
        emit "ccr_${c%%.*}"      -compression "cdf:$c,vars:none"
        emit "ccr_${c%%.*}_vgz"  -compression "cdf:$c,vars:gzip.6:1"
    done
    # Majority, sparseness, and encoding. Encodings are attempted, not assumed: NASA's Linux
    # build may decline to write some of them, and a cell that silently vanishes would look
    # like coverage we do not have — so failures are reported and recorded as gaps.
    emit row    -row
    emit column -column
    # Sparseness is set PER VARIABLE, not with `vars:`. Asking for it on every variable fails
    # outright (CANNOT_SPARSERECORDS) because some variables cannot carry sparse records at all,
    # and one refusal aborts the whole conversion. "Temp" is record-variant and uncompressed,
    # which is what the setting needs.
    local s e
    for s in no pad prev; do emit "sparse_$s" -sparseness "var:\"Temp\":srecords.$s"; done

    # Encodings. `mac` is deliberately absent: this build rejects it ("Unknown encoding"), and
    # nothing is lost — MAC_ENCODING is big-endian IEEE, byte-identical to the `sun`/`network`
    # cells that are here. Every VAX variant IS produced, which matters more than it looks:
    # VAX F/D/G floats are not IEEE-754, so these are the only files that exercise that
    # conversion end to end rather than through hand-transcribed unit vectors.
    for e in network ibmpc sun decstation sgi ibmrs hp next alphaosf1 alphavmsd alphavmsg alphavmsi vax; do
        emit "enc_$e" -encoding "$e"
    done

    bold "tier3: $made file(s) generated${failed:+, $failed cell(s) unavailable}"
    info "manifest: $OUT/manifest.tsv"
    if ! grep -q '^enc_vax' "$OUT/manifest.tsv" 2>/dev/null; then
        warn "NAMED GAP: no VAX-encoded file could be produced, and none exists in the archive."
        warn "VAX F/D/G float decoding is therefore covered only by the IFD's own conversion"
        warn "tables (Appendices A/B) as unit vectors — not by any end-to-end differential test."
    fi
}

do_status() {
    local have=0 total=0 bytes=0
    while IFS=$'\t' read -r tier file _url sha bytes_ _tags; do
        total=$((total+1))
        [ -f "$CORPUS/$tier/$file" ] && { have=$((have+1)); bytes=$((bytes+bytes_)); }
    done < <(rows)
    printf 'manifest : %s rows\n' "$total"
    printf 'fetched  : %s (%s)\n' "$have" "$(numfmt --to=iec "$bytes" 2>/dev/null || echo "$bytes B")"
    local t3; t3=$(ls "$CORPUS/tier3"/*.cdf 2>/dev/null | wc -l)
    printf 'tier3    : %s generated file(s)\n' "$t3"
}

do_clean() {
    if [ "${1:-}" = "--all" ]; then bold "Removing the whole corpus …"; rm -rf "$CORPUS"
    else bold "Removing generated tier3 only …"; rm -rf "$CORPUS/tier3"; fi
    info "done."
}

case "${1:-status}" in
    fetch)  shift; do_fetch "$@" ;;
    verify) do_verify ;;
    probe)  do_probe ;;
    derive) do_derive ;;
    status) do_status ;;
    clean)  do_clean "${2:-}" ;;
    *) die "unknown subcommand '$1' — one of: fetch, verify, probe, derive, status, clean" ;;
esac
