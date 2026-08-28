#!/usr/bin/env bash
# The NASA CDF library, fetched and built as cheatah-space's DEV-ONLY differential oracle.
#
#   scripts/cdf-oracle.sh fetch           # download + verify the pinned tarball, extract it
#   scripts/cdf-oracle.sh build           # build both zlib flavours, install into vendor/
#   scripts/cdf-oracle.sh status          # what is present, and whether it matches the pin
#   scripts/cdf-oracle.sh check-upstream  # has NASA published a newer release? (network, opt-in)
#   scripts/cdf-oracle.sh clean           # remove everything except the verified tarball
#   scripts/cdf-oracle.sh clean --all     # remove the tarball too
#
# WHY THIS EXISTS. space.cdf is written from the published CDF Internal Format Description, not
# from NASA's source. But "we decode this file the same way NASA does" is a claim that has to be
# MEASURED, and "we are faster" is a claim that needs something to be faster THAN. This script
# produces that reference: cdfdump, cdfconvert, cdfcompare, cdfvalidate, cdfirsdump, and
# libcdf.a, all under space/cdf/vendor/.
#
# WHY FROM SOURCE, not the prebuilt binaries. Three reasons, each load-bearing:
#   1. The prebuilt Linux tarball is glibc-coupled and may not run here at all.
#   2. Benchmarking someone else's optimization flags against ours measures the flags, not the
#      code. Building it ourselves puts both sides on -O2 -DNDEBUG with the same compiler.
#   3. CDF 3.9.2 bundles BOTH classic zlib and zlib-ng (-DZLIB=orig|ng, defaulting to orig).
#      Our from-scratch inflate has to beat zlib-ng or the comparison is a strawman, so we build
#      both and report against each.
#
# NOTHING THIS SCRIPT DOWNLOADS IS EVER COMMITTED. space/cdf/vendor/ is .gitignore'd and
# scripts/check_no_vendored_nasa.sh fails the QA gate if anything from it is ever staged.
# See scripts/cdf-oracle.pins for the licensing terms — the CDF distribution is NOT public domain.
#
# The QA gate NEVER runs this. It needs network on first use; everything after is offline and
# idempotent, short-circuited on a stamp recording (url, sha256, compiler, flags).
set -uo pipefail
cd "$(git rev-parse --show-toplevel)" || { echo "cdf-oracle: not in a git checkout"; exit 1; }
REPO_ROOT="$(pwd)"

PINS="$REPO_ROOT/scripts/cdf-oracle.pins"
VENDOR="$REPO_ROOT/space/cdf/vendor"
STAMP="$VENDOR/.oracle-stamp"

bold() { printf '\n\033[1m[cdf-oracle] %s\033[0m\n' "$*"; }
info() { printf '[cdf-oracle] %s\n' "$*"; }
warn() { printf '\033[33m[cdf-oracle] %s\033[0m\n' "$*"; }
die()  { printf '\n\033[31m[cdf-oracle] FAILED: %s\033[0m\n' "$*" >&2; exit 1; }

[ -f "$PINS" ] || die "missing $PINS"
# shellcheck disable=SC1090
. "$PINS"
for v in CDF_ORACLE_VERSION CDF_ORACLE_URL CDF_ORACLE_SHA256 CDF_ORACLE_MD5 CDF_ORACLE_TARROOT; do
    [ -n "${!v:-}" ] || die "$PINS does not define $v"
done

TARBALL="$VENDOR/$(basename "$CDF_ORACLE_URL")"
SRCDIR="$VENDOR/$CDF_ORACLE_TARROOT"

need() { command -v "$1" >/dev/null 2>&1 || die "$1 is required but not on PATH"; }

# The identity of a built oracle: change any of these and the build is stale. The compiler
# version is in here deliberately — a benchmark number produced by a different compiler is a
# different number, and silently reusing the old tree would hide that.
stamp_now() {
    printf 'url=%s\nsha256=%s\ncc=%s\ncflags=%s\n' \
        "$CDF_ORACLE_URL" "$CDF_ORACLE_SHA256" \
        "$(${CC:-cc} --version 2>/dev/null | head -1)" "-O2 -DNDEBUG"
}

# --- fetch ------------------------------------------------------------------------------------
# Verify first, extract second. A tarball that fails its digest is never unpacked.
do_fetch() {
    need curl; need sha256sum; need tar
    mkdir -p "$VENDOR"

    if [ -f "$TARBALL" ] && [ "$(sha256sum "$TARBALL" | cut -d' ' -f1)" = "$CDF_ORACLE_SHA256" ]; then
        info "tarball already present and matches the pin — not refetching."
    else
        bold "Fetching $CDF_ORACLE_VERSION …"
        info "$CDF_ORACLE_URL"
        curl -fsSL --max-time 900 -o "$TARBALL.part" "$CDF_ORACLE_URL" \
            || die "download failed (no network? check-upstream needs one too)"
        mv -f "$TARBALL.part" "$TARBALL"
        verify_or_explain
    fi

    bold "Extracting …"
    rm -rf "$SRCDIR"
    tar xzf "$TARBALL" -C "$VENDOR" || die "extract failed"
    [ -d "$SRCDIR" ] || die "expected $CDF_ORACLE_TARROOT/ inside the tarball; the pin's TARROOT is wrong"
    info "source tree: $SRCDIR"

    _patch_upstream
}

# The whole point of the pin. On a digest mismatch, say WHICH failure this is: NASA re-rolling a
# release in place (expected, needs a reviewed pin bump) or bytes that are simply wrong
# (corruption, a proxy, an attacker). Never auto-accept either.
verify_or_explain() {
    local got_sha got_md5 up_md5
    got_sha="$(sha256sum "$TARBALL" | cut -d' ' -f1)"
    [ "$got_sha" = "$CDF_ORACLE_SHA256" ] && { info "sha256 OK — matches the pin."; return 0; }

    got_md5="$(md5sum "$TARBALL" | cut -d' ' -f1)"
    up_md5="$(curl -fsSL --max-time 120 "$CDF_ORACLE_URL.md5" 2>/dev/null | awk '{print $1}')"

    printf '\n\033[31m[cdf-oracle] PINNED DIGEST MISMATCH\033[0m\n\n'
    printf '  expected sha256 : %s\n  got sha256      : %s\n\n' "$CDF_ORACLE_SHA256" "$got_sha"

    if [ -n "$up_md5" ] && [ "$up_md5" = "$CDF_ORACLE_MD5" ]; then
        printf '  NASA still publishes MD5 %s, which is what the pin records,\n' "$up_md5"
        printf '  but the bytes we received hash differently (md5 %s).\n\n' "$got_md5"
        rm -f "$TARBALL"
        die "the artifact does not match what NASA says it is publishing — corruption or interception. Downloaded file deleted; NOT bumping the pin."
    fi

    if [ -n "$up_md5" ]; then
        printf '  NASA now publishes MD5 %s (the pin records %s),\n' "$up_md5" "$CDF_ORACLE_MD5"
        printf '  and the bytes we received hash to md5 %s — consistent with an upstream re-roll.\n\n' "$got_md5"
    else
        printf '  Could not fetch %s.md5 to tell a re-roll from corruption.\n\n' "$CDF_ORACLE_URL"
    fi

    printf '  Review the change, then update scripts/cdf-oracle.pins:\n\n'
    printf '    CDF_ORACLE_SHA256=%s\n    CDF_ORACLE_MD5=%s\n\n' "$got_sha" "${up_md5:-$got_md5}"
    printf '  Then re-run. Bumping a @perf tag is likely too — scripts/perf_tags.sh names the\n'
    printf '  pinned oracle version, so measurements do not silently outlive their reference.\n\n'
    die "pinned digest mismatch — the pin is REVIEWED, never auto-accepted"
}

# --- patch: two upstream build bugs -----------------------------------------------------------
# Applied to the EXTRACTED tree only, after every fetch, idempotently. Both are genuine defects
# in NASA's shipped CMake for the C-only distribution (cdf39_2-dist-cdf.tar.gz).
#
#   BUG 1  liblib/CMakeLists.txt references the SHARED target `cdf` unconditionally
#          (TARGET_LINK_LIBRARIES(cdf ${CMAKE_DL_LIBS}), and an MSVC set_target_properties),
#          outside any if(BUILD_SHARED_CDF). Configuring with -DBUILD_SHARED_CDF=OFF therefore
#          dies with "Cannot specify link libraries for target cdf which is not built by this
#          project". Fix: guard both with if(TARGET cdf).
#
#   BUG 2  The if(MSVC) guard around the Java jar install rules is COMMENTED OUT ("#if(MSVC)" /
#          "#endif()"), so `cmake --install` unconditionally tries to install cdfjava.jar and
#          four siblings. Those jars ship only in the -java tarball, so installing the C-only
#          distribution fails outright. Fix: restore a real guard, keyed on the jar existing.
#
#   BUG 3  Same shape as BUG 1, in the BUNDLED zlib-ng. NASA wrapped its
#          add_library(zlib-ng OBJECT ...) in if(BUILD_SHARED_CDF), but the block that sets
#          properties on that target is still guarded only by upstream zlib-ng's own
#          "if(NOT DEFINED BUILD_SHARED_LIBS OR BUILD_SHARED_LIBS)" — which is TRUE when
#          BUILD_SHARED_LIBS is simply unset. With -DBUILD_SHARED_CDF=OFF the target does not
#          exist yet five set_target_properties calls run against it. Fix: require the target.
#
# CRITICAL: this patches BUILD PLUMBING ONLY. Not one line of NASA's C source is touched, ever.
# The whole value of this tree is that it is an untouched reference implementation — patching
# the library would corrupt the very oracle we measure against. If a future pin needs a source
# change to build, that is a STOP-and-think, not a patch to add here.
#
# Nothing patched here is ever redistributed: the tree lives in the .gitignore'd vendor dir.
_patch_upstream() {
    local f="$SRCDIR/liblib/CMakeLists.txt"
    [ -f "$f" ] || die "expected $f — the pin's TARROOT or layout is wrong"

    if grep -q 'cheatah-space: patched' "$f"; then
        info "upstream build patches already applied."
        return 0
    fi

    # BUG 1 — guard the two unconditional references to the shared `cdf` target.
    sed -i \
      -e 's|^      TARGET_LINK_LIBRARIES(cdf \${CMAKE_DL_LIBS})$|      if(TARGET cdf) # cheatah-space: patched (BUG 1)\n        TARGET_LINK_LIBRARIES(cdf ${CMAKE_DL_LIBS})\n      endif()|' \
      -e 's|^  set_target_properties(cdf PROPERTIES$|  if(TARGET cdf) # cheatah-space: patched (BUG 1)\n  set_target_properties(cdf PROPERTIES|' \
      "$f" || die "BUG 1 patch failed"

    # Close the if(TARGET cdf) opened around the MSVC set_target_properties block.
    awk '
      /# cheatah-space: patched \(BUG 1\)/ && /if\(TARGET cdf\)/ { inblk = 1 }
      { print }
      inblk && /^  \)$/ { print "  endif()"; inblk = 0 }
    ' "$f" > "$f.tmp" && mv "$f.tmp" "$f"

    # BUG 2 — restore the commented-out guard, keyed on the jar actually being present.
    awk '
      $0 == "#if(MSVC)" && !done {
          print "if(EXISTS ${CDF_SOURCE_DIR}/cdfjava/classes/cdfjava.jar) # cheatah-space: patched (BUG 2)"
          injar = 1; next
      }
      injar && $0 == "#endif()" { print "endif()"; injar = 0; done = 1; next }
      { print }
    ' "$f" > "$f.tmp" && mv "$f.tmp" "$f"

    # BUG 3 — the bundled zlib-ng sets properties on a target that only exists for shared builds.
    local z="$SRCDIR/src/lib/zlib-ng/CMakeLists.txt"
    [ -f "$z" ] || die "expected $z — the pin's layout changed"
    sed -i 's|^if(NOT DEFINED BUILD_SHARED_LIBS OR BUILD_SHARED_LIBS)$|if((NOT DEFINED BUILD_SHARED_LIBS OR BUILD_SHARED_LIBS) AND TARGET zlib-ng) # cheatah-space: patched (BUG 3)|' "$z" \
        || die "BUG 3 patch failed"
    grep -q 'cheatah-space: patched (BUG 3)' "$z" || die "BUG 3 patch did not apply — upstream zlib-ng layout changed"

    grep -q 'cheatah-space: patched (BUG 1)' "$f" || die "BUG 1 patch did not apply — upstream layout changed; re-check liblib/CMakeLists.txt against the new pin"
    grep -q 'cheatah-space: patched (BUG 2)' "$f" || die "BUG 2 patch did not apply — upstream layout changed; re-check liblib/CMakeLists.txt against the new pin"
    info "applied 3 upstream build patches (build plumbing only; C source untouched)."
}

# --- build ------------------------------------------------------------------------------------
# Both zlib flavours, installed side by side. The `ng` install is the one benchmarks quote.
do_build() {
    need cmake
    [ -d "$SRCDIR" ] || die "no source tree — run 'scripts/cdf-oracle.sh fetch' first"

    if [ -f "$STAMP" ] && [ "$(cat "$STAMP")" = "$(stamp_now)" ]; then
        info "already built from this pin with this compiler — nothing to do."
        info "(force a rebuild with 'scripts/cdf-oracle.sh clean')"
        return 0
    fi

    # Static only. We link libcdf.a, and _patch_upstream fixed the upstream defect that made
    # -DBUILD_SHARED_CDF=OFF fail, so there is no reason to build a shared object nobody loads.
    local flavour
    for flavour in ng orig; do
        bold "Building NASA CDF $CDF_ORACLE_VERSION (zlib=$flavour) …"
        cmake -S "$SRCDIR" -B "$VENDOR/build-$flavour" \
              -DCMAKE_BUILD_TYPE=Release \
              -DCMAKE_C_FLAGS="-O2 -DNDEBUG" \
              -DZLIB="$flavour" \
              -DBUILD_SHARED_CDF=OFF -DBUILD_STATIC_CDF=ON \
              -DBUILD_UTILITIES=ON -DBUILD_TESTS=OFF \
              -DCMAKE_INSTALL_PREFIX="$VENDOR/install-$flavour" \
              >"$VENDOR/build-$flavour.log" 2>&1 \
            || { tail -30 "$VENDOR/build-$flavour.log"; die "configure failed (zlib=$flavour) — see $VENDOR/build-$flavour.log"; }
        cmake --build "$VENDOR/build-$flavour" --parallel "$(nproc)" \
              >>"$VENDOR/build-$flavour.log" 2>&1 \
            || { tail -40 "$VENDOR/build-$flavour.log"; die "build failed (zlib=$flavour) — see $VENDOR/build-$flavour.log"; }
        cmake --install "$VENDOR/build-$flavour" >>"$VENDOR/build-$flavour.log" 2>&1 \
            || { tail -30 "$VENDOR/build-$flavour.log"; die "install failed (zlib=$flavour)"; }
        info "installed -> $VENDOR/install-$flavour"
    done

    # Assert the tools we actually depend on exist. A build that "succeeded" without cdfirsdump
    # would fail much later, in the differential harness, with a far worse error message.
    local missing=() t
    for t in cdfdump cdfconvert cdfcompare cdfvalidate cdfirsdump; do
        [ -x "$VENDOR/install-ng/bin/$t" ] || missing+=("$t")
    done
    [ ${#missing[@]} -eq 0 ] || die "built, but these expected tools are missing: ${missing[*]}"

    stamp_now > "$STAMP"
    bold "Oracle ready."
    "$VENDOR/install-ng/bin/cdfdump" -help >/dev/null 2>&1 || true
    info "tools: $VENDOR/install-ng/bin"
    info "lib:   $(find "$VENDOR/install-ng" -name 'libcdf.a' 2>/dev/null | head -1)"
}

# --- status -----------------------------------------------------------------------------------
do_status() {
    printf 'pin      : %s\n' "$CDF_ORACLE_VERSION"
    printf 'url      : %s\n' "$CDF_ORACLE_URL"
    printf 'sha256   : %s\n' "$CDF_ORACLE_SHA256"
    if [ -f "$TARBALL" ]; then
        local got; got="$(sha256sum "$TARBALL" | cut -d' ' -f1)"
        if [ "$got" = "$CDF_ORACLE_SHA256" ]; then printf 'tarball  : present, matches the pin\n'
        else printf 'tarball  : present, \033[31mDOES NOT MATCH\033[0m (%s)\n' "$got"; fi
    else
        printf 'tarball  : absent\n'
    fi
    [ -d "$SRCDIR" ] && printf 'source   : %s\n' "$SRCDIR" || printf 'source   : not extracted\n'
    local f
    for f in ng orig; do
        if [ -x "$VENDOR/install-$f/bin/cdfdump" ]; then printf 'zlib=%-5s: built\n' "$f"
        else printf 'zlib=%-5s: not built\n' "$f"; fi
    done
    if [ -f "$STAMP" ] && [ "$(cat "$STAMP")" = "$(stamp_now)" ]; then
        printf 'stamp    : current\n'
    elif [ -f "$STAMP" ]; then
        printf 'stamp    : \033[33mstale\033[0m (pin or compiler changed since the build)\n'
    else
        printf 'stamp    : none\n'
    fi
}

# --- check-upstream ---------------------------------------------------------------------------
# Opt-in, human- or CI-run, never part of the gate. Answers one question: has NASA published a
# release directory newer than the one we pin? Finding out from a weekly job beats finding out
# from a reviewer asking why our numbers cite a two-year-old reference.
do_check_upstream() {
    need curl
    bold "Checking for newer CDF releases …"
    local listing newer
    listing="$(curl -fsSL --max-time 120 https://spdf.gsfc.nasa.gov/pub/software/cdf/dist/ 2>/dev/null)" \
        || die "could not reach the distribution index"
    newer="$(printf '%s' "$listing" | grep -oE 'href="cdf[0-9]+_[0-9]+/"' | sed 's/href="//;s/\/"//' \
             | sort -V | awk -v pin="$CDF_ORACLE_VERSION" '$0 > pin')"
    if [ -z "$newer" ]; then
        info "none — $CDF_ORACLE_VERSION is still the newest published release."
        return 0
    fi
    warn "newer release directories exist:"
    printf '%s\n' "$newer" | sed 's/^/    /'
    warn "Bumping is a REVIEWED change: update scripts/cdf-oracle.pins, re-run fetch + build,"
    warn "then re-measure — scripts/perf_tags.sh requires every @perf tag to name the pinned version."
}

do_clean() {
    if [ "${1:-}" = "--all" ]; then
        bold "Removing the entire vendor oracle tree …"; rm -rf "$VENDOR"
    else
        bold "Removing build + install trees (keeping the verified tarball) …"
        rm -rf "$SRCDIR" "$VENDOR"/build-* "$VENDOR"/install-* "$VENDOR"/build-*.log "$STAMP"
    fi
    info "done."
}

case "${1:-status}" in
    fetch)          do_fetch ;;
    build)          do_build ;;
    status)         do_status ;;
    check-upstream) do_check_upstream ;;
    clean)          do_clean "${2:-}" ;;
    *) die "unknown subcommand '${1}' — one of: fetch, build, status, check-upstream, clean" ;;
esac
