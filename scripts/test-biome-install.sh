#!/usr/bin/env bash
# test-biome-install.sh — sandbox the EXACT experience of someone with a standard cheatah install who
# runs `biome add cheatah-space`. Nothing from this working tree (build/, .git, dev env vars) is allowed
# to make it falsely pass: we copy ONLY the installable package to a throwaway location (as biome fetches
# it), compile a fresh user project against it with cheatah env vars CLEARED, and the only wiring is the
# extension on CHEATAH_MODULE_PATH — precisely what biome's EXTENSIONS support sets.
#
# If a module header's .sha512 sidecar, the module layout, or the namespaces are wrong, this fails
# exactly as a real user's install would. cheatah-space is PURE cheatah — no GPU, no SDK, no system
# dependency — so there is nothing else to wire: the clean env + module path IS the whole install.
set -uo pipefail
cd "$(git rev-parse --show-toplevel)"
CHEATAH_DIR="${CHEATAH_DIR:-$PWD/../cheatah}"

# The user's installed toolchain (purrc/cheatah). Its stdlib root is BAKED into the binary, so io/math/
# ndarray resolve without any CHEATAH_ROOT env — just like a real install.
find_tool() {
    local n="$1"
    for c in release debug asan; do
        [ -x "$CHEATAH_DIR/build/$c/bin/$n" ] && { echo "$CHEATAH_DIR/build/$c/bin/$n"; return 0; }
    done
    command -v "$n" 2>/dev/null
}
PURRC="$(find_tool purrc)"; CHEATAH="$(find_tool cheatah)"
[ -x "$PURRC" ] && [ -x "$CHEATAH" ] || { echo "biome-install: no cheatah toolchain (set CHEATAH_DIR)"; exit 2; }

# 1. Simulate biome fetching cheatah-space to a fresh dir OUTSIDE this tree. A consumer gets the
#    `space/` package (headers + the .sha512 sidecars) and nothing else — copy exactly that.
INSTALL="$(mktemp -d)"; PROJ="$(mktemp -d)"
trap 'rm -rf "$INSTALL" "$PROJ"' EXIT
cp -r space "$INSTALL/space"
[ -f "$INSTALL/space/space.hpp.sha512" ] || { echo "biome-install: the fetched copy has no space/space.hpp.sha512 — run scripts/sign-modules.sh"; exit 1; }
[ -f "$INSTALL/space/time/time.hpp.sha512" ] || { echo "biome-install: the fetched copy has no space/time/time.hpp.sha512 — run scripts/sign-modules.sh"; exit 1; }

# 2. A brand-new user project that just imports the module and round-trips a Julian Date.
mkdir -p "$PROJ/src"
cat > "$PROJ/src/main.purr" <<'PURR'
# Seconds after `biome add cheatah-space`. This user never saw a git repo.
import io
import space.time as time

# J2000: 2000-01-01T12:00:00 is Unix 946728000 and, by definition, JD 2451545.0.
let jd = time.unix_to_jd(946728000.0)
let back = time.jd_to_unix(jd)
io.print("space.time maps J2000 to JD", jd, "and back to", back)
if jd == time.jd_j2000() {
    if back == 946728000.0 {
        io.print("RESULT: PASS")
    } else {
        io.print("RESULT: FAIL", back)
    }
} else {
    io.print("RESULT: FAIL", jd)
}
PURR

# 3. Compile from the user's project with a CLEAN environment: every cheatah env var cleared, and the
#    ONLY wiring is CHEATAH_MODULE_PATH -> the fetched package (what cheatah_add_program EXTENSIONS
#    does). cheatah-space is pure cheatah, so no SDK path (no CPATH, no Vulkan) is ever needed.
clean_env=(env -u CHEATAH_ROOT -u CHEATAH_LIB_DIR -u CHEATAH_TRUST -u CHEATAH_DIR
           CHEATAH_MODULE_PATH="$INSTALL")
if ! out="$(cd "$PROJ" && "${clean_env[@]}" "$PURRC" src/main.purr -o app.so 2>&1)"; then
    echo "biome-install: FAILED — a fresh user could not compile 'import space.time':"
    echo "$out" | sed 's/^/    /'
    exit 1
fi
run="$("$CHEATAH" "$PROJ/app.so" 2>&1)"
echo "$run" | sed 's/^/    /'
echo "$run" | grep -q "RESULT: PASS" \
    || { echo "biome-install: FAILED — the user program did not pass"; exit 1; }
echo "biome-install: PASS — fresh project, package fetched to a throwaway dir on CHEATAH_MODULE_PATH, clean env."
