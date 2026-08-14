#!/usr/bin/env bash
# sign-modules.sh — write the .sha512 sidecars that mark cheatah-space's module headers as VERIFIED
# cheatah modules. This is what makes the extension biome-installable: with the sidecar present, purrc
# resolves `import space` / `import space.*` on the extension path (CHEATAH_MODULE_PATH, which
# `biome add cheatah-space` sets) and the runtime verifies each header against its checksum on load — so
# a user with a standard cheatah install never touches git or --import-root.
#
# The sidecar is sha512sum format ("<hex>  <basename>\n"): identical to what `purrc --emit-library`
# emits for generated module headers (compiler/purrc.cpp writes hex + "  " + base_name), and what
# `sha512sum -c` validates. The runtime verifier reads only the first whitespace token, so the filename
# suffix is cosmetic — but keeping it matches purrc's own output and stays `-c`-checkable.
#
# EVERY space module header is hand-authored C++ (purrc never emits any of them — that is this
# extension's design: concept-templated C++20 the transpiler could not express), so ALL of them are
# signed here: the space/space.hpp umbrella and each space/<mod>/<mod>.hpp submodule. Re-run whenever a
# header changes; the QA gate checks the sidecars are in sync.
set -euo pipefail
cd "$(git rev-parse --show-toplevel)"

shopt -s nullglob
for h in space/space.hpp space/*/*.hpp; do
    [ -f "$h" ] || continue
    # Sign from the header's own directory so the sidecar records just the basename (matching purrc).
    ( cd "$(dirname "$h")" && sha512sum "$(basename "$h")" > "$(basename "$h").sha512" )
    echo "signed: $h.sha512"
done
