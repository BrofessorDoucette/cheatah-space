#!/usr/bin/env bash
# scripts/bench_run.sh — build, run, publish and gate space.irbem's benchmark suite.
#
# Every performance number in space/irbem/BENCHMARKS.md, and every performance number quoted in a
# header comment, comes from here. Nothing in the published tables is hand-typed: `publish`
# regenerates the region between a suite's `<!-- BENCH:<suite> begin -->` and `end` markers from the
# JSON this run produced, stamps it with the host, the commit and the exact command that produced
# it, and leaves the prose around it alone.
#
# Usage:
#   scripts/bench_run.sh build                 compile Google Benchmark (pinned) and the suite
#   scripts/bench_run.sh run <suite>           run one suite, write docs/bench/<suite>.json
#   scripts/bench_run.sh publish <suite>|all   run and regenerate that suite's table
#   scripts/bench_run.sh baseline              run every suite and rewrite bench/baseline.csv
#   scripts/bench_run.sh check                 run every suite and fail on a regression
#   scripts/bench_run.sh manifest              print the routine manifest (no timing)
#   scripts/bench_run.sh intensity             print the arithmetic-intensity model (no timing)
#
# Suites:
#   irbem-per-routine   every routine, on every lane it has, against the IRBEM -O2 oracle
#   irbem-crossover     trace throughput against batch size — the module's central argument
#   irbem-intensity     arithmetic intensity against measured verdict
#
# Environment:
#   CHEATAH_DIR                 cheatah checkout             (default ../cheatah)
#   CHEATAH_GPU_LINALG_DIR      cheatah-gpu-linalg checkout  (default ../cheatah-gpu-linalg)
#   SPACE_IRBEM_ORACLE          the -O2 IRBEM build          (default /tmp/irbem-builds/libirbem-O2.so)
#   CHEATAH_GPU_LINALG_VK_DEVICE  device name substring, e.g. NVIDIA / Intel / llvmpipe
#   BENCH_REPS                  repetitions                  (default 5)
#   BENCH_MIN_TIME              per-repetition minimum       (default 0.2s)
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
WORK="${CHEATAH_SPACE_BENCH_DIR:-$ROOT/build/bench}"
BIN="$WORK/irbem_bench"
JSON_DIR="$ROOT/docs/bench"
PAGE="$ROOT/space/irbem/BENCHMARKS.md"
BASELINE="$ROOT/bench/baseline.csv"

CHEATAH_DIR="${CHEATAH_DIR:-$ROOT/../cheatah}"
CHEATAH_GPU_LINALG_DIR="${CHEATAH_GPU_LINALG_DIR:-$ROOT/../cheatah-gpu-linalg}"
SPACE_IRBEM_ORACLE="${SPACE_IRBEM_ORACLE:-/tmp/irbem-builds/libirbem-O2.so}"
export SPACE_IRBEM_ORACLE
BENCH_REPS="${BENCH_REPS:-5}"
BENCH_MIN_TIME="${BENCH_MIN_TIME:-0.2s}"

# The pin. Same commit as cheatah-gpu-linalg/bench/CMakeLists.txt and cheatah's tests/benchmarks:
# a harness that drifts between repos produces numbers that cannot be compared between repos.
BENCHMARK_PIN=192ef10025eb2c4cdd392bc502f0c852196baa48

# -O3 -march=native, and -ffp-contract=off so the host lane's arithmetic is the arithmetic the
# golden tests pinned. NEVER -ffast-math — it would let the compiler reassociate the spherical
# harmonic sums and the RK4 accumulation, and the binary would then be measuring different
# mathematics from the one the correctness suite verified. bench/CMakeLists.txt carries the same
# list; `check` compares the two so they cannot drift apart silently.
BENCH_CXXFLAGS=(-std=c++20 -O3 -march=native -ffp-contract=off)

die() { printf 'bench_run: %s\n' "$*" >&2; exit 1; }
note() { printf '  %s\n' "$*" >&2; }

# --- the suites -------------------------------------------------------------------------------
# name : benchmark_filter : layout
suite_filter() {
    case "$1" in
        irbem-per-routine) echo '^BM_' ;;
        irbem-crossover)   echo 'trace_batch' ;;
        irbem-intensity)   echo '(dipole_field_(host_)?batch|igrf_batch|trace_batch/65536)' ;;
        *) die "unknown suite '$1' (irbem-per-routine, irbem-crossover, irbem-intensity)" ;;
    esac
}
ALL_SUITES=(irbem-per-routine irbem-crossover irbem-intensity)

# =============================================================================================
# build
# =============================================================================================

benchmark_src() {
    # Reuse a sibling checkout at the pin before cloning one: cheatah-gpu-linalg already fetches
    # this exact commit, and a second clone of it is a minute of network for no new information.
    local candidate
    for candidate in \
        "$WORK/benchmark-src" \
        "$CHEATAH_GPU_LINALG_DIR/build/_deps/benchmark-src" \
        "$CHEATAH_DIR/build/release/_deps/benchmark-src"
    do
        if [[ -f "$candidate/src/benchmark.cc" ]] &&
           [[ "$(git -C "$candidate" rev-parse HEAD 2>/dev/null || true)" == "$BENCHMARK_PIN" ]]; then
            echo "$candidate"; return 0
        fi
    done
    mkdir -p "$WORK"
    note "fetching google/benchmark at $BENCHMARK_PIN"
    rm -rf "$WORK/benchmark-src"
    git clone --quiet https://github.com/google/benchmark.git "$WORK/benchmark-src"
    git -C "$WORK/benchmark-src" checkout --quiet "$BENCHMARK_PIN"
    echo "$WORK/benchmark-src"
}

build_benchmark_lib() {
    local src="$1" out="$WORK/libbenchmark.a"
    if [[ -f "$out" && "$out" -nt "$src/src/benchmark.cc" ]]; then echo "$out"; return 0; fi
    note "building Google Benchmark (release; the vendored debug builds warn and mis-time)"
    local objdir="$WORK/benchmark-obj"; mkdir -p "$objdir"
    local f
    for f in "$src"/src/*.cc; do
        [[ "$(basename "$f")" == benchmark_main.cc ]] && continue
        g++ -std=c++17 -O3 -DNDEBUG -DHAVE_POSIX_REGEX -DBENCHMARK_STATIC_DEFINE \
            -I"$src/include" -I"$src/src" -c "$f" -o "$objdir/$(basename "${f%.cc}").o" &
    done
    wait
    ar rcs "$out" "$objdir"/*.o
    echo "$out"
}

build_shaders() {
    local spvdir="$ROOT/build/gpu/shaders"
    local slangc="${CHEATAH_SPACE_SLANGC:-$(command -v slangc || true)}"
    [[ -z "$slangc" ]] && for c in "$HOME"/Tools/vulkan-sdk/*/x86_64/bin/slangc; do
        [[ -x "$c" ]] && slangc="$c"
    done
    if [[ -z "$slangc" || ! -x "$slangc" ]]; then
        note "no slangc — device lanes will report themselves skipped"
        return 0
    fi
    mkdir -p "$spvdir"
    # entry point : -DKERNEL_ define. Kept in step with irbem.slang's own guards.
    local pair
    for pair in irbem_dipole_f32:DIPOLE irbem_igrf_f32:IGRF irbem_trace_i_f32:TRACE; do
        local entry="${pair%%:*}" def="${pair##*:}"
        local spv="$spvdir/$entry.spv"
        if [[ -f "$spv" && "$spv" -nt "$ROOT/space/irbem/gpu/irbem.slang" ]]; then continue; fi
        note "slangc $entry.spv"
        "$slangc" "$ROOT/space/irbem/gpu/irbem.slang" "-DKERNEL_$def" \
            -target spirv -entry "$entry" -stage compute -o "$spv"
    done
}

cmd_build() {
    mkdir -p "$WORK" "$JSON_DIR"
    build_shaders
    local bsrc blib
    bsrc="$(benchmark_src)"
    blib="$(build_benchmark_lib "$bsrc")"

    local inc=(-I"$ROOT"
               -I"$CHEATAH_DIR/stdlib/ndarray" -I"$CHEATAH_DIR/stdlib/builtins"
               -I"$CHEATAH_DIR/stdlib/fixarray" -I"$bsrc/include")
    local libs=("$blib" -ldl -lpthread)
    local defs=(-DCHEATAH_SPACE_IRBEM_SPV_DIR="\"$ROOT/build/gpu/shaders\"")

    # The device lane is opt-in BY INCLUDE PATH, exactly as the headers' __has_include expects.
    # Without cheatah-gpu-linalg the same source compiles host-only and the GPU rows say so.
    local gl="$CHEATAH_GPU_LINALG_DIR"
    if [[ -f "$gl/include/cheatah_gpu_linalg/context.hpp" && -f "$gl/build/libcheatah_gpu_linalg_volk.a" ]]; then
        local vk_inc=""
        for c in "$HOME"/Tools/vulkan-sdk/*/x86_64/include; do [[ -d "$c" ]] && vk_inc="$c"; done
        inc+=(-I"$gl/include" -I"$gl/../cheatah-gpu" -I"$gl/../cheatah-gpu/build/vk/_deps/volk-src")
        [[ -n "$vk_inc" ]] && inc+=(-I"$vk_inc")
        defs+=(-DCHEATAH_GPU_LINALG_SPV_DIR="\"$ROOT/build/gpu/shaders\"")
        libs=("$blib" "$gl/build/libcheatah_gpu_linalg_volk.a" -ldl -lpthread)
    else
        note "no cheatah-gpu-linalg build at $gl — host lanes only"
    fi

    note "compiling bench/irbem_bench.cpp"
    g++ "${BENCH_CXXFLAGS[@]}" "$ROOT/bench/irbem_bench.cpp" "${inc[@]}" "${defs[@]}" \
        "${libs[@]}" -o "$BIN"
    note "built $BIN"
}

need_bin() { [[ -x "$BIN" ]] || cmd_build; }

# =============================================================================================
# run
# =============================================================================================

cmd_run() {
    local suite="$1"
    local filter; filter="$(suite_filter "$suite")"
    need_bin
    mkdir -p "$JSON_DIR"
    note "running $suite (reps=$BENCH_REPS, min_time=$BENCH_MIN_TIME)"
    "$BIN" \
        --benchmark_filter="$filter" \
        --benchmark_repetitions="$BENCH_REPS" \
        --benchmark_min_time="$BENCH_MIN_TIME" \
        --benchmark_enable_random_interleaving=true \
        --benchmark_out_format=json \
        --benchmark_out="$JSON_DIR/$suite.json" \
        --benchmark_format=console
}

# Emit `run_name<TAB>aggregate<TAB>real_time<TAB>unit<TAB>items_per_second` for every aggregate row.
# A skipped benchmark (no device, no oracle) carries error_occurred and is emitted with a zero
# throughput, so a missing lane shows as a dash in the table rather than as a silent omission.
json_rows() {
    awk '
    /"run_name":/       { split($0, a, "\""); rn = a[4] }
    /"run_type":/       { split($0, a, "\""); rt = a[4] }
    /"aggregate_name":/ { split($0, a, "\""); an = a[4] }
    /"time_unit":/      { split($0, a, "\""); tu = a[4] }
    /"error_occurred":/ { err = 1 }
    /"real_time":/        { v = $2; sub(/,$/, "", v); rtime = v + 0 }
    /"items_per_second":/ { v = $2; sub(/,$/, "", v); ips = v + 0 }
    /^    }/ {
        if (rt == "aggregate" && rn != "")
            printf "%s\t%s\t%.10g\t%s\t%.10g\n", rn, an, rtime, tu, (err ? 0 : ips)
        rn = ""; rt = ""; an = ""; tu = ""; rtime = 0; ips = 0; err = 0
    }' "$1"
}

# =============================================================================================
# the tables
# =============================================================================================

# ns per point, from items_per_second — the one metric that is comparable across a scalar routine,
# a 65536-point batch and a device dispatch, because every benchmark reports how many points it
# processed via SetItemsProcessed.
NS_AWK='function ns(ips) { return ips > 0 ? 1e9 / ips : 0 }
        function fmt(x) { return x >= 1e6 ? sprintf("%.0f", x/1000) "\xc2\xb5s" \
                                          : (x >= 1000 ? sprintf("%.2f", x/1000) " \xc2\xb5s" \
                                                       : sprintf("%.2f", x) " ns") }'

table_per_routine() {
    local json="$1" manifest="$2"
    { json_rows "$json"; echo "@@"; cat "$manifest"; } | awk -F'\t' "$NS_AWK"'
    /^@@$/ { phase = 1; next }
    phase == 0 {
        if ($2 == "median") { med[$1] = ns($5) }
        if ($2 == "stddev") { sd[$1] = $5 }
        next
    }
    /^#/ { next }
    {
        label = $1; cpu = $2; gpu = $3; irb = $4; status = $5; primary = $6
        if (primary != "yes") next
        c = (cpu in med) ? med[cpu] : 0
        g = (gpu in med) ? med[gpu] : 0
        i = (irb in med) ? med[irb] : 0
        printf "| %s | %s | %s | %s | %s | %s | %s |\n",
            label,
            (c > 0 ? fmt(c) : "—"),
            (g > 0 ? fmt(g) : "—"),
            (c > 0 && g > 0 ? sprintf("%.2f\xc3\x97", c / g) : "—"),
            (i > 0 ? fmt(i) : "—"),
            (c > 0 && i > 0 ? sprintf("%.2f\xc3\x97", i / c) : "—"),
            status
    }'
}

table_crossover() {
    local json="$1"
    json_rows "$1" | awk -F'\t' "$NS_AWK"'
    $2 != "median" { next }
    {
        name = $1
        if (name !~ /trace_batch/) next
        n = name; sub(/^.*trace_batch\//, "", n); sub(/\/.*$/, "", n)
        if (name ~ /^BM_cpu_/) host[n + 0] = ns($5)
        if (name ~ /^BM_gpu_/) dev[n + 0]  = ns($5)
        seen[n + 0] = 1
    }
    END {
        m = 0
        for (n in seen) sizes[++m] = n + 0
        for (a = 1; a < m; a++) for (b = a + 1; b <= m; b++)
            if (sizes[a] > sizes[b]) { t = sizes[a]; sizes[a] = sizes[b]; sizes[b] = t }
        for (a = 1; a <= m; a++) {
            n = sizes[a]
            h = host[n]; d = dev[n]
            if (h <= 0 || d <= 0) { verdict = "—"; ratio = "—" }
            else {
                r = h / d
                ratio = sprintf("%.2f\xc3\x97", r)
                verdict = (r >= 1.0) ? "device wins" : "**host wins**"
            }
            printf "| %d | %s | %s | %s | %s |\n", n,
                (h > 0 ? fmt(h) : "—"), (d > 0 ? fmt(d) : "—"), ratio, verdict
        }
    }'
}

table_intensity() {
    local json="$1" model="$2"
    { json_rows "$json"; echo "@@"; cat "$model"; } | awk -F'\t' "$NS_AWK"'
    /^@@$/ { phase = 1; next }
    phase == 0 { if ($2 == "median") med[$1] = ns($5); next }
    /^#/ { next }
    {
        kernel = $1; label = $2; bytes = $3 + 0; flops = $4 + 0; fpb = $5 + 0
        h = med[$6]; d = med[$7]
        if (h <= 0 || d <= 0) { ratio = "—"; verdict = "not measured here" }
        else {
            r = h / d
            ratio = sprintf("%.2f\xc3\x97", r)
            verdict = (r >= 1.0) ? "GPU **wins**" : "GPU **loses**"
        }
        printf "| `%s` | %s | %d | %d | **%.1f** | %s | %s | %s | %s |\n",
            kernel, label, bytes, flops, fpb,
            (h > 0 ? fmt(h) : "—"), (d > 0 ? fmt(d) : "—"), ratio, verdict
    }'
}

# =============================================================================================
# publish — replace one BENCH region, stamp and all
# =============================================================================================

stamp() {
    local suite="$1" layout="$2" produced_by="$3"
    local competitors="none — the IRBEM oracle is dlopen'ed, see below"
    local oracle_note="absent"
    [[ -f "$SPACE_IRBEM_ORACLE" ]] && oracle_note="$SPACE_IRBEM_ORACLE"
    local device="${CHEATAH_GPU_LINALG_VK_DEVICE:-platform default}"
    cat <<EOF
<!-- cheatah-bench-stamp v1
     suite:        $suite
     generated:    $(date -u +%Y-%m-%d)
     commit:       $(git -C "$ROOT" rev-parse --short HEAD 2>/dev/null || echo unknown)
     host:         $(uname -n), $(nproc) CPUs @ $(awk '/cpu MHz/ {printf "%.0f", $4; exit}' /proc/cpuinfo 2>/dev/null || echo "?") MHz
     cpu-scaling:  $(cat /sys/devices/system/cpu/cpu0/cpufreq/scaling_governor 2>/dev/null || echo unknown)
     build:        $(g++ --version | head -1), Google Benchmark @ $BENCHMARK_PIN
     device:       $device
     competitors:  IRBEM $oracle_note (the -O2 REBUILD, never the shipped no-\`-O\` library)
     harness:      reps=$BENCH_REPS, min_time=$BENCH_MIN_TIME, random-interleaving=on
     statistic:    median over repetitions; every rate is REAL time, never CPU time — a device
                   lane is blocked, not running, for nearly all of its wall clock
     publishable:  true
     layout:       $layout
     watch:        space/irbem/, bench/irbem_bench.cpp

     PRODUCED BY:
       $produced_by
-->
EOF
}

replace_region() {
    local suite="$1" body_file="$2"
    grep -q "<!-- BENCH:$suite begin -->" "$PAGE" || die "no BENCH:$suite region in $PAGE"
    awk -v suite="$suite" -v body="$body_file" '
        $0 ~ "<!-- BENCH:" suite " begin -->" { print; while ((getline l < body) > 0) print l; skip = 1; next }
        $0 ~ "<!-- BENCH:" suite " end -->"   { skip = 0 }
        !skip { print }
    ' "$PAGE" > "$PAGE.tmp"
    mv "$PAGE.tmp" "$PAGE"
}

cmd_publish() {
    local suite="$1"
    if [[ "$suite" == all ]]; then for s in "${ALL_SUITES[@]}"; do cmd_publish "$s"; done; return; fi
    cmd_run "$suite"
    local json="$JSON_DIR/$suite.json"
    local body; body="$(mktemp)"
    local cmdline="SPACE_IRBEM_ORACLE='$SPACE_IRBEM_ORACLE' CHEATAH_GPU_LINALG_VK_DEVICE='${CHEATAH_GPU_LINALG_VK_DEVICE:-}' \\
           scripts/bench_run.sh publish $suite
       which runs: build/bench/irbem_bench --benchmark_filter='$(suite_filter "$suite")' \\
           --benchmark_repetitions=$BENCH_REPS --benchmark_min_time=$BENCH_MIN_TIME \\
           --benchmark_enable_random_interleaving=true \\
           --benchmark_out_format=json --benchmark_out=docs/bench/$suite.json"

    case "$suite" in
    irbem-per-routine)
        local man; man="$(mktemp)"; "$BIN" --manifest > "$man"
        { stamp "$suite" routines "$cmdline"; echo
          echo '| routine | CPU ns/call | GPU ns/call | GPU speedup | IRBEM `-O2` ns/call | vs IRBEM | on GPU? |'
          echo '|---|--:|--:|--:|--:|--:|---|'
          table_per_routine "$json" "$man"
        } > "$body"
        rm -f "$man" ;;
    irbem-crossover)
        { stamp "$suite" crossover "$cmdline"; echo
          echo '| field lines | host per line | device per line | speedup | verdict |'
          echo '|--:|--:|--:|--:|---|'
          table_crossover "$json"
        } > "$body" ;;
    irbem-intensity)
        local mod; mod="$(mktemp)"; "$BIN" --intensity > "$mod"
        { stamp "$suite" intensity "$cmdline"; echo
          echo '| kernel | computes | bytes/point | flops/point | flops/byte | host per point | device per point | speedup | verdict |'
          echo '|---|---|--:|--:|--:|--:|--:|--:|---|'
          table_intensity "$json" "$mod"
          echo
          grep '^# mean RK4' "$mod" | sed 's/^# /> /'
        } > "$body"
        rm -f "$mod" ;;
    esac
    replace_region "$suite" "$body"
    rm -f "$body"
    note "published $suite into $PAGE"
}

# =============================================================================================
# baseline + check
# =============================================================================================

# suite,key,ns_per_point — one row per measured lane. Keys are benchmark run_names, which are
# stable by construction: `crossover_sizes` deliberately does NOT set ->MinTime, because Google
# Benchmark writes that into the name and a name that moves when the harness is retuned would make
# every baseline row unjoinable.
emit_baseline() {
    local suite="$1"
    json_rows "$JSON_DIR/$suite.json" | awk -F'\t' -v s="$suite" '
        $2 == "median" && $5 > 0 { printf "%s,%s,%.6g\n", s, $1, 1e9 / $5 }'
}

cmd_baseline() {
    local s
    { echo "# suite,benchmark,ns_per_point — regenerate with scripts/bench_run.sh baseline"
      echo "# host: $(uname -n), $(g++ --version | head -1), device ${CHEATAH_GPU_LINALG_VK_DEVICE:-default}"
      echo "# commit: $(git -C "$ROOT" rev-parse --short HEAD 2>/dev/null || echo unknown), $(date -u +%Y-%m-%d)"
      for s in "${ALL_SUITES[@]}"; do cmd_run "$s" >/dev/null; emit_baseline "$s"; done
    } | sort -u -t, -k1,2 > "$BASELINE"
    note "wrote $BASELINE ($(grep -vc '^#' "$BASELINE") rows)"
}

# The regression gate. Tolerant BY RATIO and by an ABSOLUTE floor, because a routine measured in
# single nanoseconds moves by 30% on load alone and gating that would make the gate a coin flip.
# A row that regresses is re-measured once before it is called a regression, for the same reason.
REGRESS_RATIO="${REGRESS_RATIO:-1.40}"
REGRESS_FLOOR_NS="${REGRESS_FLOOR_NS:-2.0}"

cmd_check() {
    [[ -f "$BASELINE" ]] || die "no baseline at $BASELINE — run 'scripts/bench_run.sh baseline'"
    # The flag lists in this script and in bench/CMakeLists.txt must agree, or the published table
    # and a CMake-built binary are measuring two different programs.
    local f
    for f in "${BENCH_CXXFLAGS[@]:1}"; do
        grep -q -- "$f" "$ROOT/bench/CMakeLists.txt" ||
            die "bench/CMakeLists.txt does not carry '$f' — the two build paths have drifted"
    done
    if grep -q -- '-ffast-math' "$ROOT/bench/CMakeLists.txt"; then
        die "-ffast-math appears in bench/CMakeLists.txt; it is banned (see the flag comment)"
    fi

    local s now; now="$(mktemp)"
    for s in "${ALL_SUITES[@]}"; do cmd_run "$s" >/dev/null; emit_baseline "$s" >> "$now"; done

    local failures; failures="$(mktemp)"
    awk -F, -v ratio="$REGRESS_RATIO" -v floor="$REGRESS_FLOOR_NS" '
        NR == FNR { if ($0 !~ /^#/) base[$1 "," $2] = $3 + 0; next }
        {
            k = $1 "," $2; v = $3 + 0
            if (!(k in base)) next
            b = base[k]
            if (v > b * ratio && v - b > floor)
                printf "%s\t%.4g\t%.4g\t%.2f\n", k, b, v, v / b
        }' "$BASELINE" "$now" > "$failures"

    if [[ -s "$failures" ]]; then
        note "re-measuring $(wc -l < "$failures") suspected regressions before calling them one"
        local recheck; recheck="$(mktemp)"
        for s in "${ALL_SUITES[@]}"; do cmd_run "$s" >/dev/null; emit_baseline "$s" >> "$recheck"; done
        awk -F, -v ratio="$REGRESS_RATIO" -v floor="$REGRESS_FLOOR_NS" '
            NR == FNR { if ($0 !~ /^#/) base[$1 "," $2] = $3 + 0; next }
            {
                k = $1 "," $2; v = $3 + 0
                if (!(k in base)) next
                b = base[k]
                if (v > b * ratio && v - b > floor)
                    printf "%s\t%.4g\t%.4g\t%.2f\n", k, b, v, v / b
            }' "$BASELINE" "$recheck" > "$failures"
        rm -f "$recheck"
    fi

    if [[ -s "$failures" ]]; then
        printf 'bench_run: REGRESSION (>%sx and >%s ns slower than bench/baseline.csv)\n' \
            "$REGRESS_RATIO" "$REGRESS_FLOOR_NS" >&2
        printf '  %-64s %10s %10s %7s\n' "suite,benchmark" "baseline" "now" "ratio" >&2
        while IFS=$'\t' read -r k b v r; do
            printf '  %-64s %10s %10s %6sx\n' "$k" "$b" "$v" "$r" >&2
        done < "$failures"
        rm -f "$now" "$failures"; exit 1
    fi
    note "no regression against $BASELINE"
    rm -f "$now" "$failures"
}

# =============================================================================================

case "${1:-}" in
    build)     cmd_build ;;
    run)       cmd_run "${2:?usage: bench_run.sh run <suite>}" ;;
    publish)   cmd_publish "${2:-all}" ;;
    baseline)  cmd_baseline ;;
    check)     cmd_check ;;
    manifest)  need_bin; "$BIN" --manifest ;;
    intensity) need_bin; "$BIN" --intensity ;;
    *) sed -n '2,30p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; exit 1 ;;
esac
