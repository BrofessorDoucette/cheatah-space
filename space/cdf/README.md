# space.cdf

NASA **Common Data Format** I/O for cheatah — the format space-physics data ships in
(variables, attributes, records, 0–10-D arrays, optional compression). Written **from scratch
in C++ with zero dependencies**: NASA's CDF library is *not* linked and *not* required to
build, install, or use `space.cdf`. Goal: the fastest CDF reader/writer we can build.

> Reference: <https://cdf.gsfc.nasa.gov/>, and the authoritative
> [CDF 3.9 Internal Format Description](https://spdf.gsfc.nasa.gov/pub/software/cdf/doc/cdf39/cdf39ifd.pdf).
> We implement from that specification; we do not read or copy NASA's source.
>
> NASA's CDF distribution is **not public domain** — per
> [CDF_copyright.txt](https://spdf.gsfc.nasa.gov/pub/software/cdf/dist/CDF_copyright.txt) it "may
> be copied or redistributed as long as it is not sold for profit", may be incorporated into other
> products, and carries modification-notice requirements. So we use it **only as an optional,
> dev-only oracle** — to verify our output byte-for-byte and to benchmark against — fetched on
> demand by `scripts/cdf-oracle.sh` into `space/cdf/vendor/` (git-ignored), never committed, never
> shipped, and gated behind a flag. `scripts/check_no_vendored_nasa.sh` enforces that on every
> push. A normal `biome add cheatah-space` never touches any of it.

## Plan

**Phase 0 — bridge (done in `space.time`).** `unix_to_cdf_epoch` / `cdf_epoch_to_unix` convert
the CDF_EPOCH (ms since year 0). TT2000 (ns since J2000 + leap seconds) and EPOCH16
(picoseconds) land here with their leap-second table.

**Phase 1 — native reader.** Parse the single-file CDF layout directly: the descriptor records
(CDR/GDR), variable descriptors (zVDR/rVDR), variable index/value records (VXR/VVR), and
attribute records (ADR/AEDR). Decompress (RLE, Huffman, adaptive Huffman, GZIP — own
implementations, no zlib dependency) only the records actually touched. Decode straight into
cheatah `ndarray`; memory-map + zero-copy slice where the layout allows.

**Phase 2 — native writer + MD5 checksum + compression.**

**Verification & benchmark (optional).** A dev harness diffs our reader against the NASA C lib
over a fixed corpus (scalar-record, large multi-dim, compressed) and records ns/record + GB/s.
Each function's measured speedup vs the reference becomes its `@perf` tag; we ship only when we
win. The harness is skipped automatically when the reference isn't present.

## API sketch (subject to change)

```purr
import space.cdf
let f = cdf.open("mms1_fgm.cdf")          # parse descriptors lazily
let t = f.var("Epoch").records(0, 1000)    # ndarray, decoded once
```

## House rules

Every function documents `@complexity`, `@alloc`, truthful `@systest` tags, and — because a
reference exists — `@perf` stating the speedup over NASA's C lib.
