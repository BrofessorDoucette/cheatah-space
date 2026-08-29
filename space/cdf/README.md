# space.cdf

NASA **Common Data Format** I/O for cheatah — the format space-physics data ships in
(variables, attributes, records, 0–10-D arrays, optional compression). Written **from scratch
in C++ with zero dependencies**: NASA's CDF library is *not* linked and *not* required to
build, install, or use `space.cdf`. Goal: the fastest CDF reader/writer we can build.

> **Status:** working — reads CDF 3.x, including GZIP-compressed variables, verified against real mission files. Writing is not implemented.

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

## What works today

`space.cdf` opens a CDF 3.x file, lists its variables, and hands their records back as cheatah
`ndarray`s — including GZIP-compressed variables, which is nearly everything in NASA's archive.

```purr
import io
import ndarray
import space.cdf as cdf

let f = cdf.open("omni_hro2_1min_20150101_v01.cdf")
io.print(cdf.var_names(f))                  # 47 variables
io.print(cdf.data_type_name(f, "F"))        # CDF_REAL4
io.print(cdf.record_count(f, "F"))          # 44640 — a month of one-minute samples

let b = cdf.values(f, "F")                  # ndarray[float64], IMF magnitude in nT
io.print(b[0], b[1], b[2])                  # 6.92 5.84 5.71
```

| function | what it gives you |
| --- | --- |
| `open(path)` | maps the file, parses every descriptor |
| `var_names(f)` | every variable, rVariables first |
| `record_count(f, name)` | how many records the variable has |
| `data_type_name(f, name)` | `CDF_REAL4`, `CDF_EPOCH`, `CDF_TIME_TT2000`, … |
| `shape(f, name)` | records, then each varying dimension |
| `values(f, name)` | `ndarray[float64]` — what you plot |
| `values_i64(f, name)` | `ndarray[int]` — exact for `CDF_INT8` and `CDF_TIME_TT2000` |

`values` is lossless for every type a double can hold exactly, which is all of them except
`CDF_INT8` and `CDF_TIME_TT2000`; asking for those is an error naming `values_i64` rather than a
silent rounding. `values_i64` refuses the float types for the same reason.

### Plotting

The arrays go straight into [cheatah-plot](https://github.com/BrofessorDoucette/cheatah-plot).
`examples/purr_space/03_plot_omni_imf.purr` is the worked version — open a CDF, decode `Epoch`
and `F`, drop the fill values, and render a PNG:

```sh
scripts/cdf-corpus.sh fetch --tier 1
examples/purr_space/run_examples.sh        # writes examples/purr_space/out/*.png
```

Two things that example exists to show. cheatah-plot has no time axis — its ticks format as
`%.6g` — so a raw `CDF_EPOCH` of 6.3e13 ms renders as one repeated label; rebase onto hours and
name the epoch in the axis label. And OMNI marks a missing sample as ~1e31, which would flatten
every real variation into a line at the bottom of the frame, so the fill values are dropped
rather than plotted.

## What is deferred

The reader refuses what it does not yet handle, by name, rather than guessing:

| refused | error |
| --- | --- |
| CDF 2.x and earlier | `UnsupportedPreV26` |
| whole-file (CCR) compression | `UnsupportedCompression` |
| multi-file CDFs (`.v0`/`.z0` companions) | `UnsupportedMultiFile` |
| VAX F/D/G float encodings | `UnsupportedEncoding` |
| sparse records with real gaps | `UnsupportedSparse` |
| RLE / Huffman / adaptive Huffman | `UnsupportedCompression` |
| N-D variables in column-major files | `UnsupportedLayout` |

Also still to come: the attributes API, the writer, the MD5 checksum, cryptographic signing over
a canonical content digest, and the byte-for-byte differential harness against NASA's own
library. RLE, Huffman and adaptive Huffman are last on purpose — the format specification names
those algorithms but never documents their bitstreams, and no file in the public archive uses
any of them, so they have to be reverse-engineered from files generated on demand.

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
