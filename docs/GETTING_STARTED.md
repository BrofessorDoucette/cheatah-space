# Getting started with cheatah-space

Install the package, then read a real mission file, compute a magnetic coordinate, and plot the
result — the three modules and cheatah-plot working as one.

`space.time`, `space.cdf` and `space.irbem` are useful on their own, but the reason they ship
together is the workflow below: a file of spacecraft data goes in, and a figure of something
physical comes out. Nothing here is pseudocode — every block is taken from a program in this
repository that runs.

```sh
biome add cheatah-space
biome add cheatah-plot     # only for the plotting step
```

## 1. Open a real file and get arrays back

`space.cdf` reads NASA's Common Data Format — the format essentially all space-physics data ships
in — from scratch, with NASA's own library neither linked nor required. Variables come back as
cheatah `ndarray`s, so they go straight into the rest of the numeric stack.

```purr
import io
import ndarray
import space.cdf as cdf

let f = cdf.open("omni_hro2_1min_20150101_v01.cdf")
io.print(cdf.var_names(f))                  # every variable in the file
io.print(cdf.record_count(f, "F"))          # 44640 — a month of one-minute samples

let b = cdf.values(f, "F")                  # ndarray[float64], IMF magnitude in nT
```

GZIP-compressed variables — nearly everything in the archive — decode transparently. `values` is
lossless for every type a double holds exactly; for `CDF_INT8` and `CDF_TIME_TT2000` it refuses and
names `values_i64` rather than rounding quietly.

## 2. Put the timestamps on a scale you can reason about

`space.time` is the bridge between CDF's epochs and the Julian dates the rest of astronomy uses.

```purr
import space.time as st
import ndarray

let jd  = st.unix_to_jd(1420070400.0)                   # one timestamp -> Julian Date
let jds = st.unix_to_jd(ndarray.array([0.0, 86400.0]))  # or a whole array, vectorized
```

## 3. Compute a magnetic coordinate

`space.irbem` is a from-scratch reimplementation of the radiation-belt library
[PRBEM/IRBEM](https://github.com/PRBEM/IRBEM), written to the published papers. It gives you the
coordinates the belts are actually organised by — L\*, McIlwain's L, the mirror field — against
IGRF-14.

Build the epoch once and reuse it; it carries the internal field model and the frame rotations for
that instant.

```purr
import io
import fixarray
import space.irbem.purr as irbem

let e = irbem.epoch_at(2015, 182, 43200.0)      # 2015 day 182, 12:00 UT

# 6.6 Rₑ on the GEO x-axis — geosynchronous — for a 90° (equatorially mirroring) particle.
# sysaxes is IRBEM's frame code: 0 GDZ, 1 GEO, 2 GSM, 3 GSE, 4 SM, 5 GEI, 6 MAG, 7 SPH, 8 RLL.
let c = irbem.make_lstar(e, 6.6, 0.0, 0.0, 1, 90.0)

io.print("Lm      =", c[0])        # 6.66772
io.print("L*      =", c[1])        # 6.66194
io.print("Blocal  =", c[2], "nT")
io.print("status  =", irbem.status_name(int(c[6])))
```

Frames convert by the same integer codes, and the transform is exact both ways:

```purr
let gsm = irbem.coord_trans(e, 6.6, 0.0, 0.0, 1, 2)     # GEO -> GSM
let mlt = irbem.get_mlt(e, 6.6, 0.0, 0.0, 1)            # magnetic local time, hours
```

Nothing here returns a sentinel. The status slot names why a call declined — an open field line, a
point outside a model's fitted envelope, a shell that would not close — so you branch on a reason
rather than recognising `-1e31`.

## 4. Plot it

The arrays go straight into [cheatah-plot](https://github.com/BrofessorDoucette/cheatah-plot).
This is `examples/purr_space/03_plot_omni_imf.purr` in full — a month of interplanetary magnetic
field, read from the CDF and rendered to a PNG:

```purr
import plot
import plot.figure as figure

let fig = figure.line(figure.new_figure(), x, y)
fig = figure.title(fig, "OMNI interplanetary magnetic field, January 2015")
fig = figure.xlabel(fig, "hours since 2015-01-01T00:00:00Z")
fig = figure.ylabel(fig, "|B| (nT)")
fig = figure.size(fig, 1280, 720)
plot.save(fig, "out/03_plot_omni_imf.png")
```

Two things in that program are worth copying. Time is rebased to hours since the first sample,
because the plotter has no time axis and a raw `CDF_EPOCH` of 6.3587e13 ms would render as one
indistinguishable label repeated across the axis. And the fill values are dropped rather than
plotted: OMNI marks a missing sample as ~1e31, which would flatten every real variation into a line
at the bottom of the frame.

## Where to go next

| you want | read |
| --- | --- |
| what each module does | the module pages in the sidebar |
| which models are verified, and over what part of the domain | `space.irbem`'s verification page |
| how fast it is, and against what | `space.irbem`'s benchmarks page |
| the runnable programs | `examples/purr_space/` in the repository |

The OMNI file the examples use is fetched on demand — `scripts/cdf-corpus.sh fetch --tier 1` — and
never committed; each example prints what to run if it is absent rather than failing obscurely.
