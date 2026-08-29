# Getting started with cheatah-space

@brief Install the package, then read a real mission file, compute a magnetic coordinate, and plot
the result — the three modules and cheatah-plot working as one.

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
magnetic coordinates the belts are actually organised by — L\*, McIlwain's L, the mirror field —
against IGRF-14 and a choice of external field models.

> **Status:** this module is **C++-only today**. Its purr surface is not built yet, so the block
> below is C++ while the rest of this page is purr. Everything it shows compiles and runs.

```cpp
#include "space/irbem/irbem.hpp"
namespace ib = cheatah::space::irbem;

const auto model = ib::Igrf<13>::at(2015.5).value();
const auto rot   = ib::api::rotations_at(2015, 182, 43200.0, model);

const ib::Position<ib::Frame::GEO> p{cheatah::fixarray::vec3d{6.6, 0.0, 0.0}};
const auto shell = ib::make_lstar(model, rot.value, p, 90.0,
                                  ib::DriftShellOptions::from_irbem(0, 0));

if (shell.status == ib::Status::Ok) {
    shell.value.lstar;   // Roederer's L*, Earth radii
    shell.value.b_min;   // |B| at the magnetic equator of that field line, nT
}
```

Every result carries a named @ref cheatah::space::irbem::Status rather than a sentinel: an open
field line, a point outside a model's fitted envelope, and a shell that would not close are each
something you branch on, never a `-1e31` you have to recognise.

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
