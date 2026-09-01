// Copyright (c) 2026 BigBrain LLC. MIT-licensed (see LICENSE).

/**
 * @file irbem_bench.cpp
 * @brief space.irbem — the Google Benchmark suite, and the manifest that ties each routine to its
 *        GPU kernel and its IRBEM counterpart.
 *
 * Every performance claim this module makes has to come from here. Before this file existed the
 * headers carried prose ("306 ns/eval", "8.37x the host") that nothing regenerated and nothing
 * checked, and a prose number is a number that silently rots the first time somebody reorders a
 * loop. `scripts/bench_run.sh` runs this binary and GENERATES the tables of
 * `space/irbem/BENCHMARKS.md` from its JSON, so a number in the documentation is by construction a
 * number this binary measured on the machine that printed it.
 *
 * ### Three lanes, one row
 *
 * A row of the generated table is one routine measured on up to three lanes:
 *
 *  - **CPU** — this module's own header, fp64, `-O3 -march=native -ffp-contract=off`.
 *  - **GPU** — the same physics through `gpu/dispatch.hpp`, fp32, transfers INCLUDED. A routine has
 *    a GPU number only when `dispatch.hpp`'s `registered_kernels` actually carries a launchable row
 *    for it; see @ref gpu_status, which derives that rather than asserting it.
 *  - **IRBEM** — the vendored Fortran library, `dlopen`ed as a black box (the clean-room rule: run
 *    it, never read it) through the entry points documented in `matlab/libirbem.h` and
 *    `docs/source/api/` .rst files.
 *
 * The IRBEM lane quotes the `-O2` REBUILD, never the shipped binary. IRBEM's own makefile passes no
 * `-O` at all and `docs/ERROR_BUDGET.md` §5 measures the shipped library at 2.7x the `-O2` one:
 * benchmarking against the shipped artefact would inflate every speedup in the table by that factor
 * and the resulting claims would be about GNU Fortran's default optimization level, not about this
 * code. The default oracle is therefore `/tmp/irbem-builds/libirbem-O2.so`, overridable with
 * `SPACE_IRBEM_ORACLE`.
 *
 * ### Hot and cold, because the API has a seam IRBEM's does not
 *
 * `Rotations::at` spends ~14 transcendentals ONCE per epoch and every point of the ephemeris then
 * costs a 3x3 product; `helio_geometry` is the same shape. IRBEM's `geo2gsm1_` cannot amortize —
 * it takes a date with every point. Reporting only the amortized number would flatter this module
 * and reporting only the per-point rebuild would hide its actual design advantage, so every
 * transform that has that seam is benchmarked BOTH ways and appears as two rows: `(hot)` reuses one
 * `Rotations`, `(cold)` rebuilds it per point and is the like-for-like against the oracle.
 *
 * ### Elision-proofing
 *
 * These routines are small, pure and `constexpr`; a naive loop over literal inputs measures the
 * optimizer, not the code. Every benchmark therefore reads its inputs through a pointer laundered
 * by `benchmark::DoNotOptimize` — the compiler cannot fold what it cannot see — and sinks every
 * result through `DoNotOptimize` as well, so nothing may be hoisted out of the timed loop. The
 * laundering is done once per outer iteration and amortized over @ref kPoints inputs, so its own
 * cost is a fraction of a nanosecond per point rather than a per-point store.
 *
 * ### Inputs
 *
 * Never zeros, never the origin, never a pole, never an axis. A benchmark at a special case
 * measures the special case: an origin position takes `dipole_field_at`'s `r2 <= 0` early return,
 * a pole takes the shortest path through the geodetic iteration, and a date at a year boundary
 * takes the carry branch in `date_and_time_from_decimal_year`. @ref inputs builds a small
 * pseudo-ephemeris sweeping L = 1.5 to 8 Re over a full turn of longitude and +-58 deg of latitude,
 * and dates spread across the whole IGRF-valid span. It is deterministic — the same points on every
 * run, so a regression is a regression and not a different sample.
 *
 * @note No allocation happens inside any HOST timed loop. The input set is built once at start-up
 *       and the batch lanes write into buffers sized once, so a longer run must not allocate more.
 *       Verified with `valgrind --tool=memcheck`, reading `total heap usage`, at two iteration
 *       counts a hundred times apart: `BM_cpu_julian_day_number` 1213/1213 and
 *       `BM_cpu_geo_to_gdz` 1156/1156 and `BM_cpu_igrf_evaluate_deg13` 1251/1251 at 200x and
 *       20000x; `BM_cpu_trace_invariant` 1202/1202 and `BM_cpu_igrf_batch` 1162/1162 at 3x and
 *       60x. Identical in every case, so nothing in a timed loop reaches the heap.
 *
 *       The DEVICE lanes are excluded from that claim and do allocate per launch — staging vectors
 *       and pooled device buffers — which is a property of the seam, not of this file, and is
 *       inside the timed region because a caller pays it.
 */

#include <benchmark/benchmark.h>
#include <dlfcn.h>

#include <array>
#include <cctype>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <optional>
#include <span>
#include <sstream>
#include <string>
#include <string_view>
#include <vector>

#include "space/irbem/coords_geodetic.hpp"
#include "space/irbem/coords_helio.hpp"
#include "space/irbem/coords_rotations.hpp"
#include "space/irbem/ext_t89.hpp"
#include "space/time/calendar.hpp"
#include "space/irbem/gpu/dispatch.hpp"
#include "space/irbem/batch_soa.hpp"
#include "space/irbem/igrf.hpp"
#include "space/irbem/lstar.hpp"

namespace {

using namespace cheatah::space::irbem;  // NOLINT(google-build-using-namespace) — a bench TU
// Named rather than a whole-namespace using: `space::time` and `space::irbem` both have a
// `detail`, and pulling both in wholesale makes every `detail::` in this file ambiguous.
using cheatah::space::time::CalendarDate;
using cheatah::space::time::DateTime;
using cheatah::space::time::calendar_date;
using cheatah::space::time::date_and_time_from_decimal_year;
using cheatah::space::time::date_and_time_from_doy_and_ut;
using cheatah::space::time::day_of_year;
using cheatah::space::time::decimal_year;
using cheatah::space::time::julian_day_number;
namespace fx = cheatah::fixarray;

// ---------------------------------------------------------------------------------------------
// Inputs
// ---------------------------------------------------------------------------------------------

/// How many distinct inputs one timed iteration walks.
///
/// Large enough that the `DoNotOptimize` laundering and the loop counter are a small fraction of
/// a ~2 ns routine, small enough that every input stays in L1 (64 points x 24 bytes = 1.5 KB) so
/// the number measures arithmetic and not the memory system.
constexpr std::size_t kPoints = 64;

/// How many points the batch lanes (host-fp32 and device) move per iteration.
///
/// 2^18 is where `dispatch.hpp`'s registry records the device lane peaking on this hardware; below
/// it the launch latency dominates and above it nothing changes, so it is the size at which the
/// device number is a statement about the kernel rather than about the driver.
constexpr std::size_t kBatch = 1U << 18U;

/// A calendar instant, in the broken-down form the datetime routines take.
struct DateInput {
    int year;         ///< Astronomical year numbering.
    int month;        ///< 1..12.
    int day;          ///< 1-based day of month.
    int hour;         ///< UT hour.
    int minute;       ///< UT minute.
    int second;       ///< UT second.
    int doy;          ///< Day of year, precomputed so the day-of-year lanes do not re-derive it.
    double ut;        ///< UT seconds of day.
    double decy;      ///< The same instant as a decimal year.
    std::int64_t jdn; ///< The same instant's Julian Day Number.
    double jd_ut1;    ///< The same instant as a UT1 Julian Date (JDN - 0.5 + ut/86400).
};

/// The whole input set, built once.
struct Inputs {
    std::array<DateInput, kPoints> dates{};              ///< Instants spanning 1905..2029.
    std::array<Position<Frame::GEO>, kPoints> geo{};     ///< GEO Cartesian, Re.
    std::array<Position<Frame::GSM>, kPoints> gsm{};     ///< GSM Cartesian, Re.
    std::array<Position<Frame::SPH>, kPoints> sph{};     ///< Geographic spherical.
    std::array<Position<Frame::GDZ>, kPoints> gdz{};     ///< Geodetic (alt km, lat, lon).
    std::array<Position<Frame::RLL>, kPoints> rll{};     ///< Radius/lat/lon.
    std::array<Position<Frame::GSE>, kPoints> gse{};     ///< GSE Cartesian, Re.
    std::array<Position<Frame::HAE>, kPoints> hae{};     ///< HAE, AU.
    std::array<Position<Frame::HEE>, kPoints> hee{};     ///< HEE, AU.
    std::array<Position<Frame::HEEQ>, kPoints> heeq{};   ///< HEEQ, AU.
    std::vector<float> batch_pos;                        ///< 3*kBatch MAG-frame floats, Re.
    std::vector<float> batch_out;                        ///< 3*kBatch floats, written by the lanes.
    DipoleCoefficients dipole{};                         ///< Degree-1 IGRF at the reference epoch.
    Rotations rotations;                                 ///< Prebuilt rotations for the hot lanes.
    HelioGeometry helio{};                               ///< Prebuilt heliospheric geometry.
    /// IGRF-14 at full degree, for the reference epoch. `std::optional` because `Igrf`'s default
    /// constructor is private — the only way to obtain one is `Igrf::at`, which is the point.
    std::optional<Igrf<13>> igrf13{};
    /// The same model truncated to degree 10, so the table shows what the truncation buys.
    std::optional<Igrf<10>> igrf10{};
};

/// The reference epoch every "hot" lane's prebuilt geometry belongs to: 2015-180 12:00 UT, the
/// epoch `tools/oracle/convergence.cpp` measures the error budget at, so the two studies describe
/// the same instant.
constexpr int kRefYear = 2015;
constexpr int kRefDoy = 180;      ///< Day of year of the reference epoch.
constexpr double kRefUt = 43200.0; ///< UT seconds of the reference epoch (12:00).

/// Build the input set.
///
/// The sweep is deliberate rather than random: `k/kPoints` walks radius from 1.5 to 8 Re (the inner
/// belt out past geosynchronous), longitude through a full turn, and latitude through +-58 deg, so
/// no lane sees a repeated value, an axis, a pole or the origin. Dates step 2 years from 1903 so
/// the IGRF interpolation bracket moves, leap years and century non-leap years are both crossed,
/// and every instant stays inside the model's `[1900, 2030]` validity — which the 1905 start this
/// began with did NOT: it put the last of the 64 dates at 2031 and the oracle said so, once per
/// run, in a warning nobody was reading.
///
/// @return the populated set; it is a function-local static, built on first use.
/// @complexity O(kPoints + kBatch).
/// @alloc two `std::vector<float>` of `3*kBatch` floats, once, at start-up — never inside a timed
///        loop.
const Inputs& inputs() {
    static const Inputs* built = [] {
        // `Rotations` has no default constructor — every rotation set belongs to an epoch — so it
        // has to be initialised here; the two batch vectors are named alongside it only to keep
        // -Wmissing-field-initializers quiet, and are resized properly further down.
        auto* in = new Inputs{
            .batch_pos = {}, .batch_out = {},
            .rotations = Rotations::at(2451545.0, DipoleCoefficients{-29404.8, -1450.9, 4652.5})};
        for (std::size_t k = 0; k < kPoints; ++k) {
            const double t = static_cast<double>(k) / static_cast<double>(kPoints);
            const double r = 1.5 + (6.5 * t);
            const double lat_deg = 58.0 * std::sin(6.0 * t);
            const double lon_deg = -180.0 + (359.0 * t);
            const double lat = lat_deg * detail::kDegToRad;
            const double lon = lon_deg * detail::kDegToRad;
            const double cl = std::cos(lat);
            in->geo[k] = Position<Frame::GEO>{
                fx::vec3d{r * cl * std::cos(lon), r * cl * std::sin(lon), r * std::sin(lat)}};
            in->gsm[k] = Position<Frame::GSM>{in->geo[k].v};
            in->gse[k] = Position<Frame::GSE>{in->geo[k].v};
            in->sph[k] = Position<Frame::SPH>{fx::vec3d{r, lat_deg, lon_deg}};
            in->rll[k] = Position<Frame::RLL>{fx::vec3d{r, lat_deg, lon_deg}};
            in->gdz[k] = Position<Frame::GDZ>{fx::vec3d{200.0 + (35000.0 * t), lat_deg, lon_deg}};
            // Heliospheric inputs are in AU and sit near the Earth's own orbit, which is the only
            // place HEE/HAE/HEEQ are ever evaluated in practice.
            const double au = 0.983 + (0.034 * t);
            in->hae[k] = Position<Frame::HAE>{
                fx::vec3d{au * std::cos(lon), au * std::sin(lon), 0.004 * std::sin(lat)}};
            in->hee[k] = Position<Frame::HEE>{in->hae[k].v};
            in->heeq[k] = Position<Frame::HEEQ>{in->hae[k].v};

            // 1903 + 2k, so the last of kPoints = 64 dates is 2029 and every instant is
            // inside IGRF's [1900, 2030] span. It used to start at 1905, which put the
            // last date at 2031 and made the oracle print
            //   *** WARNING -- Input year = 2031.50 is out of valid range 1900-2030
            // on every run: one point in sixty-four was timing IRBEM's clamp path.
            const int year = 1903 + static_cast<int>(2 * k);
            const int month = 1 + static_cast<int>(k % 12);
            const int day = 1 + static_cast<int>(k % 28);
            const int hour = static_cast<int>(k % 24);
            const int minute = static_cast<int>((7 * k) % 60);
            const int second = static_cast<int>((13 * k) % 60);
            in->dates[k] = DateInput{
                .year = year,
                .month = month,
                .day = day,
                .hour = hour,
                .minute = minute,
                .second = second,
                .doy = day_of_year(year, month, day),
                .ut = (3600.0 * hour) + (60.0 * minute) + second,
                .decy = decimal_year(year, month, day, hour, minute, second),
                .jdn = julian_day_number(year, month, day),
                .jd_ut1 = static_cast<double>(julian_day_number(year, month, day)) - 0.5 +
                          (((3600.0 * hour) + (60.0 * minute) + second) / 86400.0),
            };
        }

        in->igrf13 = Igrf<13>::at(2015.5);
        in->igrf10 = Igrf<10>::at(2015.5);
        in->dipole = DipoleCoefficients{static_cast<double>(in->igrf13->g(1, 0)),
                                        static_cast<double>(in->igrf13->g(1, 1)),
                                        static_cast<double>(in->igrf13->h(1, 1))};
        const double jd_ref = static_cast<double>(julian_day_number(kRefYear, 1, 1)) +
                              (kRefDoy - 1) - 0.5 + (kRefUt / 86400.0);
        in->rotations = Rotations::at(jd_ref, in->dipole);
        in->helio = helio_geometry(jd_ref - 2400000.5);

        in->batch_pos.resize(3 * kBatch);
        in->batch_out.resize(3 * kBatch);
        for (std::size_t i = 0; i < kBatch; ++i) {
            const Position<Frame::GEO>& p = in->geo[i % kPoints];
            // Perturb per point so the device cannot be reading one cache line for the whole
            // batch, and so no point repeats bit-for-bit across the batch.
            const auto jitter = static_cast<float>(1.0 + (1e-6 * static_cast<double>(i % 1024)));
            in->batch_pos[(3 * i) + 0] = static_cast<float>(p.v[0]) * jitter;
            in->batch_pos[(3 * i) + 1] = static_cast<float>(p.v[1]) * jitter;
            in->batch_pos[(3 * i) + 2] = static_cast<float>(p.v[2]) * jitter;
        }
        return in;
    }();
    return *built;
}

/// Launder a pointer so the optimizer cannot see what it points at.
/// @tparam T the pointee type.
/// @param p the pointer.
/// @return the same address, through an opaque asm barrier.
/// @complexity O(1) — no instruction is emitted, only a scheduling and aliasing barrier.
/// @alloc none.
template <class T>
[[nodiscard]] const T* opaque(const T* p) {
    const T* q = p;
    benchmark::DoNotOptimize(q);
    return q;
}

/// Record how many points a benchmark processed, so `bench_run.sh` can derive ns/point from
/// `items_per_second` uniformly across scalar, batched and device lanes.
/// @param state the benchmark state, after its loop has run.
/// @param per_iteration how many points one iteration processed.
/// @complexity O(1).
/// @alloc none.
void set_points(benchmark::State& state, std::size_t per_iteration) {
    state.SetItemsProcessed(state.iterations() * static_cast<std::int64_t>(per_iteration));
}

// ---------------------------------------------------------------------------------------------
// The IRBEM oracle — dlopen'ed, never linked, never read
// ---------------------------------------------------------------------------------------------

/// Every IRBEM entry point this suite times, resolved once.
///
/// The signatures are the documented C bindings from `matlab/libirbem.h` and the `:callseq FORTRAN:`
/// lines of `docs/source/api/*.rst`; each was confirmed against the shipped library by evaluating a
/// known value before being used to time anything (`julday(2000,1,1) == 2451545`,
/// `get_doy(2000,3,1) == 61`, `date_and_time2decy(2000,7,2,0,0,0) == 2000.5`). Fortran passes
/// everything by reference, so every parameter is a pointer.
struct Oracle {
    void* handle = nullptr;  ///< The `dlopen` handle, or null when the library was not found.

    int (*julday)(int*, int*, int*) = nullptr;                        ///< `JULDAY(year, month, day)`.
    void (*caldat)(int*, int*, int*, int*) = nullptr;                 ///< `CALDAT(jdn, y, m, d)`.
    int (*get_doy)(int*, int*, int*) = nullptr;                       ///< `GET_DOY(year, month, day)`.
    void (*date2decy)(int*, int*, int*, int*, int*, int*, double*) = nullptr;  ///< `DATE_AND_TIME2DECY`.
    void (*decy2date)(double*, int*, int*, int*, int*, int*, int*, int*, double*) = nullptr;  ///< `DECY2DATE_AND_TIME`.
    void (*doyut2date)(int*, int*, double*, int*, int*, int*, int*, int*) = nullptr;  ///< `DOY_AND_UT2DATE_AND_TIME`.

    void (*geo2gsm)(int*, int*, double*, double*, double*, double*) = nullptr;  ///< `geo2gsm1`.
    void (*gsm2geo)(int*, int*, double*, double*, double*, double*) = nullptr;  ///< `gsm2geo1`.
    void (*geo2gse)(int*, int*, double*, double*, double*) = nullptr;           ///< `geo2gse1`.
    void (*geo2sm)(int*, int*, double*, double*, double*) = nullptr;            ///< `geo2sm1`.
    void (*geo2gei)(int*, int*, double*, double*, double*) = nullptr;           ///< `geo2gei1`.
    void (*gsm2sm)(int*, int*, double*, double*, double*) = nullptr;            ///< `gsm2sm1`.
    void (*geo2mag)(int*, double*, double*) = nullptr;                          ///< `geo2mag1`.

    void (*gdz_geo)(double*, double*, double*, double*, double*, double*) = nullptr;  ///< `gdz_geo`.
    void (*geo_gdz)(double*, double*, double*, double*, double*, double*) = nullptr;  ///< `geo_gdz`.
    void (*sph_car)(double*, double*, double*, double*) = nullptr;                    ///< `SPH_CAR`.
    void (*car_sph)(double*, double*, double*, double*) = nullptr;                    ///< `CAR_SPH`.
    void (*rll_gdz)(double*, double*, double*, double*) = nullptr;                    ///< `RLL_GDZ`.

    void (*hae2hee)(int*, int*, double*, double*, double*) = nullptr;   ///< `hae2hee1`.
    void (*hee2hae)(int*, int*, double*, double*, double*) = nullptr;   ///< `hee2hae1`.
    void (*hae2heeq)(int*, int*, double*, double*, double*) = nullptr;  ///< `hae2heeq1`.
    void (*heeq2hae)(int*, int*, double*, double*, double*) = nullptr;  ///< `heeq2hae1`.
    void (*gse2hee)(int*, int*, double*, double*, double*) = nullptr;   ///< `gse2hee1`.
    void (*hee2gse)(int*, int*, double*, double*, double*) = nullptr;   ///< `hee2gse1`.

    void (*get_field)(int*, int*, int*, int*, int*, double*, double*, double*, double*, double*,
                      double*, double*) = nullptr;  ///< `get_field1`.
    /// `make_lstar1(ntime, kext, options, sysaxes, iyear, idoy, ut, x1, x2, x3, maginput,
    /// Lm, Lstar, Blocal, Bmin, XJ, MLT)` — the trace-and-invariants entry point, and with
    /// `options(1) = 1` the drift-shell one too.
    void (*make_lstar)(int*, int*, int*, int*, int*, int*, double*, double*, double*, double*,
                       double*, double*, double*, double*, double*, double*, double*) = nullptr;
    /// `make_lstar_shell_splitting1(ntime, Nipa, kext, options, sysaxes, iyear, idoy, ut,
    /// x1, x2, x3, alpha, maginput, Lm, Lstar, Bmirr, Bmin, XJ, MLT)` — the SAME calculation for an
    /// arbitrary local pitch angle, and therefore the one that is actually the counterpart of
    /// @ref BM_cpu_trace_invariant. `make_lstar` above is documented as computing "the L\* parameter
    /// for locally mirroring particles (local pitch angle of 90 degrees)"
    /// (`docs/source/api/magnetic_coordinates.rst`), which integrates a different arc.
    void (*make_lstar_ss)(int*, int*, int*, int*, int*, int*, int*, double*, double*, double*,
                          double*, double*, double*, double*, double*, double*, double*, double*,
                          double*) = nullptr;
};

/// Where the oracle is expected. The `-O2` REBUILD, never the shipped binary — see the file brief.
constexpr const char* kDefaultOraclePath = "/tmp/irbem-builds/libirbem-O2.so";

/// Resolve the oracle once.
/// @return the resolved entry points; `handle` is null and every pointer is null when the library
///         is absent, which the oracle benchmarks turn into a skip rather than a crash.
/// @complexity O(1) after the first call.
/// @alloc `dlopen`'s, once.
const Oracle& oracle() {
    static const Oracle* o = [] {
        auto* out = new Oracle{};
        const char* env = std::getenv("SPACE_IRBEM_ORACLE");
        out->handle = dlopen(env != nullptr ? env : kDefaultOraclePath, RTLD_NOW);
        if (out->handle == nullptr) return out;
        auto sym = [h = out->handle](const char* name) { return dlsym(h, name); };
        out->julday = reinterpret_cast<decltype(out->julday)>(sym("julday_"));
        out->caldat = reinterpret_cast<decltype(out->caldat)>(sym("caldat_"));
        out->get_doy = reinterpret_cast<decltype(out->get_doy)>(sym("get_doy_"));
        out->date2decy = reinterpret_cast<decltype(out->date2decy)>(sym("date_and_time2decy_"));
        out->decy2date = reinterpret_cast<decltype(out->decy2date)>(sym("decy2date_and_time_"));
        out->doyut2date =
            reinterpret_cast<decltype(out->doyut2date)>(sym("doy_and_ut2date_and_time_"));
        out->geo2gsm = reinterpret_cast<decltype(out->geo2gsm)>(sym("geo2gsm1_"));
        out->gsm2geo = reinterpret_cast<decltype(out->gsm2geo)>(sym("gsm2geo1_"));
        out->geo2gse = reinterpret_cast<decltype(out->geo2gse)>(sym("geo2gse1_"));
        out->geo2sm = reinterpret_cast<decltype(out->geo2sm)>(sym("geo2sm1_"));
        out->geo2gei = reinterpret_cast<decltype(out->geo2gei)>(sym("geo2gei1_"));
        out->gsm2sm = reinterpret_cast<decltype(out->gsm2sm)>(sym("gsm2sm1_"));
        out->geo2mag = reinterpret_cast<decltype(out->geo2mag)>(sym("geo2mag1_"));
        out->gdz_geo = reinterpret_cast<decltype(out->gdz_geo)>(sym("gdz_geo_"));
        out->geo_gdz = reinterpret_cast<decltype(out->geo_gdz)>(sym("geo_gdz_"));
        out->sph_car = reinterpret_cast<decltype(out->sph_car)>(sym("sph_car_"));
        out->car_sph = reinterpret_cast<decltype(out->car_sph)>(sym("car_sph_"));
        out->rll_gdz = reinterpret_cast<decltype(out->rll_gdz)>(sym("rll_gdz_"));
        out->hae2hee = reinterpret_cast<decltype(out->hae2hee)>(sym("hae2hee1_"));
        out->hee2hae = reinterpret_cast<decltype(out->hee2hae)>(sym("hee2hae1_"));
        out->hae2heeq = reinterpret_cast<decltype(out->hae2heeq)>(sym("hae2heeq1_"));
        out->heeq2hae = reinterpret_cast<decltype(out->heeq2hae)>(sym("heeq2hae1_"));
        out->gse2hee = reinterpret_cast<decltype(out->gse2hee)>(sym("gse2hee1_"));
        out->hee2gse = reinterpret_cast<decltype(out->hee2gse)>(sym("hee2gse1_"));
        out->get_field = reinterpret_cast<decltype(out->get_field)>(sym("get_field1_"));
        out->make_lstar = reinterpret_cast<decltype(out->make_lstar)>(sym("make_lstar1_"));
        out->make_lstar_ss =
            reinterpret_cast<decltype(out->make_lstar_ss)>(sym("make_lstar_shell_splitting1_"));
        return out;
    }();
    return *o;
}

/// Skip the benchmark, with a named reason, when the oracle is missing.
/// @param state the benchmark state to skip.
/// @param fn the entry point that had to resolve.
/// @return true when the benchmark must return immediately.
/// @complexity O(1).
/// @alloc the message string on the failure path only.
bool oracle_missing(benchmark::State& state, const void* fn) {
    if (fn != nullptr) return false;
    state.SkipWithError("IRBEM oracle not available — build it and point SPACE_IRBEM_ORACLE at the "
                        "-O2 library (default /tmp/irbem-builds/libirbem-O2.so)");
    return true;
}

// ---------------------------------------------------------------------------------------------
// datetime.hpp
// ---------------------------------------------------------------------------------------------

/// `julian_day_number` — the integer primitive every other calendar quantity is a difference of.
/// @param state the benchmark state.
/// @complexity O(kPoints) per iteration.
/// @alloc none.
void BM_cpu_julian_day_number(benchmark::State& state) {
    const Inputs& in = inputs();
    for (auto _ : state) {
        const DateInput* d = opaque(in.dates.data());
        for (std::size_t k = 0; k < kPoints; ++k) {
            benchmark::DoNotOptimize(julian_day_number(d[k].year, d[k].month, d[k].day));
        }
    }
    set_points(state, kPoints);
}
BENCHMARK(BM_cpu_julian_day_number);

/// IRBEM's `JULDAY`, the counterpart of @ref BM_cpu_julian_day_number.
/// @param state the benchmark state.
/// @complexity O(kPoints) per iteration.
/// @alloc none.
void BM_irbem_julday(benchmark::State& state) {
    const Oracle& o = oracle();
    if (oracle_missing(state, reinterpret_cast<const void*>(o.julday))) return;
    const Inputs& in = inputs();
    for (auto _ : state) {
        const DateInput* d = opaque(in.dates.data());
        for (std::size_t k = 0; k < kPoints; ++k) {
            int y = d[k].year;
            int m = d[k].month;
            int dd = d[k].day;
            benchmark::DoNotOptimize(o.julday(&y, &m, &dd));
        }
    }
    set_points(state, kPoints);
}
BENCHMARK(BM_irbem_julday);

/// `calendar_date` — the inverse, which this module makes total where IRBEM's is not.
/// @param state the benchmark state.
/// @complexity O(kPoints) per iteration.
/// @alloc none.
void BM_cpu_calendar_date(benchmark::State& state) {
    const Inputs& in = inputs();
    for (auto _ : state) {
        const DateInput* d = opaque(in.dates.data());
        for (std::size_t k = 0; k < kPoints; ++k) {
            benchmark::DoNotOptimize(calendar_date(d[k].jdn));
        }
    }
    set_points(state, kPoints);
}
BENCHMARK(BM_cpu_calendar_date);

/// IRBEM's `CALDAT`, the counterpart of @ref BM_cpu_calendar_date.
/// @param state the benchmark state.
/// @complexity O(kPoints) per iteration.
/// @alloc none.
void BM_irbem_caldat(benchmark::State& state) {
    const Oracle& o = oracle();
    if (oracle_missing(state, reinterpret_cast<const void*>(o.caldat))) return;
    const Inputs& in = inputs();
    for (auto _ : state) {
        const DateInput* d = opaque(in.dates.data());
        for (std::size_t k = 0; k < kPoints; ++k) {
            int jd = static_cast<int>(d[k].jdn);
            int y = 0;
            int m = 0;
            int dd = 0;
            o.caldat(&jd, &y, &m, &dd);
            benchmark::DoNotOptimize(y);
            benchmark::DoNotOptimize(m);
            benchmark::DoNotOptimize(dd);
        }
    }
    set_points(state, kPoints);
}
BENCHMARK(BM_irbem_caldat);

/// `day_of_year`.
/// @param state the benchmark state.
/// @complexity O(kPoints) per iteration.
/// @alloc none.
void BM_cpu_day_of_year(benchmark::State& state) {
    const Inputs& in = inputs();
    for (auto _ : state) {
        const DateInput* d = opaque(in.dates.data());
        for (std::size_t k = 0; k < kPoints; ++k) {
            benchmark::DoNotOptimize(day_of_year(d[k].year, d[k].month, d[k].day));
        }
    }
    set_points(state, kPoints);
}
BENCHMARK(BM_cpu_day_of_year);

/// IRBEM's `GET_DOY`, the counterpart of @ref BM_cpu_day_of_year.
/// @param state the benchmark state.
/// @complexity O(kPoints) per iteration.
/// @alloc none.
void BM_irbem_get_doy(benchmark::State& state) {
    const Oracle& o = oracle();
    if (oracle_missing(state, reinterpret_cast<const void*>(o.get_doy))) return;
    const Inputs& in = inputs();
    for (auto _ : state) {
        const DateInput* d = opaque(in.dates.data());
        for (std::size_t k = 0; k < kPoints; ++k) {
            int y = d[k].year;
            int m = d[k].month;
            int dd = d[k].day;
            benchmark::DoNotOptimize(o.get_doy(&y, &m, &dd));
        }
    }
    set_points(state, kPoints);
}
BENCHMARK(BM_irbem_get_doy);

/// `decimal_year` — the form the IGRF coefficient interpolation consumes.
/// @param state the benchmark state.
/// @complexity O(kPoints) per iteration.
/// @alloc none.
void BM_cpu_decimal_year(benchmark::State& state) {
    const Inputs& in = inputs();
    for (auto _ : state) {
        const DateInput* d = opaque(in.dates.data());
        for (std::size_t k = 0; k < kPoints; ++k) {
            benchmark::DoNotOptimize(decimal_year(d[k].year, d[k].month, d[k].day, d[k].hour,
                                                  d[k].minute, d[k].second));
        }
    }
    set_points(state, kPoints);
}
BENCHMARK(BM_cpu_decimal_year);

/// IRBEM's `DATE_AND_TIME2DECY`, the counterpart of @ref BM_cpu_decimal_year.
/// @param state the benchmark state.
/// @complexity O(kPoints) per iteration.
/// @alloc none.
void BM_irbem_date_and_time2decy(benchmark::State& state) {
    const Oracle& o = oracle();
    if (oracle_missing(state, reinterpret_cast<const void*>(o.date2decy))) return;
    const Inputs& in = inputs();
    for (auto _ : state) {
        const DateInput* d = opaque(in.dates.data());
        for (std::size_t k = 0; k < kPoints; ++k) {
            int y = d[k].year;
            int m = d[k].month;
            int dd = d[k].day;
            int h = d[k].hour;
            int mi = d[k].minute;
            int s = d[k].second;
            double decy = 0.0;
            o.date2decy(&y, &m, &dd, &h, &mi, &s, &decy);
            benchmark::DoNotOptimize(decy);
        }
    }
    set_points(state, kPoints);
}
BENCHMARK(BM_irbem_date_and_time2decy);

/// `date_and_time_from_decimal_year` — the millisecond-snapped inverse.
/// @param state the benchmark state.
/// @complexity O(kPoints) per iteration.
/// @alloc none.
void BM_cpu_date_and_time_from_decimal_year(benchmark::State& state) {
    const Inputs& in = inputs();
    for (auto _ : state) {
        const DateInput* d = opaque(in.dates.data());
        for (std::size_t k = 0; k < kPoints; ++k) {
            benchmark::DoNotOptimize(date_and_time_from_decimal_year(d[k].decy));
        }
    }
    set_points(state, kPoints);
}
BENCHMARK(BM_cpu_date_and_time_from_decimal_year);

/// IRBEM's `DECY2DATE_AND_TIME`, the counterpart of
/// @ref BM_cpu_date_and_time_from_decimal_year.
/// @param state the benchmark state.
/// @complexity O(kPoints) per iteration.
/// @alloc none.
void BM_irbem_decy2date_and_time(benchmark::State& state) {
    const Oracle& o = oracle();
    if (oracle_missing(state, reinterpret_cast<const void*>(o.decy2date))) return;
    const Inputs& in = inputs();
    for (auto _ : state) {
        const DateInput* d = opaque(in.dates.data());
        for (std::size_t k = 0; k < kPoints; ++k) {
            double decy = d[k].decy;
            int y = 0;
            int m = 0;
            int dd = 0;
            int doy = 0;
            int h = 0;
            int mi = 0;
            int s = 0;
            double ut = 0.0;
            o.decy2date(&decy, &y, &m, &dd, &doy, &h, &mi, &s, &ut);
            benchmark::DoNotOptimize(ut);
            benchmark::DoNotOptimize(doy);
        }
    }
    set_points(state, kPoints);
}
BENCHMARK(BM_irbem_decy2date_and_time);

/// `date_and_time_from_doy_and_ut`.
/// @param state the benchmark state.
/// @complexity O(kPoints) per iteration.
/// @alloc none.
void BM_cpu_date_and_time_from_doy_and_ut(benchmark::State& state) {
    const Inputs& in = inputs();
    for (auto _ : state) {
        const DateInput* d = opaque(in.dates.data());
        for (std::size_t k = 0; k < kPoints; ++k) {
            benchmark::DoNotOptimize(date_and_time_from_doy_and_ut(d[k].year, d[k].doy, d[k].ut));
        }
    }
    set_points(state, kPoints);
}
BENCHMARK(BM_cpu_date_and_time_from_doy_and_ut);

/// IRBEM's `DOY_AND_UT2DATE_AND_TIME`, the counterpart of
/// @ref BM_cpu_date_and_time_from_doy_and_ut.
/// @param state the benchmark state.
/// @complexity O(kPoints) per iteration.
/// @alloc none.
void BM_irbem_doy_and_ut2date_and_time(benchmark::State& state) {
    const Oracle& o = oracle();
    if (oracle_missing(state, reinterpret_cast<const void*>(o.doyut2date))) return;
    const Inputs& in = inputs();
    for (auto _ : state) {
        const DateInput* d = opaque(in.dates.data());
        for (std::size_t k = 0; k < kPoints; ++k) {
            int y = d[k].year;
            int doy = d[k].doy;
            double ut = d[k].ut;
            int m = 0;
            int dd = 0;
            int h = 0;
            int mi = 0;
            int s = 0;
            o.doyut2date(&y, &doy, &ut, &m, &dd, &h, &mi, &s);
            benchmark::DoNotOptimize(m);
            benchmark::DoNotOptimize(s);
        }
    }
    set_points(state, kPoints);
}
BENCHMARK(BM_irbem_doy_and_ut2date_and_time);

// ---------------------------------------------------------------------------------------------
// coords_rotations.hpp
// ---------------------------------------------------------------------------------------------

/// `Rotations::at` — the once-per-epoch build every "hot" transform amortizes over its ephemeris.
/// @param state the benchmark state.
/// @complexity O(kPoints) builds per iteration, each a fixed ~14 transcendentals.
/// @alloc none.
void BM_cpu_rotations_build(benchmark::State& state) {
    const Inputs& in = inputs();
    for (auto _ : state) {
        const DateInput* d = opaque(in.dates.data());
        for (std::size_t k = 0; k < kPoints; ++k) {
            benchmark::DoNotOptimize(Rotations::at(d[k].jd_ut1, in.dipole));
        }
    }
    set_points(state, kPoints);
}
BENCHMARK(BM_cpu_rotations_build);

/// A frame transform with the epoch's @ref Rotations already built — the per-point cost of an
/// ephemeris pass.
/// @tparam To the destination frame. @tparam From the source frame.
/// @param state the benchmark state.
/// @complexity O(kPoints) 3x3 products per iteration.
/// @alloc none.
template <Frame To, Frame From>
void bench_transform_hot(benchmark::State& state) {
    const Inputs& in = inputs();
    const auto& src = [&]() -> const std::array<Position<From>, kPoints>& {
        if constexpr (From == Frame::GEO) return in.geo;
        else return in.gsm;
    }();
    for (auto _ : state) {
        const Position<From>* p = opaque(src.data());
        for (std::size_t k = 0; k < kPoints; ++k) {
            benchmark::DoNotOptimize(transform<To>(p[k], in.rotations));
        }
    }
    set_points(state, kPoints);
}

/// The same transform with the @ref Rotations rebuilt for every point — the like-for-like against
/// IRBEM, whose entry points take a date with every call and cannot amortize.
/// @tparam To the destination frame. @tparam From the source frame.
/// @param state the benchmark state.
/// @complexity O(kPoints) epoch builds plus O(kPoints) products per iteration.
/// @alloc none.
template <Frame To, Frame From>
void bench_transform_cold(benchmark::State& state) {
    const Inputs& in = inputs();
    const auto& src = [&]() -> const std::array<Position<From>, kPoints>& {
        if constexpr (From == Frame::GEO) return in.geo;
        else return in.gsm;
    }();
    for (auto _ : state) {
        const Position<From>* p = opaque(src.data());
        const DateInput* d = opaque(in.dates.data());
        for (std::size_t k = 0; k < kPoints; ++k) {
            benchmark::DoNotOptimize(
                transform<To>(p[k], Rotations::at(d[k].jd_ut1, in.dipole)));
        }
    }
    set_points(state, kPoints);
}

BENCHMARK(bench_transform_hot<Frame::GSM, Frame::GEO>)->Name("BM_cpu_geo_to_gsm_hot");
BENCHMARK(bench_transform_cold<Frame::GSM, Frame::GEO>)->Name("BM_cpu_geo_to_gsm_cold");
BENCHMARK(bench_transform_hot<Frame::GSE, Frame::GEO>)->Name("BM_cpu_geo_to_gse_hot");
BENCHMARK(bench_transform_cold<Frame::GSE, Frame::GEO>)->Name("BM_cpu_geo_to_gse_cold");
BENCHMARK(bench_transform_hot<Frame::SM, Frame::GEO>)->Name("BM_cpu_geo_to_sm_hot");
BENCHMARK(bench_transform_cold<Frame::SM, Frame::GEO>)->Name("BM_cpu_geo_to_sm_cold");
BENCHMARK(bench_transform_hot<Frame::MAG, Frame::GEO>)->Name("BM_cpu_geo_to_mag_hot");
BENCHMARK(bench_transform_cold<Frame::MAG, Frame::GEO>)->Name("BM_cpu_geo_to_mag_cold");
BENCHMARK(bench_transform_hot<Frame::GEI, Frame::GEO>)->Name("BM_cpu_geo_to_gei_hot");
BENCHMARK(bench_transform_cold<Frame::GEI, Frame::GEO>)->Name("BM_cpu_geo_to_gei_cold");
BENCHMARK(bench_transform_hot<Frame::SM, Frame::GSM>)->Name("BM_cpu_gsm_to_sm_hot");
BENCHMARK(bench_transform_cold<Frame::SM, Frame::GSM>)->Name("BM_cpu_gsm_to_sm_cold");
// The transposed direction: GSM->GEO is not stored, it is `geo_to_gsm` read back through
// `fixarray::transpose`, so this row is what that costs.
BENCHMARK(bench_transform_hot<Frame::GEO, Frame::GSM>)->Name("BM_cpu_gsm_to_geo_hot");
BENCHMARK(bench_transform_cold<Frame::GEO, Frame::GSM>)->Name("BM_cpu_gsm_to_geo_cold");

/// An IRBEM `xIN -> xOUT` transform whose signature carries a `psi` out-parameter
/// (`geo2gsm1`, `gsm2geo1`).
/// @tparam Fn a pointer to the resolved entry point in @ref Oracle.
/// @param state the benchmark state.
/// @param fn the resolved entry point.
/// @param src the input positions.
/// @complexity O(kPoints) per iteration.
/// @alloc none.
void bench_irbem_psi(benchmark::State& state,
                     void (*fn)(int*, int*, double*, double*, double*, double*),
                     const Position<Frame::GEO>* src) {
    if (oracle_missing(state, reinterpret_cast<const void*>(fn))) return;
    const Inputs& in = inputs();
    for (auto _ : state) {
        const Position<Frame::GEO>* p = opaque(src);
        const DateInput* d = opaque(in.dates.data());
        for (std::size_t k = 0; k < kPoints; ++k) {
            int y = d[k].year;
            int doy = d[k].doy;
            double ut = d[k].ut;
            double psi = 0.0;
            std::array<double, 3> xi{p[k].v[0], p[k].v[1], p[k].v[2]};
            std::array<double, 3> xo{};
            fn(&y, &doy, &ut, &psi, xi.data(), xo.data());
            benchmark::DoNotOptimize(xo);
        }
    }
    set_points(state, kPoints);
}

/// An IRBEM `xIN -> xOUT` transform with no `psi` (`geo2gse1`, `geo2sm1`, `geo2gei1`, `gsm2sm1`).
/// @param state the benchmark state.
/// @param fn the resolved entry point.
/// @param src the input positions.
/// @complexity O(kPoints) per iteration.
/// @alloc none.
void bench_irbem_plain(benchmark::State& state,
                       void (*fn)(int*, int*, double*, double*, double*),
                       const Position<Frame::GEO>* src) {
    if (oracle_missing(state, reinterpret_cast<const void*>(fn))) return;
    const Inputs& in = inputs();
    for (auto _ : state) {
        const Position<Frame::GEO>* p = opaque(src);
        const DateInput* d = opaque(in.dates.data());
        for (std::size_t k = 0; k < kPoints; ++k) {
            int y = d[k].year;
            int doy = d[k].doy;
            double ut = d[k].ut;
            std::array<double, 3> xi{p[k].v[0], p[k].v[1], p[k].v[2]};
            std::array<double, 3> xo{};
            fn(&y, &doy, &ut, xi.data(), xo.data());
            benchmark::DoNotOptimize(xo);
        }
    }
    set_points(state, kPoints);
}

/// IRBEM's `geo2gsm1`.
/// @param state the benchmark state.
/// @complexity O(kPoints) per iteration.
/// @alloc none.
void BM_irbem_geo2gsm(benchmark::State& state) {
    bench_irbem_psi(state, oracle().geo2gsm, inputs().geo.data());
}
BENCHMARK(BM_irbem_geo2gsm);

/// IRBEM's `gsm2geo1` — its side of the transposed direction.
/// @param state the benchmark state.
/// @complexity O(kPoints) per iteration.
/// @alloc none.
void BM_irbem_gsm2geo(benchmark::State& state) {
    bench_irbem_psi(state, oracle().gsm2geo,
                    reinterpret_cast<const Position<Frame::GEO>*>(inputs().gsm.data()));
}
BENCHMARK(BM_irbem_gsm2geo);

/// IRBEM's `geo2gse1`.
/// @param state the benchmark state.
/// @complexity O(kPoints) per iteration.
/// @alloc none.
void BM_irbem_geo2gse(benchmark::State& state) {
    bench_irbem_plain(state, oracle().geo2gse, inputs().geo.data());
}
BENCHMARK(BM_irbem_geo2gse);

/// IRBEM's `geo2sm1`.
/// @param state the benchmark state.
/// @complexity O(kPoints) per iteration.
/// @alloc none.
void BM_irbem_geo2sm(benchmark::State& state) {
    bench_irbem_plain(state, oracle().geo2sm, inputs().geo.data());
}
BENCHMARK(BM_irbem_geo2sm);

/// IRBEM's `geo2gei1`.
/// @param state the benchmark state.
/// @complexity O(kPoints) per iteration.
/// @alloc none.
void BM_irbem_geo2gei(benchmark::State& state) {
    bench_irbem_plain(state, oracle().geo2gei, inputs().geo.data());
}
BENCHMARK(BM_irbem_geo2gei);

/// IRBEM's `gsm2sm1`.
/// @param state the benchmark state.
/// @complexity O(kPoints) per iteration.
/// @alloc none.
void BM_irbem_gsm2sm(benchmark::State& state) {
    bench_irbem_plain(state, oracle().gsm2sm,
                      reinterpret_cast<const Position<Frame::GEO>*>(inputs().gsm.data()));
}
BENCHMARK(BM_irbem_gsm2sm);

/// IRBEM's `geo2mag1`, whose signature takes only the year — its MAG frame is epoch-of-year, not
/// epoch-of-instant, which is why its number is lower than the other transforms'.
/// @param state the benchmark state.
/// @complexity O(kPoints) per iteration.
/// @alloc none.
void BM_irbem_geo2mag(benchmark::State& state) {
    const Oracle& o = oracle();
    if (oracle_missing(state, reinterpret_cast<const void*>(o.geo2mag))) return;
    const Inputs& in = inputs();
    for (auto _ : state) {
        const Position<Frame::GEO>* p = opaque(in.geo.data());
        const DateInput* d = opaque(in.dates.data());
        for (std::size_t k = 0; k < kPoints; ++k) {
            int y = d[k].year;
            std::array<double, 3> xi{p[k].v[0], p[k].v[1], p[k].v[2]};
            std::array<double, 3> xo{};
            o.geo2mag(&y, xi.data(), xo.data());
            benchmark::DoNotOptimize(xo);
        }
    }
    set_points(state, kPoints);
}
BENCHMARK(BM_irbem_geo2mag);

// ---------------------------------------------------------------------------------------------
// coords_geodetic.hpp — no epoch, so no hot/cold split
// ---------------------------------------------------------------------------------------------

/// `gdz_to_geo` — the ellipsoid forward map.
/// @param state the benchmark state.
/// @complexity O(kPoints) per iteration.
/// @alloc none.
void BM_cpu_gdz_to_geo(benchmark::State& state) {
    const Inputs& in = inputs();
    for (auto _ : state) {
        const Position<Frame::GDZ>* p = opaque(in.gdz.data());
        for (std::size_t k = 0; k < kPoints; ++k) benchmark::DoNotOptimize(gdz_to_geo(p[k]));
    }
    set_points(state, kPoints);
}
BENCHMARK(BM_cpu_gdz_to_geo);

/// IRBEM's `gdz_geo`, the counterpart of @ref BM_cpu_gdz_to_geo.
/// @param state the benchmark state.
/// @complexity O(kPoints) per iteration.
/// @alloc none.
void BM_irbem_gdz_geo(benchmark::State& state) {
    const Oracle& o = oracle();
    if (oracle_missing(state, reinterpret_cast<const void*>(o.gdz_geo))) return;
    const Inputs& in = inputs();
    for (auto _ : state) {
        const Position<Frame::GDZ>* p = opaque(in.gdz.data());
        for (std::size_t k = 0; k < kPoints; ++k) {
            double alt = p[k].v[0];
            double lat = p[k].v[1];
            double lon = p[k].v[2];
            double x = 0.0;
            double y = 0.0;
            double z = 0.0;
            o.gdz_geo(&lat, &lon, &alt, &x, &y, &z);
            benchmark::DoNotOptimize(x);
            benchmark::DoNotOptimize(y);
            benchmark::DoNotOptimize(z);
        }
    }
    set_points(state, kPoints);
}
BENCHMARK(BM_irbem_gdz_geo);

/// `geo_to_gdz` — the Bowring inverse, four fixed iterations.
/// @param state the benchmark state.
/// @complexity O(kPoints) per iteration.
/// @alloc none.
void BM_cpu_geo_to_gdz(benchmark::State& state) {
    const Inputs& in = inputs();
    for (auto _ : state) {
        const Position<Frame::GEO>* p = opaque(in.geo.data());
        for (std::size_t k = 0; k < kPoints; ++k) benchmark::DoNotOptimize(geo_to_gdz(p[k]));
    }
    set_points(state, kPoints);
}
BENCHMARK(BM_cpu_geo_to_gdz);

/// IRBEM's `geo_gdz`, the counterpart of @ref BM_cpu_geo_to_gdz.
/// @param state the benchmark state.
/// @complexity O(kPoints) per iteration.
/// @alloc none.
void BM_irbem_geo_gdz(benchmark::State& state) {
    const Oracle& o = oracle();
    if (oracle_missing(state, reinterpret_cast<const void*>(o.geo_gdz))) return;
    const Inputs& in = inputs();
    for (auto _ : state) {
        const Position<Frame::GEO>* p = opaque(in.geo.data());
        for (std::size_t k = 0; k < kPoints; ++k) {
            double x = p[k].v[0];
            double y = p[k].v[1];
            double z = p[k].v[2];
            double lat = 0.0;
            double lon = 0.0;
            double alt = 0.0;
            o.geo_gdz(&x, &y, &z, &lat, &lon, &alt);
            benchmark::DoNotOptimize(lat);
            benchmark::DoNotOptimize(alt);
        }
    }
    set_points(state, kPoints);
}
BENCHMARK(BM_irbem_geo_gdz);

/// `sph_to_car`.
/// @param state the benchmark state.
/// @complexity O(kPoints) per iteration.
/// @alloc none.
void BM_cpu_sph_to_car(benchmark::State& state) {
    const Inputs& in = inputs();
    for (auto _ : state) {
        const Position<Frame::SPH>* p = opaque(in.sph.data());
        for (std::size_t k = 0; k < kPoints; ++k) benchmark::DoNotOptimize(sph_to_car(p[k]));
    }
    set_points(state, kPoints);
}
BENCHMARK(BM_cpu_sph_to_car);

/// IRBEM's `SPH_CAR`, the counterpart of @ref BM_cpu_sph_to_car.
/// @param state the benchmark state.
/// @complexity O(kPoints) per iteration.
/// @alloc none.
void BM_irbem_sph_car(benchmark::State& state) {
    const Oracle& o = oracle();
    if (oracle_missing(state, reinterpret_cast<const void*>(o.sph_car))) return;
    const Inputs& in = inputs();
    for (auto _ : state) {
        const Position<Frame::SPH>* p = opaque(in.sph.data());
        for (std::size_t k = 0; k < kPoints; ++k) {
            double r = p[k].v[0];
            double lat = p[k].v[1];
            double lon = p[k].v[2];
            std::array<double, 3> x{};
            o.sph_car(&r, &lat, &lon, x.data());
            benchmark::DoNotOptimize(x);
        }
    }
    set_points(state, kPoints);
}
BENCHMARK(BM_irbem_sph_car);

/// `car_to_sph`.
/// @param state the benchmark state.
/// @complexity O(kPoints) per iteration.
/// @alloc none.
void BM_cpu_car_to_sph(benchmark::State& state) {
    const Inputs& in = inputs();
    for (auto _ : state) {
        const Position<Frame::GEO>* p = opaque(in.geo.data());
        for (std::size_t k = 0; k < kPoints; ++k) benchmark::DoNotOptimize(car_to_sph(p[k]));
    }
    set_points(state, kPoints);
}
BENCHMARK(BM_cpu_car_to_sph);

/// IRBEM's `CAR_SPH`, the counterpart of @ref BM_cpu_car_to_sph.
/// @param state the benchmark state.
/// @complexity O(kPoints) per iteration.
/// @alloc none.
void BM_irbem_car_sph(benchmark::State& state) {
    const Oracle& o = oracle();
    if (oracle_missing(state, reinterpret_cast<const void*>(o.car_sph))) return;
    const Inputs& in = inputs();
    for (auto _ : state) {
        const Position<Frame::GEO>* p = opaque(in.geo.data());
        for (std::size_t k = 0; k < kPoints; ++k) {
            std::array<double, 3> x{p[k].v[0], p[k].v[1], p[k].v[2]};
            double r = 0.0;
            double lat = 0.0;
            double lon = 0.0;
            o.car_sph(x.data(), &r, &lat, &lon);
            benchmark::DoNotOptimize(r);
            benchmark::DoNotOptimize(lat);
        }
    }
    set_points(state, kPoints);
}
BENCHMARK(BM_irbem_car_sph);

/// `rll_to_gdz`.
/// @param state the benchmark state.
/// @complexity O(kPoints) per iteration.
/// @alloc none.
void BM_cpu_rll_to_gdz(benchmark::State& state) {
    const Inputs& in = inputs();
    for (auto _ : state) {
        const Position<Frame::RLL>* p = opaque(in.rll.data());
        for (std::size_t k = 0; k < kPoints; ++k) benchmark::DoNotOptimize(rll_to_gdz(p[k]));
    }
    set_points(state, kPoints);
}
BENCHMARK(BM_cpu_rll_to_gdz);

/// IRBEM's `RLL_GDZ`, the counterpart of @ref BM_cpu_rll_to_gdz.
/// @param state the benchmark state.
/// @complexity O(kPoints) per iteration.
/// @alloc none.
void BM_irbem_rll_gdz(benchmark::State& state) {
    const Oracle& o = oracle();
    if (oracle_missing(state, reinterpret_cast<const void*>(o.rll_gdz))) return;
    const Inputs& in = inputs();
    for (auto _ : state) {
        const Position<Frame::RLL>* p = opaque(in.rll.data());
        for (std::size_t k = 0; k < kPoints; ++k) {
            double r = p[k].v[0];
            double lat = p[k].v[1];
            double lon = p[k].v[2];
            double alt = 0.0;
            o.rll_gdz(&r, &lat, &lon, &alt);
            benchmark::DoNotOptimize(alt);
        }
    }
    set_points(state, kPoints);
}
BENCHMARK(BM_irbem_rll_gdz);

// ---------------------------------------------------------------------------------------------
// coords_helio.hpp
// ---------------------------------------------------------------------------------------------

/// `helio_geometry` — the once-per-epoch heliospheric build.
/// @param state the benchmark state.
/// @complexity O(kPoints) builds per iteration.
/// @alloc none.
void BM_cpu_helio_geometry(benchmark::State& state) {
    const Inputs& in = inputs();
    for (auto _ : state) {
        const DateInput* d = opaque(in.dates.data());
        for (std::size_t k = 0; k < kPoints; ++k) {
            benchmark::DoNotOptimize(helio_geometry(d[k].jd_ut1 - 2400000.5));
        }
    }
    set_points(state, kPoints);
}
BENCHMARK(BM_cpu_helio_geometry);

/// A heliospheric transform with the epoch's @ref HelioGeometry already built.
/// @tparam In the input position type. @tparam Out the output position type.
/// @tparam Fn the transform.
/// @param state the benchmark state.
/// @param src the input positions.
/// @complexity O(kPoints) per iteration.
/// @alloc none.
template <class In, class Out, Out (*Fn)(const In&, const HelioGeometry&)>
void bench_helio_hot(benchmark::State& state, const In* src) {
    const Inputs& in = inputs();
    for (auto _ : state) {
        const In* p = opaque(src);
        for (std::size_t k = 0; k < kPoints; ++k) benchmark::DoNotOptimize(Fn(p[k], in.helio));
    }
    set_points(state, kPoints);
}

/// The same transform with the @ref HelioGeometry rebuilt per point — the like-for-like against
/// IRBEM's `hae2hee1` and friends, which take a date with every call.
/// @tparam In the input position type. @tparam Out the output position type.
/// @tparam Fn the transform.
/// @param state the benchmark state.
/// @param src the input positions.
/// @complexity O(kPoints) geometry builds plus O(kPoints) products per iteration.
/// @alloc none.
template <class In, class Out, Out (*Fn)(const In&, const HelioGeometry&)>
void bench_helio_cold(benchmark::State& state, const In* src) {
    const Inputs& in = inputs();
    for (auto _ : state) {
        const In* p = opaque(src);
        const DateInput* d = opaque(in.dates.data());
        for (std::size_t k = 0; k < kPoints; ++k) {
            benchmark::DoNotOptimize(Fn(p[k], helio_geometry(d[k].jd_ut1 - 2400000.5)));
        }
    }
    set_points(state, kPoints);
}

/// The `Position` instantiation of `HAE2HEE`, as a plain function pointer.
constexpr auto kHae2Hee = static_cast<Position<Frame::HEE> (*)(const Position<Frame::HAE>&,
                                                               const HelioGeometry&)>(
    &HAE2HEE<Position>);
/// The `Position` instantiation of `HEE2HAE`.
constexpr auto kHee2Hae = static_cast<Position<Frame::HAE> (*)(const Position<Frame::HEE>&,
                                                               const HelioGeometry&)>(
    &HEE2HAE<Position>);
/// The `Position` instantiation of `HAE2HEEQ`.
constexpr auto kHae2Heeq = static_cast<Position<Frame::HEEQ> (*)(const Position<Frame::HAE>&,
                                                                 const HelioGeometry&)>(
    &HAE2HEEQ<Position>);
/// The `Position` instantiation of `HEEQ2HAE`.
constexpr auto kHeeq2Hae = static_cast<Position<Frame::HAE> (*)(const Position<Frame::HEEQ>&,
                                                                const HelioGeometry&)>(
    &HEEQ2HAE<Position>);
/// The `Position` overload of `GSE2HEE` — the one that also shifts the origin.
constexpr auto kGse2Hee = static_cast<Position<Frame::HEE> (*)(const Position<Frame::GSE>&,
                                                               const HelioGeometry&)>(&GSE2HEE);
/// The `Position` overload of `HEE2GSE`.
constexpr auto kHee2Gse = static_cast<Position<Frame::GSE> (*)(const Position<Frame::HEE>&,
                                                               const HelioGeometry&)>(&HEE2GSE);

/// `HAE2HEE` for a position, with the epoch's HelioGeometry prebuilt.
/// @param state the benchmark state.
/// @complexity O(kPoints) per iteration. @alloc none.
void BM_cpu_hae_to_hee_hot(benchmark::State& state) {
    bench_helio_hot<Position<Frame::HAE>, Position<Frame::HEE>, kHae2Hee>(state, inputs().hae.data());
}
BENCHMARK(BM_cpu_hae_to_hee_hot);

/// `HAE2HEE` for a position, rebuilding the HelioGeometry per point.
/// @param state the benchmark state.
/// @complexity O(kPoints) per iteration. @alloc none.
void BM_cpu_hae_to_hee_cold(benchmark::State& state) {
    bench_helio_cold<Position<Frame::HAE>, Position<Frame::HEE>, kHae2Hee>(state, inputs().hae.data());
}
BENCHMARK(BM_cpu_hae_to_hee_cold);

/// `HEE2HAE` for a position, with the epoch's HelioGeometry prebuilt.
/// @param state the benchmark state.
/// @complexity O(kPoints) per iteration. @alloc none.
void BM_cpu_hee_to_hae_hot(benchmark::State& state) {
    bench_helio_hot<Position<Frame::HEE>, Position<Frame::HAE>, kHee2Hae>(state, inputs().hee.data());
}
BENCHMARK(BM_cpu_hee_to_hae_hot);

/// `HEE2HAE` for a position, rebuilding the HelioGeometry per point.
/// @param state the benchmark state.
/// @complexity O(kPoints) per iteration. @alloc none.
void BM_cpu_hee_to_hae_cold(benchmark::State& state) {
    bench_helio_cold<Position<Frame::HEE>, Position<Frame::HAE>, kHee2Hae>(state, inputs().hee.data());
}
BENCHMARK(BM_cpu_hee_to_hae_cold);

/// `HAE2HEEQ` for a position, with the epoch's HelioGeometry prebuilt.
/// @param state the benchmark state.
/// @complexity O(kPoints) per iteration. @alloc none.
void BM_cpu_hae_to_heeq_hot(benchmark::State& state) {
    bench_helio_hot<Position<Frame::HAE>, Position<Frame::HEEQ>, kHae2Heeq>(state, inputs().hae.data());
}
BENCHMARK(BM_cpu_hae_to_heeq_hot);

/// `HAE2HEEQ` for a position, rebuilding the HelioGeometry per point.
/// @param state the benchmark state.
/// @complexity O(kPoints) per iteration. @alloc none.
void BM_cpu_hae_to_heeq_cold(benchmark::State& state) {
    bench_helio_cold<Position<Frame::HAE>, Position<Frame::HEEQ>, kHae2Heeq>(state, inputs().hae.data());
}
BENCHMARK(BM_cpu_hae_to_heeq_cold);

/// `HEEQ2HAE` for a position, with the epoch's HelioGeometry prebuilt.
/// @param state the benchmark state.
/// @complexity O(kPoints) per iteration. @alloc none.
void BM_cpu_heeq_to_hae_hot(benchmark::State& state) {
    bench_helio_hot<Position<Frame::HEEQ>, Position<Frame::HAE>, kHeeq2Hae>(state, inputs().heeq.data());
}
BENCHMARK(BM_cpu_heeq_to_hae_hot);

/// `HEEQ2HAE` for a position, rebuilding the HelioGeometry per point.
/// @param state the benchmark state.
/// @complexity O(kPoints) per iteration. @alloc none.
void BM_cpu_heeq_to_hae_cold(benchmark::State& state) {
    bench_helio_cold<Position<Frame::HEEQ>, Position<Frame::HAE>, kHeeq2Hae>(state, inputs().heeq.data());
}
BENCHMARK(BM_cpu_heeq_to_hae_cold);

/// `GSE2HEE` for a position, with the epoch's HelioGeometry prebuilt.
/// @param state the benchmark state.
/// @complexity O(kPoints) per iteration. @alloc none.
void BM_cpu_gse_to_hee_hot(benchmark::State& state) {
    bench_helio_hot<Position<Frame::GSE>, Position<Frame::HEE>, kGse2Hee>(state, inputs().gse.data());
}
BENCHMARK(BM_cpu_gse_to_hee_hot);

/// `GSE2HEE` for a position, rebuilding the HelioGeometry per point.
/// @param state the benchmark state.
/// @complexity O(kPoints) per iteration. @alloc none.
void BM_cpu_gse_to_hee_cold(benchmark::State& state) {
    bench_helio_cold<Position<Frame::GSE>, Position<Frame::HEE>, kGse2Hee>(state, inputs().gse.data());
}
BENCHMARK(BM_cpu_gse_to_hee_cold);

/// `HEE2GSE` for a position, with the epoch's HelioGeometry prebuilt.
/// @param state the benchmark state.
/// @complexity O(kPoints) per iteration. @alloc none.
void BM_cpu_hee_to_gse_hot(benchmark::State& state) {
    bench_helio_hot<Position<Frame::HEE>, Position<Frame::GSE>, kHee2Gse>(state, inputs().hee.data());
}
BENCHMARK(BM_cpu_hee_to_gse_hot);

/// `HEE2GSE` for a position, rebuilding the HelioGeometry per point.
/// @param state the benchmark state.
/// @complexity O(kPoints) per iteration. @alloc none.
void BM_cpu_hee_to_gse_cold(benchmark::State& state) {
    bench_helio_cold<Position<Frame::HEE>, Position<Frame::GSE>, kHee2Gse>(state, inputs().hee.data());
}
BENCHMARK(BM_cpu_hee_to_gse_cold);

/// One of IRBEM's `(iyr, idoy, UT, xIN, xOUT)` heliospheric transforms.
/// @param state the benchmark state.
/// @param fn the resolved entry point.
/// @param src the input components (AU or Re, per the routine).
/// @complexity O(kPoints) per iteration.
/// @alloc none.
void bench_irbem_helio(benchmark::State& state, void (*fn)(int*, int*, double*, double*, double*),
                       const fx::vec3d* src) {
    if (oracle_missing(state, reinterpret_cast<const void*>(fn))) return;
    const Inputs& in = inputs();
    for (auto _ : state) {
        const fx::vec3d* p = opaque(src);
        const DateInput* d = opaque(in.dates.data());
        for (std::size_t k = 0; k < kPoints; ++k) {
            int y = d[k].year;
            int doy = d[k].doy;
            double ut = d[k].ut;
            std::array<double, 3> xi{p[k][0], p[k][1], p[k][2]};
            std::array<double, 3> xo{};
            fn(&y, &doy, &ut, xi.data(), xo.data());
            benchmark::DoNotOptimize(xo);
        }
    }
    set_points(state, kPoints);
}

/// IRBEM's `hae2hee1`.
/// @param state the benchmark state.
/// @complexity O(kPoints) per iteration. @alloc none.
void BM_irbem_hae2hee(benchmark::State& state) {
    bench_irbem_helio(state, oracle().hae2hee,
                      reinterpret_cast<const fx::vec3d*>(inputs().hae.data()));
}
BENCHMARK(BM_irbem_hae2hee);

/// IRBEM's `hee2hae1`.
/// @param state the benchmark state.
/// @complexity O(kPoints) per iteration. @alloc none.
void BM_irbem_hee2hae(benchmark::State& state) {
    bench_irbem_helio(state, oracle().hee2hae,
                      reinterpret_cast<const fx::vec3d*>(inputs().hee.data()));
}
BENCHMARK(BM_irbem_hee2hae);

/// IRBEM's `hae2heeq1`.
/// @param state the benchmark state.
/// @complexity O(kPoints) per iteration. @alloc none.
void BM_irbem_hae2heeq(benchmark::State& state) {
    bench_irbem_helio(state, oracle().hae2heeq,
                      reinterpret_cast<const fx::vec3d*>(inputs().hae.data()));
}
BENCHMARK(BM_irbem_hae2heeq);

/// IRBEM's `heeq2hae1`.
/// @param state the benchmark state.
/// @complexity O(kPoints) per iteration. @alloc none.
void BM_irbem_heeq2hae(benchmark::State& state) {
    bench_irbem_helio(state, oracle().heeq2hae,
                      reinterpret_cast<const fx::vec3d*>(inputs().heeq.data()));
}
BENCHMARK(BM_irbem_heeq2hae);

/// IRBEM's `gse2hee1`.
/// @param state the benchmark state.
/// @complexity O(kPoints) per iteration. @alloc none.
void BM_irbem_gse2hee(benchmark::State& state) {
    bench_irbem_helio(state, oracle().gse2hee,
                      reinterpret_cast<const fx::vec3d*>(inputs().gse.data()));
}
BENCHMARK(BM_irbem_gse2hee);

/// IRBEM's `hee2gse1`.
/// @param state the benchmark state.
/// @complexity O(kPoints) per iteration. @alloc none.
void BM_irbem_hee2gse(benchmark::State& state) {
    bench_irbem_helio(state, oracle().hee2gse,
                      reinterpret_cast<const fx::vec3d*>(inputs().hee.data()));
}
BENCHMARK(BM_irbem_hee2gse);

// ---------------------------------------------------------------------------------------------
// igrf.hpp
// ---------------------------------------------------------------------------------------------

/// `Igrf<13>::at` — the per-epoch coefficient interpolation, one pass over 105 slots.
/// @param state the benchmark state.
/// @complexity O(kPoints) interpolations per iteration, each O(NMAX^2).
/// @alloc none.
void BM_cpu_igrf_at_deg13(benchmark::State& state) {
    const Inputs& in = inputs();
    for (auto _ : state) {
        const DateInput* d = opaque(in.dates.data());
        for (std::size_t k = 0; k < kPoints; ++k) {
            benchmark::DoNotOptimize(Igrf<13>::at(1901.0 + (d[k].decy - std::floor(d[k].decy)) +
                                                  static_cast<double>(k)));
        }
    }
    set_points(state, kPoints);
}
BENCHMARK(BM_cpu_igrf_at_deg13);

/// `Igrf<NMAX>::evaluate` on a GEO Cartesian position — the kernel of every field-line trace.
/// @tparam NMAX the truncation degree.
/// @param state the benchmark state.
/// @complexity O(kPoints * NMAX^2) per iteration.
/// @alloc none.
template <int NMAX>
void bench_igrf_evaluate(benchmark::State& state) {
    const Inputs& in = inputs();
    const auto& model = [&]() -> const Igrf<NMAX>& {
        if constexpr (NMAX == 13) return *in.igrf13;
        else return *in.igrf10;
    }();
    for (auto _ : state) {
        const Position<Frame::GEO>* p = opaque(in.geo.data());
        for (std::size_t k = 0; k < kPoints; ++k) benchmark::DoNotOptimize(model.evaluate(p[k]));
    }
    set_points(state, kPoints);
}
BENCHMARK(bench_igrf_evaluate<13>)->Name("BM_cpu_igrf_evaluate_deg13");
BENCHMARK(bench_igrf_evaluate<10>)->Name("BM_cpu_igrf_evaluate_deg10");

/// `Igrf<13>::evaluate` on a SPHERICAL position — the same harmonic sum with the spherical-to-
/// Cartesian conversion in front of it.
/// @param state the benchmark state.
/// @complexity O(kPoints * 13^2) per iteration.
/// @alloc none.
void BM_cpu_igrf_evaluate_deg13_sph(benchmark::State& state) {
    const Inputs& in = inputs();
    for (auto _ : state) {
        const Position<Frame::SPH>* p = opaque(in.sph.data());
        for (std::size_t k = 0; k < kPoints; ++k) {
            benchmark::DoNotOptimize(in.igrf13->evaluate(p[k]));
        }
    }
    set_points(state, kPoints);
}
BENCHMARK(BM_cpu_igrf_evaluate_deg13_sph);

/// IRBEM's internal-field evaluation through `get_field1` — `options(5)` selects which model.
///
/// This is NOT a pure harmonic-sum comparison and the generated table says so: `get_field1` also
/// resolves `sysaxes`, converts the epoch, and (every `options(2)` days) re-interpolates the IGRF
/// coefficients. It is nevertheless the only published entry point that evaluates the internal
/// field at a point, so it is the honest counterpart — an upper bound on IRBEM's cost for the same
/// answer, not a lower one.
/// @param state the benchmark state.
/// @param internal_model the `options(5)` code: 0 = IGRF, 5 = centred dipole.
/// @complexity O(kPoints) per iteration.
/// @alloc none.
void bench_irbem_get_field(benchmark::State& state, int internal_model) {
    const Oracle& o = oracle();
    if (oracle_missing(state, reinterpret_cast<const void*>(o.get_field))) return;
    const Inputs& in = inputs();
    std::array<double, 25> maginput{};
    for (auto _ : state) {
        const Position<Frame::GEO>* p = opaque(in.geo.data());
        const DateInput* d = opaque(in.dates.data());
        for (std::size_t k = 0; k < kPoints; ++k) {
            int kext = 0;  // internal field only
            std::array<int, 5> options{0, 0, 0, 0, internal_model};
            int sysaxes = 1;  // GEO Cartesian, Re
            int y = d[k].year;
            int doy = d[k].doy;
            double ut = d[k].ut;
            double x1 = p[k].v[0];
            double x2 = p[k].v[1];
            double x3 = p[k].v[2];
            std::array<double, 3> bgeo{};
            double bl = 0.0;
            o.get_field(&kext, options.data(), &sysaxes, &y, &doy, &ut, &x1, &x2, &x3,
                        maginput.data(), bgeo.data(), &bl);
            benchmark::DoNotOptimize(bgeo);
            benchmark::DoNotOptimize(bl);
        }
    }
    set_points(state, kPoints);
}
BENCHMARK_CAPTURE(bench_irbem_get_field, igrf, 0)->Name("BM_irbem_get_field_igrf");
BENCHMARK_CAPTURE(bench_irbem_get_field, dipole, 5)->Name("BM_irbem_get_field_dipole");

// ---------------------------------------------------------------------------------------------
// gpu/dispatch.hpp — the centred dipole, on both lanes
// ---------------------------------------------------------------------------------------------

/// `dipole_field_at` — the fp64 scalar reference the device lane is measured against.
/// @param state the benchmark state.
/// @complexity O(kPoints) per iteration.
/// @alloc none.
void BM_cpu_dipole_field_at(benchmark::State& state) {
    const Inputs& in = inputs();
    const double g10 = in.dipole.g10;
    for (auto _ : state) {
        const Position<Frame::GEO>* p = opaque(in.geo.data());
        for (std::size_t k = 0; k < kPoints; ++k) {
            benchmark::DoNotOptimize(
                gpu::dipole_field_at(Position<Frame::MAG>{p[k].v}, g10));
        }
    }
    set_points(state, kPoints);
}
BENCHMARK(BM_cpu_dipole_field_at);

/// `dipole_field_host` — the fp32 batch lane, the twin of the kernel and the lane a machine with
/// no device actually runs.
/// @param state the benchmark state.
/// @complexity O(kBatch) per iteration.
/// @alloc none; both spans are sized once at start-up.
void BM_cpu_dipole_field_host_batch(benchmark::State& state) {
    const Inputs& in = inputs();
    auto& out = const_cast<std::vector<float>&>(in.batch_out);
    const auto g10 = static_cast<float>(in.dipole.g10);
    for (auto _ : state) {
        gpu::dipole_field_host(std::span<const float>(in.batch_pos), std::span<float>(out), g10);
        benchmark::DoNotOptimize(out.data());
        benchmark::ClobberMemory();
    }
    set_points(state, kBatch);
}
BENCHMARK(BM_cpu_dipole_field_host_batch)->Unit(benchmark::kMillisecond)->UseRealTime();

/// `dipole_field_gpu` — the device lane through `dispatch.hpp`, transfers INCLUDED.
///
/// Transfers are inside the timed region on purpose. A kernel time that excludes the bus is a
/// number about the ALU and not about the operation a caller can actually perform; `dispatch.hpp`'s
/// own registry comment records the device losing this particular race precisely because the bus is
/// where the dipole's twelve flops per point are spent.
/// @param state the benchmark state.
/// @complexity O(kBatch) device work, plus 2*3*kBatch floats over the bus, per iteration.
/// @alloc `dispatch_batch`'s four pooled device buffers per launch; nothing per point.
void BM_gpu_dipole_field_batch(benchmark::State& state) {
    const Inputs& in = inputs();
    auto& out = const_cast<std::vector<float>&>(in.batch_out);
    const auto g10 = static_cast<float>(in.dipole.g10);
    if (!gpu::available()) {
        state.SkipWithError(gpu::unavailable_reason().c_str());
        return;
    }
    try {
        // Warm the pipeline cache and the buffer pool so the first timed iteration is not the
        // shader compile.
        gpu::dipole_field_gpu(std::span<const float>(in.batch_pos), std::span<float>(out), g10);
    } catch (const std::exception& e) {
        state.SkipWithError(e.what());
        return;
    }
    for (auto _ : state) {
        gpu::dipole_field_gpu(std::span<const float>(in.batch_pos), std::span<float>(out), g10);
        benchmark::DoNotOptimize(out.data());
        benchmark::ClobberMemory();
    }
    set_points(state, kBatch);
}
BENCHMARK(BM_gpu_dipole_field_batch)->Unit(benchmark::kMillisecond)->UseRealTime();

// ---------------------------------------------------------------------------------------------
// ext_t89.hpp — the first routine this module has that is PARAMETERIZED BY ACTIVITY
// ---------------------------------------------------------------------------------------------

/// The dipole tilt every T89 lane is measured at, radians.
///
/// ~20°, which is a real northern-summer tilt rather than the ψ = 0 special case: at ψ = 0 the
/// `sin ψ` terms of Tsyganenko (1989) eq. (2)–(10) drop out and the model evaluates a strictly
/// cheaper expression than it ever does in flight. Benchmarking the tilt-free case would be
/// benchmarking a branch of the model no orbit sits on.
constexpr double kT89TiltRad = 0.35;

/// `t89_field` at one Kp — the storm sweep, and the reason this section exists.
///
/// Every other row of this suite is measured at quiet conditions, because until `ext_t89.hpp`
/// landed there was nothing in the module that took an activity index at all. There is now, and a
/// benchmark suite for a radiation-belt library that only ever measures `maginput = 0` is measuring
/// the one condition the library does not exist for. `state.range(0)` is Kp in IRBEM's `maginput`
/// slot-1 scaling — Kp × 10, 0…90 (`docs/source/api/general_information.rst`) — swept across the
/// whole published envelope, so the extreme-storm bin is measured, not assumed to cost the same.
///
/// T89 selects one of seven coefficient sets by Kp and then evaluates the same closed form, so the
/// sweep is expected to be FLAT in time. That is the point: a flat sweep is a measured result about
/// this model, published, rather than an assumption that let nine tenths of the parameter space go
/// untimed. A model whose cost does vary with activity — a T96 or a T01 root-find — would show it
/// here.
/// @param state the benchmark state; `range(0)` is Kp × 10.
/// @complexity O(kPoints) per iteration.
/// @alloc none.
void BM_cpu_t89_field(benchmark::State& state) {
    const Inputs& in = inputs();
    const auto kp10 = static_cast<double>(state.range(0));
    for (auto _ : state) {
        const Position<Frame::GSM>* p = opaque(in.gsm.data());
        for (std::size_t k = 0; k < kPoints; ++k) {
            benchmark::DoNotOptimize(t89_field(p[k], kT89TiltRad, kp10));
        }
    }
    set_points(state, kPoints);
}
BENCHMARK(BM_cpu_t89_field)->DenseRange(0, 90, 10);

/// IRBEM's `get_field1` with `kext = 4` — Tsyganenko (1989c), the counterpart of
/// @ref BM_cpu_t89_field, at the same Kp.
///
/// `kext = 4` is the documented selector for T89 and it "uses 0 ≤ Kp ≤ 9, valid for rGEO ≤ 70 Re"
/// (`docs/source/api/general_information.rst`); `maginput(1)` carries Kp × 10. The internal field is
/// left at IGRF (`options(5) = 0`) on both sides, so the difference between the two rows is the
/// external model and the frame work around it, which is what a caller pays.
/// @param state the benchmark state; `range(0)` is Kp × 10.
/// @complexity O(kPoints) per iteration.
/// @alloc none inside the loop.
void BM_irbem_get_field_t89(benchmark::State& state) {
    const Oracle& o = oracle();
    if (oracle_missing(state, reinterpret_cast<const void*>(o.get_field))) return;
    const Inputs& in = inputs();
    std::array<double, 25> maginput{};
    maginput[0] = static_cast<double>(state.range(0));   // Kp x 10, slot 1
    for (auto _ : state) {
        const Position<Frame::GSM>* p = opaque(in.gsm.data());
        const DateInput* d = opaque(in.dates.data());
        for (std::size_t k = 0; k < kPoints; ++k) {
            int kext = 4;  // Tsyganenko [1989c]
            std::array<int, 5> options{0, 0, 0, 0, 0};
            int sysaxes = 2;  // GSM Cartesian, Re
            int y = d[k].year;
            int doy = d[k].doy;
            double ut = d[k].ut;
            double x1 = p[k].v[0];
            double x2 = p[k].v[1];
            double x3 = p[k].v[2];
            std::array<double, 3> bgeo{};
            double bl = 0.0;
            o.get_field(&kext, options.data(), &sysaxes, &y, &doy, &ut, &x1, &x2, &x3,
                        maginput.data(), bgeo.data(), &bl);
            benchmark::DoNotOptimize(bgeo);
            benchmark::DoNotOptimize(bl);
        }
    }
    set_points(state, kPoints);
}
BENCHMARK(BM_irbem_get_field_t89)->DenseRange(0, 90, 10);

// ---------------------------------------------------------------------------------------------
// lstar.hpp — the trace, the second invariant, and the CROSSOVER the module argues from
// ---------------------------------------------------------------------------------------------

/// The largest batch the crossover sweep runs, and therefore how many distinct field lines the
/// input set holds.
///
/// 2^16 rather than something rounder because the sweep is a geometric one — 64, 256, 1024, 4096,
/// 16384, 65536 — and the top of it is where the device's advantage is measured. The curve is still
/// climbing there, which is a statement about the device not yet being saturated and is why the
/// sweep stops at a size rather than at a plateau.
constexpr std::size_t kTraceMax = 1U << 16U;

/// The starting points and pitch angles of the field lines every trace lane walks.
///
/// L = 2…8 Re and pitch 30…80 deg, which is the shape a real drift-shell batch has: the inner belt
/// through the outer belt, and mirror points from deep in the loss cone to near the equator. Never
/// 90 deg (`I` is identically zero there, so the trace returns immediately and the benchmark would
/// measure the early return) and never on the dipole axis.
struct TraceInputs {
    std::vector<Position<Frame::GEO>> starts;  ///< kTraceMax GEO start points, Earth radii.
    std::vector<double> pitch_deg;             ///< kTraceMax local pitch angles, degrees.
    std::vector<FieldLine> out;                ///< kTraceMax results, written by every lane.
    std::vector<Status> statuses;              ///< kTraceMax per-line statuses.
};

/// The van der Corput radical inverse in base 2 — bit-reversal, read as a fraction of 1.
///
/// This is what makes the crossover sweep honest. The obvious `i/N` sweep would give the 64-line
/// batch nothing but L = 2 shells and the 65536-line batch the whole L = 2…8 range, so the two
/// batches would be tracing DIFFERENT physics and the ratio between them would be part input-set
/// artefact. The radical inverse is stratified in every prefix: the first 64 values already cover
/// [0, 1) evenly, and so do the first 65536. Every point of the crossover curve therefore traces
/// the same distribution of field lines, and the only thing that changes along it is the batch
/// size — which is the whole claim the curve makes.
///
/// @param i the index.
/// @return the radical inverse of @p i in `[0, 1)`.
/// @complexity O(1) — one bit reversal.
/// @alloc none.
constexpr double radical_inverse_base2(std::uint32_t i) {
    i = (i << 16U) | (i >> 16U);
    i = ((i & 0x00FF00FFU) << 8U) | ((i & 0xFF00FF00U) >> 8U);
    i = ((i & 0x0F0F0F0FU) << 4U) | ((i & 0xF0F0F0F0U) >> 4U);
    i = ((i & 0x33333333U) << 2U) | ((i & 0xCCCCCCCCU) >> 2U);
    i = ((i & 0x55555555U) << 1U) | ((i & 0xAAAAAAAAU) >> 1U);
    return static_cast<double>(i) * 2.3283064365386963e-10;   // 2^-32
}

/// Build the trace input set.
///
/// Deterministic, and stratified rather than swept: see @ref radical_inverse_base2 for why any
/// prefix of this set has to look like the whole of it. Longitude advances by a golden-ratio turn
/// on top, so no two lines share a meridian and the set never degenerates into a handful of
/// repeated traces the cache would flatter.
///
/// @return the populated set; a function-local static, built on first use.
/// @complexity O(kTraceMax).
/// @alloc four vectors of kTraceMax elements, once, at start-up — never inside a timed loop.
const TraceInputs& trace_inputs() {
    static const TraceInputs* built = [] {
        auto* in = new TraceInputs{};
        in->starts.resize(kTraceMax);
        in->pitch_deg.resize(kTraceMax);
        in->out.resize(kTraceMax);
        in->statuses.resize(kTraceMax);
        constexpr double kGoldenTurn = 2.399963229728653;  // 2 pi (2 - phi), radians
        for (std::size_t i = 0; i < kTraceMax; ++i) {
            const double t = radical_inverse_base2(static_cast<std::uint32_t>(i));
            const double l_shell = 2.0 + (6.0 * t);
            const double lon = static_cast<double>(i) * kGoldenTurn;
            // Start on the dipole equator of the shell: r = L, latitude 0. The trace's first stage
            // walks to the true (IGRF) minimum-B point from there, which is the work a real caller
            // hands it.
            in->starts[i] = Position<Frame::GEO>{
                fx::vec3d{l_shell * std::cos(lon), l_shell * std::sin(lon), 0.0}};
            in->pitch_deg[i] = 30.0 + (50.0 * t);
        }
        return in;
    }();
    return *built;
}

/// `trace_invariant` — ONE field line, fp64, the reference lane.
///
/// The unit is one traced line, not one point of an ephemeris: this is the routine whose cost the
/// crossover curve is a multiple of. It is not the fast path and is not meant to be — see the
/// header's own note on why a per-point loop cannot be accelerated.
/// @param state the benchmark state.
/// @complexity O(steps) IGRF evaluations per line, ~4 per step.
/// @alloc none — the trace stores no path.
void BM_cpu_trace_invariant(benchmark::State& state) {
    const Inputs& in = inputs();
    const TraceInputs& tr = trace_inputs();
    std::size_t k = 0;
    for (auto _ : state) {
        const Position<Frame::GEO>* p = opaque(tr.starts.data());
        const double* a = opaque(tr.pitch_deg.data());
        benchmark::DoNotOptimize(trace_invariant(*in.igrf13, p[k], a[k]));
        k = (k + 1) % kTraceMax;
    }
    set_points(state, 1);
}
BENCHMARK(BM_cpu_trace_invariant)->Unit(benchmark::kMicrosecond);

/// `mcilwain_l` — Hilton's closed form, the arithmetic that turns `I` and `B_m` into `L_m`.
/// @param state the benchmark state.
/// @complexity O(kPoints) per iteration, two cube roots each.
/// @alloc none.
void BM_cpu_mcilwain_l(benchmark::State& state) {
    const Inputs& in = inputs();
    const double m = dipole_moment(*in.igrf13);
    for (auto _ : state) {
        const Position<Frame::GEO>* p = opaque(in.geo.data());
        for (std::size_t k = 0; k < kPoints; ++k) {
            // I and B_m are derived from the point so the inputs vary and none is a special case.
            const double r = fx::norm(p[k].v);
            benchmark::DoNotOptimize(mcilwain_l(0.1 * r, 3.0e4 / (r * r * r), m));
        }
    }
    set_points(state, kPoints);
}
BENCHMARK(BM_cpu_mcilwain_l);

/// `dipole_moment` — `M` from the epoch's own degree-1 coefficients.
/// @param state the benchmark state.
/// @complexity O(kPoints) per iteration.
/// @alloc none.
void BM_cpu_dipole_moment(benchmark::State& state) {
    const Inputs& in = inputs();
    for (auto _ : state) {
        const Igrf<13>* m = opaque(&*in.igrf13);
        for (std::size_t k = 0; k < kPoints; ++k) benchmark::DoNotOptimize(dipole_moment(*m));
    }
    set_points(state, kPoints);
}
BENCHMARK(BM_cpu_dipole_moment);

/// The batch sizes the crossover sweep visits, and the units both its lanes report in.
///
/// REAL time, not CPU time. A device lane spends nearly all of its wall clock blocked in
/// `vkQueueWaitIdle` and almost none of it on a core, so Google Benchmark's default CPU-time
/// accounting would credit the device with a throughput several times its actual one — the exact
/// error that makes a GPU claim worthless. The host lane reports real time too, so the ratio
/// between them is a ratio of the same quantity.
///
/// The repetition count and minimum time are deliberately NOT set here: Google Benchmark writes
/// `->MinTime` into the benchmark's NAME, and a name that moves when the harness is retuned breaks
/// every join between the manifest and the JSON. `scripts/bench_run.sh` passes them on the command
/// line instead.
/// @param b the benchmark to attach them to.
/// @complexity O(1).
/// @alloc Google Benchmark's argument vector.
void crossover_sizes(benchmark::Benchmark* b) {
    for (std::int64_t n = 64; n <= static_cast<std::int64_t>(kTraceMax); n *= 4) b->Arg(n);
    b->Unit(benchmark::kMicrosecond)->UseRealTime();
}

/// The HOST lane of the batch tracer: exactly the loop `trace_invariant_batch` falls back to.
///
/// Written out here rather than called through `trace_invariant_batch` because that entry point
/// consults the crossover and would silently take the device above 512 lines — which is the right
/// behaviour for a caller and useless for measuring the two lanes against each other.
/// @param state the benchmark state; `range(0)` is the batch size.
/// @complexity O(lines x steps) IGRF evaluations, serial.
/// @alloc none — the output vectors are sized once at start-up.
void BM_cpu_trace_batch(benchmark::State& state) {
    const Inputs& in = inputs();
    const TraceInputs& tr = trace_inputs();
    const auto n = static_cast<std::size_t>(state.range(0));
    auto& out = const_cast<std::vector<FieldLine>&>(tr.out);
    for (auto _ : state) {
        for (std::size_t i = 0; i < n; ++i) {
            out[i] = trace_invariant(*in.igrf13, tr.starts[i], tr.pitch_deg[i]).value;
        }
        benchmark::DoNotOptimize(out.data());
        benchmark::ClobberMemory();
    }
    set_points(state, n);
}
BENCHMARK(BM_cpu_trace_batch)->Apply(crossover_sizes);

/// The DEVICE lane of the batch tracer, transfers and coefficient staging INCLUDED.
///
/// `detail::trace_batch_on_device` is what `trace_invariant_batch` calls above the crossover, so
/// this measures the seam a caller actually reaches, not a bare kernel time. Everything the device
/// lane has to do to be usable — interpolating the coefficients to the epoch, staging them, the
/// dispatch, and reading 4N floats plus N status words back — is inside the timed region.
/// @param state the benchmark state; `range(0)` is the batch size.
/// @complexity One dispatch over ceil(N/256) workgroups; O(N x steps) field evaluations, concurrent.
/// @alloc the lane's own staging vectors and seven pooled device buffers, per launch.
void BM_gpu_trace_batch(benchmark::State& state) {
#ifdef CHEATAH_SPACE_IRBEM_LSTAR_GPU
    const Inputs& in = inputs();
    const TraceInputs& tr = trace_inputs();
    const auto n = static_cast<std::size_t>(state.range(0));
    auto& out = const_cast<std::vector<FieldLine>&>(tr.out);
    auto& st = const_cast<std::vector<Status>&>(tr.statuses);
    if (!gpu::available()) {
        state.SkipWithError(gpu::unavailable_reason().c_str());
        return;
    }
    const auto run = [&] {
        return detail::trace_batch_on_device(
            *in.igrf13, std::span<const Position<Frame::GEO>>(tr.starts.data(), n),
            std::span<const double>(tr.pitch_deg.data(), n), std::span<FieldLine>(out.data(), n),
            std::span<Status>(st.data(), n), TraceOptions{});
    };
    // Warm the pipeline cache and the buffer pool: the first launch of a kernel compiles it, and
    // timing that once per size would be timing the driver.
    if (const Result<bool> warm = run(); !warm.value) {
        state.SkipWithError("device lane declined the batch (no SPIR-V, or no device)");
        return;
    }
    for (auto _ : state) {
        benchmark::DoNotOptimize(run());
        benchmark::ClobberMemory();
    }
    set_points(state, n);
#else
    state.SkipWithError("built without cheatah-gpu-linalg on the include path");
#endif
}
BENCHMARK(BM_gpu_trace_batch)->Apply(crossover_sizes);

/// IRBEM's `NALPHA_MAX` — the fixed second dimension of every per-pitch-angle output array.
/// `docs/source/api/general_information.rst` fixes it at 25, and the arrays are that long whatever
/// `Nipa` is, so under-sizing one is a stack smash rather than a wrong number.
constexpr std::size_t kNalphaMax = 25;

/// IRBEM's `make_lstar1` with `options(1) = 0` — trace, `L_m`, `B_min` and `XJ`, and NO drift
/// shell, for a LOCALLY MIRRORING particle (α = 90°, see the note on
/// @ref BM_irbem_make_lstar_pitch). Kept as a variant rather than as the trace comparison, because
/// it does not integrate the arc this suite's lines integrate.
///
/// `options` is the documented five-element control array (`docs/source/api/make_lstar.rst`):
/// `(1)` whether to compute L\*, `(2)` the IGRF re-initialisation cadence, `(3)`/`(4)` the drift
/// and radial resolutions, `(5)` the internal field selector. `kext = 0` is IGRF-internal only,
/// which is the field this module implements, so the two sides are computing the same physics.
/// @param state the benchmark state.
/// @param compute_lstar the `options(1)` value: 0 for the trace alone, 1 for the drift shell too.
/// @complexity O(steps) field evaluations per call, serial.
/// @alloc none inside the loop.
void bench_irbem_make_lstar(benchmark::State& state, int compute_lstar) {
    const Oracle& o = oracle();
    if (oracle_missing(state, reinterpret_cast<const void*>(o.make_lstar))) return;
    const TraceInputs& tr = trace_inputs();
    std::size_t k = 0;
    for (auto _ : state) {
        int ntime = 1;
        int kext = 0;
        int sysaxes = 1;  // GEO Cartesian, Re
        std::array<int, 5> options{compute_lstar, 0, 0, 0, 0};
        std::array<int, 1> iyear{kRefYear};
        std::array<int, 1> idoy{kRefDoy};
        std::array<double, 1> ut{kRefUt};
        std::array<double, 1> x1{tr.starts[k].v[0]};
        std::array<double, 1> x2{tr.starts[k].v[1]};
        std::array<double, 1> x3{tr.starts[k].v[2]};
        std::array<double, 25> maginput{};
        std::array<double, 1> lm{};
        std::array<double, 1> lstar{};
        std::array<double, 1> blocal{};
        std::array<double, 1> bmin{};
        std::array<double, 1> xj{};
        std::array<double, 1> mlt{};
        o.make_lstar(&ntime, &kext, options.data(), &sysaxes, iyear.data(), idoy.data(), ut.data(),
                     x1.data(), x2.data(), x3.data(), maginput.data(), lm.data(), lstar.data(),
                     blocal.data(), bmin.data(), xj.data(), mlt.data());
        benchmark::DoNotOptimize(lm);
        benchmark::DoNotOptimize(xj);
        k = (k + 1) % kTraceMax;
    }
    set_points(state, 1);
}
BENCHMARK_CAPTURE(bench_irbem_make_lstar, trace, 0)
    ->Name("BM_irbem_make_lstar_trace")
    ->Unit(benchmark::kMicrosecond);
BENCHMARK_CAPTURE(bench_irbem_make_lstar, full, 1)
    ->Name("BM_irbem_make_lstar_full")
    ->Unit(benchmark::kMicrosecond);

/// IRBEM's `make_lstar_shell_splitting1` at THIS suite's own pitch angles — the like-for-like
/// counterpart of @ref BM_cpu_trace_invariant, and the row the trace comparison must be read from.
///
/// `make_lstar1` is not that counterpart, and the table used to say it was. IRBEM's own
/// documentation (`docs/source/api/magnetic_coordinates.rst`) states that `MAKE_LSTAR` "computes
/// the L\* parameter for locally mirroring particles (local pitch angle of 90 degrees)" and refers
/// the reader to `MAKE_LSTAR_SHELL_SPLITTING` "to compute L\* for arbitrary pitch angles". A
/// locally mirroring particle mirrors AT the spacecraft, so its `I` integral spans the short arc
/// between the point and its conjugate; this suite's lines mirror at `B_local / sin²α` for
/// α = 30…80°, which is up to four times the local field and a far longer arc. Measured at the same
/// eight start points, `XJ` comes back 1.85 Re from the shell-splitting call against 5.5e-03 Re
/// from `make_lstar1` at L = 2 — two different integrals, not two implementations of one.
///
/// Called with `Nipa = 1`. At α = 90° this entry point reproduces `make_lstar1`'s `Lm` and `XJ` to
/// every printed digit and costs about twice as much (30.8 µs against 15.7 µs on this machine), so
/// roughly a 2× constant of its own is inside this row; everything above that is the longer arc.
///
/// @param state the benchmark state.
/// @complexity O(steps) field evaluations per call, serial.
/// @alloc none inside the loop.
void BM_irbem_make_lstar_pitch(benchmark::State& state) {
    const Oracle& o = oracle();
    if (oracle_missing(state, reinterpret_cast<const void*>(o.make_lstar_ss))) return;
    const TraceInputs& tr = trace_inputs();
    std::size_t k = 0;
    for (auto _ : state) {
        int ntime = 1;
        int nipa = 1;
        int kext = 0;
        int sysaxes = 1;  // GEO Cartesian, Re
        std::array<int, 5> options{0, 0, 0, 0, 0};
        std::array<int, 1> iyear{kRefYear};
        std::array<int, 1> idoy{kRefDoy};
        std::array<double, 1> ut{kRefUt};
        std::array<double, 1> x1{tr.starts[k].v[0]};
        std::array<double, 1> x2{tr.starts[k].v[1]};
        std::array<double, 1> x3{tr.starts[k].v[2]};
        // NALPHA_MAX is 25 (docs/source/api/general_information.rst) and the outputs are
        // [ntime, NALPHA_MAX], so every per-pitch-angle array is sized 25 even at Nipa = 1.
        std::array<double, kNalphaMax> alpha{};
        alpha[0] = tr.pitch_deg[k];
        std::array<double, 25> maginput{};
        std::array<double, kNalphaMax> lm{};
        std::array<double, kNalphaMax> lstar{};
        std::array<double, kNalphaMax> bmirr{};
        std::array<double, kNalphaMax> xj{};
        std::array<double, 1> bmin{};
        std::array<double, 1> mlt{};
        o.make_lstar_ss(&ntime, &nipa, &kext, options.data(), &sysaxes, iyear.data(), idoy.data(),
                        ut.data(), x1.data(), x2.data(), x3.data(), alpha.data(), maginput.data(),
                        lm.data(), lstar.data(), bmirr.data(), bmin.data(), xj.data(), mlt.data());
        benchmark::DoNotOptimize(lm);
        benchmark::DoNotOptimize(xj);
        k = (k + 1) % kTraceMax;
    }
    set_points(state, 1);
}
BENCHMARK(BM_irbem_make_lstar_pitch)->Unit(benchmark::kMicrosecond);

// ---------------------------------------------------------------------------------------------
// igrf.hpp on the device — reached PAST the seam, which is the porting gap the table reports
// ---------------------------------------------------------------------------------------------

/// How many points the IGRF batch lanes move per iteration.
///
/// Well past the registry's measured crossover of 128, and small enough that the fp64 host lane
/// (~0.3 us/eval) still finishes an iteration in ~20 ms, so a repetition is a measurement rather
/// than a wait.
constexpr std::size_t kIgrfBatch = 1U << 16U;

/// The IGRF batch inputs, GEO, Earth radii, fp32 for the device and fp64 for the host.
struct IgrfBatch {
    std::vector<Position<Frame::GEO>> geo64;  ///< kIgrfBatch fp64 points for the host lane.
    std::vector<float> pos32;                 ///< 3*kIgrfBatch fp32 floats for the device lane.
    std::vector<float> out32;                 ///< 3*kIgrfBatch floats, written by the device lane.
    std::vector<float> coef;                  ///< g then h, interpolated to the epoch.
    std::vector<float> norm;                  ///< The Legendre normalisation, e, f, diagonal.
};

/// Build the IGRF batch inputs and stage the epoch's coefficients.
///
/// The coefficient interpolation and the Legendre normalisation happen ONCE, here, exactly as
/// `lstar.hpp`'s device lane does it: the 26-epoch table stays on the host because interpolating it
/// per thread would be 10^5 redundant copies of a calculation done once per batch. Staging it
/// outside the timed loop is therefore not flattering the device — it is what the seam does.
///
/// @return the populated set; a function-local static, built on first use.
/// @complexity O(kIgrfBatch).
/// @alloc five vectors, once, at start-up.
const IgrfBatch& igrf_batch() {
    static const IgrfBatch* built = [] {
        const Inputs& in = inputs();
        auto* b = new IgrfBatch{};
        b->geo64.resize(kIgrfBatch);
        b->pos32.resize(3 * kIgrfBatch);
        b->out32.resize(3 * kIgrfBatch);
        constexpr double kGoldenTurn = 2.399963229728653;
        for (std::size_t i = 0; i < kIgrfBatch; ++i) {
            const double t = static_cast<double>(i) / static_cast<double>(kIgrfBatch);
            const double r = 1.05 + (7.0 * t);
            const double lat = 1.0 - (2.0 * t);   // sin(latitude), never a pole
            const double lon = static_cast<double>(i) * kGoldenTurn;
            const double cl = std::sqrt(std::max(0.0, 1.0 - (lat * lat)));
            b->geo64[i] = Position<Frame::GEO>{
                fx::vec3d{r * cl * std::cos(lon), r * cl * std::sin(lon), r * lat}};
            b->pos32[(3 * i) + 0] = static_cast<float>(b->geo64[i].v[0]);
            b->pos32[(3 * i) + 1] = static_cast<float>(b->geo64[i].v[1]);
            b->pos32[(3 * i) + 2] = static_cast<float>(b->geo64[i].v[2]);
        }
        constexpr int kNmax = 13;
        constexpr int kSlots = ((kNmax + 1) * (kNmax + 2)) / 2;
        b->coef.assign(2 * kSlots, 0.0F);
        for (int deg = 1; deg <= kNmax; ++deg) {
            for (int m = 0; m <= deg; ++m) {
                const auto k = static_cast<std::size_t>((deg * (deg + 1)) / 2 + m);
                b->coef[k] = static_cast<float>(in.igrf13->g(deg, m));
                b->coef[static_cast<std::size_t>(kSlots) + k] =
                    static_cast<float>(in.igrf13->h(deg, m));
            }
        }
        constexpr auto kNorm = detail::make_legendre_normalisation<kNmax, double>();
        b->norm.assign(static_cast<std::size_t>(2 * kSlots + kNmax + 1), 0.0F);
        for (int k = 0; k < kSlots; ++k) {
            b->norm[static_cast<std::size_t>(k)] =
                static_cast<float>(kNorm.e[static_cast<std::size_t>(k)]);
            b->norm[static_cast<std::size_t>(kSlots + k)] =
                static_cast<float>(kNorm.f[static_cast<std::size_t>(k)]);
        }
        for (int deg = 0; deg <= kNmax; ++deg) {
            b->norm[static_cast<std::size_t>((2 * kSlots) + deg)] =
                static_cast<float>(kNorm.diagonal[static_cast<std::size_t>(deg)]);
        }
        return b;
    }();
    return *built;
}

/// `Igrf<13>::evaluate` over a whole batch, fp64 — the host lane of the IGRF comparison.
///
/// The model reference and the point array are hoisted out of the timed loop. They are loop
/// invariants and a benchmark should hoist its invariants anyway, but there is a sharper reason
/// here: with GCC 13.3 at `-O2` and above, dereferencing the `std::optional<Igrf<13>>` INSIDE the
/// loop segfaulted, `evaluate` receiving a null point array. Not diagnosed past that — the loop is
/// correct either way and this form is the one a benchmark should have had from the start — but it
/// is recorded because "it worked at -O1" is exactly the shape of finding that gets forgotten.
/// @param state the benchmark state.
/// @complexity O(kIgrfBatch) evaluations per iteration, each O(NMAX^2).
/// @alloc none inside the loop.
void BM_cpu_igrf_batch(benchmark::State& state) {
    const Inputs& in = inputs();
    const IgrfBatch& b = igrf_batch();
    const Igrf<13>& model = *in.igrf13;
    const Position<Frame::GEO>* const points = b.geo64.data();
    for (auto _ : state) {
        const Position<Frame::GEO>* p = opaque(points);
        for (std::size_t i = 0; i < kIgrfBatch; ++i) {
            benchmark::DoNotOptimize(model.evaluate(p[i]));
        }
    }
    set_points(state, kIgrfBatch);
}
BENCHMARK(BM_cpu_igrf_batch)->Unit(benchmark::kMillisecond)->UseRealTime();

/// `igrf_batch_host` over the same batch — the SIMD strip lane of batch_soa.hpp, i.e. what
/// `field_batch`'s host path now runs. Same points, same model, same fp64 answers to the bit
/// (IrbemBatchSimd.BatchLaneIsBitIdenticalToTheScalarLane); the ratio to `BM_cpu_igrf_batch` is
/// the whole claim of that header, and it is measured here rather than quoted. The output spans
/// are sized once at start-up; the timed region allocates nothing.
/// @param state the benchmark state.
/// @complexity O(kIgrfBatch) evaluations per iteration, each O(NMAX^2), eight per strip.
/// @alloc none inside the loop.
void BM_cpu_igrf_batch_soa(benchmark::State& state) {
    const Inputs& in = inputs();
    const IgrfBatch& b = igrf_batch();
    const Igrf<13>& model = *in.igrf13;
    static std::vector<FieldVector<Frame::GEO>> out(kIgrfBatch);
    static std::vector<double> mag(kIgrfBatch);
    for (auto _ : state) {
        const Position<Frame::GEO>* p = opaque(b.geo64.data());
        Status s = igrf_batch_host(model, std::span<const Position<Frame::GEO>>(p, kIgrfBatch),
                                   std::span<FieldVector<Frame::GEO>>(out), std::span<double>(mag));
        benchmark::DoNotOptimize(s);
        benchmark::DoNotOptimize(out.data());
        benchmark::DoNotOptimize(mag.data());
        benchmark::ClobberMemory();
    }
    set_points(state, kIgrfBatch);
}
BENCHMARK(BM_cpu_igrf_batch_soa)->Unit(benchmark::kMillisecond)->UseRealTime();

#ifdef CHEATAH_SPACE_IRBEM_LSTAR_GPU
/// Launch `irbem_igrf_f32` over the staged batch, through `dispatch.hpp`'s own `launch_igrf`.
///
/// This used to be a launcher the benchmark owned, because the seam genuinely had none: the IGRF
/// kernel binds five buffers (pos, out, coef, norm, dims) where `dispatch_batch` speaks four, and
/// `launch_trace` is the tracer's own seven-binding one. `launch_igrf` since landed in the seam, so
/// the duplicate is gone — a benchmark that reaches past the seam is measuring a program no caller
/// can run, and once the seam exists, keeping a second copy of it here would let the two drift and
/// the table would be quoting the copy.
///
/// @param b the staged batch.
/// @param n how many points to launch.
/// @return false when there is no device or no compiled SPIR-V.
/// @complexity One dispatch over ceil(n/256) workgroups, plus 6n floats over the bus.
/// @alloc five pooled device buffers per launch, released on return.
bool launch_igrf_bench(const IgrfBatch& b, std::size_t n) {
    auto& out = const_cast<std::vector<float>&>(b.out32);
    const std::array<std::uint32_t, 2> dims{static_cast<std::uint32_t>(n), 13U};
    return gpu::launch_igrf(std::span<const float>(b.pos32.data(), 3 * n),
                            std::span<const float>(b.coef), std::span<const float>(b.norm),
                            std::span<const std::uint32_t>(dims),
                            std::span<float>(out.data(), 3 * n));
}
#endif

/// `irbem_igrf_f32` over a whole batch, transfers INCLUDED.
/// @param state the benchmark state.
/// @complexity O(kIgrfBatch) device work, plus 6*kIgrfBatch floats over the bus, per iteration.
/// @alloc five pooled device buffers per launch.
void BM_gpu_igrf_batch(benchmark::State& state) {
#ifdef CHEATAH_SPACE_IRBEM_LSTAR_GPU
    const IgrfBatch& b = igrf_batch();
    if (!gpu::available()) {
        state.SkipWithError(gpu::unavailable_reason().c_str());
        return;
    }
    if (!launch_igrf_bench(b, kIgrfBatch)) {   // warm the pipeline cache and the pool
        state.SkipWithError("no compiled SPIR-V for irbem_igrf_f32");
        return;
    }
    for (auto _ : state) {
        benchmark::DoNotOptimize(launch_igrf_bench(b, kIgrfBatch));
        benchmark::ClobberMemory();
    }
    set_points(state, kIgrfBatch);
#else
    state.SkipWithError("built without cheatah-gpu-linalg on the include path");
#endif
}
BENCHMARK(BM_gpu_igrf_batch)->Unit(benchmark::kMillisecond)->UseRealTime();

// ---------------------------------------------------------------------------------------------
// The manifest — which lane exists for which routine, DERIVED, never asserted
// ---------------------------------------------------------------------------------------------

/// One row of the generated table: a routine and the benchmarks that measure it on each lane.
struct Routine {
    const char* label;        ///< The routine, as the table's first column names it.
    const char* cpu_bench;    ///< The CPU benchmark's name; never null.
    const char* gpu_bench;    ///< The device benchmark's name, or null when there is no device lane.
    const char* irbem_bench;  ///< The oracle benchmark's name, or null when IRBEM has no counterpart.
    const char* kernel;       ///< The Slang entry point that would implement it, or null.
};

/// Every routine the table reports, in the order it reports them.
constexpr std::array<Routine, 48> kRoutines{{
    {"datetime: julian_day_number", "BM_cpu_julian_day_number", nullptr, "BM_irbem_julday", nullptr},
    {"datetime: calendar_date", "BM_cpu_calendar_date", nullptr, "BM_irbem_caldat", nullptr},
    {"datetime: day_of_year", "BM_cpu_day_of_year", nullptr, "BM_irbem_get_doy", nullptr},
    {"datetime: decimal_year", "BM_cpu_decimal_year", nullptr, "BM_irbem_date_and_time2decy", nullptr},
    {"datetime: date_and_time_from_decimal_year", "BM_cpu_date_and_time_from_decimal_year", nullptr,
     "BM_irbem_decy2date_and_time", nullptr},
    {"datetime: date_and_time_from_doy_and_ut", "BM_cpu_date_and_time_from_doy_and_ut", nullptr,
     "BM_irbem_doy_and_ut2date_and_time", nullptr},

    {"rotations: Rotations::at (per epoch)", "BM_cpu_rotations_build", nullptr, nullptr, nullptr},
    {"transform: GEO->GSM (hot)", "BM_cpu_geo_to_gsm_hot", nullptr, "BM_irbem_geo2gsm", nullptr},
    {"transform: GEO->GSM (cold)", "BM_cpu_geo_to_gsm_cold", nullptr, "BM_irbem_geo2gsm", nullptr},
    {"transform: GEO->GSE (hot)", "BM_cpu_geo_to_gse_hot", nullptr, "BM_irbem_geo2gse", nullptr},
    {"transform: GEO->GSE (cold)", "BM_cpu_geo_to_gse_cold", nullptr, "BM_irbem_geo2gse", nullptr},
    {"transform: GEO->SM (hot)", "BM_cpu_geo_to_sm_hot", nullptr, "BM_irbem_geo2sm", nullptr},
    {"transform: GEO->SM (cold)", "BM_cpu_geo_to_sm_cold", nullptr, "BM_irbem_geo2sm", nullptr},
    {"transform: GEO->MAG (hot)", "BM_cpu_geo_to_mag_hot", nullptr, "BM_irbem_geo2mag", nullptr},
    {"transform: GEO->MAG (cold)", "BM_cpu_geo_to_mag_cold", nullptr, "BM_irbem_geo2mag", nullptr},
    {"transform: GEO->GEI (hot)", "BM_cpu_geo_to_gei_hot", nullptr, "BM_irbem_geo2gei", nullptr},
    {"transform: GEO->GEI (cold)", "BM_cpu_geo_to_gei_cold", nullptr, "BM_irbem_geo2gei", nullptr},
    {"transform: GSM->SM (hot)", "BM_cpu_gsm_to_sm_hot", nullptr, "BM_irbem_gsm2sm", nullptr},
    {"transform: GSM->SM (cold)", "BM_cpu_gsm_to_sm_cold", nullptr, "BM_irbem_gsm2sm", nullptr},
    {"transform: GSM->GEO, transposed (hot)", "BM_cpu_gsm_to_geo_hot", nullptr, "BM_irbem_gsm2geo",
     nullptr},
    {"transform: GSM->GEO, transposed (cold)", "BM_cpu_gsm_to_geo_cold", nullptr,
     "BM_irbem_gsm2geo", nullptr},

    {"geodetic: gdz_to_geo", "BM_cpu_gdz_to_geo", nullptr, "BM_irbem_gdz_geo", nullptr},
    {"geodetic: geo_to_gdz (Bowring)", "BM_cpu_geo_to_gdz", nullptr, "BM_irbem_geo_gdz", nullptr},
    {"geodetic: sph_to_car", "BM_cpu_sph_to_car", nullptr, "BM_irbem_sph_car", nullptr},
    {"geodetic: car_to_sph", "BM_cpu_car_to_sph", nullptr, "BM_irbem_car_sph", nullptr},
    {"geodetic: rll_to_gdz", "BM_cpu_rll_to_gdz", nullptr, "BM_irbem_rll_gdz", nullptr},

    {"helio: helio_geometry (per epoch)", "BM_cpu_helio_geometry", nullptr, nullptr, nullptr},
    {"helio: HAE->HEE (hot)", "BM_cpu_hae_to_hee_hot", nullptr, "BM_irbem_hae2hee", nullptr},
    {"helio: HAE->HEE (cold)", "BM_cpu_hae_to_hee_cold", nullptr, "BM_irbem_hae2hee", nullptr},
    {"helio: HEE->HAE (hot)", "BM_cpu_hee_to_hae_hot", nullptr, "BM_irbem_hee2hae", nullptr},
    {"helio: HEE->HAE (cold)", "BM_cpu_hee_to_hae_cold", nullptr, "BM_irbem_hee2hae", nullptr},
    {"helio: HAE->HEEQ (hot)", "BM_cpu_hae_to_heeq_hot", nullptr, "BM_irbem_hae2heeq", nullptr},
    {"helio: HAE->HEEQ (cold)", "BM_cpu_hae_to_heeq_cold", nullptr, "BM_irbem_hae2heeq", nullptr},
    {"helio: HEEQ->HAE (hot)", "BM_cpu_heeq_to_hae_hot", nullptr, "BM_irbem_heeq2hae", nullptr},
    {"helio: HEEQ->HAE (cold)", "BM_cpu_heeq_to_hae_cold", nullptr, "BM_irbem_heeq2hae", nullptr},
    {"helio: GSE->HEE, position (hot)", "BM_cpu_gse_to_hee_hot", nullptr, "BM_irbem_gse2hee",
     nullptr},
    {"helio: HEE->GSE, position (hot)", "BM_cpu_hee_to_gse_hot", nullptr, "BM_irbem_hee2gse",
     nullptr},

    {"igrf: Igrf<13>::at (per epoch)", "BM_cpu_igrf_at_deg13", nullptr, nullptr, nullptr},
    {"igrf: evaluate, degree 13, GEO", "BM_cpu_igrf_evaluate_deg13", nullptr,
     "BM_irbem_get_field_igrf", "irbem_igrf_f32"},
    {"igrf: evaluate, batch of 65536", "BM_cpu_igrf_batch/real_time",
     "BM_gpu_igrf_batch/real_time",
     "BM_irbem_get_field_igrf", "irbem_igrf_f32"},
    {"igrf: batch of 65536, SIMD strips (field_batch host lane)",
     "BM_cpu_igrf_batch_soa/real_time", nullptr, "BM_irbem_get_field_igrf", nullptr},
    {"dipole: batch, fp32", "BM_cpu_dipole_field_host_batch/real_time", "BM_gpu_dipole_field_batch/real_time",
     "BM_irbem_get_field_dipole", "irbem_dipole_f32"},

    {"lstar: trace_invariant (one line)", "BM_cpu_trace_invariant", nullptr,
     "BM_irbem_make_lstar_pitch", "irbem_trace_i_f32"},
    {"lstar: trace batch of 65536", "BM_cpu_trace_batch/65536/real_time",
     "BM_gpu_trace_batch/65536/real_time", "BM_irbem_make_lstar_pitch", "irbem_trace_i_f32"},
    {"T89: t89_field, Kp = 0 (quiet)", "BM_cpu_t89_field/0", nullptr,
     "BM_irbem_get_field_t89/0", "irbem_t89_f32"},
    {"T89: t89_field, Kp = 9- (extreme storm)", "BM_cpu_t89_field/90", nullptr,
     "BM_irbem_get_field_t89/90", "irbem_t89_f32"},

    {"lstar: mcilwain_l (Hilton)", "BM_cpu_mcilwain_l", nullptr, nullptr, nullptr},
    {"lstar: dipole_moment", "BM_cpu_dipole_moment", nullptr, nullptr, nullptr},
}};

/// The extra routines that exist but do not head a table row of their own, kept here so the
/// manifest is a complete listing of what the binary measures rather than a filtered one.
constexpr std::array<Routine, 14> kExtraRoutines{{
    {"igrf: evaluate, degree 10, GEO", "BM_cpu_igrf_evaluate_deg10", nullptr, nullptr,
     "irbem_igrf_f32"},
    {"igrf: evaluate, degree 13, spherical", "BM_cpu_igrf_evaluate_deg13_sph", nullptr, nullptr,
     "irbem_igrf_f32"},
    {"dipole: scalar, fp64", "BM_cpu_dipole_field_at", nullptr, "BM_irbem_get_field_dipole",
     "irbem_dipole_f32"},
    {"helio: GSE->HEE, position (cold)", "BM_cpu_gse_to_hee_cold", nullptr, "BM_irbem_gse2hee",
     nullptr},
    {"T89: t89_field, Kp = 1.0", "BM_cpu_t89_field/10", nullptr,
     "BM_irbem_get_field_t89/10", "irbem_t89_f32"},
    {"T89: t89_field, Kp = 2.0", "BM_cpu_t89_field/20", nullptr,
     "BM_irbem_get_field_t89/20", "irbem_t89_f32"},
    {"T89: t89_field, Kp = 3.0", "BM_cpu_t89_field/30", nullptr,
     "BM_irbem_get_field_t89/30", "irbem_t89_f32"},
    {"T89: t89_field, Kp = 4.0", "BM_cpu_t89_field/40", nullptr,
     "BM_irbem_get_field_t89/40", "irbem_t89_f32"},
    {"T89: t89_field, Kp = 5.0", "BM_cpu_t89_field/50", nullptr,
     "BM_irbem_get_field_t89/50", "irbem_t89_f32"},
    {"T89: t89_field, Kp = 6.0", "BM_cpu_t89_field/60", nullptr,
     "BM_irbem_get_field_t89/60", "irbem_t89_f32"},
    {"T89: t89_field, Kp = 7.0", "BM_cpu_t89_field/70", nullptr,
     "BM_irbem_get_field_t89/70", "irbem_t89_f32"},
    {"T89: t89_field, Kp = 8.0", "BM_cpu_t89_field/80", nullptr,
     "BM_irbem_get_field_t89/80", "irbem_t89_f32"},
    {"lstar: trace at alpha = 90 deg only (locally mirroring)", "-", nullptr,
     "BM_irbem_make_lstar_trace", nullptr},
    // `driftshell.hpp::make_lstar` landed while this suite was being written, so the row no longer
    // says NOT PORTED — but nothing here times it yet, and an empty CPU cell beside IRBEM's cost is
    // the honest rendering of "ported, not yet measured". Owed: a CPU (and device) lane for it.
    {"lstar: L* drift shell (ported; NOT YET BENCHMARKED here)", "-", nullptr,
     "BM_irbem_make_lstar_full", nullptr},
}};

/// Whether `irbem.slang` defines @p entry as a compute entry point.
///
/// Read from `gpu::shader_source_path()` — the same file the registry-completeness test parses —
/// so the answer is the source's, not a second list that could drift from it.
/// @param entry the Slang entry-point name.
/// @return whether the source declares `void <entry>(`.
/// @complexity O(size of irbem.slang), once per query.
/// @alloc the file contents.
bool slang_defines(std::string_view entry) {
    std::ifstream f(gpu::shader_source_path());
    if (!f) return false;
    std::stringstream buffer;
    buffer << f.rdbuf();
    return buffer.str().find("void " + std::string(entry) + "(") != std::string::npos;
}

/// Where `dispatch.hpp` lives, derived from the header that declares the kernels rather than from
/// the working directory — the same trick @ref gpu::shader_source_path plays, for the same reason.
/// @return the path to `gpu/dispatch.hpp`.
/// @complexity O(1).
/// @alloc two short-lived paths, one of them returned.
std::filesystem::path dispatch_header_path() {
    return gpu::shader_source_path().parent_path() / "dispatch.hpp";
}

/// The named launcher in `dispatch.hpp` that dispatches @p entry, if one exists.
///
/// This used to be a hand-written rule — "four bindings means `dispatch_batch`, and
/// `irbem_trace_i_f32` means `launch_trace`, and everything else is unreachable" — which is not a
/// derivation at all: it is a second list, in a second file, that goes stale the moment somebody
/// writes a launcher. It did. `launch_igrf` was added to the seam and this column went on
/// publishing "registered, NO LAUNCHER" about a kernel the seam had learned to launch. So the
/// answer is now READ from the seam: find `dispatch_1d("<entry>"` and name the enclosing
/// `launch_*` function.
///
/// @param entry the Slang entry-point name.
/// @return the launcher's identifier, or an empty string when no launcher names @p entry.
/// @complexity O(size of dispatch.hpp), once per query.
/// @alloc the file contents and the returned name.
std::string launcher_for(std::string_view entry) {
    // Every header of the module, not just dispatch.hpp. A launcher does not have to live beside
    // the registry and two do not: trace_path_on_device is in trace_api.hpp and the igrf batch
    // lane's is in field.hpp. Reading only dispatch.hpp reported "registered, NO LAUNCHER" about a
    // kernel a caller can reach, which is a worse error than the one this column exists to prevent
    // — it understates the module rather than overstating it, but it is equally untrue.
    const std::filesystem::path dir = gpu::shader_source_path().parent_path().parent_path();
    std::vector<std::filesystem::path> headers;
    for (const std::filesystem::directory_entry& e : std::filesystem::directory_iterator(dir)) {
        if (e.path().extension() == ".hpp") headers.push_back(e.path());
    }
    headers.push_back(dispatch_header_path());
    std::sort(headers.begin(), headers.end());

    const std::string quoted = "\"" + std::string(entry) + "\"";
    for (const std::filesystem::path& h : headers) {
        std::ifstream f(h);
        if (!f) continue;
        std::stringstream buffer;
        buffer << f.rdbuf();
        const std::string src = buffer.str();

        // The kernel's name in quotes, close after a `dispatch_1d(`. Matched in two steps rather
        // than as one literal because the call spells the name two ways: bare in the older form,
        // and wrapped as `qualified("name").c_str()` since kernels started being addressed by
        // their directory. A single-literal needle silently stopped matching when that changed.
        std::size_t at = std::string::npos;
        for (std::size_t i = src.find(quoted); i != std::string::npos; i = src.find(quoted, i + 1)) {
            const std::size_t call = src.rfind("dispatch_1d(", i);
            if (call == std::string::npos || (i - call) > 48) continue;
            at = i;
            break;
        }
        if (at == std::string::npos) continue;

        // The enclosing function, read off the last signature line that starts at column 0 —
        // which is how every launcher in this module is written. Its name is the identifier
        // immediately before that line's first `(`, and is what a caller would write.
        std::size_t line = src.rfind("\n", at);
        while (line != std::string::npos && line > 0) {
            const std::size_t prev = src.rfind("\n", line - 1);
            const std::size_t bol = prev == std::string::npos ? 0 : prev + 1;
            const std::size_t paren = src.find('(', bol);
            const bool at_margin = bol < src.size() && src[bol] != ' ' && src[bol] != '\t' &&
                                   src[bol] != '\n';
            if (at_margin && paren != std::string::npos && paren < line) {
                std::size_t e = paren;
                while (e > bol && (std::isalnum(static_cast<unsigned char>(src[e - 1])) != 0 ||
                                   src[e - 1] == '_')) {
                    --e;
                }
                if (e < paren) return src.substr(e, paren - e);
            }
            line = prev;
        }
    }
    return {};
}

/// What the "on GPU?" column says for a routine, DERIVED from `dispatch.hpp`'s registry, from
/// `dispatch.hpp`'s own launchers, and from `irbem.slang` itself.
///
/// The four answers are the four real states, and the third one is the whole reason this column
/// exists: a kernel can be written, compiled and even measured in isolation while remaining
/// unreachable through the seam every caller goes through, and that is a porting gap, not a
/// success. `dispatch_batch` binds exactly four buffers — pos, out, params, dims — so a four-slot
/// registry row is reachable through the generic seam; anything wider needs a launcher of its own,
/// and @ref launcher_for goes and looks for one instead of assuming.
/// @param r the routine.
/// @return the column text.
/// @complexity O(registry size) plus one read of `dispatch.hpp`, and of `irbem.slang` when the
///             kernel is unregistered.
/// @alloc the returned string.
std::string gpu_status(const Routine& r) {
    if (r.kernel == nullptr) return "no - host only";
    for (const gpu::KernelInfo& k : gpu::registered_kernels) {
        if (std::string_view(k.name) != std::string_view(r.kernel)) continue;
        // Four bindings is what `dispatch_batch` can express generically; a wider kernel is
        // reachable only if somebody wrote it a launcher, so go and read whether somebody did.
        if (k.bindings == 4) return std::string("yes - ") + r.kernel + " via dispatch_batch";
        if (const std::string via = launcher_for(k.name); !via.empty()) {
            return std::string("yes - ") + r.kernel + " via " + via;
        }
        return std::string("registered, NO LAUNCHER - ") + r.kernel + " binds " +
               std::to_string(k.bindings) +
               " buffers; dispatch_batch expresses 4 and no launcher in dispatch.hpp names it";
    }
    if (slang_defines(r.kernel)) {
        return std::string("kernel only - ") + r.kernel +
               " compiles in irbem.slang but dispatch.hpp registers no row for it";
    }
    return "no - host only";
}

// ---------------------------------------------------------------------------------------------
// Arithmetic intensity — the one number that predicts every verdict in the table above
// ---------------------------------------------------------------------------------------------

/// Flops charged to ONE spherical-harmonic coefficient slot of the field sum.
///
/// The term model, stated so nobody has to guess what was counted: per `(n, m)` slot the kernel
/// forms `G = g·cos mφ + h·sin mφ` (two multiplies and an add), multiplies it by the associated
/// Legendre value (one), and accumulates it (one). Five. This is a TERM count, the way the cost of
/// a harmonic model is conventionally quoted — not an instruction count of the Slang inner loop,
/// which also carries the Legendre recursion and the two horizontal-component accumulations and
/// comes out roughly six times higher.
///
/// The convention is stated rather than hidden because the argument this table makes is about the
/// RATIO between rows — dipole against IGRF against trace — and a ratio is unchanged by a constant
/// factor applied to every row. Pick the other convention and every intensity below scales by ~6
/// and every verdict is identical.
constexpr double kFlopsPerHarmonicSlot = 5.0;

/// Coefficient slots in a degree-13 expansion: the `(n, m)` pairs with `0 ≤ m ≤ n ≤ 13`, which is
/// `(N+1)(N+2)/2 = 105`. That is the triangle the kernel actually loops over. It carries
/// `N(N+2) = 195` real coefficients — each slot holds a `g` and, for `m > 0`, an `h`, and the term
/// model charges both together — and one of its slots, `n = 0`, is not part of the expansion at
/// all, so this is a ~1% overcount of the loop trip count and is deliberately not corrected: the
/// table's argument is the ratio between rows, and the same constant sits in the trace row too.
constexpr double kIgrfSlots = 105.0;

/// Flops in one IGRF-14 evaluation at degree 13, under @ref kFlopsPerHarmonicSlot.
constexpr double kFlopsPerIgrfEval = kIgrfSlots * kFlopsPerHarmonicSlot;

/// Flops in one centred-dipole evaluation — `dispatch.hpp`'s own count for `irbem_dipole_f32`,
/// which is the same twelve operations `dipole_field_at` performs.
constexpr double kFlopsPerDipoleEval = 12.0;

/// IGRF evaluations one RK4 step costs. Four stages; `rk4_step` carries the fifth — the field at
/// the arrival point — forward into the next step's `k1` rather than re-evaluating it, which is
/// where the reference implementation's fifth evaluation went.
constexpr double kIgrfEvalsPerRk4Step = 4.0;

/// One row of the intensity study.
struct Intensity {
    const char* kernel;      ///< The Slang entry point.
    const char* label;       ///< How the table names it.
    double bytes_per_point;  ///< Compulsory DRAM traffic per point, in and out, bytes.
    double flops_per_point;  ///< Flops per point under the term model; 0 when measured at run time.
    const char* host_bench;  ///< The host lane's benchmark name.
    const char* gpu_bench;   ///< The device lane's benchmark name.
};

/// Measure the mean RK4 step count of a trace over the benchmark's own input set.
///
/// The trace's intensity is the only one that is not a property of the source alone: it is
/// `4 × steps × flops-per-evaluation`, and `steps` depends on the shell, the pitch angle and the
/// field. So it is MEASURED — over a stride through the same input set the crossover sweep uses —
/// rather than assumed, which is also why this row's flops/byte is the only one that can move when
/// the input set changes.
///
/// @param sample how many lines to trace.
/// @return the mean step count over the closed lines, or 0 when none closed.
/// @complexity O(sample × steps) IGRF evaluations.
/// @alloc none.
double mean_trace_steps(std::size_t sample) {
    const Inputs& in = inputs();
    const TraceInputs& tr = trace_inputs();
    const std::size_t stride = kTraceMax / sample;
    double total = 0.0;
    std::size_t counted = 0;
    for (std::size_t i = 0; i < kTraceMax; i += stride) {
        const Result<FieldLine> r = trace_invariant(*in.igrf13, tr.starts[i], tr.pitch_deg[i]);
        if (r.status != Status::Ok) continue;
        total += r.value.steps;
        ++counted;
    }
    return counted == 0 ? 0.0 : total / static_cast<double>(counted);
}

/// Print the intensity study, one tab-separated row per kernel, and exit.
///
/// Byte counts are exact — they are the kernel's own buffer layout. Flop counts are the term model
/// above, except the tracer's, which is derived from a MEASURED mean step count. Nothing here is a
/// throughput: `scripts/bench_run.sh` joins these rows against the JSON the same binary emits, so
/// the verdict column in the published table is a measurement and the intensity column is a count.
/// @complexity O(sample × steps) for the one measured row.
/// @alloc none beyond the input sets.
void print_intensity() {
    const double steps = mean_trace_steps(256);
    const std::array<Intensity, 3> rows{{
        // 3 floats in, 3 floats out.
        {"irbem_dipole_f32", "centred dipole", 24.0, kFlopsPerDipoleEval,
         "BM_cpu_dipole_field_host_batch/real_time", "BM_gpu_dipole_field_batch/real_time"},
        // 3 floats in, 3 floats out. The coefficient and normalisation tables are 1736 bytes for
        // the WHOLE batch — 0.03 bytes/point at 65536 — so they do not appear here.
        {"irbem_igrf_f32", "IGRF-14, degree 13", 24.0, kFlopsPerIgrfEval,
         "BM_cpu_igrf_batch/real_time", "BM_gpu_igrf_batch/real_time"},
        // 3 position floats + 1 pitch float in = 16; 4 result floats + 1 status word out = 20.
        {"irbem_trace_i_f32", "trace + second invariant", 36.0,
         kIgrfEvalsPerRk4Step * steps * kFlopsPerIgrfEval, "BM_cpu_trace_batch/65536/real_time",
         "BM_gpu_trace_batch/65536/real_time"},
    }};
    std::printf("# kernel\tlabel\tbytes_per_point\tflops_per_point\tflops_per_byte\thost_bench\t"
                "gpu_bench\n");
    for (const Intensity& r : rows) {
        std::printf("%s\t%s\t%.0f\t%.0f\t%.3f\t%s\t%s\n", r.kernel, r.label, r.bytes_per_point,
                    r.flops_per_point, r.flops_per_point / r.bytes_per_point, r.host_bench,
                    r.gpu_bench);
    }
    std::printf("# mean RK4 steps per traced line, measured over 256 lines of the input set: %.1f\n",
                steps);
}

/// Print the manifest, one tab-separated row per routine, and exit.
///
/// Consumed by `scripts/bench_run.sh`, which joins it against the JSON this same binary emits.
/// Nothing in the output is a measurement — every number in the table comes from the JSON — and
/// nothing in it is hand-maintained prose about the device: @ref gpu_status derives that.
/// @complexity O(routines), plus one read of `irbem.slang`.
/// @alloc the status strings.
void print_manifest() {
    std::printf("# label\tcpu_bench\tgpu_bench\tirbem_bench\tgpu_status\tprimary\n");
    auto emit = [](const Routine& r, const char* primary) {
        std::printf("%s\t%s\t%s\t%s\t%s\t%s\n", r.label, r.cpu_bench,
                    r.gpu_bench != nullptr ? r.gpu_bench : "-",
                    r.irbem_bench != nullptr ? r.irbem_bench : "-", gpu_status(r).c_str(), primary);
    };
    for (const Routine& r : kRoutines) emit(r, "yes");
    for (const Routine& r : kExtraRoutines) emit(r, "no");

    // Every registered kernel that NO row claims, derived from the registry rather than listed.
    //
    // The `kernel` field is the one hand-maintained link between a routine and the device, and it
    // is therefore the one place this table can quietly under-report itself. It did: both T89 rows
    // carried a null kernel, so gpu_status short-circuited to "no - host only" about a kernel the
    // registry has always held with four bindings, and the published page said T89 could not reach
    // the device. That is the same failure the "on GPU?" column was rewritten to prevent — a second
    // list, in a second place, going stale — reappearing one field lower down.
    //
    // A missing row cannot be derived away (a routine's kernel is genuinely per-row knowledge), but
    // an UNCLAIMED kernel can be, so it is. These lines are the porting gap stated where the table
    // is generated, and they are what a reader should compare against the models the module says it
    // implements.
    for (const gpu::KernelInfo& k : gpu::registered_kernels) {
        bool claimed = false;
        for (const Routine& r : kRoutines) {
            claimed = claimed || (r.kernel != nullptr && std::string_view(r.kernel) == k.name);
        }
        for (const Routine& r : kExtraRoutines) {
            claimed = claimed || (r.kernel != nullptr && std::string_view(r.kernel) == k.name);
        }
        if (claimed) continue;
        const Routine probe{k.name, "-", nullptr, nullptr, k.name};
        std::printf("# unclaimed\t%s\t%s\n", k.name, gpu_status(probe).c_str());
    }
}

// ---------------------------------------------------------------------------------------------
// --verify — a fast lane that computes the WRONG answer is not a fast lane
// ---------------------------------------------------------------------------------------------

/// One lane comparison: what was compared, over how many points, and how far apart the lanes got.
struct Verification {
    const char* lane;      ///< The device lane's name.
    const char* quantity;  ///< The quantity compared.
    std::size_t n;         ///< How many values were compared.
    double max_rel;        ///< The largest relative deviation seen.
    double budget;         ///< The budget from `docs/ERROR_BUDGET.md` §5.
};

/// The relative deviation of @p dev from @p host, guarded against a zero reference.
/// @param host the fp64 reference value. @param dev the fp32 device value.
/// @return `|dev − host| / max(|host|, tiny)`.
/// @complexity O(1).
/// @alloc none.
double rel_dev(double host, double dev) {
    const double d = std::abs(dev - host);
    const double s = std::abs(host);
    return s > 1e-30 ? d / s : d;
}

/// Check every device lane the table quotes against its own host lane, and report the spread.
///
/// The suite measured three device lanes for a long time without ever asking whether any of them
/// computed the right answer. A benchmark cannot: `launch_igrf` returning `true` says a dispatch
/// was submitted, not that the coefficient buffer was packed the way the kernel unpacks it, and a
/// kernel that reads the normalisation table at the wrong stride is exactly as fast as one that
/// reads it correctly. So the speedups in `BENCHMARKS.md` were, strictly, throughput claims about
/// an unverified computation. This mode closes that: it is the reason a published ratio may be
/// quoted at all, and `bench_run.sh check` runs it before it grades a single row.
///
/// Budgets are `docs/ERROR_BUDGET.md` §5 — `Blocal` 1e-6 relative, `Bmin` 1e-5, `XJ` 1e-4. The
/// device lanes are fp32 and the host lanes fp64, so these are fp32-limited agreements between two
/// implementations, not accuracies against the physics.
///
/// @return 0 when every lane is inside its budget, 1 otherwise; 0 with a note when there is no
///         device, because "no GPU here" is not a failure of the code.
/// @complexity O(kBatch + kIgrfBatch + kVerifyLines × steps) on each of two lanes.
/// @alloc the comparison buffers, once.
int verify_device_lanes() {
#ifndef CHEATAH_SPACE_IRBEM_LSTAR_GPU
    std::printf("# no device lane in this build (no cheatah-gpu-linalg on the include path)\n");
    return 0;
#else
    if (!gpu::available()) {
        std::printf("# no device: %s\n", gpu::unavailable_reason().c_str());
        return 0;
    }
    std::vector<Verification> rows;

    // --- irbem_dipole_f32: the fp32 host twin against the kernel, component by component.
    {
        const Inputs& in = inputs();
        const auto g10 = static_cast<float>(in.dipole.g10);
        std::vector<float> host(3 * kBatch);
        std::vector<float> dev(3 * kBatch);
        gpu::dipole_field_host(std::span<const float>(in.batch_pos), std::span<float>(host), g10);
        gpu::dipole_field_gpu(std::span<const float>(in.batch_pos), std::span<float>(dev), g10);
        double worst = 0.0;
        for (std::size_t i = 0; i < kBatch; ++i) {
            // Compare the magnitude, not a component: a component near a zero crossing has no
            // relative scale of its own and would report a meaningless ratio.
            const auto mag = [](const std::vector<float>& v, std::size_t k) {
                return std::sqrt((double(v[3 * k]) * v[3 * k]) +
                                 (double(v[(3 * k) + 1]) * v[(3 * k) + 1]) +
                                 (double(v[(3 * k) + 2]) * v[(3 * k) + 2]));
            };
            worst = std::max(worst, rel_dev(mag(host, i), mag(dev, i)));
        }
        rows.push_back({"irbem_dipole_f32", "|B| vs the fp32 host twin", kBatch, worst, 1e-6});
    }

    // --- irbem_igrf_f32: the kernel against igrf.hpp's own fp64 evaluate.
    {
        const Inputs& in = inputs();
        const IgrfBatch& b = igrf_batch();
        if (launch_igrf_bench(b, kIgrfBatch)) {
            double worst = 0.0;
            for (std::size_t i = 0; i < kIgrfBatch; ++i) {
                const FieldVector<Frame::GEO> h = in.igrf13->evaluate(b.geo64[i]);
                const double hm = std::sqrt((h.v[0] * h.v[0]) + (h.v[1] * h.v[1]) +
                                            (h.v[2] * h.v[2]));
                const double dm =
                    std::sqrt((double(b.out32[3 * i]) * b.out32[3 * i]) +
                              (double(b.out32[(3 * i) + 1]) * b.out32[(3 * i) + 1]) +
                              (double(b.out32[(3 * i) + 2]) * b.out32[(3 * i) + 2]));
                worst = std::max(worst, rel_dev(hm, dm));
            }
            rows.push_back({"irbem_igrf_f32", "|B| vs the fp64 host lane", kIgrfBatch, worst, 1e-6});
        }
    }

    // --- irbem_trace_i_f32: the whole RK4 chain, compared on the invariant it exists to produce.
    {
        constexpr std::size_t kVerifyLines = 4096;
        const Inputs& in = inputs();
        const TraceInputs& tr = trace_inputs();
        std::vector<FieldLine> dev(kVerifyLines);
        std::vector<Status> st(kVerifyLines);
        const Result<bool> ran = detail::trace_batch_on_device(
            *in.igrf13, std::span<const Position<Frame::GEO>>(tr.starts.data(), kVerifyLines),
            std::span<const double>(tr.pitch_deg.data(), kVerifyLines),
            std::span<FieldLine>(dev.data(), kVerifyLines),
            std::span<Status>(st.data(), kVerifyLines), TraceOptions{});
        if (ran.value) {
            double worst_i = 0.0;
            double worst_b = 0.0;
            std::size_t compared = 0;
            for (std::size_t i = 0; i < kVerifyLines; ++i) {
                if (st[i] != Status::Ok) continue;
                const Result<FieldLine> h =
                    trace_invariant(*in.igrf13, tr.starts[i], tr.pitch_deg[i]);
                if (h.status != Status::Ok) continue;
                worst_i = std::max(worst_i, rel_dev(h.value.invariant_i, dev[i].invariant_i));
                worst_b = std::max(worst_b, rel_dev(h.value.b_min, dev[i].b_min));
                ++compared;
            }
            rows.push_back({"irbem_trace_i_f32", "I vs the fp64 host lane", compared, worst_i,
                            1e-4});
            // `Bmin` is a root-find on top of B, so §5 gives it 1e-5 rather than `Blocal`'s
            // 1e-6. The measured deviation comes in an order better than that; the gate still
            // quotes the documented budget rather than a tighter one invented here.
            rows.push_back({"irbem_trace_i_f32", "Bmin vs the fp64 host lane", compared, worst_b,
                            1e-5});
        }
    }

    std::printf("# lane\tquantity\tn\tmax_rel_dev\tbudget\tverdict\n");
    int bad = 0;
    for (const Verification& v : rows) {
        const bool ok = v.max_rel <= v.budget;
        if (!ok) ++bad;
        std::printf("%s\t%s\t%zu\t%.3g\t%.3g\t%s\n", v.lane, v.quantity, v.n, v.max_rel, v.budget,
                    ok ? "inside budget" : "OUTSIDE BUDGET");
    }
    if (rows.empty()) {
        std::printf("# no device lane ran — nothing verified\n");
        return 1;
    }
    return bad == 0 ? 0 : 1;
#endif
}

}  // namespace

/**
 * The suite's entry point.
 *
 * `--manifest` prints @ref print_manifest and exits, `--intensity` prints @ref print_intensity and
 * exits, `--verify` runs @ref verify_device_lanes and exits with its status; everything else is
 * Google Benchmark's own argument handling. All three are modes of this binary rather than separate
 * ones, so the routine list, the intensity model, the correctness check and the code that measures
 * them can never be built from different sources.
 *
 * @param argc the argument count. @param argv the arguments.
 * @return 0 on success, 1 when Google Benchmark rejects its arguments.
 * @complexity the suite's.
 * @alloc Google Benchmark's, plus the input set.
 */
int main(int argc, char** argv) {
    for (int i = 1; i < argc; ++i) {
        if (std::strcmp(argv[i], "--manifest") == 0) {
            print_manifest();
            return 0;
        }
        if (std::strcmp(argv[i], "--intensity") == 0) {
            print_intensity();
            return 0;
        }
        if (std::strcmp(argv[i], "--verify") == 0) return verify_device_lanes();
    }
    benchmark::Initialize(&argc, argv);
    if (benchmark::ReportUnrecognizedArguments(argc, argv)) return 1;
    benchmark::RunSpecifiedBenchmarks();
    benchmark::Shutdown();
    return 0;
}
