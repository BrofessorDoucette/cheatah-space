// Tests for space/irbem/trace_api.hpp — the five path-returning tracing and finder routines.
//
// Two halves, and the split matters.
//
// The FIRST half is self-contained: exact-arithmetic checks, structural invariants, and the
// degenerate cases (a start on the magnetic equator, a particle in the loss cone, a buffer too
// small). Those run everywhere and are the regression net.
//
// The SECOND half is DIFFERENTIAL against the compiled IRBEM oracle, loaded with `dlopen` as a
// black box at `/tmp/irbem-builds/libirbem-O2.so` and never linked. It is skipped when the oracle
// is absent, which is every machine but a developer's — the library is LGPL-3.0 and this repo is
// MIT, so nothing here may become a link-time dependency. Override the path with
// `CHEATAH_SPACE_IRBEM_ORACLE`.
//
// THE MODEL IS MATCHED, WHICH IS WHY THE DIFFERENTIAL NUMBERS MEAN ANYTHING. Two facts about the
// oracle were measured, not assumed, and both were surprises:
//
//   1. IRBEM's internal field is evaluated at decimal year `iyear + 0.5` REGARDLESS of `idoy` and
//      `UT` — it initialises IGRF once per year, at mid-year. Scanning decimal year against its
//      `get_field1_` output puts the minimum at exactly `y + 0.5` for every year tried.
//   2. It truncates IGRF at degree 10, as `igrf.hpp`'s brief says.
//
// With `Igrf<10>::at(year + 0.5)` our `|B|` agrees with the oracle's `Blocal` to **1.3e-13
// relative**, i.e. the two libraries evaluate the same field to round-off. Every deviation the
// differential tests below report is therefore an ALGORITHM difference — the tracer, the
// quadrature, the refinements — and not a coefficient difference dressed up as one. Getting the
// epoch wrong (2015 + doy/365) puts a 3.6e-6 model error underneath every trace and makes all of
// these numbers meaningless.

#include "space/irbem/trace_api.hpp"

#include <gtest/gtest.h>

#include <dlfcn.h>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <string>
#include <vector>

namespace {

using cheatah::fixarray::vec3d;
using namespace cheatah::space::irbem;
namespace fixarray = cheatah::fixarray;

/// The epoch every test runs at, chosen so the oracle's own IGRF initialisation lands on it
/// exactly. See the file header.
constexpr int kYear = 2015;
constexpr double kDecimalYear = kYear + 0.5;
/// IRBEM's `idoy`/`UT`. They do NOT affect its internal field (finding 1 above); they are passed
/// so the calls are well formed and so a future external field would see a real timestamp.
constexpr int kDoy = 180;
constexpr double kUt = 43200.0;

/// The model at the matched epoch, truncated where the oracle truncates.
Igrf<10> model() {
    const auto m = Igrf<10>::at(kDecimalYear);
    EXPECT_TRUE(m.has_value());
    return *m;
}

Position<Frame::GEO> geo(double x, double y, double z) { return Position<Frame::GEO>{vec3d{x, y, z}}; }

double separation(const Position<Frame::GEO>& a, const Position<Frame::GEO>& b) {
    return fixarray::norm(a.v - b.v);
}

/// A root predicate for @ref detail::refine_in_step, spelled as a STRUCT rather than a lambda so
/// every test below drives the SAME template instantiation. Two lambdas are two closure types and
/// therefore two instantiations, each covering only the branches its own call happened to take;
/// one type means the converging case and the degenerate case together exercise the whole
/// function.
struct RadiusMinus {
    double target;
    double operator()(const Position<Frame::GEO>& p, const vec3d&) const {
        return p.radius() - target;
    }
};

/// A spread of start points: L from 1.5 to 8, four longitudes, three latitudes each — 84 points,
/// including the South Atlantic anomaly's longitude and the (−7.492, −2.806, 0) start that sits
/// 0.006 R_E from the magnetic equator and broke the one-probe direction rule.
std::vector<Position<Frame::GEO>> sweep_points() {
    std::vector<Position<Frame::GEO>> out;
    for (const double l : {1.5, 2.0, 3.0, 4.0, 5.0, 6.6, 8.0}) {
        for (const double lon : {0.0, 1.3, 3.5, 5.0}) {
            for (const double zf : {0.0, 0.25, -0.4}) {
                const double r = l * (1.0 - (0.35 * std::abs(zf)));
                out.push_back(geo(r * std::cos(lon), r * std::sin(lon), l * zf));
            }
        }
    }
    return out;
}

// ---------------------------------------------------------------------------------------------
// The oracle, as a black box
// ---------------------------------------------------------------------------------------------

using TraceFieldLine2Fn = void (*)(int*, int*, int*, int*, int*, double*, double*, double*, double*,
                                   double*, double*, double*, double*, double*, double*, double*,
                                   int*);
using TraceTowardEarthFn = void (*)(int*, int*, int*, int*, int*, double*, double*, double*,
                                    double*, double*, double*, double*, int*);
using FindMirrorFn = void (*)(int*, int*, int*, int*, int*, double*, double*, double*, double*,
                              double*, double*, double*, double*, double*);
using FindMagEquatorFn = void (*)(int*, int*, int*, int*, int*, double*, double*, double*, double*,
                                  double*, double*, double*);
using FindFootFn = void (*)(int*, int*, int*, int*, int*, double*, double*, double*, double*,
                            double*, int*, double*, double*, double*, double*);

/// The dev-only oracle handle. Opened once; absent on any machine that has not built IRBEM, in
/// which case every differential test skips rather than fails.
struct Oracle {
    void* handle = nullptr;
    TraceFieldLine2Fn trace_field_line2 = nullptr;
    TraceTowardEarthFn trace_toward_earth = nullptr;
    FindMirrorFn find_mirror = nullptr;
    FindMagEquatorFn find_magequator = nullptr;
    FindFootFn find_foot = nullptr;

    [[nodiscard]] bool ready() const { return handle != nullptr; }
};

const Oracle& oracle() {
    static const Oracle o = [] {
        Oracle out;
        const char* env = std::getenv("CHEATAH_SPACE_IRBEM_ORACLE");
        const std::string path = env != nullptr ? env : "/tmp/irbem-builds/libirbem-O2.so";
        out.handle = dlopen(path.c_str(), RTLD_NOW);
        if (out.handle == nullptr) return out;
        out.trace_field_line2 =
            reinterpret_cast<TraceFieldLine2Fn>(dlsym(out.handle, "trace_field_line2_1_"));
        out.trace_toward_earth = reinterpret_cast<TraceTowardEarthFn>(
            dlsym(out.handle, "trace_field_line_towards_earth1_"));
        out.find_mirror = reinterpret_cast<FindMirrorFn>(dlsym(out.handle, "find_mirror_point1_"));
        out.find_magequator =
            reinterpret_cast<FindMagEquatorFn>(dlsym(out.handle, "find_magequator1_"));
        out.find_foot = reinterpret_cast<FindFootFn>(dlsym(out.handle, "find_foot_point1_"));
        return out;
    }();
    return o;
}

/// IRBEM's refusal sentinel; anything at or below it is "no answer", not a number.
bool is_baddata(double v) { return v <= -0.99e31; }

/// The arguments every oracle entry point shares. `options` is all zeros: no L\*, IGRF internal,
/// default resolutions — the settings these five routines do not read anyway.
struct OracleCall {
    int kext = 0;
    std::array<int, 5> options{0, 0, 0, 0, 0};
    int sysaxes = 1;   // GEO Cartesian, Earth radii
    int year = kYear;
    int doy = kDoy;
    double ut = kUt;
    std::vector<double> maginput = std::vector<double>(25, 0.0);
};

// ---------------------------------------------------------------------------------------------
// Self-contained tests
// ---------------------------------------------------------------------------------------------

TEST(IrbemTraceApi, ParabolicVertexIsExactOnAParabola) {
    // b(s) = 3 + 5 (s - 0.25)^2 sampled at s = -1, 0, +1 with ds = 1. The vertex is at s = 0.25
    // with value 3, and the three-point formula must recover both EXACTLY — every term is a
    // dyadic rational, so this is `==`, not a tolerance.
    const auto b = [](double s) { return 3.0 + (5.0 * (s - 0.25) * (s - 0.25)); };
    double arc = -1.0;
    const double value = detail::parabolic_minimum(b(-1.0), b(0.0), b(1.0), 1.0, arc);
    EXPECT_EQ(arc, 0.25);
    EXPECT_EQ(value, 3.0);

    // A flat bracket has no vertex: the coarse sample stands and the arc is exactly zero, so a
    // caller never takes a partial step on a fit that did not happen.
    double flat_arc = 7.0;
    EXPECT_EQ(detail::parabolic_minimum(2.0, 2.0, 2.0, 1.0, flat_arc), 2.0);
    EXPECT_EQ(flat_arc, 0.0);

    // A CONCAVE bracket (a maximum, not a minimum) is rejected the same way rather than being
    // extrapolated to a vertex that is nowhere near the samples.
    double concave_arc = 7.0;
    EXPECT_EQ(detail::parabolic_minimum(1.0, 5.0, 1.0, 1.0, concave_arc), 5.0);
    EXPECT_EQ(concave_arc, 0.0);
}

TEST(IrbemTraceApi, RefineInStepFindsTheRadiusToRoundoff) {
    const Igrf<10> m = model();
    // Walk down a real field line to the step that actually crosses 1 R_E, then refine inside it —
    // the exact situation the trace's end caps face.
    Position<Frame::GEO> p = geo(4.0, 0.0, 0.0);
    vec3d b = m.evaluate(p).v;
    const double ds = detail::increasing_field_step(m, p, b, 0.08);
    Position<Frame::GEO> beyond = p;
    bool crossed = false;
    for (int i = 0; i < 400 && !crossed; ++i) {
        vec3d b_next{};
        const Position<Frame::GEO> q = detail::rk4_step(m, p, b, ds, b_next);
        if (q.radius() < 1.0) {
            beyond = q;
            crossed = true;
            break;
        }
        p = q;
        b = b_next;
    }
    ASSERT_TRUE(crossed);
    ASSERT_GT(p.radius(), 1.0);

    const detail::RefinedPoint hit = detail::refine_in_step(
        m, p, b, ds, p.radius() - 1.0, beyond.radius() - 1.0, RadiusMinus{1.0}, 8);
    // Four Illinois iterations reach an exact zero here; eight is the shipped setting and leaves
    // nothing. A conditional halving bug held this at 2.5e-6 and this assertion is what found it.
    EXPECT_EQ(hit.position.radius(), 1.0);
    // The root is strictly inside the step, and on the same side of the origin.
    EXPECT_GT(std::abs(hit.arc), 0.0);
    EXPECT_LT(std::abs(hit.arc), std::abs(ds));
    EXPECT_EQ(hit.arc < 0.0, ds < 0.0);
}

TEST(IrbemTraceApi, RefineInStepStaysInsideTheStepWithoutABracket) {
    // No sign change: the secant would extrapolate arbitrarily far and return a point that is not
    // on the traced step at all. Clamping keeps the answer inside the step it was handed.
    const Igrf<10> m = model();
    const Position<Frame::GEO> start = geo(1.6, 0.0, 0.9);
    const vec3d b = m.evaluate(start).v;
    const double ds = 0.4;
    vec3d b_end{};
    const Position<Frame::GEO> end = detail::rk4_step(m, start, b, ds, b_end);
    ASSERT_GT(end.radius(), 1.0) << "premise: neither endpoint is below 1 R_E";
    const detail::RefinedPoint hit = detail::refine_in_step(
        m, start, b, ds, start.radius() - 1.0, end.radius() - 1.0, RadiusMinus{1.0}, 8);
    EXPECT_GE(hit.arc, 0.0);
    EXPECT_LE(hit.arc, ds);
}

TEST(IrbemTraceApi, RefineInStepStopsOnADegenerateBracket) {
    // A predicate with no sign change gives a zero denominator on the first trial; the routine has
    // to return the origin rather than divide by it. This is the branch that keeps a caller who
    // brackets wrongly from getting a NaN position.
    const Igrf<10> m = model();
    const Position<Frame::GEO> start = geo(3.0, 0.0, 0.0);
    const vec3d b = m.evaluate(start).v;
    // A target radius of -1 makes the predicate positive everywhere, so `f_hi - f_lo` is zero on
    // the first trial and the routine must return the origin rather than divide by it.
    const detail::RefinedPoint hit =
        detail::refine_in_step(m, start, b, 0.1, 1.0, 1.0, RadiusMinus{-1.0}, 8);
    EXPECT_EQ(hit.arc, 0.0);
    EXPECT_EQ(hit.position, start);
}

TEST(IrbemTraceApi, StepSizeFollowsTheDipoleLUnlessOverridden) {
    // On the equator the dipole L is the radius, so ds = r / steps_per_l exactly.
    PathTraceOptions opt;
    opt.steps_per_l = 50.0;
    EXPECT_EQ(detail::path_step(geo(4.0, 0.0, 0.0), opt), 4.0 / 50.0);

    // An absolute step wins outright, whatever L is.
    opt.step_size = 0.02;
    EXPECT_EQ(detail::path_step(geo(8.0, 0.0, 0.0), opt), 0.02);

    // A nonsensical steps_per_l falls back to the documented 50 rather than dividing by zero.
    PathTraceOptions bad;
    bad.steps_per_l = 0.0;
    EXPECT_EQ(detail::path_step(geo(4.0, 0.0, 0.0), bad), 4.0 / 50.0);
}

TEST(IrbemTraceApi, IncreasingFieldStepPointsAwayFromTheEquator) {
    const Igrf<10> m = model();
    // A point well north of the magnetic equator: the field rises going further north, which is
    // +B-hat, and the field line's own direction.
    const Position<Frame::GEO> north = geo(3.0, 0.0, 1.5);
    const vec3d bn = m.evaluate(north).v;
    const double ds_n = detail::increasing_field_step(m, north, bn, 0.06);
    vec3d probe{};
    const Position<Frame::GEO> q = detail::rk4_step(m, north, bn, ds_n, probe);
    EXPECT_GT(fixarray::norm(probe), fixarray::norm(bn));
    EXPECT_GT(std::abs(q.v[2]), std::abs(north.v[2]));   // further from the equatorial plane

    // A point well SOUTH must go the other way. Same magnitude, opposite sign.
    const Position<Frame::GEO> south = geo(2.0, 1.0, -1.2);
    const vec3d bs = m.evaluate(south).v;
    const double ds_s = detail::increasing_field_step(m, south, bs, 0.05);
    EXPECT_EQ(std::abs(ds_s), 0.05);
    vec3d probe_s{};
    (void)detail::rk4_step(m, south, bs, ds_s, probe_s);
    EXPECT_GT(fixarray::norm(probe_s), fixarray::norm(bs));
}

TEST(IrbemTraceApi, IncreasingFieldStepIsRightBesideTheMagneticEquator) {
    // REGRESSION. At this start the magnetic equator is 0.006 R_E away, so BOTH neighbours are
    // higher than the start and the natural one-probe rule ("is the +ds neighbour higher than
    // here?") answers yes and walks the wrong way. The consequence was not subtle: FIND_FOOT_POINT
    // returned the northern foot for Hemisphere::Same when the oracle — and the physics — say
    // southern, 130 degrees of latitude away. Comparing the TWO neighbours fixes it, because |B|
    // falls monotonically toward the equator, so the larger neighbour is always the far one.
    const Igrf<10> m = model();
    const Position<Frame::GEO> start = geo(-7.492, -2.806, 0.0);
    const vec3d b = m.evaluate(start).v;

    const auto eq = find_magequator(m, start);
    ASSERT_TRUE(eq.ok());
    ASSERT_LT(separation(eq.value.position, start), 0.01) << "this start must be beside the equator";
    ASSERT_GT(eq.value.position.v[2], start.v[2]) << "and the equator must be to its north";

    // Both neighbours higher — the degeneracy the one-probe rule cannot see.
    vec3d up{};
    vec3d down{};
    (void)detail::rk4_step(m, start, b, 0.02, up);
    (void)detail::rk4_step(m, start, b, -0.02, down);
    EXPECT_GT(fixarray::norm(up), fixarray::norm(b));
    EXPECT_GT(fixarray::norm(down), fixarray::norm(b));

    // The start is SOUTH of the equator, so "away from the equator" is the −B-hat direction.
    EXPECT_LT(detail::increasing_field_step(m, start, b, 0.02), 0.0);

    const auto foot = find_foot_point(m, start, 100.0, Hemisphere::Same);
    ASSERT_TRUE(foot.ok());
    EXPECT_LT(foot.value.position.latitude(), 0.0);
}

TEST(IrbemTraceApi, MirrorPointAtNinetyDegreesIsTheInputPoint) {
    const Igrf<10> m = model();
    const Position<Frame::GEO> start = geo(3.0, 0.0, 1.5);
    const auto mp = find_mirror_point(m, start, 90.0);
    ASSERT_TRUE(mp.ok());
    // Exactly, not nearly: a locally mirroring particle mirrors where it is, and the routine
    // returns the caller's own point rather than a refinement of it.
    EXPECT_EQ(mp.value.position, start);
    EXPECT_EQ(mp.value.b_mirror, mp.value.b_local);
    EXPECT_EQ(mp.value.b_local, m.evaluate(start).magnitude());
}

TEST(IrbemTraceApi, MirrorFieldIsTheAdiabaticInvariantValue) {
    const Igrf<10> m = model();
    const Position<Frame::GEO> start = geo(4.0, 0.0, 0.0);
    const double b_local = m.evaluate(start).magnitude();
    for (const double alpha : {80.0, 60.0, 45.0, 30.0}) {
        const auto mp = find_mirror_point(m, start, alpha);
        ASSERT_TRUE(mp.ok()) << alpha;
        const double sin_a = std::sin(alpha * (std::numbers::pi / 180.0));
        // B_m = B_local / sin^2(alpha) is the first adiabatic invariant, not a search result.
        EXPECT_DOUBLE_EQ(mp.value.b_mirror, b_local / (sin_a * sin_a));
        // ... and the point the trace found must actually carry that field.
        // The refinement lands ON the mirror surface: measured residual 3e-11 nT out of ~10^3.
        EXPECT_NEAR(m.evaluate(mp.value.position).magnitude(), mp.value.b_mirror,
                    1e-12 * mp.value.b_mirror);
    }
}

TEST(IrbemTraceApi, MirrorPointInTheLossConeIsReported) {
    const Igrf<10> m = model();
    // A low-altitude start at a small pitch angle: the mirror field is reached below the surface,
    // so the particle precipitates instead of mirroring. Physics, hence OpenFieldLine.
    const auto mp = find_mirror_point(m, geo(1.2, 0.0, 0.4), 30.0);
    EXPECT_EQ(mp.status, Status::OpenFieldLine);
    EXPECT_GT(mp.value.b_mirror, mp.value.b_local);
    // The value is still populated — a partial answer is diagnostic.
    EXPECT_GT(mp.value.b_local, 0.0);
}

TEST(IrbemTraceApi, MagEquatorIsTheMinimumAlongTheLine) {
    const Igrf<10> m = model();
    std::vector<PathPoint> path(irbem_max_path_points);
    for (const Position<Frame::GEO>& start : sweep_points()) {
        if (start.radius() <= 1.05) continue;
        const auto eq = find_magequator(m, start);
        ASSERT_TRUE(eq.ok());
        const auto line = trace_field_line(m, start, path);
        ASSERT_TRUE(line.ok());
        // No sample on the whole line may be below the fitted minimum, and the fit must sit at or
        // below the smallest sample — that is what "sub-step refinement" means, checked rather
        // than asserted.
        for (std::size_t i = 0; i < line.value.point_count; ++i) {
            EXPECT_GE(path[i].b_magnitude, eq.value.b_min * (1.0 - 1e-9));
        }
        EXPECT_LE(eq.value.b_min, path[line.value.equator_index].b_magnitude);
        // The two routines must agree with each other to the sub-step level.
        EXPECT_NEAR(eq.value.b_min, line.value.b_min, 1e-9 * eq.value.b_min);
        EXPECT_LT(separation(eq.value.position, line.value.equator), 1e-6);
        // And the field really is smallest there. The parabola's vertex VALUE and `|B|` at its
        // vertex POSITION are two different estimates of the same minimum — the fit extrapolates
        // off the sample grid and the position comes from a partial RK4 step — so they agree to
        // the Bmin budget (1e-5, ERROR_BUDGET §4) and not to round-off. Measured: 2.5e-6.
        EXPECT_NEAR(m.evaluate(eq.value.position).magnitude(), eq.value.b_min,
                    1e-5 * eq.value.b_min);
    }
}

TEST(IrbemTraceApi, MagEquatorRefinesAStartThatIsAlreadyTheMinimum) {
    const Igrf<10> m = model();
    // Start ON the equator by construction: run once, then restart from the answer. The second
    // call takes the both-neighbours-higher branch, which is the one a caller hits when it feeds
    // an equatorial ephemeris back in.
    const Position<Frame::GEO> start = geo(5.0, 0.0, 0.0);
    const auto first = find_magequator(m, start);
    ASSERT_TRUE(first.ok());
    const auto second = find_magequator(m, first.value.position);
    ASSERT_TRUE(second.ok());
    EXPECT_NEAR(second.value.b_min, first.value.b_min, 1e-7 * first.value.b_min);
    EXPECT_LT(separation(second.value.position, first.value.position), 1e-3);
}

TEST(IrbemTraceApi, TraceFieldLineEndsOnTheReferenceSurfaceAtBothFeet) {
    const Igrf<10> m = model();
    std::vector<PathPoint> path(irbem_max_path_points);
    for (const double r0 : {1.0, 1.05, 0.95}) {
        PathTraceOptions opt;
        opt.r0 = r0;
        const auto line = trace_field_line(m, geo(4.0, 0.0, 0.0), path, opt);
        ASSERT_TRUE(line.ok()) << r0;
        // Both end caps are regula-falsi refinements, so they land ON the surface rather than one
        // whole step past it. Without the refinement these would be ~0.08 R_E low.
        EXPECT_NEAR(path[0].position.radius(), r0, 1e-11);
        EXPECT_NEAR(path[line.value.point_count - 1].position.radius(), r0, 1e-11);
        // The input point is stored verbatim, and its recorded field is the model's.
        EXPECT_EQ(path[line.value.start_index].position, geo(4.0, 0.0, 0.0));
        EXPECT_EQ(path[line.value.start_index].b_magnitude, line.value.b_local);
    }
}

TEST(IrbemTraceApi, PathIsOrderedAlongTheField) {
    const Igrf<10> m = model();
    std::vector<PathPoint> path(irbem_max_path_points);
    const auto line = trace_field_line(m, geo(2.0, 1.0, -1.2), path);
    ASSERT_TRUE(line.ok());
    ASSERT_GT(line.value.point_count, 10U);
    // Index increases along +B-hat everywhere: the documented ordering, checked at every sample
    // rather than at the ends, because a reversal in the middle is exactly the bug the in-place
    // reverse could introduce.
    for (std::size_t i = 0; i + 1 < line.value.point_count; ++i) {
        const vec3d step = path[i + 1].position.v - path[i].position.v;
        const vec3d b = m.evaluate(path[i].position).v;
        EXPECT_GT(fixarray::dot(step, b), 0.0) << i;
    }
    EXPECT_LT(path[0].position.v[2], path[line.value.point_count - 1].position.v[2])
        << "the geomagnetic field points from the southern foot to the northern one";
}

TEST(IrbemTraceApi, XjAgreesWithTheInvariantTracer) {
    const Igrf<10> m = model();
    std::vector<PathPoint> path(irbem_max_path_points);
    double worst = 0.0;
    for (const Position<Frame::GEO>& start : sweep_points()) {
        if (start.radius() <= 1.05) continue;
        const auto line = trace_field_line(m, start, path);
        if (!line.ok()) continue;
        // lstar.hpp's tracer, at the same resolution, for a locally mirroring particle: the SAME
        // quadrature on the SAME grid. What separates them is only that RK4 is not exactly
        // reversible, so the two walks drift apart at O(ds^5) per step.
        const auto ref = trace_invariant(m, start, 90.0, TraceOptions{});
        if (!ref.ok()) continue;
        if (ref.value.invariant_i > 0.05) {
            worst = std::max(worst, std::abs(line.value.invariant_i - ref.value.invariant_i) /
                                        ref.value.invariant_i);
        }
        EXPECT_NEAR(line.value.invariant_i, ref.value.invariant_i,
                    1e-4 * std::max(1e-3, ref.value.invariant_i));
        EXPECT_NEAR(line.value.b_local, ref.value.b_local, 1e-12 * ref.value.b_local);
    }
    // Measured: 8.6e-6 over the 84-point sweep. The bound is where a real divergence would show.
    EXPECT_LT(worst, 1e-4);
}

TEST(IrbemTraceApi, TraceFieldLineReportsATruncatedPath) {
    const Igrf<10> m = model();
    // Twenty samples cannot hold an L=6 line, and the routine must say so rather than return a
    // "field line" that is a fifth of one with plausible-looking scalars.
    std::vector<PathPoint> tiny(20);
    const auto line = trace_field_line(m, geo(6.0, 0.0, 0.0), tiny);
    EXPECT_EQ(line.status, Status::NotConverged);
    EXPECT_TRUE(line.value.truncated);
    EXPECT_LE(line.value.point_count, tiny.size());

    // The pathological case: the FIRST half alone fills the buffer, so there is no room even for
    // the input point. Reported, not overrun.
    std::vector<PathPoint> tinier(3);
    const auto stub = trace_field_line(m, geo(6.0, 0.0, 0.0), tinier);
    EXPECT_EQ(stub.status, Status::NotConverged);
    EXPECT_TRUE(stub.value.truncated);
    EXPECT_EQ(stub.value.point_count, tinier.size());
}

TEST(IrbemTraceApi, TraceFieldLineTruncatesInTheForwardHalfToo) {
    // BOTH cases in the test above take the SAME exit: the backward half alone fills the buffer,
    // so `trace_field_line` returns before the input point is even stored and every scalar is
    // zero. The other truncation shape — the backward half FITS, the forward half does not — runs
    // the whole quadrature over a partial line and was not reached by any test. It is the shape a
    // caller actually hits, because a buffer that holds one wing usually holds most of the other.
    const Igrf<10> m = model();
    std::vector<PathPoint> full(irbem_max_path_points);
    const auto whole = trace_field_line(m, geo(6.0, 0.0, 0.0), full);
    ASSERT_EQ(whole.status, Status::Ok);
    ASSERT_GT(whole.value.start_index, 1U);

    // One slot past the backward half plus the input point: the forward half gets nothing.
    std::vector<PathPoint> clipped(whole.value.start_index + 1);
    const auto part = trace_field_line(m, geo(6.0, 0.0, 0.0), clipped);
    EXPECT_EQ(part.status, Status::NotConverged);
    EXPECT_TRUE(part.value.truncated);
    EXPECT_EQ(part.value.point_count, clipped.size());
    EXPECT_EQ(part.value.start_index, whole.value.start_index);
    // Unlike the backward-half exit, THESE scalars are computed — over the partial line. Bmin is
    // real (the equator is in the wing that fitted) but the parabola guard declines the fit,
    // because the minimum now sits within two samples of the buffer's end.
    EXPECT_GT(part.value.b_min, 0.0);
    EXPECT_NEAR(part.value.b_min, whole.value.b_min, 1e-3 * whole.value.b_min);
    EXPECT_GT(part.value.invariant_i, 0.0);

    // And the backward-half exit, stated as what it is: nothing was integrated, so nothing is
    // reported. A caller that read these as a short field line would be reading zeros.
    std::vector<PathPoint> stub(whole.value.start_index);
    const auto none = trace_field_line(m, geo(6.0, 0.0, 0.0), stub);
    EXPECT_EQ(none.status, Status::NotConverged);
    EXPECT_TRUE(none.value.truncated);
    EXPECT_EQ(none.value.b_min, 0.0);
    EXPECT_EQ(none.value.invariant_i, 0.0);
    EXPECT_EQ(none.value.mcilwain_l, 0.0);
}

TEST(IrbemTraceApi, TraceFieldLineRefusesADegenerateStart) {
    const Igrf<10> m = model();
    std::vector<PathPoint> path(irbem_max_path_points);
    EXPECT_EQ(trace_field_line(m, geo(4.0, 0.0, 0.0), std::span<PathPoint>{}).status,
              Status::DomainError);
    EXPECT_EQ(trace_field_line(m, geo(0.5, 0.0, 0.0), path).status, Status::DomainError);
    EXPECT_EQ(trace_field_line(m, geo(std::nan(""), 0.0, 0.0), path).status, Status::DomainError);
    EXPECT_EQ(find_magequator(m, geo(0.5, 0.0, 0.0)).status, Status::DomainError);
    EXPECT_EQ(find_mirror_point(m, geo(0.5, 0.0, 0.0), 45.0).status, Status::DomainError);
    EXPECT_EQ(find_mirror_point(m, geo(4.0, 0.0, 0.0), 0.0).status, Status::DomainError);
    EXPECT_EQ(find_mirror_point(m, geo(4.0, 0.0, 0.0), 181.0).status, Status::DomainError);
    EXPECT_EQ(find_mirror_point(m, geo(4.0, 0.0, 0.0), 180.0).status, Status::DomainError)
        << "no perpendicular velocity, so no mirror point - and sin(180 deg) is 1.2e-16, not 0";
    EXPECT_EQ(find_foot_point(m, geo(0.5, 0.0, 0.0), 100.0, Hemisphere::Same).status,
              Status::DomainError);
    EXPECT_EQ(find_foot_point(m, geo(4.0, 0.0, 0.0), std::nan(""), Hemisphere::Same).status,
              Status::DomainError);
    EXPECT_EQ(trace_field_line_toward_earth(m, geo(4.0, 0.0, 0.0), std::span<PathPoint>{}).status,
              Status::DomainError);
    EXPECT_EQ(trace_field_line_toward_earth(m, geo(0.5, 0.0, 0.0), path).status,
              Status::DomainError);
}

TEST(IrbemTraceApi, FootPointRefusesAStartBelowTheRequestedAltitude) {
    const Igrf<10> m = model();
    // 1.05 R_E is ~331 km up; asking for the 1000 km surface from below it is not a trace that
    // terminates, it is a question with no answer on this side of the line.
    EXPECT_EQ(find_foot_point(m, geo(0.55, 0.0, 0.9), 1000.0, Hemisphere::Same).status,
              Status::DomainError);
}

TEST(IrbemTraceApi, FootPointLandsOnTheRequestedGeodeticAltitude) {
    const Igrf<10> m = model();
    double worst = 0.0;
    for (const Position<Frame::GEO>& start : sweep_points()) {
        if (start.radius() <= 1.05) continue;
        for (const double alt : {100.0, 500.0, 1000.0}) {
            for (const Hemisphere h : {Hemisphere::Same, Hemisphere::Opposite, Hemisphere::North,
                                       Hemisphere::South}) {
                const auto foot = find_foot_point(m, start, alt, h);
                if (!foot.ok()) continue;
                worst = std::max(worst, std::abs(foot.value.position.radius() - alt));
                // The termination condition is GEODETIC altitude, so this is the definition of
                // the answer, not a derived quantity. Measured worst residual: 1.9e-10 km, i.e.
                // 0.2 micrometres. The oracle's own iteration stops up to 1.0 km short.
                EXPECT_NEAR(foot.value.position.radius(), alt, 1e-6);
                EXPECT_NEAR(foot.value.b_magnitude, fixarray::norm(foot.value.field.v),
                            1e-12 * foot.value.b_magnitude);
                EXPECT_NEAR(m.evaluate(gdz_to_geo(foot.value.position)).magnitude(),
                            foot.value.b_magnitude, 1e-6 * foot.value.b_magnitude);
            }
        }
    }
    EXPECT_LT(worst, 1e-6);
}

TEST(IrbemTraceApi, FootPointHemisphereFlagsSelectTheTwoFeet) {
    const Igrf<10> m = model();
    for (const Position<Frame::GEO>& start : sweep_points()) {
        if (start.radius() <= 1.05) continue;
        const auto same = find_foot_point(m, start, 100.0, Hemisphere::Same);
        const auto opposite = find_foot_point(m, start, 100.0, Hemisphere::Opposite);
        const auto north = find_foot_point(m, start, 100.0, Hemisphere::North);
        const auto south = find_foot_point(m, start, 100.0, Hemisphere::South);
        if (!same.ok() || !opposite.ok() || !north.ok() || !south.ok()) continue;

        // The two named hemispheres really are north and south...
        EXPECT_GT(north.value.position.latitude(), 0.0);
        EXPECT_LT(south.value.position.latitude(), 0.0);
        // ... Same and Opposite are the same two feet under different names ...
        const bool same_is_north = same.value.position.latitude() > 0.0;
        const Position<Frame::GDZ>& expect_same =
            same_is_north ? north.value.position : south.value.position;
        const Position<Frame::GDZ>& expect_opp =
            same_is_north ? south.value.position : north.value.position;
        EXPECT_NEAR(same.value.position.latitude(), expect_same.latitude(), 1e-9);
        EXPECT_NEAR(opposite.value.position.latitude(), expect_opp.latitude(), 1e-9);
        // ... and they are genuinely different points.
        EXPECT_GT(std::abs(north.value.position.latitude() - south.value.position.latitude()), 1.0);
    }
}

TEST(IrbemTraceApi, TowardEarthSamplesAreOneStepApart) {
    const Igrf<10> m = model();
    std::vector<PathPoint> path(irbem_max_path_points);
    PathTraceOptions opt;
    opt.step_size = 0.02;
    const auto n = trace_field_line_toward_earth(m, geo(4.0, 0.0, 0.0), path, opt);
    ASSERT_TRUE(n.ok());
    ASSERT_GT(n.value, 100U);
    EXPECT_EQ(path[0].position, geo(4.0, 0.0, 0.0));
    // The last sample is the first one INSIDE the surface — the reference's behaviour, and what a
    // plot wants, since a path that stops short of the Earth looks like a bug.
    EXPECT_LT(path[n.value - 1].position.radius(), opt.r0);
    EXPECT_GT(path[n.value - 2].position.radius(), opt.r0);
    for (std::size_t i = 0; i + 1 < n.value; ++i) {
        // The CHORD is slightly shorter than the arc; on a line of this curvature that is a few
        // parts in 10^5, and it must never exceed the step.
        const double chord = separation(path[i + 1].position, path[i].position);
        EXPECT_LE(chord, opt.step_size);
        EXPECT_GT(chord, 0.999 * opt.step_size);
    }
    // It walks toward the Earth, not away from it.
    EXPECT_LT(path[n.value - 1].position.radius(), path[0].position.radius());
}

TEST(IrbemTraceApi, TowardEarthReportsATruncatedPath) {
    const Igrf<10> m = model();
    std::vector<PathPoint> tiny(10);
    PathTraceOptions opt;
    opt.step_size = 0.02;
    const auto n = trace_field_line_toward_earth(m, geo(6.0, 0.0, 0.0), tiny, opt);
    EXPECT_EQ(n.status, Status::NotConverged);
    EXPECT_EQ(n.value, tiny.size());
}

TEST(IrbemTraceApi, TowardEarthStopsAtTheStepCap) {
    const Igrf<10> m = model();
    std::vector<PathPoint> path(irbem_max_path_points);
    PathTraceOptions opt;
    opt.step_size = 0.02;
    opt.max_steps = 5;
    const auto n = trace_field_line_toward_earth(m, geo(6.0, 0.0, 0.0), path, opt);
    EXPECT_EQ(n.status, Status::OpenFieldLine);
    EXPECT_EQ(n.value, 6U);   // the start plus five steps
}

TEST(IrbemTraceApi, TracesStopAtTheStepCap) {
    const Igrf<10> m = model();
    std::vector<PathPoint> path(irbem_max_path_points);
    PathTraceOptions opt;
    opt.max_steps = 4;
    // Four steps cannot reach either foot, so the line does not close and both halves say so.
    const auto line = trace_field_line(m, geo(6.0, 0.0, 0.0), path, opt);
    EXPECT_EQ(line.status, Status::OpenFieldLine);
    EXPECT_EQ(line.value.point_count, 9U);   // four each way plus the start
    EXPECT_FALSE(line.value.truncated);      // room was not the problem
    EXPECT_EQ(find_magequator(m, geo(6.0, 0.0, 3.0), opt).status, Status::OpenFieldLine);
    EXPECT_EQ(find_mirror_point(m, geo(6.0, 0.0, 3.0), 5.0, opt).status, Status::OpenFieldLine);
    EXPECT_EQ(find_foot_point(m, geo(6.0, 0.0, 3.0), 100.0, Hemisphere::Same, opt).status,
              Status::OpenFieldLine);
}

TEST(IrbemTraceApi, MagEquatorReportsAnOpenLine) {
    const Igrf<10> m = model();
    // A start ~60 km up over the SOUTH ATLANTIC ANOMALY (8 deg S, 300 deg E, r = 1.01 R_E). Field
    // lines are arches, so walking toward the field minimum almost always RAISES the geocentric
    // radius — but in the anomaly the surface field is weak and the arch is tilted enough that the
    // line dips as much as 0.012 R_E below its own start before it climbs. Scanning 45 000 starts
    // over the whole globe, that is the only regime where it happens at all, and it is exactly the
    // regime a LEO dosimetry study lives in. The routine reports it instead of integrating through
    // the Earth.
    const auto eq = find_magequator(m, geo(0.500085, -0.866173, -0.140565));
    EXPECT_EQ(eq.status, Status::OpenFieldLine);
    EXPECT_GT(eq.value.b_min, 0.0);
}

TEST(IrbemTraceApi, StepSizeIsClampedNearThePoles) {
    // The dipole L of a near-polar point is r/cos^2(lat), which diverges: at 80 deg S and r = 2
    // R_E it is 66, so the unclamped step would be 1.3 R_E and one RK4 step could land inside the
    // Earth. The clamp holds it to an eighth of the radius.
    PathTraceOptions opt;
    const Position<Frame::GEO> polar = geo(0.347, 0.0, 1.970);   // ~80 deg N, r = 2 R_E
    ASSERT_GT(detail::dipole_l(polar) / opt.steps_per_l, polar.radius() / 8.0);
    EXPECT_EQ(detail::path_step(polar, opt),
              polar.radius() * detail::max_step_fraction_of_radius);
    // ... and it does not bind anywhere in the belts.
    const Position<Frame::GEO> belt = geo(4.0, 0.0, 1.0);
    EXPECT_EQ(detail::path_step(belt, opt), detail::dipole_l(belt) / opt.steps_per_l);
}

TEST(IrbemTraceApi, FootPointReportsAnOpenLine) {
    const Igrf<10> m = model();
    // An altitude below r0 can never be reached: the trace stops at the surface first.
    PathTraceOptions opt;
    opt.r0 = 1.2;
    const auto foot = find_foot_point(m, geo(4.0, 0.0, 0.0), 100.0, Hemisphere::Same, opt);
    EXPECT_EQ(foot.status, Status::OpenFieldLine);
    EXPECT_GT(foot.value.b_magnitude, 0.0);
    // A named hemisphere propagates the same refusal rather than silently trying the other side.
    EXPECT_EQ(find_foot_point(m, geo(4.0, 0.0, 0.0), 100.0, Hemisphere::North, opt).status,
              Status::OpenFieldLine);
}

TEST(IrbemTraceApi, ResultsAreDeterministic) {
    // Bit-exact self-golden: the same inputs must produce the same bits, so a refactor that
    // perturbs the last place is caught even where agreement with the oracle is a tolerance.
    const Igrf<10> m = model();
    std::vector<PathPoint> a(irbem_max_path_points);
    std::vector<PathPoint> b(irbem_max_path_points);
    const Position<Frame::GEO> start = geo(3.0, 0.0, 1.5);
    const auto first = trace_field_line(m, start, a);
    const auto second = trace_field_line(m, start, b);
    ASSERT_EQ(first.value.point_count, second.value.point_count);
    EXPECT_EQ(first.value.invariant_i, second.value.invariant_i);
    EXPECT_EQ(first.value.b_min, second.value.b_min);
    EXPECT_EQ(first.value.mcilwain_l, second.value.mcilwain_l);
    for (std::size_t i = 0; i < first.value.point_count; ++i) {
        EXPECT_EQ(a[i].position, b[i].position) << i;
        EXPECT_EQ(a[i].b_magnitude, b[i].b_magnitude) << i;
    }
}

TEST(IrbemTraceApi, PathBatchValidatesItsArguments) {
    const Igrf<10> m = model();
    const std::vector<Position<Frame::GEO>> starts{geo(4.0, 0.0, 0.0), geo(3.0, 0.0, 1.5)};
    std::vector<PathPoint> paths(2 * 64);
    std::vector<std::uint32_t> counts(2);
    std::vector<Status> statuses(2);
    PathTraceOptions opt;
    opt.step_size = 0.02;

    std::vector<std::uint32_t> short_counts(1);
    EXPECT_EQ(
        trace_field_line_toward_earth_batch(m, starts, paths, short_counts, statuses, opt).status,
        Status::DomainError);
    std::vector<Status> short_status(1);
    EXPECT_EQ(
        trace_field_line_toward_earth_batch(m, starts, paths, counts, short_status, opt).status,
        Status::DomainError);
    // The device lane takes a FIXED step, so an L-proportional one is refused outright rather
    // than silently running a different algorithm on one lane than on the other.
    PathTraceOptions no_step;
    EXPECT_EQ(trace_field_line_toward_earth_batch(m, starts, paths, counts, statuses, no_step)
                  .status,
              Status::DomainError);
    // `paths` must divide evenly into per-line strides, and the stride cannot be zero.
    std::vector<PathPoint> ragged(2 * 64 + 1);
    EXPECT_EQ(
        trace_field_line_toward_earth_batch(m, starts, ragged, counts, statuses, opt).status,
        Status::DomainError);
    EXPECT_EQ(trace_field_line_toward_earth_batch(m, starts, std::span<PathPoint>{}, counts,
                                                  statuses, opt)
                  .status,
              Status::DomainError);
    // An empty batch is a no-op, not an error, and touches no device.
    const auto empty = trace_field_line_toward_earth_batch(
        m, std::span<const Position<Frame::GEO>>{}, std::span<PathPoint>{},
        std::span<std::uint32_t>{}, std::span<Status>{}, opt);
    EXPECT_EQ(empty.status, Status::Ok);
    EXPECT_FALSE(empty.value);
}

TEST(IrbemTraceApi, PathBatchOnTheHostMatchesTheSingleLineRoutine) {
    const Igrf<10> m = model();
    std::vector<Position<Frame::GEO>> starts;
    for (const Position<Frame::GEO>& p : sweep_points()) {
        if (p.radius() > 1.05) starts.push_back(p);
    }
    constexpr std::size_t kStride = 600;
    std::vector<PathPoint> paths(starts.size() * kStride);
    std::vector<std::uint32_t> counts(starts.size());
    std::vector<Status> statuses(starts.size());
    PathTraceOptions opt;
    opt.step_size = 0.02;
    opt.max_steps = 100000;

    const auto batch =
        trace_field_line_toward_earth_batch(m, starts, paths, counts, statuses, opt);
    EXPECT_FALSE(batch.value) << "no device lane is expected in the host-only build";

    // Bit-for-bit: the batch is a loop over the single-line routine, not a second implementation.
    std::vector<PathPoint> one(kStride);
    for (std::size_t i = 0; i < starts.size(); ++i) {
        const auto r = trace_field_line_toward_earth(m, starts[i], one, opt);
        ASSERT_EQ(statuses[i], r.status) << i;
        ASSERT_EQ(counts[i], r.value) << i;
        for (std::size_t k = 0; k < r.value; ++k) {
            EXPECT_EQ(paths[(i * kStride) + k].position, one[k].position) << i << ":" << k;
            EXPECT_EQ(paths[(i * kStride) + k].b_magnitude, one[k].b_magnitude) << i << ":" << k;
        }
    }
}

// ---------------------------------------------------------------------------------------------
// Differential tests against the IRBEM oracle
// ---------------------------------------------------------------------------------------------

#define SKIP_WITHOUT_ORACLE()                                                                 \
    do {                                                                                      \
        if (!oracle().ready())                                                                \
            GTEST_SKIP() << "IRBEM oracle not present; set CHEATAH_SPACE_IRBEM_ORACLE";       \
    } while (false)

/// The resolution the differential comparisons run at.
///
/// NOT the default. `PathTraceOptions::steps_per_l` is 50 to match `lstar.hpp`, and at 50 the
/// first-order quadrature behind `XJ` — and therefore `Lm` — is 8% and 5.7e-3 relative from the
/// oracle's. Both implementations are first-order and NEITHER is converged (ERROR_BUDGET §2(a):
/// "a 0.01 agreement target is only meaningful at matched settings"), so the honest comparison is
/// at a resolution where our discretization is small against the difference being measured. 400
/// steps per L is where our own value stops moving faster than the oracle's residual. Measured,
/// over the 84-point × 4-epoch sweep:
///
///     steps_per_l :        50 :       100 :       200 :       400
///     max |dLm|/Lm:   5.7e-03 :   2.8e-03 :   1.3e-03 :   6.0e-04
///     max |dXJ|/XJ:   8.2e-02 :   4.7e-02 :   2.2e-02 :   9.2e-03
///     max |dBmin|/B:  3.9e-06 :   1.7e-06 :   2.0e-06 :   1.9e-06
///
/// `Bmin` is flat because it is not a quadrature — it converged before 50.
constexpr double kDifferentialStepsPerL = 400.0;

PathTraceOptions differential_options() {
    PathTraceOptions opt;
    opt.steps_per_l = kDifferentialStepsPerL;
    opt.max_steps = 100000;
    return opt;
}

TEST(IrbemTraceApi, OracleEvaluatesTheSameFieldWeDo) {
    SKIP_WITHOUT_ORACLE();
    // The premise of every test below. If this fails, the differential numbers are measuring a
    // coefficient difference and nothing else.
    const Igrf<10> m = model();
    OracleCall c;
    double worst = 0.0;
    for (const Position<Frame::GEO>& p : sweep_points()) {
        if (p.radius() <= 1.05) continue;
        double x1 = p.v[0];
        double x2 = p.v[1];
        double x3 = p.v[2];
        double b_local = 0.0;
        double b_mirror = 0.0;
        double alpha = 90.0;
        std::array<double, 3> xgeo{};
        oracle().find_mirror(&c.kext, c.options.data(), &c.sysaxes, &c.year, &c.doy, &c.ut, &x1,
                             &x2, &x3, &alpha, c.maginput.data(), &b_local, &b_mirror, xgeo.data());
        ASSERT_FALSE(is_baddata(b_local));
        worst = std::max(worst, std::abs(m.evaluate(p).magnitude() - b_local) / b_local);
    }
    // Measured: 1.2e-15 — the two libraries evaluate IGRF-10 at year+0.5 to round-off.
    EXPECT_LT(worst, 1e-13) << "matched-model premise broken; check the epoch convention";
}

TEST(IrbemTraceApi, MagEquatorMatchesTheOracle) {
    SKIP_WITHOUT_ORACLE();
    const Igrf<10> m = model();
    const PathTraceOptions opt = differential_options();
    OracleCall c;
    double worst_b = 0.0;
    double worst_pos = 0.0;
    int compared = 0;
    for (const Position<Frame::GEO>& p : sweep_points()) {
        if (p.radius() <= 1.05) continue;
        double x1 = p.v[0];
        double x2 = p.v[1];
        double x3 = p.v[2];
        double bmin = 0.0;
        std::array<double, 3> xgeo{};
        oracle().find_magequator(&c.kext, c.options.data(), &c.sysaxes, &c.year, &c.doy, &c.ut, &x1,
                                 &x2, &x3, c.maginput.data(), &bmin, xgeo.data());
        const auto eq = find_magequator(m, p, opt);
        if (is_baddata(bmin) || !eq.ok()) continue;
        ++compared;
        worst_b = std::max(worst_b, std::abs(eq.value.b_min - bmin) / bmin);
        worst_pos = std::max(worst_pos,
                             separation(eq.value.position, Position<Frame::GEO>{vec3d{
                                                               xgeo[0], xgeo[1], xgeo[2]}}));
    }
    EXPECT_GT(compared, 80);
    // Measured: Bmin 1.9e-6 relative, position 4.6e-5 R_E. The Bmin residual is the ORACLE's: our
    // value at this resolution is within 1e-8 of our own converged value, and the oracle sits the
    // same 1.9e-6 away from it. Budget for Bmin is 1e-5 (ERROR_BUDGET §4).
    EXPECT_LT(worst_b, 1e-5);
    EXPECT_LT(worst_pos, 2e-4);
}

TEST(IrbemTraceApi, MirrorPointMatchesTheOracle) {
    SKIP_WITHOUT_ORACLE();
    const Igrf<10> m = model();
    const PathTraceOptions opt = differential_options();
    OracleCall c;
    double worst_pos = 0.0;
    double worst_bmirr_here = 0.0;
    int compared = 0;
    for (const Position<Frame::GEO>& p : sweep_points()) {
        if (p.radius() <= 1.05) continue;
        for (const double alpha : {90.0, 75.0, 60.0, 45.0, 30.0, 15.0}) {
            double x1 = p.v[0];
            double x2 = p.v[1];
            double x3 = p.v[2];
            double a = alpha;
            double b_local = 0.0;
            double b_mirror = 0.0;
            std::array<double, 3> xgeo{};
            oracle().find_mirror(&c.kext, c.options.data(), &c.sysaxes, &c.year, &c.doy, &c.ut,
                                 &x1, &x2, &x3, &a, c.maginput.data(), &b_local, &b_mirror,
                                 xgeo.data());
            const auto mp = find_mirror_point(m, p, alpha, opt);
            if (is_baddata(b_mirror) || !mp.ok()) continue;
            ++compared;
            const Position<Frame::GEO> theirs{vec3d{xgeo[0], xgeo[1], xgeo[2]}};
            worst_pos = std::max(worst_pos, separation(mp.value.position, theirs));
            // THE ORACLE'S `Bmirr` IS NOT THE ADIABATIC VALUE. Evaluating our field at its own
            // returned position reproduces its Bmirr to 1.3e-15, so what it reports is |B| where
            // its root-find stopped, and the 1.8e-4 spread against B_local/sin^2(alpha) is that
            // root-find's residual. Ours is the invariant value exactly. Measuring the two claims
            // separately is the only way to say which is which.
            worst_bmirr_here = std::max(
                worst_bmirr_here, std::abs(m.evaluate(theirs).magnitude() - b_mirror) / b_mirror);
        }
    }
    EXPECT_GT(compared, 400);
    // Measured: position 2.1e-4 R_E (~1.3 km along a line whose step is 0.02 R_E).
    EXPECT_LT(worst_pos, 1e-3);
    EXPECT_LT(worst_bmirr_here, 1e-12);
}

TEST(IrbemTraceApi, FootPointMatchesTheOracle) {
    SKIP_WITHOUT_ORACLE();
    const Igrf<10> m = model();
    const PathTraceOptions opt = differential_options();
    OracleCall c;
    double worst_lat = 0.0;
    double worst_lon = 0.0;
    double worst_b = 0.0;
    double worst_their_alt = 0.0;
    int compared = 0;
    for (const Position<Frame::GEO>& p : sweep_points()) {
        if (p.radius() <= 1.05) continue;
        for (const double alt : {100.0, 500.0}) {
            for (const int flag : {0, 1, -1, 2}) {
                double x1 = p.v[0];
                double x2 = p.v[1];
                double x3 = p.v[2];
                double stop_alt = alt;
                int hemi = flag;
                std::array<double, 3> xfoot{};
                std::array<double, 3> bfoot{};
                double bfootmag = 0.0;
                oracle().find_foot(&c.kext, c.options.data(), &c.sysaxes, &c.year, &c.doy, &c.ut,
                                   &x1, &x2, &x3, &stop_alt, &hemi, c.maginput.data(),
                                   xfoot.data(), bfoot.data(), &bfootmag);
                if (is_baddata(bfootmag)) continue;
                worst_their_alt = std::max(worst_their_alt, std::abs(xfoot[0] - alt));
                // MATCHED-ALTITUDE COMPARISON. The oracle's own iteration stops up to 1.0 km from
                // the altitude it was asked for, so comparing its foot with ours at the REQUESTED
                // altitude measures its termination error, not the two tracers. Asking for the
                // altitude it actually reached puts both feet on the same surface.
                const auto foot =
                    find_foot_point(m, p, xfoot[0], static_cast<Hemisphere>(flag), opt);
                if (!foot.ok()) continue;
                ++compared;
                worst_lat =
                    std::max(worst_lat, std::abs(foot.value.position.latitude() - xfoot[1]));
                double dlon = std::abs(foot.value.position.longitude() - xfoot[2]);
                if (dlon > 180.0) dlon = 360.0 - dlon;
                worst_lon = std::max(worst_lon, dlon);
                worst_b = std::max(worst_b, std::abs(foot.value.b_magnitude - bfootmag) / bfootmag);
            }
        }
    }
    EXPECT_GT(compared, 500);
    // Measured, at matched altitude: latitude 5.8e-4 deg (inside the 1e-3 deg footpoint budget),
    // longitude 3.6e-3 deg — which at |lat| ~ 60 is ~0.2 km of ground distance, i.e. about twice
    // the budget's ~100 m and dominated by the two tracers' along-line disagreement rather than by
    // either one's termination. |B| agrees to 4.8e-6, inside the 1e-5 B-at-a-point budget.
    EXPECT_LT(worst_lat, 2e-3);
    EXPECT_LT(worst_lon, 1e-2);
    EXPECT_LT(worst_b, 1e-5);
    // Ours lands within 6.3e-4 km of the altitude asked for; the oracle's within 1.0 km.
    EXPECT_GT(worst_their_alt, 0.1) << "if the oracle got tighter, retune the matched comparison";
}

TEST(IrbemTraceApi, TraceFieldLineMatchesTheOracle) {
    SKIP_WITHOUT_ORACLE();
    const Igrf<10> m = model();
    const PathTraceOptions opt = differential_options();
    OracleCall c;
    std::vector<PathPoint> path(200000);
    std::vector<double> their_b(irbem_max_path_points);
    std::vector<double> their_pos(3 * irbem_max_path_points);
    double worst_lm = 0.0;
    double worst_lm_abs = 0.0;
    double worst_xj = 0.0;
    double worst_bmin = 0.0;
    double worst_path = 0.0;
    int compared = 0;

    for (const Position<Frame::GEO>& p : sweep_points()) {
        if (p.radius() <= 1.05) continue;
        double x1 = p.v[0];
        double x2 = p.v[1];
        double x3 = p.v[2];
        double r0 = 1.0;
        double lm = 0.0;
        double bmin = 0.0;
        double xj = 0.0;
        int count = 0;
        oracle().trace_field_line2(&c.kext, c.options.data(), &c.sysaxes, &c.year, &c.doy, &c.ut,
                                   &x1, &x2, &x3, c.maginput.data(), &r0, &lm, their_b.data(),
                                   &bmin, &xj, their_pos.data(), &count);
        const auto line = trace_field_line(m, p, path, opt);
        if (is_baddata(lm) || !line.ok() || count < 4) continue;
        ++compared;
        worst_lm = std::max(worst_lm, std::abs(line.value.mcilwain_l - lm) / lm);
        worst_lm_abs = std::max(worst_lm_abs, std::abs(line.value.mcilwain_l - lm));
        worst_bmin = std::max(worst_bmin, std::abs(line.value.b_min - bmin) / bmin);
        if (xj > 0.05) worst_xj = std::max(worst_xj, std::abs(line.value.invariant_i - xj) / xj);

        // GEOMETRY. The two implementations sample the line at different spacings, so comparing
        // index for index is meaningless; what must agree is the CURVE. Every fortieth of our
        // samples is projected onto the oracle's polyline and the distance recorded.
        const std::size_t stride = (line.value.point_count / 40) + 1;
        for (std::size_t i = 0; i < line.value.point_count; i += stride) {
            double best = 1e30;
            for (int k = 0; k + 1 < count; ++k) {
                const vec3d a{their_pos[3 * k], their_pos[(3 * k) + 1], their_pos[(3 * k) + 2]};
                const vec3d bpt{their_pos[(3 * k) + 3], their_pos[(3 * k) + 4],
                                their_pos[(3 * k) + 5]};
                const vec3d u = bpt - a;
                const vec3d w = path[i].position.v - a;
                const double uu = fixarray::dot(u, u);
                double t = uu > 0.0 ? fixarray::dot(w, u) / uu : 0.0;
                t = std::clamp(t, 0.0, 1.0);
                best = std::min(best, fixarray::norm(w - (u * t)));
            }
            worst_path = std::max(worst_path, best);
        }
    }

    EXPECT_GT(compared, 80);
    // Measured over 84 points x 4 epochs, at steps_per_l = 400:
    //   Lm      6.0e-4 relative, 6.1e-3 absolute (worst at Lm ~ 10.7; 3.4e-3 absolute for Lm<=6.6)
    //   Bmin    6.3e-8 relative
    //   XJ      9.2e-3 relative for XJ > 0.05; 1.6e-2 R_E absolute at worst
    //   path    1.8e-4 R_E point-to-polyline
    //
    // The Lm ABSOLUTE budget of 1e-3 (ERROR_BUDGET §4) is MET only below Lm ~ 3; at Lm ~ 10.7 the
    // gap is 6.1e-3. It is inherited whole from XJ: Lm is Hilton's closed form in (I, B_m), and
    // both libraries integrate I with a first-order rule whose error at their own step sizes is
    // ~1e-2 relative. Neither is converged, so this is a RESOLUTION disagreement and not an
    // algorithm one — our own value moves by the same order between steps_per_l 200 and 400 (see
    // kDifferentialStepsPerL's table). Closing it needs a higher-order quadrature that handles the
    // square-root endpoint singularity, which is `Compat::Improved`'s job and not this header's.
    EXPECT_LT(worst_lm, 1e-3);
    EXPECT_LT(worst_lm_abs, 1e-2);
    EXPECT_LT(worst_bmin, 1e-6);
    EXPECT_LT(worst_xj, 2e-2);
    EXPECT_LT(worst_path, 1e-3);
}

TEST(IrbemTraceApi, TowardEarthMatchesTheOracle) {
    SKIP_WITHOUT_ORACLE();
    const Igrf<10> m = model();
    OracleCall c;
    std::vector<PathPoint> path(irbem_max_path_points);
    std::vector<double> their_pos(3 * irbem_max_path_points);
    PathTraceOptions opt;
    opt.step_size = 0.02;
    opt.max_steps = 100000;
    double worst = 0.0;
    long worst_count = 0;
    int compared = 0;

    for (const Position<Frame::GEO>& p : sweep_points()) {
        if (p.radius() <= 1.05) continue;
        double x1 = p.v[0];
        double x2 = p.v[1];
        double x3 = p.v[2];
        double ds = 0.02;   // IRBEM WRITES THROUGH THIS POINTER: it flips the sign to go earthward
        int count = 0;
        oracle().trace_toward_earth(&c.kext, c.options.data(), &c.sysaxes, &c.year, &c.doy, &c.ut,
                                    &x1, &x2, &x3, c.maginput.data(), &ds, their_pos.data(),
                                    &count);
        const auto n = trace_field_line_toward_earth(m, p, path, opt);
        if (count < 3 || !n.ok()) continue;
        ++compared;
        worst_count = std::max(worst_count, std::abs(static_cast<long>(n.value) - count));
        for (std::size_t i = 0; i < n.value && static_cast<int>(i) < count; ++i) {
            worst = std::max(worst, separation(path[i].position,
                                               Position<Frame::GEO>{vec3d{their_pos[3 * i],
                                                                          their_pos[(3 * i) + 1],
                                                                          their_pos[(3 * i) + 2]}}));
        }
    }
    EXPECT_GT(compared, 80);
    // Measured: sample counts IDENTICAL in all 336 cases, and the whole path agrees to 3.2e-15
    // R_E — over ~500 steps. At a matched fixed step, with a matched field, this reproduces the
    // reference's tracer to round-off, which is the strongest statement in this file: the RK4
    // scheme, the step convention and the termination rule are the same algorithm and not merely
    // close ones.
    EXPECT_EQ(worst_count, 0);
    EXPECT_LT(worst, 1e-12);
}

// ---------------------------------------------------------------------------------------------
// The device lane
//
// Compiled only when cheatah-gpu-linalg was on the include path (the same __has_include seam
// lstar.hpp uses), and SKIPPED at runtime when no device comes up. A silent host fallback is the
// failure mode that makes a GPU performance claim worthless, so the lane a call actually took is
// asserted rather than assumed.
// ---------------------------------------------------------------------------------------------

#ifdef CHEATAH_SPACE_IRBEM_LSTAR_GPU

/// Starts spread over L = 2..8 and a range of latitudes — the shape a drift-shell sweep has, not a
/// row of identical lines that would let every warp retire together.
std::vector<Position<Frame::GEO>> device_starts(std::size_t n) {
    std::vector<Position<Frame::GEO>> out;
    out.reserve(n);
    for (std::size_t i = 0; i < n; ++i) {
        const double t = static_cast<double>(i) / static_cast<double>(n);
        const double l = 2.0 + (6.0 * t);
        const double lon = 6.2831853071795864 * std::fmod(t * 7.0, 1.0);
        const double lat = 0.6 * std::sin(t * 31.0);
        const double r = l * std::cos(lat) * std::cos(lat);
        out.push_back(geo(r * std::cos(lat) * std::cos(lon), r * std::cos(lat) * std::sin(lon),
                          r * std::sin(lat)));
    }
    return out;
}

#define SKIP_WITHOUT_DEVICE()                                                                  \
    do {                                                                                       \
        if (!gpu::available()) GTEST_SKIP() << "no device: " << gpu::unavailable_reason();      \
    } while (false)

TEST(IrbemTraceApiGpu, PathBatchUsesTheDeviceWhenOneIsAvailable) {
    SKIP_WITHOUT_DEVICE();
    const Igrf<10> m = model();
    const std::size_t n = gpu::gpu_crossover("irbem_trace_path_f32");
    const std::vector<Position<Frame::GEO>> starts = device_starts(n);
    constexpr std::size_t kMax = 512;
    std::vector<PathPoint> paths(n * kMax);
    std::vector<std::uint32_t> counts(n);
    std::vector<Status> statuses(n);
    PathTraceOptions opt;
    opt.step_size = 0.02;
    const auto r = trace_field_line_toward_earth_batch(m, starts, paths, counts, statuses, opt);
    EXPECT_TRUE(r.value) << "the batch fell back to the host at its own crossover";
    EXPECT_EQ(r.status, Status::Ok);
}

TEST(IrbemTraceApiGpu, PathKernelAgreesWithTheHostLane) {
    SKIP_WITHOUT_DEVICE();
    const Igrf<10> m = model();
    const std::size_t n = std::max<std::size_t>(gpu::gpu_crossover("irbem_trace_path_f32"), 512);
    const std::vector<Position<Frame::GEO>> starts = device_starts(n);
    constexpr std::size_t kMax = 512;
    PathTraceOptions opt;
    opt.step_size = 0.02;

    std::vector<PathPoint> device_paths(n * kMax);
    std::vector<std::uint32_t> device_counts(n);
    std::vector<Status> device_status(n);
    const auto on_device =
        trace_field_line_toward_earth_batch(m, starts, device_paths, device_counts, device_status,
                                            opt);
    ASSERT_TRUE(on_device.value);

    double worst_first = 0.0;
    double worst_last = 0.0;
    double worst_b = 0.0;
    long worst_count = 0;
    std::vector<PathPoint> host(kMax);
    for (std::size_t i = 0; i < n; ++i) {
        const auto h = trace_field_line_toward_earth(m, starts[i], host, opt);
        if (!h.ok() || device_status[i] != Status::Ok) continue;
        worst_count = std::max(worst_count,
                               std::abs(static_cast<long>(device_counts[i]) -
                                        static_cast<long>(h.value)));
        const std::size_t common = std::min<std::size_t>(device_counts[i], h.value);
        for (std::size_t k = 0; k < common; ++k) {
            const double d = separation(device_paths[(i * kMax) + k].position, host[k].position);
            if (k <= 1) worst_first = std::max(worst_first, d);
            if (k + 1 == common) worst_last = std::max(worst_last, d);
            worst_b = std::max(worst_b, std::abs(device_paths[(i * kMax) + k].b_magnitude -
                                                 host[k].b_magnitude) /
                                            host[k].b_magnitude);
        }
    }
    // fp32 POSITION accumulated over ~250 steps: the error grows along the path, which is the
    // whole reason this kernel's precision story differs from irbem_trace_i_f32's. The invariant
    // tracer accumulates one scalar; this one accumulates a position. Measured numbers are quoted
    // in the report the run prints.
    std::printf("[ path kernel ] fp32 vs fp64: first sample %.3e Re, last sample %.3e Re, "
                "|B| %.3e rel, sample-count delta %ld\n",
                worst_first, worst_last, worst_b, worst_count);
    EXPECT_LT(worst_first, 1e-6);
    EXPECT_LT(worst_last, 1e-3);
    EXPECT_LE(worst_count, 1);
}

/// Set an environment variable for the duration of a scope and put it back.
class ScopedEnv {
public:
    ScopedEnv(const char* name, const char* value) : name_(name) {
        if (const char* prev = std::getenv(name)) {
            had_ = true;
            prev_ = prev;
        }
        setenv(name, value, 1);
    }
    ScopedEnv(const ScopedEnv&) = delete;
    ScopedEnv& operator=(const ScopedEnv&) = delete;
    ScopedEnv(ScopedEnv&&) = delete;
    ScopedEnv& operator=(ScopedEnv&&) = delete;
    ~ScopedEnv() {
        if (had_) {
            setenv(name_, prev_.c_str(), 1);
            return;
        }
        unsetenv(name_);
    }

private:
    const char* name_;
    bool had_ = false;
    std::string prev_;
};

TEST(IrbemTraceApiGpu, PathBatchFallsBackToTheHostLane) {
    SKIP_WITHOUT_DEVICE();
    const Igrf<10> m = model();
    const std::size_t n = 512;
    const std::vector<Position<Frame::GEO>> starts = device_starts(n);
    constexpr std::size_t kMax = 512;
    std::vector<PathPoint> paths(n * kMax);
    std::vector<std::uint32_t> counts(n);
    std::vector<Status> statuses(n);
    PathTraceOptions opt;
    opt.step_size = 0.02;
    {
        // The operator switch: the same machine, the same batch, forced onto the host. This is
        // what makes the differential lane-vs-lane comparison possible at all.
        const ScopedEnv off("CHEATAH_SPACE_IRBEM_NO_GPU", "1");
        const auto r = trace_field_line_toward_earth_batch(m, starts, paths, counts, statuses, opt);
        EXPECT_FALSE(r.value);
        EXPECT_GT(counts[0], 0U);
    }
    {
        // A missing .spv is a DEPLOYMENT problem, not a reason to refuse to compute: the batch
        // must quietly run on the host rather than throw or return zeros.
        const ScopedEnv nowhere("CHEATAH_SPACE_IRBEM_SPV_DIR",
                                "/nonexistent-space-irbem-shader-dir");
        const ScopedEnv always("CHEATAH_SPACE_IRBEM_GPU_CROSSOVER", "1");
        const auto r = trace_field_line_toward_earth_batch(m, starts, paths, counts, statuses, opt);
        EXPECT_FALSE(r.value);
        EXPECT_GT(counts[0], 0U);
    }
}

TEST(IrbemTraceApiGpu, PathBatchHonoursTheStepCapOnEitherLane) {
    SKIP_WITHOUT_DEVICE();
    // A REGRESSION, not a hypothetical. `PathTraceOptions::max_steps` bounds the host loop, and
    // the kernel's dims buffer originally did not carry it: the device loop ran to `max_points`
    // instead. At max_steps = 10 every one of 1 024 lines came back different — device 80 samples
    // and `Ok`, host 11 and `OpenFieldLine` — so the ANSWER, not just the timing, depended on
    // whether a device happened to be present. That is the one thing a batch entry point with a
    // fallback must never do.
    const Igrf<10> m = model();
    const std::size_t n = std::max<std::size_t>(gpu::gpu_crossover("irbem_trace_path_f32"), 512);
    const std::vector<Position<Frame::GEO>> starts = device_starts(n);
    constexpr std::size_t kMax = 512;
    PathTraceOptions opt;
    opt.step_size = 0.02;
    opt.max_steps = 10;

    std::vector<PathPoint> dev(n * kMax);
    std::vector<std::uint32_t> dev_counts(n);
    std::vector<Status> dev_status(n);
    const auto on_device =
        trace_field_line_toward_earth_batch(m, starts, dev, dev_counts, dev_status, opt);
    ASSERT_TRUE(on_device.value);

    std::vector<PathPoint> hst(n * kMax);
    std::vector<std::uint32_t> hst_counts(n);
    std::vector<Status> hst_status(n);
    {
        const ScopedEnv off("CHEATAH_SPACE_IRBEM_NO_GPU", "1");
        const auto on_host =
            trace_field_line_toward_earth_batch(m, starts, hst, hst_counts, hst_status, opt);
        ASSERT_FALSE(on_host.value);
    }
    for (std::size_t i = 0; i < n; ++i) {
        ASSERT_EQ(dev_counts[i], hst_counts[i]) << "line " << i << " ran a different length";
        ASSERT_EQ(dev_status[i], hst_status[i]) << "line " << i;
    }
    // Eleven samples: the input point plus ten steps, then the cap.
    EXPECT_EQ(dev_counts[0], 11U);
    EXPECT_EQ(dev_status[0], Status::OpenFieldLine);
}

TEST(IrbemTraceApiGpu, PathKernelCrossover) {
    SKIP_WITHOUT_DEVICE();
    if (std::getenv("CHEATAH_SPACE_IRBEM_BENCH") == nullptr) {
        GTEST_SKIP() << "set CHEATAH_SPACE_IRBEM_BENCH=1 to run the crossover measurement";
    }
    const Igrf<10> m = model();
    constexpr std::size_t kMax = 512;
    PathTraceOptions opt;
    opt.step_size = 0.02;
    std::printf("%10s %14s %14s %10s\n", "lines", "device us/line", "host us/line", "speedup");
    const std::size_t sizes[] = {64, 128, 256, 512, 1024, 4096, 16384, 65536};
    for (const std::size_t n : sizes) {
        const std::vector<Position<Frame::GEO>> starts = device_starts(n);
        std::vector<PathPoint> paths(n * kMax);
        std::vector<std::uint32_t> counts(n);
        std::vector<Status> statuses(n);

        double best_device = 1e30;
        double best_host = 1e30;
        for (int rep = 0; rep < 5; ++rep) {
            {
                setenv("CHEATAH_SPACE_IRBEM_GPU_CROSSOVER", "1", 1);
                const auto t0 = std::chrono::steady_clock::now();
                const auto r =
                    trace_field_line_toward_earth_batch(m, starts, paths, counts, statuses, opt);
                const auto t1 = std::chrono::steady_clock::now();
                unsetenv("CHEATAH_SPACE_IRBEM_GPU_CROSSOVER");
                if (r.value) {
                    best_device = std::min(best_device,
                                           std::chrono::duration<double>(t1 - t0).count());
                }
            }
            {
                setenv("CHEATAH_SPACE_IRBEM_NO_GPU", "1", 1);
                const auto t0 = std::chrono::steady_clock::now();
                const auto r =
                    trace_field_line_toward_earth_batch(m, starts, paths, counts, statuses, opt);
                const auto t1 = std::chrono::steady_clock::now();
                unsetenv("CHEATAH_SPACE_IRBEM_NO_GPU");
                ASSERT_FALSE(r.value);
                best_host = std::min(best_host, std::chrono::duration<double>(t1 - t0).count());
            }
        }
        const double dn = static_cast<double>(n);
        std::printf("%10zu %14.3f %14.3f %10.2fx\n", n, best_device * 1e6 / dn,
                    best_host * 1e6 / dn, best_host / best_device);
    }
}

#endif  // CHEATAH_SPACE_IRBEM_LSTAR_GPU

}  // namespace
