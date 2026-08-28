// Host-only coverage for space.irbem's field-query guards and float-step helper.
//
// These paths are pure host arithmetic — an input validation guard, a degenerate-field fall-through
// and a float-rounding helper — but their existing tests sit behind a GPU availability check they
// do not actually need. scripts/coverage.sh configures WITHOUT the device lane, so those tests skip
// there and the lines read as uncovered, which is how a genuinely untested guard and a merely
// device-gated one become indistinguishable.
//
// Requiring a GPU to cover code that never touches one is the wrong dependency. These run anywhere.
#include <gtest/gtest.h>

#include <cmath>
#include <limits>

#include "space/irbem/api.hpp"
#include "space/irbem/driftshell.hpp"
#include "space/irbem/field.hpp"
#include "space/irbem/lstar.hpp"
#include "space/irbem/igrf.hpp"

namespace ib = cheatah::space::irbem;
namespace fx = cheatah::fixarray;

namespace {
/// A model with a switchable field: identically zero, or a plain axial dipole.
///
/// The zero case is the degenerate branch no real internal field can produce away from the origin,
/// and a fall-through nobody can reach from real data is still code the compiler emits and a future
/// edit can break. But a model that ONLY returns zero can never reach the routine's main body, so
/// this instantiation would be half-tested — and llvm-cov counts it separately from every other
/// model. One switchable model drives both halves of its own instantiation.
struct SwitchableField {
    bool zero = false;

    [[nodiscard]] ib::FieldVector<ib::Frame::GEO> evaluate(
        const ib::Position<ib::Frame::GEO>& p) const {
        if (zero) return ib::FieldVector<ib::Frame::GEO>{fx::vec3d{0.0, 0.0, 0.0}};
        // B = g10 (a/r)^3 (3 (z_hat . r_hat) r_hat - z_hat), the centred axial dipole.
        const double r2 = fx::squared_norm(p.v);
        if (r2 <= 0.0) return ib::FieldVector<ib::Frame::GEO>{fx::vec3d{0.0, 0.0, 0.0}};
        const double r = std::sqrt(r2);
        const double inv_r5 = 1.0 / (r2 * r2 * r);
        constexpr double g10 = -30000.0;
        return ib::FieldVector<ib::Frame::GEO>{
            fx::vec3d{3.0 * g10 * p.v[0] * p.v[2] * inv_r5, 3.0 * g10 * p.v[1] * p.v[2] * inv_r5,
                      g10 * ((3.0 * p.v[2] * p.v[2]) - r2) * inv_r5}};
    }
};
}  // namespace

TEST(IrbemFieldHostPaths, BderivsRefusesWhatItCannotAnswer) {
    const auto model = ib::Igrf<13>::at(2015.0);
    ASSERT_TRUE(model.has_value());
    const ib::Position<ib::Frame::GEO> good{fx::vec3d{3.0, 0.0, 0.0}};

    // The origin: r == 0, so there is no direction to difference along.
    EXPECT_EQ(ib::Status::DomainError,
              ib::bderivs(*model, ib::Position<ib::Frame::GEO>{fx::vec3d{0.0, 0.0, 0.0}}).status);
    // A non-finite coordinate cannot produce a finite derivative, and must say so rather than
    // propagating a NaN into a trace where it would surface hundreds of steps later.
    const double nan = std::numeric_limits<double>::quiet_NaN();
    EXPECT_EQ(ib::Status::DomainError,
              ib::bderivs(*model, ib::Position<ib::Frame::GEO>{fx::vec3d{nan, 0.0, 0.0}}).status);
    // A negative or non-finite step is a caller error, not a small step.
    EXPECT_EQ(ib::Status::DomainError, ib::bderivs(*model, good, -1e-3).status);
    EXPECT_EQ(ib::Status::DomainError, ib::bderivs(*model, good, nan).status);
    // ...and the same call with a valid step must succeed, so the guards are shown to be rejecting
    // the input rather than the routine being broken.
    EXPECT_EQ(ib::Status::Ok, ib::bderivs(*model, good, 1e-3).status);
}

TEST(IrbemFieldHostPaths, HemisphereRefusesWhatItCannotAnswer) {
    const auto model = ib::Igrf<13>::at(2015.0);
    ASSERT_TRUE(model.has_value());
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const ib::Position<ib::Frame::GEO> good{fx::vec3d{3.0, 0.0, 0.0}};

    EXPECT_EQ(ib::Status::DomainError,
              ib::hemisphere(*model, ib::Position<ib::Frame::GEO>{fx::vec3d{0.0, 0.0, 0.0}}).status);
    EXPECT_EQ(ib::Status::DomainError,
              ib::hemisphere(*model, ib::Position<ib::Frame::GEO>{fx::vec3d{nan, 1.0, 0.0}}).status);
    EXPECT_EQ(ib::Status::DomainError, ib::hemisphere(*model, good, -1.0).status);
    EXPECT_EQ(ib::Status::Ok, ib::hemisphere(*model, good).status);
}

TEST(IrbemFieldHostPaths, AZeroFieldHasNoHemisphere) {
    // |B| == 0 gives no direction to step along and no slope to read, so the answer is Invalid —
    // reported as Ok, because a null field is a legitimate model output, not a caller error. That
    // distinction is the whole reason Status and the value are separate.
    const SwitchableField null{true};
    const auto r = ib::hemisphere(null, ib::Position<ib::Frame::GEO>{fx::vec3d{3.0, 0.0, 0.0}});
    EXPECT_EQ(ib::Status::Ok, r.status);
    EXPECT_EQ(ib::Hemisphere::Invalid, r.value);

    // The SAME instantiation, now with a real field, so its main body is compiled AND run rather
    // than only its degenerate branch.
    const SwitchableField live{false};
    const auto north = ib::hemisphere(live, ib::Position<ib::Frame::GEO>{fx::vec3d{2.0, 0.0, 2.0}});
    const auto south = ib::hemisphere(live, ib::Position<ib::Frame::GEO>{fx::vec3d{2.0, 0.0, -2.0}});
    EXPECT_EQ(ib::Status::Ok, north.status);
    EXPECT_EQ(ib::Status::Ok, south.status);
    EXPECT_NE(north.value, south.value) << "opposite sides of the equator must not agree";
    EXPECT_EQ(ib::Status::Ok, ib::hemisphere(live, {fx::vec3d{2.0, 0.0, 2.0}}, 1e-3).status);
    EXPECT_EQ(ib::Status::DomainError,
              ib::hemisphere(live, {fx::vec3d{0.0, 0.0, 0.0}}).status);
}

TEST(IrbemFieldHostPaths, InterleavePacksPositionsForTheDeviceUpload) {
    // detail::interleave is the doubles-to-floats packing every device dispatch feeds from. It is
    // NOT compiled out without a GPU — it is ordinary host code that simply has no caller there,
    // which is how a function reads as 0% covered while looking device-specific. It is worth
    // testing directly: a transposed index here would corrupt every position the device sees, and
    // the symptom would appear as a physics error hundreds of RK4 steps downstream.
    const std::array<ib::Position<ib::Frame::GEO>, 3> pts{
        ib::Position<ib::Frame::GEO>{fx::vec3d{1.5, -2.25, 3.0}},
        ib::Position<ib::Frame::GEO>{fx::vec3d{-4.0, 5.5, -6.75}},
        ib::Position<ib::Frame::GEO>{fx::vec3d{0.0, 0.25, -0.5}}};
    std::array<float, 9> packed{};
    ib::detail::interleave(pts, packed);

    // Exactly-representable values throughout, so the double-to-float narrowing is exact and the
    // comparison can be == rather than approximate.
    for (std::size_t i = 0; i < pts.size(); ++i) {
        EXPECT_FLOAT_EQ(static_cast<float>(pts[i].v[0]), packed[(3 * i) + 0]) << i;
        EXPECT_FLOAT_EQ(static_cast<float>(pts[i].v[1]), packed[(3 * i) + 1]) << i;
        EXPECT_FLOAT_EQ(static_cast<float>(pts[i].v[2]), packed[(3 * i) + 2]) << i;
    }
    // Empty input must be a no-op, not an out-of-bounds write into a zero-length span.
    std::array<float, 0> none{};
    ib::detail::interleave({}, none);
}

TEST(IrbemFieldHostPaths, TheRealisedStepIsTheOneTheFloatLaneCanTake) {
    // The device lane differences in fp32, so the step it ACTUALLY takes is whatever survives
    // rounding x and x+h to float. Dividing by the intended h rather than the realised one is a
    // silent error in every fp32 derivative, growing as h shrinks — which is exactly when a caller
    // believes they are being more accurate.
    EXPECT_DOUBLE_EQ(0.0625, ib::detail::realised_step(1.0, 0.0625));  // exactly representable: unchanged

    // A step small enough to vanish into x's last float bit: the two floats collapse, the realised
    // step is zero, and the helper returns the intended h so the caller's division stays finite
    // rather than producing an infinity.
    EXPECT_DOUBLE_EQ(1e-30, ib::detail::realised_step(1.0, 1e-30));

    // In between, the realised step differs from the intended one and is what must be divided by.
    const double h = 1e-6;
    const double taken = ib::detail::realised_step(3.0, h);
    EXPECT_NE(0.0, taken);
    EXPECT_NEAR(h, taken, h) << "the realised step should be the same order as the intended one";
}

// ---------------------------------------------------------------------------------------------
// The same guards, for the OTHER instantiation.
//
// llvm-cov counts a template per INSTANTIATION, so a guard exercised only through `Igrf<13>` is an
// untested guard in `Igrf<1>` — and they are separately compiled code that a future edit can break
// independently. The truncated-degree model is not a test fixture either: `Igrf<1>` IS the centred
// dipole, which the drift-shell tests use as their closed-form reference, so its guards are on a
// real path rather than a hypothetical one.
// ---------------------------------------------------------------------------------------------

/// Drive every guard for one truncation degree. Called once per degree the library instantiates,
/// because llvm-cov counts a template PER INSTANTIATION — a guard exercised only through `Igrf<13>`
/// is untested code in `Igrf<10>`, and they are separately compiled bodies a future edit can break
/// independently. The degrees are not arbitrary: 13 is IGRF-14's full published degree, 10 is the
/// truncation IRBEM uses (so it is the one the differential tests run through), and 1 is the
/// centred dipole the drift-shell tests use as their closed-form reference.
template <int NMAX>
void exercise_guards_at_degree() {
    const auto model = ib::Igrf<NMAX>::at(2015.0);
    ASSERT_TRUE(model.has_value()) << "degree " << NMAX;
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const ib::Position<ib::Frame::GEO> good{fx::vec3d{3.0, 0.0, 0.0}};
    const ib::Position<ib::Frame::GEO> origin{fx::vec3d{0.0, 0.0, 0.0}};

    EXPECT_EQ(ib::Status::DomainError, ib::bderivs(*model, origin).status) << NMAX;
    EXPECT_EQ(ib::Status::DomainError, ib::bderivs(*model, good, nan).status) << NMAX;
    EXPECT_EQ(ib::Status::DomainError, ib::bderivs(*model, good, -1.0).status) << NMAX;
    EXPECT_EQ(ib::Status::Ok, ib::bderivs(*model, good, 1e-3).status) << NMAX;

    EXPECT_EQ(ib::Status::DomainError, ib::hemisphere(*model, origin).status) << NMAX;
    EXPECT_EQ(ib::Status::DomainError, ib::hemisphere(*model, good, -1.0).status) << NMAX;
    EXPECT_EQ(ib::Status::DomainError, ib::hemisphere(*model, good, nan).status) << NMAX;
    EXPECT_EQ(ib::Status::Ok, ib::hemisphere(*model, good).status) << NMAX;
    // The fp32 step path: passing 0 asks the routine to choose, which is the branch the device lane
    // takes and the one a host-only run would otherwise never compile a caller for.
    EXPECT_EQ(ib::Status::Ok, ib::hemisphere(*model, good, 0.0).status) << NMAX;

    // The trace's step cap, per degree.
    ib::TraceOptions opt;
    opt.max_steps = 2;
    EXPECT_EQ(ib::Status::OpenFieldLine,
              ib::trace_invariant(*model, ib::Position<ib::Frame::GEO>{fx::vec3d{6.0, 0.0, 0.0}},
                                  30.0, opt)
                  .status)
        << NMAX;

    // Mismatched batch spans, per degree.
    std::array<ib::Position<ib::Frame::GEO>, 2> starts{good, good};
    std::array<double, 1> pitch{45.0};
    std::array<ib::FieldLine, 2> out{};
    std::array<ib::Status, 2> sts{};
    EXPECT_EQ(ib::Status::DomainError,
              ib::trace_invariant_batch(*model, starts, pitch, out, sts).status)
        << NMAX;

    // An epoch outside IGRF's definition, per degree.
    EXPECT_EQ(ib::Status::DomainError, ib::api::rotations_at(1850, 1, 0.0, *model).status) << NMAX;
    EXPECT_EQ(ib::Status::Ok, ib::api::rotations_at(2015, 180, 43200.0, *model).status) << NMAX;
}

TEST(IrbemFieldHostPaths, TheGuardsHoldAtEveryTruncationDegree) {
    exercise_guards_at_degree<13>();  // IGRF-14's full degree
    exercise_guards_at_degree<10>();  // IRBEM's truncation — the differential path
    exercise_guards_at_degree<1>();   // the centred dipole reference
}

TEST(IrbemFieldHostPaths, ATraceThatRunsOutOfStepsSaysSo) {
    // The step cap is not an error condition — it is how an OPEN field line is detected, since a
    // line that never closes would otherwise integrate forever. Forcing it with a tiny cap makes
    // the path deterministic rather than waiting for a pathological geometry to produce it.
    const auto model = ib::Igrf<13>::at(2015.0);
    ASSERT_TRUE(model.has_value());
    ib::TraceOptions opt;
    opt.max_steps = 2;  // far too few to reach either mirror point
    const auto r = ib::trace_invariant(*model, ib::Position<ib::Frame::GEO>{fx::vec3d{6.0, 0.0, 0.0}},
                                       30.0, opt);
    EXPECT_EQ(ib::Status::OpenFieldLine, r.status);
    // The partial integral is still returned: a truncated trace is diagnostic, and silently
    // zeroing it would hide how far the trace got.
    EXPECT_GE(r.value.invariant_i, 0.0);

    const auto dipole = ib::Igrf<1>::at(2015.0);
    ASSERT_TRUE(dipole.has_value());
    EXPECT_EQ(ib::Status::OpenFieldLine,
              ib::trace_invariant(*dipole, ib::Position<ib::Frame::GEO>{fx::vec3d{6.0, 0.0, 0.0}},
                                  30.0, opt)
                  .status);
}

TEST(IrbemFieldHostPaths, MismatchedBatchSpansAreRefused) {
    // Spans of disagreeing length cannot be a partial answer — there is no correspondence between
    // input and output to preserve — so this is refused outright rather than processed to the
    // shortest, which would silently drop points.
    const auto model = ib::Igrf<13>::at(2015.0);
    ASSERT_TRUE(model.has_value());
    std::array<ib::Position<ib::Frame::GEO>, 2> starts{
        ib::Position<ib::Frame::GEO>{fx::vec3d{3.0, 0.0, 0.0}},
        ib::Position<ib::Frame::GEO>{fx::vec3d{4.0, 0.0, 0.0}}};
    std::array<double, 1> pitch{45.0};  // deliberately one short
    std::array<ib::FieldLine, 2> out{};
    std::array<ib::Status, 2> sts{};
    EXPECT_EQ(ib::Status::DomainError,
              ib::trace_invariant_batch(*model, starts, pitch, out, sts).status);
}

TEST(IrbemFieldHostPaths, RotationsRefuseAnEpochTheFieldModelCannotAnswer) {
    // IGRF is defined from 1900 to 2030. Outside that the coefficients would be an unbounded
    // linear extrapolation of a five-year secular-variation prediction, which is not a field —
    // so the epoch is refused rather than the rotations being built from nonsense.
    const auto model = ib::Igrf<13>::at(2015.0);
    ASSERT_TRUE(model.has_value());
    EXPECT_EQ(ib::Status::DomainError, ib::api::rotations_at(1850, 1, 0.0, *model).status);
    EXPECT_EQ(ib::Status::DomainError, ib::api::rotations_at(2100, 1, 0.0, *model).status);
    EXPECT_EQ(ib::Status::Ok, ib::api::rotations_at(2015, 180, 43200.0, *model).status);
}

// ---------------------------------------------------------------------------------------------
// The drift shell at IRBEM's own truncation degree.
//
// make_lstar_batch<13> sits at 98.82% lines while make_lstar_batch<10> sits at 89.35% — the same
// source, separately compiled, one exercised and one not. Degree 10 is not an incidental
// instantiation: it is the truncation IRBEM uses, so it is the code path every DIFFERENTIAL claim
// in this module actually runs through. Leaving it thinner-covered than the degree the closed-form
// tests use is exactly backwards.
// ---------------------------------------------------------------------------------------------

TEST(IrbemFieldHostPaths, TheDriftShellWorksAtIrbemsTruncationDegree) {
    const auto model = ib::Igrf<10>::at(2015.0);
    ASSERT_TRUE(model.has_value());
    const auto rot = ib::api::rotations_at(2015, 180, 43200.0, *model);
    ASSERT_EQ(ib::Status::Ok, rot.status);

    // A batch spanning the inner belt to beyond geosynchronous, equatorial and off, so the shell
    // walk and the flux integral both run rather than just the entry guards.
    const std::array<ib::Position<ib::Frame::GEO>, 4> starts{
        ib::Position<ib::Frame::GEO>{fx::vec3d{3.0, 0.0, 0.0}},
        ib::Position<ib::Frame::GEO>{fx::vec3d{4.5, 0.0, 0.0}},
        ib::Position<ib::Frame::GEO>{fx::vec3d{0.0, 6.6, 0.0}},
        ib::Position<ib::Frame::GEO>{fx::vec3d{5.0, 0.0, 1.0}}};
    const std::array<double, 4> pitch{90.0, 60.0, 90.0, 45.0};
    std::array<ib::DriftShell, 4> out{};
    std::array<ib::Status, 4> sts{};

    const auto r = ib::make_lstar_batch(*model, rot.value, starts, pitch, out, sts);
    EXPECT_NE(ib::Status::DomainError, r.status);
    for (std::size_t i = 0; i < starts.size(); ++i) {
        if (sts[i] != ib::Status::Ok) continue;
        // L* must at least be a plausible shell parameter: greater than 1 (inside the Earth is not
        // a drift shell) and not wildly beyond the starting radius.
        EXPECT_GT(out[i].lstar, 1.0) << "point " << i;
        EXPECT_LT(out[i].lstar, 20.0) << "point " << i;
    }

    // The single-point form, which wraps the batch — same degree, so its wrapper is compiled and
    // exercised too rather than only its degree-13 twin.
    const auto one = ib::make_lstar(*model, rot.value,
                                    ib::Position<ib::Frame::GEO>{fx::vec3d{4.0, 0.0, 0.0}}, 90.0);
    EXPECT_NE(ib::Status::DomainError, one.status);

    // And BOTH batch guards at this degree — the span-length check and the options check are
    // separate code in each instantiation, and this is the degree the differential tests use.
    std::array<double, 2> short_pitch{90.0, 90.0};
    EXPECT_EQ(ib::Status::DomainError,
              ib::make_lstar_batch(*model, rot.value, starts, short_pitch, out, sts).status);

    ib::DriftShellOptions too_few;
    too_few.azimuths = 2;  // cannot enclose a contour
    sts.fill(ib::Status::Ok);  // poison, so a missed per-point write is visible
    EXPECT_EQ(ib::Status::DomainError,
              ib::make_lstar_batch(*model, rot.value, starts, pitch, out, sts, too_few).status);
    for (const ib::Status st : sts) EXPECT_EQ(ib::Status::DomainError, st);
}

TEST(IrbemFieldHostPaths, NonsensicalDriftShellOptionsAreRefusedPerPoint) {
    // The options guard: fewer than three azimuths cannot enclose a contour, fewer than two
    // bracket trials cannot bracket a root, and a non-positive colatitude step cannot advance the
    // flux quadrature. None of these is a small value — each makes the algorithm meaningless — so
    // they are refused rather than clamped to something workable, which would silently return an
    // answer the caller did not ask for.
    //
    // Every per-point status is set too, not just the aggregate: a caller iterating statuses must
    // not read a stale or default-constructed Ok for a batch that was refused outright.
    const auto model = ib::Igrf<13>::at(2015.0).value();
    const auto rot = ib::api::rotations_at(2015, 180, 43200.0, model);
    ASSERT_EQ(ib::Status::Ok, rot.status);

    const std::array<ib::Position<ib::Frame::GEO>, 2> starts{
        ib::Position<ib::Frame::GEO>{fx::vec3d{4.0, 0.0, 0.0}},
        ib::Position<ib::Frame::GEO>{fx::vec3d{5.0, 0.0, 0.0}}};
    const std::array<double, 2> pitch{90.0, 90.0};
    std::array<ib::DriftShell, 2> out{};
    std::array<ib::Status, 2> sts{};

    for (const ib::DriftShellOptions bad : {
             [] { ib::DriftShellOptions o; o.azimuths = 2; return o; }(),
             [] { ib::DriftShellOptions o; o.bracket_trials = 1; return o; }(),
             [] { ib::DriftShellOptions o; o.colatitude_step_deg = 0.0; return o; }(),
             [] { ib::DriftShellOptions o; o.colatitude_step_deg = -1.0; return o; }(),
         }) {
        sts = {ib::Status::Ok, ib::Status::Ok};  // poison, so a missed write is visible
        const auto r = ib::make_lstar_batch(model, rot.value, starts, pitch, out, sts, bad);
        EXPECT_EQ(ib::Status::DomainError, r.status);
        EXPECT_FALSE(r.value) << "a refused batch cannot have used the device";
        for (const ib::Status s : sts) EXPECT_EQ(ib::Status::DomainError, s);
    }

    // The same options at their defaults must succeed, so the guard is shown to reject the OPTIONS
    // rather than the routine being broken for this input.
    const auto ok = ib::make_lstar_batch(model, rot.value, starts, pitch, out, sts);
    EXPECT_NE(ib::Status::DomainError, ok.status);
}
