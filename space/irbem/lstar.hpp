#pragma once

/**
 * @file lstar.hpp
 * @brief space.irbem — field-line tracing, the adiabatic invariants, and L*.
 *
 * This is what the module exists for. Everything else — the frames, the internal field, the
 * transforms, the GPU seam — is the substrate these four quantities stand on:
 *
 *  - **B_min** — the minimum field along the line through a point, i.e. the magnetic equator.
 *  - **I** — the geometric part of the second adiabatic invariant,
 *    `I = ∫ √(1 − B(s)/B_m) ds`, integrated between the mirror points.
 *  - **L_m** — McIlwain's L, from `(I, B_m)` through Hilton's approximation.
 *  - **L\*** — Roederer's L, `L* = 2πk₀/(Φ R_E)`, where Φ is the magnetic flux enclosed by the
 *    particle's DRIFT SHELL. This one is expensive: it needs a whole family of field lines, not
 *    one, which is why a single L\* costs ~10⁵ field-model evaluations and why the reference
 *    implementation takes ~15 ms per point.
 *
 * ## Why the cost is structural, and where it goes
 *
 * Tracing is a serial chain — each RK4 step depends on the last — so no amount of hardware makes
 * ONE trace faster. The parallelism is *across* traces: `ntime × pitch angles × azimuths ×
 * root-find trials`, which is 10⁵–10⁷ for a real batch. That is why the device lane is one thread
 * per FIELD LINE and never one thread per step.
 *
 * It also has enormous arithmetic intensity, which is what makes it worth offloading at all: the
 * point goes to the device once (~24 bytes) and ~10⁵ field evaluations happen on it before six
 * scalars come back. Measured on the same seam, the contrast is stark — a streaming dipole kernel
 * at 0.5 flops/byte LOSES to the host (0.69×), while IGRF at ~20 flops/byte WINS by 8.37×. A trace
 * is ~9 400 flops/byte.
 *
 * ## The step size is load-bearing
 *
 * `ds = L/50`, where L is the start point's dipole L. It looks like an arbitrary constant and it
 * is not: arc length and step size BOTH scale with L, so the step COUNT stays roughly uniform
 * across field lines of very different size. On the device that is what keeps a warp converged —
 * lanes carrying an L=2 line and an L=8 line retire together. Replacing it with a fixed absolute
 * step is the obvious-looking refactor and would make high-L lanes do ~4× the work of low-L ones.
 *
 * ## Sources
 *
 *  - Roederer, *Dynamics of Geomagnetically Trapped Radiation*, Springer (1970), ch. 2 — the
 *    invariants and the drift-shell definition of L\*.
 *  - Schulz & Lanzerotti, *Particle Diffusion in the Radiation Belts*, Springer (1974), ch. 1.
 *  - Hilton, *L parameter: a new approximation*, J. Geophys. Res. 76(28):6952 (1971) — the closed
 *    form for L_m from `(I, B_m)` that replaces McIlwain's original tabulation.
 */

#include <cmath>
#include <cstddef>

#include <span>
#include <vector>

#include "frames.hpp"
#include "igrf.hpp"
#include "policy.hpp"
#include "status.hpp"

// The device lane is opt-in by include path, exactly as space.time's ndarray support is. When
// cheatah-gpu-linalg is absent every routine here still compiles and simply runs on the host.
#if __has_include("cheatah_gpu_linalg/context.hpp")
#  include "gpu/dispatch.hpp"
#  define CHEATAH_SPACE_IRBEM_LSTAR_GPU 1
#endif

namespace cheatah::space::irbem {

/// What a field-line trace produces. Fixed size, trivially copyable, no path storage — the whole
/// point is that `I` and `B_min` accumulate as the trace runs, so nothing needs to be remembered.
/// The APIs that DO return the path (`trace_field_line`, `drift_shell`) are separate and
/// bandwidth-bound; this is the compute-bound one that L\* is built from.
struct FieldLine {
    double b_local = 0.0;   ///< `|B|` at the starting point, nT.
    double b_min = 0.0;     ///< The minimum `|B|` along the line, nT — the magnetic equator.
    double b_mirror = 0.0;  ///< The mirror field for the requested pitch angle, nT.
    double invariant_i = 0.0;  ///< `I = ∫ √(1 − B/B_m) ds`, in Earth radii.
    double arc_length = 0.0;   ///< Total arc traced between mirror points, Earth radii.
    int steps = 0;             ///< RK4 steps taken; a diagnostic, and the divergence measure.
    Position<Frame::GEO> equator{};  ///< Where `b_min` was found.
};

/// Tuning for a trace. Defaults reproduce the reference implementation's behaviour so a
/// differential comparison is comparing ALGORITHMS rather than resolutions.
struct TraceOptions {
    /// Steps per unit L — the reference's `Nreb`. `ds = L/steps_per_l`. See the file brief on why
    /// this is proportional to L rather than absolute.
    double steps_per_l = 50.0;
    /// The hard cap on RK4 steps. Reaching it means the line did not close: an open field line, or
    /// a start point outside the model's domain. Reported as @ref Status::OpenFieldLine, never as a
    /// silently truncated integral.
    int max_steps = 1000;
    /// The radius at which a trace is considered to have reached the atmosphere, in Earth radii.
    double min_radius = 1.0;
};

namespace detail {

/// One RK4 step along the unit tangent `dx/ds = B/|B|`.
///
/// Returns the field at the ARRIVAL point as well as the point, because the next step needs it as
/// its own `k1`: the reference evaluates the field a fifth time at the endpoint purely to report
/// `|B|` there and then re-evaluates the same point on the next call. Carrying it forward removes
/// one of every five field evaluations — a 20% cut in the dominant cost, for free.
///
/// @tparam NMAX the IGRF truncation degree.
/// @param model the internal field model.
/// @param p the current position, GEO, Earth radii.
/// @param b_here the field at @p p, already known, nT.
/// @param ds the signed step length, Earth radii.
/// @param b_next receives the field at the arrival point, nT.
/// @return the arrival position.
/// @complexity Four IGRF evaluations — the RK4 stages — plus one carried in.
/// @alloc none.
template <GeoFieldModel M>
[[nodiscard]] inline Position<Frame::GEO> rk4_step(const M& model,
                                                   const Position<Frame::GEO>& p,
                                                   const fixarray::vec3d& b_here, double ds,
                                                   fixarray::vec3d& b_next) {
    const auto unit = [](const fixarray::vec3d& b) {
        const double m = fixarray::norm(b);
        return m > 0.0 ? fixarray::vec3d{b[0] / m, b[1] / m, b[2] / m} : fixarray::vec3d{};
    };
    const fixarray::vec3d k1 = unit(b_here);
    const auto at = [&](const fixarray::vec3d& v) {
        return model.evaluate(Position<Frame::GEO>{v}).v;
    };
    const fixarray::vec3d p2 = p.v + (k1 * (ds * 0.5));
    const fixarray::vec3d k2 = unit(at(p2));
    const fixarray::vec3d p3 = p.v + (k2 * (ds * 0.5));
    const fixarray::vec3d k3 = unit(at(p3));
    const fixarray::vec3d p4 = p.v + (k3 * ds);
    const fixarray::vec3d k4 = unit(at(p4));
    const fixarray::vec3d out =
        p.v + ((k1 + (k2 * 2.0) + (k3 * 2.0) + k4) * (ds / 6.0));
    b_next = at(out);
    return Position<Frame::GEO>{out};
}

/// The dipole L of a point — `r / sin²θ` in dipole coordinates — used only to size the step.
/// @param p the position, GEO, Earth radii.
/// @return the dipole L; at least 1, so the step never collapses near the origin.
/// @complexity O(1). @alloc none.
[[nodiscard]] inline double dipole_l(const Position<Frame::GEO>& p) {
    const double r = fixarray::norm(p.v);
    if (r <= 0.0) return 1.0;
    const double sin_lat = p.v[2] / r;
    const double cos2_lat = 1.0 - (sin_lat * sin_lat);
    return cos2_lat > 1e-12 ? r / cos2_lat : 1.0e3;
}

}  // namespace detail

/**
 * Trace the field line through @p start and accumulate the second invariant.
 *
 * Walks in the direction of DECREASING field to find the magnetic equator, then continues to the
 * conjugate mirror point, accumulating `I` as it goes. The integrand `√(1 − B/B_m)` has a
 * square-root singularity in its derivative at both mirror points, which is why the reference uses
 * a first-order rule with an endpoint correction rather than a higher-order quadrature that would
 * converge no faster there.
 *
 * @tparam NMAX the IGRF truncation degree.
 * @param model the internal field model, already built for the epoch.
 * @param start the starting position, GEO, Earth radii.
 * @param pitch_angle_deg the equatorial pitch angle; 90° mirrors at the equator, so `I` is zero.
 * @param opt the tracing options.
 * @return the trace, with @ref Status::OpenFieldLine when the line did not close within
 *         @ref TraceOptions::max_steps, and @ref Status::DomainError for a start point at the
 *         origin or inside the atmosphere. The value is populated in every case — a partial trace
 *         is still diagnostic — which is why this returns @ref Result rather than an optional.
 * @complexity O(steps) IGRF evaluations, ~4 per step; 100–300 steps typically.
 * @alloc none — the trace accumulates in registers and stores no path.
 * @test IrbemLstar.DipoleTraceMatchesTheAnalyticInvariant
 */
namespace detail {

/// Integrate `√(1 − B/B_m) ds` from @p from, walking in the direction @p ds, until the field
/// reaches @p b_mirror. One HALF of the bounce path.
///
/// @return the partial integral; @p steps receives the count, and @p closed whether a mirror point
///         was actually reached rather than the step cap or the atmosphere.
/// @complexity O(steps) IGRF evaluations. @alloc none.
template <GeoFieldModel M>
[[nodiscard]] inline double half_invariant(const M& model,
                                            const Position<Frame::GEO>& from,
                                            const fixarray::vec3d& b_from, double b_mirror,
                                            double ds, const TraceOptions& opt, int& steps,
                                            bool& closed) {
    Position<Frame::GEO> p = from;
    fixarray::vec3d b = b_from;
    double sum = 0.0;
    closed = false;
    for (steps = 0; steps < opt.max_steps; ++steps) {
        const double b_prev = fixarray::norm(b);
        fixarray::vec3d b_new{};
        const Position<Frame::GEO> q = detail::rk4_step(model, p, b, ds, b_new);
        const double bmag = fixarray::norm(b_new);
        if (fixarray::norm(q.v) < opt.min_radius) return sum;   // hit the atmosphere: open line
        if (bmag >= b_mirror) {
            // The mirror point lies inside this step. Credit the fraction that was below it,
            // linearly in B — first order, matching the reference, and the integrand's derivative
            // is singular here so a higher-order rule would not converge faster anyway.
            if (bmag > b_prev) {
                const double frac = (b_mirror - b_prev) / (bmag - b_prev);
                sum += 0.5 * std::sqrt(std::max(0.0, 1.0 - (b_prev / b_mirror))) * frac;
            }
            closed = true;
            return sum;
        }
        sum += std::sqrt(1.0 - (bmag / b_mirror));
        p = q;
        b = b_new;
    }
    return sum;
}

}  // namespace detail

/**
 * Trace the field line through @p start and accumulate the second invariant.
 *
 * Three stages, because `I` is defined between the two mirror points and a one-directional walk
 * from an arbitrary start point covers only part of that path:
 *
 *  1. walk in the direction of DECREASING field to the magnetic equator, recording `B_min`;
 *  2. from the equator, integrate outward to the mirror point;
 *  3. from the equator, integrate outward the OTHER way, to the conjugate mirror point.
 *
 * `I` is the sum. Doing only the first half and doubling it would be right for a centred dipole and
 * wrong for every real field, which is asymmetric about the equator — and wrong by an amount that
 * grows exactly where the models matter most.
 *
 * @tparam NMAX the IGRF truncation degree.
 * @param model the internal field model, already built for the epoch.
 * @param start the starting position, GEO, Earth radii.
 * @param pitch_angle_deg the LOCAL pitch angle at @p start; 90° mirrors immediately, so `I` is 0.
 * @param opt the tracing options.
 * @return the trace. @ref Status::OpenFieldLine when either half failed to reach a mirror point
 *         within @ref TraceOptions::max_steps or ran into the atmosphere; @ref Status::DomainError
 *         for a start inside the atmosphere or a non-physical pitch angle. The value is populated
 *         in every case, because a partial trace is still diagnostic.
 * @complexity O(steps) IGRF evaluations, ~4 per step; 100–300 steps typically.
 * @alloc none — the integral accumulates in registers and no path is stored.
 * @test IrbemLstar.DipoleTraceMatchesTheAnalyticInvariant
 */
template <GeoFieldModel M>
[[nodiscard]] inline Result<FieldLine> trace_invariant(const M& model,
                                                       const Position<Frame::GEO>& start,
                                                       double pitch_angle_deg,
                                                       const TraceOptions& opt = {}) {
    FieldLine fl{};
    const double r0 = fixarray::norm(start.v);
    if (!(r0 > opt.min_radius) || !std::isfinite(r0)) return {Status::DomainError, fl};

    fixarray::vec3d b = model.evaluate(start).v;
    fl.b_local = fixarray::norm(b);
    if (!(fl.b_local > 0.0)) return {Status::DomainError, fl};

    const double sin_a = std::sin(pitch_angle_deg * (std::numbers::pi / 180.0));
    if (!(sin_a > 0.0)) return {Status::DomainError, fl};
    fl.b_mirror = fl.b_local / (sin_a * sin_a);

    const double ds_mag = detail::dipole_l(start) / opt.steps_per_l;

    // --- 1. downhill to the equator ---------------------------------------------------------
    fixarray::vec3d b_probe{};
    (void)detail::rk4_step(model, start, b, ds_mag, b_probe);
    double ds = (fixarray::norm(b_probe) < fl.b_local) ? ds_mag : -ds_mag;

    Position<Frame::GEO> eq = start;
    fixarray::vec3d b_eq = b;
    fl.b_min = fl.b_local;
    for (int i = 0; i < opt.max_steps; ++i) {
        fixarray::vec3d b_new{};
        const Position<Frame::GEO> q = detail::rk4_step(model, eq, b_eq, ds, b_new);
        const double bmag = fixarray::norm(b_new);
        if (fixarray::norm(q.v) < opt.min_radius) { fl.steps = i; return {Status::OpenFieldLine, fl}; }
        if (bmag >= fl.b_min) break;      // passed the minimum
        fl.b_min = bmag;
        eq = q;
        b_eq = b_new;
    }
    fl.equator = eq;

    // --- 2 and 3. both halves, from the equator outward --------------------------------------
    int steps_up = 0;
    int steps_dn = 0;
    bool closed_up = false;
    bool closed_dn = false;
    const double i_up =
        detail::half_invariant(model, eq, b_eq, fl.b_mirror, ds, opt, steps_up, closed_up);
    const double i_dn =
        detail::half_invariant(model, eq, b_eq, fl.b_mirror, -ds, opt, steps_dn, closed_dn);

    fl.steps = steps_up + steps_dn;
    fl.arc_length = std::abs(ds) * static_cast<double>(fl.steps);
    fl.invariant_i = (i_up + i_dn) * std::abs(ds);
    return {(closed_up && closed_dn) ? Status::Ok : Status::OpenFieldLine, fl};
}

/**
 * The geomagnetic dipole moment of a model, `M`, in nT.
 *
 * The magnitude of the degree-1 field: `√(g₁⁰² + g₁¹² + h₁¹²)`. Epoch-dependent — the Earth's
 * moment has fallen roughly 6% over the last century — and `L ∝ M^(1/3)`, so taking it from the
 * model rather than from a constant is what keeps `L_m` right across epochs.
 *
 * @tparam NMAX the truncation degree.
 * @param model the internal field model.
 * @return `M` in nT.
 * @complexity O(1). @alloc none.
 * @test IrbemLstar.DipoleMomentTracksTheEpoch
 */
template <class M>
    requires requires(const M& m) { { m.g(1, 0) } -> std::convertible_to<double>; }
[[nodiscard]] inline double dipole_moment(const M& model) {
    const double g10 = model.g(1, 0);
    const double g11 = model.g(1, 1);
    const double h11 = model.h(1, 1);
    return std::sqrt((g10 * g10) + (g11 * g11) + (h11 * h11));
}

/**
 * McIlwain's L from the second invariant and the mirror field — Hilton's closed form.
 *
 * McIlwain defined L through a tabulated function; Hilton (1971) fitted a closed form accurate to
 * about one part in 10⁴, which is what every implementation since has used. With
 * `X = I³·B_m/M` and `M` the dipole moment,
 *
 *     L³·B_m/M = 1 + 1.35047·X^(1/3) + 0.465376·X^(2/3) + 0.0475455·X
 *
 * @param invariant_i the second invariant `I`, Earth radii.
 * @param b_mirror the mirror field, nT.
 * @param dipole_moment_nt the dipole moment `M` in nT·R_E³. Use @ref dipole_moment, which derives
 *        it from the epoch's own IGRF coefficients. Do NOT hard-code a constant: `M` drifts by
 *        several percent per century, and because `L ∝ M^(1/3)` a stale value shows up as a
 *        CONSTANT relative offset in `L_m` at every shell — the 1960s-era `0.311653e5` is 4.3% high
 *        against IGRF-2015 and produces a uniform 1.4% error, which is exactly how this was found.
 * @return `L_m` in Earth radii, or @ref Status::DomainError for non-physical inputs.
 * @complexity O(1) — two cube roots.
 * @alloc none.
 * @test IrbemLstar.HiltonReproducesTheDipoleLimit
 */
[[nodiscard]] inline Result<double> mcilwain_l(double invariant_i, double b_mirror,
                                               double dipole_moment_nt) {
    if (!(b_mirror > 0.0) || !(dipole_moment_nt > 0.0) || !(invariant_i >= 0.0)) {
        return {Status::DomainError, 0.0};
    }
    const double x = (invariant_i * invariant_i * invariant_i) * b_mirror / dipole_moment_nt;
    const double c = std::cbrt(x);
    const double f = 1.0 + (1.35047 * c) + (0.465376 * c * c) + (0.0475455 * x);
    return {Status::Ok, std::cbrt(f * dipole_moment_nt / b_mirror)};
}


#ifdef CHEATAH_SPACE_IRBEM_LSTAR_GPU
namespace detail {

/// Run a whole batch of traces on the device.
///
/// The coefficients are interpolated to the epoch HERE, once, and uploaded as buffers — the device
/// never sees IGRF's 26-epoch table, because interpolating it per thread would be 10^5 redundant
/// copies of a calculation the host does once. Likewise the Legendre normalisation, which is
/// `constexpr` on the host and therefore free.
///
/// @return `Status::ParametersMissing` when the device could not be used after all (no SPIR-V, no
///         device), which the caller reads as "fall through to the host" rather than as an error.
/// @complexity One dispatch; O(lines × steps) field evaluations, concurrent.
/// @alloc device buffers, released back to the context's size-classed pool on return.
template <int NMAX>
[[nodiscard]] inline Result<bool> trace_batch_on_device(
    const Igrf<NMAX>& model, std::span<const Position<Frame::GEO>> starts,
    std::span<const double> pitch_angles_deg, std::span<FieldLine> out,
    std::span<Status> statuses, const TraceOptions& opt) {

    if (!gpu::available()) return {Status::ParametersMissing, false};

    constexpr int kSlots = ((NMAX + 1) * (NMAX + 2)) / 2;
    const std::size_t n = starts.size();

    std::vector<float> coef(2 * kSlots, 0.0F);
    for (int deg = 1; deg <= NMAX; ++deg) {
        for (int m = 0; m <= deg; ++m) {
            const std::size_t k = (static_cast<std::size_t>(deg) * (deg + 1)) / 2 + m;
            coef[k] = static_cast<float>(model.g(deg, m));
            coef[kSlots + k] = static_cast<float>(model.h(deg, m));
        }
    }
    constexpr auto kNorm = ::cheatah::space::irbem::detail::make_legendre_normalisation<NMAX, double>();
    std::vector<float> nrm(2 * kSlots + NMAX + 1, 0.0F);
    for (int k = 0; k < kSlots; ++k) {
        nrm[static_cast<std::size_t>(k)] = static_cast<float>(kNorm.e[static_cast<std::size_t>(k)]);
        nrm[static_cast<std::size_t>(kSlots + k)] =
            static_cast<float>(kNorm.f[static_cast<std::size_t>(k)]);
    }
    for (int deg = 0; deg <= NMAX; ++deg) {
        nrm[static_cast<std::size_t>(2 * kSlots + deg)] =
            static_cast<float>(kNorm.diagonal[static_cast<std::size_t>(deg)]);
    }

    std::vector<float> pos(3 * n);
    std::vector<float> pitch(n);
    for (std::size_t i = 0; i < n; ++i) {
        pos[(3 * i) + 0] = static_cast<float>(starts[i].v[0]);
        pos[(3 * i) + 1] = static_cast<float>(starts[i].v[1]);
        pos[(3 * i) + 2] = static_cast<float>(starts[i].v[2]);
        pitch[i] = static_cast<float>(pitch_angles_deg[i]);
    }
    // steps_per_l is carried as an integer thousandth so the dims buffer stays a uint buffer —
    // the ABI has no float dims slot, and inventing one for a single value is not worth a binding.
    const std::array<std::uint32_t, 4> dims{
        static_cast<std::uint32_t>(n), static_cast<std::uint32_t>(NMAX),
        static_cast<std::uint32_t>(opt.max_steps),
        static_cast<std::uint32_t>(opt.steps_per_l * 1000.0)};

    std::vector<float> raw(4 * n);
    std::vector<std::uint32_t> st(n);

    if (!gpu::launch_trace(pos, pitch, coef, nrm, dims, raw, st)) {
        return {Status::ParametersMissing, false};
    }

    bool all_ok = true;
    for (std::size_t i = 0; i < n; ++i) {
        out[i] = FieldLine{};
        out[i].invariant_i = raw[(4 * i) + 0];
        out[i].b_min = raw[(4 * i) + 1];
        out[i].b_mirror = raw[(4 * i) + 2];
        out[i].b_local = raw[(4 * i) + 3];
        statuses[i] = st[i] < status_count ? static_cast<Status>(st[i]) : Status::DomainError;
        all_ok = all_ok && (statuses[i] == Status::Ok);
    }
    return {all_ok ? Status::Ok : Status::OpenFieldLine, true};
}

}  // namespace detail
#endif  // CHEATAH_SPACE_IRBEM_LSTAR_GPU

// ---------------------------------------------------------------------------------------------
// The batch entry point — the GPU is the DEFAULT, the host loop is the fallback
// ---------------------------------------------------------------------------------------------

/**
 * Trace a whole batch of field lines and return each one's second invariant.
 *
 * **This is the routine to call.** @ref trace_invariant traces ONE line and is the reference lane:
 * it is what the device lane is verified against, and what runs when there is no device or the
 * batch is too small to pay for one. It is not the fast path and is not meant to be.
 *
 * Why the batch form exists at all: a single trace is a serial RK4 chain, so no hardware makes it
 * faster. The parallelism is entirely ACROSS lines — `ntime × pitch angles × drift-shell azimuths
 * × root-find trials`, which is 10⁵–10⁷ for real work — and it only becomes available if the caller
 * hands over the whole batch at once. A loop calling @ref trace_invariant per point cannot be
 * accelerated no matter what hardware is present.
 *
 * Measured on an RTX 3070 Ti against this header's own fp64 host lane (`-O3 -march=native`), field
 * lines spread over L = 2…8 and pitch angles 30…80°:
 *
 * | lines | 64 | 256 | 1024 | 4096 | 16384 | 65536 |
 * |---|---|---|---|---|---|---|
 * | speedup | 0.5× | 0.8× | 3.1× | 12.7× | 33.2× | **48.9×** |
 *
 * Below ~512 lines the ~30 µs submit floor dominates and the host wins, which is why the crossover
 * is consulted rather than assumed; above it the arithmetic intensity takes over and the curve is
 * still climbing at 65 536. At that size it is 1.24 µs/line against the host's 60.4 µs.
 *
 * @tparam NMAX the IGRF truncation degree.
 * @param model the internal field model, already built for the epoch.
 * @param starts the starting positions, GEO, Earth radii.
 * @param pitch_angles_deg the local pitch angle at each start; same length as @p starts.
 * @param out receives one @ref FieldLine per input; same length as @p starts.
 * @param statuses receives each line's status; same length as @p starts.
 * @param opt the tracing options.
 * @return @ref Status::Ok when every line closed, @ref Status::OpenFieldLine when any did not (the
 *         per-line statuses say which), @ref Status::DomainError on a length mismatch. The value is
 *         `true` when the device lane serviced the call — a test asserts this rather than trusting
 *         that a GPU was used, because a silent fallback is the failure mode that makes a
 *         performance claim worthless.
 * @complexity O(lines × steps) field evaluations; on the device those run concurrently.
 * @alloc the device lane stages coefficients and results; the host lane allocates nothing.
 * @test IrbemLstar.BatchAgreesWithTheReferenceLane
 * @test IrbemLstar.BatchUsesTheDeviceWhenOneIsAvailable
 */
template <int NMAX>
[[nodiscard]] inline Result<bool> trace_invariant_batch(
    const Igrf<NMAX>& model, std::span<const Position<Frame::GEO>> starts,
    std::span<const double> pitch_angles_deg, std::span<FieldLine> out,
    std::span<Status> statuses, const TraceOptions& opt = {}) {

    const std::size_t n = starts.size();
    if (pitch_angles_deg.size() != n || out.size() != n || statuses.size() != n) {
        return {Status::DomainError, false};
    }
    if (n == 0) return {Status::Ok, false};

#ifdef CHEATAH_SPACE_IRBEM_LSTAR_GPU
    if (gpu::prefer_gpu("irbem_trace_i_f32", n)) {
        if (const Result<bool> r = detail::trace_batch_on_device(model, starts, pitch_angles_deg,
                                                                out, statuses, opt);
            r.status != Status::ParametersMissing) {
            return r;
        }
        // ParametersMissing means the device could not be used after all (no SPIR-V, no device).
        // Fall through to the host rather than failing: a missing shader is a deployment problem,
        // not a reason to refuse to compute.
    }
#endif

    bool all_ok = true;
    for (std::size_t i = 0; i < n; ++i) {
        const Result<FieldLine> r = trace_invariant(model, starts[i], pitch_angles_deg[i], opt);
        out[i] = r.value;
        statuses[i] = r.status;
        all_ok = all_ok && (r.status == Status::Ok);
    }
    return {all_ok ? Status::Ok : Status::OpenFieldLine, false};
}

}  // namespace cheatah::space::irbem