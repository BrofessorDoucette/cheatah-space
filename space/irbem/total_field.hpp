#pragma once

/**
 * @file total_field.hpp
 * @brief space.irbem — the internal field plus an external model, as one field a tracer can follow.
 *
 * A particle does not experience the internal and external fields separately, and neither should a
 * field-line trace. This is the superposition: IGRF plus a Tsyganenko external model, summed at
 * every point, presented as a single @ref GeoFieldModel so that @ref trace_invariant,
 * @ref make_lstar and everything built on them follow the *total* field without knowing there are
 * two of them.
 *
 * ## Why this type has to exist
 *
 * Until it did, `trace_invariant` took a `const Igrf<NMAX>&` — so nothing but the internal field
 * could ever reach a trace, and three consequences followed that all look unrelated and are not:
 *
 *  - **No storm testing was possible.** Every Tsyganenko model is *parameterized by* geomagnetic
 *    activity, and `maginput` had nowhere to go. A corpus sweeping Kp, Dst, Pdyn and southward Bz
 *    could not change a single number, so the claim that this library reproduces the models under
 *    disturbed conditions was untestable rather than merely untested.
 *  - **`Status::NotConverged` was unreachable.** It fires when a drift shell has a gap, and a pure
 *    internal field is dipole-like everywhere: every shell closes, at every L out to 40 and every
 *    step cap — measured, not assumed. Open field lines need a magnetopause and a tail.
 *  - **No mock could be injected**, so those paths could not even be tested synthetically.
 *
 * ## The frames, which is where this is easy to get wrong
 *
 * IGRF is defined in GEO; every Tsyganenko model is defined in GSM. The rotation between them
 * depends only on the epoch, so it is built ONCE in @ref Rotations and reused for every point —
 * which is also what keeps the transcendentals out of the inner loop. A `FieldVector` transform is
 * a pure rotation with no origin shift, unlike a `Position` transform; the frame-tagged types make
 * that distinction a compile error rather than a plausible wrong answer.
 *
 * ## The device lane
 *
 * @ref trace_invariant_batch below routes a `TotalFieldT89` batch to `irbem_trace_total_f32`,
 * which evaluates BOTH models per RK4 stage on the device — `igrf_eval` and `t89_eval` are the
 * same shared Slang functions the batched field kernels call, so each piece of physics exists on
 * the device exactly once. The host stages three things per batch: the epoch-interpolated IGRF
 * coefficients, the ONE Kp bin's T89 parameter block (a batch shares one epoch and one Kp — a
 * per-thread bin branch would diverge the warp), and the epoch's `gsm_from_geo` rotation. Single
 * traces and small batches stay on the host, below the same measured crossover the internal
 * tracer uses.
 */

#include <cmath>
#include <numbers>

#include "coords_rotations.hpp"
#include "lstar.hpp"
#include <span>
#include <vector>

#include "ext_t89.hpp"
#include "frames.hpp"
#include "igrf.hpp"
#include "status.hpp"

namespace cheatah::space::irbem {

/**
 * IGRF plus Tsyganenko 1989, as a single field.
 *
 * Satisfies @ref GeoFieldModel, so it drops into @ref trace_invariant and everything above it.
 *
 * @tparam NMAX the internal field's truncation degree. 10 reproduces IRBEM's own choice, which is
 *         what the differential tests run through; 13 is IGRF-14's full published degree.
 */
template <int NMAX = 10>
class TotalFieldT89 {
  public:
    /// The internal part's truncation degree — what generic staging and buffer sizing read, on
    /// both this type and @ref Igrf, so `M::degree` means the same thing for either.
    static constexpr int degree = NMAX;

    /**
     * @param internal the internal field, already built for the epoch.
     * @param rotations the epoch's frame rotations — built once, reused for every point.
     * @param kp_times_ten Kp in IRBEM's slot-1 scaling (Kp x 10, nominally 0..90). T89 is
     *        Kp-BINNED, not a continuous function of Kp: seven coefficient sets, so values inside
     *        one bin give identical fields by construction rather than by approximation.
     */
    constexpr TotalFieldT89(const Igrf<NMAX>& internal, const Rotations& rotations,
                            double kp_times_ten)
        : internal_(&internal), rotations_(&rotations), kp_times_ten_(kp_times_ten) {}

    /**
     * The total field at a geographic point.
     *
     * @param p the position, GEO, Earth radii.
     * @return `B_internal + B_external`, in GEO, nT. When the external model refuses the point —
     *         outside its published validity envelope, or a non-finite input — the INTERNAL field
     *         is returned alone rather than a zero or a NaN: the internal field is still the best
     *         available answer there, and a trace that hit a NaN would fail hundreds of steps later
     *         with no indication of where. Callers who need to know whether the external model
     *         contributed ask @ref external_status for the point.
     * @complexity One IGRF evaluation, one T89 evaluation, two 3x3 rotations.
     * @alloc none.
     * @test IrbemTotalField.SuperposesInternalAndExternal
     */
    [[nodiscard]] FieldVector<Frame::GEO> evaluate(const Position<Frame::GEO>& p) const {
        const FieldVector<Frame::GEO> b_int = internal_->evaluate(p);
        const Position<Frame::GSM> p_gsm = transform<Frame::GSM>(p, *rotations_);
        const double tilt_rad = rotations_->dipole_tilt_deg * (std::numbers::pi / 180.0);
        const Result<FieldVector<Frame::GSM>> b_ext = t89_field(p_gsm, tilt_rad, kp_times_ten_);
        if (b_ext.status == Status::DomainError) return b_int;
        // A FieldVector transform is a pure rotation — no origin shift, unlike a Position.
        const FieldVector<Frame::GEO> b_ext_geo = transform<Frame::GEO>(b_ext.value, *rotations_);
        return FieldVector<Frame::GEO>{b_int.v + b_ext_geo.v};
    }

    /**
     * Whether the external model answered at @p p, and if not, why.
     * @param p the position, GEO, Earth radii.
     * @return the external model's status; @ref Status::Ok when it contributed.
     * @complexity One T89 evaluation and one rotation.
     * @alloc none.
     * @test IrbemTotalField.ReportsWhenTheExternalModelDeclines
     */
    [[nodiscard]] Status external_status(const Position<Frame::GEO>& p) const {
        const Position<Frame::GSM> p_gsm = transform<Frame::GSM>(p, *rotations_);
        const double tilt_rad = rotations_->dipole_tilt_deg * (std::numbers::pi / 180.0);
        return t89_field(p_gsm, tilt_rad, kp_times_ten_).status;
    }

    /// The activity level this field was built for, in IRBEM's Kp x 10 scaling.
    /// @return the value passed to the constructor.
    /// @complexity O(1). @alloc none.
    /// @test IrbemTotalField.SuperposesInternalAndExternal
    [[nodiscard]] constexpr double kp_times_ten() const { return kp_times_ten_; }

    /// The epoch's frame rotations — what the device staging and any caller mapping frames needs.
    /// @return the rotations this field was built with.
    /// @complexity O(1). @alloc none.
    /// @test IrbemTotalField.SuperposesInternalAndExternal
    [[nodiscard]] constexpr const Rotations& rotations() const { return *rotations_; }

    /// The internal part's Gauss coefficient `g(n, m)`, in nT.
    ///
    /// A superposition has no spherical-harmonic expansion of its own — the external field is not
    /// current-free, so no scalar potential exists to expand. What a caller asking `g(1, 0)` of a
    /// total field means, in every use this module has (the dipole moment `k0`, the trace step
    /// sizing, the device staging of the internal part), is the INTERNAL field's coefficient, and
    /// that is what this forwards to. Documented here precisely because silently answering a
    /// question the physics cannot pose is how a wrong number acquires authority.
    /// @param n the degree. @param m the order. @return the internal part's coefficient.
    /// @complexity O(1). @alloc none.
    /// @test IrbemTotalField.SuperposesInternalAndExternal
    [[nodiscard]] constexpr double g(int n, int m) const { return internal_->g(n, m); }

    /// The internal part's `h(n, m)`, in nT — see @ref g for why this forwards.
    /// @param n the degree. @param m the order. @return the internal part's coefficient.
    /// @complexity O(1). @alloc none.
    /// @test IrbemTotalField.SuperposesInternalAndExternal
    [[nodiscard]] constexpr double h(int n, int m) const { return internal_->h(n, m); }

    /// The internal field alone — what `dipole_moment` and the device staging need, since those
    /// are questions about the internal field specifically and a superposition cannot answer them.
    /// @return the internal model.
    /// @complexity O(1). @alloc none.
    /// @test IrbemTotalField.SuperposesInternalAndExternal
    [[nodiscard]] constexpr const Igrf<NMAX>& internal() const { return *internal_; }

  private:
    const Igrf<NMAX>* internal_;
    const Rotations* rotations_;
    double kp_times_ten_;
};


/**
 * Trace a batch of field lines through the TOTAL field — the GPU-default entry point.
 *
 * The shape mirrors the internal field's @ref trace_invariant_batch exactly: device above the
 * measured crossover, the fp64 host loop below it or when no device answers, and the returned
 * value says which lane actually served the call — a silent CPU fallback is the failure mode that
 * makes a performance claim worthless, so it is observable by construction.
 *
 * @tparam NMAX the internal field's truncation degree.
 * @param field the superposed model; carries the epoch's rotations and the batch's Kp.
 * @param starts the starting positions, GEO, Earth radii.
 * @param pitch_angles_deg the local pitch angle at each start; same length as @p starts.
 * @param out one @ref FieldLine per input. @param statuses one @ref Status per input.
 * @param opt the tracing options.
 * @return @ref Status::Ok when every line closed; the value is `true` when the DEVICE served the
 *         call.
 * @complexity O(lines x steps) total-field evaluations (~900 flops each); concurrent on device.
 * @alloc staging vectors for the device lane; the host lane allocates nothing.
 * @test IrbemTotalField.BatchAgreesWithTheReferenceLane
 * @test IrbemTotalField.BatchUsesTheDeviceWhenOneIsAvailable
 */
template <int NMAX>
[[nodiscard]] inline Result<bool> trace_invariant_batch(
    const TotalFieldT89<NMAX>& field, std::span<const Position<Frame::GEO>> starts,
    std::span<const double> pitch_angles_deg, std::span<FieldLine> out, std::span<Status> statuses,
    const TraceOptions& opt = {}) {

    const std::size_t n = starts.size();
    if (pitch_angles_deg.size() != n || out.size() != n || statuses.size() != n) {
        return {Status::DomainError, false};
    }
    if (n == 0) return {Status::Ok, false};

#ifdef CHEATAH_SPACE_IRBEM_LSTAR_GPU
    if (gpu::prefer_gpu("irbem_trace_total_f32", n)) {
        // --- stage: IGRF coefficients + normalisation, epoch-interpolated ONCE on the host ------
        constexpr int kSlots = ((NMAX + 1) * (NMAX + 2)) / 2;
        const Igrf<NMAX>& igrf = field.internal();
        std::vector<float> coef(2 * kSlots, 0.0F);
        for (int deg = 1; deg <= NMAX; ++deg) {
            for (int m = 0; m <= deg; ++m) {
                const std::size_t k = (static_cast<std::size_t>(deg) * (deg + 1)) / 2 + m;
                coef[k] = static_cast<float>(igrf.g(deg, m));
                coef[kSlots + k] = static_cast<float>(igrf.h(deg, m));
            }
        }
        constexpr auto kNorm = detail::make_legendre_normalisation<NMAX, double>();
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

        // --- stage: the ext block — one Kp bin's T89 parameters + gsm_from_geo -----------------
        const double tilt_rad = field.rotations().dipole_tilt_deg * (std::numbers::pi / 180.0);
        const std::array<float, t89_param_count> par = t89_param_block(
            static_cast<float>(std::sin(tilt_rad)), static_cast<float>(std::cos(tilt_rad)),
            t89_kp_bin(field.kp_times_ten()));
        std::vector<float> ext(t89_param_count + 9, 0.0F);
        for (std::size_t k = 0; k < t89_param_count; ++k) ext[k] = par[k];
        const fixarray::mat3d r = rotation_matrix<Frame::GSM, Frame::GEO>(field.rotations());
        // fixarray matrices are column-major, which is precisely the layout the kernel reads.
        const double* rd = r.data();
        for (std::size_t k = 0; k < 9; ++k) ext[t89_param_count + k] = static_cast<float>(rd[k]);

        // --- stage: positions and pitches -------------------------------------------------------
        std::vector<float> pos(3 * n);
        std::vector<float> pitch(n);
        for (std::size_t i = 0; i < n; ++i) {
            pos[(3 * i) + 0] = static_cast<float>(starts[i].v[0]);
            pos[(3 * i) + 1] = static_cast<float>(starts[i].v[1]);
            pos[(3 * i) + 2] = static_cast<float>(starts[i].v[2]);
            pitch[i] = static_cast<float>(pitch_angles_deg[i]);
        }
        const std::array<std::uint32_t, 4> dims{
            static_cast<std::uint32_t>(n), static_cast<std::uint32_t>(NMAX),
            static_cast<std::uint32_t>(opt.max_steps),
            static_cast<std::uint32_t>(opt.steps_per_l * 1000.0)};

        std::vector<float> raw(4 * n);
        std::vector<std::uint32_t> st(n);
        if (gpu::launch_trace_total(pos, pitch, coef, nrm, ext, dims, raw, st)) {
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
        // No device / no SPIR-V after all: fall through to the host, which is a deployment
        // situation and not a reason to refuse to compute.
    }
#endif

    bool all_ok = true;
    for (std::size_t i = 0; i < n; ++i) {
        const Result<FieldLine> r = trace_invariant(field, starts[i], pitch_angles_deg[i], opt);
        out[i] = r.value;
        statuses[i] = r.status;
        all_ok = all_ok && (r.status == Status::Ok);
    }
    return {all_ok ? Status::Ok : Status::OpenFieldLine, false};
}

}  // namespace cheatah::space::irbem
