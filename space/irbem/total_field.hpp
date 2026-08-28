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
 * ## What this is NOT
 *
 * This is the **host** superposition. The device lane's kernels upload one model's coefficients and
 * evaluate one model; a GPU trace through the total field needs a kernel that evaluates both, which
 * is a separate piece of work. So a batch traced through a `TotalField` runs on the CPU today, and
 * the GPU speedup quoted elsewhere is for the internal field. Saying otherwise would be the kind of
 * claim that is true of the benchmark and false of the library.
 */

#include <cmath>
#include <numbers>

#include "coords_rotations.hpp"
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

}  // namespace cheatah::space::irbem
