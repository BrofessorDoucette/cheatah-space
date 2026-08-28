#pragma once

/**
 * @file context.hpp
 * @brief space.irbem — the immutable per-epoch evaluation state, in place of mutable globals.
 *
 * IRBEM keeps the state of an evaluation in roughly twenty mutable Fortran `COMMON` blocks: the
 * epoch, the dipole tilt, the frame rotations for that instant, the solar-wind drivers, and — the
 * expensive one — a warm start in which each time point seeds its drift-shell search from the
 * previous point's converged answer. That warm start is not an implementation detail. It is why the
 * library cannot be threaded and cannot be offloaded: point `k`'s answer depends on point `k-1`
 * having run first, in the same process, through memory neither of them names in its signature.
 *
 * @ref FieldContext removes the possibility. It is a value: built once per timestamp, validated at
 * construction, `const` thereafter, with no setters, no statics, and no seat for a previous point's
 * result. Every field evaluation for that timestamp takes it by `const&`, so N threads or N GPU
 * invocations read the same bytes and cannot influence each other. Answers stop depending on
 * evaluation order — which also makes them reproducible, and makes a disagreement with the oracle a
 * real disagreement rather than a difference in what ran first.
 *
 * It is trivially copyable and about a kilobyte, so the same object that a CPU worker holds by
 * reference is a `memcpy` into a uniform block — Vulkan guarantees at least 16 KiB of
 * `maxUniformBufferRange`, so it fits everywhere with room to spare. Nothing here allocates,
 * throws, or touches the clock.
 *
 * The layout is split hot/cold on purpose. @ref HotState is what an inner-loop field call touches
 * on every one of the ~10^5 evaluations a single L* costs: the tilt and its precomputed sine and
 * cosine, the five drivers every published external model reads, and the GEO<->GSM rotation pair.
 * Those are the first eight doubles — exactly one 64-byte cache line — followed by the two
 * matrices. @ref ColdState carries the epoch bookkeeping, the full 25-slot driver vector and the
 * rotations to the remaining frames, all of which a trace reads at most once.
 *
 * Precomputing the trigonometry here is the whole point of a per-epoch object: the transcendentals
 * that a frame rotation needs depend on the timestamp and not on the position, so they are paid
 * once per timestamp instead of once per field call.
 *
 * @note The rotation matrices are an INPUT to this file, not an output of it: the coordinate
 *       transform module computes them from the epoch and hands them over. See @ref RotationTable
 *       for the exact contract.
 */

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <numbers>
#include <optional>
#include <span>
#include <string_view>
#include <type_traits>

#include "frames.hpp"

namespace cheatah::space::irbem {

// ---- the driver (`maginput`) vector -------------------------------------------------------------

/**
 * A solar-wind / geomagnetic driver — an index into the `maginput` vector.
 *
 * The names, order and units are IRBEM's `maginput` table (IRBEM `docs/source/api/`
 * `general_information.rst`, "Magnetic field inputs"), because that is the vector every caller
 * already has and every external field model is parameterized by. IRBEM numbers the table from 1;
 * these enumerators are **zero-based**, so an enumerator is the C array subscript directly and the
 * off-by-one lives in exactly one place — this comment.
 *
 * Slots 18..25 of the table are reserved and have no enumerator; they are still carried, still
 * validated, and reachable through @ref FieldContext::drivers.
 */
enum class Driver : std::uint8_t {
    Kp = 0,  ///< Kp as OMNI2 stores it: **Kp x 10**, an integer-valued double in 0..90.
    Dst,     ///< Disturbance storm-time index, nT.
    Dsw,     ///< Solar-wind proton density, cm^-3.
    Vsw,     ///< Solar-wind bulk speed, km/s.
    Pdyn,    ///< Solar-wind dynamic pressure, nPa.
    ByIMF,   ///< Interplanetary magnetic field, GSM y component, nT.
    BzIMF,   ///< Interplanetary magnetic field, GSM z component, nT.
    G1,      ///< Tsyganenko 2001 driver G1 — a one-hour average of the IMF/velocity coupling.
    G2,      ///< Tsyganenko 2001 driver G2.
    G3,      ///< Tsyganenko 2001 driver G3.
    W1,      ///< Tsyganenko & Sitnov 2005 (doi:10.1029/2004JA010798) driver W1.
    W2,      ///< Tsyganenko & Sitnov 2005 driver W2.
    W3,      ///< Tsyganenko & Sitnov 2005 driver W3.
    W4,      ///< Tsyganenko & Sitnov 2005 driver W4.
    W5,      ///< Tsyganenko & Sitnov 2005 driver W5.
    W6,      ///< Tsyganenko & Sitnov 2005 driver W6.
    AL,      ///< Auroral electrojet lower index, nT.
};

/// The length of the driver vector — IRBEM's `maginput` is 25 doubles wide, reserved slots included.
inline constexpr std::size_t driver_count = 25;

/// How many of the 25 slots have a published meaning; the rest are reserved for future use.
inline constexpr std::size_t named_driver_count = 17;

/// The driver vector itself. Fixed size, inline storage, trivially copyable — it rides inside
/// @ref FieldContext and therefore into a GPU uniform block.
using DriverSet = std::array<double, driver_count>;

/**
 * The driver's short name, as the `maginput` table and the model papers spell it.
 * @param d the driver.
 * @return a static string such as `"Pdyn"`; `"?"` for a value outside the enumerator list, which
 *         the reserved slots 18..25 have no enumerator for.
 * @complexity O(1).
 * @alloc none.
 */
constexpr std::string_view name_of(Driver d) {
    switch (d) {
        case Driver::Kp:
            return "Kp";
        case Driver::Dst:
            return "Dst";
        case Driver::Dsw:
            return "Dsw";
        case Driver::Vsw:
            return "Vsw";
        case Driver::Pdyn:
            return "Pdyn";
        case Driver::ByIMF:
            return "ByIMF";
        case Driver::BzIMF:
            return "BzIMF";
        case Driver::G1:
            return "G1";
        case Driver::G2:
            return "G2";
        case Driver::G3:
            return "G3";
        case Driver::W1:
            return "W1";
        case Driver::W2:
            return "W2";
        case Driver::W3:
            return "W3";
        case Driver::W4:
            return "W4";
        case Driver::W5:
            return "W5";
        case Driver::W6:
            return "W6";
        case Driver::AL:
            return "AL";
    }
    return "?";  // a reserved slot, or a value outside the enumerator list; keeps this total
}

// ---- the rotation table -------------------------------------------------------------------------

/// How many frames are Cartesian, and therefore have a rotation: every @ref Frame but GDZ, SPH and
/// RLL, whose components are an angle pair over a radius and are not related to anything by a
/// rotation at all.
inline constexpr std::size_t cartesian_frame_count = frame_count - 3;

/**
 * The slot @p f occupies in a @ref RotationTable.
 *
 * Computed rather than tabulated: @ref Frame's values run GDZ=0, GEO..MAG=1..6, SPH=7, RLL=8,
 * HEE..HEEQ=9..11, so removing the three angular frames is two contiguous shifts. A switch would
 * need a default arm that no correct call can reach, and unreachable code is code that cannot be
 * tested.
 *
 * @param f a Cartesian frame; the result is meaningless (and, for GDZ, wraps) otherwise, which is
 *          why every caller in this header is constrained on @ref CartesianFrame.
 * @return the slot, `0..8`, with GEO at 0.
 * @complexity O(1).
 * @alloc none.
 */
constexpr std::size_t cartesian_slot(Frame f) {
    const auto raw = static_cast<std::size_t>(f);
    return raw <= static_cast<std::size_t>(Frame::MAG) ? raw - 1 : raw - 3;
}

/**
 * The per-epoch rotations, hub-and-spoke: element `cartesian_slot(F)` is the matrix `M_F` for which
 * `v_GEO = M_F * v_F`. GEO's own slot is the identity.
 *
 * **This is the seam with the coordinate-transform module.** That module owes exactly one function,
 * `RotationTable rotations_for(const Epoch&)`, plus the tilt angle its GEO->GSM construction already
 * produces as a by-product; @ref make_field_context takes both and its signature does not change
 * when that function lands. Until then the table is a caller-supplied parameter, which is also what
 * lets this file be tested against exact synthetic rotations instead of against transcendentals.
 *
 * Storing nine matrices rather than all eighty-one ordered pairs is not only a size decision: GEO is
 * the hub because the internal (IGRF) field is evaluated in GEO, so the two hops an inner loop
 * actually makes — GSM->GEO and back — are stored outright, and only the rare cross pairs
 * (say GSE->SM) cost a matrix product.
 */
using RotationTable = std::array<fixarray::mat3d, cartesian_frame_count>;

// ---- the epoch ----------------------------------------------------------------------------------

/// The earliest epoch accepted: 1900 is the first IAGA DGRF epoch, so nothing before it can be
/// evaluated by any internal field model this library will carry.
inline constexpr std::int32_t epoch_min_year = 1900;

/// The latest epoch accepted. A deliberate sanity bound, not a model-validity bound — a model
/// still reports for itself when an epoch is outside its coefficients' span.
inline constexpr std::int32_t epoch_max_year = 2100;

/// The largest accepted seconds-of-day. 86400 is the `24:00:00` spelling of midnight and 86401
/// admits a positive UTC leap second; both are things a real ephemeris contains.
inline constexpr double max_seconds_of_day = 86401.0;

/**
 * The instant an evaluation is for.
 *
 * Both spellings are carried on purpose. `year`/`day_of_year`/`seconds_ut` is what IRBEM's C API
 * takes and what a CDF ephemeris stores; `decimal_year` is what IGRF's coefficient interpolation
 * consumes, and deriving it per point would be arithmetic repeated ~10^5 times for one L*. They are
 * checked against each other at construction rather than trusted (see @ref ContextError).
 */
struct Epoch {
    /// The epoch as a fractional year, e.g. 2015.5 — the argument IGRF's interpolation takes.
    double decimal_year;
    /// Seconds of day, UT, in `[0, 86401]`.
    double seconds_ut;
    /// Calendar year, in `[epoch_min_year, epoch_max_year]`.
    std::int32_t year;
    /// Day of year, 1-based, in `[1, 366]`.
    std::int32_t day_of_year;
};

// ---- the hot / cold split -----------------------------------------------------------------------

/**
 * What the inner loop touches — the field call and the GEO<->GSM hop, on every one of the ~10^5
 * evaluations a single L* costs.
 *
 * The first eight doubles are exactly one 64-byte cache line, and they are the ones a Tsyganenko
 * external field call needs: the tilt with its precomputed trigonometry, and the five drivers the
 * published models read (Kp for Mead & Fairfield 1975, Tsyganenko 1987/1989 and Mead-Tsyganenko;
 * Dst, Pdyn, By and Bz for Tsyganenko 1996 onwards). The two matrices follow.
 *
 * @ref geo_to_gsm and @ref gsm_to_geo duplicate what @ref ColdState::to_geo already holds for GSM.
 * That is deliberate: the duplicate keeps the inner loop off the cold block entirely and spares it
 * a transpose per call. A test pins the two copies to each other so the duplicate cannot go stale.
 */
struct HotState {
    /// The geodipole tilt angle psi, radians — the angle between the dipole axis and the GSM z axis.
    double tilt_rad;
    /// `sin(tilt_rad)`, computed once per epoch so no field call computes it.
    double sin_tilt;
    /// `cos(tilt_rad)`, likewise.
    double cos_tilt;
    /// @ref Driver::Kp, copied out of the driver vector (OMNI2 scaling: Kp x 10).
    double kp;
    /// @ref Driver::Dst, nT.
    double dst;
    /// @ref Driver::Pdyn, nPa.
    double pdyn;
    /// @ref Driver::ByIMF, nT.
    double by_imf;
    /// @ref Driver::BzIMF, nT.
    double bz_imf;
    /// Takes GEO components to GSM.
    fixarray::mat3d geo_to_gsm;
    /// Takes GSM components to GEO — the transpose of @ref geo_to_gsm, since a rotation is
    /// orthogonal and construction has already checked that this one is.
    fixarray::mat3d gsm_to_geo;
};

/**
 * What an evaluation reads at most once: the epoch bookkeeping, the whole driver vector, and the
 * rotations to the frames an inner loop does not visit.
 */
struct ColdState {
    /// The instant this context is for.
    Epoch epoch;
    /// All 25 `maginput` slots, in IRBEM's order, reserved entries included.
    DriverSet drivers;
    /// The rotations, indexed by @ref cartesian_slot; see @ref RotationTable for the convention.
    RotationTable to_geo;
};

// ---- failure ------------------------------------------------------------------------------------

/**
 * Why a set of inputs cannot form a @ref FieldContext.
 *
 * IRBEM's `baddata = -1e31` sentinel is kept at the outer API boundary for compatibility and
 * nowhere else: inside the library a failure carries a reason, so a wrong answer cannot be mistaken
 * for a number and a rejected build says which input was wrong.
 */
enum class ContextError : std::uint8_t {
    None = 0,                 ///< No failure — the inputs are usable.
    YearOutOfRange,           ///< @ref Epoch::year outside `[epoch_min_year, epoch_max_year]`.
    DayOfYearOutOfRange,      ///< @ref Epoch::day_of_year outside `[1, 366]`.
    SecondsOfDayOutOfRange,   ///< @ref Epoch::seconds_ut not finite, or outside `[0, 86401]`.
    DecimalYearInconsistent,  ///< @ref Epoch::decimal_year not finite, or not inside its own year.
    TiltNotFinite,            ///< The dipole tilt is NaN or infinite.
    TiltOutOfRange,           ///< `|tilt| > pi/2`, which no angle-to-an-axis can be.
    DriverNotFinite,          ///< Some `maginput` slot is NaN or infinite.
    RotationNotFinite,        ///< Some rotation matrix has a NaN or infinite element.
    RotationNotOrthogonal,    ///< Some rotation matrix fails `M^T M = I` to within tolerance.
    RotationImproper,         ///< Some rotation matrix is orthogonal but has a negative determinant.
};

/**
 * A human-readable reason, for a log line or a test failure message.
 * @param e the failure.
 * @return a static string such as `"driver not finite"`; `"?"` for a value outside the enumerator
 *         list.
 * @complexity O(1).
 * @alloc none.
 */
constexpr std::string_view describe(ContextError e) {
    switch (e) {
        case ContextError::None:
            return "ok";
        case ContextError::YearOutOfRange:
            return "epoch year out of range";
        case ContextError::DayOfYearOutOfRange:
            return "epoch day of year out of range";
        case ContextError::SecondsOfDayOutOfRange:
            return "epoch seconds of day out of range";
        case ContextError::DecimalYearInconsistent:
            return "decimal year disagrees with the calendar year";
        case ContextError::TiltNotFinite:
            return "dipole tilt not finite";
        case ContextError::TiltOutOfRange:
            return "dipole tilt out of range";
        case ContextError::DriverNotFinite:
            return "driver not finite";
        case ContextError::RotationNotFinite:
            return "rotation matrix not finite";
        case ContextError::RotationNotOrthogonal:
            return "rotation matrix not orthogonal";
        case ContextError::RotationImproper:
            return "rotation matrix is improper (a reflection)";
    }
    return "?";  // unreachable for a valid enumerator; keeps the function total
}

/// How far `M^T M` may stray from the identity before a matrix is refused as a rotation. Slack
/// enough that a matrix assembled from sines and cosines (error a few ulp, ~1e-16) never trips it,
/// tight enough that a scaled, sheared or garbage matrix always does.
inline constexpr double rotation_orthogonality_tolerance = 1e-9;

/// The largest possible dipole tilt. The tilt is the angle between the dipole axis and an axis, so
/// it is bounded by pi/2 by definition; the physical excursion over a year is about +-35 degrees.
inline constexpr double max_tilt_rad = std::numbers::pi / 2.0;

// ---- validation ---------------------------------------------------------------------------------

namespace detail {

/**
 * Whether every element of @p xs is finite.
 * @param xs the values; a matrix passes `std::span(m.data(), mat3d::size)`, a driver set passes
 *           itself.
 * @return false as soon as one element is NaN or infinite.
 * @complexity O(n) in the span length, and it short-circuits.
 * @alloc none.
 */
inline bool all_finite(std::span<const double> xs) {
    return std::ranges::all_of(xs, [](double x) { return std::isfinite(x); });
}

/**
 * Whether @p m really is a rotation, and if not, which way it fails.
 *
 * The finiteness scan must come first and cannot be folded into the tolerance test: every
 * comparison against a NaN is false, so `|NaN - 1| > tol` does not fire and a NaN matrix would be
 * accepted as orthogonal. Properness is a separate question from orthogonality — `M^T M = I` is
 * equally true of a reflection, which would mirror every trace and still produce plausible-looking
 * numbers.
 *
 * @param m the candidate rotation.
 * @return ContextError::None when @p m is a proper rotation, otherwise the specific defect.
 * @complexity O(1) — one 3x3 product, nine comparisons and one determinant.
 * @alloc none.
 */
inline ContextError rotation_defect(const fixarray::mat3d& m) {
    if (!all_finite(std::span<const double>(m.data(), fixarray::mat3d::size))) {
        return ContextError::RotationNotFinite;
    }
    const fixarray::mat3d gram = fixarray::transpose(m) * m;
    const fixarray::mat3d unit = fixarray::mat3d::identity();
    for (std::size_t i = 0; i < fixarray::mat3d::size; ++i) {
        if (std::fabs(gram.data()[i] - unit.data()[i]) > rotation_orthogonality_tolerance) {
            return ContextError::RotationNotOrthogonal;
        }
    }
    if (fixarray::determinant(m) <= 0.0) return ContextError::RotationImproper;
    return ContextError::None;
}

/**
 * Whether @p e is a usable epoch, and if not, which field is wrong.
 *
 * The consistency check between @ref Epoch::decimal_year and @ref Epoch::year is the one that earns
 * its keep: carrying both spellings is a performance decision, and a caller who fills one from a
 * different timestamp than the other produces field coefficients for the wrong year with no other
 * symptom.
 *
 * @param e the epoch.
 * @return ContextError::None when @p e is usable, otherwise the specific defect.
 * @complexity O(1).
 * @alloc none.
 */
inline ContextError epoch_defect(const Epoch& e) {
    if (e.year < epoch_min_year || e.year > epoch_max_year) return ContextError::YearOutOfRange;
    if (e.day_of_year < 1 || e.day_of_year > 366) return ContextError::DayOfYearOutOfRange;
    if (!std::isfinite(e.seconds_ut) || e.seconds_ut < 0.0 || e.seconds_ut > max_seconds_of_day) {
        return ContextError::SecondsOfDayOutOfRange;
    }
    const auto year = static_cast<double>(e.year);
    if (!std::isfinite(e.decimal_year) || e.decimal_year < year || e.decimal_year > year + 1.0) {
        return ContextError::DecimalYearInconsistent;
    }
    return ContextError::None;
}

}  // namespace detail

class ContextResult;

// ---- the context --------------------------------------------------------------------------------

/**
 * Everything a field evaluation needs to know about *when* it is, and nothing about *where*.
 *
 * Immutable by construction: the data are private, every accessor is `const`, there are no setters
 * and the constructor is private, so the only way to obtain one is @ref make_field_context and the
 * only contexts that exist are validated ones. Nothing in it can record what the previous
 * evaluation found, which is precisely the affordance IRBEM's `COMMON` blocks provide and precisely
 * why they cannot be threaded.
 *
 * Immutability is spelled "private data plus const accessors" rather than `const` data members on
 * purpose: `const` members would delete copy assignment and make the type useless in an ephemeris
 * buffer, while preventing nothing that is not already prevented.
 *
 * Trivially copyable and standard layout, so `memcpy` into a mapped uniform buffer is a valid way
 * to hand it to a GPU, and passing it by value between CPU workers costs a copy of about a
 * kilobyte and no indirection.
 */
class FieldContext {
  public:
    /**
     * The inner loop's block — tilt, precomputed trigonometry, the five universal drivers, and the
     * GEO<->GSM rotation pair.
     * @return a reference into this context; it lives exactly as long as the context does.
     * @complexity O(1).
     * @alloc none.
     */
    [[nodiscard]] const HotState& hot() const noexcept { return hot_; }

    /**
     * The instant this context is for.
     * @return the validated epoch.
     * @complexity O(1).
     * @alloc none.
     */
    [[nodiscard]] const Epoch& epoch() const noexcept { return cold_.epoch; }

    /**
     * The whole `maginput` vector, in IRBEM's order, reserved slots included.
     * @return all 25 drivers.
     * @complexity O(1).
     * @alloc none.
     */
    [[nodiscard]] const DriverSet& drivers() const noexcept { return cold_.drivers; }

    /**
     * One named driver.
     * @param d which driver; @ref Driver's enumerators are the zero-based `maginput` subscripts.
     * @return its value, in the unit @ref Driver documents.
     * @complexity O(1).
     * @alloc none.
     */
    [[nodiscard]] double driver(Driver d) const noexcept {
        return cold_.drivers[static_cast<std::size_t>(d)];
    }

    /**
     * The stored matrix taking @p F's components to GEO — no arithmetic, just the table entry.
     * @tparam F the source frame; Cartesian, since an angular frame is not related to anything by a
     *           rotation.
     * @return a reference to the matrix `M_F` for which `v_GEO = M_F * v_F`.
     * @complexity O(1).
     * @alloc none.
     */
    template <Frame F>
        requires CartesianFrame<F>
    [[nodiscard]] const fixarray::mat3d& rotation_to_geo() const noexcept {
        return cold_.to_geo[cartesian_slot(F)];
    }

    /**
     * The rotation taking @p From's components to @p To's, composed through the GEO hub.
     *
     * The three cases that cost nothing — a frame to itself, anything to GEO, GEO to anything — are
     * chosen at compile time, so only a genuine cross pair such as GSE->SM pays a matrix product.
     * A rotation's inverse is its transpose, which is why the GEO-to-anything case is a transpose
     * and not a solve.
     *
     * @tparam From the source frame; Cartesian.
     * @tparam To the destination frame; Cartesian.
     * @return the matrix `R` for which `v_To = R * v_From`.
     * @complexity O(1); at most one 3x3 product.
     * @alloc none.
     */
    template <Frame From, Frame To>
        requires CartesianFrame<From> && CartesianFrame<To>
    [[nodiscard]] fixarray::mat3d rotation() const noexcept {
        if constexpr (From == To) {
            return fixarray::mat3d::identity();
        } else if constexpr (To == Frame::GEO) {
            return cold_.to_geo[cartesian_slot(From)];
        } else if constexpr (From == Frame::GEO) {
            return fixarray::transpose(cold_.to_geo[cartesian_slot(To)]);
        } else {
            return fixarray::transpose(cold_.to_geo[cartesian_slot(To)]) *
                   cold_.to_geo[cartesian_slot(From)];
        }
    }

    /**
     * Re-express a position in another frame.
     *
     * This applies a rotation and only a rotation, so it is defined between Cartesian frames alone.
     * Converting to or from GDZ, SPH or RLL is not a rotation — it involves the reference ellipsoid
     * or a Cartesian-to-angular change of variables — and belongs to the coordinate transform
     * module, not here.
     *
     * @tparam To the destination frame; Cartesian.
     * @tparam From the source frame; Cartesian, deduced from @p p.
     * @param p the position, in Earth radii.
     * @return the same point, tagged with and expressed in @p To.
     * @complexity O(1).
     * @alloc none.
     */
    template <Frame To, Frame From>
        requires CartesianFrame<From> && CartesianFrame<To>
    [[nodiscard]] Position<To> rotate(Position<From> p) const noexcept {
        return Position<To>{rotation<From, To>() * p.v};
    }

    /**
     * Re-express a field vector in another frame.
     *
     * Separate from the position overload because the types are: a field in nanotesla and a position
     * in Earth radii transform by the same matrix but must never be interchangeable.
     *
     * @tparam To the destination frame; Cartesian.
     * @tparam From the source frame; Cartesian, deduced from @p b.
     * @param b the field vector, in nanotesla.
     * @return the same field, tagged with and expressed in @p To.
     * @complexity O(1).
     * @alloc none.
     */
    template <Frame To, Frame From>
        requires CartesianFrame<From> && CartesianFrame<To>
    [[nodiscard]] FieldVector<To> rotate(FieldVector<From> b) const noexcept {
        return FieldVector<To>{rotation<From, To>() * b.v};
    }

    /// The only function permitted to build one, so every context that exists has been validated.
    friend ContextResult make_field_context(const Epoch& epoch, double tilt_rad,
                                            const RotationTable& to_geo, const DriverSet& drivers);

  private:
    /**
     * Assemble the state. Private: the inputs are assumed already validated by
     * @ref make_field_context, which is the only caller.
     * @param epoch the instant. @param tilt_rad the geodipole tilt, radians.
     * @param to_geo the rotations, per @ref RotationTable. @param drivers all 25 `maginput` slots.
     * @complexity O(1) — two transcendentals and a fixed number of copies, once per epoch.
     * @alloc none.
     */
    FieldContext(const Epoch& epoch, double tilt_rad, const RotationTable& to_geo,
                 const DriverSet& drivers)
        : hot_{tilt_rad,
               std::sin(tilt_rad),
               std::cos(tilt_rad),
               drivers[static_cast<std::size_t>(Driver::Kp)],
               drivers[static_cast<std::size_t>(Driver::Dst)],
               drivers[static_cast<std::size_t>(Driver::Pdyn)],
               drivers[static_cast<std::size_t>(Driver::ByIMF)],
               drivers[static_cast<std::size_t>(Driver::BzIMF)],
               fixarray::transpose(to_geo[cartesian_slot(Frame::GSM)]),
               to_geo[cartesian_slot(Frame::GSM)]},
          cold_{epoch, drivers, to_geo} {}

    /// The inner loop's cache line and rotation pair.
    HotState hot_;
    /// Everything read at most once per evaluation.
    ColdState cold_;
};

static_assert(std::is_trivially_copyable_v<FieldContext>,
              "FieldContext is memcpy'd into a GPU uniform block; it must stay a plain value");
static_assert(std::is_standard_layout_v<FieldContext>,
              "a GPU-visible struct needs a layout the host and device agree on");
static_assert(sizeof(HotState) == 8 * sizeof(double) + 2 * sizeof(fixarray::mat3d),
              "the hot block must stay exactly one cache line of scalars plus the two matrices");

/**
 * The outcome of a context build: a context, or the named reason there is none.
 *
 * There is deliberately no default state and no public constructor. A struct of
 * `{optional<FieldContext>, ContextError}` would have a value-initialized form carrying no context
 * and no error — a success with nothing in it, which is exactly the sort of half-valid state this
 * file exists to make unrepresentable. Both constructors here are private and each establishes the
 * invariant that a context is present if and only if @ref error is ContextError::None, so
 * @ref make_field_context is the only source of either half.
 */
class ContextResult {
  public:
    /// Whether a context was built.
    /// @return true exactly when @ref error is ContextError::None.
    /// @complexity O(1).
    /// @alloc none.
    [[nodiscard]] bool has_value() const noexcept { return context_.has_value(); }

    /// Why no context was built.
    /// @return the failure, or ContextError::None when one was built.
    /// @complexity O(1).
    /// @alloc none.
    [[nodiscard]] ContextError error() const noexcept { return error_; }

    /// The context that was built.
    /// @return a reference to it; it lives as long as this result does.
    /// @throws std::bad_optional_access when the build failed — consult @ref error first.
    /// @complexity O(1).
    /// @alloc none.
    [[nodiscard]] const FieldContext& value() const { return context_.value(); }

    /// The only function permitted to build one.
    friend ContextResult make_field_context(const Epoch& epoch, double tilt_rad,
                                            const RotationTable& to_geo, const DriverSet& drivers);

  private:
    /// The success case.
    /// @param built the validated context.
    /// @complexity O(1). @alloc none.
    explicit ContextResult(const FieldContext& built) : context_(built), error_(ContextError::None) {}

    /// The failure case; @ref context_ stays disengaged.
    /// @param why the failure, never ContextError::None.
    /// @complexity O(1). @alloc none.
    explicit ContextResult(ContextError why) : error_(why) {}

    /// The context, engaged only on success.
    std::optional<FieldContext> context_;
    /// Why there is none; ContextError::None on success.
    ContextError error_;
};

/**
 * Validate the inputs and, if they are usable, build the context for this epoch.
 *
 * The checks run in a fixed order — epoch, tilt, drivers, rotations — so a caller with several
 * problems is told about them in that sequence rather than in an order that depends on the build.
 *
 * Every driver is checked for finiteness and **none is range-checked**. That is deliberate: a caller
 * legitimately leaves the slots its chosen external model ignores filled with a placeholder, so
 * refusing a negative solar-wind speed here would reject a perfectly valid internal-field-only
 * evaluation. Whether a driver is inside the range its model was fitted over is that model's
 * question, asked where the driver is actually consumed.
 *
 * @param epoch the instant; see @ref detail::epoch_defect for what makes one usable.
 * @param tilt_rad the geodipole tilt psi, radians, `|psi| <= pi/2`.
 * @param to_geo the rotations for this epoch, per @ref RotationTable — each is checked to be a
 *               finite, orthogonal, proper rotation, which is the guard on the seam with the
 *               coordinate transform module.
 * @param drivers all 25 `maginput` slots; each must be finite.
 * @return the context, or the first defect found.
 * @complexity O(1) — bounded by the 25 driver checks and the nine 3x3 rotation checks, all paid
 *             once per epoch rather than once per field evaluation.
 * @alloc none.
 */
inline ContextResult make_field_context(const Epoch& epoch, double tilt_rad,
                                        const RotationTable& to_geo, const DriverSet& drivers) {
    const ContextError when = detail::epoch_defect(epoch);
    if (when != ContextError::None) return ContextResult{when};
    if (!std::isfinite(tilt_rad)) return ContextResult{ContextError::TiltNotFinite};
    if (std::fabs(tilt_rad) > max_tilt_rad) return ContextResult{ContextError::TiltOutOfRange};
    if (!detail::all_finite(std::span<const double>(drivers))) {
        return ContextResult{ContextError::DriverNotFinite};
    }
    for (const fixarray::mat3d& m : to_geo) {
        const ContextError defect = detail::rotation_defect(m);
        if (defect != ContextError::None) return ContextResult{defect};
    }
    return ContextResult{FieldContext(epoch, tilt_rad, to_geo, drivers)};
}

}  // namespace cheatah::space::irbem
