#pragma once

/**
 * @file purr.hpp
 * @brief space.irbem's cheatah-facing surface — the physics, callable from a `.purr` program.
 *
 * The rest of this module is written as a C++20 library: the frame lives in the type
 * (@ref Position "Position<Frame::GEO>"), the truncation degree and precision are template
 * parameters (@ref Igrf "Igrf<13, Exact>"), batches take `std::span`, and @ref Igrf's constructor is
 * private behind a static factory. Every one of those is a correctness win in C++ and a wall in
 * cheatah, because cheatah's interop is namespace-path joining and its type-argument grammar cannot
 * spell a qualified enum template argument — `Position<Frame::GSM>` is unwritable, and so is
 * `Igrf<13>::at`.
 *
 * This header is the boundary @ref frames.hpp already describes: *"IRBEM's C API selects a frame
 * with a runtime `sysaxes` integer. That integer enters the typed world through
 * @ref frame_from_sysaxes at the API boundary and nowhere else."* Here, that boundary is real code.
 * Everything below takes integers, doubles and `ndarray`s and returns the same, so a cheatah
 * program can call it; internally each function is a thin call into the templated core, which is
 * unchanged.
 *
 * The shape deliberately mirrors IRBEM's own C API — `kext`, `sysaxes`, `options(3)`/`options(4)` —
 * because that is what a researcher moving off the Fortran already knows.
 *
 * @code{.purr}
 * import io
 * import fixarray
 * import space.irbem.purr as irbem
 *
 * let e = irbem.epoch_at(2015, 182, 43200.0)          # IGRF-14 + the frame rotations, once
 * let c = irbem.make_lstar(e, 6.6, 0.0, 0.0, 1, 90.0) # geosynchronous, 90 deg pitch angle
 * # c[0] Lm, c[1] L*, c[2] Blocal, c[3] Bmin, c[4] I, c[5] MLT
 * @endcode
 */

#include <cmath>
#include <optional>
#include <string>

#include "ndarray.hpp"

#include "api.hpp"
#include "driftshell.hpp"
#include "frames.hpp"
#include "igrf.hpp"
#include "status.hpp"

namespace cheatah::space::irbem::purr {

/// The internal-field truncation the facade pins. Degree 13 is IGRF-14 as IAGA publishes it, and
/// fixing it here is what removes @ref Igrf's non-type template parameter from the surface.
using Model = Igrf<tables::igrf14_max_degree, Exact>;

/**
 * An epoch: the internal field model and the frame rotations for one instant.
 *
 * A cheatah program cannot call `Igrf<13>::at` — a static member on a class template, returning
 * `std::optional` — so this absorbs it. Build one, hand it to everything else. It is a plain
 * aggregate, so a `.purr` program reads @ref ok directly.
 */
struct Epoch {
    /// IGRF at this instant. Held in an `optional` because @ref Igrf's constructor is PRIVATE —
    /// only its static factory builds one — so a default-constructible aggregate cannot hold it
    /// directly. That privacy is exactly what a cheatah program cannot work around, and absorbing
    /// it is this type's whole reason to exist.
    std::optional<Model> model;
    Rotations rotations;  ///< The frame rotations for it.
    bool ok = false;      ///< Whether the epoch was inside IGRF's definition; see @ref epoch_at.
};

/**
 * The epoch for a date and time.
 * @param year the calendar year. @param doy the day of year, 1-366.
 * @param ut_seconds seconds since midnight UT.
 * @return the epoch; @ref Epoch::ok is false when the date lies outside IGRF's published range,
 *         in which case every routine taking it returns @ref Status::OutOfValidityRange rather
 *         than a plausible wrong number.
 * @complexity O(NMAX²) once, then reused by every call that takes it.
 * @alloc none.
 * @test IrbemPurr.EpochMatchesTheTypedCore
 * @systest systests/test_irbem.purr
 *
 * @code{.purr}
 * import fixarray
 * import space.irbem.purr as irbem
 * let e = irbem.epoch_at(2015, 182, 43200.0)
 * @endcode
 */
[[nodiscard]] inline Epoch epoch_at(int year, int doy, double ut_seconds) {
    Epoch e;
    const DateTime dt = date_and_time_from_doy_and_ut(year, doy, ut_seconds);
    const double decy = decimal_year(dt.year, dt.month, dt.day, dt.hour, dt.minute, dt.second);
    const std::optional<Model> m = Model::at(decy);
    if (!m.has_value()) return e;
    e.model = m;
    const Result<Rotations> r = api::rotations_at(year, doy, ut_seconds, *e.model);
    if (r.status != Status::Ok) return e;
    e.rotations = r.value;
    e.ok = true;
    return e;
}

/**
 * Whether @p e names an instant this model can evaluate.
 * @param e the epoch. @return @ref Epoch::ok.
 * @complexity O(1). @alloc none.
 * @test IrbemPurr.EpochMatchesTheTypedCore
 * @systest systests/test_irbem.purr
 *
 * @code{.purr}
 * import fixarray
 * import space.irbem.purr as irbem
 * let e = irbem.epoch_at(1750, 1, 0.0)
 * let usable = irbem.epoch_ok(e)          # false — before IGRF begins
 * @endcode
 */
[[nodiscard]] inline bool epoch_ok(const Epoch& e) { return e.ok; }

/**
 * The name of a @ref Status code, for a program that wants to report why a call declined.
 * @param code the status, as returned in a result array's status slot.
 * @return a short name such as `"Ok"` or `"OutOfValidityRange"`; `"Unknown"` for an unrecognised
 *         code, so this never throws on a value from an older build.
 * @complexity O(1). @alloc the returned string.
 * @test IrbemPurr.StatusNamesEveryCode
 * @systest systests/test_irbem.purr
 *
 * @code{.purr}
 * import io
 * import fixarray
 * import space.irbem.purr as irbem
 * io.print(irbem.status_name(0))          # -> Ok
 * @endcode
 */
[[nodiscard]] inline std::string status_name(int code) {
    switch (static_cast<Status>(code)) {
        case Status::Ok: return "Ok";
        case Status::OutOfValidityRange: return "OutOfValidityRange";
        case Status::OpenFieldLine: return "OpenFieldLine";
        case Status::NotConverged: return "NotConverged";
        case Status::ParametersMissing: return "ParametersMissing";
        case Status::DomainError: return "DomainError";
    }
    return "Unknown";
}

/// How many values @ref make_lstar returns, and the slot each occupies.
inline constexpr long long lstar_slots = 7;

/**
 * The magnetic coordinates of one point — IRBEM's `MAKE_LSTAR`.
 *
 * @param e the epoch, from @ref epoch_at.
 * @param x1 first position component. @param x2 second. @param x3 third — read per @p sysaxes.
 * @param sysaxes IRBEM's frame code: 0 GDZ, 1 GEO, 2 GSM, 3 GSE, 4 SM, 5 GEI, 6 MAG, 7 SPH, 8 RLL.
 * @param pitch_angle_deg the local pitch angle; 90 is IRBEM's `MAKE_LSTAR` convention.
 * @return an `ndarray` of @ref lstar_slots values: `Lm`, `L*`, `Blocal`, `Bmin`, `I`, `MLT`, and
 *         the @ref Status code as the last slot. A declined call fills the physical slots with
 *         IRBEM's `baddata` and names the reason in the status slot rather than returning nothing.
 * @complexity ~10⁵ field evaluations — the drift-shell integral.
 * @alloc one array of @ref lstar_slots doubles.
 * @test IrbemPurr.MakeLstarMatchesTheTypedCore
 * @systest systests/test_irbem.purr
 *
 * @code{.purr}
 * import io
 * import fixarray
 * import space.irbem.purr as irbem
 * let e = irbem.epoch_at(2015, 182, 43200.0)
 * let c = irbem.make_lstar(e, 6.6, 0.0, 0.0, 1, 90.0)
 * io.print("L* =", c[1])
 * @endcode
 */
[[nodiscard]] inline ::cheatah::ndarray::NDArray make_lstar(
    const Epoch& e, double x1, double x2, double x3, int sysaxes, double pitch_angle_deg) {
    auto out = ::cheatah::ndarray::zeros({lstar_slots});
    auto fill = [&out](const api::MagneticCoordinates& c, Status st) {
        out[0] = c.lm;
        out[1] = c.lstar;
        out[2] = c.blocal;
        out[3] = c.bmin;
        out[4] = c.xj;
        out[5] = c.mlt;
        out[6] = static_cast<double>(static_cast<int>(st));
    };
    if (!e.ok) {
        fill(api::MagneticCoordinates{}, Status::OutOfValidityRange);
        return out;
    }
    const Result<fixarray::vec3d> geo =
        api::detail::to_geo(frame_from_sysaxes(sysaxes).value_or(Frame::GEO),
                    fixarray::vec3d{x1, x2, x3}, e.rotations);
    if (!frame_from_sysaxes(sysaxes).has_value() || geo.status != Status::Ok) {
        fill(api::MagneticCoordinates{}, Status::DomainError);
        return out;
    }
    const Result<api::MagneticCoordinates> r =
        api::make_lstar(*e.model, e.rotations, Position<Frame::GEO>{geo.value},
                        ExternalModel::None, pitch_angle_deg, DriftShellOptions{});
    fill(r.value, r.status);
    return out;
}

/**
 * Magnetic local time at a geographic point — IRBEM's `GET_MLT`.
 * @param e the epoch, from @ref epoch_at.
 * @param x1 first position component.
 * @param x2 second position component.
 * @param x3 third position component — read per @p sysaxes.
 * @param sysaxes the frame code the components are in.
 * @return MLT in hours, folded into `[0, 24)`; a negative value signals a declined call, which
 *         cannot collide with a real MLT.
 * @complexity O(1) beyond the frame transform. @alloc none.
 * @test IrbemPurr.MltMatchesTheTypedCore
 * @systest systests/test_irbem.purr
 *
 * @code{.purr}
 * import io
 * import fixarray
 * import space.irbem.purr as irbem
 * let e = irbem.epoch_at(2015, 182, 43200.0)
 * io.print("MLT =", irbem.get_mlt(e, 6.6, 0.0, 0.0, 1))
 * @endcode
 */
[[nodiscard]] inline double get_mlt(const Epoch& e, double x1, double x2, double x3, int sysaxes) {
    if (!e.ok) return -1.0;
    const std::optional<Frame> f = frame_from_sysaxes(sysaxes);
    if (!f.has_value()) return -1.0;
    const Result<fixarray::vec3d> geo = api::detail::to_geo(*f, fixarray::vec3d{x1, x2, x3}, e.rotations);
    if (geo.status != Status::Ok) return -1.0;
    const Result<double> r = api::get_mlt(Position<Frame::GEO>{geo.value}, e.rotations);
    if (r.status != Status::Ok) return -1.0;
    return r.value;
}

/**
 * Convert a position between any two of IRBEM's frames — its `COORD_TRANS`.
 * @param e the epoch, whose rotations the transform needs.
 * @param x1 first input component.
 * @param x2 second input component.
 * @param x3 third input component — read per @p sysaxes_in.
 * @param sysaxes_in the frame the input is in.
 * @param sysaxes_out the frame to convert to.
 * @return an `ndarray` of four values: the three output components, then the @ref Status code.
 * @complexity O(1) — at most two 3×3 products. @alloc one array of four doubles.
 * @test IrbemPurr.CoordTransRoundTrips
 * @systest systests/test_irbem.purr
 *
 * @code{.purr}
 * import fixarray
 * import space.irbem.purr as irbem
 * let e = irbem.epoch_at(2015, 182, 43200.0)
 * let gsm = irbem.coord_trans(e, 6.6, 0.0, 0.0, 1, 2)   # GEO -> GSM
 * @endcode
 */
[[nodiscard]] inline ::cheatah::ndarray::NDArray coord_trans(
    const Epoch& e, double x1, double x2, double x3, int sysaxes_in, int sysaxes_out) {
    auto out = ::cheatah::ndarray::zeros({4});
    out[3] = static_cast<double>(static_cast<int>(Status::DomainError));
    if (!e.ok) {
        out[3] = static_cast<double>(static_cast<int>(Status::OutOfValidityRange));
        return out;
    }
    const Result<fixarray::vec3d> r =
        api::coord_trans(sysaxes_in, sysaxes_out, fixarray::vec3d{x1, x2, x3}, e.rotations);
    out[0] = r.value[0];
    out[1] = r.value[1];
    out[2] = r.value[2];
    out[3] = static_cast<double>(static_cast<int>(r.status));
    return out;
}

}  // namespace cheatah::space::irbem::purr
