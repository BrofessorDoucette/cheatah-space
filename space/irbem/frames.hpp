#pragma once

/**
 * @file frames.hpp
 * @brief space.irbem — the reference frames, and the position/field types tagged by them.
 *
 * Magnetospheric physics is done in a dozen coordinate frames at once: a spacecraft ephemeris
 * arrives in geodetic or geographic coordinates, the internal field is evaluated in geographic,
 * every external field model is defined in solar-magnetospheric, and drift shells are traced in
 * solar-magnetic. Silently handing one frame's numbers to a routine expecting another is THE
 * classic defect of this domain — the components are all plausible, so the answer is merely wrong
 * rather than obviously broken.
 *
 * So the frame lives in the type. @ref Position and @ref FieldVector are zero-cost wrappers around
 * a `cheatah::fixarray::vec3d` whose template parameter is the frame, which makes a mismatch a
 * COMPILE error and puts the transform graph itself under compile-time check. There is no runtime
 * cost and no runtime branch: the wrapper holds exactly one vector and nothing else.
 *
 * The geometry itself is `cheatah::fixarray` — inline storage, compile-time trip counts, nothing
 * allocated, and matrices stored column-major so a frame rotation `m * v` is a sum of contiguous
 * scaled columns. This header adds the vocabulary, not the arithmetic.
 *
 * IRBEM's C API selects a frame with a runtime `sysaxes` integer. That integer enters the typed
 * world through @ref frame_from_sysaxes at the API boundary and nowhere else.
 *
 * @note Frames differ in what their three components MEAN, which is why @ref FrameKind exists:
 *       a Cartesian frame carries (x, y, z) in Earth radii, a spherical one carries
 *       (radius, latitude, east longitude), and the geodetic one carries
 *       (altitude, latitude, east longitude) with altitude in kilometres above the reference
 *       ellipsoid. Mixing those is the same class of bug one level down.
 */

#include <cstdint>
#include <optional>
#include <string_view>

#include "cheatah.hpp"
#include "fixarray.hpp"

namespace cheatah::space::irbem {

/**
 * How a frame's three components are to be read.
 *
 * The distinction is load-bearing: `radius()` means something for a spherical frame and nothing
 * for a Cartesian one, so the accessors below are constrained on it rather than left to a comment.
 */
enum class FrameKind : std::uint8_t {
    Cartesian,  ///< (x, y, z) in Earth radii.
    Spherical,  ///< (radius Re, latitude deg, east longitude deg).
    Geodetic,   ///< (altitude km, latitude deg, east longitude deg), on the reference ellipsoid.
};

/**
 * A reference frame.
 *
 * The first nine enumerators carry IRBEM's `sysaxes` codes as their values, so the boundary
 * conversion in @ref frame_from_sysaxes is a range check rather than a table. The heliospheric
 * three have no `sysaxes` code — IRBEM reaches them through dedicated routines — so they are
 * numbered above the `sysaxes` range and @ref sysaxes_of reports them as absent.
 */
enum class Frame : std::uint8_t {
    GDZ = 0,   ///< Geodetic: altitude (km), latitude, east longitude.
    GEO = 1,   ///< Geographic (Earth-fixed) Cartesian.
    GSM = 2,   ///< Geocentric Solar Magnetospheric — the frame the external field models live in.
    GSE = 3,   ///< Geocentric Solar Ecliptic.
    SM = 4,    ///< Solar Magnetic — the frame drift shells are natural in.
    GEI = 5,   ///< Geocentric Equatorial Inertial (true equator and equinox of date).
    MAG = 6,   ///< Geomagnetic (centred dipole).
    SPH = 7,   ///< Geographic spherical: radius (Re), latitude, east longitude.
    RLL = 8,   ///< Radius/latitude/longitude — as @ref Frame::SPH; IRBEM prefers this spelling.
    HEE = 9,   ///< Heliocentric Earth Ecliptic.
    HAE = 10,  ///< Heliocentric Aries Ecliptic.
    HEEQ = 11, ///< Heliocentric Earth Equatorial.
};

/// The number of frames — the bound of any loop over @ref Frame.
inline constexpr std::size_t frame_count = 12;

/**
 * How to read @p f's three components.
 * @param f the frame.
 * @return the component convention; @ref FrameKind::Cartesian for everything but GDZ/SPH/RLL.
 * @complexity O(1).
 * @alloc none.
 */
constexpr FrameKind kind_of(Frame f) {
    switch (f) {
        case Frame::GDZ:
            return FrameKind::Geodetic;
        case Frame::SPH:
        case Frame::RLL:
            return FrameKind::Spherical;
        default:
            return FrameKind::Cartesian;
    }
}

/**
 * The frame's short name, as the literature and IRBEM's own documentation spell it.
 * @param f the frame.
 * @return a static string such as `"GSM"`; never empty.
 * @complexity O(1).
 * @alloc none.
 */
constexpr std::string_view name_of(Frame f) {
    switch (f) {
        case Frame::GDZ:
            return "GDZ";
        case Frame::GEO:
            return "GEO";
        case Frame::GSM:
            return "GSM";
        case Frame::GSE:
            return "GSE";
        case Frame::SM:
            return "SM";
        case Frame::GEI:
            return "GEI";
        case Frame::MAG:
            return "MAG";
        case Frame::SPH:
            return "SPH";
        case Frame::RLL:
            return "RLL";
        case Frame::HEE:
            return "HEE";
        case Frame::HAE:
            return "HAE";
        case Frame::HEEQ:
            return "HEEQ";
    }
    return "?";  // unreachable for a valid enumerator; keeps the function total
}

/**
 * The IRBEM `sysaxes` code naming @p f, when one exists.
 * @param f the frame.
 * @return the code `0..8`, or `std::nullopt` for the heliospheric frames, which `sysaxes` cannot
 *         name (IRBEM reaches those through dedicated transform routines instead).
 * @complexity O(1).
 * @alloc none.
 */
constexpr std::optional<int> sysaxes_of(Frame f) {
    const auto raw = static_cast<int>(f);
    if (raw > static_cast<int>(Frame::RLL)) return std::nullopt;
    return raw;
}

/**
 * The frame an IRBEM `sysaxes` code names — the ONE place a runtime frame selector enters the
 * typed world.
 * @param sysaxes the IRBEM code, `0..8`.
 * @return the frame, or `std::nullopt` when @p sysaxes is outside that range. Callers turn the
 *         empty case into a named error; it is never silently defaulted to a frame.
 * @complexity O(1).
 * @alloc none.
 */
constexpr std::optional<Frame> frame_from_sysaxes(int sysaxes) {
    if (sysaxes < 0 || sysaxes > static_cast<int>(Frame::RLL)) return std::nullopt;
    return static_cast<Frame>(static_cast<std::uint8_t>(sysaxes));
}

/// A frame whose components are `(x, y, z)` in Earth radii.
template <Frame F>
concept CartesianFrame = (kind_of(F) == FrameKind::Cartesian);

/// A frame whose components are an angular pair over a radius or altitude.
template <Frame F>
concept AngularFrame = (kind_of(F) != FrameKind::Cartesian);

/**
 * A position, in frame @p F.
 *
 * An aggregate holding exactly one `vec3d`, so it is the same size, trivially copyable, and
 * passes in registers — the frame tag costs nothing at runtime and everything at compile time.
 *
 * @tparam F the frame the components are expressed in; see @ref kind_of for what they mean.
 */
template <Frame F>
struct Position {
    /// The frame these components are expressed in.
    static constexpr Frame frame = F;

    /// The components — (x, y, z) in Re, or (radius|altitude, latitude, east longitude).
    /// Zero by default: `fixarray::Fixed`'s defaulted constructor value-initializes, which the
    /// `DefaultsAreZero` test pins down.
    fixarray::vec3d v;

    /// Geocentric distance, for a Cartesian frame.
    /// @return `|v|` in Earth radii.
    /// @complexity O(1).
    /// @alloc none.
    /// @note Not `constexpr`: it takes a square root, and `fixarray::norm` is a runtime function.
    ///       Use `fixarray::squared_norm(p.v)` when a compile-time comparison will do.
    [[nodiscard]] double radius() const
        requires CartesianFrame<F>
    {
        return fixarray::norm(v);
    }

    /// Geocentric distance, for a spherical frame, where it is simply the first component.
    /// @return the radius in Earth radii (or the altitude in km, for @ref Frame::GDZ).
    /// @complexity O(1).
    /// @alloc none.
    [[nodiscard]] constexpr double radius() const
        requires AngularFrame<F>
    {
        return v[0];
    }

    /// Latitude, for an angular frame.
    /// @return degrees, north positive.
    /// @complexity O(1).
    /// @alloc none.
    [[nodiscard]] constexpr double latitude() const
        requires AngularFrame<F>
    {
        return v[1];
    }

    /// East longitude, for an angular frame.
    /// @return degrees, east positive.
    /// @complexity O(1).
    /// @alloc none.
    [[nodiscard]] constexpr double longitude() const
        requires AngularFrame<F>
    {
        return v[2];
    }

    /// Two positions are equal when their components are.
    /// @param a the left position. @param b the right position.
    /// @return whether every component matches exactly.
    /// @complexity O(1).
    /// @alloc none.
    friend constexpr bool operator==(const Position& a, const Position& b) { return a.v == b.v; }
};

/**
 * A magnetic field vector in nanotesla, in frame @p F.
 *
 * Superposition is the operation that matters: a total field is an internal model plus one or
 * more external contributions, and they must all be in the same frame. Making `+` frame-preserving
 * means that sum cannot be written wrongly.
 *
 * @tparam F the frame the components are expressed in; always Cartesian in practice, since the
 *           field models are defined on Cartesian frames.
 */
template <Frame F>
struct FieldVector {
    /// The frame these components are expressed in.
    static constexpr Frame frame = F;

    /// The components, in nanotesla. Zero by default — see @ref Position::v.
    fixarray::vec3d v;

    /// Field magnitude — the `B` of `Blocal`, `Bmin`, `Bmirr`.
    /// @return `|v|` in nanotesla.
    /// @complexity O(1).
    /// @alloc none.
    /// @note Not `constexpr`, for the same reason as @ref Position::radius.
    [[nodiscard]] double magnitude() const { return fixarray::norm(v); }

    /// Superposition: internal plus external, in one frame.
    /// @param a the left field. @param b the right field.
    /// @return their componentwise sum.
    /// @complexity O(1).
    /// @alloc none.
    friend constexpr FieldVector operator+(FieldVector a, FieldVector b) {
        return FieldVector{a.v + b.v};
    }

    /// Componentwise difference — used to form finite-difference field gradients.
    /// @param a the left field. @param b the right field.
    /// @return their componentwise difference.
    /// @complexity O(1).
    /// @alloc none.
    friend constexpr FieldVector operator-(FieldVector a, FieldVector b) {
        return FieldVector{a.v - b.v};
    }

    /// Scale the field.
    /// @param a the field. @param s the scale factor.
    /// @return @p a with every component multiplied by @p s.
    /// @complexity O(1).
    /// @alloc none.
    friend constexpr FieldVector operator*(FieldVector a, double s) { return FieldVector{a.v * s}; }

    /// Two field vectors are equal when their components are.
    /// @param a the left field. @param b the right field.
    /// @return whether every component matches exactly.
    /// @complexity O(1).
    /// @alloc none.
    friend constexpr bool operator==(const FieldVector& a, const FieldVector& b) {
        return a.v == b.v;
    }
};

}  // namespace cheatah::space::irbem
