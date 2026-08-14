#pragma once

/**
 * @file time.hpp
 * @brief space.time — Julian dates & epoch conversions (the pun is intended: spacetime).
 *
 * An extension *for* the cheatah standard time library: it adds the time scales astronomy
 * and space-physics actually use — Julian Date, Modified Julian Date, the J2000 epoch, and
 * NASA's CDF_EPOCH (the bridge to `space.cdf`). `import space.time` resolves this header.
 *
 * Every conversion is a single concept-constrained template that accepts EITHER a scalar
 * `Value` that makes sense for the function (a number of seconds/days/ms), a numeric
 * `ndarray` (a whole mission's worth of timestamps, vectorized over the ndarray SIMD path),
 * or — once `cheatah::datetime` defines its struct — an ndarray of those datetime structs.
 * The constraints are explicit so a researcher who passes the wrong thing gets a crisp
 * compiler error instead of a deep template spew.
 *
 * Cross-platform, header-only, allocation-free for scalars: pure C++20 arithmetic with no
 * platform headers, no I/O, and no global state. The current clock is intentionally NOT a
 * dependency here — compose with the standard `datetime`: `unix_to_jd(datetime.timestamp())`.
 *
 * @note `n` below is the element count of an ndarray input (1 for a scalar).
 */

#include "cheatah.hpp"

// ndarray support is OPT-IN by include path: a program that passes ndarrays imports `ndarray`,
// which puts this header on purrc's C++ include path; a scalar-only program never needs it.
#if __has_include("ndarray.hpp")
#  include "ndarray.hpp"
#  define CHEATAH_SPACE_HAVE_NDARRAY 1
#endif

#include <type_traits>

namespace cheatah::space::time {

// ---- input concepts ---------------------------------------------------------------------

/// A scalar that makes sense as a continuous time/epoch quantity — any arithmetic type.
template <class T>
concept Numeric = std::is_arithmetic_v<std::remove_cvref_t<T>>;

/**
 * Trait hook for cheatah's forthcoming datetime struct. It defaults to `false`; when
 * `cheatah::datetime` defines its broken-down datetime type, this is specialized to `true`
 * for it (here or there), which lights up @ref DatetimeScalar / @ref DatetimeArray and the
 * datetime overloads — without this header taking a dependency on that type today.
 */
template <class T>
inline constexpr bool is_datetime_v = false;

/// A single broken-down datetime value (see @ref is_datetime_v).
template <class T>
concept DatetimeScalar = is_datetime_v<std::remove_cvref_t<T>>;

#ifdef CHEATAH_SPACE_HAVE_NDARRAY
namespace detail {
/// Whether @p A is a `cheatah::ndarray::basic_ndarray<…>` (any element type).
template <class>
inline constexpr bool is_ndarray_v = false;
template <class T>
inline constexpr bool is_ndarray_v<::cheatah::ndarray::basic_ndarray<T>> = true;

/// The element (`value_type`) of an ndarray input.
template <class A>
using element_t = typename std::remove_cvref_t<A>::value_type;
}  // namespace detail

/// An ndarray with a numeric (arithmetic) entry — the vectorized continuous-time input.
template <class A>
concept NumericArray =
    detail::is_ndarray_v<std::remove_cvref_t<A>> && std::is_arithmetic_v<detail::element_t<A>>;

/// An ndarray of datetime structs (see @ref is_datetime_v) — lit up when that struct lands.
template <class A>
concept DatetimeArray =
    detail::is_ndarray_v<std::remove_cvref_t<A>> && is_datetime_v<detail::element_t<A>>;
#else
template <class A> concept NumericArray = false;
template <class A> concept DatetimeArray = false;
#endif

/// The full accepted domain: a sensible scalar, a numeric ndarray, or a datetime-struct ndarray.
template <class T>
concept TimeInput = Numeric<T> || DatetimeScalar<T> || NumericArray<T> || DatetimeArray<T>;

// ---- reference epochs -------------------------------------------------------------------
// Exposed as functions so callers never hard-code the magic numbers.

/**
 * Julian Date of the Unix epoch, 1970-01-01T00:00:00Z.
 * @return `2440587.5`.
 * @complexity O(1).
 * @alloc none.
 * @systest systests/test_jd.purr
 * @par Example
 * @code{.purr}
 * import io
 * import space.time as st
 *
 * # The Unix epoch's Julian Date is the anchor of every unix<->jd conversion.
 * io.print(st.jd_unix_epoch())                       # 2440587.5
 * io.print(st.unix_to_jd(0.0) == st.jd_unix_epoch()) # true
 * @endcode
 */
constexpr double jd_unix_epoch() { return 2440587.5; }

/**
 * Julian Date of the J2000.0 epoch, 2000-01-01T12:00:00 TT.
 * @return `2451545.0`.
 * @complexity O(1).
 * @alloc none.
 * @systest systests/test_jd.purr
 * @par Example
 * @code{.purr}
 * import io
 * import space.time as st
 *
 * # J2000.0 — the reference epoch of modern astronomical series.
 * io.print(st.jd_j2000())                            # 2451545.0
 * io.print(st.jd_to_j2000_seconds(st.jd_j2000()))    # 0.0 — J2000 is its own origin
 * @endcode
 */
constexpr double jd_j2000() { return 2451545.0; }

/**
 * CDF_EPOCH of the Unix epoch: milliseconds from 0000-01-01T00:00:00.000 to 1970-01-01 —
 * the constant that bridges Unix time and NASA CDF's millisecond epoch.
 * @return `62167219200000.0`.
 * @complexity O(1).
 * @alloc none.
 * @systest systests/test_cdf_epoch.purr
 * @par Example
 * @code{.purr}
 * import io
 * import space.time as st
 *
 * # The Unix epoch expressed as a CDF_EPOCH (ms since year 0).
 * io.print(st.cdf_epoch_unix_offset_ms())                             # 62167219200000.0
 * io.print(st.unix_to_cdf_epoch(0.0) == st.cdf_epoch_unix_offset_ms())  # true
 * @endcode
 */
constexpr double cdf_epoch_unix_offset_ms() { return 62167219200000.0; }

// ---- conversions ------------------------------------------------------------------------

/**
 * Unix epoch seconds → Julian Date.
 * @param seconds scalar or numeric ndarray of seconds since 1970-01-01T00:00:00Z.
 * @return the corresponding Julian Date(s).
 * @complexity O(n).
 * @alloc one result array for an ndarray; none for a scalar.
 * @systest systests/test_jd.purr
 * @systest systests/test_vectorized.purr
 * @par Example
 * @code{.purr}
 * import io
 * import ndarray
 * import space.time as st
 *
 * io.print(st.unix_to_jd(946728000.0))               # 2451545.0 — J2000, exactly
 * let jds = st.unix_to_jd(ndarray.array([0.0, 86400.0]))   # a whole array, vectorized
 * io.print(jds[0], jds[1])                           # 2440587.5 2440588.5
 * @endcode
 */
auto unix_to_jd(TimeInput auto&& seconds) {
    return seconds / 86400.0 + 2440587.5;
}

/**
 * Julian Date → Unix epoch seconds.
 * @param jd scalar or numeric ndarray of Julian Dates.
 * @return seconds since 1970-01-01T00:00:00Z.
 * @complexity O(n).
 * @alloc one result array for an ndarray; none for a scalar.
 * @systest systests/test_jd.purr
 * @systest systests/test_vectorized.purr
 * @par Example
 * @code{.purr}
 * import io
 * import space.time as st
 *
 * io.print(st.jd_to_unix(2451545.0))                 # 946728000.0 — J2000 in Unix seconds
 * io.print(st.jd_to_unix(st.unix_to_jd(86400.0)))    # 86400.0 — the round trip
 * @endcode
 */
auto jd_to_unix(TimeInput auto&& jd) {
    return (jd - 2440587.5) * 86400.0;
}

/**
 * Julian Date → Modified Julian Date (MJD = JD − 2400000.5).
 * @param jd scalar or numeric ndarray of Julian Dates.
 * @return the Modified Julian Date(s).
 * @complexity O(n).
 * @alloc one result array for an ndarray; none for a scalar.
 * @systest systests/test_mjd.purr
 * @par Example
 * @code{.purr}
 * import io
 * import space.time as st
 *
 * io.print(st.jd_to_mjd(st.jd_j2000()))              # 51544.5 — J2000 as an MJD
 * @endcode
 */
auto jd_to_mjd(TimeInput auto&& jd) {
    return jd - 2400000.5;
}

/**
 * Modified Julian Date → Julian Date.
 * @param mjd scalar or numeric ndarray of Modified Julian Dates.
 * @return the Julian Date(s).
 * @complexity O(n).
 * @alloc one result array for an ndarray; none for a scalar.
 * @systest systests/test_mjd.purr
 * @par Example
 * @code{.purr}
 * import io
 * import space.time as st
 *
 * io.print(st.mjd_to_jd(51544.5))                    # 2451545.0 — back to a full JD
 * @endcode
 */
auto mjd_to_jd(TimeInput auto&& mjd) {
    return mjd + 2400000.5;
}

/**
 * Unix epoch seconds → Modified Julian Date.
 * @param seconds scalar or numeric ndarray of seconds since the Unix epoch.
 * @return the Modified Julian Date(s).
 * @complexity O(n).
 * @alloc one result array for an ndarray; none for a scalar.
 * @systest systests/test_mjd.purr
 * @systest systests/test_vectorized.purr
 * @par Example
 * @code{.purr}
 * import io
 * import space.time as st
 *
 * io.print(st.unix_to_mjd(0.0))                      # 40587.0 — the Unix epoch's MJD
 * @endcode
 */
auto unix_to_mjd(TimeInput auto&& seconds) {
    return jd_to_mjd(unix_to_jd(seconds));
}

/**
 * Modified Julian Date → Unix epoch seconds.
 * @param mjd scalar or numeric ndarray of Modified Julian Dates.
 * @return seconds since the Unix epoch.
 * @complexity O(n).
 * @alloc one result array for an ndarray; none for a scalar.
 * @systest systests/test_mjd.purr
 * @par Example
 * @code{.purr}
 * import io
 * import space.time as st
 *
 * io.print(st.mjd_to_unix(40587.0))                  # 0.0 — the Unix epoch again
 * @endcode
 */
auto mjd_to_unix(TimeInput auto&& mjd) {
    return jd_to_unix(mjd_to_jd(mjd));
}

/**
 * Seconds elapsed since the J2000.0 epoch for a given Julian Date.
 * @param jd scalar or numeric ndarray of Julian Dates.
 * @return seconds since 2000-01-01T12:00:00.
 * @complexity O(n).
 * @alloc one result array for an ndarray; none for a scalar.
 * @systest systests/test_j2000.purr
 * @par Example
 * @code{.purr}
 * import io
 * import space.time as st
 *
 * io.print(st.jd_to_j2000_seconds(2451546.0))        # 86400.0 — one day after J2000
 * @endcode
 */
auto jd_to_j2000_seconds(TimeInput auto&& jd) {
    return (jd - 2451545.0) * 86400.0;
}

/**
 * Julian centuries (36525 days) elapsed since J2000.0 — the time argument for most
 * astronomical polynomial series (mean sidereal time, precession, …).
 * @param jd scalar or numeric ndarray of Julian Dates.
 * @return Julian centuries since J2000.0.
 * @complexity O(n).
 * @alloc one result array for an ndarray; none for a scalar.
 * @systest systests/test_j2000.purr
 * @par Example
 * @code{.purr}
 * import io
 * import space.time as st
 *
 * # The `T` fed to precession / sidereal-time polynomial series.
 * io.print(st.jd_to_j2000_centuries(2451545.0 + 36525.0))   # 1.0 — one Julian century
 * @endcode
 */
auto jd_to_j2000_centuries(TimeInput auto&& jd) {
    return (jd - 2451545.0) / 36525.0;
}

/**
 * Unix epoch seconds → NASA CDF_EPOCH (milliseconds since year 0). Feeds `space.cdf`
 * epoch variables directly. CDF_TT2000 (ns since J2000 with leap seconds) needs a
 * leap-second table and lands with `space.cdf` — see the cdf roadmap.
 * @param seconds scalar or numeric ndarray of seconds since the Unix epoch.
 * @return the CDF_EPOCH value(s) in milliseconds.
 * @complexity O(n).
 * @alloc one result array for an ndarray; none for a scalar.
 * @systest systests/test_cdf_epoch.purr
 * @systest systests/test_vectorized.purr
 * @par Example
 * @code{.purr}
 * import io
 * import space.time as st
 *
 * # J2000 as a CDF_EPOCH — ready to write into a CDF epoch variable.
 * io.print(st.unix_to_cdf_epoch(946728000.0))        # 63113947200000.0
 * @endcode
 */
auto unix_to_cdf_epoch(TimeInput auto&& seconds) {
    return seconds * 1000.0 + 62167219200000.0;
}

/**
 * NASA CDF_EPOCH (milliseconds since year 0) → Unix epoch seconds.
 * @param ms scalar or numeric ndarray of CDF_EPOCH values in milliseconds.
 * @return seconds since the Unix epoch.
 * @complexity O(n).
 * @alloc one result array for an ndarray; none for a scalar.
 * @systest systests/test_cdf_epoch.purr
 * @par Example
 * @code{.purr}
 * import io
 * import space.time as st
 *
 * # An epoch read from a CDF file, back onto the Unix scale.
 * io.print(st.cdf_epoch_to_unix(63113947200000.0))   # 946728000.0 — J2000
 * @endcode
 */
auto cdf_epoch_to_unix(TimeInput auto&& ms) {
    return (ms - 62167219200000.0) / 1000.0;
}

}  // namespace cheatah::space::time
