#pragma once

/**
 * @file library_info.hpp
 * @brief space.irbem — what this implementation is, for the routines IRBEM exposes to ask.
 *
 * IRBEM lets a caller interrogate the library it linked against: which IGRF generation the internal
 * field implements, how large a batch it accepts, and which build it is. Those questions are about
 * THIS implementation, so they live here.
 *
 * The calendar conversions that used to share this file moved to @ref space/time/calendar.hpp,
 * where they belong: Julian day numbers and decimal years are generic time arithmetic with no
 * space-physics content, and keeping a second copy here duplicated `is_leap_year` against
 * `space/time/civil.hpp`.
 */

#include "../time/calendar.hpp"

namespace cheatah::space::irbem {

/// The calendar conversions, re-exported so `space.irbem` callers need not name `space.time` for
/// something IRBEM presents as part of its own surface. One definition, two spellings.
using space::time::CalendarDate;
using space::time::DateTime;
using space::time::calendar_date;
using space::time::date_and_time_from_decimal_year;
using space::time::date_and_time_from_doy_and_ut;
using space::time::date_from_day_of_year;
using space::time::day_of_year;
using space::time::days_in_year;
using space::time::decimal_year;
using space::time::gregorian_reform_jdn;
using space::time::is_leap_year;
using space::time::julian_day_number;

/**
 * The batch cap on the time dimension of IRBEM's array entry points — IRBEM's
 * `GET_IRBEM_NTIME_MAX`.
 *
 * @return `100000`. That figure is IRBEM's own dimensioning constant, taken from the generated
 *         `source/ntime_max.inc` (`PARAMETER (NTIME_MAX = 100000)`) and confirmed by calling
 *         `get_irbem_ntime_max1_` on the shipped shared library; it is a number, not an algorithm.
 * @complexity O(1).
 * @alloc none.
 *
 * @note Nothing in this implementation is actually limited to that many epochs — the batch lanes
 *       take a caller-provided span of any length. The value exists so a port of an IRBEM-shaped
 *       program that sizes its own arrays by this number keeps working unchanged.
 */
constexpr int max_batch_times() { return 100000; }

/**
 * The IGRF generation this module implements — IRBEM's `GET_IGRF_VERSION`.
 *
 * @return `14`, i.e. IGRF-14 (IAGA, released 2024, valid 1900.0–2030.0). The shipped IRBEM library
 *         returns the same, measured.
 * @complexity O(1).
 * @alloc none.
 */
constexpr int igrf_generation() { return 14; }

/**
 * Implementation version — the compatibility shim standing in for IRBEM's
 * `IRBEM_FORTRAN_VERSION`.
 *
 * @return a monotonically increasing integer identifying **this** C++ implementation, starting at
 *         `1` and bumped only when the behaviour of this shim pair changes.
 * @complexity O(1).
 * @alloc none.
 *
 * @warning There is no Fortran here. IRBEM's routine reports the revision of its own Fortran
 *          sources, and returning a plausible-looking IRBEM revision would let a caller
 *          feature-detect against a value that means nothing. Comparing this number against an
 *          IRBEM revision is therefore meaningless; use @ref implementation_release, which says in
 *          words what is running.
 */
constexpr int implementation_version() { return 1; }

/**
 * Implementation release tag — the compatibility shim standing in for
 * `IRBEM_FORTRAN_RELEASE`.
 *
 * @return a static, never-empty string naming this implementation. It says "not IRBEM Fortran"
 *         in so many words, so a log line carrying it cannot be mistaken for the real library's
 *         release tag (which is a git short hash, measured).
 * @complexity O(1).
 * @alloc none — the returned view refers to a string literal with static storage duration.
 *
 * @note IRBEM's C entry point writes an 80-character space-padded buffer. This value is short
 *       enough to fit; padding it belongs in the C boundary layer, not here, so that the ordinary
 *       C++ caller gets a plain view with no trailing blanks to strip.
 */
constexpr std::string_view implementation_release() {
    return "cheatah-space space.irbem (not IRBEM Fortran)";
}

}  // namespace cheatah::space::irbem
