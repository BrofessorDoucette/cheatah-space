#pragma once

/**
 * @file calendar.hpp
 * @brief space.time — Julian day numbers, day-of-year, and decimal years, in exact arithmetic.
 *
 * The conversions between a calendar date and the numeric time scales physics actually computes
 * in: the Julian day number, the day of year, and the decimal year that every geomagnetic field
 * model takes as its epoch argument.
 *
 * **These live in `space.time`, not in a consumer.** They are generic calendar arithmetic with no
 * space-physics content — `space.irbem` needs them, `space.cdf` needs them, and anything that
 * timestamps an observation needs them. They were briefly written inside `space.irbem`, which put
 * the calendar in a module named after radiation belts and duplicated the leap rule against
 * @ref civil.hpp.
 *
 * @ref civil.hpp is the neighbouring layer and the split is deliberate: `civil.hpp` is the
 * proleptic Gregorian calendar in exact INTEGER arithmetic, because CDF_TIME_TT2000 needs exact
 * nanoseconds. This header sits above it and adds the scales that carry a fraction — a decimal year
 * is a `double` by nature, being a position within a year of varying length.
 *
 * @note @ref cheatah::space::time::is_leap_year is NOT redefined here; it is @ref civil.hpp's.
 *       Two definitions of the leap rule in one namespace is an ODR violation waiting to fire, and
 *       worse, two places for the century rule to be got wrong independently.
 *
 * The Julian day arithmetic is Fliegel & Van Flandern (1968), *A machine algorithm for processing
 * calendar dates*, Comm. ACM 11(10):657 — integer throughout, so there is no floating-point
 * intermediate and therefore no rounding, valid across the Gregorian reform and into negative
 * years under astronomical numbering.
 */

#include <cstdint>
#include <string_view>

#include "cheatah.hpp"
#include "civil.hpp"

namespace cheatah::space::time {

/// A calendar date, with no time of day.
struct CalendarDate {
    int year = 0;   ///< The year; astronomical numbering, so year 0 is 1 BC.
    int month = 1;  ///< Month, 1 = January … 12 = December.
    int day = 1;    ///< Day of the month, 1-based.

    /// Two dates are equal when all three fields are.
    /// @param a the left date. @param b the right date.
    /// @return whether they name the same day.
    /// @complexity O(1). @alloc none.
    friend constexpr bool operator==(const CalendarDate& a, const CalendarDate& b) = default;
};

/// A calendar date with a UT time of day — what IRBEM's date routines hand back.
struct DateTime {
    int year = 0;         ///< The year; astronomical numbering.
    int month = 1;        ///< Month, 1 = January … 12 = December.
    int day = 1;          ///< Day of the month, 1-based.
    int day_of_year = 1;  ///< Day of year, 1 = January 1.
    int hour = 0;         ///< UT hour of day, 0…23.
    int minute = 0;       ///< UT minute, 0…59.
    int second = 0;       ///< UT second, 0…59, truncated; @ref ut_seconds keeps the remainder.
    double ut_seconds = 0.0;  ///< UT time of day in seconds since midnight, `[0, 86400)`.

    /// Field-wise equality.
    /// @param a the left value. @param b the right value.
    /// @return whether every field matches.
    /// @complexity O(1). @alloc none.
    friend constexpr bool operator==(const DateTime& a, const DateTime& b) = default;
};

namespace detail {

/// Microseconds in a day — the unit the sub-day arithmetic carries, chosen so a day is an exact
/// integer count and no time of day can drift.
inline constexpr std::int64_t microseconds_per_day = 86'400'000'000;

/// Microseconds in a second.
inline constexpr std::int64_t microseconds_per_second = 1'000'000;

/// The rounding grain of a decimal-year round trip: a decimal year is a `double`, so it cannot
/// carry a full year at microsecond resolution. Round trips land within a millisecond — far below
/// any geomagnetic model's time resolution, and asserted rather than assumed.
inline constexpr std::int64_t decimal_year_resolution_us = 1'000;

/// Floor division for signed integers — C++ truncates toward zero, which is the wrong direction
/// for calendar arithmetic before year 0 and produces an off-by-one day there.
/// @param a the dividend. @param b the divisor; must be positive.
/// @return `floor(a / b)`.
/// @complexity O(1). @alloc none.
constexpr std::int64_t floor_div(std::int64_t a, std::int64_t b) {
    const std::int64_t q = a / b;
    return ((a % b != 0) && ((a < 0) != (b < 0))) ? q - 1 : q;
}

/// `floor` for a double, as an integer, without going through `<cmath>` so it stays `constexpr`.
/// @param x the value. @return the largest integer not greater than @p x.
/// @complexity O(1). @alloc none.
constexpr std::int64_t floor_to_int(double x) {
    const auto t = static_cast<std::int64_t>(x);
    return (x < 0.0 && static_cast<double>(t) != x) ? t - 1 : t;
}

/// Round-half-away-from-zero, as an integer, `constexpr`.
/// @param x the value. @return the nearest integer.
/// @complexity O(1). @alloc none.
constexpr std::int64_t round_to_int(double x) {
    return x >= 0.0 ? floor_to_int(x + 0.5) : -floor_to_int(-x + 0.5);
}

}  // namespace detail

/// The Julian day number at which the Gregorian calendar takes over: 1582-10-15, the day Pope
/// Gregory XIII's reform declared should follow 1582-10-04.
inline constexpr std::int64_t gregorian_reform_jdn = 2299161;

/**
 * The Julian day number of a calendar date.
 *
 * **Reform-aware, not proleptic.** Dates on or after the Gregorian reform use the Gregorian rule;
 * earlier dates use the JULIAN rule, which is what they were actually recorded in. This differs
 * deliberately from @ref civil.hpp, which is proleptic Gregorian throughout because that is what
 * the CDF format specifies — the two conventions disagree by ten days at the reform and by a
 * growing amount before it, so they are separate functions rather than one with a flag.
 *
 * Which convention a caller wants depends on what the date MEANS: a CDF timestamp is proleptic by
 * definition, whereas a historical observation is in the calendar its observer used. Geomagnetic
 * work reaches back to the 19th century and IGRF to 1900, both comfortably post-reform, so in
 * practice the two agree everywhere this module is used — but agreeing by accident is not the same
 * as agreeing by construction, and a silent ten-day error is exactly the kind that survives review.
 *
 * Integer arithmetic throughout (Fliegel & Van Flandern 1968 for the Gregorian branch), so there is
 * no floating-point intermediate and no rounding, valid into negative years under astronomical
 * numbering: JDN 0 is -4712-01-01 in the Julian calendar.
 *
 * @param year the year, astronomical numbering. @param month 1-12. @param day 1-based.
 * @return the Julian day number — a count of days where the DAY BEGINS AT NOON. A caller wanting a
 *         Julian DATE for a sidereal-time series must subtract 0.5 and add the UT fraction; getting
 *         that half-day wrong shifts GMST by twelve hours, which rotates GSM by 180 degrees and
 *         yields a field that looks plausible and is inverted.
 * @complexity O(1). @alloc none.
 */
constexpr std::int64_t julian_day_number(int year, int month, int day) {
    const std::int64_t y = year;
    const std::int64_t m = month;
    const std::int64_t d = day;

    // The Gregorian rule, evaluated first so the reform can be tested against its own result.
    const std::int64_t a = detail::floor_div(14 - m, 12);
    const std::int64_t yy = y + 4800 - a;
    const std::int64_t mm = m + (12 * a) - 3;
    const std::int64_t greg = d + detail::floor_div((153 * mm) + 2, 5) + (365 * yy) +
                              detail::floor_div(yy, 4) - detail::floor_div(yy, 100) +
                              detail::floor_div(yy, 400) - 32045;
    if (greg >= gregorian_reform_jdn) return greg;

    // Before the reform, the Julian rule — no century exception, so a leap year every fourth.
    // `(m - 9) / 7` is TRUNCATING division on purpose, not floor: the classic form assumes C
    // semantics here, and `m` is 1-12 so the numerator is -8…3 and never ambiguous. Substituting
    // floor_div (as the Gregorian branch correctly uses elsewhere) shifts JDN 0 by two days.
    const std::int64_t jm = (m - 9) / 7;
    return (367 * y) - detail::floor_div(7 * (y + 5001 + jm), 4) + detail::floor_div(275 * m, 9) +
           d + 1729777;
}

/**
 * The calendar date of a Julian day number — the exact inverse of @ref julian_day_number,
 * including across the reform.
 * @param jdn the Julian day number.
 * @return the calendar date, in whichever calendar was in force at that JDN.
 * @complexity O(1). @alloc none.
 */
constexpr CalendarDate calendar_date(std::int64_t jdn) {
    if (jdn >= gregorian_reform_jdn) {
        const std::int64_t a = jdn + 32044;
        const std::int64_t b = detail::floor_div((4 * a) + 3, 146097);
        const std::int64_t c = a - detail::floor_div(146097 * b, 4);
        const std::int64_t d = detail::floor_div((4 * c) + 3, 1461);
        const std::int64_t e = c - detail::floor_div(1461 * d, 4);
        const std::int64_t m = detail::floor_div((5 * e) + 2, 153);
        return CalendarDate{static_cast<int>((100 * b) + d - 4800 + detail::floor_div(m, 10)),
                            static_cast<int>(m + 3 - (12 * detail::floor_div(m, 10))),
                            static_cast<int>(e - detail::floor_div((153 * m) + 2, 5) + 1)};
    }
    const std::int64_t j = jdn + 1402;
    const std::int64_t k = detail::floor_div(j - 1, 1461);
    const std::int64_t l = j - (1461 * k);
    const std::int64_t n = detail::floor_div(l - 1, 365) - detail::floor_div(l, 1461);
    const std::int64_t i = l - (365 * n) + 30;
    const std::int64_t jj = detail::floor_div(80 * i, 2447);
    const std::int64_t dd = i - detail::floor_div(2447 * jj, 80);
    const std::int64_t ii = detail::floor_div(jj, 11);
    return CalendarDate{static_cast<int>((4 * k) + n + ii - 4716),
                        static_cast<int>(jj + 2 - (12 * ii)), static_cast<int>(dd)};
}

/**
 * The day of year of a calendar date, leap-year correct.
 * @param year the year. @param month 1-12. @param day 1-based.
 * @return the day of year, 1 = January 1.
 * @complexity O(1). @alloc none.
 */
constexpr int day_of_year(int year, int month, int day) {
    return static_cast<int>(julian_day_number(year, month, day) -
                            julian_day_number(year, 1, 1)) +
           1;
}

/**
 * How many days a year has.
 * @param year the year.
 * @return 366 in a leap year, 365 otherwise — computed from the calendar rather than from the leap
 *         rule directly, so the two can never disagree.
 * @complexity O(1). @alloc none.
 */
constexpr int days_in_year(int year) {
    return static_cast<int>(julian_day_number(year + 1, 1, 1) - julian_day_number(year, 1, 1));
}

/**
 * The calendar date of a given day of year — the inverse of @ref day_of_year, and the piece
 * IRBEM's `DOY_AND_UT2DATE_AND_TIME` needs.
 *
 * @param year the year, astronomical numbering.
 * @param doy the day of year, 1 = January 1. Values outside `[1, days_in_year(year)]` are not
 *            rejected; they roll into the neighbouring year, which is what makes a UT offset that
 *            crosses midnight on December 31 come out right instead of producing a "day 367".
 * @return the calendar date.
 * @complexity O(1).
 * @alloc none.
 */
constexpr CalendarDate date_from_day_of_year(int year, int doy) {
    return calendar_date(julian_day_number(year, 1, 1) + doy - 1);
}

namespace detail {

/// Assemble a @ref DateTime from a Julian Day Number and a microsecond offset into that day.
/// Shared by the decimal-year and day-of-year entry points so the two cannot drift apart in how
/// they split a day into h/m/s.
/// @param jdn the Julian Day Number of the day.
/// @param us_of_day microseconds since midnight UT, in `[0, microseconds_per_day)`.
/// @return the fully populated broken-down instant.
/// @complexity O(1).
/// @alloc none.
constexpr DateTime assemble(std::int64_t jdn, std::int64_t us_of_day) {
    const CalendarDate date = calendar_date(jdn);
    const std::int64_t whole_seconds = us_of_day / microseconds_per_second;

    DateTime out{};
    out.year = date.year;
    out.month = date.month;
    out.day = date.day;
    out.day_of_year = static_cast<int>(jdn - julian_day_number(date.year, 1, 1) + 1);
    out.hour = static_cast<int>(whole_seconds / 3600);
    out.minute = static_cast<int>((whole_seconds / 60) % 60);
    out.second = static_cast<int>(whole_seconds % 60);
    out.ut_seconds =
        static_cast<double>(us_of_day) / static_cast<double>(microseconds_per_second);
    return out;
}

}  // namespace detail

/**
 * Date and time to decimal year — IRBEM's `DATE_AND_TIME2DECY`.
 *
 * `yyyy.0` is January 1 at 00:00 UT and the year's own length is the denominator, so the fraction
 * is the elapsed portion of *that* year: a mid-year instant is `.5` in a leap year and in a common
 * year alike. This is the form IGRF coefficient interpolation consumes.
 *
 * Everything above the final division is integer, so the only rounding in the whole routine is the
 * one unavoidable division — and when the fraction is exactly representable the result is exact.
 * `decimal_year(2000, 7, 2, 0, 0, 0)` is `2000.5` bit-for-bit, because 183 of 366 days is one half.
 *
 * @param year the year, astronomical numbering.
 * @param month the month, 1 = January … 12 = December.
 * @param day the day of the month, 1-based.
 * @param hour the UT hour of day, 0…23.
 * @param minute the UT minute, 0…59.
 * @param second the UT second, 0…59.
 * @return the decimal year.
 * @complexity O(1).
 * @alloc none.
 */
constexpr double decimal_year(int year, int month, int day, int hour, int minute, int second) {
    const std::int64_t jan1 = julian_day_number(year, 1, 1);
    const std::int64_t elapsed_seconds =
        (julian_day_number(year, month, day) - jan1) * 86400 +
        static_cast<std::int64_t>(hour) * 3600 + static_cast<std::int64_t>(minute) * 60 +
        static_cast<std::int64_t>(second);
    const std::int64_t year_seconds = (julian_day_number(year + 1, 1, 1) - jan1) * 86400;

    return static_cast<double>(year) +
           static_cast<double>(elapsed_seconds) / static_cast<double>(year_seconds);
}

/**
 * Decimal year back to a broken-down instant — IRBEM's `DECY2DATE_AND_TIME`, and the inverse of
 * @ref decimal_year.
 *
 * The fraction is converted to an integer count of microseconds of year and then snapped to the
 * millisecond before anything is split off it — see `detail::decimal_year_resolution_us` for the
 * derivation. That snap is the whole point of this routine: a decimal year near 2000 resolves to
 * only about 7 microseconds, so a whole hour recovered from one lands a few microseconds short and
 * a plain truncation reports the previous second. Measured against the shipped IRBEM library, that
 * defect fires on roughly one round trip in six across a sweep of the year 2000; this routine
 * returns the instant it was given, exactly, for every input @ref decimal_year can produce.
 *
 * @param decy the decimal year, `yyyy.0` being January 1 at 00:00 UT.
 * @return the instant, with month, day, day of year, h:m:s and seconds-of-day all filled in.
 * @complexity O(1).
 * @alloc none.
 */
constexpr DateTime date_and_time_from_decimal_year(double decy) {
    std::int64_t year = detail::floor_to_int(decy);
    const double fraction = decy - static_cast<double>(year);

    const std::int64_t year_us =
        static_cast<std::int64_t>(days_in_year(static_cast<int>(year))) *
        detail::microseconds_per_day;
    // Snap to the resolution a decimal year can actually carry; see decimal_year_resolution_us
    // for why doing this in microseconds instead silently loses a second.
    const std::int64_t raw_us = detail::round_to_int(fraction * static_cast<double>(year_us));
    std::int64_t us = ((raw_us + detail::decimal_year_resolution_us / 2) /
                       detail::decimal_year_resolution_us) *
                      detail::decimal_year_resolution_us;

    // The snap can land exactly on the far end of the year — the instant that IS the next
    // January 1. One step is all that is ever possible, since the fraction is below one by
    // construction.
    if (us >= year_us) {
        us = 0;
        year += 1;
    }

    const std::int64_t day_offset = us / detail::microseconds_per_day;
    return detail::assemble(julian_day_number(static_cast<int>(year), 1, 1) + day_offset,
                            us - day_offset * detail::microseconds_per_day);
}

/**
 * Year, day of year and UT seconds to a broken-down instant — IRBEM's
 * `DOY_AND_UT2DATE_AND_TIME`.
 *
 * @param year the year, astronomical numbering.
 * @param doy the day of year, 1 = January 1.
 * @param ut_seconds UT time of day in seconds since midnight. Values outside `[0, 86400)` are
 *                   carried into the day count rather than rejected, so a caller stepping a
 *                   trajectory by adding seconds gets the right date across midnight — including
 *                   backwards, where a negative offset walks into the previous day.
 * @return the instant. `day_of_year` in the result is recomputed from the resulting date, so it
 *         reflects any such carry rather than echoing @p doy back.
 * @complexity O(1).
 * @alloc none.
 */
constexpr DateTime date_and_time_from_doy_and_ut(int year, int doy, double ut_seconds) {
    const std::int64_t us = detail::round_to_int(
        ut_seconds * static_cast<double>(detail::microseconds_per_second));
    const std::int64_t day_offset = detail::floor_div(us, detail::microseconds_per_day);

    return detail::assemble(julian_day_number(year, 1, 1) + doy - 1 + day_offset,
                            us - day_offset * detail::microseconds_per_day);
}

}  // namespace cheatah::space::time
