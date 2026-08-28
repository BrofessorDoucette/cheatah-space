#pragma once

/**
 * @file civil.hpp
 * @brief space.time — the proleptic Gregorian calendar, in exact integer arithmetic.
 *
 * The calendar layer the CDF epoch types are built on. `import space.time` resolves it
 * through `time.hpp`; it is separated only because it is a self-contained, exhaustively
 * testable unit that both `tt2000.hpp` and `cdfepoch.hpp` consume.
 *
 * Everything here is `constexpr` integer arithmetic — no floating point anywhere. That is
 * the point: `time.hpp`'s Julian Date conversions carry a `double`, which near the present
 * epoch resolves to roughly 40 microseconds, and CDF_TIME_TT2000 needs exact nanoseconds.
 * A calendar that rounds is a calendar that puts an event on the wrong day once a decade.
 *
 * Years use astronomical numbering: year 0 is 1 BC, year -1 is 2 BC. The calendar is
 * PROLEPTIC — Gregorian rules are projected back before 1582 rather than switching to
 * Julian — because that is what CDF does, and agreeing with the format matters more here
 * than agreeing with history.
 *
 * @note Two independent day-number algorithms live here, and that is deliberate rather than
 *       redundant. days_from_civil() is the correct one. nasa_julian_day_at_noon() is a
 *       literal transcription of the formula NASA's CDF library uses. They are cross-checked
 *       against each other across the whole TT2000 range in the unit tests, which is what
 *       lets us claim our epoch conversions agree with NASA's *because we checked*, rather
 *       than because both look plausible.
 *
 * Cross-platform, header-only, allocation-free: no platform headers, no I/O, no global state.
 */

#include "cheatah.hpp"

namespace cheatah::space::time {

/**
 * A calendar date in the proleptic Gregorian calendar, astronomical year numbering.
 *
 * A plain aggregate: no invariant is enforced on construction, so a `CivilDate` can hold
 * 2026-02-30. Ask is_valid_civil() when the value came from outside; the conversions here
 * are total functions over whatever they are handed.
 */
struct CivilDate {
    int year{};   ///< Astronomical year numbering: 0 is 1 BC, -1 is 2 BC.
    int month{};  ///< 1–12 for a valid date.
    int day{};    ///< 1–31 for a valid date, bounded by days_in_month().
};

/**
 * Whether @p year is a leap year under proleptic Gregorian rules.
 * @param year astronomical year number.
 * @return true when @p year has 366 days.
 * @complexity O(1).
 * @alloc none.
 * @systest systests/test_civil.purr
 * @par Example
 * @code{.purr}
 * import io
 * import space.time as st
 *
 * io.print(st.is_leap_year(2000))   # True  — divisible by 400
 * io.print(st.is_leap_year(1900))   # False — divisible by 100 but not 400
 * @endcode
 */
constexpr bool is_leap_year(int year) noexcept {
    return (year % 4 == 0 && year % 100 != 0) || year % 400 == 0;
}

/**
 * The number of days in a given month.
 * @param year astronomical year number (decides February's length).
 * @param month 1–12; anything outside that range yields 0.
 * @return days in the month, or 0 when @p month is out of range.
 * @complexity O(1).
 * @alloc none.
 * @systest systests/test_civil.purr
 * @par Example
 * @code{.purr}
 * import io
 * import space.time as st
 *
 * io.print(st.days_in_month(2024, 2))   # 29
 * io.print(st.days_in_month(2026, 2))   # 28
 * @endcode
 */
constexpr int days_in_month(int year, int month) noexcept {
    if (month < 1 || month > 12) { return 0; }
    // Indexed 1–12; entry 0 is unused padding so the month number indexes directly.
    constexpr int lengths[13] = {0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
    return (month == 2 && is_leap_year(year)) ? 29 : lengths[month];
}

/**
 * Whether @p date names a day that exists.
 * @param date the date to check.
 * @return true when the month is 1–12 and the day is within that month's length.
 * @complexity O(1).
 * @alloc none.
 * @systest systests/test_civil.purr
 * @par Example
 * @code{.purr}
 * import io
 * import space.time as st
 *
 * io.print(st.is_valid_civil(st.civil_date(2024, 2, 29)))   # True  — 2024 is a leap year
 * io.print(st.is_valid_civil(st.civil_date(2026, 2, 29)))   # False — 2026 is not
 * @endcode
 */
constexpr bool is_valid_civil(const CivilDate& date) noexcept {
    return date.month >= 1 && date.month <= 12 && date.day >= 1
           && date.day <= days_in_month(date.year, date.month);
}

/**
 * Build a CivilDate from its three components.
 *
 * Exists because a `.purr` caller cannot write a C++ aggregate initializer; from C++ prefer
 * `CivilDate{y, m, d}` directly. No validation is performed — see is_valid_civil().
 *
 * @param year astronomical year number.
 * @param month 1–12 for a valid date.
 * @param day 1–31 for a valid date.
 * @return the assembled date.
 * @complexity O(1).
 * @alloc none.
 * @systest systests/test_civil.purr
 * @par Example
 * @code{.purr}
 * import io
 * import space.time as st
 *
 * let d = st.civil_date(2026, 8, 28)
 * io.print(d.year, d.month, d.day)   # 2026 8 28
 * @endcode
 */
constexpr CivilDate civil_date(int year, int month, int day) noexcept {
    return CivilDate{year, month, day};
}

/**
 * Days since the Unix epoch (1970-01-01) for a proleptic Gregorian date.
 *
 * Howard Hinnant's `days_from_civil`, which is exact for every year representable in `int`
 * and is the inverse of civil_from_days(). Shifting the year to start in March is what makes
 * the leap day fall at the end of the year, so no special case for February is needed.
 *
 * @param year astronomical year number.
 * @param month 1–12.
 * @param day 1–31.
 * @return days before (negative) or after (positive) 1970-01-01.
 * @complexity O(1).
 * @alloc none.
 * @systest systests/test_civil.purr
 * @par Example
 * @code{.purr}
 * import io
 * import space.time as st
 *
 * io.print(st.days_from_civil(1970, 1, 1))   # 0
 * io.print(st.days_from_civil(2000, 1, 1))   # 10957 — the J2000 day
 * @endcode
 */
constexpr long long days_from_civil(int year, int month, int day) noexcept {
    long long y = year;
    y -= (month <= 2) ? 1 : 0;                       // years start in March
    const long long era = (y >= 0 ? y : y - 399) / 400;
    const long long yoe = y - era * 400;             // year of era, [0, 399]
    const long long mp = (month + (month > 2 ? -3 : 9));  // March-based month, [0, 11]
    const long long doy = (153 * mp + 2) / 5 + day - 1;   // day of year, [0, 365]
    const long long doe = yoe * 365 + yoe / 4 - yoe / 100 + doy;  // day of era, [0, 146096]
    return era * 146097 + doe - 719468;              // shift the era origin to 1970-01-01
}

/**
 * The proleptic Gregorian date for a count of days since the Unix epoch.
 *
 * The exact inverse of days_from_civil() over the whole representable range.
 *
 * @param days days before (negative) or after (positive) 1970-01-01.
 * @return the calendar date.
 * @complexity O(1).
 * @alloc none.
 * @systest systests/test_civil.purr
 * @par Example
 * @code{.purr}
 * import io
 * import space.time as st
 *
 * let d = st.civil_from_days(10957)
 * io.print(d.year, d.month, d.day)   # 2000 1 1
 * @endcode
 */
constexpr CivilDate civil_from_days(long long days) noexcept {
    const long long z = days + 719468;
    const long long era = (z >= 0 ? z : z - 146096) / 146097;
    const long long doe = z - era * 146097;                                  // [0, 146096]
    const long long yoe = (doe - doe / 1460 + doe / 36524 - doe / 146096) / 365;  // [0, 399]
    const long long doy = doe - (365 * yoe + yoe / 4 - yoe / 100);           // [0, 365]
    const long long mp = (5 * doy + 2) / 153;                                // [0, 11]
    const long long day = doy - (153 * mp + 2) / 5 + 1;                      // [1, 31]
    const long long month = mp + (mp < 10 ? 3 : -9);                         // [1, 12]
    return CivilDate{static_cast<int>(yoe + era * 400 + (month <= 2 ? 1 : 0)),
                     static_cast<int>(month), static_cast<int>(day)};
}

/**
 * The Julian Day number at noon, computed exactly the way NASA's CDF library computes it.
 *
 * A literal transcription of the `_JulianDay` formula in NASA's CDF distribution, truncating
 * integer division and all. It is here so our CDF epoch conversions can be shown to agree
 * with NASA's rather than merely be believed to: the unit tests assert this equals
 * `days_from_civil(...) + 2440588` for every day in the TT2000 range, which retires the whole
 * class of off-by-one calendar bugs in one loop.
 *
 * @warning It is EXACT for every year >= 0 and WRONG for negative years — off by one day at
 *          year -1, by two by year -401, and drifting further back. The formula needs floor
 *          division and C truncates toward zero, which only differs once the numerator goes
 *          negative. That costs CDF nothing, and the bound is the reason: CDF_EPOCH counts from
 *          year 0 and CDF_TIME_TT2000 reaches back only to 1707, so no CDF-representable instant
 *          falls in the broken region. Verified day by day in tests/space_civil_test.cpp, both
 *          that they agree for year >= 0 and that they diverge below it — the divergence is
 *          pinned deliberately, so a future "fix" to this transcription fails the build.
 *
 * Prefer days_from_civil() for new code. This exists for provenance, not for use.
 *
 * @param year astronomical year number.
 * @param month 1–12.
 * @param day 1–31.
 * @return the Julian Day number of the noon that begins this date's Julian day.
 * @complexity O(1).
 * @alloc none.
 * @systest systests/test_civil.purr
 * @par Example
 * @code{.purr}
 * import io
 * import space.time as st
 *
 * # J2000: 2000-01-01 has Julian Day 2451545 at noon.
 * io.print(st.nasa_julian_day_at_noon(2000, 1, 1))   # 2451545
 * @endcode
 */
constexpr long long nasa_julian_day_at_noon(int year, int month, int day) noexcept {
    const long long y = year;
    const long long m = month;
    return 367 * y - (7 * (y + ((m + 9) / 12))) / 4
           - (3 * (((y + ((m - 9) / 7)) / 100) + 1)) / 4
           + (275 * m) / 9 + day + 1721029;
}

}  // namespace cheatah::space::time
