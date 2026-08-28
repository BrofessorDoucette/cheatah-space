#pragma once

/**
 * @file leapseconds.hpp
 * @brief space.cdf — NASA's leap-second table, as exact integers.
 *
 * TAI-UTC: the quantity CDF_TIME_TT2000 needs, and the only reason TT2000 is harder than
 * counting nanoseconds. `import space.cdf` resolves it through `cdf.hpp`.
 *
 * The table is hand-authored here rather than read from a file at run time, and that is a
 * deliberate trade. A header-only numeric library that reaches for `getenv` or opens a file
 * during a conversion has global state, and the same CDF then decodes to different timestamps
 * on two machines — a reproducibility failure far worse for science data than being six months
 * behind on a leap second. Currency is handled OUT of band by `scripts/check_leapseconds.sh`,
 * which re-fetches NASA's table and fails when it no longer matches leap_seconds_sha256().
 *
 * ARITHMETIC. Every column is an exact scaled integer, so no `double` ever touches the
 * leap-second path:
 *
 *   - `delta_at_e7`    TAI-UTC in units of 1e-7 s   (column 4 x 1e7)
 *   - `drift_mjd`      the MJD the drift references (column 5; integral, 0 from 1972)
 *   - `drift_rate_e7`  drift in 1e-7 s/day          (column 6 x 1e7; 0 from 1972)
 *
 * All 42 rows are exactly representable that way — NASA prints at most 7 decimals — which is
 * what makes tai_minus_utc_ns() exact rather than merely close. See its docs for where that
 * makes us differ from NASA's own answer, and by exactly how much.
 *
 * THE PRE-1972 ERA. The first 14 rows are the "rubber second" era, when UTC ran at a rate
 * slightly different from TAI and the offset drifted CONTINUOUSLY with the day. Those rows are
 * why the drift columns exist; from 1972-01-01 the offset is a whole number of seconds.
 *
 * NEGATIVE LEAP SECONDS ARE NOT HYPOTHETICAL. Two rows in this very table DECREASE —
 * 1961-08-01 and 1968-02-01 — so code assuming the offset only ever grows is already wrong on
 * real data. It matters going forward too: there has been no leap second since 2017-01-01,
 * Earth's rotation has been running fast, and the next one may well be negative.
 *
 * Cross-platform, header-only, allocation-free: no platform headers, no I/O, no global state.
 */

#include "cheatah.hpp"

#include <cstddef>
#include <cstdint>
#include <string_view>

// The calendar. space.cdf depends on space.time; never the reverse.
//
// The drift term below needs the Julian Day Number at noon, and it deliberately uses
// nasa_julian_day_at_noon() — the literal transcription of the formula NASA's own CDF library
// uses — rather than a mathematically cleaner one. That makes our leap-second arithmetic
// bit-compatible with the reference implementation by construction rather than by coincidence.
// The transcription is wrong for NEGATIVE years (it needs floor division and C truncates), which
// costs nothing here: CDF_EPOCH counts from year 0 and TT2000 reaches back only to 1707, so no
// CDF-representable date falls in the broken region. Asserted day by day in the unit tests.
#include "../time/civil.hpp"

namespace cheatah::space::cdf {

/**
 * One row of NASA's `CDFLeapSeconds.txt`, held as exact scaled integers.
 *
 * A row takes effect at 00:00:00 UTC on its date and holds until the next row's date.
 */
struct LeapSecondEntry {
    int year{};                 ///< UTC year the entry takes effect.
    int month{};                ///< UTC month the entry takes effect.
    int day{};                  ///< UTC day the entry takes effect.
    long long delta_at_e7{};    ///< TAI-UTC in units of 1e-7 s.
    long long drift_mjd{};      ///< MJD the drift is referenced to; 0 from 1972 on.
    long long drift_rate_e7{};  ///< Drift in 1e-7 s per day; 0 from 1972 on.
};

namespace detail {
/// Rows in the shipped table.
inline constexpr std::size_t kLeapSecondCount = 42;
/// Rows in the pre-1972 "rubber second" era. Matches `NERA1` in NASA's implementation.
inline constexpr std::size_t kLeapSecondDriftRows = 14;
/// The table itself, ascending by date. Reached publicly through leap_second_entry().
inline constexpr LeapSecondEntry kLeapSeconds[kLeapSecondCount] = {
    {1960,  1,  1,     14178180,  37300,    12960},
    {1961,  1,  1,     14228180,  37300,    12960},
    {1961,  8,  1,     13728180,  37300,    12960},
    {1962,  1,  1,     18458580,  37665,    11232},
    {1963, 11,  1,     19458580,  37665,    11232},
    {1964,  1,  1,     32401300,  38761,    12960},
    {1964,  4,  1,     33401300,  38761,    12960},
    {1964,  9,  1,     34401300,  38761,    12960},
    {1965,  1,  1,     35401300,  38761,    12960},
    {1965,  3,  1,     36401300,  38761,    12960},
    {1965,  7,  1,     37401300,  38761,    12960},
    {1965,  9,  1,     38401300,  38761,    12960},
    {1966,  1,  1,     43131700,  39126,    25920},
    {1968,  2,  1,     42131700,  39126,    25920},
    {1972,  1,  1,    100000000,      0,        0},
    {1972,  7,  1,    110000000,      0,        0},
    {1973,  1,  1,    120000000,      0,        0},
    {1974,  1,  1,    130000000,      0,        0},
    {1975,  1,  1,    140000000,      0,        0},
    {1976,  1,  1,    150000000,      0,        0},
    {1977,  1,  1,    160000000,      0,        0},
    {1978,  1,  1,    170000000,      0,        0},
    {1979,  1,  1,    180000000,      0,        0},
    {1980,  1,  1,    190000000,      0,        0},
    {1981,  7,  1,    200000000,      0,        0},
    {1982,  7,  1,    210000000,      0,        0},
    {1983,  7,  1,    220000000,      0,        0},
    {1985,  7,  1,    230000000,      0,        0},
    {1988,  1,  1,    240000000,      0,        0},
    {1990,  1,  1,    250000000,      0,        0},
    {1991,  1,  1,    260000000,      0,        0},
    {1992,  7,  1,    270000000,      0,        0},
    {1993,  7,  1,    280000000,      0,        0},
    {1994,  7,  1,    290000000,      0,        0},
    {1996,  1,  1,    300000000,      0,        0},
    {1997,  7,  1,    310000000,      0,        0},
    {1999,  1,  1,    320000000,      0,        0},
    {2006,  1,  1,    330000000,      0,        0},
    {2009,  1,  1,    340000000,      0,        0},
    {2012,  7,  1,    350000000,      0,        0},
    {2015,  7,  1,    360000000,      0,        0},
    {2017,  1,  1,    370000000,      0,        0},
};

/// A date packed so that ordering by integer value orders by date. Cheaper than a three-way
/// lexicographic compare, and it has no sub-condition that valid input cannot reach — the
/// obvious `y > ey || (y == ey && (m > em || ...))` form contains a `day >= entry.day` test that
/// can never be false, because every entry in the table falls on the first of a month.
constexpr long long packed_date(int year, int month, int day) noexcept {
    return (static_cast<long long>(year) * 10000) + (month * 100) + day;
}

/// Index of the last entry on or before the given date, or -1 when it precedes the table
/// (before 1960-01-01, where UTC as we now define it did not yet exist).
constexpr int leap_index_for(int year, int month, int day) noexcept {
    const long long want = packed_date(year, month, day);
    int found = -1;
    for (std::size_t i = 0; i < kLeapSecondCount; ++i) {
        const LeapSecondEntry& e = kLeapSeconds[i];
        if (packed_date(e.year, e.month, e.day) > want) { break; }
        found = static_cast<int>(i);
    }
    return found;
}
}  // namespace detail

/**
 * The number of entries in the leap-second table.
 * @return 42 as shipped; grows only when a leap second is added upstream.
 * @complexity O(1).
 * @alloc none.
 * @systest systests/test_leapseconds.purr
 * @par Example
 * @code{.purr}
 * import io
 * import space.cdf as cdf
 *
 * io.print(cdf.leap_second_count())   # 42
 * @endcode
 */
constexpr long long leap_second_count() noexcept {
    return static_cast<long long>(detail::kLeapSecondCount);
}

/**
 * How many leading entries lie in the pre-1972 "rubber second" era, where TAI-UTC drifts
 * continuously with the day rather than stepping by whole seconds.
 * @return 14, matching `NERA1` in NASA's own implementation.
 * @complexity O(1).
 * @alloc none.
 * @systest systests/test_leapseconds.purr
 * @par Example
 * @code{.purr}
 * import io
 * import space.cdf as cdf
 *
 * io.print(cdf.leap_second_drift_rows())   # 14 — entries 0..13 drift
 * @endcode
 */
constexpr long long leap_second_drift_rows() noexcept {
    return static_cast<long long>(detail::kLeapSecondDriftRows);
}

/**
 * One entry of the leap-second table, by index.
 * @param index 0-based; anything outside [0, leap_second_count()) yields a zeroed entry.
 * @return the entry, by value.
 * @complexity O(1).
 * @alloc none.
 * @systest systests/test_leapseconds.purr
 * @par Example
 * @code{.purr}
 * import io
 * import space.cdf as cdf
 *
 * # The most recent entry: TAI-UTC became 37 s on 2017-01-01.
 * let e = cdf.leap_second_entry(41)
 * io.print(e.year, e.month, e.day)   # 2017 1 1
 * io.print(e.delta_at_e7)            # 370000000 — 37.0 s, in units of 1e-7 s
 * @endcode
 */
constexpr LeapSecondEntry leap_second_entry(long long index) noexcept {
    if (index < 0 || index >= static_cast<long long>(detail::kLeapSecondCount)) {
        return LeapSecondEntry{};
    }
    return detail::kLeapSeconds[static_cast<std::size_t>(index)];
}

/**
 * TAI-UTC in EXACT nanoseconds for a UTC calendar date.
 *
 * Pure integer arithmetic. Before 1972 the offset drifts with the day, which is what the drift
 * columns encode; the MJD at noon is a half-integer, so the drift term is evaluated at twice
 * scale and halved by folding the 2 into the 1e-7 -> 1e-9 conversion. Nothing rounds.
 *
 * @warning This is EXACT, and NASA's library is not. NASA computes the same quantity in
 *          `double` and truncates toward zero, so on 430 of the 4384 days between 1960-01-01
 *          and 1972-01-01 its answer is 1 ns lower than this one. From 1972-01-01 onward the two
 *          agree on every day, because the offset is a whole number of seconds and the drift
 *          term vanishes. Use tai_minus_utc_ns_nasa_compat() when bit-identical agreement with
 *          NASA matters more than being right; the divergence is enumerated day by day in
 *          tests/space_cdf_leapseconds_test.cpp so it stays a reviewed fact, not a surprise.
 *
 * @param year UTC year.
 * @param month UTC month, 1-12.
 * @param day UTC day of month.
 * @return TAI-UTC in nanoseconds; 0 before 1960-01-01, where the table does not reach.
 * @complexity O(leap_second_count()) — a linear scan of 42 entries.
 * @alloc none.
 * @systest systests/test_leapseconds.purr
 * @par Example
 * @code{.purr}
 * import io
 * import space.cdf as cdf
 *
 * # Since 2017-01-01, TAI runs 37 seconds ahead of UTC.
 * io.print(cdf.tai_minus_utc_ns(2026, 8, 28))   # 37000000000
 * @endcode
 */
constexpr long long tai_minus_utc_ns(int year, int month, int day) noexcept {
    const int i = detail::leap_index_for(year, month, day);
    if (i < 0) { return 0; }
    const LeapSecondEntry& e = detail::kLeapSeconds[i];
    if (e.drift_rate_e7 == 0) { return e.delta_at_e7 * 100; }
    // MJD = JD_at_noon - 2400000.5, a half-integer. Work at twice scale so it stays an integer:
    // two_mjd = 2*MJD. The drift term is (MJD - drift_mjd) * rate_e7 * 100 nanoseconds, i.e.
    // (two_mjd - 2*drift_mjd) * rate_e7 * 50 — exact, and at most ~1e11, so int64 holds it.
    const long long two_mjd =
        2 * ::cheatah::space::time::nasa_julian_day_at_noon(year, month, day) - 4800001;
    return e.delta_at_e7 * 100 + (two_mjd - 2 * e.drift_mjd) * e.drift_rate_e7 * 50;
}

/**
 * TAI-UTC in nanoseconds, computed the way NASA's CDF library computes it.
 *
 * A faithful reproduction of the reference implementation: accumulate in `double`, then
 * truncate toward zero. It exists so bit-identical agreement with NASA is available when that is
 * what is wanted — and so the difference between the two is a tested, enumerated fact rather
 * than a suspicion.
 *
 * @param year UTC year.
 * @param month UTC month, 1-12.
 * @param day UTC day of month.
 * @return TAI-UTC in nanoseconds as NASA would compute it; 0 before 1960-01-01.
 * @complexity O(leap_second_count()) — a linear scan of 42 entries.
 * @alloc none.
 * @systest systests/test_leapseconds.purr
 * @par Example
 * @code{.purr}
 * import io
 * import space.cdf as cdf
 *
 * # Identical to the exact form from 1972 onward; both are 37 s in 2026.
 * io.print(cdf.tai_minus_utc_ns_nasa_compat(2026, 8, 28))   # 37000000000
 * @endcode
 */
inline long long tai_minus_utc_ns_nasa_compat(int year, int month, int day) noexcept {
    const int i = detail::leap_index_for(year, month, day);
    if (i < 0) { return 0; }
    const LeapSecondEntry& e = detail::kLeapSeconds[i];
    const double mjd =
        static_cast<double>(::cheatah::space::time::nasa_julian_day_at_noon(year, month, day))
        - 2400000.5;
    const double delta_at = static_cast<double>(e.delta_at_e7) / 1.0e7;
    const double drift_rate = static_cast<double>(e.drift_rate_e7) / 1.0e7;
    const double seconds = delta_at + (mjd - static_cast<double>(e.drift_mjd)) * drift_rate;
    return static_cast<long long>(seconds * 1.0e9);  // truncates toward zero, as NASA's does
}

// ---- provenance --------------------------------------------------------------------------------
// The table is a snapshot of someone else's data, so it carries where it came from, when we last
// checked, and how far forward it can be trusted.

/**
 * Where the table comes from.
 * @return NASA's published leap-second table URL.
 * @complexity O(1).
 * @alloc none.
 * @systest systests/test_leapseconds.purr
 * @par Example
 * @code{.purr}
 * import io
 * import space.cdf as cdf
 *
 * io.print(cdf.leap_seconds_url())
 * @endcode
 */
constexpr std::string_view leap_seconds_url() noexcept {
    return "https://cdf.gsfc.nasa.gov/html/CDFLeapSeconds.txt";
}

/**
 * The `Updated:` stamp carried inside the upstream file itself.
 * @return the upstream stamp, `YYYYMMDD`.
 * @complexity O(1).
 * @alloc none.
 * @systest systests/test_leapseconds.purr
 * @par Example
 * @code{.purr}
 * import io
 * import space.cdf as cdf
 *
 * io.print(cdf.leap_seconds_upstream_updated())   # 20161025
 * @endcode
 */
constexpr std::string_view leap_seconds_upstream_updated() noexcept { return "20161025"; }

/**
 * SHA-256 of the NORMALIZED upstream table — comments and blank lines dropped, each remaining
 * line's whitespace collapsed and trimmed. Normalizing means a reformat or an edited comment
 * raises no false alarm, while any change to the data does.
 *
 * @return the 64-character lowercase hex digest.
 * @complexity O(1).
 * @alloc none.
 * @systest systests/test_leapseconds.purr
 * @par Example
 * @code{.purr}
 * import io
 * import space.cdf as cdf
 *
 * # scripts/check_leapseconds.sh recomputes this from the live file and compares.
 * io.print(cdf.leap_seconds_sha256())
 * @endcode
 */
constexpr std::string_view leap_seconds_sha256() noexcept {
    return "b21522ee7d9d763e281dce5c63e94b65be5de6107a59ef6e9ac355a58446264e";
}

/**
 * The date this copy of the table was last checked against upstream.
 * @return the verification date.
 * @complexity O(1).
 * @alloc none.
 * @systest systests/test_leapseconds.purr
 * @par Example
 * @code{.purr}
 * import io
 * import space.cdf as cdf
 *
 * let d = cdf.leap_seconds_verified_on()
 * io.print(d.year, d.month, d.day)   # 2026 8 28
 * @endcode
 */
constexpr ::cheatah::space::time::CivilDate leap_seconds_verified_on() noexcept {
    return ::cheatah::space::time::CivilDate{2026, 8, 28};
}

/**
 * The date through which this table is guaranteed COMPLETE.
 *
 * Derived, not guessed: IERS announces a leap second at least six months ahead in Bulletin C.
 * Upstream carried no entry after 2017-01-01 when we checked on leap_seconds_verified_on(), so
 * none can occur before six months after that date. Past this date the table may be merely
 * current rather than provably complete — ask leap_seconds_covers() before trusting a conversion
 * of a future timestamp.
 *
 * @return the last date the table is provably complete through.
 * @complexity O(1).
 * @alloc none.
 * @systest systests/test_leapseconds.purr
 * @par Example
 * @code{.purr}
 * import io
 * import space.cdf as cdf
 *
 * let d = cdf.leap_seconds_known_good_through()
 * io.print(d.year, d.month, d.day)   # 2027 2 28
 * @endcode
 */
constexpr ::cheatah::space::time::CivilDate leap_seconds_known_good_through() noexcept {
    return ::cheatah::space::time::CivilDate{2027, 2, 28};
}

/**
 * Whether a date falls in the range this table is provably complete for.
 *
 * A conversion outside it is not necessarily wrong — it is merely unverifiable, because a leap
 * second announced after leap_seconds_verified_on() would not be in this copy. Callers
 * converting future timestamps should branch on this rather than take a silently stale answer.
 *
 * @param year UTC year.
 * @param month UTC month, 1-12.
 * @param day UTC day of month.
 * @return true when the date is on or before leap_seconds_known_good_through().
 * @complexity O(1).
 * @alloc none.
 * @systest systests/test_leapseconds.purr
 * @par Example
 * @code{.purr}
 * import io
 * import space.cdf as cdf
 *
 * io.print(cdf.leap_seconds_covers(2026, 8, 28))   # True
 * io.print(cdf.leap_seconds_covers(2099, 1, 1))    # False — beyond the guarantee
 * @endcode
 */
constexpr bool leap_seconds_covers(int year, int month, int day) noexcept {
    const ::cheatah::space::time::CivilDate limit = leap_seconds_known_good_through();
    return ::cheatah::space::time::nasa_julian_day_at_noon(year, month, day)
           <= ::cheatah::space::time::nasa_julian_day_at_noon(limit.year, limit.month, limit.day);
}

}  // namespace cheatah::space::cdf
