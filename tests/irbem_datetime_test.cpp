// Unit tests for space.irbem's date/time layer — the Julian Day Number primitive, its inverse,
// day of year, decimal year, and the library-info shims.
//
// Two things drive the shape of this file.
//
// First, calendar arithmetic is exactly the kind of code that passes a handful of spot checks and
// still has a hole in it, so the load-bearing tests are SWEEPS over hundreds of thousands of dates
// asserting an invertibility property, not a table of anchors. The anchors are here too — as
// `static_assert`s, since every routine is `constexpr` — but they are the cheap half.
//
// Second, the QA gate requires 100% line AND function coverage of the header, and a `static_assert`
// contributes exactly nothing to a runtime coverage profile. So every entity is exercised twice:
// once at compile time (proof) and once at run time (coverage). The awkward branches — floor
// division of a negative, a decimal year that rounds up onto the next January 1 — each have a named
// test whose comment explains why that particular input reaches them.
//
// Arithmetic is asserted with `==` wherever the value is exactly representable: a Julian Day Number
// is an integer, 183 of 366 days is one half, and a whole number of seconds survives the
// microsecond quantization intact. The one place a tolerance is unavoidable is a round trip THROUGH
// a decimal year, where the year's magnitude eats the low bits of the fraction; that tolerance is
// derived from the double's ulp in a comment rather than picked to make the test pass.
#include <gtest/gtest.h>

#include <array>
#include <cstddef>
#include <cmath>
#include <cstdint>
#include <string_view>
#include <type_traits>

#include "alloc_counter.hpp"
#include "space/time/calendar.hpp"
#include "space/irbem/library_info.hpp"

namespace ib = cheatah::space::irbem;

using ib::CalendarDate;
using ib::DateTime;

// ---- compile-time anchors ---------------------------------------------------------------------

// The J2000 epoch. This is the single most-quoted Julian Day Number in the field, and it is also
// `space::time::jd_j2000()` — the two modules meet here, on a value neither one may move.
static_assert(ib::julian_day_number(2000, 1, 1) == 2451545);

// The origin of the Julian period: 1 January 4713 BC = astronomical year -4712, Julian calendar.
static_assert(ib::julian_day_number(-4712, 1, 1) == 0);

// The Modified Julian Date origin, 1858 November 17, is JD 2400000.5 — i.e. JDN 2400001 at noon.
static_assert(ib::julian_day_number(1858, 11, 17) == 2400001);

// The reform. 1582 October 4 (Julian) is followed directly by October 15 (Gregorian); the ten days
// between are not in the count, and the two branches join with no discontinuity in the NUMBER.
static_assert(ib::julian_day_number(1582, 10, 4) == 2299160);
static_assert(ib::julian_day_number(1582, 10, 15) == ib::gregorian_reform_jdn);
static_assert(ib::gregorian_reform_jdn == 2299161);
// A "date" inside the excised ten days is read as Julian, which is what the shipped IRBEM library
// does — measured through its C entry point, not read from its source.
static_assert(ib::julian_day_number(1582, 10, 5) == 2299161);
static_assert(ib::julian_day_number(1582, 10, 14) == 2299170);

// A year before the Julian-period origin, which is where floor division stops agreeing with C++'s
// truncating `/`. -5000-01-01 is 288 Julian years (72 of them leap) before JDN 0.
static_assert(ib::julian_day_number(-5000, 1, 1) == -(288 * 365 + 72));

// The inverse, on both sides of the reform and below zero.
static_assert(ib::calendar_date(2451545) == CalendarDate{2000, 1, 1});
static_assert(ib::calendar_date(0) == CalendarDate{-4712, 1, 1});
static_assert(ib::calendar_date(2299160) == CalendarDate{1582, 10, 4});
static_assert(ib::calendar_date(2299161) == CalendarDate{1582, 10, 15});
static_assert(ib::calendar_date(-105192) == CalendarDate{-5000, 1, 1});

// Day of year across the three leap rules that matter: divisible by 4, by 100, and by 400.
static_assert(ib::day_of_year(2000, 3, 1) == 61);  // 2000 is a leap year (divisible by 400)
static_assert(ib::day_of_year(1900, 3, 1) == 60);  // 1900 is not (century, not divisible by 400)
static_assert(ib::day_of_year(2100, 3, 1) == 60);  // nor is 2100
static_assert(ib::day_of_year(2004, 3, 1) == 61);
static_assert(ib::day_of_year(2000, 12, 31) == 366);
static_assert(ib::day_of_year(1582, 12, 31) == 355);  // the reform year lost ten days

static_assert(ib::days_in_year(1999) == 365);
static_assert(ib::days_in_year(2000) == 366);
static_assert(ib::days_in_year(1900) == 365);
static_assert(ib::days_in_year(1582) == 355);
static_assert(ib::days_in_year(1500) == 366);  // Julian: no century exception yet

static_assert(ib::is_leap_year(2000) && !ib::is_leap_year(1900) && !ib::is_leap_year(2100));
static_assert(ib::is_leap_year(2004) && !ib::is_leap_year(2003));
static_assert(ib::is_leap_year(-4) && !ib::is_leap_year(-2));  // negative years too

static_assert(ib::date_from_day_of_year(2000, 61) == CalendarDate{2000, 3, 1});
static_assert(ib::date_from_day_of_year(1900, 60) == CalendarDate{1900, 3, 1});

// Exactly representable decimal years: day 184 of a 366-day year is the halfway point, and 183 of
// 365 with twelve hours on it is too. Both are `.5` with no bits to spare, so `==` is honest here.
static_assert(ib::decimal_year(2000, 7, 2, 0, 0, 0) == 2000.5);
static_assert(ib::decimal_year(1999, 7, 2, 12, 0, 0) == 1999.5);
static_assert(ib::decimal_year(2000, 1, 1, 0, 0, 0) == 2000.0);

// ...and back again, on the values where the fraction is exact.
static_assert(ib::date_and_time_from_decimal_year(2000.5) ==
              DateTime{2000, 7, 2, 184, 0, 0, 0, 0.0});
static_assert(ib::date_and_time_from_decimal_year(1999.5) ==
              DateTime{1999, 7, 2, 183, 12, 0, 0, 43200.0});

static_assert(ib::date_and_time_from_doy_and_ut(2000, 184, 0.0) ==
              DateTime{2000, 7, 2, 184, 0, 0, 0, 0.0});
static_assert(ib::date_and_time_from_doy_and_ut(2004, 60, 21600.0) ==
              DateTime{2004, 2, 29, 60, 6, 0, 0, 21600.0});

// The library-info surface is a set of constants; pin them where they cannot drift unnoticed.
static_assert(ib::max_batch_times() == 100000);
static_assert(ib::igrf_generation() == 14);
static_assert(ib::implementation_version() == 1);
static_assert(!ib::implementation_release().empty());
// IRBEM's C entry point hands back an 80-character buffer; the tag has to fit in one.
static_assert(ib::implementation_release().size() <= 80);

// The types stay cheap and copyable — they travel by value through every model entry point.
static_assert(std::is_trivially_copyable_v<CalendarDate>);
static_assert(std::is_trivially_copyable_v<DateTime>);
static_assert(std::is_aggregate_v<CalendarDate>);

// ---- the sweeps -------------------------------------------------------------------------------

// The strongest statement available about julian_day_number/calendar_date: over a contiguous run of
// Julian Day Numbers spanning both calendars, calendar_date is injective and julian_day_number
// inverts it. A contiguous JDN range that round-trips is exactly the claim that consecutive civil
// days differ by one day — including at the reform, which is therefore not a special case here
// either.
TEST(IrbemDatetime, JulianDayNumberAndCalendarDateAreMutualInversesFrom1500To2200) {
    const std::int64_t first = ib::julian_day_number(1500, 1, 1);
    const std::int64_t last = ib::julian_day_number(2200, 12, 31);
    ASSERT_EQ(2268933, first);  // fixes the sweep's own endpoints, so a bug cannot shrink it
    ASSERT_EQ(2524958, last);

    // Inverting is not enough on its own: the month/day encoding the algorithm works in is
    // periodic, so a shifted-by-one decoding (April 0 for March 31, say) still round-trips
    // perfectly through julian_day_number. The range check is what pins the answer to a real
    // calendar date, and a deliberately perturbed build proved it is load-bearing.
    int failures = 0;
    int out_of_range = 0;
    std::int64_t first_failure = 0;
    for (std::int64_t jdn = first; jdn <= last; ++jdn) {
        const CalendarDate date = ib::calendar_date(jdn);
        if (date.month < 1 || date.month > 12 || date.day < 1 || date.day > 31) {
            if (out_of_range == 0) {
                ADD_FAILURE() << "JDN " << jdn << " decoded to " << date.year << '-' << date.month
                              << '-' << date.day;
            }
            ++out_of_range;
        }
        if (ib::julian_day_number(date.year, date.month, date.day) != jdn) {
            if (failures == 0) first_failure = jdn;
            ++failures;
        }
    }
    EXPECT_EQ(0, failures) << "first failing JDN " << first_failure << " of "
                           << (last - first + 1) << " swept";
    EXPECT_EQ(0, out_of_range);
}

// The same property below the Julian-period origin, where every intermediate in both directions is
// negative and floor division is the only thing keeping the answer right.
TEST(IrbemDatetime, TheInversionHoldsForNegativeJulianDayNumbers) {
    int failures = 0;
    for (std::int64_t jdn = -120000; jdn <= 0; ++jdn) {
        const CalendarDate date = ib::calendar_date(jdn);
        if (date.month < 1 || date.month > 12 || date.day < 1 || date.day > 31) ++failures;
        if (ib::julian_day_number(date.year, date.month, date.day) != jdn) ++failures;
    }
    EXPECT_EQ(0, failures);

    // Truncating division would put this one day late; floor division puts it here.
    EXPECT_EQ(-105192, ib::julian_day_number(-5000, 1, 1));
    EXPECT_EQ((CalendarDate{-5000, 1, 1}), ib::calendar_date(-105192));
}

// day_of_year and date_from_day_of_year invert each other, and the day of year never leaves
// [1, days_in_year]. Sweeping every year from 1500 covers the century non-leap years 1700/1800/1900
// and 2100, the divisible-by-400 leap years 1600 and 2000, the Julian-rule years before 1582 (where
// 1500 IS a leap year), and the 355-day reform year itself.
TEST(IrbemDatetime, DayOfYearInvertsAcrossEveryYearFrom1500To2200) {
    int failures = 0;
    long long checked = 0;
    for (int year = 1500; year <= 2200; ++year) {
        const int length = ib::days_in_year(year);
        for (int doy = 1; doy <= length; ++doy) {
            const CalendarDate date = ib::date_from_day_of_year(year, doy);
            if (date.year != year || ib::day_of_year(year, date.month, date.day) != doy ||
                date.month < 1 || date.month > 12 || date.day < 1 || date.day > 31) {
                ++failures;
            }
            ++checked;
        }
    }
    EXPECT_EQ(0, failures);
    EXPECT_EQ(256026, checked);  // the sweep really did run: 701 years of days
}

// days_in_year is measured (a difference of two Julian Day Numbers) while is_leap_year is the
// Gregorian RULE. From 1583 on they must agree exactly — two independent derivations of the same
// fact, which is the point of keeping both.
TEST(IrbemDatetime, TheLeapRuleAgreesWithTheMeasuredYearLengthAfterTheReform) {
    for (int year = 1583; year <= 2400; ++year) {
        EXPECT_EQ(365 + (ib::is_leap_year(year) ? 1 : 0), ib::days_in_year(year)) << year;
    }
    // ...and before it they must NOT, which is the whole reason days_in_year is the authority.
    EXPECT_FALSE(ib::is_leap_year(1500));
    EXPECT_EQ(366, ib::days_in_year(1500));  // Julian: divisible by 4, no century exception
    EXPECT_EQ(355, ib::days_in_year(1582));  // ten days struck out of it
}

// ---- decimal year -----------------------------------------------------------------------------

namespace {

/// A wall-clock time of day to sweep. Chosen for the boundaries that break a floor-based inverse:
/// midnight, an exact hour, an exact minute, one second either side of a rollover.
struct TimeOfDay {
    int hour;
    int minute;
    int second;
};

constexpr std::array<TimeOfDay, 8> kTimes{TimeOfDay{0, 0, 0},   TimeOfDay{0, 0, 1},
                                          TimeOfDay{1, 0, 0},   TimeOfDay{6, 30, 30},
                                          TimeOfDay{12, 0, 0},  TimeOfDay{12, 34, 56},
                                          TimeOfDay{23, 0, 0},  TimeOfDay{23, 59, 59}};

}  // namespace

// The decimal year round trip, over thousands of instants. Every INTEGER field must come back
// exactly — that is the assertion IRBEM's own routine fails, by a whole second, on roughly one in
// six of a comparable sweep (measured through its C entry point).
TEST(IrbemDatetime, DecimalYearRoundTripsExactlyInEveryIntegerField) {
    // Every 7th day of every year from 1600 to 2200, at each of the eight boundary times below.
    constexpr long long kExpectedInstants = 254824;

    // Every field, seconds-of-day included, comes back with `==` and no tolerance. That is only
    // true because the inverse snaps to the millisecond: a decimal year near 2200 carries about
    // 15 us of noise in its fraction (ulp(2200) times a year of seconds), which is far too coarse
    // for a microsecond count and far too fine to disturb a millisecond one.
    long long checked = 0;
    for (int year = 1600; year <= 2200; ++year) {
        const int length = ib::days_in_year(year);
        for (int doy = 1; doy <= length; doy += 7) {
            const CalendarDate date = ib::date_from_day_of_year(year, doy);
            for (const TimeOfDay& t : kTimes) {
                const double decy =
                    ib::decimal_year(year, date.month, date.day, t.hour, t.minute, t.second);
                ASSERT_GE(decy, static_cast<double>(year));
                ASSERT_LT(decy, static_cast<double>(year + 1));

                const DateTime got = ib::date_and_time_from_decimal_year(decy);
                ASSERT_EQ(year, got.year) << decy;
                ASSERT_EQ(date.month, got.month) << decy;
                ASSERT_EQ(date.day, got.day) << decy;
                ASSERT_GE(got.month, 1);
                ASSERT_LE(got.month, 12);
                ASSERT_GE(got.day, 1);
                ASSERT_LE(got.day, 31);
                ASSERT_EQ(doy, got.day_of_year) << decy;
                ASSERT_EQ(t.hour, got.hour) << decy;
                ASSERT_EQ(t.minute, got.minute) << decy;
                ASSERT_EQ(t.second, got.second) << decy;

                const double expected_ut = (t.hour * 3600.0) + (t.minute * 60.0) + t.second;
                ASSERT_EQ(expected_ut, got.ut_seconds) << decy;
                ++checked;
            }
        }
    }
    EXPECT_EQ(kExpectedInstants, checked);
}

// The two producers must agree: going through a decimal year and going through (doy, UT) are two
// routes to the same instant, and they share only the assembly step.
TEST(IrbemDatetime, TheDoyAndUtRouteAgreesWithTheDecimalYearRoute) {
    for (int year = 1990; year <= 2010; ++year) {
        const int length = ib::days_in_year(year);
        for (int doy = 1; doy <= length; doy += 17) {
            for (const TimeOfDay& t : kTimes) {
                const double ut = (t.hour * 3600.0) + (t.minute * 60.0) + t.second;
                const DateTime direct = ib::date_and_time_from_doy_and_ut(year, doy, ut);

                // The direct route has no floating-point year in it at all, so it is exact.
                EXPECT_EQ(year, direct.year);
                EXPECT_EQ(doy, direct.day_of_year);
                EXPECT_EQ(t.hour, direct.hour);
                EXPECT_EQ(t.minute, direct.minute);
                EXPECT_EQ(t.second, direct.second);
                EXPECT_EQ(ut, direct.ut_seconds);

                const CalendarDate date = ib::date_from_day_of_year(year, doy);
                const DateTime viaDecy = ib::date_and_time_from_decimal_year(
                    ib::decimal_year(year, date.month, date.day, t.hour, t.minute, t.second));
                EXPECT_EQ(direct.year, viaDecy.year);
                EXPECT_EQ(direct.month, viaDecy.month);
                EXPECT_EQ(direct.day, viaDecy.day);
                EXPECT_EQ(direct.day_of_year, viaDecy.day_of_year);
                EXPECT_EQ(direct.hour, viaDecy.hour);
                EXPECT_EQ(direct.minute, viaDecy.minute);
                EXPECT_EQ(direct.second, viaDecy.second);
            }
        }
    }
}

// The exact values, spelled out. A decimal year is a rational with a power-of-two-friendly
// denominator surprisingly often, and where it is, nothing is approximate.
TEST(IrbemDatetime, DecimalYearHasExactValuesWhereTheFractionIsExact) {
    EXPECT_EQ(2000.0, ib::decimal_year(2000, 1, 1, 0, 0, 0));
    EXPECT_EQ(2000.5, ib::decimal_year(2000, 7, 2, 0, 0, 0));    // 183 of 366 days
    EXPECT_EQ(1999.5, ib::decimal_year(1999, 7, 2, 12, 0, 0));   // 182.5 of 365 days
    EXPECT_EQ(2000.25, ib::decimal_year(2000, 4, 1, 12, 0, 0));   // 91.5 of 366 days
    EXPECT_EQ(2000.75, ib::decimal_year(2000, 10, 1, 12, 0, 0));  // 274.5 of 366 days

    // The denominator is the year's OWN length, so the same calendar instant is a different
    // fraction in a leap year than in a common one. Getting this wrong is how an IGRF secular
    // interpolation picks up a systematic quarter-day error.
    EXPECT_NE(ib::decimal_year(2000, 7, 1, 0, 0, 0), ib::decimal_year(2001, 7, 1, 0, 0, 0) - 1.0);
}

// The specific defect this implementation does not have. A whole hour recovered from a decimal year
// is a few nanoseconds short of the hour; a bare floor turns that into the previous second. The
// shipped IRBEM library returns 11:59:59 for this input — measured.
TEST(IrbemDatetime, AWholeHourSurvivesTheDecimalYearRoundTrip) {
    const DateTime noon =
        ib::date_and_time_from_decimal_year(ib::decimal_year(2001, 1, 1, 12, 0, 0));
    EXPECT_EQ(12, noon.hour);
    EXPECT_EQ(0, noon.minute);
    EXPECT_EQ(0, noon.second);
    EXPECT_EQ(1, noon.day_of_year);

    const DateTime last_second =
        ib::date_and_time_from_decimal_year(ib::decimal_year(2000, 12, 31, 23, 59, 59));
    EXPECT_EQ((CalendarDate{2000, 12, 31}),
              (CalendarDate{last_second.year, last_second.month, last_second.day}));
    EXPECT_EQ(366, last_second.day_of_year);
    EXPECT_EQ(23, last_second.hour);
    EXPECT_EQ(59, last_second.minute);
    EXPECT_EQ(59, last_second.second);
}

// ---- the awkward branches, each reached on purpose ---------------------------------------------

// A UT offset outside [0, 86400) is carried into the day count rather than rejected, in BOTH
// directions. The negative direction is also the only caller that reaches round_to_int's
// negative-argument branch and floor_div's negative-numerator branch on the microsecond path.
TEST(IrbemDatetime, AUtOffsetOutsideTheDayCarriesIntoTheDate) {
    // 1.5 s before the start of January 2 is January 1 at 23:59:58.5.
    const DateTime before = ib::date_and_time_from_doy_and_ut(2000, 2, -1.5);
    EXPECT_EQ(2000, before.year);
    EXPECT_EQ(1, before.month);
    EXPECT_EQ(1, before.day);
    EXPECT_EQ(1, before.day_of_year);
    EXPECT_EQ(23, before.hour);
    EXPECT_EQ(59, before.minute);
    EXPECT_EQ(58, before.second);
    EXPECT_EQ(86398.5, before.ut_seconds);

    // A full day past January 1 is January 2 — and the reported day of year follows the date.
    const DateTime after = ib::date_and_time_from_doy_and_ut(2000, 1, 86400.0);
    EXPECT_EQ((CalendarDate{2000, 1, 2}), (CalendarDate{after.year, after.month, after.day}));
    EXPECT_EQ(2, after.day_of_year);
    EXPECT_EQ(0, after.hour);
    EXPECT_EQ(0.0, after.ut_seconds);

    // Crossing a year boundary backwards lands in the previous year, with ITS day count.
    const DateTime new_year = ib::date_and_time_from_doy_and_ut(2000, 1, -1.0);
    EXPECT_EQ(1999, new_year.year);
    EXPECT_EQ(365, new_year.day_of_year);
    EXPECT_EQ(23, new_year.hour);
    EXPECT_EQ(59, new_year.minute);
    EXPECT_EQ(59, new_year.second);
}

// A negative decimal year — the only input that reaches floor_to_int's negative branch. Year -2 is
// 3 BC in civil numbering; -1.25 is three quarters of the way through it. The Julian calendar makes
// that year 365 days, so three quarters is day 274 at 18:00 exactly, and the value round-trips
// bit-for-bit because 3/4 is a binary fraction.
TEST(IrbemDatetime, ANegativeDecimalYearIsHandledInAstronomicalNumbering) {
    const DateTime t = ib::date_and_time_from_decimal_year(-1.25);
    EXPECT_EQ(-2, t.year);
    EXPECT_EQ(274, t.day_of_year);
    EXPECT_EQ(18, t.hour);
    EXPECT_EQ(0, t.minute);
    EXPECT_EQ(0, t.second);
    EXPECT_EQ(64800.0, t.ut_seconds);

    EXPECT_EQ(-1.25, ib::decimal_year(-2, t.month, t.day, t.hour, t.minute, t.second));
}

// The microsecond rounding can land on the far end of a year — the instant that IS the next
// January 1 — and the routine has to step the year rather than report day 366 of a 365-day one.
//
// The largest decimal year strictly below 2001 is 2001 minus one ulp, or 2.3e-13 short — which is
// 7 microseconds short of the new year, well inside the half-millisecond snap. So this is not a
// contrived epoch: it is simply the last representable instant of the year 2000.
TEST(IrbemDatetime, ADecimalYearThatRoundsOntoTheNextJanuaryFirstStepsTheYear) {
    const double just_under = std::nextafter(2001.0, 2000.0);
    ASSERT_LT(just_under, 2001.0);
    ASSERT_GT(just_under, 2000.0);

    const DateTime t = ib::date_and_time_from_decimal_year(just_under);
    EXPECT_EQ(2001, t.year);
    EXPECT_EQ(1, t.month);
    EXPECT_EQ(1, t.day);
    EXPECT_EQ(1, t.day_of_year);
    EXPECT_EQ(0, t.hour);
    EXPECT_EQ(0, t.minute);
    EXPECT_EQ(0, t.second);
    EXPECT_EQ(0.0, t.ut_seconds);
}

// The reform year, end to end. October 5 through 14 of 1582 are dates that never happened; the
// library reads them as Julian, which is what the shipped IRBEM does, so the two agree everywhere a
// caller can actually be.
TEST(IrbemDatetime, TheGregorianReformIsHandledInEveryRoutine) {
    EXPECT_EQ(2299160, ib::julian_day_number(1582, 10, 4));
    EXPECT_EQ(ib::gregorian_reform_jdn, ib::julian_day_number(1582, 10, 15));
    EXPECT_EQ(ib::gregorian_reform_jdn, ib::julian_day_number(1582, 10, 5));

    EXPECT_EQ((CalendarDate{1582, 10, 4}), ib::calendar_date(2299160));
    EXPECT_EQ((CalendarDate{1582, 10, 15}), ib::calendar_date(2299161));

    EXPECT_EQ(277, ib::day_of_year(1582, 10, 4));
    EXPECT_EQ(278, ib::day_of_year(1582, 10, 15));
    EXPECT_EQ(355, ib::day_of_year(1582, 12, 31));
    EXPECT_EQ(355, ib::days_in_year(1582));

    // The decimal year divides by the year's real length, so December 31 of 1582 is 354/355 of the
    // way through it — not 364/365.
    EXPECT_DOUBLE_EQ(1582.0 + (354.0 / 355.0), ib::decimal_year(1582, 12, 31, 0, 0, 0));

    const DateTime t = ib::date_and_time_from_doy_and_ut(1582, 278, 0.0);
    EXPECT_EQ((CalendarDate{1582, 10, 15}), (CalendarDate{t.year, t.month, t.day}));
}

// ---- the comparison operators and the library-info shims ---------------------------------------

TEST(IrbemDatetime, TheAggregatesCompareFieldByField) {
    const CalendarDate a{2000, 1, 1};
    EXPECT_TRUE(a == (CalendarDate{2000, 1, 1}));
    EXPECT_FALSE(a == (CalendarDate{2001, 1, 1}));
    EXPECT_FALSE(a == (CalendarDate{2000, 2, 1}));
    EXPECT_FALSE(a == (CalendarDate{2000, 1, 2}));
    EXPECT_TRUE(a != (CalendarDate{2000, 1, 2}));

    const DateTime t{2000, 7, 2, 184, 6, 30, 30, 23430.0};
    EXPECT_TRUE(t == (DateTime{2000, 7, 2, 184, 6, 30, 30, 23430.0}));
    EXPECT_FALSE(t == (DateTime{2001, 7, 2, 184, 6, 30, 30, 23430.0}));
    EXPECT_FALSE(t == (DateTime{2000, 8, 2, 184, 6, 30, 30, 23430.0}));
    EXPECT_FALSE(t == (DateTime{2000, 7, 3, 184, 6, 30, 30, 23430.0}));
    EXPECT_FALSE(t == (DateTime{2000, 7, 2, 185, 6, 30, 30, 23430.0}));
    EXPECT_FALSE(t == (DateTime{2000, 7, 2, 184, 7, 30, 30, 23430.0}));
    EXPECT_FALSE(t == (DateTime{2000, 7, 2, 184, 6, 31, 30, 23430.0}));
    EXPECT_FALSE(t == (DateTime{2000, 7, 2, 184, 6, 30, 31, 23430.0}));
    EXPECT_FALSE(t == (DateTime{2000, 7, 2, 184, 6, 30, 30, 23431.0}));
    EXPECT_TRUE(t != (DateTime{2000, 7, 2, 184, 6, 30, 30, 23431.0}));

    // Both aggregates default to a valid date rather than to all-zeroes, since month 0 and day 0
    // are not dates and a default-constructed value must not look like one.
    const CalendarDate defaulted{};
    EXPECT_EQ((CalendarDate{0, 1, 1}), defaulted);
    const DateTime defaulted_time{};
    EXPECT_EQ(1, defaulted_time.month);
    EXPECT_EQ(1, defaulted_time.day);
    EXPECT_EQ(1, defaulted_time.day_of_year);
    EXPECT_EQ(0.0, defaulted_time.ut_seconds);
}

// These four are IRBEM's library-info routines. Three report facts about this build; the fourth
// deliberately does NOT impersonate a Fortran revision.
TEST(IrbemDatetime, TheLibraryInfoRoutinesReportThisImplementation) {
    EXPECT_EQ(100000, ib::max_batch_times());
    EXPECT_EQ(14, ib::igrf_generation());
    EXPECT_EQ(1, ib::implementation_version());

    const std::string_view release = ib::implementation_release();
    EXPECT_FALSE(release.empty());
    EXPECT_LE(release.size(), 80U);  // IRBEM's C boundary hands back an 80-character buffer
    // The tag has to say what it is: a log line carrying it must not read as IRBEM's own release.
    EXPECT_NE(std::string_view::npos, release.find("cheatah-space"));
    EXPECT_NE(std::string_view::npos, release.find("not IRBEM Fortran"));
}

// ---- everything runs -------------------------------------------------------------------------

// A last pass that simply calls every entry point at run time. Two thirds of this header's surface
// is `constexpr` and pinned by `static_assert` above, but a compile-time evaluation contributes
// nothing to a coverage profile — this is the test that makes the profile honest.
//
// It also stands as the allocation check. Nothing in this header can touch the heap — it returns
// only scalars, two trivially copyable aggregates, and a `string_view` into a literal — and the
// shared tripwire in alloc_counter.hpp measures that rather than leaving it as a claim. The reading
// is taken AFTER a first pass, so a routine that lazily built a workspace on first use would still
// be caught on the second.
TEST(IrbemDatetime, EveryEntryPointRunsAtRunTime) {
    double sink = 0.0;
    for (int year = 1980; year <= 2020; ++year) {
        const DateTime t =
            ib::date_and_time_from_decimal_year(ib::decimal_year(year, 6, 15, 12, 34, 56));
        const DateTime u = ib::date_and_time_from_doy_and_ut(year, t.day_of_year, t.ut_seconds);
        EXPECT_EQ(t.hour, u.hour);
        EXPECT_EQ(t.minute, u.minute);

        const CalendarDate d = ib::calendar_date(ib::julian_day_number(year, 6, 15));
        EXPECT_EQ((CalendarDate{year, 6, 15}), d);

        sink += static_cast<double>(ib::julian_day_number(year, 1, 1)) +
                ib::day_of_year(year, 6, 15) + ib::days_in_year(year) +
                (ib::is_leap_year(year) ? 1 : 0) + ib::date_from_day_of_year(year, 1).month +
                ib::max_batch_times() + ib::igrf_generation() + ib::implementation_version() +
                static_cast<double>(ib::implementation_release().size()) +
                static_cast<double>(ib::gregorian_reform_jdn) + u.ut_seconds;
    }
    EXPECT_NE(0.0, sink);

    const std::size_t before = cheatah_space_test::allocation_count();
    double second_pass = 0.0;
    for (int year = 1980; year <= 2020; ++year) {
        const DateTime t =
            ib::date_and_time_from_decimal_year(ib::decimal_year(year, 6, 15, 12, 34, 56));
        second_pass += ib::date_and_time_from_doy_and_ut(year, t.day_of_year, t.ut_seconds)
                           .ut_seconds +
                       static_cast<double>(ib::julian_day_number(year, 1, 1)) +
                       ib::day_of_year(year, 6, 15) +
                       static_cast<double>(ib::calendar_date(2451545).year) +
                       static_cast<double>(ib::implementation_release().size());
    }
    EXPECT_EQ(before, cheatah_space_test::allocation_count());
    EXPECT_NE(0.0, second_pass);
}
