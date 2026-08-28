// Unit tests for space/time/civil.hpp — the proleptic Gregorian calendar in exact integers.
//
// Every function is called at RUN TIME as well as in a static_assert. That distinction is not
// pedantry: clang's source-based coverage counts a constexpr function evaluated only in a
// static_assert as never executed, so a header proved correct entirely at compile time reports
// 0% and fails the QA gate's 100%-functions bar. The static_asserts are here because a wrong
// calendar should be a BUILD failure; the runtime calls are here so the coverage is honest.
//
// The centrepiece is CivilRange below. Two independent day-number algorithms live in civil.hpp —
// ours and a literal transcription of NASA's — and this sweeps every day in the CDF_TIME_TT2000
// range asserting they agree and that the round trip is exact. 213,505 days is cheap to check
// and it is what lets space.cdf claim agreement with NASA's epoch handling rather than assume it.
#include <gtest/gtest.h>

#include "space/space.hpp"

namespace st = cheatah::space::time;

// ---- compile-time surface --------------------------------------------------------------------

static_assert(st::is_leap_year(2000), "divisible by 400 is a leap year");
static_assert(!st::is_leap_year(1900), "divisible by 100 but not 400 is not");
static_assert(st::is_leap_year(2024));
static_assert(!st::is_leap_year(2026));
static_assert(st::days_in_month(2024, 2) == 29);
static_assert(st::days_in_month(2026, 2) == 28);
static_assert(st::days_from_civil(1970, 1, 1) == 0, "the Unix epoch is day zero");
static_assert(st::days_from_civil(2000, 1, 1) == 10957, "the J2000 day");
static_assert(st::nasa_julian_day_at_noon(2000, 1, 1) == 2451545, "J2000 Julian Day");
static_assert(st::civil_from_days(0).year == 1970);
static_assert(st::is_valid_civil(st::civil_date(2024, 2, 29)));
static_assert(!st::is_valid_civil(st::civil_date(2026, 2, 29)));

// ---- leap years --------------------------------------------------------------------------------

TEST(Civil, LeapYearRules) {
    EXPECT_TRUE(st::is_leap_year(2000));   // /400
    EXPECT_FALSE(st::is_leap_year(1900));  // /100 not /400
    EXPECT_FALSE(st::is_leap_year(2100));
    EXPECT_TRUE(st::is_leap_year(2024));   // /4
    EXPECT_FALSE(st::is_leap_year(2026));
    EXPECT_TRUE(st::is_leap_year(0));      // year 0 is divisible by 400
    EXPECT_TRUE(st::is_leap_year(-4));     // the rule is projected backwards unchanged
}

// ---- month lengths -----------------------------------------------------------------------------

TEST(Civil, DaysInMonthCoversEveryMonthAndBothFebruaries) {
    constexpr int common[13] = {0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
    for (int m = 1; m <= 12; ++m) {
        EXPECT_EQ(st::days_in_month(2026, m), common[m]) << "month " << m;
    }
    EXPECT_EQ(st::days_in_month(2024, 2), 29);
}

TEST(Civil, DaysInMonthRejectsOutOfRange) {
    // Both guard branches, so neither is an uncovered line.
    EXPECT_EQ(st::days_in_month(2026, 0), 0);
    EXPECT_EQ(st::days_in_month(2026, 13), 0);
    EXPECT_EQ(st::days_in_month(2026, -1), 0);
}

// ---- validity ----------------------------------------------------------------------------------

TEST(Civil, IsValidCivilRejectsEachWayIndependently) {
    EXPECT_TRUE(st::is_valid_civil(st::civil_date(2026, 8, 28)));
    EXPECT_TRUE(st::is_valid_civil(st::civil_date(2024, 2, 29)));
    EXPECT_FALSE(st::is_valid_civil(st::civil_date(2026, 2, 29)));  // day past month length
    EXPECT_FALSE(st::is_valid_civil(st::civil_date(2026, 0, 10)));  // month too small
    EXPECT_FALSE(st::is_valid_civil(st::civil_date(2026, 13, 10))); // month too large
    EXPECT_FALSE(st::is_valid_civil(st::civil_date(2026, 8, 0)));   // day too small
    EXPECT_FALSE(st::is_valid_civil(st::civil_date(2026, 8, 32)));  // day too large
}

TEST(Civil, CivilDateCarriesItsComponents) {
    const st::CivilDate d = st::civil_date(2026, 8, 28);
    EXPECT_EQ(d.year, 2026);
    EXPECT_EQ(d.month, 8);
    EXPECT_EQ(d.day, 28);
    EXPECT_EQ(st::CivilDate{}.year, 0);  // the aggregate's defaulted members
}

// ---- anchors -----------------------------------------------------------------------------------

TEST(Civil, KnownDayNumbers) {
    EXPECT_EQ(st::days_from_civil(1970, 1, 1), 0);
    EXPECT_EQ(st::days_from_civil(1970, 1, 2), 1);
    EXPECT_EQ(st::days_from_civil(1969, 12, 31), -1);
    EXPECT_EQ(st::days_from_civil(2000, 1, 1), 10957);   // J2000
    EXPECT_EQ(st::days_from_civil(2024, 2, 29), 19782);  // a real leap day
}

TEST(Civil, KnownCivilDates) {
    const st::CivilDate epoch = st::civil_from_days(0);
    EXPECT_EQ(epoch.year, 1970);
    EXPECT_EQ(epoch.month, 1);
    EXPECT_EQ(epoch.day, 1);

    const st::CivilDate j2000 = st::civil_from_days(10957);
    EXPECT_EQ(j2000.year, 2000);
    EXPECT_EQ(j2000.month, 1);
    EXPECT_EQ(j2000.day, 1);

    const st::CivilDate before = st::civil_from_days(-1);
    EXPECT_EQ(before.year, 1969);
    EXPECT_EQ(before.month, 12);
    EXPECT_EQ(before.day, 31);
}

TEST(Civil, NasaJulianDayAnchors) {
    EXPECT_EQ(st::nasa_julian_day_at_noon(2000, 1, 1), 2451545);  // J2000
    EXPECT_EQ(st::nasa_julian_day_at_noon(1970, 1, 1), 2440588);  // the Unix epoch
}

// ---- the exhaustive sweep ----------------------------------------------------------------------

// The TT2000 representable range: 1707-09-22 .. 2292-04-11, the days an int64 of nanoseconds
// since J2000 can name. Sweeping it is the whole point — see the file header.
TEST(Civil, RoundTripAndNasaAgreementAcrossTheWholeTT2000Range) {
    const long long lo = st::days_from_civil(1707, 9, 22);
    const long long hi = st::days_from_civil(2292, 4, 11);
    ASSERT_LT(lo, hi);
    ASSERT_EQ(hi - lo + 1, 213505) << "the TT2000 day range changed — check the bounds";

    long long round_trip_failures = 0;
    long long nasa_disagreements = 0;
    long long first_bad_day = 0;
    for (long long z = lo; z <= hi; ++z) {
        const st::CivilDate c = st::civil_from_days(z);
        if (st::days_from_civil(c.year, c.month, c.day) != z) {
            if (round_trip_failures == 0) { first_bad_day = z; }
            ++round_trip_failures;
        }
        // NASA counts Julian Days; we count days from 1970. 2440588 is the offset between them.
        if (st::nasa_julian_day_at_noon(c.year, c.month, c.day) != z + 2440588) {
            if (nasa_disagreements == 0) { first_bad_day = z; }
            ++nasa_disagreements;
        }
    }
    EXPECT_EQ(round_trip_failures, 0) << "first at day " << first_bad_day;
    EXPECT_EQ(nasa_disagreements, 0) << "first at day " << first_bad_day;
}

// A green sweep proves nothing unless it could have gone red, so provoke each failure mode on a
// deliberately wrong input and assert the checks fire.
TEST(Civil, TheSweepsChecksCanActuallyFail) {
    // A date that does not exist round-trips to a DIFFERENT day, which is what the sweep catches.
    const long long bogus = st::days_from_civil(2026, 2, 30);  // February has 28 days in 2026
    const st::CivilDate normalized = st::civil_from_days(bogus);
    EXPECT_FALSE(st::is_valid_civil(st::civil_date(2026, 2, 30)));
    EXPECT_EQ(normalized.month, 3);  // it lands in March
    EXPECT_EQ(normalized.day, 2);

    // And the NASA cross-check is a real equality, not a tautology: a wrong offset must fail it.
    EXPECT_NE(st::nasa_julian_day_at_noon(2000, 1, 1), st::days_from_civil(2000, 1, 1) + 2440587);
    EXPECT_EQ(st::nasa_julian_day_at_noon(2000, 1, 1), st::days_from_civil(2000, 1, 1) + 2440588);
}

// ---- negative and boundary years -----------------------------------------------------------------

// Opaque identity: keeps a value out of the constant folder.
//
// This matters more than it looks. Everything in civil.hpp is `constexpr`, so a call with
// literal arguments is evaluated at COMPILE time and emits no runtime code at all — the branch
// is then reported as never taken, and a test that "exercises" it proves nothing about the
// generated code. Routing the inputs through a volatile read forces a real call. The negative-era
// arms of days_from_civil() and civil_from_days() are only reachable this way.
namespace {
template <class T>
T opaque(T v) {
    volatile T sink = v;
    return sink;
}
}  // namespace

TEST(Civil, ProlepticYearsBeforeTheCommonEra) {
    // Year 0 exists in astronomical numbering and is a leap year.
    const long long y0 = st::days_from_civil(opaque(0), opaque(1), opaque(1));
    EXPECT_LT(y0, 0);
    const st::CivilDate back = st::civil_from_days(opaque(y0));
    EXPECT_EQ(back.year, 0);
    EXPECT_EQ(back.month, 1);
    EXPECT_EQ(back.day, 1);

    // CDF_EPOCH counts milliseconds from year 0, so this is the origin space.cdf will use.
    EXPECT_EQ(y0, -719528);
}

TEST(Civil, NegativeEraArmsAreReachedAtRuntime) {
    // days_from_civil takes `y - 399` only for a negative year; civil_from_days takes
    // `z - 146096` only for a day count before 0000-03-01. Both need genuinely negative,
    // non-constant-folded inputs — see opaque() above.
    for (int year : {-1, -4, -401, -1200}) {
        for (int month : {1, 3, 12}) {
            const long long z = st::days_from_civil(opaque(year), opaque(month), opaque(1));
            EXPECT_LT(z, -719528) << year << "-" << month;
            const st::CivilDate c = st::civil_from_days(opaque(z));
            EXPECT_EQ(c.year, year) << year << "-" << month;
            EXPECT_EQ(c.month, month) << year << "-" << month;
            EXPECT_EQ(c.day, 1) << year << "-" << month;
        }
    }
}

// NASA's formula needs floor division; C truncates toward zero. The two only differ once the
// numerator is negative, which pins the damage to negative years. This asserts BOTH halves of
// that boundary, because only the pair is useful: agreement above it is what lets space.cdf
// trust NASA's epoch handling, and divergence below it is what stops someone "correcting" the
// transcription into something that no longer reproduces the reference implementation.
TEST(Civil, NasaFormulaIsExactFromYearZeroAndWrongBeforeIt) {
    // Exact for every day of every year from 0 through 3000 — which covers all of CDF_EPOCH
    // (milliseconds from year 0) and all of TT2000 (1707..2292).
    long long disagreements = 0;
    for (int y = 0; y <= 3000; ++y) {
        for (int m = 1; m <= 12; ++m) {
            for (int d = 1, dim = st::days_in_month(y, m); d <= dim; ++d) {
                if (st::nasa_julian_day_at_noon(y, m, d) != st::days_from_civil(y, m, d) + 2440588) {
                    ++disagreements;
                }
            }
        }
    }
    EXPECT_EQ(disagreements, 0) << "NASA's formula must be exact wherever CDF can reach";

    // And genuinely wrong below year 0. Pinned to exact offsets so a silent change is caught.
    struct Known { int year; long long delta; };
    constexpr Known known[] = {{-1, -1}, {-2, -1}, {-100, -1}, {-401, -2}, {-1200, -1}};
    for (const Known& k : known) {
        const long long got = st::nasa_julian_day_at_noon(opaque(k.year), opaque(1), opaque(1));
        const long long exact = st::days_from_civil(opaque(k.year), opaque(1), opaque(1)) + 2440588;
        EXPECT_EQ(got - exact, k.delta) << "year " << k.year;
    }
    // Year 0 itself is on the correct side of the boundary.
    EXPECT_EQ(st::nasa_julian_day_at_noon(opaque(0), opaque(1), opaque(1)),
              st::days_from_civil(opaque(0), opaque(1), opaque(1)) + 2440588);
}
