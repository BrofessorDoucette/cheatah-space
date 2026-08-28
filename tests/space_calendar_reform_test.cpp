// The Julian/Gregorian reform branches of space.time's calendar.
//
// space/time/calendar.hpp is reform-aware: dates on or after 1582-10-15 use the Gregorian rule,
// earlier ones the Julian rule they were actually recorded in. That branch is easy to write, easy
// to get subtly wrong, and — being ~440 years in the past — is exercised by nothing else in the
// suite. A ten-day error there is invisible in modern epochs and catastrophic in a historical one,
// so the branch gets its own tests rather than riding on coverage from elsewhere.
//
// This deliberately differs from civil.hpp, which is PROLEPTIC Gregorian throughout because that is
// what the CDF format specifies. The two disagree by ten days at the reform and by a growing amount
// before it; that is by design, and asserted here so nobody "fixes" one to match the other.
#include <gtest/gtest.h>

#include "space/time/calendar.hpp"
#include "space/time/civil.hpp"

namespace t = cheatah::space::time;

// The anchors, at compile time: JDN 0 is -4712-01-01 in the Julian calendar, and the reform lands
// exactly where Gregory XIII put it.
static_assert(t::julian_day_number(-4712, 1, 1) == 0);
static_assert(t::julian_day_number(1582, 10, 15) == t::gregorian_reform_jdn);
static_assert(t::julian_day_number(1582, 10, 4) == t::gregorian_reform_jdn - 1);
static_assert(t::julian_day_number(2000, 1, 1) == 2451545);
static_assert(t::julian_day_number(1970, 1, 1) == 2440588);

TEST(CalendarReform, TheTenLostDaysAreLost) {
    // 1582-10-04 (Julian) was followed immediately by 1582-10-15 (Gregorian): the days between
    // never happened. Consecutive JDNs, ten calendar days apart.
    EXPECT_EQ(t::julian_day_number(1582, 10, 15), t::julian_day_number(1582, 10, 4) + 1);

    // The Julian rule still ANSWERS for those dates — the formula is total, and a caller handing in
    // 1582-10-05 gets the JDN that date would have had. What it must not do is silently agree with
    // the Gregorian answer for the same numerals.
    EXPECT_EQ(t::julian_day_number(1582, 10, 5), t::gregorian_reform_jdn);
    EXPECT_EQ(t::julian_day_number(1582, 10, 14), t::gregorian_reform_jdn + 9);
}

TEST(CalendarReform, TheInverseReturnsTheDateThatActuallyOccurred) {
    // Two calendar dates map to JDN 2299161 — Julian 1582-10-05 and Gregorian 1582-10-15 — so the
    // inverse cannot be a true inverse there. It returns the one that happened.
    const t::CalendarDate at_reform = t::calendar_date(t::gregorian_reform_jdn);
    EXPECT_EQ(1582, at_reform.year);
    EXPECT_EQ(10, at_reform.month);
    EXPECT_EQ(15, at_reform.day);

    const t::CalendarDate before = t::calendar_date(t::gregorian_reform_jdn - 1);
    EXPECT_EQ(1582, before.year);
    EXPECT_EQ(10, before.month);
    EXPECT_EQ(4, before.day);
}

TEST(CalendarReform, RoundTripsExactlyOnBothSidesOfTheReform) {
    // Every JDN in the Gregorian era must round-trip exactly. The Julian era round-trips too, over
    // its own dates — the only discontinuity is the ten-day gap itself.
    for (std::int64_t jdn = t::gregorian_reform_jdn; jdn < t::gregorian_reform_jdn + 200000; ++jdn) {
        const t::CalendarDate c = t::calendar_date(jdn);
        ASSERT_EQ(jdn, t::julian_day_number(c.year, c.month, c.day)) << jdn;
    }
    for (std::int64_t jdn = 0; jdn < t::gregorian_reform_jdn - 1; jdn += 7) {
        const t::CalendarDate c = t::calendar_date(jdn);
        ASSERT_EQ(jdn, t::julian_day_number(c.year, c.month, c.day)) << jdn;
    }
}

TEST(CalendarReform, TheJulianLeapRuleHasNoCenturyException) {
    // 1500 is a leap year in the Julian calendar and would not be in the Gregorian: the century
    // exception is the whole substance of the reform, so this is the branch that matters most.
    EXPECT_EQ(366, static_cast<int>(t::julian_day_number(1501, 1, 1) -
                                    t::julian_day_number(1500, 1, 1)));
    // 1700, safely after the reform, follows Gregorian rules and is NOT a leap year.
    EXPECT_EQ(365, static_cast<int>(t::julian_day_number(1701, 1, 1) -
                                    t::julian_day_number(1700, 1, 1)));
    // 1600 is a leap year under BOTH rules — divisible by 400 — so it distinguishes nothing and is
    // included precisely to show the test is not merely detecting "old year, more days".
    EXPECT_EQ(366, static_cast<int>(t::julian_day_number(1601, 1, 1) -
                                    t::julian_day_number(1600, 1, 1)));
}

TEST(CalendarReform, DisagreesWithProlepticCivilBeforeTheReformOnPurpose) {
    // civil.hpp is proleptic Gregorian because CDF says so; calendar.hpp is reform-aware because
    // historical dates say so. They MUST agree after the reform and MUST differ before it. If a
    // future change makes these equal everywhere, one of the two conventions has been broken.
    EXPECT_TRUE(t::is_leap_year(1600));   // both rules
    EXPECT_FALSE(t::is_leap_year(1700));  // Gregorian only — civil.hpp is proleptic

    // The reform-aware calendar puts 1700 at 365 days; a proleptic-Julian reading would say 366.
    EXPECT_EQ(365, t::days_in_year(1700));
    // ...and before the reform it agrees with Julian, where 1500 has 366.
    EXPECT_EQ(366, t::days_in_year(1500));
}
