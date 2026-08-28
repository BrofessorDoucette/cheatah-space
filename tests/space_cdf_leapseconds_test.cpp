// Unit tests for space/cdf/leapseconds.hpp — NASA's leap-second table as exact integers.
//
// Two things are being established here, and only the pair is useful.
//
// First, that the table itself is intact: 42 entries, ascending, 14 of them in the pre-1972
// drift era, with the two DECREASING entries (1961-08-01, 1968-02-01) still present. Those two
// are not trivia — they are the reason the reverse TT2000 conversion cannot assume the offset
// only grows, and a "tidied" table that dropped them would break that code silently.
//
// Second, that our exact arithmetic and NASA's double-then-truncate arithmetic differ EXACTLY
// where we say they do: on 430 of the 4384 days before 1972, always by 1 ns, and never
// afterwards. That count is asserted rather than described. If a change makes us agree with
// NASA everywhere, we stopped being exact; if it makes us differ after 1972, we broke CDF
// compatibility. Both are build failures.
#include <gtest/gtest.h>

#include "space/space.hpp"

namespace cdf = cheatah::space::cdf;
namespace st = cheatah::space::time;

// ---- compile-time surface ----------------------------------------------------------------------

static_assert(cdf::leap_second_count() == 42);
static_assert(cdf::leap_second_drift_rows() == 14);
static_assert(cdf::leap_second_entry(41).year == 2017, "the most recent entry");
static_assert(cdf::leap_second_entry(41).delta_at_e7 == 370'000'000, "37.0 s");
static_assert(cdf::tai_minus_utc_ns(2026, 8, 28) == 37'000'000'000LL);

// ---- the table ---------------------------------------------------------------------------------

TEST(LeapSeconds, TableShape) {
    EXPECT_EQ(cdf::leap_second_count(), 42);
    EXPECT_EQ(cdf::leap_second_drift_rows(), 14);

    const cdf::LeapSecondEntry first = cdf::leap_second_entry(0);
    EXPECT_EQ(first.year, 1960);
    EXPECT_EQ(first.month, 1);
    EXPECT_EQ(first.day, 1);

    const cdf::LeapSecondEntry last = cdf::leap_second_entry(cdf::leap_second_count() - 1);
    EXPECT_EQ(last.year, 2017);
    EXPECT_EQ(last.month, 1);
    EXPECT_EQ(last.day, 1);
    EXPECT_EQ(last.delta_at_e7, 370'000'000);  // 37.0 s
    EXPECT_EQ(last.drift_mjd, 0);
    EXPECT_EQ(last.drift_rate_e7, 0);
}

TEST(LeapSeconds, EntriesAscendByDateAndDriftOnlyBefore1972) {
    // No sentinel value: day numbers are days from 1970, so the whole pre-1972 era is NEGATIVE
    // and any "impossible" initial value is a real day for some entry.
    long long previous_day = 0;
    bool have_previous = false;
    long long drift_rows = 0;
    for (long long i = 0; i < cdf::leap_second_count(); ++i) {
        const cdf::LeapSecondEntry e = cdf::leap_second_entry(i);
        const long long day = st::days_from_civil(e.year, e.month, e.day);
        if (have_previous) {
            EXPECT_GT(day, previous_day) << "entry " << i << " is not after its predecessor";
        }
        previous_day = day;
        have_previous = true;

        const bool drifts = (e.drift_rate_e7 != 0);
        if (drifts) { ++drift_rows; }
        // The drift columns are a pre-1972 phenomenon, and they travel together.
        EXPECT_EQ(drifts, e.year < 1972) << "entry " << i;
        EXPECT_EQ(e.drift_mjd != 0, drifts) << "entry " << i;
    }
    EXPECT_EQ(drift_rows, cdf::leap_second_drift_rows());
}

// The offset is NOT monotonic. Two real entries go down, and code downstream depends on that
// staying true, so it is asserted rather than left as a comment.
TEST(LeapSeconds, TheTwoDecreasingEntriesAreStillThere) {
    std::vector<std::string> decreases;
    for (long long i = 1; i < cdf::leap_second_count(); ++i) {
        const cdf::LeapSecondEntry prev = cdf::leap_second_entry(i - 1);
        const cdf::LeapSecondEntry cur = cdf::leap_second_entry(i);
        if (cur.delta_at_e7 < prev.delta_at_e7) {
            decreases.push_back(std::to_string(cur.year) + "-" + std::to_string(cur.month));
        }
    }
    ASSERT_EQ(decreases.size(), 2u) << "a negative step vanished from the table";
    EXPECT_EQ(decreases[0], "1961-8");
    EXPECT_EQ(decreases[1], "1968-2");
}

TEST(LeapSeconds, OutOfRangeIndexYieldsAZeroedEntry) {
    for (long long i : {-1LL, -100LL, cdf::leap_second_count(), cdf::leap_second_count() + 5}) {
        const cdf::LeapSecondEntry e = cdf::leap_second_entry(i);
        EXPECT_EQ(e.year, 0) << "index " << i;
        EXPECT_EQ(e.delta_at_e7, 0) << "index " << i;
    }
}

// ---- the lookup ----------------------------------------------------------------------------------

TEST(LeapSeconds, KnownOffsets) {
    // Whole-second era: the value is the table entry, exactly.
    EXPECT_EQ(cdf::tai_minus_utc_ns(2026, 8, 28), 37'000'000'000LL);
    EXPECT_EQ(cdf::tai_minus_utc_ns(2017, 1, 1), 37'000'000'000LL);
    EXPECT_EQ(cdf::tai_minus_utc_ns(2016, 12, 31), 36'000'000'000LL);
    EXPECT_EQ(cdf::tai_minus_utc_ns(1972, 1, 1), 10'000'000'000LL);

    // Before the table begins there is nothing to report.
    EXPECT_EQ(cdf::tai_minus_utc_ns(1959, 12, 31), 0);
    EXPECT_EQ(cdf::tai_minus_utc_ns(1900, 1, 1), 0);
    EXPECT_EQ(cdf::tai_minus_utc_ns_nasa_compat(1959, 12, 31), 0);
}

TEST(LeapSeconds, AnEntryTakesEffectOnItsOwnDateNotTheDayAfter) {
    // The boundary is the thing most likely to be off by one, so check both sides of every entry.
    for (long long i = 1; i < cdf::leap_second_count(); ++i) {
        const cdf::LeapSecondEntry e = cdf::leap_second_entry(i);
        const st::CivilDate before = st::civil_from_days(st::days_from_civil(e.year, e.month, e.day) - 1);
        EXPECT_NE(cdf::tai_minus_utc_ns(e.year, e.month, e.day),
                  cdf::tai_minus_utc_ns(before.year, before.month, before.day))
            << "entry " << i << " does not change anything on its own date";
        if (e.drift_rate_e7 == 0) {
            EXPECT_EQ(cdf::tai_minus_utc_ns(e.year, e.month, e.day), e.delta_at_e7 * 100)
                << "entry " << i;
        }
    }
}

TEST(LeapSeconds, PreOneNineSevenTwoOffsetDriftsWithinAnEntry) {
    // Inside a drift entry the offset must CHANGE from day to day — that is what makes the era
    // "rubber". A constant answer here would mean the drift term was being dropped.
    const long long a = cdf::tai_minus_utc_ns(1960, 1, 1);
    const long long b = cdf::tai_minus_utc_ns(1960, 6, 1);
    EXPECT_NE(a, b);
    EXPECT_GT(b, a) << "the 1960 entry drifts upward";
}

// ---- exact vs NASA ---------------------------------------------------------------------------------

TEST(LeapSeconds, ExactAndNasaCompatAgreeFrom1972Onward) {
    long long disagreements = 0;
    const long long from = st::days_from_civil(1972, 1, 1);
    const long long to = st::days_from_civil(2035, 1, 1);
    for (long long z = from; z <= to; ++z) {
        const st::CivilDate d = st::civil_from_days(z);
        if (cdf::tai_minus_utc_ns(d.year, d.month, d.day)
            != cdf::tai_minus_utc_ns_nasa_compat(d.year, d.month, d.day)) { ++disagreements; }
    }
    EXPECT_EQ(disagreements, 0)
        << "we must be bit-identical to NASA wherever the offset is a whole second";
}

// The pre-1972 divergence, pinned to its measured extent. This is the "closed, reviewed list"
// that licenses calling tai_minus_utc_ns() exact while still claiming CDF compatibility.
TEST(LeapSeconds, ExactAndNasaCompatDifferOnExactly430DaysBefore1972AlwaysByOneNanosecond) {
    long long differing = 0;
    long long total = 0;
    long long max_abs_delta = 0;
    long long wrong_direction = 0;
    const long long from = st::days_from_civil(1960, 1, 1);
    const long long to = st::days_from_civil(1972, 1, 1);
    for (long long z = from; z <= to; ++z) {
        ++total;
        const st::CivilDate d = st::civil_from_days(z);
        const long long delta = cdf::tai_minus_utc_ns(d.year, d.month, d.day)
                              - cdf::tai_minus_utc_ns_nasa_compat(d.year, d.month, d.day);
        if (delta != 0) {
            ++differing;
            max_abs_delta = std::max(max_abs_delta, delta < 0 ? -delta : delta);
            // NASA truncates toward zero, so its answer is never ABOVE the exact one here.
            if (delta < 0) { ++wrong_direction; }
        }
    }
    EXPECT_EQ(total, 4384);
    EXPECT_EQ(differing, 430) << "the divergence set changed — re-derive it before editing this";
    EXPECT_EQ(max_abs_delta, 1) << "the divergence must only ever be one nanosecond";
    EXPECT_EQ(wrong_direction, 0) << "NASA truncates, so it can only ever be lower";
}

// ---- provenance ------------------------------------------------------------------------------------

TEST(LeapSeconds, ProvenanceIsPresentAndSelfConsistent) {
    EXPECT_EQ(cdf::leap_seconds_url(), "https://cdf.gsfc.nasa.gov/html/CDFLeapSeconds.txt");
    EXPECT_EQ(cdf::leap_seconds_upstream_updated(), "20161025");
    EXPECT_EQ(cdf::leap_seconds_sha256().size(), 64u);
    for (char c : cdf::leap_seconds_sha256()) {
        EXPECT_TRUE((c >= '0' && c <= '9') || (c >= 'a' && c <= 'f')) << "not lowercase hex: " << c;
    }

    const st::CivilDate verified = cdf::leap_seconds_verified_on();
    const st::CivilDate through = cdf::leap_seconds_known_good_through();
    // Both dates must be real ones — civil.hpp already has the predicate for that.
    EXPECT_TRUE(st::is_valid_civil(verified));
    EXPECT_TRUE(st::is_valid_civil(through));
    // The guarantee is derived from the verification date plus IERS's six-month notice, so it
    // must be after it — and the last table entry must be before it, or the table is stale.
    EXPECT_GT(st::days_from_civil(through.year, through.month, through.day),
              st::days_from_civil(verified.year, verified.month, verified.day));
    const cdf::LeapSecondEntry last = cdf::leap_second_entry(cdf::leap_second_count() - 1);
    EXPECT_LT(st::days_from_civil(last.year, last.month, last.day),
              st::days_from_civil(verified.year, verified.month, verified.day));
}

TEST(LeapSeconds, CoversAnswersBothWays) {
    EXPECT_TRUE(cdf::leap_seconds_covers(2026, 8, 28));
    EXPECT_TRUE(cdf::leap_seconds_covers(1972, 1, 1));
    const st::CivilDate limit = cdf::leap_seconds_known_good_through();
    EXPECT_TRUE(cdf::leap_seconds_covers(limit.year, limit.month, limit.day));
    EXPECT_FALSE(cdf::leap_seconds_covers(2099, 1, 1));
    // One day past the guarantee is outside it — the boundary is inclusive.
    const st::CivilDate through = cdf::leap_seconds_known_good_through();
    const st::CivilDate after =
        st::civil_from_days(st::days_from_civil(through.year, through.month, through.day) + 1);
    EXPECT_FALSE(cdf::leap_seconds_covers(after.year, after.month, after.day));
}

// ---- the calendar assumption this file rests on ---------------------------------------------------
//
// tai_minus_utc_ns() feeds space/time/calendar.hpp's julian_day_number() into the drift term.
// NASA's library uses a DIFFERENT formula for the same quantity — the truncating `_JulianDay`
// from its source. The two agree for every year >= 0 and diverge below it, because NASA's needs
// floor division and C truncates toward zero. That divergence is harmless only because CDF
// cannot reach it: CDF_EPOCH counts from year 0 and TT2000 from 1707. Asserted here rather than
// assumed, since every leap-second answer before 1972 depends on it.
TEST(LeapSeconds, NasasJulianDayIsExactEverywhereCdfCanReach) {
    // tai_minus_utc_ns() feeds NASA's own truncating formula into the drift term, so that our
    // arithmetic is bit-compatible with theirs by construction. That formula is only correct for
    // year >= 0 — it needs floor division and C truncates — so this asserts the exactness holds
    // across everything CDF can represent, against the independent exact calendar next door.
    long long disagreements = 0;
    for (int y = 0; y <= 3000; ++y) {
        for (int m = 1; m <= 12; ++m) {
            for (int d : {1, 15, 28}) {
                if (st::nasa_julian_day_at_noon(y, m, d) != st::days_from_civil(y, m, d) + 2440588) {
                    ++disagreements;
                }
            }
        }
    }
    EXPECT_EQ(disagreements, 0)
        << "NASA's formula must be exact wherever CDF can represent a date";

    // Anchors, so a wholesale offset error cannot hide behind mutual agreement.
    EXPECT_EQ(st::nasa_julian_day_at_noon(2000, 1, 1), 2451545);  // J2000
    EXPECT_EQ(st::nasa_julian_day_at_noon(1970, 1, 1), 2440588);  // the Unix epoch
    EXPECT_EQ(st::nasa_julian_day_at_noon(0, 1, 1), 1721060);     // the CDF_EPOCH origin
}
