// Unit tests for space.time — Julian Date / MJD / J2000 / CDF_EPOCH conversions.
// These execute every function of the hand-authored space/time/time.hpp on every instantiation the
// test binary creates, so the QA gate's clang source-based coverage reports 100% lines + functions
// over the space package headers. The .purr systests (systests/test_*.purr) check the same surface
// end-to-end through the cheatah runtime; this checks it in-process for coverage.
//
// The conversions are concept-constrained templates (`TimeInput auto&&`), so each is driven on the
// scalar path (double rvalue + lvalue, and an int to prove Numeric spans the arithmetic types) AND
// on the vectorized path (ndarray<double>) — every instantiation the TU creates is also executed,
// which is what makes llvm-cov see each specialization covered. The header is hand-authored (purrc
// emits none of it), so there is no generated operator<</pretty-print/module_abi glue to touch.
#include <gtest/gtest.h>

#include <string>
#include <vector>

#include "space/space.hpp"

namespace st = cheatah::space::time;
namespace nd = cheatah::ndarray;

// ---- concepts (compile-time surface) ---------------------------------------------------------

// A stand-in for cheatah::datetime's forthcoming broken-down struct: specializing the trait hook
// is exactly how that type will opt in, so the test proves the hook works without the dependency.
struct FakeBrokenDownDatetime {
    int year;  ///< any payload; only the trait matters
};
template <>
inline constexpr bool cheatah::space::time::is_datetime_v<FakeBrokenDownDatetime> = true;

// Numeric spans exactly the arithmetic types; TimeInput is the union of the four input kinds.
static_assert(st::Numeric<double> && st::Numeric<float> && st::Numeric<int>);
static_assert(!st::Numeric<std::string> && !st::Numeric<FakeBrokenDownDatetime>);
static_assert(st::DatetimeScalar<FakeBrokenDownDatetime> && !st::DatetimeScalar<double>);
static_assert(st::NumericArray<nd::basic_ndarray<double>> && st::NumericArray<nd::basic_ndarray<int>>);
static_assert(!st::NumericArray<double> && !st::DatetimeArray<nd::basic_ndarray<double>>);
static_assert(st::TimeInput<double> && st::TimeInput<nd::basic_ndarray<double>>);
static_assert(st::TimeInput<FakeBrokenDownDatetime>);  // via DatetimeScalar — the trait hook
static_assert(!st::TimeInput<std::string> && !st::TimeInput<std::vector<double>>);

// cv/ref-qualified inputs must satisfy the concepts too (the `auto&&` parameters deduce these).
static_assert(st::Numeric<const double&> && st::NumericArray<nd::basic_ndarray<double>&>);

// ---- reference epochs ------------------------------------------------------------------------

// The three epoch constants — usable in constant expressions AND executed at runtime (the runtime
// calls are what coverage counts).
TEST(SpaceTime, Epochs) {
    static_assert(st::jd_unix_epoch() == 2440587.5);   // constexpr-usable
    static_assert(st::jd_j2000() == 2451545.0);
    static_assert(st::cdf_epoch_unix_offset_ms() == 62167219200000.0);

    EXPECT_DOUBLE_EQ(st::jd_unix_epoch(), 2440587.5);
    EXPECT_DOUBLE_EQ(st::jd_j2000(), 2451545.0);
    EXPECT_DOUBLE_EQ(st::cdf_epoch_unix_offset_ms(), 62167219200000.0);

    // The constants relate exactly as the conversions promise.
    EXPECT_DOUBLE_EQ(st::unix_to_jd(0.0), st::jd_unix_epoch());
    EXPECT_DOUBLE_EQ(st::unix_to_cdf_epoch(0.0), st::cdf_epoch_unix_offset_ms());
}

// ---- unix <-> jd -----------------------------------------------------------------------------

// Known instants on exact half-days convert exactly (Unix epoch, ±1 day, J2000).
TEST(SpaceTime, KnownDates) {
    EXPECT_DOUBLE_EQ(st::unix_to_jd(0.0), 2440587.5);
    EXPECT_DOUBLE_EQ(st::unix_to_jd(86400.0), 2440588.5);
    EXPECT_DOUBLE_EQ(st::unix_to_jd(-86400.0), 2440586.5);
    EXPECT_DOUBLE_EQ(st::unix_to_jd(946728000.0), 2451545.0);   // J2000
    EXPECT_DOUBLE_EQ(st::jd_to_unix(2440587.5), 0.0);
    EXPECT_DOUBLE_EQ(st::jd_to_unix(2451545.0), 946728000.0);
}

// Round trips over a spread of real timestamps; a JD double keeps ~50 µs of resolution near the
// present epoch, so the round trip is asserted within 1 ms (see the header's numerical note).
TEST(SpaceTime, UnixJdRoundTrip) {
    for (double s : {0.0, 1.0, 1000.0, 86400.0, 946728000.0, 1700000000.0}) {
        EXPECT_NEAR(st::jd_to_unix(st::unix_to_jd(s)), s, 1e-3);
    }
    // An lvalue drives the `double&` deduction of the forwarding parameter.
    double lv = 946728000.0;
    EXPECT_DOUBLE_EQ(st::unix_to_jd(lv), 2451545.0);
    double jd = 2451545.0;
    EXPECT_DOUBLE_EQ(st::jd_to_unix(jd), 946728000.0);
    // An int scalar proves Numeric spans the arithmetic types (result promotes to double).
    EXPECT_DOUBLE_EQ(st::unix_to_jd(86400), 2440588.5);
}

// ---- modified julian date --------------------------------------------------------------------

// JD <-> MJD is an exact 2400000.5 offset; the Unix legs ride the JD conversions.
TEST(SpaceTime, MjdRoundTrip) {
    EXPECT_DOUBLE_EQ(st::jd_to_mjd(2451545.0), 51544.5);
    EXPECT_DOUBLE_EQ(st::mjd_to_jd(51544.5), 2451545.0);
    EXPECT_DOUBLE_EQ(st::unix_to_mjd(0.0), 40587.0);
    EXPECT_DOUBLE_EQ(st::mjd_to_unix(40587.0), 0.0);
    EXPECT_DOUBLE_EQ(st::unix_to_mjd(946728000.0), 51544.5);

    for (double j : {2440587.5, 2451545.0, 2460000.5}) {
        EXPECT_DOUBLE_EQ(st::mjd_to_jd(st::jd_to_mjd(j)), j);
    }
    for (double s : {0.0, 86400.0, 946728000.0, 1700000000.0}) {
        EXPECT_NEAR(st::mjd_to_unix(st::unix_to_mjd(s)), s, 1e-3);
    }
}

// ---- j2000 offsets ---------------------------------------------------------------------------

// Zero at the epoch itself, exact one-day / one-century offsets, and the negative direction.
TEST(SpaceTime, J2000) {
    EXPECT_DOUBLE_EQ(st::jd_to_j2000_seconds(2451545.0), 0.0);
    EXPECT_DOUBLE_EQ(st::jd_to_j2000_centuries(2451545.0), 0.0);
    EXPECT_DOUBLE_EQ(st::jd_to_j2000_seconds(2451546.0), 86400.0);
    EXPECT_DOUBLE_EQ(st::jd_to_j2000_seconds(2451545.0 + 36525.0), 3155760000.0);
    EXPECT_DOUBLE_EQ(st::jd_to_j2000_centuries(2451545.0 + 36525.0), 1.0);
    EXPECT_DOUBLE_EQ(st::jd_to_j2000_centuries(2451545.0 + 18262.5), 0.5);
    EXPECT_DOUBLE_EQ(st::jd_to_j2000_centuries(2451545.0 - 36525.0), -1.0);
}

// ---- cdf_epoch bridge ------------------------------------------------------------------------

// The ms-since-year-0 scale: exact at the Unix epoch and at J2000's well-known 63113947200000.0.
TEST(SpaceTime, CdfEpoch) {
    EXPECT_DOUBLE_EQ(st::unix_to_cdf_epoch(0.0), 62167219200000.0);
    EXPECT_DOUBLE_EQ(st::cdf_epoch_to_unix(62167219200000.0), 0.0);
    EXPECT_DOUBLE_EQ(st::unix_to_cdf_epoch(946728000.0), 63113947200000.0);
    EXPECT_DOUBLE_EQ(st::cdf_epoch_to_unix(63113947200000.0), 946728000.0);

    // One second is 1000 CDF_EPOCH ms (the ~6.2e13 offset costs ~0.02 ms of double precision).
    EXPECT_NEAR(st::unix_to_cdf_epoch(1.0) - st::cdf_epoch_unix_offset_ms(), 1000.0, 0.1);
    for (double s : {0.0, 1.0, 946728000.0, 1700000000.0}) {
        EXPECT_NEAR(st::cdf_epoch_to_unix(st::unix_to_cdf_epoch(s)), s, 1e-3);
    }
}

// ---- vectorized (ndarray) --------------------------------------------------------------------

// EVERY conversion on the NumericArray path: the same templates broadcast elementwise over a whole
// ndarray of timestamps, and each element must agree exactly with the scalar path.
TEST(SpaceTime, Vectorized) {
    const auto xs = nd::array(std::vector<double>{0.0, 86400.0, 946728000.0});

    auto jd = st::unix_to_jd(xs);                       // const lvalue array instantiation
    EXPECT_DOUBLE_EQ(jd[0], 2440587.5);
    EXPECT_DOUBLE_EQ(jd[1], 2440588.5);
    EXPECT_DOUBLE_EQ(jd[2], 2451545.0);
    EXPECT_DOUBLE_EQ(jd[0], st::unix_to_jd(0.0));       // array path == scalar path
    EXPECT_DOUBLE_EQ(jd[2], st::unix_to_jd(946728000.0));

    auto back = st::jd_to_unix(jd);
    EXPECT_DOUBLE_EQ(back[0], 0.0);
    EXPECT_DOUBLE_EQ(back[2], 946728000.0);

    auto mjd = st::unix_to_mjd(xs);
    EXPECT_DOUBLE_EQ(mjd[0], 40587.0);
    EXPECT_DOUBLE_EQ(mjd[2], 51544.5);
    auto mback = st::mjd_to_unix(mjd);
    EXPECT_DOUBLE_EQ(mback[0], 0.0);

    auto asmjd = st::jd_to_mjd(jd);
    EXPECT_DOUBLE_EQ(asmjd[2], 51544.5);
    auto asjd = st::mjd_to_jd(asmjd);
    EXPECT_DOUBLE_EQ(asjd[2], 2451545.0);

    auto secs = st::jd_to_j2000_seconds(jd);
    EXPECT_DOUBLE_EQ(secs[2], 0.0);                     // J2000 is its own origin
    EXPECT_DOUBLE_EQ(secs[0], (2440587.5 - 2451545.0) * 86400.0);
    auto cent = st::jd_to_j2000_centuries(jd);
    EXPECT_DOUBLE_EQ(cent[2], 0.0);
    EXPECT_DOUBLE_EQ(cent[0], (2440587.5 - 2451545.0) / 36525.0);

    auto cdf = st::unix_to_cdf_epoch(xs);
    EXPECT_DOUBLE_EQ(cdf[0], 62167219200000.0);
    EXPECT_DOUBLE_EQ(cdf[2], 63113947200000.0);
    auto cback = st::cdf_epoch_to_unix(cdf);
    EXPECT_DOUBLE_EQ(cback[0], 0.0);
    EXPECT_DOUBLE_EQ(cback[2], 946728000.0);
}

// An rvalue ndarray drives the value-deduced instantiation (the `auto&&` forwarding path a chained
// call like unix_to_mjd's jd_to_mjd(unix_to_jd(xs)) takes internally).
TEST(SpaceTime, VectorizedRvalue) {
    auto jd = st::unix_to_jd(nd::array(std::vector<double>{0.0, 946728000.0}));
    EXPECT_DOUBLE_EQ(jd[1], 2451545.0);
    auto s = st::jd_to_unix(nd::array(std::vector<double>{2440587.5}));
    EXPECT_DOUBLE_EQ(s[0], 0.0);
}
