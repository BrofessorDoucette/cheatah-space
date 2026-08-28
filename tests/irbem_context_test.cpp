// Unit tests for space.irbem's FieldContext — the immutable per-epoch state that exists so the
// mutable Fortran COMMON blocks it replaces cannot be reintroduced.
//
// Two things make the arithmetic here exact rather than approximate:
//
//   * the rotation fixture is a table of SIGNED AXIS PERMUTATIONS (every element 0 or +-1). Those
//     are exactly orthogonal, exactly proper, and every product or round trip through them is
//     exact in binary, so a mismatch is a real mismatch and not a rounding story. A fixture built
//     from real sines and cosines would force a tolerance onto assertions whose whole job is to
//     catch a transposed or misindexed matrix.
//   * the drivers and positions are halves, quarters and the 3-4-12-13 triple.
//
// The one unavoidable transcendental is sin/cos of the tilt; it gets an explicit 1-ulp tolerance at
// its single use, and a second context built at tilt 0 pins the exact case.
//
// The QA gate wants 100% lines AND functions over the header, and clang source-based coverage
// counts each template INSTANTIATION separately, so the frame-pair exercises are folds over the
// nine Cartesian frames (81 ordered pairs) rather than hand-picked samples: no instantiation the
// library can produce is left unexercised, and adding a frame adds no test code.
#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <string>
#include <string_view>
#include <thread>
#include <type_traits>
#include <utility>

#include "space/irbem/context.hpp"

namespace ib = cheatah::space::irbem;
namespace fx = cheatah::fixarray;

using ib::ContextError;
using ib::Driver;
using ib::Frame;

// ---- on the allocation guard ---------------------------------------------------------------------
//
// The repo's rule is that nothing in a hot path touches the heap, "asserted by a test that counts
// global operator new". That counter has to live in exactly ONE translation unit of the test
// binary — the replaceable global allocation functions are a program-wide definition, so a second
// copy is a duplicate-symbol link error rather than a second measurement. Sibling test files in
// this directory already provide it, and it observes the whole process including this file's
// allocations, so this TU deliberately does not define its own.
//
// What is asserted here instead is the structural half, which the counter cannot see: a type that
// owns heap memory needs a non-trivial destructor, so a trivially destructible, trivially copyable
// FieldContext of fixed size cannot own an allocation no matter what the profiler sees on any one
// run. The static_asserts below pin that.

// ---- compile-time surface -----------------------------------------------------------------------

// The reason the type exists at all: it is a value, so N threads and N GPU invocations can hold the
// same bytes, and it is memcpy-able into a uniform block.
static_assert(std::is_trivially_copyable_v<ib::FieldContext>);
static_assert(std::is_standard_layout_v<ib::FieldContext>);
static_assert(std::is_trivially_copyable_v<ib::HotState>);
static_assert(std::is_trivially_copyable_v<ib::ColdState>);
static_assert(std::is_trivially_copyable_v<ib::Epoch>);
// A type that owns heap memory needs a non-trivial destructor to release it, so this is the
// compile-time half of "nothing on the heap": FieldContext cannot own an allocation at all.
static_assert(std::is_trivially_destructible_v<ib::FieldContext>);

// ...and nothing may sneak a pointer, a std::string or a vtable in: a trivially copyable type with
// no padding between its two blocks is exactly a byte image.
static_assert(sizeof(ib::FieldContext) == sizeof(ib::HotState) + sizeof(ib::ColdState));
static_assert(sizeof(ib::HotState) == 208);  // one 64-byte line of scalars + two 72-byte matrices
static_assert(sizeof(ib::FieldContext) <= 4096,
              "must fit a uniform block; Vulkan guarantees only 16 KiB of maxUniformBufferRange");
static_assert(!std::is_polymorphic_v<ib::FieldContext>);

// The only way to obtain one is the factory — construction cannot bypass validation.
static_assert(!std::is_default_constructible_v<ib::FieldContext>);
static_assert(!std::is_constructible_v<ib::FieldContext, ib::Epoch, double, ib::RotationTable,
                                       ib::DriverSet>);

// The driver vector is IRBEM's maginput: 25 wide, 17 of them named.
static_assert(ib::driver_count == 25);
static_assert(ib::named_driver_count == 17);
static_assert(std::tuple_size_v<ib::DriverSet> == 25);
static_assert(static_cast<std::size_t>(Driver::AL) == ib::named_driver_count - 1);
static_assert(ib::cartesian_frame_count == 9);
static_assert(ib::RotationTable{}.size() == 9);

// The slot map is a compile-time fact, so a rotation lookup is an offset and not a branch.
static_assert(ib::cartesian_slot(Frame::GEO) == 0);
static_assert(ib::cartesian_slot(Frame::MAG) == 5);
static_assert(ib::cartesian_slot(Frame::HEE) == 6);
static_assert(ib::cartesian_slot(Frame::HEEQ) == ib::cartesian_frame_count - 1);

// ---- fixtures -----------------------------------------------------------------------------------

namespace {

/// The nine Cartesian frames, in slot order — the loop bound for everything below.
constexpr std::array<Frame, ib::cartesian_frame_count> kCartesianFrames{
    Frame::GEO, Frame::GSM, Frame::GSE,  Frame::SM,  Frame::GEI,
    Frame::MAG, Frame::HEE, Frame::HAE,  Frame::HEEQ};

/// A rotation table of signed axis permutations: each is exactly orthogonal with determinant +1, so
/// products and round trips through them are exact. Physically meaningless on purpose — this file
/// tests the plumbing, and the coordinate transform module tests the angles.
/// @return the fixture table, indexed by ib::cartesian_slot.
ib::RotationTable permutation_table() {
    // Constructed in reading (row) order, which is how fixarray's constructor takes elements.
    return ib::RotationTable{
        fx::mat3d{1, 0, 0, 0, 1, 0, 0, 0, 1},     // GEO  — identity, as the hub must be
        fx::mat3d{0, -1, 0, 1, 0, 0, 0, 0, 1},    // GSM  — +90 deg about z
        fx::mat3d{1, 0, 0, 0, 0, -1, 0, 1, 0},    // GSE  — +90 deg about x
        fx::mat3d{0, 0, 1, 0, 1, 0, -1, 0, 0},    // SM   — +90 deg about y
        fx::mat3d{-1, 0, 0, 0, -1, 0, 0, 0, 1},   // GEI  — 180 deg about z
        fx::mat3d{1, 0, 0, 0, -1, 0, 0, 0, -1},   // MAG  — 180 deg about x
        fx::mat3d{-1, 0, 0, 0, 1, 0, 0, 0, -1},   // HEE  — 180 deg about y
        fx::mat3d{0, 1, 0, -1, 0, 0, 0, 0, 1},    // HAE  — -90 deg about z
        fx::mat3d{0, 0, 1, 1, 0, 0, 0, 1, 0}};    // HEEQ — the cyclic permutation x->y->z->x
}

/// A usable epoch: 2015 doy 100, 06:00 UT. The decimal year is a halved value so it is exact.
/// @return the epoch fixture.
constexpr ib::Epoch good_epoch() { return ib::Epoch{2015.25, 21600.0, 2015, 100}; }

/// A usable driver vector. Every entry distinct and exactly representable, so a misindexed read
/// cannot accidentally match.
/// @return the driver fixture; slot i holds i/4.
ib::DriverSet good_drivers() {
    ib::DriverSet d{};
    for (std::size_t i = 0; i < ib::driver_count; ++i) { d[i] = static_cast<double>(i) / 4.0; }
    return d;
}

/// The standard context — the permutation table and the good epoch and drivers — at a given tilt.
/// The explicit has_value/throw rather than an ASSERT_ macro is what lets a static analyzer see
/// that the dereference is guarded; gtest's early return is control flow it cannot follow into.
/// @param tilt_rad the geodipole tilt, radians.
/// @return the built context; throws if these fixtures were ever rejected, which would be a bug in
///         the fixtures themselves rather than an expected test failure.
ib::FieldContext context_at(double tilt_rad) {
    const ib::ContextResult r =
        ib::make_field_context(good_epoch(), tilt_rad, permutation_table(), good_drivers());
    if (!r.has_value()) { throw std::runtime_error(std::string(ib::describe(r.error()))); }
    return r.value();
}

/// The reason the factory refused these inputs (ContextError::None if it did not).
/// @param e the epoch. @param tilt the tilt in radians.
/// @param t the rotation table. @param d the drivers.
/// @return the reported failure.
ContextError refusal(const ib::Epoch& e, double tilt, const ib::RotationTable& t,
                     const ib::DriverSet& d) {
    const ib::ContextResult r = ib::make_field_context(e, tilt, t, d);
    // The invariant the header states: engaged iff None. Checking it here means every failure case
    // below also checks it.
    EXPECT_EQ(r.error() == ContextError::None, r.has_value());
    return r.error();
}

constexpr double kQuietNaN = std::numeric_limits<double>::quiet_NaN();
constexpr double kInfinity = std::numeric_limits<double>::infinity();

}  // namespace

// ---- the naming tables --------------------------------------------------------------------------

TEST(IrbemContext, EveryNamedDriverHasAName) {
    constexpr std::array<Driver, ib::named_driver_count> kAll{
        Driver::Kp, Driver::Dst, Driver::Dsw, Driver::Vsw, Driver::Pdyn, Driver::ByIMF,
        Driver::BzIMF, Driver::G1, Driver::G2, Driver::G3, Driver::W1, Driver::W2,
        Driver::W3, Driver::W4, Driver::W5, Driver::W6, Driver::AL};
    for (std::size_t i = 0; i < kAll.size(); ++i) {
        const std::string_view name = ib::name_of(kAll[i]);
        EXPECT_FALSE(name.empty());
        EXPECT_NE("?", name) << "driver " << i << " fell through name_of";
        // The enumerators must be the zero-based maginput subscripts, in the table's order — the
        // one fact that makes drivers()[static_cast<size_t>(d)] correct.
        EXPECT_EQ(i, static_cast<std::size_t>(kAll[i])) << name;
    }
    EXPECT_EQ("Pdyn", ib::name_of(Driver::Pdyn));
    EXPECT_EQ("W6", ib::name_of(Driver::W6));
}

// Slots 18..25 of maginput are reserved and have no enumerator. Driver's underlying type is uint8_t,
// so by [dcl.enum]/8 every uint8_t IS a valid Driver value — this is well defined, not a UB probe.
TEST(IrbemContext, NameOfIsTotalForAReservedSlot) {
    // NOLINTNEXTLINE(clang-analyzer-optin.core.EnumCastOutOfRange)
    EXPECT_EQ("?", ib::name_of(static_cast<Driver>(20)));
}

TEST(IrbemContext, EveryErrorHasADescription) {
    constexpr std::array<ContextError, 11> kAll{ContextError::None,
                                                ContextError::YearOutOfRange,
                                                ContextError::DayOfYearOutOfRange,
                                                ContextError::SecondsOfDayOutOfRange,
                                                ContextError::DecimalYearInconsistent,
                                                ContextError::TiltNotFinite,
                                                ContextError::TiltOutOfRange,
                                                ContextError::DriverNotFinite,
                                                ContextError::RotationNotFinite,
                                                ContextError::RotationNotOrthogonal,
                                                ContextError::RotationImproper};
    for (const ContextError e : kAll) {
        const std::string_view text = ib::describe(e);
        EXPECT_FALSE(text.empty());
        EXPECT_NE("?", text) << "error " << static_cast<int>(e) << " fell through describe";
    }
    EXPECT_EQ("ok", ib::describe(ContextError::None));
    // NOLINTNEXTLINE(clang-analyzer-optin.core.EnumCastOutOfRange)
    EXPECT_EQ("?", ib::describe(static_cast<ContextError>(200)));
}

TEST(IrbemContext, CartesianSlotIsABijectionOntoTheTable) {
    std::array<bool, ib::cartesian_frame_count> seen{};
    for (const Frame f : kCartesianFrames) {
        const std::size_t slot = ib::cartesian_slot(f);
        ASSERT_LT(slot, ib::cartesian_frame_count) << ib::name_of(f);
        EXPECT_FALSE(seen[slot]) << "two frames claim slot " << slot;
        seen[slot] = true;
    }
    for (std::size_t i = 0; i < seen.size(); ++i) { EXPECT_TRUE(seen[i]) << "slot " << i; }
    // GEO is the hub, so it must be slot 0 — the whole table is keyed on that.
    EXPECT_EQ(0U, ib::cartesian_slot(Frame::GEO));
}

// ---- what the context carries -------------------------------------------------------------------

TEST(IrbemContext, CarriesTheEpochAndEveryDriver) {
    const ib::FieldContext ctx = context_at(0.0);
    const ib::DriverSet expected = good_drivers();

    EXPECT_EQ(2015.25, ctx.epoch().decimal_year);
    EXPECT_EQ(21600.0, ctx.epoch().seconds_ut);
    EXPECT_EQ(2015, ctx.epoch().year);
    EXPECT_EQ(100, ctx.epoch().day_of_year);

    // All 25 slots survive, reserved ones included: a caller's vector comes back byte for byte.
    EXPECT_EQ(expected, ctx.drivers());
    for (std::size_t i = 0; i < ib::driver_count; ++i) {
        EXPECT_EQ(expected[i], ctx.drivers()[i]) << "slot " << i;
    }
    // ...and the named accessor addresses the same slots.
    EXPECT_EQ(expected[0], ctx.driver(Driver::Kp));
    EXPECT_EQ(expected[4], ctx.driver(Driver::Pdyn));
    EXPECT_EQ(expected[16], ctx.driver(Driver::AL));
}

TEST(IrbemContext, TheHotBlockMirrorsTheColdOneExactly) {
    const ib::FieldContext ctx = context_at(0.0);
    const ib::HotState& hot = ctx.hot();

    // The five hot drivers are copies, and a copy is exactly the kind of thing that goes stale.
    EXPECT_EQ(ctx.driver(Driver::Kp), hot.kp);
    EXPECT_EQ(ctx.driver(Driver::Dst), hot.dst);
    EXPECT_EQ(ctx.driver(Driver::Pdyn), hot.pdyn);
    EXPECT_EQ(ctx.driver(Driver::ByIMF), hot.by_imf);
    EXPECT_EQ(ctx.driver(Driver::BzIMF), hot.bz_imf);

    // The hot GSM pair must be the table entry and its transpose, not two independent guesses.
    const fx::mat3d gsm_to_geo = ctx.rotation_to_geo<Frame::GSM>();
    EXPECT_EQ(gsm_to_geo, hot.gsm_to_geo);
    EXPECT_EQ(fx::transpose(gsm_to_geo), hot.geo_to_gsm);
    EXPECT_EQ((ctx.rotation<Frame::GEO, Frame::GSM>()), hot.geo_to_gsm);
    // ...which is only meaningful if the two directions differ at all.
    EXPECT_FALSE(hot.geo_to_gsm == hot.gsm_to_geo);
}

TEST(IrbemContext, TheTiltTrigonometryIsPrecomputed) {
    // Tilt 0 is exact: no tolerance, no excuse.
    const ib::FieldContext flat = context_at(0.0);
    EXPECT_EQ(0.0, flat.hot().tilt_rad);
    EXPECT_EQ(0.0, flat.hot().sin_tilt);
    EXPECT_EQ(1.0, flat.hot().cos_tilt);

    // A real tilt needs a transcendental, so this is the one tolerance in the file. sin(pi/6) is
    // 0.5 and cos(pi/6) is sqrt(3)/2; a correctly rounded libm gives each to well under one ulp of
    // 1.0, which is 2^-52 ~ 2.2e-16.
    constexpr double kOneUlp = 2.3e-16;
    const double tilt = std::numbers::pi / 6.0;
    const ib::FieldContext leaning = context_at(tilt);
    const ib::HotState& hot = leaning.hot();
    EXPECT_EQ(tilt, hot.tilt_rad);
    EXPECT_NEAR(0.5, hot.sin_tilt, kOneUlp);
    EXPECT_NEAR(std::sqrt(3.0) / 2.0, hot.cos_tilt, kOneUlp);
    // sin^2 + cos^2 = 1 to a couple of ulp — the property a mixed-up sin/cos pair would still pass,
    // which is why the two exact checks above come first.
    EXPECT_NEAR(1.0, (hot.sin_tilt * hot.sin_tilt) + (hot.cos_tilt * hot.cos_tilt), 4.0 * kOneUlp);

    // A negative tilt is ordinary — the dipole leans both ways over a year.
    const ib::FieldContext back = context_at(-tilt);
    const ib::HotState& leaned = back.hot();
    EXPECT_EQ(-hot.sin_tilt, leaned.sin_tilt);
    EXPECT_EQ(hot.cos_tilt, leaned.cos_tilt);
}

// ---- validation ---------------------------------------------------------------------------------

TEST(IrbemContext, AcceptsTheBoundariesOfEveryEpochRange) {
    const ib::RotationTable table = permutation_table();
    const ib::DriverSet drivers = good_drivers();
    // Inclusive bounds are inclusive: 1900.0 is the first IAGA DGRF epoch, doy 366 is a leap year,
    // 86401 s is a day carrying a positive leap second, and +-pi/2 is the definitional tilt bound.
    EXPECT_EQ(ContextError::None,
              refusal(ib::Epoch{1900.0, 0.0, 1900, 1}, 0.0, table, drivers));
    EXPECT_EQ(ContextError::None,
              refusal(ib::Epoch{2101.0, 86401.0, 2100, 366}, 0.0, table, drivers));
    EXPECT_EQ(ContextError::None, refusal(good_epoch(), ib::max_tilt_rad, table, drivers));
    EXPECT_EQ(ContextError::None, refusal(good_epoch(), -ib::max_tilt_rad, table, drivers));
}

TEST(IrbemContext, RejectsAnOutOfRangeEpoch) {
    const ib::RotationTable table = permutation_table();
    const ib::DriverSet drivers = good_drivers();

    EXPECT_EQ(ContextError::YearOutOfRange,
              refusal(ib::Epoch{1899.5, 0.0, 1899, 1}, 0.0, table, drivers));
    EXPECT_EQ(ContextError::YearOutOfRange,
              refusal(ib::Epoch{2101.5, 0.0, 2101, 1}, 0.0, table, drivers));
    EXPECT_EQ(ContextError::DayOfYearOutOfRange,
              refusal(ib::Epoch{2015.0, 0.0, 2015, 0}, 0.0, table, drivers));
    EXPECT_EQ(ContextError::DayOfYearOutOfRange,
              refusal(ib::Epoch{2015.0, 0.0, 2015, 367}, 0.0, table, drivers));
    EXPECT_EQ(ContextError::SecondsOfDayOutOfRange,
              refusal(ib::Epoch{2015.0, -0.5, 2015, 1}, 0.0, table, drivers));
    EXPECT_EQ(ContextError::SecondsOfDayOutOfRange,
              refusal(ib::Epoch{2015.0, 86402.0, 2015, 1}, 0.0, table, drivers));
    EXPECT_EQ(ContextError::SecondsOfDayOutOfRange,
              refusal(ib::Epoch{2015.0, kQuietNaN, 2015, 1}, 0.0, table, drivers));
}

// The check that earns its keep: both spellings of the instant are carried for speed, so they have
// to be checked against each other or a caller silently evaluates IGRF for the wrong year.
TEST(IrbemContext, RejectsADecimalYearThatDisagreesWithTheCalendarYear) {
    const ib::RotationTable table = permutation_table();
    const ib::DriverSet drivers = good_drivers();

    EXPECT_EQ(ContextError::DecimalYearInconsistent,
              refusal(ib::Epoch{2010.5, 21600.0, 2015, 100}, 0.0, table, drivers));
    EXPECT_EQ(ContextError::DecimalYearInconsistent,
              refusal(ib::Epoch{2016.5, 21600.0, 2015, 100}, 0.0, table, drivers));
    EXPECT_EQ(ContextError::DecimalYearInconsistent,
              refusal(ib::Epoch{kQuietNaN, 21600.0, 2015, 100}, 0.0, table, drivers));
    // The closed ends are the year itself and the next new year — both legitimate.
    EXPECT_EQ(ContextError::None, refusal(ib::Epoch{2015.0, 0.0, 2015, 1}, 0.0, table, drivers));
    EXPECT_EQ(ContextError::None,
              refusal(ib::Epoch{2016.0, 86400.0, 2015, 365}, 0.0, table, drivers));
}

TEST(IrbemContext, RejectsAnUnusableTilt) {
    const ib::RotationTable table = permutation_table();
    const ib::DriverSet drivers = good_drivers();

    EXPECT_EQ(ContextError::TiltNotFinite, refusal(good_epoch(), kQuietNaN, table, drivers));
    EXPECT_EQ(ContextError::TiltNotFinite, refusal(good_epoch(), -kInfinity, table, drivers));
    // Beyond pi/2 is not a large tilt, it is a unit or sign error — the angle between an axis and
    // an axis cannot exceed a right angle.
    EXPECT_EQ(ContextError::TiltOutOfRange, refusal(good_epoch(), 2.0, table, drivers));
    EXPECT_EQ(ContextError::TiltOutOfRange, refusal(good_epoch(), -35.0, table, drivers));
}

TEST(IrbemContext, RejectsANonFiniteDriverInAnySlot) {
    const ib::RotationTable table = permutation_table();
    // Every slot is scanned, reserved ones included — a NaN parked in slot 25 is still a NaN that a
    // future model would read.
    for (std::size_t i = 0; i < ib::driver_count; ++i) {
        ib::DriverSet d = good_drivers();
        d[i] = kQuietNaN;
        EXPECT_EQ(ContextError::DriverNotFinite, refusal(good_epoch(), 0.0, table, d))
            << "slot " << i;
        d[i] = kInfinity;
        EXPECT_EQ(ContextError::DriverNotFinite, refusal(good_epoch(), 0.0, table, d))
            << "slot " << i;
    }
    // IRBEM's own baddata sentinel is finite, so a caller who leaves the slots its model ignores
    // filled with it is NOT refused. That is deliberate: range validity belongs to the model that
    // consumes the driver, not to construction.
    ib::DriverSet unset = good_drivers();
    unset[static_cast<std::size_t>(Driver::Vsw)] = -1e31;
    EXPECT_EQ(ContextError::None, refusal(good_epoch(), 0.0, table, unset));
}

TEST(IrbemContext, RejectsAMatrixThatIsNotAProperRotation) {
    const ib::DriverSet drivers = good_drivers();

    // Every slot is checked, not just the ones the inner loop uses.
    for (std::size_t i = 0; i < ib::cartesian_frame_count; ++i) {
        ib::RotationTable bad = permutation_table();
        bad[i] = fx::mat3d::filled(2.0);
        EXPECT_EQ(ContextError::RotationNotOrthogonal, refusal(good_epoch(), 0.0, bad, drivers))
            << "slot " << i;
    }

    // A NaN element must be caught by the finiteness scan BEFORE the tolerance test, because every
    // comparison against a NaN is false and `|NaN - 1| > tol` would quietly accept it.
    ib::RotationTable nan_table = permutation_table();
    nan_table[3](1, 2) = kQuietNaN;
    EXPECT_EQ(ContextError::RotationNotFinite, refusal(good_epoch(), 0.0, nan_table, drivers));
    ib::RotationTable inf_table = permutation_table();
    inf_table[0](0, 0) = kInfinity;
    EXPECT_EQ(ContextError::RotationNotFinite, refusal(good_epoch(), 0.0, inf_table, drivers));

    // A scale, a shear and a rotation that is off by more than the tolerance.
    ib::RotationTable scaled = permutation_table();
    scaled[2] = scaled[2] * 1.5;
    EXPECT_EQ(ContextError::RotationNotOrthogonal, refusal(good_epoch(), 0.0, scaled, drivers));
    ib::RotationTable sheared = permutation_table();
    sheared[5](0, 1) = 0.25;
    EXPECT_EQ(ContextError::RotationNotOrthogonal, refusal(good_epoch(), 0.0, sheared, drivers));

    // A reflection is orthogonal — M^T M is exactly the identity — and would mirror every trace
    // while looking entirely plausible. Only the determinant catches it.
    ib::RotationTable mirror = permutation_table();
    mirror[4] = fx::mat3d{-1, 0, 0, 0, 1, 0, 0, 0, 1};
    EXPECT_EQ(ContextError::RotationImproper, refusal(good_epoch(), 0.0, mirror, drivers));

    // ...and the tolerance is slack enough for a matrix assembled from real transcendentals: a
    // rotation built from sin/cos of an angle with no exact representation must be accepted.
    const double angle = 0.7;
    ib::RotationTable trig = permutation_table();
    trig[7] = fx::mat3d{std::cos(angle), -std::sin(angle), 0.0, std::sin(angle), std::cos(angle),
                        0.0, 0.0, 0.0, 1.0};
    EXPECT_EQ(ContextError::None, refusal(good_epoch(), 0.0, trig, drivers));
}

// The checks run in a documented order, so a caller with several problems gets a stable answer
// rather than one that depends on how the compiler ordered the tests.
TEST(IrbemContext, ReportsTheFirstDefectInTheDocumentedOrder) {
    ib::RotationTable bad = permutation_table();
    bad[0] = fx::mat3d::filled(0.0);
    ib::DriverSet drivers = good_drivers();
    drivers[0] = kQuietNaN;
    // Epoch before tilt before drivers before rotations.
    EXPECT_EQ(ContextError::YearOutOfRange,
              refusal(ib::Epoch{1800.0, 0.0, 1800, 1}, kQuietNaN, bad, drivers));
    EXPECT_EQ(ContextError::TiltNotFinite, refusal(good_epoch(), kQuietNaN, bad, drivers));
    EXPECT_EQ(ContextError::DriverNotFinite, refusal(good_epoch(), 0.0, bad, drivers));
    EXPECT_EQ(ContextError::RotationNotOrthogonal,
              refusal(good_epoch(), 0.0, bad, good_drivers()));
}

// ---- the rotations, over all 81 ordered pairs ---------------------------------------------------

namespace {

/// Exercise every rotation entry point for one ordered pair of frames.
/// @tparam From the source frame. @tparam To the destination frame.
/// @param ctx the context under test, built from permutation_table().
template <Frame From, Frame To>
void exercise_pair(const ib::FieldContext& ctx) {
    const ib::RotationTable table = permutation_table();
    const fx::mat3d forward = ctx.rotation<From, To>();
    const fx::mat3d backward = ctx.rotation<To, From>();

    // Independently of which compile-time branch the accessor took, the answer must be the one the
    // stored table defines: to GEO by the source's matrix, out of GEO by the destination's inverse,
    // which for an orthogonal matrix is its transpose.
    const fx::vec3d v{3.0, -4.0, 12.0};  // 3-4-12-13, exact in binary
    const fx::vec3d via_table = fx::transpose(table[ib::cartesian_slot(To)]) *
                                (table[ib::cartesian_slot(From)] * v);
    EXPECT_EQ(via_table, forward * v)
        << ib::name_of(From) << " -> " << ib::name_of(To);

    // A rotation composed with its reverse is the identity, exactly, for this fixture.
    EXPECT_EQ(fx::mat3d::identity(), backward * forward)
        << ib::name_of(From) << " -> " << ib::name_of(To);

    // The frame tag travels with the value, and the round trip is bit-identical.
    const ib::Position<From> p{v};
    const ib::Position<To> moved = ctx.rotate<To>(p);
    static_assert(std::is_same_v<decltype(moved), const ib::Position<To>>);
    EXPECT_EQ(via_table, moved.v) << ib::name_of(From) << " -> " << ib::name_of(To);
    EXPECT_EQ(p, ctx.rotate<From>(moved)) << ib::name_of(From) << " -> " << ib::name_of(To);

    // ...and a field vector rotates by the same matrix while staying a different type.
    const ib::FieldVector<From> b{fx::vec3d{0.5, -1.5, 2.0}};
    const ib::FieldVector<To> b_moved = ctx.rotate<To>(b);
    static_assert(std::is_same_v<decltype(b_moved), const ib::FieldVector<To>>);
    EXPECT_EQ(forward * b.v, b_moved.v) << ib::name_of(From) << " -> " << ib::name_of(To);
    EXPECT_EQ(b, ctx.rotate<From>(b_moved)) << ib::name_of(From) << " -> " << ib::name_of(To);
    // A rotation preserves length; with these fixtures it does so to the last bit.
    EXPECT_EQ(b.magnitude(), b_moved.magnitude()) << ib::name_of(From) << " -> " << ib::name_of(To);

    if constexpr (From == To) {
        // The self case must be exactly the identity — not "orthogonal to within rounding".
        EXPECT_EQ(fx::mat3d::identity(), forward) << ib::name_of(From);
        EXPECT_EQ(p.v, moved.v) << ib::name_of(From);
    }
}

/// One row of the pair matrix: @p FromIndex against every destination.
/// @tparam FromIndex the source frame's slot. @tparam To the destination slots.
/// @param ctx the context under test. @param dests the destination index pack.
template <std::size_t FromIndex, std::size_t... To>
void exercise_row(const ib::FieldContext& ctx, std::index_sequence<To...> dests) {
    (void)dests;
    (exercise_pair<kCartesianFrames[FromIndex], kCartesianFrames[To]>(ctx), ...);
}

/// Every ordered pair, so every instantiation the library can produce is driven.
/// @tparam From the source slots. @param ctx the context under test. @param sources the index pack.
template <std::size_t... From>
void exercise_all_pairs(const ib::FieldContext& ctx, std::index_sequence<From...> sources) {
    (void)sources;
    (exercise_row<From>(ctx, std::make_index_sequence<ib::cartesian_frame_count>{}), ...);
}

/// Check the stored-matrix accessor for one frame.
/// @tparam F the frame. @param ctx the context under test.
template <Frame F>
void exercise_to_geo(const ib::FieldContext& ctx) {
    const ib::RotationTable table = permutation_table();
    EXPECT_EQ(table[ib::cartesian_slot(F)], ctx.rotation_to_geo<F>()) << ib::name_of(F);
    // The stored matrix IS the frame-to-GEO rotation, so the composed accessor must agree.
    EXPECT_EQ((ctx.rotation<F, Frame::GEO>()), ctx.rotation_to_geo<F>()) << ib::name_of(F);
}

/// Drive @ref exercise_to_geo for every Cartesian frame.
/// @tparam I the frame slots. @param ctx the context under test. @param slots the index pack.
template <std::size_t... I>
void exercise_all_to_geo(const ib::FieldContext& ctx, std::index_sequence<I...> slots) {
    (void)slots;
    (exercise_to_geo<kCartesianFrames[I]>(ctx), ...);
}

}  // namespace

TEST(IrbemContext, EveryStoredRotationIsReadable) {
    const ib::FieldContext ctx = context_at(0.0);
    exercise_all_to_geo(ctx, std::make_index_sequence<ib::cartesian_frame_count>{});
}

TEST(IrbemContext, EveryOrderedFramePairRotatesConsistently) {
    const ib::FieldContext ctx = context_at(0.0);
    exercise_all_pairs(ctx, std::make_index_sequence<ib::cartesian_frame_count>{});
}

// ---- the properties the whole design exists for -------------------------------------------------

// IRBEM chains time points through mutable COMMON blocks, so point k's answer depends on point k-1
// having run first in the same process. Here the state is a value: many threads share one const
// context and get bit-identical answers, which is the property that makes the batch and GPU lanes
// possible at all.
TEST(IrbemContext, IsSafeToShareAcrossThreads) {
    const ib::FieldContext ctx = context_at(0.0);
    constexpr std::size_t kThreads = 4;
    constexpr std::size_t kIterations = 256;

    std::array<fx::vec3d, kThreads> results{};
    std::array<std::thread, kThreads> workers{};
    for (std::size_t t = 0; t < kThreads; ++t) {
        workers[t] = std::thread([&ctx, &results, t]() {
            ib::Position<Frame::GEO> p{fx::vec3d{3.0, -4.0, 12.0}};
            for (std::size_t i = 0; i < kIterations; ++i) {
                // A loop of rotations through four frames and back — the shape of a trace step.
                const ib::Position<Frame::GSM> gsm = ctx.rotate<Frame::GSM>(p);
                const ib::Position<Frame::SM> sm = ctx.rotate<Frame::SM>(gsm);
                const ib::Position<Frame::GEI> gei = ctx.rotate<Frame::GEI>(sm);
                p = ctx.rotate<Frame::GEO>(gei);
            }
            results[t] = p.v;
        });
    }
    for (std::thread& w : workers) { w.join(); }

    for (std::size_t t = 1; t < kThreads; ++t) {
        EXPECT_EQ(results[0], results[t]) << "thread " << t << " disagreed with thread 0";
    }
    // The fixture rotations are permutations of order dividing 4 and the loop count is a multiple
    // of 4, so the exact starting point comes back — evaluation order changed nothing.
    EXPECT_EQ((fx::vec3d{3.0, -4.0, 12.0}), results[0]);
}

// A copied context is byte-identical, which is what makes a memcpy into a uniform block a valid way
// to hand it to a GPU.
TEST(IrbemContext, CopiesAreByteIdentical) {
    const ib::FieldContext ctx = context_at(0.0);
    std::array<std::byte, sizeof(ib::FieldContext)> image{};
    std::memcpy(image.data(), &ctx, sizeof(ib::FieldContext));

    const ib::FieldContext copy = ctx;  // NOLINT(performance-unnecessary-copy-initialization)
    std::array<std::byte, sizeof(ib::FieldContext)> copy_image{};
    std::memcpy(copy_image.data(), &copy, sizeof(ib::FieldContext));

    EXPECT_EQ(image, copy_image);
    EXPECT_EQ(ctx.drivers(), copy.drivers());
    EXPECT_EQ(ctx.hot().gsm_to_geo, copy.hot().gsm_to_geo);
}
