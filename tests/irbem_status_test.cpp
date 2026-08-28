// Unit tests for space.irbem's status layer — the Result<T> pair, the baddata bridge that is the
// only place -1e31 is allowed to appear, and the per-model validity envelopes.
//
// The envelope table is a transcription of published numbers, so the tests transcribe them AGAIN,
// independently, from IRBEM's `docs/source/api/general_information.rst` kext table. A single typo in
// the header would otherwise be invisible: the header would agree with itself. The two
// transcriptions are compared entry by entry (EnvelopeTableMatchesTheIrbemKextTable), and every
// bounded driver is then exercised from BOTH sides of BOTH its limits with std::nextafter, which is
// the smallest possible step outside — a bound that is off by any amount at all fails.
//
// Boundary convention, asserted rather than assumed: the published limits are stated with <= and
// >=, so a driver exactly ON a bound is inside the envelope and reports Ok. `EveryBoundedDriver...`
// checks the bound itself, and the two neighbouring representable doubles either side of it.
//
// Every check is a status comparison, so every assertion here is `==` — there is no arithmetic to
// need a tolerance, which is the point of keeping the envelope table free of derived quantities.
//
// The status values themselves are read out of runtime tables rather than written as literals, so
// the constexpr functions under test are genuinely executed (and so instrumented) rather than
// folded away at compile time.
#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <optional>
#include <string>
#include <string_view>
#include <type_traits>
#include <vector>

#include "alloc_counter.hpp"
#include "space/irbem/status.hpp"

namespace ib = cheatah::space::irbem;

namespace {

/// Every status, in enumerator order, as runtime data so the functions under test are executed.
constexpr std::array<ib::Status, ib::status_count> all_statuses{
    ib::Status::Ok,        ib::Status::OutOfValidityRange, ib::Status::OpenFieldLine,
    ib::Status::NotConverged, ib::Status::ParametersMissing, ib::Status::DomainError,
};

/// Every model key the table covers, likewise.
constexpr std::array<ib::ExternalModel, ib::model_count> all_models{
    ib::ExternalModel::None,
    ib::ExternalModel::MeadFairfield1975,
    ib::ExternalModel::Tsyganenko1987Short,
    ib::ExternalModel::Tsyganenko1987Long,
    ib::ExternalModel::Tsyganenko1989,
    ib::ExternalModel::OlsonPfitzerQuiet1977,
    ib::ExternalModel::OlsonPfitzerDynamic1988,
    ib::ExternalModel::Tsyganenko1996,
    ib::ExternalModel::OstapenkoMaltsev1997,
    ib::ExternalModel::Tsyganenko2001,
    ib::ExternalModel::Tsyganenko2001Storm,
    ib::ExternalModel::Tsyganenko2004Storm,
    ib::ExternalModel::Alexeev2000,
    ib::ExternalModel::Tsyganenko2007,
    ib::ExternalModel::MeadTsyganenko,
};

/// A `kext` no IRBEM release defines — the unrecognised-key probe.
constexpr auto bogus_model = static_cast<ib::ExternalModel>(200);

/// Drivers that satisfy every published bound of @p m at once: the midpoint of each bounded
/// interval, and zero for a driver the model reads with no published bound. Built at run time from
/// the envelope, so a test that then moves ONE driver is moving it from a known-good vector.
ib::DriverSet nominal_drivers(ib::ExternalModel m) {
    ib::DriverSet d{};
    for (const ib::DriverBound& b : ib::envelope_of(m)) {
        const auto slot = static_cast<std::size_t>(b.driver);
        d[slot] = b.published() ? (b.lo + b.hi) / 2.0 : 0.0;
    }
    return d;
}

/// One published driver interval, transcribed a second time from IRBEM's kext table so that the
/// header's copy has something independent to disagree with.
struct BoundCase {
    ib::ExternalModel model;
    ib::Driver driver;
    double lo;
    double hi;
};

/// Kp is carried as Kp x 10 (OMNI2 scaling), so the published `0 <= Kp <= 9` is 0..90 here.
constexpr std::array<BoundCase, 18> bounded_cases{{
    {ib::ExternalModel::MeadFairfield1975, ib::Driver::Kp, 0.0, 90.0},
    {ib::ExternalModel::Tsyganenko1987Short, ib::Driver::Kp, 0.0, 90.0},
    {ib::ExternalModel::Tsyganenko1987Long, ib::Driver::Kp, 0.0, 90.0},
    {ib::ExternalModel::Tsyganenko1989, ib::Driver::Kp, 0.0, 90.0},
    {ib::ExternalModel::MeadTsyganenko, ib::Driver::Kp, 0.0, 90.0},
    {ib::ExternalModel::OlsonPfitzerDynamic1988, ib::Driver::Dsw, 5.0, 50.0},
    {ib::ExternalModel::OlsonPfitzerDynamic1988, ib::Driver::Vsw, 300.0, 500.0},
    {ib::ExternalModel::OlsonPfitzerDynamic1988, ib::Driver::Dst, -100.0, 20.0},
    {ib::ExternalModel::Tsyganenko1996, ib::Driver::Dst, -100.0, 20.0},
    {ib::ExternalModel::Tsyganenko1996, ib::Driver::Pdyn, 0.5, 10.0},
    {ib::ExternalModel::Tsyganenko1996, ib::Driver::ByIMF, -10.0, 10.0},
    {ib::ExternalModel::Tsyganenko1996, ib::Driver::BzIMF, -10.0, 10.0},
    {ib::ExternalModel::Tsyganenko2001, ib::Driver::Dst, -50.0, 20.0},
    {ib::ExternalModel::Tsyganenko2001, ib::Driver::Pdyn, 0.5, 5.0},
    {ib::ExternalModel::Tsyganenko2001, ib::Driver::ByIMF, -5.0, 5.0},
    {ib::ExternalModel::Tsyganenko2001, ib::Driver::BzIMF, -5.0, 5.0},
    {ib::ExternalModel::Tsyganenko2001, ib::Driver::G1, 0.0, 10.0},
    {ib::ExternalModel::Tsyganenko2001, ib::Driver::G2, 0.0, 10.0},
}};

/// The drivers each model reads, transcribed independently from the same table. A model missing
/// from this list reads nothing.
struct DriverListCase {
    ib::ExternalModel model;
    std::size_t count;
};

constexpr std::array<DriverListCase, ib::model_count> driver_counts{{
    {ib::ExternalModel::None, 0},
    {ib::ExternalModel::MeadFairfield1975, 1},
    {ib::ExternalModel::Tsyganenko1987Short, 1},
    {ib::ExternalModel::Tsyganenko1987Long, 1},
    {ib::ExternalModel::Tsyganenko1989, 1},
    {ib::ExternalModel::OlsonPfitzerQuiet1977, 0},
    {ib::ExternalModel::OlsonPfitzerDynamic1988, 3},
    {ib::ExternalModel::Tsyganenko1996, 4},
    {ib::ExternalModel::OstapenkoMaltsev1997, 4},
    {ib::ExternalModel::Tsyganenko2001, 6},
    {ib::ExternalModel::Tsyganenko2001Storm, 6},
    {ib::ExternalModel::Tsyganenko2004Storm, 10},
    {ib::ExternalModel::Alexeev2000, 5},
    {ib::ExternalModel::Tsyganenko2007, 1},
    {ib::ExternalModel::MeadTsyganenko, 1},
}};

/// The published spatial envelopes, transcribed independently. `inf` means the table publishes no
/// limit of that kind for that model.
struct SpatialCase {
    ib::ExternalModel model;
    double max_r_geo;
    double min_x_gsm;
};

constexpr double inf = std::numeric_limits<double>::infinity();

constexpr std::array<SpatialCase, ib::model_count> spatial_cases{{
    {ib::ExternalModel::None, inf, -inf},
    {ib::ExternalModel::MeadFairfield1975, 17.0, -inf},
    {ib::ExternalModel::Tsyganenko1987Short, 30.0, -inf},
    {ib::ExternalModel::Tsyganenko1987Long, 70.0, -inf},
    {ib::ExternalModel::Tsyganenko1989, 70.0, -inf},
    {ib::ExternalModel::OlsonPfitzerQuiet1977, 15.0, -inf},
    {ib::ExternalModel::OlsonPfitzerDynamic1988, 60.0, -inf},
    {ib::ExternalModel::Tsyganenko1996, 40.0, -inf},
    {ib::ExternalModel::OstapenkoMaltsev1997, inf, -inf},
    {ib::ExternalModel::Tsyganenko2001, inf, -15.0},
    {ib::ExternalModel::Tsyganenko2001Storm, inf, -15.0},
    {ib::ExternalModel::Tsyganenko2004Storm, inf, -15.0},
    {ib::ExternalModel::Alexeev2000, inf, -inf},
    {ib::ExternalModel::Tsyganenko2007, inf, -inf},
    {ib::ExternalModel::MeadTsyganenko, inf, -inf},
}};

/// The next representable double strictly below @p x — the smallest possible step outside a lower
/// bound.
double just_below(double x) { return std::nextafter(x, -inf); }
/// The next representable double strictly above @p x.
double just_above(double x) { return std::nextafter(x, inf); }

}  // namespace

// ---- Status ---------------------------------------------------------------------------------

TEST(IrbemStatus, DescribeCoversEveryEnumerator) {
    std::vector<std::string_view> seen;
    for (const ib::Status s : all_statuses) {
        const std::string_view text = ib::describe(s);
        EXPECT_FALSE(text.empty());
        EXPECT_NE(text, "?") << "enumerator " << static_cast<int>(s) << " has no description";
        seen.push_back(text);
    }
    // Distinct: a copy-paste in the switch would make two statuses indistinguishable in a log.
    for (std::size_t i = 0; i < seen.size(); ++i) {
        for (std::size_t j = i + 1; j < seen.size(); ++j) EXPECT_NE(seen[i], seen[j]);
    }
    EXPECT_EQ(ib::describe(ib::Status::Ok), "ok");
    // Total for a value outside the enumerator list — the case a switch without this arm would
    // fall off the end of.
    EXPECT_EQ(ib::describe(static_cast<ib::Status>(99)), "?");
}

TEST(IrbemStatus, IsOkIsExactlyOk) {
    for (const ib::Status s : all_statuses) {
        EXPECT_EQ(ib::is_ok(s), s == ib::Status::Ok);
    }
}

TEST(IrbemStatus, FirstFailureKeepsTheFirstNonOk) {
    EXPECT_EQ(ib::first_failure(ib::Status::Ok, ib::Status::Ok), ib::Status::Ok);
    // Second reason survives when the first check passed...
    EXPECT_EQ(ib::first_failure(ib::Status::Ok, ib::Status::NotConverged), ib::Status::NotConverged);
    // ...and is DISCARDED when the first failed, which is the whole point: the first reason is the
    // most fundamental one.
    EXPECT_EQ(ib::first_failure(ib::Status::DomainError, ib::Status::NotConverged),
              ib::Status::DomainError);
    EXPECT_EQ(ib::first_failure(ib::Status::OutOfValidityRange, ib::Status::Ok),
              ib::Status::OutOfValidityRange);
}

TEST(IrbemStatus, StatusCodeRoundTrips) {
    for (std::uint32_t code = 0; code < ib::status_count; ++code) {
        const std::optional<ib::Status> s = ib::status_from_code(code);
        ASSERT_TRUE(s.has_value()) << "code " << code << " is inside status_count and must decode";
        EXPECT_EQ(ib::status_code(*s), code);
    }
    for (const ib::Status s : all_statuses) {
        EXPECT_EQ(ib::status_from_code(ib::status_code(s)), s);
    }
    // The wire contract: Ok is zero, so a zeroed status buffer decodes as success.
    EXPECT_EQ(ib::status_code(ib::Status::Ok), 0u);
}

TEST(IrbemStatus, StatusCodeRejectsAnUnknownCode) {
    // The exact boundary of the decoder, from both sides.
    EXPECT_TRUE(ib::status_from_code(ib::status_count - 1).has_value());
    EXPECT_FALSE(ib::status_from_code(ib::status_count).has_value());
    EXPECT_FALSE(ib::status_from_code(0xFFFFFFFFu).has_value());
}

// ---- Result ---------------------------------------------------------------------------------

TEST(IrbemStatus, ResultIsTriviallyCopyableAndSmall) {
    static_assert(std::is_trivially_copyable_v<ib::Result<double>>);
    static_assert(std::is_trivially_copyable_v<ib::Result<float>>);
    static_assert(std::is_aggregate_v<ib::Result<double>>);
    // One alignment unit of tag, and no more: this is what lets a device write an array of them.
    EXPECT_EQ(sizeof(ib::Result<double>), sizeof(double) + alignof(double));
    EXPECT_EQ(sizeof(ib::Result<float>), sizeof(float) + alignof(float));
    EXPECT_EQ(offsetof(ib::Result<double>, status), 0u);
    // Trivially copyable is a claim about memcpy, so memcpy it.
    const ib::Result<double> a{ib::Status::OpenFieldLine, 6.5};
    std::array<unsigned char, sizeof(a)> raw{};
    std::memcpy(raw.data(), &a, sizeof(a));
    ib::Result<double> b{ib::Status::Ok, 0.0};
    std::memcpy(&b, raw.data(), sizeof(b));
    EXPECT_EQ(b.status, ib::Status::OpenFieldLine);
    EXPECT_EQ(b.value, 6.5);  // exactly representable; == is the right comparison
}

TEST(IrbemStatus, ResultOkTracksItsStatus) {
    for (const ib::Status s : all_statuses) {
        const ib::Result<double> r{s, 1.25};
        EXPECT_EQ(r.ok(), s == ib::Status::Ok);
        // The value is ALWAYS returned — that is the contract this whole header exists for.
        EXPECT_EQ(r.value, 1.25);
        const ib::Result<float> f{s, 1.25F};
        EXPECT_EQ(f.ok(), s == ib::Status::Ok);
        EXPECT_EQ(f.value, 1.25F);
    }
}

// ---- the baddata bridge ----------------------------------------------------------------------

TEST(IrbemStatus, BaddataIsTheIrbemSentinel) {
    EXPECT_EQ(ib::baddata, -1e31);
    EXPECT_LT(ib::baddata, 0.0);
    // 1e31 is NOT exactly 10^31 in binary64 — it is the nearest double to it — so the sentinel is
    // a rounded value. Equality against the same literal is still exact, which is what IRBEM's own
    // `.eq. baddata` test relies on.
    EXPECT_NE(ib::baddata, std::nextafter(ib::baddata, 0.0));
}

TEST(IrbemStatus, ToBaddataCollapsesEveryFailure) {
    for (const ib::Status s : all_statuses) {
        const ib::Result<double> r{s, 2.5};
        EXPECT_EQ(ib::to_baddata(r), s == ib::Status::Ok ? 2.5 : ib::baddata);
    }
    // OutOfValidityRange collapses too, even though the value was computed — the documented,
    // deliberate loss at the C boundary.
    EXPECT_EQ(ib::to_baddata(ib::Result<double>{ib::Status::OutOfValidityRange, 2.5}), ib::baddata);
    // A float payload widens exactly (2.5 is a dyadic rational).
    EXPECT_EQ(ib::to_baddata(ib::Result<float>{ib::Status::Ok, 2.5F}), 2.5);
    EXPECT_EQ(ib::to_baddata(ib::Result<float>{ib::Status::NotConverged, 2.5F}), ib::baddata);
}

TEST(IrbemStatus, IsBaddataDetectsOnlyTheSentinel) {
    EXPECT_TRUE(ib::is_baddata(ib::baddata));
    EXPECT_TRUE(ib::is_baddata(-1e31));
    EXPECT_FALSE(ib::is_baddata(0.0));
    EXPECT_FALSE(ib::is_baddata(1e31));
    // The neighbouring representable doubles are NOT the sentinel — the test is exact, not fuzzy.
    EXPECT_FALSE(ib::is_baddata(just_above(ib::baddata)));
    EXPECT_FALSE(ib::is_baddata(just_below(ib::baddata)));
    // A NaN is a different failure from a reported one, and must not be mistaken for it.
    EXPECT_FALSE(ib::is_baddata(std::numeric_limits<double>::quiet_NaN()));
    EXPECT_FALSE(ib::is_baddata(-inf));
}

// ---- the model table ---------------------------------------------------------------------------

TEST(IrbemStatus, ModelNamesCoverEveryKey) {
    std::vector<std::string_view> seen;
    for (std::size_t i = 0; i < all_models.size(); ++i) {
        // The enumerator's value IS the kext key; a caller migrating from IRBEM passes the integer.
        EXPECT_EQ(static_cast<std::size_t>(all_models[i]), i);
        const std::string_view n = ib::name_of(all_models[i]);
        EXPECT_FALSE(n.empty());
        EXPECT_NE(n, "?");
        seen.push_back(n);
    }
    for (std::size_t i = 0; i < seen.size(); ++i) {
        for (std::size_t j = i + 1; j < seen.size(); ++j) EXPECT_NE(seen[i], seen[j]);
    }
    EXPECT_EQ(ib::name_of(ib::ExternalModel::Tsyganenko1996), "T96");
    EXPECT_EQ(ib::name_of(bogus_model), "?");
}

TEST(IrbemStatus, UnrecognisedModelIsADomainError) {
    for (const ib::ExternalModel m : all_models) EXPECT_TRUE(ib::is_recognised(m));
    EXPECT_FALSE(ib::is_recognised(bogus_model));
    // The exact boundary of the table: the last key in, the next key out.
    EXPECT_TRUE(ib::is_recognised(static_cast<ib::ExternalModel>(ib::model_count - 1)));
    EXPECT_FALSE(ib::is_recognised(static_cast<ib::ExternalModel>(ib::model_count)));
    // Every entry point rejects it, and rejects it the same way.
    const ib::DriverSet zeros{};
    EXPECT_EQ(ib::check_validity(bogus_model, zeros), ib::Status::DomainError);
    EXPECT_EQ(ib::check_position(bogus_model, 6.6, 0.0), ib::Status::DomainError);
    EXPECT_EQ(ib::check_parameters(bogus_model, true), ib::Status::DomainError);
    // And the envelope lookup is total rather than out-of-bounds.
    EXPECT_EQ(ib::envelope_of(bogus_model).bound_count, 0u);
    EXPECT_EQ(ib::envelope_of(bogus_model).citation, "unrecognised kext");
}

TEST(IrbemStatus, EnvelopeTableMatchesTheIrbemKextTable) {
    for (const DriverListCase& c : driver_counts) {
        const ib::ValidityEnvelope& env = ib::envelope_of(c.model);
        EXPECT_EQ(env.bound_count, c.count) << ib::name_of(c.model) << " reads a different number "
                                            << "of drivers than the kext table lists";
        EXPECT_FALSE(env.citation.empty()) << ib::name_of(c.model) << " has no citation";
        // The iteration range and the count agree — a range-for sees exactly bound_count entries.
        std::size_t walked = 0;
        for (const ib::DriverBound& b : env) {
            EXPECT_LT(static_cast<std::size_t>(b.driver), ib::named_driver_count);
            ++walked;
        }
        EXPECT_EQ(walked, c.count);
        EXPECT_EQ(static_cast<std::size_t>(env.end() - env.begin()), c.count);
    }
    for (const BoundCase& c : bounded_cases) {
        const std::optional<ib::DriverBound> b = ib::envelope_of(c.model).bound_for(c.driver);
        ASSERT_TRUE(b.has_value()) << ib::name_of(c.model) << " does not read "
                                   << ib::name_of(c.driver);
        EXPECT_TRUE(b->published());
        EXPECT_EQ(b->lo, c.lo) << ib::name_of(c.model) << " " << ib::name_of(c.driver);
        EXPECT_EQ(b->hi, c.hi) << ib::name_of(c.model) << " " << ib::name_of(c.driver);
    }
    for (const SpatialCase& c : spatial_cases) {
        const ib::ValidityEnvelope& env = ib::envelope_of(c.model);
        EXPECT_EQ(env.max_r_geo, c.max_r_geo) << ib::name_of(c.model);
        EXPECT_EQ(env.min_x_gsm, c.min_x_gsm) << ib::name_of(c.model);
    }
    // Exactly one model needs coefficient files.
    std::size_t needing = 0;
    for (const ib::ExternalModel m : all_models) needing += ib::envelope_of(m).needs_coefficient_files ? 1 : 0;
    EXPECT_EQ(needing, 1u);
}

TEST(IrbemStatus, KpBoundsAreInOmniScaling) {
    // The published range is 0 <= Kp <= 9; slot 1 holds Kp x 10, so the stored bound is 0..90.
    EXPECT_EQ(ib::kp_bound.driver, ib::Driver::Kp);
    EXPECT_EQ(ib::kp_bound.lo, 0.0);
    EXPECT_EQ(ib::kp_bound.hi, 90.0);
    EXPECT_TRUE(ib::kp_bound.published());
}

TEST(IrbemStatus, TS05DriversAreUsedButUnbounded) {
    // IRBEM's kext table says of kext 10, 11: "there is no upper or lower limit for those inputs".
    // The header records that as an infinite interval, NOT as an invented number, so `published()`
    // is false and no W value can ever report OutOfValidityRange.
    const ib::ValidityEnvelope& env = ib::envelope_of(ib::ExternalModel::Tsyganenko2004Storm);
    for (const ib::DriverBound& b : env) {
        EXPECT_FALSE(b.published()) << ib::name_of(b.driver) << " claims a bound TS05 never published";
        EXPECT_EQ(b.lo, -inf);
        EXPECT_EQ(b.hi, inf);
    }
    ib::DriverSet d{};
    for (const double w : {-1e6, 0.0, 1e6}) {
        for (std::size_t slot = static_cast<std::size_t>(ib::Driver::W1);
             slot <= static_cast<std::size_t>(ib::Driver::W6); ++slot) {
            d[slot] = w;
        }
        EXPECT_EQ(ib::check_validity(ib::ExternalModel::Tsyganenko2004Storm, d), ib::Status::Ok);
    }
    // A non-finite W is still a domain error — "unbounded" is not "unchecked".
    d[static_cast<std::size_t>(ib::Driver::W3)] = std::numeric_limits<double>::quiet_NaN();
    EXPECT_EQ(ib::check_validity(ib::ExternalModel::Tsyganenko2004Storm, d), ib::Status::DomainError);
}

TEST(IrbemStatus, EnvelopeLookupFindsOnlyTheDriversAModelReads) {
    const ib::ValidityEnvelope& t96 = ib::envelope_of(ib::ExternalModel::Tsyganenko1996);
    EXPECT_TRUE(t96.bound_for(ib::Driver::Pdyn).has_value());
    // T96 does not read Kp, G1 or the W parameters, so a caller's garbage in those slots is
    // irrelevant to it.
    EXPECT_FALSE(t96.bound_for(ib::Driver::Kp).has_value());
    EXPECT_FALSE(t96.bound_for(ib::Driver::G1).has_value());
    EXPECT_FALSE(t96.bound_for(ib::Driver::W1).has_value());
    // ...and that is observable through check_validity, not merely through the table.
    ib::DriverSet d = nominal_drivers(ib::ExternalModel::Tsyganenko1996);
    d[static_cast<std::size_t>(ib::Driver::Kp)] = 1e300;
    d[static_cast<std::size_t>(ib::Driver::W1)] = -1e300;
    EXPECT_EQ(ib::check_validity(ib::ExternalModel::Tsyganenko1996, d), ib::Status::Ok);
    // A model that reads nothing at all finds nothing.
    EXPECT_FALSE(ib::envelope_of(ib::ExternalModel::None).bound_for(ib::Driver::Kp).has_value());
}

TEST(IrbemStatus, MakeEnvelopeFillsTheCountFromTheList) {
    // Built at run time so the builder is exercised rather than folded: the count must come from
    // the list, which is the whole reason the table is not written as a literal aggregate.
    const ib::ValidityEnvelope e = ib::make_envelope(
        "test", {{ib::Driver::Dst, -1.0, 1.0}, {ib::Driver::Pdyn, 0.0, 2.0}}, 5.0, -3.0, true);
    EXPECT_EQ(e.bound_count, 2u);
    EXPECT_EQ(e.citation, "test");
    EXPECT_EQ(e.max_r_geo, 5.0);
    EXPECT_EQ(e.min_x_gsm, -3.0);
    EXPECT_TRUE(e.needs_coefficient_files);
    EXPECT_EQ(e.begin()->driver, ib::Driver::Dst);
    EXPECT_EQ(e.bound_for(ib::Driver::Pdyn)->hi, 2.0);
    // An empty list is a model that reads nothing, not a malformed envelope.
    const ib::ValidityEnvelope none = ib::make_envelope("none", {}, inf, -inf, false);
    EXPECT_EQ(none.bound_count, 0u);
    EXPECT_EQ(none.begin(), none.end());
    // Overflow is dropped rather than written past the array — the guard that makes the fixed
    // capacity safe. Eleven entries into ten slots.
    const ib::ValidityEnvelope over = ib::make_envelope(
        "over",
        {{ib::Driver::Kp, 0.0, 1.0},   {ib::Driver::Dst, 0.0, 1.0},  {ib::Driver::Dsw, 0.0, 1.0},
         {ib::Driver::Vsw, 0.0, 1.0},  {ib::Driver::Pdyn, 0.0, 1.0}, {ib::Driver::ByIMF, 0.0, 1.0},
         {ib::Driver::BzIMF, 0.0, 1.0}, {ib::Driver::G1, 0.0, 1.0},  {ib::Driver::G2, 0.0, 1.0},
         {ib::Driver::G3, 0.0, 1.0},   {ib::Driver::W1, 0.0, 1.0}},
        inf, -inf, false);
    EXPECT_EQ(over.bound_count, ib::max_model_bounds);
    EXPECT_FALSE(over.bound_for(ib::Driver::W1).has_value());
}

// ---- check_validity, from both sides of every bound --------------------------------------------

TEST(IrbemStatus, EveryBoundedDriverIsCheckedFromBothSides) {
    for (const BoundCase& c : bounded_cases) {
        const auto slot = static_cast<std::size_t>(c.driver);
        ib::DriverSet d = nominal_drivers(c.model);
        const std::string label =
            std::string(ib::name_of(c.model)) + "/" + std::string(ib::name_of(c.driver));

        // Comfortably inside.
        d[slot] = (c.lo + c.hi) / 2.0;
        EXPECT_EQ(ib::check_validity(c.model, d), ib::Status::Ok) << label << " midpoint";

        // The DECIDED boundary convention: the published limits are closed, so the bound itself is
        // inside.
        d[slot] = c.lo;
        EXPECT_EQ(ib::check_validity(c.model, d), ib::Status::Ok) << label << " at lo";
        d[slot] = c.hi;
        EXPECT_EQ(ib::check_validity(c.model, d), ib::Status::Ok) << label << " at hi";

        // Just inside, by one representable step.
        d[slot] = just_above(c.lo);
        EXPECT_EQ(ib::check_validity(c.model, d), ib::Status::Ok) << label << " just inside lo";
        d[slot] = just_below(c.hi);
        EXPECT_EQ(ib::check_validity(c.model, d), ib::Status::Ok) << label << " just inside hi";

        // Just outside, by one representable step — the smallest step that must be reported.
        d[slot] = just_below(c.lo);
        EXPECT_EQ(ib::check_validity(c.model, d), ib::Status::OutOfValidityRange)
            << label << " just below lo";
        d[slot] = just_above(c.hi);
        EXPECT_EQ(ib::check_validity(c.model, d), ib::Status::OutOfValidityRange)
            << label << " just above hi";

        // Far outside, so the report is not an artefact of the one-ulp step.
        d[slot] = c.lo - 1.0;
        EXPECT_EQ(ib::check_validity(c.model, d), ib::Status::OutOfValidityRange)
            << label << " below lo";
        d[slot] = c.hi + 1.0;
        EXPECT_EQ(ib::check_validity(c.model, d), ib::Status::OutOfValidityRange)
            << label << " above hi";
    }
}

TEST(IrbemStatus, NominalDriversAreOkForEveryModel) {
    for (const ib::ExternalModel m : all_models) {
        EXPECT_EQ(ib::check_validity(m, nominal_drivers(m)), ib::Status::Ok) << ib::name_of(m);
    }
    // A model with no drivers accepts anything, including a vector of garbage.
    ib::DriverSet junk{};
    junk.fill(1e300);
    EXPECT_EQ(ib::check_validity(ib::ExternalModel::None, junk), ib::Status::Ok);
    EXPECT_EQ(ib::check_validity(ib::ExternalModel::OlsonPfitzerQuiet1977, junk), ib::Status::Ok);
}

TEST(IrbemStatus, NonFiniteDriverIsADomainErrorNotAnOutOfRange) {
    const auto pdyn = static_cast<std::size_t>(ib::Driver::Pdyn);
    for (const double bad : {std::numeric_limits<double>::quiet_NaN(), inf, -inf}) {
        ib::DriverSet d = nominal_drivers(ib::ExternalModel::Tsyganenko1996);
        d[pdyn] = bad;
        EXPECT_EQ(ib::check_validity(ib::ExternalModel::Tsyganenko1996, d), ib::Status::DomainError);
    }
    // Ordering: a NaN in one used slot outranks an out-of-range value in an EARLIER used slot.
    // Dst is listed before Pdyn in T96's envelope, so a per-entry check would return
    // OutOfValidityRange here; the two-pass check reports the more fundamental failure.
    ib::DriverSet d = nominal_drivers(ib::ExternalModel::Tsyganenko1996);
    d[static_cast<std::size_t>(ib::Driver::Dst)] = 1000.0;  // out of [-100, 20]
    d[pdyn] = std::numeric_limits<double>::quiet_NaN();
    EXPECT_EQ(ib::check_validity(ib::ExternalModel::Tsyganenko1996, d), ib::Status::DomainError);
    // A NaN in a slot T96 does NOT read is not T96's problem.
    ib::DriverSet e = nominal_drivers(ib::ExternalModel::Tsyganenko1996);
    e[static_cast<std::size_t>(ib::Driver::W4)] = std::numeric_limits<double>::quiet_NaN();
    EXPECT_EQ(ib::check_validity(ib::ExternalModel::Tsyganenko1996, e), ib::Status::Ok);
}

TEST(IrbemStatus, T96AndT01DifferOnTheSameDrivers) {
    // The concrete reason the envelope is per-model rather than global: T01 was fitted over a
    // narrower range than T96, so the same storm-time vector is inside one and outside the other.
    ib::DriverSet d{};
    d[static_cast<std::size_t>(ib::Driver::Dst)] = -80.0;   // inside T96's -100, outside T01's -50
    d[static_cast<std::size_t>(ib::Driver::Pdyn)] = 8.0;    // inside T96's 10, outside T01's 5
    d[static_cast<std::size_t>(ib::Driver::ByIMF)] = 7.0;   // inside T96's 10, outside T01's 5
    d[static_cast<std::size_t>(ib::Driver::BzIMF)] = -7.0;
    d[static_cast<std::size_t>(ib::Driver::G1)] = 5.0;
    d[static_cast<std::size_t>(ib::Driver::G2)] = 5.0;
    EXPECT_EQ(ib::check_validity(ib::ExternalModel::Tsyganenko1996, d), ib::Status::Ok);
    EXPECT_EQ(ib::check_validity(ib::ExternalModel::Tsyganenko2001, d),
              ib::Status::OutOfValidityRange);
    // ...and the T01-storm variant, which publishes no limits, accepts it.
    EXPECT_EQ(ib::check_validity(ib::ExternalModel::Tsyganenko2001Storm, d), ib::Status::Ok);
}

// ---- check_position ----------------------------------------------------------------------------

TEST(IrbemStatus, PositionEnvelopeIsCheckedFromBothSides) {
    for (const SpatialCase& c : spatial_cases) {
        const std::string label{ib::name_of(c.model)};
        // A point well inside every envelope: geosynchronous, on the dayside.
        EXPECT_EQ(ib::check_position(c.model, 6.6, 6.6), ib::Status::Ok) << label;

        if (std::isfinite(c.max_r_geo)) {
            EXPECT_EQ(ib::check_position(c.model, c.max_r_geo, 0.0), ib::Status::Ok)
                << label << " at the radial limit (closed interval)";
            EXPECT_EQ(ib::check_position(c.model, just_below(c.max_r_geo), 0.0), ib::Status::Ok)
                << label << " just inside the radial limit";
            EXPECT_EQ(ib::check_position(c.model, just_above(c.max_r_geo), 0.0),
                      ib::Status::OutOfValidityRange)
                << label << " one ulp outside the radial limit";
            EXPECT_EQ(ib::check_position(c.model, c.max_r_geo + 1.0, 0.0),
                      ib::Status::OutOfValidityRange)
                << label << " outside the radial limit";
        } else {
            // No published radial limit: arbitrarily far out is Ok as far as radius goes.
            EXPECT_EQ(ib::check_position(c.model, 1e6, 0.0), ib::Status::Ok) << label;
        }

        if (std::isfinite(c.min_x_gsm)) {
            EXPECT_EQ(ib::check_position(c.model, 20.0, c.min_x_gsm), ib::Status::Ok)
                << label << " at the tailward limit (closed interval)";
            EXPECT_EQ(ib::check_position(c.model, 20.0, just_above(c.min_x_gsm)), ib::Status::Ok)
                << label << " just inside the tailward limit";
            EXPECT_EQ(ib::check_position(c.model, 20.0, just_below(c.min_x_gsm)),
                      ib::Status::OutOfValidityRange)
                << label << " one ulp beyond the tailward limit";
            EXPECT_EQ(ib::check_position(c.model, 20.0, c.min_x_gsm - 1.0),
                      ib::Status::OutOfValidityRange)
                << label << " beyond the tailward limit";
        } else {
            EXPECT_EQ(ib::check_position(c.model, 6.6, -1e6), ib::Status::Ok) << label;
        }
    }
}

TEST(IrbemStatus, PositionInsideTheEarthIsADomainError) {
    const ib::ExternalModel m = ib::ExternalModel::Tsyganenko1996;
    // min_r_geo is the WGS84 polar radius over the equatorial radius: below it, a point is
    // underground at EVERY latitude.
    EXPECT_LT(ib::min_r_geo, 1.0);
    EXPECT_GT(ib::min_r_geo, 0.99);
    EXPECT_EQ(ib::check_position(m, ib::min_r_geo, 0.0), ib::Status::Ok);
    EXPECT_EQ(ib::check_position(m, just_above(ib::min_r_geo), 0.0), ib::Status::Ok);
    EXPECT_EQ(ib::check_position(m, just_below(ib::min_r_geo), 0.0), ib::Status::DomainError);
    EXPECT_EQ(ib::check_position(m, 0.5, 0.0), ib::Status::DomainError);
    EXPECT_EQ(ib::check_position(m, 0.0, 0.0), ib::Status::DomainError);
    EXPECT_EQ(ib::check_position(m, -1.0, 0.0), ib::Status::DomainError);
    // A domain error outranks an out-of-range: underground AND beyond the radial limit is still a
    // domain error, because the radius could not be believed in the first place.
    EXPECT_EQ(ib::check_position(m, std::numeric_limits<double>::quiet_NaN(), 0.0),
              ib::Status::DomainError);
    EXPECT_EQ(ib::check_position(m, 6.6, std::numeric_limits<double>::quiet_NaN()),
              ib::Status::DomainError);
    EXPECT_EQ(ib::check_position(m, inf, 0.0), ib::Status::DomainError);
    EXPECT_EQ(ib::check_position(m, 6.6, -inf), ib::Status::DomainError);
}

// ---- check_parameters ---------------------------------------------------------------------------

TEST(IrbemStatus, ParametersAreMissingOnlyForTs07d) {
    for (const ib::ExternalModel m : all_models) {
        EXPECT_EQ(ib::check_parameters(m, true), ib::Status::Ok) << ib::name_of(m);
        const ib::Status without = ib::check_parameters(m, false);
        if (m == ib::ExternalModel::Tsyganenko2007) {
            EXPECT_EQ(without, ib::Status::ParametersMissing);
        } else {
            EXPECT_EQ(without, ib::Status::Ok) << ib::name_of(m) << " should need no coefficients";
        }
    }
    // The composition a caller actually writes: files first, then drivers.
    ib::DriverSet d{};
    EXPECT_EQ(ib::first_failure(ib::check_parameters(ib::ExternalModel::Tsyganenko2007, false),
                                ib::check_validity(ib::ExternalModel::Tsyganenko2007, d)),
              ib::Status::ParametersMissing);
}

// ---- T89's Kp bins -------------------------------------------------------------------------------

TEST(IrbemStatus, T89KpBinsFollowThePublishedIntervals) {
    // Every Kp value OMNI2 can actually hold, in Kp x 10, with the bin the 1989 paper's seven
    // intervals put it in: {0,0+} {1-,1,1+} {2-,2,2+} {3-,3,3+} {4-,4,4+} {5-,5,5+} {>=6-}.
    struct KpCase {
        double kp10;
        int bin;
    };
    constexpr std::array<KpCase, 28> cases{{
        {0, 1},  {3, 1},                     // 0, 0+
        {7, 2},  {10, 2}, {13, 2},           // 1-, 1, 1+
        {17, 3}, {20, 3}, {23, 3},           // 2-, 2, 2+
        {27, 4}, {30, 4}, {33, 4},           // 3-, 3, 3+
        {37, 5}, {40, 5}, {43, 5},           // 4-, 4, 4+
        {47, 6}, {50, 6}, {53, 6},           // 5-, 5, 5+
        {57, 7}, {60, 7}, {63, 7}, {70, 7},  // 6- and up
        {73, 7}, {77, 7}, {80, 7}, {83, 7}, {87, 7}, {90, 7},
    }};
    for (const KpCase& c : cases) {
        EXPECT_EQ(ib::t89_kp_bin(c.kp10), c.bin) << "Kp x 10 = " << c.kp10;
    }
}

TEST(IrbemStatus, T89KpBinThresholdsFallBetweenRepresentableKp) {
    // The thresholds are at Kp x 10 = 5, 15, ..., 55 — BETWEEN the values Kp can take, which is why
    // no real Kp sits on a bin edge. Checked from both sides of every threshold anyway.
    constexpr std::array<double, 6> thresholds{5.0, 15.0, 25.0, 35.0, 45.0, 55.0};
    for (std::size_t i = 0; i < thresholds.size(); ++i) {
        const double t = thresholds[i];
        const int lower = static_cast<int>(i) + 1;
        EXPECT_EQ(ib::t89_kp_bin(just_below(t)), lower) << "just below " << t;
        EXPECT_EQ(ib::t89_kp_bin(t), lower) << "the threshold " << t << " itself belongs below";
        EXPECT_EQ(ib::t89_kp_bin(just_above(t)), lower + 1) << "just above " << t;
    }
    // Saturating and degenerate inputs are still a bin, because the caller gets a coefficient set
    // regardless; check_validity is what reports the input as unusable.
    EXPECT_EQ(ib::t89_kp_bin(-1.0), 1);
    EXPECT_EQ(ib::t89_kp_bin(-1e300), 1);
    EXPECT_EQ(ib::t89_kp_bin(std::numeric_limits<double>::quiet_NaN()), 1);
    EXPECT_EQ(ib::t89_kp_bin(1000.0), 7);
    EXPECT_EQ(ib::t89_kp_bin(inf), 7);
    EXPECT_EQ(ib::t89_kp_bin(-inf), 1);
    // ...and a Kp outside 0..90 IS reported, by the other half of the header.
    ib::DriverSet d{};
    d[static_cast<std::size_t>(ib::Driver::Kp)] = 91.0;
    EXPECT_EQ(ib::check_validity(ib::ExternalModel::Tsyganenko1989, d),
              ib::Status::OutOfValidityRange);
}

// ---- the no-allocation claim -----------------------------------------------------------------

TEST(IrbemStatus, CheckingCostsNoAllocation) {
    // Every check runs on a stack DriverSet against a .rodata envelope. The assertion that catches
    // a per-call workspace is the second measurement, not the first.
    const ib::DriverSet d = nominal_drivers(ib::ExternalModel::Tsyganenko2004Storm);
    volatile int sink = 0;
    sink += static_cast<int>(ib::check_validity(ib::ExternalModel::Tsyganenko2004Storm, d));

    const std::size_t before = cheatah_space_test::allocation_count();
    for (const ib::ExternalModel m : all_models) {
        sink += static_cast<int>(ib::check_validity(m, d));
        sink += static_cast<int>(ib::check_position(m, 6.6, 6.6));
        sink += static_cast<int>(ib::check_parameters(m, false));
        sink += ib::t89_kp_bin(43.0);
        sink += static_cast<int>(ib::envelope_of(m).bound_count);
        sink += static_cast<int>(ib::describe(ib::check_validity(m, d)).size());
        sink += static_cast<int>(ib::to_baddata(ib::Result<double>{ib::Status::Ok, 1.0}));
    }
    EXPECT_EQ(cheatah_space_test::allocation_count(), before);
    EXPECT_NE(sink, -1);  // keep the loop
}
