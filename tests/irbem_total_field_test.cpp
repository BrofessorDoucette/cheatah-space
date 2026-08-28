// The internal-plus-external superposition, and the storm sweep it makes possible.
//
// Until TotalFieldT89 existed, `trace_invariant` took a `const Igrf<NMAX>&`, so nothing but the
// internal field could reach a trace and a corpus sweeping Kp/Dst/Pdyn/Bz could not change a single
// number. These tests are the evidence that geomagnetic activity now actually propagates into the
// invariants — which is the claim a reviewer would challenge first, since every Tsyganenko model is
// parameterized BY activity and quiet-time agreement proves almost nothing.
#include <gtest/gtest.h>

#include <cmath>

#include "irbem_domain_corpus.hpp"
#include "space/irbem/api.hpp"
#include "space/irbem/lstar.hpp"
#include "space/irbem/total_field.hpp"

namespace ib = cheatah::space::irbem;
namespace fx = cheatah::fixarray;
namespace corpus = cheatah_space_test;

namespace {
ib::Rotations epoch_rotations(const ib::Igrf<10>& m) {
    const auto r = ib::api::rotations_at(2015, 180, 43200.0, m);
    EXPECT_EQ(ib::Status::Ok, r.status);
    return r.value;
}
}  // namespace

TEST(IrbemTotalField, SuperposesInternalAndExternal) {
    const ib::Igrf<10> igrf = ib::Igrf<10>::at(2015.0).value();
    const ib::Rotations rot = epoch_rotations(igrf);
    const ib::TotalFieldT89<10> total(igrf, rot, 30.0);

    static_assert(ib::GeoFieldModel<ib::TotalFieldT89<10>>,
                  "the whole point is that a tracer can follow it without knowing it is a sum");

    const ib::Position<ib::Frame::GEO> p{fx::vec3d{6.0, 0.0, 0.0}};
    const double internal_only = igrf.evaluate(p).magnitude();
    const double with_external = total.evaluate(p).magnitude();

    EXPECT_EQ(30.0, total.kp_times_ten());
    EXPECT_EQ(igrf.g(1, 0), total.internal().g(1, 0)) << "the internal field must be reachable";
    // The external field is a real contribution at L=6, not a rounding difference — if this were
    // tiny the superposition would be untestable and probably not wired up at all.
    EXPECT_GT(std::abs(with_external - internal_only) / internal_only, 1e-3);
    // ...but it is a PERTURBATION on the internal field out here, not a replacement. A result
    // differing by more than a factor of two would mean the frames were wrong.
    EXPECT_LT(std::abs(with_external - internal_only) / internal_only, 0.5);
}

TEST(IrbemTotalField, ActivityChangesTheInvariants) {
    // The headline: sweep the corpus's four activity regimes and confirm the traced invariants
    // actually move. A library whose answers are identical in quiet and storm conditions is not
    // using the external field, however well it is wired.
    const ib::Igrf<10> igrf = ib::Igrf<10>::at(2015.0).value();
    const ib::Rotations rot = epoch_rotations(igrf);
    const ib::Position<ib::Frame::GEO> p{fx::vec3d{6.0, 0.0, 0.0}};

    std::vector<double> bmin;
    for (const corpus::MagInput& m : corpus::regime_drivers) {
        const ib::TotalFieldT89<10> total(igrf, rot, m.kp * 10.0);
        const auto t = ib::trace_invariant(total, p, 45.0);
        ASSERT_EQ(ib::Status::Ok, t.status) << "Kp = " << m.kp;
        EXPECT_GT(t.value.invariant_i, 0.0);
        bmin.push_back(t.value.b_min);
    }
    ASSERT_EQ(corpus::regime_drivers.size(), bmin.size());

    // Quiet and storm must differ. This is the assertion that would have been impossible before.
    EXPECT_NE(bmin.front(), bmin[2]) << "storm conditions must not reproduce quiet ones";

    // A storm DEPRESSES the equatorial field: the ring current opposes the internal field there.
    // This is the physical signature of a storm, and it falls out of the model rather than being
    // put in by hand — which is why it is worth asserting rather than merely observing.
    EXPECT_LT(bmin[2], bmin.front()) << "the ring current should reduce Bmin during a storm";
}

TEST(IrbemTotalField, KpIsBinnedNotContinuous) {
    // T89 carries SEVEN coefficient sets, not a function of Kp. Two values inside one bin must give
    // bit-identical fields by construction — if they differed, the binning would be broken; if
    // values in DIFFERENT bins agreed, the coefficients would not be being selected at all.
    const ib::Igrf<10> igrf = ib::Igrf<10>::at(2015.0).value();
    const ib::Rotations rot = epoch_rotations(igrf);
    const ib::Position<ib::Frame::GEO> p{fx::vec3d{6.0, 0.0, 0.0}};

    // Kp 6.0 and 8.5 are both in bin 7 (">= 6-"), the open-ended top bin.
    const double a = ib::TotalFieldT89<10>(igrf, rot, 60.0).evaluate(p).magnitude();
    const double b = ib::TotalFieldT89<10>(igrf, rot, 85.0).evaluate(p).magnitude();
    EXPECT_DOUBLE_EQ(a, b) << "same bin must mean the same coefficients";

    // Kp 1 is bin 2; it must differ from the top bin.
    const double quiet = ib::TotalFieldT89<10>(igrf, rot, 10.0).evaluate(p).magnitude();
    EXPECT_NE(quiet, a) << "different bins must select different coefficients";
}

TEST(IrbemTotalField, ReportsWhenTheExternalModelDeclines) {
    // Outside T89's validity the external model declines and the INTERNAL field is returned alone
    // — the best available answer there — rather than a zero or a NaN that would surface as a
    // trace failure hundreds of RK4 steps downstream with no indication of the cause.
    const ib::Igrf<10> igrf = ib::Igrf<10>::at(2015.0).value();
    const ib::Rotations rot = epoch_rotations(igrf);
    const ib::TotalFieldT89<10> total(igrf, rot, 30.0);

    const ib::Position<ib::Frame::GEO> good{fx::vec3d{6.0, 0.0, 0.0}};
    EXPECT_EQ(ib::Status::Ok, total.external_status(good));

    const double nan = std::numeric_limits<double>::quiet_NaN();
    const ib::Position<ib::Frame::GEO> bad{fx::vec3d{nan, 0.0, 0.0}};
    EXPECT_NE(ib::Status::Ok, total.external_status(bad));
    // The value is still finite: internal-only, not a propagated NaN.
    const ib::Position<ib::Frame::GEO> far{fx::vec3d{200.0, 0.0, 0.0}};
    if (total.external_status(far) != ib::Status::Ok) {
        EXPECT_TRUE(std::isfinite(total.evaluate(far).magnitude()));
    }
}
