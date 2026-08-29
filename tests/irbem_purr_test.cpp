// space.irbem's cheatah-facing facade.
//
// The facade exists because a .purr program cannot spell `Igrf<13>::at` (a static factory on a
// class template, with a private constructor), cannot name `Position<Frame::GEO>` (a class
// template whose argument is an enum VALUE, which purr's type-argument grammar has no syntax for),
// and has no `std::span`. These tests hold it to being a THIN wrapper: for the same inputs it must
// return exactly what the typed core returns, bit for bit, or it is a second implementation and
// will drift.
//
// systests/test_irbem.purr is the other half — it proves the surface is reachable from cheatah at
// all, which is the check whose absence let five unusable #include paths survive.
#include <gtest/gtest.h>

#include <cmath>

#include "space/irbem/purr.hpp"

namespace ib = cheatah::space::irbem;
namespace pr = cheatah::space::irbem::purr;

namespace {

// 2015 day 182 at 12:00 UT — mid-2015, inside IGRF-14's published range by a wide margin.
constexpr int kYear = 2015;
constexpr int kDoy = 182;
constexpr double kUt = 43200.0;

pr::Epoch standard() { return pr::epoch_at(kYear, kDoy, kUt); }

}  // namespace

TEST(IrbemPurr, EpochMatchesTheTypedCore) {
    const pr::Epoch e = standard();
    ASSERT_TRUE(pr::epoch_ok(e));
    ASSERT_TRUE(e.model.has_value());

    // The handle must carry the SAME model the typed path builds, not merely a similar one.
    const ib::DateTime dt = ib::date_and_time_from_doy_and_ut(kYear, kDoy, kUt);
    const double decy = ib::decimal_year(dt.year, dt.month, dt.day, dt.hour, dt.minute, dt.second);
    const auto direct = pr::Model::at(decy);
    ASSERT_TRUE(direct.has_value());
    EXPECT_EQ(e.model->g(1, 0), direct->g(1, 0));
    EXPECT_EQ(e.model->g(1, 1), direct->g(1, 1));
    EXPECT_EQ(e.model->h(1, 1), direct->h(1, 1));
}

TEST(IrbemPurr, AnEpochOutsideIgrfIsRefusedRatherThanExtrapolated) {
    // Before the first DGRF and after the last secular-variation prediction. Both decline; neither
    // returns a plausible-looking field for a year the model says nothing about.
    EXPECT_FALSE(pr::epoch_ok(pr::epoch_at(1750, 1, 0.0)));
    EXPECT_FALSE(pr::epoch_ok(pr::epoch_at(2400, 1, 0.0)));
    EXPECT_TRUE(pr::epoch_ok(pr::epoch_at(1900, 1, 0.0)));
}

TEST(IrbemPurr, MakeLstarMatchesTheTypedCore) {
    const pr::Epoch e = standard();
    ASSERT_TRUE(pr::epoch_ok(e));
    auto got = pr::make_lstar(e, 6.6, 0.0, 0.0, 1, 90.0);

    const ib::Position<ib::Frame::GEO> p{cheatah::fixarray::vec3d{6.6, 0.0, 0.0}};
    const ib::Result<ib::api::MagneticCoordinates> want = ib::api::make_lstar(
        *e.model, e.rotations, p, ib::ExternalModel::None, 90.0, ib::DriftShellOptions{});
    ASSERT_EQ(want.status, ib::Status::Ok);

    // Bit-for-bit: the facade re-orders the fields into an array and does nothing else to them.
    EXPECT_EQ(got[0], want.value.lm);
    EXPECT_EQ(got[1], want.value.lstar);
    EXPECT_EQ(got[2], want.value.blocal);
    EXPECT_EQ(got[3], want.value.bmin);
    EXPECT_EQ(got[4], want.value.xj);
    EXPECT_EQ(got[5], want.value.mlt);
    EXPECT_EQ(static_cast<int>(got[6]), static_cast<int>(ib::Status::Ok));
}

TEST(IrbemPurr, MakeLstarReportsWhyItDeclined) {
    // An unusable epoch: the physical slots stay at baddata and the status names the reason, so a
    // caller never reads a number that was never computed.
    const pr::Epoch bad = pr::epoch_at(1750, 1, 0.0);
    auto out = pr::make_lstar(bad, 6.6, 0.0, 0.0, 1, 90.0);
    EXPECT_EQ(static_cast<int>(out[6]), static_cast<int>(ib::Status::OutOfValidityRange));

    // A frame code outside 0..8 is a domain error, not a silent fallback to GEO.
    auto badframe = pr::make_lstar(standard(), 6.6, 0.0, 0.0, 99, 90.0);
    EXPECT_EQ(static_cast<int>(badframe[6]), static_cast<int>(ib::Status::DomainError));
}

TEST(IrbemPurr, MltMatchesTheTypedCore) {
    const pr::Epoch e = standard();
    const double got = pr::get_mlt(e, 6.6, 0.0, 0.0, 1);
    const ib::Result<double> want =
        ib::api::get_mlt(ib::Position<ib::Frame::GEO>{cheatah::fixarray::vec3d{6.6, 0.0, 0.0}},
                         e.rotations);
    ASSERT_EQ(want.status, ib::Status::Ok);
    EXPECT_EQ(got, want.value);
    EXPECT_GE(got, 0.0);
    EXPECT_LT(got, 24.0);

    // A declined call returns a negative hour, which no real MLT can be.
    EXPECT_LT(pr::get_mlt(pr::epoch_at(1750, 1, 0.0), 6.6, 0.0, 0.0, 1), 0.0);
    EXPECT_LT(pr::get_mlt(e, 6.6, 0.0, 0.0, 99), 0.0);
}

TEST(IrbemPurr, CoordTransRoundTrips) {
    const pr::Epoch e = standard();
    auto gsm = pr::coord_trans(e, 6.6, 0.0, 0.0, 1, 2);
    ASSERT_EQ(static_cast<int>(gsm[3]), static_cast<int>(ib::Status::Ok));

    auto back = pr::coord_trans(e, gsm[0], gsm[1], gsm[2], 2, 1);
    ASSERT_EQ(static_cast<int>(back[3]), static_cast<int>(ib::Status::Ok));
    // The rotations compose to the identity, so the round trip closes to roundoff.
    EXPECT_NEAR(back[0], 6.6, 1e-12);
    EXPECT_NEAR(back[1], 0.0, 1e-12);
    EXPECT_NEAR(back[2], 0.0, 1e-12);

    auto refused = pr::coord_trans(e, 6.6, 0.0, 0.0, 1, 99);
    EXPECT_NE(static_cast<int>(refused[3]), static_cast<int>(ib::Status::Ok));
    auto no_epoch = pr::coord_trans(pr::epoch_at(1750, 1, 0.0), 6.6, 0.0, 0.0, 1, 2);
    EXPECT_EQ(static_cast<int>(no_epoch[3]), static_cast<int>(ib::Status::OutOfValidityRange));
}

TEST(IrbemPurr, StatusNamesEveryCode) {
    EXPECT_EQ(pr::status_name(static_cast<int>(ib::Status::Ok)), "Ok");
    EXPECT_EQ(pr::status_name(static_cast<int>(ib::Status::OutOfValidityRange)), "OutOfValidityRange");
    EXPECT_EQ(pr::status_name(static_cast<int>(ib::Status::OpenFieldLine)), "OpenFieldLine");
    EXPECT_EQ(pr::status_name(static_cast<int>(ib::Status::NotConverged)), "NotConverged");
    EXPECT_EQ(pr::status_name(static_cast<int>(ib::Status::ParametersMissing)), "ParametersMissing");
    EXPECT_EQ(pr::status_name(static_cast<int>(ib::Status::DomainError)), "DomainError");
    // A code from a future build names itself rather than throwing.
    EXPECT_EQ(pr::status_name(99), "Unknown");
}
