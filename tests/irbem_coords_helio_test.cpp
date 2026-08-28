// Unit tests for space.irbem's heliospheric transforms — HEE, HAE, HEEQ and the GSE bridge.
//
// The reference is Fränz & Harper 2002 (Planet. Space Sci. 50, 217; corrected version 2002-03-12),
// whose Table 8 tabulates ONE spacecraft vector in fifteen coordinate systems for 1996-08-28
// 16:46:00 TT — the same epoch Hapgood's own software uses as its reference set. Three of those
// rows (HAE_D, HEE_D, HEEQ_D) pin down exactly the three rotations this file implements, and a
// fourth (GSE_D) pins the half turn. They are reproduced verbatim below and are the primary
// acceptance criterion.
//
// Two things Table 8 deliberately does NOT test, because its vectors are all *geocentric* ("the
// vector is still geocentric since we did not apply a translation", F&H §A.2):
//
//   * the origin shift between GSE and HEE, and
//   * the fact that a field vector must not receive it.
//
// Those are tested here from the definitions instead, with exactly-representable values so the
// assertions are `==` rather than approximate: a spacecraft sitting on the Sun has GSE position
// (r0, 0, 0) and must land on the HEE origin, while a field measured there must not move at all.
//
// Coverage note: clang counts each template INSTANTIATION separately, so every rotational transform
// is driven with BOTH Position and FieldVector, through a helper that takes the vector kind as a
// template template parameter.
#include <gtest/gtest.h>

#include <cmath>

#include "space/irbem/coords_helio.hpp"
#include "space/irbem/frames.hpp"

namespace ib = cheatah::space::irbem;
namespace fx = cheatah::fixarray;

using ib::FieldVector;
using ib::Frame;
using ib::Position;

// ---- the Fränz & Harper reference epoch and vectors ---------------------------------------------

// 1996-08-28 16:46:00 TT = JD 2450324.19861111 = MJD 50323.69861111 (F&H Table 8 caption).
// Written as a whole day plus an exact count of seconds so the literal carries no transcription
// slip; 16h46m = 60360 s.
constexpr double kReferenceMjdTt = 50323.0 + (60360.0 / 86400.0);

// F&H Table 8, in units of Earth radii. Only the four rows this file is responsible for.
constexpr fx::vec3d kTable8Hae{-5.7864918, -3.0028771, 3.3908764};
constexpr fx::vec3d kTable8Hee{-4.0378470, -5.1182566, 3.3908764};
constexpr fx::vec3d kTable8Heeq{-4.4132668, -5.1924440, 2.7496187};
constexpr fx::vec3d kTable8Gse{4.0378470, 5.1182566, 3.3908764};

// F&H print Table 8 to 8 significant figures, so a component of magnitude ~5 is quoted to 5e-8.
// Everything here agrees with them an order better than that; the cap is set at the printing
// resolution so the test measures OUR error rather than their rounding.
constexpr double kTable8Tolerance = 2e-7;

// F&H's own text for this epoch: "for HEEQ_D we use θ = 259.89919 (eqn.17) ... we use the ecliptic
// longitude of the Earth λ_geo = -24.302838 (eqn.36)". Both are printed to 8 figures.
constexpr double kPublishedCentralMeridianDeg = 259.89919;
constexpr double kPublishedEarthLongitudeDeg = 360.0 - 24.302838;
constexpr double kPublishedAngleTolerance = 1e-5;

/// The geometry under test, built once — every test uses the same epoch unless it says otherwise.
const ib::HelioGeometry& reference_geometry() {
    static const ib::HelioGeometry geometry = ib::helio_geometry(kReferenceMjdTt);
    return geometry;
}

/// Componentwise comparison with an absolute cap, reported with the offending axis named.
/// @param actual the computed vector. @param expected the reference vector.
/// @param tolerance the absolute cap per component. @param what a label for the failure message.
void expect_vector_near(const fx::vec3d& actual, const fx::vec3d& expected, double tolerance,
                        const char* what) {
    for (std::size_t i = 0; i < 3; ++i) {
        EXPECT_NEAR(expected[i], actual[i], tolerance)
            << what << " component " << i << " (0=x, 1=y, 2=z)";
    }
}

// ---- compile-time surface -----------------------------------------------------------------------

// The frame tag is load-bearing: a transform must refuse a vector from the wrong frame. Deduction
// binds V from the argument, then V<Frame::HAE> must match it exactly, so HAE2HEE cannot be handed
// a GSE vector at all.
template <class T>
concept Hae2HeeAccepts = requires(T value, const ib::HelioGeometry& geometry) {
    ib::HAE2HEE(value, geometry);
};
static_assert(Hae2HeeAccepts<Position<Frame::HAE>>);
static_assert(Hae2HeeAccepts<FieldVector<Frame::HAE>>);
static_assert(!Hae2HeeAccepts<Position<Frame::GSE>>);
static_assert(!Hae2HeeAccepts<Position<Frame::HEE>>);
static_assert(!Hae2HeeAccepts<FieldVector<Frame::HEEQ>>);

// The concept admits exactly the two frame-tagged vector templates.
static_assert(ib::FrameTagged<Position>);
static_assert(ib::FrameTagged<FieldVector>);

// The output frame is fixed by the function, not by the caller.
static_assert(std::is_same_v<decltype(ib::HAE2HEEQ(Position<Frame::HAE>{}, reference_geometry())),
                             Position<Frame::HEEQ>>);
static_assert(
    std::is_same_v<decltype(ib::GSE2HEE(FieldVector<Frame::GSE>{}, reference_geometry())),
                   FieldVector<Frame::HEE>>);

// ---- the published constants --------------------------------------------------------------------

TEST(IrbemCoordsHelio, PublishedConstants) {
    // IAU 2012 Resolution B2 defines the astronomical unit as exactly 149 597 870 700 m.
    EXPECT_EQ(149597870.7, ib::astronomical_unit_km());
    // 7.25 deg is the conventional solar-equator inclination of F&H eqn. 14.
    EXPECT_EQ(7.25, ib::solar_equator_inclination_deg());
    // ... and the AU expressed in the module's own length unit (Re = 6371.2 km).
    EXPECT_NEAR(23480.3287763, ib::au_in_earth_radii(), 1e-7);
    EXPECT_EQ(ib::astronomical_unit_km() / 6371.2, ib::au_in_earth_radii());
}

// ---- the epoch geometry --------------------------------------------------------------------------

TEST(IrbemCoordsHelio, GeometryMatchesTheFranzHarperWorkedExample) {
    const ib::HelioGeometry& g = reference_geometry();

    // The two angles F&H print in prose for this epoch.
    EXPECT_NEAR(kPublishedEarthLongitudeDeg, g.earth_longitude_deg, kPublishedAngleTolerance);
    EXPECT_NEAR(kPublishedCentralMeridianDeg, g.solar_central_meridian_deg,
                kPublishedAngleTolerance);

    // Ω from eqn. 14 at T0 = -0.0334237..., and the distance from eqn. 36. Late August is seven
    // weeks past aphelion, so the Earth is still comfortably beyond 1 AU.
    EXPECT_NEAR(75.7133070, g.solar_node_deg, 1e-7);
    EXPECT_NEAR(1.00993402, g.sun_earth_distance_au, 1e-8);
    EXPECT_GT(g.sun_earth_distance_au, 1.0);

    // The two distance fields must describe one distance.
    EXPECT_EQ(g.sun_earth_distance_au * ib::au_in_earth_radii(), g.sun_earth_distance_re);

    // Aberration moves the apparent sub-Earth point by exactly 20 arcsec (F&H §3.2.1 / App. A.3).
    // The two longitudes are each folded into [0, 360) from a raw value ~1.1e3 deg, so their
    // difference carries the ulp of THAT magnitude (2.3e-13 deg), not of the 5.6e-3 deg answer.
    EXPECT_NEAR(-20.0 / 3600.0, g.apparent_earth_longitude_deg - g.earth_longitude_deg, 1e-12);

    // Every angle is folded into [0, 360).
    for (const double angle : {g.earth_longitude_deg, g.apparent_earth_longitude_deg,
                               g.solar_node_deg, g.solar_central_meridian_deg}) {
        EXPECT_GE(angle, 0.0);
        EXPECT_LT(angle, 360.0);
    }
}

TEST(IrbemCoordsHelio, CentralMeridianSatisfiesEquation17) {
    const ib::HelioGeometry& g = reference_geometry();
    // Independently of how it was computed, θ must satisfy tan θ = cos(i)·tan(λ_app − Ω), with the
    // quadrant of λ_app − Ω. Comparing sin/cos rather than the tangents keeps the identity finite.
    const double to_rad = std::acos(-1.0) / 180.0;
    const double from_node = (g.apparent_earth_longitude_deg - g.solar_node_deg) * to_rad;
    const double theta = g.solar_central_meridian_deg * to_rad;
    const double cos_i = std::cos(ib::solar_equator_inclination_deg() * to_rad);
    // sin θ · cos(λ−Ω) = cos i · cos θ · sin(λ−Ω) is the tangent identity, cleared of poles.
    EXPECT_NEAR(std::sin(theta) * std::cos(from_node),
                cos_i * std::cos(theta) * std::sin(from_node), 1e-15);
}

TEST(IrbemCoordsHelio, GeometryTracksTheEpoch) {
    // A stub returning fixed angles would pass every Table 8 assertion above if it were tuned to
    // this one epoch. So: sixty days later the Earth has moved ~59 deg along its orbit, the node
    // has crept by 1.397 deg/century, and the distance has changed by a per-cent-level amount.
    const ib::HelioGeometry& g = reference_geometry();
    const ib::HelioGeometry later = ib::helio_geometry(kReferenceMjdTt + 60.0);

    const double advance = std::fmod(later.earth_longitude_deg - g.earth_longitude_deg + 360.0,
                                     360.0);
    EXPECT_GT(advance, 55.0);
    EXPECT_LT(advance, 63.0);
    EXPECT_NEAR(60.0 / 36525.0 * 1.397, later.solar_node_deg - g.solar_node_deg, 1e-12);
    EXPECT_NE(g.sun_earth_distance_au, later.sun_earth_distance_au);

    // An epoch after J2000 exercises the other branch of the [0, 360) fold, where the raw mean
    // longitude is already positive.
    const ib::HelioGeometry modern = ib::helio_geometry(61000.0);  // 2025-11-21
    EXPECT_GE(modern.earth_longitude_deg, 0.0);
    EXPECT_LT(modern.earth_longitude_deg, 360.0);
    // Late November, six weeks short of perihelion: already well inside 1 AU.
    EXPECT_NEAR(0.98794, modern.sun_earth_distance_au, 1e-5);
}

// ---- the rotations, against Table 8 ---------------------------------------------------------------

TEST(IrbemCoordsHelio, HaeToHeeMatchesTable8) {
    const ib::HelioGeometry& g = reference_geometry();
    expect_vector_near(ib::HAE2HEE(Position<Frame::HAE>{kTable8Hae}, g).v, kTable8Hee,
                       kTable8Tolerance, "HAE2HEE position");
    expect_vector_near(ib::HAE2HEE(FieldVector<Frame::HAE>{kTable8Hae}, g).v, kTable8Hee,
                       kTable8Tolerance, "HAE2HEE field");

    // The rotation is about the ecliptic pole, so Z is untouched — exactly, not nearly.
    EXPECT_EQ(kTable8Hae[2], ib::HAE2HEE(Position<Frame::HAE>{kTable8Hae}, g).v[2]);
}

TEST(IrbemCoordsHelio, HeeToHaeMatchesTable8) {
    const ib::HelioGeometry& g = reference_geometry();
    expect_vector_near(ib::HEE2HAE(Position<Frame::HEE>{kTable8Hee}, g).v, kTable8Hae,
                       kTable8Tolerance, "HEE2HAE position");
    expect_vector_near(ib::HEE2HAE(FieldVector<Frame::HEE>{kTable8Hee}, g).v, kTable8Hae,
                       kTable8Tolerance, "HEE2HAE field");
}

TEST(IrbemCoordsHelio, HaeToHeeqMatchesTable8) {
    const ib::HelioGeometry& g = reference_geometry();
    expect_vector_near(ib::HAE2HEEQ(Position<Frame::HAE>{kTable8Hae}, g).v, kTable8Heeq,
                       kTable8Tolerance, "HAE2HEEQ position");
    expect_vector_near(ib::HAE2HEEQ(FieldVector<Frame::HAE>{kTable8Hae}, g).v, kTable8Heeq,
                       kTable8Tolerance, "HAE2HEEQ field");
}

TEST(IrbemCoordsHelio, HeeqToHaeMatchesTable8) {
    const ib::HelioGeometry& g = reference_geometry();
    expect_vector_near(ib::HEEQ2HAE(Position<Frame::HEEQ>{kTable8Heeq}, g).v, kTable8Hae,
                       kTable8Tolerance, "HEEQ2HAE position");
    expect_vector_near(ib::HEEQ2HAE(FieldVector<Frame::HEEQ>{kTable8Heeq}, g).v, kTable8Hae,
                       kTable8Tolerance, "HEEQ2HAE field");
}

TEST(IrbemCoordsHelio, GseToHeeIsExactlyTheHalfTurnOfTable8) {
    const ib::HelioGeometry& g = reference_geometry();
    // A half turn is a sign flip, so for a FIELD (no origin shift) Table 8's GSE row must map onto
    // its HEE row bit for bit. This is the one assertion in the file that can be exact against a
    // published number, and it is exact because no transcendental is involved.
    const FieldVector<Frame::HEE> hee = ib::GSE2HEE(FieldVector<Frame::GSE>{kTable8Gse}, g);
    EXPECT_EQ(kTable8Hee[0], hee.v[0]);
    EXPECT_EQ(kTable8Hee[1], hee.v[1]);
    EXPECT_EQ(kTable8Hee[2], hee.v[2]);

    const FieldVector<Frame::GSE> back = ib::HEE2GSE(FieldVector<Frame::HEE>{kTable8Hee}, g);
    EXPECT_EQ(kTable8Gse[0], back.v[0]);
    EXPECT_EQ(kTable8Gse[1], back.v[1]);
    EXPECT_EQ(kTable8Gse[2], back.v[2]);
}

// ---- the origin shift: the defect this file exists to prevent --------------------------------------

TEST(IrbemCoordsHelio, PositionTranslatesButFieldDoesNot) {
    const ib::HelioGeometry& g = reference_geometry();
    const double r0 = g.sun_earth_distance_re;

    // The GSE origin IS the Earth, whose HEE position is (r0, 0, 0) by construction.
    const Position<Frame::HEE> earth = ib::GSE2HEE(Position<Frame::GSE>{}, g);
    EXPECT_EQ(r0, earth.v[0]);
    EXPECT_EQ(0.0, earth.v[1]);
    EXPECT_EQ(0.0, earth.v[2]);

    // A FIELD measured at the Earth has no position to translate: zero in, zero out. If the origin
    // shift ever leaks into the field overload, this is where it shows — as 23 713 nT of nothing.
    const FieldVector<Frame::HEE> field = ib::GSE2HEE(FieldVector<Frame::GSE>{}, g);
    EXPECT_EQ(0.0, field.v[0]);
    EXPECT_EQ(0.0, field.v[1]);
    EXPECT_EQ(0.0, field.v[2]);

    // Symmetrically: the HEE origin IS the Sun, whose GSE position is (r0, 0, 0) — +X_GSE points at
    // the Sun, so both origins sit on each other's +X axis at the same distance.
    const Position<Frame::GSE> sun = ib::HEE2GSE(Position<Frame::HEE>{}, g);
    EXPECT_EQ(r0, sun.v[0]);
    EXPECT_EQ(0.0, sun.v[1]);
    EXPECT_EQ(0.0, sun.v[2]);

    const FieldVector<Frame::GSE> sun_field = ib::HEE2GSE(FieldVector<Frame::HEE>{}, g);
    EXPECT_EQ(0.0, sun_field.v[0]);
    EXPECT_EQ(0.0, sun_field.v[1]);
    EXPECT_EQ(0.0, sun_field.v[2]);

    // And a spacecraft parked ON the Sun — GSE (r0, 0, 0) — lands on the HEE origin exactly.
    const Position<Frame::HEE> at_sun = ib::GSE2HEE(Position<Frame::GSE>{fx::vec3d{r0, 0.0, 0.0}}, g);
    EXPECT_EQ(0.0, at_sun.v[0]);
    EXPECT_EQ(0.0, at_sun.v[1]);
    EXPECT_EQ(0.0, at_sun.v[2]);
}

TEST(IrbemCoordsHelio, PositionAndFieldDifferByExactlyTheSunEarthVector) {
    const ib::HelioGeometry& g = reference_geometry();
    // Quarters and halves: r0 is ~2.4e4, whose ulp is 2^-38, and 6.5 / 2.25 / 0.75 are all exact
    // multiples of that, so `r0 - x` is representable and the arithmetic below is exact.
    const fx::vec3d sample{6.5, -2.25, 0.75};

    const Position<Frame::HEE> position = ib::GSE2HEE(Position<Frame::GSE>{sample}, g);
    const FieldVector<Frame::HEE> field = ib::GSE2HEE(FieldVector<Frame::GSE>{sample}, g);

    // They differ in X by the Sun-Earth distance and nowhere else.
    EXPECT_EQ(g.sun_earth_distance_re, position.v[0] - field.v[0]);
    EXPECT_EQ(field.v[1], position.v[1]);
    EXPECT_EQ(field.v[2], position.v[2]);

    // The rotation itself is the half turn, exactly.
    EXPECT_EQ(-6.5, field.v[0]);
    EXPECT_EQ(2.25, field.v[1]);
    EXPECT_EQ(0.75, field.v[2]);
}

TEST(IrbemCoordsHelio, GseHeePositionRoundTripIsExact) {
    const ib::HelioGeometry& g = reference_geometry();
    const fx::vec3d sample{6.5, -2.25, 0.75};

    const Position<Frame::GSE> gse{sample};
    EXPECT_EQ(gse, ib::HEE2GSE(ib::GSE2HEE(gse, g), g));

    const Position<Frame::HEE> hee{sample};
    EXPECT_EQ(hee, ib::GSE2HEE(ib::HEE2GSE(hee, g), g));

    const FieldVector<Frame::GSE> gse_field{sample};
    EXPECT_EQ(gse_field, ib::HEE2GSE(ib::GSE2HEE(gse_field, g), g));

    const FieldVector<Frame::HEE> hee_field{sample};
    EXPECT_EQ(hee_field, ib::GSE2HEE(ib::HEE2GSE(hee_field, g), g));
}

TEST(IrbemCoordsHelio, GseHeePositionRoundTripHoldsForAdversarialValues) {
    const ib::HelioGeometry& g = reference_geometry();
    // Values chosen so that `r0 - x` is NOT representable and the subtraction really rounds. The
    // floor is then one ulp of r0 (~2.4e4 Re, ulp 3.6e-12), not one ulp of the coordinate — the
    // price of expressing a 6 Re position relative to an origin 23 000 Re away. Measured worst
    // component here: 7.3e-13 Re, i.e. 4.6 nm. The 1e-12 cap is that floor, not a guess.
    const fx::vec3d nasty{6.700000000000001, -2.3333333333333335, 0.1234567890123};

    const Position<Frame::GSE> gse{nasty};
    expect_vector_near(ib::HEE2GSE(ib::GSE2HEE(gse, g), g).v, nasty, 1e-12, "GSE->HEE->GSE");

    const Position<Frame::HEE> hee{nasty};
    expect_vector_near(ib::GSE2HEE(ib::HEE2GSE(hee, g), g).v, nasty, 1e-12, "HEE->GSE->HEE");

    // A field vector has no origin to subtract from, so its round trip is exact whatever the value.
    const FieldVector<Frame::GSE> field{nasty};
    EXPECT_EQ(field, ib::HEE2GSE(ib::GSE2HEE(field, g), g));
}

// ---- round trips through the rotations -------------------------------------------------------------

// Rotating out and back multiplies by M^T·M = I to within the conditioning of a 3x3 orthonormal
// product: three multiply-adds per component, so ~4 ulp of the largest component (~1e-15 here). The
// task's 1e-12 target has three orders of headroom over that.
constexpr double kRoundTripTolerance = 1e-14;

/// Drive one HAE-centred round trip for one vector kind. Templated on the kind so that both
/// instantiations of every transform are exercised.
/// @tparam V the vector kind, Position or FieldVector.
/// @param geometry the epoch geometry.
template <template <Frame> class V>
    requires ib::FrameTagged<V>
void expect_rotational_round_trips(const ib::HelioGeometry& geometry) {
    const fx::vec3d sample{1.5, -2.5, 0.25};

    const V<Frame::HAE> hae{sample};
    expect_vector_near(ib::HEE2HAE(ib::HAE2HEE(hae, geometry), geometry).v, sample,
                       kRoundTripTolerance, "HAE->HEE->HAE");
    expect_vector_near(ib::HEEQ2HAE(ib::HAE2HEEQ(hae, geometry), geometry).v, sample,
                       kRoundTripTolerance, "HAE->HEEQ->HAE");

    const V<Frame::HEE> hee{sample};
    expect_vector_near(ib::HAE2HEE(ib::HEE2HAE(hee, geometry), geometry).v, sample,
                       kRoundTripTolerance, "HEE->HAE->HEE");

    const V<Frame::HEEQ> heeq{sample};
    expect_vector_near(ib::HAE2HEEQ(ib::HEEQ2HAE(heeq, geometry), geometry).v, sample,
                       kRoundTripTolerance, "HEEQ->HAE->HEEQ");

    // A rotation preserves length; that is the one property a transposed-index bug cannot fake.
    EXPECT_NEAR(fx::norm(sample), fx::norm(ib::HAE2HEEQ(hae, geometry).v), kRoundTripTolerance);
    EXPECT_NEAR(fx::norm(sample), fx::norm(ib::HAE2HEE(hae, geometry).v), kRoundTripTolerance);
}

TEST(IrbemCoordsHelio, RotationalRoundTripsForBothVectorKinds) {
    expect_rotational_round_trips<Position>(reference_geometry());
    expect_rotational_round_trips<FieldVector>(reference_geometry());
}

TEST(IrbemCoordsHelio, RotationMatricesAreOrthonormalAndRightHanded) {
    const ib::HelioGeometry& g = reference_geometry();
    for (const fx::mat3d& m : {g.hae_to_hee, g.hae_to_heeq}) {
        const fx::mat3d product = m * fx::transpose(m);
        const fx::mat3d identity = fx::mat3d::identity();
        for (std::size_t row = 0; row < 3; ++row) {
            for (std::size_t col = 0; col < 3; ++col) {
                EXPECT_NEAR(identity(row, col), product(row, col), 1e-15) << row << "," << col;
            }
        }
        // +1, not -1: a reflection would satisfy the orthogonality check above and still mirror the
        // magnetosphere.
        EXPECT_NEAR(1.0, fx::determinant(m), 1e-15);
    }
}

// ---- independent geometric checks, not derived from Table 8 -------------------------------------------

TEST(IrbemCoordsHelio, HeeXAxisPointsAtTheEarth) {
    const ib::HelioGeometry& g = reference_geometry();
    // +X_HEE is the Sun->Earth direction, so in HAE it must lie in the ecliptic plane at the
    // Earth's heliocentric longitude. That is the definition of λ_geo, checked end to end.
    const fx::vec3d axis = ib::HEE2HAE(FieldVector<Frame::HEE>{fx::vec3d{1.0, 0.0, 0.0}}, g).v;
    const double to_deg = 180.0 / std::acos(-1.0);
    EXPECT_NEAR(g.earth_longitude_deg, std::fmod(std::atan2(axis[1], axis[0]) * to_deg + 360.0,
                                                 360.0),
                1e-12);
    EXPECT_NEAR(0.0, axis[2], 1e-16);
}

TEST(IrbemCoordsHelio, HeeqZAxisIsTheSolarRotationAxis) {
    const ib::HelioGeometry& g = reference_geometry();
    // The pole of a plane whose ascending node is Ω and whose inclination is i sits at ecliptic
    // longitude Ω − 90 deg and latitude 90 − i. Recovering both from the transform is an
    // independent check of the E(Ω, i, θ) composition, and it does not involve θ at all.
    const fx::vec3d axis = ib::HEEQ2HAE(FieldVector<Frame::HEEQ>{fx::vec3d{0.0, 0.0, 1.0}}, g).v;
    const double to_deg = 180.0 / std::acos(-1.0);
    EXPECT_NEAR(std::fmod(g.solar_node_deg - 90.0 + 360.0, 360.0),
                std::fmod(std::atan2(axis[1], axis[0]) * to_deg + 360.0, 360.0), 1e-12);
    EXPECT_NEAR(90.0 - ib::solar_equator_inclination_deg(), std::asin(axis[2]) * to_deg, 1e-12);
}

TEST(IrbemCoordsHelio, HeeqXAxisLiesUnderTheEarth) {
    const ib::HelioGeometry& g = reference_geometry();
    // +X_HEEQ is where the solar central meridian crosses the solar equator: it is the projection
    // of the Sun->Earth direction onto the solar equatorial plane. So the two must share a
    // heliographic meridian — equivalently, the Earth direction has zero Y in HEEQ.
    //
    // Build the apparent Sun->Earth unit vector in HAE from the apparent longitude, then take it
    // into HEEQ; its Y component must vanish and its X must be positive.
    const double to_rad = std::acos(-1.0) / 180.0;
    const double lambda = g.apparent_earth_longitude_deg * to_rad;
    const FieldVector<Frame::HAE> to_earth{fx::vec3d{std::cos(lambda), std::sin(lambda), 0.0}};
    const fx::vec3d in_heeq = ib::HAE2HEEQ(to_earth, g).v;
    EXPECT_NEAR(0.0, in_heeq[1], 1e-15);
    EXPECT_GT(in_heeq[0], 0.0);
    // Out of the solar equatorial plane by the heliographic latitude of the Earth, which in late
    // August is close to its +7.25 deg extreme (the Earth crosses the solar equator in early June
    // and early December).
    EXPECT_NEAR(7.0, std::asin(in_heeq[2]) * (180.0 / std::acos(-1.0)), 0.5);
}
