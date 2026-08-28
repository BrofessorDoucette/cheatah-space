// Unit tests for space.irbem's IGRF internal field — the coefficient table, the time
// interpolation, and the branch-free spherical-harmonic kernel.
//
// The kernel is validated four independent ways, because "it agrees with itself" proves nothing:
//
//   1. Against IAGA's own published field values, fetched from the British Geological Survey's
//      IGRF-14 web service (BGS maintains the IGRF synthesis program for IAGA V-MOD). Ten points
//      spanning the equator, both poles, both hemispheres, sea level to 850 km, and 1900 to 2030.
//      BGS publishes components rounded to whole nT, which is the tolerance.
//   2. Against an independent textbook implementation written in this file: the *other*
//      formulation — geocentric colatitude, `sin`/`cos`, Schmidt Legendre functions carrying their
//      `sinᵐθ` factors, and an explicit `1/sin θ` in the eastward component. It shares no line and
//      no algebraic rearrangement with the header, so agreement to 1e-12 is real evidence that the
//      pole-free Cartesian derivation is correct rather than merely self-consistent.
//   3. Against the closed-form centred dipole at degree 1, where the whole series collapses to
//      `B = r⁻³[3(m·r̂)r̂ - m]` and can be checked in one line.
//   4. Against the coefficient table itself: at every tabulated epoch the interpolation must
//      return the published number bit-for-bit, so a transcription error in `tables/igrf14.hpp`
//      cannot hide behind a tolerance.
//
// Coverage counts each template INSTANTIATION separately, so every exercise that can be a template
// is one, and `drive_every_member` walks the whole public surface of each instantiation.
#include <gtest/gtest.h>

#include "alloc_counter.hpp"

#include <array>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <new>
#include <optional>
#include <type_traits>

#include "space/irbem/igrf.hpp"

namespace ib = cheatah::space::irbem;
namespace fx = cheatah::fixarray;
namespace tb = cheatah::space::irbem::tables;

using ib::Frame;

// ---- allocation counter ------------------------------------------------------------------------
// The kernel claims zero heap traffic; the only honest way to assert that is to count. The
// counter and the replaced global operators live in tests/alloc_counter.cpp — one definition for
// the whole binary, so nothing can allocate behind our back.

namespace {

constexpr double kDegToRad = 3.14159265358979323846 / 180.0;

// ---- an independent reference implementation ---------------------------------------------------
// Deliberately the textbook formulation, written from the same published definition but sharing no
// algebra with the header: colatitude and longitude as angles, Schmidt functions that keep their
// sinᵐθ, and the singular 1/sinθ that the header exists to avoid. Slow, allocating, and obviously
// correct — which is exactly what a reference is for. Valid away from the poles only.
struct ReferenceField {
    double bx = 0.0;
    double by = 0.0;
    double bz = 0.0;
};

template <int NMAX>
ReferenceField reference_field(const ib::Igrf<NMAX>& model, double x, double y, double z) {
    const double r = std::sqrt((x * x) + (y * y) + (z * z));
    const double theta = std::acos(z / r);  // geocentric colatitude
    const double phi = std::atan2(y, x);    // east longitude
    const double st = std::sin(theta);
    const double ct = std::cos(theta);

    // Schmidt semi-normalised Pⁿₘ(cos θ) and dPⁿₘ/dθ, carrying their sinᵐθ factors.
    std::array<std::array<double, NMAX + 1>, NMAX + 1> p{};
    std::array<std::array<double, NMAX + 1>, NMAX + 1> dp{};
    p[0][0] = 1.0;
    for (int n = 1; n <= NMAX; ++n) {
        const double k =
            (n == 1) ? 1.0
                     : std::sqrt(static_cast<double>((2 * n) - 1) / static_cast<double>(2 * n));
        p[n][n] = k * st * p[n - 1][n - 1];
        dp[n][n] = k * ((st * dp[n - 1][n - 1]) + (ct * p[n - 1][n - 1]));
        for (int m = 0; m < n; ++m) {
            const double nn = n;
            const double mm = m;
            const double den = std::sqrt((nn * nn) - (mm * mm));
            const double num = std::sqrt(((nn - 1.0) * (nn - 1.0)) - (mm * mm));
            // The `n - 2 >= m` guard with `m >= 0` makes `n - 2` non-negative wherever the
            // subscript is evaluated, but cppcheck cannot see through the ternary and reads the
            // index as potentially -1. Suppressed rather than restructured: rewriting correct,
            // readable recursion arithmetic to appease a checker is how a real bug gets hidden.
            // cppcheck-suppress negativeContainerIndex
            const double back = (n - 2 >= m) ? p[n - 2][m] : 0.0;
            // cppcheck-suppress negativeContainerIndex
            const double dback = (n - 2 >= m) ? dp[n - 2][m] : 0.0;
            p[n][m] = ((((2.0 * nn) - 1.0) * ct * p[n - 1][m]) - (num * back)) / den;
            dp[n][m] = (((((2.0 * nn) - 1.0) * ((ct * dp[n - 1][m]) - (st * p[n - 1][m]))) -
                         (num * dback))) /
                       den;
        }
    }

    double br = 0.0;
    double bt = 0.0;
    double bp = 0.0;
    for (int n = 1; n <= NMAX; ++n) {
        const double rn = std::pow(1.0 / r, n + 2);  // a/r with a = 1 Re by unit choice
        for (int m = 0; m <= n; ++m) {
            const double cm = std::cos(m * phi);
            const double sm = std::sin(m * phi);
            const double gg = model.g(n, m);
            const double hh = model.h(n, m);
            br += rn * (n + 1) * (((gg * cm) + (hh * sm)) * p[n][m]);
            bt -= rn * (((gg * cm) + (hh * sm)) * dp[n][m]);
            bp -= rn * (m * (((-gg * sm) + (hh * cm)) * p[n][m])) / st;
        }
    }
    const double cp = std::cos(phi);
    const double sp = std::sin(phi);
    return ReferenceField{(br * st * cp) + (bt * ct * cp) - (bp * sp),
                          (br * st * sp) + (bt * ct * sp) + (bp * cp),
                          (br * ct) - (bt * st)};
}

// ---- helpers -----------------------------------------------------------------------------------

ib::Position<Frame::GEO> geo(double x, double y, double z) {
    return ib::Position<Frame::GEO>{fx::vec3d{x, y, z}};
}

// WGS84 geodetic (latitude deg, east longitude deg, altitude km) to geocentric Cartesian in units
// of the IGRF reference radius. NIMA TR8350.2: a = 6378.137 km, 1/f = 298.257223563.
ib::Position<Frame::GEO> from_geodetic(double lat_deg, double lon_deg, double alt_km) {
    constexpr double kSemiMajorKm = 6378.137;
    constexpr double kFlattening = 1.0 / 298.257223563;
    constexpr double kEcc2 = kFlattening * (2.0 - kFlattening);
    const double sp = std::sin(lat_deg * kDegToRad);
    const double cp = std::cos(lat_deg * kDegToRad);
    const double sl = std::sin(lon_deg * kDegToRad);
    const double cl = std::cos(lon_deg * kDegToRad);
    const double prime = kSemiMajorKm / std::sqrt(1.0 - (kEcc2 * sp * sp));
    const double re = ib::Igrf<>::reference_radius_km;
    return geo(((prime + alt_km) * cp * cl) / re, ((prime + alt_km) * cp * sl) / re,
               (((prime * (1.0 - kEcc2)) + alt_km) * sp) / re);
}

// Geographic Cartesian field to the local geodetic north/east/down triad BGS reports in.
struct LocalField {
    double north = 0.0;
    double east = 0.0;
    double down = 0.0;
};

LocalField to_local(const ib::FieldVector<Frame::GEO>& b, double lat_deg, double lon_deg) {
    const double sp = std::sin(lat_deg * kDegToRad);
    const double cp = std::cos(lat_deg * kDegToRad);
    const double sl = std::sin(lon_deg * kDegToRad);
    const double cl = std::cos(lon_deg * kDegToRad);
    return LocalField{(-b.v[0] * sp * cl) - (b.v[1] * sp * sl) + (b.v[2] * cp),
                      (-b.v[0] * sl) + (b.v[1] * cl),
                      (-b.v[0] * cp * cl) - (b.v[1] * cp * sl) - (b.v[2] * sp)};
}

double relative_difference(double a, double b) {
    const double scale = std::fmax(std::fabs(a), std::fabs(b));
    return scale == 0.0 ? 0.0 : std::fabs(a - b) / scale;
}

}  // namespace

// ---- detail: the compile-time square root ------------------------------------------------------

// The normalisation table is only a win if it is genuinely built at compile time. These are
// constant expressions, so if `const_sqrt` ever stopped folding, this file would stop compiling.
static_assert(ib::detail::const_sqrt(4.0) == 2.0);
static_assert(ib::detail::const_sqrt(0.25) == 0.5);
static_assert(ib::detail::const_sqrt(0.0) == 0.0);
static_assert(ib::detail::const_sqrt(-1.0) == 0.0);
static_assert(ib::detail::triangular_slots(13) == tb::igrf14_slots);
static_assert(ib::detail::slot_index(13, 13) == tb::igrf14_slots - 1);
static_assert(ib::detail::slot_index(0, 0) == 0);

TEST(IgrfConstSqrt, MatchesStdSqrtOnEveryRadicandTheTableUses) {
    // Every value `make_legendre_normalisation` ever roots, at full degree.
    for (int n = 1; n <= tb::igrf14_max_degree; ++n) {
        const double ratio =
            static_cast<double>((2 * n) - 1) / static_cast<double>(2 * n);
        EXPECT_EQ(ib::detail::const_sqrt(ratio), std::sqrt(ratio)) << "diagonal n=" << n;
        for (int m = 0; m < n; ++m) {
            const double nn = n;
            const double mm = m;
            const double den = (nn * nn) - (mm * mm);
            const double num = ((nn - 1.0) * (nn - 1.0)) - (mm * mm);
            EXPECT_EQ(ib::detail::const_sqrt(den), std::sqrt(den)) << "n=" << n << " m=" << m;
            EXPECT_EQ(ib::detail::const_sqrt(num), std::sqrt(num)) << "n=" << n << " m=" << m;
        }
    }
}

TEST(IgrfConstSqrt, PerfectSquaresAreExactAndDegenerateInputsAreZero) {
    EXPECT_EQ(ib::detail::const_sqrt(1.0), 1.0);
    EXPECT_EQ(ib::detail::const_sqrt(9.0), 3.0);
    EXPECT_EQ(ib::detail::const_sqrt(1.0 / 16.0), 0.25);
    EXPECT_EQ(ib::detail::const_sqrt(0.0), 0.0);
    EXPECT_EQ(ib::detail::const_sqrt(-4.0), 0.0);
    EXPECT_EQ(ib::detail::const_sqrt(std::numeric_limits<double>::quiet_NaN()), 0.0);
    // Far from 1 in both directions, where a Newton seeded with the radicand itself would need
    // hundreds of iterations and a fixed budget would silently return the wrong answer.
    EXPECT_EQ(ib::detail::const_sqrt(1e-30), std::sqrt(1e-30));
    EXPECT_EQ(ib::detail::const_sqrt(1e30), std::sqrt(1e30));
    EXPECT_EQ(ib::detail::const_sqrt(1e-300), std::sqrt(1e-300));
    EXPECT_EQ(ib::detail::const_sqrt(1e300), std::sqrt(1e300));
}

TEST(IgrfConstSqrt, IsCorrectlyRoundedOverAWideSweep) {
    // The claim in the header is *correctly rounded*, not *close*. 40 000 rationals across five
    // decades, every one bit-identical to the hardware square root.
    std::size_t checked = 0;
    for (int i = 1; i <= 200; ++i) {
        for (int j = 1; j <= 200; ++j) {
            const double x = static_cast<double>(i) / static_cast<double>(j);
            ASSERT_EQ(ib::detail::const_sqrt(x), std::sqrt(x)) << "x = " << x;
            ++checked;
        }
    }
    EXPECT_EQ(checked, 40000U);
}

TEST(IgrfTriangularPacking, IndexesEverySlotExactlyOnce) {
    std::array<int, tb::igrf14_slots> seen{};
    for (int n = 0; n <= tb::igrf14_max_degree; ++n) {
        for (int m = 0; m <= n; ++m) {
            const std::size_t k = ib::detail::slot_index(n, m);
            ASSERT_LT(k, tb::igrf14_slots);
            ++seen[k];
        }
    }
    for (std::size_t k = 0; k < tb::igrf14_slots; ++k) EXPECT_EQ(seen[k], 1) << "slot " << k;
    EXPECT_EQ(ib::detail::triangular_slots(0), 1U);
    EXPECT_EQ(ib::detail::triangular_slots(1), 3U);
}

// The table is `static constexpr` inside Igrf, so it is never *called* at run time. Calling it
// here both checks the arithmetic against `std::sqrt` and gives the instantiation coverage.
template <int NMAX, class T>
void check_normalisation_table() {
    const auto table = ib::detail::make_legendre_normalisation<NMAX, T>();
    EXPECT_EQ(table.diagonal[0], static_cast<T>(1));
    double diag = 1.0;
    for (int n = 1; n <= NMAX; ++n) {
        diag *= (n == 1) ? 1.0
                         : std::sqrt(static_cast<double>((2 * n) - 1) / static_cast<double>(2 * n));
        EXPECT_EQ(table.diagonal[static_cast<std::size_t>(n)], static_cast<T>(diag))
            << "degree " << n;
        for (int m = 0; m < n; ++m) {
            const double nn = n;
            const double mm = m;
            const double den = std::sqrt((nn * nn) - (mm * mm));
            const double num = std::sqrt(((nn - 1.0) * (nn - 1.0)) - (mm * mm));
            const std::size_t k = ib::detail::slot_index(n, m);
            EXPECT_EQ(table.e[k], static_cast<T>(((2.0 * nn) - 1.0) / den));
            EXPECT_EQ(table.f[k], static_cast<T>(num / den));
            // The aliased degree-(n-2) read is only safe because this factor is EXACTLY zero.
            if (m == n - 1) {
                EXPECT_EQ(table.f[k], static_cast<T>(0));
            }
        }
    }
}

TEST(IgrfNormalisation, IsTheSchmidtRecursionAndVanishesOnTheAliasedDiagonal) {
    check_normalisation_table<13, double>();
    check_normalisation_table<10, double>();
    check_normalisation_table<1, double>();
    check_normalisation_table<13, float>();
}

// ---- the coefficient table and its interpolation -----------------------------------------------

TEST(IgrfTable, HasTheShapeIagaPublishes) {
    EXPECT_EQ(tb::igrf14_max_degree, 13);
    EXPECT_EQ(tb::igrf14_epoch_count, 26U);
    EXPECT_EQ(tb::igrf14_slots, 105U);
    for (std::size_t i = 0; i < tb::igrf14_epoch_count; ++i) {
        EXPECT_EQ(tb::igrf14_epochs[i], 1900.0 + (5.0 * static_cast<double>(i)));
        // Degree 0 is not part of a source-free potential, and hⁿ₀ multiplies sin(0φ) ≡ 0.
        EXPECT_EQ(tb::igrf14_g[i][0], 0.0);
        for (int n = 1; n <= tb::igrf14_max_degree; ++n) {
            EXPECT_EQ(tb::igrf14_h[i][ib::detail::slot_index(n, 0)], 0.0);
        }
    }
    // Spot values, read straight off the published file, as a transcription tripwire.
    EXPECT_EQ(tb::igrf14_g[0][ib::detail::slot_index(1, 0)], -31543.0);   // 1900 g₁⁰
    EXPECT_EQ(tb::igrf14_h[0][ib::detail::slot_index(1, 1)], 5922.0);     // 1900 h₁¹
    EXPECT_EQ(tb::igrf14_g[25][ib::detail::slot_index(1, 0)], -29350.0);  // 2025 g₁⁰
    EXPECT_EQ(tb::igrf14_h[25][ib::detail::slot_index(13, 13)], -0.5);    // 2025 h₁₃¹³
    EXPECT_EQ(tb::igrf14_g_sv[ib::detail::slot_index(1, 0)], 12.6);
    EXPECT_EQ(tb::igrf14_h_sv[ib::detail::slot_index(1, 1)], -21.5);
    // Secular variation is published only to degree 8.
    for (int n = 9; n <= tb::igrf14_max_degree; ++n) {
        for (int m = 0; m <= n; ++m) {
            EXPECT_EQ(tb::igrf14_g_sv[ib::detail::slot_index(n, m)], 0.0);
            EXPECT_EQ(tb::igrf14_h_sv[ib::detail::slot_index(n, m)], 0.0);
        }
    }
}

TEST(IgrfEpoch, RejectsDatesTheModelDoesNotCover) {
    EXPECT_FALSE(ib::Igrf<>::at(1899.999).has_value());
    EXPECT_FALSE(ib::Igrf<>::at(2030.001).has_value());
    EXPECT_FALSE(ib::Igrf<>::at(std::numeric_limits<double>::quiet_NaN()).has_value());
    EXPECT_FALSE(ib::Igrf<>::at(-std::numeric_limits<double>::infinity()).has_value());
    EXPECT_FALSE(ib::Igrf<>::at(std::numeric_limits<double>::infinity()).has_value());
    // The endpoints themselves are covered: 1900.0 is the first DGRF, 2030.0 the end of the
    // published secular-variation prediction. Past 2030 the model is REJECTED, not extrapolated.
    ASSERT_TRUE(ib::Igrf<>::at(1900.0).has_value());
    ASSERT_TRUE(ib::Igrf<>::at(2030.0).has_value());
    EXPECT_EQ(ib::Igrf<>::at(1900.0)->year(), 1900.0);
    EXPECT_EQ(ib::Igrf<>::at(2030.0)->year(), 2030.0);
}

TEST(IgrfEpoch, ReproducesEveryTabulatedEpochBitForBit) {
    for (std::size_t i = 0; i < tb::igrf14_epoch_count; ++i) {
        const auto model = ib::Igrf<>::at(tb::igrf14_epochs[i]);
        ASSERT_TRUE(model.has_value()) << tb::igrf14_epochs[i];
        for (int n = 1; n <= tb::igrf14_max_degree; ++n) {
            for (int m = 0; m <= n; ++m) {
                const std::size_t k = ib::detail::slot_index(n, m);
                EXPECT_EQ(model->g(n, m), tb::igrf14_g[i][k])
                    << "epoch " << tb::igrf14_epochs[i] << " g " << n << "," << m;
                EXPECT_EQ(model->h(n, m), tb::igrf14_h[i][k])
                    << "epoch " << tb::igrf14_epochs[i] << " h " << n << "," << m;
            }
        }
    }
}

TEST(IgrfEpoch, InterpolatesLinearlyBetweenEpochs) {
    // 2002.5 is the exact midpoint of the 2000/2005 interval, and the weight 0.5 is exactly
    // representable, so the answer is the exact arithmetic mean and the assertion can be `==`.
    const auto mid = ib::Igrf<>::at(2002.5);
    ASSERT_TRUE(mid.has_value());
    for (int n = 1; n <= tb::igrf14_max_degree; ++n) {
        for (int m = 0; m <= n; ++m) {
            const std::size_t k = ib::detail::slot_index(n, m);
            const double lo = tb::igrf14_g[20][k];
            const double hi = tb::igrf14_g[21][k];
            EXPECT_EQ(mid->g(n, m), lo + (0.5 * (hi - lo)));
            const double lo_h = tb::igrf14_h[20][k];
            const double hi_h = tb::igrf14_h[21][k];
            EXPECT_EQ(mid->h(n, m), lo_h + (0.5 * (hi_h - lo_h)));
        }
    }
    // A quarter of the way into the 1955/1960 interval, likewise exact.
    const auto q = ib::Igrf<>::at(1956.25);
    ASSERT_TRUE(q.has_value());
    const std::size_t k10 = ib::detail::slot_index(1, 0);
    EXPECT_EQ(q->g(1, 0), tb::igrf14_g[11][k10] +
                              (0.25 * (tb::igrf14_g[12][k10] - tb::igrf14_g[11][k10])));
}

TEST(IgrfEpoch, ExtrapolatesPastTheLastEpochWithPublishedSecularVariation) {
    const auto model = ib::Igrf<>::at(2028.0);
    ASSERT_TRUE(model.has_value());
    for (int n = 1; n <= tb::igrf14_max_degree; ++n) {
        for (int m = 0; m <= n; ++m) {
            const std::size_t k = ib::detail::slot_index(n, m);
            EXPECT_EQ(model->g(n, m), tb::igrf14_g[25][k] + (3.0 * tb::igrf14_g_sv[k]));
            EXPECT_EQ(model->h(n, m), tb::igrf14_h[25][k] + (3.0 * tb::igrf14_h_sv[k]));
        }
    }
}

TEST(IgrfEpoch, IsAConstantExpression) {
    // The headline claim: a model for a known epoch costs nothing at run time.
    constexpr auto model = ib::Igrf<>::at(2025.0);
    static_assert(model.has_value());
    static_assert(model->degree == 13);
    static_assert(model->g(1, 0) == -29350.0);
    static_assert(model->h(1, 1) == 4545.5);
    static_assert(model->year() == 2025.0);
    static_assert(!ib::Igrf<>::at(1850.0).has_value());
    // Interpolation folds too.
    constexpr auto mid = ib::Igrf<>::at(2002.5);
    static_assert(mid.has_value());
    static_assert(mid->g(1, 0) == -29619.4 + (0.5 * (-29554.63 - -29619.4)));
    EXPECT_EQ(model->g(1, 0), -29350.0);  // and a runtime read of the same object
}

TEST(IgrfCoefficients, ReportExactZeroOutsideTheTruncation) {
    const auto full = ib::Igrf<13>::at(2025.0);
    const auto cut = ib::Igrf<10>::at(2025.0);
    ASSERT_TRUE(full.has_value());
    ASSERT_TRUE(cut.has_value());
    EXPECT_EQ(full->g(0, 0), 0.0);
    EXPECT_EQ(full->h(0, 0), 0.0);
    EXPECT_EQ(full->g(14, 0), 0.0);
    EXPECT_EQ(full->h(14, 0), 0.0);
    EXPECT_EQ(full->g(-1, 0), 0.0);
    EXPECT_EQ(full->h(-1, 0), 0.0);
    EXPECT_EQ(full->g(3, 4), 0.0);   // order above degree
    EXPECT_EQ(full->h(3, -1), 0.0);  // negative order
    EXPECT_EQ(full->h(5, 0), 0.0);   // hⁿ₀ is zero by definition, not by truncation
    // A truncated model agrees exactly where it is defined and is zero above.
    for (int n = 1; n <= 10; ++n) {
        for (int m = 0; m <= n; ++m) {
            EXPECT_EQ(cut->g(n, m), full->g(n, m));
            EXPECT_EQ(cut->h(n, m), full->h(n, m));
        }
    }
    EXPECT_EQ(cut->g(11, 0), 0.0);
    EXPECT_NE(full->g(11, 0), 0.0);
}

// ---- the field kernel --------------------------------------------------------------------------

TEST(IgrfField, DegreeOneIsTheClosedFormCentredDipole) {
    // At NMAX = 1 the entire series is B = r⁻³[3(m·r̂)r̂ - m] with m = (g₁¹, h₁¹, g₁⁰). If the
    // Cartesian rearrangement in the header has an algebra error this is where it shows.
    for (const double year : {1900.0, 1975.0, 2025.0, 2030.0}) {
        const auto model = ib::Igrf<1>::at(year);
        ASSERT_TRUE(model.has_value());
        const fx::vec3d moment{model->g(1, 1), model->h(1, 1), model->g(1, 0)};
        for (const auto& p : {geo(2.0, 0.0, 0.0), geo(0.0, 0.0, 4.0), geo(0.0, -2.5, 0.0),
                              geo(3.0, 4.0, 12.0), geo(-1.5, 0.5, -0.25)}) {
            const double r = fx::norm(p.v);
            const fx::vec3d rhat = p.v * (1.0 / r);
            const double proj = fx::dot(moment, rhat);
            const fx::vec3d want = ((rhat * (3.0 * proj)) - moment) * (1.0 / (r * r * r));
            const auto got = model->evaluate(p);
            for (std::size_t i = 0; i < 3; ++i) {
                EXPECT_LT(relative_difference(got.v[i], want[i]), 1e-14)
                    << "year " << year << " component " << i;
            }
        }
    }
}

// Away from the poles the header must agree with the textbook formulation to roundoff. 1e-12 is
// ~4500 ulp of the largest partial sum, which is the honest cost of two different summation orders
// over ~200 terms whose partials reach 1e5 nT; the agreement observed is nearer 1e-14.
template <int NMAX>
void check_against_reference(double year) {
    const auto model = ib::Igrf<NMAX>::at(year);
    ASSERT_TRUE(model.has_value());
    for (const auto& p : {geo(1.0, 0.0, 0.0), geo(0.0, 1.0, 0.0), geo(0.5, 0.5, 0.5),
                          geo(-0.75, 0.25, 1.5), geo(3.0, 4.0, 12.0), geo(-6.0, -2.0, 0.5),
                          geo(0.25, -0.125, -1.0), geo(-0.5, 0.0, -2.0)}) {
        const auto want = reference_field<NMAX>(*model, p.v[0], p.v[1], p.v[2]);
        const auto got = model->evaluate(p);
        EXPECT_LT(relative_difference(got.v[0], want.bx), 1e-12) << "Bx at year " << year;
        EXPECT_LT(relative_difference(got.v[1], want.by), 1e-12) << "By at year " << year;
        EXPECT_LT(relative_difference(got.v[2], want.bz), 1e-12) << "Bz at year " << year;
    }
}

TEST(IgrfField, MatchesAnIndependentTextbookImplementation) {
    for (const double year : {1900.0, 1937.5, 1980.0, 2002.5, 2020.0, 2025.0, 2029.75}) {
        check_against_reference<13>(year);
        check_against_reference<10>(year);
        check_against_reference<1>(year);
    }
}

TEST(IgrfField, ThePoleIsAnOrdinaryPoint) {
    const auto model = ib::Igrf<>::at(2025.0);
    ASSERT_TRUE(model.has_value());
    for (const double sign : {1.0, -1.0}) {
        const auto at_pole = model->evaluate(geo(0.0, 0.0, sign * 1.5));
        // Finite, and the reference implementation's limit as it is approached along two different
        // meridians must converge to the same vector — that is the singularity actually cancelling
        // rather than being papered over.
        for (std::size_t i = 0; i < 3; ++i) EXPECT_TRUE(std::isfinite(at_pole.v[i]));
        // 1e-2 nT is 4e-7 of the field — set by the REFERENCE, whose 1/sinθ multiplies its own
        // roundoff by 1e7 this close to the axis, not by the kernel, which has no such factor.
        for (const double azimuth : {0.0, 90.0, 217.0}) {
            const double eps = 1e-7;
            const double a = azimuth * kDegToRad;
            const auto near = reference_field<13>(*model, eps * std::cos(a), eps * std::sin(a),
                                                  sign * 1.5);
            EXPECT_NEAR(at_pole.v[0], near.bx, 1e-2) << "azimuth " << azimuth;
            EXPECT_NEAR(at_pole.v[1], near.by, 1e-2) << "azimuth " << azimuth;
            EXPECT_NEAR(at_pole.v[2], near.bz, 1e-2) << "azimuth " << azimuth;
        }
        // And the kernel itself is continuous into the pole, to its own roundoff.
        const auto just_off = model->evaluate(geo(1e-9, 1e-9, sign * 1.5));
        for (std::size_t i = 0; i < 3; ++i) EXPECT_NEAR(at_pole.v[i], just_off.v[i], 1e-3);
    }
}

TEST(IgrfField, MatchesTheValuesIagaPublishes) {
    // Reference values from the British Geological Survey's IGRF-14 web service — BGS maintains
    // the IGRF synthesis program for IAGA V-MOD, so these are the model's own published output:
    //   https://geomag.bgs.ac.uk/web_service/GMModels/igrf/14/?latitude=..&longitude=..&altitude=..&date=..
    // retrieved 2026-08-28. Positions are WGS84 geodetic, altitude in km above the ellipsoid,
    // components in the local north/east/down triad. Dates are 1 January so the decimal year is
    // exactly the integer and no day-of-year convention can creep in.
    //
    // Tolerance: BGS rounds every component to whole nT, so ±1 nT is the resolution of the
    // reference, not of this kernel. The largest observed deviation is 0.57 nT — pure rounding.
    struct Case {
        double lat;
        double lon;
        double alt_km;
        double year;
        double north;
        double east;
        double down;
    };
    constexpr std::array<Case, 10> cases{{
        {0.0, 0.0, 0.0, 2025.0, 27457, -1926, -15997},
        {0.0, 120.0, 500.0, 2020.0, 30880, 76, -8843},
        {90.0, 0.0, 0.0, 2010.0, 1888, -461, 56568},
        {-90.0, 0.0, 0.0, 1950.0, 14085, -6936, -58453},
        {45.0, -75.0, 0.0, 2027.0, 18474, -4196, 49447},
        {-30.5, 150.25, 200.0, 2003.0, 23835, 4657, -44152},
        {60.0, 45.0, 100.0, 2027.0, 13307, 3705, 50716},
        {0.0, 0.0, 0.0, 1900.0, 28028, -8560, -5590},
        {25.0, -100.0, 0.0, 2030.0, 25456, 1746, 33932},
        {-60.25, 30.5, 850.0, 1975.0, 9432, -6930, -25732},
    }};
    for (const auto& c : cases) {
        const auto model = ib::Igrf<>::at(c.year);
        ASSERT_TRUE(model.has_value());
        const auto local =
            to_local(model->evaluate(from_geodetic(c.lat, c.lon, c.alt_km)), c.lat, c.lon);
        EXPECT_NEAR(local.north, c.north, 1.0) << "X at " << c.lat << "," << c.lon << " " << c.year;
        EXPECT_NEAR(local.east, c.east, 1.0) << "Y at " << c.lat << "," << c.lon << " " << c.year;
        EXPECT_NEAR(local.down, c.down, 1.0) << "Z at " << c.lat << "," << c.lon << " " << c.year;
        // The total intensity is the quantity the reference quotes most precisely.
        const double f = std::sqrt((local.north * local.north) + (local.east * local.east) +
                                   (local.down * local.down));
        const double want_f = std::sqrt((c.north * c.north) + (c.east * c.east) +
                                        (c.down * c.down));
        EXPECT_NEAR(f, want_f, 1.0);
    }
}

TEST(IgrfField, TruncationDegreeIsAPhysicalDifferenceNotARoundingOne) {
    // IRBEM truncates the internal field at degree 10; IGRF-14 publishes 13. Near the surface the
    // omitted degrees are worth tens of nT — far above any tolerance in this file — so a
    // differential test against IRBEM must say which truncation it ran.
    const auto full = ib::Igrf<13>::at(2025.0);
    const auto cut = ib::Igrf<10>::at(2025.0);
    ASSERT_TRUE(full.has_value());
    ASSERT_TRUE(cut.has_value());
    const auto p = geo(0.6, -0.5, 0.62449979983983983);  // |r| ≈ 1 Re, mid-latitude
    const auto b13 = full->evaluate(p);
    const auto b10 = cut->evaluate(p);
    const double gap = fx::norm((b13 - b10).v);
    EXPECT_GT(gap, 1.0) << "degrees 11-13 must matter at the surface";
    EXPECT_LT(gap, 200.0) << "…but they are a correction, not the field";
    // High above, the omitted degrees fall off as r⁻¹⁵ and the two agree.
    const auto far = geo(0.0, 0.0, 8.0);
    EXPECT_LT(fx::norm((full->evaluate(far) - cut->evaluate(far)).v), 1e-9);
}

TEST(IgrfField, TheSphericalEntryPointAgreesWithTheCartesianOne) {
    const auto model = ib::Igrf<>::at(2020.0);
    ASSERT_TRUE(model.has_value());
    struct Sph {
        double r;
        double lat;
        double lon;
    };
    for (const auto& s : {Sph{1.0, 0.0, 0.0}, Sph{1.5, 45.0, 90.0}, Sph{2.0, -30.0, -120.0},
                          Sph{6.6, 12.5, 210.0}, Sph{1.25, 89.0, 45.0}}) {
        const double cl = std::cos(s.lat * kDegToRad);
        const auto cart = geo(s.r * cl * std::cos(s.lon * kDegToRad),
                              s.r * cl * std::sin(s.lon * kDegToRad),
                              s.r * std::sin(s.lat * kDegToRad));
        const auto a = model->evaluate(ib::Position<Frame::SPH>{fx::vec3d{s.r, s.lat, s.lon}});
        const auto b = model->evaluate(cart);
        for (std::size_t i = 0; i < 3; ++i) EXPECT_EQ(a.v[i], b.v[i]) << "component " << i;
    }
}

TEST(IgrfField, TheFastPrecisionPolicyTracksTheExactOne) {
    // The precision seam, wired to policy.hpp: `Fast` narrows only the recursions
    // (`integrand = float`) while `SoundPrecision` holds the reduction at `accum = double`. The
    // error budget allows ~1e-6 relative on B; what we actually see is a few times 1e-7.
    static_assert(std::is_same_v<ib::Igrf<13, ib::Fast>::integrand, float>);
    static_assert(std::is_same_v<ib::Igrf<13, ib::Fast>::accum, double>);
    static_assert(std::is_same_v<ib::Igrf<13, ib::Exact>::integrand, double>);
    using Dbl = ib::Igrf<13, ib::Exact>;
    using Flt = ib::Igrf<13, ib::Fast>;
    const auto d = Dbl::at(2025.0);
    const auto f = Flt::at(2025.0);
    ASSERT_TRUE(d.has_value());
    ASSERT_TRUE(f.has_value());
    EXPECT_EQ(f->year(), 2025.0);
    EXPECT_EQ(f->g(1, 0), -29350.0F);
    EXPECT_EQ(f->h(1, 1), 4545.5F);
    EXPECT_EQ(f->g(0, 0), 0.0F);
    EXPECT_EQ(f->h(0, 0), 0.0F);
    for (const auto& p : {geo(1.0, 0.0, 0.0), geo(0.5, 0.5, 0.75), geo(-2.0, 3.0, 1.0),
                          geo(0.0, 0.0, 1.0)}) {
        const auto bd = d->evaluate(p);
        const auto bf = f->evaluate(p);
        for (std::size_t i = 0; i < 3; ++i) {
            EXPECT_LT(std::fabs(bd.v[i] - bf.v[i]) / bd.magnitude(), 1e-5) << "component " << i;
        }
    }
    const auto sf = f->evaluate(ib::Position<Frame::SPH>{fx::vec3d{1.5, 30.0, 60.0}});
    const auto sd = d->evaluate(ib::Position<Frame::SPH>{fx::vec3d{1.5, 30.0, 60.0}});
    EXPECT_LT(std::fabs(sf.magnitude() - sd.magnitude()) / sd.magnitude(), 1e-5);
}

// Every public member of one instantiation, so no instantiation is left with a function nobody
// called. Returns a checksum so nothing can be optimised away.
template <int NMAX, ib::SoundPrecision P>
double drive_every_member() {
    using Model = ib::Igrf<NMAX, P>;
    EXPECT_FALSE(Model::at(1899.0).has_value());   // below the first DGRF
    EXPECT_FALSE(Model::at(2031.0).has_value());   // past the secular-variation prediction
    const auto model = Model::at(2011.25);
    EXPECT_TRUE(model.has_value());
    EXPECT_EQ(model->year(), 2011.25);
    EXPECT_EQ(Model::degree, NMAX);
    EXPECT_EQ(Model::reference_radius_km, 6371.2);
    EXPECT_EQ(Model::earliest_year, 1900.0);
    EXPECT_EQ(Model::latest_epoch_year, 2025.0);
    EXPECT_EQ(Model::latest_year, 2030.0);
    double sum = static_cast<double>(model->g(1, 0)) + static_cast<double>(model->h(1, 1)) +
                 static_cast<double>(model->g(99, 0)) + static_cast<double>(model->h(99, 0));
    sum += model->evaluate(geo(1.0, 2.0, 2.0)).magnitude();
    sum += model->evaluate(ib::Position<Frame::SPH>{fx::vec3d{1.75, -20.0, 33.0}}).magnitude();
    return sum;
}

TEST(IgrfInstantiations, EveryOneIsFullyExercised) {
    const double a = drive_every_member<13, ib::Exact>();
    const double b = drive_every_member<10, ib::Exact>();
    const double c = drive_every_member<1, ib::Exact>();
    const double d = drive_every_member<13, ib::Fast>();
    for (const double v : {a, b, c, d}) EXPECT_TRUE(std::isfinite(v));
    EXPECT_NE(a, b);
    EXPECT_NE(a, c);
    EXPECT_LT(relative_difference(a, d), 1e-5);
}

TEST(IgrfField, AllocatesNothing) {
    const auto model = ib::Igrf<>::at(2025.0);
    ASSERT_TRUE(model.has_value());
    const auto p = geo(1.5, -2.0, 3.25);
    volatile double sink = 0.0;
    // Warm everything up first, then count.
    sink += model->evaluate(p).magnitude();
    const std::size_t before = cheatah_space_test::allocation_count();
    for (int i = 0; i < 1000; ++i) sink += model->evaluate(p).v[0];
    sink += model->evaluate(ib::Position<Frame::SPH>{fx::vec3d{2.0, 10.0, 20.0}}).v[1];
    const std::size_t after = cheatah_space_test::allocation_count();
    EXPECT_EQ(after, before) << "the field kernel must not touch the heap";
    // And building a model is allocation-free too.
    const std::size_t before_build = cheatah_space_test::allocation_count();
    const auto other = ib::Igrf<>::at(1963.75);
    EXPECT_EQ(cheatah_space_test::allocation_count(), before_build);
    ASSERT_TRUE(other.has_value());
    EXPECT_TRUE(std::isfinite(static_cast<double>(sink)));
}

TEST(IgrfField, Throughput) {
    // A CPU baseline for the orchestrator, not a pass/fail gate. The checksum is accumulated and
    // asserted so no part of the loop can be elided.
    const auto model = ib::Igrf<>::at(2025.0);
    ASSERT_TRUE(model.has_value());
    constexpr int kPoints = 4096;
    constexpr int kRepeats = 64;
    std::array<ib::Position<Frame::GEO>, kPoints> points{};
    for (int i = 0; i < kPoints; ++i) {
        const double t = static_cast<double>(i) / kPoints;
        points[static_cast<std::size_t>(i)] =
            geo(1.0 + (5.0 * t), -3.0 + (7.0 * t), 0.5 - (4.0 * t));
    }
    double checksum = 0.0;
    const auto start = std::chrono::steady_clock::now();
    for (int rep = 0; rep < kRepeats; ++rep) {
        for (const auto& p : points) checksum += model->evaluate(p).v[2];
    }
    const auto elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - start);
    const double evals = static_cast<double>(kPoints) * kRepeats;
    std::printf("[   PERF   ] IGRF-14 degree 13: %.0f evaluations in %.3f s = %.2f M eval/s "
                "(%.1f ns each)\n",
                evals, elapsed.count(), evals / elapsed.count() / 1e6,
                elapsed.count() / evals * 1e9);
    EXPECT_TRUE(std::isfinite(checksum));
    EXPECT_NE(checksum, 0.0);
    EXPECT_GT(elapsed.count(), 0.0);
}
