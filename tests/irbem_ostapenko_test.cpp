/// @file irbem_ostapenko_test.cpp
/// @brief The suite for `space/irbem/ext_ostapenko.hpp` — Ostapenko & Maltsev (1997).
///
/// A regression model has no analytic truth to compare against, so what is checked here is every
/// question about whether the implementation IS the paper:
///
///  - **Transcription, two independent ways.** The Cartesian harmonics in the header are compared
///    to the paper's Table 1 as PRINTED (cylindrical, with `cos phi` and `sin phi`) at random
///    off-axis points, and the ten curl-free rows are compared to an independent construction of
///    the Schmidt-normalized solid harmonics the paper's text says they are — which is the check
///    that pins the normalization of row 17 the table misprints.
///  - **`div B = 0`**, per harmonic and for the whole field: a second-order stencil residual must
///    fall as `h^2`. Model-independent, oracle-independent, and impossible to satisfy with a wrong
///    sign in any row.
///  - **Continuity in every driver.** The model is a linear regression, so nearby drivers give
///    nearby fields, second differences vanish, and — the assertion a Kp-BINNED model like T89
///    would fail — nearby Kp values give DIFFERENT fields. Southward Bz is swept densely.
///  - **The physics the paper reports**: its Figure 5 profiles — the equatorial `B_z` depression
///    deepens with negative Dst, the dayside compresses with pressure.
///  - **Validity from both sides** of every bound the paper states, and the batch fold reporting
///    exactly what the scalar lane does.
///  - **The differential against the IRBEM oracle**, when present: the functional form to roundoff
///    by regression, and the field itself to the rounding of the printed tables, per corpus regime.
///
/// The oracle test `dlopen`s IRBEM at runtime rather than linking it: this binary must build and
/// pass on a machine that has never heard of IRBEM.

#include <gtest/gtest.h>

#include <dlfcn.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <limits>
#include <numbers>
#include <string>
#include <vector>

#include "alloc_counter.hpp"
#include "irbem_domain_corpus.hpp"
#include "space/irbem/api.hpp"
#include "space/irbem/ext_ostapenko.hpp"
#include "space/irbem/lstar.hpp"

namespace {

namespace ir = cheatah::space::irbem;
namespace corpus = cheatah_space_test;

using ir::Frame;
using ir::om97_amplitudes;
using ir::om97_basis;
using ir::om97_check_drivers;
using ir::om97_check_fitted_region;
using ir::om97_components;
using ir::om97_field;
using ir::om97_field_at;
using ir::om97_field_batch;
using ir::om97_field_host;
using ir::om97_fitted_region;
using ir::om97_harmonic_count;
using ir::om97_normalization_measured;
using ir::om97_normalization_published;
using ir::om97_param_block;
using ir::om97_param_count;
using ir::om97_relation_coefficients;
using ir::Om97Basis;
using ir::Om97Drivers;
using ir::Om97Normalization;
using ir::Position;
using ir::Status;

/// A sink the optimizer cannot see through, so the allocation test's calls actually happen.
volatile double sink = 0.0;

/// A GSM point, spelled so a test reads as coordinates rather than as a constructor call.
Position<Frame::GSM> at(double x, double y, double z) {
    return Position<Frame::GSM>{cheatah::fixarray::vec3d(x, y, z)};
}

/// The tilts the suite sweeps: zero, and both signs of a realistic seasonal-diurnal excursion.
constexpr std::array<double, 5> kTilts{-0.55, -0.21, 0.0, 0.21, 0.55};

/// The corpus's moderate drivers, the working driver set for tests that need just one.
Om97Drivers moderate() {
    const corpus::MagInput& m = corpus::regime_drivers[1];
    return Om97Drivers{m.dst, m.pdyn, m.kp * 10.0, m.bz_imf};
}

/// A corpus regime's drivers as the model reads them.
Om97Drivers drivers_of(const corpus::MagInput& m) {
    return Om97Drivers{m.dst, m.pdyn, m.kp * 10.0, m.bz_imf};
}

/// The drivers sitting exactly on the published means, where every normalized parameter is zero.
Om97Drivers at_means() {
    const Om97Normalization& n = om97_normalization_published;
    return Om97Drivers{n.dst_mean, n.pdyn_mean, n.kp_mean * 10.0, n.bz_mean};
}

/// `div B` at one point by central differences.
double divergence(const std::array<double, om97_harmonic_count>& amp, double ps, double x,
                  double y, double z, double h) {
    const double s = std::sin(ps);
    const double c = std::cos(ps);
    const auto b = [&](double a, double d, double e) {
        return om97_components<double>(amp, s, c, a, d, e);
    };
    return ((b(x + h, y, z)[0] - b(x - h, y, z)[0]) + (b(x, y + h, z)[1] - b(x, y - h, z)[1]) +
            (b(x, y, z + h)[2] - b(x, y, z - h)[2])) /
           (2.0 * h);
}

/// The largest `|div B|` over the sampled box at step @p h, over every tilt and corpus regime.
double worst_divergence(double h) {
    double worst = 0.0;
    for (double ps : kTilts) {
        for (const corpus::MagInput& m : corpus::regime_drivers) {
            const std::array<double, om97_harmonic_count> amp =
                om97_amplitudes<double>(drivers_of(m));
            // Integer induction, coordinates derived per iteration (cert-flp30).
            for (int ix = 0; ix <= 8; ++ix)
                for (int iy = 0; iy <= 6; ++iy)
                    for (int iz = 0; iz <= 5; ++iz) {
                        const double x = -9.5 + (2.3 * ix);
                        const double y = -8.0 + (2.7 * iy);
                        const double z = -6.0 + (2.4 * iz);
                        worst = std::max(worst, std::fabs(divergence(amp, ps, x, y, z, h)));
                    }
        }
    }
    return worst;
}

/// A deterministic scatter of GSM points inside the paper's fitted region at every tilt the batch
/// tests use (|psi| <= 0.35): r <= 8.5 and |z| <= 3.5 keeps `rho_SM <= 8.5` and
/// `|z_SM| <= 8.5 sin 0.35 + 3.5 cos 0.35 = 6.2 < 7`. A 64-bit LCG so a disagreement is
/// reproducible, and never a lattice, so no component is systematically zero.
std::vector<Position<Frame::GSM>> scatter(std::size_t n) {
    std::vector<Position<Frame::GSM>> out;
    out.reserve(n);
    std::uint64_t s = 0x9E3779B97F4A7C15ULL;
    const auto next = [&s] {
        s = (s * 6364136223846793005ULL) + 1442695040888963407ULL;
        return static_cast<double>(s >> 11) / 9007199254740992.0;
    };
    for (std::size_t i = 0; i < n; ++i) {
        const double r = 3.0 + (5.5 * next());
        const double th = std::acos(1.0 - (2.0 * next()));
        const double ph = 6.283185307179586 * next();
        const double z = r * std::cos(th);
        out.push_back(at(r * std::sin(th) * std::cos(ph), r * std::sin(th) * std::sin(ph),
                         std::fabs(z) > 3.5 ? 3.5 * (z < 0.0 ? -1.0 : 1.0) : z));
    }
    return out;
}

// ---- the paper's Table 1 as PRINTED, in cylindrical coordinates --------------------------------

/// One printed row, `(b_rho, b_phi, b_z)`, at cylindrical `(rho, phi, z)` and tilt sine @p sp.
/// Transcribed from the table WITHOUT the Cartesian rewrite, and with row 17's Schmidt factor.
std::array<double, 3> printed_row(std::size_t i, double rho, double phi, double z, double sp) {
    const double c = std::cos(phi);
    const double s = std::sin(phi);
    const double r2 = rho * rho;
    const double r4 = r2 * r2;
    const double z2 = z * z;
    const double s3 = std::numbers::sqrt3;
    const double s10 = std::sqrt(10.0);
    const double s32 = std::sqrt(1.5);
    const double s15 = std::sqrt(15.0);
    switch (i) {
        case 1: return {0.0, 0.0, 1.0};
        case 2: return {-3.0 * rho * z, 0.0, 3.0 * (z2 - (r2 / 2.0))};
        case 3: return {-5.0 * rho * z * ((2.0 * z2) - (3.0 * r2 / 2.0)), 0.0,
                        5.0 * ((z2 * z2) - (3.0 * r2 * z2) + (3.0 * r4 / 8.0))};
        case 4: return {0.0, 0.0, r2};
        case 5: return {0.0, 0.0, r4};
        case 6: return {-2.0 * rho * z2 * z, 0.0, z2 * z2};
        case 7: return {s3 * z * c, -s3 * z * s, s3 * rho * c};
        case 8: return {s10 * z * (z2 - (9.0 * r2 / 4.0)) * c, -s10 * z * (z2 - (3.0 * r2 / 4.0)) * s,
                        3.0 * s10 * rho * (z2 - (r2 / 4.0)) * c};
        case 9: return {z * c, -z * s, 0.0};
        case 10: return {0.0, 0.0, r2 * rho * c};
        case 11: return {z2 * z * c, -z2 * z * s, 0.0};
        case 12: return {0.0, r2 * z * s, -0.5 * rho * z2 * c};
        case 13: return {-rho * sp, 0.0, 2.0 * z * sp};
        case 14: return {3.0 * rho * ((r2 / 2.0) - (2.0 * z2)) * sp, 0.0,
                         2.0 * z * ((2.0 * z2) - (3.0 * r2)) * sp};
        case 15: return {c * sp, -s * sp, 0.0};
        case 16: return {s32 * ((-3.0 * r2 / 2.0) + (2.0 * z2)) * c * sp,
                         s32 * ((r2 / 2.0) - (2.0 * z2)) * s * sp, 4.0 * s32 * rho * z * c * sp};
        // Row 17 with the Schmidt factor sqrt(15)/8 in place of the printed 1/(8 sqrt(15)); the
        // polynomials are the table's.
        case 17: return {(s15 / 8.0) * ((5.0 * r4) - (36.0 * r2 * z2) + (8.0 * z2 * z2)) * c * sp,
                         -(s15 / 8.0) * (r4 - (12.0 * r2 * z2) + (8.0 * z2 * z2)) * s * sp,
                         s15 * ((4.0 * z2) - (3.0 * r2)) * rho * z * c * sp};
        default: return {0.0, 0.0, 0.0};
    }
}

// ---- an independent construction of the Schmidt-normalized solid harmonics ---------------------

/// The unnormalized associated Legendre function `P_n^m(u)` without the Condon-Shortley phase,
/// by the standard recurrences — a code path that shares nothing with the header.
double legendre(int n, int m, double u) {
    double pmm = 1.0;
    if (m > 0) {
        const double s = std::sqrt(1.0 - (u * u));
        double f = 1.0;
        for (int k = 1; k <= m; ++k) {
            pmm *= f * s;
            f += 2.0;
        }
    }
    if (n == m) return pmm;
    double pmm1 = u * (2.0 * m + 1.0) * pmm;
    if (n == m + 1) return pmm1;
    double p = 0.0;
    for (int k = m + 2; k <= n; ++k) {
        p = ((u * (2.0 * k - 1.0) * pmm1) - ((k + m - 1.0) * pmm)) / (k - m);
        pmm = pmm1;
        pmm1 = p;
    }
    return p;
}

/// `r^n S_n^m(cos theta) cos(m phi)`: the Schmidt quasi-normalized solid harmonic potential.
double schmidt_solid(int n, int m, double x, double y, double z) {
    const double r = std::sqrt((x * x) + (y * y) + (z * z));
    if (r == 0.0) return 0.0;
    double norm = 1.0;
    if (m > 0) {
        double ratio = 1.0;  // (n-m)! / (n+m)!
        for (int k = n - m + 1; k <= n + m; ++k) ratio /= k;
        norm = std::sqrt(2.0 * ratio);
    }
    return std::pow(r, n) * norm * legendre(n, m, z / r) * std::cos(m * std::atan2(y, x));
}

/// The gradient of @ref schmidt_solid by central differences.
std::array<double, 3> schmidt_gradient(int n, int m, double x, double y, double z) {
    const double h = 1e-5;
    return {(schmidt_solid(n, m, x + h, y, z) - schmidt_solid(n, m, x - h, y, z)) / (2.0 * h),
            (schmidt_solid(n, m, x, y + h, z) - schmidt_solid(n, m, x, y - h, z)) / (2.0 * h),
            (schmidt_solid(n, m, x, y, z + h) - schmidt_solid(n, m, x, y, z - h)) / (2.0 * h)};
}

// ---- the oracle, opened at runtime and never linked --------------------------------------------

/// `get_field1_` and `coord_trans_vec1_`, as the vendored `matlab/libirbem.h` documents them.
using GetField1 = void (*)(int*, int*, int*, int*, int*, double*, double*, double*, double*,
                           double*, double*, double*);
using CoordTransVec1 = void (*)(int*, int*, int*, int*, int*, double*, double*, double*);

/// The oracle handle, or nulls when IRBEM is not on this machine.
struct Oracle {
    void* handle = nullptr;
    GetField1 get_field = nullptr;
    CoordTransVec1 coord_trans = nullptr;

    /// Whether the oracle can be called.
    [[nodiscard]] bool usable() const { return get_field != nullptr && coord_trans != nullptr; }
};

/// Open the oracle once for the process. `CHEATAH_SPACE_IRBEM_ORACLE` overrides the path.
const Oracle& oracle() {
    static const Oracle* const opened = [] {
        auto* o = new Oracle;
        const char* env = std::getenv("CHEATAH_SPACE_IRBEM_ORACLE");
        const std::string path = env != nullptr ? env : "/tmp/irbem-builds/libirbem-O2.so";
        o->handle = dlopen(path.c_str(), RTLD_NOW);
        if (o->handle == nullptr) return o;
        o->get_field = reinterpret_cast<GetField1>(dlsym(o->handle, "get_field1_"));
        o->coord_trans = reinterpret_cast<CoordTransVec1>(dlsym(o->handle, "coord_trans_vec1_"));
        return o;
    }();
    return *opened;
}

/// One oracle sample inside the fitted box: GSM position and the oracle's external field there.
struct OracleSample {
    double x, y, z;
    std::array<double, 3> b;
};

/// The oracle's dipole tilt and its external field over a grid inside the fitted box, for one
/// epoch (2015, day @p doy, @p ut seconds) and one driver set.
double oracle_samples(const Oracle& o, int doy, double ut, const Om97Drivers& d,
                      std::vector<OracleSample>& out) {
    int year = 2015;
    int idoy = doy;
    double t = ut;
    int one = 1;
    double ps = 0.0;
    {
        int si = 4;
        int so = 2;
        std::array<double, 3> in{0.0, 0.0, 1.0};
        std::array<double, 3> v{};
        o.coord_trans(&one, &si, &so, &year, &idoy, &t, in.data(), v.data());
        ps = std::atan2(v[0], v[2]);
    }
    out.clear();
    for (int ix = 0; ix <= 6; ++ix)
        for (int iy = 0; iy <= 6; ++iy)
            for (int iz = 0; iz <= 4; ++iz) {
                const double x = -9.0 + (3.0 * ix);
                const double y = -9.0 + (3.0 * iy);
                const double z = -6.0 + (3.0 * iz);
                const double r = std::sqrt((x * x) + (y * y) + (z * z));
                const double xs = (x * std::cos(ps)) - (z * std::sin(ps));
                const double zs = (x * std::sin(ps)) + (z * std::cos(ps));
                if (r < 3.0 || std::sqrt((xs * xs) + (y * y)) > 10.0 || std::fabs(zs) > 7.0) continue;
                std::array<double, 3> gsm{x, y, z};
                std::array<double, 3> geo{};
                {
                    int si = 2;
                    int so = 1;
                    o.coord_trans(&one, &si, &so, &year, &idoy, &t, gsm.data(), geo.data());
                }
                std::array<int, 5> options{0, 0, 0, 0, 0};
                int sysaxes = 1;
                int k0 = 0;
                int k8 = 8;
                std::vector<double> mag(25, 0.0);
                mag[0] = d.kp_times_ten;
                mag[1] = d.dst;
                mag[2] = 5.0;
                mag[3] = 400.0;
                mag[4] = d.pdyn;
                mag[6] = d.bz_imf;
                std::array<double, 3> b0{};
                std::array<double, 3> b8{};
                double m0 = 0.0;
                double m8 = 0.0;
                double x1 = geo[0];
                double x2 = geo[1];
                double x3 = geo[2];
                o.get_field(&k0, options.data(), &sysaxes, &year, &idoy, &t, &x1, &x2, &x3,
                            mag.data(), b0.data(), &m0);
                o.get_field(&k8, options.data(), &sysaxes, &year, &idoy, &t, &x1, &x2, &x3,
                            mag.data(), b8.data(), &m8);
                std::array<double, 3> dgeo{b8[0] - b0[0], b8[1] - b0[1], b8[2] - b0[2]};
                std::array<double, 3> ext{};
                {
                    int si = 1;
                    int so = 2;
                    o.coord_trans(&one, &si, &so, &year, &idoy, &t, dgeo.data(), ext.data());
                }
                out.push_back({x, y, z, ext});
            }
    return ps;
}

/// Solve `A x = b` in place by Gaussian elimination with partial pivoting.
bool solve(std::vector<double>& a, std::vector<double>& b, int n) {
    for (int i = 0; i < n; ++i) {
        int piv = i;
        for (int r = i + 1; r < n; ++r)
            if (std::fabs(a[(r * n) + i]) > std::fabs(a[(piv * n) + i])) piv = r;
        if (std::fabs(a[(piv * n) + i]) < 1e-30) return false;
        if (piv != i) {
            for (int c = 0; c < n; ++c) std::swap(a[(i * n) + c], a[(piv * n) + c]);
            std::swap(b[i], b[piv]);
        }
        for (int r = i + 1; r < n; ++r) {
            const double f = a[(r * n) + i] / a[(i * n) + i];
            for (int c = i; c < n; ++c) a[(r * n) + c] -= f * a[(i * n) + c];
            b[r] -= f * b[i];
        }
    }
    for (int i = n - 1; i >= 0; --i) {
        double s = b[i];
        for (int c = i + 1; c < n; ++c) s -= a[(i * n) + c] * b[c];
        b[i] = s / a[(i * n) + i];
    }
    return true;
}

/// The header's 17 basis fields at a GSM point, rotated back to GSM.
Om97Basis<double> basis_gsm(double ps, double x, double y, double z) {
    const double sp = std::sin(ps);
    const double cp = std::cos(ps);
    const double k = 1.0 / ir::om97_length_scale_re;
    Om97Basis<double> b =
        om97_basis<double>(sp, ((x * cp) - (z * sp)) * k, y * k, ((x * sp) + (z * cp)) * k);
    for (auto& h : b) {
        const double bx = h[0];
        const double bz = h[2];
        h[0] = (bx * cp) + (bz * sp);
        h[2] = (-bx * sp) + (bz * cp);
    }
    return b;
}

}  // namespace

// ================================================================================================
// The published tables
// ================================================================================================

TEST(IrbemOm97, PublishedNormalizationIsTable2) {
    // Ostapenko & Maltsev (1997) Table 2: <A> and sigma_A for Dst (nT), p (nPa), Kp, IMFz (nT).
    const Om97Normalization& n = om97_normalization_published;
    EXPECT_EQ(n.dst_mean, -17.0);
    EXPECT_EQ(n.dst_sigma, 25.0);
    EXPECT_EQ(n.pdyn_mean, 2.2);
    EXPECT_EQ(n.pdyn_sigma, 1.9);
    EXPECT_EQ(n.kp_mean, 2.3);
    EXPECT_EQ(n.kp_sigma, 1.3);
    EXPECT_EQ(n.bz_mean, 0.0);
    EXPECT_EQ(n.bz_sigma, 3.7);
}

TEST(IrbemOm97, MeasuredNormalizationIsCloseToThePublishedOne) {
    // The measured set is the published one carried to more digits, so every dispersion must sit
    // within the paper's own two-figure rounding and every mean within a tenth of a sigma of the
    // printed value. A recovered value outside that would say the recovery measured something
    // else — a different dataset, or a different model.
    const Om97Normalization& p = om97_normalization_published;
    const Om97Normalization& m = om97_normalization_measured;
    EXPECT_NEAR(m.dst_sigma / p.dst_sigma, 1.0, 0.02);
    EXPECT_NEAR(m.pdyn_sigma / p.pdyn_sigma, 1.0, 0.03);
    EXPECT_NEAR(m.kp_sigma / p.kp_sigma, 1.0, 0.05);
    EXPECT_NEAR(m.bz_sigma / p.bz_sigma, 1.0, 0.02);
    EXPECT_NEAR(m.dst_mean, p.dst_mean, 0.1 * p.dst_sigma);
    EXPECT_NEAR(m.pdyn_mean, p.pdyn_mean, 0.1 * p.pdyn_sigma);
    EXPECT_NEAR(m.kp_mean, p.kp_mean, 0.1 * p.kp_sigma);
    EXPECT_NEAR(m.bz_mean, p.bz_mean, 0.1 * p.bz_sigma);
}

TEST(IrbemOm97, RelationCoefficientsAreTable4) {
    // Spot checks against the printed table — the corners and the rows the paper singles out
    // (it names the Dst dependence of the symmetric rows and the IMFz dependence of rows 9, 11
    // and 14 as the strongest effects) — plus the shape of the table.
    ASSERT_EQ(om97_relation_coefficients.size(), 17U);
    ASSERT_EQ(om97_relation_coefficients[0].size(), 5U);
    EXPECT_EQ(om97_relation_coefficients[0][0], -43.39);
    EXPECT_EQ(om97_relation_coefficients[0][4], -1.63);
    EXPECT_EQ(om97_relation_coefficients[3][1], -61.19);
    EXPECT_EQ(om97_relation_coefficients[5][0], -130.10);
    EXPECT_EQ(om97_relation_coefficients[8][4], -10.93);
    EXPECT_EQ(om97_relation_coefficients[9][1], -10.04);
    EXPECT_EQ(om97_relation_coefficients[15][1], -3.19);
    EXPECT_EQ(om97_relation_coefficients[16][0], 4.91);
    EXPECT_EQ(om97_relation_coefficients[16][4], -0.62);
    // The paper: "the coefficients of all the six terms of the azimuthally symmetric part depend
    // mostly on Dst" — for each of rows 1-6 the Dst column dominates the other three drivers.
    for (std::size_t i = 0; i < 6; ++i) {
        const ir::Om97Row& a = om97_relation_coefficients[i];
        EXPECT_GT(std::fabs(a[1]), std::fabs(a[2])) << "row " << (i + 1);
        EXPECT_GT(std::fabs(a[1]), std::fabs(a[3])) << "row " << (i + 1);
        EXPECT_GT(std::fabs(a[1]), std::fabs(a[4])) << "row " << (i + 1);
    }
}

TEST(IrbemOm97, FittedRegionIsThePapers) {
    // The abstract: inner boundary r = 3 Re, outer boundary rho = 10 Re, |z| <= 7 Re; the
    // discussion: invalid for Dst < -200 nT.
    EXPECT_EQ(om97_fitted_region.r_min, 3.0);
    EXPECT_EQ(om97_fitted_region.rho_sm_max, 10.0);
    EXPECT_EQ(om97_fitted_region.abs_z_sm_max, 7.0);
    EXPECT_EQ(om97_fitted_region.dst_min, -200.0);
    EXPECT_EQ(ir::om97_length_scale_re, 10.0);
    EXPECT_EQ(om97_harmonic_count, 17U);
    EXPECT_EQ(ir::om97_regressor_count, 5U);
}

TEST(IrbemOm97, AmplitudesAtTheMeansAreTheConstantColumn) {
    // At the published means every normalized parameter is EXACTLY zero (the subtraction is of
    // identical doubles; 23/10 rounds to the same double as the literal 2.3), so the amplitudes
    // are the a_i0 column bit for bit.
    const std::array<double, om97_harmonic_count> amp = om97_amplitudes<double>(at_means());
    for (std::size_t i = 0; i < om97_harmonic_count; ++i) {
        EXPECT_EQ(amp[i], om97_relation_coefficients[i][0]) << "row " << (i + 1);
    }
    // And the float instantiation is the same numbers rounded once.
    const std::array<float, om97_harmonic_count> f = om97_amplitudes<float>(at_means());
    for (std::size_t i = 0; i < om97_harmonic_count; ++i) {
        EXPECT_EQ(f[i], static_cast<float>(om97_relation_coefficients[i][0])) << "row " << (i + 1);
    }
}

TEST(IrbemOm97, AmplitudesAreLinearInEveryDriver) {
    // One sigma of each driver moves amplitude i by exactly a_ik (to roundoff), two sigma by
    // twice that, and the second difference vanishes. Kp is fed as Kp x 10.
    const Om97Drivers base = at_means();
    const Om97Normalization& n = om97_normalization_published;
    const std::array<double, om97_harmonic_count> a0 = om97_amplitudes<double>(base);
    struct Probe {
        std::size_t column;
        Om97Drivers one;
        Om97Drivers two;
    };
    const std::array<Probe, 4> probes{{
        {1, Om97Drivers{base.dst + n.dst_sigma, base.pdyn, base.kp_times_ten, base.bz_imf},
         Om97Drivers{base.dst + (2.0 * n.dst_sigma), base.pdyn, base.kp_times_ten, base.bz_imf}},
        {2, Om97Drivers{base.dst, base.pdyn + n.pdyn_sigma, base.kp_times_ten, base.bz_imf},
         Om97Drivers{base.dst, base.pdyn + (2.0 * n.pdyn_sigma), base.kp_times_ten, base.bz_imf}},
        {3, Om97Drivers{base.dst, base.pdyn, base.kp_times_ten + (10.0 * n.kp_sigma), base.bz_imf},
         Om97Drivers{base.dst, base.pdyn, base.kp_times_ten + (20.0 * n.kp_sigma), base.bz_imf}},
        {4, Om97Drivers{base.dst, base.pdyn, base.kp_times_ten, base.bz_imf + n.bz_sigma},
         Om97Drivers{base.dst, base.pdyn, base.kp_times_ten, base.bz_imf + (2.0 * n.bz_sigma)}},
    }};
    for (const Probe& p : probes) {
        const std::array<double, om97_harmonic_count> a1 = om97_amplitudes<double>(p.one);
        const std::array<double, om97_harmonic_count> a2 = om97_amplitudes<double>(p.two);
        for (std::size_t i = 0; i < om97_harmonic_count; ++i) {
            EXPECT_NEAR(a1[i] - a0[i], om97_relation_coefficients[i][p.column], 1e-9)
                << "row " << (i + 1) << " column " << p.column;
            EXPECT_NEAR((a2[i] - a1[i]) - (a1[i] - a0[i]), 0.0, 1e-9) << "row " << (i + 1);
        }
    }
    // A caller-supplied normalization is honoured: doubling every sigma halves every slope.
    const Om97Normalization wide{n.dst_mean, 2.0 * n.dst_sigma, n.pdyn_mean, 2.0 * n.pdyn_sigma,
                                 n.kp_mean, 2.0 * n.kp_sigma, n.bz_mean, 2.0 * n.bz_sigma};
    const std::array<double, om97_harmonic_count> w = om97_amplitudes<double>(probes[0].one, wide);
    for (std::size_t i = 0; i < om97_harmonic_count; ++i) {
        EXPECT_NEAR(w[i] - a0[i], 0.5 * om97_relation_coefficients[i][1], 1e-9) << "row " << (i + 1);
    }
}

// ================================================================================================
// The harmonics
// ================================================================================================

TEST(IrbemOm97, CartesianFormsMatchThePrintedCylindricalTable) {
    // The header's Cartesian rows against Table 1 as printed, at scattered off-axis points, for
    // every row: cos phi = x/rho and sin phi = y/rho substituted by hand cannot have slipped
    // anywhere without this noticing. Relative to the row's own magnitude, because the rows span
    // four orders of magnitude over the box.
    std::uint64_t s = 0x2545F4914F6CDD1DULL;
    const auto next = [&s] {
        s = (s * 6364136223846793005ULL) + 1442695040888963407ULL;
        return (static_cast<double>(s >> 11) / 9007199254740992.0) - 0.5;
    };
    double worst = 0.0;
    for (int trial = 0; trial < 400; ++trial) {
        const double x = 2.0 * next();
        const double y = 2.0 * next();
        const double z = 1.4 * next();
        const double sp = next();
        const double rho = std::sqrt((x * x) + (y * y));
        if (rho < 1e-3) continue;
        const double phi = std::atan2(y, x);
        const Om97Basis<double> b = om97_basis<double>(sp, x, y, z);
        for (std::size_t i = 1; i <= om97_harmonic_count; ++i) {
            const std::array<double, 3> cyl = printed_row(i, rho, phi, z, sp);
            const double c = std::cos(phi);
            const double sn = std::sin(phi);
            const std::array<double, 3> ref{(cyl[0] * c) - (cyl[1] * sn), (cyl[0] * sn) + (cyl[1] * c),
                                            cyl[2]};
            const double scale = 1.0 + std::fabs(ref[0]) + std::fabs(ref[1]) + std::fabs(ref[2]);
            for (int k = 0; k < 3; ++k) {
                const double err = std::fabs(b[i - 1][static_cast<std::size_t>(k)] - ref[static_cast<std::size_t>(k)]) / scale;
                worst = std::max(worst, err);
                EXPECT_LT(err, 1e-13) << "row " << i << " component " << k;
            }
        }
    }
    std::printf("[ MEASURED ] worst relative |Cartesian - printed cylindrical| over 17 rows: %.3e\n",
                worst);
}

TEST(IrbemOm97, CurlFreeHarmonicsAreSchmidtNormalizedSolidHarmonics) {
    // The paper: the curl-free harmonics "are expressed in terms of the associated Legendre
    // functions in the Schmidt normalization". Ten rows are curl-free — 1, 2, 3, 7, 8 and the
    // five tilt rows — and each must be the gradient of r^n S_n^m cos(m phi) for the (n, m) below,
    // computed here by an independent Legendre recurrence and a finite difference. Row 17 is the
    // one this test exists for: the table prints its factor as 1/(8 sqrt 15) where Schmidt
    // normalization gives sqrt(15)/8, and the header follows the paper's text, not the misprint.
    struct Row {
        std::size_t i;
        int n;
        int m;
        bool tilt;
    };
    const std::array<Row, 10> rows{{{1, 1, 0, false},
                                    {2, 3, 0, false},
                                    {3, 5, 0, false},
                                    {7, 2, 1, false},
                                    {8, 4, 1, false},
                                    {13, 2, 0, true},
                                    {14, 4, 0, true},
                                    {15, 1, 1, true},
                                    {16, 3, 1, true},
                                    {17, 5, 1, true}}};
    const std::array<std::array<double, 3>, 5> pts{{{0.35, -0.2, 0.15},
                                                    {-0.7, 0.4, -0.5},
                                                    {0.1, 0.9, 0.6},
                                                    {0.8, 0.05, -0.65},
                                                    {-0.45, -0.55, 0.3}}};
    const double sp = 0.6;
    for (const Row& r : rows) {
        for (const std::array<double, 3>& p : pts) {
            const Om97Basis<double> b = om97_basis<double>(sp, p[0], p[1], p[2]);
            const std::array<double, 3> g = schmidt_gradient(r.n, r.m, p[0], p[1], p[2]);
            const double f = r.tilt ? sp : 1.0;
            for (int k = 0; k < 3; ++k) {
                EXPECT_NEAR(b[r.i - 1][static_cast<std::size_t>(k)], f * g[static_cast<std::size_t>(k)], 2e-8)
                    << "row " << r.i << " (n=" << r.n << ", m=" << r.m << ") component " << k;
            }
        }
    }
    // And the seven that are NOT curl-free have a curl: rows 4, 5, 6, 9, 10, 11, 12. The table's
    // last column gives mu0 j_phi for each; the sign of the axially symmetric ring-current rows
    // (4, 5) is westward (negative) as the paper says the total current is.
    const double x = 0.6;
    const double y = 0.0;
    const double z = 0.3;
    const double h = 1e-5;
    for (std::size_t i : {4U, 5U, 6U, 9U, 10U, 11U, 12U}) {
        // curl_y = dB_x/dz - dB_z/dx, which at y = 0 (phi = 0) is -mu0 j_phi's sign convention
        // up to orientation; only non-vanishing is asserted here, the sign for rows 4 and 5.
        const double dbx_dz = (om97_basis<double>(sp, x, y, z + h)[i - 1][0] -
                               om97_basis<double>(sp, x, y, z - h)[i - 1][0]) / (2.0 * h);
        const double dbz_dx = (om97_basis<double>(sp, x + h, y, z)[i - 1][2] -
                               om97_basis<double>(sp, x - h, y, z)[i - 1][2]) / (2.0 * h);
        const double curl_y = dbx_dz - dbz_dx;
        EXPECT_GT(std::fabs(curl_y), 1e-3) << "row " << i << " should carry a current";
        if (i == 4U || i == 5U) {
            EXPECT_LT(curl_y, 0.0) << "row " << i << ": mu0 j_phi = -2 rho, -4 rho^3";
        }
    }
}

TEST(IrbemOm97, EveryHarmonicIsDivergenceFree) {
    // Each of the 17 rows alone, by a central difference at h = 1e-4 in normalized coordinates:
    // the residual is the stencil's truncation, O(h^2 x third derivative / 6) — for the degree-4
    // rows (3, 5, 6, 8, 14, 17) up to ~2e-7 — against row magnitudes of order one. A wrong sign
    // in any Cartesian rewrite makes it order one.
    const double h = 1e-4;
    const double sp = 0.45;
    for (std::size_t i = 0; i < om97_harmonic_count; ++i) {
        for (double x : {-0.8, -0.3, 0.2, 0.7}) {
            for (double y : {-0.6, 0.1, 0.5}) {
                for (double z : {-0.5, 0.0, 0.4}) {
                    const double d =
                        ((om97_basis<double>(sp, x + h, y, z)[i][0] - om97_basis<double>(sp, x - h, y, z)[i][0]) +
                         (om97_basis<double>(sp, x, y + h, z)[i][1] - om97_basis<double>(sp, x, y - h, z)[i][1]) +
                         (om97_basis<double>(sp, x, y, z + h)[i][2] - om97_basis<double>(sp, x, y, z - h)[i][2])) /
                        (2.0 * h);
                    EXPECT_LT(std::fabs(d), 1e-6) << "row " << (i + 1) << " at " << x << "," << y << "," << z;
                }
            }
        }
    }
}

TEST(IrbemOm97, DivergenceVanishesEverywhere) {
    // The whole field, every tilt and corpus regime, over the fitted box: the stencil residual
    // must fall as h^2 over three decades — the signature of an exactly divergence-free field
    // sampled by a second-order stencil. Every harmonic is a polynomial, so the only way this
    // fails is a transcription error, and it cannot be fooled by agreeing with a wrong reference.
    const double d2 = worst_divergence(1e-2);
    const double d3 = worst_divergence(1e-3);
    const double d4 = worst_divergence(1e-4);
    std::printf("[ MEASURED ] worst |div B|: h=1e-2 %.3e  h=1e-3 %.3e  h=1e-4 %.3e nT/Re\n", d2, d3, d4);
    EXPECT_GT(d2 / d3, 50.0) << "divergence is not falling as h^2 — a harmonic is wrong";
    EXPECT_GT(d3 / d4, 30.0) << "divergence is not falling as h^2 — a harmonic is wrong";
    EXPECT_LT(d4, 1e-7);
}

TEST(IrbemOm97, ZeroTiltIsMirrorSymmetricAboutTheEquator) {
    // At psi = 0 the five tilt rows vanish (sin psi = 0 EXACTLY) and every remaining row is even
    // in z for B_z and odd for B_x, B_y. Bitwise, because the two evaluations perform the same
    // operations on operands that differ only in one sign bit.
    for (const corpus::MagInput& m : corpus::regime_drivers) {
        const std::array<double, om97_harmonic_count> amp = om97_amplitudes<double>(drivers_of(m));
        for (double x : {-9.0, -4.5, -1.5, 3.25, 8.0}) {
            for (double y : {-7.5, -1.25, 0.0, 2.5, 6.0}) {
                for (double z : {0.5, 2.25, 6.0}) {
                    const std::array<double, 3> up = om97_components<double>(amp, 0.0, 1.0, x, y, z);
                    const std::array<double, 3> dn = om97_components<double>(amp, 0.0, 1.0, x, y, -z);
                    EXPECT_EQ(up[0], -dn[0]);
                    EXPECT_EQ(up[1], -dn[1]);
                    EXPECT_EQ(up[2], dn[2]);
                }
            }
        }
    }
}

TEST(IrbemOm97, DawnDuskSymmetryHoldsAtEveryTilt) {
    // The model has no dawn-dusk asymmetry: every B_x and B_z is even in y and every B_y odd,
    // at ANY tilt (the asymmetric rows are cos phi in B_rho, B_z and sin phi in B_phi). So
    // y -> -y gives (B_x, -B_y, B_z) exactly.
    for (double ps : kTilts) {
        const double s = std::sin(ps);
        const double c = std::cos(ps);
        for (const corpus::MagInput& m : corpus::regime_drivers) {
            const std::array<double, om97_harmonic_count> amp = om97_amplitudes<double>(drivers_of(m));
            for (double x : {-8.0, -3.5, 2.0, 7.75}) {
                for (double y : {0.75, 3.5, 8.25}) {
                    for (double z : {-5.5, -1.0, 0.0, 4.25}) {
                        const std::array<double, 3> dusk = om97_components<double>(amp, s, c, x, y, z);
                        const std::array<double, 3> dawn = om97_components<double>(amp, s, c, x, -y, z);
                        EXPECT_EQ(dusk[0], dawn[0]);
                        EXPECT_EQ(dusk[1], -dawn[1]);
                        EXPECT_EQ(dusk[2], dawn[2]);
                    }
                }
            }
        }
    }
}

TEST(IrbemOm97, OnTheNoonMidnightMeridianTheFieldStaysInThatPlane) {
    // With y = 0 every B_y term carries a factor y, so B_y is EXACTLY zero on the meridian, at
    // every tilt and every activity — and the dipole axis (rho = 0) is inside the domain of the
    // Cartesian form, where the printed cylindrical form would divide by zero.
    for (double ps : kTilts) {
        for (const corpus::MagInput& m : corpus::regime_drivers) {
            const std::array<double, om97_harmonic_count> amp = om97_amplitudes<double>(drivers_of(m));
            for (double x : {-9.0, -2.0, 0.0, 4.5, 9.0}) {
                for (double z : {-6.0, 0.0, 3.5}) {
                    const std::array<double, 3> b =
                        om97_components<double>(amp, std::sin(ps), std::cos(ps), x, 0.0, z);
                    EXPECT_EQ(b[1], 0.0);
                    EXPECT_TRUE(std::isfinite(b[0]) && std::isfinite(b[2]));
                }
            }
        }
    }
    // Exactly on the axis with the tilt rows switched on: still finite, still B_y = 0.
    const std::array<double, 3> axis =
        om97_components<double>(om97_amplitudes<double>(moderate()), 0.0, 1.0, 0.0, 0.0, 5.0);
    EXPECT_EQ(axis[1], 0.0);
    EXPECT_TRUE(std::isfinite(axis[0]) && std::isfinite(axis[2]));
}

// ================================================================================================
// The physics the paper reports, and the drivers
// ================================================================================================

TEST(IrbemOm97, EquatorialBzFollowsThePapersProfiles) {
    // Figure 2 and Figure 5 of the paper, at zero tilt, in the equatorial plane. Under average
    // conditions B_z is depressed to about -40..-60 nT a few Re down the tail and rises to about
    // +30 nT at 9 Re on the dayside (Fig. 2a); on the dawn-dusk meridian it is ~-40 nT near the
    // Earth and ~0 at +-9 Re (Fig. 2b). Dst deepens the whole profile (Fig. 5, top); pressure
    // compresses the dayside (Fig. 5, second row: B_z at X = +9 grows with p).
    const std::array<double, om97_harmonic_count> avg = om97_amplitudes<double>(at_means());
    const auto bz = [](const std::array<double, om97_harmonic_count>& a, double x, double y) {
        return om97_components<double>(a, 0.0, 1.0, x, y, 0.0)[2];
    };
    EXPECT_GT(bz(avg, -3.0, 0.0), -70.0);
    EXPECT_LT(bz(avg, -3.0, 0.0), -30.0);
    EXPECT_GT(bz(avg, 9.0, 0.0), 15.0);
    EXPECT_LT(bz(avg, 9.0, 0.0), 45.0);
    EXPECT_GT(bz(avg, -9.0, 0.0), -40.0);
    EXPECT_LT(bz(avg, -9.0, 0.0), 0.0);
    EXPECT_GT(bz(avg, 0.0, 9.0), -15.0);
    EXPECT_LT(bz(avg, 0.0, 9.0), 15.0);
    EXPECT_LT(bz(avg, 0.0, 0.0), -30.0);
    // Dst: the paper's headline. Sweep from +8 down to -42 (the three Fig. 5 curves) and beyond:
    // B_z at X = -3 must fall monotonically.
    double prev = std::numeric_limits<double>::infinity();
    for (int k = 0; k <= 10; ++k) {
        const double dst = 8.0 - (10.0 * k);
        const double v = bz(om97_amplitudes<double>(Om97Drivers{dst, 2.2, 23.0, 0.0}), -3.0, 0.0);
        EXPECT_LT(v, prev) << "Dst = " << dst;
        prev = v;
    }
    // Pressure: B_z at X = +9 grows with p (Fig. 5, p = 0.4 vs 4.2 nPa).
    const double low = bz(om97_amplitudes<double>(Om97Drivers{-17.0, 0.4, 23.0, 0.0}), 9.0, 0.0);
    const double high = bz(om97_amplitudes<double>(Om97Drivers{-17.0, 4.2, 23.0, 0.0}), 9.0, 0.0);
    EXPECT_GT(high, low + 10.0);
}

TEST(IrbemOm97, DriversAreContinuousAndSmooth) {
    // The regression is LINEAR in each driver, so across a continuous sweep of any one of them
    // the field moves linearly: the second difference along the sweep is roundoff, and the first
    // difference is bounded by the driver step. This is the property a Kp-binned model (T89)
    // does not have and must not be imitated here: nearby Kp values must give DIFFERENT fields.
    const Position<Frame::GSM> p = at(-5.5, 2.0, 1.5);
    const double tilt = 0.25;
    const Om97Drivers base = moderate();
    const auto field = [&](const Om97Drivers& d) { return om97_field(p, tilt, d).value; };
    struct Sweep {
        const char* name;
        double lo;
        double hi;
        int steps;
    };
    // Bz is swept most densely and asymmetrically southward — it is the driver that couples the
    // solar wind in, and the corpus brief says why the sampling is lopsided.
    const std::array<Sweep, 4> sweeps{{{"Dst", -220.0, 30.0, 50},
                                       {"Pdyn", 0.3, 12.0, 40},
                                       {"Kp", 0.0, 90.0, 45},
                                       {"Bz", -35.0, 5.0, 80}}};
    for (const Sweep& s : sweeps) {
        const double step = (s.hi - s.lo) / s.steps;
        std::vector<double> bz;
        bz.reserve(static_cast<std::size_t>(s.steps) + 1);
        double lipschitz = 0.0;
        for (int k = 0; k <= s.steps; ++k) {
            Om97Drivers d = base;
            const double v = s.lo + (step * k);
            if (s.name == std::string("Dst")) d.dst = v;
            if (s.name == std::string("Pdyn")) d.pdyn = v;
            if (s.name == std::string("Kp")) d.kp_times_ten = v;
            if (s.name == std::string("Bz")) d.bz_imf = v;
            const ir::FieldVector<Frame::GSM> b = field(d);
            bz.push_back(b.v[2]);
            if (k > 0) lipschitz = std::max(lipschitz, std::fabs(bz[static_cast<std::size_t>(k)] - bz[static_cast<std::size_t>(k) - 1]) / step);
        }
        double worst_second = 0.0;
        for (std::size_t k = 2; k < bz.size(); ++k) {
            worst_second = std::max(worst_second, std::fabs((bz[k] - bz[k - 1]) - (bz[k - 1] - bz[k - 2])));
        }
        std::printf("[ MEASURED ] %s sweep %g..%g: |dBz/d%s| <= %.4f nT per unit, worst second "
                    "difference %.2e nT\n",
                    s.name, s.lo, s.hi, s.name, lipschitz, worst_second);
        EXPECT_LT(worst_second, 1e-9) << s.name << " is not linear";
        EXPECT_GT(lipschitz, 1e-4) << s.name << " does not move the field at all";
        EXPECT_LT(lipschitz, 50.0) << s.name;
        // Nearby values give NEARBY fields, never identical ones: the T89 bin identity would be
        // a bug here.
        EXPECT_NE(bz[0], bz[1]) << s.name;
    }
}

TEST(IrbemOm97, SouthwardBzSweepReconfiguresTheNightside) {
    // Row 9's Bz column is the largest IMF sensitivity in Table 4 (-10.93): a uniform B_x
    // proportional to z. So under southward Bz the field above the equator on the nightside
    // acquires a sunward B_x component that grows linearly as Bz goes south — the tail-stretching
    // signature. Sample 60 values from +5 to -40 nT and assert the monotone trend.
    const Position<Frame::GSM> p = at(-6.0, 0.0, 3.0);
    double prev = -std::numeric_limits<double>::infinity();
    for (int k = 0; k <= 60; ++k) {
        const double bz_imf = 5.0 - (0.75 * k);
        Om97Drivers d = moderate();
        d.bz_imf = bz_imf;
        const double bx = om97_field(p, 0.0, d).value.v[0];
        EXPECT_GT(bx, prev) << "Bz = " << bz_imf;
        prev = bx;
    }
}

TEST(IrbemOm97, CorpusRegimesEvaluateAndReportTheirCaveats) {
    // The four corpus regimes, all local times, both hemispheres: every regime evaluates to a
    // finite field, and the paper's storm caveat fires exactly where the paper puts it — quiet,
    // moderate and storm (Dst = -150) are inside, extreme (Dst = -350) is reported.
    for (std::size_t r = 0; r < corpus::regime_drivers.size(); ++r) {
        const Om97Drivers d = drivers_of(corpus::regime_drivers[r]);
        double largest = 0.0;
        for (double lt : corpus::local_times) {
            const double phi = (lt / 24.0) * 2.0 * std::numbers::pi;
            for (double z : {-2.0, 0.0, 2.5}) {
                const Position<Frame::GSM> p = at(6.6 * std::cos(phi), 6.6 * std::sin(phi), z);
                const ir::Result<ir::FieldVector<Frame::GSM>> f = om97_field(p, 0.2, d);
                EXPECT_TRUE(std::isfinite(f.value.v[0]) && std::isfinite(f.value.v[1]) &&
                            std::isfinite(f.value.v[2]));
                EXPECT_EQ(f.status, r == 3 ? Status::OutOfValidityRange : Status::Ok)
                    << "regime " << r << " local time " << lt;
                largest = std::max(largest, f.value.magnitude());
            }
        }
        // Tens of nT somewhere around geosynchronous in every regime — individual points can sit
        // near a null of the external field (the quiet regime has one at 0.24 nT), which is why
        // this is the maximum over the ring and not a per-point floor.
        EXPECT_GT(largest, 10.0) << "regime " << r;
    }
    // The four named storm events too: every one is Dst < -200 except Starlink-2022, and the
    // status says so; every one still returns a number.
    for (const corpus::StormEvent& e : corpus::storm_events) {
        const ir::Result<ir::FieldVector<Frame::GSM>> f =
            om97_field(at(-6.0, 1.0, 0.5), 0.1, drivers_of(e.mag));
        EXPECT_EQ(f.status, e.peak_dst < -200.0 ? Status::OutOfValidityRange : Status::Ok) << e.name;
        EXPECT_TRUE(std::isfinite(f.value.v[2])) << e.name;
    }
}

TEST(IrbemOm97, ReferenceLaneMatchesTheComponentForm) {
    const std::array<double, om97_harmonic_count> amp = om97_amplitudes<double>(moderate());
    const std::array<double, 3> c = om97_components<double>(amp, std::sin(0.3), std::cos(0.3), -4.0, 2.5, 1.5);
    const ir::FieldVector<Frame::GSM> f = om97_field_at(at(-4.0, 2.5, 1.5), std::sin(0.3), std::cos(0.3), amp);
    EXPECT_EQ(f.v[0], c[0]);
    EXPECT_EQ(f.v[1], c[1]);
    EXPECT_EQ(f.v[2], c[2]);
}

// ================================================================================================
// Validity, from both sides
// ================================================================================================

TEST(IrbemOm97, PositionValidityIsCheckedFromBothSides) {
    const Om97Drivers d = moderate();
    // The inner boundary r = 3: on it is inside, just under is reported — and the value is
    // still returned, non-zero.
    EXPECT_EQ(om97_field(at(3.0, 0.0, 0.0), 0.0, d).status, Status::Ok);
    const ir::Result<ir::FieldVector<Frame::GSM>> inner = om97_field(at(2.999, 0.0, 0.0), 0.0, d);
    EXPECT_EQ(inner.status, Status::OutOfValidityRange);
    EXPECT_NE(inner.value.v[2], 0.0);
    // The outer boundary rho_SM = 10, at zero tilt where SM is GSM.
    EXPECT_EQ(om97_field(at(0.0, 10.0, 0.0), 0.0, d).status, Status::Ok);
    EXPECT_EQ(om97_field(at(0.0, 10.001, 0.0), 0.0, d).status, Status::OutOfValidityRange);
    // |z_SM| = 7.
    EXPECT_EQ(om97_field(at(4.0, 0.0, 7.0), 0.0, d).status, Status::Ok);
    EXPECT_EQ(om97_field(at(4.0, 0.0, -7.001), 0.0, d).status, Status::OutOfValidityRange);
    // The box is stated in SM: at a 30-degree tilt a GSM point with z = 6.5 that rotates to
    // z_SM = 6.5 cos 30 + 5 sin 30 = 8.13 is outside, while its mirror is inside.
    const double ps = std::numbers::pi / 6.0;
    EXPECT_EQ(om97_field(at(5.0, 0.0, 6.5), ps, d).status, Status::OutOfValidityRange);
    EXPECT_EQ(om97_field(at(-5.0, 0.0, 6.5), ps, d).status, Status::Ok);
    // Inside the Earth is a refusal, with a zero.
    const ir::Result<ir::FieldVector<Frame::GSM>> deep = om97_field(at(0.5, 0.0, 0.0), 0.0, d);
    EXPECT_EQ(deep.status, Status::DomainError);
    EXPECT_EQ(deep.value.v[2], 0.0);
    // The check itself, on its own.
    EXPECT_EQ(om97_check_fitted_region(3.0, 10.0, 7.0), Status::Ok);
    EXPECT_EQ(om97_check_fitted_region(2.9, 5.0, 1.0), Status::OutOfValidityRange);
    EXPECT_EQ(om97_check_fitted_region(5.0, 10.1, 1.0), Status::OutOfValidityRange);
    EXPECT_EQ(om97_check_fitted_region(5.0, 5.0, 7.1), Status::OutOfValidityRange);
}

TEST(IrbemOm97, DriverValidityIsCheckedFromBothSides) {
    const Position<Frame::GSM> p = at(-5.0, 1.0, 1.0);
    Om97Drivers d = moderate();
    d.dst = -200.0;
    EXPECT_EQ(om97_field(p, 0.1, d).status, Status::Ok);
    d.dst = -200.001;
    const ir::Result<ir::FieldVector<Frame::GSM>> storm = om97_field(p, 0.1, d);
    EXPECT_EQ(storm.status, Status::OutOfValidityRange);
    EXPECT_NE(storm.value.v[2], 0.0) << "the value is still returned";
    // A non-finite driver is a refusal, through status.hpp's check for this model's row.
    d = moderate();
    d.kp_times_ten = std::numeric_limits<double>::quiet_NaN();
    EXPECT_EQ(om97_check_drivers(d), Status::DomainError);
    EXPECT_EQ(om97_field(p, 0.1, d).status, Status::DomainError);
    d = moderate();
    d.bz_imf = std::numeric_limits<double>::infinity();
    EXPECT_EQ(om97_field(p, 0.1, d).status, Status::DomainError);
    // The IRBEM table publishes no ranges, so a wild but finite Pdyn is merely computed.
    d = moderate();
    d.pdyn = 60.0;
    EXPECT_EQ(om97_check_drivers(d), Status::Ok);
}

TEST(IrbemOm97, NonFiniteInputIsADomainError) {
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const double inf = std::numeric_limits<double>::infinity();
    const Om97Drivers d = moderate();
    EXPECT_EQ(om97_field(at(nan, 0.0, 0.0), 0.0, d).status, Status::DomainError);
    EXPECT_EQ(om97_field(at(5.0, inf, 0.0), 0.0, d).status, Status::DomainError);
    EXPECT_EQ(om97_field(at(5.0, 0.0, -inf), 0.0, d).status, Status::DomainError);
    EXPECT_EQ(om97_field(at(5.0, 0.0, 0.0), nan, d).status, Status::DomainError);
    EXPECT_EQ(om97_field(at(5.0, 0.0, 0.0), nan, d).value.v[2], 0.0);
}

TEST(IrbemOm97, RightAngleTiltIsNotRefusedButBeyondItIs) {
    // No tan(psi) anywhere: a tilt of exactly +-pi/2 is an (unphysical) angle to an axis and the
    // rotation is perfectly defined, so it is evaluated. Past it is not an angle to an axis.
    const Om97Drivers d = moderate();
    EXPECT_NE(om97_field(at(-5.0, 1.0, 1.0), std::numbers::pi / 2.0, d).status, Status::DomainError);
    EXPECT_NE(om97_field(at(-5.0, 1.0, 1.0), -std::numbers::pi / 2.0, d).status, Status::DomainError);
    EXPECT_EQ(om97_field(at(-5.0, 1.0, 1.0), 1.6, d).status, Status::DomainError);
    EXPECT_EQ(om97_field(at(-5.0, 1.0, 1.0), -1.6, d).status, Status::DomainError);
}

TEST(IrbemOm97, AnOverflowingExtrapolationIsADomainErrorNotANaN) {
    // A fourth-order polynomial in x/10 overflows binary64 near |x| ~ 1e78 R_E — nothing a trace
    // reaches, exactly what a units bug produces. The answer must be a refusal, never an inf or
    // a NaN, and the batch fold must agree.
    const Om97Drivers d = moderate();
    const ir::Result<ir::FieldVector<Frame::GSM>> far = om97_field(at(1e80, 0.0, 0.0), 0.2, d);
    EXPECT_EQ(far.status, Status::DomainError);
    EXPECT_EQ(far.value.v[0], 0.0);
    EXPECT_EQ(far.value.v[2], 0.0);
    // Merely far is a caveat with a (large, finite) value.
    const ir::Result<ir::FieldVector<Frame::GSM>> out = om97_field(at(60.0, 0.0, 0.0), 0.2, d);
    EXPECT_EQ(out.status, Status::OutOfValidityRange);
    EXPECT_TRUE(std::isfinite(out.value.v[2]));
}

TEST(IrbemOm97, ContextOverloadAgreesWithTheExplicitOne) {
    const ir::Epoch epoch{2015.5, 43200.0, 2015, 180};
    ir::RotationTable identity{};
    for (cheatah::fixarray::mat3d& m : identity) m = cheatah::fixarray::mat3d::identity();
    ir::DriverSet drivers{};
    drivers[static_cast<std::size_t>(ir::Driver::Kp)] = 47.0;
    drivers[static_cast<std::size_t>(ir::Driver::Dst)] = -63.0;
    drivers[static_cast<std::size_t>(ir::Driver::Pdyn)] = 4.5;
    drivers[static_cast<std::size_t>(ir::Driver::BzIMF)] = -8.5;
    const ir::ContextResult built = ir::make_field_context(epoch, -0.42, identity, drivers);
    ASSERT_TRUE(built.has_value()) << ir::describe(built.error());

    const Position<Frame::GSM> p = at(3.75, -1.5, 2.25);
    const ir::Result<ir::FieldVector<Frame::GSM>> via_context = om97_field(p, built.value());
    const ir::Result<ir::FieldVector<Frame::GSM>> via_scalars =
        om97_field(p, -0.42, Om97Drivers{-63.0, 4.5, 47.0, -8.5});
    EXPECT_EQ(via_context.status, via_scalars.status);
    EXPECT_EQ(via_context.value.v[0], via_scalars.value.v[0]);
    EXPECT_EQ(via_context.value.v[1], via_scalars.value.v[1]);
    EXPECT_EQ(via_context.value.v[2], via_scalars.value.v[2]);
    // And with the measured normalization threaded through — a different number, same status.
    const ir::Result<ir::FieldVector<Frame::GSM>> measured =
        om97_field(p, built.value(), om97_normalization_measured);
    EXPECT_EQ(measured.status, via_scalars.status);
    EXPECT_NE(measured.value.v[2], via_scalars.value.v[2]);
}

// ================================================================================================
// The batch lanes
// ================================================================================================

TEST(IrbemOm97, HostFloatLaneTracksTheReferenceLane) {
    const std::size_t n = 2048;
    const std::vector<Position<Frame::GSM>> pts = scatter(n);
    std::vector<float> pos(3 * n);
    for (std::size_t i = 0; i < n; ++i) {
        pos[(3 * i) + 0] = static_cast<float>(pts[i].v[0]);
        pos[(3 * i) + 1] = static_cast<float>(pts[i].v[1]);
        pos[(3 * i) + 2] = static_cast<float>(pts[i].v[2]);
    }
    const double ps = 0.31;
    const Om97Drivers d = moderate();
    std::vector<float> out(3 * n);
    ASSERT_TRUE(om97_field_host(pos, out, static_cast<float>(std::sin(ps)),
                                static_cast<float>(std::cos(ps)), om97_amplitudes<float>(d)));
    const std::array<double, om97_harmonic_count> amp = om97_amplitudes<double>(d);
    double worst = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        const ir::FieldVector<Frame::GSM> ref = om97_field_at(pts[i], std::sin(ps), std::cos(ps), amp);
        for (std::size_t c = 0; c < 3; ++c) worst = std::max(worst, std::fabs(out[(3 * i) + c] - ref.v[c]));
    }
    std::printf("[ MEASURED ] fp32 host lane vs fp64 reference, %zu points: max |dB| = %.3e nT\n", n, worst);
    EXPECT_LT(worst, 1e-3);
}

TEST(IrbemOm97, HostFloatLaneRejectsMismatchedSpans) {
    std::vector<float> pos(7, 1.0F);
    std::vector<float> out(7, 0.0F);
    const std::array<float, om97_harmonic_count> amp = om97_amplitudes<float>(moderate());
    EXPECT_FALSE(om97_field_host(pos, out, 0.0F, 1.0F, amp));
    pos.resize(6);
    EXPECT_FALSE(om97_field_host(pos, out, 0.0F, 1.0F, amp));
    out.resize(6);
    EXPECT_TRUE(om97_field_host(pos, out, 0.0F, 1.0F, amp));
    EXPECT_TRUE(om97_field_host(std::span<const float>{}, std::span<float>{}, 0.0F, 1.0F, amp));
}

TEST(IrbemOm97, ParameterBlockCarriesTheTiltThenTheAmplitudes) {
    const std::array<float, om97_harmonic_count> amp = om97_amplitudes<float>(moderate());
    const std::array<float, om97_param_count> block = om97_param_block(0.25F, 0.75F, amp);
    EXPECT_EQ(om97_param_count, 19U);
    EXPECT_EQ(block[0], 0.25F);
    EXPECT_EQ(block[1], 0.75F);
    for (std::size_t k = 0; k < om97_harmonic_count; ++k) EXPECT_EQ(block[2 + k], amp[k]);
#if CHEATAH_SPACE_IRBEM_OM97_GPU
    EXPECT_EQ(ir::gpu::kernel_info("irbem_om97_f32").params, om97_param_count);
    EXPECT_EQ(ir::gpu::kernel_info("irbem_om97_f32").bindings, 4U);
#endif
}

TEST(IrbemOm97, BatchAgreesWithTheReferenceLane) {
    // Below the crossover the batch is the host fp64 loop, and then it is BIT-identical to the
    // scalar lane.
    const std::vector<Position<Frame::GSM>> pts = scatter(16);
    std::vector<ir::FieldVector<Frame::GSM>> out(pts.size());
    const Om97Drivers d = moderate();
    const ir::Result<bool> r = om97_field_batch(pts, 0.28, d, out);
    EXPECT_EQ(r.status, Status::Ok);
    EXPECT_FALSE(r.value);
    for (std::size_t i = 0; i < pts.size(); ++i) {
        const ir::Result<ir::FieldVector<Frame::GSM>> s = om97_field(pts[i], 0.28, d);
        EXPECT_EQ(out[i].v[0], s.value.v[0]);
        EXPECT_EQ(out[i].v[1], s.value.v[1]);
        EXPECT_EQ(out[i].v[2], s.value.v[2]);
    }
    // The empty batch is the drivers' verdict and nothing else.
    const ir::Result<bool> none = om97_field_batch(std::span<const Position<Frame::GSM>>{}, 0.28, d,
                                                   std::span<ir::FieldVector<Frame::GSM>>{});
    EXPECT_EQ(none.status, Status::Ok);
    Om97Drivers deep = d;
    deep.dst = -300.0;
    EXPECT_EQ(om97_field_batch(std::span<const Position<Frame::GSM>>{}, 0.28, deep,
                               std::span<ir::FieldVector<Frame::GSM>>{}).status,
              Status::OutOfValidityRange);
}

TEST(IrbemOm97, BatchRejectsMismatchedSpans) {
    const std::vector<Position<Frame::GSM>> pts = scatter(8);
    std::vector<ir::FieldVector<Frame::GSM>> out(7);
    const Om97Drivers d = moderate();
    EXPECT_EQ(om97_field_batch(pts, 0.28, d, out).status, Status::DomainError);
    out.resize(8);
    EXPECT_EQ(om97_field_batch(pts, std::numeric_limits<double>::quiet_NaN(), d, out).status,
              Status::DomainError);
    EXPECT_EQ(om97_field_batch(pts, 1.7, d, out).status, Status::DomainError);
    Om97Drivers bad = d;
    bad.dst = std::numeric_limits<double>::infinity();
    EXPECT_EQ(om97_field_batch(pts, 0.28, bad, out).status, Status::DomainError);
}

TEST(IrbemOm97, BatchReportsTheSameEnvelopeTheScalarLaneDoes) {
    const Om97Drivers d = moderate();
    const double ps = 0.33;
    // A batch entirely inside the box is Ok.
    std::vector<Position<Frame::GSM>> pts = scatter(64);
    std::vector<ir::FieldVector<Frame::GSM>> out(pts.size());
    EXPECT_EQ(om97_field_batch(pts, ps, d, out).status, Status::Ok);
    // One point outside — each of the three bounds in turn — and the batch says so, with every
    // point still computed and equal to the scalar lane.
    // (3, 0, 7.5) rotates to z_SM = 3 sin 0.33 + 7.5 cos 0.33 = 8.06; its mirror (-3, 0, 7.5)
    // rotates to 6.12 and would be INSIDE — the box is an SM box, and this is where that shows.
    for (const Position<Frame::GSM>& stray : {at(2.5, 0.5, 0.5), at(0.0, 10.5, 0.0), at(3.0, 0.0, 7.5)}) {
        pts[17] = stray;
        const ir::Result<bool> r = om97_field_batch(pts, ps, d, out);
        EXPECT_EQ(r.status, Status::OutOfValidityRange);
        EXPECT_EQ(om97_field(stray, ps, d).status, Status::OutOfValidityRange);
        for (std::size_t i = 0; i < pts.size(); ++i) {
            EXPECT_EQ(out[i].v[2], om97_field(pts[i], ps, d).value.v[2]) << i;
        }
    }
    // One point inside the Earth, or not finite: the batch is refused and zeroed.
    for (const Position<Frame::GSM>& bad : {at(0.3, 0.1, 0.0), at(std::numeric_limits<double>::quiet_NaN(), 0.0, 0.0)}) {
        pts[17] = bad;
        std::fill(out.begin(), out.end(), ir::FieldVector<Frame::GSM>{cheatah::fixarray::vec3d{7.0, 7.0, 7.0}});
        const ir::Result<bool> r = om97_field_batch(pts, ps, d, out);
        EXPECT_EQ(r.status, Status::DomainError);
        for (const ir::FieldVector<Frame::GSM>& b : out) {
            EXPECT_EQ(b.v[0], 0.0);
            EXPECT_EQ(b.v[2], 0.0);
        }
    }
    // The fold's own extremes, spelled out.
    ir::Om97PositionFold fold{std::sin(ps), std::cos(ps)};
    fold.add(at(5.0, 0.0, 0.0));
    EXPECT_EQ(fold.verdict(), Status::Ok);
    fold.add(at(-3.0, 0.0, 7.5));
    EXPECT_EQ(fold.verdict(), Status::Ok) << "z_SM = 6.12: inside the SM box though z_GSM = 7.5";
    fold.add(at(6.0, 0.0, 6.0));
    EXPECT_EQ(fold.verdict(), Status::OutOfValidityRange) << "z_SM = 7.6: r = 8.5 but z_SM decides";
}

TEST(IrbemOm97, NothingOnTheHeapInTheHotPath) {
    const std::vector<Position<Frame::GSM>> pts = scatter(256);
    std::vector<ir::FieldVector<Frame::GSM>> out(pts.size());
    const Om97Drivers d = moderate();
    (void)om97_field_batch(pts, 0.2, d, out);
    (void)om97_field(at(5.0, 1.0, 1.0), 0.2, d);

    std::vector<float> pos(3 * pts.size(), 4.0F);
    std::vector<float> f32(3 * pts.size());

    const std::size_t before = cheatah_space_test::allocation_count();
    for (int i = 0; i < 64; ++i) {
        sink = sink + om97_field(at(5.0 + (0.01 * i), 1.0, 1.0), 0.2, d).value.v[2];
    }
    (void)om97_field_batch(pts, 0.2, d, out);  // below the crossover: the host lane
    const std::array<float, om97_harmonic_count> amp = om97_amplitudes<float>(d);
    (void)om97_field_host(pos, f32, 0.2F, 0.98F, amp);
    sink = sink + out[0].v[2] + f32[2];
    EXPECT_EQ(before, cheatah_space_test::allocation_count());
}

// ================================================================================================
// The total field
// ================================================================================================

TEST(IrbemOm97, TotalFieldSuperposesInternalAndExternal) {
    const ir::Igrf<10> igrf = ir::Igrf<10>::at(2015.0).value();
    const auto rot = ir::api::rotations_at(2015, 180, 43200.0, igrf);
    ASSERT_EQ(rot.status, Status::Ok);
    const Om97Drivers d = moderate();
    const ir::TotalFieldOm97<10> total(igrf, rot.value, d);
    static_assert(ir::GeoFieldModel<ir::TotalFieldOm97<10>>);
    EXPECT_EQ(ir::TotalFieldOm97<10>::degree, 10);

    const ir::Position<Frame::GEO> p{cheatah::fixarray::vec3d{6.0, 0.0, 0.0}};
    const double internal_only = igrf.evaluate(p).magnitude();
    const double with_external = total.evaluate(p).magnitude();
    EXPECT_EQ(total.drivers().dst, d.dst);
    EXPECT_EQ(total.g(1, 0), igrf.g(1, 0));
    EXPECT_EQ(total.h(1, 1), igrf.h(1, 1));
    EXPECT_EQ(&total.internal(), &igrf);
    EXPECT_EQ(&total.rotations(), &rot.value);
    EXPECT_GT(std::fabs(with_external - internal_only) / internal_only, 1e-3);
    EXPECT_LT(std::fabs(with_external - internal_only) / internal_only, 0.5);
    // The superposition is the sum of the two, in GEO, exactly.
    const ir::Position<Frame::GSM> p_gsm = ir::transform<Frame::GSM>(p, rot.value);
    const double tilt = rot.value.dipole_tilt_deg * (std::numbers::pi / 180.0);
    const ir::FieldVector<Frame::GEO> ext = ir::transform<Frame::GEO>(om97_field(p_gsm, tilt, d).value, rot.value);
    const ir::FieldVector<Frame::GEO> sum{igrf.evaluate(p).v + ext.v};
    EXPECT_NEAR(total.evaluate(p).v[2], sum.v[2], 1e-9);
    EXPECT_EQ(total.external_status(p), Status::Ok);

    // A trace through the total field closes, and activity moves its invariants: the moderate
    // and storm regimes give different B_min at geosynchronous.
    const ir::Result<ir::FieldLine> line = ir::trace_invariant(total, p, 45.0);
    EXPECT_EQ(line.status, Status::Ok);
    const ir::TotalFieldOm97<10> stormy(igrf, rot.value, drivers_of(corpus::regime_drivers[2]));
    const ir::Result<ir::FieldLine> storm_line = ir::trace_invariant(stormy, p, 45.0);
    EXPECT_EQ(storm_line.status, Status::Ok);
    EXPECT_NE(line.value.b_min, storm_line.value.b_min);
    // The measured normalization threads through the wrapper as well.
    const ir::TotalFieldOm97<10> measured(igrf, rot.value, d, om97_normalization_measured);
    EXPECT_NE(measured.evaluate(p).v[2], total.evaluate(p).v[2]);
}

TEST(IrbemOm97, TotalFieldReportsWhenTheExternalModelDeclines) {
    const ir::Igrf<10> igrf = ir::Igrf<10>::at(2015.0).value();
    const auto rot = ir::api::rotations_at(2015, 180, 43200.0, igrf);
    ASSERT_EQ(rot.status, Status::Ok);
    const ir::TotalFieldOm97<10> total(igrf, rot.value, moderate());
    // Inside the fitted box: Ok. At 1.5 Re — the inner belt, where a trace spends most of its
    // steps — the paper's fit does not reach, and the wrapper says so while still adding the
    // polynomial's (small) extrapolation.
    const ir::Position<Frame::GEO> low{cheatah::fixarray::vec3d{1.5, 0.0, 0.0}};
    EXPECT_EQ(total.external_status(low), Status::OutOfValidityRange);
    EXPECT_NE(total.evaluate(low).v[2], igrf.evaluate(low).v[2]);
    // Inside the Earth the external model refuses and the internal field is returned alone.
    const ir::Position<Frame::GEO> deep{cheatah::fixarray::vec3d{0.5, 0.0, 0.0}};
    EXPECT_EQ(total.external_status(deep), Status::DomainError);
    EXPECT_EQ(total.evaluate(deep).v[2], igrf.evaluate(deep).v[2]);
}

// ================================================================================================
// The differential against the IRBEM oracle
// ================================================================================================

TEST(IrbemOm97, AgreesWithTheIrbemOracleToTheTablesRounding) {
    const Oracle& o = oracle();
    if (!o.usable()) {
        GTEST_SKIP() << "IRBEM oracle not present (set CHEATAH_SPACE_IRBEM_ORACLE to its .so); "
                        "the oracle is a dev-only black box and is never linked";
    }
    // Three epochs spanning the tilt range, all four corpus regimes, inside the fitted box. Two
    // caps per regime, both MEASUREMENTS from tools/oracle/ostapenko_diff.cpp with ~1.4x headroom:
    // the published Table 2 (1-2% of the field: the table's two-figure rounding) and the measured
    // normalization (the residual rounding of the two-decimal a_ik). A regression in THIS
    // implementation — a wrong row, a wrong frame — fails both by an order of magnitude.
    struct Epoch {
        int doy;
        double ut;
    };
    const std::array<Epoch, 3> epochs{{{80, 39183.0}, {180, 43200.0}, {355, 7200.0}}};
    const std::array<double, 4> cap_published{{0.65, 0.65, 1.25, 2.7}};
    const std::array<double, 4> cap_measured{{0.032, 0.023, 0.063, 0.18}};
    std::vector<OracleSample> samples;
    for (std::size_t r = 0; r < corpus::regime_drivers.size(); ++r) {
        const Om97Drivers d = drivers_of(corpus::regime_drivers[r]);
        const std::array<double, om97_harmonic_count> amp_pub = om97_amplitudes<double>(d);
        const std::array<double, om97_harmonic_count> amp_meas =
            om97_amplitudes<double>(d, om97_normalization_measured);
        double sum2_pub = 0.0;
        double sum2_meas = 0.0;
        double sig2 = 0.0;
        std::size_t n = 0;
        for (const Epoch& e : epochs) {
            const double ps = oracle_samples(o, e.doy, e.ut, d, samples);
            for (const OracleSample& s : samples) {
                const std::array<double, 3> mp = om97_components<double>(amp_pub, std::sin(ps), std::cos(ps), s.x, s.y, s.z);
                const std::array<double, 3> mm = om97_components<double>(amp_meas, std::sin(ps), std::cos(ps), s.x, s.y, s.z);
                for (std::size_t c = 0; c < 3; ++c) {
                    sum2_pub += (mp[c] - s.b[c]) * (mp[c] - s.b[c]);
                    sum2_meas += (mm[c] - s.b[c]) * (mm[c] - s.b[c]);
                    sig2 += s.b[c] * s.b[c];
                }
                ++n;
            }
        }
        ASSERT_GT(n, 0U);
        const double rms_pub = std::sqrt(sum2_pub / static_cast<double>(n));
        const double rms_meas = std::sqrt(sum2_meas / static_cast<double>(n));
        std::printf("[ MEASURED ] regime %zu: RMS |dB| vs IRBEM kext=8 = %.4f nT with Table 2 as printed "
                    "(%.2f%% of the field), %.4f nT with the measured normalization (%.3f%%), %zu points\n",
                    r, rms_pub, 100.0 * std::sqrt(sum2_pub / sig2), rms_meas,
                    100.0 * std::sqrt(sum2_meas / sig2), n);
        EXPECT_LT(rms_pub, cap_published[r]) << "regime " << r;
        EXPECT_LT(rms_meas, cap_measured[r]) << "regime " << r;
    }

    // The functional form, to roundoff: regress the oracle's field at the moderate regime onto
    // the header's 17 basis fields at every epoch. The residual has no floor to hide behind.
    const Om97Drivers d = moderate();
    for (const Epoch& e : epochs) {
        const double ps = oracle_samples(o, e.doy, e.ut, d, samples);
        const int k = static_cast<int>(om97_harmonic_count);
        std::vector<double> ata(static_cast<std::size_t>(k) * k, 0.0);
        std::vector<double> atb(static_cast<std::size_t>(k), 0.0);
        for (const OracleSample& s : samples) {
            const Om97Basis<double> h = basis_gsm(ps, s.x, s.y, s.z);
            for (std::size_t c = 0; c < 3; ++c)
                for (int i = 0; i < k; ++i) {
                    atb[static_cast<std::size_t>(i)] += h[static_cast<std::size_t>(i)][c] * s.b[c];
                    for (int j = 0; j < k; ++j)
                        ata[(static_cast<std::size_t>(i) * om97_harmonic_count) + static_cast<std::size_t>(j)] +=
                            h[static_cast<std::size_t>(i)][c] * h[static_cast<std::size_t>(j)][c];
                }
        }
        ASSERT_TRUE(solve(ata, atb, k));
        double res2 = 0.0;
        for (const OracleSample& s : samples) {
            const Om97Basis<double> h = basis_gsm(ps, s.x, s.y, s.z);
            for (std::size_t c = 0; c < 3; ++c) {
                double f = 0.0;
                for (int i = 0; i < k; ++i) f += h[static_cast<std::size_t>(i)][c] * atb[static_cast<std::size_t>(i)];
                res2 += (f - s.b[c]) * (f - s.b[c]);
            }
        }
        const double rms = std::sqrt(res2 / static_cast<double>(3 * samples.size()));
        std::printf("[ MEASURED ] tilt %.2f deg: oracle regressed onto the header's basis, rms residual %.3e nT\n",
                    ps * 57.29577951308232, rms);
        EXPECT_LT(rms, 1e-9) << "the oracle's functional form is not this header's";
        // And the fitted harmonic-17 amplitude is the Schmidt-normalized one: within the
        // two-decimal rounding of Table 4 (plus the normalization's 1%), not 15x off.
        const double a17 = om97_amplitudes<double>(d, om97_normalization_measured)[16];
        EXPECT_NEAR(atb[16] / a17, 1.0, 0.05) << "row 17's normalization";
    }
}

#if CHEATAH_SPACE_IRBEM_OM97_GPU
namespace {

/// Point `CHEATAH_SPACE_IRBEM_SPV_DIR` somewhere for the life of the object, and put it back.
class SpvDirScope {
  public:
    /// @param dir the directory to advertise as the shader directory.
    explicit SpvDirScope(const std::string& dir) {
        if (const char* prev = std::getenv(kVar)) {
            had_ = true;
            prev_ = prev;
        }
        ::setenv(kVar, dir.c_str(), 1);
    }
    SpvDirScope(const SpvDirScope&) = delete;
    SpvDirScope& operator=(const SpvDirScope&) = delete;
    SpvDirScope(SpvDirScope&&) = delete;
    SpvDirScope& operator=(SpvDirScope&&) = delete;
    /// Restore what was there before.
    ~SpvDirScope() {
        if (had_) {
            ::setenv(kVar, prev_.c_str(), 1);
            return;
        }
        ::unsetenv(kVar);
    }

  private:
    static constexpr const char* kVar = "CHEATAH_SPACE_IRBEM_SPV_DIR";
    bool had_ = false;
    std::string prev_;
};

}  // namespace

TEST(IrbemOm97, BatchFallsBackToTheHostWhenTheShaderWasNeverBuilt) {
    if (!ir::gpu::available()) GTEST_SKIP() << "no device: " << ir::gpu::unavailable_reason();
    const std::size_t n = 4 * ir::gpu::gpu_crossover("irbem_om97_f32");
    const std::vector<Position<Frame::GSM>> pts = scatter(n);
    std::vector<ir::FieldVector<Frame::GSM>> out(n);
    const Om97Drivers d = moderate();
    {
        const SpvDirScope nowhere(std::filesystem::temp_directory_path().string() +
                                  "/cheatah-space-no-such-shader-dir");
        const ir::Result<bool> r = om97_field_batch(pts, 0.28, d, out);
        EXPECT_EQ(r.status, Status::Ok);
        EXPECT_FALSE(r.value) << "with no compiled shader the batch must run on the host";
    }
    const std::array<double, om97_harmonic_count> amp = om97_amplitudes<double>(d);
    for (std::size_t i = 0; i < n; ++i) {
        const ir::FieldVector<Frame::GSM> ref = om97_field_at(pts[i], std::sin(0.28), std::cos(0.28), amp);
        ASSERT_EQ(out[i].v[0], ref.v[0]) << "point " << i;
        ASSERT_EQ(out[i].v[1], ref.v[1]) << "point " << i;
        ASSERT_EQ(out[i].v[2], ref.v[2]) << "point " << i;
    }
}

TEST(IrbemOm97, BatchUsesTheDeviceWhenOneIsAvailable) {
    if (!ir::gpu::available()) GTEST_SKIP() << "no device: " << ir::gpu::unavailable_reason();
    const std::size_t n = 4 * ir::gpu::gpu_crossover("irbem_om97_f32");
    const std::vector<Position<Frame::GSM>> pts = scatter(n);
    std::vector<ir::FieldVector<Frame::GSM>> out(n);
    const Om97Drivers d = moderate();
    const ir::Result<bool> r = om97_field_batch(pts, 0.28, d, out);
    EXPECT_EQ(r.status, Status::Ok);
    EXPECT_TRUE(r.value) << "the batch fell back to the host with a device present";

    const std::array<double, om97_harmonic_count> amp = om97_amplitudes<double>(d);
    double worst = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        const ir::FieldVector<Frame::GSM> ref = om97_field_at(pts[i], std::sin(0.28), std::cos(0.28), amp);
        for (std::size_t c = 0; c < 3; ++c) worst = std::max(worst, std::fabs(out[i].v[c] - ref.v[c]));
    }
    std::printf("[ MEASURED ] device batch of %zu vs fp64 reference: max |dB| = %.3e nT\n", n, worst);
    EXPECT_LT(worst, 1e-2);

    const std::vector<Position<Frame::GSM>> few = scatter(16);
    std::vector<ir::FieldVector<Frame::GSM>> few_out(few.size());
    EXPECT_FALSE(om97_field_batch(few, 0.28, d, few_out).value);
}

TEST(IrbemOm97, TheDeviceLaneRefusesABadPointBeforeItDispatches) {
    if (!ir::gpu::available()) GTEST_SKIP() << "no device: " << ir::gpu::unavailable_reason();
    const std::size_t n = 4 * ir::gpu::gpu_crossover("irbem_om97_f32");
    std::vector<Position<Frame::GSM>> pts = scatter(n);
    pts[n / 3] = at(0.25, 0.0, 0.0);
    std::vector<ir::FieldVector<Frame::GSM>> out(n, ir::FieldVector<Frame::GSM>{cheatah::fixarray::vec3d{7.0, 7.0, 7.0}});
    const Om97Drivers d = moderate();
    const ir::Result<bool> r = om97_field_batch(pts, 0.28, d, out);
    EXPECT_EQ(r.status, Status::DomainError);
    EXPECT_FALSE(r.value) << "a refused batch must not have reached the device";
    for (const ir::FieldVector<Frame::GSM>& b : out) {
        ASSERT_EQ(b.v[0], 0.0);
        ASSERT_EQ(b.v[1], 0.0);
        ASSERT_EQ(b.v[2], 0.0);
    }
    std::vector<Position<Frame::GSM>> far = scatter(n);
    far[n / 2] = at(-12.0, 0.0, 0.0);
    const ir::Result<bool> f = om97_field_batch(far, 0.28, d, out);
    EXPECT_EQ(f.status, Status::OutOfValidityRange);
    EXPECT_TRUE(f.value) << "an out-of-validity batch must still be computed on the device";
}

TEST(IrbemOm97, DeviceKernelAgreesWithTheHostLane) {
    if (!ir::gpu::available()) GTEST_SKIP() << "no device: " << ir::gpu::unavailable_reason();
    const std::size_t n = 1 << 16;
    const std::vector<Position<Frame::GSM>> pts = scatter(n);
    std::vector<float> pos(3 * n);
    for (std::size_t i = 0; i < n; ++i) {
        pos[(3 * i) + 0] = static_cast<float>(pts[i].v[0]);
        pos[(3 * i) + 1] = static_cast<float>(pts[i].v[1]);
        pos[(3 * i) + 2] = static_cast<float>(pts[i].v[2]);
    }
    const double ps = 0.31;
    const float sp = static_cast<float>(std::sin(ps));
    const float cp = static_cast<float>(std::cos(ps));
    const std::array<float, om97_harmonic_count> amp = om97_amplitudes<float>(moderate());
    std::vector<float> host(3 * n);
    std::vector<float> device(3 * n);
    ASSERT_TRUE(om97_field_host(pos, host, sp, cp, amp));
    const std::array<float, om97_param_count> block = om97_param_block(sp, cp, amp);
    ir::gpu::dispatch_batch("irbem_om97_f32", pos, device, std::span<const float>(block));
    double worst = 0.0;
    for (std::size_t i = 0; i < 3 * n; ++i) {
        worst = std::max(worst, std::fabs(static_cast<double>(device[i]) - host[i]));
    }
    std::printf("[ MEASURED ] device vs host, %zu points: max |dB| = %.3e nT\n", n, worst);
    EXPECT_LT(worst, 1e-3) << "the device is not evaluating the same expressions";
}
#endif
