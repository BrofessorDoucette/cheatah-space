// Differential study of space.irbem's Ostapenko & Maltsev (1997) against the IRBEM oracle's
// `kext = 8`, and the experiments that say exactly what the remaining difference is.
//
// DEV-ONLY. Like t89_diff.cpp beside it, this is the one kind of program that touches IRBEM, it is
// never built by the QA gate, and it never ships. IRBEM is LGPL-3.0 and cheatah-space is MIT: the
// library is run here as a BLACK BOX (dlopen plus the documented C entry points), never read for
// its logic and never linked into anything we distribute.
//
// WHAT IT MEASURES.
//
// `space/irbem/ext_ostapenko.hpp` implements the model AS PUBLISHED: eq. (2) over the 17 harmonics
// of Table 1 (with the Schmidt-normalized row 17 the paper's text prescribes), the 85 coefficients
// of Table 4 and the normalization of Table 2. Four passes, in order:
//
//   1. `deviation`   — how far apart are they, per corpus regime (quiet / moderate / storm /
//                      extreme) and per tilt, inside the paper's fitted box, with the published
//                      Table 2 and with the measured normalization below.
//   2. `form`        — regress the oracle's external field onto THIS header's 17 basis fields, all
//                      17 amplitudes free. A residual at roundoff says the functional form (the
//                      Cartesian Table 1, in SM, Schmidt row 17) is the oracle's exactly. Two
//                      control variants (evaluated in GSM; tilt sign flipped) show what a wrong
//                      frame would have looked like.
//   3. `linearity`   — one-step and two-step slopes of every amplitude in every driver; the
//                      regression form is linear, so their difference must be roundoff.
//   4. `normalize`   — with the published a_ik held, recover the eight Table 2 scalars from the
//                      measured intercepts and slopes by least squares, report the residual of
//                      that 8-parameter model over all 85 measured numbers, and re-run the
//                      deviation pass with the recovered values. This is the experiment behind
//                      `om97_normalization_measured`.
//
// The external field is isolated from the oracle by DIFFERENCE: `get_field1_` with `kext = 8`
// minus the same call with `kext = 0`, both with `options(5) = 0`, so the internal IGRF term is
// bit-for-bit identical between them and cancels exactly. The dipole tilt is taken from the oracle
// too — the SM z-axis transformed into GSM is `(sin psi, 0, cos psi)`.
//
// Build (from the repository root):
//   g++ -O2 -std=c++20 tools/oracle/ostapenko_diff.cpp -I. \
//       -I$CHEATAH_DIR/stdlib/ndarray -I$CHEATAH_DIR/stdlib/builtins -I$CHEATAH_DIR/stdlib/fixarray \
//       -o /tmp/ostapenko_diff && /tmp/ostapenko_diff /tmp/irbem-builds/libirbem-O2.so
#include <dlfcn.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

#include "space/irbem/ext_ostapenko.hpp"
#include "tests/irbem_domain_corpus.hpp"

namespace {

namespace ir = cheatah::space::irbem;
using ir::om97_amplitudes;
using ir::om97_basis;
using ir::om97_components;
using ir::om97_harmonic_count;
using ir::Om97Drivers;
using ir::Om97Normalization;

/// `get_field1_`, as documented in the vendored matlab/libirbem.h.
using GetField1 = void (*)(int*, int*, int*, int*, int*, double*, double*, double*, double*,
                           double*, double*, double*);
/// `coord_trans_vec1_`, likewise — used for GSM<->GEO and to read the oracle's own dipole tilt.
using CoordTransVec1 = void (*)(int*, int*, int*, int*, int*, double*, double*, double*);

GetField1 g_get_field = nullptr;
CoordTransVec1 g_coord_trans = nullptr;

/// An epoch to sample at; the year is fixed at 2015 throughout.
struct Epoch {
    int doy;
    double ut;
};

/// Three epochs spanning the tilt range: near zero, strongly positive, strongly negative.
constexpr std::array<Epoch, 3> kEpochs{{{80, 39183.0}, {180, 43200.0}, {355, 7200.0}}};
constexpr int kYear = 2015;

/// One sample point: GSM position, its GEO image, and the oracle's tilt for the epoch.
struct Sample {
    double x, y, z;
    std::array<double, 3> geo;
};

/// The oracle's dipole tilt for an epoch, radians.
double oracle_tilt(const Epoch& e) {
    int year = kYear;
    int doy = e.doy;
    double ut = e.ut;
    int one = 1;
    int si = 4;
    int so = 2;
    std::array<double, 3> in{0.0, 0.0, 1.0};
    std::array<double, 3> out{};
    g_coord_trans(&one, &si, &so, &year, &doy, &ut, in.data(), out.data());
    return std::atan2(out[0], out[2]);
}

/// The sample grid inside the paper's fitted box for one epoch (the box is stated in SM, so the
/// tilt decides which GSM lattice points are inside it).
std::vector<Sample> grid(const Epoch& e, double ps, double step) {
    std::vector<Sample> out;
    int year = kYear;
    int doy = e.doy;
    double ut = e.ut;
    int one = 1;
    const int nx = static_cast<int>(std::lround(18.0 / step));
    const int nz = static_cast<int>(std::lround(12.0 / step));
    for (int ix = 0; ix <= nx; ++ix)
        for (int iy = 0; iy <= nx; ++iy)
            for (int iz = 0; iz <= nz; ++iz) {
                const double x = -9.0 + (step * ix);
                const double y = -9.0 + (step * iy);
                const double z = -6.0 + (step * iz);
                const double r = std::sqrt((x * x) + (y * y) + (z * z));
                const double xs = (x * std::cos(ps)) - (z * std::sin(ps));
                const double zs = (x * std::sin(ps)) + (z * std::cos(ps));
                if (r < 3.0 || std::sqrt((xs * xs) + (y * y)) > 10.0 || std::fabs(zs) > 7.0) continue;
                std::array<double, 3> gsm{x, y, z};
                std::array<double, 3> geo{};
                int si = 2;
                int so = 1;
                g_coord_trans(&one, &si, &so, &year, &doy, &ut, gsm.data(), geo.data());
                out.push_back({x, y, z, geo});
            }
    return out;
}

/// The oracle's external field at one sample, GSM, nT: `kext = 8` minus `kext = 0`.
std::array<double, 3> oracle_external(const Epoch& e, const Sample& s, const Om97Drivers& d) {
    int year = kYear;
    int doy = e.doy;
    double ut = e.ut;
    int one = 1;
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
    double x1 = s.geo[0];
    double x2 = s.geo[1];
    double x3 = s.geo[2];
    g_get_field(&k0, options.data(), &sysaxes, &year, &doy, &ut, &x1, &x2, &x3, mag.data(),
                b0.data(), &m0);
    g_get_field(&k8, options.data(), &sysaxes, &year, &doy, &ut, &x1, &x2, &x3, mag.data(),
                b8.data(), &m8);
    std::array<double, 3> dgeo{b8[0] - b0[0], b8[1] - b0[1], b8[2] - b0[2]};
    std::array<double, 3> gsm{};
    int si = 1;
    int so = 2;
    g_coord_trans(&one, &si, &so, &year, &doy, &ut, dgeo.data(), gsm.data());
    return gsm;
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

/// The 17 basis fields at a GSM point, in GSM, for the header's frame (variant 0), the same
/// polynomials evaluated in GSM without the SM rotation (1), or with the tilt sign flipped (2).
ir::Om97Basis<double> basis_gsm(int variant, double ps, double x, double y, double z) {
    double sp = std::sin(ps);
    const double cp = std::cos(ps);
    if (variant == 2) sp = -sp;
    const double k = 1.0 / ir::om97_length_scale_re;
    double xs = x * k;
    double ys = y * k;
    double zs = z * k;
    if (variant != 1) {
        xs = ((x * cp) - (z * sp)) * k;
        zs = ((x * sp) + (z * cp)) * k;
    }
    ir::Om97Basis<double> b = om97_basis<double>(sp, xs, ys, zs);
    if (variant == 1) return b;
    for (auto& h : b) {
        const double bx = h[0];
        const double bz = h[2];
        h[0] = (bx * cp) + (bz * sp);
        h[2] = (-bx * sp) + (bz * cp);
    }
    return b;
}

/// Least-squares fit of the 17 amplitudes to the oracle's field over @p samples; returns the RMS
/// residual and writes the amplitudes.
double fit_amplitudes(int variant, const Epoch& e, double ps, const std::vector<Sample>& samples,
                      const Om97Drivers& d, std::array<double, om97_harmonic_count>& amp,
                      double* signal_rms) {
    const int n = static_cast<int>(om97_harmonic_count);
    std::vector<double> ata(static_cast<std::size_t>(n) * n, 0.0);
    std::vector<double> atb(static_cast<std::size_t>(n), 0.0);
    std::vector<std::array<double, 3>> obs;
    obs.reserve(samples.size());
    double sig2 = 0.0;
    for (const Sample& s : samples) {
        const std::array<double, 3> o = oracle_external(e, s, d);
        obs.push_back(o);
        const ir::Om97Basis<double> h = basis_gsm(variant, ps, s.x, s.y, s.z);
        for (int c = 0; c < 3; ++c) {
            sig2 += o[c] * o[c];
            for (int i = 0; i < n; ++i) {
                atb[i] += h[i][c] * o[c];
                for (int j = 0; j < n; ++j) ata[(i * n) + j] += h[i][c] * h[j][c];
            }
        }
    }
    std::vector<double> a = ata;
    std::vector<double> b = atb;
    if (!solve(a, b, n)) return 1e30;
    for (int i = 0; i < n; ++i) amp[i] = b[i];
    double res2 = 0.0;
    for (std::size_t s = 0; s < samples.size(); ++s) {
        const ir::Om97Basis<double> h = basis_gsm(variant, ps, samples[s].x, samples[s].y, samples[s].z);
        for (int c = 0; c < 3; ++c) {
            double f = 0.0;
            for (int i = 0; i < n; ++i) f += h[i][c] * amp[i];
            res2 += (f - obs[s][c]) * (f - obs[s][c]);
        }
    }
    if (signal_rms != nullptr) *signal_rms = std::sqrt(sig2 / static_cast<double>(3 * samples.size()));
    return std::sqrt(res2 / static_cast<double>(3 * samples.size()));
}

/// Deviation statistics of the header against the oracle over @p samples.
struct Deviation {
    std::size_t n = 0;
    double rms_abs = 0.0;
    double p99_abs = 0.0;
    double max_abs = 0.0;
    double rms_rel = 0.0;
    double max_rel = 0.0;
};

Deviation deviation(const Epoch& e, double ps, const std::vector<Sample>& samples,
                    const Om97Drivers& d, const Om97Normalization& norm) {
    const std::array<double, om97_harmonic_count> amp = om97_amplitudes<double>(d, norm);
    std::vector<double> abs_err;
    double sum2 = 0.0;
    double sig2 = 0.0;
    double worst_rel = 0.0;
    for (const Sample& s : samples) {
        const std::array<double, 3> o = oracle_external(e, s, d);
        const std::array<double, 3> m =
            om97_components<double>(amp, std::sin(ps), std::cos(ps), s.x, s.y, s.z);
        double d2 = 0.0;
        double o2 = 0.0;
        for (int c = 0; c < 3; ++c) {
            d2 += (m[c] - o[c]) * (m[c] - o[c]);
            o2 += o[c] * o[c];
        }
        abs_err.push_back(std::sqrt(d2));
        sum2 += d2;
        sig2 += o2;
        worst_rel = std::max(worst_rel, std::sqrt(d2 / (o2 + 1e-30)));
    }
    std::sort(abs_err.begin(), abs_err.end());
    Deviation out;
    out.n = samples.size();
    out.rms_abs = std::sqrt(sum2 / static_cast<double>(samples.size()));
    out.p99_abs = abs_err[(samples.size() * 99) / 100];
    out.max_abs = abs_err.back();
    out.rms_rel = std::sqrt(sum2 / sig2);
    out.max_rel = worst_rel;
    return out;
}

/// The corpus regimes' drivers, as Om97Drivers.
Om97Drivers drivers_of(const cheatah_space_test::MagInput& m) {
    return Om97Drivers{m.dst, m.pdyn, m.kp * 10.0, m.bz_imf};
}

const char* regime_name(std::size_t i) {
    static const char* const names[4] = {"quiet", "moderate", "storm", "extreme"};
    return names[i];
}

}  // namespace

int main(int argc, char** argv) {
    const std::string lib = argc > 1 ? argv[1] : "/tmp/irbem-builds/libirbem-O2.so";
    void* h = dlopen(lib.c_str(), RTLD_NOW);
    if (h == nullptr) {
        // NOLINTNEXTLINE(cert-err33-c): a diagnostic on the way out; nothing to do if stderr is gone
        (void)std::fprintf(stderr, "cannot dlopen %s: %s\n", lib.c_str(), dlerror());
        return 1;
    }
    g_get_field = reinterpret_cast<GetField1>(dlsym(h, "get_field1_"));
    g_coord_trans = reinterpret_cast<CoordTransVec1>(dlsym(h, "coord_trans_vec1_"));
    if (g_get_field == nullptr || g_coord_trans == nullptr) {
        // NOLINTNEXTLINE(cert-err33-c): a diagnostic on the way out; nothing to do if stderr is gone
        (void)std::fprintf(stderr, "missing entry points in %s\n", lib.c_str());
        return 1;
    }
    std::printf("space.irbem OM97 (as published, 1997) vs IRBEM kext=8, oracle %s\n", lib.c_str());

    std::array<double, kEpochs.size()> tilts{};
    std::array<std::vector<Sample>, kEpochs.size()> grids;
    for (std::size_t e = 0; e < kEpochs.size(); ++e) {
        tilts[e] = oracle_tilt(kEpochs[e]);
        grids[e] = grid(kEpochs[e], tilts[e], 2.0);
    }

    // ---- 1. deviation, per regime and tilt, with the published Table 2 --------------------------
    const auto deviation_pass = [&](const Om97Normalization& norm, const char* label) {
        std::printf("\n--- deviation inside the fitted box (3 <= r, rho_SM <= 10, |z_SM| <= 7), %s ---\n",
                    label);
        std::printf("%-9s %7s %6s %10s %10s %10s %9s %9s\n", "regime", "tilt", "N", "rms|dB|",
                    "p99|dB|", "max|dB|", "rms rel", "max rel");
        for (std::size_t r = 0; r < cheatah_space_test::regime_drivers.size(); ++r) {
            const Om97Drivers d = drivers_of(cheatah_space_test::regime_drivers[r]);
            for (std::size_t e = 0; e < kEpochs.size(); ++e) {
                const Deviation v = deviation(kEpochs[e], tilts[e], grids[e], d, norm);
                std::printf("%-9s %7.2f %6zu %10.3e %10.3e %10.3e %9.2e %9.2e\n", regime_name(r),
                            tilts[e] * 57.29577951308232, v.n, v.rms_abs, v.p99_abs, v.max_abs,
                            v.rms_rel, v.max_rel);
            }
        }
    };
    deviation_pass(ir::om97_normalization_published, "Table 2 as printed");

    // ---- 2. form: regress the oracle onto the header's basis ----------------------------------
    std::printf("\n--- form: the oracle regressed onto this header's 17 basis fields ---\n");
    const Om97Drivers d_form = drivers_of(cheatah_space_test::regime_drivers[1]);
    for (std::size_t e = 0; e < kEpochs.size(); ++e) {
        for (int variant = 0; variant < 3; ++variant) {
            std::array<double, om97_harmonic_count> amp{};
            double sig = 0.0;
            const double res = fit_amplitudes(variant, kEpochs[e], tilts[e], grids[e], d_form, amp, &sig);
            std::printf("tilt %6.2f deg, %-32s: rms residual %.3e nT (signal rms %.2f nT, N=%zu)\n",
                        tilts[e] * 57.29577951308232,
                        variant == 0   ? "SM frame (the header)"
                        : variant == 1 ? "control: evaluated in GSM"
                                       : "control: tilt sign flipped",
                        res, sig, grids[e].size());
        }
    }

    // ---- 3. linearity + 4. normalization recovery -----------------------------------------------
    // All at the epoch nearest zero tilt so the tilt rows are small and the symmetric rows carry
    // the signal; the form pass has already shown the frame is right at every tilt.
    std::printf("\n--- linearity and the recovered Table 2 ---\n");
    const std::size_t e0 = 0;
    const Om97Normalization& pub = ir::om97_normalization_published;
    const Om97Drivers at_means{pub.dst_mean, pub.pdyn_mean, pub.kp_mean * 10.0, pub.bz_mean};
    std::array<double, om97_harmonic_count> a_means{};
    (void)fit_amplitudes(0, kEpochs[e0], tilts[e0], grids[e0], at_means, a_means, nullptr);

    // Per driver: the measured slope of every amplitude (one-step and two-step).
    struct DriverProbe {
        const char* name;
        double step;   // one sigma of the published Table 2
        int column;    // Table 4 column
    };
    const std::array<DriverProbe, 4> probes{{{"Dst", pub.dst_sigma, 1},
                                             {"Pdyn", pub.pdyn_sigma, 2},
                                             {"Kp", pub.kp_sigma, 3},
                                             {"Bz", pub.bz_sigma, 4}}};
    std::array<std::array<double, om97_harmonic_count>, 4> slopes{};
    double worst_nonlinear = 0.0;
    for (std::size_t k = 0; k < probes.size(); ++k) {
        const auto perturbed = [&](double delta) {
            Om97Drivers d = at_means;
            if (k == 0) d.dst += delta;
            if (k == 1) d.pdyn += delta;
            if (k == 2) d.kp_times_ten += 10.0 * delta;
            if (k == 3) d.bz_imf += delta;
            return d;
        };
        std::array<double, om97_harmonic_count> a1{};
        std::array<double, om97_harmonic_count> a2{};
        (void)fit_amplitudes(0, kEpochs[e0], tilts[e0], grids[e0], perturbed(probes[k].step), a1, nullptr);
        (void)fit_amplitudes(0, kEpochs[e0], tilts[e0], grids[e0], perturbed(2.0 * probes[k].step), a2, nullptr);
        for (std::size_t i = 0; i < om97_harmonic_count; ++i) {
            const double s1 = (a1[i] - a_means[i]) / probes[k].step;
            const double s2 = (a2[i] - a_means[i]) / (2.0 * probes[k].step);
            slopes[k][i] = s1;
            worst_nonlinear = std::max(worst_nonlinear, std::fabs(s1 - s2) * probes[k].step);
        }
    }
    std::printf("worst |one-step - two-step| amplitude difference over 4 drivers x 17 rows: %.3e nT\n",
                worst_nonlinear);

    // sigma_k: the least-squares scale between the published column and the measured slopes.
    std::array<double, 4> sigma{};
    std::array<double, 4> ratio_lo{};
    std::array<double, 4> ratio_hi{};
    for (std::size_t k = 0; k < 4; ++k) {
        double num = 0.0;
        double den = 0.0;
        ratio_lo[k] = 1e30;
        ratio_hi[k] = -1e30;
        for (std::size_t i = 0; i < om97_harmonic_count; ++i) {
            const double a = ir::om97_relation_coefficients[i][probes[k].column];
            num += a * a;
            den += a * slopes[k][i];
            if (std::fabs(a) > 2.0) {  // rows large enough that two-decimal rounding is < 0.25%
                const double ratio = slopes[k][i] * probes[k].step / a;
                ratio_lo[k] = std::min(ratio_lo[k], ratio);
                ratio_hi[k] = std::max(ratio_hi[k], ratio);
            }
        }
        sigma[k] = num / den;
        std::printf("%-5s: measured slope x published sigma / published a_ik spans %.4f..%.4f over "
                    "rows with |a| > 2  ->  sigma = %.4f (published %.1f)\n",
                    probes[k].name, ratio_lo[k], ratio_hi[k], sigma[k], probes[k].step);
    }
    // mu_k: intercepts at the published means minus a_i0 = sum_k a_ik (mu_pub_k - mu_k) / sigma_k,
    // a 17 x 4 linear least-squares problem in the four unknown means.
    {
        std::vector<double> ata(16, 0.0);
        std::vector<double> atb(4, 0.0);
        for (std::size_t i = 0; i < om97_harmonic_count; ++i) {
            std::array<double, 4> row{};
            for (std::size_t k = 0; k < 4; ++k)
                row[k] = ir::om97_relation_coefficients[i][probes[k].column] / sigma[k];
            const double y = a_means[i] - ir::om97_relation_coefficients[i][0];
            for (std::size_t p = 0; p < 4; ++p) {
                atb[p] += row[p] * y;
                for (std::size_t q = 0; q < 4; ++q) ata[(p * 4) + q] += row[p] * row[q];
            }
        }
        std::vector<double> a = ata;
        std::vector<double> b = atb;
        if (!solve(a, b, 4)) {
            std::printf("mean recovery: singular\n");
            return 1;
        }
        const std::array<double, 4> pub_means{pub.dst_mean, pub.pdyn_mean, pub.kp_mean, pub.bz_mean};
        Om97Normalization measured{};
        std::array<double, 4> mu{};
        for (std::size_t k = 0; k < 4; ++k) mu[k] = pub_means[k] - b[k];  // b = mu_pub - mu
        measured = Om97Normalization{mu[0], sigma[0], mu[1], sigma[1], mu[2], sigma[2], mu[3], sigma[3]};
        std::printf("recovered means: Dst %.3f  Pdyn %.4f  Kp %.4f  Bz %.4f (published %.1f %.1f %.1f %.1f)\n",
                    mu[0], mu[1], mu[2], mu[3], pub_means[0], pub_means[1], pub_means[2], pub_means[3]);

        // The residual of the 8-parameter model over all 85 measured numbers, in amplitude units,
        // against the rounding floor of a two-decimal table (uniform +-0.005 -> rms 0.0029).
        double res2 = 0.0;
        double worst = 0.0;
        for (std::size_t i = 0; i < om97_harmonic_count; ++i) {
            double pred = ir::om97_relation_coefficients[i][0];
            for (std::size_t k = 0; k < 4; ++k)
                pred += ir::om97_relation_coefficients[i][probes[k].column] * (pub_means[k] - mu[k]) / sigma[k];
            const double r0 = pred - a_means[i];
            res2 += r0 * r0;
            worst = std::max(worst, std::fabs(r0));
            for (std::size_t k = 0; k < 4; ++k) {
                const double rs = (ir::om97_relation_coefficients[i][probes[k].column] / sigma[k]) - slopes[k][i];
                // express in the same units as the table: slope x sigma
                res2 += (rs * sigma[k]) * (rs * sigma[k]);
                worst = std::max(worst, std::fabs(rs * sigma[k]));
            }
        }
        std::printf("8-parameter model vs the 85 measured numbers: rms %.4f, worst %.4f (a two-decimal "
                    "table's rounding floor is rms 0.0029, worst 0.005)\n",
                    std::sqrt(res2 / 85.0), worst);
        std::printf("=> om97_normalization_measured{%.2f, %.2f, %.2f, %.3f, %.2f, %.3f, %.2f, %.3f}\n",
                    measured.dst_mean, measured.dst_sigma, measured.pdyn_mean, measured.pdyn_sigma,
                    measured.kp_mean, measured.kp_sigma, measured.bz_mean, measured.bz_sigma);

        deviation_pass(measured, "recovered normalization");
        deviation_pass(ir::om97_normalization_measured, "om97_normalization_measured as shipped");
    }
    return 0;
}
