// Differential study and CALIBRATION of space.irbem's T01-storm (Tsyganenko, Singer & Kasper 2003)
// against the IRBEM oracle's `kext = 10`, and the experiment that quantifies the gap.
//
// DEV-ONLY. Like convergence.cpp and t89_diff.cpp beside it, this is the one kind of program that
// touches IRBEM, it is never built by the QA gate, and it never ships. IRBEM is LGPL-3.0 and
// cheatah-space is MIT: the library is run here as a BLACK BOX (dlopen plus the documented C entry
// points), never read for its logic and never linked into anything we distribute. Tsyganenko's own
// T01/TSK03 source is GPL-3.0 and is not read at all, by anything, ever.
//
// WHAT IT DOES, and why this one also CALIBRATES.
//
// The T89 harness compares a published coefficient table against the oracle. TSK03 publishes no
// such table: the 2003 paper states the module structure (magnetopause, cross-tail sheet, symmetric
// and partial ring currents, Birkeland systems) and the driver set (Dst, Pdyn, By, Bz, G2, G3), and
// the numerical coefficients — several hundred, most of them shielding-field expansions — exist
// only in the author's code distribution, which is GPL-3.0. So `space/irbem/ext_t01s.hpp` carries
// the published STRUCTURE as 108 exactly divergence-free spatial modes whose amplitudes are
// smooth in the six drivers, and the amplitude coefficients are determined HERE, by least squares
// against the oracle's kext=10 minus kext=0 field sampled as a black box. Three passes:
//
//   1. `calibrate`: sample the oracle over a driver x position x tilt grid inside the Shue (1998)
//      magnetopause, solve the 1 188 linear coefficients, optionally refine
//      the seven free geometry parameters by Nelder-Mead, and report the training residual.
//   2. `holdout`: the same on fresh driver states, a fresh epoch and an off-grid point set, which
//      is the honest number — a fit reported on its own training set proves nothing.
//   3. `deviation`: the SHIPPED header (whatever coefficients it currently carries) against the
//      oracle, per corpus regime and per real storm event, belts and full box.
//
// The external field is isolated from the oracle by DIFFERENCE: `get_field1_` with `kext = 10`
// minus the same call with `kext = 0`, both with `options(5) = 0` so the internal IGRF term is
// bit-for-bit identical between them and cancels exactly. The dipole tilt is taken from the oracle
// too — the SM z-axis transformed into GSM is `(sin psi, 0, cos psi)` — so a frame difference cannot
// masquerade as a model difference.
//
// Build (from the repository root):
//   g++ -O2 -std=c++20 tools/oracle/t01s_diff.cpp -I. \
//       -I$CHEATAH_DIR/stdlib/ndarray -I$CHEATAH_DIR/stdlib/builtins -I$CHEATAH_DIR/stdlib/fixarray \
//       -o /tmp/t01s_diff && /tmp/t01s_diff /tmp/irbem-builds/libirbem-O2.so [--refine] [--emit]
//
// `--refine` runs the Nelder-Mead geometry refinement (minutes); `--emit` prints the coefficient
// table and geometry as C++ source for pasting into ext_t01s.hpp.
#include <dlfcn.h>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <random>
#include <string>
#include <vector>

#include "space/irbem/ext_t01s.hpp"
#include "tests/irbem_domain_corpus.hpp"

namespace {

using namespace cheatah::space::irbem;

using GetField1 = void (*)(int*, int*, int*, int*, int*, double*, double*, double*, double*,
                           double*, double*, double*);
using CoordTransVec1 = void (*)(int*, int*, int*, int*, int*, double*, double*, double*);

GetField1 g_get_field = nullptr;
CoordTransVec1 g_ctv = nullptr;

struct OracleEpoch {
    int year;
    int doy;
    double ut;
};

/// One oracle sample.
struct Sample {
    double tilt;
    double x, y, z;
    T01sDrivers d;
    std::array<double, 3> b;  ///< oracle external field, GSM, nT
    double weight;
};

double oracle_tilt(const OracleEpoch& e) {
    int one = 1, si = 4, so = 2, iy = e.year, id = e.doy;
    double ut = e.ut;
    std::array<double, 3> in{0.0, 0.0, 1.0}, out{};
    g_ctv(&one, &si, &so, &iy, &id, &ut, in.data(), out.data());
    return std::atan2(out[0], out[2]);
}

std::array<double, 3> oracle_external(const OracleEpoch& e, double x, double y, double z,
                                      const T01sDrivers& d) {
    int one = 1, iy = e.year, id = e.doy;
    double ut = e.ut;
    std::array<double, 3> gsm{x, y, z}, geo{};
    {
        int si = 2, so = 1;
        g_ctv(&one, &si, &so, &iy, &id, &ut, gsm.data(), geo.data());
    }
    std::array<int, 5> options{0, 0, 0, 0, 0};
    int sysaxes = 1, k0 = 0, k10 = 10;
    std::vector<double> mag(25, 0.0);
    mag[1] = d.dst;
    mag[4] = d.pdyn;
    mag[5] = d.by_imf;
    mag[6] = d.bz_imf;
    mag[8] = d.g2;
    mag[9] = d.g3;
    std::array<double, 3> b0{}, b1{};
    double m0 = 0, m1 = 0;
    double x1 = geo[0], x2 = geo[1], x3 = geo[2];
    g_get_field(&k0, options.data(), &sysaxes, &iy, &id, &ut, &x1, &x2, &x3, mag.data(), b0.data(), &m0);
    g_get_field(&k10, options.data(), &sysaxes, &iy, &id, &ut, &x1, &x2, &x3, mag.data(), b1.data(), &m1);
    std::array<double, 3> dgeo{b1[0] - b0[0], b1[1] - b0[1], b1[2] - b0[2]}, out{};
    {
        int si = 1, so = 2;
        g_ctv(&one, &si, &so, &iy, &id, &ut, dgeo.data(), out.data());
    }
    return out;
}

/// Shue et al. (1998), JGR 103, 17691, eqs. (9)-(11): the magnetopause radius at angle theta from
/// the Sun-Earth line. Used only to keep the calibration inside the region where a magnetospheric
/// model means anything; the oracle evaluates its formulas on both sides of the boundary.
bool inside_magnetopause(double x, double y, double z, double pdyn, double bz) {
    const double r = std::sqrt((x * x) + (y * y) + (z * z));
    if (r < 1e-9) return true;
    const double cos_t = x / r;
    const double r0 = (10.22 + (1.29 * std::tanh(0.184 * (bz + 8.14)))) * std::pow(pdyn, -1.0 / 6.6);
    const double alpha = (0.58 - (0.007 * bz)) * (1.0 + (0.024 * std::log(pdyn)));
    const double rmp = r0 * std::pow(2.0 / (1.0 + cos_t), alpha);
    return r < 0.92 * rmp;
}

std::vector<std::array<double, 3>> point_set(bool jitter, unsigned seed) {
    std::vector<std::array<double, 3>> pts;
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> jit(-1.0, 1.0);
    for (double x = -15.0; x <= 12.0; x += 3.0)
        for (double y = -12.0; y <= 12.0; y += 3.0)
            for (double z = -6.0; z <= 6.0; z += 3.0) {
                std::array<double, 3> p{x, y, z};
                if (jitter) for (double& c : p) c += 1.2 * jit(rng);
                pts.push_back(p);
            }
    const std::array<double, 6> shells{2.5, 3.5, 4.5, 5.5, 6.6, 8.0};
    const std::array<double, 5> lats{0.0, 20.0, -20.0, 40.0, -40.0};
    for (double r : shells)
        for (int lt = 0; lt < 8; ++lt)
            for (double lat : lats) {
                const double phi = (lt * 45.0 + (jitter ? 22.5 : 0.0)) * std::numbers::pi / 180.0;
                const double la = (lat + (jitter ? 7.0 : 0.0)) * std::numbers::pi / 180.0;
                pts.push_back({r * std::cos(la) * std::cos(phi), r * std::cos(la) * std::sin(phi),
                               r * std::sin(la)});
            }
    return pts;
}

std::vector<T01sDrivers> driver_states(unsigned seed, int n_random) {
    std::vector<T01sDrivers> out;
    using cheatah_space_test::regime_drivers;
    using cheatah_space_test::storm_events;
    for (const auto& m : regime_drivers) out.push_back({m.dst, m.pdyn, m.by_imf, m.bz_imf, m.g2, m.g3});
    for (const auto& e : storm_events)
        out.push_back({e.mag.dst, e.mag.pdyn, e.mag.by_imf, e.mag.bz_imf, e.mag.g2, e.mag.g3});
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> u(0.0, 1.0);
    for (int i = 0; i < n_random; ++i) {
        // Storm-shaped: Dst and G2/G3 skew towards disturbed values, Bz towards southward.
        const double s = u(rng);
        T01sDrivers d{};
        d.dst = 20.0 - (470.0 * s * u(rng));
        d.pdyn = 0.5 + (24.5 * std::pow(u(rng), 1.5));
        d.by_imf = -20.0 + (40.0 * u(rng));
        d.bz_imf = 12.0 - (47.0 * std::pow(u(rng), 0.8));
        d.g2 = 90.0 * s * u(rng);
        d.g3 = 200.0 * s * u(rng);
        out.push_back(d);
    }
    return out;
}

std::vector<Sample> sample_oracle(const std::vector<OracleEpoch>& epochs,
                                  const std::vector<std::array<double, 3>>& pts,
                                  const std::vector<T01sDrivers>& states, double r_min,
                                  double r_max) {
    std::vector<Sample> out;
    for (const OracleEpoch& e : epochs) {
        const double ps = oracle_tilt(e);
        for (const T01sDrivers& d : states) {
            std::vector<Sample> batch;
            double sig2 = 0.0;
            for (const auto& p : pts) {
                const double r = std::sqrt((p[0] * p[0]) + (p[1] * p[1]) + (p[2] * p[2]));
                if (r < r_min || r > r_max) continue;
                if (!inside_magnetopause(p[0], p[1], p[2], d.pdyn, d.bz_imf)) continue;
                Sample s{ps, p[0], p[1], p[2], d, oracle_external(e, p[0], p[1], p[2], d), 1.0};
                if (std::fabs(s.b[0]) > 1e20) {
                    static int shown = 0;
                    if (shown++ < 8)
                        std::printf("BAD: epoch %d/%d p=(%.1f,%.1f,%.1f) dst %.0f pdyn %.2f by %.1f bz %.1f g2 %.1f g3 %.1f\n",
                                    e.year, e.doy, p[0], p[1], p[2], d.dst, d.pdyn, d.by_imf, d.bz_imf, d.g2, d.g3);
                    continue;
                }
                sig2 += (s.b[0] * s.b[0]) + (s.b[1] * s.b[1]) + (s.b[2] * s.b[2]);
                batch.push_back(s);
            }
            if (batch.empty()) continue;
            // Equalise the influence of driver states: a 300 nT extreme would otherwise drown
            // the 20 nT quiet field the belts spend most of their time in.
            const double rms = std::sqrt(sig2 / static_cast<double>(batch.size()));
            const double w = 1.0 / (rms + 25.0);
            for (Sample& s : batch) s.weight = w;
            out.insert(out.end(), batch.begin(), batch.end());
        }
    }
    return out;
}

constexpr std::size_t kM = t01s_mode_count;
constexpr std::size_t kF = t01s_feature_count;
constexpr std::size_t kN = kM * kF;

using Coeff = std::array<std::array<double, kF>, kM>;

/// Solve `A x = b` in place by Gaussian elimination with partial pivoting.
bool solve(std::vector<double>& a, std::vector<double>& b, std::size_t n) {
    for (std::size_t i = 0; i < n; ++i) {
        std::size_t piv = i;
        for (std::size_t r = i + 1; r < n; ++r)
            if (std::fabs(a[(r * n) + i]) > std::fabs(a[(piv * n) + i])) piv = r;
        if (std::fabs(a[(piv * n) + i]) < 1e-300) return false;
        if (piv != i) {
            for (std::size_t c = 0; c < n; ++c) std::swap(a[(i * n) + c], a[(piv * n) + c]);
            std::swap(b[i], b[piv]);
        }
        for (std::size_t r = i + 1; r < n; ++r) {
            const double f = a[(r * n) + i] / a[(i * n) + i];
            if (f == 0.0) continue;
            for (std::size_t c = i; c < n; ++c) a[(r * n) + c] -= f * a[(i * n) + c];
            b[r] -= f * b[i];
        }
    }
    for (std::size_t i = n; i-- > 0;) {
        double s = b[i];
        for (std::size_t c = i + 1; c < n; ++c) s -= a[(i * n) + c] * b[c];
        b[i] = s / a[(i * n) + i];
    }
    return true;
}

/// Weighted least squares for the mode x feature coefficients with @p geo held. Samples arrive
/// grouped by driver state (sample_oracle emits them that way), and every sample of a group shares
/// one feature vector, so the normal matrix is a sum of Kronecker products
/// `M_s (x) f_s f_s^T` with `M_s` the group's mode Gram matrix — which turns an O(N (M F)^2)
/// accumulation into O(N M^2 + groups (M F)^2). Returns the weighted RMS residual.
double fit_linear(const std::vector<Sample>& samples, const T01sGeometry& geo, Coeff& coeff,
                  double ridge) {
    std::vector<double> ata(kN * kN, 0.0), atb(kN, 0.0);
    std::array<std::array<double, 3>, kM> modes{};
    std::vector<double> gram(kM * kM), rhs(kM);
    double wsum = 0.0;
    std::size_t i = 0;
    while (i < samples.size()) {
        // one group: same drivers (and therefore the same feature vector and weight)
        std::size_t j = i;
        const T01sDrivers d0 = samples[i].d;
        const double w = samples[i].weight * samples[i].weight;
        std::fill(gram.begin(), gram.end(), 0.0);
        std::fill(rhs.begin(), rhs.end(), 0.0);
        while (j < samples.size() && std::memcmp(&samples[j].d, &d0, sizeof d0) == 0 &&
               samples[j].weight == samples[i].weight) {
            const Sample& s = samples[j];
            t01s_mode_fields<double>(geo, std::sin(s.tilt), std::cos(s.tilt), s.x, s.y, s.z, modes);
            for (int c = 0; c < 3; ++c) {
                for (std::size_t a = 0; a < kM; ++a) {
                    const double ma = modes[a][static_cast<std::size_t>(c)];
                    rhs[a] += ma * s.b[static_cast<std::size_t>(c)];
                    for (std::size_t b = 0; b < kM; ++b) gram[(a * kM) + b] += ma * modes[b][static_cast<std::size_t>(c)];
                }
            }
            wsum += 3.0 * w;
            ++j;
        }
        const std::array<double, kF> f = t01s_features(d0);
        for (std::size_t a = 0; a < kM; ++a) {
            for (std::size_t k = 0; k < kF; ++k) atb[(a * kF) + k] += w * rhs[a] * f[k];
            for (std::size_t b = 0; b < kM; ++b) {
                const double g = w * gram[(a * kM) + b];
                if (g == 0.0) continue;
                for (std::size_t k = 0; k < kF; ++k)
                    for (std::size_t l = 0; l < kF; ++l)
                        ata[(((a * kF) + k) * kN) + (b * kF) + l] += g * f[k] * f[l];
            }
        }
        i = j;
    }
    // Column scaling + a ridge so a nearly collinear pair of modes cannot blow the amplitudes up.
    std::vector<double> sc(kN);
    for (std::size_t q = 0; q < kN; ++q) sc[q] = ata[(q * kN) + q] > 0.0 ? 1.0 / std::sqrt(ata[(q * kN) + q]) : 0.0;
    std::vector<double> a(kN * kN), b(kN);
    for (std::size_t q = 0; q < kN; ++q) {
        b[q] = atb[q] * sc[q];
        for (std::size_t r = 0; r < kN; ++r) a[(q * kN) + r] = ata[(q * kN) + r] * sc[q] * sc[r];
        a[(q * kN) + q] += ridge;
        if (sc[q] == 0.0) a[(q * kN) + q] = 1.0;
    }
    if (!solve(a, b, kN)) return 1e30;
    for (std::size_t m = 0; m < kM; ++m)
        for (std::size_t k = 0; k < kF; ++k) coeff[m][k] = b[(m * kF) + k] * sc[(m * kF) + k];
    double res = 0.0;
    for (const Sample& s : samples) {
        t01s_mode_fields<double>(geo, std::sin(s.tilt), std::cos(s.tilt), s.x, s.y, s.z, modes);
        const std::array<double, kF> f = t01s_features(s.d);
        for (int c = 0; c < 3; ++c) {
            double v = 0.0;
            for (std::size_t m = 0; m < kM; ++m) {
                double amp = 0.0;
                for (std::size_t k = 0; k < kF; ++k) amp += coeff[m][k] * f[k];
                v += amp * modes[m][static_cast<std::size_t>(c)];
            }
            const double e = v - s.b[static_cast<std::size_t>(c)];
            res += s.weight * s.weight * e * e;
        }
    }
    return std::sqrt(res / wsum);
}

struct Stats {
    std::size_t n = 0;
    double rms_abs = 0, p99_abs = 0, max_abs = 0, rms_rel = 0, p99_rel = 0;
};

/// Statistics of a model (given coefficients + geometry) against samples.
Stats evaluate(const std::vector<Sample>& samples, const T01sGeometry& geo, const Coeff& coeff) {
    std::vector<double> abs_e, rel_e;
    double sum2 = 0, sig2 = 0;
    std::array<std::array<double, 3>, kM> modes{};
    for (const Sample& s : samples) {
        t01s_mode_fields<double>(geo, std::sin(s.tilt), std::cos(s.tilt), s.x, s.y, s.z, modes);
        const std::array<double, kF> f = t01s_features(s.d);
        double d2 = 0, o2 = 0;
        for (int c = 0; c < 3; ++c) {
            double v = 0.0;
            for (std::size_t m = 0; m < kM; ++m) {
                double amp = 0.0;
                for (std::size_t k = 0; k < kF; ++k) amp += coeff[m][k] * f[k];
                v += amp * modes[m][static_cast<std::size_t>(c)];
            }
            const double e = v - s.b[static_cast<std::size_t>(c)];
            d2 += e * e;
            o2 += s.b[static_cast<std::size_t>(c)] * s.b[static_cast<std::size_t>(c)];
        }
        abs_e.push_back(std::sqrt(d2));
        rel_e.push_back(std::sqrt(d2) / (std::sqrt(o2) + 1e-12));
        sum2 += d2;
        sig2 += o2;
    }
    Stats st;
    st.n = abs_e.size();
    if (st.n == 0) return st;
    std::sort(abs_e.begin(), abs_e.end());
    std::sort(rel_e.begin(), rel_e.end());
    const std::size_t p99 = (st.n * 99) / 100;
    st.rms_abs = std::sqrt(sum2 / static_cast<double>(st.n));
    st.p99_abs = abs_e[p99];
    st.max_abs = abs_e[st.n - 1];
    st.rms_rel = std::sqrt(sum2 / sig2);
    st.p99_rel = rel_e[p99];
    return st;
}

void print_stats(const char* tag, const Stats& s) {
    std::printf("%-34s %6zu  rms %8.3f  p99 %8.3f  max %8.3f nT | rms rel %6.3f  p99 rel %6.3f\n",
                tag, s.n, s.rms_abs, s.p99_abs, s.max_abs, s.rms_rel, s.p99_rel);
}

// The seven geometry parameters the refinement is allowed to move. `hinge2` is fixed at the
// published (4 Re)^2 of T89 eq. (11) and is not one of them; the radial scales `a` and half
// thicknesses `D` now live per mode in the mode tables rather than in the geometry, which is what
// changed when the model grew from fifteen modes to a hundred and eight.
std::array<double, 7> pack(const T01sGeometry& g) {
    return {g.delta_y, g.r_hinge, g.warp_g, g.l_y, g.x0_w, g.dx_w, g.dy_w};
}
T01sGeometry unpack(const T01sGeometry& base, const std::array<double, 7>& v) {
    T01sGeometry g = base;
    g.delta_y = v[0];
    g.r_hinge = v[1];
    g.warp_g = v[2];
    g.l_y = v[3];
    g.x0_w = v[4];
    g.dx_w = v[5];
    g.dy_w = v[6];
    return g;
}

/// Nelder-Mead over the seven free geometry parameters with the linear fit inside.
T01sGeometry refine(const std::vector<Sample>& samples, const T01sGeometry& start_geo, double ridge,
                    int iters) {
    constexpr int nd = 7;
    Coeff tmp{};
    const auto cost = [&](const std::array<double, 7>& v) {
        // delta_y >= 0 (a sheet that thins toward the flanks is not the published shape),
        // r_hinge, l_y, dx_w, dy_w all strictly positive, and |warp_g| bounded by the published
        // range of T89 Table 1.
        if (v[0] < 0.0 || v[0] > 0.2 || v[1] < 2.0 || v[1] > 20.0 || std::fabs(v[2]) > 12.0 ||
            v[3] < 3.0 || v[4] < -10.0 || v[4] > 20.0 || v[5] < 3.0 || v[6] < 5.0)
            return 1e30;
        return fit_linear(samples, unpack(start_geo, v), tmp, ridge);
    };
    std::vector<std::array<double, 7>> sim;
    std::vector<double> val;
    sim.push_back(pack(start_geo));
    val.push_back(cost(sim[0]));
    for (int i = 0; i < nd; ++i) {
        auto v = sim[0];
        v[static_cast<std::size_t>(i)] *= 1.2;
        sim.push_back(v);
        val.push_back(cost(v));
    }
    for (int it = 0; it < iters; ++it) {
        std::vector<int> idx(nd + 1);
        for (int i = 0; i <= nd; ++i) idx[static_cast<std::size_t>(i)] = i;
        std::sort(idx.begin(), idx.end(), [&](int a, int b) { return val[static_cast<std::size_t>(a)] < val[static_cast<std::size_t>(b)]; });
        std::vector<std::array<double, 7>> s2;
        std::vector<double> v2;
        for (int i : idx) {
            s2.push_back(sim[static_cast<std::size_t>(i)]);
            v2.push_back(val[static_cast<std::size_t>(i)]);
        }
        sim = s2;
        val = v2;
        if (it % 10 == 0) {
            std::printf("  refine it %3d: best %.4f  geo a_t %.2f D0 %.2f a_rc1 %.2f a_rc2 %.2f Rc %.2f x0 %.2f Dy %.2f\n",
                        it, val[0], sim[0][0], sim[0][1], sim[0][2], sim[0][3], sim[0][4], sim[0][5], sim[0][6]);
            std::fflush(stdout);
        }
        if (std::fabs(val[nd] - val[0]) < 1e-5 * val[0]) break;
        std::array<double, 7> centre{};
        for (int i = 0; i < nd; ++i)
            for (int j = 0; j < nd; ++j) centre[static_cast<std::size_t>(j)] += sim[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] / nd;
        const auto along = [&](double a) {
            std::array<double, 7> v{};
            for (int j = 0; j < nd; ++j)
                v[static_cast<std::size_t>(j)] = centre[static_cast<std::size_t>(j)] + (a * (centre[static_cast<std::size_t>(j)] - sim[nd][static_cast<std::size_t>(j)]));
            return v;
        };
        const auto vr = along(1.0);
        const double fr = cost(vr);
        if (fr < val[0]) {
            const auto ve = along(2.0);
            const double fe = cost(ve);
            sim[nd] = fe < fr ? ve : vr;
            val[nd] = std::min(fe, fr);
        } else if (fr < val[nd - 1]) {
            sim[nd] = vr;
            val[nd] = fr;
        } else {
            const auto vc = along(-0.5);
            const double fc = cost(vc);
            if (fc < val[nd]) {
                sim[nd] = vc;
                val[nd] = fc;
            } else {
                for (int i = 1; i <= nd; ++i) {
                    for (int j = 0; j < nd; ++j)
                        sim[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] = sim[0][static_cast<std::size_t>(j)] + (0.5 * (sim[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] - sim[0][static_cast<std::size_t>(j)]));
                    val[static_cast<std::size_t>(i)] = cost(sim[static_cast<std::size_t>(i)]);
                }
            }
        }
    }
    return unpack(start_geo, sim[0]);
}

void emit(const T01sGeometry& g, const Coeff& c) {
    std::printf("\n// ---- paste into ext_t01s.hpp ----\n");
    std::printf("inline constexpr T01sGeometry t01s_geometry{\n"
                "    /* delta_y */ %.6g, /* r_hinge */ %.6g, /* warp_g */ %.6g, /* l_y */ %.6g,\n"
                "    /* hinge2 */ %.6g,  /* x0_w */ %.6g,    /* dx_w */ %.6g,  /* dy_w */ %.6g};\n",
                g.delta_y, g.r_hinge, g.warp_g, g.l_y, g.hinge2, g.x0_w, g.dx_w, g.dy_w);
    std::printf("inline constexpr std::array<std::array<double, t01s_feature_count>, t01s_mode_count>\n"
                "    t01s_coefficients{{\n");
    for (std::size_t m = 0; m < kM; ++m) {
        std::printf("        {{");
        for (std::size_t k = 0; k < kF; ++k) std::printf("%s%.9g", k ? ", " : "", c[m][k]);
        std::printf("}},\n");
    }
    std::printf("    }};\n");
}

}  // namespace

int main(int argc, char** argv) {
    std::string lib = "/tmp/irbem-builds/libirbem-O2.so";
    bool do_refine = false, do_emit = false;
    double ridge_arg = 1e-3;  // the ridge the SHIPPED table was fit with; see ext_t01s.hpp
    for (int i = 1; i < argc; ++i) {
        if (std::strcmp(argv[i], "--refine") == 0) do_refine = true;
        else if (std::strcmp(argv[i], "--emit") == 0) do_emit = true;
        else if (std::strncmp(argv[i], "--ridge=", 8) == 0) ridge_arg = std::atof(argv[i] + 8);
        else lib = argv[i];
    }
    void* h = dlopen(lib.c_str(), RTLD_NOW);
    if (h == nullptr) {
        std::fprintf(stderr, "cannot dlopen %s: %s\n", lib.c_str(), dlerror());
        return 1;
    }
    g_get_field = reinterpret_cast<GetField1>(dlsym(h, "get_field1_"));
    g_ctv = reinterpret_cast<CoordTransVec1>(dlsym(h, "coord_trans_vec1_"));
    if (g_get_field == nullptr || g_ctv == nullptr) {
        std::fprintf(stderr, "missing entry points in %s\n", lib.c_str());
        return 1;
    }
    std::printf("space.irbem T01-storm vs IRBEM kext=10, oracle %s\n", lib.c_str());

    // Three epochs spanning the tilt range (same as t89_diff): near zero, strongly +, strongly -.
    const std::vector<OracleEpoch> train_epochs{{2015, 80, 39183.0}, {2015, 180, 43200.0}, {2015, 355, 7200.0}};
    const std::vector<OracleEpoch> hold_epochs{{2003, 303, 21600.0}, {2015, 240, 60000.0}};
    for (const OracleEpoch& e : train_epochs) std::printf("  train tilt %+.2f deg\n", oracle_tilt(e) * 180.0 / std::numbers::pi);
    for (const OracleEpoch& e : hold_epochs) std::printf("  holdout tilt %+.2f deg\n", oracle_tilt(e) * 180.0 / std::numbers::pi);

    const auto t0 = std::chrono::steady_clock::now();
    const std::vector<Sample> train =
        sample_oracle(train_epochs, point_set(false, 1), driver_states(12345, 40), 2.2, 40.0);
    const std::vector<Sample> hold =
        sample_oracle(hold_epochs, point_set(true, 7), driver_states(98765, 16), 2.2, 40.0);
    const auto t1 = std::chrono::steady_clock::now();
    std::printf("sampled %zu training + %zu holdout points in %.1f s\n", train.size(), hold.size(),
                std::chrono::duration<double>(t1 - t0).count());

    const double ridge = ridge_arg;
    T01sGeometry geo = t01s_geometry;
    Coeff coeff{};

    // ---- deviation pass of the SHIPPED header, first: what the file currently does --------------
    std::printf("\n--- shipped header vs oracle (holdout set) ---\n");
    print_stats("shipped, all holdout", evaluate(hold, t01s_geometry, t01s_coefficients));

    // ---- calibration -------------------------------------------------------------------------
    std::printf("\n--- calibration (%zu linear coefficients, geometry held) ---\n", kN);
    double r = fit_linear(train, geo, coeff, ridge);
    std::printf("  weighted rms residual %.4f (weighted units)\n", r);
    print_stats("fit, training set", evaluate(train, geo, coeff));
    print_stats("fit, holdout set", evaluate(hold, geo, coeff));
    {
        double biggest = 0.0;
        for (std::size_t m = 0; m < kM; ++m)
            for (std::size_t k = 0; k < kF; ++k) biggest = std::max(biggest, std::fabs(coeff[m][k]));
        const auto& st = cheatah_space_test::regime_drivers[2];
        double big_amp = 0.0;
        for (double a : (void)biggest, [&] {
                 std::array<double, kM> amps{};
                 const std::array<double, kF> f = t01s_features({st.dst, st.pdyn, st.by_imf, st.bz_imf, st.g2, st.g3});
                 for (std::size_t m = 0; m < kM; ++m) for (std::size_t k = 0; k < kF; ++k) amps[m] += coeff[m][k] * f[k];
                 return amps;
             }())
            big_amp = std::max(big_amp, std::fabs(a));
        std::printf("  ridge %g: largest |coefficient| %.3g, largest storm-regime |amplitude| %.3g\n", ridge, biggest, big_amp);
    }

    if (do_refine) {
        std::printf("\n--- geometry refinement (Nelder-Mead, linear fit inside) ---\n");
        geo = refine(train, geo, ridge, 300);
        r = fit_linear(train, geo, coeff, ridge);
        std::printf("  weighted rms residual %.4f\n", r);
        print_stats("refined, training set", evaluate(train, geo, coeff));
        print_stats("refined, holdout set", evaluate(hold, geo, coeff));
    }

    // ---- per-regime deviation of the calibrated model -------------------------------------------
    std::printf("\n--- per-regime deviation, calibrated model, holdout epochs ---\n");
    using cheatah_space_test::regime_drivers;
    using cheatah_space_test::storm_events;
    const std::array<const char*, 4> regime_names{"quiet", "moderate", "storm", "extreme"};
    const auto hold_pts = point_set(true, 7);
    for (int region = 0; region < 2; ++region) {
        const double r_hi = region == 0 ? 10.0 : 40.0;
        std::printf("%s\n", region == 0 ? " belts, 2.2 <= r <= 10 Re:" : " full box (x >= -15, inside the magnetopause):");
        for (std::size_t i = 0; i < 4; ++i) {
            const auto& m = regime_drivers[i];
            const std::vector<Sample> s = sample_oracle(hold_epochs, hold_pts, {{m.dst, m.pdyn, m.by_imf, m.bz_imf, m.g2, m.g3}}, 2.2, r_hi);
            print_stats(regime_names[i], evaluate(s, geo, coeff));
        }
        for (const auto& e : storm_events) {
            const std::vector<Sample> s = sample_oracle(hold_epochs, hold_pts, {{e.mag.dst, e.mag.pdyn, e.mag.by_imf, e.mag.bz_imf, e.mag.g2, e.mag.g3}}, 2.2, r_hi);
            print_stats(std::string(e.name).c_str(), evaluate(s, geo, coeff));
        }
    }

    // ---- spatial span per driver state: can the modes fit ONE state at all? ----------------
    std::printf("\n--- spatial span per state (%zu amplitudes free, belts r<=10, train epochs) ---\n", kM);
    {
        const auto pts = point_set(false, 1);
        std::vector<T01sDrivers> states;
        for (const auto& m : regime_drivers) states.push_back({m.dst, m.pdyn, m.by_imf, m.bz_imf, m.g2, m.g3});
        for (const auto& e : storm_events) states.push_back({e.mag.dst, e.mag.pdyn, e.mag.by_imf, e.mag.bz_imf, e.mag.g2, e.mag.g3});
        for (std::size_t si = 0; si < states.size(); ++si) {
            const std::vector<Sample> s = sample_oracle(train_epochs, pts, {states[si]}, 2.2, 10.0);
            // kM-unknown least squares
            std::vector<double> ata(kM * kM, 0.0), atb(kM, 0.0);
            std::array<std::array<double, 3>, kM> modes{};
            for (const Sample& q : s) {
                t01s_mode_fields<double>(geo, std::sin(q.tilt), std::cos(q.tilt), q.x, q.y, q.z, modes);
                for (int c = 0; c < 3; ++c)
                    for (std::size_t i = 0; i < kM; ++i) {
                        atb[i] += modes[i][static_cast<std::size_t>(c)] * q.b[static_cast<std::size_t>(c)];
                        for (std::size_t j = 0; j < kM; ++j) ata[(i * kM) + j] += modes[i][static_cast<std::size_t>(c)] * modes[j][static_cast<std::size_t>(c)];
                    }
            }
            std::vector<double> a = ata, b = atb;
            for (std::size_t i = 0; i < kM; ++i) a[(i * kM) + i] *= 1.0 + 1e-10;
            solve(a, b, kM);
            double res = 0, sig = 0;
            double bres[3] = {0, 0, 0}, bsig[3] = {0, 0, 0};
            for (const Sample& q : s) {
                t01s_mode_fields<double>(geo, std::sin(q.tilt), std::cos(q.tilt), q.x, q.y, q.z, modes);
                const double lat = std::fabs(std::atan2(q.z, std::sqrt(q.x * q.x + q.y * q.y))) * 180.0 / std::numbers::pi;
                const int band = lat < 15.0 ? 0 : lat < 35.0 ? 1 : 2;
                for (int c = 0; c < 3; ++c) {
                    double v = 0; for (std::size_t i = 0; i < kM; ++i) v += b[i] * modes[i][static_cast<std::size_t>(c)];
                    const double e2 = (v - q.b[static_cast<std::size_t>(c)]) * (v - q.b[static_cast<std::size_t>(c)]);
                    res += e2; bres[band] += e2;
                    sig += q.b[static_cast<std::size_t>(c)] * q.b[static_cast<std::size_t>(c)];
                    bsig[band] += q.b[static_cast<std::size_t>(c)] * q.b[static_cast<std::size_t>(c)];
                }
            }
            std::printf("    by |lat| band: <15: %.3f  15-35: %.3f  >35: %.3f\n", std::sqrt(bres[0] / bsig[0]), std::sqrt(bres[1] / bsig[1]), std::sqrt(bres[2] / bsig[2]));
            std::printf("  state %zu (dst %5.0f pdyn %5.1f g2 %5.1f g3 %5.1f): n %5zu  rms %8.3f nT  rel %.3f  amps:", si, states[si].dst, states[si].pdyn, states[si].g2, states[si].g3, s.size(), std::sqrt(res / (3.0 * s.size())), std::sqrt(res / sig));
            for (std::size_t i = 0; i < kM; ++i) std::printf(" %.3g", b[i]);
            std::printf("\n");
        }
    }

    if (do_emit) emit(geo, coeff);
    return 0;
}
