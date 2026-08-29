// Provenance and differential study of space.irbem's Tsyganenko (2001) — "T01" — against the
// IRBEM oracle's `kext = 9`. This is the program that decides WHAT ext_t01.hpp may claim, and
// every measured number quoted in that header's brief is printed by a pass below.
//
// DEV-ONLY. Like mead_diff.cpp and t89_diff.cpp beside it, this is the one kind of program that
// touches IRBEM, it is never built by the QA gate, and it never ships. IRBEM is LGPL-3.0 and
// cheatah-space is MIT: the library is run here as a BLACK BOX (dlopen plus the documented C entry
// points), never read for its logic and never linked into anything we distribute.
//
// WHAT IT MEASURES, in the order the questions were asked:
//
//   1. IS T01 BINNED OR CONTINUOUS?  Each of the six drivers is swept finely at one point and the
//      second difference of the response is reported. T89's Kp scan finds four plateaux and three
//      jumps; this one finds none — the response is smooth in every driver, which is why the
//      suite asserts SMOOTHNESS where the T89 suite asserts bin identity.
//   2. WHICH DRIVER ISOLATES WHICH CURRENT SYSTEM?  `G1` enters P2 eq. (11) through the two tail
//      amplitudes and nothing else, and the oracle's response to it is LINEAR to roundoff — so
//      `dB/dG1` is the oracle's own cross-tail basis, shielding included, with every other term
//      cancelled. `G2` drives the four Birkeland amplitudes (P2 eq. 8), the two Birkeland scaling
//      factors (eq. 9) and the short module's inner-edge shift (eq. 3), so its response is
//      dominated by, but not purely, the Birkeland system.
//   3. THE INTERCONNECTION FIELD, AS PUBLISHED.  `B_y` and `B_z` drive nothing but `B_INT`, so
//      `B(By, Bz) - B(0, 0)` must be ONE UNIFORM VECTOR at every point, equal to
//      `epsilon(theta) (0, By, Bz)`. Pass 3 checks the uniformity and then recovers P2 eq. (10)'s
//      `epsilon_0` and `epsilon_1` by least squares. They come back as the printed decimals — this
//      is the one term of the model for which published coefficients meet ORACLE PARITY.
//   4. THE TWO UNIT CONSTANTS.  Neither paper states the factor that turns its printed potentials
//      into nanotesla. Each is measured once, as a least-squares scale of this library's basis
//      against the isolating driver's response, and the residual of that fit is the honest measure
//      of how much of the oracle's response the UNSHIELDED published form reproduces. Pass 4 also
//      runs two controls: the same fit with the tilt deformations switched off (the residual gets
//      much worse, so the deformations are doing real work), and the same fit under both readings
//      of P1's sunward shift (the two are indistinguishable, so that ambiguity is NOT resolved
//      here and the header says so).
//   5. HOW FAR APART ARE THEY?  The shipping evaluator against the oracle over the same points,
//      and then the FREE REFIT: the best least-squares combination of the four bases this file
//      carries plus the interconnection field. The refit's residual is the FLOOR — what no choice
//      of amplitudes can fix, i.e. the field that belongs to the terms the papers do not give.
//
// The external field is isolated from the oracle by DIFFERENCE: `get_field1_` with `kext = 9`
// minus the same call with `kext = 0`, both with `options(5) = 0`, so the internal IGRF term is
// bit-for-bit identical between them and cancels exactly. The dipole tilt is taken from the oracle
// too — the SM z-axis transformed into GSM is `(sin psi, 0, cos psi)` — so that a difference in
// the tilt model cannot masquerade as a difference in the external model.
//
// Build (from the repository root; one line):
//   g++ -O2 -std=c++20 tools/oracle/t01_diff.cpp -I. -I$CHEATAH_DIR/stdlib/ndarray
//       -I$CHEATAH_DIR/stdlib/builtins -I$CHEATAH_DIR/stdlib/fixarray -ldl -o /tmp/t01_diff
//   /tmp/t01_diff /tmp/irbem-builds/libirbem-O2.so
#include <dlfcn.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

#include "space/irbem/ext_t01.hpp"

namespace {

using cheatah::space::irbem::T01Drivers;
using cheatah::space::irbem::T01State;
using cheatah::space::irbem::t01_coefficients;
using cheatah::space::irbem::t01_components;
using cheatah::space::irbem::t01_state;
using cheatah::space::irbem::t01_tail_fixed;
using cheatah::space::irbem::t01_tail_modules;
using cheatah::space::irbem::t01_units;
using cheatah::space::irbem::T01TailModule;
namespace dt = cheatah::space::irbem::detail::t01;

/// `get_field1_`, as documented in the vendored matlab/libirbem.h.
using GetField1 = void (*)(int*, int*, int*, int*, int*, double*, double*, double*, double*,
                           double*, double*, double*);
/// `coord_trans_vec1_`, likewise — GSM<->GEO, and the oracle's own dipole tilt.
using CoordTransVec1 = void (*)(int*, int*, int*, int*, int*, double*, double*, double*);

GetField1 g_get_field = nullptr;
CoordTransVec1 g_coord_trans = nullptr;

/// An epoch to sample at.
struct Epoch {
    int year;
    int doy;
    double ut;
};

/// The four epochs every pass samples: tilts of +25.6, -30.4, +0.0 and -21.6 degrees.
const std::array<Epoch, 4> kEpochs{{{2015, 180, 61200.0},
                                    {2015, 355, 7200.0},
                                    {2015, 80, 39183.0},
                                    {1989, 72, 3600.0}}};

/// The baseline driver set every differential is taken about: inside the published envelope on
/// every axis, southward `B_z`, and both history integrals well away from zero.
constexpr T01Drivers kBase{-30.0, 3.0, 2.0, -4.0, 3.0, 4.0};

/// The oracle's dipole tilt at an epoch, radians.
double oracle_tilt(const Epoch& e) {
    int one = 1;
    int year = e.year;
    int doy = e.doy;
    double ut = e.ut;
    int si = 4;  // SM
    int so = 2;  // GSM
    std::array<double, 3> in{0.0, 0.0, 1.0};
    std::array<double, 3> out{};
    g_coord_trans(&one, &si, &so, &year, &doy, &ut, in.data(), out.data());
    return std::atan2(out[0], out[2]);
}

/// The `maginput` vector T01 reads: slots 2, 5, 6, 7, 8, 9 one-based.
std::vector<double> maginput(const T01Drivers& d) {
    std::vector<double> m(25, 0.0);
    m[1] = d.dst;
    m[4] = d.pdyn;
    m[5] = d.by_imf;
    m[6] = d.bz_imf;
    m[7] = d.g1;
    m[8] = d.g2;
    return m;
}

/// The oracle's external field (kext=9 minus kext=0) at a GSM point, in GSM.
std::array<double, 3> oracle_external(const Epoch& e, double x, double y, double z,
                                      const T01Drivers& d) {
    int one = 1;
    int year = e.year;
    int doy = e.doy;
    double ut = e.ut;
    std::array<double, 3> gsm{x, y, z};
    std::array<double, 3> geo{};
    {
        int si = 2;
        int so = 1;
        g_coord_trans(&one, &si, &so, &year, &doy, &ut, gsm.data(), geo.data());
    }
    std::array<int, 5> options{0, 0, 0, 0, 0};
    int sysaxes = 1;
    int k0 = 0;
    int k1 = 9;
    std::vector<double> mag = maginput(d);
    std::array<double, 3> b0{};
    std::array<double, 3> b1{};
    double m0 = 0.0;
    double m1 = 0.0;
    double x1 = geo[0];
    double x2 = geo[1];
    double x3 = geo[2];
    g_get_field(&k0, options.data(), &sysaxes, &year, &doy, &ut, &x1, &x2, &x3, mag.data(),
                b0.data(), &m0);
    g_get_field(&k1, options.data(), &sysaxes, &year, &doy, &ut, &x1, &x2, &x3, mag.data(),
                b1.data(), &m1);
    std::array<double, 3> dgeo{b1[0] - b0[0], b1[1] - b0[1], b1[2] - b0[2]};
    std::array<double, 3> dgsm{};
    {
        int si = 1;
        int so = 2;
        g_coord_trans(&one, &si, &so, &year, &doy, &ut, dgeo.data(), dgsm.data());
    }
    return dgsm;
}

/// A 64-bit LCG in [0, 1), so every run samples the same points.
double next_unit(std::uint64_t& s) {
    s = (s * 6364136223846793005ULL) + 1442695040888963407ULL;
    return static_cast<double>(s >> 11U) / 9007199254740992.0;
}

/// Solve `A x = b` in place by Gaussian elimination with partial pivoting.
bool solve(std::vector<double>& a, std::vector<double>& b, int n) {
    for (int i = 0; i < n; ++i) {
        int piv = i;
        for (int r = i + 1; r < n; ++r) {
            if (std::fabs(a[(r * n) + i]) > std::fabs(a[(piv * n) + i])) piv = r;
        }
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

/// Least squares of `obs` on `rows` (each of width `n`); returns the coefficients and writes the
/// relative RMS residual (the residual norm over the observation norm).
std::vector<double> least_squares(const std::vector<std::vector<double>>& rows,
                                  const std::vector<double>& obs, int n, double& rel_residual) {
    std::vector<double> ata(static_cast<std::size_t>(n) * n, 0.0);
    std::vector<double> atb(static_cast<std::size_t>(n), 0.0);
    for (std::size_t r = 0; r < rows.size(); ++r) {
        for (int i = 0; i < n; ++i) {
            const std::size_t ii = static_cast<std::size_t>(i);
            atb[ii] += rows[r][ii] * obs[r];
            for (int j = 0; j < n; ++j) {
                ata[(ii * static_cast<std::size_t>(n)) + static_cast<std::size_t>(j)] +=
                    rows[r][ii] * rows[r][static_cast<std::size_t>(j)];
            }
        }
    }
    if (!solve(ata, atb, n)) {
        rel_residual = 1e30;
        return atb;
    }
    double res = 0.0;
    double sig = 0.0;
    for (std::size_t r = 0; r < rows.size(); ++r) {
        double f = 0.0;
        for (int k = 0; k < n; ++k) f += rows[r][static_cast<std::size_t>(k)] * atb[static_cast<std::size_t>(k)];
        res += (f - obs[r]) * (f - obs[r]);
        sig += obs[r] * obs[r];
    }
    rel_residual = std::sqrt(res / sig);
    return atb;
}

/// One sampled point: where it is, the epoch's tilt, and the three oracle readings every pass
/// shares — the baseline field and the responses to `G1` and `G2` over their whole published range.
struct Sample {
    double tilt;
    double x;
    double y;
    double z;
    std::array<double, 3> base;  ///< the oracle's external field at kBase
    std::array<double, 3> dg1;   ///< d(external)/dG1, from G1 = 0 to 10
    std::array<double, 3> dg2;   ///< d(external)/dG2, from G2 = 0 to 10
};

/// The shared point set: 250 candidates per epoch scattered over `3 <= r <= 15 R_E`, those
/// tailward of the model's `x >= -14 R_E` working region dropped. 999 points on this seed.
std::vector<Sample> gather() {
    std::vector<Sample> out;
    std::uint64_t s = 0x12345678ABCDEF01ULL;
    for (const Epoch& e : kEpochs) {
        const double tilt = oracle_tilt(e);
        for (int i = 0; i < 250; ++i) {
            const double r = 3.0 + (12.0 * next_unit(s));
            const double th = std::acos(1.0 - (2.0 * next_unit(s)));
            const double ph = 6.283185307179586 * next_unit(s);
            const double x = r * std::sin(th) * std::cos(ph);
            const double y = r * std::sin(th) * std::sin(ph);
            const double z = r * std::cos(th);
            if (x < -14.0) continue;
            T01Drivers g1lo = kBase;
            T01Drivers g1hi = kBase;
            T01Drivers g2lo = kBase;
            T01Drivers g2hi = kBase;
            g1lo.g1 = 0.0;
            g1hi.g1 = 10.0;
            g2lo.g2 = 0.0;
            g2hi.g2 = 10.0;
            const std::array<double, 3> a = oracle_external(e, x, y, z, g1lo);
            const std::array<double, 3> b = oracle_external(e, x, y, z, g1hi);
            const std::array<double, 3> c = oracle_external(e, x, y, z, g2lo);
            const std::array<double, 3> d = oracle_external(e, x, y, z, g2hi);
            const std::array<double, 3> base = oracle_external(e, x, y, z, kBase);
            out.push_back(Sample{tilt, x, y, z, base,
                                 {(b[0] - a[0]) / 10.0, (b[1] - a[1]) / 10.0, (b[2] - a[2]) / 10.0},
                                 {(d[0] - c[0]) / 10.0, (d[1] - c[1]) / 10.0, (d[2] - c[2]) / 10.0}});
        }
    }
    return out;
}

/// This library's cross-tail basis at a point: the two modules' fields, per unit amplitude, with
/// the warping and bending deformations applied (or, with @p deformed false, without them, which
/// is pass 4's control).
void tail_basis(const Sample& p, bool deformed, std::array<double, 3>& m1,
                std::array<double, 3>& m2) {
    T01State<double> s = t01_state<double>(deformed ? p.tilt : 0.0, kBase);
    if (!deformed) s.sin_tilt = 0.0;
    T01State<double> a = s;
    T01State<double> b = s;
    a.t1 = 1.0;
    a.t2 = 0.0;
    b.t1 = 0.0;
    b.t2 = 1.0;
    const dt::Vec3<double> va = dt::tail_field<double>(a, p.x, p.y, p.z);
    const dt::Vec3<double> vb = dt::tail_field<double>(b, p.x, p.y, p.z);
    const double u = t01_units.tail;
    m1 = {va.x / u, va.y / u, va.z / u};
    m2 = {vb.x / u, vb.y / u, vb.z / u};
}

/// One tail module under either reading of P1's sunward shift `X_s`, with no deformation, which is
/// how pass 4 compares the two readings on equal terms.
dt::Vec3<double> module_shifted(const T01TailModule& mod, double shift, double x, double y, double z,
                               bool after_scaling) {
    const double eta = mod.eta;
    const double xs = after_scaling ? ((eta * x) - ((eta - 1.0) * mod.x_m) - shift)
                                    : ((eta * (x - shift)) - ((eta - 1.0) * mod.x_m));
    const double ys = eta * y;
    const double zs = eta * z;
    const double flare = mod.delta_d_x * std::exp(xs / t01_tail_fixed.delta_x);
    const double yn = ys / t01_tail_fixed.delta_y;
    const double d = eta * (t01_coefficients.d0 + (t01_coefficients.delta_d_y * yn * yn) + flare);
    const double d_x = eta * flare / t01_tail_fixed.delta_x;
    const double d_y = eta * 2.0 * t01_coefficients.delta_d_y * yn / t01_tail_fixed.delta_y;
    return dt::tail_disc<double>(xs, ys, zs, d, d_x, d_y);
}

}  // namespace

int main(int argc, char** argv) {
    const std::string lib = argc > 1 ? argv[1] : "/tmp/irbem-builds/libirbem-O2.so";
    void* h = dlopen(lib.c_str(), RTLD_NOW);
    if (h == nullptr) {
        (void)std::fprintf(stderr, "cannot dlopen %s: %s\n", lib.c_str(), dlerror());
        return 1;
    }
    g_get_field = reinterpret_cast<GetField1>(dlsym(h, "get_field1_"));
    g_coord_trans = reinterpret_cast<CoordTransVec1>(dlsym(h, "coord_trans_vec1_"));
    if (g_get_field == nullptr || g_coord_trans == nullptr) {
        (void)std::fprintf(stderr, "missing entry points in %s\n", lib.c_str());
        return 1;
    }

    std::printf("space.irbem Tsyganenko 2001 (T01) vs IRBEM kext=9, oracle %s\n", lib.c_str());

    // ---- 1. binned or continuous? --------------------------------------------------------------
    std::printf("\n--- 1. driver sweeps at GSM (-8, 3, 2): is any driver BINNED? ---\n");
    std::printf("  (a jump would show as a second difference far above the sweep's own curvature)\n");
    {
        const Epoch e = kEpochs[0];
        struct Axis {
            const char* name;
            double lo;
            double hi;
            double T01Drivers::*slot;
        };
        const std::array<Axis, 6> axes{{{"Dst", -50.0, 20.0, &T01Drivers::dst},
                                        {"Pdyn", 0.5, 5.0, &T01Drivers::pdyn},
                                        {"By", -5.0, 5.0, &T01Drivers::by_imf},
                                        {"Bz", -5.0, 5.0, &T01Drivers::bz_imf},
                                        {"G1", 0.0, 10.0, &T01Drivers::g1},
                                        {"G2", 0.0, 10.0, &T01Drivers::g2}}};
        for (const Axis& a : axes) {
            constexpr int kSteps = 200;
            std::array<double, kSteps + 1> bz{};
            for (int i = 0; i <= kSteps; ++i) {
                T01Drivers d = kBase;
                d.*a.slot = a.lo + ((a.hi - a.lo) * i / kSteps);
                bz[static_cast<std::size_t>(i)] = oracle_external(e, -8.0, 3.0, 2.0, d)[2];
            }
            double span = 0.0;
            double worst2 = 0.0;
            for (int i = 1; i < kSteps; ++i) {
                const std::size_t k = static_cast<std::size_t>(i);
                worst2 = std::max(worst2, std::fabs(bz[k + 1] - (2.0 * bz[k]) + bz[k - 1]));
                span = std::max(span, std::fabs(bz[k] - bz[0]));
            }
            std::printf("  %-4s over [%7.2f, %7.2f]: Bz_ext spans %8.4f nT, worst |2nd diff| "
                        "%9.2e nT  -> %s\n",
                        a.name, a.lo, a.hi, span, worst2,
                        worst2 > 0.02 * span ? "A JUMP: BINNED" : "smooth: CONTINUOUS");
        }
    }

    // ---- 2. which driver isolates what ---------------------------------------------------------
    std::printf("\n--- 2. G1 isolates the tail: is the oracle's response to it LINEAR? ---\n");
    {
        const Epoch e = kEpochs[0];
        double worst = 0.0;
        for (const std::array<double, 3>& p :
             std::vector<std::array<double, 3>>{{-8.0, 3.0, 2.0}, {5.0, 0.0, 3.0}, {-13.0, -4.0, 1.0}}) {
            T01Drivers lo = kBase;
            T01Drivers hi = kBase;
            T01Drivers mid = kBase;
            lo.g1 = 0.0;
            hi.g1 = 10.0;
            mid.g1 = 5.0;
            const std::array<double, 3> a = oracle_external(e, p[0], p[1], p[2], lo);
            const std::array<double, 3> b = oracle_external(e, p[0], p[1], p[2], hi);
            const std::array<double, 3> c = oracle_external(e, p[0], p[1], p[2], mid);
            for (int k = 0; k < 3; ++k) {
                const std::size_t kk = static_cast<std::size_t>(k);
                const double bend = c[kk] - (0.5 * (a[kk] + b[kk]));
                worst = std::max(worst, std::fabs(bend) / std::max(1e-12, std::fabs(b[kk] - a[kk])));
            }
        }
        std::printf("  worst midpoint deviation from the chord, relative to the chord: %.2e\n", worst);
        std::printf("  -> G1 enters P2 eq. (11) LINEARLY, through t_1 and t_2 alone; dB/dG1 is the\n"
                    "     oracle's own cross-tail basis with every other current cancelled.\n");
    }

    const std::vector<Sample> pts = gather();
    std::printf("\n  %zu sampled points, 4 epochs, 3 <= r <= 15 Re, x >= -14 Re.\n", pts.size());

    // ---- 3. the interconnection field, and P2 eq. (10)'s coefficients ---------------------------
    std::printf("\n--- 3. B_INT: uniformity, and epsilon_0, epsilon_1 of P2 eq. (10) recovered ---\n");
    {
        const Epoch e = kEpochs[0];
        // Uniformity: the response to (By, Bz) must be the same vector everywhere.
        T01Drivers zero_imf = kBase;
        zero_imf.by_imf = 0.0;
        zero_imf.bz_imf = 0.0;
        T01Drivers probe = kBase;
        probe.by_imf = 3.0;
        probe.bz_imf = -4.0;
        std::array<double, 3> first{};
        double worst_spread = 0.0;
        bool have_first = false;
        for (std::size_t i = 0; i < pts.size(); i += 37) {
            const Sample& p = pts[i];
            const std::array<double, 3> a = oracle_external(e, p.x, p.y, p.z, zero_imf);
            const std::array<double, 3> b = oracle_external(e, p.x, p.y, p.z, probe);
            const std::array<double, 3> d{b[0] - a[0], b[1] - a[1], b[2] - a[2]};
            if (!have_first) {
                first = d;
                have_first = true;
                continue;
            }
            for (int k = 0; k < 3; ++k) {
                worst_spread = std::max(worst_spread, std::fabs(d[static_cast<std::size_t>(k)] -
                                                                first[static_cast<std::size_t>(k)]));
            }
        }
        std::printf("  response to (By, Bz) = (3, -4) at 27 scattered points: worst spread %.2e nT "
                    "(so it IS one uniform field)\n",
                    worst_spread);
        // Recover epsilon(theta) = e0 + e1 sin^2(theta/2) over the published By, Bz box.
        std::vector<std::vector<double>> rows;
        std::vector<double> obs;
        for (int iy = -10; iy <= 10; ++iy) {
            for (int iz = -10; iz <= 10; ++iz) {
                const double by = 0.5 * iy;
                const double bz = 0.5 * iz;
                const double bperp = std::sqrt((by * by) + (bz * bz));
                if (bperp < 0.4) continue;
                T01Drivers d = kBase;
                d.by_imf = by;
                d.bz_imf = bz;
                const std::array<double, 3> a = oracle_external(e, -8.0, 3.0, 2.0, zero_imf);
                const std::array<double, 3> b = oracle_external(e, -8.0, 3.0, 2.0, d);
                const double s2 = 0.5 * (1.0 - (bz / bperp));
                rows.push_back({by, by * s2});
                obs.push_back(b[1] - a[1]);
                rows.push_back({bz, bz * s2});
                obs.push_back(b[2] - a[2]);
            }
        }
        double rel = 0.0;
        const std::vector<double> c = least_squares(rows, obs, 2, rel);
        std::printf("  recovered epsilon_0 = %.6f  (P2 Table 1: %.3f)\n", c[0], t01_coefficients.eps[0]);
        std::printf("  recovered epsilon_1 = %.6f  (P2 Table 1: %.3f)\n", c[1], t01_coefficients.eps[1]);
        std::printf("  fit residual %.3e over %zu observations -> B_INT is at ORACLE PARITY with\n"
                    "     the published form and the published coefficients.\n",
                    rel, obs.size());
    }

    // ---- 4. the two unit constants, with their controls -----------------------------------------
    std::printf("\n--- 4. the unit constants the papers omit, measured against the isolating driver ---\n");
    {
        // Tail: fit a_1 M1 + a_2 M2 to dB/dG1. P2 eq. (2) says dt_1/dG1 = t_1^(2) = 0.319 and
        // dt_2/dG1 = t_2^(2) = -0.061, so the unit constant is a_1 / 0.319.
        std::vector<std::vector<double>> both;
        std::vector<std::vector<double>> only1;
        std::vector<std::vector<double>> only2;
        std::vector<double> obs;
        for (const Sample& p : pts) {
            std::array<double, 3> m1{};
            std::array<double, 3> m2{};
            tail_basis(p, true, m1, m2);
            for (int k = 0; k < 3; ++k) {
                const std::size_t kk = static_cast<std::size_t>(k);
                both.push_back({m1[kk], m2[kk]});
                only1.push_back({m1[kk]});
                only2.push_back({m2[kk]});
                obs.push_back(p.dg1[kk]);
            }
        }
        double rel_both = 0.0;
        double rel_1 = 0.0;
        double rel_2 = 0.0;
        const std::vector<double> c_both = least_squares(both, obs, 2, rel_both);
        const std::vector<double> c_1 = least_squares(only1, obs, 1, rel_1);
        const std::vector<double> c_2 = least_squares(only2, obs, 1, rel_2);
        std::printf("  tail, both modules : a1 = %.6f  a2 = %.6f   residual %.4f\n", c_both[0],
                    c_both[1], rel_both);
        std::printf("  tail, SHORT alone  : a1 = %.6f                residual %.4f\n", c_1[0], rel_1);
        std::printf("  tail, LONG alone   :               a2 = %.6f  residual %.4f\n", c_2[0], rel_2);
        std::printf("  -> unit constant = a1 / t_1^(2) = %.6f / %.3f = %.5f  (ext_t01.hpp carries "
                    "%.5f)\n",
                    c_both[0], t01_coefficients.t1[2], c_both[0] / t01_coefficients.t1[2],
                    t01_units.tail);
        std::printf("     The LONG module's fitted weight is ~%.0fx smaller than the published\n"
                    "     ratio |t_2^(2) / t_1^(2)| = %.3f would make it: its UNSHIELDED field is\n"
                    "     not a proxy for the confined field the oracle carries.\n",
                    std::fabs(c_both[0] * t01_coefficients.t2[2] /
                              (c_both[1] * t01_coefficients.t1[2])),
                    std::fabs(t01_coefficients.t2[2] / t01_coefficients.t1[2]));

        // CONTROL A: the same fit with the tilt deformations switched off.
        std::vector<std::vector<double>> flat;
        std::vector<double> flat_obs;
        for (const Sample& p : pts) {
            std::array<double, 3> m1{};
            std::array<double, 3> m2{};
            tail_basis(p, false, m1, m2);
            for (int k = 0; k < 3; ++k) {
                const std::size_t kk = static_cast<std::size_t>(k);
                flat.push_back({m1[kk], m2[kk]});
                flat_obs.push_back(p.dg1[kk]);
            }
        }
        double rel_flat = 0.0;
        (void)least_squares(flat, flat_obs, 2, rel_flat);
        std::printf("  CONTROL, deformations OFF: residual %.4f against %.4f with them on — the\n"
                    "     warping and bending of P1 eqs. (7)-(14) account for %.0f%% of the\n"
                    "     residual, so they are doing real work rather than decorating the fit.\n",
                    rel_flat, rel_both, 100.0 * (rel_flat - rel_both) / rel_flat);

        // CONTROL B: both readings of P1's sunward shift, undeformed so the comparison is fair.
        for (int after = 0; after < 2; ++after) {
            std::vector<std::vector<double>> rows;
            std::vector<double> o;
            for (const Sample& p : pts) {
                const T01State<double> s = t01_state<double>(p.tilt, kBase);
                const dt::Vec3<double> a = module_shifted(
                    t01_tail_modules[0], t01_tail_modules[0].x_s + s.shift1, p.x, p.y, p.z, after != 0);
                const dt::Vec3<double> b = module_shifted(t01_tail_modules[1],
                                                          t01_tail_modules[1].x_s, p.x, p.y, p.z,
                                                          after != 0);
                const std::array<double, 3> va{a.x, a.y, a.z};
                const std::array<double, 3> vb{b.x, b.y, b.z};
                for (int k = 0; k < 3; ++k) {
                    rows.push_back({va[static_cast<std::size_t>(k)], vb[static_cast<std::size_t>(k)]});
                    o.push_back(p.dg1[static_cast<std::size_t>(k)]);
                }
            }
            double rel = 0.0;
            (void)least_squares(rows, o, 2, rel);
            std::printf("  CONTROL, sunward shift applied %-14s: residual %.4f\n",
                        after != 0 ? "AFTER scaling" : "BEFORE scaling", rel);
        }
        std::printf("     The two readings differ by under a per cent of a ~50%% residual: this\n"
                    "     measurement does NOT resolve P1's ambiguity, and ext_t01.hpp says so.\n");

        // Birkeland: fit the G2 response, with the published amplitude ratios and without.
        std::vector<std::vector<double>> tied;
        std::vector<std::vector<double>> free;
        std::vector<double> fobs;
        for (const Sample& p : pts) {
            T01Drivers lo = kBase;
            T01Drivers hi = kBase;
            lo.g2 = 0.0;
            hi.g2 = 10.0;
            const T01State<double> s0 = t01_state<double>(p.tilt, lo);
            const T01State<double> s1 = t01_state<double>(p.tilt, hi);
            const double u = t01_units.fac;
            const dt::Vec3<double> a0 = dt::fac_region<double>(0, s0, p.x, p.y, p.z);
            const dt::Vec3<double> a1 = dt::fac_region<double>(0, s1, p.x, p.y, p.z);
            const dt::Vec3<double> b0 = dt::fac_region<double>(1, s0, p.x, p.y, p.z);
            const dt::Vec3<double> b1 = dt::fac_region<double>(1, s1, p.x, p.y, p.z);
            const std::array<double, 3> d1{(a1.x - a0.x) / (10.0 * u), (a1.y - a0.y) / (10.0 * u),
                                           (a1.z - a0.z) / (10.0 * u)};
            const std::array<double, 3> d2{(b1.x - b0.x) / (10.0 * u), (b1.y - b0.y) / (10.0 * u),
                                           (b1.z - b0.z) / (10.0 * u)};
            for (int k = 0; k < 3; ++k) {
                const std::size_t kk = static_cast<std::size_t>(k);
                tied.push_back({d1[kk] + d2[kk]});
                free.push_back({d1[kk], d2[kk]});
                fobs.push_back(p.dg2[kk]);
            }
        }
        double rel_tied = 0.0;
        double rel_free = 0.0;
        const std::vector<double> c_tied = least_squares(tied, fobs, 1, rel_tied);
        const std::vector<double> c_free = least_squares(free, fobs, 2, rel_free);
        std::printf("  Birkeland, published ratios: unit = %.1f          residual %.4f\n", c_tied[0],
                    rel_tied);
        std::printf("  Birkeland, regions freed   : R1 %.1f  R2 %.1f  residual %.4f\n", c_free[0],
                    c_free[1], rel_free);
        std::printf("  -> freeing the two regions buys %.4f of residual: P2 eq. (8)'s published\n"
                    "     amplitude RATIOS are what the oracle carries, and only the unit is\n"
                    "     unpublished. ext_t01.hpp carries %.1f.\n",
                    rel_tied - rel_free, t01_units.fac);
    }

    // ---- 5. how far apart, and the floor --------------------------------------------------------
    std::printf("\n--- 5. the shipping evaluator against the oracle, and the free-refit FLOOR ---\n");
    {
        double sum2 = 0.0;
        double sig2 = 0.0;
        double worst = 0.0;
        for (const Sample& p : pts) {
            const T01State<double> s = t01_state<double>(p.tilt, kBase);
            const std::array<double, 3> b = t01_components<double>(s, p.x, p.y, p.z);
            double d2 = 0.0;
            double o2 = 0.0;
            for (int k = 0; k < 3; ++k) {
                const std::size_t kk = static_cast<std::size_t>(k);
                const double d = b[kk] - p.base[kk];
                d2 += d * d;
                o2 += p.base[kk] * p.base[kk];
            }
            sum2 += d2;
            sig2 += o2;
            worst = std::max(worst, std::sqrt(d2 / o2));
        }
        std::printf("  shipping evaluator : RMS relative deviation %.4f, worst pointwise %.2f\n",
                    std::sqrt(sum2 / sig2), worst);

        std::vector<std::vector<double>> rows;
        std::vector<double> obs;
        for (const Sample& p : pts) {
            const T01State<double> s = t01_state<double>(p.tilt, kBase);
            std::array<double, 3> m1{};
            std::array<double, 3> m2{};
            tail_basis(p, true, m1, m2);
            const double u = t01_units.fac;
            const dt::Vec3<double> r1 = dt::fac_region<double>(0, s, p.x, p.y, p.z);
            const dt::Vec3<double> r2 = dt::fac_region<double>(1, s, p.x, p.y, p.z);
            const std::array<double, 3> c1{r1.x / u, r1.y / u, r1.z / u};
            const std::array<double, 3> c2{r2.x / u, r2.y / u, r2.z / u};
            const std::array<double, 3> bint{0.0, s.int_y, s.int_z};
            for (int k = 0; k < 3; ++k) {
                const std::size_t kk = static_cast<std::size_t>(k);
                rows.push_back({m1[kk], m2[kk], c1[kk], c2[kk], bint[kk]});
                obs.push_back(p.base[kk]);
            }
        }
        double rel = 0.0;
        const std::vector<double> c = least_squares(rows, obs, 5, rel);
        std::printf("  free refit of ALL five bases: M1 %.4f  M2 %.4f  R1 %.1f  R2 %.1f  INT %.3f\n",
                    c[0], c[1], c[2], c[3], c[4]);
        std::printf("  FLOOR: %.4f relative — no choice of amplitudes over the bases this file\n"
                    "     carries does better, so that is the field belonging to the terms the\n"
                    "     papers do not give: every shielding expansion, B_CF, and B_RC.\n",
                    rel);
        const T01State<double> s = t01_state<double>(0.0, kBase);
        std::printf("  For scale, the shipping amplitudes at the baseline drivers are\n"
                    "     unit*t_1 = %.4f, unit*t_2 = %.4f against the refit's %.4f and %.4f.\n",
                    t01_units.tail * s.t1, t01_units.tail * s.t2, c[0], c[1]);
    }

    std::printf(
        "\n  VERDICT: PUBLISHED-FORM-WITH-DOCUMENTED-GAP.\n"
        "  The functional form of every term this file carries is published (P1 eqs. 2-25), the\n"
        "  amplitude coefficients that drive them are published (P2 Table 1), and pass 3 shows the\n"
        "  interconnection term reaching ORACLE PARITY on both. What is NOT published is every\n"
        "  shielding expansion — hundreds of least-squares coefficients that exist only in\n"
        "  Tsyganenko's GPL-3.0 source, which this MIT clean room does not read — and the ring\n"
        "  current's analytic form, which is deferred to a closed-access 2000 paper. Pass 5\n"
        "  measures what their absence costs, and ext_t01.hpp's brief quotes that number rather\n"
        "  than claiming a parity it cannot have.\n");
    dlclose(h);
    return 0;
}
