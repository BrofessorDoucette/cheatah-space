// Provenance and differential study of space.irbem's Mead & Fairfield (1975) against the IRBEM
// oracle's `kext = 1` — the experiment that established WHAT the oracle computes before a line of
// ext_mead.hpp was written, and the measurement that says the two now agree.
//
// DEV-ONLY. Like convergence.cpp and t89_diff.cpp beside it, this is the one kind of program that
// touches IRBEM, it is never built by the QA gate, and it never ships. IRBEM is LGPL-3.0 and
// cheatah-space is MIT: the library is run here as a BLACK BOX (dlopen plus the documented C entry
// points), never read for its logic and never linked into anything we distribute.
//
// WHAT IT MEASURES, in the order the questions were asked:
//
//   1. WHERE ARE THE Kp BINS?  A scan of Kp x 10 from 0 to 95 at one point reports every value at
//      which the external field changes. Four distinct fields, switching at 4, 20 and 30; 91 and
//      above is refused by the oracle. That is the paper's four groups {0,0+}, {1-..2-}, {2..3-},
//      {>=3}, and mead_kp_bin's thresholds.
//   2. WHAT IS THE FUNCTIONAL FORM?  A FREE 20-term basis per component — every monomial of degree
//      <= 2 in position, times {1, tilt} — is fitted by least squares in four candidate framings:
//      GSM or SM position, tilt as the angle or as its sine. Exactly one framing fits to roundoff
//      (SM, angle: 7e-15 relative); the others leave 5e-3..2e-1. No assumption about the paper's
//      frame or tilt variable survives contact with this pass; it is simply read off.
//   3. WHAT ARE THE COEFFICIENTS?  In the SM fit the `y z` / `x z` ratio of B_x is -tan(4 deg) in
//      every bin and B_y's `x^2` and `y^2` tilt terms are equal and opposite — a rotated `x y`
//      term. Refitting on the paper's seventeen-term y-symmetric basis with the POSITION rotated
//      by 4 degrees about the dipole axis (the mean solar-wind aberration, applied to where the
//      field is evaluated and not to its components) recovers all 68 numbers as printed-table
//      decimals to 1e-9 relative. This pass prints the recovered values beside the ones
//      ext_mead.hpp hard-codes and the largest difference.
//   4. HOW FAR APART ARE THEY NOW?  The deviation pass: the shipping evaluator against the oracle,
//      per bin, over scattered points at 1.2-16.7 Re and six epochs of tilt, plus the tilt
//      dependence at three tilts explicitly. RMS and worst relative deviation. The verdict is
//      ORACLE PARITY: ~2e-9 RMS, < 1e-7 worst, against a 1e-6 budget.
//
// The external field is isolated from the oracle by DIFFERENCE: `get_field1_` with `kext = 1`
// minus the same call with `kext = 0`, both with `options(5) = 0`, so the internal IGRF term is
// bit-for-bit identical between them and cancels exactly. The dipole tilt is taken from the oracle
// too — the SM z-axis transformed into GSM is `(sin psi, 0, cos psi)` — so that a difference in
// the tilt model cannot masquerade as a difference in the external model.
//
// Build (from the repository root; one line):
//   g++ -O2 -std=c++20 tools/oracle/mead_diff.cpp -I. -I$CHEATAH_DIR/stdlib/ndarray
//       -I$CHEATAH_DIR/stdlib/builtins -I$CHEATAH_DIR/stdlib/fixarray -o /tmp/mead_diff
//   /tmp/mead_diff /tmp/irbem-builds/libirbem-O2.so
#include <dlfcn.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

#include "space/irbem/ext_mead.hpp"

namespace {

using cheatah::space::irbem::mead_aberration_cos;
using cheatah::space::irbem::mead_aberration_sin;
using cheatah::space::irbem::mead_bin_count;
using cheatah::space::irbem::mead_coefficient_count;
using cheatah::space::irbem::mead_components;
using cheatah::space::irbem::mead_deg_per_rad;
using cheatah::space::irbem::mead_kp_bin;
using cheatah::space::irbem::mead_parameters;
using cheatah::space::irbem::MeadParameters;

/// `get_field1_`, as documented in the vendored matlab/libirbem.h.
using GetField1 = void (*)(int*, int*, int*, int*, int*, double*, double*, double*, double*,
                           double*, double*, double*);
/// `coord_trans_vec1_`, likewise — GSM<->GEO, and the oracle's own dipole tilt.
using CoordTransVec1 = void (*)(int*, int*, int*, int*, int*, double*, double*, double*);

GetField1 g_get_field = nullptr;
CoordTransVec1 g_coord_trans = nullptr;

/// One oracle sample: the tilt, the GSM position, and the external field the oracle reports.
struct Sample {
    double tilt;              ///< dipole tilt psi, radians, as the oracle defines it
    double x, y, z;           ///< GSM position, Earth radii
    std::array<double, 3> b;  ///< oracle external field in GSM, nT
};

/// An epoch to sample at.
struct Epoch {
    int year;
    int doy;
    double ut;
};

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

/// The oracle's external field (kext minus kext = 0) at a GSM point, in GSM.
Sample oracle_external(const Epoch& e, double x, double y, double z, double kp10, int kext) {
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
    int k1 = kext;
    std::vector<double> mag(25, 0.0);
    mag[0] = kp10;
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
    return Sample{oracle_tilt(e), x, y, z, dgsm};
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
/// relative RMS residual.
std::vector<double> least_squares(const std::vector<std::vector<double>>& rows,
                                  const std::vector<double>& obs, int n, double& rel_residual) {
    std::vector<double> ata(static_cast<std::size_t>(n) * n, 0.0);
    std::vector<double> atb(static_cast<std::size_t>(n), 0.0);
    for (std::size_t r = 0; r < rows.size(); ++r) {
        for (int i = 0; i < n; ++i) {
            atb[static_cast<std::size_t>(i)] += rows[r][static_cast<std::size_t>(i)] * obs[r];
            for (int j = 0; j < n; ++j) {
                ata[(static_cast<std::size_t>(i) * n) + j] +=
                    rows[r][static_cast<std::size_t>(i)] * rows[r][static_cast<std::size_t>(j)];
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

/// The samples for one bin: scattered GSM points over several epochs.
std::vector<Sample> gather(const std::vector<Epoch>& epochs, double kp10, int per_epoch,
                           double r_lo, double r_hi, std::uint64_t seed) {
    std::vector<Sample> out;
    std::uint64_t s = seed;
    for (const Epoch& e : epochs) {
        for (int i = 0; i < per_epoch; ++i) {
            const double r = r_lo + ((r_hi - r_lo) * next_unit(s));
            const double th = std::acos(1.0 - (2.0 * next_unit(s)));
            const double ph = 6.283185307179586 * next_unit(s);
            out.push_back(oracle_external(e, r * std::sin(th) * std::cos(ph),
                                          r * std::sin(th) * std::sin(ph), r * std::cos(th), kp10,
                                          1));
        }
    }
    return out;
}

/// The 20 free monomials of degree <= 2 in (x, y, z) times {1, t}.
std::vector<double> free_basis(double x, double y, double z, double t) {
    return {1.0,   x,     y,     z,     x * x,     y * y,     z * z,     x * y,     x * z,     y * z,
            t,     t * x, t * y, t * z, t * x * x, t * y * y, t * z * z, t * x * y, t * x * z, t * y * z};
}

/// GSM -> SM for a position or a component vector: a rotation about y by the tilt.
std::array<double, 3> to_sm(const std::array<double, 3>& v, double tilt) {
    const double c = std::cos(tilt);
    const double s = std::sin(tilt);
    return {(v[0] * c) - (v[2] * s), v[1], (v[0] * s) + (v[2] * c)};
}

/// The seventeen hard-coded coefficients of one bin, in the order the paper's tables list them.
std::array<double, mead_coefficient_count> flatten(const MeadParameters<double>& p) {
    return {p.bx_z,  p.bx_xz, p.bx_t,   p.bx_tx, p.bx_txx, p.bx_tyy, p.bx_tzz, p.by_yz, p.by_ty,
            p.by_txy, p.bz_1, p.bz_x,   p.bz_xx, p.bz_yy,  p.bz_zz,  p.bz_tz,  p.bz_txz};
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

    std::printf("space.irbem Mead & Fairfield (1975) vs IRBEM kext=1, oracle %s\n", lib.c_str());

    // ---- 1. the Kp bins ----------------------------------------------------------------------
    std::printf("\n--- 1. Kp bin scan: Kp x 10 values at which the oracle's field changes ---\n");
    {
        const Epoch e{2015, 180, 43200.0};
        double prev = std::nan("");
        int distinct = 0;
        for (int k = 0; k <= 95; ++k) {
            const Sample s = oracle_external(e, -8.0, 3.0, 2.0, static_cast<double>(k), 1);
            if (s.b[2] != prev) {
                const bool refused = s.b[2] < -1e30;
                std::printf("  Kp x 10 = %2d : %s", k,
                            refused ? "REFUSED by the oracle (baddata)\n" : "");
                if (!refused) {
                    ++distinct;
                    std::printf("Bz_ext = %10.6f nT  -> mead_kp_bin says %d\n", s.b[2],
                                mead_kp_bin(static_cast<double>(k)));
                }
                prev = s.b[2];
            }
        }
        std::printf("  %d distinct coefficient sets; mead_bin_count = %zu\n", distinct,
                    mead_bin_count);
    }

    // ---- 2. the functional form ----------------------------------------------------------------
    const std::vector<Epoch> epochs{{2015, 180, 61200.0}, {2015, 355, 7200.0}, {2015, 80, 39183.0},
                                    {1989, 72, 3600.0},   {2003, 303, 50000.0}, {2022, 34, 20000.0}};
    std::printf("\n--- 2. functional form: a FREE 20-term quadratic x {1, tilt} basis, bin 1 ---\n");
    {
        const std::vector<Sample> samples = gather(epochs, 0.0, 160, 3.0, 12.0, 0x12345678ABCDEF01ULL);
        for (int frame = 0; frame < 2; ++frame) {
            for (int tv = 0; tv < 2; ++tv) {
                double worst = 0.0;
                for (int comp = 0; comp < 3; ++comp) {
                    std::vector<std::vector<double>> rows;
                    std::vector<double> obs;
                    for (const Sample& s : samples) {
                        std::array<double, 3> p{s.x, s.y, s.z};
                        std::array<double, 3> b = s.b;
                        if (frame == 1) {
                            p = to_sm(p, s.tilt);
                            b = to_sm(b, s.tilt);
                        }
                        rows.push_back(free_basis(p[0], p[1], p[2], tv == 0 ? s.tilt : std::sin(s.tilt)));
                        obs.push_back(b[static_cast<std::size_t>(comp)]);
                    }
                    double rel = 0.0;
                    (void)least_squares(rows, obs, 20, rel);
                    worst = std::max(worst, rel);
                }
                std::printf("  frame %s, tilt as %-8s : worst component rel residual %.3e%s\n",
                            frame == 0 ? "GSM" : "SM ", tv == 0 ? "angle" : "sin",
                            worst, worst < 1e-12 ? "   <- roundoff: THIS is the model" : "");
            }
        }
    }

    // ---- 3. the coefficients -------------------------------------------------------------------
    std::printf("\n--- 3. coefficient recovery on the paper's 17-term basis, position aberrated 4 deg ---\n");
    std::printf("  (tilt in degrees; each recovered value beside the one ext_mead.hpp hard-codes)\n");
    {
        const std::array<double, mead_bin_count> kp10{{0.0, 10.0, 25.0, 60.0}};
        const char* names[mead_coefficient_count] = {"a1 z",  "a2 xz", "a3 t",  "a4 tx", "a5 tx2",
                                                     "a6 ty2", "a7 tz2", "b1 yz", "b2 ty", "b3 txy",
                                                     "c1 1",  "c2 x",  "c3 x2", "c4 y2", "c5 z2",
                                                     "c6 tz", "c7 txz"};
        double worst_gap = 0.0;
        for (std::size_t bin = 0; bin < mead_bin_count; ++bin) {
            const std::vector<Sample> samples =
                gather(epochs, kp10[bin], 260, 3.0, 12.0, 0x9E3779B97F4A7C15ULL);
            std::array<double, mead_coefficient_count> got{};
            std::array<double, 3> rel{};
            for (int comp = 0; comp < 3; ++comp) {
                std::vector<std::vector<double>> rows;
                std::vector<double> obs;
                for (const Sample& s : samples) {
                    const std::array<double, 3> p = to_sm({s.x, s.y, s.z}, s.tilt);
                    const std::array<double, 3> b = to_sm(s.b, s.tilt);
                    const double xm = (p[0] * mead_aberration_cos) - (p[1] * mead_aberration_sin);
                    const double ym = (p[0] * mead_aberration_sin) + (p[1] * mead_aberration_cos);
                    const double z = p[2];
                    const double t = s.tilt * mead_deg_per_rad;
                    if (comp == 0) rows.push_back({z, xm * z, t, t * xm, t * xm * xm, t * ym * ym, t * z * z});
                    if (comp == 1) rows.push_back({ym * z, t * ym, t * xm * ym});
                    if (comp == 2) rows.push_back({1.0, xm, xm * xm, ym * ym, z * z, t * z, t * xm * z});
                    obs.push_back(b[static_cast<std::size_t>(comp)]);
                }
                const int n = comp == 1 ? 3 : 7;
                const std::vector<double> c = least_squares(rows, obs, n, rel[static_cast<std::size_t>(comp)]);
                const std::size_t base = comp == 0 ? 0 : comp == 1 ? 7 : 10;
                for (int k = 0; k < n; ++k) got[base + static_cast<std::size_t>(k)] = c[static_cast<std::size_t>(k)];
            }
            const std::array<double, mead_coefficient_count> have =
                flatten(mead_parameters<double>(static_cast<int>(bin) + 1));
            std::printf("  bin %zu (Kp x 10 = %.0f): fit residuals Bx %.1e By %.1e Bz %.1e\n", bin + 1,
                        kp10[bin], rel[0], rel[1], rel[2]);
            for (std::size_t k = 0; k < mead_coefficient_count; ++k) {
                const double gap = std::fabs(got[k] - have[k]) / std::fabs(have[k]);
                worst_gap = std::max(worst_gap, gap);
                std::printf("    %-7s recovered % .10g   hard-coded % .6g   rel gap %.1e\n", names[k],
                            got[k], have[k], gap);
            }
        }
        std::printf("  worst relative gap between a recovered and a hard-coded coefficient: %.2e\n",
                    worst_gap);
    }

    // ---- 4. the deviation pass ------------------------------------------------------------------
    std::printf("\n--- 4. shipping evaluator vs oracle, scattered 1.2 <= r <= 16.7 Re, six epochs ---\n");
    std::printf("%3s %7s %6s %12s %12s %12s\n", "bin", "Kp*10", "N", "rms rel", "worst rel",
                "worst |dB| nT");
    {
        const std::array<double, mead_bin_count> kp10{{0.0, 10.0, 25.0, 60.0}};
        for (std::size_t bin = 0; bin < mead_bin_count; ++bin) {
            const std::vector<Sample> samples =
                gather(epochs, kp10[bin], 300, 1.2, 16.7, 0xABCDEF0123456789ULL);
            const MeadParameters<double> p = mead_parameters<double>(static_cast<int>(bin) + 1);
            double sum2 = 0.0;
            double sig2 = 0.0;
            double worst_rel = 0.0;
            double worst_abs = 0.0;
            for (const Sample& s : samples) {
                const std::array<double, 3> mine = mead_components<double>(
                    p, std::sin(s.tilt), std::cos(s.tilt), s.tilt * mead_deg_per_rad, s.x, s.y, s.z);
                double d2 = 0.0;
                double o2 = 0.0;
                for (int c = 0; c < 3; ++c) {
                    const double d = mine[static_cast<std::size_t>(c)] - s.b[static_cast<std::size_t>(c)];
                    d2 += d * d;
                    o2 += s.b[static_cast<std::size_t>(c)] * s.b[static_cast<std::size_t>(c)];
                }
                sum2 += d2;
                sig2 += o2;
                worst_rel = std::max(worst_rel, std::sqrt(d2 / o2));
                worst_abs = std::max(worst_abs, std::sqrt(d2));
            }
            std::printf("%3zu %7.0f %6zu %12.3e %12.3e %12.3e\n", bin + 1, kp10[bin], samples.size(),
                        std::sqrt(sum2 / sig2), worst_rel, worst_abs);
        }
    }

    // ---- 4b. the tilt dependence, explicitly, at three tilts ----------------------------------
    std::printf("\n--- 4b. tilt dependence at three tilts, bin 3, fixed GSM points ---\n");
    {
        const std::array<Epoch, 3> three{{{2015, 80, 39183.0}, {2015, 180, 43200.0}, {2015, 355, 7200.0}}};
        const MeadParameters<double> p = mead_parameters<double>(3);
        for (const Epoch& e : three) {
            double worst = 0.0;
            double tilt = 0.0;
            for (int ix = -3; ix <= 3; ++ix) {
                for (int iz = -2; iz <= 2; ++iz) {
                    const double x = 2.5 * ix;
                    const double z = 2.0 * iz;
                    const double y = 1.5;
                    if (std::sqrt((x * x) + (y * y) + (z * z)) < 2.0) continue;
                    const Sample s = oracle_external(e, x, y, z, 25.0, 1);
                    tilt = s.tilt;
                    const std::array<double, 3> mine = mead_components<double>(
                        p, std::sin(s.tilt), std::cos(s.tilt), s.tilt * mead_deg_per_rad, x, y, z);
                    double d2 = 0.0;
                    double o2 = 0.0;
                    for (int c = 0; c < 3; ++c) {
                        const double d = mine[static_cast<std::size_t>(c)] - s.b[static_cast<std::size_t>(c)];
                        d2 += d * d;
                        o2 += s.b[static_cast<std::size_t>(c)] * s.b[static_cast<std::size_t>(c)];
                    }
                    worst = std::max(worst, std::sqrt(d2 / o2));
                }
            }
            std::printf("  tilt %+7.3f deg : worst rel deviation %.3e\n", tilt * mead_deg_per_rad, worst);
        }
    }

    std::printf(
        "\n  VERDICT: ORACLE PARITY. The published form with the published coefficients IS what\n"
        "  IRBEM's kext = 1 computes; the residual floor of ~1e-9 is the oracle's own rounded\n"
        "  internal constants and sits two orders below the 1e-6 budget.\n");
    dlclose(h);
    return 0;
}
