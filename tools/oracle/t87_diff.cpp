// Provenance study of IRBEM's `kext = 2` (Tsyganenko 1987 "short") and `kext = 3` (the same
// paper's "long" variant) — the experiment that had to run BEFORE a line of ext_t87.hpp could be
// written, and the measurement that says ext_t87.hpp must NOT be written yet.
//
// DEV-ONLY. Like mead_diff.cpp, t89_diff.cpp and convergence.cpp beside it, this is the one kind
// of program that touches IRBEM, it is never built by the QA gate, and it never ships. IRBEM is
// LGPL-3.0 and cheatah-space is MIT: the library is run here as a BLACK BOX (dlopen plus the
// documented C entry points), never read for its logic and never linked into anything we
// distribute.
//
// ============================================================================================
// THE VERDICT THIS PROGRAM ESTABLISHES
// ============================================================================================
//
// **space.irbem does not ship a T87.** Not "ships one with a documented gap" — does not ship one.
// Three independent findings, each measured below, each reproducible by re-running this file:
//
//  A. THE PUBLISHED MODEL IS UNAVAILABLE TO A CLEAN ROOM. Tsyganenko, *Global quantitative models
//     of the geomagnetic field in the cislunar magnetosphere for different disturbance levels*,
//     Planet. Space Sci. **35**(11):1347-1358 (1987) is the sole publication of BOTH the T87
//     functional form and its six Kp-binned coefficient tables, for both the short and the long
//     truncation. It is behind Elsevier's paywall (doi:10.1016/0032-0633(87)90046-8) and no
//     accessible secondary source reproduces its equations: not Tsyganenko's own 1990 review
//     (Space Sci. Rev. 54:75, also paywalled), not the Oulu space-physics reference that supplied
//     the T89c revision note in `ext_t89.hpp`, not SPENVIS, not the CCMC model pages — every one
//     of them names the model's three current systems in prose and prints no equation. The
//     clean-room rule forbids recovering the form from IRBEM's `source/*.f`, so there is no
//     published form to implement and no published table to transcribe. `ext_t89.hpp`'s response
//     to "form published, coefficients not" — implement the form, measure the coefficient gap —
//     has no analogue here, because the form is missing too.
//
//  B. BLACK-BOX RECOVERY DOES NOT REACH. `mead_diff.cpp` recovered Mead & Fairfield from the
//     oracle by least squares because that model is seventeen monomials. T87 is not: passes 5 and
//     6 below eliminate both families this library already owns. A FREE polynomial basis to degree
//     six — 168 terms per component, in GSM and SM, with the tilt entered as its angle and as its
//     sine — saturates at ~12% relative residual. The full published T89 form with ALL NINETEEN of
//     its linear amplitudes free and its nine non-linear parameters swept over the six published
//     Kp columns saturates at 22-28%. For contrast, that same experiment run against the model it
//     IS the family of reaches 4.0-5.7% (`ext_t89.hpp`'s free-refit floor). Two families
//     eliminated does not leave a third to guess at, and guessing one and printing Tsyganenko's
//     citation over it would be fabrication.
//
//  C. THE ORACLE ITSELF IS NOT A SOUND PARITY TARGET FOR EITHER VARIANT. Even granting a form,
//     there would be nothing well-defined to match:
//
//       - `kext = 3` IS NOT A FUNCTION OF ITS ARGUMENTS (pass 2). The same position, epoch and Kp
//         return two different fields depending on what was called before: the first evaluation
//         after the Kp bin changes differs from every later one by up to 18 nT in B_x at 8 R_E.
//         The transient is reproducible and bin-change-triggered — the signature of state that is
//         read before it is written — so a caller sweeping Kp gets one wrong sample per bin and a
//         caller tracing at fixed Kp gets one wrong sample per run.
//
//       - `kext = 2` BREAKS DAWN-DUSK SYMMETRY (pass 3). T87 reads no IMF B_y and no dawn-dusk
//         driver, so B_x must be even in y and B_y odd in y. B_y and B_z obey that to 3e-14 —
//         roundoff — and B_x does not, missing by up to 1.6 nT. `kext = 3` obeys all three to
//         roundoff at the same points, so this is not the harness: it is that variant.
//
//       - `kext = 2` HAS EIGHT Kp PLATEAUS WHERE THE PAPER PUBLISHES SIX (pass 1). Both variants
//         switch at Kp = 1-, 2-, 3-, 4-, 5- exactly, which is the paper's six groups and confirms
//         the harness reads `maginput(1)` correctly. `kext = 2` then switches AGAIN at Kp = 2o and
//         at Kp = 5+, twice, inside two of those groups.
//
//       - NEITHER VARIANT IS DIVERGENCE-FREE (pass 4). |div B| is flat under h-refinement from
//         1e-2 to 1e-3 R_E — so it is the model's, not the stencil's — at ~1e-4 nT/R_E for
//         `kext = 3` and ~1e-2 for `kext = 2`, a hundred times worse. A second-order stencil on a
//         genuinely solenoidal field falls as h^2 over that range, which is what
//         `IrbemT89.DivergenceVanishesEverywhere` measures and gets.
//
//     The three `kext = 2` anomalies point the same way and are consistent with a coefficient
//     table being indexed past its end, but this program does not diagnose IRBEM — it only
//     records that the target is not one a parity budget can be written against.
//
// So: no `space/irbem/ext_t87.hpp`, no `tests/irbem_t87_test.cpp`, no `t87_eval` in
// `gpu/irbem.slang`, and no CMake rows. `status.hpp` already carries the two envelope rows
// (`Tsyganenko1987Short`, `Tsyganenko1987Long`, r <= 30 and r <= 70 R_E) that the kext table
// publishes, and they stay: they are quoted from IRBEM's own documentation and are true whether
// or not this library evaluates the model. What is missing is the evaluator, and it is missing on
// purpose.
//
// ============================================================================================
// HOW THE EXTERNAL FIELD IS ISOLATED
// ============================================================================================
//
// Exactly as in mead_diff.cpp and t89_diff.cpp: `get_field1_` at the model's `kext` minus the same
// call at `kext = 0`, both with `options` all zero, so the internal IGRF term is bit-for-bit
// identical between them and cancels. The dipole tilt is taken from the oracle too — the SM z axis
// transformed into GSM is `(sin psi, 0, cos psi)` — so a difference in the tilt model cannot
// masquerade as a difference in the external model.
//
// ONE DEVIATION FROM THE OTHER STUDIES, AND IT IS FINDING C's FAULT: every sample is taken TWICE
// and the SECOND is kept. `kext = 3` returns a different answer on the first call after a Kp bin
// change (pass 2 measures exactly that), so a study that took one call per point would report the
// transient in some rows and the steady state in others, and could not tell the two apart. Warming
// makes `kext = 3` reproducible and does nothing at all to `kext = 2`, which is stateless.
//
// Build (from the repository root; one line):
//   g++ -O2 -std=c++20 tools/oracle/t87_diff.cpp -I. -I$CHEATAH_DIR/stdlib/ndarray
//       -I$CHEATAH_DIR/stdlib/builtins -I$CHEATAH_DIR/stdlib/fixarray -o /tmp/t87_diff -ldl
//   /tmp/t87_diff /tmp/irbem-builds/libirbem-O2.so
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

#include "space/irbem/ext_t89.hpp"

namespace {

using cheatah::space::irbem::t89_linear_count;
using cheatah::space::irbem::t89_parameters;
using cheatah::space::irbem::t89_published_set_count;
using cheatah::space::irbem::T89Parameters;
using cheatah::space::irbem::t89_components;

/// Degrees per radian, spelled once.
constexpr double kDegPerRad = 57.295779513082320876798154814105;

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

/// The six epochs every pass shares: four seasons of tilt plus the March 1989, Halloween 2003 and
/// February 2022 storm days, so the tilt spans -30 to +26 degrees rather than clustering.
const std::vector<Epoch>& epochs() {
    static const std::vector<Epoch> e{{2015, 180, 61200.0}, {2015, 355, 7200.0},
                                      {2015, 80, 39183.0},  {1989, 72, 3600.0},
                                      {2003, 303, 50000.0}, {2022, 34, 20000.0}};
    return e;
}

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
///
/// @param reps how many times to evaluate before keeping the answer. TWO is the study's default
///        and the reason is finding C: `kext = 3` answers differently on the first call after a Kp
///        bin change. Pass one to SEE that, which is what pass 2 does.
/// @return the external field, or all-NaN when the oracle refuses the point (its `baddata` fill).
std::array<double, 3> oracle_external(const Epoch& e, double x, double y, double z, double kp10,
                                      int kext, int reps = 2) {
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
    for (int r = 0; r < reps; ++r) {
        double m0 = 0.0;
        double m1 = 0.0;
        double x1 = geo[0];
        double x2 = geo[1];
        double x3 = geo[2];
        g_get_field(&k0, options.data(), &sysaxes, &year, &doy, &ut, &x1, &x2, &x3, mag.data(),
                    b0.data(), &m0);
        g_get_field(&k1, options.data(), &sysaxes, &year, &doy, &ut, &x1, &x2, &x3, mag.data(),
                    b1.data(), &m1);
    }
    const double nan = std::nan("");
    if (b1[0] < -1e30 || b0[0] < -1e30) return {nan, nan, nan};
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

/// Least squares of `obs` on `rows` (each of width `n`); returns the relative RMS residual and,
/// through @p out, the coefficients.
double least_squares(const std::vector<std::vector<double>>& rows, const std::vector<double>& obs,
                     int n, std::vector<double>& out) {
    const std::size_t un = static_cast<std::size_t>(n);
    std::vector<double> ata(un * un, 0.0);
    std::vector<double> atb(un, 0.0);
    for (std::size_t r = 0; r < rows.size(); ++r) {
        for (std::size_t i = 0; i < un; ++i) {
            atb[i] += rows[r][i] * obs[r];
            for (std::size_t j = 0; j < un; ++j) ata[(i * un) + j] += rows[r][i] * rows[r][j];
        }
    }
    if (!solve(ata, atb, n)) {
        out.assign(un, 0.0);
        return 1e30;
    }
    double res = 0.0;
    double sig = 0.0;
    for (std::size_t r = 0; r < rows.size(); ++r) {
        double f = 0.0;
        for (std::size_t k = 0; k < un; ++k) f += rows[r][k] * atb[k];
        res += (f - obs[r]) * (f - obs[r]);
        sig += obs[r] * obs[r];
    }
    out = atb;
    return std::sqrt(res / sig);
}

/// One oracle sample: the tilt, the GSM position, and the external field the oracle reports.
struct Sample {
    double tilt;
    double x, y, z;
    std::array<double, 3> b;
};

/// Scattered GSM points over the shared epochs, at one Kp and one `kext`.
std::vector<Sample> gather(double kp10, int kext, int per_epoch, double r_lo, double r_hi,
                           std::uint64_t seed) {
    std::vector<Sample> out;
    std::uint64_t s = seed;
    for (const Epoch& e : epochs()) {
        const double t = oracle_tilt(e);
        for (int i = 0; i < per_epoch; ++i) {
            const double r = r_lo + ((r_hi - r_lo) * next_unit(s));
            const double th = std::acos(1.0 - (2.0 * next_unit(s)));
            const double ph = 6.283185307179586 * next_unit(s);
            const double x = r * std::sin(th) * std::cos(ph);
            const double y = r * std::sin(th) * std::sin(ph);
            const double z = r * std::cos(th);
            const std::array<double, 3> b = oracle_external(e, x, y, z, kp10, kext);
            if (!std::isfinite(b[0])) continue;
            out.push_back(Sample{t, x, y, z, b});
        }
    }
    return out;
}

/// GSM -> SM for a position or a component vector: a rotation about y by the tilt.
std::array<double, 3> to_sm(const std::array<double, 3>& v, double tilt) {
    const double c = std::cos(tilt);
    const double s = std::sin(tilt);
    return {(v[0] * c) - (v[2] * s), v[1], (v[0] * s) + (v[2] * c)};
}

/// `|div B|` of the oracle's external field by a centred second-order stencil.
double oracle_div_b(const Epoch& e, double x, double y, double z, double kp10, int kext, double h) {
    const std::array<double, 3> px = oracle_external(e, x + h, y, z, kp10, kext);
    const std::array<double, 3> mx = oracle_external(e, x - h, y, z, kp10, kext);
    const std::array<double, 3> py = oracle_external(e, x, y + h, z, kp10, kext);
    const std::array<double, 3> my = oracle_external(e, x, y - h, z, kp10, kext);
    const std::array<double, 3> pz = oracle_external(e, x, y, z + h, kp10, kext);
    const std::array<double, 3> mz = oracle_external(e, x, y, z - h, kp10, kext);
    return std::fabs(((px[0] - mx[0]) + (py[1] - my[1]) + (pz[2] - mz[2])) / (2.0 * h));
}

/// The variants this study covers, and what IRBEM's kext table publishes about each.
struct Variant {
    int kext;
    const char* name;
    double max_r_geo;
};
constexpr std::array<Variant, 2> kVariants{{{2, "T87 short", 30.0}, {3, "T87 long", 70.0}}};

// -----------------------------------------------------------------------------------------
// pass 1 — where are the Kp bins?
// -----------------------------------------------------------------------------------------
void pass_kp_bins() {
    std::printf(
        "\n--- 1. Kp binning: every Kp x 10 at which the external field changes, step 0.001 ---\n");
    std::printf(
        "  The paper fits six Kp groups: {0,0+}, {1-,1,1+}, {2-,2,2+}, {3-,3,3+}, {4-,4,4+},\n"
        "  {>= 5-}. Their edges are Kp x 10 = 7, 17, 27, 37, 47 exactly -- the values 1-, 2-,\n"
        "  3-, 4-, 5- take in IRBEM's slot-1 scaling.\n");
    const Epoch e{2015, 180, 43200.0};
    for (const Variant& v : kVariants) {
        std::printf("  kext = %d (%s):\n", v.kext, v.name);
        double prev = oracle_external(e, -8.0, 3.0, 2.0, 0.0, v.kext)[2];
        int edges = 0;
        for (int i = 1; i <= 90000; ++i) {
            const double kp10 = static_cast<double>(i) * 0.001;
            const double bz = oracle_external(e, -8.0, 3.0, 2.0, kp10, v.kext)[2];
            if (bz != prev) {
                const bool published =
                    i == 7000 || i == 17000 || i == 27000 || i == 37000 || i == 47000;
                std::printf("    switches at Kp x 10 = %8.3f   %s\n", kp10,
                            published ? "(a published group edge)"
                                      : "<-- NOT A PUBLISHED GROUP EDGE");
                prev = bz;
                ++edges;
            }
        }
        std::printf("    %d edges, i.e. %d coefficient plateaus; the paper publishes 6\n", edges,
                    edges + 1);
    }
}

// -----------------------------------------------------------------------------------------
// pass 2 — is the oracle a function of its arguments?
// -----------------------------------------------------------------------------------------
void pass_statefulness() {
    std::printf("\n--- 2. is the oracle a FUNCTION of its arguments? (one call per sample) ---\n");
    std::printf(
        "  Identical position, epoch and Kp, evaluated three times in a row, then again after\n"
        "  the Kp bin has been changed and changed back. A pure function prints one row.\n");
    const Epoch e{2015, 180, 43200.0};
    for (const Variant& v : kVariants) {
        std::printf("  kext = %d (%s):\n", v.kext, v.name);
        const std::array<double, 6> sweep{20.0, 20.0, 20.0, 60.0, 20.0, 20.0};
        double worst = 0.0;
        std::array<double, 3> first{};
        for (std::size_t i = 0; i < sweep.size(); ++i) {
            const std::array<double, 3> b = oracle_external(e, -8.0, 3.0, 2.0, sweep[i], v.kext, 1);
            if (i == 0) first = b;
            const bool same_bin = sweep[i] == 20.0;
            if (same_bin) {
                for (int c = 0; c < 3; ++c) {
                    worst = std::max(worst, std::fabs(b[static_cast<std::size_t>(c)] -
                                                      first[static_cast<std::size_t>(c)]));
                }
            }
            std::printf("    call %zu at Kp x 10 = %4.0f : (%12.6f, %12.6f, %12.6f)%s\n", i + 1,
                        sweep[i], b[0], b[1], b[2], same_bin ? "" : "   <- bin changed");
        }
        std::printf("    worst spread over IDENTICAL arguments: %.6f nT%s\n", worst,
                    worst > 1e-9 ? "   <-- NOT A FUNCTION OF ITS ARGUMENTS" : "   (stateless)");
    }
}

// -----------------------------------------------------------------------------------------
// pass 3 — dawn-dusk parity
// -----------------------------------------------------------------------------------------
void pass_parity() {
    std::printf("\n--- 3. dawn-dusk parity: T87 reads no B_y and no dawn-dusk driver ---\n");
    std::printf(
        "  B_x and B_z must be EVEN in y and B_y ODD in y. Printed: the residual of each.\n");
    const Epoch e{2015, 180, 43200.0};
    const std::array<std::array<double, 3>, 4> pts{
        {{-8.0, 3.0, 2.0}, {5.0, 4.0, -1.0}, {-15.0, 6.0, 3.0}, {2.0, -5.0, 4.0}}};
    for (const Variant& v : kVariants) {
        double worst = 0.0;
        std::printf("  kext = %d (%s):\n", v.kext, v.name);
        for (const std::array<double, 3>& p : pts) {
            const std::array<double, 3> a = oracle_external(e, p[0], p[1], p[2], 20.0, v.kext);
            const std::array<double, 3> m = oracle_external(e, p[0], -p[1], p[2], 20.0, v.kext);
            const double dx = m[0] - a[0];
            const double dy = m[1] + a[1];
            const double dz = m[2] - a[2];
            worst = std::max(worst, std::max(std::fabs(dx), std::max(std::fabs(dy), std::fabs(dz))));
            std::printf("    (%6.1f,%5.1f,%5.1f) : dBx %11.3e  dBy %11.3e  dBz %11.3e\n", p[0],
                        p[1], p[2], dx, dy, dz);
        }
        std::printf("    worst parity residual %.3e nT%s\n", worst,
                    worst > 1e-9 ? "   <-- SYMMETRY BROKEN" : "   (roundoff: symmetry holds)");
    }
}

// -----------------------------------------------------------------------------------------
// pass 4 — is the oracle's field divergence-free?
// -----------------------------------------------------------------------------------------
void pass_divergence() {
    std::printf("\n--- 4. |div B| of the ORACLE's own external field, h = 1e-2 then 1e-3 R_E ---\n");
    std::printf(
        "  A second-order stencil on a solenoidal field falls as h^2 -- a hundredfold over this\n"
        "  refinement. A residual that does NOT move is the model's own, not the stencil's.\n");
    const Epoch e{2015, 180, 43200.0};
    const std::array<std::array<double, 3>, 3> pts{
        {{-8.0, 3.0, 2.0}, {5.0, 2.0, 1.0}, {-15.0, -4.0, 3.0}}};
    for (const Variant& v : kVariants) {
        std::printf("  kext = %d (%s):\n", v.kext, v.name);
        for (int b = 1; b <= 6; ++b) {
            const double kp10 = (b == 1) ? 0.0 : (static_cast<double>(b) * 10.0) - 3.0;
            std::printf("    Kp bin %d (Kp x 10 = %4.0f) :", b, kp10);
            for (const std::array<double, 3>& p : pts) {
                const double coarse = oracle_div_b(e, p[0], p[1], p[2], kp10, v.kext, 1e-2);
                const double fine = oracle_div_b(e, p[0], p[1], p[2], kp10, v.kext, 1e-3);
                std::printf("  %8.2e -> %8.2e", coarse, fine);
            }
            std::printf("\n");
        }
    }
}

// -----------------------------------------------------------------------------------------
// pass 5 — form elimination A: is it a polynomial?
// -----------------------------------------------------------------------------------------
void pass_polynomial() {
    std::printf("\n--- 5. form elimination A: a FREE polynomial basis x {1, tilt} ---\n");
    std::printf(
        "  This is the pass that recovered Mead & Fairfield outright (mead_diff.cpp step 2,\n"
        "  7e-15 relative in the right framing). Here it is run to degree six, 168 free terms\n"
        "  per component, in both frames and with the tilt entered both ways.\n");
    for (const Variant& v : kVariants) {
        const std::vector<Sample> samples = gather(20.0, v.kext, 400, 3.0, 12.0,
                                                   0x12345678ABCDEF01ULL);
        std::printf("  kext = %d (%s), Kp x 10 = 20, 3 <= r <= 12 R_E, N = %zu:\n", v.kext, v.name,
                    samples.size());
        for (int frame = 0; frame < 2; ++frame) {
            for (int tv = 0; tv < 2; ++tv) {
                for (int deg = 2; deg <= 6; deg += 2) {
                    double worst = 0.0;
                    int width = 0;
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
                            const double t = (tv == 0) ? s.tilt : std::sin(s.tilt);
                            std::vector<double> row;
                            for (int i = 0; i <= deg; ++i) {
                                for (int j = 0; i + j <= deg; ++j) {
                                    for (int k = 0; i + j + k <= deg; ++k) {
                                        const double m = std::pow(p[0], i) * std::pow(p[1], j) *
                                                         std::pow(p[2], k);
                                        row.push_back(m);
                                        row.push_back(t * m);
                                    }
                                }
                            }
                            width = static_cast<int>(row.size());
                            rows.push_back(std::move(row));
                            obs.push_back(b[static_cast<std::size_t>(comp)]);
                        }
                        std::vector<double> c;
                        worst = std::max(worst, least_squares(rows, obs, width, c));
                    }
                    std::printf(
                        "    frame %s  tilt as %-5s  degree %d (%3d terms) : worst rel resid "
                        "%.3e%s\n",
                        frame == 0 ? "GSM" : "SM ", tv == 0 ? "angle" : "sin", deg, width, worst,
                        worst < 1e-11 ? "   <- ROUNDOFF: this IS the model" : "");
                }
            }
        }
    }
}

// -----------------------------------------------------------------------------------------
// pass 6 — form elimination B: is it in the T89 family?
// -----------------------------------------------------------------------------------------
void pass_t89_family() {
    std::printf("\n--- 6. form elimination B: the PUBLISHED T89 form, all 19 amplitudes free ---\n");
    std::printf(
        "  T89 is T87's direct successor and shares its construction (disc-potential tail and\n"
        "  ring current, planar closure sheets, an exp(x/dx) polynomial). t89_components is\n"
        "  LINEAR in C_1..C_19, so the 19 field patterns are extracted with unit coefficient\n"
        "  vectors and their amplitudes solved exactly; the nine non-linear parameters are swept\n"
        "  over the six published Kp columns. The same experiment against the model this form\n"
        "  really is reaches 4.0-5.7%% (ext_t89.hpp's measured free-refit floor).\n");
    const int nc = static_cast<int>(t89_linear_count);
    for (const Variant& v : kVariants) {
        const double r_hi = (v.kext == 2) ? 12.0 : 25.0;
        const std::vector<Sample> samples =
            gather(20.0, v.kext, 300, 3.0, r_hi, 0x9E3779B97F4A7C15ULL);
        std::printf("  kext = %d (%s), Kp x 10 = 20, 3 <= r <= %.0f R_E, N = %zu:\n", v.kext,
                    v.name, r_hi, samples.size());
        double best = 1e30;
        for (std::size_t nb = 1; nb <= t89_published_set_count; ++nb) {
            const T89Parameters<double> base = t89_parameters<double>(static_cast<int>(nb));
            std::vector<std::vector<double>> rows;
            std::vector<double> obs;
            for (const Sample& s : samples) {
                std::array<std::array<double, 3>, t89_linear_count> pat{};
                for (std::size_t k = 0; k < t89_linear_count; ++k) {
                    T89Parameters<double> p = base;
                    for (std::size_t j = 0; j < t89_linear_count; ++j) {
                        p.c[j] = (j == k) ? 1.0 : 0.0;
                    }
                    pat[k] = t89_components<double>(p, std::sin(s.tilt), std::cos(s.tilt), s.x, s.y,
                                                    s.z);
                }
                for (std::size_t c = 0; c < 3; ++c) {
                    std::vector<double> row(t89_linear_count);
                    for (std::size_t k = 0; k < t89_linear_count; ++k) row[k] = pat[k][c];
                    rows.push_back(std::move(row));
                    obs.push_back(s.b[c]);
                }
            }
            std::vector<double> amp;
            const double rel = least_squares(rows, obs, nc, amp);
            best = std::min(best, rel);
            std::printf("    T89 non-linear column %zu : free 19-amplitude refit residual %.4e\n",
                        nb, rel);
        }
        std::printf("    best over all six columns: %.4e%s\n", best,
                    best > 0.10 ? "   <-- NOT THIS FAMILY" : "");
    }
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

    std::printf("space.irbem Tsyganenko (1987) provenance study vs IRBEM kext = 2 and 3\n");
    std::printf("oracle %s\n", lib.c_str());
    for (const Epoch& e : epochs()) {
        std::printf("  epoch %4d/%03d %8.0f s : oracle tilt %+8.4f deg\n", e.year, e.doy, e.ut,
                    oracle_tilt(e) * kDegPerRad);
    }

    pass_kp_bins();
    pass_statefulness();
    pass_parity();
    pass_divergence();
    pass_polynomial();
    pass_t89_family();

    std::printf(
        "\n  VERDICT: NEITHER VARIANT IS IMPLEMENTABLE AS A CLEAN ROOM, AND NEITHER IS A SOUND\n"
        "  PARITY TARGET. The 1987 paper is the sole publication of the functional form and the\n"
        "  Kp coefficient tables for both truncations and is paywalled; no accessible secondary\n"
        "  source reproduces its equations; the clean-room rule forbids reading IRBEM's Fortran.\n"
        "  Passes 5 and 6 eliminate the two families this library already owns, so black-box\n"
        "  recovery has nothing to fit. Independently, pass 2 shows kext = 3 is not a function\n"
        "  of its arguments, pass 3 shows kext = 2 breaks a symmetry the model's own driver set\n"
        "  requires, pass 1 shows kext = 2 carries eight Kp plateaus where six are published,\n"
        "  and pass 4 shows neither variant is divergence-free under h-refinement. space.irbem\n"
        "  therefore ships NO T87 evaluator; status.hpp keeps the two published envelope rows.\n");
    dlclose(h);
    return 0;
}
