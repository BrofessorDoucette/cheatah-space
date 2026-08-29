// Differential study of space.irbem's Tsyganenko (1996) against the IRBEM oracle's `kext = 7`:
// the provenance probe that decided the model's structure, the fit that produced its tables, and
// the report that measures what the published form cannot reach.
//
// DEV-ONLY. Like t89_diff.cpp beside it, this is the one kind of program that touches IRBEM, it is
// never built by the QA gate, and it never ships. IRBEM is LGPL-3.0 and cheatah-space is MIT: the
// library is run here as a BLACK BOX (dlopen plus the documented C entry points), never read for
// its logic and never linked into anything we distribute. Tsyganenko's own T96 source is GPL-3.0
// and is not read at all, by anything, ever.
//
// THREE PASSES, in the order the work was done:
//
//   probe   The measurements in ext_t96.hpp's file brief: that the oracle is exactly affine in
//           Dst, exactly homogeneous of degree one in the IMF, spanned by three clock-angle
//           functions, continuous in Pdyn and refused past 40 Re. These fix eq. (1) of the header.
//
//   fit     Sample the oracle's five driver-response fields at five tilts and three pressures,
//           refine the basis geometry by Nelder-Mead with the 3 840 linear amplitudes solved
//           exactly at every step (variable projection), and write the header's GENERATED block.
//           The basis is the header's own t96_basis: what ships is what was fitted.
//
//   report  With the header's CURRENT tables, on an INDEPENDENT sample (different scatter seed,
//           different pressures, an extra tilt): the per-family floor, the whole-field differential
//           over the four corpus regimes, and the dense southward-Bz storm sweep.
//
// The external field is isolated from the oracle by DIFFERENCE: `get_field1_` with `kext = 7`
// minus the same call with `kext = 0`, both with `options(5) = 0`, so the internal IGRF term is
// bit-for-bit identical and cancels exactly. The dipole tilt is taken from the oracle too.
//
// Build (from the repository root):
//   g++ -O3 -march=native -std=c++20 tools/oracle/t96_diff.cpp -I. \
//       -I$CHEATAH_DIR/stdlib/ndarray -I$CHEATAH_DIR/stdlib/builtins -I$CHEATAH_DIR/stdlib/fixarray \
//       -o /tmp/t96_diff
//   /tmp/t96_diff probe  [oracle.so]
//   /tmp/t96_diff fit    [oracle.so] [nm_iterations] [generated_block_out]
//   /tmp/t96_diff report [oracle.so]
//
// Quote the -O2 build, never the as-shipped one (docs/ERROR_BUDGET.md section 5).
#include <dlfcn.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

#include "space/irbem/ext_t96.hpp"

namespace {

namespace ir = cheatah::space::irbem;

using GetField1 = void (*)(int*, int*, int*, int*, int*, double*, double*, double*, double*,
                           double*, double*, double*);
using CoordTransVec1 = void (*)(int*, int*, int*, int*, int*, double*, double*, double*);

GetField1 g_get_field = nullptr;
CoordTransVec1 g_coord_trans = nullptr;
int g_year = 2015;
int g_doy = 180;
double g_ut = 43200.0;

/// The oracle's external field at a GSM point, in GSM nT: kext = 7 minus kext = 0.
std::array<double, 3> oracle_ext(double x, double y, double z, double dst, double pdyn, double by,
                                 double bz) {
    int one = 1;
    std::array<double, 3> gsm{x, y, z};
    std::array<double, 3> geo{};
    {
        int si = 2;
        int so = 1;
        g_coord_trans(&one, &si, &so, &g_year, &g_doy, &g_ut, gsm.data(), geo.data());
    }
    std::array<int, 5> options{0, 0, 0, 0, 0};
    int sysaxes = 1;
    int k0 = 0;
    int k7 = 7;
    std::vector<double> mag(25, 0.0);
    mag[1] = dst;
    mag[4] = pdyn;
    mag[5] = by;
    mag[6] = bz;
    std::array<double, 3> b0{};
    std::array<double, 3> b7{};
    double m0 = 0.0;
    double m7 = 0.0;
    double x1 = geo[0];
    double x2 = geo[1];
    double x3 = geo[2];
    g_get_field(&k0, options.data(), &sysaxes, &g_year, &g_doy, &g_ut, &x1, &x2, &x3, mag.data(),
                b0.data(), &m0);
    g_get_field(&k7, options.data(), &sysaxes, &g_year, &g_doy, &g_ut, &x1, &x2, &x3, mag.data(),
                b7.data(), &m7);
    std::array<double, 3> dgeo{b7[0] - b0[0], b7[1] - b0[1], b7[2] - b0[2]};
    std::array<double, 3> out{};
    {
        int si = 1;
        int so = 2;
        g_coord_trans(&one, &si, &so, &g_year, &g_doy, &g_ut, dgeo.data(), out.data());
    }
    return out;
}

/// The oracle's own dipole tilt for the current epoch: the SM z axis in GSM is (sin psi, 0, cos psi).
double oracle_tilt() {
    int one = 1;
    int si = 4;
    int so = 2;
    std::array<double, 3> in{0.0, 0.0, 1.0};
    std::array<double, 3> out{};
    g_coord_trans(&one, &si, &so, &g_year, &g_doy, &g_ut, in.data(), out.data());
    return std::atan2(out[0], out[2]);
}

bool is_baddata(const std::array<double, 3>& b) { return b[0] < -1e30 || b[1] < -1e30 || b[2] < -1e30; }

struct Pt {
    double x;
    double y;
    double z;
};

/// A deterministic scatter inside a magnetopause-shaped envelope: log-uniform in radius over
/// 1.2..35 Re, uniform in direction, kept when x <= 9 and y^2 + z^2 <= 180 - 15 x. The envelope is
/// where T96 means anything — outside the magnetopause the oracle's box harmonics grow without a
/// magnetosphere to shield, and fitting them there would fit nothing physical.
std::vector<Pt> make_points(std::size_t n, std::uint64_t seed) {
    std::vector<Pt> out;
    std::uint64_t s = seed;
    const auto next = [&s] {
        s = (s * 6364136223846793005ULL) + 1442695040888963407ULL;
        return static_cast<double>(s >> 11) / 9007199254740992.0;
    };
    while (out.size() < n) {
        const double r = 1.2 * std::pow(35.0 / 1.2, next());
        const double c = 1.0 - (2.0 * next());
        const double snt = std::sqrt(1.0 - (c * c));
        const double ph = 6.283185307179586 * next();
        const Pt p{r * snt * std::cos(ph), r * snt * std::sin(ph), r * c};
        if (p.x < -36.0 || p.x > 9.0) continue;
        if ((p.y * p.y) + (p.z * p.z) > 180.0 - (15.0 * p.x)) continue;
        out.push_back(p);
    }
    return out;
}

/// One sample: where, at what tilt and pressure, and the five driver-response fields there.
struct Sample {
    double tilt;
    double sqrt_pd;
    Pt p;
    std::array<std::array<double, 3>, ir::t96_family_count> tgt;  // B0, FD, FZ, FA, FY
};

/// Extract the five response fields at one point from six oracle calls — the decomposition of
/// eq. (1) in the header, using the measured facts: affine in Dst, `B_t F(theta)` in the IMF.
Sample sample_at(const Pt& p, double tilt, double pd) {
    const std::array<double, 3> b00 = oracle_ext(p.x, p.y, p.z, 0.0, pd, 0.0, 0.0);
    const std::array<double, 3> bd = oracle_ext(p.x, p.y, p.z, -100.0, pd, 0.0, 0.0);
    const std::array<double, 3> bn = oracle_ext(p.x, p.y, p.z, 0.0, pd, 0.0, 10.0);
    const std::array<double, 3> bs = oracle_ext(p.x, p.y, p.z, 0.0, pd, 0.0, -10.0);
    const std::array<double, 3> by = oracle_ext(p.x, p.y, p.z, 0.0, pd, 10.0, 0.0);
    Sample s{};
    s.tilt = tilt;
    s.sqrt_pd = std::sqrt(pd);
    s.p = p;
    for (int c = 0; c < 3; ++c) {
        s.tgt[0][c] = b00[c];
        s.tgt[1][c] = (bd[c] - b00[c]) / -100.0;
        s.tgt[2][c] = (bn[c] - b00[c]) / 10.0;                      // north: 10 F_Z
        s.tgt[3][c] = ((bs[c] - b00[c]) + (10.0 * s.tgt[2][c])) / 10.0;  // south: -10 F_Z + 10 F_A
        s.tgt[4][c] = ((by[c] - b00[c]) - ((10.0 / std::sqrt(2.0)) * s.tgt[3][c])) / 10.0;
    }
    return s;
}

struct Epoch {
    int doy;
    double ut;
};

std::vector<Sample> gather(const std::vector<Epoch>& epochs, const std::vector<double>& pds,
                           const std::vector<Pt>& pts) {
    std::vector<Sample> out;
    for (const Epoch& e : epochs) {
        g_doy = e.doy;
        g_ut = e.ut;
        const double tilt = oracle_tilt();
        std::printf("  epoch doy %3d ut %6.0f: tilt %+.4f rad\n", e.doy, e.ut, tilt);
        for (double pd : pds) {
            for (const Pt& p : pts) out.push_back(sample_at(p, tilt, pd));
        }
    }
    return out;
}

constexpr std::size_t kNB = ir::t96_basis_count;
constexpr std::size_t kNT = ir::t96_tilt_factor_count;
constexpr std::size_t kNP = ir::t96_pressure_factor_count;
constexpr std::size_t kNCOL = ir::t96_column_count;
constexpr std::size_t kNF = ir::t96_family_count;

struct Fit {
    std::array<std::vector<double>, kNF> coef;
    std::array<double, kNF> res{};
    std::array<double, kNF> sig{};
};

/// The joint linear least-squares solve: ONE normal matrix over the 768 tilt-and-pressure-expanded
/// columns (it is the same for every family) and five right-hand sides.
Fit joint_fit(const ir::T96Geometry& g, const std::vector<Sample>& S) {
    std::vector<double> ata(kNCOL * kNCOL, 0.0);
    std::vector<std::array<double, kNF>> atb(kNCOL, std::array<double, kNF>{});
    std::vector<double> row(kNCOL);
    ir::T96Basis<double> bas{};
    std::array<double, kNF> sig{};
    for (const Sample& s : S) {
        ir::t96_basis<double>(g, std::sin(s.tilt), std::cos(s.tilt), s.p.x, s.p.y, s.p.z, bas);
        const double sp = std::sin(s.tilt);
        const std::array<double, kNT> tf{1.0, sp, sp * sp};
        const std::array<double, kNP> pf{1.0, s.sqrt_pd};
        for (int c = 0; c < 3; ++c) {
            for (std::size_t j = 0; j < kNB; ++j)
                for (std::size_t t = 0; t < kNT; ++t)
                    for (std::size_t p = 0; p < kNP; ++p)
                        row[(((j * kNT) + t) * kNP) + p] = bas[j][static_cast<std::size_t>(c)] * tf[t] * pf[p];
            for (std::size_t i = 0; i < kNCOL; ++i) {
                const double ri = row[i];
                if (ri == 0.0) continue;
                for (std::size_t f = 0; f < kNF; ++f) atb[i][f] += ri * s.tgt[f][static_cast<std::size_t>(c)];
                double* arow = &ata[i * kNCOL];
                for (std::size_t k = i; k < kNCOL; ++k) arow[k] += ri * row[k];
            }
            for (std::size_t f = 0; f < kNF; ++f) sig[f] += s.tgt[f][static_cast<std::size_t>(c)] * s.tgt[f][static_cast<std::size_t>(c)];
        }
    }
    for (std::size_t i = 0; i < kNCOL; ++i)
        for (std::size_t k = 0; k < i; ++k) ata[(i * kNCOL) + k] = ata[(k * kNCOL) + i];
    double tr = 0.0;
    for (std::size_t i = 0; i < kNCOL; ++i) tr += ata[(i * kNCOL) + i];
    // A whisper of ridge: the box scales are a ladder and some columns are nearly collinear; this
    // keeps the Cholesky honest without moving anything the data constrains.
    for (std::size_t i = 0; i < kNCOL; ++i) ata[(i * kNCOL) + i] += 1e-9 * tr / static_cast<double>(kNCOL);
    std::vector<double> L(kNCOL * kNCOL, 0.0);
    for (std::size_t i = 0; i < kNCOL; ++i) {
        for (std::size_t k = 0; k <= i; ++k) {
            double sum = ata[(i * kNCOL) + k];
            for (std::size_t m = 0; m < k; ++m) sum -= L[(i * kNCOL) + m] * L[(k * kNCOL) + m];
            if (i == k) {
                L[(i * kNCOL) + i] = std::sqrt(sum > 0.0 ? sum : 1e-300);
            } else {
                L[(i * kNCOL) + k] = sum / L[(k * kNCOL) + k];
            }
        }
    }
    Fit out;
    for (std::size_t f = 0; f < kNF; ++f) {
        std::vector<double> yv(kNCOL);
        std::vector<double> xv(kNCOL);
        for (std::size_t i = 0; i < kNCOL; ++i) {
            double sum = atb[i][f];
            for (std::size_t m = 0; m < i; ++m) sum -= L[(i * kNCOL) + m] * yv[m];
            yv[i] = sum / L[(i * kNCOL) + i];
        }
        for (std::size_t ii = kNCOL; ii-- > 0;) {
            double sum = yv[ii];
            for (std::size_t m = ii + 1; m < kNCOL; ++m) sum -= L[(m * kNCOL) + ii] * xv[m];
            xv[ii] = sum / L[(ii * kNCOL) + ii];
        }
        out.coef[f] = xv;
        out.sig[f] = sig[f];
    }
    for (const Sample& s : S) {
        ir::t96_basis<double>(g, std::sin(s.tilt), std::cos(s.tilt), s.p.x, s.p.y, s.p.z, bas);
        const double sp = std::sin(s.tilt);
        const std::array<double, kNT> tf{1.0, sp, sp * sp};
        const std::array<double, kNP> pf{1.0, s.sqrt_pd};
        for (std::size_t f = 0; f < kNF; ++f) {
            for (std::size_t c = 0; c < 3; ++c) {
                double v = 0.0;
                for (std::size_t j = 0; j < kNB; ++j)
                    for (std::size_t t = 0; t < kNT; ++t)
                        for (std::size_t p = 0; p < kNP; ++p)
                            v += out.coef[f][(((j * kNT) + t) * kNP) + p] * bas[j][c] * tf[t] * pf[p];
                out.res[f] += (v - s.tgt[f][c]) * (v - s.tgt[f][c]);
            }
        }
    }
    return out;
}

const char* family_name(std::size_t f) {
    static const char* const names[kNF]{"B_0 (quiet)", "B_D (per nT Dst)", "B_Z (per nT Bz)",
                                        "B_A (per nT h)", "B_Y (per nT By)"};
    return names[f];
}

/// The geometry as a flat vector for the simplex, and back.
constexpr std::size_t kGeomDims = 5 + ir::t96_tail_truncation_count + ir::t96_tail_scale_count +
                                  ir::t96_ring_scale_count + ir::t96_fac_shell_count +
                                  ir::t96_fac_width_count + (2 * ir::t96_box_scale_count);

std::vector<double> pack(const ir::T96Geometry& g) {
    std::vector<double> v{g.r_c, g.g, g.d_0, g.gamma_rc, g.d_y};
    for (double d : g.x_0) v.push_back(d);
    for (double d : g.a_t) v.push_back(d);
    for (double d : g.a_rc) v.push_back(d);
    for (double d : g.fac_shell) v.push_back(d);
    for (double d : g.fac_width) v.push_back(d);
    for (double d : g.box_p) v.push_back(d);
    for (double d : g.box_r) v.push_back(d);
    return v;
}

ir::T96Geometry unpack(const std::vector<double>& v) {
    ir::T96Geometry g{};
    std::size_t i = 0;
    g.r_c = v[i++];
    g.g = v[i++];
    g.d_0 = v[i++];
    g.gamma_rc = v[i++];
    g.d_y = v[i++];
    for (double& d : g.x_0) d = v[i++];
    for (double& d : g.a_t) d = v[i++];
    for (double& d : g.a_rc) d = v[i++];
    for (double& d : g.fac_shell) d = v[i++];
    for (double& d : g.fac_width) d = v[i++];
    for (double& d : g.box_p) d = v[i++];
    for (double& d : g.box_r) d = v[i++];
    return g;
}

/// Whether a candidate geometry is one the basis can be evaluated at without nonsense: positive
/// thicknesses and scales, widths that are neither a step nor a smear.
bool plausible(const ir::T96Geometry& g) {
    if (g.d_0 < 0.3 || g.d_y < 3.0) return false;
    for (double d : g.a_t) if (d < 2.0) return false;
    for (double d : g.a_rc) if (d < 1.5) return false;
    for (double d : g.fac_shell) if (d < 1.5) return false;
    for (double d : g.fac_width) if (d < 0.03 || d > 0.8) return false;
    for (double d : g.box_p) if (d < 3.0) return false;
    for (double d : g.box_r) if (d < 2.0) return false;
    return true;
}

double cost_of(const ir::T96Geometry& g, const std::vector<Sample>& S) {
    if (!plausible(g)) return 1e30;
    const Fit f = joint_fit(g, S);
    double c = 0.0;
    for (std::size_t k = 0; k < kNF; ++k) c += f.res[k] / f.sig[k];
    return c;
}

/// Nelder-Mead over the geometry with the amplitudes solved exactly at every evaluation.
ir::T96Geometry refine(const ir::T96Geometry& start, const std::vector<Sample>& S, int iters) {
    const std::size_t nd = kGeomDims;
    std::vector<std::vector<double>> sx;
    std::vector<double> fv;
    sx.push_back(pack(start));
    fv.push_back(cost_of(start, S));
    for (std::size_t i = 0; i < nd; ++i) {
        std::vector<double> v = sx[0];
        v[i] += 0.15 * std::fabs(v[i]);
        sx.push_back(v);
        fv.push_back(cost_of(unpack(v), S));
    }
    for (int it = 0; it < iters; ++it) {
        std::vector<std::size_t> idx(nd + 1);
        for (std::size_t i = 0; i <= nd; ++i) idx[i] = i;
        std::sort(idx.begin(), idx.end(), [&](std::size_t a, std::size_t b) { return fv[a] < fv[b]; });
        std::vector<std::vector<double>> s2;
        std::vector<double> f2;
        for (std::size_t i : idx) {
            s2.push_back(sx[i]);
            f2.push_back(fv[i]);
        }
        sx = s2;
        fv = f2;
        if (it % 25 == 0) {
            std::printf("  NM iteration %4d: cost %.5f (worst vertex %.5f)\n", it, fv[0], fv[nd]);
            std::fflush(stdout);
        }
        if (std::fabs(fv[nd] - fv[0]) < 1e-5 * fv[0]) break;
        std::vector<double> cen(nd, 0.0);
        for (std::size_t i = 0; i < nd; ++i)
            for (std::size_t j = 0; j < nd; ++j) cen[j] += sx[i][j] / static_cast<double>(nd);
        const auto along = [&](double a) {
            std::vector<double> v(nd);
            for (std::size_t j = 0; j < nd; ++j) v[j] = cen[j] + (a * (cen[j] - sx[nd][j]));
            return v;
        };
        const std::vector<double> vr = along(1.0);
        const double fr = cost_of(unpack(vr), S);
        if (fr < fv[0]) {
            const std::vector<double> ve = along(2.0);
            const double fe = cost_of(unpack(ve), S);
            if (fe < fr) {
                sx[nd] = ve;
                fv[nd] = fe;
            } else {
                sx[nd] = vr;
                fv[nd] = fr;
            }
        } else if (fr < fv[nd - 1]) {
            sx[nd] = vr;
            fv[nd] = fr;
        } else {
            const std::vector<double> vc = along(-0.5);
            const double fc = cost_of(unpack(vc), S);
            if (fc < fv[nd]) {
                sx[nd] = vc;
                fv[nd] = fc;
            } else {
                for (std::size_t i = 1; i <= nd; ++i) {
                    for (std::size_t j = 0; j < nd; ++j) sx[i][j] = sx[0][j] + (0.5 * (sx[i][j] - sx[0][j]));
                    fv[i] = cost_of(unpack(sx[i]), S);
                }
            }
        }
    }
    return unpack(sx[0]);
}

void print_family_floors(const ir::T96Geometry& g, const Fit& f, const std::vector<Sample>& S) {
    const double edges[]{0.0, 2.0, 3.0, 5.0, 8.0, 12.0, 20.0, 40.0};
    ir::T96Basis<double> bas{};
    std::printf("\n  %-18s %10s %10s %8s   | per radius band, rel RMS: [1.2,2) [2,3) [3,5) [5,8) [8,12) [12,20) [20,35]\n",
                "family", "rms err", "rms signal", "rel");
    for (std::size_t fam = 0; fam < kNF; ++fam) {
        std::array<double, 7> rr{};
        std::array<double, 7> ss{};
        for (const Sample& s : S) {
            ir::t96_basis<double>(g, std::sin(s.tilt), std::cos(s.tilt), s.p.x, s.p.y, s.p.z, bas);
            const double sp = std::sin(s.tilt);
            const std::array<double, kNT> tf{1.0, sp, sp * sp};
            const std::array<double, kNP> pf{1.0, s.sqrt_pd};
            const double rad = std::sqrt((s.p.x * s.p.x) + (s.p.y * s.p.y) + (s.p.z * s.p.z));
            std::size_t e = 0;
            while (e < 6 && rad >= edges[e + 1]) ++e;
            for (std::size_t c = 0; c < 3; ++c) {
                double v = 0.0;
                for (std::size_t j = 0; j < kNB; ++j)
                    for (std::size_t t = 0; t < kNT; ++t)
                        for (std::size_t p = 0; p < kNP; ++p)
                            v += f.coef[fam][(((j * kNT) + t) * kNP) + p] * bas[j][c] * tf[t] * pf[p];
                rr[e] += (v - s.tgt[fam][c]) * (v - s.tgt[fam][c]);
                ss[e] += s.tgt[fam][c] * s.tgt[fam][c];
            }
        }
        std::printf("  %-18s %10.4f %10.4f %7.1f%%   |", family_name(fam),
                    std::sqrt(f.res[fam] / (3.0 * static_cast<double>(S.size()))),
                    std::sqrt(f.sig[fam] / (3.0 * static_cast<double>(S.size()))),
                    100.0 * std::sqrt(f.res[fam] / f.sig[fam]));
        for (std::size_t e = 0; e < 7; ++e) std::printf(" %5.1f%%", ss[e] > 0.0 ? 100.0 * std::sqrt(rr[e] / ss[e]) : 0.0);
        std::printf("\n");
    }
}

void write_generated(const ir::T96Geometry& g, const Fit& f, const std::string& path) {
    FILE* out = std::fopen(path.c_str(), "w");
    if (out == nullptr) {
        std::fprintf(stderr, "cannot write %s\n", path.c_str());
        return;
    }
    std::fprintf(out, "// ---- BEGIN GENERATED: written by tools/oracle/t96_diff.cpp `fit`; do not edit by hand ----------\n\n");
    std::fprintf(out, "/**\n * The frozen geometry — the starting ladders refined by the harness's Nelder-Mead pass.\n *\n * @test IrbemT96.GeometryIsTheFrozenFitGeometry\n */\n");
    std::fprintf(out, "inline constexpr T96Geometry t96_geometry{\n");
    std::fprintf(out, "    /* r_c */ %.9g, /* g */ %.9g, /* d_0 */ %.9g, /* gamma_rc */ %.9g, /* d_y */ %.9g,\n", g.r_c, g.g, g.d_0, g.gamma_rc, g.d_y);
    const auto list = [&](const char* name, const double* v, std::size_t n) {
        std::fprintf(out, "    /* %s */ {", name);
        for (std::size_t i = 0; i < n; ++i) std::fprintf(out, "%.9g%s", v[i], i + 1 < n ? ", " : "");
        std::fprintf(out, "}");
    };
    list("x_0", g.x_0.data(), g.x_0.size()); std::fprintf(out, ",\n");
    list("a_t", g.a_t.data(), g.a_t.size()); std::fprintf(out, ",\n");
    list("a_rc", g.a_rc.data(), g.a_rc.size()); std::fprintf(out, ",\n");
    list("fac_shell", g.fac_shell.data(), g.fac_shell.size()); std::fprintf(out, ",\n");
    list("fac_width", g.fac_width.data(), g.fac_width.size()); std::fprintf(out, ",\n");
    list("box_p", g.box_p.data(), g.box_p.size()); std::fprintf(out, ",\n");
    list("box_r", g.box_r.data(), g.box_r.size()); std::fprintf(out, "};\n\n");
    std::fprintf(out, "/**\n * The fitted amplitude tables: `t96_tables[family][basis * 6 + tilt_factor * 2 + pressure_factor]`.\n *\n"
                      " * Fitted to the IRBEM oracle's black-box output by the harness; see the file brief for what that\n"
                      " * means and what it does not. Families in @ref T96Family order.\n *\n * @test IrbemT96.TablesAreNotEmpty\n */\n");
    std::fprintf(out, "inline constexpr std::array<std::array<double, t96_column_count>, t96_family_count> t96_tables{{\n");
    for (std::size_t fam = 0; fam < kNF; ++fam) {
        std::fprintf(out, "    // %s\n    {{", family_name(fam));
        for (std::size_t i = 0; i < kNCOL; ++i) {
            if (i % 6 == 0) std::fprintf(out, "\n     ");
            std::fprintf(out, "%.9g%s", f.coef[fam][i], i + 1 < kNCOL ? ", " : "");
        }
        std::fprintf(out, "}},\n");
    }
    std::fprintf(out, "}};\n\n// ---- END GENERATED ---------------------------------------------------------------------------\n");
    std::fclose(out);
    std::printf("\n  wrote %s\n", path.c_str());
}

// ---------------------------------------------------------------------------------------------
// probe
// ---------------------------------------------------------------------------------------------

double rel_dev(const std::array<double, 3>& a, const std::array<double, 3>& b) {
    double e2 = 0.0;
    double m2 = 0.0;
    for (int c = 0; c < 3; ++c) {
        e2 += (a[c] - b[c]) * (a[c] - b[c]);
        m2 += b[c] * b[c];
    }
    return std::sqrt(e2 / (m2 + 1e-300));
}

void probe() {
    g_doy = 180;
    g_ut = 43200.0;
    std::printf("PROBE: the structure of kext = 7's driver dependence (tilt %+.4f rad)\n\n", oracle_tilt());
    const std::vector<Pt> pts{{-6.6, 0, 0}, {6, 0, 1}, {0, -6.6, 1}, {4, 4, -2}, {8, 2, 2}, {-10, 3, 2}, {-5, -5, 3}};

    double worst = 0.0;
    for (const Pt& p : pts) {
        const auto b0 = oracle_ext(p.x, p.y, p.z, 0, 3, 0, 0);
        const auto b1 = oracle_ext(p.x, p.y, p.z, -100, 3, 0, 0);
        const auto bt = oracle_ext(p.x, p.y, p.z, -63, 3, 0, 0);
        std::array<double, 3> pr{};
        for (int c = 0; c < 3; ++c) pr[c] = b0[c] + (0.63 * (b1[c] - b0[c]));
        worst = std::max(worst, rel_dev(pr, bt));
    }
    std::printf("  affine in Dst (predict -63 from 0 and -100):            worst rel dev %.2e\n", worst);

    worst = 0.0;
    for (const Pt& p : pts) {
        const auto b0 = oracle_ext(p.x, p.y, p.z, -20, 3, 0, 0);
        const auto b1 = oracle_ext(p.x, p.y, p.z, -20, 3, 8, -8);
        for (double l : {0.25, 0.5, 0.75}) {
            const auto bl = oracle_ext(p.x, p.y, p.z, -20, 3, 8 * l, -8 * l);
            std::array<double, 3> pr{};
            for (int c = 0; c < 3; ++c) pr[c] = b0[c] + (l * (b1[c] - b0[c]));
            worst = std::max(worst, rel_dev(pr, bl));
        }
    }
    std::printf("  homogeneous of degree 1 in (By, Bz) (rays at 1/4, 1/2, 3/4): worst rel dev %.2e\n", worst);

    worst = 0.0;
    for (const Pt& p : pts) {
        const auto b11 = oracle_ext(p.x, p.y, p.z, -100, 3, 0, -8);
        const auto b10 = oracle_ext(p.x, p.y, p.z, -100, 3, 0, 0);
        const auto b01 = oracle_ext(p.x, p.y, p.z, 0, 3, 0, -8);
        const auto b00 = oracle_ext(p.x, p.y, p.z, 0, 3, 0, 0);
        for (int c = 0; c < 3; ++c) worst = std::max(worst, std::fabs(b11[c] - b10[c] - b01[c] + b00[c]));
    }
    std::printf("  Dst x Bz cross term (Dst -100 x Bz -8):                 worst |r| %.2e nT\n", worst);

    double worst3 = 0.0;
    double worst_by_by = 0.0;
    for (const Pt& p : pts) {
        const auto b0 = oracle_ext(p.x, p.y, p.z, -20, 3, 0, 0);
        const int M = 24;
        std::vector<std::array<double, 3>> G(M);
        std::vector<double> th(M);
        for (int k = 0; k < M; ++k) {
            th[k] = 2.0 * std::numbers::pi * k / M;
            const auto b = oracle_ext(p.x, p.y, p.z, -20, 3, 10 * std::sin(th[k]), 10 * std::cos(th[k]));
            for (int c = 0; c < 3; ++c) G[k][c] = (b[c] - b0[c]) / 10.0;
        }
        for (int c = 0; c < 3; ++c) {
            // fit a sin + b cos + d sin(theta/2) by normal equations
            double a[9]{};
            double rhs[3]{};
            double sig = 0.0;
            for (int k = 0; k < M; ++k) {
                const double f[3]{std::sin(th[k]), std::cos(th[k]), std::sin(th[k] / 2.0)};
                for (int i = 0; i < 3; ++i) {
                    rhs[i] += f[i] * G[k][c];
                    for (int j = 0; j < 3; ++j) a[(i * 3) + j] += f[i] * f[j];
                }
                sig += G[k][c] * G[k][c];
            }
            // solve 3x3
            double m[9];
            std::memcpy(m, a, sizeof m);
            double x[3]{rhs[0], rhs[1], rhs[2]};
            for (int i = 0; i < 3; ++i) {
                for (int r = i + 1; r < 3; ++r) {
                    const double fct = m[(r * 3) + i] / m[(i * 3) + i];
                    for (int cc = i; cc < 3; ++cc) m[(r * 3) + cc] -= fct * m[(i * 3) + cc];
                    x[r] -= fct * x[i];
                }
            }
            for (int i = 2; i >= 0; --i) {
                double s = x[i];
                for (int cc = i + 1; cc < 3; ++cc) s -= m[(i * 3) + cc] * x[cc];
                x[i] = s / m[(i * 3) + i];
            }
            double res = 0.0;
            for (int k = 0; k < M; ++k) {
                const double v = (x[0] * std::sin(th[k])) + (x[1] * std::cos(th[k])) + (x[2] * std::sin(th[k] / 2.0));
                res += (v - G[k][c]) * (v - G[k][c]);
            }
            if (sig > 1e-12) worst3 = std::max(worst3, std::sqrt(res / sig));
        }
        // evenness in By: B(+By) + B(-By) - 2 B(0)
        const auto bp = oracle_ext(p.x, p.y, p.z, -20, 3, 10, 0);
        const auto bm = oracle_ext(p.x, p.y, p.z, -20, 3, -10, 0);
        for (int c = 0; c < 3; ++c) worst_by_by = std::max(worst_by_by, std::fabs(bp[c] + bm[c] - (2.0 * b0[c])));
    }
    std::printf("  clock angle spanned by {sin, cos, sin(theta/2)}:        worst rel RMS residual %.2e\n", worst3);
    std::printf("  By-even part B(+10) + B(-10) - 2 B(0):                   up to %.2f nT (the sin(theta/2) term)\n", worst_by_by);

    {
        const auto b1 = oracle_ext(-6.6, 0, 0, -20, 1.0, 0, -2);
        const auto b4 = oracle_ext(-6.6, 0, 0, -20, 4.0, 0, -2);
        for (double pd : {2.25, 9.0}) {
            const auto bt = oracle_ext(-6.6, 0, 0, -20, pd, 0, -2);
            const double s = (std::sqrt(pd) - 1.0) / (2.0 - 1.0);
            std::array<double, 3> pr{};
            for (int c = 0; c < 3; ++c) pr[c] = b1[c] + (s * (b4[c] - b1[c]));
            std::printf("  affine in sqrt(Pdyn) (predict %.2f nPa from 1 and 4):    rel dev %.2e\n", pd, rel_dev(pr, bt));
        }
        const auto a = oracle_ext(-6.6, 0, 0, -50, 2.999, 0, -5);
        const auto b = oracle_ext(-6.6, 0, 0, -50, 3.001, 0, -5);
        std::printf("  continuity in Pdyn across 3 nPa (step 0.002):             |dB| %.2e nT\n", rel_dev(a, b) * std::sqrt((b[0] * b[0]) + (b[1] * b[1]) + (b[2] * b[2])));
    }
    {
        const auto a = oracle_ext(-6.6, 0, 0, -20, 2, 0, -4.999);
        const auto b = oracle_ext(-6.6, 0, 0, -20, 2, 0, -5.001);
        std::printf("  continuity in Bz across -5 nT (step 0.002):              |dB| %.2e nT\n", std::sqrt(((a[0] - b[0]) * (a[0] - b[0])) + ((a[2] - b[2]) * (a[2] - b[2]))));
        const auto n = oracle_ext(-6.6, 0, 0, -20, 2, 0, 5);
        const auto z = oracle_ext(-6.6, 0, 0, -20, 2, 0, 0);
        const auto s = oracle_ext(-6.6, 0, 0, -20, 2, 0, -5);
        std::printf("  Bx slope per nT Bz at (-6.6,0,0): northward %+.4f, southward %+.4f — the kink\n", (n[0] - z[0]) / 5.0, (z[0] - s[0]) / 5.0);
    }
    for (double x : {-39.0, -41.0}) {
        const auto b = oracle_ext(x, 0, 0, -20, 3, 0, -5);
        std::printf("  at x = %.0f Re the oracle %s\n", x, is_baddata(b) ? "returns baddata (refused)" : "answers");
    }
}

// ---------------------------------------------------------------------------------------------
// fit
// ---------------------------------------------------------------------------------------------

void fit(int nm_iters, const std::string& out_path) {
    std::printf("FIT: sampling the oracle's five response fields\n");
    const std::vector<Epoch> epochs{{172, 61200.0}, {172, 43200.0}, {80, 39183.0}, {355, 43200.0}, {355, 20500.0}};
    const std::vector<double> pds{1.0, 3.0, 8.0};
    const std::vector<Pt> pts = make_points(360, 0x9E3779B97F4A7C15ULL);
    const std::vector<Sample> S = gather(epochs, pds, pts);
    std::printf("  %zu samples (%zu points x %zu tilts x %zu pressures)\n", S.size(), pts.size(), epochs.size(), pds.size());

    ir::T96Geometry g = ir::t96_geometry;
    std::printf("\n  starting geometry: cost %.5f\n", cost_of(g, S));
    if (nm_iters > 0) {
        // The simplex runs on a thinned sample so an evaluation costs ~1 s; the final solve uses all.
        std::vector<Sample> thin;
        for (std::size_t i = 0; i < S.size(); i += 3) thin.push_back(S[i]);
        std::printf("  Nelder-Mead over %zu geometric parameters on %zu samples, %d iterations\n", kGeomDims, thin.size(), nm_iters);
        g = refine(g, thin, nm_iters);
    }
    const Fit f = joint_fit(g, S);
    std::printf("\n  final geometry: cost %.5f\n", cost_of(g, S));
    std::printf("  geometry = {r_c %.4f, g %.4f, d_0 %.4f, gamma_rc %.4f, d_y %.4f, x_0 {%.3f %.3f}, a_t {%.3f %.3f %.3f %.3f},\n"
                "              a_rc {%.3f %.3f %.3f %.3f}, L {%.3f %.3f %.3f %.3f %.3f %.3f}, w {%.4f %.4f}, p {%.3f %.3f %.3f}, r {%.3f %.3f %.3f}}\n",
                g.r_c, g.g, g.d_0, g.gamma_rc, g.d_y, g.x_0[0], g.x_0[1], g.a_t[0], g.a_t[1], g.a_t[2], g.a_t[3],
                g.a_rc[0], g.a_rc[1], g.a_rc[2], g.a_rc[3], g.fac_shell[0], g.fac_shell[1], g.fac_shell[2], g.fac_shell[3], g.fac_shell[4], g.fac_shell[5],
                g.fac_width[0], g.fac_width[1], g.box_p[0], g.box_p[1], g.box_p[2], g.box_r[0], g.box_r[1], g.box_r[2]);
    print_family_floors(g, f, S);
    write_generated(g, f, out_path);
}

// ---------------------------------------------------------------------------------------------
// report
// ---------------------------------------------------------------------------------------------

/// The shipped model's field at a GSM point through the shipped entry point, GSM nT.
std::array<double, 3> mine(double tilt, const Pt& p, double dst, double pd, double by, double bz) {
    const ir::T96Amplitudes<double> a = ir::t96_amplitudes<double>(tilt, dst, pd, by, bz);
    return ir::t96_components<double>(a, p.x, p.y, p.z);
}

struct Stats {
    double sum2 = 0.0;
    double sig2 = 0.0;
    double worst = 0.0;
    std::size_t n = 0;
    void add(const std::array<double, 3>& m, const std::array<double, 3>& o) {
        double d2 = 0.0;
        double o2 = 0.0;
        for (int c = 0; c < 3; ++c) {
            d2 += (m[c] - o[c]) * (m[c] - o[c]);
            o2 += o[c] * o[c];
        }
        sum2 += d2;
        sig2 += o2;
        worst = std::max(worst, std::sqrt(d2));
        ++n;
    }
    double rms() const { return n > 0 ? std::sqrt(sum2 / static_cast<double>(n)) : 0.0; }
    double rel() const { return sig2 > 0.0 ? std::sqrt(sum2 / sig2) : 0.0; }
};

void report() {
    std::printf("REPORT: the shipped tables against the oracle on an INDEPENDENT sample\n");
    const std::vector<Epoch> epochs{{172, 61200.0}, {100, 30000.0}, {80, 39183.0}, {300, 50000.0}, {355, 20500.0}};
    const std::vector<double> pds{0.7, 2.0, 5.0, 9.5};
    const std::vector<Pt> pts = make_points(300, 0x2545F4914F6CDD1DULL);
    std::printf("  per-family floors (fitting decomposition, independent points/pressures/tilts):\n");
    const std::vector<Sample> S = gather(epochs, pds, pts);
    Fit f;
    for (std::size_t fam = 0; fam < kNF; ++fam) {
        f.coef[fam].assign(ir::t96_tables[fam].begin(), ir::t96_tables[fam].end());
        f.res[fam] = 0.0;
        f.sig[fam] = 0.0;
    }
    ir::T96Basis<double> bas{};
    for (const Sample& s : S) {
        ir::t96_basis<double>(ir::t96_geometry, std::sin(s.tilt), std::cos(s.tilt), s.p.x, s.p.y, s.p.z, bas);
        const double sp = std::sin(s.tilt);
        const std::array<double, kNT> tf{1.0, sp, sp * sp};
        const std::array<double, kNP> pf{1.0, s.sqrt_pd};
        for (std::size_t fam = 0; fam < kNF; ++fam)
            for (std::size_t c = 0; c < 3; ++c) {
                double v = 0.0;
                for (std::size_t j = 0; j < kNB; ++j)
                    for (std::size_t t = 0; t < kNT; ++t)
                        for (std::size_t p = 0; p < kNP; ++p)
                            v += f.coef[fam][(((j * kNT) + t) * kNP) + p] * bas[j][c] * tf[t] * pf[p];
                f.res[fam] += (v - s.tgt[fam][c]) * (v - s.tgt[fam][c]);
                f.sig[fam] += s.tgt[fam][c] * s.tgt[fam][c];
            }
    }
    print_family_floors(ir::t96_geometry, f, S);

    // Whole-field differentials over the corpus regimes, drivers clamped to the envelope.
    struct Regime {
        const char* name;
        double dst;
        double pd;
        double by;
        double bz;
    };
    const std::vector<Regime> regimes{
        {"quiet    (Dst -8,   Pdyn 1.8, By +1,  Bz +2)", -8.0, 1.8, 1.0, 2.0},
        {"moderate (Dst -42,  Pdyn 3.2, By -4,  Bz -5)", -42.0, 3.2, -4.0, -5.0},
        {"storm    (Dst -100, Pdyn 9.0, By +8,  Bz -10) [clamped]", -100.0, 9.0, 8.0, -10.0},
        {"extreme  (Dst -100, Pdyn 10,  By -10, Bz -10) [clamped]", -100.0, 10.0, -10.0, -10.0},
    };
    std::printf("\n  whole-field differential, vector RMS |dB| over %zu points x %zu tilts:\n", pts.size(), epochs.size());
    std::printf("  %-58s %9s %9s %9s | belts 3-8 Re: %9s %9s\n", "regime", "rms|dB|", "rms|B|", "rel", "rms|dB|", "rel");
    for (const Regime& rg : regimes) {
        Stats all;
        Stats belts;
        for (const Epoch& e : epochs) {
            g_doy = e.doy;
            g_ut = e.ut;
            const double tilt = oracle_tilt();
            for (const Pt& p : pts) {
                const auto o = oracle_ext(p.x, p.y, p.z, rg.dst, rg.pd, rg.by, rg.bz);
                if (is_baddata(o)) continue;
                const auto m = mine(tilt, p, rg.dst, rg.pd, rg.by, rg.bz);
                all.add(m, o);
                const double r = std::sqrt((p.x * p.x) + (p.y * p.y) + (p.z * p.z));
                if (r >= 3.0 && r <= 8.0) belts.add(m, o);
            }
        }
        std::printf("  %-58s %9.3f %9.3f %8.1f%% |               %9.3f %8.1f%%\n", rg.name, all.rms(),
                    std::sqrt(all.sig2 / static_cast<double>(all.n)), 100.0 * all.rel(), belts.rms(), 100.0 * belts.rel());
    }

    // The storm sweep: southward Bz densely, at the storm regime's other drivers, belts only.
    std::printf("\n  southward-Bz storm sweep at Dst -100, Pdyn 6, By +5, belts 3-8 Re, 3 tilts:\n");
    std::printf("  %6s %9s %9s %8s\n", "Bz", "rms|dB|", "rms|B|", "rel");
    for (double bz : {2.0, 0.0, -1.0, -2.0, -3.0, -4.0, -5.0, -6.0, -7.0, -8.0, -9.0, -10.0}) {
        Stats st;
        for (const Epoch& e : {epochs[0], epochs[2], epochs[4]}) {
            g_doy = e.doy;
            g_ut = e.ut;
            const double tilt = oracle_tilt();
            for (const Pt& p : pts) {
                const double r = std::sqrt((p.x * p.x) + (p.y * p.y) + (p.z * p.z));
                if (r < 3.0 || r > 8.0) continue;
                const auto o = oracle_ext(p.x, p.y, p.z, -100.0, 6.0, 5.0, bz);
                if (is_baddata(o)) continue;
                st.add(mine(tilt, p, -100.0, 6.0, 5.0, bz), o);
            }
        }
        std::printf("  %6.1f %9.3f %9.3f %7.1f%%\n", bz, st.rms(), std::sqrt(st.sig2 / static_cast<double>(st.n)), 100.0 * st.rel());
    }
}

}  // namespace

int main(int argc, char** argv) {
    const std::string pass = argc > 1 ? argv[1] : "report";
    const std::string lib = argc > 2 ? argv[2] : "/tmp/irbem-builds/libirbem-O2.so";
    void* h = dlopen(lib.c_str(), RTLD_NOW);
    if (h == nullptr) {
        std::fprintf(stderr, "cannot dlopen %s: %s\n", lib.c_str(), dlerror());
        return 1;
    }
    g_get_field = reinterpret_cast<GetField1>(dlsym(h, "get_field1_"));
    g_coord_trans = reinterpret_cast<CoordTransVec1>(dlsym(h, "coord_trans_vec1_"));
    if (g_get_field == nullptr || g_coord_trans == nullptr) {
        std::fprintf(stderr, "missing entry points in %s\n", lib.c_str());
        return 1;
    }
    std::printf("space.irbem T96 vs IRBEM kext=7, oracle %s\n\n", lib.c_str());
    if (pass == "probe") {
        probe();
    } else if (pass == "fit") {
        const int iters = argc > 3 ? std::atoi(argv[3]) : 300;
        const std::string out = argc > 4 ? argv[4] : "/tmp/t96_generated.inc";
        fit(iters, out);
    } else {
        report();
    }
    return 0;
}
