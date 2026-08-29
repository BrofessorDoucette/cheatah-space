// Differential study of space.irbem's Olson-Pfitzer dynamic model against the IRBEM oracle's
// `kext = 6`, and the black-box experiments that say what the oracle IS before asking how far away
// it is.
//
// DEV-ONLY. Like t89_diff.cpp beside it, this is the one kind of program that touches IRBEM: it is
// never built by the QA gate and never ships. IRBEM is LGPL-3.0 and cheatah-space is MIT, so the
// library is run here as a BLACK BOX (dlopen plus the documented C entry points), never read for
// its logic and never linked into anything distributed.
//
// WHAT IT MEASURES, in order:
//
//   1. STRUCTURE. Four black-box facts about the oracle, each a numerical identity that either
//      holds to roundoff or does not: (a) `n` and `V` enter only through `n V^2`; (b) the field is
//      affine in Dst; (c) the Dst gradient does not depend on the pressure; (d) it refuses (baddata)
//      strictly outside `5..50`, `300..500`, `-100..20`, `r <= 60`. `space/irbem/ext_opd.hpp`
//      reproduces (a)-(c) by construction and reports (d) instead of refusing.
//   2. NON-SIMILARITY. The rank of the oracle's `n V^2` family over a point set, and the best fit
//      of `a B(s r; P_ref) + g G(r)` per pressure. A pure self-similar compression of one field
//      would give rank ~1-2 and a roundoff residual; the measurement says otherwise, which is why
//      the published-structure model cannot reach parity and the numbers below are a model-family
//      gap, not a bug.
//   3. DEVIATION per regime, sweeping the drivers CONTINUOUSLY (this is a continuous-driver model)
//      through the corpus's four regimes clipped to the model's documented envelope, plus dense
//      samples along the Dst axis, at three dipole tilts, in two regions.
//   4. THE FLOOR. How close the published-structure form can get with its 20 linear amplitudes
//      free (the 19 quiet coefficients under the compression, plus the ring), driver by driver.
//      That is the gap the STRUCTURE leaves; the deviation in (3) is the gap the published
//      constants leave on top of it.
//
// Build (from the repository root):
//   g++ -O2 -std=c++20 tools/oracle/opd_diff.cpp -I. \
//       -I$CHEATAH_DIR/stdlib/ndarray -I$CHEATAH_DIR/stdlib/builtins -I$CHEATAH_DIR/stdlib/fixarray \
//       -o /tmp/opd_diff && /tmp/opd_diff /tmp/irbem-builds/libirbem-O2.so
#include <dlfcn.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

#include "space/irbem/ext_opd.hpp"

namespace {

namespace ir = cheatah::space::irbem;

/// `get_field1_`, as documented in the vendored matlab/libirbem.h.
using GetField1 = void (*)(int*, int*, int*, int*, int*, double*, double*, double*, double*,
                           double*, double*, double*);
/// `coord_trans_vec1_`, likewise.
using CoordTransVec1 = void (*)(int*, int*, int*, int*, int*, double*, double*, double*);

GetField1 g_get_field = nullptr;
CoordTransVec1 g_ctv = nullptr;
int g_year = 2015;
int g_doy = 80;
double g_ut = 39183.0;

/// The oracle's own dipole tilt at the current epoch: the SM z axis in GSM is (sin psi, 0, cos psi).
double oracle_tilt() {
    int one = 1;
    int si = 4;
    int so = 2;
    std::array<double, 3> in{0.0, 0.0, 1.0};
    std::array<double, 3> out{};
    g_ctv(&one, &si, &so, &g_year, &g_doy, &g_ut, in.data(), out.data());
    return std::atan2(out[0], out[2]);
}

/// The oracle's external field for @p kext at a GSM point, in GSM nT, isolated as kext minus 0.
/// A component below -1e30 is the oracle's baddata refusal.
std::array<double, 3> oracle_ext(int kext, double x, double y, double z, double dsw, double vsw,
                                 double dst) {
    int one = 1;
    int si = 2;
    int so = 1;
    std::array<double, 3> gsm{x, y, z};
    std::array<double, 3> geo{};
    g_ctv(&one, &si, &so, &g_year, &g_doy, &g_ut, gsm.data(), geo.data());
    std::array<int, 5> options{0, 0, 0, 0, 0};
    int sysaxes = 1;
    int k0 = 0;
    int k = kext;
    std::vector<double> mag(25, 0.0);
    mag[1] = dst;
    mag[2] = dsw;
    mag[3] = vsw;
    std::array<double, 3> b0{};
    std::array<double, 3> bk{};
    double m0 = 0.0;
    double mk = 0.0;
    double x1 = geo[0];
    double x2 = geo[1];
    double x3 = geo[2];
    g_get_field(&k0, options.data(), &sysaxes, &g_year, &g_doy, &g_ut, &x1, &x2, &x3, mag.data(),
                b0.data(), &m0);
    g_get_field(&k, options.data(), &sysaxes, &g_year, &g_doy, &g_ut, &x1, &x2, &x3, mag.data(),
                bk.data(), &mk);
    if (bk[0] < -1e30) return {bk[0], bk[1], bk[2]};
    std::array<double, 3> dgeo{bk[0] - b0[0], bk[1] - b0[1], bk[2] - b0[2]};
    std::array<double, 3> out{};
    si = 1;
    so = 2;
    g_ctv(&one, &si, &so, &g_year, &g_doy, &g_ut, dgeo.data(), out.data());
    return out;
}

bool is_bad(const std::array<double, 3>& b) { return b[0] < -1e30; }

double norm(const std::array<double, 3>& v) {
    return std::sqrt((v[0] * v[0]) + (v[1] * v[1]) + (v[2] * v[2]));
}

double dist(const std::array<double, 3>& a, const std::array<double, 3>& b) {
    return std::sqrt(((a[0] - b[0]) * (a[0] - b[0])) + ((a[1] - b[1]) * (a[1] - b[1])) +
                     ((a[2] - b[2]) * (a[2] - b[2])));
}

/// A deterministic scatter inside a shell, so the sweeps are reproducible.
std::vector<std::array<double, 3>> scatter(std::size_t n, double r_lo, double r_hi,
                                           std::uint64_t seed) {
    std::vector<std::array<double, 3>> out;
    std::uint64_t s = seed;
    const auto next = [&s] {
        s = (s * 6364136223846793005ULL) + 1442695040888963407ULL;
        return static_cast<double>(s >> 11) / 9007199254740992.0;
    };
    for (std::size_t i = 0; i < n; ++i) {
        const double r = r_lo + ((r_hi - r_lo) * next());
        const double th = std::acos(1.0 - (2.0 * next()));
        const double ph = 6.283185307179586 * next();
        out.push_back({r * std::sin(th) * std::cos(ph), r * std::sin(th) * std::sin(ph),
                       r * std::cos(th)});
    }
    return out;
}

/// Whether a GSM point is inside the magnetopause for the given solar wind, by the published
/// Shue et al. (1997, JGR 102, 9497) shape `r = r_0 (2 / (1 + cos theta))^alpha` with
/// `alpha = 0.58` and `r_0` the pressure-balance standoff `11.15 (P / P_0)^(-1/6)` Re — the
/// Kivelson & Russell (1995) `107.4 (n V^2)^(-1/6)` at the reference (5, 400). A 10% margin
/// keeps both models away from a boundary neither is meaningful beyond: the oracle's quiet-model
/// exponentials diverge outside it, and comparing two models where neither means anything would
/// measure nothing. The clipping is stated in the output.
bool inside_magnetopause(double x, double y, double z, double dsw, double vsw) {
    const double r = std::sqrt((x * x) + (y * y) + (z * z));
    if (r < 1e-9) return true;
    const double r0 = 11.15 * std::pow(ir::opd_pressure_npa(dsw, vsw) / ir::opd_reference_pressure_npa(), -1.0 / 6.0);
    const double cos_th = x / r;
    const double mp = r0 * std::pow(2.0 / (1.0 + cos_th), 0.58);
    return r < 0.9 * mp;
}

/// One oracle sample for the floor experiment.
struct Sample {
    double tilt;
    double x, y, z;
    double dsw, vsw, dst;
    std::array<double, 3> b;
};

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

/// The 20 basis fields of the published-structure form at one sample: the 19 quiet coefficients
/// one-hot under the compression, and the unit ring. Built through the SHIPPING evaluator so the
/// floor is a property of the header and not of a copy of it.
constexpr int kBasis = static_cast<int>(ir::t89_linear_count) + 1;
void basis(const Sample& s, std::array<std::array<double, 3>, kBasis>& out) {
    const double sn = std::sin(s.tilt);
    const double cs = std::cos(s.tilt);
    ir::OpdParameters<double> p = ir::opd_parameters<double>(s.dsw, s.vsw, s.dst);
    for (int k = 0; k < kBasis - 1; ++k) {
        ir::OpdParameters<double> one = p;
        one.quiet.c.fill(0.0);
        one.quiet.c[static_cast<std::size_t>(k)] = 1.0;
        one.ring = 0.0;
        out[static_cast<std::size_t>(k)] = ir::opd_components<double>(one, sn, cs, s.x, s.y, s.z);
    }
    ir::OpdParameters<double> ring = p;
    ring.quiet.c.fill(0.0);
    ring.ring = 1.0;
    out[kBasis - 1] = ir::opd_components<double>(ring, sn, cs, s.x, s.y, s.z);
}

/// RMS residual of the best 20-amplitude fit of the structure to @p samples, as a fraction of the
/// oracle's RMS external field. The amplitudes are shared across the whole sample set — one model
/// for all drivers — which is what "the structure" means.
double floor_residual(const std::vector<Sample>& samples, double& signal_rms) {
    const int n = kBasis;
    std::vector<double> ata(static_cast<std::size_t>(n) * n, 0.0);
    std::vector<double> atb(static_cast<std::size_t>(n), 0.0);
    std::vector<std::array<double, kBasis>> rows;
    std::vector<double> obs;
    std::array<std::array<double, 3>, kBasis> bs{};
    double sig2 = 0.0;
    for (const Sample& s : samples) {
        basis(s, bs);
        for (int c = 0; c < 3; ++c) {
            std::array<double, kBasis> row{};
            for (int k = 0; k < n; ++k) row[static_cast<std::size_t>(k)] = bs[static_cast<std::size_t>(k)][static_cast<std::size_t>(c)];
            for (int i = 0; i < n; ++i) {
                atb[static_cast<std::size_t>(i)] += row[static_cast<std::size_t>(i)] * s.b[static_cast<std::size_t>(c)];
                for (int j = 0; j < n; ++j)
                    ata[(static_cast<std::size_t>(i) * n) + j] += row[static_cast<std::size_t>(i)] * row[static_cast<std::size_t>(j)];
            }
            rows.push_back(row);
            obs.push_back(s.b[static_cast<std::size_t>(c)]);
            sig2 += s.b[static_cast<std::size_t>(c)] * s.b[static_cast<std::size_t>(c)];
        }
    }
    signal_rms = std::sqrt(sig2 / static_cast<double>(obs.size()));
    if (!solve(ata, atb, n)) return 1e30;
    double res = 0.0;
    for (std::size_t r = 0; r < obs.size(); ++r) {
        double f = 0.0;
        for (int k = 0; k < n; ++k) f += rows[r][static_cast<std::size_t>(k)] * atb[static_cast<std::size_t>(k)];
        res += (f - obs[r]) * (f - obs[r]);
    }
    return std::sqrt(res / static_cast<double>(obs.size()));
}

/// The four corpus regimes' (n, V, Dst), CLIPPED to the model's documented envelope — the storm
/// and extreme regimes lie outside it, and the oracle refuses there, so what is compared is the
/// model at the edge of its own envelope, which is where it is used in practice.
struct Regime {
    const char* name;
    double dsw, vsw, dst;
};
constexpr std::array<Regime, 4> kRegimes{{
    {"quiet    (n=5,  V=380, Dst=-8)", 5.0, 380.0, -8.0},
    {"moderate (n=8,  V=450, Dst=-42)", 8.0, 450.0, -42.0},
    {"storm    (n=20, V=500, Dst=-100) [Vsw, Dst clipped to envelope]", 20.0, 500.0, -100.0},
    {"extreme  (n=45, V=500, Dst=-100) [Vsw, Dst clipped to envelope]", 45.0, 500.0, -100.0},
}};

}  // namespace

int main(int argc, char** argv) {
    const std::string lib = argc > 1 ? argv[1] : "/tmp/irbem-builds/libirbem-O2.so";
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
    std::printf("space.irbem Olson-Pfitzer dynamic (documented structure) vs IRBEM kext=6, oracle %s\n",
                lib.c_str());

    // ---- 1. STRUCTURE: what the oracle is ---------------------------------------------------
    std::printf("\n=== 1. black-box structure of the oracle (epoch doy 80, tilt %.3f deg) ===\n",
                oracle_tilt() * 57.29577951308232);
    // Inside the magnetopause at EVERY pressure the probes use (standoff 7.6 Re at n = 50), so no
    // sample sits where the oracle's quiet-model exponentials have run away.
    const std::vector<std::array<double, 3>> probe = scatter(40, 2.5, 6.5, 0x51ED270693C1A2B7ULL);
    {
        // (a) n and V only through n V^2: (20, 400) vs (32, sqrt(1e5)) have the same product.
        double worst = 0.0;
        double worst_rel = 0.0;
        for (const auto& p : probe) {
            const auto a = oracle_ext(6, p[0], p[1], p[2], 20.0, 400.0, -30.0);
            const auto b = oracle_ext(6, p[0], p[1], p[2], 32.0, 316.22776601683794, -30.0);
            worst = std::max(worst, dist(a, b));
            worst_rel = std::max(worst_rel, dist(a, b) / (norm(a) + 1e-12));
        }
        std::printf("  (a) equal n V^2, different (n, V): max |dB| = %.3e nT (rel %.3e)\n", worst,
                    worst_rel);
        // (b) affine in Dst: the midpoint identity.
        double worst_aff = 0.0;
        for (const auto& p : probe) {
            const auto b0 = oracle_ext(6, p[0], p[1], p[2], 10.0, 450.0, 0.0);
            const auto b1 = oracle_ext(6, p[0], p[1], p[2], 10.0, 450.0, -100.0);
            const auto bm = oracle_ext(6, p[0], p[1], p[2], 10.0, 450.0, -50.0);
            const std::array<double, 3> lin{0.5 * (b0[0] + b1[0]), 0.5 * (b0[1] + b1[1]),
                                            0.5 * (b0[2] + b1[2])};
            worst_aff = std::max(worst_aff, dist(bm, lin));
        }
        std::printf("  (b) affine in Dst (midpoint identity): max |dB| = %.3e nT\n", worst_aff);
        // (c) the Dst gradient is pressure-independent.
        double worst_g = 0.0;
        double ref_g = 0.0;
        for (const auto& p : probe) {
            const auto l0 = oracle_ext(6, p[0], p[1], p[2], 5.0, 400.0, 0.0);
            const auto l1 = oracle_ext(6, p[0], p[1], p[2], 5.0, 400.0, -100.0);
            const auto h0 = oracle_ext(6, p[0], p[1], p[2], 50.0, 400.0, 0.0);
            const auto h1 = oracle_ext(6, p[0], p[1], p[2], 50.0, 400.0, -100.0);
            std::array<double, 3> gl{};
            std::array<double, 3> gh{};
            for (int c = 0; c < 3; ++c) {
                gl[static_cast<std::size_t>(c)] = (l1[static_cast<std::size_t>(c)] - l0[static_cast<std::size_t>(c)]) / -100.0;
                gh[static_cast<std::size_t>(c)] = (h1[static_cast<std::size_t>(c)] - h0[static_cast<std::size_t>(c)]) / -100.0;
            }
            worst_g = std::max(worst_g, dist(gl, gh));
            ref_g = std::max(ref_g, norm(gl));
        }
        std::printf("  (c) dB/dDst at n=5 vs n=50: max |dG| = %.3e nT/nT (|G| up to %.3f)\n", worst_g,
                    ref_g);
        // (d) the refusals.
        const auto bad = [](std::array<double, 3> b) { return is_bad(b) ? "baddata" : "value"; };
        std::printf("  (d) oracle at n=50: %s, 50.0001: %s | n=5: %s, 4.9999: %s\n",
                    bad(oracle_ext(6, 6.6, 0, 0, 50.0, 400, 0)), bad(oracle_ext(6, 6.6, 0, 0, 50.0001, 400, 0)),
                    bad(oracle_ext(6, 6.6, 0, 0, 5.0, 400, 0)), bad(oracle_ext(6, 6.6, 0, 0, 4.9999, 400, 0)));
        std::printf("      V=500: %s, 500.001: %s | V=300: %s, 299.999: %s\n",
                    bad(oracle_ext(6, 6.6, 0, 0, 10, 500.0, 0)), bad(oracle_ext(6, 6.6, 0, 0, 10, 500.001, 0)),
                    bad(oracle_ext(6, 6.6, 0, 0, 10, 300.0, 0)), bad(oracle_ext(6, 6.6, 0, 0, 10, 299.999, 0)));
        std::printf("      Dst=20: %s, 20.001: %s | Dst=-100: %s, -100.001: %s\n",
                    bad(oracle_ext(6, 6.6, 0, 0, 10, 400, 20.0)), bad(oracle_ext(6, 6.6, 0, 0, 10, 400, 20.001)),
                    bad(oracle_ext(6, 6.6, 0, 0, 10, 400, -100.0)), bad(oracle_ext(6, 6.6, 0, 0, 10, 400, -100.001)));
        std::printf("      r=59.9: %s, r=60.05: %s\n", bad(oracle_ext(6, -59.9, 0, 0, 10, 400, 0)),
                    bad(oracle_ext(6, -60.05, 0, 0, 10, 400, 0)));
    }

    // ---- 2. NON-SIMILARITY ------------------------------------------------------------------
    std::printf("\n=== 2. is the oracle a similarity transform of one field? ===\n");
    {
        const std::vector<double> dvals{5, 7, 10, 14, 20, 28, 40, 50};
        const int n = static_cast<int>(dvals.size());
        std::vector<std::vector<double>> rows;
        for (double d : dvals) {
            std::vector<double> row;
            for (const auto& p : probe) {
                const auto b = oracle_ext(6, p[0], p[1], p[2], d, 400.0, 0.0);
                row.push_back(b[0]);
                row.push_back(b[1]);
                row.push_back(b[2]);
            }
            rows.push_back(row);
        }
        // Gram matrix, Jacobi eigenvalues -> singular values of the 8-field family.
        std::vector<double> a(static_cast<std::size_t>(n) * n, 0.0);
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j) {
                double s = 0.0;
                for (std::size_t k = 0; k < rows[0].size(); ++k) s += rows[static_cast<std::size_t>(i)][k] * rows[static_cast<std::size_t>(j)][k];
                a[(static_cast<std::size_t>(i) * n) + j] = s;
            }
        for (int sweep = 0; sweep < 100; ++sweep) {
            double off = 0.0;
            for (int i = 0; i < n; ++i)
                for (int j = i + 1; j < n; ++j) off += a[(static_cast<std::size_t>(i) * n) + j] * a[(static_cast<std::size_t>(i) * n) + j];
            if (off < 1e-20) break;
            for (int p = 0; p < n; ++p)
                for (int q = p + 1; q < n; ++q) {
                    const double apq = a[(static_cast<std::size_t>(p) * n) + q];
                    if (std::fabs(apq) < 1e-300) continue;
                    const double th = 0.5 * std::atan2(2.0 * apq, a[(static_cast<std::size_t>(q) * n) + q] - a[(static_cast<std::size_t>(p) * n) + p]);
                    const double c = std::cos(th);
                    const double s = std::sin(th);
                    for (int k = 0; k < n; ++k) {
                        const double akp = a[(static_cast<std::size_t>(k) * n) + p];
                        const double akq = a[(static_cast<std::size_t>(k) * n) + q];
                        a[(static_cast<std::size_t>(k) * n) + p] = (c * akp) - (s * akq);
                        a[(static_cast<std::size_t>(k) * n) + q] = (s * akp) + (c * akq);
                    }
                    for (int k = 0; k < n; ++k) {
                        const double apk = a[(static_cast<std::size_t>(p) * n) + k];
                        const double aqk = a[(static_cast<std::size_t>(q) * n) + k];
                        a[(static_cast<std::size_t>(p) * n) + k] = (c * apk) - (s * aqk);
                        a[(static_cast<std::size_t>(q) * n) + k] = (s * apk) + (c * aqk);
                    }
                }
        }
        std::vector<double> sv;
        for (int i = 0; i < n; ++i) sv.push_back(std::sqrt(std::max(0.0, a[(static_cast<std::size_t>(i) * n) + i])));
        std::sort(sv.rbegin(), sv.rend());
        std::printf("  singular values of the 8-pressure family (n = 5..50 at V = 400, Dst = 0):\n   ");
        for (double v : sv) std::printf(" %.3e", v / sv[0]);
        std::printf("  (relative to the largest)\n");
        std::printf("  a rank-1 family (pure amplitude scaling) would show one value; the family needs\n"
                    "  ~4-5 to reach 1e-3, so no single field's similarity transform reproduces it.\n");

        // The best per-pressure fit of a B(s r; P_ref) + g G(r), s scanned.
        std::printf("  best fit B(r; n) = a B(s r; n=5) + g G(r), s and (a, g) free per pressure:\n");
        for (int which = 0; which < 2; ++which) {
        // Two probe sets: the inner magnetosphere (2.5..6.5 Re, all local times), and the
        // NIGHTSIDE out to 15 Re, inside the magnetopause at every pressure used, to see whether
        // the geometric rescaling the published idea predicts shows up where the tail is.
        std::vector<std::array<double, 3>> probe = which == 0 ? scatter(40, 2.5, 6.5, 0x51ED270693C1A2B7ULL)
                                                              : scatter(60, 3.0, 15.0, 0x7C1A2B351ED27069ULL);
        if (which == 1) for (auto& p : probe) p[0] = -std::fabs(p[0]);
        std::vector<std::array<double, 3>> G(probe.size());
        for (std::size_t i = 0; i < probe.size(); ++i) {
            const auto b0 = oracle_ext(6, probe[i][0], probe[i][1], probe[i][2], 5.0, 400.0, 0.0);
            const auto b1 = oracle_ext(6, probe[i][0], probe[i][1], probe[i][2], 5.0, 400.0, -100.0);
            for (int c = 0; c < 3; ++c) G[i][static_cast<std::size_t>(c)] = (b1[static_cast<std::size_t>(c)] - b0[static_cast<std::size_t>(c)]) / -100.0;
        }
        std::printf("   %s:\n", which == 0 ? "inner magnetosphere, 2.5..6.5 Re" : "nightside, 3..15 Re");
        for (double d : {10.0, 20.0, 40.0, 50.0}) {
            std::vector<std::array<double, 3>> tgt(probe.size());
            for (std::size_t i = 0; i < probe.size(); ++i) tgt[i] = oracle_ext(6, probe[i][0], probe[i][1], probe[i][2], d, 400.0, 0.0);
            double best = 1e30;
            double bs = 0.0;
            double ba = 0.0;
            double bg = 0.0;
            for (int is = 0; is <= 160; ++is) {
                const double s = 0.85 + (0.005 * is);
                double a11 = 0.0, a12 = 0.0, a22 = 0.0, r1 = 0.0, r2 = 0.0;
                std::vector<std::array<double, 3>> q(probe.size());
                for (std::size_t i = 0; i < probe.size(); ++i) {
                    q[i] = oracle_ext(6, probe[i][0] * s, probe[i][1] * s, probe[i][2] * s, 5.0, 400.0, 0.0);
                    for (int c = 0; c < 3; ++c) {
                        const std::size_t cc = static_cast<std::size_t>(c);
                        a11 += q[i][cc] * q[i][cc];
                        a12 += q[i][cc] * G[i][cc];
                        a22 += G[i][cc] * G[i][cc];
                        r1 += q[i][cc] * tgt[i][cc];
                        r2 += G[i][cc] * tgt[i][cc];
                    }
                }
                const double det = (a11 * a22) - (a12 * a12);
                const double aa = ((a22 * r1) - (a12 * r2)) / det;
                const double gg = ((-a12 * r1) + (a11 * r2)) / det;
                double res = 0.0;
                double sig = 0.0;
                for (std::size_t i = 0; i < probe.size(); ++i)
                    for (int c = 0; c < 3; ++c) {
                        const std::size_t cc = static_cast<std::size_t>(c);
                        const double f = (aa * q[i][cc]) + (gg * G[i][cc]);
                        res += (tgt[i][cc] - f) * (tgt[i][cc] - f);
                        sig += tgt[i][cc] * tgt[i][cc];
                    }
                if (res / sig < best) {
                    best = res / sig;
                    bs = s;
                    ba = aa;
                    bg = gg;
                }
            }
            std::printf("    n=%3.0f: best s = %.3f (P^(1/6) = %.3f), a = %.3f (s^3 = %.3f), g = %+.1f nT, "
                        "residual %.1f%% of the field\n",
                        d, bs, std::pow(d / 5.0, 1.0 / 6.0), ba, bs * bs * bs, bg, 100.0 * std::sqrt(best));
        }
        }
    }

    // ---- 3. DEVIATION per regime, drivers swept continuously ----------------------------------
    std::printf("\n=== 3. deviation of the published-structure model, per regime ===\n");
    struct Epoch {
        int doy;
        double ut;
    };
    const std::array<Epoch, 3> epochs{{{80, 39183.0}, {180, 43200.0}, {355, 7200.0}}};
    std::vector<Sample> floor_samples;
    for (int region = 0; region < 2; ++region) {
        const double r_hi = region == 0 ? 10.0 : 35.0;
        std::printf("\n--- %s, inside 0.9 x the Shue (1997) magnetopause ---\n",
                    region == 0 ? "radiation belts, 3 <= r <= 10 Re" : "full box, 3 <= r <= 35 Re");
        std::printf("%-66s %6s %9s %9s %9s %8s\n", "regime", "N", "rms|dB|", "p99|dB|", "max|dB|",
                    "rms rel");
        for (const Regime& rg : kRegimes) {
            std::vector<double> abs_err;
            double sum2 = 0.0;
            double sig2 = 0.0;
            for (const Epoch& e : epochs) {
                g_doy = e.doy;
                g_ut = e.ut;
                const double ps = oracle_tilt();
                // The drivers are swept CONTINUOUSLY around the regime's representative point:
                // nine (n, V, Dst) triples per regime spread over +-20% in n, +-8% in V and
                // +-15 nT in Dst, all clipped to the envelope.
                for (int k = 0; k < 9; ++k) {
                    const double fn = 0.8 + (0.05 * k);
                    const double fv = 0.92 + (0.02 * k);
                    const double dd = -15.0 + (3.75 * k);
                    const double dsw = std::clamp(rg.dsw * fn, 5.0, 50.0);
                    const double vsw = std::clamp(rg.vsw * fv, 300.0, 500.0);
                    const double dst = std::clamp(rg.dst + dd, -100.0, 20.0);
                    const std::vector<std::array<double, 3>> pts =
                        scatter(60, 3.0, r_hi, 0xA5A5A5A5DEADBEEFULL + static_cast<std::uint64_t>(k));
                    for (const auto& p : pts) {
                        if (!inside_magnetopause(p[0], p[1], p[2], dsw, vsw)) continue;
                        const auto ora = oracle_ext(6, p[0], p[1], p[2], dsw, vsw, dst);
                        if (is_bad(ora)) continue;
                        const ir::FieldVector<ir::Frame::GSM> mine = ir::opd_field_at(
                            ir::Position<ir::Frame::GSM>{cheatah::fixarray::vec3d{p[0], p[1], p[2]}},
                            std::sin(ps), std::cos(ps), dsw, vsw, dst);
                        double d2 = 0.0;
                        double o2 = 0.0;
                        for (int c = 0; c < 3; ++c) {
                            const std::size_t cc = static_cast<std::size_t>(c);
                            d2 += (mine.v[c] - ora[cc]) * (mine.v[c] - ora[cc]);
                            o2 += ora[cc] * ora[cc];
                        }
                        abs_err.push_back(std::sqrt(d2));
                        sum2 += d2;
                        sig2 += o2;
                        if (region == 1 && (k % 4) == 0) floor_samples.push_back({ps, p[0], p[1], p[2], dsw, vsw, dst, ora});
                    }
                }
            }
            std::sort(abs_err.begin(), abs_err.end());
            const std::size_t n = abs_err.size();
            std::printf("%-66s %6zu %9.3f %9.3f %9.3f %8.4f\n", rg.name, n,
                        std::sqrt(sum2 / static_cast<double>(n)), abs_err[(n * 99) / 100],
                        abs_err[n - 1], std::sqrt(sum2 / sig2));
        }
    }
    // Dense Dst axis at fixed pressure, geosynchronous shell: where the ring term is the story.
    g_doy = 80;
    g_ut = 39183.0;
    {
        const double ps = oracle_tilt();
        std::printf("\n--- dense Dst sweep, r = 6.6 Re shell, n = 10, V = 450, tilt %.2f deg ---\n",
                    ps * 57.29577951308232);
        std::printf("%8s %9s %8s | %s\n", "Dst", "rms|dB|", "rms rel", "dBz/dDst at noon: model vs oracle");
        const std::vector<std::array<double, 3>> pts = scatter(80, 6.6, 6.6, 0x0123456789ABCDEFULL);
        for (int k = 0; k <= 12; ++k) {
            const double dst = -100.0 + (10.0 * k);
            double sum2 = 0.0;
            double sig2 = 0.0;
            for (const auto& p : pts) {
                const auto ora = oracle_ext(6, p[0], p[1], p[2], 10.0, 450.0, dst);
                const ir::FieldVector<ir::Frame::GSM> mine = ir::opd_field_at(
                    ir::Position<ir::Frame::GSM>{cheatah::fixarray::vec3d{p[0], p[1], p[2]}},
                    std::sin(ps), std::cos(ps), 10.0, 450.0, dst);
                for (int c = 0; c < 3; ++c) {
                    const std::size_t cc = static_cast<std::size_t>(c);
                    sum2 += (mine.v[c] - ora[cc]) * (mine.v[c] - ora[cc]);
                    sig2 += ora[cc] * ora[cc];
                }
            }
            const auto o0 = oracle_ext(6, 6.6, 0, 0, 10.0, 450.0, 0.0);
            const auto o1 = oracle_ext(6, 6.6, 0, 0, 10.0, 450.0, -100.0);
            const ir::Position<ir::Frame::GSM> noon{cheatah::fixarray::vec3d{6.6, 0.0, 0.0}};
            const double m0 = ir::opd_field_at(noon, std::sin(ps), std::cos(ps), 10.0, 450.0, 0.0).v[2];
            const double m1 = ir::opd_field_at(noon, std::sin(ps), std::cos(ps), 10.0, 450.0, -100.0).v[2];
            std::printf("%8.1f %9.3f %8.4f | %.4f vs %.4f\n", dst, std::sqrt(sum2 / static_cast<double>(pts.size())),
                        std::sqrt(sum2 / sig2), (m1 - m0) / -100.0, (o1[2] - o0[2]) / -100.0);
        }
    }

    // ---- 4. THE FLOOR -------------------------------------------------------------------------
    std::printf("\n=== 4. the floor: the published structure with all 20 amplitudes free ===\n");
    {
        double signal = 0.0;
        const double res = floor_residual(floor_samples, signal);
        // The shipped constants' own residual over the SAME samples, for the comparison.
        double ship2 = 0.0;
        for (const Sample& s : floor_samples) {
            const ir::FieldVector<ir::Frame::GSM> mine = ir::opd_field_at(
                ir::Position<ir::Frame::GSM>{cheatah::fixarray::vec3d{s.x, s.y, s.z}}, std::sin(s.tilt),
                std::cos(s.tilt), s.dsw, s.vsw, s.dst);
            for (int c = 0; c < 3; ++c) ship2 += (mine.v[c] - s.b[static_cast<std::size_t>(c)]) * (mine.v[c] - s.b[static_cast<std::size_t>(c)]);
        }
        const double ship = std::sqrt(ship2 / static_cast<double>(3 * floor_samples.size()));
        std::printf("  %zu samples over all regimes, tilts and the 3..35 Re box; oracle rms field %.3f nT\n",
                    floor_samples.size(), signal);
        std::printf("  published constants, nothing free     : rms residual %.3f nT (%.1f%% of signal)\n",
                    ship, 100.0 * ship / signal);
        std::printf("  one shared set of 20 amplitudes free  : rms residual %.3f nT (%.1f%% of signal)\n",
                    res, 100.0 * res / signal);
        std::printf("\n  The second line is the gap the STRUCTURE leaves: no choice of amplitudes for the\n"
                    "  compressed-quiet-field-plus-ring form gets closer than that to the distributed\n"
                    "  model. The first minus the second is what the published constants cost on top.\n"
                    "  Neither is a defect of the arithmetic; the div-B and symmetry suites establish that.\n");
    }
    return 0;
}
