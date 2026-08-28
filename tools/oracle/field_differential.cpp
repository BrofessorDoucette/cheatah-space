// Differential study of space/irbem/field.hpp against the IRBEM oracle — the measurement behind
// every number in that header's @file block, and the source of the golden vectors in
// tests/irbem_field_test.cpp.
//
// DEV-ONLY, exactly like convergence.cpp: IRBEM is LGPL-3.0 and cheatah-space is MIT, so the
// library is run here as a BLACK BOX (dlopen + the documented C entry points), never linked into
// anything we ship and never read for its logic. Nothing in the QA gate builds this file.
//
// What it establishes, in the order the header cites it:
//
//   1. GET_FIELD_MULTI is the internal field at the point, and with options(2)=0 the epoch is
//      year.5 exactly (the published options table says "initialize IGRF field once per year
//      (year.5)"). Sweeping the truncation degree 8..13 collapses the residual at 10, which is
//      how IRBEM's internal truncation was identified without reading a line of Fortran.
//   2. GET_BDERIVS is a FORWARD difference of GET_FIELD_MULTI at the caller's dX -- reproduced
//      from the oracle's own field values, bit for bit.
//   3. COMPUTE_GRAD_CURV_CURL's eight outputs, each reproduced from field.hpp's algebra, under
//      STORM conditions: Kp 0..9, Dst 0..-400 nT, Pdyn 0.5..40 nPa, and Bz swept densely
//      southward, across kext 4 (T89), 7 (T96) and 10 (T01 storm). The algebra is validated
//      independently of which field model produced B, which is what makes a storm sweep
//      meaningful for a module that does not yet implement an external field.
//   4. GET_HEMI_MULTI is the sign of d|B|/ds along +Bhat.
//
// Build (needs no CMake):
//   g++ -std=c++20 -O2 -ffp-contract=off tools/oracle/field_differential.cpp -I. \
//       -I$CHEATAH_DIR/stdlib/ndarray -I$CHEATAH_DIR/stdlib/builtins -I$CHEATAH_DIR/stdlib/fixarray \
//       -ldl -o /tmp/field_differential && /tmp/field_differential
#include <dlfcn.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <random>
#include <string>
#include <vector>

#include "space/irbem/field.hpp"

using namespace cheatah::space::irbem;
namespace fx = cheatah::fixarray;

namespace {

using GetFieldMulti = void (*)(int*, int*, int*, int*, int*, int*, double*, double*, double*,
                               double*, double*, double*, double*);
using GetBderivs = void (*)(int*, int*, int*, int*, double*, int*, int*, double*, double*, double*,
                            double*, double*, double*, double*, double*, double*);
using GetHemiMulti = void (*)(int*, int*, int*, int*, int*, int*, double*, double*, double*,
                              double*, double*, int*);
using ComputeGradCurvCurl = void (*)(int*, double*, double*, double*, double*, double*, double*,
                                     double*, double*, double*, double*, double*, double*);

GetFieldMulti g_field = nullptr;
GetBderivs g_bderivs = nullptr;
GetHemiMulti g_hemi = nullptr;
ComputeGradCurvCurl g_gcc = nullptr;

/// 2015 day 180, 12:00 UT. With options(2)=0 IRBEM initializes IGRF at year.5, so this is the
/// epoch every comparison below builds the model at.
constexpr double kEpoch = 2015.5;

struct Storm {
    const char* name;
    int kext;
    std::array<double, 25> maginput;
};

/// Fill the maginput slots by IRBEM's 1-based index (Kp is 1, Dst 2, Pdyn 5, By 6, Bz 7).
std::array<double, 25> drivers(double kp10, double dst, double pdyn, double by, double bz) {
    std::array<double, 25> m{};
    m[0] = kp10;
    m[1] = dst;
    m[2] = 5.0;    // Dsw
    m[3] = 400.0;  // Vsw
    m[4] = pdyn;
    m[5] = by;
    m[6] = bz;
    m[7] = 5.0;   // G1
    m[8] = 5.0;   // G2
    m[9] = 5.0;   // G3
    return m;
}

double rel(double a, double b) { return std::abs(a - b) / std::max(1e-12, std::abs(b)); }

}  // namespace

int main(int argc, char** argv) {
    const std::string lib = argc > 1 ? argv[1] : "/tmp/irbem-builds/libirbem-O2.so";
    void* h = dlopen(lib.c_str(), RTLD_NOW);
    if (h == nullptr) {
        std::fprintf(stderr, "cannot dlopen %s: %s\n", lib.c_str(), dlerror());
        return 1;
    }
    g_field = reinterpret_cast<GetFieldMulti>(dlsym(h, "get_field_multi_"));
    g_bderivs = reinterpret_cast<GetBderivs>(dlsym(h, "get_bderivs_"));
    g_hemi = reinterpret_cast<GetHemiMulti>(dlsym(h, "get_hemi_multi_"));
    g_gcc = reinterpret_cast<ComputeGradCurvCurl>(dlsym(h, "compute_grad_curv_curl_"));
    if (g_field == nullptr || g_bderivs == nullptr || g_hemi == nullptr || g_gcc == nullptr) {
        std::fprintf(stderr, "missing entry points in %s\n", lib.c_str());
        return 1;
    }

    std::mt19937 rng(20260828);
    std::uniform_real_distribution<double> u(-1.0, 1.0);
    const std::size_t n = 300;
    std::vector<double> x1(n);
    std::vector<double> x2(n);
    std::vector<double> x3(n);
    std::vector<Position<Frame::GEO>> pts(n);
    for (std::size_t i = 0; i < n; ++i) {
        const double a = u(rng) * 3.14159265358979;
        const double b = u(rng) * 1.45;
        const double r = 1.5 + (std::abs(u(rng)) * 7.0);
        x1[i] = r * std::cos(b) * std::cos(a);
        x2[i] = r * std::cos(b) * std::sin(a);
        x3[i] = r * std::sin(b);
        pts[i] = Position<Frame::GEO>{fx::vec3d{x1[i], x2[i], x3[i]}};
    }

    std::vector<int> iyear(n, 2015);
    std::vector<int> idoy(n, 180);
    std::vector<double> ut(n, 43200.0);
    std::array<int, 5> options{0, 0, 0, 0, 0};
    int ntime = static_cast<int>(n);
    int sysaxes = 1;

    auto oracle_field = [&](int kext, const std::array<double, 25>& drv, std::vector<double>& bgeo,
                            std::vector<double>& bl, const std::vector<double>& a1,
                            const std::vector<double>& a2, const std::vector<double>& a3) {
        std::vector<double> mag(25 * n);
        for (std::size_t i = 0; i < n; ++i) {
            for (std::size_t k = 0; k < 25; ++k) mag[(25 * i) + k] = drv[k];
        }
        bgeo.assign(3 * n, 0.0);
        bl.assign(n, 0.0);
        std::vector<double> p1 = a1;
        std::vector<double> p2 = a2;
        std::vector<double> p3 = a3;
        int k = kext;
        g_field(&ntime, &k, options.data(), &sysaxes, iyear.data(), idoy.data(), ut.data(),
                p1.data(), p2.data(), p3.data(), mag.data(), bgeo.data(), bl.data());
    };

    // ---- 1. GET_FIELD_MULTI, and which truncation IRBEM uses ---------------------------------
    std::printf("=== 1. GET_FIELD_MULTI vs field_batch, kext=0, epoch %.1f ===\n", kEpoch);
    std::vector<double> bgeo;
    std::vector<double> bl;
    oracle_field(0, drivers(0, 0, 2, 0, 0), bgeo, bl, x1, x2, x3);
    const auto score = [&](auto model, int degree) {
        std::vector<FieldVector<Frame::GEO>> b(n);
        std::vector<double> bm(n);
        const Result<bool> r = field_batch(model, pts, b, bm);
        double worst_mag = 0.0;
        double worst_vec = 0.0;
        for (std::size_t i = 0; i < n; ++i) {
            worst_mag = std::max(worst_mag, rel(bm[i], bl[i]));
            const fx::vec3d o{bgeo[3 * i], bgeo[(3 * i) + 1], bgeo[(3 * i) + 2]};
            worst_vec = std::max(worst_vec, fx::norm(b[i].v - o) / fx::norm(o));
        }
        std::printf("  degree %2d (status %d, gpu %d): max rel |B| %.3e   max rel vector %.3e\n",
                    degree, static_cast<int>(r.status), static_cast<int>(r.value), worst_mag,
                    worst_vec);
    };
    score(*Igrf<8>::at(kEpoch), 8);
    score(*Igrf<9>::at(kEpoch), 9);
    score(*Igrf<10>::at(kEpoch), 10);
    score(*Igrf<11>::at(kEpoch), 11);
    score(*Igrf<13>::at(kEpoch), 13);

    const Igrf<10> model = *Igrf<10>::at(kEpoch);

    // ---- 2. GET_BDERIVS is a forward difference of GET_FIELD_MULTI ----------------------------
    std::printf("\n=== 2. GET_BDERIVS reconstructed from the oracle's OWN field values ===\n");
    for (const double dx : {1.0e-1, 1.0e-2, 1.0e-3}) {
        std::vector<double> mag(25 * n, 0.0);
        std::vector<double> ob(3 * n);
        std::vector<double> obm(n);
        std::vector<double> ogr(3 * n);
        std::vector<double> odf(9 * n);
        double step = dx;
        int kext = 0;
        std::vector<double> p1 = x1;
        std::vector<double> p2 = x2;
        std::vector<double> p3 = x3;
        g_bderivs(&ntime, &kext, options.data(), &sysaxes, &step, iyear.data(), idoy.data(),
                  ut.data(), p1.data(), p2.data(), p3.data(), mag.data(), ob.data(), obm.data(),
                  ogr.data(), odf.data());
        double worst_g = 0.0;
        double worst_d = 0.0;
        for (int j = 0; j < 3; ++j) {
            std::vector<double> s1 = x1;
            std::vector<double> s2 = x2;
            std::vector<double> s3 = x3;
            std::vector<double>& target = (j == 0) ? s1 : (j == 1 ? s2 : s3);
            for (double& e : target) e += dx;
            std::vector<double> sb;
            std::vector<double> sbl;
            oracle_field(0, drivers(0, 0, 2, 0, 0), sb, sbl, s1, s2, s3);
            for (std::size_t i = 0; i < n; ++i) {
                worst_g = std::max(worst_g,
                                   rel((sbl[i] - obm[i]) / dx, ogr[(3 * i) + static_cast<std::size_t>(j)]));
                for (int c = 0; c < 3; ++c) {
                    const std::size_t k = (9 * i) + static_cast<std::size_t>((3 * j) + c);
                    worst_d = std::max(worst_d,
                                       rel((sb[(3 * i) + static_cast<std::size_t>(c)] -
                                            ob[(3 * i) + static_cast<std::size_t>(c)]) / dx,
                                           odf[k]));
                }
            }
        }
        std::printf("  dX = %.0e : forward-difference reconstruction  gradBmag %.3e  diffB %.3e\n",
                    dx, worst_g, worst_d);

        // ... and field.hpp's own bderivs against the same oracle outputs.
        std::vector<BDerivatives> mine(n);
        (void)bderivs_batch(model, pts, mine, dx);
        double mg = 0.0;
        double md = 0.0;
        for (std::size_t i = 0; i < n; ++i) {
            const double sc = fx::norm(fx::vec3d{ogr[3 * i], ogr[(3 * i) + 1], ogr[(3 * i) + 2]});
            for (std::size_t j = 0; j < 3; ++j) {
                mg = std::max(mg, std::abs(mine[i].grad_b_mag[j] - ogr[(3 * i) + j]) / sc);
                for (std::size_t c = 0; c < 3; ++c) {
                    md = std::max(md, std::abs(mine[i].diff_b(c, j) -
                                               odf[(9 * i) + (3 * j) + c]) / sc);
                }
            }
        }
        std::printf("             field.hpp bderivs_batch vs oracle  gradBmag %.3e  diffB %.3e\n",
                    mg, md);
    }

    // ---- 3. COMPUTE_GRAD_CURV_CURL under storm conditions -------------------------------------
    std::printf("\n=== 3. COMPUTE_GRAD_CURV_CURL: field.hpp's algebra vs the oracle's, in storms ===\n");
    std::vector<Storm> storms;
    for (double kp = 0; kp <= 9; kp += 1.0) {
        storms.push_back({"T89 Kp", 4, drivers(kp * 10.0, -20.0, 2.0, 0.0, 0.0)});
    }
    for (double dst : {0.0, -50.0, -100.0, -200.0, -300.0, -400.0}) {
        storms.push_back({"T96 Dst", 7, drivers(50.0, dst, 4.0, 0.0, -5.0)});
        storms.push_back({"T01s Dst", 10, drivers(50.0, dst, 4.0, 0.0, -5.0)});
    }
    for (double pdyn : {0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 40.0}) {
        storms.push_back({"T01s Pdyn", 10, drivers(60.0, -150.0, pdyn, 0.0, -10.0)});
    }
    for (double bz = 0.0; bz >= -30.0; bz -= 2.0) {
        storms.push_back({"T01s Bz south", 10, drivers(70.0, -200.0, 8.0, 3.0, bz)});
    }
    storms.push_back({"quiet internal", 0, drivers(0.0, 0.0, 2.0, 0.0, 0.0)});

    std::array<double, 8> worst{};
    const char* names[8] = {"grad_par",  "grad_perp", "grad_drift", "curvature",
                            "r_curv",    "curv_drift", "curl_b",    "div_b"};
    std::size_t compared = 0;
    for (const Storm& s : storms) {
        std::vector<double> mag(25 * n);
        for (std::size_t i = 0; i < n; ++i) {
            for (std::size_t k = 0; k < 25; ++k) mag[(25 * i) + k] = s.maginput[k];
        }
        std::vector<double> ob(3 * n);
        std::vector<double> obm(n);
        std::vector<double> ogr(3 * n);
        std::vector<double> odf(9 * n);
        double step = 1.0e-3;
        int kext = s.kext;
        std::vector<double> p1 = x1;
        std::vector<double> p2 = x2;
        std::vector<double> p3 = x3;
        g_bderivs(&ntime, &kext, options.data(), &sysaxes, &step, iyear.data(), idoy.data(),
                  ut.data(), p1.data(), p2.data(), p3.data(), mag.data(), ob.data(), obm.data(),
                  ogr.data(), odf.data());

        std::vector<double> gp(n);
        std::vector<double> gperp(3 * n);
        std::vector<double> gdrift(3 * n);
        std::vector<double> curv(3 * n);
        std::vector<double> rc(n);
        std::vector<double> cdrift(3 * n);
        std::vector<double> curlb(3 * n);
        std::vector<double> divb(n);
        g_gcc(&ntime, ob.data(), obm.data(), ogr.data(), odf.data(), gp.data(), gperp.data(),
              gdrift.data(), curv.data(), rc.data(), cdrift.data(), curlb.data(), divb.data());

        for (std::size_t i = 0; i < n; ++i) {
            if (!(obm[i] > 0.0) || obm[i] > 1e30) continue;   // the oracle's baddata sentinel
            BDerivatives d{};
            d.b = FieldVector<Frame::GEO>{fx::vec3d{ob[3 * i], ob[(3 * i) + 1], ob[(3 * i) + 2]}};
            d.b_mag = obm[i];
            d.grad_b_mag = fx::vec3d{ogr[3 * i], ogr[(3 * i) + 1], ogr[(3 * i) + 2]};
            for (std::size_t j = 0; j < 3; ++j) {
                for (std::size_t c = 0; c < 3; ++c) d.diff_b(c, j) = odf[(9 * i) + (3 * j) + c];
            }
            const Result<GradCurvCurl> g = grad_curv_curl(d);
            if (g.status != Status::Ok) continue;
            ++compared;
            worst[0] = std::max(worst[0], rel(g.value.grad_par, gp[i]));
            for (std::size_t k = 0; k < 3; ++k) {
                worst[1] = std::max(worst[1], rel(g.value.grad_perp[k], gperp[(3 * i) + k]));
                worst[2] = std::max(worst[2], rel(g.value.grad_drift[k], gdrift[(3 * i) + k]));
                worst[3] = std::max(worst[3], rel(g.value.curvature[k], curv[(3 * i) + k]));
                worst[5] = std::max(worst[5], rel(g.value.curv_drift[k], cdrift[(3 * i) + k]));
                worst[6] = std::max(worst[6], rel(g.value.curl_b[k], curlb[(3 * i) + k]));
            }
            worst[4] = std::max(worst[4], rel(g.value.r_curv, rc[i]));
            worst[7] = std::max(worst[7], rel(g.value.div_b, divb[i]));
        }
    }
    std::printf("  %zu storm configurations x %zu points = %zu comparisons\n", storms.size(), n,
                compared);
    for (int k = 0; k < 8; ++k) std::printf("    %-11s max rel %.3e\n", names[k], worst[k]);

    // ---- 4. GET_HEMI_MULTI --------------------------------------------------------------------
    std::printf("\n=== 4. GET_HEMI_MULTI vs hemisphere_batch, kext=0 ===\n");
    {
        std::vector<double> mag(25 * n, 0.0);
        std::vector<int> oh(n, 0);
        int kext = 0;
        std::vector<double> p1 = x1;
        std::vector<double> p2 = x2;
        std::vector<double> p3 = x3;
        g_hemi(&ntime, &kext, options.data(), &sysaxes, iyear.data(), idoy.data(), ut.data(),
               p1.data(), p2.data(), p3.data(), mag.data(), oh.data());
        std::vector<Hemisphere> mine(n);
        (void)hemisphere_batch(model, pts, mine);
        std::size_t agree = 0;
        for (std::size_t i = 0; i < n; ++i) {
            if (static_cast<int>(mine[i]) == oh[i]) ++agree;
        }
        std::printf("  agree %zu / %zu\n", agree, n);
    }

    // ---- 5. golden vectors for the unit test --------------------------------------------------
    std::printf("\n=== 5. golden vectors (kext=0, epoch %.1f, degree 10, dX = 1e-3) ===\n", kEpoch);
    const std::array<fx::vec3d, 4> golden{fx::vec3d{2.0, 0.0, 0.0}, fx::vec3d{3.0, 1.0, 0.5},
                                          fx::vec3d{-4.0, 2.0, 1.5}, fx::vec3d{5.5, 0.0, 2.0}};
    {
        const int m = static_cast<int>(golden.size());
        std::vector<double> a1(m);
        std::vector<double> a2(m);
        std::vector<double> a3(m);
        for (int i = 0; i < m; ++i) {
            a1[i] = golden[static_cast<std::size_t>(i)][0];
            a2[i] = golden[static_cast<std::size_t>(i)][1];
            a3[i] = golden[static_cast<std::size_t>(i)][2];
        }
        std::vector<int> yy(m, 2015);
        std::vector<int> dd(m, 180);
        std::vector<double> tt(m, 43200.0);
        std::vector<double> mag(25 * m, 0.0);
        std::vector<double> ob(3 * m);
        std::vector<double> obm(m);
        std::vector<double> ogr(3 * m);
        std::vector<double> odf(9 * m);
        std::vector<int> oh(m);
        int mt = m;
        int kext = 0;
        double step = 1.0e-3;
        std::vector<double> q1 = a1;
        std::vector<double> q2 = a2;
        std::vector<double> q3 = a3;
        g_bderivs(&mt, &kext, options.data(), &sysaxes, &step, yy.data(), dd.data(), tt.data(),
                  q1.data(), q2.data(), q3.data(), mag.data(), ob.data(), obm.data(), ogr.data(),
                  odf.data());
        q1 = a1;
        q2 = a2;
        q3 = a3;
        g_hemi(&mt, &kext, options.data(), &sysaxes, yy.data(), dd.data(), tt.data(), q1.data(),
               q2.data(), q3.data(), mag.data(), oh.data());
        std::vector<double> gp(m);
        std::vector<double> gperp(3 * m);
        std::vector<double> gdrift(3 * m);
        std::vector<double> curv(3 * m);
        std::vector<double> rc(m);
        std::vector<double> cdrift(3 * m);
        std::vector<double> curlb(3 * m);
        std::vector<double> divb(m);
        g_gcc(&mt, ob.data(), obm.data(), ogr.data(), odf.data(), gp.data(), gperp.data(),
              gdrift.data(), curv.data(), rc.data(), cdrift.data(), curlb.data(), divb.data());
        for (int i = 0; i < m; ++i) {
            std::printf("  {%.1f, %.1f, %.1f}\n", a1[i], a2[i], a3[i]);
            std::printf("    B      = {%.9g, %.9g, %.9g}   |B| = %.10g\n", ob[3 * i],
                        ob[(3 * i) + 1], ob[(3 * i) + 2], obm[i]);
            std::printf("    gradB  = {%.9g, %.9g, %.9g}\n", ogr[3 * i], ogr[(3 * i) + 1],
                        ogr[(3 * i) + 2]);
            std::printf("    diffB  = {");
            for (int k = 0; k < 9; ++k) std::printf("%.9g%s", odf[(9 * i) + k], k == 8 ? "}\n" : ", ");
            std::printf("    gpar %.9g  rcurv %.9g  divB %.9g  hemi %+d\n", gp[i], rc[i], divb[i],
                        oh[i]);
            std::printf("    curl   = {%.9g, %.9g, %.9g}\n", curlb[3 * i], curlb[(3 * i) + 1],
                        curlb[(3 * i) + 2]);
            std::printf("    curvat = {%.9g, %.9g, %.9g}\n", curv[3 * i], curv[(3 * i) + 1],
                        curv[(3 * i) + 2]);
            // grad_perp, grad_drift and curv_drift were absent from this print until a review
            // showed the consequence: with no golden for them, DOUBLING any of the three still
            // passed the whole unit suite, because the only assertions on them were the
            // orthogonality relations, which a scale factor leaves intact. They are printed here
            // and frozen in tests/irbem_field_test.cpp so a wrong 1/Bmag can no longer hide.
            std::printf("    gperp  = {%.9g, %.9g, %.9g}\n", gperp[3 * i], gperp[(3 * i) + 1],
                        gperp[(3 * i) + 2]);
            std::printf("    gdrift = {%.9g, %.9g, %.9g}\n", gdrift[3 * i], gdrift[(3 * i) + 1],
                        gdrift[(3 * i) + 2]);
            std::printf("    cdrift = {%.9g, %.9g, %.9g}\n", cdrift[3 * i], cdrift[(3 * i) + 1],
                        cdrift[(3 * i) + 2]);
        }
    }
    return 0;
}
