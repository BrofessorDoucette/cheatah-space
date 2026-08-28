// L* against the IRBEM oracle — the harness behind every number `space/irbem/driftshell.hpp`
// quotes, and the one that regenerates the golden tables `tests/irbem_driftshell_test.cpp`
// transcribes.
//
// DEV-ONLY, exactly like convergence.cpp beside it. IRBEM is LGPL-3.0 and cheatah-space is MIT, so
// the library is run here as a BLACK BOX (dlopen plus the documented C entry points), never read
// for its logic, never built by the QA gate, and never linked into anything we ship. The unit
// suite carries the numbers this program prints as constants with provenance, so the gate needs no
// oracle at all.
//
// It answers three questions, and keeps them apart because they mean different things:
//
//   1. ACCURACY at matched resolution. `docs/ERROR_BUDGET.md` §2(a) shows a 0.01 L* target is only
//      meaningful at matched `options(3)`/`options(4)`, IRBEM's own default-resolution error being
//      0.010-0.017 at L ≈ 6. So both matched settings are reported, and so is each side's distance
//      from IRBEM's most converged answer — which is the comparison that says which
//      implementation's default resolution is actually closer to the truth.
//   2. COVERAGE of the parameter space: pitch angles 15-90 degrees (shell splitting, where I ranges
//      over three orders of magnitude) and epochs 1900-2029 (where the dipole moment moves 7.8 %,
//      the control on k0 being read from the model rather than hard-coded).
//   3. COST, in milliseconds per L* point: our host lane, our device lane at several batch sizes,
//      and the oracle. The oracle is timed HERE rather than quoted, so the three numbers come off
//      the same machine in the same run.
//
// Build:
//   g++ -std=c++20 -O3 -march=native -ffp-contract=off tools/oracle/lstar_diff.cpp \
//     -I. -I$CHEATAH_DIR/stdlib/{ndarray,builtins,fixarray} [gpu includes] -ldl -o /tmp/lstar_diff
// Run:  /tmp/lstar_diff [path-to-libirbem-O2.so]
#include <dlfcn.h>

#include <array>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

#include "space/irbem/api.hpp"
#include "space/irbem/driftshell.hpp"

namespace ib = cheatah::space::irbem;
namespace fx = cheatah::fixarray;
using ib::Frame;

namespace {

/// make_lstar1_, as documented in the vendored matlab/libirbem.h.
using MakeLstar1 = void (*)(int*, int*, int*, int*, int*, int*, double*, double*, double*, double*,
                            double*, double*, double*, double*, double*, double*, double*);
/// make_lstar_shell_splitting1_ — the same, for an arbitrary pitch angle.
using SplitLstar1 = void (*)(int*, int*, int*, int*, int*, int*, int*, double*, double*, double*,
                             double*, double*, double*, double*, double*, double*, double*,
                             double*, double*);

/// IRBEM's `NALPHA_MAX`: the shell-splitting outputs are (ntime, NALPHA_MAX) whatever Nipa is.
constexpr int kAlphaMax = 25;

struct Point {
    const char* name;
    double x;
    double y;
    double z;
};

/// The twelve geometries the golden tables carry: equatorial from the inner belt through
/// geosynchronous to the outer belt, two off the Greenwich meridian, one on the diagonal, and two
/// off the equator where the line is longer and the trace works harder.
constexpr std::array<Point, 12> kPoints{{
    {"eq2.0", 2.0, 0.0, 0.0},   {"eq3.0", 3.0, 0.0, 0.0},
    {"eq4.0", 4.0, 0.0, 0.0},   {"eq5.0", 5.0, 0.0, 0.0},
    {"eq6.0", 6.0, 0.0, 0.0},   {"eq6.6", 6.6, 0.0, 0.0},
    {"eq8.0", 8.0, 0.0, 0.0},   {"y4.0", 0.0, 4.0, 0.0},
    {"y6.0", 0.0, 6.0, 0.0},
    {"xy45", 3.5355339059327378, 3.5355339059327378, 0.0},
    {"off5.5", 5.5, 0.0, 2.0},  {"off3.0", 3.0, 0.0, 1.0},
}};

struct OracleRun {
    double lm{};
    double lstar{};
    double blocal{};
    double bmin{};
    double xj{};
    double seconds{};
};

/// One `make_lstar1_` evaluation at a chosen resolution, kext = 0 (internal field only).
OracleRun run_oracle(MakeLstar1 f, const Point& p, int year, int doy, double ut, int resolution) {
    int ntime = 1;
    int kext = 0;
    int sysaxes = 1;   // GEO Cartesian, Re
    std::array<int, 5> options{1, 0, resolution, resolution, 0};
    std::array<int, 1> iyear{year};
    std::array<int, 1> idoy{doy};
    std::array<double, 1> ut_v{ut};
    std::array<double, 1> x1{p.x};
    std::array<double, 1> x2{p.y};
    std::array<double, 1> x3{p.z};
    std::vector<double> maginput(25, 0.0);
    std::array<double, 1> lm{};
    std::array<double, 1> lstar{};
    std::array<double, 1> blocal{};
    std::array<double, 1> bmin{};
    std::array<double, 1> xj{};
    std::array<double, 1> mlt{};
    const auto t0 = std::chrono::steady_clock::now();
    f(&ntime, &kext, options.data(), &sysaxes, iyear.data(), idoy.data(), ut_v.data(), x1.data(),
      x2.data(), x3.data(), maginput.data(), lm.data(), lstar.data(), blocal.data(), bmin.data(),
      xj.data(), mlt.data());
    const auto t1 = std::chrono::steady_clock::now();
    return OracleRun{lm[0],  lstar[0], blocal[0], bmin[0], xj[0],
                     std::chrono::duration<double>(t1 - t0).count()};
}

/// One `make_lstar_shell_splitting1_` evaluation at a chosen pitch angle.
double run_oracle_split(SplitLstar1 f, const Point& p, double alpha, int year, int doy, double ut,
                        int resolution) {
    int ntime = 1;
    int nipa = 1;
    int kext = 0;
    int sysaxes = 1;
    std::array<int, 5> options{1, 0, resolution, resolution, 0};
    std::array<int, 1> iyear{year};
    std::array<int, 1> idoy{doy};
    std::array<double, 1> ut_v{ut};
    std::array<double, 1> x1{p.x};
    std::array<double, 1> x2{p.y};
    std::array<double, 1> x3{p.z};
    std::array<double, 1> al{alpha};
    std::vector<double> maginput(25, 0.0);
    std::vector<double> lm(kAlphaMax, 0.0);
    std::vector<double> lstar(kAlphaMax, 0.0);
    std::vector<double> bmirr(kAlphaMax, 0.0);
    std::vector<double> bmin(kAlphaMax, 0.0);
    std::vector<double> xj(kAlphaMax, 0.0);
    std::vector<double> mlt(kAlphaMax, 0.0);
    f(&ntime, &nipa, &kext, options.data(), &sysaxes, iyear.data(), idoy.data(), ut_v.data(),
      x1.data(), x2.data(), x3.data(), al.data(), maginput.data(), lm.data(), lstar.data(),
      bmirr.data(), bmin.data(), xj.data(), mlt.data());
    return lstar[0];
}

/// Our own L*, at the same resolution, for one point.
double ours(const ib::Igrf<13>& model, const ib::Rotations& rot, const Point& p, double alpha,
            int resolution) {
    const ib::DriftShellOptions opt = ib::DriftShellOptions::from_irbem(resolution, resolution);
    const ib::Result<ib::DriftShell> s =
        ib::make_lstar(model, rot, ib::Position<Frame::GEO>{fx::vec3d{p.x, p.y, p.z}}, alpha, opt);
    return s.status == ib::Status::Ok ? s.value.lstar : std::nan("");
}

}  // namespace

int main(int argc, char** argv) {
    const std::string lib = argc > 1 ? argv[1] : "/tmp/irbem-builds/libirbem-O2.so";
    void* h = dlopen(lib.c_str(), RTLD_NOW);
    if (h == nullptr) {
        std::fprintf(stderr, "cannot dlopen %s: %s\n", lib.c_str(), dlerror());
        return 1;
    }
    auto make_lstar = reinterpret_cast<MakeLstar1>(dlsym(h, "make_lstar1_"));
    auto split = reinterpret_cast<SplitLstar1>(dlsym(h, "make_lstar_shell_splitting1_"));
    if (make_lstar == nullptr || split == nullptr) {
        std::fprintf(stderr, "missing entry points in %s\n", lib.c_str());
        return 1;
    }

    constexpr int kYear = 2015;
    constexpr int kDoy = 180;
    constexpr double kUt = 43200.0;
    const auto model = ib::Igrf<13>::at(kYear + (179.5 / 365.0));
    const ib::Result<ib::Rotations> rot = ib::api::rotations_at(kYear, kDoy, kUt, *model);

    std::printf("# space.irbem L* vs IRBEM (%s), kext=0 (IGRF internal), 2015-180 12:00 UT\n\n",
                lib.c_str());

    // --- 1. accuracy, at both matched resolutions ------------------------------------------------
    for (const int res : {0, 9}) {
        std::printf("=== matched options(3,4) = %d  (Nder = %d, dtheta = %g deg) ===\n", res,
                    25 * (res + 1), 180.0 / (720.0 * (res + 1)));
        std::printf("%-8s %14s %14s %12s %14s\n", "point", "ours", "IRBEM", "ours-IRBEM",
                    "ours-converged");
        double worst_matched = 0.0;
        double worst_converged = 0.0;
        double sum_converged = 0.0;
        double sum_oracle_own = 0.0;
        for (const Point& p : kPoints) {
            const OracleRun reference = run_oracle(make_lstar, p, kYear, kDoy, kUt, 9);
            const OracleRun matched =
                res == 9 ? reference : run_oracle(make_lstar, p, kYear, kDoy, kUt, res);
            const double mine = ours(*model, rot.value, p, 90.0, res);
            worst_matched = std::max(worst_matched, std::abs(mine - matched.lstar));
            worst_converged = std::max(worst_converged, std::abs(mine - reference.lstar));
            sum_converged += std::abs(mine - reference.lstar);
            sum_oracle_own += std::abs(matched.lstar - reference.lstar);
            std::printf("%-8s %14.9f %14.9f %+12.4f %+14.4f\n", p.name, mine, matched.lstar,
                        mine - matched.lstar, mine - reference.lstar);
        }
        const double n = static_cast<double>(kPoints.size());
        std::printf("  worst |ours - IRBEM| at this setting        : %.4f\n", worst_matched);
        std::printf("  worst |ours - converged IRBEM|             : %.4f  (mean %.4f)\n",
                    worst_converged, sum_converged / n);
        std::printf("  mean  |IRBEM at this setting - converged|  : %.4f\n\n", sum_oracle_own / n);
    }

    // --- 2. shell splitting and epochs ------------------------------------------------------------
    std::printf("=== shell splitting, matched options(3,4) = 9 ===\n");
    double worst_split = 0.0;
    for (const Point& p : {kPoints[2], kPoints[4], kPoints[10]}) {
        for (const double alpha : {90.0, 75.0, 60.0, 45.0, 30.0, 15.0}) {
            const double o = run_oracle_split(split, p, alpha, kYear, kDoy, kUt, 9);
            const double mine = ours(*model, rot.value, p, alpha, 9);
            worst_split = std::max(worst_split, std::abs(mine - o));
            std::printf("%-8s a=%5.1f  ours %12.7f  IRBEM %12.7f  %+8.4f\n", p.name, alpha, mine, o,
                        mine - o);
        }
    }
    std::printf("  worst: %.4f\n\n=== epoch sweep at (6,0,0), 90 deg, options = 9 ===\n",
                worst_split);
    double worst_epoch = 0.0;
    for (const int year : {1900, 1940, 1975, 2000, 2015, 2025, 2029}) {
        const auto m = ib::Igrf<13>::at(year + (179.5 / 365.0));
        const ib::Result<ib::Rotations> r = ib::api::rotations_at(year, kDoy, kUt, *m);
        const OracleRun o = run_oracle(make_lstar, kPoints[4], year, kDoy, kUt, 9);
        const double mine = ours(*m, r.value, kPoints[4], 90.0, 9);
        worst_epoch = std::max(worst_epoch, std::abs(mine - o.lstar));
        std::printf("%4d  k0 = %9.2f nT  ours %12.7f  IRBEM %12.7f  %+8.4f\n", year,
                    ib::dipole_moment(*m), mine, o.lstar, mine - o.lstar);
    }
    std::printf("  worst: %.4f\n\n", worst_epoch);

    // --- 3. cost -----------------------------------------------------------------------------------
    // The oracle first, at its recommended resolution, so all three numbers come off this machine.
    double oracle_seconds = 0.0;
    for (const Point& p : kPoints) oracle_seconds += run_oracle(make_lstar, p, kYear, kDoy, kUt, 0).seconds;
    std::printf("=== cost, ms per L* point, IRBEM default resolution ===\n");
    std::printf("  IRBEM %s : %8.3f\n", lib.c_str(),
                1000.0 * oracle_seconds / static_cast<double>(kPoints.size()));

    const ib::DriftShellOptions opt;
    for (const std::size_t batch : {std::size_t{1}, std::size_t{64}, std::size_t{512},
                                    std::size_t{2048}}) {
        std::vector<ib::Position<Frame::GEO>> pts(batch);
        std::vector<double> pitch(batch, 90.0);
        std::uint64_t seed = 0x9E3779B97F4A7C15ULL;
        const auto next = [&seed] {
            seed = (seed * 6364136223846793005ULL) + 1442695040888963407ULL;
            return static_cast<double>(seed >> 11U) / static_cast<double>(1ULL << 53U);
        };
        for (std::size_t i = 0; i < batch; ++i) {
            const double l = 2.0 + (6.0 * next());
            const double ph = 6.283185307179586 * next();
            pts[i] = ib::Position<Frame::GEO>{
                fx::vec3d{l * std::cos(ph), l * std::sin(ph), 0.0}};
        }
        std::vector<ib::DriftShell> out(batch);
        std::vector<ib::Status> st(batch);
        (void)ib::make_lstar_batch(*model, rot.value, pts, pitch, out, st, opt);   // warm up
        const auto t0 = std::chrono::steady_clock::now();
        const ib::Result<bool> r =
            ib::make_lstar_batch(*model, rot.value, pts, pitch, out, st, opt);
        const auto t1 = std::chrono::steady_clock::now();
        std::printf("  space.irbem batch %5zu, device=%d : %8.4f   (%d traces per point)\n", batch,
                    static_cast<int>(r.value), std::chrono::duration<double, std::milli>(t1 - t0).count() /
                                                   static_cast<double>(batch),
                    out[0].traces);
    }
    dlclose(h);
    return 0;
}
