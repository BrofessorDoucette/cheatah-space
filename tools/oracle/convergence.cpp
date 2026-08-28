// Convergence study against the IRBEM oracle — the measurement that turns space.irbem's
// ERROR_BUDGET.md from a set of plausible claims into observed numbers.
//
// DEV-ONLY. This is the one place anything links IRBEM, it is never built by the QA gate, and it
// never ships. IRBEM is LGPL-3.0 and cheatah-space is MIT: the library is run here as a BLACK BOX
// (dlopen + the documented C entry points), never read for its logic and never linked into
// anything we distribute.
//
// What it measures: IRBEM exposes two resolution knobs on L*, and its own documentation admits
// "an error of ~2% at L=6" at the recommended setting. options(3) sets the theta step used by the
// drift-shell root-find and the flux quadrature (dtheta = pi/(720*(options(3)+1))); options(4)
// sets the azimuth count (Nder = 25*(options(4)+1)). Sweeping both and watching L* move gives the
// DISCRETIZATION error directly — the term the budget claims dominates arithmetic roundoff by two
// to three orders of magnitude. If that claim survives contact with these numbers, an fp32
// integrand on the GPU is justified; if it does not, the whole precision policy changes.
//
// Build (see run.sh):
//   g++ -O2 -std=c++20 convergence.cpp -ldl -o convergence
#include <dlfcn.h>

#include <array>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

namespace {

/// IRBEM's "no data" sentinel; any output at or below this is a refusal, not a number.
constexpr double kBadData = -1.0e31;
[[nodiscard]] bool is_bad(double v) { return v <= kBadData * 0.99; }

/// make_lstar1_, as documented in the vendored matlab/libirbem.h.
using MakeLstar1 = void (*)(int*, int*, int*, int*, int*, int*, double*, double*, double*, double*,
                            double*, double*, double*, double*, double*, double*, double*);

/// One evaluation of make_lstar1 for a single point, at a chosen resolution.
struct Sample {
    double lm{}, lstar{}, blocal{}, bmin{}, xj{}, mlt{};
    double seconds{};
};

/// A test point, in GEO Cartesian (sysaxes = 1), Earth radii.
struct Point {
    const char* name;
    double x, y, z;
};

Sample run_one(MakeLstar1 make_lstar, const Point& p, int kext, int t_resol, int r_resol) {
    int ntime = 1;
    int kext_v = kext;
    int sysaxes = 1;  // GEO Cartesian, Re
    // options(1)=1 -> compute L* (0 means DO NOT, per the IRBEM options table); (2)=0 -> init IGRF
    // once per year; (3)=t_resol; (4)=r_resol; (5)=0 -> IGRF internal field.
    std::array<int, 5> options{1, 0, t_resol, r_resol, 0};
    std::array<int, 1> iyear{2015};
    std::array<int, 1> idoy{180};
    std::array<double, 1> ut{43200.0};
    std::array<double, 1> x1{p.x}, x2{p.y}, x3{p.z};
    // maginput is (25, ntime), column-major. kext 0 and 5 need none of it, but it must be valid.
    std::vector<double> maginput(25 * 1, 0.0);

    std::array<double, 1> lm{}, lstar{}, blocal{}, bmin{}, xj{}, mlt{};

    const auto t0 = std::chrono::steady_clock::now();
    make_lstar(&ntime, &kext_v, options.data(), &sysaxes, iyear.data(), idoy.data(), ut.data(),
               x1.data(), x2.data(), x3.data(), maginput.data(), lm.data(), lstar.data(),
               blocal.data(), bmin.data(), xj.data(), mlt.data());
    const auto t1 = std::chrono::steady_clock::now();

    return Sample{lm[0],  lstar[0], blocal[0],
                  bmin[0], xj[0],   mlt[0],
                  std::chrono::duration<double>(t1 - t0).count()};
}

/// Relative move from @p reference to @p value, or absolute when the reference is ~0.
[[nodiscard]] double relative_move(double value, double reference) {
    if (is_bad(value) || is_bad(reference)) return std::nan("");
    const double denom = std::abs(reference);
    return denom > 1e-12 ? std::abs(value - reference) / denom : std::abs(value - reference);
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
    if (make_lstar == nullptr) {
        std::fprintf(stderr, "no make_lstar1_ in %s\n", lib.c_str());
        return 1;
    }

    // Equatorial points spanning the inner belt through geosynchronous and out to the outer belt,
    // plus two off-equator points where the field line is longer and the trace works harder.
    const std::array<Point, 7> points{{
        {"L~2 eq", 2.0, 0.0, 0.0},
        {"L~3 eq", 3.0, 0.0, 0.0},
        {"L~4 eq", 4.0, 0.0, 0.0},
        {"L~6 eq", 6.0, 0.0, 0.0},
        {"L~6.6 GEO", 6.6, 0.0, 0.0},
        {"L~8 eq", 8.0, 0.0, 0.0},
        {"L~6 off-eq", 5.5, 0.0, 2.0},
    }};
    // kext 0 (internal field only) isolates the L* discretization from any external-model
    // parameterization; kext 5 (Olson-Pfitzer quiet) adds a realistic external field and still
    // needs no maginput, so the comparison stays clean.
    const std::array<int, 2> kexts{0, 5};

    std::printf("# IRBEM convergence study — library: %s\n", lib.c_str());
    std::printf("# options(5)=0 (IGRF internal), 2015-180 12:00 UT, GEO Cartesian\n");
    std::printf("# t_resol = options(3): dtheta = pi/(720*(t+1));  r_resol = options(4): Nder = 25*(r+1)\n\n");

    for (const int kext : kexts) {
        std::printf("=== kext=%d %s ===\n", kext,
                    kext == 0 ? "(internal field only)" : "(Olson-Pfitzer quiet)");
        std::printf("%-12s %8s %8s %14s %14s %14s %12s\n", "point", "t_resol", "r_resol", "Lstar",
                    "rel move", "XJ", "sec");

        for (const Point& p : points) {
            // The most-converged setting we will run is the reference the coarser ones move from.
            const Sample best = run_one(make_lstar, p, kext, 9, 9);

            for (const int res : {0, 1, 2, 4, 9}) {
                const Sample s = run_one(make_lstar, p, kext, res, res);
                const double move = relative_move(s.lstar, best.lstar);
                if (is_bad(s.lstar)) {
                    std::printf("%-12s %8d %8d %14s %14s %14.6g %12.4f\n", p.name, res, res,
                                "baddata", "-", s.xj, s.seconds);
                } else {
                    std::printf("%-12s %8d %8d %14.9f %14.3e %14.6g %12.4f\n", p.name, res, res,
                                s.lstar, move, s.xj, s.seconds);
                }
            }
            std::printf("\n");
        }
    }

    dlclose(h);
    return 0;
}
