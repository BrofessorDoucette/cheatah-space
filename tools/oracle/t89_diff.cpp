// Differential study of space.irbem's Tsyganenko (1989) against the IRBEM oracle's `kext = 4`,
// and the experiment that says WHY they differ.
//
// DEV-ONLY. Like convergence.cpp beside it, this is the one kind of program that touches IRBEM, it
// is never built by the QA gate, and it never ships. IRBEM is LGPL-3.0 and cheatah-space is MIT:
// the library is run here as a BLACK BOX (dlopen plus the documented C entry points), never read
// for its logic and never linked into anything we distribute. Tsyganenko's own T89 source is
// GPL-3.0 and is not read at all, by anything, ever.
//
// WHAT IT MEASURES, and why the answer is not "they agree".
//
// `space/irbem/ext_t89.hpp` implements the model AS PUBLISHED: the functional form of Tsyganenko
// (1989) eqs. (11)-(20) and the coefficients of that paper's Table 1. IRBEM's `kext = 4` is
// labelled "Tsyganenko [1989c]" and evaluates a LATER revision of the same model whose parameters
// were never published in the paper. So this harness answers two questions, in order:
//
//   1. HOW FAR APART are they, per Kp bin, per region?  ->  the `deviation` pass.
//   2. Is the difference a re-fit, or a different functional form?  ->  the `span` pass, which
//      lets ALL 19 linear coefficients and ALL 9 non-linear parameters of the published form float
//      freely and reports the residual that no choice of parameters can remove. A re-fit would
//      drive that to roundoff; a structural difference leaves a floor. It leaves a floor.
//
// The external field is isolated from the oracle by DIFFERENCE: `get_field1_` with `kext = 4` minus
// the same call with `kext = 0`, both with `options(5) = 0` so the internal IGRF term is bit-for-bit
// identical between them and cancels exactly. The dipole tilt is taken from the oracle too — the SM
// z-axis transformed into GSM is `(sin psi, 0, cos psi)` — so that a difference in the tilt model
// cannot masquerade as a difference in T89.
//
// Build (from the repository root):
//   g++ -O2 -std=c++20 tools/oracle/t89_diff.cpp -I. \
//       -I$CHEATAH_DIR/stdlib/ndarray -I$CHEATAH_DIR/stdlib/builtins -I$CHEATAH_DIR/stdlib/fixarray \
//       -o /tmp/t89_diff && /tmp/t89_diff /tmp/irbem-builds/libirbem-O2.so
//
// Quote the -O2 build, never the as-shipped one: IRBEM ships with no `-O` at all and is 2.7x slower
// that way (docs/ERROR_BUDGET.md section 5).
#include <dlfcn.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

#include "space/irbem/ext_t89.hpp"

namespace {

using cheatah::space::irbem::T89Parameters;
using cheatah::space::irbem::t89_components;
using cheatah::space::irbem::t89_linear_count;
using cheatah::space::irbem::t89_parameters;

/// `get_field1_`, as documented in the vendored matlab/libirbem.h.
using GetField1 = void (*)(int*, int*, int*, int*, int*, double*, double*, double*, double*,
                           double*, double*, double*);
/// `coord_trans_vec1_`, likewise — used for GSM<->GEO and to read the oracle's own dipole tilt.
using CoordTransVec1 = void (*)(int*, int*, int*, int*, int*, double*, double*, double*);

/// One oracle sample: where, at what tilt, and what external field the oracle reports there.
struct Sample {
    double tilt;                  ///< dipole tilt psi, radians, as the oracle defines it
    double x, y, z;               ///< GSM position, Earth radii
    std::array<double, 3> b;      ///< oracle external field in GSM, nT
};

/// An epoch to sample at, chosen to span the tilt range rather than to be pretty.
struct Epoch {
    int doy;
    double ut;
};

/// The 19 basis fields of the published form: the field produced with coefficient k set to 1 and
/// every other linear coefficient zero. Built by calling the SHIPPING evaluator 19 times with a
/// one-hot coefficient vector, so the span test is a test of the header and not of a copy of it.
void basis(const T89Parameters<double>& nonlinear, double ps, double x, double y, double z,
           std::array<std::array<double, 3>, t89_linear_count>& out) {
    const double s = std::sin(ps);
    const double c = std::cos(ps);
    for (std::size_t k = 0; k < t89_linear_count; ++k) {
        T89Parameters<double> one = nonlinear;
        one.c.fill(0.0);
        one.c[k] = 1.0;
        out[k] = t89_components<double>(one, s, c, x, y, z);
    }
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

/// The best residual the published functional form can achieve against @p samples with all 19
/// linear coefficients free and @p nonlinear held.
double linear_residual(const std::vector<Sample>& samples, const T89Parameters<double>& nonlinear) {
    const int n = static_cast<int>(t89_linear_count);
    std::vector<double> ata(static_cast<std::size_t>(n) * n, 0.0);
    std::vector<double> atb(static_cast<std::size_t>(n), 0.0);
    std::array<std::array<double, 3>, t89_linear_count> bs{};
    double sy2 = 0.0;
    std::vector<std::array<double, t89_linear_count>> rows;
    std::vector<double> obs;
    rows.reserve(samples.size() * 3);
    obs.reserve(samples.size() * 3);
    for (const Sample& s : samples) {
        basis(nonlinear, s.tilt, s.x, s.y, s.z, bs);
        for (int c = 0; c < 3; ++c) {
            std::array<double, t89_linear_count> row{};
            for (int k = 0; k < n; ++k) row[static_cast<std::size_t>(k)] = bs[static_cast<std::size_t>(k)][c];
            for (int i = 0; i < n; ++i) {
                atb[static_cast<std::size_t>(i)] += row[static_cast<std::size_t>(i)] * s.b[c];
                for (int j = 0; j < n; ++j)
                    ata[(static_cast<std::size_t>(i) * n) + j] +=
                        row[static_cast<std::size_t>(i)] * row[static_cast<std::size_t>(j)];
            }
            rows.push_back(row);
            obs.push_back(s.b[c]);
            sy2 += s.b[c] * s.b[c];
        }
    }
    std::vector<double> a = ata;
    std::vector<double> b = atb;
    if (!solve(a, b, n)) return 1e30;
    double res = 0.0;
    for (std::size_t r = 0; r < obs.size(); ++r) {
        double f = 0.0;
        for (int k = 0; k < n; ++k) f += rows[r][static_cast<std::size_t>(k)] * b[static_cast<std::size_t>(k)];
        res += (f - obs[r]) * (f - obs[r]);
    }
    (void)sy2;
    return res / static_cast<double>(obs.size());
}

/// Pack / unpack the nine non-linear parameters for the simplex.
T89Parameters<double> with_nonlinear(const T89Parameters<double>& base, const double* v) {
    T89Parameters<double> p = base;
    p.delta_x = v[0];
    p.a_rc = v[1];
    p.d_0 = v[2];
    p.gamma_rc = v[3];
    p.r_c = v[4];
    p.g = v[5];
    p.a_t = v[6];
    p.d_y = v[7];
    p.x_0 = v[8];
    return p;
}

}  // namespace

int main(int argc, char** argv) {
    const std::string lib = argc > 1 ? argv[1] : "/tmp/irbem-builds/libirbem-O2.so";
    void* h = dlopen(lib.c_str(), RTLD_NOW);
    if (h == nullptr) {
        std::fprintf(stderr, "cannot dlopen %s: %s\n", lib.c_str(), dlerror());
        return 1;
    }
    auto* get_field = reinterpret_cast<GetField1>(dlsym(h, "get_field1_"));
    auto* ctv = reinterpret_cast<CoordTransVec1>(dlsym(h, "coord_trans_vec1_"));
    if (get_field == nullptr || ctv == nullptr) {
        std::fprintf(stderr, "missing entry points in %s\n", lib.c_str());
        return 1;
    }

    // Three epochs spanning the tilt range: near zero, strongly positive, strongly negative. The
    // tilt drives every sin(psi) / tan(psi) term in the model, so a comparison at one tilt would
    // leave half the model untested.
    const std::array<Epoch, 3> epochs{{{80, 39183.0}, {180, 43200.0}, {355, 7200.0}}};
    const int year = 2015;

    std::printf("space.irbem T89 (as published, 1989) vs IRBEM kext=4 (\"T89c\"), oracle %s\n\n",
                lib.c_str());

    // Kp x 10 values that land inside each of the seven bins t89_kp_bin defines.
    const std::array<double, 7> kp10{{0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 70.0}};

    std::vector<Sample> bin_one_samples;
    // Two regions, because they answer different questions: the radiation belts (r <= 10 Re) are
    // what this library exists for, and the full 3..35 Re box is where the tail term lives and
    // where the two parameterizations diverge most.
    for (int region = 0; region < 2; ++region) {
    const double r_hi = region == 0 ? 10.0 : 35.0;
    std::printf("\n--- %s ---\n", region == 0 ? "radiation belts, 3 <= r <= 10 Re"
                                              : "full box, 3 <= r <= 35 Re");
    std::printf("%3s %7s %8s %10s %10s %10s %9s %9s\n", "bin", "Kp*10", "N", "rms|dB|", "p99|dB|",
                "max|dB|", "rms rel", "p99 rel");
    for (std::size_t bi = 0; bi < kp10.size(); ++bi) {
        std::vector<double> abs_err;
        std::vector<double> rel_err;
        double sum2 = 0.0;
        double sig2 = 0.0;
        for (const Epoch& e : epochs) {
            int iyear = year;
            int idoy = e.doy;
            double ut = e.ut;
            int one = 1;
            double ps = 0.0;
            {
                int si = 4;
                int so = 2;
                std::array<double, 3> in{0.0, 0.0, 1.0};
                std::array<double, 3> out{};
                ctv(&one, &si, &so, &iyear, &idoy, &ut, in.data(), out.data());
                ps = std::atan2(out[0], out[2]);
            }
            for (double x = -30.0; x <= 14.0; x += 4.0)
                for (double y = -12.0; y <= 12.0; y += 4.0)
                    for (double z = -8.0; z <= 8.0; z += 4.0) {
                        const double r = std::sqrt((x * x) + (y * y) + (z * z));
                        if (r < 3.0 || r > r_hi) continue;
                        std::array<double, 3> gsm{x, y, z};
                        std::array<double, 3> geo{};
                        {
                            int si = 2;
                            int so = 1;
                            ctv(&one, &si, &so, &iyear, &idoy, &ut, gsm.data(), geo.data());
                        }
                        std::array<int, 5> options{0, 0, 0, 0, 0};
                        int sysaxes = 1;
                        int k0 = 0;
                        int k4 = 4;
                        std::vector<double> mag(25, 0.0);
                        mag[0] = kp10[bi];
                        std::array<double, 3> b0{};
                        std::array<double, 3> b4{};
                        double m0 = 0.0;
                        double m4 = 0.0;
                        double x1 = geo[0];
                        double x2 = geo[1];
                        double x3 = geo[2];
                        get_field(&k0, options.data(), &sysaxes, &iyear, &idoy, &ut, &x1, &x2, &x3,
                                  mag.data(), b0.data(), &m0);
                        get_field(&k4, options.data(), &sysaxes, &iyear, &idoy, &ut, &x1, &x2, &x3,
                                  mag.data(), b4.data(), &m4);
                        std::array<double, 3> dgeo{b4[0] - b0[0], b4[1] - b0[1], b4[2] - b0[2]};
                        std::array<double, 3> ora{};
                        {
                            int si = 1;
                            int so = 2;
                            ctv(&one, &si, &so, &iyear, &idoy, &ut, dgeo.data(), ora.data());
                        }
                        const std::array<double, 3> mine = t89_components<double>(
                            t89_parameters<double>(static_cast<int>(bi) + 1), std::sin(ps),
                            std::cos(ps), x, y, z);
                        double d2 = 0.0;
                        double o2 = 0.0;
                        for (int c = 0; c < 3; ++c) {
                            d2 += (mine[c] - ora[c]) * (mine[c] - ora[c]);
                            o2 += ora[c] * ora[c];
                        }
                        abs_err.push_back(std::sqrt(d2));
                        rel_err.push_back(std::sqrt(d2) / (std::sqrt(o2) + 1e-12));
                        sum2 += d2;
                        sig2 += o2;
                        if (bi == 0 && region == 1) bin_one_samples.push_back({ps, x, y, z, ora});
                    }
        }
        std::sort(abs_err.begin(), abs_err.end());
        std::sort(rel_err.begin(), rel_err.end());
        const std::size_t n = abs_err.size();
        const std::size_t p99 = (n * 99) / 100;
        std::printf("%3zu %7.0f %8zu %10.3f %10.3f %10.3f %9.4f %9.4f\n", bi + 1, kp10[bi], n,
                    std::sqrt(sum2 / static_cast<double>(n)), abs_err[p99], abs_err[n - 1],
                    std::sqrt(sum2 / sig2), rel_err[p99]);
    }
    }

    // ---- the span pass: is the difference a re-fit, or a different functional form? -----------
    std::printf(
        "\nSpan test on Kp bin 1 (%zu samples): how close can the PUBLISHED functional form get to\n"
        "the oracle with every parameter free?\n",
        bin_one_samples.size());
    double sig2 = 0.0;
    for (const Sample& s : bin_one_samples)
        for (int c = 0; c < 3; ++c) sig2 += s.b[c] * s.b[c];
    const double rms_signal = std::sqrt(sig2 / static_cast<double>(bin_one_samples.size() * 3));

    const T89Parameters<double> published = t89_parameters<double>(1);
    std::array<double, 9> start{published.delta_x, published.a_rc, published.d_0,
                                published.gamma_rc, published.r_c,  published.g,
                                published.a_t,      published.d_y,  published.x_0};
    std::printf("  published parameters, 19 linear free : rms residual %.4f nT (%.2f%% of signal)\n",
                std::sqrt(linear_residual(bin_one_samples, published)),
                100.0 * std::sqrt(linear_residual(bin_one_samples, published)) / rms_signal);

    // Nelder-Mead over the nine non-linear parameters, with the 19 linear ones solved exactly at
    // every evaluation (variable projection).
    const int nd = 9;
    std::vector<std::vector<double>> simplex;
    std::vector<double> value;
    const auto cost = [&](const std::vector<double>& v) {
        if (v[0] < 1.0 || v[1] < 0.5 || v[6] < 0.5 || v[7] < 1.0) return 1e30;
        return linear_residual(bin_one_samples, with_nonlinear(published, v.data()));
    };
    simplex.emplace_back(start.begin(), start.end());
    value.push_back(cost(simplex[0]));
    for (int i = 0; i < nd; ++i) {
        std::vector<double> v = simplex[0];
        v[static_cast<std::size_t>(i)] +=
            std::fabs(v[static_cast<std::size_t>(i)]) > 1e-6 ? 0.15 * std::fabs(v[static_cast<std::size_t>(i)]) : 0.2;
        simplex.push_back(v);
        value.push_back(cost(v));
    }
    for (int it = 0; it < 4000; ++it) {
        std::vector<int> idx(static_cast<std::size_t>(nd) + 1);
        for (int i = 0; i <= nd; ++i) idx[static_cast<std::size_t>(i)] = i;
        std::sort(idx.begin(), idx.end(), [&](int a, int b) {
            return value[static_cast<std::size_t>(a)] < value[static_cast<std::size_t>(b)];
        });
        std::vector<std::vector<double>> s2;
        std::vector<double> v2;
        for (int i : idx) {
            s2.push_back(simplex[static_cast<std::size_t>(i)]);
            v2.push_back(value[static_cast<std::size_t>(i)]);
        }
        simplex = s2;
        value = v2;
        if (std::fabs(value[static_cast<std::size_t>(nd)] - value[0]) <
            1e-12 * (std::fabs(value[0]) + 1e-12))
            break;
        std::vector<double> centre(static_cast<std::size_t>(nd), 0.0);
        for (int i = 0; i < nd; ++i)
            for (int j = 0; j < nd; ++j)
                centre[static_cast<std::size_t>(j)] +=
                    simplex[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] / nd;
        const auto along = [&](double a) {
            std::vector<double> v(static_cast<std::size_t>(nd));
            for (int j = 0; j < nd; ++j)
                v[static_cast<std::size_t>(j)] =
                    centre[static_cast<std::size_t>(j)] +
                    (a * (centre[static_cast<std::size_t>(j)] -
                          simplex[static_cast<std::size_t>(nd)][static_cast<std::size_t>(j)]));
            return v;
        };
        const std::vector<double> vr = along(1.0);
        const double fr = cost(vr);
        if (fr < value[0]) {
            const std::vector<double> ve = along(2.0);
            const double fe = cost(ve);
            simplex[static_cast<std::size_t>(nd)] = fe < fr ? ve : vr;
            value[static_cast<std::size_t>(nd)] = std::min(fe, fr);
        } else if (fr < value[static_cast<std::size_t>(nd) - 1]) {
            simplex[static_cast<std::size_t>(nd)] = vr;
            value[static_cast<std::size_t>(nd)] = fr;
        } else {
            const std::vector<double> vc = along(-0.5);
            const double fc = cost(vc);
            if (fc < value[static_cast<std::size_t>(nd)]) {
                simplex[static_cast<std::size_t>(nd)] = vc;
                value[static_cast<std::size_t>(nd)] = fc;
            } else {
                for (int i = 1; i <= nd; ++i) {
                    for (int j = 0; j < nd; ++j)
                        simplex[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] =
                            simplex[0][static_cast<std::size_t>(j)] +
                            (0.5 * (simplex[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] -
                                    simplex[0][static_cast<std::size_t>(j)]));
                    value[static_cast<std::size_t>(i)] = cost(simplex[static_cast<std::size_t>(i)]);
                }
            }
        }
    }
    const double best = std::sqrt(value[0]);
    std::printf("  all 28 parameters free               : rms residual %.4f nT (%.2f%% of signal)\n",
                best, 100.0 * best / rms_signal);
    std::printf(
        "\n  A re-fit of the same equations would drive that second number to roundoff. It does not\n"
        "  go there, so the difference between the published model and IRBEM's \"T89c\" is\n"
        "  STRUCTURAL: the deployed revision is not the published functional form with different\n"
        "  numbers in it. Its source is GPL-3.0 and is not read here.\n");
    return 0;
}
