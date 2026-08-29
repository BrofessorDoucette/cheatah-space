// The SIMD batch lane of space.irbem — space/irbem/batch_soa.hpp.
//
// The contract under test is unusual and worth stating: the strip evaluator is a SCHEDULING change
// only. Every assertion of value equality here is therefore memcmp — bit identity with the scalar
// kernel, not a tolerance — because a tolerance would accept precisely the reassociation and
// contraction bugs this lane must not have. The other half of the contract is that the lane stays
// VECTORISED: a regression that quietly decays it to scalar code changes no bits and fails no
// numeric test, so the emitted machine code itself is under test (StripCodegenCarriesPackedVectorOps
// disassembles this binary and counts packed mul/add/fma in the strip symbol).
//
// Perturbation controls (run during development, then reverted):
//  - re-associating the strip's reduction to `rpow * (accum(n + 1) * ap0)` — the same value,
//    rounded differently — flipped all four memcmp tests to FAIL at the first strip;
//  - replacing the vector rows with a scalar per-lane struct of the same arithmetic (a shadow copy
//    of the header, `-fno-tree-vectorize`) kept every memcmp test green — it IS the same
//    arithmetic — and flipped StripCodegenCarriesPackedVectorOps to FAIL, which is the division
//    of labour between the two halves of this suite;
//  - deleting `[[gnu::noinline]]` changed nothing observable at -O3 (GCC keeps the symbol
//    out-of-line on size grounds anyway), so the attribute is a guarantee against a future
//    inliner heuristic, not something this suite can currently distinguish.

#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <cstring>
#include <random>
#include <string>
#include <unistd.h>
#include <vector>

#include "space/irbem/batch_soa.hpp"
#include "space/irbem/field.hpp"
#include "space/irbem/igrf.hpp"

#include "alloc_counter.hpp"

namespace ib = cheatah::space::irbem;
namespace fx = cheatah::fixarray;

namespace {

/// A deterministic batch that starts with the places a dipole-ish field is least generic — both
/// poles, the equator at all four local times, off-equator points — and fills the rest with a
/// seeded pseudo-random shell sweep over r in [1, 10] Re. Bit-identity needs no exactly
/// representable inputs (both lanes read the same doubles), so the sweep favours breadth. The
/// seed is a parameter with a fixed default (the house shape, cf. sample_points in
/// irbem_field_test.cpp): reproducible, and not the constant-seed pattern cert-msc51 rejects.
std::vector<ib::Position<ib::Frame::GEO>> corpus(std::size_t n, unsigned long long seed = 20260828ULL) {
    std::vector<ib::Position<ib::Frame::GEO>> pts;
    pts.reserve(n);
    const std::array<fx::vec3d, 10> anchors{{
        fx::vec3d{0.0, 0.0, 1.5},    // north polar axis: the point the Cartesian kernel needs no branch for
        fx::vec3d{0.0, 0.0, -1.5},   // south polar axis
        fx::vec3d{1.5, 0.0, 0.0},    // equator, noon meridian
        fx::vec3d{-1.5, 0.0, 0.0},   // equator, midnight
        fx::vec3d{0.0, 1.5, 0.0},    // dusk
        fx::vec3d{0.0, -1.5, 0.0},   // dawn
        fx::vec3d{1.0, 0.0, 0.0},    // the r-validity floor of the series, on the reference sphere
        fx::vec3d{3.0, 4.0, 5.0},    // off-equator, all components active
        fx::vec3d{-6.6, 1.0, -0.5},  // geosynchronous-ish, southern
        fx::vec3d{0.25, 0.25, 9.5},  // near-axis at large r
    }};
    for (const fx::vec3d& v : anchors) {
        if (pts.size() == n) break;
        pts.push_back(ib::Position<ib::Frame::GEO>{v});
    }
    std::mt19937_64 rng(seed);  // a fixed default seed, as sample_points in irbem_field_test.cpp
    std::uniform_real_distribution<double> radius(1.0, 10.0);
    std::uniform_real_distribution<double> comp(-1.0, 1.0);
    while (pts.size() < n) {
        const double x = comp(rng);
        const double y = comp(rng);
        const double z = comp(rng);
        const double s = std::sqrt((x * x) + (y * y) + (z * z));
        if (!(s > 0.1)) continue;  // reject near-origin directions rather than dividing by ~0
        const double r = radius(rng);
        pts.push_back(ib::Position<ib::Frame::GEO>{fx::vec3d{r * x / s, r * y / s, r * z / s}});
    }
    return pts;
}

/// The scalar reference: the exact per-point loop the batch lane replaced.
template <int NMAX, ib::SoundPrecision P>
void scalar_reference(const ib::Igrf<NMAX, P>& model,
                      const std::vector<ib::Position<ib::Frame::GEO>>& pts,
                      std::vector<ib::FieldVector<ib::Frame::GEO>>& b, std::vector<double>& mag) {
    b.resize(pts.size());
    mag.resize(pts.size());
    for (std::size_t i = 0; i < pts.size(); ++i) {
        b[i] = model.evaluate(pts[i]);
        mag[i] = b[i].magnitude();
    }
}

}  // namespace

// The core contract: over a large batch spanning both truncations, two epochs, poles, all local
// times and 1-10 Re, every byte of every field vector and every magnitude is identical to the
// scalar kernel's. 4099 points is deliberately 3 (mod 8): the last strip boundary and a non-empty
// tail are inside the comparison, not a separate code path trusted on faith.
TEST(IrbemBatchSimd, BatchLaneIsBitIdenticalToTheScalarLane) {
    const auto pts = corpus(4099);
    std::vector<ib::FieldVector<ib::Frame::GEO>> want;
    std::vector<double> want_mag;
    std::vector<ib::FieldVector<ib::Frame::GEO>> got(pts.size());
    std::vector<double> got_mag(pts.size());

    {
        const auto model = ib::Igrf<13>::at(2015.5);
        ASSERT_TRUE(model.has_value());
        scalar_reference(*model, pts, want, want_mag);
        ASSERT_EQ(ib::Status::Ok, ib::igrf_batch_host(*model, std::span(pts), std::span(got),
                                                      std::span(got_mag)));
        EXPECT_EQ(0, std::memcmp(want.data(), got.data(), pts.size() * sizeof(want[0])));
        EXPECT_EQ(0, std::memcmp(want_mag.data(), got_mag.data(), pts.size() * sizeof(double)));
    }
    {
        // The other truncation, at an interpolated historical epoch: a second complete
        // instantiation of the strip template, same bits.
        const auto model = ib::Igrf<10>::at(1987.25);
        ASSERT_TRUE(model.has_value());
        scalar_reference(*model, pts, want, want_mag);
        ASSERT_EQ(ib::Status::Ok, ib::igrf_batch_host(*model, std::span(pts), std::span(got),
                                                      std::span(got_mag)));
        EXPECT_EQ(0, std::memcmp(want.data(), got.data(), pts.size() * sizeof(want[0])));
        EXPECT_EQ(0, std::memcmp(want_mag.data(), got_mag.data(), pts.size() * sizeof(double)));
    }
}

// The fp32-integrand policy drives the OTHER LaneVec specialisation and the integrand-to-
// accumulator re-rowing (8-lane float rows carried into pairs of 4-lane double rows). Same
// contract, same memcmp.
TEST(IrbemBatchSimd, FastPolicyStripIsBitIdenticalToo) {
    static_assert(ib::igrf_strip_points<ib::Fast> == 16);
    const auto pts = corpus(1013);  // 5 (mod 16): strips plus a tail
    const auto model = ib::Igrf<10, ib::Fast>::at(2003.0);
    ASSERT_TRUE(model.has_value());
    std::vector<ib::FieldVector<ib::Frame::GEO>> want;
    std::vector<double> want_mag;
    scalar_reference(*model, pts, want, want_mag);
    std::vector<ib::FieldVector<ib::Frame::GEO>> got(pts.size());
    std::vector<double> got_mag(pts.size());
    ASSERT_EQ(ib::Status::Ok,
              ib::igrf_batch_host(*model, std::span(pts), std::span(got), std::span(got_mag)));
    EXPECT_EQ(0, std::memcmp(want.data(), got.data(), pts.size() * sizeof(want[0])));
    EXPECT_EQ(0, std::memcmp(want_mag.data(), got_mag.data(), pts.size() * sizeof(double)));
}

// Every batch length from empty through two whole strips: each length exercises a different
// strip/tail split, and the empty batch must be Ok-and-untouched rather than a special case.
TEST(IrbemBatchSimd, EveryTailLengthIsBitIdentical) {
    static_assert(ib::igrf_strip_points<ib::Exact> == 8);
    constexpr std::size_t kMax = 2 * ib::igrf_strip_points<ib::Exact>;
    const auto model = ib::Igrf<13>::at(2020.0);
    ASSERT_TRUE(model.has_value());
    const auto all = corpus(kMax);
    for (std::size_t n = 0; n <= kMax; ++n) {
        const std::vector<ib::Position<ib::Frame::GEO>> pts(all.begin(),
                                                            all.begin() + static_cast<long>(n));
        std::vector<ib::FieldVector<ib::Frame::GEO>> want;
        std::vector<double> want_mag;
        scalar_reference(*model, pts, want, want_mag);
        std::vector<ib::FieldVector<ib::Frame::GEO>> got(n);
        std::vector<double> got_mag(n);
        ASSERT_EQ(ib::Status::Ok,
                  ib::igrf_batch_host(*model, std::span(pts), std::span(got), std::span(got_mag)))
            << "length " << n;
        if (n == 0) continue;  // an empty vector's data() is null, and memcmp's pointers are nonnull
        EXPECT_EQ(0, std::memcmp(want.data(), got.data(), n * sizeof(want[0]))) << "length " << n;
        EXPECT_EQ(0, std::memcmp(want_mag.data(), got_mag.data(), n * sizeof(double)))
            << "length " << n;
    }
}

// field_batch's host lane IS this lane now: the same memcmp, through the public entry point the
// rest of the library calls. This TU is compiled without the GPU seam, so the host path is the
// only path and the assertion cannot be satisfied by a device answer.
TEST(IrbemBatchSimd, FieldBatchHostLaneIsBitIdentical) {
    const auto pts = corpus(1027);
    const auto model = ib::Igrf<13>::at(2015.5);
    ASSERT_TRUE(model.has_value());
    std::vector<ib::FieldVector<ib::Frame::GEO>> want;
    std::vector<double> want_mag;
    scalar_reference(*model, pts, want, want_mag);
    std::vector<ib::FieldVector<ib::Frame::GEO>> got(pts.size());
    std::vector<double> got_mag(pts.size());
    const auto r =
        ib::field_batch(*model, std::span(pts), std::span(got), std::span(got_mag));
    ASSERT_EQ(ib::Status::Ok, r.status);
    EXPECT_FALSE(r.value);  // the host lane answered, not a device
    EXPECT_EQ(0, std::memcmp(want.data(), got.data(), pts.size() * sizeof(want[0])));
    EXPECT_EQ(0, std::memcmp(want_mag.data(), got_mag.data(), pts.size() * sizeof(double)));
}

// The refusal branch: mismatched span lengths are a caller bug, answered with DomainError before
// anything is written.
TEST(IrbemBatchSimd, BatchRefusesMismatchedSpans) {
    const auto model = ib::Igrf<13>::at(2015.5);
    ASSERT_TRUE(model.has_value());
    const auto pts = corpus(8);
    std::vector<ib::FieldVector<ib::Frame::GEO>> b(7);   // one short
    std::vector<double> mag(8);
    EXPECT_EQ(ib::Status::DomainError,
              ib::igrf_batch_host(*model, std::span(pts), std::span(b), std::span(mag)));
    b.resize(8);
    mag.resize(9);  // one long
    EXPECT_EQ(ib::Status::DomainError,
              ib::igrf_batch_host(*model, std::span(pts), std::span(b), std::span(mag)));
}

// "No heap in the hot path" is a claim, so it gets the process-wide counter: the second call must
// allocate exactly nothing (the first would too — the workspace is the stack — but the second call
// is the one that catches a per-call workspace, per alloc_counter.hpp's own brief).
TEST(IrbemBatchSimd, BatchLaneNeverTouchesTheHeap) {
    const auto model = ib::Igrf<13>::at(2015.5);
    ASSERT_TRUE(model.has_value());
    const auto pts = corpus(64);
    std::vector<ib::FieldVector<ib::Frame::GEO>> b(pts.size());
    std::vector<double> mag(pts.size());
    ASSERT_EQ(ib::Status::Ok,
              ib::igrf_batch_host(*model, std::span(pts), std::span(b), std::span(mag)));
    const std::size_t before = cheatah_space_test::allocation_count();
    ASSERT_EQ(ib::Status::Ok,
              ib::igrf_batch_host(*model, std::span(pts), std::span(b), std::span(mag)));
    EXPECT_EQ(before, cheatah_space_test::allocation_count());
}

// The machine-code half of the contract. Bit-identity cannot distinguish the vector lane from the
// scalar one — that is its point — so the packed instructions are asserted directly: disassemble
// THIS binary, find the noinline strip symbol, and count packed-double mul/add/sub/fma against
// scalar ones. GCC's typed vectors emit packed instructions at every optimisation level (they are
// vector-typed IR, not an autovectorisation outcome), so the floor assertion runs even in debug
// builds; the dominance assertion (packed outnumbers scalar) is meaningful only for the optimised
// build this lane exists for, and is gated accordingly.
TEST(IrbemBatchSimd, StripCodegenCarriesPackedVectorOps) {
    // Force the instantiation whose symbol is inspected (the suite already instantiates it, but a
    // test should not depend on its neighbours for its own preconditions).
    auto* strip = &ib::detail::igrf_evaluate_strip<13, ib::Exact>;
    ASSERT_NE(nullptr, strip);

    std::array<char, 4096> exe{};
    const ssize_t len = readlink("/proc/self/exe", exe.data(), exe.size() - 1);
    ASSERT_GT(len, 0) << "cannot resolve /proc/self/exe";
    const std::string cmd =
        std::string("objdump -d --no-show-raw-insn '") + exe.data() + "' 2>/dev/null";
    FILE* pipe = popen(cmd.c_str(), "r");  // NOLINT(cert-env33-c) — fixed binary, quoted argument
    ASSERT_NE(nullptr, pipe);

    std::size_t packed = 0;
    std::size_t scalar = 0;
    bool inside = false;
    bool found = false;
    std::array<char, 1024> line{};
    while (std::fgets(line.data(), static_cast<int>(line.size()), pipe) != nullptr) {
        const std::string s(line.data());
        if (s.find(">:") != std::string::npos) {  // a new symbol starts
            inside = s.find("igrf_evaluate_strip") != std::string::npos;
            found = found || inside;
            continue;
        }
        if (!inside) continue;
        // Packed double arithmetic: mulpd/addpd/subpd (SSE2 and their AVX v- forms) and any
        // vfmadd/vfmsub/vfnmadd variant on packed doubles.
        if (s.find("mulpd") != std::string::npos || s.find("addpd") != std::string::npos ||
            s.find("subpd") != std::string::npos ||
            (s.find("vfm") != std::string::npos && s.find("pd ") != std::string::npos) ||
            (s.find("vfnm") != std::string::npos && s.find("pd ") != std::string::npos)) {
            ++packed;
        }
        if (s.find("mulsd") != std::string::npos || s.find("addsd") != std::string::npos ||
            s.find("subsd") != std::string::npos) {
            ++scalar;
        }
    }
    const int rc = pclose(pipe);
    if (rc != 0) GTEST_SKIP() << "objdump unavailable (exit " << rc << ")";
    ASSERT_TRUE(found) << "igrf_evaluate_strip symbol not found — [[gnu::noinline]] removed?";

    // The recursion alone is ~24 packed statements per loop body; anything under this floor means
    // the vector types were replaced with scalars.
    EXPECT_GE(packed, 24U) << "packed-double ops in the strip collapsed — the SIMD lane is gone";
#if defined(__OPTIMIZE__) && defined(__AVX__)
    EXPECT_GT(packed, scalar)
        << "the strip is majority-scalar again (packed " << packed << ", scalar " << scalar << ")";
#endif
}
