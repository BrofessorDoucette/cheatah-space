// Unit tests for space.irbem's GPU seam — the kernel registry, the availability probes, the SPIR-V
// resolution, the launcher's four named failures, and the centred-dipole kernel itself.
//
// This file is compiled in TWO mutually exclusive configurations, and both are part of the gate:
//
//   (a) WITHOUT cheatah-gpu-linalg on the include path — `__has_include` is false, the header's
//       device half is preprocessed away, and every routine must still compile, still answer, and
//       report itself unavailable by name instead of crashing;
//   (b) WITH the consumer seam from cheatah-gpu-linalg/build/consumer.cmake, where the RTX
//       actually runs the kernel.
//
// Every test here is written to be meaningful in both. The device-only tests ask `available()` and
// GTEST_SKIP with `unavailable_reason()` — and a skip is reported as a skip, never as a pass.
//
// Arithmetic: the dipole formula is a ratio of polynomials, so there ARE exactly-representable
// cases, and they are used rather than a tolerance. The 3-4-5 triangle gives r = 5 and r^5 = 3125
// exactly; feeding g10 = 3125 nT makes every intermediate a small integer and the answer an exact
// integer in both binary32 and binary64 (Bx = 3*3125*3*4/3125 = 36, Bz = 3125*(3*16-25)/3125 = 23),
// so host and device are compared with `==`. A tolerance is used only where a genuine irrational
// enters — the pseudo-random lattice — and there it is the repo's own `Blocal` budget.
#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <set>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#include "space/irbem/gpu/dispatch.hpp"

namespace ib = cheatah::space::irbem;
namespace ig = cheatah::space::irbem::gpu;
namespace fx = cheatah::fixarray;

using ib::Frame;

// ---- test-local helpers -------------------------------------------------------------------

namespace {

/// Sets an environment variable for the life of the object and puts the old value (or its
/// absence) back afterwards, so no test can leak an override into the next one.
class EnvGuard {
public:
    /// @param name the variable. @param value the value, or nullptr to unset it.
    EnvGuard(const char* name, const char* value) : name_(name) {
        if (const char* prev = std::getenv(name)) {
            had_ = true;
            prev_ = prev;
        }
        if (value != nullptr) {
            ::setenv(name, value, 1);
        } else {
            ::unsetenv(name);
        }
    }
    EnvGuard(const EnvGuard&) = delete;
    EnvGuard& operator=(const EnvGuard&) = delete;
    EnvGuard(EnvGuard&&) = delete;
    EnvGuard& operator=(EnvGuard&&) = delete;
    ~EnvGuard() {
        if (had_) {
            ::setenv(name_, prev_.c_str(), 1);
        } else {
            ::unsetenv(name_);
        }
    }

private:
    const char* name_;
    bool had_ = false;
    std::string prev_;
};

/// The entry points the Slang source actually defines: every `void NAME(` that follows a
/// `[shader("compute")]` attribute, ignoring comment lines.
std::set<std::string> slang_entry_points(const std::filesystem::path& src) {
    std::set<std::string> names;
    std::ifstream in(src);
    EXPECT_TRUE(in.good()) << "cannot read " << src;
    std::string line;
    bool pending = false;
    while (std::getline(in, line)) {
        const std::size_t first = line.find_first_not_of(" \t");
        if (first == std::string::npos) continue;
        if (line.compare(first, 2, "//") == 0) continue;
        if (line.find("[shader(\"compute\")]") != std::string::npos) {
            pending = true;
            continue;
        }
        if (!pending) continue;
        const std::size_t v = line.find("void ");
        if (v == std::string::npos) continue;
        const std::size_t open = line.find('(', v);
        if (open == std::string::npos) continue;
        names.insert(line.substr(v + 5, open - (v + 5)));
        pending = false;
    }
    return names;
}

/// A deterministic scatter of points over the inner magnetosphere, never at the origin: a 64-bit
/// LCG (Knuth's MMIX multiplier) mapped onto radii of 1.5-8 Re. Deterministic so a device/host
/// disagreement is reproducible, and never a lattice, so no component is systematically zero.
std::vector<float> scatter(std::size_t n) {
    std::vector<float> pos(3 * n);
    std::uint64_t s = 0x9E3779B97F4A7C15ULL;
    const auto next = [&s] {
        s = (s * 6364136223846793005ULL) + 1442695040888963407ULL;
        return static_cast<double>(s >> 11U) / static_cast<double>(1ULL << 53U);
    };
    for (std::size_t i = 0; i < n; ++i) {
        const double r = 1.5 + (6.5 * next());
        const double ct = -1.0 + (2.0 * next());
        const double ph = 2.0 * M_PI * next();
        const double st = std::sqrt(1.0 - (ct * ct));
        pos[(3 * i) + 0] = static_cast<float>(r * st * std::cos(ph));
        pos[(3 * i) + 1] = static_cast<float>(r * st * std::sin(ph));
        pos[(3 * i) + 2] = static_cast<float>(r * ct);
    }
    return pos;
}

/// g10 = 3125 nT = 5^5 makes the 3-4-5 geometry below exact; the SIGN and MAGNITUDE of the real
/// coefficient are irrelevant to the algebra, and the physical value is used where realism matters.
constexpr double kG10Exact = 3125.0;
/// IGRF-13, epoch 2020.0, degree 1 order 0 (Alken et al. 2021, Table 2).
constexpr double kG10Igrf2020 = -29404.8;

}  // namespace

// ---- the registry -----------------------------------------------------------------------------

TEST(IrbemGpu, RegistryIsNonEmptyAndSelfConsistent) {
    ASSERT_FALSE(ig::registered_kernels.empty());
    for (const ig::KernelInfo& k : ig::registered_kernels) {
        EXPECT_NE(k.name, nullptr);
        EXPECT_FALSE(k.brief.empty()) << k.name;
        EXPECT_GE(k.bindings, 2U) << k.name;   // at least one data buffer plus the dims buffer
        EXPECT_LE(k.bindings, 8U) << k.name;   // Context::kMaxBindings
        EXPECT_EQ(&ig::kernel_info(k.name), &k);
    }
}

// The negative control this whole seam exists for: a name nothing implements must be refused
// BEFORE a device is touched, both when asked directly and when launched.
TEST(IrbemGpu, UnknownKernelIsRefusedByNameAndByLaunch) {
    EXPECT_THROW((void)ig::kernel_info("irbem_no_such_kernel"), ig::UnknownKernel);
    const std::array<float, 3> pos{1.0F, 0.0F, 0.0F};
    std::array<float, 3> out{};
    const std::array<float, 1> par{1.0F};
    EXPECT_THROW(ig::dispatch_batch("irbem_no_such_kernel", pos, out, par), ig::UnknownKernel);
    try {
        (void)ig::kernel_info("irbem_no_such_kernel");
        FAIL() << "expected UnknownKernel";
    } catch (const ig::UnknownKernel& e) {
        EXPECT_NE(std::string(e.what()).find("irbem_no_such_kernel"), std::string::npos);
    }
}

// The registry and the Slang source must be the same list. This is what catches a kernel added to
// the shader with no host entry (unreachable) or a host entry with no shader (ShaderMissing at a
// customer's first launch).
TEST(IrbemGpu, RegistryMatchesTheSlangEntryPoints) {
    const std::set<std::string> in_slang = slang_entry_points(ig::shader_source_path());
    ASSERT_FALSE(in_slang.empty()) << "no [shader(\"compute\")] entry points parsed out of "
                                   << ig::shader_source_path();
    std::set<std::string> in_registry;
    for (const ig::KernelInfo& k : ig::registered_kernels) in_registry.insert(k.name);
    EXPECT_EQ(in_slang, in_registry);
}

// ---- the fp64 host reference --------------------------------------------------------------

TEST(IrbemGpu, DipoleReferenceIsExactOnIntegerGeometry) {
    // 3-4-5 in the x-z plane: r = 5, r^5 = 3125, and g10 = 3125 cancels the denominator exactly.
    const ib::FieldVector<Frame::MAG> b =
        ig::dipole_field_at(ib::Position<Frame::MAG>{fx::vec3d{3.0, 0.0, 4.0}}, kG10Exact);
    EXPECT_EQ(b.v[0], 36.0);   // 3*g10*x*z/r^5 = 3*3125*3*4/3125
    EXPECT_EQ(b.v[1], 0.0);
    EXPECT_EQ(b.v[2], 23.0);   // g10*(3z^2 - r^2)/r^5 = 3125*(48-25)/3125

    // The same triangle rotated into y-z: x and y swap roles, z is untouched.
    const ib::FieldVector<Frame::MAG> byz =
        ig::dipole_field_at(ib::Position<Frame::MAG>{fx::vec3d{0.0, 3.0, 4.0}}, kG10Exact);
    EXPECT_EQ(byz.v[0], 0.0);
    EXPECT_EQ(byz.v[1], 36.0);
    EXPECT_EQ(byz.v[2], 23.0);

    // On the dipole equator z = 0, so B is purely -z-hat: B = -g10 z-hat / r^3.
    const ib::FieldVector<Frame::MAG> eq =
        ig::dipole_field_at(ib::Position<Frame::MAG>{fx::vec3d{3.0, 4.0, 0.0}}, kG10Exact);
    EXPECT_EQ(eq.v[0], 0.0);
    EXPECT_EQ(eq.v[1], 0.0);
    EXPECT_EQ(eq.v[2], -25.0);   // g10*(0-25)/3125
}

TEST(IrbemGpu, DipoleReferenceOnTheAxisAndAtThePole) {
    // On the axis at r = 2 the field is 2*g10/r^3 * z-hat = g10/4 * z-hat.
    const ib::FieldVector<Frame::MAG> up =
        ig::dipole_field_at(ib::Position<Frame::MAG>{fx::vec3d{0.0, 0.0, 2.0}}, -32768.0);
    EXPECT_EQ(up.v[0], 0.0);
    EXPECT_EQ(up.v[1], 0.0);
    EXPECT_EQ(up.v[2], -8192.0);

    // ...and the mirror point below is the SAME vector, not its negation: a dipole field is
    // symmetric about the equator in Bz and antisymmetric in Bx/By, which is exactly the property a
    // sign slip in the z-term would break.
    const ib::FieldVector<Frame::MAG> down =
        ig::dipole_field_at(ib::Position<Frame::MAG>{fx::vec3d{0.0, 0.0, -2.0}}, -32768.0);
    EXPECT_EQ(down.v[2], up.v[2]);

    // The origin is a pole of the field; the header defines it as zero rather than a NaN.
    const ib::FieldVector<Frame::MAG> o =
        ig::dipole_field_at(ib::Position<Frame::MAG>{fx::vec3d{0.0, 0.0, 0.0}}, kG10Igrf2020);
    EXPECT_EQ(o.v[0], 0.0);
    EXPECT_EQ(o.v[1], 0.0);
    EXPECT_EQ(o.v[2], 0.0);
}

// The physical check the exact-arithmetic cases cannot make: with the REAL (negative) g10, the
// field at the northern dipole pole points DOWN — into the ground — and is twice the equatorial
// field at the same radius. Getting either wrong is the classic dipole sign error.
TEST(IrbemGpu, DipoleHasThePhysicalSignAndTheFactorOfTwo) {
    const ib::FieldVector<Frame::MAG> pole =
        ig::dipole_field_at(ib::Position<Frame::MAG>{fx::vec3d{0.0, 0.0, 4.0}}, kG10Igrf2020);
    const ib::FieldVector<Frame::MAG> eq =
        ig::dipole_field_at(ib::Position<Frame::MAG>{fx::vec3d{4.0, 0.0, 0.0}}, kG10Igrf2020);
    EXPECT_LT(pole.v[2], 0.0);              // northward on the axis => into the Earth
    EXPECT_GT(eq.v[2], 0.0);                // and out of it on the equator
    // On the axis B = 2 g10 z-hat / r^3; on the equator B = -g10 z-hat / r^3. So at equal r the
    // polar field is twice the equatorial one in MAGNITUDE and OPPOSITE in sign — the field lines
    // leave one pole and come back at the other. Both halves of that are exact here.
    EXPECT_EQ(pole.v[2], -2.0 * eq.v[2]);
    // Magnitude sanity against the textbook surface value: at r = 1 Re on the equator the centred
    // dipole gives |B| = |g10| = 29404.8 nT.
    EXPECT_DOUBLE_EQ(
        ig::dipole_field_at(ib::Position<Frame::MAG>{fx::vec3d{1.0, 0.0, 0.0}}, kG10Igrf2020)
            .magnitude(),
        29404.8);
}

// ---- the fp32 host batch lane -------------------------------------------------------------

TEST(IrbemGpu, HostBatchReproducesTheReferenceExactly) {
    const std::array<float, 12> pos{3.0F, 0.0F, 4.0F, 0.0F, 3.0F, 4.0F,
                                    3.0F, 4.0F, 0.0F, 0.0F, 0.0F, 0.0F};
    std::array<float, 12> out{};
    ig::dipole_field_host(pos, out, static_cast<float>(kG10Exact));
    const std::array<float, 12> want{36.0F, 0.0F,  23.0F, 0.0F, 36.0F, 23.0F,
                                     0.0F,  0.0F, -25.0F, 0.0F, 0.0F,  0.0F};
    EXPECT_EQ(out, want);
}

TEST(IrbemGpu, HostBatchRejectsMalformedSpans) {
    const std::array<float, 4> pos{};
    std::array<float, 4> out{};
    EXPECT_THROW(ig::dipole_field_host(pos, out, 1.0F), std::invalid_argument);
    const std::array<float, 3> pos3{};
    std::array<float, 6> out6{};
    EXPECT_THROW(ig::dipole_field_host(pos3, out6, 1.0F), std::invalid_argument);
}

// ---- shader resolution ----------------------------------------------------------------------

TEST(IrbemGpu, ShaderDirHonoursTheEnvironmentThenFallsBack) {
    {
        const EnvGuard g("CHEATAH_SPACE_IRBEM_SPV_DIR", "/some/where/else");
        EXPECT_EQ(ig::shader_dir(), std::filesystem::path("/some/where/else"));
        EXPECT_EQ(ig::shader_path("k"), std::filesystem::path("/some/where/else/k.spv"));
    }
    const EnvGuard g("CHEATAH_SPACE_IRBEM_SPV_DIR", nullptr);
    const std::filesystem::path d = ig::shader_dir();
    // The contract is "<some build tree>/shaders", NOT a directory literally named "build": which
    // build tree the shaders land in is a configuration choice (this repo configures build/gpu for
    // the GPU lane, cheatah-gpu-linalg uses build/), and asserting the name would make an
    // unrelated -B choice fail the suite. What must hold is the leaf and that it sits inside the
    // checkout rather than somewhere arbitrary on the machine.
    EXPECT_EQ(d.filename(), "shaders");
    EXPECT_NE(d.parent_path().filename(), "") << d;
    EXPECT_NE(d.string().find("cheatah-space"), std::string::npos) << d;
}

TEST(IrbemGpu, SlangSourceIsWhereTheHeaderSaysItIs) {
    const std::filesystem::path src = ig::shader_source_path();
    ASSERT_TRUE(std::filesystem::exists(src)) << src;
    EXPECT_EQ(src.filename(), "irbem.slang");
}

// ---- availability ------------------------------------------------------------------------------

TEST(IrbemGpu, AvailabilityAndReasonAgree) {
    if (ig::available()) {
        EXPECT_TRUE(ig::unavailable_reason().empty());
    } else {
        EXPECT_FALSE(ig::unavailable_reason().empty());
    }
    EXPECT_EQ(ig::compiled_with_gpu, CHEATAH_SPACE_IRBEM_HAVE_GPU != 0);

    // The operator override must be able to switch the lane off on a machine that HAS a device —
    // this is what makes the no-device failure paths testable everywhere.
    const EnvGuard g("CHEATAH_SPACE_IRBEM_NO_GPU", "1");
    EXPECT_FALSE(ig::available());
    EXPECT_FALSE(ig::unavailable_reason().empty());
}

// An override that is SET must still be able to mean "on": an exported-but-empty variable and an
// explicit "0" are both how a shell says no, and reading either as "yes, disable" would silently
// strand a fleet on the CPU lane. (Compiled without the device stack the answer is false either
// way, and this test then only pins that the strings do not crash the parse.)
TEST(IrbemGpu, TheOverrideOnlyDisablesWhenItSaysSo) {
    const bool baseline = ig::available();
    {
        const EnvGuard empty("CHEATAH_SPACE_IRBEM_NO_GPU", "");
        EXPECT_EQ(ig::available(), baseline);
    }
    {
        const EnvGuard zero("CHEATAH_SPACE_IRBEM_NO_GPU", "0");
        EXPECT_EQ(ig::available(), baseline);
    }
    {
        const EnvGuard zeroish("CHEATAH_SPACE_IRBEM_NO_GPU", "00");   // not exactly "0"
        EXPECT_FALSE(ig::available());
    }
}

// The dipole kernel is MEASURED never to beat the host lane (the sweep is recorded next to the
// registry row), so Auto must never send it to the device however large the batch gets. A
// crossover invented rather than measured is exactly how a "GPU-accelerated" path ends up 3x
// slower than the loop it replaced, which is what the earlier placeholder of 2^14 would have done.
TEST(IrbemGpu, TheDipoleKernelHasNoMeasuredCrossover) {
    EXPECT_EQ(ig::gpu_crossover("irbem_dipole_f32"), ig::never_faster_on_device);
    EXPECT_FALSE(ig::prefer_gpu("irbem_dipole_f32", 0));
    EXPECT_FALSE(ig::prefer_gpu("irbem_dipole_f32", 1U << 24U));
    EXPECT_THROW((void)ig::prefer_gpu("irbem_nope", 1), ig::UnknownKernel);
}

TEST(IrbemGpu, TheCrossoverOverrideIsHonouredAndTypoProof) {
    {
        const EnvGuard g("CHEATAH_SPACE_IRBEM_GPU_CROSSOVER", "1000");
        EXPECT_EQ(ig::gpu_crossover("irbem_dipole_f32"), 1000U);
        EXPECT_FALSE(ig::prefer_gpu("irbem_dipole_f32", 999));
        EXPECT_EQ(ig::prefer_gpu("irbem_dipole_f32", 1000), ig::available());
    }
    {
        // A typo must fall back to the measured value, not be read as "always use the device".
        const EnvGuard g("CHEATAH_SPACE_IRBEM_GPU_CROSSOVER", "1000x");
        EXPECT_EQ(ig::gpu_crossover("irbem_dipole_f32"), ig::never_faster_on_device);
    }
    {
        const EnvGuard g("CHEATAH_SPACE_IRBEM_GPU_CROSSOVER", "not-a-number");
        EXPECT_EQ(ig::gpu_crossover("irbem_dipole_f32"), ig::never_faster_on_device);
    }
    const EnvGuard c("CHEATAH_SPACE_IRBEM_GPU_CROSSOVER", "0");
    const EnvGuard g("CHEATAH_SPACE_IRBEM_NO_GPU", "1");
    EXPECT_FALSE(ig::prefer_gpu("irbem_dipole_f32", 1U << 20U));   // no device beats any crossover
}

// ---- the launcher's named failures ---------------------------------------------------------

TEST(IrbemGpu, MalformedBatchesAreRefused) {
    std::array<float, 4> four{};
    std::array<float, 3> three{};
    std::array<float, 6> six{};
    const std::array<float, 1> par{1.0F};
    const std::array<float, 2> par2{1.0F, 2.0F};
    EXPECT_THROW(ig::dispatch_batch("irbem_dipole_f32", four, four, par), std::invalid_argument);
    EXPECT_THROW(ig::dispatch_batch("irbem_dipole_f32", three, six, par), std::invalid_argument);
    EXPECT_THROW(ig::dispatch_batch("irbem_dipole_f32", three, three, par2), std::invalid_argument);
}

TEST(IrbemGpu, EmptyBatchIsANoOpEvenWithNoDeviceAndNoShader) {
    const EnvGuard sd("CHEATAH_SPACE_IRBEM_SPV_DIR", "/definitely/not/here");
    const EnvGuard ng("CHEATAH_SPACE_IRBEM_NO_GPU", "1");
    const std::span<const float> none;
    std::span<float> nout;
    const std::array<float, 1> par{1.0F};
    EXPECT_NO_THROW(ig::dispatch_batch("irbem_dipole_f32", none, nout, par));
}

// A missing .spv is a BUILD failure and is diagnosed before availability, so it surfaces on a
// machine with no GPU too — which is exactly the machine where the shader step gets skipped.
TEST(IrbemGpu, MissingShaderIsNamedAndIsDiagnosedBeforeAvailability) {
    const std::filesystem::path empty =
        std::filesystem::temp_directory_path() / "irbem_gpu_test_no_shaders";
    std::filesystem::create_directories(empty);
    const EnvGuard sd("CHEATAH_SPACE_IRBEM_SPV_DIR", empty.string().c_str());
    const EnvGuard ng("CHEATAH_SPACE_IRBEM_NO_GPU", "1");   // no device either: shader still wins
    const std::array<float, 3> pos{1.0F, 0.0F, 0.0F};
    std::array<float, 3> out{};
    const std::array<float, 1> par{1.0F};
    try {
        ig::dispatch_batch("irbem_dipole_f32", pos, out, par);
        FAIL() << "expected ShaderMissing";
    } catch (const ig::ShaderMissing& e) {
        EXPECT_NE(std::string(e.what()).find(empty.string()), std::string::npos);
        EXPECT_NE(std::string(e.what()).find("irbem_dipole_f32"), std::string::npos);
    }
}

// With the shader present but no device, the launcher must say so by name and carry the reason.
TEST(IrbemGpu, NoDeviceIsNamedAndCarriesTheReason) {
    if (!std::filesystem::exists(ig::shader_path("irbem_dipole_f32")))
        GTEST_SKIP() << "shader not built at " << ig::shader_path("irbem_dipole_f32");
    const EnvGuard ng("CHEATAH_SPACE_IRBEM_NO_GPU", "1");
    const std::array<float, 3> pos{1.0F, 0.0F, 0.0F};
    std::array<float, 3> out{};
    try {
        ig::dipole_field_gpu(pos, out, 1.0F);
        FAIL() << "expected GpuUnavailable";
    } catch (const ig::GpuUnavailable& e) {
        const std::string msg = e.what();
        EXPECT_NE(msg.find("irbem_dipole_f32"), std::string::npos) << msg;
        // The reason must be the ACTUAL one. Compiled against the device stack, that is the
        // operator override above; compiled without it, it is the absent header — and saying the
        // wrong one would send a user hunting the wrong problem.
        EXPECT_NE(msg.find(ig::compiled_with_gpu ? "CHEATAH_SPACE_IRBEM_NO_GPU"
                                                 : "cheatah-gpu-linalg"),
                  std::string::npos)
            << msg;
    }
    // ...and a FORCED device lane must throw rather than silently fall back to the CPU, or a test
    // that meant to exercise the device would pass without ever reaching one.
    EXPECT_THROW(ig::dipole_field(pos, out, 1.0F, ig::Lane::Gpu), ig::GpuUnavailable);
}

// ---- lane selection --------------------------------------------------------------------------

TEST(IrbemGpu, AutoAndForcedHostLanesComputeTheSameThing) {
    const std::array<float, 6> pos{3.0F, 0.0F, 4.0F, 0.0F, 3.0F, 4.0F};
    std::array<float, 6> a{};
    std::array<float, 6> b{};
    ig::dipole_field(pos, a, static_cast<float>(kG10Exact));                  // Auto, small batch
    ig::dipole_field(pos, b, static_cast<float>(kG10Exact), ig::Lane::Host);  // forced
    EXPECT_EQ(a, b);
    const std::array<float, 6> want{36.0F, 0.0F, 23.0F, 0.0F, 36.0F, 23.0F};
    EXPECT_EQ(a, want);
}

// ---- the device lane ---------------------------------------------------------------------------

TEST(IrbemGpu, KernelIsExactOnIntegerGeometry) {
    if (!ig::available()) GTEST_SKIP() << "no device: " << ig::unavailable_reason();
    const std::array<float, 12> pos{3.0F, 0.0F, 4.0F, 0.0F, 3.0F, 4.0F,
                                    3.0F, 4.0F, 0.0F, 0.0F, 0.0F, 0.0F};
    std::array<float, 12> out{};
    ig::dipole_field(pos, out, static_cast<float>(kG10Exact), ig::Lane::Gpu);
    // Every intermediate is a small integer, so contraction and reassociation cannot change the
    // answer: this is exact even on a driver free to fuse.
    const std::array<float, 12> want{36.0F, 0.0F,  23.0F, 0.0F, 36.0F, 23.0F,
                                     0.0F,  0.0F, -25.0F, 0.0F, 0.0F,  0.0F};
    EXPECT_EQ(out, want);
}

TEST(IrbemGpu, KernelMatchesTheHostLaneToTheBlocalBudget) {
    if (!ig::available()) GTEST_SKIP() << "no device: " << ig::unavailable_reason();
    const std::size_t n = 1U << 20U;   // above gpu_crossover_points, so Lane::Auto picks the device
    const std::vector<float> pos = scatter(n);
    std::vector<float> host(3 * n);
    std::vector<float> dev(3 * n);
    ig::dipole_field_host(pos, host, static_cast<float>(kG10Igrf2020));
    // Auto, with the crossover overridden down so the automatic policy — not a forced lane — is
    // what puts this batch on the device.
    const EnvGuard over("CHEATAH_SPACE_IRBEM_GPU_CROSSOVER", "1024");
    ASSERT_TRUE(ig::prefer_gpu("irbem_dipole_f32", n));
    ig::dipole_field(pos, dev, static_cast<float>(kG10Igrf2020));

    double worst = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        const double hx = host[(3 * i) + 0];
        const double hy = host[(3 * i) + 1];
        const double hz = host[(3 * i) + 2];
        const double mag = std::sqrt((hx * hx) + (hy * hy) + (hz * hz));
        ASSERT_GT(mag, 0.0) << "point " << i << " landed on the origin";
        for (std::size_t c = 0; c < 3; ++c)
            worst = std::max(worst, std::abs(dev[(3 * i) + c] - host[(3 * i) + c]) / mag);
    }
    // docs/ERROR_BUDGET.md §4: Bgeo/Blocal are budgeted at 1e-6 relative. §6: the device lane
    // cannot carry bit-goldens because FMA contraction is the driver's choice, so it is held to
    // that budget against the host reference instead.
    EXPECT_LT(worst, 1e-6) << "max relative device-vs-host deviation";
    std::printf("[ measured ] max relative device-vs-host deviation over %zu points: %.3e\n", n,
                worst);
}

TEST(IrbemGpu, KernelIsBitIdenticalRunToRun) {
    if (!ig::available()) GTEST_SKIP() << "no device: " << ig::unavailable_reason();
    const std::size_t n = 4096;
    const std::vector<float> pos = scatter(n);
    std::vector<float> a(3 * n);
    std::vector<float> b(3 * n);
    ig::dipole_field_gpu(pos, a, static_cast<float>(kG10Igrf2020));
    ig::dipole_field_gpu(pos, b, static_cast<float>(kG10Igrf2020));
    EXPECT_EQ(std::memcmp(a.data(), b.data(), a.size() * sizeof(float)), 0);
}

// The workgroup size is 256 and the grid is rounded up to it, so the tail of a batch that is not a
// multiple of 256 is where an unguarded kernel writes out of bounds or leaves stale data.
TEST(IrbemGpu, RaggedBatchesAreComputedInFull) {
    if (!ig::available()) GTEST_SKIP() << "no device: " << ig::unavailable_reason();
    for (const std::size_t n : {std::size_t{1}, std::size_t{255}, std::size_t{257},
                                std::size_t{1000}}) {
        const std::vector<float> pos = scatter(n);
        std::vector<float> host(3 * n);
        std::vector<float> dev(3 * n, -12345.0F);   // poison, so an unwritten element is visible
        ig::dipole_field_host(pos, host, static_cast<float>(kG10Igrf2020));
        ig::dipole_field_gpu(pos, dev, static_cast<float>(kG10Igrf2020));
        for (std::size_t i = 0; i < 3 * n; ++i) {
            ASSERT_NE(dev[i], -12345.0F) << "n=" << n << " element " << i << " never written";
            ASSERT_NEAR(dev[i], host[i], 1e-3 * std::abs(host[i]) + 1e-6)
                << "n=" << n << " element " << i;
        }
    }
}

// Both destructor paths of the SPIR-V directory scope: the variable was already set, and it was
// not. Getting the "was already set" path wrong would leave cheatah-gpu-linalg's own kernels
// pointed at our directory for the rest of the process.
TEST(IrbemGpu, TheSpvDirScopeIsRestoredEitherWay) {
    if (!ig::available()) GTEST_SKIP() << "no device: " << ig::unavailable_reason();
    const std::array<float, 3> pos{3.0F, 0.0F, 4.0F};
    std::array<float, 3> out{};
    {
        const EnvGuard g("CHEATAH_GPU_LINALG_SPV_DIR", nullptr);
        ig::dipole_field_gpu(pos, out, static_cast<float>(kG10Exact));
        EXPECT_EQ(std::getenv("CHEATAH_GPU_LINALG_SPV_DIR"), nullptr);
    }
    {
        const std::string mine = ig::shader_dir().string();
        const EnvGuard g("CHEATAH_GPU_LINALG_SPV_DIR", mine.c_str());
        ig::dipole_field_gpu(pos, out, static_cast<float>(kG10Exact));
        ASSERT_NE(std::getenv("CHEATAH_GPU_LINALG_SPV_DIR"), nullptr);
        EXPECT_EQ(std::string(std::getenv("CHEATAH_GPU_LINALG_SPV_DIR")), mine);
    }
    EXPECT_EQ(out[0], 36.0F);
}

TEST(IrbemGpu, ThroughputIsReported) {
    if (!ig::available()) GTEST_SKIP() << "no device: " << ig::unavailable_reason();
    const std::size_t n = 1U << 20U;
    const std::vector<float> pos = scatter(n);
    std::vector<float> out(3 * n);
    using clock = std::chrono::steady_clock;
    ig::dipole_field_gpu(pos, out, static_cast<float>(kG10Igrf2020));   // warm the pipeline cache
    const clock::time_point t0 = clock::now();
    ig::dipole_field_gpu(pos, out, static_cast<float>(kG10Igrf2020));
    const clock::time_point t1 = clock::now();
    ig::dipole_field_host(pos, out, static_cast<float>(kG10Igrf2020));
    const clock::time_point t2 = clock::now();
    const double gpu_s = std::chrono::duration<double>(t1 - t0).count();
    const double cpu_s = std::chrono::duration<double>(t2 - t1).count();
    std::printf("[ measured ] %zu points: device %.3f ms (%.1f Mpts/s, transfers included), "
                "host %.3f ms (%.1f Mpts/s)\n",
                n, gpu_s * 1e3, static_cast<double>(n) / gpu_s / 1e6, cpu_s * 1e3,
                static_cast<double>(n) / cpu_s / 1e6);
    EXPECT_GT(gpu_s, 0.0);
    EXPECT_GT(cpu_s, 0.0);
}
