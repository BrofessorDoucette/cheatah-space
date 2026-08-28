// Unit tests for space.irbem's frame layer — Frame/FrameKind, the sysaxes boundary, and the
// frame-tagged Position<F> / FieldVector<F>.
//
// The QA gate requires 100% lines AND functions over space/**.hpp, and clang source-based coverage
// counts each template INSTANTIATION separately. So the per-frame exercises are themselves
// templates, invoked once per frame through a fold over the whole Frame enum: adding a frame adds
// no test code, and no instantiation can be left unexercised by accident.
//
// Arithmetic uses exactly-representable values (halves, and the 3-4-5 triangle) so every assertion
// is an exact `==`. A tolerance here would hide precisely the component-swap and frame-mix bugs
// this layer exists to prevent.
#include <gtest/gtest.h>

#include <string_view>
#include <type_traits>
#include <utility>

#include "space/irbem/frames.hpp"

namespace ib = cheatah::space::irbem;
namespace fx = cheatah::fixarray;

using ib::Frame;
using ib::FrameKind;

// ---- compile-time surface ---------------------------------------------------------------------

// The frame tag is the whole point: two frames are two types, so a mismatch cannot compile.
static_assert(!std::is_same_v<ib::Position<Frame::GEO>, ib::Position<Frame::GSM>>);
static_assert(!std::is_same_v<ib::FieldVector<Frame::GSM>, ib::Position<Frame::GSM>>);

// ...and it must cost nothing. A Position is exactly its vector, and trivially copyable, so it
// passes in registers and can be memcpy'd into a GPU upload buffer.
static_assert(sizeof(ib::Position<Frame::GEO>) == sizeof(fx::vec3d));
static_assert(sizeof(ib::FieldVector<Frame::GSM>) == sizeof(fx::vec3d));
static_assert(std::is_trivially_copyable_v<ib::Position<Frame::GEO>>);
static_assert(std::is_trivially_copyable_v<ib::FieldVector<Frame::GSM>>);

// The kind concepts partition the frames.
static_assert(ib::CartesianFrame<Frame::GEO> && ib::CartesianFrame<Frame::HEEQ>);
static_assert(!ib::CartesianFrame<Frame::GDZ> && !ib::CartesianFrame<Frame::RLL>);
static_assert(ib::AngularFrame<Frame::GDZ> && ib::AngularFrame<Frame::SPH>);
static_assert(!ib::AngularFrame<Frame::SM>);

// The boundary conversion is usable at compile time, which is what lets a compile-time model
// selection skip the runtime switch entirely.
static_assert(ib::frame_from_sysaxes(2) == Frame::GSM);
static_assert(!ib::frame_from_sysaxes(9).has_value());
static_assert(ib::sysaxes_of(Frame::RLL) == 8);
static_assert(!ib::sysaxes_of(Frame::HEE).has_value());

// ---- the frame table ---------------------------------------------------------------------------

/// Every enumerator, so a loop over frames cannot silently miss one.
constexpr std::array<Frame, ib::frame_count> kAllFrames{
    Frame::GDZ, Frame::GEO, Frame::GSM, Frame::GSE, Frame::SM,  Frame::GEI,
    Frame::MAG, Frame::SPH, Frame::RLL, Frame::HEE, Frame::HAE, Frame::HEEQ};

TEST(IrbemFrames, EveryFrameHasAKindAndAName) {
    for (const Frame f : kAllFrames) {
        const std::string_view name = ib::name_of(f);
        EXPECT_FALSE(name.empty()) << static_cast<int>(f);
        EXPECT_NE("?", name) << "frame " << static_cast<int>(f) << " fell through name_of";

        // kind_of must be total; the three non-Cartesian frames are exactly GDZ, SPH and RLL.
        const FrameKind kind = ib::kind_of(f);
        const bool angular = (f == Frame::GDZ || f == Frame::SPH || f == Frame::RLL);
        EXPECT_EQ(angular, kind != FrameKind::Cartesian) << name;
    }
    EXPECT_EQ(FrameKind::Geodetic, ib::kind_of(Frame::GDZ));
    EXPECT_EQ(FrameKind::Spherical, ib::kind_of(Frame::SPH));
    EXPECT_EQ("GSM", ib::name_of(Frame::GSM));
}

// name_of stays total even for a value outside the enumerator set. Frame's underlying type is
// uint8_t, so every uint8_t IS a valid Frame value — this is well-defined, not a UB probe, and it
// covers the fall-through that keeps the function from running off the end.
TEST(IrbemFrames, NameOfIsTotalForAnUnnamedValue) {
    // Frame fixes its underlying type to uint8_t, so by [dcl.enum]/8 every uint8_t value IS a
    // valid Frame value and this cast is well-defined — the analyzer's check is about the
    // ENUMERATOR list, which is deliberately narrower than the value range here.
    // NOLINTNEXTLINE(clang-analyzer-optin.core.EnumCastOutOfRange)
    EXPECT_EQ("?", ib::name_of(static_cast<Frame>(99)));
}

TEST(IrbemFrames, SysaxesRoundTripsForTheCodedFramesAndRejectsTheRest) {
    for (const Frame f : kAllFrames) {
        const auto code = ib::sysaxes_of(f);
        if (code.has_value()) {
            EXPECT_EQ(f, ib::frame_from_sysaxes(*code)) << ib::name_of(f);
        } else {
            // The heliospheric frames have no sysaxes spelling — that must be reported, not faked.
            EXPECT_TRUE(f == Frame::HEE || f == Frame::HAE || f == Frame::HEEQ) << ib::name_of(f);
        }
    }
    // Out-of-range codes are refused rather than defaulted to a frame: a bad sysaxes must surface
    // as a named error at the API boundary, never as silently-GDZ coordinates.
    EXPECT_FALSE(ib::frame_from_sysaxes(-1).has_value());
    EXPECT_FALSE(ib::frame_from_sysaxes(9).has_value());
    EXPECT_FALSE(ib::frame_from_sysaxes(1000).has_value());
}

// ---- Position / FieldVector, once per frame -----------------------------------------------------

/// Exercises every member a Cartesian-frame Position has.
template <Frame F>
void exercise_cartesian_position() {
    const ib::Position<F> p{fx::vec3d{3.0, 4.0, 12.0}};
    EXPECT_EQ(F, ib::Position<F>::frame);
    EXPECT_EQ(13.0, p.radius()) << ib::name_of(F);  // 3-4-12-13, exact in binary
    EXPECT_EQ(p, (ib::Position<F>{fx::vec3d{3.0, 4.0, 12.0}}));
    EXPECT_FALSE(p == (ib::Position<F>{fx::vec3d{3.0, 4.0, 12.5}}));
}

/// Exercises every member an angular-frame Position has.
template <Frame F>
void exercise_angular_position() {
    const ib::Position<F> p{fx::vec3d{6.5, -45.0, 120.0}};
    EXPECT_EQ(F, ib::Position<F>::frame);
    EXPECT_EQ(6.5, p.radius()) << ib::name_of(F);
    EXPECT_EQ(-45.0, p.latitude()) << ib::name_of(F);
    EXPECT_EQ(120.0, p.longitude()) << ib::name_of(F);
    EXPECT_EQ(p, (ib::Position<F>{fx::vec3d{6.5, -45.0, 120.0}}));
    EXPECT_FALSE(p == (ib::Position<F>{fx::vec3d{6.5, -45.0, 121.0}}));
}

/// Exercises every member a FieldVector has, in frame @p F.
template <Frame F>
void exercise_field_vector() {
    const ib::FieldVector<F> internal{fx::vec3d{0.0, 3.0, 4.0}};
    const ib::FieldVector<F> external{fx::vec3d{2.0, 1.0, 0.5}};

    EXPECT_EQ(F, ib::FieldVector<F>::frame);
    EXPECT_EQ(5.0, internal.magnitude()) << ib::name_of(F);

    // Superposition — the operation the whole type exists to protect.
    const ib::FieldVector<F> total = internal + external;
    EXPECT_EQ((ib::FieldVector<F>{fx::vec3d{2.0, 4.0, 4.5}}), total) << ib::name_of(F);

    // ...and its inverse, which is how finite-difference gradients are formed.
    EXPECT_EQ(internal, total - external) << ib::name_of(F);

    EXPECT_EQ((ib::FieldVector<F>{fx::vec3d{0.0, 1.5, 2.0}}), internal * 0.5) << ib::name_of(F);
    EXPECT_FALSE(internal == external);
}

/// Drive the right exercise for each frame, chosen by its kind at compile time.
template <Frame F>
void exercise_frame() {
    if constexpr (ib::CartesianFrame<F>) {
        exercise_cartesian_position<F>();
    } else {
        exercise_angular_position<F>();
    }
    exercise_field_vector<F>();
}

/// Fold over the whole enum, so every instantiation the library can produce is exercised.
template <std::size_t... I>
void exercise_all_frames(std::index_sequence<I...> /*frames*/) {
    (exercise_frame<static_cast<Frame>(I)>(), ...);
}

TEST(IrbemFrames, PositionAndFieldVectorWorkInEveryFrame) {
    exercise_all_frames(std::make_index_sequence<ib::frame_count>{});
}

// A default-constructed value is the zero vector — the state a caller gets before a routine fills
// it in, and the one an uninitialized-memory bug would perturb.
TEST(IrbemFrames, DefaultsAreZero) {
    const ib::Position<Frame::GEO> p{};
    const ib::FieldVector<Frame::GSM> b{};
    EXPECT_EQ(0.0, p.radius());
    EXPECT_EQ(0.0, b.magnitude());
}
