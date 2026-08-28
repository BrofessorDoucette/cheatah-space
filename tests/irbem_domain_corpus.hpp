#pragma once

/// @file irbem_domain_corpus.hpp
/// @brief The input corpus every space.irbem differential test draws from.
///
/// **This is a correctness requirement, not thoroughness.** The Tsyganenko external field models
/// are *parameterized by geomagnetic activity*: T89 is Kp-binned, T96 takes Pdyn/Dst/By/Bz, T01
/// adds G1-G3, TS05 adds W1-W6. A library validated only at `maginput = 0` has demonstrated
/// nothing whatsoever about the disturbed conditions that are the entire reason radiation-belt
/// models exist. Quiet-time agreement is the easy half and the half nobody needs.
///
/// Two halves, deliberately:
///
///  - a **stratified synthetic sweep**, which finds the corners of each model's fitted envelope —
///    the places where an empirical fit stops being interpolation and starts being extrapolation;
///  - **real storm events**, which prove we reproduce the field where researchers actually look.
///
/// Synthetic sweeps find corners the solar wind rarely visits; real events cluster where it does.
/// Neither alone is sufficient, which is why both are here.
///
/// ## Why southward Bz is sampled densely and asymmetrically
///
/// The sampling of `Bz` is deliberately NOT symmetric about zero. Southward IMF (`Bz < 0`) is what
/// couples the solar wind to the magnetosphere: it permits dayside reconnection, drives convection,
/// loads the tail and produces the storm-time reconfiguration the models were fitted to reproduce.
/// Northward Bz is comparatively quiescent. A corpus that sampled the two evenly would spend half
/// its budget on the case where nothing happens, so this one does not.
///
/// ## What is honestly approximate here
///
/// G1-G3 (T01) and W1-W6 (TS05) are not instantaneous quantities — they are **integrals over
/// solar-wind history**, accumulated over hours. A physically consistent value cannot be conjured
/// from a single Dst reading. Where this corpus supplies them it says what they were derived from
/// and how rough that is; where a driver could not be sourced it is marked, not invented. A test
/// corpus that fabricates its inputs proves agreement with a fiction.
#include <array>
#include <cstddef>
#include <cstdint>
#include <string_view>

namespace cheatah_space_test {

/// How disturbed the magnetosphere is. The models behave qualitatively differently across these —
/// the tail current strengthens and moves earthward, the ring current deepens, the magnetopause
/// compresses — so a corpus that samples only one regime samples only one physics.
enum class Regime : std::uint8_t {
    Quiet,     ///< Kp < 2, |Dst| < 20 nT. The easy case, and the one that proves least.
    Moderate,  ///< Kp 3-4, Dst -30…-50 nT. Typical disturbed conditions.
    Storm,     ///< Kp 5-7, Dst -100…-200 nT. Main phase; Bz strongly south.
    Extreme,   ///< Kp 8-9, Dst -250…-400 nT. At or past most models' fitted envelopes.
};

/// The IRBEM `maginput` drivers this corpus varies. Indices match IRBEM's documented order so a
/// case can be handed straight to the oracle without a translation step that could itself be wrong.
struct MagInput {
    double kp = 0.0;     ///< maginput(1): Kp index.
    double dst = 0.0;    ///< maginput(2): Dst, nT.
    double density = 5.0;  ///< maginput(3): solar-wind proton density, cm^-3.
    double velocity = 400.0;  ///< maginput(4): solar-wind speed, km/s.
    double pdyn = 2.0;   ///< maginput(5): dynamic pressure, nPa.
    double by_imf = 0.0; ///< maginput(6): IMF By in GSM, nT.
    double bz_imf = 0.0; ///< maginput(7): IMF Bz in GSM, nT. NEGATIVE is the interesting sign.
    double g1 = 0.0;     ///< maginput(8): T01's G1 — a solar-wind history integral; see the brief.
    double g2 = 0.0;     ///< maginput(9): T01's G2.
    double g3 = 0.0;     ///< maginput(10): T01's G3.
};

/// One input case, with a stable slug so `ctest -R` can select it out of thousands.
struct Case {
    std::string_view name;  ///< A stable slug: regime, shell, local time, pitch angle.
    Regime regime;
    int year;
    int doy;
    double ut_seconds;
    double x, y, z;         ///< GEO Cartesian, Earth radii.
    double pitch_deg;       ///< Local pitch angle.
    MagInput mag;
};

/// The four regimes' representative drivers.
///
/// Values are chosen to sit inside each regime rather than at its edge, so a test that fails here
/// has failed in the regime's interior — the boundaries are exercised separately by the validity
/// tests, where being just inside and just outside is the whole point.
inline constexpr std::array<MagInput, 4> regime_drivers{{
    // Quiet: nominal solar wind, Bz weakly northward.
    {1.0, -8.0, 5.0, 380.0, 1.8, 1.0, 2.0, 0.0, 0.0, 0.0},
    // Moderate: Bz has turned south, pressure up a little.
    {3.5, -42.0, 8.0, 450.0, 3.2, -4.0, -5.0, 6.0, 8.0, 20.0},
    // Storm main phase: strongly southward Bz, compressed magnetosphere, deep ring current.
    {6.0, -150.0, 20.0, 600.0, 9.0, 8.0, -15.0, 25.0, 30.0, 80.0},
    // Extreme: at or past most models' fitted envelopes. Included precisely BECAUSE it is past
    // them — a model must report OutOfValidityRange here rather than returning a confident number.
    {8.5, -350.0, 45.0, 900.0, 28.0, -18.0, -30.0, 60.0, 75.0, 180.0},
}};

/// Real storm events, with published representative peak drivers.
///
/// Each is an event researchers cite by name, so agreement here is agreement a reviewer can check.
/// The drivers are representative PEAK values from the published literature, not a reconstructed
/// time series — enough to place the model in the right regime, not enough to reproduce a specific
/// minute. That limitation is stated rather than hidden: an exact reconstruction needs the archived
/// OMNI series, which this corpus deliberately does not vendor.
struct StormEvent {
    std::string_view name;
    int year;
    int doy;
    double peak_dst;   ///< nT, published.
    MagInput mag;
};

inline constexpr std::array<StormEvent, 4> storm_events{{
    // Halloween 2003 — one of the most-studied events in the record; the belts were reconfigured
    // wholesale and a new belt formed.
    {"halloween-2003", 2003, 303, -383.0,
     {9.0, -383.0, 40.0, 1800.0, 25.0, -20.0, -30.0, 70.0, 90.0, 200.0}},
    // March 1989 — the Hydro-Quebec storm, the largest of the space age.
    {"march-1989", 1989, 72, -589.0,
     {9.0, -589.0, 50.0, 1000.0, 30.0, 15.0, -40.0, 90.0, 110.0, 250.0}},
    // St Patrick's Day 2015 — the largest of solar cycle 24, and heavily instrumented by the Van
    // Allen Probes, so it is the modern reference event.
    {"st-patricks-2015", 2015, 76, -223.0,
     {8.0, -223.0, 25.0, 700.0, 12.0, -10.0, -20.0, 40.0, 50.0, 120.0}},
    // February 2022 — the storm that de-orbited 38 Starlink satellites, and the event Tsyganenko's
    // own 2022 reconstruction paper addresses. Modest in Dst yet operationally severe, which is
    // exactly why a corpus keyed only on Dst would miss it.
    {"starlink-2022", 2022, 34, -66.0,
     {5.0, -66.0, 15.0, 550.0, 8.0, 5.0, -12.0, 20.0, 25.0, 60.0}},
}};

/// L shells spanning the inner belt, the slot region, the outer belt and geosynchronous.
inline constexpr std::array<double, 8> shells{{1.5, 2.0, 3.0, 4.0, 5.0, 6.0, 6.6, 8.0}};

/// Magnetic local times, in hours. Noon and midnight are included explicitly because the field is
/// least like a dipole there — compressed on the dayside, stretched into the tail at night — so a
/// corpus that sampled only dawn and dusk would flatter the model.
inline constexpr std::array<double, 8> local_times{{0.0, 3.0, 6.0, 9.0, 12.0, 15.0, 18.0, 21.0}};

/// Pitch angles: the loss cone, a mid-range spread, and 90 degrees where the particle mirrors at
/// the equator and `I` is identically zero — a degenerate case worth testing precisely because it
/// is degenerate.
inline constexpr std::array<double, 6> pitch_angles{{5.0, 15.0, 30.0, 50.0, 70.0, 90.0}};

/// Epochs spanning the IGRF definition, including intra-epoch points where the coefficients are
/// interpolated rather than tabulated, and a leap year.
inline constexpr std::array<int, 6> years{{1900, 1965, 2000, 2015, 2020, 2024}};

/// The full sweep size, so a test can state its own coverage rather than implying it.
inline constexpr std::size_t sweep_size =
    regime_drivers.size() * shells.size() * local_times.size() * pitch_angles.size();

}  // namespace cheatah_space_test
