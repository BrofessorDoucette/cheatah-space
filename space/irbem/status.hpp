#pragma once

/**
 * @file status.hpp
 * @brief space.irbem — the failure vocabulary: a value and an inspectable reason, in place of a
 *        magic number that means "something went wrong, guess what".
 *
 * IRBEM answers every question with a `double`. When the question cannot be answered — the field
 * line never closed, the mirror-point search ran out of iterations, the drivers describe a
 * magnetosphere the empirical fit has never seen — it returns `baddata = -1e31` and the caller is
 * left to infer why. Two of those three cases are *physics* the caller would act on differently:
 * an open field line at geosynchronous during a storm is a result, a non-converged root-find is a
 * bug report, and a driver vector outside T96's fitting envelope is neither — the number is
 * perfectly computable and merely unsupported by the data the model was fitted to.
 *
 * So a routine here returns @ref Result: the value it computed AND why, if anything, it is not to
 * be trusted. Three properties are load-bearing:
 *
 * - **The value is always returned.** @ref Status::OutOfValidityRange in particular does not
 *   suppress the answer; extrapolating an empirical fit is a decision only the caller can make, and
 *   a library that silently substitutes a sentinel has made it for them. What this header removes
 *   is the *silence*, not the number.
 * - **@ref Result is trivially copyable and no larger than its value plus one alignment unit**
 *   (`static_assert`ed below). The GPU lane writes one `uint` status per point into a buffer beside
 *   the values, so a single point that leaves the model's envelope is reported for that point
 *   instead of branching the whole workgroup out — see @ref status_code / @ref status_from_code for
 *   the wire encoding, which is the enumerator's own value and is fixed by that ABI.
 * - **`baddata` survives only at the outermost C-compatible boundary**, through
 *   @ref to_baddata / @ref is_baddata. Inside the library nothing compares a `double` against
 *   -1e31, because that comparison is exactly the defect this file exists to delete.
 *
 * ### Validity envelopes
 *
 * Half of this header is a table. Every empirical magnetospheric field model is a fit to a data
 * set, and every data set has edges: T96 saw no hour with `Pdyn > 10 nPa`, so what it returns there
 * is an extrapolation of a functional form, not a measurement-backed field. IRBEM knows those edges
 * — its `docs/source/api/general_information.rst` `kext` table states them — and does nothing with
 * them. @ref envelope_of turns the same table into `constexpr` data and @ref check_validity reports
 * against it.
 *
 * The bounds below are quoted from that table, which is IRBEM's own restatement of each model's
 * published range; the per-model citation in each @ref ValidityEnvelope names the paper the range
 * comes from. **Where a bound is genuinely unpublished the entry says so rather than guessing**: an
 * entry with infinite `lo`/`hi` means "this model reads this driver and no upper or lower limit has
 * been published for it" — which is the literal wording of the IRBEM table's rows for the T01-storm
 * and TS05 storm models, whose W- and G-parameters are unbounded by construction. Absence of a
 * bound is a documented fact here, never an invented number.
 *
 * @note Nothing in this header allocates, throws, or has a non-constant-time path. Every function
 *       is `constexpr` and total: a value outside an enumerator list returns a defined answer
 *       ("?", or @ref Status::DomainError) rather than falling off the end.
 */

#include <array>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <initializer_list>
#include <limits>
#include <optional>
#include <string_view>
#include <type_traits>

#include "context.hpp"

namespace cheatah::space::irbem {

// -------------------------------------------------------------------------------------------
// The status enumeration
// -------------------------------------------------------------------------------------------

/**
 * Why a computed value may not be trustworthy — the whole failure vocabulary of this library.
 *
 * The enumerators are ordered so that `Ok == 0` and the underlying values are the wire encoding the
 * GPU lane writes (@ref status_code). Appending is safe; renumbering is an ABI break.
 *
 * The distinction that matters is between the *physical* conditions and the *input* one.
 * @ref OpenFieldLine and @ref NotConverged are outcomes of a computation on well-formed input;
 * @ref OutOfValidityRange says the input is well-formed and outside the fit; @ref DomainError says
 * the input was never usable. A caller triaging a batch of points wants those three piles separate.
 */
enum class Status : std::uint8_t {
    Ok = 0,             ///< The value is the model's answer, with no caveat.
    OutOfValidityRange, ///< Outside the model's PUBLISHED fitting envelope; the value is an
                        ///< extrapolation of the functional form, and is still returned.
    OpenFieldLine,      ///< The trace left the domain without closing — physics, not a defect.
    NotConverged,       ///< A root-find or drift-shell iteration hit its cap.
    ParametersMissing,  ///< TS07D / RBF coefficient files were not provisioned; see
                        ///< @ref check_parameters.
    DomainError,        ///< NaN or infinite input, a radius inside the Earth, an epoch outside the
                        ///< model's definition, or an unrecognised model key.
};

/// How many enumerators @ref Status has — the bound the wire decoder checks against.
inline constexpr std::uint32_t status_count = 6;

/**
 * A human-readable reason, for a log line or a test failure message.
 *
 * @param s the status.
 * @return a static string such as `"outside the model's published validity range"`; `"?"` for a
 *         value outside the enumerator list, which keeps this function total.
 * @complexity O(1).
 * @alloc none.
 * @test IrbemStatus.DescribeCoversEveryEnumerator
 */
[[nodiscard]] constexpr std::string_view describe(Status s) {
    switch (s) {
        case Status::Ok:
            return "ok";
        case Status::OutOfValidityRange:
            return "outside the model's published validity range";
        case Status::OpenFieldLine:
            return "field line left the domain without closing";
        case Status::NotConverged:
            return "iteration limit reached without convergence";
        case Status::ParametersMissing:
            return "model coefficient files not provisioned";
        case Status::DomainError:
            return "input outside the domain of definition";
    }
    return "?";  // a value outside the enumerator list; keeps this total
}

/**
 * Whether @p s is the no-caveat status.
 *
 * Spelled as a function rather than left to `== Status::Ok` at every call site so that "success"
 * has one definition — if a future status ever has to be treated as benign, it is added here.
 *
 * @param s the status.
 * @return `true` only for @ref Status::Ok.
 * @complexity O(1).
 * @alloc none.
 * @test IrbemStatus.IsOkIsExactlyOk
 */
[[nodiscard]] constexpr bool is_ok(Status s) { return s == Status::Ok; }

/**
 * Compose two checks: the first non-`Ok` of @p a and @p b.
 *
 * A caller that must satisfy several envelopes at once (drivers, position, coefficient files) wants
 * one status out, and wants the *first* reason rather than the last, because the checks are ordered
 * cheapest-and-most-fundamental first.
 *
 * @param a the first check's result.
 * @param b the second check's result.
 * @return @p a when it is not @ref Status::Ok, otherwise @p b.
 * @complexity O(1).
 * @alloc none.
 * @test IrbemStatus.FirstFailureKeepsTheFirstNonOk
 */
[[nodiscard]] constexpr Status first_failure(Status a, Status b) { return is_ok(a) ? b : a; }

// -------------------------------------------------------------------------------------------
// The wire encoding — one uint per point, beside the value
// -------------------------------------------------------------------------------------------

/**
 * The status as the `uint` a compute kernel writes into its per-point status buffer.
 *
 * It is the enumerator's own value, so the kernel needs no table: the device-side code assigns the
 * same small integers by hand and this function is the host-side statement of that contract.
 *
 * @param s the status.
 * @return the code, `0..status_count-1` for every enumerator.
 * @complexity O(1).
 * @alloc none.
 * @test IrbemStatus.StatusCodeRoundTrips
 */
[[nodiscard]] constexpr std::uint32_t status_code(Status s) { return static_cast<std::uint32_t>(s); }

/**
 * The inverse of @ref status_code, for reading a device-written status buffer back.
 *
 * Validated rather than cast, because the bytes come from a buffer a driver wrote: an uninitialised
 * slot, a kernel that never ran, or a device-side bug produces an integer with no meaning, and
 * `static_cast`ing it into a `Status` would be undefined behaviour dressed up as a result.
 *
 * @param code the raw value read out of the status buffer.
 * @return the corresponding status, or `std::nullopt` when @p code names no enumerator.
 * @complexity O(1).
 * @alloc none.
 * @test IrbemStatus.StatusCodeRejectsAnUnknownCode
 */
[[nodiscard]] constexpr std::optional<Status> status_from_code(std::uint32_t code) {
    if (code >= status_count) return std::nullopt;
    return static_cast<Status>(code);
}

// -------------------------------------------------------------------------------------------
// Result
// -------------------------------------------------------------------------------------------

/**
 * A computed value and the reason it may not be trustworthy.
 *
 * An aggregate on purpose: no constructors, no invariants, no accessors beyond @ref ok. That is
 * what keeps it trivially copyable — the property the GPU lane needs, since a `Result<float>` array
 * is memcpy'd out of a device buffer — and what lets it be built with `{Status::Ok, x}` at every
 * one of the several hundred call sites this library will eventually have.
 *
 * The `status` member is first so that the layout is the same one the device writes: a small
 * integer tag followed by the payload at its natural alignment.
 *
 * @tparam T the value type; anything trivially copyable, in practice `double`, `float`, or a
 *           `cheatah::fixarray` vector.
 * @test IrbemStatus.ResultIsTriviallyCopyableAndSmall
 */
template <class T>
struct Result {
    /// Why the value may not be trustworthy; @ref Status::Ok when it may.
    Status status;
    /// The computed value — **always** populated, including when `status` is
    /// @ref Status::OutOfValidityRange. See the file brief.
    T value;

    /**
     * Whether the value carries no caveat.
     * @return `true` only when `status` is @ref Status::Ok.
     * @complexity O(1).
     * @alloc none.
     * @test IrbemStatus.ResultOkTracksItsStatus
     */
    [[nodiscard]] constexpr bool ok() const { return is_ok(status); }
};

/// A `Result` is exactly its payload plus one alignment unit of tag — the size claim the GPU lane
/// rests on, checked here for the three payloads that actually cross the device boundary.
static_assert(std::is_trivially_copyable_v<Result<double>>);
static_assert(std::is_trivially_copyable_v<Result<float>>);
static_assert(std::is_trivially_copyable_v<Result<Status>>);
static_assert(sizeof(Result<double>) == sizeof(double) + alignof(double));
static_assert(sizeof(Result<float>) == sizeof(float) + alignof(float));
static_assert(std::is_standard_layout_v<Result<double>>);
static_assert(offsetof(Result<double>, status) == 0,
              "the device writes the tag first; see the wire encoding above");

// -------------------------------------------------------------------------------------------
// The baddata bridge — the ONLY place -1e31 appears
// -------------------------------------------------------------------------------------------

/**
 * IRBEM's "no answer" sentinel, `-1e31`.
 *
 * Kept for one reason: a caller coming from IRBEM's C, Fortran, IDL or Python bindings tests
 * against this number, and the C-compatible boundary of this library must not break them. Note that
 * `1e31` is not exactly representable in binary64 (10^31 = 2^31 * 5^31, and 5^31 needs 73
 * significant bits), so this constant is the nearest `double` to -10^31 — the same rounding
 * IRBEM's own `-1d31` literal gets, which is why an equality test against it is nonetheless exact.
 *
 * @test IrbemStatus.BaddataIsTheIrbemSentinel
 */
inline constexpr double baddata = -1e31;

/**
 * Collapse a @ref Result down to IRBEM's convention, for the C-compatible boundary only.
 *
 * This is a lossy, one-way bridge: six statuses map onto one number and the reason is gone. There
 * is deliberately no inverse — a caller who wants to know *why* keeps the @ref Result.
 *
 * Note that @ref Status::OutOfValidityRange collapses to the sentinel here even though the value
 * was computed, because IRBEM's convention has no way to say "here is the number, with a caveat".
 * That loss is precisely the cost of the boundary, and precisely why nothing inside this library
 * crosses it.
 *
 * @param r the result to collapse.
 * @return `r.value` widened to `double` when `r.ok()`, otherwise @ref baddata.
 * @tparam T the payload type; must convert to `double`.
 * @complexity O(1).
 * @alloc none.
 * @test IrbemStatus.ToBaddataCollapsesEveryFailure
 */
template <std::convertible_to<double> T>
[[nodiscard]] constexpr double to_baddata(const Result<T>& r) {
    return r.ok() ? static_cast<double>(r.value) : baddata;
}

/**
 * Whether @p v is the sentinel, for code reading an IRBEM-convention array back.
 *
 * Exact equality, matching IRBEM's own `.eq. baddata` test. A NaN is *not* the sentinel — it
 * compares unequal to everything, this function included — which is correct: a NaN that escaped a
 * computation is a different failure from one the library deliberately reported.
 *
 * @param v the value read back.
 * @return `true` only for the bit pattern @ref baddata names.
 * @complexity O(1).
 * @alloc none.
 * @test IrbemStatus.IsBaddataDetectsOnlyTheSentinel
 */
[[nodiscard]] constexpr bool is_baddata(double v) { return v == baddata; }

// -------------------------------------------------------------------------------------------
// The external field models
// -------------------------------------------------------------------------------------------

/**
 * An external magnetospheric field model — IRBEM's `kext` key, with the same numbering.
 *
 * The numbering is IRBEM's because it is a published API: a caller migrating from IRBEM passes the
 * integer it already has, and @ref envelope_of indexes the envelope table with it directly. The
 * table is IRBEM `docs/source/api/general_information.rst`, "External magnetic field model".
 *
 * @test IrbemStatus.ModelNamesCoverEveryKey
 */
enum class ExternalModel : std::uint8_t {
    None = 0,                 ///< No external field — the internal model alone.
    MeadFairfield1975,        ///< Mead & Fairfield, JGR 80, 523 (1975).
    Tsyganenko1987Short,      ///< Tsyganenko, Planet. Space Sci. 35, 1347 (1987), short version.
    Tsyganenko1987Long,       ///< Tsyganenko, Planet. Space Sci. 35, 1347 (1987), long version.
    Tsyganenko1989,           ///< Tsyganenko, Planet. Space Sci. 37, 5 (1989) — "T89c".
    OlsonPfitzerQuiet1977,    ///< Olson & Pfitzer (1977), quiet-time model.
    OlsonPfitzerDynamic1988,  ///< Olson & Pfitzer (1988), dynamic model.
    Tsyganenko1996,           ///< Tsyganenko & Stern, JGR 101, 27187 (1996) — "T96".
    OstapenkoMaltsev1997,     ///< Ostapenko & Maltsev, Geomagn. Aeron. 37 (1997).
    Tsyganenko2001,           ///< Tsyganenko, JGR 107, 1179 (2002) — "T01".
    Tsyganenko2001Storm,      ///< Tsyganenko et al., JGR 108, 1209 (2003) — "T01s".
    Tsyganenko2004Storm,      ///< Tsyganenko & Sitnov (2005), doi:10.1029/2004JA010798 — "TS05".
    Alexeev2000,              ///< Alexeev et al. (2000), the paraboloid model.
    Tsyganenko2007,           ///< Tsyganenko & Sitnov (2007), doi:10.1029/2007JA012260 — "TS07D".
    MeadTsyganenko,           ///< IRBEM's own Mead-form fit to T89 (`kext = 14`).
};

/// How many `kext` keys the envelope table covers — `0..14`, which is IRBEM's whole range.
inline constexpr std::size_t model_count = 15;

/**
 * Whether @p m names one of the @ref model_count keys the envelope table covers.
 *
 * @param m the model key, possibly a raw integer cast in from a C caller.
 * @return `true` when `m` is a table index.
 * @complexity O(1).
 * @alloc none.
 * @test IrbemStatus.UnrecognisedModelIsADomainError
 */
[[nodiscard]] constexpr bool is_recognised(ExternalModel m) {
    return static_cast<std::size_t>(m) < model_count;
}

/**
 * The model's short name, as its papers and IRBEM's table spell it.
 * @param m the model.
 * @return a static string such as `"T96"`; `"?"` for an unrecognised key.
 * @complexity O(1).
 * @alloc none.
 * @test IrbemStatus.ModelNamesCoverEveryKey
 */
[[nodiscard]] constexpr std::string_view name_of(ExternalModel m) {
    switch (m) {
        case ExternalModel::None:
            return "none";
        case ExternalModel::MeadFairfield1975:
            return "MF75";
        case ExternalModel::Tsyganenko1987Short:
            return "T87short";
        case ExternalModel::Tsyganenko1987Long:
            return "T87long";
        case ExternalModel::Tsyganenko1989:
            return "T89";
        case ExternalModel::OlsonPfitzerQuiet1977:
            return "OP77quiet";
        case ExternalModel::OlsonPfitzerDynamic1988:
            return "OP88dyn";
        case ExternalModel::Tsyganenko1996:
            return "T96";
        case ExternalModel::OstapenkoMaltsev1997:
            return "OM97";
        case ExternalModel::Tsyganenko2001:
            return "T01";
        case ExternalModel::Tsyganenko2001Storm:
            return "T01storm";
        case ExternalModel::Tsyganenko2004Storm:
            return "TS05";
        case ExternalModel::Alexeev2000:
            return "A2000";
        case ExternalModel::Tsyganenko2007:
            return "TS07D";
        case ExternalModel::MeadTsyganenko:
            return "MeadT";
    }
    return "?";  // an unrecognised key; keeps this total
}

// -------------------------------------------------------------------------------------------
// Validity envelopes
// -------------------------------------------------------------------------------------------

/// The largest number of drivers any one model reads — TS05 reads ten (Dst, Pdyn, By, Bz, W1..W6).
inline constexpr std::size_t max_model_bounds = 10;

/**
 * One driver a model reads, and the inclusive interval the model was fitted over.
 *
 * **The interval is closed.** The published tables state their limits with `<=` and `>=`
 * ("-100 <= Dst <= 20"), so a driver sitting exactly on a bound is inside the envelope and
 * @ref check_validity returns @ref Status::Ok for it. Only a value strictly beyond a bound is
 * reported. The bounds are decimal values that are mostly not exactly representable in binary64
 * (0.5, 10 and 20 are; -100 is; 5 is), but that is irrelevant to the comparison: the same rounded
 * `double` is on both sides of it.
 *
 * An entry whose `lo` is `-inf` and `hi` is `+inf` means **the model reads this driver and no
 * bound has been published for it**. That is a fact, not a placeholder: IRBEM's `kext` table says
 * of the T01-storm and TS05 rows, in those words, "there is no upper or lower limit for those
 * inputs". Such an entry still gets the finiteness check, and never reports
 * @ref Status::OutOfValidityRange.
 *
 * @test IrbemStatus.EnvelopeTableMatchesTheIrbemKextTable
 */
struct DriverBound {
    /// Which `maginput` slot; @ref Driver's enumerators are the zero-based subscripts.
    Driver driver;
    /// Inclusive lower limit, or `-infinity` when unpublished.
    double lo;
    /// Inclusive upper limit, or `+infinity` when unpublished.
    double hi;

    /**
     * Whether a finite fitting range has been published for this driver.
     * @return `false` for the "used, unbounded" entries described above.
     * @complexity O(1).
     * @alloc none.
     * @test IrbemStatus.TS05DriversAreUsedButUnbounded
     */
    [[nodiscard]] constexpr bool published() const {
        return std::isfinite(lo) && std::isfinite(hi);
    }
};

/**
 * Everything published about where one model may be believed: which drivers it reads and over what
 * range, how far out in space it was fitted, and whether it needs coefficient files on disk.
 *
 * Aggregate, trivially copyable, and about 300 bytes — it lives in `.rodata` and is reached by
 * reference (@ref envelope_of), never copied into a hot loop.
 *
 * @test IrbemStatus.EnvelopeTableMatchesTheIrbemKextTable
 */
struct ValidityEnvelope {
    /// The source the ranges below are quoted from, so a reader can check them.
    std::string_view citation;
    /// The drivers this model reads; only the first `bound_count` entries are meaningful.
    std::array<DriverBound, max_model_bounds> bound_storage;
    /// How many entries of @ref bound_storage are meaningful.
    std::size_t bound_count;
    /// Largest geocentric radius, in Earth radii, the model was fitted out to; `+infinity` when
    /// the published sources state no radial limit.
    double max_r_geo;
    /// Smallest GSM x, in Earth radii, the model was fitted to; `-infinity` when none is published.
    /// The Tsyganenko 2001 and later models state this instead of a radial limit, because their
    /// tail current sheet is what runs out.
    double min_x_gsm;
    /// Whether the model needs per-interval coefficient files provisioned on disk before it can be
    /// evaluated at all; see @ref check_parameters.
    bool needs_coefficient_files;

    /**
     * The meaningful part of @ref bound_storage.
     * @return a pointer to the first bound; never null, and valid for @ref bound_count entries.
     * @complexity O(1).
     * @alloc none.
     * @test IrbemStatus.EnvelopeTableMatchesTheIrbemKextTable
     */
    [[nodiscard]] constexpr const DriverBound* begin() const { return bound_storage.data(); }

    /**
     * One past the meaningful part of @ref bound_storage, so an envelope is range-`for`-able.
     * @return the end iterator.
     * @complexity O(1).
     * @alloc none.
     * @test IrbemStatus.EnvelopeTableMatchesTheIrbemKextTable
     */
    [[nodiscard]] constexpr const DriverBound* end() const {
        return bound_storage.data() + bound_count;
    }

    /**
     * The bound this model publishes for @p d, if it reads @p d at all.
     * @param d the driver to look up.
     * @return the entry, or `std::nullopt` when the model does not read that driver.
     * @complexity O(@ref bound_count) — at most ten compares, over data in one cache line pair.
     * @alloc none.
     * @test IrbemStatus.EnvelopeLookupFindsOnlyTheDriversAModelReads
     */
    [[nodiscard]] constexpr std::optional<DriverBound> bound_for(Driver d) const {
        for (const DriverBound& b : *this) {
            if (b.driver == d) return b;
        }
        return std::nullopt;
    }
};

/// Positive infinity, spelled once so the table below reads as a table.
inline constexpr double unbounded_above = std::numeric_limits<double>::infinity();
/// Negative infinity, likewise.
inline constexpr double unbounded_below = -std::numeric_limits<double>::infinity();

/**
 * Build a @ref ValidityEnvelope from a braced list of bounds, so the table below cannot get its
 * count out of step with its contents.
 *
 * @param citation the source the ranges are quoted from.
 * @param bounds the drivers the model reads, at most @ref max_model_bounds of them; entries beyond
 *               that are dropped, which the table's own `static_assert` makes unreachable.
 * @param max_r_geo the published radial limit in Earth radii, or @ref unbounded_above.
 * @param min_x_gsm the published GSM-x limit in Earth radii, or @ref unbounded_below.
 * @param needs_coefficient_files whether coefficient files must be provisioned first.
 * @return the envelope.
 * @complexity O(`bounds.size()`).
 * @alloc none — the storage is the returned object's own inline array.
 * @test IrbemStatus.MakeEnvelopeFillsTheCountFromTheList
 */
[[nodiscard]] constexpr ValidityEnvelope make_envelope(std::string_view citation,
                                                       std::initializer_list<DriverBound> bounds,
                                                       double max_r_geo, double min_x_gsm,
                                                       bool needs_coefficient_files) {
    ValidityEnvelope e{citation, {}, 0, max_r_geo, min_x_gsm, needs_coefficient_files};
    for (const DriverBound& b : bounds) {
        if (e.bound_count == max_model_bounds) break;
        e.bound_storage[e.bound_count] = b;
        ++e.bound_count;
    }
    return e;
}

/**
 * Kp's bounds as the `maginput` vector carries it.
 *
 * Every Kp-driven model publishes its range as `0 <= Kp <= 9`, but IRBEM's slot 1 holds **Kp x 10**
 * (`general_information.rst`, "Magnetic field inputs": "consistent with OMNI2, this is Kp*10, and
 * it is in the range 0 to 90"). The scaling is applied once, here, so no model row has to remember
 * it.
 *
 * @test IrbemStatus.KpBoundsAreInOmniScaling
 */
inline constexpr DriverBound kp_bound{Driver::Kp, 0.0, 90.0};

/**
 * The envelope table, indexed by `kext`.
 *
 * Every number is quoted from IRBEM `docs/source/api/general_information.rst`, "External magnetic
 * field model", which restates each model's published range; the per-row citation names the paper
 * that range comes from. Rows whose drivers are listed with infinite bounds are the ones that table
 * explicitly declares unbounded, or for which the row states a driver list and no limits at all.
 *
 * **Known gaps, stated rather than filled in:**
 * - `OM97` (kext 8) and `A2000` (kext 12): the table names the drivers each reads and publishes no
 *   ranges, and no ranges are given in the sources this library is allowed to read. Their entries
 *   are therefore all-unbounded.
 * - `T01storm` (kext 10) and `TS05` (kext 11): the table says in so many words that "there is no
 *   upper or lower limit for those inputs".
 * - `TS07D` (kext 13): the table gives neither drivers nor limits. The model's own paper
 *   (Tsyganenko & Sitnov 2007) parameterizes it by `Pdyn` plus a per-interval coefficient set, so
 *   `Pdyn` is listed unbounded and `needs_coefficient_files` is set; IRBEM provisions those files
 *   with its `setup_ts07d_files.sh`, which fetches a `Coeffs/` and a `TAIL_PAR/` directory.
 * - `MeadT` (kext 14): an IRBEM-specific refit with no separate publication, so it inherits only
 *   the Kp range its Mead form is defined over and states no spatial limit.
 *
 * @test IrbemStatus.EnvelopeTableMatchesTheIrbemKextTable
 */
inline constexpr std::array<ValidityEnvelope, model_count> validity_envelopes{{
    // 0 — no external field at all; nothing to be outside of.
    make_envelope("no external field", {}, unbounded_above, unbounded_below, false),
    // 1 — Mead & Fairfield [1975]: uses 0 <= Kp <= 9, valid for rGEO <= 17 Re.
    make_envelope("Mead & Fairfield, JGR 80, 523 (1975); IRBEM kext table", {kp_bound}, 17.0,
                  unbounded_below, false),
    // 2 — Tsyganenko short [1987]: uses 0 <= Kp <= 9, valid for rGEO <= 30 Re.
    make_envelope("Tsyganenko, Planet. Space Sci. 35, 1347 (1987); IRBEM kext table", {kp_bound},
                  30.0, unbounded_below, false),
    // 3 — Tsyganenko long [1987]: uses 0 <= Kp <= 9, valid for rGEO <= 70 Re.
    make_envelope("Tsyganenko, Planet. Space Sci. 35, 1347 (1987); IRBEM kext table", {kp_bound},
                  70.0, unbounded_below, false),
    // 4 — Tsyganenko [1989c]: uses 0 <= Kp <= 9, valid for rGEO <= 70 Re. The model is Kp-BINNED
    //     rather than continuous in Kp; see t89_kp_bin.
    make_envelope("Tsyganenko, Planet. Space Sci. 37, 5 (1989); IRBEM kext table", {kp_bound}, 70.0,
                  unbounded_below, false),
    // 5 — Olson & Pfitzer quiet [1977]: no drivers, valid for rGEO <= 15 Re.
    make_envelope("Olson & Pfitzer (1977); IRBEM kext table", {}, 15.0, unbounded_below, false),
    // 6 — Olson & Pfitzer dynamic [1988]: 5 <= Dsw <= 50, 300 <= Vsw <= 500, -100 <= Dst <= 20,
    //     valid for rGEO <= 60 Re.
    make_envelope("Olson & Pfitzer (1988); IRBEM kext table",
                  {{Driver::Dsw, 5.0, 50.0}, {Driver::Vsw, 300.0, 500.0}, {Driver::Dst, -100.0, 20.0}},
                  60.0, unbounded_below, false),
    // 7 — Tsyganenko [1996]: -100 <= Dst <= 20, 0.5 <= Pdyn <= 10, |By| <= 10, |Bz| <= 10,
    //     valid for rGEO <= 40 Re.
    make_envelope("Tsyganenko & Stern, JGR 101, 27187 (1996); IRBEM kext table",
                  {{Driver::Dst, -100.0, 20.0},
                   {Driver::Pdyn, 0.5, 10.0},
                   {Driver::ByIMF, -10.0, 10.0},
                   {Driver::BzIMF, -10.0, 10.0}},
                  40.0, unbounded_below, false),
    // 8 — Ostapenko & Maltsev [1997]: uses Dst, Pdyn, Bz, Kp; NO published ranges (see the gap list).
    make_envelope("Ostapenko & Maltsev, Geomagn. Aeron. 37 (1997); IRBEM kext table lists the "
                  "drivers and publishes no ranges",
                  {{Driver::Dst, unbounded_below, unbounded_above},
                   {Driver::Pdyn, unbounded_below, unbounded_above},
                   {Driver::BzIMF, unbounded_below, unbounded_above},
                   {Driver::Kp, unbounded_below, unbounded_above}},
                  unbounded_above, unbounded_below, false),
    // 9 — Tsyganenko [2001]: -50 <= Dst <= 20, 0.5 <= Pdyn <= 5, |By| <= 5, |Bz| <= 5,
    //     0 <= G1 <= 10, 0 <= G2 <= 10, valid for xGSM >= -15 Re.
    make_envelope("Tsyganenko, JGR 107, 1179 (2002); IRBEM kext table",
                  {{Driver::Dst, -50.0, 20.0},
                   {Driver::Pdyn, 0.5, 5.0},
                   {Driver::ByIMF, -5.0, 5.0},
                   {Driver::BzIMF, -5.0, 5.0},
                   {Driver::G1, 0.0, 10.0},
                   {Driver::G2, 0.0, 10.0}},
                  unbounded_above, -15.0, false),
    // 10 — Tsyganenko [2001] storm: uses Dst, Pdyn, By, Bz, G2, G3; "there is no upper or lower
    //      limit for those inputs"; valid for xGSM >= -15 Re.
    make_envelope("Tsyganenko et al., JGR 108, 1209 (2003); IRBEM kext table states no limits",
                  {{Driver::Dst, unbounded_below, unbounded_above},
                   {Driver::Pdyn, unbounded_below, unbounded_above},
                   {Driver::ByIMF, unbounded_below, unbounded_above},
                   {Driver::BzIMF, unbounded_below, unbounded_above},
                   {Driver::G2, unbounded_below, unbounded_above},
                   {Driver::G3, unbounded_below, unbounded_above}},
                  unbounded_above, -15.0, false),
    // 11 — Tsyganenko [2004] storm (TS05): uses Dst, Pdyn, By, Bz, W1..W6; "there is no upper or
    //      lower limit for those inputs"; valid for xGSM >= -15 Re.
    make_envelope("Tsyganenko & Sitnov, doi:10.1029/2004JA010798 (2005); IRBEM kext table states "
                  "no limits",
                  {{Driver::Dst, unbounded_below, unbounded_above},
                   {Driver::Pdyn, unbounded_below, unbounded_above},
                   {Driver::ByIMF, unbounded_below, unbounded_above},
                   {Driver::BzIMF, unbounded_below, unbounded_above},
                   {Driver::W1, unbounded_below, unbounded_above},
                   {Driver::W2, unbounded_below, unbounded_above},
                   {Driver::W3, unbounded_below, unbounded_above},
                   {Driver::W4, unbounded_below, unbounded_above},
                   {Driver::W5, unbounded_below, unbounded_above},
                   {Driver::W6, unbounded_below, unbounded_above}},
                  unbounded_above, -15.0, false),
    // 12 — Alexeev [2000] paraboloid: uses Dsw, Vsw, Dst, Bz, AL; no published ranges.
    make_envelope("Alexeev et al. (2000); IRBEM kext table lists the drivers and publishes no "
                  "ranges",
                  {{Driver::Dsw, unbounded_below, unbounded_above},
                   {Driver::Vsw, unbounded_below, unbounded_above},
                   {Driver::Dst, unbounded_below, unbounded_above},
                   {Driver::BzIMF, unbounded_below, unbounded_above},
                   {Driver::AL, unbounded_below, unbounded_above}},
                  unbounded_above, unbounded_below, false),
    // 13 — Tsyganenko [2007] (TS07D): coefficient files required; Pdyn is the only maginput driver.
    make_envelope("Tsyganenko & Sitnov, doi:10.1029/2007JA012260 (2007); IRBEM kext table gives no "
                  "ranges; coefficients per IRBEM setup_ts07d_files.sh",
                  {{Driver::Pdyn, unbounded_below, unbounded_above}}, unbounded_above,
                  unbounded_below, true),
    // 14 — Mead-Tsyganenko: IRBEM's Mead-form refit of T89; uses Kp, no separate publication.
    make_envelope("IRBEM kext table: \"onera model where the Tsyganenko 89 model is best fitted by "
                  "a Mead model\"",
                  {kp_bound}, unbounded_above, unbounded_below, false),
}};

/// The envelope handed back for a `kext` outside `0..14`: no drivers, no limits, and never reached
/// through @ref check_validity, which rejects an unrecognised key before it looks anything up.
inline constexpr ValidityEnvelope unknown_envelope =
    make_envelope("unrecognised kext", {}, unbounded_above, unbounded_below, false);

/**
 * The published envelope for @p m.
 *
 * @param m the model key.
 * @return a reference to the table row, or to @ref unknown_envelope for an unrecognised key. The
 *         referent has static storage duration and outlives every caller.
 * @complexity O(1) — one bounds check and an index.
 * @alloc none.
 * @test IrbemStatus.EnvelopeTableMatchesTheIrbemKextTable
 */
[[nodiscard]] constexpr const ValidityEnvelope& envelope_of(ExternalModel m) {
    if (!is_recognised(m)) return unknown_envelope;
    return validity_envelopes[static_cast<std::size_t>(m)];
}

// -------------------------------------------------------------------------------------------
// The checks
// -------------------------------------------------------------------------------------------

/**
 * Whether @p drivers put @p m outside the data it was fitted to.
 *
 * **This never suppresses a value.** It is called beside the evaluation, not instead of it: the
 * caller gets the number the functional form produces AND the knowledge that the number is an
 * extrapolation. Deciding what to do about that — clamp, refuse, or proceed and annotate — is a
 * scientific judgement this library will not make on the caller's behalf.
 *
 * Only the drivers the model actually reads are examined, so a `maginput` vector full of fill
 * values in the slots a model ignores is still `Ok` for that model. That matters in practice: a
 * caller running the same vector through T89 and T96 has, at most, the union of what both need.
 *
 * The two failures are ordered, most fundamental first: **every** driver the model reads is checked
 * for finiteness before **any** is checked against its range, so a NaN anywhere in the used set
 * reports @ref Status::DomainError rather than whichever came first in the list.
 *
 * @param m the model.
 * @param drivers the whole 25-slot `maginput` vector; the reserved slots are never read.
 * @return @ref Status::DomainError for an unrecognised key or a non-finite used driver,
 *         @ref Status::OutOfValidityRange for a used driver strictly outside its published closed
 *         interval, otherwise @ref Status::Ok.
 * @complexity O(@ref max_model_bounds) — at most twenty compares, no branches on data.
 * @alloc none.
 * @test IrbemStatus.EveryBoundedDriverIsCheckedFromBothSides
 */
[[nodiscard]] inline Status check_validity(ExternalModel m, const DriverSet& drivers) {
    if (!is_recognised(m)) return Status::DomainError;
    const ValidityEnvelope& env = envelope_of(m);
    for (const DriverBound& b : env) {
        if (!std::isfinite(drivers[static_cast<std::size_t>(b.driver)])) return Status::DomainError;
    }
    for (const DriverBound& b : env) {
        const double v = drivers[static_cast<std::size_t>(b.driver)];
        if (v < b.lo || v > b.hi) return Status::OutOfValidityRange;
    }
    return Status::Ok;
}

/**
 * The smallest geocentric radius that is unambiguously above the ground, in Earth radii: the WGS84
 * polar semi-minor axis over the equatorial semi-major axis, `6356.752314245 / 6378.137`.
 *
 * A point below this radius is inside the solid Earth at *every* latitude, so rejecting it needs no
 * geodetic conversion and cannot reject a legitimate low-altitude polar point by mistake. Between
 * this and 1 Re a point may be above ground (near the poles) or below it (near the equator);
 * deciding which is the geodetic module's job, and this check deliberately does not attempt it.
 *
 * WGS84 axes per NIMA TR8350.2, 3rd ed. (2000), §3.2.
 *
 * @test IrbemStatus.PositionInsideTheEarthIsADomainError
 */
inline constexpr double min_r_geo = 6356.752314245 / 6378.137;

/**
 * Whether a point is inside the region @p m was fitted over.
 *
 * The spatial envelope is stated two different ways by the two model families, and both are
 * checked: the pre-2001 models publish a maximum geocentric radius (they stop having a tail),
 * while T01 and later publish a minimum GSM x (their tail current sheet is what runs out).
 * A model publishing neither accepts every finite point.
 *
 * @param m the model.
 * @param r_geo the geocentric radius, in Earth radii.
 * @param x_gsm the GSM x coordinate, in Earth radii.
 * @return @ref Status::DomainError for an unrecognised key, a non-finite coordinate, or a radius
 *         below @ref min_r_geo; @ref Status::OutOfValidityRange outside the published envelope;
 *         otherwise @ref Status::Ok.
 * @complexity O(1).
 * @alloc none.
 * @test IrbemStatus.PositionEnvelopeIsCheckedFromBothSides
 */
[[nodiscard]] inline Status check_position(ExternalModel m, double r_geo, double x_gsm) {
    if (!is_recognised(m)) return Status::DomainError;
    if (!std::isfinite(r_geo) || !std::isfinite(x_gsm)) return Status::DomainError;
    if (r_geo < min_r_geo) return Status::DomainError;
    const ValidityEnvelope& env = envelope_of(m);
    if (r_geo > env.max_r_geo || x_gsm < env.min_x_gsm) return Status::OutOfValidityRange;
    return Status::Ok;
}

/**
 * Whether @p m can be evaluated at all with the coefficient files currently provisioned.
 *
 * Only TS07D needs them: it is not a closed-form fit but a per-interval expansion whose radial
 * basis function coefficients are downloaded per six-minute interval (IRBEM ships
 * `setup_ts07d_files.sh` to fetch a `Coeffs/` and a `TAIL_PAR/` tree for exactly this reason).
 * Without them the model has no parameters, which is a different failure from having bad ones —
 * hence its own status rather than @ref Status::DomainError.
 *
 * @param m the model.
 * @param coefficients_present whether the caller has located the files; the *locating* is the
 *                             caller's business, so that this header stays free of filesystem I/O.
 * @return @ref Status::DomainError for an unrecognised key, @ref Status::ParametersMissing when the
 *         model needs files the caller does not have, otherwise @ref Status::Ok.
 * @complexity O(1).
 * @alloc none.
 * @test IrbemStatus.ParametersAreMissingOnlyForTs07d
 */
[[nodiscard]] constexpr Status check_parameters(ExternalModel m, bool coefficients_present) {
    if (!is_recognised(m)) return Status::DomainError;
    if (envelope_of(m).needs_coefficient_files && !coefficients_present) {
        return Status::ParametersMissing;
    }
    return Status::Ok;
}

/**
 * The Kp bin T89 uses, from Kp in IRBEM's OMNI2 scaling.
 *
 * T89 is not continuous in Kp. Tsyganenko (Planet. Space Sci. 37, 5, 1989) sorted the data into
 * seven Kp intervals and fitted a separate coefficient set to each: `{0, 0+}`, `{1-, 1, 1+}`,
 * `{2-, 2, 2+}`, `{3-, 3, 3+}`, `{4-, 4, 4+}`, `{5-, 5, 5+}`, and `{>= 6-}`. The returned index is
 * the 1-based bin number, which is the `iopt` argument the published T89 interface takes.
 *
 * The thresholds therefore fall at Kp x 10 = 5, 15, 25, 35, 45, 55 — **between** the values Kp can
 * actually take, since the third-of-a-unit steps land on 0, 3, 7, 10, 13, 17, ... in this scaling.
 * There is consequently no "exact boundary" case for a real Kp to sit on, and the bin edges are
 * chosen midway rather than at a bin's own extreme so that the mapping is robust to a caller who
 * has already rounded. The convention where it cannot matter is nonetheless decided and uniform:
 * **a threshold value belongs to the bin below it**, so the bins are the half-open intervals
 * `(-inf, 5]`, `(5, 15]`, ... `(55, +inf)` in this scaling. A negative Kp — which
 * @ref check_validity reports separately — still yields bin 1, and any Kp above 5.5 yields bin 7.
 *
 * @param kp_times_ten Kp in IRBEM's slot-1 scaling, i.e. Kp x 10, nominally 0..90.
 * @return the bin, `1..7`.
 * @complexity O(1).
 * @alloc none.
 * @test IrbemStatus.T89KpBinsFollowThePublishedIntervals
 */
[[nodiscard]] inline int t89_kp_bin(double kp_times_ten) {
    if (!(kp_times_ten > 5.0)) return 1;  // the negation also catches NaN and every negative Kp
    // ceil, not floor: a value sitting exactly on a threshold belongs to the bin below it. The
    // divisions are exact — every threshold minus 5 is a multiple of 10, and 10 is not a power of
    // two, but 10/10, 20/10, ... 50/10 are all exactly representable quotients.
    const double index = std::ceil((kp_times_ten - 5.0) / 10.0);
    if (!(index < 6.0)) return 7;  // saturates, and catches +infinity
    return 1 + static_cast<int>(index);
}

}  // namespace cheatah::space::irbem
