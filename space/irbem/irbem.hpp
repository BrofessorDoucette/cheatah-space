#pragma once

/**
 * @file irbem.hpp
 * @brief space.irbem — the package umbrella; `import space.irbem` resolves here.
 *
 * Radiation-belt and magnetic-field models: a from-scratch C++20 reimplementation of
 * [PRBEM/IRBEM](https://github.com/PRBEM/IRBEM), written to the published papers rather than
 * derived from its source. IRBEM is LGPL-3.0 and Tsyganenko's own reference code is GPL-3.0;
 * cheatah-space is MIT, so neither is a safe thing to derive from. IRBEM is used only as an
 * optional, dev-time black-box oracle, never linked and never shipped.
 *
 * The reason this module exists is speed. IRBEM is single-threaded, unvectorized, and ships
 * compiled with no `-O` flag at all; one L\* evaluation costs ~10^5 magnetic-field model calls, and
 * the field-line and drift-shell integrals that dominate it are embarrassingly parallel. They run
 * on the GPU here — see gpu/dispatch.hpp — with the reductions ordered on the host.
 *
 * Layers, innermost first:
 *
 *  - @ref frames.hpp — the reference frames, with the frame in the TYPE, so handing a GSM vector to
 *    a routine expecting GEO is a compile error rather than a plausible wrong answer;
 *  - @ref policy.hpp — `Precision` and `Compat` as compile-time policies;
 *  - @ref context.hpp — `FieldContext`, immutable, replacing IRBEM's ~20 mutable `COMMON` blocks
 *    (that global state is precisely what makes the reference impossible to thread or offload);
 *  - @ref datetime.hpp, @ref coords_rotations.hpp, @ref coords_geodetic.hpp, @ref coords_helio.hpp;
 *  - @ref igrf.hpp — the internal field, and the hottest kernel in the library;
 *  - @ref ext_t89.hpp — Tsyganenko (1989), the first EXTERNAL field model: the magnetosphere's own
 *    currents, without which a drift shell traced through IGRF alone is a dipole exercise;
 *  - @ref gpu/dispatch.hpp — the device seam, present only when the GPU stack is on the include
 *    path and reporting itself unavailable rather than failing when it is not.
 *
 * Accuracy is a published budget, not an accident: see docs/ERROR_BUDGET.md, whose numbers are
 * measured against the oracle rather than asserted.
 */

#include "frames.hpp"
#include "policy.hpp"
#include "context.hpp"
#include "library_info.hpp"
#include "coords_rotations.hpp"
#include "coords_geodetic.hpp"
#include "coords_helio.hpp"
#include "igrf.hpp"
#include "ext_t89.hpp"
#include "ext_mead.hpp"
#include "ext_opq.hpp"
#include "ext_opd.hpp"
#include "ext_ostapenko.hpp"
#include "field.hpp"
#include "batch_soa.hpp"
#include "total_field.hpp"
#include "lstar.hpp"
#include "driftshell.hpp"
#include "trace_api.hpp"
#include "purr.hpp"
#include "api.hpp"

// The device lane is opt-in by include path, exactly as space.time's ndarray support is: a program
// built without cheatah-gpu-linalg still compiles every routine here and simply has no GPU.
#if __has_include("cheatah_gpu_linalg/context.hpp")
#  include "gpu/dispatch.hpp"
#endif
