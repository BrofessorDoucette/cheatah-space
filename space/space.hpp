#pragma once
// cheatah-deps: ndarray
//
// space.cdf hands variables back as ndarrays, so the package always needs it. purrc reads
// this marker from the header it RESOLVES — `import space.cdf` resolves module `space`,
// i.e. this file — and puts ndarray on the include and link paths for every consumer.

/**
 * @file space.hpp
 * @brief The `space` package — cheatah-space's astronomy / space-physics extension.
 *
 * `import space` pulls in every space.* submodule; `import space.time` (etc.) pulls in just
 * one. The package header simply includes the submodules so the `cheatah::space::*`
 * namespaces are defined, mirroring the first-party `parsers` package.
 *
 * Submodules:
 *   - space.time  — Julian dates & epoch conversions (the CDF_EPOCH bridge).   [working]
 *   - space.cdf   — NASA Common Data Format I/O.                                [working]
 *   - space.irbem — radiation-belt / magnetic-field models (IRBEM reimpl).      [working]
 */

#include "time/time.hpp"
#include "cdf/cdf.hpp"

// space.irbem's geometry is `cheatah::fixarray`, so it is included ONLY when fixarray is on the
// include path — exactly the opt-in-by-include-path idiom space/time/time.hpp already uses for
// ndarray. Without the guard, `import space.time` would fail for every program that did not also
// import fixarray, which is a dependency space.time does not have and should not acquire by being
// in the same package as something that does.
//
// A program wanting the physics writes `import fixarray` alongside `import space.irbem`, or simply
// includes `space/irbem/irbem.hpp` directly.
#if __has_include("fixarray.hpp")
#  include "irbem/irbem.hpp"
#endif

namespace cheatah::space {}
