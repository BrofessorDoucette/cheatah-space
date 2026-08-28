#pragma once

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
 *   - space.cdf   — NASA Common Data Format I/O.                                [roadmap]
 *   - space.irbem — radiation-belt / magnetic-field models (IRBEM reimpl).      [building]
 */

#include "time/time.hpp"
#include "cdf/cdf.hpp"
#include "irbem/irbem.hpp"

namespace cheatah::space {}
