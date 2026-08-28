#pragma once

/**
 * @file cdf.hpp
 * @brief space.cdf — NASA Common Data Format I/O, written from scratch.
 *
 * The format essentially all heliophysics data ships in — RBSP, MMS, THEMIS, OMNI, Voyager.
 * `import space.cdf` resolves this header, which pulls in the submodule headers so the
 * `cheatah::space::cdf` namespace is complete.
 *
 * ZERO DEPENDENCIES. NASA's CDF library is neither linked nor required to build, install or
 * use space.cdf. We implement from the published
 * [Internal Format Description](https://spdf.gsfc.nasa.gov/pub/software/cdf/doc/cdf39/cdf39ifd.pdf);
 * their distribution is used only as an optional, dev-only differential oracle, fetched by
 * `scripts/cdf-oracle.sh` into the git-ignored `space/cdf/vendor/`. See the module README.
 *
 * What is here so far — this module is being built up in order, and the file lists what
 * actually exists rather than what is planned:
 *
 *   - `types.hpp`        the format vocabulary: data types, encodings, record types, errors.
 *   - `bytes.hpp`        the one bounds-checked, byte-swapping reader over a file's bytes.
 *   - `mapping.hpp`      read-only memory mapping, with the syscalls behind a testable seam.
 *   - `records.hpp`      the internal records (CDR/GDR/VDR/VXR/VVR/CVVR/CPR), parsed and validated.
 *   - `index.hpp`        a variable's VXR tree flattened, iteratively and capped, into extents.
 *   - `leapseconds.hpp`  NASA's leap-second table as exact integers, plus TAI-UTC lookup.
 *
 * The reader, writer, compression codecs, checksum and signing land on top of this.
 */

#include "bytes.hpp"
#include "leapseconds.hpp"
#include "index.hpp"
#include "mapping.hpp"
#include "records.hpp"
#include "types.hpp"

namespace cheatah::space::cdf {}
