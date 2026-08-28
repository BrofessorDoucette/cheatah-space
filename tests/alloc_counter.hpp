#pragma once

/// @file alloc_counter.hpp
/// @brief A process-wide heap tripwire for the space.irbem tests.
///
/// "Nothing on the heap in a hot path" is a claim, so it gets a counter rather than a comment.
/// The global `operator new` family is replaced in alloc_counter.cpp so nothing can allocate
/// behind our back — not the routine under test, not a container it constructs internally.
///
/// **This lives in exactly ONE translation unit on purpose.** Replacing `operator new` is a
/// whole-program act: three test files independently rolled their own counter, and linking them
/// together produced `multiple definition of operator new`. One definition, one counter, shared.
///
/// Usage — the assertion that actually catches per-call allocation is the SECOND call, not the
/// first. A routine that allocates a workspace on every invocation passes a single-call check and
/// fails this one:
/// @code
/// run_once(...);                                        // may populate a workspace
/// const std::size_t before = allocation_count();
/// run_once(...);                                        // must reuse it
/// EXPECT_EQ(before, allocation_count());
/// @endcode

#include <cstddef>

namespace cheatah_space_test {

/// Total global allocations since process start. Monotonic; compare two reads to bound a region.
[[nodiscard]] std::size_t allocation_count() noexcept;

}  // namespace cheatah_space_test
