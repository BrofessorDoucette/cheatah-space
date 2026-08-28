// The single definition of the test suite's heap tripwire. See alloc_counter.hpp for why this is
// one translation unit and not one per test.
//
// EVERY global allocation form is replaced, not just the common one: if `operator new(size_t)` is
// replaced but the aligned or array form is not, an allocation can slip past the counter, and —
// worse — a deallocation can be routed to a different implementation than its allocation. The
// deletes are replaced for exactly that pairing reason, and deliberately do not count: only
// allocations are of interest, and counting frees would make the assertion order-dependent.
#include "alloc_counter.hpp"

#include <cstdlib>
#include <new>

namespace {

std::size_t g_allocations = 0;  // NOLINT(cppcoreguidelines-avoid-non-const-global-variables)

/// Count, then allocate. `bytes == 0` is bumped to 1 because a zero-size request must still return
/// a unique pointer.
void* tracked_allocate(std::size_t bytes) noexcept {
    ++g_allocations;
    return std::malloc(bytes == 0 ? 1 : bytes);
}

void* tracked_allocate_aligned(std::size_t bytes, std::size_t alignment) noexcept {
    ++g_allocations;
    // aligned_alloc requires the size to be a multiple of the alignment.
    const std::size_t rounded = ((bytes == 0 ? 1 : bytes) + alignment - 1) / alignment * alignment;
    return std::aligned_alloc(alignment, rounded);
}

}  // namespace

namespace cheatah_space_test {

std::size_t allocation_count() noexcept { return g_allocations; }

}  // namespace cheatah_space_test

void* operator new(std::size_t bytes) {
    void* p = tracked_allocate(bytes);
    if (p == nullptr) throw std::bad_alloc();
    return p;
}
void* operator new[](std::size_t bytes) {
    void* p = tracked_allocate(bytes);
    if (p == nullptr) throw std::bad_alloc();
    return p;
}
void* operator new(std::size_t bytes, const std::nothrow_t& /*tag*/) noexcept {
    return tracked_allocate(bytes);
}
void* operator new[](std::size_t bytes, const std::nothrow_t& /*tag*/) noexcept {
    return tracked_allocate(bytes);
}
void* operator new(std::size_t bytes, std::align_val_t alignment) {
    void* p = tracked_allocate_aligned(bytes, static_cast<std::size_t>(alignment));
    if (p == nullptr) throw std::bad_alloc();
    return p;
}
void* operator new[](std::size_t bytes, std::align_val_t alignment) {
    void* p = tracked_allocate_aligned(bytes, static_cast<std::size_t>(alignment));
    if (p == nullptr) throw std::bad_alloc();
    return p;
}

void operator delete(void* p) noexcept { std::free(p); }
void operator delete[](void* p) noexcept { std::free(p); }
void operator delete(void* p, std::size_t /*bytes*/) noexcept { std::free(p); }
void operator delete[](void* p, std::size_t /*bytes*/) noexcept { std::free(p); }
void operator delete(void* p, const std::nothrow_t& /*tag*/) noexcept { std::free(p); }
void operator delete[](void* p, const std::nothrow_t& /*tag*/) noexcept { std::free(p); }
void operator delete(void* p, std::align_val_t /*alignment*/) noexcept { std::free(p); }
void operator delete[](void* p, std::align_val_t /*alignment*/) noexcept { std::free(p); }
void operator delete(void* p, std::size_t /*bytes*/, std::align_val_t /*alignment*/) noexcept {
    std::free(p);
}
void operator delete[](void* p, std::size_t /*bytes*/, std::align_val_t /*alignment*/) noexcept {
    std::free(p);
}
