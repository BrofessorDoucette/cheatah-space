#pragma once

/**
 * @file mapping.hpp
 * @brief space.cdf — memory-mapping a file, with the syscalls behind a testable seam.
 *
 * The only file in space.cdf that includes a platform header. Everything above it works on a
 * `Bytes` span and neither knows nor cares where the bytes came from, which is what lets the
 * hostile-input tests run entirely on heap buffers — see below for why that matters.
 *
 * **Why mmap rather than read().** NASA's library reads through a `FILE*` and an internal pool of
 * cache buffers (that is what `CDFsetCacheSize` and the `-cachesize` flags on their tools tune),
 * so every record access is a seek, a read into the pool, and a copy out of it. Mapping the file
 * once makes the operating system's page cache the buffer: no syscall per record, no pool, no
 * LRU, and no copy on the paths where the bytes are already in the layout we want.
 *
 * **Why the syscall seam.** `mmap` failing is not reachable from any input a test can craft — you
 * cannot write a malformed CDF that makes the kernel run out of address space. Templating on a
 * `SysOps` policy lets the tests instantiate `BasicFileMapping<FailingSysOps>` and drive each
 * failure branch deliberately. In production the default policy is a stateless empty struct whose
 * members are one-line wrappers around the real calls, so it inlines away to nothing: the seam
 * costs an empty base and zero instructions, and buys real coverage of the error paths that the
 * QA gate's 100%-line bar would otherwise make unreachable.
 *
 * **A note the fuzzing depends on.** ASan cannot see an over-read that stays inside an mmap'd
 * region — the pages are mapped, so nothing traps. It sees every byte of a heap buffer. That is
 * why `open_from_memory()` exists on the reader and why every hostile-input test goes through it
 * rather than through a file: a bounds bug that a mapped read would silently tolerate becomes a
 * hard ASan failure. Mapping is the fast path for trusted files; the heap is the path that is
 * actually being audited.
 */

#include "cheatah.hpp"

#include <cstddef>
#include <cstdint>
#include <string>
#include <utility>

#include "bytes.hpp"
#include "types.hpp"

#if defined(_WIN32)
#  include <windows.h>
#else
#  include <fcntl.h>
#  include <sys/mman.h>
#  include <sys/stat.h>
#  include <unistd.h>
#endif

namespace cheatah::space::cdf::detail {

/**
 * The syscalls a mapping needs, as a policy so the tests can fail them on demand.
 *
 * Stateless: every member is static, so `BasicFileMapping` stores nothing for it and the calls
 * inline straight through to the platform.
 */
struct PosixSysOps {
#if defined(_WIN32)
    /// Open @p path for reading. @return a handle, or nullptr on failure.
    static void* open_file(const std::string& path) {
        void* h = ::CreateFileA(path.c_str(), GENERIC_READ, FILE_SHARE_READ, nullptr,
                                OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL, nullptr);
        return h == INVALID_HANDLE_VALUE ? nullptr : h;
    }
    /// @param handle an open file. @return its size in bytes, or 0 on failure.
    static std::uint64_t file_size(void* handle) {
        LARGE_INTEGER size{};
        return ::GetFileSizeEx(handle, &size) ? static_cast<std::uint64_t>(size.QuadPart) : 0;
    }
    /// Map @p size bytes of @p handle read-only. @return the base address, or nullptr.
    static void* map_file(void* handle, std::uint64_t size) {
        void* m = ::CreateFileMappingA(handle, nullptr, PAGE_READONLY, 0, 0, nullptr);
        if (m == nullptr) { return nullptr; }
        void* view = ::MapViewOfFile(m, FILE_MAP_READ, 0, 0, static_cast<SIZE_T>(size));
        ::CloseHandle(m);
        return view;
    }
    /// Unmap a previously mapped view. @param addr base. @param size length.
    static void unmap_file(void* addr, std::uint64_t size) {
        (void)size;
        ::UnmapViewOfFile(addr);
    }
    /// Close an open file. @param handle the handle.
    static void close_file(void* handle) { ::CloseHandle(handle); }
    /// Advise the kernel that the whole mapping will be read. @param addr base. @param size length.
    static void advise_willneed(void* addr, std::uint64_t size) { (void)addr; (void)size; }
#else
    /// Open @p path for reading. @return a file descriptor, or -1 on failure.
    static int open_file(const std::string& path) { return ::open(path.c_str(), O_RDONLY); }

    /// @param fd an open descriptor. @return the file's size in bytes, or 0 on failure.
    static std::uint64_t file_size(int fd) {
        struct stat st {};
        if (::fstat(fd, &st) != 0 || st.st_size < 0) { return 0; }
        return static_cast<std::uint64_t>(st.st_size);
    }

    /// Map @p size bytes of @p fd read-only. @return the base address, or nullptr on failure.
    static void* map_file(int fd, std::uint64_t size) {
        void* addr = ::mmap(nullptr, static_cast<std::size_t>(size), PROT_READ, MAP_PRIVATE, fd, 0);
        return addr == MAP_FAILED ? nullptr : addr;
    }

    /// Unmap a mapping. @param addr base. @param size length.
    static void unmap_file(void* addr, std::uint64_t size) {
        ::munmap(addr, static_cast<std::size_t>(size));
    }

    /// Close a descriptor. @param fd the descriptor.
    static void close_file(int fd) { ::close(fd); }

    /// Tell the kernel the descriptors and index will be read, so it can fault them in ahead of
    /// use. Advisory only — a failure is not an error and is deliberately ignored.
    /// @param addr base. @param size length.
    static void advise_willneed(void* addr, std::uint64_t size) {
#  if defined(MADV_WILLNEED)
        ::madvise(addr, static_cast<std::size_t>(size), MADV_WILLNEED);
#  else
        (void)addr;
        (void)size;
#  endif
    }
#endif
};

/**
 * A read-only memory mapping of a file, released when it goes out of scope.
 *
 * Move-only: two owners would unmap twice. Templated on the syscall policy purely so the failure
 * branches are testable — see the file docs.
 *
 * @tparam Sys the syscall policy; PosixSysOps in production.
 */
template <class Sys = PosixSysOps>
class BasicFileMapping {
  public:
#if defined(_WIN32)
    using handle_t = void*;                              ///< Platform file handle.
    static constexpr handle_t kNoHandle = nullptr;       ///< Value meaning "not open".
#else
    using handle_t = int;                                ///< Platform file descriptor.
    static constexpr handle_t kNoHandle = -1;            ///< Value meaning "not open".
#endif

    /// An empty mapping, owning nothing.
    BasicFileMapping() = default;

    /**
     * Map @p path read-only.
     *
     * @param path the file to map.
     * @throws CdfError ErrorCode::CannotOpen when the file will not open, ErrorCode::EmptyFile
     *         when it is too small to be a CDF, ErrorCode::CannotMap when mapping fails.
     */
    explicit BasicFileMapping(const std::string& path) {
        handle_ = Sys::open_file(path);
        detail::require(handle_ != kNoHandle, ErrorCode::CannotOpen, 0);

        size_ = Sys::file_size(handle_);
        // Below the two magic numbers there is nothing that could be a CDF, and a zero-length
        // mmap is an error on POSIX — so this check also keeps the mapping call well-defined.
        if (size_ < kCdrOffset) {
            Sys::close_file(handle_);
            handle_ = kNoHandle;
            detail::throw_cdf(ErrorCode::EmptyFile, 0);
        }

        addr_ = Sys::map_file(handle_, size_);
        if (addr_ == nullptr) {
            Sys::close_file(handle_);
            handle_ = kNoHandle;
            size_ = 0;
            detail::throw_cdf(ErrorCode::CannotMap, 0);
        }

        // The descriptors and the variable index live near the front and are walked immediately;
        // asking for them up front turns a scatter of minor faults into one readahead.
        Sys::advise_willneed(addr_, size_ < kAdviseBytes ? size_ : kAdviseBytes);
    }

    /// Release the mapping.
    ~BasicFileMapping() { reset(); }

    BasicFileMapping(const BasicFileMapping&) = delete;
    BasicFileMapping& operator=(const BasicFileMapping&) = delete;

    /// Take ownership of @p other's mapping. @param other the mapping to move from.
    BasicFileMapping(BasicFileMapping&& other) noexcept
        : addr_(other.addr_), size_(other.size_), handle_(other.handle_) {
        other.addr_ = nullptr;
        other.size_ = 0;
        other.handle_ = kNoHandle;
    }

    /// Take ownership of @p other's mapping, releasing ours. @param other the mapping to move
    /// from. @return this mapping.
    BasicFileMapping& operator=(BasicFileMapping&& other) noexcept {
        if (this != &other) {
            reset();
            addr_ = other.addr_;
            size_ = other.size_;
            handle_ = other.handle_;
            other.addr_ = nullptr;
            other.size_ = 0;
            other.handle_ = kNoHandle;
        }
        return *this;
    }

    /// @return a bounds-checked view over the whole file.
    [[nodiscard]] Bytes bytes() const noexcept {
        return Bytes(static_cast<const std::byte*>(addr_), size_);
    }

    /// @return the mapped size in bytes.
    [[nodiscard]] std::uint64_t size() const noexcept { return size_; }

    /// @return whether this object currently owns a mapping.
    [[nodiscard]] bool is_open() const noexcept { return addr_ != nullptr; }

  private:
    /// How much of the file to ask the kernel to fault in ahead of use. The descriptors and the
    /// variable index sit near the front; reading the whole of a multi-gigabyte file here would
    /// defeat the point of mapping it lazily.
    static constexpr std::uint64_t kAdviseBytes = 1U << 20;

    /// Release whatever is held, leaving the object empty.
    void reset() noexcept {
        if (addr_ != nullptr) {
            Sys::unmap_file(addr_, size_);
            addr_ = nullptr;
        }
        if (handle_ != kNoHandle) {
            Sys::close_file(handle_);
            handle_ = kNoHandle;
        }
        size_ = 0;
    }

    void* addr_ = nullptr;
    std::uint64_t size_ = 0;
    handle_t handle_ = kNoHandle;
};

/// The mapping used in production.
using FileMapping = BasicFileMapping<PosixSysOps>;

}  // namespace cheatah::space::cdf::detail
