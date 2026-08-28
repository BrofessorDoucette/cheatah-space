#pragma once

/**
 * @file bytes.hpp
 * @brief space.cdf — the one place raw bytes are read, bounds-checked, and byte-swapped.
 *
 * Every multi-byte load in space.cdf goes through this file. That is a deliberate chokepoint
 * rather than a convenience, for three reasons that each independently justify it:
 *
 * **1. Alignment.** CDF's fields are not naturally aligned. `GDRoffset` is an 8-byte integer at
 * offset 12 of the CDR; `VXRnext` is 8 bytes at offset 12 of a VXR; a VXR's `Last[]` array starts
 * at `28 + 4*Nentries`, which is odd-aligned whenever `Nentries` is odd. And the mapping itself
 * offers no alignment guarantee for variable data. So every load here is `std::memcpy` into a
 * local followed by `__builtin_bswap`, and **never** a `reinterpret_cast` to a wider pointer.
 * That is not caution: a wide-pointer read of an unaligned address is undefined behaviour, and
 * UBSan is a hard QA-gate stage. Compilers fold memcpy+bswap into a single `movbe`/`ldr+rev`, so
 * the correct form is also the fast one.
 *
 * **2. Bounds.** space.cdf parses untrusted files. Confining reads to one type means there is
 * exactly one place a bounds check can be forgotten, and it is checked.
 *
 * **3. Coverage.** The `require()` below is the reason a parser with hundreds of validations can
 * hit the 100%-line bar. Each check is a single call whose line executes on every happy-path
 * parse; the `throw` lives in one function, covered once. Written the obvious way —
 * `if (bad) { throw CdfError{...}; }` at each site — every one of those throws would be a
 * separate uncovered line needing its own crafted malformed file just to satisfy the gate. That
 * is why parsers famously cannot reach 100%, and it is the trap this design steps around. The
 * errors are still each provoked individually by the unit tests, because SECURITY.md requires it;
 * coverage simply stops being the reason they must be.
 *
 * All control integers in a CDF are **big-endian regardless of the file's data encoding** — the
 * encoding applies to variable and attribute VALUES only. Readers that miss that decode the
 * record structure of a little-endian file into nonsense.
 *
 * Header-only and allocation-free. No platform headers: this reads a span someone else mapped.
 */

#include "cheatah.hpp"

#include <bit>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <stdexcept>
#include <string>

#include "types.hpp"

namespace cheatah::space::cdf {

/**
 * The exception every malformed-file path raises.
 *
 * Carries the byte offset as well as the code, because "bad record type" without a location is
 * not actionable on an 8 MB file.
 */
class CdfError : public std::runtime_error {
  public:
    /**
     * Construct from a code and the offset the trouble was found at.
     * @param code what went wrong.
     * @param offset byte offset into the file.
     */
    CdfError(ErrorCode code, std::uint64_t offset)
        : std::runtime_error(std::string(detail::error_message(code)) + " (at byte "
                             + std::to_string(offset) + ")"),
          code_(code), offset_(offset) {}

    /// @return the error code, for programmatic handling.
    [[nodiscard]] ErrorCode code() const noexcept { return code_; }
    /// @return the byte offset the error was detected at.
    [[nodiscard]] std::uint64_t offset() const noexcept { return offset_; }

  private:
    ErrorCode code_;
    std::uint64_t offset_;
};

namespace detail {

/// Raise a CdfError. The ONLY place space.cdf throws — see the file docs for why that matters.
/// @param code what went wrong. @param offset where. @throws CdfError always.
[[noreturn]] inline void throw_cdf(ErrorCode code, std::uint64_t offset) {
    throw CdfError(code, offset);
}

/// Assert a parser invariant. One line per check, and the line runs on every successful parse.
/// @param ok the invariant. @param code the failure to report. @param offset where.
inline void require(bool ok, ErrorCode code, std::uint64_t offset) {
    if (!ok) [[unlikely]] { throw_cdf(code, offset); }
}

/**
 * A bounds-checked view over the bytes of a CDF, with big-endian accessors.
 *
 * Copyable and cheap — it is a pointer and a length. It does not own the memory; the mapping
 * (or buffer) outlives every Bytes made from it.
 */
class Bytes {
  public:
    /// An empty view.
    Bytes() = default;

    /// View @p size bytes at @p data.
    /// @param data first byte. @param size how many bytes are readable.
    Bytes(const std::byte* data, std::uint64_t size) noexcept : data_(data), size_(size) {}

    /// @return how many bytes the view spans.
    [[nodiscard]] std::uint64_t size() const noexcept { return size_; }
    /// @return a pointer to the first byte.
    [[nodiscard]] const std::byte* data() const noexcept { return data_; }

    /// Whether [offset, offset+count) lies inside the view. Overflow-safe: the addition is done
    /// as a subtraction against the size so a hostile offset near 2^64 cannot wrap into range.
    /// @param offset start. @param count length.
    /// @return true when the range is readable.
    [[nodiscard]] bool contains(std::uint64_t offset, std::uint64_t count) const noexcept {
        return offset <= size_ && count <= size_ - offset;
    }

    /// Bounds-check a range, raising ErrorCode::TruncatedFile when it escapes the view.
    /// @param offset start. @param count length.
    void require_range(std::uint64_t offset, std::uint64_t count) const {
        require(contains(offset, count), ErrorCode::TruncatedFile, offset);
    }

    /// A sub-view, bounds-checked.
    /// @param offset start. @param count length.
    /// @return the sub-view.
    [[nodiscard]] Bytes subspan(std::uint64_t offset, std::uint64_t count) const {
        require_range(offset, count);
        return Bytes(data_ + offset, count);
    }

    /// One byte, bounds-checked.
    /// @param offset which byte.
    /// @return its value.
    [[nodiscard]] std::uint8_t u8(std::uint64_t offset) const {
        require_range(offset, 1);
        return static_cast<std::uint8_t>(data_[offset]);
    }

    /// A big-endian unsigned 16-bit integer.
    /// @param offset where the field starts.
    /// @return the decoded value.
    [[nodiscard]] std::uint16_t be_u16(std::uint64_t offset) const {
        return load_be<std::uint16_t>(offset);
    }

    /// A big-endian unsigned 32-bit integer.
    /// @param offset where the field starts.
    /// @return the decoded value.
    [[nodiscard]] std::uint32_t be_u32(std::uint64_t offset) const {
        return load_be<std::uint32_t>(offset);
    }

    /// A big-endian unsigned 64-bit integer.
    /// @param offset where the field starts.
    /// @return the decoded value.
    [[nodiscard]] std::uint64_t be_u64(std::uint64_t offset) const {
        return load_be<std::uint64_t>(offset);
    }

    /// A big-endian signed 32-bit integer — the width of every record type, data type and count
    /// field in the format.
    /// @param offset where the field starts.
    /// @return the decoded value.
    [[nodiscard]] std::int32_t be_i32(std::uint64_t offset) const {
        return static_cast<std::int32_t>(load_be<std::uint32_t>(offset));
    }

    /// A big-endian signed 64-bit integer — the width of every file offset and record size in
    /// CDF 3.x.
    /// @param offset where the field starts.
    /// @return the decoded value.
    [[nodiscard]] std::int64_t be_i64(std::uint64_t offset) const {
        return static_cast<std::int64_t>(load_be<std::uint64_t>(offset));
    }

    /**
     * A file offset that must land inside the file.
     *
     * Offsets are stored signed, and a negative one is meaningless; checking here rather than at
     * each of the ~30 call sites is the point of the chokepoint.
     *
     * @param offset where the field starts.
     * @return the offset it holds.
     */
    [[nodiscard]] std::uint64_t be_offset(std::uint64_t offset) const {
        const std::int64_t raw = be_i64(offset);
        require(raw >= 0, ErrorCode::BadOffset, offset);
        const auto value = static_cast<std::uint64_t>(raw);
        require(value <= size_, ErrorCode::BadOffset, offset);
        return value;
    }

    /**
     * A fixed-width, NUL-padded name field.
     *
     * CDF pads names into a fixed field (256 bytes in 3.x, 64 pre-2.5) and the spec calls the
     * bytes after the terminator undefined, so the trailing padding is stripped rather than
     * trusted. The result points into the mapping; it does not own its bytes.
     *
     * @param offset where the field starts.
     * @param width the field's full width in bytes.
     * @return the name, without padding.
     */
    [[nodiscard]] std::string_view name(std::uint64_t offset, std::uint64_t width) const {
        require_range(offset, width);
        const auto* first = reinterpret_cast<const char*>(data_ + offset);
        std::uint64_t length = 0;
        while (length < width && first[length] != '\0') { ++length; }
        return {first, static_cast<std::size_t>(length)};
    }

  private:
    /// Load a big-endian integer of type T.
    ///
    /// memcpy, then byteswap. Never a wide-pointer cast — CDF fields are unaligned (see the file
    /// docs). `std::endian` decides the swap at compile time, so exactly one form is emitted and
    /// the other leaves no uncovered branch behind.
    template <class T>
    [[nodiscard]] T load_be(std::uint64_t offset) const {
        require_range(offset, sizeof(T));
        T value{};
        std::memcpy(&value, data_ + offset, sizeof(T));
        if constexpr (std::endian::native == std::endian::little) {
            if constexpr (sizeof(T) == 2) { value = __builtin_bswap16(value); }
            else if constexpr (sizeof(T) == 4) { value = __builtin_bswap32(value); }
            else { value = __builtin_bswap64(value); }
        }
        return value;
    }

    const std::byte* data_ = nullptr;
    std::uint64_t size_ = 0;
};

}  // namespace detail

}  // namespace cheatah::space::cdf
