#pragma once

/**
 * @file encoding.hpp
 * @brief space.cdf — turning a run of stored values into a run of host values, in one pass.
 *
 * The decoder is chosen once per request, not once per value: `decode_run` switches on the
 * variable's storage type a single time and then runs a tight loop that loads, byte-swaps if the
 * file's encoding differs from the host's, and converts. The swap decision is a runtime bool
 * because it depends on the FILE, but the swap itself is selected by width at compile time.
 *
 * Seventeen data types collapse to nine storage types (an EPOCH is a double, a TT2000 is an
 * int64, a CHAR is an unsigned byte, an EPOCH16 is two doubles). Converting through the storage
 * type is what keeps this file small — and what makes the lossless rule checkable: a `double`
 * holds every storage type exactly except int64, and a `long long` holds every integer type but
 * no float, so those two refusals are the whole of the policy.
 *
 * Only the two IEEE decode classes are handled. VAX floats are refused at open, before any value
 * is read — see reader.hpp.
 */

#include "cheatah.hpp"

#include <cstdint>
#include <cstring>
#include <type_traits>

#include "bytes.hpp"
#include "types.hpp"

namespace cheatah::space::cdf::detail {

/// How a data type is physically laid out. Seventeen CDF types share nine layouts, and this is
/// the only thing the decoder switches on — which is what keeps decode_run() free of arms that
/// nothing can reach, and therefore free of lines no test can cover.
enum class StorageKind : std::uint8_t { I8, I16, I32, I64, U8, U16, U32, F32, F64, Unknown };

/// The layout a data type is stored in. EPOCH16 is two F64 per value; the caller accounts for
/// that by counting raw elements rather than values.
/// @param type a CDF data type. @return its layout, or Unknown for a value the format does not define.
constexpr StorageKind storage_kind(DataType type) noexcept {
    switch (type) {
        case DataType::Int1:
        case DataType::Byte: return StorageKind::I8;
        case DataType::Int2: return StorageKind::I16;
        case DataType::Int4: return StorageKind::I32;
        case DataType::Int8:
        case DataType::TimeTt2000: return StorageKind::I64;
        case DataType::Uint1:
        case DataType::Char:
        case DataType::Uchar: return StorageKind::U8;
        case DataType::Uint2: return StorageKind::U16;
        case DataType::Uint4: return StorageKind::U32;
        case DataType::Real4:
        case DataType::Float: return StorageKind::F32;
        case DataType::Real8:
        case DataType::Double:
        case DataType::Epoch:
        case DataType::Epoch16: return StorageKind::F64;
    }
    return StorageKind::Unknown;
}

/// Whether a layout holds a floating-point value — the half of the lossless rule that decides
/// which of values()/values_i64() may read a variable.
/// @param kind a layout. @return true for F32 and F64.
constexpr bool is_float_kind(StorageKind kind) noexcept {
    return kind == StorageKind::F32 || kind == StorageKind::F64;
}

/// Bytes in one raw element of a layout.
/// @param kind a layout. @return its size, or 0 for Unknown.
constexpr std::uint64_t storage_bytes(StorageKind kind) noexcept {
    switch (kind) {
        case StorageKind::I8:
        case StorageKind::U8: return 1;
        case StorageKind::I16:
        case StorageKind::U16: return 2;
        case StorageKind::I32:
        case StorageKind::U32:
        case StorageKind::F32: return 4;
        case StorageKind::I64:
        case StorageKind::F64: return 8;
        case StorageKind::Unknown: return 0;
    }
    return 0;
}

/// Bytes per raw element of a data type: its element size, except for EPOCH16, whose two
/// doubles are decoded as two elements.
/// @param type a CDF data type. @return bytes per raw element.
constexpr std::uint64_t raw_element_bytes(DataType type) noexcept {
    return storage_bytes(storage_kind(type));
}

/// Whether values in this decode class must be byte-swapped on this host.
/// @param cls the file's decode class. @return true when a swap is needed.
constexpr bool needs_swap(EncodingClass cls) noexcept { return !is_host_native(cls); }

/// Load one raw value, swapping bytes if asked. memcpy in and out — never a wide-pointer read
/// of possibly unaligned file bytes.
/// @tparam Raw the storage type. @param p the bytes. @param swap whether to reverse them.
/// @return the value.
template <class Raw>
inline Raw load_raw(const std::byte* p, bool swap) noexcept {
    if constexpr (sizeof(Raw) == 1) {
        Raw v{};
        std::memcpy(&v, p, 1);
        return v;
    } else {
        using U = std::conditional_t<sizeof(Raw) == 2, std::uint16_t,
                  std::conditional_t<sizeof(Raw) == 4, std::uint32_t, std::uint64_t>>;
        U u{};
        std::memcpy(&u, p, sizeof(U));
        if (swap) {
            if constexpr (sizeof(U) == 2) { u = __builtin_bswap16(u); }
            else if constexpr (sizeof(U) == 4) { u = __builtin_bswap32(u); }
            else { u = __builtin_bswap64(u); }
        }
        Raw v{};
        std::memcpy(&v, &u, sizeof(Raw));
        return v;
    }
}

/// Convert @p n raw elements into @p dst.
/// @tparam Out the host type. @tparam Raw the storage type.
/// @param src the stored bytes. @param dst where to write. @param n how many. @param swap byte order.
template <class Out, class Raw>
inline void convert_run(const std::byte* src, Out* dst, std::uint64_t n, bool swap) noexcept {
    for (std::uint64_t i = 0; i < n; ++i) {
        dst[i] = static_cast<Out>(load_raw<Raw>(src + i * sizeof(Raw), swap));
    }
}

/**
 * Decode @p n raw elements of @p type into host values.
 *
 * Refuses the two lossy combinations with ErrorCode::LossyConversion rather than rounding
 * silently: int64 storage into a double, and any float storage into an integer.
 *
 * @tparam Out double or long long.
 * @param src the stored bytes, at least n * raw_element_bytes(type) of them.
 * @param type the variable's data type.
 * @param swap whether the file's byte order differs from the host's.
 * @param n how many raw elements.
 * @param dst where to write n values.
 * @param offset the file offset, for the error.
 */
template <class Out>
inline void decode_run(const std::byte* src, DataType type, bool swap, std::uint64_t n, Out* dst,
                       std::uint64_t offset) {
    constexpr bool kToFloat = std::is_floating_point_v<Out>;
    // The lossless rule is enforced inside the arms rather than ahead of the switch: a guard
    // before it would make the very arms it guards unreachable, and an arm no instantiation can
    // enter is a line no test can cover. A double holds every layout exactly except int64; an
    // int64 holds every integer layout and no float.
    switch (storage_kind(type)) {
        case StorageKind::I8: convert_run<Out, std::int8_t>(src, dst, n, swap); return;
        case StorageKind::I16: convert_run<Out, std::int16_t>(src, dst, n, swap); return;
        case StorageKind::I32: convert_run<Out, std::int32_t>(src, dst, n, swap); return;
        case StorageKind::U8: convert_run<Out, std::uint8_t>(src, dst, n, swap); return;
        case StorageKind::U16: convert_run<Out, std::uint16_t>(src, dst, n, swap); return;
        case StorageKind::U32: convert_run<Out, std::uint32_t>(src, dst, n, swap); return;
        case StorageKind::I64:
            if constexpr (kToFloat) { throw_cdf(ErrorCode::LossyConversion, offset); }
            else { convert_run<Out, std::int64_t>(src, dst, n, swap); return; }
        case StorageKind::F32:
            if constexpr (kToFloat) { convert_run<Out, float>(src, dst, n, swap); return; }
            else { throw_cdf(ErrorCode::LossyConversion, offset); }
        case StorageKind::F64:
            if constexpr (kToFloat) { convert_run<Out, double>(src, dst, n, swap); return; }
            else { throw_cdf(ErrorCode::LossyConversion, offset); }
        case StorageKind::Unknown: break;
    }
    throw_cdf(ErrorCode::BadDataType, offset);
}

}  // namespace cheatah::space::cdf::detail
