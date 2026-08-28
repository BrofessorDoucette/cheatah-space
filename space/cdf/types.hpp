#pragma once

/**
 * @file types.hpp
 * @brief space.cdf — the format's vocabulary: data types, encodings, record types, errors.
 *
 * The closed enumerations the CDF Internal Format Description defines, plus the small total
 * functions over them. No I/O and no parsing happen here; this is the layer everything else
 * agrees on. `import space.cdf` resolves it through `cdf.hpp`.
 *
 * Every `switch` over an enumeration in this file deliberately has **no `default:`**. Clang then
 * warns when an enumerator is added without a case, so the tables cannot silently rot — and
 * there is no unreachable default line for the coverage gate to complain about.
 *
 * A note on the numbering, because it looks arbitrary and is not: the values are the ones
 * written into CDF files, so they are fixed by the format and not by us. CDF_REAL4 is 21 and
 * CDF_FLOAT is 44 even though both are a 4-byte IEEE single; CDF_REAL8 (22) and CDF_DOUBLE (45)
 * are likewise the same layout. The pairs exist for historical reasons and both appear in real
 * files, so both are handled everywhere rather than normalized away at the boundary.
 *
 * Cross-platform, header-only, allocation-free: no platform headers, no I/O, no global state.
 */

#include "cheatah.hpp"

#include <cstddef>
#include <cstdint>
#include <string_view>

namespace cheatah::space::cdf {

// ---- file identification -------------------------------------------------------------------

/// First magic number of a CDF 3.x file, at byte offset 0.
inline constexpr std::uint32_t kMagicV3 = 0xCDF30001U;
/// First magic number of a CDF 2.6/2.7 file — the 32-bit-offset generation.
inline constexpr std::uint32_t kMagicV26 = 0xCDF26002U;
/// First magic number of a pre-2.6 file. Not supported; see ErrorCode::UnsupportedPreV26.
inline constexpr std::uint32_t kMagicPreV26 = 0x0000FFFFU;
/// Second magic number: the file's records are stored uncompressed.
inline constexpr std::uint32_t kMagicUncompressed = 0x0000FFFFU;
/// Second magic number: the whole file is compressed and a CCR follows the magic.
inline constexpr std::uint32_t kMagicCompressed = 0xCCCC0001U;

/// Byte offset of the CDR, which always follows the two 4-byte magic numbers.
inline constexpr std::uint64_t kCdrOffset = 8;
/// Bytes a checksummed file carries past `GDR.eof` — an MD5 digest, not a CDF record.
inline constexpr std::uint64_t kChecksumBytes = 16;

// ---- record types ---------------------------------------------------------------------------

/**
 * The internal record types, as written in the 4-byte big-endian field at record offset 8.
 *
 * @note UnusedInternal is -1 and marks a hole left by a rewrite; such records are skipped, never
 *       parsed. Every other value is a real record with a defined layout.
 */
enum class RecordType : std::int32_t {
    Cdr = 1,                 ///< CDF Descriptor Record — always at offset 8.
    Gdr = 2,                 ///< Global Descriptor Record — the entry point to everything else.
    RVariableDescriptor = 3, ///< rVDR: an rVariable's descriptor.
    Adr = 4,                 ///< Attribute Descriptor Record.
    AgrEdr = 5,              ///< Attribute g/rEntry Descriptor Record.
    Vxr = 6,                 ///< Variable Index Record — may index VVRs, CVVRs, or more VXRs.
    Vvr = 7,                 ///< Variable Values Record — uncompressed record data.
    ZVariableDescriptor = 8, ///< zVDR: a zVariable's descriptor.
    AzEdr = 9,               ///< Attribute zEntry Descriptor Record.
    Ccr = 10,                ///< Compressed CDF Record — the whole file, compressed.
    Cpr = 11,                ///< Compression Parameters Record.
    Spr = 12,                ///< Sparseness Parameters Record.
    Cvvr = 13,               ///< Compressed Variable Values Record.
    UnusedInternal = -1,     ///< A hole in the file. Skipped, never parsed.
};

namespace detail {
// Total functions over the format's closed enumerations. They live in `detail` because
// nothing user-facing names them yet: the purr-visible way to ask a type question is to
// ask a VARIABLE (`cdf.type_name(f, "Epoch")`), which arrives with the reader. They get
// promoted, with worked examples, when that lands.

/**
 * The name of a record type, for diagnostics.
 * @param type the record type read from a file.
 * @return the spec's abbreviation, or "unknown" for a value the format does not define.
 * @complexity O(1).
 * @alloc none.
 */
constexpr std::string_view record_type_name(RecordType type) noexcept {
    switch (type) {
        case RecordType::Cdr: return "CDR";
        case RecordType::Gdr: return "GDR";
        case RecordType::RVariableDescriptor: return "rVDR";
        case RecordType::Adr: return "ADR";
        case RecordType::AgrEdr: return "AgrEDR";
        case RecordType::Vxr: return "VXR";
        case RecordType::Vvr: return "VVR";
        case RecordType::ZVariableDescriptor: return "zVDR";
        case RecordType::AzEdr: return "AzEDR";
        case RecordType::Ccr: return "CCR";
        case RecordType::Cpr: return "CPR";
        case RecordType::Spr: return "SPR";
        case RecordType::Cvvr: return "CVVR";
        case RecordType::UnusedInternal: return "UIR";
    }
    return "unknown";
}

}  // namespace detail

// ---- data types -----------------------------------------------------------------------------

/**
 * The CDF data types. The numbering is the format's, not ours — see the file docs for why
 * REAL4/FLOAT and REAL8/DOUBLE are distinct values with identical layouts.
 */
enum class DataType : std::int32_t {
    Int1 = 1,      ///< 1-byte signed.
    Int2 = 2,      ///< 2-byte signed.
    Int4 = 4,      ///< 4-byte signed.
    Int8 = 8,      ///< 8-byte signed. Not representable in a double without loss.
    Uint1 = 11,    ///< 1-byte unsigned.
    Uint2 = 12,    ///< 2-byte unsigned.
    Uint4 = 14,    ///< 4-byte unsigned.
    Real4 = 21,    ///< 4-byte IEEE single.
    Real8 = 22,    ///< 8-byte IEEE double.
    Epoch = 31,    ///< CDF_EPOCH: a double, milliseconds since year 0.
    Epoch16 = 32,  ///< CDF_EPOCH16: two doubles — seconds since year 0, then picoseconds.
    TimeTt2000 = 33,  ///< CDF_TIME_TT2000: int64 nanoseconds since J2000, with leap seconds.
    Byte = 41,     ///< 1-byte signed. Distinct from Int1 only by name.
    Float = 44,    ///< 4-byte IEEE single. Same layout as Real4.
    Double = 45,   ///< 8-byte IEEE double. Same layout as Real8.
    Char = 51,     ///< Character data. NEVER NUL-terminated in file; NumElems is the length.
    Uchar = 52,    ///< Unsigned character data. Same storage as Char.
};

namespace detail {

/**
 * The size in bytes of one element of a data type.
 * @param type a CDF data type.
 * @return the element size, or 0 when @p type is not a value the format defines.
 * @complexity O(1).
 * @alloc none.
 */
constexpr long long element_size(DataType type) noexcept {
    switch (type) {
        case DataType::Int1:
        case DataType::Uint1:
        case DataType::Byte:
        case DataType::Char:
        case DataType::Uchar: return 1;
        case DataType::Int2:
        case DataType::Uint2: return 2;
        case DataType::Int4:
        case DataType::Uint4:
        case DataType::Real4:
        case DataType::Float: return 4;
        case DataType::Int8:
        case DataType::Real8:
        case DataType::Epoch:
        case DataType::Double:
        case DataType::TimeTt2000: return 8;
        case DataType::Epoch16: return 16;
    }
    return 0;
}

/**
 * Whether a data type is one of the format's defined values.
 *
 * Worth having as its own predicate: a corrupt or hostile file can put any 32-bit integer in a
 * data-type field, and element_size() returning 0 must not be confused with a real zero-size type.
 *
 * @param type the value read from a file.
 * @return true when @p type is a data type the format defines.
 * @complexity O(1).
 * @alloc none.
 */
constexpr bool is_known_data_type(DataType type) noexcept { return element_size(type) != 0; }

/**
 * Whether a data type holds character data rather than numbers.
 *
 * The distinction is load-bearing: for CHAR and UCHAR the variable's `NumElems` is the STRING
 * LENGTH, not a count of numeric elements, and the bytes are not NUL-terminated in the file.
 *
 * @param type a CDF data type.
 * @return true for Char and Uchar.
 * @complexity O(1).
 * @alloc none.
 */
constexpr bool is_char_type(DataType type) noexcept {
    return type == DataType::Char || type == DataType::Uchar;
}

/**
 * Whether a data type is one of the three CDF time types.
 * @param type a CDF data type.
 * @return true for Epoch, Epoch16 and TimeTt2000.
 * @complexity O(1).
 * @alloc none.
 */
constexpr bool is_time_type(DataType type) noexcept {
    return type == DataType::Epoch || type == DataType::Epoch16
           || type == DataType::TimeTt2000;
}

/**
 * The spec's name for a data type, for diagnostics and for listing a file's contents.
 * @param type a CDF data type.
 * @return the `CDF_*` name, or "unknown" for a value the format does not define.
 * @complexity O(1).
 * @alloc none.
 */
constexpr std::string_view data_type_name(DataType type) noexcept {
    switch (type) {
        case DataType::Int1: return "CDF_INT1";
        case DataType::Int2: return "CDF_INT2";
        case DataType::Int4: return "CDF_INT4";
        case DataType::Int8: return "CDF_INT8";
        case DataType::Uint1: return "CDF_UINT1";
        case DataType::Uint2: return "CDF_UINT2";
        case DataType::Uint4: return "CDF_UINT4";
        case DataType::Real4: return "CDF_REAL4";
        case DataType::Real8: return "CDF_REAL8";
        case DataType::Epoch: return "CDF_EPOCH";
        case DataType::Epoch16: return "CDF_EPOCH16";
        case DataType::TimeTt2000: return "CDF_TIME_TT2000";
        case DataType::Byte: return "CDF_BYTE";
        case DataType::Float: return "CDF_FLOAT";
        case DataType::Double: return "CDF_DOUBLE";
        case DataType::Char: return "CDF_CHAR";
        case DataType::Uchar: return "CDF_UCHAR";
    }
    return "unknown";
}

}  // namespace detail

// ---- encodings ------------------------------------------------------------------------------

/**
 * The byte orders and floating-point formats a CDF's values may be stored in.
 *
 * @note "NETWORK" is not XDR despite the name — it is plain big-endian IEEE-754. That is worth
 *       stating because assuming otherwise invents a whole conversion layer that is not needed.
 */
enum class Encoding : std::int32_t {
    Network = 1,      ///< Big-endian IEEE. The most common encoding in the public archive.
    Sun = 2,          ///< Big-endian IEEE.
    Vax = 3,          ///< Little-endian, VAX F/D_FLOAT — not IEEE.
    DecStation = 4,   ///< Little-endian IEEE.
    Sgi = 5,          ///< Big-endian IEEE.
    IbmPc = 6,        ///< Little-endian IEEE. The other common archive encoding.
    IbmRs = 7,        ///< Big-endian IEEE.
    Host = 8,         ///< Write-time selector meaning "this machine"; resolved when reading.
    Ppc = 9,          ///< Big-endian IEEE.
    Mac = 10,         ///< Big-endian IEEE.
    Hp = 11,          ///< Big-endian IEEE.
    NeXT = 12,        ///< Big-endian IEEE.
    AlphaOsf1 = 13,   ///< Little-endian IEEE.
    AlphaVmsD = 14,   ///< Little-endian, VAX D_FLOAT.
    AlphaVmsG = 15,   ///< Little-endian, VAX G_FLOAT.
    AlphaVmsI = 16,   ///< Little-endian IEEE.
    ArmLittle = 17,   ///< Little-endian IEEE.
    ArmBig = 18,      ///< Big-endian IEEE.
    Ia64VmsI = 19,    ///< Little-endian IEEE.
    Ia64VmsD = 20,    ///< Little-endian, VAX D_FLOAT.
    Ia64VmsG = 21,    ///< Little-endian, VAX G_FLOAT.
};

/**
 * The four ways a CDF's bytes can actually need decoding.
 *
 * Twenty-one encodings collapse to four decode behaviours, which is the single simplification
 * that keeps the decoder small: what matters is integer byte order and whether floats are IEEE
 * or one of the two VAX formats.
 */
enum class EncodingClass : std::uint8_t {
    BigIeee,     ///< Big-endian integers, IEEE-754 floats.
    LittleIeee,  ///< Little-endian integers, IEEE-754 floats.
    LittleVaxD,  ///< Little-endian integers, VAX F_FLOAT / D_FLOAT.
    LittleVaxG,  ///< Little-endian integers, VAX F_FLOAT / G_FLOAT.
};

/// True when this build targets a little-endian host. CDF stores no encoding meaning "native",
/// so the Host selector has to be resolved against the machine actually doing the reading.
inline constexpr bool kHostIsLittleEndian = (__BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__);

namespace detail {

/**
 * Which decode behaviour an encoding requires.
 *
 * Encoding::Host is a write-time selector, so it resolves here against the host's own byte order
 * rather than naming a fixed class.
 *
 * @param encoding the encoding read from a file's CDR.
 * @return the decode class to use.
 * @complexity O(1).
 * @alloc none.
 */
constexpr EncodingClass encoding_class_of(Encoding encoding) noexcept {
    switch (encoding) {
        case Encoding::Network:
        case Encoding::Sun:
        case Encoding::Sgi:
        case Encoding::IbmRs:
        case Encoding::Ppc:
        case Encoding::Mac:
        case Encoding::Hp:
        case Encoding::NeXT:
        case Encoding::ArmBig: return EncodingClass::BigIeee;
        case Encoding::DecStation:
        case Encoding::IbmPc:
        case Encoding::AlphaOsf1:
        case Encoding::AlphaVmsI:
        case Encoding::ArmLittle:
        case Encoding::Ia64VmsI: return EncodingClass::LittleIeee;
        case Encoding::Vax:
        case Encoding::AlphaVmsD:
        case Encoding::Ia64VmsD: return EncodingClass::LittleVaxD;
        case Encoding::AlphaVmsG:
        case Encoding::Ia64VmsG: return EncodingClass::LittleVaxG;
        case Encoding::Host:
            return kHostIsLittleEndian ? EncodingClass::LittleIeee : EncodingClass::BigIeee;
    }
    return EncodingClass::BigIeee;
}

/**
 * Whether an encoding value is one the format defines. Every code 1..21 is assigned, so this is
 * a range check — but a corrupt CDR can hold anything, and encoding_class_of() must not be the
 * thing that discovers it.
 * @param encoding the value read from a CDR.
 * @return true when the format defines it.
 * @complexity O(1).
 * @alloc none.
 */
constexpr bool is_known_encoding(Encoding encoding) noexcept {
    const auto v = static_cast<std::int32_t>(encoding);
    return v >= static_cast<std::int32_t>(Encoding::Network)
           && v <= static_cast<std::int32_t>(Encoding::Ia64VmsG);
}

/**
 * The name of a decode class, for diagnostics.
 * @param cls a decode class.
 * @return a short lowercase name.
 * @complexity O(1).
 * @alloc none.
 */
constexpr std::string_view encoding_class_name(EncodingClass cls) noexcept {
    switch (cls) {
        case EncodingClass::BigIeee: return "big-ieee";
        case EncodingClass::LittleIeee: return "little-ieee";
        case EncodingClass::LittleVaxD: return "little-vax-d";
        case EncodingClass::LittleVaxG: return "little-vax-g";
    }
    return "unknown";
}

/**
 * Whether values in this class can be read straight from the file with no conversion at all.
 *
 * True only when the stored integer byte order matches the host's and the floats are IEEE — the
 * case where a reader may hand back a span into the mapping instead of decoding.
 *
 * @param cls a decode class.
 * @return true when no byte swapping or float conversion is needed on this host.
 * @complexity O(1).
 * @alloc none.
 */
constexpr bool is_host_native(EncodingClass cls) noexcept {
    return cls == (kHostIsLittleEndian ? EncodingClass::LittleIeee : EncodingClass::BigIeee);
}

}  // namespace detail

// ---- majority, sparseness, compression --------------------------------------------------------

/// How a multi-dimensional variable's values are laid out within one record.
enum class Majority : std::uint8_t {
    Row,     ///< Last dimension varies fastest — C order.
    Column,  ///< First dimension varies fastest — Fortran order.
};

/// What a variable's unwritten records contain.
enum class SparseRecords : std::int32_t {
    None = 0,      ///< Every record 0..MaxRec is physically present.
    Pad = 1,       ///< Missing records read as the variable's pad value.
    Previous = 2,  ///< Missing records repeat the last written record.
};

/// The compression algorithms CDF defines. Only Gzip occurs in the public archive.
enum class CompressionType : std::int32_t {
    None = 0,   ///< Stored uncompressed.
    Rle = 1,    ///< Run-length encoding of zeros.
    Huff = 2,   ///< Static Huffman.
    Ahuff = 3,  ///< Adaptive Huffman.
    Gzip = 5,   ///< DEFLATE. Note 4 is unassigned by the format.
};

namespace detail {

/**
 * The name of a compression type, for diagnostics.
 * @param type a compression type.
 * @return a short lowercase name, or "unknown" for an undefined value.
 * @complexity O(1).
 * @alloc none.
 */
constexpr std::string_view compression_name(CompressionType type) noexcept {
    switch (type) {
        case CompressionType::None: return "none";
        case CompressionType::Rle: return "rle";
        case CompressionType::Huff: return "huff";
        case CompressionType::Ahuff: return "ahuff";
        case CompressionType::Gzip: return "gzip";
    }
    return "unknown";
}

}  // namespace detail

// ---- CDR flag bits ------------------------------------------------------------------------------

/// CDR flags bit 0: set when the file's majority is row-major.
inline constexpr std::uint32_t kCdrFlagRowMajority = 1U << 0;
/// CDR flags bit 1: set when all data lives in this one file rather than `.v*`/`.z*` companions.
inline constexpr std::uint32_t kCdrFlagSingleFile = 1U << 1;
/// CDR flags bit 2: set when the file carries a trailing checksum.
inline constexpr std::uint32_t kCdrFlagChecksum = 1U << 2;
/// CDR flags bit 3: set when that checksum is MD5. Meaningful only with kCdrFlagChecksum.
inline constexpr std::uint32_t kCdrFlagChecksumMd5 = 1U << 3;

// ---- errors ---------------------------------------------------------------------------------

/**
 * Everything that can be wrong with a CDF we are asked to read.
 *
 * One flat enumeration rather than exception subtypes: a caller almost always wants to report
 * what happened and where, not to branch on the kind, and a closed set can be table-tested for
 * reachability — every code below is provoked by a crafted file in the unit tests, so none of
 * them is aspirational.
 */
enum class ErrorCode : std::uint8_t {
    None = 0,             ///< No error.
    NotCdf,               ///< The first magic number matches no CDF generation.
    UnsupportedPreV26,    ///< A pre-2.6 file: documented-undefined fields make it unparseable.
    TruncatedFile,        ///< A read ran past the end of the mapping.
    BadMagic,             ///< The second magic number is neither uncompressed nor compressed.
    BadRecordType,        ///< A record's type field is not one the format defines.
    UnexpectedRecordType, ///< A record is well-formed but not the type the layout requires here.
    BadRecordSize,        ///< A record's size is smaller than its own header, or runs past EOF.
    BadOffset,            ///< An internal offset points outside the file.
    BadDataType,          ///< A variable or attribute entry names an undefined data type.
    BadEncoding,          ///< The CDR names an encoding the format does not define.
    BadDimensions,        ///< A dimension count or size is negative, or the product overflows.
    VxrTreeTooDeep,       ///< The variable index nests deeper than the cap — likely a cycle.
    VxrTreeTooLarge,      ///< The variable index has more nodes than the cap — likely a cycle.
    RecordOutOfRange,     ///< A requested record number does not exist in the variable.
    VariableNotFound,     ///< No variable in the file carries the requested name.
    DuplicateVariable,    ///< Two variables share a name, so a lookup by name is ambiguous.
    UnsupportedMultiFile, ///< Data lives in companion files; not supported yet.
    UnsupportedSparse,    ///< The variable has sparse records; not supported yet.
    UnsupportedCompression, ///< The variable uses a compression this build cannot decode yet.
    DecompressionFailed,  ///< The compressed stream is malformed.
    DecompressedSizeMismatch, ///< Decompression produced a different size than the record claims.
    CannotOpen,           ///< The file could not be opened.
    CannotMap,            ///< The file could not be memory-mapped.
    EmptyFile,            ///< The file is too small to contain even the magic numbers.
    UnsupportedLayout,    ///< An N-D variable in a column-major file; the transpose is not in yet.
    UnsupportedEncoding,  ///< A VAX float encoding; the conversion is not in yet.
};

/// How many ErrorCode enumerators exist, including None. The unit tests assert one crafted
/// failure per code against this, so a new code cannot be added without a test that provokes it.
inline constexpr std::size_t kErrorCodeCount = 27;

namespace detail {

/**
 * A human-readable explanation of an error code.
 * @param code the error.
 * @return a one-line description, or "unknown error" for a value outside the enumeration.
 * @complexity O(1).
 * @alloc none.
 */
constexpr std::string_view error_message(ErrorCode code) noexcept {
    switch (code) {
        case ErrorCode::None: return "no error";
        case ErrorCode::NotCdf: return "not a CDF file (bad magic number)";
        case ErrorCode::UnsupportedPreV26: return "CDF 2.x files are not supported yet";
        case ErrorCode::TruncatedFile: return "file ends in the middle of a record";
        case ErrorCode::BadMagic: return "second magic number is neither uncompressed nor compressed";
        case ErrorCode::BadRecordType: return "undefined record type";
        case ErrorCode::UnexpectedRecordType: return "record is not the type expected here";
        case ErrorCode::BadRecordSize: return "record size is impossible";
        case ErrorCode::BadOffset: return "internal offset points outside the file";
        case ErrorCode::BadDataType: return "undefined data type";
        case ErrorCode::BadEncoding: return "undefined encoding";
        case ErrorCode::BadDimensions: return "impossible dimension count or size";
        case ErrorCode::VxrTreeTooDeep: return "variable index nests too deeply (cycle?)";
        case ErrorCode::VxrTreeTooLarge: return "variable index has too many nodes (cycle?)";
        case ErrorCode::RecordOutOfRange: return "record number does not exist in this variable";
        case ErrorCode::VariableNotFound: return "no variable with that name";
        case ErrorCode::DuplicateVariable: return "two variables share that name";
        case ErrorCode::UnsupportedMultiFile: return "multi-file CDFs are not supported yet";
        case ErrorCode::UnsupportedSparse: return "sparse records are not supported yet";
        case ErrorCode::UnsupportedCompression: return "this compression is not supported yet";
        case ErrorCode::DecompressionFailed: return "compressed data is malformed";
        case ErrorCode::DecompressedSizeMismatch: return "decompressed size does not match the record";
        case ErrorCode::CannotOpen: return "cannot open the file";
        case ErrorCode::CannotMap: return "cannot memory-map the file";
        case ErrorCode::EmptyFile: return "file is too small to be a CDF";
        case ErrorCode::UnsupportedLayout: return "multi-dimensional column-major variables are not supported yet";
        case ErrorCode::UnsupportedEncoding: return "VAX floating-point encodings are not supported yet";
    }
    return "unknown error";
}

}  // namespace detail

}  // namespace cheatah::space::cdf
