#pragma once

/**
 * @file records.hpp
 * @brief space.cdf — the internal records of a CDF 3.x file, parsed into plain aggregates.
 *
 * One `parse_*` per record type the reader needs: the file header, CDR, GDR, VDR (r and z),
 * VXR, VVR, CVVR and CPR. Each returns a small struct by value; names are `string_view`s into
 * the mapping, and nothing is copied at parse time. Everything is validated with one
 * `require()` line per invariant — see bytes.hpp for why that shape is load-bearing.
 *
 * The offsets below are the CDF 3.x Internal Format Description's, and every one of them was
 * checked against real files before being trusted here. Two are worth calling out because a
 * first reading of the spec gets them wrong: a VDR's `NumElems` is at +64 (not +56), and
 * `CPRorSPRoffset` is **-1**, not 0, when there is no compression record — so it must be read
 * as a signed value and only validated as an offset when the compressed flag is set.
 *
 * Records deliberately NOT parsed here, because nothing in this tranche reads them: ADR/AEDR
 * (attributes), SPR (sparse arrays, which the spec says are unimplemented), CCR (whole-file
 * compression), UIR (holes). The reader refuses what it cannot handle with a typed error rather
 * than parsing a record it will then ignore.
 */

#include "cheatah.hpp"

#include <array>
#include <cstddef>
#include <cstdint>
#include <string_view>

#include "bytes.hpp"
#include "types.hpp"

namespace cheatah::space::cdf::detail {

/// Every internal record begins with an 8-byte size and a 4-byte type.
inline constexpr std::uint64_t kRecordHeaderBytes = 12;
/// The most dimensions a variable may have. The format's own limit.
inline constexpr std::int32_t kMaxDims = 10;
/// Width of a name field in CDF 3.x.
inline constexpr std::uint64_t kNameBytes = 256;
/// Byte offset of a VDR's name field.
inline constexpr std::uint64_t kVdrNameOffset = 84;
/// Where a zVDR's dimension count sits: right after the 256-byte name.
inline constexpr std::uint64_t kVdrDimsOffset = kVdrNameOffset + kNameBytes;
/// A sane cap on VXR entries. The library writes 7; anything wildly larger is corruption.
inline constexpr std::int32_t kMaxVxrEntries = 65536;
/// A sane cap on CPR parameters. The format uses at most 1.
inline constexpr std::int32_t kMaxCprParams = 8;

// ---- headers ---------------------------------------------------------------------------------

/// The 12 bytes every record starts with, plus where it was read from.
struct RecordHeader {
    std::uint64_t offset{};  ///< Where the record starts in the file.
    std::uint64_t size{};    ///< Its total size in bytes, header included.
    RecordType type{};       ///< What it is.
};

/// Read and validate a record header: the size must cover the header and stay inside the file,
/// and the type must be one the format defines.
/// @param b the file. @param offset where the record starts.
/// @return the header.
inline RecordHeader read_record_header(const Bytes& b, std::uint64_t offset) {
    b.require_range(offset, kRecordHeaderBytes);
    const std::int64_t raw_size = b.be_i64(offset);
    require(raw_size >= static_cast<std::int64_t>(kRecordHeaderBytes), ErrorCode::BadRecordSize, offset);
    const auto size = static_cast<std::uint64_t>(raw_size);
    require(b.contains(offset, size), ErrorCode::BadRecordSize, offset);
    const auto type = static_cast<RecordType>(b.be_i32(offset + 8));
    require(record_type_name(type) != "unknown", ErrorCode::BadRecordType, offset + 8);
    return RecordHeader{offset, size, type};
}

/// Read a record header and insist on a particular type.
/// @param b the file. @param offset where the record starts. @param want the type required.
/// @return the header.
inline RecordHeader expect_record(const Bytes& b, std::uint64_t offset, RecordType want) {
    const RecordHeader h = read_record_header(b, offset);
    require(h.type == want, ErrorCode::UnexpectedRecordType, offset + 8);
    return h;
}

// ---- the file header -------------------------------------------------------------------------

/// The two magic numbers.
struct FileHeader {
    std::uint32_t magic1{};    ///< Identifies the format generation.
    std::uint32_t magic2{};    ///< Uncompressed, or whole-file compressed.
    bool compressed{};         ///< True when a CCR follows and the real file is inside it.
};

/// Validate the magic numbers. Generations other than 3.x are refused here, before anything
/// else is read, with the error that names why.
/// @param b the file.
/// @return the header.
inline FileHeader read_file_header(const Bytes& b) {
    require(b.size() >= kCdrOffset, ErrorCode::EmptyFile, 0);
    const std::uint32_t m1 = b.be_u32(0);
    const std::uint32_t m2 = b.be_u32(4);
    require(m1 != kMagicV26 && m1 != kMagicPreV26, ErrorCode::UnsupportedPreV26, 0);
    require(m1 == kMagicV3, ErrorCode::NotCdf, 0);
    require(m2 == kMagicUncompressed || m2 == kMagicCompressed, ErrorCode::BadMagic, 4);
    return FileHeader{m1, m2, m2 == kMagicCompressed};
}

// ---- CDR -------------------------------------------------------------------------------------

/// The CDF Descriptor Record: format version, encoding, and the layout flags.
struct Cdr {
    std::uint64_t gdr_offset{};  ///< Where the GDR is.
    std::int32_t version{};      ///< Major version (3).
    std::int32_t release{};      ///< Minor version.
    std::int32_t increment{};    ///< Patch level.
    Encoding encoding{};         ///< How VALUES are stored. Control ints are always big-endian.
    std::uint32_t flags{};       ///< The kCdrFlag* bits.

    /// @return whether multi-dimensional records are laid out last-dimension-fastest.
    [[nodiscard]] bool row_major() const noexcept { return (flags & kCdrFlagRowMajority) != 0; }
    /// @return whether all data lives in this one file.
    [[nodiscard]] bool single_file() const noexcept { return (flags & kCdrFlagSingleFile) != 0; }
    /// @return whether 16 checksum bytes follow `GDR.eof`.
    [[nodiscard]] bool has_checksum() const noexcept { return (flags & kCdrFlagChecksum) != 0; }
};

/// Bytes the CDR's fixed fields occupy before the 256-byte copyright text.
inline constexpr std::uint64_t kCdrFixedBytes = 56;

/// Parse the CDR, which is always at kCdrOffset.
/// @param b the file.
/// @return the record.
inline Cdr parse_cdr(const Bytes& b) {
    const RecordHeader h = expect_record(b, kCdrOffset, RecordType::Cdr);
    require(h.size >= kCdrFixedBytes, ErrorCode::BadRecordSize, h.offset);
    Cdr c;
    c.gdr_offset = b.be_offset(h.offset + 12);
    c.version = b.be_i32(h.offset + 20);
    c.release = b.be_i32(h.offset + 24);
    c.encoding = static_cast<Encoding>(b.be_i32(h.offset + 28));
    require(is_known_encoding(c.encoding), ErrorCode::BadEncoding, h.offset + 28);
    c.flags = b.be_u32(h.offset + 32);
    c.increment = b.be_i32(h.offset + 44);
    return c;
}

// ---- GDR -------------------------------------------------------------------------------------

/// The Global Descriptor Record: the heads of every chain, and the rVariable dimensions.
struct Gdr {
    std::uint64_t offset{};       ///< Where it was read from.
    std::uint64_t rvdr_head{};    ///< First rVDR, or 0.
    std::uint64_t zvdr_head{};    ///< First zVDR, or 0.
    std::uint64_t adr_head{};     ///< First ADR, or 0.
    std::uint64_t eof{};          ///< End of the last record. Excludes any checksum.
    std::int32_t nr_vars{};       ///< Number of rVariables.
    std::int32_t num_attr{};      ///< Number of attributes.
    std::int32_t r_max_rec{};     ///< Highest record written across all rVariables.
    std::int32_t r_num_dims{};    ///< Dimensions shared by every rVariable.
    std::int32_t nz_vars{};       ///< Number of zVariables.
    std::uint64_t uir_head{};     ///< First unused record, or 0.
    std::array<std::int32_t, kMaxDims> r_dim_sizes{};  ///< rVariable dimension sizes.
};

/// Bytes of the GDR before the rDimSizes array.
inline constexpr std::uint64_t kGdrFixedBytes = 84;

/// Parse the GDR at the offset the CDR gave.
/// @param b the file. @param offset from Cdr::gdr_offset.
/// @return the record.
inline Gdr parse_gdr(const Bytes& b, std::uint64_t offset) {
    const RecordHeader h = expect_record(b, offset, RecordType::Gdr);
    require(h.size >= kGdrFixedBytes, ErrorCode::BadRecordSize, h.offset);
    Gdr g;
    g.offset = h.offset;
    g.rvdr_head = b.be_offset(h.offset + 12);
    g.zvdr_head = b.be_offset(h.offset + 20);
    g.adr_head = b.be_offset(h.offset + 28);
    g.eof = b.be_offset(h.offset + 36);
    g.nr_vars = b.be_i32(h.offset + 44);
    g.num_attr = b.be_i32(h.offset + 48);
    g.r_max_rec = b.be_i32(h.offset + 52);
    g.r_num_dims = b.be_i32(h.offset + 56);
    g.nz_vars = b.be_i32(h.offset + 60);
    g.uir_head = b.be_offset(h.offset + 64);
    require(g.nr_vars >= 0 && g.nz_vars >= 0 && g.num_attr >= 0, ErrorCode::BadDimensions, h.offset + 44);
    require(g.r_num_dims >= 0 && g.r_num_dims <= kMaxDims, ErrorCode::BadDimensions, h.offset + 56);
    require(h.size >= kGdrFixedBytes + 4 * static_cast<std::uint64_t>(g.r_num_dims),
            ErrorCode::BadRecordSize, h.offset);
    for (std::int32_t i = 0; i < g.r_num_dims; ++i) {
        g.r_dim_sizes[static_cast<std::size_t>(i)] = b.be_i32(h.offset + kGdrFixedBytes + 4 * static_cast<std::uint64_t>(i));
        require(g.r_dim_sizes[static_cast<std::size_t>(i)] >= 1, ErrorCode::BadDimensions,
                h.offset + kGdrFixedBytes + 4 * static_cast<std::uint64_t>(i));
    }
    return g;
}

// ---- VDR -------------------------------------------------------------------------------------

/// VDR flag bit 0: the variable has one record per record number (vs. one shared record).
inline constexpr std::uint32_t kVdrFlagRecordVariance = 1U << 0;
/// VDR flag bit 1: a pad value follows the dimension arrays.
inline constexpr std::uint32_t kVdrFlagPadValue = 1U << 1;
/// VDR flag bit 2: the variable's records are compressed; `cpr_offset` names the CPR.
inline constexpr std::uint32_t kVdrFlagCompressed = 1U << 2;

/**
 * A variable descriptor — rVDR and zVDR share this shape.
 *
 * The two differ only in where the dimensions come from: a zVDR carries its own, an rVDR uses
 * the GDR's. Both are resolved at parse time so nothing downstream has to know which it was.
 */
struct Vdr {
    std::uint64_t offset{};         ///< Where it was read from.
    bool is_z{};                    ///< zVDR (true) or rVDR (false).
    std::uint64_t next{};           ///< Next VDR in the chain, or 0.
    DataType data_type{};           ///< Element type.
    std::int32_t max_rec{};         ///< Highest record number written; -1 when none.
    std::uint64_t vxr_head{};       ///< Root of the index tree, or 0 when no records.
    std::uint64_t vxr_tail{};       ///< Last VXR in the top-level chain, or 0.
    std::uint32_t flags{};          ///< The kVdrFlag* bits.
    SparseRecords s_records{};      ///< What unwritten records read as.
    std::int32_t num_elems{};       ///< Elements per value; for CHAR, the string length.
    std::int32_t num{};             ///< The variable's number within its kind.
    std::uint64_t cpr_offset{};     ///< The CPR when compressed(), else 0.
    std::int32_t blocking_factor{}; ///< Records per allocation, advisory.
    std::string_view name{};        ///< Points into the mapping.
    std::int32_t num_dims{};        ///< Dimensions, 0..kMaxDims.
    std::array<std::int32_t, kMaxDims> dim_sizes{};  ///< Size of each dimension.
    std::array<bool, kMaxDims> dim_varys{};          ///< Whether each dimension varies per record.
    std::uint64_t pad_offset{};     ///< Where the pad value sits when has_pad(), else 0.

    /// @return whether record number selects a distinct record.
    [[nodiscard]] bool record_varies() const noexcept { return (flags & kVdrFlagRecordVariance) != 0; }
    /// @return whether a pad value is present.
    [[nodiscard]] bool has_pad() const noexcept { return (flags & kVdrFlagPadValue) != 0; }
    /// @return whether the records are compressed.
    [[nodiscard]] bool compressed() const noexcept { return (flags & kVdrFlagCompressed) != 0; }
    /// @return how many records exist: max_rec + 1.
    [[nodiscard]] std::int64_t record_count() const noexcept { return static_cast<std::int64_t>(max_rec) + 1; }

    /// @return the number of VALUES in one record — the product of the varying dimensions.
    [[nodiscard]] std::uint64_t values_per_record() const noexcept {
        std::uint64_t n = 1;
        for (std::int32_t i = 0; i < num_dims; ++i) {
            if (dim_varys[static_cast<std::size_t>(i)]) { n *= static_cast<std::uint64_t>(dim_sizes[static_cast<std::size_t>(i)]); }
        }
        return n;
    }
    /// @return bytes in one stored value: element size times NumElems.
    [[nodiscard]] std::uint64_t value_bytes() const noexcept {
        return static_cast<std::uint64_t>(element_size(data_type)) * static_cast<std::uint64_t>(num_elems);
    }
    /// @return bytes in one record.
    [[nodiscard]] std::uint64_t record_bytes() const noexcept { return values_per_record() * value_bytes(); }
};

/// Parse an rVDR or zVDR. The kind is taken from the record type, not asked for.
/// @param b the file. @param offset where the VDR starts. @param gdr the GDR, for rVDR dimensions.
/// @return the descriptor.
inline Vdr parse_vdr(const Bytes& b, std::uint64_t offset, const Gdr& gdr) {
    const RecordHeader h = read_record_header(b, offset);
    require(h.type == RecordType::ZVariableDescriptor || h.type == RecordType::RVariableDescriptor,
            ErrorCode::UnexpectedRecordType, offset + 8);
    Vdr v;
    v.offset = h.offset;
    v.is_z = (h.type == RecordType::ZVariableDescriptor);
    require(h.size >= kVdrDimsOffset, ErrorCode::BadRecordSize, h.offset);
    v.next = b.be_offset(h.offset + 12);
    v.data_type = static_cast<DataType>(b.be_i32(h.offset + 20));
    require(is_known_data_type(v.data_type), ErrorCode::BadDataType, h.offset + 20);
    v.max_rec = b.be_i32(h.offset + 24);
    require(v.max_rec >= -1, ErrorCode::BadDimensions, h.offset + 24);
    v.vxr_head = b.be_offset(h.offset + 28);
    v.vxr_tail = b.be_offset(h.offset + 36);
    v.flags = b.be_u32(h.offset + 44);
    v.s_records = static_cast<SparseRecords>(b.be_i32(h.offset + 48));
    require(v.s_records == SparseRecords::None || v.s_records == SparseRecords::Pad
                || v.s_records == SparseRecords::Previous,
            ErrorCode::UnsupportedSparse, h.offset + 48);
    v.num_elems = b.be_i32(h.offset + 64);
    require(v.num_elems >= 1, ErrorCode::BadDimensions, h.offset + 64);
    v.num = b.be_i32(h.offset + 68);
    // -1 when absent — NOT 0 — so only treat it as an offset when the flag says there is one.
    if (v.compressed()) { v.cpr_offset = b.be_offset(h.offset + 72); }
    v.blocking_factor = b.be_i32(h.offset + 80);
    v.name = b.name(h.offset + kVdrNameOffset, kNameBytes);

    std::uint64_t cursor = h.offset + kVdrDimsOffset;
    if (v.is_z) {
        v.num_dims = b.be_i32(cursor);
        cursor += 4;
        require(v.num_dims >= 0 && v.num_dims <= kMaxDims, ErrorCode::BadDimensions, cursor - 4);
        // Size and vary words for every dimension must fit inside THIS record before any is
        // read — the bytes past a short record belong to whatever follows it, not to us.
        require(cursor + 8 * static_cast<std::uint64_t>(v.num_dims) <= h.offset + h.size,
                ErrorCode::BadRecordSize, cursor - 4);
        for (std::int32_t i = 0; i < v.num_dims; ++i) {
            v.dim_sizes[static_cast<std::size_t>(i)] = b.be_i32(cursor);
            require(v.dim_sizes[static_cast<std::size_t>(i)] >= 1, ErrorCode::BadDimensions, cursor);
            cursor += 4;
        }
    } else {
        v.num_dims = gdr.r_num_dims;
        v.dim_sizes = gdr.r_dim_sizes;
        require(cursor + 4 * static_cast<std::uint64_t>(v.num_dims) <= h.offset + h.size,
                ErrorCode::BadRecordSize, cursor);
    }
    for (std::int32_t i = 0; i < v.num_dims; ++i) {
        v.dim_varys[static_cast<std::size_t>(i)] = (b.be_i32(cursor) != 0);  // -1 is true
        cursor += 4;
    }
    if (v.has_pad()) {
        v.pad_offset = cursor;
        require(cursor + v.value_bytes() <= h.offset + h.size, ErrorCode::BadRecordSize, cursor);
    }
    require(cursor <= h.offset + h.size, ErrorCode::BadRecordSize, cursor);
    return v;
}

// ---- VXR -------------------------------------------------------------------------------------

/// Bytes of a VXR before its three entry arrays.
inline constexpr std::uint64_t kVxrFixedBytes = 28;

/**
 * A Variable Index Record. The entry arrays are read on demand rather than copied.
 *
 * An entry's offset points at a VVR, a CVVR, or another VXR — the record carries no tag saying
 * which, so the reader dispatches on the RecordType found at the target. See index.hpp.
 */
struct Vxr {
    std::uint64_t offset{};    ///< Where it was read from.
    std::uint64_t next{};      ///< Next sibling VXR, or 0.
    std::int32_t nentries{};   ///< Capacity of the arrays.
    std::int32_t nused{};      ///< How many entries are live.

    /// @param b the file. @param i entry index, < nused. @return the first record it covers.
    [[nodiscard]] std::int32_t first(const Bytes& b, std::int32_t i) const {
        return b.be_i32(offset + kVxrFixedBytes + 4 * static_cast<std::uint64_t>(i));
    }
    /// @param b the file. @param i entry index. @return the last record it covers, inclusive.
    [[nodiscard]] std::int32_t last(const Bytes& b, std::int32_t i) const {
        return b.be_i32(offset + kVxrFixedBytes + 4 * static_cast<std::uint64_t>(nentries)
                        + 4 * static_cast<std::uint64_t>(i));
    }
    /// @param b the file. @param i entry index. @return where the entry's record lives.
    [[nodiscard]] std::uint64_t child_offset(const Bytes& b, std::int32_t i) const {
        return b.be_offset(offset + kVxrFixedBytes + 8 * static_cast<std::uint64_t>(nentries)
                           + 8 * static_cast<std::uint64_t>(i));
    }
};

/// Parse a VXR.
/// @param b the file. @param offset where the VXR starts.
/// @return the record.
inline Vxr parse_vxr(const Bytes& b, std::uint64_t offset) {
    const RecordHeader h = expect_record(b, offset, RecordType::Vxr);
    require(h.size >= kVxrFixedBytes, ErrorCode::BadRecordSize, h.offset);
    Vxr x;
    x.offset = h.offset;
    x.next = b.be_offset(h.offset + 12);
    x.nentries = b.be_i32(h.offset + 20);
    x.nused = b.be_i32(h.offset + 24);
    require(x.nentries >= 1 && x.nentries <= kMaxVxrEntries, ErrorCode::BadRecordSize, h.offset + 20);
    require(x.nused >= 0 && x.nused <= x.nentries, ErrorCode::BadRecordSize, h.offset + 24);
    require(h.size >= kVxrFixedBytes + 16 * static_cast<std::uint64_t>(x.nentries),
            ErrorCode::BadRecordSize, h.offset);
    return x;
}

// ---- VVR / CVVR -------------------------------------------------------------------------------

/// The payload of a VVR: everything after the 12-byte header.
/// @param b the file. @param h a header already known to be a VVR.
/// @return the record data.
inline Bytes vvr_data(const Bytes& b, const RecordHeader& h) {
    return b.subspan(h.offset + kRecordHeaderBytes, h.size - kRecordHeaderBytes);
}

/// Bytes of a CVVR before the compressed payload.
inline constexpr std::uint64_t kCvvrFixedBytes = 24;

/// A Compressed Variable Values Record: where its payload is and how long it is.
struct Cvvr {
    std::uint64_t offset{};    ///< Where it was read from.
    std::uint64_t csize{};     ///< Compressed payload length.
    Bytes payload{};           ///< The compressed bytes.
};

/// Parse a CVVR.
/// @param b the file. @param offset where the CVVR starts.
/// @return the record.
inline Cvvr parse_cvvr(const Bytes& b, std::uint64_t offset) {
    const RecordHeader h = expect_record(b, offset, RecordType::Cvvr);
    require(h.size >= kCvvrFixedBytes, ErrorCode::BadRecordSize, h.offset);
    const std::int64_t raw = b.be_i64(h.offset + 16);
    require(raw >= 0, ErrorCode::BadRecordSize, h.offset + 16);
    Cvvr c;
    c.offset = h.offset;
    c.csize = static_cast<std::uint64_t>(raw);
    require(kCvvrFixedBytes + c.csize <= h.size, ErrorCode::BadRecordSize, h.offset + 16);
    c.payload = b.subspan(h.offset + kCvvrFixedBytes, c.csize);
    return c;
}

// ---- CPR -------------------------------------------------------------------------------------

/// Bytes of a CPR before its parameter array.
inline constexpr std::uint64_t kCprFixedBytes = 24;

/// A Compression Parameters Record.
struct Cpr {
    CompressionType type{};    ///< The algorithm.
    std::int32_t pcount{};     ///< How many parameters follow.
    std::int32_t level{};      ///< The first parameter — GZIP's level — or 0.
};

/// Parse a CPR.
/// @param b the file. @param offset from Vdr::cpr_offset.
/// @return the record.
inline Cpr parse_cpr(const Bytes& b, std::uint64_t offset) {
    const RecordHeader h = expect_record(b, offset, RecordType::Cpr);
    require(h.size >= kCprFixedBytes, ErrorCode::BadRecordSize, h.offset);
    Cpr c;
    c.type = static_cast<CompressionType>(b.be_i32(h.offset + 12));
    require(compression_name(c.type) != "unknown", ErrorCode::UnsupportedCompression, h.offset + 12);
    c.pcount = b.be_i32(h.offset + 20);
    require(c.pcount >= 0 && c.pcount <= kMaxCprParams, ErrorCode::BadRecordSize, h.offset + 20);
    require(h.size >= kCprFixedBytes + 4 * static_cast<std::uint64_t>(c.pcount), ErrorCode::BadRecordSize, h.offset);
    if (c.pcount > 0) { c.level = b.be_i32(h.offset + kCprFixedBytes); }
    return c;
}

}  // namespace cheatah::space::cdf::detail
