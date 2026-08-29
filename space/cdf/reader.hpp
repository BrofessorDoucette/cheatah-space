#pragma once

/**
 * @file reader.hpp
 * @brief space.cdf — opening a CDF and reading its variables into ndarrays.
 *
 * The purr-facing surface of the module. A `File` is opened once — magic, CDR, GDR and every
 * variable descriptor are parsed then, so listing variables is free — and is immutable after
 * that: it can be copied, shared between threads, and read from concurrently with no locking.
 *
 * Reading a variable is one pass. The index is flattened (index.hpp), and each extent's bytes
 * are decoded (encoding.hpp) straight into the ndarray's buffer — no intermediate vector, no
 * zero-fill. The array is a **copy** of the file's data, deliberately: an ndarray owns its
 * buffer and cannot yet view foreign memory, and CDF's dominant encoding is big-endian, so on
 * an x86 host a conversion pass was mandatory regardless.
 *
 * What this reader refuses, with a typed error rather than a guess: CDF 2.x, whole-file
 * compression, multi-file CDFs, VAX float encodings, sparse records with gaps, and
 * multi-dimensional variables in column-major files. Each names why in its ErrorCode.
 *
 * @note The index describes allocation, not truth — a VVR is allocated at blocking-factor
 *       granularity and may extend past the last written record. Every read here clamps to
 *       the variable's record count. See index.hpp.
 */

// cheatah-deps: ndarray

#include "cheatah.hpp"
#include "ndarray.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

#include "bytes.hpp"
#include "encoding.hpp"
#include "index.hpp"
#include "inflate.hpp"
#include "mapping.hpp"
#include "records.hpp"
#include "types.hpp"

namespace cheatah::space::cdf {

namespace detail {

/// A ceiling on how many records one compressed extent may claim, so a corrupt index cannot ask
/// for an arbitrarily large scratch buffer before a single byte is decompressed.
inline constexpr std::uint64_t kMaxRecordsPerExtent = 1U << 24;

/// Everything an opened file knows. Built once, then only read.
struct FileImpl {
    std::unique_ptr<FileMapping> mapping;          ///< The mapping, when opened from a path.
    std::vector<std::byte> owned;                  ///< The bytes, when opened from memory.
    Bytes bytes;                                   ///< A view over whichever of the two is live.
    Cdr cdr;                                       ///< The descriptor record.
    Gdr gdr;                                       ///< The global descriptor.
    std::vector<Vdr> vars;                         ///< rVariables first, then zVariables.
    std::unordered_map<std::string, std::size_t> by_name;  ///< Name → index into vars.
    bool swap{};                                   ///< Whether values need byte-swapping here.
};

/// Parse the descriptors. Refuses everything this tranche does not read.
/// @param f the file to fill; `bytes` must already be set.
inline void load_descriptors(FileImpl& f) {
    const FileHeader fh = read_file_header(f.bytes);
    require(!fh.compressed, ErrorCode::UnsupportedCompression, 4);
    f.cdr = parse_cdr(f.bytes);
    require(f.cdr.single_file(), ErrorCode::UnsupportedMultiFile, kCdrOffset + 32);
    const EncodingClass cls = encoding_class_of(f.cdr.encoding);
    require(cls == EncodingClass::BigIeee || cls == EncodingClass::LittleIeee,
            ErrorCode::UnsupportedEncoding, kCdrOffset + 28);
    f.swap = needs_swap(cls);
    f.gdr = parse_gdr(f.bytes, f.cdr.gdr_offset);

    // Both chains, capped at the count the GDR declares — that is what stops a `next` cycle.
    const auto walk = [&](std::uint64_t head, std::int32_t declared) {
        std::int32_t seen = 0;
        for (std::uint64_t off = head; off != 0;) {
            require(seen < declared, ErrorCode::BadDimensions, off);
            Vdr v = parse_vdr(f.bytes, off, f.gdr);
            std::string name(v.name);
            require(f.by_name.emplace(name, f.vars.size()).second, ErrorCode::DuplicateVariable, off);
            off = v.next;
            f.vars.push_back(v);
            ++seen;
        }
    };
    walk(f.gdr.rvdr_head, f.gdr.nr_vars);
    walk(f.gdr.zvdr_head, f.gdr.nz_vars);
}

}  // namespace detail

/**
 * An opened CDF.
 *
 * Cheap to copy — copies share the underlying mapping — and immutable, so any number of
 * readers may use one concurrently. Obtain one from open(); a default-constructed File is
 * empty and every query on it fails with ErrorCode::EmptyFile.
 */
class File {
  public:
    /// An empty file. Every query fails until one is assigned from open().
    File() = default;

    /// Wrap a fully loaded implementation. Internal; use open().
    /// @param impl the loaded state.
    explicit File(std::shared_ptr<const detail::FileImpl> impl) noexcept : impl_(std::move(impl)) {}

    /// @return the loaded state. Throws ErrorCode::EmptyFile on a default-constructed File.
    [[nodiscard]] const detail::FileImpl& state() const {
        detail::require(impl_ != nullptr, ErrorCode::EmptyFile, 0);
        return *impl_;
    }

  private:
    std::shared_ptr<const detail::FileImpl> impl_;
};

/**
 * Open a CDF by path.
 *
 * Maps the file and parses every descriptor. Only CDF 3.x single-file, uncompressed-container,
 * IEEE-encoded files are accepted; anything else is refused here with the ErrorCode that says
 * why, so a later read cannot fail for a reason that was knowable at open.
 *
 * @param path the `.cdf` file.
 * @return the opened file.
 * @complexity O(v) in the number of variables; no record data is read.
 * @alloc one descriptor list and name index; the file itself is memory-mapped, not read.
 * @systest systests/test_cdf_read.purr
 * @par Example
 * @code{.purr}
 * import io
 * import space.cdf as cdf
 *
 * let f = cdf.open("omni_hro2_1min_20150101_v01.cdf")
 * io.print(len(cdf.var_names(f)))   # 47
 * @endcode
 */
inline File open(const std::string& path) {
    auto impl = std::make_shared<detail::FileImpl>();
    impl->mapping = std::make_unique<detail::FileMapping>(path);
    impl->bytes = impl->mapping->bytes();
    detail::load_descriptors(*impl);
    return File(std::move(impl));
}

/**
 * Open a CDF from bytes already in memory.
 *
 * The bytes are copied and owned by the File. This is the path every hostile-input test uses:
 * ASan sees each byte of a heap buffer, but cannot see an over-read that stays inside a mapped
 * page, so a bounds bug is only reliably caught here.
 *
 * @param data the file's bytes.
 * @param size how many.
 * @return the opened file.
 * @complexity O(size) to copy, then O(v) in the number of variables.
 * @alloc one copy of the bytes, plus the descriptor list and name index.
 * @systest systests/test_cdf_read.purr
 * @par Example
 * @code{.purr}
 * import io
 * import space.cdf as cdf
 *
 * # From purr, prefer cdf.open(path); this entry point exists for tests and embedders.
 * let f = cdf.open("omni_hro2_1min_20150101_v01.cdf")
 * io.print(cdf.record_count(f, "Epoch"))   # 44640
 * @endcode
 */
inline File open_from_memory(const std::byte* data, std::uint64_t size) {
    auto impl = std::make_shared<detail::FileImpl>();
    impl->owned.assign(data, data + size);
    impl->bytes = detail::Bytes(impl->owned.data(), impl->owned.size());
    detail::load_descriptors(*impl);
    return File(std::move(impl));
}

namespace detail {
/// Look a variable up by name.
/// @param f the file. @param name the variable. @return its descriptor.
inline const Vdr& find_var(const File& f, const std::string& name) {
    const FileImpl& s = f.state();
    const auto it = s.by_name.find(name);
    require(it != s.by_name.end(), ErrorCode::VariableNotFound, 0);
    return s.vars[it->second];
}

/// The shape values() returns for a variable: records, then each varying dimension, then the
/// per-value multiplicity a CHAR string or an EPOCH16 pair adds.
/// @param v the descriptor. @return the shape.
inline std::vector<std::size_t> shape_of(const Vdr& v) {
    std::vector<std::size_t> s;
    s.push_back(static_cast<std::size_t>(v.record_count()));
    for (std::int32_t i = 0; i < v.num_dims; ++i) {
        if (v.dim_varys[static_cast<std::size_t>(i)]) {
            s.push_back(static_cast<std::size_t>(v.dim_sizes[static_cast<std::size_t>(i)]));
        }
    }
    if (is_char_type(v.data_type) && v.num_elems > 1) { s.push_back(static_cast<std::size_t>(v.num_elems)); }
    if (v.data_type == DataType::Epoch16) { s.push_back(2); }
    return s;
}

/// Read every record of a variable into a freshly allocated array of Out.
/// @tparam Out double or long long. @param f the file. @param name the variable.
/// @return the array, shaped by shape_of().
template <class Out>
inline ::cheatah::ndarray::basic_ndarray<Out> read_values(const File& f, const std::string& name) {
    const FileImpl& s = f.state();
    const Vdr& v = find_var(f, name);
    require(v.num_dims <= 1 || s.cdr.row_major(), ErrorCode::UnsupportedLayout, v.offset);
    // Refuse the lossy pairings before allocating anything, so a caller that asked for the
    // wrong accessor gets the error rather than an array it then has to throw away.
    const StorageKind kind = storage_kind(v.data_type);
    if constexpr (std::is_floating_point_v<Out>) {
        require(kind != StorageKind::I64, ErrorCode::LossyConversion, v.offset);
    } else {
        require(!is_float_kind(kind), ErrorCode::LossyConversion, v.offset);
    }

    const std::int64_t n = v.record_count();
    const std::uint64_t rec_bytes = v.record_bytes();
    const std::uint64_t elems = rec_bytes / raw_element_bytes(v.data_type);
    auto out = ::cheatah::ndarray::basic_ndarray<Out>::uninitialized(shape_of(v));
    if (n == 0) { return out; }
    Out* dst = out.buffer()->data();

    const std::vector<RecordExtent> extents = build_index(s.bytes, v.vxr_head);
    const ErrorCode gap = (v.s_records == SparseRecords::None) ? ErrorCode::RecordOutOfRange
                                                               : ErrorCode::UnsupportedSparse;
    std::int64_t cursor = 0;
    for (const RecordExtent& e : extents) {
        if (e.first >= n) { break; }
        require(e.first == cursor, gap, e.offset);
        const std::int64_t hi = std::min<std::int64_t>(e.last, n - 1);
        const auto count = static_cast<std::uint64_t>(hi - e.first + 1);
        if (e.compressed) {
            // A CVVR holds the WHOLE extent compressed, so the full block is inflated even when
            // only its first records are wanted — DEFLATE has no random access, and the trailing
            // records of the last extent are allocation slack anyway.
            const auto whole = static_cast<std::uint64_t>(e.last - e.first + 1);
            require(whole <= kMaxRecordsPerExtent, ErrorCode::DecompressionFailed, e.offset);
            const Cvvr c = parse_cvvr(s.bytes, e.offset);
            std::vector<std::byte> plain(whole * rec_bytes);
            inflate(c.payload, plain.data(), plain.size(), e.offset);
            decode_run<Out>(plain.data(), v.data_type, s.swap, count * elems,
                            dst + static_cast<std::uint64_t>(cursor) * elems, e.offset);
        } else {
            const RecordHeader h = expect_record(s.bytes, e.offset, RecordType::Vvr);
            const Bytes data = vvr_data(s.bytes, h);
            require(data.size() >= count * rec_bytes, ErrorCode::BadRecordSize, e.offset);
            decode_run<Out>(data.data(), v.data_type, s.swap, count * elems,
                            dst + static_cast<std::uint64_t>(cursor) * elems, e.offset);
        }
        cursor = hi + 1;
    }
    require(cursor == n, gap, v.vxr_head);
    return out;
}
}  // namespace detail

/**
 * The names of every variable in the file, rVariables first, in file order.
 * @param f an opened file.
 * @return the names.
 * @complexity O(v).
 * @alloc one list of v strings.
 * @systest systests/test_cdf_read.purr
 * @par Example
 * @code{.purr}
 * import io
 * import space.cdf as cdf
 *
 * let f = cdf.open("omni_hro2_1min_20150101_v01.cdf")
 * io.print(cdf.var_names(f)[0])   # Epoch
 * @endcode
 */
inline std::vector<std::string> var_names(const File& f) {
    const detail::FileImpl& s = f.state();
    std::vector<std::string> names;
    names.reserve(s.vars.size());
    for (const detail::Vdr& v : s.vars) { names.emplace_back(v.name); }
    return names;
}

/**
 * How many records a variable has.
 * @param f an opened file.
 * @param name the variable.
 * @return the record count; 0 for a variable with nothing written.
 * @complexity O(1) after the name lookup.
 * @alloc none.
 * @systest systests/test_cdf_read.purr
 * @par Example
 * @code{.purr}
 * import io
 * import space.cdf as cdf
 *
 * let f = cdf.open("omni_hro2_1min_20150101_v01.cdf")
 * io.print(cdf.record_count(f, "F"))   # 44640 — one day at one-minute cadence
 * @endcode
 */
inline long long record_count(const File& f, const std::string& name) {
    return detail::find_var(f, name).record_count();
}

/**
 * The CDF data type of a variable, by its format name.
 * @param f an opened file.
 * @param name the variable.
 * @return e.g. "CDF_REAL4", "CDF_EPOCH", "CDF_TIME_TT2000".
 * @complexity O(1) after the name lookup.
 * @alloc one short string.
 * @systest systests/test_cdf_read.purr
 * @par Example
 * @code{.purr}
 * import io
 * import space.cdf as cdf
 *
 * let f = cdf.open("omni_hro2_1min_20150101_v01.cdf")
 * io.print(cdf.data_type_name(f, "Epoch"))   # CDF_EPOCH
 * @endcode
 */
inline std::string data_type_name(const File& f, const std::string& name) {
    return std::string(detail::data_type_name(detail::find_var(f, name).data_type));
}

/**
 * The shape values() will return for a variable.
 *
 * Records first, then each dimension that varies per record, then one trailing axis when a
 * value is itself several elements — the string length of a CHAR variable, or 2 for EPOCH16.
 *
 * @param f an opened file.
 * @param name the variable.
 * @return the shape.
 * @complexity O(d) in the number of dimensions.
 * @alloc one short list.
 * @systest systests/test_cdf_read.purr
 * @par Example
 * @code{.purr}
 * import io
 * import space.cdf as cdf
 *
 * let f = cdf.open("omni_hro2_1min_20150101_v01.cdf")
 * io.print(cdf.shape(f, "F"))   # [44640] — a scalar per record
 * @endcode
 */
inline std::vector<long long> shape(const File& f, const std::string& name) {
    std::vector<long long> out;
    for (std::size_t d : detail::shape_of(detail::find_var(f, name))) { out.push_back(static_cast<long long>(d)); }
    return out;
}

/**
 * Every record of a variable, as float64.
 *
 * Lossless for every CDF type except CDF_INT8 and CDF_TIME_TT2000, which a double cannot hold
 * exactly; those are refused with ErrorCode::LossyConversion — use values_i64(). CHAR data
 * comes back as byte values with the string length as the last axis; EPOCH16 as a trailing
 * axis of 2. This is the array cheatah-plot draws directly.
 *
 * @param f an opened file.
 * @param name the variable.
 * @return an ndarray shaped as shape() describes. A copy of the file's data.
 * @complexity O(n) in the number of values — one decode pass.
 * @alloc the result array, exactly once; nothing intermediate.
 * @systest systests/test_cdf_read.purr
 * @par Example
 * @code{.purr}
 * import io
 * import ndarray
 * import space.cdf as cdf
 *
 * let f = cdf.open("omni_hro2_1min_20150101_v01.cdf")
 * let b = cdf.values(f, "F")            # IMF magnitude, nT, one value per minute
 * io.print(b[0], b[1], b[2])            # 6.92 5.84 5.71
 * @endcode
 */
inline ::cheatah::ndarray::basic_ndarray<double> values(const File& f, const std::string& name) {
    return detail::read_values<double>(f, name);
}

/**
 * Every record of an integer variable, as int64 — the exact path for CDF_INT8 and TT2000.
 *
 * Accepts every integer type; refuses REAL4/REAL8/EPOCH/EPOCH16 with ErrorCode::LossyConversion
 * because truncating a float is not a conversion anyone asked for.
 *
 * @param f an opened file.
 * @param name the variable.
 * @return an ndarray of int shaped as shape() describes. A copy of the file's data.
 * @complexity O(n) in the number of values — one decode pass.
 * @alloc the result array, exactly once; nothing intermediate.
 * @systest systests/test_cdf_read.purr
 * @par Example
 * @code{.purr}
 * import io
 * import ndarray
 * import space.cdf as cdf
 *
 * # TT2000 is nanoseconds since J2000; an int64 holds it exactly, a double does not.
 * let f = cdf.open("rbsp-a_density_emfisis-l4_20130101_v1.5.17.cdf")
 * let t = cdf.values_i64(f, "Epoch")
 * io.print(t[0])
 * @endcode
 */
inline ::cheatah::ndarray::basic_ndarray<long long> values_i64(const File& f, const std::string& name) {
    return detail::read_values<long long>(f, name);
}

}  // namespace cheatah::space::cdf
