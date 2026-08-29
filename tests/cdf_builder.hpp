// A tiny in-memory CDF 3.x writer for crafted-bytes tests.
//
// Emits the smallest well-formed single-file CDF that has one zVariable with real records —
// magic, CDR, GDR, zVDR, one VXR, one VVR — with every record's offset recorded so a test can
// reach in and corrupt exactly one field. It is deliberately not a general writer: it exists
// so that each ErrorCode has a file that provokes it, and so hostile-input tests run over heap
// memory where ASan can see every byte (it cannot see an over-read inside an mmap'd page).
//
// tests/ is outside the coverage, doc and sidecar globs, so this file carries no such burden.
#pragma once

#include <cstddef>
#include <cstdint>
#include <cstring>
#include <string>
#include <vector>

#include "space/cdf/types.hpp"

namespace cdftest {

struct Blob {
    std::vector<std::byte> bytes;

    std::uint64_t size() const { return bytes.size(); }

    void put_be32(std::uint64_t at, std::uint32_t v) {
        for (int i = 3; i >= 0; --i) { bytes[at + static_cast<std::uint64_t>(3 - i)] = std::byte{static_cast<unsigned char>(v >> (8 * i))}; }
    }
    void put_be64(std::uint64_t at, std::uint64_t v) {
        for (int i = 7; i >= 0; --i) { bytes[at + static_cast<std::uint64_t>(7 - i)] = std::byte{static_cast<unsigned char>(v >> (8 * i))}; }
    }
    void put_i32(std::uint64_t at, std::int32_t v) { put_be32(at, static_cast<std::uint32_t>(v)); }
    void put_i64(std::uint64_t at, std::int64_t v) { put_be64(at, static_cast<std::uint64_t>(v)); }
    void put_f64_be(std::uint64_t at, double v) {
        std::uint64_t u{};
        std::memcpy(&u, &v, 8);
        put_be64(at, u);
    }
    std::uint64_t append(std::uint64_t n) {
        const std::uint64_t at = bytes.size();
        bytes.resize(bytes.size() + n, std::byte{0});
        return at;
    }
};

// The layout of a minimal file, so tests can name fields by offset.
struct Minimal {
    Blob blob;
    std::uint64_t cdr{}, gdr{}, zvdr{}, vxr{}, vvr{};
    std::int32_t records{};

    // Field offsets a test may want to corrupt.
    std::uint64_t cdr_gdr_offset() const { return cdr + 12; }
    std::uint64_t cdr_encoding() const { return cdr + 28; }
    std::uint64_t cdr_flags() const { return cdr + 32; }
    std::uint64_t gdr_zvdr_head() const { return gdr + 20; }
    std::uint64_t gdr_eof() const { return gdr + 36; }
    std::uint64_t gdr_nz_vars() const { return gdr + 60; }
    std::uint64_t gdr_r_num_dims() const { return gdr + 56; }
    std::uint64_t zvdr_next() const { return zvdr + 12; }
    std::uint64_t zvdr_data_type() const { return zvdr + 20; }
    std::uint64_t zvdr_max_rec() const { return zvdr + 24; }
    std::uint64_t zvdr_vxr_head() const { return zvdr + 28; }
    std::uint64_t zvdr_flags() const { return zvdr + 44; }
    std::uint64_t zvdr_s_records() const { return zvdr + 48; }
    std::uint64_t zvdr_num_elems() const { return zvdr + 64; }
    std::uint64_t zvdr_cpr_offset() const { return zvdr + 72; }
    std::uint64_t zvdr_num_dims() const { return zvdr + 340; }
    std::uint64_t vxr_next() const { return vxr + 12; }
    std::uint64_t vxr_nentries() const { return vxr + 20; }
    std::uint64_t vxr_nused() const { return vxr + 24; }
    std::uint64_t vxr_first0() const { return vxr + 28; }
    std::uint64_t vxr_last0() const { return vxr + 32; }
    std::uint64_t vxr_child0() const { return vxr + 36; }
    std::uint64_t vvr_data() const { return vvr + 12; }
};

// One zVariable of any type, its records given as raw bytes already in the file's byte order.
// `record_bytes` is the concatenation of all records; `records` is how many there are.
inline Minimal minimal_raw(const std::string& name, cheatah::space::cdf::DataType type,
                           std::int32_t num_elems, std::int32_t records,
                           const std::vector<std::byte>& record_bytes,
                           cheatah::space::cdf::Encoding enc = cheatah::space::cdf::Encoding::Network) {
    using namespace cheatah::space::cdf;
    Minimal m;
    Blob& b = m.blob;
    m.records = records;

    b.append(8);
    b.put_be32(0, kMagicV3);
    b.put_be32(4, kMagicUncompressed);

    m.cdr = b.append(312);
    b.put_i64(m.cdr + 0, 312);
    b.put_i32(m.cdr + 8, static_cast<std::int32_t>(RecordType::Cdr));
    b.put_i32(m.cdr + 20, 3);
    b.put_i32(m.cdr + 24, 9);
    b.put_i32(m.cdr + 28, static_cast<std::int32_t>(enc));
    b.put_be32(m.cdr + 32, kCdrFlagRowMajority | kCdrFlagSingleFile);

    m.gdr = b.append(84);
    b.put_i64(m.cdr + 12, static_cast<std::int64_t>(m.gdr));
    b.put_i64(m.gdr + 0, 84);
    b.put_i32(m.gdr + 8, static_cast<std::int32_t>(RecordType::Gdr));
    b.put_i32(m.gdr + 52, -1);
    b.put_i32(m.gdr + 60, 1);

    m.zvdr = b.append(344);
    b.put_i64(m.gdr + 20, static_cast<std::int64_t>(m.zvdr));
    b.put_i64(m.zvdr + 0, 344);
    b.put_i32(m.zvdr + 8, static_cast<std::int32_t>(RecordType::ZVariableDescriptor));
    b.put_i32(m.zvdr + 20, static_cast<std::int32_t>(type));
    b.put_i32(m.zvdr + 24, records - 1);
    b.put_be32(m.zvdr + 44, 1U);                       // record variance, no pad, uncompressed
    b.put_i32(m.zvdr + 64, num_elems);
    b.put_i32(m.zvdr + 68, 0);
    b.put_i64(m.zvdr + 72, -1);
    std::memcpy(b.bytes.data() + m.zvdr + 84, name.data(), name.size());

    if (records > 0) {
        m.vxr = b.append(44);
        b.put_i64(m.zvdr + 28, static_cast<std::int64_t>(m.vxr));
        b.put_i64(m.zvdr + 36, static_cast<std::int64_t>(m.vxr));
        b.put_i64(m.vxr + 0, 44);
        b.put_i32(m.vxr + 8, static_cast<std::int32_t>(RecordType::Vxr));
        b.put_i32(m.vxr + 20, 1);
        b.put_i32(m.vxr + 24, 1);
        b.put_i32(m.vxr + 28, 0);
        b.put_i32(m.vxr + 32, records - 1);

        const std::uint64_t vvr_size = 12 + record_bytes.size();
        m.vvr = b.append(vvr_size);
        b.put_i64(m.vxr + 36, static_cast<std::int64_t>(m.vvr));
        b.put_i64(m.vvr + 0, static_cast<std::int64_t>(vvr_size));
        b.put_i32(m.vvr + 8, static_cast<std::int32_t>(RecordType::Vvr));
        std::memcpy(b.bytes.data() + m.vvr + 12, record_bytes.data(), record_bytes.size());
    }
    b.put_i64(m.gdr + 36, static_cast<std::int64_t>(b.size()));
    return m;
}

// Encode a host value into `bytes` in the given byte order (big-endian for Network).
template <class T>
inline void push_value(std::vector<std::byte>& bytes, T v, bool big_endian) {
    unsigned char raw[sizeof(T)];
    std::memcpy(raw, &v, sizeof(T));
    if (big_endian) {
        for (std::size_t i = 0; i < sizeof(T); ++i) { bytes.push_back(std::byte{raw[sizeof(T) - 1 - i]}); }
    } else {
        for (std::size_t i = 0; i < sizeof(T); ++i) { bytes.push_back(std::byte{raw[i]}); }
    }
}

// One CDF_REAL8 zVariable named `name` holding values[0..n).
inline Minimal minimal_real8(const std::string& name, const std::vector<double>& values,
                             cheatah::space::cdf::Encoding enc = cheatah::space::cdf::Encoding::Network) {
    std::vector<std::byte> raw;
    for (double v : values) { push_value(raw, v, enc == cheatah::space::cdf::Encoding::Network); }
    return minimal_raw(name, cheatah::space::cdf::DataType::Real8, 1,
                       static_cast<std::int32_t>(values.size()), raw, enc);
}

}  // namespace cdftest
