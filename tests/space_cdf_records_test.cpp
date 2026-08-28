// Unit tests for space/cdf/records.hpp and index.hpp — parsing the internal records of a CDF
// and flattening a variable's index tree.
//
// The real-file assertions are against values decoded by hand from the raw bytes before this
// code existed (OMNI: 47 zVariables, Epoch's index has exactly 24 leaves in a two-level tree;
// test_alltypes: `volume` is [2,4,2], `Longitude` is GZIP level 9). The crafted-bytes tests use
// the in-memory builder so that every ErrorCode the parsers can raise has a file that raises it
// — and so the walk's cycle caps are shown to fire, not assumed to.
#include <gtest/gtest.h>

#include <cstring>
#include <fstream>
#include <functional>
#include <string>
#include <vector>

#include "cdf_builder.hpp"
#include "space/space.hpp"

namespace cdf = cheatah::space::cdf;
namespace det = cheatah::space::cdf::detail;

namespace {

std::string corpus_path(const char* rel) {
    std::string prefix;
    for (int up = 0; up < 5; ++up) {
        const std::string candidate = prefix + "space/cdf/vendor/corpus/" + rel;
        if (std::ifstream(candidate, std::ios::binary).good()) { return candidate; }
        prefix += "../";
    }
    return {};
}

std::vector<std::byte> read_all(const std::string& path) {
    std::ifstream in(path, std::ios::binary);
    const std::vector<char> raw((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
    std::vector<std::byte> out(raw.size());
    std::memcpy(out.data(), raw.data(), raw.size());
    return out;
}

det::Bytes view(const cdftest::Minimal& m) { return {m.blob.bytes.data(), m.blob.size()}; }

cdf::ErrorCode code_of(const std::function<void()>& f) {
    try { f(); } catch (const cdf::CdfError& e) { return e.code(); }
    return cdf::ErrorCode::None;
}

}  // namespace

// ---- a real file: OMNI --------------------------------------------------------------------------

TEST(CdfRecords, OmniDescriptorsDecodeAsHandDecoded) {
    const std::string path = corpus_path("tier1/omni_hro2_1min_20150101_v01.cdf");
    if (path.empty()) { GTEST_SKIP() << "run: scripts/cdf-corpus.sh fetch --tier 1"; }
    const auto buf = read_all(path);
    const det::Bytes b(buf.data(), buf.size());

    const det::FileHeader fh = det::read_file_header(b);
    EXPECT_FALSE(fh.compressed);

    const det::Cdr cdr = det::parse_cdr(b);
    EXPECT_EQ(cdr.gdr_offset, 320u);
    EXPECT_EQ(cdr.version, 3);
    EXPECT_EQ(cdr.release, 9);
    EXPECT_EQ(cdr.encoding, cdf::Encoding::Network);
    EXPECT_TRUE(cdr.row_major());
    EXPECT_TRUE(cdr.single_file());
    EXPECT_FALSE(cdr.has_checksum());

    const det::Gdr gdr = det::parse_gdr(b, cdr.gdr_offset);
    EXPECT_EQ(gdr.zvdr_head, 21601u);
    EXPECT_EQ(gdr.rvdr_head, 0u);
    EXPECT_EQ(gdr.nz_vars, 47);
    EXPECT_EQ(gdr.nr_vars, 0);
    EXPECT_EQ(gdr.num_attr, 47);
    EXPECT_EQ(gdr.r_num_dims, 0);
    EXPECT_EQ(gdr.eof, b.size());  // no checksum => eof is the whole file

    // Walk the whole zVDR chain: exactly NzVars descriptors, first is Epoch.
    int count = 0;
    std::uint64_t off = gdr.zvdr_head;
    det::Vdr first{};
    while (off != 0) {
        const det::Vdr v = det::parse_vdr(b, off, gdr);
        if (count == 0) { first = v; }
        EXPECT_TRUE(v.is_z);
        EXPECT_EQ(v.num, count);
        ++count;
        off = v.next;
    }
    EXPECT_EQ(count, 47);
    EXPECT_EQ(first.name, "Epoch");
    EXPECT_EQ(first.data_type, cdf::DataType::Epoch);
    EXPECT_EQ(first.max_rec, 44639);
    EXPECT_EQ(first.record_count(), 44640);
    EXPECT_EQ(first.num_elems, 1);
    EXPECT_EQ(first.num_dims, 0);
    EXPECT_EQ(first.values_per_record(), 1u);
    EXPECT_EQ(first.record_bytes(), 8u);
    EXPECT_TRUE(first.record_varies());
    EXPECT_FALSE(first.compressed());
    EXPECT_EQ(first.cpr_offset, 0u);
    EXPECT_EQ(first.vxr_head, 5607364u);
}

TEST(CdfIndex, OmniEpochIndexIsATwoLevelTreeWith24Leaves) {
    const std::string path = corpus_path("tier1/omni_hro2_1min_20150101_v01.cdf");
    if (path.empty()) { GTEST_SKIP() << "run: scripts/cdf-corpus.sh fetch --tier 1"; }
    const auto buf = read_all(path);
    const det::Bytes b(buf.data(), buf.size());
    const det::Gdr gdr = det::parse_gdr(b, det::parse_cdr(b).gdr_offset);
    const det::Vdr epoch = det::parse_vdr(b, gdr.zvdr_head, gdr);

    const std::vector<det::RecordExtent> ext = det::build_index(b, epoch.vxr_head);
    ASSERT_EQ(ext.size(), 24u);
    // Sorted, contiguous, and covering exactly 0..maxRec with no gaps or overlaps.
    std::int64_t expect_first = 0;
    for (const det::RecordExtent& e : ext) {
        EXPECT_EQ(e.first, expect_first);
        EXPECT_GE(e.last, e.first);
        EXPECT_FALSE(e.compressed);
        expect_first = e.last + 1;
    }
    // The index describes ALLOCATION, not truth: the last VVR is a full 1024-record block
    // (44032..45055) although only records through maxRec 44639 were written. A reader must
    // clamp to Vdr::record_count(); the index alone would hand back 416 records of garbage.
    EXPECT_EQ(ext.back().first, 44032);
    EXPECT_EQ(ext.back().last, 45055);
    EXPECT_GE(ext.back().last, epoch.max_rec);
    EXPECT_EQ(ext.front().offset, 81272u);   // the first VVR, as walked by hand
    EXPECT_EQ(ext.front().count(), 1024);    // 8204-byte VVR = 12 + 1024 * 8

    // Lookups land on the right extent, including both ends of an extent and past the end.
    EXPECT_EQ(det::find_extent(ext, 0), ext.data());
    EXPECT_EQ(det::find_extent(ext, 1023), ext.data());
    EXPECT_EQ(det::find_extent(ext, 1024), &ext[1]);
    EXPECT_EQ(det::find_extent(ext, 44639), &ext[23]);
    EXPECT_EQ(det::find_extent(ext, 45055), &ext[23]);
    EXPECT_EQ(det::find_extent(ext, 45056), nullptr);
    EXPECT_EQ(det::find_extent(ext, -1), nullptr);

    // The leaf's VVR holds what the index says it holds.
    const det::RecordHeader h = det::expect_record(b, ext[0].offset, cdf::RecordType::Vvr);
    EXPECT_EQ(det::vvr_data(b, h).size(), 1024u * 8u);
}

// ---- a real file: test_alltypes ---------------------------------------------------------------------

TEST(CdfRecords, TestAlltypesDimensionsPadsAndCompression) {
    const std::string path = corpus_path("tier0/test_alltypes.cdf");
    if (path.empty()) { GTEST_SKIP() << "run: scripts/cdf-corpus.sh fetch"; }
    const auto buf = read_all(path);
    const det::Bytes b(buf.data(), buf.size());
    const det::Cdr cdr = det::parse_cdr(b);
    EXPECT_EQ(cdr.encoding, cdf::Encoding::IbmPc);
    EXPECT_TRUE(cdr.has_checksum());
    const det::Gdr gdr = det::parse_gdr(b, cdr.gdr_offset);
    EXPECT_EQ(gdr.eof + cdf::kChecksumBytes, b.size());  // checksum => 16 bytes past eof

    bool saw_volume = false, saw_longitude = false, saw_name = false, saw_temp = false;
    for (std::uint64_t off = gdr.zvdr_head; off != 0;) {
        const det::Vdr v = det::parse_vdr(b, off, gdr);
        if (v.name == "volume") {
            saw_volume = true;
            EXPECT_EQ(v.num_dims, 3);
            EXPECT_EQ(v.dim_sizes[0], 2); EXPECT_EQ(v.dim_sizes[1], 4); EXPECT_EQ(v.dim_sizes[2], 2);
            EXPECT_TRUE(v.dim_varys[0] && v.dim_varys[1] && v.dim_varys[2]);
            EXPECT_EQ(v.values_per_record(), 16u);
            EXPECT_EQ(v.record_bytes(), 64u);   // INT4
            EXPECT_TRUE(v.has_pad());
            EXPECT_NE(v.pad_offset, 0u);
        } else if (v.name == "Longitude") {
            saw_longitude = true;
            EXPECT_TRUE(v.compressed());
            EXPECT_EQ(v.cpr_offset, 11098u);
            EXPECT_EQ(v.blocking_factor, 130);
            const det::Cpr cpr = det::parse_cpr(b, v.cpr_offset);
            EXPECT_EQ(cpr.type, cdf::CompressionType::Gzip);
            EXPECT_EQ(cpr.pcount, 1);
            EXPECT_EQ(cpr.level, 9);
            const auto ext = det::build_index(b, v.vxr_head);
            ASSERT_FALSE(ext.empty());
            EXPECT_TRUE(ext[0].compressed);
            const det::Cvvr c = det::parse_cvvr(b, ext[0].offset);
            EXPECT_EQ(c.csize, c.payload.size());
            EXPECT_EQ(c.payload.u8(0), 0x1Fu);  // a gzip container, as the corpus probe found
            EXPECT_EQ(c.payload.u8(1), 0x8Bu);
        } else if (v.name == "Name") {
            saw_name = true;
            EXPECT_EQ(v.data_type, cdf::DataType::Char);
            EXPECT_EQ(v.num_elems, 10);       // the string length
            EXPECT_EQ(v.value_bytes(), 10u);
        } else if (v.name == "Temp") {
            saw_temp = true;
            EXPECT_EQ(v.s_records, cdf::SparseRecords::Pad);
        }
        off = v.next;
    }
    EXPECT_TRUE(saw_volume && saw_longitude && saw_name && saw_temp);
}

// ---- crafted bytes: the happy path of the builder itself --------------------------------------------

TEST(CdfRecords, MinimalBuiltFileParsesEndToEnd) {
    const cdftest::Minimal m = cdftest::minimal_real8("x", {1.0, 2.0, 3.0});
    const det::Bytes b = view(m);
    EXPECT_FALSE(det::read_file_header(b).compressed);
    const det::Cdr cdr = det::parse_cdr(b);
    EXPECT_EQ(cdr.gdr_offset, m.gdr);
    const det::Gdr gdr = det::parse_gdr(b, cdr.gdr_offset);
    EXPECT_EQ(gdr.nz_vars, 1);
    EXPECT_EQ(gdr.eof, b.size());
    const det::Vdr v = det::parse_vdr(b, gdr.zvdr_head, gdr);
    EXPECT_EQ(v.name, "x");
    EXPECT_EQ(v.record_count(), 3);
    EXPECT_EQ(v.next, 0u);
    const auto ext = det::build_index(b, v.vxr_head);
    ASSERT_EQ(ext.size(), 1u);
    EXPECT_EQ(ext[0].first, 0);
    EXPECT_EQ(ext[0].last, 2);
    EXPECT_EQ(ext[0].offset, m.vvr);

    // Zero records: no VXR at all, and the index is empty rather than an error.
    const cdftest::Minimal empty = cdftest::minimal_real8("e", {});
    const det::Bytes eb = view(empty);
    const det::Vdr ev = det::parse_vdr(eb, det::parse_gdr(eb, det::parse_cdr(eb).gdr_offset).zvdr_head, det::parse_gdr(eb, det::parse_cdr(eb).gdr_offset));
    EXPECT_EQ(ev.record_count(), 0);
    EXPECT_EQ(ev.vxr_head, 0u);
    EXPECT_TRUE(det::build_index(eb, ev.vxr_head).empty());
}

// ---- crafted bytes: every error the parsers can raise --------------------------------------------------

TEST(CdfRecords, FileHeaderErrors) {
    cdftest::Minimal m = cdftest::minimal_real8("x", {1.0});
    m.blob.put_be32(0, cdf::kMagicV26);
    EXPECT_EQ(code_of([&] { det::read_file_header(view(m)); }), cdf::ErrorCode::UnsupportedPreV26);
    m.blob.put_be32(0, cdf::kMagicPreV26);
    EXPECT_EQ(code_of([&] { det::read_file_header(view(m)); }), cdf::ErrorCode::UnsupportedPreV26);
    m.blob.put_be32(0, 0x12345678u);
    EXPECT_EQ(code_of([&] { det::read_file_header(view(m)); }), cdf::ErrorCode::NotCdf);
    m.blob.put_be32(0, cdf::kMagicV3);
    m.blob.put_be32(4, 0xDEADBEEFu);
    EXPECT_EQ(code_of([&] { det::read_file_header(view(m)); }), cdf::ErrorCode::BadMagic);
    m.blob.put_be32(4, cdf::kMagicCompressed);
    EXPECT_TRUE(det::read_file_header(view(m)).compressed);

    const std::byte tiny[4] = {};
    EXPECT_EQ(code_of([&] { det::read_file_header(det::Bytes(tiny, 4)); }), cdf::ErrorCode::EmptyFile);
}

TEST(CdfRecords, RecordHeaderErrors) {
    cdftest::Minimal m = cdftest::minimal_real8("x", {1.0});
    // Size smaller than the header itself.
    m.blob.put_i64(m.cdr, 4);
    EXPECT_EQ(code_of([&] { det::parse_cdr(view(m)); }), cdf::ErrorCode::BadRecordSize);
    // Size running past the end of the file.
    m.blob.put_i64(m.cdr, static_cast<std::int64_t>(m.blob.size()));
    EXPECT_EQ(code_of([&] { det::parse_cdr(view(m)); }), cdf::ErrorCode::BadRecordSize);
    // Negative size.
    m.blob.put_i64(m.cdr, -1);
    EXPECT_EQ(code_of([&] { det::parse_cdr(view(m)); }), cdf::ErrorCode::BadRecordSize);
    m.blob.put_i64(m.cdr, 312);
    // Undefined record type.
    m.blob.put_i32(m.cdr + 8, 99);
    EXPECT_EQ(code_of([&] { det::parse_cdr(view(m)); }), cdf::ErrorCode::BadRecordType);
    // Defined, but not the type expected here.
    m.blob.put_i32(m.cdr + 8, static_cast<std::int32_t>(cdf::RecordType::Gdr));
    EXPECT_EQ(code_of([&] { det::parse_cdr(view(m)); }), cdf::ErrorCode::UnexpectedRecordType);
    m.blob.put_i32(m.cdr + 8, static_cast<std::int32_t>(cdf::RecordType::Cdr));
    // CDR too short for its fixed fields — but a legal record size.
    m.blob.put_i64(m.cdr, 40);
    EXPECT_EQ(code_of([&] { det::parse_cdr(view(m)); }), cdf::ErrorCode::BadRecordSize);
    m.blob.put_i64(m.cdr, 312);
    // Truncated mid-record.
    const det::Bytes cut(m.blob.bytes.data(), m.cdr + 30);
    EXPECT_EQ(code_of([&] { det::parse_cdr(cut); }), cdf::ErrorCode::BadRecordSize);
}

TEST(CdfRecords, CdrAndGdrFieldErrors) {
    cdftest::Minimal m = cdftest::minimal_real8("x", {1.0});
    m.blob.put_i32(m.cdr_encoding(), 99);
    EXPECT_EQ(code_of([&] { det::parse_cdr(view(m)); }), cdf::ErrorCode::BadEncoding);
    m.blob.put_i32(m.cdr_encoding(), static_cast<std::int32_t>(cdf::Encoding::Network));
    m.blob.put_i64(m.cdr_gdr_offset(), static_cast<std::int64_t>(m.blob.size() + 100));
    EXPECT_EQ(code_of([&] { det::parse_cdr(view(m)); }), cdf::ErrorCode::BadOffset);
    m.blob.put_i64(m.cdr_gdr_offset(), static_cast<std::int64_t>(m.gdr));

    const det::Bytes b = view(m);
    m.blob.put_i32(m.gdr_nz_vars(), -3);
    EXPECT_EQ(code_of([&] { det::parse_gdr(b, m.gdr); }), cdf::ErrorCode::BadDimensions);
    m.blob.put_i32(m.gdr_nz_vars(), 1);
    m.blob.put_i32(m.gdr_r_num_dims(), 11);
    EXPECT_EQ(code_of([&] { det::parse_gdr(b, m.gdr); }), cdf::ErrorCode::BadDimensions);
    // A dim count the 84-byte record cannot hold.
    m.blob.put_i32(m.gdr_r_num_dims(), 2);
    EXPECT_EQ(code_of([&] { det::parse_gdr(b, m.gdr); }), cdf::ErrorCode::BadRecordSize);
    m.blob.put_i32(m.gdr_r_num_dims(), 0);
    m.blob.put_i64(m.gdr_eof(), -5);
    EXPECT_EQ(code_of([&] { det::parse_gdr(b, m.gdr); }), cdf::ErrorCode::BadOffset);
    m.blob.put_i64(m.gdr_eof(), static_cast<std::int64_t>(m.blob.size()));
    // An rVDR with GDR dims: build a GDR carrying one dim of size 0 to hit the size check.
    cdftest::Minimal r = cdftest::minimal_real8("r", {1.0});
    r.blob.bytes.resize(r.blob.size() + 4, std::byte{0});          // room for rDimSizes[0]
    r.blob.put_i64(r.gdr, 88);                                        // GDR grows by 4 …
    // … but the zVDR now starts 4 bytes later than the GDR says: rebuild that link.
    r.blob.put_i32(r.gdr_r_num_dims(), 1);
    r.blob.put_i32(r.gdr + 84, 0);                                    // size 0 is invalid
    EXPECT_EQ(code_of([&] { det::parse_gdr(view(r), r.gdr); }), cdf::ErrorCode::BadDimensions);
}

TEST(CdfRecords, VdrFieldErrors) {
    cdftest::Minimal m = cdftest::minimal_real8("x", {1.0, 2.0});
    const det::Bytes b = view(m);
    const det::Gdr gdr = det::parse_gdr(b, det::parse_cdr(b).gdr_offset);
    auto parse = [&] { (void)det::parse_vdr(view(m), m.zvdr, gdr); };

    m.blob.put_i32(m.zvdr + 8, static_cast<std::int32_t>(cdf::RecordType::Vxr));
    EXPECT_EQ(code_of(parse), cdf::ErrorCode::UnexpectedRecordType);
    m.blob.put_i32(m.zvdr + 8, static_cast<std::int32_t>(cdf::RecordType::ZVariableDescriptor));

    m.blob.put_i32(m.zvdr_data_type(), 7777);
    EXPECT_EQ(code_of(parse), cdf::ErrorCode::BadDataType);
    m.blob.put_i32(m.zvdr_data_type(), static_cast<std::int32_t>(cdf::DataType::Real8));

    m.blob.put_i32(m.zvdr_max_rec(), -2);
    EXPECT_EQ(code_of(parse), cdf::ErrorCode::BadDimensions);
    m.blob.put_i32(m.zvdr_max_rec(), 1);

    m.blob.put_i32(m.zvdr_s_records(), 9);
    EXPECT_EQ(code_of(parse), cdf::ErrorCode::UnsupportedSparse);
    m.blob.put_i32(m.zvdr_s_records(), 0);

    m.blob.put_i32(m.zvdr_num_elems(), 0);
    EXPECT_EQ(code_of(parse), cdf::ErrorCode::BadDimensions);
    m.blob.put_i32(m.zvdr_num_elems(), 1);

    m.blob.put_i32(m.zvdr_num_dims(), 11);
    EXPECT_EQ(code_of(parse), cdf::ErrorCode::BadDimensions);
    m.blob.put_i32(m.zvdr_num_dims(), 0);

    // Compressed flag set but the CPR offset is the -1 sentinel: that is a bad offset, not a CPR.
    m.blob.put_be32(m.zvdr_flags(), det::kVdrFlagRecordVariance | det::kVdrFlagCompressed);
    EXPECT_EQ(code_of(parse), cdf::ErrorCode::BadOffset);
    // Pad flag set but the record has no room for the pad value.
    m.blob.put_be32(m.zvdr_flags(), det::kVdrFlagRecordVariance | det::kVdrFlagPadValue);
    EXPECT_EQ(code_of(parse), cdf::ErrorCode::BadRecordSize);
    m.blob.put_be32(m.zvdr_flags(), det::kVdrFlagRecordVariance);

    // A zVDR claiming one dimension in a record too short to hold its size and vary words. The
    // bytes past the record are the VXR's — inside the file, so this is a record-size error, not
    // a truncation, and the parser must refuse before reading them.
    m.blob.put_i32(m.zvdr_num_dims(), 1);
    EXPECT_EQ(code_of(parse), cdf::ErrorCode::BadRecordSize);
    m.blob.put_i32(m.zvdr_num_dims(), 0);

    // The VDR record too short for even its fixed fields.
    m.blob.put_i64(m.zvdr, 100);
    EXPECT_EQ(code_of(parse), cdf::ErrorCode::BadRecordSize);
    m.blob.put_i64(m.zvdr, 344);
    EXPECT_EQ(code_of(parse), cdf::ErrorCode::None);
}

TEST(CdfRecords, VdrWithDimensionsAndPadParses) {
    // Grow the minimal zVDR into one with a [3] dimension and a pad value, then parse it.
    cdftest::Minimal m = cdftest::minimal_real8("d", {1.0});
    // 344 + 4 (size) + 4 (vary) + 8 (pad) = 360; move the VXR/VVR by rebuilding is overkill —
    // instead point the zVDR's next at 0 and let the record simply be larger than the fields
    // that follow it need; the VXR is reachable by offset, not by adjacency.
    std::vector<std::byte>& bytes = m.blob.bytes;
    bytes.insert(bytes.begin() + static_cast<std::ptrdiff_t>(m.zvdr + 344), 16, std::byte{0});
    m.blob.put_i64(m.zvdr, 360);
    m.blob.put_i32(m.zvdr_num_dims(), 1);
    m.blob.put_i32(m.zvdr + 344, 3);          // zDimSizes[0]
    m.blob.put_i32(m.zvdr + 348, -1);         // DimVarys[0] = true
    m.blob.put_be32(m.zvdr_flags(), det::kVdrFlagRecordVariance | det::kVdrFlagPadValue);
    m.blob.put_f64_be(m.zvdr + 352, -1.0e31);  // the pad
    // Every record after the zVDR moved by 16: fix the offsets that point at them.
    m.blob.put_i64(m.zvdr_vxr_head(), static_cast<std::int64_t>(m.vxr + 16));
    m.blob.put_i64(m.zvdr + 36, static_cast<std::int64_t>(m.vxr + 16));
    m.blob.put_i64(m.vxr + 16 + 36, static_cast<std::int64_t>(m.vvr + 16));
    m.blob.put_i64(m.gdr_eof(), static_cast<std::int64_t>(m.blob.size()));

    const det::Bytes b = view(m);
    const det::Vdr v = det::parse_vdr(b, m.zvdr, det::parse_gdr(b, det::parse_cdr(b).gdr_offset));
    EXPECT_EQ(v.num_dims, 1);
    EXPECT_EQ(v.dim_sizes[0], 3);
    EXPECT_TRUE(v.dim_varys[0]);
    EXPECT_EQ(v.values_per_record(), 3u);
    EXPECT_EQ(v.record_bytes(), 24u);
    EXPECT_TRUE(v.has_pad());
    EXPECT_EQ(v.pad_offset, m.zvdr + 352);
    EXPECT_EQ(det::build_index(b, v.vxr_head).size(), 1u);

    // A non-varying dimension does not multiply the record size.
    m.blob.put_i32(m.zvdr + 348, 0);
    EXPECT_EQ(det::parse_vdr(view(m), m.zvdr, det::parse_gdr(view(m), m.gdr)).values_per_record(), 1u);
    // A zero dimension size is rejected.
    m.blob.put_i32(m.zvdr + 344, 0);
    EXPECT_EQ(code_of([&] { det::parse_vdr(view(m), m.zvdr, det::parse_gdr(view(m), m.gdr)); }),
              cdf::ErrorCode::BadDimensions);
}

TEST(CdfRecords, RVdrTakesItsDimensionsFromTheGdr) {
    // Turn the minimal zVDR into an rVDR: type 3, and dims come from the GDR (0 here), so the
    // 4 bytes at +340 are read as DimVarys for zero dims — i.e. not read at all.
    cdftest::Minimal m = cdftest::minimal_real8("r", {5.0});
    m.blob.put_i32(m.zvdr + 8, static_cast<std::int32_t>(cdf::RecordType::RVariableDescriptor));
    const det::Bytes b = view(m);
    const det::Gdr gdr = det::parse_gdr(b, det::parse_cdr(b).gdr_offset);
    const det::Vdr v = det::parse_vdr(b, m.zvdr, gdr);
    EXPECT_FALSE(v.is_z);
    EXPECT_EQ(v.num_dims, 0);
    EXPECT_EQ(v.name, "r");
    EXPECT_EQ(v.record_bytes(), 8u);
}

TEST(CdfRecords, VxrCvvrAndCprErrors) {
    cdftest::Minimal m = cdftest::minimal_real8("x", {1.0, 2.0});
    auto vxr = [&] { (void)det::parse_vxr(view(m), m.vxr); };

    m.blob.put_i32(m.vxr_nentries(), 0);
    EXPECT_EQ(code_of(vxr), cdf::ErrorCode::BadRecordSize);
    m.blob.put_i32(m.vxr_nentries(), 100000);
    EXPECT_EQ(code_of(vxr), cdf::ErrorCode::BadRecordSize);
    m.blob.put_i32(m.vxr_nentries(), 2);                   // more entries than the 44-byte record holds
    EXPECT_EQ(code_of(vxr), cdf::ErrorCode::BadRecordSize);
    m.blob.put_i32(m.vxr_nentries(), 1);
    m.blob.put_i32(m.vxr_nused(), 2);                      // used > capacity
    EXPECT_EQ(code_of(vxr), cdf::ErrorCode::BadRecordSize);
    m.blob.put_i32(m.vxr_nused(), -1);
    EXPECT_EQ(code_of(vxr), cdf::ErrorCode::BadRecordSize);
    m.blob.put_i32(m.vxr_nused(), 1);
    m.blob.put_i64(m.vxr, 20);
    EXPECT_EQ(code_of(vxr), cdf::ErrorCode::BadRecordSize);
    m.blob.put_i64(m.vxr, 44);
    EXPECT_EQ(code_of(vxr), cdf::ErrorCode::None);

    // A CVVR: reinterpret the VVR as one.
    m.blob.put_i32(m.vvr + 8, static_cast<std::int32_t>(cdf::RecordType::Cvvr));
    m.blob.put_i64(m.vvr + 16, 4);                        // cSize
    auto cvvr = [&] { (void)det::parse_cvvr(view(m), m.vvr); };
    EXPECT_EQ(code_of(cvvr), cdf::ErrorCode::None);
    EXPECT_EQ(det::parse_cvvr(view(m), m.vvr).payload.size(), 4u);
    m.blob.put_i64(m.vvr + 16, 1000);                     // payload past the record
    EXPECT_EQ(code_of(cvvr), cdf::ErrorCode::BadRecordSize);
    m.blob.put_i64(m.vvr + 16, -1);
    EXPECT_EQ(code_of(cvvr), cdf::ErrorCode::BadRecordSize);
    m.blob.put_i64(m.vvr, 20);                            // record too short for the CVVR header
    EXPECT_EQ(code_of(cvvr), cdf::ErrorCode::BadRecordSize);
    m.blob.put_i64(m.vvr, 28);
    m.blob.put_i32(m.vvr + 8, static_cast<std::int32_t>(cdf::RecordType::Vvr));
    EXPECT_EQ(code_of(cvvr), cdf::ErrorCode::UnexpectedRecordType);

    // A CPR: reinterpret the VVR as one (28 bytes: 24 fixed + 1 param).
    m.blob.put_i32(m.vvr + 8, static_cast<std::int32_t>(cdf::RecordType::Cpr));
    m.blob.put_i32(m.vvr + 12, static_cast<std::int32_t>(cdf::CompressionType::Gzip));
    m.blob.put_i32(m.vvr + 20, 1);
    m.blob.put_i32(m.vvr + 24, 6);
    auto cpr = [&] { (void)det::parse_cpr(view(m), m.vvr); };
    EXPECT_EQ(code_of(cpr), cdf::ErrorCode::None);
    EXPECT_EQ(det::parse_cpr(view(m), m.vvr).level, 6);
    m.blob.put_i32(m.vvr + 20, 0);
    EXPECT_EQ(det::parse_cpr(view(m), m.vvr).level, 0);   // no params => level 0
    m.blob.put_i32(m.vvr + 20, 2);                        // 2 params need 32 bytes
    EXPECT_EQ(code_of(cpr), cdf::ErrorCode::BadRecordSize);
    m.blob.put_i32(m.vvr + 20, 99);
    EXPECT_EQ(code_of(cpr), cdf::ErrorCode::BadRecordSize);
    m.blob.put_i32(m.vvr + 20, 1);
    m.blob.put_i32(m.vvr + 12, 4);                        // the unassigned code
    EXPECT_EQ(code_of(cpr), cdf::ErrorCode::UnsupportedCompression);
    m.blob.put_i32(m.vvr + 12, static_cast<std::int32_t>(cdf::CompressionType::Gzip));
    m.blob.put_i64(m.vvr, 20);
    EXPECT_EQ(code_of(cpr), cdf::ErrorCode::BadRecordSize);
}

// ---- the index walk's caps actually fire -----------------------------------------------------------

TEST(CdfIndex, CyclesAreCaughtByTheCaps) {
    cdftest::Minimal m = cdftest::minimal_real8("x", {1.0, 2.0});
    // A VXR whose only entry points back at itself: depth grows without bound.
    m.blob.put_i64(m.vxr_child0(), static_cast<std::int64_t>(m.vxr));
    EXPECT_EQ(code_of([&] { det::build_index(view(m), m.vxr); }), cdf::ErrorCode::VxrTreeTooDeep);

    // A VXR whose `next` is itself: same depth forever, so the node cap is what stops it.
    m.blob.put_i64(m.vxr_child0(), static_cast<std::int64_t>(m.vvr));
    m.blob.put_i64(m.vxr_next(), static_cast<std::int64_t>(m.vxr));
    EXPECT_EQ(code_of([&] { det::build_index(view(m), m.vxr); }), cdf::ErrorCode::VxrTreeTooLarge);
    m.blob.put_i64(m.vxr_next(), 0);

    // An entry pointing at a record that is neither VXR, VVR nor CVVR.
    m.blob.put_i64(m.vxr_child0(), static_cast<std::int64_t>(m.gdr));
    EXPECT_EQ(code_of([&] { det::build_index(view(m), m.vxr); }), cdf::ErrorCode::UnexpectedRecordType);
    // An entry pointing past the end of the file.
    m.blob.put_i64(m.vxr_child0(), static_cast<std::int64_t>(m.blob.size() + 8));
    EXPECT_EQ(code_of([&] { det::build_index(view(m), m.vxr); }), cdf::ErrorCode::BadOffset);
    m.blob.put_i64(m.vxr_child0(), static_cast<std::int64_t>(m.vvr));
    // An entry whose record range is inverted or negative.
    m.blob.put_i32(m.vxr_first0(), 5);
    m.blob.put_i32(m.vxr_last0(), 2);
    EXPECT_EQ(code_of([&] { det::build_index(view(m), m.vxr); }), cdf::ErrorCode::RecordOutOfRange);
    m.blob.put_i32(m.vxr_first0(), -1);
    EXPECT_EQ(code_of([&] { det::build_index(view(m), m.vxr); }), cdf::ErrorCode::RecordOutOfRange);
    m.blob.put_i32(m.vxr_first0(), 0);
    m.blob.put_i32(m.vxr_last0(), 1);
    EXPECT_EQ(det::build_index(view(m), m.vxr).size(), 1u);

    // A VXR with zero used entries and no next contributes nothing.
    m.blob.put_i32(m.vxr_nused(), 0);
    EXPECT_TRUE(det::build_index(view(m), m.vxr).empty());
}

TEST(CdfIndex, ExtentsAreSortedRegardlessOfFileOrder) {
    // Two VXR entries listed in descending record order must come back ascending.
    cdftest::Minimal m = cdftest::minimal_real8("x", {1.0, 2.0, 3.0, 4.0});
    // Grow the VXR to 2 entries: 28 + 32 = 60 bytes. Insert 16 bytes after it, shifting the VVR.
    std::vector<std::byte>& bytes = m.blob.bytes;
    bytes.insert(bytes.begin() + static_cast<std::ptrdiff_t>(m.vxr + 44), 16, std::byte{0});
    const std::uint64_t vvr = m.vvr + 16;
    m.blob.put_i64(m.vxr, 60);
    m.blob.put_i32(m.vxr + 20, 2);
    m.blob.put_i32(m.vxr + 24, 2);
    // entry 0 covers 2..3 (later records), entry 1 covers 0..1 — both point at the same VVR;
    // the offsets are what is being sorted, the payload is not read here.
    m.blob.put_i32(m.vxr + 28, 2); m.blob.put_i32(m.vxr + 32, 0);          // First[0], First[1]
    m.blob.put_i32(m.vxr + 36, 3); m.blob.put_i32(m.vxr + 40, 1);          // Last[0], Last[1]
    m.blob.put_i64(m.vxr + 44, static_cast<std::int64_t>(vvr));            // Offset[0]
    m.blob.put_i64(m.vxr + 52, static_cast<std::int64_t>(vvr));            // Offset[1]
    m.blob.put_i64(m.gdr_eof(), static_cast<std::int64_t>(m.blob.size()));
    const auto ext = det::build_index(view(m), m.vxr);
    ASSERT_EQ(ext.size(), 2u);
    EXPECT_EQ(ext[0].first, 0);
    EXPECT_EQ(ext[1].first, 2);
    EXPECT_EQ(det::find_extent(ext, 1), ext.data());
    EXPECT_EQ(det::find_extent(ext, 2), &ext[1]);
}
