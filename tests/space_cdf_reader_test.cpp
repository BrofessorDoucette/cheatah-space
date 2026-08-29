// Unit tests for space/cdf/reader.hpp and encoding.hpp — opening a CDF and reading variables
// into ndarrays.
//
// The real-file values are the ones decoded by hand from the raw bytes before any of this code
// existed: OMNI's `F` opens with 6.92, 5.84, 5.71 nT and its `Epoch` with 63587289600000.0 ms.
// The builder-based tests cover every storage type in both byte orders offline, so the decode
// kernels are fully exercised even on a checkout with no corpus. The final test flips every byte
// of a small file one at a time, through open_from_memory so ASan sees each read: the reader
// must either succeed or throw CdfError, and never do anything else.
#include <gtest/gtest.h>

#include <cmath>
#include <cstring>
#include <fstream>
#include <functional>
#include <string>
#include <vector>

#include "cdf_builder.hpp"
#include "space/space.hpp"

namespace cdf = cheatah::space::cdf;
namespace det = cheatah::space::cdf::detail;
namespace nd = cheatah::ndarray;

namespace {

// Flat element access. `a[i]` needs as many indices as the array has dimensions, and most CDF
// variables are shaped [records, dims...]; the reader fills the buffer in flat row-major order,
// which is exactly what these tests want to compare.
template <class T>
const T* flat(const cheatah::ndarray::basic_ndarray<T>& a) {
    return a.buffer()->data() + a.offset();
}

std::string corpus_path(const char* rel) {
    std::string prefix;
    for (int up = 0; up < 5; ++up) {
        const std::string candidate = prefix + "space/cdf/vendor/corpus/" + rel;
        if (std::ifstream(candidate, std::ios::binary).good()) { return candidate; }
        prefix += "../";
    }
    return {};
}

cdf::ErrorCode code_of(const std::function<void()>& f) {
    try { f(); } catch (const cdf::CdfError& e) { return e.code(); }
    return cdf::ErrorCode::None;
}

cdf::File open_built(const cdftest::Minimal& m) {
    return cdf::open_from_memory(m.blob.bytes.data(), m.blob.size());
}

}  // namespace

// ---- a real file: OMNI, the hand-decoded numbers ------------------------------------------------

TEST(CdfReader, OmniOpensAndDecodesTheHandVerifiedValues) {
    const std::string path = corpus_path("tier1/omni_hro2_1min_20150101_v01.cdf");
    if (path.empty()) { GTEST_SKIP() << "run: scripts/cdf-corpus.sh fetch --tier 1"; }
    const cdf::File f = cdf::open(path);

    const std::vector<std::string> names = cdf::var_names(f);
    ASSERT_EQ(names.size(), 47u);
    EXPECT_EQ(names[0], "Epoch");
    EXPECT_EQ(cdf::record_count(f, "F"), 44640);
    EXPECT_EQ(cdf::data_type_name(f, "F"), "CDF_REAL4");
    EXPECT_EQ(cdf::data_type_name(f, "Epoch"), "CDF_EPOCH");
    EXPECT_EQ(cdf::shape(f, "F"), (std::vector<long long>{44640}));

    const nd::basic_ndarray<double> b = cdf::values(f, "F");
    ASSERT_EQ(b.size(), 44640u);
    EXPECT_NEAR(flat(b)[0], 6.92, 0.005);
    EXPECT_NEAR(flat(b)[1], 5.84, 0.005);
    EXPECT_NEAR(flat(b)[2], 5.71, 0.005);

    const nd::basic_ndarray<double> t = cdf::values(f, "Epoch");
    ASSERT_EQ(t.size(), 44640u);
    EXPECT_EQ(flat(t)[0], 63587289600000.0);   // 2015-01-01T00:00:00Z as CDF_EPOCH ms
    EXPECT_EQ(flat(t)[1] - flat(t)[0], 60000.0);              // one-minute cadence
    EXPECT_EQ(flat(t)[44639] - flat(t)[0], 44639 * 60000.0);

    EXPECT_EQ(code_of([&] { (void)cdf::values(f, "no_such_variable"); }), cdf::ErrorCode::VariableNotFound);

    // A copy shares the mapping and reads the same bytes.
    const cdf::File g = f;  // NOLINT(performance-unnecessary-copy-initialization) — the copy IS the test
    const nd::basic_ndarray<double> gb = cdf::values(g, "F");
    EXPECT_EQ(flat(gb)[0], flat(b)[0]);
}

// ---- a real file: every type test_alltypes has ---------------------------------------------------

TEST(CdfReader, TestAlltypesReadsEveryUncompressedVariable) {
    const std::string path = corpus_path("tier0/test_alltypes.cdf");
    if (path.empty()) { GTEST_SKIP() << "run: scripts/cdf-corpus.sh fetch"; }
    const cdf::File f = cdf::open(path);
    const det::FileImpl& s = f.state();
    ASSERT_EQ(cdf::var_names(f).size(), 21u);

    int read_as_float = 0, read_as_int = 0, refused_lossy = 0, compressed = 0, sparse = 0;
    for (const det::Vdr& v : s.vars) {
        const std::string name(v.name);
        const std::vector<long long> shp = cdf::shape(f, name);
        ASSERT_FALSE(shp.empty());
        EXPECT_EQ(shp[0], cdf::record_count(f, name));
        if (v.s_records != cdf::SparseRecords::None) {
            // `Temp` is srecords.PAD with real gaps in the index. Filling them needs the pad
            // value and the previous-record rule, which is deferred — so the reader must refuse
            // rather than hand back whatever bytes happen to sit in the allocation slack.
            EXPECT_EQ(code_of([&] { (void)cdf::values(f, name); }), cdf::ErrorCode::UnsupportedSparse) << name;
            ++sparse;
            continue;
        }
        if (v.compressed()) {
            // GZIP lands in M8; until then the reader must refuse, not guess.
            EXPECT_EQ(code_of([&] { (void)cdf::values(f, name); }), cdf::ErrorCode::UnsupportedCompression) << name;
            ++compressed;
            continue;
        }
        const det::StorageKind kind = det::storage_kind(v.data_type);
        if (kind == det::StorageKind::I64) {
            EXPECT_EQ(code_of([&] { (void)cdf::values(f, name); }), cdf::ErrorCode::LossyConversion) << name;
            auto a = cdf::values_i64(f, name);
            EXPECT_EQ(static_cast<long long>(a.size()), cdf::record_count(f, name) * static_cast<long long>(v.values_per_record()));
            ++read_as_int;
        } else if (det::is_float_kind(kind)) {
            EXPECT_EQ(code_of([&] { (void)cdf::values_i64(f, name); }), cdf::ErrorCode::LossyConversion) << name;
            auto a = cdf::values(f, name);
            std::size_t expect = 1;
            for (long long d : shp) { expect *= static_cast<std::size_t>(d); }
            EXPECT_EQ(a.size(), expect) << name;
            ++read_as_float;
            ++refused_lossy;
        } else {
            // Every integer type reads both ways, and the two agree value for value.
            const auto d = cdf::values(f, name);
            const auto i = cdf::values_i64(f, name);
            ASSERT_EQ(d.size(), i.size()) << name;
            for (std::size_t k = 0; k < d.size(); ++k) {
                EXPECT_EQ(flat(d)[k], static_cast<double>(flat(i)[k])) << name << " [" << k << "]";
            }
            ++read_as_float;
            ++read_as_int;
        }
    }
    EXPECT_EQ(compressed, 2);            // Longitude and longitude_dup
    EXPECT_EQ(sparse, 1);                // Temp
    EXPECT_GT(read_as_float, 10);
    EXPECT_GT(read_as_int, 5);
    EXPECT_GT(refused_lossy, 3);

    // Shapes carry dimensions, string lengths and the EPOCH16 pair.
    EXPECT_EQ(cdf::shape(f, "volume"), (std::vector<long long>{1, 2, 4, 2}));
    EXPECT_EQ(cdf::shape(f, "Name"), (std::vector<long long>{2, 2, 10}));
    bool saw_epoch16 = false;
    for (const std::string& n : cdf::var_names(f)) {
        if (cdf::data_type_name(f, n) == "CDF_EPOCH16") {
            saw_epoch16 = true;
            EXPECT_EQ(cdf::shape(f, n).back(), 2);
            const auto e = cdf::values(f, n);
            EXPECT_EQ(e.size() % 2, 0u);
        }
    }
    EXPECT_TRUE(saw_epoch16);
    // A CHAR variable comes back as bytes: printable ASCII in this file.
    const auto chars = cdf::values(f, "Name");
    EXPECT_EQ(chars.size(), 40u);
    EXPECT_GE(flat(chars)[0], 32.0);
    EXPECT_LE(flat(chars)[0], 126.0);
}

// ---- crafted files: every storage type in both byte orders ------------------------------------------

namespace {
template <class Raw>
void check_type(cdf::DataType type, const std::vector<Raw>& vals, cdf::Encoding enc) {
    std::vector<std::byte> raw;
    for (Raw v : vals) { cdftest::push_value(raw, v, enc == cdf::Encoding::Network); }
    const cdftest::Minimal m = cdftest::minimal_raw("v", type, 1, static_cast<std::int32_t>(vals.size()), raw, enc);
    const cdf::File f = open_built(m);
    EXPECT_EQ(cdf::record_count(f, "v"), static_cast<long long>(vals.size()));
    if constexpr (std::is_floating_point_v<Raw>) {
        auto d = cdf::values(f, "v");
        ASSERT_EQ(d.size(), vals.size());
        for (std::size_t i = 0; i < vals.size(); ++i) { EXPECT_EQ(d[i], static_cast<double>(vals[i])); }
        EXPECT_EQ(code_of([&] { (void)cdf::values_i64(f, "v"); }), cdf::ErrorCode::LossyConversion);
    } else {
        const auto i64 = cdf::values_i64(f, "v");
        ASSERT_EQ(i64.size(), vals.size());
        for (std::size_t i = 0; i < vals.size(); ++i) { EXPECT_EQ(flat(i64)[i], static_cast<long long>(vals[i])); }
        if constexpr (sizeof(Raw) == 8) {
            EXPECT_EQ(code_of([&] { (void)cdf::values(f, "v"); }), cdf::ErrorCode::LossyConversion);
        } else {
            const auto d = cdf::values(f, "v");
            for (std::size_t i = 0; i < vals.size(); ++i) { EXPECT_EQ(flat(d)[i], static_cast<double>(vals[i])); }
        }
    }
}
}  // namespace

TEST(CdfReader, EveryStorageTypeDecodesInBothByteOrders) {
    for (cdf::Encoding enc : {cdf::Encoding::Network, cdf::Encoding::IbmPc}) {
        check_type<std::int8_t>(cdf::DataType::Int1, {-128, -1, 0, 1, 127}, enc);
        check_type<std::int8_t>(cdf::DataType::Byte, {-5, 5}, enc);
        check_type<std::int16_t>(cdf::DataType::Int2, {-32768, -2, 0, 2, 32767}, enc);
        check_type<std::int32_t>(cdf::DataType::Int4, {-2147483647 - 1, -3, 0, 3, 2147483647}, enc);
        check_type<std::int64_t>(cdf::DataType::Int8, {-9007199254740993LL, 0, 9007199254740993LL}, enc);
        check_type<std::int64_t>(cdf::DataType::TimeTt2000, {473385668184000000LL, 0}, enc);   // 2015-01-01
        check_type<std::uint8_t>(cdf::DataType::Uint1, {0, 1, 255}, enc);
        check_type<std::uint8_t>(cdf::DataType::Uchar, {65, 66}, enc);
        check_type<std::uint16_t>(cdf::DataType::Uint2, {0, 1, 65535}, enc);
        check_type<std::uint32_t>(cdf::DataType::Uint4, {0, 1, 4294967295u}, enc);
        check_type<float>(cdf::DataType::Real4, {-1.5f, 0.0f, 6.92f, 1e30f}, enc);
        check_type<float>(cdf::DataType::Float, {2.5f}, enc);
        check_type<double>(cdf::DataType::Real8, {-1.0e31, 0.0, 63587289600000.0}, enc);
        check_type<double>(cdf::DataType::Double, {3.25}, enc);
        check_type<double>(cdf::DataType::Epoch, {63587289600000.0, 63587289660000.0}, enc);
    }
    // -0.0 and NaN survive the byte path bit-for-bit — a decoder that went through text would not.
    std::vector<std::byte> raw;
    cdftest::push_value(raw, -0.0, true);
    cdftest::push_value(raw, std::nan(""), true);
    const cdftest::Minimal m = cdftest::minimal_raw("v", cdf::DataType::Real8, 1, 2, raw);
    const auto d = cdf::values(open_built(m), "v");
    EXPECT_TRUE(std::signbit(flat(d)[0]));
    EXPECT_TRUE(std::isnan(flat(d)[1]));
}

TEST(CdfReader, CharStringsAndEpoch16PairsShapeCorrectly) {
    // Two records of a CHAR[4] string: shape [2, 4], byte values.
    std::vector<std::byte> raw;
    for (char c : std::string("abcdWXYZ")) { raw.push_back(std::byte{static_cast<unsigned char>(c)}); }
    const cdftest::Minimal m = cdftest::minimal_raw("s", cdf::DataType::Char, 4, 2, raw);
    const cdf::File f = open_built(m);
    EXPECT_EQ(cdf::shape(f, "s"), (std::vector<long long>{2, 4}));
    const auto d = cdf::values(f, "s");
    ASSERT_EQ(d.size(), 8u);
    EXPECT_EQ(flat(d)[0], 'a');
    EXPECT_EQ(flat(d)[7], 'Z');
    const auto i = cdf::values_i64(f, "s");
    EXPECT_EQ(flat(i)[4], 'W');

    // EPOCH16: one record is two doubles, shape [1, 2].
    std::vector<std::byte> e16;
    cdftest::push_value(e16, 1.0, true);
    cdftest::push_value(e16, 2.0, true);
    const cdftest::Minimal m2 = cdftest::minimal_raw("e", cdf::DataType::Epoch16, 1, 1, e16);
    const cdf::File g = open_built(m2);
    EXPECT_EQ(cdf::shape(g, "e"), (std::vector<long long>{1, 2}));
    const auto v = cdf::values(g, "e");
    ASSERT_EQ(v.size(), 2u);
    EXPECT_EQ(flat(v)[0], 1.0);
    EXPECT_EQ(flat(v)[1], 2.0);
}

TEST(CdfReader, EmptyVariableAndEmptyFileHandle) {
    const cdftest::Minimal m = cdftest::minimal_real8("none", {});
    const cdf::File f = open_built(m);
    EXPECT_EQ(cdf::record_count(f, "none"), 0);
    EXPECT_EQ(cdf::values(f, "none").size(), 0u);
    EXPECT_EQ(cdf::shape(f, "none"), (std::vector<long long>{0}));

    const cdf::File empty;
    EXPECT_EQ(code_of([&] { (void)cdf::var_names(empty); }), cdf::ErrorCode::EmptyFile);
    EXPECT_EQ(code_of([&] { cdf::open("/nonexistent/x.cdf"); }), cdf::ErrorCode::CannotOpen);
    const std::byte tiny[3] = {};
    EXPECT_EQ(code_of([&] { cdf::open_from_memory(tiny, 3); }), cdf::ErrorCode::EmptyFile);
}

// ---- crafted files: every refusal --------------------------------------------------------------------

TEST(CdfReader, RefusesWhatItDoesNotReadYetByName) {
    cdftest::Minimal m = cdftest::minimal_real8("x", {1.0, 2.0});
    auto open = [&] { (void)open_built(m); };

    m.blob.put_be32(4, cdf::kMagicCompressed);                            // whole-file CCR
    EXPECT_EQ(code_of(open), cdf::ErrorCode::UnsupportedCompression);
    m.blob.put_be32(4, cdf::kMagicUncompressed);

    m.blob.put_be32(m.cdr_flags(), cdf::kCdrFlagRowMajority);              // single-file bit cleared
    EXPECT_EQ(code_of(open), cdf::ErrorCode::UnsupportedMultiFile);
    m.blob.put_be32(m.cdr_flags(), cdf::kCdrFlagRowMajority | cdf::kCdrFlagSingleFile);

    for (cdf::Encoding vax : {cdf::Encoding::Vax, cdf::Encoding::AlphaVmsD, cdf::Encoding::AlphaVmsG}) {
        m.blob.put_i32(m.cdr_encoding(), static_cast<std::int32_t>(vax));
        EXPECT_EQ(code_of(open), cdf::ErrorCode::UnsupportedEncoding);
    }
    m.blob.put_i32(m.cdr_encoding(), static_cast<std::int32_t>(cdf::Encoding::Network));
    EXPECT_EQ(code_of(open), cdf::ErrorCode::None);

    // The rVDR chain pointing at the zVDR: the same name arrives twice.
    m.blob.put_i64(m.gdr + 12, static_cast<std::int64_t>(m.zvdr));
    m.blob.put_i32(m.gdr + 44, 1);
    EXPECT_EQ(code_of(open), cdf::ErrorCode::DuplicateVariable);
    m.blob.put_i64(m.gdr + 12, 0);
    m.blob.put_i32(m.gdr + 44, 0);

    // A zVDR chain that loops: more descriptors than the GDR declares.
    m.blob.put_i64(m.zvdr_next(), static_cast<std::int64_t>(m.zvdr));
    EXPECT_EQ(code_of(open), cdf::ErrorCode::BadDimensions);
    m.blob.put_i64(m.zvdr_next(), 0);

    // A gap: the only extent starts at record 1.
    m.blob.put_i32(m.vxr_first0(), 1);
    EXPECT_EQ(code_of([&] { (void)cdf::values(open_built(m), "x"); }), cdf::ErrorCode::RecordOutOfRange);
    m.blob.put_i32(m.zvdr_s_records(), static_cast<std::int32_t>(cdf::SparseRecords::Pad));
    EXPECT_EQ(code_of([&] { (void)cdf::values(open_built(m), "x"); }), cdf::ErrorCode::UnsupportedSparse);
    m.blob.put_i32(m.zvdr_s_records(), 0);
    m.blob.put_i32(m.vxr_first0(), 0);

    // An extent that ends early: records 0..0 only, but two records exist.
    m.blob.put_i32(m.vxr_last0(), 0);
    EXPECT_EQ(code_of([&] { (void)cdf::values(open_built(m), "x"); }), cdf::ErrorCode::RecordOutOfRange);
    m.blob.put_i32(m.vxr_last0(), 1);

    // A VVR too short for the records the index says it holds.
    m.blob.put_i64(m.vvr, 12 + 8);
    EXPECT_EQ(code_of([&] { (void)cdf::values(open_built(m), "x"); }), cdf::ErrorCode::BadRecordSize);
    m.blob.put_i64(m.vvr, 12 + 16);

    // The extent points at a CVVR: compression is M8's job, so refuse for now.
    m.blob.put_i32(m.vvr + 8, static_cast<std::int32_t>(cdf::RecordType::Cvvr));
    m.blob.put_i64(m.vvr + 16, 4);
    EXPECT_EQ(code_of([&] { (void)cdf::values(open_built(m), "x"); }), cdf::ErrorCode::UnsupportedCompression);
    m.blob.put_i32(m.vvr + 8, static_cast<std::int32_t>(cdf::RecordType::Vvr));

    // An extent past the record count is ignored, not an error: allocation slack is normal.
    m.blob.put_i32(m.zvdr_max_rec(), 0);
    EXPECT_EQ(cdf::values(open_built(m), "x").size(), 1u);
    m.blob.put_i32(m.zvdr_max_rec(), 1);
    EXPECT_EQ(cdf::values(open_built(m), "x").size(), 2u);
}

TEST(CdfReader, ColumnMajorMultiDimensionalIsRefusedButOneDIsFine) {
    // Grow the zVDR to carry two dimensions [2,1], both varying, and clear the row-major flag.
    cdftest::Minimal m = cdftest::minimal_real8("d", {1.0, 2.0});
    std::vector<std::byte>& bytes = m.blob.bytes;
    bytes.insert(bytes.begin() + static_cast<std::ptrdiff_t>(m.zvdr + 344), 16, std::byte{0});
    m.blob.put_i64(m.zvdr, 360);
    m.blob.put_i32(m.zvdr_num_dims(), 2);
    m.blob.put_i32(m.zvdr + 344, 2); m.blob.put_i32(m.zvdr + 348, 1);      // sizes
    m.blob.put_i32(m.zvdr + 352, -1); m.blob.put_i32(m.zvdr + 356, -1);    // varys
    m.blob.put_i32(m.zvdr_max_rec(), 0);                                    // one record of 2x1
    m.blob.put_i64(m.zvdr_vxr_head(), static_cast<std::int64_t>(m.vxr + 16));
    m.blob.put_i64(m.zvdr + 36, static_cast<std::int64_t>(m.vxr + 16));
    m.blob.put_i32(m.vxr + 16 + 32, 0);                                     // Last[0]
    m.blob.put_i64(m.vxr + 16 + 36, static_cast<std::int64_t>(m.vvr + 16));
    m.blob.put_i64(m.gdr_eof(), static_cast<std::int64_t>(m.blob.size()));

    EXPECT_EQ(cdf::shape(open_built(m), "d"), (std::vector<long long>{1, 2, 1}));
    EXPECT_EQ(cdf::values(open_built(m), "d").size(), 2u);                 // row-major: fine

    m.blob.put_be32(m.cdr_flags(), cdf::kCdrFlagSingleFile);               // column-major
    EXPECT_EQ(code_of([&] { (void)cdf::values(open_built(m), "d"); }), cdf::ErrorCode::UnsupportedLayout);

    // One dimension is the same in either majority, so it is not refused.
    m.blob.put_i32(m.zvdr_num_dims(), 1);
    m.blob.put_i32(m.zvdr + 348, -1);                                       // varys[0] now at +348
    EXPECT_EQ(cdf::values(open_built(m), "d").size(), 2u);
}

// ---- hostile input: every single-byte corruption of a whole file ---------------------------------------

TEST(CdfReader, EverySingleByteFlipEitherReadsOrThrowsCdfError) {
    const cdftest::Minimal m = cdftest::minimal_real8("x", {1.0, 2.0, 3.0});
    const std::vector<std::byte> good = m.blob.bytes;
    int ok = 0, refused = 0;
    for (std::size_t i = 0; i < good.size(); ++i) {
        for (unsigned char flip : {0xFFu, 0x80u, 0x01u}) {
            std::vector<std::byte> bad = good;
            bad[i] = std::byte{static_cast<unsigned char>(static_cast<unsigned char>(bad[i]) ^ flip)};
            try {
                const cdf::File f = cdf::open_from_memory(bad.data(), bad.size());
                for (const std::string& n : cdf::var_names(f)) {
                    (void)cdf::shape(f, n);
                    (void)cdf::data_type_name(f, n);
                    // A read may legitimately refuse (a flipped type byte, a bent offset); it
                    // may not do anything else. Counting rather than swallowing keeps that a
                    // measurement instead of a silence.
                    if (code_of([&] { (void)cdf::values(f, n); }) != cdf::ErrorCode::None) { ++refused; }
                    if (code_of([&] { (void)cdf::values_i64(f, n); }) != cdf::ErrorCode::None) { ++refused; }
                }
                ++ok;
            } catch (const cdf::CdfError&) {
                ++refused;
            }
            // Anything else — a std::exception, a crash, an ASan report — fails the test.
        }
    }
    EXPECT_GT(ok, 0);
    EXPECT_GT(refused, 0);
    // Truncation at every length, too. A short file must be refused or read what survives —
    // never over-read, which is what running this through open_from_memory lets ASan prove.
    int truncated_ok = 0;
    for (std::size_t len = 0; len < good.size(); ++len) {
        const cdf::ErrorCode c = code_of([&] {
            const cdf::File f = cdf::open_from_memory(good.data(), len);
            for (const std::string& n : cdf::var_names(f)) {
                (void)code_of([&] { (void)cdf::values(f, n); });
            }
        });
        if (c == cdf::ErrorCode::None) { ++truncated_ok; }
    }
    EXPECT_LT(truncated_ok, static_cast<int>(good.size()));  // most truncations must be refused
}

// ---- the decode layer's own edges ------------------------------------------------------------------
//
// reader.hpp validates the data type when it parses the VDR, so a bogus type can never reach
// decode_run() through the public API. The layer is still tested directly: a switch arm that no
// caller can reach is exactly the kind of thing that rots, and the alternative — deleting the
// arm — would cost the compiler warning that catches a newly added enumerator.
TEST(CdfEncoding, StorageKindAndTheUnknownArms) {
    // Every defined type maps to a real layout whose size matches the format's element size,
    // except EPOCH16, which is deliberately two F64 elements per value.
    struct Case { cdf::DataType type; det::StorageKind kind; };
    constexpr Case kCases[] = {
        {cdf::DataType::Int1, det::StorageKind::I8},   {cdf::DataType::Byte, det::StorageKind::I8},
        {cdf::DataType::Int2, det::StorageKind::I16},  {cdf::DataType::Int4, det::StorageKind::I32},
        {cdf::DataType::Int8, det::StorageKind::I64},  {cdf::DataType::TimeTt2000, det::StorageKind::I64},
        {cdf::DataType::Uint1, det::StorageKind::U8},  {cdf::DataType::Char, det::StorageKind::U8},
        {cdf::DataType::Uchar, det::StorageKind::U8},  {cdf::DataType::Uint2, det::StorageKind::U16},
        {cdf::DataType::Uint4, det::StorageKind::U32}, {cdf::DataType::Real4, det::StorageKind::F32},
        {cdf::DataType::Float, det::StorageKind::F32}, {cdf::DataType::Real8, det::StorageKind::F64},
        {cdf::DataType::Double, det::StorageKind::F64},{cdf::DataType::Epoch, det::StorageKind::F64},
        {cdf::DataType::Epoch16, det::StorageKind::F64},
    };
    for (const Case& c : kCases) {
        EXPECT_EQ(det::storage_kind(c.type), c.kind) << det::data_type_name(c.type);
        EXPECT_EQ(det::raw_element_bytes(c.type), det::storage_bytes(c.kind));
        const bool is_float = det::is_float_kind(c.kind);
        EXPECT_EQ(is_float, c.kind == det::StorageKind::F32 || c.kind == det::StorageKind::F64);
        if (c.type != cdf::DataType::Epoch16) {
            EXPECT_EQ(det::raw_element_bytes(c.type), static_cast<std::uint64_t>(det::element_size(c.type)));
        } else {
            EXPECT_EQ(det::raw_element_bytes(c.type), 8u);   // two of these per EPOCH16 value
            EXPECT_EQ(det::element_size(c.type), 16);
        }
    }

    // A type the format does not define maps to Unknown, sizes to 0, and decodes to an error
    // rather than reading whatever the switch happened to fall through to.
    const auto bogus = static_cast<cdf::DataType>(4242);
    EXPECT_EQ(det::storage_kind(bogus), det::StorageKind::Unknown);
    EXPECT_EQ(det::storage_bytes(det::StorageKind::Unknown), 0u);
    EXPECT_EQ(det::raw_element_bytes(bogus), 0u);
    EXPECT_FALSE(det::is_float_kind(det::StorageKind::Unknown));
    EXPECT_EQ(det::storage_bytes(static_cast<det::StorageKind>(200)), 0u);

    const std::byte src[8] = {};
    double d[1]{};
    long long i[1]{};
    EXPECT_EQ(code_of([&] { det::decode_run<double>(src, bogus, false, 1, d, 0); }),
              cdf::ErrorCode::BadDataType);
    EXPECT_EQ(code_of([&] { det::decode_run<long long>(src, bogus, false, 1, i, 0); }),
              cdf::ErrorCode::BadDataType);
    // And the lossy pairings, refused at this layer too — not only at the reader's front door.
    EXPECT_EQ(code_of([&] { det::decode_run<double>(src, cdf::DataType::Int8, false, 1, d, 0); }),
              cdf::ErrorCode::LossyConversion);
    EXPECT_EQ(code_of([&] { det::decode_run<long long>(src, cdf::DataType::Real8, false, 1, i, 0); }),
              cdf::ErrorCode::LossyConversion);
}
