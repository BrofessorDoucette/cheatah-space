// Unit tests for the space.cdf format layer — types.hpp, bytes.hpp, mapping.hpp.
//
// Three things are being established, and the third is the one that matters most.
//
// 1. The vocabulary is complete and total: every enumerator has a size, a name and a decode
//    class, and every "unknown" fallback is reachable, because a hostile file can put any 32-bit
//    integer in a type field.
//
// 2. The byte reader decodes a REAL NASA CDF correctly. The expected values below were obtained
//    by decoding the OMNI file by hand from its raw bytes before any of this code existed, so
//    this is a check against the file, not against our own parser's opinion of the file.
//
// 3. Every guard actually bites. A bounds check that has never been observed to reject anything
//    is not a bounds check. Reads past the end, offsets past the end and negative offsets each
//    get their own assertion, and the syscall failures — which no crafted input can provoke —
//    are driven through the SysOps seam.
#include <gtest/gtest.h>

#include <cstdio>
#include <cstring>
#include <fstream>
#include <string>
#include <vector>

#include "space/space.hpp"

namespace cdf = cheatah::space::cdf;
namespace det = cheatah::space::cdf::detail;

namespace {

// The file decoded by hand in the planning notes. Absent on a fresh checkout, so every test that
// needs it skips with the command that fetches it — the QA gate must stay green offline.
//
// The path is searched UPWARD from the working directory rather than taken as-is, because the
// suite is run from at least two places: ctest runs it from the repo root, and scripts/coverage.sh
// runs it from build/cov. A bare relative path silently resolves in one and not the other, which
// does not fail the tests — it SKIPS them, so the coverage run quietly stopped exercising the
// real-file paths at all. That is a much worse failure than a broken test.
const std::string& omni_path() {
    static const std::string path = [] {
        const char* rel = "space/cdf/vendor/corpus/tier1/omni_hro2_1min_20150101_v01.cdf";
        std::string prefix;
        for (int up = 0; up < 5; ++up) {
            std::string candidate = prefix + rel;
            if (std::ifstream(candidate, std::ios::binary).good()) { return candidate; }
            prefix += "../";
        }
        return std::string{};
    }();
    return path;
}

bool have_omni() { return !omni_path().empty(); }

std::vector<std::byte> read_all(const char* path) {
    std::ifstream in(path, std::ios::binary);
    std::vector<char> raw((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
    const auto* p = reinterpret_cast<const std::byte*>(raw.data());
    // The two-pointer braced init selects the iterator-range constructor, which COPIES the bytes
    // into the returned vector before `raw` dies — cppcheck reads it as returning a pointer into
    // the local and flags a dangling lifetime that is not there.
    // cppcheck-suppress returnDanglingLifetime
    return {p, p + raw.size()};
}

// Every syscall fails. Drives the error paths in BasicFileMapping that no input can reach —
// you cannot craft a CDF that makes the kernel refuse to mmap.
struct FailingSysOps {
    static int open_file(const std::string& /*path*/) { return -1; }
    static std::uint64_t file_size(int /*fd*/) { return 0; }
    static void* map_file(int /*fd*/, std::uint64_t /*size*/) { return nullptr; }
    static void unmap_file(void* /*addr*/, std::uint64_t /*size*/) {}
    static void close_file(int /*fd*/) {}
    static void advise_willneed(void* /*addr*/, std::uint64_t /*size*/) {}
};

// Opens, reports a plausible size, then fails to map — the branch that must still close the
// descriptor it already opened.
struct MapFailsSysOps {
    static int open_file(const std::string& /*path*/) { return 3; }
    static std::uint64_t file_size(int /*fd*/) { return 4096; }
    static void* map_file(int /*fd*/, std::uint64_t /*size*/) { return nullptr; }
    static void unmap_file(void* /*addr*/, std::uint64_t /*size*/) {}
    static void close_file(int /*fd*/) {}
    static void advise_willneed(void* /*addr*/, std::uint64_t /*size*/) {}
};

// Opens but reports an empty file — the EmptyFile branch, again with a descriptor to close.
struct EmptySysOps {
    static int open_file(const std::string& /*path*/) { return 3; }
    static std::uint64_t file_size(int /*fd*/) { return 0; }
    static void* map_file(int /*fd*/, std::uint64_t /*size*/) { return nullptr; }
    static void unmap_file(void* /*addr*/, std::uint64_t /*size*/) {}
    static void close_file(int /*fd*/) {}
    static void advise_willneed(void* /*addr*/, std::uint64_t /*size*/) {}
};

}  // namespace

// ---- types: the vocabulary is total ------------------------------------------------------------

static_assert(det::element_size(cdf::DataType::Real4) == 4);
static_assert(det::element_size(cdf::DataType::Epoch16) == 16);
static_assert(det::is_time_type(cdf::DataType::TimeTt2000));
static_assert(det::encoding_class_of(cdf::Encoding::Network) == cdf::EncodingClass::BigIeee);

TEST(CdfTypes, EverySizeMatchesTheFormat) {
    struct Case { cdf::DataType type; long long size; };
    constexpr Case kCases[] = {
        {cdf::DataType::Int1, 1},   {cdf::DataType::Int2, 2},    {cdf::DataType::Int4, 4},
        {cdf::DataType::Int8, 8},   {cdf::DataType::Uint1, 1},   {cdf::DataType::Uint2, 2},
        {cdf::DataType::Uint4, 4},  {cdf::DataType::Real4, 4},   {cdf::DataType::Real8, 8},
        {cdf::DataType::Epoch, 8},  {cdf::DataType::Epoch16, 16},
        {cdf::DataType::TimeTt2000, 8},
        {cdf::DataType::Byte, 1},   {cdf::DataType::Float, 4},   {cdf::DataType::Double, 8},
        {cdf::DataType::Char, 1},   {cdf::DataType::Uchar, 1},
    };
    for (const Case& c : kCases) {
        EXPECT_EQ(det::element_size(c.type), c.size) << det::data_type_name(c.type);
        EXPECT_TRUE(det::is_known_data_type(c.type)) << det::data_type_name(c.type);
        EXPECT_NE(det::data_type_name(c.type), "unknown");
    }
    // REAL4/FLOAT and REAL8/DOUBLE are distinct format codes with identical layouts. Both occur
    // in real files, so both must be handled rather than normalized away.
    EXPECT_EQ(det::element_size(cdf::DataType::Real4), det::element_size(cdf::DataType::Float));
    EXPECT_EQ(det::element_size(cdf::DataType::Real8), det::element_size(cdf::DataType::Double));
    EXPECT_NE(det::data_type_name(cdf::DataType::Real4), det::data_type_name(cdf::DataType::Float));
}

TEST(CdfTypes, UndefinedValuesAreRejectedNotGuessed) {
    // A corrupt file can put any int32 in a type field. element_size() returning 0 must be
    // distinguishable from a real type, which is what is_known_data_type() is for.
    // Casting an out-of-range value is the point of the test, not an oversight: a corrupt file
    // supplies exactly this, and the code must cope rather than read past a table.
    // NOLINTNEXTLINE(clang-analyzer-optin.core.EnumCastOutOfRange)
    const auto bogus = static_cast<cdf::DataType>(9999);
    EXPECT_EQ(det::element_size(bogus), 0);
    EXPECT_FALSE(det::is_known_data_type(bogus));
    EXPECT_EQ(det::data_type_name(bogus), "unknown");
    EXPECT_FALSE(det::is_char_type(bogus));
    EXPECT_FALSE(det::is_time_type(bogus));

    EXPECT_EQ(det::record_type_name(static_cast<cdf::RecordType>(77)), "unknown");
    EXPECT_EQ(det::compression_name(static_cast<cdf::CompressionType>(4)), "unknown");
    EXPECT_EQ(det::encoding_class_name(static_cast<cdf::EncodingClass>(200)), "unknown");
    EXPECT_EQ(det::error_message(static_cast<cdf::ErrorCode>(200)), "unknown error");
    // Compression type 4 is genuinely unassigned by the format, not an oversight here.
    constexpr cdf::CompressionType kComp[] = {
        cdf::CompressionType::None, cdf::CompressionType::Rle, cdf::CompressionType::Huff,
        cdf::CompressionType::Ahuff, cdf::CompressionType::Gzip,
    };
    for (cdf::CompressionType c : kComp) { EXPECT_NE(det::compression_name(c), "unknown"); }
    EXPECT_EQ(det::compression_name(cdf::CompressionType::Gzip), "gzip");
    EXPECT_EQ(det::compression_name(cdf::CompressionType::None), "none");

    // encoding_class_of()'s switch covers every enumerator, so its trailing return is reachable
    // only through a value the format does not define — which is exactly what a corrupt CDR
    // supplies. It must fall back rather than read uninitialized memory.
    EXPECT_EQ(det::encoding_class_of(static_cast<cdf::Encoding>(999)),
              cdf::EncodingClass::BigIeee);
}

TEST(CdfTypes, CharAndTimePredicates) {
    EXPECT_TRUE(det::is_char_type(cdf::DataType::Char));
    EXPECT_TRUE(det::is_char_type(cdf::DataType::Uchar));
    EXPECT_FALSE(det::is_char_type(cdf::DataType::Int1));  // same size, different meaning
    EXPECT_TRUE(det::is_time_type(cdf::DataType::Epoch));
    EXPECT_TRUE(det::is_time_type(cdf::DataType::Epoch16));
    EXPECT_TRUE(det::is_time_type(cdf::DataType::TimeTt2000));
    EXPECT_FALSE(det::is_time_type(cdf::DataType::Real8));  // same storage as Epoch
}

TEST(CdfTypes, EveryRecordTypeIsNamed) {
    constexpr cdf::RecordType kAll[] = {
        cdf::RecordType::Cdr, cdf::RecordType::Gdr, cdf::RecordType::RVariableDescriptor,
        cdf::RecordType::Adr, cdf::RecordType::AgrEdr, cdf::RecordType::Vxr,
        cdf::RecordType::Vvr, cdf::RecordType::ZVariableDescriptor, cdf::RecordType::AzEdr,
        cdf::RecordType::Ccr, cdf::RecordType::Cpr, cdf::RecordType::Spr,
        cdf::RecordType::Cvvr, cdf::RecordType::UnusedInternal,
    };
    for (cdf::RecordType t : kAll) { EXPECT_NE(det::record_type_name(t), "unknown"); }
    EXPECT_EQ(det::record_type_name(cdf::RecordType::ZVariableDescriptor), "zVDR");
    EXPECT_EQ(det::record_type_name(cdf::RecordType::Vxr), "VXR");
}

TEST(CdfTypes, AllTwentyOneEncodingsCollapseToFourClasses) {
    struct Case { cdf::Encoding enc; cdf::EncodingClass cls; };
    constexpr Case kCases[] = {
        {cdf::Encoding::Network, cdf::EncodingClass::BigIeee},
        {cdf::Encoding::Sun, cdf::EncodingClass::BigIeee},
        {cdf::Encoding::Sgi, cdf::EncodingClass::BigIeee},
        {cdf::Encoding::IbmRs, cdf::EncodingClass::BigIeee},
        {cdf::Encoding::Ppc, cdf::EncodingClass::BigIeee},
        {cdf::Encoding::Mac, cdf::EncodingClass::BigIeee},
        {cdf::Encoding::Hp, cdf::EncodingClass::BigIeee},
        {cdf::Encoding::NeXT, cdf::EncodingClass::BigIeee},
        {cdf::Encoding::ArmBig, cdf::EncodingClass::BigIeee},
        {cdf::Encoding::DecStation, cdf::EncodingClass::LittleIeee},
        {cdf::Encoding::IbmPc, cdf::EncodingClass::LittleIeee},
        {cdf::Encoding::AlphaOsf1, cdf::EncodingClass::LittleIeee},
        {cdf::Encoding::AlphaVmsI, cdf::EncodingClass::LittleIeee},
        {cdf::Encoding::ArmLittle, cdf::EncodingClass::LittleIeee},
        {cdf::Encoding::Ia64VmsI, cdf::EncodingClass::LittleIeee},
        {cdf::Encoding::Vax, cdf::EncodingClass::LittleVaxD},
        {cdf::Encoding::AlphaVmsD, cdf::EncodingClass::LittleVaxD},
        {cdf::Encoding::Ia64VmsD, cdf::EncodingClass::LittleVaxD},
        {cdf::Encoding::AlphaVmsG, cdf::EncodingClass::LittleVaxG},
        {cdf::Encoding::Ia64VmsG, cdf::EncodingClass::LittleVaxG},
    };
    for (const Case& c : kCases) {
        EXPECT_EQ(det::encoding_class_of(c.enc), c.cls) << static_cast<int>(c.enc);
        EXPECT_NE(det::encoding_class_name(c.cls), "unknown");
    }
    // Host is a write-time selector, so it must resolve to whatever this machine is.
    EXPECT_EQ(det::encoding_class_of(cdf::Encoding::Host),
              cdf::kHostIsLittleEndian ? cdf::EncodingClass::LittleIeee
                                       : cdf::EncodingClass::BigIeee);
    EXPECT_TRUE(det::is_host_native(det::encoding_class_of(cdf::Encoding::Host)));
    // VAX floats are never IEEE, so they always need conversion regardless of byte order.
    EXPECT_FALSE(det::is_host_native(cdf::EncodingClass::LittleVaxD));
    EXPECT_FALSE(det::is_host_native(cdf::EncodingClass::LittleVaxG));
    // Exactly one of the two IEEE classes is free on any given host.
    EXPECT_NE(det::is_host_native(cdf::EncodingClass::BigIeee),
              det::is_host_native(cdf::EncodingClass::LittleIeee));
}

TEST(CdfTypes, EveryErrorCodeHasAMessage) {
    for (std::size_t i = 0; i < cdf::kErrorCodeCount; ++i) {
        const auto code = static_cast<cdf::ErrorCode>(i);
        EXPECT_NE(det::error_message(code), "unknown error") << "code " << i;
    }
    // The count is the ratchet: adding an enumerator without extending the message table, or
    // without bumping kErrorCodeCount, fails here.
    EXPECT_EQ(det::error_message(static_cast<cdf::ErrorCode>(cdf::kErrorCodeCount)),
              "unknown error");
}

// ---- bytes: decoding a real file ---------------------------------------------------------------

TEST(CdfBytes, DecodesTheOmniFileExactlyAsHandDecoded) {
    if (!have_omni()) { GTEST_SKIP() << "run: scripts/cdf-corpus.sh fetch --tier 1"; }
    const std::vector<std::byte> buf = read_all(omni_path().c_str());
    const det::Bytes b(buf.data(), buf.size());

    // Magic, then the CDR at a fixed offset.
    EXPECT_EQ(b.be_u32(0), cdf::kMagicV3);
    EXPECT_EQ(b.be_u32(4), cdf::kMagicUncompressed);
    EXPECT_EQ(b.be_i64(cdf::kCdrOffset), 312);
    EXPECT_EQ(b.be_i32(cdf::kCdrOffset + 8), static_cast<std::int32_t>(cdf::RecordType::Cdr));

    // The GDR sits immediately after the CDR — a self-consistency check the format guarantees.
    const std::uint64_t gdr = b.be_offset(cdf::kCdrOffset + 12);
    EXPECT_EQ(gdr, 320u);
    EXPECT_EQ(gdr, cdf::kCdrOffset + 312);

    EXPECT_EQ(b.be_i32(cdf::kCdrOffset + 28), static_cast<std::int32_t>(cdf::Encoding::Network));
    const std::uint32_t flags = b.be_u32(cdf::kCdrOffset + 32);
    EXPECT_TRUE((flags & cdf::kCdrFlagRowMajority) != 0);
    EXPECT_TRUE((flags & cdf::kCdrFlagSingleFile) != 0);
    EXPECT_FALSE((flags & cdf::kCdrFlagChecksum) != 0);

    EXPECT_EQ(b.be_i32(gdr + 8), static_cast<std::int32_t>(cdf::RecordType::Gdr));
    EXPECT_EQ(b.be_offset(gdr + 12), 0u);       // rVDRhead: this file has no rVariables
    EXPECT_EQ(b.be_i32(gdr + 44), 0);           // NrVars agrees
    EXPECT_EQ(b.be_i32(gdr + 60), 47);          // NzVars
    // No checksum, so eof is the whole file. With one it would be filesize - 16.
    EXPECT_EQ(static_cast<std::uint64_t>(b.be_i64(gdr + 36)), b.size());

    const std::uint64_t zvdr = b.be_offset(gdr + 20);
    EXPECT_EQ(zvdr, 21601u);
    EXPECT_EQ(b.be_i32(zvdr + 8),
              static_cast<std::int32_t>(cdf::RecordType::ZVariableDescriptor));
    EXPECT_EQ(b.name(zvdr + 84, 256), "Epoch");
    EXPECT_EQ(b.be_i32(zvdr + 20), static_cast<std::int32_t>(cdf::DataType::Epoch));
    EXPECT_EQ(b.be_i32(zvdr + 24), 44639);      // maxRec: one day at 1-minute cadence
    EXPECT_EQ(b.be_i32(zvdr + 64), 1);          // NumElems — at offset 64, not 56
}

TEST(CdfBytes, AccessorsAgreeOnKnownBytes) {
    const std::byte raw[] = {
        std::byte{0x01}, std::byte{0x02}, std::byte{0x03}, std::byte{0x04},
        std::byte{0x05}, std::byte{0x06}, std::byte{0x07}, std::byte{0x08},
        std::byte{0xFF}, std::byte{0xFF}, std::byte{0xFF}, std::byte{0xFE},
    };
    const det::Bytes b(raw, sizeof(raw));
    EXPECT_EQ(b.size(), sizeof(raw));
    EXPECT_EQ(b.data(), raw);
    EXPECT_EQ(b.u8(0), 0x01u);
    EXPECT_EQ(b.be_u16(0), 0x0102u);
    EXPECT_EQ(b.be_u32(0), 0x01020304u);
    EXPECT_EQ(b.be_u64(0), 0x0102030405060708ull);
    EXPECT_EQ(b.be_i32(0), 0x01020304);
    EXPECT_EQ(b.be_i64(0), 0x0102030405060708LL);
    EXPECT_EQ(b.be_i32(8), -2);  // sign extension, not a huge unsigned

    const det::Bytes sub = b.subspan(4, 4);
    EXPECT_EQ(sub.size(), 4u);
    EXPECT_EQ(sub.be_u32(0), 0x05060708u);

    // An empty view is legal and reads nothing.
    const det::Bytes empty;
    EXPECT_EQ(empty.size(), 0u);
    EXPECT_EQ(empty.data(), nullptr);
    EXPECT_TRUE(empty.contains(0, 0));
    EXPECT_FALSE(empty.contains(0, 1));
}

TEST(CdfBytes, NamesAreStrippedOfPadding) {
    std::vector<std::byte> raw(300, std::byte{0});
    const char* name = "Epoch";
    std::memcpy(raw.data(), name, 5);
    const det::Bytes b(raw.data(), raw.size());
    EXPECT_EQ(b.name(0, 256), "Epoch");
    // A name that fills the field has no terminator; the width is the limit.
    std::vector<std::byte> full(8, std::byte{'A'});
    const det::Bytes fb(full.data(), full.size());
    EXPECT_EQ(fb.name(0, 8), "AAAAAAAA");
    EXPECT_EQ(fb.name(0, 8).size(), 8u);
}

TEST(CdfBytes, EveryGuardActuallyBites) {
    const std::byte raw[16] = {};
    const det::Bytes b(raw, sizeof(raw));

    // Reading past the end.
    EXPECT_THROW((void)b.be_u32(14), cdf::CdfError);
    EXPECT_THROW((void)b.be_u64(9), cdf::CdfError);
    EXPECT_THROW((void)b.u8(16), cdf::CdfError);
    EXPECT_THROW((void)b.subspan(8, 9), cdf::CdfError);
    EXPECT_THROW((void)b.name(12, 8), cdf::CdfError);

    // An offset near 2^64 must not wrap into range. This is the check that a naive
    // `offset + count <= size` would fail.
    EXPECT_FALSE(b.contains(0xFFFFFFFFFFFFFFF0ull, 32));
    EXPECT_THROW((void)b.subspan(0xFFFFFFFFFFFFFFF0ull, 32), cdf::CdfError);

    // Exactly-at-the-edge reads are legal; one past is not.
    EXPECT_NO_THROW((void)b.be_u64(8));
    EXPECT_THROW((void)b.be_u64(9), cdf::CdfError);
    EXPECT_TRUE(b.contains(16, 0));
    EXPECT_FALSE(b.contains(17, 0));
}

TEST(CdfBytes, OffsetsAreValidatedNotTrusted) {
    std::vector<std::byte> raw(64, std::byte{0});
    const det::Bytes b(raw.data(), raw.size());

    EXPECT_EQ(b.be_offset(0), 0u);  // zero means "no such record" and is legal

    // Negative offset.
    std::fill(raw.begin(), raw.begin() + 8, std::byte{0xFF});
    EXPECT_THROW((void)b.be_offset(0), cdf::CdfError);

    // Offset past the end of the file.
    std::fill(raw.begin(), raw.begin() + 8, std::byte{0});
    raw[7] = std::byte{0xFF};  // 255, well past a 64-byte view
    EXPECT_THROW((void)b.be_offset(0), cdf::CdfError);

    // The largest in-range offset is the size itself (a zero-length record at EOF).
    std::fill(raw.begin(), raw.begin() + 8, std::byte{0});
    raw[7] = std::byte{64};
    EXPECT_EQ(b.be_offset(0), 64u);
}

TEST(CdfBytes, ErrorsCarryCodeAndOffset) {
    const std::byte raw[4] = {};
    const det::Bytes b(raw, sizeof(raw));
    try {
        (void)b.be_u64(2);
        FAIL() << "expected a CdfError";
    } catch (const cdf::CdfError& e) {
        EXPECT_EQ(e.code(), cdf::ErrorCode::TruncatedFile);
        EXPECT_EQ(e.offset(), 2u);
        // The message names the offset — "bad record" alone is not actionable on an 8 MB file.
        EXPECT_NE(std::string(e.what()).find("byte 2"), std::string::npos) << e.what();
    }
}

// ---- mapping ---------------------------------------------------------------------------------

TEST(CdfMapping, MapsARealFile) {
    if (!have_omni()) { GTEST_SKIP() << "run: scripts/cdf-corpus.sh fetch --tier 1"; }
    det::FileMapping m(omni_path());
    EXPECT_TRUE(m.is_open());
    EXPECT_EQ(m.size(), 8774208u);
    EXPECT_EQ(m.bytes().be_u32(0), cdf::kMagicV3);
    EXPECT_EQ(m.bytes().name(21601 + 84, 256), "Epoch");
}

TEST(CdfMapping, MoveTransfersOwnershipExactlyOnce) {
    if (!have_omni()) { GTEST_SKIP() << "run: scripts/cdf-corpus.sh fetch --tier 1"; }
    det::FileMapping a(omni_path());
    const std::uint64_t size = a.size();

    det::FileMapping b(std::move(a));
    EXPECT_TRUE(b.is_open());
    EXPECT_EQ(b.size(), size);
    // Reading the moved-from object IS the assertion: a move must leave it owning nothing, or
    // the destructor unmaps twice.
    // NOLINTNEXTLINE(bugprone-use-after-move,clang-analyzer-cplusplus.Move)
    EXPECT_FALSE(a.is_open());
    // NOLINTNEXTLINE(bugprone-use-after-move,clang-analyzer-cplusplus.Move)
    EXPECT_EQ(a.size(), 0u);

    det::FileMapping c;
    EXPECT_FALSE(c.is_open());
    c = std::move(b);
    EXPECT_TRUE(c.is_open());
    EXPECT_EQ(c.bytes().be_u32(0), cdf::kMagicV3);

    // Self-move must not unmap what it is about to keep. Valgrind and ASan are the real judges.
    det::FileMapping* alias = &c;
    c = std::move(*alias);
    EXPECT_TRUE(c.is_open());  // NOLINT(clang-analyzer-cplusplus.Move) — self-move
    EXPECT_EQ(c.size(), size);  // NOLINT(clang-analyzer-cplusplus.Move)
}

TEST(CdfMapping, RealFilesystemFailures) {
    try {
        det::FileMapping m("/nonexistent/definitely/not/here.cdf");
        FAIL() << "expected a CdfError";
    } catch (const cdf::CdfError& e) {
        EXPECT_EQ(e.code(), cdf::ErrorCode::CannotOpen);
    }

    // A file too small to hold even the two magic numbers, through the REAL syscalls.
    const std::string tiny = std::string(std::tmpnam(nullptr)) + ".cdf";
    { std::ofstream out(tiny, std::ios::binary); out << "abc"; }
    try {
        det::FileMapping m(tiny);
        FAIL() << "expected a CdfError";
    } catch (const cdf::CdfError& e) {
        EXPECT_EQ(e.code(), cdf::ErrorCode::EmptyFile);
    }
    // Best-effort teardown; a stale temp file is not an error a test can act on.
    (void)std::remove(tiny.c_str());  // NOLINT(cert-err33-c): see above

    // A directory opens and stats to a plausible size, then refuses to map — the one way to
    // reach the real mmap-failure path without a fault-injecting policy.
    try {
        det::FileMapping m("/tmp");
        FAIL() << "expected a CdfError";
    } catch (const cdf::CdfError& e) {
        EXPECT_EQ(e.code(), cdf::ErrorCode::CannotMap);
    }
}

// The syscall failure paths. No CDF, however malformed, can make mmap fail — so without this
// seam these branches would be permanently uncovered and, worse, never exercised at all.
TEST(CdfMapping, SyscallFailuresAreHandledAndLeakNothing) {
    using Failing = det::BasicFileMapping<FailingSysOps>;
    try {
        Failing m("anything");
        FAIL() << "expected a CdfError";
    } catch (const cdf::CdfError& e) {
        EXPECT_EQ(e.code(), cdf::ErrorCode::CannotOpen);
    }

    using Empty = det::BasicFileMapping<EmptySysOps>;
    try {
        Empty m("anything");
        FAIL() << "expected a CdfError";
    } catch (const cdf::CdfError& e) {
        EXPECT_EQ(e.code(), cdf::ErrorCode::EmptyFile);
    }

    using MapFails = det::BasicFileMapping<MapFailsSysOps>;
    try {
        MapFails m("anything");
        FAIL() << "expected a CdfError";
    } catch (const cdf::CdfError& e) {
        EXPECT_EQ(e.code(), cdf::ErrorCode::CannotMap);
    }

    // A default-constructed mapping owns nothing and destructs cleanly.
    const Failing none;
    EXPECT_FALSE(none.is_open());
    EXPECT_EQ(none.size(), 0u);
    EXPECT_EQ(none.bytes().size(), 0u);
}
