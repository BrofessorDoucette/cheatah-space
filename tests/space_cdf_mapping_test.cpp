// The real-syscall lane of space.cdf's file mapping, plus the small enum tails coverage showed
// were never driven.
//
// mapping.hpp puts the mmap syscalls behind a SysOps seam so the mapping logic can be tested with
// a fake — and then only the fake got tested: PosixSysOps and the BasicFileMapping<PosixSysOps>
// instantiation (the one every real CDF read uses) sat at 0% while the suite read green. A seam
// whose REAL side is never exercised has inverted its purpose: the test double is covered and the
// production code is not. These tests map an actual temporary file through actual mmap.
#include <gtest/gtest.h>

#include <cstdint>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <string>
#include <utility>

#include "space/cdf/mapping.hpp"
#include "space/cdf/types.hpp"

namespace cdf = cheatah::space::cdf;

namespace {
/// A real file on disk with known bytes, removed on destruction.
struct TempFile {
    std::string path;
    explicit TempFile(const std::string& contents)
        : path(std::string(::testing::TempDir()) + "cdf_mapping_probe.bin") {
        std::ofstream out(path, std::ios::binary | std::ios::trunc);
        out.write(contents.data(), static_cast<std::streamsize>(contents.size()));
    }
    ~TempFile() {
        // Best-effort teardown: a temp file that outlives a failed remove is a stale byte-blob in
        // the test tempdir, not an error a destructor could act on — and throwing from here would
        // turn one test's cleanup into the whole run's abort.
        (void)std::remove(path.c_str());  // NOLINT(cert-err33-c): see above
    }
};
}  // namespace

TEST(CdfMappingPosix, MapsARealFileThroughRealSyscalls) {
    const std::string bytes = "CDF mapping probe: 32 real bytes";
    const TempFile f(bytes);

    cdf::detail::FileMapping map(f.path);
    ASSERT_TRUE(map.is_open()) << "mapping a readable file must succeed";
    ASSERT_EQ(bytes.size(), map.size());
    // The mapped memory IS the file: byte-for-byte, through a real mmap, not a fake.
    EXPECT_EQ(0, std::memcmp(map.bytes().data(), bytes.data(), bytes.size()));
}

TEST(CdfMappingPosix, MoveTransfersOwnershipExactlyOnce) {
    const TempFile f("move semantics probe");
    cdf::detail::FileMapping a(f.path);
    ASSERT_TRUE(a.is_open());
    const std::uint64_t size = a.size();

    // Move construction: the source must be left closed, or the destructor unmaps twice.
    cdf::detail::FileMapping b(std::move(a));
    EXPECT_TRUE(b.is_open());
    EXPECT_EQ(size, b.size());
    EXPECT_FALSE(a.is_open());  // NOLINT(bugprone-use-after-move): the post-move state IS the test

    // Move assignment over an open mapping: the target's old mapping must be released first.
    cdf::detail::FileMapping c(f.path);
    ASSERT_TRUE(c.is_open());
    c = std::move(b);
    EXPECT_TRUE(c.is_open());
    EXPECT_FALSE(b.is_open());  // NOLINT(bugprone-use-after-move)
}

TEST(CdfMappingPosix, AMissingFileThrowsAndATinyFileIsRefused) {
    // The constructor's contract is to THROW, not to hand back a closed object a caller must
    // remember to test: a mapping either exists or the reason it does not is in the exception.
    EXPECT_THROW(cdf::detail::FileMapping("/nonexistent/definitely/not/here.cdf"),
                 std::exception);

    // A file below the two magic numbers cannot be a CDF, and a zero-length mmap is UB on POSIX —
    // so the too-small path must close the descriptor it already opened and then throw. The close
    // is the part a leak would hide in, which is why this real-file case exists rather than
    // trusting the fake-SysOps test alone.
    const TempFile tiny("x");
    EXPECT_THROW(cdf::detail::FileMapping{tiny.path}, std::exception);
}

TEST(CdfTypes, EveryEncodingClassAndCompressionNameIsReachable) {
    using cdf::EncodingClass;
    using cdf::CompressionType;
    // The enum tails the sweep never hit: an out-of-enum encoding falls back to BigIeee (the
    // conservative read for CDF, whose default heritage encoding is big-endian), and every
    // compression enumerator plus the unknown fall-through names itself.
    EXPECT_EQ(EncodingClass::BigIeee,
              cdf::detail::encoding_class_of(static_cast<cdf::Encoding>(255)));
    EXPECT_EQ("none", cdf::detail::compression_name(CompressionType::None));
    EXPECT_EQ("rle", cdf::detail::compression_name(CompressionType::Rle));
    EXPECT_EQ("huff", cdf::detail::compression_name(CompressionType::Huff));
    EXPECT_EQ("ahuff", cdf::detail::compression_name(CompressionType::Ahuff));
    EXPECT_EQ("gzip", cdf::detail::compression_name(CompressionType::Gzip));
    EXPECT_EQ("unknown", cdf::detail::compression_name(static_cast<CompressionType>(200)));
}
