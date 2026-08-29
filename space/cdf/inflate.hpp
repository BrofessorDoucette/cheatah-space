#pragma once

/**
 * @file inflate.hpp
 * @brief space.cdf — DEFLATE decompression (RFC 1951), from scratch, with no zlib.
 *
 * CDF's GZIP mode is the only compression that occurs in NASA's public archive — every
 * compressed variable in every corpus file uses it — so this is what turns "opens OMNI" into
 * "opens RBSP, MMS and THEMIS". The payload inside a CVVR is a full gzip container (1f 8b 08 00),
 * verified on both `test_alltypes.cdf` and an RBSP HOPE file; zlib and raw-deflate streams are
 * accepted as well, because the format does not promise which a writer produces.
 *
 * The decoder is a plain canonical-Huffman implementation: a 64-bit bit reservoir refilled a
 * byte at a time, and per-symbol decode walking code lengths. It is not the fastest possible
 * shape — a multi-level lookup table would beat it — but it is small enough to read, and the
 * benchmark work that would justify the table is deliberately deferred. What matters here is
 * that it is correct and that it cannot be talked into running away: every output write is
 * bounds-checked against the size the caller already knows from the record index, and any
 * disagreement is an error rather than a resize.
 *
 * Everything is bounds-checked because this parses untrusted files. There is no allocation
 * beyond the caller's output buffer and two small fixed tables.
 */

#include "cheatah.hpp"

#include <array>
#include <cstddef>
#include <cstdint>
#include <cstring>

#include "bytes.hpp"
#include "types.hpp"

namespace cheatah::space::cdf::detail {

/// Largest DEFLATE code length, and so the number of length buckets a canonical table needs.
inline constexpr int kMaxCodeBits = 15;
/// Literal/length alphabet size.
inline constexpr int kLitLenSymbols = 288;
/// Distance alphabet size.
inline constexpr int kDistSymbols = 32;
/// Code-length alphabet size.
inline constexpr int kCodeLenSymbols = 19;

/// A canonical Huffman decoding table: counts per bit length, plus symbols in canonical order.
struct HuffTable {
    std::array<std::uint16_t, kMaxCodeBits + 1> count{};  ///< How many codes have each length.
    std::array<std::uint16_t, kLitLenSymbols> symbol{};   ///< Symbols, shortest code first.
};

/// Build a canonical table from code lengths.
/// @param lengths one length per symbol; 0 means the symbol is unused.
/// @param n how many symbols. @param out the table to fill. @param offset for errors.
inline void build_huffman(const std::uint8_t* lengths, int n, HuffTable& out, std::uint64_t offset) {
    out.count.fill(0);
    for (int i = 0; i < n; ++i) {
        require(lengths[i] <= kMaxCodeBits, ErrorCode::DecompressionFailed, offset);
        ++out.count[lengths[i]];
    }
    out.count[0] = 0;
    // An over-subscribed set would make two symbols share a code; an incomplete one leaves a
    // code that decodes to nothing. Both mean the stream is malformed, not merely unusual.
    int left = 1;
    for (int len = 1; len <= kMaxCodeBits; ++len) {
        left <<= 1;
        left -= out.count[len];
        require(left >= 0, ErrorCode::DecompressionFailed, offset);
    }
    std::array<std::uint16_t, kMaxCodeBits + 1> offs{};
    for (int len = 1; len < kMaxCodeBits; ++len) {
        offs[len + 1] = static_cast<std::uint16_t>(offs[len] + out.count[len]);
    }
    for (int i = 0; i < n; ++i) {
        if (lengths[i] != 0) { out.symbol[offs[lengths[i]]++] = static_cast<std::uint16_t>(i); }
    }
}

/// A DEFLATE stream being decoded: the input, a bit reservoir, and the output window.
class Inflater {
  public:
    /// @param in the compressed bytes. @param out where to write. @param out_size how many bytes
    /// the caller already knows to expect — the record index says so, and a stream that disagrees
    /// is an error rather than a reason to grow a buffer.
    Inflater(Bytes in, std::byte* out, std::uint64_t out_size) noexcept
        : in_(in), out_(out), out_size_(out_size) {}

    /// Decode until the final block. @return how many bytes were written.
    std::uint64_t run() {
        bool final_block = false;
        while (!final_block) {
            final_block = (bits(1) != 0);
            const std::uint32_t type = bits(2);
            if (type == 0) { stored(); }
            else if (type == 1) { fixed(); }
            else if (type == 2) { dynamic(); }
            else { throw_cdf(ErrorCode::DecompressionFailed, pos_); }
        }
        return written_;
    }

  private:
    /// Read @p n bits, refilling the reservoir a byte at a time.
    /// @param n 0..24. @return the bits, LSB-first as DEFLATE stores them.
    std::uint32_t bits(int n) {
        while (nbits_ < n) {
            require(pos_ < in_.size(), ErrorCode::DecompressionFailed, pos_);
            reservoir_ |= static_cast<std::uint64_t>(in_.u8(pos_)) << nbits_;
            ++pos_;
            nbits_ += 8;
        }
        const auto v = static_cast<std::uint32_t>(reservoir_ & ((1ULL << n) - 1));
        reservoir_ >>= n;
        nbits_ -= n;
        return v;
    }

    /// Decode one symbol with @p t.
    /// @param t a built table. @return the symbol.
    int decode(const HuffTable& t) {
        int code = 0;
        int first = 0;
        int index = 0;
        for (int len = 1; len <= kMaxCodeBits; ++len) {
            code |= static_cast<int>(bits(1));
            const int count = t.count[len];
            if (code - first < count) { return t.symbol[static_cast<std::size_t>(index + (code - first))]; }
            index += count;
            first = (first + count) << 1;
            code <<= 1;
        }
        throw_cdf(ErrorCode::DecompressionFailed, pos_);
    }

    /// Copy an uncompressed block.
    void stored() {
        reservoir_ = 0;
        nbits_ = 0;   // stored blocks resume on a byte boundary
        require(pos_ + 4 <= in_.size(), ErrorCode::DecompressionFailed, pos_);
        const std::uint32_t len = static_cast<std::uint32_t>(in_.u8(pos_))
                                  | (static_cast<std::uint32_t>(in_.u8(pos_ + 1)) << 8);
        const std::uint32_t nlen = static_cast<std::uint32_t>(in_.u8(pos_ + 2))
                                   | (static_cast<std::uint32_t>(in_.u8(pos_ + 3)) << 8);
        require((len ^ 0xFFFFU) == nlen, ErrorCode::DecompressionFailed, pos_);
        pos_ += 4;
        require(pos_ + len <= in_.size(), ErrorCode::DecompressionFailed, pos_);
        require(written_ + len <= out_size_, ErrorCode::DecompressedSizeMismatch, pos_);
        std::memcpy(out_ + written_, in_.data() + pos_, len);
        pos_ += len;
        written_ += len;
    }

    /// Decode a block with the fixed tables RFC 1951 section 3.2.6 defines.
    void fixed() {
        std::array<std::uint8_t, kLitLenSymbols> ll{};
        for (int i = 0; i < 144; ++i) { ll[i] = 8; }
        for (int i = 144; i < 256; ++i) { ll[i] = 9; }
        for (int i = 256; i < 280; ++i) { ll[i] = 7; }
        for (int i = 280; i < kLitLenSymbols; ++i) { ll[i] = 8; }
        std::array<std::uint8_t, kDistSymbols> d{};
        d.fill(5);
        HuffTable lt;
        HuffTable dt;
        build_huffman(ll.data(), kLitLenSymbols, lt, pos_);
        build_huffman(d.data(), kDistSymbols, dt, pos_);
        block(lt, dt);
    }

    /// Read the per-block tables, then decode with them.
    void dynamic() {
        const int hlit = static_cast<int>(bits(5)) + 257;
        const int hdist = static_cast<int>(bits(5)) + 1;
        const int hclen = static_cast<int>(bits(4)) + 4;
        require(hlit <= kLitLenSymbols && hdist <= kDistSymbols, ErrorCode::DecompressionFailed, pos_);
        // The code-length alphabet is itself Huffman-coded, in this fixed permutation.
        static constexpr std::array<std::uint8_t, kCodeLenSymbols> kOrder = {
            16, 17, 18, 0, 8, 7, 9, 6, 10, 5, 11, 4, 12, 3, 13, 2, 14, 1, 15};
        std::array<std::uint8_t, kCodeLenSymbols> cl{};
        for (int i = 0; i < hclen; ++i) { cl[kOrder[static_cast<std::size_t>(i)]] = static_cast<std::uint8_t>(bits(3)); }
        HuffTable ct;
        build_huffman(cl.data(), kCodeLenSymbols, ct, pos_);

        std::array<std::uint8_t, kLitLenSymbols + kDistSymbols> lengths{};
        int i = 0;
        while (i < hlit + hdist) {
            const int sym = decode(ct);
            if (sym < 16) {
                lengths[static_cast<std::size_t>(i++)] = static_cast<std::uint8_t>(sym);
            } else {
                std::uint8_t repeat_value = 0;
                int repeat = 0;
                if (sym == 16) {
                    require(i > 0, ErrorCode::DecompressionFailed, pos_);
                    repeat_value = lengths[static_cast<std::size_t>(i - 1)];
                    repeat = 3 + static_cast<int>(bits(2));
                } else if (sym == 17) {
                    repeat = 3 + static_cast<int>(bits(3));
                } else {
                    repeat = 11 + static_cast<int>(bits(7));
                }
                require(i + repeat <= hlit + hdist, ErrorCode::DecompressionFailed, pos_);
                for (int k = 0; k < repeat; ++k) { lengths[static_cast<std::size_t>(i++)] = repeat_value; }
            }
        }
        HuffTable lt;
        HuffTable dt;
        build_huffman(lengths.data(), hlit, lt, pos_);
        build_huffman(lengths.data() + hlit, hdist, dt, pos_);
        block(lt, dt);
    }

    /// Decode literals and back-references until the end-of-block symbol.
    /// @param lt the literal/length table. @param dt the distance table.
    void block(const HuffTable& lt, const HuffTable& dt) {
        static constexpr std::array<std::uint16_t, 29> kLenBase = {
            3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 15, 17, 19, 23, 27, 31, 35, 43, 51, 59,
            67, 83, 99, 115, 131, 163, 195, 227, 258};
        static constexpr std::array<std::uint8_t, 29> kLenExtra = {
            0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 0};
        static constexpr std::array<std::uint16_t, 30> kDistBase = {
            1, 2, 3, 4, 5, 7, 9, 13, 17, 25, 33, 49, 65, 97, 129, 193, 257, 385, 513, 769,
            1025, 1537, 2049, 3073, 4097, 6145, 8193, 12289, 16385, 24577};
        static constexpr std::array<std::uint8_t, 30> kDistExtra = {
            0, 0, 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, 10, 11, 11, 12, 12, 13, 13};

        for (;;) {
            const int sym = decode(lt);
            if (sym < 256) {
                require(written_ < out_size_, ErrorCode::DecompressedSizeMismatch, pos_);
                out_[written_++] = static_cast<std::byte>(sym);
            } else if (sym == 256) {
                return;
            } else {
                const auto li = static_cast<std::size_t>(sym - 257);
                require(li < kLenBase.size(), ErrorCode::DecompressionFailed, pos_);
                const std::uint64_t length = kLenBase[li] + bits(kLenExtra[li]);
                const int dsym = decode(dt);
                require(static_cast<std::size_t>(dsym) < kDistBase.size(), ErrorCode::DecompressionFailed, pos_);
                const std::uint64_t dist = kDistBase[static_cast<std::size_t>(dsym)]
                                           + bits(kDistExtra[static_cast<std::size_t>(dsym)]);
                // A distance reaching before the start of the output is the classic malformed
                // stream; it would otherwise read whatever memory precedes the buffer.
                require(dist <= written_, ErrorCode::DecompressionFailed, pos_);
                require(written_ + length <= out_size_, ErrorCode::DecompressedSizeMismatch, pos_);
                // Byte at a time on purpose: overlapping copies are legal and load-bearing in
                // DEFLATE (a distance of 1 is how a run of one byte is encoded), so memcpy or
                // memmove would both be wrong.
                for (std::uint64_t k = 0; k < length; ++k) {
                    out_[written_] = out_[written_ - dist];
                    ++written_;
                }
            }
        }
    }

    Bytes in_;
    std::byte* out_;
    std::uint64_t out_size_;
    std::uint64_t pos_ = 0;
    std::uint64_t written_ = 0;
    std::uint64_t reservoir_ = 0;
    int nbits_ = 0;
};

/**
 * Decompress a DEFLATE stream, in whichever container it arrives in.
 *
 * gzip (1f 8b) and zlib (0x78 and friends) headers are stripped; anything else is treated as a
 * raw stream. CDF's own writer emits gzip, but the format does not promise it, so all three are
 * accepted.
 *
 * @param in the compressed bytes.
 * @param out where to write.
 * @param out_size exactly how many bytes are expected — from the record index.
 * @param offset the file offset, for errors.
 */
inline void inflate(Bytes in, std::byte* out, std::uint64_t out_size, std::uint64_t offset) {
    require(in.size() >= 2, ErrorCode::DecompressionFailed, offset);
    std::uint64_t start = 0;
    if (in.u8(0) == 0x1F && in.u8(1) == 0x8B) {
        // gzip: 10-byte header, then the optional extra/name/comment/CRC fields.
        require(in.size() >= 10, ErrorCode::DecompressionFailed, offset);
        require(in.u8(2) == 8, ErrorCode::DecompressionFailed, offset);   // DEFLATE, the only method
        const std::uint8_t flags = in.u8(3);
        start = 10;
        if ((flags & 0x04) != 0) {                                        // FEXTRA
            require(start + 2 <= in.size(), ErrorCode::DecompressionFailed, offset);
            const std::uint64_t xlen = in.u8(start) | (static_cast<std::uint64_t>(in.u8(start + 1)) << 8);
            start += 2 + xlen;
        }
        if ((flags & 0x08) != 0) {                                        // FNAME
            while (start < in.size() && in.u8(start) != 0) { ++start; }
            ++start;
        }
        if ((flags & 0x10) != 0) {                                        // FCOMMENT
            while (start < in.size() && in.u8(start) != 0) { ++start; }
            ++start;
        }
        if ((flags & 0x02) != 0) { start += 2; }                          // FHCRC
        require(start <= in.size(), ErrorCode::DecompressionFailed, offset);
    } else if ((in.u8(0) & 0x0F) == 8 && ((in.u8(0) << 8) | in.u8(1)) % 31 == 0) {
        start = 2;                                                        // zlib
    }

    Inflater z(in.subspan(start, in.size() - start), out, out_size);
    const std::uint64_t got = z.run();
    require(got == out_size, ErrorCode::DecompressedSizeMismatch, offset);
}

}  // namespace cheatah::space::cdf::detail
