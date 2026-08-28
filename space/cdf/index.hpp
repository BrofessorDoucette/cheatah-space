#pragma once

/**
 * @file index.hpp
 * @brief space.cdf — flattening a variable's VXR tree into a sorted list of record extents.
 *
 * A variable's records are found through its VXRs, and the VXRs form a tree: an entry may point
 * at a VVR, a CVVR, or another VXR, and each VXR may also have a `next` sibling. NASA's library
 * walks that structure from the root on every access. This walks it once per variable into a
 * flat vector sorted by first record, so every later lookup is a binary search.
 *
 * The walk is **iterative, with an explicit stack, and capped** on both depth and total node
 * count. That is not defensive decoration: a hostile file makes `VXRnext` point at itself in
 * one byte, and a recursive walker then either overflows the C stack or spins forever. Both are
 * QA-gate failures (ASan, or the test timeout) before they are anything else.
 *
 * The tree is real, not theoretical. The OMNI `Epoch` variable in the corpus has a two-level
 * index — a root whose three entries each point at another VXR — with a `next` chain on top,
 * for 24 leaf records in all. It was the first real file this reader opened.
 */

#include "cheatah.hpp"

#include <algorithm>
#include <cstdint>
#include <vector>

#include "bytes.hpp"
#include "records.hpp"
#include "types.hpp"

namespace cheatah::space::cdf::detail {

/// One contiguous run of records and where their bytes live.
///
/// An extent may run PAST the variable's last written record: VVRs are allocated at
/// blocking-factor granularity, so the final leaf of the OMNI `Epoch` index covers records
/// 44032..45055 while `maxRec` is 44639. The index describes allocation; `Vdr::record_count()`
/// describes truth, and every reader must clamp to it.
struct RecordExtent {
    std::int64_t first{};      ///< First record number covered.
    std::int64_t last{};       ///< Last record number covered, inclusive.
    std::uint64_t offset{};    ///< The VVR or CVVR holding them.
    bool compressed{};         ///< CVVR (true) or VVR (false).

    /// @return how many records this extent holds.
    [[nodiscard]] std::int64_t count() const noexcept { return last - first + 1; }
};

/// How deep an index tree may nest before it is judged a cycle.
inline constexpr std::int32_t kMaxVxrDepth = 64;
/// How many VXR nodes a single variable may have before it is judged a cycle.
inline constexpr std::uint64_t kMaxVxrNodes = 1'000'000;

/**
 * Flatten a variable's index into extents sorted by first record.
 *
 * @param b the file.
 * @param vxr_head Vdr::vxr_head; 0 yields an empty index.
 * @return every leaf extent, sorted ascending by `first`.
 */
inline std::vector<RecordExtent> build_index(const Bytes& b, std::uint64_t vxr_head) {
    std::vector<RecordExtent> extents;
    if (vxr_head == 0) { return extents; }

    struct Pending { std::uint64_t offset; std::int32_t depth; };
    std::vector<Pending> stack;
    stack.push_back({vxr_head, 0});
    std::uint64_t nodes = 0;

    while (!stack.empty()) {
        const Pending cur = stack.back();
        stack.pop_back();
        require(cur.depth < kMaxVxrDepth, ErrorCode::VxrTreeTooDeep, cur.offset);
        require(++nodes <= kMaxVxrNodes, ErrorCode::VxrTreeTooLarge, cur.offset);

        const Vxr x = parse_vxr(b, cur.offset);
        for (std::int32_t i = 0; i < x.nused; ++i) {
            const std::int32_t first = x.first(b, i);
            const std::int32_t last = x.last(b, i);
            require(first >= 0 && last >= first, ErrorCode::RecordOutOfRange, x.offset);
            const std::uint64_t child = x.child_offset(b, i);
            const RecordHeader h = read_record_header(b, child);
            if (h.type == RecordType::Vxr) {
                stack.push_back({child, cur.depth + 1});
            } else {
                require(h.type == RecordType::Vvr || h.type == RecordType::Cvvr,
                        ErrorCode::UnexpectedRecordType, child + 8);
                extents.push_back({first, last, child, h.type == RecordType::Cvvr});
            }
        }
        // A sibling is at the same depth as the node it follows; a cycle through `next` is
        // caught by the node cap rather than the depth cap.
        if (x.next != 0) { stack.push_back({x.next, cur.depth}); }
    }

    std::sort(extents.begin(), extents.end(),
              [](const RecordExtent& lhs, const RecordExtent& rhs) { return lhs.first < rhs.first; });
    return extents;
}

/**
 * The extent containing @p record, if any.
 *
 * @param extents a sorted index from build_index().
 * @param record the record number wanted.
 * @return a pointer into @p extents, or nullptr when no extent covers @p record (a sparse gap).
 */
inline const RecordExtent* find_extent(const std::vector<RecordExtent>& extents, std::int64_t record) noexcept {
    auto it = std::upper_bound(extents.begin(), extents.end(), record,
                               [](std::int64_t r, const RecordExtent& e) { return r < e.first; });
    if (it == extents.begin()) { return nullptr; }
    --it;
    return (record <= it->last) ? &*it : nullptr;
}

}  // namespace cheatah::space::cdf::detail
