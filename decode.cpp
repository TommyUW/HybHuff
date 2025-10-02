// ============================================================================
// File: decode_blocked.cpp  (self-contained helpers + paged decoder)
// Functions:
//   • extract_block_slices: copy the minimal per-block subarrays and bit windows
//   • decode_block_partial: decode entries [i0, i1) using those slices
//   • decodeSide (paged)  : iterate blocks → extract → decode → (optional) verify
// Key idea:
//   Work in block-aligned windows, then locally advance bit pointers within the
//   copied slices to decode each entry without touching unrelated data.
// ============================================================================

#include <vector>
#include <cstdint>
#include <cstring>
#include <algorithm>
#include <iostream>

#include "decode.h"       // readBits(...), readBits64(...), bytes_for_bits(...)
#include "huffman_tree.h" // struct HuffmanNode { HuffmanNode* l,*r; int v; }
#include "misc.h"

using std::vector;
using std::uint8_t;
using std::uint16_t;
using std::uint64_t;

/**
 * extract_block_slices
 * --------------------
 * Extract per-entry metadata and byte-aligned bit-slices covering entries [i0, i1).
 *
 * Inputs:
 *   i0, i1           : half-open entry range [i0, i1)
 *   BLOCK            : block size in entries
 *   deg, bitCount,
 *   huffCount        : per-entry metadata arrays (global)
 *   neighHi, neighLo : global bitstreams (hi = Huffman section, lo = fallback section)
 *   blockHi, blockLo : bit offsets at block starts; block k begins at entry k*BLOCK
 *   fallbackBitsOpp  : fixed width (bits) for fallback-coded values
 *
 * Outputs (sized/filled by the function):
 *   deg_blk, bit_blk, huff_blk : slices of per-entry metadata for [i0, i1)
 *   hiSlice, loSlice           : byte windows copied from neighHi/neighLo
 *   hiBitOffset, loBitOffset   : starting bit offsets *inside the first byte*
 *                                of hiSlice/loSlice (for the entry i0)
 *
 * Notes:
 *   • The pair (hiSlice, hiBitOffset) and (loSlice, loBitOffset) together define
 *     bit pointers where local decoding should start for entry i0.
 *   • The caller advances local bit pointers per entry using bit_blk/huff_blk/deg_blk.
 */
void extract_block_slices(
    int i0, int i1,
    int BLOCK,
    const int*       deg,
    const uint16_t*  bitCount,
    const uint16_t*  huffCount,
    const uint8_t*   neighHi,
    const uint8_t*   neighLo,
    const uint64_t*  blockHi,       // bit offsets at block starts
    const uint64_t*  blockLo,       // bit offsets at block starts
    int              fallbackBitsOpp,
    // outputs
    std::vector<int>&       deg_blk,
    std::vector<uint16_t>&  bit_blk,
    std::vector<uint16_t>&  huff_blk,
    std::vector<uint8_t>&   hiSlice,
    std::vector<uint8_t>&   loSlice,
    uint32_t&               hiBitOffset,
    uint32_t&               loBitOffset
) {
    const int len = std::max(0, i1 - i0);
    deg_blk.resize(len);
    bit_blk.resize(len);
    huff_blk.resize(len);
    for (int k = 0; k < len; ++k) {
        deg_blk[k]  = deg[i0 + k];
        bit_blk[k]  = bitCount[i0 + k];
        huff_blk[k] = huffCount[i0 + k];
    }

    // Block containing entry i0
    const int blk         = (BLOCK > 0) ? (i0 / BLOCK) : 0;
    const int blkStartIdx = blk * BLOCK;

    // Bits within the block prior to i0 (to position the local window)
    uint64_t hiBitsBefore = 0, loBitsBefore = 0;
    for (int t = blkStartIdx; t < i0; ++t) {
        hiBitsBefore += bitCount[t];
        loBitsBefore += (uint64_t)(deg[t] - huffCount[t]) * (uint64_t)fallbackBitsOpp;
    }

    // Absolute bit starts for entry i0, then split into byte-aligned copy + intra-byte offset
    const uint64_t hiStartBit = blockHi[blk] + hiBitsBefore;
    const uint64_t loStartBit = blockLo[blk] + loBitsBefore;

    // Bits needed to cover [i0, i1)
    uint64_t hiBitsThis = 0, loBitsThis = 0;
    for (int i = i0; i < i1; ++i) {
        hiBitsThis += bitCount[i];
        loBitsThis += (uint64_t)(deg[i] - huffCount[i]) * (uint64_t)fallbackBitsOpp;
    }

    // Copy byte windows; return the intra-byte starting offsets for i0
    const size_t hiSrcByte0 = (size_t)(hiStartBit >> 3);
    const size_t loSrcByte0 = (size_t)(loStartBit >> 3);
    hiBitOffset = (uint32_t)(hiStartBit & 7ULL);
    loBitOffset = (uint32_t)(loStartBit & 7ULL);

    const size_t hiBytes = bytes_for_bits(hiBitOffset + hiBitsThis);
    const size_t loBytes = bytes_for_bits(loBitOffset + loBitsThis);

    hiSlice.resize(hiBytes);
    loSlice.resize(loBytes);

    if (hiBytes) std::memcpy(hiSlice.data(), neighHi + hiSrcByte0, hiBytes);
    if (loBytes) std::memcpy(loSlice.data(), neighLo + loSrcByte0, loBytes);
}

/**
 * decode_single_entry_from_slices
 * -------------------------------
 * Decode exactly one entry using local slices and intra-slice bit offsets.
 *
 * Inputs:
 *   deg_blk[0], bit_blk[0], huff_blk[0] : per-entry metadata for this entry
 *   hiSlice/loSlice                      : local byte windows
 *   hiBitOffset/loBitOffset              : starting bit positions within slices
 *   treeOpp, fallbackBitsOpp             : Huffman tree + fixed-width fallback size
 *
 * Output:
 *   Vector of decoded integers (neighbors).
 *
 * Steps:
 *   1) Decode h = huff_blk[0] symbols from hiSlice (bounded by hb = bit_blk[0] bits).
 *   2) Skip to the end of hb bits (consume padding if any).
 *   3) Decode the remaining (deg - h) values from loSlice using fixed width.
 */
std::vector<int> decode_single_entry_from_slices(
    const int*       deg_blk,        // length >= 1
    const uint16_t*  bit_blk,        // length >= 1
    const uint16_t*  huff_blk,       // length >= 1
    const uint8_t*   hiSlice,
    const uint8_t*   loSlice,
    uint32_t         hiBitOffset,
    uint32_t         loBitOffset,
    HuffmanNode*     treeOpp,
    int              fallbackBitsOpp
) {
    std::vector<int> out;
    const int d  = deg_blk[0];
    if (d <= 0) return out;

    out.reserve(d);

    const int h  = huff_blk[0];
    const int hb = bit_blk[0];

    uint64_t ptrHi_bits = hiBitOffset;  // bit pointer into hiSlice
    uint64_t ptrLo_bits = loBitOffset;  // bit pointer into loSlice

    // Huffman section: decode up to h symbols within hb bits
    const uint64_t hiStart = ptrHi_bits;
    int usedH   = 0;
    int usedBit = 0;
    while (usedH < h && usedBit < hb) {
        HuffmanNode* cur = treeOpp;
        // Walk down to a leaf
        do {
            bool bit = (hiSlice[ptrHi_bits >> 3] >> (ptrHi_bits & 7)) & 1;
            cur = bit ? cur->r : cur->l;
            ++ptrHi_bits; ++usedBit;
        } while (cur->l || cur->r);
        out.push_back(cur->v);
        ++usedH;
    }
    // Align to hb bits (skip any remainder/padding)
    const uint64_t shouldBe = hiStart + (uint64_t)hb;
    if (ptrHi_bits < shouldBe) ptrHi_bits = shouldBe;

    // Fallback section: fixed-width decode for the rest
    for (int t = usedH; t < d; ++t) {
        int x = readBits64(loSlice, ptrLo_bits, fallbackBitsOpp);
        out.push_back(x);
    }

    return out;
}

/**
 * decode_block_partial
 * --------------------
 * Decode entries [i0, i1) using pre-extracted slices and per-entry metadata.
 *
 * Inputs:
 *   deg_blk/bit_blk/huff_blk  : metadata slices for [i0, i1)
 *   hiSlice/loSlice           : byte windows covering the whole [i0, i1) region
 *   hiBitOffset/loBitOffset   : bit positions for entry i0 inside the slices
 *   treeOpp, fallbackBitsOpp  : decoding parameters
 *   edges (optional)          : if non-null, used for sum-check verification
 *   ePtr, totalError          : running verification counters (updated if edges != nullptr)
 *
 * Side effects:
 *   adj[i] is filled with decoded neighbors for every i in [i0, i1).
 *   ePtr/totalError advanced if verifying.
 */
void decode_block_partial(
    int i0, int i1,
    const int*       deg_blk,
    const uint16_t*  bit_blk,
    const uint16_t*  huff_blk,
    const uint8_t*   hiSlice,
    const uint8_t*   loSlice,
    uint32_t         hiBitOffset,
    uint32_t         loBitOffset,
    HuffmanNode*     treeOpp,
    int              fallbackBitsOpp,
    std::vector<std::vector<int>>& adj,
    const int*       edges,      // may be nullptr
    long&            ePtr,       // running pointer into edges (global)
    long&            totalError  // accumulated sum-check difference (global)
) {
    if (i0 >= i1) return;

    // Local bit cursors within the slices; advance per entry
    uint64_t curHiBits = hiBitOffset;
    uint64_t curLoBits = loBitOffset;

    for (int i = i0; i < i1; ++i) {
        const int k  = i - i0;

        // Decode one entry using local slices (no global bitstream access)
        std::vector<int> out = decode_single_entry_from_slices(
            /*deg_blk=*/      deg_blk  + k,
            /*bit_blk=*/      bit_blk  + k,
            /*huff_blk=*/     huff_blk + k,
            /*hiSlice=*/      hiSlice,
            /*loSlice=*/      loSlice,
            /*hiBitOffset=*/  static_cast<uint32_t>(curHiBits),
            /*loBitOffset=*/  static_cast<uint32_t>(curLoBits),
            /*treeOpp=*/      treeOpp,
            /*fallbackBits=*/ fallbackBitsOpp
        );

        adj[i] = std::move(out);

        // Advance to the next entry’s starting bits within the same slices
        const int d  = deg_blk[k];
        const int h  = huff_blk[k];
        const int hb = bit_blk[k];
        curHiBits += static_cast<uint64_t>(hb);
        curLoBits += static_cast<uint64_t>(d - h) * static_cast<uint64_t>(fallbackBitsOpp);

        // Optional verification (sum-check)
        if (edges) {
            long sumOut = 0, sumOrig = 0;
            for (int j = 0; j < d; ++j, ++ePtr) {
                sumOut  += adj[i][j];
                sumOrig += edges[ePtr];
            }
            totalError += (sumOut - sumOrig);
        }
    }
}

/**
 * decodeSide
 * ----------
 * Top-level paged decoder.
 * Splits [0, n) into blocks of size BLOCK, extracts per-block slices, then decodes.
 * If 'edges' is provided, computes a sum-check for the decoded neighborhood lists.
 *
 * Returns:
 *   adj: vector of length n where adj[i] holds decoded neighbors of entry i.
 */
std::vector<std::vector<int>> decodeSide(
    int n,
    const int *deg,
    const int *edges,                 // may be nullptr to disable verification
    const uint16_t *huffCount,
    const uint16_t *bitCount,
    const uint8_t *neighHi,
    const uint8_t *neighLo,
    HuffmanNode* treeOpp,
    int fallbackBitsOpp,
    const uint64_t* blockHi,
    const uint64_t* blockLo,
    int BLOCK
) {
    std::vector<std::vector<int>> adj(n);
    if (n <= 0) return adj;

    const int nBlocks = (n + BLOCK - 1) / BLOCK;
    long ePtr = 0;        // cursor for verification against edges[]
    long totalError = 0;  // accumulated sum-check difference

    for (int b = 0; b < nBlocks; ++b) {
        const int i0 = b * BLOCK;
        const int i1 = std::min(n, i0 + BLOCK);

        std::vector<int>       deg_blk;
        std::vector<uint16_t>  bit_blk, huff_blk;
        std::vector<uint8_t>   hiSlice, loSlice;
        uint32_t hiBitOffset = 0, loBitOffset = 0;

        // Copy minimal metadata/slices for this block range
        extract_block_slices(
            i0, i1, BLOCK,
            deg, bitCount, huffCount,
            neighHi, neighLo,
            blockHi, blockLo,
            fallbackBitsOpp,
            deg_blk, bit_blk, huff_blk,
            hiSlice, loSlice,
            hiBitOffset, loBitOffset
        );

        // Decode entries covered by this block’s slices
        decode_block_partial(
            i0, i1,
            deg_blk.data(), bit_blk.data(), huff_blk.data(),
            hiSlice.data(), loSlice.data(),
            hiBitOffset, loBitOffset,
            treeOpp, fallbackBitsOpp,
            adj,
            edges, ePtr, totalError
        );
    }

    if (edges) {
        if (totalError == 0)
            std::cout << "decodeSide(paged): verified by sum-check.\n";
        else
            std::cout << "decodeSide(paged): totalError=" << totalError << "\n";
    }
    return adj;
}
