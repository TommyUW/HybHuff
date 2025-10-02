// ===== File: decode.h =====
//
// Description
// -----------
// Interface for the block-paged hybrid decoder used in compressed hypergraph/graph
// processing. A two-tier format stores each entry’s neighbors as:
//   • neighHi : Huffman-encoded symbols for frequent values
//   • neighLo : fixed-width (fallback) values for the remainder
//
// This interface supports decoding individual entries from pre-extracted slices
// as well as a top-level, block-paged decode of an entire side. Optionally, a
// lightweight sum-check can verify decoded results against the original edges[].
//
// Usage
// -----
// Include this header to call `decodeSide(...)` when reconstructing adjacency
// lists for algorithms like BFS, k-core, PageRank, etc.
//
// Dependencies
// ------------
// • Requires Huffman tree types from "huffman_tree.h"
//
#ifndef DECODE_H
#define DECODE_H

#include <cstdint>
#include <vector>
#include "huffman_tree.h"
using namespace std;  // (Header-wide 'using' is generally discouraged; kept for compatibility.)

// -----------------------------------------------------------------------------
// decodeSide
// -----------------------------------------------------------------------------
// Top-level paged decoder.
// Splits [0, n) into blocks of size BLOCK. For each block, it extracts minimal
// metadata and bit windows (from neighHi/neighLo) using blockHi/blockLo start
// offsets, then decodes entries into adj[i]. If edges != nullptr, performs a
// sum-check (running difference of decoded vs. original values).
//
// Parameters:
//   n                 : number of entries to decode (vertices or hyperedges)
//   deg[n]            : per-entry degrees
//   edges[sum(deg)]   : original flat adjacency for optional verification
//   huffCount[n]      : number of Huffman symbols per entry
//   bitCount[n]       : number of bits occupied by the Huffman section per entry
//   neighHi, neighLo  : global bitstreams for Huffman/fallback sections
//   treeOpp           : Huffman tree used for decoding the Huffman section
//   fallbackBitsOpp   : bit width for each fallback value
//   blockHi, blockLo  : bit offsets at block starts (entry k*BLOCK)
//   BLOCK             : block size (in entries)
//
// Returns:
//   Vector< Vector<int> > of length n with decoded neighbors per entry.
//
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
    const uint64_t* blockHi,          // start bit offset in 'hi' per block
    const uint64_t* blockLo,          // start bit offset in 'lo' per block
    int BLOCK                          // block size used during encoding/decoding
);

// -----------------------------------------------------------------------------
// extract_block_slices
// -----------------------------------------------------------------------------
// Build minimal per-block views for entries [i0, i1):
//   • deg_blk/bit_blk/huff_blk : slices of metadata for that range
//   • hiSlice/loSlice          : byte windows copied from the global bitstreams
//   • hiBitOffset/loBitOffset  : starting bit positions inside those slices for i0
//
// Notes:
//   • The pair (hiSlice, hiBitOffset) and (loSlice, loBitOffset) define where
//     local decoding of entry i0 begins inside the copied windows.
//   • The caller advances bit pointers between entries using bit_blk/huff_blk/deg_blk.
//
void extract_block_slices(
    int i0, int i1,
    int BLOCK,
    const int*       deg,
    const uint16_t*  bitCount,
    const uint16_t*  huffCount,
    const uint8_t*   neighHi,
    const uint8_t*   neighLo,
    const uint64_t*  blockHi,       // bit offsets at block starts (hi)
    const uint64_t*  blockLo,       // bit offsets at block starts (lo)
    int              fallbackBitsOpp,
    // outputs
    std::vector<int>&       deg_blk,
    std::vector<uint16_t>&  bit_blk,
    std::vector<uint16_t>&  huff_blk,
    std::vector<uint8_t>&   hiSlice,
    std::vector<uint8_t>&   loSlice,
    uint32_t&               hiBitOffset,  // starting bit inside first byte of hiSlice
    uint32_t&               loBitOffset   // starting bit inside first byte of loSlice
);

// -----------------------------------------------------------------------------
// decode_block_partial
// -----------------------------------------------------------------------------
// Decode entries [i0, i1) using pre-extracted slices and per-entry metadata.
//
// Parameters:
//   deg_blk/bit_blk/huff_blk : metadata slices for [i0, i1)
//   hiSlice/loSlice          : byte windows covering [i0, i1)
//   hiBitOffset/loBitOffset  : bit positions for entry i0 within those slices
//   treeOpp/fallbackBitsOpp  : decoding parameters
//   adj                      : destination adjacency vectors (size >= n)
//   edges/ePtr/totalError    : optional verification state; updated if edges != nullptr
//
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
);

// -----------------------------------------------------------------------------
// decode_single_entry_from_slices
// -----------------------------------------------------------------------------
// Decode exactly one entry (deg_blk[0]) from local slices, starting from the given
// intra-slice bit offsets. First consume huff_blk[0] symbols within bit_blk[0] bits
// from hiSlice, then decode the remaining (deg - huff) values from loSlice using
// fixed-width fallbackBitsOpp.
//
std::vector<int> decode_single_entry_from_slices(
    const int*       deg_blk,
    const uint16_t*  bit_blk,
    const uint16_t*  huff_blk,
    const uint8_t*   hiSlice,
    const uint8_t*   loSlice,
    uint32_t         hiBitOffset,
    uint32_t         loBitOffset,
    HuffmanNode*     treeOpp,
    int              fallbackBitsOpp
);

#endif  // DECODE_H
