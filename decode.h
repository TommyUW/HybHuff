// ===== File: decode.h =====
//
// Description:
// Interface for the hybrid decoding routine used in compressed hypergraph processing.
// This decoder reconstructs the adjacency list of vertices or hyperedges from a two-tier
// encoded format:
//   - `neighHi`: Encoded using Huffman tree compression for frequent elements.
//   - `neighLo`: Encoded using fixed-width fallback encoding for rare elements.
// It also performs a checksum comparison (sum-check) against the original uncompressed
// values to ensure decoding correctness.
//
// Usage:
// Include this header when calling `decodeSide()` to reconstruct the decoded graph
// for tasks such as BFS or other graph traversals.
//
// Dependencies:
// - Requires the Huffman tree structure defined in "huffman_tree.h"
//

#ifndef DECODE_H
#define DECODE_H

#include <cstdint>
#include <vector>
#include "huffman_tree.h"

// Function: decodeSide
//
// Parameters:
//   - n:              Number of vertices (or hyperedges) to decode
//   - deg:            Array of size n, degree of each vertex (or size of each hyperedge)
//   - edges:          Flat array of original uncompressed neighbor values (for verification)
//   - bitCount:       Array of size n, total Huffman bits assigned to each entry
//   - neighHi:        Bitstream storing Huffman-encoded values (packed into bytes)
//   - neighLo:        Bitstream storing fallback-encoded values (fixed width per value)
//   - treeOpp:        Root pointer to Huffman tree used for decoding
//   - fallbackBitsOpp:Number of bits per fallback value (typically log2 of range)
//
// Returns:
//   A vector of n vectors, where each inner vector contains the decoded neighbors
//   of the corresponding vertex/hyperedge.
//
std::vector<std::vector<int>> decodeSide(
    int n,
    const int *deg,
    const int *edges,
    const int *bitCount,
    const uint8_t *neighHi,
    const uint8_t *neighLo,
    HuffmanNode* treeOpp,
    int fallbackBitsOpp
);

#endif  // DECODE_H
