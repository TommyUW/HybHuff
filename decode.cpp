// ===== File: decode.cpp =====
//
// Description:
// This function performs hybrid decoding of an encoded hypergraph representation.
// The input is compressed using a two-part strategy:
//   - A Huffman tree is used to decode frequent values (encoded in `neighHi`).
//   - A fixed-width fallback encoding is used for the remaining less frequent values (`neighLo`).
// The function decodes the compressed format and reconstructs the original adjacency lists
// (stored in a 2D vector `adj`). It also verifies correctness via a sum-check against the
// original values (`edges[]`), printing diagnostics in case of mismatch.
//
// Inputs:
//   - n: number of vertices (or hyperedges) to decode
//   - deg: degree of each vertex (or size of each hyperedge)
//   - edges: original flat array of values for comparison
//   - bitCount: how many bits are allocated to Huffman codes for each entry
//   - neighHi: bit-packed Huffman-encoded data
//   - neighLo: bit-packed fallback data (fixed-width encoding)
//   - treeOpp: root of the Huffman decoding tree
//   - fallbackBitsOpp: number of bits per value in fallback encoding
//
// Output:
//   - A vector of vectors (adjacency list), where each entry `adj[i]` holds the decoded values
//     for vertex (or hyperedge) `i`.
//
// =======================================

#include "decode.h"
#include <iostream>
#include <vector>
#include <cassert>

// Helper function to read `n` bits from array A starting at bit position `pos` (big-endian)
inline int readBits(const uint8_t *A, long &pos, int n) {
    int x = 0;
    for (int i = 0; i < n; i++) {
        // Extract one bit from byte array A at bit offset pos
        x = (x << 1) | ((A[pos / 8] >> (pos % 8)) & 1);
        pos++;
    }
    return x;
}

// Main decode function for one side of the bipartite hypergraph
std::vector<std::vector<int>> decodeSide(
    int n,
    const int *deg,
    const int *edges,
    const int *bitCount,
    const uint8_t *neighHi,
    const uint8_t *neighLo,
    HuffmanNode* treeOpp,
    int fallbackBitsOpp
) {
    std::vector<std::vector<int>> adj(n);  // Output adjacency list
    adj.reserve(n);

    long ptrHi = 0;    // Bit pointer into Huffman bitstream
    long ptrLo = 0;    // Bit pointer into fallback bitstream
    long ePtr  = 0;    // Pointer into original edges array

    long totalError = 0;  // Used to accumulate total sum-check differences

    for (int i = 0; i < n; ++i) {
        int d = deg[i];  // Number of elements to decode for this entry
        std::vector<int> out;
        out.reserve(d);

        // Step 1: Decode Huffman-encoded values using bitCount[i] bits
        int collected = 0;
        for (int b = 0; b < bitCount[i]; ) {
            HuffmanNode *cur = treeOpp;

            // Traverse the Huffman tree bit by bit until a leaf is reached
            while (cur->l || cur->r) {
                if (ptrHi < 0) {
                    std::cerr << "ERROR: negative ptrHi\n";
                    break;
                }

                // Read a single bit from neighHi
                bool bit = ((neighHi[ptrHi / 8] >> (ptrHi % 8)) & 1);
                cur = bit ? cur->r : cur->l;
                ++ptrHi;
                ++b;
            }

            // Store the decoded value from the Huffman leaf
            out.push_back(cur->v);
            ++collected;
        }

        // Step 2: Decode any remaining values using fixed-width fallback encoding
        while (collected < d) {
            int x = readBits(neighLo, ptrLo, fallbackBitsOpp);
            out.push_back(x);
            ++collected;
        }

        // Step 3: Perform sum-check to compare with original edges[]
        long sumOut  = 0;
        long sumOrig = 0;
        for (int idx = 0; idx < d; idx++, ++ePtr) {
            sumOut  += out[idx];
            sumOrig += edges[ePtr];
        }

        long diff = sumOut - sumOrig;
        totalError += diff;

        // Save decoded values into output adjacency list
        adj[i] = std::move(out);
    }

    // Final report on decoding accuracy
    if (totalError == 0) {
        std::cout << "decodeSide: all " << n << " entries verified by sum‐check\n";
    } else {
        std::cout << "decodeSide: totalError=" << totalError << "\n";
    }

    return adj;
}
