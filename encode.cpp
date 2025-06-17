// ===== File: encode.cpp =====
//
// Description:
// This file implements the hybrid encoding routine for hypergraph compression.
// Each neighbor (or edge) is either encoded using:
//   - Huffman coding (via `codesOpp`) for high-frequency values
//   - Fixed-width binary encoding (`fallbackBitsOpp`) for others
//
// Two separate bitstreams are produced:
//   - `neighHi`: Huffman-encoded bits (packed into a byte array)
//   - `neighLo`: Fallback-encoded bits (packed into a byte array)
//
// The function also computes and returns:
//   - `bitCount[i]`: number of bits used by Huffman codes for entry `i`
//   - `(bitsH, bitsL)`: total bits used for Huffman and fallback encoding respectively
//
// Inputs:
//   - n:            Number of vertices (or hyperedges) to encode
//   - deg:          Degree of each entry (array of size n)
//   - edges:        Flat array of neighbor indices
//   - codesOpp:     Huffman code table mapping values to binary strings
//   - fallbackBitsOpp: Number of bits for fixed-width fallback encoding
//
// Outputs (allocated within function):
//   - bitCount:     Array of size n, storing total Huffman bits used per entry
//   - neighHi:      Byte array of packed Huffman-encoded bits
//   - neighLo:      Byte array of packed fallback-encoded bits
//
// Return value:
//   - Pair (bitsH, bitsL): total number of bits in neighHi and neighLo respectively
//

#include "encode.h"
#include <vector>
#include <unordered_map>
#include <sys/resource.h>
#include <unistd.h>
#include <chrono>
#include <malloc.h>

using namespace std;

pair<long, long> encodeSide(
    int n,
    const int *deg,
    const int *edges,
    const unordered_map<int, std::string> &codesOpp,
    int fallbackBitsOpp,
    int *&bitCount,
    uint8_t *&neighHi,
    uint8_t *&neighLo
) {
    vector<bool> HF; // Huffman bit stream (as a vector of bits)
    vector<bool> LF; // Fallback bit stream (as a vector of bits)

    bitCount = new int[n]();  // Allocate and initialize to 0
    long ptr = 0;             // Index into the flat `edges` array
    long bitsH = 0, bitsL = 0;

    // Loop through each entry
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < deg[i]; j++, ptr++) {
            int nei = edges[ptr];

            // Try Huffman encoding first
            auto it = codesOpp.find(nei);
            if (it != codesOpp.end()) {
                // Encode using Huffman (push each bit)
                for (char c : it->second) {
                    HF.push_back(c == '1');
                    bitsH++;
                }
                bitCount[i] += it->second.size();  // Track bits used
            } else {
                // Fallback to fixed-width binary encoding
                int x = nei;
                for (int b = fallbackBitsOpp - 1; b >= 0; b--) {
                    LF.push_back((x >> b) & 1);  // Extract each bit (MSB-first)
                    bitsL++;
                }
            }
        }
    }

    // Compute number of bytes needed for each stream (round up to byte boundary)
    long bytesH = (bitsH + 7) / 8;
    long bytesL = (bitsL + 7) / 8;

    // Allocate output byte arrays and initialize to 0
    neighHi = bytesH ? new uint8_t[bytesH]() : nullptr;
    neighLo = bytesL ? new uint8_t[bytesL]() : nullptr;

    // Pack bits into bytes for Huffman stream
    for (long i = 0; i < bitsH; i++) {
        if (HF[i])
            neighHi[i / 8] |= 1 << (i % 8);  // Set the i-th bit in the byte array
    }

    // Pack bits into bytes for fallback stream
    for (long i = 0; i < bitsL; i++) {
        if (LF[i])
            neighLo[i / 8] |= 1 << (i % 8);
    }

    return {bitsH, bitsL};  // Return total bits in each stream
}
