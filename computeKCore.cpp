#include "computeKCore.h"
#include <queue>
#include <iostream>
#include <chrono>
#include <sys/resource.h>
#include <unistd.h>
#include <vector>
#include <algorithm>
#include <omp.h>
#include <cstdint>
#include "decode.h"
#include "huffman_tree.h"
#include "misc.h"
#include "malloc.h"
using namespace std;

// --------------------------------------------------------------
// K-core on the encoded representation (on-demand decoding).
// Strategy:
//   • Maintain mutable degrees (curDeg) and a removal frontier R.
//   • Iteratively peel vertices with degree < k_thresh.
//   • For each vertex v to peel, decode only v’s neighbors from
//     compressed slices extracted for [v, v+1).
//   • Updates to curDeg are atomic; next round candidates deduped.
// Inputs:
//   n                 : number of vertices
//   deg               : original degrees (encoded-side semantics)
//   huffCount,bitCount: per-entry symbol/bit counts for hi/lo streams
//   neighHi,neighLo   : compressed bitstreams
//   treeOpp           : Huffman tree for opposite side
//   fallbackBitsOpp   : fixed width for fallback-coded items
//   blockHi,blockLo   : per-block base offsets into hi/lo
//   BLOCK             : vertex block size
//   k_thresh          : K threshold
// Output:
//   inCore mask (1 = remains in k-core, 0 = peeled)
// --------------------------------------------------------------
std::vector<char> computeKCore_onDemand_decodeRandom(
    int n,
    const int* deg,
    const uint16_t* huffCount,
    const uint16_t* bitCount,
    const uint8_t* neighHi,
    const uint8_t* neighLo,
    HuffmanNode* treeOpp,
    int fallbackBitsOpp,
    // block index built at encode-time
    const uint64_t* blockHi,
    const uint64_t* blockLo,
    int BLOCK,
    int k_thresh
) {
    using clock = std::chrono::high_resolution_clock;
    auto t0 = clock::now();

    if (n <= 0) return {};

    // Working degree, removal flags, and initial frontier
    std::vector<int>  curDeg(n);
    std::vector<char> removed(n, 0);
    std::vector<int>  R; R.reserve(1024);

    for (int v = 0; v < n; ++v) {
        curDeg[v] = deg[v];
        if (curDeg[v] < k_thresh) { removed[v] = 1; R.push_back(v); }
    }

    // Placeholder for any per-slice metric if tracked externally (unused)
    size_t peak_slice_rss_bytes = 0;

    // Peeling rounds
    while (!R.empty()) {
        std::vector<int> nextR;

        #pragma omp parallel
        {
            // Thread-local buffers
            std::vector<int> localNext; localNext.reserve(256);
            std::vector<int>       deg_blk;
            std::vector<uint16_t>  bit_blk, huff_blk;
            std::vector<uint8_t>   hiSlice, loSlice;

            // Each thread extracts the minimal slice needed for vertex v
            #pragma omp for schedule(dynamic, 1024)
            for (int idx = 0; idx < (int)R.size(); ++idx) {
                const int v = R[idx];
                if (deg[v] == 0) continue;  // isolated

                // Extract contiguous slices for [v, v+1)
                uint32_t hiBitOffset = 0, loBitOffset = 0;

                // Guarded extraction; also ensures consistent shared buffers if any
                #pragma omp critical(__kcore_block_rss)
                {
                    extract_block_slices(
                        /*i0=*/v, /*i1=*/v+1,
                        BLOCK,
                        deg, bitCount, huffCount,
                        neighHi, neighLo,
                        blockHi, blockLo,
                        fallbackBitsOpp,
                        deg_blk, bit_blk, huff_blk,
                        hiSlice, loSlice,
                        hiBitOffset, loBitOffset
                    );
                }

                // Decode exactly one entry (v) from the local slices
                std::vector<int> nbrs = decode_single_entry_from_slices(
                    deg_blk.data(), bit_blk.data(), huff_blk.data(),
                    hiSlice.data(),  loSlice.data(),
                    hiBitOffset,     loBitOffset,
                    treeOpp, fallbackBitsOpp
                );

                // Decrement degrees of neighbors; enqueue threshold-crossers
                for (int u : nbrs) {
                    if ((unsigned)u >= (unsigned)n) continue;
                    if (removed[u]) continue;

                    int old;
                    #pragma omp atomic capture
                    { old = curDeg[u]; curDeg[u] = old - 1; }

                    if (old == k_thresh) {
                        localNext.push_back(u);
                    }
                }
            } // omp for

            // Merge thread-local candidates
            #pragma omp critical
            nextR.insert(nextR.end(), localNext.begin(), localNext.end());
        } // omp parallel

        if (nextR.empty()) break;

        // Deduplicate and finalize removals
        std::sort(nextR.begin(), nextR.end());
        nextR.erase(std::unique(nextR.begin(), nextR.end()), nextR.end());

        int w = 0;
        for (int u : nextR) {
            if (!removed[u]) {
                removed[u] = 1;
                nextR[w++] = u;
            }
        }
        nextR.resize(w);

        R.swap(nextR);
    }

    // Build in-core mask
    std::vector<char> inCore(n, 0);
    int remain = 0;
    for (int v = 0; v < n; ++v) { inCore[v] = !removed[v]; remain += inCore[v]; }

    auto t1 = clock::now();
    std::cout << "KCore(Hybuff) time: " << std::chrono::duration<double>(t1 - t0).count() << " s\n";
    return inCore;
}

// --------------------------------------------------------------
// K-core on the raw (uncompressed) adjacency.
// Strategy:
//   • Build CSR offsets from deg[] for O(1) neighbor range lookup.
//   • Standard peeling loop identical to encoded version, but
//     neighbors come directly from edges[] ranges.
// Inputs:
//   n, deg[], edges[], k_thresh
// Output:
//   inCore mask (1 = remains in k-core, 0 = peeled)
// --------------------------------------------------------------
std::vector<char> computeKCore_raw(
    int n,
    const int* deg,
    const int* edges,
    int k_thresh
) {
    using clock = std::chrono::high_resolution_clock;
    auto t0 = clock::now();

    if (n <= 0) return {};

    // Optional process-metric placeholder (unused)
    size_t r_entry = rss_bytes();

    // CSR offsets to index neighbors quickly
    std::vector<uint64_t> off(n + 1, 0);
    for (int v = 0; v < n; ++v) off[v + 1] = off[v] + (uint64_t)deg[v];

    std::vector<int>  curDeg(n);
    std::vector<char> removed(n, 0);
    for (int v = 0; v < n; ++v) curDeg[v] = deg[v];

    // Initial removals
    std::vector<int> R; R.reserve(1024);
    for (int v = 0; v < n; ++v) {
        if (curDeg[v] < k_thresh) { removed[v] = 1; R.push_back(v); }
    }

    // Peeling rounds
    while (!R.empty()) {
        std::vector<int> nextR;

        #pragma omp parallel
        {
            std::vector<int> localNext;
            localNext.reserve(256);

            #pragma omp for schedule(dynamic, 1024)
            for (int i = 0; i < (int)R.size(); ++i) {
                int v = R[i];
                if (deg[v] == 0) continue;

                const uint64_t s = off[v];
                const uint64_t e = off[v + 1];

                for (uint64_t p = s; p < e; ++p) {
                    int u = edges[p];
                    if ((unsigned)u >= (unsigned)n) continue;
                    if (removed[u]) continue;

                    int old;
                    #pragma omp atomic capture
                    { old = curDeg[u]; curDeg[u] = old - 1; }

                    if (old == k_thresh) localNext.push_back(u);
                }
            }

            // Merge thread-local candidates
            #pragma omp critical
            nextR.insert(nextR.end(), localNext.begin(), localNext.end());
        } // parallel

        if (nextR.empty()) break;

        // Deduplicate and finalize removals
        std::sort(nextR.begin(), nextR.end());
        nextR.erase(std::unique(nextR.begin(), nextR.end()), nextR.end());

        int w = 0;
        for (int u : nextR) {
            if (!removed[u]) { removed[u] = 1; nextR[w++] = u; }
        }
        nextR.resize(w);

        R.swap(nextR);
    }

    // Build in-core mask
    std::vector<char> inCore(n, 0);
    int remain = 0;
    for (int v = 0; v < n; ++v) { inCore[v] = !removed[v]; remain += inCore[v]; }

    auto t1 = clock::now();
    std::cout << "KCore(raw) time: " << std::chrono::duration<double>(t1 - t0).count() << " s\n";
    return inCore;
}
