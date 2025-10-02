// ==========================================================
// File: BFS_block_memprofile.cpp
// Purpose: Level-synchronous single-source BFS for two modes:
//   (1) Encoded graph with block-indexed, on-demand decoding.
//   (2) Raw CSR baseline (deg[] + edges[]).
// Notes:
//   • Encoded mode caches one block per thread and decodes a
//     single vertex’s neighbors from local slices.
//   • Raw mode builds CSR offsets and performs a standard BFS.
// ==========================================================

#include "BFS.h"
#include <iostream>
#include <queue>
#include <chrono>
#include <vector>
#include <climits>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <algorithm>
#include <atomic>
#include <unistd.h>       // sysconf
#include <sys/resource.h>
#include "decode.h"          // extract_block_slices(...), decode_single_entry_from_slices(...)
#include "huffman_tree.h"    // struct HuffmanNode { HuffmanNode *l,*r; int v; }
#include "misc.h"
#include <malloc.h>
using namespace std;

// =============================================================
// Encoded BFS (block-indexed, on-demand decode; parallel-friendly).
// Inputs:
//   n                : #vertices
//   deg              : per-vertex degree (encoded-side semantics)
//   edges            : unused in encoded mode (kept for symmetry)
//   huffCount        : #Huffman-coded items per vertex in block
//   bitCount         : #bits contributed by Huffman-coded items
//   hi, lo           : compressed bitstreams
//   treeOpp          : Huffman tree (opposite side)
//   fallbackBitsOpp  : fixed width for fallback-coded items
//   src              : BFS source
//   blockHi, blockLo : per-block base bit offsets
//   BLOCK            : block size (vertex granularity)
// Behavior:
//   • Processes current frontier; per thread, lazily loads the
//     block containing u and decodes neighbors of u only.
//   • Thread-local vectors avoid cross-thread contention.
//   • Next frontier is deduplicated before the next level.
// =============================================================
void runBFSFromSingleSource(
    int n, const int* deg, const int* edges,
    const uint16_t* huffCount, const uint16_t* bitCount,
    const uint8_t* hi, const uint8_t* lo,
    HuffmanNode* treeOpp, int fallbackBitsOpp,
    int src,
    const uint64_t* blockHi, const uint64_t* blockLo, int BLOCK
) {
    using clock = std::chrono::high_resolution_clock;
    auto t0 = clock::now();

    // Basic input checks (fail fast).
    if (n <= 0 || src < 0 || src >= n) { std::cout << "Single-source BFS: invalid n/src\n"; return; }
    if (!deg || !huffCount || !bitCount || !hi || !lo || !treeOpp ||
        !blockHi || !blockLo || BLOCK <= 0) {
        std::cerr << "Single-source BFS: missing/invalid inputs\n";
        return;
    }

    const int INF = INT_MAX/4;
    std::vector<int> dist(n, INF);
    dist[src] = 0;

    // Level-synchronous frontier
    std::vector<int> frontier, nextFrontier;
    frontier.push_back(src);
    int maxLevel = 0, reached = 1;

    // Peak block metric placeholder (intentionally unused here).
    size_t peak_block_rss_bytes = 0;

    while (!frontier.empty()) {
        const int nextLevel = maxLevel + 1;
        nextFrontier.clear();

        // Parallelize across current frontier.
        #pragma omp parallel
        {
            // -------- Thread-local block cache --------
            // The following buffers/slices are re-populated when
            // the thread switches to a different block.
            int t_curBlk = -1, t_blk_i0 = 0, t_blk_i1 = 0;
            std::vector<int>       t_deg_blk;
            std::vector<uint16_t>  t_bit_blk, t_huff_blk;
            std::vector<uint8_t>   t_hiSlice, t_loSlice;
            uint32_t               t_baseHiOffset = 0, t_baseLoOffset = 0;
            std::vector<uint64_t>  t_offHiBits, t_offLoBits;

            // Accumulate thread-local discoveries and merge once.
            std::vector<int> t_localNext;  t_localNext.reserve(1024);

            // Load/calc slices for block 'blk' if not already active.
            // Populates per-entry offset arrays used for decode.
            auto ensure_block_loaded = [&](int blk) {
                if (blk == t_curBlk) return;
                t_curBlk = blk;
                t_blk_i0 = blk * BLOCK;
                t_blk_i1 = std::min(n, t_blk_i0 + BLOCK);
                t_baseHiOffset = t_baseLoOffset = 0;

                // Extract contiguous slices for [t_blk_i0, t_blk_i1).
                #pragma omp critical(__bfs_block_rss)
                {
                    extract_block_slices(
                        t_blk_i0, t_blk_i1,
                        BLOCK,
                        deg, bitCount, huffCount,
                        hi, lo,
                        blockHi, blockLo,
                        fallbackBitsOpp,
                        t_deg_blk, t_bit_blk, t_huff_blk,
                        t_hiSlice, t_loSlice,
                        t_baseHiOffset, t_baseLoOffset
                    );
                }

                // Precompute per-entry bit offsets within the slices.
                const int L = t_blk_i1 - t_blk_i0;
                t_offHiBits.assign(L, 0);
                t_offLoBits.assign(L, 0);
                uint64_t accH = 0, accL = 0;
                for (int k = 0; k < L; ++k) {
                    t_offHiBits[k] = accH;
                    t_offLoBits[k] = accL;
                    accH += (uint64_t)t_bit_blk[k];
                    accL += (uint64_t)(t_deg_blk[k] - t_huff_blk[k]) * (uint64_t)fallbackBitsOpp;
                }
            };

            // Work over the current frontier; each u may trigger a block load.
            #pragma omp for schedule(dynamic, 1024)
            for (int i = 0; i < (int)frontier.size(); ++i) {
                int u = frontier[i];
                if (deg[u] == 0) continue; // isolated vertex

                const int blk = u / BLOCK;
                ensure_block_loaded(blk);

                // Local index inside the active block.
                const int k = u - t_blk_i0;
                const uint32_t hiOff = t_baseHiOffset + (uint32_t)t_offHiBits[k];
                const uint32_t loOff = t_baseLoOffset + (uint32_t)t_offLoBits[k];

                // Decode neighbors of u from thread-local slices only.
                std::vector<int> nbrs = decode_single_entry_from_slices(
                    &t_deg_blk[k], &t_bit_blk[k], &t_huff_blk[k],
                    t_hiSlice.data(), t_loSlice.data(),
                    hiOff, loOff,
                    treeOpp, fallbackBitsOpp
                );

                // Standard BFS relax.
                for (int v : nbrs) {
                    if ((unsigned)v >= (unsigned)n) continue;
                    if (dist[v] != INF) continue;
                    #pragma omp critical(__bfs_visit_dist)
                    {
                        if (dist[v] == INF) {
                            dist[v] = nextLevel;
                            t_localNext.push_back(v);
                        }
                    }
                }
            } // frontier loop

            // Merge thread-local next frontier.
            #pragma omp critical(__bfs_merge_next)
            {
                reached += (int)t_localNext.size();
                nextFrontier.insert(nextFrontier.end(),
                                    t_localNext.begin(), t_localNext.end());
            }
        } // parallel region

        // Remove duplicates before proceeding to the next level.
        std::sort(nextFrontier.begin(), nextFrontier.end());
        nextFrontier.erase(std::unique(nextFrontier.begin(), nextFrontier.end()), nextFrontier.end());

        frontier.swap(nextFrontier);
        if (!frontier.empty()) maxLevel = nextLevel;
    }

    auto t1 = clock::now();
    std::cout << "BFS(Hybhuff) time: "
              << std::chrono::duration<double>(t1 - t0).count() << " s\n";
}

// =============================================================
// Raw CSR BFS (no decoding).
// Inputs:
//   n, deg[], edges[], src
// Behavior:
//   • Builds CSR offsets from deg[].
//   • Standard queue-based BFS over edges[].
//   • Returns distance array (INF for unreachable).
// =============================================================
std::vector<int> bfs_single_source_raw(
    int n,
    const int* deg,
    const int* edges,
    int src
) {
    using clock = std::chrono::high_resolution_clock;
    auto t0 = clock::now();

    std::vector<int> dist(n, INT_MAX/4);
    if (n <= 0 || src < 0 || src >= n) {
        return dist;
    }

    // CSR offsets from degrees.
    std::vector<long long> off(n + 1, 0);
    for (int i = 0; i < n; ++i) off[i + 1] = off[i] + (long long)deg[i];

    std::queue<int> q;
    dist[src] = 0;
    q.push(src);

    int reached = 1, maxLevel = 0;

    while (!q.empty()) {
        int u = q.front(); q.pop();
        int nextLevel = dist[u] + 1;

        const long long b = off[u], e = off[u + 1];
        for (long long p = b; p < e; ++p) {
            int v = edges[p];
            if ((unsigned)v >= (unsigned)n) continue;
            if (dist[v] != INT_MAX/4) continue;
            dist[v] = nextLevel;
            q.push(v);
            ++reached;
        }
        if (e > b && nextLevel > maxLevel) maxLevel = nextLevel;
    }

    auto t1 = clock::now();
    std::cout << "BFS(raw) time: "
              << std::chrono::duration<double>(t1 - t0).count() << " s\n";
    return dist;
}
