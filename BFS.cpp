// ==========================================================
// File: BFS_block_memprofile.cpp (RSS-only profiling)
// - Encoded BFS (block-indexed, on-demand decode) with OS-level
//   Resident Set Size (RSS) measurements in BYTES.
//   We report: peak block RSS Δ observed while loading a block.
// - Raw BFS prints RSS at entry/exit and Δ.
// - No allocator/capacity math; no external profilers.
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
#include "decode.h"          // extract_block_slices(...), readBits(...)
#include "huffman_tree.h"    // struct HuffmanNode { HuffmanNode *l,*r; int v; }
#include "misc.h"
#include <malloc.h>
using namespace std;

static inline size_t rss_bytes_strict() {
#if defined(__linux__)
    if (FILE* f = fopen("/proc/self/smaps_rollup", "r")) {
        char key[64]; unsigned long kb;
        while (fscanf(f, "%63s %lu", key, &kb) == 2) {
            if (strcmp(key, "Rss:") == 0) { fclose(f); return (size_t)kb * 1024ull; }
        }
        fclose(f);
    }
    if (FILE* f = fopen("/proc/self/status", "r")) {
        char key[64]; unsigned long kb;
        while (fscanf(f, "%63s %lu", key, &kb) == 2) {
            if (strncmp(key, "VmRSS:", 6) == 0) { fclose(f); return (size_t)kb * 1024ull; }
        }
        fclose(f);
    }
    if (FILE* f = fopen("/proc/self/statm", "r")) {
        unsigned long pages=0, dummy=0;
        if (fscanf(f, "%lu %lu", &dummy, &pages) == 2) {
            fclose(f);
            long ps = sysconf(_SC_PAGESIZE); if (ps <= 0) ps = 4096;
            return (size_t)pages * (size_t)ps;
        }
        fclose(f);
    }
#endif
    // portable fallback
    rusage ru; getrusage(RUSAGE_SELF, &ru);
#if defined(__APPLE__) && defined(__MACH__)
    return (size_t)(ru.ru_maxrss / 1024) * 1024ull;
#else
    return (size_t)ru.ru_maxrss * 1024ull;
#endif
}

// =============================================================
// Single-source BFS (encoded graph, block-indexed, parallel),
// with RSS-only accounting of block loads.
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

    // RSS-based accounting for encoded BFS (OS-level resident bytes)
    size_t peak_block_rss_bytes = 0; // max RSS delta seen on a single block load

    while (!frontier.empty()) {
        const int nextLevel = maxLevel + 1;
        nextFrontier.clear();

        // Parallel over the current frontier
        #pragma omp parallel
        {
            // Thread-local block cache (avoid contention)
            int t_curBlk = -1, t_blk_i0 = 0, t_blk_i1 = 0;
            std::vector<int>       t_deg_blk;
            std::vector<uint16_t>  t_bit_blk, t_huff_blk;
            std::vector<uint8_t>   t_hiSlice, t_loSlice;
            uint32_t               t_baseHiOffset = 0, t_baseLoOffset = 0;
            std::vector<uint64_t>  t_offHiBits, t_offLoBits; // helper arrays (not measured specifically)

            std::vector<int> t_localNext;  t_localNext.reserve(1024); // thread-local buffer

            auto ensure_block_loaded = [&](int blk) {
                if (blk == t_curBlk) return;
                t_curBlk = blk;
                t_blk_i0 = blk * BLOCK;
                t_blk_i1 = std::min(n, t_blk_i0 + BLOCK);

                // Release previous block buffers to reduce allocator noise
                std::vector<uint8_t>().swap(t_hiSlice);
                std::vector<uint8_t>().swap(t_loSlice);
                std::vector<int>().swap(t_deg_blk);
                std::vector<uint16_t>().swap(t_bit_blk);
                std::vector<uint16_t>().swap(t_huff_blk);
                malloc_trim(0);  // try to return free memory to OS
                t_baseHiOffset = t_baseLoOffset = 0;

                // Measure RSS delta around the block load.
                size_t r0 = 0, r1 = 0;
                #pragma omp critical(__bfs_block_rss)
                {
                    r0 = rss_bytes_strict();

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

                    r1 = rss_bytes_strict();
                }
                const size_t loaded_bytes = (t_hiSlice.size() + t_loSlice.size());
                const size_t my_block_rss = (loaded_bytes == 0) ? 0 : (r1 > r0 ? (r1 - r0) : 0);
                #pragma omp critical(__bfs_peak_block_rss)
                peak_block_rss_bytes = std::max(peak_block_rss_bytes, my_block_rss);


                // Precompute per-entry bit offsets
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

            #pragma omp for schedule(dynamic, 1024)
            for (int i = 0; i < (int)frontier.size(); ++i) {
                int u = frontier[i];
                if (deg[u] == 0) continue; // isolated

                const int blk = u / BLOCK;
                ensure_block_loaded(blk);

                const int k = u - t_blk_i0;
                const uint32_t hiOff = t_baseHiOffset + (uint32_t)t_offHiBits[k];
                const uint32_t loOff = t_baseLoOffset + (uint32_t)t_offLoBits[k];

                // Decode neighbors of u from the thread-local cached slices
                std::vector<int> nbrs = decode_single_entry_from_slices(
                    &t_deg_blk[k], &t_bit_blk[k], &t_huff_blk[k],
                    t_hiSlice.data(), t_loSlice.data(),
                    hiOff, loOff,
                    treeOpp, fallbackBitsOpp
                );

                // Visit neighbors
                for (int v : nbrs) {
                    if ((unsigned)v >= (unsigned)n) continue;
                    if (dist[v] != INF) continue;
                    // mark visited once
                    #pragma omp critical(__bfs_visit_dist)
                    {
                        if (dist[v] == INF) {
                            dist[v] = nextLevel;
                            t_localNext.push_back(v);
                        }
                    }
                }
            } // for frontier

            // Merge thread-local next frontier
            #pragma omp critical(__bfs_merge_next)
            {
                reached += (int)t_localNext.size();
                nextFrontier.insert(nextFrontier.end(),
                                    t_localNext.begin(), t_localNext.end());
            }
        } // parallel

        // Deduplicate nextFrontier
        std::sort(nextFrontier.begin(), nextFrontier.end());
        nextFrontier.erase(std::unique(nextFrontier.begin(), nextFrontier.end()), nextFrontier.end());

        frontier.swap(nextFrontier);
        if (!frontier.empty()) maxLevel = nextLevel;
    }

    auto t1 = clock::now();

    // Report RSS-based block memory usage (bytes)
    std::cout << std::fixed;
    std::cout << "BFS(Hybhuff): peak block RSS \xCE\x94 = " << peak_block_rss_bytes << " bytes"
              << "Final RSS: " << rss_bytes()
                << " | time = " << std::chrono::duration<double>(t1 - t0).count() << " s\n";
    
}

// =============================================================
// BFS without decoding: uses deg[] + edges[] directly.
// Prints RSS at entry/exit and delta in BYTES.
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

    // Build CSR offsets once (prefix sum)
    std::vector<long long> off(n + 1, 0);
    for (int i = 0; i < n; ++i) off[i + 1] = off[i] + (long long)deg[i];

    size_t r_before = rss_bytes();
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
    size_t r_after = rss_bytes();
    
    std::cout << "BFS(raw): reached=" << reached
              << " | RSS \xCE\x94 = " << (r_after > r_before ? (r_after - r_before) : 0) << " bytes"
              << " | time=" << std::chrono::duration<double>(t1 - t0).count() << " s\n";
    return dist;
}
