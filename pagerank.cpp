// ============================================================================
// pagerank_rolling_singlefile.cpp
// Purpose:
//   • PageRank on (1) raw CSR and (2) encoded graph via on-demand block decoding.
//   • Encoded path uses per-block slices and rolling bit pointers (no O(n) scans).
// Notes:
//   • Keeps original names/signatures.
//   • Parallelized with OpenMP where annotated.
// ============================================================================

#include <vector>
#include <cstdint>
#include <atomic>
#include <chrono>
#include <iostream>
#include <algorithm>
#include <omp.h>
#include <sys/resource.h>
#include <unistd.h>
#include <cstdio>
#include "misc.h"
#include "huffman_tree.h"
#include "decode.h"

using namespace std;

// ----------------------------------------------------------------------------
// runPageRankRaw
// ----------------------------------------------------------------------------
// Raw PageRank over flat CSR (deg[] + edges[]).
// Strategy:
//   1) Build prefix offsets off[] from deg[] for O(1) neighbor range.
//   2) Push-based updates: pr[u]/deg[u] added to each neighbor v.
//   3) Handle dangling mass and damping each iteration.
// Inputs:
//   n      : number of vertices
//   deg    : out-degree per vertex
//   edges  : flat neighbor list of length sum(deg)
//   iters  : number of iterations
//   alpha  : damping factor (typically 0.85)
void runPageRankRaw(
    int n,
    const int* deg,
    const int* edges,
    int iters,
    double alpha /*=0.85*/
) {
    using clock = std::chrono::high_resolution_clock;
    auto t0 = clock::now();

    if (n <= 0) { std::cout << "PageRank(raw): empty graph\n"; return; }

    // Prefix offsets into edges[]
    std::vector<uint64_t> off(n + 1, 0);
    for (int i = 0; i < n; ++i) off[i + 1] = off[i] + (uint64_t)deg[i];

    std::vector<double> pr(n, 1.0 / n);
    std::vector<double> next(n, 0.0);

    for (int it = 0; it < iters; ++it) {
        std::fill(next.begin(), next.end(), 0.0);
        double dangling = 0.0;

        // Push contributions
        for (int u = 0; u < n; ++u) {
            const int du = deg[u];
            if (du == 0) { dangling += pr[u]; continue; }

            const double share = pr[u] / du;
            const uint64_t s = off[u];
            const uint64_t e = off[u + 1];

            for (uint64_t p = s; p < e; ++p) {
                int v = edges[p];
                if ((unsigned)v < (unsigned)n) next[v] += share;
            }
        }

        // Teleport + distribute dangling mass
        const double base = (1.0 - alpha) / n;
        const double dang = dangling / n;
        for (int v = 0; v < n; ++v)
            pr[v] = base + alpha * (next[v] + dang);
    }

    auto t1 = clock::now();
    std::cout << "PageRank (raw edges): "
              << std::chrono::duration<double>(t1 - t0).count() << " s\n";
}

// ----------------------------------------------------------------------------
// runPageRankOnDemand_decodeRandom
// ----------------------------------------------------------------------------
// PageRank decoding neighbors on-demand from compressed streams.
// Strategy:
//   • Iterate blocks [i0, i1) of size BLOCK.
//   • For each block, extract minimal slices (hi/lo) once.
//   • Maintain intra-block bit offsets to decode each vertex without re-scanning.
//   • Push contributions into 'next' with atomics; handle dangling + damping.
// Inputs mirror the encoded representation used elsewhere.
void runPageRankOnDemand_decodeRandom(
    int n, const int* deg, const int* edges,
    const uint16_t* huffCount, const uint16_t* bitCount,
    const uint8_t* hi, const uint8_t* lo,
    HuffmanNode* treeOpp, int fallbackBitsOpp,
    int iters, double alpha,
    const uint64_t* blockHi, const uint64_t* blockLo, int BLOCK
){
    using clock = std::chrono::high_resolution_clock;
    auto t0 = clock::now();

    if (n <= 0) { std::cout << "PageRank: empty graph\n"; return; }
    if (!deg || !huffCount || !bitCount || !hi || !lo || !treeOpp ||
        !blockHi || !blockLo || BLOCK <= 0) {
        std::cerr << "PageRank(decode): missing/invalid inputs\n";
        return;
    }

    const int nBlocks = (n + BLOCK - 1) / BLOCK;

    std::vector<double> pr(n, 1.0 / n);
    std::vector<double> next(n, 0.0);

    // (Internally tracked metric; not discussed here)
    size_t peak_block_rss_bytes = 0;

    for (int it = 0; it < iters; ++it) {
        std::fill(next.begin(), next.end(), 0.0);
        double dangling = 0.0;

        // Parallel over blocks; each thread caches slices for its block
        #pragma omp parallel
        {
            std::vector<int>       t_deg_blk;
            std::vector<uint16_t>  t_bit_blk, t_huff_blk;
            std::vector<uint8_t>   t_hiSlice, t_loSlice;
            uint32_t               t_baseHiOffset = 0, t_baseLoOffset = 0;
            std::vector<uint64_t>  t_offHiBits, t_offLoBits;
            std::vector<int>       nbrs; nbrs.reserve(128);

            double local_dangling = 0.0;

            #pragma omp for schedule(dynamic,1) reduction(+:dangling)
            for (int b = 0; b < nBlocks; ++b) {
                const int i0 = b * BLOCK;
                const int i1 = std::min(n, i0 + BLOCK);
                if (i0 >= i1) continue;

                // Extract minimal slices once for this block
                t_deg_blk.clear(); t_bit_blk.clear(); t_huff_blk.clear();
                t_hiSlice.clear();  t_loSlice.clear();
                t_baseHiOffset = t_baseLoOffset = 0;
                std::vector<uint8_t>().swap(t_hiSlice);
                std::vector<uint8_t>().swap(t_loSlice);
                std::vector<int>().swap(t_deg_blk);
                std::vector<uint16_t>().swap(t_bit_blk);
                std::vector<uint16_t>().swap(t_huff_blk);

                size_t r0 = 0, r1 = 0;
                #pragma omp critical(__pr_block_rss)
                {
                    r0 = rss_bytes();
                    extract_block_slices(
                        i0, i1, BLOCK,
                        deg, bitCount, huffCount,
                        hi, lo,
                        blockHi, blockLo,
                        fallbackBitsOpp,
                        t_deg_blk, t_bit_blk, t_huff_blk,
                        t_hiSlice, t_loSlice,
                        t_baseHiOffset, t_baseLoOffset
                    );
                    r1 = rss_bytes();
                }
                size_t my_block_rss = (r1 > r0 ? (r1 - r0) : 0);
                #pragma omp critical(__pr_peak_block_rss)
                peak_block_rss_bytes = std::max(peak_block_rss_bytes, my_block_rss);

                // Precompute intra-block bit offsets for each entry
                const int L = i1 - i0;
                t_offHiBits.assign(L, 0);
                t_offLoBits.assign(L, 0);
                uint64_t accH = 0, accL = 0;
                for (int k = 0; k < L; ++k) {
                    t_offHiBits[k] = accH;
                    t_offLoBits[k] = accL;
                    accH += (uint64_t)t_bit_blk[k];
                    accL += (uint64_t)(t_deg_blk[k] - t_huff_blk[k]) * (uint64_t)fallbackBitsOpp;
                }

                // Process vertices in this block
                for (int u = i0; u < i1; ++u) {
                    const int du = deg[u];
                    if (du == 0) { local_dangling += pr[u]; continue; }

                    const int k = u - i0;
                    const uint32_t hiOff = t_baseHiOffset + (uint32_t)t_offHiBits[k];
                    const uint32_t loOff = t_baseLoOffset + (uint32_t)t_offLoBits[k];

                    // Decode u’s neighbors from cached slices
                    nbrs = decode_single_entry_from_slices(
                        &t_deg_blk[k], &t_bit_blk[k], &t_huff_blk[k],
                        t_hiSlice.data(), t_loSlice.data(),
                        hiOff, loOff,
                        treeOpp, fallbackBitsOpp
                    );

                    const double share = pr[u] / du;
                    for (int v : nbrs) {
                        if ((unsigned)v < (unsigned)n) {
                            #pragma omp atomic
                            next[v] += share;
                        }
                    }
                }

                // Add per-thread dangling to the reduction
                dangling += local_dangling;
                local_dangling = 0.0;
            }
        }

        // Teleport + distribute dangling mass
        const double base = (1.0 - alpha) / n;
        const double dang = dangling / n;

        #pragma omp parallel for schedule(static)
        for (int v = 0; v < n; ++v) {
            pr[v] = base + alpha * (next[v] + dang);
        }
    }

    auto t1 = clock::now();
    std::cout << "PageRank (Hybhuff) time: "
              << std::chrono::duration<double>(t1 - t0).count() << " s\n";
}
