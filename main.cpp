// ============================================================================
// File: main.cpp (complete, buildable)
// Driver for HybHuff encode → estimate footprint → run BFS / K-core / PageRank
// on both compressed (on-demand decode) and raw baselines, plus a full decode.
// ============================================================================

#include <chrono>
#include <iostream>
#include <vector>
#include <queue>
#include <unordered_map>
#include <unordered_set>
#include <numeric>
#include <algorithm>
#include <sys/resource.h>
#include <unistd.h>
#include <cstdlib>
#include <cstdint>     // uint8_t, uint16_t, uint64_t

#include "hypergraph.h"
#include "huffman_tree.h"
#include "huffman_code.h"
#include "encode.h"
#include "decode.h"
#include "BFS.h"
#include "computeKCore.h"
#include "pagerank.h"
#include "misc.h"
using namespace std;

int main(int argc, char** argv) {
    if (argc != 3) {
        cerr << "Usage: " << argv[0] << " <pct> <input.hyper>\n";
        return 1;
    }

    // -------------------------------------------------------------------------
    // Parse params and load hypergraph
    //   pct   : fraction (0–1) of items to encode via Huffman (passed as <pct>%)
    //   fname : input hypergraph file
    //   encodeH: choose which side to encode (encode hyperedges if nh >= nv)
    // -------------------------------------------------------------------------
    double pct = atof(argv[1]) * 0.01;
    const char* fname = argv[2];

    Hypergraph* G = readHypergraph(fname);
    int nv = G->nv, nh = G->nh;
    bool encodeH = nv <= nh;

    // -------------------------------------------------------------------------
    // Encoding: build Huffman tree on the opposite side and encode the chosen side
    // Outputs:
    //   huffCount[i]  : #Huffman symbols in entry i
    //   bitCount[i]   : #bits occupied by the Huffman section of entry i
    //   hi / lo       : concatenated bitstreams (Huffman / fallback)
    //   blockHi/Lo    : start bit offsets for each block (BLOCK entries)
    //   bitsH / bitsL : total bits across the encoded side in hi/lo
    // -------------------------------------------------------------------------
    cout << "Encoding\n";

    int fbV = 0, fbH = 0;                    // fallback bit-widths (V, H)
    HuffmanNode *tV = nullptr, *tH = nullptr;
    unordered_map<int,string> cV, cH;        // codes (not needed at runtime)
    uint16_t *bitCount = nullptr;
    uint16_t *huffCount = nullptr;
    uint8_t  *hi = nullptr, *lo = nullptr;
    int bitsH = 0, bitsL = 0;

    uint64_t *blockHi = nullptr, *blockLo = nullptr;
    int nBlocks = 0;
    const int BLOCK = 1024;                   // block size in entries (tunable)

    auto t_encode0 = chrono::high_resolution_clock::now();

    if (encodeH) {
        tV = buildTree(G->degreeV, nv, pct, fbV);
        genCodes(tV, cV);
        tie(bitsH, bitsL) = encodeSide(
            nh, G->degreeH, G->edgesH, cV, fbV,
            huffCount, bitCount, hi, lo,
            BLOCK, blockHi, blockLo, &nBlocks
        );
    } else {
        tH = buildTree(G->degreeH, nh, pct, fbH);
        genCodes(tH, cH);
        tie(bitsH, bitsL) = encodeSide(
            nv, G->degreeV, G->edgesV, cH, fbH,
            huffCount, bitCount, hi, lo,
            BLOCK, blockHi, blockLo, &nBlocks
        );
    }

    auto t_encode1 = chrono::high_resolution_clock::now();
    cout << "Encoding Time: " << chrono::duration<double>(t_encode1 - t_encode0).count() << " s\n";

    // -------------------------------------------------------------------------
    // Footprint estimate: content (streams + metadata) + Huffman tree
    //   Choose integer widths (uint16_t vs int) based on max deg/bitcount.
    // -------------------------------------------------------------------------
    size_t treeSize = encodeH ? treeMem(tV) : treeMem(tH);
    int n = encodeH ? nh : nv;

    int maxDeg = 0, maxBC = 0;
    for (int i = 0; i < n; ++i) {
        maxDeg = max(maxDeg, encodeH ? G->degreeH[i] : G->degreeV[i]);
        maxBC  = max(maxBC, (int)bitCount[i]);
    }

    size_t contentSize = 0;
    if (maxDeg < (1 << 16) && maxBC < (1 << 16)) {
        contentSize = footprint_both(n, bitsH, bitsL);
    } else if (maxDeg < (1 << 16) || maxBC < (1 << 16)) {
        contentSize = footprint_one(n, bitsH, bitsL);
    } else {
        contentSize = footprint(n, bitsH, bitsL);
    }

    cout << "Content size: " << (contentSize / 1024) << " KB\n";
    cout << "Tree    size: " << (treeSize   / 1024) << " KB\n";
    cout << "Total   size: " << ((contentSize + treeSize) / 1024) << " KB\n";

    // -------------------------------------------------------------------------
    // BFS: compressed (on-demand decode) vs raw
    //   Compressed: neighbors decoded per-vertex from block slices.
    //   Raw: standard CSR BFS (deg[] + edges[]).
    // -------------------------------------------------------------------------
    auto t_bfs0 = chrono::high_resolution_clock::now();
    for (int i = 0; i < 5; i++) {
        runBFSFromSingleSource(
            encodeH ? nh : nv,
            encodeH ? G->degreeH : G->degreeV,
            encodeH ? G->edgesH  : G->edgesV,
            huffCount, bitCount, hi, lo,
            encodeH ? tV : tH, encodeH ? fbV : fbH,
            i,
            blockHi, blockLo, BLOCK
        );
    }
    auto t_bfs1 = chrono::high_resolution_clock::now();
    cout << "BFS(compressed) Time: "
         << chrono::duration<double>(t_bfs1 - t_bfs0).count() << " s\n";

    auto t_bfs2 = chrono::high_resolution_clock::now();
    for (int i = 0; i < 5; i++) {
        auto dist_raw = bfs_single_source_raw(
            encodeH ? nh : nv,
            encodeH ? G->degreeH : G->degreeV,
            encodeH ? G->edgesH  : G->edgesV,
            i
        );
        (void)dist_raw; // used only for timing in this driver
    }
    auto t_bfs3 = chrono::high_resolution_clock::now();
    cout << "BFS(raw) Time: "
         << chrono::duration<double>(t_bfs3 - t_bfs2).count() << " s\n";

    // -------------------------------------------------------------------------
    // K-core: compressed (decode-on-demand) vs raw
    //   We sweep k=1..5 on compressed, and k=1..2 on raw for timing.
    // -------------------------------------------------------------------------
    auto t_kcore0 = chrono::high_resolution_clock::now();
    for (int i = 1; i <= 5; i++) {
        auto inCore_dec = computeKCore_onDemand_decodeRandom(
            n,
            encodeH ? G->degreeH : G->degreeV,
            huffCount, bitCount,
            hi, lo,
            encodeH ? tV : tH, encodeH ? fbV : fbH,
            blockHi, blockLo, BLOCK,
            i
        );
        (void)inCore_dec;
    }
    auto t_kcore1 = chrono::high_resolution_clock::now();
    cout << "K-core(compressed) Time: "
         << chrono::duration<double>(t_kcore1 - t_kcore0).count() << " s\n";

    auto t_kcore2 = chrono::high_resolution_clock::now();
    for (int i = 1; i <= 2; i++) {
        auto inCore_raw = computeKCore_raw(
            n,
            encodeH ? G->degreeH : G->degreeV,
            encodeH ? G->edgesH  : G->edgesV,
            i
        );
        (void)inCore_raw;
    }
    auto t_kcore3 = chrono::high_resolution_clock::now();
    cout << "K-core(raw) Time: "
         << chrono::duration<double>(t_kcore3 - t_kcore2).count() << " s\n";

    // -------------------------------------------------------------------------
    // PageRank: compressed (decode-on-demand) and raw
    //   pr_iters is provided by pagerank.h / build system; alpha = 0.85.
    //   Pass edges[] as needed by your verification settings.
    // -------------------------------------------------------------------------
    runPageRankOnDemand_decodeRandom(
        encodeH ? nh : nv,
        encodeH ? G->degreeH : G->degreeV,
        encodeH ? G->edgesH  : G->edgesV,    // pass nullptr if per-vertex verification is disabled
        huffCount, bitCount,
        hi, lo,
        encodeH ? tV : tH, encodeH ? fbV : fbH,
        20, 0.85,
        blockHi, blockLo, BLOCK
    );

    runPageRankRaw(
        encodeH ? nh : nv,
        encodeH ? G->degreeH : G->degreeV,
        encodeH ? G->edgesH  : G->edgesV,
        20,
        0.85
    );

    // -------------------------------------------------------------------------
    // Full paged decode of the encoded side (optional verification if edges!=nullptr)
    // -------------------------------------------------------------------------
    using clock = std::chrono::high_resolution_clock;
    auto t0 = clock::now();
    auto nbrs = decodeSide(
        encodeH ? nh : nv,
        encodeH ? G->degreeH : G->degreeV,
        encodeH ? G->edgesH  : G->edgesV,
        huffCount, bitCount,
        hi, lo,
        encodeH ? tV : tH,
        encodeH ? fbV : fbH,
        blockHi, blockLo, BLOCK
    );
    (void)nbrs; // not used further in this driver
    auto t1 = clock::now();
    cout << "Decode Time: " << std::chrono::duration<double>(t1 - t0).count() << " s\n";

    // -------------------------------------------------------------------------
    // Cleanup
    // -------------------------------------------------------------------------
    freeTree(tV);
    freeTree(tH);
    delete[] bitCount;
    delete[] huffCount;
    delete[] hi;
    delete[] lo;
    delete[] blockHi;
    delete[] blockLo;
    delete G;

    return 0;
}
