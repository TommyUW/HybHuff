// ============================================================================
// File: main.cpp
//
// Description:
// This program demonstrates a memory-efficient hypergraph compression algorithm
// using a hybrid Huffman + bitwise fallback encoding scheme. The hypergraph is 
// encoded by building a Huffman tree on the more frequent side (vertices or 
// hyperedges), and encoding the other side using a mix of Huffman and fixed-width
// binary codes.
//
// The encoded structure is then validated with decoding and used to perform 
// analyses such as BFS, connected components, k-core decomposition, label 
// propagation, and PageRank.
//
// The memory footprint is estimated based on whether the per-element bitCount
// and degree arrays can be stored as uint16_t.
//
// Developed by:
//   Tianyu Zhao, University of Washington
//
// Research Advisors:
//   Prof. Dongfang Zhao  
//   Dr. Luanzheng Guo  
//   Dr. Nathan Tallent
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

#include "hypergraph.h"
#include "huffman_tree.h"
#include "huffman_code.h"
#include "encode.h"
#include "decode.h"
#include "BFS.h"
#include "computeKCore.h"
#include "label_propagation.h"
#include "pagerank.h"

using namespace std;

// -----------------------------------------------------------------------------
// Memory footprint estimation functions
// -----------------------------------------------------------------------------

static size_t footprint(int n, long bitsH, long bitsL) {
    long bH = (bitsH + 7) / 8, bL = (bitsL + 7) / 8;
    return 2 * sizeof(int*)
         + 2 * sizeof(uint8_t*)
         + sizeof(int) * n
         + sizeof(int) * n
         + bH * sizeof(uint8_t)
         + bL * sizeof(uint8_t);
}

static size_t footprint_one(int n, long bitsH, long bitsL) {
    long bH = (bitsH + 7) / 8, bL = (bitsL + 7) / 8;
    return 2 * sizeof(int*)
         + 2 * sizeof(uint8_t*)
         + sizeof(uint16_t) * n
         + sizeof(int) * n
         + bH * sizeof(uint8_t)
         + bL * sizeof(uint8_t);
}

static size_t footprint_both(int n, long bitsH, long bitsL) {
    long bH = (bitsH + 7) / 8, bL = (bitsL + 7) / 8;
    return 2 * sizeof(int*)
         + 2 * sizeof(uint8_t*)
         + sizeof(uint16_t) * n
         + sizeof(uint16_t) * n
         + bH * sizeof(uint8_t)
         + bL * sizeof(uint8_t);
}

int main(int argc, char** argv) {
    if (argc != 6) {
        cerr << "Usage: " << argv[0]
             << " <pct> <input.hyper> <k_thresh> <lp_iters> <pr_iters>\n";
        return 1;
    }

    // -------------------------------------------------------------------------
    // Parse input parameters
    // -------------------------------------------------------------------------
    double pct = atof(argv[1]) * 0.01;
    const char* fname = argv[2];
    int k_thresh  = atoi(argv[3]);
    int lp_iters  = atoi(argv[4]);
    int pr_iters  = atoi(argv[5]);

    Hypergraph* G = readHypergraph(fname);
    int nv = G->nv, nh = G->nh;
    bool encodeH = nv <= nh;

    cout << "Compressing " << (encodeH ? "hyperedges" : "vertices")
         << " (tree on " << (encodeH ? "vertices" : "hyperedges") << ")\n";

    // -------------------------------------------------------------------------
    // Encoding
    // -------------------------------------------------------------------------
    int fbV=0, fbH=0;
    HuffmanNode *tV=nullptr, *tH=nullptr;
    unordered_map<int,string> cV, cH;
    int *bitCount = nullptr;
    uint8_t *hi=nullptr, *lo=nullptr;
    int bitsH = 0, bitsL = 0;

    auto t_encode0 = chrono::high_resolution_clock::now();

    if (encodeH) {
        tV = buildTree(G->degreeV, nv, pct, fbV);
        genCodes(tV, cV);
        tie(bitsH,bitsL) = encodeSide(nh, G->degreeH, G->edgesH, cV, fbV, bitCount, hi, lo);
    } else {
        tH = buildTree(G->degreeH, nh, pct, fbH);
        genCodes(tH, cH);
        tie(bitsH,bitsL) = encodeSide(nv, G->degreeV, G->edgesV, cH, fbH, bitCount, hi, lo);
    }

    auto t_encode1 = chrono::high_resolution_clock::now();

    // -------------------------------------------------------------------------
    // Decoding
    // -------------------------------------------------------------------------
    auto t_decode0 = chrono::high_resolution_clock::now();
    auto adj = encodeH
      ? decodeSide(nh, G->degreeH, G->edgesH, bitCount, hi, lo, tV, fbV)
      : decodeSide(nv, G->degreeV, G->edgesV, bitCount, hi, lo, tH, fbH);
    auto t_decode1 = chrono::high_resolution_clock::now();

    // -------------------------------------------------------------------------
    // Reconstruct bipartite structure
    // -------------------------------------------------------------------------
    vector<vector<int>> v2he, he2v;
    if (encodeH) {
        he2v = adj;
        v2he.assign(nv, {});
        for (int e = 0; e < nh; ++e)
            for (int v : he2v[e])
                v2he[v].push_back(e);
    } else {
        v2he = adj;
        he2v.assign(nh, {});
        for (int v = 0; v < nv; ++v)
            for (int e : v2he[v])
                he2v[e].push_back(v);
    }

    // -------------------------------------------------------------------------
    // Estimate memory footprint
    // -------------------------------------------------------------------------
    size_t treeSize = encodeH ? treeMem(tV) : treeMem(tH);
    int n = encodeH ? nh : nv;

    int maxDeg = 0, maxBC = 0;
    for (int i = 0; i < n; ++i) {
        maxDeg = max(maxDeg, encodeH ? G->degreeH[i] : G->degreeV[i]);
        maxBC  = max(maxBC, bitCount[i]);
    }

    size_t contentSize = 0;
    if (maxDeg < (1 << 16) && maxBC < (1 << 16)) {
        contentSize = footprint_both(n, bitsH, bitsL);
        cout << "[Footprint] Using uint16_t for both degree and bitCount\n";
    } else if (maxDeg < (1 << 16) || maxBC < (1 << 16)) {
        contentSize = footprint_one(n, bitsH, bitsL);
        cout << "[Footprint] Using uint16_t for one of degree or bitCount\n";
    } else {
        contentSize = footprint(n, bitsH, bitsL);
        cout << "[Footprint] Using full int32_t for both degree and bitCount\n";
    }

    cout << "Content size: " << (contentSize / 1024) << " KB\n";
    cout << "Tree    size: " << (treeSize / 1024) << " KB\n";
    cout << "Total   size: " << ((contentSize + treeSize) / 1024) << " KB\n";

    // -------------------------------------------------------------------------
    // Analytic workloads: BFS, k-core, label propagation, PageRank
    // -------------------------------------------------------------------------
    auto t_bfs0 = chrono::high_resolution_clock::now();

    vector<int> offsetV(nv + 1), offsetH(nh + 1);
    runConnectedComponents(G, nv, nh, adj, offsetV, offsetH);

    runBFS(nv, nh, G->degreeV, G->edgesV, G->degreeH, G->edgesH,
           offsetV.data(), offsetH.data());

    auto t_bfs1 = chrono::high_resolution_clock::now();

    //runFullBFSAllVertices(adj);

    auto t_bfs2 = chrono::high_resolution_clock::now();

    auto t_k0 = chrono::high_resolution_clock::now();

    auto inCore = computeKCore(v2he, he2v, k_thresh);
    runLabelPropagation(v2he, he2v, inCore, lp_iters);

    auto t_k1 = chrono::high_resolution_clock::now();

    auto t_p0 = chrono::high_resolution_clock::now();
    runPageRank(v2he, he2v, pr_iters);
    auto t_p1 = chrono::high_resolution_clock::now();

    // -------------------------------------------------------------------------
    // Timing summary
    // -------------------------------------------------------------------------
    cout << "Encoding Time: " << chrono::duration<double>(t_encode1 - t_encode0).count() << " s\n";
    cout << "Decoding Time: " << chrono::duration<double>(t_decode1 - t_decode0).count() << " s\n";
    cout << "Total (Enc+Dec): " << chrono::duration<double>(t_decode1 - t_encode0).count() << " s\n";

    cout << "BFS Time: " << chrono::duration<double>(t_bfs1 - t_bfs0).count() << " s\n"; 
    cout << "BFS All Vertices Time: " << chrono::duration<double>(t_bfs2 - t_bfs0).count() << " s\n"; 
    cout << "K-Core + Label Propagation Time: " << chrono::duration<double>(t_k1 - t_k0).count() << " s\n"; 
    cout << "PageRank Time: " << chrono::duration<double>(t_p1 - t_p0).count() << " s\n"; 

    // -------------------------------------------------------------------------
    // Cleanup
    // -------------------------------------------------------------------------
    freeTree(tV);
    freeTree(tH);
    delete[] bitCount;
    delete[] hi;
    delete[] lo;
    delete G;

    return 0;
}
