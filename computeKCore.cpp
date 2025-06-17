#include "computeKCore.h"
#include <queue>
#include <iostream>
#include <chrono>
#include <sys/resource.h>
#include <unistd.h>

using namespace std;

/**
 * @brief Compute the k-core of a bipartite hypergraph.
 *
 * The algorithm performs iterative pruning of vertices whose degree
 * (number of incident hyperedges) is less than `k_thresh`. Once such
 * a vertex is removed, any hyperedge it participates in is marked
 * removed, and its incident vertices are re-evaluated.
 *
 * This process continues until all remaining vertices meet the
 * minimum degree requirement.
 *
 * @param v2he     Vertex-to-hyperedge adjacency list.
 * @param he2v     Hyperedge-to-vertex adjacency list.
 * @param k_thresh Minimum degree threshold (k).
 * @return         Binary vector of size nv. `1` = in k-core, `0` = removed.
 */
vector<char> computeKCore(const vector<vector<int>>& v2he,
                          const vector<vector<int>>& he2v,
                          int k_thresh) {
    int nv = (int)v2he.size();  // number of vertices
    int nh = (int)he2v.size();  // number of hyperedges

    vector<int> degV(nv);           // degree of each vertex
    vector<int> removedHE(nh, 0);   // whether each hyperedge is removed
    vector<char> removedV(nv, 0);   // whether each vertex is removed
    queue<int> Q;                   // vertices queued for removal

    // Step 1: Initialize vertex degrees and mark those below threshold
    for (int v = 0; v < nv; ++v) {
        degV[v] = (int)v2he[v].size();
        if (degV[v] < k_thresh) {
            removedV[v] = 1;
            Q.push(v);
        }
    }

    // Step 2: Iteratively prune vertices and affected neighbors
    while (!Q.empty()) {
        int v = Q.front(); Q.pop();
        for (int e : v2he[v]) {
            // Only remove each hyperedge once
            if (++removedHE[e] == 1) {
                for (int u : he2v[e]) {
                    if (!removedV[u] && --degV[u] < k_thresh) {
                        removedV[u] = 1;
                        Q.push(u);
                    }
                }
            }
        }
    }

    // Step 3: Mark remaining vertices as in-core
    vector<char> inCore(nv);
    int remain = 0;
    for (int v = 0; v < nv; ++v) {
        inCore[v] = !removedV[v];
        if (inCore[v]) ++remain;
    }

    cout << "K-Core: " << remain << " vertices remain in core\n";
    return inCore;
}
