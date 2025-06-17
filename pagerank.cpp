#include "pagerank.h"
#include <vector>
#include <iostream>
#include <chrono>
#include <sys/resource.h>
#include <unistd.h>

using namespace std;

/**
 * @brief Run PageRank on a bipartite hypergraph.
 *
 * This function computes PageRank scores over a bipartite hypergraph using
 * alternating updates between vertices and hyperedges.
 * 
 * The algorithm performs the following for a fixed number of iterations:
 *   1. Distribute scores from vertices to hyperedges (V → E).
 *   2. Distribute scores from hyperedges to vertices (E → V).
 *   3. Apply damping factor to simulate random restarts.
 *
 * @param v2he  Vertex-to-hyperedge adjacency list.
 * @param he2v  Hyperedge-to-vertex adjacency list.
 * @param iters Number of iterations to run.
 */
void runPageRank(const vector<vector<int>>& v2he,
                 const vector<vector<int>>& he2v,
                 int iters) {
    int nv = (int)v2he.size();
    int nh = (int)he2v.size();

    vector<double> prV(nv, 1.0 / nv);  // PageRank scores for vertices
    vector<double> prE(nh, 1.0 / nh);  // PageRank scores for hyperedges
    vector<double> nextV(nv), nextE(nh);  // Buffers for next iteration
    const double alpha = 0.85;  // Damping factor

    for (int it = 0; it < iters; ++it) {
        fill(nextV.begin(), nextV.end(), 0.0);
        fill(nextE.begin(), nextE.end(), 0.0);

        // Step 1: V → E
        for (int v = 0; v < nv; ++v) {
            if (!v2he[v].empty()) {
                double share = prV[v] / v2he[v].size();
                for (int e : v2he[v]) nextE[e] += share;
            }
        }

        // Step 2: E → V
        for (int e = 0; e < nh; ++e) {
            if (!he2v[e].empty()) {
                double share = prE[e] / he2v[e].size();
                for (int v : he2v[e]) nextV[v] += share;
            }
        }

        // Step 3: Apply damping
        for (int v = 0; v < nv; ++v)
            prV[v] = (1.0 - alpha) / nv + alpha * nextV[v];
        for (int e = 0; e < nh; ++e)
            prE[e] = (1.0 - alpha) / nh + alpha * nextE[e];
    }

    // (Optional) Print top-k scores or validate results here if needed
}
