#include "label_propagation.h"
#include <unordered_map>
#include <unordered_set>
#include <numeric>
#include <iostream>
#include <chrono>
#include <sys/resource.h>
#include <unistd.h>

using namespace std;

/**
 * @brief Runs label propagation on the k-core of a hypergraph.
 *
 * This implementation operates on the in-core subset of vertices and builds a reduced
 * bipartite representation (`v2he2` and `he2v2`). It then iteratively propagates labels
 * based on majority voting among neighbors. The process repeats for `maxIter` rounds.
 *
 * @param v2he     Vertex-to-hyperedge adjacency list.
 * @param he2v     Hyperedge-to-vertex adjacency list.
 * @param inCore   Vector indicating whether a vertex is in the k-core (char: 0 or 1).
 * @param maxIter  Maximum number of iterations to perform.
 */
void runLabelPropagation(const vector<vector<int>>& v2he,
                         const vector<vector<int>>& he2v,
                         const vector<char>& inCore,
                         int maxIter) {
    int nv = (int)v2he.size();
    int nh = (int)he2v.size();

    // Step 1: Remap in-core vertices to new indices (compact form)
    vector<int> old2new(nv, -1);
    int nv2 = 0;
    for (int v = 0; v < nv; ++v) {
        if (inCore[v]) old2new[v] = nv2++;
    }

    // Step 2: Construct compact bipartite structure for in-core subgraph
    vector<vector<int>> v2he2(nv2);
    vector<vector<int>> he2v2;
    he2v2.reserve(nh);
    for (int e = 0; e < nh; ++e) {
        vector<int> memb2;
        for (int v : he2v[e]) {
            if (inCore[v]) memb2.push_back(old2new[v]);
        }
        if (!memb2.empty()) {
            int e2 = (int)he2v2.size();
            he2v2.push_back(memb2);
            for (int v2 : memb2)
                v2he2[v2].push_back(e2);
        }
    }

    // Step 3: Initialize labels (each vertex starts with its own ID)
    vector<int> label(nv2);
    iota(label.begin(), label.end(), 0); // label[i] = i

    // Step 4: Iterative label propagation
    for (int it = 0; it < maxIter; ++it) {
        bool changed = false;
        for (int v = 0; v < nv2; ++v) {
            unordered_map<int, int> freq;
            for (int e : v2he2[v]) {
                for (int u : he2v2[e]) {
                    // Boost label signal by repeating the vote (optional, here x10)
                    for (int k = 0; k < 10; ++k)
                        ++freq[label[u]];
                }
            }
            // Select the most frequent label (break ties by choosing smaller label)
            int bestLab = label[v], bestCnt = 0;
            for (auto& p : freq) {
                if (p.second > bestCnt || (p.second == bestCnt && p.first < bestLab)) {
                    bestLab = p.first;
                    bestCnt = p.second;
                }
            }
            if (bestLab != label[v]) {
                label[v] = bestLab;
                changed = true;
            }
        }
        // Optional: Early termination if no labels changed
        // if (!changed) break;
    }

    // Step 5: Summarize number of unique labels (communities)
    unordered_set<int> uniq(label.begin(), label.end());
    cout << "Label Propagation: " << uniq.size() << " communities found in core\n";
}
