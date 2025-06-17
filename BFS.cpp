// ==========================================================
// File: BFS.cpp
//
// Description:
// Implements breadth-first search (BFS) and connected components
// traversal on bipartite hypergraphs represented using flat arrays
// and adjacency lists. Designed to work with compressed hypergraph
// representations.
//
// Algorithms:
// - Connected Component Detection using alternating vertex ↔ hyperedge BFS
// - Vertex-centric BFS from a fixed root
// - All-pairs BFS from all vertices to estimate global reachability
//
// Author: Tianyu Zhao, University of Washington
// Advisors: Prof. Dongfang Zhao, Dr. Luanzheng Guo, Dr. Nathan Tallent
// ==========================================================

#include "BFS.h"
#include <iostream>
#include <queue>
#include <chrono>
#include <vector>
#include <cassert>
#include <algorithm>
using namespace std;

// ------------------------------------------------------------
// Connected Component Detection (Bipartite BFS)
//
// Given the hypergraph in flattened CSR format (degree + edge list),
// this function performs BFS alternately across vertices and
// hyperedges. Each connected component is labeled with an ID.
// This traversal is used to build component-aware offsets for BFS.
// ------------------------------------------------------------
void runConnectedComponents(Hypergraph* G, int nv, int nh,
                            vector<vector<int>>& adj,
                            vector<int>& offsetV,
                            vector<int>& offsetH) {
    auto t_cc0 = chrono::high_resolution_clock::now();

    // Compute CSR-style offsets for vertex and hyperedge neighbors
    offsetV.assign(nv + 1, 0);
    for (int v = 0; v < nv; ++v)
        offsetV[v + 1] = offsetV[v] + G->degreeV[v];
    assert(offsetV[nv] == G->mv);

    offsetH.assign(nh + 1, 0);
    for (int e = 0; e < nh; ++e)
        offsetH[e + 1] = offsetH[e] + G->degreeH[e];
    assert(offsetH[nh] == G->mh);

    // Initialize visitation and component ID markers
    vector<char> visitedV(nv, 0), visitedE(nh, 0);
    vector<int> compIDV(nv, -1), compIDE(nh, -1);
    queue<pair<bool, int>> bfsQ;

    int compLabel = 0;

    // Start BFS from each unvisited vertex
    for (int v0 = 0; v0 < nv; ++v0) {
        if (visitedV[v0]) continue;

        visitedV[v0] = 1;
        compIDV[v0] = compLabel;
        bfsQ.push({false, v0});  // (isHyperedge=false, vertexID)

        // Perform bipartite BFS traversal (vertex ↔ hyperedge)
        while (!bfsQ.empty()) {
            auto [isHE, id] = bfsQ.front();
            bfsQ.pop();

            if (!isHE) {  // vertex → hyperedge neighbors
                int base = offsetV[id];
                for (int i = 0; i < G->degreeV[id]; ++i) {
                    int e = G->edgesV[base + i];
                    if (e < 0 || e >= nh || visitedE[e]) continue;
                    visitedE[e] = 1;
                    compIDE[e] = compLabel;
                    bfsQ.push({true, e});
                }
            } else {  // hyperedge → vertex neighbors
                int base = offsetH[id];
                for (int i = 0; i < G->degreeH[id]; ++i) {
                    int u = G->edgesH[base + i];
                    if (u < 0 || u >= nv || visitedV[u]) continue;
                    visitedV[u] = 1;
                    compIDV[u] = compLabel;
                    bfsQ.push({false, u});
                }
            }
        }

        compLabel++;  // new component
    }

    auto t_cc1 = chrono::high_resolution_clock::now();
    cout << "Connected Components time: "
         << chrono::duration<double>(t_cc1 - t_cc0).count() << " s\n";
}

// ------------------------------------------------------------
// Single-source BFS from fixed vertex (e.g., vertex 5050)
//
// Traverses bipartite hypergraph starting from a vertex,
// alternates vertex → hyperedge → vertex expansions.
// Measures distance from the root to all reachable vertices.
// ------------------------------------------------------------
void runBFS(int nv, int nh, const int* degreeV, const int* edgesV,
            const int* degreeH, const int* edgesH,
            const int* offsetV, const int* offsetH) {

    auto t_bfs0 = chrono::high_resolution_clock::now();
    vector<int> distB(nv, -1);
    queue<int> qB;

    int startV = 5050;  // source vertex
    if (startV >= 0 && startV < nv) {
        distB[startV] = 0;
        qB.push(startV);

        while (!qB.empty()) {
            int v = qB.front(); qB.pop();
            int d = distB[v];

            // vertex → hyperedges
            for (int i = 0; i < degreeV[v]; ++i) {
                int e = edgesV[offsetV[v] + i];
                if (e < 0 || e >= nh) continue;

                // hyperedge → vertices
                for (int j = 0; j < degreeH[e]; ++j) {
                    int u = edgesH[offsetH[e] + j];
                    if (u < 0 || u >= nv || distB[u] != -1) continue;
                    distB[u] = d + 1;
                    qB.push(u);
                }
            }
        }

        // Report stats
        int reached = count_if(distB.begin(), distB.end(), [](int d) { return d != -1; });
        cout << "\nBFS from vertex " << startV << " reached " << reached << " / " << nv << " vertices\n";

        int toShow = min(nv, 10);
        cout << "Distances (vertex → distance):\n";
        for (int v = 0; v < toShow; ++v) {
            cout << "  v=" << v << (distB[v] == -1 ? " unreachable\n"
                                                   : (" → " + to_string(distB[v]) + "\n"));
        }
    }

    auto t_bfs1 = chrono::high_resolution_clock::now();
    cout << "BFS runtime: " << chrono::duration<double>(t_bfs1 - t_bfs0).count() << " s\n";
}

// ------------------------------------------------------------
// Full BFS from every vertex
//
// Runs BFS starting from each vertex in the decoded adjacency list.
// Aggregates the total number of reachable vertices across runs.
// Used to simulate dense traversal workloads for comparison.
// ------------------------------------------------------------
void runFullBFSAllVertices(const vector<vector<int>>& adj) {
    auto t0 = chrono::high_resolution_clock::now();
    size_t N = adj.size();
    int totalReach = 0;

    for (size_t src = 0; src < N; ++src) {
        vector<char> vis(N, 0);
        queue<int> q;
        vis[src] = 1;
        q.push(src);
        int count = 1;

        // Classic unweighted BFS from `src`
        while (!q.empty()) {
            int u = q.front(); q.pop();
            for (int v : adj[u]) {
                if ((unsigned)v < N && !vis[v]) {
                    vis[v] = 1;
                    q.push(v);
                    ++count;
                }
            }
        }

        totalReach += count;
    }

    auto t1 = chrono::high_resolution_clock::now();
    cout << "All-vertex BFS finished. Total reach count: " << totalReach
         << " | Time: " << chrono::duration<double>(t1 - t0).count() << " s\n";
}
