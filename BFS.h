#ifndef BFS_H
#define BFS_H

#include <vector>
#include "hypergraph.h"

/**
 * @file BFS.h
 * @brief Header for BFS-based graph algorithms on hypergraphs.
 *
 * This module provides utilities to run breadth-first search (BFS),
 * connected component analysis, and full-graph traversals on hypergraphs
 * represented either as flat arrays (CSR-style) or as adjacency lists.
 */

/**
 * @brief Run a BFS traversal starting from a fixed vertex (e.g., vertex 5050).
 *
 * Traverses the hypergraph using alternating vertex-hyperedge hops.
 *
 * @param nv       Number of vertices.
 * @param nh       Number of hyperedges.
 * @param degreeV  Vertex degrees.
 * @param edgesV   Vertex→hyperedge flattened edge list.
 * @param degreeH  Hyperedge degrees.
 * @param edgesH   Hyperedge→vertex flattened edge list.
 * @param offsetV  CSR-style offset array for vertex→hyperedges.
 * @param offsetH  CSR-style offset array for hyperedge→vertices.
 */
void runBFS(int nv, int nh, const int* degreeV, const int* edgesV,
            const int* degreeH, const int* edgesH,
            const int* offsetV, const int* offsetH);

/**
 * @brief Run BFS from every vertex in the decoded adjacency list.
 *
 * Used to evaluate the full traversal cost and reachability.
 *
 * @param adj Adjacency list representation of the hypergraph (bipartite).
 */
void runFullBFSAllVertices(const std::vector<std::vector<int>>& adj);

/**
 * @brief Compute connected components of the bipartite hypergraph.
 *
 * Uses alternating BFS on the vertex-hyperedge bipartite structure.
 * Outputs the CSR offsets used by BFS routines.
 *
 * @param G        Pointer to the hypergraph structure.
 * @param nv       Number of vertices.
 * @param nh       Number of hyperedges.
 * @param adj      Output adjacency list (decoded bipartite graph).
 * @param offsetV  Output CSR offset for vertex→hyperedges.
 * @param offsetH  Output CSR offset for hyperedge→vertices.
 */
void runConnectedComponents(Hypergraph* G, int nv, int nh,
                            std::vector<std::vector<int>>& adj,
                            std::vector<int>& offsetV,
                            std::vector<int>& offsetH);

#endif // BFS_H
