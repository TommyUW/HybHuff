#ifndef CORE_DECOMPOSITION_H
#define CORE_DECOMPOSITION_H

#include <vector>

/**
 * @file core_decomposition.h
 * @brief Header file for k-core decomposition on bipartite hypergraphs.
 *
 * The k-core of a graph is the maximal subgraph in which all vertices
 * have degree at least k. In the case of a hypergraph (represented as
 * bipartite), this function computes the vertices that belong to the
 * k-core by iteratively removing vertices and updating incident
 * hyperedges.
 *
 * The result is used for downstream algorithms such as label propagation.
 */

/**
 * @brief Compute the k-core of a hypergraph.
 *
 * This function identifies all vertices that belong to the k-core,
 * defined as the subgraph where every vertex has at least `k_thresh`
 * incident hyperedges.
 *
 * @param v2he     Vertex-to-hyperedge adjacency list.
 * @param he2v     Hyperedge-to-vertex adjacency list.
 * @param k_thresh Minimum degree threshold (k).
 * @return         A binary vector (`char` type), where 1 means in-core and 0 means pruned.
 */
std::vector<char> computeKCore(const std::vector<std::vector<int>>& v2he,
                               const std::vector<std::vector<int>>& he2v,
                               int k_thresh);

#endif // CORE_DECOMPOSITION_H
