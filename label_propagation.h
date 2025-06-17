#ifndef LABEL_PROPAGATION_H
#define LABEL_PROPAGATION_H

#include <vector>

/**
 * @file label_propagation.h
 * @brief Header file for running label propagation on bipartite hypergraphs.
 *
 * This function performs label propagation on a hypergraph represented by
 * vertex-to-hyperedge and hyperedge-to-vertex adjacency lists.
 * 
 * Label propagation is a simple iterative algorithm used for community detection.
 * Initially, each node starts with a unique label. In each iteration, a node adopts
 * the most frequent label among its neighbors. The process is repeated for a
 * specified number of iterations or until convergence.
 *
 * The function restricts label updates to vertices and hyperedges marked as "inCore".
 */

/**
 * @brief Runs label propagation algorithm for community detection.
 *
 * @param v2he     Vertex-to-hyperedge adjacency list.
 * @param he2v     Hyperedge-to-vertex adjacency list.
 * @param inCore   Binary flag vector indicating whether a vertex/hyperedge is in the k-core.
 * @param maxIter  Maximum number of label propagation iterations.
 */
void runLabelPropagation(const std::vector<std::vector<int>>& v2he,
                         const std::vector<std::vector<int>>& he2v,
                         const std::vector<char>& inCore,
                         int maxIter);

#endif // LABEL_PROPAGATION_H
