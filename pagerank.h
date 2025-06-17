#ifndef PAGERANK_H
#define PAGERANK_H

#include <vector>

/**
 * @file pagerank.h
 * @brief Header file for hypergraph PageRank implementation.
 *
 * This module defines the interface for running the PageRank algorithm
 * on a bipartite hypergraph. The hypergraph is represented as a pair of
 * adjacency lists:
 *   - v2he: vertex-to-hyperedge connectivity
 *   - he2v: hyperedge-to-vertex connectivity
 *
 * PageRank is iteratively computed over both sides, with score updates
 * alternating between vertices and hyperedges.
 *
 */

/**
 * @brief Run PageRank on a bipartite hypergraph.
 *
 * @param v2he  Adjacency list from vertices to hyperedges
 * @param he2v  Adjacency list from hyperedges to vertices
 * @param iters Number of PageRank iterations to perform
 */
void runPageRank(const std::vector<std::vector<int>>& v2he,
                 const std::vector<std::vector<int>>& he2v,
                 int iters);

#endif // PAGERANK_H
