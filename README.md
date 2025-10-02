# HybHuff: Compressing Hypergraphs of Skewed Topology through Huffman–Bitwise Coordination

This repository contains the code accompanying the paper:
HybHuff: Compressing Hypergraphs of Skewed Topology through Huffman–Bitwise Coordination

## Abstract

Hypergraphs provide a natural representation for many-to-many relationships in data-intensive applications, yet their practicality is limited by high memory consumption from explicitly storing vertices, hyperedges, and all incidences. This challenge is exacerbated by the skewed topology common in real-world hypergraphs. We propose HybHuff, a hybrid compression framework that coordinates Huffman encoding for frequent values with fixed-width bitwise encoding for the remainder. This coordination avoids large Huffman-tree overheads while exploiting frequency skew that bitwise packing alone would miss. We provide a theoretical analysis proving that, under mild conditions, an optimal allocation between the two sub-compressors exists, and we design a practical algorithm to approximate this allocation. We further introduce a blocked decoding scheme that enables random access to compressed hypergraphs and keeps runtime memory low. Experiments on real-world datasets show that HybHuff outperforms general-purpose compressors—achieving up to 2.3× smaller size than Zip and 1.9× than ZFP—and hypergraph-specialized methods—1.6× smaller than Hygra and 1.4× than HyperCSA. To assess real-world utility, we run Breadth-First Search, PageRank, and k-core directly over compressed data; across four benchmark datasets, HybHuff delivers competitive runtime while markedly reducing memory footprint.

## Quick Start
make clean && make
./main <percentage> Hypergraphs/<dataset>
For example:
./main 1.27 Hypergraphs/com-lj.all.cmty-hygra 
