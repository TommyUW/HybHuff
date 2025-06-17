
CXX = g++
CXXFLAGS = -std=c++17 -O3
OBJS = main.o hypergraph.o huffman_tree.o huffman_code.o encode.o decode.o BFS.o computeKCore.o label_propagation.o pagerank.o

all: main

main: $(OBJS)
	$(CXX) $(CXXFLAGS) -o main $(OBJS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f *.o main