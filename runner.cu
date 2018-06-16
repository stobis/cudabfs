#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <queue>
#include <vector>
using namespace std;

#include "bfs-mgpu.cuh"

int main(int argc, char* argv[]) {
  if (argc != 2) {
    cerr << "Usage: " << argv[0] << " GRAPH" << endl;
    exit(1);
  }

  ifstream in(argv[1], ios::binary);  
  assert(in.is_open());
  int n, m;
  in.read((char*)&n, sizeof(int));
  in.read((char*)&m, sizeof(int));
  vector<int> nodes(n + 1), edges(m);
  in.read((char*)nodes.data(), nodes.size() * sizeof(int));
  in.read((char*)edges.data(), edges.size() * sizeof(int));

  for (int i = 0; i < 5; ++i) {
    int source = rand() % n;
    uint64_t seqsum = SequentialBFS(nodes, edges, source);
    uint64_t parsum = ParallelBFS(nodes, edges, source);
    assert(seqsum == parsum);
  }
}

