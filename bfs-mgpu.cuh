#ifndef BFS_MGPU_CUH
#define BFS_MGPU_CUH

#include <moderngpu/context.hxx>
#include <moderngpu/memory.hxx>
using namespace mgpu;

#include <vector>
using namespace std;

void ParallelBFS(int n, int m, mem_t<int>& nodes, mem_t<int>& edges, int source,
                 mem_t<int>& distance, context_t& context);

uint64_t ParallelBFS(const vector<int>& nodes, const vector<int>& edges,
                     int source);

uint64_t SequentialBFS(const vector<int>& nodes, const vector<int>& edges, 
                        int source);

#endif  // BFS_MGPU_CUH
