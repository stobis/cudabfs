#ifndef BFS_MGPU_CUH
#define BFS_MGPU_CUH

#include <moderngpu/context.hxx>
#include <moderngpu/memory.hxx>

#include <vector>

namespace bfs_mgpu {

void ParallelBFS(int, int, mgpu::mem_t<int>&, mgpu::mem_t<int>&, int,
                 mgpu::mem_t<int>&, mgpu::context_t&);

uint64_t ParallelBFS(const std::vector<int>&, const std::vector<int>&, int);

uint64_t SequentialBFS(const std::vector<int>&, const std::vector<int>&, int);

}  // namespace bfs_mgpu

#endif  // BFS_MGPU_CUH
