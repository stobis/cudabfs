NVCC=/usr/local/cuda/bin/nvcc
NVCCFLAGS=-arch sm_50 -O3 -std=c++11 --expt-extended-lambda

MGPU=3rdparty/moderngpu
MGPUFLAGS=-I $(MGPU)/src/

CXX=g++
CXXFLAGS=-O3 -std=c++11

all: bfs-mgpu.e dimacs-parser.e

bfs-mgpu.e: bfs-mgpu.cu
	$(NVCC) $(NVCCFLAGS) $(MGPUFLAGS) $< -o $@

dimacs-parser.e: dimacs-parser.cpp
	$(CXX) $(CXXFLAGS) $< -o $@

clean:
	rm -f *.e

.PHONY: all clean
