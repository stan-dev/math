#ifndef STAN_MATH_PRIM_MAT_KERN_GPU_MATRIX_MULTIPLY_KERNEL_HPP
#define STAN_MATH_PRIM_MAT_KERN_GPU_MATRIX_MULTIPLY_KERNEL_HPP

#include <string>

#if defined(__APPLE__) || defined(__MACOSX)
#include <OpenCL/cl.hpp>
#else
#include <CL/cl.hpp>
#endif


/** @file stan/math/prim/mat/kern_gpu/matrix_multiply_kernels.hpp
    @brief 
*/

namespace stan {
  namespace math {
    namespace kernel_sources{

      static std::string scalar_mul_diagonal= 
      "__kernel void scalar_mul_diagonal(\n"
      "          __global double *a,\n"
      "          double scalar,\n"
      "          unsigned int M,\n"
      "          unsigned int N) \n"
      "{ \n"
      "    int i = get_global_id(0); \n"
      "    if (i < M && i < N){ \n"
      "     a[i*N+i] *= scalar  ;\n"
      "    } \n"
      "};\n";
      
      static std::string scalar_mul= 
      "__kernel void scalar_mul(\n"
      "          __global double *a,\n"
      "          double scalar,\n"
      "          unsigned int M,\n"
      "          unsigned int N) \n"
      "{ \n"
      "    int i = get_global_id(0); \n"
      "    int j = get_global_id(1); \n"
      "    if (i < M && j < N){ \n"
      "     a[i*N+j]*= scalar; \n"
      "    } \n"
      "};\n";

      static std::string basic_multiply=
      " #define TS 16 \n"
      "__kernel void basic_multiply(const int M, const int N, const int K, \n"
      "                      const __global double* A, \n"
      "                      const __global double* B, \n"
      "                      __global double* C) { \n"
      "      \n"
      "    // Thread identifiers \n"
      "    const int row = get_local_id(0); // Local row ID (max: TS)  \n"
      "    const int col = get_local_id(1); // Local col ID (max: TS)  \n"
      "    const int globalRow = TS*get_group_id(0) + row; // Row ID of C (0..M) \n"
      "    const int globalCol = TS*get_group_id(1) + col; // Col ID of C (0..N) \n"
      "  \n"
      "    // Local memory to fit a tile of TS*TS elements of A and B  \n"
      "    __local double Asub[TS][TS];  \n"
      "    __local double Bsub[TS][TS];  \n"
      "  \n"
      "    // Initialise the accumulation register \n"
      "    double acc = 0.0; \n"
      "    \n"
      "    // Loop over all tiles \n"
      "    const int numTiles = (K+TS-1)/TS;\n"
      "    for (int t=0; t<numTiles; t++) {\n"
      " \n"
      "        // Load one tile of A and B into local memory \n"
      "        const int tiledRow = TS*t + row; \n"
      "        const int tiledCol = TS*t + col; \n"
      "        if (tiledRow<K && globalCol<N){  \n"
      "         Asub[col][row] = B[tiledRow*N + globalCol];  \n"
      "        }else{\n"
      "         Asub[col][row] = 0.0;  \n"
      "        }\n"
      "        if (tiledCol<K && globalRow<M){  \n"
      "         Bsub[col][row] = A[globalRow*K + tiledCol];  \n"
      "        }else{\n"
      "         Bsub[col][row] = 0.0;  \n"
      "        }\n"
      " \n"
      "        // Synchronise to make sure the tile is loaded \n"
      "        barrier(CLK_LOCAL_MEM_FENCE);  \n"
      " \n"
      "        // Perform the computation for a single tile \n"
      "        for (int k=0; k<TS; k++) { \n"
      "            acc += Bsub[k][row] * Asub[col][k];  \n"
      "        }  \n"
      " \n"
      "        // Synchronise before loading the next tile  \n"
      "        barrier(CLK_LOCAL_MEM_FENCE);  \n"
      "    }  \n"
      " \n"
      "    // Store the final result in C \n"
      "    if (globalCol<N && globalRow<M){  \n"
      "     C[globalRow*N + globalCol] = acc;  \n"
      "    } \n"
      " } \n";

    }
  }
}

#endif

