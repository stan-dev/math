#ifndef STAN_MATH_PRIM_MAT_KERN_GPU_BASIC_MATRIX_KERNEL_HPP
#define STAN_MATH_PRIM_MAT_KERN_GPU_BASIC_MATRIX_KERNEL_HPP

#include <string>

#if defined(__APPLE__) || defined(__MACOSX)
#include <OpenCL/cl.hpp>
#else
#include <CL/cl.hpp>
#endif


/** @file stanmathcl/basic_matrix_kernels.hpp
    @brief 
*/

namespace stan {
  namespace math {
    namespace kernel_sources{

      std::string transpose= 
      "__kernel void transpose(\n"
      "          __global double *a,\n"
      "          __global double *b,\n"
      "          unsigned int M,\n"
      "          unsigned int N) \n"
      "{ \n"
      "    int i = get_global_id(0); \n"
      "    int j = get_global_id(1); \n"
      "    if( i < M && j < N ){ \n"
      "     a[j*M+i] = b[i*N+j];\n"
      "    } \n"
      "};\n";
      
      std::string copy= 
      "__kernel void copy(\n"
      "          __global double *a,\n"
      "          __global double *b,\n"
      "          unsigned int M,\n"
      "          unsigned int N) \n"
      "{ \n"
      "    int i = get_global_id(0); \n"
      "    int j = get_global_id(1); \n"
      "    if( i < M && j < N ){ \n"
      "     a[i*N+j] = b[i*N+j];\n"
      "    } \n"
      "};\n";  

      std::string zeros= 
      "__kernel void zeros(\n"
      "          __global double *a,\n"
      "          unsigned int M,\n"
      "          unsigned int N, \n"
      "          unsigned int part) \n"
      "{ \n"
      "    int i = get_global_id(0); \n"
      "    int j = get_global_id(1); \n"
      "    if( i < M && j < N ){ \n"
      "     if( part==0 && i<j ){ \n"
      "       a[i*N+j] = 0;\n"
      "     }else if( part==1 && i>j ){ \n"  
      "       a[i*N+j] = 0;\n"
      "     }else{ \n"
      "     } \n"
      "    } \n"
      "};\n";  
      
      std::string identity= 
      "__kernel void identity(\n"
      "          __global double *a,\n"
      "          unsigned int M,\n"
      "          unsigned int N) \n"
      "{ \n"
      "    int i = get_global_id(0); \n"
      "    int j = get_global_id(1); \n"
      "    if( i < M && j < N ){ \n"
      "       if( i==j ){ \n"
      "         a[i*N+j] = 1.0;\n"
      "       }else{ \n"
      "         a[i*N+j] = 0.0;\n"
      "       } \n"
      "    } \n"
      "};\n"; 
      
      std::string copy_triangular= 
      "__kernel void copy_triangular(\n"
      "          __global double *a,\n"
      "          __global double *b,\n"
      "          unsigned int M,\n"
      "          unsigned int N, \n"
      "          unsigned int lower_upper) \n"
      "{ \n"
      "    int i = get_global_id(0); \n"
      "    int j = get_global_id(1); \n"
      "    if( i < M && j < N){ \n"
      "     if( !lower_upper && i < M && j <=i ){ \n"
      "       a[i*N+j] = b[i*N+j];\n"
      "     }else if(!lower_upper && i<M && j<N){ \n"
      "       a[i*N+j]=0; \n"
      "     }else if(lower_upper && i<M && j>=i && j<N){ \n"
      "       a[i*N+j]=b[i*N+j];\n"
      "     }else if(lower_upper && i<M && j<i ){ \n"
      "       a[i*N+j]=0;\n"
      "     } \n"
      "    } \n"
      "};\n";  

      std::string copy_triangular_transposed= 
      "__kernel void copy_triangular_transposed(\n"
      "          __global double *a,\n"
      "          unsigned int M,\n"
      "          unsigned int N, \n"
      "          unsigned int lower_to_upper) \n"
      "{ \n"
      "    int i = get_global_id(0); \n"
      "    int j = get_global_id(1); \n"
      "    if( i < M && j < N ){ \n"
      "     if( lower_to_upper && j>i){ \n"
      "       a[j*N+i]=a[i*N+j];  \n"
      "     }else if( !lower_to_upper && j<i){ \n"
      "       a[j*N+i]=a[i*N+j];  \n"
      "     } \n"
      "    } \n"
      "};\n";
      
      std::string add= 
      "__kernel void add(\n"
      "          __global double *c,\n"
      "          __global double *a,\n"
      "          __global double *b,\n"
      "          unsigned int M,\n"
      "          unsigned int N) \n"
      "{ \n"
      "    int i = get_global_id(0); \n"
      "    int j = get_global_id(1); \n"
      "    if( i < M && j < N ){ \n"
      "     c[i*N+j] = a[i*N+j]+b[i*N+j];\n"
      "    } \n"
      "};\n"; 
      
      std::string subtract= 
      "__kernel void subtract(\n"
      "          __global double *c,\n"
      "          __global double *a,\n"
      "          __global double *b,\n"
      "          unsigned int M,\n"
      "          unsigned int N) \n"
      "{ \n"
      "    int i = get_global_id(0); \n"
      "    int j = get_global_id(1); \n"
      "    if( i < M && j < N ){ \n"
      "     c[i*N+j] = a[i*N+j]-b[i*N+j];\n"
      "    } \n"
      "};\n"; 
  
    }  
  }
}
#endif

