#ifndef STAN_MATH_PRIM_MAT_FUN_VIENNACL_HPP
#define STAN_MATH_PRIM_MAT_FUN_VIENNACL_HPP

#ifndef VIENNACL_WITH_EIGEN
  #define VIENNACL_WITH_EIGEN 1
#endif

#include "viennacl/vector.hpp"
#include "viennacl/matrix.hpp"
#include "viennacl/linalg/direct_solve.hpp"
#include "viennacl/linalg/lu.hpp"
#include "viennacl/linalg/ilu_operations.hpp"

//a custom openCL kernel for copying the lower triangle to the upper
static const char * custom_kernel_lower_upper =
"__kernel void copy_lower_tri_upper_tri(\n"
"          __global double * vec1,\n"
"          unsigned int M,\n"
"          unsigned int N, \n"
"          unsigned int Npad) \n"
"{ \n"
"    int i = get_global_id(0); \n"
"    int j = get_global_id(1); \n"
"    vec1[i*Npad+j] = vec1[j*Npad+i] * (i < M && j < N && i < j);\n"
"};\n";
				
				
viennacl::ocl::program & my_prog = viennacl::ocl::current_context().add_program(custom_kernel_lower_upper, "custom_kernel_lower_upper");

viennacl::ocl::kernel & my_kernel_mul = my_prog.get_kernel("copy_lower_tri_upper_tri");

static const char * custom_kernel_lower_tri_mult = 
"__kernel void lower_triangular_matmul(const unsigned int L, const unsigned int M, \n"
"__global unsigned int *a, __global unsigned int *b, __global unsigned int *c) { \n"
"   int i, j, k, bl, di; \n"
"   i = get_group_id(1) * get_local_size(1) + get_local_id(1); \n"
"   j = get_group_id(0) * get_local_size(0) + get_local_id(0); \n"
"   // This gets the bottom left of the work group? \n"
"   bl = get_group_id(1) * get_local_size(1) + L; \n"
"   di = j * (L + 1); \n"
"   unsigned int temp = 0; \n"
"   // If bottom left of block is not in lower tri, exit \n"
"   if (bl > di) { \n"
"   } else { \n"
"   // normal matmul \n"
"       for(k = 0; k<M; k++) { \n"
"           temp += a[j * M + k]* b[k * L + i] * (i < L && j < M); \n"
"       } \n"
"       c[j * L + i] = temp; \n"
"   } \n"
"} \n";

viennacl::ocl::program & my_lower_tri_mult = viennacl::ocl::current_context().add_program(custom_kernel_lower_tri_mult, "custom_kernel_lower_tri_mult");

viennacl::ocl::kernel & my_kernel_lower_tri_mul = my_lower_tri_mult.get_kernel("lower_triangular_matmul");


 
#endif
