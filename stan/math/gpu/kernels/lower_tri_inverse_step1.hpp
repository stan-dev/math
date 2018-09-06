#ifndef STAN_MATH_GPU_KERNELS_LOWER_TRI_INVERSE_STEP1_HPP
#define STAN_MATH_GPU_KERNELS_LOWER_TRI_INVERSE_STEP1_HPP
#ifdef STAN_OPENCL

#include <stan/math/gpu/kernel_cl.hpp>

namespace stan {
namespace math {
namespace opencl_kernels {
// \cond
const char* lower_tri_inverse_step1_kernel_code = STRINGIFY(
    // \endcond
    /**
     * Calculates 
     *
     * @param[in,out] A The identity matrix output.
     * @param rows The number of rows for A.
     * @param cols The number of cols for A.
     * @note Code is a <code>const char*</code> held in
     * <code>identity_kernel_code.</code>
     *  Used in math/gpu/identity_opencl.hpp.
     *  This kernel uses the helper macros available in helpers.cl.
     */
    __kernel void lower_tri_inverse_step1(
        __global double* A,
        __global double* V,
        int remainder,
        int part_size_fixed,
        int rows) {
      int indeks = get_global_id(0);
      int i = indeks*part_size_fixed;
      int part_size;
      double faktor;

      if ( indeks < remainder ) {
        i += indeks;
        part_size = part_size_fixed+1;
      } else {
        i += remainder;
        part_size = part_size_fixed;
      }
      int offset = indeks*(part_size_fixed+1)*(part_size_fixed+1);
      for (int p = 0; p < part_size; p++) {
        for (int r = 0; r < part_size; r++) {
          if ( p == r ) {
            V(p,r) = 1;
          } else {
            V(p,r) = 0;
          }
        }
      }
      for (unsigned int ii = 0; ii < part_size; ii++) {
        if ( ii > 0 ) {
          for (unsigned int j = ii; j < part_size; j++) {
            faktor = A((j+i),(i+ii-1));
            for (unsigned int k = 0; k < part_size; k++) {
              V(j,k) -= faktor*V((ii-1),k);
            }
          }
        }
        faktor = A((ii+i),(ii+i));
        for (unsigned int k = 0; k < part_size; k++) {
          V(ii,k) /= faktor;
        }
      }
      for (int p = 0; p < part_size; p++) {
        for (int r=0; r < part_size; r++) {
          A((p+i),(i+r)) = V(p,r);
        }
      }
  }
    // \cond
);
// \endcond

/**
 * See the docs for \link kernels/matrix_multiply.hpp add() \endlink
 */
const local_range_kernel<cl::Buffer, cl::Buffer, int, int, int>
    lower_tri_inverse_step1("lower_tri_inverse_step1", lower_tri_inverse_step1_kernel_code);

}  // namespace opencl_kernels
}  // namespace math
}  // namespace stan
#endif
#endif
