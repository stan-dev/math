#define STRINGIFY(src) #src

STRINGIFY(
    /**
     * Check if the <code>matrix_gpu</code> has zeros on the diagonal
     *
     * @param[in] A Matrix to check.
     * @param[out] flag the flag to be written to if any diagonal is zero.
     * @param rows The number of rows for A.
     * @param cols The number of cols of A.
     *
     * @note Kernel for stan/math/gpu/err/check_diagonal_zeros.hpp.
     *  This kernel uses the helper macros available in helpers.cl.
     */
    __kernel void is_zero_on_diagonal(__global read_only double *A,
                                      __global int *flag,
                                      read_only unsigned int rows,
                                      write_only unsigned int cols) {
      const int i = get_global_id(0);
      if (i < rows && i < cols) {
        if (A(i, i) == 0) {
          flag[0] = 1;
        }
      }
    })
