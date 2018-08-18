#define STRINGIFY(src) #src

STRINGIFY(
    /**
     * Takes the transpose of the matrix on the GPU.
     *
     * @param[out] B The output matrix to hold transpose of A.
     * @param[in] A The input matrix to transpose into B.
     * @param rows The number of rows for A.
     * @param cols The number of columns for A.
     *
     * @note This kernel uses the helper macros available in helpers.cl.
     */
    __kernel void transpose(__global read_write double *B,
                            __global read_only double *A,
                            read_only unsigned int rows,
                            read_only unsigned int cols) {
      int i = get_global_id(0);
      int j = get_global_id(1);
      if (i < rows && j < cols) {
        BT(j, i) = A(i, j);
      }
    })
