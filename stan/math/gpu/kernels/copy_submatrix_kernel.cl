R"(
#ifndef src
#define src(i,j) src[j * src_rows + i]
#endif
#ifndef dst
#define dst(i,j) dst[j * dst_rows + i]
#endif
/**
 * Copies a submatrix of the source matrix to
 * the destination matrix. The submatrix to copy
 * starts at (src_offset_rows, src_offset_cols)
 * and is of size size_rows x size_cols.
 * The submatrix is copied to the
 * destination matrix starting at
 * (dst_offset_rows, dst_offset_cols)
 *
 * @param src (read) The source matrix.
 * @param dst (write) The destination submatrix.
 * @param src_offset_i The offset row in src.
 * @param src_offset_j The offset column in src.
 * @param dst_offset_i The offset row in dst.
 * @param dst_offset_j The offset column in dst.
 * @param size_i The number of rows in the submatrix.
 * @param size_j The number of columns in the submatrix.
 * @param src_rows The number of rows in the source matrix.
 * @param src_cols The number of cols in the source matrix.
 * @param src_rows The number of rows in the destination matrix.
 * @param dst_cols The number of cols in the destination matrix.
 *
 * @note used in math/gpu/copy_submatrix_opencl.hpp
 */
__kernel void copy_submatrix(__global double *src, __global double *dst,
  unsigned int src_offset_i, unsigned int src_offset_j,
  unsigned int dst_offset_i, unsigned int dst_offset_j,
  unsigned int size_i, unsigned int size_j,
  unsigned int src_rows, unsigned int src_cols,
  unsigned int dst_rows, unsigned int dst_cols) {
  int i = get_global_id(0);
  int j = get_global_id(1);
  if ((i + src_offset_i) < src_rows && (j + src_offset_j) < src_cols &&
   (i + dst_offset_i) < dst_rows && (j + dst_offset_j) < dst_cols) {
    dst((dst_offset_i + i), (dst_offset_j + j)) =
      src((src_offset_i + i),(src_offset_j + j));
  }
};)"
