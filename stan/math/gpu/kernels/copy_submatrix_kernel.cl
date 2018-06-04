R"(
#define src(i,j) src[j*src_rows+i]
#define dst(i,j) dst[j*dst_rows+i]
__kernel void copy_submatrix(
      __global double *src,
      __global double *dst,
      unsigned int src_offset_i,
      unsigned int src_offset_j,
      unsigned int dst_offset_i,
      unsigned int dst_offset_j,
      unsigned int size_i,
      unsigned int size_j,
      unsigned int src_rows,
      unsigned int src_cols,
      unsigned int dst_rows,
      unsigned int dst_cols) {
  int i = get_global_id(0);
  int j = get_global_id(1);
  if ((i+src_offset_i) < src_rows && (j+src_offset_j) < src_cols
    && (i+dst_offset_i) < dst_rows && (j+dst_offset_j) < dst_cols) {
    dst((dst_offset_i+i), (dst_offset_j+j)) =
      src((src_offset_i+i),(src_offset_j+j));
  }
};)"
