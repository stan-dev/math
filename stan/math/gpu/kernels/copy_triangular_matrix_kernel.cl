R"(
#ifndef A
#define A(i, j)  A[j*rows+i]
#endif

#ifndef B
#define B(i, j)  B[j*rows+i]
#endif
__kernel void copy_triangular(
      __global double *A,
      __global double *B,
      unsigned int rows,
      unsigned int cols,
      unsigned int lower_upper) {
  int i = get_global_id(0);
  int j = get_global_id(1);
  if ( i < rows && j < cols ) {
    if ( !lower_upper && j <= i ) {
      A(i,j) = B(i,j);
    } else if ( !lower_upper ) {
      A(i,j) = 0;
    } else if ( lower_upper && j >= i ) {
      A(i,j) = B(i,j);
    } else if ( lower_upper && j < i ) {
      A(i,j) = 0;
    }
  }
};)"
