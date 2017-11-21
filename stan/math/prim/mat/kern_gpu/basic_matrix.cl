R"=====(
#define A(i,j)  A[i*cols+j]
#define AT(i,j)  A[i*rows+j]
#define B(i,j)  B[i*cols+j]
#define BT(i,j)  B[i*rows+j]
#define C(i,j)  C[i*cols+j]

#ifdef cl_khr_fp64
    #pragma OPENCL EXTENSION cl_khr_fp64 : enable
#elif defined(cl_amd_fp64)
    #pragma OPENCL EXTENSION cl_amd_fp64 : enable
#else
    #error "Double not supported by OpenCL implementation."
#endif


__kernel void transpose(
          __global double *B,
          __global double *A,
          int rows,
          int cols ) {
    int i = get_global_id(0);
    int j = get_global_id(1);
    if ( i < rows && j < cols ) {
     BT(j,i) = A(i,j);     
    }
}

__kernel void copy(
          __global double *A,
          __global double *B,
          unsigned int rows,
          unsigned int cols) {
    int i = get_global_id(0);
    int j = get_global_id(1);
    if ( i < rows && j < cols ) {
     B(i,j) = A(i,j);
    }
}

__kernel void zeros(
          __global double *A,
          unsigned int rows,
          unsigned int cols,
          unsigned int part) {
    int i = get_global_id(0);
    int j = get_global_id(1);
    if ( i < rows && j < cols ) {
     if ( part == 0 && j < i ) {
       A(i,j) = 0;
     } else if ( part == 1 && j > i ) {
       A(i,j) = 0;
     } else if ( part == 2 ) {
       A(i,j) = 0;
     }
    }
}

__kernel void identity(
          __global double *A,
          unsigned int rows,
          unsigned int cols) {
    int i = get_global_id(0);
    int j = get_global_id(1);
    if ( i < rows && j < cols ) {
       if ( i == j ) {
         A(i,j) = 1.0;
       } else {
         A(i,j) = 0.0;
       }
    }
}

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
}

__kernel void copy_triangular_transposed(
          __global double *A,
          unsigned int rows,
          unsigned int cols,
          unsigned int lower_to_upper) {
    int i = get_global_id(0);
    int j = get_global_id(1);
    if ( i < rows && j < cols ) {
     if ( lower_to_upper && j > i ) {
       AT(j,i) = A(i,j);
     } else if ( !lower_to_upper && j > i ) {
       A(i,j) = AT(j,i);
     }
    }
}

__kernel void add(
          __global double *C,
          __global double *A,
          __global double *B,
          unsigned int rows,
          unsigned int cols) {
    int i = get_global_id(0);
    int j = get_global_id(1);
    if ( i < rows && j < cols ) {
      C(i,j) = A(i,j) + B(i,j);
    }
}

__kernel void subtract(
          __global double *C,
          __global double *A,
          __global double *B,
          unsigned int rows,
          unsigned int cols) {
    int i = get_global_id(0);
    int j = get_global_id(1);
    if ( i < rows && j < cols ) {
     C(i,j) = A(i,j)-B(i,j);
    }
}

#define src(i,j) src[i*src_cols+j]
#define dst(i,j) dst[i*dst_cols+j]

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
    if ( (i+src_offset_i) < src_rows && (j+src_offset_j) < src_cols
      && (i+dst_offset_i) < dst_rows && (j+dst_offset_j) < dst_cols ) {
      dst((dst_offset_i+i), (dst_offset_j+j)) =
          src((src_offset_i+i),(src_offset_j+j));
    }
}


)====="
