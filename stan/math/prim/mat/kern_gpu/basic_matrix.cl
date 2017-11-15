R"=====(
#ifdef cl_khr_fp64
    #pragma OPENCL EXTENSION cl_khr_fp64 : enable
#elif defined(cl_amd_fp64)
    #pragma OPENCL EXTENSION cl_amd_fp64 : enable
#else
    #error "Double not supported by OpenCL implementation."
#endif

__kernel void transpose(
          __global double *a,
          __global double *b,
          int M,
          int N ) {
    int i = get_global_id(0);
    int j = get_global_id(1);
    if ( i < M && j < N ) {
     a[j*M+i] = b[i*N+j];     
    }
}

__kernel void copy(
          __global double *a,
          __global double *b,
          unsigned int M,
          unsigned int N) {
    int i = get_global_id(0);
    int j = get_global_id(1);
    if ( i < M && j < N ) {
     b[i*N+j] = a[i*N+j];
    }
}

__kernel void zeros(
          __global double *a,
          unsigned int M,
          unsigned int N,
          unsigned int part) {
    int i = get_global_id(0);
    int j = get_global_id(1);
    if ( i < M && j < N ) {
     if ( part == 0 && j < i ) {
       a[i*N+j] = 0;
     } else if ( part == 1 && j > i ) {
       a[i*N+j] = 0;
     } else if ( part == 2 ) {
       a[i*N+j] = 0;
     }
    }
}

__kernel void identity(
          __global double *a,
          unsigned int M,
          unsigned int N) {
    int i = get_global_id(0);
    int j = get_global_id(1);
    if ( i < M && j < N ) {
       if ( i == j ) {
         a[i*N+j] = 1.0;
       } else {
         a[i*N+j] = 0.0;
       }
    }
}

__kernel void copy_triangular(
          __global double *a,
          __global double *b,
          unsigned int M,
          unsigned int N,
          unsigned int lower_upper) {
    int i = get_global_id(0);
    int j = get_global_id(1);
    if ( i < M && j < N ) {
     if ( !lower_upper && j <= i ) {
       a[i*N+j] = b[i*N+j];
     } else if ( !lower_upper ) {
       a[i*N+j] = 0;
     } else if ( lower_upper && j >= i ) {
       a[i*N+j] = b[i*N+j];
     } else if ( lower_upper && j < i ) {
       a[i*N+j] = 0;
     }
    }
}

__kernel void copy_triangular_transposed(
          __global double *a,
          unsigned int M,
          unsigned int N,
          unsigned int lower_to_upper) {
    int i = get_global_id(0);
    int j = get_global_id(1);
    if ( i < M && j < N ) {
     if ( lower_to_upper && j > i ) {
       a[j*N+i] = a[i*N+j];
     } else if ( !lower_to_upper && j > i ) {
       a[i*N+j] = a[j*N+i];
     }
    }
}

__kernel void add(
          __global double *c,
          __global double *a,
          __global double *b,
          unsigned int M,
          unsigned int N) {
    int i = get_global_id(0);
    int j = get_global_id(1);
    if ( i < M && j < N ) {
     c[i*N+j] = a[i*N+j]+b[i*N+j];
    }
}

__kernel void subtract(
          __global double *c,
          __global double *a,
          __global double *b,
          unsigned int M,
          unsigned int N) {
    int i = get_global_id(0);
    int j = get_global_id(1);
    if ( i < M && j < N ) {
     c[i*N+j] = a[i*N+j] - b[i*N+j];
    }
}

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
      dst[(dst_offset_i+i)*dst_cols+dst_offset_j+j] =
          src[(src_offset_i+i)*src_cols+src_offset_j+j];
    }
}


)====="
