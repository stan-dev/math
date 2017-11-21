R"=====(
#ifdef cl_khr_fp64
    #pragma OPENCL EXTENSION cl_khr_fp64 : enable
#elif defined(cl_amd_fp64)
    #pragma OPENCL EXTENSION cl_amd_fp64 : enable
#else
    #error "Double not supported by OpenCL implementation."
#endif

#define A(i,j)  A[i*cols+j]
#define AT(i,j)  A[i*rows+j]
#define B(i,j)  B[i*cols+j]
#define BT(i,j)  B[i*rows+j]
#define C(i,j)  C[i*cols+j]

__kernel void check_nan(
            __global double *A,
            int rows,
            int cols,
            __global int *flag) {
  const int i = get_global_id(0);
  const int j = get_global_id(1);
  if (isnan(A(i,j))) {
    flag[0] = 1;
  }
}

__kernel void check_diagonal_zeros(
              __global double *A,
              int rows,
              int cols,
              __global int *flag) {
  const int i = get_global_id(0);
  if (A(i,i) == 0) {
    flag[0] = 1;
  }
}

__kernel void check_symmetric(
            __global double *A,
            int rows,
            int cols,
            __global int *flag,
            double tolerance) {
  const int i = get_global_id(0);
  const int j = get_global_id(1);
  double diff = fabs(A(i,j)-A(j,i));
  if ( diff > tolerance ) {
    flag[0] = 1;
  }
}

)====="
