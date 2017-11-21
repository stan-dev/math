R"=====(
#ifdef cl_khr_fp64
    #pragma OPENCL EXTENSION cl_khr_fp64 : enable
#elif defined(cl_amd_fp64)
    #pragma OPENCL EXTENSION cl_amd_fp64 : enable
#else
    #error "Double not supported by OpenCL implementation."
#endif

#define V(i,j) V[i*n+j]
#define d(i,j) d[i*n+j]

__kernel void cholesky_block(
                __global double *V,
                __global double *d,
                int n) {
  int f = get_local_id(0);
  for (int j = 0; j <n; j++) { 
      double ss = 0;
      for (int k = 0; k < j; k++) {
          ss += V(j,k) * V(j,k);
      }
      V(j,j) = sqrt(d(j,j) - ss);
      barrier(CLK_LOCAL_MEM_FENCE);
      int i=f;
      if ( i >= (j+1) && i < n ){
          double s = 0;
          for (int k = 0; k < j; k++) {
              s += V(i,k) * V(j,k);
          }
          V(i,j) = (1.0 / V(j,j) * (d(i,j) - s));
      }
  }
}
)====="
