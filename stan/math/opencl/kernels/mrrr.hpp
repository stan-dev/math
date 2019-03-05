#ifndef STAN_MATH_GPU_KERNELS_MRRR_HPP
#define STAN_MATH_GPU_KERNELS_MRRR_HPP
#ifndef STAN_OPENCL
#error "NO STAN_OPENCL"
#endif
#ifdef STAN_OPENCL

#include <stan/math/opencl/kernel_cl.hpp>

namespace stan {
namespace math {
namespace opencl_kernels {

// \cond
const char* eigenvals_bisect_kernel_code = STRINGIFY(
// \endcond
        int getSturmCountLdl(__global double* l, const __global double* d, double shift, int n){
          double s = -shift;
          double l_plus;
          double d_plus;
          int count = 0;
          for (int i = 0; i < n; i++) {
            d_plus = s + d[i];
            count += d_plus >= 0;
            if (isinf(d_plus) && isinf(s)) { // this happens if d_plus==0 -> in next iteration d_plus==inf and s==inf
              s = l[i] * l[i] * d[i] - shift;
            }
            else {
              s = l[i] * l[i] * s * (d[i] / d_plus) - shift;
            }
          }
          d_plus = s + d[n];
          count += d_plus >= 0;
          return count;
        }

        __kernel void eigenvals_bisect(__global double* l, const __global double* d, __global double* low_res, __global double* high_res,
                const double min_eigval, const double max_eigval, const int n) {
          const int lid = get_local_id(0);
          const int gid = get_global_id(0);
          const int gsize = get_global_size(0);
          const int lsize = get_local_size(0);
          const int ngroups = get_num_groups(0);
          const int wgid = get_group_id(0);

          double eps = 3e-16;
          double min_norm = 3e-308; //(approximately) smallest normalized double, larger than 0

          int i=gid;
          double low = min_eigval;
          double high = max_eigval;

          while (fabs((high - low) / (high + low)) > eps && fabs(high - low) > min_norm) {
            double mid = (high + low) * 0.5;
            int count = getSturmCountLdl(l, d, mid, n-1);
            if (count > i) {
              low = mid;
            }
            else {
              high = mid;
            }
          }
          low_res[i]=low;
          high_res[i]=high;
        }
// \cond
);
// \endcond

const global_range_kernel<cl::Buffer, cl::Buffer, cl::Buffer, cl::Buffer, double, double, int >
        eigenvals_bisect("eigenvals_bisect", eigenvals_bisect_kernel_code);

}  // namespace opencl_kernels
}  // namespace math
}  // namespace stan
#endif
#endif