#ifndef STAN_MATH_OPENCL_PRIM_CUMULATIVE_SUM_HPP
#define STAN_MATH_OPENCL_PRIM_CUMULATIVE_SUM_HPP
#ifdef STAN_OPENCL
#include <stan/math/opencl/prim/size.hpp>
#include <stan/math/opencl/kernels/scan.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>

namespace stan {
namespace math {

/**
 * Returns a Matrix with values added along the main diagonal
 *
 * @tparam T_m type of input kernel generator expression for the input matrix
 * @tparam T_a type of input kernel generator expression to add along the
 * diagonal
 *
 * @param mat input kernel generator expression
 * @param to_add scalar value or input kernel generator expression to add along
 * the diagonal
 * @return a kernel generator expressio with to_add added along main diagonal
 */
template <typename T_vec,
          require_all_kernel_expressions_and_none_scalar_t<T_vec>* = nullptr>
inline auto cumulative_sum(T_vec&& v) {
  check_vector("cumulative_sum(OpenCL)", "v", v);
  const int local_size
      = opencl_kernels::scan_double_sum.get_option("LOCAL_SIZE_");
  int work_groups
      = opencl_context.device()[0].getInfo<CL_DEVICE_MAX_COMPUTE_UNITS>() * 16;
  matrix_cl<double> tmp_threads(local_size * work_groups, 1);
  matrix_cl<double> tmp_wgs(work_groups, 1);
  try {
    opencl_kernels::scan_double_sum(cl::NDRange(local_size * work_groups),
                                    cl::NDRange(local_size), tmp_wgs,
                                    tmp_threads, v, v.size());
  } catch (const cl::Error& e) {
    check_opencl_error("cumulative_sum", e);
  }
  return tmp_wgs;
}
}  // namespace math
}  // namespace stan

#endif
#endif
