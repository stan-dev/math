#ifndef STAN_MATH_OPENCL_PRIM_SORT_ASC_HPP
#define STAN_MATH_OPENCL_PRIM_SORT_ASC_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/kernels/mergesort.hpp>
#include <algorithm>

#include <stan/math/opencl/copy.hpp>
#include <stan/math/opencl/kernel_generator.hpp>

namespace stan {
namespace math {

/**
 * Sort the (row)vector in ascending order.
 * @param input vector of row vector to sort
 */
template <typename T,
          require_all_nonscalar_prim_or_rev_kernel_expression_t<T>* = nullptr>
inline auto sort_asc(T&& input) {
  using T_val = value_type_t<T>;
  check_vector("sort_asc(OpenCL)", "input", input);
  matrix_cl<T_val> a = std::forward<T>(input);
  matrix_cl<T_val> b(a.rows(), a.cols());
  matrix_cl<T_val>* in_ptr = &a;
  matrix_cl<T_val>* out_ptr = &b;
  int compute_units
      = opencl_context.device()[0].getInfo<CL_DEVICE_MAX_COMPUTE_UNITS>();
  int local = 64;
  for (int run_len = 1; run_len < input.size(); run_len *= 2) {
    int runs = (input.size() + run_len - 1) / run_len;
    opencl_kernels::sort_asc<T_val>::merge_step(
        cl::NDRange(4 * compute_units * local), cl::NDRange(local), *out_ptr,
        *in_ptr, run_len, input.size(), (runs + 1) / 2);
    std::swap(in_ptr, out_ptr);
  }
  return *in_ptr;
}
}  // namespace math
}  // namespace stan
#endif
#endif
