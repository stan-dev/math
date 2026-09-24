#ifndef STAN_MATH_OPENCL_REV_SD_HPP
#define STAN_MATH_OPENCL_REV_SD_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/rev/arena_matrix_cl.hpp>
#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/opencl/prim/mean.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/value_of.hpp>
#include <stan/math/opencl/rev/scalar_cl.hpp>

namespace stan {
namespace math {

/**
 * Return the sample standard deviation of the var_value matrix
 *
 * @tparam T Input type
 * @param[in] A input matrix
 * @return sample standard deviation of specified matrix
 * @throw domain error size is not greater than zero.
 */
template <typename T,
          require_all_kernel_expressions_and_none_scalar_t<T>* = nullptr>
inline opencl::ScalarCl<var> sd(const var_value<T>& A) {
  if (A.size() == 1) {
    return opencl::ScalarCl<var>();
  }
  opencl::ScalarCl<double> A_mean = mean(A.val());
  arena_matrix_cl<double> diff;
  matrix_cl<double> sq_norm;
  auto diff_expr = A.val() - A_mean;
  results(diff, sq_norm) = expressions(diff_expr, sum_2d(square(diff_expr)));
  return opencl::make_callback_scalar_cl(
      sqrt(sum(sq_norm) / (A.size() - 1.0)),
      [A, diff](const auto& res_adj, const auto& res_val) mutable {
        opencl::ScalarCl<double> factor
            = res_adj / (res_val * (A.size() - 1.0));
        A.adj() += factor * diff;
      });
}

}  // namespace math
}  // namespace stan

#endif
#endif
