#ifndef STAN_MATH_OPENCL_REV_SD_HPP
#define STAN_MATH_OPENCL_REV_SD_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/rev/arena_matrix_cl.hpp>
#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/opencl/prim/mean.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/value_of.hpp>

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
inline var sd(const var_value<T>& A) {
  if (A.size() == 1) {
    return 0.0;
  }
  double A_mean = mean(A.val());
  arena_matrix_cl<double> diff;
  matrix_cl<double> sq_norm;
  auto diff_expr = A.val() - A_mean;
  auto sq_norm_expr = sum_2d(square(diff));
  results(diff, sq_norm) = expressions(diff_expr, sq_norm_expr);

  return make_callback_var(
      sqrt(from_matrix_cl(sq_norm).sum() / (A.size() - 1.0)),
      [A, diff](vari& res) mutable {
        A.adj() += res.adj() / (res.val() * (A.size() - 1.0)) * diff;
      });
}

}  // namespace math
}  // namespace stan

#endif
#endif
