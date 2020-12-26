#ifndef STAN_MATH_OPENCL_REV_FLOOR_HPP
#define STAN_MATH_OPENCL_REV_FLOOR_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/value_of.hpp>

namespace stan {
namespace math {

/**
 * Returns the elementwise `floor()` of the input
 * `var_value<matrix_cl<double>>`.
 *
 * @param A input `var_value<matrix_cl<double>>`.
 * @return Elementwise `floor()` of the input argument.
 */
inline var_value<matrix_cl<double>> floor(
    const var_value<matrix_cl<double>>& A) {
  var_value<matrix_cl<double>> res = floor(A.val());

  reverse_pass_callback([A, res]() mutable {
    A.adj() = select(isnan(A.val()), constant(NOT_A_NUMBER, A.rows(), A.cols()),
                     A.adj());
  });

  return res;
}

}  // namespace math
}  // namespace stan

#endif
#endif
