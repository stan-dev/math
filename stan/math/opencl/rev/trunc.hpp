#ifndef STAN_MATH_OPENCL_REV_TRUNC_HPP
#define STAN_MATH_OPENCL_REV_TRUNC_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/value_of.hpp>

namespace stan {
namespace math {

/**
 * Returns the elementwise `trunc()` of the
 * input `var_value<matrix_cl<double>>`.
 *
 * @param A input `var_value<matrix_cl<double>>`.
 * @return Elementwise `trunc()` of the input argument.
 */
inline var_value<matrix_cl<double>> trunc(
    const var_value<matrix_cl<double>>& A) {
  var_value<matrix_cl<double>> res = trunc(A.val());

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
