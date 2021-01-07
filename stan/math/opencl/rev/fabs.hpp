#ifndef STAN_MATH_OPENCL_REV_FABS_HPP
#define STAN_MATH_OPENCL_REV_FABS_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/value_of.hpp>

namespace stan {
namespace math {

/**
 * Returns the elementwise `fabs()` of the
 * input `var_value<matrix_cl<double>>`.
 *
 * @param A input `var_value<matrix_cl<double>>`.
 * @return Elementwise `fabs()` of the input argument.
 */
inline var_value<matrix_cl<double>> fabs(
    const var_value<matrix_cl<double>>& A) {
  var_value<matrix_cl<double>> res = fabs(A.val());

  reverse_pass_callback([A, res]() mutable {
    A.adj() = select(
        isnan(A.val()), constant(NOT_A_NUMBER, A.rows(), A.cols()),
        select(A.val() < 0.0, A.adj() - res.adj(), A.adj() + res.adj()));
  });

  return res;
}

}  // namespace math
}  // namespace stan

#endif
#endif
