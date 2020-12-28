#ifndef STAN_MATH_OPENCL_REV_LOGIT_HPP
#define STAN_MATH_OPENCL_REV_LOGIT_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/value_of.hpp>

namespace stan {
namespace math {

/**
 * Returns the elementwise `logit()` of a var_value<matrix_cl<double>>.
 *
 * @param A argument
 * @return Elementwise `logit()` of the input.
 */
inline var_value<matrix_cl<double>> logit(
    const var_value<matrix_cl<double>>& A) {
  var_value<matrix_cl<double>> res = logit(A.val());

  reverse_pass_callback([A, res]() mutable {
    A.adj() = A.adj() + elt_divide(res.adj(), A.val() * LOG_TEN);
  });

  return res;
}

}  // namespace math
}  // namespace stan

#endif
#endif
