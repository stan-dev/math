#ifndef STAN_MATH_OPENCL_REV_INV_CLOGLOG_HPP
#define STAN_MATH_OPENCL_REV_INV_CLOGLOG_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/value_of.hpp>

namespace stan {
namespace math {

/**
 * Returns the elementwise `inv_cloglog()` of a var_value<matrix_cl<double>>.
 *
 * @param A argument
 * @return Elementwise `inv_cloglog()` of the input.
 */
inline var_value<matrix_cl<double>> inv_cloglog(
    const var_value<matrix_cl<double>>& A) {
  var_value<matrix_cl<double>> res = inv_cloglog(A.val());

  reverse_pass_callback([A, res]() mutable {
    A.adj() = A.adj() + elt_multiply(res.adj(), exp(A.val() - exp(A.val())));
  });

  return res;
}

}  // namespace math
}  // namespace stan

#endif
#endif
