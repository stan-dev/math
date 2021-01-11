#ifndef STAN_MATH_OPENCL_REV_LGAMMA_HPP
#define STAN_MATH_OPENCL_REV_LGAMMA_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/value_of.hpp>

namespace stan {
namespace math {

/**
 * Returns the elementwise `lgamma()` of a var_value<matrix_cl<double>>.
 *
 * @param A argument
 * @return Elementwise `lgamma()` of the input.
 */
inline var_value<matrix_cl<double>> lgamma(
    const var_value<matrix_cl<double>>& A) {
  var_value<matrix_cl<double>> res = lgamma(A.val());
  reverse_pass_callback([A, res]() mutable {
    A.adj() = A.adj() + elt_multiply(res.adj(), digamma(A.val()));
  });
  return res;
}

}  // namespace math
}  // namespace stan

#endif
#endif
