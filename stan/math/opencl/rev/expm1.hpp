#ifndef STAN_MATH_OPENCL_REV_EXPM1_HPP
#define STAN_MATH_OPENCL_REV_EXPM1_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/value_of.hpp>

namespace stan {
namespace math {

/**
 * Returns the elementwise `expm1()` of a var_value<matrix_cl<double>>.
 *
 * @param A argument
 * @return Elementwise `expm1()` of the input.
 */
inline var_value<matrix_cl<double>> expm1(
    const var_value<matrix_cl<double>>& A) {
  var_value<matrix_cl<double>> res = expm1(A.val());

  reverse_pass_callback([A, res]() mutable {
    A.adj() = A.adj() + elt_multiply(res.adj(), res.val() + 1.0);
  });

  return res;
}

}  // namespace math
}  // namespace stan

#endif
#endif
