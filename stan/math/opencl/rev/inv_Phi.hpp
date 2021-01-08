#ifndef STAN_MATH_OPENCL_REV_INV_PHI_HPP
#define STAN_MATH_OPENCL_REV_INV_PHI_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/value_of.hpp>

namespace stan {
namespace math {

/**
 * Returns the elementwise `inv_Phi()` of a var_value<matrix_cl<double>>.
 *
 * @param A argument
 * @return Elementwise `inv_Phi()` of the input.
 */
inline var_value<matrix_cl<double>> inv_Phi(
    const var_value<matrix_cl<double>>& A) {
  var_value<matrix_cl<double>> res = inv_Phi(A.val());

  reverse_pass_callback([A, res]() mutable {
    A.adj() = A.adj()
              + elt_multiply(
                    res.adj(),
                    elt_divide(SQRT_TWO_PI,
                               exp(-0.5 * elt_multiply(res.val(), res.val()))));
  });

  return res;
}

}  // namespace math
}  // namespace stan

#endif
#endif
