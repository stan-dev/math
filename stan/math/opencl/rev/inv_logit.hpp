#ifndef STAN_MATH_OPENCL_REV_INV_LOGIT_HPP
#define STAN_MATH_OPENCL_REV_INV_LOGIT_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/value_of.hpp>

namespace stan {
namespace math {

/**
 * Returns the elementwise `inv_logit()` of a var_value<matrix_cl<double>>.
 *
 * @param A argument
 * @return Elementwise `inv_logit()` of the input.
 */
inline var_value<matrix_cl<double>> inv_logit(
    const var_value<matrix_cl<double>>& A) {
  var_value<matrix_cl<double>> res = inv_logit(A.val());

  reverse_pass_callback([A, res]() mutable {
    A.adj()
        = A.adj()
          + elt_multiply(res.adj(), elt_multiply(res.val(), 1.0 - res.val()));
  });
  return res;
}

}  // namespace math
}  // namespace stan

#endif
#endif
