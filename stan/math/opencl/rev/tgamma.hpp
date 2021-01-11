#ifndef STAN_MATH_OPENCL_REV_TGAMMA_HPP
#define STAN_MATH_OPENCL_REV_TGAMMA_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/value_of.hpp>

namespace stan {
namespace math {

/**
 * Returns the elementwise `tgamma()` of a var_value<matrix_cl<double>>.
 *
 * @param A argument
 * @return Elementwise `tgamma()` of the input.
 */
inline var_value<matrix_cl<double>> tgamma(
    const var_value<matrix_cl<double>>& A) {
  var_value<matrix_cl<double>> res = tgamma(A.val());

  reverse_pass_callback([A, res]() mutable {
    A.adj()
        = A.adj()
          + elt_multiply(res.adj(), elt_multiply(res.val(), digamma(A.val())));
  });

  return res;
}

}  // namespace math
}  // namespace stan

#endif
#endif
